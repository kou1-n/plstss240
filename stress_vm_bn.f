      subroutine stress_vm_bn(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens, g_val,
     &                  ierror, itr, histi )
c
c     ============================================================
c     Block Newton法によるvon Mises弾塑性応力更新
c     山本ら(2021)の論文に基づく実装
c     ============================================================
c
      implicit double precision(a-h,o-z)
c
      dimension prope(20)
      dimension str(3,3),stry(3,3),seta(3,3)
      dimension sig(3,3),oun(3,3),sd(3,3),eel(3,3)
      dimension plstrg(3,3),betaeg(3,3),stri(3,3)
      dimension sxi(3,3)  ! 相対応力ξ
      dimension ctens(3,3,3,3)
      dimension ehist(20), histi(50)
c
c     --- Block Newton用変数 (論文のBox 1準拠) ---
      REAL*8 deltag, alpeg_new, theta
      REAL*8 deltagi, alpegi, ftri
      REAL*8 deltag_prev, delta_gamma_inc
      REAL*8 g_val, xa(3,3)
c     --- L, M, N行列 (論文式48-50) ---
      dimension L_mat(3,3), M_mat(3,3)
      REAL*8 N_scalar, N_inv, stno_new, hard_new
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),
     &               DTENS(3,3,3,3)
c
c **********************************************************************
c
c ***** 履歴データの読み込み *******************************************
c     === 内部変数の復元（前ステップから） ===
      alpeg = ehist(1)           ! 相当塑性ひずみ
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          plstrg(jj,ii) = ehist(kk)  ! 塑性ひずみテンソル
        enddo
      enddo
      do ii=1,3
        do jj=1,3
          kk = kk+1
          betaeg(jj,ii) = ehist(kk)  ! 移動硬化背応力
        enddo
      enddo
c
c     === 反復履歴の処理（Block Newton法） ===
      if(itr.gt.1) then
c       --- phexa8.fで計算されたδγ増分を取得 ---
        delta_gamma_inc = histi(1)  ! δγ増分（phexa8から）
        deltag_prev = histi(34)     ! 前回の累積Δγ（位置34）
c       累積値を更新: Δγ^(k+1) = Δγ^(k) + δγ
        deltagi = deltag_prev + delta_gamma_inc
        alpegi  = histi(3)          ! 前回のalpeg
        ftri    = histi(4)          ! 前回の降伏関数値
        kk = 4
        do jj=1,3
          do ii=1,3
            kk = kk+1
            stri(ii,jj) = histi(kk)  ! 前回のひずみ
          enddo
        enddo
        do jj=1,3
          do ii=1,3
            kk = kk+1
            xa(ii,jj) = histi(kk)    ! 感度行列
          enddo
        enddo
      else
c       === 初回反復：BOX 2初期化 ===
        deltagi = 0.d0
        deltag_prev = 0.d0
        delta_gamma_inc = 0.d0
        alpegi  = alpeg
        ftri    = 0.d0
        stri    = str
        xa      = 0.d0
c       初回反復でもhisti(34)を初期化
        histi(34) = 0.d0            ! 初回の累積Δγ=0
      endif
c
c ***** 材料定数の設定 **************************************************
c     === 弾性定数 ===
      yng = prope(1)        ! ヤング率
      poi = prope(2)        ! ポアソン比
      vmu = yng/(2.d0*(1.d0 +poi))               ! せん断弾性係数
      vlm = poi*yng/((1.d0 +poi)*(1.d0 -2.d0*poi))  ! ラメ第1定数
      vkp = yng/(3.d0*(1.d0 -2.d0*poi))          ! 体積弾性係数
c
c     === 塑性パラメータ ===
      yld = prope(10)       ! 初期降伏応力
      hk  = prope(11)       ! 線形硬化係数
      hpa = prope(12)       ! Voce硬化飽和応力
      hpb = prope(13)       ! Voce硬化指数
      hpd = prope(15)       ! 硬化スケール係数
c
c ***** 初期化 **********************************************************
      deltag = 0.d0
      ctens  = 0.d0
      ierror = 0
c
c  === BOX 1 STEP 1: 試行弾性応力の計算 ===
c  弾性予測子による応力計算
      etrs = str(1,1) + str(2,2) + str(3,3)    ! 体積ひずみ
      emean = etrs/3.d0                         ! 平均ひずみ
c
c  試行偏差応力: s^tr = 2μ{dev[ε] - ε^p_n}
      stry = 2.d0*vmu*(str - emean*DELTA - plstrg)
c
c  相対応力テンソル: ξ^tr = s^tr - β_n
      seta = stry - betaeg
c
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno + seta(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)    ! ||ξ^tr||（相対応力のノルム）
c
c  流動則の方向ベクトル n = ξ^tr/||ξ^tr||
      if(stno.gt.1.0d-16) then
        oun(:,:) = seta(:,:)/stno
      else
        oun(:,:) = 0.d0
      endif
c
c  === 硬化関数と降伏関数の評価（試行応力での） ===
c     Voce型硬化則: K = hpd*(hk*α + (hpa-yld)*(1-exp(-hpb*α)))
      hard = hpd*(hk*alpeg + (hpa - yld)*(1.d0 - dexp(-hpb*alpeg)))
      dhard = hpd*(hk + hpb*(hpa - yld)*dexp(-hpb*alpeg))  ! dK/dα
c
c     降伏関数: f = ||ξ|| - √(2/3)(σ_Y + K)
      ftreg = stno - dsqrt(2.d0/3.d0)*(yld + hard)
c
c  ===== 塑性状態：Block Newton法による応力更新 =====
      if(ftreg.gt.-ctol) then
        idepg = 1
c
c       === BOX 1 STEP 3: Update consistency parameter ===
        if(itr.eq.1) then
c         === BOX 2: Initialization (k=0) ===
c
c         --- BOX 2 STEP 1: Set initial consistency parameter ---
          deltag = 0.d0
c
c         --- BOX 2 STEP 2: Compute initial stress (elastic predictor) ---
c         σ_{n+1}^(0) = s^tr + κtr(ε)I
          sig(:,:) = stry(:,:) + vkp*etrs*DELTA(:,:)
c
c         --- BOX 2 STEP 3: Evaluate initial yield function ---
c         s^(0) = dev(σ_{n+1}^(0))
          smean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
          sd(:,:) = sig(:,:) - smean*DELTA(:,:)
c
c         ξ^(0) = s^(0) - β_n
          kk = 10
          do ii=1,3
            do jj=1,3
              kk = kk+1
              betaeg(jj,ii) = ehist(kk)  ! β_n from previous step
            enddo
          enddo
          sxi(:,:) = sd(:,:) - betaeg(:,:)
c
c         ||ξ^(0)||
          stno_new = 0.d0
          do ii=1,3
            do jj=1,3
              stno_new = stno_new + sxi(ii,jj)**2
            enddo
          enddo
          stno_new = dsqrt(stno_new)
c
c         K(α_n)
          hard = hpd*(hk*alpeg + (hpa - yld)*(1.d0 - dexp(-hpb*alpeg)))
c
c         g^(0) = ||ξ^(0)|| - √(2/3)(σ_Y + K(α_n))
          g_val = stno_new - dsqrt(2.d0/3.d0)*(yld + hard)
c
c         --- BOX 2 STEP 4: Set initial plastic variables ---
c         ε^p_{n+1}^(0) = ε^p_n (already in plstrg from ehist)
c         α_{n+1}^(0) = α_n
          alpeg_new = alpeg
c         β_{n+1}^(0) = β_n (already in betaeg)
c
c         --- BOX 2 STEP 5: Compute initial flow direction ---
c         n^(0) = ξ^(0)/||ξ^(0)||
          if(stno_new.gt.1.d-16) then
            do jj=1,3
              do ii=1,3
                oun(ii,jj) = sxi(ii,jj)/stno_new
              enddo
            enddo
          else
            oun(:,:) = 0.d0
          endif
c
c         --- BOX 2 STEP 6: Initialize matrices for tangent computation ---
c         (L, M, N matrices will be computed below in Step 6)
c
        else
c         === BOX 1 STEP 3: Update consistency parameter ===
c         Δγ^(k+1) = Δγ^(k) + δγ (phexa8から渡されたδγ増分を使用)
          deltag = deltagi
          if(deltag.lt.0.d0) deltag = 0.d0
c
c         === BOX 1 STEP 4: Update state variables ===
c         Step 4.1: 相当塑性ひずみの更新 α_{n+1}^(k) = α_n + √(2/3)Δγ^(k)
          alpeg_new = alpeg + dsqrt(2.d0/3.d0)*deltag
c
c         Step 4.2: 硬化係数の計算 K'(α_{n+1}^(k))
          dhard = hpd*(hk + hpb*(hpa - yld)*dexp(-hpb*alpeg_new))
c
c         Step 4.3: 仮の応力から相対応力を計算してn^(k)を更新
c         まず仮の応力を計算: σ_temp = s^tr - 2μΔγn_old + κtr(ε)I
          sig(:,:) = stry(:,:) - 2.d0*vmu*deltag*oun(:,:)
     &             + vkp*etrs*DELTA(:,:)
c
c         偏差応力を計算
          smean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
          sd(:,:) = sig(:,:) - smean*DELTA(:,:)
c
c         背応力の更新: β_{n+1}^(k) = β_n + (2/3)K'Δγ^(k)n^(k-1)
          kk = 10
          do ii=1,3
            do jj=1,3
              kk = kk+1
              betaeg(jj,ii) = ehist(kk)
     &                      + (2.d0/3.d0)*dhard*deltag*oun(jj,ii)
            enddo
          enddo
c
c         相対応力ξ^(k) = s^(k) - β_{n+1}^(k)
          sxi(:,:) = sd(:,:) - betaeg(:,:)
          stno_new = 0.d0
          do ii=1,3
            do jj=1,3
              stno_new = stno_new + sxi(ii,jj)**2
            enddo
          enddo
          stno_new = dsqrt(stno_new)
c
c         Step 4.4: 流動方向の更新 n^(k) = ξ^(k)/||ξ^(k)||
          if(stno_new.gt.1.d-16) then
            do jj=1,3
              do ii=1,3
                oun(ii,jj) = sxi(ii,jj)/stno_new
              enddo
            enddo
          else
            oun(:,:) = 0.d0
          endif
c
c         Step 4.5: 塑性ひずみの更新 ε^p_{n+1}^(k) = ε^p_n + Δγ^(k)n^(k)
c         NOTE: 最新のn^(k)を使用
          kk = 1
          do ii=1,3
            do jj=1,3
              kk = kk+1
              plstrg(jj,ii) = ehist(kk)
     &                      + deltag*oun(jj,ii)
            enddo
          enddo
c
c         Step 4.6: 背応力を最新のn^(k)で再計算
          kk = 10
          do ii=1,3
            do jj=1,3
              kk = kk+1
              betaeg(jj,ii) = ehist(kk)
     &                      + (2.d0/3.d0)*dhard*deltag*oun(jj,ii)
            enddo
          enddo
c
c         === BOX 1 STEP 5: Update stress and compute yield function ===
c         Step 5.1: 応力更新 σ_{n+1}^(k) = s^tr - 2μΔγ^(k)n^(k) + κtr(ε)I
c         NOTE: 最新のn^(k)を使用して応力を再計算
          sig(:,:) = stry(:,:) - 2.d0*vmu*deltag*oun(:,:)
     &             + vkp*etrs*DELTA(:,:)
c
c         Step 5.2: 偏差応力と相対応力の計算
          smean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
          sd(:,:) = sig(:,:) - smean*DELTA(:,:)
c
c         Step 5.3: ξ^(k) = dev(σ_{n+1}^(k)) - β_{n+1}^(k)
          sxi(:,:) = sd(:,:) - betaeg(:,:)
          stno_new = 0.d0
          do ii=1,3
            do jj=1,3
              stno_new = stno_new + sxi(ii,jj)**2
            enddo
          enddo
          stno_new = dsqrt(stno_new)
c
c         Step 5.4: 硬化関数K(α_{n+1}^(k))
          hard_new = hpd*(hk*alpeg_new
     &             + (hpa - yld)*(1.d0 - dexp(-hpb*alpeg_new)))
c
c         Step 5.5: 降伏関数 g^(k) = ||ξ^(k)|| - √(2/3)(σ_Y + K)
          g_val = stno_new - dsqrt(2.d0/3.d0)*(yld + hard_new)
        endif
c
c
c       === BOX 1 STEP 6: Compute tangent modulus ===
c       L, M, N行列の計算（接線剛性用）
c
c       L = -2μn (von Mises)
        do jj=1,3
          do ii=1,3
            L_mat(ii,jj) = -2.d0*vmu*oun(ii,jj)
          enddo
        enddo
c
c       M = 2μn (von Mises)
        do jj=1,3
          do ii=1,3
            M_mat(ii,jj) = 2.d0*vmu*oun(ii,jj)
          enddo
        enddo
c
c       N = -(2μ + (2/3)K')
        N_scalar = -(2.d0*vmu + (2.d0/3.d0)*dhard)
c       数値安定性の確保
        if(dabs(N_scalar).lt.1.d-10) then
          N_scalar = -2.d0*vmu
        endif
c
c       θ = 1（簡単化のため幾何学的軟化を無視）
        theta = 1.d0
c
c       === C^ep = C - N^-1 L ⊗ M（論文式51） ===
c       N_scalarが小さすぎる場合は弾性接線を使用
        if(dabs(N_scalar).lt.1.d-10) then
          N_inv = 0.d0
        else
          N_inv = 1.d0 / N_scalar
        endif
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
c               C行列（弾性）
                C_ijkl = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &                 + 2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &                 - (1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
c               塑性補正項を適用（N_inv=0の時は弾性のみ）
                ctens(ii,jj,kk,ll) = C_ijkl
     &                      - N_inv * L_mat(ii,jj) * M_mat(kk,ll)
              enddo
            enddo
          enddo
        enddo
c
c       === 次回反復用の値を保存 ===
c       NOTE: histi(1)はphexa8からのdelta_gamma_incで使用されるため上書きしない
        histi(34) = deltag      ! 累積Δγを位置34に保存（次回のdeltag_prev）
        histi(3) = alpeg_new    ! 更新された相当塑性ひずみ
        histi(4) = g_val        ! 降伏関数値
c       === N_scalarとg_valをδγ計算用に保存 ===
        histi(5) = N_scalar
        histi(6) = g_val
        kk = 6
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = str(ii,jj)  ! 現在のひずみ
          enddo
        enddo
c       感度行列xa = (1/N) * (∂g/∂ε)を保存
        if(dabs(N_scalar).gt.1.d-16) then
          xa(:,:) = (1.d0/N_scalar)*oun(:,:)
        else
          xa(:,:) = 0.d0
        endif
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = xa(ii,jj)
          enddo
        enddo
c       === M行列をphexa8でのδγ計算用に保存 ===
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = M_mat(ii,jj)
          enddo
        enddo
c
c       履歴変数の更新
        alpeg = alpeg_new
c
c  ===== 弾性状態 =====
      else
        idepg = 0
        deltag = 0.d0
c
c       === 弾性状態：降伏関数をゼロとする（論文式84） ===
c       山本らBox 1: g ≡ 0 for elastic state
        g_val = 0.d0
c
c       === 弾性応力 ===
        sig(:,:) = stry(:,:) + vkp*etrs*DELTA(:,:)
c
c       === 弾性接線剛性 ===
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll) = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &            + 2.d0*vmu*( FIT(ii,jj,kk,ll)
     &                       - (1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
              enddo
            enddo
          enddo
        enddo
c
c       弾性状態の値を保存
c       NOTE: histi(1)はphexa8からのdelta_gamma_incで使用されるため上書きしない
        histi(34) = 0.d0        ! 弾性状態ではΔγ=0（位置34）
        histi(3) = alpeg        ! 現在の相当塑性ひずみ（変化なし）
        histi(4) = 0.d0         ! 弾性状態では降伏関数=0
        histi(5) = -2.d0*vmu    ! N_scalar（弾性値）
        histi(6) = 0.d0         ! g_val = 0（弾性）
      endif
c
c ***** 内部変数の保存 **************************************************
      ehist(1) = alpeg
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          ehist(kk) = plstrg(jj,ii)
        enddo
      enddo
      do ii=1,3
        do jj=1,3
          kk = kk+1
          ehist(kk) = betaeg(jj,ii)
        enddo
      enddo
c
c ***** その他の出力値 **************************************************
c     von Mises応力の計算
      smean = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
      sd(:,:) = sig(:,:) -smean*DELTA(:,:)
      vons = 0.d0
      do jj=1,3
        do ii=1,3
          vons = vons +sd(ii,jj)*sd(ii,jj)
        enddo
      enddo
      vons = dsqrt(1.5d0*vons)
c
c     エネルギー密度の計算
      eel(:,:) = str(:,:) -plstrg(:,:)
      e_dns = 0.d0
      do jj=1,3
        do ii=1,3
          e_dns = e_dns +eel(ii,jj)*sig(ii,jj)
        enddo
      enddo
      e_dns = e_dns*0.5d0
c
      p_dns = 0.d0
      do jj=1,3
        do ii=1,3
          p_dns = p_dns+deltag*oun(ii,jj)*sig(ii,jj)
        enddo
      enddo
c
      return
      end