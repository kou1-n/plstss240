      subroutine stress_dp_bn(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens, g_val,
     &                  ierror,  itr , histi )
c
c **********************************************************************
c     Block Newton法によるDrucker-Prager弾塑性構成則
c
c     参考文献: Yamamoto et al. (2021) "Simultaneously iterative procedure
c              based on block Newton method for elastoplastic problems"
c              Int. J. Numer. Methods Eng. 122(9), 2145-2178
c
c     特徴:
c     - 局所反復を行わない（大域的Newton反復のみ）
c     - 平衡方程式と降伏条件を同時に解く
c     - consistency parameter Δγを代数的に更新
c
c     処理の流れ:
c     1. 試行応力の計算（弾性予測子）
c     2. 降伏判定
c     3. 塑性の場合：BOX 1アルゴリズムで応力更新
c     4. 接線剛性の計算（C^ep = C - N^-1 L ⊗ M）
c **********************************************************************
c
      implicit double precision(a-h,o-z)
c
      dimension prope(20)
c
c     --- 共通ブロック（将来の拡張用） ---
      
      dimension str(3,3),stry(3,3),
     &          sig(3,3),oun(3,3),sd(3,3),eel(3,3),
     &          plstrg(3,3), psig(3,3), stri(3,3)
      dimension ctens(3,3,3,3)
      dimension ehist(20), histi(50)
c
c     --- Block Newton variables (1-variable system for D-P) ---
      REAL*8 deltag, alpeg_new, Dg
      REAL*8 deltagi, alpegi, ftri
      REAL*8 deltag_prev, delta_gamma_inc
      REAL*8 gg, xa(3,3)
c     --- L, M, N matrices for Box 2 formulation ---
      dimension L_mat(3,3), M_mat(3,3)
      REAL*8 N_scalar
c     --- Drucker-Prager parameters ---
      REAL*8 eta_dp, xi_dp, etabar_dp, phi_dp, psi_dp
c     --- 作業変数 ---
      REAL*8 g_val, N_inv
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
c
c ***** 履歴変数の読み込み **********************************************
c     ehist: 前ステップからの履歴（塑性ひずみ等）
      alpeg = ehist(1)        ! 等価塑性ひずみ
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          plstrg(jj,ii) = ehist(kk)  ! 塑性ひずみテンソル
        enddo
      enddo
c
c     === 反復履歴の処理（Block Newton法） ===
      if(itr.gt.1) then
c       --- phexa8.fで計算されたδγ増分を取得 ---
        delta_gamma_inc = histi(1)  ! δγ増分（phexa8から）
        deltag_prev = histi(34)     ! 前回の累積Δγ（位置34に移動）
c       累積値を更新: Δγ^(k+1) = Δγ^(k) + δγ
        deltagi = deltag_prev + delta_gamma_inc
        alpegi  = histi(3)          ! 前回のalpeg
        ftri    = histi(4)          ! 前回の降伏関数値
        kk = 4
        do jj=1,3
          do ii=1,3
            kk = kk+1
            stri(ii,jj) = histi(kk)  ! Previous strain
          enddo
        enddo
        do jj=1,3
          do ii=1,3
            kk = kk+1
            xa(ii,jj) = histi(kk)  ! Previous xa (sensitivity)
          enddo
        enddo
      else
c       === 初回反復：BOX 2初期化 ===
        deltagi = 0.d0
        deltag_prev = 0.d0          ! 初回なので累積Δγ=0
        delta_gamma_inc = 0.d0      ! 初回なのでδγ=0
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
      vmu = yng/(2.d0*(1.d0 +poi))                       ! せん断弾性係数
      vlm = poi*yng/((1.d0 +poi)*(1.d0 -2.d0*poi))     ! ラメ第1定数
      vkp = yng/(3.d0*(1.d0 -2.d0*poi))                 ! 体積弾性係数
c
c     === 塑性パラメータ ===
      yld = prope(10)       ! 初期降伏応力
      hk  = prope(11)       ! 硬化係数1
      hpa = prope(12)       ! 硬化係数2
      hpb = prope(13)       ! 硬化係数3
      hpc = prope(14)       ! 硬化係数4
      hpd = prope(15)       ! 硬化係数5
      phi_dp = prope(16)    ! 内部摩擦角
      psi_dp = prope(17)    ! ダイレイタンシー角
c
c     === Drucker-Pragerパラメータの計算 ===
c     Mohr-Coulombからの変換式
      eta_dp   = 6.d0 * sin(phi_dp) / ( dsqrt(3.d0) *
     &                                 (3.d0 - sin(phi_dp)) )
      xi_dp    = 6.d0 * cos(phi_dp) / ( dsqrt(3.d0) *
     &                                 (3.d0 - sin(phi_dp)) )
      etabar_dp = 6.d0 * sin(psi_dp) / ( dsqrt(3.d0) *
     &                                  (3.d0 - sin(psi_dp)) )
c     Drucker-Pragerパラメータは計算済み
c
c ***** 初期化 **********************************************************
      deltag = 0.d0
      psig   = 0.d0
      ctens  = 0.d0
      ierror = 0
c
c   === BOX 1 STEP 1: 試行弾性応力の計算 ===
c   弾性予測子による応力計算
      etrs = str(1,1) +str(2,2) +str(3,3)    ! 体積ひずみ
      emean = etrs/3.d0                       ! 平均ひずみ
c
c   試行偏差応力: s^tr = 2μ{dev[ε] - ε^p}
      stry = 2.d0*vmu*(str -emean*DELTA -plstrg)
c
c   相対応力テンソル: ξ^tr = s^tr - β
c   （現在はβ=0、移動硬化への拡張可能）
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +stry(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)    ! ||s^tr||（偏差応力のノルム）
c
c   流動則の方向ベクトル n = s^tr/||s^tr||
      if(stno.gt.1.0d-16) then
        oun(:,:) = stry(:,:)/stno
      else
        oun(:,:) = 0.d0
      endif
c
c  === 硬化関数と降伏関数の評価 ===
c     Voce型硬化則: H = hpd*(hk*α + (hpa-yld)*(1-exp(-hpb*α)))
      hard = hpd*(hk*alpeg
     &     +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg)))
      dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg))  ! dH/dα
c
c     降伏関数: f = √(1/2)||s|| + η*p - ξ(σ_Y + H)
      ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*vkp*etrs
     &      - xi_dp*(yld +hard)
c
c  ===== 塑性状態：Block Newton法による応力更新 =====
      if(ftreg.gt.-ctol) then
        idepg = 1
c
        if(itr.eq.1) then
c         === BOX 2初期化（論文のBOX 2 STEP 1） ===
c
c         1. Δγ⁽⁰⁾ = 0で初期化
          deltag = 0.d0
c
c         2. 初期降伏関数値の評価
          hard = hpd*(hk*alpeg
     &         +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg)))
          gg = dsqrt(1.d0/2.d0)*stno
     &       + eta_dp*vkp*etrs
     &       - xi_dp*(yld + hard)
c
c         3. ε^p(u_{n+1}^{(0)}, γ_{n+1}^{(0)}) = ε^p(u_n, γ_n)
c            (plstrg already loaded from ehist - no change needed)
c         
c         4. α(γ_{n+1}^{(0)}) = α(γ_n)  
          alpeg_new = alpeg
c         
c         5. β(u_{n+1}^{(0)}, γ_{n+1}^{(0)}) = β(u_n, γ_n)
c            (For D-P: kinematic hardening handled through α relationship)
c         
c         6. n(u_{n+1}^{(0)}, γ_{n+1}^{(0)}) = n(u_n, γ_n)
c            (oun computed from current trial stress - keep existing logic)
c         
c         === BOX 2 STEP 2: Compute tangent moduli for first calculation ===
c         Compute L, M, N matrices according to paper equations (48)-(50)
c         
c         L = -C^e : (n + α/3*I) = -2μn - κα/3*I (Drucker-Prager)
          do jj=1,3
            do ii=1,3
              L_mat(ii,jj) = -2.d0*vmu*oun(ii,jj)
     &                       - eta_dp*vkp*DELTA(ii,jj)/3.d0
            enddo
          enddo
c         
c         M = 2μ(n + β/3*I) = 2μn + 2μβ/3*I (Drucker-Prager)
          do jj=1,3
            do ii=1,3
              M_mat(ii,jj) = 2.d0*vmu*oun(ii,jj)
     &                       + 2.d0*vmu*etabar_dp*DELTA(ii,jj)/3.d0
            enddo
          enddo
c
c         N = -{2μ + κη̄η + ξ²K'} (Drucker-Prager specific)
          N_scalar = -(2.d0*vmu + vkp*etabar_dp*eta_dp
     &              + xi_dp*xi_dp*dhard)
c         Prevent singular matrix (ensure N_scalar is not too small)
          if(dabs(N_scalar).lt.1.d-12) then
            N_scalar = -2.d0*vmu    ! Use elastic value as fallback
          endif
c         
c         Keep Dg for backward compatibility
          Dg = N_scalar
        else
c         === BOX 1 STEP 2: Block Newton (no local iteration) ===
c         Following Yamamoto et al.: direct state update without local Newton
c
c         === Evaluate with current deltagi (from global system) ===
c         No local Newton iteration - use global coupling result
c
c         === Update L, M, N matrices for current iteration ===
c         L = -C^e : (n + α/3*I) = -2μn - κα/3*I (Drucker-Prager)
          do jj=1,3
            do ii=1,3
              L_mat(ii,jj) = -2.d0*vmu*oun(ii,jj)
     &                       - eta_dp*vkp*DELTA(ii,jj)/3.d0
            enddo
          enddo
c
c         M = 2μ(n + β/3*I) = 2μn + 2μβ/3*I (Drucker-Prager)
          do jj=1,3
            do ii=1,3
              M_mat(ii,jj) = 2.d0*vmu*oun(ii,jj)
     &                       + 2.d0*vmu*etabar_dp*DELTA(ii,jj)/3.d0
            enddo
          enddo
c
c         === BOX 1 STEP 2: Check consistency parameter Δγ^(k+1) ===
          deltag = deltagi
          if(deltag.lt.0.d0) deltag = 0.d0
c         Removed debug code for minimal implementation

c         === BOX 1: Update state variables ===
c         Update plastic strain: εᵖ = εᵖₙ + Δγ·n
          plstrg(:,:) = plstrg(:,:) + deltag*oun(:,:)*dsqrt(1.d0/2.d0)
     &                + deltag*etabar_dp*DELTA(:,:)/3.d0
c         Update equivalent plastic strain: α = αₙ + √(2/3)·Δγ
          alpeg_new = alpeg + xi_dp*deltag
c         Update hardening derivative for updated state
          dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg_new))
c         Update N_scalar for current state (should be negative)
c         N = -(2μ + κη̄η + ξ²K')
          N_scalar = -(2.d0*vmu + vkp*etabar_dp*eta_dp
     &              + xi_dp*xi_dp*dhard)
c         Ensure N_scalar is numerically stable
          if(dabs(N_scalar).lt.1.d-10) then
            N_scalar = -2.d0*vmu  ! Use elastic value as fallback
          endif
c         Update stress: σ = κtr[ε]I + s^tr - √2μΔγn
          sig(:,:) = stry(:,:) - dsqrt(2.d0)*vmu*deltag*oun(:,:)
     &             + (vkp*etrs - vkp*etabar_dp*deltag)*DELTA(:,:)

c       === BOX 1: Evaluate yield function with updated state ===
        smean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
        stno = 0.d0
        do ii=1,3
          do jj=1,3
            stno = stno +(sig(ii,jj) - smean*DELTA(ii,jj))**2
          enddo
        enddo
        stno = dsqrt(stno)
        hard = hpd*(hk*alpeg_new
     &       +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg_new)))
c         Removed debug code for minimal implementation
        g_val = dsqrt(1.d0/2.d0)*stno + eta_dp*smean
     &        - xi_dp*(yld + hard)

c         Removed debug code for minimal implementation
c
c       === BOX 1 STEP 3: Compute tangent moduli ===
c       Use already computed L, M, N matrices
c
c       === Compute C matrix (with geometric softening) ===
        if(deltag.gt.1.d-16 .and. stno.gt.1.d-16) then
          theta = 1.d0 - (dsqrt(2.d0)*vmu*deltag)/stno
        else
          theta = 1.d0
        endif
c
c       === Compute C^ep = C - N^-1 L ⊗ M (Paper equation) ===
c       N_scalar should be negative, so N_inv should be negative
        if(dabs(N_scalar).gt.1.d-10) then
          N_inv = 1.d0 / N_scalar
        else
c         Use elastic tangent if N_scalar is too small
          N_inv = 0.d0
        endif
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
c               C matrix (elastic + geometric softening)
                C_ijkl = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &                 + 2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &                 - (1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
c               
c               C^ep = C - N^-1 L ⊗ M (Paper formulation)
                ctens(ii,jj,kk,ll) = C_ijkl 
     &                             - N_inv * L_mat(ii,jj) * M_mat(kk,ll)
              enddo
            enddo
          enddo
        enddo
c
c       === Store values for next iteration ===
c       NOTE: histi(1)はphexa8からのdelta_gamma_incで使用されるため上書きしない
        histi(34) = deltag      ! 累積Δγを位置34に保存（次回のdeltag_prev）
        histi(3) = alpeg_new    ! 更新された相当塑性ひずみ
        histi(4) = ftreg        ! 降伏関数値
c       === NEW: Store N_scalar and g_val for delta_gamma calculation ===
        histi(5) = N_scalar
        histi(6) = g_val
        kk = 6
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = str(ii,jj)
          enddo
        enddo
c       Store xa (sensitivity ∂Δγ/∂ε) for next iteration
c       xa = (1/Dg) * (∂gg/∂ε)
        if(dabs(N_scalar).gt.1.d-16) then
          xa(:,:) = (1.d0/N_scalar)*(dsqrt(1.d0/2.d0)*oun(:,:)
     &                        + eta_dp*DELTA(:,:)/3.d0)
        else
          xa(:,:) = 0.d0
        endif
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = xa(ii,jj)
          enddo
        enddo
c       === NEW: Store M_mat for delta_gamma calculation in phexa8 ===
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = M_mat(ii,jj)
          enddo
        enddo
c
c       Update alpeg for storage
        alpeg = alpeg_new
      endif
c
c  ===== ELASTIC CASE =====
      else
        idepg = 0
        deltag = 0.d0
c
c       === Elastic case: set yield function to zero (Paper eq. 84) ===
c       Following Yamamoto et al. Box 1: g ≡ 0 for elastic state
        gg = 0.d0
c
c       === Elastic stress ===
        sig(:,:) = stry(:,:) +vkp*etrs*DELTA(:,:)
c         Removed debug code for minimal implementation
c
c       === Elastic tangent moduli ===
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &              +2.d0*vmu*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
              enddo
            enddo
          enddo
        enddo
c
c       Store elastic values
c       NOTE: histi(1)はphexa8からのdelta_gamma_incで使用されるため上書きしない
        histi(34) = 0.d0        ! 弾性状態ではΔγ=0（位置34）
        histi(3) = alpeg        ! 現在の相当塑性ひずみ（変化なし）
        histi(4) = 0.d0         ! 弾性状態では降伏関数=0
c       Store yield function value: 0 for elastic case (Paper eq. 84)
        g_val = gg
      endif
c
c   === Compute von Mises Stress ===
      smean = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
      sd(:,:) = sig(:,:) -smean*DELTA(:,:)
c
      vons = 0.d0
      do jj=1,3
        do ii=1,3
          vons = vons +sd(ii,jj)*sd(ii,jj)
        enddo
      enddo
      vons = dsqrt(1.5d0*vons)
c
c   === Compute the Energy Density ===
      eel(:,:) = str(:,:) -plstrg(:,:)
c
c     --- compute the elastic strain energy density
      e_dns = 0.d0
      do jj=1,3
        do ii=1,3
          e_dns = e_dns +eel(ii,jj)*sig(ii,jj)
        enddo
      enddo
      e_dns = e_dns*0.5d0
c
c     --- compute the plastic strain energy density
      p_dns = 0.d0
      do jj=1,3
        do ii=1,3
          p_dns = p_dns+deltag*oun(ii,jj)*sig(ii,jj)
        enddo
      enddo
c
c ***** Store Deformation Histories ************************************
      ehist(1) = alpeg
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          ehist(kk) = plstrg(jj,ii)
        enddo
      enddo
c
c **********************************************************************
      RETURN
      END