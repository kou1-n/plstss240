      subroutine phexa8(   nel,  nelx,   ndf,  node, ngaus,
     &                    idep,MATYPE, prope,    xe,    ue,
     &                     sts,   stn, finte,  vone,  epse,
     &                  dhist0,dhist1,
     &                   e_ene, p_ene,tene_e,dene_p, tempe,
     &                   dtemp,dtempe,
     &                  dndx_g, det_g,ctensg, g_vals,
     &                  ierror, itr, histi0 )
c
c **********************************************************************
c     8節点六面体要素の要素計算サブルーチン
c
c     機能:
c     1. ガウス積分点での応力・ひずみ計算
c     2. 要素剛性行列の計算
c     3. 内力ベクトルの計算
c     4. Block Newton法でのδγ計算（MATYPE=5の場合）
c
c     入力:
c     - nel: 要素番号
c     - xe: 要素節点座標
c     - ue: 要素節点変位
c     - prope: 材料定数
c     - MATYPE: 材料モデルタイプ
c     - itr: 反復回数
c
c     出力:
c     - finte: 内力ベクトル
c     - sts, stn: 要素平均応力・ひずみ
c     - ctensg: 接線剛性テンソル
c     - g_vals: 降伏関数値（Block Newton用）
c **********************************************************************
c
      implicit double precision (a-h,o-z)
c
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
      dimension dhist0(20,ngaus,nelx),dhist1(20,ngaus,nelx)
      dimension histi0(50,ngaus,nelx)
      dimension g_vals(ngaus)
c
      dimension idep(8)
c
      dimension prope(20)
      dimension xe(3,8),dn(3,8),dndx(3,8),ue(3,8)
      dimension wg(2),xg(2)
      dimension dxdl(3,3),dldx(3,3),
     &          sig(3,3),sts(3,3),str(3,3),stn(3,3),
     &          sd(3,3),stry(3,3),oun(3,3),seta(3,3),
     &          str_prev(3,3), delta_str(3,3), M_mat(3,3)
      dimension ctens(3,3,3,3)
      REAL*8 N_scalar, g_val_prev, M_eps, delta_gamma_inc
      dimension finte(24)
      dimension ehist(20)
      dimension histi(50)
c
      common /tvalu/ ctol,stol
      common /debug_info/ nel_current, ig_current, lstep_current
c **********************************************************************
c　収束が悪いのでitrmaxを増やした2025年8月23日
      itrmax = 50
c
      wg(1) = 1.d0
      wg(2) = 1.d0
      xg(1) = -1.d0/dsqrt(3.d0)
      xg(2) =  1.d0/dsqrt(3.d0)
c
c ***** Initilization **************************************************
      finte = 0.d0 ! finte(:) = 0.d0
c
      area = 0.d0
      vone = 0.d0
      epse = 0.d0
      p_ene = 0.d0
      e_ene = 0.d0
c
      sts = 0.d0 ! sts(:,:) = 0.d0
      stn = 0.d0 ! stn(:,:) = 0.d0
c     Initialize yield function values for BN method
      g_vals = 0.d0
c
c ***** Set Material Properties ****************************************
c    --- thermal parameters
c     row = prope(3)
c     ccc = prope(4)
c     aaa = prope(5)
c
c ***** ガウス積分開始（2×2×2=8点積分） ****************************
c     各積分点で以下を実行:
c     1. 形状関数の微分を計算
c     2. ひずみテンソルを計算
c     3. 構成則を呼び出して応力と接線剛性を取得
c     4. 内力ベクトルと剛性行列に積分
      ig = 0
      do 100 ix=1,2
        xl = xg(ix)
        wx = wg(ix)
        do 110 iy=1,2
          yl = xg(iy)
          wy = wg(iy)
          do 120 iz=1,2
            zl = xg(iz)
            wz = wg(iz)
            ig = ig +1
c
c         === Initilization of Local Tensors ===
            sig = 0.d0 ! sig(:,:) = 0.d0
            str = 0.d0 ! str(:,:) = 0.d0
c
c         === Derivative of Shape Function r.w.t. Local Coordinate ===
c             ( dn_ij = d(N^j)/d(l_i) )
c           dn(1,1) = -0.125d0*(1.d0 -yl)*(1.d0 -zl)
c           dn(2,1) = -0.125d0*(1.d0 -xl)*(1.d0 -zl)
c           dn(3,1) = -0.125d0*(1.d0 -xl)*(1.d0 -yl)
c           dn(1,2) =  0.125d0*(1.d0 -yl)*(1.d0 -zl)
c           dn(2,2) = -0.125d0*(1.d0 +xl)*(1.d0 -zl)
c           dn(3,2) = -0.125d0*(1.d0 +xl)*(1.d0 -yl)
c           dn(1,3) =  0.125d0*(1.d0 +yl)*(1.d0 -zl)
c           dn(2,3) =  0.125d0*(1.d0 +xl)*(1.d0 -zl)
c           dn(3,3) = -0.125d0*(1.d0 +xl)*(1.d0 +yl)
c           dn(1,4) = -0.125d0*(1.d0 +yl)*(1.d0 -zl)
c           dn(2,4) =  0.125d0*(1.d0 -xl)*(1.d0 -zl)
c           dn(3,4) = -0.125d0*(1.d0 -xl)*(1.d0 +yl)
c           dn(1,5) = -0.125d0*(1.d0 -yl)*(1.d0 +zl)
c           dn(2,5) = -0.125d0*(1.d0 -xl)*(1.d0 +zl)
c           dn(3,5) =  0.125d0*(1.d0 -xl)*(1.d0 -yl)
c           dn(1,6) =  0.125d0*(1.d0 -yl)*(1.d0 +zl)
c           dn(2,6) = -0.125d0*(1.d0 +xl)*(1.d0 +zl)
c           dn(3,6) =  0.125d0*(1.d0 +xl)*(1.d0 -yl)
c           dn(1,7) =  0.125d0*(1.d0 +yl)*(1.d0 +zl)
c           dn(2,7) =  0.125d0*(1.d0 +xl)*(1.d0 +zl)
c           dn(3,7) =  0.125d0*(1.d0 +xl)*(1.d0 +yl)
c           dn(1,8) = -0.125d0*(1.d0 +yl)*(1.d0 +zl)
c           dn(2,8) =  0.125d0*(1.d0 -xl)*(1.d0 +zl)
c           dn(3,8) =  0.125d0*(1.d0 -xl)*(1.d0 +yl)
c
c       === Derivative of Current Position r.w.t. Loacal Coordinate ===
c           ( dxdl_ij = d(x_i)/d(l_j) )
c           do ndia=1,ndf
c             do ndjb=1,ndf
c               dxdlj = 0.d0
c               do ia=1,node
c                 dxdlj = dxdlj +xe(ndia,ia)*dn(ndjb,ia)
c               enddo
c               dxdl(ndia,ndjb) = dxdlj
c             enddo
c           enddo
c
c         === Derivative of Local Coordinate r.w.t. Current Position ===
c             ( dldx_ij = d(li)/d(xj) = ( d(xi)/d(lj) )^{-1} )
c           CALL inv_33(  dxdl,   det,  dldx,ierror )
c           if(det.le.0.d0) CALL zerodt(   nel,    ig,   det)
c
            det = det_g(ig,nel)
c           det_g(ig,nel) = det
            detwxyz = det*wx*wy*wz
            area = area +detwxyz
c
c       === Derivative of Local Coordinate r.w.t. Current Position ===
c           ( dldx_ij = (dxdl_ij)^{-1} )
c           deti = 1.d0/(det*2.d0)
c           do jj=1,3
c             do ii=1,3
c               dldx(ii,jj) = dldx(ii,jj)*deti
c             enddo
c           enddo
c
c       === Derivative of Shape Function r.w.t. Current Position ===
c           ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
            dndx(:,:) = dndx_g(:,:,ig,nel)
c           do ndia=1,ndf
c             do ia=1,node
c               dndxj = 0.d0
c               do ndjb=1,ndf
c                 dndxj = dndxj +dn(ndjb,ia)*dldx(ndjb,ndia)
c               enddo
c               dndx(ndia,ia) = dndxj
c               dndx_g(ndia,ia,ig,nel) = dndxj
c             enddo
c           enddo
c
c       === ひずみテンソルの計算 ===
c           ε_ij = 1/2(u_i,j + u_j,i)
c           節点変位と形状関数微分から計算
            do jj=1,ndf
              do ii=1,ndf
                eij = 0.d0
                do no=1,node
                  eij = eij +0.5d0*(
     &                  ue(ii,no)*dndx(jj,no) +ue(jj,no)*dndx(ii,no) )
                enddo
                str(ii,jj) = eij
              enddo
            enddo
c
c         === Set Deformation Histtories ===
            ehist(:) = dhist0(:,ig,nel)
            histi(:) = histi0(:,ig,nel)
c
c         === 構成則の呼び出し（材料モデルに応じて分岐） ===
          if( MATYPE.eq.1) then
            CALL stress(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror, itr, histi )
          elseif(MATYPE.eq.2) then
            CALL st_gtn(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror, itr, histi )
          elseif(MATYPE.eq.3) then
            CALL stress_vm(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror, itr, histi )
          elseif(MATYPE.eq.4) then
            CALL stress_dp_rm(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror )
          elseif(MATYPE.eq.5) then
c           ===================================================
c           Block Newton法の実装（山本ら2021年の論文に基づく）
c           ===================================================
            nel_current = nel
            ig_current = ig
c
c           === 反復2回目以降：δγの計算 ===
c           論文のBOX 2 Step 3.2: δγ = -N^{-1}(M:ε(δu) + g)
            if(itr.gt.1) then
c             前回反復から保存された値を取得
              N_scalar = histi(5)     ! N係数（負値）
              g_val_prev = histi(6)   ! 前回の降伏関数値
c
c             === ひずみ増分の計算：δε = ε^(k) - ε^(k-1) ===
              kk = 6
              do jj=1,3
                do ii=1,3
                  kk = kk+1
                  str_prev(ii,jj) = histi(kk)    ! 前回のひずみ
                  delta_str(ii,jj) = str(ii,jj) - str_prev(ii,jj)
                enddo
              enddo
c
c             === M行列の取得（histi配列から） ===
              kk = 24  ! xa値をスキップ（位置16-24）
              do jj=1,3
                do ii=1,3
                  kk = kk+1
                  M_mat(ii,jj) = histi(kk)  ! M = 2μ(n + β/3*I)
                enddo
              enddo
c
c             === M:ε(δu)の計算（テンソル内積） ===
              M_eps = 0.d0
              do jj=1,3
                do ii=1,3
                  M_eps = M_eps + M_mat(ii,jj) * delta_str(ii,jj)
                enddo
              enddo
c
c             === δγの計算（論文式） ===
c             δγ = -N^{-1}(M:ε(δu) + g)
              if(dabs(N_scalar).gt.1.d-12) then
                delta_gamma_inc = -(M_eps + g_val_prev) / N_scalar
              else
                delta_gamma_inc = 0.d0
              endif
c
c             stress_dp_bnに渡すためhisti(1)に保存
              histi(1) = delta_gamma_inc
c
            else
c             === 初回反復：δγ = 0で初期化（BOX 2 Step 1） ===
              histi(1) = 0.d0
            endif
c
            CALL stress_dp_bn(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens, g_vals(ig),
     &                  ierror,  itr , histi )     
          else
            STOP 'Something wrong in phexa8'
          endif
            if(ierror.ne.0) RETURN
c
c ***** Store Deformation Histories (at current state) *****************
            dhist1(:,ig,nel) = ehist(:)
            histi0(:,ig,nel) = histi(:)
            ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
            idep(ig) = idepg
c
            alpeg = ehist(1)
c ***** Compute the Elemental Value (Sum-up for Element ) **************
            vone = vone +vons*detwxyz
            epse = epse +alpeg*detwxyz
            sts(:,:) = sts(:,:) +sig(:,:)*detwxyz
            stn(:,:) = stn(:,:) +str(:,:)*detwxyz
            e_ene = e_ene +e_dns*detwxyz
            p_ene = p_ene +p_dns*detwxyz
c
c ***** 内力ベクトルの計算 *********************************************
c       f_int = ∫ B^T σ dV
c       各節点の各自由度に対して応力の寄与を積分
            do no=1,node
              do ii=1,ndf
                ndof = ndf*(no -1) +ii
                qol = 0.d0
                do jj=1,ndf
                  qol = qol +sig(jj,ii)*dndx(jj,no)
                enddo
                finte(ndof) = finte(ndof) +qol*detwxyz
              enddo
            enddo
c
  120     continue
  110   continue
  100 continue
c
c ***** Compute the Elemental Valume (Volume Average) ******************
      vone = vone/area
      epse = epse/area
c
      tene_e = tene_e +e_ene
      dene_p = dene_p +p_ene
c
      e_ene = e_ene/area
      p_ene = p_ene/area
c
      sts = sts/area ! sts(:,:) = sts(:,:)/area
      stn = stn/area ! sts(:,:) = stn(:,:)/area
c
c ***** Convert the Plastic ene. into Temperature *****
c
c     bbb = aaa*row*ccc
c     if( bbb.eq.0.d0) then
c       dtempe = 0.d0
c     else
c       dtempe = (p_ene*10.d0**9)/(aaa*row*ccc)
c     endif
c     dtemp = dtemp + dtempe*area
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
