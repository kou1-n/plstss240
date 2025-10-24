       subroutine stress_dp_rm(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror )
c
      implicit double precision(a-h,o-z)
c
      dimension prope(20)
c
c     --- common for tracking loading step and element info ---
      common /debug_info/ nel_current, ig_current, lstep_current
      dimension sts(3,3),stn(3,3)
      dimension str(3,3),stry(3,3),
     &          sig(3,3),oun(3,3),sd(3,3),eel(3,3),
     &          plstrg(3,3)
      dimension ctens(3,3,3,3)
      dimension ehist(20)
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
c     write(*,*) 'stress'
c
c ***** Load Deformation Histories *************************************
      alpeg = ehist(1)
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          plstrg(jj,ii) = ehist(kk)!塑性ひずみテンソルε^p
        enddo
      enddo
c     移動硬化は使用しない（等方硬化のみ）
c
c     DP塑性でのtr(plstrg)確認（体積塑性ひずみ）
c      tr_plstrg = plstrg(1,1) + plstrg(2,2) + plstrg(3,3)
c      write(*,'(A,E12.5)') '  DP tr(plstrg) before = ', tr_plstrg
c
c ***** Set Material Properties ****************************************
c    --- elastic parameters
c     --- Young modulus
      yng = prope(1)
c     --- Poisson's ratio
      poi = prope(2)
c     --- Lame const. 'mu' = shear modulus
      vmu = yng/(2.d0*(1.d0 +poi))
c     --- Lame const. 'lamuda'
      vlm = poi*yng/((1.d0 +poi)*(1.d0 -2.d0*poi))
c     --- Bulk modulus 'kappa'
      vkp = yng/(3.d0*(1.d0 -2.d0*poi))
c    --- plastic parameters
      yld = prope(10)
      hk  = prope(11)
      hpa = prope(12)
      hpb = prope(13)
      hpc = prope(14)
      hpd = prope(15)
      phi_dp = prope(16)
      psi_dp = prope(17)
      eta_dp   = 6.d0 * sin(phi_dp) / ( dsqrt(3.d0) * 
     &                                 (3.d0 - sin(phi_dp)) )
      xi_dp    = 6.d0 * cos(phi_dp) / ( dsqrt(3.d0) *
     &                                 (3.d0 - sin(phi_dp)) )
      etabar_dp = 6.d0 * sin(psi_dp) / ( dsqrt(3.d0) *
     &                                  (3.d0 - sin(psi_dp)) )
c
c    --- thermal parameters
      row = prope(3)
      ccc = prope(4)
      aaa = prope(5)
c
c ***** Initialization *************************************************
      deltag = 0.d0
      ctens = 0.d0 ! ctens(:,:,:,:) = 0.d0
c
c   === Deviatoric Stress => Trial Stress ===
c   DP FIX: 平均（体積）ひずみ
      emean  = (str(1,1)+str(2,2)+str(3,3)) / 3.d0
      epmean = (plstrg(1,1)+plstrg(2,2)+plstrg(3,3)) / 3.d0
c
c   体積ひずみ（全ひずみのトレース）
      etrs = str(1,1) +str(2,2) +str(3,3)  !ε_v_n+1^tr = tr(ε_n+1)
c
c   DP FIX: 試行平均応力の計算（塑性体積ひずみを考慮）
c   epv0 = tr(ε^p_n) = 3*epmean (既存変数を活用)
      epv0 = 3.d0 * epmean  ! 塑性体積ひずみ
      p_try = vkp * (etrs - epv0)  ! 試行平均応力 p_try = K*(ε_v - ε_v^p)
c      write(*,'(A,E12.5)') '  DP p_try = K(εv-εvp) = ', p_try
c
c   DP FIX: 試行偏差応力（塑性ひずみの体積成分除外）
      stry = 2.d0*vmu*( (str -emean*DELTA)
     &                 -(plstrg -epmean*DELTA) )
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +stry(ii,jj)**2 !試行のstryで作成したノルム
        enddo
      enddo
      stno = dsqrt(stno)
c
c  === Compute Hardening Function & Yield Function ===
      hard = hpd*(hk*alpeg
     &     +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg)))
c
c     DP FIX: 降伏関数の体積項をp_tryに置換
      ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*p_try
     &      - xi_dp*(yld +hard)
c
c          write(*,'(A,E12.5)') '  ftreg     = ', ftreg
c  ===== PLASTIC CASE
c          determine the Lagrange multiplier by N.R. iteration =====
      if(ftreg.gt.0.d0) then
c
c       === DEBUG: Material parameters and initial state ===
c        write(*,'(A)') '  --- Entering Plastic Correction ---'
c        write(*,'(A,E14.6,A)') '    phi_dp   = ', phi_dp*180.d0/
c     &    3.14159265358979d0, ' deg'
c        write(*,'(A,E14.6,A)') '    psi_dp   = ', psi_dp*180.d0/
c     &    3.14159265358979d0, ' deg'
c        write(*,'(A,E14.6)') '    eta_dp   = ', eta_dp
c        write(*,'(A,E14.6)') '    etabar_dp= ', etabar_dp
c        write(*,'(A,E14.6)') '    xi_dp    = ', xi_dp
c        write(*,'(A,E14.6)') '    G (vmu)  = ', vmu
c        write(*,'(A,E14.6)') '    K (vkp)  = ', vkp
c        write(*,'(A,E14.6)') '    K/G      = ', vkp/vmu
c        write(*,'(A,E14.6)') '    p_try    = ', p_try
c        write(*,'(A,E14.6)') '    ||s_try||= ', stno
c        write(*,'(A,E14.6)') '    alpeg    = ', alpeg
c        write(*,'(A,E14.6)') '    hk       = ', hk
c        write(*,'(A,E14.6)') '    hpa      = ', hpa
c        write(*,'(A,E14.6)') '    hpb      = ', hpb
c       === END DEBUG ===
c
c             write(*,*) 1
        idepg = 1
c     --- initilization ( Box 3.1. step 1 )
        deltag = 0.d0
c 前回の降伏曲面上にいた状態からのスタートなので，deltag=0でよい
        alptmp = alpeg
c
c     --- compute "\Delta gamma" by N.R. iteration
c                       ( Box 3.1. step 2 )
        do 200 it=1,itrmax
c
c         --- K(\alpha^{(n)}_{n+1}) -\sigma_Y
          hrdtmp = hpd*(hk*alptmp
     &            +(hpa -yld)*(1.d0 -dexp(-hpb*alptmp)))
c         --- K'(\alpha^{(n)}_{n+1})
          dhdtmp = hpd*(hk
     &            +hpb*(hpa -yld)*dexp(-hpb*alptmp))
c
c         ---gg(8.117 p363)
c         DP FIX: Newton残差の体積項をp_tryに置換
          gg = dsqrt(1.d0/2.d0)*stno
     &       - vmu*deltag
     &       + eta_dp*(p_try - vkp*etabar_dp*deltag)
     &       - xi_dp*(yld +hrdtmp)
c         ---ggのdeltagによる偏微分
c         CORRECTED: H = ξ*(∂K/∂α)/√3, NOT ξ²*(∂K/∂α)
          sqrt3 = dsqrt(3.d0)
          Dg = -vmu
     &         -eta_dp*vkp*etabar_dp
     &         -xi_dp*dhdtmp/sqrt3
c
          deltag = deltag -gg/Dg
c         --- check convergence
c         CORRECTED: Δα = Δγ/√3, NOT ξ*Δγ
          alptmp = alpeg + deltag/sqrt3
c
c         DP塑性でのNewton反復中のdeltag確認（全反復で出力）
c          write(*,'(A,I3,A,E12.5)') '    Iter[', it, '] deltag = ',
c     &                               deltag
c
c         --- デバッグ用出力（iteration毎に整理）
c         if(it.eq.1) then
c           write(*,*)
c           write(*,'(A)') 
c    &        '========== Plastic Return Mapping Debug Info =========='
c           write(*,'(A,I3,A,I5,A,I3)') 
c    &        '  Loading Step: ', lstep_current,
c    &        '  Element: ', nel_current, 
c    &        '  Gauss Point: ', ig_current
c           write(*,'(A)') 
c    &        '--------------------------------------------------------'
c           write(*,'(A,E12.5)') '  Initial stno     = ', stno
c           write(*,'(A,E12.5)') '  Convergence tol  = ', ctol
c           write(*,'(A,E12.5)') '  vmu (shear mod)  = ', vmu
c           write(*,'(A,E12.5)') '  vkp (bulk mod)   = ', vkp
c           write(*,'(A,E12.5)') '  eta_dp           = ', eta_dp
c           write(*,'(A,E12.5)') '  etabar_dp        = ', etabar_dp
c           write(*,'(A,E12.5)') '  xi_dp            = ', xi_dp
c           write(*,'(A)') 
c    &        '--------------------------------------------------------'
c         endif
c
c          write(*,'(A,I3,A)') '  Iter[', it, ']:'
c          write(*,'(A,E12.5,A,E12.5)')
c     &      '    gg=', gg, '  Dg=', Dg
c          write(*,'(A,E12.5,A,E12.5)')
c     &      '    gg/Dg=', gg/Dg, '  |gg/Dg|=', dabs(gg/Dg)
c          write(*,'(A,E12.5,A,E12.5)')
c     &      '    deltag=', deltag, '  alptmp=', alptmp
c          write(*,'(A,E12.5,A,E12.5)')
c     &      '    hrdtmp=', hrdtmp, '  dhdtmp=', dhdtmp
c          write(*,'(A,E12.5,A,E12.5)')
c     &      '    H_correct=', xi_dp*dhdtmp/sqrt3,
c     &      '  H_old(wrong)=', xi_dp*xi_dp*dhdtmp
c
          if( dabs(gg/Dg).lt.ctol) then
c           write(*,'(A,I3,A)') 
c    &        '  --> Converged at iteration ', it, ' !'
c           write(*,'(A)') 
c    &        '======================================================='
c           write(*,*)
            GOTO 210
          endif
  200   continue
c
c     --- error section: Failed to compute "Delta gamma"
c       write(*,*)
c       write(*,'(A)') 
c    &    '********************************************************'
c       write(*,'(A)') 
c    &    '**** ERROR: Plastic Return Mapping Failed (ierror=17) *'
c       write(*,'(A)') 
c    &    '********************************************************'
c       write(*,'(A,I3,A,I5,A,I3)') 
c    &    '  Loading Step: ', lstep_current,
c    &    '  Element: ', nel_current, 
c    &    '  Gauss Point: ', ig_current
c       write(*,'(A,I3)') '  Total iterations attempted: ', itrmax
c       write(*,'(A)') '  Final convergence status:'
c       write(*,'(A,E12.5)') '    Final gg       = ', gg
c       write(*,'(A,E12.5)') '    Final Dg       = ', Dg
c       if(dabs(Dg).gt.1.0d-20) then
c         write(*,'(A,E12.5)') '    Final gg/Dg    = ', gg/Dg
c         write(*,'(A,E12.5)') '    Final |gg/Dg|  = ', dabs(gg/Dg)
c       else
c         write(*,'(A)') '    WARNING: Dg is nearly zero!'
c       endif
c       write(*,'(A,E12.5)') '    Required tol   = ', ctol
c       write(*,'(A,E12.5)') '    Final deltag   = ', deltag
c       write(*,'(A,E12.5)') '    Final alptmp   = ', alptmp
c       write(*,'(A)') 
c    &    '********************************************************'
c       write(*,*)
        ierror = 17
        RETURN
c
  210   CONTINUE
c
c       === DEBUG: Return mapping converged ===
c        write(*,'(A)') '  --- Return Mapping Converged ---'
c        write(*,'(A,E14.6)') '    Final deltag = ', deltag
c        write(*,'(A,E14.6)') '    sqrt(2)*G*deltag = ',
c     &    dsqrt(2.d0)*vmu*deltag
c        write(*,'(A,E14.6)') '    ||s_try||    = ', stno
c        write(*,'(A,E14.6)') '    Ratio sqrt(2)*G*deltag/||s_try|| = ',
c     &    (dsqrt(2.d0)*vmu*deltag)/stno
c       === END DEBUG ===
c
c     --- update equivalent plastic strain "alpha"
c       CORRECTED: Δα = Δγ/√3, NOT ξ*Δγ
        alptmp = alpeg
        alpeg = alpeg +deltag/sqrt3
c
c       DP塑性でのΔγ（塑性乗数）確認
c        write(*,'(A,E12.5)') '  DP deltag (Δγ)       = ', deltag
c        write(*,'(A,E12.5)') '  DP alpeg (α_n+1)     = ', alpeg
c
c     --- outward unit normal vector in stress space
        oun(:,:) = stry(:,:)/stno
c     ounもNRによって求められたΔγの値に依存するのでは...？？stryとstnoでは違和感がある
c     →試行状態の流れベクトルと偏差方向の流れベクトルは同じ方向を向くので問題ない
c
c
c       DP塑性でのtr(oun)確認（偏差方向）
c        tr_oun = oun(1,1) + oun(2,2) + oun(3,3)
c        write(*,'(A,E12.5)') '  DP tr(oun)           = ', tr_oun
c
c     --- compute the plastic strain et.al.
c       etrs already computed at the beginning
c    8.109に基づく偏差ひずみのアップデート        
        plstrg(:,:) = plstrg(:,:)
     &              + deltag*(
     &              + oun(:,:)*dsqrt(1.d0/2.d0)
     &              + etabar_dp*DELTA(:,:)/3.d0
     &              ) !ok
c       DP塑性でのtr(plstrg)確認（更新後、体積変化あり）
c        tr_plstrg = plstrg(1,1) + plstrg(2,2) + plstrg(3,3)
c        write(*,'(A,E12.5)') '  DP tr(plstrg) after  = ', tr_plstrg
c       DP塑性の体積ひずみ増分
c        tr_delta = deltag * etabar_dp
c        write(*,'(A,E12.5)') '  DP Δε_v (=Δγ*etabar) = ', tr_delta


c       DP FIX: 応力更新（塑性時）の体積項をp_tryに置換
        sig(:,:) = stry(:,:) -dsqrt(2.d0)*vmu*deltag*oun(:,:)
     &            +(p_try - vkp*etabar_dp*deltag)*DELTA(:,:)
c       p_n+1 = p_try - K*etabar*Δγ (体積応力更新)


c     --- update consititutive tensor: "ctens"
        dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg))
c       CORRECTED: H = ξ*(∂K/∂α)/√3, NOT ξ²*(∂K/∂α)
        A = 1.d0/(vmu + vkp*etabar_dp*eta_dp
     &          + xi_dp*dhard/sqrt3)

c       === DEBUG: Check consistency tangent coefficients ===
c        write(*,'(A)') '  --- Consistent Tangent Diagnostics ---'
c        write(*,'(A,E14.6)') '    dhard            = ', dhard
c        write(*,'(A,E14.6)') '    A (corrected)    = ', A
c        write(*,'(A,E14.6)') '    Denominator (corrected) = ',
c          &        vmu + vkp*etabar_dp*eta_dp + xi_dp*dhard/sqrt3
c        write(*,'(A,E14.6)') '    H_correct = ξ*(∂K/∂α)/√3 = ',
c          &        xi_dp*dhard/sqrt3
c        write(*,'(A,E14.6)') '    H_old(wrong) = ξ²*(∂K/∂α) = ',
c          &        xi_dp*xi_dp*dhard
c        write(*,'(A,E14.6)') '    K*eta*etabar     = ',
c          &        vkp*eta_dp*etabar_dp
c        write(*,'(A,E14.6)') '    K*eta*etabar*A   = ',
c          &        vkp*eta_dp*etabar_dp*A

        theta = 1.d0 -(dsqrt(2.d0)*vmu*deltag)/stno !ok
        thetab= (dsqrt(2.d0)*vmu*deltag)/stno - vmu*A !ok

c        write(*,'(A,E14.6)') '    theta            = ', theta
c        write(*,'(A,E14.6)') '    thetab           = ', thetab
c        write(*,'(A,E14.6)') '    bulk_coef (1-K*eta*etabar*A) = ',
c          &        1.d0-vkp*eta_dp*etabar_dp*A
c
        if(theta .lt. 0.d0) then
c          write(*,'(A)') '    *** WARNING: theta < 0 ***'
        endif
        if((1.d0-vkp*eta_dp*etabar_dp*A) .lt. 0.01d0) then
c          write(*,'(A)') '    *** WARNING: bulk stiffness too small ***'
        endif
c       === END DEBUG ===
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            =  2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
     &
     &              +2.d0*vmu*thetab*oun(ii,jj)*oun(kk,ll)
     &
     &              - dsqrt(2.d0)*vmu*A*vkp
     &              *(eta_dp*oun(ii,jj)*DELTA(kk,ll)
     &              + etabar_dp*DELTA(ii,jj)*oun(kk,ll)) 
     & 
     &              + vkp*(1.d0-vkp*eta_dp*etabar_dp*A)
     &                         *DELTA(ii,jj)*DELTA(kk,ll)
     
              enddo
            enddo
          enddo
        enddo
c
c  ===== ELASTIC CASE:
c          Identify the Trial Stress as Actual One =====
      else
        idepg = 0
c
c       DP塑性：弾性ケースのdeltag確認（ゼロのはず）
c        write(*,'(A,E12.5)') '  DP deltag (elastic)  = ', deltag
c
c       DP FIX: 応力更新（弾性時）の体積項をp_tryに置換
c       塑性体積ひずみを考慮してp_try = K*(ε_v - ε_v^p)を使用
        sig(:,:) = stry(:,:) +p_try*DELTA(:,:)
c       write(*,*) etrs,ctol
c       write(*,*) vkp*etrs,'aaa'
c
c     --- donot update "ctens" ---
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
c               ctens(ii,jj,kk,ll) = vlm*DELTA(ii,jj)*DELTA(kk,ll)
c    &                              +2.d0*vmu*FIT(ii,jj,kk,ll)
                ctens(ii,jj,kk,ll)
     &            = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &              +2.d0*vmu*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
              enddo
            enddo
          enddo
        enddo
c
      endif
c
c   === Compute von Mieses Stress  ===
      smean = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
c     write(*,*) smean
c     write(*,*) sig
c
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
c     --- separate the strain into the elastic & plastic
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
c     移動硬化変数は保存しない（使用しないため）
c
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
