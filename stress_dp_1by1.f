       subroutine stress_dp_1by1(itrmax, idepg,
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
c   emean:静水圧成分（スカラー）
      emean = (str(1,1) +str(2,2) +str(3,3))/3.d0
c
      stry = 2.d0*vmu*(str -emean*DELTA -plstrg) !p366 s_n+1^tr
c
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
      ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*emean 
     &      - xi_dp*(yld +hard)
c
c  ===== PLASTIC CASE
c          determine the Lagrange multiplier by N.R. iteration =====
      if(ftreg.gt.0.d0) then
c             write(*,*) 1
        idepg = 1
c     --- initilization ( Box 3.1. step 1 )
        deltag = 0.d0
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
c         ---ggに降伏関数の前のステップの値を入れる
          gg = dsqrt(1.d0/2.d0)*stno + eta_dp*emean
     &      - xi_dp*(yld +hrdtmp)
c         ---Dgに降伏関数の微分の前のステップの値を入れる
          Dg = -vmu 
     &         -vkp*eta_dp*etabar_dp
     &         -xi_dp*xi_dp*dhdtmp
c
          deltag = deltag -gg/Dg
c         --- check convergence
          alptmp = alpeg + xi_dp*deltag
c
c         --- デバッグ用出力（iteration毎に整理）
          if(it.eq.1) then
            write(*,*)
            write(*,'(A)') 
     &        '========== Plastic Return Mapping Debug Info =========='
            write(*,'(A,I3,A,I5,A,I3)') 
     &        '  Loading Step: ', lstep_current,
     &        '  Element: ', nel_current, 
     &        '  Gauss Point: ', ig_current
            write(*,'(A)') 
     &        '--------------------------------------------------------'
            write(*,'(A,E12.5)') '  Initial stno     = ', stno
            write(*,'(A,E12.5)') '  Convergence tol  = ', ctol
            write(*,'(A,E12.5)') '  vmu (shear mod)  = ', vmu
            write(*,'(A,E12.5)') '  vkp (bulk mod)   = ', vkp
            write(*,'(A,E12.5)') '  eta_dp           = ', eta_dp
            write(*,'(A,E12.5)') '  etabar_dp        = ', etabar_dp
            write(*,'(A,E12.5)') '  xi_dp            = ', xi_dp
            write(*,'(A)') 
     &        '--------------------------------------------------------'
          endif
c
          write(*,'(A,I3,A)') '  Iter[', it, ']:'
          write(*,'(A,E12.5,A,E12.5)') 
     &      '    gg=', gg, '  Dg=', Dg
          write(*,'(A,E12.5,A,E12.5)') 
     &      '    gg/Dg=', gg/Dg, '  |gg/Dg|=', dabs(gg/Dg)
          write(*,'(A,E12.5,A,E12.5)') 
     &      '    deltag=', deltag, '  alptmp=', alptmp
          write(*,'(A,E12.5,A,E12.5)') 
     &      '    hrdtmp=', hrdtmp, '  dhdtmp=', dhdtmp
c
          if( dabs(gg/Dg).lt.ctol) then
            write(*,'(A,I3,A)') 
     &        '  --> Converged at iteration ', it, ' !'
            write(*,'(A)') 
     &        '========================================================'
            write(*,*)
            GOTO 210
          endif
  200   continue
c
c     --- error section: Failed to compute "Delta gamma"
        write(*,*)
        write(*,'(A)') 
     &    '********************************************************'
        write(*,'(A)') 
     &    '**** ERROR: Plastic Return Mapping Failed (ierror=17) *'
        write(*,'(A)') 
     &    '********************************************************'
        write(*,'(A,I3,A,I5,A,I3)') 
     &    '  Loading Step: ', lstep_current,
     &    '  Element: ', nel_current, 
     &    '  Gauss Point: ', ig_current
        write(*,'(A,I3)') '  Total iterations attempted: ', itrmax
        write(*,'(A)') '  Final convergence status:'
        write(*,'(A,E12.5)') '    Final gg       = ', gg
        write(*,'(A,E12.5)') '    Final Dg       = ', Dg
        if(dabs(Dg).gt.1.0d-20) then
          write(*,'(A,E12.5)') '    Final gg/Dg    = ', gg/Dg
          write(*,'(A,E12.5)') '    Final |gg/Dg|  = ', dabs(gg/Dg)
        else
          write(*,'(A)') '    WARNING: Dg is nearly zero!'
        endif
        write(*,'(A,E12.5)') '    Required tol   = ', ctol
        write(*,'(A,E12.5)') '    Final deltag   = ', deltag
        write(*,'(A,E12.5)') '    Final alptmp   = ', alptmp
        write(*,'(A)') 
     &    '********************************************************'
        write(*,*)
        ierror = 17
        RETURN
c
  210   CONTINUE
c     --- update equivalent plastic strain "alpha"
        alptmp = alpeg
        alpeg = alpeg +xi_dp*deltag
c
c     --- outward unit normal vector in stress space
        oun(:,:) = stry(:,:)/stno
c
c     --- compute the plastic strain et.al.
        etrs = str(1,1) +str(2,2) +str(3,3)
c    8.109に基づく偏差ひずみのアップデート        
        plstrg(:,:) = plstrg(:,:) 
     &              + deltag*(
     &              + oun(:,:)*dsqrt(1.d0/2.d0)     
     &              + etabar_dp*DELTA(:,:)/3.d0
     &              ) !ok


        sig(:,:) = stry(:,:) -dsqrt(2.d0)*vmu*deltag*oun(:,:) !ok
     &            +((etrs/3.d0) - vkp*etabar_dp*deltag)*DELTA(:,:)!ok
c


c     --- update consititutive tensor: "ctens"
        dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg))
        A = 1.d0/(vmu + vkp*etabar_dp*eta_dp 
     &          + xi_dp*xi_dp*dhard)

        theta = 1.d0 -(dsqrt(2.d0)*vmu*deltag)/stno !ok
        thetab= (dsqrt(2.d0)*vmu*deltag)/stno - vmu*A !ok
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            = +2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
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
c 以上は大丈夫0715_1745
c  ===== ELASTIC CASE:
c          Identify the Trial Stress as Actual One =====
      else
        idepg = 0
c
        etrs = str(1,1) +str(2,2) +str(3,3)
c       if(dabs(etrs).le.1.0d-16) etrs = 0.d0
        sig(:,:) = stry(:,:) +vkp*etrs*DELTA(:,:)
c       write(*,*) etrs,ctol
c       write(*,*) vkp*etrs,'aaa'
c
c     --- donot update "ctens" ---
c     --- ここctens,von-misesと同じにしているけどいいのか・・・
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
