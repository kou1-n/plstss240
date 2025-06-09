      subroutine stress(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
c    &                   prope,!plstrg,betaeg, alpeg,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror )
c
      implicit double precision(a-h,o-z)
c     External hardening functions
      external H_iso, dH_iso_dk, H_kin, dH_kin_dk
      double precision H_iso, dH_iso_dk, H_kin, dH_kin_dk
c
      dimension prope(20)
      dimension sts(3,3),stn(3,3) !使われていない
      dimension str(3,3),stry(3,3),seta(3,3),
     &          sig(3,3),oun(3,3),sd(3,3),eel(3,3),
     &          plstrg(3,3),betaeg(3,3)
      dimension ctens(3,3,3,3)
      dimension ehist(20)
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)


c ***** ここで A を定義する (コンシステントタンジェントで用いるスカラー) *****
      double precision A
c ***** ここで N を定義する (非関連流れ則で用いる2階のテンソルN(;,;)) *****
      double precision N(3,3)
c **********************************************************************
c     write(*,*) 'stress'
c
c ***** Load Deformation Histories *************************************
      alpeg = ehist(1)
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          plstrg(jj,ii) = ehist(kk) !塑性ひずみテンソル（Plastic Strain Tensor）
        enddo
      enddo
      do ii=1,3
        do jj=1,3
          kk = kk+1
          betaeg(jj,ii) = 0.d0 !移動硬化は考えないのでehist(kk)=0
        enddo
      enddo
c
c ***** Set Material Properties ****************************************
c    --- elastic parameters
c     --- Young modulus
      yng = prope(1)
c     --- Poisson's ratio
      poi = prope(2)
c       --- Lame const. 'mu' = shear modulus
      vmu = yng/(2.d0*(1.d0 +poi))
c       --- Lame const. 'lamuda'
      vlm = poi*yng/((1.d0 +poi)*(1.d0 -2.d0*poi))
c       --- Bulk modulus 'kappa'
      vkp = yng/(3.d0*(1.d0 -2.d0*poi))
c    --- plastic parameters
      yld = prope(10) !c0の値に対応
      hk  = prope(11) !H_barの値に対応
      hpa = prope(12) !H_inf_barに対応
      hpb = prope(13) !xi_barの減衰率に対応
c     hpc = prope(14)
      hpd = prope(15) !:等方硬化(1)と移動硬化(0)関数の割引を表すtheta[0,1]
      phi_dp = prope(16) !FRICTION_ANGLE
      psi_dp = prope(17) !DILATANCY_ANGLE
c Drucker-Prager の η, ξ, および 非関連の場合の η̄ の外側一致の式に合わせて再定義
      eta_dp=  sin(phi_dp) / ( dsqrt(3.d0) * (3.d0 - sin(phi_dp)) )
      xi_dp = 6.d0 * cos(phi_dp) / ( dsqrt(3.d0) 
     &                              * (3.d0 - sin(phi_dp)) )
      etabar_dp = 6.d0 * sin(psi_dp) / ( dsqrt(3.d0)  
     &                                          * (3.d0 - sin(psi_dp)) )


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
c   === Deviatoric Stress => Trial Stress ===BOX 8.8 step1
      emean = (str(1,1) +str(2,2) +str(3,3))/3.d0 !スカラー
      ptry = vkp*(str(1,1) +str(2,2) +str(3,3)) !試行静水圧成分 K*epsilon_v(tr(epsilon))
c
      stry = 2.d0*vmu*(str -emean*DELTA -plstrg)
      !stry(:,:) = 2.d0*vmu*(str(:,:) -emean*DELTA(:,:) -plstrg(:,:)) !試行偏差応力
      seta = stry -betaeg !試行相当応力　＝　試行偏差応力 - 背応力
      !seta(:,:) = stry(:,:) -betaeg(:,:)
c
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +seta(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)
c
c  === Compute Hardening Function & Yield Function === BOX 8.8 step2
      hard = hpd* hk*alpeg
     &     +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg))
c
      ftreg =dsqrt(1.d0/2.d0)*stno 
     &      +eta_dp*ptry 
     &      -xi_dp*(yld +hard)
c      write(*,*) 'ftreg = ', ftreg
c
c  ===== PLASTIC CASE 
c          determine the Lagrange multiplier by N.R. iteration =====
      if(ftreg.gt.0.d0) then
c             write(*,*) 1
        idepg = 1
c        print *,"PLASTIC"
c     --- initilization ( Box 3.1. step 1 ) box 8.9 step 1
        deltag = 0.d0
        alptmp = alpeg
c
c     --- compute "\Delta gamma" by N.R. iteration
c                       ( Box 3.1. step 2 )
        do 200 it=1,itrmax
          alpd = alptmp -alpeg!alpdは"alpha"の増分（= α^{(n+1)} - α^{(n)}）
c
c         --- K(\alpha^{(n)}_{n+1}) -\sigma_Y
          hrdtmp = H_iso(alptmp, yld, hk, hpa, hpb) - yld
c         --- K'(\alpha^{(n)}_{n+1})
          dhdtmp = dH_iso_dk(alptmp, yld, hk, hpa, hpb)
c
c         --- H(\alpha^{(n)}_{n+1}) -H(\alpha_{n})
          tmpkrd = H_kin(alptmp, hk, hpd) - H_kin(alpeg, hk, hpd)
          ! tmpkrd は運動硬化項 (H) の増分:
          !   H(α^{(n+1)}) - H(α^n) = (1 - hpd)*hk*(α^{(n+1)} - α^n)
c          
c         --- H'(\alpha^{(n)}_{n+1}) -H(\alpha_{n})
          dkdtmp = dH_kin_dk(hk, hpd)
c         ! dkdtmp は dH/dα：線形のため定数 (1 - hpd)*hk
c
c
          gg = dsqrt(1.d0/2.d0)*stno
     &         -vmu*deltag
     &         +eta_dp*ptry
     &         -eta_dp*vkp*etabar_dp*deltag
     &         -xi_dp*(yld +hrdtmp +tmpkrd) 
c             
          Dg = -vmu 
     &         -vkp*etabar_dp*eta_dp
     &         -xi_dp*xi_dp*(dhdtmp +dkdtmp)
c
          deltag = deltag -gg/Dg
c          print '("Update deltag=",e25.15)',deltag
          alptmp = alpeg + xi_dp*deltag !第二項にsqrt(2/3)をかける必要がある?
c
c
          if( dabs(gg/Dg).lt.ctol) then
            GOTO 210
          endif
  200   continue
c
c     --- error section: Failed to compute "Delta gamma"
        ierror = 17
        RETURN
c
  210   CONTINUE
c     --- update equivalent plastic strain "alpha"
        alptmp = alpeg
        alpeg = alpeg +xi_dp*deltag
c
c     --- outward unit normal vector in stress space
C     
        oun(:,:) = seta(:,:)/stno

        N(:,:)   = seta(:,:)/ (dsqrt(2.d0)*stno) 
     &           + etabar_dp/ 3.d0 * DELTA(:,:)
c
c     --- compute the plastic strain et.al.
        alpd = alpeg -alptmp
c
        tmpkrd = (1.d0 -hpd)* hk*alpd
        etrs = str(1,1) +str(2,2) +str(3,3)
c===Neto  BOX7.5 p290 step
        betaeg(:,:) = betaeg(:,:) +deltag*dsqrt(2.d0)*tmpkrd*N(:,:)

        plstrg(:,:) = plstrg(:,:) +deltag*N(:,:)
        sig(:,:) = stry(:,:) -dsqrt(2.d0)*vmu*deltag*oun(:,:)
     &            +(ptry-vkp*etabar_dp*deltag)*DELTA(:,:)
c
c     --- update consititutive tensor: "ctens"
        dhard = dH_iso_dk(alpeg, yld, hk, hpa, hpb)
        dkard = dH_kin_dk(hk, hpd)
c    &        +(hpa -yld)*(1.d0 -dexp(-hpb*alpha  )) コメントアウトして移動硬化を線形化している
            A = 1.d0 / ( vmu + vkp*etabar_dp*eta_dp 
     &               + xi_dp*xi_dp*(dhdtmp + dkdtmp) )

c            D(:,:) = seta(:,:) / stno
        theta = 1.d0 -(dsqrt(2.d0)*vmu*deltag)/stno
        thetab= (dsqrt(2.d0)*vmu*deltag)/stno -vmu*A
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            = vkp*(1.d0 -vkp*eta_dp*etabar_dp*A)
     &                 *DELTA(ii,jj)*DELTA(kk,ll)
     &              +2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
     &              +2.d0*vmu*thetab*oun(ii,jj)*oun(kk,ll)
     &              -dsqrt(2.d0)*vmu*A*vkp*
     &                         ( eta_dp*oun(ii,jj)*DELTA(kk,ll)
     &                        + etabar_dp*DELTA(ii,jj)*oun(kk,ll) )
              enddo
            enddo
          enddo
        enddo
c
c  ===== ELASTIC CASE:
c          Identify the Trial Stress as Actual One =====
      else
        idepg = 0
c        print *,"ELASTIC"
c
        etrs = str(1,1) +str(2,2) +str(3,3)
c       if(dabs(etrs).le.1.0d-16) etrs = 0.d0
        sig(:,:) = stry(:,:) +vkp*etrs*DELTA(:,:)
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
      do ii=1,3
        do jj=1,3
          kk = kk+1
          ehist(kk) = betaeg(jj,ii)
        enddo
      enddo
c      write(*,*) plstrg
c
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
