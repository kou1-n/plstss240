      subroutine stress(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
c    &                   prope,!plstrg,betaeg, alpeg,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror )
c
      implicit double precision(a-h,o-z)
c
      dimension prope(20)
      dimension sts(3,3),stn(3,3)
      dimension str(3,3),stry(3,3),seta(3,3),
     &          sig(3,3),oun(3,3),sd(3,3),eel(3,3),
     &          plstrg(3,3),betaeg(3,3)
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
          plstrg(jj,ii) = ehist(kk)
        enddo
      enddo
      do ii=1,3
        do jj=1,3
          kk = kk+1
          betaeg(jj,ii) = ehist(kk)
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
      yld = prope(10)
      hk  = prope(11)
      hpa = prope(12)
      hpb = prope(13)
c     hpc = prope(14)
      hpd = prope(15)
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
      emean = (str(1,1) +str(2,2) +str(3,3))/3.d0
c
      stry = 2.d0*vmu*(str -emean*DELTA -plstrg)
      !stry(:,:) = 2.d0*vmu*(str(:,:) -emean*DELTA(:,:) -plstrg(:,:))
      seta = stry -betaeg
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
c  === Compute Hardening Function & Yield Function ===
      hard = hpd* hk*alpeg
     &     +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg))
c
      ftreg = stno -dsqrt(2.d0/3.d0)*(yld +hard)
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
          alpd = alptmp -alpeg
c
c         --- K(\alpha^{(n)}_{n+1}) -\sigma_Y
          hrdtmp = hpd*hk*alptmp
     &            +(hpa -yld)*(1.d0 -dexp(-hpb*alptmp))
c         --- K'(\alpha^{(n)}_{n+1})
          dhdtmp = hpd*hk
     &            +hpb*(hpa -yld)*dexp(-hpb*alptmp)
c
c         --- H(\alpha^{(n)}_{n+1}) -H(\alpha_{n})
          tmpkrd = (1.d0 -hpd)* hk*alpd
c
c         --- H'(\alpha^{(n)}_{n+1}) -H(\alpha_{n})
          dkdtmp = (1.d0 -hpd)* hk
c
          gg = -dsqrt(2.d0/3.d0)*(yld +hrdtmp +tmpkrd)
     &         +stno -2.d0*vmu*deltag
          Dg = -2.d0*vmu -(2.d0/3.d0)*(dhdtmp+dkdtmp)
c
          deltag = deltag -gg/Dg
          alptmp = alpeg +dsqrt(2.d0/3.d0)*deltag
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
        alpeg = alpeg +dsqrt(2.d0/3.d0)*deltag
c
c     --- outward unit normal vector in stress space
        oun(:,:) = seta(:,:)/stno
c
c     --- compute the plastic strain et.al.
        alpd = alpeg -alptmp
c
        tmpkrd = (1.d0 -hpd)* hk*alpd
        etrs = str(1,1) +str(2,2) +str(3,3)
c
        betaeg(:,:) = betaeg(:,:) +dsqrt(2.d0/3.d0)*tmpkrd*oun(:,:)
        plstrg(:,:) = plstrg(:,:) +deltag*oun(:,:)
        sig(:,:) = stry(:,:) -2.d0*vmu*deltag*oun(:,:)
     &                                 +vkp*etrs*DELTA(:,:)
c
c     --- update consititutive tensor: "ctens"
        dhard = hpd* hk +hpb*(hpa -yld)*dexp(-hpb*alpeg)
        dkard = (1.d0 -hpd)* hk
c    &        +(hpa -yld)*(1.d0 -dexp(-hpb*alpha  ))

        theta = 1.d0 -(2.d0*vmu*deltag)/stno
        thetab= 1.d0/(1.d0 +(dhard+dkard)/(3.d0*vmu)) -(1.d0 -theta)
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &              +2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
     &              -2.d0*vmu*thetab*oun(ii,jj)*oun(kk,ll)
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
c
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
