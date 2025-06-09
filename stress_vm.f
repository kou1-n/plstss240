      subroutine stress_vm(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror )
      implicit double precision(a-h,o-z)

      dimension prope(20)
      dimension str(3,3), sig(3,3), stry(3,3), seta(3,3)
      dimension plstrg(3,3), betaeg(3,3), oun(3,3)
      dimension sd(3,3), eel(3,3)
      dimension ctens(3,3,3,3)
      dimension ehist(20)

      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)

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

c ***** Set Material Properties ****************************************
      yng = prope(1)
      poi = prope(2)
      vmu = yng/(2.d0*(1.d0+poi))
      vlm = poi*yng/((1.d0+poi)*(1.d0-2.d0*poi))
      vkp = yng/(3.d0*(1.d0-2.d0*poi))
      yld = prope(10)
      hk  = prope(11)
      hpd = prope(15)

c ***** Initialization *************************************************
      ctens = 0.d0
      deltag = 0.d0

c === Deviatoric Stress => Trial Stress ===
      emean = (str(1,1)+str(2,2)+str(3,3))/3.d0
      stry(:,:) = 2.d0*vmu*(str(:,:) -emean*DELTA(:,:) -plstrg(:,:))
      seta(:,:) = stry(:,:) -betaeg(:,:)
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +seta(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)

c === Yield Function ===
      hard = hpd*hk*alpeg
      ftreg = stno -dsqrt(2.d0/3.d0)*(yld +hard)

c ===== PLASTIC CASE =====
      if(ftreg.gt.0.d0) then
        idepg = 1
        deltag = ftreg /(2.d0*vmu +(2.d0/3.d0)*hk)
        alpeg = alpeg +dsqrt(2.d0/3.d0)*deltag
        oun(:,:) = seta(:,:)/stno
      betaeg(:,:) = betaeg(:,:)
     &              +(2.d0/3.d0)*(1.d0-hpd)*hk*deltag*oun(:,:)
        plstrg(:,:) = plstrg(:,:) +deltag*oun(:,:)
        etrs = str(1,1)+str(2,2)+str(3,3)
        sig(:,:) = stry(:,:) -2.d0*vmu*deltag*oun(:,:)
     &            +vkp*etrs*DELTA(:,:)
        theta = 1.d0 -(2.d0*vmu*deltag)/stno
        thetab= 1.d0/(1.d0 +hk/(3.d0*vmu)) -(1.d0 -theta)
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
      else
        idepg = 0
        etrs = str(1,1)+str(2,2)+str(3,3)
        sig(:,:) = stry(:,:) +vkp*etrs*DELTA(:,:)
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
      endif

c === von Mises Stress ===
      smean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
      sd(:,:) = sig(:,:) -smean*DELTA(:,:)
      vons = 0.d0
      do jj=1,3
        do ii=1,3
          vons = vons +sd(ii,jj)*sd(ii,jj)
        enddo
      enddo
      vons = dsqrt(1.5d0*vons)

c === Energy Density ===
      eel(:,:) = str(:,:) -plstrg(:,:)
      e_dns = 0.d0
      do jj=1,3
        do ii=1,3
          e_dns = e_dns +eel(ii,jj)*sig(ii,jj)
        enddo
      enddo
      e_dns = 0.5d0*e_dns
      p_dns = 0.d0
      do jj=1,3
        do ii=1,3
          p_dns = p_dns +deltag*oun(ii,jj)*sig(ii,jj)
        enddo
      enddo

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

      RETURN
      END
