      subroutine stress_dp_rm(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror )
c
      implicit double precision(a-h,o-z)
c
      dimension prope(20)
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
c ***** Load Deformation Histories *************************************
      alpeg = ehist(1)
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          plstrg(jj,ii) = ehist(kk)
        enddo
      enddo
c
c ***** Set Material Properties ****************************************
c    --- elastic parameters
      yng = prope(1)
      poi = prope(2)
      vmu = yng/(2.d0*(1.d0 +poi))
      vlm = poi*yng/((1.d0 +poi)*(1.d0 -2.d0*poi))
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
c     Drucker–Pragerの係数
      eta_dp    = 2.d0*sin(phi_dp)/(3.d0 - sin(phi_dp))
      xi_dp     = 3.d0*cos(phi_dp)/(3.d0 - sin(phi_dp))
      etabar_dp = 2.d0*sin(psi_dp)/(3.d0 - sin(psi_dp))
c
c ***** Initialization *************************************************
      deltag = 0.d0
      ctens  = 0.d0
c
c   === Trial quantities ==============================================
      emean  = (str(1,1)+str(2,2)+str(3,3)) / 3.d0
      epmean = (plstrg(1,1)+plstrg(2,2)+plstrg(3,3)) / 3.d0
      etrs   = str(1,1) +str(2,2) +str(3,3)
      epv0   = 3.d0*epmean
      p_try  = vkp * (etrs - epv0)
c     deviatoric trial stress
      stry = 2.d0*vmu*( (str -emean*DELTA)
     &                 -(plstrg -epmean*DELTA) )
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno + stry(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)
c
c   === Hardening & yield function ====================================
      hard = hpd*(hk*alpeg + (hpa -yld)*(1.d0 -dexp(-hpb*alpeg)))
      ftreg = dsqrt(3.d0/2.d0)*stno + eta_dp*p_try
     &      - xi_dp*(yld +hard)
c
c  ===== PLASTIC CASE ==================================================
      if(ftreg.gt.0.d0) then
        idepg = 1
c
c       --- ここで「関連流れ」に基づく dα/dγ を先に定義 ---
c           N_f = d f / d σ = sqrt(3/2) n + (eta/3) I
c           ||N_f||^2 = 3/2 + (eta^2)/3
        nf_norm  = dsqrt(1.5d0 + (eta_dp*eta_dp)/3.d0)
        dalp_dgam= dsqrt(2.d0/3.d0) * nf_norm   ! dα/dγ（定数）
c
        deltag = 0.d0
        alptmp = alpeg
c
c       --- Newton iteration for Δγ -----------------------------------
        do 200 it=1,itrmax
c         hard(alptmp) and derivative
          hrdtmp = hpd*( hk*alptmp
     &            +(hpa -yld)*(1.d0 -dexp(-hpb*alptmp)) )
          dhdtmp = hpd*( hk
     &            +hpb*(hpa -yld)*dexp(-hpb*alptmp) )
c
c         residual g(Δγ)=0  （q形式）
          gg = dsqrt(3.d0/2.d0)*stno
     &       - 3.d0*vmu*deltag
     &       + eta_dp*(p_try - vkp*etabar_dp*deltag)
     &       - xi_dp*(yld + hrdtmp)
c
c         tangent d g / d(Δγ)
          Dg = -3.d0*vmu
     &         - eta_dp*vkp*etabar_dp
     &         - xi_dp*dhdtmp*dalp_dgam
c
          deltag = deltag - gg/Dg
c         α^(trial) の更新は関連流れの規準に従う
          alptmp = alpeg + dalp_dgam*deltag
          if(dabs(gg/Dg).lt.ctol) goto 210
  200   continue
        ierror = 17
        return
c
  210   continue
c       --- 確定更新 ---------------------------------------------------
c       有効塑性ひずみ（関連流れに基づく）
        alpeg = alpeg + dalp_dgam*deltag
c
c       偏差方向単位テンソル
        oun(:,:) = stry(:,:)/stno
c
c       塑性ひずみ（非関連流れ g：ψ を使用）
        plstrg(:,:) = plstrg(:,:)
     &              + deltag*( dsqrt(3.d0/2.d0)*oun(:,:)
     &              + (etabar_dp/3.d0)*DELTA(:,:) )
c
c       応力更新
        sig(:,:) = stry(:,:) - dsqrt(6.d0)*vmu*deltag*oun(:,:)
     &            + (p_try - vkp*etabar_dp*deltag)*DELTA(:,:)
c
c       一貫接線（A の分母に dα/dγ を反映）
        dhard = hpd*( hk + hpb*(hpa -yld)*dexp(-hpb*alpeg) )
        A = 1.d0 / ( 3.d0*vmu + vkp*etabar_dp*eta_dp
     &               + xi_dp*dhard*dalp_dgam )
c
        theta = 1.d0 - (dsqrt(6.d0)*vmu*deltag)/stno
        thetab= (dsqrt(6.d0)*vmu*deltag)/stno - 3.d0*vmu*A
c
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            =  2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
     &              +2.d0*vmu*thetab*oun(ii,jj)*oun(kk,ll)
     &              - dsqrt(6.d0)*vmu*A*vkp
     &              *(eta_dp*oun(ii,jj)*DELTA(kk,ll)
     &              + etabar_dp*DELTA(ii,jj)*oun(kk,ll))
     &              + vkp*(1.d0 - vkp*eta_dp*etabar_dp*A)
     &                         *DELTA(ii,jj)*DELTA(kk,ll)
              enddo
            enddo
          enddo
        enddo
c
c  ===== ELASTIC CASE ==================================================
      else
        idepg = 0
        sig(:,:) = stry(:,:) + p_try*DELTA(:,:)
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
c
c   === Von Mises stress, energies ====================================
      smean = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
      sd(:,:) = sig(:,:) - smean*DELTA(:,:)
      vons = 0.d0
      do jj=1,3
        do ii=1,3
          vons = vons + sd(ii,jj)*sd(ii,jj)
        enddo
      enddo
      vons = dsqrt(1.5d0*vons)
c
      eel(:,:) = str(:,:) - plstrg(:,:)
      e_dns = 0.d0
      do jj=1,3
        do ii=1,3
          e_dns = e_dns + eel(ii,jj)*sig(ii,jj)
        enddo
      enddo
      e_dns = 0.5d0*e_dns
c
      p_dns = 0.d0
      do jj=1,3
        do ii=1,3
          p_dns = p_dns + deltag*oun(ii,jj)*sig(ii,jj)
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
      return
      end
