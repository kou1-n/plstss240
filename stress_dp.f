      SUBROUTINE stress_dp(itrmax, idepg,
     &                     prope, sig, str, ehist,
     &                     ctol, vons, e_dns, p_dns,
     &                     ctens,
     &                     ierror)
c *********************************************************************
c *                                                                   *
c *   STRESS Update Routine for Drucker-Prager Model                  *
c *   with Combined Nonlinear Isotropic/Kinematic Hardening           *
c *   Using Return-Mapping Algorithm with Two Unknowns                *
c *                                                                   *
c *   Input Variables:                                                *
c *     itrmax    : maximum number of Newton iterations               *
c *     idepg     : plastic flag (returned)                           *
c *     prope(20) : material constants (unused)                       *
c *     sig(3,3)  : stress tensor                                     *
c *     str(3,3)  : strain tensor                                     *
c *     ehist(20) : history variables                                 *
c *     ctol      : convergence tolerance for Newton iteration        *
c *                                                                   *
c *   Output Variables:                                               *
c *     vons      : von Mises stress                                  *
c *     e_dns     : elastic energy density                            *
c *     p_dns     : plastic work density                              *
c *     ctens(3,3,3,3) : constitutive tensor                          *
c *     ehist(20) : updated history variables                         *
c *     ierror    : error flag                                        *
c *                                                                   *
c **********************************************************************
      implicit none

c     External hardening functions
      external H_iso, dH_iso_dk, H_kin, dH_kin_dk
      double precision H_iso, dH_iso_dk, H_kin, dH_kin_dk

c     Material parameters are taken from the MATER block
c     (Young modulus, Poisson ratio, etc.)
c
c ====================== Local Variables ===============================
      integer i,j,k,l,it
      integer itrmax, idepg, ierror

      real*8  prope(20)
      real*8  ctol
      real*8  str(3,3),sig(3,3)
      real*8  ctens(3,3,3,3)
      real*8  ehist(20)
      real*8  vons,e_dns,p_dns
c
      real*8  DELTA(3,3),FIT(3,3,3,3)
      real*8  plstrg(3,3),betaeg(3,3)
      real*8  stry(3,3),seta(3,3)
      real*8  eel(3,3)
      real*8  N(3,3),oun(3,3)
c     Material constants
      real*8  yng,poi,vmu,vlm,vkp
      real*8  yld,hk,hpa,hpb,hpd
      real*8  phi_dp,psi_dp,eta_dp,xi_dp,etabar_dp
c
      real*8  alpeg,alptmp,alpd
      real*8  deltag,stno
      real*8  hrdtmp,dhdtmp
      real*8  tmpkrd,dkdtmp
      real*8  dhard,dkard
      real*8  etrs,smean
      real*8  ptry,ftreg
      real*8  theta,thetab,A
c
c --- Variables for 2x2 Newton-Raphson system
      real*8  R(2),dR(2,2),dx(2)
      real*8  sig_eq,sig_eq_trial
      real*8  det
c     tolerance for Jacobian singularity check
      real*8, parameter :: TOL_DET = 1.0d-12
c
c ===================== Initialize Variables ==========================
      ierror = 0
      idepg = 0
c
c --- Set material parameters from the prope array --------------------
      yng    = prope(1)
      poi    = prope(2)
      vmu    = yng/(2.d0*(1.d0 + poi))
      vlm    = poi*yng/((1.d0 + poi)*(1.d0 - 2.d0*poi))
      vkp    = yng/(3.d0*(1.d0 - 2.d0*poi))
      yld    = prope(10)
      hk     = prope(11)
      hpa    = prope(12)
      hpb    = prope(13)
      hpd    = prope(15)
      phi_dp = prope(16)
      psi_dp = prope(17)

      eta_dp   =  sin(phi_dp) / ( dsqrt(3.d0) * (3.d0 - sin(phi_dp)) )
      xi_dp    = 6.d0 * cos(phi_dp) / ( dsqrt(3.d0) *
     &                                   (3.d0 - sin(phi_dp)) )
      etabar_dp = 6.d0 * sin(psi_dp) / ( dsqrt(3.d0) *
     &                                   (3.d0 - sin(psi_dp)) )
c
c --- Initialize strain/stress tensors
      deltag = 0.d0
      do j=1,3
        do i=1,3
          plstrg(i,j) = 0.d0
          betaeg(i,j) = 0.d0
          stry(i,j) = 0.d0
          seta(i,j) = 0.d0
          eel(i,j) = 0.d0
          N(i,j) = 0.d0
          oun(i,j) = 0.d0
        enddo
      enddo
c
c --- Initialize history variables
      alpeg = ehist(1)
      k = 1
      do i=1,3
        do j=1,3
          k = k+1
          plstrg(j,i) = ehist(k)
        enddo
      enddo
      do i=1,3
        do j=1,3
          k = k+1
          betaeg(j,i) = ehist(k)
        enddo
      enddo
c
c --- Initialize unit tensors
      do j=1,3
        do i=1,3
          DELTA(i,j) = 0.d0
        enddo
      enddo
      do i=1,3
        DELTA(i,i) = 1.d0
      enddo
      do l=1,3
        do k=1,3
          do j=1,3
            do i=1,3
              FIT(i,j,k,l) = 0.5d0*( DELTA(i,k)*DELTA(j,l)
     &                             +DELTA(i,l)*DELTA(j,k) )
            enddo
          enddo
        enddo
      enddo
c
c ================= Elastic Predictor Step ===========================
c --- Compute trial stress
      etrs = str(1,1) + str(2,2) + str(3,3)
      do j=1,3
        do i=1,3
          stry(i,j) = 2.d0*vmu*(str(i,j) 
     &               - (1.d0/3.d0)*etrs*DELTA(i,j))
     &               + vkp*etrs*DELTA(i,j)
        enddo
      enddo
c
c --- Compute trial deviatoric stress
      ptry = (stry(1,1) + stry(2,2) + stry(3,3))/3.d0
      do j=1,3
        do i=1,3
          seta(i,j) = stry(i,j) - ptry*DELTA(i,j)
        enddo
      enddo
c
c --- Compute trial stress norm
      stno = 0.d0
      do j=1,3
        do i=1,3
          stno = stno + seta(i,j)*seta(i,j)
        enddo
      enddo
      stno = dsqrt(stno)
      sig_eq_trial = dsqrt(1.5d0)*stno
c
c --- Check yield condition
      ftreg = sig_eq_trial
     &      - dsqrt(2.d0/3.d0)
     &        *( H_iso(alpeg, yld, hk, hpa, hpb)
     &          + H_kin(alpeg, hk, hpd) )
     &      - eta_dp*ptry
c
c ================= Plastic Corrector Step ==========================
      if(ftreg.gt.0.d0) then
        idepg = 1
c     --- Plastic state detected (output managed at step level)
c
c --- Initialize unknowns
        deltag = 0.d0
        sig_eq = sig_eq_trial
        alptmp = alpeg
c
c --- Newton-Raphson iteration for two unknowns
        do it=1,itrmax
          alpd = alptmp - alpeg
c
c --- Compute hardening terms
          hrdtmp = H_iso(alptmp, yld, hk, hpa, hpb) - yld
          dhdtmp = dH_iso_dk(alptmp, yld, hk, hpa, hpb)
c
          tmpkrd = H_kin(alptmp, hk, hpd) - H_kin(alpeg, hk, hpd)
          dkdtmp = dH_kin_dk(hk, hpd)
c
c --- Compute residuals
          R(1) = sig_eq 
     &         - sig_eq_trial 
     &         + dsqrt(2.d0)*vmu*deltag
          R(2) = sig_eq
     &         - dsqrt(2.d0/3.d0)
     &           *( H_iso(alptmp, yld, hk, hpa, hpb)
     &              + H_kin(alptmp, hk, hpd) )
     &         - eta_dp*ptry
     &         + eta_dp*vkp*etabar_dp*deltag
c
c --- Compute Jacobian
          dR(1,1) = dsqrt(2.d0)*vmu
          dR(1,2) = 1.d0
          dR(2,1) = eta_dp*vkp*etabar_dp
     &               - dsqrt(2.d0/3.d0)*(dhdtmp + dkdtmp)*xi_dp
          dR(2,2) = 1.d0
c
c --- Solve 2x2 system
          det = dR(1,1)*dR(2,2) - dR(1,2)*dR(2,1)
          if(dabs(det).lt.TOL_DET) then
            ierror = 21
            return
          endif
          dx(1) = (-R(1)*dR(2,2) + R(2)*dR(1,2))/det
          dx(2) = (-R(2)*dR(1,1) + R(1)*dR(2,1))/det
c
c --- Update unknowns
          deltag = deltag + dx(1)
          sig_eq = sig_eq + dx(2)
          alptmp = alpeg + xi_dp*deltag
c
c --- Check convergence
          if(dabs(dx(1)).lt.ctol .and. dabs(dx(2)).lt.ctol) then
            goto 210
          endif
        enddo
c
c --- Error section: Failed to converge
        ierror = 17
        return
c
  210   continue
c
c --- Update equivalent plastic strain
        alptmp = alpeg
        alpeg = alpeg + xi_dp*deltag
c
c --- Compute outward normal
        do j=1,3
          do i=1,3
            oun(i,j) = seta(i,j)/stno
            N(i,j) = seta(i,j)/(dsqrt(2.d0)*stno)
     &             + etabar_dp/3.d0*DELTA(i,j)
          enddo
        enddo
c
c --- Update plastic strain and backstress
        alpd = alpeg - alptmp
        tmpkrd = (1.d0-hpd)*hk*alpd
        etrs = str(1,1) + str(2,2) + str(3,3)
c
        do j=1,3
          do i=1,3
            betaeg(i,j) = betaeg(i,j) 
     &                  + deltag*dsqrt(2.d0)*tmpkrd*N(i,j)
            plstrg(i,j) = plstrg(i,j) + deltag*N(i,j)
            sig(i,j) = stry(i,j) 
     &               - dsqrt(2.d0)*vmu*deltag*oun(i,j)
     &               + (ptry-vkp*etabar_dp*deltag)*DELTA(i,j)
          enddo
        enddo
c
c --- Update constitutive tensor
        dhard = dH_iso_dk(alpeg, yld, hk, hpa, hpb)
        dkard = dH_kin_dk(hk, hpd)
c
        A = 1.d0/(vmu + vkp*etabar_dp*eta_dp 
     &          + xi_dp*xi_dp*(dhdtmp + dkdtmp))
c
        theta = 1.d0 - (dsqrt(2.d0)*vmu*deltag)/stno
        thetab = (dsqrt(2.d0)*vmu*deltag)/stno - vmu*A
c
        do l=1,3
          do k=1,3
            do j=1,3
              do i=1,3
                ctens(i,j,k,l) = vkp*(1.d0-vkp*eta_dp*etabar_dp*A)
     &                         *DELTA(i,j)*DELTA(k,l)
     &                         + 2.d0*vmu*theta*(FIT(i,j,k,l)
     &                         - (1.d0/3.d0)*DELTA(i,j)*DELTA(k,l))
     &                         + 2.d0*vmu*thetab*oun(i,j)*oun(k,l)
     &                         - dsqrt(2.d0)*vmu*A*vkp
     &                         *(eta_dp*oun(i,j)*DELTA(k,l)
     &                         + etabar_dp*DELTA(i,j)*oun(k,l))
              enddo
            enddo
          enddo
        enddo
c
c ================= Elastic Case ====================================
      else
        idepg = 0
c
        etrs = str(1,1) + str(2,2) + str(3,3)
        do j=1,3
          do i=1,3
            sig(i,j) = stry(i,j) + vkp*etrs*DELTA(i,j)
          enddo
        enddo
c
c --- Update constitutive tensor
        do l=1,3
          do k=1,3
            do j=1,3
              do i=1,3
                ctens(i,j,k,l) = vkp*DELTA(i,j)*DELTA(k,l)
     &                         + 2.d0*vmu*(FIT(i,j,k,l)
     &                         - (1.d0/3.d0)*DELTA(i,j)*DELTA(k,l))
              enddo
            enddo
          enddo
        enddo
      endif
c
c ================= Post-Processing =================================
c --- Compute von Mises stress
      smean = (sig(1,1) + sig(2,2) + sig(3,3))/3.d0
      do j=1,3
        do i=1,3
          seta(i,j) = sig(i,j) - smean*DELTA(i,j)
        enddo
      enddo
c
      vons = 0.d0
      do j=1,3
        do i=1,3
          vons = vons + seta(i,j)*seta(i,j)
        enddo
      enddo
      vons = dsqrt(1.5d0*vons)
c
c --- Compute energy density
      do j=1,3
        do i=1,3
          eel(i,j) = str(i,j) - plstrg(i,j)
        enddo
      enddo
c
      e_dns = 0.d0
      do j=1,3
        do i=1,3
          e_dns = e_dns + eel(i,j)*sig(i,j)
        enddo
      enddo
      e_dns = 0.5d0*e_dns
c
      p_dns = 0.d0
      do j=1,3
        do i=1,3
          p_dns = p_dns + deltag*oun(i,j)*sig(i,j)
        enddo
      enddo
c
c --- Store history variables
      ehist(1) = alpeg
      k = 1
      do i=1,3
        do j=1,3
          k = k+1
          ehist(k) = plstrg(j,i)
        enddo
      enddo
      do i=1,3
        do j=1,3
          k = k+1
          ehist(k) = betaeg(j,i)
        enddo
      enddo
c
      return
      end 
