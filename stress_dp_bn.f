      subroutine stress_dp_bn(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens, g_val,
     &                  ierror,  itr , histi )
c
c **********************************************************************
c     Block Newton method implementation for Drucker-Prager plasticity
c     Based on: Yamamoto et al. (2021) Box 2 formulation
c     
c     Features:
c     - 1-variable system (deltag only) since alpha = alpha0 + xi*deltag
c     - No local iterations - uses global Block Newton approach
c     - Pseudo stress correction for yield condition residual
c **********************************************************************
c
      implicit double precision(a-h,o-z)
c
      dimension prope(20)
c
c     --- common for tracking loading step and element info ---
      common /debug_info/ nel_current, ig_current, lstep_current
      
      dimension str(3,3),stry(3,3),
     &          sig(3,3),oun(3,3),sd(3,3),eel(3,3),
     &          plstrg(3,3), psig(3,3), stri(3,3)
      dimension ctens(3,3,3,3)
      dimension ehist(20), histi(50)
c
c     --- Block Newton variables (1-variable system for D-P) ---
      REAL*8 deltag, alpeg_new, Dg
      REAL*8 deltagi, alpegi, ftri
      REAL*8 gg, xa(3,3)
c     --- L, M, N matrices for Box 2 formulation ---
      dimension L_mat(3,3), M_mat(3,3)
      REAL*8 N_scalar
c     --- Drucker-Prager parameters ---
      REAL*8 eta_dp, xi_dp, etabar_dp, phi_dp, psi_dp
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
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
c
c     --- Load iteration history from previous global iteration ---
      if(itr.gt.1) then
        deltagi = histi(1)
        alpegi  = histi(2)
        ftri    = histi(3)  ! Previous yield function value
        kk = 3
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
        deltagi = 0.d0
        alpegi  = alpeg
        ftri    = 0.d0
        stri    = str
        xa      = 0.d0
      endif
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
c     DEBUG: Check input angles before calculation
      if(nel_current.eq.1 .and. ig_current.eq.1
     &   .and. lstep_current.eq.1) then
        write(*,'(A,2E12.4)') ' INPUT phi_dp,psi_dp:', phi_dp, psi_dp
      endif
c     === Drucker-Prager parameters ===
      eta_dp   = 6.d0 * sin(phi_dp) / ( dsqrt(3.d0) *
     &                                 (3.d0 - sin(phi_dp)) )
      xi_dp    = 6.d0 * cos(phi_dp) / ( dsqrt(3.d0) *
     &                                 (3.d0 - sin(phi_dp)) )
      etabar_dp = 6.d0 * sin(psi_dp) / ( dsqrt(3.d0) *
     &                                  (3.d0 - sin(psi_dp)) )
c     DEBUG: Print calculated parameters (temporarily disabled)
c      if(nel_current.eq.1 .and. ig_current.eq.1
c     &   .and. lstep_current.eq.1) then
c        write(*,'(A)') ' === CALCULATED D-P PARAMETERS ==='
c        write(*,'(A,3E12.4)') ' eta,xi,etabar:', eta_dp, xi_dp, etabar_dp
c        write(*,'(A,2E12.4)') ' vmu,vkp:', vmu, vkp
c      endif
c
c ***** Initialization *************************************************
      deltag = 0.d0
      psig   = 0.d0
      ctens  = 0.d0
      ierror = 0
c
c   === BOX 1 STEP 1: Compute Trial Elastic Stress ===
      etrs = str(1,1) +str(2,2) +str(3,3)
      emean = etrs/3.d0
c
c   Trial deviatoric stress: s^tr = 2μ{dev[ε] - ε^p}
      stry = 2.d0*vmu*(str -emean*DELTA -plstrg)
c
c   Trial relative stress: ξ^tr = s^tr - β (back stress from history)
c   For now, assuming β = 0 (can be extended for kinematic hardening)
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +stry(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)
c
c   Unit normal
      if(stno.gt.1.0d-16) then
        oun(:,:) = stry(:,:)/stno
      else
        oun(:,:) = 0.d0
      endif
c
c  === Compute Hardening Function & Yield Function ===
c     === Complete yield function evaluation for plastic check ===
      hard = hpd*(hk*alpeg
     &     +(hpa -yld) *(1.d0 -dexp(-hpb*alpeg)))
      dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg))
c     DEBUG: Check hardening derivative
      if(nel_current.eq.1 .and. ig_current.eq.1) then
        write(*,'(A,I3,A,5E12.4)') ' Step',lstep_current,
     &    ' hpd,hk,hpb,hpa,yld:', hpd, hk, hpb, hpa, yld
        write(*,'(A,I3,A,2E12.4)') ' Step',lstep_current,
     &    ' alpeg,dhard:', alpeg, dhard
      endif
      ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*vkp*etrs
     &      - xi_dp*(yld +hard)
c
c  ===== PLASTIC CASE: Block Newton Method =====
      if(ftreg.gt.-ctol) then
        idepg = 1
c
        if(itr.eq.1) then
c         === BOX 2 STEP 1: INITIALIZE ===
c         Following Yamamoto et al. (2021) Box 2 initialization exactly
c
c         1. Δγ⁽⁰⁾ = 0
          deltag = 0.d0
c
c         2. For Block Newton: evaluate initial yield function
c            Initially deltag=0, so similar to trial state evaluation
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
c         DEBUG: Analyze N_scalar components in detail
          term1 = 2.d0*vmu
          term2 = vkp*etabar_dp*eta_dp
          term3 = xi_dp*xi_dp*dhard
          if(nel_current.eq.1 .and. ig_current.eq.1) then
            write(*,'(A,I3)') ' === N_scalar DEBUG Step:', lstep_current
            write(*,'(A,E12.4)') '   2*mu        =', term1
            write(*,'(A,E12.4)') '   kappa*eta^2 =', term2
            write(*,'(A,E12.4)') '   xi^2*dhard  =', term3
            write(*,'(A,E12.4)') '   vmu         =', vmu
            write(*,'(A,E12.4)') '   vkp         =', vkp
            write(*,'(A,E12.4)') '   etabar_dp   =', etabar_dp
            write(*,'(A,E12.4)') '   eta_dp      =', eta_dp
            write(*,'(A,E12.4)') '   xi_dp       =', xi_dp
            write(*,'(A,E12.4)') '   dhard       =', dhard
          endif
c         N = -{2μ + κβα + 硬化項} (Drucker-Prager specific)
          N_scalar = -(term1 + term2 + term3)
          if(nel_current.eq.1 .and. ig_current.eq.1) then
            write(*,'(A,E12.4)') '   N_total     =', N_scalar
            write(*,'(A,E12.4)') '   |N_total|   =', dabs(N_scalar)
            if(dabs(N_scalar).lt.1.d-6) then
              write(*,'(A)') ' *** CRITICAL: N_scalar very small! ***'
            endif
          endif
c         Prevent singular matrix (ensure N_scalar is not too small)
          if(dabs(N_scalar).lt.1.d-12) then
            write(*,'(A,E12.4)') ' WARNING: N_scalar too small:', N_scalar
            N_scalar = -2.d0*vmu
            write(*,'(A,E12.4)') ' Reset to -2*mu:', N_scalar
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
c         DEBUG: Print deltag value
          if(nel_current.eq.1 .and. ig_current.eq.1) then
            write(*,'(A,I3,A,2E12.4)') ' Step',lstep_current,
     &        ' deltagi,deltag:', deltagi, deltag
          endif

c         === BOX 1: Update state variables ===
c         Update plastic strain: εᵖ = εᵖₙ + Δγ·n
          plstrg(:,:) = plstrg(:,:) + deltag*oun(:,:)*dsqrt(1.d0/2.d0)
     &                + deltag*etabar_dp*DELTA(:,:)/3.d0
c         Update equivalent plastic strain: α = αₙ + √(2/3)·Δγ
          alpeg_new = alpeg + xi_dp*deltag
c         Update hardening derivative for updated state
          dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg_new))
c         Update N_scalar for current state
          N_scalar = -(2.d0*vmu + vkp*etabar_dp*eta_dp
     &              + xi_dp*xi_dp*dhard)
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
c       DEBUG: Print yield function components
        if(nel_current.eq.1 .and. ig_current.eq.1) then
          write(*,'(A,I3,A,3E12.4)') ' Step',lstep_current,
     &      ' q,p,hard:', dsqrt(1.d0/2.d0)*stno, eta_dp*smean, hard
        endif
        g_val = dsqrt(1.d0/2.d0)*stno + eta_dp*smean
     &        - xi_dp*(yld + hard)
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
        if(dabs(N_scalar).gt.1.d-16) then
          N_inv = 1.d0 / N_scalar
        else
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
        histi(1) = deltag
        histi(2) = alpeg_new
        histi(3) = ftreg
        kk = 3
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
        histi(1) = 0.d0
        histi(2) = alpeg
        histi(3) = 0.d0
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