      double precision function H_iso(kappa, yld, hk, hpa, hpb)
      implicit none
      double precision kappa, yld, hk, hpa, hpb
c     Isotropic hardening law K(\kappa) as in CMLformat07.md
      H_iso = yld + hk*kappa + (hpa - yld)*(1.d0 - dexp(-hpb*kappa))
      return
      end

c      Derivative of isotropic hardening with respect to kappa
      double precision function dH_iso_dk(kappa, yld, hk, hpa, hpb)
      implicit none
      double precision kappa, yld, hk, hpa, hpb
c     Derivative dK/d\kappa for isotropic hardening
      dH_iso_dk = hk + (hpa - yld)*hpb*dexp(-hpb*kappa)
      return
      end

      double precision function H_kin(kappa, hk, theta)
      implicit none
      double precision kappa, hk, theta
c     Linear kinematic hardening component
      H_kin = (1.d0 - theta)*hk*kappa
      return
      end

      double precision function dH_kin_dk(hk, theta)
      implicit none
      double precision hk, theta
c     Derivative dH_kin/d\kappa is constant
      dH_kin_dk = (1.d0 - theta)*hk
      return
      end
