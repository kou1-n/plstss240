      subroutine plastc( prope, ctens,
     &                     vmu,   vlm,   vkp,
     &                    ftry, alpha,deltag,   sig,    bb,
     &                  ierror )
c
c ----------------------------------------------------------------------
c Purpose: Make elastic-plastic consistent tangent tensors
c ----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension prope(20)
      dimension ctens(3,3,3,3)
      dimension sig(3,3),sd(3,3),oun(3,3),bb(3,3)
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
c
c ***** Initialization *************************************************
      ctens = 0.d0 ! ctens(:,:,:,:) = 0.d0
c
c ***** Set Material Properties ****************************************
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
c     --- Plastic Parameters
      yld = prope(10)
      hk  = prope(11)
      hpa = prope(12)
      hpb = prope(13)
      hpc = prope(14)
      hpd = prope(15)
c
c ***** Set Deformation Histories **************************************
c   === Deviatoric Stress ===
      smean = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
c
      do jj=1,3
        do ii=1,3
          sd(ii,jj) = sig(ii,jj) -smean*DELTA(ii,jj) -bb(ii,jj)
        enddo
      enddo
c
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +sig(ii,jj)*sig(ii,jj)
        enddo
      enddo
      stno = dsqrt(stno)
c
c   === Outward Unit Normal Vector in Stress Space ===
      do jj=1,3
        do ii=1,3
          oun(ii,jj) = sig(ii,jj)/stno
        enddo
      enddo
c
c   === Set Plastic Parameters ===
      dhard = hpd* hk +hpb*(hpa -yld)*dexp(-hpb*alpha)
      dkard = (1.d0 -hpd)* hk

      theta = 1.d0 -(2.d0*vmu*deltag)/stno
      thetab= 1.d0/(1.d0 +(dhard+dkard)/(3.d0*vmu)) -(1.d0 -theta)
c
c ***** Set Elastic Tensor *********************************************
      do ll=1,3
        do kk=1,3
          do jj=1,3
            do ii=1,3
              ctens(ii,jj,kk,ll) 
     &          = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &            +2.d0*vmu*theta*( FIT(ii,jj,kk,ll) 
     &                          -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
     &            -2.d0*vmu*thetab*oun(ii,jj)*oun(kk,ll)
            enddo
          enddo
        enddo
      enddo
c
c **********************************************************************
c     do ii=1,3
c       do jj=1,3
c         do kk=1,3
c           do ll=1,3
c             write(*,'(4i2,e25.17)') ii,jj,kk,ll,ctens(ii,jj,kk,ll)
c           enddo
c         enddo
c       enddo
c     enddo
c **********************************************************************
c     write(*,*) ' '
c     write(*,*) ctens
c **********************************************************************
      RETURN
      END
