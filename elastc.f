      subroutine elastc( prope, ctens,
     &                     vmu,   vlm,   vkp,
     &                  ierror )
c
c ----------------------------------------------------------------------
c Purpose: Make elastic material property tensors
c ----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension prope(20)
      dimension ctens(3,3,3,3)
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
c
c ***** Set Elastic Tensor *********************************************
      do ll=1,3
        do kk=1,3
          do jj=1,3
            do ii=1,3
c             ctens(ii,jj,kk,ll) = vlm*DELTA(ii,jj)*DELTA(kk,ll)
c    &                            +2.d0*vmu*FIT(ii,jj,kk,ll)
              ctens(ii,jj,kk,ll)
     &          = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &            +2.d0*vmu*( FIT(ii,jj,kk,ll)
     &                     -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
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
c **********************************************************************
      RETURN
      END
