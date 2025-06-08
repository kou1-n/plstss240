      subroutine inv_22(    aa,   det,    bb,ierror )
c
      implicit double precision(a-h,o-z)
c
      dimension aa(2,2),bb(2,2)
c
c **********************************************************************
c   ----- Determinant -----
      det = aa(1,1)*aa(2,2) -aa(1,2)*aa(2,1)
c
      if(det.le.0.d0) then
        ierror = 99
        RETURN
      endif
c
c   ----- Inverse -----
      deti = 1.d0/det
      bb(1,1) =  aa(2,2)*deti
      bb(2,2) =  aa(1,1)*deti
      bb(1,2) = -aa(1,2)*deti
      bb(2,1) = -aa(2,1)*deti
c
c
c   ----- for Checking -----
c     do ii=1,3
c       do jj=1,3
c         tmp = 0.d0
c         do kk=1,3
c           tmp = tmp +aa(ii,kk)*bb(kk,jj)
c         enddo
c         write(*,*) tmp,ii,jj
c       enddo
c     enddo
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
