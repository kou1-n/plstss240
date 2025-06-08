      subroutine inv_33(    aa,   det,    bb,
     &                  ierror )
c
      implicit double precision(a-h,o-z)
c
      dimension aa(3,3),bb(3,3)
c
c **********************************************************************
c   ----- Determinant -----
      det = aa(1,1)*aa(2,2)*aa(3,3)
     &     +aa(2,1)*aa(3,2)*aa(1,3)
     &     +aa(3,1)*aa(1,2)*aa(2,3)
     &     -aa(1,3)*aa(2,2)*aa(3,1)
     &     -aa(1,1)*aa(3,2)*aa(2,3)
     &     -aa(1,2)*aa(2,1)*aa(3,3)
c
      if(det.le.0.d0) then
        ierror = 99
        RETURN
      endif
c
c   ----- Inverse -----
      bb(1,1) = aa(2,2)*aa(3,3) -aa(2,3)*aa(3,2)
      bb(2,1) = aa(2,3)*aa(3,1) -aa(2,1)*aa(3,3)
      bb(3,1) = aa(2,1)*aa(3,2) -aa(2,2)*aa(3,1)
      bb(1,2) = aa(1,3)*aa(3,2) -aa(1,2)*aa(3,3)
      bb(2,2) = aa(1,1)*aa(3,3) -aa(1,3)*aa(3,1)
      bb(3,2) = aa(1,2)*aa(3,1) -aa(1,1)*aa(3,2)
      bb(1,3) = aa(1,2)*aa(2,3) -aa(1,3)*aa(2,2)
      bb(2,3) = aa(1,3)*aa(2,1) -aa(1,1)*aa(2,3)
      bb(3,3) = aa(1,1)*aa(2,2) -aa(1,2)*aa(2,1)
c
      deti = 1.d0/det
      bb = bb*deti ! bb(:,:) = bb(:,:)*deti
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
