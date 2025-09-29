      subroutine pres3d(   nel, ijk_e,ijkl_e,
     &                      xe,    fe, vsl_e,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk_e(8)
      dimension ijkl_e(4)
c
      dimension xe(3,8)
      dimension fe(32)
      dimension vsl_e(4)
c
      dimension sh(4)
      dimension dn(2,4)
      dimension dxdl(3,2)
      dimension gs(3)
      dimension wg(2),xg(2)
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
      wg(1) = 1.d0
      wg(2) = 1.d0
      xg(1) = -1.d0/dsqrt(3.d0)
      xg(2) =  1.d0/dsqrt(3.d0)
c
c ****** Initilization *************************************************
      dxdl = 0.d0 ! dxdl(:,:) = 0.d0
c
c ****** Numerical Integration (Surface Int.) **************************
      do 100 ix=1,2
        xl = xg(ix)
        wx = wg(ix)
        do 110 iy=1,2
          yl = xg(iy)
          wy = wg(iy)
c
c     ===== Shape Functions & Their Derivatives in Local Coord. =====
c         ( dn_ij = d(N^j)/d(l_i) )
          sh(1) = 0.25d0*(1.d0-xl)*(1.d0-yl)
          sh(2) = 0.25d0*(1.d0+xl)*(1.d0-yl)
          sh(3) = 0.25d0*(1.d0+xl)*(1.d0+yl)
          sh(4) = 0.25d0*(1.d0-xl)*(1.d0+yl)
c
          dn(1,1) = -0.25d0*(1.d0 -yl)
          dn(1,2) =  0.25d0*(1.d0 -yl)
          dn(1,3) =  0.25d0*(1.d0 +yl)
          dn(1,4) = -0.25d0*(1.d0 +yl)
c
          dn(2,1) = -0.25d0*(1.d0 -xl)
          dn(2,2) = -0.25d0*(1.d0 +xl)
          dn(2,3) =  0.25d0*(1.d0 +xl)
          dn(2,4) =  0.25d0*(1.d0 -xl)
c
c       === Derivative of Current Position r.w.t. Loacal Coordinate ===
c           ( dxdl_ij = d(x_i)/d(l_j) )
          do ndia=1,3
            do ndjb=1,2
              dxdlj = 0.d0
              do ia=1,4
                inn = ijkl_e(ia)
                do no=1,8
                  if(inn.eq.ijk_e(no)) then
c                   write(*,*) ia,no,ijkl_e(ia),ijk_e(no)
                    dxdlj = dxdlj +xe(ndia,no)*dn(ndjb,ia)
                  endif
                enddo
              enddo
              dxdl(ndia,ndjb) = dxdlj
            enddo
          enddo
c
c         do ii=1,3
c           write(*,'(3e15.5)') (dxdl(ii,jj),jj=1,2)
c         enddo
c
c       === G: Corresponde to det(J) in Volume Int. ===
          do kc=1,3
            ggss = 0.d0
            do ia=1,3
              do jb=1,3
                ggss = ggss +EPSLN(ia,jb,kc)*dxdl(ia,1)*dxdl(jb,2)
              enddo
            enddo
            gs(kc) = ggss
          enddo
c
          gsurf = dsqrt(gs(1)**2 +gs(2)**2 +gs(3)**2)*wx*wy
c
c         write(*,*) gsurf
c
          do 200 in=1,4
            inn = ijkl_e(in)
            do no=1,8
              if(inn.eq.ijk_e(no)) then
                ia = no
              endif
            enddo
c
            do nd=1,3
              iai = 3*(ia -1) +nd
              piai = vsl_e(nd+1)
              do jb=1,4
                fe(iai) = fe(iai) +sh(in)*sh(jb)*piai*gsurf
              enddo
            enddo
  200     continue
c
  110   continue
  100 continue
c
c     do ii=1,8
c       write(*,'(4e15.5)') (fe(3*(ii-1)+jj),jj=1,3)
c     enddo
c **********************************************************************
c **********************************************************************
      RETURN
      END
