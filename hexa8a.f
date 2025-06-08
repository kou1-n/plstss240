      subroutine hexa8a(   nel,  nelx,   ndf,  node, ngaus,
     &                   prope,   ske,  ijke,    xe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
c
      dimension ijke(8),idep(8)
c
      dimension prope(20)
      dimension xe(3,8),dndx(3,8)
      dimension wg(2),xg(2)
      dimension ctens(3,3,3,3)
      dimension ske(32,32)
c
c **********************************************************************
      wg(1) = 1.d0
      wg(2) = 1.d0
      xg(1) = -1.d0/dsqrt(3.d0)
      xg(2) =  1.d0/dsqrt(3.d0)
c
c ***** Initialization *************************************************
      pi = 4.d0*datan(1.d0)
c
c ***** Start Numerical Integration ************************************
      ig = 0
      do 100 ix=1,2
        xl = xg(ix)
        wx = wg(ix)
        do 110 iy=1,2
          yl = xg(iy)
          wy = wg(iy)
          do 120 iz=1,2
            zl = xg(iz)
            wz = wg(iz)
            ig = ig +1

c
c         === Derivative of Shape Function r.w.t. Local Coordinate ===
c             ( dn_ij = d(N^j)/d(l_i) )
c
c         === Derivative of Current Position r.w.t. Loacal Coordinate ===
c             ( dxdl_ij = d(x_i)/d(l_j) )
c
c           ( Determinant => Jacobian )
            det = det_g(ig,nel)
            detwxyz = det*wx*wy*wz
c
c         === Derivative of Local Coordinate r.w.t. Current Position ===
c           ( dldx_ij = (dxdl_ij)^{-1} )
c
c         === Derivative of Shape Function r.w.t. Current Position ===
c             ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
            dndx(:,:) = dndx_g(:,:,ig,nel)
c
c
c         === Set Elastic-Plastic Consistent Tangent ===
            ctens(:,:,:,:) = ctensg(:,:,:,:,ig,nel)
c
c         === Constract Element Stiffness Matrix: SKE_uu ===
            do ia=1,node
              do ii=1,ndf
                iai = ndf*(ia-1) +ii
                do jb=1,node
                  do jj=1,ndf
                    jbj = ndf*(jb-1) +jj
c
                    sum = 0.d0
                    do kk=1,ndf
                      do ll=1,ndf
                        sum = sum
     &                      +dndx(kk,ia)*ctens(ii,kk,jj,ll)*dndx(ll,jb)
                      enddo
                    enddo
                    ske(iai,jbj) = ske(iai,jbj) +sum*detwxyz
c
                  enddo
                enddo
              enddo
            enddo
c
  120     continue
  110   continue
  100 continue
c
c **********************************************************************
 9000   FORMAT('************************************************',/,
     &         '************************************************',/,
     &         '      Determinant becomes ZERO!! in :',i8,i5,e13.5)
c **********************************************************************
c **********************************************************************
      RETURN
      END
