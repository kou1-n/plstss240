      subroutine initt4(   nel,  nelx,   ndf,  node, ngaus,
     &                   prope,  ijke,    xe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
c
      dimension ijke(8)
c
      dimension prope(20)
      dimension xe(3,8),dn(3,8),dndx(3,8)
      dimension dxdl(3,3),dldx(3,3)
      dimension wg(2), xg(2)
      dimension ctens(3,3,3,3)
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
      wg(1) = 1.d0
      wg(2) = 1.d0
      xg(1) = -1.d0/dsqrt(3.d0)
      xg(2) =  1.d0/dsqrt(3.d0)
c
c ***** Initialization *************************************************
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
            dn(1,1) = -0.125d0*(1.d0 -yl)*(1.d0 -zl)
            dn(2,1) = -0.125d0*(1.d0 -xl)*(1.d0 -zl)
            dn(3,1) = -0.125d0*(1.d0 -xl)*(1.d0 -yl)
c
            dn(1,2) =  0.125d0*(1.d0 -yl)*(1.d0 -zl)
            dn(2,2) = -0.125d0*(1.d0 +xl)*(1.d0 -zl)
            dn(3,2) = -0.125d0*(1.d0 +xl)*(1.d0 -yl)
c
            dn(1,3) =  0.125d0*(1.d0 +yl)*(1.d0 -zl)
            dn(2,3) =  0.125d0*(1.d0 +xl)*(1.d0 -zl)
            dn(3,3) = -0.125d0*(1.d0 +xl)*(1.d0 +yl)
c
            dn(1,4) = -0.125d0*(1.d0 +yl)*(1.d0 -zl)
            dn(2,4) =  0.125d0*(1.d0 -xl)*(1.d0 -zl)
            dn(3,4) = -0.125d0*(1.d0 -xl)*(1.d0 +yl)
c
            dn(1,5) = -0.125d0*(1.d0 -yl)*(1.d0 +zl)
            dn(2,5) = -0.125d0*(1.d0 -xl)*(1.d0 +zl)
            dn(3,5) =  0.125d0*(1.d0 -xl)*(1.d0 -yl)
c
            dn(1,6) =  0.125d0*(1.d0 -yl)*(1.d0 +zl)
            dn(2,6) = -0.125d0*(1.d0 +xl)*(1.d0 +zl)
            dn(3,6) =  0.125d0*(1.d0 +xl)*(1.d0 -yl)
c
            dn(1,7) =  0.125d0*(1.d0 +yl)*(1.d0 +zl)
            dn(2,7) =  0.125d0*(1.d0 +xl)*(1.d0 +zl)
            dn(3,7) =  0.125d0*(1.d0 +xl)*(1.d0 +yl)
c
            dn(1,8) = -0.125d0*(1.d0 +yl)*(1.d0 +zl)
            dn(2,8) =  0.125d0*(1.d0 -xl)*(1.d0 +zl)
            dn(3,8) =  0.125d0*(1.d0 -xl)*(1.d0 +yl)
c
c         === Derivative of Current Position r.w.t. Loacal Coordinate ===
c             ( dxdl_ij = d(x_i)/d(l_j) )
            do ndia=1,ndf
              do ndjb=1,ndf
                dxdlj = 0.d0
                do ia=1,node
                  dxdlj = dxdlj +xe(ndia,ia)*dn(ndjb,ia)
                enddo
                dxdl(ndia,ndjb) = dxdlj
              enddo
            enddo
c
c           ( Determinant => Jacobian )
            det = 0.d0
            do 200 ia=1,3
              do ii=1,3
                ddj = 0.d0
                do 210 jb=1,3
                  do jj=1,3
                    dxdlj = dxdl(jb,jj)
                    do kc=1,3
                      pijk = EPSLN(ia,jb,kc)
                      do kk=1,3
                        sss = pijk*EPSLN(ii,jj,kk)*dxdlj*dxdl(kc,kk)
                        det = det +sss*dxdl(ia,ii)
                        ddj = ddj +sss
                      enddo
                    enddo
                  enddo
  210           continue
                dldx(ii,ia) = ddj
              enddo
  200       continue

c
            det = det/6.d0
            det_g(ig,nel) = det
            detwxyz = det*wx*wy*wz
c
            if(det.le.0.d0) then
              WRITE(*,9000) nel,ig,det
              ierror = 15
              RETURN
            endif
c
c         === Derivative of Local Coordinate r.w.t. Current Position ===
c           ( dldx_ij = (dxdl_ij)^{-1} )
            deti = 1.d0/(det*2.d0)
            do jj=1,3
              do ii=1,3
                dldx(ii,jj) = dldx(ii,jj)*deti
              enddo
            enddo
c
c         do ii=1,3
c           write(*,'(3e15.5)') (dldx(ii,jj),jj=1,3)
c         enddo
c
c         === Derivative of Shape Function r.w.t. Current Position ===
c             ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
            do ndia=1,ndf
              do ia=1,node
                dndxj = 0.d0
                do ndjb=1,ndf
                  dndxj = dndxj +dn(ndjb,ia)*dldx(ndjb,ndia)
                enddo
                dndx(ndia,ia) = dndxj
                dndx_g(ndia,ia,ig,nel) = dndxj
              enddo
            enddo
c
c
c       === Set Elastic-Plastic Consistent Tangent ===
            CALL elastc( prope, ctens,
     &                     vmu,   vlm,   vkp,
     &                  ierror )
c
          ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
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
