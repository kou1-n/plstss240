      subroutine initq4(   nel,  nelx,   ndf,  node, ngaus,
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
      dimension xe(3,8)
      dimension wg(2),xg(2)
      dimension dn(2,4),dndx(2,4)
      dimension dxdl(2,2),dldx(2,2)
      dimension ctens(3,3,3,3)
c
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
          ig = ig +1
c
c       === Derivative of Shape Function r.w.t. Local Coordinate ===
c           ( dn_ij = d(N^j)/d(l_i) )
          dn(1,1) = -0.25d0*(1.d0 -yl)
          dn(2,1) = -0.25d0*(1.d0 -xl)
c
          dn(1,2) =  0.25d0*(1.d0 -yl)
          dn(2,2) = -0.25d0*(1.d0 +xl)
c
          dn(1,3) =  0.25d0*(1.d0 +yl)
          dn(2,3) =  0.25d0*(1.d0 +xl)
c
          dn(1,4) = -0.25d0*(1.d0 +yl)
          dn(2,4) =  0.25d0*(1.d0 -xl)
c
c       === Derivative of Current Position r.w.t. Loacal Coordinate ===
c           ( dxdl_ij = d(x_i)/d(l_j) )
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
c       === derivative of LOCAL COORDIANTE r.w.t. current position ===
c             ( dldx_ij = d(li)/d(xj) = ( d(xi)/d(lj) )^{-1} )
          CALL inv_22(  dxdl,   det,  dldx,ierror )
          if(det.le.0.d0) CALL zerodt(   nel,    ig,   det)
c
          det_g(ig,nel) = det
          detwxyz = det*wx*wy
c
c       === Derivative of Shape Function r.w.t. Current Position ===
c           ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
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
c       === Set Elastic-Plastic Consistent Tangent ===
          CALL elastc( prope, ctens,
     &                   vmu,   vlm,   vkp,
     &                ierror )
c
          ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
c
c
  110   continue
  100 continue
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
