      subroutine initt3(   nel,  nelx,   ndf,  node, ngaus,
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
c     dimension wg(2),xg(2)
      dimension dn(2,4),dndx(2,3)
c     dimension dn(2,4),dndx(2,4)
c     dimension dxdl(2,2),dldx(2,2)
      dimension ctens(3,3,3,3)
c
c **********************************************************************
c
c ***** Initialization *************************************************
      ig = 1
c
c ***** Set Global Coordinate ******************************************
      x1 = xe(1,1)
      y1 = xe(2,1)
      x2 = xe(1,2)
      y2 = xe(2,2)
      x3 = xe(1,3)
      y3 = xe(2,3)
c
c   === Area of this element ===
c     det = ( x1*(y2 -y3) +x2*(y3 -y1) +x3*(y1 -y2) )
      det  = x1*y2 +x2*y3 +x3*y1 -(x1*y3 +x2*y1 +x3*y2)
      area = 0.5d0*det
      deti = 1.d0/det
      detwxyz = area
      det_g(ig,nel) = det
c
c ***** Start Numerical Integration ************************************
c
c   === Derivative of Shape Function r.w.t. Current Position ===
c       ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
      dndx(1,1) = (y2 -y3)*deti
      dndx(2,1) = (x3 -x2)*deti
      dndx(1,2) = (y3 -y1)*deti
      dndx(2,2) = (x1 -x3)*deti
      dndx(1,3) = (y1 -y2)*deti
      dndx(2,3) = (x2 -x1)*deti
c
      dndx_g(:,:,ig,nel) = dndx(:,:)
c
c   === Set Elastic-Plastic Consistent Tangent ===
      CALL elastc( prope, ctens,
     &               vmu,   vlm,   vkp,
     &            ierror )
c
      ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
c
c **********************************************************************
 9000 FORMAT('************************************************',/,
     &       '************************************************',/,
     &       '      Determinant becomes ZERO!! in :',i8,i5,e13.5)
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
