      subroutine initt6(   nel,  nelx,   ndf,  node, ngaus,
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
      dimension wg(3),xg(3),yg(3)
      dimension dn(2,6),dndx(2,6)
      dimension dxdl(2,2),dldx(2,2)
      dimension ctens(3,3,3,3)
c
c **********************************************************************
      wg(1) = 1.d0/6.d0
      wg(2) = 1.d0/6.d0
      wg(3) = 1.d0/6.d0
c
      xg(1) = 1.d0/6.d0
      xg(2) = 2.d0/3.d0
      xg(3) = 1.d0/6.d0
c
      yg(1) = 1.d0/6.d0
      yg(2) = 1.d0/6.d0
      yg(3) = 2.d0/3.d0
c
c ***** Initialization *************************************************
c
c ***** Start Numerical Integration ************************************
      ig = 0
      do 100 ix =1,3
       xl = xg(ix)
       yl = yg(ix)
       zl = 1.d0 - xl -yl
       wx= wg(ix)
       ig= ig +1
c
c       === Derivative of Shape Function r.w.t. Local Coordinate ===
c           ( dn_ij = d(N^j)/d(l_i) )
          dn(1,1) =  4.d0*xl -1.d0
          dn(2,1) =  0.d0
          dn(1,2) =  0.d0
          dn(2,2) =  4.d0*yl -1.d0
          dn(1,3) =  1.d0-4.d0*zl 
          dn(2,3) =  1.d0-4.d0*zl
          dn(1,4) =  4.d0*yl
          dn(2,4) =  4.d0*xl
          dn(1,5) =  -4.d0*yl
          dn(2,5) =  4.d0*(zl-yl)
          dn(1,6) =  4.d0*(zl-xl)
          dn(2,6) =  -4.d0*xl
c
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
c         ( Determinant => Jacobian )
          CALL inv_22(  dxdl,   det,  dldx,ierror )
          if(det.le.0.d0) CALL zerodt(   nel,    ig,   det)
          det_g(ig,nel) = det
c
          if(det.le.0.d0) then
            WRITE(*,9000) nel,ig,det
            ierror = 13
            RETURN
          endif
c
c       === Derivative of Local Coordinate r.w.t. Current Position ===
c           ( dldx_ij = (dxdl_ij)^{-1} )
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
  100 continue
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
