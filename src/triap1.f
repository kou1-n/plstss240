      subroutine triap1(   nel,  nelx,   ndf,  node, ngaus,
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
      dimension ijke(8)
      dimension ijk_t(3,4)
c
      dimension prope(20)
      dimension xe(3,8)
c     dimension wg(2),xg(2)
      dimension dndx(2,3)
      dimension dndx_m(2,6)
c     dimension dn(2,4),dndx(2,4)
c     dimension dxdl(2,2),dldx(2,2)
      dimension ctens(3,3,3,3),c_dev(3,3,3,3)
      dimension ske(32,32)
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c
c **********************************************************************
      ske(:,:) = 0.d0
c
c ***** Set Local Conectivities ****************************************
      ijk_t(1,1) = 1
      ijk_t(2,1) = 4
      ijk_t(3,1) = 6
      ijk_t(1,2) = 4
      ijk_t(2,2) = 2
      ijk_t(3,2) = 5
      ijk_t(1,3) = 6
      ijk_t(2,3) = 4
      ijk_t(3,3) = 5
      ijk_t(1,4) = 6
      ijk_t(2,4) = 5
      ijk_t(3,4) = 3
c
c ***** Initialization *************************************************
      area = 0.d0
      c_volm = 0.d0
      dndx_m(:,:) = 0.d0
c
c **********************************************************************
c ***** Start Numerical Integration (Volumetric part) ******************
c **********************************************************************
      do 100 ig=1,4
c
        det = det_g(ig,nel)
        detwxyz = 0.5d0*det
        area = area +detwxyz
c
c   === Derivative of Shape Function r.w.t. Current Position ===
c       ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
c       dndx(:,:) = dndx_g(:,:,ig,nel)
        do no=1,3
          nom = ijk_t(no,ig)
          do nd=1,2
            dndx_m(nd,nom) = dndx_m(nd,nom) 
     &                      +dndx_g(nd,no,ig,nel)*detwxyz
          enddo
        enddo
c
c   === Set Elastic-Plastic Consistent Tangent ===
c       ctens(:,:,:,:) = ctensg(:,:,:,:,ig,nel)
c     --- for volumetric part ---
        c_vol = 0.d0
        do ii=1,3
          do jj=1,3
            do kk=1,3
              do ll=1,3
                c_vol = c_vol +DELTA(ii,jj)*ctensg(ii,jj,kk,ll,ig,nel)
     &                                     *DELTA(kk,ll)
              enddo
            enddo
          enddo
        enddo
        c_vol = c_vol/9.d0
c
        c_volm = c_volm +c_vol*detwxyz
c
c **********************************************************************
  100 CONTINUE
c
c   === Volume average ===
      c_volm = c_volm/area
      dndx_m(:,:) = dndx_m(:,:)/area
c
c   === Constract Element Stiffness Matrix: SKE_vol ===
      do ia=1,node
        do ii=1,ndf
          iai = ndf*(ia-1) +ii
          do jb=1,node
            do jj=1,ndf
              jbj = ndf*(jb-1) +jj
              ske(iai,jbj) = ske(iai,jbj) +dndx_m(ii,ia)*dndx_m(jj,jb)
c    &                                      *area*c_volm
            enddo
          enddo
        enddo
      enddo
c
      ske(:,:) = ske(:,:)*area*c_volm
c
c **********************************************************************
c ***** Start Numerical Integration (Deviatoric part) ******************
c **********************************************************************
      do 200 ig=1,4
c
        det = det_g(ig,nel)
        detwxyz = 0.5d0*det
c
c   === Derivative of Shape Function r.w.t. Current Position ===
c       ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
        dndx(:,:) = dndx_g(:,:,ig,nel)
c
c   === Set Elastic-Plastic Consistent Tangent ===
c     --- for deviatoric part ---
        ctens(:,:,:,:) = ctensg(:,:,:,:,ig,nel)
        do ii=1,3
          do jj=1,3
            do kk=1,3
              do ll=1,3
c
                sum = 0.d0
                do mm=1,3
                  do nn=1,3
                    do mo=1,3
                      do no=1,3
                    sum = sum +DTENS(mm,nn,ii,jj)*ctens(mm,nn,mo,no)
     &                                            *DTENS(mo,no,kk,ll)
                      enddo
                    enddo
                  enddo
                enddo
                c_dev(ii,jj,kk,ll) = sum
c
              enddo
            enddo
          enddo
        enddo
c
c     === Constract Element Stiffness Matrix: SKE_dev ===
        do ia=1,3
          do ii=1,ndf
c           iai = ndf*(ia-1) +ii
            iai  =ndf*(ijk_t(ia,ig)-1)+ii
            do jb=1,3
              do jj=1,ndf
c               jbj = ndf*(jb-1) +jj
                jbj  =ndf*(ijk_t(jb,ig)-1)+jj
c
                sum = 0.d0
                do kk=1,ndf
                  do ll=1,ndf
                    sum = sum
     &                  +dndx(kk,ia)*dndx(ll,jb)*c_dev(ii,kk,jj,ll)
                  enddo
                enddo
                ske(iai,jbj) = ske(iai,jbj) +sum*detwxyz
c
              enddo
            enddo
          enddo
        enddo
c
  200 CONTINUE
c
c **********************************************************************
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
