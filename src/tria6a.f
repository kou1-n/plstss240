      subroutine tria6a(   nel,  nelx,   ndf,  node, ngaus,
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
c
      dimension prope(20)
      dimension xe(3,8)
      dimension wg(3),xg(3),yg(3)
      dimension dn(2,6),dndx(2,6)
      dimension dxdl(2,2),dldx(2,2)
      dimension ctens(3,3,3,3)
      dimension ske(32,32)
c
c **********************************************************************
c

      ske=0.d0
      wg(1) = 1.d0/6.d0
      wg(2) = 1.d0/6.d0
      wg(3) = 1.d0/6.d0
c       
      xg(1) = 1.d0/6.d0
      xg(2) = 2.d0/3.d0
      xg(3) = 1.d0/6.d0
c       
      yg(1) =  1.d0/6.d0
      yg(2) =  1.d0/6.d0
      yg(3) =  2.d0/3.d0
c ***** Initialization *************************************************
      ig = 0
c ***** Set Global Coordinate ******************************************
      do 100 ix = 1,3
       xl = xg(ix)
       yl = xg(ix)
       zl = 1-xl-yl
         wx = wg(ix)
         ig = ig +1
c   === Area of this element ===
      det = det_g(ig,nel)
      detwxyz = det*wx
c
c ***** Start Numerical Integration ************************************
c
c   === Derivative of Shape Function r.w.t. Current Position ===
c       ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
      dndx(:,:) = dndx_g(:,:,ig,nel)
c
c   === Set Elastic-Plastic Consistent Tangent ===
      ctens(:,:,:,:) = ctensg(:,:,:,:,ig,nel)
c
c   === Constract Element Stiffness Matrix: SKE_uu ===
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
     &                +dndx(kk,ia)*dndx(ll,jb)*ctens(ii,kk,jj,ll)
                enddo
              enddo
              ske(iai,jbj) = ske(iai,jbj) +sum*detwxyz
c
            enddo
          enddo
        enddo
      enddo
c
  100 continue
c     write(*,*) 'ske'
c     do ii=1,12
c     write(*,'(12f10.5)') (ske(ii,jj),jj=1,12)
c     enddo
c **********************************************************************
 9000 FORMAT('************************************************',/,
     &       '************************************************',/,
     &       '      Determinant becomes ZERO!! in :',i8,i5,e13.5)
c
c **********************************************************************
c     write(*,*) ske(1,1)
c     write(*,'(/,i5)') nel
c     do no=1,12
c       do mo=1,12
c         write(*,'(2i3,2e25.17,e15.5)') no,mo,ske(no,mo),ske(mo,no),
c    &                            ske(no,mo)-ske(mo,no)
c       enddo
c     enddo
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
