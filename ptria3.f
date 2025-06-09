      subroutine ptria3(   nel,  nelx,   ndf,  node, ngaus,
     &                    idep, prope,    xe,    ue,
     &                     sts,   stn, finte,  vone,  epse,
     &                  dhist0,dhist1,
     &                   p_ene, e_ene,tene_e,dene_p, tempe,
     &                   dtemp,dtempe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
      dimension dhist0(20,ngaus,nelx),dhist1(20,ngaus,nelx)
c
      dimension idep(8)
c
      dimension prope(20)
      dimension xe(3,8),ue(3,8)
      dimension wg(2),xg(2)
      dimension dn(2,4),dndx(2,3)
      dimension dxdl(2,2),dldx(2,2)
      dimension sts(3,3),stn(3,3)
      dimension sig(3,3),str(3,3)
      dimension ctens(3,3,3,3)
      dimension finte(24)
      dimension ehist(20)
c
      common /tvalu/ ctol,stol
c **********************************************************************
      itrmax = 20
c
c ***** Initialization *************************************************
      finte = 0.d0 ! finte(:) = 0.d0
c
      vone = 0.d0
      epse = 0.d0
      p_ene = 0.d0
      e_ene = 0.d0
c
      sts = 0.d0 ! sts(:,:) = 0.d0
      stn = 0.d0 ! stn(:,:) = 0.d0
c
c ***** Set Global Coordinate ******************************************
c   === Area of this element ===
c     det = ( x1*(y2 -y3) +x2*(y3 -y1) +x3*(y1 -y2) )
c     det_g(1,nel) = det
      det = det_g(1,nel)
c
      area = 0.5d0*det
      detwxyz = area
      deti = 1.d0/det
c
c ***** Set Material Properties ****************************************
c    --- thermal parameters
c     row = prope(3)
c     ccc = prope(4)
c     aaa = prope(5)
c
c ***** Start Numerical Integration ************************************
      ig = 1
c
c   === Initilization of Local Tensors ===
c     sig = 0.d0 ! sig(:,:) = 0.d0
      str = 0.d0 ! str(:,:) = 0.d0
c
c   === Derivative of Shape Function r.w.t. Current Position ===
c       ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
      dndx(:,:) = dndx_g(:,:,ig,nel)
c
c   === Compute the Local Strain Tensor ===
      do jj=1,2
        do ii=1,2
          eij = 0.d0
          do no=1,node
            eij = eij +0.5d0*(
     &            ue(ii,no)*dndx(jj,no) +ue(jj,no)*dndx(ii,no) )
          enddo
          str(ii,jj) = eij
        enddo
      enddo
c
c   === Set Deformation Histtories ===
      ehist(:) = dhist0(:,ig,nel)
c
c   === Compute Local Stresses (Drucker-Prager) ===
      CALL stress_dp(itrmax, idepg,
     &             prope,   sig,   str, ehist,
     &              ctol,  vons, e_dns, p_dns,
     &             ctens,
     &            ierror )
      if(ierror.ne.0) RETURN
c
c ***** Store Deformation Histories (at current state) *****************
      dhist1(:,ig,nel) = ehist(:)
      ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
      idep(ig) = idepg
c
      alpeg = ehist(1)
c ***** Compute the Elemental Value (Sum-up for Element ) **************
      vone = vone +vons*detwxyz
      epse = epse +alpeg*detwxyz
      sts = sts +sig*detwxyz ! sts(:,:) = sts(:,:) +sig(:,:)*detwxyz
      stn = stn +str*detwxyz ! stn(:,:) = stn(:,:) +str(:,:)*detwxyz
      e_ene = e_ene +e_dns*detwxyz
      p_ene = p_ene +p_dns*detwxyz
c
c ***** Compute the Internal Force Vector ******************************
      do no=1,node
        do ii=1,ndf
          ndof = ndf*(no -1) +ii
          qol = 0.d0
          do jj=1,ndf
            qol = qol +sig(jj,ii)*dndx(jj,no)
          enddo
          finte(ndof) = finte(ndof) +qol*detwxyz
        enddo
      enddo
c
c ***** Compute the Elemental Value (Volume Average) *******************
      vone = vone/area
      epse = epse/area
c
      tene_e = tene_e +e_ene
      dene_p = dene_p +p_ene
c
      e_ene=e_ene/area
      p_ene=p_ene/area
c
      sts = sts/area ! sts(:,:) = sts(:,:)/area
      stn = stn/area ! sts(:,:) = sts(:,:)/area
c
c ***** Convert the Plastic ene. into Temperature *****
c
c     bbb = aaa*row*ccc
c     if( bbb.eq.0.d0) then
c       dtempe = 0.d0
c     else
c       dtempe = (p_ene*10.d0**9)/(aaa*row*ccc)
c     endif
c     dtemp = dtemp + dtempe*area
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
