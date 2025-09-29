      subroutine ptrip1(   nel,  nelx,   ndf,  node, ngaus,
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
      dimension ijk_t(3,4)
c
      dimension prope(20)
      dimension xe(3,8),ue(3,8)
c     dimension wg(2),xg(2)
      dimension dn(2,4),dndx(2,4)
      dimension dndx_m(2,6)
c     dimension dxdl(2,2),dldx(2,2)
      dimension sts(3,3),stn(3,3)
      dimension sig(3,3),str(3,3),sd(3,3)
      dimension ctens(3,3,3,3)
      dimension finte(24)
      dimension ehist(20)
c
      dimension uel(3,3)
c
      common /tvalu/ ctol,stol
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
      itrmax = 20
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
      finte = 0.d0 ! finte(:) = 0.d0
c
      area = 0.d0
      vone = 0.d0
      epse = 0.d0
      p_ene = 0.d0
      e_ene = 0.d0
c
      sts = 0.d0 ! sts(:,:) = 0.d0
      stn = 0.d0 ! stn(:,:) = 0.d0
c
      emean = 0.d0
      dndx_m = 0.d0 ! dndx_m(:,:) = 0.d0
c
c **********************************************************************
c ***** Start Numerical Integration (Volumetric part) ******************
c **********************************************************************
      do 100 ig=1,4
c
c ***** Set Global Coordinate ******************************************
c     --- Define local conectivities ---
        do ii=1,3
          jj = ijk_t(ii,ig)
          do ll=1,2
            uel(ll,ii) = ue(ll,jj)
          enddo
        enddo
c
c   === Area of this element ===
        det = det_g(ig,nel)
        detwxyz = 0.5d0* det
        area = area +detwxyz
c
c   === Initilization of Local Tensors ===
c       sig = 0.d0 ! sig(:,:) = 0.d0
        str = 0.d0 ! str(:,:) = 0.d0
c
c ***** Compute Local Strain *******************************************
c   === Volume Average of dndx ===
        dndx(:,:) = dndx_g(:,:,ig,nel)
        do no=1,3
          nom = ijk_t(no,ig)
          do nd=1,2
            dndx_m(nd,nom) = dndx_m(nd,nom)
     &                      +dndx(nd,no)*detwxyz
          enddo
        enddo
c
c   === Compute the Local Strain Tensor ===
        do jj=1,2
          do ii=1,2
            eij = 0.d0
            do no=1,3
              eij = eij +0.5d0*(
     &            uel(ii,no)*dndx(jj,no)+uel(jj,no)*dndx(ii,no))
            enddo
            str(ii,jj) = eij
          enddo
        enddo
c
c   === Volumetric part ===
        emean = emean +( str(1,1) +str(2,2) +str(3,3) )/3.d0*detwxyz
c
  100 CONTINUE
c
c ***** Compute Local Volumetirc Stress ********************************
      emean = emean/area
      dndx_m = dndx_m/area ! dndx_m(:,:) = dndx_m(:,:)/area
c
c   === Set Deformation Histories ===
      ehist = 0.d0 ! ehist(:) = 0.d0
      str = 0.d0 ! str(:,:) = 0.d0
c
c   === Compute the Averaged Volumetric Strain ===
      do jj=1,2
        do ii=1,2
          eij = 0.d0
          do no=1,6
            eij = eij +0.5d0*(
     &          ue(ii,no)*dndx_m(jj,no)+ue(jj,no)*dndx_m(ii,no))
          enddo
          str(ii,jj) = eij
        enddo
      enddo
      emstr = (str(1,1)+str(2,2)+str(3,3))/3.d0
c
c   === Compute Averaged Volumetoric Stresses ===
      CALL stress_dp_bn(itrmax, idepg,
     &             prope,  sig,   str, ehist,
     &              ctol,  vons, e_dns, p_dns,
     &             ctens,
     &            ierror,  itr , histi )
      if(ierror.ne.0) RETURN
c
      smean = ( sig(1,1) +sig(2,2) +sig(3,3) )/3.d0
      sig = 0.d0 ! sig(:,:) = 0.d0
      do ii=1,3
        sig(ii,ii) = smean
c       sts(ii,ii) = smean
      enddo
c
c ***** Compute the Internal Force Vector ******************************
        do no=1,6
          do ii=1,ndf
            ndof = ndf*(no -1) +ii
            qol = 0.d0
            do jj=1,ndf
              qol = qol +sig(jj,ii)*dndx_m(jj,no)
            enddo
            finte(ndof) = finte(ndof) +qol *area
          enddo
        enddo
c
c **********************************************************************
c ***** Start Numerical Integration (Deviatoric part) ******************
c **********************************************************************
c     write(*,*) 'Deviatoric'
      do 200 ig=1,4
c
c     --- Set local displacements ---
        do ii=1,3
          jj=ijk_t(ii,ig)
          do ll=1,2
            uel(ll,ii)=ue(ll,jj)
          enddo
        enddo
c
        str = 0.d0 ! str(:,:) = 0.d0
        sig = 0.d0 ! sig(:,:) = 0.d0
c
c   === Area of this element ===
        det = det_g(ig,nel)
        detwxyz = 0.5d0*det
c
c ***** Compute Local Strain *******************************************
c   === Derivative of Shape Function r.w.t. Current Position ===
c       ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
        dndx(:,:) = dndx_g(:,:,ig,nel)
c
c   === Compute Local Strain Tensor ===
        do jj=1,2
          do ii=1,2
            eij = 0.d0
            do no=1,3
              eij = eij +0.5d0*(
     &            uel(ii,no)*dndx(jj,no)+uel(jj,no)*dndx(ii,no))
            enddo
            str(ii,jj) = eij
          enddo
        enddo
c
c   === Deviatoric Part ===
        em = ( str(1,1) +str(2,2) +str(3,3) )/3.d0
        do ii=1,3
         str(ii,ii) = str(ii,ii) -em
        enddo
c
c ***** Compute Local Deviatoric Stresses ******************************
c   === Set Deformation Histtories ===
        ehist(:) = dhist0(:,ig,nel)
c
c   === Compute Local Stresses (Drucker-Prager) ===
        CALL stress_dp_bn(itrmax, idepg,
     &               prope,  sig,   str, ehist,
     &                ctol,  vons, e_dns, p_dns,
     &               ctens,
     &              ierror,  itr , histi )
        if(ierror.ne.0) RETURN
c
c   === Clean up Volumetric Part ===
        sm = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
c
c ***** Store Deformation Histories (at current state) *****************
        dhist1(:,ig,nel) = ehist(:)
        ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
        idep(ig) = idepg
c
c ***** Compute the Elemental Value (Sum-up for Element ) **************
        alpeg = ehist(1)
        vone = vone +vons*detwxyz
        epse = epse +alpeg*detwxyz
        sts(:,:) = sts(:,:) +sig(:,:)*detwxyz
        stn(:,:) = stn(:,:) +str(:,:)*detwxyz
        e_ene = e_ene +e_dns*detwxyz
        p_ene = p_ene +p_dns*detwxyz
c
c ***** Compute the Internal Force Vector ******************************
        do no=1,3
          mp = ijk_t(no,ig)
          do ii=1,ndf
            ndof  =ndf*(mp -1) +ii
c           ndof = ndf*(no -1) +ii
            qol = 0.d0
            do jj=1,ndf
              qol = qol +sig(jj,ii)*dndx(jj,no)
            enddo
            finte(ndof) = finte(ndof) +qol *detwxyz
          enddo
        enddo
c
c **********************************************************************
  200 CONTINUE
c
c
c ***** Compute the Elemental Value (Volume Average) *******************
      vone = vone/area
      epse = epse/area

      tene_e = tene_e +e_ene
      dene_p = dene_p +p_ene
c
      e_ene=e_ene/area
      p_ene=p_ene/area
c
      sts = sts/area ! sts(:,:) = sts(:,:)/area
      stn = stn/area ! sts(:,:) = sts(:,:)/area
c
      do ii=1,3
        sts(ii,ii) = sts(ii,ii) +smean
      enddo
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
