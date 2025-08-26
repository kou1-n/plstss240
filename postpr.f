      subroutine postpr(   neq,
     *                     ijk,   mpe,  mang, idep1, melem,
     &                   matid,
     &                     xyz,  prop, angle,    u1, sigma,
     &                   epsln,   von,  fint,dhist0,dhist1,
     &                     eps,  pene,  eene,
     &                  tene_e,dene_p,  temp, dtemp, tempd,
     &                  dndx_g, det_g,ctensg, g_norm,
     &                  ierror, itr , histi0  )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension mpe(nelx),mang(nelx),melem(nelx)
      dimension idep1(ngaus,nelx)
      dimension matid(lmat)
c
      dimension xyz(3,nx)
      dimension prop(20,lmat)
      dimension angle(3,lang)
      dimension u1(neq),fint(neq)
      dimension sigma(6,nelx),epsln(6,nelx)
      dimension von(nelx),eps(nelx),pene(nelx),eene(nelx),temp(nelx),
     &          tempd(nelx)
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
      dimension psogg(3,3,ngaus,nelx)
      dimension dhist0(20,ngaus,nelx),dhist1(20,ngaus,nelx)
      dimension histi0(50,ngaus,nelx)
c
      dimension ijke(8),idep(ngaus)
c     dimension ijke(8),idep(8)
c
      dimension sts(3,3),stn(3,3)
      dimension xe(3,8),ue(3,8)
      dimension prope(20)
      dimension finte(24)
      dimension g_vals(ngaus)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
c
c **********************************************************************
c
c ****** Initialization ************************************************
      sigma = 0.d0 ! sigma(:,:) = 0.d0
      epsln = 0.d0 ! epsln(:,:) = 0.d0
c
      fint = 0.d0 ! fint(:) = 0.d0
c
      tene_e = 0.d0
      dtemp = 0.d0
c     dtempe = 0.d0
      dene_p = 0.d0
c
c     Initialize yield residual norm
      g_norm = 0.d0
c
c ****** Loop for ALL Elements *****************************************
      do 100 nel=1,nelx
c
c ===== Define Elelental Basic Data =====
c     --- Material Properties ---
        mat = mpe(nel)
        man = mang(nel)
        prope(:) = prop(:,mat)
        MATYPE = matid(mat)
c
c ===== Computation of Elemental Internal Forces =====
c     --- /QUAD4/ 4-node Isoparametric Element ---
        if( melem(nel).eq.24 ) then
          do no=1,4
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            ue(1,no) = u1(ndf*(ijkia-1) +1)
            ue(2,no) = u1(ndf*(ijkia-1) +2)
          enddo
c
          CALL pquad4(   nel,  nelx,   ndf,  node, ngaus,
     &                  idep, prope,    xe,    ue,
     &                   sts,   stn, finte,  vone,  epse,
     &                dhist0,dhist1,
     &                 p_ene, e_ene,tene_e,dene_p, tempe,
     &                 dtemp,dtempe,
     &                dndx_g, det_g,ctensg,
     &                ierror )
          if(ierror.ne.0) RETURN
c
c     --- /TRIA3/ 3-node Constant Strain Triangle
        elseif( melem(nel).eq.23 ) then
          do no=1,3
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            ue(1,no) = u1(ndf*(ijkia-1) +1)
            ue(2,no) = u1(ndf*(ijkia-1) +2)
          enddo
c
          CALL ptria3(   nel,  nelx,   ndf,  node, ngaus,
     &                  idep, prope,    xe,    ue,
     &                   sts,   stn, finte,  vone,  epse,
     &                dhist0,dhist1,
     &                 p_ene, e_ene,tene_e,dene_p, tempe,
     &                 dtemp,dtempe,
     &                dndx_g, det_g,ctensg,
     &                ierror )
          if(ierror.ne.0) RETURN
c
c     --- /QUAD8/ 8-node Isoparametric Element (2nd order)
        elseif( melem(nel).eq.28 ) then
          do no=1,8
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            ue(1,no) = u1(ndf*(ijkia-1) +1)
            ue(2,no) = u1(ndf*(ijkia-1) +2)
          enddo
c
          CALL pquad8(   nel,  nelx,   ndf,  node, ngaus,
     &                  idep, prope,    xe,    ue,
     &                   sts,   stn, finte,  vone,  epse,
     &                dhist0,dhist1,
     &                 p_ene, e_ene,tene_e,dene_p, tempe,
     &                 dtemp,dtempe,
     &                dndx_g, det_g,ctensg,
     &                ierror )
          if(ierror.ne.0) RETURN
c
c     --- /TRIA6/ 6-node Isoparametric Element (2nd order)
        elseif( melem(nel).eq.26 ) then
          do no=1,6
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            ue(1,no) = u1(ndf*(ijkia-1) +1)
            ue(2,no) = u1(ndf*(ijkia-1) +2)
          enddo
c
          if(incomp.eq.0) then
c         === Usual 6-node parabolic element ===
            CALL ptria6(   nel,  nelx,   ndf,  node, ngaus,
     &                    idep, prope,    xe,    ue,
     &                     sts,   stn, finte,  vone,  epse,
     &                  dhist0,dhist1,
     &                   p_ene, e_ene,tene_e,dene_p, tempe,
     &                   dtemp,dtempe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
c
          elseif(incomp.eq.2) then
c         === Usual 6-node parabolic element ===
            CALL ptrip1(   nel,  nelx,   ndf,  node, ngaus,
     &                    idep, prope,    xe,    ue,
     &                     sts,   stn, finte,  vone,  epse,
     &                  dhist0,dhist1,
     &                   p_ene, e_ene,tene_e,dene_p, tempe,
     &                   dtemp,dtempe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
          endif
          if(ierror.ne.0) RETURN
c
c     --- /HEXA8/ 8-node Isoparametric Element ---
        elseif( melem(nel).eq.38 ) then
          do no=1,8
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            xe(3,no) = xyz(3,ijkia)
            ue(1,no) = u1(ndf*(ijkia-1) +1)
            ue(2,no) = u1(ndf*(ijkia-1) +2)
            ue(3,no) = u1(ndf*(ijkia-1) +3)
          enddo
c
          CALL phexa8(   nel,  nelx,   ndf,  node, ngaus,
     &                  idep,MATYPE, prope,    xe,    ue,
     &                   sts,   stn, finte,  vone,  epse,
     &                dhist0,dhist1,
     &                 p_ene, e_ene,tene_e,dene_p, tempe,
     &                 dtemp,dtempe,
     &                dndx_g, det_g,ctensg, g_vals,
     &                ierror, itr , histi0 )
          if(ierror.ne.0) RETURN
c
c       --- Accumulate yield residual norm for BN method ---
          if(MATYPE.eq.5) then
            do ig=1,ngaus
              g_norm = g_norm + g_vals(ig)**2
            enddo
          endif
c
c     --- /TTRA4/ 4-node Constant Strain Triangle
        elseif( melem(nel).eq.34 ) then
c
        endif
c
c ===== Compute Internal Force Vector =====
        do no=1,node
          ijkia = ijk(no,nel)
          if(ijkia.ne.0) then
            do ii=1,ndf
              ndof = ndf*(ijkia -1) +ii
              nd   = ndf*(no -1) +ii
              fint(ndof) = fint(ndof) +finte(nd)
            enddo
          endif
        enddo
c       do nn=1,nx
c         write(*,*) nn,(fint(ndf*(nn-1)+jj),jj=1,2)
c       enddo
c
c ===== Store Deformation Histories =====
        do ig=1,ngaus
          idep1(ig,nel) = idep(ig)
        enddo
c
c
c ===== Store Elemental Data =====
c     --- Stress ---
        sigma(1,nel) = sts(1,1)
        sigma(2,nel) = sts(2,1)
        sigma(3,nel) = sts(3,1)
        sigma(4,nel) = sts(2,2)
        sigma(5,nel) = sts(3,2)
        sigma(6,nel) = sts(3,3)
c
        von(nel) = vone
        eps(nel) = epse
c
c     --- Strain ---
        epsln(1,nel) = stn(1,1)
        epsln(2,nel) = stn(2,1)
        epsln(3,nel) = stn(3,1)
        epsln(4,nel) = stn(2,2)
        epsln(5,nel) = stn(3,2)
        epsln(6,nel) = stn(3,3)
c
c     --- Energy ---
        pene(nel) = p_ene
        eene(nel) = e_ene
c
c       temp(nel) = tempe
        tempd(nel) = dtempe
c       write(*,*) dtempe
c
  100 continue
c
c ***** Compute the L2 norm of yield residuals for BN method *********
      if(g_norm.gt.0.d0) then
        g_norm = dsqrt(g_norm)
      endif
c
c **********************************************************************
c     do ne=1,ndf*nx
c       write(*,'(i5,3e15.5)') ne,fint(ne)
c     enddo
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
