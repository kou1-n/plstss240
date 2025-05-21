      subroutine initia(  NGSK,   neq,  neqm, NGSKo,  neqo,
     &                     ijk,   mpe,  mang,  mdof,  idof,
     &                   jdiag,index0,index1,jcolmn, melem,
     &                     xyz,  prop, angle,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension mpe(nelx),mang(nelx),melem(nelx)
      dimension mdof(neq),idof(neq)
      dimension jdiag(neq+1)
      dimension index0(neqo+1)
      dimension index1(NGSKo)
      dimension jcolmn(nsolvr+1)
c
      dimension xyz(3,nx)
      dimension prop(20,lmat)
      dimension angle(3,lang)
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
c
      dimension ijke(8)
c
      dimension xe(3,8)
      dimension prope(20)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
c
c **********************************************************************
      NGSKm = NGSK -NGSKo
c
c ****** Loop for All Elements *****************************************
      do 100 nel=1,nelx
c
c ===== Define Elemental Basic Data =====
c     --- Material Properties ---
        mat = mpe(nel)
        man = mang(nel)
        prope(:) = prop(:,mat)
c
c ===== Construct Elemental Stiffness Matrix "ske" =====
c     --- /QUAD4/ 4-node Isoparametric Element
        if( melem(nel).eq.24 ) then
          do no=1,4
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
          enddo
c
          CALL initq4(   nel,  nelx,   ndf,  node, ngaus,
     &                 prope,  ijke,    xe,
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
          enddo
c
          CALL initt3(   nel,  nelx,   ndf,  node, ngaus,
     &                 prope,  ijke,    xe,
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
          enddo
c
          CALL initq8(   nel,  nelx,   ndf,  node, ngaus,
     &                 prope,  ijke,    xe,
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
          enddo
c
c       === Usual 6-node parabolic element ===
          if(incomp.eq.0) then
            CALL initt6(   nel,  nelx,   ndf,  node, ngaus,
     &                   prope,  ijke,    xe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
c
          elseif(incomp.eq.2) then
c       === P1-iso-P2/P0 element ===
            CALL initp1(   nel,  nelx,   ndf,  node, ngaus,
     &                   prope,  ijke,    xe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror )
          endif
          if(ierror.ne.0) RETURN
c
c     --- /HEXA8/ 8-node Isoparametric Element
        elseif( melem(nel).eq.38 ) then
          do no=1,8
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            xe(3,no) = xyz(3,ijkia)
          enddo
c
          CALL inith8(   nel,  nelx,   ndf,  node, ngaus,
     &                 prope,  ijke,    xe,
     &                dndx_g, det_g,ctensg,
     &                ierror )
          if(ierror.ne.0) RETURN
c
c     --- /TTRA4/ 4-node Constant Strain Triangle
        elseif( melem(nel).eq.34 ) then
          do no=1,4
            ijkia = ijk(no,nel)
            ijke(no) = ijkia
            xe(1,no) = xyz(1,ijkia)
            xe(2,no) = xyz(2,ijkia)
            xe(3,no) = xyz(3,ijkia)
          enddo
c
          CALL initt4(   nel,  nelx,   ndf,  node, ngaus,
     &                 prope,  ijke,    xe,
     &                dndx_g, det_g,ctensg,
     &                ierror )
          if(ierror.ne.0) RETURN
c
        endif
c
  100 continue
c     stop 'after initia'
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
