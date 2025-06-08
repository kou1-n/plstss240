      subroutine forces(   neq,
     &                     ijk,  nfix,  ndir,
     &                  npload,nsload,  ijkl, 
     &                     xyz,  vfix,vpload,
     &                  vsload,   foc,  disp,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension nfix(nspc)
      dimension ndir(6,nspc)
      dimension npload(npoin)
      dimension nsload(npres)
      dimension ijkl(nsn,npres)
c
      dimension xyz(3,nx)
      dimension vfix(6,nspc)
      dimension vpload(6,npoin)
      dimension vsload(4,npres)
      dimension foc(neq),disp(neq)
c
      dimension ijk_e(8)
      dimension ijkl_e(4)
c
      dimension xe(3,8)
      dimension fe(32)
      dimension vsl_e(4)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
c
c **********************************************************************
c
c ****** Initialization ************************************************
      foc = 0.d0 ! foc(:) = 0.d0
      disp = 0.d0 ! disp(:) = 0.d0
c
c ****** Mechanical External Forces (Distributed Force) ****************
      do 100 np=1,npres
c
c       --- initialization ---
        ijk_e = 0 ! ijk_e(:) = 0
        ijkl_e = 0 ! ijkl_e(:) = 0
c
        fe = 0.d0 ! fe(:) = 0.d0
        xe = 0.d0 ! xe(:,:) = 0.d0
        vsl_e = 0.d0 ! vsl_e(:) = 0.d0
c
c       --- elemental basic settings ---
        nel = nsload(np)
        do ia=1,node
          ij = ijk(ia,nel)
          ijk_e(ia) = ij
          if(ij.ne.0) then
            do kk=1,3
              xe(kk,ia) = xyz(kk,ij)
            enddo
            nod_el = ia
          endif
        enddo
c
        do ns=1,nsn
          ijkl_e(ns) = ijkl(ns,np)
        enddo
c
        do ii=1,4
          vsl_e(ii) = vsload(ii,np)
        enddo
c
c       --- compute equivalent nodal force ---
c     ( for 2-D problem )
        if(ndf.eq.2) then
          CALL pres2d(   nel, ijk_e,ijkl_e,
     &                    xe,    fe, vsl_e,
     &                incomp,
     &                ierror )
c
c     ( for 3-D problem )
        elseif(ndf.eq.3) then
          CALL pres3d(   nel, ijk_e,ijkl_e,
     &                    xe,    fe, vsl_e,
     &                ierror )
        endif
c
c       --- assemble global force vector ---
        do ia=1,node
          ij = ijk(ia,nel)
          if(ij.ne.0) then
            nd0 = ndf*(ij -1)
            nd1 = ndf*(ia -1)
            do nd=1,ndf
              niai = nd0 +nd
              iai  = nd1 +nd
              foc(niai) = foc(niai) +fe(iai)
            enddo
          endif
        enddo
c
  100 continue
c
c ****** Mechanical External Forces (Point Force) **********************
      do np=1,npoin
        npl = npload(np)
        nd0 = ndf*(npl -1)
        do nd=1,ndf
          niai = nd0 +nd
          foc(niai) = foc(niai) +vpload(nd,np)
        enddo
      enddo
c
c ****** Enforced Displacement *****************************************
      do ns=1,nspc
        nf = nfix(ns)
        nd0 = ndf*(nf -1)
        do nd=1,ndf
          if(ndir(nd,ns).eq.1) then
            ne = nd0 +nd
            disp(ne) = vfix(nd,ns)
          endif
        enddo
      enddo
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
