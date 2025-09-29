      subroutine redcml(   neq,  neqm,
     &                     ijk,   mpe,  mang,  nfix,  ndir,
     &                  npload,nsload,nbload,  ijkl, melem,
     &                     ncp,  mdof,  idof, mpstp, muprt,
     &                   mfprt, mnprt, msprt, mbprt, matid,
     &                  arcini,
     &                     xyz,  prop, angle,  vfix,vpload,
     &                  vsload,vbload,  finc,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension mpe(nelx),mang(nelx),melem(nelx)
      dimension ijk(node,nelx)
      dimension nfix(nspc),ndir(6,nspc)
      dimension npload(npoin)
      dimension nsload(npres)
      dimension nbload(nbody)
      dimension ijkl(nsn,npres)
      dimension ncp(neq),mdof(neq),idof(neq)
      dimension mpstp(lpstp)
      dimension muprt(2,luprt)
      dimension mfprt(2,lfprt)
      dimension mnprt(2,lnprt)
      dimension msprt(2,lsprt)
      dimension mbprt(2,lbprt)
      dimension matid(lmat)
c
      dimension xyz(3,nx)
      dimension prop(20,lmat)
      dimension angle(3,lang)
      dimension vfix(6,nspc)
      dimension vpload(6,npoin)
      dimension vsload(4,npres)
      dimension vbload(3,nbody)
      dimension finc(nstep)
c
      character*80 title
      character*7  charin
      character*7  keyLAS,keyEND,keySOL,
     &             keyPST,keyPDI,keyPFO,keyPSN,keyPSS,keyPTE,
     &             key000,key010,
     &             key100,key233,key236,key244,key248,key264,key268,
     &             key500,key550,key600,key610,key620,key650
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
      common /print/ lpstp,luprt,lfprt,lnprt,lsprt,lbprt
c
c **********************************************************************
      REWIND lra
      WRITE(*,100)
c
c ***** Initialization *************************************************
      ijk = 0 ! ijk(:,:) = 0
      ncp = 0 ! ncp(:) = 0
c
c ***** Definition of CML-formatted Formats ****************************
  100 FORMAT(/,5x,'READing the FEM data ',$)
  101 FORMAT('.',$)
  102 FORMAT('. done!!',/)
c
 9000 FORMAT(a7)
 9001 FORMAT(a80)
c
 9201 FORMAT(i8,3e15.5)
 9202 FORMAT(i8,4i5,8i8)
 9203 FORMAT(i8,3i5,8i8)
 9204 FORMAT(5e12.5)
 9205 FORMAT(i5,3e13.5)
 9206 FORMAT(i8,i5,i2,5i1,1p6e12.5)
 9207 FORMAT(i8,6e12.5)
 9208 FORMAT(3i8,4e12.5)
 9209 FORMAT(5i8,4e12.5)
 9210 FORMAT(i8,4e12.5)
c
 9500 FORMAT(5i8)
 9501 FORMAT(i8,i5)
 9502 FORMAT(5i5)
 9503 FORMAT(i5,i8,i5)
 9504 FORMAT(f12.5)
 9505 FORMAT(2i5,f12.5)
 9506 FORMAT(i5,i8)
c
c ***** Definition of CML-formatted Headers ****************************
      key000 = '/TITLE/'
      key010 = '/SOLVR/'
      keyLAS = '/LASTD/'
      keyEND = '/ENDOF/'
      keySOL = '/SOLUT/'
      keyPST = '/PSTEP/'
      keyPDI = '/PDISP/'
      keyPFO = '/PFOCE/'
      keyPSN = '/PSTRN/'
      keyPSS = '/PSTRS/'
      keyPTE = '/PTEMP/'
c
      key100 = '/COORD/'
      key233 = '/TRIA3/'
      key236 = '/TRIA6/'
      key244 = '/QUAD4/'
      key248 = '/QUAD8/'
      key264 = '/TTRA4/'
      key268 = '/HEXA8/'
c
      key500 = '/MATER/'
      key550 = '/EULER/'
      key600 = '/CONST/'
      key610 = '/PCNST/'
      key620 = '/INFTY/'
      key650 = '/LOADC/'
c
c ***** Read Some Data of the FE-model *********************************
 1000 CONTINUE
      READ(lra,9000,err=9998, end=9999) charin
c
c /TITLE/ ***** Title of this FE model *****
      if(charin.eq.key000) then
      WRITE(*,101)
      READ(lra,*) title
c
c /COORD/ ***** Coordinates of each node *****
      elseif(charin.eq.key100) then
        WRITE(*,101)
        READ(lra,9500) ndum
        do nn=1,nx
          READ(lra,9201) ndum,(xyz(kk,nn), kk=1,3)
        enddo
c
c /TRIA3/ ***** 3-nodes TRIANGLE element (2D) *****
      elseif(charin.eq.key233) then
        WRITE(*,101)
        READ(lra,9501) nel_t3,ndum
        do nel=1,nel_t3
          READ(lra,9202) kk,mpe(kk),mang(kk),ndum,ndum,
     &                   (ijk(ll,kk),ll=1,3)
          melem(kk) = 23
        enddo
c
c /TRIA6/ ***** 6-nodes TRIANGLE element (2D) *****
      elseif(charin.eq.key236) then
        WRITE(*,101)
        READ(lra,9501) nel_t6,ndum
        do nel=1,nel_t6
          READ(lra,9202) kk,mpe(kk),mang(kk),ndum,ndum,
     &                   (ijk(ll,kk),ll=1,6)
          melem(kk) = 26
        enddo
c
c /QUAD4/ ***** 4-nodes Isoparametric element (2D) *****
      elseif(charin.eq.key244) then
        WRITE(*,101)
        READ(lra,9501) nel_q4,ndum
        do nel=1,nel_q4
          READ(lra,9202) kk,mpe(kk),mang(kk),ndum,ndum,
     &                   (ijk(ll,kk),ll=1,4)
          melem(kk) = 24
        enddo
c
c /QUAD8/ ***** 8-nodes Isoparametric element (2D) *****
      elseif(charin.eq.key248) then
        WRITE(*,101)
        READ(lra,9501) nel_q8,ndum
        do nel=1,nel_q8
          READ(lra,9202) kk,mpe(kk),mang(kk),ndum,ndum,
     &                   (ijk(ll,kk),ll=1,8)
          melem(kk) = 28
        enddo
c
c /TTRA4/ ***** 4-nodes Tetrahedral element (3D) *****
      elseif(charin.eq.key264) then
        WRITE(*,101)
        READ(lra,9501) nel_t4,ndum
        do nel=1,nel_t4
          READ(lra,9203) kk,mpe(kk),mang(kk),ndum,
     &                   (ijk(ll,kk),ll=1,4)
          melem(kk) = 34
        enddo
c
c
c /HEXA8/ ***** 8-nodes Hexahedral element (3D) *****
      elseif(charin.eq.key268) then
        WRITE(*,101)
        READ(lra,9501) nel_h8,ndum
        do nel=1,nel_h8
          READ(lra,9203) kk,mpe(kk),mang(kk),ndum,
     &                   (ijk(ll,kk),ll=1,8)
          melem(kk) = 38
        enddo
c
c /MATER/ ***** Material Properties *****
      elseif(charin.eq.key500) then
        WRITE(*,101)
        READ(lra,9502) ndum
        do lm=1,lmat
          READ(lra,9502) nn,matid(lm)
          READ(lra,9204) (prop(kk,nn),kk= 1,20)
        enddo
c
c /EULER/ ***** Euler angle for material *****
      elseif(charin.eq.key550) then
        WRITE(*,101)
        READ(lra,9502) ndum
        do la=1,lang
          READ(lra,9205) nn,(angle(kk,la),kk=1,3)
        enddo
c
c /CONST/ ***** 1st Boundary Condition: Constraint *****
      elseif(charin.eq.key600) then
        WRITE(*,101)
        READ(lra,9502) ndum,ndum,ndum
c
c       --- Single Point Constraint ---
        do ns=1,nspc
          READ(lra,9206) nfix(ns),mm,
     &                   (ndir(kk,ns),kk=1,6),
     &                   (vfix(kk,ns),kk=1,6)
        enddo
c
c /LOADC/ ***** 2nd Boundary Condition: Loading *****
      elseif(charin.eq.key650) then
        WRITE(*,101)
        READ(lra,9502) ndum
        do l=1,lnum
          READ(lra,9502) ndum,ndum,ndum
c
c         --- Point Load ---
          do np=1,npoin
            READ(lra,9207) npload(np),(vpload(kk,np),kk=1,6)
          enddo
c
c         --- Distributed (surface) Load ---
          if(ndf.eq.2) then
            do np=1,npres
              READ(lra,9208) nsload(np),
     &                       (ijkl(kk,np),kk=1,nsn),
     &                       (vsload(kk,np),kk=1,4)
            enddo
          elseif(ndf.eq.3) then
            do np=1,npres
              READ(lra,9209) nsload(np),
     &                       (ijkl(kk,np),kk=1,nsn),
     &                       (vsload(kk,np),kk=1,4)
            enddo
          endif
c
c         --- Body Force ---
          do nb=1,nbody
            READ(lra,9210) nbload(nb),dum,(vbload(kk,np),kk=1,3)
          enddo
c
        enddo
c
c /SOLVR/ ***** Solvers & Their Options *****
      elseif(charin.eq.key010) then
      WRITE(*,101)
      READ(lra,*)
      READ(lra,*)
c
c /SOLUT/ ***** Solution Procedures *****
      elseif(charin.eq.keySOL) then
        WRITE(*,101)
        READ(lra,9502) nstep,mdum,mdum,istart,iarc
        READ(lra,9502) mdum,itrobj,mdum,mdum,mdum
c    === for Load/Displacement Control ===
        if(iarc.eq.0) then
c       --- Load increments are defined EXPLICITLY
          if(itrobj.ge.1) then
            do it=1,itrobj
              READ(lra,9505) mst,med,df
              do mm=mst,med
                finc(mm) = df
              enddo
            enddo
c       --- EQUAL Load increments
          else
            do ns=1,nstep
              finc(ns) = 1.d0/dble(nstep)
            enddo
          endif
c    === for Arc-Length Method ===
        else
          READ(lra,9504) arcini
        endif
c
c /PSTEP/ ***** Print control section *****
      elseif(charin.eq.keyPST) then
        WRITE(*,101)
        READ(lra,9502) mm
        if(mm.eq.0) then
          WRITE(*,640)
          mm = -1
        endif
c
        if(mm.lt.0) then
          nsect = -mm
          do lp=1,lpstp
            mpstp(lp) = nsect*lp
          enddo
        elseif(mm.gt.0) then
          do lp=1,lpstp
            READ(lra,*) mdum,mpstp(lp)
          enddo
        endif
c
c /PDISP/ ***** ID & NDOF where nodal displacements are needed *****
      elseif(charin.eq.keyPDI) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,luprt
          READ(lra,9503) mdum,muprt(1,nn),muprt(2,nn)
        enddo
c
c /PFOCE/ ***** ID & NDOF where nodal forces are needed *****
      elseif(charin.eq.keyPFO) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,lfprt
          READ(lra,9503) mdum,mfprt(1,nn),mfprt(2,nn)
        enddo
c
c /PSTRN/ ***** ID & DOF where elemental STRAIN are needed *****
      elseif(charin.eq.keyPSN) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,lnprt
          READ(lra,9503) mdum,mnprt(1,nn),mnprt(2,nn)
        enddo
c
c /PSTRS/ ***** ID & DOF where elemental STRESS are needed *****
      elseif(charin.eq.keyPSS) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,lsprt
          READ(lra,9503) mdum,msprt(1,nn),msprt(2,nn)
        enddo
c
c /PTEMP/ ***** ID where elemental TEMPERATURE are needed *****
       elseif(charin.eq.keyPTE) then
         WRITE(*,101)
         read(lra,*)
         do nn=1,lbprt
           read(lra,9506) mdum,mbprt(1,nn)
         enddo
c
c /PENER/ ***** ID where elemental Energy are needed *****
c       elseif(charin.eq.keyPEN) then
c         write(*,101)
c         read(lra,*)
c         do nn=1,laprt
c           read(lra,9502) mdum,maprt(1,nn)
c         enddo
c
c /LASTD/ ***** Last of input data *****
      elseif(charin.eq.keyLAS) then
        WRITE(*,101)
c
c /ENDOF/ ***** End of this file *****
      elseif(charin.eq.keyEND) then
        WRITE(*,102)
c       WRITE(*,101)
        GOTO 2000
c
c /*****/ ***** NOT Supported Indicator ****
      else
        WRITE(*,*) charin
      endif
c
      GOTO 1000
c
c ***** Error Section **************************************************
 9998 CONTINUE
      WRITE(*,*) 'Something WRONG in your CML-format!!'
      ierror = 8
      RETURN
c
 9999 CONTINUE
      WRITE(*,*) 'You forgot /ENDOF/ indicator in your CML-file'
      ierror = 9
      RETURN
c
c **********************************************************************
 2000 CONTINUE
c
c ****** Set ON/OFF Flags for Constraint *******************************
c     --- constraint for displacement value (single point const.)
      do ns=1,nspc
        do nd=1,ndf
          if(ndir(nd,ns).eq.1) then
            nd0 = ndf*(nfix(ns) -1)
            ncp(nd0 +nd) = 1
          endif
        enddo
      enddo
c
c ****** Renumbering the DOF *******************************************
      kk = 0
      do ne=1,neq
        if(ncp(ne).eq.0) then
          kk = kk +1
          mdof(ne) = kk
          idof(kk) = ne
        endif
      enddo
      neqm = kk
c
      do ne=1,neq
        if(ncp(ne).ne.0) then
          kk = kk +1
          mdof(ne) = kk
          idof(kk) = ne
        endif
      enddo
c
c **********************************************************************
  640 FORMAT(/,3x,
     &       '***** Parameter "lpstp" is NOT defined *****',
     &       /,10x,'Default Number is used "lpstp = 1" ',/)
c **********************************************************************
c **********************************************************************
      RETURN
      END
