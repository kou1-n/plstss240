      subroutine estima(   neq,  neqm,  NGSK, NGSKm,  iRCM,
     &                     ijk,   ncp, jdiag,  mdof,  idof,
     &                  iw_neq,  list,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension ncp(neq),mdof(neq),idof(neq),iw_neq(neq)
      dimension jdiag(neq+1)
      dimension list(nmax,nx)
c
      dimension ndir(6)
c
      character*80 title
      character*7  charin
      character*7  key000,keyLAS,keyEND,keySOL,
     &             keyPST,keyPDI,keyPFO,keyPSN,keyPSS,keyPTE,
     &             key010,
     &             key100,key233,key236,key244,key248,key264,key268,
     &             key500,key550,key600,key610,key620,key650
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
      common /print/ lpstp,luprt,lfprt,lnprt,lsprt,lbprt
c
c **********************************************************************
      REWIND lra
c
c ***** Definition of CML-formatted Formats ****************************
  101 FORMAT('.',$)
  102 FORMAT('. done!!')
c
 9000 FORMAT(a7)
 9001 FORMAT(a80)
 9202 FORMAT(i8,4i5,8i8)
 9203 FORMAT(i8,3i5,8i8)
 9206 FORMAT(i8,i5,i2,5i1,1p6e12.5)
 9500 FORMAT(5i8)
 9501 FORMAT(i8,i5)
 9502 FORMAT(5i5)
 9503 FORMAT(i5,i8,i5)
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
c ***** Initialization of Arrays ***************************************
      ijk = 0 ! ijk(:,:) = 0
c
      ncp = 0 ! ncp(:) = 0
      jdiag = 0 ! jdiag(:) = 0
      mdof = 0 ! mdof(:) = 0
      idof = 0 ! idof(:) = 0
c
c ***** Read Necessary Data from CML-formatted File ********************
 1000 CONTINUE
      READ(lra,9000,err=9998, end=9999) charin
c
c /TITLE/ ***** Title of this FE model *****
      if(charin.eq.key000) then
      WRITE(*,101)
      READ(lra,*)
c
c /COORD/ ***** Coordinates of each node *****
      elseif(charin.eq.key100) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,nx
          READ(lra,*)
        enddo
c
c /TRIA3/ ***** 3-nodes TRIANGLE element (2D) *****
      elseif(charin.eq.key233) then
        WRITE(*,101)
        READ(lra,9501) nel_t3,ndum
        do nel=1,nel_t3
          READ(lra,9202) nn,ndum,ndum,ndum,ndum,
     &                   (ijk(kk,nn),kk=1,3)
        enddo
c
c /TRIA6/ ***** 6-nodes TRIANGLE element (2D) *****
      elseif(charin.eq.key236) then
        WRITE(*,101)
        READ(lra,9501) nel_t6,ndum
        do nel=1,nel_t6
          READ(lra,9202) nn,ndum,ndum,ndum,ndum,
     &                   (ijk(kk,nn),kk=1,6)
        enddo
c
c /QUAD4/ ***** 4-nodes Isoparametric element (2D) *****
      elseif(charin.eq.key244) then
        WRITE(*,101)
        READ(lra,9501) nel_q4,ndum
        do nel=1,nel_q4
          READ(lra,9202) nn,ndum,ndum,ndum,ndum,
     &                   (ijk(kk,nn),kk=1,4)
        enddo
c
c /QUAD8/ ***** 8-nodes Isoparametric element (2D) *****
      elseif(charin.eq.key248) then
        WRITE(*,101)
        READ(lra,9501) nel_q8,ndum
        do nel=1,nel_q8
          READ(lra,9202) nn,ndum,ndum,ndum,ndum,
     &                   (ijk(kk,nn),kk=1,8)
        enddo
c
c /TTRA4/ ***** 4-nodes Tetrahedral element (3D) *****
      elseif(charin.eq.key264) then
        WRITE(*,101)
        READ(lra,9501) nel_t4,ndum
        do nel=1,nel_t4
          READ(lra,9203) nn,ndum,ndum,ndum,
     &                   (ijk(kk,nn),kk=1,4)
        enddo
c
c /HEXA8/ ***** 8-nodes Hexahedral element (3D) *****
      elseif(charin.eq.key268) then
        WRITE(*,101)
        READ(lra,9501) nel_h8,ndum
        do nel=1,nel_h8
          READ(lra,9203) nn,ndum,ndum,ndum,
     &                   (ijk(kk,nn),kk=1,8)
        enddo
c
c /MATER/ ***** Material Properties *****
      elseif(charin.eq.key500) then
        WRITE(*,101)
        READ(lra,*)
        do lm=1,lmat
          READ(lra,*)
          do ii=1,4
            READ(lra,*)
          enddo
        enddo
c
c /EULER/ ***** Euler angles for materials *****
      elseif(charin.eq.key550) then
        WRITE(*,101)
        READ(lra,*)
        do la=1,lang
          READ(lra,*)
        enddo
c
c /CONST/ ***** 1st boundary condition: Constraint *****
      elseif(charin.eq.key600) then
        WRITE(*,101)
        READ(lra,*)
        do lo=1,locc
          READ(lra,*)
        enddo
c       --- In the case of "Single Point Constraint"
        do ns=1,nspc
          READ(lra,9206) nf,mm,(ndir(k),k=1,6),
     &                   dum,dum,dum,dum,dum,dum
          nd0 = ndf*(nf-1)
          do nd=1,ndf
            if(ndir(nd).eq.1) then
              ncp(nd0 +nd) = 1
            endif
          enddo
        enddo
        do mp=1,mpc
          READ(lra,*)
        enddo
c
c /LOADC/ ***** 2nd boundary condition: Loading condition *****
      elseif(charin.eq.key650) then
        WRITE(*,101)
        READ(lra,*)
        do ln=1,lnum
          READ(lra,9502) npoin,npres,nbody
          do np=1,npoin
            READ(lra,*)
          enddo
          do np=1,npres
            READ(lra,*)
          enddo
          do nb=1,nbody
            READ(lra,*)
          enddo
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
        READ(lra,9502) nstep,mdum,mdum,mdum,iarc
        READ(lra,9502) mdum,itrobj,mdum,mdum,mdum
c    === for Load/Displacement Control ===
        if(iarc.eq.0) then
c       --- Load increments are defined EXPLICITLY
          if(itrobj.ge.1) then
            do it=1,itrobj
              READ(lra,*)
            enddo
          endif
c    === for Arc-Length Method ===
        else
          READ(lra,*)
        endif
c
c /PSTEP/ ***** Print control section *****
      elseif(charin.eq.keyPST) then
        WRITE(*,101)
        READ(lra,*) mm
        if(mm.gt.0) then
          do nn=1,lpstp
            READ(lra,*)
          enddo
        endif
c
c /PDISP/ ***** ID & DOF where nodal displacements are needed *****
      elseif(charin.eq.keyPDI) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,luprt
          READ(lra,*)
        enddo
c
c /PFOCE/ ***** ID & NDOF where nodal forces are needed *****
      elseif(charin.eq.keyPFO) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,lfprt
          READ(lra,*)
        enddo
c
c /PSTRN/ ***** ID & DOF where elemental STRAIN are needed *****
      elseif(charin.eq.keyPSN) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,lnprt
          READ(lra,*)
        enddo
c
c /PSTRS/ ***** ID & DOF where elemental STRESS are needed *****
      elseif(charin.eq.keyPSS) then
        WRITE(*,101)
        READ(lra,*)
        do nn=1,lsprt
          READ(lra,*)
        enddo
c
c /PTEMP/ ***** ID where elemental Temperature are needed *****
       elseif(charin.eq.keyPTE) then
       WRITE(*,101)
       read(lra,*)
       do nn=1,lbprt
         read(lra,*)
       enddo
c
c /PENER/ ***** ID where elemental Energy are needed *****
c       elseif(charin.eq.keyPEN) then
c         write(*,101)
c         read(lra,*)
c         do nn=1,laprt
c           read(lra,*)
c         enddo
c
c /LASTD/ ***** Last of input data *****
      elseif(charin.eq.keyLAS) then
c
c /ENDOF/ ***** End of this file *****
      elseif(charin.eq.keyEND) then
        WRITE(*,102)
        GOTO 2000
c
      else
        WRITE(*,*) charin
      endif
c
      GOTO 1000
c
c ***** Error Section **************************************************
 9998 CONTINUE
      WRITE(*,*) 'Something WRONG in your CML-format!!'
      ierror = 5
      RETURN
c
 9999 CONTINUE
      WRITE(*,*) 'You forgot /ENDOF/ indicator in your CML-file'
      ierror = 6
      RETURN
c
c **********************************************************************
 2000 CONTINUE
c
c     do nel=1,nelx
c       write(*,*) nel,(ijk(ii,nel),ii=1,node)
c     enddo
c
c     do nn=1,neq
c       write(*,*) nn,ncp(nn)
c     enddo
c
c ***** Renumbering the DOF ********************************************
c     ===== Default Numbering =====
      if(iRCM.eq.0) then
        kk = 0
        do ne=1,neq
          if(ncp(ne).eq.0) then
            kk = kk +1
            mdof(ne) = kk
            idof(kk) = ne
          endif
        enddo
        neqm = kk
        do ne=1,neq
          if(ncp(ne).ne.0) then
            kk = kk +1
            mdof(ne) = kk
            idof(kk) = ne
          endif
        enddo
c
c     ===== Use RNM_***.cml File =====
      else
        READ(lrb,*) 
        READ(lrb,*) 
        READ(lrb,*) 
        READ(lrb,*) 
        READ(lrb,'(2i8)') nn,neqm
        do ne=1,neq
          READ(lrb,'(3i8)') nn,idof(ne),mdof(ne)
        enddo
        REWIND(lrb)
      endif
c
c     do nn=1,neq
c       write(*,'(4i5)') nn,ncp(nn),mdof(nn),idof(nn)
c     enddo
c
c ***** Compute the Hight of Stiffness Matrix "SK" *********************
c   ===== PARDISO Solver =====
      if( (isolvr.eq.1).or.(isolvr.eq.2) ) then
        CALL prepar(    nx,  nelx,   ndf,  node,
     &                 neq,  neqm,  NGSK, NGSKm,  nmax,
     &               nodsk,
     &                 ijk,  mdof,  idof, jdiag,
     &              iw_neq,  list,
     &              ierror )
        if(ierror.ne.0) RETURN
c
      elseif(isolvr.eq.0) then
c   ===== Conventional SKYLINE Solver =====
        CALL presky(    nx,  nelx,   ndf,  node,
     &                 neq,  neqm,  NGSK, NGSKm,
     &                 ijk,  mdof,  idof, jdiag,
     &              iw_neq,
     &              ierror )
        nodsk = 1
        if(ierror.ne.0) RETURN
      endif
c     WRITE(*,*) NGSKm, NGSK
c     stop
c
c 998 FORMAT(/,3x,'***** Reduced SKYLINE Hight Successfully!! *****',/,
c    &       8x,'Original Size of "SK"             :',i10,/,
c    &       8x,' Size of "SK" to be Solved        :',i10,/,
c    &       8x,'Size of "SK" Reordered by RCM     :',i10,/,
c    &       8x,' Size of Reduced "SK" to be Solved:',i10,/)
c
c     do ne=1,neq
c       write(*,*) ne,jdiag(ne)
c     enddo
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
