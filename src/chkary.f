      subroutine chkary(     a,
     &                    NDIM,  NGSK, NGSKo,   neq,  neqo,
     &                    iRCM,
     &                      m1,    m2,    m3,    m4,    m5,
     &                      m6,    m7,    m8,    m9,   m10,
     &                     m11,   m12,   m13,   m14,   m15,
     &                     m16,   m17,   m18,   m19,   m20,
     &                     m21,   m22,   m23,   m24,   m25,
     &                     m26,   m27,   m28,   m29,
c
     &                      n1,    n2,    n3,    n4,    n5,
     &                      n6,    n7,    n8,    n9,   n10,
     &                     n11,   n12,   n13,   n14,   n15,
     &                     n16,   n17,   n18,   n19,   n20,
     &                     n21,   n22,   n23,   n24,   n25,
     &                     n26,   n27,   n28,   n29,   n30,
     &                     n31,   n32,   n33,   n34,   n35,
     &                     n36,   n37,   n38,   n39,   n40,
     &                     n41,   n42,   n43,   n44,   n45,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension a(NDIM)
c
      character*80 title
      character*7  charin
      character*7  key000,keyLAS,keyEND,keySOL,
     &             keyPST,keyPDI,keyPFO,keyPSN,keyPSS,keyPTE,
     &             key010,
     &             key100,key233,key236,key244,key248,key264,key268,
     &             key500,key550,key600,key610,key620,key650
      character*5  w_dim
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
      common /print/ lpstp,luprt,lfprt,lnprt,lsprt,lbprt
      common /tvalu/ ctol,stol
c
c **********************************************************************
c  isolvr: Define the Solver
c          = -1: Solver is NOT Defined
c          =  0: Conventional SKYLINE solver (Default)
c          =  1: PARDISO Solver
c          =  2: RCI Conjugate Gradient Solver (Ver. 1.1--)
c
      isolvr = -1
c
      nstep = -1
      lpstp = 0
c
      nmax = 10  !hsugiyama
c     nmax = 8
c
      REWIND(lra)
      WRITE(*,100)
c
c ***** Definition of CML-formatted Formats ****************************
  100 FORMAT(/,5x,'CHECKing the FEM data ',$)
  101 FORMAT('.',$)
  102 FORMAT('. done!!')
c
 9000 FORMAT(a7)
 9001 FORMAT(a80)
 9500 FORMAT(5i8)
 9501 FORMAT(i8,i5)
 9502 FORMAT(5i5)
 9503 FORMAT(i5,i8,i5)
 9504 FORMAT(f12.5)
 9505 FORMAT(2i5,f12.5)
 9506 FORMAT(i5,i8)
c
c ***** Initialization of Basic Parameters *****************************
      nx = 0
      nelx = 0
      ndf = 0
      node = 0
      nsn = 0
      ngaus = 0
      nel_t3 = 0
      nel_t6 = 0
      nel_q4 = 0
      nel_q8 = 0
      nel_t4 = 0
      nel_h8 = 0
      lnum = 1
      incomp = 0 !hsugiyama
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
c ***** Check Dimensions of the FE-model *******************************
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
        READ(lra,9500) nx
        do nn=1,nx
          READ(lra,*)
        enddo
c
c /TRIA3/ ***** 3-nodes TRIANGLE element (2D) *****
      elseif(charin.eq.key233) then
        WRITE(*,101)
        READ(lra,9501) nel_t3,ndum
        do nel=1,nel_t3
          READ(lra,*)
        enddo
        ndf   = max0(ndf,2)
        node  = max0(node,3)
        nsn   = max0(nsn,2)
        ngaus = max0(ngaus,1)
c
c /TRIA6/ ***** 6-nodes TRIANGLE element (2D) *****
      elseif(charin.eq.key236) then
        WRITE(*,101)
        READ(lra,9501) nel_t6,incomp
        do nel=1,nel_t6
          READ(lra,*)
        enddo
        ndf   = max0(ndf,2)
        node  = max0(node,6)
        nsn   = max0(nsn,2)
        if(incomp.eq.2) then
          ngaus = max0(ngaus,4)
        else
          ngaus = max0(ngaus,7)
        endif
c
c /QUAD4/ ***** 4-nodes Isoparametric element (2D) *****
      elseif(charin.eq.key244) then
        WRITE(*,101)
        READ(lra,9501) nel_q4,ndum
        do nel=1,nel_q4
          READ(lra,*)
        enddo
        ndf   = max0(ndf,2)
        node  = max0(node,4)
        nsn   = max0(nsn,2)
        ngaus = max0(ngaus,4)
c
c /QUAD8/ ***** 8-nodes Isoparametric element (2D) *****
      elseif(charin.eq.key248) then
        WRITE(*,101)
        READ(lra,9501) nel_q8,ndum
        do nel=1,nel_q8
          READ(lra,*)
        enddo
        ndf   = max0(ndf,2)
        node  = max0(node,8)
        nsn   = max0(nsn,2)
        ngaus = max0(ngaus,9)
c
c /TTRA4/ ***** 4-nodes Tetrahedral element (3D) *****
      elseif(charin.eq.key264) then
        WRITE(*,101)
        READ(lra,9501) nel_t4,ndum
        do nel=1,nel_t4
          READ(lra,*)
        enddo
        ndf   = max0(ndf,3)
        node  = max0(node,4)
        nsn   = max0(nsn,3)
        ngaus = max0(ngaus,1)
c
c /HEXA8/ ***** 8-nodes Hexahedral element (3D) *****
      elseif(charin.eq.key268) then
        WRITE(*,101)
        READ(lra,9501) nel_h8,ndum
        do nel=1,nel_h8
          READ(lra,*)
        enddo
        ndf   = max0(ndf,3)
        node  = max0(node,8)
        nsn   = max0(nsn,4)
        ngaus = max0(ngaus,8)
c
c /MATER/ ***** Material Properties *****
      elseif(charin.eq.key500) then
        WRITE(*,101)
        READ(lra,9502) lmat
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
        READ(lra,9502) lang
        do la=1,lang
          READ(lra,*)
        enddo
c
c /CONST/ ***** 1st boundary condition: Constraint *****
      elseif(charin.eq.key600) then
        WRITE(*,101)
        READ(lra,9502) locc,nspc,mpc
        do lo=1,locc
          READ(lra,*)
        enddo
        do ns=1,nspc
          READ(lra,*)
        enddo
        do mp=1,mpc
          READ(lra,*)
        enddo
c
c /PCNST/ ***** Constraint condition for the Potential values *****
c     elseif(charin.eq.key610) then
c       WRITE(*,101)
c       READ(lra,9502) ntn,ndum,ndum
c       do nt=1,ntn
c         READ(lra,*)
c       enddo
c
c /LOADC/ ***** 2nd boundary condition: Loading condition *****
      elseif(charin.eq.key650) then
        WRITE(*,101)
        READ(lra,9502) lnum
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
        READ(lra,9502) isolvr,msol1,msol2,idum,idum
        READ(lra,9502) idum,idum,idum,idum,mtol
c   === TOLERANCE for RCI CG Solver ===
        if(isolvr.eq.2) then
          stol = 10.d0**(-1*mtol)
          if(stol.gt.1.d-5) then
            WRITE(*,9995)
            ierror = 32
            RETURN
          endif
        endif
c     msol1 (for PARDISO Solver)
c           = 0: NO output from PARDISO Solver
c           = 1: Prints Statistical Information 
c     msol1 (for RCI CG Solver)
c           = 0: DONOT Use "Preconditioner"
c           = 1: USE "Jacobi Preconditioner"
c
c     --- flag for solution coltrol ---
      jsol = 0
c
c /SOLUT/ ***** Solution Procedures *****
      elseif(charin.eq.keySOL) then
        WRITE(*,101)
        READ(lra,9502) nstep,mdum,mdum,mdum,iarc
        READ(lra,9502) mdum,itrobj,mdum,mdum,mtol
c    === TOLERANCE of this computation ===
        ctol = 10.d0**(-1*mtol)
        if(ctol.gt.1.d-5) then
          WRITE(*,9995)
          ierror = 31
          RETURN
        endif
c
c    === for Load/Displacement Control ===
        if(iarc.eq.0) then
c       --- Load increments are defined EXPLICITLY
          if(itrobj.ge.1) then
            do it=1,itrobj
              READ(lra,9505) mdum,nstped,dum
            enddo
            if(nstped.ne.nstep) then
              WRITE(*,9996)
              ierror = 19
              RETURN
            endif
          endif
c    === for Arc-Length Method ===
        else
          READ(lra,9504) arcini
        endif
c
c /PSTEP/ ***** Print control section *****
      elseif(charin.eq.keyPST) then
        WRITE(*,101)
        READ(lra,9502) lpstp
        if(lpstp.eq.0) then
          WRITE(*,640)
          lpstp = -1
        endif
c
        if(lpstp.lt.0) then
          if(nstep.le.0) then
            WRITE(*,9997)
            ierror = 12
            RETURN
          endif
          nsect = -lpstp
          lpstp = nstep/nsect
        elseif(lpstp.gt.0) then
          do nn=1,lpstp
            READ(lra,*)
          enddo
        endif
c
c /PDISP/ ***** ID & DOF where nodal displacements are needed *****
      elseif(charin.eq.keyPDI) then
        WRITE(*,101)
        READ(lra,9502) luprt
        do nn=1,luprt
          READ(lra,*)
        enddo
c
c /PFOCE/ ***** ID & DOF where nodal forces are needed *****
      elseif(charin.eq.keyPFO) then
        WRITE(*,101)
        READ(lra,9502) lfprt
        do nn=1,lfprt
          READ(lra,*)
        enddo
c
c /PSTRN/ ***** ID & DOF where elemental STRAIN are needed *****
      elseif(charin.eq.keyPSN) then
        WRITE(*,101)
        READ(lra,9502) lnprt
        do nn=1,lnprt
          READ(lra,*)
        enddo
c
c /PSTRS/ ***** ID & DOF where elemental STRESS are needed *****
      elseif(charin.eq.keyPSS) then
        WRITE(*,101)
        READ(lra,9502) lsprt
        do nn=1,lsprt
          READ(lra,*)
        enddo
c
c /PTEMP/ ***** ID where elemental Temperature are needed *****
       elseif(charin.eq.keyPTE) then
       WRITE(*,101)
       read(lra,9506) lbprt
       do nn=1,lbprt
         read(lra,*)
       enddo
c
c /PENER/ ***** ID where elemental Energy are needed *****
c       elseif(charin.eq.keyPEN) then
c         write(*,101)
c         read(lra,9502) laprt
c         do nn=1,laprt
c           read(lra,*)
c         enddo
c
c /LASTD/ ***** Last of input data *****
      elseif(charin.eq.keyLAS) then
c
c /ENDOF/ ***** End of this file *****
      elseif(charin.eq.keyEND) then
c       WRITE(*,102)
        GOTO 2000
c
      else
c       write(*,*) 'charin'
        WRITE(*,*) charin
      endif
c
      GOTO 1000
c
c ***** Error Section **************************************************
 9994 FORMAT(/,/,/,'**********',
     &       /,'The Telerance defined in /SOLVR/ seems too LARGE! ',
     &       /,'     Check the parameters in /SOLVR/ section',/)
c
 9995 FORMAT(/,/,/,'**********',
     &       /,'The Telerance defined in /SOLUT/ seems too LARGE! ',
     &       /,'     Check the parameters in /SOLUT/ section',/)
c
 9996 FORMAT(/,/,/,'**********',
     &       /,'Number of Loading Steps (nstep) is NOT consistent ',
     &       /,'     Check the parameters in /SOLUT/ section',/)
c
 9997 FORMAT(/,'/PSTEP/ statement must be defined after ',
     &         '/SOLUT/ statement in CML-format')
c
 9998 CONTINUE
      WRITE(*,*) 'Something WRONG in your CML-format!!'
      ierror = 1
      RETURN
c
 9999 CONTINUE
      WRITE(*,*) 'You forgot /ENDOF/ indicator in your CML-file'
      ierror = 2
      RETURN
c
c ***** Compute Total Number of Finite Elements ************************
 2000 CONTINUE
      nelx = nel_t3 +nel_t6 +nel_q4 +nel_q8
     &      +nel_t4 +nel_h8
c
      if((nx.le.0).or.(nelx.le.0)) then
        WRITE(*,*) 'You forget the most IMPORTANT informations!!'
        ierror = 3
        RETURN
      endif
c
c ***** Check Section **************************************************
c /SOLUT/: Number of Loading Steps
      if(nstep.eq.-1) then
        WRITE(*,9996)
        ierror = 19
        RETURN
      endif
c
c /PSTEP/: Number of Printing Steps
      if(lpstp.eq.0) then
        WRITE(*,640)
        lpstp = -1
c
        nsect = -lpstp
        lpstp = nstep/nsect
      endif
c
c ***** Define the Dimension of the Problem ****************************S
      mdf = ndf
c     neq = (ndf +1)*nx
      neq = mdf*nx
c
c ***** Estimate Size of Stiffness Matrix "SK" *************************
c   === Define the Pointers for Estimation ===
      nelxh = nelx/2 +mod(nelx,2)
      nxh   = nx/2 +mod(nx,2)
      neqh  = neq/2 +mod(neq,2)
c
      neq1h = (neq+1)/2 +mod(neq+1,2)
      iRCMh = iRCM/2 +mod(iRCM,2)
c
      k1 = 1
      k2 = k1 +node*nelxh
      k3 = k2 +neqh
      k4 = k3 +neqh
      k5 = k4 +neqh
      k6 = k5 +neqh
c
      k7 = k6 +neqh
      k8 = k7 +nmax*nxh
c
      kend = k8 -1
c
c   === Insufficient Size Parameter "NDIM" ===
      if(kend.gt.NDIM) then
        emem = 8.d0*dble(NDIM)
        w_dim = ' Byte'
        if(emem.gt.1024.d0) then
          emem = emem/1024.d0
          w_dim = 'KByte'
          if(emem.gt.1024.d0) then
            emem = emem/1024.d0
            w_dim = 'MByte'
          endif
        endif
c
        WRITE(*,601)   kend,  NDIM,  emem, w_dim
        ierror = 4
        RETURN
      endif
c
c   === In the Case, the Solver is NOT Defined ===
      if( (isolvr.ne.0).and.(isolvr.ne.1).and.
     &    (isolvr.ne.2) ) then
        WRITE(*,629)
        isolvr = 0
      endif
c
c   === Estimate the Size of "SK" ===
      CALL estima(   neq,  neqm,  NGSK, NGSKm,  iRCM,
     &             a(k1), a(k2), a(k3), a(k4), a(k5),
     &             a(k6), a(k7),
     &            ierror )
      if(ierror.ne.0) RETURN
c
c   === Compute the Additional Parameters ===
      NGSKo = NGSK -NGSKm
      neqo = neq -neqm
c
c     neqsol: dimension for working array in solvers
c     dw_sol(neqsol)
      if( isolvr.eq.1) then
        nsolvr = NGSKm
        neqsol = neq
      elseif(isolvr.eq.2) then
        nsolvr = NGSKm
        neqsol = 6*neq
      elseif(isolvr.eq.0) then
        nsolvr = 1
        neqsol = 1
      endif
c
c     write(*,*) nsolvr
c
c     write(*,*) NGSKo,NGSK,NGSKm
c     write(*,*) neqo,neq,neqm
c
c ****** Define Address Pointers for Each Array ************************
      CALL addres(  nend,  NGSK, NGSKo,   neq,  neqo,
     &                m1,    m2,    m3,    m4,    m5,
     &                m6,    m7,    m8,    m9,   m10,
     &               m11,   m12,   m13,   m14,   m15,
     &               m16,   m17,   m18,   m19,   m20,
     &               m21,   m22,   m23,   m24,   m25,
     &               m26,   m27,   m28,   m29,
c
     &                n1,    n2,    n3,    n4,    n5,
     &                n6,    n7,    n8,    n9,   n10,
     &               n11,   n12,   n13,   n14,   n15,
     &               n16,   n17,   n18,   n19,   n20,
     &               n21,   n22,   n23,   n24,   n25,
     &               n26,   n27,   n28,   n29,   n30,
     &               n31,   n32,   n33,   n34,   n35,
     &               n36,   n37,   n38,   n39,   n40,
     &               n41,   n42,   n43,   n44,   n45,
     &            ierror )
      if(ierror.ne.0) RETURN
c
c     write(*,'(10i5)') m1,m2,m3,m4,m5
c     write(*,'(10i5)') m6,m7,m8,m9,m10
c     write(*,'(10i5)') m11,m12,m13,m14,m15
c     write(*,'(10i5)') m16,m17,m18,m19
c,m20
c
c     write(*,'(10i5)') n1,n2,n3,n4,n5
c     write(*,'(10i5)') n6,n7,n8,n9,n10
c     write(*,'(10i5)') n11,n12,n13,n14,n15
c,n16,n17
c,n18,n19,n20   
c
c ***** Display the Name of Solver which are Choosen *******************
      if(isolvr.eq.0) then
c   ===== SKYLINE Solver =====
        WRITE(*,630)
c
      elseif(isolvr.eq.1) then
c   ===== PARDISO Solver =====
        WRITE(*,631)
c
      elseif(isolvr.eq.2) then
c   ===== RCI CG Solver =====
        if(msol1.eq.0) then
          WRITE(*,632) stol
        elseif( msol1.eq.1) then
          WRITE(*,633) stol
        else
          WRITE(*,634) stol
          msol1 = 0
        endif
      endif
c
c ***** Print-out Parameters *******************************************
      WRITE(*,610)     nx,  nelx,   neq,  neqm,  NGSK, NGSKm, ndf,
     &               node,  lmat,  lang,  lnum, npoin,
     &              npres, nbody,  nspc,
     &               ctol
c
      emem = 8.d0*dble(NDIM)
      w_dim = ' Byte'
      if(emem.gt.1024.d0) then
        emem = emem/1024.d0
        w_dim = 'KByte'
        if(emem.gt.1024.d0) then
          emem = emem/1024.d0
          w_dim = 'MByte'
        endif
      endif
c
      if(nend.gt.NDIM) then
        emem = 8.d0*dble(NDIM)
        w_dim = ' Byte'
        if(emem.gt.1024.d0) then
          emem = emem/1024.d0
          w_dim = 'KByte'
          if(emem.gt.1024.d0) then
            emem = emem/1024.d0
            w_dim = 'MByte'
          endif
        endif
c
c       WRITE(*,601) nend,NDIM
        WRITE(*,601)   nend,  NDIM,  emem, w_dim
        ierror = 7
        RETURN
      else
        WRITE(*,620)   nend,  NDIM,  emem, w_dim
      endif
c
c **********************************************************************
  641 FORMAT(/,3x,
     &       '***** Parameter "nstep" is NOT defined *****',
     &       /,10x,'Default Number is used "nstep = 1" ',/)
c
  640 FORMAT(/,3x,
     &       '***** Parameter "lpstp" is NOT defined *****',
     &       /,10x,'Default Number is used "lpstp = 1" ',/)
c
  634 FORMAT(/,3x,
     &       '***** RCI CG Solver will be Used for the Analysis *****',
     &       /,
     &       '       Default Settings (w/o Preconditioner) is applied',
     &       /,'       *) TOLERANCE for the CG Solver:',1pe15.5 )
c
  633 FORMAT(/,3x,
     &       '***** RCI CG Solver will be Used for the Analysis *****',
     &       /,
     &       '       with "Jacobi Preconditioner"',
     &       /,'       *) TOLERANCE for the PCG Solver:',1pe15.5 )
c
  632 FORMAT(/,3x,
     &       '***** RCI CG Solver will be Used for the Analysis *****',
     &       /,
     &       '       without "Preconditioner"',
     &       /,'       *) TOLERANCE for the CG Solver:',1pe15.5 )
c
  631 FORMAT(/,3x,
     &       '***** PARDISO Solver will be Used for the Analysis *****',
     &       /,
     &       '       with Compressed Row Storage format for Upper Part')
c
  630 FORMAT(/,3x,
     &       '***** SKYLINE Solver for Symmetric Matrix *****',/,
     &       '       with profile storage format for Lower Part')
c
  629 FORMAT(///,1x,'***** WARNNING!! Solver is NOT Defined *********',
     &         /,1x,'  As a Default Setting, ',
     &         /,1x,'    SKYLINE Solver will be Used in the Following',
     &         /,1x,'************************************************',
     &       /// )
c
  620 FORMAT(/,
     &       3x,'***** PLSTss Array Size for Computation *****',
     &     /,5x,'*) Total Array Size Required for Computation: ',i10,
     &     /,5x,'*) Array Size Declared in "main.f" (NDIM):    ',i10,
     &     /,5x,'*) Approximately Memory Size:                 ',
     &                    f10.3,a6,
     &     / )
c
  610 FORMAT(/,
     &       3x,'***** PLSTss Parameters for FE-model *****',
     &     /,5x,' 1) Total Number of Nodes:             ',i10,
     &     /,5x,' 2) Total Number of Elements:          ',i10,
     &     /,5x,' 3) Total Number of Degrees of Freedom:',i10,
     &     /,5x,' 4) Total DOF to be solved:            ',i10,
     &     /,5x,' 5) Total Stiffness Array Size:        ',i10,
     &     /,5x,' 6) Total SAS to be solved:            ',i10,
     &     /,5x,' 7) Number of DOF/node:                ',i10,
     &     /,5x,' 8) Number of Nodes/Element:           ',i10,
     &     /,5x,' 9) Number of Material Properties:     ',i10,
     &     /,5x,'10) Number of Material Angles:         ',i10,
     &     /,5x,'11) Number of Loading Sets:            ',i10,
     &     /,5x,'12) Number of Point Loads:             ',i10,
     &     /,5x,'13) Number of Distributed Loads:       ',i10,
     &     /,5x,'14) Number of Body Forces:             ',i10,
     &     /,5x,'15) Number of Prescribed Displacements:',i10,
     &     /,
     &     /,5x,' *) TOLERANCE of this computation:',1pe15.5,
     &     / )
c
  601 FORMAT(///,1x,'***** Out of Memory *****',
     &         /,1x,'Total array size required is  : ',i15,
     &         /,1x,'but NDIM is currently defined : ',i15,
     &         /,5x,'(Approximately Memory Size: ',f10.3,a6,')',
     &         /,1x,'Increase "NDIM" in parameter statement in main.f',
     &         /// )
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
