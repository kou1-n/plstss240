      subroutine rowupr(    nx,  nelx,   mdf,  node,
     &                     neq,  neqm,  NGSK, NGSKm,
     &                     ijk,  mdof,  idof, jdiag,jcolmn,
     &                  iw_neq,  list,ncolmn,index0,index1,
     &                  ierror )
c
c ----------------------------------------------------------------------
c Purpose: Define the Index Arrays for the PARDISO solver
c          values:   A real or complex array that contain the non-zeo
c          (sk)      entries of "A". The non-zero values of "A" are
c                    mapped into the "values" array using the ROW MAJOR,
c                    UPPER TRIANGULAR storage mapping.
c          column:   Element "i" of the integer array "column" contain
c          (jcolmn)  the number of the column in "A" that contained the
c                    value in "values".
c          rowIndex: Element "j" of the integer array "rowIndex" gives
c          (jdiag)   the index into the "values" array that contains the
c                    first non-zero element in a row "j" of "A".
c ----------------------------------------------------------------------
c Non-Zero Elements are Stored in ROW Major UPPER Triangular Format
c ----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension mdof(neq),idof(neq)
      dimension jdiag(neq+1)
      dimension jcolmn(nsolvr+1)
      dimension ncolmn(nodsk)
c
      dimension index0(neq-neqm +1)
      dimension index1(NGSK-NGSKm)
c
c --- working array ---
      dimension iw_neq(neq)
      dimension list(nmax,nx)
c
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
c **********************************************************************
c
c ***** Initilization **************************************************
      jcc = 0
c
      jdiag = 0 ! jdiag(:) = 0
      jdiag(1) = 1
c
      ncolmn = 0 ! ncolmn(:) = 0
c
c   ----- Index Arrays for Constrained DOF -----
      index0 = 0 ! index0(:) = 0
      index1 = 0 ! index1(:) = 0
c
c ***** Check the Dimension of the Arrasy for the List *****************
c         ( Already checked in the subroutine "prepar" )
c
c ***** Make List of Nodes which Elements are Associated to the Node ***
c              iw_neq(nx):   Number of the Elements
c              list(nmax,nx) ID of the Elements
c     ----- for usual elements
      do nel=1,nelx
        do no=1,node
          ijkia = ijk(no,nel)
          iw_neq(ijkia) = iw_neq(ijkia) +1
          list(iw_neq(ijkia),ijkia) = nel
        enddo
      enddo
c
c     write(*,*) ' '
c     do nn=1,nx
c       write(*,*) nn,iw_neq(nn),(list(kk,nn),kk=1,iw_neq(nn))
c     enddo
c
c ***** Make Index for Constrained DOF (index0, index1) ****************
      index0(1) = 1
      do ne=neqm+1,neq		! LOOP for Const. DOF's in modified order
        nd = idof(ne)		! Original DOF
        nn = (nd -1)/mdf +1	! Node ID of the DOF
        ns = index0(ne -neqm)	! Start Pointer for the DOF
        jcc = 0
c
        do jd=1,iw_neq(nn)	! LOOP for the Elements
          nc = list(jd,nn)	! ID of the Element, asso. with the NODE
c
c       --- for usual elements
          do no=1,node		! LOOP for the nodes
            ijkia = ijk(no,nc)	! ID of the node
            do mdl=1,mdf	! LOOP for dof
              me = mdf*(ijkia-1) +mdl	! Original DOF
              md = mdof(me)		! modified DOF
c
c           ----- Only for Lower Triangle -----
              if(md.le.ne) then
c               write(*,*) nc,no,ijkia,me,md,ns
c
c             ----- Check the Douplicity the Node ID -----
                ichk =0
                do jc=1,jcc
                  if(index1(ns-1 +jc).eq.md) ichk =1
                enddo
c             ----- store as the associated dof -----
                if(ichk.eq.0) then
                  jcc = jcc +1
                  index1(ns-1 +jcc) = md
                  do jci=jcc,2,-1
                    ir = index1(ns-1 +jci)
                    il = index1(ns-1 +jci-1)
                    if(ir.lt.il) then
                      index1(ns-1 +jci) = il
                      index1(ns-1 +jci-1) = ir
                    endif
                  enddo
                endif
              endif
c
            enddo
          enddo
        enddo
c
c     ----- Update the Diagonal Index -----
        index0(ne-neqm +1) = ns +jcc
c
      enddo
c
c ***** Find the Nodes which are Associated with the Target Node *******
      do nn=1,nx
        jcc = 0
        jd = jdiag(nn)
        ncolmn(jd) = nn
c
        do iw=1,iw_neq(nn)
          nel = list(iw,nn)
c
c         *** for Usual Elements ***
          do no=1,node
            ijkia = ijk(no,nel)
c
c         ----- Only for Lower Triangle -----
            if(ijkia.gt.nn) then
              ichk = 0
c
c           --- Check the Douplicity the Node ID ---
              do jc=1,jcc
                if(ncolmn(jd+jc).eq.ijkia) ichk = 1
              enddo
c
              if(ichk.eq.0) then
c
c             ----- Sort the ID in Increasing Order -----
                do jj=1,jcc
                  if(ncolmn(jd+jj).gt.ijkia) then
                    do kk=jcc,jj,-1
                      ncolmn(jd+kk+1) = ncolmn(jd+kk)
                    enddo
                    ncolmn(jd+jj) = ijkia
                    jcc = jcc +1
                    GOTO 111
                  endif
                enddo
c
c             ----- In the case the ID is maximum -----
                jcc = jcc +1
                ncolmn(jd+jcc) = ijkia
c
  111           CONTINUE
              endif
c
            endif
          enddo
        enddo
c
        jdiag(nn+1) = jd +jcc +1
c
      enddo
c
c ===== Copy "jdiag" to "iw_neq" =====
      do nn=1,nx+1
        iw_neq(nn) = jdiag(nn)
      enddo
c
      jdiag = 0 ! jdiag(:) = 0
c
c     write(*,*) jcc
c     write(*,*) (ncolmn(jj),jj=1,jcc)
c
c     write(*,*) ' '
c     write(*,*) 'iw_neq (jdiag)'
c     write(*,*) (iw_neq(nn),nn=1,nx+1)
c
c     write(*,*) ' '
c     write(*,*) 'ncolmn'
c     do nn=1,nx
c       write(*,'(10i3)') nn,(ncolmn(mm),mm=iw_neq(nn),iw_neq(nn+1)-1)
c     enddo
c
c ***** Expand the Nodal Index to the Index for DOFs *******************
c         ( ONLY for NON-Constrained DOF )
c     write(*,*) neq,neqm
      jdiag(1) = 1
c
c   ----- for Non-Constrained DOF -----
      do ne=1,neqm
        nd = idof(ne)		! Original ID of DOF
        nn = (nd -1)/mdf +1	! ID of the node
        np = iw_neq(nn)		! Pointer for the Diagonal of the node
        ns = jdiag(ne)		! Pointer for Diagonal of the dof
c
        jcc = 0
        jcolmn(ns+jcc) = ne
        do jd=iw_neq(nn),iw_neq(nn+1)-1
          nc = ncolmn(jd)	! ID of the associated node
c       ----- dispacement DOF's
          do mdl=1,mdf
            me = mdf*(nc-1) +mdl	! Original DOF
            md = mdof(me)		! modified DOF
            if((md.gt.ne).and.(md.le.neqm)) then
              jcc = jcc +1
              jcolmn(ns+jcc) = md
            endif
          enddo
        enddo
c
        jdiag(ne+1) = ns +jcc +1
c
      enddo
c
c **********************************************************************
 9000 FORMAT(/,/,
     &       3x,'**************************************************',/,
     &       3x,' Increase the Parameter "nmax" in sub. "chkary"!!',/,
     &       3x,'   Current value:',i5, ', but desired value:',i5,/,
     &       3x,'**************************************************')
c **********************************************************************
c     write(*,*) ' '
c     do ns=1,nsolvr
c       write(*,*) ns,jcolmn(ns)
c     enddo
c
c     write(*,*) ' '
c     do ne=1,neq+1
c       write(*,'(3i5)') ne,jdiag(ne)
c     enddo
c
c **********************************************************************
c **********************************************************************
      END
