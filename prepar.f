      subroutine prepar(    nx,  nelx,   mdf,  node,
     &                     neq,  neqm,  NGSK, NGSKm,  nmax,
     &                   nodsk,
     &                     ijk,  mdof,  idof, jdiag,
     &                  iw_neq,  list,
     &                  ierror )
c
c ----------------------------------------------------------------------
c Purpose: Conut the Number of Non-Zero Element in "SK"
c          To Prepare the Arrays for the PARDISO software
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
c
      implicit double precision (a-h,o-z)
c
      dimension mdof(neq),idof(neq)
      dimension ijk(node,nelx)
c
c --- working array ---
      dimension iw_neq(neq)
      dimension jdiag(neq+1)
      dimension list(nmax,nx)
c
c **********************************************************************
c
c ***** Initilization **************************************************
      NGSK  = 0
      NGSKm = 0
c
      iw_neq = 0 ! iw_neq(:) = 0
c
c     nmax=9  !hsugiyama
c ***** Check the Dimension of the Arrasy for the List *****************
c   ( count the number of element associated to the node )
c     --- for usual elements
      do nel=1,nelx
        do no=1,node
          ijkia = ijk(no,nel)
          iw_neq(ijkia) = iw_neq(ijkia) +1
        enddo
      enddo
c
c     write(*,*) iw_neq
      nnmax =0
      do nn=1,nx
        if(iw_neq(nn).gt.nnmax) nnmax = iw_neq(nn)
      enddo
c
      if(nnmax.gt.nmax) then
        ierror = 10
        WRITE(*,9000) nmax,nnmax
        RETURN
      endif
c     write(*,*) nmax,nnmax
c
c     ( reset the working array )
      iw_neq = 0 ! iw_neq(:) = 0
c
c ***** Make List of Nodes which Elements are Associated to the Node ***
c     --- for usual elements
      do nel=1,nelx
        do no=1,node
          ijkia = ijk(no,nel)
          iw_neq(ijkia) = iw_neq(ijkia) +1
          list(iw_neq(ijkia),ijkia) = nel
        enddo
      enddo
c
c     do nn=1,nx
c       write(*,*) nn,iw_neq(nn),(list(kk,nn),kk=1,iw_neq(nn))
c     enddo
c
c
c ***** Count the Number of Nodes Associated with a Node ***************
c   --- from DISP to DISP with usual elements
      do nn=1,nx
        do mm=nn,nx
          jdiag(mm) = 0
        enddo
c
        do iw=1,iw_neq(nn)
          nel = list(iw,nn)
          do no=1,node
            ijkia = ijk(no,nel)
            if(ijkia.ge.nn) jdiag(ijkia) = 1
          enddo
        enddo
c       write(*,'(i5,3x,1000i1)') nn,(jdiag(mm),mm=1,nx)
c
        do mm=nn,nx
          NGSK = NGSK +jdiag(mm)
        enddo
c
      enddo
c
c   ===== Count the Non-Zero Element for Full Matrix =====
      nodsk = NGSK
c     write(*,*) nodsk
c     write(*,*) NGSK
      NGSK = (NGSK -nx)*2 +nx
c     write(*,*) NGSK
      NGSK = NGSK*mdf*mdf
c     write(*,*) NGSK
c
c   ===== Count the Non-Zero Element in Upper Triangle =====
c         ( for ALL DOF )
      NGSK = (NGSK -neq) /2 +neq
c     write(*,*) NGSK
c
c ***** Count the Number of Nodes Associated with a Node: "idof(ne)" ***
c       ( for Constrained DOF, Lower Triangle )
      do ne=neqm+1,neq
        nd = idof(ne)
        nn = (nd -1)/mdf +1
c
        do mm=1,ne
          jdiag(mm) = 0
        enddo
c
c       write(*,*) ne,nd,nn,iw_neq(nn)
        do iw=1,iw_neq(nn)
          nel = list(iw,nn)
          do no=1,node
            ijkia = ijk(no,nel)
            do md=1,mdf
              kdf = mdf*(ijkia -1) +md
              mm = mdof(kdf)
c             write(*,*) nel,no,ijkia,md,kdf,mm
              if(mm.le.ne) jdiag(mm) = 1
            enddo
          enddo
c
        enddo
c       write(*,'(i5,3x,1000i1)') nn,(jdiag(mm),mm=1,ne)
c
        do mm=1,ne
          NGSKm = NGSKm +jdiag(mm)
        enddo
c
      enddo
c
c **********************************************************************
c **********************************************************************
c
c ***** Compute the Toltal Number of Non-Zero Elements *****************
c   --- Non-Zero Elements for CONSTRAINED DOF ---
      NGSKm = NGSK -NGSKm
c
c   --- Non-Zero Elements for NON-CONSTRAINED DOF ---
      NGSKo = NGSK -NGSKm
c
c **********************************************************************
 9000 FORMAT(/,/,
     &       3x,'**************************************************',/,
     &       3x,' Increase the Parameter "nmax" in sub. "chkary"!!',/,
     &       3x,'   Current value:',i5, ', but desired value:',i5,/,
     &       3x,'**************************************************')
c **********************************************************************
c     write(*,*) NGSKm, NGSK
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
