      subroutine presky(    nx,  nelx,   mdf,  node,
     &                     neq,  neqm,  NGSK, NGSKm,
     &                     ijk,  mdof,  idof, jdiag,
     &                  iw_neq,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension mdof(neq),idof(neq),iw_neq(neq)
      dimension jdiag(neq+1)
c **********************************************************************
c
c ***** Compute the Skyline Hight for "SK" to be solved ****************
c     --- from DISP to DISP with usual elements
      do 100 nel=1,nelx
c
        do 120 ia=1,node
          ijkia = ijk(ia,nel)
          if(ijkia.eq.0) goto 120
          do ii=1,mdf
            niai = mdf*(ijkia -1) +ii
            m_ia = mdof(niai)
c
            if(m_ia.le.neqm) then
c             write(*,'(6i5)') nel,ia,ii, niai, m_ia,n_ia
c
              do 140 jb=1,node
                ijkjb = ijk(jb,nel)
                if(ijkjb.eq.0) goto 140
                do jj=1,mdf
                  njbj = mdf*(ijkjb -1) +jj
                  m_jb = mdof(njbj)
c
                  kij = m_ia -m_jb +1
                  if(kij.ge.1) then
                    jdiag(m_ia) = max0(jdiag(m_ia),kij)
                  endif
                enddo
  140         continue
            endif
c
          enddo
  120   continue
c
  100 continue
c
c     do ne=1,neqm
c       write(*,'(3i8)') ne,idof(ne),jdiag(ne)
c     enddo
c
c ***** Define the Diagonal Pointer for "SK" to be Solved **************
c     do ne=2,neq
      do ne=2,neqm
        jdiag(ne) = jdiag(ne-1) +jdiag(ne)
      enddo
c
c     do ne=1,neqm
c       write(*,'(3i8)') ne,idof(ne),jdiag(ne)
c     enddo
c
c     NGSK = jdiag(neq)
      if(neqm.eq.0) then
        NGSKm = 0
      else
        NGSKm = jdiag(neqm)
      endif
c
c ***** Count the # of Nodes Associated with a Constrainted Node *******
      do 200 ne=neqm+1,neq
        nd = idof(ne)
        nn = (nd -1)/mdf +1
c
c ===== Initilization of Array =====
        iw_neq = 0 ! iw_neq(:) = 0
c
c ===== Check the Eligibility of Each Element =====
        do 210 nel=1,nelx
          lap = 0
          do ia=1,node
            iai = ijk(ia,nel)
            if(iai.eq.nn) lap = 1
          enddo
c
          if(lap.eq.1) then
c           write(*,*) ne,nd,nn,nel
            do ia=1,node
              iai = ijk(ia,nel)
              do mm=1,mdf
                nd2 = mdf*(iai -1) +mm
                nd2m = mdof(nd2)
c               write(*,*) '        ',iai,nd2,nd2m
                iw_neq(nd2m) = 1
              enddo
            enddo
          endif
c
  210   continue
c       write(*,'(1000i1)') (iw_neq(mm),mm=1,neq)
c
c ===== Count the # of Nodes for Lower Part of "SK" =====
        ksum = 0
        do mm=1,ne
          if(iw_neq(mm).eq.1) then
            ksum = ksum +1
          endif
        enddo
c
        if(ne.eq.1) then
          jdiag(1) = 1
        else
          jdiag(ne) = jdiag(ne-1) +ksum
        endif
  200 continue
c
      NGSK = jdiag(neq)
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
