      subroutine mapcrs(   nel,  nelx,   mdf,  node,
     &                    NGSK,   neq,  neqm, NGSKo,  neqo,
     &                     ijk,  mdof, jdiag,index0,index1,
     &                  jcolmn,
     &                      sk,   ske,
     &                  ierror )
c
c ----------------------------------------------------------------------
c Purpose: Make the global stiffness matrix "SK" in CRS-format
c            from the element stiffness "ske"
c          CRS: Compressed Row-major Storage for PARDISO solver
c ----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension mdof(neq)
      dimension jdiag(neq+1)
      dimension index0(neqo+1)
      dimension index1(NGSKo)
      dimension jcolmn(nsolvr+1)
c
      dimension sk(NGSK)
c
      dimension ske(32,32)
c
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
c **********************************************************************
      NGSKm = NGSK -NGSKo
c
c **********************************************************************
c ***** Map the element Stiffness Matrix (ske) to the Global ***********
c          ( Non-Constrained DOF )
c **********************************************************************
      do 200 ia=1,node
        ijkia = ijk(ia,nel)
        if(ijkia.eq.0) GOTO 200
c
        do ii=1,mdf
          iai = mdf*(ia-1) +ii
          ijkiai = mdf*(ijkia-1) +ii
          m_ia = mdof(ijkiai)
c
c         ===== for Non-Constrained DOF (DISP) =====
          if(m_ia.le.neqm) then
c
            inds = jdiag(m_ia)
            inde = jdiag(m_ia +1) -1
c
            do 220 jb=1,node
              ijkjb = ijk(jb,nel)
              if(ijkjb.eq.0) GOTO 220
              do jj=1,mdf
                jbj = mdf*(jb-1) +jj
                ijkjbj = mdf*(ijkjb-1) +jj
                m_jb = mdof(ijkjbj)
c
c               --- for Non-constrained DOF's (disp)
                if(m_jb.le.neqm) then
                  do ind=inds,inde
                    mm = jcolmn(ind)
                    if(mm.eq.m_jb) then
c                     write(*,*) mm,ind,iai,jbj,ske(iai,jbj),ind
                      sk(ind) = sk(ind) +ske(iai,jbj)
                    endif
                  enddo
                endif
c
              enddo
c
  220       CONTINUE
c
c         ===== for Constrained DOF (DISP) =====
          else
c           write(*,'(/,10i5)') nel,ia,ii,iai,ijkia,ijkiai,m_ia
            mo_ia = m_ia -neqm
            inds = index0(mo_ia)
            inde = index0(mo_ia +1) -1
c           write(*,'(30x,4i5)') m_ia,inds,inde
c
            do 225 jb=1,node
              ijkjb = ijk(jb,nel)
              if(ijkjb.eq.0) GOTO 225
              do jj=1,mdf
                jbj = mdf*(jb-1) +jj
                ijkjbj = mdf*(ijkjb-1) +jj
                m_jb = mdof(ijkjbj)
c
c             --- Only for Lower Parts ---
                if(m_jb.le.m_ia) then
                  do ind=inds,inde
                    mm = index1(ind)
                    if(mm.eq.m_jb) then
                      nd = NGSKm +ind
c                   write(*,*) mm,ind,iai,jbj,ske(iai,jbj),nd
                      sk(nd) = sk(nd) +ske(iai,jbj)
                    endif
                  enddo
                endif
c
              enddo
  225       CONTINUE
c
          endif
        enddo
  200 continue
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
