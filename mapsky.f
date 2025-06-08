      subroutine mapsky(   nel, ijdim,   mdf,  node,
     &                    NGSK,   neq,  neqm, NGSKo,  neqo,
     &                     ijk,  mdof, jdiag,index0,index1,
     &                      sk,   ske,
     &                  ierror )
c
c ----------------------------------------------------------------------
c Purpose: Make the global stiffness matrix "SK" in SKYLINE-format
c            from the element stiffness "ske"
c ----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,ijdim)
      dimension mdof(neq)
      dimension jdiag(neq+1)
      dimension index0(neqo+1)
      dimension index1(NGSKo)
c
      dimension sk(NGSK)
c
      dimension icdof(32)
c
      dimension ske(32,32)
c
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
c **********************************************************************
      NGSKm = NGSK -NGSKo
      iconst = 0
c
c ***** Map the element Stiffness Matrix (ske) to the Global ***********
c          ( Non-Constrained DOF )
      do 200 ia=1,node
        ijkia = ijk(ia,nel)
        if(ijkia.eq.0) GOTO 200
c
c     ***** For displacements DOF's *****
        do ii=1,mdf
          iai = mdf*(ia-1) +ii
          ijkiai = mdf*(ijkia-1) +ii
          m_ia = mdof(ijkiai)
c
c       ===== for Non-Constrained DOF =====
          if(m_ia.le.neqm) then
c
            do 220 jb=1,node
              ijkjb = ijk(jb,nel)
              if(ijkjb.eq.0) GOTO 220
              do jj=1,mdf
                jbj = mdf*(jb-1) +jj
                ijkjbj = mdf*(ijkjb-1) +jj
                m_jb = mdof(ijkjbj)
c
c               --- for Non-constrained DOF's
                if(m_jb.le.neqm) then
                  kij = m_jb -m_ia
                  if(kij.ge.0) then
                    nd = jdiag(m_jb) -kij
                    sk(nd) = sk(nd) +ske(iai,jbj)
                  endif
                endif
              enddo
  220       continue
c
c       ===== for Constrained DOF: Store Only Non-Zero Terms =====
          else
            iconst = iconst +1
            icdof(iconst) = m_ia
          endif
c
        enddo
  200 continue
c
c ***** Map the element Stiffness Matrix (ske) to the Global ***********
c          ( Constrained DOF )
c     write(*,'(2i3,9i5)') nel,iconst,(icdof(ic),ic=1,iconst)
      do 300 ic=1,iconst
        id = icdof(ic)
        ido = id -neqm
c
        do 310 ia=1,node
          ijkia = ijk(ia,nel)
          if(ijkia.eq.0) GOTO 310
          do ii=1,mdf
            ijkiai = mdf*(ijkia-1) +ii
            m_ia = mdof(ijkiai)
            if(m_ia.eq.id) then
              iai = mdf*(ia-1) +ii
              GOTO 315
            endif
          enddo
  310   continue
c
  315   CONTINUE
c       write(*,*) ia,ii,iai
c
        do 320 jb=1,node
          ijkjb = ijk(jb,nel)
          if(ijkjb.eq.0) GOTO 320
          do jj=1,mdf
            jbj = mdf*(jb-1) +jj
            ijkjbj = mdf*(ijkjb-1) +jj
            m_jb = mdof(ijkjbj)
            if(m_jb.le.id) then
c             write(*,'(10i5)') nel,jb,jj,jbj,ijkjb,ijkjbj,m_jb
c
              inds = index0(ido)
              inde = index0(ido+1) -1
c             write(*,'(3i3)') inds,inde
              do ind=inds,inde
                m_ia = index1(ind)
                if(m_ia.eq.m_jb) then
                  nd = NGSKm +ind 
                  sk(nd) = sk(nd) +ske(jbj,iai)
                endif
              enddo
            endif
          enddo
  320   continue
c
  300 continue
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
