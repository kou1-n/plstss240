      subroutine inform(    nx,  nelx,   mdf,  node,
     &                     neq,  neqm, NGSKo,  neqo,
     &                     ijk,  mdof,  idof,
     &                  iw_neq,index0,index1,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk(node,nelx)
      dimension mdof(neq),idof(neq),iw_neq(neq)
      dimension index0(neqo+1)
      dimension index1(NGSKo)
c
c **********************************************************************
      index0(1) = 1
c
c ***** Non-Zero Information of"SK" for Constrainted Nodes *************
      do 200 ne=neqm+1,neq    ! LOOP for Const. DOF's in modified order
        nd = idof(ne)         ! Original DOF
        nn = (nd -1)/mdf +1   ! Node ID of the org-DOF
c       write(*,*) ne,nd,nn
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
            do ia=1,node
              iai = ijk(ia,nel)
              do mm=1,mdf
                nd2 = mdf*(iai -1) +mm
                nd2m = mdof(nd2)
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
            id0 = index0(ne -neqm)
            index1(id0 +ksum) = mm
c           write(*,*) id0,ne,neqm,ksum,mm,index1(id0+ksum)
            ksum = ksum +1
          endif
        enddo
        index0(ne -neqm +1) = index0(ne -neqm) +ksum
  200 continue
c
c **********************************************************************
c     do ne=neqm+1,neq
c       write(*,*) index0(ne-neqm)
c     enddo
c     write(*,*) index0
c     write(*,*) index1
c
c     do ne=1,neqo
c       write(*,*) ' '
c       write(*,*) ne
c       write(*,'(15i4)') (index1(mm),mm=index0(ne),index0(ne+1)-1)
c     enddo
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
