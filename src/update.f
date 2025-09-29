      subroutine update(   neq,
     &                    idof,
     &                     del,    du,    u0,    u1,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension idof(neq)
c
      dimension del(neq),du(neq),u0(neq),u1(neq)
c
c **********************************************************************
c     del: modification of displacements during this N-R iteration
c     du:  displacement increments during this load step
c     u0:  total displacement at previous equilibrium state
c     u1:  total displacement at current load step
c **********************************************************************
c
c     do ne=1,neq
c       write(*,'(i5,3e15.5)') ne,del(ne)
c     enddo
c
c ***** Re-order the Computed Displacements and Potentials *************
c         ( Here, 'u1' is used for working array )
      do ne=1,neq
        id = idof(ne)
        u1(id) = del(ne)
      enddo
      del = u1 ! del(:) = u1(:)
c
c ***** Update Displacement Increment & Total Displacement *************
c     write(*,*) ' '
      du = du +del ! du(:) = du(:) +del(:)
      u1 = u0 +du ! u1(:) = u0(:) +du(:)
c
c **********************************************************************
c     write(*,*) ' '
c     do ne=1,neq
c       write(*,'(i5,4e15.5)') ne, u1(ne),u0(ne),du(ne),del(ne)
c       write(*,'(i5, e15.5)') ne, u1(ne)
c       write(*,'(a5, e15.5)') 'u0',      u0(ne)
c       write(*,'(a5, e15.5)') 'du',             du(ne)
c       write(*,'(a5, e15.5)') 'del',                    del(ne)
c     write(*,*) ' '
c     enddo
c     stop 'update'
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
