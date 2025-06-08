c     *********************************
c     *                               *
c     *  Create /SOLUT/ for CML file  *
c     *                               *
c     *********************************
c
      implicit double precision(a-h,o-z)
c
      9502 FORMAT(5i5)
      9505 FORMAT(2i5,f12.5)
c
      write(*,*) 'input the number of step on first loading'
      read(*,*) dff
      write(*,*) 'Input the load increments'
      read(*,*) df
      write(*,*)  'Input times of loading'
      read(*,*) nn
c
        write(lra,9502) nstep,mdum,mdum,mdum,mdum
        write(lra,9502) mdum,nn,mdum,mdum,mdum
c
      aaa1=0.d0
      aaa2=0.d0
c
      do i=1,nn
        if(i==1) then
          aaa1 = 1.d0
          aaa2 = dff
        else
          aaa1 = aaa1+(2.d0*dff)
          aaa2 = aaa2+(2.d0*dff)
        end if
c
        write(lra,9505) aaa1,aaa2,df
     enddo
c
c 