      subroutine writeboth(fmt, val1, val2, val3, val4, val5)
c
c     Write to both terminal (unit 6) and log file (unit 99)
c     Supports up to 5 values
c
      implicit none
      character*(*) fmt
      double precision val1, val2, val3, val4, val5
      logical file_opened
      integer nargs
c
c     Check if log file is opened
      inquire(unit=99, opened=file_opened)
c
c     Count number of arguments based on format string
      nargs = 1
      if (index(fmt, 'E14.6,A,E14.6') .gt. 0) nargs = 2
c
c     Write to terminal
      if (nargs .eq. 1) then
        write(*,fmt) val1
      else if (nargs .eq. 2) then
        write(*,fmt) val1, val2
      endif
c
c     Write to log file if opened
      if (file_opened) then
        if (nargs .eq. 1) then
          write(99,fmt) val1
        else if (nargs .eq. 2) then
          write(99,fmt) val1, val2
        endif
      endif
c
      return
      end