      subroutine write_log(message)
c
c     Subroutine to write output to both terminal and log file
c     Unit 99 is used for log file
c
      implicit none
      character*(*) message
      logical file_opened
c
c     Check if log file is already opened
      inquire(unit=99, opened=file_opened)
c
c     Write to terminal
      write(*,'(A)') message
c
c     Write to log file if opened
      if (file_opened) then
        write(99,'(A)') message
        call flush(99)
      endif
c
      return
      end
c
c======================================================================
c
      subroutine write_log_fmt(fmt_str, value)
c
c     Formatted write to both terminal and log file
c
      implicit none
      character*(*) fmt_str
      double precision value
      logical file_opened
c
c     Check if log file is already opened
      inquire(unit=99, opened=file_opened)
c
c     Write to terminal
      write(*,fmt_str) value
c
c     Write to log file if opened
      if (file_opened) then
        write(99,fmt_str) value
        call flush(99)
      endif
c
      return
      end
c
c======================================================================
c
      subroutine open_log_file(infem, ni)
c
c     Open log file for analysis output
c
      implicit none
      character*50 infem, flname
      integer ni
c
c     Close if already opened
      close(99, err=100)
c
  100 continue
c     Open new log file
      flname = 'output/LOG_'//infem(1:ni)//'.txt'
      open(99, file=flname, status='unknown')
      write(99,'(A)') 'PLSTss Analysis Log'
      write(99,'(A)') '==================='
      write(99,*)
c
      return
      end
c
c======================================================================
c
      subroutine close_log_file()
c
c     Close log file
c
      implicit none
      logical file_opened
c
      inquire(unit=99, opened=file_opened)
      if (file_opened) then
        write(99,*)
        write(99,'(A)') '==================='
        write(99,'(A)') 'End of Analysis Log'
        close(99)
      endif
c
      return
      end