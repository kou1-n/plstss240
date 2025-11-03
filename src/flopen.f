      subroutine flopen(ierror,  iRCM)
c
      implicit double precision(a-h,o-z)
      intrinsic date_and_time
c
      character*30 infem,flname
      character*1 fl
      character*10 date_str,time_str
      character*5 zone
      character*19 timestamp_line
      integer date_vals(8)
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
c
      call date_and_time(date_str,time_str,zone,date_vals)
      write(timestamp_line,
     &     '(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)')
     &     date_vals(1),date_vals(2),date_vals(3),
     &     date_vals(5),date_vals(6),date_vals(7)
      WRITE(*,90) timestamp_line
      WRITE(*,100)
c
   90 FORMAT(/,5x,'[Analysis Timestamp] ',A,/)
  100 FORMAT(/,
     &   5x,'************************************************',
     & /,5x,'**                                            **',
     & /,5x,'**     Welcome to "PLSTss"  Version 2.2       **',
     & /,5x,'**                          ( 2012.11.20 )    **',
     & /,5x,'**                                            **',
     & /,5x,'************************************************')
c
 1000 CONTINUE
c
      WRITE(*,200)
  200 FORMAT(/,/,'*) Enter the input file name without extension:  ',$)
c
      READ(*,*) flname
c
      do 10 i=1,29
        fl=flname(30-i:30-i+1)
        if(fl.ne.' ') then
          infem=flname(1:30-i)
          ni=30-i
          goto 20
        endif
   10 continue
c
   20 CONTINUE
c
c ***** Open all the needed files *****
c
c     === Input file in CML-format ===
      flname = infem(1:ni)//'.cml'
      OPEN(lra,file=flname,status='old',err=1000)
c
c     === Output File in CML-format ===
      flname = 'output/RES_'//infem(1:ni)//'.cml'
      OPEN(lwa,file=flname,status='unknown')
      WRITE(lwa,'(A)') '# Analysis Timestamp: '//timestamp_line
      WRITE(lwa,'(a7)') '/TITLE/'
      WRITE(lwa,*) 'PLSTss ver. 2.2'
      WRITE(lwa,'(a7)') '/LASTD/'
c
c     === Output File for Stress-Strain Curve ===
      flname = 'output/STS_'//infem(1:ni)//'.txt'
      OPEN(lwb,file=flname,status='unknown')
      WRITE(lwb,'(A)') '# Analysis Timestamp: '//timestamp_line
c
c     === Output File for Nodal Displacements & Forces ===
      flname = 'output/DIS_'//infem(1:ni)//'.txt'
      OPEN(lwc,file=flname,status='unknown')
      WRITE(lwc,'(A)') '# Analysis Timestamp: '//timestamp_line
c
c     === Output File for Some Norms ===
      flname = 'output/NOR_'//infem(1:ni)//'.txt'
      OPEN(lwd,file=flname,status='unknown')
      WRITE(lwd,'(A)') '# Analysis Timestamp: '//timestamp_line
c
c     === Output File for Plastic & Elastic Energy ===
      flname = 'output/ENE_'//infem(1:ni)//'.txt'
      OPEN(lwe,file=flname,status='unknown')
      WRITE(lwe,'(A)') '# Analysis Timestamp: '//timestamp_line
c
c     === Output File for Temperature ===
      flname = 'output/TMP_'//infem(1:ni)//'.txt'
      OPEN(lwf,file=flname,status='unknown')
      WRITE(lwf,'(A)') '# Analysis Timestamp: '//timestamp_line
c
c     === Input File for Renumbered DOF ===
      flname = 'RNM_'//infem(1:ni)//'.txt'
      OPEN(lrb,file=flname,status='old',err=1100)
      WRITE(*,1150)
      iRCM = 1
c
      RETURN
c
c ***** In the Case "RNM_***.cml" file does not exist ******************
 1100 CONTINUE
      WRITE(*,1151)
      iRCM = 0
c
      RETURN
c
c **********************************************************************
 1150 FORMAT(/,/,3x,'***** Information for the Renumbered DOF',
     &         ' exist! *****',/,
     &         9x,'Use the Renumbered DOF for Assemble',/ )
c
 1151 FORMAT(/,/,3x,'***** Information for the Renumbered DOF',
     &         ' does NOT exist! *****',/,
     &         9x,'Use the Default Numbering for Assemble',/ )
c
c **********************************************************************
c **********************************************************************
      END
