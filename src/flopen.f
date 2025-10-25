      subroutine flopen(ierror,  iRCM)
c
      implicit double precision(a-h,o-z)
c
      character*30 infem,flname
      character*1 fl
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
c
      WRITE(*,100)
c
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
      WRITE(lwa,'(a7)') '/TITLE/'
      WRITE(lwa,*) 'PLSTss ver. 2.2'
      WRITE(lwa,'(a7)') '/LASTD/'
c
c     === Output File for Stress-Strain Curve ===
      flname = 'output/STS_'//infem(1:ni)//'.txt'
      OPEN(lwb,file=flname,status='unknown')
c
c     === Output File for Nodal Displacements & Forces ===
      flname = 'output/DIS_'//infem(1:ni)//'.txt'
      OPEN(lwc,file=flname,status='unknown')
c
c     === Output File for Some Norms ===
      flname = 'output/NOR_'//infem(1:ni)//'.txt'
      OPEN(lwd,file=flname,status='unknown')
c
c     === Output File for Plastic & Elastic Energy ===
      flname = 'output/ENE_'//infem(1:ni)//'.txt'
      OPEN(lwe,file=flname,status='unknown')
c
c     === Output File for Temperature ===
      flname = 'output/TMP_'//infem(1:ni)//'.txt'
      OPEN(lwf,file=flname,status='unknown')
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
