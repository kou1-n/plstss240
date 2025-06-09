      subroutine output(    in,   itr,   neq,
     &                   mpstp, muprt, mfprt, mnprt, msprt,
     &                   mbprt, dfact, unorm, tnorm,
     &                    disp, sigma, epsln,   von,  fint,
     &                     eps,  pene,  eene,tene_e,tene_p,
     &                    temp,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      character*7 keyLAS,keyEND,
     &            key800,key810,key820,key830
c
      dimension mpstp(lpstp)
      dimension muprt(2,luprt)
      dimension mfprt(2,lfprt)
      dimension mnprt(2,lnprt)
      dimension msprt(2,lsprt)
      dimension mbprt(2,lbprt)
c
      dimension angle(3,lang)
c
      dimension disp(neq),fint(neq)
      dimension sigma(6,nelx),epsln(6,nelx)
      dimension von(nelx),eps(nelx),
     &          pene(nelx),eene(nelx),temp(nelx)
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /print/ lpstp,luprt,lfprt,lnprt,lsprt,lbprt
c
      save tene_p_prev
      data tene_p_prev/0.d0/
c
c **********************************************************************
      keyLAS = '/LASTD/'
      keyEND = '/ENDOF/'
      key800 = '/NODAL/'
      key810 = '/ELMTL/'
      key820 = '/NOPZT/'
      key830 = '/ELPZT/'
c
      m1 = 1
      dum = 0.d0
c
 9000 FORMAT(a7)
 9001 FORMAT(3i5)
 9002 FORMAT(i8)
 9003 FORMAT(i8,6e13.5)
 9004 FORMAT(8x,6e13.5)
c
 9101 FORMAT(i8,4e16.8)
c
 9201 FORMAT('%Step',$)
 9202 FORMAT(i8,'-',i1,'(disp)',$)
 9203 FORMAT(i8,'-',i1,'(foce)',$)
 9204 FORMAT(i5,1x,$)
 9205 FORMAT(e16.8,$)
 9211 FORMAT('%Step',$)
 9212 FORMAT(i8,'-',i1,'(strn)',$)
 9213 FORMAT(i8,'-',i1,'(strs)',$)
 9214 FORMAT(i8,'-',i1,'(temp)',$)
 9215 FORMAT('% for 2D Problem... 1: xx, 2: yy, 3: xy')
 9216 FORMAT('% for 3D Problem... 1: xx, 2: yy, 3: zz ',
     &                           '4: yz, 5: zx, 6: xy')
c
c ***** Output for Stress-Strain Curve (STS_***.txt) *******************
c     sigma(1,nel) = sts(1,1)
c     sigma(2,nel) = sts(2,1)
c     sigma(3,nel) = sts(3,1)
c     sigma(4,nel) = sts(2,2)
c     sigma(5,nel) = sts(3,2)
c     sigma(6,nel) = sts(3,3)
c     WRITE(lwb,9101) in,epsln(1,1),sigma(1,1)
      if( (lnprt+lsprt).gt.1) then
c   === Hearders ===
        if(in.eq.1) then
          if(ndf.eq.2) then
            WRITE(lwb,9215)
          elseif(ndf.eq.3) then
            WRITE(lwb,9216)
          endif
          WRITE(lwb,9211)
          do ln=1,lnprt
            WRITE(lwb,9212) (mnprt(i,ln),i=1,2)
          enddo
          do ls=1,lsprt
            WRITE(lwb,9213) (msprt(i,ls),i=1,2)
          enddo
          WRITE(lwb,*) ' '
        endif
c
c   === Loading Step ===
        WRITE(lwb,9204) in
c
c   === Strain ===
        do ln=1,lnprt
          mm = mnprt(1,ln)
          nn = mnprt(2,ln)
          WRITE(lwb,9205) epsln(nn,mm)
        enddo
c
c   === Stress ===
        do ls=1,lsprt
          mm = msprt(1,ls)
          nn = msprt(2,ls)
          WRITE(lwb,9205) sigma(nn,mm)
        enddo
c
        WRITE(lwb,*) ' '
      endif
c
c ***** Output of Nodal Displacements & Forces (DIS_***.txt) ***********
      if( (luprt+lfprt).gt.1) then
c   === Hearders ===
        if(in.eq.1) then
          WRITE(lwc,9201)
          do lu=1,luprt
            WRITE(lwc,9202) (muprt(i,lu),i=1,2)
          enddo
          do lf=1,lfprt
            WRITE(lwc,9203) (mfprt(i,lf),i=1,2)
          enddo
          WRITE(lwc,*) ' '
        endif
c
c   === Loading Step ===
        WRITE(lwc,9204) in
c
c   === Displacement ===
        do lu=1,luprt
          mm = ndf*(muprt(1,lu)-1) +muprt(2,lu)
          WRITE(lwc,9205) disp(mm)
        enddo
c
c   === Forces ===
        do lf=1,lfprt
          mm = ndf*(mfprt(1,lf)-1) +mfprt(2,lf)
          WRITE(lwc,9205) fint(mm)
        enddo
c
        WRITE(lwc,*) ' '
      endif
c
c ***** Output of Some NORMs (NOR_***.txt) *****************************
c   === Headers ===
      if(in.eq.1) then
        WRITE(lwd,9301)
        WRITE(lwd,9302) 0, 0, 0.d0, 0.d0, 0.d0
      endif
 9301 FORMAT('%Step Iter     dfact          unorm          tnorm')
c
c   === Step, # of Iterations, and some Norms ===
      WRITE(lwd,9302) in,itr,dfact,unorm,tnorm
 9302 FORMAT(2i5,3e15.7)
c
c   === Check for first plasticity ===
      if(tene_p.gt.1.d-12 .and. tene_p_prev.le.1.d-12) then
        WRITE(lwd,'(a,i3)') 'First plastic deformation at step: ',in
      endif
      tene_p_prev = tene_p
c
c ***** Output of Energy (ENE_***.txt)*********************************
c   === Headers ===
      if(in.eq.1) then
        WRITE(lwe,9401)
        WRITE(lwe,9402) 0, 0.d0, 0.d0, 0.d0
      endif
c
 9401 FORMAT('%Step      TotalEne.   Elastic Ene.   Plastic Ene.')
c
c   === Step, Total, Elastic, Plastic energy ===
      WRITE(lwe,9402) in,tene_e+tene_p,tene_e,tene_p
c
 9402 FORMAT(i5,3e15.7)
c
c ***** Output of Temperature (TMP_***.txt)*****************************
c   === Hearders ===
        if(in.eq.1) then
          WRITE(lwf,9211)
          do lb=1,lbprt
            WRITE(lwf,9214) (mbprt(i,lb),i=1,2)
          enddo
          WRITE(lwf,*) ' '
        endif
c
c   === Loading Step ===
        WRITE(lwf,9204) in
c
c   === Element Temperature ===
        do lb=1,lbprt
          mm = mbprt(1,lb)
          WRITE(lwf,9205) temp(mm)
        enddo
c   === Total Temperature ===
        WRITE(lwf,9205) (tene_p*10.d0**9)
     &        /(0.427d0*9.8d0/4.18d0*7800.d0*0.11d0*4.18d0*10.d0**3)
c
        WRITE(lwf,*) ' '
c      endif
c
c ***** Output to CML-formatted File (RES_***.cml) *********************
      ip = 0
      do lp=1,lpstp
        if(in.eq.mpstp(lp)) ip = 1
      enddo
c     WRITE(lwa,9000) keyLAS
c
      if(ip.eq.1) then
c ===== /NODAL/ Nodal Value (Displacement and Potential) =====
        WRITE(lwa,9000) key800
        WRITE(lwa,9001) in,m1,m1
        WRITE(lwa,9002) nx
        if(ndf.eq.2) then
          do nn=1,nx
            ux = disp( ndf*(nn-1) +1)
            uy = disp( ndf*(nn-1) +2)
            WRITE(lwa,9003) nn,ux,uy,dum,dum,dum,dum
          enddo
        else
          do nn=1,nx
            ux = disp( ndf*(nn-1) +1)
            uy = disp( ndf*(nn-1) +2)
            uz = disp( ndf*(nn-1) +3)
            WRITE(lwa,9003) nn,ux,uy,uz,dum,dum,dum
          enddo
        endif
c
c ===== /ELMTL/ Elemental Value (Stress, Strain ...etc) =====
c     sigma(1,nel) = sts(1,1)
c     sigma(2,nel) = sts(2,1)
c     sigma(3,nel) = sts(3,1)
c     sigma(4,nel) = sts(2,2)
c     sigma(5,nel) = sts(3,2)
c     sigma(6,nel) = sts(3,3)
        WRITE(lwa,9000) key810
        WRITE(lwa,9001) in,m1,m1
        WRITE(lwa,9002) nelx
        if(ndf.eq.2) then
          do nel=1,nelx
            WRITE(lwa,9003) nel,
     &                      sigma(1,nel),sigma(4,nel),sigma(2,nel),
     &                      epsln(1,nel),epsln(4,nel),epsln(2,nel)
            WRITE(lwa,9004) von(nel),eps(nel),eene(nel)+pene(nel),
     &                      eene(nel),pene(nel),temp(nel)
          enddo
        else
          do nel=1,nelx
            WRITE(lwa,9003) nel,
     &                      sigma(1,nel),sigma(4,nel),sigma(6,nel),
     &                      sigma(5,nel),sigma(3,nel),sigma(2,nel)
            WRITE(lwa,9004) epsln(1,nel),epsln(4,nel),epsln(6,nel),
     &                      epsln(5,nel),epsln(3,nel),epsln(2,nel)
            WRITE(lwa,9004) von(nel),eps(nel), eene(nel)+pene(nel),
     &                      eene(nel),pene(nel),temp(nel)
          enddo
        endif
c
c ===== /ENDOF/ Indicator of the End Of File =====
c     WRITE(lwa,9000) keyEND
c
      endif
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
