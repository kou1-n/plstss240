      subroutine st_gtn(itrmax, idepg,
     &                   prope,  sig,   str, ehist,
c    &                   prope,!plstrg,betaeg, alpeg,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror,  itr , histi )
c
      implicit double precision(a-h,o-z)
c     implicit integer(i-n)
c     implicit none
c
      dimension prope(20)
      dimension sts(3,3),stn(3,3)
      dimension str(3,3),stry(3,3),etri(3,3),
     &          sig(3,3),oun(3,3),sd(3,3),eel(3,3),
     &          plstrg(3,3),edev(3,3)
      dimension ctens(3,3,3,3),DEV(3,3,3,3)
      dimension ehist(20),r(3),test(3,3)
      dimension histi(50)
c   for subroutine
      integer n,m,NP,MP
      PARAMETER(MP=20,NP=20)
      REAL*8 a(NP,NP),b(NP,MP),ai(NP,NP),c(NP,NP),x(NP,MP),bi(NP,MP)

c   for st_gtn      
      REAL*8 dmg,dmged,kkk,kkked,hyds,eqst,omega,aaa,bbb,ccc,ddd,eee
      REAL*8 deltag
c   for Block Newton method
      dimension psig(3,3),xa1(3,3),xa2(3,3),stri(3,3)
      dimension xa1i(3,3),xa2i(3,3),plstrgi(3,3),delstr(3,3)
      REAL*8    deltagi,hydsi,dmgi,kkki,eqpsig
c
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)

      do ll=1,3
        do kk=1,3
          do jj=1,3
            do ii=1,3
              DEV(ii,jj,kk,ll)
     &           =FIT(ii,jj,kk,ll)-(1.d0/3.d0)*
     &             DELTA(ii,jj)*DELTA(kk,ll)
            enddo
          enddo
        enddo
      enddo
c **********************************************************************
c     write(*,*) 'stress'
c
c ***** Load Deformation Histories *************************************
      alpeg = ehist(1)
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          plstrg(jj,ii) = ehist(kk)
        enddo
      enddo
      dmg    = ehist(19)
      kkk    = ehist(20)
c
      deltagi = histi(1)
      hydsi   = histi(2)
      dmgi    = histi(3)
      kkki    = histi(4)
      kk=8
      do jj=1,3
        do ii=1,3
          kk = kk +1
          stri(ii,jj) = histi(kk)
        enddo
      enddo      
      do jj=1,3
        do ii=1,3
          kk = kk +1
          xa1i(ii,jj) = histi(kk)
        enddo
      enddo
      do jj=1,3
        do ii=1,3
          kk = kk+1
          xa2i(ii,jj) = histi(kk)
        enddo
      enddo
c      do jj=1,3
c        do ii=1,3
c          kk = kk+1
c          plstrgi(:,:) = histi(kk)
c        enddo
c      enddo
c ***** Set Material Properties ****************************************
c    --- elastic parameters
c     --- Young modulus
      yng = prope(1)
c     --- Poisson's ratio
      poi = prope(2)
c       --- Lame const. 'mu' = shear modulus
        vmu = yng/(2.d0*(1.d0 +poi))
c       --- Lame const. 'lamuda'
        vlm = poi*yng/((1.d0 +poi)*(1.d0 -2.d0*poi))
c       --- Bulk modulus 'kappa'
        vkp = yng/(3.d0*(1.d0 -2.d0*poi))
c    --- plastic parameters
      yld = prope(10)
      hk  = prope(11)
      hpa = prope(12)
      hpb = prope(13)
c     hpc = prope(14)
      hpd = prope(15)
      fracI = prope(7)
c
c ***** Initialization *************************************************
      deltag=0
      psig   = 0.d0
      ctens = 0.d0 ! ctens(:,:,:,:) = 0.d0
      dmg   = dmg + fracI
      dmged = dmg
      kkked = kkk
c      print '("fracI=",e15.5)',fracI
c      print '("Initial dmg=",e15.5)',dmg
c
c   === Deviatoric Stress => Trial Stress ===
      emean = (str(1,1) +str(2,2) +str(3,3))/3.d0
c
      stry = 2.d0*vmu*(str -emean*DELTA -plstrg)
      !stry(:,:) = 2.d0*vmu*(str(:,:) -emean*DELTA(:,:) -plstrg(:,:))
      stno = 0.d0
      do jj=1,3
        do ii=1,3
          stno = stno +stry(ii,jj)**2
        enddo
      enddo
      stno = dsqrt(stno)
c      
      edev = str -emean*DELTA -plstrg
c
      edevno = 0.d0
      do jj=1,3
        do ii=1,3
          edevno = edevno +edev(ii,jj)**2
        enddo
      enddo
      edevno = dsqrt(edevno)
c      print '("edevno=",e15.5)',edevno
c
      etrs = str(1,1) +str(2,2) +str(3,3)
      hyds = vkp*etrs/3.d0
c      write(*,*) etrs,hyds
c      print '("hyds=",e15.5)',hyds
      oun(:,:) = stry(:,:)/stno
c
      eqst = 0.d0
      do jj=1,3
        do ii=1,3
          eqst = eqst + stry(ii,jj)**2
        enddo
      enddo
      eqst = dsqrt(1.5d0*eqst)
c      print '("eqst",e15.5)',eqst
c
c  === Compute Hardening Function & Yield Function ===
      hard = yld + hpd*hk*kkk
     &     +(hpa -yld) *(1.d0 -dexp(-hpb*kkk))
c      print '("Initial hard",e15.5)',hard
     
      omega = 1.5d0*hyds/hard
      aaa   = dsqrt(1.d0+dmg**2-2.d0*dmg*dcosh(omega))
c      print '("aaa=",e15.5)',aaa
c
      ftreg = eqst-aaa*hard
c     print '("Trial ftreg=",e25.15)',ftreg,dmg
c      write(*,*) 'Trial ftreg=',ftreg
c      write(*,*) 'Trial omega=',omega
c      write(*,*) 'Trial aaa=',aaa
c      write(*,*) eqst,aaa*hard
c      write(*,*) aaa,hard
c
c  ===== PLASTIC CASE
c     if(ftreg.gt.0.d0) then
      if(ftreg.gt.-ctol) then
c             write(*,*) 1
        idepg = 1
c        print *,"PLASTIC"
        if(itr==1) then
          deltag = 0.d0
          hyds   = vkp*etrs/3.d0
          dmg    = ehist(19)
          kkk    = ehist(20)

        else
c     --- Withdrow R_q          
          b(1,1) = histi(5)
          b(2,1) = histi(6)
          b(3,1) = histi(7)
          b(4,1) = histi(8)   
c         
          a=0.d0
          x=0.d0
          ai=0.d0
          c=0.d0
          n=4
          m=1
c          
          hard = yld + hpd*hk*kkki
     &       +(hpa -yld) *(1.d0 -dexp(-hpb*kkki))
          omega  = 1.5d0*hydsi/hard
          aaa    = dsqrt(1.d0+dmgi**2-2.d0*dmgi*dcosh(omega))
          dhard  = hk + (hpa -yld) *hpb*dexp(-hpb*kkki)
          ddhard = -(hpa -yld) *hpb**2*dexp(-hpb*kkki)
c 
          a(1,1) = -3.d0*vmu
          a(1,2) = (3.d0*dmgi*dsinh(omega))/(2.d0*aaa)
          a(1,3) = -hard/aaa*(dmgi-dcosh(omega))
          a(1,4) = -(aaa+omega/aaa*dmgi*dsinh(omega))*dhard
c
          a(2,1) = dmgi/(2.d0*aaa)*dsinh(omega)
          a(2,2) = 1.d0+3.d0*dmgi*deltagi/(4.d0*aaa*hard)*(dcosh(omega)
     &            +(dmgi/aaa*dsinh(omega))**2)
          a(2,3) = deltagi/(2.d0*aaa)*(1.d0-dmgi/aaa**2*
     &            (dcosh(omega)-dmgi))*dsinh(omega)
          a(2,4) = -0.75d0*deltagi*hydsi*dmgi/(aaa*hard**2)*dhard*
     &            ((dmgi/aaa*dsinh(omega))**2+dcosh(omega))
c
          a(3,1) = -1.5d0*(dmgi-dmgi**2)/aaa*dsinh(omega)
          a(3,2) = -2.25d0*deltagi*(dmgi-dmgi**2)/(aaa*hard)*
     &            (dmgi*(dsinh(omega)/aaa)**2+dcosh(omega))
          a(3,3) = 1.d0+1.5d0*deltagi/aaa*((dmgi-dmgi**2)*
     &            (dmgi-dcosh(omega))/aaa**2
     &            -(1.d0-2.d0*dmgi))*dsinh(omega)
          a(3,4) = 2.25d0*deltagi*(dmgi-dmgi*2)*hyds/(aaa*hard**2)*
     &            dhard*(dcosh(omega)+dmgi*(dsinh(omega)/aaa)**2)
c
          a(4,1) = -dhard*(omega*dmgi/aaa*dsinh(omega)+aaa)
          a(4,2) = -deltagi*dhard*2.25d0*hydsi*dmgi/(aaa*hard**2)*
     &            (dcosh(omega)+dmgi/aaa**2*(dsinh(omega))**2)
          a(4,3) = -deltagi/aaa*dhard*
     &            (dmgi-dcosh(omega)-(omega*dmgi*(dmgi-dcosh(omega))/
     &            aaa**2-omega)*dsinh(omega))
          a(4,4) = 1.d0-deltagi*ddhard*(aaa+omega*dmgi/aaa*dsinh(omega))
     &            +deltagi*dhard**2*omega**2*dmgi/(aaa*hard)*
     &            (dcosh(omega)+dmgi/aaa**2*(dsinh(omega))**2)   
c
          do 203 l=1,n
            do 201 k=1,n
              ai(k,l)=a(k,l)
  201       continue
            do 202 k=1,m
              x(l,k)=b(l,k)
  202       continue
  203     continue

c     --- Use Subroutine "gaussj"          
          call dgaussj(ai,n,NP,x,m,MP) 
c          print '("ai(1,1)=",e15.5)',ai(1,1) 
c          print '("ai(2,1)=",e15.5)',ai(2,1) 
c          print '("ai(1,2)=",e15.5)',ai(1,2) 
c          print '("ai(2,2)=",e15.5)',ai(2,2)       

c     --- Update variables
          c = matmul(ai,b)
c          delstr(:,:)=str(:,:)-stri(:,:)
c          write(*,*) delstr
          do jj=1,3
            do ii=1,3
              c(1,1) = c(1,1)-(ai(1,1)*xa1i(ii,jj)+ai(1,2)*xa2i(ii,jj))*
     &              (str(ii,jj)-stri(ii,jj))
              c(2,1) = c(2,1)-(ai(2,1)*xa1i(ii,jj)+ai(2,2)*xa2i(ii,jj))*
     &              (str(ii,jj)-stri(ii,jj))
            enddo
          enddo
          deltag  = deltagi + c(1,1)
          hyds    = hydsi   + c(2,1)
          dmg     = dmgi    + c(3,1)
          kkk     = kkki    + c(4,1)
        endif  
c        print '("c(1)=",e25.15)',c(1,1)
c        print '("c(2)=",e25.15)',c(2,1)
c        print '("c(3)=",e25.15)',c(3,1)
c        print '("c(3)=",e25.15)',c(4,1)
c     --- Store Updated variables for next iteration
        histi(1) = deltag
        histi(2) = hyds
        histi(3) = dmg
        histi(4) = kkk
c        print '("c(1,1)=",e25.15)',c(1,1)
c        print '("histi(2)=",e25.15)',histi(2)
c        print '("histi(3)=",e25.15)',histi(3)
c        print '("histi(4)=",e25.15)',histi(4)
c     --- Finish  Update variables     
c
        hard    = yld + hk*kkk
     &             +(hpa -yld) *(1.d0 -dexp(-hpb*kkk))
        dhard  = hk + (hpa -yld) *hpb*dexp(-hpb*kkk)
        ddhard = -(hpa -yld) *hpb**2*dexp(-hpb*kkk)
c          print '("Update hard",e15.5)',hard
        omega   = 1.5d0*hyds/hard
        aaa   = dsqrt(1+dmg**2-2.d0*dmg*dcosh(omega))

        ftreg = (eqst-3.d0*vmu*deltag)-aaa*hard
c        write(*,*) it
c        print '("Update ftreg=",e25.15)',ftreg
c        r(1) = hyds-vkp*emean+deltag/aaa*dmg*dsinh(omega)
        r(1) = hyds-vkp*etrs/3+deltag/aaa*dmg*dsinh(omega)
        r(2) = dmg - dmged -deltag*(dmg-dmg**2)*
     &             1.5d0/aaa*dsinh(omega)
        r(3) = kkk - kkked -deltag*(omega/aaa* 
     &             dmg*dsinh(omega)+aaa)*dhard 
        histi(5) = -ftreg
        histi(6) = -r(1)
        histi(7) = -r(2)
        histi(8) = -r(3)
c        print '("r(1)=",e25.15)',r(1)
c        print '("r(2)=",e25.15)',r(2)
c        print '("r(3)=",e25.15)',r(3)
c     --- update equivalent plastic strain "alpha"
c     --- compute the plastic strain etc.
        alpeg  = alpeg +dsqrt(2.d0/3.d0)*deltag
        plstrg(:,:)=plstrg(:,:) +320.d0*deltag*(oun(:,:)
     &                 +1.5d0/aaa*dmg*dsinh(omega)*DELTA(:,:))
c   
c     --- update consititutive tensor: "ctens"
        bbb = 2.d0*vmu*(1.d0 -3.d0*vmu*deltag/eqst)
        ccc = 12.d0*vmu**3*deltag/(eqst**2*edevno)
        ddd = dsqrt(6.d0)/edevno*(1.d0-3.d0*vmu*deltag/eqst)
        eee = dsqrt(1.5d0)*12.d0*vmu**3*deltag/eqst**2

        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
                ctens(ii,jj,kk,ll)
     &            = bbb*DEV(ii,jj,kk,ll)
     &              +ccc*edev(ii,jj)*edev(kk,ll)     
     &              -6.d0*vmu**2/eqst*edev(ii,jj)*
     &              (-ai(1,1)*(ddd+eee)*edev(kk,ll)
     &              +ai(1,2)*vkp*DELTA(kk,ll))
     &              +DELTA(ii,jj)*(-ai(2,1)*(ddd+eee)*edev(kk,ll)
     &              +ai(2,2)*vkp*DELTA(kk,ll))
              enddo
            enddo
          enddo
        enddo
        xa1(:,:) = (ddd+eee)*edev(:,:)
        xa2(:,:) = -vkp*DELTA(:,:)
        kk=8
        do jj=1,3
          do ii=1,3
            kk = kk +1
            histi(kk) = str(ii,jj)
          enddo
        enddo      
        do jj=1,3
          do ii=1,3
            kk = kk +1
            histi(kk) = xa1(ii,jj)
          enddo
        enddo
        do jj=1,3
          do ii=1,3
            kk = kk+1
            histi(kk) = xa2(ii,jj)
          enddo
        enddo
c        do  l=1,n
c          do  k=1,n
c            ai(k,l)=a(k,l)
c          enddo
c          do  k=1,m
c            x(l,k)=b(l,k)
c          enddo
c        enddo
c     --- compute the Stress and Pseud Stress
c        psig(:,:) = -ftreg*(-ai(1,1)*6.d0*vmu**2/eqst*edev(:,:)
c      &                          +ai(2,1)*DELTA(:,:))
c     &              -r(1)*(-ai(1,2)*6.d0*vmu**2/eqst*edev(:,:)
c     &                          +ai(2,2)*DELTA(:,:))
c        print '("rftr=",e25.15)',ftreg        
c        print '("r(1)=",e25.15)',r(1)
c        print '("r(2)=",e25.15)',r(2)
c        print '("r(3)=",e25.15)',r(3)
        psig(:,:) = -ftreg*(-ai(1,1)*6.d0*vmu**2/eqst*edev(:,:)
     &                          +ai(2,1)*DELTA(:,:))     
        eqpsig = 0.d0
        do jj=1,3
          do ii=1,3
            eqpsig = eqpsig + psig(ii,jj)**2
          enddo
        enddo
        eqpsig = dsqrt(1.5d0*eqpsig)
c        print '("eq_psig=",e25.15)',eqpsig
        sig(:,:)  = stry(:,:) -6.d0*vmu**2*deltag*edev(:,:)
     &                    /eqst+vkp*etrs*DELTA(:,:) + psig(:,:)
c        print '("vkp*etrs   =",e25.15)',vkp*etrs
c        print '("vkp*etrs/3 =",e25.15)',vkp*etrs/3
c        print '("hyds       =",e25.15)',hyds
c        sig(:,:)  = stry(:,:) -6.d0*vmu**2*deltag*edev(:,:)
c     &                    /eqst+hyds*DELTA(:,:) + psig(:,:)
c  ===== ELASTIC CASE:
c          Identify the Trial Stress as Actual One =====
      else
        idepg = 0
c        print *,"ELASTIC"
c
        etrs = str(1,1) +str(2,2) +str(3,3)
c       if(dabs(etrs).le.1.0d-16) etrs = 0.d0
        sig(:,:) = stry(:,:) +vkp*etrs*DELTA(:,:)
c       write(*,*) etrs,ctol
c       write(*,*) vkp*etrs,'aaa'
c
c     --- donot update "ctens" ---
        do ll=1,3
          do kk=1,3
            do jj=1,3
              do ii=1,3
c               ctens(ii,jj,kk,ll) = vlm*DELTA(ii,jj)*DELTA(kk,ll)
c    &                              +2.d0*vmu*FIT(ii,jj,kk,ll)
                ctens(ii,jj,kk,ll)
     &            = vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &              +2.d0*vmu*( FIT(ii,jj,kk,ll)
     &                       -(1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
              enddo
            enddo
          enddo
        enddo
c
      endif
c
c   === Compute von Mieses Stress  ===
      smean = (sig(1,1) +sig(2,2) +sig(3,3))/3.d0
c      print '("smean=",e15.5)',smean
c     write(*,*) smean
c     write(*,*) sig
c
      sd(:,:) = sig(:,:) -smean*DELTA(:,:)
c
      vons = 0.d0
      do jj=1,3
        do ii=1,3
          vons = vons +sd(ii,jj)*sd(ii,jj)
        enddo
      enddo
      vons = dsqrt(1.5d0*vons)
c
c   === Compute the Energy Density ===
c     --- separate the strain into the elastic & plastic
      eel(:,:) = str(:,:) -plstrg(:,:)
c
c     --- compute the elastic strain energy density
      e_dns = 0.d0
      do jj=1,3
        do ii=1,3
          e_dns = e_dns +eel(ii,jj)*sig(ii,jj)
        enddo
      enddo
      e_dns = e_dns*0.5d0
c
c     --- compute the plastic strain energy density
      p_dns = 0.d0
      do jj=1,3
        do ii=1,3
          p_dns = p_dns+deltag*oun(ii,jj)*sig(ii,jj)
        enddo
      enddo
c
c ***** Store Deformation Histories ************************************
      ehist(1) = alpeg
      kk = 1
      do ii=1,3
        do jj=1,3
          kk = kk+1
          ehist(kk) = plstrg(jj,ii)
        enddo
      enddo
c      write(*,*) plstrg
c      print '("final dmg=",e25.15)',dmg
c      print '("KAKUNOU dmg=",e25.15)',dmg - fracI
      ehist(19) = dmg - fracI
      ehist(20) = kkk
c      print '("final hard=",e25.15)',hard
c      print '("final kkk=",e25.15)',kkk
c
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
