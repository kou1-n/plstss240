      subroutine analys(  NGSK, NGSKo,   neq,  neqo,  iRCM,
     &                     ijk,   mpe,  mang,  nfix,  ndir,
     &                  npload,nsload,nbload,  ijkl,   ncp,
     &                    mdof,  idof, jdiag,iw_neq,index0,
     &                  index1,jcolmn,  list,ncolmn, mpstp,
     &                   muprt, mfprt, mnprt, msprt, mbprt,
     &                   idep0, idep1, melem, matid,
c
     &                     xyz,  prop, angle,  vfix,vpload,
     &                  vsload,vbload,   foc,  disp,   del,
     &                   sigma, epsln,   von,dw_sol,    du,
     &                      u0,    u1,   res,  fint,   eps,
     &                    finc,  pene,  eene,  temp, tempd,
     &                  dndx_g, det_g,ctensg,dhist0,dhist1,
     &                      sk,
     &                  secslv,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      real*4 sec,sec0,sec1
c
      dimension ijk(node,nelx)
      dimension mpe(nelx),mang(nelx),melem(nelx)
      dimension idep0(ngaus,nelx),idep1(ngaus,nelx)
      dimension nfix(nspc)
      dimension ndir(6,nspc)
      dimension npload(npoin)
      dimension nsload(npres)
      dimension nbload(nbody)
      dimension ijkl(nsn,npres)
      dimension ncp(neq),mdof(neq),idof(neq),
     &          iw_neq(neq)
      dimension jdiag(neq+1)
      dimension index0(neqo+1)
      dimension index1(NGSKo)
      dimension jcolmn(nsolvr+1),ncolmn(nodsk)
      dimension list(nmax,nx)
      dimension mpstp(lpstp)
      dimension muprt(2,luprt)
      dimension mfprt(2,lfprt)
      dimension mnprt(2,lnprt)
      dimension msprt(2,lsprt)
      dimension mbprt(2,lbprt)
      dimension matid(lmat)
c
      dimension xyz(3,nx)
      dimension prop(20,lmat)
      dimension angle(3,lang)
      dimension vfix(6,nspc)
      dimension vpload(6,npoin)
      dimension vsload(4,npres)
      dimension vbload(3,nbody)
      dimension dw_sol(neqsol)
      dimension foc(neq),disp(neq),del(neq),
     &          du(neq),u0(neq),u1(neq),res(neq),fint(neq)
      dimension sigma(6,nelx),epsln(6,nelx)
      dimension von(nelx),eps(nelx),pene(nelx),eene(nelx),temp(nelx),
     &          tempd(nelx)
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
      dimension finc(nstep)
      dimension sk(NGSK)
      dimension dhist0(20,ngaus,nelx),dhist1(20,ngaus,nelx)
      dimension histi0(50,ngaus,nelx)
c
c     --- arrays and parameters for PARDISO solver ---
      integer*8 pt(64)
      integer iparm(64)
      integer maxfct,  mnum, mtype, phase,  nrhs
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
      common /tvalu/ ctol,stol,switch
      common /print/ lpstp,luprt,lfprt,lnprt,lsprt,lbprt
      common /debug_info/ nel_current, ig_current, lstep_current
c
c **********************************************************************
      arc = 0.d0
c
c ****** Initialization ************************************************
      index0 = 0 ! index0(:) = 0
      index1 = 0 ! index1(:) = 0
c
      u0 = 0.d0 ! u0(:) = 0.d0
c
c ****** Read All Data from CML-formatted File *************************
c     CALL CPU_TIME(  sec0)
      CALL redcml(   neq,  neqm,
     &               ijk,   mpe,  mang,  nfix,  ndir,
     &            npload,nsload,nbload,  ijkl, melem,
     &               ncp,  mdof,  idof, mpstp, muprt,
     &             mfprt, mnprt, msprt, mbprt, matid,
     &            arcini,
     &               xyz,  prop, angle,  vfix,vpload,
     &            vsload,vbload,  finc,
     &            ierror )
      if(ierror.ne.0) RETURN
c     CALL CPU_TIME(  sec1)
c     WRITE(*,102) sec1 -sec0
c
c     do ne=1,neq
c       write(*,'(5i5)') ne, ncp(ne), mdof(ne), idof(ne)
c     enddo
c
c ****** Re-Ordering for Reduction of Skyline Hight ********************
      if(iRCM.eq.1) then
        do nn=1,5
          READ(lrb,*)
        enddo
        do ne=1,neq
          READ(lrb,'(3i8)') nn,idof(ne),mdof(ne)
        enddo
        CLOSE(lrb)
      endif
c
c ****** Pre-Processing for Solvers for Un-Constrained DOF *************
c     WRITE(*,105)
c     CALL CPU_TIME(  sec0)
      if( (isolvr.eq.1).or.(isolvr.eq.2) ) then
c   ===== PARDISO Solver & RCI CG Solver =====
        NGSKm = NGSK -NGSKo
c     --- Define the Index Arrays ---
c           ( ROW major UPPER triangle & LOWER Triangle )
        CALL rowupr(    nx,  nelx,   ndf,  node,
     &                 neq,  neqm,  NGSK, NGSKm,
     &                 ijk,  mdof,  idof, jdiag,jcolmn,
     &              iw_neq,  list,ncolmn,index0,index1,
     &              ierror )
        if(ierror.ne.0) RETURN
c
      elseif(isolvr.eq.0) then
c   ===== SKYLINE Solver =====
c     --- Define the Index Arrays for NON-Constrained DOF ---
        CALL presky(    nx,  nelx,   ndf,  node,
     &                 neq,  neqm,  NGSK, NGSKm,
     &                 ijk,  mdof,  idof, jdiag,
     &              iw_neq,
     &              ierror )
        if(ierror.ne.0) RETURN
c
c     --- Define the Index Arrays for Constrained DOF ---
         CALL inform(    nx,  nelx,   ndf,  node,
     &                  neq,  neqm, NGSKo,  neqo,
     &                  ijk,  mdof,  idof,
     &               iw_neq,index0,index1,
     &               ierror )
        if(ierror.ne.0) RETURN
      endif
c      write(*,*) ' ',NGSKo
c      write(*,*) index1
c      write(*,*) ' ',neqo
c      write(*,*) index0
c      write(*,*) ' ',nsolvr
c      write(*,*) jcolmn
c     CALL CPU_TIME(  sec1)
c     WRITE(*,102) sec1-sec0
c
c ****** Compute Total External Force Vector etc. **********************
      CALL forces(   neq,
     &               ijk,  nfix,  ndir,
     &            npload,nsload,  ijkl,
     &               xyz,  vfix,vpload,
     &            vsload,   foc,  disp,
     &            ierror )
      if(ierror.ne.0) RETURN
c
c     do ne=1,neq
c       write(*,*) ne,disp(ne),foc(ne)
c     enddo
      l_beg = 1
      l_end = nstep
c
      dfact = 0.d0
      dfact0= 0.d0
c
c ***** Initial Settings for Computation *******************************
      CALL initia(  NGSK,   neq,  neqm, NGSKo,  neqo,
     &               ijk,   mpe,  mang,  mdof,  idof,
     &             jdiag,index0,index1,jcolmn, melem,
     &               xyz,  prop, angle,
     &            dndx_g, det_g,ctensg,
     &            ierror )
      if(ierror.ne.0) RETURN
c
c ****** Initialize Deformation Histories ******************************
      tene_p = 0.d0
      dhist0 = 0.d0 ! dhist0(:,:,:) = 0.d0
      dhist1 = 0.d0 ! dhist1(:,:,:) = 0.d0
      histi0 = 0.d0
c
      idep0 = 0 ! idep0(:,:) = 0
      idep1 = 0 ! idep1(:,:) = 0
c
c     ----- Set Initial Temperature -----
      temp = 0.d0 ! temp(:) = 0.d0
c
c     write(*,*) l_beg,l_end
c ***************************************
c ***** START Incremental Procedure *****
c ***************************************
      do 1000 in=l_beg,l_end
c
        WRITE(*,8001) in,l_end
c       --- Set current loading step for debug ---
        lstep_current = in
c
c     ----- Control Loading Parameter -----
        df = finc(in)
c
 1400 CONTINUE
c
c     ----- Initilize Displacement Increment -----
        du = 0.d0 ! du(:) = 0.d0
c
c     ----- Set Iinitial Load Increment in this Loading Step -----
        res = df*foc ! res(:) = df*foc(:)
c       write(*,*) res
c
        itrmax = 20
c       /////////////////////////////////////////
c      //// START Newton-Raphson Iteration /////
c     /////////////////////////////////////////
        do 2000 itr=1,itrmax
c
          WRITE(*,8002) itr
c
c ****** Assemble the Global Stiffness Matrix: SK **********************
c         WRITE(*,101)
c         CALL CPU_TIME(  sec0)
          CALL assemb(  NGSK,   neq,  neqm, NGSKo,  neqo,
     &                   ijk,   mpe,  mang,  mdof,  idof,
     &                 jdiag,index0,index1,jcolmn, melem,
     &                   xyz,  prop, angle,    sk,
     &                dndx_g, det_g,ctensg,
     &                ierror )
          if(ierror.ne.0) RETURN
c         write(*,*) sk
c         stop
c         CALL CPU_TIME(  sec1)
c         WRITE(*,102) sec1 -sec0
c       write(*,*) ''
c       write(*,*) dndx_g
c       write(*,*) det_g(1,1)
c       write(*,*) ''
c       write(*,*) 'CLEAR ASSEMB'
c
c ****** Convert Enforced Displacement to the Equivalent Force *********
          if(itr.eq.1) then
            CALL constr(  NGSK,   neq,  neqm, NGSKo,  neqo,
     &                    mdof,  idof, jdiag,index0,index1,
     &                      df,
     &                     res,  disp,   del,    sk,
     &                  ierror )
            if(ierror.ne.0) RETURN
c
          else
c       ( Re-order the Total Force Vector )
            do ne=1,neq
              m_ne = mdof(ne)
              del(m_ne) = res(ne)
c             res(ne) = 0.d0
            enddo
          endif
c
          res = 0.d0 ! res(:) = 0.d0
c
c         write(*,*) m_ne
c         do ng=1,NGSK
c           write(*,'(i5,e15.5)') ng,sk(ng)
c         enddo
c         do ne=1,neq
c           res(ne) = 0.d0
c         enddo
c         do ne=1,neq
c           write(*,*) ne,del(ne)
c         enddo
c
c ****** Solve Linear Equations ****************************************
c         WRITE(*,103)
c         CALL CPU_TIME(  sec0)
c         sec_10 = dble(sec0)
c
          if(isolvr.eq.1) then
c       ===== PARDISO Solver =====
            if(jsol.eq.0) then
c           --- Just for First Computation 
c                 Analysis, Symbolic fact., Numerical fact. & Solve
c                 Controling parameters are set in the following
              CALL pars00(  neqm,  NGSK,nsolvr, msol1,
     &                     jdiag,jcolmn,iw_neq,
     &                        sk,   del,dw_sol,
     &                        pt, iparm,
     &                    maxfct,  mnum, mtype, phase,  nrhs,
     &                    ierror )
              jsol = 1
            else
c           --- for the latter computation
c                 Numerical factorization & Solve
              phase = 23
              CALL parsol(  neqm,  NGSK,nsolvr, msol1,
     &                     jdiag,jcolmn,iw_neq,
     &                        sk,   del,dw_sol,
     &                        pt, iparm,
     &                    maxfct,  mnum, mtype, phase,  nrhs,
     &                    ierror )
            endif
c
          elseif(isolvr.eq.2) then
c       ===== RCI CG Solver =====
            tolsol = 1.d-10
            CALL pcgsol(  neqm,  NGSK,nsolvr, msol1,
     &                   jdiag,jcolmn,
     &                    stol,    sk,   del,dw_sol(1),
     &                  dw_sol(neq+1), dw_sol(2*neq+1),
     &                  ierror )
c
          elseif(isolvr.eq.0) then
c       ===== SKYLINE Solver =====
            ntt = 0
            CALL skylin(    sk,   del, jdiag,     1,  neqm,   ntt)
          endif
c
c         CALL CPU_TIME(  sec1)
c         sec_11 = dble(sec1)
c         secslv = sec_11-sec_10
c         WRITE(*,102),sec1 -sec0
c
c ****** Update the Displacement and Potentials ************************
          CALL update(   neq,
     &                  idof,
     &                   del,    du,    u0,    u1,
     &                ierror )
          if(ierror.ne.0) RETURN
c         write(*,*) u1
c         write(*,*) '  Total displacement '      !hs
c         do nn=1,nx
c           write(*,'(i5,1p3e12.4)') nn,(u1(ndf*(nn-1)+kk),kk=1,ndf) !hs
c           write(*,*) nn,(u1(ndf*(nn-1)+kk),kk=1,ndf) !hs
c         enddo
c         stop
c         write(*,*) '  Total displacement Modification' !hs
c         do nn=1,nx
c           write(*,'(i5,1p3e12.4)') nn,(del(ndf*(nn-1)+kk),kk=1,ndf)
c         enddo
c
c ****** Compute the Internal Forces etc *******************************
c         WRITE(*,104)
c         CALL CPU_TIME(  sec0)
          CALL postpr(   neq,
     &                   ijk,   mpe,  mang, idep1, melem,
     &                 matid,
     &                   xyz,  prop, angle,    u1, sigma,
     &                 epsln,   von,  fint,dhist0,dhist1,
     &                   eps,  pene,  eene,
     &                tene_e,dene_p,  temp, dtemp, tempd,
     &                dndx_g, det_g,ctensg, g_norm,
     &                ierror, itr , histi0)
          if(ierror.ne.0) RETURN
c         CALL CPU_TIME(  sec1)
c         WRITE(*,102) sec1 -sec0
c
c ****** Check Occurrence of Unloading *********************************
          munld = 0
          do nel=1,nelx
            do ig=1,ngaus
              if( (idep0(ig,nel).eq.1).and.(idep1(ig,nel).eq.0) )then
                munld = munld +1
                idep0(ig,nel) = -1
c               write(*,*) munld,ig,nel
              endif
            enddo
          enddo

C          if(munld.gt.0) then
C            itres = itres +1
C            CALL restor(   neq,
C     &                      u0,    u1,dhist0,dhist1,
C     &                  ierror )
C            WRITE(*,8005) itres
C            GOTO 1400
C          endif
C          itres = 0
c
c       ===== Compute the Residual Forces & Their Norms =====
          fnorm = 0.d0
          rnorm = 0.d0
          tnorm = 0.d0
          unorm = 0.d0
          anorm = 0.d0
          dfact = dfact0 +df
c         write(*,*) dfact
          do ne=1,neq
            fnorm = fnorm +(df*foc(ne))**2
            unorm = unorm +(u1(ne))**2
            tnorm = tnorm +(dfact*foc(ne))**2
c           write(*,'(i5,3e20.9)') ne,fint(ne),dfact*foc(ne),res(ne)
            if(ncp(ne).eq.0) then
              res(ne) = dfact*foc(ne) -fint(ne)
              rnorm = rnorm +res(ne)**2
c           write(*,'(i5,3e20.9)') ne,fint(ne),dfact*foc(ne),res(ne)
            endif
          anorm = anorm +fint(ne)**2
          enddo
          unorm = dsqrt(unorm)
          tnorm = dsqrt(tnorm)
c
c         do ne=1,neq
c           write(*,'(2i5,3e25.16)') ne,ncp(ne),res(ne)
c         enddo
c
c       ===== Check the Convergence =====
          if(dsqrt(fnorm).lt.ctol) then
c           fnorm = 1.d0
            fnorm = anorm
          endif
c
          rbf = dsqrt(rnorm/fnorm)
          WRITE(*,8004) df,dfact,arc
c
c       --- Check if using Block Newton method (MATYPE=5) ---
          isbnm = 0
          do imat=1,lmat
            if(matid(imat).eq.5) isbnm = 1
          enddo
c
          if(isbnm.eq.1) then
c         --- Extended output for Block Newton method ---
            WRITE(*,8006) dsqrt(rnorm),dsqrt(fnorm),rbf,g_norm
          else
c         --- Standard output for other methods ---
            WRITE(*,8003) dsqrt(rnorm),dsqrt(fnorm),rbf
          endif
c
          if(rbf.lt.ctol) then
            GOTO 1500
          endif
c
 2000   CONTINUE
        WRITE(*,*) 'Iteration Does NOT CONVERGED!!'
        ierror = 20
        RETURN
c       //////////////////////////////////////////
c      //// END of Newton-Raphson Iteration /////
c     //////////////////////////////////////////
c
 1500   CONTINUE
c
c ****** Update & Store the Deformation Histories **********************
        histi0 = 0.d0        
        tene_p = tene_p +dene_p
        temp = temp +tempd ! temp(:) = temp(:) +tempd(:)
c
        CALL stored(   neq,
     &              dfact0, dfact,
     &               idep0, idep1,
     &                  u0,    u1,dhist0,dhist1,
     &              ierror )
c
c ****** Output the Computed Results ***********************************
c       WRITE(*,106)
c       CALL CPU_TIME(  sec0)
        CALL output(    in,   itr,   neq,
     &               mpstp, muprt, mfprt, mnprt, msprt,
     &               mbprt, dfact, unorm, tnorm,
     &                  u1, sigma, epsln,   von,  fint,
     &                 eps,  pene,  eene,tene_e,tene_p,
     &                temp,
     &              ierror )
        if(ierror.ne.0) RETURN
c       CALL CPU_TIME(  sec1)
c       WRITE(*,102) sec1 -sec0
c
c ****************************************
c ***** END of Incremental Procedure *****
c ****************************************
 1000 continue
c
c ****** Release the Internal Memory used by PARDISO Solver ************
      CALL pars99(  neqm,  NGSK,nsolvr, msol1,
     &             jdiag,jcolmn,iw_neq,
     &                sk,   del,dw_sol,
     &                pt, iparm,
     &            maxfct,  mnum, mtype, phase,  nrhs,
     &            ierror )
c
c **********************************************************************
  101 FORMAT(5x,'Assembling ....................',$)
c 102 FORMAT(' done!!',/)
  102 FORMAT(' done!!',3x,'(',e12.5,' sec )',/)
  103 FORMAT(5x,'Solving Linear Equations ......',$)
  104 FORMAT(5x,'Post-Processing ...............',$)
  105 FORMAT(5x,'Making Indexed for "SK" .......',$)
  106 FORMAT(5x,'Writing the Results ...........',$)
c
 8001 FORMAT (/,
     &       '***** Loading Step : ',i5,'/',i5,' *****')
c
 8002 FORMAT('  Iteration : ',i3)
c
 8003 FORMAT('      rnorm :',e12.5,', fnorm : ',e12.5,
     &       ', RESIDUAL : ',e12.5)
c
 8004 FORMAT('         df :',e12.5,', dfact : ',e12.5,
     &       ',      arc : ',e12.5 )
c
 8005 FORMAT(5x,'##### Unloading has happened !! ',i3,' #####',/,
     &       8x,'Computation will REstart with elastic stiffness')
c
 8006 FORMAT('      rnorm :',e12.5,', fnorm : ',e12.5,/,
     &       '   ||Rf||/||F|| : ',e12.5,', ||Rg|| : ',e12.5,
     &       ' [Block Newton]')
c **********************************************************************
c **********************************************************************
      RETURN
      END

