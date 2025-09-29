      subroutine pres2d(   nel, ijk_e,ijkl_e,
     &                      xe,    fe, vsl_e,
     &                  incomp,               !hsugiyama
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension ijk_e(8)
      dimension ijkl_e(4)
      dimension ijkl_t(2,2)
      dimension ijkl_n(3)
c
      dimension xe(3,8)
      dimension fe(32)
      dimension fe1(32)
      dimension vsl_e(4)
c
c     dimension wg(2),xg(2),sh(2),dn(2),gs(2)
      dimension wg(3),xg(3),sh(3),dn(3),gs(2),inp(2)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
c **********************************************************************
c
c ***** Initialization *************************************************
      fe = 0.d0 !fe(:) = 0.d0
c
      if(node.le.4) then
c **********************************************************************
c ***** Line Integration for 2D, Linear Elements ***********************
c **********************************************************************
c          <1> ------------------ <2>
c     ----- Set Weight  for Numerical Integration -----
        wg(1) = 1.d0
        wg(2) = 1.d0
        xg(1) = -1.d0/dsqrt(3.d0)
        xg(2) =  1.d0/dsqrt(3.d0)
c
c     ----- Find local connectivities -----
        do ii=1,2
          inn = ijkl_e(ii)
          do no=1,node
            if(inn.eq.ijk_e(no)) then
              ijkl_n(ii) = no
            endif
          enddo
        enddo
c       write(*,*) ijkl_n
c
c     ----- Numerical Integration (line Int.) -----
        do 100 ig=1,2
          xl = xg(ig)
          wx = wg(ig)
c
c       === Shape Functions and Their Derivative ===
          sh(1) = 0.5d0*(1.d0 -xl)
          sh(2) = 0.5d0*(1.d0 +xl)
          dn(1) = -0.5d0
          dn(2) =  0.5d0
c
c       === Compute the Length of THIS Edge ===
          gs = 0.d0  !gs(:) = 0.d0
          do in=1,2
            no = ijkl_n(in)
            do ii=1,2
              gs(ii) = gs(ii) +dn(in)*xe(ii,no)
            enddo
          enddo
          gsurf = dsqrt( gs(1)**2 +gs(2)**2 )*wx
c
c       === Compute the Elemental Force Vector ===
          do in=1,2
            ia = ijkl_n(in)
            do nd=1,2
              iai = 2*(ia -1) +nd
              piai = vsl_e(nd+1)
              do jb=1,2
                fe(iai) = fe(iai) +sh(in)*sh(jb)*piai*gsurf
              enddo
            enddo
c
          enddo
c
  100   CONTINUE
c
      else
c **********************************************************************
c ***** Line Integration for 2D, Parabolic Elements ********************
c **********************************************************************
c          <1> ---------  <3> ---------- <2>
        if(incomp.eq.2) then
c     ===== P1-iso-P2/P0 elements =====
c           (Combination of Linear Elements)
c
c       ----- Set Weight  for Numerical Integration -----
          wg(1) = 1.d0
          wg(2) = 1.d0
          xg(1) = -1.d0/dsqrt(3.d0)
          xg(2) =  1.d0/dsqrt(3.d0)
c
c       ----- Update local connectivity -----
          inn = ijkl_e(1)
          do no=1,node
            if(inn.eq.ijk_e(no)) then
              mp = no +node/2
              jmm = ijk_e(mp)
            endif
          enddo
          ijkl_t(1,1) = ijkl_e(1)
          ijkl_t(2,1) = jmm
          ijkl_t(1,2) = jmm
          ijkl_t(2,2) = ijkl_e(2)
c
c       ----- Loop for Segments -----
          do 300 igg=1,2
c
c         ----- Find local connectivities -----
            do ii=1,2
              inn = ijkl_t(ii,igg)
              do no=1,node
                if(inn.eq.ijk_e(no)) then
                  ijkl_n(ii) = no
                endif
              enddo
            enddo
c           write(*,*) ijkl_n
c
c         ----- Numerical Integration (line Int.) -----
            do 350 ig=1,2
              xl = xg(ig)
              wx = wg(ig)
c
c           === Shape Functions and Their Derivative ===
              sh(1) = 0.5d0*(1.d0 -xl)
              sh(2) = 0.5d0*(1.d0 +xl)
              dn(1) = -0.5d0
              dn(2) =  0.5d0
c
c           === Compute the Length of THIS Edge ===
              gs = 0.d0  !gs(:) = 0.d0
              do in=1,2
                no = ijkl_n(in)
                do ii=1,2
                  gs(ii) = gs(ii) +dn(in)*xe(ii,no)
                enddo
              enddo
c
              gsurf = dsqrt( gs(1)**2 +gs(2)**2 )*wx
c
c           === Compute the Elemental Force Vector ===
              do in=1,2
                ia = ijkl_n(in)
                do nd=1,2
                  iai = 2*(ia -1) +nd
                  piai = vsl_e(nd+1)
                  do jb=1,2
                    fe(iai) = fe(iai) +sh(in)*sh(jb)*piai*gsurf
                  enddo
                enddo
              enddo
c
 350        CONTINUE
c
 300      CONTINUE
c
        else
c     ===== Usual P2 parabolic elements =====
c
c     ----- Find local connectivities -----
          do ii=1,2
            inn = ijkl_e(ii)
            do no=1,node
              if(inn.eq.ijk_e(no)) then
                ijkl_n(ii) = no
              endif
            enddo
          enddo
          ijkl_n(3) = ijkl_n(1) +node/2
c         write(*,*) ijkl_n
c
c       ----- Set Weight  for Numerical Integration -----
          wg(1) = 5.d0/9.d0
          wg(2) = 8.d0/9.d0
          wg(3) = 5.d0/9.d0
          xg(1) = -dsqrt(3.d0/5.d0)
          xg(2) = 0.d0
          xg(3) =  dsqrt(3.d0/5.d0)
c
c       ----- Numerical Integration (line Int.) -----
          do 200 ig=1,3
            xl = xg(ig)
            wx = wg(ig)
c
c         === Shape Functions and Their Derivative ===
            sh(1) = (xl**2 - xl)/2.d0
            sh(2) = (xl**2 + xl)/2.d0
            sh(3) = 1.d0 - xl**2
            dn(1) = xl -0.5d0
            dn(2) = xl +0.5d0
            dn(3) = -2.d0*xl
c
c       === Compute the Length of THIS Edge ===
            gs = 0.d0  !gs(:) = 0.d0
            do in=1,3
              ia = ijkl_n(in)
              do ii=1,ndf
                gs(ii) = gs(ii) +dn(in)*xe(ii,ia)
              enddo
            enddo
            gsurf = dsqrt( gs(1)**2 +gs(2)**2 )*wx
c           write(*,*) gsurf
c
c         === Compute the Elemental Force Vector ===
            do in=1,3
              ia = ijkl_n(in)
              do nd=1,2
                iai = 2*(ia -1) +nd
                piai = vsl_e(nd+1)
                do jb=1,3
                  fe(iai) = fe(iai) +sh(in)*sh(jb)*piai*gsurf
                enddo
              enddo
            enddo
c
  200     CONTINUE
c
        endif
c
      endif
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
