      subroutine addres(  nend,  NGSK, NGSKo,   neq,  neqo,
     &                      m1,    m2,    m3,    m4,    m5,
     &                      m6,    m7,    m8,    m9,   m10,
     &                     m11,   m12,   m13,   m14,   m15,
     &                     m16,   m17,   m18,   m19,   m20,
     &                     m21,   m22,   m23,   m24,   m25,
     &                     m26,   m27,   m28,   m29,
c
     &                      n1,    n2,    n3,    n4,    n5,
     &                      n6,    n7,    n8,    n9,   n10,
     &                     n11,   n12,   n13,   n14,   n15,
     &                     n16,   n17,   n18,   n19,   n20,
     &                     n21,   n22,   n23,   n24,   n25,
     &                     n26,   n27,   n28,   n29,   n30,
     &                     n31,   n32,   n33,   n34,   n35,
     &                     n36,   n37,   n38,   n39,   n40,
     &                     n41,   n42,   n43,   n44,   n45,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
      common /bound/ lnum,locc,nspc,mpc,npoin,npres,nbody,ntn
      common /solvr/ isolvr,nsolvr,msol1,nmax,nodsk,neqsol,jsol
      common /cntrl/ nstep,istart,iarc,itrobj,incomp
      common /print/ lpstp,luprt,lfprt,lnprt,lsprt,lbprt
c
c **********************************************************************
      nxh    = nx/2 +mod(nx,2)
      nelxh  = nelx/2 +mod(nelx,2)
      nspch  = nspc/2 +mod(nspc,2)
      npoinh = npoin/2 +mod(npoin,2)
      npresh = npres/2 +mod(npres,2)
      nbodyh = nbody/2 +mod(nbody,2)
      ntnh   = ntn/2 +mod(ntn,2)
      neqo1h = (neqo +1)/2 +mod(neqo+1,2)
      ngskoh = NGSKo/2 +mod(NGSKo,2)
      neq1h  = (neq+1)/2 +mod(neq+1,2)
      nsolvh = nsolvr/2 +mod(nsolvr,2)
      nsolv1h= (nsolvr+1)/2 +mod(nsolvr+1,2)
      nodskh = nodsk/2 +mod(nodsk,2)
      lpstph = lpstp/2 +mod(lpstp,2)
      lmath = lmat/2 +mod(lmat,2)
c
      neqh   = neq/2 +mod(neq,2)
c
c ****** Integer Type Arrays *******************************************
      m1  = 1
c     a( m1):  m1 --  m2 -1: ijk(node,nelx)
c     a( m2):  m2 --  m3 -1: mpe(nelx)
c     a( m3):  m3 --  m4 -1: mang(nelx)
c     a( m4):  m4 --  m5 -1: nfix(nspc)
c     a( m5):  m5 --  m6 -1: ndir(6,nspc)
      m2  = m1  +node*nelxh
      m3  = m2  +nelxh
      m4  = m3  +nelxh
      m5  = m4  +nspch
                             if(nspch.eq.0) m5 = m4 +1
      m6  = m5  +6*nspch
                             if(nspch.eq.0) m6 = m5 +1
c
c     a( m6):  m6 --  m7 -1: npload(npoin)
c     a( m7):  m7 --  m8 -1: nsload(npres)
c     a( m8):  m8 --  m9 -1: nbload(nbody)
c     a( m9):  m9 -- m10 -1: ijkl(nsn,npres)
c     a(m10): m10 -- m11 -1: ncp(neq)
      m7  = m6  +npoinh
                             if(npoinh.eq.0) m7 = m6 +1
      m8  = m7  +npresh
                             if(npresh.eq.0) m8 = m7 +1
      m9  = m8  +nbodyh
                             if(nbodyh.eq.0) m9 = m8 +1
      m10 = m9  +nsn*npresh
                             if(npresh.eq.0) m10 = m9 +1
      m11 = m10 +neqh
c
c     a(m11): m11 -- m12 -1: mdof(neq)
c     a(m12): m12 -- m13 -1: idof(neq)
c     a(m13): m13 -- m14 -1: jdiag(neq+1)
c     a(m14): m14 -- m15 -1: iw_neq(neq)
c     a(m15): m15 -- m16 -1: index0(neqo+1)
      m12 = m11 +neqh
      m13 = m12 +neqh
      m14 = m13 +neq1h
      m15 = m14 +neqh
      m16 = m15 +neqo1h
c
c     a(m16): m16 -- m17 -1: index1(NGSKo)
c     a(m17): m17 -- m18 -1: jcolmn(nsolvr)     : for PARDISO
c     a(m18): m18 -- m19 -1: list(nmax,nx)      : for PARDISO
c     a(m19): m19 -- m20 -1: ncolmn(nodsk)      : for PARDISO
c     a(m20): m20 -- m21 -1: mpstp(lpstp)
      m17 = m16 +ngskoh
      m18 = m17 +nsolv1h
      m19 = m18 +nmax*nxh
      m20 = m19 +nodskh
      m21 = m20 +lpstph
c
c     a(m21): m21 -- m22 -1: muprt(luprt)
c     a(m22): m22 -- m23 -1: mfprt(lfprt)
c     a(m23): m23 -- m24 -1: mnprt(lnprt)
c     a(m24): m24 -- m25 -1: msprt(lsprt)
c     a(m25): m25 -- m26 -1: mbprt(lbprt)
      m22 = m21 +luprt
      m23 = m22 +lfprt
      m24 = m23 +lnprt
      m25 = m24 +lsprt
      m26 = m25 +lbprt
c
c     a(m26): m26 -- m27 -1: idep0(ngaus,nelx)
c     a(m27): m27 -- m28 -1: idep1(ngaus,nelx)
c     a(m28): m28 -- m29 -1: matid(lmat)
      m27 = m26 +ngaus*nelxh
      m28 = m27 +ngaus*nelxh
      m29 = m28 +nelxh
      m30 = m29 +lmath

      mend = m30 -1
c
c ****** Real Type Arrays **********************************************
      n1  = mend +1
c     a( n1):  n1 --  n2 -1: xyz(3,nelx)
c     a( n2):  n2 --  n3 -1: prop(20,lmat)
c     a( n3):  n3 --  n4 -1: angle(3,lang)
c     a( n4):  n4 --  n5 -1: vfix(6,nspc)
c     a( n5):  n5 --  n6 -1: vpload(6,npoin)
      n2  = n1  +3*nx
      n3  = n2  +20*lmat
      n4  = n3  +3*lang
                             if(lang.eq.0) n4 = n3 +1
      n5  = n4  +6*nspc
                             if(nspc.eq.0) n5 = n4 +1
      n6  = n5  +6*npoin
                             if(npoin.eq.0) n6 = n5 +1
c
c     a( n6):  n6 --  n7 -1: vsload(4,npres)
c     a( n7):  n7 --  n8 -1: vbload(3,nbody)
c     a( n8):  n8 --  n9 -1: foc(neq)
c     a( n9):  n9 -- n10 -1: disp(neq)
c     a(n10): n10 -- n11 -1: del(neq)
      n7  = n6  +4*npres
                             if(npres.eq.0) n7 = n6 +1
      n8  = n7  +6*nbody
                             if(nbody.eq.0) n8 = n7 +1
      n9  = n8  +neq
      n10 = n9  +neq
      n11 = n10 +neq
c
c     a(n11): n11 -- n12 -1: sigma(6,nelx)
c     a(n12): n12 -- n13 -1: epsln(6,nelx)
c     a(n13): n13 -- n14 -1: von(nelx)
c     a(n14): n14 -- n15 -1: dw_sol(neqsol)
c     a(n15): n15 -- n16 -1: du(neq)
      n12 = n11 +6*nelx
      n13 = n12 +6*nelx
      n14 = n13 +nelx
      n15 = n14 +neqsol
      n16 = n15 +neq
c
c     a(n16): n16 -- n17 -1: u0(neq)
c     a(n17): n17 -- n18 -1: u1(neq)
c     a(n18): n18 -- n19 -1: res(neq)
c     a(n19): n19 -- n20 -1: fint(neq)
c     a(n20): n20 -- n21 -1: eps(nelx)
      n17 = n16 +neq
      n18 = n17 +neq
      n19 = n18 +neq
      n20 = n19 +neq
      n21 = n20 +nelx
c
c     a(n21): n21 -- n22 -1: finc(nstep)
c     a(n22): n22 -- n23 -1: pene(nelx)
c     a(n23): n23 -- n24 -1: eene(nelx)
c     a(n24): n24 -- n25 -1: temp(nelx)
c     a(n25): n25 -- n26 -1: tempd(nelx)
      n22 = n21 +nstep
      n23 = n22 +nelx
      n24 = n23 +nelx
      n25 = n24 +nelx
      n26 = n25 +nelx
c
c     a(n26): n26 -- n27 -1: dndx_g(ndf,node,ngaus,nelx)
c     a(n27): n27 -- n28 -1: det_g(ngaus,nelx)
c     a(n28): n28 -- n29 -1: ctensg(3,3,3,3,ngaus,nelx)
c     a(n29): n29 -- n30 -1: dhist0(20,ngaus,nelx)
c     a(n30): n30 -- n31 -1: dhist1(20,ngaus,nelx)
      n27 = n26 +ndf*node*ngaus*nelx
      n28 = n27 +ngaus*nelx
      n29 = n28 +3*3*3*3*ngaus*nelx
      n30 = n29 +20*ngaus*nelx
      n31 = n30 +20*ngaus*nelx
c
c     a(n31): n31 -- n32 -1: sk(NGSK)
c     a(n32): n32 -- n33 -1: 
c     a(n33): n33 -- n34 -1: 
c     a(n34): n34 -- n35 -1: 
c     a(n35): n35 -- n36 -1: 
      n32 = n31 +NGSK
      n33 = n32 +0
      n34 = n33 +0
      n35 = n34 +0
      n36 = n35 +0
c
c     a(n36): n36 -- n37 -1: 
c     a(n37): n37 -- n38 -1: 
c     a(n38): n38 -- n39 -1: 
c     a(n39): n39 -- n40 -1: 
c     a(n40): n40 -- n41 -1: 
      n37 = n36 +0
      n38 = n37 +0
      n39 = n38 +0
      n40 = n39 +0
      n41 = n40 +0
c
c     a(n41): n39 -- n41 -1: 
c     a(n42): n41 -- n42 -1: 
c     a(n43): n42 -- n43 -1: 
c     a(n44): n43 -- n44 -1: 
c     a(n45): n44 -- n45 -1: 
      n42 = n41 +0
      n43 = n42 +0
      n44 = n43 +0
      n45 = n43 +0
      n46 = n43 +0
c
      nend = n32 -1
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
