      subroutine pars99(    nn,  NGSK,nsolvr, msol1,
     &                      ia,    ja,  perm,
     &                      aa,    bb,    xx,
     &                      pt, iparm,
     &                  maxfct,  mnum, mtype, phase,  nrhs,
     &                  ierror )
c
      implicit none
c
      integer*8 pt(64)
c
      integer nn,mm,ierror,NGSK,nsolvr,idum,msol1
      integer maxfct,mnum,mtype,phase,nrhs,msglvl
c
      integer ia(nn+1)
      integer ja(nsolvr),perm(nn)
c
      integer iparm(64)
c
      double precision dum
c
      double precision aa(NGSK)
      double precision bb(nn),xx(nn)
c
c **********************************************************************
CC
CC    --> Already Set in the subroutine "pars00"
CC
c ****** Termination and Release of Memory *****************************
      phase = -1
      CALL pardiso(    pt,maxfct,   mnum, mtype, phase,
     &                 nn,   dum,   idum,  idum,  idum,
     &               nrhs, iparm, msglvl,   dum,   dum,
     &             ierror )
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
