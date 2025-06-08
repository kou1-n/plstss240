      subroutine pars00(    nn,  NGSK,nsolvr, msol1,
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
c ----- For PARDISO Solver -----
      integer omp_get_max_threads
c **********************************************************************
c
c ----------------------------------------------------------------------
c Purpose: Pre-Condition for PARDISO and Solve
c ----------------------------------------------------------------------
c Parameters:
c  pt: address pointer for solver internal variables
c  maxfct: max. number of factors which are kept at the same time
c          typical value = 1
c  mnum: actual matrix for the solution phase.
c          typical value = 1
c  mtype: Definition of the Matrix Type
c            2: Real and Symmetric POSITIVE DEFINE Matrix
c           -2: Real and Symmetric INDEFINITE Matrix
c           11: Real and Unsymmetric Matrix
c                               (see, Reference Manual for other values)
c phase: Cotrol Parameter for the Solver
c           11: Analysis, Symbolic Fact.
c           12: Analysis, Symbolic Fact., Numerical Fact.
c           13: Analysis, Symbolic Fact., Numerical Fact., Solve
c           22: Numerical Factorization
c           23: Numerical Factorization, Solve
c           33: Solve
c            0: Release Internal Memory for "L" and "U"
c           -1: Release ALL Internal Memory
c nn: Number of Equations
c aa: Non-Zero values of the "SK" (values)
c ia: Index of Row (rowIndex)
c ja: Index of Column (column)
c perm: Permutation Vector of Size nn
c nrhs: Number of Right-Hand Sides that Need to be Solved for
c iparm: Various Parameters for PARDISO routine
c        iparm(1) = 0 then, default parameters are set
c msglvl: Message Level Information
c            0: No output (printout) from PARDISO
c            1: Prints Statistical Information in file
c                       "pardiso.stat.(#nproc)"
c bb: Right Hand Side Vector/Matrix w/t dimension(n,nrhs)
c      On output, the Array is Replaced with the Solution if iparm(6)=1
c xx: On output, the Array Contains Solution if iparm(6)=0
c error: Error Indicator
c            0: No Error
c           -1: Input Inconsistent                             >> 21
c           -2: Not Enough Memory                              >> 22
c           -3: Reordering Problem                             >> 23
c           -4: ZERO PIVOT, Numerical Factorization Problem    >> 24
c           -5: Unclassified (Internal) Error                  >> 25
c           -6: Preordering Failed (matrix types 11, 13 only)  >> 26
c           -7: Diagonal Matrix Problem                        >> 27
c                               (ierror = 21, 22, 23, 24, 25, 26, 27)
c
c ****** Set Parameters ************************************************
      maxfct = 1
      mnum = 1
      mtype = 2
      phase = 13
      nrhs = 1
      msglvl = msol1
c
c ----------------------------------------------------------------------
c Parameter "iparm(64)": various paramters for PARDISO routine
c                        based on 630813-033US/October 2009
c ----------------------------------------------------------------------
c   --- (1) if = 0, then default values are used for (2) and (4)--(64)
c           * NOT solver default =1
c     iparm( 1) = 1
      iparm( 1) = 0
c   --- (2) Control the Fill-in reducing ordering for the input matrix
c             0: the minimum degree algorithm
c             2: the nested dissection algorithm from the METIS package
c                (default)
      iparm( 2) = 2
c
c
c   --- (3) Number of Processors (No default value!!)
c           * The number must be equal to the OpenMP environment
c            variable "OMP_NUM_THREADS"
c     iparm( 3) = 2
c     iparm( 3) = omp_get_max_threads()
c
c
c   --- (4) Control parameter for preconditioned CGS for unsymmetric or
c          structural symmetric matrices and CG for symmetric matrices
c             0: Always computed as required by "phase"
c             1: or 2: (see ref. manual for detail)
      iparm( 4) = 0
c
c   --- (5) Use my own fill reducting permutation vector
c             0: no (default) / 1: yes
      iparm( 5) = 0
c   --- (6) Solution of the equations
c             0: The value of "xx" and "bb" is not changed
      iparm( 6) = 0
c   --- (7) Not in use
      iparm( 7) = 0
c   --- (8) Maximum number of iterative refinement steps
c             default = 0
      iparm( 8) = 0
c   --- (9) Not in use
      iparm( 9) = 0
c
c
c   --- (10) Perturbe the pivot elements with 1.0E-13
c             define the power of 10^?? (default = 13)
      iparm(10) = 13
c
c
c   --- (11) Parameter for "UNsymmetric matrices" (default = 1)
      iparm(11) = 1
c   --- (12) Not in use
      iparm(12) = 0
c   --- (13) Not in use
      iparm(13) = 0
c   --- (14) [OUTPUT]: Purturbed pivots (for unsymmetric families)
      iparm(14) = 0
c
c
c   --- (15) [OUTPUT]: Total Peak Memory in KBytes
      iparm(15) = 0
c   --- (16) [OUTPUT]: Permanent Memory in KBytes
      iparm(16) = 0
c   --- (17) [OUTPUT]: Total Double Precision Memory consumption(KBytes)
      iparm(17) = 0
c
c
c   --- (18) [OUTPUT]: Number of nonzeros on the factors
      iparm(18) = 0
c   --- (19) [OUTPUT]: Number of Operations in MFlops
      iparm(19) = 0
c   --- (20) [OUTPUT]: CG/CGS diagnostics
      iparm(20) = 0
c ----------------------------------------------------------------------
c
c
c ****** Initialization ************************************************
      do mm=1,64
        pt(mm) = 0
      enddo
c
c ****** Solve Linear Equations by "PARDISO" ***************************
      CALL pardiso(    pt,maxfct,   mnum, mtype, phase,
     &                 nn,    aa,     ia,    ja,  perm,
     &               nrhs, iparm, msglvl,    bb,    xx,
     &             ierror )
      if(ierror.ne.0) then
        ierror = 20 -ierror
        RETURN
      endif
c
c ****** Update the Solution Vector "del" ******************************
      bb = xx ! bb(:) = xx(:)
c     do mm=1,nn
c       bb(mm) = xx(mm)
c     enddo
c
c ****** Output from PARDISO Solver ************************************
      WRITE(*,9001) iparm(15), iparm(16),iparm(17)
c
 9001 FORMAT(/,10x,'Total Peak Memory during factorization:   ',
     &              i8,' (KB)',/,
     &         10x,'Parmanet Memory from phase1 to phase3:    ',
     *              i8,' (KB)',/,
     &         10x,'Total Double Precision Memory for solving:',
     &              i8,' (KB)',/ )
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
