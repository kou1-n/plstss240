      subroutine pcgsol(    nn,  NGSK,nsolvr, msol1,
     &                      ia,    ja,
     &                    stol,    aa,    bb,    xx,
     &                    temp,   tmp,
     &                  ierror )
c
      implicit none
c
      character MATDES(3)
      integer, parameter:: length = 128
      integer ipar(length)
      integer RCI_request
c
      double precision tnorm,fnorm,DNRM2
      double precision dpar(length)
c
      integer nn,NGSK,nsolvr,ierror,msol1
c     msol1 (for RCI CG Solver)
c           = 0: DONOT Use "Preconditioner"
c           = 1: USE "Jacobi Preconditioner"
      integer mm
c
      integer ia(nn+1)
      integer ja(nsolvr)
c
      double precision aa(NGSK)
      double precision bb(nn),xx(nn),temp(nn)
      double precision tmp(nn,4)
      double precision stol
c
      external DNRM2
c
c **********************************************************************
c     write(*,*) ' '
c     write(*,*) nn
c     do mm=1,nn
c       write(*,*) mm,ia(mm),aa(ia(mm))
c     enddo
c **********************************************************************
c
c ----------------------------------------------------------------------
c Purpose: Solve the Equation by Using RCI CG Solver (Intel MKL)
c                                            Ver. 9.1.023 2008/11/11
c ----------------------------------------------------------------------
c Parameters:
c  nn: Number of Equations
c  aa: Non-Zero values of the "SK" (values)
c  ia: Index of Row (rowIndex)
c  ja: Index of Column (column)
c  bb: Right Hand Side Vector/Matrix w/ dimension (nn,nrhs)
c    nrhs: Number of Right Hand Sides that Need to be Solved for
c          nrhs = 1: (Default, Single Right Hand Side; SRHS)
c          nrhs = 2--: (NOT USED!! Multiple Right Hand Side; MRHS)
c  xx: Array Contains Solution
c
c  ipar: Parameters for RCI CG Solver (integer)
c           ipar(128): for SRHS (Default)
c           ipar(128 +2*nrhs): for MRHS
c  dpar: Parameters for RCI CG Solver (double precision)
c           dpar(128): for SRHS (Default)
c           dpar(128 +2*nrhs): for MRHS
c  tmp: Temporary Space for the RCI CG routine
c           tmp(nn, 4): for SRHS
c           tmp(nn, 3+nrhs): for MRHS
c              tmp(:,1): current search direction
c              tmp(:,2): matrix multiplied by the current search dir.
c              tmp(:,3): current residual
c              tmp(:,4): The Inverse of the Preconditioner
c                        applied to the current residual
c                        tmp(:,4--nrhs) for MRHS
c
c  RCI_request : Information about the Result of the RCI CG routines
c         negative: the routine is completed w/ error or warnings
c                0: the task has been successfully completed
c                1: WE (users) should perform the following action;
c                   Mutiply the matrix by "tmp(1:n,1)", put the result
c                   in "tmp(1:n, 2)", and return the control to the
c                   "dcg/dcgmrhs" routine.
c                2: WE (users) should perform the following action;
c                   Perform the stopping test. If they fail, return
c                   the control to the "dcg/dcgmrhs" routine. Otherwise
c                   the solution if found and stored in the "x".
c                3: WE (users) should perform the following action;
c                   (for SRHS)
c                   apply the preconditioner to "tmp(1:n, 3)", put the
c                   results in "tmp(1:n,4)", and return the control to
c                   the "dcg" routine.
c                   (for MRHS)
c                   apply the preconditioner to "tmp(1:n, 3+ipar(3))",
c                   put the results in "tmp(1:n,3)", and return the
c                   control to the "dcgmrhs" routine.
c
c ierror: Error Indicator
c            0: No Error
c         ( Solver internal error )
c           21: RCI_request is NOT equal 0 in "dcg_init"
c           22: RCI_request is NOT equal 0 in "dcg_check"
c           23: RCI_request is NOT equal 0 in "dcg"
c           24: PCG Solver does NOT converged during the maximum
c                number of iterations; ( max(150, nn) )
c
c ----------------------------------------------------------------------
c Parameter "ipar(128)": various paramters for RCI CG Solver
c ----------------------------------------------------------------------
c   --- (1) Specify the Size of the Problem
c
c   --- (2) Specify the Type of Output for Error and Warning Messages
c             6: (Default) ALL messages are Displayed on the Screen
c
c   --- (3) for SRHS;
c             Contains the Current Stage of the RCI CG computations
c             Initial value = 1
c           for MRHS;
c             Contains the Right Hand Side for which the Calculations
c              are Currently Performed
c
c   --- (4) the Current Iteration Number
c             Initial value = 0
c
c   --- (6) Control the Error Message;
c             0: DONOT Output error message at all
c             else: Follow the Specification at ipar(2)
c             Default Value = 1
c
c   --- (7) Control the Warning Message;
c             0: DONOT Output warning message at all
c             else: Follow the Specification at ipar(2)
c             Default Value = 1
c
c   --- (8) Switch for "Stopping Test"
c             0: DONOT Perform the Stopping Test
c             else: DO the Stopping Test for the Maximum Number of
c                   Iterations; ipar(4).le.ipar(5)
c             Default Value = 1
c
c   --- (12:128) for SRHS          Reserved for Future Use
c   --- (12:128+ 2*nrhs) for MRHS  Reserved for Future Use
c
c ----------------------------------------------------------------------
c Parameter "dpar(128)": various paramters for RCI CG Solver
c ----------------------------------------------------------------------
c   --- (2) Absolute Tolerance
c             Default Value = 0.0E-00
c
c   --- (3) ( the square norm of initial residual )
c           ( computed in the "dcg/dcgmrhs" routines )
c             Default Value = 0.
c
c   --- (4) ( service variable; dpar(4) = dpar(1)*dpar(3) +dpar(2) )
c           ( computed in the "dcg/dcgmrhs" routines )
c             Default Value = 0.
c
c   --- (5) ( square norm of current residual)
c             Default Value = 0.
c
c   --- (6) ( square norm of residual from the previous iteration step )
c             Default Value = 0.
c
c   --- (7) ( "alpha" parameter of CG method )
c             Default Value = 0.
c
c   --- (8) ( "beta" parameter of CG method; dpar(8)=dpar(5)/dpar(6) )
c             Default Value = 0.
c
c   --- (9:128) for SRHS          Reserved for Future Use
c   --- (9:128+ 2*nrhs) for MRHS  Reserved for Future Use
c
c ****** Initialization the Solution ***********************************
      do mm=1,nn
        xx(mm) = 1.d0
      enddo
      MATDES(1)='D'
      MATDES(2)='L'
      MATDES(3)='N'
c
c ****** Initilize the Solver ******************************************
      CALL dcg_init( nn, xx, bb, RCI_request, ipar, dpar, tmp)
      if( RCI_request.ne.0) then
        WRITE(*,*) 'PCG Solver Returned Error Code 1', RCI_request
        ierror = 21
        RETURN
      endif
c
c ****** Set Parameters ************************************************
c     ipar(5) Specify the Maximum Number of Iterations
c               Default Value = min(150, nn)
      ipar(5) = max(150,nn)
c
c     ipar(9) Switch for "Residual Stopping Test"
c               0: DONOT Perform the Residual Stopping Test
c               else: Perform the Residual Stopping Test;
c                     dpar(5).le.dpar(4)=dpar(1)*dpar(3)+dpar(2)
c               Default Value = 0
c     ipar(9) = 1
c
c     ipar(10) Switch for "User Defined Stopping Test"
c               0: DONOT Perform the Use Defined Stopping Test
c               else: Perform the User Defined Stopping Test by
c                     setting "RCI_request = 2"
c               Default Value = 1
c     ipar(10) = 1
c
c     ipar(11) Switch for "Preconditioner"
c               0: DONOT use the "Preconditioner"
c               else: Preconditioned Version of the CG method by
c                     setting and asking "RCI_request = 3"
c               Default Value = 0 (Non-precondtioned version)
      ipar(11) = msol1
c
c     dpar(1) Relative Tolerance
c               Default Value = 1.0E-06
c     dpar(1) = 1.d-5
c
c
c ****** Check the Correctness & Consistency of the Parameters *********
      CALL dcg_check( nn, xx, bb, RCI_request, ipar, dpar, tmp)
      if( RCI_request.ne.0) then
        WRITE(*,*) 'PCG Solver Returned Error Code 2', RCI_request
        ierror = 22
        RETURN
      endif
c
c ****** Compute the Solution by RCI PCG Solver ************************
c        ( Reverse Communication Starts Here )
  100 CONTINUE
      CALL dcg( nn, xx, bb, RCI_request, ipar, dpar, tmp)
c
c     /* Solution was Found with the required precision */
c     =>> /* Solution has NOR Found with the precision */
c            ( CG Solver has been stopped by Iteration Number )
      if( RCI_request.eq.0 ) then
        WRITE(*,*) 'PCG Solver DOES NOT CONVERGED !!'
        ierror = 24
        RETURN
c
c     /* Compute the Vector AA*tmp(:,1), and puthe result in tmp(:,2) */
      elseif( RCI_request.eq.1 ) then
c       write(*,*) RCI_request, xx
        CALL MKL_DCSRSYMV('U', nn, aa, ia, ja, tmp(1,1), tmp(1,2))
        GOTO 100
c
c     /* Compute and Apply the Preconditioner Matrix C^{-1} on  */
c     /*  Vector tmp(:,3) and put the result in tmp(:,4)        */
      elseif( RCI_request.eq.3 ) then
        CALL MKL_DCSRSV('N', nn, 1.d0, MATDES,
     &             aa, ja, ia, ia(2), tmp(1,3), tmp(1,4))
        GOTO 100
c
c     /* DO the User-Defined Stopping Test */
      elseif( RCI_request.eq.2 ) then
        CALL MKL_DCSRSYMV('U', nn, aa, ia, ja, xx, temp)
        CALL DAXPY( nn, -1.d0, bb, 1, temp, 1)
        tnorm = DNRM2( nn, bb, 1)
        if(tnorm.lt.dsqrt(stol) ) tnorm = 1.d0
        fnorm = DNRM2( nn, temp, 1)
c       write(*,*) RCI_request, temp
c       write(*,*) 'norm ',tnorm,fnorm,fnorm/tnorm
c
c       --- The solution has been found
        if( fnorm/tnorm.lt.stol) then
c       if( fnorm.lt.stol) then
c         write(*,*) xx
c         GOTO 200
c
c       --- The solution has not been found yet
        else
          GOTO 100
        endif
c
c     /* The dcg routine Failed to Compute the Solution */
      else
        WRITE(*,*) 'PCG Solver Returned Error Code 3', RCI_request
        ierror = 23
        RETURN
      endif
c
c **********************************************************************
c ***** CG Solver Successfully Found the Solution **********************
c **********************************************************************
c
c ****** Update the Solution Vector "del" ******************************
      do mm=1,nn
        bb(mm) = xx(mm)
      enddo
c
c **********************************************************************
c     write(*,*) xx
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
