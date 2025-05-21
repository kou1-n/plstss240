      program PLSTss
c
c     *******************************************************
c     *                                                     *
c     *    PLSTss version 2.0   Aug. 2012                   *
c     *                                                     *
c     *    (c) copyright 2005   by K. Matsui & S. Suzuki    *
c     *                                                     *
c     *******************************************************
c
c   --------------------------------------------------------------------
c   Modification of this Program                               Kazumi M.
c   --------------------------------------------------------------------
c
c     Version0.1: First Version of elastic-PLaSTic analysis
c                  under Small Strain conditions.
c                 ( based on the "FEM_ss" )                     28/11/05
c
c     Version0.2: Change the Solver SKYLINE -->> PARDISO
c                 ( based on the "FEM_ss ver. 0.3" )            10/12/05
c
c     Version0.3: Modify the Routine assemble "SK"
c                                                               10/12/05
c
c     Version0.4: *** BUG FIXED !! in assembling process ***
c                 1. Add Controling section for analysis
c                    ( Explicit definition of loading parameter )
c                 2. OUTPUT some nodal/elemental data
c                 3. Tolerance of the computation is defined in 
c                    input CML-formatted file (/SOLUT/)
c                                                               01/01/06
c
c     Version0.5: Some brush-up on PARDISO solver
c                                                               05/05/06
c            0.51: Broken down by (S.Suzuki)
c                                                               10/01/06
c            0.52: Compute & Output Energy (S.Suzuki)
c                                                               10/01/06
c            0.53: Brushed-up by kzm
c                                                               18/01/06
c
c     Version0.6: ADDed loading/unloading conditions
c                 * Unloading case: use elastic tangent stiffness
c                                                               20/01/06
c      0.61-0.64: Try to add kinematic hardning, but failed. (S.Suzuki)
c                                                               01/02/06
c
c     Version0.7: Kinematic Hardening Law
c                 1. parameter for Isotropic/Kinematic combined
c                    hardening law is defined at "prop(15)" in /MATER/
c                    prop(15) = 0 --> pure kinematic
c                    prop(15) = 1 --> pure isotropic
c                                                               02/02/06
c           0.71: Extend to 3D and skyline solver.(S.Suzuki)
c           0.72: Fix the energy section
c           0.73: Fix the temperature section
c           0.74: Brushed up
c
c     Version1.0: STABLE VERSION
c                 * BUG-FIXED for dimension statement of "jdiag"
c                 * The others are equal to v.0.74
c                                                              01/Feb/08
c
c     Version1.1: Added the solver option
c                 * RCI CG Solver (Intel MKL)
c                   * isolvr = 0: Skyline Solver
c                   * isolvr = 1: PARDISO Solver (Intel MKL)
c                   * isolvr = 2: RCI CG Solver (Intel MKL)
c                                                              17/Nov/08
c
c     Version1.2: Some BUG-FIXEDs and Packed the subroutines
c      (=1.16)    * BUG-FIXED for kinematic hardening
c                 * BUG-FIXED in re-initializing res() after constr
c                 * Packed subroutine "stress" ( for both Q4/H8 )
c                                                              06/Dec/08
c
c     Version1.3: Some BUG-FIXEDs and Packed the subroutines
c      (=1.21)    * Added TRIA3 element
c                                                              06/Dec/08
c
c     Version1.4: Some BUG-FIXEDs and Packed the subroutines
c      (=1.31)    * In PARDISO solver, 
c                   "Analysis & Simbolic Factorization" process
c                     are executed at JUST first computation
c                                                              06/Dec/08
c
c     Version1.5: Some BUG-FIXEDs and Packed the subroutines
c                   for assembling "assemb" & "asmbpd"
c                                                              06/Jan/09
c
c     Version2.0: RE-build version
c      (=1.57)    * Major modifications on material routines
c                                                              26/Aug/12
c
c     Version2.3: Added some elements	by H.Sugiyama
c      (=2.29)    * TRIA6: 6-node parablic element
c                 * P1-ISO-P2/P0 element
c                                                              20/Nov/12
c
c   --------------------------------------------------------------------
c   --------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      real*4 sec
c
      character infem*30
c
c     dimension a(NDIM)
      real*8, allocatable::a(:)
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
      common /fnctn/ DELTA(3,3),EPSLN(3,3,3),FIT(3,3,3,3),DTENS(3,3,3,3)
c **********************************************************************
c
c     NDIM = 30000
      NDIM =100000000
      allocate( a(NDIM) )
c
      CALL CPU_Time(   sec)
      sec_00 = dble(sec)
c
c ***** Initializations ***********************************************
      a = 0.d0 ! a(:) = 0.d0
c
      ierror = 0
      iRCM = 0
c
c ***** Define Basic Tnsors *******************************************
c     ( Kronecker's Delta )
      DELTA = 0.d0 ! DELTA(:,:) = 0.d0
      do ii=1,3
        DELTA(ii,ii) = 1.d0
      enddo
c
c     ( Permutation Tensor )
      EPSLN(:,:,:) = 0.d0
      do ii=1,3
        jj = ii +1
        kk = ii +2
        if(jj.gt.3) jj = jj -3
        if(kk.gt.3) kk = kk -3
        EPSLN(ii,jj,kk) = 1.d0
        EPSLN(kk,jj,ii) = -1.d0
      enddo
c
c     ( Forth order identity tensor )
      do ll=1,3
        do kk=1,3
          do jj=1,3
            do ii=1,3
              FIT(ii,jj,kk,ll) = 0.5d0*( DELTA(ii,kk)*DELTA(jj,ll)
     &                                  +DELTA(ii,ll)*DELTA(jj,kk) )
            enddo
          enddo
        enddo
      enddo
c
c     ( Forth order Deviatoric operator )
      do ll=1,3
        do kk=1,3
          do jj=1,3
            do ii=1,3
              DTENS(ii,jj,kk,ll) = DELTA(ii,kk)*DELTA(jj,ll)
     &                            -DELTA(ii,jj)*DELTA(kk,ll)/3.d0
            enddo
          enddo
        enddo
      enddo
c
c ***** Open some files ***********************************************
      CALL flopen(ierror,  iRCM)
      if(ierror.ne.0) GOTO 9999
c
c ***** Check size of each array **************************************
      CALL chkary(     a,
     &              NDIM,  NGSK, NGSKo,   neq,  neqo,
     &              iRCM,
     &                m1,    m2,    m3,    m4,    m5,
     &                m6,    m7,    m8,    m9,   m10,
     &               m11,   m12,   m13,   m14,   m15,
     &               m16,   m17,   m18,   m19,   m20,
     &               m21,   m22,   m23,   m24,   m25,
     &               m26,   m27,   m28,   m29,
c
     &                n1,    n2,    n3,    n4,    n5,
     &                n6,    n7,    n8,    n9,   n10,
     &               n11,   n12,   n13,   n14,   n15,
     &               n16,   n17,   n18,   n19,   n20,
     &               n21,   n22,   n23,   n24,   n25,
     &               n26,   n27,   n28,   n29,   n30,
     &               n31,   n32,   n33,   n34,   n35,
     &               n36,   n37,   n38,   n39,   n40,
     &               n41,   n42,   n43,   n44,   n45,
     &            ierror )
      if(ierror.ne.0) GOTO 9999
c
c ***** Re-Initialize the whole Array *********************************
      a = 0.d0 ! a(:) = 0.d0
c
c ***** Main routine for Finite Element Analysis **********************
      CALL analys(  NGSK, NGSKo,   neq,  neqo,  iRCM,
     &            a( m1),a( m2),a( m3),a( m4),a( m5),
     &            a( m6),a( m7),a( m8),a( m9),a(m10),
     &            a(m11),a(m12),a(m13),a(m14),a(m15),
     &            a(m16),a(m17),a(m18),a(m19),a(m20),
     &            a(m21),a(m22),a(m23),a(m24),a(m25),
     &            a(m26),a(m27),a(m28),a(m29),
c
     &            a( n1),a( n2),a( n3),a( n4),a( n5),
     &            a( n6),a( n7),a( n8),a( n9),a(n10),
     &            a(n11),a(n12),a(n13),a(n14),a(n15),
     &            a(n16),a(n17),a(n18),a(n19),a(n20),
     &            a(n21),a(n22),a(n23),a(n24),a(n25),
     &            a(n26),a(n27),a(n28),a(n29),a(n30),
     &            a(n31),!a(n32),a(n33),a(n34),a(n35),
c    &            a(n36),a(n37),a(n38),a(n39),a(n40),
c    &            a(n41),a(n42),a(n43),a(n44),a(n45),
     &            secslv,
     &            ierror )
      if(ierror.ne.0) GOTO 9999
c
c ***** Post Process ***************************************************
      WRITE(lwa,'(a7)') '/ENDOF/'
      CLOSE(lwa)
      CLOSE(lwb)
      CLOSE(lwc)
      CLOSE(lwd)
      CLOSE(lwe)
      CLOSE(lwf)
c
      CALL CPU_Time(   sec)
      sec_99 = dble(sec)
      WRITE(*,100) sec_99,secslv
c
      deallocate( a )
c
      STOP 'Program "PLSTss" has done !!'
c
c ***** Error Section **************************************************
 9999 CONTINUE
c
      WRITE(*,9100) ierror
      deallocate( a )
      STOP 'Emergency STOP!'
c
c **********************************************************************
  100 FORMAT(/,5x,'The Computation Time:',e15.7,' sec',
     &       /,5x,'(for solving linear equation',e15.7,' sec'//)
 9100 FORMAT(/,'Some error has occured Check the error ID',
     &       /,2x,'error ID : ',i5)
c
c **********************************************************************
c **********************************************************************
      END
c
c *********************************************************************
c **********************   Block Data Statement    ********************
c *********************************************************************
c
      BLOCK DATA
c
      common /iodev/ lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
c     common /tvalu/ ctol,stol
c
c     double precision ctol
      integer lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf
c
c     ( I/O device )
      data lra,lrb,lwa,lwb,lwc,lwd,lwe,lwf/10,11,20,21,22,23,24,25/
c
c     ( Tolerance of this computation )
c   ==>> Defined in /SOLUT/ section in CML-formatted file
c     data ctol / 1.d-10 /
c
      END
