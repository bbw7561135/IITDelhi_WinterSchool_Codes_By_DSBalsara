! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

      SUBROUTINE EIGENVECMHD

     1             ( nriembegin, nriemend, iprimitive,

     1               gamma, smallnum, prsfloor, rhofloor, pi, 
     1               debug_stop,

     1               gasconst,
     1               cp_specie, molewt_specie, heat_of_formation,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Returns normalized left and right eigenvectors for a specified set
! of mean states for adiabatic MHD. 

! When "iprimitive == 1" the eigenvectors expressed in terms of primitive
! variables are returned.

! When "iprimitive == 0" the eigenvectors in terms of conserved variables
! are returned.

! When "iprimitive == -1" the eigenvectors in terms of "conserved" variables
! are returned. Here we use the equation for entropy advection instead of
! the equation for total energy. This formulation is useful for maintaining
! pressure positivity.

! When primitive variables are used the vector of variables is in the
! order given by ( rho, vlx,y,z, prs, bfy,z, species_fractions)^T.

! When conserved variables are used the vector of variables is in the
! order given by ( rho, momx,y,z, eng, bfy,z, species_concentration)^T.

! When entropy-based "conserved" variables are used the vector of variables
! is in the order given by :
! ( rho, momx,y,z, prs/rho**(gamma-1), bfy,z, species_concentration)^T.

! Properly normalized eigenvectors as given in Roe and Balsara (1996)
! are generated.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL RHO_PRS_IN_ENG_T_GAMA_DERVS_OUT,
     1         RHO_PRS_IN_ENG_T_GAMA_OUT

! ----------------------------------------------------------------------

      INTEGER, INTENT ( IN) :: nriembegin, nriemend, iprimitive

      INTEGER, INTENT ( INOUT) :: debug_stop

      REAL, INTENT ( IN) :: gamma, smallnum, prsfloor, rhofloor,
     1      pi, gasconst

      REAL, DIMENSION ( 0 : NFLUID), INTENT ( IN) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1     INTENT ( INOUT) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     bfxbar, bfybar, bfzbar, xbar

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1     INTENT ( INOUT) ::
     1      speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID), INTENT ( INOUT) ::
     1     eigenval

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID), INTENT ( INOUT) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

      INTEGER i, ii, iii, j, k, ifluid

      REAL fourpi, small_speed


! Variables for making wave speeds and eigenvectors.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     va_x, c_s, v_a, v_f, v_s, beta_y, beta_z,
     1     alpha_f, alpha_s, sgnbx, sq_4pi_rho, va_perp,
     1     eps_zone, one_by_4pi_rho,

     1     eng, temperature,
     1     mean_mol_wt, gamma_eos, gamma_soundspeed,
     1     dprs_drho, dprs_deng, deng_drho, deng_dprs,

     1     temparr_01, temparr_02, temparr_03

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     dprs_dspecie, deng_dspecie

! Temporary usage of left and right eigenvectors in the primitive variables.
! Transformation matrices to go from conserved to primitive and back.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID) ::

     1     lf_eigenvec1, rt_eigenvec1, dudw, dwdu

! ----------------------------------------------------------------------

! Make some settings.

      fourpi = 4.0 * pi

! If the alfven speeds in the x or ( y, z) directions fall below "small_speed"
! times the local sound speed, manoeuvers have to be done to prevent the
! slow and fast eigenvectors from developing pathologies.

      small_speed = 0.000001

! ----------------------------------------------------------------------

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)"EIGENVECMHD : NFIELD .NE. 7", NFIELD, 7
      debug_stop = 1
      STOP
      END IF

      IF ( ( nriembegin .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriembegin .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EIGENVECMHD : nriembegin out of bounds",
     1       nriembegin, nriembegin
      debug_stop = 1
      STOP
      END IF

      IF ( ( nriemend .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriemend .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EIGENVECMHD : nriemend out of bounds",
     1       nriemend, nriemend
      debug_stop = 1
      STOP
      END IF

! ----------------------------------------------------------------------

* Bound the density and pressure variables. Divide "xbar" by 4 * pi.

      DO i = nriembegin, nriemend

      rhobar ( i) = AMAX1 ( rhobar ( i), rhofloor)
      prsbar ( i) = AMAX1 ( prsbar ( i), prsfloor)

      xbar ( i) = xbar ( i) / fourpi

      END DO

! ----------------------------------------------------------------------

! Fix up "gamma_soundspeed" and partial derivatives of various
! thermodynamic quantities wrt one another.

! The below subroutine also automatically checks to make sure that the
! species fractions are individually and cumulatively within reasonable
! ranges. 

      CALL RHO_PRS_IN_ENG_T_GAMA_DERVS_OUT
     1        ( nriembegin, nriemend, 0,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rhobar, prsbar, speciebar,

     1          eng, temperature,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed,

     1          dprs_drho, dprs_deng, dprs_dspecie,
     1          deng_drho, deng_dprs, deng_dspecie)

! ----------------------------------------------------------------------

! Make eigenvalues:

      DO i = nriembegin, nriemend

! ----------

      sgnbx ( i) = SIGN ( 1.0, bfxbar ( i) )

      sq_4pi_rho ( i) = SQRT ( fourpi * rhobar ( i) )

! ----------

      one_by_4pi_rho ( i) = 1.0 / ( fourpi * rhobar ( i) )

      va_x ( i) = bfxbar ( i)**2 * one_by_4pi_rho ( i)
      va_x ( i) = SQRT ( va_x ( i) )

      c_s ( i) = gamma_soundspeed ( i) * prsbar ( i) / rhobar ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      v_a ( i) = ( bfxbar ( i)**2 + bfybar ( i)**2
     1           + bfzbar ( i)**2) * one_by_4pi_rho ( i)
      v_a ( i) = SQRT ( v_a ( i) )

      temparr_01 ( i) = v_a ( i)**2 + c_s ( i)**2
      temparr_02 ( i) = temparr_01 ( i)**2
     1                - 4.0 * va_x ( i)**2 * c_s ( i)**2
      temparr_02 ( i) = AMAX1 ( temparr_02 ( i), 0.0)
      temparr_02 ( i) = SQRT ( temparr_02 ( i))

      v_f ( i) = 0.5 * ( temparr_01 ( i) + temparr_02 ( i))
      v_f ( i) = SQRT ( v_f ( i) )

      v_s ( i) = 0.5 * ( temparr_01 ( i) - temparr_02 ( i))
      v_s ( i) = AMAX1 ( v_s ( i), 0.0)
      v_s ( i) = SQRT ( v_s ( i) )

! ----------

! Now make quantities that are essential for making the eigenvectors ( though
! inessential for just the eigenvalues).

      eps_zone ( i) = small_speed * c_s ( i)


      va_perp ( i) = SQRT ( bfybar ( i)**2 + bfzbar ( i)**2)
     1             / sq_4pi_rho ( i)


      IF ( va_perp ( i) .GT. eps_zone ( i)) THEN

      beta_y ( i) = bfybar ( i)
     1            / ( va_perp ( i) * sq_4pi_rho ( i) )
      beta_z ( i) = bfzbar ( i)
     1            / ( va_perp ( i) * sq_4pi_rho ( i) )

      ELSE

      beta_y ( i) = 1.0 / SQRT ( 2.0)
      beta_z ( i) = beta_y ( i)

      END IF

! ----------

      eigenval ( i, 1) = vlxbar ( i) - v_f ( i)
      eigenval ( i, 2) = vlxbar ( i) - va_x ( i)
      eigenval ( i, 3) = vlxbar ( i) - v_s ( i)
      eigenval ( i, 4) = vlxbar ( i)
      eigenval ( i, 5) = vlxbar ( i) + v_s ( i)
      eigenval ( i, 6) = vlxbar ( i) + va_x ( i)
      eigenval ( i, 7) = vlxbar ( i) + v_f ( i)

! ----------

      END DO

! ----------

#if ( NFLUID >= 1)
      DO i = nriembegin, nriemend

      eigenval ( i, NFIELD + 1: NFIELD + NFLUID) = vlxbar ( i)

      END DO
#endif

! ----------------------------------------------------------------------

! Make a consistent set of "alpha_f/s" for all the cases.

      DO i = nriembegin, nriemend

! ----------

      IF ( v_f ( i) .GE. v_s ( i) + eps_zone ( i) ) THEN

! Straightforward evaluation is adequate here. We are far away from
! the triple umbilic point.

      temparr_01 ( i) = 1.0 / ( v_f ( i)**2 - v_s ( i)**2)

      alpha_f ( i) = AMAX1 ( ( c_s ( i)**2 - v_s ( i)**2), 0.0 )
      alpha_f ( i) = SQRT ( temparr_01 ( i) * alpha_f ( i))

      alpha_s ( i) = AMAX1 ( ( v_f ( i)**2 - c_s ( i)**2), 0.0 )
      alpha_s ( i) = SQRT ( temparr_01 ( i) * alpha_s ( i))

      ELSE IF ( ABS ( va_x ( i) - c_s ( i)) .GT.
     1          smallnum * c_s ( i) ) THEN

! We are very close to the triple umbilic point but not on it.

      temparr_01 ( i) = va_perp ( i)
     1       / ( va_x ( i) - c_s ( i))
      temparr_01 ( i) = ATAN ( temparr_01 ( i)) / 2.0

      alpha_f ( i) = ABS ( SIN ( temparr_01 ( i)))
      alpha_s ( i) = ABS ( COS ( temparr_01 ( i)))

      ELSE

! We are bang on the triple umbilic point ( indistinguishable
! upto machine accuracy). Now which way we resolve this could
! only depend on the structure of the waves in the adjoining zones
! and we don't address that issue here.

      temparr_01 ( i) = 0.25 * pi
     1                * SIGN ( 1.0, va_x ( i) - c_s ( i))

      alpha_f ( i) = ABS ( SIN ( temparr_01 ( i)))
      alpha_s ( i) = ABS ( COS ( temparr_01 ( i)))

      END IF

! ----------

      END DO

! ----------------------------------------------------------------------

! First zero out the eigenvectors. This is only needed because we can
! then fill in the species eigenvectors by putting 1's at suitable
! places on the diagonals.

      rt_eigenvec1 ( :, :, :) = 0.0
      lf_eigenvec1 ( :, :, :) = 0.0

! ----------

! Make right & left eigenvectors for alfven and entropy waves.

      DO i = nriembegin, nriemend

! ----------

! Alfven waves

! ----------


!      rt_eigenvec1 ( i, 1, 2) = 0.0
!      rt_eigenvec1 ( i, 1, 6) = 0.0

!      rt_eigenvec1 ( i, 2, 2) = 0.0
!      rt_eigenvec1 ( i, 2, 6) = 0.0

      rt_eigenvec1 ( i, 3, 2) = - beta_z ( i)
      rt_eigenvec1 ( i, 3, 6) = - rt_eigenvec1 ( i, 3, 2)

      rt_eigenvec1 ( i, 4, 2) = beta_y ( i)
      rt_eigenvec1 ( i, 4, 6) = - rt_eigenvec1 ( i, 4, 2)

!      rt_eigenvec1 ( i, 5, 2) = 0.0
!      rt_eigenvec1 ( i, 5, 6) = 0.0

      rt_eigenvec1 ( i, 6, 2) = - sgnbx ( i)
     1                        * sq_4pi_rho ( i) * beta_z ( i)
      rt_eigenvec1 ( i, 6, 6) = rt_eigenvec1 ( i, 6, 2)

      rt_eigenvec1 ( i, 7, 2) = sgnbx ( i)
     1                        * sq_4pi_rho ( i) * beta_y ( i)
      rt_eigenvec1 ( i, 7, 6) = rt_eigenvec1 ( i, 7, 2)


! ----------


!      lf_eigenvec1 ( i, 2, 1) = 0.0
!      lf_eigenvec1 ( i, 6, 1) = 0.0

!      lf_eigenvec1 ( i, 2, 2) = 0.0
!      lf_eigenvec1 ( i, 6, 2) = 0.0

      lf_eigenvec1 ( i, 2, 3) = - 0.5 * beta_z ( i)
      lf_eigenvec1 ( i, 6, 3) = - lf_eigenvec1 ( i, 2, 3)

      lf_eigenvec1 ( i, 2, 4) = 0.5 * beta_y ( i)
      lf_eigenvec1 ( i, 6, 4) = - lf_eigenvec1 ( i, 2, 4)

!      lf_eigenvec1 ( i, 2, 5) = 0.0
!      lf_eigenvec1 ( i, 6, 5) = 0.0

      lf_eigenvec1 ( i, 2, 6) = - 0.5 * sgnbx ( i)
     1                        * beta_z ( i) / sq_4pi_rho ( i)
      lf_eigenvec1 ( i, 6, 6) = lf_eigenvec1 ( i, 2, 6)

      lf_eigenvec1 ( i, 2, 7) = 0.5 * sgnbx ( i)
     1                        * beta_y ( i) / sq_4pi_rho ( i)
      lf_eigenvec1 ( i, 6, 7) = lf_eigenvec1 ( i, 2, 7)


! ----------

! Entropy wave

! ----------

      rt_eigenvec1 ( i, 1, 4) = 1.0
!      rt_eigenvec1 ( i, 2, 4) = 0.0
!      rt_eigenvec1 ( i, 3, 4) = 0.0
!      rt_eigenvec1 ( i, 4, 4) = 0.0
      rt_eigenvec1 ( i, 5, 4) = - xbar ( i)
!      rt_eigenvec1 ( i, 6, 4) = 0.0
!      rt_eigenvec1 ( i, 7, 4) = 0.0

! ----------

      lf_eigenvec1 ( i, 4, 1) = 1.0 - xbar ( i) / c_s ( i)**2
!      lf_eigenvec1 ( i, 4, 2) = 0.0
!      lf_eigenvec1 ( i, 4, 3) = 0.0
!      lf_eigenvec1 ( i, 4, 4) = 0.0
      lf_eigenvec1 ( i, 4, 5) = - 1.0 / c_s ( i)**2
!      lf_eigenvec1 ( i, 4, 6) = 0.0
!      lf_eigenvec1 ( i, 4, 7) = 0.0

! ----------

      END DO

! ----------------------------------------------------------------------

! Make right & left eigenvectors for slow and fast waves.

      DO i = nriembegin, nriemend

! ----------

! Fast waves

! ----------

      rt_eigenvec1 ( i, 1, 1) = alpha_f ( i) * rhobar ( i)
      rt_eigenvec1 ( i, 1, 7) = rt_eigenvec1 ( i, 1, 1)

      rt_eigenvec1 ( i, 2, 1) = - alpha_f ( i) * v_f ( i)
      rt_eigenvec1 ( i, 2, 7) = - rt_eigenvec1 ( i, 2, 1)

      rt_eigenvec1 ( i, 3, 1) = alpha_s ( i) * v_s ( i)
     1                        * beta_y ( i) * sgnbx ( i)
      rt_eigenvec1 ( i, 3, 7) = - rt_eigenvec1 ( i, 3, 1)

      rt_eigenvec1 ( i, 4, 1) = alpha_s ( i) * v_s ( i)
     1                        * beta_z ( i) * sgnbx ( i)
      rt_eigenvec1 ( i, 4, 7) = - rt_eigenvec1 ( i, 4, 1)

      rt_eigenvec1 ( i, 5, 1) = alpha_f ( i) * rhobar ( i)
     1                        * ( c_s ( i)**2 - xbar ( i))
      rt_eigenvec1 ( i, 5, 7) = rt_eigenvec1 ( i, 5, 1)

      rt_eigenvec1 ( i, 6, 1) = alpha_s ( i) * sq_4pi_rho ( i)
     1                        * c_s ( i) * beta_y ( i)
      rt_eigenvec1 ( i, 6, 7) = rt_eigenvec1 ( i, 6, 1)

      rt_eigenvec1 ( i, 7, 1) = alpha_s ( i) * sq_4pi_rho ( i)
     1                        * c_s ( i) * beta_z ( i)
      rt_eigenvec1 ( i, 7, 7) = rt_eigenvec1 ( i, 7, 1)

! ----------

      temparr_01 ( i) = 0.5 / c_s ( i)**2

      lf_eigenvec1 ( i, 1, 1) = alpha_f ( i) * xbar ( i)
     1                        * temparr_01 ( i) / rhobar ( i)
      lf_eigenvec1 ( i, 7, 1) = lf_eigenvec1 ( i, 1, 1)

      lf_eigenvec1 ( i, 1, 2) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 2, 1)
      lf_eigenvec1 ( i, 7, 2) = - lf_eigenvec1 ( i, 1, 2)

      lf_eigenvec1 ( i, 1, 3) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 3, 1)
      lf_eigenvec1 ( i, 7, 3) = - lf_eigenvec1 ( i, 1, 3)

      lf_eigenvec1 ( i, 1, 4) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 4, 1)
      lf_eigenvec1 ( i, 7, 4) = - lf_eigenvec1 ( i, 1, 4)

      lf_eigenvec1 ( i, 1, 5) = temparr_01 ( i)
     1                        * alpha_f ( i) / rhobar ( i)
      lf_eigenvec1 ( i, 7, 5) = lf_eigenvec1 ( i, 1, 5)

      lf_eigenvec1 ( i, 1, 6) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 6, 1)
     1                        * one_by_4pi_rho ( i)
      lf_eigenvec1 ( i, 7, 6) = lf_eigenvec1 ( i, 1, 6)

      lf_eigenvec1 ( i, 1, 7) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 7, 1)
     1                        * one_by_4pi_rho ( i)
      lf_eigenvec1 ( i, 7, 7) = lf_eigenvec1 ( i, 1, 7)

! ----------

! Slow waves

! ----------

      rt_eigenvec1 ( i, 1, 3) = alpha_s ( i) * rhobar ( i)
      rt_eigenvec1 ( i, 1, 5) = rt_eigenvec1 ( i, 1, 3)

      rt_eigenvec1 ( i, 2, 3) = - alpha_s ( i) * v_s ( i)
      rt_eigenvec1 ( i, 2, 5) = - rt_eigenvec1 ( i, 2, 3)

      rt_eigenvec1 ( i, 3, 3) = - alpha_f ( i) * v_f ( i)
     1                        * beta_y ( i) * sgnbx ( i)
      rt_eigenvec1 ( i, 3, 5) = - rt_eigenvec1 ( i, 3, 3)

      rt_eigenvec1 ( i, 4, 3) = - alpha_f ( i) * v_f ( i)
     1                        * beta_z ( i) * sgnbx ( i)
      rt_eigenvec1 ( i, 4, 5) = - rt_eigenvec1 ( i, 4, 3)

      rt_eigenvec1 ( i, 5, 3) = alpha_s ( i) * rhobar ( i)
     1                        * ( c_s ( i)**2 - xbar ( i))
      rt_eigenvec1 ( i, 5, 5) = rt_eigenvec1 ( i, 5, 3)

      rt_eigenvec1 ( i, 6, 3) = - alpha_f ( i) * sq_4pi_rho ( i)
     1                        * c_s ( i) * beta_y ( i)
      rt_eigenvec1 ( i, 6, 5) = rt_eigenvec1 ( i, 6, 3)

      rt_eigenvec1 ( i, 7, 3) = - alpha_f ( i) * sq_4pi_rho ( i)
     1                        * c_s ( i) * beta_z ( i)
      rt_eigenvec1 ( i, 7, 5) = rt_eigenvec1 ( i, 7, 3)

! ----------

      temparr_01 ( i) = 0.5 / c_s ( i)**2

      lf_eigenvec1 ( i, 3, 1) = alpha_s ( i) * xbar ( i)
     1                        * temparr_01 ( i) / rhobar ( i)
      lf_eigenvec1 ( i, 5, 1) = lf_eigenvec1 ( i, 3, 1)

      lf_eigenvec1 ( i, 3, 2) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 2, 3)
      lf_eigenvec1 ( i, 5, 2) = - lf_eigenvec1 ( i, 3, 2)

      lf_eigenvec1 ( i, 3, 3) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 3, 3)
      lf_eigenvec1 ( i, 5, 3) = - lf_eigenvec1 ( i, 3, 3)

      lf_eigenvec1 ( i, 3, 4) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 4, 3)
      lf_eigenvec1 ( i, 5, 4) = - lf_eigenvec1 ( i, 3, 4)

      lf_eigenvec1 ( i, 3, 5) = temparr_01 ( i)
     1                        * alpha_s ( i) / rhobar ( i)
      lf_eigenvec1 ( i, 5, 5) = lf_eigenvec1 ( i, 3, 5)

      lf_eigenvec1 ( i, 3, 6) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 6, 3)
     1                        * one_by_4pi_rho ( i)
      lf_eigenvec1 ( i, 5, 6) = lf_eigenvec1 ( i, 3, 6)

      lf_eigenvec1 ( i, 3, 7) = temparr_01 ( i)
     1                        * rt_eigenvec1 ( i, 7, 3)
     1                        * one_by_4pi_rho ( i)
      lf_eigenvec1 ( i, 5, 7) = lf_eigenvec1 ( i, 3, 7)


! ----------

      END DO

! ----------

! For the multiple species, fix up the eigenvectors that carry the
! contributions from the species fractions. The zero components of
! the eigenvectors have already been previously assigned.

#if ( NFLUID >= 1)

      DO j = NFIELD + 1, NFIELD + NFLUID

      DO i = nriembegin, nriemend

      rt_eigenvec1 ( i, j, j) = 1.0
      lf_eigenvec1 ( i, j, j) = 1.0

      END DO

      END DO

#endif

! ----------------------------------------------------------------------

! At this point the eigenvectors in the primitive variables are made
! and if only that much is desired, we can return at this point.

      IF ( iprimitive .EQ. 1) THEN

      DO j = 1, NFIELD + NFLUID
      DO k = 1, NFIELD + NFLUID
      DO i = nriembegin, nriemend
      rt_eigenvec ( i, j, k) = rt_eigenvec1 ( i, j, k)
      lf_eigenvec ( i, j, k) = lf_eigenvec1 ( i, j, k)
      END DO
      END DO
      END DO

      RETURN

      END IF

! ----------------------------------------------------------------------

      IF ( iprimitive .EQ. 0) THEN

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. First for
! the right eigenvectors. Don't make those matrix elements that are zero.

      DO i = nriembegin, nriemend


      dudw ( i, 1, 1) = 1.0
!      dudw ( i, 1, 2) = 0.0
!      dudw ( i, 1, 3) = 0.0
!      dudw ( i, 1, 4) = 0.0
!      dudw ( i, 1, 5) = 0.0
!      dudw ( i, 1, 6) = 0.0
!      dudw ( i, 1, 7) = 0.0


      dudw ( i, 2, 1) = vlxbar ( i)
      dudw ( i, 2, 2) = rhobar ( i)
!      dudw ( i, 2, 3) = 0.0
!      dudw ( i, 2, 4) = 0.0
!      dudw ( i, 2, 5) = 0.0
!      dudw ( i, 2, 6) = 0.0
!      dudw ( i, 2, 7) = 0.0


      dudw ( i, 3, 1) = vlybar ( i)
!      dudw ( i, 3, 2) = 0.0
      dudw ( i, 3, 3) = rhobar ( i)
!      dudw ( i, 3, 4) = 0.0
!      dudw ( i, 3, 5) = 0.0
!      dudw ( i, 3, 6) = 0.0
!      dudw ( i, 3, 7) = 0.0


      dudw ( i, 4, 1) = vlzbar ( i)
!      dudw ( i, 4, 2) = 0.0
!      dudw ( i, 4, 3) = 0.0
      dudw ( i, 4, 4) = rhobar ( i)
!      dudw ( i, 4, 5) = 0.0
!      dudw ( i, 4, 6) = 0.0
!      dudw ( i, 4, 7) = 0.0


      dudw ( i, 5, 1) = 0.5 * ( vlxbar ( i)**2
     1                        + vlybar ( i)**2 + vlzbar ( i)**2)
     1                + xbar ( i) + deng_drho ( i)
      dudw ( i, 5, 2) = rhobar ( i) * vlxbar ( i)
      dudw ( i, 5, 3) = rhobar ( i) * vlybar ( i)
      dudw ( i, 5, 4) = rhobar ( i) * vlzbar ( i)
      dudw ( i, 5, 5) = deng_dprs ( i)
      dudw ( i, 5, 6) = bfybar ( i) / fourpi
      dudw ( i, 5, 7) = bfzbar ( i) / fourpi


!      dudw ( i, 6, 1) = 0.0
!      dudw ( i, 6, 2) = 0.0
!      dudw ( i, 6, 3) = 0.0
!      dudw ( i, 6, 4) = 0.0
!      dudw ( i, 6, 5) = 0.0
      dudw ( i, 6, 6) = 1.0
!      dudw ( i, 6, 7) = 0.0


!      dudw ( i, 7, 1) = 0.0
!      dudw ( i, 7, 2) = 0.0
!      dudw ( i, 7, 3) = 0.0
!      dudw ( i, 7, 4) = 0.0
!      dudw ( i, 7, 5) = 0.0
!      dudw ( i, 7, 6) = 0.0
      dudw ( i, 7, 7) = 1.0


      END DO

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. Then for
! the left eigenvectors. Don't make those matrix elements that are zero.

      DO i = nriembegin, nriemend

      temparr_01 ( i) = 1.0 / rhobar ( i)


      dwdu ( i, 1, 1) = 1.0
      dwdu ( i, 2, 1) = - temparr_01 ( i) * vlxbar ( i)
      dwdu ( i, 3, 1) = - temparr_01 ( i) * vlybar ( i)
      dwdu ( i, 4, 1) = - temparr_01 ( i) * vlzbar ( i)
      dwdu ( i, 5, 1) = dprs_drho ( i) + dprs_deng ( i)
     1                * ( 0.5 * ( vlxbar ( i)**2 + vlybar ( i)**2
     1                          + vlzbar ( i)**2) - xbar ( i) )
!      dwdu ( i, 6, 1) = 0.0
!      dwdu ( i, 7, 1) = 0.0


!      dwdu ( i, 1, 2) = 0.0
      dwdu ( i, 2, 2) = temparr_01 ( i)
!      dwdu ( i, 3, 2) = 0.0
!      dwdu ( i, 4, 2) = 0.0
      dwdu ( i, 5, 2) = - dprs_deng ( i) * vlxbar ( i)
!      dwdu ( i, 6, 2) = 0.0
!      dwdu ( i, 7, 2) = 0.0


!      dwdu ( i, 1, 3) = 0.0
!      dwdu ( i, 2, 3) = 0.0
      dwdu ( i, 3, 3) = temparr_01 ( i)
!      dwdu ( i, 4, 3) = 0.0
      dwdu ( i, 5, 3) = - dprs_deng ( i) * vlybar ( i)
!      dwdu ( i, 6, 3) = 0.0
!      dwdu ( i, 7, 3) = 0.0


!      dwdu ( i, 1, 4) = 0.0
!      dwdu ( i, 2, 4) = 0.0
!      dwdu ( i, 3, 4) = 0.0
      dwdu ( i, 4, 4) = temparr_01 ( i)
      dwdu ( i, 5, 4) = - dprs_deng ( i) * vlzbar ( i)
!      dwdu ( i, 6, 4) = 0.0
!      dwdu ( i, 7, 4) = 0.0


!      dwdu ( i, 1, 5) = 0.0
!      dwdu ( i, 2, 5) = 0.0
!      dwdu ( i, 3, 5) = 0.0
!      dwdu ( i, 4, 5) = 0.0
      dwdu ( i, 5, 5) = dprs_deng ( i)
!      dwdu ( i, 6, 5) = 0.0
!      dwdu ( i, 7, 5) = 0.0


!      dwdu ( i, 1, 6) = 0.0
!      dwdu ( i, 2, 6) = 0.0
!      dwdu ( i, 3, 6) = 0.0
!      dwdu ( i, 4, 6) = 0.0
      dwdu ( i, 5, 6) = - dprs_deng ( i) * bfybar ( i)
     1                / fourpi
      dwdu ( i, 6, 6) = 1.0
!      dwdu ( i, 7, 6) = 0.0


!      dwdu ( i, 1, 7) = 0.0
!      dwdu ( i, 2, 7) = 0.0
!      dwdu ( i, 3, 7) = 0.0
!      dwdu ( i, 4, 7) = 0.0
      dwdu ( i, 5, 7) = - dprs_deng ( i) * bfzbar ( i)
     1                / fourpi
!      dwdu ( i, 6, 7) = 0.0
      dwdu ( i, 7, 7) = 1.0


      END DO

! ----------------------------------------------------------------------

! Do the matrix multiplications. Then save the eigenvalues. This
! gives eigenvectors in the conserved variables.
! Don't multiply those matrix elements that are zero.


      DO k = 1, NFIELD

      DO i = nriembegin, nriemend


      rt_eigenvec ( i, 1, k) =
     1     dudw ( i, 1, 1) * rt_eigenvec1 ( i, 1, k)

      rt_eigenvec ( i, 2, k) =
     1     dudw ( i, 2, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 2, 2) * rt_eigenvec1 ( i, 2, k)

      rt_eigenvec ( i, 3, k) =
     1     dudw ( i, 3, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 3, 3) * rt_eigenvec1 ( i, 3, k)

      rt_eigenvec ( i, 4, k) =
     1     dudw ( i, 4, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 4, 4) * rt_eigenvec1 ( i, 4, k)

      rt_eigenvec ( i, 5, k) =
     1     dudw ( i, 5, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 5, 2) * rt_eigenvec1 ( i, 2, k)
     1   + dudw ( i, 5, 3) * rt_eigenvec1 ( i, 3, k)
     1   + dudw ( i, 5, 4) * rt_eigenvec1 ( i, 4, k)
     1   + dudw ( i, 5, 5) * rt_eigenvec1 ( i, 5, k)
     1   + dudw ( i, 5, 6) * rt_eigenvec1 ( i, 6, k)
     1   + dudw ( i, 5, 7) * rt_eigenvec1 ( i, 7, k)

      rt_eigenvec ( i, 6, k) =
     1     dudw ( i, 6, 6) * rt_eigenvec1 ( i, 6, k)

      rt_eigenvec ( i, 7, k) =
     1     dudw ( i, 7, 7) * rt_eigenvec1 ( i, 7, k)


      lf_eigenvec ( i, k, 1) =
     1     lf_eigenvec1 ( i, k, 1) * dwdu ( i, 1, 1)
     1   + lf_eigenvec1 ( i, k, 2) * dwdu ( i, 2, 1)
     1   + lf_eigenvec1 ( i, k, 3) * dwdu ( i, 3, 1)
     1   + lf_eigenvec1 ( i, k, 4) * dwdu ( i, 4, 1)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 1)

      lf_eigenvec ( i, k, 2) =
     1     lf_eigenvec1 ( i, k, 2) * dwdu ( i, 2, 2)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 2)

      lf_eigenvec ( i, k, 3) =
     1     lf_eigenvec1 ( i, k, 3) * dwdu ( i, 3, 3)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 3)

      lf_eigenvec ( i, k, 4) =
     1     lf_eigenvec1 ( i, k, 4) * dwdu ( i, 4, 4)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 4)

      lf_eigenvec ( i, k, 5) =
     1     lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 5)

      lf_eigenvec ( i, k, 6) =
     1     lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 6)
     1   + lf_eigenvec1 ( i, k, 6) * dwdu ( i, 6, 6)

      lf_eigenvec ( i, k, 7) =
     1     lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 7)
     1   + lf_eigenvec1 ( i, k, 7) * dwdu ( i, 7, 7)


      END DO

      END DO

! At this point the eigenvectors in the conserved variables are made
! for a single species gas.

! ------------------------------

! For the situation where multiple species are present, fix up the
! eigenvectors. This consists of putting in the additional terms in
! the already existing seven MHD eigenvectors. Furthermore, it also
! consists of constructing the remaining eigenvectors that carry the
! species contributions.

#if ( NFLUID >= 1)

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

! Extra parts of the first seven right eigenvectors for MHD.
      rt_eigenvec ( i, k, 1) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 1)
      rt_eigenvec ( i, k, 2) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 2)
      rt_eigenvec ( i, k, 3) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 3)
      rt_eigenvec ( i, k, 4) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 4)
      rt_eigenvec ( i, k, 5) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 5)
      rt_eigenvec ( i, k, 6) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 6)
      rt_eigenvec ( i, k, 7) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 7)

! Do the entire right eigenvectors for each of the species.
      rt_eigenvec ( i, 1:4, k) = 0.0
      rt_eigenvec ( i, 5, k) = deng_dspecie ( ifluid, i)
      rt_eigenvec ( i, 6:k-1, k) = 0.0
      rt_eigenvec ( i, k, k) = rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      rt_eigenvec ( i, k+1:NFIELD + NFLUID, k) = 0.0
      END IF

! Extra parts of the first seven left eigenvectors for MHD.
      lf_eigenvec ( i, 1, k) = lf_eigenvec1 ( i, 1, 5)
     1                       * dprs_dspecie ( ifluid, i)
      lf_eigenvec ( i, 2, k) = lf_eigenvec1 ( i, 2, 5) 
     1                       * dprs_dspecie ( ifluid, i)
      lf_eigenvec ( i, 3, k) = lf_eigenvec1 ( i, 3, 5)
     1                       * dprs_dspecie ( ifluid, i)
      lf_eigenvec ( i, 4, k) = lf_eigenvec1 ( i, 4, 5)
     1                       * dprs_dspecie ( ifluid, i)
      lf_eigenvec ( i, 5, k) = lf_eigenvec1 ( i, 5, 5)
     1                       * dprs_dspecie ( ifluid, i)
      lf_eigenvec ( i, 6, k) = lf_eigenvec1 ( i, 6, 5)
     1                       * dprs_dspecie ( ifluid, i)
      lf_eigenvec ( i, 7, k) = lf_eigenvec1 ( i, 7, 5)
     1                       * dprs_dspecie ( ifluid, i)

! Do the entire left eigenvectors for each of the species.
      lf_eigenvec ( i, k, 1) = - speciebar ( ifluid, i)
     1                       / rhobar ( i)
      lf_eigenvec ( i, k, 2:k-1) = 0.0
      lf_eigenvec ( i, k, k) = 1.0 / rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      lf_eigenvec ( i, k, k+1:NFIELD + NFLUID) = 0.0
      END IF

      END DO

      END DO

#endif

! ----------------------------------------------------------------------

      END IF ! End of "iprimitive .EQ. 0"

! ----------------------------------------------------------------------

      IF ( iprimitive .EQ. -1) THEN

! ----------------------------------------------------------------------

! Make transformation matrices to the fake "conserved" variables. These
! are the variables for maintaining pressure positivity. First for
! the right eigenvectors. Don't make those matrix elements that are zero.

! Even when multiple species are used, these eigenvectors are based
! on a "frozen local gamma approximation". Thus, the ratio of specific
! heats, gamma, can vary from zone to zone but is assumed constant within
! each zone. This is justified because gamma is a very slowly varying
! function of the species fractions so it is ok to think of it as being
! (almost) exactly a constant.

      DO i = nriembegin, nriemend

      temparr_01 ( i) = 1.0 / rhobar ( i)
      temparr_02 ( i) = gamma_soundspeed ( i) - 1.0
      temparr_03 ( i) = rhobar ( i)**temparr_02 ( i)

      dudw ( i, 1, 1) = 1.0
!      dudw ( i, 1, 2) = 0.0
!      dudw ( i, 1, 3) = 0.0
!      dudw ( i, 1, 4) = 0.0
!      dudw ( i, 1, 5) = 0.0
!      dudw ( i, 1, 6) = 0.0
!      dudw ( i, 1, 7) = 0.0


      dudw ( i, 2, 1) = vlxbar ( i)
      dudw ( i, 2, 2) = rhobar ( i)
!      dudw ( i, 2, 3) = 0.0
!      dudw ( i, 2, 4) = 0.0
!      dudw ( i, 2, 5) = 0.0
!      dudw ( i, 2, 6) = 0.0
!      dudw ( i, 2, 7) = 0.0


      dudw ( i, 3, 1) = vlybar ( i)
!      dudw ( i, 3, 2) = 0.0
      dudw ( i, 3, 3) = rhobar ( i)
!      dudw ( i, 3, 4) = 0.0
!      dudw ( i, 3, 5) = 0.0
!      dudw ( i, 3, 6) = 0.0
!      dudw ( i, 3, 7) = 0.0


      dudw ( i, 4, 1) = vlzbar ( i)
!      dudw ( i, 4, 2) = 0.0
!      dudw ( i, 4, 3) = 0.0
      dudw ( i, 4, 4) = rhobar ( i)
!      dudw ( i, 4, 5) = 0.0
!      dudw ( i, 4, 6) = 0.0
!      dudw ( i, 4, 7) = 0.0


      dudw ( i, 5, 1) = - temparr_02 ( i) * prsbar ( i)
     1                * temparr_01 ( i) / temparr_03 ( i)
      dudw ( i, 5, 2) = 0.0
      dudw ( i, 5, 3) = 0.0
      dudw ( i, 5, 4) = 0.0
      dudw ( i, 5, 5) = 1.0 / temparr_03 ( i)
      dudw ( i, 5, 6) = 0.0
      dudw ( i, 5, 7) = 0.0


!      dudw ( i, 6, 1) = 0.0
!      dudw ( i, 6, 2) = 0.0
!      dudw ( i, 6, 3) = 0.0
!      dudw ( i, 6, 4) = 0.0
!      dudw ( i, 6, 5) = 0.0
      dudw ( i, 6, 6) = 1.0
!      dudw ( i, 6, 7) = 0.0


!      dudw ( i, 7, 1) = 0.0
!      dudw ( i, 7, 2) = 0.0
!      dudw ( i, 7, 3) = 0.0
!      dudw ( i, 7, 4) = 0.0
!      dudw ( i, 7, 5) = 0.0
!      dudw ( i, 7, 6) = 0.0
      dudw ( i, 7, 7) = 1.0


      END DO

! ----------------------------------------------------------------------

! Make transformation matrices to the fake "conserved" variables. These
! are the variables for maintaining pressure positivity. Then for
! the left eigenvectors. Don't make those matrix elements that are zero.

      DO i = nriembegin, nriemend

      temparr_01 ( i) = 1.0 / rhobar ( i)
      temparr_02 ( i) = gamma_soundspeed ( i) - 1.0
      temparr_03 ( i) = rhobar ( i)**temparr_02 ( i)


      dwdu ( i, 1, 1) = 1.0
      dwdu ( i, 2, 1) = - temparr_01 ( i) * vlxbar ( i)
      dwdu ( i, 3, 1) = - temparr_01 ( i) * vlybar ( i)
      dwdu ( i, 4, 1) = - temparr_01 ( i) * vlzbar ( i)
      dwdu ( i, 5, 1) = temparr_02 ( i) * prsbar ( i)
     1                / rhobar ( i)
!      dwdu ( i, 6, 1) = 0.0
!      dwdu ( i, 7, 1) = 0.0


!      dwdu ( i, 1, 2) = 0.0
      dwdu ( i, 2, 2) = temparr_01 ( i)
!      dwdu ( i, 3, 2) = 0.0
!      dwdu ( i, 4, 2) = 0.0
      dwdu ( i, 5, 2) = 0.0
!      dwdu ( i, 6, 2) = 0.0
!      dwdu ( i, 7, 2) = 0.0


!      dwdu ( i, 1, 3) = 0.0
!      dwdu ( i, 2, 3) = 0.0
      dwdu ( i, 3, 3) = temparr_01 ( i)
!      dwdu ( i, 4, 3) = 0.0
      dwdu ( i, 5, 3) = 0.0
!      dwdu ( i, 6, 3) = 0.0
!      dwdu ( i, 7, 3) = 0.0


!      dwdu ( i, 1, 4) = 0.0
!      dwdu ( i, 2, 4) = 0.0
!      dwdu ( i, 3, 4) = 0.0
      dwdu ( i, 4, 4) = temparr_01 ( i)
      dwdu ( i, 5, 4) = 0.0
!      dwdu ( i, 6, 4) = 0.0
!      dwdu ( i, 7, 4) = 0.0


!      dwdu ( i, 1, 5) = 0.0
!      dwdu ( i, 2, 5) = 0.0
!      dwdu ( i, 3, 5) = 0.0
!      dwdu ( i, 4, 5) = 0.0
      dwdu ( i, 5, 5) = temparr_03 ( i)
!      dwdu ( i, 6, 5) = 0.0
!      dwdu ( i, 7, 5) = 0.0


!      dwdu ( i, 1, 6) = 0.0
!      dwdu ( i, 2, 6) = 0.0
!      dwdu ( i, 3, 6) = 0.0
!      dwdu ( i, 4, 6) = 0.0
      dwdu ( i, 5, 6) = 0.0
      dwdu ( i, 6, 6) = 1.0
!      dwdu ( i, 7, 6) = 0.0


!      dwdu ( i, 1, 7) = 0.0
!      dwdu ( i, 2, 7) = 0.0
!      dwdu ( i, 3, 7) = 0.0
!      dwdu ( i, 4, 7) = 0.0
      dwdu ( i, 5, 7) = 0.0
!      dwdu ( i, 6, 7) = 0.0
      dwdu ( i, 7, 7) = 1.0


      END DO

! ----------------------------------------------------------------------

! Do the matrix multiplications. Then save the eigenvalues. This
! gives eigenvectors in the fake "conserved" variables. These
! are needed to maintain the pressure positivity.
! Don't multiply those matrix elements that are zero.


      DO k = 1, NFIELD

      DO i = nriembegin, nriemend


      rt_eigenvec ( i, 1, k) =
     1     dudw ( i, 1, 1) * rt_eigenvec1 ( i, 1, k)

      rt_eigenvec ( i, 2, k) =
     1     dudw ( i, 2, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 2, 2) * rt_eigenvec1 ( i, 2, k)

      rt_eigenvec ( i, 3, k) =
     1     dudw ( i, 3, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 3, 3) * rt_eigenvec1 ( i, 3, k)

      rt_eigenvec ( i, 4, k) =
     1     dudw ( i, 4, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 4, 4) * rt_eigenvec1 ( i, 4, k)

      rt_eigenvec ( i, 5, k) =
     1     dudw ( i, 5, 1) * rt_eigenvec1 ( i, 1, k)
     1   + dudw ( i, 5, 2) * rt_eigenvec1 ( i, 2, k)
     1   + dudw ( i, 5, 3) * rt_eigenvec1 ( i, 3, k)
     1   + dudw ( i, 5, 4) * rt_eigenvec1 ( i, 4, k)
     1   + dudw ( i, 5, 5) * rt_eigenvec1 ( i, 5, k)
     1   + dudw ( i, 5, 6) * rt_eigenvec1 ( i, 6, k)
     1   + dudw ( i, 5, 7) * rt_eigenvec1 ( i, 7, k)

      rt_eigenvec ( i, 6, k) =
     1     dudw ( i, 6, 6) * rt_eigenvec1 ( i, 6, k)

      rt_eigenvec ( i, 7, k) =
     1     dudw ( i, 7, 7) * rt_eigenvec1 ( i, 7, k)


      lf_eigenvec ( i, k, 1) =
     1     lf_eigenvec1 ( i, k, 1) * dwdu ( i, 1, 1)
     1   + lf_eigenvec1 ( i, k, 2) * dwdu ( i, 2, 1)
     1   + lf_eigenvec1 ( i, k, 3) * dwdu ( i, 3, 1)
     1   + lf_eigenvec1 ( i, k, 4) * dwdu ( i, 4, 1)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 1)

      lf_eigenvec ( i, k, 2) =
     1     lf_eigenvec1 ( i, k, 2) * dwdu ( i, 2, 2)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 2)

      lf_eigenvec ( i, k, 3) =
     1     lf_eigenvec1 ( i, k, 3) * dwdu ( i, 3, 3)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 3)

      lf_eigenvec ( i, k, 4) =
     1     lf_eigenvec1 ( i, k, 4) * dwdu ( i, 4, 4)
     1   + lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 4)

      lf_eigenvec ( i, k, 5) =
     1     lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 5)

      lf_eigenvec ( i, k, 6) =
     1     lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 6)
     1   + lf_eigenvec1 ( i, k, 6) * dwdu ( i, 6, 6)

      lf_eigenvec ( i, k, 7) =
     1     lf_eigenvec1 ( i, k, 5) * dwdu ( i, 5, 7)
     1   + lf_eigenvec1 ( i, k, 7) * dwdu ( i, 7, 7)


      END DO

      END DO

! At this point the eigenvectors in the fake "conserved" variables are made
! for a single species gas.

! ------------------------------

! For the situation where multiple species are present, fix up the
! eigenvectors. This consists of putting in the additional terms in
! the already existing seven MHD eigenvectors. Furthermore, it also
! consists of constructing the remaining eigenvectors that carry the
! species contributions.

#if ( NFLUID >= 1)

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

! Extra parts of the first seven right eigenvectors for MHD.
      rt_eigenvec ( i, k, 1) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 1)
      rt_eigenvec ( i, k, 2) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 2)
      rt_eigenvec ( i, k, 3) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 3)
      rt_eigenvec ( i, k, 4) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 4)
      rt_eigenvec ( i, k, 5) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 5)
      rt_eigenvec ( i, k, 6) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 6)
      rt_eigenvec ( i, k, 7) =
     1     speciebar ( ifluid, i) * rt_eigenvec1 ( i, 1, 7)

! Do the entire right eigenvectors for each of the species.
      rt_eigenvec ( i, 1:k-1, k) = 0.0
      rt_eigenvec ( i, k, k) = rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      rt_eigenvec ( i, k+1:NFIELD + NFLUID, k) = 0.0
      END IF

! Do the entire left eigenvectors for each of the species.
      lf_eigenvec ( i, k, 1) = - speciebar ( ifluid, i)
     1                       / rhobar ( i)
      lf_eigenvec ( i, k, 2:k-1) = 0.0
      lf_eigenvec ( i, k, k) = 1.0 / rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      lf_eigenvec ( i, k, k+1:NFIELD + NFLUID) = 0.0
      END IF

      END DO

      END DO

#endif

! ----------------------------------------------------------------------

      END IF ! End of "iprimitive .EQ. -1".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE EIGENVECMHD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EIGENVECMHD_PTWISE_8

     1             ( iprimitive, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, pi, 
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Returns normalized left and right eigenvectors for a specified set
! of mean states for adiabatic MHD. 

! When "iprimitive == 1" the eigenvectors expressed in terms of primitive
! variables are returned.

! When "iprimitive == 0" the eigenvectors in terms of conserved variables
! are returned.

! When "idirection = 1,2,3" we build eigenvectors in the x,y,z-directions.
! Velocity and magnetic field directions are permuted going in.
! The permutation is reversed in the eigenvector components, velocities
! and magnetic fields when leaving the subroutine.

! When primitive variables are used, the vector of variables inside this
! subroutine is in the order given by
! ( rho, vlx,y,z, prs, bfx,y,z, species_fractions)^T.

! When conserved variables are used, the vector of variables inside this
! subroutine is in the order given by
! ( rho, momx,y,z, eng, bfx,y,z, species_concentration)^T.

! Suitable permutations of the eigenvector components is done before exit
! when "idirection = 2, 3".

! Properly normalized eigenvectors as given in Roe and Balsara (1996)
! are generated.

! Ideal, non-reactive gas is assumed.

! ----------------------------------------------------------------------

      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER, INTENT ( IN) :: iprimitive, idirection

      INTEGER, INTENT ( INOUT) :: debug_stop

      REAL, INTENT ( IN) :: gamma, smallnum, prsfloor, rhofloor,
     1      pi, gasconst

      REAL, INTENT ( INOUT) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     bfxbar, bfybar, bfzbar, xbar

      REAL, DIMENSION ( 0: NFLUID), INTENT ( INOUT) :: speciebar

      REAL, DIMENSION ( NFIELD + NFLUID + 1), INTENT ( INOUT) ::
     1     eigenval

      REAL, DIMENSION ( NFIELD + NFLUID + 1, NFIELD + NFLUID + 1),
     1      INTENT ( INOUT) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

      INTEGER j, k, l

      REAL fourpi, small_speed


! Variables for making wave speeds and eigenvectors.

      REAL :: va_x, c_s, v_a, v_f, v_s, beta_y, beta_z,
     1     alpha_f, alpha_s, sgnbx, sq_4pi_rho, va_perp,
     1     eps_zone, one_by_4pi_rho,

     1     eng,
     1     gamma_eos, gamma_soundspeed,
     1     dprs_drho, dprs_deng, deng_drho, deng_dprs,

     1     temparr_01, temparr_02, temparr_03

! Temporary usage of left and right eigenvectors in the primitive variables.
! Transformation matrices to go from conserved to primitive and back.

      REAL, DIMENSION ( NFIELD + NFLUID + 1, NFIELD + NFLUID + 1) ::

     1     lf_eigenvec1, rt_eigenvec1, dudw, dwdu

! ----------------------------------------------------------------------

! Make some settings.

      fourpi = 4.0 * pi

! If the alfven speeds in the x or ( y, z) directions fall below "small_speed"
! times the local sound speed, manoeuvers have to be done to prevent the
! slow and fast eigenvectors from developing pathologies.

      small_speed = 0.000001

! ----------------------------------------------------------------------

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)"EIGENVECMHD_PTWISE_8 : NFIELD .NE. 7", NFIELD, 7
      debug_stop = 1
      STOP
      END IF

      IF ( ( iprimitive < 0) .OR. ( iprimitive > 1)) THEN
      WRITE ( 6, *)"EIGENVECMHD_PTWISE_8 : iprimitive = ", iprimitive
      debug_stop = 1
      STOP
      END IF

      IF ( ( idirection < 0) .OR. ( idirection > 3)) THEN
      WRITE ( 6, *)"EIGENVECMHD_PTWISE_8 : idirection = ", idirection
      debug_stop = 1
      STOP
      END IF

! ----------------------------------------------------------------------

! Rotate/permute variables for "idirection = 2 or 3".

      IF ( idirection == 2) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_02
      vlybar = temparr_03
      vlzbar = temparr_01

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_02
      bfybar = temparr_03
      bfzbar = temparr_01

      END IF

! ----------

      IF ( idirection == 3) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_03
      vlybar = temparr_01
      vlzbar = temparr_02

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_03
      bfybar = temparr_01
      bfzbar = temparr_02

      END IF

! ----------------------------------------------------------------------

! Bound the density and pressure variables. Divide "xbar" by 4 * pi.
! Set "xbar = 0.0" for pure eigenvectors.

      rhobar = AMAX1 ( rhobar, rhofloor)
      prsbar = AMAX1 ( prsbar, prsfloor)

      xbar = xbar / fourpi

! ----------------------------------------------------------------------

! Fix up "gamma_soundspeed" and partial derivatives of various
! thermodynamic quantities wrt one another.

      eng = prsbar / (gamma - 1.0)
      gamma_eos = gamma
      gamma_soundspeed = gamma

      dprs_drho = 0.0
      dprs_deng = gamma - 1.0
      deng_drho = 0.0
      deng_dprs = 1.0 / ( gamma - 1.0)

! ----------------------------------------------------------------------

! Make eigenvalues:

! ----------

      sgnbx = SIGN ( 1.0, bfxbar )

      sq_4pi_rho = SQRT ( fourpi * rhobar )

! ----------

      one_by_4pi_rho = 1.0 / ( fourpi * rhobar )

      va_x = bfxbar**2 * one_by_4pi_rho
      va_x = SQRT ( va_x )

      c_s = gamma_soundspeed * prsbar / rhobar
      c_s = SQRT ( c_s )

      v_a = ( bfxbar**2 + bfybar**2 + bfzbar**2) * one_by_4pi_rho
      v_a = SQRT ( v_a )

      temparr_01 = v_a**2 + c_s**2
      temparr_02 = temparr_01**2 - 4.0 * va_x**2 * c_s**2
      temparr_02 = AMAX1 ( temparr_02, 0.0)
      temparr_02 = SQRT ( temparr_02)

      v_f = 0.5 * ( temparr_01 + temparr_02 )
      v_f = SQRT ( v_f )

      v_s = 0.5 * ( temparr_01 - temparr_02 )
      v_s = AMAX1 ( v_s, 0.0)
      v_s = SQRT ( v_s )

! ----------

! Now make quantities that are essential for making the eigenvectors ( though
! inessential for just the eigenvalues).

      eps_zone = small_speed * c_s


      va_perp = SQRT ( bfybar**2 + bfzbar**2) / sq_4pi_rho


      IF ( va_perp .GT. eps_zone) THEN

      beta_y = bfybar / ( va_perp * sq_4pi_rho )
      beta_z = bfzbar / ( va_perp * sq_4pi_rho )

      ELSE

      beta_y = 1.0 / SQRT ( 2.0)
      beta_z = beta_y

      END IF

! ----------

      eigenval ( 1) = vlxbar - v_f
      eigenval ( 2) = vlxbar - va_x
      eigenval ( 3) = vlxbar - v_s
      eigenval ( 4) = vlxbar
      eigenval ( 5) = vlxbar + v_s
      eigenval ( 6) = vlxbar + va_x
      eigenval ( 7) = vlxbar + v_f
      eigenval ( 8) = 0.0

! ----------

#if ( NFLUID >= 1)
      eigenval ( NFIELD + 2: NFIELD + NFLUID + 1) = vlxbar
#endif

! ----------------------------------------------------------------------

! Make a consistent set of "alpha_f/s" for all the cases.

! ----------

      IF ( v_f .GE. v_s + eps_zone ) THEN

! Straightforward evaluation is adequate here. We are far away from
! the triple umbilic point.

      temparr_01 = 1.0 / ( v_f**2 - v_s**2)

      alpha_f = AMAX1 ( ( c_s**2 - v_s**2), 0.0 )
      alpha_f = SQRT ( temparr_01 * alpha_f)

      alpha_s = AMAX1 ( ( v_f**2 - c_s**2), 0.0 )
      alpha_s = SQRT ( temparr_01 * alpha_s)

      ELSE IF ( ABS ( va_x - c_s) .GT. smallnum * c_s ) THEN

! We are very close to the triple umbilic point but not on it.

      temparr_01 = va_perp / ( va_x - c_s)
      temparr_01 = ATAN ( temparr_01) / 2.0

      alpha_f = ABS ( SIN ( temparr_01))
      alpha_s = ABS ( COS ( temparr_01))

      ELSE

! We are bang on the triple umbilic point ( indistinguishable
! upto machine accuracy). Now which way we resolve this could
! only depend on the structure of the waves in the adjoining zones
! and we don't address that issue here.

      temparr_01 = 0.25 * pi * SIGN ( 1.0, va_x - c_s)

      alpha_f = ABS ( SIN ( temparr_01))
      alpha_s = ABS ( COS ( temparr_01))

      END IF

! ----------------------------------------------------------------------

! First zero out the eigenvectors. This is only needed because we can
! then fill in the species eigenvectors by putting 1's at suitable
! places on the diagonals.

      rt_eigenvec1 ( :, :) = 0.0
      lf_eigenvec1 ( :, :) = 0.0

! ----------

! Make right & left eigenvectors for alfven and entropy waves.

! ----------

! Alfven waves

! ----------

!      rt_eigenvec1 ( 1, 2) = 0.0
!      rt_eigenvec1 ( 1, 6) = 0.0

!      rt_eigenvec1 ( 2, 2) = 0.0
!      rt_eigenvec1 ( 2, 6) = 0.0

      rt_eigenvec1 ( 3, 2) = - beta_z
      rt_eigenvec1 ( 3, 6) = - rt_eigenvec1 ( 3, 2)

      rt_eigenvec1 ( 4, 2) = beta_y
      rt_eigenvec1 ( 4, 6) = - rt_eigenvec1 ( 4, 2)

!      rt_eigenvec1 ( 5, 2) = 0.0
!      rt_eigenvec1 ( 5, 6) = 0.0

!      rt_eigenvec1 ( 6, 2) = 0.0
!      rt_eigenvec1 ( 6, 6) = 0.0

      rt_eigenvec1 ( 7, 2) = - sgnbx * sq_4pi_rho * beta_z
      rt_eigenvec1 ( 7, 6) = rt_eigenvec1 ( 7, 2)

      rt_eigenvec1 ( 8, 2) = sgnbx * sq_4pi_rho * beta_y
      rt_eigenvec1 ( 8, 6) = rt_eigenvec1 ( 8, 2)

! ----------

!      lf_eigenvec1 ( 2, 1) = 0.0
!      lf_eigenvec1 ( 6, 1) = 0.0

!      lf_eigenvec1 ( 2, 2) = 0.0
!      lf_eigenvec1 ( 6, 2) = 0.0

      lf_eigenvec1 ( 2, 3) = - 0.5 * beta_z
      lf_eigenvec1 ( 6, 3) = - lf_eigenvec1 ( 2, 3)

      lf_eigenvec1 ( 2, 4) = 0.5 * beta_y
      lf_eigenvec1 ( 6, 4) = - lf_eigenvec1 ( 2, 4)

!      lf_eigenvec1 ( 2, 5) = 0.0
!      lf_eigenvec1 ( 6, 5) = 0.0

!      lf_eigenvec1 ( 2, 6) = 0.0
!      lf_eigenvec1 ( 6, 6) = 0.0

      lf_eigenvec1 ( 2, 7) = - 0.5 * sgnbx * beta_z / sq_4pi_rho
      lf_eigenvec1 ( 6, 7) = lf_eigenvec1 ( 2, 7)

      lf_eigenvec1 ( 2, 8) = 0.5 * sgnbx * beta_y / sq_4pi_rho
      lf_eigenvec1 ( 6, 8) = lf_eigenvec1 ( 2, 8)

! ----------

! Entropy wave

! ----------

      rt_eigenvec1 ( 1, 4) = 1.0
!      rt_eigenvec1 ( 2, 4) = 0.0
!      rt_eigenvec1 ( 3, 4) = 0.0
!      rt_eigenvec1 ( 4, 4) = 0.0
      rt_eigenvec1 ( 5, 4) = - xbar
!      rt_eigenvec1 ( 6, 4) = 0.0
!      rt_eigenvec1 ( 7, 4) = 0.0
!      rt_eigenvec1 ( 8, 4) = 0.0

! ----------

      lf_eigenvec1 ( 4, 1) = 1.0 - xbar / c_s**2
!      lf_eigenvec1 ( 4, 2) = 0.0
!      lf_eigenvec1 ( 4, 3) = 0.0
!      lf_eigenvec1 ( 4, 4) = 0.0
      lf_eigenvec1 ( 4, 5) = - 1.0 / c_s**2
!      lf_eigenvec1 ( 4, 6) = 0.0
!      lf_eigenvec1 ( 4, 7) = 0.0
!      lf_eigenvec1 ( 4, 8) = 0.0

! ----------

! Bx component of field

! ----------

!      rt_eigenvec1 ( 1, 8) = 0.0
!      rt_eigenvec1 ( 2, 8) = 0.0
!      rt_eigenvec1 ( 3, 8) = 0.0
!      rt_eigenvec1 ( 4, 8) = 0.0
!      rt_eigenvec1 ( 5, 8) = 0.0
      rt_eigenvec1 ( 6, 8) = 1.0
!      rt_eigenvec1 ( 7, 8) = 0.0
!      rt_eigenvec1 ( 8, 8) = 0.0

! ----------

!      lf_eigenvec1 ( 8, 1) = 0.0
!      lf_eigenvec1 ( 8, 2) = 0.0
!      lf_eigenvec1 ( 8, 3) = 0.0
!      lf_eigenvec1 ( 8, 4) = 0.0
!      lf_eigenvec1 ( 8, 5) = 0.0
      lf_eigenvec1 ( 8, 6) = 1.0
!      lf_eigenvec1 ( 8, 7) = 0.0
!      lf_eigenvec1 ( 8, 8) = 0.0

! ----------------------------------------------------------------------

! Make right & left eigenvectors for slow and fast waves.

! ----------

! Fast waves

! ----------

      rt_eigenvec1 ( 1, 1) = alpha_f * rhobar
      rt_eigenvec1 ( 1, 7) = rt_eigenvec1 ( 1, 1)

      rt_eigenvec1 ( 2, 1) = - alpha_f * v_f
      rt_eigenvec1 ( 2, 7) = - rt_eigenvec1 ( 2, 1)

      rt_eigenvec1 ( 3, 1) = alpha_s * v_s * beta_y * sgnbx
      rt_eigenvec1 ( 3, 7) = - rt_eigenvec1 ( 3, 1)

      rt_eigenvec1 ( 4, 1) = alpha_s * v_s * beta_z * sgnbx
      rt_eigenvec1 ( 4, 7) = - rt_eigenvec1 ( 4, 1)

      rt_eigenvec1 ( 5, 1) = alpha_f * rhobar * ( c_s**2 - xbar)
      rt_eigenvec1 ( 5, 7) = rt_eigenvec1 ( 5, 1)

!      rt_eigenvec1 ( 6, 1) = 0.0
!      rt_eigenvec1 ( 6, 7) = 0.0

      rt_eigenvec1 ( 7, 1) = alpha_s * sq_4pi_rho * c_s * beta_y
      rt_eigenvec1 ( 7, 7) = rt_eigenvec1 ( 7, 1)

      rt_eigenvec1 ( 8, 1) = alpha_s * sq_4pi_rho * c_s * beta_z
      rt_eigenvec1 ( 8, 7) = rt_eigenvec1 ( 8, 1)

! ----------

      temparr_01 = 0.5 / c_s**2

      lf_eigenvec1 ( 1, 1) = alpha_f * xbar * temparr_01 / rhobar
      lf_eigenvec1 ( 7, 1) = lf_eigenvec1 ( 1, 1)

      lf_eigenvec1 ( 1, 2) = temparr_01 * rt_eigenvec1 ( 2, 1)
      lf_eigenvec1 ( 7, 2) = - lf_eigenvec1 ( 1, 2)

      lf_eigenvec1 ( 1, 3) = temparr_01 * rt_eigenvec1 ( 3, 1)
      lf_eigenvec1 ( 7, 3) = - lf_eigenvec1 ( 1, 3)

      lf_eigenvec1 ( 1, 4) = temparr_01 * rt_eigenvec1 ( 4, 1)
      lf_eigenvec1 ( 7, 4) = - lf_eigenvec1 ( 1, 4)

      lf_eigenvec1 ( 1, 5) = temparr_01 * alpha_f / rhobar
      lf_eigenvec1 ( 7, 5) = lf_eigenvec1 ( 1, 5)

!      lf_eigenvec1 ( 1, 6) = 0.0
!      lf_eigenvec1 ( 7, 6) = 0.0

      lf_eigenvec1 ( 1, 7) = temparr_01 * rt_eigenvec1 ( 7, 1)
     1                     * one_by_4pi_rho
      lf_eigenvec1 ( 7, 7) = lf_eigenvec1 ( 1, 7)

      lf_eigenvec1 ( 1, 8) = temparr_01 * rt_eigenvec1 ( 8, 1)
     1                     * one_by_4pi_rho
      lf_eigenvec1 ( 7, 8) = lf_eigenvec1 ( 1, 8)

! ----------

! Slow waves

! ----------

      rt_eigenvec1 ( 1, 3) = alpha_s * rhobar
      rt_eigenvec1 ( 1, 5) = rt_eigenvec1 ( 1, 3)

      rt_eigenvec1 ( 2, 3) = - alpha_s * v_s
      rt_eigenvec1 ( 2, 5) = - rt_eigenvec1 ( 2, 3)

      rt_eigenvec1 ( 3, 3) = - alpha_f * v_f * beta_y * sgnbx
      rt_eigenvec1 ( 3, 5) = - rt_eigenvec1 ( 3, 3)

      rt_eigenvec1 ( 4, 3) = - alpha_f * v_f * beta_z * sgnbx
      rt_eigenvec1 ( 4, 5) = - rt_eigenvec1 ( 4, 3)

      rt_eigenvec1 ( 5, 3) = alpha_s * rhobar * ( c_s**2 - xbar)
      rt_eigenvec1 ( 5, 5) = rt_eigenvec1 ( 5, 3)

!      rt_eigenvec1 ( 6, 3) = 0.0
!      rt_eigenvec1 ( 6, 5) = 0.0

      rt_eigenvec1 ( 7, 3) = - alpha_f * sq_4pi_rho * c_s * beta_y
      rt_eigenvec1 ( 7, 5) = rt_eigenvec1 ( 7, 3)

      rt_eigenvec1 ( 8, 3) = - alpha_f * sq_4pi_rho * c_s * beta_z
      rt_eigenvec1 ( 8, 5) = rt_eigenvec1 ( 8, 3)

! ----------

      temparr_01 = 0.5 / c_s**2

      lf_eigenvec1 ( 3, 1) = alpha_s * xbar * temparr_01 / rhobar
      lf_eigenvec1 ( 5, 1) = lf_eigenvec1 ( 3, 1)

      lf_eigenvec1 ( 3, 2) = temparr_01 * rt_eigenvec1 ( 2, 3)
      lf_eigenvec1 ( 5, 2) = - lf_eigenvec1 ( 3, 2)

      lf_eigenvec1 ( 3, 3) = temparr_01 * rt_eigenvec1 ( 3, 3)
      lf_eigenvec1 ( 5, 3) = - lf_eigenvec1 ( 3, 3)

      lf_eigenvec1 ( 3, 4) = temparr_01 * rt_eigenvec1 ( 4, 3)
      lf_eigenvec1 ( 5, 4) = - lf_eigenvec1 ( 3, 4)

      lf_eigenvec1 ( 3, 5) = temparr_01 * alpha_s / rhobar
      lf_eigenvec1 ( 5, 5) = lf_eigenvec1 ( 3, 5)

!      lf_eigenvec1 ( 3, 6) = 0.0
!      lf_eigenvec1 ( 5, 6) = 0.0

      lf_eigenvec1 ( 3, 7) = temparr_01 * rt_eigenvec1 ( 7, 3)
     1                     * one_by_4pi_rho
      lf_eigenvec1 ( 5, 7) = lf_eigenvec1 ( 3, 7)

      lf_eigenvec1 ( 3, 8) = temparr_01 * rt_eigenvec1 ( 8, 3)
     1                     * one_by_4pi_rho
      lf_eigenvec1 ( 5, 8) = lf_eigenvec1 ( 3, 8)


! ----------

! For the multiple species, fix up the eigenvectors that carry the
! contributions from the species fractions. The zero components of
! the eigenvectors have already been previously assigned.

#if ( NFLUID >= 1)

      DO j = NFIELD + 2, NFIELD + NFLUID + 1

      rt_eigenvec1 ( j, j) = 1.0
      lf_eigenvec1 ( j, j) = 1.0

      END DO

#endif

! ----------------------------------------------------------------------

! At this point the eigenvectors in the primitive variables are made
! and if only that much is desired, we can permute and return at this point.

      IF ( iprimitive .EQ. 1) THEN

! ----------------------------------------------------------------------

      rt_eigenvec ( :, :) = rt_eigenvec1 ( :, :)
      lf_eigenvec ( :, :) = lf_eigenvec1 ( :, :)

! ----------------------------------------------------------------------

      ELSE  ! End of building eigenvectors in primitive variables.

! ----------------------------------------------------------------------

      dudw ( :, :) = 0.0
      dwdu ( :, :) = 0.0

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. First for
! the right eigenvectors. Don't make those matrix elements that are zero.

! Here we have "dwdu * dudw = I" and also "dudw * dwdu = I" where I is the
! identity matrix. Consequently, if "lf_eigenvec1 * rt_eigenvec1 = I" for
! the eigenvectors in primitive form then we also have
! "(lf_eigenvec1 * dwdu) * (dudw * rt_eigenvec1) = I".

! Derivatives of U(1) w.r.t. V(1:8)
      dudw ( 1, 1) = 1.0
!      dudw ( 1, 2) = 0.0
!      dudw ( 1, 3) = 0.0
!      dudw ( 1, 4) = 0.0
!      dudw ( 1, 5) = 0.0
!      dudw ( 1, 6) = 0.0
!      dudw ( 1, 7) = 0.0
!      dudw ( 1, 8) = 0.0


! Derivatives of U(2) w.r.t. V(1:8)
      dudw ( 2, 1) = vlxbar
      dudw ( 2, 2) = rhobar
!      dudw ( 2, 3) = 0.0
!      dudw ( 2, 4) = 0.0
!      dudw ( 2, 5) = 0.0
!      dudw ( 2, 6) = 0.0
!      dudw ( 2, 7) = 0.0
!      dudw ( 2, 8) = 0.0


      dudw ( 3, 1) = vlybar
!      dudw ( 3, 2) = 0.0
      dudw ( 3, 3) = rhobar
!      dudw ( 3, 4) = 0.0
!      dudw ( 3, 5) = 0.0
!      dudw ( 3, 6) = 0.0
!      dudw ( 3, 7) = 0.0
!      dudw ( 3, 8) = 0.0


      dudw ( 4, 1) = vlzbar
!      dudw ( 4, 2) = 0.0
!      dudw ( 4, 3) = 0.0
      dudw ( 4, 4) = rhobar
!      dudw ( 4, 5) = 0.0
!      dudw ( 4, 6) = 0.0
!      dudw ( 4, 7) = 0.0
!      dudw ( 4, 8) = 0.0


      dudw ( 5, 1) = 0.5 * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             + xbar + deng_drho
      dudw ( 5, 2) = rhobar * vlxbar
      dudw ( 5, 3) = rhobar * vlybar
      dudw ( 5, 4) = rhobar * vlzbar
      dudw ( 5, 5) = deng_dprs
      dudw ( 5, 6) = bfxbar / fourpi
      dudw ( 5, 7) = bfybar / fourpi
      dudw ( 5, 8) = bfzbar / fourpi


!      dudw ( 6, 1) = 0.0
!      dudw ( 6, 2) = 0.0
!      dudw ( 6, 3) = 0.0
!      dudw ( 6, 4) = 0.0
!      dudw ( 6, 5) = 0.0
      dudw ( 6, 6) = 1.0
!      dudw ( 6, 7) = 0.0
!      dudw ( 6, 8) = 0.0


!      dudw ( 7, 1) = 0.0
!      dudw ( 7, 2) = 0.0
!      dudw ( 7, 3) = 0.0
!      dudw ( 7, 4) = 0.0
!      dudw ( 7, 5) = 0.0
!      dudw ( 7, 6) = 0.0
      dudw ( 7, 7) = 1.0
!      dudw ( 7, 8) = 0.0


!      dudw ( 8, 1) = 0.0
!      dudw ( 8, 2) = 0.0
!      dudw ( 8, 3) = 0.0
!      dudw ( 8, 4) = 0.0
!      dudw ( 8, 5) = 0.0
!      dudw ( 8, 6) = 0.0
!      dudw ( 8, 7) = 0.0
      dudw ( 8, 8) = 1.0


#if ( NFLUID >= 1)
      DO k = 9, NFIELD + NFLUID + 1

      dudw ( k, 1) = speciebar ( k - 8)
      dudw ( k, k) = rhobar

      END DO
#endif

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. Then for
! the left eigenvectors. Don't make those matrix elements that are zero.

      temparr_01 = 1.0 / rhobar


! Derivatives of V(1:8) w.r.t. U(1)
      dwdu ( 1, 1) = 1.0
      dwdu ( 2, 1) = - temparr_01 * vlxbar
      dwdu ( 3, 1) = - temparr_01 * vlybar
      dwdu ( 4, 1) = - temparr_01 * vlzbar
      dwdu ( 5, 1) = dprs_drho + dprs_deng * ( 0.5 * ( vlxbar**2
     1                                 + vlybar**2 + vlzbar**2) - xbar)
!      dwdu ( 6, 1) = 0.0
!      dwdu ( 7, 1) = 0.0
!      dwdu ( 8, 1) = 0.0


! Derivatives of V(1:8) w.r.t. U(2)
!      dwdu ( 1, 2) = 0.0
      dwdu ( 2, 2) = temparr_01
!      dwdu ( 3, 2) = 0.0
!      dwdu ( 4, 2) = 0.0
      dwdu ( 5, 2) = - dprs_deng * vlxbar
!      dwdu ( 6, 2) = 0.0
!      dwdu ( 7, 2) = 0.0
!      dwdu ( 8, 2) = 0.0


!      dwdu ( 1, 3) = 0.0
!      dwdu ( 2, 3) = 0.0
      dwdu ( 3, 3) = temparr_01
!      dwdu ( 4, 3) = 0.0
      dwdu ( 5, 3) = - dprs_deng * vlybar
!      dwdu ( 6, 3) = 0.0
!      dwdu ( 7, 3) = 0.0
!      dwdu ( 8, 3) = 0.0


!      dwdu ( 1, 4) = 0.0
!      dwdu ( 2, 4) = 0.0
!      dwdu ( 3, 4) = 0.0
      dwdu ( 4, 4) = temparr_01
      dwdu ( 5, 4) = - dprs_deng * vlzbar
!      dwdu ( 6, 4) = 0.0
!      dwdu ( 7, 4) = 0.0
!      dwdu ( 8, 4) = 0.0


!      dwdu ( 1, 5) = 0.0
!      dwdu ( 2, 5) = 0.0
!      dwdu ( 3, 5) = 0.0
!      dwdu ( 4, 5) = 0.0
      dwdu ( 5, 5) = dprs_deng
!      dwdu ( 6, 5) = 0.0
!      dwdu ( 7, 5) = 0.0
!      dwdu ( 8, 5) = 0.0


!      dwdu ( 1, 6) = 0.0
!      dwdu ( 2, 6) = 0.0
!      dwdu ( 3, 6) = 0.0
!      dwdu ( 4, 6) = 0.0
      dwdu ( 5, 6) = - dprs_deng * bfxbar / fourpi
      dwdu ( 6, 6) = 1.0
!      dwdu ( 7, 6) = 0.0
!      dwdu ( 8, 6) = 0.0


!      dwdu ( 1, 7) = 0.0
!      dwdu ( 2, 7) = 0.0
!      dwdu ( 3, 7) = 0.0
!      dwdu ( 4, 7) = 0.0
      dwdu ( 5, 7) = - dprs_deng * bfybar / fourpi
!      dwdu ( 6, 7) = 0.0
      dwdu ( 7, 7) = 1.0
!      dwdu ( 8, 7) = 0.0


!      dwdu ( 1, 8) = 0.0
!      dwdu ( 2, 8) = 0.0
!      dwdu ( 3, 8) = 0.0
!      dwdu ( 4, 8) = 0.0
      dwdu ( 5, 8) = - dprs_deng * bfzbar / fourpi
!      dwdu ( 6, 8) = 0.0
!      dwdu ( 7, 8) = 0.0
      dwdu ( 8, 8) = 1.0


#if ( NFLUID >= 1)
      DO k = 9, NFIELD + NFLUID + 1

      dwdu ( k, 1) = - temparr_01 * speciebar ( k - 8)
      dwdu ( k, k) = temparr_01

      END DO
#endif

! ----------------------------------------------------------------------

! Do the matrix multiplications. Then save the eigenvalues. This
! gives eigenvectors in the conserved variables.
! Don't multiply those matrix elements that are zero.

      rt_eigenvec ( :, :) = 0.0
      lf_eigenvec ( :, :) = 0.0

! Loop over k^th right eigenvector or k^th left eigenvector.
      DO k = 1, NFIELD + NFLUID + 1 

      DO j = 1, NFIELD + NFLUID + 1

      DO l = 1, NFIELD + NFLUID + 1

      rt_eigenvec ( j, k) = rt_eigenvec ( j, k)
     1                    + dudw ( j, l) * rt_eigenvec1 ( l, k)

      lf_eigenvec ( k, j) = lf_eigenvec ( k, j)
     1                    + lf_eigenvec1 ( k, l) * dwdu ( l, j)

      END DO

      END DO

      END DO

! ----------------------------------------------------------------------

      END IF  ! End of building eigenvectors in conserved variables.

! ----------------------------------------------------------------------
! The permutation is reversed in the eigenvector components, velocities
! and magnetic fields when leaving the subroutine.

      IF ( idirection == 2) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_03
      vlybar = temparr_01
      vlzbar = temparr_02

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_03
      bfybar = temparr_01
      bfzbar = temparr_02

! -----

      DO k = 1, NFIELD + NFLUID + 1

! -----

      temparr_01 = rt_eigenvec ( 2, k)
      temparr_02 = rt_eigenvec ( 3, k)
      temparr_03 = rt_eigenvec ( 4, k)

      rt_eigenvec ( 2, k) = temparr_03
      rt_eigenvec ( 3, k) = temparr_01
      rt_eigenvec ( 4, k) = temparr_02

! -----

      temparr_01 = rt_eigenvec ( 6, k)
      temparr_02 = rt_eigenvec ( 7, k)
      temparr_03 = rt_eigenvec ( 8, k)

      rt_eigenvec ( 6, k) = temparr_03
      rt_eigenvec ( 7, k) = temparr_01
      rt_eigenvec ( 8, k) = temparr_02

! -----

      temparr_01 = lf_eigenvec ( k, 2)
      temparr_02 = lf_eigenvec ( k, 3)
      temparr_03 = lf_eigenvec ( k, 4)

      lf_eigenvec ( k, 2) = temparr_03
      lf_eigenvec ( k, 3) = temparr_01
      lf_eigenvec ( k, 4) = temparr_02

! -----

      temparr_01 = lf_eigenvec ( k, 6)
      temparr_02 = lf_eigenvec ( k, 7)
      temparr_03 = lf_eigenvec ( k, 8)

      lf_eigenvec ( k, 6) = temparr_03
      lf_eigenvec ( k, 7) = temparr_01
      lf_eigenvec ( k, 8) = temparr_02

! -----

      END DO

! -----

      END IF

! ----------

      IF ( idirection == 3) THEN

! -----

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_02
      vlybar = temparr_03
      vlzbar = temparr_01

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_02
      bfybar = temparr_03
      bfzbar = temparr_01

! -----

      DO k = 1, NFIELD + NFLUID + 1

! -----

      temparr_01 = rt_eigenvec ( 2, k)
      temparr_02 = rt_eigenvec ( 3, k)
      temparr_03 = rt_eigenvec ( 4, k)

      rt_eigenvec ( 2, k) = temparr_02
      rt_eigenvec ( 3, k) = temparr_03
      rt_eigenvec ( 4, k) = temparr_01

! -----

      temparr_01 = rt_eigenvec ( 6, k)
      temparr_02 = rt_eigenvec ( 7, k)
      temparr_03 = rt_eigenvec ( 8, k)

      rt_eigenvec ( 6, k) = temparr_02
      rt_eigenvec ( 7, k) = temparr_03
      rt_eigenvec ( 8, k) = temparr_01

! -----

      temparr_01 = lf_eigenvec ( k, 2)
      temparr_02 = lf_eigenvec ( k, 3)
      temparr_03 = lf_eigenvec ( k, 4)

      lf_eigenvec ( k, 2) = temparr_02
      lf_eigenvec ( k, 3) = temparr_03
      lf_eigenvec ( k, 4) = temparr_01

! -----

      temparr_01 = lf_eigenvec ( k, 6)
      temparr_02 = lf_eigenvec ( k, 7)
      temparr_03 = lf_eigenvec ( k, 8)

      lf_eigenvec ( k, 6) = temparr_02
      lf_eigenvec ( k, 7) = temparr_03
      lf_eigenvec ( k, 8) = temparr_01

! -----

      END DO

! -----

      END IF

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE EIGENVECMHD_PTWISE_8

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EIGENVECMHD_PTWISE_LINDEG

     1             ( iprimitive, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, pi, 
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Returns normalized left and right eigenvectors for a specified set
! of mean states for adiabatic MHD. Does this only for the linearly degenerate
! waves. Rest of the eigenvectors are zeroed. This is useful in endowing
! substructure in multid. Riemann problems.

! When "iprimitive == 1" the eigenvectors expressed in terms of primitive
! variables are returned.

! When "iprimitive == 0" the eigenvectors in terms of conserved variables
! are returned.

! When "idirection = 1,2,3" we build eigenvectors in the x,y,z-directions.
! Velocity and magnetic field directions are permuted going in.
! The permutation is reversed in the eigenvector components, velocities
! and magnetic fields when leaving the subroutine.

! When primitive variables are used, the vector of variables inside this
! subroutine is in the order given by
! ( rho, vlx,y,z, prs, bfx,y,z, species_fractions)^T.

! When conserved variables are used, the vector of variables inside this
! subroutine is in the order given by
! ( rho, momx,y,z, eng, bfx,y,z, species_concentration)^T.

! Suitable permutations of the eigenvector components is done before exit
! when "idirection = 2, 3".

! Properly normalized eigenvectors as given in Roe and Balsara (1996)
! are generated.

! Ideal, non-reactive gas is assumed.

! ----------------------------------------------------------------------

      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER, INTENT ( IN) :: iprimitive, idirection

      INTEGER, INTENT ( INOUT) :: debug_stop

      REAL, INTENT ( IN) :: gamma, smallnum, prsfloor, rhofloor,
     1      pi, gasconst

      REAL, INTENT ( INOUT) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     bfxbar, bfybar, bfzbar, xbar

      REAL, DIMENSION ( 0: NFLUID), INTENT ( INOUT) :: speciebar

      REAL, DIMENSION ( NFIELD + NFLUID + 1), INTENT ( INOUT) ::
     1     eigenval

      REAL, DIMENSION ( NFIELD + NFLUID + 1, NFIELD + NFLUID + 1),
     1      INTENT ( INOUT) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

      INTEGER j, k, l

      REAL fourpi, small_speed


! Variables for making wave speeds and eigenvectors.

      REAL :: va_x, c_s, v_a, v_f, v_s, beta_y, beta_z,
     1     alpha_f, alpha_s, sgnbx, sq_4pi_rho, va_perp,
     1     eps_zone, one_by_4pi_rho,

     1     eng,
     1     gamma_eos, gamma_soundspeed,
     1     dprs_drho, dprs_deng, deng_drho, deng_dprs,

     1     temparr_01, temparr_02, temparr_03

! Temporary usage of left and right eigenvectors in the primitive variables.
! Transformation matrices to go from conserved to primitive and back.

      REAL, DIMENSION ( NFIELD + NFLUID + 1, NFIELD + NFLUID + 1) ::

     1     lf_eigenvec1, rt_eigenvec1, dudw, dwdu

! ----------------------------------------------------------------------

! Make some settings.

      fourpi = 4.0 * pi

! If the alfven speeds in the x or ( y, z) directions fall below "small_speed"
! times the local sound speed, manoeuvers have to be done to prevent the
! slow and fast eigenvectors from developing pathologies.

      small_speed = 0.000001

! ----------------------------------------------------------------------

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1      "EIGENVECMHD_PTWISE_LINDEG : NFIELD .NE. 7", NFIELD, 7
      debug_stop = 1
      STOP
      END IF

      IF ( ( iprimitive < 0) .OR. ( iprimitive > 1)) THEN
      WRITE ( 6, *)
     1      "EIGENVECMHD_PTWISE_LINDEG : iprimitive = ", iprimitive
      debug_stop = 1
      STOP
      END IF

      IF ( ( idirection < 0) .OR. ( idirection > 3)) THEN
      WRITE ( 6, *)
     1      "EIGENVECMHD_PTWISE_LINDEG : idirection = ", idirection
      debug_stop = 1
      STOP
      END IF

! ----------------------------------------------------------------------

! Rotate/permute variables for "idirection = 2 or 3".

      IF ( idirection == 2) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_02
      vlybar = temparr_03
      vlzbar = temparr_01

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_02
      bfybar = temparr_03
      bfzbar = temparr_01

      END IF

! ----------

      IF ( idirection == 3) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_03
      vlybar = temparr_01
      vlzbar = temparr_02

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_03
      bfybar = temparr_01
      bfzbar = temparr_02

      END IF

! ----------------------------------------------------------------------

! Bound the density and pressure variables. Divide "xbar" by 4 * pi.
! Set "xbar = 0.0" for pure eigenvectors.

      rhobar = AMAX1 ( rhobar, rhofloor)
      prsbar = AMAX1 ( prsbar, prsfloor)

      xbar = xbar / fourpi

! ----------------------------------------------------------------------

! Fix up "gamma_soundspeed" and partial derivatives of various
! thermodynamic quantities wrt one another.

      eng = prsbar / (gamma - 1.0)
      gamma_eos = gamma
      gamma_soundspeed = gamma

      dprs_drho = 0.0
      dprs_deng = gamma - 1.0
      deng_drho = 0.0
      deng_dprs = 1.0 / ( gamma - 1.0)

! ----------------------------------------------------------------------

! Make eigenvalues:

! ----------

      sgnbx = SIGN ( 1.0, bfxbar )

      sq_4pi_rho = SQRT ( fourpi * rhobar )

! ----------

      one_by_4pi_rho = 1.0 / ( fourpi * rhobar )

      va_x = bfxbar**2 * one_by_4pi_rho
      va_x = SQRT ( va_x )

      c_s = gamma_soundspeed * prsbar / rhobar
      c_s = SQRT ( c_s )

      v_a = ( bfxbar**2 + bfybar**2 + bfzbar**2) * one_by_4pi_rho
      v_a = SQRT ( v_a )

      temparr_01 = v_a**2 + c_s**2
      temparr_02 = temparr_01**2 - 4.0 * va_x**2 * c_s**2
      temparr_02 = AMAX1 ( temparr_02, 0.0)
      temparr_02 = SQRT ( temparr_02)

      v_f = 0.5 * ( temparr_01 + temparr_02 )
      v_f = SQRT ( v_f )

      v_s = 0.5 * ( temparr_01 - temparr_02 )
      v_s = AMAX1 ( v_s, 0.0)
      v_s = SQRT ( v_s )

! ----------

! Now make quantities that are essential for making the eigenvectors ( though
! inessential for just the eigenvalues).

      eps_zone = small_speed * c_s


      va_perp = SQRT ( bfybar**2 + bfzbar**2) / sq_4pi_rho


      IF ( va_perp .GT. eps_zone) THEN

      beta_y = bfybar / ( va_perp * sq_4pi_rho )
      beta_z = bfzbar / ( va_perp * sq_4pi_rho )

      ELSE

      beta_y = 1.0 / SQRT ( 2.0)
      beta_z = beta_y

      END IF

! ----------

      eigenval ( 1) = vlxbar - v_f
      eigenval ( 2) = vlxbar - va_x
      eigenval ( 3) = vlxbar - v_s
      eigenval ( 4) = vlxbar
      eigenval ( 5) = vlxbar + v_s
      eigenval ( 6) = vlxbar + va_x
      eigenval ( 7) = vlxbar + v_f
      eigenval ( 8) = 0.0

! ----------

#if ( NFLUID >= 1)
      eigenval ( NFIELD + 2: NFIELD + NFLUID + 1) = vlxbar
#endif

! ----------------------------------------------------------------------

! Make a consistent set of "alpha_f/s" for all the cases.

! ----------

      IF ( v_f .GE. v_s + eps_zone ) THEN

! Straightforward evaluation is adequate here. We are far away from
! the triple umbilic point.

      temparr_01 = 1.0 / ( v_f**2 - v_s**2)

      alpha_f = AMAX1 ( ( c_s**2 - v_s**2), 0.0 )
      alpha_f = SQRT ( temparr_01 * alpha_f)

      alpha_s = AMAX1 ( ( v_f**2 - c_s**2), 0.0 )
      alpha_s = SQRT ( temparr_01 * alpha_s)

      ELSE IF ( ABS ( va_x - c_s) .GT. smallnum * c_s ) THEN

! We are very close to the triple umbilic point but not on it.

      temparr_01 = va_perp / ( va_x - c_s)
      temparr_01 = ATAN ( temparr_01) / 2.0

      alpha_f = ABS ( SIN ( temparr_01))
      alpha_s = ABS ( COS ( temparr_01))

      ELSE

! We are bang on the triple umbilic point ( indistinguishable
! upto machine accuracy). Now which way we resolve this could
! only depend on the structure of the waves in the adjoining zones
! and we don't address that issue here.

      temparr_01 = 0.25 * pi * SIGN ( 1.0, va_x - c_s)

      alpha_f = ABS ( SIN ( temparr_01))
      alpha_s = ABS ( COS ( temparr_01))

      END IF

! ----------------------------------------------------------------------

! First zero out the eigenvectors. This is only needed because we can
! then fill in the species eigenvectors by putting 1's at suitable
! places on the diagonals.

      rt_eigenvec1 ( :, :) = 0.0
      lf_eigenvec1 ( :, :) = 0.0

! ----------

! Make right & left eigenvectors for alfven and entropy waves.

! ----------

! Alfven waves

! If the slow waves are retained alongside with the Alfven waves then
! they form a complete eigenvector basis into which the full system can
! be projected even in the limit where the longitudinal magnetic field
! goes to zero.

! However, that is not what we are doing here. For that reason, when we
! discover that limit, we just revert to the hydrodynamical shear waves.

! ----------

      IF ( va_x > eps_zone) THEN

! ----------

!      rt_eigenvec1 ( 1, 2) = 0.0
!      rt_eigenvec1 ( 1, 6) = 0.0

!      rt_eigenvec1 ( 2, 2) = 0.0
!      rt_eigenvec1 ( 2, 6) = 0.0

      rt_eigenvec1 ( 3, 2) = - beta_z
      rt_eigenvec1 ( 3, 6) = - rt_eigenvec1 ( 3, 2)

      rt_eigenvec1 ( 4, 2) = beta_y
      rt_eigenvec1 ( 4, 6) = - rt_eigenvec1 ( 4, 2)

!      rt_eigenvec1 ( 5, 2) = 0.0
!      rt_eigenvec1 ( 5, 6) = 0.0

!      rt_eigenvec1 ( 6, 2) = 0.0
!      rt_eigenvec1 ( 6, 6) = 0.0

      rt_eigenvec1 ( 7, 2) = - sgnbx * sq_4pi_rho * beta_z
      rt_eigenvec1 ( 7, 6) = rt_eigenvec1 ( 7, 2)

      rt_eigenvec1 ( 8, 2) = sgnbx * sq_4pi_rho * beta_y
      rt_eigenvec1 ( 8, 6) = rt_eigenvec1 ( 8, 2)

! ----------

!      lf_eigenvec1 ( 2, 1) = 0.0
!      lf_eigenvec1 ( 6, 1) = 0.0

!      lf_eigenvec1 ( 2, 2) = 0.0
!      lf_eigenvec1 ( 6, 2) = 0.0

      lf_eigenvec1 ( 2, 3) = - 0.5 * beta_z
      lf_eigenvec1 ( 6, 3) = - lf_eigenvec1 ( 2, 3)

      lf_eigenvec1 ( 2, 4) = 0.5 * beta_y
      lf_eigenvec1 ( 6, 4) = - lf_eigenvec1 ( 2, 4)

!      lf_eigenvec1 ( 2, 5) = 0.0
!      lf_eigenvec1 ( 6, 5) = 0.0

!      lf_eigenvec1 ( 2, 6) = 0.0
!      lf_eigenvec1 ( 6, 6) = 0.0

      lf_eigenvec1 ( 2, 7) = - 0.5 * sgnbx * beta_z / sq_4pi_rho
      lf_eigenvec1 ( 6, 7) = lf_eigenvec1 ( 2, 7)

      lf_eigenvec1 ( 2, 8) = 0.5 * sgnbx * beta_y / sq_4pi_rho
      lf_eigenvec1 ( 6, 8) = lf_eigenvec1 ( 2, 8)

! ----------

      ELSE

! ----------

!      rt_eigenvec1 ( 1, 2) = 0.0
!      rt_eigenvec1 ( 1, 6) = 0.0

!      rt_eigenvec1 ( 2, 2) = 0.0
!      rt_eigenvec1 ( 2, 6) = 0.0

      rt_eigenvec1 ( 3, 2) = 1.0
      rt_eigenvec1 ( 3, 6) = 0.0

      rt_eigenvec1 ( 4, 2) = 0.0
      rt_eigenvec1 ( 4, 6) = 1.0

!      rt_eigenvec1 ( 5, 2) = 0.0
!      rt_eigenvec1 ( 5, 6) = 0.0

!      rt_eigenvec1 ( 6, 2) = 0.0
!      rt_eigenvec1 ( 6, 6) = 0.0

!      rt_eigenvec1 ( 7, 2) = 0.0
!      rt_eigenvec1 ( 7, 6) = 0.0

!      rt_eigenvec1 ( 8, 2) = 0.0
!      rt_eigenvec1 ( 8, 6) = 0.0

! ----------

!      lf_eigenvec1 ( 2, 1) = 0.0
!      lf_eigenvec1 ( 6, 1) = 0.0

!      lf_eigenvec1 ( 2, 2) = 0.0
!      lf_eigenvec1 ( 6, 2) = 0.0

      lf_eigenvec1 ( 2, 3) = 1.0
      lf_eigenvec1 ( 6, 3) = 0.0

      lf_eigenvec1 ( 2, 4) = 0.0
      lf_eigenvec1 ( 6, 4) = 1.0

!      lf_eigenvec1 ( 2, 5) = 0.0
!      lf_eigenvec1 ( 6, 5) = 0.0

!      lf_eigenvec1 ( 2, 6) = 0.0
!      lf_eigenvec1 ( 6, 6) = 0.0

!      lf_eigenvec1 ( 2, 7) = 0.0
!      lf_eigenvec1 ( 6, 7) = 0.0

!      lf_eigenvec1 ( 2, 8) = 0.0
!      lf_eigenvec1 ( 6, 8) = 0.0

! ----------

      END IF

! ----------

! Entropy wave

! ----------

      rt_eigenvec1 ( 1, 4) = 1.0
!      rt_eigenvec1 ( 2, 4) = 0.0
!      rt_eigenvec1 ( 3, 4) = 0.0
!      rt_eigenvec1 ( 4, 4) = 0.0
      rt_eigenvec1 ( 5, 4) = - xbar
!      rt_eigenvec1 ( 6, 4) = 0.0
!      rt_eigenvec1 ( 7, 4) = 0.0
!      rt_eigenvec1 ( 8, 4) = 0.0

! ----------

      lf_eigenvec1 ( 4, 1) = 1.0 - xbar / c_s**2
!      lf_eigenvec1 ( 4, 2) = 0.0
!      lf_eigenvec1 ( 4, 3) = 0.0
!      lf_eigenvec1 ( 4, 4) = 0.0
      lf_eigenvec1 ( 4, 5) = - 1.0 / c_s**2
!      lf_eigenvec1 ( 4, 6) = 0.0
!      lf_eigenvec1 ( 4, 7) = 0.0
!      lf_eigenvec1 ( 4, 8) = 0.0

! ----------

! Bx component of field

! ----------

!      rt_eigenvec1 ( 1, 8) = 0.0
!      rt_eigenvec1 ( 2, 8) = 0.0
!      rt_eigenvec1 ( 3, 8) = 0.0
!      rt_eigenvec1 ( 4, 8) = 0.0
!      rt_eigenvec1 ( 5, 8) = 0.0
      rt_eigenvec1 ( 6, 8) = 1.0
!      rt_eigenvec1 ( 7, 8) = 0.0
!      rt_eigenvec1 ( 8, 8) = 0.0

! ----------

!      lf_eigenvec1 ( 8, 1) = 0.0
!      lf_eigenvec1 ( 8, 2) = 0.0
!      lf_eigenvec1 ( 8, 3) = 0.0
!      lf_eigenvec1 ( 8, 4) = 0.0
!      lf_eigenvec1 ( 8, 5) = 0.0
      lf_eigenvec1 ( 8, 6) = 1.0
!      lf_eigenvec1 ( 8, 7) = 0.0
!      lf_eigenvec1 ( 8, 8) = 0.0

! ----------------------------------------------------------------------

! For the multiple species, fix up the eigenvectors that carry the
! contributions from the species fractions. The zero components of
! the eigenvectors have already been previously assigned.

#if ( NFLUID >= 1)

      DO j = NFIELD + 2, NFIELD + NFLUID + 1

      rt_eigenvec1 ( j, j) = 1.0
      lf_eigenvec1 ( j, j) = 1.0

      END DO

#endif

! ----------------------------------------------------------------------

! At this point the eigenvectors in the primitive variables are made
! and if only that much is desired, we can return at this point.

      IF ( iprimitive .EQ. 1) THEN

! ----------------------------------------------------------------------

      rt_eigenvec ( :, :) = rt_eigenvec1 ( :, :)
      lf_eigenvec ( :, :) = lf_eigenvec1 ( :, :)

! ----------------------------------------------------------------------

      ELSE  ! End of building eigenvectors in primitive variables.

! ----------------------------------------------------------------------

      dudw ( :, :) = 0.0
      dwdu ( :, :) = 0.0

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. First for
! the right eigenvectors. Don't make those matrix elements that are zero.

! Here we have "dwdu * dudw = I" and also "dudw * dwdu = I" where I is the
! identity matrix. Consequently, if "lf_eigenvec1 * rt_eigenvec1 = I" for
! the eigenvectors in primitive form then we also have
! "(lf_eigenvec1 * dwdu) * (dudw * rt_eigenvec1) = I".

! Derivatives of U(1) w.r.t. V(1:8)
      dudw ( 1, 1) = 1.0
!      dudw ( 1, 2) = 0.0
!      dudw ( 1, 3) = 0.0
!      dudw ( 1, 4) = 0.0
!      dudw ( 1, 5) = 0.0
!      dudw ( 1, 6) = 0.0
!      dudw ( 1, 7) = 0.0
!      dudw ( 1, 8) = 0.0


! Derivatives of U(2) w.r.t. V(1:8)
      dudw ( 2, 1) = vlxbar
      dudw ( 2, 2) = rhobar
!      dudw ( 2, 3) = 0.0
!      dudw ( 2, 4) = 0.0
!      dudw ( 2, 5) = 0.0
!      dudw ( 2, 6) = 0.0
!      dudw ( 2, 7) = 0.0
!      dudw ( 2, 8) = 0.0


      dudw ( 3, 1) = vlybar
!      dudw ( 3, 2) = 0.0
      dudw ( 3, 3) = rhobar
!      dudw ( 3, 4) = 0.0
!      dudw ( 3, 5) = 0.0
!      dudw ( 3, 6) = 0.0
!      dudw ( 3, 7) = 0.0
!      dudw ( 3, 8) = 0.0


      dudw ( 4, 1) = vlzbar
!      dudw ( 4, 2) = 0.0
!      dudw ( 4, 3) = 0.0
      dudw ( 4, 4) = rhobar
!      dudw ( 4, 5) = 0.0
!      dudw ( 4, 6) = 0.0
!      dudw ( 4, 7) = 0.0
!      dudw ( 4, 8) = 0.0


      dudw ( 5, 1) = 0.5 * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             + xbar + deng_drho
      dudw ( 5, 2) = rhobar * vlxbar
      dudw ( 5, 3) = rhobar * vlybar
      dudw ( 5, 4) = rhobar * vlzbar
      dudw ( 5, 5) = deng_dprs
      dudw ( 5, 6) = bfxbar / fourpi
      dudw ( 5, 7) = bfybar / fourpi
      dudw ( 5, 8) = bfzbar / fourpi


!      dudw ( 6, 1) = 0.0
!      dudw ( 6, 2) = 0.0
!      dudw ( 6, 3) = 0.0
!      dudw ( 6, 4) = 0.0
!      dudw ( 6, 5) = 0.0
      dudw ( 6, 6) = 1.0
!      dudw ( 6, 7) = 0.0
!      dudw ( 6, 8) = 0.0


!      dudw ( 7, 1) = 0.0
!      dudw ( 7, 2) = 0.0
!      dudw ( 7, 3) = 0.0
!      dudw ( 7, 4) = 0.0
!      dudw ( 7, 5) = 0.0
!      dudw ( 7, 6) = 0.0
      dudw ( 7, 7) = 1.0
!      dudw ( 7, 8) = 0.0


!      dudw ( 8, 1) = 0.0
!      dudw ( 8, 2) = 0.0
!      dudw ( 8, 3) = 0.0
!      dudw ( 8, 4) = 0.0
!      dudw ( 8, 5) = 0.0
!      dudw ( 8, 6) = 0.0
!      dudw ( 8, 7) = 0.0
      dudw ( 8, 8) = 1.0


#if ( NFLUID >= 1)
      DO k = 9, NFIELD + NFLUID + 1

      dudw ( k, 1) = speciebar ( k - 8)
      dudw ( k, k) = rhobar

      END DO
#endif

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. Then for
! the left eigenvectors. Don't make those matrix elements that are zero.

      temparr_01 = 1.0 / rhobar


! Derivatives of V(1:8) w.r.t. U(1)
      dwdu ( 1, 1) = 1.0
      dwdu ( 2, 1) = - temparr_01 * vlxbar
      dwdu ( 3, 1) = - temparr_01 * vlybar
      dwdu ( 4, 1) = - temparr_01 * vlzbar
      dwdu ( 5, 1) = dprs_drho + dprs_deng * ( 0.5 * ( vlxbar**2
     1                                 + vlybar**2 + vlzbar**2) - xbar)
!      dwdu ( 6, 1) = 0.0
!      dwdu ( 7, 1) = 0.0
!      dwdu ( 8, 1) = 0.0


! Derivatives of V(1:8) w.r.t. U(2)
!      dwdu ( 1, 2) = 0.0
      dwdu ( 2, 2) = temparr_01
!      dwdu ( 3, 2) = 0.0
!      dwdu ( 4, 2) = 0.0
      dwdu ( 5, 2) = - dprs_deng * vlxbar
!      dwdu ( 6, 2) = 0.0
!      dwdu ( 7, 2) = 0.0
!      dwdu ( 8, 2) = 0.0


!      dwdu ( 1, 3) = 0.0
!      dwdu ( 2, 3) = 0.0
      dwdu ( 3, 3) = temparr_01
!      dwdu ( 4, 3) = 0.0
      dwdu ( 5, 3) = - dprs_deng * vlybar
!      dwdu ( 6, 3) = 0.0
!      dwdu ( 7, 3) = 0.0
!      dwdu ( 8, 3) = 0.0


!      dwdu ( 1, 4) = 0.0
!      dwdu ( 2, 4) = 0.0
!      dwdu ( 3, 4) = 0.0
      dwdu ( 4, 4) = temparr_01
      dwdu ( 5, 4) = - dprs_deng * vlzbar
!      dwdu ( 6, 4) = 0.0
!      dwdu ( 7, 4) = 0.0
!      dwdu ( 8, 4) = 0.0


!      dwdu ( 1, 5) = 0.0
!      dwdu ( 2, 5) = 0.0
!      dwdu ( 3, 5) = 0.0
!      dwdu ( 4, 5) = 0.0
      dwdu ( 5, 5) = dprs_deng
!      dwdu ( 6, 5) = 0.0
!      dwdu ( 7, 5) = 0.0
!      dwdu ( 8, 5) = 0.0


!      dwdu ( 1, 6) = 0.0
!      dwdu ( 2, 6) = 0.0
!      dwdu ( 3, 6) = 0.0
!      dwdu ( 4, 6) = 0.0
      dwdu ( 5, 6) = - dprs_deng * bfxbar / fourpi
      dwdu ( 6, 6) = 1.0
!      dwdu ( 7, 6) = 0.0
!      dwdu ( 8, 6) = 0.0


!      dwdu ( 1, 7) = 0.0
!      dwdu ( 2, 7) = 0.0
!      dwdu ( 3, 7) = 0.0
!      dwdu ( 4, 7) = 0.0
      dwdu ( 5, 7) = - dprs_deng * bfybar / fourpi
!      dwdu ( 6, 7) = 0.0
      dwdu ( 7, 7) = 1.0
!      dwdu ( 8, 7) = 0.0


!      dwdu ( 1, 8) = 0.0
!      dwdu ( 2, 8) = 0.0
!      dwdu ( 3, 8) = 0.0
!      dwdu ( 4, 8) = 0.0
      dwdu ( 5, 8) = - dprs_deng * bfzbar / fourpi
!      dwdu ( 6, 8) = 0.0
!      dwdu ( 7, 8) = 0.0
      dwdu ( 8, 8) = 1.0


#if ( NFLUID >= 1)
      DO k = 9, NFIELD + NFLUID + 1

      dwdu ( k, 1) = - temparr_01 * speciebar ( k - 8)
      dwdu ( k, k) = temparr_01

      END DO
#endif

! ----------------------------------------------------------------------

! Do the matrix multiplications. Then save the eigenvalues. This
! gives eigenvectors in the conserved variables.
! Don't multiply those matrix elements that are zero.

      rt_eigenvec ( :, :) = 0.0
      lf_eigenvec ( :, :) = 0.0

! Loop over k^th right eigenvector or k^th left eigenvector.
      DO k = 1, NFIELD + NFLUID + 1 

      IF ( ( k == 2) .OR. ( k == 4) .OR. ( k == 6) ) THEN

      DO j = 1, NFIELD + NFLUID + 1

      DO l = 1, NFIELD + NFLUID + 1

      rt_eigenvec ( j, k) = rt_eigenvec ( j, k)
     1                    + dudw ( j, l) * rt_eigenvec1 ( l, k)

      lf_eigenvec ( k, j) = lf_eigenvec ( k, j)
     1                    + lf_eigenvec1 ( k, l) * dwdu ( l, j)

      END DO

      END DO

      END IF  ! End of evaluating just for the entropy and Alfven waves.

      END DO

! ----------------------------------------------------------------------

      END IF  ! End of building eigenvectors in conserved variables.

! ----------------------------------------------------------------------
! The permutation is reversed in the eigenvector components, velocities
! and magnetic fields when leaving the subroutine.

      IF ( idirection == 2) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_03
      vlybar = temparr_01
      vlzbar = temparr_02

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_03
      bfybar = temparr_01
      bfzbar = temparr_02

! -----

      DO k = 1, NFIELD + NFLUID + 1

! -----

      temparr_01 = rt_eigenvec ( 2, k)
      temparr_02 = rt_eigenvec ( 3, k)
      temparr_03 = rt_eigenvec ( 4, k)

      rt_eigenvec ( 2, k) = temparr_03
      rt_eigenvec ( 3, k) = temparr_01
      rt_eigenvec ( 4, k) = temparr_02

! -----

      temparr_01 = rt_eigenvec ( 6, k)
      temparr_02 = rt_eigenvec ( 7, k)
      temparr_03 = rt_eigenvec ( 8, k)

      rt_eigenvec ( 6, k) = temparr_03
      rt_eigenvec ( 7, k) = temparr_01
      rt_eigenvec ( 8, k) = temparr_02

! -----

      temparr_01 = lf_eigenvec ( k, 2)
      temparr_02 = lf_eigenvec ( k, 3)
      temparr_03 = lf_eigenvec ( k, 4)

      lf_eigenvec ( k, 2) = temparr_03
      lf_eigenvec ( k, 3) = temparr_01
      lf_eigenvec ( k, 4) = temparr_02

! -----

      temparr_01 = lf_eigenvec ( k, 6)
      temparr_02 = lf_eigenvec ( k, 7)
      temparr_03 = lf_eigenvec ( k, 8)

      lf_eigenvec ( k, 6) = temparr_03
      lf_eigenvec ( k, 7) = temparr_01
      lf_eigenvec ( k, 8) = temparr_02

! -----

      END DO

! -----

      END IF

! ----------

      IF ( idirection == 3) THEN

! -----

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_02
      vlybar = temparr_03
      vlzbar = temparr_01

! -----

      temparr_01 = bfxbar
      temparr_02 = bfybar
      temparr_03 = bfzbar

      bfxbar = temparr_02
      bfybar = temparr_03
      bfzbar = temparr_01

! -----

      DO k = 1, NFIELD + NFLUID + 1

! -----

      temparr_01 = rt_eigenvec ( 2, k)
      temparr_02 = rt_eigenvec ( 3, k)
      temparr_03 = rt_eigenvec ( 4, k)

      rt_eigenvec ( 2, k) = temparr_02
      rt_eigenvec ( 3, k) = temparr_03
      rt_eigenvec ( 4, k) = temparr_01

! -----

      temparr_01 = rt_eigenvec ( 6, k)
      temparr_02 = rt_eigenvec ( 7, k)
      temparr_03 = rt_eigenvec ( 8, k)

      rt_eigenvec ( 6, k) = temparr_02
      rt_eigenvec ( 7, k) = temparr_03
      rt_eigenvec ( 8, k) = temparr_01

! -----

      temparr_01 = lf_eigenvec ( k, 2)
      temparr_02 = lf_eigenvec ( k, 3)
      temparr_03 = lf_eigenvec ( k, 4)

      lf_eigenvec ( k, 2) = temparr_02
      lf_eigenvec ( k, 3) = temparr_03
      lf_eigenvec ( k, 4) = temparr_01

! -----

      temparr_01 = lf_eigenvec ( k, 6)
      temparr_02 = lf_eigenvec ( k, 7)
      temparr_03 = lf_eigenvec ( k, 8)

      lf_eigenvec ( k, 6) = temparr_02
      lf_eigenvec ( k, 7) = temparr_03
      lf_eigenvec ( k, 8) = temparr_01

! -----

      END DO

! -----

      END IF

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE EIGENVECMHD_PTWISE_LINDEG

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)
      PROGRAM TEST_EIGENVECMHD_PTWISE_8

! ----------------------------------------------------------------------
! Variables used for testing "EIGENVECMHD_PTWISE_8".

      INTEGER :: nriembegin, nriemend, iprimitive, idirection,
     1           i, j, k, m

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor, pi, gasconst

      REAL :: rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     bfxbar, bfybar, bfzbar, xbar,

     1     ob_Rho, ob_4pi, gm1, gm2, Prs_Tot_Star, Eng_Tot_Star,
     1     V_Dot_B_Star

      REAL, DIMENSION ( 0: NFLUID) :: speciebar

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: eigenval,

     1      dPrs_Tot_dU

      REAL, DIMENSION ( NFIELD + NFLUID + 1, NFIELD + NFLUID + 1) ::
     1     lf_eigenvec, rt_eigenvec,

     1     Itest, Abar, Bbar, AbarEig, BbarEig

! ----------------------------------------------------------------------
! Variables used for comparing with "EIGENVECMHD".

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar_1d, prsbar_1d, vlxbar_1d, vlybar_1d, vlzbar_1d,
     1     bfxbar_1d, bfybar_1d, bfzbar_1d, xbar_1d

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      speciebar_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID) ::
     1     eigenval_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID) ::
     1     lf_eigenvec_1d, rt_eigenvec_1d

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "test_eigenvecmhd.out")

      iprimitive = 0
      idirection = 1

      nriembegin = 1
      nriemend = 10

      gamma = 1.4
      gasconst = 1.0

      cp_specie ( :) = 2.0
      molewt_specie ( :) = 3.0
      heat_of_formation ( :) = 4.0

      smallnum = 1.0e-12
      rhofloor = 1.0e-8
      prsfloor = 1.0e-12

      pi = 4.0 * ATAN ( 1.0)
      ob_4pi = 1.0 / (4.0 * pi)

      debug_stop = 0

      rhobar = 0.1
      ob_Rho = 1.0 / rhobar
      prsbar = 2.2

      vlxbar = 1.0
      vlybar = 2.0
      vlzbar = -3.0
      
      bfxbar = -2.0
      bfybar = 1.0
      bfzbar = -1.5
      xbar = 0.0

      Prs_Tot_Star = prsbar + ( bfxbar**2 + bfybar**2 + bfzbar**2)
     1                        / ( 8.0 * pi)

      Eng_Tot_Star = prsbar / ( gamma - 1.0)
     1             + 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             + ( bfxbar**2 + bfybar**2 + bfzbar**2) / ( 8.0 * pi)

      V_Dot_B_Star = vlxbar * bfxbar + vlybar * bfybar + vlzbar * bfzbar

#if ( NFLUID >= 1)
      speciebar ( :) = 1.0 / NFLUID
#endif

! ----------------------------------------------------------------------

      rhobar_1d ( nriembegin: nriemend) = rhobar
      prsbar_1d ( nriembegin: nriemend) = prsbar

      vlxbar_1d ( nriembegin: nriemend) = vlxbar
      vlybar_1d ( nriembegin: nriemend) = vlybar
      vlzbar_1d ( nriembegin: nriemend) = vlzbar

      bfxbar_1d ( nriembegin: nriemend) = bfxbar
      bfybar_1d ( nriembegin: nriemend) = bfybar
      bfzbar_1d ( nriembegin: nriemend) = bfzbar
      xbar_1d ( nriembegin: nriemend) = xbar

#if ( NFLUID >= 1)
      DO i = nriembegin, nriemend
      speciebar_1d ( :, i) = speciebar ( :)
      END DO
#endif

! ----------------------------------------------------------------------

      CALL EIGENVECMHD_PTWISE_8

     1             ( iprimitive, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! ------------------------------

      WRITE ( 6, *)"eigenval = ", eigenval
      WRITE ( 6, *)" "

      Itest ( :, :) = 0.0

      DO i = 1, NFIELD + NFLUID + 1
      DO j = 1, NFIELD + NFLUID + 1
      DO k = 1, NFIELD + NFLUID + 1

      Itest ( i, j) = Itest ( i, j)
     1              + lf_eigenvec ( i, k) * rt_eigenvec ( k, j)

      END DO
      END DO
      END DO

      DO i = 1, NFIELD + NFLUID + 1
      WRITE ( 6, *)"i = ", i
      WRITE ( 6, *)( Itest ( i, j), j = 1, NFIELD + NFLUID + 1)
      WRITE ( 6, *)" "
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"-----------------------------------------------"
      WRITE ( 6, *)"Done with orthonormality test for eigenvectors"
      WRITE ( 6, *)"-----------------------------------------------"

! ----------------------------------------------------------------------

      CALL EIGENVECMHD

     1             ( nriembegin, nriemend, iprimitive,

     1               gamma, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,
     1               cp_specie, molewt_specie, heat_of_formation,

     1               rhobar_1d, prsbar_1d,
     1               vlxbar_1d, vlybar_1d, vlzbar_1d,
     1               bfxbar_1d, bfybar_1d, bfzbar_1d, xbar_1d,
     1               speciebar_1d,

     1               eigenval_1d,
     1               lf_eigenvec_1d, rt_eigenvec_1d )

! ------------------------------

      i = 1 ! Only compare the first zone for correctness; they're all the same

      DO k = 1, 7  ! Index "k" runs over waves

      DO j = 1, 7  ! Index "j" runs over eigenvector components.

      m = j
      IF ( m > 5) m = j + 1

      IF ( ABS ( rt_eigenvec_1d ( i, j, k) - rt_eigenvec ( m, k) ) >
     1     smallnum )
     1   WRITE ( 6, *)"Disagreement rt : j, k = ", j, k," values = ",
     1   rt_eigenvec_1d ( i, j, k), rt_eigenvec ( m, k)

      IF ( ABS ( lf_eigenvec_1d ( i, k, j) - lf_eigenvec ( k, m) ) >
     1     smallnum )
     1   WRITE ( 6, *)"Disagreement lf : j, k = ", j, k," values = ",
     1   lf_eigenvec_1d ( i, k, j), lf_eigenvec ( k, m)

      END DO  ! Index "j" runs over eigenvector components.

      END DO  ! Index "k" runs over waves

      WRITE ( 6, *)" "
      WRITE ( 6, *)"-----------------------------------------------"
      WRITE ( 6, *)"Done comparing old and new eigenvector routines"
      WRITE ( 6, *)"-----------------------------------------------"

! ----------------------------------------------------------------------

! Now test whether the characteristic matrices are exactly reproduced.

! ----------------------------------------------------------------------

! Zero out characteristic matrices.

      Abar ( :, :) = 0.0
      Bbar ( :, :) = 0.0

      AbarEig ( :, :) = 0.0
      BbarEig ( :, :) = 0.0

! ------------------------------

! Evaluate the derivatives of the total pressure.

      dPrs_Tot_dU ( :) = 0.0

      gm1 = gamma - 1.0
      gm2 = gamma - 2.0

      dPrs_Tot_dU ( 1) = 0.5 * gm1
     1                 * (vlxbar**2 + vlybar**2 + vlzbar**2)

      dPrs_Tot_dU ( 2) = - gm1 * vlxbar
      dPrs_Tot_dU ( 3) = - gm1 * vlybar
      dPrs_Tot_dU ( 4) = - gm1 * vlzbar

      dPrs_Tot_dU ( 5) = gm1

      dPrs_Tot_dU ( 6) = - ob_4pi * gm2 * bfxbar
      dPrs_Tot_dU ( 7) = - ob_4pi * gm2 * bfybar
      dPrs_Tot_dU ( 8) = - ob_4pi * gm2 * bfzbar

! ------------------------------

! Build the characteristic matrix "Abar".

      Abar ( 1, 2) = 1.0

      Abar ( 2, 1) = - vlxbar**2 + dPrs_Tot_dU ( 1)
      Abar ( 2, 2) = 2.0 * vlxbar + dPrs_Tot_dU ( 2)
      Abar ( 2, 3) = dPrs_Tot_dU ( 3)
      Abar ( 2, 4) = dPrs_Tot_dU ( 4)
      Abar ( 2, 5) = dPrs_Tot_dU ( 5)
      Abar ( 2, 6) = - 2.0 * ob_4pi * bfxbar + dPrs_Tot_dU ( 6)
      Abar ( 2, 7) = dPrs_Tot_dU ( 7)
      Abar ( 2, 8) = dPrs_Tot_dU ( 8)

      Abar ( 3, 1) = - vlxbar * vlybar
      Abar ( 3, 2) = vlybar
      Abar ( 3, 3) = vlxbar
      Abar ( 3, 6) = - ob_4pi * bfybar
      Abar ( 3, 7) = - ob_4pi * bfxbar

      Abar ( 4, 1) = - vlxbar * vlzbar
      Abar ( 4, 2) = vlzbar
      Abar ( 4, 4) = vlxbar
      Abar ( 4, 6) = - ob_4pi * bfzbar
      Abar ( 4, 8) = - ob_4pi * bfxbar

! -----

      Abar ( 5, 1) = -(Eng_Tot_Star + Prs_Tot_Star) * ob_Rho * vlxbar
     1             + dPrs_Tot_dU ( 1) * vlxbar
     1             + ob_4pi * ob_Rho * bfxbar * V_Dot_B_Star
      Abar ( 5, 2) = (Eng_Tot_Star + Prs_Tot_Star) * ob_Rho
     1             + dPrs_Tot_dU ( 2) * vlxbar
     1             - ob_4pi * ob_Rho * bfxbar**2
      Abar ( 5, 3) = dPrs_Tot_dU ( 3) * vlxbar
     1             - ob_4pi * ob_Rho * bfxbar * bfybar
      Abar ( 5, 4) = dPrs_Tot_dU ( 4) * vlxbar
     1             - ob_4pi * ob_Rho * bfxbar * bfzbar
      Abar ( 5, 5) = ( 1.0 + dPrs_Tot_dU ( 5) ) * vlxbar
      Abar ( 5, 6) = dPrs_Tot_dU ( 6) * vlxbar
     1             - ob_4pi * vlxbar * bfxbar
     1             - ob_4pi * V_Dot_B_Star
      Abar ( 5, 7) = dPrs_Tot_dU ( 7) * vlxbar
     1             - ob_4pi * vlybar * bfxbar
      Abar ( 5, 8) = dPrs_Tot_dU ( 8) * vlxbar
     1             - ob_4pi * vlzbar * bfxbar

! -----

      Abar ( 7, 1) = - ob_Rho
     1             * ( vlxbar * bfybar - vlybar * bfxbar) 
      Abar ( 7, 2) = ob_Rho * bfybar
      Abar ( 7, 3) = - ob_Rho * bfxbar
      Abar ( 7, 6) = - vlybar
      Abar ( 7, 7) = vlxbar

      Abar ( 8, 1) = - ob_Rho 
     1             * ( vlxbar * bfzbar - vlzbar * bfxbar)
      Abar ( 8, 2) = ob_Rho * bfzbar
      Abar ( 8, 4) = - ob_Rho * bfxbar
      Abar ( 8, 6) = - vlzbar
      Abar ( 8, 8) = vlxbar

! ----------

! If we acknowledge that "bxbar" cannot change, then the 6th term's
! contributions will be zero. So we have a change that is more suited
! for comparison with eigenvectors.

      Abar ( :, 6) = 0.0

! ------------------------------

! Build the characteristic matrix "Bbar".

      Bbar ( 1, 3) = 1.0

      Bbar ( 2, 1) = - vlxbar * vlybar
      Bbar ( 2, 2) = vlybar
      Bbar ( 2, 3) = vlxbar
      Bbar ( 2, 6) = - ob_4pi * bfybar
      Bbar ( 2, 7) = - ob_4pi * bfxbar

      Bbar ( 3, 1) = - vlybar**2 + dPrs_Tot_dU ( 1)
      Bbar ( 3, 2) = dPrs_Tot_dU ( 2)
      Bbar ( 3, 3) = 2.0 * vlybar + dPrs_Tot_dU ( 3)
      Bbar ( 3, 4) = dPrs_Tot_dU ( 4)
      Bbar ( 3, 5) = dPrs_Tot_dU ( 5)
      Bbar ( 3, 6) = dPrs_Tot_dU ( 6)
      Bbar ( 3, 7) = - 2.0 * ob_4pi * bfybar + dPrs_Tot_dU ( 7)
      Bbar ( 3, 8) = dPrs_Tot_dU ( 8)

      Bbar ( 4, 1) = - vlybar * vlzbar
      Bbar ( 4, 3) = vlzbar
      Bbar ( 4, 4) = vlybar
      Bbar ( 4, 7) = - ob_4pi * bfzbar
      Bbar ( 4, 8) = - ob_4pi * bfybar

! -----

      Bbar ( 5, 1) = -(Eng_Tot_Star + Prs_Tot_Star) * ob_Rho * vlybar
     1             + dPrs_Tot_dU ( 1) * vlybar
     1             + ob_4pi * ob_Rho * bfybar * V_Dot_B_Star
      Bbar ( 5, 2) = dPrs_Tot_dU ( 2) * vlybar
     1             - ob_4pi * ob_Rho * bfxbar * bfybar
      Bbar ( 5, 3) = (Eng_Tot_Star + Prs_Tot_Star) * ob_Rho
     1             + dPrs_Tot_dU ( 3) * vlybar
     1             - ob_4pi * ob_Rho * bfybar**2
      Bbar ( 5, 4) = dPrs_Tot_dU ( 4) * vlybar
     1             - ob_4pi * ob_Rho * bfybar * bfzbar
      Bbar ( 5, 5) = ( 1.0 + dPrs_Tot_dU ( 5) ) * vlybar
      Bbar ( 5, 6) = dPrs_Tot_dU ( 6) * vlybar
     1             - ob_4pi * vlxbar * bfybar
      Bbar ( 5, 7) = dPrs_Tot_dU ( 7) * vlybar
     1             - ob_4pi * vlybar * bfybar
     1             - ob_4pi * V_Dot_B_Star
      Bbar ( 5, 8) = dPrs_Tot_dU ( 8) * vlybar
     1             - ob_4pi * vlzbar * bfybar

! -----

      Bbar ( 6, 1) = - ob_Rho
     1             * ( vlybar * bfxbar - vlxbar * bfybar)
      Bbar ( 6, 2) = - ob_Rho * bfybar
      Bbar ( 6, 3) = ob_Rho * bfxbar
      Bbar ( 6, 6) = vlybar
      Bbar ( 6, 7) = - vlxbar

      Bbar ( 8, 1) = - ob_Rho
     1             * ( vlybar * bfzbar - vlzbar * bfybar)
      Bbar ( 8, 3) = ob_Rho * bfzbar
      Bbar ( 8, 4) = - ob_Rho * bfybar
      Bbar ( 8, 7) = - vlzbar
      Bbar ( 8, 8) = vlybar

! ----------

! If we acknowledge that "bybar" cannot change, then the 7th term's
! contributions will be zero. So we have a change that is more suited
! for comparison with eigenvectors.

      Bbar ( :, 7) = 0.0

! ----------------------------------------------------------------------
! Build alternative version of "AbarEig" based on eigenvectors and test
! against "Abar".

      iprimitive = 0
      idirection = 1

      CALL EIGENVECMHD_PTWISE_8

     1             ( iprimitive, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! --------------------

      AbarEig ( :, :) = 0.0

      DO k = 1, NFIELD + NFLUID + 1

        DO i = 1, NFIELD + NFLUID + 1

          DO j = 1, NFIELD + NFLUID + 1

          AbarEig ( i, j) = AbarEig ( i, j)
     1+ rt_eigenvec ( i, k) * eigenval ( k) * lf_eigenvec ( k, j)

          END DO

        END DO

      END DO

! --------------------
! B_x is prevented from varying. Except for the "AbarEig ( i, 6)" component,
! the other components should agree.

      DO i = 1, NFIELD + NFLUID + 1

        DO j = 1, NFIELD + NFLUID + 1

        IF ( ABS ( AbarEig ( i, j) - Abar ( i, j) ) > smallnum)
     1  WRITE ( 6, *)"Problem : AbarEig, Abar ( i, j) = ", i, j

        END DO

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"-----------------------------------------------"
      WRITE ( 6, *)"Done comparing AbarEig with Abar"
      WRITE ( 6, *)"-----------------------------------------------"

! ----------------------------------------------------------------------
! Build alternative version of "BbarEig" based on eigenvectors and test
! against "Bbar".

      iprimitive = 0
      idirection = 2

      CALL EIGENVECMHD_PTWISE_8

     1             ( iprimitive, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! --------------------

      BbarEig ( :, :) = 0.0

      DO k = 1, NFIELD + NFLUID + 1

        DO i = 1, NFIELD + NFLUID + 1

          DO j = 1, NFIELD + NFLUID + 1

          BbarEig ( i, j) = BbarEig ( i, j)
     1+ rt_eigenvec ( i, k) * eigenval ( k) * lf_eigenvec ( k, j)

          END DO

        END DO

      END DO

! --------------------
! B_y is prevented from varying. Except for the "BbarEig ( i, 7)" component,
! the other components should agree.

      DO i = 1, NFIELD + NFLUID + 1

        DO j = 1, NFIELD + NFLUID + 1

        IF ( ABS ( BbarEig ( i, j) - Bbar ( i, j) ) > smallnum)
     1  WRITE ( 6, *)"Problem : BbarEig, Bbar ( i, j) = ", i, j

        END DO

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"-----------------------------------------------"
      WRITE ( 6, *)"Done comparing BbarEig with Bbar"
      WRITE ( 6, *)"-----------------------------------------------"

! ----------------------------------------------------------------------





! ----------------------------------------------------------------------

      STOP

      END PROGRAM TEST_EIGENVECMHD_PTWISE_8
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------















































