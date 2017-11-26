! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

      SUBROUTINE MHD_RIEM_HLLEM_PTWISE

     1     ( gamma, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vlx_L, vlx_R,
     1        vly_L, vly_R,
     1        vlz_L, vlz_R,
     1        bfxre,
     1        bfy_L, bfy_R,
     1        bfz_L, bfz_R,
     1        specie_L, specie_R,

     1        gamma_eos_L, gamma_eos_R,
     1        gamma_soundspeed_L, gamma_soundspeed_R,
     1        mean_mol_wt_L, mean_mol_wt_R,

     1        blend_HLL,

     1        S_L, S_R, 
     1        F_L, F_R, U_HLL, F_HLL,
     1        U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM, 

     1        flux_re, entropy_flux_re)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds a 1D HLLEM Riemann solver for MHD based on the work of
! Dumbser & Balsara (2016) JCP and Balsara et al. (2016) JCP.

! Exact preservation of stationary intermediate waves is achieved.
! This includes contact discontinuities and also Alfven waves.

! ----------------------------------------------------------------------

! Implements the Dumbser and Balsara (2016) JCP paper (DB16) and the
! Balsara et al (2016) JCP paper (Betal16).

! Inputs :-

! "iHLLEM_choice = 1" gives original Dumbser-Balsara HLLEM; it has
! slighly less dissipation for fast-moving waves. "iHLLEM_choice = 2" gives
! Balsara et al HLLEM; it is more stable for strong shocks.
! Both will capture stationary intermediate waves. "iHLLEM_choice = 0" gives
! HLL.

! split_prsflux_re : Splits the momentum flux into a pressure part and an
!                advected part if set to 1. Useful in certain geometries.

! gasconst, cp_specie, molewt_specie, heat_of_formation : Thermodynamic inputs.
! rho_L, rho_R : left and right state densities, all other flow variables are
!                similarly also inputs.
! gamma_eos_L, gamma_eos_R : left and right state ratio of specific heats.
!                Similarly for other variables.
! blend_HLL : Gives the amount of blend_HLLing, i.e. reduction in
!                "dU_HLLEM, dF_HLLEM" required to stabilize strong shocks.
!                "blend_HLL = 0" means we get full HLLEM flux;
!                "blend_HLL = 1" means we get HLL flux.

! Outputs :-

! S_L, S_R : Extremal speeds of the HLL Riemann solver.
! F_L, F_R : Left and right fluxes, used as input for HLL RS.
! U_HLL, F_HLL : State and flux coming out of HLL RS.
! U_HLLEM, F_HLLEM : Final HLLEM state and flux that overlies the time axis
!            in the Riemann fan.
! dU_HLLEM, dF_HLLEM : These are the terms that have tobe suitably added
!            to the HLL state and flux in order to get the HLLEM state and flux.
! dU_HLLEM_Min : "dU_HLLEM" has a compression built in. But we cannot use that
!            to build the linear profile in the multid RS (because it would be
!            a double counting of the slope. "dU_HLLEM_Min" removes the
!            compression which would have doubled the steepening of the slopes.
! flux_re : Seven component flux, with pressure split, specific to certain
!            my MHD code.
! entropy_flux_re : Flux for entropy variable.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER use_approx_eos, split_prsflux_re, iHLLEM_choice

      REAL gamma, smallnum, prsfloor, rhofloor, pi, gasconst,

     1     cp_specie ( 0 : NFLUID), molewt_specie ( 0 : NFLUID),
     1     heat_of_formation ( 0 : NFLUID),

     1     rho_L, rho_R, prs_L, prs_R,
     1     vlx_L, vlx_R, vly_L, vly_R, vlz_L, vlz_R,

     1     bfxre, bfy_L, bfy_R, bfz_L, bfz_R,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     blend_HLL, S_L, S_R, entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: flux_re,
     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

! ----------------------------------------------------------------------

! Variables past this point are not used in the subroutine interface.

      INTEGER i, j, iprim, idirection, debug_stop

      REAL fourpi, ob_4pi, ob_8pi,
     1     rhofloor_plus_delta, prsfloor_plus_delta,

     1     eng_L, eng_R, v2_L, v2_R, b2_L, b2_R, vdb_L, vdb_R,
     1     total_enthalpy_L, total_enthalpy_R,

     1     entropy_L, entropy_R, entropy_flux_L, entropy_flux_R,

     1     gamma_eos, gamma_soundspeed, mean_mol_wt,
     1     temp_01, temp_02, temp_03,
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, total_enthalpybar,
     1     bfxbar, bfybar, bfzbar, xbar, eng, temperature,

     1     ob_4pi_rho, va_x, c_s, v_a, v_f, v_s

      REAL, DIMENSION ( 0: NFLUID) :: speciebar

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: U_L, U_R,
     1     Del_U, eigenval, delta, eigenwt

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NFIELD + 1 + NFLUID) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 * ob_4pi

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)"MHD_RIEM_HLLEM_PTWISE : NFIELD .NE. 7", NFIELD, 7
      STOP
      END IF

! ----------
! Because floor values of density and pressure have already been taken,
! we set a bottom value that is a little above that floor.

      rhofloor_plus_delta = 1.000001 * rhofloor
      prsfloor_plus_delta = 1.000001 * prsfloor

      rho_L = AMAX1 ( rho_L, rhofloor)
      rho_R = AMAX1 ( rho_R, rhofloor)

      prs_L = AMAX1 ( prs_L, prsfloor)
      prs_R = AMAX1 ( prs_R, prsfloor)

! ----------
! Relying on approximate EOS's here to make thermal energy density.

      eng_L = prs_L / ( gamma_eos_L - 1.0)
      eng_R = prs_R / ( gamma_eos_R - 1.0)

! ----------
! Make "total_enthalpy_L" and "total_enthalpy_R". These are
! useful in making Roe-averages.

      v2_L = vlx_L**2 + vly_L**2 + vlz_L**2
      b2_L = bfxre**2 + bfy_L**2 + bfz_L**2
      vdb_L = vlx_L * bfxre + vly_L * bfy_L + vlz_L * bfz_L

      total_enthalpy_L = ( eng_L + prs_L
     1    + 0.5 * rho_L * v2_L + ob_4pi * b2_L ) / rho_L

! -----

      v2_R = vlx_R**2 + vly_R**2 + vlz_R**2
      b2_R = bfxre**2 + bfy_R**2 + bfz_R**2
      vdb_R = vlx_R * bfxre + vly_R * bfy_R + vlz_R * bfz_R

      total_enthalpy_R = ( eng_R + prs_R
     1    + 0.5 * rho_R * v2_R + ob_4pi * b2_R ) / rho_R

! ------------------------------
! Make "U_L, U_R".

      U_L ( 1) = rho_L
      U_L ( 2) = rho_L * vlx_L
      U_L ( 3) = rho_L * vly_L
      U_L ( 4) = rho_L * vlz_L
      U_L ( 5) = eng_L + 0.5 * ( rho_L * v2_L + ob_4pi * b2_L )
      U_L ( 6) = bfxre
      U_L ( 7) = bfy_L
      U_L ( 8) = bfz_L
#if ( NFLUID >= 1)
      U_L ( 9: NFIELD + 1 + NFLUID) = rho_L * specie_L ( 1:  NFLUID)
#endif

! ----------

      U_R ( 1) = rho_R
      U_R ( 2) = rho_R * vlx_R
      U_R ( 3) = rho_R * vly_R
      U_R ( 4) = rho_R * vlz_R
      U_R ( 5) = eng_R + 0.5 * ( rho_R * v2_R + ob_4pi * b2_R )
      U_R ( 6) = bfxre
      U_R ( 7) = bfy_R
      U_R ( 8) = bfz_R
#if ( NFLUID >= 1)
      U_R ( 9: NFIELD + 1 + NFLUID) = rho_R * specie_R ( 1:  NFLUID)
#endif

! ------------------------------
! Make "F_L, F_R".

      F_L ( :) = vlx_L * U_L ( :)


      F_L ( 2) = F_L ( 2) + prs_L
     1         + 0.5 * ( - bfxre**2 + bfy_L**2 + bfz_L**2) * ob_4pi

      F_L ( 3) = F_L ( 3) - bfxre * bfy_L * ob_4pi

      F_L ( 4) = F_L ( 4) - bfxre * bfz_L * ob_4pi

      F_L ( 5) = F_L ( 5) + vlx_L * ( prs_L + 0.5 * b2_L * ob_4pi )
     1         - bfxre * vdb_L * ob_4pi

      F_L ( 6) = 0.0

      F_L ( 7) = F_L ( 7) - bfxre * vly_L

      F_L ( 8) = F_L ( 8) - bfxre * vlz_L

! ----------

      F_R ( :) = vlx_R * U_R ( :)


      F_R ( 2) = F_R ( 2) + prs_R
     1         + 0.5 * ( - bfxre**2 + bfy_R**2 + bfz_R**2) * ob_4pi

      F_R ( 3) = F_R ( 3) - bfxre * bfy_R * ob_4pi

      F_R ( 4) = F_R ( 4) - bfxre * bfz_R * ob_4pi

      F_R ( 5) = F_R ( 5) + vlx_R * ( prs_R + 0.5 * b2_R * ob_4pi )
     1         - bfxre * vdb_R * ob_4pi

      F_R ( 6) = 0.0

      F_R ( 7) = F_R ( 7) - bfxre * vly_R

      F_R ( 8) = F_R ( 8) - bfxre * vlz_R

! ------------------------------
! Make variables that are used for obtaining an entropy flux.

      entropy_L = prs_L / rho_L**( gamma_soundspeed_L - 1.0)

      entropy_R = prs_R / rho_R**( gamma_soundspeed_R - 1.0)

      entropy_flux_L = entropy_L * vlx_L

      entropy_flux_R = entropy_R * vlx_R

! ----------------------------------------------------------------------
! Here we use Roe averaging and make only one call to the
! eigenvector subroutine. The "eng, temperature" etc. that we make here
! can even serve as good initial guesses for an exact EOS.

      gamma_eos = 0.5 * ( gamma_eos_L + gamma_eos_R)
      gamma_soundspeed = 0.5 * ( gamma_soundspeed_L
     1                         + gamma_soundspeed_R)
      mean_mol_wt = 0.5 * ( mean_mol_wt_L + mean_mol_wt_R)

! ----------

      temp_01 = SQRT ( rho_R / rho_L)
      temp_02 = 1.0 / ( 1.0 + temp_01)

      rhobar = SQRT ( rho_L * rho_R)

      vlxbar = ( vlx_L + vlx_R * temp_01) * temp_02
      vlybar = ( vly_L + vly_R * temp_01) * temp_02
      vlzbar = ( vlz_L + vlz_R * temp_01) * temp_02

      total_enthalpybar = ( total_enthalpy_L
     1                    + total_enthalpy_R * temp_01) * temp_02

      bfxbar = bfxre
      bfybar = ( bfy_L * temp_01 + bfy_R) * temp_02
      bfzbar = ( bfz_L * temp_01 + bfz_R) * temp_02

      xbar = 0.5 * ( ( bfy_R - bfy_L)**2 + ( bfz_R - bfz_L)**2)
     1     / ( SQRT ( rho_L) + SQRT ( rho_R) )**2

      prsbar = total_enthalpybar
     1       - 0.5 * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1       - ( bfxre**2 + bfybar**2 + bfzbar**2) / ( fourpi * rhobar)
      prsbar = ( gamma_eos - 1.0) * rhobar * prsbar / gamma_eos
      prsbar = prsbar - ( gamma_eos - 2.0) * rhobar * xbar
     1                / ( fourpi * gamma_eos )
      prsbar = AMAX1 ( prsbar, prsfloor)

! ----------
! Use the option of using arithmetic averages if needed.

      IF ( ( rhobar .LE. rhofloor_plus_delta) .OR.
     1     ( prsbar .LE. prsfloor_plus_delta) ) THEN

      rhobar = 0.5 * ( rho_L + rho_R)
      prsbar = 0.5 * ( prs_L + prs_R)

      vlxbar = 0.5 * ( vlx_L + vlx_R)
      vlybar = 0.5 * ( vly_L + vly_R)
      vlzbar = 0.5 * ( vlz_L + vlz_R)

      bfxbar = bfxre
      bfybar = 0.5 * ( bfy_L + bfy_R)
      bfzbar = 0.5 * ( bfz_L + bfz_R)

      xbar = 0.0

      END IF

      eng = prsbar / ( gamma_eos - 1.0)
      temperature = prsbar * mean_mol_wt / ( gasconst * rhobar)

! ----------

#if ( NFLUID >= 1)
      temp_01 = SQRT ( rho_R / rho_L)
      temp_02 = 1.0 / ( 1.0 + temp_01)

      speciebar ( 1: NFLUID) = temp_02 * ( specie_L ( 1: NFLUID)
     1                       + specie_R ( 1: NFLUID) * temp_01)
#endif

! ----------------------------------------------------------------------
! Make extremal speeds, "S_L, S_R", using the Roe averaged state.

      ob_4pi_rho = 1.0 / ( fourpi * rhobar )

      va_x = bfxre**2 * ob_4pi_rho
      va_x = SQRT ( va_x )

      c_s = gamma_soundspeed * prsbar / rhobar
      c_s = SQRT ( c_s )

      v_a = ( bfxre**2 + bfybar**2 + bfzbar**2) * ob_4pi_rho
      v_a = SQRT ( v_a )

      temp_01 = v_a**2 + c_s**2
      temp_02 = temp_01**2 - 4.0 * va_x**2 * c_s**2
      temp_02 = AMAX1 ( temp_02, 0.0)
      temp_02 = SQRT ( temp_02)

      v_f = 0.5 * ( temp_01 + temp_02)
      v_f = SQRT ( v_f )

!      v_s = 0.5 * ( temp_01 - temp_02)
!      v_s = AMAX1 ( v_s, 0.0)
!      v_s = SQRT ( v_s )

! ----------

      S_L = vlxbar - v_f
      S_R = vlxbar + v_f

! ----------------------------------------------------------------------
! Expand the extremal speeds, "S_L, S_R", using "rho_L" etc.

      ob_4pi_rho = 1.0 / ( fourpi * rho_L )

      va_x = bfxre**2 * ob_4pi_rho
      va_x = SQRT ( va_x )

      c_s = gamma_soundspeed_L * prs_L / rho_L
      c_s = SQRT ( c_s )

      v_a = ( bfxre**2 + bfy_L**2 + bfz_L**2) * ob_4pi_rho
      v_a = SQRT ( v_a )

      temp_01 = v_a**2 + c_s**2
      temp_02 = temp_01**2 - 4.0 * va_x**2 * c_s**2
      temp_02 = AMAX1 ( temp_02, 0.0)
      temp_02 = SQRT ( temp_02)

      v_f = 0.5 * ( temp_01 + temp_02)
      v_f = SQRT ( v_f )

!      v_s = 0.5 * ( temp_01 - temp_02)
!      v_s = AMAX1 ( v_s, 0.0)
!      v_s = SQRT ( v_s )

! ----------

      S_L = AMIN1 ( S_L, vlx_L - v_f)
! While it is ok to uncomment the next line, in principle, the left state
! does not contribute to "S_R".
!      S_R = AMAX1 ( S_R, vlx_L + v_f)

! ----------------------------------------------------------------------
! Expand the extremal speeds, "S_L, S_R", using "rho_R" etc.

      ob_4pi_rho = 1.0 / ( fourpi * rho_R )

      va_x = bfxre**2 * ob_4pi_rho
      va_x = SQRT ( va_x )

      c_s = gamma_soundspeed_R * prs_R / rho_R
      c_s = SQRT ( c_s )

      v_a = ( bfxre**2 + bfy_R**2 + bfz_R**2) * ob_4pi_rho
      v_a = SQRT ( v_a )

      temp_01 = v_a**2 + c_s**2
      temp_02 = temp_01**2 - 4.0 * va_x**2 * c_s**2
      temp_02 = AMAX1 ( temp_02, 0.0)
      temp_02 = SQRT ( temp_02)

      v_f = 0.5 * ( temp_01 + temp_02)
      v_f = SQRT ( v_f )

!      v_s = 0.5 * ( temp_01 - temp_02)
!      v_s = AMAX1 ( v_s, 0.0)
!      v_s = SQRT ( v_s )

! ----------

! While it is ok to uncomment the next line, in principle, the right state
! does not contribute to "S_L".
!      S_L = AMIN1 ( S_L, vlx_R - v_f)
      S_R = AMAX1 ( S_R, vlx_R + v_f)

! ----------------------------------------------------------------------
! Now that "S_L, S_R" and "U_L, U_R, F_L, F_R" are built, we
! can build the HLL state and flux "U_HLL, F_HLL". Do similarly for
! entropy flux. The entropy flux does not take HLLEM-based contributions,
! so it is fully evaluated here.

! "U_HLLEM, F_HLLEM", and everything that contributes to them, are also
! set up here for the supersonic cases. For the subsonic cases, we
! do that later.

      IF ( S_L >= 0.0) THEN

        U_HLL ( :) = U_L ( :)
        F_HLL ( :) = F_L ( :)
        entropy_flux_re = entropy_flux_L

        U_HLLEM ( :) = U_L ( :)
        F_HLLEM ( :) = F_L ( :)

        dU_HLLEM ( :) = 0.0
        dU_HLLEM_Min ( :) = 0.0
        dF_HLLEM ( :) = 0.0

        delta ( :) = 0.0
        eigenwt ( :) = 0.0

      ELSE IF ( S_R <= 0.0) THEN

        U_HLL ( :) = U_R ( :)
        F_HLL ( :) = F_R ( :)
        entropy_flux_re = entropy_flux_R

        U_HLLEM ( :) = U_R ( :)
        F_HLLEM ( :) = F_R ( :)

        dU_HLLEM ( :) = 0.0
        dU_HLLEM_Min ( :) = 0.0
        dF_HLLEM ( :) = 0.0

        delta ( :) = 0.0
        eigenwt ( :) = 0.0

      ELSE

        temp_01 = 1.0 / ( S_R - S_L)
        temp_02 = S_R * S_L * temp_01

        U_HLL ( :) = ( S_R * U_R ( :) - S_L * U_L ( :)) * temp_01
     1             - ( F_R ( :) - F_L ( :)) * temp_01

        F_HLL ( :) = ( S_R * F_L ( :) - S_L * F_R ( :)) * temp_01
     1             + ( U_R ( :) - U_L ( :)) * temp_02

        entropy_flux_re = ( S_R * entropy_flux_L
     1                    - S_L * entropy_flux_R) * temp_01
     1                  + ( entropy_R - entropy_L) * temp_02

        U_HLLEM ( :) = U_HLL ( :)  ! Ensure that "U_HLLEM ( :)" always initialized.
        F_HLLEM ( :) = F_HLL ( :)  ! Ensure that "F_HLLEM ( :)" always initialized.

        dU_HLLEM ( :) = 0.0
        dU_HLLEM_Min ( :) = 0.0
        dF_HLLEM ( :) = 0.0

        delta ( :) = 0.0
        eigenwt ( :) = 0.0

      END IF

! ----------------------------------------------------------------------
! Add the HLLEM parts to the pre-existing HLL RS. Do this only for
! the subsonic case.

      IF ( ( S_L < 0.0) .AND. ( S_R > 0.0) .AND.
     1     ( iHLLEM_choice > 0) ) THEN

! ------------------------------------------------------------
! First call "EIGENVECMHD_PTWISE_LINDEG". Load up "eigenval,
! lf_eigenvec, rt_eigenvec".

      iprim = 0  ! Need eigenvectors in conserved variables
      idirection = 1  ! in the x-direction.
      debug_stop = 0

#if ( COMPLETE_EIGENSYSTEM == 0)
      CALL EIGENVECMHD_PTWISE_LINDEG
#else
      CALL EIGENVECMHD_PTWISE_8
#endif

     1             ( iprim, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! ------------------------------------------------------------
! Now build "delta, eigenwt" for each of the waves of interest.
! The contributions from other waves are zeroed.

      delta ( :) = 0.0
      eigenwt ( :) = 0.0

      Del_U ( :) = U_R ( :) - U_L ( :)

! ------------------------------

      DO j = 1, NFIELD + 1 + NFLUID

! ----------

#if ( COMPLETE_EIGENSYSTEM == 0)

! Entropy wave and Alfven waves only.
        IF ( ( j == 2) .OR. ( j == 4) .OR. ( j == 6) ) THEN

#elif ( COMPLETE_EIGENSYSTEM == 1)

! Entropy wave, Alfven waves and slow magnetosonic waves.
! The Alfven + slow msonics carry most of the shear in the flow.
        IF ( ( j >= 2) .AND. ( j <= 6) ) THEN

#else

! Fast magnetosonic waves also included. The fast msonics carry most of the
! pressure jump.
        IF ( ( j >= 1) .AND. ( j <= 7) ) THEN

#endif

! ----------

        eigenval ( j) = AMIN1 ( AMAX1 ( eigenval ( j), S_L), S_R)

        IF ( iHLLEM_choice == 1) THEN

! Eqn. (28) of DB16.

          delta ( j) = 1.0
     1               - AMIN1 ( eigenval ( j), 0.0) / ( S_L - smallnum)
     1               - AMAX1 ( eigenval ( j), 0.0) / ( S_R + smallnum)

        ELSE

! Eqn. (2.11) of Betal16.

          temp_01 = - 3.0 * S_R * S_L / ( S_R - S_L)**2

          temp_02 = ( S_R - S_L)**2 / 3.0 - eigenval ( j) * ( S_R + S_L)

          temp_03 = S_R * AMIN1 ( eigenval ( j), 0.0)
     1            + S_L * AMAX1 ( eigenval ( j), 0.0) - S_R * S_L

          delta ( j) = temp_01

          IF ( temp_02 > 1.0e-6 * ( S_R - S_L)**2 ) delta ( j) =
     1         AMAX1 ( 0.0, AMIN1 ( delta ( j), temp_03 / temp_02) )

        END IF

        delta ( j) = delta ( j) * ( 1.0 - blend_HLL)

! ----------

        DO i = 1, NFIELD + 1 + NFLUID

        eigenwt ( j) = eigenwt ( j)
     1               + lf_eigenvec ( j, i) * Del_U ( i)

        END DO  ! End of "DO i = 1, NFIELD + 1 + NFLUID" loop.

! ----------

        END IF  ! End of "IF ( ( j == 2) .OR. ( j == 4) .OR. ( j == 6) )".

      END DO  ! End of "DO j = 1, NFIELD + 1 + NFLUID" loop.

! ------------------------------------------------------------
! Lastly, build "dU_HLLEM, dF_HLLEM". If the linear profile from "dU_HLLEM"
! gives physical states, add in the contributions to the fluxes.

      dU_HLLEM ( :) = 0.0
      dU_HLLEM_Min ( :) = 0.0
      dF_HLLEM ( :) = 0.0

! ------------------------------

      DO j = 1, NFIELD + 1 + NFLUID

! ----------

#if ( COMPLETE_EIGENSYSTEM == 0)

! Entropy wave and Alfven waves only.
        IF ( ( j == 2) .OR. ( j == 4) .OR. ( j == 6) ) THEN

#elif ( COMPLETE_EIGENSYSTEM == 1)

! Entropy wave, Alfven waves and slow magnetosonic waves. The Alfven + slow msonics
! cary most of the shear in the flow.
        IF ( ( j >= 2) .AND. ( j <= 6) ) THEN

#else

! Fast magnetosonic waves also included. The fast msonics carry most of the
! pressure jump.
        IF ( ( j >= 1) .AND. ( j <= 7) ) THEN

#endif

! ----------

! Eqn. (2.8) of Betal16.
        temp_01 = 2.0 * delta ( j) * eigenwt ( j)

        IF ( iHLLEM_choice == 1) THEN

! Eqn. (30) of DB16.
          temp_02 = - eigenwt ( j) * delta ( j) * S_R * S_L
     1            / ( S_R - S_L)

! "dU_HLLEM_Min" gets a contribution that is bounded by unity times "eigenwt".

          temp_03 = eigenwt ( j) * delta ( j)

! Another option is "eigenwt" itself. This is the steepest profile that will be physical
! when it comes to evaluating transverse fluxes.  With this profile the multid HLLEM
! will, in principle, revert to the 1D HLLEM for mesh-aligned flows.

!          temp_03 = eigenwt ( j)

        ELSE

! Eqn. (2.9b) of Betal16.
          temp_02 = 2.0 * eigenwt ( j) * delta ( j) * eigenval ( j)

! "dU_HLLEM_Min" gets a contribution that is bounded by unity times "eigenwt".

          temp_03 = eigenwt ( j) * 4.0 * delta ( j) / 3.0

! Another option is "eigenwt" itself. This is the steepest profile that will be physical
! when it comes to evaluating transverse fluxes.  With this profile the multid HLLEM
! will, in principle, revert to the 1D HLLEM for mesh-aligned flows.

!          temp_03 = eigenwt ( j)

        END IF

! ----------

        dU_HLLEM ( :) = dU_HLLEM ( :) + temp_01 * rt_eigenvec ( :, j)

        dU_HLLEM_Min ( :) = dU_HLLEM_Min ( :)
     1                                + temp_03 * rt_eigenvec ( :, j)

        dF_HLLEM ( :) = dF_HLLEM ( :) + temp_02 * rt_eigenvec ( :, j)

! ----------

        END IF  ! End of "IF ( ( j == 2) .OR. ( j == 4) .OR. ( j == 6) )".

      END DO  ! End of "DO j = 1, NFIELD + 1 + NFLUID" loop.

! ------------------------------------------------------------

! Density and pressure from HLL state should be bounded by densities and
! pressures from either side.

      rhobar = AMAX1 ( U_HLL ( 1), rhofloor)
      temp_01 = 1.0 / rhobar

      vlxbar = U_HLL ( 2) * temp_01
      vlybar = U_HLL ( 3) * temp_01
      vlzbar = U_HLL ( 4) * temp_01

      bfxbar = U_HLL ( 6)
      bfybar = U_HLL ( 7)
      bfzbar = U_HLL ( 8)

      prsbar = U_HLL ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             - ob_8pi * ( bfxbar**2 + bfybar**2 + bfzbar**2)

      prsbar = AMAX1 ( prsbar * ( gamma_eos - 1.0), prsfloor)

      IF ( ( rhobar < 0.9 * AMIN1 ( rho_L, rho_R) ) .OR.
     1     ( rhobar > 1.1 * AMAX1 ( rho_L, rho_R) ) .OR.
     1     ( prsbar < 0.9 * AMIN1 ( prs_L, prs_R) ) .OR.
     1     ( prsbar > 1.1 * AMAX1 ( prs_L, prs_R) ) ) THEN

        dU_HLLEM ( :) = 0.0
        dF_HLLEM ( :) = 0.0

      END IF

! ------------------------------
! Examine right state within the Riemann fan for positivity.

      Del_U ( :) = U_HLL ( :) + 0.5 * dU_HLLEM ( :)


      rhobar = AMAX1 ( Del_U ( 1), rhofloor)
      temp_01 = 1.0 / rhobar

      vlxbar = Del_U ( 2) * temp_01
      vlybar = Del_U ( 3) * temp_01
      vlzbar = Del_U ( 4) * temp_01

      bfxbar = Del_U ( 6)
      bfybar = Del_U ( 7)
      bfzbar = Del_U ( 8)

      prsbar = Del_U ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             - ob_8pi * ( bfxbar**2 + bfybar**2 + bfzbar**2)

      prsbar = AMAX1 ( prsbar * ( gamma_eos - 1.0), prsfloor)

      IF ( ( rhobar <= rhofloor_plus_delta) .OR.
     1     ( prsbar <= prsfloor_plus_delta) ) THEN

        dU_HLLEM ( :) = 0.0
        dF_HLLEM ( :) = 0.0

      END IF

! ------------------------------
! Examine left state within the Riemann fan for positivity.

      Del_U ( :) = U_HLL ( :) - 0.5 * dU_HLLEM ( :)


      rhobar = AMAX1 ( Del_U ( 1), rhofloor)
      temp_01 = 1.0 / rhobar

      vlxbar = Del_U ( 2) * temp_01
      vlybar = Del_U ( 3) * temp_01
      vlzbar = Del_U ( 4) * temp_01

      bfxbar = Del_U ( 6)
      bfybar = Del_U ( 7)
      bfzbar = Del_U ( 8)

      prsbar = Del_U ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             - ob_8pi * ( bfxbar**2 + bfybar**2 + bfzbar**2)

      prsbar = AMAX1 ( prsbar * ( gamma_eos - 1.0), prsfloor)

      IF ( ( rhobar <= rhofloor_plus_delta) .OR.
     1     ( prsbar <= prsfloor_plus_delta) ) THEN

        dU_HLLEM ( :) = 0.0
        dF_HLLEM ( :) = 0.0

      END IF

! ------------------------------
! Add in contributions from "dU_HLLEM" and "dF_HLLEM" to the state
! and fluxes. This gives us "U_HLLEM, F_HLLEM".

        temp_01 = - 0.5 * ( S_R + S_L) / ( S_R - S_L)

! Eqns. (2.1) and (2.3) of Betal16.
        U_HLLEM ( :) = U_HLL ( :) + temp_01 * dU_HLLEM ( :)

! ----------

        IF ( iHLLEM_choice == 1) THEN

! Eqn. (30) of DB16.
          F_HLLEM ( :) = F_HLL ( :) + dF_HLLEM ( :)

        ELSE

! Eqn. (2.9b) of Betal16.
          temp_02 = ( S_R - S_L) / 6.0

          F_HLLEM ( :) = F_HLL ( :) + temp_02 * dU_HLLEM ( :)
     1                 + temp_01 * dF_HLLEM ( :)

        END IF

! ------------------------------------------------------------

      END IF  ! End of "IF ( ( S_L < 0.0) .AND. ( S_R > 0.0) )" loop.

! ----------------------------------------------------------------------
! Wrap up and exit subroutine.

      flux_re ( :) = F_HLLEM ( :)

      flux_re ( 6) = F_HLLEM ( 7)
      flux_re ( 7) = F_HLLEM ( 8)
      flux_re ( 8) = 0.0

! ----------
! Also split the momentum flux_re into advected and pressure parts and
! store the pressure part separately. Do this only if it is called for.

      IF ( split_prsflux_re .EQ. 1) THEN

      rhobar = AMAX1 ( U_HLLEM ( 1), rhofloor)
      temp_01 = 1.0 / rhobar

      vlxbar = U_HLLEM ( 2) * temp_01
      vlybar = U_HLLEM ( 3) * temp_01
      vlzbar = U_HLLEM ( 4) * temp_01

      bfxbar = U_HLLEM ( 6)
      bfybar = U_HLLEM ( 7)
      bfzbar = U_HLLEM ( 8)

      prsbar = U_HLLEM ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             - ob_8pi * ( bfxbar**2 + bfybar**2 + bfzbar**2)

      prsbar = AMAX1 ( prsbar * ( gamma_eos - 1.0), prsfloor)

! This is the part of the x-momentum flux that acts like a pressure.
      flux_re ( 8) = flux_re ( 2) - rhobar * vlxbar**2
     1             + ob_4pi * bfxre**2

      flux_re ( 2) = flux_re ( 2) - flux_re ( 8)

      END IF

! ----------------------------------------------------------------------

      END SUBROUTINE MHD_RIEM_HLLEM_PTWISE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MHD_TRANSVERSE_FLUXES_PTWISE

     1     ( gammabar, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     S_L, S_R, S_Mid, neta_x, neta_y,

     1     U_HLL, F_HLL, G_HLL,
     1     U_HLLEM, F_HLLEM, G_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM,

     1     U_Minus, F_Minus, G_Minus, U_Plus, F_Plus, G_Plus)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes the outputs from the 1D HLLEM solver for MHD and makes all the
! transverse fluxes. Also does internal rotations in the components of
! the state and flux vectors so as to produce the same in a global space.

! ----------------------------------------------------------------------

! Inputs:-

! gammabar : Averaged gamma.

! "iHLLEM_choice = 1" gives original Dumbser-Balsara HLLEM; it has
! slighly less dissipation for fast-moving waves. "iHLLEM_choice = 2" gives
! Balsara et al HLLEM; it is more stable for strong shocks.
! Both will capture stationary intermediate waves. "iHLLEM_choice = 0" gives
! HLL. Keep this consistent across subroutines.

! S_L, S_R : Extremal speeds of the HLL Riemann solver.

! neta_x, neta_y : Outward normal to the face of the multidimensional wave
!            model. (tau, z, neta) form a right-handed triplet.

! U_HLL, F_HLL : HLL state and flux. On input, the 2,3,4 and 6,7,8
!            components are (tau, z, neta)-directional velocities and
!            magnetic fields. Besides, F_HLL is flux in the tau-direction.

! U_HLLEM, F_HLLEM : HLLEM state and flux. Same logic as above.

! dU_HLLEM, dF_HLLEM : These are the additions that turn an HLL state or
!            flux into an HLLEM flux state or flux.
! dU_HLLEM_Min : Reduced slope in the linear profile to avoid double steepening
!               of the multidimensional profile.

! Outputs:-

! S_Mid : Middle of the two extremal speeds.

! U_HLL, F_HLL, G_HLL : HLL state and both fluxes. On output, the
!            2,3,4 and 6,7,8 components are x,y,z-directional velocities
!            and magnetic fields. Besides F_HLL, G_HLL are x- and y-fluxes.

! U_HLLEM, F_HLLEM, G_HLLEM : HLLEM state and both fluxes. On output, the
!            2,3,4 and 6,7,8 components are x,y,z-directional velocities
!            and magnetic fields. Besides F_HLLEM, G_HLLEM are x- and
!            y-fluxes.

! U_Minus, U_Plus, F_Minus, F_Plus, G_Minus, G_Plus : For HLLEM,
!        we fit a linear profile in the Riemann fan. As a result,
!        we split the fan [S_L, S_R] into two equal halves [S_L, S_Mid],
!        [S_Mid, S_R]. These are the states and fluxes in each of those halves.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER use_approx_eos, iHLLEM_choice

      REAL gammabar, smallnum, prsfloor, rhofloor, pi, gasconst,

     1     cp_specie ( 0 : NFLUID), molewt_specie ( 0 : NFLUID),
     1     heat_of_formation ( 0 : NFLUID),

     1     S_L, S_R, S_Mid, neta_x, neta_y

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: U_HLL, F_HLL, G_HLL,
     1     U_HLLEM, F_HLLEM, G_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM,

     1     U_Minus, F_Minus, G_Minus, U_Plus, F_Plus, G_Plus

! ----------------------------------------------------------------------

! Variables past this point are not used in the subroutine interface.

      REAL temp_01, temp_02, temp_03, fourpi, ob_4pi, ob_8pi,
     1     rho, prs, vlx, vly, vlz, bfx, bfy, bfz, v2, b2, vdb,
     1     vlneta, vltau, bfneta, bftau

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: U_temp, F_temp, G_temp

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 * ob_4pi

      S_Mid = 0.5 * ( S_L + S_R)  ! Analogue of contact discontinuity in HLLEM.

! ----------------------------------------------------------------------

! First build "U_Minus, U_Plus, F_Minus, F_Plus" for the choice of
! "iHLLEM_choice".

! Realize that 1D HLLEM steepens the profile already. We have two choices:-

! 1) Either go out of this subroutine with the steepened profile and don't
! steepen it again later. This is the choice used here.
! In that case "U_Minus, U_Plus" get updated with "dU_HLLEM".

! 2) Alternatively, don't oversteepen the profile here. But then you have to
! steepen it later. This choice does not work as well in certain instances.
! In that case "U_Minus, U_Plus" get updated with "dU_HLLEM_Min".

! ------------------------------

      U_Minus ( :) = U_HLL ( :) - 0.25 * dU_HLLEM_Min ( :)

      IF ( iHLLEM_choice == 1) THEN

! Eqn. (2.14) of Betal16.
        temp_01 = 0.75 * S_L + 0.25 * S_R
        temp_02 = 0.5 * ( temp_01**2 - S_R * S_L) / ( S_R - S_L)

        F_Minus ( :) = F_HLL ( :) + temp_02 * dU_HLLEM ( :)

      ELSE

! Eqn. (2.12) of Betal16.
        temp_01 = ( S_R - S_L) / 6.0

        F_Minus ( :) = F_HLL ( :) + temp_01 * dU_HLLEM ( :)
     1                            - 0.25 * dF_HLLEM ( :)

      END IF

! ------------------------------

      U_Plus ( :) = U_HLL ( :) + 0.25 * dU_HLLEM_Min ( :)

      IF ( iHLLEM_choice == 1) THEN

! Eqn. (2.14) of Betal16.
        temp_01 = 0.25 * S_L + 0.75 * S_R
        temp_02 = 0.5 * ( temp_01**2 - S_R * S_L) / ( S_R - S_L)

        F_Plus ( :) = F_HLL ( :) + temp_02 * dU_HLLEM ( :)

      ELSE

! Eqn. (2.12) of Betal16.
        temp_01 = ( S_R - S_L) / 6.0

        F_Plus ( :) = F_HLL ( :) + temp_01 * dU_HLLEM ( :)
     1                           + 0.25 * dF_HLLEM ( :)

      END IF

! ----------------------------------------------------------------------
! Build primitives from "U_Minus ( :), F_Minus ( :)".
! Then build "G_Minus ( :)" in the neta-direction.
! Make external rotations of the flux vectors "F_Minus ( :), G_Minus ( :)".
! Make reverse rotations within the 2,4 and 6,8 components of state & fluxes.

! ------------------------------
! Make primitives.

      rho = AMAX1 ( U_Minus ( 1), rhofloor)
      temp_01 = 1.0 / rho

      vltau  = U_Minus ( 2) * temp_01
      vlz    = U_Minus ( 3) * temp_01
      vlneta = U_Minus ( 4) * temp_01

      bftau  = U_Minus ( 6)
      bfz    = U_Minus ( 7)
      bfneta = U_Minus ( 8)

      vlx = neta_x * vlneta - neta_y * vltau
      vly = neta_y * vlneta + neta_x * vltau

      bfx = neta_x * bfneta - neta_y * bftau
      bfy = neta_y * bfneta + neta_x * bftau

      v2 = vlx**2 + vly**2 + vlz**2
      b2 = bfx**2 + bfy**2 + bfz**2
      vdb = vlx * bfx + vly * bfy + vlz * bfz

#if ( 1 == 1)
      prs = U_Minus ( 5) - 0.5 * rho * v2 - ob_8pi * b2
      prs = AMAX1 ( prs * ( gammabar - 1.0), prsfloor)
#else
      prs = F_Minus ( 2) - rho * vltau**2 - ob_8pi * b2
     1                                    + ob_4pi * bftau**2
      prs = AMAX1 ( prs, prsfloor)
#endif

! ------------------------------
! "G_Minus ( :)" is flux in the neta-direction, still in the rotated coordinates.
! Primitives from above are needed to make this flux.

      G_Minus ( :) = vlneta * U_Minus ( :)

      G_Minus ( 2) = G_Minus ( 2) - ob_4pi * bfneta * bftau
      G_Minus ( 3) = G_Minus ( 3) - ob_4pi * bfneta * bfz
      G_Minus ( 4) = G_Minus ( 4) + prs + ob_8pi * b2
     1                                  - ob_4pi * bfneta**2

      G_Minus ( 5) = G_Minus ( 5) + ( prs + ob_8pi * b2 ) * vlneta
     1             - ob_4pi * bfneta * vdb

      G_Minus ( 6) = G_Minus ( 6) - bfneta * vltau
      G_Minus ( 7) = G_Minus ( 7) - bfneta * vlz
      G_Minus ( 8) = G_Minus ( 8) - bfneta * vlneta  ! This should be zero.

! ------------------------------
! External rotations of the flux vectors "F_Minus ( :), G_Minus ( :)" to get
! "F_temp ( :), G_temp ( :)".

      U_temp ( :) = U_Minus ( :)
      F_temp ( :) = neta_x * G_Minus ( :) - neta_y * F_Minus ( :)
      G_temp ( :) = neta_y * G_Minus ( :) + neta_x * F_Minus ( :)

! ------------------------------
! Internal rotations applied to "U_Minus ( :), F_Minus ( :), G_Minus ( :)".

      U_Minus ( 1) = U_temp ( 1)

      U_Minus ( 2) = neta_x * U_temp ( 4) - neta_y * U_temp ( 2)
      U_Minus ( 3) = neta_y * U_temp ( 4) + neta_x * U_temp ( 2)
      U_Minus ( 4) = U_temp ( 3)

      U_Minus ( 5) = U_temp ( 5)

      U_Minus ( 6) = neta_x * U_temp ( 8) - neta_y * U_temp ( 6)
      U_Minus ( 7) = neta_y * U_temp ( 8) + neta_x * U_temp ( 6)
      U_Minus ( 8) = U_temp ( 7)

#if ( NFLUID >= 1)
      U_Minus ( 9: ) = U_temp ( 9: )
#endif

! ----------

      F_Minus ( 1) = F_temp ( 1)

      F_Minus ( 2) = neta_x * F_temp ( 4) - neta_y * F_temp ( 2)
      F_Minus ( 3) = neta_y * F_temp ( 4) + neta_x * F_temp ( 2)
      F_Minus ( 4) = F_temp ( 3)

      F_Minus ( 5) = F_temp ( 5)

      F_Minus ( 6) = neta_x * F_temp ( 8) - neta_y * F_temp ( 6)
      F_Minus ( 7) = neta_y * F_temp ( 8) + neta_x * F_temp ( 6)
      F_Minus ( 8) = F_temp ( 7)

#if ( NFLUID >= 1)
      F_Minus ( 9: ) = F_temp ( 9: )
#endif

! ----------

      G_Minus ( 1) = G_temp ( 1)

      G_Minus ( 2) = neta_x * G_temp ( 4) - neta_y * G_temp ( 2)
      G_Minus ( 3) = neta_y * G_temp ( 4) + neta_x * G_temp ( 2)
      G_Minus ( 4) = G_temp ( 3)

      G_Minus ( 5) = G_temp ( 5)

      G_Minus ( 6) = neta_x * G_temp ( 8) - neta_y * G_temp ( 6)
      G_Minus ( 7) = neta_y * G_temp ( 8) + neta_x * G_temp ( 6)
      G_Minus ( 8) = G_temp ( 7)

#if ( NFLUID >= 1)
      G_Minus ( 9: ) = G_temp ( 9: )
#endif

! ------------------------------
! We usually need "G_Minus ( :)" in the neta direction and "F_Minus ( :)"
! in the tau direction. However, we want the Cartesian ordering of the
! flux components. That is done here.

      F_temp ( :) = F_Minus ( :)
      G_temp ( :) = G_Minus ( :)
      F_Minus ( :) = neta_x * G_temp ( :) - neta_y * F_temp ( :)
      G_Minus ( :) = neta_y * G_temp ( :) + neta_x * F_temp ( :)

! ----------------------------------------------------------------------
! Build primitives from "U_Plus ( :), F_Plus ( :)".
! Then build "G_Plus ( :)" in the neta-direction.
! Make external rotations of the flux vectors "F_Plus ( :), G_Plus ( :)".
! Make reverse rotations within the 2,4 and 6,8 components of state & fluxes.

! ------------------------------
! Make primitives.

      rho = AMAX1 ( U_Plus ( 1), rhofloor)
      temp_01 = 1.0 / rho

      vltau  = U_Plus ( 2) * temp_01
      vlz    = U_Plus ( 3) * temp_01
      vlneta = U_Plus ( 4) * temp_01

      bftau  = U_Plus ( 6)
      bfz    = U_Plus ( 7)
      bfneta = U_Plus ( 8)

      vlx = neta_x * vlneta - neta_y * vltau
      vly = neta_y * vlneta + neta_x * vltau

      bfx = neta_x * bfneta - neta_y * bftau
      bfy = neta_y * bfneta + neta_x * bftau

      v2 = vlx**2 + vly**2 + vlz**2
      b2 = bfx**2 + bfy**2 + bfz**2
      vdb = vlx * bfx + vly * bfy + vlz * bfz

#if ( 1 == 1)
      prs = U_Plus ( 5) - 0.5 * rho * v2 - ob_8pi * b2
      prs = AMAX1 ( prs * ( gammabar - 1.0), prsfloor)
#else
      prs = F_Plus ( 2) - rho * vltau**2 - ob_8pi * b2
     1                                    + ob_4pi * bftau**2
      prs = AMAX1 ( prs, prsfloor)
#endif

! ------------------------------
! "G_Plus ( :)" is flux in the neta-direction, still in the rotated coordinates.
! Primitives from above are needed to make this flux.

      G_Plus ( :) = vlneta * U_Plus ( :)

      G_Plus ( 2) = G_Plus ( 2) - ob_4pi * bfneta * bftau
      G_Plus ( 3) = G_Plus ( 3) - ob_4pi * bfneta * bfz
      G_Plus ( 4) = G_Plus ( 4) + prs + ob_8pi * b2
     1                                  - ob_4pi * bfneta**2

      G_Plus ( 5) = G_Plus ( 5) + ( prs + ob_8pi * b2 ) * vlneta
     1             - ob_4pi * bfneta * vdb

      G_Plus ( 6) = G_Plus ( 6) - bfneta * vltau
      G_Plus ( 7) = G_Plus ( 7) - bfneta * vlz
      G_Plus ( 8) = G_Plus ( 8) - bfneta * vlneta  ! This should be zero.

! ------------------------------
! External rotations of the flux vectors "F_Plus ( :), G_Plus ( :)" to get
! "F_temp ( :), G_temp ( :)".

      U_temp ( :) = U_Plus ( :)
      F_temp ( :) = neta_x * G_Plus ( :) - neta_y * F_Plus ( :)
      G_temp ( :) = neta_y * G_Plus ( :) + neta_x * F_Plus ( :)

! ------------------------------
! Internal rotations applied to "U_Plus ( :), F_Plus ( :), G_Plus ( :)".

      U_Plus ( 1) = U_temp ( 1)

      U_Plus ( 2) = neta_x * U_temp ( 4) - neta_y * U_temp ( 2)
      U_Plus ( 3) = neta_y * U_temp ( 4) + neta_x * U_temp ( 2)
      U_Plus ( 4) = U_temp ( 3)

      U_Plus ( 5) = U_temp ( 5)

      U_Plus ( 6) = neta_x * U_temp ( 8) - neta_y * U_temp ( 6)
      U_Plus ( 7) = neta_y * U_temp ( 8) + neta_x * U_temp ( 6)
      U_Plus ( 8) = U_temp ( 7)

#if ( NFLUID >= 1)
      U_Plus ( 9: ) = U_temp ( 9: )
#endif

! ----------

      F_Plus ( 1) = F_temp ( 1)

      F_Plus ( 2) = neta_x * F_temp ( 4) - neta_y * F_temp ( 2)
      F_Plus ( 3) = neta_y * F_temp ( 4) + neta_x * F_temp ( 2)
      F_Plus ( 4) = F_temp ( 3)

      F_Plus ( 5) = F_temp ( 5)

      F_Plus ( 6) = neta_x * F_temp ( 8) - neta_y * F_temp ( 6)
      F_Plus ( 7) = neta_y * F_temp ( 8) + neta_x * F_temp ( 6)
      F_Plus ( 8) = F_temp ( 7)

#if ( NFLUID >= 1)
      F_Plus ( 9: ) = F_temp ( 9: )
#endif

! ----------

      G_Plus ( 1) = G_temp ( 1)

      G_Plus ( 2) = neta_x * G_temp ( 4) - neta_y * G_temp ( 2)
      G_Plus ( 3) = neta_y * G_temp ( 4) + neta_x * G_temp ( 2)
      G_Plus ( 4) = G_temp ( 3)

      G_Plus ( 5) = G_temp ( 5)

      G_Plus ( 6) = neta_x * G_temp ( 8) - neta_y * G_temp ( 6)
      G_Plus ( 7) = neta_y * G_temp ( 8) + neta_x * G_temp ( 6)
      G_Plus ( 8) = G_temp ( 7)

#if ( NFLUID >= 1)
      G_Plus ( 9: ) = G_temp ( 9: )
#endif

! ------------------------------
! We usually need "G_Plus ( :)" in the neta direction and "F_Plus ( :)"
! in the tau direction. However, we want the Cartesian ordering of the
! flux components. That is done here.

      F_temp ( :) = F_Plus ( :)
      G_temp ( :) = G_Plus ( :)
      F_Plus ( :) = neta_x * G_temp ( :) - neta_y * F_temp ( :)
      G_Plus ( :) = neta_y * G_temp ( :) + neta_x * F_temp ( :)

! ----------------------------------------------------------------------
! Build primitives from "U_HLL ( :), F_HLL ( :)".
! Then build "G_HLL ( :)" in the neta-direction.
! Make external rotations of the flux vectors "F_HLL ( :), G_HLL ( :)".
! Make reverse rotations within the 2,4 and 6,8 components of state & fluxes.

! ------------------------------
! Make primitives.

      rho = AMAX1 ( U_HLL ( 1), rhofloor)
      temp_01 = 1.0 / rho

      vltau  = U_HLL ( 2) * temp_01
      vlz    = U_HLL ( 3) * temp_01
      vlneta = U_HLL ( 4) * temp_01

      bftau  = U_HLL ( 6)
      bfz    = U_HLL ( 7)
      bfneta = U_HLL ( 8)

      vlx = neta_x * vlneta - neta_y * vltau
      vly = neta_y * vlneta + neta_x * vltau

      bfx = neta_x * bfneta - neta_y * bftau
      bfy = neta_y * bfneta + neta_x * bftau

      v2 = vlx**2 + vly**2 + vlz**2
      b2 = bfx**2 + bfy**2 + bfz**2
      vdb = vlx * bfx + vly * bfy + vlz * bfz

#if ( 1 == 1)
      prs = U_HLL ( 5) - 0.5 * rho * v2 - ob_8pi * b2
      prs = AMAX1 ( prs * ( gammabar - 1.0), prsfloor)
#else
      prs = F_HLL ( 2) - rho * vltau**2 - ob_8pi * b2
     1                                    + ob_4pi * bftau**2
      prs = AMAX1 ( prs, prsfloor)
#endif

! ------------------------------
! "G_HLL ( :)" is flux in the neta-direction, still in the rotated coordinates.
! Primitives from above are needed to make this flux.

      G_HLL ( :) = vlneta * U_HLL ( :)

      G_HLL ( 2) = G_HLL ( 2) - ob_4pi * bfneta * bftau
      G_HLL ( 3) = G_HLL ( 3) - ob_4pi * bfneta * bfz
      G_HLL ( 4) = G_HLL ( 4) + prs + ob_8pi * b2
     1                                  - ob_4pi * bfneta**2

      G_HLL ( 5) = G_HLL ( 5) + ( prs + ob_8pi * b2 ) * vlneta
     1             - ob_4pi * bfneta * vdb

      G_HLL ( 6) = G_HLL ( 6) - bfneta * vltau
      G_HLL ( 7) = G_HLL ( 7) - bfneta * vlz
      G_HLL ( 8) = G_HLL ( 8) - bfneta * vlneta  ! This should be zero.

! ------------------------------
! External rotations of the flux vectors "F_HLL ( :), G_HLL ( :)" to get
! "F_temp ( :), G_temp ( :)".

      U_temp ( :) = U_HLL ( :)
      F_temp ( :) = neta_x * G_HLL ( :) - neta_y * F_HLL ( :)
      G_temp ( :) = neta_y * G_HLL ( :) + neta_x * F_HLL ( :)

! ------------------------------
! Internal rotations applied to "U_HLL ( :), F_HLL ( :), G_HLL ( :)".

      U_HLL ( 1) = U_temp ( 1)

      U_HLL ( 2) = neta_x * U_temp ( 4) - neta_y * U_temp ( 2)
      U_HLL ( 3) = neta_y * U_temp ( 4) + neta_x * U_temp ( 2)
      U_HLL ( 4) = U_temp ( 3)

      U_HLL ( 5) = U_temp ( 5)

      U_HLL ( 6) = neta_x * U_temp ( 8) - neta_y * U_temp ( 6)
      U_HLL ( 7) = neta_y * U_temp ( 8) + neta_x * U_temp ( 6)
      U_HLL ( 8) = U_temp ( 7)

#if ( NFLUID >= 1)
      U_HLL ( 9: ) = U_temp ( 9: )
#endif

! ----------

      F_HLL ( 1) = F_temp ( 1)

      F_HLL ( 2) = neta_x * F_temp ( 4) - neta_y * F_temp ( 2)
      F_HLL ( 3) = neta_y * F_temp ( 4) + neta_x * F_temp ( 2)
      F_HLL ( 4) = F_temp ( 3)

      F_HLL ( 5) = F_temp ( 5)

      F_HLL ( 6) = neta_x * F_temp ( 8) - neta_y * F_temp ( 6)
      F_HLL ( 7) = neta_y * F_temp ( 8) + neta_x * F_temp ( 6)
      F_HLL ( 8) = F_temp ( 7)

#if ( NFLUID >= 1)
      F_HLL ( 9: ) = F_temp ( 9: )
#endif

! ----------

      G_HLL ( 1) = G_temp ( 1)

      G_HLL ( 2) = neta_x * G_temp ( 4) - neta_y * G_temp ( 2)
      G_HLL ( 3) = neta_y * G_temp ( 4) + neta_x * G_temp ( 2)
      G_HLL ( 4) = G_temp ( 3)

      G_HLL ( 5) = G_temp ( 5)

      G_HLL ( 6) = neta_x * G_temp ( 8) - neta_y * G_temp ( 6)
      G_HLL ( 7) = neta_y * G_temp ( 8) + neta_x * G_temp ( 6)
      G_HLL ( 8) = G_temp ( 7)

#if ( NFLUID >= 1)
      G_HLL ( 9: ) = G_temp ( 9: )
#endif

! ------------------------------
! We usually need "G_HLL ( :)" in the neta direction and "F_HLL ( :)"
! in the tau direction. However, we want the Cartesian ordering of the
! flux components. That is done here.

      F_temp ( :) = F_HLL ( :)
      G_temp ( :) = G_HLL ( :)
      F_HLL ( :) = neta_x * G_temp ( :) - neta_y * F_temp ( :)
      G_HLL ( :) = neta_y * G_temp ( :) + neta_x * F_temp ( :)

! ----------------------------------------------------------------------
! Build primitives from "U_HLLEM ( :), F_HLLEM ( :)".
! Then build "G_HLLEM ( :)" in the neta-direction.
! Make external rotations of the flux vectors "F_HLLEM ( :), G_HLLEM ( :)".
! Make reverse rotations within the 2,4 and 6,8 components of state & fluxes.

! ------------------------------
! Make primitives.

      rho = AMAX1 ( U_HLLEM ( 1), rhofloor)
      temp_01 = 1.0 / rho

      vltau  = U_HLLEM ( 2) * temp_01
      vlz    = U_HLLEM ( 3) * temp_01
      vlneta = U_HLLEM ( 4) * temp_01

      bftau  = U_HLLEM ( 6)
      bfz    = U_HLLEM ( 7)
      bfneta = U_HLLEM ( 8)

      vlx = neta_x * vlneta - neta_y * vltau
      vly = neta_y * vlneta + neta_x * vltau

      bfx = neta_x * bfneta - neta_y * bftau
      bfy = neta_y * bfneta + neta_x * bftau

      v2 = vlx**2 + vly**2 + vlz**2
      b2 = bfx**2 + bfy**2 + bfz**2
      vdb = vlx * bfx + vly * bfy + vlz * bfz

#if ( 1 == 1)
      prs = U_HLLEM ( 5) - 0.5 * rho * v2 - ob_8pi * b2
      prs = AMAX1 ( prs * ( gammabar - 1.0), prsfloor)
#else
      prs = F_HLLEM ( 2) - rho * vltau**2 - ob_8pi * b2
     1                                    + ob_4pi * bftau**2
      prs = AMAX1 ( prs, prsfloor)
#endif

! ------------------------------
! "G_HLLEM ( :)" is flux in the neta-direction, still in the rotated coordinates.
! Primitives from above are needed to make this flux.

      G_HLLEM ( :) = vlneta * U_HLLEM ( :)

      G_HLLEM ( 2) = G_HLLEM ( 2) - ob_4pi * bfneta * bftau
      G_HLLEM ( 3) = G_HLLEM ( 3) - ob_4pi * bfneta * bfz
      G_HLLEM ( 4) = G_HLLEM ( 4) + prs + ob_8pi * b2
     1                                  - ob_4pi * bfneta**2

      G_HLLEM ( 5) = G_HLLEM ( 5) + ( prs + ob_8pi * b2 ) * vlneta
     1             - ob_4pi * bfneta * vdb

      G_HLLEM ( 6) = G_HLLEM ( 6) - bfneta * vltau
      G_HLLEM ( 7) = G_HLLEM ( 7) - bfneta * vlz
      G_HLLEM ( 8) = G_HLLEM ( 8) - bfneta * vlneta  ! This should be zero.

! ------------------------------
! External rotations of the flux vectors "F_HLLEM ( :), G_HLLEM ( :)" to get
! "F_temp ( :), G_temp ( :)".

      U_temp ( :) = U_HLLEM ( :)
      F_temp ( :) = neta_x * G_HLLEM ( :) - neta_y * F_HLLEM ( :)
      G_temp ( :) = neta_y * G_HLLEM ( :) + neta_x * F_HLLEM ( :)

! ------------------------------
! Internal rotations applied to "U_HLLEM ( :), F_HLLEM ( :), G_HLLEM ( :)".

      U_HLLEM ( 1) = U_temp ( 1)

      U_HLLEM ( 2) = neta_x * U_temp ( 4) - neta_y * U_temp ( 2)
      U_HLLEM ( 3) = neta_y * U_temp ( 4) + neta_x * U_temp ( 2)
      U_HLLEM ( 4) = U_temp ( 3)

      U_HLLEM ( 5) = U_temp ( 5)

      U_HLLEM ( 6) = neta_x * U_temp ( 8) - neta_y * U_temp ( 6)
      U_HLLEM ( 7) = neta_y * U_temp ( 8) + neta_x * U_temp ( 6)
      U_HLLEM ( 8) = U_temp ( 7)

#if ( NFLUID >= 1)
      U_HLLEM ( 9: ) = U_temp ( 9: )
#endif

! ----------

      F_HLLEM ( 1) = F_temp ( 1)

      F_HLLEM ( 2) = neta_x * F_temp ( 4) - neta_y * F_temp ( 2)
      F_HLLEM ( 3) = neta_y * F_temp ( 4) + neta_x * F_temp ( 2)
      F_HLLEM ( 4) = F_temp ( 3)

      F_HLLEM ( 5) = F_temp ( 5)

      F_HLLEM ( 6) = neta_x * F_temp ( 8) - neta_y * F_temp ( 6)
      F_HLLEM ( 7) = neta_y * F_temp ( 8) + neta_x * F_temp ( 6)
      F_HLLEM ( 8) = F_temp ( 7)

#if ( NFLUID >= 1)
      F_HLLEM ( 9: ) = F_temp ( 9: )
#endif

! ----------

      G_HLLEM ( 1) = G_temp ( 1)

      G_HLLEM ( 2) = neta_x * G_temp ( 4) - neta_y * G_temp ( 2)
      G_HLLEM ( 3) = neta_y * G_temp ( 4) + neta_x * G_temp ( 2)
      G_HLLEM ( 4) = G_temp ( 3)

      G_HLLEM ( 5) = G_temp ( 5)

      G_HLLEM ( 6) = neta_x * G_temp ( 8) - neta_y * G_temp ( 6)
      G_HLLEM ( 7) = neta_y * G_temp ( 8) + neta_x * G_temp ( 6)
      G_HLLEM ( 8) = G_temp ( 7)

#if ( NFLUID >= 1)
      G_HLLEM ( 9: ) = G_temp ( 9: )
#endif

! ------------------------------
! We usually need "G_HLLEM ( :)" in the neta direction and "F_HLLEM ( :)"
! in the tau direction. However, we want the Cartesian ordering of the
! flux components. That is done here.

      F_temp ( :) = F_HLLEM ( :)
      G_temp ( :) = G_HLLEM ( :)
      F_HLLEM ( :) = neta_x * G_temp ( :) - neta_y * F_temp ( :)
      G_HLLEM ( :) = neta_y * G_temp ( :) + neta_x * F_temp ( :)

! ----------------------------------------------------------------------

      END SUBROUTINE MHD_TRANSVERSE_FLUXES_PTWISE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MHD_RIEM_HLLEM_DB

     1   ( nriembegin, nriemend,

     1     double_dissipate_efz, double_dissipate_efy,

     1     gamma, smallnum, prsfloor, rhofloor, pi, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        bfxre,
     1        bfyremi, bfyrepl,
     1        bfzremi, bfzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, bfy0, bfz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Packages the pointwise HLLEM RS in a form that works with the
! dimension-by-dimension Riemann solver's interface. This is just
! useful for testing.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL RHO_PRS_IN_ENG_T_GAMA_OUT

! ----------------------------------------------------------------------

      INTEGER :: nriembegin, nriemend, 
     1        double_dissipate_efz, double_dissipate_efy,
     1        ientropyfix, iprimitive, use_approx_eos,
     1        split_prsflux_re

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor, 
     1     pi, fracdiff, gasconst, divratio_einfeldt_rs

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     bfxre, bfyremi, bfyrepl, bfzremi, bfzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID) ::
     1    lf_eigenvec_re, rt_eigenvec_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID) ::
     1     eigenval_re, eigenwt_re, conserved_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD + NFLUID + 1) ::
     1     flux_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD + NFLUID) ::
     1      zonebdy_flux, zonebdy_conserved

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho0, prs0, vlx0, vly0, vlz0, bfy0, bfz0,
     1     entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Past this point we only declare variables that are internal to the
! subroutine.

      INTEGER i, iHLLEM_choice

      REAL rho_L, rho_R, prs_L, prs_R,
     1     vlx_L, vlx_R, vly_L, vly_R, vlz_L, vlz_R,

     1     bfxpt, bfy_L, bfy_R, bfz_L, bfz_R,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     blend_HLL, S_L, S_R, entropy_flux_pt

      REAL, DIMENSION ( 0: NFLUID) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: flux_pt,
     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

! ----------------------------------------------------------------------

      blend_HLL = 0.0  ! We want full-fledged HLLEM Riemann solver.

      iHLLEM_choice = 1

      DO i = nriembegin, nriemend

! ------------------------------

      rho_L = rhoremi ( i)
      prs_L = prsremi ( i)

      vlx_L = vlxremi ( i)
      vly_L = vlyremi ( i)
      vlz_L = vlzremi ( i)

      bfy_L = bfyremi ( i)
      bfz_L = bfzremi ( i)

      mean_mol_wt_L = mean_mol_wt_remi ( i)
      gamma_eos_L = gamma_eos_remi ( i)
      gamma_soundspeed_L = gamma_soundspeed_remi ( i)


      bfxpt = bfxre ( i)


      rho_R = rhorepl ( i)
      prs_R = prsrepl ( i)

      vlx_R = vlxrepl ( i)
      vly_R = vlyrepl ( i)
      vlz_R = vlzrepl ( i)

      bfy_R = bfyrepl ( i)
      bfz_R = bfzrepl ( i)

      mean_mol_wt_R = mean_mol_wt_repl ( i)
      gamma_eos_R = gamma_eos_repl ( i)
      gamma_soundspeed_R = gamma_soundspeed_repl ( i)

! ------------------------------

      CALL MHD_RIEM_HLLEM_PTWISE

     1     ( gamma, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vlx_L, vlx_R,
     1        vly_L, vly_R,
     1        vlz_L, vlz_R,
     1        bfxpt,
     1        bfy_L, bfy_R,
     1        bfz_L, bfz_R,
     1        specie_L, specie_R,

     1        gamma_eos_L, gamma_eos_R,
     1        gamma_soundspeed_L, gamma_soundspeed_R,
     1        mean_mol_wt_L, mean_mol_wt_R,

     1        blend_HLL,

     1        S_L, S_R,
     1        F_L, F_R, U_HLL, F_HLL,
     1        U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM,

     1        flux_pt, entropy_flux_pt)

! ------------------------------

      entropy_flux_re ( i) = entropy_flux_pt
      flux_re ( i, :) = flux_pt ( :)

! ------------------------------

      END DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MHD_RIEM_HLLEM_DB

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)
      PROGRAM TEST_MHD_RIEM_HLLEM_DB

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: i, nriembegin, nriemend,
     1        double_dissipate_efz, double_dissipate_efy,
     1        ientropyfix, iprimitive, use_approx_eos,
     1        split_prsflux_re

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor,
     1     pi, fracdiff, gasconst, divratio_einfeldt_rs

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     bfxre, bfyremi, bfyrepl, bfzremi, bfzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID) ::
     1    lf_eigenvec_re, rt_eigenvec_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID) ::
     1     eigenval_re, eigenwt_re, conserved_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD + NFLUID + 1) ::
     1     flux_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD + NFLUID) ::
     1      zonebdy_flux, zonebdy_conserved

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho0, prs0, vlx0, vly0, vlz0, bfy0, bfz0,
     1     entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------

      nriembegin = 1
      nriemend = 1
      double_dissipate_efz = 0
      double_dissipate_efy = 0
      ientropyfix = 1
      iprimitive = 0
      use_approx_eos = 1
      split_prsflux_re = 0

      gamma = 1.666

      smallnum = 1.0e-12
      rhofloor = 1.0e-8
      prsfloor = 1.0e-12
      pi = 4.0 * ATAN ( 1.0)
      fracdiff = 0.1
      gasconst = 1.0
      divratio_einfeldt_rs = 0.1

! ----------------------------------------------------------------------

      DO i = nriembegin, nriemend

      rhoremi ( i) = 1.0
      prsremi ( i) = 20.0

      vlxremi ( i) = 10.0
      vlyremi ( i) = 0.0
      vlzremi ( i) = 0.0

      bfyremi ( i) = 5.0
      bfzremi ( i) = 0.0

      bfxre ( i) = 5.0

      rhorepl ( i) = 1.0
      prsrepl ( i) = 1.0

      vlxrepl ( i) = -10.0
      vlyrepl ( i) = 0.0
      vlzrepl ( i) = 0.0

      bfyrepl ( i) = 5.0
      bfzrepl ( i) = 0.0


      gamma_eos_remi ( i) = gamma
      gamma_eos_repl ( i) = gamma

      gamma_soundspeed_remi ( i) = gamma
      gamma_soundspeed_repl ( i) = gamma

      END DO

      CALL MHD_RIEM_HLLEM_DB

     1   ( nriembegin, nriemend,

     1     double_dissipate_efz, double_dissipate_efy,

     1     gamma, smallnum, prsfloor, rhofloor, pi, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        bfxre,
     1        bfyremi, bfyrepl,
     1        bfzremi, bfzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, bfy0, bfz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

      END PROGRAM TEST_MHD_RIEM_HLLEM_DB
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


