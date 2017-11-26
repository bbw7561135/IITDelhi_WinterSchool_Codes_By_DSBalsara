************************************************************************
************************************************************************

#include "directives"

      SUBROUTINE RHO_ENG_IN_PRS_T_GAMA_OUT_EULER
     1        ( oned_begin, oned_end, oned_off,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rho, eng, specie,

     1          prs, temperature, entropy,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed)

************************************************************************

* Based on the formulation of:
* J-S Shuen, M-S Liou and B. van Leer, J. Comput. Phys., vol. 90,
* pg. 371, (1990)

* INPUTS:
* gasconst : Universal gas const.
* prsfloor : Floor value for the pressure.
* cp_specie : Specific heat (per mole) at constant pressure.
* heat_of_formation : heat of formation (per mole).

* rho : density
* eng : thermal energy density
* specie : Species fraction (ranges from 0 to 1).

* OUTPUTS:
* prs : Pressure; temperature ; entropy : Acutally a pseudo-entropy
* mean_mol_wt : mean molecular weight of the gas mixture.
* gamma_eos : gas gamma based on pressure and energy density.
* gamma_soundspeed : gas gamma based on sound speed.

****************

* This subroutine would be extremely useful for initializing pressure
* and pseudo-entropy at the beginning of each fractional timestep.

****************

* Takes density, energy density and species fractions (ranging from
* 0 to 1.0) and uses the atomic properties to produce the pressure,
* temperature and pseudo-entropy.

* Also make mean molecular weight, "mean_mol_wt".

* Also makes "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?

* We assume constant "cp_specie" for each specie. For a real gas it may
* become necessary to fit polynomial functions (w.r.t. temperature)
* to "cp_specie". In that case, one will have to construct not just
* one "mean_cp" but also mean values for the higher polynomial coefficients.

* This subroutine is valid for ideal gases for which "gamma_soundspeed"
* can be calculated without an evaluation of derivatives of thermodynamic
* variables. All the details for the non-ideal case have also been
* programmed in though. They have, however, been deleted from this
* subroutine (though they are all there in the previous one.

* Written by: Dinshaw Balsara

************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: oned_begin, oned_end, oned_off

      REAL, INTENT (IN) :: gasconst, gamma,
     1      smallnum, prsfloor, rhofloor

      REAL, DIMENSION ( 0 : NFLUID_EULER), INTENT (IN) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1      INTENT (INOUT) ::
     1      rho, eng, prs, temperature, entropy,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1      INTENT (INOUT) ::
     1      specie

************************************************************************

* Variables that are internal to this subroutine.

      INTEGER i, ifluid

      REAL engfloor

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      mean_cp, mean_cv, mean_heat_of_formation,
     1      temparr_01, temparr_02

************************************************************************

      engfloor = prsfloor / ( gamma - 1.0)

* Provide a floor for the energy density.

* If multispecies is not called for, make a quick exit.

#if (NFLUID_EULER == 0)

      DO i = oned_begin - oned_off, oned_end + oned_off

      IF ( eng ( i) .GT. engfloor) THEN
      prs ( i) = eng ( i) * (gamma - 1.0)
      ELSE
      eng ( i) = engfloor
      prs ( i) = prsfloor
      END IF

      mean_mol_wt ( i) = molewt_specie ( 0)
      gamma_eos ( i) = gamma
      gamma_soundspeed ( i) = gamma

      temperature ( i) = prs ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rho ( i))

      entropy ( i) = prs ( i) / rho ( i)**gamma

      END DO

      RETURN

#endif

************************************************************************

* Build the mass fraction for the inert specie given that all the
* others are satisfied.
* Enforce physical constraint that each of the mass fractions should
* be bounded between zero and one.
* Species that serve as tags or tracers are only kept between 0 and 1 without
* further change.

      DO ifluid = 1, NFLUID_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      specie ( ifluid, i) = MAX ( specie ( ifluid, i), 0.0)
      specie ( ifluid, i) = MIN ( specie ( ifluid, i), 1.0)

      END DO

      END DO

      specie ( 0:0, :) = 1.0

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      specie ( 0, i) = specie ( 0, i) - specie ( ifluid, i)

      END DO

      END DO

* Enforce physical constraint that the sum of mass fractions for
* the active species should be bounded between zero and one.

      DO i = oned_begin - oned_off, oned_end + oned_off

      IF ( specie ( 0, i) .LT. - smallnum) THEN

      specie ( 0, i) = - specie ( 0, i) + 1.0

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      specie ( ifluid, i) = specie ( ifluid, i) / specie ( 0, i)
      specie ( ifluid, i) = MAX ( specie ( ifluid, i), 0.0)
      specie ( ifluid, i) = MIN ( specie ( ifluid, i), 1.0)

      END DO

      specie ( 0, i) = 0.0

      END IF

      specie ( 0, i) = MAX ( specie ( 0, i), 0.0)
      specie ( 0, i) = MIN ( specie ( 0, i), 1.0)

      END DO

****************

* Build the mean molecular weight in each zone.

* Also build "mean_cp" and "mean_cv" the specific heats
* at constant pressure and volume. Here we have assumed that
* "cp_specie" is a constant. If "cp_specie" has a polynomial
* dependence on temperature then the present "mean_cp" only
* gives us the constant part of the polynomial. In that case
* a similar process has to be followed for all the coefficients
* of the polynomial.

      mean_mol_wt ( :) = 1.0 / molewt_specie ( 0)
      mean_cp ( :) = cp_specie ( 0)

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_mol_wt ( i) = mean_mol_wt ( i)
     1        + specie ( ifluid, i)
     1        * ( 1.0 / molewt_specie ( ifluid)
     1          - 1.0 / molewt_specie ( 0))

      mean_cp ( i) = mean_cp ( i) + specie ( ifluid, i)
     1                  * ( cp_specie ( ifluid)
     1                    - cp_specie ( 0))

      END DO

      END DO


      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_mol_wt ( i) = 1.0 / mean_mol_wt ( i)

      mean_cv ( i) = mean_cp ( i) - gasconst / mean_mol_wt ( i)

      END DO

****************

* Find the latent heat, i.e. the heat of formation of the mixture.

      mean_heat_of_formation ( :) = 0.0

      DO ifluid = 0, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_heat_of_formation ( i) = mean_heat_of_formation ( i)
     1          + rho ( i) * specie ( ifluid, i)
     1          * heat_of_formation ( ifluid)

      END DO

      END DO

****************

* Provide a floor for the energy density.

* Use that to find the temperature and pressure.

* If "cp_specie" has a polynomial dependence on temperature then 
* one has to use a root solver.
* Also provide a floor for the pressure.

      DO i = oned_begin - oned_off, oned_end + oned_off

      eng ( i) = MAX ( eng ( i), engfloor)

      temparr_01 ( i) = eng ( i) - mean_heat_of_formation ( i)
      temparr_02 ( i) = rho ( i)
     1                * ( mean_cp ( i) - gasconst / mean_mol_wt ( i))

      temperature ( i) = temparr_01 ( i) / temparr_02 ( i)
      prs ( i) = gasconst * rho ( i) * temperature ( i)
     1         / mean_mol_wt ( i)

      IF ( eng ( i) .LE. engfloor) THEN

      prs ( i) = prsfloor
      temperature ( i) = prs ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rho ( i))
      eng ( i) = mean_heat_of_formation ( i)
     1         + temperature ( i) * temparr_02 ( i)

      END IF

      END DO

****************

* Make "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?
* Also make the pseudo-entropy.

* For ideal gases, the following is a quick way of evaluating
* "gamma_soundspeed" without evaluation of any computationally
* costly partial derivatives. As a result, the more computationally
* costly portion (which even works for non-ideal gases) is commented out.

      DO i = oned_begin - oned_off, oned_end + oned_off

      gamma_eos ( i) = 1.0 + prs ( i) / eng ( i)

      gamma_soundspeed ( i) = 1.0 + gasconst
     1                      / ( mean_mol_wt ( i) * mean_cv ( i))

      entropy ( i) = prs ( i) / rho ( i)**gamma_soundspeed ( i)

      END DO

************************************************************************

      RETURN

      END SUBROUTINE RHO_ENG_IN_PRS_T_GAMA_OUT_EULER

************************************************************************
************************************************************************

      SUBROUTINE RHO_PRS_IN_ENG_T_GAMA_OUT_EULER
     1        ( oned_begin, oned_end, oned_off,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rho, prs, specie,

     1          eng, temperature,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed)

************************************************************************

* Based on the formulation of:
* J-S Shuen, M-S Liou and B. van Leer, J. Comput. Phys., vol. 90,
* pg. 371, (1990)

* INPUTS:
* gasconst : Universal gas const.
* prsfloor : Floor value for the pressure.
* cp_specie : Specific heat (per mole) at constant pressure.
* heat_of_formation : heat of formation (per mole).

* rho : density
* prs : Pressure;
* specie : Species fraction (ranges from 0 to 1).

* OUTPUTS:
* eng : thermal energy density
* temperature ;
* mean_mol_wt : mean molecular weight of the gas mixture.
* gamma_eos : gas gamma based on pressure and energy density.
* gamma_soundspeed : gas gamma based on sound speed.

****************

* This subroutine would be extremely useful in the flux computation
* subroutines.

****************

* Takes density, energy density and species fractions (ranging from
* 0 to 1.0) and uses the atomic properties to produce the pressure,
* temperature and pseudo-entropy.

* Also make mean molecular weight, "mean_mol_wt".

* Also makes "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?

* We assume constant "cp_specie" for each specie. For a real gas it may
* become necessary to fit polynomial functions (w.r.t. temperature)
* to "cp_specie". In that case, one will have to construct not just
* one "mean_cp" but also mean values for the higher polynomial coefficients.

* This subroutine is valid for ideal gases for which "gamma_soundspeed"
* can be calculated without an evaluation of derivatives of thermodynamic
* variables. All the details for the non-ideal case have also been
* programmed in though. They have, however, been eliminated here
* thought they have been retained in the previous subroutine.

* Written by: Dinshaw Balsara

************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: oned_begin, oned_end, oned_off

      REAL, INTENT (IN) :: gasconst, gamma,
     1      smallnum, prsfloor, rhofloor

      REAL, DIMENSION ( 0 : NFLUID_EULER), INTENT (IN) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1      INTENT (INOUT) ::
     1      rho, eng, prs, temperature,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1      INTENT (INOUT) ::
     1      specie

************************************************************************

* Variables that are internal to this subroutine.

      INTEGER i, ifluid

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      mean_cp, mean_cv, mean_heat_of_formation,
     1      temparr_01, temparr_02

************************************************************************

* If multispecies is not called for, make a quick exit.

#if (NFLUID_EULER == 0)

      DO i = oned_begin - oned_off, oned_end + oned_off

      eng ( i) = prs ( i) / (gamma - 1.0)
      mean_mol_wt ( i) = molewt_specie ( 0)
      gamma_eos ( i) = gamma
      gamma_soundspeed ( i) = gamma

      temperature ( i) = prs ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rho ( i))

      END DO

      RETURN

#endif

************************************************************************

* Build the mass fraction for the inert specie given that all the
* others are satisfied.
* Enforce physical constraint that each of the mass fractions should
* be bounded between zero and one.
* Species that serve as tags or tracers are only kept between 0 and 1 without
* further change.

      DO ifluid = 1, NFLUID_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      specie ( ifluid, i) = MAX ( specie ( ifluid, i), 0.0)
      specie ( ifluid, i) = MIN ( specie ( ifluid, i), 1.0)

      END DO

      END DO

      specie ( 0:0, :) = 1.0

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      specie ( 0, i) = specie ( 0, i) - specie ( ifluid, i)

      END DO

      END DO

* Enforce physical constraint that the sum of mass fractions for
* the active species should be bounded between zero and one.

      DO i = oned_begin - oned_off, oned_end + oned_off

      IF ( specie ( 0, i) .LT. - smallnum) THEN

      specie ( 0, i) = - specie ( 0, i) + 1.0

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      specie ( ifluid, i) = specie ( ifluid, i) / specie ( 0, i)
      specie ( ifluid, i) = MAX ( specie ( ifluid, i), 0.0)
      specie ( ifluid, i) = MIN ( specie ( ifluid, i), 1.0)

      END DO

      specie ( 0, i) = 0.0

      END IF

      specie ( 0, i) = MAX ( specie ( 0, i), 0.0)
      specie ( 0, i) = MIN ( specie ( 0, i), 1.0)

      END DO

****************

* Build the mean molecular weight in each zone.

* Also build "mean_cp" and "mean_cv" the specific heats
* at constant pressure and volume. Here we have assumed that
* "cp_specie" is a constant. If "cp_specie" has a polynomial
* dependence on temperature then the present "mean_cp" only
* gives us the constant part of the polynomial. In that case
* a similar process has to be followed for all the coefficients
* of the polynomial.

      mean_mol_wt ( :) = 1.0 / molewt_specie ( 0)
      mean_cp ( :) = cp_specie ( 0)

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_mol_wt ( i) = mean_mol_wt ( i)
     1        + specie ( ifluid, i)
     1        * ( 1.0 / molewt_specie ( ifluid)
     1          - 1.0 / molewt_specie ( 0))

      mean_cp ( i) = mean_cp ( i) + specie ( ifluid, i)
     1                  * ( cp_specie ( ifluid)
     1                    - cp_specie ( 0))

      END DO

      END DO


      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_mol_wt ( i) = 1.0 / mean_mol_wt ( i)

      mean_cv ( i) = mean_cp ( i) - gasconst / mean_mol_wt ( i)

      END DO

****************

* Find the latent heat, i.e. the heat of formation of the mixture.

      mean_heat_of_formation ( :) = 0.0

      DO ifluid = 0, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_heat_of_formation ( i) = mean_heat_of_formation ( i)
     1          + rho ( i) * specie ( ifluid, i)
     1          * heat_of_formation ( ifluid)

      END DO

      END DO

****************

* Use that to find the temperature and thermal energy density.

* If "cp_specie" has a polynomial dependence on temperature then 
* one has to use a root solver.
* Also provide a floor for the pressure.

      DO i = oned_begin - oned_off, oned_end + oned_off

      temparr_02 ( i) = rho ( i)
     1                * ( mean_cp ( i) - gasconst / mean_mol_wt ( i))

      IF ( prs ( i) .LE. prsfloor) prs ( i) = prsfloor

      temperature ( i) = prs ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rho ( i))
      eng ( i) = mean_heat_of_formation ( i)
     1         + temperature ( i) * temparr_02 ( i)

      END DO

****************

* Make "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?
* Also make the pseudo-entropy.

* For ideal gases, the following is a quick way of evaluating
* "gamma_soundspeed" without evaluation of any computationally
* costly partial derivatives. As a result, the more computationally
* costly portion (which even works for non-ideal gases) is commented out.

      DO i = oned_begin - oned_off, oned_end + oned_off

      gamma_eos ( i) = 1.0 + prs ( i) / eng ( i)

      gamma_soundspeed ( i) = 1.0 + gasconst
     1                      / ( mean_mol_wt ( i) * mean_cv ( i))

*      entropy ( i) = prs ( i) / rho ( i)**gamma_soundspeed ( i)

      END DO

************************************************************************

      RETURN

      END SUBROUTINE RHO_PRS_IN_ENG_T_GAMA_OUT_EULER

************************************************************************
************************************************************************

      SUBROUTINE RHO_PRS_IN_ENG_T_GAMA_DERVS_O_EU
     1        ( oned_begin, oned_end, oned_off,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rho, prs, specie,

     1          eng, temperature,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed,

     1          dprs_drho, dprs_deng, dprs_dspecie,
     1          deng_drho, deng_dprs, deng_dspecie)

************************************************************************

* Based on the formulation of:
* J-S Shuen, M-S Liou and B. van Leer, J. Comput. Phys., vol. 90,
* pg. 371, (1990)

* INPUTS:
* gasconst : Universal gas const.
* prsfloor : Floor value for the pressure.
* cp_specie : Specific heat (per mole) at constant pressure.
* heat_of_formation : heat of formation (per mole).

* rho : density
* prs : Pressure;
* specie : Species fraction (ranges from 0 to 1).

* OUTPUTS:
* eng : thermal energy density
* temperature ;
* mean_mol_wt : mean molecular weight of the gas mixture.
* gamma_eos : gas gamma based on pressure and energy density.
* gamma_soundspeed : gas gamma based on sound speed.

* dprs_drho, dprs_deng, dprs_dspecie : Partial derivatives of the
*       pressure w.r.t. density; thermal energy density and species mass
*       concentration ( not mass fraction).
* The variables to be held fixed are two of the following three:
* density, energy density and species mass concentrations.

* deng_drho, deng_dprs, deng_dspecie : Partial derivatives of the
*       thermal energy density w.r.t. density; pressure and
*       species fractions.
* The variables to be held fixed are two of the following three:
* density, pressure and species fractions.

****************

* This subroutine would be extremely useful in the eigenvector subroutines
* when going from eigenvectors in primitive variables to eigenvectors
* in conserved variables.

****************

* Takes density, energy density and species fractions (ranging from
* 0 to 1.0) and uses the atomic properties to produce the pressure,
* temperature and pseudo-entropy.

* Also make mean molecular weight, "mean_mol_wt".

* Also makes "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?

* Also produces partial derivatives of the pressure w.r.t.: a) the density;
* b) the thermal energy density and c) The species densities.
* (Useful for transforming the left eigenvector from primitive
* to conserved variables.)

* Also produces partial derivatives of the energy density w.r.t.:
* a) the density; b) the pressure and c) The species fractions.
* (Useful for transforming the right eigenvector from primitive
* to conserved variables.)

* We assume constant "cp_specie" for each specie. For a real gas it may
* become necessary to fit polynomial functions (w.r.t. temperature)
* to "cp_specie". In that case, one will have to construct not just
* one "mean_cp" but also mean values for the higher polynomial coefficients.

* This subroutine is valid for ideal gases for which "gamma_soundspeed"
* can be calculated without an evaluation of derivatives of thermodynamic
* variables. All the details for the non-ideal case have also been
* programmed in though. They have, however, been commented out.

* For ideal gases, when derivatives are not needed, the whole
* portion that evaluates derivatives of thermodynamic variables can
* be eliminated.

* Written by: Dinshaw Balsara

************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: oned_begin, oned_end, oned_off

      REAL, INTENT (IN) :: gasconst, gamma,
     1      smallnum, prsfloor, rhofloor

      REAL, DIMENSION ( 0 : NFLUID_EULER), INTENT (IN) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1      INTENT (INOUT) ::
     1      rho, eng, prs, temperature,
     1      dprs_drho, dprs_deng, deng_drho, deng_dprs,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1      INTENT (INOUT) ::
     1      specie,
     1      dprs_dspecie, deng_dspecie

************************************************************************

* Variables that are internal to this subroutine.

      INTEGER i, ifluid

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      mean_cp, mean_cv, mean_heat_of_formation,
     1      temparr_01, temparr_02

************************************************************************

* If multispecies is not called for, make a quick exit.

#if (NFLUID_EULER == 0)

      DO i = oned_begin - oned_off, oned_end + oned_off

      eng ( i) = prs ( i) / (gamma - 1.0)
      mean_mol_wt ( i) = molewt_specie ( 0)
      gamma_eos ( i) = gamma
      gamma_soundspeed ( i) = gamma

      temperature ( i) = prs ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rho ( i))

      dprs_drho ( i) = 0.0
      dprs_deng ( i) = gamma - 1.0
      deng_drho ( i) = 0.0
      deng_dprs ( i) = 1.0 / ( gamma - 1.0)

      END DO

      RETURN

#endif

************************************************************************

* Build the mass fraction for the inert specie given that all the
* others are satisfied.
* Enforce physical constraint that each of the mass fractions should
* be bounded between zero and one.
* Species that serve as tags or tracers are only kept between 0 and 1 without
* further change.

      DO ifluid = 1, NFLUID_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      specie ( ifluid, i) = MAX ( specie ( ifluid, i), 0.0)
      specie ( ifluid, i) = MIN ( specie ( ifluid, i), 1.0)

      END DO

      END DO

      specie ( 0:0, :) = 1.0

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      specie ( 0, i) = specie ( 0, i) - specie ( ifluid, i)

      END DO

      END DO

* Enforce physical constraint that the sum of mass fractions for
* the active species should be bounded between zero and one.

      DO i = oned_begin - oned_off, oned_end + oned_off

      IF ( specie ( 0, i) .LT. - smallnum) THEN

      specie ( 0, i) = - specie ( 0, i) + 1.0

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      specie ( ifluid, i) = specie ( ifluid, i) / specie ( 0, i)
      specie ( ifluid, i) = MAX ( specie ( ifluid, i), 0.0)
      specie ( ifluid, i) = MIN ( specie ( ifluid, i), 1.0)

      END DO

      specie ( 0, i) = 0.0

      END IF

      specie ( 0, i) = MAX ( specie ( 0, i), 0.0)
      specie ( 0, i) = MIN ( specie ( 0, i), 1.0)

      END DO

****************

* Build the mean molecular weight in each zone.

* Also build "mean_cp" and "mean_cv" the specific heats
* at constant pressure and volume. Here we have assumed that
* "cp_specie" is a constant. If "cp_specie" has a polynomial
* dependence on temperature then the present "mean_cp" only
* gives us the constant part of the polynomial. In that case
* a similar process has to be followed for all the coefficients
* of the polynomial.

      mean_mol_wt ( :) = 1.0 / molewt_specie ( 0)
      mean_cp ( :) = cp_specie ( 0)

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_mol_wt ( i) = mean_mol_wt ( i)
     1        + specie ( ifluid, i)
     1        * ( 1.0 / molewt_specie ( ifluid)
     1          - 1.0 / molewt_specie ( 0))

      mean_cp ( i) = mean_cp ( i) + specie ( ifluid, i)
     1                  * ( cp_specie ( ifluid)
     1                    - cp_specie ( 0))

      END DO

      END DO


      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_mol_wt ( i) = 1.0 / mean_mol_wt ( i)

      mean_cv ( i) = mean_cp ( i) - gasconst / mean_mol_wt ( i)

      END DO

****************

* Find the latent heat, i.e. the heat of formation of the mixture.

      mean_heat_of_formation ( :) = 0.0

      DO ifluid = 0, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      mean_heat_of_formation ( i) = mean_heat_of_formation ( i)
     1          + rho ( i) * specie ( ifluid, i)
     1          * heat_of_formation ( ifluid)

      END DO

      END DO

****************

* Use that to find the temperature and thermal energy density.

* If "cp_specie" has a polynomial dependence on temperature then 
* one has to use a root solver.
* Also provide a floor for the pressure.

      DO i = oned_begin - oned_off, oned_end + oned_off

      temparr_02 ( i) = rho ( i)
     1                * ( mean_cp ( i) - gasconst / mean_mol_wt ( i))

      IF ( prs ( i) .LE. prsfloor) prs ( i) = prsfloor

      temperature ( i) = prs ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rho ( i))
      eng ( i) = mean_heat_of_formation ( i)
     1         + temperature ( i) * temparr_02 ( i)

      END DO

****************

* Make "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?
* Also make the pseudo-entropy.

* For ideal gases, the following is a quick way of evaluating
* "gamma_soundspeed" without evaluation of any computationally
* costly partial derivatives. As a result, the more computationally
* costly portion (which even works for non-ideal gases) is commented out.

      DO i = oned_begin - oned_off, oned_end + oned_off

      gamma_eos ( i) = 1.0 + prs ( i) / eng ( i)

      gamma_soundspeed ( i) = 1.0 + gasconst
     1                      / ( mean_mol_wt ( i) * mean_cv ( i))

*      entropy ( i) = prs ( i) / rho ( i)**gamma_soundspeed ( i)

      END DO

****************

* Calculate the derivatives of pressure w.r.t. a) the density;
* b) the thermal energy density and c) The species densities.

* Also calculate partial derivatives of the energy density w.r.t.:
* a) the density; b) the pressure and c) The species fractions.

      DO i = oned_begin - oned_off, oned_end + oned_off

      temparr_01 ( i) = gasconst / mean_mol_wt ( i)

      dprs_drho ( i) = temparr_01 ( i)
     1               * ( temperature ( i)
     1                 - eng ( i) / ( rho ( i) * mean_cv ( i)) )

      dprs_deng ( i) = temparr_01 ( i) / mean_cv ( i)

      deng_drho ( i) = eng ( i) / rho ( i)
     1               - prs ( i) * mean_cv ( i)
     1               / ( temparr_01 ( i) * rho ( i))

      deng_dprs ( i) = mean_cv ( i) / temparr_01 ( i)

      END DO


* If "cp_specie" is a polynomial function of temperature then we have
* to integrate it.

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      temparr_01 ( i) = cp_specie ( ifluid) * temperature ( i)
     1                + heat_of_formation ( ifluid)
     1                - gasconst * temperature ( i)
     1                / molewt_specie ( ifluid)

      temparr_02 ( i) = cp_specie ( 0) * temperature ( i)
     1                + heat_of_formation ( 0)
     1                - gasconst * temperature ( i)
     1                / molewt_specie ( 0)

      dprs_dspecie ( ifluid, i) = temperature ( i)
     1                * ( 1.0 / molewt_specie ( ifluid)
     1                  - 1.0 / molewt_specie ( 0))
     1                - ( temparr_01 ( i) - temparr_02 ( i))
     1                / ( mean_mol_wt ( i) * mean_cv ( i))

      dprs_dspecie ( ifluid, i) = dprs_dspecie ( ifluid, i)
     1                          * gasconst * rho ( i)

      deng_dspecie ( ifluid, i) =
     1                  ( temparr_01 ( i) - temparr_02 ( i))
     1                - mean_mol_wt ( i) * mean_cv ( i)
     1                * temperature ( i)
     1                * ( 1.0 / molewt_specie ( ifluid)
     1                  - 1.0 / molewt_specie ( 0))

      deng_dspecie ( ifluid, i) = rho ( i)
     1                          * deng_dspecie ( ifluid, i)

      END DO

      END DO

****************

* Make "gamma_eos, gamma_soundspeed". These say that if the gas were
* a gamma law gas, what is the effective gamma we would deduce from
* either the equation of state or from the sound speed?
* Also make the pseudo-entropy.

* This is a more computationally costly way of evaluating
* "gamma_soundspeed" which is even valid for non-ideal gases.
* For ideal gases, which we do here, it is not needed and, therefore,
* commented out.

*      DO i = oned_begin - oned_off, oned_end + oned_off
*
*      gamma_eos ( i) = 1.0 + prs ( i) / eng ( i)
*
*      gamma_soundspeed ( i) = ( eng ( i) + prs ( i)) * dprs_deng ( i)
*     1                      + dprs_drho ( i)
*      gamma_soundspeed ( i) = gamma_soundspeed ( i) / prs ( i)
*
*      entropy ( i) = prs ( i) / rho ( i)**gamma_soundspeed ( i)
*
*      END DO

****************

* Go from species fractions to species densities as the independent
* variables for the pressure derivatives. This involves adding new terms
* to "dprs_drho".
* Note: Once that is done, "dprs_drho" can't any more be used
* to obtain the soundspeed. To obtain the sound speed we can
* only use "gamma_soundspeed".

      DO ifluid = 1, NFLUID_EULER - NTRACERS_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      dprs_drho ( i) = dprs_drho ( i) 
     1               - dprs_dspecie ( ifluid, i)
     1               * specie ( ifluid, i) / rho ( i)

      dprs_dspecie ( ifluid, i) = dprs_dspecie ( ifluid, i) 
     1                          / rho ( i)

      END DO

      END DO

****************

* Compensate for the species fractions that are just tracers.

      DO ifluid = NFLUID_EULER - NTRACERS_EULER + 1, NFLUID_EULER

      DO i = oned_begin - oned_off, oned_end + oned_off

      dprs_dspecie ( ifluid, i) = 0.0
      deng_dspecie ( ifluid, i) = 0.0

      END DO

      END DO

************************************************************************

      RETURN

      END SUBROUTINE RHO_PRS_IN_ENG_T_GAMA_DERVS_O_EU

************************************************************************
************************************************************************

#if ( 0 == 1)

      PROGRAM TESTEOS

* cpp -traditional eos.src eos.f
* ifc -autodouble -w -xW -O2 -o xeos eos.f

************************************************************************
************************************************************************

      IMPLICIT NONE

      INTEGER oned_begin, oned_end, oned_off, i, ifluid

      REAL gasconst, gamma, smallnum, prsfloor, rhofloor

      REAL, DIMENSION ( 0 : NFLUID_EULER) ::
     1      cp_specie, molewt_specie, heat_of_formation,
     1      cv_specie, gamma_specie

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      rho, eng, prs, temperature, entropy,
     1      dprs_drho, dprs_deng, deng_drho, deng_dprs,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie,
     1      dprs_dspecie, deng_dspecie

************************************************************************

      OPEN ( unit = 6, file = "testeos.out")


      oned_begin = ONEDFLUIDBEGIN
      oned_end = ONEDFLUIDEND
      oned_off = 0

      gamma = 5.0 / 3.0
      gasconst = 1.0
      smallnum = 1.0e-12
      prsfloor = 1.0e-16
      rhofloor = 1.0e-8


      DO ifluid = 0, NFLUID_EULER

      gamma_specie ( ifluid) = gamma

      molewt_specie ( ifluid) = 1.0

      cp_specie ( ifluid) = ( gamma_specie ( ifluid)
     1                      / ( gamma_specie ( ifluid) - 1.0) )
     1                    * gasconst
     1                    / molewt_specie ( ifluid)

      heat_of_formation ( ifluid) = 0.0

      cv_specie ( ifluid) = cp_specie ( ifluid)
     1                    / gamma_specie ( ifluid)

      END DO


      DO i = oned_begin - oned_off, oned_end + oned_off

      rho ( i) = 1.0
      prs ( i) = 1.0
      eng ( i) = prs ( i) / ( gamma - 1.0)
      specie ( :, i) = 1.0 / AMAX1 ( 1.0, REAL ( NFLUID_EULER))

      END DO


      CALL RHO_ENG_IN_PRS_T_GAMA_DERVS_OUT
     1        ( oned_begin, oned_end, oned_off,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rho, eng, specie,

     1          prs, temperature, entropy,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed,

     1          dprs_drho, dprs_deng, dprs_dspecie,
     1          deng_drho, deng_dprs, deng_dspecie)


      DO i = oned_begin - oned_off, oned_end + oned_off

      WRITE ( 6, *)"****************************************"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"i = ", i
      WRITE ( 6, *)"rho, eng, specie = ", rho ( i), eng ( i),
     1              specie ( 1, i)
      WRITE ( 6, *)"prs, temperature, entropy = ", prs ( i),
     1              temperature ( i), entropy ( i)
      WRITE ( 6, *)"mean_mol_wt, gamma_eos, gamma_soundspeed = ",
     1          mean_mol_wt ( i), gamma_eos ( i), gamma_soundspeed ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"dprs_drho, dprs_deng = ",
     1          dprs_drho ( i), dprs_deng ( i)
      WRITE ( 6, *)"dprs_dspecie = ",
     1  ( dprs_dspecie ( ifluid, i), ifluid = 1, NFLUID_EULER)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"deng_drho, deng_dprs = ",
     1          deng_drho ( i), deng_dprs ( i)
      WRITE ( 6, *)"deng_dspecie = ",
     1  ( deng_dspecie ( ifluid, i), ifluid = 1, NFLUID_EULER)
      WRITE ( 6, *)" "

      END DO

************************************************************************
************************************************************************

      STOP

      END

#endif

************************************************************************
************************************************************************






