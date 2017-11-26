! ----------------------------------------------------------------------

      PROGRAM EULER_RK_APPROX_RS

! Two-stage Runge-Kutta scheme for updating Euler equations with various
! approximate Riemann solvers incorporated in the same package. This enables
! one to intercompare Riemann solvers. Linearized, HLL and HLLC and HLLEM Riemann
! solvers for the Euler equations included here. The linearized Riemann
! solver includes an HLL fixup.

! GNU Fortran compilation lines:-
! cp euler_RK_approxRS.src euler_RK_approxRS.F
! gfortran -fdefault-integer-8 -fdefault-real-8 -O2 -o xeuler_1d euler_RK_approxRS.F

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix, ifield, j, k, islope, icharlim, iriem_solvr,
     1        nzones, ntstep, nsys, nriembegin, nriemend,
     1        ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,
     1        istep, iproblem, ifracstep, iprimitive, rieiter,
     1        ientropyfix, use_approx_eos, split_prsflux_re, debug_stop

      REAL mu, t, dtcur, dtprev, dx, tstop, tempa, div_v_coef, lap_visc,
     1     prsfac, rhofloor, prsfloor, smallnum, gamma, gasconst,
     1     fracdiff, divratio_einfeldt_rs

      PARAMETER ( nzones = 200, NFIELD = 5, NFLUID = 0,
     1            nsys = NFIELD + NFLUID,
     1            ONEDFLUIDBEGIN = 1 - 3, ONEDFLUIDEND = nzones + 3)

      REAL, DIMENSION ( 0 : NFLUID) :: cv_specie, cp_specie,
     1      molewt_specie, heat_of_formation

! Big data variables  used for RK timestep.

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND ) :: divflag

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND ) :: x,
     1      sound_spd, visc_coef

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, nsys) ::
     1      u, u0, dudt, v, del_v, del_eigwt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, nsys, -2: 2) ::
     1      eigwt

! ----------------------------------------------------------------------

! Big data variables passed through eigenvector subroutine interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, nsys) ::
     1     eigenval

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, nsys, nsys) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

! Big data variables passed through Riemann Solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      rhoremi, prsremi, vlxremi, vlyremi, vlzremi,
     1      rhorepl, prsrepl, vlxrepl, vlyrepl, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr,

     1     rhoflux, engflux, prsflux, vlxflux, vlyflux, vlzflux

      REAL, DIMENSION ( 0 : NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD + NFLUID, NFIELD + NFLUID) ::
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
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION
     1    ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------

! Declare some info. for TVD limiters.

      REAL MC_BETA_LIMITER, MIN_MOD, a, b, mcb

      MIN_MOD ( a, b) = AMIN1 ( ABS (a), ABS (b))
     1                  * ( SIGN ( 0.5, a) + SIGN ( 0.5, b))

      MC_BETA_LIMITER ( a, b, mcb) =
     1              ( SIGN ( 0.5, a) + SIGN ( 0.5, b))
     1            * AMIN1 ( 0.5 * ( ABS ( a + b)),
     1                      mcb * ABS ( a), mcb * ABS ( b) )

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "euler_RK_approxRS.out")

! iproblem = 1, 2, 3, 4, 5 for Sod, Lax, WC, Noh, rarefaction
      iproblem = 5

! islope controls how the slopes are built. islope = 1,2 for min-mod & MC-beta.
      islope = 2

! Decide on whether we want to reconstruct characteristically or componentwise;
! i.e. set to 1 or 0 respectively.
      icharlim = 1

! Choose Riemann solver; iriem_solvr == 1, 2, 3, 4 for linearized, HLLC, HLL, HLLEM RS.
      iriem_solvr = 4

      mu = 0.8            ! CFL number

      dx = 1.0 / REAL ( nzones)   ! Zone size

      t = 0.0             ! Set initial time to zero.

      IF ( iproblem .EQ. 1) THEN
        tstop = 0.2    ! Stopping time.
        gamma = 1.4    ! Set polytropic index
      ELSE IF ( iproblem .EQ. 2) THEN
        tstop = 0.13
        gamma = 1.4
      ELSE IF ( iproblem .EQ. 3) THEN
        tstop = 0.01
        gamma = 1.4
      ELSE IF ( iproblem .EQ. 4) THEN
        tstop = 0.4
        gamma = 5.0 / 3.0
      ELSE IF ( iproblem .EQ. 5) THEN
        tstop = 0.08
        gamma = 5.0 / 3.0
      ELSE
        STOP'bad iproblem'
      END IF

! gas constant can be scaled to any number as long as temperature is not used.
      gasconst = 1.0

      ntstep = 35000 ! max. # of timesteps.

      mcb = 1.5    ! ranges between 1 and 2.

! If the divergence of the velocity at a zone boundary is negative and exceeds
! this coef times the sound speed then we have to apply more stringent limiter
! to that zone.
      div_v_coef = 0.4

! If a strong shock is indeed detected then put in a small amount of
! artificial viscosity at that zone boundary.
      lap_visc = 0.1

      rhofloor = 1.0e-8   ! Set density floor

      prsfloor = 1.0e-12 ! Set pressure floor

      smallnum = 1.0e-12  ! Set a small number for the code.

! Riemann Solver variables.

      fracdiff = 0.1  !Set fractional change needed for entropy fix.(deprecated)

      ientropyfix = 1 ! Set choice of entropy fix 1,2 or 3.

      use_approx_eos = 1 ! Use an approximate EOS.

      split_prsflux_re = 0 ! Don't split momentum flux into pressure + other.

! For divergences that are larger than this factor, use HLL.
      divratio_einfeldt_rs = 0.4

! Any values for specific heats will do for polytropic gas as long as the
! ratio is preserved.
      cv_specie ( :) = 1.0
      cp_specie ( :) = gamma * cv_specie ( :)

      molewt_specie ( :) = 1.0
      heat_of_formation ( :) = 0.0

      mean_mol_wt_remi ( :) = 1.0
      gamma_eos_remi ( :) = gamma
      gamma_soundspeed_remi ( :) = gamma

      mean_mol_wt_repl ( :) = 1.0
      gamma_eos_repl ( :) = gamma
      gamma_soundspeed_repl ( :) = gamma

      debug_stop = 0

! ------------------------------

      WRITE ( 6, *)"iproblem = ", iproblem
      WRITE ( 6, *)" "
      WRITE ( 6, *)"islope, icharlim, mcb = ", islope, icharlim, mcb
      WRITE ( 6, *)" "
      WRITE ( 6, *)"div_v_coef, lap_visc = ", div_v_coef, lap_visc
      WRITE ( 6, *)" "
      WRITE ( 6, *)"mu, dx, tstop, ntstep = ", mu, dx, tstop, ntstep
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma, rhofloor, prsfloor, smallnum = ",
     1              gamma, rhofloor, prsfloor, smallnum
      WRITE ( 6, *)" "
      WRITE ( 6, *)"prsfac, rieiter = ", prsfac, rieiter
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

! Set up coordinate value at zone centers.

      DO ix = ONEDFLUIDBEGIN, ONEDFLUIDEND
      x ( ix) = - 0.5 + ( ix - 0.5) * dx
      END DO

! Set up primitive variables : v(ix, 1) = density; v(ix, 2) = x-velocity;
! v(ix, 3) = y-velocity; v(ix, 4) = y-velocity; v(ix, 5) = pressure.
! Do this for all the problems, "iproblem" , being set up here.
! Also set up timestep "dtcur".

      dtcur = 1.0e30

      DO ix = ONEDFLUIDBEGIN, ONEDFLUIDEND

! ----------

        IF ( iproblem .EQ. 1) THEN

! Sod problem.

! ----------

           IF ( x ( ix) .LT. 0.0) THEN
             v ( ix, 1) = 1.0
             v ( ix, 2) = 0.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 1.0
           ELSE
             v ( ix, 1) = 0.125
             v ( ix, 2) = 0.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 0.1
           END IF

! ----------

        ELSE IF ( iproblem .EQ. 2) THEN

! Lax problem.

! ----------

           IF ( x ( ix) .LT. 0.0) THEN
             v ( ix, 1) = 0.445
             v ( ix, 2) = 0.698
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 3.528
           ELSE
             v ( ix, 1) = 0.5
             v ( ix, 2) = 0.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 0.571
           END IF

! ----------

        ELSE IF ( iproblem .EQ. 3) THEN

! Half of the Woodward & Colella (1984) problem.

! ----------

           IF ( x ( ix) .LT. 0.0) THEN
             v ( ix, 1) = 1.0
             v ( ix, 2) = 0.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 0.01
           ELSE
             v ( ix, 1) = 1.0
             v ( ix, 2) = 0.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 1000.0
           END IF

! ----------

        ELSE IF ( iproblem .EQ. 4) THEN

! Noh-like problem with colliding streams.

! ----------

           IF ( x ( ix) .LT. 0.0) THEN
             v ( ix, 1) = 1.0
             v ( ix, 2) = 2.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 0.2
           ELSE
             v ( ix, 1) = 1.5
             v ( ix, 2) = - 2.0
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 0.2
           END IF

! ----------

        ELSE IF ( iproblem .EQ. 5) THEN

! Two strong rarefaction fans.

! ----------

           IF ( x ( ix) .LT. 0.0) THEN
             v ( ix, 1) = 1.0
             v ( ix, 2) = - 2.5
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 2.0
           ELSE
             v ( ix, 1) = 1.5
             v ( ix, 2) = 2.5
             v ( ix, 3) = 0.0
             v ( ix, 4) = 0.0
             v ( ix, 5) = 4.0
           END IF

! ----------

        ELSE

        STOP'bad iproblem'

        END IF

! ----------

! Make "dtcur" here.

      tempa = ABS ( v ( ix, 2))
     1      + SQRT ( gamma * v ( ix, 5) / v ( ix, 1))

      dtcur = AMIN1 ( dtcur, mu * dx / tempa)

! ----------

      END DO ! End of "DO ix ="


      dtprev = dtcur
      WRITE ( 6, *)"dtcur = ", dtcur
      WRITE ( 6, *)" "


! ------------------------------

! Obtain conserved variables from primitive variables: u(ix, 1) = density;
! u(ix, 2) = x-momentum density; u(ix, 3) = y-momentum density;
! u(ix, 4) = z-momentum density; u(ix, 5) = total energy density.

      DO ix = ONEDFLUIDBEGIN, ONEDFLUIDEND

      u ( ix, 1) =  v ( ix, 1)

      u ( ix, 2) = v ( ix, 1) * v ( ix, 2)
      u ( ix, 3) = v ( ix, 1) * v ( ix, 3)
      u ( ix, 4) = v ( ix, 1) * v ( ix, 4)

      u ( ix, 5) = 0.5 * v ( ix, 1) * ( v ( ix, 2)**2 + v ( ix, 3)**2
     1                                + v ( ix, 4)**2)
     1           + v ( ix, 5) / ( gamma - 1.0)

      END DO ! End of "DO ix ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Take "ntstep" timesteps. Main time-integration loop begins here.

      DO istep = 1, ntstep

      DO ifracstep = 1, 2

      WRITE ( 6, *)"istep, ifracstep, t, dtcur = ",
     1              istep, ifracstep, t, dtcur
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

! Obtain primitive variables from current conserved variables for this
! fractional timestep.

      DO ix = ONEDFLUIDBEGIN, ONEDFLUIDEND

      v ( ix, 1) = AMAX1 ( u ( ix, 1), rhofloor)

      v ( ix, 2) = u ( ix, 2) / v ( ix, 1)
      v ( ix, 3) = u ( ix, 3) / v ( ix, 1)
      v ( ix, 4) = u ( ix, 4) / v ( ix, 1)

      tempa = 0.5 * v ( ix, 1) * ( v ( ix, 2)**2 + v ( ix, 3)**2
     1                           + v ( ix, 4)**2)

      v ( ix, 5) = ( u ( ix, 5) - tempa) * ( gamma - 1.0)
      v ( ix, 5) = AMAX1 ( v ( ix, 5), prsfloor)

      sound_spd ( ix) = SQRT ( gamma * v ( ix, 5) / v ( ix, 1))

      IF ( v ( ix, 5) .LE. prsfloor) u ( ix, 5) = tempa
     1                             + prsfloor / ( gamma - 1.0)

      END DO ! End of "DO ix ="

! ------------------------------

! Set a face-centered flag that tracks regions of strong shock.
! Also make the artificial viscosity coefficient.
! Also save divergence of velocity at zone boundaries.

      DO ix = ONEDFLUIDBEGIN + 1, ONEDFLUIDEND - 1

      tempa = v ( ix + 1, 2) - v ( ix, 2)

      divvel_arr ( ix) = ABS ( tempa)
      msonic_arr ( ix) = AMIN1 ( sound_spd ( ix), sound_spd ( ix + 1))

      IF ( tempa .LT. - div_v_coef
     1          * AMIN1 ( sound_spd ( ix), sound_spd ( ix + 1)) ) THEN

        divflag ( ix) = 1
        visc_coef ( ix) = lap_visc * ABS ( tempa)

      ELSE

        divflag ( ix) = 0
        visc_coef ( ix) = 0.0

      END IF

      END DO ! End of "DO ix ="

        divflag ( ONEDFLUIDBEGIN) = 0
        divflag ( ONEDFLUIDEND) = 0

! ----------------------------------------------------------------------

! Apply limiters to obtain "del_v". Do this in the space of primitive
! variables or in the space of characteristic variables depending on
! "icharlim". Choice of limiters according to "islope".

      IF ( icharlim .EQ. 0) THEN

! Limit the variables in physical space.

! ------------------------------

        IF ( islope .EQ. 1) THEN

! MinMod limiter in physical space.

        DO k = 1, nsys

          DO ix = ONEDFLUIDBEGIN + 1, ONEDFLUIDEND - 1

          del_v ( ix, k) = MIN_MOD ( v ( ix + 1, k) - v ( ix, k),
     1                               v ( ix, k) - v ( ix - 1, k) )

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

        ELSE

! MC-beta limiter in physical space.

        DO k = 1, nsys

          DO ix = ONEDFLUIDBEGIN + 1, ONEDFLUIDEND - 1

          IF ( ( divflag ( ix) .EQ. 1) .OR.
     1         ( divflag ( ix - 1) .EQ. 1) ) THEN

          del_v ( ix, k) = MIN_MOD ( v ( ix + 1, k) - v ( ix, k),
     1                               v ( ix, k) - v ( ix - 1, k) )

          ELSE

          del_v ( ix, k)=MC_BETA_LIMITER ( v ( ix + 1, k) - v ( ix, k),
     1                               v ( ix, k) - v ( ix - 1, k), mcb )

          END IF

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

        END IF ! End of "IF ( islope .EQ. 1)".

! ------------------------------

      ELSE

! Limit the variables in characteristic space.

! ------------------------------

! Load up variables that are needed for eigenvector evaluation. Then call the
! eigenvector routine.

      nriembegin = ONEDFLUIDBEGIN
      nriemend = ONEDFLUIDEND
      iprimitive = 1

      rhobar ( :) = v ( :, 1)
      vlxbar ( :) = v ( :, 2)
      vlybar ( :) = v ( :, 3)
      vlzbar ( :) = v ( :, 4)
      prsbar ( :) = v ( :, 5)

      CALL EIGENVEC_EULER

     1             ( nriembegin, nriemend, iprimitive,
     1               ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1               gamma, smallnum, rhofloor, prsfloor,

     1               gasconst,
     1               cp_specie, molewt_specie, heat_of_formation,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

! ------------------------------

! Now project into the space of eigenweights.

      eigwt ( :, :, :) = 0.0

      DO j = - 2, 2

      DO k = 1, nsys

        DO ifield = 1, nsys

          DO ix = ONEDFLUIDBEGIN + 2, ONEDFLUIDEND - 2

          eigwt ( ix, ifield, j) = eigwt ( ix, ifield, j)
     1  + lf_eigenvec ( ix, ifield, k) * v ( ix + j, k)

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

      END DO ! End of "DO k ="

      END DO ! End of "DO j ="

! ------------------------------


        IF ( islope .EQ. 1) THEN

! MinMod limiter in characteristic space.

        DO k = 1, nsys

          DO ix = ONEDFLUIDBEGIN + 2, ONEDFLUIDEND - 2

          del_eigwt ( ix, k) =
     1    MIN_MOD ( eigwt ( ix, k, 1) - eigwt ( ix, k, 0),
     1              eigwt ( ix, k, 0) - eigwt ( ix, k, -1) )

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

        ELSE

! MC-beta limiter in characteristic space.

        DO k = 1, nsys

          DO ix = ONEDFLUIDBEGIN + 1, ONEDFLUIDEND - 1

          IF ( ( divflag ( ix) .EQ. 1) .OR.
     1         ( divflag ( ix - 1) .EQ. 1) ) THEN

          del_eigwt ( ix, k) =
     1    MIN_MOD ( eigwt ( ix, k, 1) - eigwt ( ix, k, 0),
     1              eigwt ( ix, k, 0) - eigwt ( ix, k, -1) )

          ELSE

          del_eigwt ( ix, k) =
     1    MC_BETA_LIMITER ( eigwt ( ix, k, 1) - eigwt ( ix, k, 0),
     1              eigwt ( ix, k, 0) - eigwt ( ix, k, -1), mcb )

          END IF

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

        END IF ! End of "IF ( islope .EQ. 1)".

! ------------------------------

! Now project back into physical space.

      del_v ( :, :) = 0.0

      DO j = 1, nsys

        DO ifield = 1, nsys

          DO ix = ONEDFLUIDBEGIN + 2, ONEDFLUIDEND - 2

          del_v ( ix, j) = del_v ( ix, j) + del_eigwt ( ix, ifield)
     1                               * rt_eigenvec ( ix, j, ifield)

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

      END DO ! End of "DO j ="

! ------------------------------

      END IF ! End of "IF ( icharlim .EQ. 0)".

! ----------------------------------------------------------------------

! Use "v" and "del_v" to obtain left and right states at the zone boundaries.

      DO ix = 1 - 1, nzones + 1

      rhoremi ( ix) = v ( ix, 1) + 0.5 * del_v ( ix, 1)
      vlxremi ( ix) = v ( ix, 2) + 0.5 * del_v ( ix, 2)
      vlyremi ( ix) = v ( ix, 3) + 0.5 * del_v ( ix, 3)
      vlzremi ( ix) = v ( ix, 4) + 0.5 * del_v ( ix, 4)
      prsremi ( ix) = v ( ix, 5) + 0.5 * del_v ( ix, 5)

      rhorepl ( ix - 1) = v ( ix, 1) - 0.5 * del_v ( ix, 1)
      vlxrepl ( ix - 1) = v ( ix, 2) - 0.5 * del_v ( ix, 2)
      vlyrepl ( ix - 1) = v ( ix, 3) - 0.5 * del_v ( ix, 3)
      vlzrepl ( ix - 1) = v ( ix, 4) - 0.5 * del_v ( ix, 4)
      prsrepl ( ix - 1) = v ( ix, 5) - 0.5 * del_v ( ix, 5)

      END DO ! End of "DO ix =".

! ----------------------------------------------------------------------

! Invoke the Riemann solver. Uncomment one of the three choices, as desired.
! There are no multispecies eqns here.
! Choose Riemann solver; iriem_solvr == 1, 2, 3 for linearized, HLLC, HLL

      nriembegin = 0
      nriemend = nzones
      iprimitive = 0

      IF ( iriem_solvr == 1) THEN

      CALL EULER_RIEM_ROE

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

      ELSE IF ( iriem_solvr == 2) THEN

      CALL EULER_RIEM_HLLC

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

      ELSE IF ( iriem_solvr == 3) THEN

      CALL EULER_RIEM_HLL

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

      ELSE

      CALL EULER_RIEM_HLLEM_DB

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

      END IF

! ----------------------------------------------------------------------

! Unpack "flux_re". Add in an extra Lapidus viscosity if it is warranted.

      DO ix = 0, nzones

      rhoflux ( ix) = flux_re ( ix, 1)

      vlxflux ( ix) = flux_re ( ix, 2)
      vlyflux ( ix) = flux_re ( ix, 3)
      vlzflux ( ix) = flux_re ( ix, 4)

      engflux ( ix) = flux_re ( ix, 5)

      prsflux ( ix) = flux_re ( ix, 6)


      IF ( divflag ( ix) .EQ. 1) THEN

        rhoflux ( ix) = rhoflux ( ix) - visc_coef ( ix)
     1                                * ( u ( ix + 1, 1) - u ( ix, 1) )

        vlxflux ( ix) = vlxflux ( ix) - visc_coef ( ix)
     1                                * ( u ( ix + 1, 2) - u ( ix, 2) )

        vlyflux ( ix) = vlyflux ( ix) - visc_coef ( ix)
     1                                * ( u ( ix + 1, 3) - u ( ix, 3) )

        vlzflux ( ix) = vlzflux ( ix) - visc_coef ( ix)
     1                                * ( u ( ix + 1, 4) - u ( ix, 4) )

        engflux ( ix) = engflux ( ix) - visc_coef ( ix)
     1                                * ( u ( ix + 1, 5) - u ( ix, 5) )

      END IF

      END DO ! End of "DO ix =".

! ----------------------------------------------------------------------

! Use the fluxes to obtain "dudt", the time rate of update.

      DO ix = 1, nzones

      dudt ( ix, 1) = - ( rhoflux ( ix) - rhoflux ( ix - 1)) / dx

      dudt ( ix, 2) = - ( vlxflux ( ix) - vlxflux ( ix - 1)) / dx
     1                - ( prsflux ( ix) - prsflux ( ix - 1)) / dx

      dudt ( ix, 3) = - ( vlyflux ( ix) - vlyflux ( ix - 1)) / dx

      dudt ( ix, 4) = - ( vlzflux ( ix) - vlzflux ( ix - 1)) / dx

      dudt ( ix, 5) = - ( engflux ( ix) - engflux ( ix - 1)) / dx

      END DO ! End of "DO ix =".

! ----------------------------------------------------------------------

! Make the predictor or corrector step updates here.

      IF ( ifracstep .EQ. 1) THEN

        DO ix = 1, nzones

        u0 ( ix, :) = u ( ix, :)

        u ( ix, :) = u ( ix, :) + 0.5 * dtcur * dudt ( ix, :)

        END DO ! End of "DO ix =".

      ELSE

        DO ix = 1, nzones

        u ( ix, :) = u0 ( ix, :) + dtcur * dudt ( ix, :)

        END DO ! End of "DO ix =".

      END IF ! End of "IF ( ifracstep .EQ. 1)".

! ----------------------------------------------------------------------

! Obtain primitive variables from recently updated conserved variables
! and use that to obtain the new timestep "dtcur".

      IF ( ifracstep .EQ. 2) THEN

      t = t + dtcur

      dtprev = dtcur
      dtcur = 1.0e30

      DO ix = ONEDFLUIDBEGIN, ONEDFLUIDEND

! ----------

      v ( ix, 1) = AMAX1 ( u ( ix, 1), rhofloor)

      v ( ix, 2) = u ( ix, 2) / v ( ix, 1)
      v ( ix, 3) = u ( ix, 3) / v ( ix, 1)
      v ( ix, 4) = u ( ix, 4) / v ( ix, 1)

      tempa = 0.5 * v ( ix, 1) * ( v ( ix, 2)**2 + v ( ix, 3)**2
     1                           + v ( ix, 4)**2)

      v ( ix, 5) = ( u ( ix, 5) - tempa) * ( gamma - 1.0)
      v ( ix, 5) = AMAX1 ( v ( ix, 5), prsfloor)

      IF ( v ( ix, 5) .LE. prsfloor) u ( ix, 5) = tempa
     1                             + prsfloor / ( gamma - 1.0)

! ----------

! Make "dtcur" here.

      tempa = ABS ( v ( ix, 2))
     1      + SQRT ( gamma * v ( ix, 5) / v ( ix, 1))

      dtcur = AMIN1 ( dtcur, mu * dx / tempa)

! ----------

      END DO ! End of "DO ix ="

      IF ( dtcur .GT. 1.1 * dtprev) dtcur = 1.1 * dtprev

      END IF ! End of "IF ( ifracstep .EQ. 2)" where "dtcur" is evaluated.

! ----------------------------------------------------------------------

! Write out final data for imaging and wrap up the timestep loop if
! the final time has been reached.

      IF ( ( ifracstep .EQ. 2) .AND. ( t .GE. tstop) ) THEN

      OPEN ( unit = 11, file = "rh010001")
#if ( 0 == 1)
      DO ix = 1, nzones
      WRITE ( 11, *) v ( ix, 1)
      END DO
      DO ix = 1, nzones
      WRITE ( 11, *) x ( ix)
      END DO
#else
      WRITE ( 11, *) ( v ( ix, 1), ix = 1, nzones)
#endif
      CLOSE ( 11)

      OPEN ( unit = 11, file = "vx010001")
#if ( 0 == 1)
      DO ix = 1, nzones
      WRITE ( 11, *) v ( ix, 2)
      END DO
      DO ix = 1, nzones
      WRITE ( 11, *) x ( ix)
      END DO
#else
      WRITE ( 11, *) ( v ( ix, 2), ix = 1, nzones)
#endif
      CLOSE ( 11)

      OPEN ( unit = 11, file = "vy010001")
#if ( 0 == 1)
      DO ix = 1, nzones
      WRITE ( 11, *) v ( ix, 3)
      END DO
      DO ix = 1, nzones
      WRITE ( 11, *) x ( ix)
      END DO
#else
      WRITE ( 11, *) ( v ( ix, 3), ix = 1, nzones)
#endif
      CLOSE ( 11)

      OPEN ( unit = 11, file = "vz010001")
#if ( 0 == 1)
      DO ix = 1, nzones
      WRITE ( 11, *) v ( ix, 4)
      END DO
      DO ix = 1, nzones
      WRITE ( 11, *) x ( ix)
      END DO
#else
      WRITE ( 11, *) ( v ( ix, 4), ix = 1, nzones)
#endif
      CLOSE ( 11)

      OPEN ( unit = 11, file = "pr010001")
#if ( 0 == 1)
      DO ix = 1, nzones
      WRITE ( 11, *) v ( ix, 5)
      END DO
      DO ix = 1, nzones
      WRITE ( 11, *) x ( ix)
      END DO
#else
      WRITE ( 11, *) ( v ( ix, 5), ix = 1, nzones)
#endif
      CLOSE ( 11)

      GO TO 100

      END IF ! End of "IF ( ( ifracstep .EQ. 2) .AND ( t .GT. tstop) )".

! ----------------------------------------------------------------------

      END DO ! End of "DO ifracstep =".

      END DO ! End of "DO istep =".

100   CONTINUE

! ----------------------------------------------------------------------

      STOP

      END PROGRAM EULER_RK_APPROX_RS

! ----------------------------------------------------------------------


*******************************************************************
*******************************************************************

      SUBROUTINE EIGENVEC_EULER

     1             ( nriembegin, nriemend, iprimitive,
     1               ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1               gamma, smallnum, rhofloor, prsfloor,

     1               gasconst,
     1               cp_specie, molewt_specie, heat_of_formation,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

*******************************************************************

* Written by : Dinshaw S. Balsara

* Returns normalized left and right eigenvectors for a specified set
* of mean states for adiabatic Euler flow.

* When "iprimitive == 1" the eigenvectors expressed in terms of primitive
* variables are returned.

* When "iprimitive == 0" the eigenvectors in terms of conserved variables
* are returned.

* When "iprimitive == -1" the eigenvectors in terms of "conserved" variables
* are returned. Here we use the equation for entropy advection instead of
* the equation for total energy. This formulation is useful for maintaining
* pressure positivity.

* When primitive variables are used the vector of variables is in the
* order given by ( rho, vlx,y,z, prs, species_fractions)^T.

* When conserved variables are used the vector of variables is in the
* order given by ( rho, momx,y,z, eng, species_concentration)^T.

* When entropy-based "conserved" variables are used the vector of variables
* is in the order given by :
* ( rho, momx,y,z, prs/rho**(gamma-1), species_concentration)^T.

* Properly normalized eigenvectors as given in Roe and Balsara (1994)
* are generated.

*******************************************************************

      IMPLICIT NONE

************************************************************************
************************************************************************

      INTEGER, INTENT ( IN) :: nriembegin, nriemend, iprimitive,
     1      ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID

      REAL, INTENT ( IN) :: gamma, smallnum, prsfloor, rhofloor,
     1      gasconst

      REAL, DIMENSION ( 0 : NFLUID), INTENT ( IN) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1     INTENT ( INOUT) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1     INTENT ( INOUT) ::
     1      speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID), INTENT ( INOUT) ::
     1     eigenval

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID),
     1      INTENT ( INOUT) ::
     1     lf_eigenvec, rt_eigenvec

*****************************************************************
*****************************************************************

      INTEGER i, ii, iii, j, k, ifluid

* Variables for making wave speeds and eigenvectors.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     c_s, eng, temperature,
     1     mean_mol_wt, gamma_eos, gamma_soundspeed,
     1     dprs_drho, dprs_deng, deng_drho, deng_dprs,

     1     temparr_01, temparr_02, temparr_03

      REAL, DIMENSION ( 0: NFLUID,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     dprs_dspecie, deng_dspecie

* Temporary usage of left and right eigenvectors in the primitive variables.
* Transformation matrices to go from conserved to primitive and back.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID,
     1      NFIELD + NFLUID) ::

     1     lf_eigenvec1, rt_eigenvec1, dudw, dwdu

******************************************************************

!      COMMON /EIGENVEC_EULER_COM/

!     1     c_s, eng, temperature,
!     1     mean_mol_wt, gamma_eos, gamma_soundspeed,
!     1     dprs_drho, dprs_deng, deng_drho, deng_dprs,

!     1     dprs_dspecie, deng_dspecie,

!     1     temparr_01, temparr_02, temparr_03,

!     1     lf_eigenvec1, rt_eigenvec1, dudw, dwdu

#if ( MAKEOPENMP == 1)
!   !$OMP THREADPRIVATE (/EIGENVEC_EULER_COM/)
#endif

******************************************************************

      IF ( NFIELD .NE. 5) THEN
      WRITE ( 6, *)"EIGENVEC_EULER: NFIELD .NE. 5",NFIELD,5
      STOP
      END IF

      IF ( ( nriembegin .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriembegin .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EIGENVEC_EULER : nriembegin out of bounds",
     1       nriembegin, nriembegin
      STOP
      END IF

      IF ( ( nriemend .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriemend .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EIGENVEC_EULER : nriemend out of bounds",
     1       nriemend, nriemend
      STOP
      END IF


*****************************************************************

* Bound the density and pressure variables.

      DO i = nriembegin, nriemend

      rhobar ( i) = AMAX1 ( rhobar ( i), rhofloor)
      prsbar ( i) = AMAX1 ( prsbar ( i), prsfloor)

      END DO

******************************************************************

* Fix up "gamma_soundspeed" and partial derivatives of various
* thermodynamic quantities wrt one another.

* Because this is a demonstrator code, we just provide gamma-law gas formulae.

      DO i = nriembegin, nriemend

      eng ( i) = prsbar ( i) / (gamma - 1.0)
      mean_mol_wt ( i) = molewt_specie ( 0)
      gamma_eos ( i) = gamma
      gamma_soundspeed ( i) = gamma

      temperature ( i) = prsbar ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rhobar ( i))

      dprs_drho ( i) = 0.0
      dprs_deng ( i) = gamma - 1.0
      deng_drho ( i) = 0.0
      deng_dprs ( i) = 1.0 / ( gamma - 1.0)

      END DO

******************************************************************

* Make eigenvalues:

      DO i = nriembegin, nriemend

      c_s ( i) = gamma_soundspeed ( i) * prsbar ( i) / rhobar ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval ( i, 1) = vlxbar ( i) - c_s ( i)
      eigenval ( i, 2) = vlxbar ( i)
      eigenval ( i, 3) = vlxbar ( i)
      eigenval ( i, 4) = vlxbar ( i)
      eigenval ( i, 5) = vlxbar ( i) + c_s ( i)

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval ( i, NFIELD + 1: NFIELD + NFLUID) =
     1         vlxbar ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

*******************************************************************

* First zero out the eigenvectors. This is only needed because we can
* then fill in the species eigenvectors by putting 1's at suitable
* places on the diagonals.

      rt_eigenvec1 ( :, :, :) = 0.0
      lf_eigenvec1 ( :, :, :) = 0.0

*******************

* Make right & left eigenvectors for entropy waves and transverse
* velocity waves.

      DO i = nriembegin, nriemend

*******************

* Entropy wave

      rt_eigenvec1 ( i, 1, 2) = 1.0
      rt_eigenvec1 ( i, 2, 2) = 0.0
      rt_eigenvec1 ( i, 3, 2) = 0.0
      rt_eigenvec1 ( i, 4, 2) = 0.0
      rt_eigenvec1 ( i, 5, 2) = 0.0

      lf_eigenvec1 ( i, 2, 1) = 1.0
      lf_eigenvec1 ( i, 2, 2) = 0.0
      lf_eigenvec1 ( i, 2, 3) = 0.0
      lf_eigenvec1 ( i, 2, 4) = 0.0
      lf_eigenvec1 ( i, 2, 5) = - 1.0 / c_s ( i)**2

*******************

* Transverse y-velocity wave

      rt_eigenvec1 ( i, 1, 3) = 0.0
      rt_eigenvec1 ( i, 2, 3) = 0.0
      rt_eigenvec1 ( i, 3, 3) = 1.0
      rt_eigenvec1 ( i, 4, 3) = 0.0
      rt_eigenvec1 ( i, 5, 3) = 0.0

      lf_eigenvec1 ( i, 3, 1) = 0.0
      lf_eigenvec1 ( i, 3, 2) = 0.0
      lf_eigenvec1 ( i, 3, 3) = 1.0
      lf_eigenvec1 ( i, 3, 4) = 0.0
      lf_eigenvec1 ( i, 3, 5) = 0.0

*******************

* Transverse z-velocity wave

      rt_eigenvec1 ( i, 1, 4) = 0.0
      rt_eigenvec1 ( i, 2, 4) = 0.0
      rt_eigenvec1 ( i, 3, 4) = 0.0
      rt_eigenvec1 ( i, 4, 4) = 1.0
      rt_eigenvec1 ( i, 5, 4) = 0.0

      lf_eigenvec1 ( i, 4, 1) = 0.0
      lf_eigenvec1 ( i, 4, 2) = 0.0
      lf_eigenvec1 ( i, 4, 3) = 0.0
      lf_eigenvec1 ( i, 4, 4) = 1.0
      lf_eigenvec1 ( i, 4, 5) = 0.0

*******************

      END DO

*******************************************************************

* Make right & left eigenvectors for sound waves.

      DO i = nriembegin, nriemend

*******************

* Sound waves

*******************

      rt_eigenvec1 ( i, 1, 1) = rhobar ( i)
      rt_eigenvec1 ( i, 1, 5) = rhobar ( i)

      rt_eigenvec1 ( i, 2, 1) = - c_s ( i)
      rt_eigenvec1 ( i, 2, 5) = c_s ( i)

      rt_eigenvec1 ( i, 3, 1) = 0.0
      rt_eigenvec1 ( i, 3, 5) = 0.0

      rt_eigenvec1 ( i, 4, 1) = 0.0
      rt_eigenvec1 ( i, 4, 5) = 0.0

      rt_eigenvec1 ( i, 5, 1) = rhobar ( i) * c_s ( i)**2
      rt_eigenvec1 ( i, 5, 5) = rt_eigenvec1 ( i, 5, 1)

*******************

      temparr_01 ( i) = 0.5 / c_s ( i)**2

      lf_eigenvec1 ( i, 1, 1) = 0.0
      lf_eigenvec1 ( i, 5, 1) = 0.0

      lf_eigenvec1 ( i, 1, 2) = - temparr_01 ( i) * c_s ( i)
      lf_eigenvec1 ( i, 5, 2) = - lf_eigenvec1 ( i, 1, 2)

      lf_eigenvec1 ( i, 1, 3) = 0.0
      lf_eigenvec1 ( i, 5, 3) = 0.0

      lf_eigenvec1 ( i, 1, 4) = 0.0
      lf_eigenvec1 ( i, 5, 4) = 0.0

      lf_eigenvec1 ( i, 1, 5) = temparr_01 ( i) / rhobar ( i)
      lf_eigenvec1 ( i, 5, 5) = lf_eigenvec1 ( i, 1, 5)

*******************

      END DO

*******************

* For the multiple species, fix up the eigenvectors that carry the
* contributions from the species fractions. The zero components of
* the eigenvectors have already been previously assigned.

      IF ( NFLUID .GE. 1) THEN

      DO j = NFIELD + 1, NFIELD + NFLUID

      DO i = nriembegin, nriemend

      rt_eigenvec1 ( i, j, j) = 1.0
      lf_eigenvec1 ( i, j, j) = 1.0

      END DO

      END DO

      END IF
! End of "NFLUID .GE. 1".

*******************************************************************

* At this point the eigenvectors in the primitive variables are made
* and if only that much is desired, we can return at this point.

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

******************************************************************

      IF ( iprimitive .EQ. 0) THEN

******************************************************************

* Make transformation matrices to the conserved variables. First for
* the right eigenvectors. Don't make those matrix elements that are zero.

      DO i = nriembegin, nriemend


      dudw ( i, 1, 1) = 1.0
*      dudw ( i, 1, 2) = 0.0
*      dudw ( i, 1, 3) = 0.0
*      dudw ( i, 1, 4) = 0.0
*      dudw ( i, 1, 5) = 0.0


      dudw ( i, 2, 1) = vlxbar ( i)
      dudw ( i, 2, 2) = rhobar ( i)
*      dudw ( i, 2, 3) = 0.0
*      dudw ( i, 2, 4) = 0.0
*      dudw ( i, 2, 5) = 0.0


      dudw ( i, 3, 1) = vlybar ( i)
*      dudw ( i, 3, 2) = 0.0
      dudw ( i, 3, 3) = rhobar ( i)
*      dudw ( i, 3, 4) = 0.0
*      dudw ( i, 3, 5) = 0.0


      dudw ( i, 4, 1) = vlzbar ( i)
*      dudw ( i, 4, 2) = 0.0
*      dudw ( i, 4, 3) = 0.0
      dudw ( i, 4, 4) = rhobar ( i)
*      dudw ( i, 4, 5) = 0.0


      dudw ( i, 5, 1) = 0.5 * ( vlxbar ( i)**2
     1                        + vlybar ( i)**2 + vlzbar ( i)**2)
     1                + deng_drho ( i)
      dudw ( i, 5, 2) = rhobar ( i) * vlxbar ( i)
      dudw ( i, 5, 3) = rhobar ( i) * vlybar ( i)
      dudw ( i, 5, 4) = rhobar ( i) * vlzbar ( i)
      dudw ( i, 5, 5) = deng_dprs ( i)

      END DO

******************************************************************

* Make transformation matrices to the conserved variables. Then for
* the left eigenvectors. Don't make those matrix elements that are zero.

      DO i = nriembegin, nriemend

      temparr_01 ( i) = 1.0 / rhobar ( i)


      dwdu ( i, 1, 1) = 1.0
      dwdu ( i, 2, 1) = - temparr_01 ( i) * vlxbar ( i)
      dwdu ( i, 3, 1) = - temparr_01 ( i) * vlybar ( i)
      dwdu ( i, 4, 1) = - temparr_01 ( i) * vlzbar ( i)
      dwdu ( i, 5, 1) = dprs_drho ( i) + 0.5 * dprs_deng ( i)
     1                * ( vlxbar ( i)**2
     1                  + vlybar ( i)**2 + vlzbar ( i)**2)


*      dwdu ( i, 1, 2) = 0.0
      dwdu ( i, 2, 2) = temparr_01 ( i)
*      dwdu ( i, 3, 2) = 0.0
*      dwdu ( i, 4, 2) = 0.0
      dwdu ( i, 5, 2) = - dprs_deng ( i) * vlxbar ( i)


*      dwdu ( i, 1, 3) = 0.0
*      dwdu ( i, 2, 3) = 0.0
      dwdu ( i, 3, 3) = temparr_01 ( i)
*      dwdu ( i, 4, 3) = 0.0
      dwdu ( i, 5, 3) = - dprs_deng ( i) * vlybar ( i)


*      dwdu ( i, 1, 4) = 0.0
*      dwdu ( i, 2, 4) = 0.0
*      dwdu ( i, 3, 4) = 0.0
      dwdu ( i, 4, 4) = temparr_01 ( i)
      dwdu ( i, 5, 4) = - dprs_deng ( i) * vlzbar ( i)


*      dwdu ( i, 1, 5) = 0.0
*      dwdu ( i, 2, 5) = 0.0
*      dwdu ( i, 3, 5) = 0.0
*      dwdu ( i, 4, 5) = 0.0
      dwdu ( i, 5, 5) = dprs_deng ( i)

      END DO

******************************************************************

* Do the matrix multiplications. Then save the eigenvalues. This
* gives eigenvectors in the conserved variables.
* Don't multiply those matrix elements that are zero.


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


      END DO

      END DO

* At this point the eigenvectors in the conserved variables are made
* for a single species gas.

*******************

* For the situation where multiple species are present, fix up the
* eigenvectors. This consists of putting in the additional terms in
* the already existing five Euler flow eigenvectors. Furthermore, it also
* consists of constructing the remaining eigenvectors that carry the
* species contributions.

      IF ( NFLUID .GE. 1) THEN

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

* Extra parts of the first five right eigenvectors for Euler flow.
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

* Do the entire right eigenvectors for each of the species.
      rt_eigenvec ( i, 1:4, k) = 0.0
      rt_eigenvec ( i, 5, k) = deng_dspecie ( ifluid, i)
      rt_eigenvec ( i, 6:k-1, k) = 0.0
      rt_eigenvec ( i, k, k) = rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      rt_eigenvec ( i, k+1:NFIELD + NFLUID, k) = 0.0
      END IF

* Extra parts of the first five left eigenvectors for Euler flow.
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

* Do the entire left eigenvectors for each of the species.
      lf_eigenvec ( i, k, 1) = - speciebar ( ifluid, i)
     1                       / rhobar ( i)
      lf_eigenvec ( i, k, 2:k-1) = 0.0
      lf_eigenvec ( i, k, k) = 1.0 / rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      lf_eigenvec ( i, k, k+1:NFIELD + NFLUID) = 0.0
      END IF

      END DO

      END DO

      END IF
! End of "NFLUID .GE. 1".

******************************************************************

      END IF

******************************************************************

      IF ( iprimitive .EQ. -1) THEN

******************************************************************

* Make transformation matrices to the fake "conserved" variables. These
* are the variables for maintaining pressure positivity. First for
* the right eigenvectors. Don't make those matrix elements that are zero.

* Even when multiple species are used, these eigenvectors are based
* on a "frozen local gamma approximation". Thus, the ratio of specific
* heats, gamma, can vary from zone to zone but is assumed constant within
* each zone. This is justified because gamma is a very slowly varying
* function of the species fractions so it is ok to think of it as being
* (almost) exactly a constant.

      DO i = nriembegin, nriemend

      temparr_01 ( i) = 1.0 / rhobar ( i)
      temparr_02 ( i) = gamma_soundspeed ( i) - 1.0
      temparr_03 ( i) = rhobar ( i)**temparr_02 ( i)

      dudw ( i, 1, 1) = 1.0
*      dudw ( i, 1, 2) = 0.0
*      dudw ( i, 1, 3) = 0.0
*      dudw ( i, 1, 4) = 0.0
*      dudw ( i, 1, 5) = 0.0


      dudw ( i, 2, 1) = vlxbar ( i)
      dudw ( i, 2, 2) = rhobar ( i)
*      dudw ( i, 2, 3) = 0.0
*      dudw ( i, 2, 4) = 0.0
*      dudw ( i, 2, 5) = 0.0


      dudw ( i, 3, 1) = vlybar ( i)
*      dudw ( i, 3, 2) = 0.0
      dudw ( i, 3, 3) = rhobar ( i)
*      dudw ( i, 3, 4) = 0.0
*      dudw ( i, 3, 5) = 0.0


      dudw ( i, 4, 1) = vlzbar ( i)
*      dudw ( i, 4, 2) = 0.0
*      dudw ( i, 4, 3) = 0.0
      dudw ( i, 4, 4) = rhobar ( i)
*      dudw ( i, 4, 5) = 0.0


      dudw ( i, 5, 1) = - temparr_02 ( i) * prsbar ( i)
     1                * temparr_01 ( i) / temparr_03 ( i)
      dudw ( i, 5, 2) = 0.0
      dudw ( i, 5, 3) = 0.0
      dudw ( i, 5, 4) = 0.0
      dudw ( i, 5, 5) = 1.0 / temparr_03 ( i)


      END DO

******************************************************************

* Make transformation matrices to the fake "conserved" variables. These
* are the variables for maintaining pressure positivity. Then for
* the left eigenvectors. Don't make those matrix elements that are zero.

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


*      dwdu ( i, 1, 2) = 0.0
      dwdu ( i, 2, 2) = temparr_01 ( i)
*      dwdu ( i, 3, 2) = 0.0
*      dwdu ( i, 4, 2) = 0.0
      dwdu ( i, 5, 2) = 0.0


*      dwdu ( i, 1, 3) = 0.0
*      dwdu ( i, 2, 3) = 0.0
      dwdu ( i, 3, 3) = temparr_01 ( i)
*      dwdu ( i, 4, 3) = 0.0
      dwdu ( i, 5, 3) = 0.0


*      dwdu ( i, 1, 4) = 0.0
*      dwdu ( i, 2, 4) = 0.0
*      dwdu ( i, 3, 4) = 0.0
      dwdu ( i, 4, 4) = temparr_01 ( i)
      dwdu ( i, 5, 4) = 0.0


*      dwdu ( i, 1, 5) = 0.0
*      dwdu ( i, 2, 5) = 0.0
*      dwdu ( i, 3, 5) = 0.0
*      dwdu ( i, 4, 5) = 0.0
      dwdu ( i, 5, 5) = temparr_03 ( i)


      END DO

******************************************************************

* Do the matrix multiplications. Then save the eigenvalues. This
* gives eigenvectors in the fake "conserved" variables. These
* are needed to maintain the pressure positivity.
* Don't multiply those matrix elements that are zero.


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


      END DO

      END DO

* At this point the eigenvectors in the fake "conserved" variables are made
* for a single species gas.

*******************

* For the situation where multiple species are present, fix up the
* eigenvectors. This consists of putting in the additional terms in
* the already existing five Euler fluid eigenvectors. Furthermore, it also
* consists of constructing the remaining eigenvectors that carry the
* species contributions.

      IF ( NFLUID .GE. 1) THEN

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

* Extra parts of the first five right eigenvectors for Euler flow.
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

* Do the entire right eigenvectors for each of the species.
      rt_eigenvec ( i, 1:k-1, k) = 0.0
      rt_eigenvec ( i, k, k) = rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      rt_eigenvec ( i, k+1:NFIELD + NFLUID, k) = 0.0
      END IF

* Do the entire left eigenvectors for each of the species.
      lf_eigenvec ( i, k, 1) = - speciebar ( ifluid, i)
     1                       / rhobar ( i)
      lf_eigenvec ( i, k, 2:k-1) = 0.0
      lf_eigenvec ( i, k, k) = 1.0 / rhobar ( i)
      IF ( k + 1 .LE. NFIELD + NFLUID) THEN
      lf_eigenvec ( i, k, k+1:NFIELD + NFLUID) = 0.0
      END IF

      END DO

      END DO

      END IF
! End of "NFLUID .GE. 1".

******************************************************************

      END IF

******************************************************************
******************************************************************

      RETURN

      END SUBROUTINE EIGENVEC_EULER

******************************************************************
*******************************************************************

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! We have the option of minimizing the number of EOS calls
! (which can be quite expensive) by setting "use_approx_eos = 1".
! This can be a risky step if "gamma_eos" is rapidly varying.

      SUBROUTINE EULER_RIEM_ROE

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

* Written by : Dinshaw S. Balsara

* Does the riemann solver based on Balsara's modification of the Roe-type
* riemann solver. See Balsara, Astrophysical J. Supp., vol. 132, pg. 1, (2001)
* and Roe, J. Comput. Phys., vol. 43, pg. 357, (1981) and
* Harten & Hymnan, J. Comput. Phys., vol. 50, pg. 297, (1983).

* The x-momentum flux_re is split into an advected part and a ( thermal +
* magnetic) pressure part.

* ientropyfix = 1 will give the maximal dissipation by using the LLF flux
* even for characteristic fields that go into eachother.

* ientropyfix = 2 will provide an entropy fix when characteristic fields are
* going apart but not when they are going into eachother.

* ientropyfix = 3 will provide the basic Harten and Hyman entropy fix.

* In those zone boundaries where the pressure is negative, an Einfeld
* or LLF flux is used. A few lines at the beginning of the last
* loop turn this into an LLF flux. Since that is more stable, we keep
* those lines uncommented so that the default flux is the LLF flux.

* While the Roe type riemann solver is little more accurate than the
* HLL solver, it can become very expensive when several species
* are involved.

*******************************************************************

      IMPLICIT NONE

      EXTERNAL EIGENVECMHD, RHO_PRS_IN_ENG_T_GAMA_OUT_EULER

*******************************************************************

      INTEGER :: nriembegin, nriemend, 
     1        ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,
     1        ientropyfix, iprimitive, use_approx_eos,
     1        split_prsflux_re

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor, 
     1     fracdiff, gasconst, divratio_einfeldt_rs

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr

      REAL, DIMENSION
     1   ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD + NFLUID, NFIELD + NFLUID) ::
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
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION
     1    ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

************************************************************************
************************************************************************

* Past this point we only declare variables that are internal to the
* subroutine.

      INTEGER i, ii, j, k, ifield, ifluid, num_ein

      REAL small_speed,
     1     rhofloor_plus_delta, prsfloor_plus_delta

* Variables used for enoized riemann solver.

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      iwave_1, iwave_2, make_2, make_entropyfix

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                     NFIELD + NFLUID) ::

     1        use_eigenvec, do_roe


      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID) ::

     1      zonebdy_eigenval, eigenval_lf, eigenval_rt, del_cons


      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD + NFLUID) ::

     1      eigenval_db


      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID) ::

     1      eigenwt_cons, eigenwt_flux,
     1      eigenval_1, eigenval_2, conserved_2

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     momxremi, momxrepl, momyremi, momyrepl,
     1     momzremi, momzrepl, engl, engr, eps_zone,

     1     proj_flux_plus, proj_flux_minus,
     1     proj_cons_plus, proj_cons_minus,

     1     coef_flux_plus_f, coef_flux_plus_c,
     1     coef_flux_minus_f, coef_flux_minus_c,

     1     coef_cons_plus, coef_cons_minus, wavspd_1, wavspd_2,

     1     temparr_01, temparr_02, temparr_03, ratio_ein


******************************************************************
******************************************************************

* Variables used for Einfeldt fix.

* "conserved_ein, flux_ein" contain the resolved fields and fluxes
* on either side of the riemann problem.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID) ::

     1     conserved_ein_l, conserved_ein_r,
     1     flux_ein_l, flux_ein_r, conserved_ein, flux_ein


* These should be loaded with the min and max speeds from the roe solver.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     wavespd_ein_l, wavespd_ein_r

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     eng, temperature,
     1     mean_mol_wt, gamma_eos, gamma_soundspeed,
     1     engremi, temperatureremi, total_enthalpy_remi,
     1     engrepl, temperaturerepl, total_enthalpy_repl,
     1     mean_mol_wt_re, gamma_eos_re, gamma_soundspeed_re,

     1     entropy_ein_l, entropy_flux_ein_l,
     1     entropy_ein_r, entropy_flux_ein_r

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1        indx_ein

******************************************************************
******************************************************************

* Variables used for EIGENVEC_EULER


      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     c_s,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, total_enthalpy_bar,

     1     dprs_drho, dprs_deng, deng_drho, deng_dprs


      REAL, DIMENSION
     1    ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      speciebar, dprs_dspecie, deng_dspecie

! ----------------------------------------------------------------------

* Make some settings.
* Because floor values of density and pressure have already been taken,
* we set a bottom value that is a little above that floor.


      rhofloor_plus_delta = 1.000001 * rhofloor
      prsfloor_plus_delta = 1.000001 * prsfloor


      IF ( NFIELD .NE. 5) THEN
      WRITE ( 6, *)"EULER_RIEM_ROE : NFIELD .NE. 5", NFIELD,5
      debug_stop = 1
      RETURN
      END IF

      IF ( ( nriembegin .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriembegin .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EULER_RIEM_ROE : nriembegin out of bounds",
     1       nriembegin, ONEDFLUIDBEGIN, ONEDFLUIDEND
      debug_stop = 1
      RETURN
      END IF

      IF ( ( nriemend .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriemend .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EULER_RIEM_ROE : nriemend out of bounds",
     1       nriemend, ONEDFLUIDBEGIN, ONEDFLUIDEND
      debug_stop = 1
      RETURN
      END IF

      IF ( ( iprimitive .LT. -1) .OR.
     1     ( iprimitive .GT. 0) ) THEN
      WRITE ( 6, *)"EULER_RIEM_ROE : iprimitive out of bounds",
     1  iprimitive
      debug_stop = 1
      RETURN
      END IF

      IF ( ( ientropyfix .LT. 1) .OR.
     1     ( ientropyfix .GT. 3) ) THEN
      WRITE ( 6, *)"EULER_RIEM_ROE : ientropyfix out of bounds",
     1  ientropyfix
      debug_stop = 1
      RETURN
      END IF

      IF ( ( split_prsflux_re .LT. 0) .OR.
     1     ( split_prsflux_re .GT. 1)) THEN
      WRITE ( 6, *)"EULER_RIEM_ROE : split_prsflux_re ",
     1             " out of bounds", split_prsflux_re
      debug_stop = 1
      RETURN
      END IF


******************************************************************
******************************************************************

* Bound the density and pressure variables.
* For starters none of the zones need Einfeldt fix.
* Also "rhobar 2" and "prsbar 2" may not be evaluated in all cases
* later on so we give them a dummy initialization.

      DO i = nriembegin, nriemend

      do_ein ( i) = 0

      rhoremi ( i) = AMAX1 ( rhoremi ( i), rhofloor)
      rhorepl ( i) = AMAX1 ( rhorepl ( i), rhofloor)

      prsremi ( i) = AMAX1 ( prsremi ( i), prsfloor)
      prsrepl ( i) = AMAX1 ( prsrepl ( i), prsfloor)
 
      END DO

******************************************************************
******************************************************************

* Evaluate the thermal energy density and the local gammas. In doing so, we
* also make sure the species fractions are individually and cumulatively
* within reasonable ranges. Do this only if exact EOS's are to be used.

      IF ( use_approx_eos .NE. 1) THEN

      WRITE ( 6, *)
     1      "Stopping because EOS routines not present in demo. code."
      STOP

      ELSE

* Relying on approximate EOS's here.

      DO i = nriembegin, nriemend

      engremi ( i) = prsremi ( i) / ( gamma_eos_remi ( i) - 1.0)
      temperatureremi ( i) = prsremi ( i) * mean_mol_wt_remi ( i)
     1                     / ( gasconst * rhoremi ( i))

      engrepl ( i) = prsrepl ( i) / ( gamma_eos_repl ( i) - 1.0)
      temperaturerepl ( i) = prsrepl ( i) * mean_mol_wt_repl ( i)
     1                     / ( gasconst * rhorepl ( i))

      END DO

      END IF

***************************

* Make "total_enthalpy_remi" and "total_enthalpy_repl". These are
* useful in making Roe-averages.

      DO i = nriembegin, nriemend

      total_enthalpy_remi ( i) = ( engremi ( i) + prsremi ( i)
     1    + 0.5 * rhoremi ( i) * ( vlxremi ( i)**2 + vlyremi ( i)**2
     1                           + vlzremi ( i)**2) ) / rhoremi ( i)

      total_enthalpy_repl ( i) = ( engrepl ( i) + prsrepl ( i)
     1    + 0.5 * rhorepl ( i) * ( vlxrepl ( i)**2 + vlyrepl ( i)**2
     1                           + vlzrepl ( i)**2) ) / rhorepl ( i)

      END DO

******************************************************************
******************************************************************

* Make "zonebdy_conserved" and "zonebdy_flux". Uncommenting this part would
* make this a regular linearized Riemann solver. We only need to comment
* it out when using ADER schemes because they may produce higher order fluxes
* and conserved variables for use in the Riemann solver.

      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 1) = rhoremi ( i)
      zonebdy_conserved ( i, 0, 2) = rhoremi ( i)
     1                             * vlxremi ( i)
      zonebdy_conserved ( i, 0, 3) = rhoremi ( i)
     1                             * vlyremi ( i)
      zonebdy_conserved ( i, 0, 4) = rhoremi ( i)
     1                             * vlzremi ( i)


      zonebdy_conserved ( i, 1, 1) = rhorepl ( i)
      zonebdy_conserved ( i, 1, 2) = rhorepl ( i)
     1                             * vlxrepl ( i)
      zonebdy_conserved ( i, 1, 3) = rhorepl ( i)
     1                             * vlyrepl ( i)
      zonebdy_conserved ( i, 1, 4) = rhorepl ( i)
     1                             * vlzrepl ( i)


      zonebdy_flux ( i, 0, 1) = zonebdy_conserved ( i, 0, 2)
      zonebdy_flux ( i, 0, 2) =
     1     zonebdy_conserved ( i, 0, 2) * vlxremi ( i)
     1   + prsremi ( i)
      zonebdy_flux ( i, 0, 3) =
     1     zonebdy_conserved ( i, 0, 3) * vlxremi ( i)
      zonebdy_flux ( i, 0, 4) =
     1     zonebdy_conserved ( i, 0, 4) * vlxremi ( i)


      zonebdy_flux ( i, 1, 1) = zonebdy_conserved ( i, 1, 2)
      zonebdy_flux ( i, 1, 2) =
     1     zonebdy_conserved ( i, 1, 2) * vlxrepl ( i)
     1   + prsrepl ( i)
      zonebdy_flux ( i, 1, 3) =
     1     zonebdy_conserved ( i, 1, 3) * vlxrepl ( i)
      zonebdy_flux ( i, 1, 4) =
     1     zonebdy_conserved ( i, 1, 4) * vlxrepl ( i)


      END DO


* Make the fifth component differently depending on which flux
* is desired.

      IF ( iprimitive .EQ. 0) THEN

      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 5) =
     1      engremi ( i)
     1    + 0.5 * rhoremi ( i)
     1    * ( vlxremi ( i)**2 + vlyremi ( i)**2
     1                        + vlzremi ( i)**2)

      zonebdy_conserved ( i, 1, 5) =
     1      engrepl ( i)
     1    + 0.5 * rhorepl ( i)
     1    * ( vlxrepl ( i)**2 + vlyrepl ( i)**2
     1                        + vlzrepl ( i)**2)

      zonebdy_flux ( i, 0, 5) =
     1     ( zonebdy_conserved ( i, 0, 5) + prsremi ( i) )
     1   * vlxremi ( i)

      zonebdy_flux ( i, 1, 5) =
     1     ( zonebdy_conserved ( i, 1, 5) + prsrepl ( i) )
     1   * vlxrepl ( i)

      END DO

      END IF


      IF ( iprimitive .EQ. - 1) THEN

* EOS
      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 5) = prsremi ( i)
     1        / rhoremi ( i)**( gamma_soundspeed_remi ( i) - 1.0)

      zonebdy_conserved ( i, 1, 5) = prsrepl ( i)
     1        / rhorepl ( i)**( gamma_soundspeed_repl ( i) - 1.0)

      zonebdy_flux ( i, 0, 5) = zonebdy_conserved ( i, 0, 5)
     1                        * vlxremi ( i)

      zonebdy_flux ( i, 1, 5) = zonebdy_conserved ( i, 1, 5)
     1                        * vlxrepl ( i)

      END DO

      END IF


* If multiple species are present then initialize their fluxes etc.
* appropriately also.

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        zonebdy_conserved ( i, 0, k) = rhoremi ( i)
     1                               * specieremi ( ifluid, i)

        zonebdy_conserved ( i, 1, k) = rhorepl ( i)
     1                               * specierepl ( ifluid, i)

        zonebdy_flux ( i, 0, k) = zonebdy_conserved ( i, 0, k)
     1                          * vlxremi ( i)

        zonebdy_flux ( i, 1, k) = zonebdy_conserved ( i, 1, k)
     1                          * vlxrepl ( i)

        conserved_re ( i, k) = 0.0

        flux_re ( i, k + 1) = 0.0

      END DO

      END DO

******************************************************************
******************************************************************

! Zero out "flux_re" here.

      DO i = nriembegin, nriemend

      conserved_re ( i, 1) = 0.0
      conserved_re ( i, 2) = 0.0
      conserved_re ( i, 3) = 0.0
      conserved_re ( i, 4) = 0.0
      conserved_re ( i, 5) = 0.0

      flux_re ( i, 1) = 0.0
      flux_re ( i, 2) = 0.0
      flux_re ( i, 3) = 0.0
      flux_re ( i, 4) = 0.0
      flux_re ( i, 5) = 0.0
      flux_re ( i, 6) = 0.0

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        conserved_re ( i, k) = 0.0

        flux_re ( i, k + 1) = 0.0

      END DO

      END DO


******************************************************************
******************************************************************

* Here we use Roe averaging and make only one call to the
* eigenvector subroutine. The "eng, temperature" etc. that we make here
* can even serve as good initial guesses for an exact EOS.


      DO i = nriembegin, nriemend

! ----------

      gamma_eos ( i) = 0.5 * ( gamma_eos_remi ( i)
     1                       + gamma_eos_repl ( i))
      gamma_soundspeed ( i) = 0.5 * ( gamma_soundspeed_remi ( i)
     1                              + gamma_soundspeed_repl ( i))
      mean_mol_wt ( i) = 0.5 * ( mean_mol_wt_remi ( i)
     1                         + mean_mol_wt_repl ( i))

! ----------

      temparr_01 ( i) = SQRT ( rhorepl ( i) / rhoremi ( i))
      temparr_02 ( i) = 1.0 / ( 1.0 + temparr_01 ( i))

      rhobar ( i) = SQRT ( rhoremi ( i) * rhorepl ( i))

      vlxbar ( i) = ( vlxremi ( i) + vlxrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)
      vlybar ( i) = ( vlyremi ( i) + vlyrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)
      vlzbar ( i) = ( vlzremi ( i) + vlzrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)

      total_enthalpy_bar ( i) = ( total_enthalpy_remi ( i)
     1  + total_enthalpy_repl ( i) * temparr_01 ( i)) * temparr_02 ( i)

      prsbar ( i) = total_enthalpy_bar ( i) - 0.5 * ( vlxbar ( i)**2
     1                     + vlybar ( i)**2 + vlzbar ( i)**2)
      prsbar ( i) = ( gamma_eos ( i) - 1.0) * rhobar ( i) * prsbar ( i)
     1            / gamma_eos ( i)
      prsbar ( i) = AMAX1 ( prsbar ( i), prsfloor)

! ----------

! Give yourself the option of using arithmetic averages if needed.

      IF ( ( rhobar ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prsbar ( i) .LE. prsfloor_plus_delta) ) THEN

      rhobar ( i) = 0.5 * ( rhoremi ( i) + rhorepl ( i))
      prsbar ( i) = 0.5 * ( prsremi ( i) + prsrepl ( i))

      vlxbar ( i) = 0.5 * ( vlxremi ( i) + vlxrepl ( i))
      vlybar ( i) = 0.5 * ( vlyremi ( i) + vlyrepl ( i))
      vlzbar ( i) = 0.5 * ( vlzremi ( i) + vlzrepl ( i))

      END IF

! ----------

      eng ( i) = prsbar ( i) / ( gamma_eos ( i) - 1.0)
      temperature ( i) = prsbar ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rhobar ( i))

! ----------

      END DO

! ------------------------------

      DO ifluid = 1, NFLUID

      DO i = nriembegin, nriemend

      temparr_01 ( i) = SQRT ( rhorepl ( i) / rhoremi ( i))
      temparr_02 ( i) = 1.0 / ( 1.0 + temparr_01 ( i))

      speciebar ( ifluid, i) = ( specieremi ( ifluid, i)
     1 + specierepl ( ifluid, i) * temparr_01 ( i)) * temparr_02 ( i)

      END DO

      END DO


* Using exact EOS's here. If that is not needed, then "eng, temperature"
* etc. have already been made above and the present call is not needed.

      IF ( use_approx_eos .NE. 1) THEN

      WRITE ( 6, *)
     1      "Stopping because EOS routines not present in demo. code."

      END IF

*********************************************************************
*********************************************************************

* Make eigenvectors and eigenvalues for the arithmetic averaged
* state of the riemann problem.

* The RHO_PRS_IN_ENG_T_GAMA_OUT_EULER subroutine that has been called before
* will automatically check to make sure that "speciebar" is within the
* limits.

      CALL EIGENVEC_EULER

     1           ( nriembegin, nriemend, iprimitive,
     1             ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1             gamma, smallnum, rhofloor, prsfloor,

     1             gasconst,
     1             cp_specie, molewt_specie, heat_of_formation,

     1             rhobar, prsbar,
     1             vlxbar, vlybar, vlzbar,
     1             speciebar,

     1             eigenval_re,
     1             lf_eigenvec_re, rt_eigenvec_re )

*********************************************************************
*********************************************************************

      IF ( debug_stop .EQ. 1) THEN
      RETURN
      END IF


* And store the eigenvalues on either side of
* the zone boundary.

      DO i = nriembegin, nriemend

      zonebdy_eigenval ( i, 1) = eigenval_re ( i, 1)
      eigenval_db ( i, 0, 1) = eigenval_re ( i, 1)
      eigenval_db ( i, 1, 1) = eigenval_re ( i, 1)

      zonebdy_eigenval ( i, 2) = eigenval_re ( i, 2)
      eigenval_db ( i, 0, 2) = eigenval_re ( i, 2)
      eigenval_db ( i, 1, 2) = eigenval_re ( i, 2)

      zonebdy_eigenval ( i, 3) = eigenval_re ( i, 3)
      eigenval_db ( i, 0, 3) = eigenval_re ( i, 3)
      eigenval_db ( i, 1, 3) = eigenval_re ( i, 3)

      zonebdy_eigenval ( i, 4) = eigenval_re ( i, 4)
      eigenval_db ( i, 0, 4) = eigenval_re ( i, 4)
      eigenval_db ( i, 1, 4) = eigenval_re ( i, 4)

      zonebdy_eigenval ( i, 5) = eigenval_re ( i, 5)
      eigenval_db ( i, 0, 5) = eigenval_re ( i, 5)
      eigenval_db ( i, 1, 5) = eigenval_re ( i, 5)

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        zonebdy_eigenval ( i, k) = eigenval_re ( i, k)
        eigenval_db ( i, 0, k) = eigenval_re ( i, k)
        eigenval_db ( i, 1, k) = eigenval_re ( i, k)

      END DO

      END DO

*********************************************************************
*********************************************************************

* Here we explicitly need to get the eigenvalues on either
* side of the zone boundary so that we can maximize / minimize w.r.t.
* the eigenvalues of the mean state. We reset "eigenval_db" here
* so that "eigenval_db ( 0, ifield, i)" is the left-most eigenvalue
* in that characteristic field and "eigenval_db ( 1, ifield, i)" is
* the right-most.

* Make "eigenval_lf":

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed_remi ( i) * prsremi ( i)
     1         / rhoremi ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_lf ( i, 1) = vlxremi ( i) - c_s ( i)
      eigenval_lf ( i, 2) = vlxremi ( i)
      eigenval_lf ( i, 3) = vlxremi ( i)
      eigenval_lf ( i, 4) = vlxremi ( i)
      eigenval_lf ( i, 5) = vlxremi ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_lf ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxremi ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_lf"

*********************************************************************
*********************************************************************

      DO i = nriembegin, nriemend


      eigenval_db ( i, 0, 1) =
     1      AMIN1 ( eigenval_db ( i, 0, 1), eigenval_lf ( i, 1) )

      eigenval_db ( i, 1, 1) =
     1      AMAX1 ( eigenval_db ( i, 1, 1), eigenval_lf ( i, 1) )


      eigenval_db ( i, 0, 2) =
     1      AMIN1 ( eigenval_db ( i, 0, 2), eigenval_lf ( i, 2) )

      eigenval_db ( i, 1, 2) =
     1      AMAX1 ( eigenval_db ( i, 1, 2), eigenval_lf ( i, 2) )


      eigenval_db ( i, 0, 3) =
     1      AMIN1 ( eigenval_db ( i, 0, 3), eigenval_lf ( i, 3) )

      eigenval_db ( i, 1, 3) =
     1      AMAX1 ( eigenval_db ( i, 1, 3), eigenval_lf ( i, 3) )


      eigenval_db ( i, 0, 4) =
     1      AMIN1 ( eigenval_db ( i, 0, 4), eigenval_lf ( i, 4) )

      eigenval_db ( i, 1, 4) =
     1      AMAX1 ( eigenval_db ( i, 1, 4), eigenval_lf ( i, 4) )


      eigenval_db ( i, 0, 5) =
     1      AMIN1 ( eigenval_db ( i, 0, 5), eigenval_lf ( i, 5) )

      eigenval_db ( i, 1, 5) =
     1      AMAX1 ( eigenval_db ( i, 1, 5), eigenval_lf ( i, 5) )


      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        eigenval_db ( i, 0, k) =
     1        AMIN1 ( eigenval_db ( i, 0, k), eigenval_lf ( i, k) )

        eigenval_db ( i, 1, k) =
     1        AMAX1 ( eigenval_db ( i, 1, k), eigenval_lf ( i, k) )

      END DO

      END DO

*********************************************************************
*********************************************************************

* Make "eigenval_rt":

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed_repl ( i) * prsrepl ( i)
     1         / rhorepl ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_rt ( i, 1) = vlxrepl ( i) - c_s ( i)
      eigenval_rt ( i, 2) = vlxrepl ( i)
      eigenval_rt ( i, 3) = vlxrepl ( i)
      eigenval_rt ( i, 4) = vlxrepl ( i)
      eigenval_rt ( i, 5) = vlxrepl ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_rt ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxrepl ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_rt"

*********************************************************************
*********************************************************************

      DO i = nriembegin, nriemend


      eigenval_db ( i, 0, 1) =
     1      AMIN1 ( eigenval_db ( i, 0, 1), eigenval_rt ( i, 1) )

      eigenval_db ( i, 1, 1) =
     1      AMAX1 ( eigenval_db ( i, 1, 1), eigenval_rt ( i, 1) )


      eigenval_db ( i, 0, 2) =
     1      AMIN1 ( eigenval_db ( i, 0, 2), eigenval_rt ( i, 2) )

      eigenval_db ( i, 1, 2) =
     1      AMAX1 ( eigenval_db ( i, 1, 2), eigenval_rt ( i, 2) )


      eigenval_db ( i, 0, 3) =
     1      AMIN1 ( eigenval_db ( i, 0, 3), eigenval_rt ( i, 3) )

      eigenval_db ( i, 1, 3) =
     1      AMAX1 ( eigenval_db ( i, 1, 3), eigenval_rt ( i, 3) )


      eigenval_db ( i, 0, 4) =
     1      AMIN1 ( eigenval_db ( i, 0, 4), eigenval_rt ( i, 4) )

      eigenval_db ( i, 1, 4) =
     1      AMAX1 ( eigenval_db ( i, 1, 4), eigenval_rt ( i, 4) )


      eigenval_db ( i, 0, 5) =
     1      AMIN1 ( eigenval_db ( i, 0, 5), eigenval_rt ( i, 5) )

      eigenval_db ( i, 1, 5) =
     1      AMAX1 ( eigenval_db ( i, 1, 5), eigenval_rt ( i, 5) )

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        eigenval_db ( i, 0, k) =
     1        AMIN1 ( eigenval_db ( i, 0, k), eigenval_rt ( i, k) )

        eigenval_db ( i, 1, k) =
     1        AMAX1 ( eigenval_db ( i, 1, k), eigenval_rt ( i, k) )

      END DO

      END DO

* End of making the eigenvectors and appropriate wave speeds.

************************************************************************
************************************************************************

* This is a simpler and easier entropy fix strategy.

* No Einfeldt fix is attempted.
* Likewise, no attempt is made to split the momentum flux_re into pressure and
* advected parts.
* The emphasis here is on speed!

      IF ( ientropyfix .EQ. 3) THEN

******************************************************************
******************************************************************

* Make the differences in conserved_re.

      DO k = 1, NFIELD + NFLUID

        DO i = nriembegin, nriemend

        del_cons ( i, k) = zonebdy_conserved ( i, 1, k)
     1                   - zonebdy_conserved ( i, 0, k)

        END DO

      END DO

! ------------------------------

* Project the conserved_re differences onto the right eigenvectors.

      eigenwt_re ( nriembegin : nriemend, :) = 0.0

* Left-multiply the kth eigenvector.
      DO k = 1, NFIELD + NFLUID

* dot product loop.
        DO j = 1, NFIELD + NFLUID

          DO i = nriembegin, nriemend

          eigenwt_re ( i, k) = eigenwt_re ( i, k)
     1                   + lf_eigenvec_re ( i, k, j)
     1                            * del_cons ( i, j)

          END DO

        END DO

      END DO

! ------------------------------

* Calculate the mean fluxes.

* kth flux vector component
      DO k = 1, NFIELD + NFLUID

        DO i = nriembegin, nriemend

        flux_re ( i, k) = 0.5 * ( zonebdy_flux ( i, 0, k)
     1                          + zonebdy_flux ( i, 1, k) )

        conserved_re ( i, k) =
     1                    0.5 * ( zonebdy_conserved ( i, 0, k)
     1                          + zonebdy_conserved ( i, 1, k) )

        END DO

      END DO

! ------------------------------

* Make entropy fix for the wavespeeds.

* kth wavespeed being considered here.
      DO k = 1, NFIELD + NFLUID

        DO i = nriembegin, nriemend

! --------------------

        IF ( ( eigenval_db ( i, 0, k) .LT. - smallnum) .AND.
     1       ( eigenval_db ( i, 1, k) .GT. smallnum) ) THEN

! ----------

! Rarefaction.


        temparr_01 ( i) = 1.0 / ( eigenval_db ( i, 1, k)
     1                          - eigenval_db ( i, 0, k))

! From Eqs. (6.41) and (6.42) of book; this is "0.5*(alpha^minus - alpha^plus)".
        eigenwt_cons ( i, k) =
     1   0.5 * ( ( eigenval_db ( i, 1, k) + eigenval_db ( i, 0, k))
     1         - 2.0 * eigenval_re ( i, k) )
     1 * temparr_01 ( i) * eigenwt_re ( i, k)

! Eqs. (6.41) and (6.43) of book; this is "lambda^minus".
        temparr_02 ( i) = eigenval_db ( i, 0, k)
     1 * ( eigenval_db ( i, 1, k) - eigenval_re ( i, k))
     1 * temparr_01 ( i)

! Eqs. (6.41) and (6.43) of book; this is "lambda^plus".
        temparr_03 ( i) = eigenval_db ( i, 1, k)
     1 * ( eigenval_re ( i, k) - eigenval_db ( i, 0, k))
     1 * temparr_01 ( i)

! From Eqn. (6.43) of book; this is "0.5*( lambda^minus - lambda^plus)*alpha^m".
        eigenwt_flux ( i, k) =
     1   0.5 * ( temparr_02 ( i) - temparr_03 ( i))
     1 * eigenwt_re ( i, k)

! ----------

        ELSE IF ( eigenval_re ( i, k) .LT. 0.0) THEN

! ----------

! Left-going wave, no rarefaction.


        eigenwt_cons ( i, k) = 0.5 * eigenwt_re ( i, k)

        eigenwt_flux ( i, k) = 0.5 * eigenval_re ( i, k)
     1                       * eigenwt_re ( i, k)

! ----------

        ELSE

! ----------

! Right-going wave, no rarefaction.


        eigenwt_cons ( i, k) = - 0.5 * eigenwt_re ( i, k)

        eigenwt_flux ( i, k) = - 0.5 * eigenval_re ( i, k)
     1                       * eigenwt_re ( i, k)

! ----------

        END IF

! --------------------

        END DO

      END DO

! ------------------------------

* Add in the contribution from the entropy fix to the mean fluxes
* and conserved variables.

* kth term in the vector of conserved variables.
      DO k = 1, NFIELD + NFLUID

* jth wave's contribution being considered here.
        DO j = 1, NFIELD + NFLUID

          DO i = nriembegin, nriemend

          flux_re ( i, k) = flux_re ( i, k)
     1               + eigenwt_flux ( i, j)
     1             * rt_eigenvec_re ( i, k, j)

          conserved_re ( i, k) = conserved_re ( i, k)
     1                         + eigenwt_cons ( i, j)
     1                       * rt_eigenvec_re ( i, k, j) 

          END DO

        END DO

      END DO


******************************************************************
******************************************************************

      END IF

* End of ientropyfix = 3 case.

******************************************************************
******************************************************************

* This entropy fix strategy goes all the way from here to the end of
* the subroutine.

      IF ( ( ientropyfix .EQ. 1) .OR.
     1     ( ientropyfix .EQ. 2) ) THEN

******************************************************************
******************************************************************

* Make "use_eigenvec" and "do_roe". "do_roe = 1" tells us that this zone
* boundary may use a roe construction. "use_eigenvec" tells us
* which eigensystem to use ( left or right) when roe construction
* is used. If roe construction is not used we use entropy fix construction.
*
* Note too that situations where the left and right wavespeeds are
* both 0 are treated as shocks. This prevents them from generating
* divide by zero in the entropy fix ( were they to be treated as
* rarefactions).

******************************************************************

      IF ( ientropyfix .EQ. 1) THEN

* Beginning of region that is specific to "ientropyfix = 1".

******************************************************************

* kth wavespeeds being considered here.
      DO k = 1, NFIELD + NFLUID

        DO i = nriembegin, nriemend

        IF ( ( eigenval_db ( i, 0, k) .GE. 0.0) .AND.
     1       ( eigenval_db ( i, 1, k) .GE. 0.0) ) THEN

        do_roe ( i, k) = 1
        use_eigenvec ( i, k) = 0

        ELSE IF ( ( eigenval_db ( i, 0, k) .LT. 0.0) .AND.
     1            ( eigenval_db ( i, 1, k) .LT. 0.0) ) THEN

        do_roe ( i, k) = 1
        use_eigenvec ( i, k) = 1

        ELSE

        do_roe ( i, k) = 0
        use_eigenvec ( i, k) = 1

        END IF

        END DO

      END DO


* For the case where the characteristic fields are going into each other
* and both characteristic fields are going into the zone boundary it may
* help to increase the local dissipation. If so, uncomment the below
* loop.

* kth wavespeeds being considered here.
      DO k = 1, NFIELD + NFLUID

        DO i = nriembegin, nriemend

        IF ( ( do_roe ( i, k) .EQ. 0) .AND.
     1       ( eigenval_lf ( i, k) .GT. 0.0) .AND.
     1       ( eigenval_rt ( i, k) .LT. 0.0) ) THEN

        eigenval_db ( i, 1, k) =
     1         AMAX1 ( ABS ( eigenval_db ( i, 1, k)),
     1                 ABS ( eigenval_db ( i, 0, k)),
     1                 ABS ( eigenval_lf ( i, k)),
     1                 ABS ( eigenval_rt ( i, k)) )

        eigenval_db ( i, 0, k) = - eigenval_db ( i, 1, k)

        END IF

        END DO

      END DO

******************************************************************

      END IF

* End of region that is specific to "ientropyfix = 1".

******************************************************************

      IF ( ientropyfix .EQ. 2) THEN

* Beginning of region that is specific to "ientropyfix = 2".

******************************************************************

* kth wavespeeds being considered here.
      DO k = 1, NFIELD + NFLUID

        DO i = nriembegin, nriemend

        IF ( ( eigenval_lf ( i, k) .LT. 0.0) .AND.
     1       ( eigenval_rt ( i, k) .GT. 0.0) ) THEN

        do_roe ( i, k) = 0
        use_eigenvec ( i, k) = 1

        ELSE IF ( eigenval_re ( i, k) .GE. 0.0) THEN

        do_roe ( i, k) = 1
        use_eigenvec ( i, k) = 0

        ELSE

        do_roe ( i, k) = 1
        use_eigenvec ( i, k) = 1

        END IF

        END DO

      END DO

******************************************************************

      END IF

* End of region that is specific to "ientropyfix = 2".

******************************************************************

* Loop over the characteristic fields.

      DO ifield = 1, NFIELD + NFLUID

******************************************************************

* Do the Roe fluxes and conserved variables here.


* Zero out dot products.

      DO i = nriembegin, nriemend

      proj_flux_plus ( i) = 0.0
      proj_cons_plus ( i) = 0.0

      END DO


* Fill in the dot products.

      DO i = nriembegin, nriemend


      IF ( do_roe ( i, ifield) .EQ. 1 ) THEN

      k = use_eigenvec ( i, ifield)

        DO j = 1, NFIELD + NFLUID

        proj_flux_plus ( i) = proj_flux_plus ( i)
     1                      + lf_eigenvec_re ( i, ifield, j)
     1                             * zonebdy_flux ( i, k, j)

        proj_cons_plus ( i) = proj_cons_plus ( i)
     1                      + lf_eigenvec_re ( i, ifield, j)
     1                        * zonebdy_conserved ( i, k, j)

        END DO

      END IF


      END DO


* Add in the suitable eigenvectors to the fluxes and conserved variables.

      DO i = nriembegin, nriemend


      IF ( do_roe ( i, ifield) .EQ. 1 ) THEN

      flux_re ( i, 1 : NFIELD + NFLUID) =
     1            flux_re ( i, 1 : NFIELD + NFLUID)
     1   + proj_flux_plus ( i)
     1   * rt_eigenvec_re ( i, 1 : NFIELD + NFLUID, ifield)

      conserved_re ( i, 1 : NFIELD + NFLUID) =
     1       conserved_re ( i, 1 : NFIELD + NFLUID)
     1   + proj_cons_plus ( i)
     1   * rt_eigenvec_re ( i, 1 : NFIELD + NFLUID, ifield)

      END IF


      END DO
*
******************************************************************

* Do the entropy fix for the fluxes and conserved variables here.

* Zero out dot products.

      DO i = nriembegin, nriemend

      proj_flux_plus ( i) = 0.0
      proj_cons_plus ( i) = 0.0

      proj_flux_minus ( i) = 0.0
      proj_cons_minus ( i) = 0.0

      END DO


* Make sure left and right eigenvalues are at least a little separated.

      DO i = nriembegin, nriemend

      IF ( ( do_roe ( i, ifield) .EQ. 0 ) .AND.
     1     ( eigenval_db ( i, 1, ifield)
     1     - eigenval_db ( i, 0, ifield) .LE. smallnum )
     1   ) 

     1 eigenval_db ( i, 0, ifield) = eigenval_db ( i, 1, ifield)
     1                            - smallnum

      END DO


* Make coefficients for the flux and conserved variable projections here.

      DO i = nriembegin, nriemend

      IF ( do_roe ( i, ifield) .EQ. 0 ) THEN


      temparr_01 ( i) = 1.0 / ( eigenval_db ( i, 1, ifield)
     1                        - eigenval_db ( i, 0, ifield) )


      coef_flux_plus_f ( i) = eigenval_db ( i, 1, ifield)
     1                      * temparr_01 ( i)

      coef_flux_plus_c ( i) = - eigenval_db ( i, 0, ifield)
     1                      * coef_flux_plus_f ( i)

      coef_flux_minus_f ( i) = - eigenval_db ( i, 0, ifield)
     1                      * temparr_01 ( i)

      coef_flux_minus_c ( i) = - coef_flux_plus_c ( i)


      coef_cons_plus ( i) = ( zonebdy_eigenval ( i, ifield)
     1                      - eigenval_db ( i, 0, ifield) )
     1                    * temparr_01 ( i)

      coef_cons_minus ( i) = ( eigenval_db ( i, 1, ifield)
     1                       - zonebdy_eigenval ( i, ifield) )
     1                    * temparr_01 ( i)

      END IF

      END DO


* Fill in the dot products.

      DO i = nriembegin, nriemend

      IF ( do_roe ( i, ifield) .EQ. 0 ) THEN

        DO j = 1, NFIELD + NFLUID

        proj_flux_plus ( i) = proj_flux_plus ( i)
     1                      + lf_eigenvec_re ( i, ifield, j)
     1                        * zonebdy_flux ( i, 0, j)

        proj_cons_plus ( i) = proj_cons_plus ( i)
     1                      + lf_eigenvec_re ( i, ifield, j)
     1                   * zonebdy_conserved ( i, 0, j)

        proj_flux_minus ( i) = proj_flux_minus ( i)
     1                        + lf_eigenvec_re ( i, ifield, j)
     1                          * zonebdy_flux ( i, 1, j)

        proj_cons_minus ( i) = proj_cons_minus ( i)
     1                        + lf_eigenvec_re ( i, ifield, j)
     1                     * zonebdy_conserved ( i, 1, j)

        END DO

      END IF

      END DO


* Add in the suitable eigenvectors to the fluxes and conserved variables.

      DO i = nriembegin, nriemend


      IF ( do_roe ( i, ifield) .EQ. 0 ) THEN


      temparr_01 ( i) =

     1       coef_flux_plus_f ( i) * proj_flux_plus ( i)
     1     + coef_flux_plus_c ( i) * proj_cons_plus ( i)

     1     + coef_flux_minus_f ( i) * proj_flux_minus ( i)
     1     + coef_flux_minus_c ( i) * proj_cons_minus ( i)

      flux_re ( i, 1 : NFIELD + NFLUID) =
     1            flux_re ( i, 1 : NFIELD + NFLUID)
     1       + temparr_01 ( i)
     1   * rt_eigenvec_re ( i, 1 : NFIELD + NFLUID, ifield)


      temparr_02 ( i) =

     1       coef_cons_plus ( i) * proj_cons_plus ( i)

     1     + coef_cons_minus ( i) * proj_cons_minus ( i)

      conserved_re ( i, 1 : NFIELD + NFLUID) =
     1       conserved_re ( i, 1 : NFIELD + NFLUID)
     1       + temparr_02 ( i)
     1   * rt_eigenvec_re ( i, 1 : NFIELD + NFLUID, ifield)

      END IF


      END DO


******************************************************************

      END DO

* End of ifield loop.

******************************************************************
******************************************************************

* This way of obtaining the entropy flux ensures that even during a
* normal call to the Riemann solver we get an HLL-based entropy
* flux with which to enforce pressure positivity.

      DO i = nriembegin, nriemend

      entropy_ein_l ( i) = prsremi ( i)
     1        / rhoremi ( i)**( gamma_soundspeed_remi ( i) - 1.0)
      entropy_ein_r ( i) = prsrepl ( i)
     1        / rhorepl ( i)**( gamma_soundspeed_repl ( i) - 1.0)

      entropy_flux_ein_l ( i) = entropy_ein_l ( i) * vlxremi ( i)
      entropy_flux_ein_r ( i) = entropy_ein_r ( i) * vlxrepl ( i)


      wavespd_ein_l ( i) = AMIN1 ( eigenval_db ( i, 0, 1),
     1                             eigenval_db ( i, 0, 2) )

      wavespd_ein_r ( i) = AMAX1 ( eigenval_db ( i, 1, 4),
     1                             eigenval_db ( i, 1, 5) )

      wavespd_ein_l ( i) = AMIN1 ( wavespd_ein_l ( i), 0.0)
      wavespd_ein_r ( i) = AMAX1 ( wavespd_ein_r ( i), 0.0)

      temparr_01 ( i) = 1.0 / ( wavespd_ein_r ( i)
     1                        - wavespd_ein_l ( i))

      temparr_02 ( i) = wavespd_ein_r ( i) * wavespd_ein_l ( i)
     1       / ( wavespd_ein_r ( i) - wavespd_ein_l ( i))


      entropy_flux_re ( i) = 

     1     ( wavespd_ein_r ( i) * entropy_flux_ein_l ( i)
     1     - wavespd_ein_l ( i) * entropy_flux_ein_r ( i) )
     1     * temparr_01 ( i)

     1     + ( entropy_ein_r ( i) - entropy_ein_l ( i) )
     1     * temparr_02 ( i)

      END DO

******************************************************************
******************************************************************

      END IF

* End of ientropyfix = 1 or 2 case.

******************************************************************
******************************************************************

! Build "do_ein" and "ratio_ein" so that one can treat strong
! shock situations by blending in parts of an HLLE riemann solver.
! There are two alternatives, one which treats strong rarefactions
! on the same footing as shocks; the other which only goes after shocks.

      DO i = nriembegin, nriemend

      temparr_01 ( i) = AMIN1 ( msonic_arr ( i),
     1                          eigenval_lf ( i, 5) - vlxremi ( i),
     1                          eigenval_rt ( i, 5) - vlxrepl ( i) )

#if ( 1 == 1)
! Here we treat strong rarefactions on an equal footing with strong shocks.
      temparr_02 ( i) = AMAX1 ( ABS ( divvel_arr ( i) ),
     1                          ABS ( vlxrepl ( i) - vlxremi ( i) ),
     1            ABS ( eigenval_rt ( i, 5) - eigenval_lf ( i, 5) ),
     1            ABS ( eigenval_rt ( i, 1) - eigenval_lf ( i, 1) ) )

      temparr_03 ( i) = - temparr_02 ( i)
     1                + divratio_einfeldt_rs * temparr_01 ( i)
#else
! Here we only try to detect strong shocks. We don't detect strong rarefactions.
      temparr_02 ( i) = AMIN1 ( divvel_arr ( i),
     1                          vlxrepl ( i) - vlxremi ( i),
     1            eigenval_rt ( i, 5) - eigenval_lf ( i, 5),
     1            eigenval_rt ( i, 1) - eigenval_lf ( i, 1) )

      temparr_03 ( i) = temparr_02 ( i)
     1                + divratio_einfeldt_rs * temparr_01 ( i)
#endif

      IF ( temparr_03 ( i) .GE. 0.0) THEN
        do_ein ( i) = 0
        ratio_ein ( i) = 0.0
      ELSE
        do_ein ( i) = 1
        ratio_ein ( i) = ABS ( temparr_03 ( i))
     1                   / ( divratio_einfeldt_rs * temparr_01 ( i))
        ratio_ein ( i) = AMIN1 ( ratio_ein ( i), 1.0)
      END IF

      END DO

******************************************************************
******************************************************************

! Now check to see if each zone boundary has positive density
! and pressure. If not, give that zone an Einfeldt fix.

! For multispecies we have to also approximate the value of "gamma".
! This is done in the most reasonable way that is relatively inexpensive.
! i.e. it is not the scientifically best choice, but it is
! the most expedient one.

! For a more traditional interpretation of the Einfeldt fix, uncomment
! the ensuing do loop.
!      DO i = nriembegin, nriemend
!      do_ein ( i) = 0
!      ratio_ein ( i) = 0.0
!      END DO


      DO i = nriembegin, nriemend

! Form the ratio of specific heats etc.

      IF ( ( eigenval_db ( i, 0, 2) .GT. 0.0) .AND.
     1     ( eigenval_db ( i, 1, 2) .GT. 0.0) ) THEN

      mean_mol_wt_re ( i) = mean_mol_wt_remi ( i)
      gamma_eos_re ( i) = gamma_eos_remi ( i)
      gamma_soundspeed_re ( i) = gamma_soundspeed_remi ( i)

      ELSE IF ( ( eigenval_db ( i, 0, 2) .LT. 0.0) .AND.
     1          ( eigenval_db ( i, 1, 2) .LT. 0.0) ) THEN

      mean_mol_wt_re ( i) = mean_mol_wt_repl ( i)
      gamma_eos_re ( i) = gamma_eos_repl ( i)
      gamma_soundspeed_re ( i) = gamma_soundspeed_repl ( i)

      ELSE

      mean_mol_wt_re ( i) = mean_mol_wt ( i)
      gamma_eos_re ( i) = gamma_eos ( i)
      gamma_soundspeed_re ( i) = gamma_soundspeed ( i)

      END IF


! Form the flow variables at the resolved state.

      rho0 ( i) = AMAX1 ( conserved_re ( i, 1), rhofloor)

      temparr_01 ( i) = 1.0 / rho0 ( i)
      vlx0 ( i) = conserved_re ( i, 2) * temparr_01 ( i)
      vly0 ( i) = conserved_re ( i, 3) * temparr_01 ( i)
      vlz0 ( i) = conserved_re ( i, 4) * temparr_01 ( i)


      prs0 ( i) = conserved_re ( i, 5)
     1   - 0.5 * rho0 ( i) * ( vlx0 ( i)**2 + vly0 ( i)**2
     1                       + vlz0 ( i)**2 )
      prs0 ( i) = ( gamma_eos_re ( i) - 1.0) * prs0 ( i)


      IF ( NFLUID .GE. 1) THEN

      temparr_01 ( i) = 1.0 / rho0 ( i)
      specie0 ( 1: NFLUID, i) =
     1 conserved_re ( i, NFIELD + 1: NFIELD + NFLUID)
     1 * temparr_01 ( i)

      END IF
! End of "NFLUID .GE. 1".


! Flag zones for Einfeldt fix, if that is needed.

      IF ( ( rho0 ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prs0 ( i) .LE. prsfloor_plus_delta) ) THEN
        do_ein ( i) = 1
        ratio_ein ( i) = 1.0
      END IF

      IF ( ( rhoremi ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prsremi ( i) .LE. prsfloor_plus_delta) ) THEN
        do_ein ( i) = 1
        ratio_ein ( i) = 1.0
      END IF

      IF ( ( rhorepl ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prsrepl ( i) .LE. prsfloor_plus_delta) ) THEN
        do_ein ( i) = 1
        ratio_ein ( i) = 1.0
      END IF

      IF ( ( rhobar ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prsbar ( i) .LE. prsfloor_plus_delta) ) THEN
        do_ein ( i) = 1
        ratio_ein ( i) = 1.0
      END IF

      END DO


* If the species fractions in the resolved state are less than 0.0
* or greater than 1.0 then too we need to invoke the HLL riemann solver.
* Because it is based on a linear, convex arithmetic sum of the left
* and right states, the HLL riemann solver can be relied on to
* produce species fractions that satisfy this physical requirement
* of lying between 0 and 1. This property also ensures that if both
* the right and left states of the riemann problem have the property
* that the sum of the active species fractions remains less than one then 
* the resolved state from the riemann solver also has that property.

      DO ifluid = 1, NFLUID

        DO i = nriembegin, nriemend

        IF ( ( specie0 ( ifluid, i) .LT. - smallnum) .OR.
     1       ( specie0 ( ifluid, i) .GT. 1.0 + smallnum) ) THEN
          do_ein ( i) = 1
          ratio_ein ( i) = 1.0
        END IF

        END DO

      END DO
*
*
* Now make index array with which we will access those zones that need
* the Einfeldt fix.
*
      num_ein = nriembegin - 1

      DO i = nriembegin, nriemend

      IF ( do_ein ( i) .EQ. 1) THEN

      num_ein = num_ein + 1

      indx_ein ( num_ein) = i

      END IF

      END DO

! Note that we have commented out the warning when Einfeldt fix is invoked.

!      IF ( num_ein .GE. nriembegin )
!     1      WRITE ( 6, *)" ntstep = ", ntstep,
!     1       " gives num_ein = ", num_ein - nriembegin + 1

*******************************************************************

      IF ( num_ein .GE. nriembegin) THEN

*******************************************************************

! First gather zones that need Einfeldt fix. 
! The bounding eigenvalues are those produced from the
! Roe solver.
! Since the Roe solver always produces a good
! value for the Alfven speed we include it in here. 

      DO i = nriembegin, num_ein

****************

      ii = indx_ein ( i)

****************

      conserved_ein_l ( i, 1: NFIELD + NFLUID) =
     1    zonebdy_conserved ( ii, 0, 1: NFIELD + NFLUID)

      conserved_ein_r ( i, 1: NFIELD + NFLUID) =
     1    zonebdy_conserved ( ii, 1, 1: NFIELD + NFLUID)


      flux_ein_l ( i, 1: NFIELD + NFLUID) =
     1    zonebdy_flux ( ii, 0, 1: NFIELD + NFLUID)

      flux_ein_r ( i, 1: NFIELD + NFLUID) =
     1    zonebdy_flux ( ii, 1, 1: NFIELD + NFLUID)


      wavespd_ein_l ( i) = AMIN1 ( eigenval_db ( ii, 0, 1),
     1                             eigenval_db ( ii, 0, 2) )

      wavespd_ein_r ( i) = AMAX1 ( eigenval_db ( ii, 1, 5),
     1                             eigenval_db ( ii, 1, 4) )

****************

! The lines below make the flux into an LLF flux (which is more stabilizing).
! Comment them out to get an Einfeldt flux.

      temparr_01 ( i) = AMAX1 ( ABS ( eigenval_db ( ii, 0, 1)), 
     1                          ABS ( eigenval_db ( ii, 0, 2)), 
     1                          ABS ( eigenval_db ( ii, 1, 5)),
     1                          ABS ( eigenval_db ( ii, 1, 4)))
      wavespd_ein_l ( i) = - 1.0 * temparr_01 ( i)
      wavespd_ein_r ( i) = 1.0 * temparr_01 ( i)

****************

! Here we make the Einfeldt fix directly without invoking a subroutine call.

      IF ( wavespd_ein_l ( i) .GT. 0.0) THEN
! Pick left conserved state.

      conserved_ein ( i, 1: NFIELD + NFLUID) =
     1      conserved_ein_l ( i, 1: NFIELD + NFLUID)

      ELSE IF ( wavespd_ein_r ( i) .LT. 0.0) THEN
! Pick right conserved state.

      conserved_ein ( i, 1: NFIELD + NFLUID) =
     1      conserved_ein_r ( i, 1: NFIELD + NFLUID)

      ELSE
! Pick HLL with one mean, conserved state.

      temparr_01 ( i) = 1.0 / ( wavespd_ein_r ( i)
     1                        - wavespd_ein_l ( i))

      conserved_ein ( i, 1: NFIELD + NFLUID) =

     1  ( wavespd_ein_r ( i)
     1  * conserved_ein_r ( i, 1: NFIELD + NFLUID)
     1  -  wavespd_ein_l ( i)
     1  * conserved_ein_l ( i, 1: NFIELD + NFLUID) )
     1  * temparr_01 ( i)

     1  - ( flux_ein_r ( i, 1: NFIELD + NFLUID)
     1    - flux_ein_l ( i, 1: NFIELD + NFLUID) )
     1  * temparr_01 ( i)

      END IF

****************

! Introduce just one mean state or, alternatively, pick either the
! left or right states.

      wavespd_ein_l ( i) = AMIN1 ( wavespd_ein_l ( i), 0.0)
      wavespd_ein_r ( i) = AMAX1 ( wavespd_ein_r ( i), 0.0)

      temparr_01 ( i) = 1.0 / ( wavespd_ein_r ( i)
     1                        - wavespd_ein_l ( i))

      temparr_02 ( i) = wavespd_ein_r ( i) * wavespd_ein_l ( i)
     1       / ( wavespd_ein_r ( i) - wavespd_ein_l ( i))


      flux_ein ( i, 1: NFIELD + NFLUID) =

     1 ( wavespd_ein_r ( i)
     1 * flux_ein_l ( i, 1: NFIELD + NFLUID)
     1 - wavespd_ein_l ( i)
     1 * flux_ein_r ( i, 1: NFIELD + NFLUID) )
     1 * temparr_01 ( i)

     1 + ( conserved_ein_r ( i, 1: NFIELD + NFLUID)
     1   - conserved_ein_l ( i, 1: NFIELD + NFLUID) )
     1 * temparr_02 ( i)

****************

! And hand them back. This is where we blend "1.0 - ratio_ein ( ii)" of the
! Roe type Riemann solver with "ratio_ein ( ii)" of the HLLE Riemann solver.

      temparr_02 ( i) = 1.0 - ratio_ein ( ii)

      conserved_re ( ii, 1: NFIELD + NFLUID) =
     1   temparr_02 ( i)
     1 * conserved_re ( ii, 1: NFIELD + NFLUID)
     1 + ratio_ein ( ii)
     1 * conserved_ein ( i, 1: NFIELD + NFLUID)

      flux_re ( ii, 1: NFIELD + NFLUID) = 
     1   temparr_02 ( i)
     1 * flux_re ( ii, 1: NFIELD + NFLUID)
     1 + ratio_ein ( ii)
     1 * flux_ein ( i, 1: NFIELD + NFLUID)

      END DO

******************************************************************

      END IF

! End of making Einfeldt fix.

******************************************************************
******************************************************************

! Wind-up operations done past this point:

******************************************************************

! Move the species fluxes by one to make way for the split
! in the advected and pressure parts of the x-momentum flux.

! To do that, we realize first that the 6th component of the flux will
! have to contain the pressure. Right now, the fluxes are all contiguous
! even when multiple species are present. So we first move the fluxes
! over by one when multiple species are present.

      DO ifluid = NFLUID, 1, -1

      k = NFIELD + ifluid

        DO i = nriembegin, nriemend

        flux_re ( i, k + 1) = flux_re ( i, k)

        END DO

      END DO

! Zero out the sixth component of "flux_re" because we have not
! yet made the pressure split.

      DO i = nriembegin, nriemend

      flux_re ( i, 6) = 0.0

      END DO

******************************************************************
******************************************************************

! Form the flow variables at the resolved state.
! Also split the momentum flux_re into advected and pressure parts and
! store the pressure part separately.

      IF ( split_prsflux_re .EQ. 1) THEN

      DO i = nriembegin, nriemend

      rho0 ( i) = AMAX1 ( conserved_re ( i, 1), rhofloor)

      temparr_01 ( i) = 1.0 / rho0 ( i)
      vlx0 ( i) = conserved_re ( i, 2) * temparr_01 ( i)
      vly0 ( i) = conserved_re ( i, 3) * temparr_01 ( i)
      vlz0 ( i) = conserved_re ( i, 4) * temparr_01 ( i)


      prs0 ( i) = conserved_re ( i, 5)
     1   - 0.5 * rho0 ( i) * ( vlx0 ( i)**2 + vly0 ( i)**2
     1                       + vlz0 ( i)**2 )
      prs0 ( i) = ( gamma_eos_re ( i) - 1.0) * prs0 ( i)


      IF ( NFLUID .GE. 1) THEN

      temparr_01 ( i) = 1.0 / rho0 ( i)
      specie0 ( 1: NFLUID, i) =
     1 conserved_re ( i, NFIELD + 1: NFIELD + NFLUID)
     1 * temparr_01 ( i)

      END IF
! End of "NFLUID .GE. 1".


! Make the part of the momentum flux that acts like a pressure.

      flux_re ( i, 6) = prs0 ( i)

      flux_re ( i, 2) = flux_re ( i, 2) - flux_re ( i, 6)

      END DO

      END IF
! End of "split_prsflux_re .EQ. 1" part.

**********************************************************************
**********************************************************************

      RETURN

      END SUBROUTINE EULER_RIEM_ROE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! We have the option of minimizing the number of EOS calls
! (which can be quite expensive) by setting "use_approx_eos = 1".
! This can be a risky step if "gamma_eos" is rapidly varying.

      SUBROUTINE EULER_RIEM_HLL

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

* Written by : Dinshaw S. Balsara

* Does the HLL riemann solver. The interface and inputs/outputs are
* intentionally kept the same as the linearized Riemann solver so that
* the subroutines can be used interchangeably.

* The x-momentum flux_re is split into an advected part and a ( thermal +
* magnetic) pressure part.

* ientropyfix = 1 will give the maximal dissipation by using the LLF flux
* even for characteristic fields that go into eachother.

* ientropyfix = 2 will provide an entropy fix when characteristic fields are
* going apart but not when they are going into eachother.

* ientropyfix = 3 will provide the basic Harten and Hyman entropy fix.

* While the Roe type riemann solver is little more accurate than the
* HLL solver, it can become very expensive when several species
* are involved.

*******************************************************************

      IMPLICIT NONE

      EXTERNAL EIGENVECMHD, RHO_PRS_IN_ENG_T_GAMA_OUT_EULER

*******************************************************************

      INTEGER :: nriembegin, nriemend, 
     1        ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,
     1        ientropyfix, iprimitive, use_approx_eos,
     1        split_prsflux_re

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor, 
     1     fracdiff, gasconst, divratio_einfeldt_rs

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr

      REAL, DIMENSION
     1   ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD + NFLUID, NFIELD + NFLUID) ::
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
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION
     1    ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

************************************************************************
************************************************************************

* Past this point we only declare variables that are internal to the
* subroutine.

      INTEGER i, ii, j, k, ifield, ifluid, num_ein

      REAL small_speed,
     1     rhofloor_plus_delta, prsfloor_plus_delta

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID) ::

     1      zonebdy_eigenval, eigenval_lf, eigenval_rt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD + NFLUID) ::

     1      eigenval_db

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     temparr_01, temparr_02, temparr_03


******************************************************************
******************************************************************

* Variables used for Einfeldt fix.

* These should be loaded with the min and max speeds from the roe solver.
      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     wavespd_ein_l, wavespd_ein_r

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     eng, temperature,
     1     mean_mol_wt, gamma_eos, gamma_soundspeed,
     1     engremi, temperatureremi, total_enthalpy_remi,
     1     engrepl, temperaturerepl, total_enthalpy_repl,
     1     mean_mol_wt_re, gamma_eos_re, gamma_soundspeed_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     c_s,
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, total_enthalpy_bar

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      speciebar

! ----------------------------------------------------------------------

* Make some settings.
* Because floor values of density and pressure have already been taken,
* we set a bottom value that is a little above that floor.


      rhofloor_plus_delta = 1.000001 * rhofloor
      prsfloor_plus_delta = 1.000001 * prsfloor


      IF ( NFIELD .NE. 5) THEN
      WRITE ( 6, *)"EULER_RIEM_HLL : NFIELD .NE. 5", NFIELD,5
      debug_stop = 1
      RETURN
      END IF

      IF ( ( nriembegin .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriembegin .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLL : nriembegin out of bounds",
     1       nriembegin, ONEDFLUIDBEGIN, ONEDFLUIDEND
      debug_stop = 1
      RETURN
      END IF

      IF ( ( nriemend .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriemend .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLL : nriemend out of bounds",
     1       nriemend, ONEDFLUIDBEGIN, ONEDFLUIDEND
      debug_stop = 1
      RETURN
      END IF

      IF ( ( iprimitive .LT. -1) .OR.
     1     ( iprimitive .GT. 0) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLL : iprimitive out of bounds",
     1  iprimitive
      debug_stop = 1
      RETURN
      END IF

      IF ( ( ientropyfix .LT. 1) .OR.
     1     ( ientropyfix .GT. 3) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLL : ientropyfix out of bounds",
     1  ientropyfix
      debug_stop = 1
      RETURN
      END IF

      IF ( ( split_prsflux_re .LT. 0) .OR.
     1     ( split_prsflux_re .GT. 1)) THEN
      WRITE ( 6, *)"EULER_RIEM_HLL : split_prsflux_re ",
     1             " out of bounds", split_prsflux_re
      debug_stop = 1
      RETURN
      END IF


******************************************************************
******************************************************************

* Bound the density and pressure variables.
* For starters none of the zones need Einfeldt fix.
* Also "rhobar 2" and "prsbar 2" may not be evaluated in all cases
* later on so we give them a dummy initialization.

      DO i = nriembegin, nriemend

      do_ein ( i) = 0

      rhoremi ( i) = AMAX1 ( rhoremi ( i), rhofloor)
      rhorepl ( i) = AMAX1 ( rhorepl ( i), rhofloor)

      prsremi ( i) = AMAX1 ( prsremi ( i), prsfloor)
      prsrepl ( i) = AMAX1 ( prsrepl ( i), prsfloor)
 
      END DO

******************************************************************
******************************************************************

* Evaluate the thermal energy density and the local gammas. In doing so, we
* also make sure the species fractions are individually and cumulatively
* within reasonable ranges. Do this only if exact EOS's are to be used.

      IF ( use_approx_eos .NE. 1) THEN

      WRITE ( 6, *)
     1      "Stopping because EOS routines not present in demo. code."
      STOP

      ELSE

* Relying on approximate EOS's here.

      DO i = nriembegin, nriemend

      engremi ( i) = prsremi ( i) / ( gamma_eos_remi ( i) - 1.0)
      temperatureremi ( i) = prsremi ( i) * mean_mol_wt_remi ( i)
     1                     / ( gasconst * rhoremi ( i))

      engrepl ( i) = prsrepl ( i) / ( gamma_eos_repl ( i) - 1.0)
      temperaturerepl ( i) = prsrepl ( i) * mean_mol_wt_repl ( i)
     1                     / ( gasconst * rhorepl ( i))

      END DO

      END IF

***************************

* Make "total_enthalpy_remi" and "total_enthalpy_repl". These are
* useful in making Roe-averages.

      DO i = nriembegin, nriemend

      total_enthalpy_remi ( i) = ( engremi ( i) + prsremi ( i)
     1    + 0.5 * rhoremi ( i) * ( vlxremi ( i)**2 + vlyremi ( i)**2
     1                           + vlzremi ( i)**2) ) / rhoremi ( i)

      total_enthalpy_repl ( i) = ( engrepl ( i) + prsrepl ( i)
     1    + 0.5 * rhorepl ( i) * ( vlxrepl ( i)**2 + vlyrepl ( i)**2
     1                           + vlzrepl ( i)**2) ) / rhorepl ( i)

      END DO

******************************************************************
******************************************************************

* Make "zonebdy_conserved" and "zonebdy_flux". Uncommenting this part would
* make this a regular linearized Riemann solver. We only need to comment
* it out when using ADER schemes because they may produce higher order fluxes
* and conserved variables for use in the Riemann solver.

      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 1) = rhoremi ( i)
      zonebdy_conserved ( i, 0, 2) = rhoremi ( i)
     1                             * vlxremi ( i)
      zonebdy_conserved ( i, 0, 3) = rhoremi ( i)
     1                             * vlyremi ( i)
      zonebdy_conserved ( i, 0, 4) = rhoremi ( i)
     1                             * vlzremi ( i)


      zonebdy_conserved ( i, 1, 1) = rhorepl ( i)
      zonebdy_conserved ( i, 1, 2) = rhorepl ( i)
     1                             * vlxrepl ( i)
      zonebdy_conserved ( i, 1, 3) = rhorepl ( i)
     1                             * vlyrepl ( i)
      zonebdy_conserved ( i, 1, 4) = rhorepl ( i)
     1                             * vlzrepl ( i)


      zonebdy_flux ( i, 0, 1) = zonebdy_conserved ( i, 0, 2)
      zonebdy_flux ( i, 0, 2) =
     1     zonebdy_conserved ( i, 0, 2) * vlxremi ( i)
     1   + prsremi ( i)
      zonebdy_flux ( i, 0, 3) =
     1     zonebdy_conserved ( i, 0, 3) * vlxremi ( i)
      zonebdy_flux ( i, 0, 4) =
     1     zonebdy_conserved ( i, 0, 4) * vlxremi ( i)


      zonebdy_flux ( i, 1, 1) = zonebdy_conserved ( i, 1, 2)
      zonebdy_flux ( i, 1, 2) =
     1     zonebdy_conserved ( i, 1, 2) * vlxrepl ( i)
     1   + prsrepl ( i)
      zonebdy_flux ( i, 1, 3) =
     1     zonebdy_conserved ( i, 1, 3) * vlxrepl ( i)
      zonebdy_flux ( i, 1, 4) =
     1     zonebdy_conserved ( i, 1, 4) * vlxrepl ( i)


      END DO


* Make the fifth component differently depending on which flux
* is desired.

      IF ( iprimitive .EQ. 0) THEN

      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 5) =
     1      engremi ( i)
     1    + 0.5 * rhoremi ( i)
     1    * ( vlxremi ( i)**2 + vlyremi ( i)**2
     1                        + vlzremi ( i)**2)

      zonebdy_conserved ( i, 1, 5) =
     1      engrepl ( i)
     1    + 0.5 * rhorepl ( i)
     1    * ( vlxrepl ( i)**2 + vlyrepl ( i)**2
     1                        + vlzrepl ( i)**2)

      zonebdy_flux ( i, 0, 5) =
     1     ( zonebdy_conserved ( i, 0, 5) + prsremi ( i) )
     1   * vlxremi ( i)

      zonebdy_flux ( i, 1, 5) =
     1     ( zonebdy_conserved ( i, 1, 5) + prsrepl ( i) )
     1   * vlxrepl ( i)

      END DO

      END IF


      IF ( iprimitive .EQ. - 1) THEN

* EOS
      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 5) = prsremi ( i)
     1        / rhoremi ( i)**( gamma_soundspeed_remi ( i) - 1.0)

      zonebdy_conserved ( i, 1, 5) = prsrepl ( i)
     1        / rhorepl ( i)**( gamma_soundspeed_repl ( i) - 1.0)

      zonebdy_flux ( i, 0, 5) = zonebdy_conserved ( i, 0, 5)
     1                        * vlxremi ( i)

      zonebdy_flux ( i, 1, 5) = zonebdy_conserved ( i, 1, 5)
     1                        * vlxrepl ( i)

      END DO

      END IF


* If multiple species are present then initialize their fluxes etc.
* appropriately also.

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        zonebdy_conserved ( i, 0, k) = rhoremi ( i)
     1                               * specieremi ( ifluid, i)

        zonebdy_conserved ( i, 1, k) = rhorepl ( i)
     1                               * specierepl ( ifluid, i)

        zonebdy_flux ( i, 0, k) = zonebdy_conserved ( i, 0, k)
     1                          * vlxremi ( i)

        zonebdy_flux ( i, 1, k) = zonebdy_conserved ( i, 1, k)
     1                          * vlxrepl ( i)

        conserved_re ( i, k) = 0.0

        flux_re ( i, k + 1) = 0.0

      END DO

      END DO

******************************************************************
******************************************************************

! Zero out "flux_re" here.

      DO i = nriembegin, nriemend

      conserved_re ( i, 1) = 0.0
      conserved_re ( i, 2) = 0.0
      conserved_re ( i, 3) = 0.0
      conserved_re ( i, 4) = 0.0
      conserved_re ( i, 5) = 0.0

      flux_re ( i, 1) = 0.0
      flux_re ( i, 2) = 0.0
      flux_re ( i, 3) = 0.0
      flux_re ( i, 4) = 0.0
      flux_re ( i, 5) = 0.0
      flux_re ( i, 6) = 0.0

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        conserved_re ( i, k) = 0.0

        flux_re ( i, k + 1) = 0.0

      END DO

      END DO


******************************************************************
******************************************************************

* Here we use Roe averaging and make only one call to the
* eigenvector subroutine. The "eng, temperature" etc. that we make here
* can even serve as good initial guesses for an exact EOS.


      DO i = nriembegin, nriemend

! ----------

      gamma_eos ( i) = 0.5 * ( gamma_eos_remi ( i)
     1                       + gamma_eos_repl ( i))
      gamma_soundspeed ( i) = 0.5 * ( gamma_soundspeed_remi ( i)
     1                              + gamma_soundspeed_repl ( i))
      mean_mol_wt ( i) = 0.5 * ( mean_mol_wt_remi ( i)
     1                         + mean_mol_wt_repl ( i))

! ----------

      temparr_01 ( i) = SQRT ( rhorepl ( i) / rhoremi ( i))
      temparr_02 ( i) = 1.0 / ( 1.0 + temparr_01 ( i))

      rhobar ( i) = SQRT ( rhoremi ( i) * rhorepl ( i))

      vlxbar ( i) = ( vlxremi ( i) + vlxrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)
      vlybar ( i) = ( vlyremi ( i) + vlyrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)
      vlzbar ( i) = ( vlzremi ( i) + vlzrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)

      total_enthalpy_bar ( i) = ( total_enthalpy_remi ( i)
     1  + total_enthalpy_repl ( i) * temparr_01 ( i)) * temparr_02 ( i)

      prsbar ( i) = total_enthalpy_bar ( i) - 0.5 * ( vlxbar ( i)**2
     1                     + vlybar ( i)**2 + vlzbar ( i)**2)
      prsbar ( i) = ( gamma_eos ( i) - 1.0) * rhobar ( i) * prsbar ( i)
     1            / gamma_eos ( i)
      prsbar ( i) = AMAX1 ( prsbar ( i), prsfloor)

! ----------

! Give yourself the option of using arithmetic averages if needed.

      IF ( ( rhobar ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prsbar ( i) .LE. prsfloor_plus_delta) ) THEN

      rhobar ( i) = 0.5 * ( rhoremi ( i) + rhorepl ( i))
      prsbar ( i) = 0.5 * ( prsremi ( i) + prsrepl ( i))

      vlxbar ( i) = 0.5 * ( vlxremi ( i) + vlxrepl ( i))
      vlybar ( i) = 0.5 * ( vlyremi ( i) + vlyrepl ( i))
      vlzbar ( i) = 0.5 * ( vlzremi ( i) + vlzrepl ( i))

      END IF

! ----------

      eng ( i) = prsbar ( i) / ( gamma_eos ( i) - 1.0)
      temperature ( i) = prsbar ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rhobar ( i))

! ----------

      END DO

! ------------------------------

      DO ifluid = 1, NFLUID

      DO i = nriembegin, nriemend

        speciebar ( ifluid, i) = 0.5 * ( specieremi ( ifluid, i)
     1                                 + specierepl ( ifluid, i))

      END DO

      END DO


* Using exact EOS's here. If that is not needed, then "eng, temperature"
* etc. have already been made above and the present call is not needed.

      IF ( use_approx_eos .NE. 1) THEN

      WRITE ( 6, *)
     1      "Stopping because EOS routines not present in demo. code."

      END IF

*********************************************************************
*********************************************************************

* Make "eigenval_re" by using the Roe-averaged values:

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed ( i) * prsbar ( i) / rhobar ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_re ( i, 1) = vlxbar ( i) - c_s ( i)
      eigenval_re ( i, 2) = vlxbar ( i)
      eigenval_re ( i, 3) = vlxbar ( i)
      eigenval_re ( i, 4) = vlxbar ( i)
      eigenval_re ( i, 5) = vlxbar ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_re ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxbar ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_rt"

*********************************************************************
*********************************************************************

      IF ( debug_stop .EQ. 1) THEN
      RETURN
      END IF


* And store the eigenvalues on either side of
* the zone boundary.

      DO i = nriembegin, nriemend

      zonebdy_eigenval ( i, 1) = eigenval_re ( i, 1)
      eigenval_db ( i, 0, 1) = eigenval_re ( i, 1)
      eigenval_db ( i, 1, 1) = eigenval_re ( i, 1)

      zonebdy_eigenval ( i, 2) = eigenval_re ( i, 2)
      eigenval_db ( i, 0, 2) = eigenval_re ( i, 2)
      eigenval_db ( i, 1, 2) = eigenval_re ( i, 2)

      zonebdy_eigenval ( i, 3) = eigenval_re ( i, 3)
      eigenval_db ( i, 0, 3) = eigenval_re ( i, 3)
      eigenval_db ( i, 1, 3) = eigenval_re ( i, 3)

      zonebdy_eigenval ( i, 4) = eigenval_re ( i, 4)
      eigenval_db ( i, 0, 4) = eigenval_re ( i, 4)
      eigenval_db ( i, 1, 4) = eigenval_re ( i, 4)

      zonebdy_eigenval ( i, 5) = eigenval_re ( i, 5)
      eigenval_db ( i, 0, 5) = eigenval_re ( i, 5)
      eigenval_db ( i, 1, 5) = eigenval_re ( i, 5)

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        zonebdy_eigenval ( i, k) = eigenval_re ( i, k)
        eigenval_db ( i, 0, k) = eigenval_re ( i, k)
        eigenval_db ( i, 1, k) = eigenval_re ( i, k)

      END DO

      END DO

*********************************************************************
*********************************************************************

* Here we explicitly need to get the eigenvalues on either
* side of the zone boundary so that we can maximize / minimize w.r.t.
* the eigenvalues of the mean state. We reset "eigenval_db" here
* so that "eigenval_db ( 0, ifield, i)" is the left-most eigenvalue
* in that characteristic field and "eigenval_db ( 1, ifield, i)" is
* the right-most.

* Make "eigenval_lf":

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed_remi ( i) * prsremi ( i)
     1         / rhoremi ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_lf ( i, 1) = vlxremi ( i) - c_s ( i)
      eigenval_lf ( i, 2) = vlxremi ( i)
      eigenval_lf ( i, 3) = vlxremi ( i)
      eigenval_lf ( i, 4) = vlxremi ( i)
      eigenval_lf ( i, 5) = vlxremi ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_lf ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxremi ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_lf"

*********************************************************************
*********************************************************************

      DO i = nriembegin, nriemend


      eigenval_db ( i, 0, 1) =
     1      AMIN1 ( eigenval_db ( i, 0, 1), eigenval_lf ( i, 1) )

      eigenval_db ( i, 1, 1) =
     1      AMAX1 ( eigenval_db ( i, 1, 1), eigenval_lf ( i, 1) )


      eigenval_db ( i, 0, 2) =
     1      AMIN1 ( eigenval_db ( i, 0, 2), eigenval_lf ( i, 2) )

      eigenval_db ( i, 1, 2) =
     1      AMAX1 ( eigenval_db ( i, 1, 2), eigenval_lf ( i, 2) )


      eigenval_db ( i, 0, 3) =
     1      AMIN1 ( eigenval_db ( i, 0, 3), eigenval_lf ( i, 3) )

      eigenval_db ( i, 1, 3) =
     1      AMAX1 ( eigenval_db ( i, 1, 3), eigenval_lf ( i, 3) )


      eigenval_db ( i, 0, 4) =
     1      AMIN1 ( eigenval_db ( i, 0, 4), eigenval_lf ( i, 4) )

      eigenval_db ( i, 1, 4) =
     1      AMAX1 ( eigenval_db ( i, 1, 4), eigenval_lf ( i, 4) )


      eigenval_db ( i, 0, 5) =
     1      AMIN1 ( eigenval_db ( i, 0, 5), eigenval_lf ( i, 5) )

      eigenval_db ( i, 1, 5) =
     1      AMAX1 ( eigenval_db ( i, 1, 5), eigenval_lf ( i, 5) )


      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        eigenval_db ( i, 0, k) =
     1        AMIN1 ( eigenval_db ( i, 0, k), eigenval_lf ( i, k) )

        eigenval_db ( i, 1, k) =
     1        AMAX1 ( eigenval_db ( i, 1, k), eigenval_lf ( i, k) )

      END DO

      END DO

*********************************************************************
*********************************************************************

* Make "eigenval_rt":

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed_repl ( i) * prsrepl ( i)
     1         / rhorepl ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_rt ( i, 1) = vlxrepl ( i) - c_s ( i)
      eigenval_rt ( i, 2) = vlxrepl ( i)
      eigenval_rt ( i, 3) = vlxrepl ( i)
      eigenval_rt ( i, 4) = vlxrepl ( i)
      eigenval_rt ( i, 5) = vlxrepl ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_rt ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxrepl ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_rt"

*********************************************************************
*********************************************************************

      DO i = nriembegin, nriemend


      eigenval_db ( i, 0, 1) =
     1      AMIN1 ( eigenval_db ( i, 0, 1), eigenval_rt ( i, 1) )

      eigenval_db ( i, 1, 1) =
     1      AMAX1 ( eigenval_db ( i, 1, 1), eigenval_rt ( i, 1) )


      eigenval_db ( i, 0, 2) =
     1      AMIN1 ( eigenval_db ( i, 0, 2), eigenval_rt ( i, 2) )

      eigenval_db ( i, 1, 2) =
     1      AMAX1 ( eigenval_db ( i, 1, 2), eigenval_rt ( i, 2) )


      eigenval_db ( i, 0, 3) =
     1      AMIN1 ( eigenval_db ( i, 0, 3), eigenval_rt ( i, 3) )

      eigenval_db ( i, 1, 3) =
     1      AMAX1 ( eigenval_db ( i, 1, 3), eigenval_rt ( i, 3) )


      eigenval_db ( i, 0, 4) =
     1      AMIN1 ( eigenval_db ( i, 0, 4), eigenval_rt ( i, 4) )

      eigenval_db ( i, 1, 4) =
     1      AMAX1 ( eigenval_db ( i, 1, 4), eigenval_rt ( i, 4) )


      eigenval_db ( i, 0, 5) =
     1      AMIN1 ( eigenval_db ( i, 0, 5), eigenval_rt ( i, 5) )

      eigenval_db ( i, 1, 5) =
     1      AMAX1 ( eigenval_db ( i, 1, 5), eigenval_rt ( i, 5) )

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        eigenval_db ( i, 0, k) =
     1        AMIN1 ( eigenval_db ( i, 0, k), eigenval_rt ( i, k) )

        eigenval_db ( i, 1, k) =
     1        AMAX1 ( eigenval_db ( i, 1, k), eigenval_rt ( i, k) )

      END DO

      END DO

* End of making the eigenvalues and appropriate wave speeds.

************************************************************************
************************************************************************

! First gather zones that need Einfeldt fix. 
! The bounding eigenvalues are those produced from the
! Roe solver.
! Since the Roe solver always produces a good
! value for the Alfven speed we include it in here. 

      DO i = nriembegin, nriemend

****************

      wavespd_ein_l ( i) = AMIN1 ( eigenval_db ( i, 0, 1),
     1                             eigenval_db ( i, 0, 2) )

      wavespd_ein_r ( i) = AMAX1 ( eigenval_db ( i, 1, 5),
     1                             eigenval_db ( i, 1, 4) )

****************

! The lines below make the flux into an LLF flux (which is more stabilizing).
! Comment them out to get an Einfeldt flux.

!      temparr_01 ( i) = AMAX1 ( ABS ( eigenval_db ( i, 0, 1)), 
!     1                          ABS ( eigenval_db ( i, 0, 2)), 
!     1                          ABS ( eigenval_db ( i, 1, 5)),
!     1                          ABS ( eigenval_db ( i, 1, 4)))
!      wavespd_ein_l ( i) = - 1.0 * temparr_01 ( i)
!      wavespd_ein_r ( i) = 1.0 * temparr_01 ( i)

****************

! Here we make the Einfeldt fix directly without invoking a subroutine call.

      IF ( wavespd_ein_l ( i) .GT. 0.0) THEN
! Pick left conserved state.

      conserved_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_conserved ( i, 0, 1: NFIELD + NFLUID)

      flux_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_flux ( i, 0, 1: NFIELD + NFLUID)

      ELSE IF ( wavespd_ein_r ( i) .LT. 0.0) THEN
! Pick right conserved state.

      conserved_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_conserved ( i, 1, 1: NFIELD + NFLUID)

      flux_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_flux ( i, 1, 1: NFIELD + NFLUID)

      ELSE
! Pick HLL with one mean, conserved state.

      temparr_01 ( i) = 1.0 / ( wavespd_ein_r ( i)
     1                        - wavespd_ein_l ( i))

      temparr_02 ( i) = wavespd_ein_r ( i) * wavespd_ein_l ( i)
     1       / ( wavespd_ein_r ( i) - wavespd_ein_l ( i))

      conserved_re ( i, 1: NFIELD + NFLUID) =

     1  ( wavespd_ein_r ( i)
     1  * zonebdy_conserved ( i, 1, 1: NFIELD + NFLUID)
     1  - wavespd_ein_l ( i)
     1  * zonebdy_conserved ( i, 0, 1: NFIELD + NFLUID) )
     1  * temparr_01 ( i)

     1  - ( zonebdy_flux ( i, 1, 1: NFIELD + NFLUID)
     1    - zonebdy_flux ( i, 0, 1: NFIELD + NFLUID) )
     1  * temparr_01 ( i)

      flux_re ( i, 1: NFIELD + NFLUID) =

     1 ( wavespd_ein_r ( i)
     1 * zonebdy_flux ( i, 0, 1: NFIELD + NFLUID)
     1 - wavespd_ein_l ( i)
     1 * zonebdy_flux ( i, 1, 1: NFIELD + NFLUID) )
     1 * temparr_01 ( i)

     1 + ( zonebdy_conserved ( i, 1, 1: NFIELD + NFLUID)
     1   - zonebdy_conserved ( i, 0, 1: NFIELD + NFLUID) )
     1 * temparr_02 ( i)

      END IF

      END DO

! End of making Einfeldt fix.

******************************************************************
******************************************************************

! Wind-up operations done past this point:

******************************************************************

! Move the species fluxes by one to make way for the split
! in the advected and pressure parts of the x-momentum flux.

! To do that, we realize first that the 6th component of the flux will
! have to contain the pressure. Right now, the fluxes are all contiguous
! even when multiple species are present. So we first move the fluxes
! over by one when multiple species are present.

      DO ifluid = NFLUID, 1, -1

      k = NFIELD + ifluid

        DO i = nriembegin, nriemend

        flux_re ( i, k + 1) = flux_re ( i, k)

        END DO

      END DO

! Zero out the sixth component of "flux_re" because we have not
! yet made the pressure split.

      DO i = nriembegin, nriemend

      flux_re ( i, 6) = 0.0

      END DO

******************************************************************
******************************************************************

! Form the flow variables at the resolved state.
! Also split the momentum flux_re into advected and pressure parts and
! store the pressure part separately.

      IF ( split_prsflux_re .EQ. 1) THEN

      DO i = nriembegin, nriemend

      rho0 ( i) = AMAX1 ( conserved_re ( i, 1), rhofloor)

      temparr_01 ( i) = 1.0 / rho0 ( i)
      vlx0 ( i) = conserved_re ( i, 2) * temparr_01 ( i)
      vly0 ( i) = conserved_re ( i, 3) * temparr_01 ( i)
      vlz0 ( i) = conserved_re ( i, 4) * temparr_01 ( i)


      prs0 ( i) = conserved_re ( i, 5)
     1   - 0.5 * rho0 ( i) * ( vlx0 ( i)**2 + vly0 ( i)**2
     1                       + vlz0 ( i)**2 )
      prs0 ( i) = ( gamma_eos_re ( i) - 1.0) * prs0 ( i)


      IF ( NFLUID .GE. 1) THEN

      temparr_01 ( i) = 1.0 / rho0 ( i)
      specie0 ( 1: NFLUID, i) =
     1 conserved_re ( i, NFIELD + 1: NFIELD + NFLUID)
     1 * temparr_01 ( i)

      END IF
! End of "NFLUID .GE. 1".


! Make the part of the momentum flux that acts like a pressure.

      flux_re ( i, 6) = prs0 ( i)

      flux_re ( i, 2) = flux_re ( i, 2) - flux_re ( i, 6)

      END DO

      END IF
! End of "split_prsflux_re .EQ. 1" part.

**********************************************************************
**********************************************************************

      RETURN

      END SUBROUTINE EULER_RIEM_HLL

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! We have the option of minimizing the number of EOS calls
! (which can be quite expensive) by setting "use_approx_eos = 1".
! This can be a risky step if "gamma_eos" is rapidly varying.

      SUBROUTINE EULER_RIEM_HLLC

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

* Written by : Dinshaw S. Balsara

* Does the HLLD riemann solver. The interface and inputs/outputs are
* intentionally kept the same as the linearized Riemann solver so that
* the subroutines can be used interchangeably. It is based on
* Toro, Spruce and Speares (1994) and Batten et al. (1997).

* The x-momentum flux_re is split into an advected part and a ( thermal +
* magnetic) pressure part.

* ientropyfix = 1 will give the maximal dissipation by using the LLF flux
* even for characteristic fields that go into eachother.

* ientropyfix = 2 will provide an entropy fix when characteristic fields are
* going apart but not when they are going into eachother.

* ientropyfix = 3 will provide the basic Harten and Hyman entropy fix.

* While the Roe type riemann solver is little more accurate than the
* HLLD solver, it can become very expensive when several species
* are involved.

*******************************************************************

      IMPLICIT NONE

      EXTERNAL EIGENVECMHD, RHO_PRS_IN_ENG_T_GAMA_OUT_EULER

*******************************************************************

      INTEGER :: nriembegin, nriemend, 
     1        ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,
     1        ientropyfix, iprimitive, use_approx_eos,
     1        split_prsflux_re

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor, 
     1     fracdiff, gasconst, divratio_einfeldt_rs

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr

      REAL, DIMENSION
     1   ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD + NFLUID, NFIELD + NFLUID) ::
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
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION
     1    ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

************************************************************************
************************************************************************

* Past this point we only declare variables that are internal to the
* subroutine.

      INTEGER i, ii, j, k, ifield, ifluid, num_ein

      REAL small_speed,
     1     rhofloor_plus_delta, prsfloor_plus_delta

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID) ::

     1      zonebdy_eigenval, eigenval_lf, eigenval_rt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD + NFLUID) ::

     1      eigenval_db

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     temparr_01, temparr_02, temparr_03


******************************************************************
******************************************************************

* Variables used for Einfeldt fix.

* These should be loaded with the min and max speeds from the roe solver.
      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     wavespd_ein_l, wavespd_ein_r

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     eng, temperature, eng0,
     1     mean_mol_wt, gamma_eos, gamma_soundspeed,
     1     engremi, temperatureremi, total_enthalpy_remi,
     1     engrepl, temperaturerepl, total_enthalpy_repl,
     1     mean_mol_wt_re, gamma_eos_re, gamma_soundspeed_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     c_s,
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, total_enthalpy_bar

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      speciebar

! ----------------------------------------------------------------------

* Make some settings.
* Because floor values of density and pressure have already been taken,
* we set a bottom value that is a little above that floor.


      rhofloor_plus_delta = 1.000001 * rhofloor
      prsfloor_plus_delta = 1.000001 * prsfloor


      IF ( NFIELD .NE. 5) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLC : NFIELD .NE. 5", NFIELD,5
      debug_stop = 1
      RETURN
      END IF

      IF ( ( nriembegin .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriembegin .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLC : nriembegin out of bounds",
     1       nriembegin, ONEDFLUIDBEGIN, ONEDFLUIDEND
      debug_stop = 1
      RETURN
      END IF

      IF ( ( nriemend .LT. ONEDFLUIDBEGIN) .OR.
     1     ( nriemend .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLC : nriemend out of bounds",
     1       nriemend, ONEDFLUIDBEGIN, ONEDFLUIDEND
      debug_stop = 1
      RETURN
      END IF

      IF ( ( iprimitive .LT. -1) .OR.
     1     ( iprimitive .GT. 0) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLC : iprimitive out of bounds",
     1  iprimitive
      debug_stop = 1
      RETURN
      END IF

      IF ( ( ientropyfix .LT. 1) .OR.
     1     ( ientropyfix .GT. 3) ) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLC : ientropyfix out of bounds",
     1  ientropyfix
      debug_stop = 1
      RETURN
      END IF

      IF ( ( split_prsflux_re .LT. 0) .OR.
     1     ( split_prsflux_re .GT. 1)) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLC : split_prsflux_re ",
     1             " out of bounds", split_prsflux_re
      debug_stop = 1
      RETURN
      END IF


******************************************************************
******************************************************************

* Bound the density and pressure variables.
* For starters none of the zones need Einfeldt fix.
* Also "rhobar 2" and "prsbar 2" may not be evaluated in all cases
* later on so we give them a dummy initialization.

      DO i = nriembegin, nriemend

      do_ein ( i) = 0

      rhoremi ( i) = AMAX1 ( rhoremi ( i), rhofloor)
      rhorepl ( i) = AMAX1 ( rhorepl ( i), rhofloor)

      prsremi ( i) = AMAX1 ( prsremi ( i), prsfloor)
      prsrepl ( i) = AMAX1 ( prsrepl ( i), prsfloor)
 
      END DO

******************************************************************
******************************************************************

* Evaluate the thermal energy density and the local gammas. In doing so, we
* also make sure the species fractions are individually and cumulatively
* within reasonable ranges. Do this only if exact EOS's are to be used.

      IF ( use_approx_eos .NE. 1) THEN

      WRITE ( 6, *)
     1      "Stopping because EOS routines not present in demo. code."
      STOP

      ELSE

* Relying on approximate EOS's here.

      DO i = nriembegin, nriemend

      engremi ( i) = prsremi ( i) / ( gamma_eos_remi ( i) - 1.0)
      temperatureremi ( i) = prsremi ( i) * mean_mol_wt_remi ( i)
     1                     / ( gasconst * rhoremi ( i))

      engrepl ( i) = prsrepl ( i) / ( gamma_eos_repl ( i) - 1.0)
      temperaturerepl ( i) = prsrepl ( i) * mean_mol_wt_repl ( i)
     1                     / ( gasconst * rhorepl ( i))

      END DO

      END IF

***************************

* Make "total_enthalpy_remi" and "total_enthalpy_repl". These are
* useful in making Roe-averages.

      DO i = nriembegin, nriemend

      total_enthalpy_remi ( i) = ( engremi ( i) + prsremi ( i)
     1    + 0.5 * rhoremi ( i) * ( vlxremi ( i)**2 + vlyremi ( i)**2
     1                           + vlzremi ( i)**2) ) / rhoremi ( i)

      total_enthalpy_repl ( i) = ( engrepl ( i) + prsrepl ( i)
     1    + 0.5 * rhorepl ( i) * ( vlxrepl ( i)**2 + vlyrepl ( i)**2
     1                           + vlzrepl ( i)**2) ) / rhorepl ( i)

      END DO

******************************************************************
******************************************************************

* Make "zonebdy_conserved" and "zonebdy_flux". Uncommenting this part would
* make this a regular linearized Riemann solver. We only need to comment
* it out when using ADER schemes because they may produce higher order fluxes
* and conserved variables for use in the Riemann solver.

      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 1) = rhoremi ( i)
      zonebdy_conserved ( i, 0, 2) = rhoremi ( i)
     1                             * vlxremi ( i)
      zonebdy_conserved ( i, 0, 3) = rhoremi ( i)
     1                             * vlyremi ( i)
      zonebdy_conserved ( i, 0, 4) = rhoremi ( i)
     1                             * vlzremi ( i)


      zonebdy_conserved ( i, 1, 1) = rhorepl ( i)
      zonebdy_conserved ( i, 1, 2) = rhorepl ( i)
     1                             * vlxrepl ( i)
      zonebdy_conserved ( i, 1, 3) = rhorepl ( i)
     1                             * vlyrepl ( i)
      zonebdy_conserved ( i, 1, 4) = rhorepl ( i)
     1                             * vlzrepl ( i)


      zonebdy_flux ( i, 0, 1) = zonebdy_conserved ( i, 0, 2)
      zonebdy_flux ( i, 0, 2) =
     1     zonebdy_conserved ( i, 0, 2) * vlxremi ( i)
     1   + prsremi ( i)
      zonebdy_flux ( i, 0, 3) =
     1     zonebdy_conserved ( i, 0, 3) * vlxremi ( i)
      zonebdy_flux ( i, 0, 4) =
     1     zonebdy_conserved ( i, 0, 4) * vlxremi ( i)


      zonebdy_flux ( i, 1, 1) = zonebdy_conserved ( i, 1, 2)
      zonebdy_flux ( i, 1, 2) =
     1     zonebdy_conserved ( i, 1, 2) * vlxrepl ( i)
     1   + prsrepl ( i)
      zonebdy_flux ( i, 1, 3) =
     1     zonebdy_conserved ( i, 1, 3) * vlxrepl ( i)
      zonebdy_flux ( i, 1, 4) =
     1     zonebdy_conserved ( i, 1, 4) * vlxrepl ( i)


      END DO


* Make the fifth component differently depending on which flux
* is desired.

      IF ( iprimitive .EQ. 0) THEN

      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 5) =
     1      engremi ( i)
     1    + 0.5 * rhoremi ( i)
     1    * ( vlxremi ( i)**2 + vlyremi ( i)**2
     1                        + vlzremi ( i)**2)

      zonebdy_conserved ( i, 1, 5) =
     1      engrepl ( i)
     1    + 0.5 * rhorepl ( i)
     1    * ( vlxrepl ( i)**2 + vlyrepl ( i)**2
     1                        + vlzrepl ( i)**2)

      zonebdy_flux ( i, 0, 5) =
     1     ( zonebdy_conserved ( i, 0, 5) + prsremi ( i) )
     1   * vlxremi ( i)

      zonebdy_flux ( i, 1, 5) =
     1     ( zonebdy_conserved ( i, 1, 5) + prsrepl ( i) )
     1   * vlxrepl ( i)

      END DO

      END IF


      IF ( iprimitive .EQ. - 1) THEN

* EOS
      DO i = nriembegin, nriemend

      zonebdy_conserved ( i, 0, 5) = prsremi ( i)
     1        / rhoremi ( i)**( gamma_soundspeed_remi ( i) - 1.0)

      zonebdy_conserved ( i, 1, 5) = prsrepl ( i)
     1        / rhorepl ( i)**( gamma_soundspeed_repl ( i) - 1.0)

      zonebdy_flux ( i, 0, 5) = zonebdy_conserved ( i, 0, 5)
     1                        * vlxremi ( i)

      zonebdy_flux ( i, 1, 5) = zonebdy_conserved ( i, 1, 5)
     1                        * vlxrepl ( i)

      END DO

      END IF


* If multiple species are present then initialize their fluxes etc.
* appropriately also.

      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        zonebdy_conserved ( i, 0, k) = rhoremi ( i)
     1                               * specieremi ( ifluid, i)

        zonebdy_conserved ( i, 1, k) = rhorepl ( i)
     1                               * specierepl ( ifluid, i)

        zonebdy_flux ( i, 0, k) = zonebdy_conserved ( i, 0, k)
     1                          * vlxremi ( i)

        zonebdy_flux ( i, 1, k) = zonebdy_conserved ( i, 1, k)
     1                          * vlxrepl ( i)

        conserved_re ( i, k) = 0.0

        flux_re ( i, k + 1) = 0.0

      END DO

      END DO

******************************************************************
******************************************************************

! Zero out "flux_re" here.

      DO i = nriembegin, nriemend

      conserved_re ( i, 1) = 0.0
      conserved_re ( i, 2) = 0.0
      conserved_re ( i, 3) = 0.0
      conserved_re ( i, 4) = 0.0
      conserved_re ( i, 5) = 0.0

      flux_re ( i, 1) = 0.0
      flux_re ( i, 2) = 0.0
      flux_re ( i, 3) = 0.0
      flux_re ( i, 4) = 0.0
      flux_re ( i, 5) = 0.0
      flux_re ( i, 6) = 0.0

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        conserved_re ( i, k) = 0.0

        flux_re ( i, k + 1) = 0.0

      END DO

      END DO


******************************************************************
******************************************************************

* Here we use Roe averaging and make only one call to the
* eigenvector subroutine. The "eng, temperature" etc. that we make here
* can even serve as good initial guesses for an exact EOS.


      DO i = nriembegin, nriemend

! ----------

      gamma_eos ( i) = 0.5 * ( gamma_eos_remi ( i)
     1                       + gamma_eos_repl ( i))
      gamma_soundspeed ( i) = 0.5 * ( gamma_soundspeed_remi ( i)
     1                              + gamma_soundspeed_repl ( i))
      mean_mol_wt ( i) = 0.5 * ( mean_mol_wt_remi ( i)
     1                         + mean_mol_wt_repl ( i))

! ----------

      temparr_01 ( i) = SQRT ( rhorepl ( i) / rhoremi ( i))
      temparr_02 ( i) = 1.0 / ( 1.0 + temparr_01 ( i))

      rhobar ( i) = SQRT ( rhoremi ( i) * rhorepl ( i))

      vlxbar ( i) = ( vlxremi ( i) + vlxrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)
      vlybar ( i) = ( vlyremi ( i) + vlyrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)
      vlzbar ( i) = ( vlzremi ( i) + vlzrepl ( i) * temparr_01 ( i))
     1            * temparr_02 ( i)

      total_enthalpy_bar ( i) = ( total_enthalpy_remi ( i)
     1  + total_enthalpy_repl ( i) * temparr_01 ( i)) * temparr_02 ( i)

      prsbar ( i) = total_enthalpy_bar ( i) - 0.5 * ( vlxbar ( i)**2
     1                     + vlybar ( i)**2 + vlzbar ( i)**2)
      prsbar ( i) = ( gamma_eos ( i) - 1.0) * rhobar ( i) * prsbar ( i)
     1            / gamma_eos ( i)
      prsbar ( i) = AMAX1 ( prsbar ( i), prsfloor)

! ----------

! Give yourself the option of using arithmetic averages if needed.

      IF ( ( rhobar ( i) .LE. rhofloor_plus_delta) .OR.
     1     ( prsbar ( i) .LE. prsfloor_plus_delta) ) THEN

      rhobar ( i) = 0.5 * ( rhoremi ( i) + rhorepl ( i))
      prsbar ( i) = 0.5 * ( prsremi ( i) + prsrepl ( i))

      vlxbar ( i) = 0.5 * ( vlxremi ( i) + vlxrepl ( i))
      vlybar ( i) = 0.5 * ( vlyremi ( i) + vlyrepl ( i))
      vlzbar ( i) = 0.5 * ( vlzremi ( i) + vlzrepl ( i))

      END IF

! ----------

      eng ( i) = prsbar ( i) / ( gamma_eos ( i) - 1.0)
      temperature ( i) = prsbar ( i) * mean_mol_wt ( i)
     1                 / ( gasconst * rhobar ( i))

! ----------

      END DO

! ------------------------------

      DO ifluid = 1, NFLUID

      DO i = nriembegin, nriemend

        speciebar ( ifluid, i) = 0.5 * ( specieremi ( ifluid, i)
     1                                 + specierepl ( ifluid, i))

      END DO

      END DO


* Using exact EOS's here. If that is not needed, then "eng, temperature"
* etc. have already been made above and the present call is not needed.

      IF ( use_approx_eos .NE. 1) THEN

      WRITE ( 6, *)
     1      "Stopping because EOS routines not present in demo. code."

      END IF

*********************************************************************
*********************************************************************

* Make "eigenval_re" by using the Roe-averaged values:

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed ( i) * prsbar ( i) / rhobar ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_re ( i, 1) = vlxbar ( i) - c_s ( i)
      eigenval_re ( i, 2) = vlxbar ( i)
      eigenval_re ( i, 3) = vlxbar ( i)
      eigenval_re ( i, 4) = vlxbar ( i)
      eigenval_re ( i, 5) = vlxbar ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_re ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxbar ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_rt"

*********************************************************************
*********************************************************************

      IF ( debug_stop .EQ. 1) THEN
      RETURN
      END IF


* And store the eigenvalues on either side of
* the zone boundary.

      DO i = nriembegin, nriemend

      zonebdy_eigenval ( i, 1) = eigenval_re ( i, 1)
      eigenval_db ( i, 0, 1) = eigenval_re ( i, 1)
      eigenval_db ( i, 1, 1) = eigenval_re ( i, 1)

      zonebdy_eigenval ( i, 2) = eigenval_re ( i, 2)
      eigenval_db ( i, 0, 2) = eigenval_re ( i, 2)
      eigenval_db ( i, 1, 2) = eigenval_re ( i, 2)

      zonebdy_eigenval ( i, 3) = eigenval_re ( i, 3)
      eigenval_db ( i, 0, 3) = eigenval_re ( i, 3)
      eigenval_db ( i, 1, 3) = eigenval_re ( i, 3)

      zonebdy_eigenval ( i, 4) = eigenval_re ( i, 4)
      eigenval_db ( i, 0, 4) = eigenval_re ( i, 4)
      eigenval_db ( i, 1, 4) = eigenval_re ( i, 4)

      zonebdy_eigenval ( i, 5) = eigenval_re ( i, 5)
      eigenval_db ( i, 0, 5) = eigenval_re ( i, 5)
      eigenval_db ( i, 1, 5) = eigenval_re ( i, 5)

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        zonebdy_eigenval ( i, k) = eigenval_re ( i, k)
        eigenval_db ( i, 0, k) = eigenval_re ( i, k)
        eigenval_db ( i, 1, k) = eigenval_re ( i, k)

      END DO

      END DO

*********************************************************************
*********************************************************************

* Here we explicitly need to get the eigenvalues on either
* side of the zone boundary so that we can maximize / minimize w.r.t.
* the eigenvalues of the mean state. We reset "eigenval_db" here
* so that "eigenval_db ( 0, ifield, i)" is the left-most eigenvalue
* in that characteristic field and "eigenval_db ( 1, ifield, i)" is
* the right-most.

* Make "eigenval_lf":

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed_remi ( i) * prsremi ( i)
     1         / rhoremi ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_lf ( i, 1) = vlxremi ( i) - c_s ( i)
      eigenval_lf ( i, 2) = vlxremi ( i)
      eigenval_lf ( i, 3) = vlxremi ( i)
      eigenval_lf ( i, 4) = vlxremi ( i)
      eigenval_lf ( i, 5) = vlxremi ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_lf ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxremi ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_lf"

*********************************************************************
*********************************************************************

      DO i = nriembegin, nriemend


      eigenval_db ( i, 0, 1) =
     1      AMIN1 ( eigenval_db ( i, 0, 1), eigenval_lf ( i, 1) )

      eigenval_db ( i, 1, 1) =
     1      AMAX1 ( eigenval_db ( i, 1, 1), eigenval_lf ( i, 1) )


      eigenval_db ( i, 0, 2) =
     1      AMIN1 ( eigenval_db ( i, 0, 2), eigenval_lf ( i, 2) )

      eigenval_db ( i, 1, 2) =
     1      AMAX1 ( eigenval_db ( i, 1, 2), eigenval_lf ( i, 2) )


      eigenval_db ( i, 0, 3) =
     1      AMIN1 ( eigenval_db ( i, 0, 3), eigenval_lf ( i, 3) )

      eigenval_db ( i, 1, 3) =
     1      AMAX1 ( eigenval_db ( i, 1, 3), eigenval_lf ( i, 3) )


      eigenval_db ( i, 0, 4) =
     1      AMIN1 ( eigenval_db ( i, 0, 4), eigenval_lf ( i, 4) )

      eigenval_db ( i, 1, 4) =
     1      AMAX1 ( eigenval_db ( i, 1, 4), eigenval_lf ( i, 4) )


      eigenval_db ( i, 0, 5) =
     1      AMIN1 ( eigenval_db ( i, 0, 5), eigenval_lf ( i, 5) )

      eigenval_db ( i, 1, 5) =
     1      AMAX1 ( eigenval_db ( i, 1, 5), eigenval_lf ( i, 5) )


      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        eigenval_db ( i, 0, k) =
     1        AMIN1 ( eigenval_db ( i, 0, k), eigenval_lf ( i, k) )

        eigenval_db ( i, 1, k) =
     1        AMAX1 ( eigenval_db ( i, 1, k), eigenval_lf ( i, k) )

      END DO

      END DO

*********************************************************************
*********************************************************************

* Make "eigenval_rt":

      DO i = nriembegin, nriemend

*****************

      c_s ( i) = gamma_soundspeed_repl ( i) * prsrepl ( i)
     1         / rhorepl ( i)
      c_s ( i) = SQRT ( c_s ( i) )

      eigenval_rt ( i, 1) = vlxrepl ( i) - c_s ( i)
      eigenval_rt ( i, 2) = vlxrepl ( i)
      eigenval_rt ( i, 3) = vlxrepl ( i)
      eigenval_rt ( i, 4) = vlxrepl ( i)
      eigenval_rt ( i, 5) = vlxrepl ( i) + c_s ( i)

*******************

      END DO

*******************

      IF ( NFLUID .GE. 1) THEN

      DO i = nriembegin, nriemend

      eigenval_rt ( i, NFIELD + 1: NFIELD + NFLUID) =
     1            vlxrepl ( i)

      END DO

      END IF
! End of "NFLUID .GE. 1".

* End of making "eigenval_rt"

*********************************************************************
*********************************************************************

      DO i = nriembegin, nriemend


      eigenval_db ( i, 0, 1) =
     1      AMIN1 ( eigenval_db ( i, 0, 1), eigenval_rt ( i, 1) )

      eigenval_db ( i, 1, 1) =
     1      AMAX1 ( eigenval_db ( i, 1, 1), eigenval_rt ( i, 1) )


      eigenval_db ( i, 0, 2) =
     1      AMIN1 ( eigenval_db ( i, 0, 2), eigenval_rt ( i, 2) )

      eigenval_db ( i, 1, 2) =
     1      AMAX1 ( eigenval_db ( i, 1, 2), eigenval_rt ( i, 2) )


      eigenval_db ( i, 0, 3) =
     1      AMIN1 ( eigenval_db ( i, 0, 3), eigenval_rt ( i, 3) )

      eigenval_db ( i, 1, 3) =
     1      AMAX1 ( eigenval_db ( i, 1, 3), eigenval_rt ( i, 3) )


      eigenval_db ( i, 0, 4) =
     1      AMIN1 ( eigenval_db ( i, 0, 4), eigenval_rt ( i, 4) )

      eigenval_db ( i, 1, 4) =
     1      AMAX1 ( eigenval_db ( i, 1, 4), eigenval_rt ( i, 4) )


      eigenval_db ( i, 0, 5) =
     1      AMIN1 ( eigenval_db ( i, 0, 5), eigenval_rt ( i, 5) )

      eigenval_db ( i, 1, 5) =
     1      AMAX1 ( eigenval_db ( i, 1, 5), eigenval_rt ( i, 5) )

      END DO


      DO ifluid = 1, NFLUID

      k = NFIELD + ifluid

      DO i = nriembegin, nriemend

        eigenval_db ( i, 0, k) =
     1        AMIN1 ( eigenval_db ( i, 0, k), eigenval_rt ( i, k) )

        eigenval_db ( i, 1, k) =
     1        AMAX1 ( eigenval_db ( i, 1, k), eigenval_rt ( i, k) )

      END DO

      END DO

* End of making the eigenvalues and appropriate wave speeds.

************************************************************************
************************************************************************

! First gather zones that need Einfeldt fix. 
! The bounding eigenvalues are those produced from the
! Roe solver.
! Since the Roe solver always produces a good
! value for the Alfven speed we include it in here. 

      DO i = nriembegin, nriemend

! ------------------------------

      wavespd_ein_l ( i) = AMIN1 ( eigenval_db ( i, 0, 1),
     1                             eigenval_db ( i, 0, 2) )

      wavespd_ein_r ( i) = AMAX1 ( eigenval_db ( i, 1, 5),
     1                             eigenval_db ( i, 1, 4) )

! ------------------------------

! The lines below make the flux into an LLF flux (which is more stabilizing).
! Comment them out to get an Einfeldt flux.

!      temparr_01 ( i) = AMAX1 ( ABS ( eigenval_db ( i, 0, 1)), 
!     1                          ABS ( eigenval_db ( i, 0, 2)), 
!     1                          ABS ( eigenval_db ( i, 1, 5)),
!     1                          ABS ( eigenval_db ( i, 1, 4)))
!      wavespd_ein_l ( i) = - 1.0 * temparr_01 ( i)
!      wavespd_ein_r ( i) = 1.0 * temparr_01 ( i)

! ------------------------------

! Here we make the Einfeldt fix directly without invoking a subroutine call.

      IF ( wavespd_ein_l ( i) .GT. 0.0) THEN
! Pick left conserved state.

! ------------------------------

      conserved_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_conserved ( i, 0, 1: NFIELD + NFLUID)

      flux_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_flux ( i, 0, 1: NFIELD + NFLUID)


      rho0 ( i) = rhoremi ( i)
      prs0 ( i) = prsremi ( i)

      vlx0 ( i) = vlxremi ( i)
      vly0 ( i) = vlyremi ( i)
      vlz0 ( i) = vlzremi ( i)

      specie0 ( i, 1: NFLUID) = specieremi ( i, 1: NFLUID)

! ------------------------------

      ELSE IF ( wavespd_ein_r ( i) .LT. 0.0) THEN
! Pick right conserved state.

! ------------------------------

      conserved_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_conserved ( i, 1, 1: NFIELD + NFLUID)

      flux_re ( i, 1: NFIELD + NFLUID) =
     1      zonebdy_flux ( i, 1, 1: NFIELD + NFLUID)


      rho0 ( i) = rhorepl ( i)
      prs0 ( i) = prsrepl ( i)

      vlx0 ( i) = vlxrepl ( i)
      vly0 ( i) = vlyrepl ( i)
      vlz0 ( i) = vlzrepl ( i)

      specie0 ( i, 1: NFLUID) = specierepl ( i, 1: NFLUID)

! ------------------------------

      ELSE
! Pick HLLD with two mean, conserved states.

! ------------------------------

      temparr_01 ( i) = 1.0 / ( wavespd_ein_r ( i)
     1                        - wavespd_ein_l ( i))

      conserved_re ( i, 1: 2) =

     1  ( wavespd_ein_r ( i) * zonebdy_conserved ( i, 1, 1: 2)
     1  - wavespd_ein_l ( i) * zonebdy_conserved ( i, 0, 1: 2) )
     1  * temparr_01 ( i)

     1  - ( zonebdy_flux ( i, 1, 1: 2) - zonebdy_flux ( i, 0, 1: 2) )
     1  * temparr_01 ( i)

! Make x-velocity using Eqn. (6.60) of book.
      vlx0 ( i) = conserved_re ( i, 2)
     1          / AMAX1 ( conserved_re ( i, 1), rhofloor)

! Make pressure using Eqn. (6.66) of book.
      prs0 ( i) = prsremi ( i) + rhoremi ( i)
     1          * ( wavespd_ein_l ( i) - vlxremi ( i))
     1          * ( vlx0 ( i) - vlxremi ( i))

        IF ( vlx0 ( i) .GE. 0.0) THEN

! Make density using Eqn. (6.65) of book.
        rho0 ( i) = rhoremi ( i)
     1            * ( wavespd_ein_l ( i) - vlxremi ( i))
     1            / ( wavespd_ein_l ( i) - vlx0 ( i))
        rho0 ( i) = AMAX1 ( rho0 ( i), rhofloor)

! Make transverse velocities using Eqn. (6.67) of book.
        vly0 ( i) = vlyremi ( i)
        vlz0 ( i) = vlzremi ( i)

! Make energy density using Eqn. (6.68) of book.
        eng0 ( i) = ( ( wavespd_ein_l ( i) - vlxremi ( i))
     1              * zonebdy_conserved ( i, 0, 5)
     1        - prsremi ( i) * vlxremi ( i) + prs0 ( i) * vlx0 ( i) )
     1            / ( wavespd_ein_l ( i) - vlx0 ( i))

        specie0 ( i, 1: NFLUID) = specieremi ( i, 1: NFLUID)

! Now make conserved variables and fluxes.

        conserved_re ( i, 1) = rho0 ( i)

        conserved_re ( i, 2) = rho0 ( i) * vlx0 ( i)
        conserved_re ( i, 3) = rho0 ( i) * vly0 ( i)
        conserved_re ( i, 4) = rho0 ( i) * vlz0 ( i)

        conserved_re ( i, 5) = eng0 ( i)

        conserved_re ( i, NFIELD + 1: NFIELD + NFLUID) =
     1                         rho0 ( i) * specie0 ( i, 1: NFLUID)

        flux_re ( i, 1) = conserved_re ( i, 2)

        flux_re ( i, 2) = conserved_re ( i, 2) * vlx0 ( i) + prs0 ( i)
        flux_re ( i, 3) = conserved_re ( i, 2) * vly0 ( i)
        flux_re ( i, 4) = conserved_re ( i, 2) * vlz0 ( i)

        flux_re ( i, 5) = ( eng0 ( i) + prs0 ( i)) * vlx0 ( i)

        flux_re ( i, NFIELD + 1: NFIELD + NFLUID) =
     1  conserved_re ( i, NFIELD + 1: NFIELD + NFLUID) * vlx0 ( i)

        ELSE

! Make density using Eqn. (6.65) of book.
        rho0 ( i) = rhorepl ( i)
     1            * ( wavespd_ein_r ( i) - vlxrepl ( i))
     1            / ( wavespd_ein_r ( i) - vlx0 ( i))
        rho0 ( i) = AMAX1 ( rho0 ( i), rhofloor)

! Make transverse velocities using Eqn. (6.67) of book.
        vly0 ( i) = vlyrepl ( i)
        vlz0 ( i) = vlzrepl ( i)

! Make energy density using Eqn. (6.68) of book.
        eng0 ( i) = ( ( wavespd_ein_r ( i) - vlxrepl ( i))
     1              * zonebdy_conserved ( i, 1, 5)
     1        - prsrepl ( i) * vlxrepl ( i) + prs0 ( i) * vlx0 ( i) )
     1            / ( wavespd_ein_r ( i) - vlx0 ( i))

        specie0 ( i, 1: NFLUID) = specierepl ( i, 1: NFLUID)

! Now make conserved variables and fluxes.

        conserved_re ( i, 1) = rho0 ( i)

        conserved_re ( i, 2) = rho0 ( i) * vlx0 ( i)
        conserved_re ( i, 3) = rho0 ( i) * vly0 ( i)
        conserved_re ( i, 4) = rho0 ( i) * vlz0 ( i)

        conserved_re ( i, 5) = eng0 ( i)

        conserved_re ( i, NFIELD + 1: NFIELD + NFLUID) =
     1                         rho0 ( i) * specie0 ( i, 1: NFLUID)

        flux_re ( i, 1) = conserved_re ( i, 2)

        flux_re ( i, 2) = conserved_re ( i, 2) * vlx0 ( i) + prs0 ( i)
        flux_re ( i, 3) = conserved_re ( i, 2) * vly0 ( i)
        flux_re ( i, 4) = conserved_re ( i, 2) * vlz0 ( i)

        flux_re ( i, 5) = ( eng0 ( i) + prs0 ( i)) * vlx0 ( i)

        flux_re ( i, NFIELD + 1: NFIELD + NFLUID) =
     1  conserved_re ( i, NFIELD + 1: NFIELD + NFLUID) * vlx0 ( i)

        END IF

! ------------------------------

      END IF

! End of making HLLC Riemann solver.

! ------------------------------

      END DO

! End of making all HLLC Riemann solvers.

******************************************************************
******************************************************************

! Wind-up operations done past this point:

******************************************************************

! Move the species fluxes by one to make way for the split
! in the advected and pressure parts of the x-momentum flux.

! To do that, we realize first that the 6th component of the flux will
! have to contain the pressure. Right now, the fluxes are all contiguous
! even when multiple species are present. So we first move the fluxes
! over by one when multiple species are present.

      DO ifluid = NFLUID, 1, -1

      k = NFIELD + ifluid

        DO i = nriembegin, nriemend

        flux_re ( i, k + 1) = flux_re ( i, k)

        END DO

      END DO

! Zero out the sixth component of "flux_re" because we have not
! yet made the pressure split.

      DO i = nriembegin, nriemend

      flux_re ( i, 6) = 0.0

      END DO

******************************************************************
******************************************************************

! Form the flow variables at the resolved state.
! Also split the momentum flux_re into advected and pressure parts and
! store the pressure part separately.

      IF ( split_prsflux_re .EQ. 1) THEN

      DO i = nriembegin, nriemend

! Make the part of the momentum flux that acts like a pressure.

      flux_re ( i, 6) = prs0 ( i)

      flux_re ( i, 2) = flux_re ( i, 2) - flux_re ( i, 6)

      END DO

      END IF
! End of "split_prsflux_re .EQ. 1" part.

**********************************************************************
**********************************************************************

      RETURN

      END SUBROUTINE EULER_RIEM_HLLC

**********************************************************************
**********************************************************************


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EULER_RIEM_HLLEM_PTWISE

     1     ( NFIELD_EULER, NFLUID_EULER,
     1       gamma, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vlx_L, vlx_R,
     1        vly_L, vly_R,
     1        vlz_L, vlz_R,
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

! Builds a 1D HLLEM Riemann solver for Euler based on the work of
! Dumbser & Balsara (2016) JCP and Balsara et al. (2016) JCP.

! Exact preservation of stationary intermediate waves is achieved.

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
! flux_re : Five component flux, with pressure split, specific to certain
!            my MHD code.
! entropy_flux_re : Flux for entropy variable.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NFIELD_EULER, NFLUID_EULER,
     1        use_approx_eos, split_prsflux_re, iHLLEM_choice

      REAL gamma, smallnum, prsfloor, rhofloor, pi, gasconst,

     1     cp_specie ( 0 : NFLUID_EULER),
     1     molewt_specie ( 0 : NFLUID_EULER),
     1     heat_of_formation ( 0 : NFLUID_EULER),

     1     rho_L, rho_R, prs_L, prs_R,
     1     vlx_L, vlx_R, vly_L, vly_R, vlz_L, vlz_R,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     blend_HLL, S_L, S_R, entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID_EULER) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD_EULER + NFLUID_EULER) :: 
     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

      REAL, DIMENSION ( NFIELD_EULER + 1 + NFLUID_EULER) :: flux_re

! ----------------------------------------------------------------------

! Variables past this point are not used in the subroutine interface.

      INTEGER i, j, iprim, idirection, debug_stop

      REAL rhofloor_plus_delta, prsfloor_plus_delta,

     1     eng_L, eng_R, v2_L, v2_R,
     1     total_enthalpy_L, total_enthalpy_R,

     1     entropy_L, entropy_R, entropy_flux_L, entropy_flux_R,

     1     gamma_eos, gamma_soundspeed, mean_mol_wt,
     1     temp_01, temp_02, temp_03,
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, total_enthalpybar,
     1     eng, temperature,

     1     c_s

      REAL, DIMENSION ( 0: NFLUID_EULER) :: speciebar

      REAL, DIMENSION ( NFIELD_EULER + NFLUID_EULER) :: U_L, U_R,
     1     Del_U, eigenval, delta, eigenwt

      REAL, DIMENSION ( NFIELD_EULER + NFLUID_EULER,
     1                  NFIELD_EULER + NFLUID_EULER) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

      IF ( NFIELD_EULER .NE. 5) THEN
      WRITE ( 6, *)"EULER_RIEM_HLLEM_PTWISE : NFIELD","_EULER .NE. 5",
     1              NFIELD_EULER, 5
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

      total_enthalpy_L = ( eng_L + prs_L
     1    + 0.5 * rho_L * v2_L ) / rho_L

! -----

      v2_R = vlx_R**2 + vly_R**2 + vlz_R**2

      total_enthalpy_R = ( eng_R + prs_R
     1    + 0.5 * rho_R * v2_R ) / rho_R

! ------------------------------
! Make "U_L, U_R".

      U_L ( 1) = rho_L
      U_L ( 2) = rho_L * vlx_L
      U_L ( 3) = rho_L * vly_L
      U_L ( 4) = rho_L * vlz_L
      U_L ( 5) = eng_L + 0.5 * rho_L * v2_L
#if ( NFLUID_EULER >= 1)
      U_L ( 6: NFIELD_EULER + 1 + NFLUID_EULER) =
     1     rho_L * specie_L ( 1:  NFLUID_EULER)
#endif

! ----------

      U_R ( 1) = rho_R
      U_R ( 2) = rho_R * vlx_R
      U_R ( 3) = rho_R * vly_R
      U_R ( 4) = rho_R * vlz_R
      U_R ( 5) = eng_R + 0.5 * rho_R * v2_R
#if ( NFLUID_EULER >= 1)
      U_R ( 6: NFIELD_EULER + 1 + NFLUID_EULER) =
     1     rho_R * specie_R ( 1:  NFLUID_EULER)
#endif

! ------------------------------
! Make "F_L, F_R".

      F_L ( :) = vlx_L * U_L ( :)


      F_L ( 2) = F_L ( 2) + prs_L

      F_L ( 3) = F_L ( 3)

      F_L ( 4) = F_L ( 4)

      F_L ( 5) = F_L ( 5) + vlx_L * prs_L

! ----------

      F_R ( :) = vlx_R * U_R ( :)


      F_R ( 2) = F_R ( 2) + prs_R

      F_R ( 3) = F_R ( 3)

      F_R ( 4) = F_R ( 4)

      F_R ( 5) = F_R ( 5) + vlx_R * prs_R

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

      prsbar = total_enthalpybar
     1       - 0.5 * ( vlxbar**2 + vlybar**2 + vlzbar**2)
      prsbar = ( gamma_eos - 1.0) * rhobar * prsbar / gamma_eos
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

      END IF

      eng = prsbar / ( gamma_eos - 1.0)
      temperature = prsbar * mean_mol_wt / ( gasconst * rhobar)

! ----------

#if ( NFLUID_EULER >= 1)
      temp_01 = SQRT ( rho_R / rho_L)
      temp_02 = 1.0 / ( 1.0 + temp_01)

      speciebar ( 1: NFLUID_EULER) = temp_02
     1                     * ( specie_L ( 1: NFLUID_EULER)
     1                       + specie_R ( 1: NFLUID_EULER) * temp_01)
#endif

! ----------------------------------------------------------------------
! Make extremal speeds, "S_L, S_R", using the Roe averaged state.

      c_s = gamma_soundspeed * prsbar / rhobar
      c_s = SQRT ( c_s )

! ----------

      S_L = vlxbar - c_s
      S_R = vlxbar + c_s

! ----------------------------------------------------------------------
! Expand the extremal speeds, "S_L, S_R", using "rho_L" etc.

      c_s = gamma_soundspeed_L * prs_L / rho_L
      c_s = SQRT ( c_s )

! ----------

      S_L = AMIN1 ( S_L, vlx_L - c_s)
! While it is ok to uncomment the next line, in principle, the left state
! does not contribute to "S_R".
!      S_R = AMAX1 ( S_R, vlx_L + c_s)

! ----------------------------------------------------------------------
! Expand the extremal speeds, "S_L, S_R", using "rho_R" etc.

      c_s = gamma_soundspeed_R * prs_R / rho_R
      c_s = SQRT ( c_s )

! ----------

! While it is ok to uncomment the next line, in principle, the right state
! does not contribute to "S_L".
!      S_L = AMIN1 ( S_L, vlx_R - c_s)
      S_R = AMAX1 ( S_R, vlx_R + c_s)

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
! First call "EIGENVEC_EULER_PTWISE_LINDEG". Load up "eigenval,
! lf_eigenvec, rt_eigenvec".

      iprim = 0  ! Need eigenvectors in conserved variables
      idirection = 1  ! in the x-direction.
      debug_stop = 0

      CALL EIGENVEC_EULER_PTWISE_LINDEG

     1             ( iprim, idirection, NFIELD_EULER, NFLUID_EULER,

     1               gamma, smallnum, prsfloor, rhofloor,
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
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

      DO j = 1, NFIELD_EULER + NFLUID_EULER

! ----------

! Just include intermediate waves for contact and shear discontinuities.
        IF ( ( j >= 2) .AND. ( j <= 4) ) THEN

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

        DO i = 1, NFIELD_EULER + NFLUID_EULER

        eigenwt ( j) = eigenwt ( j)
     1               + lf_eigenvec ( j, i) * Del_U ( i)

        END DO  ! End of "DO i = 1, NFIELD_EULER + NFLUID_EULER" loop.

! ----------

        END IF  ! End of "IF ( ( j >= 2) .AND. ( j <= 4) )".

      END DO  ! End of "DO j = 1, NFIELD_EULER + NFLUID_EULER" loop.

! ------------------------------------------------------------
! Lastly, build "dU_HLLEM, dF_HLLEM". If the linear profile from "dU_HLLEM"
! gives physical states, add in the contributions to the fluxes.

      dU_HLLEM ( :) = 0.0
      dU_HLLEM_Min ( :) = 0.0
      dF_HLLEM ( :) = 0.0

! ------------------------------

      DO j = 1, NFIELD_EULER + NFLUID_EULER

! ----------

! Just include intermediate waves for contact and shear discontinuities.
        IF ( ( j >= 2) .AND. ( j <= 4) ) THEN

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

        END IF  ! End of "IF ( ( j >= 2) .AND. ( j <= 4) )".

      END DO  ! End of "DO j = 1, NFIELD_EULER + NFLUID_EULER" loop.

! ------------------------------------------------------------

! Density and pressure from HLL state should be bounded by densities and
! pressures from either side.

      rhobar = AMAX1 ( U_HLL ( 1), rhofloor)
      temp_01 = 1.0 / rhobar

      vlxbar = U_HLL ( 2) * temp_01
      vlybar = U_HLL ( 3) * temp_01
      vlzbar = U_HLL ( 4) * temp_01

      prsbar = U_HLL ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)

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

      prsbar = Del_U ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)

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

      prsbar = Del_U ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)

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

      flux_re ( 1: NFIELD_EULER) = F_HLLEM ( 1: NFIELD_EULER)

! The pressure that is split off from the momentum flux may (optionally)
! be stored in the sixth component of the flux that is returned.
      flux_re ( 6) = 0.0

      DO j = 6, NFIELD_EULER + NFLUID_EULER
      flux_re ( j + 1) = F_HLLEM ( j)
      END DO

! ----------
! Also split the momentum flux_re into advected and pressure parts and
! store the pressure part separately. Do this only if it is called for.

      IF ( split_prsflux_re .EQ. 1) THEN

      rhobar = AMAX1 ( U_HLLEM ( 1), rhofloor)
      temp_01 = 1.0 / rhobar

      vlxbar = U_HLLEM ( 2) * temp_01
      vlybar = U_HLLEM ( 3) * temp_01
      vlzbar = U_HLLEM ( 4) * temp_01

      prsbar = U_HLLEM ( 5)
     1       - 0.5 * rhobar * ( vlxbar**2 + vlybar**2 + vlzbar**2)

      prsbar = AMAX1 ( prsbar * ( gamma_eos - 1.0), prsfloor)

! This is the part of the x-momentum flux that acts like a pressure.
      flux_re ( 6) = flux_re ( 2) - rhobar * vlxbar**2

      flux_re ( 2) = flux_re ( 2) - flux_re ( 6)

      END IF

! ----------------------------------------------------------------------

      END SUBROUTINE EULER_RIEM_HLLEM_PTWISE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EULER_RIEM_HLLEM_DB

     1   ( nriembegin, nriemend,
     1     ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,

     1     gamma, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Packages the pointwise HLLEM RS in a form that works with the
! dimension-by-dimension Riemann solver's interface. This is just
! useful for testing.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: nriembegin, nriemend,
     1        ONEDFLUIDBEGIN, ONEDFLUIDEND, NFIELD, NFLUID,
     1        ientropyfix, iprimitive, use_approx_eos,
     1        split_prsflux_re

      INTEGER :: debug_stop

      REAL :: gamma, smallnum, prsfloor, rhofloor,
     1     fracdiff, gasconst, divratio_einfeldt_rs

      REAL, DIMENSION ( 0 : NFLUID) ::
     1      cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr

      REAL, DIMENSION
     1   ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD + NFLUID, NFIELD + NFLUID) ::
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
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION
     1    ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Past this point we only declare variables that are internal to the
! subroutine.

      INTEGER i, iHLLEM_choice

      REAL rho_L, rho_R, prs_L, prs_R, pi,
     1     vlx_L, vlx_R, vly_L, vly_R, vlz_L, vlz_R,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     blend_HLL, S_L, S_R, entropy_flux_pt

      REAL, DIMENSION ( 0: NFLUID) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD + NFLUID) ::

     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: flux_pt

! ----------------------------------------------------------------------

      pi = 4.0 * ATAN ( 1.0)

      blend_HLL = 0.0

      iHLLEM_choice = 1

      DO i = nriembegin, nriemend

! ------------------------------

      rho_L = rhoremi ( i)
      prs_L = prsremi ( i)

      vlx_L = vlxremi ( i)
      vly_L = vlyremi ( i)
      vlz_L = vlzremi ( i)

      mean_mol_wt_L = mean_mol_wt_remi ( i)
      gamma_eos_L = gamma_eos_remi ( i)
      gamma_soundspeed_L = gamma_soundspeed_remi ( i)


      rho_R = rhorepl ( i)
      prs_R = prsrepl ( i)

      vlx_R = vlxrepl ( i)
      vly_R = vlyrepl ( i)
      vlz_R = vlzrepl ( i)

      mean_mol_wt_R = mean_mol_wt_repl ( i)
      gamma_eos_R = gamma_eos_repl ( i)
      gamma_soundspeed_R = gamma_soundspeed_repl ( i)

! ------------------------------

      CALL EULER_RIEM_HLLEM_PTWISE

     1     ( NFIELD, NFLUID,
     1       gamma, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vlx_L, vlx_R,
     1        vly_L, vly_R,
     1        vlz_L, vlz_R,
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

      END SUBROUTINE EULER_RIEM_HLLEM_DB

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EIGENVEC_EULER_PTWISE_LINDEG

     1             ( iprimitive, idirection, NFIELD_EULER, NFLUID_EULER,

     1               gamma, smallnum, prsfloor, rhofloor, 
     1               debug_stop,

     1               gasconst,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
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

      INTEGER, INTENT ( IN) :: iprimitive, idirection,
     1         NFIELD_EULER, NFLUID_EULER

      INTEGER, INTENT ( INOUT) :: debug_stop

      REAL, INTENT ( IN) :: gamma, smallnum, prsfloor, rhofloor,
     1      gasconst

      REAL, INTENT ( INOUT) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar

      REAL, DIMENSION ( 0: NFLUID_EULER), INTENT ( INOUT) :: speciebar

      REAL, DIMENSION ( NFIELD_EULER + NFLUID_EULER), INTENT ( INOUT) ::
     1     eigenval

      REAL, DIMENSION ( NFIELD_EULER + NFLUID_EULER,
     1                  NFIELD_EULER + NFLUID_EULER),
     1      INTENT ( INOUT) ::
     1     lf_eigenvec, rt_eigenvec

! ----------------------------------------------------------------------

      INTEGER j, k, l

! Variables for making wave speeds and eigenvectors.

      REAL :: c_s, eps_zone,

     1     eng,
     1     gamma_eos, gamma_soundspeed,
     1     dprs_drho, dprs_deng, deng_drho, deng_dprs,

     1     temparr_01, temparr_02, temparr_03

! Temporary usage of left and right eigenvectors in the primitive variables.
! Transformation matrices to go from conserved to primitive and back.

      REAL, DIMENSION ( NFIELD_EULER + NFLUID_EULER,
     1                  NFIELD_EULER + NFLUID_EULER) ::

     1     lf_eigenvec1, rt_eigenvec1, dudw, dwdu

! ----------------------------------------------------------------------

      IF ( NFIELD_EULER .NE. 5) THEN
      WRITE ( 6, *)
     1      "EIGENVEC_EULER_PTWISE_LINDEG : NFIELD_EULER .NE. 5",
     1      NFIELD_EULER, 5
      debug_stop = 1
      STOP
      END IF

      IF ( ( iprimitive < 0) .OR. ( iprimitive > 1)) THEN
      WRITE ( 6, *)
     1      "EIGENVEC_EULER_PTWISE_LINDEG : iprimitive = ", iprimitive
      debug_stop = 1
      STOP
      END IF

      IF ( ( idirection < 0) .OR. ( idirection > 3)) THEN
      WRITE ( 6, *)
     1      "EIGENVEC_EULER_PTWISE_LINDEG : idirection = ", idirection
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

      END IF

! ----------

      IF ( idirection == 3) THEN

      temparr_01 = vlxbar
      temparr_02 = vlybar
      temparr_03 = vlzbar

      vlxbar = temparr_03
      vlybar = temparr_01
      vlzbar = temparr_02

      END IF

! ----------------------------------------------------------------------

! Bound the density and pressure variables.

      rhobar = AMAX1 ( rhobar, rhofloor)
      prsbar = AMAX1 ( prsbar, prsfloor)

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

      c_s = gamma_soundspeed * prsbar / rhobar
      c_s = SQRT ( c_s )

! ----------

      eigenval ( 1) = vlxbar - c_s
      eigenval ( 2) = vlxbar
      eigenval ( 3) = vlxbar
      eigenval ( 4) = vlxbar
      eigenval ( 5) = vlxbar + c_s

! ----------

#if ( NFLUID_EULER >= 1)
      eigenval ( NFIELD_EULER + 1: NFIELD_EULER + NFLUID_EULER) = vlxbar
#endif

! ----------------------------------------------------------------------

! First zero out the eigenvectors. This is only needed because we can
! then fill in the species eigenvectors by putting 1's at suitable
! places on the diagonals.

      rt_eigenvec1 ( :, :) = 0.0
      lf_eigenvec1 ( :, :) = 0.0

! ----------

! Make right & left eigenvectors for shear and entropy waves.

! ----------

! Shear waves

! ----------

      rt_eigenvec1 ( 3, 3) = 1.0

      rt_eigenvec1 ( 4, 4) = 1.0

! ----------

      lf_eigenvec1 ( 3, 3) = 1.0

      lf_eigenvec1 ( 4, 4) = 1.0

! ----------

! Entropy wave

! ----------

      rt_eigenvec1 ( 1, 2) = 1.0

! ----------

      lf_eigenvec1 ( 2, 1) = 1.0
      lf_eigenvec1 ( 2, 5) = - 1.0 / c_s**2

! ----------------------------------------------------------------------

! For the multiple species, fix up the eigenvectors that carry the
! contributions from the species fractions. The zero components of
! the eigenvectors have already been previously assigned.

#if ( NFLUID_EULER >= 1)

      DO j = NFIELD_EULER + 1, NFIELD_EULER + NFLUID_EULER

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


! Derivatives of U(2) w.r.t. V(1:8)
      dudw ( 2, 1) = vlxbar
      dudw ( 2, 2) = rhobar


      dudw ( 3, 1) = vlybar
      dudw ( 3, 3) = rhobar


      dudw ( 4, 1) = vlzbar
      dudw ( 4, 4) = rhobar


      dudw ( 5, 1) = 0.5 * ( vlxbar**2 + vlybar**2 + vlzbar**2)
     1             + deng_drho
      dudw ( 5, 2) = rhobar * vlxbar
      dudw ( 5, 3) = rhobar * vlybar
      dudw ( 5, 4) = rhobar * vlzbar
      dudw ( 5, 5) = deng_dprs


#if ( NFLUID_EULER >= 1)
      DO k = 6, NFIELD_EULER + NFLUID_EULER

      dudw ( k, 1) = speciebar ( k - 5)
      dudw ( k, k) = rhobar

      END DO
#endif

! ----------------------------------------------------------------------

! Make transformation matrices to the conserved variables. Then for
! the left eigenvectors. Don't make those matrix elements that are zero.

      temparr_01 = 1.0 / rhobar


      dwdu ( 1, 1) = 1.0
      dwdu ( 2, 1) = - temparr_01 * vlxbar
      dwdu ( 3, 1) = - temparr_01 * vlybar
      dwdu ( 4, 1) = - temparr_01 * vlzbar
      dwdu ( 5, 1) = dprs_drho + dprs_deng * 0.5 * ( vlxbar**2
     1                                 + vlybar**2 + vlzbar**2)


      dwdu ( 2, 2) = temparr_01
      dwdu ( 5, 2) = - dprs_deng * vlxbar


      dwdu ( 3, 3) = temparr_01
      dwdu ( 5, 3) = - dprs_deng * vlybar


      dwdu ( 4, 4) = temparr_01
      dwdu ( 5, 4) = - dprs_deng * vlzbar


      dwdu ( 5, 5) = dprs_deng


#if ( NFLUID_EULER >= 1)
      DO k = 6, NFIELD_EULER + NFLUID_EULER

      dwdu ( k, 1) = - temparr_01 * speciebar ( k - 5)
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
      DO k = 1, NFIELD_EULER + NFLUID_EULER 

      IF ( ( k == 2) .OR. ( k == 3) .OR. ( k == 4) ) THEN

      DO j = 1, NFIELD_EULER + NFLUID_EULER

      DO l = 1, NFIELD_EULER + NFLUID_EULER

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

      DO k = 1, NFIELD_EULER + NFLUID_EULER

! -----

      temparr_01 = rt_eigenvec ( 2, k)
      temparr_02 = rt_eigenvec ( 3, k)
      temparr_03 = rt_eigenvec ( 4, k)

      rt_eigenvec ( 2, k) = temparr_03
      rt_eigenvec ( 3, k) = temparr_01
      rt_eigenvec ( 4, k) = temparr_02

! -----

      temparr_01 = lf_eigenvec ( k, 2)
      temparr_02 = lf_eigenvec ( k, 3)
      temparr_03 = lf_eigenvec ( k, 4)

      lf_eigenvec ( k, 2) = temparr_03
      lf_eigenvec ( k, 3) = temparr_01
      lf_eigenvec ( k, 4) = temparr_02

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

      DO k = 1, NFIELD_EULER + NFLUID_EULER

! -----

      temparr_01 = rt_eigenvec ( 2, k)
      temparr_02 = rt_eigenvec ( 3, k)
      temparr_03 = rt_eigenvec ( 4, k)

      rt_eigenvec ( 2, k) = temparr_02
      rt_eigenvec ( 3, k) = temparr_03
      rt_eigenvec ( 4, k) = temparr_01

! -----

      temparr_01 = lf_eigenvec ( k, 2)
      temparr_02 = lf_eigenvec ( k, 3)
      temparr_03 = lf_eigenvec ( k, 4)

      lf_eigenvec ( k, 2) = temparr_02
      lf_eigenvec ( k, 3) = temparr_03
      lf_eigenvec ( k, 4) = temparr_01

! -----

      END DO

! -----

      END IF

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE EIGENVEC_EULER_PTWISE_LINDEG

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
















