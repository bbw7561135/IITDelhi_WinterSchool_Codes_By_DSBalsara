! ----------------------------------------------------------------------

#define BURGERS 1

      PROGRAM NONLIN_SCALAR_RK

! Solves the Burgers and Buckley-Leverett equations or the general scalar
! hyp eqn. u_t + f(u)_x = 0 in 1d.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix, ixp1, ixm1, islope, nzones, ntstep,
     1        istep, iproblem, ifracstep

      REAL mu, t, dt, t_break, dx, tstop, umin, umax,
     1     tempa, tempb, tempc, speed_lf, speed_rt

      PARAMETER ( nzones = 500)

      REAL, DIMENSION ( -1: nzones + 2) :: x, un, un0, dun,
     1         un_rt, un_lf, flx, tu

! ----------------------------------------------------------------------

! Declare some info. for TVD limiters.

      REAL MC_BETA_LIMITER, MIN_MOD, VAN_LEER, SUPERBEE,
     1     F, DF, DDF, a, b, mcb

      MIN_MOD ( a, b) = AMIN1 ( ABS (a), ABS (b))
     1                  * ( SIGN ( 0.5, a) + SIGN ( 0.5, b))

      MC_BETA_LIMITER ( a, b, mcb) =
     1              ( SIGN ( 0.5, a) + SIGN ( 0.5, b))
     1            * AMIN1 ( 0.5 * ( ABS ( a + b)),
     1                      mcb * ABS ( a), mcb * ABS ( b) )

      VAN_LEER ( a, b) = a * b * SIGN ( 1.0, a)
     1                    * ( SIGN ( 1.0, a * b) + 1.0)
     1                    / ( ABS (a) + ABS (b) + 1.0e-30)

      SUPERBEE ( a, b) = ( SIGN ( 0.5, a) + SIGN ( 0.5, b))
     1               * AMAX1 ( AMIN1 ( ABS ( a), 2.0 * ABS ( b)),
     1                         AMIN1 ( ABS ( b), 2.0 * ABS ( a)) )

#if ( BURGERS == 1)
! Use Burgers flux.

      F ( a) = 0.5 * a**2

      DF ( a) = a

      DDF ( a) = 1.0

#else
! Use Buckley-Leverett flux.

      F ( a) = 4.0 * a**2 / ( 4.0 * a**2 + (1.0 - a)**2 )

      DF ( a) = 8.0 * a * (1.0 - a) / ( 4.0 * a**2 + (1.0 - a)**2 )**2

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "rhoa0001")

      OPEN ( unit = 5, file = "tracechar.in")

! islope controls how the slopes are built. islope = 1,2 for min-mod & MC-beta.
      islope = 2

! iproblems 1, 2, 3, 5 for Burgers ; 4 for Buckley-Leverett
      iproblem = 3

      mu = 0.6 ! CFL #

      dx = 1.0 / REAL ( nzones)   ! Zone size

      t = 0.0

! Stopping time.
! Set to 0.0, 0.08, 0.1116, 0.6 for iproblem = 1; the shock forms at 0.1166.
! Set to 0.0, 0.3, 0.6 for iproblem = 2; the shock propagates with unit speed.
! Set to 0.0, 0.5, 1.0 for iproblem = 3; the rarefaction propagates to the right.
! Set to 0.5 for iproblem = 4; the compound wave travels to the right.
! Set to 0.35 for iproblem = 5
      tstop = 1.0

      ntstep = 35000 ! # of timesteps.

      mcb = 2.0    ! ranges between 1 and 2.

!      WRITE ( 6, *)"islope = ", islope," iproblem = ", iproblem,
!     1             " mu = ", mu," mcb = ", mcb
!      WRITE ( 6, *)"dx = ", dx," tstop = ", tstop,
!     1             " ntstep = ", ntstep

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Initialize the data. Set ghost zone values too.

      umin = 1.0e32
      umax = -1.0e32

      dt = 1.0e32

      DO ix = 1 - 2, nzones + 2

      x ( ix) = - 0.5 + ( ix - 0.5) * dx

      IF ( iproblem .EQ. 1) THEN
! Gaussian being initialized

      un ( ix) = 0.5 + EXP ( - ( x ( ix) + 0.25)**2 / 0.1**2 )

      ELSE IF ( iproblem .EQ. 2) THEN
! Shock being initialized.

        IF ( x ( ix) .LT. -0.25) THEN
        un ( ix) = 2.0
        ELSE
        un ( ix) = 0.0
        END IF

      ELSE IF ( iproblem .EQ. 3) THEN
! Rarefaction fan being initialized.

        IF ( x ( ix) .LT. -0.25) THEN
        un ( ix) = 0.0
        ELSE
        un ( ix) = 0.5
        END IF

      ELSE IF ( iproblem .EQ. 4) THEN
! Compound wave from Buckley-Leverett being initialized.

        IF ( x ( ix) .LT. -0.4) THEN
        un ( ix) = 1.0
        ELSE
        un ( ix) = 0.0
        END IF

      ELSE IF ( iproblem .EQ. 5) THEN
! Rarefaction fan being initialized with -1 and 1

        IF ( x ( ix) .LT. 0.0) THEN
        un ( ix) = -0.5
        ELSE
        un ( ix) = 1.0
        END IF

      END IF

#if ( BURGERS == 1)

      dt = AMIN1 ( dt, mu * dx / AMAX1 ( ABS ( DF ( un ( ix)) ),
     1                                   1.0e-16)  )

#else

      tempa = 0.5 * ( un ( ix) + un ( ix + 1) )
      tempa = ABS ( DF ( tempa))

      dt = AMIN1 ( dt, mu * dx / AMAX1 ( ABS ( DF ( un ( ix)) ),
     1                            tempa, 1.0e-16)  )

#endif

      umin = AMIN1 ( umin, un ( ix))
      umax = AMAX1 ( umax, un ( ix))

      END DO ! End of "DO ix ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Find the wave breaking time, i.e. the time at which the first shock forms.

      t_break = 1.0e5

      DO ix = 1, nzones

#if ( BURGERS == 1)
! For Burgers flux.

      tempa = DDF ( un ( ix))
     1      * 0.5 * ( un ( ix + 1) - un ( ix - 1) ) / dx

#else
! For Buckley-Leverett problem we use numerical differentiation.

      tempa =DF ( un ( ix + 1)) - DF ( un ( ix - 1))
      tempa = tempa / (2.0 * dx)
      tempa = tempa
     1      * 0.5 * ( un ( ix + 1) - un ( ix - 1) ) / dx

#endif

      IF ( tempa .LT. 0.0) THEN
      t_break = AMIN1 ( t_break, - 1.0 / tempa)
      END IF

      END DO ! End of "DO ix ="

! ----------

!      WRITE ( 6, *)"dt, t_break, umin, umax = ",
!     1              dt, t_break, umin, umax

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Take "ntstep" timesteps.

      DO istep = 1, ntstep

      DO ifracstep = 1, 2

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Apply the limiter.

        dun ( :) = 0.0

        DO ix = 1 - 1, nzones + 1   ! ix runs over zone centers.

        IF ( islope .EQ. 1) THEN

        dun ( ix) = MIN_MOD ( un ( ix + 1) - un ( ix),
     1                        un ( ix) - un ( ix - 1) )

        ELSE

        dun ( ix) = MC_BETA_LIMITER ( un ( ix + 1) - un ( ix),
     1                                un ( ix) - un ( ix - 1), mcb)

        END IF

        END DO ! End of "DO ix ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Build the fluxes using the HLL Riemann solver.

      un_rt ( :) = un ( :) + 0.5 * dun ( :)
      un_lf ( :) = un ( :) - 0.5 * dun ( :)

      flx ( :) = 0.0

      DO ix = 1 - 1, nzones    ! ix runs over zone boundaries.

        tempa = DF ( un_rt ( ix))
        tempb = DF ( un_lf ( ix + 1))
        tempc = DF ( 0.5 * ( un_rt ( ix) + un_lf ( ix + 1)) )

        speed_lf = AMIN1 ( AMIN1 ( tempa, tempc), 0.0)
        speed_rt = AMAX1 ( AMAX1 ( tempb, tempc), 0.0)

        IF ( ABS ( speed_rt - speed_lf) .LT. 1.0e-16)
     1   speed_rt = speed_lf + 1.0e-16

        tempa = speed_rt / ( speed_rt - speed_lf)
        tempb = speed_lf / ( speed_rt - speed_lf)
        tempc = speed_rt * speed_lf / ( speed_rt - speed_lf)

        flx ( ix) = tempa * F ( un_rt ( ix))
     1            - tempb * F ( un_lf ( ix + 1))
     1            + tempc * ( un_lf ( ix + 1) - un_rt ( ix))

      END DO ! End of "DO ix ="

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Two-step Runge-Kutta update.

      IF ( ifracstep .EQ. 1) THEN
! Predictor part of RK.

      un0 ( :) = un ( :)

        DO ix = 1, nzones    ! ix runs over zone centers.

        un ( ix) = un ( ix)
     1     - ( flx ( ix) - flx ( ix - 1) ) * 0.5 * dt / dx

        END DO ! End of "DO ix ="

      ELSE
! Corrector part of RK.

        DO ix = 1, nzones    ! ix runs over zone centers.

        un ( ix) = un0 ( ix)
     1     - ( flx ( ix) - flx ( ix - 1) ) * dt / dx

        END DO ! End of "DO ix ="

      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 1) THEN
! Impose periodic boundary values on "un" here if desired.

      un ( 0) = un ( nzones)
      un ( -1) = un ( nzones - 1)

      un ( nzones + 1) = un ( 1)
      un ( nzones + 2) = un ( 2)

      ELSE
! Impose continuitive boundary values on "un" here if desired.

      un ( 0) = un ( 1)
      un ( -1) = un ( 1)

      un ( nzones + 1) = un ( nzones)
      un ( nzones + 2) = un ( nzones)

      END IF

! ----------------------------------------------------------------------

#if ( 0 == 1)
! Write out diagnostics.

        umin = 1.0e32
        umax = -1.0e32

        DO ix = 1, nzones
        umin = AMIN1 ( umin, un ( ix))
        umax = AMAX1 ( umax, un ( ix))
        END DO

!        WRITE ( 6, *)"from output : time, umin, umax = ",
!     1                t, umin, umax

#endif

! ----------------------------------------------------------------------

#if ( 1 == 1)
! Write out data for constructing the plot for characteristic fields.
      IF ( ifracstep .EQ. 1) THEN

        DO ix = 1 - 2, nzones + 2
        WRITE ( 5, *) un ( ix)
        END DO

      END IF
#endif

! ----------------------------------------------------------------------

      IF ( ifracstep .EQ. 2) t = t + dt  ! Increment time.

        IF ( ( t .GT. tstop) .AND. ( ifracstep .EQ. 2)) THEN
! Write out data if we have reached the final time and exit.

        umin = 1.0e32
        umax = -1.0e32

        DO ix = 1, nzones
        umin = AMIN1 ( umin, un ( ix))
        umax = AMAX1 ( umax, un ( ix))
        END DO

!        WRITE ( 6, *)"from output : umin, umax = ", umin, umax
!        WRITE ( 6, *)"istep = ", istep

!        WRITE ( 6, *)"-------------- Output --------------------"

#if ( 0 == 1)
        DO ix = 1, nzones
        WRITE ( 6, *) un ( ix)
        END DO

        DO ix = 1, nzones
        WRITE ( 6, *) x ( ix)
        END DO
#else
        WRITE ( 6, *) ( un ( ix), ix = 1, nzones)
#endif

        GO TO 100

        END IF

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO ! End of "DO ifracstep =".

      END DO ! End of "DO istep =".

100   CONTINUE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      STOP

      END PROGRAM NONLIN_SCALAR_RK

! ----------------------------------------------------------------------



