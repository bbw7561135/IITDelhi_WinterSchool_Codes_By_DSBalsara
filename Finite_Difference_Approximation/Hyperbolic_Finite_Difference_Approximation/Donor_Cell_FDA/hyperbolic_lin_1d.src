! ----------------------------------------------------------------------

      PROGRAM ONED_SCALAR_ADVECTION

! Solves linear advection problem in 1d. Gaussian or top-hat profiles possible.
! Uses explicit Lax-Friedrichs, Lax-Wendroff, 2nd order Runge-Kutta or Donor
! Cell schemes. 

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzones, ntstep, istep, ix, ixp1, ixp2, ixm1, ixm2,
     1        iproblem, ischeme

      PARAMETER ( nzones = 100)

      REAL x, dx, t, dt, tstop, mu, c

      REAL, DIMENSION ( nzones) :: u, rhs

! ----------------------------------------------------------------------

! Declare some scalar variables.

      OPEN ( unit = 6, file = "rhoa0001")

      mu = 0.4 ! CFL number; limited by unity.

      iproblem = 2  ! iproblem=1 ==> Gaussian; iproblem=2 ==> Square pulse

! ischeme=1,2,3,4 ==> Lax-Friedrichs, Lax-Wendroff, Heun, 1st order upwind.
      ischeme = 4

      c = 1.0  ! Advection speed, assumed positive

! Stopping time. 0.0, 1.0, 2.0 for Gaussian; 0.0, 0.25 and 0.75 for Sq. pulse.
      tstop = 1.0

      dx = 1.0 / REAL ( nzones)   ! Make zone size on a unit 1d domain.

      dt = mu * dx / c
      t = 0.0

      ntstep = 35000 ! # of timesteps.

!      WRITE ( 6, *)"mu = ", mu," c = ", c," iproblem = ", iproblem,
!     1             " ischeme = ", ischeme
!      WRITE ( 6, *)"tstop = ", tstop," dx = ", dx," dt = ", dt,
!     1             " ntstep = ", ntstep

! ----------------------------------------------------------------------

! Initialize the data with a Gaussian or square wave pulse.
! Periodic boundary conditions.

      DO ix = 1, nzones
      x = - 0.5 + ( ix - 0.5) * dx

      IF ( iproblem .EQ. 1) THEN
! Gaussian pulse set up here.

        u ( ix) = EXP ( -( x / 0.1)**2)

      ELSE
! Square pulse set up here.

        IF ( ( -0.05 .LE. x) .AND. ( x .LE. 0.05) ) THEN
        u ( ix) = 1.0
        ELSE
        u ( ix) = 0.0
        END IF

      END IF

      END DO

! ----------------------------------------------------------------------

! Take "ntstep" timesteps.

      DO istep = 1, ntstep

! ----------------------------------------------------------------------

        DO ix = 1, nzones

        ixp2 = ix + 2
        if ( ixp2 .GT. nzones) ixp2 = ixp2 - nzones
        ixp1 = ix + 1
        if ( ixp1 .GT. nzones) ixp1 = ixp1 - nzones
        ixm1 = ix - 1
        if ( ixm1 .LT. 1) ixm1 = ixm1 + nzones
        ixm2 = ix - 2
        if ( ixm2 .LT. 1) ixm2 = ixm2 + nzones

        IF ( ischeme .EQ. 1) THEN
! Lax-Friedrichs

          rhs ( ix) = 0.5 * ( u (ixp1) + u (ixm1))
     1            - 0.5 * mu * ( u (ixp1) - u (ixm1))


        ELSE IF ( ischeme .EQ. 2) THEN
! Lax-Wendroff

          rhs ( ix) = u (ix) - 0.5 * mu * ( u (ixp1) - u (ixm1))
     1            + 0.5 * mu**2 * ( u (ixp1) - 2.0 * u (ix) + u (ixm1))

        ELSE IF ( ischeme .EQ. 3) THEN
! Two-step Runge-Kutta (Heun)

          rhs ( ix) = u (ix) - 0.5 * mu * ( u (ixp1) - u (ixm1))
     1         + 0.125 * mu**2 * ( u (ixp2) - 2.0 * u (ix) + u (ixm2))

        ELSE
! First order upwind

          rhs ( ix) = u (ix) - mu * ( u (ix) - u (ixm1))

        END IF

        END DO ! End of "DO ix =".

      u ( :) = rhs ( :)

! ------------------------------

      t = t + dt  ! Increment time.

        IF ( t .GT. tstop) THEN
! Write out data if we have reached the final time and exit.

!          WRITE ( 6, *)"-------------- Output --------------------"

#if ( 0 == 1)
          DO ix = 1, nzones
          WRITE ( 6, *) u ( ix)
          END DO

          DO ix = 1, nzones
          WRITE ( 6, *) ( - 0.5 + ( ix - 0.5) * dx)
          END DO
#else
          WRITE ( 6, *) ( u ( ix), ix = 1, nzones)
#endif

        EXIT
        END IF

! ----------------------------------------------------------------------

      END DO ! End of "DO istep =".

! ----------------------------------------------------------------------

      STOP

      END PROGRAM ONED_SCALAR_ADVECTION

! ----------------------------------------------------------------------

