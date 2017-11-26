! ----------------------------------------------------------------------

      PROGRAM SLOPE_LIMITER_SCHEME

! Solves the advection equation in periodic geometry using slopes that are
! provided to each zone. The slopes may or may not be limited slopes.
! The Gaussian pulse and top-hat profile are treated.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix, ixp1, ixm1, islope, nzones, ntstep, istep, iproblem

      REAL mu, c, t, dt, dx, tstop, timage

      PARAMETER ( nzones = 100)

      REAL, DIMENSION ( nzones) :: x, un, dun, unp1

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

      OPEN ( unit = 6, file = "rhoa0001")

! islope controls how the slopes are built.
! islope = 1(right-biased); 2(left-biased), 3(minmod); 4(MC-beta)
      islope = 4

      iproblem = 1 ! iproblem = 1,2 for Gaussian, square pulse.

      mu = 0.8 ! CFL #

      c = 1.0  ! Advection speed.

      dx = 1.0 / REAL ( nzones)   ! Zone size

      dt = mu * dx / c  ! Timestep

      t = 0.0

! Stopping time. 0.0, 1.0, 2.0 for Gaussian; 0.0, 0.25 and 0.75 for Sq. pulse.
      tstop = 1.0

      timage = 1.0

      ntstep = 35000 ! # of timesteps.

      mcb = 2.0    ! ranges between 1 and 2.

!      WRITE ( 6, *)"islope = ", islope," iproblem = ", iproblem,
!     1             " c = ", c," mu = ", mu," mcb = ", mcb
!      WRITE ( 6, *)"dx = ", dx," dt = ", dt," tstop = ", tstop,
!     1             " ntstep = ", ntstep

! ----------------------------------------------------------------------

! Initialize the data with a Gaussian or square wave pulse.
! Periodic boundary conditions.

      DO ix = 1, nzones

      x ( ix) = - 0.5 + ( ix - 0.5) * dx

      IF ( iproblem .EQ. 1) THEN
! Gaussian pulse set up here.

        un ( ix) = EXP ( -( x ( ix) / 0.1)**2)

      ELSE
! Square pulse set up here.

        IF ( ( -0.05 .LE. x ( ix)) .AND. ( x ( ix) .LE. 0.05) ) THEN
        un ( ix) = 1.0
        ELSE
        un ( ix) = 0.0
        END IF

      END IF

      END DO ! End of "DO ix ="

! ----------

! Write out the initial profile.

!          WRITE ( 6, *)"--------- Initial Output --------------"

!          DO ix = 1, nzones
!          WRITE ( 6, *) un ( ix)
!          END DO

!          DO ix = 1, nzones
!          WRITE ( 6, *) x ( ix)
!          END DO

! ----------------------------------------------------------------------

! Take "ntstep" timesteps.

      DO istep = 1, ntstep

! ----------------------------------------------------------------------

! Make slopes using choice of limiter. Various choices are possible based
! on "islope". Periodic geometry is enforced.

        DO ix = 1, nzones

        ixp1 = ix + 1
        if ( ixp1 .GT. nzones) ixp1 = ixp1 - nzones
        ixm1 = ix - 1
        if ( ixm1 .LT. 1) ixm1 = ixm1 + nzones

        IF ( islope == 1) THEN

          dun ( ix) = un ( ixp1) - un ( ix)

        ELSE IF ( islope == 2) THEN

          dun ( ix) = un ( ix) - un ( ixm1)

        ELSE IF ( islope == 3) THEN
! Min-mod limiter.

          dun ( ix) = MIN_MOD ( un ( ixp1) - un ( ix),
     1                          un ( ix) - un ( ixm1) )

        ELSE
! MC-beta limiter

          dun ( ix) = MC_BETA_LIMITER ( un ( ixp1) - un ( ix),
     1                                  un ( ix) - un ( ixm1), mcb )

        END IF

        END DO ! End of "DO ix ="

! ------------------------------

! Make update of "un" and store in "unp1". Apply slope limiter scheme.
! This is just eqn. (3.7) of text with slopes
! that may or may not be limited. Periodic geometry is enforced.

        DO ix = 1, nzones

        ixp1 = ix + 1
        if ( ixp1 .GT. nzones) ixp1 = ixp1 - nzones
        ixm1 = ix - 1
        if ( ixm1 .LT. 1) ixm1 = ixm1 + nzones

        unp1 ( ix) = un ( ix) - mu * ( un ( ix) - un ( ixm1) )
     1        - 0.5 * mu * ( 1.0 - mu) * ( dun ( ix) - dun ( ixm1) )

        END DO ! End of "DO ix ="

! ------------------------------

! Lastly, hand back the new timestep to "un".

        un ( :) = unp1 ( :)

! ------------------------------

      t = t + dt  ! Increment time.

! ----------

        IF ( ( t .GE. timage) .AND. ( t .LT. timage + dt) ) THEN

!          WRITE ( 6, *)"--------- Interim Output --------------"

!          DO ix = 1, nzones
!          WRITE ( 6, *) un ( ix)
!          END DO

!          DO ix = 1, nzones
!          WRITE ( 6, *) x ( ix)
!          END DO

        END IF

! ----------

        IF ( t .GE. tstop) THEN
! Write out data if we have reached the final time and exit.

!          WRITE ( 6, *)"--------- Final Output --------------"

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
        EXIT
        END IF

! ----------------------------------------------------------------------

      END DO ! End of "DO istep =".

! ----------------------------------------------------------------------

      STOP

      END PROGRAM SLOPE_LIMITER_SCHEME

! ----------------------------------------------------------------------

