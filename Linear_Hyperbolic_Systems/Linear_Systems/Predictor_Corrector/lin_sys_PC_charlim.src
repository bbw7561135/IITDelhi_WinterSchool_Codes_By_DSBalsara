! ----------------------------------------------------------------------

      PROGRAM LIN_SYS_PC_CHARLIM

! Single-step, predictor-corrector strategy for updating linear systems.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix, ixp1, ixm1, j, k, l, islope, nzones, ntstep, nsys,
     1        istep, iproblem

      REAL mu, t, dt, dx, tstop, vx0, c0, rho0, umin, umax,
     1     tempa, tempb

      PARAMETER ( nzones = 70, nsys = 3)

      REAL, DIMENSION ( nsys) :: eigenval,
     1      eigenval_pl, eigenval_mi, eigenval_abs

      REAL, DIMENSION ( nsys, nsys) :: rt_eigenvec, lf_eigenvec,
     1      a_tot, a_pl, a_mi, a_abs, m_1, m_2, m_3

      REAL, DIMENSION ( -1: nzones + 2) :: x

      REAL, DIMENSION ( -1: nzones + 2, nsys) :: un, dun,
     1      char_un, dchar_un, un_rt, un_lf, flx, tu

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

! Set up the vectors and matrices for eigenvalues and eigenvectors.
! Eigenvalues must be ordered from smallest to largest.

      IF ( nsys .EQ. 2) THEN
! We solve the 2x2 linear system given by {{0,1},{1,0}}.

      a_tot ( 1, 1) = 0.0
      a_tot ( 1, 2) = 1.0
      a_tot ( 2, 1) = 1.0
      a_tot ( 2, 2) = 0.0

      eigenval ( 1) = -1
      eigenval ( 2) = 1

      rt_eigenvec ( 1, 1) = 1.0 / SQRT ( 2.0)
      rt_eigenvec ( 2, 1) = -1.0 / SQRT ( 2.0)
      rt_eigenvec ( 1, 2) = 1.0 / SQRT ( 2.0)
      rt_eigenvec ( 2, 2) = 1.0 / SQRT ( 2.0)

      lf_eigenvec ( 1, 1) = 1.0 / SQRT ( 2.0)
      lf_eigenvec ( 1, 2) = -1.0 / SQRT ( 2.0)
      lf_eigenvec ( 2, 1) = 1.0 / SQRT ( 2.0)
      lf_eigenvec ( 2, 2) = 1.0 / SQRT ( 2.0)

      ELSE IF ( nsys .EQ. 3) THEN

#if ( 0 == 1)

! We solve the 3x3 linear acoustics system given by
! {{vx0, rho0 * c0**2, 0},{1/rho0,vx0,0},{0,0,vx0}}.

      vx0 = 0.25
      c0 = 0.75
      rho0 = 1.0

      a_tot ( 1, 1) = vx0
      a_tot ( 1, 2) = rho0 * c0**2
      a_tot ( 1, 3) = 0.0

      a_tot ( 2, 1) = 1.0 / rho0
      a_tot ( 2, 2) = vx0
      a_tot ( 2, 3) = 0.0

      a_tot ( 3, 1) = 0.0
      a_tot ( 3, 2) = 0.0
      a_tot ( 3, 3) = vx0

      eigenval ( 1) = vx0 - c0
      eigenval ( 2) = vx0
      eigenval ( 3) = vx0 + c0

      rt_eigenvec ( 1, 1) = - rho0 * c0**2
      rt_eigenvec ( 2, 1) = 1.0
      rt_eigenvec ( 3, 1) = 0.0

      rt_eigenvec ( 1, 2) = 0.0
      rt_eigenvec ( 2, 2) = 0.0
      rt_eigenvec ( 3, 2) = 1.0

      rt_eigenvec ( 1, 3) = rho0 * c0**2
      rt_eigenvec ( 2, 3) = 1.0
      rt_eigenvec ( 3, 3) = 0.0

      lf_eigenvec ( 1, 1) = - 1.0 / ( 2.0 * rho0 * c0**2)
      lf_eigenvec ( 1, 2) = 0.5
      lf_eigenvec ( 1, 3) = 0.0

      lf_eigenvec ( 2, 1) = 0.0
      lf_eigenvec ( 2, 2) = 0.0
      lf_eigenvec ( 2, 3) = 1.0

      lf_eigenvec ( 3, 1) = 1.0 / ( 2.0 * rho0 * c0**2)
      lf_eigenvec ( 3, 2) = 0.5
      lf_eigenvec ( 3, 3) = 0.0

#else

! We solve the 3x3 linearized Euler system given by
! {{vx0, rho0, 0},{0,vx0,1/rho0},{0,rho0*c0**2,vx0}}.

      vx0 = 0.25
      c0 = 0.75
      rho0 = 1.0

      a_tot ( 1, 1) = vx0
      a_tot ( 1, 2) = rho0
      a_tot ( 1, 3) = 0.0

      a_tot ( 2, 1) = 0.0
      a_tot ( 2, 2) = vx0
      a_tot ( 2, 3) = 1.0 / rho0

      a_tot ( 3, 1) = 0.0
      a_tot ( 3, 2) = rho0 * c0**2
      a_tot ( 3, 3) = vx0

      eigenval ( 1) = vx0 - c0
      eigenval ( 2) = vx0
      eigenval ( 3) = vx0 + c0

      rt_eigenvec ( 1, 1) = rho0
      rt_eigenvec ( 2, 1) = - c0
      rt_eigenvec ( 3, 1) = rho0 * c0**2

      rt_eigenvec ( 1, 2) = 1.0
      rt_eigenvec ( 2, 2) = 0.0
      rt_eigenvec ( 3, 2) = 0.0

      rt_eigenvec ( 1, 3) = rho0
      rt_eigenvec ( 2, 3) = c0
      rt_eigenvec ( 3, 3) = rho0 * c0**2

      lf_eigenvec ( 1, 1) = 0.0
      lf_eigenvec ( 1, 2) = - 0.5 / c0
      lf_eigenvec ( 1, 3) = 0.5 / ( rho0 * c0**2)

      lf_eigenvec ( 2, 1) = 1.0
      lf_eigenvec ( 2, 2) = 0.0
      lf_eigenvec ( 2, 3) = - 1.0 / c0**2

      lf_eigenvec ( 3, 1) = 0.0
      lf_eigenvec ( 3, 2) = 0.5 / c0
      lf_eigenvec ( 3, 3) = 0.5 / ( rho0 * c0**2)

#endif

      ELSE

      STOP'nsys out of range'

      END IF

! ------------------------------

! Build A^+ i.e. a_pl; A^- i.e. a_mi; |A| i.e. a_abs; A itself i.e. a_tot.

      DO j = 1, nsys
      eigenval_pl ( j) = AMAX1 ( eigenval ( j), 0.0)
      eigenval_mi ( j) = AMIN1 ( eigenval ( j), 0.0)
      eigenval_abs ( j) = ABS ( eigenval ( j))
      END DO

      m_1 ( :, :) = 0.0
      m_2 ( :, :) = 0.0
      m_3 ( :, :) = 0.0

      DO j = 1, nsys
        DO k = 1, nsys
        m_1 ( j, k) = rt_eigenvec ( j, k) * eigenval_pl ( k)
        m_2 ( j, k) = rt_eigenvec ( j, k) * eigenval_mi ( k)
        m_3 ( j, k) = rt_eigenvec ( j, k) * eigenval_abs ( k)
        END DO
      END DO

      a_pl ( :, :) = 0.0
      a_mi ( :, :) = 0.0
      a_abs ( :, :) = 0.0

      DO j = 1, nsys
        DO k = 1, nsys
          DO l = 1, nsys
          a_pl ( j, k) = a_pl ( j, k)
     1                 + m_1 ( j, l) * lf_eigenvec ( l, k)
          a_mi ( j, k) = a_mi ( j, k)
     1                 + m_2 ( j, l) * lf_eigenvec ( l, k)
          a_abs ( j, k) = a_abs ( j, k)
     1                  + m_3 ( j, l) * lf_eigenvec ( l, k)
          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "rhoa0001")

! islope controls how the slopes are built. islope = 1,2 for min-mod & MC-beta.
      islope = 1

      iproblem = 1

      mu = 0.6 ! CFL #

      dx = 1.0 / REAL ( nzones)   ! Zone size

      dt = mu * dx / AMAX1 ( ABS ( eigenval ( 1)),
     1                       ABS ( eigenval ( nsys)) )

      t = 0.0

      tstop = 0.4    ! Stopping time.

      ntstep = 35000 ! # of timesteps.

      mcb = 2.0    ! ranges between 1 and 2.

!      WRITE ( 6, *)"islope = ", islope," iproblem = ", iproblem,
!     1             " mu = ", mu," mcb = ", mcb
!      WRITE ( 6, *)"dx = ", dx," dt = ", dt," tstop = ", tstop,
!     1             " ntstep = ", ntstep

! ----------------------------------------------------------------------

! Initialize the data. Set ghost zone values too.

      umin = 1.0e32
      umax = -1.0e32

      DO ix = 1 - 2, nzones + 2

      x ( ix) = - 0.5 + ( ix - 0.5) * dx

      IF ( x ( ix) .LE. 0.0) THEN

        un ( ix, 1) = 0.2
        un ( ix, 2) = 0.3
        IF ( nsys .EQ. 3) un ( ix, 3) = 0.1

      ELSE

        un ( ix, 1) = - 0.1
        un ( ix, 2) = 0.4
        IF ( nsys .EQ. 3) un ( ix, 3) = 0.3

      END IF

      umin = AMIN1 ( umin, un ( ix, 1))
      umin = AMIN1 ( umin, un ( ix, 2))
      IF ( nsys .EQ. 3) umin = AMIN1 ( umin, un ( ix, 3))
      umax = AMAX1 ( umax, un ( ix, 1))
      umax = AMAX1 ( umax, un ( ix, 2))
      IF ( nsys .EQ. 3) umax = AMAX1 ( umax, un ( ix, 3))

      END DO ! End of "DO ix ="

!      WRITE ( 6, *)"umin, umax = ", umin, umax

! ----------------------------------------------------------------------

! Take "ntstep" timesteps.

      DO istep = 1, ntstep

! ----------------------------------------------------------------------

! Make slopes using choice of limiters. Also decide whether one wants to
! limit characteristically or componentwise.

! Limit on the characteristic variables.

! ------------------------------

! Project to make characteristic variables "char_un".

      char_un ( :, :) = 0.0

      DO j = 1, nsys     ! jth characteristic field being built.

        DO k = 1, nsys   ! k is a matrix multiplication index.

          DO ix = 1 - 2, nzones + 2   ! ix runs over zone centers.

          char_un ( ix, j) = char_un ( ix, j)
     1                     + lf_eigenvec ( j, k) * un ( ix, k)

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

      END DO ! End of "DO j ="

! ------------------------------

! Limit the undivided differences in the characteristic variables to make
! "dchar_un".

      DO j = 1, nsys     ! jth characteristic field being limited.

        DO ix = 1 - 1, nzones + 1   ! ix runs over zone centers.

        IF ( islope .EQ. 1) THEN

        dchar_un ( ix, j) =
     1    MIN_MOD ( char_un ( ix + 1, j) - char_un ( ix, j),
     1              char_un ( ix, j) - char_un ( ix - 1, j) )

        ELSE

        dchar_un ( ix, j) =
     1  MC_BETA_LIMITER ( char_un ( ix + 1, j) - char_un ( ix, j),
     1                    char_un ( ix, j) - char_un ( ix - 1, j), mcb)

        END IF

        END DO ! End of "DO ix ="

      END DO ! End of "DO j ="

! ------------------------------

! Project the differences in the characteristic variables back to physical
! space using "rt_eigenvec". This gives "dun". This step not needed in LW.

      dun ( :, :) = 0.0


      DO j = 1, nsys     ! jth component of undivided difference "dun".

        DO k = 1, nsys   ! k is a matrix multiplication index.

          DO ix = 1 - 1, nzones + 1   ! ix runs over zone centers.

          dun ( ix, j) = dun ( ix, j)
     1                 + rt_eigenvec ( j, k) * dchar_un ( ix, k)

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

      END DO ! End of "DO j ="

! ----------------------------------------------------------------------

! Build the fluxes using the Roe type Riemann solver.

! Predictor part where we add in the time-dependence to get the left
! and right states at dt/2.

      tempa = 1.0 - AMAX1 ( eigenval ( nsys), 0.0) * dt / dx
      tempb = 1.0 + AMIN1 ( eigenval ( 1), 0.0) * dt / dx

!      tempa = 1.0
!      tempb = 1.0

      DO j = 1, nsys     ! jth component of "un_rt" and "un_lf"

        DO ix = 1 - 2, nzones + 2

        un_rt ( ix, j) = un ( ix, j) + 0.5 * tempb * dun ( ix, j)
        un_lf ( ix, j) = un ( ix, j) - 0.5 * tempa * dun ( ix, j)

        END DO ! End of "DO ix ="

      END DO ! End of "DO j ="


      DO j = 1, nsys     ! jth component of "un_rt" and "un_lf"

        DO k = 1, nsys   ! k loops over all the rt. eigenvectors being added in.

          DO ix = 1 - 1, nzones + 1    ! ix runs over zone centers.

          tempa = - (0.5 * eigenval ( k) * dt / dx) * dchar_un ( ix, k)
     1          * rt_eigenvec ( j, k)

          un_rt ( ix, j) = un_rt ( ix, j) - tempa

          un_lf ( ix, j) = un_lf ( ix, j) - tempa

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

      END DO ! End of "DO j ="


! ----------------------------------------------------------------------

! Corrector part of Predictor-Corrector method.

      flx ( :, :) = 0.0

      DO j = 1, nsys     ! jth component of flux "flx"

        DO k = 1, nsys   ! k is a matrix multiplication index.

          DO ix = 1 - 1, nzones    ! ix runs over zone boundaries.

          flx ( ix, j) = flx ( ix, j)
     1  + 0.5 * a_tot ( j, k) * ( un_lf ( ix + 1, k) + un_rt ( ix, k) )
     1  - 0.5 * a_abs ( j, k) * ( un_lf ( ix + 1, k) - un_rt ( ix, k) )

          END DO ! End of "DO ix ="

        END DO ! End of "DO k ="

      END DO ! End of "DO j ="

! ------------------------------

! One-step update.

      DO j = 1, nsys     ! jth component of "un" being updated.

        DO ix = 1, nzones    ! ix runs over zone centers.

        un ( ix, j) = un ( ix, j)
     1     - ( flx ( ix, j) - flx ( ix - 1, j) ) * dt / dx

        END DO ! End of "DO ix ="

      END DO ! End of "DO j ="

! ------------------------------

! Impose boundary values on "un" here if desired.

      un ( 0, :) = un ( 1, :)
      un ( -1, :) = un ( 1, :)

      un ( nzones + 1, :) = un ( nzones, :)
      un ( nzones + 2, :) = un ( nzones, :)

! ----------------------------------------------------------------------

#if ( 0 == 1)
! Write out diagnostics.

        umin = 1.0e32
        umax = -1.0e32

        DO j = 1, nsys
        DO ix = 1, nzones
        umin = AMIN1 ( umin, un ( ix, j))
        umax = AMAX1 ( umax, un ( ix, j))
        END DO
        END DO

!        WRITE ( 6, *)"from output : time, umin, umax = ",
!     1                t, umin, umax

#endif

! ----------------------------------------------------------------------

      t = t + dt  ! Increment time.

        IF ( t .GT. tstop) THEN
! Write out data if we have reached the final time and exit.

        umin = 1.0e32
        umax = -1.0e32

        DO j = 1, nsys
        DO ix = 1, nzones
        umin = AMIN1 ( umin, un ( ix, j))
        umax = AMAX1 ( umax, un ( ix, j))
        END DO
        END DO

!        WRITE ( 6, *)"from output : umin, umax = ", umin, umax

!        WRITE ( 6, *)"-------------- Output --------------------"

#if ( 0 == 1)
        DO j = 1, nsys
          DO ix = 1, nzones
          WRITE ( 6, *) un ( ix, j)
          END DO
        END DO

        DO ix = 1, nzones
        WRITE ( 6, *) x ( ix)
        END DO
#else
          WRITE ( 6, *) ( ( un ( ix, j), ix = 1, nzones), j = 1, nsys)
#endif

        GO TO 100

        END IF

! ----------------------------------------------------------------------

      END DO ! End of "DO istep =".

100   CONTINUE

! ----------------------------------------------------------------------

      STOP

      END PROGRAM LIN_SYS_PC_CHARLIM

! ----------------------------------------------------------------------


