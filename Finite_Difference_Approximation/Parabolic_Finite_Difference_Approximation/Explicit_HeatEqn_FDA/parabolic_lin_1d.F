! ----------------------------------------------------------------------

      PROGRAM ONED_HEAT_CONDUCTION

! Solves linear heat conduction problem in 1d. Can range from fully explicit
! to fully implicit depending on "alpha".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzones, ntstep, istep, ix

      PARAMETER ( nzones = 64)

      REAL x, dx, t, dt, tstop, mu, sigma, alpha,
     1     a_lf, b_lf, c_lf, a_rt, b_rt, c_rt

      REAL, DIMENSION ( 0: nzones) :: u, rhs, a, b, c

! ----------------------------------------------------------------------

! Declare some scalar variables.

      OPEN ( unit = 6, file = "rhoa0001")

      mu = 0.4  ! Controls the parabolic timestep. 0.5 is limit for explicit.

      alpha = 1.0  ! alpha=1.0 ==> fully explicit; alpha=0.0 ==> fully implicit.

      sigma = 1.0  ! heat conduction coefficient.

      tstop = 0.01  ! Stopping time.

      dx = 1.0 / REAL ( nzones)   ! Make zone size on a unit 1d domain.

      dt = mu * dx**2 / sigma
      t = 0.0

      ntstep = 35000 ! # of timesteps.

!      WRITE ( 6, *)"mu = ", mu," alpha = ", alpha," sigma = ", sigma
!      WRITE ( 6, *)"tstop = ", tstop," dx = ", dx," dt = ", dt,
!     1             " ntstep = ", ntstep

! ----------------------------------------------------------------------

! Initialize the data with a step function.
! Dirichlet boundary conditions with zero at boundaries.

      DO ix = 0, nzones
      x = ix * dx

        IF ( ( 0.4 .LE. x) .AND. ( x .LE. 0.6) ) THEN
        u ( ix) = 1.0
        ELSE
        u ( ix) = 0.0
        END IF

      END DO

! ------------------------------

! Left boundary condition with "a_lf * du_dx + b_lf * u = c_lf".
      a_lf = 0.0; b_lf = 1.0; c_lf = 0.0

! Right boundary condition with "a_rt * du_dx + b_rt * u = c_rt".
      a_rt = 0.0; b_rt = 1.0; c_rt = 0.0

! ----------------------------------------------------------------------

! Take "ntstep" timesteps.

      DO istep = 1, ntstep

! ----------------------------------------------------------------------

        DO ix = 1, nzones - 1
        rhs ( ix) = u ( ix) + alpha * mu
     1            * ( u ( ix + 1) - 2.0 * u ( ix) + u ( ix - 1) )
        END DO ! End of "DO ix =".

! ------------------------------

        IF ( ABS ( alpha - 1.0) .LT. 1.0e-8) THEN
! Fully explicit, no need for matrix inversion.
! Dirichlet boundary conditions assumed.

        u ( 1: nzones - 1) = rhs ( 1: nzones - 1)

        ELSE
! Semi-implicit or fully implicit.
! We use the band matrix solver "tridag" from Numerical Recipes.

        b ( 0) = b_lf * dx - a_lf
        a ( 0) = 0.0
        c ( 0) = a_lf

        b ( nzones) = b_rt * dx + a_rt
        a ( nzones) = - a_rt
        c ( nzones) = 0.0

        b ( 1: nzones - 1) = 1.0 + 2.0 * ( 1.0 - alpha) * mu ! diagonal terms
        a ( 1: nzones - 1) = - ( 1.0 - alpha) * mu         ! subdiagonal terms
        c ( 1: nzones - 1) = - ( 1.0 - alpha) * mu        ! superdiagonal terms

        rhs ( 0) = c_lf * dx
        rhs ( nzones) = c_rt * dx

        call tridag ( a, b, c, rhs, u, nzones + 1)

        END IF

! ------------------------------

      t = t + dt  ! Increment time.

        IF ( t .GT. tstop) THEN
! Write out data if we have reached the final time and exit.

!          WRITE ( 6, *)"-------------- Output --------------------"

#if ( 0 == 1)
          DO ix = 0, nzones
          WRITE ( 6, *) u ( ix)
          END DO

          DO ix = 0, nzones
          WRITE ( 6, *) ix * dx
          END DO
#else
        WRITE ( 6, *) ( u ( ix), ix = 0, nzones )
#endif

        EXIT
        END IF

! ----------------------------------------------------------------------

      END DO ! End of "DO istep =".

! ----------------------------------------------------------------------

      STOP

      END PROGRAM ONED_HEAT_CONDUCTION

! ----------------------------------------------------------------------

      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      REAL bet,gam(NMAX)
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

! ----------------------------------------------------------------------

