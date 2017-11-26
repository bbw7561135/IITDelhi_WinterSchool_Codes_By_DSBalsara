! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_BLAST_3D_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 3d Hydrodynamic blast problem.
! Can be collapsed to a 2d problem by switching off one of the directions.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,
     1     vlz1, vlz2, 
     1     radius, xcenter, ycenter, zcenter,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u


      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! Stop this problem at a time of "timestop = 0.01".

      rho1 = 1.0
      rho2 = 1.0

      prs1 = 1000.0
      prs2 = 0.1

      vlx1 = 0.0
      vlx2 = 0.0

      vly1 = 0.0
      vly2 = 0.0

      vlz1 = 0.0
      vlz2 = 0.0

      radius = 0.1

! Make the center of the domain so that the explosion can be set up there.

      xcenter = 0.5 * ( base_grid_xmin + base_grid_xmax)
      ycenter = 0.5 * ( base_grid_ymin + base_grid_ymax)
      zcenter = 0.5 * ( base_grid_zmin + base_grid_zmax)

! ----------------------------------------------------------------------


! ------------------------------

! Now set up interior. Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          tempaa = 0.0

          IF ( ix1 .GT. 1) tempaa = tempaa
     1                            + ( xc_save ( ixx) - xcenter)**2
          IF ( iy1 .GT. 1) tempaa = tempaa
     1                            + ( yc_save ( iyy) - ycenter)**2
          IF ( iz1 .GT. 1) tempaa = tempaa
     1                            + ( zc_save ( izz) - zcenter)**2

          tempaa = SQRT ( tempaa)

          IF ( tempaa .LE. radius) THEN

          u ( ixx, iyy, izz, 1, 1) = rho1

          u ( ixx, iyy, izz, 2, 1) = rho1 * vlx1
          u ( ixx, iyy, izz, 3, 1) = rho1 * vly1
          u ( ixx, iyy, izz, 4, 1) = rho1 * vlz1

          u ( ixx, iyy, izz, 5, 1) = prs1 / (gamma_euler - 1.0)
     1         + 0.5 * rho1 * ( vlx1**2 + vly1**2 + vlz1**2 )

#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz, 6: NFIELD_EULER + NFLUID_EULER, 1) = rho1
#endif

          tempab = SQRT ( vlx1**2 + vly1**2 + vlz1**2)
          tempac = gamma_euler * prs1
          tempac = SQRT ( tempac / rho1) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

          ELSE

          u ( ixx, iyy, izz, 1, 1) = rho2

          u ( ixx, iyy, izz, 2, 1) = rho2 * vlx2
          u ( ixx, iyy, izz, 3, 1) = rho2 * vly2
          u ( ixx, iyy, izz, 4, 1) = rho2 * vlz2

          u ( ixx, iyy, izz, 5, 1) = prs2 / (gamma_euler - 1.0)
     1         + 0.5 * rho2 * ( vlx2**2 + vly2**2 + vlz2**2 )

#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz, 6: NFIELD_EULER + NFLUID_EULER, 1) = 0.0
#endif

          tempab = SQRT ( vlx2**2 + vly2**2 + vlz2**2)
          tempac = gamma_euler * prs2
          tempac = SQRT ( tempac / rho2) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

          END IF

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_BLAST_3D_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_JET_3D

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 3d Euler jet problem.

! Can be collapsed to a 2d problem by switching off one of the directions.
! When so collapsed, it yields the 2d Euler slab jet problem. We contrive
! to have the sound speed in the jet to be unity. The jet flows in the
! z-direction, so we can set the z-velocity equal to the desired Mach number.

! To obtain that problem on a 200x1x400 zone mesh on the
! interval [-.5,.5]x[-.5,.5]x[-.5,1.5], set "gamma_euler = 1.4",
! "rho1 = gamma_euler",
! "rho2 = 0.1 * gamma_euler", "prs1 = 1.0, prs2 = 1.0", "vlz1 = 80.0",
! "radius = 0.05".

! Run till "timestop = 0.030" for Mach 80 jet.
! Run till "timestop = 0.0030" for Mach 800 jet.
! Run till "timestop = 0.0012" for Mach 2000 jet.

! The timestep is controlled by letting the jet pop a few zones into the
! computational domain. This ensures that the CFL number of the interior
! solution (without considering boundaries) properly represents the desired
! CFL number.

! Put this line when setting the lower boundary of "bcarr_z". The other
! boundaries can remain continuitive.
!      IF ( ( SQRT ( xc_save ( ixx)**2 + yc_save ( iyy)**2) .LE. 0.05)
!     1   .AND. ( ABS ( zmin - base_grid_zmin) .LT. smallnum) )
!     1   bcarr_z ( ixx, iyy, 1) = 1

! In "PAD_EULER_BOUNDARY_Z", for the inflow at lower z-boundary, set:
!      rho_tmp = gamma_euler
!      momx_tmp = 0.0
!      momy_tmp = 0.0
!      momz_tmp = 80.0 * rho_tmp ! Reset this var to change jet's inflow speed.
!      eng_tmp = 1.0 / ( gamma_euler - 1.0)
!     1        + 0.5 * momz_tmp**2 / rho_tmp

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,
     1     vlz1, vlz2,
     1     radius, xcenter, ycenter, zcenter,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! This choice ensures that we have a dense jet and the sound speed in
! the jet is unity.

      rho1 = gamma_euler
      rho2 = 0.1 * gamma_euler

      prs1 = 1.0
      prs2 = 1.0

      vlx1 = 0.0
      vlx2 = 0.0

      vly1 = 0.0
      vly2 = 0.0

      vlz1 = 80.0
      vlz2 = 0.0

      radius = 0.05

! Make the center of the domain so that the explosion can be set up there.

      xcenter = 0.5 * ( base_grid_xmin + base_grid_xmax)
      ycenter = 0.5 * ( base_grid_ymin + base_grid_ymax)
      zcenter = 0.5 * ( base_grid_zmin + base_grid_zmax)

! ----------------------------------------------------------------------

! Now set up interior. Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          tempaa = 0.0

          IF ( ix1 .GT. 1) tempaa = tempaa
     1                            + ( xc_save ( ixx) - xcenter)**2
          IF ( iy1 .GT. 1) tempaa = tempaa
     1                            + ( yc_save ( iyy) - ycenter)**2

          tempaa = SQRT ( tempaa)

          IF ( ( tempaa .LE. radius) .AND. ( izz .LE. 5) ) THEN

          u ( ixx, iyy, izz, 1, 1) = rho1

          u ( ixx, iyy, izz, 2, 1) = rho1 * vlx1
          u ( ixx, iyy, izz, 3, 1) = rho1 * vly1
          u ( ixx, iyy, izz, 4, 1) = rho1 * vlz1

          u ( ixx, iyy, izz, 5, 1) = prs1 / (gamma_euler - 1.0)
     1         + 0.5 * rho1 * ( vlx1**2 + vly1**2 + vlz1**2 )

          tempab = SQRT ( vlx1**2 + vly1**2 + vlz1**2)
          tempac = gamma_euler * prs1
          tempac = SQRT ( tempac / rho1) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz,
     1        NFIELD_EULER + 1: NFIELD_EULER + NFLUID_EULER, 1) = rho1
#endif

          ELSE

          u ( ixx, iyy, izz, 1, 1) = rho2

          u ( ixx, iyy, izz, 2, 1) = rho2 * vlx2
          u ( ixx, iyy, izz, 3, 1) = rho2 * vly2
          u ( ixx, iyy, izz, 4, 1) = rho2 * vlz2

          u ( ixx, iyy, izz, 5, 1) = prs2 / (gamma_euler - 1.0)
     1         + 0.5 * rho2 * ( vlx2**2 + vly2**2 + vlz2**2 )

          tempab = SQRT ( vlx2**2 + vly2**2 + vlz2**2)
          tempac = gamma_euler * prs2
          tempac = SQRT ( tempac / rho2) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz,
     1        NFIELD_EULER + 1: NFIELD_EULER + NFLUID_EULER, 1) = 0.0
#endif

          END IF

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_JET_3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ISENTROPIC_VORTEX_XY_EULER
     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up an isentropic hydrodynamical vortex in 2d. Default is xy-plane
! with iz1 = 1. The problem is catalogued in Balsara & Shu (2000) JCP.
! Use a domain of [-5,5]x[-5,5]x[-5,5] on a range of meshes
! going from 64x64x1 to 512x512x1 zones for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma_euler = 1.4" suggested. Suggested stopping time is
! "timestop = 10.0". The vortex propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Can analyze the
! accuracy of the problem by comparing the numerical solution with
! the analytical solution.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho0, prs0, entropy0, vmag, vlx0, vly0, delta_t,

     1     radius, tempaa, tempab, tempac,
     1     costheta, sintheta, velrms, bfldrms,
     1     VLXFN, VLYFN, TEMPFN,
     1     x, y,
     1     xlo, xhi, ylo, yhi,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     prsmin, prsmax

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------

      VLXFN ( vmag, x, y) = - vmag * y * EXP ( 0.5*(1.0-x**2-y**2))
      VLYFN ( vmag, x, y) = vmag * x * EXP ( 0.5*(1.0-x**2-y**2))

      TEMPFN ( pi, vmag, gamma_euler, x, y) =
     1               - ( gamma_euler - 1.0) * vmag**2
     1               * EXP (1.0-x**2-y**2) / ( 2.0 * gamma_euler)

! ----------------------------------------------------------------------

! Default parameters for isentropic hydrodynamical vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0
      entropy0 = prs0 / rho0**gamma_euler

      vlx0 = 1.0
      vly0 = 1.0

      vmag = 5.0 / ( 2.0 * pi)

! ----------------------------------------------------------------------

      dtcur = 1.0e30

      prsmin = 1.0e30
      prsmax = -1.0e30

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = xb_save ( ixx - 1)
          xhi = xb_save ( ixx)
          ylo = yb_save ( iyy - 1)
          yhi = yb_save ( iyy)

          vlxuse = vlx0 + 0.25 * VLXFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLXFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLXFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLXFN ( vmag, xhi, yhi)

          vlyuse = vly0 + 0.25 * VLYFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLYFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLYFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLYFN ( vmag, xhi, yhi)

          vlzuse = 0.0

          delta_t = 0.25 * TEMPFN ( pi, vmag, gamma_euler, xlo, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma_euler, xhi, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma_euler, xlo, yhi)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma_euler, xhi, yhi)

          rhouse = ( ( ( prs0 / rho0) + delta_t)
     1               / entropy0)**(1.0/(gamma_euler-1.0) )

          prsuse = rhouse * ( ( prs0 / rho0) + delta_t)

          prsmin = MIN ( prsmin, prsuse)
          prsmax = MAX ( prsmax, prsuse)

          enguse = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

          u ( ixx, iyy, izz, 1, 1) = rhouse
          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse
          u ( ixx, iyy, izz, 5, 1) = enguse

! ----------


          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2

          tempab = 
     1            gamma_euler * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ISENTROPIC_VORTEX_XY_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ERROR_ISENTROPIC_VORTEX_XY_EULER
     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Evaluates the error in the isentropic hydrodynamical vortex in 2d. This
! subroutine should be placed in the timestep loop & it will
! keep writing out the erorr. Default is xy-plane
! with iz1 = 1. The problem is catalogued in Balsara & Shu (2000) JCP.
! Use a domain of [-5,5]x[-5,5]x[-5,5] on a 120x120x1 zone mesh for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma_euler = 1.4" suggested. Suggested stopping time is
! "timestop = 10.0". The vortex propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Can analyze the
! accuracy of the problem by comparing the numerical solution with
! the analytical solution.

! Some accuracy data is provided below. Notice that high quality reconstruction
! + good RS yields substantially better results on smaller meshes and
! a much faster approach to the theoretical accuracy.

! ----------------------------------------------------------------------
! WENO reconstruction with multid. RS and CFL no. 0.4.
! zones     rhoerr_L1     accuracy_L1     rhoerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64     1.1922E-3                     2.5498E-2
! 128x128   2.7608E-4     2.11            6.0024E-3     2.08
! 256x256   6.6438E-5     2.05            1.3497E-3     2.15
! 512x512   1.6459E-5     2.01            3.1559E-4     2.10
! 1024x1024 4.1068E-6     2.00            7.6091E-5     2.05
! 2048x2048 1.0263E-6     2.01            1.8676E-5     2.03
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz, shift_pattern

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     rho0, prs0, entropy0, vmag, vlx0, vly0, delta_t,

     1     radius, tempaa, tempab, tempac,
     1     costheta, sintheta, velrms, bfldrms,
     1     VLXFN, VLYFN, TEMPFN,
     1     x, y,
     1     xlo, xhi, ylo, yhi,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,

     1     rhoerr_L1, rhoerr_Linf, engerr_L1, engerr_Linf,
     1     vlxerr_L1, vlxerr_Linf, vlyerr_L1, vlyerr_Linf,
     1     vlzerr_L1, vlzerr_Linf

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------

      VLXFN ( vmag, x, y) = - vmag * y * EXP ( 0.5*(1.0-x**2-y**2))
      VLYFN ( vmag, x, y) = vmag * x * EXP ( 0.5*(1.0-x**2-y**2))

      TEMPFN ( pi, vmag, gamma_euler, x, y) =
     1               - ( gamma_euler - 1.0) * vmag**2 
     1                 * EXP (1.0-x**2-y**2) / ( 2.0 * gamma_euler)

! ----------------------------------------------------------------------

! Default parameters for isentropic vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0
      entropy0 = prs0 / rho0**gamma_euler

      vlx0 = 1.0
      vly0 = 1.0

      vmag = 5.0 / ( 2.0 * pi)

! ----------------------------------------------------------------------

! Shift the rotor pattern along till it returns to its original location
! at time of 10.0 units. When that happens, there is no need to shift it.

      rhoerr_L1 = 0.0
      rhoerr_Linf = -1.0e32
      engerr_L1 = 0.0
      engerr_Linf = -1.0e32

      vlxerr_L1 = 0.0
      vlxerr_Linf = -1.0e32
      vlyerr_L1 = 0.0
      vlyerr_Linf = -1.0e32
      vlzerr_L1 = 0.0
      vlzerr_Linf = -1.0e32

      shift_pattern = 1
      IF ( ABS ( time - 10.0) .LT. smallnum) shift_pattern = 0

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          IF ( shift_pattern .EQ. 1) THEN

          xlo = xb_save ( ixx - 1) - vlx0 * time
          IF ( xlo .GE. base_grid_xmax) xlo = xlo
     1                       - ( base_grid_xmax - base_grid_xmin)
          IF ( xlo .LT. base_grid_xmin) xlo = xlo
     1                       + ( base_grid_xmax - base_grid_xmin)
          xhi = xlo + dx_save ( ixx)

          ylo = yb_save ( iyy - 1) - vly0 * time
          IF ( ylo .GE. base_grid_ymax) ylo = ylo
     1                       - ( base_grid_ymax - base_grid_ymin)
          IF ( ylo .LT. base_grid_ymin) ylo = ylo
     1                       + ( base_grid_ymax - base_grid_ymin)
          yhi = ylo + dy_save ( iyy)

          ELSE

          xlo = xb_save ( ixx - 1)
          xhi = xb_save ( ixx)
          ylo = yb_save ( iyy - 1)
          yhi = yb_save ( iyy)

          END IF

! ----------

          vlxuse = vlx0 + 0.25 * VLXFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLXFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLXFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLXFN ( vmag, xhi, yhi)

          vlyuse = vly0 + 0.25 * VLYFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLYFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLYFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLYFN ( vmag, xhi, yhi)

          vlzuse = 0.0

          delta_t = 0.25 * TEMPFN ( pi, vmag, gamma_euler, xlo, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma_euler, xhi, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma_euler, xlo, yhi)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma_euler, xhi, yhi)

          rhouse = ( ( ( prs0 / rho0) + delta_t)
     1               / entropy0)**(1.0/(gamma_euler-1.0) )

          prsuse = rhouse * ( ( prs0 / rho0) + delta_t)

          enguse = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          tempaa = ABS ( u ( ixx, iyy, izz, 1, 1) - rhouse)
          rhoerr_L1 = rhoerr_L1 + tempaa
          rhoerr_Linf = MAX ( rhoerr_Linf, tempaa)

          tempaa = ABS ( u ( ixx, iyy, izz, 5, 1) - enguse)
          engerr_L1 = engerr_L1 + tempaa
          engerr_Linf = MAX ( engerr_Linf, tempaa)


          tempaa = ABS ( u ( ixx, iyy, izz, 2, 1) - rhouse * vlxuse)
          vlxerr_L1 = vlxerr_L1 + tempaa
          vlxerr_Linf = MAX ( vlxerr_Linf, tempaa)

          tempaa = ABS ( u ( ixx, iyy, izz, 3, 1) - rhouse * vlyuse)
          vlyerr_L1 = vlyerr_L1 + tempaa
          vlyerr_Linf = MAX ( vlyerr_Linf, tempaa)

          tempaa = ABS ( u ( ixx, iyy, izz, 4, 1) - rhouse * vlzuse)
          vlzerr_L1 = vlzerr_L1 + tempaa
          vlzerr_Linf = MAX ( vlzerr_Linf, tempaa)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      rhoerr_L1 = rhoerr_L1 / FLOAT ( ix1 * iy1 * iz1)
      engerr_L1 = engerr_L1 / FLOAT ( ix1 * iy1 * iz1)

      vlxerr_L1 = vlxerr_L1 / FLOAT ( ix1 * iy1 * iz1)
      vlyerr_L1 = vlyerr_L1 / FLOAT ( ix1 * iy1 * iz1)
      vlzerr_L1 = vlzerr_L1 / FLOAT ( ix1 * iy1 * iz1)

! ----------------------------------------------------------------------

      WRITE ( 6, *)"ERROR_ISENTROPIC_VORTEX_XY gives :-"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"shift_pattern = ", shift_pattern
      WRITE ( 6, *)" "
      WRITE ( 6, *)"rhoerr_L1, rhoerr_Linf = ",
     1              rhoerr_L1, rhoerr_Linf
      WRITE ( 6, *)"engerr_L1, engerr_Linf = ",
     1              engerr_L1, engerr_Linf
      WRITE ( 6, *)" "
      WRITE ( 6, *)"vlxerr_L1, vlxerr_Linf = ",
     1              vlxerr_L1, vlxerr_Linf
      WRITE ( 6, *)"vlyerr_L1, vlyerr_Linf = ",
     1              vlyerr_L1, vlyerr_Linf
      WRITE ( 6, *)"vlzerr_L1, vlzerr_Linf = ",
     1              vlzerr_L1, vlzerr_Linf
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ERROR_ISENTROPIC_VORTEX_XY_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_DOUBLE_MACH_2D_XZ_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler, 

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up the 2d hydrodynamical double Mach reflection problem in the xz-plane
! with "iy1 = 1, ioffy = 0". 480x1x120 zones is the traditional resolution.
! The problem is catalogued in Woodward & Colella (1984) JCP.
! Set up on a [0,4]x[-1,1]x[0,1] rectangular domain. Set "gamma_euler=1.4"
! and "timestop = 0.2". The boundaries are tricky so we give the steps:

! Step I) Set left x-boundary to inflow (no. 1) with primitive
! values given by "rho1 = 8.0, vlx1 = 7.14470958, vly1 = 0,
! vlz1 = -4.125, prs1 = 116.5". These values will have to be typed
! in by hand and the conserved variables made in PAD_BOUNDARY_X.
! Therefore, at the lower x-boundary of PAD_BOUNDARY_X insert:
!      rho_tmp = 8.0
!      momx_tmp = 8.0 * 4.125 * SQRT ( 3.0)
!      momy_tmp = 0.0
!      momz_tmp = - 8.0 * 4.125
!      eng_tmp = 116.5 / ( gamma_euler - 1.0) + 0.5 * 8.0 * 8.25 ** 2


! Step II) Set right x-boundary to outflow (no. 2).

! Step III) It is inflow (no. 1) for "xc_save ( ixx) .LE. 0.16666666666667";
! reflective (no. 3) for all other values. Set up the same values as Step I
! for the inflow. These values will have to be typed
! in by hand and the conserved variables made in PAD_BOUNDARY_Z.
! The lines that enforce the boundary conditions are:

!        IF ( xc_save ( ixx) .LE. 0.16666666666667) THEN
!          bcarr_z ( ixx, iyy, 1) = 1
!        ELSE
!          bcarr_z ( ixx, iyy, 1) = 3
!        END IF

! Step IV) Upper z-boundary has to be set via call to ADD_DOUBLE_MACH_2D_XZ.
! This call should be made once in each timestep loop and the specific
! location is specified in that subroutine's notes.
! Set upper z-boundary condition to something non-interfering (like 17)
! because it is set exclusively by the ADD_DOUBLE_MACH_2D_XZ_EULER subroutine.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler, 

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, prs1, vlx1, vly1, vlz1, rho2, prs2, vlx2, vly2, vlz2,
     1     mach21, slope12, xline12, zline12, velshk, xbottom, xtop,
     1     xtopmove, x1, z1, x2, z2, dst1, dst2, dst3, dst4,
     1     rhouse, prsuse, vlxuse, vlyuse, vlzuse,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! Set up the parameters for the Mach 10 oblique shock in this problem.

      rho2 = 1.4
      prs2 = 1.0
      mach21 = 10.0

      slope12 = SQRT ( 3.0)
      xline12 = 0.16666666666667
      zline12 = 0.0

      vlx2 = 0.0
      vly2 = 0.0
      vlz2 = 0.0

! The shock speed is given by (67.06) Courant and Friedrichs SFSW :
      tempaa = mach21 * SQRT ( gamma_euler * prs2 / rho2)

      velshk = tempaa

! And (40.08) of SFSW gives:
      tempab = ( gamma_euler - 1.0) / ( gamma_euler + 1.0)

! (67.07) gives:
      prs1 = prs2 * ( ( 1.0 + tempab) * mach21**2 - tempab )

! (67.02) gives:
      rho1 = rho2 * ( prs1 + tempab * prs2) / ( prs2 + tempab * prs1)

! (67.10) gives ( we expressly use the fact that vlx2, vlz2 = 0.0) :
      vlx1 = ( 1.0 - tempab) * ( mach21 - 1.0 / mach21)
      vlx1 = vlx1 * SQRT ( gamma_euler * prs2 / rho2)

      vly1 = 0.0

      tempaa = vlx1
      vlx1 = tempaa * slope12 / SQRT ( 1.0 + slope12**2 )
      vlz1 = - tempaa * 1.0 / SQRT ( 1.0 + slope12**2 )

      xbottom = xline12 + ( zmin - zline12) / slope12
      xtop = xline12 + ( zmax - zline12) / slope12

! Note that the shock moves obliquely wrt the wall so there is a
! sine factor that needs to be divided to the speed of its motion
! to get its speed of motion wrt the wall.

      velshk = velshk * SQRT ( 1.0 + slope12**2 ) / slope12

! "xtopmove" is the location of the shock at the top z-boundary of
! the domain. At the initial time it remains constant.

      xtopmove = xtop

! ----------------------------------------------------------------------


! ------------------------------

! Now set up interior. The top of the shock will always pass through the point
! ( xtopmove, 0, zmax). Hence we key off that to decide whether a zone
! is in the unshocked or shocked region and set it accordingly.

! Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          x1 = xb_save ( ixx - 1)
          x2 = xb_save ( ixx)
          z1 = zb_save ( izz - 1)
          z2 = zb_save ( izz)

          dst1 = -slope12 * ( x1 - xtopmove) + ( z1 - zmax)
          dst2 = -slope12 * ( x2 - xtopmove) + ( z1 - zmax)
          dst3 = -slope12 * ( x1 - xtopmove) + ( z2 - zmax)
          dst4 = -slope12 * ( x2 - xtopmove) + ( z2 - zmax)

! ----------

          IF ( ( dst1 .LT. 0.0) .AND. ( dst2 .LT. 0.0) .AND.
     1         ( dst3 .LT. 0.0) .AND. ( dst4 .LT. 0.0) ) THEN

! We are squarely in the unshocked region (region 2).

          rhouse = rho2
          vlxuse = vlx2
          vlyuse = vly2
          vlzuse = vlz2
          prsuse = prs2

          ELSE IF ( ( dst1 .GT. 0.0) .AND. ( dst2 .GT. 0.0) .AND. 
     1              ( dst3 .GT. 0.0) .AND. ( dst4 .GT. 0.0) ) THEN

! We are squarely in the shocked region (region 1).

          rhouse = rho1
          vlxuse = vlx1
          vlyuse = vly1
          vlzuse = vlz1
          prsuse = prs1

          ELSE

! The zone straddles the shock so average its values.

          rhouse = 0.5 * ( rho1 + rho2)
          vlxuse = 0.5 * ( vlx1 + vlx2)
          vlyuse = 0.5 * ( vly1 + vly2)
          vlzuse = 0.5 * ( vlz1 + vlz2)
          prsuse = 0.5 * ( prs1 + prs2)

          END IF

! ----------

! Initialize the interior values.

          u ( ixx, iyy, izz, 1, 1) = rhouse
          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse
          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          tempaa = SQRT ( vlxuse**2 + vlyuse**2 + vlzuse**2 )

          tempab = SQRT ( gamma_euler * prsuse / rhouse )

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_DOUBLE_MACH_2D_XZ_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_DOUBLE_MACH_2D_XZ_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler, 

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Upper z-boundary of the double Mach reflection problem
! has to be set via call to this subroutine -- ADD_DOUBLE_MACH_2D_XZ_EULER.
! This call should be made once in each timestep loop. Best place
! to put it is after the padding of the ghost zones but before the
! call to the limiter.

! Note that the interface of this subroutine is different from that of
! INIT_DOUBLE_MACH_2D_XZ_EULER. The "dtcur" has been replaced with "time".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     rho1, prs1, vlx1, vly1, vlz1, rho2, prs2, vlx2, vly2, vlz2,
     1     mach21, slope12, xline12, zline12, velshk, xbottom, xtop,
     1     xtopmove, x1, z1, x2, z2, dst1, dst2, dst3, dst4,
     1     rhouse, prsuse, vlxuse, vlyuse, vlzuse,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! Set up the parameters for the Mach 10 oblique shock in this problem.

      rho2 = 1.4
      prs2 = 1.0
      mach21 = 10.0

      slope12 = SQRT ( 3.0)
      xline12 = 0.16666666666667
      zline12 = 0.0

      vlx2 = 0.0
      vly2 = 0.0
      vlz2 = 0.0

! The shock speed is given by (67.06) Courant and Friedrichs SFSW :
      tempaa = mach21 * SQRT ( gamma_euler * prs2 / rho2)

      velshk = tempaa

! And (40.08) of SFSW gives:
      tempab = ( gamma_euler - 1.0) / ( gamma_euler + 1.0)

! (67.07) gives:
      prs1 = prs2 * ( ( 1.0 + tempab) * mach21**2 - tempab )

! (67.02) gives:
      rho1 = rho2 * ( prs1 + tempab * prs2) / ( prs2 + tempab * prs1)

! (67.10) gives ( we expressly use the fact that vlx2, vlz2 = 0.0) :
      vlx1 = ( 1.0 - tempab) * ( mach21 - 1.0 / mach21)
      vlx1 = vlx1 * SQRT ( gamma_euler * prs2 / rho2)

      vly1 = 0.0

      tempaa = vlx1
      vlx1 = tempaa * slope12 / SQRT ( 1.0 + slope12**2 )
      vlz1 = - tempaa * 1.0 / SQRT ( 1.0 + slope12**2 )

      xbottom = xline12 + ( zmin - zline12) / slope12
      xtop = xline12 + ( zmax - zline12) / slope12

! Note that the shock moves obliquely wrt the wall so there is a
! sine factor that needs to be divided to the speed of its motion
! to get its speed of motion wrt the wall.

      velshk = velshk * SQRT ( 1.0 + slope12**2 ) / slope12

! "xtopmove" is the location of the shock at the top z-boundary of
! the domain. It moves as time progresses.

      xtopmove = xtop + velshk * time

! ----------------------------------------------------------------------

! Now set up interior. The top of the shock will always pass through the point
! ( xtopmove, 0, zmax). Hence we key off that to decide whether a zone
! is in the unshocked or shocked region and set it accordingly.

      DO izz = iz1 + 1, iz1 + ioffz  ! Notice, this is the top z-boundary.

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          x1 = xb_save ( ixx - 1)
          x2 = xb_save ( ixx)
          z1 = zb_save ( izz - 1)
          z2 = zb_save ( izz)

          dst1 = -slope12 * ( x1 - xtopmove) + ( z1 - zmax)
          dst2 = -slope12 * ( x2 - xtopmove) + ( z1 - zmax)
          dst3 = -slope12 * ( x1 - xtopmove) + ( z2 - zmax)
          dst4 = -slope12 * ( x2 - xtopmove) + ( z2 - zmax)

! ----------

          IF ( ( dst1 .LT. 0.0) .AND. ( dst2 .LT. 0.0) .AND.
     1         ( dst3 .LT. 0.0) .AND. ( dst4 .LT. 0.0) ) THEN

! We are squarely in the unshocked region (region 2).

          rhouse = rho2
          vlxuse = vlx2
          vlyuse = vly2
          vlzuse = vlz2
          prsuse = prs2

          ELSE IF ( ( dst1 .GT. 0.0) .AND. ( dst2 .GT. 0.0) .AND. 
     1              ( dst3 .GT. 0.0) .AND. ( dst4 .GT. 0.0) ) THEN

! We are squarely in the shocked region (region 1).

          rhouse = rho1
          vlxuse = vlx1
          vlyuse = vly1
          vlzuse = vlz1
          prsuse = prs1

          ELSE

! The zone straddles the shock so average its values.

          rhouse = 0.5 * ( rho1 + rho2)
          vlxuse = 0.5 * ( vlx1 + vlx2)
          vlyuse = 0.5 * ( vly1 + vly2)
          vlzuse = 0.5 * ( vlz1 + vlz2)
          prsuse = 0.5 * ( prs1 + prs2)

          END IF

! ----------

! Set ghost zone values at the top z-boundary.

          u ( ixx, iyy, izz, 1, 1) = rhouse
          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse
          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADD_DOUBLE_MACH_2D_XZ_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_FF_STEP_2D_XY_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler, 

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up the 2d hydrodynamical forward facing problem in the
! xy-plane with iz1 = 1. 240x80x1 zones is the traditional resolution.
! At the above resolution, the corner of the step is at the right upper
! corner of the zone given by "ixx = 48, iyy = 16". 
! The problem is catalogued in Woodward & Colella (1984) JCP.
! Set up on a [0,3]x[0,1]x[-1,1] rectangular domain. Set "gamma_euler=1.4"
! and "timestop = 4.0".

! Declare and set "nxstep = 48 * ( ix1 / 240); nystep = 16 * ( iy1 / 80)". 
! The zone "( nxstep, nystep, :)" is the most important zone and everything 
! keys off it.

! The boundaries are tricky so we give the steps:

! Step I) Set the x and y boundary conditions.

!      IF ( iyy .LE. nystep) x_indx_limits ( 2, iyy, izz) = nxstep
!      bcarr_x ( 1, iyy, izz) = 1  ! And remember to set inflow data.
!      bcarr_x ( 2, iyy, izz) = 2  ! Continuitive outflow at top x-boundary.
!      IF ( iyy .LE. nystep) bcarr_x ( 2, iyy, izz) = 3 ! Reflective at step.

!      IF ( ixx .GT. nxstep) y_indx_limits ( ixx, 1, izz) = nystep + 1
!      bcarr_y ( ixx, 1, izz) = 3  ! Reflective at bottom y-boundary.
!      bcarr_y ( ixx, 2, izz) = 3  ! Reflective at top y-boundary.

! Inflow at lower x-boundary is set as follows:

!      rho_tmp = 1.4
!      momx_tmp = 1.4 * 3.0
!      momy_tmp = 0.0
!      momz_tmp = 0.0
!      eng_tmp = 1.0/(gamma_euler-1.0) + 0.5 * 1.4 * 3.0**2

! Step II) Because the reflection at the top x-face of the step is most
! important, call PAD_BOUNDARY_Y before PAD_BOUNDARY_X.

! Step III) After padding boundary and before limiter, call
! ADD_FF_STEP_2D_XY_EULER.

! STEP IV) We can cover up the flow variables in the step with:
!      izz = 1
!      DO iyy = 1, nystep
!      DO ixx = nxstep + 1, ix1
!      r4gr ( ixx, iyy, izz) = r4gr ( 1, 1, izz)
!      END DO
!      END DO

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler, 

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rhouse, prsuse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u


      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! Set up the parameters for the Mach 3 flow in this problem.

      rhouse = 1.4
      prsuse = 1.0

      vlxuse = 3.0
      vlyuse = 0.0
      vlzuse = 0.0


! ------------------------------

! Now set up interior. The top of the shock will always pass through the point
! ( xtopmove, 0, zmax). Hence we key off that to decide whether a zone
! is in the unshocked or shocked region and set it accordingly.

! Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

! Initialize the interior values.

          u ( ixx, iyy, izz, 1, 1) = rhouse
          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse
          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          tempaa = SQRT ( vlxuse**2 + vlyuse**2 + vlzuse**2 )

          tempab = SQRT ( gamma_euler * prsuse / rhouse )

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_FF_STEP_2D_XY_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_FF_STEP_2D_XY_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Provides the entropy and enthalpy fixes to 6 zones at the corner of the step.
! This fix has been described in Woodward & Colella (1984) JCP.

! Place this subroutine after padding boundaries but before call to limiter.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz, nxstep, nystep

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     rho_step, prs_step, vlx_step, vly_step, vlz_step,
     1     entropy_step, enthalpy_step,
     1     rhouse, prsuse, vlxuse, vlyuse, vlzuse, td_enthalpyuse,
     1     bfxuse, bfyuse, bfzuse, tempaa

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! Set up the free stream specific enthalpy and entropy for the first zone
! that abuts the step in the x-direction.

      nxstep = 48 * ( ix1 / 240)
      nystep = 16 * ( iy1 / 80)
      izz = 1

      rho_step = AMAX1 ( u ( nxstep, nystep, izz, 1, 1), rhofloor)

      vlx_step = u ( nxstep, nystep, izz, 2, 1) / rho_step
      vly_step = u ( nxstep, nystep, izz, 3, 1) / rho_step
      vlz_step = u ( nxstep, nystep, izz, 4, 1) / rho_step

      prs_step = (gamma_euler - 1.0) * ( u ( nxstep, nystep, izz, 5, 1)
     1 - 0.5 * rho_step * ( vlx_step**2 + vly_step**2 + vlz_step**2) )
      prs_step = AMAX1 ( prs_step, prsfloor)

      entropy_step = prs_step / rho_step**gamma_euler
      enthalpy_step = ( gamma_euler * prs_step)
     1              / ( rho_step * ( gamma_euler - 1.0))
     1 + 0.5 * ( vlx_step**2 + vly_step**2 + vlz_step**2)

! ----------------------------------------------------------------------

! Entropy fix, zone 1.

      ixx = nxstep + 1
      iyy = nystep + 1
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma_euler - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma_euler)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma_euler * prsuse)
     1               / ( rhouse * ( gamma_euler - 1.0))

      tempaa = AMAX1 ( enthalpy_step - td_enthalpyuse, smallnum)
      tempaa = 2.0 * tempaa
     1       / ( vlxuse**2 + vlyuse**2 + vlzuse**2 + 2.0 * smallnum)
      tempaa = SQRT ( tempaa)

      vlxuse = tempaa * vlxuse
      vlyuse = tempaa * vlyuse
      vlzuse = tempaa * vlzuse

      u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
      u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
      u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma_euler - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

! Entropy fix, zone 2.

      ixx = nxstep + 2
      iyy = nystep + 1
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma_euler - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma_euler)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma_euler * prsuse)
     1               / ( rhouse * ( gamma_euler - 1.0))

      tempaa = AMAX1 ( enthalpy_step - td_enthalpyuse, smallnum)
      tempaa = 2.0 * tempaa
     1       / ( vlxuse**2 + vlyuse**2 + vlzuse**2 + 2.0 * smallnum)
      tempaa = SQRT ( tempaa)

      vlxuse = tempaa * vlxuse
      vlyuse = tempaa * vlyuse
      vlzuse = tempaa * vlzuse

      u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
      u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
      u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma_euler - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

! Entropy fix, zone 3.

      ixx = nxstep + 3
      iyy = nystep + 1
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma_euler - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma_euler)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma_euler * prsuse)
     1    / ( rhouse * ( gamma_euler - 1.0))

      tempaa = AMAX1 ( enthalpy_step - td_enthalpyuse, smallnum)
      tempaa = 2.0 * tempaa
     1       / ( vlxuse**2 + vlyuse**2 + vlzuse**2 + 2.0 * smallnum)
      tempaa = SQRT ( tempaa)

      vlxuse = tempaa * vlxuse
      vlyuse = tempaa * vlyuse
      vlzuse = tempaa * vlzuse

      u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
      u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
      u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma_euler - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

! Entropy fix, zone 4.

      ixx = nxstep + 4
      iyy = nystep + 1
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma_euler - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma_euler)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma_euler * prsuse)
     1    / ( rhouse * ( gamma_euler - 1.0))

      tempaa = AMAX1 ( enthalpy_step - td_enthalpyuse, smallnum)
      tempaa = 2.0 * tempaa
     1       / ( vlxuse**2 + vlyuse**2 + vlzuse**2 + 2.0 * smallnum)
      tempaa = SQRT ( tempaa)

      vlxuse = tempaa * vlxuse
      vlyuse = tempaa * vlyuse
      vlzuse = tempaa * vlzuse

      u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
      u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
      u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma_euler - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

! Entropy fix, zone 5.

      ixx = nxstep + 1
      iyy = nystep + 2
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma_euler - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma_euler)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma_euler * prsuse)
     1               / ( rhouse * ( gamma_euler - 1.0))

      tempaa = AMAX1 ( enthalpy_step - td_enthalpyuse, smallnum)
      tempaa = 2.0 * tempaa
     1       / ( vlxuse**2 + vlyuse**2 + vlzuse**2 + 2.0 * smallnum)
      tempaa = SQRT ( tempaa)

      vlxuse = tempaa * vlxuse
      vlyuse = tempaa * vlyuse
      vlzuse = tempaa * vlzuse

      u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
      u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
      u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma_euler - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

! Entropy fix, zone 6.

      ixx = nxstep + 2
      iyy = nystep + 2
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma_euler - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma_euler)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma_euler * prsuse)
     1               / ( rhouse * ( gamma_euler - 1.0))

      tempaa = AMAX1 ( enthalpy_step - td_enthalpyuse, smallnum)
      tempaa = 2.0 * tempaa
     1       / ( vlxuse**2 + vlyuse**2 + vlzuse**2 + 2.0 * smallnum)
      tempaa = SQRT ( tempaa)

      vlxuse = tempaa * vlxuse
      vlyuse = tempaa * vlyuse
      vlzuse = tempaa * vlzuse

      u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
      u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
      u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma_euler - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADD_FF_STEP_2D_XY_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_TWO_D_RIEM_PROBLEM_XY_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds two variants of two-dimensional Riemann problems drawn from
! Brio, Zakharian & Webb (2001) JCP. The problems are set up in the XY
! plane with 400x400x1 zones spanning [-1,1]x[-1,1]x[-1,1].
! Set "gamma_euler = 1.4" and keep the boundaries at continuitive (no. 2).
! Such 2d Riemann problems consist of four constant states in the 4 quadrants.
! The first problem has two initial weak shocks and two slip lines and
! is stopped at "timestop = 0.52". The second problem produces double
! Mach reflection and a shock propagating at an angle to the grid lines
! and is stopped at "timestop = 0.6".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, prs1, vlx1, vly1, rho2, prs2, vlx2, vly2,
     1     rho3, prs3, vlx3, vly3, rho4, prs4, vlx4, vly4,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

! Numbers 1, 2, 3 & 4 in the variable names pertain to quadrants I, II, III & IV.

#if ( 1 == 1)
! First 2d Riemann problem.

      rho1 = 0.5313
      prs1 = 0.4
      vlx1 = 0.0
      vly1 = 0.0

      rho2 = 1.0
      prs2 = 1.0
      vlx2 = 0.7276
      vly2 = 0.0

      rho3 = 0.8
      prs3 = 1.0
      vlx3 = 0.0
      vly3 = 0.0

      rho4 = 1.0
      prs4 = 1.0
      vlx4 = 0.0
      vly4 = 0.7276

#else
! Second 2d Riemann problem.

      rho1 = 1.5
      prs1 = 1.5
      vlx1 = 0.0
      vly1 = 0.0

      rho2 = 0.5323
      prs2 = 0.3
      vlx2 = 1.206
      vly2 = 0.0

      rho3 = 0.1379
      prs3 = 0.029
      vlx3 = 1.206
      vly3 = 1.206

      rho4 = 0.5323
      prs4 = 0.3
      vlx4 = 0.0
      vly4 = 1.206

#endif

! ----------------------------------------------------------------------


! Set up interior. The top of the shock will always pass through the point
! ( xtopmove, 0, zmax). Hence we key off that to decide whether a zone
! is in the unshocked or shocked region and set it accordingly.

! Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          IF ( ( xc_save ( ixx) .GE. 0.0) .AND.
     1         ( yc_save ( iyy) .GE. 0.0) ) THEN
! Quadrant I

            rhouse = rho1
            prsuse = prs1
            vlxuse = vlx1
            vlyuse = vly1
            vlzuse = 0.0

          ELSE IF ( ( xc_save ( ixx) .LT. 0.0) .AND.
     1              ( yc_save ( iyy) .GE. 0.0) ) THEN
! Quadrant II

            rhouse = rho2
            prsuse = prs2
            vlxuse = vlx2
            vlyuse = vly2
            vlzuse = 0.0

          ELSE IF ( ( xc_save ( ixx) .LT. 0.0) .AND.
     1              ( yc_save ( iyy) .LT. 0.0) ) THEN
! Quadrant III

            rhouse = rho3
            prsuse = prs3
            vlxuse = vlx3
            vlyuse = vly3
            vlzuse = 0.0

          ELSE
! Quadrant IV

            rhouse = rho4
            prsuse = prs4
            vlxuse = vlx4
            vlyuse = vly4
            vlzuse = 0.0

          END IF

! ----------

! Initialize the interior values.

          u ( ixx, iyy, izz, 1, 1) = rhouse
          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse
          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          tempaa = SQRT ( vlxuse**2 + vlyuse**2 + vlzuse**2 )

          tempab = SQRT ( gamma_euler * prsuse / rhouse )

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_TWO_D_RIEM_PROBLEM_XY_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_SHKTUBE_1D_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, 

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 1d Euler shocktube problem.
! Problem can be done in x, y or z-direction but only one dimension
! can be dynamically active.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        ixx, iyy, izz, iproblem, nzones

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,

     1     tstop, tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) :: u

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

      REAL v ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1         NFIELD_EULER + NFLUID_EULER)

! ----------------------------------------------------------------------

! Pick a 1d hydrodynamical test problem.
! "gamma_euler" will be appropriately reset,
! but "timestop" has to be set by hand. Note too that this suite of
! Riemann problems are assumed to be set up on a unit interval.

      iproblem = 1

      nzones = MAX0 ( ix1, iy1, iz1)

      IF ( ( ioffx * ioffy .GT. 0) .OR. ( ioffy * ioffz .GT. 0) .OR.
     1     ( ioffx * ioffz .GT. 0) ) THEN
      WRITE ( 6, *)"Stopping in INIT_SHKTUBE_1D:",
     1             " must be doing 1d problem"
      STOP
      END IF

      CALL SET_RP_EULER ( iproblem, nzones, gamma_euler, tstop, v)

      WRITE ( 6, *)"Chose iproblem = ", iproblem," on unit interval"
      WRITE ( 6, *)"Suggested stopping time for this problem = ", tstop

! ----------------------------------------------------------------------

! Now set up interior. Depending on which dimension is active, the
! problem is cyclically rotated in a suitable way. Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          IF ( ix1 .GT. 1) THEN

          rhouse = v ( ixx, 1)

          vlxuse = v ( ixx, 2)
          vlyuse = v ( ixx, 3)
          vlzuse = v ( ixx, 4)

          prsuse = v ( ixx, 5)

          ELSE IF ( iy1 .GT. 1) THEN

          rhouse = v ( iyy, 1)

          vlxuse = v ( iyy, 4)
          vlyuse = v ( iyy, 2)
          vlzuse = v ( iyy, 3)

          prsuse = v ( iyy, 5)

          ELSE

          rhouse = v ( izz, 1)

          vlxuse = v ( izz, 3)
          vlyuse = v ( izz, 4)
          vlzuse = v ( izz, 2)

          prsuse = v ( izz, 5)

          END IF

          enguse = prsuse / (gamma_euler - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          u ( ixx, iyy, izz, 1, 1) = rhouse
          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse
          u ( ixx, iyy, izz, 5, 1) = enguse

! ----------

          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2

          tempab = 
     1            gamma_euler * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_SHKTUBE_1D_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE SET_RP_EULER( iproblem, nzones, gamma_euler, tstop, v)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iproblem, nzones, izonecenter, ix

      REAL gamma_euler, tstop,
     1     v ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1         NFIELD_EULER + NFLUID_EULER),

     1 rhore_l, prsre_l, vlxre_l, vlyre_l, vlzre_l, 
     1 rhore_r, prsre_r, vlxre_r, vlyre_r, vlzre_r

! ----------------------------------------------------------------------

      IF ( NFIELD_EULER + NFLUID_EULER .ne. 5) THEN
         WRITE(6,*) "Stopping in SET_RP. Check v dimensioning."
         STOP
      END IF
      IF (iproblem.lt.1 .or. iproblem .gt. 5) THEN
         WRITE (6,*) "Stopping in SET_RP.  iproblem = ", iproblem
         STOP
      END IF

*****************
*
      IF ( iproblem .EQ. 1) THEN
*
* Two hydrodynamical rarefaction fans going off to either side.
*
      gamma_euler = 1.666
*
      izonecenter = nzones / 2
*
      tstop = 0.08
*
      rhore_l = 1.0
      prsre_l = 2.0
*
      vlxre_l = -2.5
      vlyre_l = 0.0
      vlzre_l = 0.0
*
*
      rhore_r = 1.5
      prsre_r = 4.0
*
      vlxre_r = 2.5
      vlyre_r = 0.0
      vlzre_r = 0.0
*
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 2) THEN
*
* Noh-like shock tube with colliding streams.
*
      gamma_euler = 1.666
*
      izonecenter = nzones / 2
*
      tstop = 0.4
*
      rhore_l = 1.0
      prsre_l = 2.0e-1
*
      vlxre_l = 2.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      rhore_r = 1.5
      prsre_r = 2.0e-1
*
      vlxre_r = -2.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 3) THEN
*
* WC shock tube
*
      gamma_euler = 1.4
*
      izonecenter = nzones / 2
*
      tstop = 0.01
*
      rhore_l = 1.0
      prsre_l = 1000.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      rhore_r = 1.0
      prsre_r = 1.0e-2
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 4) THEN
*
* Lax shock tube
*
      gamma_euler = 1.4
*
      izonecenter = nzones / 2
*
      tstop = 0.13
*
      rhore_l = 0.445
      prsre_l = 3.528
*
      vlxre_l = 0.698
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      rhore_r = 0.5
      prsre_r = 0.571
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 5) THEN
*
* Sod shock tube
*
      gamma_euler = 1.4
*
      izonecenter = nzones / 2
*
      tstop = 0.2
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      rhore_r = 0.125
      prsre_r = 0.1
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      END IF
*
*****************
*

*
*****************

! Now set up all the primitive variables.

      DO ix = - 3, nzones + 4

      IF ( ix .LE. izonecenter) THEN

      v ( ix, 1) = rhore_l
      v ( ix, 2) = vlxre_l
      v ( ix, 3) = vlyre_l
      v ( ix, 4) = vlzre_l
      v ( ix, 5) = prsre_l

      ELSE

      v ( ix, 1) = rhore_r
      v ( ix, 2) = vlxre_r
      v ( ix, 3) = vlyre_r
      v ( ix, 4) = vlzre_r
      v ( ix, 5) = prsre_r

      END IF

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma_euler, tstop = ", gamma_euler, tstop 
      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"rhore_l, prsre_l, vlxre_l, vlyre_l, vlzre_l"
      WRITE ( 6, *)
     1 rhore_l, prsre_l, vlxre_l, vlyre_l, vlzre_l
      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"rhore_r, prsre_r, vlxre_r, vlyre_r, vlzre_r"
      WRITE ( 6, *)
     1 rhore_r, prsre_r, vlxre_r, vlyre_r, vlzre_r

      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE SET_RP_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------





