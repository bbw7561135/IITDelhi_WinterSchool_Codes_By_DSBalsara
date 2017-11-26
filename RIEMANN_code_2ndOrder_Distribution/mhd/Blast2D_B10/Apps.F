! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_BLAST_3D

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 3d MHD blast problem. The problem is catalogued in
! Balsara et al. (2009) JCP.

! Can be collapsed to a 2d problem by switching off one of the directions.
! When so collapsed, it yields the 2d MHD blast problem of Balsara & Spicer
! (1999) JCP. To obtain that problem on a 200x200x1 zone mesh on the
! interval [-.5,.5]x[-.5,.5]x[-.5,.5], set "gamma = 1.4", "rho1/2 = 1.0",
! "prs1 = 1000.0, prs2 = 0.1", "bfx1/2 = 100.0", "bfy1/2 = 0.0",
! "bfz1/2 = 0.0", "radius = 0.1". Run till "timestop = 0.01".
! "bfx1/2 = 100.0" does not require pressure-positivity fix, if it is
! run with "MULTID_LIMITING 1", "cfl_coef_sdsp/_advct = 0.7/0.45",
! "blend_multid = 0.7", "blend_mdhll_to_mdhllc = 0.2", "lapidus_coef = 0.05"
! in this code. A minimum of "blend_mdhll_to_mdhllc = 0.2" has to be supplied
! to all Riemann solvers in "APPLY_MULTID_RS_?_EDGE" routines.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,
     1     vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,
     1     radius, xcenter, ycenter, zcenter,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

! When using "PREDICT_WITH_RS 1" and "EVAL_EFIELDS_AT_EDGES 1" we can run
! with magnetic fields as large as 40.0 in 3d for this problem without
! encountering pressure positivity issues.

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

      bfx1 = 10.0 / SQRT ( 2.0)
      bfx2 = bfx1

      bfy1 = 10.0 / SQRT ( 2.0)
      bfy2 = bfy1

      bfz1 = 0.0
      bfz2 = bfz1

      radius = 0.1

! Make the center of the domain so that the explosion can be set up there.

      xcenter = 0.5 * ( base_grid_xmin + base_grid_xmax)
      ycenter = 0.5 * ( base_grid_ymin + base_grid_ymax)
      zcenter = 0.5 * ( base_grid_zmin + base_grid_zmax)

! ----------------------------------------------------------------------

! Make magnetic fields in active directions.

      IF ( ix1 .GT. 1) THEN

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1
          bfxbdy ( ixx, iyy, izz, 1) = bfx2
          END DO
        END DO
      END DO

      ELSE

          bfx1 = 0.0
          bfx2 = 0.0
          bfxbdy ( :, :, :, 1) = 0.0

      END IF


      IF ( iy1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 0, iy1
          DO izz = 1, iz1
          bfybdy ( ixx, iyy, izz, 1) = bfy2
          END DO
        END DO
      END DO

      ELSE

          bfy1 = 0.0
          bfy2 = 0.0
          bfybdy ( :, :, :, 1) = 0.0

      END IF

      IF ( iz1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 0, iz1
          bfzbdy ( ixx, iyy, izz, 1) = bfz2
          END DO
        END DO
      END DO

      ELSE

          bfz1 = 0.0
          bfz2 = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

      END IF

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

          u ( ixx, iyy, izz, 5, 1) = prs1 / (gamma - 1.0)
     1         + 0.5 * rho1 * ( vlx1**2 + vly1**2 + vlz1**2 )
     1         + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx1
          u ( ixx, iyy, izz, 7, 1) = bfy1
          u ( ixx, iyy, izz, 8, 1) = bfz1

          tempab = SQRT ( vlx1**2 + vly1**2 + vlz1**2)
          tempac = gamma * prs1
     1           + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 4.0 * pi)
          tempac = SQRT ( tempac / rho1) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef_advct * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) = rho1
#endif

          ELSE

          u ( ixx, iyy, izz, 1, 1) = rho2

          u ( ixx, iyy, izz, 2, 1) = rho2 * vlx2
          u ( ixx, iyy, izz, 3, 1) = rho2 * vly2
          u ( ixx, iyy, izz, 4, 1) = rho2 * vlz2

          u ( ixx, iyy, izz, 5, 1) = prs2 / (gamma - 1.0)
     1         + 0.5 * rho2 * ( vlx2**2 + vly2**2 + vlz2**2 )
     1         + ( bfx2**2 + bfy2**2 + bfz2**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx2
          u ( ixx, iyy, izz, 7, 1) = bfy2
          u ( ixx, iyy, izz, 8, 1) = bfz2

          tempab = SQRT ( vlx2**2 + vly2**2 + vlz2**2)
          tempac = gamma * prs2
     1           + ( bfx2**2 + bfy2**2 + bfz2**2) / ( 4.0 * pi)
          tempac = SQRT ( tempac / rho2) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef_advct * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) = 0.0
#endif

          END IF

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_BLAST_3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ROTOR_2D_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 2d MHD rotor problem in the xy-plane with iz1 = 1. The problem
! is catalogued in Balsara & Spicer (1999) JCP and Balsara (2004) ApJS.
! Set up on a unit square with continuitive (no. 2) boundary conditions.
! 200x200 zones is ample. Stop at "timestop = 0.29". Use "gamma = 5.0/3.0".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, rho2, prs1, bfx1, bfy1, bfz1, 
     1     radius1, radius2, vrot1, xcenter, ycenter, zcenter,
     1     costheta, sintheta, radius, tempvrot,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

! The taper is designed so that it spreads over about 6 zones on
! whichever single_grid level one uses.

      tempaa = AMIN1 ( (xmax-xmin)/ix1, (ymax-ymin)/iy1,
     1                 (zmax-zmin)/iz1)

      rho1 = 1.0
      rho2 = 10.0

      prs1 = 0.5

      bfx1 = 2.5
      bfy1 = 0.0
      bfz1 = 0.0

      vrot1 = 1.0

      radius1 = 0.1
      radius2 = radius1 + 6.0 * tempaa

! Make the center of the domain so that the rotor can be set up there.

      xcenter = 0.5 * ( base_grid_xmin + base_grid_xmax)
      ycenter = 0.5 * ( base_grid_ymin + base_grid_ymax)
      zcenter = 0.5 * ( base_grid_zmin + base_grid_zmax)

! ----------------------------------------------------------------------

! Make magnetic fields in active directions.

      IF ( ix1 .GT. 1) THEN

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1
          bfxbdy ( ixx, iyy, izz, 1) = bfx1
          END DO
        END DO
      END DO

      ELSE

          bfx1 = 0.0
          bfxbdy ( :, :, :, 1) = 0.0

      END IF


      IF ( iy1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 0, iy1
          DO izz = 1, iz1
          bfybdy ( ixx, iyy, izz, 1) = bfy1
          END DO
        END DO
      END DO

      ELSE

          bfy1 = 0.0
          bfybdy ( :, :, :, 1) = 0.0

      END IF

      IF ( iz1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 0, iz1
          bfzbdy ( ixx, iyy, izz, 1) = bfz1
          END DO
        END DO
      END DO

      ELSE

          bfz1 = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

      END IF

! ------------------------------

! Now set up interior. Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          radius = ( xc_save ( ixx) - xcenter)**2
     1           + ( yc_save ( iyy) - ycenter)**2

          radius = SQRT ( radius)
          costheta = ( xc_save ( ixx) - xcenter) / radius
          sintheta = ( yc_save ( iyy) - ycenter) / radius

          IF ( radius .LE. radius1) THEN

            tempaa = 1.0
            tempvrot = vrot1 * radius / radius1

          ELSE IF ( ( radius .LE. radius2) .AND.
     1              ( radius2 .GT. radius1) ) THEN

            tempaa = ( radius2 - radius) / ( radius2 - radius1)
            tempvrot = vrot1 * tempaa

          ELSE

            tempaa = 0.0
            tempvrot = 0.0

          END IF

          u ( ixx, iyy, izz, 1, 1) = rho1 * ( 1.0 - tempaa)
     1                       + rho2 * tempaa

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 1, 1)
     1  * ( - tempvrot * sintheta)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 1, 1)
     1  * ( tempvrot * costheta)

          u ( ixx, iyy, izz, 4, 1) = 0.0

          u ( ixx, iyy, izz, 5, 1) = prs1 / (gamma - 1.0)
     1  + 0.5 * u ( ixx, iyy, izz, 1, 1)
     1      * ( ( - tempvrot * sintheta)**2
     1        + ( tempvrot * costheta)**2)
     1         + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx1
          u ( ixx, iyy, izz, 7, 1) = bfy1
          u ( ixx, iyy, izz, 8, 1) = bfz1

          tempab = SQRT ( ( - tempvrot * sintheta)**2
     1                  + ( tempvrot * costheta)**2 )
          tempac = gamma * prs1
     1           + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 4.0 * pi)
          tempac = SQRT ( tempac / u ( ixx, iyy, izz, 1, 1)) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef_advct * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ROTOR_2D_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ROTOR_2D_XZ

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 2d MHD rotor problem in the xz-plane with iy1 = 1. The problem
! is catalogued in Balsara & Spicer (1999) JCP and Balsara (2004) ApJS.
! Set up on a unit square with continuitive (no. 2) boundary conditions.
! 200x200 zones is ample. Stop at "timestop = 0.29". Use "gamma = 5.0/3.0".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, rho2, prs1, bfx1, bfy1, bfz1, 
     1     radius1, radius2, vrot1, xcenter, ycenter, zcenter,
     1     costheta, sintheta, radius, tempvrot,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

! The taper is designed so that it spreads over about 6 zones on
! whichever single_grid level one uses.

      tempaa = AMIN1 ( (xmax-xmin)/ix1, (ymax-ymin)/iy1,
     1                 (zmax-zmin)/iz1)

      rho1 = 1.0
      rho2 = 10.0

      prs1 = 0.5

      bfx1 = 2.5
      bfy1 = 0.0
      bfz1 = 0.0

      vrot1 = 1.0

      radius1 = 0.1
      radius2 = radius1 + 6.0 * tempaa

! Make the center of the domain so that the rotor can be set up there.

      xcenter = 0.5 * ( base_grid_xmin + base_grid_xmax)
      ycenter = 0.5 * ( base_grid_ymin + base_grid_ymax)
      zcenter = 0.5 * ( base_grid_zmin + base_grid_zmax)

! ----------------------------------------------------------------------

! Make magnetic fields in active directions.

      IF ( ix1 .GT. 1) THEN

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1
          bfxbdy ( ixx, iyy, izz, 1) = bfx1
          END DO
        END DO
      END DO

      ELSE

          bfx1 = 0.0
          bfxbdy ( :, :, :, 1) = 0.0

      END IF


      IF ( iy1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 0, iy1
          DO izz = 1, iz1
          bfybdy ( ixx, iyy, izz, 1) = bfy1
          END DO
        END DO
      END DO

      ELSE

          bfy1 = 0.0
          bfybdy ( :, :, :, 1) = 0.0

      END IF

      IF ( iz1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 0, iz1
          bfzbdy ( ixx, iyy, izz, 1) = bfz1
          END DO
        END DO
      END DO

      ELSE

          bfz1 = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

      END IF

! ------------------------------

! Now set up interior. Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          radius = ( xc_save ( ixx) - xcenter)**2
     1           + ( zc_save ( izz) - zcenter)**2

          radius = SQRT ( radius)
          costheta = ( xc_save ( ixx) - xcenter) / radius
          sintheta = ( zc_save ( izz) - zcenter) / radius

          IF ( radius .LE. radius1) THEN

            tempaa = 1.0
            tempvrot = vrot1 * radius / radius1

          ELSE IF ( ( radius .LE. radius2) .AND.
     1              ( radius2 .GT. radius1) ) THEN

            tempaa = ( radius2 - radius) / ( radius2 - radius1)
            tempvrot = vrot1 * tempaa

          ELSE

            tempaa = 0.0
            tempvrot = 0.0

          END IF

          u ( ixx, iyy, izz, 1, 1) = rho1 * ( 1.0 - tempaa)
     1                       + rho2 * tempaa

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 1, 1)
     1  * ( - tempvrot * sintheta)

          u ( ixx, iyy, izz, 3, 1) = 0.0

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 1, 1)
     1  * ( tempvrot * costheta)

          u ( ixx, iyy, izz, 5, 1) = prs1 / (gamma - 1.0)
     1  + 0.5 * u ( ixx, iyy, izz, 1, 1)
     1      * ( ( - tempvrot * sintheta)**2
     1        + ( tempvrot * costheta)**2)
     1         + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx1
          u ( ixx, iyy, izz, 7, 1) = bfy1
          u ( ixx, iyy, izz, 8, 1) = bfz1

          tempab = SQRT ( ( - tempvrot * sintheta)**2
     1                  + ( tempvrot * costheta)**2 )
          tempac = gamma * prs1
     1           + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 4.0 * pi)
          tempac = SQRT ( tempac / u ( ixx, iyy, izz, 1, 1)) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef_advct * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ROTOR_2D_XZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ROTOR_2D_YZ

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 2d MHD rotor problem in the yz-plane with ix1 = 1. The problem
! is catalogued in Balsara & Spicer (1999) JCP and Balsara (2004) ApJS.
! Set up on a unit square with continuitive (no. 2) boundary conditions.
! 200x200 zones is ample. Stop at "timestop = 0.29". Use "gamma = 5.0/3.0".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, rho2, prs1, bfx1, bfy1, bfz1, 
     1     radius1, radius2, vrot1, xcenter, ycenter, zcenter,
     1     costheta, sintheta, radius, tempvrot,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

! The taper is designed so that it spreads over about 6 zones on
! whichever single_grid level one uses.

      tempaa = AMIN1 ( (xmax-xmin)/ix1, (ymax-ymin)/iy1,
     1                 (zmax-zmin)/iz1)

      rho1 = 1.0
      rho2 = 10.0

      prs1 = 0.5

      bfx1 = 0.0

      bfy1 = 2.5

      bfz1 = 0.0

      vrot1 = 1.0

      radius1 = 0.1
      radius2 = radius1 + 6.0 * tempaa

! Make the center of the domain so that the rotor can be set up there.

      xcenter = 0.5 * ( base_grid_xmin + base_grid_xmax)
      ycenter = 0.5 * ( base_grid_ymin + base_grid_ymax)
      zcenter = 0.5 * ( base_grid_zmin + base_grid_zmax)

! ----------------------------------------------------------------------

! Make magnetic fields in active directions.

      IF ( ix1 .GT. 1) THEN

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1
          bfxbdy ( ixx, iyy, izz, 1) = bfx1
          END DO
        END DO
      END DO

      ELSE

          bfx1 = 0.0
          bfxbdy ( :, :, :, 1) = 0.0

      END IF


      IF ( iy1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 0, iy1
          DO izz = 1, iz1
          bfybdy ( ixx, iyy, izz, 1) = bfy1
          END DO
        END DO
      END DO

      ELSE

          bfy1 = 0.0
          bfybdy ( :, :, :, 1) = 0.0

      END IF

      IF ( iz1 .GT. 1) THEN

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 0, iz1
          bfzbdy ( ixx, iyy, izz, 1) = bfz1
          END DO
        END DO
      END DO

      ELSE

          bfz1 = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

      END IF

! ------------------------------

! Now set up interior. Also set up the timestep.

      dtcur = 1.0e32

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          radius = ( yc_save ( iyy) - ycenter)**2
     1           + ( zc_save ( izz) - zcenter)**2

          radius = SQRT ( radius)
          costheta = ( yc_save ( iyy) - ycenter) / radius
          sintheta = ( zc_save ( izz) - zcenter) / radius

          IF ( radius .LE. radius1) THEN

            tempaa = 1.0
            tempvrot = vrot1 * radius / radius1

          ELSE IF ( ( radius .LE. radius2) .AND.
     1              ( radius2 .GT. radius1) ) THEN

            tempaa = ( radius2 - radius) / ( radius2 - radius1)
            tempvrot = vrot1 * tempaa

          ELSE

            tempaa = 0.0
            tempvrot = 0.0

          END IF

          u ( ixx, iyy, izz, 1, 1) = rho1 * ( 1.0 - tempaa)
     1                       + rho2 * tempaa

          u ( ixx, iyy, izz, 2, 1) = 0.0

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 1, 1)
     1  * ( - tempvrot * sintheta)

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 1, 1)
     1  * ( tempvrot * costheta)

          u ( ixx, iyy, izz, 5, 1) = prs1 / (gamma - 1.0)
     1  + 0.5 * u ( ixx, iyy, izz, 1, 1)
     1      * ( ( - tempvrot * sintheta)**2
     1        + ( tempvrot * costheta)**2)
     1         + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx1
          u ( ixx, iyy, izz, 7, 1) = bfy1
          u ( ixx, iyy, izz, 8, 1) = bfz1

          tempab = SQRT ( ( - tempvrot * sintheta)**2
     1                  + ( tempvrot * costheta)**2 )
          tempac = gamma * prs1
     1           + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 4.0 * pi)
          tempac = SQRT ( tempac / u ( ixx, iyy, izz, 1, 1)) + tempab
          tempad = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))
          tempab = cfl_coef_advct * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ROTOR_2D_YZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ALFVEN_2D_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up an obliquely propagating torsional Alfven wave in the 2d xy-plane
! with iz1 = 1. The problem is catalogued in Balsara (2004) ApJS.
! Use a domain of [-3,3]x[-3,3]x[-3,3] on a 120x120x1 zone mesh for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 1.4" suggested. Suggested stopping time is
! "timestop = 129.0". Measure the decay in the amplitude of the
! Alfven wave with time. This is most easily done by plotting out the
! evolution of the maximum z-velocity or the maximum z-component of
! magnetic field with time.

! When used with multi-dimensional Riemann solvers, realize that this isn't
! a "TRUE_2D 1" problem.

! Poorer interpolation schemes along with lower-grade Riemann solvers
! cause very rapid decay in the z-velocity or z-magnetic field. Very good
! interpolation schemes along with high-grade Riemann solvers
! reduce this decay by several orders of magnitude!

! The final amplitude at 129 time units for the magnetic field is given below:


! TVD-MC reconstruction with HLL RS: 6.5386E-3
! TVD-MC reconstruction with HLLC RS: 7.7403E-2
! TVD-MC reconstruction with HLLD RS: 0.21814
! TVD-MC reconstruction with ROE RS: 0.21825
! TVD-MC reconstruction with Multi-d HLLC RS: 0.08945


! WENO reconstruction with HLL RS: 2.1753E-2

! WENO reconstruction with HLLC RS: 0.12504

! WENO reconstruction with HLLD RS: 0.3918

! WENO reconstruction with ROE RS: 0.3918

! WENO reconstruction with Multi-d HLLC RS: 0.1913

! ----------

! WENO reconstruction with 1d HLL RS (CFL 0.4, ADER predictor,
!       PREDICT_WITH_RS 0, EVAL_EFIELDS_AT_EDGES 0, MULTID_LIMITING 0) : 0.02169
!                           The GS variant of double dissipation gives : 0.02168 

! WENO reconstruction with 1d HLLC RS (CFL 0.4, ADER predictor,
!       PREDICT_WITH_RS 0, EVAL_EFIELDS_AT_EDGES 0, MULTID_LIMITING 0) : 0.1235
!                           The GS variant of double dissipation gives : 0.1235

! WENO reconstruction with 1d HLLD RS (CFL 0.4, ADER predictor,
!       PREDICT_WITH_RS 0, EVAL_EFIELDS_AT_EDGES 0, MULTID_LIMITING 0) : 0.3861
!                           The GS variant of double dissipation gives : 0.3861

! ----------

! WENO reconstruction with 1d HLL RS (CFL 0.4, ADER predictor,
!       PREDICT_WITH_RS 0, EVAL_EFIELDS_AT_EDGES 1, MULTID_LIMITING 0) : 0.02951
!                           The GS variant of double dissipation gives :

! WENO reconstruction with 1d HLLC RS (CFL 0.4, ADER predictor,
!       PREDICT_WITH_RS 0, EVAL_EFIELDS_AT_EDGES 1, MULTID_LIMITING 0) : 0.1171
!                           The GS variant of double dissipation gives :

! WENO reconstruction with 1d HLLD RS (CFL 0.4, ADER predictor,
!       PREDICT_WITH_RS 0, EVAL_EFIELDS_AT_EDGES 1, MULTID_LIMITING 0) : 0.3495
!                           The GS variant of double dissipation gives :

! ----------

! WENO reconstruction with Multi-d HLL RS (CFL 0.4, ADER predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 1.0, use_blend_HLL = 1.0,
!                                          multid_wavemodel_num = 0): 0.03635

! WENO reconstruction with Multi-d HLLC RS (CFL 0.4, ADER predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 0.25,
!                                          multid_wavemodel_num = 0):

! WENO reconstruction with Multi-d HLLD RS (CFL 0.4, ADER predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 0.25,
!                                          multid_wavemodel_num = 0):

! ----------

! WENO reconstruction with Multi-d HLL RS (CFL 0.7, multid RS predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 1.0, use_blend_HLL = 1.0,
!                                          multid_wavemodel_num = 0): 0.03423

! WENO reconstruction with Multi-d HLLC RS (CFL 0.7, multid RS predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 0.25,
!                                          multid_wavemodel_num = 0): 

! WENO reconstruction with Multi-d HLLD RS (CFL 0.7, multid RS predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 0.25,
!                                          multid_wavemodel_num = 0): 

! ----------

! WENO reconstruction with Multi-d HLL RS (CFL 0.9, multid RS predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 1.0, use_blend_HLL = 1.0,
!                                          multid_wavemodel_num = 2): unstable

! WENO reconstruction with Multi-d HLLC+HLL RS (CFL 0.9, multid RS predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 0.0,
!                                          multid_wavemodel_num = 2): unstable

! WENO reconstruction with Multi-d HLLD+HLL RS (CFL 0.9, multid RS predictor,
!       blend_multid_flux = 0.33, blend_mdhll_to_mdhllc = 0.0,
!                                          multid_wavemodel_num = 2): unstable


! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     vlx_tmp, vly_tmp, vlz_tmp, bfx_tmp, bfy_tmp, bfz_tmp, 
     1     tempaa, tempab, tempac,
     1     rho0, prs0, vln, bfn, vperp0, tempk, vax, tempphase,
     1     tempsin, tempcos, tempvy1, tempvz1, tempby1, tempbz1, time,
     1     rr, nx, ny, x, y, z, kx, ky, phase, x1, x2, y1, y2

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

! Stop this MHD torsional Alfven wave problem at a time of "timestop = 129.0".

! Set up the direction of propagation.

      rr = base_grid_xmax - base_grid_xmin
      nx = 1.0 / SQRT ( rr**2 + 1.0)
      ny = rr / SQRT ( rr**2 + 1.0)

! Set up density, pressure, magnetic field and velocity along (nx,ny).

      rho0 = 1.0
      prs0 = 1.0
      vln = 0.0
      bfn = 1.0

! Also set up the amplitude of the velocity perturbation in the
! transverse direction.

      vperp0 = 0.2

! Find the Alfven speed in the direction of propagation.

      vax = ABS ( bfn) / SQRT ( 4.0 * pi * rho0)

! Make wave numbers (kx,ky) and "phase".

      time = 0.0

      kx = nx * ( 2.0 * pi / ny)
      ky = ny * ( 2.0 * pi / ny)
      phase = - vax * time * ( 2.0 * pi / ny)

! ----------------------------------------------------------------------

! Make magnetic fields in active directions using a vector potential approach.
! That ensures divergence-free magnetic fields.

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1

            x = xb_save ( ixx)
            y = yb_save ( iyy)
            z = zc_save ( izz)
            tempphase = kx * x + ky * y + phase
            tempsin = SIN ( tempphase)
            tempcos = COS ( tempphase)
            tempaa = - bfn * ny * x + bfn * nx * y
     1             + vperp0 * ny * SQRT ( 4.0 * pi * rho0)
     1             * tempsin / ( 2.0 * pi)

            x = xb_save ( ixx)
            y = yb_save ( iyy - 1)
            z = zc_save ( izz)
            tempphase = kx * x + ky * y + phase
            tempsin = SIN ( tempphase)
            tempcos = COS ( tempphase)
            tempab = - bfn * ny * x + bfn * nx * y
     1             + vperp0 * ny * SQRT ( 4.0 * pi * rho0)
     1             * tempsin / ( 2.0 * pi)

            bfxbdy ( ixx, iyy, izz, 1) = ( tempaa - tempab)
     1                                / dy_save ( iyy)

          END DO
        END DO
      END DO

! ----------

      DO ixx = 1, ix1
        DO iyy = 0, iy1
          DO izz = 1, iz1

            x = xb_save ( ixx)
            y = yb_save ( iyy)
            z = zc_save ( izz)
            tempphase = kx * x + ky * y + phase
            tempsin = SIN ( tempphase)
            tempcos = COS ( tempphase)
            tempaa = - bfn * ny * x + bfn * nx * y
     1             + vperp0 * ny * SQRT ( 4.0 * pi * rho0)
     1             * tempsin / ( 2.0 * pi)

            x = xb_save ( ixx - 1)
            y = yb_save ( iyy)
            z = zc_save ( izz)
            tempphase = kx * x + ky * y + phase
            tempsin = SIN ( tempphase)
            tempcos = COS ( tempphase)
            tempab = - bfn * ny * x + bfn * nx * y
     1             + vperp0 * ny * SQRT ( 4.0 * pi * rho0)
     1             * tempsin / ( 2.0 * pi)

            bfybdy ( ixx, iyy, izz, 1) = - ( tempaa - tempab)
     1                                / dx_save ( ixx)

          END DO
        END DO
      END DO

! ----------

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 0, iz1

            x = xc_save ( ixx)
            y = yb_save ( iyy)
            z = zb_save ( izz)
            tempphase = kx * x + ky * y + phase
            tempsin = SIN ( tempphase)
            tempcos = COS ( tempphase)
            tempaa = - vperp0 * SQRT ( 4.0 * pi * rho0) * tempcos
     1             / ( 2.0 * pi)

            x = xc_save ( ixx)
            y = yb_save ( iyy - 1)
            z = zb_save ( izz)
            tempphase = kx * x + ky * y + phase
            tempsin = SIN ( tempphase)
            tempcos = COS ( tempphase)
            tempab = - vperp0 * SQRT ( 4.0 * pi * rho0) * tempcos
     1             / ( 2.0 * pi)

            bfzbdy ( ixx, iyy, izz, 1) = - ( tempaa - tempab)
     1                                / dy_save ( iyy)

          END DO
        END DO
      END DO

! ------------------------------

! Now set up interior zone-centered variables. The zone-centered velocities
! and (zone-averaged) magnetic fields are initialized by an exact, analytical
! value of "Sin[ kx x + ky y + phase]" and "Cos[ kx x + ky y + phase]"
! averaged over the domain [x1,x2]X[y1,y2]. As a result, this set-up
! routine is accurate for all orders!

      dtcur = ( vln + vperp0 + vax + SQRT (gamma * prs0 / rho0))
     1      / AMIN1( dx_save ( 1), dy_save ( 1), dz_save (1) )
      dtcur = cfl_coef_advct / dtcur

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1


          x = xc_save ( ixx)
          y = yc_save ( iyy)
          z = zc_save ( izz)

          x1 = x - 0.5 * dx_save ( ixx)
          x2 = x + 0.5 * dx_save ( ixx)

          y1 = y - 0.5 * dy_save ( iyy)
          y2 = y + 0.5 * dy_save ( iyy)

          tempac = kx * ky * ( x2 - x1) * ( y2 - y1)
          tempac = 1.0 / tempac

          tempsin = ( - SIN ( phase + kx * x1 + ky * y1)
     1                + SIN ( phase + kx * x2 + ky * y1)
     1                + SIN ( phase + kx * x1 + ky * y2)
     1                - SIN ( phase + kx * x2 + ky * y2) ) * tempac

          tempcos = ( - COS ( phase + kx * x1 + ky * y1)
     1                + COS ( phase + kx * x2 + ky * y1)
     1                + COS ( phase + kx * x1 + ky * y2)
     1                - COS ( phase + kx * x2 + ky * y2) ) * tempac

          vlx_tmp = vln * nx - vperp0 * ny * tempcos
          vly_tmp = vln * ny + vperp0 * nx * tempcos
          vlz_tmp = vperp0 * tempsin

          bfx_tmp = bfn * nx + vperp0 * ny * tempcos
     1                                  * SQRT ( 4.0 * pi * rho0)
          bfy_tmp = bfn * ny - vperp0 * nx * tempcos
     1                                  * SQRT ( 4.0 * pi * rho0)
          bfz_tmp = - vperp0 * tempsin
     1                         * SQRT ( 4.0 * pi * rho0)

          u ( ixx, iyy, izz, 1, 1) = rho0

          u ( ixx, iyy, izz, 2, 1) = rho0 * vlx_tmp
          u ( ixx, iyy, izz, 3, 1) = rho0 * vly_tmp
          u ( ixx, iyy, izz, 4, 1) = rho0 * vlz_tmp

          u ( ixx, iyy, izz, 5, 1) = prs0 / (gamma - 1.0)
     1  + 0.5 * rho0 * ( vlx_tmp**2 + vly_tmp**2 + vlz_tmp**2 )
     1  + ( bfx_tmp**2 + bfy_tmp**2 + bfz_tmp**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx_tmp
          u ( ixx, iyy, izz, 7, 1) = bfy_tmp
          u ( ixx, iyy, izz, 8, 1) = bfz_tmp

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ALFVEN_2D_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_MAGNETIZED_VORTEX_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a magnetized vortex in 2d. Default is xy-plane
! with iz1 = 1. The problem is catalogued in Balsara (2004) ApJS.
! Use a domain of [-5,5]x[-5,5]x[-5,5] on a 120x120x1 zone mesh for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 5.0/3.0" suggested. Suggested stopping time is
! "timestop = 10.0". The vortex propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Can analyze the
! accuracy of the problem by comparing the numerical solution with
! the analytical solution.

! Because this problem has several orders of magnitude variation in the
! magnetic field, the magnetic field divergence, if it is being diagnosed, 
! should be based on an unscaled, undivided divergence.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho0, prs0, vmag, bmag, vlx0, vly0,

     1     radius, tempaa, tempab, tempac,
     1     costheta, sintheta, velrms, bfldrms,
     1     VLXFN, VLYFN, BFXFN, BFYFN, VECPOTFN, PRSFN,
     1     x, y,
     1     xlo, xhi, ylo, yhi,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, bfxbdy_xp1, bfxbdy_xm1,
     1     bfybdy_yp1, bfybdy_ym1, bfzbdy_zp1, bfzbdy_zm1,
     1     prsmin, prsmax, maxdivb

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      BFXFN ( bmag, x, y) = - bmag * y * EXP ( 0.5*(1.0-x**2-y**2))
      BFYFN ( bmag, x, y) = bmag * x * EXP ( 0.5*(1.0-x**2-y**2))

      VECPOTFN ( bmag, x, y) = bmag * EXP ( 0.5*(1.0-x**2-y**2))

      PRSFN ( pi, rho0, prs0, vmag, bmag, x, y) = prs0
     1      - rho0 * 0.5 * vmag**2 * EXP ( 1.0-x**2-y**2)
     1      + bmag**2 * ( 1.0-x**2-y**2) * EXP ( 1.0-x**2-y**2)
     1      / (8.0*pi)

! ----------------------------------------------------------------------

! Default parameters for MHD magnetic vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0

      vlx0 = 1.0
      vly0 = 1.0

      vmag = 1.0 / ( 2.0 * pi)
      bmag = SQRT ( 4.0 * pi * rho0) / ( 2.0 * pi)

! ----------------------------------------------------------------------

! Make magnetic fields in active directions using a vector potential approach.
! That ensures divergence-free magnetic fields.

      dtcur = 1.0e30

      prsmin = 1.0e30
      prsmax = -1.0e30

      maxdivb = 0.0

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = xb_save ( ixx - 1)
          xhi = xb_save ( ixx)
          ylo = yb_save ( iyy - 1)
          yhi = yb_save ( iyy)

          rhouse = rho0

          vlxuse = vlx0 + 0.25 * VLXFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLXFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLXFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLXFN ( vmag, xhi, yhi)

          vlyuse = vly0 + 0.25 * VLYFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLYFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLYFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLYFN ( vmag, xhi, yhi)

          vlzuse = 0.0

          bfxuse = 0.25 * BFXFN ( bmag, xlo, ylo)
     1           + 0.25 * BFXFN ( bmag, xhi, ylo)
     1           + 0.25 * BFXFN ( bmag, xlo, yhi)
     1           + 0.25 * BFXFN ( bmag, xhi, yhi)

          bfyuse = 0.25 * BFYFN ( bmag, xlo, ylo)
     1           + 0.25 * BFYFN ( bmag, xhi, ylo)
     1           + 0.25 * BFYFN ( bmag, xlo, yhi)
     1           + 0.25 * BFYFN ( bmag, xhi, yhi)

          bfzuse = 0.0

          prsuse = 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xlo, ylo)
     1           + 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xhi, ylo)
     1           + 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xlo, yhi)
     1           + 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xhi, yhi)
          prsmin = MIN ( prsmin, prsuse)
          prsmax = MAX ( prsmax, prsuse)

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)

          u ( ixx, iyy, izz, 1, 1) = rhouse

          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

          u ( ixx, iyy, izz, 5, 1) = enguse

          u ( ixx, iyy, izz, 6, 1) = bfxuse
          u ( ixx, iyy, izz, 7, 1) = bfyuse
          u ( ixx, iyy, izz, 8, 1) = bfzuse

! ----------

          bfxbdy_xp1 = ( VECPOTFN ( bmag, xhi, yhi)
     1                 - VECPOTFN ( bmag, xhi, ylo) ) / (yhi - ylo)
          bfxbdy_xm1 = ( VECPOTFN ( bmag, xlo, yhi)
     1                 - VECPOTFN ( bmag, xlo, ylo) ) / (yhi - ylo)

          bfybdy_yp1 = - ( VECPOTFN ( bmag, xhi, yhi)
     1                   - VECPOTFN ( bmag, xlo, yhi) ) / (xhi - xlo)
          bfybdy_ym1 = - ( VECPOTFN ( bmag, xhi, ylo)
     1                   - VECPOTFN ( bmag, xlo, ylo) ) / (xhi - xlo)

          bfzbdy_zp1 = 0.0
          bfzbdy_zm1 = 0.0


          tempaa = ( bfxbdy_xp1 - bfxbdy_xm1) / (xhi - xlo)
     1           + ( bfybdy_yp1 - bfybdy_ym1) / (yhi - ylo)
          tempaa = ABS ( tempaa)
          maxdivb = MAX ( maxdivb, tempaa)


          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy_xp1
          bfxbdy ( ixx - 1, iyy, izz, 1) = bfxbdy_xm1

          bfybdy ( ixx, iyy, izz, 1) = bfybdy_yp1
          bfybdy ( ixx, iyy - 1, izz, 1) = bfybdy_ym1

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy_zp1
          bfzbdy ( ixx, iyy, izz - 1, 1) = bfzbdy_zm1

! ----------

          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2
          tempab = bfxuse**2 + bfyuse**2 + bfzuse**2

          tempab = tempab / ( 4.0 * pi * rhouse)
     1           + gamma * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_MAGNETIZED_VORTEX_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ERROR_MAGNETIZED_VORTEX_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Evaluates the error in the magnetized vortex in 2d. This
! subroutine should be placed in the timestep loop & it will
! keep writing out the erorr. Default is xy-plane
! with iz1 = 1. The problem is catalogued in Balsara (2004) ApJS.
! Use a domain of [-5,5]x[-5,5]x[-5,5] on a 120x120x1 zone mesh for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 5.0/3.0" suggested. Suggested stopping time is
! "timestop = 10.0". The vortex propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Can analyze the
! accuracy of the problem by comparing the numerical solution with
! the analytical solution.

! Because this problem has several orders of magnitude variation in the
! magnetic field, the magnetic field divergence, if it is being diagnosed, 
! should be based on an unscaled, undivided divergence.

! Some accuracy data is provided below. Notice that high quality reconstruction
! + good RS yields substantially better results on smaller meshes and
! a much faster approach to the theoretical accuracy.



! magvortex_ADER_WENO2_CFL0.65


! zones       bfxbdyerr_L1     accuracy_L1     bfxbdyerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64       3.1688E-3                        2.5430E-2
! 128x128     7.7893E-4        2.02            6.1421E-3        2.05
! 256x256     1.9341E-4        2.01            1.5527E-3        1.98
! ----------------------------------------------------------------------







! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz, shift_pattern

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     rho0, prs0, vmag, bmag, vlx0, vly0,

     1     radius, tempaa, tempab, tempac,
     1     costheta, sintheta, velrms, bfldrms,
     1     VLXFN, VLYFN, BFXFN, BFYFN, VECPOTFN, PRSFN,
     1     x, y,
     1     xlo, xhi, ylo, yhi,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, bfxbdy_xp1, bfxbdy_xm1,
     1     bfybdy_yp1, bfybdy_ym1, bfzbdy_zp1, bfzbdy_zm1,

     1     rhoerr_L1, rhoerr_Linf, engerr_L1, engerr_Linf,
     1     vlxerr_L1, vlxerr_Linf, vlyerr_L1, vlyerr_Linf,
     1     vlzerr_L1, vlzerr_Linf,
     1     bfxbdyerr_L1, bfxbdyerr_Linf,
     1     bfybdyerr_L1, bfybdyerr_Linf,
     1     bfzbdyerr_L1, bfzbdyerr_Linf

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      BFXFN ( bmag, x, y) = - bmag * y * EXP ( 0.5*(1.0-x**2-y**2))
      BFYFN ( bmag, x, y) = bmag * x * EXP ( 0.5*(1.0-x**2-y**2))

      VECPOTFN ( bmag, x, y) = bmag * EXP ( 0.5*(1.0-x**2-y**2))

      PRSFN ( pi, rho0, prs0, vmag, bmag, x, y) = prs0
     1      - rho0 * 0.5 * vmag**2 * EXP ( 1.0-x**2-y**2)
     1      + bmag**2 * ( 1.0-x**2-y**2) * EXP ( 1.0-x**2-y**2)
     1      / (8.0*pi)

! ----------------------------------------------------------------------

! Default parameters for MHD magnetic vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0

      vlx0 = 1.0
      vly0 = 1.0

      vmag = 1.0 / ( 2.0 * pi)
      bmag = SQRT ( 4.0 * pi * rho0) / ( 2.0 * pi)

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

      bfxbdyerr_L1 = 0.0
      bfxbdyerr_Linf = -1.0e32
      bfybdyerr_L1 = 0.0
      bfybdyerr_Linf = -1.0e32
      bfzbdyerr_L1 = 0.0
      bfzbdyerr_Linf = -1.0e32

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

          rhouse = rho0

          vlxuse = vlx0 + 0.25 * VLXFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLXFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLXFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLXFN ( vmag, xhi, yhi)

          vlyuse = vly0 + 0.25 * VLYFN ( vmag, xlo, ylo)
     1                  + 0.25 * VLYFN ( vmag, xhi, ylo)
     1                  + 0.25 * VLYFN ( vmag, xlo, yhi)
     1                  + 0.25 * VLYFN ( vmag, xhi, yhi)

          vlzuse = 0.0

          bfxuse = 0.25 * BFXFN ( bmag, xlo, ylo)
     1           + 0.25 * BFXFN ( bmag, xhi, ylo)
     1           + 0.25 * BFXFN ( bmag, xlo, yhi)
     1           + 0.25 * BFXFN ( bmag, xhi, yhi)

          bfyuse = 0.25 * BFYFN ( bmag, xlo, ylo)
     1           + 0.25 * BFYFN ( bmag, xhi, ylo)
     1           + 0.25 * BFYFN ( bmag, xlo, yhi)
     1           + 0.25 * BFYFN ( bmag, xhi, yhi)

          bfzuse = 0.0

          prsuse = 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xlo, ylo)
     1           + 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xhi, ylo)
     1           + 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xlo, yhi)
     1           + 0.25 * PRSFN ( pi, rho0, prs0, vmag, bmag, xhi, yhi)

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)


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

          bfxbdy_xp1 = ( VECPOTFN ( bmag, xhi, yhi)
     1                 - VECPOTFN ( bmag, xhi, ylo) ) / (yhi - ylo)
          bfxbdy_xm1 = ( VECPOTFN ( bmag, xlo, yhi)
     1                 - VECPOTFN ( bmag, xlo, ylo) ) / (yhi - ylo)

          bfybdy_yp1 = - ( VECPOTFN ( bmag, xhi, yhi)
     1                   - VECPOTFN ( bmag, xlo, yhi) ) / (xhi - xlo)
          bfybdy_ym1 = - ( VECPOTFN ( bmag, xhi, ylo)
     1                   - VECPOTFN ( bmag, xlo, ylo) ) / (xhi - xlo)

          bfzbdy_zp1 = 0.0
          bfzbdy_zm1 = 0.0


          tempaa = ABS ( bfxbdy ( ixx, iyy, izz, 1) - bfxbdy_xp1)
          bfxbdyerr_L1 = bfxbdyerr_L1 + tempaa
          bfxbdyerr_Linf = MAX ( bfxbdyerr_Linf, tempaa)

          tempaa = ABS ( bfybdy ( ixx, iyy, izz, 1) - bfybdy_yp1)
          bfybdyerr_L1 = bfybdyerr_L1 + tempaa
          bfybdyerr_Linf = MAX ( bfybdyerr_Linf, tempaa)

          tempaa = ABS ( bfzbdy ( ixx, iyy, izz, 1) - bfzbdy_zp1)
          bfzbdyerr_L1 = bfzbdyerr_L1 + tempaa
          bfzbdyerr_Linf = MAX ( bfzbdyerr_Linf, tempaa)

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

      bfxbdyerr_L1 = bfxbdyerr_L1 / FLOAT ( ix1 * iy1 * iz1)
      bfybdyerr_L1 = bfybdyerr_L1 / FLOAT ( ix1 * iy1 * iz1)
      bfzbdyerr_L1 = bfzbdyerr_L1 / FLOAT ( ix1 * iy1 * iz1)

! ----------------------------------------------------------------------

      WRITE ( 6, *)"ERROR_MAGNETIZED_VORTEX_XY gives :-"
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
      WRITE ( 6, *)"bfxbdyerr_L1, bfxbdyerr_Linf = ",
     1              bfxbdyerr_L1, bfxbdyerr_Linf
      WRITE ( 6, *)"bfybdyerr_L1, bfybdyerr_Linf = ",
     1              bfybdyerr_L1, bfybdyerr_Linf
      WRITE ( 6, *)"bfzbdyerr_L1, bfzbdyerr_Linf = ",
     1              bfzbdyerr_L1, bfzbdyerr_Linf
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ERROR_MAGNETIZED_VORTEX_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ISENTROPIC_VORTEX_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

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
! "gamma = 1.4" suggested. Suggested stopping time is
! "timestop = 10.0". The vortex propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Can analyze the
! accuracy of the problem by comparing the numerical solution with
! the analytical solution.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

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
     1     bfxuse, bfyuse, bfzuse, bfxbdy_xp1, bfxbdy_xm1,
     1     bfybdy_yp1, bfybdy_ym1, bfzbdy_zp1, bfzbdy_zm1,
     1     prsmin, prsmax

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      TEMPFN ( pi, vmag, gamma, x, y) = - ( gamma - 1.0) * vmag**2
     1                          * EXP (1.0-x**2-y**2) / ( 2.0 * gamma)

! ----------------------------------------------------------------------

! Default parameters for isentropic hydrodynamical vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0
      entropy0 = prs0 / rho0**gamma

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

          bfxuse = 0.0
          bfyuse = 0.0
          bfzuse = 0.0

          delta_t = 0.25 * TEMPFN ( pi, vmag, gamma, xlo, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma, xhi, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma, xlo, yhi)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma, xhi, yhi)

          rhouse = ( ( ( prs0 / rho0) + delta_t)
     1               / entropy0)**(1.0/(gamma-1.0) )

          prsuse = rhouse * ( ( prs0 / rho0) + delta_t)

          prsmin = MIN ( prsmin, prsuse)
          prsmax = MAX ( prsmax, prsuse)

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)

          u ( ixx, iyy, izz, 1, 1) = rhouse

          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

          u ( ixx, iyy, izz, 5, 1) = enguse

          u ( ixx, iyy, izz, 6, 1) = 0.0
          u ( ixx, iyy, izz, 7, 1) = 0.0
          u ( ixx, iyy, izz, 8, 1) = 0.0

! ----------

          bfxbdy ( ixx, iyy, izz, 1) = 0.0
          bfxbdy ( ixx - 1, iyy, izz, 1) = 0.0

          bfybdy ( ixx, iyy, izz, 1) = 0.0
          bfybdy ( ixx, iyy - 1, izz, 1) = 0.0

          bfzbdy ( ixx, iyy, izz, 1) = 0.0
          bfzbdy ( ixx, iyy, izz - 1, 1) = 0.0

! ----------

          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2
          tempab = bfxuse**2 + bfyuse**2 + bfzuse**2

          tempab = tempab / ( 4.0 * pi * rhouse)
     1           + gamma * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ISENTROPIC_VORTEX_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ERROR_ISENTROPIC_VORTEX_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,

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
! "gamma = 1.4" suggested. Suggested stopping time is
! "timestop = 10.0". The vortex propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Can analyze the
! accuracy of the problem by comparing the numerical solution with
! the analytical solution.

! Some accuracy data is provided below. Notice that high quality reconstruction
! + good RS yields substantially better results on smaller meshes and
! a much faster approach to the theoretical accuracy.


! hydrovortex_ADER_WENO2_CFL0.8

! zones       rhoerr_L1     accuracy_L1     rhoerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64       9.5949E-4                     1.5694E-2
! 128x128     1.9390E-4     2.31            3.4514E-3     2.18
! 256x256     4.1801E-5     2.21            7.3805E-4     2.23
! ----------------------------------------------------------------------






! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

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
     1     bfxuse, bfyuse, bfzuse, bfxbdy_xp1, bfxbdy_xm1,
     1     bfybdy_yp1, bfybdy_ym1, bfzbdy_zp1, bfzbdy_zm1,

     1     rhoerr_L1, rhoerr_Linf, engerr_L1, engerr_Linf,
     1     vlxerr_L1, vlxerr_Linf, vlyerr_L1, vlyerr_Linf,
     1     vlzerr_L1, vlzerr_Linf

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      TEMPFN ( pi, vmag, gamma, x, y) = - ( gamma - 1.0) * vmag**2 
     1                          * EXP (1.0-x**2-y**2) / ( 2.0 * gamma)

! ----------------------------------------------------------------------

! Default parameters for isentropic vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0
      entropy0 = prs0 / rho0**gamma

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

          bfxuse = 0.0
          bfyuse = 0.0
          bfzuse = 0.0

          delta_t = 0.25 * TEMPFN ( pi, vmag, gamma, xlo, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma, xhi, ylo)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma, xlo, yhi)
     1            + 0.25 * TEMPFN ( pi, vmag, gamma, xhi, yhi)

          rhouse = ( ( ( prs0 / rho0) + delta_t)
     1               / entropy0)**(1.0/(gamma-1.0) )

          prsuse = rhouse * ( ( prs0 / rho0) + delta_t)

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)

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

      END SUBROUTINE ERROR_ISENTROPIC_VORTEX_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_FIELD_LOOP_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a magnetized field loop in 2d. Default is xy-plane
! with iz1 = 1. The problem is catalogued in Gardiner & Stone (2005) JCP.
! Use a domain of [-1,1]x[-0.5,0.5]x[-1,1] on a 128x64x1 zone mesh for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 1.4" suggested. Suggested stopping time is
! "timestop = 1.0". The field loop propagates along the diagonal of
! the computing mesh in a form-preserving fashion. Because of the
! discontinuity in the field at the loop's boundary and the loop's center,
! this problem is not suitable for accuracy analysis.

! Because this problem has several orders of magnitude variation in the
! magnetic field, the magnetic field divergence, if it is being diagnosed, 
! should be based on an unscaled, undivided divergence.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho0, prs0, amag, vlx0, vly0,

     1     radius, tempaa, tempab, tempac,
     1     costheta, sintheta, velrms, bfldrms,
     1     VECPOTFN,
     1     x, y,
     1     xlo, xhi, ylo, yhi,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, bfxbdy_xp1, bfxbdy_xm1,
     1     bfybdy_yp1, bfybdy_ym1, bfzbdy_zp1, bfzbdy_zm1,
     1     prsmin, prsmax, maxdivb

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      VECPOTFN ( amag, x, y) = amag * ( 0.3 - SQRT ( x**2 + y**2))
     1 * ( 0.5 + SIGN ( 0.5, 0.3 - SQRT ( x**2 + y**2) ) )

! ----------------------------------------------------------------------

! Default parameters for MHD magnetic vortex problem.
! Stop this one at a time of 10.0.

      rho0 = 1.0
      prs0 = 1.0

      vlx0 = SQRT ( 5.0) * ( 2.0 / SQRT ( 5.0))
      vly0 = SQRT ( 5.0) * ( 1.0 / SQRT ( 5.0))

      amag = SQRT ( 4.0 * pi * rho0) * 1.0e-3

! ----------------------------------------------------------------------

! Make magnetic fields in active directions using a vector potential approach.
! That ensures divergence-free magnetic fields.

      dtcur = 1.0e30

      prsmin = 1.0e30
      prsmax = -1.0e30

      maxdivb = 0.0

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = xb_save ( ixx - 1)
          xhi = xb_save ( ixx)
          ylo = yb_save ( iyy - 1)
          yhi = yb_save ( iyy)

! ----------

          bfxbdy_xp1 = ( VECPOTFN ( amag, xhi, yhi)
     1                 - VECPOTFN ( amag, xhi, ylo) ) / (yhi - ylo)
          bfxbdy_xm1 = ( VECPOTFN ( amag, xlo, yhi)
     1                 - VECPOTFN ( amag, xlo, ylo) ) / (yhi - ylo)

          bfybdy_yp1 = - ( VECPOTFN ( amag, xhi, yhi)
     1                   - VECPOTFN ( amag, xlo, yhi) ) / (xhi - xlo)
          bfybdy_ym1 = - ( VECPOTFN ( amag, xhi, ylo)
     1                   - VECPOTFN ( amag, xlo, ylo) ) / (xhi - xlo)

          bfzbdy_zp1 = 0.0
          bfzbdy_zm1 = 0.0


          tempaa = ( bfxbdy_xp1 - bfxbdy_xm1) / (xhi - xlo)
     1           + ( bfybdy_yp1 - bfybdy_ym1) / (yhi - ylo)
          tempaa = ABS ( tempaa)
          maxdivb = MAX ( maxdivb, tempaa)


          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy_xp1
          bfxbdy ( ixx - 1, iyy, izz, 1) = bfxbdy_xm1

          bfybdy ( ixx, iyy, izz, 1) = bfybdy_yp1
          bfybdy ( ixx, iyy - 1, izz, 1) = bfybdy_ym1

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy_zp1
          bfzbdy ( ixx, iyy, izz - 1, 1) = bfzbdy_zm1

          bfxuse = 0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                   + bfxbdy ( ixx - 1, iyy, izz, 1) )

          bfyuse = 0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                   + bfybdy ( ixx, iyy - 1, izz, 1) )

          bfzuse = 0.0

! ----------

          rhouse = rho0

          vlxuse = vlx0

          vlyuse = vly0

          vlzuse = 0.0

          prsuse = prs0
          prsmin = MIN ( prsmin, prsuse)
          prsmax = MAX ( prsmax, prsuse)

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)

          u ( ixx, iyy, izz, 1, 1) = rhouse

          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

          u ( ixx, iyy, izz, 5, 1) = enguse

          u ( ixx, iyy, izz, 6, 1) = bfxuse
          u ( ixx, iyy, izz, 7, 1) = bfyuse
          u ( ixx, iyy, izz, 8, 1) = bfzuse

! ----------

          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2
          tempab = bfxuse**2 + bfyuse**2 + bfzuse**2

          tempab = tempab / ( 4.0 * pi * rhouse)
     1           + gamma * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_FIELD_LOOP_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_FIELD_LOOP_3D

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a magnetized field loop in 3d. The problem is catalogued in
! Gardiner & Stone (2008) JCP. A domain of [-0.5,0.5]x[-0.5,0.5]x[-1,1]
! is used with meshes ranging from 32x32x64 to 128x128x256 having periodic
! boundaries (boundary condition # 4). Use "gamma = 1.4" and "timestop = 1.0".
! Because of the discontinuity in the field at the loop's boundary and
! the loop's center, this problem is not suitable for accuracy analysis.

! Because this problem has several orders of magnitude variation in the
! magnetic field, the magnetic field divergence, if it is being diagnosed,
! should be based on an unscaled, undivided divergence.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0, div_max,

     1     x, y, z, Ax, Ay, Az,

     1     A_x_pl, A_x_mi, A_y_pl, A_y_mi, A_z_pl, A_z_mi,

     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      rho0 = 1.0
      prs0 = 1.0

      vlx0 = 1.0
      vly0 = 1.0
      vlz0 = 2.0

! ----------------------------------------------------------------------

! Step a: We make a numerical integraton of the vector potential's
! components along the zone edges. We then differentiate this vector
! potential to get the magnetic field.

! ------------------------------

! Make "bfxbdy" using a magnetic vector potential approach.

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1

! -----

          x = xb_save ( ixx)

          y = yb_save ( iyy)

          z = zc_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_z_pl = Az

! -----

          x = xb_save ( ixx)

          y = yb_save ( iyy - 1)

          z = zc_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_z_mi = Az

! -----

          x = xb_save ( ixx)

          y = yc_save ( iyy)

          z = zb_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_y_pl = Ay

! -----

          x = xb_save ( ixx)

          y = yc_save ( iyy)

          z = zb_save ( izz - 1)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_y_mi = Ay

! -----

          bfxbdy ( ixx, iyy, izz, 1) =
     1    ( A_z_pl - A_z_mi ) / dy_save ( iyy)
     1  - ( A_y_pl - A_y_mi ) / dz_save ( izz)

! -----

          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

! Make "bfybdy" using a magnetic vector potential approach.

      DO ixx = 1, ix1
        DO iyy = 0, iy1
          DO izz = 1, iz1

! -----

          x = xb_save ( ixx)

          y = yb_save ( iyy)

          z = zc_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_z_pl = Az

! -----

          x = xb_save ( ixx - 1)

          y = yb_save ( iyy)

          z = zc_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_z_mi = Az

! -----

          x = xc_save ( ixx)

          y = yb_save ( iyy)

          z = zb_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_x_pl = Ax

! -----

          x = xc_save ( ixx)

          y = yb_save ( iyy)

          z = zb_save ( izz - 1)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_x_mi = Ax

! -----

          bfybdy ( ixx, iyy, izz, 1) =
     1    ( A_x_pl - A_x_mi ) / dz_save ( izz)
     1  - ( A_z_pl - A_z_mi ) / dx_save ( ixx)

! -----

          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

! Make "bfzbdy" using a magnetic vector potential approach.

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 0, iz1

! -----

          x = xc_save ( ixx)

          y = yb_save ( iyy)

          z = zb_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_x_pl = Ax

! -----

          x = xc_save ( ixx)

          y = yb_save ( iyy - 1)

          z = zb_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_x_mi = Ax

! -----

          x = xb_save ( ixx)

          y = yc_save ( iyy)

          z = zb_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_y_pl = Ay

! -----

          x = xb_save ( ixx - 1)

          y = yc_save ( iyy)

          z = zb_save ( izz)

          CALL FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

          A_y_mi = Ay

! -----

          bfzbdy ( ixx, iyy, izz, 1) =
     1    ( A_y_pl - A_y_mi ) / dx_save ( ixx)
     1  - ( A_x_pl - A_x_mi ) / dy_save ( iyy)

! -----

          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

! Step b: Test for div-free aspect in the magnetic field.
! Since this test is based on the integral form of Stokes law it
! is accurate at all orders.

! ------------------------------

      div_max = 0.0

      DO ixx = 1, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1

! ----------

          tempaa = ( bfxbdy ( ixx, iyy, izz, 1)
     1             - bfxbdy ( ixx-1, iyy, izz, 1))
     1           * dy_save ( iyy) * dz_save ( izz)
     1           + ( bfybdy ( ixx, iyy, izz, 1)
     1             - bfybdy ( ixx, iyy-1, izz, 1))
     1           * dx_save ( ixx) * dz_save ( izz)
     1           + ( bfzbdy ( ixx, iyy, izz, 1)
     1             - bfzbdy ( ixx, iyy, izz-1, 1))
     1           * dx_save ( ixx) * dy_save ( iyy)

          div_max = AMAX1 ( div_max, tempaa)

! ----------

          END DO
        END DO
      END DO

      IF ( div_max .GT. smallnum) THEN
      WRITE ( 6, *)"WARNING : div_max .GT. smallnum, div_max = ",
     1              div_max
      END IF

      WRITE ( 6, *)"From b-field setup step: div_max = ",
     1                             div_max

! ----------------------------------------------------------------------

! Make all the primitive (and therefrom conserved) variables at the zone center.

      dtcur = 1.0e30

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          bfx0 = 0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                 + bfxbdy ( ixx-1, iyy, izz, 1))

          bfy0 = 0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                 + bfybdy ( ixx, iyy-1, izz, 1))

          bfz0 = 0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                 + bfzbdy ( ixx, iyy, izz-1, 1))

          u ( ixx, iyy, izz, 1, 1) = rho0

          u ( ixx, iyy, izz, 2, 1) = rho0 * vlx0
          u ( ixx, iyy, izz, 3, 1) = rho0 * vly0
          u ( ixx, iyy, izz, 4, 1) = rho0 * vlz0

          u ( ixx, iyy, izz, 5, 1) = prs0 / ( gamma - 1.0)
     1    + 0.5 * rho0 * ( vlx0**2 + vly0**2 + vlz0**2)
     1    + ( bfx0**2 + bfy0**2 + bfz0**2) / ( 8.0 * pi)

          u ( ixx, iyy, izz, 6, 1) = bfx0
          u ( ixx, iyy, izz, 7, 1) = bfy0
          u ( ixx, iyy, izz, 8, 1) = bfz0

          tempaa = SQRT ( vlx0**2 + vly0**2 + vlz0**2)
          tempab = gamma * prs0
     1           + ( bfx0**2 + bfy0**2 + bfz0**2) / ( 4.0 * pi)
          tempac = SQRT ( tempab / rho0)

          tempad = cfl_coef_advct
     1      * AMIN1 ( dx_save ( ixx), dy_save ( iyy), dz_save ( izz))
     1           / ( tempaa + tempac)

          dtcur = AMIN1 ( dtcur, tempad)

          END DO

        END DO

      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE INIT_FIELD_LOOP_3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ERROR_FIELD_LOOP_3D

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Measures the total magnetic energy as a function of time. Also measures
! the absolute value of the magnetic field along the axis of the vortex.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     bfx0, bfy0, bfz0, ob5, tb5, b2, b3_abs

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      ob5 = 1.0 / SQRT ( 5.0)
      tb5 = 2.0 / SQRT ( 5.0)

! ----------------------------------------------------------------------

      b2 = 0.0
      b3_abs = 0.0

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          bfx0 = 0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                 + bfxbdy ( ixx-1, iyy, izz, 1))

          bfy0 = 0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                 + bfybdy ( ixx, iyy-1, izz, 1))

          bfz0 = 0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                 + bfzbdy ( ixx, iyy, izz-1, 1))

          b2 = b2 + ( bfx0**2 + bfy0**2 + bfz0**2 )
     1            * dx_save ( ixx) * dy_save ( iyy) * dz_save ( izz)

!          b3_abs = b3_abs + ABS ( - ob5 * bfx0 + tb5 * bfz0)
          b3_abs = b3_abs + ( ( - ob5 * bfx0 + tb5 * bfz0)**2 )
     1            * dx_save ( ixx) * dy_save ( iyy) * dz_save ( izz)

          END DO

        END DO

      END DO

      WRITE ( 7, *)time, b2, b3_abs

! ----------------------------------------------------------------------

      END SUBROUTINE ERROR_FIELD_LOOP_3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#if ( 0 == 1)
! This little code takes input from previous subroutine and turns it into
! properly normalized output.

      INTEGER i, npoints
      PARAMETER ( npoints = 203)
      REAL time ( npoints), b2 ( npoints), b3_abs ( npoints), b0

      OPEN ( unit = 7, file = "fort.7")
      OPEN ( unit = 8, file = "bmag0001")
      OPEN ( unit = 9, file = "bx3c0001")

      DO i = 1, npoints
      READ ( 7, *) time ( i), b2 ( i), b3_abs ( i)
      IF ( i == 1) b0 = b2 ( 1)
      b2 ( i) = b2 ( i) / b0
      b3_abs ( i) = b3_abs ( i) / b0
      END DO

      DO i = 1, npoints
      WRITE ( 8, *) b2 ( i)
      END DO

      DO i = 1, npoints
!      WRITE ( 9, *) b3_abs ( i)
      WRITE ( 9, *) ALOG10 ( b3_abs ( i) )
      END DO

      DO i = 1, npoints
      WRITE ( 8, *) time ( i)
      WRITE ( 9, *) time ( i)
      END DO

      STOP
      END

#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FIELD_LOOP_VECPOT ( x, y, z, Ax, Ay, Az)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER i

      REAL x, y, z, Ax, Ay, Az,
     1     x1, x2, x3, Ax1, Ax2, Ax3, ob5, tb5, r

! ----------------------------------------------------------------------

      ob5 = 1.0 / SQRT ( 5.0)
      tb5 = 2.0 / SQRT ( 5.0)

      x1 = tb5 * x + ob5 * z
      x2 = y
      x3 = - ob5 * x + tb5 * z

      IF ( x1 < -0.5 * tb5) THEN

        DO i = 1, 10
        x1 = x1 + tb5
        IF ( x1 >= -0.5 * tb5) GO TO 1010
        END DO
 
      END IF

      IF ( x1 > 0.5 * tb5) THEN

        DO i = 1, 10
        x1 = x1 - tb5
        IF ( x1 <= 0.5 * tb5) GO TO 1010
        END DO

      END IF

1010  CONTINUE

      Ax1 = 0.0
      Ax2 = 0.0

      r = SQRT ( x1**2 + x2**2)

      IF ( r <= 0.3) THEN
        Ax3 = 1.0e-3 * ( 0.3 - r)
      ELSE
        Ax3 = 0.0
      END IF

      Ax = - ob5 * Ax3
      Ay = 0.0
      Az = tb5 * Ax3

! ----------------------------------------------------------------------

      END SUBROUTINE FIELD_LOOP_VECPOT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ORZAG_TANG_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up the Orzag-Tang problem in 2d. Default is xy-plane with iz1 = 1.
! Use a domain of [0,2]x[0,2]x[-1,1] on a 200x200x1 zone mesh for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 5.0/3.0" suggested. Suggested stopping time is
! "timestop = 1.0". The problem is thought to mimic the transition to
! MHD turbulence, see Dahlburg & Picone (1989) Phys. of Fluids A.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho0, prs0,

     1     radius, tempaa, tempab, tempac,
     1     costheta, sintheta, velrms, bfldrms,
     1     VECPOTFN,
     1     x, y,
     1     xlo, xhi, ylo, yhi,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, bfxbdy_xp1, bfxbdy_xm1,
     1     bfybdy_yp1, bfybdy_ym1, bfzbdy_zp1, bfzbdy_zm1,
     1     prsmin, prsmax, maxdivb

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      VECPOTFN ( pi, x, y) = - SQRT ( 4.0 * pi)
     1    * ( COS ( 2.0 * pi * x) + 2.0 * COS ( pi * y) ) / ( 2.0 * pi)

! ----------------------------------------------------------------------

! Default parameters for MHD magnetic vortex problem.
! Stop this one at a time of 10.0.

      rho0 = gamma**2
      prs0 = gamma

! ----------------------------------------------------------------------

! Make magnetic fields in active directions using a vector potential approach.
! That ensures divergence-free magnetic fields.

      dtcur = 1.0e30

      prsmin = 1.0e30
      prsmax = -1.0e30

      maxdivb = 0.0

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = xb_save ( ixx - 1)
          xhi = xb_save ( ixx)
          ylo = yb_save ( iyy - 1)
          yhi = yb_save ( iyy)

! ----------

          bfxbdy_xp1 = ( VECPOTFN ( pi, xhi, yhi)
     1                 - VECPOTFN ( pi, xhi, ylo) ) / (yhi - ylo)
          bfxbdy_xm1 = ( VECPOTFN ( pi, xlo, yhi)
     1                 - VECPOTFN ( pi, xlo, ylo) ) / (yhi - ylo)

          bfybdy_yp1 = - ( VECPOTFN ( pi, xhi, yhi)
     1                   - VECPOTFN ( pi, xlo, yhi) ) / (xhi - xlo)
          bfybdy_ym1 = - ( VECPOTFN ( pi, xhi, ylo)
     1                   - VECPOTFN ( pi, xlo, ylo) ) / (xhi - xlo)

          bfzbdy_zp1 = 0.0
          bfzbdy_zm1 = 0.0


          tempaa = ( bfxbdy_xp1 - bfxbdy_xm1) / (xhi - xlo)
     1           + ( bfybdy_yp1 - bfybdy_ym1) / (yhi - ylo)
          tempaa = ABS ( tempaa)
          maxdivb = MAX ( maxdivb, tempaa)


          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy_xp1
          bfxbdy ( ixx - 1, iyy, izz, 1) = bfxbdy_xm1

          bfybdy ( ixx, iyy, izz, 1) = bfybdy_yp1
          bfybdy ( ixx, iyy - 1, izz, 1) = bfybdy_ym1

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy_zp1
          bfzbdy ( ixx, iyy, izz - 1, 1) = bfzbdy_zm1

          bfxuse = 0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                   + bfxbdy ( ixx - 1, iyy, izz, 1) )

          bfyuse = 0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                   + bfybdy ( ixx, iyy - 1, izz, 1) )

          bfzuse = 0.0

! ----------

          rhouse = rho0

          vlxuse = - SIN ( pi * yc_save ( iyy) )

          vlyuse = SIN ( pi * xc_save ( ixx) )

          vlzuse = 0.0

          prsuse = prs0
          prsmin = MIN ( prsmin, prsuse)
          prsmax = MAX ( prsmax, prsuse)

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)

          u ( ixx, iyy, izz, 1, 1) = rhouse

          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

          u ( ixx, iyy, izz, 5, 1) = enguse

          u ( ixx, iyy, izz, 6, 1) = bfxuse
          u ( ixx, iyy, izz, 7, 1) = bfyuse
          u ( ixx, iyy, izz, 8, 1) = bfzuse

! ----------

          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2
          tempab = bfxuse**2 + bfyuse**2 + bfzuse**2

          tempab = tempab / ( 4.0 * pi * rhouse)
     1           + gamma * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ORZAG_TANG_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_DOUBLE_MACH_2D_XZ

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up the 2d hydrodynamical double Mach reflection problem in the xz-plane
! with "iy1 = 1, ioffy = 0". 480x1x120 zones is the traditional resolution.
! The problem is catalogued in Woodward & Colella (1984) JCP.
! Set up on a [0,4]x[-1,1]x[0,1] rectangular domain. Set "gamma=1.4"
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
!      eng_tmp = 116.5 / ( gamma - 1.0) + 0.5 * 8.0 * 8.25 ** 2
!      bfx_tmp = 0.0
!      bfy_tmp = 0.0
!      bfz_tmp = 0.0

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

! At the lower z-boundary of PAD_BOUNDARY_Z, insert the same boundary
! conditions as in Step I.

! Step IV) Upper z-boundary has to be set via call to ADD_DOUBLE_MACH_2D_XZ.
! This call should be made once in each timestep loop and the specific
! location is specified in that subroutine's notes.
! Set upper z-boundary condition to something non-interfering (like 17)
! because it is set exclusively by the ADD_DOUBLE_MACH_2D_XZ subroutine.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

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
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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
      tempaa = mach21 * SQRT ( gamma * prs2 / rho2)

      velshk = tempaa

! And (40.08) of SFSW gives:
      tempab = ( gamma - 1.0) / ( gamma + 1.0)

! (67.07) gives:
      prs1 = prs2 * ( ( 1.0 + tempab) * mach21**2 - tempab )

! (67.02) gives:
      rho1 = rho2 * ( prs1 + tempab * prs2) / ( prs2 + tempab * prs1)

! (67.10) gives ( we expressly use the fact that vlx2, vlz2 = 0.0) :
      vlx1 = ( 1.0 - tempab) * ( mach21 - 1.0 / mach21)
      vlx1 = vlx1 * SQRT ( gamma * prs2 / rho2)

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

! Set magnetic fields to zero.

          bfxbdy ( :, :, :, 1) = 0.0
          bfybdy ( :, :, :, 1) = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

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

          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

          u ( ixx, iyy, izz, 6, 1) = 0.0
          u ( ixx, iyy, izz, 7, 1) = 0.0
          u ( ixx, iyy, izz, 8, 1) = 0.0

! ----------

          tempaa = SQRT ( vlxuse**2 + vlyuse**2 + vlzuse**2 )

          tempab = SQRT ( gamma * prsuse / rhouse )

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_DOUBLE_MACH_2D_XZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_DOUBLE_MACH_2D_XZ

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Upper z-boundary of the double Mach reflection problem
! has to be set via call to this subroutine -- ADD_DOUBLE_MACH_2D_XZ.
! This call should be made once in each timestep loop. Best place
! to put it is after the padding of the ghost zones but before the
! call to the limiter.

! Note that the interface of this subroutine is different from that of
! INIT_DOUBLE_MACH_2D_XZ. The "dtcur" has been replaced with "time".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

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
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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
      tempaa = mach21 * SQRT ( gamma * prs2 / rho2)

      velshk = tempaa

! And (40.08) of SFSW gives:
      tempab = ( gamma - 1.0) / ( gamma + 1.0)

! (67.07) gives:
      prs1 = prs2 * ( ( 1.0 + tempab) * mach21**2 - tempab )

! (67.02) gives:
      rho1 = rho2 * ( prs1 + tempab * prs2) / ( prs2 + tempab * prs1)

! (67.10) gives ( we expressly use the fact that vlx2, vlz2 = 0.0) :
      vlx1 = ( 1.0 - tempab) * ( mach21 - 1.0 / mach21)
      vlx1 = vlx1 * SQRT ( gamma * prs2 / rho2)

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

          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADD_DOUBLE_MACH_2D_XZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_FF_STEP_2D_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

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
! Set up on a [0,3]x[0,1]x[-1,1] rectangular domain. Set "gamma=1.4"
! and "timestop = 4.0".

! Declare and set "nxstep = (48 * ix1) / 240; nystep = (16 * iy1) / 80". 
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
!      eng_tmp = 1.0/(gamma-1.0) + 0.5 * 1.4 * 3.0**2
!      bfx_tmp = 0.0
!      bfy_tmp = 0.0
!      bfz_tmp = 0.0

! Step II) Because the reflection at the top x-face of the step is most
! important, call PAD_BOUNDARY_Y before PAD_BOUNDARY_X.

! Step III) After padding boundary and before limiter, call ADD_FF_STEP_2D_XY.

! STEP IV) We can cover up the flow variables in the step with:
!      r4gr ( nxstep + 1: ix1, 1: nystep, :) = r4gr ( 1, 1, 1)

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rhouse, prsuse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      bfxuse = 0.0
      bfyuse = 0.0
      bfzuse = 0.0

! ----------------------------------------------------------------------

! Set magnetic fields to zero.

          bfxbdy ( :, :, :, 1) = 0.0
          bfybdy ( :, :, :, 1) = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

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

          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

          u ( ixx, iyy, izz, 6, 1) = 0.0
          u ( ixx, iyy, izz, 7, 1) = 0.0
          u ( ixx, iyy, izz, 8, 1) = 0.0

! ----------

          tempaa = SQRT ( vlxuse**2 + vlyuse**2 + vlzuse**2 )

          tempab = SQRT ( gamma * prsuse / rhouse )

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_FF_STEP_2D_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_FF_STEP_2D_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,

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
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

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
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      prs_step = (gamma - 1.0) * ( u ( nxstep, nystep, izz, 5, 1)
     1 - 0.5 * rho_step * ( vlx_step**2 + vly_step**2 + vlz_step**2) )
      prs_step = AMAX1 ( prs_step, prsfloor)

      entropy_step = prs_step / rho_step**gamma
      enthalpy_step = ( gamma * prs_step) / ( rho_step * ( gamma - 1.0))
     1 + 0.5 * ( vlx_step**2 + vly_step**2 + vlz_step**2)

! ----------------------------------------------------------------------

! Zero out any stray magnetic field buildup.

      bfxbdy ( :, :, :, 1) = 0.0

      bfybdy ( :, :, :, 1) = 0.0

      bfzbdy ( :, :, :, 1) = 0.0

! ----------------------------------------------------------------------

! Entropy fix, zone 1.

      ixx = nxstep + 1
      iyy = nystep + 1
      izz = 1

      rhouse = AMAX1 ( u ( ixx, iyy, izz, 1, 1), rhofloor)

      vlxuse = u ( ixx, iyy, izz, 2, 1) / rhouse
      vlyuse = u ( ixx, iyy, izz, 3, 1) / rhouse
      vlzuse = u ( ixx, iyy, izz, 4, 1) / rhouse

      prsuse = ( gamma - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma * prsuse) / ( rhouse * ( gamma - 1.0))

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

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma - 1.0)
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

      prsuse = ( gamma - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma * prsuse) / ( rhouse * ( gamma - 1.0))

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

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma - 1.0)
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

      prsuse = ( gamma - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma * prsuse) / ( rhouse * ( gamma - 1.0))

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

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma - 1.0)
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

      prsuse = ( gamma - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma * prsuse) / ( rhouse * ( gamma - 1.0))

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

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma - 1.0)
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

      prsuse = ( gamma - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma * prsuse) / ( rhouse * ( gamma - 1.0))

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

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma - 1.0)
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

      prsuse = ( gamma - 1.0) * ( u ( ixx, iyy, izz, 5, 1)
     1 - 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2) )
      prsuse = AMAX1 ( prsuse, prsfloor)

! Reset density based on matching entropy in this zone to that in front of 
! the step.

      rhouse = ( prsuse / entropy_step)**( 1.0/gamma)
      u ( ixx, iyy, izz, 1, 1) = rhouse

! Reset velocities in this zone based on matching the kinetic part of
! the free stream specific enthalpy in this zone to that in front of 
! the step.

      td_enthalpyuse = ( gamma * prsuse) / ( rhouse * ( gamma - 1.0))

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

      u ( ixx, iyy, izz, 5, 1) = prsuse / ( gamma - 1.0)
     1 + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

! ----------------------------------------------------------------------

! At a distance of "ioffx, ioffy" zones inside the step, the values will
! certainly not be reset by the padbdy routines. Therefore, give all those
! zones completely inoffensive values. Keep resetting those values at
! every timestep. The values from zone ( 1, iy1, 1) are used.

      izz = 1

      DO iyy = 1 - ioffy, nystep - ioffy
      DO ixx = nxstep + ioffx, ix1 + ioffx

      u ( ixx, iyy, izz, :, 1) = u ( 1, iy1, 1, :, 1)
      u ( ixx, iyy, izz, :, 2:) = 0.0

      END DO
      END DO

! ----------------------------------------------------------------------

! Do the same initialization for the two square ghost regions that lie
! immediately outside the inflow boundary.

      izz = 1

      DO iyy = 1 - ioffy, 0
      DO ixx = 1 - ioffx, 0

      u ( ixx, iyy, izz, :, 1) = u ( 1, iy1, 1, :, 1)
      u ( ixx, iyy, izz, :, 2:) = 0.0

      END DO
      END DO

! ------------------------------

      izz = 1

      DO iyy = iy1 + 1, iy1 + ioffy
      DO ixx = 1 - ioffx, 0

      u ( ixx, iyy, izz, :, 1) = u ( 1, iy1, 1, :, 1)
      u ( ixx, iyy, izz, :, 2:) = 0.0

      END DO
      END DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADD_FF_STEP_2D_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_TWO_D_RIEM_PROBLEM_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds three variants of two-dimensional Riemann problems. The first two
! are drawn from Brio, Zakharian & Webb (2001) JCP, the third one from
! Schulz-Rinne, Collins & Galz (1993) SIAM J. Sci. Comput.
! The problems are set up in the XY
! plane with 400x400x1 zones spanning [-1,1]x[-1,1]x[-1,1].
! Set "gamma = 1.4" and keep the boundaries at continuitive (no. 2).
! Such 2d Riemann problems consist of four constant states in the 4 quadrants.
! The first problem has two initial weak shocks and two slip lines and
! is stopped at "timestop = 0.62". The second problem produces double
! Mach reflection and a shock propagating at an angle to the grid lines
! and is stopped at "timestop = 1.1". The third problem creates interacting
! shock structures which result in three-shock configurations. The subsonic
! area is bounded by the joining Mach shocks and the reflected shocks so that
! it has an oval shape. It should be stopped at "timestop = 0.55".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz, iproblem

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, prs1, vlx1, vly1, rho2, prs2, vlx2, vly2,
     1     rho3, prs3, vlx3, vly3, rho4, prs4, vlx4, vly4,
     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse,
     1     tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      iproblem = 1

      IF ( iproblem == 1) THEN
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

      ELSE IF ( iproblem == 2) THEN
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

      ELSE

      prs1 = 1.1
      rho1 = 1.1
      vlx1 = 0.0
      vly1 = 0.0

! Set up a pure left-going shock between states 1 and 2. I.e. the shock
! propagates from Quadrant I into Quadrant II.

      prs2 = 0.35
      rho2 = rho1 * ( prs2 / prs1 + (gamma - 1.0)/(gamma + 1.0) )
     1     / ( 1.0 + ( (gamma - 1.0) * prs2) / ( (gamma + 1.0) * prs1))
      vlx2 = vlx1 + SQRT ( ( prs2 - prs1) * ( rho2 - rho1)
     1                   / (rho2 * rho1) ) 
      vly2 = vly1

! There is also a pure shock that propagates in the y-direction from
! Quadrant III into Quadrant II.
! The same state in Quadrant III also causes a pure shock that propagates
! in the x-direction from Quadrant III into Quadrant IV.

      prs3 = prs1
      rho3 = rho1
      vlx3 = vlx2
      vly3 = vlx2

      prs4 = prs2
      rho4 = rho2
      vlx4 = vlx1
      vly4 = vlx2

      END IF

! ----------------------------------------------------------------------

! Set magnetic fields to zero.

          bfxbdy ( :, :, :, 1) = 0.0
          bfybdy ( :, :, :, 1) = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

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

          u ( ixx, iyy, izz, 5, 1) = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)

          u ( ixx, iyy, izz, 6, 1) = 0.0
          u ( ixx, iyy, izz, 7, 1) = 0.0
          u ( ixx, iyy, izz, 8, 1) = 0.0

! ----------

          tempaa = SQRT ( vlxuse**2 + vlyuse**2 + vlzuse**2 )

          tempab = SQRT ( gamma * prsuse / rhouse )

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_TWO_D_RIEM_PROBLEM_XY

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_SHKTUBE_1D ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Sets up a 1d MHD shocktube problem. The problems are catalogued in
! Ryu & Jones (1995) ApJ and Dai & Woodward (1994).
! Problem can be done in x, y or z-direction but only one dimension
! can be dynamically active.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz, iproblem, nzones

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rhouse, prsuse, enguse, vlxuse, vlyuse, vlzuse,
     1     bfxuse, bfyuse, bfzuse, spuse,

     1     tstop, tempaa, tempab, tempac, tempad

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

      REAL bfxbar ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1     v ( ONEDFLUIDBEGIN : ONEDFLUIDEND, NFIELD + NFLUID)

! ----------------------------------------------------------------------

! Pick a 1d MHD test problem. "gamma" will be appropriately reset,
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

      CALL SET_RP ( iproblem, nzones, gamma, tstop, v, bfxbar)

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

          bfxuse = bfxbar ( ixx)
          bfyuse = v ( ixx, 6)
          bfzuse = v ( ixx, 7)

#if ( NFLUID >= 1)
          spuse = v ( ixx, 8)
#endif

          ELSE IF ( iy1 .GT. 1) THEN

          rhouse = v ( iyy, 1)

          vlxuse = v ( iyy, 4)
          vlyuse = v ( iyy, 2)
          vlzuse = v ( iyy, 3)

          prsuse = v ( iyy, 5)

          bfxuse = v ( iyy, 7)
          bfyuse = bfxbar ( iyy)
          bfzuse = v ( iyy, 6)

#if ( NFLUID >= 1)
          spuse = v ( iyy, 8)
#endif

          ELSE

          rhouse = v ( izz, 1)

          vlxuse = v ( izz, 3)
          vlyuse = v ( izz, 4)
          vlzuse = v ( izz, 2)

          prsuse = v ( izz, 5)

          bfxuse = v ( izz, 6)
          bfyuse = v ( izz, 7)
          bfzuse = bfxbar ( izz)

#if ( NFLUID >= 1)
          spuse = v ( izz, 8)
#endif

          END IF

          enguse = prsuse / (gamma - 1.0)
     1       + 0.5 * rhouse * ( vlxuse**2 + vlyuse**2 + vlzuse**2)
     1       + ( bfxuse**2 + bfyuse**2 + bfzuse**2) / (8.0 * pi)

! ----------

          u ( ixx, iyy, izz, 1, 1) = rhouse

          u ( ixx, iyy, izz, 2, 1) = rhouse * vlxuse
          u ( ixx, iyy, izz, 3, 1) = rhouse * vlyuse
          u ( ixx, iyy, izz, 4, 1) = rhouse * vlzuse

          u ( ixx, iyy, izz, 5, 1) = enguse

          u ( ixx, iyy, izz, 6, 1) = bfxuse
          u ( ixx, iyy, izz, 7, 1) = bfyuse
          u ( ixx, iyy, izz, 8, 1) = bfzuse

#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9:, 1) = rhouse * spuse
#endif

          bfxbdy ( ixx, iyy, izz, 1) = bfxuse
          bfxbdy ( ixx - 1, iyy, izz, 1) = bfxuse

          bfybdy ( ixx, iyy, izz, 1) = bfyuse
          bfybdy ( ixx, iyy - 1, izz, 1) = bfyuse

          bfzbdy ( ixx, iyy, izz, 1) = bfzuse
          bfzbdy ( ixx, iyy, izz - 1, 1) = bfzuse

! ----------

          tempaa = vlxuse**2 + vlyuse**2 + vlzuse**2
          tempab = bfxuse**2 + bfyuse**2 + bfzuse**2

          tempab = tempab / ( 4.0 * pi * rhouse)
     1           + gamma * prsuse / rhouse
          tempab = SQRT ( tempab)
          tempaa = SQRT ( tempaa)

          tempac = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                     dz_save ( izz))

          tempac = cfl_coef_advct * tempac / ( tempaa + tempab)

          dtcur = AMIN1 ( dtcur, tempac)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_SHKTUBE_1D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE SET_RP ( iproblem, nzones, gamma, tstop, v, bfxbar)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iproblem, nzones, izonecenter, ix

      REAL gamma, tstop, 
     1     bfxbar ( ONEDFLUIDBEGIN : ONEDFLUIDEND),
     1     v ( ONEDFLUIDBEGIN : ONEDFLUIDEND, NFIELD + NFLUID),

     1 rhore_l, prsre_l, vlxre_l, vlyre_l, vlzre_l, bfyre_l, bfzre_l,
     1 rhore_r, prsre_r, vlxre_r, vlyre_r, vlzre_r, bfyre_r, bfzre_r,
     1 bfxre_0, fourpi

! ----------------------------------------------------------------------

      fourpi = 16.0 * ATAN ( 1.0)

*****************
*
      IF ( iproblem .EQ. 1) THEN
*
* Brio and Wu shock tube
*
      gamma = 2.0
*
      izonecenter = nzones / 2
*
      tstop = 0.1
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 1.0 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 0.75 * SQRT ( fourpi)
*
*
      rhore_r = 0.125
      prsre_r = 0.1
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = -1.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 2) THEN
*
* Supersonic intersecting streams.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.08
*
*
      rhore_l = 1.0
      prsre_l = 20.0
*
      vlxre_l = 10.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 5.0
      bfzre_l = 0.0
*
*
      bfxre_0 = 5.0
*
*
      rhore_r = 1.0
      prsre_r = 1.0
*
      vlxre_r = -10.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 5.0
      bfzre_r = 0.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 3) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.03
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 5.0
      bfzre_l = 0.0
*
*
      bfxre_0 = 3.0
*
*
      rhore_r = 0.1
      prsre_r = 10.0
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 2.0
      bfzre_r = 0.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 4) THEN
*
* Non-coplanar problem with all 7 waves.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.2
*
*
      rhore_l = 1.08
      prsre_l = 0.95
*
      vlxre_l = 1.2
      vlyre_l = 0.01
      vlzre_l = 0.5
*
      bfyre_l = 3.6
      bfzre_l = 2.0
*
*
      bfxre_0 = 2.0
*
*
      rhore_r = 1.0
      prsre_r = 1.0
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 4.0
      bfzre_r = 2.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 5) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.035
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 6.0
      bfzre_l = 0.0
*
*
      bfxre_0 = 3.0
*
*
      rhore_r = 0.1
      prsre_r = 10.0
*
      vlxre_r = 0.0
      vlyre_r = 2.0
      vlzre_r = 1.0
*
      bfyre_r = 1.0
      bfzre_r = 0.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 6) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.01
*
*
      rhore_l = 0.1
      prsre_l = 0.4
*
      vlxre_l = 50.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = -1.0
      bfzre_l = -2.0
*
*
      bfxre_0 = 0.0
*
*
      rhore_r = 0.1
      prsre_r = 0.2
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 1.0
      bfzre_r = 2.0
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 7) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.1
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = -1.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 1.0 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 0.0 * SQRT ( fourpi)
*
*
      rhore_r = 1.0
      prsre_r = 1.0
*
      vlxre_r = 1.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 1.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 8) THEN
*
* Switch-on fast shock.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.15
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 1.0 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 1.0 * SQRT ( fourpi)
*
*
      rhore_r = 0.2
      prsre_r = 0.1
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 0.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 9) THEN
*
* Switch-off fast rarefaction.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.15
*
*
      rhore_l = 0.4
      prsre_l = 0.52467
*
      vlxre_l = -0.66991
      vlyre_l = 0.98263
      vlzre_l = 0.0
*
      bfyre_l = 0.0025293 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 1.3 * SQRT ( fourpi)
*
*
      rhore_r = 1.0
      prsre_r = 1.0
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 1.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 10) THEN
*
* Switch-off slow shock.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.15
*
*
      rhore_l = 0.65
      prsre_l = 0.5
*
      vlxre_l = 0.667
      vlyre_l = -0.257
      vlzre_l = 0.0
*
      bfyre_l = 0.55 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 0.75 * SQRT ( fourpi)
*
*
      rhore_r = 1.0
      prsre_r = 0.75
*
      vlxre_r = 0.4
      vlyre_r = -0.94
      vlzre_r = 0.0
*
      bfyre_r = 0.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 11) THEN
*
* Switch-on slow rarefaction.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.16
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 0.0 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 0.7 * SQRT ( fourpi)
*
*
      rhore_r = 0.3
      prsre_r = 0.2
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 1.0
*
      bfyre_r = 1.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*****************
*
      IF ( iproblem .EQ. 12) THEN
*
* Compound wave formed by fast rarefaction and Alfven wave.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.16
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 1.0 * SQRT ( fourpi)
      bfzre_l = 0.0 * SQRT ( fourpi)
*
*
      bfxre_0 = 1.3 * SQRT ( fourpi)
*
*
      rhore_r = 0.4
      prsre_r = 0.4
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = -1.0 * SQRT ( fourpi)
      bfzre_r = 0.0 * SQRT ( fourpi)
*
      END IF
*
*
*****************
*
      IF ( iproblem .EQ. 13) THEN
*
* Noh-like problem with colliding supersonic streams.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.04
*
*
      rhore_l = 0.15
      prsre_l = 0.28
*
      vlxre_l = 21.55
      vlyre_l = 1.0
      vlzre_l = 1.0
*
      bfyre_l = -2.0
      bfzre_l = -1.0
*
*
      bfxre_0 = 0.0
*
*
      rhore_r = 0.1
      prsre_r = 0.1
*
      vlxre_r = -26.45
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 2.0
      bfzre_r = 1.0
*
      END IF
*
*
*****************
*
      IF ( iproblem .EQ. 14) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 5
*
      tstop = 0.16
*
*
      rhore_l = 0.32
      prsre_l = 14.83
*
      vlxre_l = 0.0
      vlyre_l = 1.91
      vlzre_l = 0.9551
*
      bfyre_l = 2.239
      bfzre_l = 1.12
*
*
      bfxre_0 = 2.0
*
*
      rhore_r = 0.1
      prsre_r = 0.1
*
      vlxre_r = -10.41
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 2.0
      bfzre_r = 1.0
*
      END IF
*
*
*****************
*
      IF ( iproblem .EQ. 15) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 10
*
      tstop = 0.08
*
*
      rhore_l = 0.2
      prsre_l = 10.79
*
      vlxre_l = 4.589
      vlyre_l = -2.061
      vlzre_l = 3.317
*
      bfyre_l = 3.8
      bfzre_l = 3.8
*
*
      bfxre_0 = 2.0
*
*
      rhore_r = 0.1
      prsre_r = 0.2
*
      vlxre_r = -4.589
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 1.0
      bfzre_r = 1.0
*
      END IF
*
*
*****************
*
      IF ( iproblem .EQ. 16) THEN
*
* Noh-like problem with very fast colliding supersonic streams.
*
      gamma = 1.6666666666
*
      izonecenter = nzones / 2
*
      tstop = 0.03
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 36.87
      vlyre_l = -0.155
      vlzre_l = -0.0386
*
      bfyre_l = 4.0
      bfzre_l = 1.0
*
*
      bfxre_0 = 4.0
*
*
      rhore_r = 1.0
      prsre_r = 1.0
*
      vlxre_r = -36.87
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 4.0
      bfzre_r = 1.0
*
      END IF
*
*
*****************
*
      IF ( iproblem .EQ. 17) THEN
*
*
      gamma = 1.6666666666
*
      izonecenter = 6 * nzones / 10
*
      tstop = 0.4
*
*
      rhore_l = 0.8129
      prsre_l = 0.4809
*
      vlxre_l = 1.801
      vlyre_l = 0.3672
      vlzre_l = 0.1836
*
      bfyre_l = 1.7856
      bfzre_l = 0.8928
*
*
      bfxre_0 = 4.0
*
*
      rhore_r = 1.0
      prsre_r = 0.1
*
      vlxre_r = -1.7942
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 2.0
      bfzre_r = 1.0
*
      END IF
*
*****************
* Stationary Contact discontinuity.
*
      IF ( iproblem .EQ. 18) THEN
*
*
      gamma = 1.4
*
      izonecenter = nzones / 2
*
      tstop = 0.25
*
*
      rhore_l = 1.0
      prsre_l = 1.0
*
      vlxre_l = 0.0
      vlyre_l = 0.0
      vlzre_l = 0.0
*
      bfyre_l = 0.0
      bfzre_l = 0.0
*
*
      bfxre_0 = 1.0
*
*
      rhore_r = 0.1
      prsre_r = 1.0
*
      vlxre_r = 0.0
      vlyre_r = 0.0
      vlzre_r = 0.0
*
      bfyre_r = 0.0
      bfzre_r = 0.0
*
      END IF
*
*****************
* Stationary right-going Alfven wave.
*
      IF ( iproblem .EQ. 19) THEN
*
*
      gamma = 1.4
*
      izonecenter = nzones / 2
*
      tstop = 0.25
*
*
      rhore_l = 1.0 / fourpi
      prsre_l = 1.0
*
      vlxre_l = - 1.0
      vlyre_l = 1.0
      vlzre_l = - 1.0
*
      bfyre_l = - 1.0
      bfzre_l = 1.0
*
*
      bfxre_0 = 1.0
*
*
      rhore_r = 1.0 / fourpi
      prsre_r = 1.0
*
      vlxre_r = - 1.0
      vlyre_r = - 1.0
      vlzre_r = - 1.0
*
      bfyre_r = 1.0
      bfzre_r = 1.0
*
      END IF
*
*****************

! Now set up all the primitive variables.

      bfxbar ( :) = bfxre_0

      DO ix = - 3, nzones + 4

      IF ( ix .LE. izonecenter) THEN

      v ( ix, 1) = rhore_l
      v ( ix, 2) = vlxre_l
      v ( ix, 3) = vlyre_l
      v ( ix, 4) = vlzre_l
      v ( ix, 5) = prsre_l
      v ( ix, 6) = bfyre_l
      v ( ix, 7) = bfzre_l

#if ( NFLUID >= 1)
      v ( ix, 8:) = 1.0
#endif

      ELSE

      v ( ix, 1) = rhore_r
      v ( ix, 2) = vlxre_r
      v ( ix, 3) = vlyre_r
      v ( ix, 4) = vlzre_r
      v ( ix, 5) = prsre_r
      v ( ix, 6) = bfyre_r
      v ( ix, 7) = bfzre_r

#if ( NFLUID >= 1)
      v ( ix, 8:) = 0.0
#endif

      END IF

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma, tstop, bfxre_0 = ", gamma, tstop, bfxre_0
      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"rhore_l, prsre_l, vlxre_l, vlyre_l, vlzre_l, bfyre_l, bfzre_l"
      WRITE ( 6, *)
     1 rhore_l, prsre_l, vlxre_l, vlyre_l, vlzre_l, bfyre_l, bfzre_l
      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"rhore_r, prsre_r, vlxre_r, vlyre_r, vlzre_r, bfyre_r, bfzre_r"
      WRITE ( 6, *)
     1 rhore_r, prsre_r, vlxre_r, vlyre_r, vlzre_r, bfyre_r, bfzre_r

      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE SET_RP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------




