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
! "bfx1/2 = 100.0" requires pressure-positivity fix, which is not there
! in this code. So it is best to simplify with "bfx1/2 = 20.0".

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

      bfx1 = 10.0 / SQRT ( 3.0)
      bfx2 = bfx1

      bfy1 = 10.0 / SQRT ( 3.0)
      bfy2 = bfy1

      bfz1 = 10.0 / SQRT ( 3.0)
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
          tempab = cfl_coef * tempad / tempac
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
          tempab = cfl_coef * tempad / tempac
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
          tempab = cfl_coef * tempad / tempac
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
          tempab = cfl_coef * tempad / tempac
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
          tempab = cfl_coef * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ROTOR_2D_YZ

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_DECAYING_TURBULENCE_3D

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

! We assume that the second, third and fourth components are loaded with
! a velocity that scales to an rms value of 1. We also initialize
! a density of unity to keep it all very simple. The magnetic field can then
! be set to some constant value. Pressure is also "1/gamma" so that
! sound speed is unity to start with.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     rho1, prs1, vlx1, vly1,
     1     vlz1, bfx1, bfy1, bfz1, v_mach, va_mach,
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

! Stop this problem at a time of "timestop = 0.25".


      v_mach = 2.0  ! Set up a Mach 2 turbulence. (rms vel / c_s)

      va_mach = 2.0  ! Set up an Alfvenic Mach number of 2. (alfven speed / c_s)

      rho1 = 1.0

      prs1 = 1.0 / gamma

      vlx1 = 0.0

      vly1 = 0.0

      vlz1 = 0.0

      bfx1 = va_mach * SQRT ( 16.0 * ATAN ( 1.0) ) ! This gives V_A = 2 * c_s

      bfy1 = 0.0

      bfz1 = 0.0

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

          u ( ixx, iyy, izz, 1, 1) = rho1

          u ( ixx, iyy, izz, 2, 1) = v_mach * rho1
     1                             * u ( ixx, iyy, izz, 2, 1) 
          u ( ixx, iyy, izz, 3, 1) = v_mach * rho1
     1                             * u ( ixx, iyy, izz, 3, 1)
          u ( ixx, iyy, izz, 4, 1) = v_mach * rho1
     1                             * u ( ixx, iyy, izz, 4, 1)

          vlx1 = u ( ixx, iyy, izz, 2, 1) / rho1
          vly1 = u ( ixx, iyy, izz, 3, 1) / rho1
          vlz1 = u ( ixx, iyy, izz, 4, 1) / rho1

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
          tempab = cfl_coef * tempad / tempac
          dtcur = AMIN1 ( dtcur, tempab)

#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) = rho1
#endif

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_DECAYING_TURBULENCE_3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_FORCING ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,

     1        dtcur, vrms_target, temp_target, temp_time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr,

     1        accngr_x, accngr_y, accngr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Uses the accelerations in "accngr_x, accngr_y, accngr_z" to figure
! out how much of a push is required to bring the rms velocity up to
! "vrms_target". If "vrms_target <= 0" it is a decaying turbulence so
! we exit right away.

! We also want the "temperature" i.e. (pressure/density) to reach a
! value of "temp_target" in an amount of time "temp_time".

! Then provides the push to the velocities (not momenta). Similarly,
! temperature is given a small change to bring it up to "temp_target"

! Finally, new momenta and total energy are rebuilt based on the push.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe

      REAL dtcur, vrms_target, temp_target, temp_time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1      bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr,

     1        accngr_x, accngr_y, accngr_z

! ----------

! Arrays that specify geometry.

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

      INTEGER ixx, iyy, izz, i_comp,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2

      REAL ob_8pi, tempaa, tempab, tempac, tempad,
     1     v_dot_v, v_dot_a, a_dot_a, accn_factor

! ----------------------------------------------------------------------

      IF ( vrms_target <= 0.0) RETURN

      ob_8pi = 1.0 / ( 8.0 * pi)

! ----------------------------------------------------------------------

! Build primitive variables.

! Build accumulation variables for "v * v", "v * a", "a * a".

      v_dot_v = 0.0
      v_dot_a = 0.0
      a_dot_a = 0.0

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, tempaa, tempab, tempac)                    &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ob_8pi, gamma, rhofloor, prsfloor,                            &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy,                                       &
!$OMP&    temperaturegr, gamma_eos_gr, gamma_soundspeed_gr,             &
!$OMP&    mean_mol_wt_gr, accngr_x, accngr_y, accngr_z)                 &
!$OMP&    REDUCTION ( +: v_dot_v)                                       &
!$OMP&    REDUCTION ( +: v_dot_a)                                       &
!$OMP&    REDUCTION ( +: a_dot_a)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ------------------------------
! Build primitive variables.

          u ( ixx, iyy, izz, 1, 1) = AMAX1 ( rhofloor,
     1                               u ( ixx, iyy, izz, 1, 1))

          rhogr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 1, 1)

          tempaa = 1.0 / rhogr ( ixx, iyy, izz)

          vlxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 2, 1) * tempaa
          vlygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 3, 1) * tempaa
          vlzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 4, 1) * tempaa

          bfxgr ( ixx, iyy, izz) =0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                                  + bfxbdy ( ixxm1, iyy, izz, 1))
          bfygr ( ixx, iyy, izz) =0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                                  + bfybdy ( ixx, iyym1, izz, 1))
          bfzgr ( ixx, iyy, izz) =0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                                  + bfzbdy ( ixx, iyy, izzm1, 1))

          u ( ixx, iyy, izz, 6, 1) = bfxgr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 7, 1) = bfygr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 8, 1) = bfzgr ( ixx, iyy, izz)

          tempab = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          tempac = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          prsgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 5, 1)
     1                           - tempab - tempac
          prsgr ( ixx, iyy, izz) = AMAX1 ( prsfloor, ( gamma - 1.0)
     1                             * prsgr ( ixx, iyy, izz) )

          temperaturegr ( ixx, iyy, izz) = prsgr ( ixx, iyy, izz)
     1                                   * tempaa

! ------------------------------
! Build accumulation variables.

      v_dot_v = v_dot_v
     1        + vlxgr ( ixx, iyy, izz) * vlxgr ( ixx, iyy, izz)
     1        + vlygr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)
     1        + vlzgr ( ixx, iyy, izz) * vlzgr ( ixx, iyy, izz)

      v_dot_a = v_dot_a
     1        + vlxgr ( ixx, iyy, izz) * accngr_x ( ixx, iyy, izz)
     1        + vlygr ( ixx, iyy, izz) * accngr_y ( ixx, iyy, izz)
     1        + vlzgr ( ixx, iyy, izz) * accngr_z ( ixx, iyy, izz)

      a_dot_a = a_dot_a
     1        + accngr_x ( ixx, iyy, izz) * accngr_x ( ixx, iyy, izz)
     1        + accngr_y ( ixx, iyy, izz) * accngr_y ( ixx, iyy, izz)
     1        + accngr_z ( ixx, iyy, izz) * accngr_z ( ixx, iyy, izz)

! ------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------------------------------------

! Finally, normalize the accumulation variables and write them out.
! Will need to be done differently in parallel.

      tempaa = 1.0 / REAL ( ix1 * iy1 * iz1)

      v_dot_v = v_dot_v * tempaa
      v_dot_a = v_dot_a * tempaa
      a_dot_a = a_dot_a * tempaa

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"v_dot_v, v_dot_a, a_dot_a = ",
     1              v_dot_v, v_dot_a, a_dot_a
      WRITE ( 6, *)"vrms_target**2, temp_target = ",
     1              vrms_target**2, temp_target

! ------------------------------

! Try to increase the "v_dot_v" by ~ 1-5% in the next timestep if it dips below
! 95% of "vrms_target**2".

! This is the coefficient in front of the expression
! for "tempac = - 0.01 * AMAX1 ( v_dot_v, 0.1)". I.e. 0.01 means 1%; 0.05
! means 5%; etc.

! 0.1 is internal floor value for velocity in 
! "tempac = - 0.01 * AMAX1 ( v_dot_v, 0.1)".

! "accn_factor" is the factor by which "accngr_x,y,z" should be multiplied
! to get the intended velocity push.

      IF ( v_dot_v >= vrms_target**2) THEN

        tempaa = 0.0
        tempab = 0.0
        tempac = 0.0
        tempad = 0.0
        accn_factor = 0.0  ! No acceleration if target rms velocity is exceeded.

      ELSE IF ( v_dot_v < 0.95 * vrms_target**2) THEN

        tempac = - 0.01 * AMAX1 ( v_dot_v, 0.1)  ! 0.1 is floor value of rms-vel
        tempab = dtcur * v_dot_a
        tempaa = dtcur**2 * a_dot_a

        tempad = AMAX1 ( 0.0, tempab**2 - 4.0 * tempaa * tempac)
        tempad = SQRT ( tempad)

        accn_factor = 0.5 * ( - tempab + tempad) / tempaa

      ELSE

        tempac = v_dot_v - vrms_target**2  ! We are very close to target.
        tempab = dtcur * v_dot_a
        tempaa = dtcur**2 * a_dot_a

        tempad = AMAX1 ( 0.0, tempab**2 - 4.0 * tempaa * tempac)
        tempad = SQRT ( tempad)

        accn_factor = 0.5 * ( - tempab + tempad) / tempaa

      END IF

      WRITE ( 6, *)" "
      WRITE ( 6, *)"tempaa, tempab, tempac, tempad = ",
     1              tempaa, tempab, tempac, tempad
      WRITE ( 6, *)"Accelerating with accn_factor = ", accn_factor
      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"

! ----------------------------------------------------------------------

! Now we are ready to add in the accelerations to the velocities.
! We will also apply cooling/heating to the flow.

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, tempaa, tempab, tempac, tempad)            &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ob_8pi, gamma, rhofloor, prsfloor, accn_factor, dtcur,        &
!$OMP&    temp_target, temp_time,                                       &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy,                                       &
!$OMP&    temperaturegr, gamma_eos_gr, gamma_soundspeed_gr,             &
!$OMP&    mean_mol_wt_gr, accngr_x, accngr_y, accngr_z)                 &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ------------------------------
! Rebuild velocity. This is an explicit push with the accelerations.

          vlxgr ( ixx, iyy, izz) = vlxgr ( ixx, iyy, izz)
     1          + accn_factor * dtcur * accngr_x ( ixx, iyy, izz)

          vlygr ( ixx, iyy, izz) = vlygr ( ixx, iyy, izz)
     1          + accn_factor * dtcur * accngr_y ( ixx, iyy, izz)

          vlzgr ( ixx, iyy, izz) = vlzgr ( ixx, iyy, izz)
     1          + accn_factor * dtcur * accngr_z ( ixx, iyy, izz)

! ------------------------------
! Rebuild "temperature", which is solved in implicit fashion, and
! then the pressure.

          tempad = dtcur / temp_time

          temperaturegr ( ixx, iyy, izz) =
     1  ( temperaturegr ( ixx, iyy, izz) + tempad * temp_target)
     1 / ( 1.0 + tempad)

          prsgr ( ixx, iyy, izz) = temperaturegr ( ixx, iyy, izz)
     1                           * rhogr ( ixx, iyy, izz)
          prsgr ( ixx, iyy, izz) = AMAX1 ( prsfloor,
     1                                     prsgr ( ixx, iyy, izz) )

! ------------------------------
! Rebuild energy density and momenta.

          tempab = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          tempac = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          u ( ixx, iyy, izz, 5, 1) =
     1    prsgr ( ixx, iyy, izz) / (gamma - 1.0) + tempab + tempac

          u ( ixx, iyy, izz, 2, 1) = rhogr ( ixx, iyy, izz)
     1                             * vlxgr ( ixx, iyy, izz)

          u ( ixx, iyy, izz, 3, 1) = rhogr ( ixx, iyy, izz)
     1                             * vlygr ( ixx, iyy, izz)

          u ( ixx, iyy, izz, 4, 1) = rhogr ( ixx, iyy, izz)
     1                             * vlzgr ( ixx, iyy, izz)

! ------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADD_FORCING

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

! When this application is run with ADER with HLL Riemann solver, it helps
! to reduce the CFL to 0.25 or so. The Roe type Riemann solver runs fine
! with CFL of 0.45 in 2d. This is most likely because the HLL does not resolve
! Alfven waves while the Roe RS does.

! Poorer interpolation schemes along with lower-grade Riemann solvers
! cause very rapid decay in the z-velocity or z-magnetic field. Very good
! interpolation schemes along with high-grade Riemann solvers
! reduce this decay by several orders of magnitude!

! The final amplitude at 129 time units for the magnetic field is given below:

! 2nd order ADER-MC-HLL, CFL 0.25 : vlzmax = 3.2448E-4, bfzmax = 1.1123E-3
! 2nd order ADER-WENO-HLL, CFL 0.25 : vlzmax = 9.5630E-4, bfzmax = 3.4044E-3
! 2nd order ADER-WENO-Roe, CFL 0.45 : vlzmax = 6.4084E-2, bfzmax = 0.2259

! 3rd order ADER-WENO-HLL, CFL 0.25 : vlzmax = 5.9188E-3, bfzmax = 2.1087E-2
! 3rd order ADER-WENO-Roe, CFL 0.45 : vlzmax = 9.1566E-2, bfzmax = 0.3248
! 3rd order ADER-PPM-HLL, CFL 0.2 : vlzmax = , bfzmax = 
! 3rd order ADER-PPM-Roe, CFL 0.2 : vlzmax = , bfzmax = 

! 4th order ADER-WENO-HLL, CFL 0.25 : vlzmax = 0.1822, bfzmax = 0.6459
! 4th order ADER-WENO-Roe, CFL 0.45 : vlzmax = 0.1953, bfzmax = 0.6922

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
     1     rr, nx, ny, mx, my, x, y, z, kx, ky, phase, x1, x2, y1, y2

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
      mx = - ny
      my = nx

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
      dtcur = cfl_coef / dtcur

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

      SUBROUTINE INIT_VORTEX

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

! Sets up a magnetized vortex in 2d. Decide which problem by setting
! "iproblem = 1,...,6". The same choice should be reflected in the
! error routine.

! The Euler problem is catalogued in Balsara & Shu (2000) JCP.
! Use a domain of [-5,5]x[-5,5]x[-5,5] on a range of meshes
! going from 64x64x1 to 512x512x1 zones for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 1.4" suggested. Suggested stopping time is
! "timestop = 10.0".

! The MHD problem is catalogued in Balsara (2004) ApJS.
! Use a domain of [-5,5]x[-5,5]x[-5,5] on a range of meshes
! going from 64x64x1 to 512x512x1 zones for this
! problem along with periodic boundaries (boundary condition # 4).
! "gamma = 5.0/3.0" suggested. Suggested stopping time is
! "timestop = 10.0".

! For fourth order & higher, double domain and stopping time.

! The vortex propagates along the diagonal of
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

     1        ixx, iyy, izz, iproblem

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     time, bmag, fourpi, max_signal_speed,

     1     tempaa, tempab, tempac, tempad,
     1     VECPOTFN, x, y,
     1     xlo, xhi, ylo, yhi, maxdivb,
     1     bfxbdy_xp1, bfxbdy_xm1, bfybdy_yp1, bfybdy_ym1,
     1     bfzbdy_zp1, bfzbdy_zm1

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

      VECPOTFN ( bmag, x, y) = bmag * EXP ( 0.5*(1.0-x**2-y**2))

! ----------------------------------------------------------------------

! Variables that contain the quadrature point info. as well as the
! information about making nodal to modal transcription.

      INTEGER iquadpt, i_euler, ivar_lo, ivar

      REAL xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( 1:3, NUM_ZONECEN_QUADPTS) :: zonecen_quadpt,
     1      zc_qp

      REAL, DIMENSION ( 1:3, NUM_FACECEN_QUADPTS) :: xface_quadpt,
     1      yface_quadpt, zface_quadpt, xf_qp, yf_qp, zf_qp

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_QUADPTS) ::
     1      cons_qp, prim_qp, t

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_BASES) ::
     1      cons_modal

! ----------------------------------------------------------------------

! Select a problem number 1,...,6.
      IF ( iz1 .EQ. 1) iproblem = 1    ! 1 or 2 for MHD/Euler cases
      IF ( ix1 .EQ. 1) iproblem = 3    ! 3 or 4 for MHD/Euler cases
      IF ( iy1 .EQ. 1) iproblem = 5    ! 5 or 6 for MHD/Euler cases

! Define some temporary variables needed for this problem.

      fourpi = 4.0 * pi

      time = 0.0

      dtcur = 1.0e32
      max_signal_speed = 1.0e-32

! This gives the magnitude of the magnetic vortex assuming a unit density.
      bmag = SQRT ( 4.0 * pi * 1.0) / ( 2.0 * pi)

! ----------------------------------------------------------------------

! The zone-centered quadrature points in the unit cube are stored
! in "zonecen_quadpt".
! The x-face centered quadrature points in the unit cube are stored
! in "xface_quadpt".
! The y-face centered quadrature points in the unit cube are stored
! in "yface_quadpt".
! The z-face centered quadrature points in the unit cube are stored
! in "zface_quadpt".

! ------------------------------

#if ( NUMORDER == 2)

      zonecen_quadpt ( 1:3, 01) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.0 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      xface_quadpt ( 1:3, 2) = (/ 0.5,-0.5, 0.0 /)
      xface_quadpt ( 1:3, 3) = (/ 0.5, 0.0, 0.5 /)
      xface_quadpt ( 1:3, 4) = (/ 0.5, 0.0,-0.5 /)

      yface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 2) = (/-0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      yface_quadpt ( 1:3, 4) = (/ 0.0, 0.5,-0.5 /)

      zface_quadpt ( 1:3, 1) = (/ 0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 2) = (/-0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      zface_quadpt ( 1:3, 4) = (/ 0.0,-0.5, 0.5 /)

#endif
! End "NUMORDER == 2".

! ------------------------------

#if ( NUMORDER == 3)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/ -0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/ -0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/ -0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 09) = (/ -0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/  0.5,  0.0,  0.0 /)
      xface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      xface_quadpt ( 1:3, 3) = (/  0.5, -0.5,  0.0 /)
      xface_quadpt ( 1:3, 4) = (/  0.5,  0.0,  0.5 /)
      xface_quadpt ( 1:3, 5) = (/  0.5,  0.0, -0.5 /)
      xface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      xface_quadpt ( 1:3, 7) = (/  0.5, -0.5,  0.5 /)
      xface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      xface_quadpt ( 1:3, 9) = (/  0.5, -0.5, -0.5 /)

      yface_quadpt ( 1:3, 1) = (/  0.0,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 3) = (/ -0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 5) = (/  0.0,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 9) = (/ -0.5,  0.5, -0.5 /)

      zface_quadpt ( 1:3, 1) = (/  0.0,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 2) = (/  0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 3) = (/ -0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 5) = (/  0.0, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 8) = (/  0.5, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 9) = (/ -0.5, -0.5,  0.5 /)

#endif
! End "NUMORDER == 3".

! ------------------------------

#if ( NUMORDER == 4)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/ -0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/  0.0 ,  0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/  0.0 , -0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 09) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 ,  0.0 ,  0.25 /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.0 , -0.25 /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 14) = (/  0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 15) = (/  0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 16) = (/  0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 17) = (/  0.5 , -0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 18) = (/ -0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 19) = (/ -0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 20) = (/ -0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 21) = (/ -0.5 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/   0.5,   0.0,   0.0 /)
      xface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      xface_quadpt ( 1:3, 3) = (/   0.5,  0.25,   0.0 /)
      xface_quadpt ( 1:3, 4) = (/   0.5, -0.25,   0.0 /)
      xface_quadpt ( 1:3, 5) = (/   0.5,  -0.5,   0.0 /)
      xface_quadpt ( 1:3, 6) = (/   0.5,   0.0,   0.5 /)
      xface_quadpt ( 1:3, 7) = (/   0.5,   0.0,  0.25 /)
      xface_quadpt ( 1:3, 8) = (/   0.5,   0.0, -0.25 /)
      xface_quadpt ( 1:3, 9) = (/   0.5,   0.0,  -0.5 /)
      xface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      xface_quadpt ( 1:3,11) = (/   0.5,  -0.5,   0.5 /)
      xface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      xface_quadpt ( 1:3,13) = (/   0.5,  -0.5,  -0.5 /)

      yface_quadpt ( 1:3, 1) = (/   0.0,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 3) = (/  0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 4) = (/ -0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 5) = (/  -0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      yface_quadpt ( 1:3, 7) = (/   0.0,   0.5,  0.25 /)
      yface_quadpt ( 1:3, 8) = (/   0.0,   0.5, -0.25 /)
      yface_quadpt ( 1:3, 9) = (/   0.0,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,13) = (/  -0.5,   0.5,  -0.5 /)

      zface_quadpt ( 1:3, 1) = (/   0.0,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 2) = (/   0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 3) = (/  0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 4) = (/ -0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 5) = (/  -0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      zface_quadpt ( 1:3, 7) = (/   0.0,  0.25,   0.5 /)
      zface_quadpt ( 1:3, 8) = (/   0.0, -0.25,   0.5 /)
      zface_quadpt ( 1:3, 9) = (/   0.0,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,12) = (/   0.5,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,13) = (/  -0.5,  -0.5,   0.5 /)

#endif
! End "NUMORDER == 4".

! ----------------------------------------------------------------------

! Make magnetic fields in active directions using a vector potential approach.
! That ensures divergence-free magnetic fields.

      maxdivb = 0.0

! ------------------------------
! Set up B-field in the xy-direction.

      IF ( iproblem .EQ. 1) THEN

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = xb_save ( ixx - 1)
          xhi = xb_save ( ixx)
          ylo = yb_save ( iyy - 1)
          yhi = yb_save ( iyy)

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

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

      END IF
! End of "iproblem .EQ. 1"

! ------------------------------
! Set up B-field in the yz-direction.

      IF ( iproblem .EQ. 3) THEN

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = yb_save ( iyy - 1)
          xhi = yb_save ( iyy)
          ylo = zb_save ( izz - 1)
          yhi = zb_save ( izz)

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


          bfybdy ( ixx, iyy, izz, 1) = bfxbdy_xp1
          bfybdy ( ixx, iyy - 1, izz, 1) = bfxbdy_xm1

          bfzbdy ( ixx, iyy, izz, 1) = bfybdy_yp1
          bfzbdy ( ixx, iyy, izz - 1, 1) = bfybdy_ym1

          bfxbdy ( ixx, iyy, izz, 1) = bfzbdy_zp1
          bfxbdy ( ixx - 1, iyy, izz, 1) = bfzbdy_zm1

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

      END IF
! End of "iproblem .EQ. 3"

! ------------------------------
! Set up B-field in the zx-direction.

      IF ( iproblem .EQ. 5) THEN

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------

          xlo = zb_save ( izz - 1)
          xhi = zb_save ( izz)
          ylo = xb_save ( ixx - 1)
          yhi = xb_save ( ixx)

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


          bfzbdy ( ixx, iyy, izz, 1) = bfxbdy_xp1
          bfzbdy ( ixx, iyy, izz - 1, 1) = bfxbdy_xm1

          bfxbdy ( ixx, iyy, izz, 1) = bfybdy_yp1
          bfxbdy ( ixx - 1, iyy, izz, 1) = bfybdy_ym1

          bfybdy ( ixx, iyy, izz, 1) = bfzbdy_zp1
          bfybdy ( ixx, iyy - 1, izz, 1) = bfzbdy_zm1

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

      END IF
! End of "iproblem .EQ. 5"

! ------------------------------

! Zero out magnetic fields for those problems where they are not needed.

      IF ( ( iproblem .EQ. 2) .OR. ( iproblem .EQ. 4) .OR.
     1     ( iproblem .EQ. 6)) THEN

          bfxbdy ( :, :, :, 1) = 0.0
          bfybdy ( :, :, :, 1) = 0.0
          bfzbdy ( :, :, :, 1) = 0.0

      END IF

! ----------------------------------------------------------------------

! Make all the primitive (and therefrom conserved) variables at the
! zone-centered quadrature points.
! The zone-centered quadrature points in the computational domain are stored
! in "zc_qp".
! The primitive variables at the quadrature points are stored in "prim_qp".
! The conserved variables at the quadrature points are stored in "cons_qp".
! The conserved variables in modal basis are stored in "cons_modal".

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------------------------------------------------------------------

          xc = xc_save ( ixx)
          dx = dx_save ( ixx)
          yc = yc_save ( iyy)
          dy = dy_save ( iyy)
          zc = zc_save ( izz)
          dz = dz_save ( izz)

          zc_qp ( 1, :) = zonecen_quadpt ( 1, :) * dx + xc
          zc_qp ( 2, :) = zonecen_quadpt ( 2, :) * dy + yc
          zc_qp ( 3, :) = zonecen_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! --------------------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, zc_qp ( 1, iquadpt), zc_qp ( 2, iquadpt),
     1    zc_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! --------------------

! Only the MHD fluid is used for constraining the timestep. If the
! Euler fluids can also constrain the timestep then include them here.

! This is the flow speed.
          tempaa = SQRT ( prim_qp ( 2, iquadpt)**2
     1                  + prim_qp ( 3, iquadpt)**2
     1                  + prim_qp ( 4, iquadpt)**2)
          tempab = gamma * prim_qp ( 5, iquadpt)
     1           + (1.0/fourpi) * ( prim_qp ( 6, iquadpt)**2
     1                            + prim_qp ( 7, iquadpt)**2
     1                            + prim_qp ( 8, iquadpt)**2)
! This is the magnetosonic signal speed.
          tempab = SQRT ( tempab / prim_qp ( 1, iquadpt))
! This is the flow + magnetosonic signal speed.
          tempac = tempaa + tempab

          tempad = cfl_coef * AMIN1 ( dx, dy, dz) / tempac
          dtcur = AMIN1 ( dtcur, tempad)

          max_signal_speed = AMAX1 ( max_signal_speed, tempac)

          cons_qp ( 1, iquadpt) = prim_qp ( 1, iquadpt)

          cons_qp ( 2, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 2, iquadpt)
          cons_qp ( 3, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 3, iquadpt)
          cons_qp ( 4, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 4, iquadpt)

          cons_qp ( 5, iquadpt) = prim_qp ( 5, iquadpt) / (gamma - 1.0)
     1     + 0.5 * prim_qp ( 1, iquadpt) * ( prim_qp ( 2, iquadpt)**2
     1           + prim_qp ( 3, iquadpt)**2 + prim_qp ( 4, iquadpt)**2)
     1     + ( prim_qp ( 6, iquadpt)**2 + prim_qp ( 7, iquadpt)**2
     1       + prim_qp ( 8, iquadpt)**2) / ( 8.0 * pi)

          cons_qp ( 6, iquadpt) = prim_qp ( 6, iquadpt)
          cons_qp ( 7, iquadpt) = prim_qp ( 7, iquadpt)
          cons_qp ( 8, iquadpt) = prim_qp ( 8, iquadpt)

#if ( NFLUID >= 1)
          cons_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1             prim_qp ( 1, iquadpt)
     1           * prim_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
#endif

          END DO
! End of "DO iquadpt ="

! ------------------------------

! Now that the conserved variables are available at all the quadrature
! points, we make the zone-centered conserved variables.

! For good measure, we make all their moments too,
! though they may not be needed.

! ----------

#if ( NUMORDER == 2)

          cons_modal ( :, 1) = ( cons_qp ( :, 01) + cons_qp ( :, 02)
     1                         + cons_qp ( :, 03) + cons_qp ( :, 04)
     1                         + cons_qp ( :, 05) + cons_qp ( :, 06))
     1                         / 6.0

          cons_modal ( :, 2) = cons_qp ( :, 01) - cons_qp ( :, 02)
          cons_modal ( :, 3) = cons_qp ( :, 03) - cons_qp ( :, 04)
          cons_modal ( :, 4) = cons_qp ( :, 05) - cons_qp ( :, 06)

#endif
! End of "NUMORDER == 2"
! ----------

#if ( NUMORDER == 3)

      t ( :, 01) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 04) + cons_qp ( :, 05)

      t ( :, 02) = cons_qp ( :, 06) + cons_qp ( :, 07)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 03) = cons_qp ( :, 10) + cons_qp ( :, 11)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 04) = cons_qp ( :, 02) + cons_qp ( :, 06)
     1           + cons_qp ( :, 10) + cons_qp ( :, 12)

      t ( :, 05) = cons_qp ( :, 03) + cons_qp ( :, 07)
     1           + cons_qp ( :, 11) + cons_qp ( :, 13)

      t ( :, 06) = cons_qp ( :, 04) + cons_qp ( :, 05)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 07) = cons_qp ( :, 04) + cons_qp ( :, 08)
     1           + cons_qp ( :, 10) + cons_qp ( :, 11)

      t ( :, 08) = cons_qp ( :, 05) + cons_qp ( :, 09)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 09) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 02) = 0.25 * ( t ( :, 01) - t ( :, 02) )

      cons_modal ( :, 03) = 0.25 * ( t ( :, 04) - t ( :, 05) )

      cons_modal ( :, 04) = 0.25 * ( t ( :, 07) - t ( :, 08) )

      cons_modal ( :, 05) = 0.5 * ( t ( :, 01) + t ( :, 02)
     1 - t ( :, 03) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 06) = 0.5 * ( t ( :, 04) + t ( :, 05)
     1 - t ( :, 06) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 07) = 0.5 * ( t ( :, 07) + t ( :, 08)
     1 - t ( :, 09) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 08) = cons_qp ( :, 02) - cons_qp ( :, 03)
     1 - cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 09) = cons_qp ( :, 10) - cons_qp ( :, 11)
     1 - cons_qp ( :, 12) + cons_qp ( :, 13)

      cons_modal ( :, 10) = cons_qp ( :, 04) - cons_qp ( :, 05)
     1 - cons_qp ( :, 08) + cons_qp ( :, 09)

      cons_modal ( :, 01) = cons_qp ( :, 01) + ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 12.0

#endif
! End of "NUMORDER == 3".
! ----------

#if ( NUMORDER == 4)

! Note that in this formulation "cons_modal ( :, 14:19)" require
! space-time nodal points. They, therefore, can't be built here.

      cons_modal ( :, 05) = 2.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 05) )

      cons_modal ( :, 06) = 2.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 09) )

      cons_modal ( :, 07) = 2.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 13) )

      cons_modal ( :, 08) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 15)
     1 + cons_qp ( :, 16) - cons_qp ( :, 17) - cons_qp ( :, 18)
     1 + cons_qp ( :, 19) - cons_qp ( :, 20) + cons_qp ( :, 21) )

      cons_modal ( :, 09) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 16)
     1 + cons_qp ( :, 18) - cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 17) - cons_qp ( :, 19) + cons_qp ( :, 21) )

      cons_modal ( :, 10) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 18)
     1 + cons_qp ( :, 15) - cons_qp ( :, 19) - cons_qp ( :, 16)
     1 + cons_qp ( :, 20) - cons_qp ( :, 17) + cons_qp ( :, 21) )

      cons_modal ( :, 11) = 16.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 03) + 2.0 * cons_qp ( :, 04)
     1 - cons_qp ( :, 05) ) / 3.0

      cons_modal ( :, 12) = 16.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 07) + 2.0 * cons_qp ( :, 08)
     1 - cons_qp ( :, 09) ) / 3.0

      cons_modal ( :, 13) = 16.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 11) + 2.0 * cons_qp ( :, 12)
     1 - cons_qp ( :, 13) ) / 3.0

      t ( :, 01) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 16) + cons_qp ( :, 17)

      t ( :, 02) = cons_qp ( :, 18) + cons_qp ( :, 19)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      t ( :, 03) = cons_qp ( :, 14) + cons_qp ( :, 16)
     1 + cons_qp ( :, 18) + cons_qp ( :, 20)

      t ( :, 04) = cons_qp ( :, 15) + cons_qp ( :, 17)
     1 + cons_qp ( :, 19) + cons_qp ( :, 21)

      t ( :, 05) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19)

      t ( :, 06) = cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      cons_modal ( :, 02) = 2.0 * ( cons_qp ( :, 02)
     1 - cons_qp ( :, 05) ) / 3.0 - 0.1 * cons_modal ( :, 11)
     1 + 0.25 * ( t ( :, 01) - t ( :, 02) ) / 3.0

      cons_modal ( :, 03) = 2.0 * ( cons_qp ( :, 06)
     1 - cons_qp ( :, 09) ) / 3.0 - 0.1 * cons_modal ( :, 12)
     1 + 0.25 * ( t ( :, 03) - t ( :, 04) ) / 3.0

      cons_modal ( :, 04) = 2.0 * ( cons_qp ( :, 10)
     1 - cons_qp ( :, 13) ) / 3.0 - 0.1 * cons_modal ( :, 13)
     1 + 0.25 * ( t ( :, 05) - t ( :, 06) ) / 3.0

      cons_modal ( :, 20) = cons_qp ( :, 14) - cons_qp ( :, 18)
     1 - cons_qp ( :, 16) + cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 19) + cons_qp ( :, 17) - cons_qp ( :, 21)

      cons_modal ( :, 01) = 0.125 * ( cons_qp ( :, 14)
     1 + cons_qp ( :, 15) + cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19) + cons_qp ( :, 20)
     1 + cons_qp ( :, 21) ) - 0.5 * ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 3.0

#endif
! End of "NUMORDER == 4".
! ----------
! ------------------------------

! Now store "cons_modal ( :, 1)" in the code's conserved variable meshes.

          u ( ixx, iyy, izz, :, 1) = cons_modal ( :, 1)

! ----------------------------------------------------------------------

! Make all the x-magnetic field variables at the top x-face-centered
! quadrature points.
! The top x-face-centered quadrature points in the computational domain
! are stored in "xf_qp".
! The x-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 6, :)".

! --------------------

          xf_qp ( 1, :) = xface_quadpt ( 1, :) * dx + xc
          xf_qp ( 2, :) = xface_quadpt ( 2, :) * dy + yc
          xf_qp ( 3, :) = xface_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, xf_qp ( 1, iquadpt), xf_qp ( 2, iquadpt),
     1    xf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 6, 01) = 0.25 * ( prim_qp ( 6, 01)
     1       + prim_qp ( 6, 02) + prim_qp ( 6, 03) + prim_qp ( 6, 04))

          cons_modal ( 6, 02) = prim_qp ( 6, 01) - prim_qp ( 6, 02)

          cons_modal ( 6, 03) = prim_qp ( 6, 03) - prim_qp ( 6, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 03)

          cons_modal ( 6, 03) = prim_qp ( 6, 04) - prim_qp ( 6, 05)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 03)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 04) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 06) = prim_qp ( 6, 06) - prim_qp ( 6, 07)
     1       - prim_qp ( 6, 08) + prim_qp ( 6, 09)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 06) + 2.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 06) = prim_qp ( 6, 10) - prim_qp ( 6, 11)
     1       - prim_qp ( 6, 12) + prim_qp ( 6, 13)

          cons_modal ( 6, 07) = ( 16.0 * prim_qp ( 6, 02)
     1       - 16.0 * prim_qp ( 6, 05) - 32.0 * prim_qp ( 6, 03)
     1       + 32.0 * prim_qp ( 6, 04) ) / 3.0

          cons_modal ( 6, 08) = ( 16.0 * prim_qp ( 6, 06)
     1       - 16.0 * prim_qp ( 6, 09) - 32.0 * prim_qp ( 6, 07)
     1       + 32.0 * prim_qp ( 6, 08) ) / 3.0

          cons_modal ( 6, 09) = 2.0 * prim_qp ( 6, 10)
     1       + 2.0 * prim_qp ( 6, 11) - 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 06)
     1       + 4.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 10) = 2.0 * prim_qp ( 6, 10)
     1       - 2.0 * prim_qp ( 6, 11) + 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 02)
     1       + 4.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 05)
     1       - 0.1 * cons_modal ( 6, 07) + cons_modal ( 6, 10) / 12.0

          cons_modal ( 6, 03) = prim_qp ( 6, 06) - prim_qp ( 6, 09)
     1       - 0.1 * cons_modal ( 6, 08) + cons_modal ( 6, 09) / 12.0

#endif
! ----------

!          bfxbdy ( ixx, iyy, izz, 1) = cons_modal ( 6, 1)

! ----------------------------------------------------------------------

! Make all the x-magnetic field variables at the bottom x-face-centered
! quadrature points.
! The bottom x-face-centered quadrature points in the computational domain
! are stored in "xf_qp".
! The x-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 6, :)".

! --------------------

          xf_qp ( 1, :) = - xface_quadpt ( 1, :) * dx + xc
          xf_qp ( 2, :) = xface_quadpt ( 2, :) * dy + yc
          xf_qp ( 3, :) = xface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, xf_qp ( 1, iquadpt), xf_qp ( 2, iquadpt),
     1    xf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 6, 01) = 0.25 * ( prim_qp ( 6, 01)
     1       + prim_qp ( 6, 02) + prim_qp ( 6, 03) + prim_qp ( 6, 04))

          cons_modal ( 6, 02) = prim_qp ( 6, 01) - prim_qp ( 6, 02)

          cons_modal ( 6, 03) = prim_qp ( 6, 03) - prim_qp ( 6, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 03)

          cons_modal ( 6, 03) = prim_qp ( 6, 04) - prim_qp ( 6, 05)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 03)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 04) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 06) = prim_qp ( 6, 06) - prim_qp ( 6, 07)
     1       - prim_qp ( 6, 08) + prim_qp ( 6, 09)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 06) + 2.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 06) = prim_qp ( 6, 10) - prim_qp ( 6, 11)
     1       - prim_qp ( 6, 12) + prim_qp ( 6, 13)

          cons_modal ( 6, 07) = ( 16.0 * prim_qp ( 6, 02)
     1       - 16.0 * prim_qp ( 6, 05) - 32.0 * prim_qp ( 6, 03)
     1       + 32.0 * prim_qp ( 6, 04) ) / 3.0

          cons_modal ( 6, 08) = ( 16.0 * prim_qp ( 6, 06)
     1       - 16.0 * prim_qp ( 6, 09) - 32.0 * prim_qp ( 6, 07)
     1       + 32.0 * prim_qp ( 6, 08) ) / 3.0

          cons_modal ( 6, 09) = 2.0 * prim_qp ( 6, 10)
     1       + 2.0 * prim_qp ( 6, 11) - 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 06)
     1       + 4.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 10) = 2.0 * prim_qp ( 6, 10)
     1       - 2.0 * prim_qp ( 6, 11) + 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 02)
     1       + 4.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 05)
     1       - 0.1 * cons_modal ( 6, 07) + cons_modal ( 6, 10) / 12.0

          cons_modal ( 6, 03) = prim_qp ( 6, 06) - prim_qp ( 6, 09)
     1       - 0.1 * cons_modal ( 6, 08) + cons_modal ( 6, 09) / 12.0

#endif
! ----------

!          bfxbdy ( ixx - 1, iyy, izz, 1) = cons_modal ( 6, 1)

! ----------------------------------------------------------------------

! Make all the y-magnetic field variables at the top y-face-centered
! quadrature points.
! The top y-face-centered quadrature points in the computational domain
! are stored in "yf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 7, :)".

! --------------------

          yf_qp ( 1, :) = yface_quadpt ( 1, :) * dx + xc
          yf_qp ( 2, :) = yface_quadpt ( 2, :) * dy + yc
          yf_qp ( 3, :) = yface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, yf_qp ( 1, iquadpt), yf_qp ( 2, iquadpt),
     1    yf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 7, 01) = 0.25 * ( prim_qp ( 7, 01)
     1       + prim_qp ( 7, 02) + prim_qp ( 7, 03) + prim_qp ( 7, 04))

          cons_modal ( 7, 02) = prim_qp ( 7, 01) - prim_qp ( 7, 02)

          cons_modal ( 7, 03) = prim_qp ( 7, 03) - prim_qp ( 7, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 03)

          cons_modal ( 7, 03) = prim_qp ( 7, 04) - prim_qp ( 7, 05)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 03)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 04) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 06) = prim_qp ( 7, 06) - prim_qp ( 7, 07)
     1       - prim_qp ( 7, 08) + prim_qp ( 7, 09)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 06) + 2.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 06) = prim_qp ( 7, 10) - prim_qp ( 7, 11)
     1       - prim_qp ( 7, 12) + prim_qp ( 7, 13)

          cons_modal ( 7, 07) = ( 16.0 * prim_qp ( 7, 02)
     1       - 16.0 * prim_qp ( 7, 05) - 32.0 * prim_qp ( 7, 03)
     1       + 32.0 * prim_qp ( 7, 04) ) / 3.0

          cons_modal ( 7, 08) = ( 16.0 * prim_qp ( 7, 06)
     1       - 16.0 * prim_qp ( 7, 09) - 32.0 * prim_qp ( 7, 07)
     1       + 32.0 * prim_qp ( 7, 08) ) / 3.0

          cons_modal ( 7, 09) = 2.0 * prim_qp ( 7, 10)
     1       + 2.0 * prim_qp ( 7, 11) - 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 06)
     1       + 4.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 10) = 2.0 * prim_qp ( 7, 10)
     1       - 2.0 * prim_qp ( 7, 11) + 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 02)
     1       + 4.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 05)
     1       - 0.1 * cons_modal ( 7, 07) + cons_modal ( 7, 10) / 12.0

          cons_modal ( 7, 03) = prim_qp ( 7, 06) - prim_qp ( 7, 09)
     1       - 0.1 * cons_modal ( 7, 08) + cons_modal ( 7, 09) / 12.0

#endif
! ----------

!          bfybdy ( ixx, iyy, izz, 1) = cons_modal ( 7, 1)

! ----------------------------------------------------------------------

! Make all the y-magnetic field variables at the bottom y-face-centered
! quadrature points.
! The bottom y-face-centered quadrature points in the computational domain
! are stored in "yf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 7, :)".

! --------------------

          yf_qp ( 1, :) = yface_quadpt ( 1, :) * dx + xc
          yf_qp ( 2, :) = - yface_quadpt ( 2, :) * dy + yc
          yf_qp ( 3, :) = yface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, yf_qp ( 1, iquadpt), yf_qp ( 2, iquadpt),
     1    yf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 7, 01) = 0.25 * ( prim_qp ( 7, 01)
     1       + prim_qp ( 7, 02) + prim_qp ( 7, 03) + prim_qp ( 7, 04))

          cons_modal ( 7, 02) = prim_qp ( 7, 01) - prim_qp ( 7, 02)

          cons_modal ( 7, 03) = prim_qp ( 7, 03) - prim_qp ( 7, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 03)

          cons_modal ( 7, 03) = prim_qp ( 7, 04) - prim_qp ( 7, 05)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 03)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 04) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 06) = prim_qp ( 7, 06) - prim_qp ( 7, 07)
     1       - prim_qp ( 7, 08) + prim_qp ( 7, 09)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 06) + 2.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 06) = prim_qp ( 7, 10) - prim_qp ( 7, 11)
     1       - prim_qp ( 7, 12) + prim_qp ( 7, 13)

          cons_modal ( 7, 07) = ( 16.0 * prim_qp ( 7, 02)
     1       - 16.0 * prim_qp ( 7, 05) - 32.0 * prim_qp ( 7, 03)
     1       + 32.0 * prim_qp ( 7, 04) ) / 3.0

          cons_modal ( 7, 08) = ( 16.0 * prim_qp ( 7, 06)
     1       - 16.0 * prim_qp ( 7, 09) - 32.0 * prim_qp ( 7, 07)
     1       + 32.0 * prim_qp ( 7, 08) ) / 3.0

          cons_modal ( 7, 09) = 2.0 * prim_qp ( 7, 10)
     1       + 2.0 * prim_qp ( 7, 11) - 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 06)
     1       + 4.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 10) = 2.0 * prim_qp ( 7, 10)
     1       - 2.0 * prim_qp ( 7, 11) + 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 02)
     1       + 4.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 05)
     1       - 0.1 * cons_modal ( 7, 07) + cons_modal ( 7, 10) / 12.0

          cons_modal ( 7, 03) = prim_qp ( 7, 06) - prim_qp ( 7, 09)
     1       - 0.1 * cons_modal ( 7, 08) + cons_modal ( 7, 09) / 12.0

#endif
! ----------

!          bfybdy ( ixx, iyy - 1, izz, 1) = cons_modal ( 7, 1)

! ----------------------------------------------------------------------

! Make all the z-magnetic field variables at the top z-face-centered
! quadrature points.
! The top z-face-centered quadrature points in the computational domain
! are stored in "zf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 8, :)".

! --------------------

          zf_qp ( 1, :) = zface_quadpt ( 1, :) * dx + xc
          zf_qp ( 2, :) = zface_quadpt ( 2, :) * dy + yc
          zf_qp ( 3, :) = zface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, zf_qp ( 1, iquadpt), zf_qp ( 2, iquadpt),
     1    zf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 8, 01) = 0.25 * ( prim_qp ( 8, 01)
     1       + prim_qp ( 8, 02) + prim_qp ( 8, 03) + prim_qp ( 8, 04))

          cons_modal ( 8, 02) = prim_qp ( 8, 01) - prim_qp ( 8, 02)

          cons_modal ( 8, 03) = prim_qp ( 8, 03) - prim_qp ( 8, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 03)

          cons_modal ( 8, 03) = prim_qp ( 8, 04) - prim_qp ( 8, 05)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 03)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 04) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 06) = prim_qp ( 8, 06) - prim_qp ( 8, 07)
     1       - prim_qp ( 8, 08) + prim_qp ( 8, 09)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 06) + 2.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 06) = prim_qp ( 8, 10) - prim_qp ( 8, 11)
     1       - prim_qp ( 8, 12) + prim_qp ( 8, 13)

          cons_modal ( 8, 07) = ( 16.0 * prim_qp ( 8, 02)
     1       - 16.0 * prim_qp ( 8, 05) - 32.0 * prim_qp ( 8, 03)
     1       + 32.0 * prim_qp ( 8, 04) ) / 3.0

          cons_modal ( 8, 08) = ( 16.0 * prim_qp ( 8, 06)
     1       - 16.0 * prim_qp ( 8, 09) - 32.0 * prim_qp ( 8, 07)
     1       + 32.0 * prim_qp ( 8, 08) ) / 3.0

          cons_modal ( 8, 09) = 2.0 * prim_qp ( 8, 10)
     1       + 2.0 * prim_qp ( 8, 11) - 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 06)
     1       + 4.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 10) = 2.0 * prim_qp ( 8, 10)
     1       - 2.0 * prim_qp ( 8, 11) + 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 02)
     1       + 4.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 05)
     1       - 0.1 * cons_modal ( 8, 07) + cons_modal ( 8, 10) / 12.0

          cons_modal ( 8, 03) = prim_qp ( 8, 06) - prim_qp ( 8, 09)
     1       - 0.1 * cons_modal ( 8, 08) + cons_modal ( 8, 09) / 12.0

#endif
! ----------

!          bfzbdy ( ixx, iyy, izz, 1) = cons_modal ( 8, 1)

! ----------------------------------------------------------------------

! Make all the z-magnetic field variables at the bottom z-face-centered
! quadrature points.
! The bottom z-face-centered quadrature points in the computational domain
! are stored in "zf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 8, :)".

! --------------------

          zf_qp ( 1, :) = zface_quadpt ( 1, :) * dx + xc
          zf_qp ( 2, :) = zface_quadpt ( 2, :) * dy + yc
          zf_qp ( 3, :) = - zface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, zf_qp ( 1, iquadpt), zf_qp ( 2, iquadpt),
     1    zf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 8, 01) = 0.25 * ( prim_qp ( 8, 01)
     1       + prim_qp ( 8, 02) + prim_qp ( 8, 03) + prim_qp ( 8, 04))

          cons_modal ( 8, 02) = prim_qp ( 8, 01) - prim_qp ( 8, 02)

          cons_modal ( 8, 03) = prim_qp ( 8, 03) - prim_qp ( 8, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 03)

          cons_modal ( 8, 03) = prim_qp ( 8, 04) - prim_qp ( 8, 05)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 03)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 04) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 06) = prim_qp ( 8, 06) - prim_qp ( 8, 07)
     1       - prim_qp ( 8, 08) + prim_qp ( 8, 09)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 06) + 2.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 06) = prim_qp ( 8, 10) - prim_qp ( 8, 11)
     1       - prim_qp ( 8, 12) + prim_qp ( 8, 13)

          cons_modal ( 8, 07) = ( 16.0 * prim_qp ( 8, 02)
     1       - 16.0 * prim_qp ( 8, 05) - 32.0 * prim_qp ( 8, 03)
     1       + 32.0 * prim_qp ( 8, 04) ) / 3.0

          cons_modal ( 8, 08) = ( 16.0 * prim_qp ( 8, 06)
     1       - 16.0 * prim_qp ( 8, 09) - 32.0 * prim_qp ( 8, 07)
     1       + 32.0 * prim_qp ( 8, 08) ) / 3.0

          cons_modal ( 8, 09) = 2.0 * prim_qp ( 8, 10)
     1       + 2.0 * prim_qp ( 8, 11) - 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 06)
     1       + 4.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 10) = 2.0 * prim_qp ( 8, 10)
     1       - 2.0 * prim_qp ( 8, 11) + 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 02)
     1       + 4.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 05)
     1       - 0.1 * cons_modal ( 8, 07) + cons_modal ( 8, 10) / 12.0

          cons_modal ( 8, 03) = prim_qp ( 8, 06) - prim_qp ( 8, 09)
     1       - 0.1 * cons_modal ( 8, 08) + cons_modal ( 8, 09) / 12.0

#endif
! ----------

!          bfzbdy ( ixx, iyy, izz - 1, 1) = cons_modal ( 8, 1)

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_VORTEX

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ERROR_VORTEX

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

! Evaluates the error in the magnetized or unmagnetized vortex in 2d.
! This subroutine should be placed in the timestep loop (ideally, it
! should be placed after the call to LIMIT_ALL_MOMENTS) & it will
! keep writing out the error. See further usage notes in the set-up code.
! Can analyze the accuracy of the problem by comparing the numerical
! solution with the analytical solution. Keep "iproblem" consistent.

! This is only done for serial setting. For parallel setting
! we would have to accumulate "L_1_error and L_inf_error" from all the
! processors.

! Everything is in place for P1Pm schemes also, i.e. we also
! evaluate the slopes.

! ----------------------------------------------------------------------



! hydrovortex_ADER_WENO3_CFL0.8

! zones     rhoerr_L1     accuracy_L1     rhoerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64     5.4996E-4                     9.5739E-3
! 128x128   8.0649E-5     2.77            1.3067E-3     2.87
! 256x256   1.0361E-5     2.96            1.7063E-4     2.94
! ----------------------------------------------------------------------



! hydrovortex_ADER_WENO4_CFL0.4

! zones     rhoerr_L1     accuracy_L1     rhoerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64     4.7883E-4                     3.2613E-2
! 128x128   2.5250E-5     4.25            1.6253E-3     4.33
! 256x256   1.0724E-6     4.56            7.9549E-5     4.35
! ----------------------------------------------------------------------



! magvortex_ADER_WENO3_CFL0.6

! zones     Bxerr_L1     accuracy_L1     Bxerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64     5.4912E-4                    1.1205E-2
! 128x128   6.6134E-5    3.05            1.7339E-3    2.69
! 256x256   8.0960E-6    3.03            2.8449E-4    2.61
! ----------------------------------------------------------------------



! magvortex_ADER_WENO4_CFL0.4

! zones     Bxerr_L1     accuracy_L1     Bxerr_Linf   accuracy_Linf
! ----------------------------------------------------------------------
! 64x64     4.6499E-4                    3.0675E-2
! 128x128   2.3300E-5    4.32            1.3958E-3    4.46
! 256x256   1.3189E-6    4.14            6.2268E-5    4.49
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz, iproblem

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     tempaa, tempab, tempac, fourpi

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

! Variables that contain the quadrature point info. as well as the
! information about making nodal to modal transcription.

      INTEGER iquadpt, i_euler, ivar_lo, ivar

      REAL xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( 1:3, NUM_ZONECEN_QUADPTS) :: zonecen_quadpt,
     1      zc_qp

      REAL, DIMENSION ( 1:3, NUM_FACECEN_QUADPTS) :: xface_quadpt,
     1      yface_quadpt, zface_quadpt, xf_qp, yf_qp, zf_qp

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_QUADPTS) ::
     1      cons_qp, prim_qp, t

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_BASES) ::
     1      cons_modal

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: diff_cons,
     1      L_1_error, L_inf_error

! ----------------------------------------------------------------------

! Select a problem number 1,...,6.
      IF ( iz1 .EQ. 1) iproblem = 1    ! 1 or 2 for MHD/Euler cases
      IF ( ix1 .EQ. 1) iproblem = 3    ! 3 or 4 for MHD/Euler cases
      IF ( iy1 .EQ. 1) iproblem = 5    ! 5 or 6 for MHD/Euler cases

! Define some temporary variables needed for this problem.

      fourpi = 4.0 * pi

      L_1_error ( :) = 0.0
      L_inf_error ( :) = 0.0

! ----------------------------------------------------------------------

! The zone-centered quadrature points in the unit cube are stored
! in "zonecen_quadpt".
! The x-face centered quadrature points in the unit cube are stored
! in "xface_quadpt".
! The y-face centered quadrature points in the unit cube are stored
! in "yface_quadpt".
! The z-face centered quadrature points in the unit cube are stored
! in "zface_quadpt".

! ------------------------------

#if ( NUMORDER == 2)

      zonecen_quadpt ( 1:3, 01) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.0 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      xface_quadpt ( 1:3, 2) = (/ 0.5,-0.5, 0.0 /)
      xface_quadpt ( 1:3, 3) = (/ 0.5, 0.0, 0.5 /)
      xface_quadpt ( 1:3, 4) = (/ 0.5, 0.0,-0.5 /)

      yface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 2) = (/-0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      yface_quadpt ( 1:3, 4) = (/ 0.0, 0.5,-0.5 /)

      zface_quadpt ( 1:3, 1) = (/ 0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 2) = (/-0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      zface_quadpt ( 1:3, 4) = (/ 0.0,-0.5, 0.5 /)

#endif
! End "NUMORDER == 2".

! ------------------------------

#if ( NUMORDER == 3)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/ -0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/ -0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/ -0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 09) = (/ -0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/  0.5,  0.0,  0.0 /)
      xface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      xface_quadpt ( 1:3, 3) = (/  0.5, -0.5,  0.0 /)
      xface_quadpt ( 1:3, 4) = (/  0.5,  0.0,  0.5 /)
      xface_quadpt ( 1:3, 5) = (/  0.5,  0.0, -0.5 /)
      xface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      xface_quadpt ( 1:3, 7) = (/  0.5, -0.5,  0.5 /)
      xface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      xface_quadpt ( 1:3, 9) = (/  0.5, -0.5, -0.5 /)

      yface_quadpt ( 1:3, 1) = (/  0.0,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 3) = (/ -0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 5) = (/  0.0,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 9) = (/ -0.5,  0.5, -0.5 /)

      zface_quadpt ( 1:3, 1) = (/  0.0,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 2) = (/  0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 3) = (/ -0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 5) = (/  0.0, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 8) = (/  0.5, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 9) = (/ -0.5, -0.5,  0.5 /)

#endif
! End "NUMORDER == 3".

! ------------------------------

#if ( NUMORDER == 4)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/ -0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/  0.0 ,  0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/  0.0 , -0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 09) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 ,  0.0 ,  0.25 /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.0 , -0.25 /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 14) = (/  0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 15) = (/  0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 16) = (/  0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 17) = (/  0.5 , -0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 18) = (/ -0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 19) = (/ -0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 20) = (/ -0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 21) = (/ -0.5 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/   0.5,   0.0,   0.0 /)
      xface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      xface_quadpt ( 1:3, 3) = (/   0.5,  0.25,   0.0 /)
      xface_quadpt ( 1:3, 4) = (/   0.5, -0.25,   0.0 /)
      xface_quadpt ( 1:3, 5) = (/   0.5,  -0.5,   0.0 /)
      xface_quadpt ( 1:3, 6) = (/   0.5,   0.0,   0.5 /)
      xface_quadpt ( 1:3, 7) = (/   0.5,   0.0,  0.25 /)
      xface_quadpt ( 1:3, 8) = (/   0.5,   0.0, -0.25 /)
      xface_quadpt ( 1:3, 9) = (/   0.5,   0.0,  -0.5 /)
      xface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      xface_quadpt ( 1:3,11) = (/   0.5,  -0.5,   0.5 /)
      xface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      xface_quadpt ( 1:3,13) = (/   0.5,  -0.5,  -0.5 /)

      yface_quadpt ( 1:3, 1) = (/   0.0,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 3) = (/  0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 4) = (/ -0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 5) = (/  -0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      yface_quadpt ( 1:3, 7) = (/   0.0,   0.5,  0.25 /)
      yface_quadpt ( 1:3, 8) = (/   0.0,   0.5, -0.25 /)
      yface_quadpt ( 1:3, 9) = (/   0.0,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,13) = (/  -0.5,   0.5,  -0.5 /)

      zface_quadpt ( 1:3, 1) = (/   0.0,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 2) = (/   0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 3) = (/  0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 4) = (/ -0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 5) = (/  -0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      zface_quadpt ( 1:3, 7) = (/   0.0,  0.25,   0.5 /)
      zface_quadpt ( 1:3, 8) = (/   0.0, -0.25,   0.5 /)
      zface_quadpt ( 1:3, 9) = (/   0.0,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,12) = (/   0.5,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,13) = (/  -0.5,  -0.5,   0.5 /)

#endif
! End "NUMORDER == 4".

! ----------------------------------------------------------------------

! Make all the primitive (and therefrom conserved) variables at the
! zone-centered quadrature points.
! The zone-centered quadrature points in the computational domain are stored
! in "zc_qp".
! The primitive variables at the quadrature points are stored in "prim_qp".
! The conserved variables at the quadrature points are stored in "cons_qp".
! The conserved variables in modal basis are stored in "cons_modal".

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------------------------------------------------------------------

          xc = xc_save ( ixx)
          dx = dx_save ( ixx)
          yc = yc_save ( iyy)
          dy = dy_save ( iyy)
          zc = zc_save ( izz)
          dz = dz_save ( izz)

          zc_qp ( 1, :) = zonecen_quadpt ( 1, :) * dx + xc
          zc_qp ( 2, :) = zonecen_quadpt ( 2, :) * dy + yc
          zc_qp ( 3, :) = zonecen_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! --------------------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, zc_qp ( 1, iquadpt), zc_qp ( 2, iquadpt),
     1    zc_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! --------------------

          cons_qp ( 1, iquadpt) = prim_qp ( 1, iquadpt)

          cons_qp ( 2, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 2, iquadpt)
          cons_qp ( 3, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 3, iquadpt)
          cons_qp ( 4, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 4, iquadpt)

          cons_qp ( 5, iquadpt) = prim_qp ( 5, iquadpt) / (gamma - 1.0)
     1     + 0.5 * prim_qp ( 1, iquadpt) * ( prim_qp ( 2, iquadpt)**2
     1           + prim_qp ( 3, iquadpt)**2 + prim_qp ( 4, iquadpt)**2)
     1     + ( prim_qp ( 6, iquadpt)**2 + prim_qp ( 7, iquadpt)**2
     1       + prim_qp ( 8, iquadpt)**2) / ( 8.0 * pi)

          cons_qp ( 6, iquadpt) = prim_qp ( 6, iquadpt)
          cons_qp ( 7, iquadpt) = prim_qp ( 7, iquadpt)
          cons_qp ( 8, iquadpt) = prim_qp ( 8, iquadpt)

#if ( NFLUID >= 1)
          cons_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1             prim_qp ( 1, iquadpt)
     1           * prim_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
#endif

          END DO
! End of "DO iquadpt ="

! ------------------------------

! Now that the conserved variables are available at all the quadrature
! points, we make the zone-centered conserved variables.

! For good measure, we make all their moments too,
! though they may not be needed.

! ----------

#if ( NUMORDER == 2)

          cons_modal ( :, 1) = ( cons_qp ( :, 01) + cons_qp ( :, 02)
     1                         + cons_qp ( :, 03) + cons_qp ( :, 04)
     1                         + cons_qp ( :, 05) + cons_qp ( :, 06))
     1                         / 6.0

          cons_modal ( :, 2) = cons_qp ( :, 01) - cons_qp ( :, 02)
          cons_modal ( :, 3) = cons_qp ( :, 03) - cons_qp ( :, 04)
          cons_modal ( :, 4) = cons_qp ( :, 05) - cons_qp ( :, 06)

#endif

! ----------

#if ( NUMORDER == 3)

      t ( :, 01) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 04) + cons_qp ( :, 05)

      t ( :, 02) = cons_qp ( :, 06) + cons_qp ( :, 07)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 03) = cons_qp ( :, 10) + cons_qp ( :, 11)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 04) = cons_qp ( :, 02) + cons_qp ( :, 06)
     1           + cons_qp ( :, 10) + cons_qp ( :, 12)

      t ( :, 05) = cons_qp ( :, 03) + cons_qp ( :, 07)
     1           + cons_qp ( :, 11) + cons_qp ( :, 13)

      t ( :, 06) = cons_qp ( :, 04) + cons_qp ( :, 05)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 07) = cons_qp ( :, 04) + cons_qp ( :, 08)
     1           + cons_qp ( :, 10) + cons_qp ( :, 11)

      t ( :, 08) = cons_qp ( :, 05) + cons_qp ( :, 09)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 09) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 02) = 0.25 * ( t ( :, 01) - t ( :, 02) )

      cons_modal ( :, 03) = 0.25 * ( t ( :, 04) - t ( :, 05) )

      cons_modal ( :, 04) = 0.25 * ( t ( :, 07) - t ( :, 08) )

      cons_modal ( :, 05) = 0.5 * ( t ( :, 01) + t ( :, 02)
     1 - t ( :, 03) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 06) = 0.5 * ( t ( :, 04) + t ( :, 05)
     1 - t ( :, 06) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 07) = 0.5 * ( t ( :, 07) + t ( :, 08)
     1 - t ( :, 09) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 08) = cons_qp ( :, 02) - cons_qp ( :, 03)
     1 - cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 09) = cons_qp ( :, 10) - cons_qp ( :, 11)
     1 - cons_qp ( :, 12) + cons_qp ( :, 13)

      cons_modal ( :, 10) = cons_qp ( :, 04) - cons_qp ( :, 05)
     1 - cons_qp ( :, 08) + cons_qp ( :, 09)

      cons_modal ( :, 01) = cons_qp ( :, 01) + ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 12.0

#endif

! ----------

#if ( NUMORDER == 4)

! Note that in this formulation "cons_modal ( :, 14:19)" require
! space-time nodal points. They, therefore, can't be built here.

      cons_modal ( :, 05) = 2.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 05) )

      cons_modal ( :, 06) = 2.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 09) )

      cons_modal ( :, 07) = 2.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 13) )

      cons_modal ( :, 08) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 15)
     1 + cons_qp ( :, 16) - cons_qp ( :, 17) - cons_qp ( :, 18)
     1 + cons_qp ( :, 19) - cons_qp ( :, 20) + cons_qp ( :, 21) )

      cons_modal ( :, 09) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 16)
     1 + cons_qp ( :, 18) - cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 17) - cons_qp ( :, 19) + cons_qp ( :, 21) )

      cons_modal ( :, 10) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 18)
     1 + cons_qp ( :, 15) - cons_qp ( :, 19) - cons_qp ( :, 16)
     1 + cons_qp ( :, 20) - cons_qp ( :, 17) + cons_qp ( :, 21) )

      cons_modal ( :, 11) = 16.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 03) + 2.0 * cons_qp ( :, 04)
     1 - cons_qp ( :, 05) ) / 3.0

      cons_modal ( :, 12) = 16.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 07) + 2.0 * cons_qp ( :, 08)
     1 - cons_qp ( :, 09) ) / 3.0

      cons_modal ( :, 13) = 16.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 11) + 2.0 * cons_qp ( :, 12)
     1 - cons_qp ( :, 13) ) / 3.0

      t ( :, 01) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 16) + cons_qp ( :, 17)

      t ( :, 02) = cons_qp ( :, 18) + cons_qp ( :, 19)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      t ( :, 03) = cons_qp ( :, 14) + cons_qp ( :, 16)
     1 + cons_qp ( :, 18) + cons_qp ( :, 20)

      t ( :, 04) = cons_qp ( :, 15) + cons_qp ( :, 17)
     1 + cons_qp ( :, 19) + cons_qp ( :, 21)

      t ( :, 05) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19)

      t ( :, 06) = cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      cons_modal ( :, 02) = 2.0 * ( cons_qp ( :, 02)
     1 - cons_qp ( :, 05) ) / 3.0 - 0.1 * cons_modal ( :, 11)
     1 + 0.25 * ( t ( :, 01) - t ( :, 02) ) / 3.0

      cons_modal ( :, 03) = 2.0 * ( cons_qp ( :, 06)
     1 - cons_qp ( :, 09) ) / 3.0 - 0.1 * cons_modal ( :, 12)
     1 + 0.25 * ( t ( :, 03) - t ( :, 04) ) / 3.0

      cons_modal ( :, 04) = 2.0 * ( cons_qp ( :, 10)
     1 - cons_qp ( :, 13) ) / 3.0 - 0.1 * cons_modal ( :, 13)
     1 + 0.25 * ( t ( :, 05) - t ( :, 06) ) / 3.0

      cons_modal ( :, 20) = cons_qp ( :, 14) - cons_qp ( :, 18)
     1 - cons_qp ( :, 16) + cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 19) + cons_qp ( :, 17) - cons_qp ( :, 21)

      cons_modal ( :, 01) = 0.125 * ( cons_qp ( :, 14)
     1 + cons_qp ( :, 15) + cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19) + cons_qp ( :, 20)
     1 + cons_qp ( :, 21) ) - 0.5 * ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 3.0

#endif

! ------------------------------

! Now store "diff_cons ( :)" between "cons_modal ( :, 1)" and the code's
! conserved variable mesh function "u". While these are zone-centered
! components, the magnetic field components will be overwritten.

          diff_cons ( :) = ABS ( u ( ixx, iyy, izz, :, 1)
     1                         - cons_modal ( :, 1) )

! ----------------------------------------------------------------------

! Make all the x-magnetic field variables at the top x-face-centered
! quadrature points.
! The top x-face-centered quadrature points in the computational domain
! are stored in "xf_qp".
! The x-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 6, :)".

! --------------------

          xf_qp ( 1, :) = xface_quadpt ( 1, :) * dx + xc
          xf_qp ( 2, :) = xface_quadpt ( 2, :) * dy + yc
          xf_qp ( 3, :) = xface_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, xf_qp ( 1, iquadpt), xf_qp ( 2, iquadpt),
     1    xf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 6, 01) = 0.25 * ( prim_qp ( 6, 01)
     1       + prim_qp ( 6, 02) + prim_qp ( 6, 03) + prim_qp ( 6, 04))

          cons_modal ( 6, 02) = prim_qp ( 6, 01) - prim_qp ( 6, 02)

          cons_modal ( 6, 03) = prim_qp ( 6, 03) - prim_qp ( 6, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 03)

          cons_modal ( 6, 03) = prim_qp ( 6, 04) - prim_qp ( 6, 05)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 03)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 04) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 06) = prim_qp ( 6, 06) - prim_qp ( 6, 07)
     1       - prim_qp ( 6, 08) + prim_qp ( 6, 09)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 06) + 2.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 06) = prim_qp ( 6, 10) - prim_qp ( 6, 11)
     1       - prim_qp ( 6, 12) + prim_qp ( 6, 13)

          cons_modal ( 6, 07) = ( 16.0 * prim_qp ( 6, 02)
     1       - 16.0 * prim_qp ( 6, 05) - 32.0 * prim_qp ( 6, 03)
     1       + 32.0 * prim_qp ( 6, 04) ) / 3.0

          cons_modal ( 6, 08) = ( 16.0 * prim_qp ( 6, 06)
     1       - 16.0 * prim_qp ( 6, 09) - 32.0 * prim_qp ( 6, 07)
     1       + 32.0 * prim_qp ( 6, 08) ) / 3.0

          cons_modal ( 6, 09) = 2.0 * prim_qp ( 6, 10)
     1       + 2.0 * prim_qp ( 6, 11) - 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 06)
     1       + 4.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 10) = 2.0 * prim_qp ( 6, 10)
     1       - 2.0 * prim_qp ( 6, 11) + 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 02)
     1       + 4.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 05)
     1       - 0.1 * cons_modal ( 6, 07) + cons_modal ( 6, 10) / 12.0

          cons_modal ( 6, 03) = prim_qp ( 6, 06) - prim_qp ( 6, 09)
     1       - 0.1 * cons_modal ( 6, 08) + cons_modal ( 6, 09) / 12.0

#endif
! ----------

          diff_cons ( 6) = ABS ( bfxbdy ( ixx, iyy, izz, 1)
     1                         - cons_modal ( 6, 1) )

! ----------------------------------------------------------------------

! Make all the y-magnetic field variables at the top y-face-centered
! quadrature points.
! The top y-face-centered quadrature points in the computational domain
! are stored in "yf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 7, :)".

! --------------------

          yf_qp ( 1, :) = yface_quadpt ( 1, :) * dx + xc
          yf_qp ( 2, :) = yface_quadpt ( 2, :) * dy + yc
          yf_qp ( 3, :) = yface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, yf_qp ( 1, iquadpt), yf_qp ( 2, iquadpt),
     1    yf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 7, 01) = 0.25 * ( prim_qp ( 7, 01)
     1       + prim_qp ( 7, 02) + prim_qp ( 7, 03) + prim_qp ( 7, 04))

          cons_modal ( 7, 02) = prim_qp ( 7, 01) - prim_qp ( 7, 02)

          cons_modal ( 7, 03) = prim_qp ( 7, 03) - prim_qp ( 7, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 03)

          cons_modal ( 7, 03) = prim_qp ( 7, 04) - prim_qp ( 7, 05)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 03)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 04) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 06) = prim_qp ( 7, 06) - prim_qp ( 7, 07)
     1       - prim_qp ( 7, 08) + prim_qp ( 7, 09)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 06) + 2.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 06) = prim_qp ( 7, 10) - prim_qp ( 7, 11)
     1       - prim_qp ( 7, 12) + prim_qp ( 7, 13)

          cons_modal ( 7, 07) = ( 16.0 * prim_qp ( 7, 02)
     1       - 16.0 * prim_qp ( 7, 05) - 32.0 * prim_qp ( 7, 03)
     1       + 32.0 * prim_qp ( 7, 04) ) / 3.0

          cons_modal ( 7, 08) = ( 16.0 * prim_qp ( 7, 06)
     1       - 16.0 * prim_qp ( 7, 09) - 32.0 * prim_qp ( 7, 07)
     1       + 32.0 * prim_qp ( 7, 08) ) / 3.0

          cons_modal ( 7, 09) = 2.0 * prim_qp ( 7, 10)
     1       + 2.0 * prim_qp ( 7, 11) - 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 06)
     1       + 4.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 10) = 2.0 * prim_qp ( 7, 10)
     1       - 2.0 * prim_qp ( 7, 11) + 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 02)
     1       + 4.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 05)
     1       - 0.1 * cons_modal ( 7, 07) + cons_modal ( 7, 10) / 12.0

          cons_modal ( 7, 03) = prim_qp ( 7, 06) - prim_qp ( 7, 09)
     1       - 0.1 * cons_modal ( 7, 08) + cons_modal ( 7, 09) / 12.0

#endif
! ----------

          diff_cons ( 7) = ABS ( bfybdy ( ixx, iyy, izz, 1)
     1                         - cons_modal ( 7, 1) )

! ----------------------------------------------------------------------

! Make all the z-magnetic field variables at the top z-face-centered
! quadrature points.
! The top z-face-centered quadrature points in the computational domain
! are stored in "zf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 8, :)".

! --------------------

          zf_qp ( 1, :) = zface_quadpt ( 1, :) * dx + xc
          zf_qp ( 2, :) = zface_quadpt ( 2, :) * dy + yc
          zf_qp ( 3, :) = zface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL VORTEX_GENERATOR_MAKE_PRIM
     1  ( iproblem, zf_qp ( 1, iquadpt), zf_qp ( 2, iquadpt),
     1    zf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 8, 01) = 0.25 * ( prim_qp ( 8, 01)
     1       + prim_qp ( 8, 02) + prim_qp ( 8, 03) + prim_qp ( 8, 04))

          cons_modal ( 8, 02) = prim_qp ( 8, 01) - prim_qp ( 8, 02)

          cons_modal ( 8, 03) = prim_qp ( 8, 03) - prim_qp ( 8, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 03)

          cons_modal ( 8, 03) = prim_qp ( 8, 04) - prim_qp ( 8, 05)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 03)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 04) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 06) = prim_qp ( 8, 06) - prim_qp ( 8, 07)
     1       - prim_qp ( 8, 08) + prim_qp ( 8, 09)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 06) + 2.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 06) = prim_qp ( 8, 10) - prim_qp ( 8, 11)
     1       - prim_qp ( 8, 12) + prim_qp ( 8, 13)

          cons_modal ( 8, 07) = ( 16.0 * prim_qp ( 8, 02)
     1       - 16.0 * prim_qp ( 8, 05) - 32.0 * prim_qp ( 8, 03)
     1       + 32.0 * prim_qp ( 8, 04) ) / 3.0

          cons_modal ( 8, 08) = ( 16.0 * prim_qp ( 8, 06)
     1       - 16.0 * prim_qp ( 8, 09) - 32.0 * prim_qp ( 8, 07)
     1       + 32.0 * prim_qp ( 8, 08) ) / 3.0

          cons_modal ( 8, 09) = 2.0 * prim_qp ( 8, 10)
     1       + 2.0 * prim_qp ( 8, 11) - 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 06)
     1       + 4.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 10) = 2.0 * prim_qp ( 8, 10)
     1       - 2.0 * prim_qp ( 8, 11) + 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 02)
     1       + 4.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 05)
     1       - 0.1 * cons_modal ( 8, 07) + cons_modal ( 8, 10) / 12.0

          cons_modal ( 8, 03) = prim_qp ( 8, 06) - prim_qp ( 8, 09)
     1       - 0.1 * cons_modal ( 8, 08) + cons_modal ( 8, 09) / 12.0

#endif
! ----------

          diff_cons ( 8) = ABS ( bfzbdy ( ixx, iyy, izz, 1)
     1                         - cons_modal ( 8, 1) )

! ----------------------------------------------------------------------

! Now use "diff_cons ( :)" to build the "L_1_error and L_inf_error" arrays.

          L_1_error ( :) = L_1_error ( :) + diff_cons ( :)
          L_inf_error ( :) = AMAX1 ( L_inf_error ( :), diff_cons ( :))

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

! Write out the accuracy information.

      L_1_error ( :) = L_1_error ( :) / ( ix1 * iy1 * iz1)

      WRITE ( 6, *)
     1"---------------------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Accuracy analysis for iproblem = ", iproblem,
     1             " at time = ", time
      WRITE ( 6, *)"# of zones : ix1, iy1, iz1 = ", ix1, iy1, iz1


      DO ixx = 1, NFIELD + 1 + NFLUID
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Error in variable no. = ", ixx
      WRITE ( 6, *)"L_1_error, L_inf_error = ",
     1              L_1_error ( ixx), L_inf_error ( ixx)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"---------------------------------------------------------------"

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ERROR_VORTEX

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE VORTEX_GENERATOR_MAKE_PRIM

     1   ( iproblem, x, y, z, time,
     1     xmin, xmax, ymin, ymax, zmin, zmax,
     1     gamma, prim)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! For higher order schemes, it is most convenient to just get a
! vector of primitive variables at any general space or time point 
! in the computational domain. This information can then be used to
! set the variables at zone-centers, faces etc.

! To have a div-free set up of magnetic field in multi-dimensions, it
! is still best to start with a vector potential approach.

! The present code sets up either a magnetized, isodensity vortex or
! an unmagnetized, isentropic vortex.
! The magnetized, isodensity vortex is from Balsara (2004) Ap.J.Supp,
! vol. 151, pg. 168.
! The unmagnetized, isentropic vortex is from Balsara & Shu (2000) JCP,
! vol. 160, pg. 424.

! The problem and its direction depends on "iproblem".
! "iproblem=1" magnetized vortex in xy-plane;
! "iproblem=2" unmagnetized, isentropic vortex in xy-plane;
! "iproblem=3" magnetized vortex in yz-plane;
! "iproblem=4" unmagnetized, isentropic vortex in yz-plane;
! "iproblem=5" magnetized vortex in xz-plane;
! "iproblem=6" unmagnetized, isentropic vortex in xz-plane;

! The variable names and meanings (and their input/output status) are given
! below:

! "iproblem" : Problem number that is being set up. (input)
! "x, y, z, time" : space-time position at which variables are desired. (input)
! "xmin, xmax, ymin, ymax, zmin, zmax" : Domain boundary in 3d. (input)
! "gamma" : gamma for MHD and Euler fluids. (input)

! "prim ( NFIELD + 1 + NFLUID)" : Vector of primitive variables. (output)
!                  ( Order : mhd density, mhd x,y,z-velocities, mhd pressure,
!                    bfx,y,z, mhd species fraction;
!                    euler density, euler x,y,z-velocities, euler pressure,
!                    euler species fraction; repeat for more euler families)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iproblem

      REAL x, y, z, time, gamma, prim ( NFIELD + 1 + NFLUID),
     1     xmin, xmax, ymin, ymax, zmin, zmax

! ----------------------------------------------------------------------

      REAL rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0, pi,
     1     kappa, mu, tt, xx, yy, zz, rr2, exp_rr2, tempaa, tempab

! ----------------------------------------------------------------------

      pi = 4.0 * ATAN ( 1.0)

      IF ( NFIELD + 1 + NFLUID .NE. 8) THEN
      WRITE ( 6, *)"Stopping in VORTEX_GENERATOR_MAKE_PRIM, no. vars= ",
     1              NFIELD + 1 + NFLUID
      STOP
      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 1) THEN

! This is the magnetized vortex problem in the xy-plane from Balsara (2004)
! Ap.J.Supp.

! ----------

      rho0 = 1.0
      prs0 = 1.0
      vlx0 = 1.0
      vly0 = 1.0
      vlz0 = 0.0
      bfx0 = 0.0
      bfy0 = 0.0
      bfz0 = 0.0

! Strength of the fluid vorticity.
      kappa = 1.0

! Strength of the curl of the magnetic field.
      mu = SQRT ( 4.0 * pi)

! ----------

! "tt" is the time. However, if the domain has shifted over by exactly
! one period, we set "tt = 0.0".
      tt = time
      IF ( ( ABS ( xmax - vlx0 * tt - xmin) .LT. 1.0e-14) .OR.
     1     ( ABS ( ymax - vly0 * tt - ymin) .LT. 1.0e-14) ) tt = 0.0

      xx = x - vlx0 * tt
      yy = y - vly0 * tt

      IF ( xx .LT. xmin) xx = xx + ( xmax - xmin)
      IF ( yy .LT. ymin) yy = yy + ( ymax - ymin)

      rr2 = xx**2 + yy**2
      exp_rr2 = EXP ( 0.5 * ( 1.0 - rr2) )

! ----------

      prim ( 1) = rho0

      prim ( 2) = vlx0 - yy * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 3) = vly0 + xx * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 4) = vlz0

      prim ( 5) = prs0 + ( exp_rr2**2 / ( 4.0 * pi**2) )
     1                 * ( mu**2 * (1.0 - rr2) / (8.0 * pi)
     1                   - 0.5 * kappa**2)

      prim ( 6) = bfx0 - yy * mu * exp_rr2 / ( 2.0 * pi)
      prim ( 7) = bfy0 + xx * mu * exp_rr2 / ( 2.0 * pi)
      prim ( 8) = bfz0

! ----------

      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 2) THEN

! This is the unmagnetized vortex problem in the xy-plane from
! Balsara & Shu (2000) JCP.

! ----------

      rho0 = 1.0
      prs0 = 1.0
      vlx0 = 1.0
      vly0 = 1.0
      vlz0 = 0.0
      bfx0 = 0.0
      bfy0 = 0.0
      bfz0 = 0.0

! Strength of the fluid vorticity. This "kappa" is the same as "epsilon"
! from Balsara & Shu (2000) JCP.
      kappa = 5.0

! ----------

! "tt" is the time. However, if the domain has shifted over by exactly
! one period, we set "tt = 0.0".
      tt = time
      IF ( ( ABS ( xmax - vlx0 * tt - xmin) .LT. 1.0e-14) .OR.
     1     ( ABS ( ymax - vly0 * tt - ymin) .LT. 1.0e-14) ) tt = 0.0

      xx = x - vlx0 * tt
      yy = y - vly0 * tt

      IF ( xx .LT. xmin) xx = xx + ( xmax - xmin)
      IF ( yy .LT. ymin) yy = yy + ( ymax - ymin)

      rr2 = xx**2 + yy**2
      exp_rr2 = EXP ( 0.5 * ( 1.0 - rr2) )

! ----------

! Evaluate density (along with its fluctuation) for an isentropic flow
! with a temperature fluctuation specified by "tempaa" below.
      tempaa = - exp_rr2**2 * kappa**2 * ( gamma - 1.0)
     1       / ( 8.0 * gamma * pi**2)
      tempab = tempaa + prs0 / rho0
      tempab = tempab * rho0**gamma / prs0

      prim ( 1) = tempab**( 1.0 / ( gamma - 1.0))

      prim ( 2) = vlx0 - yy * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 3) = vly0 + xx * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 4) = vlz0

      prim ( 5) = prim ( 1) * ( tempaa + prs0 / rho0)

      prim ( 6) = bfx0
      prim ( 7) = bfy0
      prim ( 8) = bfz0

! ----------

      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 3) THEN

! This is the magnetized vortex problem in the yz-plane from Balsara (2004)
! Ap.J.Supp.

! ----------

      rho0 = 1.0
      prs0 = 1.0
      vlx0 = 0.0
      vly0 = 1.0
      vlz0 = 1.0
      bfx0 = 0.0
      bfy0 = 0.0
      bfz0 = 0.0

! Strength of the fluid vorticity.
      kappa = 1.0

! Strength of the curl of the magnetic field.
      mu = SQRT ( 4.0 * pi)

! ----------

! "tt" is the time. However, if the domain has shifted over by exactly
! one period, we set "tt = 0.0".
      tt = time
      IF ( ( ABS ( ymax - vly0 * tt - ymin) .LT. 1.0e-14) .OR.
     1     ( ABS ( zmax - vlz0 * tt - zmin) .LT. 1.0e-14) ) tt = 0.0

      yy = y - vly0 * tt
      zz = z - vlz0 * tt


      IF ( yy .LT. ymin) yy = yy + ( ymax - ymin)
      IF ( zz .LT. zmin) zz = zz + ( zmax - zmin)

      rr2 = yy**2 + zz**2
      exp_rr2 = EXP ( 0.5 * ( 1.0 - rr2) )

! ----------

      prim ( 1) = rho0

      prim ( 2) = vlx0 
      prim ( 3) = vly0 - zz * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 4) = vlz0 + yy * kappa * exp_rr2 / ( 2.0 * pi)

      prim ( 5) = prs0 + ( exp_rr2**2 / ( 4.0 * pi**2) )
     1                 * ( mu**2 * (1.0 - rr2) / (8.0 * pi)
     1                   - 0.5 * kappa**2)

      prim ( 6) = bfx0 
      prim ( 7) = bfy0 - zz * mu * exp_rr2 / ( 2.0 * pi)
      prim ( 8) = bfz0 + yy * mu * exp_rr2 / ( 2.0 * pi)

! ----------

      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 4) THEN


! This is the unmagnetized vortex problem in the yz-plane from
! Balsara & Shu (2000) JCP.

! ----------

      rho0 = 1.0
      prs0 = 1.0
      vlx0 = 0.0
      vly0 = 1.0
      vlz0 = 1.0
      bfx0 = 0.0
      bfy0 = 0.0
      bfz0 = 0.0

! Strength of the fluid vorticity. This "kappa" is the same as "epsilon"
! from Balsara & Shu (2000) JCP.
      kappa = 5.0

! ----------

! "tt" is the time. However, if the domain has shifted over by exactly
! one period, we set "tt = 0.0".
      tt = time
      IF ( ( ABS ( ymax - vly0 * tt - ymin) .LT. 1.0e-14) .OR.
     1     ( ABS ( zmax - vlz0 * tt - zmin) .LT. 1.0e-14) ) tt = 0.0

      yy = y - vly0 * tt
      zz = z - vlz0 * tt

      IF ( yy .LT. ymin) yy = yy + ( ymax - ymin)
      IF ( zz .LT. zmin) zz = zz + ( zmax - zmin)

      rr2 =  yy**2 + zz**2
      exp_rr2 = EXP ( 0.5 * ( 1.0 - rr2) )

! ----------

! Evaluate density (along with its fluctuation) for an isentropic flow
! with a temperature fluctuation specified by "tempaa" below.
      tempaa = - exp_rr2**2 * kappa**2 * ( gamma - 1.0)
     1       / ( 8.0 * gamma * pi**2)
      tempab = tempaa + prs0 / rho0
      tempab = tempab * rho0**gamma / prs0

      prim ( 1) = tempab**( 1.0 / ( gamma - 1.0))

      prim ( 2) = vlx0
      prim ( 3) = vly0 - zz * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 4) = vlz0 + yy * kappa * exp_rr2 / ( 2.0 * pi)

      prim ( 5) = prim ( 1) * ( tempaa + prs0 / rho0)

      prim ( 6) = bfx0
      prim ( 7) = bfy0
      prim ( 8) = bfz0

! ----------

      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 5) THEN

! This is the magnetized vortex problem in the xz-plane from Balsara (2004)
! Ap.J.Supp.

! ----------

      rho0 = 1.0
      prs0 = 1.0
      vlx0 = 1.0
      vly0 = 0.0
      vlz0 = 1.0
      bfx0 = 0.0
      bfy0 = 0.0
      bfz0 = 0.0

! Strength of the fluid vorticity.
      kappa = 1.0

! Strength of the curl of the magnetic field.
      mu = SQRT ( 4.0 * pi)

! ----------

! "tt" is the time. However, if the domain has shifted over by exactly
! one period, we set "tt = 0.0".
      tt = time
      IF ( ( ABS ( xmax - vlx0 * tt - xmin) .LT. 1.0e-14) .OR.
     1     ( ABS ( zmax - vlz0 * tt - zmin) .LT. 1.0e-14) ) tt = 0.0

      xx = x - vlx0 * tt
      zz = z - vlz0 * tt


      IF ( xx .LT. xmin) xx = xx + ( xmax - xmin)
      IF ( zz .LT. zmin) zz = zz + ( zmax - zmin)

      rr2 = xx**2 + zz**2
      exp_rr2 = EXP ( 0.5 * ( 1.0 - rr2) )

! ----------

      prim ( 1) = rho0

      prim ( 2) = vlx0 + zz * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 3) = vly0 
      prim ( 4) = vlz0 - xx * kappa * exp_rr2 / ( 2.0 * pi)

      prim ( 5) = prs0 + ( exp_rr2**2 / ( 4.0 * pi**2) )
     1                 * ( mu**2 * (1.0 - rr2) / (8.0 * pi)
     1                   - 0.5 * kappa**2)

      prim ( 6) = bfx0 + zz * mu * exp_rr2 / ( 2.0 * pi)
      prim ( 7) = bfy0
      prim ( 8) = bfz0 - xx * mu * exp_rr2 / ( 2.0 * pi)

! ----------


      END IF

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 6) THEN

! This is the unmagnetized vortex problem in the xz-plane from
! Balsara & Shu (2000) JCP.

! ----------

      rho0 = 1.0
      prs0 = 1.0
      vlx0 = 1.0
      vly0 = 0.0
      vlz0 = 1.0
      bfx0 = 0.0
      bfy0 = 0.0
      bfz0 = 0.0

! Strength of the fluid vorticity. This "kappa" is the same as "epsilon"
! from Balsara & Shu (2000) JCP.
      kappa = 5.0

! ----------

! "tt" is the time. However, if the domain has shifted over by exactly
! one period, we set "tt = 0.0".
      tt = time
      IF ( ( ABS ( xmax - vlx0 * tt - xmin) .LT. 1.0e-14) .OR.
     1     ( ABS ( zmax - vlz0 * tt - zmin) .LT. 1.0e-14) ) tt = 0.0

      xx = x - vlx0 * tt
      zz = z - vlz0 * tt

      IF ( xx .LT. xmin) xx = xx + ( xmax - xmin)
      IF ( zz .LT. xmin) zz = zz + ( zmax - zmin)

      rr2 =  xx**2 + zz**2
      exp_rr2 = EXP ( 0.5 * ( 1.0 - rr2) )

! ----------

! Evaluate density (along with its fluctuation) for an isentropic flow
! with a temperature fluctuation specified by "tempaa" below.
      tempaa = - exp_rr2**2 * kappa**2 * ( gamma - 1.0)
     1       / ( 8.0 * gamma * pi**2)
      tempab = tempaa + prs0 / rho0
      tempab = tempab * rho0**gamma / prs0

      prim ( 1) = tempab**( 1.0 / ( gamma - 1.0))

      prim ( 2) = vlx0 - zz * kappa * exp_rr2 / ( 2.0 * pi)
      prim ( 3) = vly0
      prim ( 4) = vlz0 + xx * kappa * exp_rr2 / ( 2.0 * pi)

      prim ( 5) = prim ( 1) * ( tempaa + prs0 / rho0)

      prim ( 6) = bfx0
      prim ( 7) = bfy0
      prim ( 8) = bfz0

! ----------


      END IF

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE VORTEX_GENERATOR_MAKE_PRIM

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INIT_ALFVEN3D

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

! Setting up the variables for a higher order scheme can be slightly tricky.
! The present subroutine is provided to facilitate that process.
! It takes the information about a given "time" for a given "iproblem"
! and makes conserved variables at the desired order of accuracy for all the
! flow variables.

! Set "gamma = 1.4" and unit cube [-0.5,0.5]^3 with periodic geometry.
! Stop Alfven problem at "SQRT ( 3.0) / 2.0". Stop entropy wave problem
! at "SQRT ( 3.0)".

! We have to provide a subroutine that generates a vector of primitive
! variables at any space-time point within the computational domain.

! Everything is in place for P1Pm schemes also, i.e. we also
! evaluate the slopes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz, iproblem

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     dtcur,

     1     time, fourpi, max_signal_speed,

     1     tempaa, tempab, tempac, tempad, div_max

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

! Variables that contain the quadrature point info. as well as the
! information about making nodal to modal transcription.

      INTEGER iquadpt, i_euler, ivar_lo, ivar

      REAL xc, yc, zc, dx, dy, dz, xb, yb, zb, xb_m1, yb_m1, zb_m1,
     1     A_x_1, A_x_2, A_x_3, A_x_4, A_x_pl, A_x_mi,
     1     A_y_1, A_y_2, A_y_3, A_y_4, A_y_pl, A_y_mi,
     1     A_z_1, A_z_2, A_z_3, A_z_4, A_z_pl, A_z_mi,
     1     ALFVEN3D_A_x, ALFVEN3D_A_y, ALFVEN3D_A_z

      REAL, DIMENSION ( 1:3, NUM_ZONECEN_QUADPTS) :: zonecen_quadpt,
     1      zc_qp

      REAL, DIMENSION ( 1:3, NUM_FACECEN_QUADPTS) :: xface_quadpt,
     1      yface_quadpt, zface_quadpt, xf_qp, yf_qp, zf_qp

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_QUADPTS) ::
     1      cons_qp, prim_qp, t

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_BASES) ::
     1      cons_modal

! ----------------------------------------------------------------------

! Select a problem number 1 or 2 for MHD/Euler cases.
      iproblem = 1

      IF ( ( ix1 .EQ. 1) .OR. ( iy1 .EQ. 1) .OR. ( iz1 .EQ. 1) ) THEN
      WRITE ( 6, *)"Stopping because INIT_ALFVEN3D must be 3d"
      STOP
      END IF

! Define some temporary variables needed for this problem.

      fourpi = 4.0 * pi

      time = 0.0

      dtcur = 1.0e32
      max_signal_speed = 1.0e-32

! ----------------------------------------------------------------------

! The zone-centered quadrature points in the unit cube are stored
! in "zonecen_quadpt".
! The x-face centered quadrature points in the unit cube are stored
! in "xface_quadpt".
! The y-face centered quadrature points in the unit cube are stored
! in "yface_quadpt".
! The z-face centered quadrature points in the unit cube are stored
! in "zface_quadpt".

! ------------------------------

#if ( NUMORDER == 2)

      zonecen_quadpt ( 1:3, 01) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.0 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      xface_quadpt ( 1:3, 2) = (/ 0.5,-0.5, 0.0 /)
      xface_quadpt ( 1:3, 3) = (/ 0.5, 0.0, 0.5 /)
      xface_quadpt ( 1:3, 4) = (/ 0.5, 0.0,-0.5 /)

      yface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 2) = (/-0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      yface_quadpt ( 1:3, 4) = (/ 0.0, 0.5,-0.5 /)

      zface_quadpt ( 1:3, 1) = (/ 0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 2) = (/-0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      zface_quadpt ( 1:3, 4) = (/ 0.0,-0.5, 0.5 /)

#endif
! End "NUMORDER == 2".

! ------------------------------

#if ( NUMORDER == 3)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/ -0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/ -0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/ -0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 09) = (/ -0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/  0.5,  0.0,  0.0 /)
      xface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      xface_quadpt ( 1:3, 3) = (/  0.5, -0.5,  0.0 /)
      xface_quadpt ( 1:3, 4) = (/  0.5,  0.0,  0.5 /)
      xface_quadpt ( 1:3, 5) = (/  0.5,  0.0, -0.5 /)
      xface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      xface_quadpt ( 1:3, 7) = (/  0.5, -0.5,  0.5 /)
      xface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      xface_quadpt ( 1:3, 9) = (/  0.5, -0.5, -0.5 /)

      yface_quadpt ( 1:3, 1) = (/  0.0,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 3) = (/ -0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 5) = (/  0.0,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 9) = (/ -0.5,  0.5, -0.5 /)

      zface_quadpt ( 1:3, 1) = (/  0.0,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 2) = (/  0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 3) = (/ -0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 5) = (/  0.0, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 8) = (/  0.5, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 9) = (/ -0.5, -0.5,  0.5 /)

#endif
! End "NUMORDER == 3".

! ------------------------------

#if ( NUMORDER == 4)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/ -0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/  0.0 ,  0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/  0.0 , -0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 09) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 ,  0.0 ,  0.25 /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.0 , -0.25 /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 14) = (/  0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 15) = (/  0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 16) = (/  0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 17) = (/  0.5 , -0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 18) = (/ -0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 19) = (/ -0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 20) = (/ -0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 21) = (/ -0.5 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/   0.5,   0.0,   0.0 /)
      xface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      xface_quadpt ( 1:3, 3) = (/   0.5,  0.25,   0.0 /)
      xface_quadpt ( 1:3, 4) = (/   0.5, -0.25,   0.0 /)
      xface_quadpt ( 1:3, 5) = (/   0.5,  -0.5,   0.0 /)
      xface_quadpt ( 1:3, 6) = (/   0.5,   0.0,   0.5 /)
      xface_quadpt ( 1:3, 7) = (/   0.5,   0.0,  0.25 /)
      xface_quadpt ( 1:3, 8) = (/   0.5,   0.0, -0.25 /)
      xface_quadpt ( 1:3, 9) = (/   0.5,   0.0,  -0.5 /)
      xface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      xface_quadpt ( 1:3,11) = (/   0.5,  -0.5,   0.5 /)
      xface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      xface_quadpt ( 1:3,13) = (/   0.5,  -0.5,  -0.5 /)

      yface_quadpt ( 1:3, 1) = (/   0.0,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 3) = (/  0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 4) = (/ -0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 5) = (/  -0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      yface_quadpt ( 1:3, 7) = (/   0.0,   0.5,  0.25 /)
      yface_quadpt ( 1:3, 8) = (/   0.0,   0.5, -0.25 /)
      yface_quadpt ( 1:3, 9) = (/   0.0,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,13) = (/  -0.5,   0.5,  -0.5 /)

      zface_quadpt ( 1:3, 1) = (/   0.0,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 2) = (/   0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 3) = (/  0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 4) = (/ -0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 5) = (/  -0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      zface_quadpt ( 1:3, 7) = (/   0.0,  0.25,   0.5 /)
      zface_quadpt ( 1:3, 8) = (/   0.0, -0.25,   0.5 /)
      zface_quadpt ( 1:3, 9) = (/   0.0,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,12) = (/   0.5,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,13) = (/  -0.5,  -0.5,   0.5 /)

#endif
! End "NUMORDER == 4".

! ----------------------------------------------------------------------

! Set the magnetic field vector to be minus of the curl of the magnetic
! vector potential on the mesh.
! This is the only way for setting up a div-free magnetic field.

! ----------------------------------------------------------------------

      IF ( iproblem .EQ. 1) THEN

! ----------------------------------------------------------------------

! Step a: We make a numerical integraton of the vector potential's
! components along the zone edges with 4th order accuracy. We then
! differentiate this vector potential to get the magnetic field.
! Note though that the differentiation does not drop one order since
! it is being done in a direction that is transverse to the integration.
! This set up should be accurate up to 4th order.

! ------------------------------

! Make "bfxbdy" using a magnetic vector potential approach.

      DO ixx = 0, ix1
        DO iyy = 1, iy1
          DO izz = 1, iz1

! -----

          xb = xb_save ( ixx)

          yc = 0.5 * ( yb_save ( iyy) + yb_save ( iyy - 1) )
          yb = yb_save ( iyy)
          yb_m1 = yb_save ( iyy - 1)
          dy = dy_save ( iyy)

          zc = 0.5 * ( zb_save ( izz) + zb_save ( izz - 1) )
          zb = zb_save ( izz)
          zb_m1 = zb_save ( izz - 1)
          dz = dz_save ( izz)

! -----

          A_y_1 = ALFVEN3D_A_y ( xb, yc - dy/2.0, zb)
          A_y_2 = ALFVEN3D_A_y ( xb, yc - dy/6.0, zb)
          A_y_3 = ALFVEN3D_A_y ( xb, yc + dy/6.0, zb)
          A_y_4 = ALFVEN3D_A_y ( xb, yc + dy/2.0, zb)
          A_y_pl = ( A_y_1 + 3.0 * A_y_2 + 3.0 * A_y_3 + A_y_4 ) / 8.0

          A_y_1 = ALFVEN3D_A_y ( xb, yc - dy/2.0, zb_m1)
          A_y_2 = ALFVEN3D_A_y ( xb, yc - dy/6.0, zb_m1)
          A_y_3 = ALFVEN3D_A_y ( xb, yc + dy/6.0, zb_m1)
          A_y_4 = ALFVEN3D_A_y ( xb, yc + dy/2.0, zb_m1)
          A_y_mi = ( A_y_1 + 3.0 * A_y_2 + 3.0 * A_y_3 + A_y_4 ) / 8.0

! -----

          A_z_1 = ALFVEN3D_A_z ( xb, yb, zc - dz/2.0)
          A_z_2 = ALFVEN3D_A_z ( xb, yb, zc - dz/6.0)
          A_z_3 = ALFVEN3D_A_z ( xb, yb, zc + dz/6.0)
          A_z_4 = ALFVEN3D_A_z ( xb, yb, zc + dz/2.0)
          A_z_pl = ( A_z_1 + 3.0 * A_z_2 + 3.0 * A_z_3 + A_z_4 ) / 8.0

          A_z_1 = ALFVEN3D_A_z ( xb, yb_m1, zc - dz/2.0)
          A_z_2 = ALFVEN3D_A_z ( xb, yb_m1, zc - dz/6.0)
          A_z_3 = ALFVEN3D_A_z ( xb, yb_m1, zc + dz/6.0)
          A_z_4 = ALFVEN3D_A_z ( xb, yb_m1, zc + dz/2.0)
          A_z_mi = ( A_z_1 + 3.0 * A_z_2 + 3.0 * A_z_3 + A_z_4 ) / 8.0

! -----

          bfxbdy ( ixx, iyy, izz, 1) =
     1    ( A_z_pl - A_z_mi ) / ( yb - yb_m1 )
     1  - ( A_y_pl - A_y_mi ) / ( zb - zb_m1 )

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

          xc = 0.5 * ( xb_save ( ixx) + xb_save ( ixx - 1) )
          xb = xb_save ( ixx)
          xb_m1 = xb_save ( ixx - 1)
          dx = dx_save ( ixx)

          yb = yb_save ( iyy)

          zc = 0.5 * ( zb_save ( izz) + zb_save ( izz - 1) )
          zb = zb_save ( izz)
          zb_m1 = zb_save ( izz - 1)
          dz = dz_save ( izz)

! -----

          A_x_1 = ALFVEN3D_A_x ( xc - dx/2.0, yb, zb)
          A_x_2 = ALFVEN3D_A_x ( xc - dx/6.0, yb, zb)
          A_x_3 = ALFVEN3D_A_x ( xc + dx/6.0, yb, zb)
          A_x_4 = ALFVEN3D_A_x ( xc + dx/2.0, yb, zb)
          A_x_pl = ( A_x_1 + 3.0 * A_x_2 + 3.0 * A_x_3 + A_x_4 ) / 8.0

          A_x_1 = ALFVEN3D_A_x ( xc - dx/2.0, yb, zb_m1)
          A_x_2 = ALFVEN3D_A_x ( xc - dx/6.0, yb, zb_m1)
          A_x_3 = ALFVEN3D_A_x ( xc + dx/6.0, yb, zb_m1)
          A_x_4 = ALFVEN3D_A_x ( xc + dx/2.0, yb, zb_m1)
          A_x_mi = ( A_x_1 + 3.0 * A_x_2 + 3.0 * A_x_3 + A_x_4 ) / 8.0

! -----

          A_z_1 = ALFVEN3D_A_z ( xb, yb, zc - dz/2.0)
          A_z_2 = ALFVEN3D_A_z ( xb, yb, zc - dz/6.0)
          A_z_3 = ALFVEN3D_A_z ( xb, yb, zc + dz/6.0)
          A_z_4 = ALFVEN3D_A_z ( xb, yb, zc + dz/2.0)
          A_z_pl = ( A_z_1 + 3.0 * A_z_2 + 3.0 * A_z_3 + A_z_4 ) / 8.0

          A_z_1 = ALFVEN3D_A_z ( xb_m1, yb, zc - dz/2.0)
          A_z_2 = ALFVEN3D_A_z ( xb_m1, yb, zc - dz/6.0)
          A_z_3 = ALFVEN3D_A_z ( xb_m1, yb, zc + dz/6.0)
          A_z_4 = ALFVEN3D_A_z ( xb_m1, yb, zc + dz/2.0)
          A_z_mi = ( A_z_1 + 3.0 * A_z_2 + 3.0 * A_z_3 + A_z_4 ) / 8.0

! -----

          bfybdy ( ixx, iyy, izz, 1) =
     1    ( A_x_pl - A_x_mi ) / ( zb - zb_m1 )
     1  - ( A_z_pl - A_z_mi ) / ( xb - xb_m1 )

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

          xc = 0.5 * ( xb_save ( ixx) + xb_save ( ixx - 1) )
          xb = xb_save ( ixx)
          xb_m1 = xb_save ( ixx - 1)
          dx = dx_save ( ixx)

          yc = 0.5 * ( yb_save ( iyy) + yb_save ( iyy - 1) )
          yb = yb_save ( iyy)
          yb_m1 = yb_save ( iyy - 1)
          dy = dy_save ( iyy)

          zb = zb_save ( izz)

! -----

          A_x_1 = ALFVEN3D_A_x ( xc - dx/2.0, yb, zb)
          A_x_2 = ALFVEN3D_A_x ( xc - dx/6.0, yb, zb)
          A_x_3 = ALFVEN3D_A_x ( xc + dx/6.0, yb, zb)
          A_x_4 = ALFVEN3D_A_x ( xc + dx/2.0, yb, zb)
          A_x_pl = ( A_x_1 + 3.0 * A_x_2 + 3.0 * A_x_3 + A_x_4 ) / 8.0

          A_x_1 = ALFVEN3D_A_x ( xc - dx/2.0, yb_m1, zb)
          A_x_2 = ALFVEN3D_A_x ( xc - dx/6.0, yb_m1, zb)
          A_x_3 = ALFVEN3D_A_x ( xc + dx/6.0, yb_m1, zb)
          A_x_4 = ALFVEN3D_A_x ( xc + dx/2.0, yb_m1, zb)
          A_x_mi = ( A_x_1 + 3.0 * A_x_2 + 3.0 * A_x_3 + A_x_4 ) / 8.0

! -----

          A_y_1 = ALFVEN3D_A_y ( xb, yc - dy/2.0, zb)
          A_y_2 = ALFVEN3D_A_y ( xb, yc - dy/6.0, zb)
          A_y_3 = ALFVEN3D_A_y ( xb, yc + dy/6.0, zb)
          A_y_4 = ALFVEN3D_A_y ( xb, yc + dy/2.0, zb)
          A_y_pl = ( A_y_1 + 3.0 * A_y_2 + 3.0 * A_y_3 + A_y_4 ) / 8.0

          A_y_1 = ALFVEN3D_A_y ( xb_m1, yc - dy/2.0, zb)
          A_y_2 = ALFVEN3D_A_y ( xb_m1, yc - dy/6.0, zb)
          A_y_3 = ALFVEN3D_A_y ( xb_m1, yc + dy/6.0, zb)
          A_y_4 = ALFVEN3D_A_y ( xb_m1, yc + dy/2.0, zb)
          A_y_mi = ( A_y_1 + 3.0 * A_y_2 + 3.0 * A_y_3 + A_y_4 ) / 8.0

! -----

          bfzbdy ( ixx, iyy, izz, 1) =
     1    ( A_y_pl - A_y_mi ) / ( xb - xb_m1 )
     1  - ( A_x_pl - A_x_mi ) / ( yb - yb_m1)

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

      ELSE

! If we get here then we are just setting up a density wave along the diagonal.

      bfxbdy ( :, :, :, 1) = 0.0
      bfybdy ( :, :, :, 1) = 0.0
      bfzbdy ( :, :, :, 1) = 0.0

      END IF
! End of "iproblem .EQ. 1"

! ----------------------------------------------------------------------

! Make all the primitive (and therefrom conserved) variables at the
! zone-centered quadrature points.
! The zone-centered quadrature points in the computational domain are stored
! in "zc_qp".
! The primitive variables at the quadrature points are stored in "prim_qp".
! The conserved variables at the quadrature points are stored in "cons_qp".
! The conserved variables in modal basis are stored in "cons_modal".

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------------------------------------------------------------------

          xc = xc_save ( ixx)
          dx = dx_save ( ixx)
          yc = yc_save ( iyy)
          dy = dy_save ( iyy)
          zc = zc_save ( izz)
          dz = dz_save ( izz)

          zc_qp ( 1, :) = zonecen_quadpt ( 1, :) * dx + xc
          zc_qp ( 2, :) = zonecen_quadpt ( 2, :) * dy + yc
          zc_qp ( 3, :) = zonecen_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! --------------------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, zc_qp ( 1, iquadpt), zc_qp ( 2, iquadpt),
     1    zc_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! --------------------

! Only the MHD fluid is used for constraining the timestep. If the
! Euler fluids can also constrain the timestep then include them here.

! This is the flow speed.
          tempaa = SQRT ( prim_qp ( 2, iquadpt)**2
     1                  + prim_qp ( 3, iquadpt)**2
     1                  + prim_qp ( 4, iquadpt)**2)
          tempab = gamma * prim_qp ( 5, iquadpt)
     1           + (1.0/fourpi) * ( prim_qp ( 6, iquadpt)**2
     1                            + prim_qp ( 7, iquadpt)**2
     1                            + prim_qp ( 8, iquadpt)**2)
! This is the magnetosonic signal speed.
          tempab = SQRT ( tempab / prim_qp ( 1, iquadpt))
! This is the flow + magnetosonic signal speed.
          tempac = tempaa + tempab

          tempad = cfl_coef * AMIN1 ( dx, dy, dz) / tempac
          dtcur = AMIN1 ( dtcur, tempad)

          max_signal_speed = AMAX1 ( max_signal_speed, tempac)

          cons_qp ( 1, iquadpt) = prim_qp ( 1, iquadpt)

          cons_qp ( 2, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 2, iquadpt)
          cons_qp ( 3, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 3, iquadpt)
          cons_qp ( 4, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 4, iquadpt)

          cons_qp ( 5, iquadpt) = prim_qp ( 5, iquadpt) / (gamma - 1.0)
     1     + 0.5 * prim_qp ( 1, iquadpt) * ( prim_qp ( 2, iquadpt)**2
     1           + prim_qp ( 3, iquadpt)**2 + prim_qp ( 4, iquadpt)**2)
     1     + ( prim_qp ( 6, iquadpt)**2 + prim_qp ( 7, iquadpt)**2
     1       + prim_qp ( 8, iquadpt)**2) / ( 8.0 * pi)

          cons_qp ( 6, iquadpt) = prim_qp ( 6, iquadpt)
          cons_qp ( 7, iquadpt) = prim_qp ( 7, iquadpt)
          cons_qp ( 8, iquadpt) = prim_qp ( 8, iquadpt)

#if ( NFLUID >= 1)
          cons_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1             prim_qp ( 1, iquadpt)
     1           * prim_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
#endif

          END DO
! End of "DO iquadpt ="

! ------------------------------

! Now that the conserved variables are available at all the quadrature
! points, we make the zone-centered conserved variables.

! For good measure, we make all their moments too,
! though they may not be needed.

! ----------

#if ( NUMORDER == 2)

          cons_modal ( :, 1) = ( cons_qp ( :, 01) + cons_qp ( :, 02)
     1                         + cons_qp ( :, 03) + cons_qp ( :, 04)
     1                         + cons_qp ( :, 05) + cons_qp ( :, 06))
     1                         / 6.0

          cons_modal ( :, 2) = cons_qp ( :, 01) - cons_qp ( :, 02)
          cons_modal ( :, 3) = cons_qp ( :, 03) - cons_qp ( :, 04)
          cons_modal ( :, 4) = cons_qp ( :, 05) - cons_qp ( :, 06)

#endif
! End of "NUMORDER == 2"
! ----------

#if ( NUMORDER == 3)

      t ( :, 01) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 04) + cons_qp ( :, 05)

      t ( :, 02) = cons_qp ( :, 06) + cons_qp ( :, 07)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 03) = cons_qp ( :, 10) + cons_qp ( :, 11)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 04) = cons_qp ( :, 02) + cons_qp ( :, 06)
     1           + cons_qp ( :, 10) + cons_qp ( :, 12)

      t ( :, 05) = cons_qp ( :, 03) + cons_qp ( :, 07)
     1           + cons_qp ( :, 11) + cons_qp ( :, 13)

      t ( :, 06) = cons_qp ( :, 04) + cons_qp ( :, 05)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 07) = cons_qp ( :, 04) + cons_qp ( :, 08)
     1           + cons_qp ( :, 10) + cons_qp ( :, 11)

      t ( :, 08) = cons_qp ( :, 05) + cons_qp ( :, 09)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 09) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 02) = 0.25 * ( t ( :, 01) - t ( :, 02) )

      cons_modal ( :, 03) = 0.25 * ( t ( :, 04) - t ( :, 05) )

      cons_modal ( :, 04) = 0.25 * ( t ( :, 07) - t ( :, 08) )

      cons_modal ( :, 05) = 0.5 * ( t ( :, 01) + t ( :, 02)
     1 - t ( :, 03) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 06) = 0.5 * ( t ( :, 04) + t ( :, 05)
     1 - t ( :, 06) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 07) = 0.5 * ( t ( :, 07) + t ( :, 08)
     1 - t ( :, 09) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 08) = cons_qp ( :, 02) - cons_qp ( :, 03)
     1 - cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 09) = cons_qp ( :, 10) - cons_qp ( :, 11)
     1 - cons_qp ( :, 12) + cons_qp ( :, 13)

      cons_modal ( :, 10) = cons_qp ( :, 04) - cons_qp ( :, 05)
     1 - cons_qp ( :, 08) + cons_qp ( :, 09)

      cons_modal ( :, 01) = cons_qp ( :, 01) + ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 12.0

#endif
! End of "NUMORDER == 3".
! ----------

#if ( NUMORDER == 4)

! Note that in this formulation "cons_modal ( :, 14:19)" require
! space-time nodal points. They, therefore, can't be built here.

      cons_modal ( :, 05) = 2.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 05) )

      cons_modal ( :, 06) = 2.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 09) )

      cons_modal ( :, 07) = 2.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 13) )

      cons_modal ( :, 08) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 15)
     1 + cons_qp ( :, 16) - cons_qp ( :, 17) - cons_qp ( :, 18)
     1 + cons_qp ( :, 19) - cons_qp ( :, 20) + cons_qp ( :, 21) )

      cons_modal ( :, 09) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 16)
     1 + cons_qp ( :, 18) - cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 17) - cons_qp ( :, 19) + cons_qp ( :, 21) )

      cons_modal ( :, 10) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 18)
     1 + cons_qp ( :, 15) - cons_qp ( :, 19) - cons_qp ( :, 16)
     1 + cons_qp ( :, 20) - cons_qp ( :, 17) + cons_qp ( :, 21) )

      cons_modal ( :, 11) = 16.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 03) + 2.0 * cons_qp ( :, 04)
     1 - cons_qp ( :, 05) ) / 3.0

      cons_modal ( :, 12) = 16.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 07) + 2.0 * cons_qp ( :, 08)
     1 - cons_qp ( :, 09) ) / 3.0

      cons_modal ( :, 13) = 16.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 11) + 2.0 * cons_qp ( :, 12)
     1 - cons_qp ( :, 13) ) / 3.0

      t ( :, 01) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 16) + cons_qp ( :, 17)

      t ( :, 02) = cons_qp ( :, 18) + cons_qp ( :, 19)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      t ( :, 03) = cons_qp ( :, 14) + cons_qp ( :, 16)
     1 + cons_qp ( :, 18) + cons_qp ( :, 20)

      t ( :, 04) = cons_qp ( :, 15) + cons_qp ( :, 17)
     1 + cons_qp ( :, 19) + cons_qp ( :, 21)

      t ( :, 05) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19)

      t ( :, 06) = cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      cons_modal ( :, 02) = 2.0 * ( cons_qp ( :, 02)
     1 - cons_qp ( :, 05) ) / 3.0 - 0.1 * cons_modal ( :, 11)
     1 + 0.25 * ( t ( :, 01) - t ( :, 02) ) / 3.0

      cons_modal ( :, 03) = 2.0 * ( cons_qp ( :, 06)
     1 - cons_qp ( :, 09) ) / 3.0 - 0.1 * cons_modal ( :, 12)
     1 + 0.25 * ( t ( :, 03) - t ( :, 04) ) / 3.0

      cons_modal ( :, 04) = 2.0 * ( cons_qp ( :, 10)
     1 - cons_qp ( :, 13) ) / 3.0 - 0.1 * cons_modal ( :, 13)
     1 + 0.25 * ( t ( :, 05) - t ( :, 06) ) / 3.0

      cons_modal ( :, 20) = cons_qp ( :, 14) - cons_qp ( :, 18)
     1 - cons_qp ( :, 16) + cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 19) + cons_qp ( :, 17) - cons_qp ( :, 21)

      cons_modal ( :, 01) = 0.125 * ( cons_qp ( :, 14)
     1 + cons_qp ( :, 15) + cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19) + cons_qp ( :, 20)
     1 + cons_qp ( :, 21) ) - 0.5 * ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 3.0

#endif
! End of "NUMORDER == 4".
! ----------
! ------------------------------

! Now store "cons_modal ( :, 1)" in the code's conserved variable meshes.

          u ( ixx, iyy, izz, :, 1) = cons_modal ( :, 1)

! ----------------------------------------------------------------------

! Make all the x-magnetic field variables at the top x-face-centered
! quadrature points.
! The top x-face-centered quadrature points in the computational domain
! are stored in "xf_qp".
! The x-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 6, :)".

! --------------------

          xf_qp ( 1, :) = xface_quadpt ( 1, :) * dx + xc
          xf_qp ( 2, :) = xface_quadpt ( 2, :) * dy + yc
          xf_qp ( 3, :) = xface_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, xf_qp ( 1, iquadpt), xf_qp ( 2, iquadpt),
     1    xf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 6, 01) = 0.25 * ( prim_qp ( 6, 01)
     1       + prim_qp ( 6, 02) + prim_qp ( 6, 03) + prim_qp ( 6, 04))

          cons_modal ( 6, 02) = prim_qp ( 6, 01) - prim_qp ( 6, 02)

          cons_modal ( 6, 03) = prim_qp ( 6, 03) - prim_qp ( 6, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 03)

          cons_modal ( 6, 03) = prim_qp ( 6, 04) - prim_qp ( 6, 05)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 03)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 04) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 06) = prim_qp ( 6, 06) - prim_qp ( 6, 07)
     1       - prim_qp ( 6, 08) + prim_qp ( 6, 09)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 06) + 2.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 06) = prim_qp ( 6, 10) - prim_qp ( 6, 11)
     1       - prim_qp ( 6, 12) + prim_qp ( 6, 13)

          cons_modal ( 6, 07) = ( 16.0 * prim_qp ( 6, 02)
     1       - 16.0 * prim_qp ( 6, 05) - 32.0 * prim_qp ( 6, 03)
     1       + 32.0 * prim_qp ( 6, 04) ) / 3.0

          cons_modal ( 6, 08) = ( 16.0 * prim_qp ( 6, 06)
     1       - 16.0 * prim_qp ( 6, 09) - 32.0 * prim_qp ( 6, 07)
     1       + 32.0 * prim_qp ( 6, 08) ) / 3.0

          cons_modal ( 6, 09) = 2.0 * prim_qp ( 6, 10)
     1       + 2.0 * prim_qp ( 6, 11) - 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 06)
     1       + 4.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 10) = 2.0 * prim_qp ( 6, 10)
     1       - 2.0 * prim_qp ( 6, 11) + 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 02)
     1       + 4.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 05)
     1       - 0.1 * cons_modal ( 6, 07) + cons_modal ( 6, 10) / 12.0

          cons_modal ( 6, 03) = prim_qp ( 6, 06) - prim_qp ( 6, 09)
     1       - 0.1 * cons_modal ( 6, 08) + cons_modal ( 6, 09) / 12.0

#endif
! ----------

!          bfxbdy ( ixx, iyy, izz, 1) = cons_modal ( 6, 1)

! ----------------------------------------------------------------------

! Make all the x-magnetic field variables at the bottom x-face-centered
! quadrature points.
! The bottom x-face-centered quadrature points in the computational domain
! are stored in "xf_qp".
! The x-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 6, :)".

! --------------------

          xf_qp ( 1, :) = - xface_quadpt ( 1, :) * dx + xc
          xf_qp ( 2, :) = xface_quadpt ( 2, :) * dy + yc
          xf_qp ( 3, :) = xface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, xf_qp ( 1, iquadpt), xf_qp ( 2, iquadpt),
     1    xf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 6, 01) = 0.25 * ( prim_qp ( 6, 01)
     1       + prim_qp ( 6, 02) + prim_qp ( 6, 03) + prim_qp ( 6, 04))

          cons_modal ( 6, 02) = prim_qp ( 6, 01) - prim_qp ( 6, 02)

          cons_modal ( 6, 03) = prim_qp ( 6, 03) - prim_qp ( 6, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 03)

          cons_modal ( 6, 03) = prim_qp ( 6, 04) - prim_qp ( 6, 05)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 03)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 04) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 06) = prim_qp ( 6, 06) - prim_qp ( 6, 07)
     1       - prim_qp ( 6, 08) + prim_qp ( 6, 09)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 06) + 2.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 06) = prim_qp ( 6, 10) - prim_qp ( 6, 11)
     1       - prim_qp ( 6, 12) + prim_qp ( 6, 13)

          cons_modal ( 6, 07) = ( 16.0 * prim_qp ( 6, 02)
     1       - 16.0 * prim_qp ( 6, 05) - 32.0 * prim_qp ( 6, 03)
     1       + 32.0 * prim_qp ( 6, 04) ) / 3.0

          cons_modal ( 6, 08) = ( 16.0 * prim_qp ( 6, 06)
     1       - 16.0 * prim_qp ( 6, 09) - 32.0 * prim_qp ( 6, 07)
     1       + 32.0 * prim_qp ( 6, 08) ) / 3.0

          cons_modal ( 6, 09) = 2.0 * prim_qp ( 6, 10)
     1       + 2.0 * prim_qp ( 6, 11) - 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 06)
     1       + 4.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 10) = 2.0 * prim_qp ( 6, 10)
     1       - 2.0 * prim_qp ( 6, 11) + 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 02)
     1       + 4.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 05)
     1       - 0.1 * cons_modal ( 6, 07) + cons_modal ( 6, 10) / 12.0

          cons_modal ( 6, 03) = prim_qp ( 6, 06) - prim_qp ( 6, 09)
     1       - 0.1 * cons_modal ( 6, 08) + cons_modal ( 6, 09) / 12.0

#endif
! ----------

!          bfxbdy ( ixx - 1, iyy, izz, 1) = cons_modal ( 6, 1)

! ----------------------------------------------------------------------

! Make all the y-magnetic field variables at the top y-face-centered
! quadrature points.
! The top y-face-centered quadrature points in the computational domain
! are stored in "yf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 7, :)".

! --------------------

          yf_qp ( 1, :) = yface_quadpt ( 1, :) * dx + xc
          yf_qp ( 2, :) = yface_quadpt ( 2, :) * dy + yc
          yf_qp ( 3, :) = yface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, yf_qp ( 1, iquadpt), yf_qp ( 2, iquadpt),
     1    yf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 7, 01) = 0.25 * ( prim_qp ( 7, 01)
     1       + prim_qp ( 7, 02) + prim_qp ( 7, 03) + prim_qp ( 7, 04))

          cons_modal ( 7, 02) = prim_qp ( 7, 01) - prim_qp ( 7, 02)

          cons_modal ( 7, 03) = prim_qp ( 7, 03) - prim_qp ( 7, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 03)

          cons_modal ( 7, 03) = prim_qp ( 7, 04) - prim_qp ( 7, 05)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 03)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 04) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 06) = prim_qp ( 7, 06) - prim_qp ( 7, 07)
     1       - prim_qp ( 7, 08) + prim_qp ( 7, 09)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 06) + 2.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 06) = prim_qp ( 7, 10) - prim_qp ( 7, 11)
     1       - prim_qp ( 7, 12) + prim_qp ( 7, 13)

          cons_modal ( 7, 07) = ( 16.0 * prim_qp ( 7, 02)
     1       - 16.0 * prim_qp ( 7, 05) - 32.0 * prim_qp ( 7, 03)
     1       + 32.0 * prim_qp ( 7, 04) ) / 3.0

          cons_modal ( 7, 08) = ( 16.0 * prim_qp ( 7, 06)
     1       - 16.0 * prim_qp ( 7, 09) - 32.0 * prim_qp ( 7, 07)
     1       + 32.0 * prim_qp ( 7, 08) ) / 3.0

          cons_modal ( 7, 09) = 2.0 * prim_qp ( 7, 10)
     1       + 2.0 * prim_qp ( 7, 11) - 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 06)
     1       + 4.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 10) = 2.0 * prim_qp ( 7, 10)
     1       - 2.0 * prim_qp ( 7, 11) + 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 02)
     1       + 4.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 05)
     1       - 0.1 * cons_modal ( 7, 07) + cons_modal ( 7, 10) / 12.0

          cons_modal ( 7, 03) = prim_qp ( 7, 06) - prim_qp ( 7, 09)
     1       - 0.1 * cons_modal ( 7, 08) + cons_modal ( 7, 09) / 12.0

#endif
! ----------

!          bfybdy ( ixx, iyy, izz, 1) = cons_modal ( 7, 1)

! ----------------------------------------------------------------------

! Make all the y-magnetic field variables at the bottom y-face-centered
! quadrature points.
! The bottom y-face-centered quadrature points in the computational domain
! are stored in "yf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 7, :)".

! --------------------

          yf_qp ( 1, :) = yface_quadpt ( 1, :) * dx + xc
          yf_qp ( 2, :) = - yface_quadpt ( 2, :) * dy + yc
          yf_qp ( 3, :) = yface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, yf_qp ( 1, iquadpt), yf_qp ( 2, iquadpt),
     1    yf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 7, 01) = 0.25 * ( prim_qp ( 7, 01)
     1       + prim_qp ( 7, 02) + prim_qp ( 7, 03) + prim_qp ( 7, 04))

          cons_modal ( 7, 02) = prim_qp ( 7, 01) - prim_qp ( 7, 02)

          cons_modal ( 7, 03) = prim_qp ( 7, 03) - prim_qp ( 7, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 03)

          cons_modal ( 7, 03) = prim_qp ( 7, 04) - prim_qp ( 7, 05)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 03)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 04) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 06) = prim_qp ( 7, 06) - prim_qp ( 7, 07)
     1       - prim_qp ( 7, 08) + prim_qp ( 7, 09)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 06) + 2.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 06) = prim_qp ( 7, 10) - prim_qp ( 7, 11)
     1       - prim_qp ( 7, 12) + prim_qp ( 7, 13)

          cons_modal ( 7, 07) = ( 16.0 * prim_qp ( 7, 02)
     1       - 16.0 * prim_qp ( 7, 05) - 32.0 * prim_qp ( 7, 03)
     1       + 32.0 * prim_qp ( 7, 04) ) / 3.0

          cons_modal ( 7, 08) = ( 16.0 * prim_qp ( 7, 06)
     1       - 16.0 * prim_qp ( 7, 09) - 32.0 * prim_qp ( 7, 07)
     1       + 32.0 * prim_qp ( 7, 08) ) / 3.0

          cons_modal ( 7, 09) = 2.0 * prim_qp ( 7, 10)
     1       + 2.0 * prim_qp ( 7, 11) - 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 06)
     1       + 4.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 10) = 2.0 * prim_qp ( 7, 10)
     1       - 2.0 * prim_qp ( 7, 11) + 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 02)
     1       + 4.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 05)
     1       - 0.1 * cons_modal ( 7, 07) + cons_modal ( 7, 10) / 12.0

          cons_modal ( 7, 03) = prim_qp ( 7, 06) - prim_qp ( 7, 09)
     1       - 0.1 * cons_modal ( 7, 08) + cons_modal ( 7, 09) / 12.0

#endif
! ----------

!          bfybdy ( ixx, iyy - 1, izz, 1) = cons_modal ( 7, 1)

! ----------------------------------------------------------------------

! Make all the z-magnetic field variables at the top z-face-centered
! quadrature points.
! The top z-face-centered quadrature points in the computational domain
! are stored in "zf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 8, :)".

! --------------------

          zf_qp ( 1, :) = zface_quadpt ( 1, :) * dx + xc
          zf_qp ( 2, :) = zface_quadpt ( 2, :) * dy + yc
          zf_qp ( 3, :) = zface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, zf_qp ( 1, iquadpt), zf_qp ( 2, iquadpt),
     1    zf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 8, 01) = 0.25 * ( prim_qp ( 8, 01)
     1       + prim_qp ( 8, 02) + prim_qp ( 8, 03) + prim_qp ( 8, 04))

          cons_modal ( 8, 02) = prim_qp ( 8, 01) - prim_qp ( 8, 02)

          cons_modal ( 8, 03) = prim_qp ( 8, 03) - prim_qp ( 8, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 03)

          cons_modal ( 8, 03) = prim_qp ( 8, 04) - prim_qp ( 8, 05)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 03)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 04) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 06) = prim_qp ( 8, 06) - prim_qp ( 8, 07)
     1       - prim_qp ( 8, 08) + prim_qp ( 8, 09)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 06) + 2.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 06) = prim_qp ( 8, 10) - prim_qp ( 8, 11)
     1       - prim_qp ( 8, 12) + prim_qp ( 8, 13)

          cons_modal ( 8, 07) = ( 16.0 * prim_qp ( 8, 02)
     1       - 16.0 * prim_qp ( 8, 05) - 32.0 * prim_qp ( 8, 03)
     1       + 32.0 * prim_qp ( 8, 04) ) / 3.0

          cons_modal ( 8, 08) = ( 16.0 * prim_qp ( 8, 06)
     1       - 16.0 * prim_qp ( 8, 09) - 32.0 * prim_qp ( 8, 07)
     1       + 32.0 * prim_qp ( 8, 08) ) / 3.0

          cons_modal ( 8, 09) = 2.0 * prim_qp ( 8, 10)
     1       + 2.0 * prim_qp ( 8, 11) - 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 06)
     1       + 4.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 10) = 2.0 * prim_qp ( 8, 10)
     1       - 2.0 * prim_qp ( 8, 11) + 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 02)
     1       + 4.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 05)
     1       - 0.1 * cons_modal ( 8, 07) + cons_modal ( 8, 10) / 12.0

          cons_modal ( 8, 03) = prim_qp ( 8, 06) - prim_qp ( 8, 09)
     1       - 0.1 * cons_modal ( 8, 08) + cons_modal ( 8, 09) / 12.0

#endif
! ----------

!          bfzbdy ( ixx, iyy, izz, 1) = cons_modal ( 8, 1)

! ----------------------------------------------------------------------

! Make all the z-magnetic field variables at the bottom z-face-centered
! quadrature points.
! The bottom z-face-centered quadrature points in the computational domain
! are stored in "zf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 8, :)".

! --------------------

          zf_qp ( 1, :) = zface_quadpt ( 1, :) * dx + xc
          zf_qp ( 2, :) = zface_quadpt ( 2, :) * dy + yc
          zf_qp ( 3, :) = - zface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, zf_qp ( 1, iquadpt), zf_qp ( 2, iquadpt),
     1    zf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 8, 01) = 0.25 * ( prim_qp ( 8, 01)
     1       + prim_qp ( 8, 02) + prim_qp ( 8, 03) + prim_qp ( 8, 04))

          cons_modal ( 8, 02) = prim_qp ( 8, 01) - prim_qp ( 8, 02)

          cons_modal ( 8, 03) = prim_qp ( 8, 03) - prim_qp ( 8, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 03)

          cons_modal ( 8, 03) = prim_qp ( 8, 04) - prim_qp ( 8, 05)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 03)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 04) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 06) = prim_qp ( 8, 06) - prim_qp ( 8, 07)
     1       - prim_qp ( 8, 08) + prim_qp ( 8, 09)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 06) + 2.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 06) = prim_qp ( 8, 10) - prim_qp ( 8, 11)
     1       - prim_qp ( 8, 12) + prim_qp ( 8, 13)

          cons_modal ( 8, 07) = ( 16.0 * prim_qp ( 8, 02)
     1       - 16.0 * prim_qp ( 8, 05) - 32.0 * prim_qp ( 8, 03)
     1       + 32.0 * prim_qp ( 8, 04) ) / 3.0

          cons_modal ( 8, 08) = ( 16.0 * prim_qp ( 8, 06)
     1       - 16.0 * prim_qp ( 8, 09) - 32.0 * prim_qp ( 8, 07)
     1       + 32.0 * prim_qp ( 8, 08) ) / 3.0

          cons_modal ( 8, 09) = 2.0 * prim_qp ( 8, 10)
     1       + 2.0 * prim_qp ( 8, 11) - 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 06)
     1       + 4.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 10) = 2.0 * prim_qp ( 8, 10)
     1       - 2.0 * prim_qp ( 8, 11) + 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 02)
     1       + 4.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 05)
     1       - 0.1 * cons_modal ( 8, 07) + cons_modal ( 8, 10) / 12.0

          cons_modal ( 8, 03) = prim_qp ( 8, 06) - prim_qp ( 8, 09)
     1       - 0.1 * cons_modal ( 8, 08) + cons_modal ( 8, 09) / 12.0

#endif
! ----------

!          bfzbdy ( ixx, iyy, izz - 1, 1) = cons_modal ( 8, 1)

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INIT_ALFVEN3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ERROR_ALFVEN3D

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

! Takes the data at a given "time" for a given "iproblem" and makes
! and writes out L_1 and L_infinity error norms for all the flow variables
! for it. This is only done for serial setting. For parallel setting
! we would have to accumulate "L_1_error and L_inf_error" from all the
! processors.

! We have to provide a subroutine that generates a vector of primitive
! variables at any space-time point within the computational domain.

! Everything is in place for P1Pm schemes also, i.e. we also
! evaluate the slopes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        ixx, iyy, izz, iproblem

      REAL xmin, xmax, ymin, ymax, zmin, zmax,

     1     time,

     1     tempaa, tempab, tempac, fourpi

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

! Variables that contain the quadrature point info. as well as the
! information about making nodal to modal transcription.

      INTEGER iquadpt, i_euler, ivar_lo, ivar

      REAL xc, yc, zc, dx, dy, dz

      REAL, DIMENSION ( 1:3, NUM_ZONECEN_QUADPTS) :: zonecen_quadpt,
     1      zc_qp

      REAL, DIMENSION ( 1:3, NUM_FACECEN_QUADPTS) :: xface_quadpt,
     1      yface_quadpt, zface_quadpt, xf_qp, yf_qp, zf_qp

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_QUADPTS) ::
     1      cons_qp, prim_qp, t

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NUM_ZONECEN_BASES) ::
     1      cons_modal

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: diff_cons,
     1      L_1_error, L_inf_error

! ----------------------------------------------------------------------

! Select a problem number 1 or 2 for MHD/Euler cases.
      iproblem = 1

! Define some temporary variables needed for this problem.

      fourpi = 4.0 * pi

      L_1_error ( :) = 0.0
      L_inf_error ( :) = 0.0

! ----------------------------------------------------------------------

! The zone-centered quadrature points in the unit cube are stored
! in "zonecen_quadpt".
! The x-face centered quadrature points in the unit cube are stored
! in "xface_quadpt".
! The y-face centered quadrature points in the unit cube are stored
! in "yface_quadpt".
! The z-face centered quadrature points in the unit cube are stored
! in "zface_quadpt".

! ------------------------------

#if ( NUMORDER == 2)

      zonecen_quadpt ( 1:3, 01) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.0 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      xface_quadpt ( 1:3, 2) = (/ 0.5,-0.5, 0.0 /)
      xface_quadpt ( 1:3, 3) = (/ 0.5, 0.0, 0.5 /)
      xface_quadpt ( 1:3, 4) = (/ 0.5, 0.0,-0.5 /)

      yface_quadpt ( 1:3, 1) = (/ 0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 2) = (/-0.5, 0.5, 0.0 /)
      yface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      yface_quadpt ( 1:3, 4) = (/ 0.0, 0.5,-0.5 /)

      zface_quadpt ( 1:3, 1) = (/ 0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 2) = (/-0.5, 0.0, 0.5 /)
      zface_quadpt ( 1:3, 3) = (/ 0.0, 0.5, 0.5 /)
      zface_quadpt ( 1:3, 4) = (/ 0.0,-0.5, 0.5 /)

#endif
! End "NUMORDER == 2".

! ------------------------------

#if ( NUMORDER == 3)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/  0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 05) = (/  0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 06) = (/ -0.5 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/ -0.5 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/ -0.5 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 09) = (/ -0.5 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/  0.5,  0.0,  0.0 /)
      xface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      xface_quadpt ( 1:3, 3) = (/  0.5, -0.5,  0.0 /)
      xface_quadpt ( 1:3, 4) = (/  0.5,  0.0,  0.5 /)
      xface_quadpt ( 1:3, 5) = (/  0.5,  0.0, -0.5 /)
      xface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      xface_quadpt ( 1:3, 7) = (/  0.5, -0.5,  0.5 /)
      xface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      xface_quadpt ( 1:3, 9) = (/  0.5, -0.5, -0.5 /)

      yface_quadpt ( 1:3, 1) = (/  0.0,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 2) = (/  0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 3) = (/ -0.5,  0.5,  0.0 /)
      yface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 5) = (/  0.0,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      yface_quadpt ( 1:3, 8) = (/  0.5,  0.5, -0.5 /)
      yface_quadpt ( 1:3, 9) = (/ -0.5,  0.5, -0.5 /)

      zface_quadpt ( 1:3, 1) = (/  0.0,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 2) = (/  0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 3) = (/ -0.5,  0.0,  0.5 /)
      zface_quadpt ( 1:3, 4) = (/  0.0,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 5) = (/  0.0, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 6) = (/  0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 7) = (/ -0.5,  0.5,  0.5 /)
      zface_quadpt ( 1:3, 8) = (/  0.5, -0.5,  0.5 /)
      zface_quadpt ( 1:3, 9) = (/ -0.5, -0.5,  0.5 /)

#endif
! End "NUMORDER == 3".

! ------------------------------

#if ( NUMORDER == 4)

      zonecen_quadpt ( 1:3, 01) = (/  0.0 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 02) = (/  0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 03) = (/  0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 04) = (/ -0.25,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 05) = (/ -0.5 ,  0.0 ,  0.0  /)
      zonecen_quadpt ( 1:3, 06) = (/  0.0 ,  0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 07) = (/  0.0 ,  0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 08) = (/  0.0 , -0.25,  0.0  /)
      zonecen_quadpt ( 1:3, 09) = (/  0.0 , -0.5 ,  0.0  /)
      zonecen_quadpt ( 1:3, 10) = (/  0.0 ,  0.0 ,  0.5  /)
      zonecen_quadpt ( 1:3, 11) = (/  0.0 ,  0.0 ,  0.25 /)
      zonecen_quadpt ( 1:3, 12) = (/  0.0 ,  0.0 , -0.25 /)
      zonecen_quadpt ( 1:3, 13) = (/  0.0 ,  0.0 , -0.5  /)
      zonecen_quadpt ( 1:3, 14) = (/  0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 15) = (/  0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 16) = (/  0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 17) = (/  0.5 , -0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 18) = (/ -0.5 ,  0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 19) = (/ -0.5 , -0.5 ,  0.5  /)
      zonecen_quadpt ( 1:3, 20) = (/ -0.5 ,  0.5 , -0.5  /)
      zonecen_quadpt ( 1:3, 21) = (/ -0.5 , -0.5 , -0.5  /)

      xface_quadpt ( 1:3, 1) = (/   0.5,   0.0,   0.0 /)
      xface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      xface_quadpt ( 1:3, 3) = (/   0.5,  0.25,   0.0 /)
      xface_quadpt ( 1:3, 4) = (/   0.5, -0.25,   0.0 /)
      xface_quadpt ( 1:3, 5) = (/   0.5,  -0.5,   0.0 /)
      xface_quadpt ( 1:3, 6) = (/   0.5,   0.0,   0.5 /)
      xface_quadpt ( 1:3, 7) = (/   0.5,   0.0,  0.25 /)
      xface_quadpt ( 1:3, 8) = (/   0.5,   0.0, -0.25 /)
      xface_quadpt ( 1:3, 9) = (/   0.5,   0.0,  -0.5 /)
      xface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      xface_quadpt ( 1:3,11) = (/   0.5,  -0.5,   0.5 /)
      xface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      xface_quadpt ( 1:3,13) = (/   0.5,  -0.5,  -0.5 /)

      yface_quadpt ( 1:3, 1) = (/   0.0,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 2) = (/   0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 3) = (/  0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 4) = (/ -0.25,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 5) = (/  -0.5,   0.5,   0.0 /)
      yface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      yface_quadpt ( 1:3, 7) = (/   0.0,   0.5,  0.25 /)
      yface_quadpt ( 1:3, 8) = (/   0.0,   0.5, -0.25 /)
      yface_quadpt ( 1:3, 9) = (/   0.0,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      yface_quadpt ( 1:3,12) = (/   0.5,   0.5,  -0.5 /)
      yface_quadpt ( 1:3,13) = (/  -0.5,   0.5,  -0.5 /)

      zface_quadpt ( 1:3, 1) = (/   0.0,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 2) = (/   0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 3) = (/  0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 4) = (/ -0.25,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 5) = (/  -0.5,   0.0,   0.5 /)
      zface_quadpt ( 1:3, 6) = (/   0.0,   0.5,   0.5 /)
      zface_quadpt ( 1:3, 7) = (/   0.0,  0.25,   0.5 /)
      zface_quadpt ( 1:3, 8) = (/   0.0, -0.25,   0.5 /)
      zface_quadpt ( 1:3, 9) = (/   0.0,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,10) = (/   0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,11) = (/  -0.5,   0.5,   0.5 /)
      zface_quadpt ( 1:3,12) = (/   0.5,  -0.5,   0.5 /)
      zface_quadpt ( 1:3,13) = (/  -0.5,  -0.5,   0.5 /)

#endif
! End "NUMORDER == 4".

! ----------------------------------------------------------------------

! Make all the primitive (and therefrom conserved) variables at the
! zone-centered quadrature points.
! The zone-centered quadrature points in the computational domain are stored
! in "zc_qp".
! The primitive variables at the quadrature points are stored in "prim_qp".
! The conserved variables at the quadrature points are stored in "cons_qp".
! The conserved variables in modal basis are stored in "cons_modal".

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------------------------------------------------------------------

          xc = xc_save ( ixx)
          dx = dx_save ( ixx)
          yc = yc_save ( iyy)
          dy = dy_save ( iyy)
          zc = zc_save ( izz)
          dz = dz_save ( izz)

          zc_qp ( 1, :) = zonecen_quadpt ( 1, :) * dx + xc
          zc_qp ( 2, :) = zonecen_quadpt ( 2, :) * dy + yc
          zc_qp ( 3, :) = zonecen_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_ZONECEN_QUADPTS

! --------------------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, zc_qp ( 1, iquadpt), zc_qp ( 2, iquadpt),
     1    zc_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! --------------------

          cons_qp ( 1, iquadpt) = prim_qp ( 1, iquadpt)

          cons_qp ( 2, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 2, iquadpt)
          cons_qp ( 3, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 3, iquadpt)
          cons_qp ( 4, iquadpt) = prim_qp ( 1, iquadpt)
     1                          * prim_qp ( 4, iquadpt)

          cons_qp ( 5, iquadpt) = prim_qp ( 5, iquadpt) / (gamma - 1.0)
     1     + 0.5 * prim_qp ( 1, iquadpt) * ( prim_qp ( 2, iquadpt)**2
     1           + prim_qp ( 3, iquadpt)**2 + prim_qp ( 4, iquadpt)**2)
     1     + ( prim_qp ( 6, iquadpt)**2 + prim_qp ( 7, iquadpt)**2
     1       + prim_qp ( 8, iquadpt)**2) / ( 8.0 * pi)

          cons_qp ( 6, iquadpt) = prim_qp ( 6, iquadpt)
          cons_qp ( 7, iquadpt) = prim_qp ( 7, iquadpt)
          cons_qp ( 8, iquadpt) = prim_qp ( 8, iquadpt)

#if ( NFLUID >= 1)
          cons_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1             prim_qp ( 1, iquadpt)
     1           * prim_qp ( NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
#endif

          END DO
! End of "DO iquadpt ="

! ------------------------------

! Now that the conserved variables are available at all the quadrature
! points, we make the zone-centered conserved variables.

! For good measure, we make all their moments too,
! though they may not be needed.

! ----------

#if ( NUMORDER == 2)

          cons_modal ( :, 1) = ( cons_qp ( :, 01) + cons_qp ( :, 02)
     1                         + cons_qp ( :, 03) + cons_qp ( :, 04)
     1                         + cons_qp ( :, 05) + cons_qp ( :, 06))
     1                         / 6.0

          cons_modal ( :, 2) = cons_qp ( :, 01) - cons_qp ( :, 02)
          cons_modal ( :, 3) = cons_qp ( :, 03) - cons_qp ( :, 04)
          cons_modal ( :, 4) = cons_qp ( :, 05) - cons_qp ( :, 06)

#endif

! ----------

#if ( NUMORDER == 3)

      t ( :, 01) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 04) + cons_qp ( :, 05)

      t ( :, 02) = cons_qp ( :, 06) + cons_qp ( :, 07)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 03) = cons_qp ( :, 10) + cons_qp ( :, 11)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 04) = cons_qp ( :, 02) + cons_qp ( :, 06)
     1           + cons_qp ( :, 10) + cons_qp ( :, 12)

      t ( :, 05) = cons_qp ( :, 03) + cons_qp ( :, 07)
     1           + cons_qp ( :, 11) + cons_qp ( :, 13)

      t ( :, 06) = cons_qp ( :, 04) + cons_qp ( :, 05)
     1           + cons_qp ( :, 08) + cons_qp ( :, 09)

      t ( :, 07) = cons_qp ( :, 04) + cons_qp ( :, 08)
     1           + cons_qp ( :, 10) + cons_qp ( :, 11)

      t ( :, 08) = cons_qp ( :, 05) + cons_qp ( :, 09)
     1           + cons_qp ( :, 12) + cons_qp ( :, 13)

      t ( :, 09) = cons_qp ( :, 02) + cons_qp ( :, 03)
     1           + cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 02) = 0.25 * ( t ( :, 01) - t ( :, 02) )

      cons_modal ( :, 03) = 0.25 * ( t ( :, 04) - t ( :, 05) )

      cons_modal ( :, 04) = 0.25 * ( t ( :, 07) - t ( :, 08) )

      cons_modal ( :, 05) = 0.5 * ( t ( :, 01) + t ( :, 02)
     1 - t ( :, 03) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 06) = 0.5 * ( t ( :, 04) + t ( :, 05)
     1 - t ( :, 06) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 07) = 0.5 * ( t ( :, 07) + t ( :, 08)
     1 - t ( :, 09) ) - 2.0 * cons_qp ( :, 01)

      cons_modal ( :, 08) = cons_qp ( :, 02) - cons_qp ( :, 03)
     1 - cons_qp ( :, 06) + cons_qp ( :, 07)

      cons_modal ( :, 09) = cons_qp ( :, 10) - cons_qp ( :, 11)
     1 - cons_qp ( :, 12) + cons_qp ( :, 13)

      cons_modal ( :, 10) = cons_qp ( :, 04) - cons_qp ( :, 05)
     1 - cons_qp ( :, 08) + cons_qp ( :, 09)

      cons_modal ( :, 01) = cons_qp ( :, 01) + ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 12.0

#endif

! ----------

#if ( NUMORDER == 4)

! Note that in this formulation "cons_modal ( :, 14:19)" require
! space-time nodal points. They, therefore, can't be built here.

      cons_modal ( :, 05) = 2.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 05) )

      cons_modal ( :, 06) = 2.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 09) )

      cons_modal ( :, 07) = 2.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 01) + cons_qp ( :, 13) )

      cons_modal ( :, 08) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 15)
     1 + cons_qp ( :, 16) - cons_qp ( :, 17) - cons_qp ( :, 18)
     1 + cons_qp ( :, 19) - cons_qp ( :, 20) + cons_qp ( :, 21) )

      cons_modal ( :, 09) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 16)
     1 + cons_qp ( :, 18) - cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 17) - cons_qp ( :, 19) + cons_qp ( :, 21) )

      cons_modal ( :, 10) = 0.5 * ( cons_qp ( :, 14) - cons_qp ( :, 18)
     1 + cons_qp ( :, 15) - cons_qp ( :, 19) - cons_qp ( :, 16)
     1 + cons_qp ( :, 20) - cons_qp ( :, 17) + cons_qp ( :, 21) )

      cons_modal ( :, 11) = 16.0 * ( cons_qp ( :, 02)
     1 - 2.0 * cons_qp ( :, 03) + 2.0 * cons_qp ( :, 04)
     1 - cons_qp ( :, 05) ) / 3.0

      cons_modal ( :, 12) = 16.0 * ( cons_qp ( :, 06)
     1 - 2.0 * cons_qp ( :, 07) + 2.0 * cons_qp ( :, 08)
     1 - cons_qp ( :, 09) ) / 3.0

      cons_modal ( :, 13) = 16.0 * ( cons_qp ( :, 10)
     1 - 2.0 * cons_qp ( :, 11) + 2.0 * cons_qp ( :, 12)
     1 - cons_qp ( :, 13) ) / 3.0

      t ( :, 01) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 16) + cons_qp ( :, 17)

      t ( :, 02) = cons_qp ( :, 18) + cons_qp ( :, 19)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      t ( :, 03) = cons_qp ( :, 14) + cons_qp ( :, 16)
     1 + cons_qp ( :, 18) + cons_qp ( :, 20)

      t ( :, 04) = cons_qp ( :, 15) + cons_qp ( :, 17)
     1 + cons_qp ( :, 19) + cons_qp ( :, 21)

      t ( :, 05) = cons_qp ( :, 14) + cons_qp ( :, 15)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19)

      t ( :, 06) = cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 20) + cons_qp ( :, 21)

      cons_modal ( :, 02) = 2.0 * ( cons_qp ( :, 02)
     1 - cons_qp ( :, 05) ) / 3.0 - 0.1 * cons_modal ( :, 11)
     1 + 0.25 * ( t ( :, 01) - t ( :, 02) ) / 3.0

      cons_modal ( :, 03) = 2.0 * ( cons_qp ( :, 06)
     1 - cons_qp ( :, 09) ) / 3.0 - 0.1 * cons_modal ( :, 12)
     1 + 0.25 * ( t ( :, 03) - t ( :, 04) ) / 3.0

      cons_modal ( :, 04) = 2.0 * ( cons_qp ( :, 10)
     1 - cons_qp ( :, 13) ) / 3.0 - 0.1 * cons_modal ( :, 13)
     1 + 0.25 * ( t ( :, 05) - t ( :, 06) ) / 3.0

      cons_modal ( :, 20) = cons_qp ( :, 14) - cons_qp ( :, 18)
     1 - cons_qp ( :, 16) + cons_qp ( :, 20) - cons_qp ( :, 15)
     1 + cons_qp ( :, 19) + cons_qp ( :, 17) - cons_qp ( :, 21)

      cons_modal ( :, 01) = 0.125 * ( cons_qp ( :, 14)
     1 + cons_qp ( :, 15) + cons_qp ( :, 16) + cons_qp ( :, 17)
     1 + cons_qp ( :, 18) + cons_qp ( :, 19) + cons_qp ( :, 20)
     1 + cons_qp ( :, 21) ) - 0.5 * ( cons_modal ( :, 05)
     1 + cons_modal ( :, 06) + cons_modal ( :, 07) ) / 3.0

#endif

! ------------------------------

! Now store "diff_cons ( :)" between "cons_modal ( :, 1)" and the code's
! conserved variable mesh function "u". While these are zone-centered
! components, the magnetic field components will be overwritten.

          diff_cons ( :) = ABS ( u ( ixx, iyy, izz, :, 1)
     1                         - cons_modal ( :, 1) )

! ----------------------------------------------------------------------

! Make all the x-magnetic field variables at the top x-face-centered
! quadrature points.
! The top x-face-centered quadrature points in the computational domain
! are stored in "xf_qp".
! The x-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 6, :)".

! --------------------

          xf_qp ( 1, :) = xface_quadpt ( 1, :) * dx + xc
          xf_qp ( 2, :) = xface_quadpt ( 2, :) * dy + yc
          xf_qp ( 3, :) = xface_quadpt ( 3, :) * dz + zc

! ------------------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, xf_qp ( 1, iquadpt), xf_qp ( 2, iquadpt),
     1    xf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 6, 01) = 0.25 * ( prim_qp ( 6, 01)
     1       + prim_qp ( 6, 02) + prim_qp ( 6, 03) + prim_qp ( 6, 04))

          cons_modal ( 6, 02) = prim_qp ( 6, 01) - prim_qp ( 6, 02)

          cons_modal ( 6, 03) = prim_qp ( 6, 03) - prim_qp ( 6, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 03)

          cons_modal ( 6, 03) = prim_qp ( 6, 04) - prim_qp ( 6, 05)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 03)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 04) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 06) = prim_qp ( 6, 06) - prim_qp ( 6, 07)
     1       - prim_qp ( 6, 08) + prim_qp ( 6, 09)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 6, 04) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 02) + 2.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 05) = - 4.0 * prim_qp ( 6, 01)
     1       + 2.0 * prim_qp ( 6, 06) + 2.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 06) = prim_qp ( 6, 10) - prim_qp ( 6, 11)
     1       - prim_qp ( 6, 12) + prim_qp ( 6, 13)

          cons_modal ( 6, 07) = ( 16.0 * prim_qp ( 6, 02)
     1       - 16.0 * prim_qp ( 6, 05) - 32.0 * prim_qp ( 6, 03)
     1       + 32.0 * prim_qp ( 6, 04) ) / 3.0

          cons_modal ( 6, 08) = ( 16.0 * prim_qp ( 6, 06)
     1       - 16.0 * prim_qp ( 6, 09) - 32.0 * prim_qp ( 6, 07)
     1       + 32.0 * prim_qp ( 6, 08) ) / 3.0

          cons_modal ( 6, 09) = 2.0 * prim_qp ( 6, 10)
     1       + 2.0 * prim_qp ( 6, 11) - 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 06)
     1       + 4.0 * prim_qp ( 6, 09)

          cons_modal ( 6, 10) = 2.0 * prim_qp ( 6, 10)
     1       - 2.0 * prim_qp ( 6, 11) + 2.0 * prim_qp ( 6, 12)
     1       - 2.0 * prim_qp ( 6, 13) - 4.0 * prim_qp ( 6, 02)
     1       + 4.0 * prim_qp ( 6, 05)

          cons_modal ( 6, 01) = prim_qp ( 6, 01)
     1       + ( cons_modal ( 6, 04) + cons_modal ( 6, 05) ) / 12.0

          cons_modal ( 6, 02) = prim_qp ( 6, 02) - prim_qp ( 6, 05)
     1       - 0.1 * cons_modal ( 6, 07) + cons_modal ( 6, 10) / 12.0

          cons_modal ( 6, 03) = prim_qp ( 6, 06) - prim_qp ( 6, 09)
     1       - 0.1 * cons_modal ( 6, 08) + cons_modal ( 6, 09) / 12.0

#endif
! ----------

          diff_cons ( 6) = ABS ( bfxbdy ( ixx, iyy, izz, 1)
     1                         - cons_modal ( 6, 1) )

! ----------------------------------------------------------------------

! Make all the y-magnetic field variables at the top y-face-centered
! quadrature points.
! The top y-face-centered quadrature points in the computational domain
! are stored in "yf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 7, :)".

! --------------------

          yf_qp ( 1, :) = yface_quadpt ( 1, :) * dx + xc
          yf_qp ( 2, :) = yface_quadpt ( 2, :) * dy + yc
          yf_qp ( 3, :) = yface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, yf_qp ( 1, iquadpt), yf_qp ( 2, iquadpt),
     1    yf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 7, 01) = 0.25 * ( prim_qp ( 7, 01)
     1       + prim_qp ( 7, 02) + prim_qp ( 7, 03) + prim_qp ( 7, 04))

          cons_modal ( 7, 02) = prim_qp ( 7, 01) - prim_qp ( 7, 02)

          cons_modal ( 7, 03) = prim_qp ( 7, 03) - prim_qp ( 7, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 03)

          cons_modal ( 7, 03) = prim_qp ( 7, 04) - prim_qp ( 7, 05)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 03)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 04) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 06) = prim_qp ( 7, 06) - prim_qp ( 7, 07)
     1       - prim_qp ( 7, 08) + prim_qp ( 7, 09)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 7, 04) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 02) + 2.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 05) = - 4.0 * prim_qp ( 7, 01)
     1       + 2.0 * prim_qp ( 7, 06) + 2.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 06) = prim_qp ( 7, 10) - prim_qp ( 7, 11)
     1       - prim_qp ( 7, 12) + prim_qp ( 7, 13)

          cons_modal ( 7, 07) = ( 16.0 * prim_qp ( 7, 02)
     1       - 16.0 * prim_qp ( 7, 05) - 32.0 * prim_qp ( 7, 03)
     1       + 32.0 * prim_qp ( 7, 04) ) / 3.0

          cons_modal ( 7, 08) = ( 16.0 * prim_qp ( 7, 06)
     1       - 16.0 * prim_qp ( 7, 09) - 32.0 * prim_qp ( 7, 07)
     1       + 32.0 * prim_qp ( 7, 08) ) / 3.0

          cons_modal ( 7, 09) = 2.0 * prim_qp ( 7, 10)
     1       + 2.0 * prim_qp ( 7, 11) - 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 06)
     1       + 4.0 * prim_qp ( 7, 09)

          cons_modal ( 7, 10) = 2.0 * prim_qp ( 7, 10)
     1       - 2.0 * prim_qp ( 7, 11) + 2.0 * prim_qp ( 7, 12)
     1       - 2.0 * prim_qp ( 7, 13) - 4.0 * prim_qp ( 7, 02)
     1       + 4.0 * prim_qp ( 7, 05)

          cons_modal ( 7, 01) = prim_qp ( 7, 01)
     1       + ( cons_modal ( 7, 04) + cons_modal ( 7, 05) ) / 12.0

          cons_modal ( 7, 02) = prim_qp ( 7, 02) - prim_qp ( 7, 05)
     1       - 0.1 * cons_modal ( 7, 07) + cons_modal ( 7, 10) / 12.0

          cons_modal ( 7, 03) = prim_qp ( 7, 06) - prim_qp ( 7, 09)
     1       - 0.1 * cons_modal ( 7, 08) + cons_modal ( 7, 09) / 12.0

#endif
! ----------

          diff_cons ( 7) = ABS ( bfybdy ( ixx, iyy, izz, 1)
     1                         - cons_modal ( 7, 1) )

! ----------------------------------------------------------------------

! Make all the z-magnetic field variables at the top z-face-centered
! quadrature points.
! The top z-face-centered quadrature points in the computational domain
! are stored in "zf_qp".
! The y-magnetic fields at the facial quadrature points are stored
! in "prim_qp ( 8, :)".

! --------------------

          zf_qp ( 1, :) = zface_quadpt ( 1, :) * dx + xc
          zf_qp ( 2, :) = zface_quadpt ( 2, :) * dy + yc
          zf_qp ( 3, :) = zface_quadpt ( 3, :) * dz + zc

! --------------------

          DO iquadpt = 1, NUM_FACECEN_QUADPTS

! ----------

          CALL ALFVEN3D_GENERATOR_MAKE_PRIM
     1  ( iproblem, zf_qp ( 1, iquadpt), zf_qp ( 2, iquadpt),
     1    zf_qp ( 3, iquadpt), time,
     1    base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1    base_grid_ymax, base_grid_zmin, base_grid_zmax,
     1    gamma, prim_qp ( 1, iquadpt) )

! ----------

          END DO
! End of "DO iquadpt ="

! --------------------
#if ( NUMORDER == 2)

          cons_modal ( 8, 01) = 0.25 * ( prim_qp ( 8, 01)
     1       + prim_qp ( 8, 02) + prim_qp ( 8, 03) + prim_qp ( 8, 04))

          cons_modal ( 8, 02) = prim_qp ( 8, 01) - prim_qp ( 8, 02)

          cons_modal ( 8, 03) = prim_qp ( 8, 03) - prim_qp ( 8, 04)

#endif
! ----------
#if ( NUMORDER == 3)

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 03)

          cons_modal ( 8, 03) = prim_qp ( 8, 04) - prim_qp ( 8, 05)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 03)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 04) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 06) = prim_qp ( 8, 06) - prim_qp ( 8, 07)
     1       - prim_qp ( 8, 08) + prim_qp ( 8, 09)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

#endif
! ----------
#if ( NUMORDER == 4)

          cons_modal ( 8, 04) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 02) + 2.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 05) = - 4.0 * prim_qp ( 8, 01)
     1       + 2.0 * prim_qp ( 8, 06) + 2.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 06) = prim_qp ( 8, 10) - prim_qp ( 8, 11)
     1       - prim_qp ( 8, 12) + prim_qp ( 8, 13)

          cons_modal ( 8, 07) = ( 16.0 * prim_qp ( 8, 02)
     1       - 16.0 * prim_qp ( 8, 05) - 32.0 * prim_qp ( 8, 03)
     1       + 32.0 * prim_qp ( 8, 04) ) / 3.0

          cons_modal ( 8, 08) = ( 16.0 * prim_qp ( 8, 06)
     1       - 16.0 * prim_qp ( 8, 09) - 32.0 * prim_qp ( 8, 07)
     1       + 32.0 * prim_qp ( 8, 08) ) / 3.0

          cons_modal ( 8, 09) = 2.0 * prim_qp ( 8, 10)
     1       + 2.0 * prim_qp ( 8, 11) - 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 06)
     1       + 4.0 * prim_qp ( 8, 09)

          cons_modal ( 8, 10) = 2.0 * prim_qp ( 8, 10)
     1       - 2.0 * prim_qp ( 8, 11) + 2.0 * prim_qp ( 8, 12)
     1       - 2.0 * prim_qp ( 8, 13) - 4.0 * prim_qp ( 8, 02)
     1       + 4.0 * prim_qp ( 8, 05)

          cons_modal ( 8, 01) = prim_qp ( 8, 01)
     1       + ( cons_modal ( 8, 04) + cons_modal ( 8, 05) ) / 12.0

          cons_modal ( 8, 02) = prim_qp ( 8, 02) - prim_qp ( 8, 05)
     1       - 0.1 * cons_modal ( 8, 07) + cons_modal ( 8, 10) / 12.0

          cons_modal ( 8, 03) = prim_qp ( 8, 06) - prim_qp ( 8, 09)
     1       - 0.1 * cons_modal ( 8, 08) + cons_modal ( 8, 09) / 12.0

#endif
! ----------

          diff_cons ( 8) = ABS ( bfzbdy ( ixx, iyy, izz, 1)
     1                         - cons_modal ( 8, 1) )

! ----------------------------------------------------------------------

! Now use "diff_cons ( :)" to build the "L_1_error and L_inf_error" arrays.

          L_1_error ( :) = L_1_error ( :) + diff_cons ( :)
          L_inf_error ( :) = AMAX1 ( L_inf_error ( :), diff_cons ( :))

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

! Write out the accuracy information.

      L_1_error ( :) = L_1_error ( :) / ( ix1 * iy1 * iz1)

      WRITE ( 6, *)
     1"---------------------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Accuracy analysis for iproblem = ", iproblem,
     1             " at time = ", time
      WRITE ( 6, *)"# of zones : ix1, iy1, iz1 = ", ix1, iy1, iz1


      DO ixx = 1, NFIELD + 1 + NFLUID
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Error in variable no. = ", ixx
      WRITE ( 6, *)"L_1_error, L_inf_error = ",
     1              L_1_error ( ixx), L_inf_error ( ixx)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"---------------------------------------------------------------"

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ERROR_ALFVEN3D

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ALFVEN3D_GENERATOR_MAKE_PRIM

     1   ( iproblem, x, y, z, time,
     1     xmin, xmax, ymin, ymax, zmin, zmax,
     1     gamma, prim)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! For higher order schemes, it is most convenient to just get a
! vector of primitive variables at any general space or time point 
! in the computational domain. This information can then be used to
! set the variables at zone-centers, faces etc.

! The variable names and meanings (and their input/output status) are given
! below:

! "iproblem" : Problem number that is being set up. (1: Alfven; 2:density pulse)
! "x, y, z, time" : space-time position at which variables are desired. (input)
! "xmin, xmax, ymin, ymax, zmin, zmax" : Domain boundary in 3d. (input)
! "gamma" : gamma for MHD and Euler fluids. (input)

! "prim ( NFIELD + 1 + NFLUID)" : Vector of primitive variables. (output)
!                  ( Order : mhd density, mhd x,y,z-velocities, mhd pressure,
!                    bfx,y,z, mhd species fraction;
!                    euler density, euler x,y,z-velocities, euler pressure,
!                    euler species fraction; repeat for more euler families)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iproblem

      REAL x, y, z, time, gamma, prim ( NFIELD + 1 + NFLUID),
     1     xmin, xmax, ymin, ymax, zmin, zmax

! ----------------------------------------------------------------------

      INTEGER i, j

      REAL pi, b0, v0, rho0, prs0, eps, va0, lambda, phase,
     1     vx, vy, vz, bx, by, bz, tempaa,
     1     A_x, A_y, A_z, a0_x, a0_y, a0_z,
     1     phi, theta, psi, a ( 3, 3), ainv ( 3, 3)

! ----------------------------------------------------------------------

! Set up some variables here.

      pi = 4.0 * ATAN ( 1.0)

      IF ( NFIELD + 1 + NFLUID .NE. 8) THEN
      WRITE ( 6, *)"Stopping in ALFVEN3D_GENERATOR_MAKE_PRIM, NVAR = ",
     1              NFIELD + 1 + NFLUID
      STOP
      END IF

! ------------------------------

      IF ( iproblem .EQ. 1) THEN

! Alfven wave problem set up here.

! Set density and pressure.
      rho0 = 1.0
      prs0 = 10.0

! Advection speed.
      v0 = 1.0

! Set unperturbed magnetic field so that Alfven speed equals advection speed.
      b0 = v0 * SQRT ( 4.0 * pi * rho0)

! Amplitude of the perturbation.
      eps = 0.02

! ----------

      ELSE

! ----------

! Density wave problem set up here.

! Set density and pressure.
      rho0 = 1.0
      prs0 = 1.0

! Advection speed.
      v0 = 1.0

! Magnetic field set to zero.
      b0 = 0.0

! Amplitude of the perturbation.
      eps = 0.2

      END IF

! ------------------------------

! Alfven speed
      va0 = ABS ( b0) / SQRT ( 4.0 * pi * rho0)

! Wavelength "lambda" of the wave is set to equal "1.0 / SQRT (3.0)" for the
! unit cube.
! This choice of wavelength makes the wave periodic in all directions.
      lambda = 1.0 / SQRT (3.0)

! ----------------------------------------------------------------------

! Make the angles for the rotation matrix. This enables us to rotate the
! vectorial entities that are needed in setting up the wave. Thus an
! Alfven wave that is originally propagating along the x-axis is taken and the
! computational domain is rotated so that the wave now propagates along
! the diagonal.
      phi = - pi / 4.0
      theta = ASIN ( - SQRT ( 2.0) / SQRT ( 3.0))
      psi = ASIN ( 0.25 * ( SQRT ( 2.0) - SQRT ( 6.0) ) )

! Now make the rotation matrix and its inverse.

      a ( 1, 1) = COS ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * SIN ( psi)
      a ( 1, 2) = COS ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * SIN ( psi)
      a ( 1, 3) = SIN ( psi) * SIN ( theta)

      a ( 2, 1) = - SIN ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * COS ( psi)
      a ( 2, 2) = - SIN ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * COS ( psi)
      a ( 2, 3) = COS ( psi) * SIN ( theta)

      a ( 3, 1) = SIN ( theta) * SIN ( phi)
      a ( 3, 2) = - SIN ( theta) * COS ( phi)
      a ( 3, 3) = COS ( theta)

! Take the transpose of the above rotation matrix. That is also its inverse.
      DO i = 1, 3
      DO j = 1, 3
      ainv ( i, j) = a ( j, i)
      END DO
      END DO

! ----------------------------------------------------------------------

      phase = (2.0 * pi / lambda)
     1      * ( ainv ( 1, 1) * x + ainv ( 1, 2) * y + ainv ( 1, 3) * z
     1        - ( v0 + va0) * time )

! Set up the velocity components along the original x-axis, (before rotation).
      vx = v0
      vy = eps * COS ( phase)
      vz = eps * SIN ( phase)

! Set up the magnetic field components along the original x-axis
! (before rotation).
      tempaa = eps * SQRT ( 4.0 * pi * rho0)
      bx = b0
      by = - tempaa * COS ( phase)
      bz = - tempaa * SIN ( phase)

! Now make the vector of primitive variables. Remember to rotate the
! vectorial entities into the diagonal of the unit cube using the
! rotation matrix. Then return.

      IF ( iproblem .EQ. 1) THEN

      prim ( 1) = rho0
      prim ( 2) = a ( 1, 1) * vx + a ( 1, 2) * vy + a ( 1, 3) * vz
      prim ( 3) = a ( 2, 1) * vx + a ( 2, 2) * vy + a ( 2, 3) * vz
      prim ( 4) = a ( 3, 1) * vx + a ( 3, 2) * vy + a ( 3, 3) * vz
      prim ( 5) = prs0
      prim ( 6) = a ( 1, 1) * bx + a ( 1, 2) * by + a ( 1, 3) * bz
      prim ( 7) = a ( 2, 1) * bx + a ( 2, 2) * by + a ( 2, 3) * bz
      prim ( 8) = a ( 3, 1) * bx + a ( 3, 2) * by + a ( 3, 3) * bz

      ELSE

      prim ( 1) = rho0 * ( 1.0 + eps * SIN ( phase) )
      prim ( 2) = a ( 1, 1) * v0
      prim ( 3) = a ( 2, 1) * v0
      prim ( 4) = a ( 3, 1) * v0
      prim ( 5) = prs0
      prim ( 6) = 0.0
      prim ( 7) = 0.0
      prim ( 8) = 0.0

      END IF

! ----------------------------------------------------------------------

! For good measure, we also build the magnetic vector potential here, though
! it is not passed back through the subroutine interface of this routine.
! The original magnetic vector potential is given by "( a0_x, a0_y, a0_z)"
! and its rotated counterpart is given by "( A_x, A_y, A_z)".

      tempaa = lambda * eps * SQRT ( 4.0 * pi * rho0) / ( 2.0 * pi)
      a0_x = 0.0
      a0_y = tempaa * COS ( phase)
      a0_z = b0
     1     * ( ainv ( 2, 1) * x + ainv ( 2, 2) * y + ainv ( 2, 3) * z )
     1     + tempaa * SIN ( phase)

      A_x = a ( 1, 1) * a0_x + a ( 1, 2) * a0_y + a ( 1, 3) * a0_z
      A_y = a ( 2, 1) * a0_x + a ( 2, 2) * a0_y + a ( 2, 3) * a0_z
      A_z = a ( 3, 1) * a0_x + a ( 3, 2) * a0_y + a ( 3, 3) * a0_z

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ALFVEN3D_GENERATOR_MAKE_PRIM

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      REAL FUNCTION ALFVEN3D_A_x ( x, y, z)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! For higher order schemes, it is most convenient to just get a
! vector of primitive variables at any general space or time point 
! in the computational domain. This information can then be used to
! set the variables at zone-centers, faces etc.

! The variable names and meanings (and their input/output status) are given
! below:

! "x, y, z, time" : space-time position at which variables are desired. (input)
! "ALFVEN3D_A_x": Value of the rotated vector potential along the cube's x-axis.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL x, y, z

! ----------------------------------------------------------------------

      INTEGER i, j

      REAL pi, b0, v0, rho0, prs0, eps, va0, lambda, phase,
     1     vx, vy, vz, bx, by, bz, tempaa, time,
     1     A_x, A_y, A_z, a0_x, a0_y, a0_z,
     1     phi, theta, psi, a ( 3, 3), ainv ( 3, 3)

! ----------------------------------------------------------------------

! Set up some variables here.

      time = 0.0

      pi = 4.0 * ATAN ( 1.0)

      IF ( NFIELD + 1 + NFLUID .NE. 8) THEN
      WRITE ( 6, *)"Stopping in function ALFVEN3D_A_x, NVAR = ",
     1              NFIELD + 1 + NFLUID
      STOP
      END IF

! ------------------------------

! Alfven wave problem set up here.

! Set density and pressure.
      rho0 = 1.0
      prs0 = 10.0

! Advection speed.
      v0 = 1.0

! Set unperturbed magnetic field so that Alfven speed equals advection speed.
      b0 = v0 * SQRT ( 4.0 * pi * rho0)

! Amplitude of the perturbation.
      eps = 0.02

! ------------------------------

! Alfven speed
      va0 = ABS ( b0) / SQRT ( 4.0 * pi * rho0)

! Wavelength "lambda" of the wave is set to equal "1.0 / SQRT (3.0)" for the
! unit cube.
! This choice of wavelength makes the wave periodic in all directions.
      lambda = 1.0 / SQRT (3.0)

! ----------------------------------------------------------------------

! Make the angles for the rotation matrix. This enables us to rotate the
! vectorial entities that are needed in setting up the wave. Thus an
! Alfven wave that is originally propagating along the x-axis is taken and the
! computational domain is rotated so that the wave now propagates along
! the diagonal.
      phi = - pi / 4.0
      theta = ASIN ( - SQRT ( 2.0) / SQRT ( 3.0))
      psi = ASIN ( 0.25 * ( SQRT ( 2.0) - SQRT ( 6.0) ) )

! Now make the rotation matrix and its inverse.

      a ( 1, 1) = COS ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * SIN ( psi)
      a ( 1, 2) = COS ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * SIN ( psi)
      a ( 1, 3) = SIN ( psi) * SIN ( theta)

      a ( 2, 1) = - SIN ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * COS ( psi)
      a ( 2, 2) = - SIN ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * COS ( psi)
      a ( 2, 3) = COS ( psi) * SIN ( theta)

      a ( 3, 1) = SIN ( theta) * SIN ( phi)
      a ( 3, 2) = - SIN ( theta) * COS ( phi)
      a ( 3, 3) = COS ( theta)

! Take the transpose of the above rotation matrix. That is also its inverse.
      DO i = 1, 3
      DO j = 1, 3
      ainv ( i, j) = a ( j, i)
      END DO
      END DO

! ----------------------------------------------------------------------

      phase = (2.0 * pi / lambda)
     1      * ( ainv ( 1, 1) * x + ainv ( 1, 2) * y + ainv ( 1, 3) * z
     1        - ( v0 + va0) * time )

! ----------------------------------------------------------------------

! We build the magnetic vector potential here.
! The original magnetic vector potential is given by "( a0_x, a0_y, a0_z)"
! and its rotated counterpart is given by "( A_x, A_y, A_z)".

      tempaa = lambda * eps * SQRT ( 4.0 * pi * rho0) / ( 2.0 * pi)
      a0_x = 0.0
      a0_y = tempaa * COS ( phase)
      a0_z = b0
     1     * ( ainv ( 2, 1) * x + ainv ( 2, 2) * y + ainv ( 2, 3) * z )
     1     + tempaa * SIN ( phase)

      A_x = a ( 1, 1) * a0_x + a ( 1, 2) * a0_y + a ( 1, 3) * a0_z
      A_y = a ( 2, 1) * a0_x + a ( 2, 2) * a0_y + a ( 2, 3) * a0_z
      A_z = a ( 3, 1) * a0_x + a ( 3, 2) * a0_y + a ( 3, 3) * a0_z

      ALFVEN3D_A_x = A_x

! ----------------------------------------------------------------------

      RETURN

      END FUNCTION ALFVEN3D_A_x

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      REAL FUNCTION ALFVEN3D_A_y ( x, y, z)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! For higher order schemes, it is most convenient to just get a
! vector of primitive variables at any general space or time point 
! in the computational domain. This information can then be used to
! set the variables at zone-centers, faces etc.

! The variable names and meanings (and their input/output status) are given
! below:

! "x, y, z, time" : space-time position at which variables are desired. (input)
! "ALFVEN3D_A_y": Value of the rotated vector potential along the cube's y-axis.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL x, y, z

! ----------------------------------------------------------------------

      INTEGER i, j

      REAL pi, b0, v0, rho0, prs0, eps, va0, lambda, phase,
     1     vx, vy, vz, bx, by, bz, tempaa, time,
     1     A_x, A_y, A_z, a0_x, a0_y, a0_z,
     1     phi, theta, psi, a ( 3, 3), ainv ( 3, 3)

! ----------------------------------------------------------------------

! Set up some variables here.

      time = 0.0

      pi = 4.0 * ATAN ( 1.0)

      IF ( NFIELD + 1 + NFLUID .NE. 8) THEN
      WRITE ( 6, *)"Stopping in function ALFVEN3D_A_y, NVAR = ",
     1              NFIELD + 1 + NFLUID
      STOP
      END IF

! ------------------------------

! Alfven wave problem set up here.

! Set density and pressure.
      rho0 = 1.0
      prs0 = 10.0

! Advection speed.
      v0 = 1.0

! Set unperturbed magnetic field so that Alfven speed equals advection speed.
      b0 = v0 * SQRT ( 4.0 * pi * rho0)

! Amplitude of the perturbation.
      eps = 0.02

! ------------------------------

! Alfven speed
      va0 = ABS ( b0) / SQRT ( 4.0 * pi * rho0)

! Wavelength "lambda" of the wave is set to equal "1.0 / SQRT (3.0)" for the
! unit cube.
! This choice of wavelength makes the wave periodic in all directions.
      lambda = 1.0 / SQRT (3.0)

! ----------------------------------------------------------------------

! Make the angles for the rotation matrix. This enables us to rotate the
! vectorial entities that are needed in setting up the wave. Thus an
! Alfven wave that is originally propagating along the x-axis is taken and the
! computational domain is rotated so that the wave now propagates along
! the diagonal.
      phi = - pi / 4.0
      theta = ASIN ( - SQRT ( 2.0) / SQRT ( 3.0))
      psi = ASIN ( 0.25 * ( SQRT ( 2.0) - SQRT ( 6.0) ) )

! Now make the rotation matrix and its inverse.

      a ( 1, 1) = COS ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * SIN ( psi)
      a ( 1, 2) = COS ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * SIN ( psi)
      a ( 1, 3) = SIN ( psi) * SIN ( theta)

      a ( 2, 1) = - SIN ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * COS ( psi)
      a ( 2, 2) = - SIN ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * COS ( psi)
      a ( 2, 3) = COS ( psi) * SIN ( theta)

      a ( 3, 1) = SIN ( theta) * SIN ( phi)
      a ( 3, 2) = - SIN ( theta) * COS ( phi)
      a ( 3, 3) = COS ( theta)

! Take the transpose of the above rotation matrix. That is also its inverse.
      DO i = 1, 3
      DO j = 1, 3
      ainv ( i, j) = a ( j, i)
      END DO
      END DO

! ----------------------------------------------------------------------

      phase = (2.0 * pi / lambda)
     1      * ( ainv ( 1, 1) * x + ainv ( 1, 2) * y + ainv ( 1, 3) * z
     1        - ( v0 + va0) * time )

! ----------------------------------------------------------------------

! We build the magnetic vector potential here.
! The original magnetic vector potential is given by "( a0_x, a0_y, a0_z)"
! and its rotated counterpart is given by "( A_x, A_y, A_z)".

      tempaa = lambda * eps * SQRT ( 4.0 * pi * rho0) / ( 2.0 * pi)
      a0_x = 0.0
      a0_y = tempaa * COS ( phase)
      a0_z = b0
     1     * ( ainv ( 2, 1) * x + ainv ( 2, 2) * y + ainv ( 2, 3) * z )
     1     + tempaa * SIN ( phase)

      A_x = a ( 1, 1) * a0_x + a ( 1, 2) * a0_y + a ( 1, 3) * a0_z
      A_y = a ( 2, 1) * a0_x + a ( 2, 2) * a0_y + a ( 2, 3) * a0_z
      A_z = a ( 3, 1) * a0_x + a ( 3, 2) * a0_y + a ( 3, 3) * a0_z

      ALFVEN3D_A_y = A_y

! ----------------------------------------------------------------------

      RETURN

      END FUNCTION ALFVEN3D_A_y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      REAL FUNCTION ALFVEN3D_A_z ( x, y, z)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! For higher order schemes, it is most convenient to just get a
! vector of primitive variables at any general space or time point 
! in the computational domain. This information can then be used to
! set the variables at zone-centers, faces etc.

! The variable names and meanings (and their input/output status) are given
! below:

! "x, y, z, time" : space-time position at which variables are desired. (input)
! "ALFVEN3D_A_z": Value of the rotated vector potential along the cube's z-axis.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL x, y, z

! ----------------------------------------------------------------------

      INTEGER i, j

      REAL pi, b0, v0, rho0, prs0, eps, va0, lambda, phase,
     1     vx, vy, vz, bx, by, bz, tempaa, time,
     1     A_x, A_y, A_z, a0_x, a0_y, a0_z,
     1     phi, theta, psi, a ( 3, 3), ainv ( 3, 3)

! ----------------------------------------------------------------------

! Set up some variables here.

      time = 0.0

      pi = 4.0 * ATAN ( 1.0)

      IF ( NFIELD + 1 + NFLUID .NE. 8) THEN
      WRITE ( 6, *)"Stopping in function ALFVEN3D_A_z, NVAR = ",
     1              NFIELD + 1 + NFLUID
      STOP
      END IF

! ------------------------------

! Alfven wave problem set up here.

! Set density and pressure.
      rho0 = 1.0
      prs0 = 10.0

! Advection speed.
      v0 = 1.0

! Set unperturbed magnetic field so that Alfven speed equals advection speed.
      b0 = v0 * SQRT ( 4.0 * pi * rho0)

! Amplitude of the perturbation.
      eps = 0.02

! ------------------------------

! Alfven speed
      va0 = ABS ( b0) / SQRT ( 4.0 * pi * rho0)

! Wavelength "lambda" of the wave is set to equal "1.0 / SQRT (3.0)" for the
! unit cube.
! This choice of wavelength makes the wave periodic in all directions.
      lambda = 1.0 / SQRT (3.0)

! ----------------------------------------------------------------------

! Make the angles for the rotation matrix. This enables us to rotate the
! vectorial entities that are needed in setting up the wave. Thus an
! Alfven wave that is originally propagating along the x-axis is taken and the
! computational domain is rotated so that the wave now propagates along
! the diagonal.
      phi = - pi / 4.0
      theta = ASIN ( - SQRT ( 2.0) / SQRT ( 3.0))
      psi = ASIN ( 0.25 * ( SQRT ( 2.0) - SQRT ( 6.0) ) )

! Now make the rotation matrix and its inverse.

      a ( 1, 1) = COS ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * SIN ( psi)
      a ( 1, 2) = COS ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * SIN ( psi)
      a ( 1, 3) = SIN ( psi) * SIN ( theta)

      a ( 2, 1) = - SIN ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * COS ( psi)
      a ( 2, 2) = - SIN ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * COS ( psi)
      a ( 2, 3) = COS ( psi) * SIN ( theta)

      a ( 3, 1) = SIN ( theta) * SIN ( phi)
      a ( 3, 2) = - SIN ( theta) * COS ( phi)
      a ( 3, 3) = COS ( theta)

! Take the transpose of the above rotation matrix. That is also its inverse.
      DO i = 1, 3
      DO j = 1, 3
      ainv ( i, j) = a ( j, i)
      END DO
      END DO

! ----------------------------------------------------------------------

      phase = (2.0 * pi / lambda)
     1      * ( ainv ( 1, 1) * x + ainv ( 1, 2) * y + ainv ( 1, 3) * z
     1        - ( v0 + va0) * time )

! ----------------------------------------------------------------------

! We build the magnetic vector potential here.
! The original magnetic vector potential is given by "( a0_x, a0_y, a0_z)"
! and its rotated counterpart is given by "( A_x, A_y, A_z)".

      tempaa = lambda * eps * SQRT ( 4.0 * pi * rho0) / ( 2.0 * pi)
      a0_x = 0.0
      a0_y = tempaa * COS ( phase)
      a0_z = b0
     1     * ( ainv ( 2, 1) * x + ainv ( 2, 2) * y + ainv ( 2, 3) * z )
     1     + tempaa * SIN ( phase)

      A_x = a ( 1, 1) * a0_x + a ( 1, 2) * a0_y + a ( 1, 3) * a0_z
      A_y = a ( 2, 1) * a0_x + a ( 2, 2) * a0_y + a ( 2, 3) * a0_z
      A_z = a ( 3, 1) * a0_x + a ( 3, 2) * a0_y + a ( 3, 3) * a0_z

      ALFVEN3D_A_z = A_z

! ----------------------------------------------------------------------

      RETURN

      END FUNCTION ALFVEN3D_A_z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)

! This is just a tester code for testing whether our rotation matrices
! are working well. We need such a rotation matrix to rotate the waves
! into the diagonal of the unit cube.

      IMPLICIT NONE

! Define the rotation matrix from eq.(4-46) of Goldstein, 2nd edition.
      REAL phi, theta, psi

      REAL pi, a ( 3, 3), ainv ( 3, 3), achk ( 3, 3)

      INTEGER i, j


      pi = 4.0 * ATAN ( 1.0)

      phi = - pi / 4.0

      theta = ASIN ( - SQRT ( 2.0) / SQRT ( 3.0))

      psi = ASIN ( 0.25 * ( SQRT ( 2.0) - SQRT ( 6.0) ) )

      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "

      WRITE ( 6, *)"phi, theta, psi = ", phi, theta, psi
      WRITE ( 6, *)"phi, theta, psi (in degrees) = ",
     1    phi*(180.0/pi), theta*(180.0/pi), psi*(180.0/pi)

      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "

      a ( 1, 1) = COS ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * SIN ( psi)
      a ( 1, 2) = COS ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * SIN ( psi)
      a ( 1, 3) = SIN ( psi) * SIN ( theta)

      a ( 2, 1) = - SIN ( psi) * COS ( phi)
     1          - COS ( theta) * SIN ( phi) * COS ( psi)
      a ( 2, 2) = - SIN ( psi) * SIN ( phi)
     1          + COS ( theta) * COS ( phi) * COS ( psi)
      a ( 2, 3) = COS ( psi) * SIN ( theta)

      a ( 3, 1) = SIN ( theta) * SIN ( phi)
      a ( 3, 2) = - SIN ( theta) * COS ( phi)
      a ( 3, 3) = COS ( theta)

! Take the transpose of the above rotation matrix. That is also its inverse.
      DO i = 1, 3
      DO j = 1, 3
      ainv ( i, j) = a ( j, i)
      END DO
      END DO

      achk ( :, :) = 0.0

      DO i = 1, 3
      DO j = 1, 3
      achk ( i, j) = ainv ( i, 1) * a ( 1, j) + ainv ( i, 2) * a ( 2, j)
     1             + ainv ( i, 3) * a ( 3, j)
      WRITE ( 6, *)"achk ( ",i,", ",j," ) = ", achk ( i, j)
      END DO
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "

      DO i = 1, 3
      DO j = 1, 3
      WRITE ( 6, *)"a ( ",i,", ",j," ) = ", a ( i, j)
      END DO
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "

      DO i = 1, 3
      DO j = 1, 3
      WRITE ( 6, *)"ainv ( ",i,", ",j," ) = ", ainv ( i, j)
      END DO
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "


      STOP

      END

#endif

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

          tempac = cfl_coef * tempac / ( tempaa + tempab)

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

          tempad = cfl_coef
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

          tempac = cfl_coef * tempac / ( tempaa + tempab)

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
!      eng_tmp = 116.5 / ( gamma - 1.0) + 0.5 * 8.0 * 8.25**2
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
! This call should be made once in each timestep loop. Best place
! to put it is after the padding of the ghost zones but before the
! call to the limiter.
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

          tempac = cfl_coef * tempac / ( tempaa + tempab)

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

! Zero out any stray magnetic field buildup.

      bfxbdy ( :, :, :, 1) = 0.0

      bfybdy ( :, :, :, 1) = 0.0

      bfzbdy ( :, :, :, 1) = 0.0

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

          tempac = cfl_coef * tempac / ( tempaa + tempab)

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

! Builds two variants of two-dimensional Riemann problems drawn from
! Brio, Zakharian & Webb (2001) JCP. The problems are set up in the XY
! plane with 400x400x1 zones spanning [-1,1]x[-1,1]x[-1,1].
! Set "gamma = 1.4" and keep the boundaries at continuitive (no. 2).
! Such 2d Riemann problems consist of four constant states in the 4 quadrants.
! The first problem has two initial weak shocks and two slip lines and
! is stopped at "timestop = 0.52". The second problem produces double
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

      iproblem = 2

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

          tempac = cfl_coef * tempac / ( tempaa + tempab)

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
     1     bfxuse, bfyuse, bfzuse,

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

          ELSE IF ( iy1 .GT. 1) THEN

          rhouse = v ( iyy, 1)

          vlxuse = v ( iyy, 4)
          vlyuse = v ( iyy, 2)
          vlzuse = v ( iyy, 3)

          prsuse = v ( iyy, 5)

          bfxuse = v ( iyy, 7)
          bfyuse = bfxbar ( iyy)
          bfzuse = v ( iyy, 6)

          ELSE

          rhouse = v ( izz, 1)

          vlxuse = v ( izz, 3)
          vlyuse = v ( izz, 4)
          vlzuse = v ( izz, 2)

          prsuse = v ( izz, 5)

          bfxuse = v ( izz, 6)
          bfyuse = v ( izz, 7)
          bfzuse = bfxbar ( izz)

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

          tempac = cfl_coef * tempac / ( tempaa + tempab)

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

      ELSE

      v ( ix, 1) = rhore_r
      v ( ix, 2) = vlxre_r
      v ( ix, 3) = vlyre_r
      v ( ix, 4) = vlzre_r
      v ( ix, 5) = prsre_r
      v ( ix, 6) = bfyre_r
      v ( ix, 7) = bfzre_r

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




