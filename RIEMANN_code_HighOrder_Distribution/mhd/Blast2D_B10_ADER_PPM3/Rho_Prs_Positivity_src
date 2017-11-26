! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! If one also wants to bracket the velocity and magnetic field range, then
! set it to 1.
#define REDUCE_VEL_BF 0

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE POSITIVE_RHO_PRS ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Restores density and pressure positivity to the reconstruction that
! was carried out in "LIMIT_ALL_MOMENTS". For many problems this is not
! needed at all. For some problems with very strong shocks and rarefactions,
! this may be needed once in a while in a few zones for an extraordinarily
! small number of timesteps.

! It is best not to put any restrictions on the velocities or magnetic
! fields, even though "RECONST_POS_RHO_PRS_MHD_CUBIC" has the means to
! do so. In other words, it is best to set "REDUCE_VEL_BF" to zero.

! This is the subroutine where the "entropygr" variable is initialized
! for this timestep. Its entire space-time variation within the zone
! (i.e. the equivalent of a second order accurate ADER step) is also
! built here.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe

      REAL dtcur

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
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

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

! Variables that are internal to this subroutine.

      INTEGER ixx, iyy, izz, ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1

! These variables go through the interface of the 1d subroutine.

      INTEGER oned_begin, oned_end, oned_off

      REAL dtmove

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1      xb_use, yb_use, yb_use_m1, zb_use,
     1      flatten, gamma_zonecen, molewt_zonecen

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR) :: vmin, vmax

! ----------------------------------------------------------------------

      oned_begin = 1
      oned_end = ix1
      oned_off = ioffx

      dtmove = dtcur

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, xc_use, dx_use, xb_use, yc_use, dy_use,    &
!$OMP&    yb_use, yb_use_m1, zc_use, dz_use, zb_use,                    &
!$OMP&    flatten, w, entropy_1d, gamma_zonecen, molewt_zonecen,        &
!$OMP&    vmin, vmax, corrupted_zone)                                   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    oned_begin, oned_end, oned_off, dtmove, igeom, pi,            &
!$OMP&    rhofloor, prsfloor, smallnum, gamma, gasconst,                &
!$OMP&    u, entropygr,                                                 &
!$OMP&    rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,       &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,            &
!$OMP&    temperaturegr, divvelgr, msonicgr, flattengr,                 &
!$OMP&    xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,         &
!$OMP&    zb_save, dz_save, zc_save)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

      zc_use ( 1 - ioffx : ix1 + ioffx) = zc_save ( izz)
      dz_use ( 1 - ioffx : ix1 + ioffx) = dz_save ( izz)
      zb_use ( 1 - ioffx : ix1 + ioffx) = zb_save ( izz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

        yc_use ( 1 - ioffx : ix1 + ioffx) = yc_save ( iyy)
        dy_use ( 1 - ioffx : ix1 + ioffx) = dy_save ( iyy)
        yb_use ( 1 - ioffx : ix1 + ioffx) = yb_save ( iyy)
        yb_use_m1 ( 1 - ioffx : ix1 + ioffx) = yb_save ( iyy - 1)

! ------------------------------

          xc_use ( 1 - ioffx : ix1 + ioffx) =
     1       xc_save ( 1 - ioffx : ix1 + ioffx)
          dx_use ( 1 - ioffx : ix1 + ioffx) =
     1       dx_save ( 1 - ioffx : ix1 + ioffx)
          xb_use ( 0 - ioffx : ix1 + ioffx) =
     1       xb_save ( 0 - ioffx : ix1 + ioffx)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

! Find the local min and max of the density variable.
          vmin ( ixx, 1) = AMIN1 (
     1  rhogr ( ixxm1, iyy, izz), rhogr ( ixx, iyy, izz),
     1  rhogr ( ixxp1, iyy, izz),
     1  rhogr ( ixxm1, iyyp1, izz), rhogr ( ixx, iyyp1, izz), 
     1  rhogr ( ixxp1, iyyp1, izz),
     1  rhogr ( ixxm1, iyym1, izz), rhogr ( ixx, iyym1, izz), 
     1  rhogr ( ixxp1, iyym1, izz),

     1  rhogr ( ixxm1, iyy, izzp1), rhogr ( ixx, iyy, izzp1),
     1  rhogr ( ixxp1, iyy, izzp1),
     1  rhogr ( ixxm1, iyyp1, izzp1), rhogr ( ixx, iyyp1, izzp1), 
     1  rhogr ( ixxp1, iyyp1, izzp1),
     1  rhogr ( ixxm1, iyym1, izzp1), rhogr ( ixx, iyym1, izzp1), 
     1  rhogr ( ixxp1, iyym1, izzp1),

     1  rhogr ( ixxm1, iyy, izzm1), rhogr ( ixx, iyy, izzm1),
     1  rhogr ( ixxp1, iyy, izzm1),
     1  rhogr ( ixxm1, iyyp1, izzm1), rhogr ( ixx, iyyp1, izzm1), 
     1  rhogr ( ixxp1, iyyp1, izzm1),
     1  rhogr ( ixxm1, iyym1, izzm1), rhogr ( ixx, iyym1, izzm1), 
     1  rhogr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 1) = AMAX1 (
     1  rhogr ( ixxm1, iyy, izz), rhogr ( ixx, iyy, izz),
     1  rhogr ( ixxp1, iyy, izz),
     1  rhogr ( ixxm1, iyyp1, izz), rhogr ( ixx, iyyp1, izz), 
     1  rhogr ( ixxp1, iyyp1, izz),
     1  rhogr ( ixxm1, iyym1, izz), rhogr ( ixx, iyym1, izz), 
     1  rhogr ( ixxp1, iyym1, izz),

     1  rhogr ( ixxm1, iyy, izzp1), rhogr ( ixx, iyy, izzp1),
     1  rhogr ( ixxp1, iyy, izzp1),
     1  rhogr ( ixxm1, iyyp1, izzp1), rhogr ( ixx, iyyp1, izzp1), 
     1  rhogr ( ixxp1, iyyp1, izzp1),
     1  rhogr ( ixxm1, iyym1, izzp1), rhogr ( ixx, iyym1, izzp1), 
     1  rhogr ( ixxp1, iyym1, izzp1),

     1  rhogr ( ixxm1, iyy, izzm1), rhogr ( ixx, iyy, izzm1),
     1  rhogr ( ixxp1, iyy, izzm1),
     1  rhogr ( ixxm1, iyyp1, izzm1), rhogr ( ixx, iyyp1, izzm1), 
     1  rhogr ( ixxp1, iyyp1, izzm1),
     1  rhogr ( ixxm1, iyym1, izzm1), rhogr ( ixx, iyym1, izzm1), 
     1  rhogr ( ixxp1, iyym1, izzm1) )

! ----------

! Find the local min and max of the pressure variable.
          vmin ( ixx, 5) = AMIN1 (
     1  prsgr ( ixxm1, iyy, izz), prsgr ( ixx, iyy, izz),
     1  prsgr ( ixxp1, iyy, izz),
     1  prsgr ( ixxm1, iyyp1, izz), prsgr ( ixx, iyyp1, izz), 
     1  prsgr ( ixxp1, iyyp1, izz),
     1  prsgr ( ixxm1, iyym1, izz), prsgr ( ixx, iyym1, izz), 
     1  prsgr ( ixxp1, iyym1, izz),

     1  prsgr ( ixxm1, iyy, izzp1), prsgr ( ixx, iyy, izzp1),
     1  prsgr ( ixxp1, iyy, izzp1),
     1  prsgr ( ixxm1, iyyp1, izzp1), prsgr ( ixx, iyyp1, izzp1), 
     1  prsgr ( ixxp1, iyyp1, izzp1),
     1  prsgr ( ixxm1, iyym1, izzp1), prsgr ( ixx, iyym1, izzp1), 
     1  prsgr ( ixxp1, iyym1, izzp1),

     1  prsgr ( ixxm1, iyy, izzm1), prsgr ( ixx, iyy, izzm1),
     1  prsgr ( ixxp1, iyy, izzm1),
     1  prsgr ( ixxm1, iyyp1, izzm1), prsgr ( ixx, iyyp1, izzm1), 
     1  prsgr ( ixxp1, iyyp1, izzm1),
     1  prsgr ( ixxm1, iyym1, izzm1), prsgr ( ixx, iyym1, izzm1), 
     1  prsgr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 5) = AMAX1 (
     1  prsgr ( ixxm1, iyy, izz), prsgr ( ixx, iyy, izz),
     1  prsgr ( ixxp1, iyy, izz),
     1  prsgr ( ixxm1, iyyp1, izz), prsgr ( ixx, iyyp1, izz), 
     1  prsgr ( ixxp1, iyyp1, izz),
     1  prsgr ( ixxm1, iyym1, izz), prsgr ( ixx, iyym1, izz), 
     1  prsgr ( ixxp1, iyym1, izz),

     1  prsgr ( ixxm1, iyy, izzp1), prsgr ( ixx, iyy, izzp1),
     1  prsgr ( ixxp1, iyy, izzp1),
     1  prsgr ( ixxm1, iyyp1, izzp1), prsgr ( ixx, iyyp1, izzp1), 
     1  prsgr ( ixxp1, iyyp1, izzp1),
     1  prsgr ( ixxm1, iyym1, izzp1), prsgr ( ixx, iyym1, izzp1), 
     1  prsgr ( ixxp1, iyym1, izzp1),

     1  prsgr ( ixxm1, iyy, izzm1), prsgr ( ixx, iyy, izzm1),
     1  prsgr ( ixxp1, iyy, izzm1),
     1  prsgr ( ixxm1, iyyp1, izzm1), prsgr ( ixx, iyyp1, izzm1), 
     1  prsgr ( ixxp1, iyyp1, izzm1),
     1  prsgr ( ixxm1, iyym1, izzm1), prsgr ( ixx, iyym1, izzm1), 
     1  prsgr ( ixxp1, iyym1, izzm1) )

! ----------
#if ( REDUCE_VEL_BF == 1)

! Find the local min and max of the x-velocity variable.
          vmin ( ixx, 2) = AMIN1 (
     1  vlxgr ( ixxm1, iyy, izz), vlxgr ( ixx, iyy, izz),
     1  vlxgr ( ixxp1, iyy, izz),
     1  vlxgr ( ixxm1, iyyp1, izz), vlxgr ( ixx, iyyp1, izz),
     1  vlxgr ( ixxp1, iyyp1, izz),
     1  vlxgr ( ixxm1, iyym1, izz), vlxgr ( ixx, iyym1, izz),
     1  vlxgr ( ixxp1, iyym1, izz),

     1  vlxgr ( ixxm1, iyy, izzp1), vlxgr ( ixx, iyy, izzp1),
     1  vlxgr ( ixxp1, iyy, izzp1),
     1  vlxgr ( ixxm1, iyyp1, izzp1), vlxgr ( ixx, iyyp1, izzp1),
     1  vlxgr ( ixxp1, iyyp1, izzp1),
     1  vlxgr ( ixxm1, iyym1, izzp1), vlxgr ( ixx, iyym1, izzp1),
     1  vlxgr ( ixxp1, iyym1, izzp1),

     1  vlxgr ( ixxm1, iyy, izzm1), vlxgr ( ixx, iyy, izzm1),
     1  vlxgr ( ixxp1, iyy, izzm1),
     1  vlxgr ( ixxm1, iyyp1, izzm1), vlxgr ( ixx, iyyp1, izzm1),
     1  vlxgr ( ixxp1, iyyp1, izzm1),
     1  vlxgr ( ixxm1, iyym1, izzm1), vlxgr ( ixx, iyym1, izzm1),
     1  vlxgr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 2) = AMAX1 (
     1  vlxgr ( ixxm1, iyy, izz), vlxgr ( ixx, iyy, izz),
     1  vlxgr ( ixxp1, iyy, izz),
     1  vlxgr ( ixxm1, iyyp1, izz), vlxgr ( ixx, iyyp1, izz),
     1  vlxgr ( ixxp1, iyyp1, izz),
     1  vlxgr ( ixxm1, iyym1, izz), vlxgr ( ixx, iyym1, izz),
     1  vlxgr ( ixxp1, iyym1, izz),

     1  vlxgr ( ixxm1, iyy, izzp1), vlxgr ( ixx, iyy, izzp1),
     1  vlxgr ( ixxp1, iyy, izzp1),
     1  vlxgr ( ixxm1, iyyp1, izzp1), vlxgr ( ixx, iyyp1, izzp1),
     1  vlxgr ( ixxp1, iyyp1, izzp1),
     1  vlxgr ( ixxm1, iyym1, izzp1), vlxgr ( ixx, iyym1, izzp1),
     1  vlxgr ( ixxp1, iyym1, izzp1),

     1  vlxgr ( ixxm1, iyy, izzm1), vlxgr ( ixx, iyy, izzm1),
     1  vlxgr ( ixxp1, iyy, izzm1),
     1  vlxgr ( ixxm1, iyyp1, izzm1), vlxgr ( ixx, iyyp1, izzm1),
     1  vlxgr ( ixxp1, iyyp1, izzm1),
     1  vlxgr ( ixxm1, iyym1, izzm1), vlxgr ( ixx, iyym1, izzm1),
     1  vlxgr ( ixxp1, iyym1, izzm1) )

! ----------

! Find the local min and max of the y-velocity variable.
          vmin ( ixx, 3) = AMIN1 (
     1  vlygr ( ixxm1, iyy, izz), vlygr ( ixx, iyy, izz),
     1  vlygr ( ixxp1, iyy, izz),
     1  vlygr ( ixxm1, iyyp1, izz), vlygr ( ixx, iyyp1, izz),
     1  vlygr ( ixxp1, iyyp1, izz),
     1  vlygr ( ixxm1, iyym1, izz), vlygr ( ixx, iyym1, izz),
     1  vlygr ( ixxp1, iyym1, izz),

     1  vlygr ( ixxm1, iyy, izzp1), vlygr ( ixx, iyy, izzp1),
     1  vlygr ( ixxp1, iyy, izzp1),
     1  vlygr ( ixxm1, iyyp1, izzp1), vlygr ( ixx, iyyp1, izzp1),
     1  vlygr ( ixxp1, iyyp1, izzp1),
     1  vlygr ( ixxm1, iyym1, izzp1), vlygr ( ixx, iyym1, izzp1),
     1  vlygr ( ixxp1, iyym1, izzp1),

     1  vlygr ( ixxm1, iyy, izzm1), vlygr ( ixx, iyy, izzm1),
     1  vlygr ( ixxp1, iyy, izzm1),
     1  vlygr ( ixxm1, iyyp1, izzm1), vlygr ( ixx, iyyp1, izzm1),
     1  vlygr ( ixxp1, iyyp1, izzm1),
     1  vlygr ( ixxm1, iyym1, izzm1), vlygr ( ixx, iyym1, izzm1),
     1  vlygr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 3) = AMAX1 (
     1  vlygr ( ixxm1, iyy, izz), vlygr ( ixx, iyy, izz),
     1  vlygr ( ixxp1, iyy, izz),
     1  vlygr ( ixxm1, iyyp1, izz), vlygr ( ixx, iyyp1, izz),
     1  vlygr ( ixxp1, iyyp1, izz),
     1  vlygr ( ixxm1, iyym1, izz), vlygr ( ixx, iyym1, izz),
     1  vlygr ( ixxp1, iyym1, izz),

     1  vlygr ( ixxm1, iyy, izzp1), vlygr ( ixx, iyy, izzp1),
     1  vlygr ( ixxp1, iyy, izzp1),
     1  vlygr ( ixxm1, iyyp1, izzp1), vlygr ( ixx, iyyp1, izzp1),
     1  vlygr ( ixxp1, iyyp1, izzp1),
     1  vlygr ( ixxm1, iyym1, izzp1), vlygr ( ixx, iyym1, izzp1),
     1  vlygr ( ixxp1, iyym1, izzp1),

     1  vlygr ( ixxm1, iyy, izzm1), vlygr ( ixx, iyy, izzm1),
     1  vlygr ( ixxp1, iyy, izzm1),
     1  vlygr ( ixxm1, iyyp1, izzm1), vlygr ( ixx, iyyp1, izzm1),
     1  vlygr ( ixxp1, iyyp1, izzm1),
     1  vlygr ( ixxm1, iyym1, izzm1), vlygr ( ixx, iyym1, izzm1),
     1  vlygr ( ixxp1, iyym1, izzm1) )

! ----------

! Find the local min and max of the z-velocity variable.
          vmin ( ixx, 4) = AMIN1 (
     1  vlzgr ( ixxm1, iyy, izz), vlzgr ( ixx, iyy, izz),
     1  vlzgr ( ixxp1, iyy, izz),
     1  vlzgr ( ixxm1, iyyp1, izz), vlzgr ( ixx, iyyp1, izz),
     1  vlzgr ( ixxp1, iyyp1, izz),
     1  vlzgr ( ixxm1, iyym1, izz), vlzgr ( ixx, iyym1, izz),
     1  vlzgr ( ixxp1, iyym1, izz),

     1  vlzgr ( ixxm1, iyy, izzp1), vlzgr ( ixx, iyy, izzp1),
     1  vlzgr ( ixxp1, iyy, izzp1),
     1  vlzgr ( ixxm1, iyyp1, izzp1), vlzgr ( ixx, iyyp1, izzp1),
     1  vlzgr ( ixxp1, iyyp1, izzp1),
     1  vlzgr ( ixxm1, iyym1, izzp1), vlzgr ( ixx, iyym1, izzp1),
     1  vlzgr ( ixxp1, iyym1, izzp1),

     1  vlzgr ( ixxm1, iyy, izzm1), vlzgr ( ixx, iyy, izzm1),
     1  vlzgr ( ixxp1, iyy, izzm1),
     1  vlzgr ( ixxm1, iyyp1, izzm1), vlzgr ( ixx, iyyp1, izzm1),
     1  vlzgr ( ixxp1, iyyp1, izzm1),
     1  vlzgr ( ixxm1, iyym1, izzm1), vlzgr ( ixx, iyym1, izzm1),
     1  vlzgr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 4) = AMAX1 (
     1  vlzgr ( ixxm1, iyy, izz), vlzgr ( ixx, iyy, izz),
     1  vlzgr ( ixxp1, iyy, izz),
     1  vlzgr ( ixxm1, iyyp1, izz), vlzgr ( ixx, iyyp1, izz),
     1  vlzgr ( ixxp1, iyyp1, izz),
     1  vlzgr ( ixxm1, iyym1, izz), vlzgr ( ixx, iyym1, izz),
     1  vlzgr ( ixxp1, iyym1, izz),

     1  vlzgr ( ixxm1, iyy, izzp1), vlzgr ( ixx, iyy, izzp1),
     1  vlzgr ( ixxp1, iyy, izzp1),
     1  vlzgr ( ixxm1, iyyp1, izzp1), vlzgr ( ixx, iyyp1, izzp1),
     1  vlzgr ( ixxp1, iyyp1, izzp1),
     1  vlzgr ( ixxm1, iyym1, izzp1), vlzgr ( ixx, iyym1, izzp1),
     1  vlzgr ( ixxp1, iyym1, izzp1),

     1  vlzgr ( ixxm1, iyy, izzm1), vlzgr ( ixx, iyy, izzm1),
     1  vlzgr ( ixxp1, iyy, izzm1),
     1  vlzgr ( ixxm1, iyyp1, izzm1), vlzgr ( ixx, iyyp1, izzm1),
     1  vlzgr ( ixxp1, iyyp1, izzm1),
     1  vlzgr ( ixxm1, iyym1, izzm1), vlzgr ( ixx, iyym1, izzm1),
     1  vlzgr ( ixxp1, iyym1, izzm1) )

! ----------

! Find the local min and max of the x-magnetic field variable.
          vmin ( ixx, 6) = AMIN1 (
     1  bfxgr ( ixxm1, iyy, izz), bfxgr ( ixx, iyy, izz),
     1  bfxgr ( ixxp1, iyy, izz),
     1  bfxgr ( ixxm1, iyyp1, izz), bfxgr ( ixx, iyyp1, izz),
     1  bfxgr ( ixxp1, iyyp1, izz),
     1  bfxgr ( ixxm1, iyym1, izz), bfxgr ( ixx, iyym1, izz),
     1  bfxgr ( ixxp1, iyym1, izz),

     1  bfxgr ( ixxm1, iyy, izzp1), bfxgr ( ixx, iyy, izzp1),
     1  bfxgr ( ixxp1, iyy, izzp1),
     1  bfxgr ( ixxm1, iyyp1, izzp1), bfxgr ( ixx, iyyp1, izzp1),
     1  bfxgr ( ixxp1, iyyp1, izzp1),
     1  bfxgr ( ixxm1, iyym1, izzp1), bfxgr ( ixx, iyym1, izzp1),
     1  bfxgr ( ixxp1, iyym1, izzp1),

     1  bfxgr ( ixxm1, iyy, izzm1), bfxgr ( ixx, iyy, izzm1),
     1  bfxgr ( ixxp1, iyy, izzm1),
     1  bfxgr ( ixxm1, iyyp1, izzm1), bfxgr ( ixx, iyyp1, izzm1),
     1  bfxgr ( ixxp1, iyyp1, izzm1),
     1  bfxgr ( ixxm1, iyym1, izzm1), bfxgr ( ixx, iyym1, izzm1),
     1  bfxgr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 6) = AMAX1 (
     1  bfxgr ( ixxm1, iyy, izz), bfxgr ( ixx, iyy, izz),
     1  bfxgr ( ixxp1, iyy, izz),
     1  bfxgr ( ixxm1, iyyp1, izz), bfxgr ( ixx, iyyp1, izz),
     1  bfxgr ( ixxp1, iyyp1, izz),
     1  bfxgr ( ixxm1, iyym1, izz), bfxgr ( ixx, iyym1, izz),
     1  bfxgr ( ixxp1, iyym1, izz),

     1  bfxgr ( ixxm1, iyy, izzp1), bfxgr ( ixx, iyy, izzp1),
     1  bfxgr ( ixxp1, iyy, izzp1),
     1  bfxgr ( ixxm1, iyyp1, izzp1), bfxgr ( ixx, iyyp1, izzp1),
     1  bfxgr ( ixxp1, iyyp1, izzp1),
     1  bfxgr ( ixxm1, iyym1, izzp1), bfxgr ( ixx, iyym1, izzp1),
     1  bfxgr ( ixxp1, iyym1, izzp1),

     1  bfxgr ( ixxm1, iyy, izzm1), bfxgr ( ixx, iyy, izzm1),
     1  bfxgr ( ixxp1, iyy, izzm1),
     1  bfxgr ( ixxm1, iyyp1, izzm1), bfxgr ( ixx, iyyp1, izzm1),
     1  bfxgr ( ixxp1, iyyp1, izzm1),
     1  bfxgr ( ixxm1, iyym1, izzm1), bfxgr ( ixx, iyym1, izzm1),
     1  bfxgr ( ixxp1, iyym1, izzm1) )

! ----------

! Find the local min and max of the y-magnetic field variable.
          vmin ( ixx, 7) = AMIN1 (
     1  bfygr ( ixxm1, iyy, izz), bfygr ( ixx, iyy, izz),
     1  bfygr ( ixxp1, iyy, izz),
     1  bfygr ( ixxm1, iyyp1, izz), bfygr ( ixx, iyyp1, izz),
     1  bfygr ( ixxp1, iyyp1, izz),
     1  bfygr ( ixxm1, iyym1, izz), bfygr ( ixx, iyym1, izz),
     1  bfygr ( ixxp1, iyym1, izz),

     1  bfygr ( ixxm1, iyy, izzp1), bfygr ( ixx, iyy, izzp1),
     1  bfygr ( ixxp1, iyy, izzp1),
     1  bfygr ( ixxm1, iyyp1, izzp1), bfygr ( ixx, iyyp1, izzp1),
     1  bfygr ( ixxp1, iyyp1, izzp1),
     1  bfygr ( ixxm1, iyym1, izzp1), bfygr ( ixx, iyym1, izzp1),
     1  bfygr ( ixxp1, iyym1, izzp1),

     1  bfygr ( ixxm1, iyy, izzm1), bfygr ( ixx, iyy, izzm1),
     1  bfygr ( ixxp1, iyy, izzm1),
     1  bfygr ( ixxm1, iyyp1, izzm1), bfygr ( ixx, iyyp1, izzm1),
     1  bfygr ( ixxp1, iyyp1, izzm1),
     1  bfygr ( ixxm1, iyym1, izzm1), bfygr ( ixx, iyym1, izzm1),
     1  bfygr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 7) = AMAX1 (
     1  bfygr ( ixxm1, iyy, izz), bfygr ( ixx, iyy, izz),
     1  bfygr ( ixxp1, iyy, izz),
     1  bfygr ( ixxm1, iyyp1, izz), bfygr ( ixx, iyyp1, izz),
     1  bfygr ( ixxp1, iyyp1, izz),
     1  bfygr ( ixxm1, iyym1, izz), bfygr ( ixx, iyym1, izz),
     1  bfygr ( ixxp1, iyym1, izz),

     1  bfygr ( ixxm1, iyy, izzp1), bfygr ( ixx, iyy, izzp1),
     1  bfygr ( ixxp1, iyy, izzp1),
     1  bfygr ( ixxm1, iyyp1, izzp1), bfygr ( ixx, iyyp1, izzp1),
     1  bfygr ( ixxp1, iyyp1, izzp1),
     1  bfygr ( ixxm1, iyym1, izzp1), bfygr ( ixx, iyym1, izzp1),
     1  bfygr ( ixxp1, iyym1, izzp1),

     1  bfygr ( ixxm1, iyy, izzm1), bfygr ( ixx, iyy, izzm1),
     1  bfygr ( ixxp1, iyy, izzm1),
     1  bfygr ( ixxm1, iyyp1, izzm1), bfygr ( ixx, iyyp1, izzm1),
     1  bfygr ( ixxp1, iyyp1, izzm1),
     1  bfygr ( ixxm1, iyym1, izzm1), bfygr ( ixx, iyym1, izzm1),
     1  bfygr ( ixxp1, iyym1, izzm1) )

! ----------

! Find the local min and max of the z-magnetic field variable.
          vmin ( ixx, 8) = AMIN1 (
     1  bfzgr ( ixxm1, iyy, izz), bfzgr ( ixx, iyy, izz),
     1  bfzgr ( ixxp1, iyy, izz),
     1  bfzgr ( ixxm1, iyyp1, izz), bfzgr ( ixx, iyyp1, izz),
     1  bfzgr ( ixxp1, iyyp1, izz),
     1  bfzgr ( ixxm1, iyym1, izz), bfzgr ( ixx, iyym1, izz),
     1  bfzgr ( ixxp1, iyym1, izz),

     1  bfzgr ( ixxm1, iyy, izzp1), bfzgr ( ixx, iyy, izzp1),
     1  bfzgr ( ixxp1, iyy, izzp1),
     1  bfzgr ( ixxm1, iyyp1, izzp1), bfzgr ( ixx, iyyp1, izzp1),
     1  bfzgr ( ixxp1, iyyp1, izzp1),
     1  bfzgr ( ixxm1, iyym1, izzp1), bfzgr ( ixx, iyym1, izzp1),
     1  bfzgr ( ixxp1, iyym1, izzp1),

     1  bfzgr ( ixxm1, iyy, izzm1), bfzgr ( ixx, iyy, izzm1),
     1  bfzgr ( ixxp1, iyy, izzm1),
     1  bfzgr ( ixxm1, iyyp1, izzm1), bfzgr ( ixx, iyyp1, izzm1),
     1  bfzgr ( ixxp1, iyyp1, izzm1),
     1  bfzgr ( ixxm1, iyym1, izzm1), bfzgr ( ixx, iyym1, izzm1),
     1  bfzgr ( ixxp1, iyym1, izzm1) )

          vmax ( ixx, 8) = AMAX1 (
     1  bfzgr ( ixxm1, iyy, izz), bfzgr ( ixx, iyy, izz),
     1  bfzgr ( ixxp1, iyy, izz),
     1  bfzgr ( ixxm1, iyyp1, izz), bfzgr ( ixx, iyyp1, izz),
     1  bfzgr ( ixxp1, iyyp1, izz),
     1  bfzgr ( ixxm1, iyym1, izz), bfzgr ( ixx, iyym1, izz),
     1  bfzgr ( ixxp1, iyym1, izz),

     1  bfzgr ( ixxm1, iyy, izzp1), bfzgr ( ixx, iyy, izzp1),
     1  bfzgr ( ixxp1, iyy, izzp1),
     1  bfzgr ( ixxm1, iyyp1, izzp1), bfzgr ( ixx, iyyp1, izzp1),
     1  bfzgr ( ixxp1, iyyp1, izzp1),
     1  bfzgr ( ixxm1, iyym1, izzp1), bfzgr ( ixx, iyym1, izzp1),
     1  bfzgr ( ixxp1, iyym1, izzp1),

     1  bfzgr ( ixxm1, iyy, izzm1), bfzgr ( ixx, iyy, izzm1),
     1  bfzgr ( ixxp1, iyy, izzm1),
     1  bfzgr ( ixxm1, iyyp1, izzm1), bfzgr ( ixx, iyyp1, izzm1),
     1  bfzgr ( ixxp1, iyyp1, izzm1),
     1  bfzgr ( ixxm1, iyym1, izzm1), bfzgr ( ixx, iyym1, izzm1),
     1  bfzgr ( ixxp1, iyym1, izzm1) )

#endif
! End of "REDUCE_VEL_BF == 1"
! ----------

          gamma_zonecen ( ixx) = gamma_eos_gr ( ixx, iyy, izz)
          molewt_zonecen ( ixx) = mean_mol_wt_gr ( ixx, iyy, izz)

          flatten ( ixx) = flattengr ( ixx, iyy, izz)

          w ( ixx, 1: NFIELD + 1 + NFLUID, 1: NUM_ZONECEN_BASES) =
     1    u ( ixx, iyy, izz, 1: NFIELD + 1 + NFLUID,
     1                                     1: NUM_ZONECEN_BASES)

          entropy_1d ( ixx, 1:4) = entropygr ( ixx, iyy, izz, 1:4)

! ----------

          END DO ! End of "DO ixx =".

! ------------------------------

! Make call to "RECONST_POS_RHO_PRS_MHD_CUBIC".

          CALL RECONST_POS_RHO_PRS_MHD_CUBIC

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove, igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gasconst,

     1             xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1             xb_use, yb_use, yb_use_m1, zb_use,

     1             flatten, w, entropy_1d,
     1             gamma_zonecen, molewt_zonecen,
     1             vmin, vmax, corrupted_zone)

! ------------------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

          u ( ixx, iyy, izz, 1: NFIELD + 1 + NFLUID,
     1                                     1: NUM_ZONECEN_BASES) =
     1    w ( ixx, 1: NFIELD + 1 + NFLUID, 1: NUM_ZONECEN_BASES)

          entropygr ( ixx, iyy, izz, 5) = entropy_1d ( ixx, 5)

          END DO ! End of "DO ixx =".

! ------------------------------

        END DO ! End of "DO iyy =".

      END DO ! End of "DO izz =".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE POSITIVE_RHO_PRS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE RECONST_POS_RHO_PRS_MHD_CUBIC

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove, igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gasconst,

     1             xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1             xb_use, yb_use, yb_use_m1, zb_use,

     1             flatten, w, entropy_1d,
     1             gamma_zonecen, molewt_zonecen,
     1             vmin, vmax, corrupted_zone)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Checks the reconstruction for density and pressure positivity. Without
! these the fluxes become non-convex and unphysical. Comes in with
! the following 1d arrays loaded:

! xc_use, yc_use, zc_use, dx_use, dy_use, dz_use : geometric 1d variables.
! flatten : Has the value from "flattengr".
! w : modal vector of conserved variables.
! entropy_1d : Comes in with the four spatial modes and exits with one time
! mode of a second order accurate entropy function.
! gamma_zonecen, molewt_zonecen : zone-centered gamma and molewt.
! vmin, vmax : Min and max values of the primitive variables from all the
! neighboring zones. Done outside.
! corrupted_zone : If the zone has taken any density or pressure positivity
! fix-ups, it is set to 1; otherwise it is 0.

! The fundamental idea is that we construct nodal values in space (not
! time). We then test them for positivity. If the positivity is not obtained,
! we reduce the spatially varying part of the reconstructed profile.
! This is done in two stages. First, the density variation is shrunk,
! if needed. Second, all the variations (except certain necessary components
! of the magnetic field) are shrunk consistent with obtaining
! positive pressures, if that is needed.

! In this version, a cubic is produced for the MHD equations and solved
! via iteration. It takes more time than the quadratic version. However, it
! also allows much (though not all) of the variation in the magnetic field
! modes to be reduced. Note though that the cubic is only invoked in those
! zones that have a negative pressure. Such zones are usually nonexistent
! or minimal in most real simulations.

! This is one of ONLY two places where we may perhaps trade
! conservation in return for physical consistency. The other place is
! in the "EVALUATE_TIMESTEP" routine.

! ----------------------------------------------------------------------

      IMPLICIT NONE

! These variables go through the subroutine interface.

      INTEGER oned_begin, oned_end, oned_off, igeom

      REAL dtmove, pi, rhofloor, prsfloor, smallnum, gamma,
     1     gasconst

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1      xb_use, yb_use, yb_use_m1, zb_use,
     1      flatten, gamma_zonecen, molewt_zonecen

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR) :: vmin, vmax

! ----------------------------------------------------------------------

! Variables that are internal to the subroutine.

      INTEGER ix, iquadpt, imode, ivar, num_pp_quadpts, iter, itermax

      REAL fourpi, ob4pi, ob8pi, ob_root_4pi, expandfac,
     1     del_rho, del_eng, del_mx, del_my, del_mz,
     1     bx_c, by_c, bz_c, bx_v, by_v, bz_v,
     1     del_bx, del_by, del_bz, m_bar2, del_m2, m_dot_del_m,
     1     b_bar2, del_b2, b_dot_del_b,
     1     ax1, ax2, ay1, ay2, az1, az2, ob_dx, ob_dy, ob_dz,
     1     tempaa, tempab, tempac, tempad

! Since we are interested in extremal variations within a zone, we need
! to define nodal sets that have more quadrature points than the smallest
! symmetrical set of quadrature points.
#if ( NUMORDER == 2)
      PARAMETER ( num_pp_quadpts = 8)
#elif ( NUMORDER == 3)
      PARAMETER ( num_pp_quadpts = 21)
#elif ( NUMORDER == 4)
      PARAMETER ( num_pp_quadpts = 33)
#endif

! These are the flags for flagging various situations and coefficients
! for curing various pathologies.

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: pp_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: reduce_modes

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                     num_pp_quadpts) :: pp_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                     num_pp_quadpts) :: reduce_quadpt,
     1                     a0, a1, a2, a3, flip_sign,
     1                     x_low, x_high, f_low, f_high,
     1                     x_root, f_root, df_root

! These are the bounds on the primitive variables from neighboring zones
! with some inbuilt expansions.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR) :: vminuse, vmaxuse, tempvec

! These are the flow variables at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  num_pp_quadpts) ::

     1      rho, prs, eng_t, ent, vlx, vly, vlz, bfx, bfy, bfz,
     1      bfx_const, bfy_const, bfz_const, bfx_var, bfy_var, bfz_var,
     1      gamma_eos_quadpt, molewt_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0: NFLUID, num_pp_quadpts) :: specie

! These are the conserved variables at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, num_pp_quadpts) ::

     1      u_quadpt, tu

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob4pi = 1.0 / fourpi
      ob8pi = 0.5 * ob4pi
      ob_root_4pi = 1.0 / SQRT ( fourpi)

! Maximum number of iterations for the cubic root solver. Ok to do just a few.
      itermax = 4

! Set the factor by which we want to expand the range of each variable.
! Setting this to 0 gives us the TVD range in multi-dimensions.
      expandfac = 0.4
! ----------------------------------------------------------------------

! STEP 0) Start by saying that none of the zones are corrupted, none have
! negative pressures and that none of them need to have their moments reduced.
! Make the same assertion at all the quadrature points within each zone.
! By the end of this subroutine, most zones should continue to be that way.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      corrupted_zone ( ix) = 0
      pp_zone ( ix) = 0
      reduce_modes ( ix) = 1.0

      pp_quadpt ( ix, :) = 0
      reduce_quadpt ( ix, :) = 1.0

      gamma_eos_quadpt ( ix, :) = gamma_zonecen ( ix)
      molewt_quadpt ( ix, :) = molewt_zonecen ( ix)

      END DO

! ----------------------------------------------------------------------

! STEP 1) Set "corrupted_zone" to 0 and expand "vminuse, vmaxuse" if possible.
! The bounds on variables that should remain positive are set
! multiplicatively. The bounds on variables that aren't positive definite
! are set based on expanding the range.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      corrupted_zone ( ix) = 0.0

! Expand density range:
      vminuse ( ix, 1) = vmin ( ix, 1)
     1      * ( 1.0 - expandfac + expandfac * flatten ( ix) )
      vminuse ( ix, 1) = AMAX1 ( vminuse ( ix, 1), rhofloor)
      vmaxuse ( ix, 1) = vmax ( ix, 1)
     1      * ( 1.0 + expandfac - expandfac * flatten ( ix) )

! Expand pressure range:
      vminuse ( ix, 5) = vmin ( ix, 5)
     1      * ( 1.0 - expandfac + expandfac * flatten ( ix) )
      vminuse ( ix, 5) = AMAX1 ( vminuse ( ix, 5), prsfloor)
      vmaxuse ( ix, 5) = vmax ( ix, 5)
     1      * ( 1.0 + expandfac - expandfac * flatten ( ix) )

#if ( REDUCE_VEL_BF == 1)
! Expand velocity range:
      tempvec ( ix, 2:4) = expandfac * ( 1.0 - flatten ( ix) )
     1       * ( vmax ( ix, 2:4) - vmin ( ix, 2:4) )
      vminuse ( ix, 2:4) = vmin ( ix, 2:4) - tempvec ( ix, 2:4)
      vmaxuse ( ix, 2:4) = vmax ( ix, 2:4) + tempvec ( ix, 2:4)

! Expand magnetic field range:
      tempvec ( ix, 6:8) = expandfac * ( 1.0 - flatten ( ix) )
     1       * ( vmax ( ix, 6:8) - vmin ( ix, 6:8) )
      vminuse ( ix, 6:8) = vmin ( ix, 6:8) - tempvec ( ix, 6:8)
      vmaxuse ( ix, 6:8) = vmax ( ix, 6:8) + tempvec ( ix, 6:8)
#endif

#if ( NFLUID >= 1)
! Expand the specie fraction range:
      vminuse ( ix, 9:8 + NFLUID) = vmin ( ix, 9:8 + NFLUID)
     1      * ( 1.0 - expandfac + expandfac * flatten ( ix) )
      vminuse ( ix, 9:8 + NFLUID) =
     1                      AMAX1 ( vminuse ( ix, 9:8 + NFLUID), 0.0)
      vmaxuse ( ix, 9:8 + NFLUID) = vmax ( ix, 9:8 + NFLUID)
     1      * ( 1.0 + expandfac - expandfac * flatten ( ix) )
      vmaxuse ( ix, 9:8 + NFLUID) =
     1                      AMIN1 ( vmaxuse ( ix, 9:8 + NFLUID), 1.0)
#endif

      END DO ! End of "DO ix =".

! ----------------------------------------------------------------------

! STEP 2) Use only the spatial modes to obtain the conserved variables
! at the nodes of interest. I.e., use the spatial modes to make 
! "u_quadpt ( ix, :, 1: num_pp_quadpts)"
! from "w ( ix, :, 1: NUM_ZONECEN_BASES)" at time t = 0.

! Certain magnetic field modes cannot be varied without losing consistency
! of the normal component of the magnetic field at the faces. Their
! contributions are built and stored in "bfx_const" etc.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal ( at zone vertices).

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 02) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 03) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 04) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 05) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 06) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 07) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 08) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

! ----------

      bfx_const ( ix, 01) = w ( ix, 6, 01) + 0.5 * w ( ix, 6, 02)

      bfx_const ( ix, 02) = w ( ix, 6, 01) - 0.5 * w ( ix, 6, 02)

      bfx_const ( ix, 03) = bfx_const ( ix, 01)

      bfx_const ( ix, 04) = bfx_const ( ix, 02)

      bfx_const ( ix, 05) = bfx_const ( ix, 01)

      bfx_const ( ix, 06) = bfx_const ( ix, 02)

      bfx_const ( ix, 07) = bfx_const ( ix, 01)

      bfx_const ( ix, 08) = bfx_const ( ix, 02)

! ----------

      bfy_const ( ix, 01) = w ( ix, 7, 01) + 0.5 * w ( ix, 7, 03)

      bfy_const ( ix, 02) = bfy_const ( ix, 01)

      bfy_const ( ix, 03) = w ( ix, 7, 01) - 0.5 * w ( ix, 7, 03)

      bfy_const ( ix, 04) = bfy_const ( ix, 03)

      bfy_const ( ix, 05) = bfy_const ( ix, 01)

      bfy_const ( ix, 06) = bfy_const ( ix, 01)

      bfy_const ( ix, 07) = bfy_const ( ix, 03)

      bfy_const ( ix, 08) = bfy_const ( ix, 03)

! ----------

      bfz_const ( ix, 01) = w ( ix, 8, 01) + 0.5 * w ( ix, 8, 04)

      bfz_const ( ix, 02) = bfz_const ( ix, 01)

      bfz_const ( ix, 03) = bfz_const ( ix, 01)

      bfz_const ( ix, 04) = bfz_const ( ix, 01)

      bfz_const ( ix, 05) = w ( ix, 8, 01) - 0.5 * w ( ix, 8, 04)

      bfz_const ( ix, 06) = bfz_const ( ix, 05)

      bfz_const ( ix, 07) = bfz_const ( ix, 05)

      bfz_const ( ix, 08) = bfz_const ( ix, 05)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal ( including zone vertices).

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = w ( ix, :, 01) - ( w ( ix, :, 05)
     1 + w ( ix, :, 06) + w ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( w ( ix, :, 02) + w ( ix, :, 03) )
     1 + 0.25 * ( w ( ix, :, 05) + w ( ix, :, 06) + w ( ix, :, 08) )

      u_quadpt ( ix, :, 03) = u_quadpt ( ix, :, 02) - w ( ix, :, 03)
     1 - 0.5 * w ( ix, :, 08)

      u_quadpt ( ix, :, 04) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( w ( ix, :, 02) + w ( ix, :, 04) )
     1 + 0.25 * ( w ( ix, :, 05) + w ( ix, :, 07) + w ( ix, :, 10) )

      u_quadpt ( ix, :, 05) = u_quadpt ( ix, :, 04) - w ( ix, :, 04)
     1 - 0.5 * w ( ix, :, 10)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 03) - w ( ix, :, 02)
     1 + w ( ix, :, 03)

      u_quadpt ( ix, :, 07) = u_quadpt ( ix, :, 02) - w ( ix, :, 02)
     1 - w ( ix, :, 03)

      u_quadpt ( ix, :, 08) = u_quadpt ( ix, :, 05) - w ( ix, :, 02)
     1 + w ( ix, :, 04)

      u_quadpt ( ix, :, 09) = u_quadpt ( ix, :, 04) - w ( ix, :, 02)
     1 - w ( ix, :, 04)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( w ( ix, :, 03) + w ( ix, :, 04) )
     1 + 0.25 * ( w ( ix, :, 06) + w ( ix, :, 07) + w ( ix, :, 09) )

      u_quadpt ( ix, :, 11) = u_quadpt ( ix, :, 10) - w ( ix, :, 03)
     1 - 0.5 * w ( ix, :, 09)

      u_quadpt ( ix, :, 12) = u_quadpt ( ix, :, 11) + w ( ix, :, 03)
     1 - w ( ix, :, 04)

      u_quadpt ( ix, :, 13) = u_quadpt ( ix, :, 10) - w ( ix, :, 03)
     1 - w ( ix, :, 04)

! -----
! Additional nodes at vertices.

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 04) + u_quadpt ( ix, :, 10)
     1 - 2.0 * u_quadpt ( ix, :, 01)
     1 - 0.5 * ( w ( ix, :, 02) + w ( ix, :, 03) + w ( ix, :, 04) )
     1 - 0.25 * ( w ( ix, :, 05) + w ( ix, :, 06) + w ( ix, :, 07) )

      u_quadpt ( ix, :, 15) = 2.0 * u_quadpt ( ix, :, 10)
     1 - u_quadpt ( ix, :, 14) + 0.5 * w ( ix, :, 05)

      u_quadpt ( ix, :, 16) = 2.0 * u_quadpt ( ix, :, 04)
     1 - u_quadpt ( ix, :, 14) + 0.5 * w ( ix, :, 06)

      u_quadpt ( ix, :, 17) = 2.0 * u_quadpt ( ix, :, 11)
     1 - u_quadpt ( ix, :, 16) + 0.5 * w ( ix, :, 05)

      u_quadpt ( ix, :, 18) = 2.0 * u_quadpt ( ix, :, 02)
     1 - u_quadpt ( ix, :, 14) + 0.5 * w ( ix, :, 07)

      u_quadpt ( ix, :, 19) = 2.0 * u_quadpt ( ix, :, 12)
     1 - u_quadpt ( ix, :, 18) + 0.5 * w ( ix, :, 05)

      u_quadpt ( ix, :, 20) = 2.0 * u_quadpt ( ix, :, 05)
     1 - u_quadpt ( ix, :, 18) + 0.5 * w ( ix, :, 06)

      u_quadpt ( ix, :, 21) = 2.0 * u_quadpt ( ix, :, 13)
     1 - u_quadpt ( ix, :, 20) + 0.5 * w ( ix, :, 05)

! ----------

      bfx_const ( ix, 01) = w ( ix, 6, 01) - w ( ix, 6, 05) / 12.0

      bfx_const ( ix, 02) = w ( ix, 6, 01) + 0.5 * w ( ix, 6, 02)
     1 + w ( ix, 6, 05) / 6.0

      bfx_const ( ix, 03) = bfx_const ( ix, 02)

      bfx_const ( ix, 04) = bfx_const ( ix, 02)

      bfx_const ( ix, 05) = bfx_const ( ix, 02)

      bfx_const ( ix, 06) = w ( ix, 6, 01) - 0.5 * w ( ix, 6, 02)
     1 + w ( ix, 6, 05) / 6.0

      bfx_const ( ix, 07) = bfx_const ( ix, 06)

      bfx_const ( ix, 08) = bfx_const ( ix, 06)

      bfx_const ( ix, 09) = bfx_const ( ix, 06)

      bfx_const ( ix, 10) = bfx_const ( ix, 01)

      bfx_const ( ix, 11) = bfx_const ( ix, 01)

      bfx_const ( ix, 12) = bfx_const ( ix, 01)

      bfx_const ( ix, 13) = bfx_const ( ix, 01)


      bfx_const ( ix, 14) = bfx_const ( ix, 02)

      bfx_const ( ix, 15) = bfx_const ( ix, 06)

      bfx_const ( ix, 16) = bfx_const ( ix, 02)

      bfx_const ( ix, 17) = bfx_const ( ix, 06)

      bfx_const ( ix, 18) = bfx_const ( ix, 02)

      bfx_const ( ix, 19) = bfx_const ( ix, 06)

      bfx_const ( ix, 20) = bfx_const ( ix, 02)

      bfx_const ( ix, 21) = bfx_const ( ix, 06)

! ----------

      bfy_const ( ix, 01) = w ( ix, 7, 01) - w ( ix, 7, 06) / 12.0

      bfy_const ( ix, 02) = w ( ix, 7, 01) + 0.5 * w ( ix, 7, 03)
     1 + w ( ix, 7, 06) / 6.0

      bfy_const ( ix, 03) = w ( ix, 7, 01) - 0.5 * w ( ix, 7, 03)
     1 + w ( ix, 7, 06) / 6.0

      bfy_const ( ix, 04) = bfy_const ( ix, 01)

      bfy_const ( ix, 05) = bfy_const ( ix, 01)

      bfy_const ( ix, 06) = bfy_const ( ix, 02)

      bfy_const ( ix, 07) = bfy_const ( ix, 03)

      bfy_const ( ix, 08) = bfy_const ( ix, 01)

      bfy_const ( ix, 09) = bfy_const ( ix, 01)

      bfy_const ( ix, 10) = bfy_const ( ix, 02)

      bfy_const ( ix, 11) = bfy_const ( ix, 03)

      bfy_const ( ix, 12) = bfy_const ( ix, 02)

      bfy_const ( ix, 13) = bfy_const ( ix, 03)


      bfy_const ( ix, 14) = bfy_const ( ix, 02)

      bfy_const ( ix, 15) = bfy_const ( ix, 02)

      bfy_const ( ix, 16) = bfy_const ( ix, 03)

      bfy_const ( ix, 17) = bfy_const ( ix, 03)

      bfy_const ( ix, 18) = bfy_const ( ix, 02)

      bfy_const ( ix, 19) = bfy_const ( ix, 02)

      bfy_const ( ix, 20) = bfy_const ( ix, 03)

      bfy_const ( ix, 21) = bfy_const ( ix, 03)

! ----------

      bfz_const ( ix, 01) = w ( ix, 8, 01) - w ( ix, 8, 07) / 12.0

      bfz_const ( ix, 02) = bfz_const ( ix, 01)

      bfz_const ( ix, 03) = bfz_const ( ix, 01)

      bfz_const ( ix, 04) = w ( ix, 8, 01) + 0.5 * w ( ix, 8, 04)
     1 + w ( ix, 8, 07) / 6.0

      bfz_const ( ix, 05) = w ( ix, 8, 01) - 0.5 * w ( ix, 8, 04)
     1 + w ( ix, 8, 07) / 6.0

      bfz_const ( ix, 06) = bfz_const ( ix, 01)

      bfz_const ( ix, 07) = bfz_const ( ix, 01)

      bfz_const ( ix, 08) = bfz_const ( ix, 04)

      bfz_const ( ix, 09) = bfz_const ( ix, 05)

      bfz_const ( ix, 10) = bfz_const ( ix, 04)

      bfz_const ( ix, 11) = bfz_const ( ix, 04)

      bfz_const ( ix, 12) = bfz_const ( ix, 05)

      bfz_const ( ix, 13) = bfz_const ( ix, 05)


      bfz_const ( ix, 14) = bfz_const ( ix, 04)

      bfz_const ( ix, 15) = bfz_const ( ix, 04)

      bfz_const ( ix, 16) = bfz_const ( ix, 04)

      bfz_const ( ix, 17) = bfz_const ( ix, 04)

      bfz_const ( ix, 18) = bfz_const ( ix, 05)

      bfz_const ( ix, 19) = bfz_const ( ix, 05)

      bfz_const ( ix, 20) = bfz_const ( ix, 05)

      bfz_const ( ix, 21) = bfz_const ( ix, 05)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal (zone vertices are already included).

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = w ( ix, :, 01) - ( w ( ix, :, 05)
     1 + w ( ix, :, 06) + w ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * w ( ix, :, 02) + 30.0 * w ( ix, :, 05)
     1 + 6.0 * w ( ix, :, 11) - 5.0 * w ( ix, :, 16)
     1 - 5.0 * w ( ix, :, 18) ) / 120.0

      u_quadpt ( ix, :, 03) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) ) - 0.015625 * ( 4 * w ( ix, :, 05)
     1 + 3.0 * w ( ix, :, 11) )

      u_quadpt ( ix, :, 04) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) - 3.0 * u_quadpt ( ix, :, 03)
     1 - 0.09375 * w ( ix, :, 11)

      u_quadpt ( ix, :, 05) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 02) - 8.0 * u_quadpt ( ix, :, 03)
     1 - 0.375 * w ( ix, :, 11)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * w ( ix, :, 03) + 30.0 * w ( ix, :, 06)
     1 + 6.0 * w ( ix, :, 12) - 5.0 * w ( ix, :, 14)
     1 - 5.0 * w ( ix, :, 19) ) / 120.0

      u_quadpt ( ix, :, 07) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) ) - 0.015625 * ( 4.0 * w ( ix, :, 06)
     1 + 3.0 * w ( ix, :, 12) )

      u_quadpt ( ix, :, 08) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) - 3.0 * u_quadpt ( ix, :, 07)
     1 - 0.09375 * w ( ix, :, 12)

      u_quadpt ( ix, :, 09) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 06) - 8.0 * u_quadpt ( ix, :, 07)
     1 - 0.375 * w ( ix, :, 12)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * w ( ix, :, 04) + 30.0 * w ( ix, :, 07)
     1 + 6.0 * w ( ix, :, 13) - 5.0 * w ( ix, :, 15)
     1 - 5.0 * w ( ix, :, 17) ) / 120.0

      u_quadpt ( ix, :, 11) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) ) - 0.015625 * ( 4.0 * w ( ix, :, 07)
     1 + 3.0 * w ( ix, :, 13) )

      u_quadpt ( ix, :, 12) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) - 3.0 * u_quadpt ( ix, :, 11)
     1 - 0.09375 * w ( ix, :, 13)

      u_quadpt ( ix, :, 13) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 10) - 8.0 * u_quadpt ( ix, :, 11)
     1 - 0.375 * w ( ix, :, 13)

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 07) + u_quadpt ( ix, :, 10)
     1 - 2.0 * u_quadpt ( ix, :, 01) + 0.25 * ( w ( ix, :, 03)
     1 + w ( ix, :, 08) + w ( ix, :, 09) + w ( ix, :, 10) )
     1 + 0.125 * ( w ( ix, :, 15) + w ( ix, :, 16) + w ( ix, :, 17)
     1 + w ( ix, :, 18) + w ( ix, :, 20) ) + 0.1875 * w ( ix, :, 06)
     1 + 5.0 * ( w ( ix, :, 14) + w ( ix, :, 19) ) / 48.0
     1 + 0.071875 * w ( ix, :, 12)

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 - 0.5 * ( w ( ix, :, 08) + w ( ix, :, 09) )
     1 - 0.25 * ( w ( ix, :, 14) + w ( ix, :, 19) + w ( ix, :, 20) )

      u_quadpt ( ix, :, 16) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( w ( ix, :, 09) + w ( ix, :, 10) )
     1 - 0.25 * ( w ( ix, :, 15) + w ( ix, :, 17) + w ( ix, :, 20) )

      u_quadpt ( ix, :, 17) = u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 16) - u_quadpt ( ix, :, 14) + w ( ix, :, 09)
     1 + 0.5 * w ( ix, :, 20)

      u_quadpt ( ix, :, 18) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * ( w ( ix, :, 08) + w ( ix, :, 10) )
     1 - 0.25 * ( w ( ix, :, 16) + w ( ix, :, 18) + w ( ix, :, 20) )

      u_quadpt ( ix, :, 19) = u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 + 0.5 * ( w ( ix, :, 08) - w ( ix, :, 09) )
     1 - 0.25 * ( w ( ix, :, 14) + w ( ix, :, 19) - w ( ix, :, 20) )

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( w ( ix, :, 09) - w ( ix, :, 10) )
     1 - 0.25 * ( w ( ix, :, 15) + w ( ix, :, 17) - w ( ix, :, 20) )

      u_quadpt ( ix, :, 21) = u_quadpt ( ix, :, 19)
     1 + u_quadpt ( ix, :, 20) - u_quadpt ( ix, :, 18) + w ( ix, :, 09)
     1 - 0.5 * w ( ix, :, 20)

! -----
! Additional nodes at zone edges.

      u_quadpt ( ix, :, 22) = 0.5 * ( u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 16) ) - 0.25 * w ( ix, :, 07)
     1 - 0.125 * ( w ( ix, :, 18) + w ( ix, :, 19) )

      u_quadpt ( ix, :, 23) = 0.5 * ( u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 17) ) - 0.25 * w ( ix, :, 07)
     1 - 0.125 * ( w ( ix, :, 18) - w ( ix, :, 19) )

      u_quadpt ( ix, :, 24) = 0.5 * ( u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 15) ) - 0.25 * w ( ix, :, 06)
     1 - 0.125 * ( w ( ix, :, 16) + w ( ix, :, 17) )

      u_quadpt ( ix, :, 25) = 0.5 * ( u_quadpt ( ix, :, 16)
     1 + u_quadpt ( ix, :, 17) ) - 0.25 * w ( ix, :, 06)
     1 - 0.125 * ( w ( ix, :, 16) - w ( ix, :, 17) )

      u_quadpt ( ix, :, 26) = 0.5 * ( u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 20) ) - 0.25 * w ( ix, :, 07)
     1 + 0.125 * ( w ( ix, :, 18) - w ( ix, :, 19) )

      u_quadpt ( ix, :, 27) = 0.5 * ( u_quadpt ( ix, :, 19)
     1 + u_quadpt ( ix, :, 21) ) - 0.25 * w ( ix, :, 07)
     1 + 0.125 * ( w ( ix, :, 18) + w ( ix, :, 19) )

      u_quadpt ( ix, :, 28) = 0.5 * ( u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 19) ) - 0.25 * w ( ix, :, 06)
     1 + 0.125 * ( w ( ix, :, 16) - w ( ix, :, 17) )

      u_quadpt ( ix, :, 29) = 0.5 * ( u_quadpt ( ix, :, 20)
     1 + u_quadpt ( ix, :, 21) ) - 0.25 * w ( ix, :, 06)
     1 + 0.125 * ( w ( ix, :, 16) + w ( ix, :, 17) )

      u_quadpt ( ix, :, 30) = 0.5 * ( u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 18) ) - 0.25 * w ( ix, :, 05)
     1 - 0.125 * ( w ( ix, :, 14) + w ( ix, :, 15) )

      u_quadpt ( ix, :, 31) = 0.5 * ( u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 19) ) - 0.25 * w ( ix, :, 05)
     1 + 0.125 * ( w ( ix, :, 14) - w ( ix, :, 15) )

      u_quadpt ( ix, :, 32) = 0.5 * ( u_quadpt ( ix, :, 16)
     1 + u_quadpt ( ix, :, 20) ) - 0.25 * w ( ix, :, 05)
     1 - 0.125 * ( w ( ix, :, 14) - w ( ix, :, 15) )

      u_quadpt ( ix, :, 33) = 0.5 * ( u_quadpt ( ix, :, 17)
     1 + u_quadpt ( ix, :, 21) ) - 0.25 * w ( ix, :, 05)
     1 + 0.125 * ( w ( ix, :, 14) + w ( ix, :, 15) )

! ----------

      bfx_const ( ix, 01) = w ( ix, 6, 01) - w ( ix, 6, 05) / 12.0

      bfx_const ( ix, 02) = w ( ix, 6, 01) + 0.5 * w ( ix, 6, 02)
     1 + 0.05 * w ( ix, 6, 11) + w ( ix, 6, 05) / 6.0

      bfx_const ( ix, 03) = w ( ix, 6, 01) + 0.25 * w ( ix, 6, 02)
     1 - 0.021875 * w ( ix, 6, 11) - w ( ix, 6, 05) / 48.0

      bfx_const ( ix, 04) = w ( ix, 6, 01) - 0.25 * w ( ix, 6, 02)
     1 + 0.021875 * w ( ix, 6, 11) - w ( ix, 6, 05) / 48.0

      bfx_const ( ix, 05) = w ( ix, 6, 01) - 0.5 * w ( ix, 6, 02)
     1 - 0.05 * w ( ix, 6, 11) + w ( ix, 6, 05) / 6.0

      bfx_const ( ix, 06) = bfx_const ( ix, 01)

      bfx_const ( ix, 07) = bfx_const ( ix, 01)

      bfx_const ( ix, 08) = bfx_const ( ix, 01)

      bfx_const ( ix, 09) = bfx_const ( ix, 01)

      bfx_const ( ix, 10) = bfx_const ( ix, 01)

      bfx_const ( ix, 11) = bfx_const ( ix, 01)

      bfx_const ( ix, 12) = bfx_const ( ix, 01)

      bfx_const ( ix, 13) = bfx_const ( ix, 01)

      bfx_const ( ix, 14) = bfx_const ( ix, 02)

      bfx_const ( ix, 15) = bfx_const ( ix, 02)

      bfx_const ( ix, 16) = bfx_const ( ix, 02)

      bfx_const ( ix, 17) = bfx_const ( ix, 02)

      bfx_const ( ix, 18) = bfx_const ( ix, 05)

      bfx_const ( ix, 19) = bfx_const ( ix, 05)

      bfx_const ( ix, 20) = bfx_const ( ix, 05)

      bfx_const ( ix, 21) = bfx_const ( ix, 05)


      bfx_const ( ix, 22) = bfx_const ( ix, 02)

      bfx_const ( ix, 23) = bfx_const ( ix, 02)

      bfx_const ( ix, 24) = bfx_const ( ix, 02)

      bfx_const ( ix, 25) = bfx_const ( ix, 02)

      bfx_const ( ix, 26) = bfx_const ( ix, 05)

      bfx_const ( ix, 27) = bfx_const ( ix, 05)

      bfx_const ( ix, 28) = bfx_const ( ix, 05)

      bfx_const ( ix, 29) = bfx_const ( ix, 05)

      bfx_const ( ix, 30) = bfx_const ( ix, 01)

      bfx_const ( ix, 31) = bfx_const ( ix, 01)

      bfx_const ( ix, 32) = bfx_const ( ix, 01)

      bfx_const ( ix, 33) = bfx_const ( ix, 01)

! ----------

      bfy_const ( ix, 01) = w ( ix, 7, 01) - w ( ix, 7, 06) / 12.0

      bfy_const ( ix, 02) = bfy_const ( ix, 01)

      bfy_const ( ix, 03) = bfy_const ( ix, 01)

      bfy_const ( ix, 04) = bfy_const ( ix, 01)

      bfy_const ( ix, 05) = bfy_const ( ix, 01)

      bfy_const ( ix, 06) = w ( ix, 7, 01) + 0.5 * w ( ix, 7, 03)
     1 + 0.05 * w ( ix, 7, 12) + w ( ix, 7, 06) / 6.0

      bfy_const ( ix, 07) = w ( ix, 7, 01) + 0.25 * w ( ix, 7, 03)
     1 - 0.021875 * w ( ix, 7, 12) - w ( ix, 7, 06) / 48.0

      bfy_const ( ix, 08) = w ( ix, 7, 01) - 0.25 * w ( ix, 7, 03)
     1 + 0.021875 * w ( ix, 7, 12) - w ( ix, 7, 06) / 48.0

      bfy_const ( ix, 09) = w ( ix, 7, 01) - 0.5 * w ( ix, 7, 03)
     1 - 0.05 * w ( ix, 7, 12) + w ( ix, 7, 06) / 6.0

      bfy_const ( ix, 10) = bfy_const ( ix, 01)

      bfy_const ( ix, 11) = bfy_const ( ix, 01)

      bfy_const ( ix, 12) = bfy_const ( ix, 01)

      bfy_const ( ix, 13) = bfy_const ( ix, 01)

      bfy_const ( ix, 14) = bfy_const ( ix, 06)

      bfy_const ( ix, 15) = bfy_const ( ix, 09)

      bfy_const ( ix, 16) = bfy_const ( ix, 06)

      bfy_const ( ix, 17) = bfy_const ( ix, 09)

      bfy_const ( ix, 18) = bfy_const ( ix, 06)

      bfy_const ( ix, 19) = bfy_const ( ix, 09)

      bfy_const ( ix, 20) = bfy_const ( ix, 06)

      bfy_const ( ix, 21) = bfy_const ( ix, 09)


      bfy_const ( ix, 22) = bfy_const ( ix, 06)

      bfy_const ( ix, 23) = bfy_const ( ix, 09)

      bfy_const ( ix, 24) = bfy_const ( ix, 01)

      bfy_const ( ix, 25) = bfy_const ( ix, 01)

      bfy_const ( ix, 26) = bfy_const ( ix, 06)

      bfy_const ( ix, 27) = bfy_const ( ix, 09)

      bfy_const ( ix, 28) = bfy_const ( ix, 01)

      bfy_const ( ix, 29) = bfy_const ( ix, 01)

      bfy_const ( ix, 30) = bfy_const ( ix, 06)

      bfy_const ( ix, 31) = bfy_const ( ix, 09)

      bfy_const ( ix, 32) = bfy_const ( ix, 06)

      bfy_const ( ix, 33) = bfy_const ( ix, 09)

! ----------

      bfz_const ( ix, 01) = w ( ix, 8, 01) - w ( ix, 8, 07) / 12.0

      bfz_const ( ix, 02) = bfz_const ( ix, 01)

      bfz_const ( ix, 03) = bfz_const ( ix, 01)

      bfz_const ( ix, 04) = bfz_const ( ix, 01)

      bfz_const ( ix, 05) = bfz_const ( ix, 01)

      bfz_const ( ix, 06) = bfz_const ( ix, 01)

      bfz_const ( ix, 07) = bfz_const ( ix, 01)

      bfz_const ( ix, 08) = bfz_const ( ix, 01)

      bfz_const ( ix, 09) = bfz_const ( ix, 01)

      bfz_const ( ix, 10) = w ( ix, 8, 01) + 0.5 * w ( ix, 8, 04)
     1 + 0.05 * w ( ix, 8, 13) + w ( ix, 8, 07) / 6.0

      bfz_const ( ix, 11) = w ( ix, 8, 01) + 0.25 * w ( ix, 8, 04)
     1 - 0.021875 * w ( ix, 8, 13) - w ( ix, 8, 07) / 48.0

      bfz_const ( ix, 12) = w ( ix, 8, 01) - 0.25 * w ( ix, 8, 04)
     1 + 0.021875 * w ( ix, 8, 13) - w ( ix, 8, 07) / 48.0

      bfz_const ( ix, 13) = w ( ix, 8, 01) - 0.5 * w ( ix, 8, 04)
     1 - 0.05 * w ( ix, 8, 13) + w ( ix, 8, 07) / 6.0

      bfz_const ( ix, 14) = bfz_const ( ix, 10)

      bfz_const ( ix, 15) = bfz_const ( ix, 10)

      bfz_const ( ix, 16) = bfz_const ( ix, 13)

      bfz_const ( ix, 17) = bfz_const ( ix, 13)

      bfz_const ( ix, 18) = bfz_const ( ix, 10)

      bfz_const ( ix, 19) = bfz_const ( ix, 10)

      bfz_const ( ix, 20) = bfz_const ( ix, 13)

      bfz_const ( ix, 21) = bfz_const ( ix, 13)


      bfz_const ( ix, 22) = bfz_const ( ix, 01)

      bfz_const ( ix, 23) = bfz_const ( ix, 01)

      bfz_const ( ix, 24) = bfz_const ( ix, 10)

      bfz_const ( ix, 25) = bfz_const ( ix, 13)

      bfz_const ( ix, 26) = bfz_const ( ix, 01)

      bfz_const ( ix, 27) = bfz_const ( ix, 01)

      bfz_const ( ix, 28) = bfz_const ( ix, 10)

      bfz_const ( ix, 29) = bfz_const ( ix, 13)

      bfz_const ( ix, 30) = bfz_const ( ix, 10)

      bfz_const ( ix, 31) = bfz_const ( ix, 10)

      bfz_const ( ix, 32) = bfz_const ( ix, 13)

      bfz_const ( ix, 33) = bfz_const ( ix, 13)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------

! STEP 3) Get the primitive variables from the conserved variables,
! but only do this for density, velocities and magnetic fields; not
! pressure. These are easy to bring within bounds, so we do that first.

! Actually, it is best to focus on bringing only the density within bounds.

! Here we build the parts of the magnetic field that can be varied without
! losing consistency. These parts are stored in "bfx_var" etc.

! --------------------

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      bfx_var ( ix, iquadpt) = bfx ( ix, iquadpt)
     1                       - bfx_const ( ix, iquadpt)
      bfy_var ( ix, iquadpt) = bfy ( ix, iquadpt)
     1                       - bfy_const ( ix, iquadpt)
      bfz_var ( ix, iquadpt) = bfz ( ix, iquadpt)
     1                       - bfz_const ( ix, iquadpt)

! ----------

! Bring density within bounds.

      IF ( rho ( ix, iquadpt) .GT. vmaxuse ( ix, 1) * 1.000001) THEN
      corrupted_zone ( ix) = 1
      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 1) - w ( ix, 1, 1) )
     1                   / ( rho ( ix, iquadpt) - w ( ix, 1, 1) ) )
      END IF

      IF ( rho ( ix, iquadpt) .LT. vminuse ( ix, 1) * 0.999999) THEN
      corrupted_zone ( ix) = 1
      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 1) - w ( ix, 1, 1) )
     1                   / ( rho ( ix, iquadpt) - w ( ix, 1, 1) ) )
      END IF

! -----
#if ( REDUCE_VEL_BF == 1)

! Bring x-velocity within bounds.

      IF ( vlx ( ix, iquadpt) .GT. vmaxuse ( ix, 2) + smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 2, 1) - w ( ix, 1, 1) * vmaxuse ( ix, 2)
      tempac = u_quadpt ( ix, 2, iquadpt)
     1       - rho ( ix, iquadpt) * vmaxuse ( ix, 2)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

      IF ( vlx ( ix, iquadpt) .LT. vminuse ( ix, 2) - smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 2, 1) - w ( ix, 1, 1) * vminuse ( ix, 2)
      tempac = u_quadpt ( ix, 2, iquadpt)
     1       - rho ( ix, iquadpt) * vminuse ( ix, 2)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

! -----

! Bring y-velocity within bounds.

      IF ( vly ( ix, iquadpt) .GT. vmaxuse ( ix, 3) + smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 3, 1) - w ( ix, 1, 1) * vmaxuse ( ix, 3)
      tempac = u_quadpt ( ix, 3, iquadpt)
     1       - rho ( ix, iquadpt) * vmaxuse ( ix, 3)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

      IF ( vly ( ix, iquadpt) .LT. vminuse ( ix, 3) - smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 3, 1) - w ( ix, 1, 1) * vminuse ( ix, 3)
      tempac = u_quadpt ( ix, 3, iquadpt)
     1       - rho ( ix, iquadpt) * vminuse ( ix, 3)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

! -----

! Bring z-velocity within bounds.

      IF ( vlz ( ix, iquadpt) .GT. vmaxuse ( ix, 4) + smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 4, 1) - w ( ix, 1, 1) * vmaxuse ( ix, 4)
      tempac = u_quadpt ( ix, 4, iquadpt)
     1       - rho ( ix, iquadpt) * vmaxuse ( ix, 4)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

      IF ( vlz ( ix, iquadpt) .LT. vminuse ( ix, 4) - smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 4, 1) - w ( ix, 1, 1) * vminuse ( ix, 4)
      tempac = u_quadpt ( ix, 4, iquadpt)
     1       - rho ( ix, iquadpt) * vminuse ( ix, 4)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

! -----

! Bring Bx within bounds.

      IF ( bfx ( ix, iquadpt) .GT. vmaxuse ( ix, 6) + smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 6) - w ( ix, 6, 1) )
     1                   / ( bfx ( ix, iquadpt) - w ( ix, 6, 1) ) )

      END IF

      IF ( bfx ( ix, iquadpt) .LT. vminuse ( ix, 6) - smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 6) - w ( ix, 6, 1) )
     1                   / ( bfx ( ix, iquadpt) - w ( ix, 6, 1) ) )

      END IF

! -----

! Bring By within bounds.

      IF ( bfy ( ix, iquadpt) .GT. vmaxuse ( ix, 7) + smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 7) - w ( ix, 7, 1) )
     1                   / ( bfy ( ix, iquadpt) - w ( ix, 7, 1) ) )

      END IF

      IF ( bfy ( ix, iquadpt) .LT. vminuse ( ix, 7) - smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 7) - w ( ix, 7, 1) )
     1                   / ( bfy ( ix, iquadpt) - w ( ix, 7, 1) ) )

      END IF

! -----

! Bring Bz within bounds.

      IF ( bfz ( ix, iquadpt) .GT. vmaxuse ( ix, 8) + smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 8) - w ( ix, 8, 1) )
     1                   / ( bfz ( ix, iquadpt) - w ( ix, 8, 1) ) )

      END IF

      IF ( bfz ( ix, iquadpt) .LT. vminuse ( ix, 8) - smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 8) - w ( ix, 8, 1) )
     1                   / ( bfz ( ix, iquadpt) - w ( ix, 8, 1) ) )

      END IF

#endif
! End of "REDUCE_VEL_BF == 1"
! -----

! Lastly, just make sure "reduce_modes" is inside [0,1].

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix), 1.0)
      reduce_modes ( ix) = AMAX1 ( reduce_modes ( ix), 0.0)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 4) Now if we have found a suitable "reduce_modes", we reduce ALL the
! modal variations by that amount. We also reduce the variation in the
! conserved variables by that amount. This leaves us free to rebuild the
! nodal density, velocities and magnetic fields that are within bounds.
! All that remains is to consider the pressure.

! The reduction in velocities and fields also
! gives the pressure variable a second chance to become positive without
! further touch up. We will, eventually, touch it up if it is still
! non-positive.

! --------------------

      DO imode = 2, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      w ( ix, 1: 5, imode) = reduce_modes ( ix) * w ( ix, 1: 5, imode)

#if ( NFLUID >= 1)
      w ( ix, 9: 8 + NFLUID, imode) = reduce_modes ( ix)
     1                              * w ( ix, 9: 8 + NFLUID, imode)
#endif

! -----
#if ( NUMORDER == 2)

      IF ( imode .NE. 2) w ( ix, 6, imode) = reduce_modes ( ix)
     1                 * w ( ix, 6, imode)

      IF ( imode .NE. 3) w ( ix, 7, imode) = reduce_modes ( ix)
     1                 * w ( ix, 7, imode)

      IF ( imode .NE. 4) w ( ix, 8, imode) = reduce_modes ( ix)
     1                 * w ( ix, 8, imode)

#endif
! End of "NUMORDER == 2".
! -----
#if ( NUMORDER == 3)

      IF ( ( imode .NE. 2) .AND. ( imode .NE. 5))
     1                   w ( ix, 6, imode) = reduce_modes ( ix)
     1                 * w ( ix, 6, imode)

      IF ( ( imode .NE. 3) .AND. ( imode .NE. 6))
     1                   w ( ix, 7, imode) = reduce_modes ( ix)
     1                 * w ( ix, 7, imode)

      IF ( ( imode .NE. 4) .AND. ( imode .NE. 7))
     1                   w ( ix, 8, imode) = reduce_modes ( ix)
     1                 * w ( ix, 8, imode)

#endif
! End of "NUMORDER == 3".
! -----
#if ( NUMORDER == 4)

      IF ( ( imode .NE. 2) .AND. ( imode .NE. 5) .AND. ( imode .NE. 11))
     1                   w ( ix, 6, imode) = reduce_modes ( ix)
     1                 * w ( ix, 6, imode)

      IF ( ( imode .NE. 3) .AND. ( imode .NE. 6) .AND. ( imode .NE. 12))
     1                   w ( ix, 7, imode) = reduce_modes ( ix)
     1                 * w ( ix, 7, imode)

      IF ( ( imode .NE. 4) .AND. ( imode .NE. 7) .AND. ( imode .NE. 13))
     1                   w ( ix, 8, imode) = reduce_modes ( ix)
     1                 * w ( ix, 8, imode)

#endif
! End of "NUMORDER == 4".
! -----

      END IF ! End of "corrupted_zone ( ix) .EQ. 1".

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO imode =" loop.

! --------------------

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      tempaa = 1.0 - reduce_modes ( ix)

      u_quadpt ( ix, 1: 5, iquadpt) = tempaa * w ( ix, 1: 5, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 1: 5, iquadpt)

#if ( NFLUID >= 1)
      u_quadpt ( ix, 9: 8 + NFLUID, iquadpt) =
     1                       tempaa * w ( ix, 9: 8 + NFLUID, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 9: 8 + NFLUID, iquadpt)
#endif

! -----

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempab = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempab
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempab
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempab

! -----

      bfx_var ( ix, iquadpt) = reduce_modes ( ix)
     1                       * bfx_var ( ix, iquadpt)
      bfy_var ( ix, iquadpt) = reduce_modes ( ix)
     1                       * bfy_var ( ix, iquadpt)
      bfz_var ( ix, iquadpt) = reduce_modes ( ix)
     1                       * bfz_var ( ix, iquadpt)

      bfx ( ix, iquadpt) = bfx_const ( ix, iquadpt)
     1                   + bfx_var ( ix, iquadpt)
      bfy ( ix, iquadpt) = bfy_const ( ix, iquadpt)
     1                   + bfy_var ( ix, iquadpt)
      bfz ( ix, iquadpt) = bfz_const ( ix, iquadpt)
     1                   + bfz_var ( ix, iquadpt)

      u_quadpt ( ix, 6, iquadpt) = bfx ( ix, iquadpt)
      u_quadpt ( ix, 7, iquadpt) = bfy ( ix, iquadpt)
      u_quadpt ( ix, 8, iquadpt) = bfz ( ix, iquadpt)

! -----

      END IF ! End of "corrupted_zone ( ix) .EQ. 1".

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 5) The zones with "corrupted_zone ( ix) == 1" now have taken on
! corrections, the rest have not. Build pressure for all zones here.
! Flag "pp_zone" for zones with pressures that go below floor values.

! Build the coefficients of the cubic.

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

! Now build pressure and species fractions and flag if needed.

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      prs ( ix, iquadpt) = ( u_quadpt ( ix, 5, iquadpt) - tempab
     1         - tempac) * ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      IF ( prs ( ix, iquadpt) .LT. vminuse ( ix, 5) ) THEN

      corrupted_zone ( ix) = 1
      pp_zone ( ix) = 1
      pp_quadpt ( ix, iquadpt) = 1

      bx_c = bfx_const ( ix, iquadpt) * ob_root_4pi
      by_c = bfy_const ( ix, iquadpt) * ob_root_4pi
      bz_c = bfz_const ( ix, iquadpt) * ob_root_4pi

      bx_v = bfx_var ( ix, iquadpt) * ob_root_4pi
      by_v = bfy_var ( ix, iquadpt) * ob_root_4pi
      bz_v = bfz_var ( ix, iquadpt) * ob_root_4pi

      del_rho = rho ( ix, iquadpt) - w ( ix, 1, 1)
      del_eng = eng_t ( ix, iquadpt) - w ( ix, 5, 1)

      del_mx = u_quadpt ( ix, 2, iquadpt) - w ( ix, 2, 1)
      del_my = u_quadpt ( ix, 3, iquadpt) - w ( ix, 3, 1)
      del_mz = u_quadpt ( ix, 4, iquadpt) - w ( ix, 4, 1)

      m_bar2 = w ( ix, 2, 1)**2 + w ( ix, 3, 1)**2 + w ( ix, 4, 1)**2
      del_m2 = del_mx**2 + del_my**2 + del_mz**2
      m_dot_del_m = w ( ix, 2, 1) * del_mx + w ( ix, 3, 1) * del_my
     1            + w ( ix, 4, 1) * del_mz

      b_bar2 = bx_c**2 + by_c**2 + bz_c**2
      del_b2 = bx_v**2 + by_v**2 + bz_v**2
      b_dot_del_b = bx_c * bx_v + by_c * by_v + bz_c * bz_v

      tempad = 1.0 / ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

      a3 ( ix, iquadpt) = - del_rho * del_b2

      a2 ( ix, iquadpt) = 2.0 * del_rho * del_eng - del_m2
     1                  - w ( ix, 1, 1) * del_b2
     1                  - 2.0 * del_rho * b_dot_del_b

      a1 ( ix, iquadpt) = 2.0 * w ( ix, 1, 1) * del_eng
     1                  + 2.0 * w ( ix, 5, 1) * del_rho
     1   - 2.0 * m_dot_del_m - 2.0 * w ( ix, 1, 1) * b_dot_del_b
     1   - del_rho * b_bar2 - 2.0 * vminuse ( ix, 5) * tempad * del_rho

      a0 ( ix, iquadpt) = 2.0 * w ( ix, 1, 1) * w ( ix, 5, 1) - m_bar2
     1   - w ( ix, 1, 1) * b_bar2
     1   - 2.0 * vminuse ( ix, 5) * tempad * w ( ix, 1, 1)

      x_low ( ix, iquadpt) = 0.0
      x_high ( ix, iquadpt) = 1.0
      x_root ( ix, iquadpt) = 0.5

      f_low ( ix, iquadpt) = a0 ( ix, iquadpt)
     1                 + a1 ( ix, iquadpt) * x_low ( ix, iquadpt)
     1                 + a2 ( ix, iquadpt) * x_low ( ix, iquadpt)**2
     1                 + a3 ( ix, iquadpt) * x_low ( ix, iquadpt)**3

      f_high ( ix, iquadpt) = a0 ( ix, iquadpt)
     1                 + a1 ( ix, iquadpt) * x_high ( ix, iquadpt)
     1                 + a2 ( ix, iquadpt) * x_high ( ix, iquadpt)**2
     1                 + a3 ( ix, iquadpt) * x_high ( ix, iquadpt)**3

      f_root ( ix, iquadpt) = a0 ( ix, iquadpt)
     1                 + a1 ( ix, iquadpt) * x_root ( ix, iquadpt)
     1                 + a2 ( ix, iquadpt) * x_root ( ix, iquadpt)**2
     1                 + a3 ( ix, iquadpt) * x_root ( ix, iquadpt)**3

      df_root ( ix, iquadpt) = a1 ( ix, iquadpt)
     1           + 2.0 * a2 ( ix, iquadpt) * x_root ( ix, iquadpt)
     1           + 3.0 * a3 ( ix, iquadpt) * x_root ( ix, iquadpt)**2

      flip_sign ( ix, iquadpt) = 1.0
      IF ( f_high ( ix, iquadpt) < f_low ( ix, iquadpt) )
     1                             flip_sign ( ix, iquadpt) = - 1.0

      END IF ! End of "prs ( ix, iquadpt) .LT. vminuse ( ix, 5)".

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 6) Now iterate the roots of the cubic to convergence, or at least
! iterate them a few times.

      DO iter = 1, itermax

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_quadpt ( ix, iquadpt) .EQ. 1) THEN

! ----------

      tempaa = ( x_root ( ix, iquadpt) - x_low ( ix, iquadpt) )
     1      * df_root ( ix, iquadpt) - f_root ( ix, iquadpt)
      tempab = ( x_high ( ix, iquadpt) - x_root ( ix, iquadpt) )
     1      * df_root ( ix, iquadpt) + f_root ( ix, iquadpt)
      tempac = df_root ( ix, iquadpt)
     1      * ( x_high ( ix, iquadpt) - x_low ( ix, iquadpt) )

      IF ( ( tempaa * tempab > 0.0 ) .AND.
     1   ( ABS ( 2.0 * f_root ( ix, iquadpt) ) < ABS ( tempac) ) ) THEN
! Choosing Newton here.
          x_root ( ix, iquadpt) = x_root ( ix, iquadpt)
     1  - f_root ( ix, iquadpt) / df_root ( ix, iquadpt)
      ELSE
! Choosing Bisection here.
          x_root ( ix, iquadpt) = 0.5 * ( x_low ( ix, iquadpt)
     1                                  + x_high ( ix, iquadpt) )
      END IF

! -----

! Evaluate the function & its derivativative.

      f_root ( ix, iquadpt) = a0 ( ix, iquadpt)
     1                 + a1 ( ix, iquadpt) * x_root ( ix, iquadpt)
     1                 + a2 ( ix, iquadpt) * x_root ( ix, iquadpt)**2
     1                 + a3 ( ix, iquadpt) * x_root ( ix, iquadpt)**3

      df_root ( ix, iquadpt) = a1 ( ix, iquadpt)
     1           + 2.0 * a2 ( ix, iquadpt) * x_root ( ix, iquadpt)
     1           + 3.0 * a3 ( ix, iquadpt) * x_root ( ix, iquadpt)**2

! -----

! Shrink the domain.

      IF ( flip_sign ( ix, iquadpt) * f_root ( ix, iquadpt) > 0) THEN
          x_high ( ix, iquadpt) = x_root ( ix, iquadpt)
          f_high ( ix, iquadpt) = f_root ( ix, iquadpt)
      ELSE
          x_low ( ix, iquadpt) = x_root ( ix, iquadpt)
          f_low ( ix, iquadpt) = f_root ( ix, iquadpt)
      END IF

! ----------

      END IF ! End of "pp_quadpt ( ix, iquadpt) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

      END DO
! End of "DO iter =" loop.

! ----------------------------------------------------------------------

! STEP 7) By this point we know how much further we need to shrink the modes.
! Poll all the nodes of a flagged zone to find the smallest value of
! "x_root ( ix, iquadpt)" over all the quadrature points in that zone.

! To begin with, reset "reduce_modes". Then minimize "reduce_modes".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      reduce_modes ( ix) = 1.0
      END DO


      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_quadpt ( ix, iquadpt) .EQ. 1) THEN

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                             0.95 * x_root ( ix, iquadpt) )

      END IF ! End of "pp_quadpt ( ix, iquadpt) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 8) This is almost a repeat of STEP 4, where we reduce all the
! possible modes that can be reduced by a factor of "reduce_modes".
! This time we key off "pp_zone".

! --------------------

      DO imode = 2, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( pp_zone ( ix) .EQ. 1) THEN

! -----

      w ( ix, 1: 5, imode) = reduce_modes ( ix) * w ( ix, 1: 5, imode)

#if ( NFLUID >= 1)
      w ( ix, 9: 8 + NFLUID, imode) = reduce_modes ( ix)
     1                              * w ( ix, 9: 8 + NFLUID, imode)
#endif

! -----
#if ( NUMORDER == 2)

      IF ( imode .NE. 2) w ( ix, 6, imode) = reduce_modes ( ix)
     1                 * w ( ix, 6, imode)

      IF ( imode .NE. 3) w ( ix, 7, imode) = reduce_modes ( ix)
     1                 * w ( ix, 7, imode)

      IF ( imode .NE. 4) w ( ix, 8, imode) = reduce_modes ( ix)
     1                 * w ( ix, 8, imode)

#endif
! End of "NUMORDER == 2".
! -----
#if ( NUMORDER == 3)

      IF ( ( imode .NE. 2) .AND. ( imode .NE. 5))
     1                   w ( ix, 6, imode) = reduce_modes ( ix)
     1                 * w ( ix, 6, imode)

      IF ( ( imode .NE. 3) .AND. ( imode .NE. 6))
     1                   w ( ix, 7, imode) = reduce_modes ( ix)
     1                 * w ( ix, 7, imode)

      IF ( ( imode .NE. 4) .AND. ( imode .NE. 7))
     1                   w ( ix, 8, imode) = reduce_modes ( ix)
     1                 * w ( ix, 8, imode)

#endif
! End of "NUMORDER == 3".
! -----
#if ( NUMORDER == 4)

      IF ( ( imode .NE. 2) .AND. ( imode .NE. 5) .AND. ( imode .NE. 11))
     1                   w ( ix, 6, imode) = reduce_modes ( ix)
     1                 * w ( ix, 6, imode)

      IF ( ( imode .NE. 3) .AND. ( imode .NE. 6) .AND. ( imode .NE. 12))
     1                   w ( ix, 7, imode) = reduce_modes ( ix)
     1                 * w ( ix, 7, imode)

      IF ( ( imode .NE. 4) .AND. ( imode .NE. 7) .AND. ( imode .NE. 13))
     1                   w ( ix, 8, imode) = reduce_modes ( ix)
     1                 * w ( ix, 8, imode)

#endif
! End of "NUMORDER == 4".
! -----

      END IF ! End of "pp_zone ( ix) .EQ. 1".

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO imode =" loop.

! ----------------------------------------------------------------------

! STEP 9) Now make sure that the nodal values in the zones that originally
! lost pressure positivity at least have a pressure that exceeds 
! "vminuse ( ix, 5)".

! Rebuild the total energy density at those nodes, i.e. rebuild
! "u_quadpt ( ix, 5, iquadpt)" for the flagged nodes.

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( pp_zone ( ix) .EQ. 1) THEN

! -----

      tempaa = 1.0 - reduce_modes ( ix)

      u_quadpt ( ix, 1: 5, iquadpt) = tempaa * w ( ix, 1: 5, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 1: 5, iquadpt)

#if ( NFLUID >= 1)
      u_quadpt ( ix, 9: 8 + NFLUID, iquadpt) =
     1                       tempaa * w ( ix, 9: 8 + NFLUID, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 9: 8 + NFLUID, iquadpt)
#endif

! -----

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempab = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempab
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempab
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempab

! -----

      bfx_var ( ix, iquadpt) = reduce_modes ( ix)
     1                       * bfx_var ( ix, iquadpt)
      bfy_var ( ix, iquadpt) = reduce_modes ( ix)
     1                       * bfy_var ( ix, iquadpt)
      bfz_var ( ix, iquadpt) = reduce_modes ( ix)
     1                       * bfz_var ( ix, iquadpt)

      bfx ( ix, iquadpt) = bfx_const ( ix, iquadpt)
     1                   + bfx_var ( ix, iquadpt)
      bfy ( ix, iquadpt) = bfy_const ( ix, iquadpt)
     1                   + bfy_var ( ix, iquadpt)
      bfz ( ix, iquadpt) = bfz_const ( ix, iquadpt)
     1                   + bfz_var ( ix, iquadpt)

      u_quadpt ( ix, 6, iquadpt) = bfx ( ix, iquadpt)
      u_quadpt ( ix, 7, iquadpt) = bfy ( ix, iquadpt)
      u_quadpt ( ix, 8, iquadpt) = bfz ( ix, iquadpt)

! -----

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      prs ( ix, iquadpt) = ( u_quadpt ( ix, 5, iquadpt) - tempab
     1         - tempac) * ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

! This ensuing line is the only place where we might cause a loss of
! total energy conservation.
      prs ( ix, iquadpt)= AMAX1 ( prs ( ix, iquadpt), vminuse ( ix, 5))

      u_quadpt ( ix, 5, iquadpt) = tempab + tempac
     1  + prs ( ix, iquadpt) / ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

! -----

      END IF ! End of "pp_zone ( ix) .EQ. 1".

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 10) Now use the nodal values of just the total energy density in
! the affected zone to rebuild ALL the moments of the energy density
! (and just the energy density) in that zone.

! Notice that we also rebuild the zeroth moment.
! If "AMAX1 ( prs ( ix, iquadpt), vminuse ( ix, 5))" had caused any changes
! in the pressure in the previous step then we will also be changing
! the total energy density in the zone. In that case, we choose positivity
! over and above total energy conservation.

! We also build the entropy variable and its first moments for all the zones.

! --------------------------------------------------

      ivar = 5

! --------------------------------------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_zone ( ix) .EQ. 1) THEN

      w ( ix, ivar, 01) = 0.125 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 02) + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 04) + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 06) + u_quadpt ( ix, ivar, 07)
     1 + u_quadpt ( ix, ivar, 08) )

      w ( ix, ivar, 02) =  0.25 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 03) + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 07) - u_quadpt ( ix, ivar, 02)
     1 - u_quadpt ( ix, ivar, 04) - u_quadpt ( ix, ivar, 06)
     1 - u_quadpt ( ix, ivar, 08) )

      w ( ix, ivar, 03) =  0.25 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 02) + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 06) - u_quadpt ( ix, ivar, 03)
     1 - u_quadpt ( ix, ivar, 04) - u_quadpt ( ix, ivar, 07)
     1 - u_quadpt ( ix, ivar, 08) )

      w ( ix, ivar, 04) =  0.25 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 02) + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 04) - u_quadpt ( ix, ivar, 05)
     1 - u_quadpt ( ix, ivar, 06) - u_quadpt ( ix, ivar, 07)
     1 - u_quadpt ( ix, ivar, 08) )

      END IF ! End of "pp_zone ( ix) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ent ( ix, 01) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 02) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 03) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 04) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 05) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 04)

      ent ( ix, 06) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 04)

! ----------

      IF ( igeom .EQ. 1) THEN
      ax1 = 1.0
      ax2 = 1.0
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove / dx_use ( ix)
      ob_dy = dtmove / dy_use ( ix)
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 2) THEN
      ax1 = xb_use ( ix)
      ax2 = xb_use ( ix - 1)
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 2.0 / ( xb_use ( ix)**2 - xb_use ( ix - 1)**2)
      ob_dy = dtmove / ( xc_use ( ix) * dy_use ( ix))
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 3) THEN
      ax1 = xb_use ( ix)**2
      ax2 = xb_use ( ix - 1)**2
      ay1 = SIN ( yb_use ( ix))
      ay2 = SIN ( yb_use_m1 ( ix))
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 3.0 / ( xb_use ( ix)**3 - xb_use ( ix - 1)**3)
      ob_dy = dtmove / ( COS ( yb_use_m1 ( ix)) - COS ( yb_use ( ix)) )
      ob_dz = dtmove / ( xc_use ( ix) * SIN ( yc_use ( ix))
     1                                * dz_use ( ix))
      END IF

      entropy_1d ( ix, 05) = - 0.25 * ax1 * ob_dx * ent ( ix, 01)
     1*( vlx ( ix, 01) + vlx ( ix, 03) + vlx ( ix, 05) + vlx ( ix, 07))
     1                       + 0.25 * ax2 * ob_dx * ent ( ix, 02)
     1*( vlx ( ix, 02) + vlx ( ix, 04) + vlx ( ix, 06) + vlx ( ix, 08))

     1                       - 0.25 * ay1 * ob_dy * ent ( ix, 03)
     1*( vly ( ix, 01) + vly ( ix, 02) + vly ( ix, 05) + vly ( ix, 06))
     1                       + 0.25 * ay2 * ob_dy * ent ( ix, 04)
     1*( vly ( ix, 03) + vly ( ix, 04) + vly ( ix, 07) + vly ( ix, 08))

     1                       - 0.25 * az1 * ob_dz * ent ( ix, 05)
     1*( vlz ( ix, 01) + vlz ( ix, 02) + vlz ( ix, 03) + vlz ( ix, 04))
     1                       + 0.25 * az2 * ob_dz * ent ( ix, 06)
     1*( vlz ( ix, 05) + vlz ( ix, 06) + vlz ( ix, 07) + vlz ( ix, 08))

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".
! --------------------------------------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_zone ( ix) .EQ. 1) THEN

      tu ( ix, ivar, 01) = u_quadpt ( ix, ivar, 02)
     1 + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 04) + u_quadpt ( ix, ivar, 05)

      tu ( ix, ivar, 02) = u_quadpt ( ix, ivar, 06)
     1 + u_quadpt ( ix, ivar, 07)
     1 + u_quadpt ( ix, ivar, 08) + u_quadpt ( ix, ivar, 09)

      tu ( ix, ivar, 03) = u_quadpt ( ix, ivar, 10)
     1 + u_quadpt ( ix, ivar, 11)
     1 + u_quadpt ( ix, ivar, 12) + u_quadpt ( ix, ivar, 13)

      tu ( ix, ivar, 04) = u_quadpt ( ix, ivar, 02)
     1 + u_quadpt ( ix, ivar, 06)
     1 + u_quadpt ( ix, ivar, 10) + u_quadpt ( ix, ivar, 12)

      tu ( ix, ivar, 05) = u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 07)
     1 + u_quadpt ( ix, ivar, 11) + u_quadpt ( ix, ivar, 13)

      tu ( ix, ivar, 06) = u_quadpt ( ix, ivar, 04)
     1 + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 08) + u_quadpt ( ix, ivar, 09)

      tu ( ix, ivar, 07) = u_quadpt ( ix, ivar, 04)
     1 + u_quadpt ( ix, ivar, 08)
     1 + u_quadpt ( ix, ivar, 10) + u_quadpt ( ix, ivar, 11)

      tu ( ix, ivar, 08) = u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 09)
     1 + u_quadpt ( ix, ivar, 12) + u_quadpt ( ix, ivar, 13)

      tu ( ix, ivar, 09) = u_quadpt ( ix, ivar, 02)
     1 + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 06) + u_quadpt ( ix, ivar, 07)

      w ( ix, ivar, 02) = 0.25 * ( tu ( ix, ivar, 01)
     1 - tu ( ix, ivar, 02) )

      w ( ix, ivar, 03) = 0.25 * ( tu ( ix, ivar, 04)
     1 - tu ( ix, ivar, 05) )

      w ( ix, ivar, 04) = 0.25 * ( tu ( ix, ivar, 07)
     1 - tu ( ix, ivar, 08) )

      w ( ix, ivar, 05) = 0.5 * ( tu ( ix, ivar, 01)
     1 + tu ( ix, ivar, 02)
     1 - tu ( ix, ivar, 03) ) - 2.0 * u_quadpt ( ix, ivar, 01)

      w ( ix, ivar, 06) = 0.5 * ( tu ( ix, ivar, 04)
     1 + tu ( ix, ivar, 05)
     1 - tu ( ix, ivar, 06) ) - 2.0 * u_quadpt ( ix, ivar, 01)

      w ( ix, ivar, 07) = 0.5 * ( tu ( ix, ivar, 07)
     1 + tu ( ix, ivar, 08)
     1 - tu ( ix, ivar, 09) ) - 2.0 * u_quadpt ( ix, ivar, 01)

      w ( ix, ivar, 08) = u_quadpt ( ix, ivar, 02)
     1 - u_quadpt ( ix, ivar, 03)
     1 - u_quadpt ( ix, ivar, 06) + u_quadpt ( ix, ivar, 07)

      w ( ix, ivar, 09) = u_quadpt ( ix, ivar, 10)
     1 - u_quadpt ( ix, ivar, 11)
     1 - u_quadpt ( ix, ivar, 12) + u_quadpt ( ix, ivar, 13)

      w ( ix, ivar, 10) = u_quadpt ( ix, ivar, 04)
     1 - u_quadpt ( ix, ivar, 05)
     1 - u_quadpt ( ix, ivar, 08) + u_quadpt ( ix, ivar, 09)

      w ( ix, ivar, 01) = u_quadpt ( ix, ivar, 01)
     1 + ( w ( ix, ivar, 05)
     1 + w ( ix, ivar, 06) + w ( ix, ivar, 07) ) / 12.0

      END IF ! End of "pp_zone ( ix) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ent ( ix, 01) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 02) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 03) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 04) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 05) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 04)

      ent ( ix, 06) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 04)

! ----------

      IF ( igeom .EQ. 1) THEN
      ax1 = 1.0
      ax2 = 1.0
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove / dx_use ( ix)
      ob_dy = dtmove / dy_use ( ix)
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 2) THEN
      ax1 = xb_use ( ix)
      ax2 = xb_use ( ix - 1)
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 2.0 / ( xb_use ( ix)**2 - xb_use ( ix - 1)**2)
      ob_dy = dtmove / ( xc_use ( ix) * dy_use ( ix))
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 3) THEN
      ax1 = xb_use ( ix)**2
      ax2 = xb_use ( ix - 1)**2
      ay1 = SIN ( yb_use ( ix))
      ay2 = SIN ( yb_use_m1 ( ix))
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 3.0 / ( xb_use ( ix)**3 - xb_use ( ix - 1)**3)
      ob_dy = dtmove / ( COS ( yb_use_m1 ( ix)) - COS ( yb_use ( ix)) )
      ob_dz = dtmove / ( xc_use ( ix) * SIN ( yc_use ( ix))
     1                                * dz_use ( ix))
      END IF

      entropy_1d ( ix, 05) = - 0.25 * ax1 * ob_dx * ent ( ix, 01)
     1*( vlx ( ix, 14) + vlx ( ix, 16) + vlx ( ix, 18) + vlx ( ix, 20))
     1                       + 0.25 * ax2 * ob_dx * ent ( ix, 02)
     1*( vlx ( ix, 15) + vlx ( ix, 17) + vlx ( ix, 19) + vlx ( ix, 21))

     1                       - 0.25 * ay1 * ob_dy * ent ( ix, 03)
     1*( vly ( ix, 14) + vly ( ix, 15) + vly ( ix, 18) + vly ( ix, 19))
     1                       + 0.25 * ay2 * ob_dy * ent ( ix, 04)
     1*( vly ( ix, 16) + vly ( ix, 17) + vly ( ix, 20) + vly ( ix, 21))

     1                       - 0.25 * az1 * ob_dz * ent ( ix, 05)
     1*( vlz ( ix, 14) + vlz ( ix, 15) + vlz ( ix, 16) + vlz ( ix, 17))
     1                       + 0.25 * az2 * ob_dz * ent ( ix, 06)
     1*( vlz ( ix, 18) + vlz ( ix, 19) + vlz ( ix, 20) + vlz ( ix, 21))

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".
! --------------------------------------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_zone ( ix) .EQ. 1) THEN

      w ( ix, ivar, 05) = 2.0 * ( u_quadpt ( ix, ivar, 02)
     1 - 2.0 * u_quadpt ( ix, ivar, 01) + u_quadpt ( ix, ivar, 05) )

      w ( ix, ivar, 06) = 2.0 * ( u_quadpt ( ix, ivar, 06)
     1 - 2.0 * u_quadpt ( ix, ivar, 01) + u_quadpt ( ix, ivar, 09) )

      w ( ix, ivar, 07) = 2.0 * ( u_quadpt ( ix, ivar, 10)
     1 - 2.0 * u_quadpt ( ix, ivar, 01) + u_quadpt ( ix, ivar, 13) )

      w ( ix, ivar, 08) = 0.5 * ( u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 16)
     1 - u_quadpt ( ix, ivar, 17) - u_quadpt ( ix, ivar, 18)
     1 + u_quadpt ( ix, ivar, 19) - u_quadpt ( ix, ivar, 20)
     1 + u_quadpt ( ix, ivar, 21) )

      w ( ix, ivar, 09) = 0.5 * ( u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 18)
     1 - u_quadpt ( ix, ivar, 20) - u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 17) - u_quadpt ( ix, ivar, 19)
     1 + u_quadpt ( ix, ivar, 21) )

      w ( ix, ivar, 10) = 0.5 * ( u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 15)
     1 - u_quadpt ( ix, ivar, 19) - u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 20) - u_quadpt ( ix, ivar, 17)
     1 + u_quadpt ( ix, ivar, 21) )

      w ( ix, ivar, 11) = 16 * ( u_quadpt ( ix, ivar, 02)
     1 - 2.0 * u_quadpt ( ix, ivar, 03)
     1 + 2.0 * u_quadpt ( ix, ivar, 04)
     1 - u_quadpt ( ix, ivar, 05) ) / 3.0

      w ( ix, ivar, 12) = 16.0 * ( u_quadpt ( ix, ivar, 06)
     1 - 2.0 * u_quadpt ( ix, ivar, 07)
     1 + 2.0 * u_quadpt ( ix, ivar, 08)
     1 - u_quadpt ( ix, ivar, 09) ) / 3.0

      w ( ix, ivar, 13) = 16.0 * ( u_quadpt ( ix, ivar, 10)
     1 - 2.0 * u_quadpt ( ix, ivar, 11)
     1 + 2.0 * u_quadpt ( ix, ivar, 12)
     1 - u_quadpt ( ix, ivar, 13) ) / 3.0

      tu ( ix, ivar, 01) = u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 17)

      tu ( ix, ivar, 02) = u_quadpt ( ix, ivar, 18)
     1 + u_quadpt ( ix, ivar, 19)
     1 + u_quadpt ( ix, ivar, 20) + u_quadpt ( ix, ivar, 21)

      tu ( ix, ivar, 03) = u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 20)

      tu ( ix, ivar, 04) = u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 17)
     1 + u_quadpt ( ix, ivar, 19) + u_quadpt ( ix, ivar, 21)

      tu ( ix, ivar, 05) = u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 19)

      tu ( ix, ivar, 06) = u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 17)
     1 + u_quadpt ( ix, ivar, 20) + u_quadpt ( ix, ivar, 21)

      w ( ix, ivar, 02) = 2.0 * ( u_quadpt ( ix, ivar, 02)
     1 - u_quadpt ( ix, ivar, 05) ) / 3.0 - 0.1 * w ( ix, ivar, 11)
     1 + 0.25 * ( tu ( ix, ivar, 01) - tu ( ix, ivar, 02) ) / 3.0

      w ( ix, ivar, 03) = 2.0 * ( u_quadpt ( ix, ivar, 06)
     1 - u_quadpt ( ix, ivar, 09) ) / 3.0 - 0.1 * w ( ix, ivar, 12)
     1 + 0.25 * ( tu ( ix, ivar, 03) - tu ( ix, ivar, 04) ) / 3.0

      w ( ix, ivar, 04) = 2.0 * ( u_quadpt ( ix, ivar, 10)
     1 - u_quadpt ( ix, ivar, 13) ) / 3.0 - 0.1 * w ( ix, ivar, 13)
     1 + 0.25 * ( tu ( ix, ivar, 05) - tu ( ix, ivar, 06) ) / 3.0

      w ( ix, ivar, 20) = u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 18)
     1 - u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 20)
     1 - u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 19)
     1 + u_quadpt ( ix, ivar, 17) - u_quadpt ( ix, ivar, 21)

      w ( ix, ivar, 01) = 0.125 * ( u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 17) + u_quadpt ( ix, ivar, 18)
     1 + u_quadpt ( ix, ivar, 19) + u_quadpt ( ix, ivar, 20)
     1 + u_quadpt ( ix, ivar, 21) ) - 0.5 * ( w ( ix, ivar, 05)
     1 + w ( ix, ivar, 06) + w ( ix, ivar, 07) ) / 3.0


      tu ( ix, ivar, 22) = 2 * u_quadpt ( ix, ivar, 22)
     1 - ( u_quadpt ( ix, ivar, 14) + u_quadpt ( ix, ivar, 16) )

      tu ( ix, ivar, 23) = 2 * u_quadpt ( ix, ivar, 23)
     1 - ( u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 17) )

      tu ( ix, ivar, 24) = 2 * u_quadpt ( ix, ivar, 24)
     1 - ( u_quadpt ( ix, ivar, 14) + u_quadpt ( ix, ivar, 15) )

      tu ( ix, ivar, 25) = 2 * u_quadpt ( ix, ivar, 25)
     1 - ( u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 17) )

      tu ( ix, ivar, 26) = 2 * u_quadpt ( ix, ivar, 26)
     1 - ( u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 20) )

      tu ( ix, ivar, 27) = 2 * u_quadpt ( ix, ivar, 27)
     1 - ( u_quadpt ( ix, ivar, 19) + u_quadpt ( ix, ivar, 21) )

      tu ( ix, ivar, 28) = 2 * u_quadpt ( ix, ivar, 28)
     1 - ( u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 19) )

      tu ( ix, ivar, 29) = 2 * u_quadpt ( ix, ivar, 29)
     1 - ( u_quadpt ( ix, ivar, 20) + u_quadpt ( ix, ivar, 21) )

      tu ( ix, ivar, 30) = 2 * u_quadpt ( ix, ivar, 30)
     1 - ( u_quadpt ( ix, ivar, 14) + u_quadpt ( ix, ivar, 18) )

      tu ( ix, ivar, 31) = 2 * u_quadpt ( ix, ivar, 31)
     1 - ( u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 19) )

      tu ( ix, ivar, 32) = 2 * u_quadpt ( ix, ivar, 32)
     1 - ( u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 20) )

      tu ( ix, ivar, 33) = 2 * u_quadpt ( ix, ivar, 33)
     1 - ( u_quadpt ( ix, ivar, 17) + u_quadpt ( ix, ivar, 21) )


      w ( ix, ivar, 14) = tu ( ix, ivar, 31) + tu ( ix, ivar, 33)
     1 - tu ( ix, ivar, 30) - tu ( ix, ivar, 32)

      w ( ix, ivar, 15) = tu ( ix, ivar, 32) + tu ( ix, ivar, 33)
     1 - tu ( ix, ivar, 30) - tu ( ix, ivar, 31)

      w ( ix, ivar, 16) = tu ( ix, ivar, 28) + tu ( ix, ivar, 29)
     1 - tu ( ix, ivar, 24) - tu ( ix, ivar, 25)

      w ( ix, ivar, 17) = tu ( ix, ivar, 25) + tu ( ix, ivar, 29)
     1 - tu ( ix, ivar, 24) - tu ( ix, ivar, 28)

      w ( ix, ivar, 18) = tu ( ix, ivar, 26) + tu ( ix, ivar, 27)
     1 - tu ( ix, ivar, 22) - tu ( ix, ivar, 23)

      w ( ix, ivar, 19) = tu ( ix, ivar, 27) - tu ( ix, ivar, 26)
     1 + tu ( ix, ivar, 23) - tu ( ix, ivar, 22)


      END IF ! End of "pp_zone ( ix) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ent ( ix, 01) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 02) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 03) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 04) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 05) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 04)

      ent ( ix, 06) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 04)

! ----------

      IF ( igeom .EQ. 1) THEN
      ax1 = 1.0
      ax2 = 1.0
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove / dx_use ( ix)
      ob_dy = dtmove / dy_use ( ix)
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 2) THEN
      ax1 = xb_use ( ix)
      ax2 = xb_use ( ix - 1)
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 2.0 / ( xb_use ( ix)**2 - xb_use ( ix - 1)**2)
      ob_dy = dtmove / ( xc_use ( ix) * dy_use ( ix))
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 3) THEN
      ax1 = xb_use ( ix)**2
      ax2 = xb_use ( ix - 1)**2
      ay1 = SIN ( yb_use ( ix))
      ay2 = SIN ( yb_use_m1 ( ix))
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 3.0 / ( xb_use ( ix)**3 - xb_use ( ix - 1)**3)
      ob_dy = dtmove / ( COS ( yb_use_m1 ( ix)) - COS ( yb_use ( ix)) )
      ob_dz = dtmove / ( xc_use ( ix) * SIN ( yc_use ( ix))
     1                                * dz_use ( ix))
      END IF

      entropy_1d ( ix, 05) = - 0.25 * ax1 * ob_dx * ent ( ix, 01)
     1*( vlx ( ix, 14) + vlx ( ix, 15) + vlx ( ix, 16) + vlx ( ix, 17))
     1                       + 0.25 * ax2 * ob_dx * ent ( ix, 02)
     1*( vlx ( ix, 18) + vlx ( ix, 19) + vlx ( ix, 20) + vlx ( ix, 21))

     1                       - 0.25 * ay1 * ob_dy * ent ( ix, 03)
     1*( vly ( ix, 14) + vly ( ix, 18) + vly ( ix, 16) + vly ( ix, 20))
     1                       + 0.25 * ay2 * ob_dy * ent ( ix, 04)
     1*( vly ( ix, 15) + vly ( ix, 19) + vly ( ix, 17) + vly ( ix, 21))

     1                       - 0.25 * az1 * ob_dz * ent ( ix, 05)
     1*( vlz ( ix, 14) + vlz ( ix, 18) + vlz ( ix, 15) + vlz ( ix, 19))
     1                       + 0.25 * az2 * ob_dz * ent ( ix, 06)
     1*( vlz ( ix, 16) + vlz ( ix, 20) + vlz ( ix, 17) + vlz ( ix, 21))

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".
! ----------------------------------------------------------------------

      END SUBROUTINE RECONST_POS_RHO_PRS_MHD_CUBIC

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE RECONST_POS_RHO_PRS_MHD_QUAD

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove, igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gasconst,

     1             xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1             xb_use, yb_use, yb_use_m1, zb_use,

     1             flatten, w, entropy_1d,
     1             gamma_zonecen, molewt_zonecen,
     1             vmin, vmax, corrupted_zone)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Checks the reconstruction for density and pressure positivity. Without
! these the fluxes become non-convex and unphysical. Comes in with
! the following 1d arrays loaded:

! xc_use, yc_use, zc_use, dx_use, dy_use, dz_use : geometric 1d variables.
! flatten : Has the value from "flattengr".
! w : modal vector of conserved variables.
! entropy_1d : Comes in with the four spatial modes and exits with one time 
! mode of a second order accurate entropy function.
! gamma_zonecen, molewt_zonecen : zone-centered gamma and molewt.
! vmin, vmax : Min and max values of the primitive variables from all the
! neighboring zones. Done outside.
! corrupted_zone : If the zone has taken any density or pressure positivity
! fix-ups, it is set to 1; otherwise it is 0.

! The fundamental idea is that we construct nodal values in space (not
! time). We then test them for positivity. If the positivity is not obtained,
! we reduce the spatially varying part of the reconstructed profile.
! This is done in two stages. First, the density variation is shrunk,
! if needed. Second, all the variations (except certain necessary components
! of the magnetic field) are shrunk consistent with obtaining
! positive pressures, if that is needed.

! In this version, a quadratic is produced for the MHD equations and solved
! directly. The magnetic field modes are not allowed to vary with the
! result that the pressure positivity enforcement turns into a quadratic.
! The quadratic is much easier to solve than a cubic. Besides, for pure
! Euler flow, we don't even need to solve a cubic.

! This is one of ONLY two places where we may perhaps trade
! conservation in return for physical consistency. The other place is
! in the "EVALUATE_TIMESTEP" routine.

! ----------------------------------------------------------------------

      IMPLICIT NONE

! These variables go through the subroutine interface.

      INTEGER oned_begin, oned_end, oned_off, igeom

      REAL dtmove, pi, rhofloor, prsfloor, smallnum, gamma,
     1     gasconst

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1      xb_use, yb_use, yb_use_m1, zb_use,
     1      flatten, gamma_zonecen, molewt_zonecen

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR) :: vmin, vmax

! ----------------------------------------------------------------------

! Variables that are internal to the subroutine.

      INTEGER ix, iquadpt, imode, ivar, num_pp_quadpts

      REAL fourpi, ob4pi, ob8pi, ob_root_4pi, expandfac,
     1     del_rho, del_eng, del_mx, del_my, del_mz,
     1     del_bx, del_by, del_bz, m_bar2, del_m2, m_dot_del_m,
     1     ax1, ax2, ay1, ay2, az1, az2, ob_dx, ob_dy, ob_dz,
     1     tempaa, tempab, tempac, tempad

! Since we are interested in extremal variations within a zone, we need
! to define nodal sets that have more quadrature points than the smallest
! symmetrical set of quadrature points.
#if ( NUMORDER == 2)
      PARAMETER ( num_pp_quadpts = 8)
#elif ( NUMORDER == 3)
      PARAMETER ( num_pp_quadpts = 21)
#elif ( NUMORDER == 4)
      PARAMETER ( num_pp_quadpts = 33)
#endif

! These are the flags for flagging various situations and coefficients
! for curing various pathologies.

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: pp_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: reduce_modes

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                     num_pp_quadpts) :: pp_quadpt

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                     num_pp_quadpts) :: reduce_quadpt,
     1                     a0, a1, a2

! These are the bounds on the primitive variables from neighboring zones
! with some inbuilt expansions.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR) :: vminuse, vmaxuse, tempvec

! These are the flow variables at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  num_pp_quadpts) ::

     1      rho, prs, eng_t, ent, vlx, vly, vlz, bfx, bfy, bfz, bprs,
     1      gamma_eos_quadpt, molewt_quadpt, x_root

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0: NFLUID, num_pp_quadpts) :: specie

! These are the conserved variables at the quadrature points.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, num_pp_quadpts) ::

     1      u_quadpt, tu

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob4pi = 1.0 / fourpi
      ob8pi = 0.5 * ob4pi
      ob_root_4pi = 1.0 / SQRT ( fourpi)

! Set the factor by which we want to expand the range of each variable.
! Setting this to 0 gives us the TVD range in multi-dimensions.
      expandfac = 0.4
! ----------------------------------------------------------------------

! STEP 0) Start by saying that none of the zones are corrupted, none have
! negative pressures and that none of them need to have their moments reduced.
! Make the same assertion at all the quadrature points within each zone.
! By the end of this subroutine, most zones should continue to be that way.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      corrupted_zone ( ix) = 0
      pp_zone ( ix) = 0
      reduce_modes ( ix) = 1.0

      pp_quadpt ( ix, :) = 0
      reduce_quadpt ( ix, :) = 1.0

      gamma_eos_quadpt ( ix, :) = gamma_zonecen ( ix)
      molewt_quadpt ( ix, :) = molewt_zonecen ( ix)

      END DO

! ----------------------------------------------------------------------

! STEP 1) Set "corrupted_zone" to 0 and expand "vminuse, vmaxuse" if possible.
! The bounds on variables that should remain positive are set
! multiplicatively. The bounds on variables that aren't positive definite
! are set based on expanding the range.

      DO ix = oned_begin - oned_off, oned_end + oned_off

      corrupted_zone ( ix) = 0.0

! Expand density range:
      vminuse ( ix, 1) = vmin ( ix, 1)
     1      * ( 1.0 - expandfac + expandfac * flatten ( ix) )
      vminuse ( ix, 1) = AMAX1 ( vminuse ( ix, 1), rhofloor)
      vmaxuse ( ix, 1) = vmax ( ix, 1)
     1      * ( 1.0 + expandfac - expandfac * flatten ( ix) )

! Expand pressure range:
      vminuse ( ix, 5) = vmin ( ix, 5)
     1      * ( 1.0 - expandfac + expandfac * flatten ( ix) )
      vminuse ( ix, 5) = AMAX1 ( vminuse ( ix, 5), prsfloor)
      vmaxuse ( ix, 5) = vmax ( ix, 5)
     1      * ( 1.0 + expandfac - expandfac * flatten ( ix) )

#if ( REDUCE_VEL_BF == 1)
! Expand velocity range:
      tempvec ( ix, 2:4) = expandfac * ( 1.0 - flatten ( ix) )
     1       * ( vmax ( ix, 2:4) - vmin ( ix, 2:4) )
      vminuse ( ix, 2:4) = vmin ( ix, 2:4) - tempvec ( ix, 2:4)
      vmaxuse ( ix, 2:4) = vmax ( ix, 2:4) + tempvec ( ix, 2:4)

! Expand magnetic field range:
      tempvec ( ix, 6:8) = expandfac * ( 1.0 - flatten ( ix) )
     1       * ( vmax ( ix, 6:8) - vmin ( ix, 6:8) )
      vminuse ( ix, 6:8) = vmin ( ix, 6:8) - tempvec ( ix, 6:8)
      vmaxuse ( ix, 6:8) = vmax ( ix, 6:8) + tempvec ( ix, 6:8)
#endif

#if ( NFLUID >= 1)
! Expand the specie fraction range:
      vminuse ( ix, 9:8 + NFLUID) = vmin ( ix, 9:8 + NFLUID)
     1      * ( 1.0 - expandfac + expandfac * flatten ( ix) )
      vminuse ( ix, 9:8 + NFLUID) =
     1                      AMAX1 ( vminuse ( ix, 9:8 + NFLUID), 0.0)
      vmaxuse ( ix, 9:8 + NFLUID) = vmax ( ix, 9:8 + NFLUID)
     1      * ( 1.0 + expandfac - expandfac * flatten ( ix) )
      vmaxuse ( ix, 9:8 + NFLUID) =
     1                      AMIN1 ( vmaxuse ( ix, 9:8 + NFLUID), 1.0)
#endif

      END DO ! End of "DO ix =".

! ----------------------------------------------------------------------

! STEP 2) Use only the spatial modes to obtain the conserved variables
! at the nodes of interest. I.e., use the spatial modes to make 
! "u_quadpt ( ix, :, 1: num_pp_quadpts)"
! from "w ( ix, :, 1: NUM_ZONECEN_BASES)" at time t = 0.

! No changes are made to the nodal magnetic fields.

! --------------------
#if ( NUMORDER == 2)
! cg, 2nd order : modal -> nodal ( at zone vertices).

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 02) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 03) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 04) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) + 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 05) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 06) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                + 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 07) = w ( ix, :, 01) + 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

      u_quadpt ( ix, :, 08) = w ( ix, :, 01) - 0.5 * w ( ix, :, 02)
     1                - 0.5 * w ( ix, :, 03) - 0.5 * w ( ix, :, 04)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2"
! --------------------
#if ( NUMORDER == 3)
! cg, 3rd order : modal -> nodal ( including zone vertices).

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = w ( ix, :, 01) - ( w ( ix, :, 05)
     1 + w ( ix, :, 06) + w ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( w ( ix, :, 02) + w ( ix, :, 03) )
     1 + 0.25 * ( w ( ix, :, 05) + w ( ix, :, 06) + w ( ix, :, 08) )

      u_quadpt ( ix, :, 03) = u_quadpt ( ix, :, 02) - w ( ix, :, 03)
     1 - 0.5 * w ( ix, :, 08)

      u_quadpt ( ix, :, 04) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( w ( ix, :, 02) + w ( ix, :, 04) )
     1 + 0.25 * ( w ( ix, :, 05) + w ( ix, :, 07) + w ( ix, :, 10) )

      u_quadpt ( ix, :, 05) = u_quadpt ( ix, :, 04) - w ( ix, :, 04)
     1 - 0.5 * w ( ix, :, 10)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 03) - w ( ix, :, 02)
     1 + w ( ix, :, 03)

      u_quadpt ( ix, :, 07) = u_quadpt ( ix, :, 02) - w ( ix, :, 02)
     1 - w ( ix, :, 03)

      u_quadpt ( ix, :, 08) = u_quadpt ( ix, :, 05) - w ( ix, :, 02)
     1 + w ( ix, :, 04)

      u_quadpt ( ix, :, 09) = u_quadpt ( ix, :, 04) - w ( ix, :, 02)
     1 - w ( ix, :, 04)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + 0.5 * ( w ( ix, :, 03) + w ( ix, :, 04) )
     1 + 0.25 * ( w ( ix, :, 06) + w ( ix, :, 07) + w ( ix, :, 09) )

      u_quadpt ( ix, :, 11) = u_quadpt ( ix, :, 10) - w ( ix, :, 03)
     1 - 0.5 * w ( ix, :, 09)

      u_quadpt ( ix, :, 12) = u_quadpt ( ix, :, 11) + w ( ix, :, 03)
     1 - w ( ix, :, 04)

      u_quadpt ( ix, :, 13) = u_quadpt ( ix, :, 10) - w ( ix, :, 03)
     1 - w ( ix, :, 04)

! -----
! Additional nodes at vertices.

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 04) + u_quadpt ( ix, :, 10)
     1 - 2.0 * u_quadpt ( ix, :, 01)
     1 - 0.5 * ( w ( ix, :, 02) + w ( ix, :, 03) + w ( ix, :, 04) )
     1 - 0.25 * ( w ( ix, :, 05) + w ( ix, :, 06) + w ( ix, :, 07) )

      u_quadpt ( ix, :, 15) = 2.0 * u_quadpt ( ix, :, 10)
     1 - u_quadpt ( ix, :, 14) + 0.5 * w ( ix, :, 05)

      u_quadpt ( ix, :, 16) = 2.0 * u_quadpt ( ix, :, 04)
     1 - u_quadpt ( ix, :, 14) + 0.5 * w ( ix, :, 06)

      u_quadpt ( ix, :, 17) = 2.0 * u_quadpt ( ix, :, 11)
     1 - u_quadpt ( ix, :, 16) + 0.5 * w ( ix, :, 05)

      u_quadpt ( ix, :, 18) = 2.0 * u_quadpt ( ix, :, 02)
     1 - u_quadpt ( ix, :, 14) + 0.5 * w ( ix, :, 07)

      u_quadpt ( ix, :, 19) = 2.0 * u_quadpt ( ix, :, 12)
     1 - u_quadpt ( ix, :, 18) + 0.5 * w ( ix, :, 05)

      u_quadpt ( ix, :, 20) = 2.0 * u_quadpt ( ix, :, 05)
     1 - u_quadpt ( ix, :, 18) + 0.5 * w ( ix, :, 06)

      u_quadpt ( ix, :, 21) = 2.0 * u_quadpt ( ix, :, 13)
     1 - u_quadpt ( ix, :, 20) + 0.5 * w ( ix, :, 05)

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3"
! --------------------
#if ( NUMORDER == 4)
! cg, 4th order : modal -> nodal (zone vertices are already included).

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      u_quadpt ( ix, :, 01) = w ( ix, :, 01) - ( w ( ix, :, 05)
     1 + w ( ix, :, 06) + w ( ix, :, 07) ) / 12.0

      u_quadpt ( ix, :, 02) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * w ( ix, :, 02) + 30.0 * w ( ix, :, 05)
     1 + 6.0 * w ( ix, :, 11) - 5.0 * w ( ix, :, 16)
     1 - 5.0 * w ( ix, :, 18) ) / 120.0

      u_quadpt ( ix, :, 03) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) ) - 0.015625 * ( 4 * w ( ix, :, 05)
     1 + 3.0 * w ( ix, :, 11) )

      u_quadpt ( ix, :, 04) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 02) - 3.0 * u_quadpt ( ix, :, 03)
     1 - 0.09375 * w ( ix, :, 11)

      u_quadpt ( ix, :, 05) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 02) - 8.0 * u_quadpt ( ix, :, 03)
     1 - 0.375 * w ( ix, :, 11)

      u_quadpt ( ix, :, 06) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * w ( ix, :, 03) + 30.0 * w ( ix, :, 06)
     1 + 6.0 * w ( ix, :, 12) - 5.0 * w ( ix, :, 14)
     1 - 5.0 * w ( ix, :, 19) ) / 120.0

      u_quadpt ( ix, :, 07) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) ) - 0.015625 * ( 4.0 * w ( ix, :, 06)
     1 + 3.0 * w ( ix, :, 12) )

      u_quadpt ( ix, :, 08) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 06) - 3.0 * u_quadpt ( ix, :, 07)
     1 - 0.09375 * w ( ix, :, 12)

      u_quadpt ( ix, :, 09) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 06) - 8.0 * u_quadpt ( ix, :, 07)
     1 - 0.375 * w ( ix, :, 12)

      u_quadpt ( ix, :, 10) = u_quadpt ( ix, :, 01)
     1 + ( 60.0 * w ( ix, :, 04) + 30.0 * w ( ix, :, 07)
     1 + 6.0 * w ( ix, :, 13) - 5.0 * w ( ix, :, 15)
     1 - 5.0 * w ( ix, :, 17) ) / 120.0

      u_quadpt ( ix, :, 11) = 0.5 * ( u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) ) - 0.015625 * ( 4.0 * w ( ix, :, 07)
     1 + 3.0 * w ( ix, :, 13) )

      u_quadpt ( ix, :, 12) = 3.0 * u_quadpt ( ix, :, 01)
     1 + u_quadpt ( ix, :, 10) - 3.0 * u_quadpt ( ix, :, 11)
     1 - 0.09375 * w ( ix, :, 13)

      u_quadpt ( ix, :, 13) = 6.0 * u_quadpt ( ix, :, 01)
     1 + 3.0 * u_quadpt ( ix, :, 10) - 8.0 * u_quadpt ( ix, :, 11)
     1 - 0.375 * w ( ix, :, 13)

      u_quadpt ( ix, :, 14) = u_quadpt ( ix, :, 02)
     1 + u_quadpt ( ix, :, 07) + u_quadpt ( ix, :, 10)
     1 - 2.0 * u_quadpt ( ix, :, 01) + 0.25 * ( w ( ix, :, 03)
     1 + w ( ix, :, 08) + w ( ix, :, 09) + w ( ix, :, 10) )
     1 + 0.125 * ( w ( ix, :, 15) + w ( ix, :, 16) + w ( ix, :, 17)
     1 + w ( ix, :, 18) + w ( ix, :, 20) ) + 0.1875 * w ( ix, :, 06)
     1 + 5.0 * ( w ( ix, :, 14) + w ( ix, :, 19) ) / 48.0
     1 + 0.071875 * w ( ix, :, 12)

      u_quadpt ( ix, :, 15) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 - 0.5 * ( w ( ix, :, 08) + w ( ix, :, 09) )
     1 - 0.25 * ( w ( ix, :, 14) + w ( ix, :, 19) + w ( ix, :, 20) )

      u_quadpt ( ix, :, 16) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( w ( ix, :, 09) + w ( ix, :, 10) )
     1 - 0.25 * ( w ( ix, :, 15) + w ( ix, :, 17) + w ( ix, :, 20) )

      u_quadpt ( ix, :, 17) = u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 16) - u_quadpt ( ix, :, 14) + w ( ix, :, 09)
     1 + 0.5 * w ( ix, :, 20)

      u_quadpt ( ix, :, 18) = u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 05) - u_quadpt ( ix, :, 02)
     1 - 0.5 * ( w ( ix, :, 08) + w ( ix, :, 10) )
     1 - 0.25 * ( w ( ix, :, 16) + w ( ix, :, 18) + w ( ix, :, 20) )

      u_quadpt ( ix, :, 19) = u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 09) - u_quadpt ( ix, :, 06)
     1 + 0.5 * ( w ( ix, :, 08) - w ( ix, :, 09) )
     1 - 0.25 * ( w ( ix, :, 14) + w ( ix, :, 19) - w ( ix, :, 20) )

      u_quadpt ( ix, :, 20) = u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 13) - u_quadpt ( ix, :, 10)
     1 - 0.5 * ( w ( ix, :, 09) - w ( ix, :, 10) )
     1 - 0.25 * ( w ( ix, :, 15) + w ( ix, :, 17) - w ( ix, :, 20) )

      u_quadpt ( ix, :, 21) = u_quadpt ( ix, :, 19)
     1 + u_quadpt ( ix, :, 20) - u_quadpt ( ix, :, 18) + w ( ix, :, 09)
     1 - 0.5 * w ( ix, :, 20)

! -----
! Additional nodes at zone edges.

      u_quadpt ( ix, :, 22) = 0.5 * ( u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 16) ) - 0.25 * w ( ix, :, 07)
     1 - 0.125 * ( w ( ix, :, 18) + w ( ix, :, 19) )

      u_quadpt ( ix, :, 23) = 0.5 * ( u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 17) ) - 0.25 * w ( ix, :, 07)
     1 - 0.125 * ( w ( ix, :, 18) - w ( ix, :, 19) )

      u_quadpt ( ix, :, 24) = 0.5 * ( u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 15) ) - 0.25 * w ( ix, :, 06)
     1 - 0.125 * ( w ( ix, :, 16) + w ( ix, :, 17) )

      u_quadpt ( ix, :, 25) = 0.5 * ( u_quadpt ( ix, :, 16)
     1 + u_quadpt ( ix, :, 17) ) - 0.25 * w ( ix, :, 06)
     1 - 0.125 * ( w ( ix, :, 16) - w ( ix, :, 17) )

      u_quadpt ( ix, :, 26) = 0.5 * ( u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 20) ) - 0.25 * w ( ix, :, 07)
     1 + 0.125 * ( w ( ix, :, 18) - w ( ix, :, 19) )

      u_quadpt ( ix, :, 27) = 0.5 * ( u_quadpt ( ix, :, 19)
     1 + u_quadpt ( ix, :, 21) ) - 0.25 * w ( ix, :, 07)
     1 + 0.125 * ( w ( ix, :, 18) + w ( ix, :, 19) )

      u_quadpt ( ix, :, 28) = 0.5 * ( u_quadpt ( ix, :, 18)
     1 + u_quadpt ( ix, :, 19) ) - 0.25 * w ( ix, :, 06)
     1 + 0.125 * ( w ( ix, :, 16) - w ( ix, :, 17) )

      u_quadpt ( ix, :, 29) = 0.5 * ( u_quadpt ( ix, :, 20)
     1 + u_quadpt ( ix, :, 21) ) - 0.25 * w ( ix, :, 06)
     1 + 0.125 * ( w ( ix, :, 16) + w ( ix, :, 17) )

      u_quadpt ( ix, :, 30) = 0.5 * ( u_quadpt ( ix, :, 14)
     1 + u_quadpt ( ix, :, 18) ) - 0.25 * w ( ix, :, 05)
     1 - 0.125 * ( w ( ix, :, 14) + w ( ix, :, 15) )

      u_quadpt ( ix, :, 31) = 0.5 * ( u_quadpt ( ix, :, 15)
     1 + u_quadpt ( ix, :, 19) ) - 0.25 * w ( ix, :, 05)
     1 + 0.125 * ( w ( ix, :, 14) - w ( ix, :, 15) )

      u_quadpt ( ix, :, 32) = 0.5 * ( u_quadpt ( ix, :, 16)
     1 + u_quadpt ( ix, :, 20) ) - 0.25 * w ( ix, :, 05)
     1 - 0.125 * ( w ( ix, :, 14) - w ( ix, :, 15) )

      u_quadpt ( ix, :, 33) = 0.5 * ( u_quadpt ( ix, :, 17)
     1 + u_quadpt ( ix, :, 21) ) - 0.25 * w ( ix, :, 05)
     1 + 0.125 * ( w ( ix, :, 14) + w ( ix, :, 15) )

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4"
! ----------------------------------------------------------------------

! STEP 3) Get the primitive variables from the conserved variables,
! but only do this for density, velocities and magnetic fields; not
! pressure. These are easy to bring within bounds, so we do that first.

! Actually, it is best to focus on bringing only the density within bounds.

! Here we build the parts of the magnetic field that can be varied without
! losing consistency. These parts are stored in "bfx_var" etc.

! --------------------

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempaa
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempaa
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempaa

      bfx ( ix, iquadpt) = u_quadpt ( ix, 6, iquadpt)
      bfy ( ix, iquadpt) = u_quadpt ( ix, 7, iquadpt)
      bfz ( ix, iquadpt) = u_quadpt ( ix, 8, iquadpt)

      bprs ( ix, iquadpt) = ( bfx ( ix, iquadpt)**2
     1  + bfy ( ix, iquadpt)**2 + bfz ( ix, iquadpt)**2 ) * ob8pi

! ----------

! Bring density within bounds.

      IF ( rho ( ix, iquadpt) .GT. vmaxuse ( ix, 1) * 1.000001) THEN
      corrupted_zone ( ix) = 1
      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 1) - w ( ix, 1, 1) )
     1                   / ( rho ( ix, iquadpt) - w ( ix, 1, 1) ) )
      END IF

      IF ( rho ( ix, iquadpt) .LT. vminuse ( ix, 1) * 0.999999) THEN
      corrupted_zone ( ix) = 1
      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 1) - w ( ix, 1, 1) )
     1                   / ( rho ( ix, iquadpt) - w ( ix, 1, 1) ) )
      END IF

! -----
#if ( REDUCE_VEL_BF == 1)

! Bring x-velocity within bounds.

      IF ( vlx ( ix, iquadpt) .GT. vmaxuse ( ix, 2) + smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 2, 1) - w ( ix, 1, 1) * vmaxuse ( ix, 2)
      tempac = u_quadpt ( ix, 2, iquadpt)
     1       - rho ( ix, iquadpt) * vmaxuse ( ix, 2)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

      IF ( vlx ( ix, iquadpt) .LT. vminuse ( ix, 2) - smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 2, 1) - w ( ix, 1, 1) * vminuse ( ix, 2)
      tempac = u_quadpt ( ix, 2, iquadpt)
     1       - rho ( ix, iquadpt) * vminuse ( ix, 2)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

! -----

! Bring y-velocity within bounds.

      IF ( vly ( ix, iquadpt) .GT. vmaxuse ( ix, 3) + smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 3, 1) - w ( ix, 1, 1) * vmaxuse ( ix, 3)
      tempac = u_quadpt ( ix, 3, iquadpt)
     1       - rho ( ix, iquadpt) * vmaxuse ( ix, 3)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

      IF ( vly ( ix, iquadpt) .LT. vminuse ( ix, 3) - smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 3, 1) - w ( ix, 1, 1) * vminuse ( ix, 3)
      tempac = u_quadpt ( ix, 3, iquadpt)
     1       - rho ( ix, iquadpt) * vminuse ( ix, 3)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

! -----

! Bring z-velocity within bounds.

      IF ( vlz ( ix, iquadpt) .GT. vmaxuse ( ix, 4) + smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 4, 1) - w ( ix, 1, 1) * vmaxuse ( ix, 4)
      tempac = u_quadpt ( ix, 4, iquadpt)
     1       - rho ( ix, iquadpt) * vmaxuse ( ix, 4)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

      IF ( vlz ( ix, iquadpt) .LT. vminuse ( ix, 4) - smallnum) THEN

      corrupted_zone ( ix) = 1

      tempab = w ( ix, 4, 1) - w ( ix, 1, 1) * vminuse ( ix, 4)
      tempac = u_quadpt ( ix, 4, iquadpt)
     1       - rho ( ix, iquadpt) * vminuse ( ix, 4)
      tempac = tempab - tempac

        IF ( ABS ( tempac) .GT. smallnum) THEN
          reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                                 tempab / tempac)
        ELSE
          reduce_modes ( ix) = 0.0
        END IF

      END IF

! -----

! Bring Bx within bounds.

      IF ( bfx ( ix, iquadpt) .GT. vmaxuse ( ix, 6) + smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 6) - w ( ix, 6, 1) )
     1                   / ( bfx ( ix, iquadpt) - w ( ix, 6, 1) ) )

      END IF

      IF ( bfx ( ix, iquadpt) .LT. vminuse ( ix, 6) - smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 6) - w ( ix, 6, 1) )
     1                   / ( bfx ( ix, iquadpt) - w ( ix, 6, 1) ) )

      END IF

! -----

! Bring By within bounds.

      IF ( bfy ( ix, iquadpt) .GT. vmaxuse ( ix, 7) + smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 7) - w ( ix, 7, 1) )
     1                   / ( bfy ( ix, iquadpt) - w ( ix, 7, 1) ) )

      END IF

      IF ( bfy ( ix, iquadpt) .LT. vminuse ( ix, 7) - smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 7) - w ( ix, 7, 1) )
     1                   / ( bfy ( ix, iquadpt) - w ( ix, 7, 1) ) )

      END IF

! -----

! Bring Bz within bounds.

      IF ( bfz ( ix, iquadpt) .GT. vmaxuse ( ix, 8) + smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vmaxuse ( ix, 8) - w ( ix, 8, 1) )
     1                   / ( bfz ( ix, iquadpt) - w ( ix, 8, 1) ) )

      END IF

      IF ( bfz ( ix, iquadpt) .LT. vminuse ( ix, 8) - smallnum) THEN

      corrupted_zone ( ix) = 1

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                     ( vminuse ( ix, 8) - w ( ix, 8, 1) )
     1                   / ( bfz ( ix, iquadpt) - w ( ix, 8, 1) ) )

      END IF

#endif
! End of "REDUCE_VEL_BF == 1"
! -----

! Lastly, just make sure "reduce_modes" is inside [0,1].

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix), 1.0)
      reduce_modes ( ix) = AMAX1 ( reduce_modes ( ix), 0.0)

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 4) Now if we have found a suitable "reduce_modes", we reduce ALL the
! modal variations by that amount. We also reduce the variation in the 
! conserved variables by that amount. This leaves us free to rebuild the 
! nodal density, velocities and magnetic fields that are within bounds. 
! All that remains is to consider the pressure. 

! The reduction in velocities and fields also
! gives the pressure variable a second chance to become positive without
! further touch up. We will, eventually, touch it up if it is still
! non-positive.

! --------------------

      DO imode = 2, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      w ( ix, 1: 5, imode) = reduce_modes ( ix) * w ( ix, 1: 5, imode)

#if ( NFLUID >= 1)
      w ( ix, 9: 8 + NFLUID, imode) = reduce_modes ( ix)
     1                              * w ( ix, 9: 8 + NFLUID, imode)
#endif

! -----

      END IF ! End of "corrupted_zone ( ix) .EQ. 1".

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO imode =" loop.

! --------------------

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( corrupted_zone ( ix) .EQ. 1) THEN

! -----

      tempaa = 1.0 - reduce_modes ( ix)

      u_quadpt ( ix, 1: 5, iquadpt) = tempaa * w ( ix, 1: 5, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 1: 5, iquadpt)

#if ( NFLUID >= 1)
      u_quadpt ( ix, 9: 8 + NFLUID, iquadpt) =
     1                       tempaa * w ( ix, 9: 8 + NFLUID, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 9: 8 + NFLUID, iquadpt)
#endif

! -----

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempab = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempab
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempab
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempab

! -----

      END IF ! End of "corrupted_zone ( ix) .EQ. 1".

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 5) The zones with "corrupted_zone ( ix) == 1" now have taken on
! corrections, the rest have not. Build pressure for all zones here.
! Flag "pp_zone" for zones with pressures that go below floor values.
! Solve the quadratic to decide on the factor by which we will shrink
! the variation in the modes.

! Build the coefficients of the quadratic and obtain the roots that lie
! within [0,1].

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

! Now build pressure and species fractions and flag if needed.

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)

      tempaa = 1.0 / rho ( ix, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      prs ( ix, iquadpt) = ( u_quadpt ( ix, 5, iquadpt) - tempab
     1         - tempac) * ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

#if ( NFLUID >= 1)
      specie ( ix, 1: NFLUID, iquadpt) = tempaa
     1       * u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt)
      specie ( ix, 1: NFLUID, iquadpt) = AMIN1 ( AMAX1 (
     1                    specie ( ix, 1: NFLUID, iquadpt), 0.0), 1.0 )
      u_quadpt ( ix, NFIELD + 2: NFIELD + 1 + NFLUID, iquadpt) =
     1       specie ( ix, 1: NFLUID, iquadpt) * rho ( ix, iquadpt)
#endif

! ----------

      IF ( prs ( ix, iquadpt) .LT. vminuse ( ix, 5) ) THEN

      corrupted_zone ( ix) = 1
      pp_zone ( ix) = 1
      pp_quadpt ( ix, iquadpt) = 1

      del_rho = rho ( ix, iquadpt) - w ( ix, 1, 1)
      del_eng = eng_t ( ix, iquadpt) - w ( ix, 5, 1)

      del_mx = u_quadpt ( ix, 2, iquadpt) - w ( ix, 2, 1)
      del_my = u_quadpt ( ix, 3, iquadpt) - w ( ix, 3, 1)
      del_mz = u_quadpt ( ix, 4, iquadpt) - w ( ix, 4, 1)

      m_bar2 = w ( ix, 2, 1)**2 + w ( ix, 3, 1)**2 + w ( ix, 4, 1)**2
      del_m2 = del_mx**2 + del_my**2 + del_mz**2
      m_dot_del_m = w ( ix, 2, 1) * del_mx + w ( ix, 3, 1) * del_my
     1            + w ( ix, 4, 1) * del_mz

      tempad = bprs ( ix, iquadpt) + vminuse ( ix, 5)
     1       / ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

      a2 ( ix, iquadpt) = 2.0 * del_rho * del_eng - del_m2

      a1 ( ix, iquadpt) = 2.0 * w ( ix, 1, 1) * del_eng
     1                  + 2.0 * w ( ix, 5, 1) * del_rho
     1   - 2.0 * m_dot_del_m - 2.0 * tempad * del_rho

      a0 ( ix, iquadpt) = 2.0 * w ( ix, 1, 1) * w ( ix, 5, 1) - m_bar2
     1   - 2.0 * tempad * w ( ix, 1, 1)

! Evaluate discriminant.
      tempaa = a1 ( ix, iquadpt)**2
     1       - 4.0 * a2 ( ix, iquadpt) * a0 ( ix, iquadpt)

        IF ( ABS ( a2 ( ix, iquadpt)) .LT. smallnum) THEN
! If the quadratic reduces to a linear equation.

          IF ( ABS ( a1 ( ix, iquadpt)) .GT. smallnum) THEN
            tempab = - a0 ( ix, iquadpt) / a1 ( ix, iquadpt)
          ELSE
            tempab = 2.0
          END IF

          tempac = 2.0

        ELSE IF ( tempaa .GE. 0.0) THEN
! If the quadratic has two real roots, store them in "tempab" and "tempac".

          tempab = 0.5 * ( - a1 ( ix, iquadpt) + SQRT ( tempaa) )
     1           / a2 ( ix, iquadpt)
          tempac = 0.5 * ( - a1 ( ix, iquadpt) - SQRT ( tempaa) )
     1           / a2 ( ix, iquadpt)

        ELSE
! Quadratic with complex roots, useless for our purposes.Set to absurd values.

          tempab = 2.0
          tempac = 2.0

        END IF

! Now obtain "x_root ( ix, iquadpt)" by polling the one or two roots.

        x_root ( ix, iquadpt) = 2.0

        IF ( (0.0 .LE. tempab) .AND. (tempab .LE. 1.0) )
     1  x_root ( ix, iquadpt) = AMIN1 ( x_root ( ix, iquadpt), tempab)

        IF ( (0.0 .LE. tempac) .AND. (tempac .LE. 1.0) )
     1  x_root ( ix, iquadpt) = AMIN1 ( x_root ( ix, iquadpt), tempac)

        IF ( x_root ( ix, iquadpt) .GT. 1.0)x_root ( ix, iquadpt) = 0.0

      END IF ! End of "prs ( ix, iquadpt) .LT. vminuse ( ix, 5)".

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 6) By this point we know how much further we need to shrink the modes.
! Poll all the nodes of a flagged zone to find the smallest value of
! "x_root ( ix, iquadpt)" over all the quadrature points in that zone.

! To begin with, reset "reduce_modes". Then minimize "reduce_modes".

      DO ix = oned_begin - oned_off, oned_end + oned_off
      reduce_modes ( ix) = 1.0
      END DO


      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_quadpt ( ix, iquadpt) .EQ. 1) THEN

      reduce_modes ( ix) = AMIN1 ( reduce_modes ( ix),
     1                             0.95 * x_root ( ix, iquadpt) )

      END IF ! End of "pp_quadpt ( ix, iquadpt) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 7) This is almost a repeat of STEP 4, where we reduce all the
! possible modes that can be reduced by a factor of "reduce_modes".
! This time we key off "pp_zone".

! --------------------

      DO imode = 2, NUM_ZONECEN_BASES

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( pp_zone ( ix) .EQ. 1) THEN

! -----

      w ( ix, 1: 5, imode) = reduce_modes ( ix) * w ( ix, 1: 5, imode)

#if ( NFLUID >= 1)
      w ( ix, 9: 8 + NFLUID, imode) = reduce_modes ( ix)
     1                              * w ( ix, 9: 8 + NFLUID, imode)
#endif

! -----

      END IF

      END DO
! End of "DO ix =" loop.

      END DO
! End of "DO imode =" loop.

! ----------------------------------------------------------------------

! STEP 8) Now make sure that the nodal values in the zones that originally
! lost pressure positivity at least have a pressure that exceeds
! "vminuse ( ix, 5)".

! Rebuild the total energy density at those nodes, i.e. rebuild
! "u_quadpt ( ix, 5, iquadpt)" for the flagged nodes.

      DO iquadpt = 1, num_pp_quadpts

      DO ix = oned_begin - oned_off, oned_end + oned_off

      IF ( pp_zone ( ix) .EQ. 1) THEN

! -----

      tempaa = 1.0 - reduce_modes ( ix)

      u_quadpt ( ix, 1: 5, iquadpt) = tempaa * w ( ix, 1: 5, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 1: 5, iquadpt)

#if ( NFLUID >= 1)
      u_quadpt ( ix, 9: 8 + NFLUID, iquadpt) =
     1                       tempaa * w ( ix, 9: 8 + NFLUID, 1)
     1  + reduce_modes ( ix) * u_quadpt ( ix, 9: 8 + NFLUID, iquadpt)
#endif

! -----

      rho ( ix, iquadpt) = u_quadpt ( ix, 1, iquadpt)

      tempab = 1.0 / rho ( ix, iquadpt)

      vlx ( ix, iquadpt) = u_quadpt ( ix, 2, iquadpt) * tempab
      vly ( ix, iquadpt) = u_quadpt ( ix, 3, iquadpt) * tempab
      vlz ( ix, iquadpt) = u_quadpt ( ix, 4, iquadpt) * tempab

! -----

      eng_t ( ix, iquadpt) = u_quadpt ( ix, 5, iquadpt)

      tempab = 0.5 * rho ( ix, iquadpt) * ( vlx ( ix, iquadpt)**2
     1         + vly ( ix, iquadpt)**2 + vlz ( ix, iquadpt)**2)

      tempac = ( bfx ( ix, iquadpt)**2 + bfy ( ix, iquadpt)**2
     1         + bfz ( ix, iquadpt)**2) * ob8pi

      prs ( ix, iquadpt) = ( u_quadpt ( ix, 5, iquadpt) - tempab
     1         - tempac) * ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

! This ensuing line is the only place where we might cause a loss of
! total energy conservation.
      prs ( ix, iquadpt)= AMAX1 ( prs ( ix, iquadpt), vminuse ( ix, 5))

      u_quadpt ( ix, 5, iquadpt) = tempab + tempac
     1  + prs ( ix, iquadpt) / ( gamma_eos_quadpt ( ix, iquadpt) - 1.0)

! -----

      END IF ! End of "IF ( pp_zone ( ix) .EQ. 1)".

      END DO
! End of "DO ix =" loop.

      END DO
! Finished "DO iquadpt ="

! ----------------------------------------------------------------------

! STEP 9) Now use the nodal values of just the total energy density in
! the affected zone to rebuild ALL the moments of the energy density
! (and just the energy density) in that zone.

! Notice that we also rebuild the zeroth moment.
! If "AMAX1 ( prs ( ix, iquadpt), vminuse ( ix, 5))" had caused any changes
! in the pressure in the previous step then we will also be changing
! the total energy density in the zone. In that case, we choose positivity
! over and above total energy conservation.

! We also build the entropy variable and its first moments for all the zones.

! --------------------------------------------------

      ivar = 5

! --------------------------------------------------
#if ( NUMORDER == 2)
! cg, 2nd order : nodal -> modal.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_zone ( ix) .EQ. 1) THEN

      w ( ix, ivar, 01) = 0.125 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 02) + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 04) + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 06) + u_quadpt ( ix, ivar, 07)
     1 + u_quadpt ( ix, ivar, 08) )

      w ( ix, ivar, 02) =  0.25 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 03) + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 07) - u_quadpt ( ix, ivar, 02)
     1 - u_quadpt ( ix, ivar, 04) - u_quadpt ( ix, ivar, 06)
     1 - u_quadpt ( ix, ivar, 08) )

      w ( ix, ivar, 03) =  0.25 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 02) + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 06) - u_quadpt ( ix, ivar, 03)
     1 - u_quadpt ( ix, ivar, 04) - u_quadpt ( ix, ivar, 07)
     1 - u_quadpt ( ix, ivar, 08) )

      w ( ix, ivar, 04) =  0.25 * ( u_quadpt ( ix, ivar, 01)
     1 + u_quadpt ( ix, ivar, 02) + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 04) - u_quadpt ( ix, ivar, 05)
     1 - u_quadpt ( ix, ivar, 06) - u_quadpt ( ix, ivar, 07)
     1 - u_quadpt ( ix, ivar, 08) )

      END IF ! End of "pp_zone ( ix) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ent ( ix, 01) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 02) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 03) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 04) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 05) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 04)

      ent ( ix, 06) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 04)

! ----------

      IF ( igeom .EQ. 1) THEN
      ax1 = 1.0
      ax2 = 1.0
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove / dx_use ( ix)
      ob_dy = dtmove / dy_use ( ix)
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 2) THEN
      ax1 = xb_use ( ix)
      ax2 = xb_use ( ix - 1)
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 2.0 / ( xb_use ( ix)**2 - xb_use ( ix - 1)**2)
      ob_dy = dtmove / ( xc_use ( ix) * dy_use ( ix))
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 3) THEN
      ax1 = xb_use ( ix)**2
      ax2 = xb_use ( ix - 1)**2
      ay1 = SIN ( yb_use ( ix))
      ay2 = SIN ( yb_use_m1 ( ix))
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 3.0 / ( xb_use ( ix)**3 - xb_use ( ix - 1)**3)
      ob_dy = dtmove / ( COS ( yb_use_m1 ( ix)) - COS ( yb_use ( ix)) )
      ob_dz = dtmove / ( xc_use ( ix) * SIN ( yc_use ( ix))
     1                                * dz_use ( ix))
      END IF

      entropy_1d ( ix, 05) = - 0.25 * ax1 * ob_dx * ent ( ix, 01)
     1*( vlx ( ix, 01) + vlx ( ix, 03) + vlx ( ix, 05) + vlx ( ix, 07))
     1                       + 0.25 * ax2 * ob_dx * ent ( ix, 02)
     1*( vlx ( ix, 02) + vlx ( ix, 04) + vlx ( ix, 06) + vlx ( ix, 08))

     1                       - 0.25 * ay1 * ob_dy * ent ( ix, 03)
     1*( vly ( ix, 01) + vly ( ix, 02) + vly ( ix, 05) + vly ( ix, 06))
     1                       + 0.25 * ay2 * ob_dy * ent ( ix, 04)
     1*( vly ( ix, 03) + vly ( ix, 04) + vly ( ix, 07) + vly ( ix, 08))

     1                       - 0.25 * az1 * ob_dz * ent ( ix, 05)
     1*( vlz ( ix, 01) + vlz ( ix, 02) + vlz ( ix, 03) + vlz ( ix, 04))
     1                       + 0.25 * az2 * ob_dz * ent ( ix, 06)
     1*( vlz ( ix, 05) + vlz ( ix, 06) + vlz ( ix, 07) + vlz ( ix, 08))

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 2".
! --------------------------------------------------
#if ( NUMORDER == 3)
! cg, 3rd order : nodal -> modal.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_zone ( ix) .EQ. 1) THEN

      tu ( ix, ivar, 01) = u_quadpt ( ix, ivar, 02)
     1 + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 04) + u_quadpt ( ix, ivar, 05)

      tu ( ix, ivar, 02) = u_quadpt ( ix, ivar, 06)
     1 + u_quadpt ( ix, ivar, 07)
     1 + u_quadpt ( ix, ivar, 08) + u_quadpt ( ix, ivar, 09)

      tu ( ix, ivar, 03) = u_quadpt ( ix, ivar, 10)
     1 + u_quadpt ( ix, ivar, 11)
     1 + u_quadpt ( ix, ivar, 12) + u_quadpt ( ix, ivar, 13)

      tu ( ix, ivar, 04) = u_quadpt ( ix, ivar, 02)
     1 + u_quadpt ( ix, ivar, 06)
     1 + u_quadpt ( ix, ivar, 10) + u_quadpt ( ix, ivar, 12)

      tu ( ix, ivar, 05) = u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 07)
     1 + u_quadpt ( ix, ivar, 11) + u_quadpt ( ix, ivar, 13)

      tu ( ix, ivar, 06) = u_quadpt ( ix, ivar, 04)
     1 + u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 08) + u_quadpt ( ix, ivar, 09)

      tu ( ix, ivar, 07) = u_quadpt ( ix, ivar, 04)
     1 + u_quadpt ( ix, ivar, 08)
     1 + u_quadpt ( ix, ivar, 10) + u_quadpt ( ix, ivar, 11)

      tu ( ix, ivar, 08) = u_quadpt ( ix, ivar, 05)
     1 + u_quadpt ( ix, ivar, 09)
     1 + u_quadpt ( ix, ivar, 12) + u_quadpt ( ix, ivar, 13)

      tu ( ix, ivar, 09) = u_quadpt ( ix, ivar, 02)
     1 + u_quadpt ( ix, ivar, 03)
     1 + u_quadpt ( ix, ivar, 06) + u_quadpt ( ix, ivar, 07)

      w ( ix, ivar, 02) = 0.25 * ( tu ( ix, ivar, 01)
     1 - tu ( ix, ivar, 02) )

      w ( ix, ivar, 03) = 0.25 * ( tu ( ix, ivar, 04)
     1 - tu ( ix, ivar, 05) )

      w ( ix, ivar, 04) = 0.25 * ( tu ( ix, ivar, 07)
     1 - tu ( ix, ivar, 08) )

      w ( ix, ivar, 05) = 0.5 * ( tu ( ix, ivar, 01)
     1 + tu ( ix, ivar, 02)
     1 - tu ( ix, ivar, 03) ) - 2.0 * u_quadpt ( ix, ivar, 01)

      w ( ix, ivar, 06) = 0.5 * ( tu ( ix, ivar, 04)
     1 + tu ( ix, ivar, 05)
     1 - tu ( ix, ivar, 06) ) - 2.0 * u_quadpt ( ix, ivar, 01)

      w ( ix, ivar, 07) = 0.5 * ( tu ( ix, ivar, 07)
     1 + tu ( ix, ivar, 08)
     1 - tu ( ix, ivar, 09) ) - 2.0 * u_quadpt ( ix, ivar, 01)

      w ( ix, ivar, 08) = u_quadpt ( ix, ivar, 02)
     1 - u_quadpt ( ix, ivar, 03)
     1 - u_quadpt ( ix, ivar, 06) + u_quadpt ( ix, ivar, 07)

      w ( ix, ivar, 09) = u_quadpt ( ix, ivar, 10)
     1 - u_quadpt ( ix, ivar, 11)
     1 - u_quadpt ( ix, ivar, 12) + u_quadpt ( ix, ivar, 13)

      w ( ix, ivar, 10) = u_quadpt ( ix, ivar, 04)
     1 - u_quadpt ( ix, ivar, 05)
     1 - u_quadpt ( ix, ivar, 08) + u_quadpt ( ix, ivar, 09)

      w ( ix, ivar, 01) = u_quadpt ( ix, ivar, 01)
     1 + ( w ( ix, ivar, 05)
     1 + w ( ix, ivar, 06) + w ( ix, ivar, 07) ) / 12.0

      END IF ! End of "pp_zone ( ix) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ent ( ix, 01) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 02) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 03) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 04) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 05) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 04)

      ent ( ix, 06) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 04)

! ----------

      IF ( igeom .EQ. 1) THEN
      ax1 = 1.0
      ax2 = 1.0
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove / dx_use ( ix)
      ob_dy = dtmove / dy_use ( ix)
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 2) THEN
      ax1 = xb_use ( ix)
      ax2 = xb_use ( ix - 1)
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 2.0 / ( xb_use ( ix)**2 - xb_use ( ix - 1)**2)
      ob_dy = dtmove / ( xc_use ( ix) * dy_use ( ix))
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 3) THEN
      ax1 = xb_use ( ix)**2
      ax2 = xb_use ( ix - 1)**2
      ay1 = SIN ( yb_use ( ix))
      ay2 = SIN ( yb_use_m1 ( ix))
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 3.0 / ( xb_use ( ix)**3 - xb_use ( ix - 1)**3)
      ob_dy = dtmove / ( COS ( yb_use_m1 ( ix)) - COS ( yb_use ( ix)) )
      ob_dz = dtmove / ( xc_use ( ix) * SIN ( yc_use ( ix))
     1                                * dz_use ( ix))
      END IF

      entropy_1d ( ix, 05) = - 0.25 * ax1 * ob_dx * ent ( ix, 01)
     1*( vlx ( ix, 14) + vlx ( ix, 16) + vlx ( ix, 18) + vlx ( ix, 20))
     1                       + 0.25 * ax2 * ob_dx * ent ( ix, 02)
     1*( vlx ( ix, 15) + vlx ( ix, 17) + vlx ( ix, 19) + vlx ( ix, 21))

     1                       - 0.25 * ay1 * ob_dy * ent ( ix, 03)
     1*( vly ( ix, 14) + vly ( ix, 15) + vly ( ix, 18) + vly ( ix, 19))
     1                       + 0.25 * ay2 * ob_dy * ent ( ix, 04)
     1*( vly ( ix, 16) + vly ( ix, 17) + vly ( ix, 20) + vly ( ix, 21))

     1                       - 0.25 * az1 * ob_dz * ent ( ix, 05)
     1*( vlz ( ix, 14) + vlz ( ix, 15) + vlz ( ix, 16) + vlz ( ix, 17))
     1                       + 0.25 * az2 * ob_dz * ent ( ix, 06)
     1*( vlz ( ix, 18) + vlz ( ix, 19) + vlz ( ix, 20) + vlz ( ix, 21))

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 3".
! --------------------------------------------------
#if ( NUMORDER == 4)
! cg, 4th order : nodal -> modal.

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      IF ( pp_zone ( ix) .EQ. 1) THEN

      w ( ix, ivar, 05) = 2.0 * ( u_quadpt ( ix, ivar, 02)
     1 - 2.0 * u_quadpt ( ix, ivar, 01) + u_quadpt ( ix, ivar, 05) )

      w ( ix, ivar, 06) = 2.0 * ( u_quadpt ( ix, ivar, 06)
     1 - 2.0 * u_quadpt ( ix, ivar, 01) + u_quadpt ( ix, ivar, 09) )

      w ( ix, ivar, 07) = 2.0 * ( u_quadpt ( ix, ivar, 10)
     1 - 2.0 * u_quadpt ( ix, ivar, 01) + u_quadpt ( ix, ivar, 13) )

      w ( ix, ivar, 08) = 0.5 * ( u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 16)
     1 - u_quadpt ( ix, ivar, 17) - u_quadpt ( ix, ivar, 18)
     1 + u_quadpt ( ix, ivar, 19) - u_quadpt ( ix, ivar, 20)
     1 + u_quadpt ( ix, ivar, 21) )

      w ( ix, ivar, 09) = 0.5 * ( u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 18)
     1 - u_quadpt ( ix, ivar, 20) - u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 17) - u_quadpt ( ix, ivar, 19)
     1 + u_quadpt ( ix, ivar, 21) )

      w ( ix, ivar, 10) = 0.5 * ( u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 15)
     1 - u_quadpt ( ix, ivar, 19) - u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 20) - u_quadpt ( ix, ivar, 17)
     1 + u_quadpt ( ix, ivar, 21) )

      w ( ix, ivar, 11) = 16 * ( u_quadpt ( ix, ivar, 02)
     1 - 2.0 * u_quadpt ( ix, ivar, 03)
     1 + 2.0 * u_quadpt ( ix, ivar, 04)
     1 - u_quadpt ( ix, ivar, 05) ) / 3.0

      w ( ix, ivar, 12) = 16.0 * ( u_quadpt ( ix, ivar, 06)
     1 - 2.0 * u_quadpt ( ix, ivar, 07)
     1 + 2.0 * u_quadpt ( ix, ivar, 08)
     1 - u_quadpt ( ix, ivar, 09) ) / 3.0

      w ( ix, ivar, 13) = 16.0 * ( u_quadpt ( ix, ivar, 10)
     1 - 2.0 * u_quadpt ( ix, ivar, 11)
     1 + 2.0 * u_quadpt ( ix, ivar, 12)
     1 - u_quadpt ( ix, ivar, 13) ) / 3.0

      tu ( ix, ivar, 01) = u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 17)

      tu ( ix, ivar, 02) = u_quadpt ( ix, ivar, 18)
     1 + u_quadpt ( ix, ivar, 19)
     1 + u_quadpt ( ix, ivar, 20) + u_quadpt ( ix, ivar, 21)

      tu ( ix, ivar, 03) = u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 20)

      tu ( ix, ivar, 04) = u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 17)
     1 + u_quadpt ( ix, ivar, 19) + u_quadpt ( ix, ivar, 21)

      tu ( ix, ivar, 05) = u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 15)
     1 + u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 19)

      tu ( ix, ivar, 06) = u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 17)
     1 + u_quadpt ( ix, ivar, 20) + u_quadpt ( ix, ivar, 21)

      w ( ix, ivar, 02) = 2.0 * ( u_quadpt ( ix, ivar, 02)
     1 - u_quadpt ( ix, ivar, 05) ) / 3.0 - 0.1 * w ( ix, ivar, 11)
     1 + 0.25 * ( tu ( ix, ivar, 01) - tu ( ix, ivar, 02) ) / 3.0

      w ( ix, ivar, 03) = 2.0 * ( u_quadpt ( ix, ivar, 06)
     1 - u_quadpt ( ix, ivar, 09) ) / 3.0 - 0.1 * w ( ix, ivar, 12)
     1 + 0.25 * ( tu ( ix, ivar, 03) - tu ( ix, ivar, 04) ) / 3.0

      w ( ix, ivar, 04) = 2.0 * ( u_quadpt ( ix, ivar, 10)
     1 - u_quadpt ( ix, ivar, 13) ) / 3.0 - 0.1 * w ( ix, ivar, 13)
     1 + 0.25 * ( tu ( ix, ivar, 05) - tu ( ix, ivar, 06) ) / 3.0

      w ( ix, ivar, 20) = u_quadpt ( ix, ivar, 14)
     1 - u_quadpt ( ix, ivar, 18)
     1 - u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 20)
     1 - u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 19)
     1 + u_quadpt ( ix, ivar, 17) - u_quadpt ( ix, ivar, 21)

      w ( ix, ivar, 01) = 0.125 * ( u_quadpt ( ix, ivar, 14)
     1 + u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 16)
     1 + u_quadpt ( ix, ivar, 17) + u_quadpt ( ix, ivar, 18)
     1 + u_quadpt ( ix, ivar, 19) + u_quadpt ( ix, ivar, 20)
     1 + u_quadpt ( ix, ivar, 21) ) - 0.5 * ( w ( ix, ivar, 05)
     1 + w ( ix, ivar, 06) + w ( ix, ivar, 07) ) / 3.0


      tu ( ix, ivar, 22) = 2 * u_quadpt ( ix, ivar, 22)
     1 - ( u_quadpt ( ix, ivar, 14) + u_quadpt ( ix, ivar, 16) )

      tu ( ix, ivar, 23) = 2 * u_quadpt ( ix, ivar, 23)
     1 - ( u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 17) )

      tu ( ix, ivar, 24) = 2 * u_quadpt ( ix, ivar, 24)
     1 - ( u_quadpt ( ix, ivar, 14) + u_quadpt ( ix, ivar, 15) )

      tu ( ix, ivar, 25) = 2 * u_quadpt ( ix, ivar, 25)
     1 - ( u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 17) )

      tu ( ix, ivar, 26) = 2 * u_quadpt ( ix, ivar, 26)
     1 - ( u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 20) )

      tu ( ix, ivar, 27) = 2 * u_quadpt ( ix, ivar, 27)
     1 - ( u_quadpt ( ix, ivar, 19) + u_quadpt ( ix, ivar, 21) )

      tu ( ix, ivar, 28) = 2 * u_quadpt ( ix, ivar, 28)
     1 - ( u_quadpt ( ix, ivar, 18) + u_quadpt ( ix, ivar, 19) )

      tu ( ix, ivar, 29) = 2 * u_quadpt ( ix, ivar, 29)
     1 - ( u_quadpt ( ix, ivar, 20) + u_quadpt ( ix, ivar, 21) )

      tu ( ix, ivar, 30) = 2 * u_quadpt ( ix, ivar, 30)
     1 - ( u_quadpt ( ix, ivar, 14) + u_quadpt ( ix, ivar, 18) )

      tu ( ix, ivar, 31) = 2 * u_quadpt ( ix, ivar, 31)
     1 - ( u_quadpt ( ix, ivar, 15) + u_quadpt ( ix, ivar, 19) )

      tu ( ix, ivar, 32) = 2 * u_quadpt ( ix, ivar, 32)
     1 - ( u_quadpt ( ix, ivar, 16) + u_quadpt ( ix, ivar, 20) )

      tu ( ix, ivar, 33) = 2 * u_quadpt ( ix, ivar, 33)
     1 - ( u_quadpt ( ix, ivar, 17) + u_quadpt ( ix, ivar, 21) )


      w ( ix, ivar, 14) = tu ( ix, ivar, 31) + tu ( ix, ivar, 33)
     1 - tu ( ix, ivar, 30) - tu ( ix, ivar, 32)

      w ( ix, ivar, 15) = tu ( ix, ivar, 32) + tu ( ix, ivar, 33)
     1 - tu ( ix, ivar, 30) - tu ( ix, ivar, 31)

      w ( ix, ivar, 16) = tu ( ix, ivar, 28) + tu ( ix, ivar, 29)
     1 - tu ( ix, ivar, 24) - tu ( ix, ivar, 25)

      w ( ix, ivar, 17) = tu ( ix, ivar, 25) + tu ( ix, ivar, 29)
     1 - tu ( ix, ivar, 24) - tu ( ix, ivar, 28)

      w ( ix, ivar, 18) = tu ( ix, ivar, 26) + tu ( ix, ivar, 27)
     1 - tu ( ix, ivar, 22) - tu ( ix, ivar, 23)

      w ( ix, ivar, 19) = tu ( ix, ivar, 27) - tu ( ix, ivar, 26)
     1 + tu ( ix, ivar, 23) - tu ( ix, ivar, 22)


      END IF ! End of "pp_zone ( ix) .EQ. 1".

! ----------

      END DO
! End of "DO ix =" loop.

! ------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

! ----------

      ent ( ix, 01) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 02) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 02)

      ent ( ix, 03) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 04) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 03)

      ent ( ix, 05) = entropy_1d ( ix, 01) + 0.5 * entropy_1d ( ix, 04)

      ent ( ix, 06) = entropy_1d ( ix, 01) - 0.5 * entropy_1d ( ix, 04)

! ----------

      IF ( igeom .EQ. 1) THEN
      ax1 = 1.0
      ax2 = 1.0
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove / dx_use ( ix)
      ob_dy = dtmove / dy_use ( ix)
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 2) THEN
      ax1 = xb_use ( ix)
      ax2 = xb_use ( ix - 1)
      ay1 = 1.0
      ay2 = 1.0
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 2.0 / ( xb_use ( ix)**2 - xb_use ( ix - 1)**2)
      ob_dy = dtmove / ( xc_use ( ix) * dy_use ( ix))
      ob_dz = dtmove / dz_use ( ix)
      END IF

      IF ( igeom .EQ. 3) THEN
      ax1 = xb_use ( ix)**2
      ax2 = xb_use ( ix - 1)**2
      ay1 = SIN ( yb_use ( ix))
      ay2 = SIN ( yb_use_m1 ( ix))
      az1 = 1.0
      az2 = 1.0
      ob_dx = dtmove * 3.0 / ( xb_use ( ix)**3 - xb_use ( ix - 1)**3)
      ob_dy = dtmove / ( COS ( yb_use_m1 ( ix)) - COS ( yb_use ( ix)) )
      ob_dz = dtmove / ( xc_use ( ix) * SIN ( yc_use ( ix))
     1                                * dz_use ( ix))
      END IF

      entropy_1d ( ix, 05) = - 0.25 * ax1 * ob_dx * ent ( ix, 01)
     1*( vlx ( ix, 14) + vlx ( ix, 15) + vlx ( ix, 16) + vlx ( ix, 17))
     1                       + 0.25 * ax2 * ob_dx * ent ( ix, 02)
     1*( vlx ( ix, 18) + vlx ( ix, 19) + vlx ( ix, 20) + vlx ( ix, 21))

     1                       - 0.25 * ay1 * ob_dy * ent ( ix, 03)
     1*( vly ( ix, 14) + vly ( ix, 18) + vly ( ix, 16) + vly ( ix, 20))
     1                       + 0.25 * ay2 * ob_dy * ent ( ix, 04)
     1*( vly ( ix, 15) + vly ( ix, 19) + vly ( ix, 17) + vly ( ix, 21))

     1                       - 0.25 * az1 * ob_dz * ent ( ix, 05)
     1*( vlz ( ix, 14) + vlz ( ix, 18) + vlz ( ix, 15) + vlz ( ix, 19))
     1                       + 0.25 * az2 * ob_dz * ent ( ix, 06)
     1*( vlz ( ix, 16) + vlz ( ix, 20) + vlz ( ix, 17) + vlz ( ix, 21))

! ----------

      END DO
! End of "DO ix =" loop.

#endif
! End of "NUMORDER == 4".
! ----------------------------------------------------------------------

      END SUBROUTINE RECONST_POS_RHO_PRS_MHD_QUAD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)
      PROGRAM TEST_POS_RHO_PRS_MHD

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER oned_begin, oned_end, oned_off, igeom, ix

      REAL dtmove, pi, rhofloor, prsfloor, smallnum, gamma,
     1     gasconst

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1                     corrupted_zone

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1      xb_use, yb_use, yb_use_m1, zb_use,
     1      flatten, gamma_zonecen, molewt_zonecen

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR, NUM_ZONECEN_BASES) :: w

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND, 5) :: entropy_1d

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NVAR) :: vmin, vmax

! ----------------------------------------------------------------------

      oned_begin = 1
      oned_end = 10
      oned_off = 0
      igeom = 1

      dtmove = 1.0

      pi = 4.0 * ATAN ( 1.0)
      rhofloor = 1.0e-8
      prsfloor = 1.0e-12
      smallnum = 1.0e-14
      gamma = 1.4
      gasconst = 1.0

      gamma_zonecen ( :) = gamma
      molewt_zonecen ( :) = 1.0

      flatten ( :) = 0.0
      w ( :, :, :) = 0.0

! ----------------------------------------------------------------------

      DO ix = oned_begin - oned_off, oned_end + oned_off

      dx_use ( ix) = 1.0
      dy_use ( ix) = 1.0
      dz_use ( ix) = 1.0

      xc_use ( ix) = ( ix - 0.5) * dx_use ( ix)
      yc_use ( ix) = 2.5
      zc_use ( ix) = 3.5

      xb_use  ( ix) = xc_use ( ix) + 0.5 * dx_use ( ix)
      xb_use  ( ix - 1) = xc_use ( ix) - 0.5 * dx_use ( ix)
      yb_use ( ix) = 2.5
      yb_use_m1 ( ix) = 2.0
      zb_use ( ix) = 3.5

! ----------

      w ( ix, 1, 1) = 1.0
      vmin ( ix, 1) = 0.9
      vmax ( ix, 1) = 1.1

      w ( ix, 2, 1) = 1.0
      w ( ix, 3, 1) = -2.0
      w ( ix, 4, 1) = 3.0

      w ( ix, 6, 1) = 2.0
      w ( ix, 7, 1) = 4.0
      w ( ix, 8, 1) = -1.0

      w ( ix, 5, 1) = 2.0 / ( gamma - 1.0)
     1              + 0.5 * ( w ( ix, 2, 1)**2 + w ( ix, 3, 1)**2
     1                      + w ( ix, 4, 1)**2 ) / w ( ix, 1, 1)
     1              + ( w ( ix, 6, 1)**2 + w ( ix, 7, 1)**2
     1                + w ( ix, 8, 1)**2 ) / ( 8.0 * pi)
      vmin ( ix, 5) = 1.5
      vmax ( ix, 5) = 2.5

      END DO

! ------------------------------

      ix = 5
      w ( ix, 1, 2) = 2.0
      w ( ix, 5, 2) = 4.0 / ( gamma - 1.0)

! ----------------------------------------------------------------------

      CALL RECONST_POS_RHO_PRS_MHD_CUBIC

     1           ( oned_begin, oned_end, oned_off,
     1             dtmove, igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gasconst,

     1             xc_use, yc_use, zc_use, dx_use, dy_use, dz_use,
     1             xb_use, yb_use, yb_use_m1, zb_use,

     1             flatten, w, entropy_1d,
     1             gamma_zonecen, molewt_zonecen,
     1             vmin, vmax, corrupted_zone)

! ----------------------------------------------------------------------

      ix = 5
      WRITE ( 6, *)"density slope = ", w ( ix, 1, 2)
      WRITE ( 6, *)"energy density slope = ", w ( ix, 5, 2)


! ----------------------------------------------------------------------

      END PROGRAM TEST_POS_RHO_PRS_MHD
#endif

! ----------------------------------------------------------------------



