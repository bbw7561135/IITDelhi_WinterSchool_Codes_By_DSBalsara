! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! The "TRUE_2D" is only set to 1 if the problem is truly 2d with no
! variation in the third direction. It enables the 2d RS to operate more
! efficiently in the dimensionality where it comes into its own.
#define TRUE_2D 1

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM RIEMANN_MULTID_MHD

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This is a very fast, three dimensional, second order accurate, fully
! functional MHD code that does divergence-free MHD using a
! second order accurate in time ADER timestep. The spatial interpolation
! is either TVD or WENO. This code was put together to go along with
! Balsara's text book.
! Several sample MHD test problems are provided with the code. 

! Code operates in 1, 2 or 3 dimensions. A dimension is switched off
! by putting "ioff? = 0" and "i?1 = 1" in that direction. Thus to switch
! off the z-direction, set "ioffz = 0" and "iz1 = 1". If a direction
! is dynamically active, it is a good idea to set "ioff? = 4" in that
! direction.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone, expand_by_one_zone_corrector,

     1        n_rad_components, n_rad_modes, n_super_timestep_fld,

     1        ixx, iyy, izz, i_comp, i_rad, ifluid, istep, ntstep,
     1        nimgstep, ndumpstep

! Set size of problem domain.

      PARAMETER ( ix1 = 128, iy1 = 128, iz1 = 1,
     1            ioffx = 4, ioffy = 4, ioffz = 0,
     1            n_cc_components = NFIELD + 1 + NFLUID,
     1            n_cc_modes = 5, n_fc_modes = 3,
     1            n_cc_bfld_modes = 7)

      REAL time, dtcur, dtcurp1, timestop, timeimage, dtimage,

     1     xmin, xmax, ymin, ymax, zmin, zmax,

     1     dt_shift,

     1     timing_info1, timing_info2, a,
     1     OMP_GET_WTIME, timing_info_par1, timing_info_par2

      CHARACTER*8 rhofil, vlxfil, vlyfil, vlzfil, prsfil,
     1            bfxfil, bfyfil, bfzfil, dumpfil,
     1            speciefil, temp_speciefil,
     1            radfil, temp_radfil

! ------------------------------

! Big static arrays declared past this point.

! "u" holds zone-centered solution as well as all its space-time moments.
! "u ( ixx, iyy, izz, 6:8, 1)" hold the zone-averaged magnetic fields.
! The order of variables in "u" is "density, x,y,z-momentum densities,
! energy density, x,y,z-magnetic fields, species density".

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

! "entropygr" holds the entropy density and is useful when pressure
! positivity switch is turned on.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

! "bfxbdy, bfybdy, bfzbdy" hold facial magnetic fields and their
! moments in the transverse directions. "bfxbdy" is collocated at the
! center of the x-face and so on...

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

! "bfxccmodes, bfyccmodes, bfzccmodes" hold the divergence-free,
! cell-centered modal representation of the magnetic field.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1      bfxccmodes, bfyccmodes, bfzccmodes

! "flux_x, flux_y, flux_z" hold the fluxes in each of the three directions.
! We also need one "flux_temp" array to parallelize in OpenMP.

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) ::
     1      flux_x, flux_y, flux_z, flux_temp

! "ef?_??_xbdygr, ef?_??_ybdygr, ef?_??_zbdygr" hold the electric field
! components evaluated at the faces.
! "efxgr, efygr, efzgr" hold the edge-centered electric fields.
! "prs_flux_x, prs_flux_y, prs_flux_z" hold the fluxes in each of the
! three directions.

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::

     1      efxgr, efygr, efzgr,

     1      prs_flux_x, prs_flux_y, prs_flux_z,
     1      entropy_flux_x, entropy_flux_y, entropy_flux_z,
     1      entropy_flux_temp

! Zone centered primitive variables and some information for divergence of
! velocity and flattener. These are derived variables.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, msonicgr, divvelgr, flattengr

! Real*4 variable used for imaging variables.

      REAL*4, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: r4gr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  2,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  2) ::
     1         z_indx_limits, bcarr_z

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

      OPEN ( unit = 6, file = "riemann.out")

! It is safest to expand all dynamically active domains by one zone.

      expand_by_one_zone = 1

! ----------------------------------------------------------------------

! Set up some parameters for running the code.
! These parameters control the physics.

      igeom = 1    ! Choose geometry: 1,2,3: Cartesian, Cylindrical, Spherical

      pi = 4.0 * ATAN ( 1.0)

      gamma = 1.4     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 1.4    ! Set up mole wt. of the gas.

      gamma_euler = 1.4     ! Polytropic index for the gas.
      gasconst_euler = 1.0  ! Set up scaled gas constant.
      molewt_euler = 1.4    ! Set up mole wt. of the gas.

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      isospeed = 1.0  ! Set the isothermal sound speed for this simulation.
      barotropic_stiffness = 0.0  ! Set the barotropic_stiffness.

      lightspd = 3e10  ! Speed of light. Scale to suit problem.
      stefan = 5.6697e-5  ! Stefan-Boltzmann constant.

! ------------------------------

! These parameters control the numerics.

! Courant numbers: "cfl_coef_sdsp" modulates the sound(fast msonic) speed;
! "cfl_coef_advct" modulates the advection speed.
! When multi-d Riemann solvers are used, they can be set in the ranges
! [0.3,0.45] in 3d; [0.5,0.65] in 2d and [0.6,0.8] in 1d.
      cfl_coef_sdsp = 0.8
      cfl_coef_advct = 0.8

! If the MC_beta limiter is being used, decide on how much of a compression
! factor one wants to put in. Ranges between 1.0 and 2.0.
! 1.0 makes it MinMod; 2.0 makes it full van Leer MC.
      mc_coef = 1.5

! Choose an entropy fix for the linearized Riemann solver (1, 2 or 3).
      ientropyfix = 1

! If the divergence of the velocity at a zone center exceeds
! this coef times the sound speed then we have to flatten the slopes.
      flatten_coef = 0.4

! This is the factor by which the pressure on one or the other side
! of a zone has to be larger. If this condition is met, we look to
! that side to increase the flattener in the zone being considered.
      prs_flatten_coef = 1.5

      fracdiff = 0.1  !Set fractional change needed for entropy fix.(deprecated)

! Start with a reasonable value for the number of super-timesteps for
! the radiation.
      n_super_timestep_fld = 5

! If the divergence of the velocity at a zone boundary exceeds this coef times
! the sound speed then we have to blend in some fraction of HLL RS.
      divratio_einfeldt_rs = 0.4

! Set a threshold value that triggers the Lapidus viscosity as well as
! its value. Set "lapidus_coef" between [0.0, 0.2]. Default is 0.0.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05

! Using 100% multid flux can sometimes mean we take flux contributions from
! a stencil that does not depend on the zones that contribute to that face.
! For that reason, we blend multid RS with 1d RS associated with that face.
! Blending in more multidimensionality also seems to help out with pressure
! positivity without any artificial fixes.
! This value should range between [0.4, 0.6] for good performance.

! "blend_multid_flux = 1.0 / 6.0" means we use Simpson rule with 1/6 of the
! flux from each edge.
! "blend_multid_flux = 0.0" means we use full flux from suitable 1d RS.
! The electric field is always fully evaluated by the multid. RS.

      blend_multid_flux = 0.1666

! In the fully subsonic case (just) the multi-d HLLC RS can oversteepen
! (almost or fully) stationary contacts. To prevent that we blend in a small
! fraction of the multi-d HLL RS. Ranging between [0.1, 0.2]

! While the flattener also plays a role; "blend_mdhll_to_mdhllc = 1.0"
! converts HLLC into a pure HLL. "blend_mdhll_to_mdhllc = 0.0" does the
! reverse in smooth flow.

      blend_mdhll_to_mdhllc = 0.0

! ------------------------------

! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 50000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc.
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.
      dtimage = 10.0 / 20.0

! Time at which you want the simulation to stop.
      timestop = 10.0

! ------------------------------

! For serial version, set the following to defaults.
      num_level = 0    ! level number is zero.
      num_grid = 0     ! grid number is zero.
      num_pe = 0       ! processor number is zero.

! ------------------------------

! Set default multispecies values, just to make sure they are not NANS.
! For an actual application these will have to be reset to relevant values.

      gamma_specie ( :) = gamma
      molewt_specie ( :) = molewt
      heat_of_formation ( :) = 0.0
      cp_specie ( :) = gasconst * gamma_specie ( :)
     1            / ( ( gamma_specie ( :) - 1.0) * molewt_specie ( :) )
      cv_specie ( :) = cp_specie ( :) / gamma_specie ( :)

! ----------

      gamma_specie_euler ( :) = gamma_euler
      molewt_specie_euler ( :) = molewt_euler
      heat_of_formation_euler ( :) = 0.0
      cp_specie_euler ( :) = gasconst_euler * gamma_specie_euler ( :)
     1 / ( ( gamma_specie_euler ( :) - 1.0) * molewt_specie_euler ( :))
      cv_specie_euler ( :) = cp_specie_euler ( :)
     1                     / gamma_specie_euler ( :)

! ----------------------------------------------------------------------

! Write out some data.

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting problem with the following parameters"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ix1, iy1, iz1 = ", ix1, iy1, iz1
      WRITE ( 6, *)"ioffx, ioffy, ioffz = ", ioffx, ioffy, ioffz
      WRITE ( 6, *)"igeom = ", igeom
      WRITE ( 6, *)" "
      WRITE ( 6, *)"expand_by_one_zone = ", expand_by_one_zone
      WRITE ( 6, *)" "
      WRITE ( 6, *)"pi, rhofloor, prsfloor, smallnum = ",
     1              pi, rhofloor, prsfloor, smallnum
      WRITE ( 6, *)"gamma, gasconst, molewt = ",gamma, gasconst, molewt
      WRITE ( 6, *)" "
      WRITE ( 6, *)"isospeed, barotropic_stiffness = ",
     1              isospeed, barotropic_stiffness
      WRITE ( 6, *)" "
      WRITE ( 6, *)"lightspd, stefan = ", lightspd, stefan
      WRITE ( 6, *)" "
      WRITE ( 6, *)"cfl_coef_sdsp, cfl_coef_advct = ",
     1              cfl_coef_sdsp, cfl_coef_advct
      WRITE ( 6, *)" "
      WRITE ( 6, *)"mc_coef = ", mc_coef
      WRITE ( 6, *)" "
      WRITE ( 6, *)"flatten_coef, divratio_einfeldt_rs = ",
     1              flatten_coef, divratio_einfeldt_rs
      WRITE ( 6, *)" "
      WRITE ( 6, *)"lapidus_threshold, lapidus_coef = ",
     1              lapidus_threshold, lapidus_coef
      WRITE ( 6, *)" "
      WRITE ( 6, *)"blend_multid_flux, blend_mdhll_to_mdhllc = ",
     1              blend_multid_flux, blend_mdhll_to_mdhllc
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ientropyfix, fracdiff = ", ientropyfix, fracdiff
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ntstep, ndumpstep, nimgstep, time, timestop = ",
     1              ntstep, ndumpstep, nimgstep, time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Algorithms chosen for this run : "
      WRITE ( 6, *)" "

#if ( SLOPE_LIMITER == 1)
      WRITE ( 6, *)"TVD limiter being used."
#else
      WRITE ( 6, *)"WENO limiter being used."
#endif
      WRITE ( 6, *)" "

#if ( USE_RIEMANN_SOLVER == 1)
      WRITE ( 6, *)"Linearized Riemann solver being used."
#elif ( USE_RIEMANN_SOLVER == 2)
      WRITE ( 6, *)"HLLD Riemann solver being used."
#else
      WRITE ( 6, *)"HLL Riemann solver being used."
#endif
      WRITE ( 6, *)" "

#if ( DIV_FREE_RECONSTRUCT == 1)
      WRITE ( 6, *)"Face-centered slopes used with div-free",
     1             " reconstruction for the B-field."
#else
      WRITE ( 6, *)"Zone-centered slopes used w/o div-free", 
     1             " reconstruction for the B-field."
#endif
      WRITE ( 6, *)" "

#if ( MULTID_LIMITING == 1)
      WRITE ( 6, *)"Multidimensional TVD limiting will be applied",
     1             " to all variables."
#else
      WRITE ( 6, *)"Multidimensional TVD limiting not being used."
#endif
      WRITE ( 6, *)" "

#if ( PRESSURE_POSITIVITY == 1)
      WRITE ( 6, *)"Positivity of pressure will be enforced as much",
     1             " as possible."
#else
      WRITE ( 6, *)"Positivity of pressure will not be enforced."
#endif
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

! If dimensioning is inconsistent, then write out some info. and stop.

      IF ( ( ix1 .GT. 1) .AND. ( ioffx .LT. 4) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : increase ioffx = ", ioffx
      STOP
      END IF

      IF ( ( ix1 .EQ. 1) .AND. ( ioffx .NE. 0) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : ioffx must be 0 ", ioffx
      STOP
      END IF

      IF ( ( iy1 .GT. 1) .AND. ( ioffy .LT. 4) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : increase ioffy = ", ioffy
      STOP
      END IF

      IF ( ( iy1 .EQ. 1) .AND. ( ioffy .NE. 0) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : ioffy must be 0 ", ioffy
      STOP
      END IF

      IF ( ( iz1 .GT. 1) .AND. ( ioffz .LT. 4) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : increase ioffz = ", ioffz
      STOP
      END IF

      IF ( ( iz1 .EQ. 1) .AND. ( ioffz .NE. 0) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : ioffz must be 0 ", ioffz
      STOP
      END IF

      IF ( NTRACERS > NFLUID) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : NTRACER","S > NFLUI","D",
     1              NTRACERS, NFLUID
      STOP
      END IF

      IF ( ( TRUE_2D == 1) .AND. ( ioffx * ioffy * ioffz > 0) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MULTID_MHD : TRUE _2D == 1 for 3d problem"
      STOP
      END IF

! ----------------------------------------------------------------------

! Set the size of the patch being simulated.

      xmin = -5.0
      xmax = 5.0
      ymin = -5.0
      ymax = 5.0
      zmin = -5.0
      zmax = 5.0

! Because this is a serial or shared memory version, the problem domain 
! is the same as the patch being set up.

      base_grid_xmin = xmin
      base_grid_xmax = xmax
      base_grid_ymin = ymin
      base_grid_ymax = ymax
      base_grid_zmin = zmin
      base_grid_zmax = zmax

! ----------------------------------------------------------------------

! Set file names for the image files. Note that this code just dumps
! out zone-centered variables for imaging. Can't restart from these files.

      rhofil = "rhoa0001"
      vlxfil = "vlxa0001"
      vlyfil = "vlya0001"
      vlzfil = "vlza0001"
      prsfil = "prsa0001"

      bfxfil = "bfxa0001"
      bfyfil = "bfya0001"
      bfzfil = "bfza0001"

      speciefil = "sp010001"

      radfil = "rd010001"

! Set file name for the first dump file.

      dumpfil = "dump0001"

! ----------------------------------------------------------------------

! Use the above patch/domain boundary info. to set up the mesh.
! Uniform mesh is default, or provide your own.
! All the routines of this code are designed to work with any form of
! mesh ratioing.

! ----------

      dx_save ( 1) = ( xmax - xmin) / ix1
      xb_save ( 0) = xmin

      DO ixx = 1, ix1 + ioffx
      dx_save ( ixx) = dx_save ( 1)
      xb_save ( ixx) = xb_save ( ixx - 1) + dx_save ( ixx)
      xc_save ( ixx) = 0.5 * ( xb_save ( ixx) + xb_save ( ixx - 1))
      END DO

      DO ixx = 0, 1 - ioffx, -1
      dx_save ( ixx) = dx_save ( 1)
      xb_save ( ixx - 1) = xb_save ( ixx) - dx_save ( ixx)
      xc_save ( ixx) = 0.5 * ( xb_save ( ixx) + xb_save ( ixx - 1))
      END DO

! ----------

      dy_save ( 1) = ( ymax - ymin) / iy1
      yb_save ( 0) = ymin

      DO iyy = 1, iy1 + ioffy
      dy_save ( iyy) = dy_save ( 1)
      yb_save ( iyy) = yb_save ( iyy - 1) + dy_save ( iyy)
      yc_save ( iyy) = 0.5 * ( yb_save ( iyy) + yb_save ( iyy - 1))
      END DO

      DO iyy = 0, 1 - ioffy, -1
      dy_save ( iyy) = dy_save ( 1)
      yb_save ( iyy - 1) = yb_save ( iyy) - dy_save ( iyy)
      yc_save ( iyy) = 0.5 * ( yb_save ( iyy) + yb_save ( iyy - 1))
      END DO

! ----------

      dz_save ( 1) = ( zmax - zmin) / iz1
      zb_save ( 0) = zmin

      DO izz = 1, iz1 + ioffz
      dz_save ( izz) = dz_save ( 1)
      zb_save ( izz) = zb_save ( izz - 1) + dz_save ( izz)
      zc_save ( izz) = 0.5 * ( zb_save ( izz) + zb_save ( izz - 1))
      END DO

      DO izz = 0, 1 - ioffz, -1
      dz_save ( izz) = dz_save ( 1)
      zb_save ( izz - 1) = zb_save ( izz) - dz_save ( izz)
      zc_save ( izz) = 0.5 * ( zb_save ( izz) + zb_save ( izz - 1))
      END DO

! ----------------------------------------------------------------------

! Initialize the boundary conditions as well as the domain interior here.
! See "PAD_BOUNDARY_X,Y,Z" for meanings of boundary condition arrays.
! Boundary condition number 17 is just reserved for "do nothing", i.e.
! for situations when an external driver sets the boundary conditions.

! ------------------------------

! Bottom and top yz planes.

      DO iyy = 1 - ioffy, iy1 + ioffy

        DO izz = 1 - ioffz, iz1 + ioffz

        x_indx_limits ( 1, iyy, izz) = 1
        x_indx_limits ( 2, iyy, izz) = ix1

        IF ( ABS ( xmin - base_grid_xmin) .LT. smallnum) THEN
          bcarr_x ( 1, iyy, izz) = 4
        ELSE
          bcarr_x ( 1, iyy, izz) = 17
        END IF

        IF ( ABS ( xmax - base_grid_xmax) .LT. smallnum) THEN
          bcarr_x ( 2, iyy, izz) = 4
        ELSE
          bcarr_x ( 2, iyy, izz) = 17
        END IF

        END DO ! End DO izz = loop.

      END DO ! End DO iyy = loop.

! ------------------------------

! Bottom and top xz planes.

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO izz = 1 - ioffz, iz1 + ioffz

        y_indx_limits ( ixx, 1, izz) = 1
        y_indx_limits ( ixx, 2, izz) = iy1

        IF ( ABS ( ymin - base_grid_ymin) .LT. smallnum) THEN
          bcarr_y ( ixx, 1, izz) = 4
        ELSE
          bcarr_y ( ixx, 1, izz) = 17
        END IF

        IF ( ABS ( ymax - base_grid_ymax) .LT. smallnum) THEN
          bcarr_y ( ixx, 2, izz) = 4
        ELSE
          bcarr_y ( ixx, 2, izz) = 17
        END IF

        END DO ! End DO izz = loop.

      END DO ! End DO ixx = loop.

! ------------------------------

! Bottom and top xy planes.

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO iyy = 1 - ioffy, iy1 + ioffy

        z_indx_limits ( ixx, iyy, 1) = 1
        z_indx_limits ( ixx, iyy, 2) = iz1

        IF ( ABS ( zmin - base_grid_zmin) .LT. smallnum) THEN
          bcarr_z ( ixx, iyy, 1) = 4
        ELSE
          bcarr_z ( ixx, iyy, 1) = 17
        END IF

        IF ( ABS ( zmax - base_grid_zmax) .LT. smallnum) THEN
          bcarr_z ( ixx, iyy, 2) = 4
        ELSE
          bcarr_z ( ixx, iyy, 2) = 17
        END IF

        END DO ! End DO iyy = loop.

      END DO ! End DO ixx = loop.

! ----------------------------------------------------------------------

! Initialize the domain interior here. This is usually the only part
! that undergoes significant changes from one problem to the next.
! The set-up subroutine should also initialize "dtcur" and it should
! be copied into "dtcurp1".

! For some problems involving stringent blasts, it may help to start
! with a much smaller CFL number and then to let the timestep float upwards.

#if ( RESTART == 0)

      a = cfl_coef_advct         ! Save advective CFL number.
      cfl_coef_advct = 0.05      ! Start problem with reduced CFL number.

      CALL INIT_ISENTROPIC_VORTEX_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

      cfl_coef_advct = a     ! Copy back the original advective CFL number.

      dtcurp1 = dtcur

      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Initialization complete, dtcur = ", dtcur
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

#endif
! End of "RESTART == 0"

! ----------------------------------------------------------------------

! Read in a prior dump file if one has been written out. Increment dump
! file name.

#if ( RESTART == 1)

! ----------

      OPEN ( unit = 5, file = "riemann.in")

      READ ( 5, *) time, dtcur, i_comp, i_rad

      READ ( 5, *) dumpfil

      READ ( 5, *) rhofil
      READ ( 5, *) vlxfil
      READ ( 5, *) vlyfil
      READ ( 5, *) vlzfil
      READ ( 5, *) prsfil
      READ ( 5, *) bfxfil
      READ ( 5, *) bfyfil
      READ ( 5, *) bfzfil
      READ ( 5, *) speciefil
      READ ( 5, *) radfil

      CLOSE ( 5)

      IF ( i_comp .NE. n_cc_components) THEN
      WRITE ( 6, *)"Stopping because : i_comp .NE. n_cc_components ",
     1              i_comp, n_cc_components
      STOP
      END IF

      IF ( i_rad .NE. n_rad_components) THEN
      WRITE ( 6, *)"Stopping because : i_rad .NE. n_rad_components ",
     1              i_rad, n_rad_components
      STOP
      END IF

      timeimage = time

! ----------

      OPEN ( UNIT = 98, file = dumpfil,
     1       FORM = 'UNFORMATTED', STATUS = 'UNKNOWN', ERR = 10 )

      GO TO 20

10    CONTINUE
      WRITE ( 6, *)"Stopping because nonexistent dumpfil = ", dumpfil
      STOP
20    CONTINUE


      READ ( 98) ( ( ( ( u ( ixx, iyy, izz, i_comp, 1), ixx = 1, ix1),
     1    iyy = 1, iy1), izz = 1, iz1), i_comp = 1, n_cc_components)

      READ ( 98) ( ( ( bfxbdy ( ixx, iyy, izz, 1), ixx = 0, ix1),
     1    iyy = 1, iy1), izz = 1, iz1)

      READ ( 98) ( ( ( bfybdy ( ixx, iyy, izz, 1), ixx = 1, ix1),
     1    iyy = 0, iy1), izz = 1, iz1)

      READ ( 98) ( ( ( bfzbdy ( ixx, iyy, izz, 1), ixx = 1, ix1),
     1    iyy = 1, iy1), izz = 0, iz1)

      CLOSE ( 98)

! ----------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *) dumpfil,
     1             " read in with time, dtcur = ", time, dtcur
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Image file names = ", rhofil," ", vlxfil," ",
     1            vlyfil," ", vlzfil," ", prsfil," ", bfxfil," ",
     1            bfyfil," ", bfzfil," ", speciefil," ", radfil
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

! ----------

      CALL INCNAME_8 ( dumpfil, 1)

#endif
! End of "RESTART == 1"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Main timestep loop begins here.

      DO istep = 1, ntstep

      CALL CPU_TIME ( timing_info1)
!      timing_info_par1 = OMP_GET_WTIME ()

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Write out a dump file if "ndumpstep" steps have been taken
! or stopping time "timestop" is reached or this is the last timestep.
! Increment dump file name. Re-initialize "riemann.in" from which
! we can restart the problem so that it can be restarted from the
! most recent dump file.

      IF ( ( MOD ( istep, ndumpstep) .EQ. 1) .OR.
     1     ( time .GE. timestop) .OR.
     1     ( istep .EQ. ntstep) ) THEN

! ----------

      OPEN ( UNIT = 98, file = dumpfil,
     1       FORM = 'UNFORMATTED', STATUS = 'UNKNOWN' )

      WRITE ( 98) ( ( ( ( u ( ixx, iyy, izz, i_comp, 1), ixx = 1, ix1),
     1    iyy = 1, iy1), izz = 1, iz1), i_comp = 1, n_cc_components)

      WRITE ( 98) ( ( ( bfxbdy ( ixx, iyy, izz, 1), ixx = 0, ix1),
     1    iyy = 1, iy1), izz = 1, iz1)

      WRITE ( 98) ( ( ( bfybdy ( ixx, iyy, izz, 1), ixx = 1, ix1),
     1    iyy = 0, iy1), izz = 1, iz1)

      WRITE ( 98) ( ( ( bfzbdy ( ixx, iyy, izz, 1), ixx = 1, ix1),
     1    iyy = 1, iy1), izz = 0, iz1)

      CLOSE ( 98)

! ----------

      WRITE ( 6, *) dumpfil,
     1      " written out with time, dtcur, n_cc_components = ",
     1      time, dtcur, n_cc_components
      WRITE ( 6, *)" ix1, iy1, iz1 = ", ix1, iy1, iz1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Image file names = ", rhofil," ", vlxfil," ",
     1            vlyfil," ", vlzfil," ", prsfil," ", bfxfil," ",
     1            bfyfil," ", bfzfil," ", speciefil," ", radfil
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

! ----------

      OPEN ( unit = 5, file = "riemann.in")

      WRITE ( 5, *) time, dtcur, n_cc_components, n_rad_components

      WRITE ( 5, *) dumpfil

      WRITE ( 5, *) rhofil
      WRITE ( 5, *) vlxfil
      WRITE ( 5, *) vlyfil
      WRITE ( 5, *) vlzfil
      WRITE ( 5, *) prsfil
      WRITE ( 5, *) bfxfil
      WRITE ( 5, *) bfyfil
      WRITE ( 5, *) bfzfil
      WRITE ( 5, *) speciefil
      WRITE ( 5, *) radfil

      CLOSE ( 5)

! ----------

      WRITE ( 6, *) time, dtcur, n_cc_components, n_rad_components

      WRITE ( 6, *) dumpfil

      WRITE ( 6, *) rhofil
      WRITE ( 6, *) vlxfil
      WRITE ( 6, *) vlyfil
      WRITE ( 6, *) vlzfil
      WRITE ( 6, *) prsfil
      WRITE ( 6, *) bfxfil
      WRITE ( 6, *) bfyfil
      WRITE ( 6, *) bfzfil
      WRITE ( 6, *) speciefil
      WRITE ( 6, *) radfil

      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

! ----------

      CALL INCNAME_8 ( dumpfil, 1)

      END IF ! End "MOD ( istep, ndumpstep) .EQ. 1".

! ----------------------------------------------------------------------

! Fill in the ghost cells in all directions. Do this either by initializing
! with inflowing variables or by extrapolating suitably from the interior
! solution.

      CALL PAD_BOUNDARY_X ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bcarr_x, x_indx_limits,

     1        xc_save, yc_save, zc_save)

      CALL PAD_BOUNDARY_Y ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bcarr_y, y_indx_limits,

     1        xc_save, yc_save, zc_save)

      CALL PAD_BOUNDARY_Z ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bcarr_z, z_indx_limits,

     1        xc_save, yc_save, zc_save)

! ----------------------------------------------------------------------

! Limit the slopes of all the flow variables. The time rate of update
! "u ( ixx, iyy, izz, :, 5)" is also zeroed out here. The flattener
! is also built. The primitive variables are built.

      CALL LIMIT_ALL_SLOPES ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

      CALL ERROR_ISENTROPIC_VORTEX_XY

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

! Call subroutine that enforces positivity of density and pressure. It
! can lower the slopes, so remove it if that is not desired.

#if ( PRESSURE_POSITIVITY == 1)
      CALL POSITIVE_RHO_PRS ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
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
#endif

! ----------------------------------------------------------------------
!      IF ( ( time .GE. timestop) .OR.
!     1     ( MOD ( istep, nimgstep) .EQ. 1) ) THEN

      IF ( ( time .GE. timeimage) .OR.
     1     ( time .GE. timestop) .OR.
     1     ( istep .EQ. ntstep) ) THEN
! ----------------------------------------------------------------------

! Images written out and file names incremented here. This is because the
! previous subroutine has just initialized all the primitive variables.
! For 3d data one may optionally choose to write out select slices of interest.

      timeimage = timeimage + dtimage

      WRITE ( 6, *)"Writing image files = ", rhofil," ", vlxfil," ",
     1            vlyfil," ", vlzfil," ", prsfil," ", bfxfil," ",
     1            bfyfil," ", bfzfil," ", speciefil
      WRITE ( 6, *)" "

! -----

      OPEN ( unit = 99, file = rhofil)
      r4gr ( :, :, :) = rhogr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlxfil)
      r4gr ( :, :, :) = vlxgr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlyfil)
      r4gr ( :, :, :) = vlygr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlzfil)
      r4gr ( :, :, :) = vlzgr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = prsfil)
      r4gr ( :, :, :) = prsgr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = bfxfil)
      r4gr ( :, :, :) = bfxgr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = bfyfil)
      r4gr ( :, :, :) = bfygr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = bfzfil)
      r4gr ( :, :, :) = bfzgr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----
#if ( NFLUID >= 1)

      temp_speciefil = speciefil

      DO ifluid = 1, NFLUID

      OPEN ( unit = 99, file = temp_speciefil)
      r4gr ( :, :, :) = u ( :, :, :, NFIELD + 1 + ifluid, 1)
     1                / rhogr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

      CALL INCNAME_SP ( temp_speciefil, 1)

      END DO

#endif
! ----------

      CALL INCNAME_8 ( rhofil, 1)
      CALL INCNAME_8 ( vlxfil, 1)
      CALL INCNAME_8 ( vlyfil, 1)
      CALL INCNAME_8 ( vlzfil, 1)
      CALL INCNAME_8 ( prsfil, 1)

      CALL INCNAME_8 ( bfxfil, 1)
      CALL INCNAME_8 ( bfyfil, 1)
      CALL INCNAME_8 ( bfzfil, 1)

      CALL INCNAME_8 ( speciefil, 1)

      CALL INCNAME_8 ( radfil, 1)

! ----------------------------------------------------------------------
      END IF ! End of "( time .GT. timestop)..." condition.
! ----------------------------------------------------------------------

! This is the principal stopping point in this code.

      IF ( time .GE. timestop) THEN
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Stopping in RIEMANN_MULTID_MHD: time .GE. timestop"
      WRITE ( 6, *)"time, timestop = ", time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      STOP
      END IF

! ----------------------------------------------------------------------

! Build the first, i.e. predictor, stage.

! ----------------------------------------------------------------------


      CALL MAKE_ADER_UPDT ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! The flux contribution to the predictor step was built in the previous
! subroutine. Here we add in the contribution from the sources.

      CALL ADD_PREDICTOR_SOURCE_TERMS

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Now that the time rate of update "u ( ixx, iyy, izz, :, 5)" is built
! within a zone, we are in a position to build time-centered fluxes,
! i.e. with "dt_shift = 0.5". These can be built at face centers or
! at the edges of the face. The latter give a more compact representation
! of the electric field.
! The fluxes in this corrector stage do need to be zeroed out before
! subroutine entry.

! ------------------------------

      dt_shift = 0.5

      expand_by_one_zone_corrector = 0

      flux_x ( :, :, :, :) = 0.0
      prs_flux_x ( :, :, :) = 0.0
      entropy_flux_x ( :, :, :) = 0.0

      flux_y ( :, :, :, :) = 0.0
      prs_flux_y ( :, :, :) = 0.0
      entropy_flux_y ( :, :, :) = 0.0

      flux_z ( :, :, :, :) = 0.0
      prs_flux_z ( :, :, :) = 0.0
      entropy_flux_z ( :, :, :) = 0.0

! ------------------------------

      flux_temp ( :, :, :, :) = 0.0
      entropy_flux_temp ( :, :, :) = 0.0

      CALL APPLY_MULTID_RS_X_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_y, flux_z, flux_temp, efxgr,
     1        entropy_flux_y, entropy_flux_z, entropy_flux_temp,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      flux_temp ( :, :, :, :) = 0.0
      entropy_flux_temp ( :, :, :) = 0.0

      CALL APPLY_MULTID_RS_Y_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_x, flux_z, flux_temp, efygr,
     1        entropy_flux_x, entropy_flux_z, entropy_flux_temp,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      flux_temp ( :, :, :, :) = 0.0
      entropy_flux_temp ( :, :, :) = 0.0

      CALL APPLY_MULTID_RS_Z_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_x, flux_y, flux_temp, efzgr,
     1        entropy_flux_x, entropy_flux_z, entropy_flux_temp,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

      CALL MAKE_FLUX_X ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur, dt_shift,

     1        u, entropygr, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_x, prs_flux_x, entropy_flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      CALL MAKE_FLUX_Y ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur, dt_shift,

     1        u, entropygr, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_y, prs_flux_y, entropy_flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      CALL MAKE_FLUX_Z ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur, dt_shift,

     1        u, entropygr, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_z, prs_flux_z, entropy_flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

      CALL ADD_LAPIDUS_X ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      CALL ADD_LAPIDUS_Y ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      CALL ADD_LAPIDUS_Z ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Update the flow variables in this corrector step and get new "dtcurp1".

      CALL UPDATE_VARS_CORRECTOR_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone_corrector,
     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy, flux_x, flux_y, flux_z,
     1        prs_flux_x, prs_flux_y, prs_flux_z,
     1       entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z,

     1        efxgr, efygr, efzgr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Make the next timestep that will be minimized in the ensuing subroutine.
! This is also where the primitive variables are refreshed and pressure
! positivity, if it is called for, is reasserted.

      dtcurp1 = 1.25 * dtcur

      CALL EVALUATE_TIMESTEP

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        dtcur, dtcurp1,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------

! Update "time" and write out a small amount of timestep info.
! Gather timing information.
! If the "time" is about to exceed "timestop" in the next time step,
! curtail the next timestep so as to stop exactly at "timestop".

      CALL CPU_TIME ( timing_info2)
!      timing_info_par2 = OMP_GET_WTIME ()

      time = time + dtcur

! If you plan to restart the problem then comment out this set of lines.
! They are only useful for accuracy analysis and nowhere else.
!      IF ( time + dtcurp1 .GT. timestop) THEN
!      dtcurp1 = dtcurp1 - ( time + dtcurp1 - timestop) + smallnum
!      dtcurp1 = AMAX1 ( dtcurp1, smallnum / 100.0)
!      timeimage = time
!      END IF

      WRITE ( 6, *)"Timestep istep = ", istep," time = ", time,
     1             " dtcur = ", dtcur," dtcurp1 = ", dtcurp1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed this step : ", timing_info2 - timing_info1
!      WRITE ( 6, *)"Timed thread 0: ",timing_info_par2-timing_info_par1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

      CALL FLUSH ( 6)

      dtcur = dtcurp1

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO ! End of DO istep = loop.

! Main timestep loop ends here.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      STOP

      END PROGRAM RIEMANN_MULTID_MHD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PAD_BOUNDARY_X ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bcarr_x, x_indx_limits,

     1        xc_save, yc_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! "bcarr_x" controls the boundary conditions along the x-direction:
! 1 for inflow; 2 for outflow/continuitive; 3 for reflective;
! 4 for periodic in serial setting.

! "x_indx_limits" gives the starting and ending indices in the x-direction
! at which the boundary conditions are to be applied. This makes it
! easy to set up interior reflecting wall type boundary conditions.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe

      REAL time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

      REAL xc_save ( 1 - ioffx : ix1 + ioffx)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, ilhs, irhs, oned_begin, oned_end,
     1        nspecies_min, nspecies_max

      REAL rho_tmp, momx_tmp, momy_tmp, momz_tmp, eng_tmp,
     1     bfx_tmp, bfy_tmp, bfz_tmp, species_tmp ( 0: NFLUID)

! ----------------------------------------------------------------------

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  & 
!$OMP&    eng_tmp, bfx_tmp, bfy_tmp, bfz_tmp, species_tmp)              &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    x_indx_limits, bcarr_x, u, bfxbdy, bfybdy, bfzbdy)            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

      DO iyy = 1 - ioffy, iy1 + ioffy

      oned_begin = x_indx_limits ( 1, iyy, izz)
      oned_end = x_indx_limits ( 2, iyy, izz)

! ----------------------------------------------------------------------

! Inflow at lower x-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 1.0
      momx_tmp = 1.0
      momy_tmp = 1.0
      momz_tmp = 1.0
      eng_tmp = 1.0
      bfx_tmp = 1.0
      bfy_tmp = 1.0
      bfz_tmp = 1.0
      species_tmp ( :) = 0.0

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 1) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx

      u ( ilhs, iyy, izz, 1, 1) = rho_tmp
      u ( ilhs, iyy, izz, 2, 1) = momx_tmp
      u ( ilhs, iyy, izz, 3, 1) = momy_tmp
      u ( ilhs, iyy, izz, 4, 1) = momz_tmp
      u ( ilhs, iyy, izz, 5, 1) = eng_tmp
      u ( ilhs, iyy, izz, 6, 1) = bfx_tmp
      u ( ilhs, iyy, izz, 7, 1) = bfy_tmp
      u ( ilhs, iyy, izz, 8, 1) = bfz_tmp
#if ( NFLUID >= 1)
      u ( ilhs, iyy, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID)
#endif

      bfxbdy ( ilhs - 1, iyy, izz, 1) = bfx_tmp

      bfybdy ( ilhs, iyy, izz, 1) = bfy_tmp
      bfybdy ( ilhs, iyy - 1, izz, 1) = bfy_tmp

      bfzbdy ( ilhs, iyy, izz, 1) = bfz_tmp
      bfzbdy ( ilhs, iyy, izz - 1, 1) = bfz_tmp

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at lower x-boundary.

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 2) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx
      irhs = oned_begin

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

      bfxbdy ( ilhs - 1, iyy, izz, 1) =
     1                          bfxbdy ( irhs - 1, iyy, izz, 1)

      bfybdy ( ilhs, iyy, izz, 1) = bfybdy ( irhs, iyy, izz, 1)
      bfybdy ( ilhs, iyy - 1, izz, 1) =
     1                          bfybdy ( irhs, iyy - 1, izz, 1)

      bfzbdy ( ilhs, iyy, izz, 1) = bfzbdy ( irhs, iyy, izz, 1)
      bfzbdy ( ilhs, iyy, izz - 1, 1) =
     1                          bfzbdy ( irhs, iyy, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Reflective at lower x-boundary.
! Flip longitudinal velocity and transverse magnetic fields.

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 3) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx
      irhs = oned_begin + ixx - 1

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)
      u ( ilhs, iyy, izz, 2, 1) = - u ( ilhs, iyy, izz, 2, 1)

      bfxbdy ( ilhs - 1, iyy, izz, 1) = bfxbdy ( irhs, iyy, izz, 1)

      bfybdy ( ilhs, iyy, izz, 1) = - bfybdy ( irhs, iyy, izz, 1)
      bfybdy ( ilhs, iyy - 1, izz, 1) =
     1                          - bfybdy ( irhs, iyy - 1, izz, 1)

      bfzbdy ( ilhs, iyy, izz, 1) = - bfzbdy ( irhs, iyy, izz, 1)
      bfzbdy ( ilhs, iyy, izz - 1, 1) =
     1                          - bfzbdy ( irhs, iyy, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Periodic at lower x-boundary.

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 4) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx
      irhs = oned_end + 1 - ixx

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

      bfxbdy ( ilhs - 1, iyy, izz, 1) =
     1                          bfxbdy ( irhs - 1, iyy, izz, 1)

      bfybdy ( ilhs, iyy, izz, 1) = bfybdy ( irhs, iyy, izz, 1)
      bfybdy ( ilhs, iyy - 1, izz, 1) =
     1                          bfybdy ( irhs, iyy - 1, izz, 1)

      bfzbdy ( ilhs, iyy, izz, 1) = bfzbdy ( irhs, iyy, izz, 1)
      bfzbdy ( ilhs, iyy, izz - 1, 1) =
     1                          bfzbdy ( irhs, iyy, izz - 1, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, bfx_tmp, bfy_tmp, bfz_tmp, species_tmp)              &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    x_indx_limits, bcarr_x, u, bfxbdy, bfybdy, bfzbdy)            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

      DO iyy = 1 - ioffy, iy1 + ioffy

      oned_begin = x_indx_limits ( 1, iyy, izz)
      oned_end = x_indx_limits ( 2, iyy, izz)

! ----------------------------------------------------------------------

! Inflow at upper x-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 1.0
      momx_tmp = 1.0
      momy_tmp = 1.0
      momz_tmp = 1.0
      eng_tmp = 1.0
      bfx_tmp = 1.0
      bfy_tmp = 1.0
      bfz_tmp = 1.0
      species_tmp ( :) = 0.0

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 1) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx

      u ( ilhs, iyy, izz, 1, 1) = rho_tmp
      u ( ilhs, iyy, izz, 2, 1) = momx_tmp
      u ( ilhs, iyy, izz, 3, 1) = momy_tmp
      u ( ilhs, iyy, izz, 4, 1) = momz_tmp
      u ( ilhs, iyy, izz, 5, 1) = eng_tmp
      u ( ilhs, iyy, izz, 6, 1) = bfx_tmp
      u ( ilhs, iyy, izz, 7, 1) = bfy_tmp
      u ( ilhs, iyy, izz, 8, 1) = bfz_tmp
#if ( NFLUID >= 1)
      u ( ilhs, iyy, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID)
#endif

      bfxbdy ( ilhs, iyy, izz, 1) = bfx_tmp

      bfybdy ( ilhs, iyy, izz, 1) = bfy_tmp
      bfybdy ( ilhs, iyy - 1, izz, 1) = bfy_tmp

      bfzbdy ( ilhs, iyy, izz, 1) = bfz_tmp
      bfzbdy ( ilhs, iyy, izz - 1, 1) = bfz_tmp

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at upper x-boundary.

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 2) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx
      irhs = oned_end

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

      bfxbdy ( ilhs, iyy, izz, 1) = bfxbdy ( irhs, iyy, izz, 1)

      bfybdy ( ilhs, iyy, izz, 1) = bfybdy ( irhs, iyy, izz, 1)
      bfybdy ( ilhs, iyy - 1, izz, 1) =
     1                          bfybdy ( irhs, iyy - 1, izz, 1)

      bfzbdy ( ilhs, iyy, izz, 1) = bfzbdy ( irhs, iyy, izz, 1)
      bfzbdy ( ilhs, iyy, izz - 1, 1) =
     1                          bfzbdy ( irhs, iyy, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Reflective at upper x-boundary.
! Flip longitudinal velocity and transverse magnetic fields.

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 3) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx
      irhs = oned_end - ixx + 1

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)
      u ( ilhs, iyy, izz, 2, 1) = - u ( ilhs, iyy, izz, 2, 1)

      bfxbdy ( ilhs, iyy, izz, 1) = bfxbdy ( irhs - 1, iyy, izz, 1)

      bfybdy ( ilhs, iyy, izz, 1) = - bfybdy ( irhs, iyy, izz, 1)
      bfybdy ( ilhs, iyy - 1, izz, 1) =
     1                        - bfybdy ( irhs, iyy - 1, izz, 1)

      bfzbdy ( ilhs, iyy, izz, 1) = - bfzbdy ( irhs, iyy, izz, 1)
      bfzbdy ( ilhs, iyy, izz - 1, 1) =
     1                        - bfzbdy ( irhs, iyy, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Periodic at upper x-boundary.

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 4) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx
      irhs = oned_begin + ixx - 1

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

      bfxbdy ( ilhs, iyy, izz, 1) = bfxbdy ( irhs, iyy, izz, 1)

      bfybdy ( ilhs, iyy, izz, 1) = bfybdy ( irhs, iyy, izz, 1)
      bfybdy ( ilhs, iyy - 1, izz, 1) =
     1                          bfybdy ( irhs, iyy - 1, izz, 1)

      bfzbdy ( ilhs, iyy, izz, 1) = bfzbdy ( irhs, iyy, izz, 1)
      bfzbdy ( ilhs, iyy, izz - 1, 1) =
     1                          bfzbdy ( irhs, iyy, izz - 1, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE PAD_BOUNDARY_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PAD_BOUNDARY_Y ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bcarr_y, y_indx_limits,

     1        xc_save, yc_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! "bcarr_y" controls the boundary conditions along the y-direction:
! 1 for inflow; 2 for outflow/continuitive; 3 for reflective;
! 4 for periodic in serial setting.

! "y_indx_limits" gives the starting and ending indices in the y-direction
! at which the boundary conditions are to be applied. This makes it
! easy to set up interior reflecting wall type boundary conditions.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe

      REAL time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  2,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

      REAL xc_save ( 1 - ioffx : ix1 + ioffx)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, ilhs, irhs, oned_begin, oned_end,
     1        nspecies_min, nspecies_max

      REAL rho_tmp, momx_tmp, momy_tmp, momz_tmp, eng_tmp,
     1     bfx_tmp, bfy_tmp, bfz_tmp, species_tmp ( 0: NFLUID)

! ----------------------------------------------------------------------

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, bfx_tmp, bfy_tmp, bfz_tmp, species_tmp)              &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    y_indx_limits, bcarr_y, u, bfxbdy, bfybdy, bfzbdy)            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

      DO ixx = 1 - ioffx, ix1 + ioffx

      oned_begin = y_indx_limits ( ixx, 1, izz)
      oned_end = y_indx_limits ( ixx, 2, izz)

! ----------------------------------------------------------------------

! Inflow at lower y-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 1.0
      momx_tmp = 1.0
      momy_tmp = 1.0
      momz_tmp = 1.0
      eng_tmp = 1.0
      bfx_tmp = 1.0
      bfy_tmp = 1.0
      bfz_tmp = 1.0
      species_tmp ( :) = 0.0

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 1) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy

      u ( ixx, ilhs, izz, 1, 1) = rho_tmp
      u ( ixx, ilhs, izz, 2, 1) = momx_tmp
      u ( ixx, ilhs, izz, 3, 1) = momy_tmp
      u ( ixx, ilhs, izz, 4, 1) = momz_tmp
      u ( ixx, ilhs, izz, 5, 1) = eng_tmp
      u ( ixx, ilhs, izz, 6, 1) = bfx_tmp
      u ( ixx, ilhs, izz, 7, 1) = bfy_tmp
      u ( ixx, ilhs, izz, 8, 1) = bfz_tmp
#if ( NFLUID >= 1)
      u ( ixx, ilhs, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID)
#endif

      bfxbdy ( ixx, ilhs, izz, 1) = bfx_tmp
      bfxbdy ( ixx - 1, ilhs, izz, 1) = bfx_tmp

      bfybdy ( ixx, ilhs - 1, izz, 1) = bfy_tmp

      bfzbdy ( ixx, ilhs, izz, 1) = bfz_tmp
      bfzbdy ( ixx, ilhs, izz - 1, 1) = bfz_tmp

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at lower y-boundary.

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 2) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy
      irhs = oned_begin

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

      bfxbdy ( ixx, ilhs, izz, 1) = bfxbdy ( ixx, irhs, izz, 1)
      bfxbdy ( ixx - 1, ilhs, izz, 1) =
     1                          bfxbdy ( ixx - 1, irhs, izz, 1)

      bfybdy ( ixx, ilhs - 1, izz, 1) =
     1                          bfybdy ( ixx, irhs - 1, izz, 1)

      bfzbdy ( ixx, ilhs, izz, 1) = bfzbdy ( ixx, irhs, izz, 1)
      bfzbdy ( ixx, ilhs, izz - 1, 1) =
     1                          bfzbdy ( ixx, irhs, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Reflective at lower y-boundary.
! Flip longitudinal velocity and transverse magnetic fields.

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 3) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy
      irhs = oned_begin + iyy - 1

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)
      u ( ixx, ilhs, izz, 3, 1) = - u ( ixx, ilhs, izz, 3, 1)

      bfxbdy ( ixx, ilhs, izz, 1) = - bfxbdy ( ixx, irhs, izz, 1)
      bfxbdy ( ixx - 1, ilhs, izz, 1) =
     1                        - bfxbdy ( ixx - 1, irhs, izz, 1)

      bfybdy ( ixx, ilhs - 1, izz, 1) = bfybdy ( ixx, irhs, izz, 1)

      bfzbdy ( ixx, ilhs, izz, 1) = - bfzbdy ( ixx, irhs, izz, 1)
      bfzbdy ( ixx, ilhs, izz - 1, 1) =
     1                        - bfzbdy ( ixx, irhs, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Periodic at lower y-boundary.

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 4) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy
      irhs = oned_end + 1 - iyy

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

      bfxbdy ( ixx, ilhs, izz, 1) = bfxbdy ( ixx, irhs, izz, 1)
      bfxbdy ( ixx - 1, ilhs, izz, 1) =
     1                          bfxbdy ( ixx - 1, irhs, izz, 1)

      bfybdy ( ixx, ilhs - 1, izz, 1) =
     1                          bfybdy ( ixx, irhs - 1, izz, 1)

      bfzbdy ( ixx, ilhs, izz, 1) = bfzbdy ( ixx, irhs, izz, 1)
      bfzbdy ( ixx, ilhs, izz - 1, 1) =
     1                          bfzbdy ( ixx, irhs, izz - 1, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, bfx_tmp, bfy_tmp, bfz_tmp, species_tmp)              &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    y_indx_limits, bcarr_y, u, bfxbdy, bfybdy, bfzbdy)            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

      DO ixx = 1 - ioffx, ix1 + ioffx

      oned_begin = y_indx_limits ( ixx, 1, izz)
      oned_end = y_indx_limits ( ixx, 2, izz)

! ----------------------------------------------------------------------

! Inflow at upper y-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 1.0
      momx_tmp = 1.0
      momy_tmp = 1.0
      momz_tmp = 1.0
      eng_tmp = 1.0
      bfx_tmp = 1.0
      bfy_tmp = 1.0
      bfz_tmp = 1.0
      species_tmp ( :) = 0.0

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 1) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy

      u ( ixx, ilhs, izz, 1, 1) = rho_tmp
      u ( ixx, ilhs, izz, 2, 1) = momx_tmp
      u ( ixx, ilhs, izz, 3, 1) = momy_tmp
      u ( ixx, ilhs, izz, 4, 1) = momz_tmp
      u ( ixx, ilhs, izz, 5, 1) = eng_tmp
      u ( ixx, ilhs, izz, 6, 1) = bfx_tmp
      u ( ixx, ilhs, izz, 7, 1) = bfy_tmp
      u ( ixx, ilhs, izz, 8, 1) = bfz_tmp
#if ( NFLUID >= 1)
      u ( ixx, ilhs, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID)
#endif

      bfxbdy ( ixx, ilhs, izz, 1) = bfx_tmp
      bfxbdy ( ixx - 1, ilhs, izz, 1) = bfx_tmp

      bfybdy ( ixx, ilhs, izz, 1) = bfy_tmp

      bfzbdy ( ixx, ilhs, izz, 1) = bfz_tmp
      bfzbdy ( ixx, ilhs, izz - 1, 1) = bfz_tmp

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at upper y-boundary.

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 2) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy
      irhs = oned_end

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

      bfxbdy ( ixx, ilhs, izz, 1) = bfxbdy ( ixx, irhs, izz, 1)
      bfxbdy ( ixx - 1, ilhs, izz, 1) =
     1                          bfxbdy ( ixx - 1, irhs, izz, 1)

      bfybdy ( ixx, ilhs, izz, 1) = bfybdy ( ixx, irhs, izz, 1)

      bfzbdy ( ixx, ilhs, izz, 1) = bfzbdy ( ixx, irhs, izz, 1)
      bfzbdy ( ixx, ilhs, izz - 1, 1) =
     1                          bfzbdy ( ixx, irhs, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Reflective at upper y-boundary.
! Flip longitudinal velocity and transverse magnetic fields.

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 3) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy
      irhs = oned_end - iyy + 1

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)
      u ( ixx, ilhs, izz, 3, 1) = - u ( ixx, ilhs, izz, 3, 1)

      bfxbdy ( ixx, ilhs, izz, 1) = - bfxbdy ( ixx, irhs, izz, 1)
      bfxbdy ( ixx - 1, ilhs, izz, 1) =
     1                        - bfxbdy ( ixx - 1, irhs, izz, 1)

      bfybdy ( ixx, ilhs, izz, 1) = bfybdy ( ixx, irhs - 1, izz, 1)

      bfzbdy ( ixx, ilhs, izz, 1) = - bfzbdy ( ixx, irhs, izz, 1)
      bfzbdy ( ixx, ilhs, izz - 1, 1) =
     1                        - bfzbdy ( ixx, irhs, izz - 1, 1)

      END DO

      END IF

! ------------------------------

! Periodic at upper y-boundary.

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 4) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy
      irhs = oned_begin + iyy - 1

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

      bfxbdy ( ixx, ilhs, izz, 1) = bfxbdy ( ixx, irhs, izz, 1)
      bfxbdy ( ixx - 1, ilhs, izz, 1) =
     1                          bfxbdy ( ixx - 1, irhs, izz, 1)

      bfybdy ( ixx, ilhs, izz, 1) = bfybdy ( ixx, irhs, izz, 1)

      bfzbdy ( ixx, ilhs, izz, 1) = bfzbdy ( ixx, irhs, izz, 1)
      bfzbdy ( ixx, ilhs, izz - 1, 1) =
     1                          bfzbdy ( ixx, irhs, izz - 1, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE PAD_BOUNDARY_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PAD_BOUNDARY_Z ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        time,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bcarr_z, z_indx_limits,

     1        xc_save, yc_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! "bcarr_z" controls the boundary conditions along the z-direction:
! 1 for inflow; 2 for outflow/continuitive; 3 for reflective;
! 4 for periodic in serial setting.

! "z_indx_limits" gives the starting and ending indices in the z-direction
! at which the boundary conditions are to be applied. This makes it
! easy to set up interior reflecting wall type boundary conditions.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe

      REAL time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  2) ::
     1         z_indx_limits, bcarr_z

      REAL xc_save ( 1 - ioffx : ix1 + ioffx)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, ilhs, irhs, oned_begin, oned_end,
     1        nspecies_min, nspecies_max

      REAL rho_tmp, momx_tmp, momy_tmp, momz_tmp, eng_tmp,
     1     bfx_tmp, bfy_tmp, bfz_tmp, species_tmp ( 0: NFLUID)

! ----------------------------------------------------------------------

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, bfx_tmp, bfy_tmp, bfz_tmp, species_tmp)              &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    z_indx_limits, bcarr_z, u, bfxbdy, bfybdy, bfzbdy)            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 1 - ioffy, iy1 + ioffy

      DO ixx = 1 - ioffx, ix1 + ioffx

      oned_begin = z_indx_limits ( ixx, iyy, 1)
      oned_end = z_indx_limits ( ixx, iyy, 2)

! ----------------------------------------------------------------------

! Inflow at lower z-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 1.0
      momx_tmp = 1.0
      momy_tmp = 1.0
      momz_tmp = 1.0
      eng_tmp = 1.0
      bfx_tmp = 1.0
      bfy_tmp = 1.0
      bfz_tmp = 1.0
      species_tmp ( :) = 0.0

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 1) THEN

      DO izz = 1, ioffz

      ilhs = oned_begin - izz

      u ( ixx, iyy, ilhs, 1, 1) = rho_tmp
      u ( ixx, iyy, ilhs, 2, 1) = momx_tmp
      u ( ixx, iyy, ilhs, 3, 1) = momy_tmp
      u ( ixx, iyy, ilhs, 4, 1) = momz_tmp
      u ( ixx, iyy, ilhs, 5, 1) = eng_tmp
      u ( ixx, iyy, ilhs, 6, 1) = bfx_tmp
      u ( ixx, iyy, ilhs, 7, 1) = bfy_tmp
      u ( ixx, iyy, ilhs, 8, 1) = bfz_tmp
#if ( NFLUID >= 1)
      u ( ixx, iyy, ilhs, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID)
#endif

      bfxbdy ( ixx, iyy, ilhs, 1) = bfx_tmp
      bfxbdy ( ixx - 1, iyy, ilhs, 1) = bfx_tmp

      bfybdy ( ixx, iyy, ilhs, 1) = bfy_tmp
      bfybdy ( ixx, iyy - 1, ilhs, 1) = bfy_tmp

      bfzbdy ( ixx, iyy, ilhs - 1, 1) = bfz_tmp

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at lower z-boundary.

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 2) THEN

      DO izz = 1, ioffz

      ilhs = oned_begin - izz
      irhs = oned_begin

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

      bfxbdy ( ixx, iyy, ilhs, 1) = bfxbdy ( ixx, iyy, irhs, 1)
      bfxbdy ( ixx - 1, iyy, ilhs, 1) =
     1                          bfxbdy ( ixx - 1, iyy, irhs, 1)

      bfybdy ( ixx, iyy, ilhs, 1) = bfybdy ( ixx, iyy, irhs, 1)
      bfybdy ( ixx, iyy - 1, ilhs, 1) =
     1                          bfybdy ( ixx, iyy - 1, irhs, 1)

      bfzbdy ( ixx, iyy, ilhs - 1, 1) =
     1                          bfzbdy ( ixx, iyy, irhs - 1, 1)

      END DO

      END IF

! ------------------------------

! Reflective at lower z-boundary.
! Flip longitudinal velocity and transverse magnetic fields.

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 3) THEN

      DO izz = 1, ioffz 

      ilhs = oned_begin - izz
      irhs = oned_begin + izz - 1

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)
      u ( ixx, iyy, ilhs, 4, 1) = - u ( ixx, iyy, ilhs, 4, 1)

      bfxbdy ( ixx, iyy, ilhs, 1) = - bfxbdy ( ixx, iyy, irhs, 1)
      bfxbdy ( ixx - 1, iyy, ilhs, 1) =
     1                        - bfxbdy ( ixx - 1, iyy, irhs, 1)

      bfybdy ( ixx, iyy, ilhs, 1) = - bfybdy ( ixx, iyy, irhs, 1)
      bfybdy ( ixx, iyy - 1, ilhs, 1) =
     1                        - bfybdy ( ixx, iyy - 1, irhs, 1)

      bfzbdy ( ixx, iyy, ilhs - 1, 1) = bfzbdy ( ixx, iyy, irhs, 1)

      END DO

      END IF

! ------------------------------

! Periodic at lower z-boundary.

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 4) THEN

      DO izz = 1, ioffz

      ilhs = oned_begin - izz
      irhs = oned_end + 1 - izz

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

      bfxbdy ( ixx, iyy, ilhs, 1) = bfxbdy ( ixx, iyy, irhs, 1)
      bfxbdy ( ixx - 1, iyy, ilhs, 1) =
     1                          bfxbdy ( ixx - 1, iyy, irhs, 1)

      bfybdy ( ixx, iyy, ilhs, 1) = bfybdy ( ixx, iyy, irhs, 1)
      bfybdy ( ixx, iyy - 1, ilhs, 1) =
     1                          bfybdy ( ixx, iyy - 1, irhs, 1)

      bfzbdy ( ixx, iyy, ilhs - 1, 1) =
     1                          bfzbdy ( ixx, iyy, irhs - 1, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, bfx_tmp, bfy_tmp, bfz_tmp, species_tmp)              &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    z_indx_limits, bcarr_z, u, bfxbdy, bfybdy, bfzbdy)            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 1 - ioffy, iy1 + ioffy

      DO ixx = 1 - ioffx, ix1 + ioffx

      oned_begin = z_indx_limits ( ixx, iyy, 1)
      oned_end = z_indx_limits ( ixx, iyy, 2)

! ----------------------------------------------------------------------

! Inflow at upper z-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 1.0
      momx_tmp = 1.0
      momy_tmp = 1.0
      momz_tmp = 1.0
      eng_tmp = 1.0
      bfx_tmp = 1.0
      bfy_tmp = 1.0
      bfz_tmp = 1.0
      species_tmp ( :) = 0.0

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 1) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz

      u ( ixx, iyy, ilhs, 1, 1) = rho_tmp
      u ( ixx, iyy, ilhs, 2, 1) = momx_tmp
      u ( ixx, iyy, ilhs, 3, 1) = momy_tmp
      u ( ixx, iyy, ilhs, 4, 1) = momz_tmp
      u ( ixx, iyy, ilhs, 5, 1) = eng_tmp
      u ( ixx, iyy, ilhs, 6, 1) = bfx_tmp
      u ( ixx, iyy, ilhs, 7, 1) = bfy_tmp
      u ( ixx, iyy, ilhs, 8, 1) = bfz_tmp
#if ( NFLUID >= 1)
      u ( ixx, iyy, ilhs, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID)
#endif

      bfxbdy ( ixx, iyy, ilhs, 1) = bfx_tmp
      bfxbdy ( ixx - 1, iyy, ilhs, 1) = bfx_tmp

      bfybdy ( ixx, iyy, ilhs, 1) = bfy_tmp
      bfybdy ( ixx, iyy - 1, ilhs, 1) = bfy_tmp

      bfzbdy ( ixx, iyy, ilhs, 1) = bfz_tmp

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at upper z-boundary.

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 2) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz
      irhs = oned_end

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

      bfxbdy ( ixx, iyy, ilhs, 1) = bfxbdy ( ixx, iyy, irhs, 1)
      bfxbdy ( ixx - 1, iyy, ilhs, 1) =
     1                          bfxbdy ( ixx - 1, iyy, irhs, 1)

      bfybdy ( ixx, iyy, ilhs, 1) = bfybdy ( ixx, iyy, irhs, 1)
      bfybdy ( ixx, iyy - 1, ilhs, 1) =
     1                          bfybdy ( ixx, iyy - 1, irhs, 1)

      bfzbdy ( ixx, iyy, ilhs, 1) = bfzbdy ( ixx, iyy, irhs, 1)

      END DO

      END IF

! ------------------------------

! Reflective at upper z-boundary.
! Flip longitudinal velocity and transverse magnetic fields.

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 3) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz
      irhs = oned_end - izz + 1

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)
      u ( ixx, iyy, ilhs, 4, 1) = - u ( ixx, iyy, ilhs, 4, 1)

      bfxbdy ( ixx, iyy, ilhs, 1) = - bfxbdy ( ixx, iyy, irhs, 1)
      bfxbdy ( ixx - 1, iyy, ilhs, 1) =
     1                        - bfxbdy ( ixx - 1, iyy, irhs, 1)

      bfybdy ( ixx, iyy, ilhs, 1) = - bfybdy ( ixx, iyy, irhs, 1)
      bfybdy ( ixx, iyy - 1, ilhs, 1) =
     1                        - bfybdy ( ixx, iyy - 1, irhs, 1)

      bfzbdy ( ixx, iyy, ilhs, 1) = bfzbdy ( ixx, iyy, irhs - 1, 1)

      END DO

      END IF

! ------------------------------

! Periodic at upper z-boundary.

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 4) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz
      irhs = oned_begin + izz - 1

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

      bfxbdy ( ixx, iyy, ilhs, 1) = bfxbdy ( ixx, iyy, irhs, 1)
      bfxbdy ( ixx - 1, iyy, ilhs, 1) =
     1                          bfxbdy ( ixx - 1, iyy, irhs, 1)

      bfybdy ( ixx, iyy, ilhs, 1) = bfybdy ( ixx, iyy, irhs, 1)
      bfybdy ( ixx, iyy - 1, ilhs, 1) =
     1                          bfybdy ( ixx, iyy - 1, irhs, 1)

      bfzbdy ( ixx, iyy, ilhs, 1) = bfzbdy ( ixx, iyy, irhs, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE PAD_BOUNDARY_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LIMIT_ALL_SLOPES ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Obtains primitive variables from conserved variables.

! Limits the conserved variables "u" in all three directions. Also limits
! the facial magnetic fields "bfxbdy, bfybdy, bfzbdy" in the transverse
! directions. Makes their undivided differences (slopes).
! The conserved variables in "u" are "density, x,y,z-momentum densities,
! energy density, x,y,z-magnetic fields, species density".

! Characteristic limiting can also be done in addition to limiting on
! the conserved variables. It does increase the cost but provides more
! stability. Set LIMIT_TWICE to 1 to invoke additional limiting.

! If pressure positivity switch is on, it also builds up an "entropygr"
! variable which contains the entropy density. The slopes of this variable
! are set conservatively and it is evolved in time with the other variables.
! The caveat is that it is evolved using physical pressures unless the
! pressures go negative. If that happens, the pressures are supplemented
! by the entropy variable.

! Multidimensional limiting can also be switched on. This is useful in
! keeping the flow variation inside a zone within the bounds set by the
! flow variables in the adjacent zones (viewed in a multidimensional sense).
! Amongst other things, it ensures that the density remains positive all
! through a zone and not just at the face centers. Multidimensional limiting
! helps with stability at the cost of degrading the accuracy, so use it sparingly.

! Makes a cell-centered divergence-free reconstruction of the magnetic field
! and stores it in "bfxccmodes, bfyccmodes, bfzccmodes". This is described in
! detail in Balsara (2001) JCP, (2004) ApJSupp, (2009) JCP. The zone-centered
! reconstruction of face-centered magnetic fields yields the following modes:

! bfx (x,y,z) = bfx_1 + bfx_2 x + bfx_3 y + bfx_4 z
!          + bfx_5 (x^2-1/12) + bfx_6 x y + bfx_7 x z

! bfy (x,y,z) = bfy_1 + bfy_2 x + bfy_3 y + bfy_4 z
!          + bfy_5 (y^2-1/12) + bfy_6 x y + bfy_7 y z

! bfz (x,y,z) = bfz_1 + bfz_2 x + bfz_3 y + bfz_4 z
!          + bfz_5 (z^2-1/12) + bfz_6 y z + bfz_7 x z

! SLOPE_LIMITER if set to 1, gives TVD limiting.
! Otherwise, gives the linear portion of r=3 centered WENO limiting.

! "flattengr" applies a flattener in regions of strong shock.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

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

      INTEGER ixx, iyy, izz, i_comp, ixx0, iyy0, izz0,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: workgr

      REAL ob_8pi, psi_limiter, central_zone_wt, ob_6, r13_12,
     1     MinMod, MC_beta, ta, tb,
     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax, div_min, div_max,
     1     flat_min, flat_max

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, flatten_1d,

     1     a, b, c, d, e, f, g, h,
     1     var_zone_min, var_zone_max,
     1     var_neighbor_min, var_neighbor_max,

     1     u_ip2, u_ip1, u_i, u_im1, u_im2,
     1     R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3,

     1     dx, dy, dz, ob_dx, ob_dy, ob_dz,

     1     bfx_pl_0, bfx_pl_y, bfx_pl_z, bfx_mi_0, bfx_mi_y, bfx_mi_z,
     1     bfy_pl_0, bfy_pl_x, bfy_pl_z, bfy_mi_0, bfy_mi_x, bfy_mi_z,
     1     bfz_pl_0, bfz_pl_x, bfz_pl_y, bfz_mi_0, bfz_mi_x, bfz_mi_y

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID, 2) :: u_1d

! MinMod limiter
!      MinMod ( ta, tb) = AMIN1 ( ABS ( ta), ABS ( tb))
!     1                    * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( ta, tb, mc_coef) = AMIN1 ( 0.5 * ( ABS ( ta + tb)),
     1                     mc_coef * ABS ( ta), mc_coef * ABS ( tb) )
     1                   * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! ----------------------------------------------------------------------

! Variables that go through the EOS subroutine boundary.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      rho, eng, prs, bprs, temperature, entropy,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)

      ob_6 = 1.0 / 6.0

      r13_12 = 13.0 / 12.0

      central_zone_wt = 100.0

! ----------

      IF ( n_cc_components .NE. NFIELD + 1 + NFLUID) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_SLOPES : n_cc_components = ",
     1              n_cc_components
      STOP
      END IF

      IF ( n_cc_modes .NE. 5) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_SLOPES : n_cc_modes = ",
     1              n_cc_modes
      STOP
      END IF

      IF ( n_fc_modes .NE. 3) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_SLOPES : n_fc_modes = ",
     1              n_fc_modes
      STOP
      END IF

      IF ( n_cc_bfld_modes .NE. 7) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_SLOPES : n_cc_bfld_modes = ",
     1              n_cc_bfld_modes
      STOP
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES : NFIE","LD"
      WRITE ( 6, *) NFIELD, ".NE. 7"
      STOP
      END IF

! ----------------------------------------------------------------------

! Make zone centered primitive variables. The 6, 7, 8 components of the
! conserved variables contain the zone-averaged magnetic fields.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, rho, eng, bprs, specie, prs,         &
!$OMP&   temperature, entropy, gamma_eos, gamma_soundspeed, mean_mol_wt)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ob_8pi, gasconst, gamma, smallnum, prsfloor, rhofloor, molewt,&
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, entropygr,                            &
!$OMP&    temperaturegr, gamma_eos_gr, gamma_soundspeed_gr,             &
!$OMP&    mean_mol_wt_gr, msonicgr)                                     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          u ( ixx, iyy, izz, :, 5) = 0.0 ! Zero out ADER's time rate of update.

          u ( ixx, iyy, izz, 1, 1) = AMAX1 ( rhofloor,
     1                               u ( ixx, iyy, izz, 1, 1))

          rhogr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 1, 1)

          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)

          vlxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 2, 1) * a ( ixx)
          vlygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 3, 1) * a ( ixx)
          vlzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 4, 1) * a ( ixx)

          bfxgr ( ixx, iyy, izz) =0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                                  + bfxbdy ( ixxm1, iyy, izz, 1))
          bfygr ( ixx, iyy, izz) =0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                                  + bfybdy ( ixx, iyym1, izz, 1))
          bfzgr ( ixx, iyy, izz) =0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                                  + bfzbdy ( ixx, iyy, izzm1, 1))

          u ( ixx, iyy, izz, 6, 1) = bfxgr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 7, 1) = bfygr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 8, 1) = bfzgr ( ixx, iyy, izz)

          b ( ixx) = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          bprs ( ixx) = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          rho ( ixx) = rhogr ( ixx, iyy, izz)

          eng ( ixx) = ( u ( ixx, iyy, izz, 5, 1) - b ( ixx)
     1                                            - bprs ( ixx) )

#if ( NFLUID >= 1)
          specie ( 1:NFLUID, ixx) =
     1         u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) * a ( ixx)
#endif

          END DO ! End of DO ixx = loop.

! ----------

! Choose between a simple or more complicated EOS. Species fraction 
! can be strictly restricted to lie between 0 and 1 if it is copied back.


#if ( NFLUID == 0)

          DO ixx = 1 - ioffx, ix1 + ioffx

          gamma_eos ( ixx) = gamma
          gamma_soundspeed ( ixx) = gamma
          mean_mol_wt ( ixx) = molewt

          prs ( ixx) = AMAX1 ( eng ( ixx) * ( gamma - 1.0), prsfloor)
          temperature ( ixx) = molewt * prs ( ixx)
     1                       / ( gasconst * rho ( ixx))

          END DO ! End of DO ixx = loop.

#else

      CALL RHO_ENG_IN_PRS_T_GAMA_OUT
     1        ( 1, ix1, ioffx,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rho, eng, specie,

     1          prs, temperature, entropy,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed)
          
#endif

! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)

          prsgr ( ixx, iyy, izz) = prs ( ixx)
          temperaturegr ( ixx, iyy, izz) = temperature ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
! This is the only place in the code where the "entropygr" variable is
! rebuilt for the present timestep.
          entropygr ( ixx, iyy, izz, 1) = prs ( ixx)
     1                         / rho ( ixx)**( gamma_eos ( ixx) - 1.0)
          entropygr ( ixx, iyy, izz, 2: 5) = 0.0
#endif

          gamma_eos_gr ( ixx, iyy, izz) = gamma_eos ( ixx)
          gamma_soundspeed_gr ( ixx, iyy, izz) = gamma_soundspeed ( ixx)
          mean_mol_wt_gr ( ixx, iyy, izz) = mean_mol_wt ( ixx)

#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) =
     1    specie ( 1:NFLUID, ixx) * rhogr ( ixx, iyy, izz)
#endif

! Choose either the full magnetosonic speed or sound speed for the flattener.
          msonicgr ( ixx, iyy, izz) =
!     1    ( 2.0 * bprs ( ixx) + gamma_soundspeed ( ixx)
!     1                      * prs ( ixx) ) * a ( ixx)     ! full msonic speed
     1     gamma_soundspeed ( ixx) * prs ( ixx) * a ( ixx) ! sound speed
          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! Make the undivided divergence of the velocity, "divvelgr". Use it along
! with the minimal value of "msonicgr" in a halo of zones around the
! current one to make the flattener, "flattengr".

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, dx, dy, dz, ob_dx, ob_dy, ob_dz)     &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, flatten_coef, xb_save, dx_save, xc_save,               &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    vlxgr, vlygr, vlzgr, divvelgr, msonicgr, flattengr, workgr)   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

! ----------

          IF ( igeom .EQ. 1) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 2) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 3) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1              * dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

! ----------

          divvelgr ( ixx, iyy, izz) =
     1    ( vlxgr ( ixxp1, iyy, izz) - vlxgr ( ixxm1, iyy, izz) )
     1  * ob_dx ( ixx)
     1  + ( vlygr ( ixx, iyyp1, izz) - vlygr ( ixx, iyym1, izz) )
     1  * ob_dy ( ixx)
     1  + ( vlzgr ( ixx, iyy, izzp1) - vlzgr ( ixx, iyy, izzm1) )
     1  * ob_dz ( ixx)

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)
     1  * AMIN1 ( dx ( ixx), dy ( ixx), dz ( ixx) )

          workgr ( ixx, iyy, izz) = AMIN1 (

     1  msonicgr ( ixxm1, iyy, izz), msonicgr ( ixx, iyy, izz),
     1  msonicgr ( ixxp1, iyy, izz),
     1  msonicgr ( ixxm1, iyyp1, izz), msonicgr ( ixx, iyyp1, izz),
     1  msonicgr ( ixxp1, iyyp1, izz),
     1  msonicgr ( ixxm1, iyym1, izz), msonicgr ( ixx, iyym1, izz),
     1  msonicgr ( ixxp1, iyym1, izz),

     1  msonicgr ( ixxm1, iyy, izzp1), msonicgr ( ixx, iyy, izzp1),
     1  msonicgr ( ixxp1, iyy, izzp1),
     1  msonicgr ( ixxm1, iyyp1, izzp1), msonicgr ( ixx, iyyp1, izzp1),
     1  msonicgr ( ixxp1, iyyp1, izzp1),
     1  msonicgr ( ixxm1, iyym1, izzp1), msonicgr ( ixx, iyym1, izzp1),
     1  msonicgr ( ixxp1, iyym1, izzp1),

     1  msonicgr ( ixxm1, iyy, izzm1), msonicgr ( ixx, iyy, izzm1),
     1  msonicgr ( ixxp1, iyy, izzm1),
     1  msonicgr ( ixxm1, iyyp1, izzm1), msonicgr ( ixx, iyyp1, izzm1),
     1  msonicgr ( ixxp1, iyyp1, izzm1),
     1  msonicgr ( ixxm1, iyym1, izzm1), msonicgr ( ixx, iyym1, izzm1),
     1  msonicgr ( ixxp1, iyym1, izzm1)

     1                                    )

          a ( ixx) = ABS ( divvelgr ( ixx, iyy, izz))
     1      / ( flatten_coef * workgr ( ixx, iyy, izz))

          a ( ixx) = a ( ixx) - 1.0

          IF ( a ( ixx) .LE. 0.0) THEN
            flattengr ( ixx, iyy, izz) = 0.0
          ELSE
            flattengr ( ixx, iyy, izz) = AMIN1 ( 1.0, a ( ixx))
          END IF

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Upwind the flattener in the x-direction.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ixxp1, ixxm1, flatten_1d)    &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    prs_flatten_coef, prsgr, flattengr)                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! -----

          DO ixx = 1 - ioffx, ix1 + ioffx

          flatten_1d ( ixx) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO ixx = loop.

! -----

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          IF ( prsgr ( ixxm1, iyy, izz) .GT.
     1         prs_flatten_coef * prsgr ( ixx, iyy, izz) ) THEN

          flattengr ( ixx, iyy, izz) =
     1                    AMAX1 ( flattengr ( ixx, iyy, izz),
     1                            flatten_1d ( ixxm1) )

          END IF

          IF ( prsgr ( ixxp1, iyy, izz) .GT.
     1         prs_flatten_coef * prsgr ( ixx, iyy, izz) ) THEN

          flattengr ( ixx, iyy, izz) =
     1                    AMAX1 ( flattengr ( ixx, iyy, izz),
     1                            flatten_1d ( ixxp1) )

          END IF

          END DO ! End of DO ixx = loop.

! -----

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Upwind the flattener in the y-direction.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, iyyp1, iyym1, ixx, flatten_1d)    &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    prs_flatten_coef, prsgr, flattengr)                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO ixx = 1 - ioffx, ix1 + ioffx

! -----

          DO iyy = 1 - ioffy, iy1 + ioffy

          flatten_1d ( iyy) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO iyy = loop.

! -----

          DO iyy = 1 - ioffy, iy1 + ioffy
          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          IF ( prsgr ( ixx, iyym1, izz) .GT.
     1         prs_flatten_coef * prsgr ( ixx, iyy, izz) ) THEN

          flattengr ( ixx, iyy, izz) =
     1                    AMAX1 ( flattengr ( ixx, iyy, izz),
     1                            flatten_1d ( iyym1) )

          END IF

          IF ( prsgr ( ixx, iyyp1, izz) .GT.
     1         prs_flatten_coef * prsgr ( ixx, iyy, izz) ) THEN

          flattengr ( ixx, iyy, izz) =
     1                    AMAX1 ( flattengr ( ixx, iyy, izz),
     1                            flatten_1d ( iyyp1) )

          END IF

          END DO ! End of DO iyy = loop.

! -----

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Upwind the flattener in the z-direction.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, ixx, flatten_1d)    &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    prs_flatten_coef, prsgr, flattengr)                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO iyy = 1 - ioffy, iy1 + ioffy

        DO ixx = 1 - ioffx, ix1 + ioffx

! -----

          DO izz = 1 - ioffz, iz1 + ioffz

          flatten_1d ( izz) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO izz = loop.

! -----

          DO izz = 1 - ioffz, iz1 + ioffz
          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

          IF ( prsgr ( ixx, iyy, izzm1) .GT.
     1         prs_flatten_coef * prsgr ( ixx, iyy, izz) ) THEN

          flattengr ( ixx, iyy, izz) =
     1                    AMAX1 ( flattengr ( ixx, iyy, izz),
     1                            flatten_1d ( izzm1) )

          END IF

          IF ( prsgr ( ixx, iyy, izzp1) .GT.
     1         prs_flatten_coef * prsgr ( ixx, iyy, izz) ) THEN

          flattengr ( ixx, iyy, izz) =
     1                    AMAX1 ( flattengr ( ixx, iyy, izz),
     1                            flatten_1d ( izzp1) )

          END IF

          END DO ! End of DO izz = loop.

! -----

        END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ------------------------------

! If so desired, smear out the effects of the flattener by one zone in
! all directions.

      workgr ( :, :, :) = flattengr ( :, :, :)

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1)                                            &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    workgr, flattengr)                                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

! ----------

          flattengr ( ixx, iyy, izz) = AMAX1 (

     1  workgr ( ixxm1, iyy, izz), workgr ( ixx, iyy, izz),
     1  workgr ( ixxp1, iyy, izz),
     1  workgr ( ixxm1, iyyp1, izz), workgr ( ixx, iyyp1, izz),
     1  workgr ( ixxp1, iyyp1, izz),
     1  workgr ( ixxm1, iyym1, izz), workgr ( ixx, iyym1, izz),
     1  workgr ( ixxp1, iyym1, izz),

     1  workgr ( ixxm1, iyy, izzp1), workgr ( ixx, iyy, izzp1),
     1  workgr ( ixxp1, iyy, izzp1),
     1  workgr ( ixxm1, iyyp1, izzp1), workgr ( ixx, iyyp1, izzp1),
     1  workgr ( ixxp1, iyyp1, izzp1),
     1  workgr ( ixxm1, iyym1, izzp1), workgr ( ixx, iyym1, izzp1),
     1  workgr ( ixxp1, iyym1, izzp1),

     1  workgr ( ixxm1, iyy, izzm1), workgr ( ixx, iyy, izzm1),
     1  workgr ( ixxp1, iyy, izzm1),
     1  workgr ( ixxm1, iyyp1, izzm1), workgr ( ixx, iyyp1, izzm1),
     1  workgr ( ixxp1, iyyp1, izzm1),
     1  workgr ( ixxm1, iyym1, izzm1), workgr ( ixx, iyym1, izzm1),
     1  workgr ( ixxp1, iyym1, izzm1)

     1                                    )

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! TVD/r=3 WENO limiters at second order start here.
! ----------------------------------------------------------------------

! Now make the slopes (with flattening) for the zone-centered
! conserved variables. If facial slopes are needed for the magnetic field
! components because they can be done later. If zone-centered slopes
! are needed for the magnetic field, then do them here.

! ------------------------------
#if ( SLOPE_LIMITER == 1)
! ------------------------------
! MC_beta Limiter.

      DO i_comp = 1, n_cc_components

      IF ( ( i_comp .LE. 5) .OR. ( i_comp .GE. 9) .OR.
     1     ( DIV_FREE_RECONSTRUCT .NE. 1) ) THEN

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c)                                   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    i_comp, mc_coef, u, flattengr)                                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! --------------------
! Loop split just to help with vectorization.

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

          a ( ixx) = u ( ixxp1, iyy, izz, i_comp, 1)
     1             - u ( ixx, iyy, izz, i_comp, 1)

          b ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
     1             - u ( ixxm1, iyy, izz, i_comp, 1)

          END DO

! --------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

! ----------

          c ( ixx) = 1.0 - flattengr ( ixx, iyy, izz)

! ----------

          u ( ixx, iyy, izz, i_comp, 2) =
     1    MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = u ( ixx, iyyp1, izz, i_comp, 1)
     1             - u ( ixx, iyy, izz, i_comp, 1)

          b ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
     1             - u ( ixx, iyym1, izz, i_comp, 1)

          u ( ixx, iyy, izz, i_comp, 3) =
     1    MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = u ( ixx, iyy, izzp1, i_comp, 1)
     1             - u ( ixx, iyy, izz, i_comp, 1)

          b ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
     1             - u ( ixx, iyy, izzm1, i_comp, 1)

          u ( ixx, iyy, izz, i_comp, 4) =
     1    MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

! --------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! Don't limit magnetic field components.

      END DO ! End of DO i_comp = loop.

! ------------------------------
#else
! ------------------------------
! WENO Limiter.

      DO i_comp = 1, n_cc_components

      IF ( ( i_comp .LE. 5) .OR. ( i_comp .GE. 9) .OR.
     1     ( DIV_FREE_RECONSTRUCT .NE. 1) ) THEN

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                     &
!$OMP&    R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3)               &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    i_comp, central_zone_wt, smallnum, r13_12, u, flattengr)      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! --------------------
! Loop split just to help with vectorization.

          DO ixx = 1 - ioffx, ix1 + ioffx

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

          u_ip2 ( ixx) = u ( ixxp2, iyy, izz, i_comp, 1)
          u_ip1 ( ixx) = u ( ixxp1, iyy, izz, i_comp, 1)
          u_i ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
          u_im1 ( ixx) = u ( ixxm1, iyy, izz, i_comp, 1)
          u_im2 ( ixx) = u ( ixxm2, iyy, izz, i_comp, 1)

          END DO

! --------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 1 - ioffz)

! ----------

          c = 1.0 - flattengr ( ixx, iyy, izz)

! ----------

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1        + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1        + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1        + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1        - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1        + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          u ( ixx, iyy, izz, i_comp, 2) =
     1             ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1             + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          u_ip2 ( ixx) = u ( ixx, iyyp2, izz, i_comp, 1)
          u_ip1 ( ixx) = u ( ixx, iyyp1, izz, i_comp, 1)
          u_i ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
          u_im1 ( ixx) = u ( ixx, iyym1, izz, i_comp, 1)
          u_im2 ( ixx) = u ( ixx, iyym2, izz, i_comp, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                  + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx)
     1        + u_im2 ( ixx) + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                   + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1        + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                  - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1        + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          u ( ixx, iyy, izz, i_comp, 3) =
     1        ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1        + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          u_ip2 ( ixx) = u ( ixx, iyy, izzp2, i_comp, 1)
          u_ip1 ( ixx) = u ( ixx, iyy, izzp1, i_comp, 1)
          u_i ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
          u_im1 ( ixx) = u ( ixx, iyy, izzm1, i_comp, 1)
          u_im2 ( ixx) = u ( ixx, iyy, izzm2, i_comp, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx)
     1                  + 0.5 * u_im2 ( ixx) + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx)
     1        + u_im2 ( ixx) + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                   + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1               + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                  - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          u ( ixx, iyy, izz, i_comp, 4) =
     1         ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1         + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------

      END IF ! Don't limit magnetic field components.

      END DO ! End of DO i_comp = loop.

! ------------------------------
#endif
! ----------------------------------------------------------------------
! Do this part only if facial slopes are desired for the
! divergence-free reconstruction. The alternative is to use zone-centered
! transverse slopes even for the magnetic fields.
#if ( DIV_FREE_RECONSTRUCT == 1)
! ----------------------------------------------------------------------

! Now make the slopes (with flattening) for the x-magnetic field
! in the y and z-directions and save them in that order.

! ------------------------------
#if ( SLOPE_LIMITER == 1)
! ------------------------------
! MC_beta Limiter.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c)                                   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    mc_coef, bfxbdy, flattengr)                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 0 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

! ----------

          c ( ixx) = 1.0
     1      - 0.5 * ( flattengr ( MAX0 ( ixx, 1 - ioffx), iyy, izz)
     1              + flattengr ( ixxp1, iyy, izz))

! ----------

          a ( ixx) = bfxbdy ( ixx, iyyp1, izz, 1)
     1             - bfxbdy ( ixx, iyy, izz, 1)

          b ( ixx) = bfxbdy ( ixx, iyy, izz, 1)
     1             - bfxbdy ( ixx, iyym1, izz, 1)

          bfxbdy ( ixx, iyy, izz, 2) =
     1           MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = bfxbdy ( ixx, iyy, izzp1, 1)
     1             - bfxbdy ( ixx, iyy, izz, 1)

          b ( ixx) = bfxbdy ( ixx, iyy, izz, 1)
     1             - bfxbdy ( ixx, iyy, izzm1, 1)

          bfxbdy ( ixx, iyy, izz, 3) =
     1           MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#else
! ------------------------------
! WENO Limiter.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                     &
!$OMP&    R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3)               &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    central_zone_wt, smallnum, r13_12, bfxbdy, flattengr)         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 0 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 0 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 1 - ioffz)

! ----------

          c ( ixx) = 1.0
     1      - 0.5 * ( flattengr ( MAX0 ( ixx, 1 - ioffx), iyy, izz)
     1              + flattengr ( ixxp1, iyy, izz))

! ----------

          u_ip2 ( ixx) = bfxbdy ( ixx, iyyp2, izz, 1)
          u_ip1 ( ixx) = bfxbdy ( ixx, iyyp1, izz, 1)
          u_i ( ixx) = bfxbdy ( ixx, iyy, izz, 1)
          u_im1 ( ixx) = bfxbdy ( ixx, iyym1, izz, 1)
          u_im2 ( ixx) = bfxbdy ( ixx, iyym2, izz, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                          + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                 + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          bfxbdy ( ixx, iyy, izz, 2) = 
     1           ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1           + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          u_ip2 ( ixx) = bfxbdy ( ixx, iyy, izzp2, 1)
          u_ip1 ( ixx) = bfxbdy ( ixx, iyy, izzp1, 1)
          u_i ( ixx) = bfxbdy ( ixx, iyy, izz, 1)
          u_im1 ( ixx) = bfxbdy ( ixx, iyy, izzm1, 1)
          u_im2 ( ixx) = bfxbdy ( ixx, iyy, izzm2, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                          + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                 + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          bfxbdy ( ixx, iyy, izz, 3) =
     1           ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx) 
     1           + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#endif
! ------------------------------

! ----------------------------------------------------------------------

! Now make the slopes (with flattening) for the y-magnetic field
! in the x and z-directions and save them in that order.

! ------------------------------
#if ( SLOPE_LIMITER == 1)
! ------------------------------
! MC_beta Limiter.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c)                                   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    mc_coef, bfybdy, flattengr)                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 0 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 0 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

! ----------

          c ( ixx) = 1.0
     1      - 0.5 * ( flattengr ( ixx, MAX0 ( iyy, 1 - ioffy), izz)
     1              + flattengr ( ixx, iyyp1, izz))

! ----------

          a ( ixx) = bfybdy ( ixxp1, iyy, izz, 1)
     1             - bfybdy ( ixx, iyy, izz, 1)

          b ( ixx) = bfybdy ( ixx, iyy, izz, 1)
     1             - bfybdy ( ixxm1, iyy, izz, 1)

          bfybdy ( ixx, iyy, izz, 2) =
     1           MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = bfybdy ( ixx, iyy, izzp1, 1)
     1             - bfybdy ( ixx, iyy, izz, 1)

          b ( ixx) = bfybdy ( ixx, iyy, izz, 1)
     1             - bfybdy ( ixx, iyy, izzm1, 1)

          bfybdy ( ixx, iyy, izz, 3) =
     1           MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#else
! ------------------------------
! WENO Limiter.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                     &
!$OMP&    R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3)               &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    central_zone_wt, smallnum, r13_12, bfybdy, flattengr)         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 0 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 0 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 0 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 1 - ioffz)

! ----------

          c ( ixx) = 1.0
     1      - 0.5 * ( flattengr ( ixx, MAX0 ( iyy, 1 - ioffy), izz)
     1              + flattengr ( ixx, iyyp1, izz))

! ----------

          u_ip2 ( ixx) = bfybdy ( ixxp2, iyy, izz, 1)
          u_ip1 ( ixx) = bfybdy ( ixxp1, iyy, izz, 1)
          u_i ( ixx) = bfybdy ( ixx, iyy, izz, 1)
          u_im1 ( ixx) = bfybdy ( ixxm1, iyy, izz, 1)
          u_im2 ( ixx) = bfybdy ( ixxm2, iyy, izz, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                          + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                 + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          bfybdy ( ixx, iyy, izz, 2) =
     1           ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1           + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          u_ip2 ( ixx) = bfybdy ( ixx, iyy, izzp2, 1)
          u_ip1 ( ixx) = bfybdy ( ixx, iyy, izzp1, 1)
          u_i ( ixx) = bfybdy ( ixx, iyy, izz, 1)
          u_im1 ( ixx) = bfybdy ( ixx, iyy, izzm1, 1)
          u_im2 ( ixx) = bfybdy ( ixx, iyy, izzm2, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                          + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                 + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          bfybdy ( ixx, iyy, izz, 3) =
     1           ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1           + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#endif
! ------------------------------

! ----------------------------------------------------------------------

! Now make the slopes (with flattening) for the z-magnetic field
! in the x and y-directions and save them in that order.

! ------------------------------
#if ( SLOPE_LIMITER == 1)
! ------------------------------
! MC_beta Limiter.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c)                                   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    mc_coef, bfzbdy, flattengr)                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 0 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 0 - ioffz)

! ----------

          c ( ixx) = 1.0
     1      - 0.5 * ( flattengr ( ixx, iyy, MAX0 ( izz, 1 - ioffz))
     1              + flattengr ( ixx, iyy, izzp1))

! ----------

          a ( ixx) = bfzbdy ( ixxp1, iyy, izz, 1)
     1             - bfzbdy ( ixx, iyy, izz, 1)

          b ( ixx) = bfzbdy ( ixx, iyy, izz, 1)
     1             - bfzbdy ( ixxm1, iyy, izz, 1)

          bfzbdy ( ixx, iyy, izz, 2) =
     1           MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = bfzbdy ( ixx, iyyp1, izz, 1)
     1             - bfzbdy ( ixx, iyy, izz, 1)

          b ( ixx) = bfzbdy ( ixx, iyy, izz, 1)
     1             - bfzbdy ( ixx, iyym1, izz, 1)

          bfzbdy ( ixx, iyy, izz, 3) =
     1           MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#else
! ------------------------------
! WENO Limiter.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                     &
!$OMP&    R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3)               &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    central_zone_wt, smallnum, r13_12, bfzbdy, flattengr)         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 0 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 0 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 0 - ioffz)

! ----------

          c ( ixx) = 1.0
     1      - 0.5 * ( flattengr ( ixx, iyy, MAX0 ( izz, 1 - ioffz))
     1              + flattengr ( ixx, iyy, izzp1))

! ----------

          u_ip2 ( ixx) = bfzbdy ( ixxp2, iyy, izz, 1)
          u_ip1 ( ixx) = bfzbdy ( ixxp1, iyy, izz, 1)
          u_i ( ixx) = bfzbdy ( ixx, iyy, izz, 1)
          u_im1 ( ixx) = bfzbdy ( ixxm1, iyy, izz, 1)
          u_im2 ( ixx) = bfzbdy ( ixxm2, iyy, izz, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                          + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                 + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          bfzbdy ( ixx, iyy, izz, 2) =
     1           ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1           + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          u_ip2 ( ixx) = bfzbdy ( ixx, iyyp2, izz, 1)
          u_ip1 ( ixx) = bfzbdy ( ixx, iyyp1, izz, 1)
          u_i ( ixx) = bfzbdy ( ixx, iyy, izz, 1)
          u_im1 ( ixx) = bfzbdy ( ixx, iyym1, izz, 1)
          u_im2 ( ixx) = bfzbdy ( ixx, iyym2, izz, 1)

          R1_x ( ixx) = - 2.0 * u_im1 ( ixx) + 0.5 * u_im2 ( ixx)
     1                + 1.5 * u_i ( ixx)
          Is1 ( ixx) = r13_12 * ( - 2.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                          + u_i ( ixx))**2
     1        + 0.25 * ( - 4.0 * u_im1 ( ixx) + u_im2 ( ixx)
     1                 + 3.0 * u_i ( ixx))**2

          R2_x ( ixx) = 0.5 * ( u_ip1 ( ixx) - u_im1 ( ixx) )
          Is2 ( ixx) = r13_12 * ( u_im1 ( ixx) - 2.0 * u_i ( ixx)
     1                          + u_ip1 ( ixx) )**2
     1        + 0.25 * ( u_im1 ( ixx) - u_ip1 ( ixx) )**2

          R3_x ( ixx) = - 1.5 * u_i ( ixx) + 2.0 * u_ip1 ( ixx)
     1                - 0.5 * u_ip2 ( ixx)
          Is3 ( ixx) = r13_12 * ( u_i ( ixx) - 2.0 * u_ip1 ( ixx)
     1                          + u_ip2 ( ixx))**2
     1        + 0.25 * ( 3.0 * u_i ( ixx) - 4.0 * u_ip1 ( ixx)
     1                 + u_ip2 ( ixx))**2

          Wt1 ( ixx) = 1.0 / ( Is1 ( ixx) + smallnum)**4
          Wt2 ( ixx) = central_zone_wt / ( Is2 ( ixx) + smallnum)**4
          Wt3 ( ixx) = 1.0 / ( Is3 ( ixx) + smallnum)**4

          a ( ixx) = 1.0 / ( Wt1 ( ixx) + Wt2 ( ixx) + Wt3 ( ixx))
          Wt1 ( ixx) = a ( ixx) * Wt1 ( ixx)
          Wt2 ( ixx) = a ( ixx) * Wt2 ( ixx)
          Wt3 ( ixx) = a ( ixx) * Wt3 ( ixx)

          bfzbdy ( ixx, iyy, izz, 3) =
     1           ( Wt1 ( ixx) * R1_x ( ixx) + Wt2 ( ixx) * R2_x ( ixx)
     1           + Wt3 ( ixx) * R3_x ( ixx) ) * c ( ixx)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#endif
! ------------------------------

! ----------------------------------------------------------------------
#endif
! End of "DIV_FREE_RECONSTRUCT == 1".
! ----------------------------------------------------------------------
! TVD/r=3 WENO limiters at second order end here.
! ----------------------------------------------------------------------
! The multidimensional limiter, if invoked, restores positivity in a
! multidimensional sense.
! ----------------------------------------------------------------------
! The better limiters can, on rare occasions, produce -ve densities at the
! edges/vertices. This can become a problem if limiters are invoked at those
! edges/vertices. 
#if ( MULTID_LIMITING == 1)

! ------------------------------

      DO i_comp = 1, n_cc_components

      IF ( ( i_comp .LE. 5) .OR. ( i_comp .GE. 9) .OR.
     1     ( DIV_FREE_RECONSTRUCT .NE. 1) ) THEN

      psi_limiter = 0.5    ! For tough problems, use this value for all vars.

!      IF ( ( i_comp .EQ. 1) .OR. ( i_comp .EQ. 6) .OR. ( i_comp .EQ. 7)
!     1     .OR. ( i_comp .EQ. 8) ) psi_limiter = 1.0

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d, e, f, g, h,                    &
!$OMP&    var_zone_min, var_zone_max, var_neighbor_min,var_neighbor_max)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    i_comp, smallnum, psi_limiter, u)                             &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! --------------------
! Loop split just to help with vectorization.

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)

! ----------

          var_neighbor_min ( ixx) =
     1                  AMIN1 ( u ( ixx, iyy, izz, i_comp, 1),
     1                          u ( ixxp1, iyy, izz, i_comp, 1),
     1                          u ( ixxm1, iyy, izz, i_comp, 1),
     1                          u ( ixx, iyyp1, izz, i_comp, 1),
     1                          u ( ixxp1, iyyp1, izz, i_comp, 1),
     1                          u ( ixxm1, iyyp1, izz, i_comp, 1),
     1                          u ( ixx, iyym1, izz, i_comp, 1),
     1                          u ( ixxp1, iyym1, izz, i_comp, 1),
     1                          u ( ixxm1, iyym1, izz, i_comp, 1),

     1                          u ( ixx, iyy, izzp1, i_comp, 1),
     1                          u ( ixxp1, iyy, izzp1, i_comp, 1),
     1                          u ( ixxm1, iyy, izzp1, i_comp, 1),
     1                          u ( ixx, iyyp1, izzp1, i_comp, 1),
     1                          u ( ixxp1, iyyp1, izzp1, i_comp, 1),
     1                          u ( ixxm1, iyyp1, izzp1, i_comp, 1),
     1                          u ( ixx, iyym1, izzp1, i_comp, 1),
     1                          u ( ixxp1, iyym1, izzp1, i_comp, 1),
     1                          u ( ixxm1, iyym1, izzp1, i_comp, 1),

     1                          u ( ixx, iyy, izzm1, i_comp, 1),
     1                          u ( ixxp1, iyy, izzm1, i_comp, 1),
     1                          u ( ixxm1, iyy, izzm1, i_comp, 1),
     1                          u ( ixx, iyyp1, izzm1, i_comp, 1),
     1                          u ( ixxp1, iyyp1, izzm1, i_comp, 1),
     1                          u ( ixxm1, iyyp1, izzm1, i_comp, 1),
     1                          u ( ixx, iyym1, izzm1, i_comp, 1),
     1                          u ( ixxp1, iyym1, izzm1, i_comp, 1),
     1                          u ( ixxm1, iyym1, izzm1, i_comp, 1) )


          var_neighbor_max ( ixx) =
     1                  AMAX1 ( u ( ixx, iyy, izz, i_comp, 1),
     1                          u ( ixxp1, iyy, izz, i_comp, 1),
     1                          u ( ixxm1, iyy, izz, i_comp, 1),
     1                          u ( ixx, iyyp1, izz, i_comp, 1),
     1                          u ( ixxp1, iyyp1, izz, i_comp, 1),
     1                          u ( ixxm1, iyyp1, izz, i_comp, 1),
     1                          u ( ixx, iyym1, izz, i_comp, 1),
     1                          u ( ixxp1, iyym1, izz, i_comp, 1),
     1                          u ( ixxm1, iyym1, izz, i_comp, 1),

     1                          u ( ixx, iyy, izzp1, i_comp, 1),
     1                          u ( ixxp1, iyy, izzp1, i_comp, 1),
     1                          u ( ixxm1, iyy, izzp1, i_comp, 1),
     1                          u ( ixx, iyyp1, izzp1, i_comp, 1),
     1                          u ( ixxp1, iyyp1, izzp1, i_comp, 1),
     1                          u ( ixxm1, iyyp1, izzp1, i_comp, 1),
     1                          u ( ixx, iyym1, izzp1, i_comp, 1),
     1                          u ( ixxp1, iyym1, izzp1, i_comp, 1),
     1                          u ( ixxm1, iyym1, izzp1, i_comp, 1),

     1                          u ( ixx, iyy, izzm1, i_comp, 1),
     1                          u ( ixxp1, iyy, izzm1, i_comp, 1),
     1                          u ( ixxm1, iyy, izzm1, i_comp, 1),
     1                          u ( ixx, iyyp1, izzm1, i_comp, 1),
     1                          u ( ixxp1, iyyp1, izzm1, i_comp, 1),
     1                          u ( ixxm1, iyyp1, izzm1, i_comp, 1),
     1                          u ( ixx, iyym1, izzm1, i_comp, 1),
     1                          u ( ixxp1, iyym1, izzm1, i_comp, 1),
     1                          u ( ixxm1, iyym1, izzm1, i_comp, 1) )

! ----------

          END DO

! --------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          a ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          b ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          c ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          d ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          e ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          f ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          g ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 + 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          h ( ixx) =  u ( ixx, iyy, izz, i_comp, 1)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 2)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 3)
     1 - 0.5 * u ( ixx, iyy, izz, i_comp, 4)

          var_zone_min ( ixx) = AMIN1 ( a ( ixx), b ( ixx), c ( ixx),
     1               d ( ixx), e ( ixx), f ( ixx), g ( ixx), h ( ixx))

          var_zone_max ( ixx) = AMAX1 ( a ( ixx), b ( ixx), c ( ixx),
     1               d ( ixx), e ( ixx), f ( ixx), g ( ixx), h ( ixx))

! ----------

          b ( ixx) = ( var_neighbor_max ( ixx)
     1               - u ( ixx, iyy, izz, i_comp, 1) )
     1      / ( var_zone_max ( ixx) - u ( ixx, iyy, izz, i_comp, 1)
     1        + smallnum)

          c ( ixx) = ( u ( ixx, iyy, izz, i_comp, 1)
     1               - var_neighbor_min ( ixx) )
     1      / ( u ( ixx, iyy, izz, i_comp, 1) - var_zone_min ( ixx)
     1         + smallnum)

          a ( ixx) = AMIN1 ( 1.0, psi_limiter
     1                          * AMIN1 ( b ( ixx), c ( ixx)) )

          u ( ixx, iyy, izz, i_comp, 2) = a ( ixx)
     1                                  * u ( ixx, iyy, izz, i_comp, 2)
          u ( ixx, iyy, izz, i_comp, 3) = a ( ixx)
     1                                  * u ( ixx, iyy, izz, i_comp, 3)
          u ( ixx, iyy, izz, i_comp, 4) = a ( ixx)
     1                                  * u ( ixx, iyy, izz, i_comp, 4)

! ----------

          END DO ! End of DO ixx = loop.

! --------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------

      END IF ! Don't limit magnetic field components.

      END DO ! End of DO i_comp = loop.

! ----------------------------------------------------------------------
! If we are building facial slopes then we limit them multidimensionally too.
#if ( DIV_FREE_RECONSTRUCT == 1)
! ----------------------------------------------------------------------
! Apply multidimensional limiting to the x-magnetic fields in the x-faces.

      psi_limiter = 0.5    ! For tough problems, use this value for all vars.

!      psi_limiter = 1.0

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d,                                &
!$OMP&    var_zone_min, var_zone_max, var_neighbor_min,var_neighbor_max)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    smallnum, psi_limiter, bfxbdy)                                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 0 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)

! ----------

          a ( ixx) =  bfxbdy ( ixx, iyy, izz, 1)
     1        + 0.5 * bfxbdy ( ixx, iyy, izz, 2)
     1        + 0.5 * bfxbdy ( ixx, iyy, izz, 3)

          b ( ixx) =  bfxbdy ( ixx, iyy, izz, 1)
     1        - 0.5 * bfxbdy ( ixx, iyy, izz, 2)
     1        + 0.5 * bfxbdy ( ixx, iyy, izz, 3)

          c ( ixx) =  bfxbdy ( ixx, iyy, izz, 1)
     1        + 0.5 * bfxbdy ( ixx, iyy, izz, 2)
     1        - 0.5 * bfxbdy ( ixx, iyy, izz, 3)

          d ( ixx) =  bfxbdy ( ixx, iyy, izz, 1)
     1        - 0.5 * bfxbdy ( ixx, iyy, izz, 2)
     1        - 0.5 * bfxbdy ( ixx, iyy, izz, 3)

          var_zone_min ( ixx) = AMIN1 ( a ( ixx), b ( ixx), c ( ixx),
     1                                  d ( ixx))

          var_zone_max ( ixx) = AMAX1 ( a ( ixx), b ( ixx), c ( ixx),
     1                                  d ( ixx))

! ---------- 

          var_neighbor_min ( ixx) = AMIN1 ( bfxbdy ( ixx, iyy, izz, 1),
     1                          bfxbdy ( ixx, iyyp1, izz, 1),
     1                          bfxbdy ( ixx, iyym1, izz, 1),
     1                          bfxbdy ( ixx, iyy, izzp1, 1),
     1                          bfxbdy ( ixx, iyyp1, izzp1, 1),
     1                          bfxbdy ( ixx, iyym1, izzp1, 1),
     1                          bfxbdy ( ixx, iyy, izzm1, 1),
     1                          bfxbdy ( ixx, iyyp1, izzm1, 1),
     1                          bfxbdy ( ixx, iyym1, izzm1, 1))

          var_neighbor_max ( ixx) = AMAX1 ( bfxbdy ( ixx, iyy, izz, 1),
     1                          bfxbdy ( ixx, iyyp1, izz, 1),
     1                          bfxbdy ( ixx, iyym1, izz, 1),
     1                          bfxbdy ( ixx, iyy, izzp1, 1),
     1                          bfxbdy ( ixx, iyyp1, izzp1, 1),
     1                          bfxbdy ( ixx, iyym1, izzp1, 1),
     1                          bfxbdy ( ixx, iyy, izzm1, 1),
     1                          bfxbdy ( ixx, iyyp1, izzm1, 1),
     1                          bfxbdy ( ixx, iyym1, izzm1, 1))

! ----------

          b ( ixx) = ( var_neighbor_max ( ixx)
     1               - bfxbdy ( ixx, iyy, izz, 1) )
     1      / ( var_zone_max ( ixx)
     1        - bfxbdy ( ixx, iyy, izz, 1) + smallnum)

          c ( ixx) = ( bfxbdy ( ixx, iyy, izz, 1)
     1               - var_neighbor_min ( ixx) )
     1      / ( bfxbdy ( ixx, iyy, izz, 1)
     1        - var_zone_min ( ixx) + smallnum)

          a ( ixx) = AMIN1 ( 1.0, psi_limiter
     1                          * AMIN1 ( b ( ixx), c ( ixx)) )

          bfxbdy ( ixx, iyy, izz, 2) = a ( ixx)
     1                               * bfxbdy ( ixx, iyy, izz, 2)
          bfxbdy ( ixx, iyy, izz, 3) = a ( ixx)
     1                               * bfxbdy ( ixx, iyy, izz, 3)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Apply multidimensional limiting to the magnetic fields on the y-faces.

      psi_limiter = 0.5    ! For tough problems, use this value for all vars.

!      psi_limiter = 1.0

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d,                                &
!$OMP&    var_zone_min, var_zone_max, var_neighbor_min,var_neighbor_max)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    smallnum, psi_limiter, bfybdy)                                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 0 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 0 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

          a ( ixx) =  bfybdy ( ixx, iyy, izz, 1)
     1        + 0.5 * bfybdy ( ixx, iyy, izz, 2)
     1        + 0.5 * bfybdy ( ixx, iyy, izz, 3)

          b ( ixx) =  bfybdy ( ixx, iyy, izz, 1)
     1        - 0.5 * bfybdy ( ixx, iyy, izz, 2)
     1        + 0.5 * bfybdy ( ixx, iyy, izz, 3)

          c ( ixx) =  bfybdy ( ixx, iyy, izz, 1)
     1        + 0.5 * bfybdy ( ixx, iyy, izz, 2)
     1        - 0.5 * bfybdy ( ixx, iyy, izz, 3)

          d ( ixx) =  bfybdy ( ixx, iyy, izz, 1)
     1        - 0.5 * bfybdy ( ixx, iyy, izz, 2)
     1        - 0.5 * bfybdy ( ixx, iyy, izz, 3)

          var_zone_min ( ixx) = AMIN1 ( a ( ixx), b ( ixx), c ( ixx),
     1                                  d ( ixx))

          var_zone_max ( ixx) = AMAX1 ( a ( ixx), b ( ixx), c ( ixx),
     1                                  d ( ixx))

! ---------- 

          var_neighbor_min ( ixx) = AMIN1 ( bfybdy ( ixx, iyy, izz, 1),
     1                          bfybdy ( ixxp1, iyy, izz, 1),
     1                          bfybdy ( ixxm1, iyy, izz, 1),
     1                          bfybdy ( ixx, iyy, izzp1, 1),
     1                          bfybdy ( ixxp1, iyy, izzp1, 1),
     1                          bfybdy ( ixxm1, iyy, izzp1, 1),
     1                          bfybdy ( ixx, iyy, izzm1, 1),
     1                          bfybdy ( ixxp1, iyy, izzm1, 1),
     1                          bfybdy ( ixxm1, iyy, izzm1, 1))

          var_neighbor_max ( ixx) = AMAX1 ( bfybdy ( ixx, iyy, izz, 1),
     1                          bfybdy ( ixxp1, iyy, izz, 1),
     1                          bfybdy ( ixxm1, iyy, izz, 1),
     1                          bfybdy ( ixx, iyy, izzp1, 1),
     1                          bfybdy ( ixxp1, iyy, izzp1, 1),
     1                          bfybdy ( ixxm1, iyy, izzp1, 1),
     1                          bfybdy ( ixx, iyy, izzm1, 1),
     1                          bfybdy ( ixxp1, iyy, izzm1, 1),
     1                          bfybdy ( ixxm1, iyy, izzm1, 1))

! ----------

          b ( ixx) = ( var_neighbor_max ( ixx)
     1               - bfybdy ( ixx, iyy, izz, 1) )
     1      / ( var_zone_max ( ixx)
     1        - bfybdy ( ixx, iyy, izz, 1) + smallnum)

          c ( ixx) = ( bfybdy ( ixx, iyy, izz, 1)
     1               - var_neighbor_min ( ixx) )
     1      / ( bfybdy ( ixx, iyy, izz, 1)
     1        - var_zone_min ( ixx) + smallnum)

          a ( ixx) = AMIN1 ( 1.0, psi_limiter
     1                          * AMIN1 ( b ( ixx), c ( ixx)) )

          bfybdy ( ixx, iyy, izz, 2) = a ( ixx)
     1                               * bfybdy ( ixx, iyy, izz, 2)
          bfybdy ( ixx, iyy, izz, 3) = a ( ixx)
     1                               * bfybdy ( ixx, iyy, izz, 3)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Apply multidimensional limiting to the magnetic fields on the z-faces.

      psi_limiter = 0.5    ! For tough problems, use this value for all vars.

!      psi_limiter = 1.0

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d,                                &
!$OMP&    var_zone_min, var_zone_max, var_neighbor_min,var_neighbor_max)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    smallnum, psi_limiter, bfzbdy)                                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 0 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

          a ( ixx) =  bfzbdy ( ixx, iyy, izz, 1)
     1        + 0.5 * bfzbdy ( ixx, iyy, izz, 2)
     1        + 0.5 * bfzbdy ( ixx, iyy, izz, 3)

          b ( ixx) =  bfzbdy ( ixx, iyy, izz, 1)
     1        - 0.5 * bfzbdy ( ixx, iyy, izz, 2)
     1        + 0.5 * bfzbdy ( ixx, iyy, izz, 3)

          c ( ixx) =  bfzbdy ( ixx, iyy, izz, 1)
     1        + 0.5 * bfzbdy ( ixx, iyy, izz, 2)
     1        - 0.5 * bfzbdy ( ixx, iyy, izz, 3)

          d ( ixx) =  bfzbdy ( ixx, iyy, izz, 1)
     1        - 0.5 * bfzbdy ( ixx, iyy, izz, 2)
     1        - 0.5 * bfzbdy ( ixx, iyy, izz, 3)

          var_zone_min ( ixx) = AMIN1 ( a ( ixx), b ( ixx), c ( ixx),
     1                                  d ( ixx))

          var_zone_max ( ixx) = AMAX1 ( a ( ixx), b ( ixx), c ( ixx),
     1                                  d ( ixx))

! ---------- 

          var_neighbor_min ( ixx) = AMIN1 ( bfzbdy ( ixx, iyy, izz, 1),
     1                          bfzbdy ( ixxp1, iyy, izz, 1),
     1                          bfzbdy ( ixxm1, iyy, izz, 1),
     1                          bfzbdy ( ixx, iyyp1, izz, 1),
     1                          bfzbdy ( ixxp1, iyyp1, izz, 1),
     1                          bfzbdy ( ixxm1, iyyp1, izz, 1),
     1                          bfzbdy ( ixx, iyym1, izz, 1),
     1                          bfzbdy ( ixxp1, iyym1, izz, 1),
     1                          bfzbdy ( ixxm1, iyym1, izz, 1))

          var_neighbor_max ( ixx) = AMAX1 ( bfzbdy ( ixx, iyy, izz, 1),
     1                          bfzbdy ( ixxp1, iyy, izz, 1),
     1                          bfzbdy ( ixxm1, iyy, izz, 1),
     1                          bfzbdy ( ixx, iyyp1, izz, 1),
     1                          bfzbdy ( ixxp1, iyyp1, izz, 1),
     1                          bfzbdy ( ixxm1, iyyp1, izz, 1),
     1                          bfzbdy ( ixx, iyym1, izz, 1),
     1                          bfzbdy ( ixxp1, iyym1, izz, 1),
     1                          bfzbdy ( ixxm1, iyym1, izz, 1))

! ----------

          b ( ixx) = ( var_neighbor_max ( ixx)
     1               - bfzbdy ( ixx, iyy, izz, 1) )
     1      / ( var_zone_max ( ixx)
     1        - bfzbdy ( ixx, iyy, izz, 1) + smallnum)

          c ( ixx) = ( bfzbdy ( ixx, iyy, izz, 1)
     1               - var_neighbor_min ( ixx) )
     1      / ( bfzbdy ( ixx, iyy, izz, 1)
     1        - var_zone_min ( ixx) + smallnum)

          a ( ixx) = AMIN1 ( 1.0, psi_limiter
     1                          * AMIN1 ( b ( ixx), c ( ixx)) )

          bfzbdy ( ixx, iyy, izz, 2) = a ( ixx)
     1                               * bfzbdy ( ixx, iyy, izz, 2)
          bfzbdy ( ixx, iyy, izz, 3) = a ( ixx)
     1                               * bfzbdy ( ixx, iyy, izz, 3)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------

! ----------------------------------------------------------------------
#endif
! End of "DIV_FREE_RECONSTRUCT == 1".
! ----------------------------------------------------------------------
#endif
! End of "MULTID_LIMITING == 1"
! ----------------------------------------------------------------------
! Multidimensional limiting ends here.
! ----------------------------------------------------------------------
! For problems with very strong shocks and large magnetic fields, one
! may have a loss of pressure positivity. The entropy formulation cures
! that and it is invoked by switching on the pressure positivity switch.
! It does cause a local loss of conservation in the zone where it is invoked.

! The gradients of the entropy density are built here. This is done with
! lower quality reconstruction because it is also more stable.

#if ( PRESSURE_POSITIVITY == 1)

! ----------------------------------------------------------------------

      psi_limiter = 0.5    ! For tough problems, use this value for entropy.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d, e, f, g, h,                    &
!$OMP&    var_zone_min, var_zone_max, var_neighbor_min,var_neighbor_max)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    mc_coef, smallnum, psi_limiter, entropygr, flattengr)         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

! --------------------
! Loop split just to help with vectorization.

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

          a ( ixx) = entropygr ( ixxp1, iyy, izz, 1)
     1             - entropygr ( ixx, iyy, izz, 1)

          b ( ixx) = entropygr ( ixx, iyy, izz, 1)
     1             - entropygr ( ixxm1, iyy, izz, 1)

! ----------

          var_neighbor_min ( ixx) =
     1                  AMIN1 ( entropygr ( ixx, iyy, izz, 1),
     1                          entropygr ( ixxp1, iyy, izz, 1),
     1                          entropygr ( ixxm1, iyy, izz, 1),
     1                          entropygr ( ixx, iyyp1, izz, 1),
     1                          entropygr ( ixxp1, iyyp1, izz, 1),
     1                          entropygr ( ixxm1, iyyp1, izz, 1),
     1                          entropygr ( ixx, iyym1, izz, 1),
     1                          entropygr ( ixxp1, iyym1, izz, 1),
     1                          entropygr ( ixxm1, iyym1, izz, 1),

     1                          entropygr ( ixx, iyy, izzp1, 1),
     1                          entropygr ( ixxp1, iyy, izzp1, 1),
     1                          entropygr ( ixxm1, iyy, izzp1, 1),
     1                          entropygr ( ixx, iyyp1, izzp1, 1),
     1                          entropygr ( ixxp1, iyyp1, izzp1, 1),
     1                          entropygr ( ixxm1, iyyp1, izzp1, 1),
     1                          entropygr ( ixx, iyym1, izzp1, 1),
     1                          entropygr ( ixxp1, iyym1, izzp1, 1),
     1                          entropygr ( ixxm1, iyym1, izzp1, 1),

     1                          entropygr ( ixx, iyy, izzm1, 1),
     1                          entropygr ( ixxp1, iyy, izzm1, 1),
     1                          entropygr ( ixxm1, iyy, izzm1, 1),
     1                          entropygr ( ixx, iyyp1, izzm1, 1),
     1                          entropygr ( ixxp1, iyyp1, izzm1, 1),
     1                          entropygr ( ixxm1, iyyp1, izzm1, 1),
     1                          entropygr ( ixx, iyym1, izzm1, 1),
     1                          entropygr ( ixxp1, iyym1, izzm1, 1),
     1                          entropygr ( ixxm1, iyym1, izzm1, 1) )

          var_neighbor_max ( ixx) =
     1                  AMAX1 ( entropygr ( ixx, iyy, izz, 1),
     1                          entropygr ( ixxp1, iyy, izz, 1),
     1                          entropygr ( ixxm1, iyy, izz, 1),
     1                          entropygr ( ixx, iyyp1, izz, 1),
     1                          entropygr ( ixxp1, iyyp1, izz, 1),
     1                          entropygr ( ixxm1, iyyp1, izz, 1),
     1                          entropygr ( ixx, iyym1, izz, 1),
     1                          entropygr ( ixxp1, iyym1, izz, 1),
     1                          entropygr ( ixxm1, iyym1, izz, 1),

     1                          entropygr ( ixx, iyy, izzp1, 1),
     1                          entropygr ( ixxp1, iyy, izzp1, 1),
     1                          entropygr ( ixxm1, iyy, izzp1, 1),
     1                          entropygr ( ixx, iyyp1, izzp1, 1),
     1                          entropygr ( ixxp1, iyyp1, izzp1, 1),
     1                          entropygr ( ixxm1, iyyp1, izzp1, 1),
     1                          entropygr ( ixx, iyym1, izzp1, 1),
     1                          entropygr ( ixxp1, iyym1, izzp1, 1),
     1                          entropygr ( ixxm1, iyym1, izzp1, 1),

     1                          entropygr ( ixx, iyy, izzm1, 1),
     1                          entropygr ( ixxp1, iyy, izzm1, 1),
     1                          entropygr ( ixxm1, iyy, izzm1, 1),
     1                          entropygr ( ixx, iyyp1, izzm1, 1),
     1                          entropygr ( ixxp1, iyyp1, izzm1, 1),
     1                          entropygr ( ixxm1, iyyp1, izzm1, 1),
     1                          entropygr ( ixx, iyym1, izzm1, 1),
     1                          entropygr ( ixxp1, iyym1, izzm1, 1),
     1                          entropygr ( ixxm1, iyym1, izzm1, 1) )

! ----------

          END DO

! --------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          c ( ixx) = 1.0 - flattengr ( ixx, iyy, izz)

! ----------

          entropygr ( ixx, iyy, izz, 2) =
     1         MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = entropygr ( ixx, iyyp1, izz, 1)
     1             - entropygr ( ixx, iyy, izz, 1)

          b ( ixx) = entropygr ( ixx, iyy, izz, 1)
     1             - entropygr ( ixx, iyym1, izz, 1)

          entropygr ( ixx, iyy, izz, 3) =
     1         MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) = entropygr ( ixx, iyy, izzp1, 1)
     1             - entropygr ( ixx, iyy, izz, 1)

          b ( ixx) = entropygr ( ixx, iyy, izz, 1)
     1             - entropygr ( ixx, iyy, izzm1, 1)

          entropygr ( ixx, iyy, izz, 4) =
     1         MC_beta ( a ( ixx), b ( ixx), mc_coef) * c ( ixx)

! ----------

          a ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          b ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          c ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          d ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          e ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          f ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          g ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          h ( ixx) =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          var_zone_min ( ixx) = AMIN1 ( a ( ixx), b ( ixx), c ( ixx),
     1               d ( ixx), e ( ixx), f ( ixx), g ( ixx), h ( ixx))

          var_zone_max ( ixx) = AMAX1 ( a ( ixx), b ( ixx), c ( ixx),
     1               d ( ixx), e ( ixx), f ( ixx), g ( ixx), h ( ixx))

! ----------

          b ( ixx) = ( var_neighbor_max ( ixx)
     1               - entropygr ( ixx, iyy, izz, 1) )
     1      / ( var_zone_max ( ixx) - entropygr ( ixx, iyy, izz, 1)
     1        + smallnum)

          c ( ixx) = ( entropygr ( ixx, iyy, izz, 1)
     1               - var_neighbor_min ( ixx) )
     1      / ( entropygr ( ixx, iyy, izz, 1) - var_zone_min ( ixx)
     1        + smallnum)

          a ( ixx) = AMIN1 ( 1.0, psi_limiter
     1                          * AMIN1 ( b ( ixx), c ( ixx)) )

          entropygr ( ixx, iyy, izz, 2) =
     1    a ( ixx) * entropygr ( ixx, iyy, izz, 2)
          entropygr ( ixx, iyy, izz, 3) =
     1    a ( ixx) * entropygr ( ixx, iyy, izz, 3)
          entropygr ( ixx, iyy, izz, 4) =
     1    a ( ixx) * entropygr ( ixx, iyy, izz, 4)

! ----------

          END DO ! End of DO ixx = loop.

! --------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
#endif
! End of "PRESSURE_POSITIVITY == 1".
! ----------------------------------------------------------------------

! Now make the divergence-free reconstruction within each zone.
! Use that to obtain the mean magnetic field within the zone as
! well as the slopes of the magnetic field in all three directions
! at the zone centers. Also use that to obtain the zone-centered
! pressure.

#if ( DIV_FREE_RECONSTRUCT == 1)

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, dx, dy, dz, ob_dx, ob_dy, ob_dz, a, b, c,  &
!$OMP&    bfx_pl_0, bfx_pl_y, bfx_pl_z, bfx_mi_0, bfx_mi_y, bfx_mi_z,   &
!$OMP&    bfy_pl_0, bfy_pl_x, bfy_pl_z, bfy_mi_0, bfy_mi_x, bfy_mi_z,   &
!$OMP&    bfz_pl_0, bfz_pl_x, bfz_pl_y, bfz_mi_0, bfz_mi_x, bfz_mi_y)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, prsfloor, ob_8pi, ob_6, xb_save, dx_save, xc_save,     &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,       &
!$OMP&    gamma_eos_gr)                                                 &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzm1 = izz - 1

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyym1 = iyy - 1

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxm1 = ixx - 1

! ----------

          IF ( igeom .EQ. 1) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 2) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 3) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1              * dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

! ----------

          bfx_pl_0 ( ixx) = bfxbdy ( ixx, iyy, izz, 1) * ob_dx ( ixx)
          bfx_pl_y ( ixx) = bfxbdy ( ixx, iyy, izz, 2) * ob_dx ( ixx)
          bfx_pl_z ( ixx) = bfxbdy ( ixx, iyy, izz, 3) * ob_dx ( ixx)

          bfx_mi_0 ( ixx) = bfxbdy ( ixxm1, iyy, izz, 1) * ob_dx ( ixx)
          bfx_mi_y ( ixx) = bfxbdy ( ixxm1, iyy, izz, 2) * ob_dx ( ixx)
          bfx_mi_z ( ixx) = bfxbdy ( ixxm1, iyy, izz, 3) * ob_dx ( ixx)

          bfy_pl_0 ( ixx) = bfybdy ( ixx, iyy, izz, 1) * ob_dy ( ixx)
          bfy_pl_x ( ixx) = bfybdy ( ixx, iyy, izz, 2) * ob_dy ( ixx)
          bfy_pl_z ( ixx) = bfybdy ( ixx, iyy, izz, 3) * ob_dy ( ixx)

          bfy_mi_0 ( ixx) = bfybdy ( ixx, iyym1, izz, 1) * ob_dy ( ixx)
          bfy_mi_x ( ixx) = bfybdy ( ixx, iyym1, izz, 2) * ob_dy ( ixx)
          bfy_mi_z ( ixx) = bfybdy ( ixx, iyym1, izz, 3) * ob_dy ( ixx)

          bfz_pl_0 ( ixx) = bfzbdy ( ixx, iyy, izz, 1) * ob_dz ( ixx)
          bfz_pl_x ( ixx) = bfzbdy ( ixx, iyy, izz, 2) * ob_dz ( ixx)
          bfz_pl_y ( ixx) = bfzbdy ( ixx, iyy, izz, 3) * ob_dz ( ixx)

          bfz_mi_0 ( ixx) = bfzbdy ( ixx, iyy, izzm1, 1) * ob_dz ( ixx)
          bfz_mi_x ( ixx) = bfzbdy ( ixx, iyy, izzm1, 2) * ob_dz ( ixx)
          bfz_mi_y ( ixx) = bfzbdy ( ixx, iyy, izzm1, 3) * ob_dz ( ixx)

! ----------

          bfxccmodes ( ixx, iyy, izz, 3) =
     1               0.5 * ( bfx_pl_y ( ixx) + bfx_mi_y ( ixx) )
          bfxccmodes ( ixx, iyy, izz, 6) =
     1               bfx_pl_y ( ixx) - bfx_mi_y ( ixx)
          bfxccmodes ( ixx, iyy, izz, 4) =
     1               0.5 * ( bfx_pl_z ( ixx) + bfx_mi_z ( ixx) )
          bfxccmodes ( ixx, iyy, izz, 7) =
     1               bfx_pl_z ( ixx) - bfx_mi_z ( ixx)

          bfyccmodes ( ixx, iyy, izz, 2) =
     1               0.5 * ( bfy_pl_x ( ixx) + bfy_mi_x ( ixx) )
          bfyccmodes ( ixx, iyy, izz, 6) =
     1               bfy_pl_x ( ixx) - bfy_mi_x ( ixx)
          bfyccmodes ( ixx, iyy, izz, 4) =
     1               0.5 * ( bfy_pl_z ( ixx) + bfy_mi_z ( ixx) )
          bfyccmodes ( ixx, iyy, izz, 7) =
     1               bfy_pl_z ( ixx) - bfy_mi_z ( ixx)

          bfzccmodes ( ixx, iyy, izz, 2) =
     1               0.5 * ( bfz_pl_x ( ixx) + bfz_mi_x ( ixx) )
          bfzccmodes ( ixx, iyy, izz, 7) =
     1               bfz_pl_x ( ixx) - bfz_mi_x ( ixx)
          bfzccmodes ( ixx, iyy, izz, 3) =
     1               0.5 * ( bfz_pl_y ( ixx) + bfz_mi_y ( ixx) )
          bfzccmodes ( ixx, iyy, izz, 6) =
     1               bfz_pl_y ( ixx) - bfz_mi_y ( ixx)

          bfxccmodes ( ixx, iyy, izz, 5) =
     1             - 0.5 * ( bfyccmodes ( ixx, iyy, izz, 6)
     1                     + bfzccmodes ( ixx, iyy, izz, 7) )
          bfyccmodes ( ixx, iyy, izz, 5) =
     1             - 0.5 * ( bfxccmodes ( ixx, iyy, izz, 6)
     1                     + bfzccmodes ( ixx, iyy, izz, 6) )
          bfzccmodes ( ixx, iyy, izz, 5) =
     1             - 0.5 * ( bfxccmodes ( ixx, iyy, izz, 7)
     1                     + bfyccmodes ( ixx, iyy, izz, 7) ) 

          bfxccmodes ( ixx, iyy, izz, 1) =
     1               0.5 * ( bfx_pl_0 ( ixx) + bfx_mi_0 ( ixx) )
     1             - ob_6 * bfxccmodes ( ixx, iyy, izz, 5)
          bfxccmodes ( ixx, iyy, izz, 2) =
     1               bfx_pl_0 ( ixx) - bfx_mi_0 ( ixx)
          bfyccmodes ( ixx, iyy, izz, 1) =
     1               0.5 * ( bfy_pl_0 ( ixx) + bfy_mi_0 ( ixx) )
     1             - ob_6 * bfyccmodes ( ixx, iyy, izz, 5)
          bfyccmodes ( ixx, iyy, izz, 3) =
     1               bfy_pl_0 ( ixx) - bfy_mi_0 ( ixx)
          bfzccmodes ( ixx, iyy, izz, 1) =
     1               0.5 * ( bfz_pl_0 ( ixx) + bfz_mi_0 ( ixx) )
     1             - ob_6 * bfzccmodes ( ixx, iyy, izz, 5)
          bfzccmodes ( ixx, iyy, izz, 4) =
     1               bfz_pl_0 ( ixx) - bfz_mi_0 ( ixx)

! ----------

          bfxccmodes ( ixx, iyy, izz, :) =
     1               bfxccmodes ( ixx, iyy, izz, :) * dx ( ixx)

          bfyccmodes ( ixx, iyy, izz, :) =
     1               bfyccmodes ( ixx, iyy, izz, :) * dy ( ixx)

          bfzccmodes ( ixx, iyy, izz, :) =
     1               bfzccmodes ( ixx, iyy, izz, :) * dz ( ixx)

! ----------

          u ( ixx, iyy, izz, 6, 1:4) = bfxccmodes ( ixx, iyy, izz, 1:4)

          u ( ixx, iyy, izz, 7, 1:4) = bfyccmodes ( ixx, iyy, izz, 1:4)

          u ( ixx, iyy, izz, 8, 1:4) = bfzccmodes ( ixx, iyy, izz, 1:4)

! ----------

          bfxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 6, 1)

          bfygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 7, 1)

          bfzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 8, 1)

          b ( ixx) = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          c ( ixx) = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          prsgr ( ixx, iyy, izz) =
     1          ( u ( ixx, iyy, izz, 5, 1) - b ( ixx) - c ( ixx) )
     1                         * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

          IF ( prsgr ( ixx, iyy, izz) .LT. prsfloor) THEN
          prsgr ( ixx, iyy, izz) = prsfloor
          u ( ixx, iyy, izz, 5, 1) = b ( ixx) + c ( ixx) + prsfloor
     1                         / ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
          END IF

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#else
! ------------------------------

! Here the magnetic field slopes in a zone are built using just the
! zone-centered variables in the transverse directions. The longitudinal
! slopes are still built here in a consistent fashion, i.e. at the
! face-centers the normal components match up from either side.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1)                                            &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    bfxgr, bfygr, bfzgr)                                          &
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

! ----------
! Make longitudinal slopes

          bfxccmodes ( ixx, iyy, izz, 2) =
     1                               ( bfxbdy ( ixx, iyy, izz, 1)
     1                               - bfxbdy ( ixxm1, iyy, izz, 1))
          bfxccmodes ( ixx, iyy, izz, 5:7) = 0.0

          bfyccmodes ( ixx, iyy, izz, 3) =
     1                               ( bfybdy ( ixx, iyy, izz, 1)
     1                               - bfybdy ( ixx, iyym1, izz, 1))
          bfyccmodes ( ixx, iyy, izz, 5:7) = 0.0

          bfzccmodes ( ixx, iyy, izz, 4) =
     1                               ( bfzbdy ( ixx, iyy, izz, 1)
     1                               - bfzbdy ( ixx, iyy, izzm1, 1))
          bfzccmodes ( ixx, iyy, izz, 5:7) = 0.0

! ----------
! Restore consistency between "u" and "bfxccmodes, bfyccmodes, bfzccmodes".

          u ( ixx, iyy, izz, 6, 2) = bfxccmodes ( ixx, iyy, izz, 2)
          bfxccmodes ( ixx, iyy, izz, 1) = u ( ixx, iyy, izz, 6, 1)
          bfxccmodes ( ixx, iyy, izz, 3) = u ( ixx, iyy, izz, 6, 3)
          bfxccmodes ( ixx, iyy, izz, 4) = u ( ixx, iyy, izz, 6, 4)

          u ( ixx, iyy, izz, 7, 3) = bfyccmodes ( ixx, iyy, izz, 3)
          bfyccmodes ( ixx, iyy, izz, 1) = u ( ixx, iyy, izz, 7, 1)
          bfyccmodes ( ixx, iyy, izz, 2) = u ( ixx, iyy, izz, 7, 2)
          bfyccmodes ( ixx, iyy, izz, 4) = u ( ixx, iyy, izz, 7, 4)

          u ( ixx, iyy, izz, 8, 4) = bfzccmodes ( ixx, iyy, izz, 4)
          bfzccmodes ( ixx, iyy, izz, 1) = u ( ixx, iyy, izz, 8, 1)
          bfzccmodes ( ixx, iyy, izz, 2) = u ( ixx, iyy, izz, 8, 2)
          bfzccmodes ( ixx, iyy, izz, 3) = u ( ixx, iyy, izz, 8, 3)

! ----------
! This may be redundant, but restore consistency for zone-averaged fields.

          bfxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 6, 1)

          bfygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 7, 1)

          bfzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 8, 1)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
! Average the zone-centered slopes to the face-centered slopes at x-faces.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixx0, ixxp1, ixxm1)                                      &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    u, bfxbdy)                                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 0 - ioffx, ix1 + ioffx
          ixx0 = MAX0 ( ixx, 1 - ioffx)
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          bfxbdy ( ixx, iyy, izz, 2) =
     1           0.5 * ( u ( ixx0, iyy, izz, 6, 3)
     1                 + u ( ixxp1, iyy, izz, 6, 3) )

          bfxbdy ( ixx, iyy, izz, 3) =
     1           0.5 * ( u ( ixx0, iyy, izz, 6, 4)
     1                 + u ( ixxp1, iyy, izz, 6, 4) )

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
! Average the zone-centered slopes to the face-centered slopes at y-faces.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyy0, iyyp1, iyym1, &
!$OMP&    ixx, ixxp1, ixxm1)                                            &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    u, bfybdy)                                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 0 - ioffy, iy1 + ioffy
        iyy0 = MAX0 ( iyy, 1 - ioffy)
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          bfybdy ( ixx, iyy, izz, 2) =
     1           0.5 * ( u ( ixx, iyy0, izz, 7, 2)
     1                 + u ( ixx, iyyp1, izz, 7, 2) )

          bfybdy ( ixx, iyy, izz, 3) =
     1           0.5 * ( u ( ixx, iyy0, izz, 7, 4)
     1                 + u ( ixx, iyyp1, izz, 7, 4) )

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
! Average the zone-centered slopes to the face-centered slopes at z-faces.

!$OMP PARALLEL DO PRIVATE ( izz, izz0, izzp1, izzm1, iyy, iyyp1, iyym1, &
!$OMP&    ixx, ixxp1, ixxm1)                                            &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    u, bfzbdy)                                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz, iz1 + ioffz
      izz0 = MAX0 ( izz, 1 - ioffz)
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          bfzbdy ( ixx, iyy, izz, 2) =
     1           0.5 * ( u ( ixx, iyy, izz0, 8, 2)
     1                 + u ( ixx, iyy, izzp1, 8, 2) )

          bfzbdy ( ixx, iyy, izz, 3) =
     1           0.5 * ( u ( ixx, iyy, izz0, 8, 3)
     1                 + u ( ixx, iyy, izzp1, 8, 3) )

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------
#endif
! End of "DIV_FREE_RECONSTRUCT == 1".
! ----------------------------------------------------------------------
! If we additionally want limiting on the characteristic variables then we
! do the below three loop nests.
#if ( LIMIT_TWICE == 1)
! ----------------------------------------------------------------------

! Characteristic limiting in the x-direction.

      IF ( ioffx .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,            &
!$OMP&    bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma, smallnum, prsfloor, rhofloor, pi, gasconst,            &
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, flattengr)                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx

          u_1d ( ixx, 1, 1) = u ( ixx, iyy, izz, 1, 1)
          u_1d ( ixx, 2, 1) = u ( ixx, iyy, izz, 2, 1)
          u_1d ( ixx, 3, 1) = u ( ixx, iyy, izz, 3, 1)
          u_1d ( ixx, 4, 1) = u ( ixx, iyy, izz, 4, 1)
          u_1d ( ixx, 5, 1) = u ( ixx, iyy, izz, 5, 1)
          u_1d ( ixx, 6, 1) = u ( ixx, iyy, izz, 7, 1)
          u_1d ( ixx, 7, 1) = u ( ixx, iyy, izz, 8, 1)
#if ( NFLUID >= 1)
          u_1d ( ixx, 8: NFIELD + NFLUID, 1) =
     1                       u ( ixx, iyy, izz, 9: n_cc_components, 1)
#endif

          u_1d ( ixx, 1, 2) = u ( ixx, iyy, izz, 1, 2)
          u_1d ( ixx, 2, 2) = u ( ixx, iyy, izz, 2, 2)
          u_1d ( ixx, 3, 2) = u ( ixx, iyy, izz, 3, 2)
          u_1d ( ixx, 4, 2) = u ( ixx, iyy, izz, 4, 2)
          u_1d ( ixx, 5, 2) = u ( ixx, iyy, izz, 5, 2)
          u_1d ( ixx, 6, 2) = u ( ixx, iyy, izz, 7, 2)
          u_1d ( ixx, 7, 2) = u ( ixx, iyy, izz, 8, 2)
#if ( NFLUID >= 1)
          u_1d ( ixx, 8: NFIELD + NFLUID, 2) = 
     1                       u ( ixx, iyy, izz, 9: n_cc_components, 2)
#endif

          rhobar ( ixx) = rhogr ( ixx, iyy, izz)
          prsbar ( ixx) = prsgr ( ixx, iyy, izz)

          vlxbar ( ixx) = vlxgr ( ixx, iyy, izz)
          vlybar ( ixx) = vlygr ( ixx, iyy, izz)
          vlzbar ( ixx) = vlzgr ( ixx, iyy, izz)

          bfxbar ( ixx) = bfxgr ( ixx, iyy, izz)
          bfybar ( ixx) = bfygr ( ixx, iyy, izz)
          bfzbar ( ixx) = bfzgr ( ixx, iyy, izz)

#if ( NFLUID >= 1)
          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( 1: NFLUID, ixx) =
     1           u ( ixx, iyy, izz, 9: n_cc_components, 1) * a ( ixx)
#endif

          flatten_1d ( ixx) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO ixx = loop.

! ----------

          bfxbdy_1d ( 0 - ioffx: ix1 + ioffx) =
     1    bfxbdy ( 0 - ioffx: ix1 + ioffx, iyy, izz, 1)

          CALL LIMIT_1D_MHD_CONS ( 1, ix1, ioffx,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx

          u ( ixx, iyy, izz, 1, 2) = u_1d ( ixx, 1, 2)
          u ( ixx, iyy, izz, 2, 2) = u_1d ( ixx, 2, 2)
          u ( ixx, iyy, izz, 3, 2) = u_1d ( ixx, 3, 2)
          u ( ixx, iyy, izz, 4, 2) = u_1d ( ixx, 4, 2)
          u ( ixx, iyy, izz, 5, 2) = u_1d ( ixx, 5, 2)
          u ( ixx, iyy, izz, 7, 2) = u_1d ( ixx, 6, 2)
          u ( ixx, iyy, izz, 8, 2) = u_1d ( ixx, 7, 2)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 2) =
     1                       u_1d ( ixx, 8: NFIELD + NFLUID, 2)
#endif

          END DO ! End of DO ixx = loop.

! ----------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "ioffx .GT. 0".

! ----------------------------------------------------------------------

! Characteristic limiting in the y-direction.

      IF ( ioffy .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,            &
!$OMP&    bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma, smallnum, prsfloor, rhofloor, pi, gasconst,            &
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, flattengr)                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

        DO ixx = 1 - ioffx, ix1 + ioffx

! ----------

          DO iyy = 1 - ioffy, iy1 + ioffy

          u_1d ( iyy, 1, 1) = u ( ixx, iyy, izz, 1, 1)
          u_1d ( iyy, 2, 1) = u ( ixx, iyy, izz, 3, 1)
          u_1d ( iyy, 3, 1) = u ( ixx, iyy, izz, 4, 1)
          u_1d ( iyy, 4, 1) = u ( ixx, iyy, izz, 2, 1)
          u_1d ( iyy, 5, 1) = u ( ixx, iyy, izz, 5, 1)
          u_1d ( iyy, 6, 1) = u ( ixx, iyy, izz, 8, 1)
          u_1d ( iyy, 7, 1) = u ( ixx, iyy, izz, 6, 1)
#if ( NFLUID >= 1)
          u_1d ( iyy, 8: NFIELD + NFLUID, 1) =
     1                       u ( ixx, iyy, izz, 9: n_cc_components, 1)
#endif

          u_1d ( iyy, 1, 2) = u ( ixx, iyy, izz, 1, 3)
          u_1d ( iyy, 2, 2) = u ( ixx, iyy, izz, 3, 3)
          u_1d ( iyy, 3, 2) = u ( ixx, iyy, izz, 4, 3)
          u_1d ( iyy, 4, 2) = u ( ixx, iyy, izz, 2, 3)
          u_1d ( iyy, 5, 2) = u ( ixx, iyy, izz, 5, 3)
          u_1d ( iyy, 6, 2) = u ( ixx, iyy, izz, 8, 3)
          u_1d ( iyy, 7, 2) = u ( ixx, iyy, izz, 6, 3)
#if ( NFLUID >= 1)
          u_1d ( iyy, 8: NFIELD + NFLUID, 2) =
     1                       u ( ixx, iyy, izz, 9: n_cc_components, 3)
#endif

          rhobar ( iyy) = rhogr ( ixx, iyy, izz)
          prsbar ( iyy) = prsgr ( ixx, iyy, izz)

          vlxbar ( iyy) = vlygr ( ixx, iyy, izz)
          vlybar ( iyy) = vlzgr ( ixx, iyy, izz)
          vlzbar ( iyy) = vlxgr ( ixx, iyy, izz)

          bfxbar ( iyy) = bfygr ( ixx, iyy, izz)
          bfybar ( iyy) = bfzgr ( ixx, iyy, izz)
          bfzbar ( iyy) = bfxgr ( ixx, iyy, izz)

#if ( NFLUID >= 1)
          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( 1: NFLUID, iyy) =
     1           u ( ixx, iyy, izz, 9: n_cc_components, 1) * a ( ixx)
#endif

          flatten_1d ( iyy) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO iyy = loop.

! ----------

          bfxbdy_1d ( 0 - ioffy: iy1 + ioffy) =
     1    bfybdy ( ixx, 0 - ioffy: iy1 + ioffy, izz, 1)

          CALL LIMIT_1D_MHD_CONS ( 1, iy1, ioffy,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------

          DO iyy = 1 - ioffy, iy1 + ioffy

          u ( ixx, iyy, izz, 1, 3) = u_1d ( iyy, 1, 2)
          u ( ixx, iyy, izz, 3, 3) = u_1d ( iyy, 2, 2)
          u ( ixx, iyy, izz, 4, 3) = u_1d ( iyy, 3, 2)
          u ( ixx, iyy, izz, 2, 3) = u_1d ( iyy, 4, 2)
          u ( ixx, iyy, izz, 5, 3) = u_1d ( iyy, 5, 2)
          u ( ixx, iyy, izz, 8, 3) = u_1d ( iyy, 6, 2)
          u ( ixx, iyy, izz, 6, 3) = u_1d ( iyy, 7, 2)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 3) =
     1                       u_1d ( iyy, 8: NFIELD + NFLUID, 2)
#endif

          END DO ! End of DO iyy = loop.

! ----------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "ioffy .GT. 0".

! ----------------------------------------------------------------------

! Characteristic limiting in the z-direction.

      IF ( ioffz .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,            &
!$OMP&    bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma, smallnum, prsfloor, rhofloor, pi, gasconst,            &
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, flattengr)                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO iyy = 1 - ioffy, iy1 + ioffy

! ----------

          DO izz = 1 - ioffz, iz1 + ioffz

          u_1d ( izz, 1, 1) = u ( ixx, iyy, izz, 1, 1)
          u_1d ( izz, 2, 1) = u ( ixx, iyy, izz, 4, 1)
          u_1d ( izz, 3, 1) = u ( ixx, iyy, izz, 2, 1)
          u_1d ( izz, 4, 1) = u ( ixx, iyy, izz, 3, 1)
          u_1d ( izz, 5, 1) = u ( ixx, iyy, izz, 5, 1)
          u_1d ( izz, 6, 1) = u ( ixx, iyy, izz, 6, 1)
          u_1d ( izz, 7, 1) = u ( ixx, iyy, izz, 7, 1)
#if ( NFLUID >= 1)
          u_1d ( izz, 8: NFIELD + NFLUID, 1) =
     1                       u ( ixx, iyy, izz, 9: n_cc_components, 1)
#endif

          u_1d ( izz, 1, 2) = u ( ixx, iyy, izz, 1, 4)
          u_1d ( izz, 2, 2) = u ( ixx, iyy, izz, 4, 4)
          u_1d ( izz, 3, 2) = u ( ixx, iyy, izz, 2, 4)
          u_1d ( izz, 4, 2) = u ( ixx, iyy, izz, 3, 4)
          u_1d ( izz, 5, 2) = u ( ixx, iyy, izz, 5, 4)
          u_1d ( izz, 6, 2) = u ( ixx, iyy, izz, 6, 4)
          u_1d ( izz, 7, 2) = u ( ixx, iyy, izz, 7, 4)
#if ( NFLUID >= 1)
          u_1d ( izz, 8: NFIELD + NFLUID, 2) =
     1                       u ( ixx, iyy, izz, 9: n_cc_components, 4)
#endif

          rhobar ( izz) = rhogr ( ixx, iyy, izz)
          prsbar ( izz) = prsgr ( ixx, iyy, izz)

          vlxbar ( izz) = vlzgr ( ixx, iyy, izz)
          vlybar ( izz) = vlxgr ( ixx, iyy, izz)
          vlzbar ( izz) = vlygr ( ixx, iyy, izz)

          bfxbar ( izz) = bfzgr ( ixx, iyy, izz)
          bfybar ( izz) = bfxgr ( ixx, iyy, izz)
          bfzbar ( izz) = bfygr ( ixx, iyy, izz)

#if ( NFLUID >= 1)
          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( 1: NFLUID, izz) =
     1           u ( ixx, iyy, izz, 9: n_cc_components, 1) * a ( ixx)
#endif

          flatten_1d ( izz) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO izz = loop.

! ----------

          bfxbdy_1d ( 0 - ioffz: iz1 + ioffz) =
     1    bfzbdy ( ixx, iyy, 0 - ioffz: iz1 + ioffz, 1)

          CALL LIMIT_1D_MHD_CONS ( 1, iz1, ioffz,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------

          DO izz = 1 - ioffz, iz1 + ioffz

          u ( ixx, iyy, izz, 1, 4) = u_1d ( izz, 1, 2)
          u ( ixx, iyy, izz, 4, 4) = u_1d ( izz, 2, 2)
          u ( ixx, iyy, izz, 2, 4) = u_1d ( izz, 3, 2)
          u ( ixx, iyy, izz, 3, 4) = u_1d ( izz, 4, 2)
          u ( ixx, iyy, izz, 5, 4) = u_1d ( izz, 5, 2)
          u ( ixx, iyy, izz, 6, 4) = u_1d ( izz, 6, 2)
          u ( ixx, iyy, izz, 7, 4) = u_1d ( izz, 7, 2)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 4) =
     1                       u_1d ( izz, 8: NFIELD + NFLUID, 2)
#endif

          END DO ! End of DO izz = loop.

! ----------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO ixx = loop.

!$OMP END PARALLEL DO

      END IF ! End of "ioffz .GT. 0".

! ----------------------------------------------------------------------
#endif
! End of "LIMIT_TWICE == 1".
! ----------------------------------------------------------------------

! Optionally evaluate and write out the scaled, undivided divergence here.
! In regions with very small magnetic field, it might not be wise to scale
! the undivided divergence by the magnitude of the magnetic field.
! In such situations, either don't scale the divergence at all or
! uncomment the klugey fix-up provided below.

! This is also an opportunity to evaluate and write out the min and
! max of the flow variables if that is desired.

#if ( WRITE_DIAGNOSTICS == 1)

      div_min = 1.0e30
      div_max = - 1.0e30

      flat_min = 1.0e30
      flat_max = - 1.0e30

      rhomin = 1.0e30
      rhomax = - 1.0e30

      vlxmin = 1.0e30
      vlxmax = - 1.0e30
      vlymin = 1.0e30
      vlymax = - 1.0e30
      vlzmin = 1.0e30
      vlzmax = - 1.0e30

      prsmin = 1.0e30
      prsmax = - 1.0e30

      bfxmin = 1.0e30
      bfxmax = - 1.0e30
      bfymin = 1.0e30
      bfymax = - 1.0e30
      bfzmin = 1.0e30
      bfzmax = - 1.0e30

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, dx, dy, dz,                  & 
!$OMP&    ob_dx, ob_dy, ob_dz, a, b, c)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, smallnum, xb_save, dx_save, xc_save,                   &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    & 
!$OMP&    bfxbdy, bfybdy, bfzbdy, msonicgr, flattengr)                  &
!$OMP&           REDUCTION ( MIN: div_min) REDUCTION ( MAX: div_max)    &
!$OMP&           REDUCTION ( MIN: flat_min) REDUCTION ( MAX: flat_max)  &
!$OMP&           REDUCTION ( MIN: rhomin) REDUCTION ( MAX: rhomax)      &
!$OMP&           REDUCTION ( MIN: vlxmin) REDUCTION ( MAX: vlxmax)      &
!$OMP&           REDUCTION ( MIN: vlymin) REDUCTION ( MAX: vlymax)      &
!$OMP&           REDUCTION ( MIN: vlzmin) REDUCTION ( MAX: vlzmax)      &
!$OMP&           REDUCTION ( MIN: prsmin) REDUCTION ( MAX: prsmax)      &
!$OMP&           REDUCTION ( MIN: bfxmin) REDUCTION ( MAX: bfxmax)      &
!$OMP&           REDUCTION ( MIN: bfymin) REDUCTION ( MAX: bfymax)      &
!$OMP&           REDUCTION ( MIN: bfzmin) REDUCTION ( MAX: bfzmax)      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1

! ----------

          IF ( igeom .EQ. 1) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 2) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 3) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1              * dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

! ----------

          a ( ixx) = ( bfxbdy ( ixx, iyy, izz, 1)
     1               - bfxbdy ( ixx - 1, iyy, izz, 1) ) * ob_dx ( ixx)
     1             + ( bfybdy ( ixx, iyy, izz, 1)
     1               - bfybdy ( ixx, iyy - 1, izz, 1) ) * ob_dy ( ixx)
     1             + ( bfzbdy ( ixx, iyy, izz, 1)
     1               - bfzbdy ( ixx, iyy, izz - 1, 1) ) * ob_dz ( ixx)

          a ( ixx) = a ( ixx) * AMIN1 ( dx ( ixx), dy ( ixx),
     1                                             dz ( ixx) )

          b ( ixx) = u ( ixx, iyy, izz, 6, 1)**2
     1             + u ( ixx, iyy, izz, 7, 1)**2
     1             + u ( ixx, iyy, izz, 8, 1)**2
          b ( ixx) = SQRT ( b ( ixx))
          b ( ixx) = AMAX1 ( b ( ixx), smallnum)

          a ( ixx) = a ( ixx) / b ( ixx)
!          a ( ixx) = a ( ixx) / AMAX1 ( b ( ixx),
!     1                     msonicgr ( ixx, iyy, izz)
!     1                     * SQRT ( rhogr ( ixx, iyy, izz) ) )

          div_min = AMIN1 ( div_min, a ( ixx))
          div_max = AMAX1 ( div_max, a ( ixx))

! ----------

          flat_min = AMIN1 ( flat_min, flattengr ( ixx, iyy, izz))
          flat_max = AMAX1 ( flat_max, flattengr ( ixx, iyy, izz))

          rhomin = AMIN1 ( rhomin, rhogr ( ixx, iyy, izz))
          rhomax = AMAX1 ( rhomax, rhogr ( ixx, iyy, izz))

          vlxmin = AMIN1 ( vlxmin, vlxgr ( ixx, iyy, izz))
          vlxmax = AMAX1 ( vlxmax, vlxgr ( ixx, iyy, izz))
          vlymin = AMIN1 ( vlymin, vlygr ( ixx, iyy, izz))
          vlymax = AMAX1 ( vlymax, vlygr ( ixx, iyy, izz))
          vlzmin = AMIN1 ( vlzmin, vlzgr ( ixx, iyy, izz))
          vlzmax = AMAX1 ( vlzmax, vlzgr ( ixx, iyy, izz))

          prsmin = AMIN1 ( prsmin, prsgr ( ixx, iyy, izz))
          prsmax = AMAX1 ( prsmax, prsgr ( ixx, iyy, izz))

          bfxmin = AMIN1 ( bfxmin, bfxgr ( ixx, iyy, izz))
          bfxmax = AMAX1 ( bfxmax, bfxgr ( ixx, iyy, izz))
          bfymin = AMIN1 ( bfymin, bfygr ( ixx, iyy, izz))
          bfymax = AMAX1 ( bfymax, bfygr ( ixx, iyy, izz))
          bfzmin = AMIN1 ( bfzmin, bfzgr ( ixx, iyy, izz))
          bfzmax = AMAX1 ( bfzmax, bfzgr ( ixx, iyy, izz))

! ----------

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      WRITE ( 6, *)"From LIMIT_ALL_SLOPES: ",
     1              "div_min, div_max, flat_min, flat_max = ",
     1              div_min, div_max, flat_min, flat_max
      WRITE ( 6, *)" "
      WRITE ( 6, *)"rhomin, rhomax, prsmin, prsmax = ",
     1              rhomin, rhomax, prsmin, prsmax
      WRITE ( 6, *)" "
      WRITE ( 6, *)"vlxmin, vlxmax, vlymin, vlymax, vlzmin, vlzmax = ",
     1              vlxmin, vlxmax, vlymin, vlymax, vlzmin, vlzmax
      WRITE ( 6, *)" "
      WRITE ( 6, *)"bfxmin, bfxmax, bfymin, bfymax, bfzmin, bfzmax = ",
     1              bfxmin, bfxmax, bfymin, bfymax, bfzmin, bfzmax
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE LIMIT_ALL_SLOPES

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LIMIT_1D_MHD_CONS ( oned_begin, oned_end, oned_off,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Does characteristic limiting for the MHD variables. The variables are
! assumed to come in with the appropriate limiting in physical space
! already done and supplied in "u_1d". The primitive variables are also
! supplied and are used for building eigenvectors. The slopes are
! returned in "u_1d ( :, :, 2)".

! We also check whether the boundary extrapolated states produce positive
! densities and pressures. If not, we zero out the slope in that zone.

! ----------------------------------------------------------------------

      IMPLICIT NONE

! These variables are passed across the subroutine boundary.

      INTEGER oned_begin, oned_end, oned_off

      REAL gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     mc_coef, central_zone_wt

      REAL, DIMENSION ( 0 : NFLUID) ::
     1     cp_specie, molewt_specie, heat_of_formation

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, flatten_1d

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID, 2) :: u_1d

! ----------------------------------------------------------------------

      INTEGER nriembegin, nriemend, iprimitive, nsys, ix, j, k, ifield,
     1        debug_stop

      REAL ob_8pi, r13_12, ta, tb, MinMod, MC_beta

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: xbar, a, b,
     1     u_ip2, u_ip1, u_i, u_im1, u_im2,
     1     R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3,
     1     rh, en, vx, vy, vz, bx, by, bz

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID) ::
     1     eigenval, del_eigwt, del_u

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID) ::
     1     lf_eigenvec, rt_eigenvec

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID, -2 : 2) :: eigwt

! MinMod limiter
      MinMod ( ta, tb) = AMIN1 ( ABS ( ta), ABS ( tb))
     1                    * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( ta, tb, mc_coef) = AMIN1 ( 0.5 * ( ABS ( ta + tb)),
     1                      mc_coef * ABS ( ta), mc_coef * ABS ( tb) )
     1                    * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! ----------------------------------------------------------------------

      COMMON / LIMIT_1D_MHD_CONS_COM / xbar, a, b,

     1     u_ip2, u_ip1, u_i, u_im1, u_im2,
     1     R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3,
     1     rh, en, vx, vy, vz, bx, by, bz,

     1     eigenval, del_eigwt, del_u,

     1     lf_eigenvec, rt_eigenvec,

     1     eigwt,

     1     nriembegin, nriemend, iprimitive, nsys, ix, j, k, ifield,
     1     debug_stop,
     1     ob_8pi, r13_12

!$OMP THREADPRIVATE (/LIMIT_1D_MHD_CONS_COM/)

! ----------------------------------------------------------------------

      nriembegin = oned_begin - oned_off
      nriemend = oned_end + oned_off
      iprimitive = 0

      nsys = NFIELD + NFLUID

      ob_8pi = 1.0 / ( 8.0 * pi)

      r13_12 = 13.0 / 12.0

      debug_stop = 0

! ------------------------------

      xbar ( nriembegin: nriemend) = 0.0

      CALL EIGENVECMHD

     1             ( nriembegin, nriemend, iprimitive,

     1               gamma, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,
     1               cp_specie, molewt_specie, heat_of_formation,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
     1               bfxbar, bfybar, bfzbar, xbar,
     1               speciebar,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )


! ----------------------------------------------------------------------

! Now project into the space of eigenweights.

      eigwt ( oned_begin - oned_off: oned_end + oned_off, :, :) = 0.0

      DO j = - 2, 2

      DO k = 1, nsys

        DO ifield = 1, nsys

          DO ix = oned_begin - oned_off + 2, oned_end + oned_off - 2

          eigwt ( ix, ifield, j) = eigwt ( ix, ifield, j)
     1  + lf_eigenvec ( ix, ifield, k) * u_1d ( ix + j, k, 1)

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

      END DO ! End of "DO k ="

      END DO ! End of "DO j ="

! ----------------------------------------------------------------------
#if ( SLOPE_LIMITER == 1)
! ------------------------------

! MC-beta limiter in characteristic space.

        DO ifield = 1, nsys

          DO ix = oned_begin - oned_off + 2, oned_end + oned_off - 2

          del_eigwt ( ix, ifield) =
     1    MC_beta ( eigwt ( ix, ifield, 1) - eigwt ( ix, ifield, 0),
     1              eigwt ( ix, ifield, 0) - eigwt ( ix, ifield, -1),
     1              mc_coef )

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

! ------------------------------
#else
! ------------------------------

! WENO limiter in characteristic space.

        DO ifield = 1, nsys

          DO ix = oned_begin - oned_off + 2, oned_end + oned_off - 2

          u_ip2 ( ix) = eigwt ( ix, ifield, 2)
          u_ip1 ( ix) = eigwt ( ix, ifield, 1)
          u_i ( ix) = eigwt ( ix, ifield, 0)
          u_im1 ( ix) = eigwt ( ix, ifield, -1)
          u_im2 ( ix) = eigwt ( ix, ifield, -2)

          R1_x ( ix) = - 2.0 * u_im1 ( ix) + 0.5 * u_im2 ( ix)
     1                 + 1.5 * u_i ( ix)
          Is1 ( ix) = r13_12 * ( - 2.0 * u_im1 ( ix)
     1                           + u_im2 ( ix) + u_i ( ix))**2
     1              + 0.25 * ( - 4.0 * u_im1 ( ix) + u_im2 ( ix)
     1                         + 3.0 * u_i ( ix))**2

          R2_x ( ix) = 0.5 * ( u_ip1 ( ix) - u_im1 ( ix) )
          Is2 ( ix) = r13_12 * ( u_im1 ( ix) - 2.0 * u_i ( ix)
     1              + u_ip1 ( ix) )**2
     1              + 0.25 * ( u_im1 ( ix) - u_ip1 ( ix) )**2

          R3_x ( ix) = - 1.5 * u_i ( ix) + 2.0 * u_ip1 ( ix)
     1                 - 0.5 * u_ip2 ( ix)
          Is3 ( ix) = r13_12 * ( u_i ( ix) - 2.0 * u_ip1 ( ix)
     1              + u_ip2 ( ix))**2
     1              + 0.25 * ( 3.0 * u_i ( ix) - 4.0 * u_ip1 ( ix)
     1                       + u_ip2 ( ix))**2

          Wt1 ( ix) = 1.0 / ( Is1 ( ix) + smallnum)**4
          Wt2 ( ix) = central_zone_wt / ( Is2 ( ix) + smallnum)**4
          Wt3 ( ix) = 1.0 / ( Is3 ( ix) + smallnum)**4

          a ( ix) = 1.0 / ( Wt1 ( ix) + Wt2 ( ix) + Wt3 ( ix))
          Wt1 ( ix) = a ( ix) * Wt1 ( ix)
          Wt2 ( ix) = a ( ix) * Wt2 ( ix)
          Wt3 ( ix) = a ( ix) * Wt3 ( ix)

          del_eigwt ( ix, ifield) = Wt1 ( ix) * R1_x ( ix)
     1   + Wt2 ( ix) * R2_x ( ix) + Wt3 ( ix) * R3_x ( ix)

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

! ------------------------------
#endif
! End of "SLOPE_LIMITER ==" part.
! ----------------------------------------------------------------------

! Now project back into physical space.

      del_u ( oned_begin - oned_off: oned_end + oned_off, :) = 0.0

      DO k = 1, nsys

        DO ifield = 1, nsys

          DO ix = oned_begin - oned_off + 2, oned_end + oned_off - 2

          del_u ( ix, k) = del_u ( ix, k) + del_eigwt ( ix, ifield)
     1                               * rt_eigenvec ( ix, k, ifield)

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

      END DO ! End of "DO k ="

! ----------------------------------------------------------------------

! And choose the minimum of the physical space and characteristic space
! limiters. Also introduce a flattener. One could also choose pure
! characteristic limiting.
      
      DO k = 1, nsys

        DO ix = oned_begin - oned_off, oned_end + oned_off

#if ( 1 == 1)
! Here we limit twice.
        u_1d ( ix, k, 2) = MinMod ( u_1d ( ix, k, 2),
     1                     del_u ( ix, k) * ( 1.0 - flatten_1d ( ix)) )
#else
! Here we limit on just the characteristic variables.
        u_1d ( ix, k, 2) = del_u ( ix, k) * ( 1.0 - flatten_1d ( ix))
#endif

        END DO ! End of "DO ix ="

      END DO ! End of "DO k ="

! ----------------------------------------------------------------------

! Check the boundary extrapolated states for positivity of density & pressure.
! If that is not achieved, zero the slopes in that zone.

      DO ix = oned_begin - oned_off, oned_end + oned_off

        rh ( ix) = u_1d ( ix, 1, 1) + 0.5 * u_1d ( ix, 1, 2)
        en ( ix) = u_1d ( ix, 5, 1) + 0.5 * u_1d ( ix, 5, 2)

        vx ( ix) = u_1d ( ix, 2, 1) + 0.5 * u_1d ( ix, 2, 2)
        vy ( ix) = u_1d ( ix, 3, 1) + 0.5 * u_1d ( ix, 3, 2)
        vz ( ix) = u_1d ( ix, 4, 1) + 0.5 * u_1d ( ix, 4, 2)

        bx ( ix) = bfxbdy_1d ( ix)
        by ( ix) = u_1d ( ix, 6, 1) + 0.5 * u_1d ( ix, 6, 2)
        bz ( ix) = u_1d ( ix, 7, 1) + 0.5 * u_1d ( ix, 7, 2)

        IF ( ( rh ( ix) .LE. 0.0) .OR.
     1       ( en ( ix) - 0.5 * ( vx ( ix)**2 + vy ( ix)**2
     1                          + vz ( ix)**2) / rh ( ix)
     1       - ob_8pi * ( bx ( ix)**2 + by ( ix)**2 + bz ( ix)**2)
     1       .LE. 0.0) ) THEN
        u_1d ( ix, :, 2) = 0.0 
        END IF

      END DO

! ----------

      DO ix = oned_begin - oned_off, oned_end + oned_off

        rh ( ix) = u_1d ( ix, 1, 1) - 0.5 * u_1d ( ix, 1, 2)
        en ( ix) = u_1d ( ix, 5, 1) - 0.5 * u_1d ( ix, 5, 2)

        vx ( ix) = u_1d ( ix, 2, 1) - 0.5 * u_1d ( ix, 2, 2)
        vy ( ix) = u_1d ( ix, 3, 1) - 0.5 * u_1d ( ix, 3, 2)
        vz ( ix) = u_1d ( ix, 4, 1) - 0.5 * u_1d ( ix, 4, 2)

        bx ( ix) = bfxbdy_1d ( ix - 1)
        by ( ix) = u_1d ( ix, 6, 1) - 0.5 * u_1d ( ix, 6, 2)
        bz ( ix) = u_1d ( ix, 7, 1) - 0.5 * u_1d ( ix, 7, 2)

        IF ( ( rh ( ix) .LE. 0.0) .OR.
     1       ( en ( ix) - 0.5 * ( vx ( ix)**2 + vy ( ix)**2
     1                          + vz ( ix)**2) / rh ( ix)
     1       - ob_8pi * ( bx ( ix)**2 + by ( ix)**2 + bz ( ix)**2)
     1       .LE. 0.0) ) THEN
        u_1d ( ix, :, 2) = 0.0
        END IF

      END DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE LIMIT_1D_MHD_CONS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_ADER_UPDT ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, 

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy" and produces an ADER-based time
! rate of update. This is stored in "u ( ixx, iyy, izz, :, 5)" within
! each zone. In this version the fluxes are not obtained from the
! Riemann solvers but rather by direct flux evaluation at dynamically
! active zone faces and by diffencing those fluxes.

! In other words,
! we obtain the primitive variabless at points (1/2,0,0) and (-1/2,0,0) and
! use them to construct x-fluxes at those points. Similarly, obtain
! primitive variables at (0,1/2,0) and (0,-1/2,0) and obtain the
! y-fluxes. Similarly, in the z-direction. The finite differences
! of those fluxes then gives us the time rate of update, which is
! stored in "u ( ixx, iyy, izz, :, 5)".

! If pressure positivity is called for, it will also take in "entropygr"
! and build "entropygr ( ixx, iyy, izz, 5)".

! This subroutine is effectively Stage I, i.e. the predictor stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, nspecies_min, nspecies_max,
     1        ioffx_expand, ioffy_expand, ioffz_expand

      REAL ob_8pi, ob_4pi, ob_12

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b, c

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho_pl, rho_mi, prs_pl, prs_mi,
     1     eng_pl, eng_mi, b_prs_pl, b_prs_mi,

     1     vlx_pl, vlx_mi, vly_pl, vly_mi,
     1     vlz_pl, vlz_mi,

     1     bfx_pl, bfx_mi, bfy_pl, bfy_mi, bfz_pl, bfz_mi,

     1     prs_flux_pl, prs_flux_mi, entropy_pl, entropy_mi

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD + NFLUID + 1) ::
     1     flux_pl, flux_mi

! ----------------------------------------------------------------------

      ob_12 = 1.0 / 12.0

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

!      ioffx_expand = expand_by_one_zone * ioffx
!      ioffy_expand = expand_by_one_zone * ioffy
!      ioffz_expand = expand_by_one_zone * ioffz
      ioffx_expand = MIN0 ( 1, ioffx)
      ioffy_expand = MIN0 ( 1, ioffy)
      ioffz_expand = MIN0 ( 1, ioffz)

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_ADER_UPDT : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_ADER_UPDT : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_ADER_UPDT : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

! ----------------------------------------------------------------------
! Treat x-fluxes here.

      IF ( ioffx .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    rho_pl, prs_pl, eng_pl, vlx_pl, vly_pl, vlz_pl, entropy_pl,   & 
!$OMP&    bfx_pl, bfy_pl, bfz_pl, b_prs_pl, flux_pl, prs_flux_pl,       &
!$OMP&    rho_mi, prs_mi, eng_mi, vlx_mi, vly_mi, vlz_mi, entropy_mi,   & 
!$OMP&    bfx_mi, bfy_mi, bfz_mi, b_prs_mi, flux_mi, prs_flux_mi)       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    igeom, smallnum, rhofloor, prsfloor, ob_12, ob_8pi, ob_4pi,   &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    entropygr, gamma_eos_gr)                                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

! ----------------------------------------------------------------------

! Obtain primitive variables; evaluate fluxes and add them in at x-faces.

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! ----------

          rho_pl ( ixx) = u ( ixx, iyy, izz, 1, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 1, 2)

          vlx_pl ( ixx) = u ( ixx, iyy, izz, 2, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 2, 2)

          vly_pl ( ixx) = u ( ixx, iyy, izz, 3, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 3, 2)

          vlz_pl ( ixx) = u ( ixx, iyy, izz, 4, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 4, 2)

          eng_pl ( ixx) = u ( ixx, iyy, izz, 5, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 5, 2)

          bfx_pl ( ixx) = bfxbdy ( ixx, iyy, izz, 1)

          bfy_pl ( ixx) = u ( ixx, iyy, izz, 7, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 7, 2)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          bfz_pl ( ixx) = u ( ixx, iyy, izz, 8, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 8, 2)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

! ----------

          rho_mi ( ixx) = u ( ixx, iyy, izz, 1, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 1, 2)

          vlx_mi ( ixx) = u ( ixx, iyy, izz, 2, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 2, 2)

          vly_mi ( ixx) = u ( ixx, iyy, izz, 3, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 3, 2)

          vlz_mi ( ixx) = u ( ixx, iyy, izz, 4, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 4, 2)

          eng_mi ( ixx) = u ( ixx, iyy, izz, 5, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 5, 2)

          bfx_mi ( ixx) = bfxbdy ( ixx - 1, iyy, izz, 1)

          bfy_mi ( ixx) = u ( ixx, iyy, izz, 7, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 7, 2)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          bfz_mi ( ixx) = u ( ixx, iyy, izz, 8, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 8, 2)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

! ----------

          rho_pl ( ixx) = AMAX1 ( rhofloor, rho_pl ( ixx))
          a ( ixx) = 1.0 / rho_pl ( ixx)

          vlx_pl ( ixx) = vlx_pl ( ixx) * a ( ixx)
          vly_pl ( ixx) = vly_pl ( ixx) * a ( ixx)
          vlz_pl ( ixx) = vlz_pl ( ixx) * a ( ixx)

          b ( ixx) = 0.5 * rho_pl ( ixx) * ( vlx_pl ( ixx)**2
     1      + vly_pl ( ixx)**2 + vlz_pl ( ixx)**2)

          b_prs_pl ( ixx) = ( bfx_pl ( ixx)**2 + bfy_pl ( ixx)**2
     1                      + bfz_pl ( ixx)**2) * ob_8pi

          prs_pl ( ixx) = ( eng_pl ( ixx) - b ( ixx) - b_prs_pl ( ixx))
     1                  * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_pl ( ixx) = entropygr ( ixx, iyy, izz, 1)
     1             + 0.5 * entropygr ( ixx, iyy, izz, 2)
          IF ( prs_pl ( ixx) .LE. prsfloor) prs_pl ( ixx) =
     1    entropy_pl ( ixx)
     1  * rho_pl ( ixx)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#endif
          prs_pl ( ixx) = AMAX1 ( prsfloor, prs_pl ( ixx))

! ----------

          rho_mi ( ixx) = AMAX1 ( rhofloor, rho_mi ( ixx))
          a ( ixx) = 1.0 / rho_mi ( ixx)

          vlx_mi ( ixx) = vlx_mi ( ixx) * a ( ixx)
          vly_mi ( ixx) = vly_mi ( ixx) * a ( ixx)
          vlz_mi ( ixx) = vlz_mi ( ixx) * a ( ixx)

          b ( ixx) = 0.5 * rho_mi ( ixx) * ( vlx_mi ( ixx)**2
     1      + vly_mi ( ixx)**2 + vlz_mi ( ixx)**2)

          b_prs_mi ( ixx) = ( bfx_mi ( ixx)**2 + bfy_mi ( ixx)**2
     1                      + bfz_mi ( ixx)**2) * ob_8pi

          prs_mi ( ixx) = ( eng_mi ( ixx) - b ( ixx) - b_prs_mi ( ixx))
     1                  * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_mi ( ixx) = entropygr ( ixx, iyy, izz, 1)
     1             - 0.5 * entropygr ( ixx, iyy, izz, 2)
          IF ( prs_mi ( ixx) .LE. prsfloor) prs_mi ( ixx) =
     1    entropy_mi ( ixx)
     1  * rho_mi ( ixx)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#endif
          prs_mi ( ixx) = AMAX1 ( prsfloor, prs_mi ( ixx))

! ----------

          flux_pl ( ixx, 1) = rho_pl ( ixx) * vlx_pl ( ixx)

          flux_pl ( ixx, 2) = rho_pl ( ixx) * vlx_pl ( ixx)**2
     1                      - bfx_pl ( ixx)**2 * ob_4pi

          prs_flux_pl ( ixx) = prs_pl ( ixx) + b_prs_pl ( ixx)

          flux_pl ( ixx, 3) = rho_pl ( ixx) * vlx_pl ( ixx)
     1                                      * vly_pl ( ixx)
     1  - bfx_pl ( ixx) * bfy_pl ( ixx) * ob_4pi

          flux_pl ( ixx, 4) = rho_pl ( ixx) * vlx_pl ( ixx)
     1                                      * vlz_pl ( ixx)
     1  - bfx_pl ( ixx) * bfz_pl ( ixx) * ob_4pi

          flux_pl ( ixx, 5) = ( eng_pl ( ixx) + prs_pl ( ixx)
     1  + b_prs_pl ( ixx) ) * vlx_pl ( ixx)
     1  - bfx_pl ( ixx) * ob_4pi * ( vlx_pl ( ixx) * bfx_pl ( ixx)
     1  + vly_pl ( ixx) * bfy_pl ( ixx) + vlz_pl ( ixx) * bfz_pl ( ixx))

          flux_pl ( ixx, 6) = 0.0

          flux_pl ( ixx, 7) = vlx_pl ( ixx) * bfy_pl ( ixx)
     1                      - vly_pl ( ixx) * bfx_pl ( ixx)

          flux_pl ( ixx, 8) = vlx_pl ( ixx) * bfz_pl ( ixx)
     1                      - vlz_pl ( ixx) * bfx_pl ( ixx)

#if ( NFLUID >= 1)
          flux_pl ( ixx, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 2) )
     1          * vlx_pl ( ixx)
#endif

! ----------

          flux_mi ( ixx, 1) = rho_mi ( ixx) * vlx_mi ( ixx)

          flux_mi ( ixx, 2) = rho_mi ( ixx) * vlx_mi ( ixx)**2
     1                      - bfx_mi ( ixx)**2 * ob_4pi

          prs_flux_mi ( ixx) = prs_mi ( ixx) + b_prs_mi ( ixx)

          flux_mi ( ixx, 3) = rho_mi ( ixx) * vlx_mi ( ixx)
     1                                      * vly_mi ( ixx)
     1  - bfx_mi ( ixx) * bfy_mi ( ixx) * ob_4pi

          flux_mi ( ixx, 4) = rho_mi ( ixx) * vlx_mi ( ixx)
     1                                      * vlz_mi ( ixx)
     1  - bfx_mi ( ixx) * bfz_mi ( ixx) * ob_4pi

          flux_mi ( ixx, 5) = ( eng_mi ( ixx) + prs_mi ( ixx) 
     1  + b_prs_mi ( ixx) ) * vlx_mi ( ixx)
     1  - bfx_mi ( ixx) * ob_4pi * ( vlx_mi ( ixx) * bfx_mi ( ixx)
     1  + vly_mi ( ixx) * bfy_mi ( ixx) + vlz_mi ( ixx) * bfz_mi ( ixx))

          flux_mi ( ixx, 6) = 0.0

          flux_mi ( ixx, 7) = vlx_mi ( ixx) * bfy_mi ( ixx)
     1                      - vly_mi ( ixx) * bfx_mi ( ixx)

          flux_mi ( ixx, 8) = vlx_mi ( ixx) * bfz_mi ( ixx)
     1                      - vlz_mi ( ixx) * bfx_mi ( ixx)

#if ( NFLUID >= 1)
          flux_mi ( ixx, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      - 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 2) )
     1          * vlx_mi ( ixx)
#endif

! ----------

          END DO ! End of DO ixx = loop.

! ----------
! Cartesian

          IF ( igeom .EQ. 1) THEN

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( ixx, 1: nspecies_max)
     1    - flux_mi ( ixx, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( ixx) * vlx_pl ( ixx)
     1    - entropy_mi ( ixx) * vlx_mi ( ixx) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_pl ( ixx) - prs_flux_mi ( ixx) ) * a ( ixx)

          END DO ! End of DO ixx = loop.

          END IF

! -----
! Cylindrical

          IF ( igeom .EQ. 2) THEN

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

          a ( ixx) = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = xb_save ( ixx)

          c ( ixx) = xb_save ( ixx - 1)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( ixx, 1: nspecies_max) * b ( ixx)
     1    - flux_mi ( ixx, 1: nspecies_max) * c ( ixx) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( ixx) * vlx_pl ( ixx) * b ( ixx)
     1    - entropy_mi ( ixx) * vlx_mi ( ixx) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_pl ( ixx) - prs_flux_mi ( ixx) ) * a ( ixx)

          END DO ! End of DO ixx = loop.

          END IF

! -----
! Spherical

          IF ( igeom .EQ. 3) THEN

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

          a ( ixx) = ( xb_save ( ixx)**3
     1               - xb_save ( ixx - 1)**3 ) / 3.0
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = xb_save ( ixx)**2

          c ( ixx) = xb_save ( ixx - 1)**2

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( ixx, 1: nspecies_max) * b ( ixx)
     1    - flux_mi ( ixx, 1: nspecies_max) * c ( ixx) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( ixx) * vlx_pl ( ixx) * b ( ixx)
     1    - entropy_mi ( ixx) * vlx_mi ( ixx) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_pl ( ixx) - prs_flux_mi ( ixx) ) * a ( ixx)

          END DO ! End of DO ixx = loop.

          END IF

! ----------------------------------------------------------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "ioffx .GT. 0".

! End of treating x-fluxes.
! ----------------------------------------------------------------------
! Treat y-fluxes here.

      IF ( ioffy .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    rho_pl, prs_pl, eng_pl, vlx_pl, vly_pl, vlz_pl, entropy_pl,   &
!$OMP&    bfx_pl, bfy_pl, bfz_pl, b_prs_pl, flux_pl, prs_flux_pl,       &
!$OMP&    rho_mi, prs_mi, eng_mi, vlx_mi, vly_mi, vlz_mi, entropy_mi,   &
!$OMP&    bfx_mi, bfy_mi, bfz_mi, b_prs_mi, flux_mi, prs_flux_mi)       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    igeom, smallnum, rhofloor, prsfloor, ob_12, ob_8pi, ob_4pi,   &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    entropygr, gamma_eos_gr)                                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! ----------------------------------------------------------------------

! Obtain primitive variables; evaluate fluxes and add them in at y-faces.

          DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

! ----------

          rho_pl ( iyy) = u ( ixx, iyy, izz, 1, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 1, 3)

          vlx_pl ( iyy) = u ( ixx, iyy, izz, 2, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 2, 3)

          vly_pl ( iyy) = u ( ixx, iyy, izz, 3, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 3, 3)

          vlz_pl ( iyy) = u ( ixx, iyy, izz, 4, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 4, 3)

          eng_pl ( iyy) = u ( ixx, iyy, izz, 5, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 5, 3)

          bfx_pl ( iyy) = u ( ixx, iyy, izz, 6, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 6, 3)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          bfy_pl ( iyy) = bfybdy ( ixx, iyy, izz, 1)

          bfz_pl ( iyy) = u ( ixx, iyy, izz, 8, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 8, 3)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

! ----------

          rho_mi ( iyy) = u ( ixx, iyy, izz, 1, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 1, 3)

          vlx_mi ( iyy) = u ( ixx, iyy, izz, 2, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 2, 3)

          vly_mi ( iyy) = u ( ixx, iyy, izz, 3, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 3, 3)

          vlz_mi ( iyy) = u ( ixx, iyy, izz, 4, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 4, 3)

          eng_mi ( iyy) = u ( ixx, iyy, izz, 5, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 5, 3)

          bfx_mi ( iyy) = u ( ixx, iyy, izz, 6, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 6, 3)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          bfy_mi ( iyy) = bfybdy ( ixx, iyy - 1, izz, 1)

          bfz_mi ( iyy) = u ( ixx, iyy, izz, 8, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 8, 3)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

! ----------

          rho_pl ( iyy) = AMAX1 ( rhofloor, rho_pl ( iyy))
          a ( iyy) = 1.0 / rho_pl ( iyy)

          vlx_pl ( iyy) = vlx_pl ( iyy) * a ( iyy)
          vly_pl ( iyy) = vly_pl ( iyy) * a ( iyy)
          vlz_pl ( iyy) = vlz_pl ( iyy) * a ( iyy)

          b ( iyy) = 0.5 * rho_pl ( iyy) * ( vlx_pl ( iyy)**2
     1      + vly_pl ( iyy)**2 + vlz_pl ( iyy)**2)

          b_prs_pl ( iyy) = ( bfx_pl ( iyy)**2 + bfy_pl ( iyy)**2
     1                      + bfz_pl ( iyy)**2) * ob_8pi

          prs_pl ( iyy) = ( eng_pl ( iyy) - b ( iyy) - b_prs_pl ( iyy))
     1                  * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_pl ( iyy) = entropygr ( ixx, iyy, izz, 1)
     1             + 0.5 * entropygr ( ixx, iyy, izz, 3)
          IF ( prs_pl ( iyy) .LE. prsfloor) prs_pl ( iyy) =
     1    entropy_pl ( iyy)
     1  * rho_pl ( iyy)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#endif
          prs_pl ( iyy) = AMAX1 ( prsfloor, prs_pl ( iyy))

! ----------

          rho_mi ( iyy) = AMAX1 ( rhofloor, rho_mi ( iyy))
          a ( iyy) = 1.0 / rho_mi ( iyy)

          vlx_mi ( iyy) = vlx_mi ( iyy) * a ( iyy)
          vly_mi ( iyy) = vly_mi ( iyy) * a ( iyy)
          vlz_mi ( iyy) = vlz_mi ( iyy) * a ( iyy)

          b ( iyy) = 0.5 * rho_mi ( iyy) * ( vlx_mi ( iyy)**2
     1      + vly_mi ( iyy)**2 + vlz_mi ( iyy)**2)

          b_prs_mi ( iyy) = ( bfx_mi ( iyy)**2 + bfy_mi ( iyy)**2
     1                      + bfz_mi ( iyy)**2) * ob_8pi

          prs_mi ( iyy) = ( eng_mi ( iyy) - b ( iyy) - b_prs_mi ( iyy))
     1                  * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_mi ( iyy) = entropygr ( ixx, iyy, izz, 1)
     1             - 0.5 * entropygr ( ixx, iyy, izz, 3)
          IF ( prs_mi ( iyy) .LE. prsfloor) prs_mi ( iyy) =
     1    entropy_mi ( iyy)
     1  * rho_mi ( iyy)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#endif
          prs_mi ( iyy) = AMAX1 ( prsfloor, prs_mi ( iyy))

! ----------

          flux_pl ( iyy, 1) = rho_pl ( iyy) * vly_pl ( iyy)

          flux_pl ( iyy, 2) = rho_pl ( iyy) * vly_pl ( iyy)
     1                      * vlx_pl ( iyy)
     1  - bfy_pl ( iyy) * bfx_pl ( iyy) * ob_4pi

          flux_pl ( iyy, 3) = rho_pl ( iyy) * vly_pl ( iyy)**2
     1                      - bfy_pl ( iyy)**2 * ob_4pi

          prs_flux_pl ( iyy) = prs_pl ( iyy) + b_prs_pl ( iyy)

          flux_pl ( iyy, 4) = rho_pl ( iyy) * vly_pl ( iyy)
     1                      * vlz_pl ( iyy)
     1  - bfy_pl ( iyy) * bfz_pl ( iyy) * ob_4pi

          flux_pl ( iyy, 5) = ( eng_pl ( iyy) + prs_pl ( iyy)
     1  + b_prs_pl ( iyy) ) * vly_pl ( iyy)
     1  - bfy_pl ( iyy) * ob_4pi * ( vlx_pl ( iyy) * bfx_pl ( iyy)
     1  + vly_pl ( iyy) * bfy_pl ( iyy) + vlz_pl ( iyy) * bfz_pl ( iyy))

          flux_pl ( iyy, 6) = vly_pl ( iyy) * bfx_pl ( iyy)
     1                      - vlx_pl ( iyy) * bfy_pl ( iyy)

          flux_pl ( iyy, 7) = 0.0

          flux_pl ( iyy, 8) = vly_pl ( iyy) * bfz_pl ( iyy)
     1                      - vlz_pl ( iyy) * bfy_pl ( iyy)

#if ( NFLUID >= 1)
          flux_pl ( iyy, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 3) )
     1          * vly_pl ( iyy)
#endif

! ----------

          flux_mi ( iyy, 1) = rho_mi ( iyy) * vly_mi ( iyy)

          flux_mi ( iyy, 2) = rho_mi ( iyy) * vly_mi ( iyy)
     1                      * vlx_mi ( iyy)
     1  - bfy_mi ( iyy) * bfx_mi ( iyy) * ob_4pi

          flux_mi ( iyy, 3) = rho_mi ( iyy) * vly_mi ( iyy)**2
     1                      - bfy_mi ( iyy)**2 * ob_4pi

          prs_flux_mi ( iyy) = prs_mi ( iyy) + b_prs_mi ( iyy)

          flux_mi ( iyy, 4) = rho_mi ( iyy) * vly_mi ( iyy)
     1                      * vlz_mi ( iyy)
     1  - bfy_mi ( iyy) * bfz_mi ( iyy) * ob_4pi

          flux_mi ( iyy, 5) = ( eng_mi ( iyy) + prs_mi ( iyy)
     1  + b_prs_mi ( iyy) ) * vly_mi ( iyy)
     1  - bfy_mi ( iyy) * ob_4pi * ( vlx_mi ( iyy) * bfx_mi ( iyy)
     1  + vly_mi ( iyy) * bfy_mi ( iyy) + vlz_mi ( iyy) * bfz_mi ( iyy))

          flux_mi ( iyy, 6) = vly_mi ( iyy) * bfx_mi ( iyy)
     1                      - vlx_mi ( iyy) * bfy_mi ( iyy)

          flux_mi ( iyy, 7) = 0.0

          flux_mi ( iyy, 8) = vly_mi ( iyy) * bfz_mi ( iyy)
     1                      - vlz_mi ( iyy) * bfy_mi ( iyy)

#if ( NFLUID >= 1)
          flux_mi ( iyy, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      - 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 3) )
     1          * vly_mi ( iyy)
#endif

! ----------

          END DO ! End of DO iyy = loop.

! ----------
! Cartesian

          IF ( igeom .EQ. 1) THEN

          DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          a ( iyy) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( iyy, 1: nspecies_max)
     1    - flux_mi ( iyy, 1: nspecies_max) ) * a ( iyy)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( iyy) * vly_pl ( iyy)
     1    - entropy_mi ( iyy) * vly_mi ( iyy) ) * a ( iyy)
#endif

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_pl ( iyy) - prs_flux_mi ( iyy) ) * a ( iyy)

          END DO ! End of DO iyy = loop.

          END IF

! -----
! Cylindrical

          IF ( igeom .EQ. 2) THEN

          DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          a ( iyy) = xc_save ( ixx) * dy_save ( iyy)
          a ( iyy) = dtcur / a ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( iyy, 1: nspecies_max)
     1    - flux_mi ( iyy, 1: nspecies_max) ) * a ( iyy)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( iyy) * vly_pl ( iyy)
     1    - entropy_mi ( iyy) * vly_mi ( iyy) ) * a ( iyy)
#endif

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_pl ( iyy) - prs_flux_mi ( iyy) ) * a ( iyy)

          END DO ! End of DO iyy = loop.

          END IF

! -----
! Spherical

          IF ( igeom .EQ. 3) THEN

          DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          a ( iyy) = - xc_save ( ixx) * ( COS ( yb_save ( iyy))
     1                           - COS ( yb_save ( iyy - 1)) )
          a ( iyy) = dtcur / a ( iyy)

          b ( iyy) = SIN ( yb_save ( iyy))

          c ( iyy) = SIN ( yb_save ( iyy - 1))

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( iyy, 1: nspecies_max) * b ( iyy)
     1    - flux_mi ( iyy, 1: nspecies_max) * c ( iyy) ) * a ( iyy)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( iyy) * vly_pl ( iyy) * b ( iyy)
     1    - entropy_mi ( iyy) * vly_mi ( iyy) * c ( iyy) ) * a ( iyy)
#endif

          a ( iyy) = xc_save ( ixx) * dy_save ( iyy)
          a ( iyy) = dtcur / a ( iyy)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_pl ( iyy) - prs_flux_mi ( iyy) ) * a ( iyy)

          END DO ! End of DO iyy = loop.

          END IF

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "ioffy .GT. 0".

! End of treating y-fluxes.
! ----------------------------------------------------------------------
! Treat z-fluxes here.

      IF ( ioffz .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    rho_pl, prs_pl, eng_pl, vlx_pl, vly_pl, vlz_pl, entropy_pl,   &
!$OMP&    bfx_pl, bfy_pl, bfz_pl, b_prs_pl, flux_pl, prs_flux_pl,       &
!$OMP&    rho_mi, prs_mi, eng_mi, vlx_mi, vly_mi, vlz_mi, entropy_mi,   &
!$OMP&    bfx_mi, bfy_mi, bfz_mi, b_prs_mi, flux_mi, prs_flux_mi)       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    igeom, smallnum, rhofloor, prsfloor, ob_12, ob_8pi, ob_4pi,   &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    entropygr, gamma_eos_gr)                                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

        DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! ----------------------------------------------------------------------

! Obtain primitive variables; evaluate fluxes and add them in at z-faces.

          DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

! ----------

          rho_pl ( izz) = u ( ixx, iyy, izz, 1, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 1, 4)

          vlx_pl ( izz) = u ( ixx, iyy, izz, 2, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 2, 4)

          vly_pl ( izz) = u ( ixx, iyy, izz, 3, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 3, 4)

          vlz_pl ( izz) = u ( ixx, iyy, izz, 4, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 4, 4)

          eng_pl ( izz) = u ( ixx, iyy, izz, 5, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 5, 4)

          bfx_pl ( izz) = u ( ixx, iyy, izz, 6, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 6, 4)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          bfy_pl ( izz) = u ( ixx, iyy, izz, 7, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 7, 4)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          bfz_pl ( izz) = bfzbdy ( ixx, iyy, izz, 1)

! ----------

          rho_mi ( izz) = u ( ixx, iyy, izz, 1, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 1, 4)

          vlx_mi ( izz) = u ( ixx, iyy, izz, 2, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 2, 4)

          vly_mi ( izz) = u ( ixx, iyy, izz, 3, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 3, 4)

          vlz_mi ( izz) = u ( ixx, iyy, izz, 4, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 4, 4)

          eng_mi ( izz) = u ( ixx, iyy, izz, 5, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 5, 4)

          bfx_mi ( izz) = u ( ixx, iyy, izz, 6, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 6, 4)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          bfy_mi ( izz) = u ( ixx, iyy, izz, 7, 1)
     1             - 0.5 * u ( ixx, iyy, izz, 7, 4)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1       - ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          bfz_mi ( izz) = bfzbdy ( ixx, iyy, izz - 1, 1)

! ----------

          rho_pl ( izz) = AMAX1 ( rhofloor, rho_pl ( izz))
          a ( izz) = 1.0 / rho_pl ( izz)

          vlx_pl ( izz) = vlx_pl ( izz) * a ( izz)
          vly_pl ( izz) = vly_pl ( izz) * a ( izz)
          vlz_pl ( izz) = vlz_pl ( izz) * a ( izz)

          b ( izz) = 0.5 * rho_pl ( izz) * ( vlx_pl ( izz)**2
     1      + vly_pl ( izz)**2 + vlz_pl ( izz)**2)

          b_prs_pl ( izz) = ( bfx_pl ( izz)**2 + bfy_pl ( izz)**2
     1                      + bfz_pl ( izz)**2) * ob_8pi

          prs_pl ( izz) = ( eng_pl ( izz) - b ( izz) - b_prs_pl ( izz))
     1                  * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_pl ( izz) = entropygr ( ixx, iyy, izz, 1)
     1             + 0.5 * entropygr ( ixx, iyy, izz, 4)
          IF ( prs_pl ( izz) .LE. prsfloor) prs_pl ( izz) =
     1    entropy_pl ( izz)
     1  * rho_pl ( izz)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#endif
          prs_pl ( izz) = AMAX1 ( prsfloor, prs_pl ( izz))

! ----------

          rho_mi ( izz) = AMAX1 ( rhofloor, rho_mi ( izz))
          a ( izz) = 1.0 / rho_mi ( izz)

          vlx_mi ( izz) = vlx_mi ( izz) * a ( izz)
          vly_mi ( izz) = vly_mi ( izz) * a ( izz)
          vlz_mi ( izz) = vlz_mi ( izz) * a ( izz)

          b ( izz) = 0.5 * rho_mi ( izz) * ( vlx_mi ( izz)**2
     1      + vly_mi ( izz)**2 + vlz_mi ( izz)**2)

          b_prs_mi ( izz) = ( bfx_mi ( izz)**2 + bfy_mi ( izz)**2
     1                      + bfz_mi ( izz)**2) * ob_8pi

          prs_mi ( izz) = ( eng_mi ( izz) - b ( izz) - b_prs_mi ( izz))
     1                  * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_mi ( izz) = entropygr ( ixx, iyy, izz, 1)
     1             - 0.5 * entropygr ( ixx, iyy, izz, 4)
          IF ( prs_mi ( izz) .LE. prsfloor) prs_mi ( izz) =
     1    entropy_mi ( izz)
     1  * rho_mi ( izz)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
#endif
          prs_mi ( izz) = AMAX1 ( prsfloor, prs_mi ( izz))

! ----------

          flux_pl ( izz, 1) = rho_pl ( izz) * vlz_pl ( izz)

          flux_pl ( izz, 2) = rho_pl ( izz) * vlz_pl ( izz)
     1                      * vlx_pl ( izz)
     1  - bfz_pl ( izz) * bfx_pl ( izz) * ob_4pi

          flux_pl ( izz, 3) = rho_pl ( izz) * vlz_pl ( izz)
     1                      * vly_pl ( izz)
     1  - bfz_pl ( izz) * bfy_pl ( izz) * ob_4pi

          flux_pl ( izz, 4) = rho_pl ( izz) * vlz_pl ( izz)**2
     1                      - bfz_pl ( izz)**2 * ob_4pi

          prs_flux_pl ( izz) = prs_pl ( izz) + b_prs_pl ( izz)

          flux_pl ( izz, 5) = ( eng_pl ( izz) + prs_pl ( izz)
     1  + b_prs_pl ( izz) ) * vlz_pl ( izz)
     1  - bfz_pl ( izz) * ob_4pi * ( vlx_pl ( izz) * bfx_pl ( izz)
     1  + vly_pl ( izz) * bfy_pl ( izz) + vlz_pl ( izz) * bfz_pl ( izz))

          flux_pl ( izz, 6) = vlz_pl ( izz) * bfx_pl ( izz)
     1                      - vlx_pl ( izz) * bfz_pl ( izz)

          flux_pl ( izz, 7) = vlz_pl ( izz) * bfy_pl ( izz)
     1                      - vly_pl ( izz) * bfz_pl ( izz)

          flux_pl ( izz, 8) = 0.0

#if ( NFLUID >= 1)
          flux_pl ( izz, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 4) )
     1          * vlz_pl ( izz)
#endif

! ----------

          flux_mi ( izz, 1) = rho_mi ( izz) * vlz_mi ( izz)

          flux_mi ( izz, 2) = rho_mi ( izz) * vlz_mi ( izz)
     1                      * vlx_mi ( izz)
     1  - bfz_mi ( izz) * bfx_mi ( izz) * ob_4pi

          flux_mi ( izz, 3) = rho_mi ( izz) * vlz_mi ( izz)
     1                      * vly_mi ( izz)
     1  - bfz_mi ( izz) * bfy_mi ( izz) * ob_4pi

          flux_mi ( izz, 4) = rho_mi ( izz) * vlz_mi ( izz)**2
     1                      - bfz_mi ( izz)**2 * ob_4pi

          prs_flux_mi ( izz) = prs_mi ( izz) + b_prs_mi ( izz)

          flux_mi ( izz, 5) = ( eng_mi ( izz) + prs_mi ( izz)
     1  + b_prs_mi ( izz) ) * vlz_mi ( izz)
     1  - bfz_mi ( izz) * ob_4pi * ( vlx_mi ( izz) * bfx_mi ( izz)
     1  + vly_mi ( izz) * bfy_mi ( izz) + vlz_mi ( izz) * bfz_mi ( izz))

          flux_mi ( izz, 6) = vlz_mi ( izz) * bfx_mi ( izz)
     1                      - vlx_mi ( izz) * bfz_mi ( izz)

          flux_mi ( izz, 7) = vlz_mi ( izz) * bfy_mi ( izz)
     1                      - vly_mi ( izz) * bfz_mi ( izz)

          flux_mi ( izz, 8) = 0.0

#if ( NFLUID >= 1)
          flux_mi ( izz, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      - 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 4) )
     1          * vlz_mi ( izz)
#endif

! ----------

          END DO ! End of DO izz = loop.

! ----------
! Cartesian & Cylindrical

          IF ( ( igeom .EQ. 1) .OR. ( igeom .EQ. 2)) THEN

          DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

          a ( izz) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( izz, 1: nspecies_max)
     1    - flux_mi ( izz, 1: nspecies_max) ) * a ( izz)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( izz) * vlz_pl ( izz)
     1    - entropy_mi ( izz) * vlz_mi ( izz) ) * a ( izz)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_pl ( izz) - prs_flux_mi ( izz) ) * a ( izz)

          END DO ! End of DO izz = loop.

          END IF

! -----
! Spherical

      IF ( igeom .EQ. 3) THEN

          DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

          a ( izz) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1             * dz_save ( izz)
          a ( izz) = dtcur / a ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_pl ( izz, 1: nspecies_max)
     1    - flux_mi ( izz, 1: nspecies_max) ) * a ( izz)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_pl ( izz) * vlz_pl ( izz)
     1    - entropy_mi ( izz) * vlz_mi ( izz) ) * a ( izz)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_pl ( izz) - prs_flux_mi ( izz) ) * a ( izz)

          END DO ! End of DO izz = loop.

          END IF

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

      END IF ! End of "ioffz .GT. 0".

! End of treating z-fluxes.
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_ADER_UPDT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE APPLY_MULTID_RS_X_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_y, flux_z, flux_temp, efxgr,
     1        entropy_flux_y, entropy_flux_z, entropy_flux_temp,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and evaluates the MuSIC Riemann solver at each x-edge. One resolved
! state and two fluxes are produced, so we add in the contribution from the
! x-edge to "flux_y, flux_z". Also builds "efxgr" at that edge. All fluxes must
! be zeroed before calls to this suite of three subroutines.

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

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
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1                  bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_y, flux_z, flux_temp

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) :: efxgr,
     1        entropy_flux_y, entropy_flux_z, entropy_flux_temp

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr

! ----------

      INTEGER, DIMENSION ( 2,
     1                     1 - ioffy : iy1 + ioffy,
     1                     1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                     2,
     1                     1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                     1 - ioffy : iy1 + ioffy,
     1                     2) ::
     1         z_indx_limits, bcarr_z

! ----------

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, iyym1, iyy0, iyyp1, izzm1, izz0, izzp1,
     1        nriembegin, nriemend, idirection,
     1        ioffx_expand, ioffy_expand, ioffz_expand,
     1        multid_wavemodel_num, split_prsflux_re

      REAL ob_6_multid, ob_6_oned, ob_6, ob_12

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + 1 + NFLUID) ::

     1      U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,

     1      gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,
     1      gamma_eos_Rt_Dn,
     1      gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,
     1      gamma_sdsp_Rt_Dn,

     1      blend_HLL,

     1      Prs_Star,

     1      temp_01, temp_02, temp_03, temp_04, temp_05, dy_use, dz_use

! ----------------------------------------------------------------------

! These are the arrays that actually go through the vertex-based RS.

      INTEGER nzone_vrtx

      REAL

! Amount of HLL we will have to blend in for stabilization.
     1     my_blend_HLL, my_blend_MDHLL,

! Velocity of vertex; total and gas pressure in resolved state.
     1     Vlx_ALE, Vly_ALE,
     1     my_Prs_Tot_Star, my_Prs_Gas_Star

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables that come together at a zone edge.
! We also come in with their primitive variables and associated
! x- and y-fluxes.
     1      U_vrtx, F_vrtx, G_vrtx,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

! entropy variables that go through RS interface.
     1      entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,

! Components of the unit normal.
     1      nx, ny

! Final, strongly-interacting state and Global x- and y-fluxes.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) ::
     1      U_Star, F_Star, G_Star

! ----------------------------------------------------------------------

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

      ob_6 = 1.0 / 6.0
      ob_12 = 1.0 / 12.0

! ----------

! We are interested in the fluxes for a non-moving mesh.
      Vlx_ALE = 0.0
      Vly_ALE = 0.0

! ----------

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

! In true 2d, only the x-edges contribute to a problem whose variation
! is entirely in the yz-plane. If this is truly a 2d problem with no variation
! in the x-direction then only two edges contribute fluxes to a face
! (instead of four), so "ob_6_multid and ob_6_oned" are increased.

        multid_wavemodel_num = 1  ! 1 includes supersonic; 2 is fully subsonic.

! "4 * ob_6_multid + 4 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux
        ob_6_oned = 0.25 * ( 1.0 - 4.0 * blend_multid_flux)

#if ( TRUE_2D == 1)
        IF ( ioffx .NE. 0) RETURN
! "2 * ob_6_multid + 2 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux
        ob_6_oned = 0.5 * ( 1.0 - 2.0 * blend_multid_flux)
#endif

! ----------

      idirection = 1

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in APPLY_MULTID_RS_X_EDGE : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iyym1, iyy0, iyyp1,          &
!$OMP&    izzm1, izz0, izzp1, nriembegin, nriemend,                     &
!$OMP&    U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn,                           &
!$OMP&    entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,   &

!$OMP&    gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,            &
!$OMP&    gamma_eos_Rt_Dn,                                              &
!$OMP&    gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,         &
!$OMP&    gamma_sdsp_Rt_Dn,                                             &

!$OMP&    blend_HLL, my_Prs_Gas_Star,                                   &
!$OMP&    my_blend_HLL, my_blend_MDHLL, my_Prs_Tot_Star,                &
!$OMP&    U_vrtx, F_vrtx, G_vrtx,                                       &
!$OMP&    gamma_eos, gamma_sdsp, U_Star, F_Star, G_Star,                &
!$OMP&    U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,     &
!$OMP&    entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,             &

!$OMP&    temp_01, temp_02, temp_03, temp_04, temp_05)                  &

!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    ob_6_multid, ob_6_oned, blend_mdhll_to_mdhllc,                &
!$OMP&    dt_shift, smallnum, prsfloor, rhofloor, pi, idirection,       &

!$OMP&    nzone_vrtx, Vlx_ALE, Vly_ALE, nx, ny, ob_6, ob_12,            &
!$OMP&    multid_wavemodel_num, split_prsflux_re,                       &

!$OMP&    x_indx_limits, y_indx_limits, z_indx_limits,                  &
!$OMP&    bcarr_x, bcarr_y, bcarr_z, dx_save, dy_save, dz_save,         &
!$OMP&    u, entropygr, gamma_eos_gr, gamma_soundspeed_gr, flattengr,   &
!$OMP&    bfxccmodes, bfyccmodes, bfzccmodes,                           &
!$OMP&    efxgr, flux_y, flux_z, flux_temp,                             &
!$OMP&    entropy_flux_y, entropy_flux_z, entropy_flux_temp)            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 0 - ioffy_expand, iy1 + ioffy_expand

! ----------------------------------------------------------------------
! Indices ( iyy0, izz0) always map to an interior point on the mesh including
! its ghost zones. This is needed because ( iyy, izz) map to the edges.

        izzm1 = MAX0 ( izz - 1, 0 - ioffz)
        izz0 = MAX0 ( izz, 1 - ioffz)
        izzp1 = MIN0 ( izz + 1, iz1 + ioffz)

        iyym1 = MAX0 ( iyy - 1, 0 - ioffy)
        iyy0 = MAX0 ( iyy, 1 - ioffy)
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)

        nriembegin = x_indx_limits ( 1, iyy0, izz0) - ioffx_expand
        nriemend = x_indx_limits ( 2, iyy0, izz0) + ioffx_expand

! ----------------------------------------------------------------------
! Load in the four states that surround the x-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------

          DO ixx = nriembegin, nriemend

! ------------------------------
! First, we build the four states.

! ----------
! From zone Lf_Dn == ( ixx, iyy0, izz0) go to y = 0.5, z = 0.5

          U_Lf_Dn ( ixx, :) = u ( ixx, iyy0, izz0, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 4)
     1           + dt_shift * u ( ixx, iyy0, izz0, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Lf_Dn ( ixx, 6) = U_Lf_Dn ( ixx, 6)
     1                 - ob_12 * bfxccmodes ( ixx, iyy0, izz0, 5)

          U_Lf_Dn ( ixx, 7) = U_Lf_Dn ( ixx, 7)
     1                 + ob_6 * bfyccmodes ( ixx, iyy0, izz0, 5)
     1                 + 0.25 * bfyccmodes ( ixx, iyy0, izz0, 7)

          U_Lf_Dn ( ixx, 8) = U_Lf_Dn ( ixx, 8)
     1                 + ob_6 * bfzccmodes ( ixx, iyy0, izz0, 5)
     1                 + 0.25 * bfzccmodes ( ixx, iyy0, izz0, 6)

#endif

          entropy_Lf_Dn ( ixx) = entropygr ( ixx, iyy0, izz0, 1)
     1                   + 0.5 * entropygr ( ixx, iyy0, izz0, 3)
     1                   + 0.5 * entropygr ( ixx, iyy0, izz0, 4)
     1              + dt_shift * entropygr ( ixx, iyy0, izz0, 5)

          gamma_eos_Lf_Dn ( ixx) = gamma_eos_gr ( ixx, iyy0, izz0)
          gamma_sdsp_Lf_Dn ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyy0, izz0)

! ----------
! From zone Rt_Dn == ( ixx, iyyp1, izz0) go to y = - 0.5, z = 0.5

          U_Rt_Dn ( ixx, :) = u ( ixx, iyyp1, izz0, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 4)
     1           + dt_shift * u ( ixx, iyyp1, izz0, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Rt_Dn ( ixx, 6) = U_Rt_Dn ( ixx, 6)
     1                 - ob_12 * bfxccmodes ( ixx, iyyp1, izz0, 5)

          U_Rt_Dn ( ixx, 7) = U_Rt_Dn ( ixx, 7)
     1                 + ob_6 * bfyccmodes ( ixx, iyyp1, izz0, 5)
     1                 - 0.25 * bfyccmodes ( ixx, iyyp1, izz0, 7)

          U_Rt_Dn ( ixx, 8) = U_Rt_Dn ( ixx, 8)
     1                 + ob_6 * bfzccmodes ( ixx, iyyp1, izz0, 5)
     1                 - 0.25 * bfzccmodes ( ixx, iyyp1, izz0, 6)

#endif

          entropy_Rt_Dn ( ixx) = entropygr ( ixx, iyyp1, izz0, 1)
     1                   - 0.5 * entropygr ( ixx, iyyp1, izz0, 3)
     1                   + 0.5 * entropygr ( ixx, iyyp1, izz0, 4)
     1              + dt_shift * entropygr ( ixx, iyyp1, izz0, 5)

          gamma_eos_Rt_Dn ( ixx) = gamma_eos_gr ( ixx, iyyp1, izz0)
          gamma_sdsp_Rt_Dn ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyyp1, izz0)

! ----------
! From zone Lf_Up == ( ixx, iyy0, izzp1) go to y = 0.5, z = - 0.5

          U_Lf_Up ( ixx, :) = u ( ixx, iyy0, izzp1, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyy0, izzp1, :, 4)
     1           + dt_shift * u ( ixx, iyy0, izzp1, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Lf_Up ( ixx, 6) = U_Lf_Up ( ixx, 6)
     1                 - ob_12 * bfxccmodes ( ixx, iyy0, izzp1, 5)

          U_Lf_Up ( ixx, 7) = U_Lf_Up ( ixx, 7)
     1                 + ob_6 * bfyccmodes ( ixx, iyy0, izzp1, 5)
     1                 - 0.25 * bfyccmodes ( ixx, iyy0, izzp1, 7)

          U_Lf_Up ( ixx, 8) = U_Lf_Up ( ixx, 8)
     1                 + ob_6 * bfzccmodes ( ixx, iyy0, izzp1, 5)
     1                 - 0.25 * bfzccmodes ( ixx, iyy0, izzp1, 6)

#endif

          entropy_Lf_Up ( ixx) = entropygr ( ixx, iyy0, izzp1, 1)
     1                   + 0.5 * entropygr ( ixx, iyy0, izzp1, 3)
     1                   - 0.5 * entropygr ( ixx, iyy0, izzp1, 4)
     1              + dt_shift * entropygr ( ixx, iyy0, izzp1, 5)

          gamma_eos_Lf_Up ( ixx) = gamma_eos_gr ( ixx, iyy0, izzp1)
          gamma_sdsp_Lf_Up ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyy0, izzp1)

! ----------
! From zone Rt_Up == ( ixx, iyyp1, izzp1) go to y = 0.5, z = 0.5

          U_Rt_Up ( ixx, :) = u ( ixx, iyyp1, izzp1, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 4)
     1           + dt_shift * u ( ixx, iyyp1, izzp1, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Rt_Up ( ixx, 6) = U_Rt_Up ( ixx, 6)
     1                 - ob_12 * bfxccmodes ( ixx, iyyp1, izzp1, 5)

          U_Rt_Up ( ixx, 7) = U_Rt_Up ( ixx, 7)
     1                 + ob_6 * bfyccmodes ( ixx, iyyp1, izzp1, 5)
     1                 + 0.25 * bfyccmodes ( ixx, iyyp1, izzp1, 7)

          U_Rt_Up ( ixx, 8) = U_Rt_Up ( ixx, 8)
     1                 + ob_6 * bfzccmodes ( ixx, iyyp1, izzp1, 5)
     1                 + 0.25 * bfzccmodes ( ixx, iyyp1, izzp1, 6)

#endif

          entropy_Rt_Up ( ixx) = entropygr ( ixx, iyyp1, izzp1, 1)
     1                   - 0.5 * entropygr ( ixx, iyyp1, izzp1, 3)
     1                   - 0.5 * entropygr ( ixx, iyyp1, izzp1, 4)
     1              + dt_shift * entropygr ( ixx, iyyp1, izzp1, 5)

          gamma_eos_Rt_Up ( ixx) = gamma_eos_gr ( ixx, iyyp1, izzp1)
          gamma_sdsp_Rt_Up ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyyp1, izzp1)

! ------------------------------
! Second, we use the flattener to decide how much HLL has to be blended
! into 1D And 2D HLLC/D Riemann solvers.

          blend_HLL ( ixx) = 0.25 * ( flattengr ( ixx, iyy0, izz0)
     1                              + flattengr ( ixx, iyyp1, izz0)
     1                              + flattengr ( ixx, iyy0, izzp1)
     1                              + flattengr ( ixx, iyyp1, izzp1) )

! ------------------------------

          END DO ! End of DO ixx = loop.

! ----------------------------------------------------------------------
! End of loading in the four states that surround the x-edge for the 
! multidimensional Riemann solver.
! ----------------------------------------------------------------------
! Impose boundary conditions.

! Observe that we place checks on the edge-centered indices ( iyy, izz)
! to determine whether we are at a physical boundary. However, the
! "?_indx_limits" and "bcarr_?" have to be accessed in zone-centered fashion.
! ----------------------------------------------------------------------

          DO ixx = nriembegin, nriemend

! ------------------------------
! We are at the top y-face. States that are obtained from zone "iyyp1" should
! be reset with values from zone "iyy0".

          IF ( iyy == y_indx_limits ( ixx, 2, izz0) ) THEN

! ----------

            IF ( bcarr_y ( ixx, 2, izz0) == 2) THEN

              U_Rt_Dn ( ixx, :) = U_Lf_Dn ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 2, izz0) == 3) THEN

              U_Rt_Dn ( ixx, :) = U_Lf_Dn ( ixx, :)
              U_Rt_Dn ( ixx, 3) = - U_Rt_Dn ( ixx, 3)
              U_Rt_Dn ( ixx, 6) = - U_Rt_Dn ( ixx, 6)
              U_Rt_Dn ( ixx, 8) = - U_Rt_Dn ( ixx, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixx, 2, izzp1) == 2) THEN

              U_Rt_Up ( ixx, :) = U_Lf_Up ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 2, izzp1) == 3) THEN

              U_Rt_Up ( ixx, :) = U_Lf_Up ( ixx, :)
              U_Rt_Up ( ixx, 3) = - U_Rt_Up ( ixx, 3)
              U_Rt_Up ( ixx, 6) = - U_Rt_Up ( ixx, 6)
              U_Rt_Up ( ixx, 8) = - U_Rt_Up ( ixx, 8)

            END IF

! ----------

          END IF ! End of top y-face.

! ------------------------------
! We are at the bottom y-face. States that are obtained from zone "iyy0" should
! be reset with values from zone "iyyp1".

          IF ( iyy == y_indx_limits ( ixx, 1, izz0) - 1) THEN

! ----------

            IF ( bcarr_y ( ixx, 1, izz0) == 2) THEN

              U_Lf_Dn ( ixx, :) = U_Rt_Dn ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 1, izz0) == 3) THEN

              U_Lf_Dn ( ixx, :) = U_Rt_Dn ( ixx, :)
              U_Lf_Dn ( ixx, 3) = - U_Lf_Dn ( ixx, 3)
              U_Lf_Dn ( ixx, 6) = - U_Lf_Dn ( ixx, 6)
              U_Lf_Dn ( ixx, 8) = - U_Lf_Dn ( ixx, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixx, 1, izzp1) == 2) THEN

              U_Lf_Up ( ixx, :) = U_Rt_Up ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 1, izzp1) == 3) THEN

              U_Lf_Up ( ixx, :) = U_Rt_Up ( ixx, :)
              U_Lf_Up ( ixx, 3) = - U_Lf_Up ( ixx, 3)
              U_Lf_Up ( ixx, 6) = - U_Lf_Up ( ixx, 6)
              U_Lf_Up ( ixx, 8) = - U_Lf_Up ( ixx, 8)

            END IF

! ----------

          END IF ! End of bottom y-face.

! ------------------------------
! We are at the top z-face. States that are obtained from zone "izzp1" should
! be reset with values from zone "izz0".

          IF ( izz == z_indx_limits ( ixx, iyy0, 2) ) THEN

! ----------

            IF ( bcarr_z ( ixx, iyy0, 2) == 2) THEN

              U_Lf_Up ( ixx, :) = U_Lf_Dn ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyy0, 2) == 3) THEN

              U_Lf_Up ( ixx, :) = U_Lf_Dn ( ixx, :)
              U_Lf_Up ( ixx, 4) = - U_Lf_Up ( ixx, 4)
              U_Lf_Up ( ixx, 6) = - U_Lf_Up ( ixx, 6)
              U_Lf_Up ( ixx, 7) = - U_Lf_Up ( ixx, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixx, iyyp1, 2) == 2) THEN

              U_Rt_Up ( ixx, :) = U_Rt_Dn ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyyp1, 2) == 3) THEN

              U_Rt_Up ( ixx, :) = U_Rt_Dn ( ixx, :)
              U_Rt_Up ( ixx, 4) = - U_Rt_Up ( ixx, 4)
              U_Rt_Up ( ixx, 6) = - U_Rt_Up ( ixx, 6)
              U_Rt_Up ( ixx, 7) = - U_Rt_Up ( ixx, 7)

            END IF

! ----------

          END IF ! End of top z-face.

! ------------------------------
! We are at the bottom z-face. States that are obtained from zone "izz0" should
! be reset with values from zone "izzp1".

          IF ( izz == z_indx_limits ( ixx, iyy0, 1) - 1 ) THEN

! ----------

            IF ( bcarr_z ( ixx, iyy0, 1) == 2) THEN

              U_Lf_Dn ( ixx, :) = U_Lf_Up ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyy0, 1) == 3) THEN

              U_Lf_Dn ( ixx, :) = U_Lf_Up ( ixx, :)
              U_Lf_Dn ( ixx, 4) = - U_Lf_Dn ( ixx, 4)
              U_Lf_Dn ( ixx, 6) = - U_Lf_Dn ( ixx, 6)
              U_Lf_Dn ( ixx, 7) = - U_Lf_Dn ( ixx, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixx, iyyp1, 1) == 2) THEN

              U_Rt_Dn ( ixx, :) = U_Rt_Up ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyyp1, 1) == 3) THEN

              U_Rt_Dn ( ixx, :) = U_Rt_Up ( ixx, :)
              U_Rt_Dn ( ixx, 4) = - U_Rt_Dn ( ixx, 4)
              U_Rt_Dn ( ixx, 6) = - U_Rt_Dn ( ixx, 6)
              U_Rt_Dn ( ixx, 7) = - U_Rt_Dn ( ixx, 7)

            END IF

! ----------

          END IF ! End of bottom z-face.

! ------------------------------

          END DO ! End of DO ixx = loop.

! ----------------------------------------------------------------------
! End of imposing boundary conditions.
! ----------------------------------------------------------------------
! Begin the steps for calling the Riemann solver.
! ----------------------------------------------------------------------

          DO ixx = nriembegin, nriemend

! ------------------------------

          my_blend_HLL = blend_HLL ( ixx)
          my_blend_MDHLL = blend_mdhll_to_mdhllc

          U_vrtx ( 1, :) = U_Rt_Up ( ixx, :)
          U_vrtx ( 2, :) = U_Lf_Up ( ixx, :)
          U_vrtx ( 3, :) = U_Lf_Dn ( ixx, :)
          U_vrtx ( 4, :) = U_Rt_Dn ( ixx, :)

          F_vrtx ( 1, :) = 0.0
          F_vrtx ( 2, :) = 0.0
          F_vrtx ( 3, :) = 0.0
          F_vrtx ( 4, :) = 0.0

          G_vrtx ( 1, :) = 0.0
          G_vrtx ( 2, :) = 0.0
          G_vrtx ( 3, :) = 0.0
          G_vrtx ( 4, :) = 0.0

          entropy_vrtx ( 1) = entropy_Rt_Up ( ixx)
          entropy_vrtx ( 2) = entropy_Lf_Up ( ixx)
          entropy_vrtx ( 3) = entropy_Lf_Dn ( ixx)
          entropy_vrtx ( 4) = entropy_Rt_Dn ( ixx)

          gamma_eos ( 1) = gamma_eos_Rt_Up ( ixx)
          gamma_eos ( 2) = gamma_eos_Lf_Up ( ixx)
          gamma_eos ( 3) = gamma_eos_Lf_Dn ( ixx)
          gamma_eos ( 4) = gamma_eos_Rt_Dn ( ixx)

          gamma_sdsp ( 1) = gamma_sdsp_Rt_Up ( ixx)
          gamma_sdsp ( 2) = gamma_sdsp_Lf_Up ( ixx)
          gamma_sdsp ( 3) = gamma_sdsp_Lf_Dn ( ixx)
          gamma_sdsp ( 4) = gamma_sdsp_Rt_Dn ( ixx)

! ------------------------------

          CALL FLIP_U ( nzone_vrtx, idirection, U_vrtx)


          CALL MuSIC_MHD_2DRS_1_STRUCT

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, prsfloor, rhofloor, pi,

     1     my_blend_HLL, my_blend_MDHLL,
     1     my_Prs_Tot_Star, my_Prs_Gas_Star,

     1     Vlx_ALE, Vly_ALE,

     1     U_vrtx, F_vrtx, G_vrtx, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy_vrtx,

     1     U_Star, F_Star, G_Star,
     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)


          CALL UNFLIP_FG ( nzone_vrtx, idirection, F_Star, G_Star,
     1                       F_HLLC_Global_x, G_HLLC_Global_y)

! Now "F_Star, F_HLLC_Global_x" hold y-flux and "G_Star, G_HLLC_Global_y"
! hold z-flux.

! ------------------------------

          efxgr ( ixx, iyy, izz) = 0.5 * ( G_Star ( 7) - F_Star ( 8) )

! ------------------------------
! Add y-flux to y-face that lies in the lower z-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffz_expand <= izz) .AND.
     1         ( izz <= iz1 + ioffz_expand) ) THEN

            flux_y ( ixx, iyy, izz, :) = flux_y ( ixx, iyy, izz, :)
     1                                 + ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_y ( ixx, iyy, izz, 3) =
     1      flux_y ( ixx, iyy, izz, 3) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add y-flux to y-face that lies in the upper z-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffz_expand <= izz + 1) .AND.
     1         ( izz + 1 <= iz1 + ioffz_expand) ) THEN

            flux_temp ( ixx, iyy, izz + 1, :) =
     1                                   ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1)
     1      flux_temp ( ixx, iyy, izz + 1, 3) =
     1      flux_temp ( ixx, iyy, izz + 1, 3)
     1    - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add z-flux to z-face that lies in the lower y-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffy_expand <= iyy) .AND.
     1         ( iyy <= iy1 + ioffy_expand) ) THEN

            flux_z ( ixx, iyy, izz, :) = flux_z ( ixx, iyy, izz, :)
     1                                 + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_z ( ixx, iyy, izz, 4) =
     1      flux_z ( ixx, iyy, izz, 4) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add z-flux to z-face that lies in the upper y-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffy_expand <= iyy + 1) .AND.
     1         ( iyy + 1 <= iy1 + ioffy_expand) ) THEN

            flux_z ( ixx, iyy + 1, izz, :) =
     1                             flux_z ( ixx, iyy + 1, izz, :)
     1                                 + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_z ( ixx, iyy + 1, izz, 4) =
     1    flux_z ( ixx, iyy + 1, izz, 4) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------

          END DO ! End of DO ixx = loop.

! ----------------------------------------------------------------------
! End of steps for calling the Riemann solver.
! ----------------------------------------------------------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    flux_y, flux_temp)                                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)
      DO izz = 0 - ioffz, iz1 + ioffz
        DO iyy = 0 - ioffy, iy1 + ioffy
          DO ixx = 0 - ioffx, ix1 + ioffx

            flux_y ( ixx, iyy, izz, :) = flux_y ( ixx, iyy, izz, :)
     1                              + flux_temp ( ixx, iyy, izz, :)

          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE APPLY_MULTID_RS_X_EDGE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE APPLY_MULTID_RS_Y_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_x, flux_z, flux_temp, efygr,
     1        entropy_flux_x, entropy_flux_z, entropy_flux_temp,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and evaluates the MuSIC Riemann solver at each y-edge. One resolved 
! state and two fluxes are produced, so we add in the contribution from the
! y-edge to "flux_x, flux_z". Also builds "efygr" at that edge. All fluxes must
! be zeroed before calls to this suite of three subroutines.

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

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
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1                  bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_z, flux_temp

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) :: efygr,
     1        entropy_flux_x, entropy_flux_z, entropy_flux_temp

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr

! ----------

      INTEGER, DIMENSION ( 2,
     1                     1 - ioffy : iy1 + ioffy,
     1                     1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                     2,
     1                     1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                     1 - ioffy : iy1 + ioffy,
     1                     2) ::
     1         z_indx_limits, bcarr_z

! ----------

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, ixxm1, ixx0, ixxp1, izzm1, izz0, izzp1,
     1        nriembegin, nriemend, idirection,
     1        ioffx_expand, ioffy_expand, ioffz_expand,
     1        multid_wavemodel_num, split_prsflux_re

      REAL ob_6_multid, ob_6_oned, ob_6, ob_12

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + 1 + NFLUID) ::

     1      U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,

     1      gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,
     1      gamma_eos_Rt_Dn,
     1      gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,
     1      gamma_sdsp_Rt_Dn,

     1      blend_HLL,

     1      Prs_Star,

     1      temp_01, temp_02, temp_03, temp_04, temp_05, dx_use, dz_use

! ----------------------------------------------------------------------

! These are the arrays that actually go through the vertex-based RS.

      INTEGER nzone_vrtx

      REAL

! Amount of HLL we will have to blend in for stabilization.
     1     my_blend_HLL, my_blend_MDHLL,

! Velocity of vertex; total and gas pressure in resolved state.
     1     Vlx_ALE, Vly_ALE,
     1     my_Prs_Tot_Star, my_Prs_Gas_Star

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables that come together at a zone edge.
! We also come in with their primitive variables and associated
! x- and y-fluxes.
     1      U_vrtx, F_vrtx, G_vrtx,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

! entropy variables that go through RS interface.
     1      entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,

! Components of the unit normal.
     1      nx, ny

! Final, strongly-interacting state and Global x- and y-fluxes.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) ::
     1      U_Star, F_Star, G_Star

! ----------------------------------------------------------------------

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

      ob_6 = 1.0 / 6.0
      ob_12 = 1.0 / 12.0

! ----------

! We are interested in the fluxes for a non-moving mesh.
      Vlx_ALE = 0.0
      Vly_ALE = 0.0

! ----------

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

! In true 2d, only the y-edges contribute to a problem whose variation
! is entirely in the xz-plane. If this is truly a 2d problem with no variation
! in the x-direction then only two edges contribute fluxes to a face
! (instead of four), so "ob_6_multid and ob_6_oned" are increased.

        multid_wavemodel_num = 1  ! 1 includes supersonic; 2 is fully subsonic.

! "4 * ob_6_multid + 4 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux
        ob_6_oned = 0.25 * ( 1.0 - 4.0 * blend_multid_flux)

#if ( TRUE_2D == 1)
        IF ( ioffy .NE. 0) RETURN
! "2 * ob_6_multid + 2 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux
        ob_6_oned = 0.5 * ( 1.0 - 2.0 * blend_multid_flux)
#endif

! ----------

      idirection = 2

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in APPLY_MULTID_RS_Y_EDGE : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, ixx0, ixxp1,          &
!$OMP&    izzm1, izz0, izzp1, nriembegin, nriemend,                     &
!$OMP&    U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn,                           &
!$OMP&    entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,   &

!$OMP&    gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,            &
!$OMP&    gamma_eos_Rt_Dn,                                              &
!$OMP&    gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,         &
!$OMP&    gamma_sdsp_Rt_Dn,                                             &

!$OMP&    blend_HLL, my_Prs_Gas_Star,                                   &
!$OMP&    my_blend_HLL, my_blend_MDHLL, my_Prs_Tot_Star,                &
!$OMP&    U_vrtx,  F_vrtx, G_vrtx,                                      &
!$OMP&    gamma_eos, gamma_sdsp, U_Star, F_Star, G_Star,                &
!$OMP&    U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,     &
!$OMP&    entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,             &

!$OMP&    temp_01, temp_02, temp_03, temp_04, temp_05)                  &

!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    ob_6_multid, ob_6_oned, blend_mdhll_to_mdhllc,                &
!$OMP&    dt_shift, smallnum, prsfloor, rhofloor, pi, idirection,       &

!$OMP&    nzone_vrtx, Vlx_ALE, Vly_ALE, nx, ny, ob_6, ob_12,            &
!$OMP&    multid_wavemodel_num, split_prsflux_re,                       &

!$OMP&    x_indx_limits, y_indx_limits, z_indx_limits,                  &
!$OMP&    bcarr_x, bcarr_y, bcarr_z, dx_save, dy_save, dz_save,         &
!$OMP&    u, entropygr, gamma_eos_gr, gamma_soundspeed_gr, flattengr,   &
!$OMP&    bfxccmodes, bfyccmodes, bfzccmodes,                           &
!$OMP&    efygr, flux_x, flux_z, flux_temp,                             &
!$OMP&    entropy_flux_x, entropy_flux_z, entropy_flux_temp)            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz_expand, iz1 + ioffz_expand

        DO ixx = 0 - ioffx_expand, ix1 + ioffx_expand

! ----------------------------------------------------------------------
! Indices ( ixx0, izz0) always map to an interior point on the mesh including
! its ghost zones. This is needed because ( ixx, izz) map to the edges.

        izzm1 = MAX0 ( izz - 1, 0 - ioffz)
        izz0 = MAX0 ( izz, 1 - ioffz)
        izzp1 = MIN0 ( izz + 1, iz1 + ioffz)

        ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)
        ixx0 = MAX0 ( ixx, 1 - ioffx)
        ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)

        nriembegin = y_indx_limits ( ixx0, 1, izz0) - ioffy_expand
        nriemend = y_indx_limits ( ixx0, 2, izz0) + ioffy_expand

! ----------------------------------------------------------------------
! Load in the four states that surround the y-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------

          DO iyy = nriembegin, nriemend

! ------------------------------
! First, we build the four states.

! ----------
! From zone Lf_Dn == ( ixx0, iyy, izz0) go to x = 0.5, z = 0.5

          U_Lf_Dn ( iyy, :) = u ( ixx0, iyy, izz0, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 4)
     1           + dt_shift * u ( ixx0, iyy, izz0, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Lf_Dn ( iyy, 6) = U_Lf_Dn ( iyy, 6)
     1                 + ob_6 * bfxccmodes ( ixx0, iyy, izz0, 5)
     1                 + 0.25 * bfxccmodes ( ixx0, iyy, izz0, 7)

          U_Lf_Dn ( iyy, 7) = U_Lf_Dn ( iyy, 7)
     1                 - ob_12 * bfyccmodes ( ixx0, iyy, izz0, 5)

          U_Lf_Dn ( iyy, 8) = U_Lf_Dn ( iyy, 8)
     1                 + ob_6 * bfzccmodes ( ixx0, iyy, izz0, 5)
     1                 + 0.25 * bfzccmodes ( ixx0, iyy, izz0, 7)

#endif

          entropy_Lf_Dn ( iyy) = entropygr ( ixx0, iyy, izz0, 1)
     1                   + 0.5 * entropygr ( ixx0, iyy, izz0, 2)
     1                   + 0.5 * entropygr ( ixx0, iyy, izz0, 4)
     1              + dt_shift * entropygr ( ixx0, iyy, izz0, 5)

          gamma_eos_Lf_Dn ( iyy) = gamma_eos_gr ( ixx0, iyy, izz0)
          gamma_sdsp_Lf_Dn ( iyy) =
     1                      gamma_soundspeed_gr ( ixx0, iyy, izz0)

! ----------
! From zone Lf_Up == ( ixxp1, iyy, izz0) go to x = - 0.5, z = 0.5

          U_Lf_Up ( iyy, :) = u ( ixxp1, iyy, izz0, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 4)
     1           + dt_shift * u ( ixxp1, iyy, izz0, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Lf_Up ( iyy, 6) = U_Lf_Up ( iyy, 6)
     1                 + ob_6 * bfxccmodes ( ixxp1, iyy, izz0, 5)
     1                 - 0.25 * bfxccmodes ( ixxp1, iyy, izz0, 7)

          U_Lf_Up ( iyy, 7) = U_Lf_Up ( iyy, 7)
     1                 - ob_12 * bfyccmodes ( ixxp1, iyy, izz0, 5)

          U_Lf_Up ( iyy, 8) = U_Lf_Up ( iyy, 8)
     1                 + ob_6 * bfzccmodes ( ixxp1, iyy, izz0, 5)
     1                 - 0.25 * bfzccmodes ( ixxp1, iyy, izz0, 7)

#endif

          entropy_Lf_Up ( iyy) = entropygr ( ixxp1, iyy, izz0, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyy, izz0, 2)
     1                   + 0.5 * entropygr ( ixxp1, iyy, izz0, 4)
     1              + dt_shift * entropygr ( ixxp1, iyy, izz0, 5)

          gamma_eos_Lf_Up ( iyy) = gamma_eos_gr ( ixxp1, iyy, izz0)
          gamma_sdsp_Lf_Up ( iyy) =
     1                      gamma_soundspeed_gr ( ixxp1, iyy, izz0)

! ----------
! From zone Rt_Dn == ( ixx0, iyy, izzp1) go to x = 0.5, z = - 0.5

          U_Rt_Dn ( iyy, :) = u ( ixx0, iyy, izzp1, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixx0, iyy, izzp1, :, 4)
     1           + dt_shift * u ( ixx0, iyy, izzp1, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Rt_Dn ( iyy, 6) = U_Rt_Dn ( iyy, 6)
     1                 + ob_6 * bfxccmodes ( ixx0, iyy, izzp1, 5)
     1                 - 0.25 * bfxccmodes ( ixx0, iyy, izzp1, 7)

          U_Rt_Dn ( iyy, 7) = U_Rt_Dn ( iyy, 7)
     1                 - ob_12 * bfyccmodes ( ixx0, iyy, izzp1, 5)

          U_Rt_Dn ( iyy, 8) = U_Rt_Dn ( iyy, 8)
     1                 + ob_6 * bfzccmodes ( ixx0, iyy, izzp1, 5)
     1                 - 0.25 * bfzccmodes ( ixx0, iyy, izzp1, 7)

#endif

          entropy_Rt_Dn ( iyy) = entropygr ( ixx0, iyy, izzp1, 1)
     1                   + 0.5 * entropygr ( ixx0, iyy, izzp1, 2)
     1                   - 0.5 * entropygr ( ixx0, iyy, izzp1, 4)
     1              + dt_shift * entropygr ( ixx0, iyy, izzp1, 5)

          gamma_eos_Rt_Dn ( iyy) = gamma_eos_gr ( ixx0, iyy, izzp1)
          gamma_sdsp_Rt_Dn ( iyy) =
     1                      gamma_soundspeed_gr ( ixx0, iyy, izzp1)

! ----------
! From zone Rt_Up == ( ixxp1, iyy, izzp1) go to x = - 0.5, z = - 0.5

          U_Rt_Up ( iyy, :) = u ( ixxp1, iyy, izzp1, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 4)
     1           + dt_shift * u ( ixxp1, iyy, izzp1, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Rt_Up ( iyy, 6) = U_Rt_Up ( iyy, 6)
     1                 + ob_6 * bfxccmodes ( ixxp1, iyy, izzp1, 5)
     1                 + 0.25 * bfxccmodes ( ixxp1, iyy, izzp1, 7)

          U_Rt_Up ( iyy, 7) = U_Rt_Up ( iyy, 7)
     1                 - ob_12 * bfyccmodes ( ixxp1, iyy, izzp1, 5)

          U_Rt_Up ( iyy, 8) = U_Rt_Up ( iyy, 8)
     1                 + ob_6 * bfzccmodes ( ixxp1, iyy, izzp1, 5)
     1                 + 0.25 * bfzccmodes ( ixxp1, iyy, izzp1, 7)

#endif

          entropy_Rt_Up ( iyy) = entropygr ( ixxp1, iyy, izzp1, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyy, izzp1, 2)
     1                   - 0.5 * entropygr ( ixxp1, iyy, izzp1, 4)
     1              + dt_shift * entropygr ( ixxp1, iyy, izzp1, 5)

          gamma_eos_Rt_Up ( iyy) = gamma_eos_gr ( ixxp1, iyy, izzp1)
          gamma_sdsp_Rt_Up ( iyy) =
     1                      gamma_soundspeed_gr ( ixxp1, iyy, izzp1)

! ------------------------------
! Second, we use the flattener to decide how much HLL has to be blended
! into 1D And 2D HLLC/D Riemann solvers.

          blend_HLL ( iyy) = 0.25 * ( flattengr ( ixx0, iyy, izz0)
     1                              + flattengr ( ixxp1, iyy, izz0)
     1                              + flattengr ( ixx0, iyy, izzp1)
     1                              + flattengr ( ixxp1, iyy, izzp1) )

! ------------------------------

          END DO ! End of DO iyy = loop.

! ----------------------------------------------------------------------
! End of loading in the four states that surround the y-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------
! Impose boundary conditions.

! Observe that we place checks on the edge-centered indices ( ixx, izz) 
! to determine whether we are at a physical boundary. However, the
! "?_indx_limits" and "bcarr_?" have to be accessed in zone-centered fashion.
! ----------------------------------------------------------------------

          DO iyy = nriembegin, nriemend

! ------------------------------
! We are at the top x-face. States that are obtained from zone "ixxp1" should
! be reset with values from zone "ixx0".

          IF ( ixx == x_indx_limits ( 2, iyy, izz0) ) THEN

! ----------

            IF ( bcarr_x ( 2, iyy, izz0) == 2) THEN

              U_Lf_Up ( iyy, :) = U_Lf_Dn ( iyy, :)

            ELSE IF ( bcarr_x ( 2, iyy, izz0) == 3) THEN

              U_Lf_Up ( iyy, :) = U_Lf_Dn ( iyy, :)
              U_Lf_Up ( iyy, 2) = - U_Lf_Up ( iyy, 2)
              U_Lf_Up ( iyy, 7) = - U_Lf_Up ( iyy, 7)
              U_Lf_Up ( iyy, 8) = - U_Lf_Up ( iyy, 8)

            END IF

! ----------

            IF ( bcarr_x ( 2, iyy, izzp1) == 2) THEN

              U_Rt_Up ( iyy, :) = U_Rt_Dn ( iyy, :)

            ELSE IF ( bcarr_x ( 2, iyy, izzp1) == 3) THEN

              U_Rt_Up ( iyy, :) = U_Rt_Dn ( iyy, :)
              U_Rt_Up ( iyy, 2) = - U_Rt_Up ( iyy, 2)
              U_Rt_Up ( iyy, 7) = - U_Rt_Up ( iyy, 7)
              U_Rt_Up ( iyy, 8) = - U_Rt_Up ( iyy, 8)

            END IF

! ----------

          END IF ! End of top x-face.

! ------------------------------
! We are at the bottom x-face. States that are obtained from zone "ixx0" should
! be reset with values from zone "ixxp1".

          IF ( ixx == x_indx_limits ( 1, iyy, izz0) - 1) THEN

! ----------

            IF ( bcarr_x ( 1, iyy, izz0) == 2) THEN

              U_Lf_Dn ( iyy, :) = U_Lf_Up ( iyy, :)

            ELSE IF ( bcarr_x ( 1, iyy, izz0) == 3) THEN

              U_Lf_Dn ( iyy, :) = U_Lf_Up ( iyy, :)
              U_Lf_Dn ( iyy, 2) = - U_Lf_Dn ( iyy, 2)
              U_Lf_Dn ( iyy, 7) = - U_Lf_Dn ( iyy, 7)
              U_Lf_Dn ( iyy, 8) = - U_Lf_Dn ( iyy, 8)

            END IF

! ----------

            IF ( bcarr_x ( 1, iyy, izzp1) == 2) THEN

              U_Rt_Dn ( iyy, :) = U_Rt_Up ( iyy, :)

            ELSE IF ( bcarr_x ( 1, iyy, izzp1) == 3) THEN

              U_Rt_Dn ( iyy, :) = U_Rt_Up ( iyy, :)
              U_Rt_Dn ( iyy, 2) = - U_Rt_Dn ( iyy, 2)
              U_Rt_Dn ( iyy, 7) = - U_Rt_Dn ( iyy, 7)
              U_Rt_Dn ( iyy, 8) = - U_Rt_Dn ( iyy, 8)

            END IF

! ----------

          END IF ! End of bottom x-face.

! ------------------------------
! We are at the top z-face. States that are obtained from zone "izzp1" should
! be reset with values from zone "izz0".

          IF ( izz == z_indx_limits ( ixx0, iyy, 2) ) THEN

! ----------

            IF ( bcarr_z ( ixx0, iyy, 2) == 2) THEN

              U_Rt_Dn ( iyy, :) = U_Lf_Dn ( iyy, :)

            ELSE IF ( bcarr_z ( ixx0, iyy, 2) == 3) THEN

              U_Rt_Dn ( iyy, :) = U_Lf_Dn ( iyy, :)
              U_Rt_Dn ( iyy, 4) = - U_Rt_Dn ( iyy, 4)
              U_Rt_Dn ( iyy, 6) = - U_Rt_Dn ( iyy, 6)
              U_Rt_Dn ( iyy, 7) = - U_Rt_Dn ( iyy, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixxp1, iyy, 2) == 2) THEN

              U_Rt_Up ( iyy, :) = U_Lf_Up ( iyy, :)

            ELSE IF ( bcarr_z ( ixxp1, iyy, 2) == 3) THEN

              U_Rt_Up ( iyy, :) = U_Lf_Up ( iyy, :)
              U_Rt_Up ( iyy, 4) = - U_Rt_Up ( iyy, 4)
              U_Rt_Up ( iyy, 6) = - U_Rt_Up ( iyy, 6)
              U_Rt_Up ( iyy, 7) = - U_Rt_Up ( iyy, 7)

            END IF

! ----------

          END IF ! End of top z-face.

! ------------------------------
! We are at the bottom z-face. States that are obtained from zone "izz0" should
! be reset with values from zone "izzp1".

          IF ( izz == z_indx_limits ( ixx0, iyy, 1) - 1) THEN

! ----------

            IF ( bcarr_z ( ixx0, iyy, 1) == 2) THEN

              U_Lf_Dn ( iyy, :) = U_Rt_Dn ( iyy, :)

            ELSE IF ( bcarr_z ( ixx0, iyy, 1) == 3) THEN

              U_Lf_Dn ( iyy, :) = U_Rt_Dn ( iyy, :)
              U_Lf_Dn ( iyy, 4) = - U_Lf_Dn ( iyy, 4)
              U_Lf_Dn ( iyy, 6) = - U_Lf_Dn ( iyy, 6)
              U_Lf_Dn ( iyy, 7) = - U_Lf_Dn ( iyy, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixxp1, iyy, 1) == 2) THEN

              U_Lf_Up ( iyy, :) = U_Rt_Up ( iyy, :)

            ELSE IF ( bcarr_z ( ixxp1, iyy, 1) == 3) THEN

              U_Lf_Up ( iyy, :) = U_Rt_Up ( iyy, :)
              U_Lf_Up ( iyy, 4) = - U_Lf_Up ( iyy, 4)
              U_Lf_Up ( iyy, 6) = - U_Lf_Up ( iyy, 6)
              U_Lf_Up ( iyy, 7) = - U_Lf_Up ( iyy, 7)

            END IF

! ----------

          END IF ! End of bottom z-face.

! ------------------------------

          END DO ! End of DO iyy = loop.

! ----------------------------------------------------------------------
! End of imposing boundary conditions.
! ----------------------------------------------------------------------
! Begin the steps for calling the Riemann solver.
! ----------------------------------------------------------------------

          DO iyy = nriembegin, nriemend

! ------------------------------

          my_blend_HLL = blend_HLL ( iyy)
          my_blend_MDHLL = blend_mdhll_to_mdhllc

          U_vrtx ( 1, :) = U_Rt_Up ( iyy, :)
          U_vrtx ( 2, :) = U_Lf_Up ( iyy, :)
          U_vrtx ( 3, :) = U_Lf_Dn ( iyy, :)
          U_vrtx ( 4, :) = U_Rt_Dn ( iyy, :)

          F_vrtx ( 1, :) = 0.0
          F_vrtx ( 2, :) = 0.0
          F_vrtx ( 3, :) = 0.0
          F_vrtx ( 4, :) = 0.0

          G_vrtx ( 1, :) = 0.0
          G_vrtx ( 2, :) = 0.0
          G_vrtx ( 3, :) = 0.0
          G_vrtx ( 4, :) = 0.0

          entropy_vrtx ( 1) = entropy_Rt_Up ( iyy)
          entropy_vrtx ( 2) = entropy_Lf_Up ( iyy)
          entropy_vrtx ( 3) = entropy_Lf_Dn ( iyy)
          entropy_vrtx ( 4) = entropy_Rt_Dn ( iyy)

          gamma_eos ( 1) = gamma_eos_Rt_Up ( iyy)
          gamma_eos ( 2) = gamma_eos_Lf_Up ( iyy)
          gamma_eos ( 3) = gamma_eos_Lf_Dn ( iyy)
          gamma_eos ( 4) = gamma_eos_Rt_Dn ( iyy)

          gamma_sdsp ( 1) = gamma_sdsp_Rt_Up ( iyy)
          gamma_sdsp ( 2) = gamma_sdsp_Lf_Up ( iyy)
          gamma_sdsp ( 3) = gamma_sdsp_Lf_Dn ( iyy)
          gamma_sdsp ( 4) = gamma_sdsp_Rt_Dn ( iyy)

! ------------------------------

          CALL FLIP_U ( nzone_vrtx, idirection, U_vrtx)


          CALL MuSIC_MHD_2DRS_1_STRUCT

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, prsfloor, rhofloor, pi,

     1     my_blend_HLL, my_blend_MDHLL,
     1     my_Prs_Tot_Star, my_Prs_Gas_Star,

     1     Vlx_ALE, Vly_ALE,

     1     U_vrtx, F_vrtx, G_vrtx, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy_vrtx,

     1     U_Star, F_Star, G_Star,
     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)


          CALL UNFLIP_FG ( nzone_vrtx, idirection, F_Star, G_Star,
     1                       F_HLLC_Global_x, G_HLLC_Global_y)

! Now "F_Star, F_HLLC_Global_x" hold z-flux and "G_Star, G_HLLC_Global_y"
! hold x-flux.

! ------------------------------

          efygr ( ixx, iyy, izz) = 0.5 * ( G_Star ( 8) - F_Star ( 6) )

! ------------------------------
! Add x-flux to x-face that lies in the lower z-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffz_expand <= izz) .AND.
     1         ( izz <= iz1 + ioffz_expand) ) THEN

            flux_x ( ixx, iyy, izz, :) = flux_x ( ixx, iyy, izz, :)
     1                                 + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_x ( ixx, iyy, izz, 2) =
     1      flux_x ( ixx, iyy, izz, 2) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add x-flux to x-face that lies in the upper z-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffz_expand <= izz + 1) .AND.
     1         ( izz + 1 <= iz1 + ioffz_expand) ) THEN

            flux_temp ( ixx, iyy, izz + 1, :) =
     1                                   ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1)
     1      flux_temp ( ixx, iyy, izz + 1, 2) =
     1      flux_temp ( ixx, iyy, izz + 1, 2)
     1    - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add z-flux to z-face that lies in the lower x-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffx_expand <= ixx) .AND.
     1         ( ixx <= ix1 + ioffx_expand) ) THEN

            flux_z ( ixx, iyy, izz, :) = flux_z ( ixx, iyy, izz, :)
     1                                 + ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_z ( ixx, iyy, izz, 4) =
     1      flux_z ( ixx, iyy, izz, 4) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add z-flux to z-face that lies in the upper x-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffx_expand <= ixx + 1) .AND.
     1         ( ixx + 1 <= ix1 + ioffx_expand) ) THEN

            flux_z ( ixx + 1, iyy, izz, :) =
     1                                flux_z ( ixx + 1, iyy, izz, :)
     1                                 + ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_z ( ixx + 1, iyy, izz, 4) =
     1    flux_z ( ixx + 1, iyy, izz, 4) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------

          END DO ! End of DO iyy = loop.

! ----------------------------------------------------------------------
! End of steps for calling the Riemann solver.
! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    flux_x, flux_temp)                                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)
      DO izz = 0 - ioffz, iz1 + ioffz
        DO iyy = 0 - ioffy, iy1 + ioffy
          DO ixx = 0 - ioffx, ix1 + ioffx

            flux_x ( ixx, iyy, izz, :) = flux_x ( ixx, iyy, izz, :)
     1                              + flux_temp ( ixx, iyy, izz, :)

          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE APPLY_MULTID_RS_Y_EDGE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE APPLY_MULTID_RS_Z_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_x, flux_y, flux_temp, efzgr,
     1        entropy_flux_x, entropy_flux_y, entropy_flux_temp,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and evaluates the MuSIC Riemann solver at each z-edge. One resolved 
! state and two fluxes are produced, so we add in the contribution from the
! z-edge to "flux_x, flux_y". Also builds "efzgr" at that edge. All fluxes must
! be zeroed before calls to this suite of three subroutines.

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

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
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1                  bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_y, flux_temp

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) :: efzgr,
     1        entropy_flux_x, entropy_flux_y, entropy_flux_temp

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr

! ----------

      INTEGER, DIMENSION ( 2,
     1                     1 - ioffy : iy1 + ioffy,
     1                     1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                     2,
     1                     1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                     1 - ioffy : iy1 + ioffy,
     1                     2) ::
     1         z_indx_limits, bcarr_z

! ----------

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, ixxm1, ixx0, ixxp1, iyym1, iyy0, iyyp1,
     1        nriembegin, nriemend, idirection,
     1        ioffx_expand, ioffy_expand, ioffz_expand,
     1        multid_wavemodel_num, split_prsflux_re

      REAL ob_6_multid, ob_6_oned, ob_6, ob_12

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + 1 + NFLUID) ::

     1      U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1      entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,

     1      gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,
     1      gamma_eos_Rt_Dn,
     1      gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,
     1      gamma_sdsp_Rt_Dn,

     1      blend_HLL,

     1      Prs_Star,

     1      temp_01, temp_02, temp_03, temp_04, temp_05, dx_use, dy_use

! ----------------------------------------------------------------------

! These are the arrays that actually go through the vertex-based RS.

      INTEGER nzone_vrtx

      REAL

! Amount of HLL we will have to blend in for stabilization.
     1     my_blend_HLL, my_blend_MDHLL,

! Velocity of vertex; total and gas pressure in resolved state.
     1     Vlx_ALE, Vly_ALE,
     1     my_Prs_Tot_Star, my_Prs_Gas_Star

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables that come together at a zone edge.
! We also come in with their primitive variables and associated
! x- and y-fluxes.
     1      U_vrtx, F_vrtx, G_vrtx,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

! entropy variables that go through RS interface.
     1      entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,

! Components of the unit normal.
     1      nx, ny

! Final, strongly-interacting state and Global x- and y-fluxes.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) ::
     1      U_Star, F_Star, G_Star

! ----------------------------------------------------------------------

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

      ob_6 = 1.0 / 6.0
      ob_12 = 1.0 / 12.0

! ----------

! We are interested in the fluxes for a non-moving mesh.
      Vlx_ALE = 0.0
      Vly_ALE = 0.0

! ----------

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

! In true 2d, only the z-edges contribute to a problem whose variation
! is entirely in the xy-plane. If this is truly a 2d problem with no variation
! in the x-direction then only two edges contribute fluxes to a face
! (instead of four), so "ob_6_multid and ob_6_oned" are increased.

        multid_wavemodel_num = 1  ! 1 includes supersonic; 2 is fully subsonic.

! "4 * ob_6_multid + 4 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux
        ob_6_oned = 0.25 * ( 1.0 - 4.0 * blend_multid_flux)

#if ( TRUE_2D == 1)
        IF ( ioffz .NE. 0) RETURN
! "2 * ob_6_multid + 2 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux
        ob_6_oned = 0.5 * ( 1.0 - 2.0 * blend_multid_flux)
#endif

! ----------

      idirection = 3

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in APPLY_MULTID_RS_Z_EDGE : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, ixx0, ixxp1,          &
!$OMP&    iyym1, iyy0, iyyp1, nriembegin, nriemend,                     &
!$OMP&    U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn,                           &
!$OMP&    entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,   &

!$OMP&    gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,            &
!$OMP&    gamma_eos_Rt_Dn,                                              &
!$OMP&    gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,         &
!$OMP&    gamma_sdsp_Rt_Dn,                                             &

!$OMP&    blend_HLL, my_Prs_Gas_Star,                                   &
!$OMP&    my_blend_HLL, my_blend_MDHLL, my_Prs_Tot_Star,                &
!$OMP&    U_vrtx, F_vrtx, G_vrtx,                                       &
!$OMP&    gamma_eos, gamma_sdsp, U_Star, F_Star, G_Star,                &
!$OMP&    U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,     &
!$OMP&    entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,             &

!$OMP&    temp_01, temp_02, temp_03, temp_04, temp_05)                  &

!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    ob_6_multid, ob_6_oned, blend_mdhll_to_mdhllc,                &
!$OMP&    dt_shift, smallnum, prsfloor, rhofloor, pi, idirection,       &

!$OMP&    nzone_vrtx, Vlx_ALE, Vly_ALE, nx, ny, ob_6, ob_12,            &
!$OMP&    multid_wavemodel_num, split_prsflux_re,                       &

!$OMP&    x_indx_limits, y_indx_limits, z_indx_limits,                  &
!$OMP&    bcarr_x, bcarr_y, bcarr_z, dx_save, dy_save, dz_save,         &
!$OMP&    u, entropygr, gamma_eos_gr, gamma_soundspeed_gr, flattengr,   &
!$OMP&    bfxccmodes, bfyccmodes, bfzccmodes,                           &
!$OMP&    efzgr, flux_x, flux_y, flux_temp,                             &
!$OMP&    entropy_flux_x, entropy_flux_y, entropy_flux_temp)            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 0 - ioffy_expand, iy1 + ioffy_expand

        DO ixx = 0 - ioffx_expand, ix1 + ioffx_expand

! ----------------------------------------------------------------------
! Indices ( ixx0, iyy0) always map to an interior point on the mesh including
! its ghost zones. This is needed because ( ixx, iyy) map to the edges.

        iyym1 = MAX0 ( iyy - 1, 0 - ioffy)
        iyy0 = MAX0 ( iyy, 1 - ioffy)
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)

        ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)
        ixx0 = MAX0 ( ixx, 1 - ioffx)
        ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)

        nriembegin = z_indx_limits ( ixx0, iyy0, 1) - ioffz_expand
        nriemend = z_indx_limits ( ixx0, iyy0, 2) + ioffz_expand

! ----------------------------------------------------------------------
! Load in the four states that surround the z-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------

          DO izz = nriembegin, nriemend

! ------------------------------
! First, we build the four states.

! ----------
! From zone Lf_Dn == ( ixx0, iyy0, izz) go to x = 0.5, y = 0.5

          U_Lf_Dn ( izz, :) = u ( ixx0, iyy0, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 3)
     1           + dt_shift * u ( ixx0, iyy0, izz, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Lf_Dn ( izz, 6) = U_Lf_Dn ( izz, 6)
     1                 + ob_6 * bfxccmodes ( ixx0, iyy0, izz, 5)
     1                 + 0.25 * bfxccmodes ( ixx0, iyy0, izz, 6)

          U_Lf_Dn ( izz, 7) = U_Lf_Dn ( izz, 7)
     1                 + ob_6 * bfyccmodes ( ixx0, iyy0, izz, 5)
     1                 + 0.25 * bfyccmodes ( ixx0, iyy0, izz, 6)

          U_Lf_Dn ( izz, 8) = U_Lf_Dn ( izz, 8)
     1                 - ob_12 * bfzccmodes ( ixx0, iyy0, izz, 5)

#endif

          entropy_Lf_Dn ( izz) = entropygr ( ixx0, iyy0, izz, 1)
     1                   + 0.5 * entropygr ( ixx0, iyy0, izz, 2)
     1                   + 0.5 * entropygr ( ixx0, iyy0, izz, 3)
     1              + dt_shift * entropygr ( ixx0, iyy0, izz, 5)

          gamma_eos_Lf_Dn ( izz) = gamma_eos_gr ( ixx0, iyy0, izz)
          gamma_sdsp_Lf_Dn ( izz) =
     1                      gamma_soundspeed_gr ( ixx0, iyy0, izz)

! ----------
! From zone Rt_Dn == ( ixxp1, iyy0, izz) go to x = - 0.5, y = 0.5

          U_Rt_Dn ( izz, :) = u ( ixxp1, iyy0, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 3)
     1           + dt_shift * u ( ixxp1, iyy0, izz, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Rt_Dn ( izz, 6) = U_Rt_Dn ( izz, 6)
     1                 + ob_6 * bfxccmodes ( ixxp1, iyy0, izz, 5)
     1                 - 0.25 * bfxccmodes ( ixxp1, iyy0, izz, 6)

          U_Rt_Dn ( izz, 7) = U_Rt_Dn ( izz, 7)
     1                 + ob_6 * bfyccmodes ( ixxp1, iyy0, izz, 5)
     1                 - 0.25 * bfyccmodes ( ixxp1, iyy0, izz, 6)

          U_Rt_Dn ( izz, 8) = U_Rt_Dn ( izz, 8)
     1                 - ob_12 * bfzccmodes ( ixxp1, iyy0, izz, 5)

#endif

          entropy_Rt_Dn ( izz) = entropygr ( ixxp1, iyy0, izz, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyy0, izz, 2)
     1                   + 0.5 * entropygr ( ixxp1, iyy0, izz, 3)
     1              + dt_shift * entropygr ( ixxp1, iyy0, izz, 5)

          gamma_eos_Rt_Dn ( izz) = gamma_eos_gr ( ixxp1, iyy0, izz)
          gamma_sdsp_Rt_Dn ( izz) =
     1                      gamma_soundspeed_gr ( ixxp1, iyy0, izz)

! ----------
! From zone Lf_Up == ( ixx0, iyyp1, izz) go to x = 0.5, y = - 0.5

          U_Lf_Up ( izz, :) = u ( ixx0, iyyp1, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixx0, iyyp1, izz, :, 3)
     1           + dt_shift * u ( ixx0, iyyp1, izz, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Lf_Up ( izz, 6) = U_Lf_Up ( izz, 6)
     1                 + ob_6 * bfxccmodes ( ixx0, iyyp1, izz, 5)
     1                 - 0.25 * bfxccmodes ( ixx0, iyyp1, izz, 6)

          U_Lf_Up ( izz, 7) = U_Lf_Up ( izz, 7)
     1                 + ob_6 * bfyccmodes ( ixx0, iyyp1, izz, 5)
     1                 - 0.25 * bfyccmodes ( ixx0, iyyp1, izz, 6)

          U_Lf_Up ( izz, 8) = U_Lf_Up ( izz, 8)
     1                 - ob_12 * bfzccmodes ( ixx0, iyyp1, izz, 5)

#endif

          entropy_Lf_Up ( izz) = entropygr ( ixx0, iyyp1, izz, 1)
     1                   + 0.5 * entropygr ( ixx0, iyyp1, izz, 2)
     1                   - 0.5 * entropygr ( ixx0, iyyp1, izz, 3)
     1              + dt_shift * entropygr ( ixx0, iyyp1, izz, 5)

          gamma_eos_Lf_Up ( izz) = gamma_eos_gr ( ixx0, iyyp1, izz)
          gamma_sdsp_Lf_Up ( izz) =
     1                      gamma_soundspeed_gr ( ixx0, iyyp1, izz)

! ----------
! From zone Rt_Up == ( ixxp1, iyyp1, izz) go to x = - 0.5, y = - 0.5

          U_Rt_Up ( izz, :) = u ( ixxp1, iyyp1, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 3)
     1           + dt_shift * u ( ixxp1, iyyp1, izz, :, 5)

#if ( DIV_FREE_RECONSTRUCT == 1)

          U_Rt_Up ( izz, 6) = U_Rt_Up ( izz, 6)
     1                 + ob_6 * bfxccmodes ( ixxp1, iyyp1, izz, 5)
     1                 + 0.25 * bfxccmodes ( ixxp1, iyyp1, izz, 6)

          U_Rt_Up ( izz, 7) = U_Rt_Up ( izz, 7)
     1                 + ob_6 * bfyccmodes ( ixxp1, iyyp1, izz, 5)
     1                 + 0.25 * bfyccmodes ( ixxp1, iyyp1, izz, 6)

          U_Rt_Up ( izz, 8) = U_Rt_Up ( izz, 8)
     1                 - ob_12 * bfzccmodes ( ixxp1, iyyp1, izz, 5)

#endif

          entropy_Rt_Up ( izz) = entropygr ( ixxp1, iyyp1, izz, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyyp1, izz, 2)
     1                   - 0.5 * entropygr ( ixxp1, iyyp1, izz, 3)
     1              + dt_shift * entropygr ( ixxp1, iyyp1, izz, 5)

          gamma_eos_Rt_Up ( izz) = gamma_eos_gr ( ixxp1, iyyp1, izz)
          gamma_sdsp_Rt_Up ( izz) =
     1                      gamma_soundspeed_gr ( ixxp1, iyyp1, izz)

! ------------------------------
! Second, we use the flattener to decide how much HLL has to be blended
! into 1D And 2D HLLC/D Riemann solvers.

          blend_HLL ( izz) = 0.25 * ( flattengr ( ixx0, iyy0, izz)
     1                              + flattengr ( ixxp1, iyy0, izz)
     1                              + flattengr ( ixx0, iyyp1, izz)
     1                              + flattengr ( ixxp1, iyyp1, izz) )

! ------------------------------

          END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------
! End of loading in the four states that surround the z-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------
! Impose boundary conditions.

! Observe that we place checks on the edge-centered indices ( ixx, iyy)
! to determine whether we are at a physical boundary. However, the
! "?_indx_limits" and "bcarr_?" have to be accessed in zone-centered fashion.
! ----------------------------------------------------------------------

          DO izz = nriembegin, nriemend

! ------------------------------
! We are at the top x-face. States that are obtained from zone "ixxp1" should
! be reset with values from zone "ixx0".

          IF ( ixx == x_indx_limits ( 2, iyy0, izz) ) THEN

! ----------

            IF ( bcarr_x ( 2, iyy0, izz) == 2) THEN

              U_Rt_Dn ( izz, :) = U_Lf_Dn ( izz, :)

            ELSE IF ( bcarr_x ( 2, iyy0, izz) == 3) THEN

              U_Rt_Dn ( izz, :) = U_Lf_Dn ( izz, :)
              U_Rt_Dn ( izz, 2) = - U_Rt_Dn ( izz, 2)
              U_Rt_Dn ( izz, 7) = - U_Rt_Dn ( izz, 7)
              U_Rt_Dn ( izz, 8) = - U_Rt_Dn ( izz, 8)

            END IF

! ----------

            IF ( bcarr_x ( 2, iyyp1, izz) == 2) THEN

              U_Rt_Up ( izz, :) = U_Lf_Up ( izz, :)

            ELSE IF ( bcarr_x ( 2, iyyp1, izz) == 3) THEN

              U_Rt_Up ( izz, :) = U_Lf_Up ( izz, :)
              U_Rt_Up ( izz, 2) = - U_Rt_Up ( izz, 2)
              U_Rt_Up ( izz, 7) = - U_Rt_Up ( izz, 7)
              U_Rt_Up ( izz, 8) = - U_Rt_Up ( izz, 8)

            END IF

! ----------

          END IF ! End of top x-face.

! ------------------------------
! We are at the bottom x-face. States that are obtained from zone "ixx0" should
! be reset with values from zone "ixxp1".

          IF ( ixx == x_indx_limits ( 1, iyy0, izz) - 1 ) THEN

! ----------

            IF ( bcarr_x ( 1, iyy0, izz) == 2) THEN

              U_Lf_Dn ( izz, :) = U_Rt_Dn ( izz, :)

            ELSE IF ( bcarr_x ( 1, iyy0, izz) == 3) THEN

              U_Lf_Dn ( izz, :) = U_Rt_Dn ( izz, :)
              U_Lf_Dn ( izz, 2) = - U_Lf_Dn ( izz, 2)
              U_Lf_Dn ( izz, 7) = - U_Lf_Dn ( izz, 7)
              U_Lf_Dn ( izz, 8) = - U_Lf_Dn ( izz, 8)

            END IF

! ----------

            IF ( bcarr_x ( 1, iyyp1, izz) == 2) THEN

              U_Lf_Up ( izz, :) = U_Rt_Up ( izz, :)

            ELSE IF ( bcarr_x ( 1, iyyp1, izz) == 3) THEN

              U_Lf_Up ( izz, :) = U_Rt_Up ( izz, :)
              U_Lf_Up ( izz, 2) = - U_Lf_Up ( izz, 2)
              U_Lf_Up ( izz, 7) = - U_Lf_Up ( izz, 7)
              U_Lf_Up ( izz, 8) = - U_Lf_Up ( izz, 8)

            END IF

! ----------

          END IF ! End of bottom x-face.

! ------------------------------
! We are at the top y-face. States that are obtained from zone "iyyp1" should
! be reset with values from zone "iyy0".

          IF ( iyy == y_indx_limits ( ixx0, 2, izz) ) THEN

! ----------

            IF ( bcarr_y ( ixx0, 2, izz) == 2) THEN

              U_Lf_Up ( izz, :) = U_Lf_Dn ( izz, :)

            ELSE IF ( bcarr_y ( ixx0, 2, izz) == 3) THEN

              U_Lf_Up ( izz, :) = U_Lf_Dn ( izz, :)
              U_Lf_Up ( izz, 3) = - U_Lf_Up ( izz, 3)
              U_Lf_Up ( izz, 6) = - U_Lf_Up ( izz, 6)
              U_Lf_Up ( izz, 8) = - U_Lf_Up ( izz, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixxp1, 2, izz) == 2) THEN

              U_Rt_Up ( izz, :) = U_Rt_Dn  ( izz, :)

            ELSE IF ( bcarr_y ( ixxp1, 2, izz) == 3) THEN

              U_Rt_Up ( izz, :) = U_Rt_Dn  ( izz, :)
              U_Rt_Up ( izz, 3) = - U_Rt_Up ( izz, 3)
              U_Rt_Up ( izz, 6) = - U_Rt_Up ( izz, 6)
              U_Rt_Up ( izz, 8) = - U_Rt_Up ( izz, 8)

            END IF

! ----------

          END IF ! End of top y-face.

! ------------------------------
! We are at the bottom y-face. States that are obtained from zone "iyy0" should
! be reset with values from zone "iyyp1".

          IF ( iyy == y_indx_limits ( ixx0, 1, izz) - 1 ) THEN

! ----------

            IF ( bcarr_y ( ixx0, 1, izz) == 2) THEN

              U_Lf_Dn ( izz, :) = U_Lf_Up ( izz, :)

            ELSE IF ( bcarr_y ( ixx0, 1, izz) == 3) THEN

              U_Lf_Dn ( izz, :) = U_Lf_Up ( izz, :)
              U_Lf_Dn ( izz, 3) = - U_Lf_Dn ( izz, 3)
              U_Lf_Dn ( izz, 6) = - U_Lf_Dn ( izz, 6)
              U_Lf_Dn ( izz, 8) = - U_Lf_Dn ( izz, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixxp1, 1, izz) == 2) THEN

              U_Rt_Dn  ( izz, :) = U_Rt_Up ( izz, :)

            ELSE IF ( bcarr_y ( ixxp1, 1, izz) == 3) THEN

              U_Rt_Dn  ( izz, :) = U_Rt_Up ( izz, :)
              U_Rt_Dn  ( izz, 3) = - U_Rt_Dn  ( izz, 3)
              U_Rt_Dn  ( izz, 6) = - U_Rt_Dn  ( izz, 6)
              U_Rt_Dn  ( izz, 8) = - U_Rt_Dn  ( izz, 8)

            END IF

! ----------

          END IF ! End of bottom y-face.

! ------------------------------

          END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------
! End of imposing boundary conditions.
! ----------------------------------------------------------------------
! Begin the steps for calling the Riemann solver.
! ----------------------------------------------------------------------

          DO izz = nriembegin, nriemend

! ------------------------------

          my_blend_HLL = blend_HLL ( izz)
          my_blend_MDHLL = blend_mdhll_to_mdhllc

          U_vrtx ( 1, :) = U_Rt_Up ( izz, :)
          U_vrtx ( 2, :) = U_Lf_Up ( izz, :)
          U_vrtx ( 3, :) = U_Lf_Dn ( izz, :)
          U_vrtx ( 4, :) = U_Rt_Dn ( izz, :)

          F_vrtx ( 1, :) = 0.0
          F_vrtx ( 2, :) = 0.0
          F_vrtx ( 3, :) = 0.0
          F_vrtx ( 4, :) = 0.0

          G_vrtx ( 1, :) = 0.0
          G_vrtx ( 2, :) = 0.0
          G_vrtx ( 3, :) = 0.0
          G_vrtx ( 4, :) = 0.0

          entropy_vrtx ( 1) = entropy_Rt_Up ( izz)
          entropy_vrtx ( 2) = entropy_Lf_Up ( izz)
          entropy_vrtx ( 3) = entropy_Lf_Dn ( izz)
          entropy_vrtx ( 4) = entropy_Rt_Dn ( izz)

          gamma_eos ( 1) = gamma_eos_Rt_Up ( izz)
          gamma_eos ( 2) = gamma_eos_Lf_Up ( izz)
          gamma_eos ( 3) = gamma_eos_Lf_Dn ( izz)
          gamma_eos ( 4) = gamma_eos_Rt_Dn ( izz)

          gamma_sdsp ( 1) = gamma_sdsp_Rt_Up ( izz)
          gamma_sdsp ( 2) = gamma_sdsp_Lf_Up ( izz)
          gamma_sdsp ( 3) = gamma_sdsp_Lf_Dn ( izz)
          gamma_sdsp ( 4) = gamma_sdsp_Rt_Dn ( izz)

! ------------------------------

          CALL MuSIC_MHD_2DRS_1_STRUCT

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, prsfloor, rhofloor, pi,

     1     my_blend_HLL, my_blend_MDHLL,
     1     my_Prs_Tot_Star, my_Prs_Gas_Star,

     1     Vlx_ALE, Vly_ALE,

     1     U_vrtx, F_vrtx, G_vrtx, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy_vrtx,

     1     U_Star, F_Star, G_Star,
     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)

! Now "F_Star, F_HLLC_Global_x" hold x-flux and "G_Star, G_HLLC_Global_y"
! hold y-flux.

! ------------------------------

          efzgr ( ixx, iyy, izz) = 0.5 * ( G_Star ( 6) - F_Star ( 7) )

! ------------------------------
! Add x-flux to x-face that lies in the lower y-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffy_expand <= iyy) .AND.
     1         ( iyy <= iy1 + ioffy_expand) ) THEN

            flux_x ( ixx, iyy, izz, :) = flux_x ( ixx, iyy, izz, :)
     1                                 + ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_x ( ixx, iyy, izz, 2) =
     1      flux_x ( ixx, iyy, izz, 2) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add x-flux to x-face that lies in the upper y-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffy_expand <= iyy + 1) .AND.
     1         ( iyy + 1 <= iy1 + ioffy_expand) ) THEN

            flux_temp ( ixx, iyy + 1, izz, :)=
     1                                   ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1)
     1      flux_temp ( ixx, iyy + 1, izz, 2) =
     1      flux_temp ( ixx, iyy + 1, izz, 2)
     1    - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add y-flux to y-face that lies in the lower x-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffx_expand <= ixx) .AND.
     1         ( ixx <= ix1 + ioffx_expand) ) THEN

            flux_y ( ixx, iyy, izz, :) = flux_y ( ixx, iyy, izz, :)
     1                                 + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_y ( ixx, iyy, izz, 3) =
     1      flux_y ( ixx, iyy, izz, 3) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
! Add y-flux to y-face that lies in the upper x-direction.
! If the face is reflecting, just use 1D flux.

          IF ( ( 1 - ioffx_expand <= ixx + 1) .AND.
     1         ( ixx + 1 <= ix1 + ioffx_expand) ) THEN

            flux_y ( ixx + 1, iyy, izz, :) =
     1                                flux_y ( ixx + 1, iyy, izz, :)
     1                                 + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_y ( ixx + 1, iyy, izz, 3) =
     1    flux_y ( ixx + 1, iyy, izz, 3) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------

          END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------
! End of steps for calling the Riemann solver.
! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    flux_x, flux_temp)                                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)
      DO izz = 0 - ioffz, iz1 + ioffz
        DO iyy = 0 - ioffy, iy1 + ioffy
          DO ixx = 0 - ioffx, ix1 + ioffx

            flux_x ( ixx, iyy, izz, :) = flux_x ( ixx, iyy, izz, :)
     1                              + flux_temp ( ixx, iyy, izz, :)

          END DO
        END DO
      END DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE APPLY_MULTID_RS_Z_EDGE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_FLUX_X ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, entropygr, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_x, prs_flux_x, entropy_flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes"
! and adds one-dimensional contribution to "flux_x and prs_flux_x".

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_x".

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1      bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1     prs_flux_x, entropy_flux_x

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, ixxp1, nriembegin, nriemend, one_x,
     1        nspecies_min, nspecies_max

      REAL ob_12, ob_8pi, a, b, ob_6_oned

! ----------------------------------------------------------------------

! These are the variables that go through the 1D Riemann solver interface.

      INTEGER use_approx_eos, split_prsflux_re, iHLLEM_choice

      REAL rho_L, rho_R, prs_L, prs_R,
     1     entropy_L, entropy_R,

     1     vlx_L, vlx_R, vly_L, vly_R,
     1     vlz_L, vlz_R,

     1     bfxre, bfy_L, bfy_R, bfz_L, bfz_R,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     msonic, divvel, blend_HLL, S_L, S_R, entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: flux_re,
     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

! ----------------------------------------------------------------------

! The computation is not dynamically active in the x-direction
! ( as indicated by setting "ioffx = 0" on the outside) so return.

      IF ( ioffx .EQ. 0) RETURN

      ob_6_oned = 1.0 - 4.0 * blend_multid_flux

#if ( TRUE_2D == 1)
      ob_6_oned = 1.0 - 2.0 * blend_multid_flux
#endif

      iHLLEM_choice = 1  ! 1 or 2 gives HLLEM; 0 gives HLL.

! ----------

      ob_12 = 1.0 / 12.0
      ob_8pi = 1.0 / ( 8.0 * pi)

      one_x = MIN0 ( expand_by_one_zone, ioffx)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

! ----------

      use_approx_eos = 1 ! Just use simple EOS with constant polytropic index.

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ixxp1, a, b,                 &
!$OMP&    nriembegin, nriemend,                                         &

!$OMP&    rho_L, prs_L, vlx_L, vly_L, vlz_L, bfy_L, bfz_L,              &
!$OMP&    mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,               &
!$OMP&    specie_L, entropy_L,                                          &

!$OMP&    bfxre, msonic, divvel, blend_HLL,                             &

!$OMP&    rho_R, prs_R, vlx_R, vly_R, vlz_R, bfy_R, bfz_R,              &
!$OMP&    mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,               &
!$OMP&    specie_R, entropy_R,                                          &

!$OMP&    S_L, S_R, F_L, F_R, U_HLL, F_HLL,                             &
!$OMP&    U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM,           &
!$OMP&    flux_re, entropy_flux_re)                                     &

!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, ob_12, ob_8pi, one_x, expand_by_one_zone,              &
!$OMP&    dtcur, dt_shift, ob_6_oned, iHLLEM_choice,                    &
!$OMP&    nspecies_min, nspecies_max,                                   &

!$OMP&    gamma, smallnum, rhofloor, prsfloor, pi, fracdiff,            &
!$OMP&    use_approx_eos, split_prsflux_re,                             &
!$OMP&    lapidus_threshold, lapidus_coef, flatten_coef,                &
!$OMP&    gasconst, cp_specie, molewt_specie, heat_of_formation,        &

!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &

!$OMP&    x_indx_limits, bcarr_x, u, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    flux_x, prs_flux_x, entropygr, entropy_flux_x,                &
!$OMP&    mean_mol_wt_gr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr,  &
!$OMP&    divvelgr)                                                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - MIN0 ( expand_by_one_zone, ioffz),
     1       iz1 + MIN0 ( expand_by_one_zone, ioffz)

        DO iyy = 1 - MIN0 ( expand_by_one_zone, ioffy),
     1         iy1 + MIN0 ( expand_by_one_zone, ioffy)

! ----------------------------------------------------------------------

! Load in the left and right states for the Riemann solver.

          nriembegin = x_indx_limits ( 1, iyy, izz) - 1 - one_x
          nriemend = x_indx_limits ( 2, iyy, izz) + one_x

          DO ixx = nriembegin, nriemend

          ixxp1 = ixx + 1

! ----------------------------------------------------------------------

          rho_L =   u ( ixx, iyy, izz, 1, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 1, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 1, 5)

          vlx_L =   u ( ixx, iyy, izz, 2, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 2, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 2, 5)

          vly_L =   u ( ixx, iyy, izz, 3, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 3, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 3, 5)

          vlz_L =   u ( ixx, iyy, izz, 4, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 4, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 4, 5)

          prs_L =   u ( ixx, iyy, izz, 5, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 5, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 5, 5)

          bfy_L =   u ( ixx, iyy, izz, 7, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 7, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 7, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          bfz_L =   u ( ixx, iyy, izz, 8, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 8, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 8, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

          mean_mol_wt_L = mean_mol_wt_gr ( ixx, iyy, izz)
          gamma_eos_L = gamma_eos_gr ( ixx, iyy, izz)
          gamma_soundspeed_L = gamma_soundspeed_gr ( ixx, iyy, izz)

! ----------

          a =       u ( ixx, iyy, izz, 6, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 6, 2)
     1 + dt_shift * u ( ixx, iyy, izz, 6, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 + 2.0 * ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          b =       u ( ixxp1, iyy, izz, 6, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 6, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 6, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 + 2.0 * ob_12 * bfxccmodes ( ixxp1, iyy, izz, 5)
#endif

          bfxre = 0.5 * ( a + b)

! ----------

          rho_R =   u ( ixxp1, iyy, izz, 1, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 1, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 1, 5)

          vlx_R =   u ( ixxp1, iyy, izz, 2, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 2, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 2, 5)

          vly_R =   u ( ixxp1, iyy, izz, 3, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 3, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 3, 5)

          vlz_R =   u ( ixxp1, iyy, izz, 4, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 4, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 4, 5)

          prs_R =   u ( ixxp1, iyy, izz, 5, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 5, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 5, 5)

          bfy_R =   u ( ixxp1, iyy, izz, 7, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 7, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 7, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfyccmodes ( ixxp1, iyy, izz, 5)
#endif

          bfz_R =   u ( ixxp1, iyy, izz, 8, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, 8, 2)
     1 + dt_shift * u ( ixxp1, iyy, izz, 8, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfzccmodes ( ixxp1, iyy, izz, 5)
#endif

          mean_mol_wt_R = mean_mol_wt_gr ( ixxp1, iyy, izz)
          gamma_eos_R = gamma_eos_gr ( ixxp1, iyy, izz)
          gamma_soundspeed_R = gamma_soundspeed_gr ( ixxp1, iyy, izz)

! ----------

          rho_L = AMAX1 ( rhofloor, rho_L)
          a = 1.0 / rho_L

          vlx_L = vlx_L * a
          vly_L = vly_L * a
          vlz_L = vlz_L * a

          b = 0.5 * rho_L * ( vlx_L**2 + vly_L**2 + vlz_L**2)
     1      + ( bfxre**2 + bfy_L**2 + bfz_L**2) * ob_8pi

          prs_L = ( prs_L - b) * ( gamma_eos_L - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_L = entropygr ( ixx, iyy, izz, 1)
     1          + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1     + dt_shift * entropygr ( ixx, iyy, izz, 5)
          IF ( prs_L .LE. prsfloor) prs_L = entropy_L
     1                              * rho_L**( gamma_eos_L - 1.0)
#endif
          prs_L = AMAX1 ( prsfloor, prs_L)

#if ( NFLUID >= 1)
          specie_L ( 1: NFLUID) =
     1             u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1     + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 2)
     1+ dt_shift * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 5)

          specie_L ( 1: NFLUID) = specie_L ( 1: NFLUID) * a
#endif

! ----------

          rho_R = AMAX1 ( rhofloor, rho_R)
          a = 1.0 / rho_R

          vlx_R = vlx_R * a
          vly_R = vly_R * a
          vlz_R = vlz_R * a

          b = 0.5 * rho_R * ( vlx_R**2 + vly_R**2 + vlz_R**2)
     1      + ( bfxre**2 + bfy_R**2 + bfz_R**2) * ob_8pi

          prs_R = ( prs_R - b) * ( gamma_eos_R - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_R = entropygr ( ixxp1, iyy, izz, 1)
     1          - 0.5 * entropygr ( ixxp1, iyy, izz, 2)
     1     + dt_shift * entropygr ( ixxp1, iyy, izz, 5)
          IF ( prs_R .LE. prsfloor) prs_R = entropy_R
     1                              * rho_R**( gamma_eos_R - 1.0)
#endif
          prs_R = AMAX1 ( prsfloor, prs_R)

#if ( NFLUID >= 1)
          specie_R ( 1: NFLUID) =
     1             u ( ixxp1, iyy, izz, nspecies_min: nspecies_max, 1)
     1     - 0.5 * u ( ixxp1, iyy, izz, nspecies_min: nspecies_max, 2)
     1+ dt_shift * u ( ixxp1, iyy, izz, nspecies_min: nspecies_max, 5)

          specie_R ( 1: NFLUID) = specie_R ( 1: NFLUID) * a
#endif

! ----------

          msonic = AMIN1 ( msonicgr ( ixx, iyy, izz),
     1                     msonicgr ( ixxp1, iyy, izz))
          divvel = AMAX1 ( ABS ( divvelgr ( ixx, iyy, izz)),
     1                     ABS ( divvelgr ( ixxp1, iyy, izz)) )

          a = divvel / ( flatten_coef * msonic) - 1.0

          IF ( a < 0.0 ) THEN
            blend_HLL = 0.0
          ELSE
            blend_HLL = AMIN1 ( 1.0, a)
          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Continuitive boundary conditions in the x-direction shown below.

          IF ( ( ixx == nriembegin + one_x) .AND.
     1         ( bcarr_x ( 1, iyy, izz) == 2) ) THEN

            rho_L = rho_R
            prs_L = prs_R

            vlx_L = vlx_R
            vly_L = vly_R
            vlz_L = vlz_R

            bfy_L = bfy_R
            bfz_L = bfz_R

#if ( NFLUID >= 1)
            specie_L ( 1: NFLUID) = specie_R ( 1: NFLUID)
#endif

          END IF

          IF ( ( ixx == nriemend - one_x) .AND.
     1         ( bcarr_x ( 2, iyy, izz) == 2) ) THEN

            rho_R = rho_L
            prs_R = prs_L

            vlx_R = vlx_L
            vly_R = vly_L
            vlz_R = vlz_L

            bfy_R = bfy_L
            bfz_R = bfz_L

#if ( NFLUID >= 1)
            specie_R ( 1: NFLUID) = specie_L ( 1: NFLUID)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Reflection boundary conditions in the x-direction shown below.

          IF ( ( ixx == nriembegin + one_x) .AND.
     1         ( bcarr_x ( 1, iyy, izz) == 3) ) THEN

            rho_L = rho_R
            prs_L = prs_R

            vlx_L= -vlx_R
            vly_L = vly_R
            vlz_L = vlz_R

            bfy_L= -bfy_R
            bfz_L= -bfz_R

#if ( NFLUID >= 1)
            specie_L ( 1: NFLUID) = specie_R ( 1: NFLUID)
#endif

          END IF

          IF ( ( ixx == nriemend - one_x) .AND.
     1         ( bcarr_x ( 2, iyy, izz) == 3) ) THEN

            rho_R = rho_L
            prs_R = prs_L

            vlx_R= -vlx_L
            vly_R = vly_L
            vlz_R = vlz_L

            bfy_R= -bfy_L
            bfz_R= -bfz_L

#if ( NFLUID >= 1)
            specie_R ( 1: NFLUID) = specie_L ( 1: NFLUID)
#endif

          END IF

! ----------------------------------------------------------------------

! Call the Riemann solver. HLLEM is the default.

      CALL MHD_RIEM_HLLEM_PTWISE

     1     ( gamma, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vlx_L, vlx_R,
     1        vly_L, vly_R,
     1        vlz_L, vlz_R,
     1        bfxre,
     1        bfy_L, bfy_R,
     1        bfz_L, bfz_R,
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

! Fluxes have been evaluated at face centers. Now hand them back.

          IF ( ( ixx == nriembegin + one_x) .AND.
     1         ( bcarr_x ( 1, iyy, izz) == 3) ) THEN

! ----------

            flux_x ( ixx, iyy, izz, 1: 5) = flux_re ( 1: 5)
            flux_x ( ixx, iyy, izz, 6) = 0.0
            flux_x ( ixx, iyy, izz, 7) = flux_re ( 6)
            flux_x ( ixx, iyy, izz, 8) = flux_re ( 7)

            prs_flux_x ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_x ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1                    flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_x ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          ELSE IF ( ( ixx == nriemend - one_x) .AND.
     1              ( bcarr_x ( 2, iyy, izz) == 3) ) THEN

! ----------

            flux_x ( ixx, iyy, izz, 1: 5) = flux_re ( 1: 5)
            flux_x ( ixx, iyy, izz, 6) = 0.0
            flux_x ( ixx, iyy, izz, 7) = flux_re ( 6)
            flux_x ( ixx, iyy, izz, 8) = flux_re ( 7)

            prs_flux_x ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_x ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1                    flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_x ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          ELSE

! ----------

            flux_x ( ixx, iyy, izz, 1: 5) =
     1                                flux_x ( ixx, iyy, izz, 1: 5)
     1                                + ob_6_oned * flux_re ( 1: 5)
            flux_x ( ixx, iyy, izz, 6) = 0.0
            flux_x ( ixx, iyy, izz, 7) = flux_x ( ixx, iyy, izz, 7)
     1                                   + ob_6_oned * flux_re ( 6)
            flux_x ( ixx, iyy, izz, 8) = flux_x ( ixx, iyy, izz, 8)
     1                                   + ob_6_oned * flux_re ( 7)

            prs_flux_x ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_x ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1      flux_x ( ixx, iyy, izz, nspecies_min: nspecies_max)
     1      + ob_6_oned * flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_x ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          END IF

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

! ----------------------------------------------------------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_FLUX_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_FLUX_Y ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, entropygr, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_y, prs_flux_y, entropy_flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes"
! and adds one-dimensional contribution to "flux_y and prs_flux_y".

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_y".

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1      bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_y

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1     prs_flux_y, entropy_flux_y

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  2,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, iyyp1, nriembegin, nriemend, one_y,
     1        nspecies_min, nspecies_max

      REAL ob_12, ob_8pi, a, b, ob_6_oned

! ----------------------------------------------------------------------

! These are the variables that go through the 1D Riemann solver interface.

      INTEGER use_approx_eos, split_prsflux_re, iHLLEM_choice

      REAL rho_L, rho_R, prs_L, prs_R,
     1     entropy_L, entropy_R,

     1     vlx_L, vlx_R, vly_L, vly_R,
     1     vlz_L, vlz_R,

     1     bfx_L, bfx_R, bfyre, bfz_L, bfz_R,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     msonic, divvel, blend_HLL, S_L, S_R, entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: flux_re,
     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

! ----------------------------------------------------------------------

! The computation is not dynamically active in the y-direction
! ( as indicated by setting "ioffy = 0" on the outside) so return.

      IF ( ioffy .EQ. 0) RETURN

      ob_6_oned = 1.0 - 4.0 * blend_multid_flux

#if ( TRUE_2D == 1)
      ob_6_oned = 1.0 - 2.0 * blend_multid_flux
#endif

      iHLLEM_choice = 1  ! 1 or 2 gives HLLEM; 0 gives HLL.

! ----------

      ob_12 = 1.0 / 12.0
      ob_8pi = 1.0 / ( 8.0 * pi)

      one_y = MIN0 ( expand_by_one_zone, ioffy)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

! ----------

      use_approx_eos = 1 ! Just use simple EOS with constant polytropic index.

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, iyyp1, a, b,                 &
!$OMP&    nriembegin, nriemend,                                         &

!$OMP&    rho_L, prs_L, vlx_L, vly_L, vlz_L, bfx_L, bfz_L,              &
!$OMP&    mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,               &
!$OMP&    specie_L, entropy_L,                                          &

!$OMP&    bfyre, msonic, divvel, blend_HLL,                             &

!$OMP&    rho_R, prs_R, vlx_R, vly_R, vlz_R, bfx_R, bfz_R,              &
!$OMP&    mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,               &
!$OMP&    specie_R, entropy_R,                                          &

!$OMP&    S_L, S_R, F_L, F_R, U_HLL, F_HLL,                             &
!$OMP&    U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM,           &
!$OMP&    flux_re, entropy_flux_re)                                     &

!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, ob_12, ob_8pi, one_y, expand_by_one_zone,              &
!$OMP&    dtcur, dt_shift, ob_6_oned, iHLLEM_choice,                    &
!$OMP&    nspecies_min, nspecies_max,                                   &

!$OMP&    gamma, smallnum, rhofloor, prsfloor, pi, fracdiff,            &
!$OMP&    use_approx_eos, split_prsflux_re,                             &
!$OMP&    lapidus_threshold, lapidus_coef, flatten_coef,                &
!$OMP&    gasconst, cp_specie, molewt_specie, heat_of_formation,        &

!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &

!$OMP&    y_indx_limits, bcarr_y, u, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    flux_y, prs_flux_y, entropygr, entropy_flux_y,                &
!$OMP&    mean_mol_wt_gr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr,  &
!$OMP&    divvelgr)                                                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - MIN0 ( expand_by_one_zone, ioffz),
     1       iz1 + MIN0 ( expand_by_one_zone, ioffz)

        DO ixx = 1 - MIN0 ( expand_by_one_zone, ioffx),
     1         ix1 + MIN0 ( expand_by_one_zone, ioffx)

! ----------------------------------------------------------------------

! Load in the left and right states for the Riemann solver.

          nriembegin = y_indx_limits ( ixx, 1, izz) - 1 - one_y
          nriemend = y_indx_limits ( ixx, 2, izz) + one_y

          DO iyy = nriembegin, nriemend

          iyyp1 = iyy + 1

! ----------------------------------------------------------------------

          rho_L =   u ( ixx, iyy, izz, 1, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 1, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 1, 5)

          vlx_L =   u ( ixx, iyy, izz, 2, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 2, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 2, 5)

          vly_L =   u ( ixx, iyy, izz, 3, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 3, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 3, 5)

          vlz_L =   u ( ixx, iyy, izz, 4, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 4, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 4, 5)

          prs_L =   u ( ixx, iyy, izz, 5, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 5, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 5, 5)

          bfx_L =   u ( ixx, iyy, izz, 6, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 6, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 6, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          bfz_L =   u ( ixx, iyy, izz, 8, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 8, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 8, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

          mean_mol_wt_L = mean_mol_wt_gr ( ixx, iyy, izz)
          gamma_eos_L = gamma_eos_gr ( ixx, iyy, izz)
          gamma_soundspeed_L = gamma_soundspeed_gr ( ixx, iyy, izz) 

! ----------

          a =       u ( ixx, iyy, izz, 7, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 7, 3)
     1 + dt_shift * u ( ixx, iyy, izz, 7, 5) 
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 + 2.0 * ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          b =       u ( ixx, iyyp1, izz, 7, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 7, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 7, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 + 2.0 * ob_12 * bfyccmodes ( ixx, iyyp1, izz, 5)
#endif

          bfyre = 0.5 * ( a + b)

! ----------

          rho_R =   u ( ixx, iyyp1, izz, 1, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 1, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 1, 5)

          vlx_R =   u ( ixx, iyyp1, izz, 2, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 2, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 2, 5)

          vly_R =   u ( ixx, iyyp1, izz, 3, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 3, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 3, 5)

          vlz_R =   u ( ixx, iyyp1, izz, 4, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 4, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 4, 5)

          prs_R =   u ( ixx, iyyp1, izz, 5, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 5, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 5, 5)

          bfx_R =   u ( ixx, iyyp1, izz, 6, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 6, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 6, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 - ob_12 * bfxccmodes ( ixx, iyyp1, izz, 5)
#endif

          bfz_R =   u ( ixx, iyyp1, izz, 8, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, 8, 3)
     1 + dt_shift * u ( ixx, iyyp1, izz, 8, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 - ob_12 * bfzccmodes ( ixx, iyyp1, izz, 5)
#endif

          mean_mol_wt_R = mean_mol_wt_gr ( ixx, iyyp1, izz)
          gamma_eos_R = gamma_eos_gr ( ixx, iyyp1, izz)
          gamma_soundspeed_R = gamma_soundspeed_gr ( ixx, iyyp1, izz)

! ----------

          rho_L = AMAX1 ( rhofloor, rho_L)
          a = 1.0 / rho_L

          vlx_L = vlx_L * a
          vly_L = vly_L * a
          vlz_L = vlz_L * a

          b = 0.5 * rho_L * ( vlx_L**2 + vly_L**2 + vlz_L**2)
     1      + ( bfx_L**2 + bfyre**2 + bfz_L**2) * ob_8pi

          prs_L = ( prs_L - b) * ( gamma_eos_L - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_L = entropygr ( ixx, iyy, izz, 1)
     1          + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1     + dt_shift * entropygr ( ixx, iyy, izz, 5)
          IF ( prs_L .LE. prsfloor) prs_L = entropy_L
     1                              * rho_L**( gamma_eos_L - 1.0)
#endif
          prs_L = AMAX1 ( prsfloor, prs_L)

#if ( NFLUID >= 1)
          specie_L ( 1: NFLUID) =
     1             u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1     + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 3)
     1+ dt_shift * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 5)

          specie_L ( 1: NFLUID) = specie_L ( 1: NFLUID) * a
#endif

! ----------

          rho_R = AMAX1 ( rhofloor, rho_R)
          a = 1.0 / rho_R

          vlx_R = vlx_R * a
          vly_R = vly_R * a
          vlz_R = vlz_R * a

          b = 0.5 * rho_R * ( vlx_R**2 + vly_R**2 + vlz_R**2)
     1      + ( bfx_R**2 + bfyre**2 + bfz_R**2) * ob_8pi

          prs_R = ( prs_R - b) * ( gamma_eos_R - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_R = entropygr ( ixx, iyyp1, izz, 1)
     1          - 0.5 * entropygr ( ixx, iyyp1, izz, 3)
     1     + dt_shift * entropygr ( ixx, iyyp1, izz, 5)
          IF ( prs_R .LE. prsfloor) prs_R = entropy_R
     1                              * rho_R**( gamma_eos_R - 1.0)
#endif
          prs_R = AMAX1 ( prsfloor, prs_R)

#if ( NFLUID >= 1)
          specie_R ( 1: NFLUID) =
     1             u ( ixx, iyyp1, izz, nspecies_min: nspecies_max, 1)
     1     - 0.5 * u ( ixx, iyyp1, izz, nspecies_min: nspecies_max, 3)
     1+ dt_shift * u ( ixx, iyyp1, izz, nspecies_min: nspecies_max, 5)

          specie_R ( 1: NFLUID) = specie_R ( 1: NFLUID) * a
#endif

! ----------

          msonic = AMIN1 ( msonicgr ( ixx, iyy, izz),
     1                     msonicgr ( ixx, iyyp1, izz))
          divvel = AMAX1 ( ABS ( divvelgr ( ixx, iyy, izz)),
     1                     ABS ( divvelgr ( ixx, iyyp1, izz)) )

          a = divvel / ( flatten_coef * msonic) - 1.0

          IF ( a < 0.0 ) THEN
            blend_HLL = 0.0
          ELSE
            blend_HLL = AMIN1 ( 1.0, a)
          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Continuitive boundary conditions in the y-direction shown below.

          IF ( ( iyy == nriembegin + one_y) .AND.
     1         ( bcarr_y ( ixx, 1, izz) == 2) ) THEN

            rho_L = rho_R
            prs_L = prs_R

            vlx_L = vlx_R
            vly_L = vly_R
            vlz_L = vlz_R

            bfx_L = bfx_R
            bfz_L = bfz_R

#if ( NFLUID >= 1)
            specie_L ( 1: NFLUID) = specie_R ( 1: NFLUID)
#endif

          END IF

          IF ( ( iyy == nriemend - one_y) .AND.
     1         ( bcarr_y ( ixx, 2, izz) == 2) ) THEN

            rho_R = rho_L
            prs_R = prs_L

            vlx_R = vlx_L
            vly_R = vly_L
            vlz_R = vlz_L

            bfx_R = bfx_L
            bfz_R = bfz_L

#if ( NFLUID >= 1)
            specie_R ( 1: NFLUID) = specie_L ( 1: NFLUID)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Reflection boundary conditions in the y-direction shown below.

          IF ( ( iyy == nriembegin + one_y) .AND.
     1         ( bcarr_y ( ixx, 1, izz) == 3) ) THEN

            rho_L = rho_R
            prs_L = prs_R

            vlx_L = vlx_R
            vly_L= -vly_R
            vlz_L = vlz_R

            bfx_L= -bfx_R
            bfz_L= -bfz_R

#if ( NFLUID >= 1)
            specie_L ( 1: NFLUID) = specie_R ( 1: NFLUID)
#endif

          END IF

          IF ( ( iyy == nriemend - one_y) .AND.
     1         ( bcarr_y ( ixx, 2, izz) == 3) ) THEN

            rho_R = rho_L
            prs_R = prs_L

            vlx_R = vlx_L
            vly_R= -vly_L
            vlz_R = vlz_L

            bfx_R= -bfx_L
            bfz_R= -bfz_L

#if ( NFLUID >= 1)
            specie_R ( 1: NFLUID) = specie_L ( 1: NFLUID)
#endif

          END IF

! ----------------------------------------------------------------------

! Call the Riemann solver with cyclic rotation of input variables.
! HLLEM is the default.

      CALL MHD_RIEM_HLLEM_PTWISE

     1     ( gamma, smallnum, prsfloor, rhofloor, pi,
     1     use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vly_L, vly_R,
     1        vlz_L, vlz_R,
     1        vlx_L, vlx_R,
     1        bfyre,
     1        bfz_L, bfz_R,
     1        bfx_L, bfx_R,
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

! Fluxes have been evaluated at face centers.

          IF ( ( iyy == nriembegin + one_y) .AND.
     1         ( bcarr_y ( ixx, 1, izz) == 3) ) THEN

! ----------

            flux_y ( ixx, iyy, izz, 1) = flux_re ( 1)
            flux_y ( ixx, iyy, izz, 2) = flux_re ( 4)
            flux_y ( ixx, iyy, izz, 3) = flux_re ( 2)
            flux_y ( ixx, iyy, izz, 4) = flux_re ( 3)
            flux_y ( ixx, iyy, izz, 5) = flux_re ( 5)

            flux_y ( ixx, iyy, izz, 6) = flux_re ( 7)
            flux_y ( ixx, iyy, izz, 7) = 0.0
            flux_y ( ixx, iyy, izz, 8) = flux_re ( 6)

            prs_flux_y ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_y ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1                    flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_y ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          ELSE IF ( ( iyy == nriemend - one_y) .AND.
     1              ( bcarr_y ( ixx, 2, izz) == 3) ) THEN

! ----------

            flux_y ( ixx, iyy, izz, 1) = flux_re ( 1)
            flux_y ( ixx, iyy, izz, 2) = flux_re ( 4)
            flux_y ( ixx, iyy, izz, 3) = flux_re ( 2)
            flux_y ( ixx, iyy, izz, 4) = flux_re ( 3)
            flux_y ( ixx, iyy, izz, 5) = flux_re ( 5)

            flux_y ( ixx, iyy, izz, 6) = flux_re ( 7)
            flux_y ( ixx, iyy, izz, 7) = 0.0
            flux_y ( ixx, iyy, izz, 8) = flux_re ( 6)

            prs_flux_y ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_y ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1                    flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_y ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          ELSE

! ----------

            flux_y ( ixx, iyy, izz, 1) = flux_y ( ixx, iyy, izz, 1)
     1                                   + ob_6_oned * flux_re ( 1)

            flux_y ( ixx, iyy, izz, 2) = flux_y ( ixx, iyy, izz, 2)
     1                                   + ob_6_oned * flux_re ( 4)
            flux_y ( ixx, iyy, izz, 3) = flux_y ( ixx, iyy, izz, 3)
     1                                   + ob_6_oned * flux_re ( 2)
            flux_y ( ixx, iyy, izz, 4) = flux_y ( ixx, iyy, izz, 4)
     1                                   + ob_6_oned * flux_re ( 3)

            flux_y ( ixx, iyy, izz, 5) = flux_y ( ixx, iyy, izz, 5)
     1                                   + ob_6_oned * flux_re ( 5)

            flux_y ( ixx, iyy, izz, 6) = flux_y ( ixx, iyy, izz, 6)
     1                                   + ob_6_oned * flux_re ( 7)
            flux_y ( ixx, iyy, izz, 7) = 0.0
            flux_y ( ixx, iyy, izz, 8) = flux_y ( ixx, iyy, izz, 8)
     1                                   + ob_6_oned * flux_re ( 6)

            prs_flux_y ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_y ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1      flux_y ( ixx, iyy, izz, nspecies_min: nspecies_max)
     1      + ob_6_oned * flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_y ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          END IF

! ----------------------------------------------------------------------

          END DO ! End of DO iyy = loop.

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_FLUX_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_FLUX_Z ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, entropygr, bfxbdy, bfybdy, bfzbdy,
     1        bfxccmodes, bfyccmodes, bfzccmodes,
     1        flux_z, prs_flux_z, entropy_flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy, bfxccmodes, bfyccmodes, bfzccmodes"
! and adds one-dimensional contribution to "flux_z and prs_flux_z".

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_z".

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes, n_cc_bfld_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_bfld_modes) ::
     1      bfxccmodes, bfyccmodes, bfzccmodes

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_z

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1      prs_flux_z, entropy_flux_z

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  2) ::
     1         z_indx_limits, bcarr_z

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, izzp1, nriembegin, nriemend, one_z,
     1        nspecies_min, nspecies_max

      REAL ob_12, ob_8pi, a, b, ob_6_oned

! ----------------------------------------------------------------------

! These are the variables that go through the 1D Riemann solver interface.

      INTEGER use_approx_eos, split_prsflux_re, iHLLEM_choice

      REAL rho_L, rho_R, prs_L, prs_R,
     1     entropy_L, entropy_R,

     1     vlx_L, vlx_R, vly_L, vly_R,
     1     vlz_L, vlz_R,

     1     bfx_L, bfx_R, bfy_L, bfy_R, bfzre,

     1     mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,
     1     mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,

     1     msonic, divvel, blend_HLL, S_L, S_R, entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID) :: specie_L, specie_R

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: flux_re,
     1     F_L, F_R, U_HLL, F_HLL,
     1     U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM

! ----------------------------------------------------------------------

! The computation is not dynamically active in the z-direction
! ( as indicated by setting "ioffz = 0" on the outside) so return.

      IF ( ioffz .EQ. 0) RETURN

      ob_6_oned = 1.0 - 4.0 * blend_multid_flux

#if ( TRUE_2D == 1)
      ob_6_oned = 1.0 - 2.0 * blend_multid_flux
#endif

      iHLLEM_choice = 1  ! 1 or 2 gives HLLEM; 0 gives HLL.

! ----------

      ob_12 = 1.0 / 12.0
      ob_8pi = 1.0 / ( 8.0 * pi)

      one_z = MIN0 ( expand_by_one_zone, ioffz)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

! ----------

      use_approx_eos = 1 ! Just use simple EOS with constant polytropic index.

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, izzp1, a, b,                 &
!$OMP&    nriembegin, nriemend,                                         &

!$OMP&    rho_L, prs_L, vlx_L, vly_L, vlz_L, bfx_L, bfy_L,              &
!$OMP&    mean_mol_wt_L, gamma_eos_L, gamma_soundspeed_L,               &
!$OMP&    specie_L, entropy_L,                                          &

!$OMP&    bfzre, msonic, divvel, blend_HLL,                             &

!$OMP&    rho_R, prs_R, vlx_R, vly_R, vlz_R, bfx_R, bfy_R,              &
!$OMP&    mean_mol_wt_R, gamma_eos_R, gamma_soundspeed_R,               &
!$OMP&    specie_R, entropy_R,                                          &

!$OMP&    S_L, S_R, F_L, F_R, U_HLL, F_HLL,                             &
!$OMP&    U_HLLEM, F_HLLEM, dU_HLLEM, dU_HLLEM_Min, dF_HLLEM,           &
!$OMP&    flux_re, entropy_flux_re)                                     &

!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, ob_12, ob_8pi, one_z, expand_by_one_zone,              &
!$OMP&    dtcur, dt_shift, ob_6_oned, iHLLEM_choice,                    &
!$OMP&    nspecies_min, nspecies_max,                                   &

!$OMP&    gamma, smallnum, rhofloor, prsfloor, pi, fracdiff,            &
!$OMP&    use_approx_eos, split_prsflux_re,                             &
!$OMP&    lapidus_threshold, lapidus_coef, flatten_coef,                &
!$OMP&    gasconst, cp_specie, molewt_specie, heat_of_formation,        &

!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &

!$OMP&    z_indx_limits, bcarr_z, u, bfxccmodes, bfyccmodes, bfzccmodes,&
!$OMP&    flux_z, prs_flux_z, entropygr, entropy_flux_z,                &
!$OMP&    mean_mol_wt_gr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr,  &
!$OMP&    divvelgr)                                                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO iyy = 1 - MIN0 ( expand_by_one_zone, ioffy),
     1       iy1 + MIN0 ( expand_by_one_zone, ioffy)

        DO ixx = 1 - MIN0 ( expand_by_one_zone, ioffx),
     1         ix1 + MIN0 ( expand_by_one_zone, ioffx)

! ----------------------------------------------------------------------

! Load in the left and right states for the Riemann solver.

          nriembegin = z_indx_limits ( ixx, iyy, 1) - 1 - one_z
          nriemend = z_indx_limits ( ixx, iyy, 2) + one_z

          DO izz = nriembegin, nriemend

          izzp1 = izz + 1

! ----------------------------------------------------------------------

          rho_L =   u ( ixx, iyy, izz, 1, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 1, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 1, 5)

          vlx_L =   u ( ixx, iyy, izz, 2, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 2, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 2, 5)

          vly_L =   u ( ixx, iyy, izz, 3, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 3, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 3, 5)

          vlz_L =   u ( ixx, iyy, izz, 4, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 4, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 4, 5)

          prs_L =   u ( ixx, iyy, izz, 5, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 5, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 5, 5)

          bfx_L =   u ( ixx, iyy, izz, 6, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 6, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 6, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfxccmodes ( ixx, iyy, izz, 5)
#endif

          bfy_L =   u ( ixx, iyy, izz, 7, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 7, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 7, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfyccmodes ( ixx, iyy, izz, 5)
#endif

          mean_mol_wt_L = mean_mol_wt_gr ( ixx, iyy, izz)
          gamma_eos_L = gamma_eos_gr ( ixx, iyy, izz)
          gamma_soundspeed_L = gamma_soundspeed_gr ( ixx, iyy, izz) 

! ----------

          a =       u ( ixx, iyy, izz, 8, 1)
     1      + 0.5 * u ( ixx, iyy, izz, 8, 4)
     1 + dt_shift * u ( ixx, iyy, izz, 8, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 + 2.0 * ob_12 * bfzccmodes ( ixx, iyy, izz, 5)
#endif

          b =       u ( ixx, iyy, izzp1, 8, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 8, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 8, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1 + 2.0 * ob_12 * bfzccmodes ( ixx, iyy, izzp1, 5)
#endif

          bfzre = 0.5 * ( a + b)

! ----------

          rho_R =   u ( ixx, iyy, izzp1, 1, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 1, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 1, 5)

          vlx_R =   u ( ixx, iyy, izzp1, 2, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 2, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 2, 5)

          vly_R =   u ( ixx, iyy, izzp1, 3, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 3, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 3, 5)

          vlz_R =   u ( ixx, iyy, izzp1, 4, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 4, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 4, 5)

          prs_R =   u ( ixx, iyy, izzp1, 5, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 5, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 5, 5)

          bfx_R =   u ( ixx, iyy, izzp1, 6, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 6, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 6, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfxccmodes ( ixx, iyy, izzp1, 5)
#endif

          bfy_R =   u ( ixx, iyy, izzp1, 7, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, 7, 4)
     1 + dt_shift * u ( ixx, iyy, izzp1, 7, 5)
#if ( DIV_FREE_RECONSTRUCT == 1)
     1    - ob_12 * bfyccmodes ( ixx, iyy, izzp1, 5)
#endif

          mean_mol_wt_R = mean_mol_wt_gr ( ixx, iyy, izzp1)
          gamma_eos_R = gamma_eos_gr ( ixx, iyy, izzp1)
          gamma_soundspeed_R = gamma_soundspeed_gr ( ixx, iyy, izzp1) 

! ----------

          rho_L = AMAX1 ( rhofloor, rho_L)
          a = 1.0 / rho_L

          vlx_L = vlx_L * a
          vly_L = vly_L * a
          vlz_L = vlz_L * a

          b = 0.5 * rho_L * ( vlx_L**2 + vly_L**2 + vlz_L**2)
     1      + ( bfx_L**2 + bfy_L**2 + bfzre**2) * ob_8pi

          prs_L = ( prs_L - b) * ( gamma_eos_L - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_L = entropygr ( ixx, iyy, izz, 1)
     1          + 0.5 * entropygr ( ixx, iyy, izz, 4)
     1     + dt_shift * entropygr ( ixx, iyy, izz, 5)
          IF ( prs_L .LE. prsfloor) prs_L = entropy_L
     1                              * rho_L**( gamma_eos_L - 1.0)
#endif
          prs_L = AMAX1 ( prsfloor, prs_L)

#if ( NFLUID >= 1)
          specie_L ( 1: NFLUID) =
     1             u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1     + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 4)
     1+ dt_shift * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 5)

          specie_L ( 1: NFLUID) = specie_L ( 1: NFLUID) * a
#endif

! ----------

          rho_R = AMAX1 ( rhofloor, rho_R)
          a = 1.0 / rho_R

          vlx_R = vlx_R * a
          vly_R = vly_R * a
          vlz_R = vlz_R * a

          b = 0.5 * rho_R * ( vlx_R**2 + vly_R**2 + vlz_R**2)
     1      + ( bfx_R**2 + bfy_R**2 + bfzre**2) * ob_8pi

          prs_R = ( prs_R - b) * ( gamma_eos_R - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropy_R = entropygr ( ixx, iyy, izzp1, 1)
     1          - 0.5 * entropygr ( ixx, iyy, izzp1, 4)
     1     + dt_shift * entropygr ( ixx, iyy, izzp1, 5)
          IF ( prs_R .LE. prsfloor) prs_R = entropy_R
     1                              * rho_R**( gamma_eos_R - 1.0)
#endif
          prs_R = AMAX1 ( prsfloor, prs_R)

#if ( NFLUID >= 1)
          specie_R ( 1: NFLUID) =
     1             u ( ixx, iyy, izzp1, nspecies_min: nspecies_max, 1)
     1     - 0.5 * u ( ixx, iyy, izzp1, nspecies_min: nspecies_max, 4)
     1+ dt_shift * u ( ixx, iyy, izzp1, nspecies_min: nspecies_max, 5)

          specie_R ( 1: NFLUID) = specie_R ( 1: NFLUID) * a
#endif

! ----------

          msonic = AMIN1 ( msonicgr ( ixx, iyy, izz),
     1                     msonicgr ( ixx, iyy, izzp1))
          divvel = AMAX1 ( ABS ( divvelgr ( ixx, iyy, izz)),
     1                     ABS ( divvelgr ( ixx, iyy, izzp1)) )

          a = divvel / ( flatten_coef * msonic) - 1.0

          IF ( a < 0.0 ) THEN
            blend_HLL = 0.0
          ELSE
            blend_HLL = AMIN1 ( 1.0, a)
          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Continuitive boundary conditions in the z-direction shown below.

          IF ( ( izz == nriembegin + one_z) .AND.
     1         ( bcarr_z ( ixx, iyy, 1) == 2) ) THEN

            rho_L = rho_R
            prs_L = prs_R

            vlx_L = vlx_R
            vly_L = vly_R
            vlz_L = vlz_R

            bfx_L = bfx_R
            bfy_L = bfy_R

#if ( NFLUID >= 1)
            specie_L ( 1: NFLUID) = specie_R ( 1: NFLUID)
#endif

          END IF

          IF ( ( izz == nriemend - one_z) .AND.
     1         ( bcarr_z ( ixx, iyy, 2) == 2) ) THEN

            rho_R = rho_L
            prs_R = prs_L

            vlx_R = vlx_L
            vly_R = vly_L
            vlz_R = vlz_L

            bfx_R = bfx_L
            bfy_R = bfy_L

#if ( NFLUID >= 1)
            specie_R ( 1: NFLUID) = specie_L ( 1: NFLUID)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Reflection boundary conditions in the z-direction shown below.

          IF ( ( izz == nriembegin + one_z) .AND.
     1         ( bcarr_z ( ixx, iyy, 1) == 3) ) THEN

            rho_L = rho_R
            prs_L = prs_R

            vlx_L = vlx_R
            vly_L = vly_R
            vlz_L= -vlz_R

            bfx_L= -bfx_R
            bfy_L= -bfy_R

#if ( NFLUID >= 1)
            specie_L ( 1: NFLUID) = specie_R ( 1: NFLUID)
#endif

          END IF

          IF ( ( izz == nriemend - one_z) .AND. 
     1         ( bcarr_z ( ixx, iyy, 2) == 3) ) THEN

            rho_R = rho_L
            prs_R = prs_L

            vlx_R = vlx_L
            vly_R = vly_L
            vlz_R= -vlz_L

            bfx_R= -bfx_L
            bfy_R= -bfy_L

#if ( NFLUID >= 1)
            specie_R ( 1: NFLUID) = specie_L ( 1: NFLUID)
#endif

          END IF

! ----------------------------------------------------------------------

! Call the Riemann solver with cyclic rotation of input variables.
! HLLEM is the default.

      CALL MHD_RIEM_HLLEM_PTWISE

     1     ( gamma, smallnum, prsfloor, rhofloor, pi,
     1       use_approx_eos, split_prsflux_re, iHLLEM_choice,

     1     gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1        rho_L, rho_R,
     1        prs_L, prs_R,
     1        vlz_L, vlz_R,
     1        vlx_L, vlx_R,
     1        vly_L, vly_R,
     1        bfzre,
     1        bfx_L, bfx_R,
     1        bfy_L, bfy_R,
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

! Fluxes have been evaluated at face centers.

          IF ( ( izz == nriembegin + one_z) .AND.
     1         ( bcarr_z ( ixx, iyy, 1) == 3) ) THEN

! ----------

            flux_z ( ixx, iyy, izz, 1) = flux_re ( 1)
            flux_z ( ixx, iyy, izz, 2) = flux_re ( 3)
            flux_z ( ixx, iyy, izz, 3) = flux_re ( 4)
            flux_z ( ixx, iyy, izz, 4) = flux_re ( 2)
            flux_z ( ixx, iyy, izz, 5) = flux_re ( 5)

            flux_z ( ixx, iyy, izz, 6) = flux_re ( 6)
            flux_z ( ixx, iyy, izz, 7) = flux_re ( 7)
            flux_z ( ixx, iyy, izz, 8) = 0.0

            prs_flux_z ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_z ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1                    flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_z ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          ELSE IF ( ( izz == nriemend - one_z) .AND.
     1         ( bcarr_z ( ixx, iyy, 2) == 3) ) THEN

! ----------

            flux_z ( ixx, iyy, izz, 1) = flux_re ( 1)
            flux_z ( ixx, iyy, izz, 2) = flux_re ( 3)
            flux_z ( ixx, iyy, izz, 3) = flux_re ( 4)
            flux_z ( ixx, iyy, izz, 4) = flux_re ( 2)
            flux_z ( ixx, iyy, izz, 5) = flux_re ( 5)

            flux_z ( ixx, iyy, izz, 6) = flux_re ( 6)
            flux_z ( ixx, iyy, izz, 7) = flux_re ( 7)
            flux_z ( ixx, iyy, izz, 8) = 0.0

            prs_flux_z ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_z ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1                    flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_z ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          ELSE

! ----------

            flux_z ( ixx, iyy, izz, 1) = flux_z ( ixx, iyy, izz, 1)
     1                                   + ob_6_oned * flux_re ( 1)

            flux_z ( ixx, iyy, izz, 2) = flux_z ( ixx, iyy, izz, 2)
     1                                   + ob_6_oned * flux_re ( 3)
            flux_z ( ixx, iyy, izz, 3) = flux_z ( ixx, iyy, izz, 3)
     1                                   + ob_6_oned * flux_re ( 4)
            flux_z ( ixx, iyy, izz, 4) = flux_z ( ixx, iyy, izz, 4)
     1                                   + ob_6_oned * flux_re ( 2)

            flux_z ( ixx, iyy, izz, 5) = flux_z ( ixx, iyy, izz, 5)
     1                                   + ob_6_oned * flux_re ( 5)

            flux_z ( ixx, iyy, izz, 6) = flux_z ( ixx, iyy, izz, 6)
     1                                   + ob_6_oned * flux_re ( 6)
            flux_z ( ixx, iyy, izz, 7) = flux_z ( ixx, iyy, izz, 7)
     1                                   + ob_6_oned * flux_re ( 7)
            flux_z ( ixx, iyy, izz, 8) = 0.0

            prs_flux_z ( ixx, iyy, izz) = flux_re ( 8)

#if ( NFLUID >= 1)
            flux_z ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1      flux_z ( ixx, iyy, izz, nspecies_min: nspecies_max)
     1           + ob_6_oned * flux_re ( nspecies_min: nspecies_max)
#endif

#if ( PRESSURE_POSITIVITY == 1)
            entropy_flux_z ( ixx, iyy, izz) = entropy_flux_re
#endif

! ----------

          END IF

! ----------------------------------------------------------------------

          END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_FLUX_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_LAPIDUS_X ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Adds a Lapidus flux contribution to the x-fluxes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, ixxp1, nriembegin, nriemend,
     1        ioffx_expand, ioffy_expand, ioffz_expand

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     msonic_arr, divvel_arr, lapidus_arr

! ----------------------------------------------------------------------

      IF ( lapidus_coef <= 0.0) RETURN
      IF ( ioffx .EQ. 0) RETURN

!      ioffx_expand = expand_by_one_zone * ioffx
!      ioffy_expand = expand_by_one_zone * ioffy
!      ioffz_expand = expand_by_one_zone * ioffz

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ixxp1,                       &
!$OMP&    nriembegin, nriemend,                                         &
!$OMP&    msonic_arr, divvel_arr, lapidus_arr)                          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    lapidus_threshold, lapidus_coef,                              &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    x_indx_limits, bcarr_x, u, flux_x,                            &
!$OMP&    msonicgr, divvelgr)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

! ----------------------------------------------------------------------

          nriembegin = x_indx_limits ( 1, iyy, izz) - 1 - ioffx_expand
          nriemend = x_indx_limits ( 2, iyy, izz) + ioffx_expand

          DO ixx = nriembegin, nriemend

          ixxp1 = ixx + 1

! ----------

          msonic_arr ( ixx) = 0.5 * ( msonicgr ( ixx, iyy, izz)
     1                              + msonicgr ( ixxp1, iyy, izz))
          divvel_arr ( ixx) = 0.5 * ( divvelgr ( ixx, iyy, izz)
     1                              + divvelgr ( ixxp1, iyy, izz))

! This triggers a Lapidus viscosity in regions of strong compression.
!          IF ( divvel_arr ( ixx) .LE.
!     1  - lapidus_threshold * msonic_arr ( ixx) ) THEN
!            lapidus_arr ( ixx) =
!     1      lapidus_coef * AMAX1 ( - divvel_arr ( ixx), 0.0)
!          ELSE
!            lapidus_arr ( ixx) = 0.0
!          END IF

! This triggers a Lapidus viscosity in regions of strong compression or
! rarefaction.
          IF ( ABS ( divvel_arr ( ixx) ) .GT.
     1    lapidus_threshold * msonic_arr ( ixx) ) THEN
            lapidus_arr ( ixx) =
     1      lapidus_coef * ABS ( divvel_arr ( ixx) )
          ELSE
            lapidus_arr ( ixx) = 0.0
          END IF

! ----------

          flux_x ( ixx, iyy, izz, :) = flux_x ( ixx, iyy, izz, :)
     1  + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, :, 1) - u ( ixxp1, iyy, izz, :, 1) )

! ----------

          END DO

! ----------------------------------------------------------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE ADD_LAPIDUS_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_LAPIDUS_Y ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Adds a Lapidus flux contribution to the y-fluxes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_y

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  2,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, iyyp1, nriembegin, nriemend,
     1        ioffx_expand, ioffy_expand, ioffz_expand

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     msonic_arr, divvel_arr, lapidus_arr

! ----------------------------------------------------------------------

      IF ( lapidus_coef <= 0.0) RETURN
      IF ( ioffy .EQ. 0) RETURN

!      ioffx_expand = expand_by_one_zone * ioffx
!      ioffy_expand = expand_by_one_zone * ioffy
!      ioffz_expand = expand_by_one_zone * ioffz

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, iyyp1,                       &
!$OMP&    nriembegin, nriemend,                                         &
!$OMP&    msonic_arr, divvel_arr, lapidus_arr)                          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    lapidus_threshold, lapidus_coef,                              &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    y_indx_limits, bcarr_y, u, flux_y,                            &
!$OMP&    msonicgr, divvelgr)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! ----------------------------------------------------------------------

          nriembegin = y_indx_limits ( ixx, 1, izz) - 1 - ioffy_expand
          nriemend = y_indx_limits ( ixx, 2, izz) + ioffy_expand

          DO iyy = nriembegin, nriemend

          iyyp1 = iyy + 1

! ----------

          msonic_arr ( iyy) = 0.5 * ( msonicgr ( ixx, iyy, izz)
     1                              + msonicgr ( ixx, iyyp1, izz))
          divvel_arr ( iyy) = 0.5 * ( divvelgr ( ixx, iyy, izz)
     1                              + divvelgr ( ixx, iyyp1, izz))

! This triggers a Lapidus viscosity in regions of strong compression.
!          IF ( divvel_arr ( iyy) .LE.
!     1  - lapidus_threshold * msonic_arr ( iyy) ) THEN
!            lapidus_arr ( iyy) =
!     1      lapidus_coef * AMAX1 ( - divvel_arr ( iyy), 0.0)
!          ELSE
!            lapidus_arr ( iyy) = 0.0
!          END IF

! This triggers a Lapidus viscosity in regions of strong compression or
! rarefaction.
          IF ( ABS ( divvel_arr ( iyy) ) .GT.
     1    lapidus_threshold * msonic_arr ( iyy) ) THEN
            lapidus_arr ( iyy) =
     1      lapidus_coef * ABS ( divvel_arr ( iyy) )
          ELSE
            lapidus_arr ( iyy) = 0.0
          END IF

! ----------

          flux_y ( ixx, iyy, izz, :) = flux_y ( ixx, iyy, izz, :)
     1  + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, :, 1) - u ( ixx, iyyp1, izz, :, 1) )

! ----------

          END DO

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE ADD_LAPIDUS_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_LAPIDUS_Z ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Adds a Lapidus flux contribution to the z-fluxes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_z

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  2) ::
     1         z_indx_limits, bcarr_z

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, izzp1, nriembegin, nriemend,
     1        ioffx_expand, ioffy_expand, ioffz_expand

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     msonic_arr, divvel_arr, lapidus_arr

! ----------------------------------------------------------------------

      IF ( lapidus_coef <= 0.0) RETURN
      IF ( ioffz .EQ. 0) RETURN

!      ioffx_expand = expand_by_one_zone * ioffx
!      ioffy_expand = expand_by_one_zone * ioffy
!      ioffz_expand = expand_by_one_zone * ioffz

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, izzp1,                       &
!$OMP&    nriembegin, nriemend,                                         &
!$OMP&    msonic_arr, divvel_arr, lapidus_arr)                          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    lapidus_threshold, lapidus_coef,                              &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    z_indx_limits, bcarr_z, u, flux_z,                            &
!$OMP&    msonicgr, divvelgr)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)                 &
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE ( 2)
#endif

      DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

        DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! ----------------------------------------------------------------------

          nriembegin = z_indx_limits ( ixx, iyy, 1) - 1 - ioffz_expand
          nriemend = z_indx_limits ( ixx, iyy, 2) + ioffz_expand

          DO izz = nriembegin, nriemend

          izzp1 = izz + 1

! ----------

          msonic_arr ( izz) = 0.5 * ( msonicgr ( ixx, iyy, izz)
     1                              + msonicgr ( ixx, iyy, izzp1))
          divvel_arr ( izz) = 0.5 * ( divvelgr ( ixx, iyy, izz)
     1                              + divvelgr ( ixx, iyy, izzp1))

! This triggers a Lapidus viscosity in regions of strong compression.
!          IF ( divvel_arr ( izz) .LE.
!     1  - lapidus_threshold * msonic_arr ( izz) ) THEN
!            lapidus_arr ( izz) =
!     1      lapidus_coef * AMAX1 ( - divvel_arr ( izz), 0.0)
!          ELSE
!            lapidus_arr ( izz) = 0.0
!          END IF

! This triggers a Lapidus viscosity in regions of strong compression or
! rarefaction.
          IF ( ABS ( divvel_arr ( izz) ) .GT.
     1    lapidus_threshold * msonic_arr ( izz) ) THEN
            lapidus_arr ( izz) =
     1      lapidus_coef * ABS ( divvel_arr ( izz) )
          ELSE
            lapidus_arr ( izz) = 0.0
          END IF

! ----------

          flux_z ( ixx, iyy, izz, :) = flux_z ( ixx, iyy, izz, :)
     1  + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, :, 1) - u ( ixx, iyy, izzp1, :, 1) )

! ----------

          END DO

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE ADD_LAPIDUS_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE UPDATE_VARS_PREDICTOR_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy, flux_x, flux_y, flux_z,
     1        prs_flux_x, prs_flux_y, prs_flux_z,
     1       entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z,

     1        efxgr, efygr, efzgr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This subroutine does the following tasks:

! 1) Makes a one-step update for "u" using the "flux_x, flux_y, flux_z"
! and "prs_flux_x, prs_flux_y, prs_flux_z" variables. If pressure positivity
! is called for, we also update "entropygr" using the
! "entropy_flux_x, entropy_flux_y, entropy_flux_z" variables.

! 2) Puts in the time-centered source terms in various geometries.
! External body forces have to be added in by hand.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_y, flux_z

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::

     1      efxgr, efygr, efzgr,

     1      prs_flux_x, prs_flux_y, prs_flux_z,
     1      entropy_flux_x, entropy_flux_y, entropy_flux_z

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits

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
     1        ixxp1, iyyp1, izzp1, ixx0, iyy0, izz0,
     1        ixxm1, iyym1, izzm1, nspecies_min, nspecies_max,
     1        ioffx_expand, ioffy_expand, ioffz_expand

      REAL ob_8pi, ob_4pi,

     1     div_min, div_max,
     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     a, b, c, d, e, f, vlx_accn, vly_accn, vlz_accn,

     1     rho_temp, vlx_temp, vly_temp, vlz_temp,
     1     bfx_temp, bfy_temp, bfz_temp,

     1     delx1, delx2, dely1, dely2, delz1, delz2, areaint

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

      IF ( ( expand_by_one_zone .LT. 0) .OR.
     1     ( expand_by_one_zone .GT. 1) ) THEN
      WRITE ( 6, *)
     1    "Stopping in UPDATE_VARS_PREDICTOR_EDGE ",
     1    ": expand_by_one_zone = ", expand_by_one_zone
      STOP
      END IF

!      ioffx_expand = expand_by_one_zone * ioffx
!      ioffy_expand = expand_by_one_zone * ioffy
!      ioffz_expand = expand_by_one_zone * ioffz

      IF ( ( ioffx > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffx_expand = 1
      ELSE
        ioffx_expand = 0
      END IF

      IF ( ( ioffy > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffy_expand = 1
      ELSE
        ioffy_expand = 0
      END IF

      IF ( ( ioffz > 0) .AND. ( expand_by_one_zone > 0) ) THEN
        ioffz_expand = 1
      ELSE
        ioffz_expand = 0
      END IF

! ----------------------------------------------------------------------
! Update the zone centered variables.
! ----------------------------------------------------------------------

! First, zero out the time rate of update, which we will then fill and
! use to update the conserved variables.
! The time rate of change of entropy does have to be zeroed out here.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&                     nspecies_max, u, entropygr)                  &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

          u ( ixx, iyy, izz, 1: nspecies_max, 5) = 0.0

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = 0.0
#endif

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Cartesian

      IF ( igeom .EQ. 1) THEN

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, flux_x, flux_y, flux_z, prs_flux_x, prs_flux_y, prs_flux_z,&
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! -----

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 1".

! ----------------------------------------------------------------------
! Cylindrical

      IF ( igeom .EQ. 2) THEN

! ------------------------------

! Transform the y-fluxes to an angular momentum flux.

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    flux_x, xb_save)                                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 0 - ioffx_expand, ix1 + ioffx_expand
          flux_x ( ixx, iyy, izz, 3) = flux_x ( ixx, iyy, izz, 3)
     1                               * xb_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               & 
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    flux_y, xc_save)                                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 0 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          flux_y ( ixx, iyy, izz, 3) = flux_y ( ixx, iyy, izz, 3)
     1                               * xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               & 
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    flux_z, xc_save)                                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          flux_z ( ixx, iyy, izz, 3) = flux_z ( ixx, iyy, izz, 3)
     1                               * xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, flux_x, flux_y, flux_z, prs_flux_x, prs_flux_y, prs_flux_z,&
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! -----

          a ( ixx) = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = xb_save ( ixx)

          c ( ixx) = xb_save ( ixx - 1)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Transform the y-angular momentum density time rate of update to a 
! momentum density.
! As long as the fluxes are not reused (no AMR for non-cartesian) this is ok.
! As long as the coordinate system does not touch the origin, it is easy
! to go from angular momentum fluxes to momentum fluxes.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    u, xc_save)                                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             / xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------
! Spherical

      IF ( igeom .EQ. 3) THEN

! ------------------------------

! Transform the z-momentum fluxes to an angular momentum density.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    flux_x, xb_save, yc_save)                                     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 0 - ioffx_expand, ix1 + ioffx_expand
          flux_x ( ixx, iyy, izz, 4) = flux_x ( ixx, iyy, izz, 4)
     1                       * xb_save ( ixx) * SIN ( yc_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    flux_y, xc_save, yb_save)                                     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 0 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          flux_y ( ixx, iyy, izz, 4) = flux_y ( ixx, iyy, izz, 4)
     1                       * xc_save ( ixx) * SIN ( yb_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    flux_z, xc_save, yc_save)                                     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          flux_z ( ixx, iyy, izz, 4) = flux_z ( ixx, iyy, izz, 4)
     1                       * xc_save ( ixx) * SIN ( yc_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, flux_x, flux_y, flux_z, prs_flux_x, prs_flux_y, prs_flux_z,&
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

! -----

          a ( ixx) = ( xb_save ( ixx)**3 - xb_save ( ixx - 1)**3 ) / 3.0
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = xb_save ( ixx)**2

          c ( ixx) = xb_save ( ixx - 1)**2

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = - xc_save ( ixx) * ( COS ( yb_save ( iyy))
     1                           - COS ( yb_save ( iyy - 1)) )
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = SIN ( yb_save ( iyy))

          c ( ixx) = SIN ( yb_save ( iyy - 1))

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1             * dz_save ( izz)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Transform the z-angular momentum density back to a momentum density.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    u, xc_save, yc_save)                                          &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                    / ( xc_save ( ixx) * SIN ( yc_save ( iyy)) )
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------
! For second order schemes, sources can be added in an operator split
! fashion. This is the place where one would put in the source terms.
! Note that the "rhogr" etc. variables are still from the old time
! while "u ( :, :, :, :, 1), bf?bdy ( :, :, :, 1)" have been updated.
! This has been used to advantage here to build time-centered sources.
! ----------------------------------------------------------------------
! Cartesian

#if ( 0 == 1)
      IF ( igeom .EQ. 1) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c, rho_temp,           &
!$OMP&    vlx_temp, vly_temp, vlz_temp, bfx_temp, bfy_temp, bfz_temp,   & 
!$OMP&    vlx_accn, vly_accn, vlz_accn)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1                             + vlx_accn ( ixx)
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             + vly_accn ( ixx)
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          a ( ixx) = vlx_accn ( ixx) * vlx_temp ( ixx)
     1             + vly_accn ( ixx) * vly_temp ( ixx)
     1             + vlz_accn ( ixx) * vlz_temp ( ixx)

          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + a ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 1".
#endif

! ----------------------------------------------------------------------
! Cylindrical

      IF ( igeom .EQ. 2) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, izzm1, iyym1, ixxm1,         &
!$OMP&    a, b, c, rho_temp, vlx_temp, vly_temp, vlz_temp,              &
!$OMP&    bfx_temp, bfy_temp, bfz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

          bfx_temp ( ixx) = 0.25 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                    + bfxbdy ( ixxm1, iyy, izz, 1))
     1              + 0.5 * bfxgr ( ixx, iyy, izz)
          bfy_temp ( ixx) = 0.25 * ( bfybdy ( ixx, iyy, izz, 1)
     1                    + bfybdy ( ixx, iyym1, izz, 1))
     1              + 0.5 * bfygr ( ixx, iyy, izz)
          bfz_temp ( ixx) = 0.25 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                    + bfzbdy ( ixx, iyy, izzm1, 1))
     1              + 0.5 * bfzgr ( ixx, iyy, izz)

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

            b ( ixx) =  vlx_accn ( ixx)
     1      + ( rho_temp ( ixx) * vly_temp ( ixx)**2
     1        - ob_4pi * bfy_temp ( ixx)**2 )
     1      * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + b ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             + vly_accn ( ixx)

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          c ( ixx) = vlx_accn ( ixx) * vlx_temp ( ixx)
     1             + vly_accn ( ixx) * vly_temp ( ixx)
     1             + vlz_accn ( ixx) * vlz_temp ( ixx)

          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + c ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------
! Spherical

      IF ( igeom .EQ. 3) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, izzm1, iyym1, ixxm1,         &
!$OMP&    a, b, c, d, e, f, rho_temp, vlx_temp, vly_temp, vlz_temp,     &
!$OMP&    bfx_temp, bfy_temp, bfz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ioffx_expand, ioffy_expand, ioffz_expand,                     &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz_expand, iz1 + ioffz_expand
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - ioffy_expand, iy1 + ioffy_expand
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - ioffx_expand, ix1 + ioffx_expand
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

          bfx_temp ( ixx) = 0.25 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                    + bfxbdy ( ixxm1, iyy, izz, 1))
     1              + 0.5 * bfxgr ( ixx, iyy, izz)
          bfy_temp ( ixx) = 0.25 * ( bfybdy ( ixx, iyy, izz, 1)
     1                    + bfybdy ( ixx, iyym1, izz, 1))
     1              + 0.5 * bfygr ( ixx, iyy, izz)
          bfz_temp ( ixx) = 0.25 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                    + bfzbdy ( ixx, iyy, izzm1, 1))
     1              + 0.5 * bfzgr ( ixx, iyy, izz)

          c ( ixx) = COS ( yc_save ( iyy)) / SIN ( yc_save ( iyy))

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

          d ( ixx) = + vlx_accn ( ixx)
     1    + ( rho_temp ( ixx) * ( vly_temp ( ixx)**2
     1                          + vlz_temp ( ixx)**2 )
     1        - ob_4pi * ( bfy_temp ( ixx)**2 + bfz_temp ( ixx)**2 ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + d ( ixx)

          e ( ixx) = vly_accn ( ixx)
     1    - ( rho_temp ( ixx) * ( vlx_temp ( ixx) * vly_temp ( ixx)
     1                          - vlz_temp ( ixx)**2 * c ( ixx) )
     1        - ob_4pi * ( bfx_temp ( ixx) * bfy_temp ( ixx)
     1                   - bfz_temp ( ixx)**2 * c ( ixx) ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5) + e ( ixx)

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          f ( ixx) = vlx_accn ( ixx) * vlx_temp ( ixx)
     1             + vly_accn ( ixx) * vly_temp ( ixx)
     1             + vlz_accn ( ixx) * vlz_temp ( ixx)

          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + f ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE UPDATE_VARS_PREDICTOR_EDGE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_PREDICTOR_SOURCE_TERMS

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, 

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and adds in source terms to the time rate of update.
! This is stored in "u ( ixx, iyy, izz, :, 5)" within each zone.

! It is effectively Stage I, i.e. the predictor stage, of a two-stage
! update.

! External body forces have to be added in by hand.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

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

! Variables past this point are not used in the subroutine interface.

      INTEGER ixx, iyy, izz, nspecies_min, nspecies_max

      REAL ob_8pi, ob_4pi

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b, c,
     1        vlx_accn, vly_accn, vlz_accn

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

! ----------------------------------------------------------------------

! Only external body forces need to be included in cartesian geometry.
! This part is just an instantiation of how acceleration terms contribute
! to the momentum and energy densities in the predictor part of the timestep.

#if ( 0 == 1)
      IF ( igeom .EQ. 1) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    vlx_accn, vly_accn, vlz_accn)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            & 
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0  
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

          vlx_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
     1                    * vlx_accn ( ixx)
          vly_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
     1                    * vly_accn ( ixx)
          vlz_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
     1                    * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1                             + vlx_accn ( ixx)
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             + vly_accn ( ixx)
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5)
     1                      + vlx_accn ( ixx) * vlxgr ( ixx, iyy, izz)
     1                      + vly_accn ( ixx) * vlygr ( ixx, iyy, izz)
     1                      + vlz_accn ( ixx) * vlzgr ( ixx, iyy, izz)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 1".
#endif

! ----------------------------------------------------------------------

! Now that the contribution from the flux terms has been accounted for,
! pick up the contribution from the source terms in cylindrical geometry.
! Also instantiates how body forces are to be included.

      IF ( igeom .EQ. 2) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    vlx_accn, vly_accn, vlz_accn)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            & 
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
!     1                    * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
!     1                    * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
!     1                    * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1      + vlx_accn ( ixx)
     1      + ( rhogr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)**2
     1        - ob_4pi * bfygr ( ixx, iyy, izz)**2 )
     1      * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1      + vly_accn ( ixx)
     1      - ( rhogr ( ixx, iyy, izz) * vlxgr ( ixx, iyy, izz)
     1                                 * vlygr ( ixx, iyy, izz)
     1      - ob_4pi * bfxgr ( ixx, iyy, izz) * bfygr ( ixx, iyy, izz))
     1      * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5)
     1                      + vlx_accn ( ixx) * vlxgr ( ixx, iyy, izz)
     1                      + vly_accn ( ixx) * vlygr ( ixx, iyy, izz)
     1                      + vlz_accn ( ixx) * vlzgr ( ixx, iyy, izz)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------

! Now that the contribution from the flux terms has been accounted for,
! pick up the contribution from the source terms in spherical geometry.
! Also instantiates how body forces are to be included.

      IF ( igeom .EQ. 3) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    vlx_accn, vly_accn, vlz_accn)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            & 
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

          c ( ixx) = COS ( yc_save ( iyy)) / SIN ( yc_save ( iyy))

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
!     1                    * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
!     1                    * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rhogr ( ixx, iyy, izz)
!     1                    * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1    + vlx_accn ( ixx)
     1    + ( rhogr ( ixx, iyy, izz) * ( vlygr ( ixx, iyy, izz)**2
     1                                 + vlzgr ( ixx, iyy, izz)**2 )
     1      - ob_4pi * ( bfygr ( ixx, iyy, izz)**2 
     1                 + bfzgr ( ixx, iyy, izz)**2 ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1    + vly_accn ( ixx) - ( rhogr ( ixx, iyy, izz)
     1    * ( vlxgr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)
     1      - vlzgr ( ixx, iyy, izz)**2 * c ( ixx) )
     1      - ob_4pi * ( bfxgr ( ixx, iyy, izz) * bfygr ( ixx, iyy, izz)
     1                 - bfzgr ( ixx, iyy, izz)**2 * c ( ixx) ) )
     1    * dtcur / xc_save ( ixx)

          a = rhogr ( ixx, iyy, izz) * vlxgr ( ixx, iyy, izz)
     1                               * vlzgr ( ixx, iyy, izz)
     1      - ob_4pi * bfxgr ( ixx, iyy, izz) * bfzgr ( ixx, iyy, izz)

          b ( ixx) = rhogr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)
     1                               * vlzgr ( ixx, iyy, izz)
     1      - ob_4pi * bfygr ( ixx, iyy, izz) * bfzgr ( ixx, iyy, izz)

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1    + vlz_accn ( ixx)
     1    - dtcur * ( a ( ixx) + b ( ixx) * c ( ixx)) / xc_save ( ixx)

          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5)
     1                      + vlx_accn ( ixx) * vlxgr ( ixx, iyy, izz)
     1                      + vly_accn ( ixx) * vlygr ( ixx, iyy, izz)
     1                      + vlz_accn ( ixx) * vlzgr ( ixx, iyy, izz)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE ADD_PREDICTOR_SOURCE_TERMS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE UPDATE_VARS_CORRECTOR_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy, flux_x, flux_y, flux_z,
     1        prs_flux_x, prs_flux_y, prs_flux_z,
     1       entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z,

     1        efxgr, efygr, efzgr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This subroutine does the following tasks:

! 1) Makes a one-step update for "u" using the "flux_x, flux_y, flux_z"
! and "prs_flux_x, prs_flux_y, prs_flux_z" variables. If pressure positivity
! is called for, we also update "entropygr" using the
! "entropy_flux_x, entropy_flux_y, entropy_flux_z" variables.

! 2) Updates "bfxbdy, bfybdy, bfzbdy" using the "efxgr, efygr, efzgr"
! variables.

! 3) Puts in the time-centered source terms in various geometries.
! External body forces have to be added in by hand.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_y, flux_z

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::

     1      efxgr, efygr, efzgr,

     1      prs_flux_x, prs_flux_y, prs_flux_z,
     1      entropy_flux_x, entropy_flux_y, entropy_flux_z

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits

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
     1        ixxp1, iyyp1, izzp1, ixx0, iyy0, izz0,
     1        ixxm1, iyym1, izzm1, nspecies_min, nspecies_max

      REAL ob_8pi, ob_4pi,

     1     div_min, div_max,
     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     a, b, c, d, e, f, vlx_accn, vly_accn, vlz_accn,

     1     rho_temp, vlx_temp, vly_temp, vlz_temp,
     1     bfx_temp, bfy_temp, bfz_temp,

     1     delx1, delx2, dely1, dely2, delz1, delz2, areaint

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

      IF ( ( expand_by_one_zone .LT. 0) .OR.
     1     ( expand_by_one_zone .GT. 1) ) THEN
      WRITE ( 6, *)
     1    "Stopping in UPDATE_VARS_CORRECTOR_EDGE ",
     1    ": expand_by_one_zone = ", expand_by_one_zone
      STOP
      END IF

! ----------------------------------------------------------------------
! Update the zone centered variables.
! ----------------------------------------------------------------------

! First, zero out the time rate of update, which we will then fill and
! use to update the conserved variables.
! The time rate of change of entropy does have to be zeroed out here.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, nspecies_max, u, entropygr)   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          u ( ixx, iyy, izz, 1: nspecies_max, 5) = 0.0

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = 0.0
#endif

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Cartesian

      IF ( igeom .EQ. 1) THEN

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, flux_x, flux_y, flux_z, prs_flux_x, prs_flux_y, prs_flux_z,&
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! -----

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 1".

! ----------------------------------------------------------------------
! Cylindrical

      IF ( igeom .EQ. 2) THEN

! ------------------------------

! Transform the y-fluxes to an angular momentum flux.

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, flux_x, xb_save)              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 0, ix1
          flux_x ( ixx, iyy, izz, 3) = flux_x ( ixx, iyy, izz, 3)
     1                               * xb_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, flux_y, xc_save)              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 0, iy1
          DO ixx = 1, ix1
          flux_y ( ixx, iyy, izz, 3) = flux_y ( ixx, iyy, izz, 3)
     1                               * xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, flux_z, xc_save)              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          flux_z ( ixx, iyy, izz, 3) = flux_z ( ixx, iyy, izz, 3)
     1                               * xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, flux_x, flux_y, flux_z, prs_flux_x, prs_flux_y, prs_flux_z,&
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! -----

          a ( ixx) = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = xb_save ( ixx)

          c ( ixx) = xb_save ( ixx - 1)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Transform the y-angular momentum density time rate of update to a 
! momentum density.
! As long as the fluxes are not reused (no AMR for non-cartesian) this is ok.
! As long as the coordinate system does not touch the origin, it is easy
! to go from angular momentum fluxes to momentum fluxes.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, u, xc_save)                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             / xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------
! Spherical

      IF ( igeom .EQ. 3) THEN

! ------------------------------

! Transform the z-momentum fluxes to an angular momentum density.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, flux_x, xb_save, yc_save)     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 0, ix1
          flux_x ( ixx, iyy, izz, 4) = flux_x ( ixx, iyy, izz, 4)
     1                       * xb_save ( ixx) * SIN ( yc_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, flux_y, xc_save, yb_save)     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 0, iy1
          DO ixx = 1, ix1
          flux_y ( ixx, iyy, izz, 4) = flux_y ( ixx, iyy, izz, 4)
     1                       * xc_save ( ixx) * SIN ( yb_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, flux_z, xc_save, yc_save)     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          flux_z ( ixx, iyy, izz, 4) = flux_z ( ixx, iyy, izz, 4)
     1                       * xc_save ( ixx) * SIN ( yc_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, flux_x, flux_y, flux_z, prs_flux_x, prs_flux_y, prs_flux_z,&
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! -----

          a ( ixx) = ( xb_save ( ixx)**3 - xb_save ( ixx - 1)**3 ) / 3.0
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = xb_save ( ixx)**2

          c ( ixx) = xb_save ( ixx - 1)**2

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = - xc_save ( ixx) * ( COS ( yb_save ( iyy))
     1                           - COS ( yb_save ( iyy - 1)) )
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = SIN ( yb_save ( iyy))

          c ( ixx) = SIN ( yb_save ( iyy - 1))

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1             * dz_save ( izz)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Transform the z-angular momentum density back to a momentum density.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, u, xc_save, yc_save)          &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                    / ( xc_save ( ixx) * SIN ( yc_save ( iyy)) )
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 3".

! ------------------------------
! Now, actually make the update

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, nspecies_max, u, entropygr)   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1

          u ( ixx, iyy, izz, 1: nspecies_max, 1) = 
     1          u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1        + u ( ixx, iyy, izz, 1: nspecies_max, 5)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1                                  + entropygr ( ixx, iyy, izz, 5)
#endif

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Now that the edge-centered electric fields are built, update the
! face-centered magnetic fields.
! ----------------------------------------------------------------------

      IF ( igeom .EQ. 1) THEN

! ------------------------------

! Update "bfxbdy" in Cartesian geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, dy_save, dz_save,      &
!$OMP&    bfxbdy, efygr, efzgr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 0, ix1

          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - dy_save ( iyy) * efygr ( ixx, iyy, izz)
     1              + dy_save ( iyy) * efygr ( ixx, iyy, izz - 1)
     1              + dz_save ( izz) * efzgr ( ixx, iyy, izz)
     1              - dz_save ( izz) * efzgr ( ixx, iyy - 1, izz)
     1             ) / ( dy_save ( iyy) * dz_save ( izz) )

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ----------

! Update "bfybdy" in Cartesian geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, dx_save, dz_save,      &
!$OMP&    bfybdy, efxgr, efzgr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 0, iy1

          DO ixx = 1, ix1

          bfybdy ( ixx, iyy, izz, 1) = bfybdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( dx_save ( ixx) * efxgr ( ixx, iyy, izz)
     1            - dx_save ( ixx) * efxgr ( ixx, iyy, izz - 1)
     1            - dz_save ( izz) * efzgr ( ixx, iyy, izz)
     1            + dz_save ( izz) * efzgr ( ixx - 1, iyy, izz)
     1             ) / ( dx_save ( ixx) * dz_save ( izz) )

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ----------

! Update "bfzbdy" in Cartesian geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, dx_save, dy_save,      &
!$OMP&    bfzbdy, efxgr, efygr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - dx_save ( ixx) * efxgr ( ixx, iyy, izz)
     1              + dx_save ( ixx) * efxgr ( ixx, iyy - 1, izz)
     1              + dy_save ( iyy) * efygr ( ixx, iyy, izz)
     1              - dy_save ( iyy) * efygr ( ixx - 1, iyy, izz)
     1             ) / ( dx_save ( ixx) * dy_save ( iyy) )

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 1".

! ----------------------------------------------------------------------

      IF ( igeom .EQ. 2) THEN

! ------------------------------

! Update "bfxbdy" in Cylindrical geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, dely1, dely2, delz1, delz2,  &
!$OMP&    areaint)                                                      &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, smallnum,              &
!$OMP&    xb_save, dy_save, dz_save,                                    &
!$OMP&    bfxbdy, efygr, efzgr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 0, ix1

          dely1 ( ixx) = xb_save ( ixx) * dy_save ( iyy)
          dely2 ( ixx) = xb_save ( ixx) * dy_save ( iyy)
          delz1 ( ixx) = dz_save ( izz)
          delz2 ( ixx) = dz_save ( izz)

          areaint ( ixx) = xb_save ( ixx) * dy_save ( iyy)
     1                   * dz_save ( izz)

          IF ( areaint ( ixx) .GT. smallnum) THEN

          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - dely1 ( ixx) * efygr ( ixx, iyy, izz)
     1              + dely2 ( ixx) * efygr ( ixx, iyy, izz - 1)
     1              + delz1 ( ixx) * efzgr ( ixx, iyy, izz)
     1              - delz2 ( ixx) * efzgr ( ixx, iyy - 1, izz)
     1             ) / areaint ( ixx)

          END IF

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ----------

! Update "bfybdy" in Cylindrical geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, dx_save, dz_save,      &
!$OMP&    bfybdy, efxgr, efzgr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 0, iy1

          DO ixx = 1, ix1

          bfybdy ( ixx, iyy, izz, 1) = bfybdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( dx_save ( ixx) * efxgr ( ixx, iyy, izz)
     1            - dx_save ( ixx) * efxgr ( ixx, iyy, izz - 1)
     1            - dz_save ( izz) * efzgr ( ixx, iyy, izz)
     1            + dz_save ( izz) * efzgr ( ixx - 1, iyy, izz)
     1             ) / ( dx_save ( ixx) * dz_save ( izz) )

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ----------

! Update "bfzbdy" in Cylindrical geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, delx1, delx2, dely1, dely2,  &
!$OMP&    areaint)                                                      &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur,                        &
!$OMP&    dx_save, xb_save, dy_save,                                    &
!$OMP&    bfzbdy, efxgr, efygr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          delx1 ( ixx) = dx_save ( ixx)
          delx2 ( ixx) = dx_save ( ixx)
          dely1 ( ixx) = xb_save ( ixx) * dy_save ( iyy)
          dely2 ( ixx) = xb_save ( ixx - 1) * dy_save ( iyy)

          areaint ( ixx) = 0.5
     1  * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2) * dy_save ( iyy)

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - delx1 ( ixx) * efxgr ( ixx, iyy, izz)
     1              + delx2 ( ixx) * efxgr ( ixx, iyy - 1, izz)
     1              + dely1 ( ixx) * efygr ( ixx, iyy, izz)
     1              - dely2 ( ixx) * efygr ( ixx - 1, iyy, izz)
     1             ) / areaint ( ixx)

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------

      IF ( igeom .EQ. 3) THEN

! ------------------------------

! Update "bfxbdy" in Spherical geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, dely1, dely2, delz1, delz2,  &
!$OMP&    areaint)                                                      &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, smallnum,              &
!$OMP&    xb_save, dy_save, yb_save, dz_save,                           &
!$OMP&    bfxbdy, efygr, efzgr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 0, ix1

          dely1 ( ixx) = xb_save ( ixx) * dy_save ( iyy)
          dely2 ( ixx) = xb_save ( ixx) * dy_save ( iyy)
          delz1 ( ixx) = xb_save ( ixx) * SIN ( yb_save ( iyy))
     1          * dz_save ( izz)
          delz2 ( ixx) = xb_save ( ixx) * SIN ( yb_save ( iyy - 1))
     1          * dz_save ( izz)

          areaint ( ixx) = - xb_save ( ixx)**2 * dz_save ( izz)
     1         * ( COS ( yb_save ( iyy)) - COS ( yb_save ( iyy - 1)) )

          IF ( areaint ( ixx) .GT. smallnum) THEN

          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - dely1 ( ixx) * efygr ( ixx, iyy, izz)
     1              + dely2 ( ixx) * efygr ( ixx, iyy, izz - 1)
     1              + delz1 ( ixx) * efzgr ( ixx, iyy, izz)
     1              - delz2 ( ixx) * efzgr ( ixx, iyy - 1, izz)
     1             ) / areaint ( ixx)

          END IF

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ----------

! Update "bfybdy" in Spherical geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, delx1, delx2, delz1, delz2,  &
!$OMP&    areaint)                                                      &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur, smallnum,              &
!$OMP&    dx_save, xb_save, yb_save, dz_save,                           &
!$OMP&    bfybdy, efxgr, efzgr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 0, iy1

          DO ixx = 1, ix1

          delx1 ( ixx) = dx_save ( ixx)
          delx2 ( ixx) = dx_save ( ixx)
          delz1 ( ixx) = xb_save ( ixx) * SIN ( yb_save ( iyy))
     1          * dz_save ( izz)
          delz2 ( ixx) = xb_save ( ixx - 1) * SIN ( yb_save ( iyy))
     1          * dz_save ( izz)

          areaint ( ixx) =
     1        0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
     1            * SIN ( yb_save ( iyy)) * dz_save ( izz)

          IF ( areaint ( ixx) .GT. smallnum) THEN

          bfybdy ( ixx, iyy, izz, 1) = bfybdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( delx1 ( ixx) * efxgr ( ixx, iyy, izz)
     1            - delx2 ( ixx) * efxgr ( ixx, iyy, izz - 1)
     1            - delz1 ( ixx) * efzgr ( ixx, iyy, izz)
     1            + delz2 ( ixx) * efzgr ( ixx - 1, iyy, izz)
     1             ) / areaint ( ixx)

          END IF

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ----------

! Update "bfzbdy" in Spherical geometry.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, delx1, delx2, dely1, dely2,  & 
!$OMP&    areaint)                                                      &
!$OMP&            SHARED ( ix1, iy1, iz1, dtcur,                        & 
!$OMP&    dx_save, xb_save, dy_save,                                    &
!$OMP&    bfzbdy, efxgr, efygr)                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          delx1 ( ixx) = dx_save ( ixx)
          delx2 ( ixx) = dx_save ( ixx)
          dely1 ( ixx) = xb_save ( ixx) * dy_save ( iyy)
          dely2 ( ixx) = xb_save ( ixx - 1) * dy_save ( iyy)

          areaint ( ixx) =
     1        0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
     1            * dy_save ( iyy)

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - delx1 ( ixx) * efxgr ( ixx, iyy, izz)
     1              + delx2 ( ixx) * efxgr ( ixx, iyy - 1, izz)
     1              + dely1 ( ixx) * efygr ( ixx, iyy, izz)
     1              - dely2 ( ixx) * efygr ( ixx - 1, iyy, izz)
     1             ) / areaint ( ixx)

          END DO

        END DO

      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------
! For second order schemes, sources can be added in an operator split
! fashion. This is the place where one would put in the source terms.
! Note that the "rhogr" etc. variables are still from the old time
! while "u ( :, :, :, :, 1), bf?bdy ( :, :, :, 1)" have been updated.
! This has been used to advantage here to build time-centered sources.
! ----------------------------------------------------------------------
! Cartesian

#if ( 0 == 1)
      IF ( igeom .EQ. 1) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c, rho_temp,           &
!$OMP&    vlx_temp, vly_temp, vlz_temp, bfx_temp, bfy_temp, bfz_temp,   & 
!$OMP&    vlx_accn, vly_accn, vlz_accn)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1)
     1                             + vlx_accn ( ixx)
          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1                             + vly_accn ( ixx)
          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1                             + vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5)
     1                             + vlx_accn ( ixx)
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             + vly_accn ( ixx)
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          a ( ixx) = vlx_accn ( ixx) * vlx_temp ( ixx)
     1             + vly_accn ( ixx) * vly_temp ( ixx)
     1             + vlz_accn ( ixx) * vlz_temp ( ixx)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1) + a ( ixx)
          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + a ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 1".
#endif

! ----------------------------------------------------------------------
! Cylindrical

      IF ( igeom .EQ. 2) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, izzm1, iyym1, ixxm1,         &
!$OMP&    a, b, c, rho_temp, vlx_temp, vly_temp, vlz_temp,              &
!$OMP&    bfx_temp, bfy_temp, bfz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

          bfx_temp ( ixx) = 0.25 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                    + bfxbdy ( ixxm1, iyy, izz, 1))
     1              + 0.5 * bfxgr ( ixx, iyy, izz)
          bfy_temp ( ixx) = 0.25 * ( bfybdy ( ixx, iyy, izz, 1)
     1                    + bfybdy ( ixx, iyym1, izz, 1))
     1              + 0.5 * bfygr ( ixx, iyy, izz)
          bfz_temp ( ixx) = 0.25 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                    + bfzbdy ( ixx, iyy, izzm1, 1))
     1              + 0.5 * bfzgr ( ixx, iyy, izz)

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

            b ( ixx) =  vlx_accn ( ixx)
     1      + ( rho_temp ( ixx) * vly_temp ( ixx)**2
     1        - ob_4pi * bfy_temp ( ixx)**2 )
     1      * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1) + b ( ixx)
          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + b ( ixx)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1                             + vly_accn ( ixx)
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1                             + vly_accn ( ixx)

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1                             + vlz_accn ( ixx)
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          c ( ixx) = vlx_accn ( ixx) * vlx_temp ( ixx)
     1             + vly_accn ( ixx) * vly_temp ( ixx)
     1             + vlz_accn ( ixx) * vlz_temp ( ixx)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1) + c ( ixx)
          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + c ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------
! Spherical

      IF ( igeom .EQ. 3) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, izzm1, iyym1, ixxm1,         &
!$OMP&    a, b, c, d, e, f, rho_temp, vlx_temp, vly_temp, vlz_temp,     &
!$OMP&    bfx_temp, bfy_temp, bfz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, ob_8pi, ob_4pi, nspecies_min, nspecies_max,            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

          bfx_temp ( ixx) = 0.25 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                    + bfxbdy ( ixxm1, iyy, izz, 1))
     1              + 0.5 * bfxgr ( ixx, iyy, izz)
          bfy_temp ( ixx) = 0.25 * ( bfybdy ( ixx, iyy, izz, 1)
     1                    + bfybdy ( ixx, iyym1, izz, 1))
     1              + 0.5 * bfygr ( ixx, iyy, izz)
          bfz_temp ( ixx) = 0.25 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                    + bfzbdy ( ixx, iyy, izzm1, 1))
     1              + 0.5 * bfzgr ( ixx, iyy, izz)

          c ( ixx) = COS ( yc_save ( iyy)) / SIN ( yc_save ( iyy))

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

          d ( ixx) = + vlx_accn ( ixx)
     1    + ( rho_temp ( ixx) * ( vly_temp ( ixx)**2
     1                          + vlz_temp ( ixx)**2 )
     1        - ob_4pi * ( bfy_temp ( ixx)**2 + bfz_temp ( ixx)**2 ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1) + d ( ixx)
          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + d ( ixx)

          e ( ixx) = vly_accn ( ixx)
     1    - ( rho_temp ( ixx) * ( vlx_temp ( ixx) * vly_temp ( ixx)
     1                          - vlz_temp ( ixx)**2 * c ( ixx) )
     1        - ob_4pi * ( bfx_temp ( ixx) * bfy_temp ( ixx)
     1                   - bfz_temp ( ixx)**2 * c ( ixx) ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1) + e ( ixx)
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5) + e ( ixx)

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1                             + vlz_accn ( ixx)
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1                             + vlz_accn ( ixx)

          f ( ixx) = vlx_accn ( ixx) * vlx_temp ( ixx)
     1             + vly_accn ( ixx) * vly_temp ( ixx)
     1             + vlz_accn ( ixx) * vlz_temp ( ixx)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1) + f ( ixx)
          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + f ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE UPDATE_VARS_CORRECTOR_EDGE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EVALUATE_TIMESTEP

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        dtcur, dtcurp1,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Finds a new timestep. "dtcurp1" is loaded with the next timestep.
! Pressure positivity, if it is called for, is also enforced here.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe

      REAL dtcur, dtcurp1

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr

! ----------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits

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
     1        ixxp1, iyyp1, izzp1, ixx0, iyy0, izz0,
     1        ixxm1, iyym1, izzm1, nspecies_min, nspecies_max

      REAL ob_8pi, ob_4pi

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     a, b, c, d, e

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

! ----------------------------------------------------------------------
! Evaluate "dtcurp1", i.e. the new timestep, by directly minimizing over
! the updated flow variables.
! Also, make all zone-centered primitive variables here. That way they are
! available at the appropriate location in the conserved variables as and
! when one needs them.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d, e)                             &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, ob_8pi, cfl_coef_sdsp, cfl_coef_advct,                 &
!$OMP&    rhofloor, prsfloor,                                           &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    x_indx_limits, u, bfxbdy, bfybdy, bfzbdy, entropygr,          & 
!$OMP&    rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,       &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, msonicgr)                  &
!$OMP&            REDUCTION ( MIN: dtcurp1)                             &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)

        DO iyy = 1, iy1
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)

          DO ixx = x_indx_limits ( 1, iyy, izz),
     1             x_indx_limits ( 2, iyy, izz)
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

          u ( ixx, iyy, izz, 1, 1) = AMAX1 ( rhofloor,
     1                               u ( ixx, iyy, izz, 1, 1))

          rhogr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 1, 1)

          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)

          vlxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 2, 1) * a ( ixx)
          vlygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 3, 1) * a ( ixx)
          vlzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 4, 1) * a ( ixx)

          bfxgr ( ixx, iyy, izz) =0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                                  + bfxbdy ( ixxm1, iyy, izz, 1))
          bfygr ( ixx, iyy, izz) =0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                                  + bfybdy ( ixx, iyym1, izz, 1))
          bfzgr ( ixx, iyy, izz) =0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                                  + bfzbdy ( ixx, iyy, izzm1, 1))

          u ( ixx, iyy, izz, 6, 1) = bfxgr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 7, 1) = bfygr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 8, 1) = bfzgr ( ixx, iyy, izz)

          e ( ixx) = ( vlxgr ( ixx, iyy, izz)**2
     1        + vlygr ( ixx, iyy, izz)**2 + vlzgr ( ixx, iyy, izz)**2)

          b ( ixx) = 0.5 * rhogr ( ixx, iyy, izz) * e ( ixx)

          c ( ixx) = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          prsgr ( ixx, iyy, izz) = ( u ( ixx, iyy, izz, 5, 1)
     1                             - b ( ixx) - c ( ixx) )
     1                         * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

! ----------

#if ( PRESSURE_POSITIVITY == 1)
          IF ( prsgr ( ixx, iyy, izz) .LE. prsfloor) THEN

! If the pressure that is evaluated from conserved variables is negative then,
! and then only, replace it with a pressure that is evaluated from the
! entropy density. This is the ONLY place where the energy density is
! "touched up". Thus, this is the ONLY place where we trade conservation
! in return for physical consistency.

          d ( ixx) = prsgr ( ixx, iyy, izz)

          prsgr ( ixx, iyy, izz) = entropygr ( ixx, iyy, izz, 1)
     1 * rhogr ( ixx, iyy, izz)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1)
     1    + ( prsgr ( ixx, iyy, izz) - d ( ixx) )
     1    / ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

          END IF
#endif

! Even without pressure positivity, we minimally elevate the energy density
! so that the pressure has its floor value.

          IF ( prsgr ( ixx, iyy, izz) .LT. prsfloor) THEN
          prsgr ( ixx, iyy, izz) = prsfloor
          u ( ixx, iyy, izz, 5, 1) = b ( ixx) + c ( ixx) + prsfloor
     1                         / ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
          END IF

! ----------

          msonicgr ( ixx, iyy, izz) =
     1    ( 2.0 * c ( ixx) + gamma_soundspeed_gr ( ixx, iyy, izz)
     1              * prsgr ( ixx, iyy, izz) ) * a ( ixx)
          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

          e ( ixx) = SQRT ( e ( ixx))

! ----------

          IF ( igeom .EQ. 1) THEN

            d ( ixx) = AMIN1 ( dx_save ( ixx), 
     1                         dy_save ( iyy), dz_save ( izz))

          ELSE IF ( igeom .EQ. 2) THEN

            d ( ixx) = AMIN1 ( dx_save ( ixx), 
     1                         xc_save ( ixx) * dy_save ( iyy),
     1                         dz_save ( izz) )

          ELSE IF ( igeom .EQ. 3) THEN

            d ( ixx) = AMIN1 ( dx_save ( ixx), 
     1                         xc_save ( ixx) * dy_save ( iyy),
     1        xc_save ( ixx) * SIN ( yc_save ( iyy)) * dz_save ( izz) )

          END IF

! ----------

          d ( ixx) = d ( ixx)
     1             / ( msonicgr ( ixx, iyy, izz) / cfl_coef_sdsp
     1               + e ( ixx) / cfl_coef_advct)

          dtcurp1 = AMIN1 ( dtcurp1, d ( ixx))

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE EVALUATE_TIMESTEP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INCNAME_8 ( filname, inc)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL MYSTOP

      CHARACTER*8 filname

      CHARACTER*1 unitchar ( 8)

      INTEGER unitint ( 8), zeroint, value, inc, i, j

! ----------------------------------------------------------------------

      unitchar ( 1) = filname ( 1:1)
      unitchar ( 2) = filname ( 2:2)
      unitchar ( 3) = filname ( 3:3)
      unitchar ( 4) = filname ( 4:4)
      unitchar ( 5) = filname ( 5:5)
      unitchar ( 6) = filname ( 6:6)
      unitchar ( 7) = filname ( 7:7)
      unitchar ( 8) = filname ( 8:8)

      zeroint = ichar ( '0')

      unitint ( 8) = ichar ( unitchar ( 8) ) - zeroint
      unitint ( 7) = ichar ( unitchar ( 7) ) - zeroint
      unitint ( 6) = ichar ( unitchar ( 6) ) - zeroint
      unitint ( 5) = ichar ( unitchar ( 5) ) - zeroint

      value = 1000 * unitint ( 5) + 100 * unitint ( 6)
     1              + 10 * unitint ( 7) + unitint ( 8)

      value = value + inc


      unitint ( 5) = value / 1000

      value = value - 1000 * unitint ( 5)

      unitint ( 6) = value / 100

      value = value - 100 * unitint ( 6)

      unitint ( 7) = value / 10

      value = value - 10 * unitint ( 7)

      unitint ( 8) = value

      unitint ( 5) = unitint ( 5) + zeroint
      unitint ( 6) = unitint ( 6) + zeroint
      unitint ( 7) = unitint ( 7) + zeroint
      unitint ( 8) = unitint ( 8) + zeroint

      unitchar ( 5) = char ( unitint ( 5) )
      unitchar ( 6) = char ( unitint ( 6) )
      unitchar ( 7) = char ( unitint ( 7) )
      unitchar ( 8) = char ( unitint ( 8) )

      filname ( 5:5) = unitchar ( 5)
      filname ( 6:6) = unitchar ( 6)
      filname ( 7:7) = unitchar ( 7)
      filname ( 8:8) = unitchar ( 8)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INCNAME_8

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INCNAME_SP ( filname, inc)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL MYSTOP

      CHARACTER*8 filname

      CHARACTER*1 unitchar ( 8)

      INTEGER unitint ( 8), zeroint, value, inc, i, j

! ----------------------------------------------------------------------

      unitchar ( 1) = filname ( 1:1)
      unitchar ( 2) = filname ( 2:2)
      unitchar ( 3) = filname ( 3:3)
      unitchar ( 4) = filname ( 4:4)
      unitchar ( 5) = filname ( 5:5)
      unitchar ( 6) = filname ( 6:6)
      unitchar ( 7) = filname ( 7:7)
      unitchar ( 8) = filname ( 8:8)

      zeroint = ichar ( '0')

      unitint ( 4) = ichar ( unitchar ( 4) ) - zeroint
      unitint ( 3) = ichar ( unitchar ( 3) ) - zeroint

      value = 10 * unitint ( 3) + unitint ( 4)

      value = value + inc

      unitint ( 3) = value / 10

      value = value - 10 * unitint ( 3)

      unitint ( 4) = value

      unitint ( 3) = unitint ( 3) + zeroint
      unitint ( 4) = unitint ( 4) + zeroint

      unitchar ( 3) = char ( unitint ( 3) )
      unitchar ( 4) = char ( unitint ( 4) )

      filname ( 3:3) = unitchar ( 3)
      filname ( 4:4) = unitchar ( 4)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INCNAME_SP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MESHGEN ( level, base_grid_level, finest_grid_level,
     1      ioffx, ioffy, ioffz, smallnum,

     1      base_grid_nx, base_grid_ny, base_grid_nz,
     1      base_grid_firsthalf_nx, base_grid_firsthalf_ny,
     1      base_grid_firsthalf_nz,

     1      base_grid_xmin, base_grid_xmax,
     1      base_grid_ymin, base_grid_ymax,
     1      base_grid_zmin, base_grid_zmax,
     1      base_grid_xratio, base_grid_yratio, base_grid_zratio,

     1      my_level_ixmin, my_level_ixmax,
     1      my_level_iymin, my_level_iymax,
     1      my_level_izmin, my_level_izmax,

     1      my_xmin, my_xmax, my_dx0, my_xratio,
     1      my_cut_ix0, my_cut_xb0,

     1      my_ymin, my_ymax, my_dy0, my_yratio,
     1      my_cut_iy0, my_cut_yb0,

     1      my_zmin, my_zmax, my_dz0, my_zratio,
     1      my_cut_iz0, my_cut_zb0)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! Takes the beginning and ending indices of a patch in a global index space
! ( i.e. "my_level_ixmin, my_level_ixmax" etc., at a given "level" ) and
! builds a few things about that patch:

! a) It builds the physical boundaries of that patch, "my_xmin, my_xmax" etc.

! b) It builds the zone size that one should start with on that patch,
!    i.e. "my_dx0" etc.

! c) It builds the zone ratioing at that level for the patch,
!    i.e. "my_xratio" etc.

! d) If the ratioed mesh is to be built in two halves then "my_cut_ix0" tell us
!    where to put the cut. If this number is less than 0 then don't worry
!    about it.

! e) This tells us what value to assign to the cut boundary. The logic is:
!    "this%xsav ( my_cut_ix0 - 1) = my_cut_xb0".

! The variables we key off on are: "base_grid_nx,y,z", "base_grid_x,y,zratio",
! "base_grid_x,y,zmin", "base_grid_x,y,zmax", all of which are specified
! on the base level mesh.

! It works with uniform meshes and geometrically ratioed meshes.
! For any other mesh, pl. provide your own version of "MESHGEN".

! Uniform meshes are produced by setting "base_grid_xratio == 1.0", ||larly
! for y and z. Setting "base_grid_xratio != 1.0" triggers ratioed meshes.

! Setting "base_grid_firsthalf_nx <= 0" etc. triggers ratioed meshes that
! increase in one direction.

! When "base_grid_firsthalf_nx > 0", we are telling the code that the first
! "base_grid_firsthalf_nx" zones on the base grid have ratioing that
! increases to the left while the next "base_grid_nx-base_grid_firsthalf_nx"
! on the base grid have ratioing that increases to the right. Thus the
! ratioing is symmetrical about the right boundary of zone
! "base_grid_firsthalf_nx".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER level, base_grid_level, finest_grid_level,
     1      ioffx, ioffy, ioffz,

     1      base_grid_nx, base_grid_ny, base_grid_nz,
     1      base_grid_firsthalf_nx, base_grid_firsthalf_ny,
     1      base_grid_firsthalf_nz,

     1      my_level_ixmin, my_level_ixmax,
     1      my_level_iymin, my_level_iymax,
     1      my_level_izmin, my_level_izmax,

     1      my_cut_ix0, my_cut_iy0, my_cut_iz0

      REAL smallnum,
     1     base_grid_xmin, base_grid_xmax,
     1     base_grid_ymin, base_grid_ymax,
     1     base_grid_zmin, base_grid_zmax,
     1     base_grid_xratio, base_grid_yratio, base_grid_zratio,

     1     my_cut_xb0, my_cut_yb0, my_cut_zb0,

     1     my_xmin, my_xmax, my_dx0, my_xratio,
     1     my_ymin, my_ymax, my_dy0, my_yratio,
     1     my_zmin, my_zmax, my_dz0, my_zratio

! ----------------------------------------------------------------------

      INTEGER i, level_nx, level_ny, level_nz,
     1      level_firsthalf_nx, level_firsthalf_ny, level_firsthalf_nz

      REAL tempaa, tempab

! ----------------------------------------------------------------------

      my_xratio = base_grid_xratio
      my_yratio = base_grid_yratio
      my_zratio = base_grid_zratio

      level_nx = base_grid_nx
      level_ny = base_grid_ny
      level_nz = base_grid_nz

      level_firsthalf_nx = base_grid_firsthalf_nx
      level_firsthalf_ny = base_grid_firsthalf_ny
      level_firsthalf_nz = base_grid_firsthalf_nz

      my_cut_ix0 = - base_grid_nx
      my_cut_iy0 = - base_grid_ny
      my_cut_iz0 = - base_grid_nz

      my_cut_xb0 = 0.0
      my_cut_yb0 = 0.0
      my_cut_zb0 = 0.0

! ------------------------------

! For "level > base_grid_level" we do the following:

      DO i = base_grid_level + 1, level

      IF ( ioffx > 0) THEN
        my_xratio = SQRT ( my_xratio)
        level_nx = level_nx * 2
        level_firsthalf_nx = level_firsthalf_nx * 2
      END IF

      IF ( ioffy > 0) THEN
        my_yratio = SQRT ( my_yratio)
        level_ny = level_ny * 2
        level_firsthalf_ny = level_firsthalf_ny * 2
      END IF

      IF ( ioffz > 0) THEN
        my_zratio = SQRT ( my_zratio)
        level_nz = level_nz * 2
        level_firsthalf_nz = level_firsthalf_nz * 2
      END IF

      END DO

! ------------------------------

! For "level < base_grid_level" we do the following:

      DO i = base_grid_level - 1, level, - 1

      IF ( ioffx > 0) THEN
        my_xratio = my_xratio**2
        level_nx = level_nx / 2
        level_firsthalf_nx = level_firsthalf_nx / 2
      END IF

      IF ( ioffy > 0) THEN
        my_yratio = my_yratio**2
        level_ny = level_ny / 2
        level_firsthalf_ny = level_firsthalf_ny / 2
      END IF

      IF ( ioffz > 0) THEN
        my_zratio = my_zratio**2
        level_nz = level_nz / 2
        level_firsthalf_nz = level_firsthalf_nz / 2
      END IF

      END DO

! ----------------------------------------------------------------------

! Check to make sure that "my_level_ixmin, my_level_ixmax" etc. are within
! the limits

! Ensure that starting and ending indices are well-bounded.

      IF ( ( my_level_ixmin .LT. 1) .OR.
     1     ( my_level_ixmin .GT. level_nx) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"level_nx, my_level_ixmin=", level_nx, my_level_ixmin
      STOP
      END IF

      IF ( ( my_level_ixmax .LT. 1) .OR.
     1     ( my_level_ixmax .GT. level_nx) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"level_nx, my_level_ixmax=", level_nx, my_level_ixmax
      STOP
      END IF

      IF ( my_level_ixmin .GT. my_level_ixmax) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"my_level_ixmin, my_level_ixmax = ",
     1              my_level_ixmin, my_level_ixmax
      STOP
      END IF

! Only active directions may have ratioed meshes.

      IF ( ( ioffx <= 0) .AND.
     1     ( ABS ( my_xratio - 1.0) .GT. smallnum) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"x-direction inactive but my_xratio = ", my_xratio
      STOP
      END IF

! If the mesh is ratioed, it should either have "base_grid_firsthalf_nx <= 0"
! (which produces ratioing in only one direction),
! or it should have lower and upper bounds that straddle zero.
! When a ratioed mesh straddles zero, it should satisfy some strict
! conditions associated with the zones being of the same size across zero.

      IF ( ABS ( my_xratio - 1.0) .GT. smallnum) THEN

        IF ( base_grid_firsthalf_nx .GT. 0) THEN

          IF ( ( level_firsthalf_nx .LT. 1) .OR.
     1         ( level_firsthalf_nx .GT. level_nx) ) THEN
          WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
          WRITE ( 6, *)"level_firsthalf_nx, level_nx = ",
     1                  level_firsthalf_nx, level_nx
          STOP
          END IF

        END IF

      END IF

! ------------------------------

      IF ( ( my_level_iymin .LT. 1) .OR.
     1     ( my_level_iymin .GT. level_ny) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"level_ny, my_level_iymin=", level_ny, my_level_iymin
      STOP
      END IF

      IF ( ( my_level_iymax .LT. 1) .OR.
     1     ( my_level_iymax .GT. level_ny) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"level_ny, my_level_iymax=", level_ny, my_level_iymax
      STOP
      END IF

      IF ( my_level_iymin .GT. my_level_iymax) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"my_level_iymin, my_level_iymax = ",
     1              my_level_iymin, my_level_iymax
      STOP
      END IF

      IF ( ( ioffy <= 0) .AND.
     1     ( ABS ( my_yratio - 1.0) .GT. smallnum) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"y-direction inactive but my_yratio = ", my_yratio
      STOP
      END IF

      IF ( ABS ( my_yratio - 1.0) .GT. smallnum) THEN

        IF ( base_grid_firsthalf_ny .GT. 0) THEN

          IF ( ( level_firsthalf_ny .LT. 1) .OR.
     1         ( level_firsthalf_ny .GT. level_ny) ) THEN
          WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
          WRITE ( 6, *)"level_firsthalf_ny, level_ny = ",
     1                  level_firsthalf_ny, level_ny
          STOP
          END IF

        END IF

      END IF

! ------------------------------

      IF ( ( my_level_izmin .LT. 1) .OR.
     1     ( my_level_izmin .GT. level_nz) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"level_nz, my_level_izmin=", level_nz, my_level_izmin
      STOP
      END IF

      IF ( ( my_level_izmax .LT. 1) .OR.
     1     ( my_level_izmax .GT. level_nz) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"level_nz, my_level_izmax=", level_nz, my_level_izmax
      STOP
      END IF

      IF ( my_level_izmin .GT. my_level_izmax) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"my_level_izmin, my_level_izmax = ",
     1              my_level_izmin, my_level_izmax
      STOP
      END IF

      IF ( ( ioffz <= 0) .AND.
     1     ( ABS ( my_zratio - 1.0) .GT. smallnum) ) THEN
      WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
      WRITE ( 6, *)"z-direction inactive but my_zratio = ", my_zratio
      STOP
      END IF

      IF ( ABS ( my_zratio - 1.0) .GT. smallnum) THEN

        IF ( base_grid_firsthalf_nz .GT. 0) THEN

          IF ( ( level_firsthalf_nz .LT. 1) .OR.
     1         ( level_firsthalf_nz .GT. level_nz) ) THEN
          WRITE ( 6, *)"Stopping in MESHGEN, level = ", level
          WRITE ( 6, *)"level_firsthalf_nz, level_nz = ",
     1                  level_firsthalf_nz, level_nz
          STOP
          END IF

        END IF

      END IF

! ----------------------------------------------------------------------

      IF ( ABS ( my_xratio - 1.0) .LE. smallnum) THEN

! Uniform meshes along x.

! ------------------------------

        my_dx0 = ( base_grid_xmax - base_grid_xmin) / level_nx
        my_xmin = base_grid_xmin + ( my_level_ixmin - 1.0) * my_dx0
        my_xmax = base_grid_xmin + my_level_ixmax * my_dx0

! ------------------------------

      ELSE IF ( base_grid_firsthalf_nx .LE. 0) THEN

! Ratioed meshes in one direction along x.

! ------------------------------

        my_dx0 = ( base_grid_xmax - base_grid_xmin)
     1         * ( my_xratio - 1.0) / ( my_xratio**level_nx - 1.0)

        my_xmin = base_grid_xmin
     1    + my_dx0 * ( my_xratio**my_level_ixmin / my_xratio - 1.0)
     1             / ( my_xratio - 1.0)

        my_xmax = base_grid_xmin
     1    + my_dx0 * ( my_xratio**my_level_ixmax - 1.0)
     1             / ( my_xratio - 1.0)

        my_dx0 = my_dx0 * my_xratio**my_level_ixmin / my_xratio

! ------------------------------

      ELSE

! Symmetrically ratioed meshes along x.

! ------------------------------

        tempaa = my_xratio**level_firsthalf_nx - 1.0
     1         + my_xratio**(level_nx - level_firsthalf_nx) - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)

        my_dx0 = ( base_grid_xmax - base_grid_xmin) / tempaa

        tempaa = my_xratio**level_firsthalf_nx - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_cut_xb0 = base_grid_xmin + my_dx0 * tempaa

! ----------

        IF ( my_level_ixmin .GT. level_firsthalf_nx) THEN

        tempaa = my_xratio**(my_level_ixmin - level_firsthalf_nx)
     1         / my_xratio - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_xmin = my_cut_xb0 + my_dx0 * tempaa

        tempaa = my_xratio**(my_level_ixmax - level_firsthalf_nx) - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_xmax = my_cut_xb0 + my_dx0 * tempaa

        my_dx0 = my_dx0 * my_xratio**(my_level_ixmin
     1                               - level_firsthalf_nx) / my_xratio

! ----------

        ELSE IF ( ( my_level_ixmin .LE. level_firsthalf_nx) .AND.
     1            ( my_level_ixmax .GE. level_firsthalf_nx)) THEN

        my_cut_ix0 = level_firsthalf_nx - my_level_ixmin + 2

        tempaa = my_xratio**( level_firsthalf_nx - my_level_ixmin + 1)
     1         - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_xmin = my_cut_xb0 - my_dx0 * tempaa

        tempaa = my_xratio**( my_level_ixmax - level_firsthalf_nx + 1)
     1         / my_xratio - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_xmax = my_cut_xb0 + my_dx0 * tempaa

! ----------

        ELSE

        my_cut_ix0 = my_level_ixmax - my_level_ixmin + 2

        tempaa = my_xratio**( level_firsthalf_nx - my_level_ixmin + 1)
     1         - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_xmin = my_cut_xb0 - my_dx0 * tempaa

        tempaa = my_xratio**( level_firsthalf_nx - my_level_ixmax) - 1.0
        tempaa = tempaa / ( my_xratio - 1.0)
        my_xmax = my_cut_xb0 - my_dx0 * tempaa

        my_cut_xb0 = my_xmax

        my_dx0 = my_dx0 * my_xratio**( level_firsthalf_nx
     1                               - my_level_ixmax)

        END IF

! ------------------------------

      END IF

! ----------------------------------------------------------------------

      IF ( ABS ( my_yratio - 1.0) .LE. smallnum) THEN

! Uniform meshes along y.

! ------------------------------

        my_dy0 = ( base_grid_ymax - base_grid_ymin) / level_ny
        my_ymin = base_grid_ymin + ( my_level_iymin - 1.0) * my_dy0
        my_ymax = base_grid_ymin + my_level_iymax * my_dy0

! ------------------------------

      ELSE IF ( base_grid_firsthalf_ny .LE. 0) THEN

! Ratioed meshes in one direction along y.

! ------------------------------

        my_dy0 = ( base_grid_ymax - base_grid_ymin)
     1         * ( my_yratio - 1.0) / ( my_yratio**level_ny - 1.0)

        my_ymin = base_grid_ymin
     1    + my_dy0 * ( my_yratio**my_level_iymin / my_yratio - 1.0)
     1             / ( my_yratio - 1.0)

        my_ymax = base_grid_ymin
     1    + my_dy0 * ( my_yratio**my_level_iymax - 1.0)
     1             / ( my_yratio - 1.0)

        my_dy0 = my_dy0 * my_yratio**my_level_iymin / my_yratio

! ------------------------------

      ELSE

! Symmetrically ratioed meshes along y.

! ------------------------------

        tempaa = my_yratio**level_firsthalf_ny - 1.0
     1         + my_yratio**(level_ny - level_firsthalf_ny) - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)

        my_dy0 = ( base_grid_ymax - base_grid_ymin) / tempaa

        tempaa = my_yratio**level_firsthalf_ny - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_cut_yb0 = base_grid_ymin + my_dy0 * tempaa

! ----------

        IF ( my_level_iymin .GT. level_firsthalf_ny) THEN

        tempaa = my_yratio**(my_level_iymin - level_firsthalf_ny)
     1         / my_yratio - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_ymin = my_cut_yb0 + my_dy0 * tempaa

        tempaa = my_yratio**(my_level_iymax - level_firsthalf_ny) - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_ymax = my_cut_yb0 + my_dy0 * tempaa

        my_dy0 = my_dy0 * my_yratio**(my_level_iymin
     1                               - level_firsthalf_ny) / my_yratio

! ----------

        ELSE IF ( ( my_level_iymin .LE. level_firsthalf_ny) .AND.
     1            ( my_level_iymax .GE. level_firsthalf_ny)) THEN

        my_cut_iy0 = level_firsthalf_ny - my_level_iymin + 2

        tempaa = my_yratio**( level_firsthalf_ny - my_level_iymin + 1)
     1         - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_ymin = my_cut_yb0 - my_dy0 * tempaa

        tempaa = my_yratio**( my_level_iymax - level_firsthalf_ny + 1)
     1         / my_yratio - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_ymax = my_cut_yb0 + my_dy0 * tempaa

! ----------

        ELSE

        my_cut_iy0 = my_level_iymax - my_level_iymin + 2

        tempaa = my_yratio**( level_firsthalf_ny - my_level_iymin + 1)
     1         - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_ymin = my_cut_yb0 - my_dy0 * tempaa

        tempaa = my_yratio**( level_firsthalf_ny - my_level_iymax) - 1.0
        tempaa = tempaa / ( my_yratio - 1.0)
        my_ymax = my_cut_yb0 - my_dy0 * tempaa

        my_cut_yb0 = my_ymax

        my_dy0 = my_dy0 * my_yratio**( level_firsthalf_ny
     1                               - my_level_iymax)

        END IF

! ------------------------------

      END IF

! ----------------------------------------------------------------------

      IF ( ABS ( my_zratio - 1.0) .LE. smallnum) THEN

! Uniform meshes along z.

! ------------------------------

        my_dz0 = ( base_grid_zmax - base_grid_zmin) / level_nz
        my_zmin = base_grid_zmin + ( my_level_izmin - 1.0) * my_dz0
        my_zmax = base_grid_zmin + my_level_izmax * my_dz0

! ------------------------------

      ELSE IF ( base_grid_firsthalf_nz .LE. 0) THEN

! Ratioed meshes in one direction along z.

! ------------------------------

        my_dz0 = ( base_grid_zmax - base_grid_zmin)
     1         * ( my_zratio - 1.0) / ( my_zratio**level_nz - 1.0)

        my_zmin = base_grid_zmin
     1    + my_dz0 * ( my_zratio**my_level_izmin / my_zratio - 1.0)
     1             / ( my_zratio - 1.0)

        my_zmax = base_grid_zmin
     1    + my_dz0 * ( my_zratio**my_level_izmax - 1.0)
     1             / ( my_zratio - 1.0)

        my_dz0 = my_dz0 * my_zratio**my_level_izmin / my_zratio

! ------------------------------

      ELSE

! Symmetrically ratioed meshes along z.

! ------------------------------

        tempaa = my_zratio**level_firsthalf_nz - 1.0
     1         + my_zratio**(level_nz - level_firsthalf_nz) - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)

        my_dz0 = ( base_grid_zmax - base_grid_zmin) / tempaa

        tempaa = my_zratio**level_firsthalf_nz - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_cut_zb0 = base_grid_zmin + my_dz0 * tempaa

! ----------

        IF ( my_level_izmin .GT. level_firsthalf_nz) THEN

        tempaa = my_zratio**(my_level_izmin - level_firsthalf_nz)
     1         / my_zratio - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_zmin = my_cut_zb0 + my_dz0 * tempaa

        tempaa = my_zratio**(my_level_izmax - level_firsthalf_nz) - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_zmax = my_cut_zb0 + my_dz0 * tempaa

        my_dz0 = my_dz0 * my_zratio**(my_level_izmin
     1                               - level_firsthalf_nz) / my_zratio

! ----------

        ELSE IF ( ( my_level_izmin .LE. level_firsthalf_nz) .AND.
     1            ( my_level_izmax .GE. level_firsthalf_nz)) THEN

        my_cut_iz0 = level_firsthalf_nz - my_level_izmin + 2

        tempaa = my_zratio**( level_firsthalf_nz - my_level_izmin + 1)
     1         - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_zmin = my_cut_zb0 - my_dz0 * tempaa

        tempaa = my_zratio**( my_level_izmax - level_firsthalf_nz + 1)
     1         / my_zratio - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_zmax = my_cut_zb0 + my_dz0 * tempaa

! ----------

        ELSE

        my_cut_iz0 = my_level_izmax - my_level_izmin + 2

        tempaa = my_zratio**( level_firsthalf_nz - my_level_izmin + 1)
     1         - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_zmin = my_cut_zb0 - my_dz0 * tempaa

        tempaa = my_zratio**( level_firsthalf_nz - my_level_izmax) - 1.0
        tempaa = tempaa / ( my_zratio - 1.0)
        my_zmax = my_cut_zb0 - my_dz0 * tempaa

        my_cut_zb0 = my_zmax

        my_dz0 = my_dz0 * my_zratio**( level_firsthalf_nz
     1                               - my_level_izmax)

        END IF

! ------------------------------

      END IF

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MESHGEN

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Start of debugging test code for MESHGEN.
#if ( 0 == 1)
      PROGRAM TEST_MESHGEN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER level, base_grid_level, finest_grid_level,
     1      ioffx, ioffy, ioffz,

     1      base_grid_nx, base_grid_ny, base_grid_nz,
     1      base_grid_firsthalf_nx, base_grid_firsthalf_ny,
     1      base_grid_firsthalf_nz,

     1      my_level_ixmin, my_level_ixmax,
     1      my_level_iymin, my_level_iymax,
     1      my_level_izmin, my_level_izmax,

     1      my_cut_ix0, my_cut_iy0, my_cut_iz0

      REAL smallnum,
     1     base_grid_xmin, base_grid_xmax,
     1     base_grid_ymin, base_grid_ymax,
     1     base_grid_zmin, base_grid_zmax,
     1     base_grid_xratio, base_grid_yratio, base_grid_zratio,

     1     my_cut_xb0, my_cut_yb0, my_cut_zb0,

     1     my_xmin, my_xmax, my_dx0, my_xratio,
     1     my_ymin, my_ymax, my_dy0, my_yratio,
     1     my_zmin, my_zmax, my_dz0, my_zratio

! ----------------------------------------------------------------------

      INTEGER ix1, iy1, iz1, ixx, iyy, izz,

     1        test_level_ixmin, test_level_ixmax,
     1        test_level_iymin, test_level_iymax,
     1        test_level_izmin, test_level_izmax

      PARAMETER ( ix1 = 10, iy1 = 16, iz1 = 20,
     1            ioffx = 4, ioffy = 4, ioffz = 4,

     1            test_level_ixmin = 1, test_level_ixmax = 40,
     1            test_level_iymin = 1, test_level_iymax = 64,
     1            test_level_izmin = 1, test_level_izmax = 80)

! ----------

! Arrays that specify global geometry.

      REAL global_xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL global_dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL global_xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL global_yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL global_dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL global_yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL global_zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL global_dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL global_zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------

! Arrays that specify local geometry.

      REAL my_xb_save
     1        ( test_level_ixmin - 1 - ioffx: test_level_ixmax + ioffx)
      REAL my_dx_save
     1            ( test_level_ixmin - ioffx: test_level_ixmax + ioffx)
      REAL my_xc_save
     1            ( test_level_ixmin - ioffx: test_level_ixmax + ioffx)

      REAL my_yb_save
     1        ( test_level_iymin - 1 - ioffy: test_level_iymax + ioffy)
      REAL my_dy_save
     1            ( test_level_iymin - ioffy: test_level_iymax + ioffy)
      REAL my_yc_save
     1            ( test_level_iymin - ioffy: test_level_iymax + ioffy)

      REAL my_zb_save
     1        ( test_level_izmin - 1 - ioffz: test_level_izmax + ioffz)
      REAL my_dz_save
     1            ( test_level_izmin - ioffz: test_level_izmax + ioffz)
      REAL my_zc_save
     1            ( test_level_izmin - ioffz: test_level_izmax + ioffz)


! ----------------------------------------------------------------------

      smallnum = 1.0e-12

      base_grid_level = 5
      finest_grid_level = 8

      OPEN ( unit = 6, file = "test_meshgen.out")

! ----------------------------------------------------------------------
! Set up information about the base level grid, i.e. the global grid.
! Then build the base level mesh.

! Number of global zones in each direction on the base grid.
      base_grid_nx = ix1
      base_grid_ny = iy1
      base_grid_nz = iz1

! If we want the base level grid to be symmetrically ratioed about the origin,
! then specify that here. Else, set this to a negative number in which case
! the zones will only be increasing in ratio as the mesh's index increases.
      base_grid_firsthalf_nx = ix1 / 2
      base_grid_firsthalf_ny = iy1 / 2
      base_grid_firsthalf_nz = iz1 / 2

! Set up the physical limits of the base mesh; min and max in each direction.
      base_grid_xmin = - 1.5
      base_grid_xmax = 1.5
      base_grid_ymin = - 2.0
      base_grid_ymax = 2.0
      base_grid_zmin = - 1.0
      base_grid_zmax = 1.0

! Ratio by which grid zones are inflated in each direction.
      base_grid_xratio = 1.02
      base_grid_yratio = 1.04
      base_grid_zratio = 1.06

! ------------------------------
! Since we want to build the base level mesh, specify the integer limits to
! be the global limits and the level to the base level.

      level = base_grid_level

      my_level_ixmin = 1
      my_level_ixmax = ix1
      my_level_iymin = 1
      my_level_iymax = iy1
      my_level_izmin = 1
      my_level_izmax = iz1

! ------------------------------
! Call MESHGEN for the base level mesh.

      CALL MESHGEN ( level, base_grid_level, finest_grid_level,
     1      ioffx, ioffy, ioffz, smallnum,

     1      base_grid_nx, base_grid_ny, base_grid_nz,
     1      base_grid_firsthalf_nx, base_grid_firsthalf_ny,
     1      base_grid_firsthalf_nz,

     1      base_grid_xmin, base_grid_xmax,
     1      base_grid_ymin, base_grid_ymax,
     1      base_grid_zmin, base_grid_zmax,
     1      base_grid_xratio, base_grid_yratio, base_grid_zratio,

     1      my_level_ixmin, my_level_ixmax,
     1      my_level_iymin, my_level_iymax,
     1      my_level_izmin, my_level_izmax,

     1      my_xmin, my_xmax, my_dx0, my_xratio,
     1      my_cut_ix0, my_cut_xb0,

     1      my_ymin, my_ymax, my_dy0, my_yratio,
     1      my_cut_iy0, my_cut_yb0,

     1      my_zmin, my_zmax, my_dz0, my_zratio,
     1      my_cut_iz0, my_cut_zb0)

! ------------------------------
! Now build the base level mesh.

! ------------------------------
! Build mesh in the x-direction.

      IF ( my_cut_ix0 .LE. 0) THEN

! -----

        global_xb_save ( 0) = my_xmin

! -----

        DO ixx = 1, my_level_ixmax + ioffx

        global_dx_save ( ixx) = my_dx0 * my_xratio**ixx / my_xratio
        global_xb_save ( ixx) = global_xb_save ( ixx - 1)
     1                        + global_dx_save ( ixx)
        global_xc_save ( ixx) = 0.5 * ( global_xb_save ( ixx)
     1                                + global_xb_save ( ixx - 1) )

        END DO

! -----

        DO ixx = 0, my_level_ixmin - ioffx, -1

        global_dx_save ( ixx) = global_dx_save ( my_level_ixmin - ixx)
        global_xb_save ( ixx - 1) = global_xb_save ( ixx)
     1                            - global_dx_save ( ixx)
        global_xc_save ( ixx) = 0.5 * ( global_xb_save ( ixx)
     1                                + global_xb_save ( ixx - 1) )

        END DO

! -----

      ELSE

! -----

        global_xb_save ( my_cut_ix0 - 1) = my_cut_xb0

! -----

        DO ixx = my_cut_ix0, my_level_ixmax + ioffx

        global_dx_save ( ixx) = my_dx0
     1          * my_xratio**( ixx - my_cut_ix0 + 1) / my_xratio
        global_xb_save ( ixx) = global_xb_save ( ixx - 1)
     1                        + global_dx_save ( ixx)
        global_xc_save ( ixx) = 0.5 * ( global_xb_save ( ixx)
     1                                + global_xb_save ( ixx - 1) )

        END DO

! -----

        DO ixx = my_cut_ix0 - 1, my_level_ixmin - ioffx, -1

        global_dx_save ( ixx) = my_dx0
     1          * my_xratio**( my_cut_ix0 - ixx) / my_xratio
        global_xb_save ( ixx - 1) = global_xb_save ( ixx)
     1                            - global_dx_save ( ixx)
        global_xc_save ( ixx) = 0.5 * ( global_xb_save ( ixx)
     1                                + global_xb_save ( ixx - 1) )

        END DO

! -----

      END IF

! ------------------------------
! Build mesh in the y-direction.

      IF ( my_cut_iy0 .LE. 0) THEN

! -----

        global_yb_save ( 0) = my_ymin

! -----

        DO iyy = 1, my_level_iymax + ioffy

        global_dy_save ( iyy) = my_dy0 * my_yratio**iyy / my_yratio
        global_yb_save ( iyy) = global_yb_save ( iyy - 1)
     1                        + global_dy_save ( iyy)
        global_yc_save ( iyy) = 0.5 * ( global_yb_save ( iyy)
     1                                + global_yb_save ( iyy - 1) )

        END DO

! -----

        DO iyy = 0, my_level_iymin - ioffy, -1

        global_dy_save ( iyy) = global_dy_save ( my_level_iymin - iyy)
        global_yb_save ( iyy - 1) = global_yb_save ( iyy)
     1                            - global_dy_save ( iyy)
        global_yc_save ( iyy) = 0.5 * ( global_yb_save ( iyy)
     1                                + global_yb_save ( iyy - 1) )

        END DO

! -----

      ELSE

! -----

        global_yb_save ( my_cut_iy0 - 1) = my_cut_yb0

! -----

        DO iyy = my_cut_iy0, my_level_iymax + ioffy

        global_dy_save ( iyy) = my_dy0
     1          * my_yratio**( iyy - my_cut_iy0 + 1) / my_yratio
        global_yb_save ( iyy) = global_yb_save ( iyy - 1)
     1                        + global_dy_save ( iyy)
        global_yc_save ( iyy) = 0.5 * ( global_yb_save ( iyy)
     1                                + global_yb_save ( iyy - 1) )

        END DO

! -----

        DO iyy = my_cut_iy0 - 1, my_level_iymin - ioffy, -1

        global_dy_save ( iyy) = my_dy0
     1          * my_yratio**( my_cut_iy0 - iyy) / my_yratio
        global_yb_save ( iyy - 1) = global_yb_save ( iyy)
     1                            - global_dy_save ( iyy)
        global_yc_save ( iyy) = 0.5 * ( global_yb_save ( iyy)
     1                                + global_yb_save ( iyy - 1) )

        END DO

! -----

      END IF

! ------------------------------
! Build mesh in the z-direction.

      IF ( my_cut_iz0 .LE. 0) THEN

! -----

        global_zb_save ( 0) = my_zmin

! -----

        DO izz = 1, my_level_izmax + ioffz

        global_dz_save ( izz) = my_dz0 * my_zratio**izz / my_zratio
        global_zb_save ( izz) = global_zb_save ( izz - 1)
     1                        + global_dz_save ( izz)
        global_zc_save ( izz) = 0.5 * ( global_zb_save ( izz)
     1                                + global_zb_save ( izz - 1) )

        END DO

! -----

        DO izz = 0, my_level_izmin - ioffz, -1

        global_dz_save ( izz) = global_dz_save ( my_level_izmin - izz)
        global_zb_save ( izz - 1) = global_zb_save ( izz)
     1                            - global_dz_save ( izz)
        global_zc_save ( izz) = 0.5 * ( global_zb_save ( izz)
     1                                + global_zb_save ( izz - 1) )

        END DO

! -----

      ELSE

! -----

        global_zb_save ( my_cut_iz0 - 1) = my_cut_zb0

! -----

        DO izz = my_cut_iz0, my_level_izmax + ioffz

        global_dz_save ( izz) = my_dz0
     1          * my_zratio**( izz - my_cut_iz0 + 1) / my_zratio
        global_zb_save ( izz) = global_zb_save ( izz - 1)
     1                        + global_dz_save ( izz)
        global_zc_save ( izz) = 0.5 * ( global_zb_save ( izz)
     1                                + global_zb_save ( izz - 1) )

        END DO

! -----

        DO izz = my_cut_iz0 - 1, my_level_izmin - ioffz, -1

        global_dz_save ( izz) = my_dz0
     1          * my_zratio**( my_cut_iz0 - izz) / my_zratio
        global_zb_save ( izz - 1) = global_zb_save ( izz)
     1                            - global_dz_save ( izz)
        global_zc_save ( izz) = 0.5 * ( global_zb_save ( izz)
     1                                + global_zb_save ( izz - 1) )

        END DO

! -----

      END IF

! ------------------------------
! Write out the information about the base level mesh in the x-direction.

      WRITE ( 6, *)" ---------------------------------"
      WRITE ( 6, *)" base level mesh in the x-direction"
      WRITE ( 6, *)" ---------------------------------"
      ixx = my_level_ixmin - ioffx - 1
      WRITE ( 6, *)" ixx = ", ixx
      WRITE ( 6, *)"global_xb_save = ", global_xb_save ( ixx)
      WRITE ( 6, *)" "

      DO ixx = my_level_ixmin - ioffx, my_level_ixmax + ioffx

      WRITE ( 6, *)" ixx = ", ixx
      WRITE ( 6, *)"global_dx_save = ", global_dx_save ( ixx),
     1             " global_xc_save = ", global_xc_save ( ixx)
      WRITE ( 6, *)"global_xb_save = ", global_xb_save ( ixx)
      WRITE ( 6, *)" "

      END DO

! ------------------------------
! Write out the information about the base level mesh in the y-direction.

      WRITE ( 6, *)" ---------------------------------"
      WRITE ( 6, *)" base level mesh in the y-direction"
      WRITE ( 6, *)" ---------------------------------"
      iyy = my_level_iymin - ioffy - 1
      WRITE ( 6, *)" iyy = ", iyy
      WRITE ( 6, *)"global_yb_save = ", global_yb_save ( iyy)
      WRITE ( 6, *)" "

      DO iyy = my_level_iymin - ioffy, my_level_iymax + ioffy

      WRITE ( 6, *)" iyy = ", iyy
      WRITE ( 6, *)"global_dy_save = ", global_dy_save ( iyy),
     1             " global_yc_save = ", global_yc_save ( iyy)
      WRITE ( 6, *)"global_yb_save = ", global_yb_save ( iyy)
      WRITE ( 6, *)" "

      END DO

! ------------------------------
! Write out the information about the base level mesh in the z-direction.

      WRITE ( 6, *)" ---------------------------------"
      WRITE ( 6, *)" base level mesh in the z-direction"
      WRITE ( 6, *)" ---------------------------------"
      izz = my_level_izmin - ioffz - 1
      WRITE ( 6, *)" izz = ", izz
      WRITE ( 6, *)"global_zb_save = ", global_zb_save ( izz)
      WRITE ( 6, *)" "

      DO izz = my_level_izmin - ioffz, my_level_izmax + ioffz

      WRITE ( 6, *)" izz = ", izz
      WRITE ( 6, *)"global_dz_save = ", global_dz_save ( izz),
     1             " global_zc_save = ", global_zc_save ( izz)
      WRITE ( 6, *)"global_zb_save = ", global_zb_save ( izz)
      WRITE ( 6, *)" "

      END DO

! ----------------------------------------------------------------------
! Since we want to build the finer/coarser level mesh, specify the integer
! limits to be the global limits and the level to the base level.

      level = base_grid_level + 2

      my_level_ixmin = test_level_ixmin
      my_level_ixmax = test_level_ixmax
      my_level_iymin = test_level_iymin
      my_level_iymax = test_level_iymax
      my_level_izmin = test_level_izmin
      my_level_izmax = test_level_izmax

! ------------------------------
! Call MESHGEN for the finer/coarser level mesh.

      CALL MESHGEN ( level, base_grid_level, finest_grid_level,
     1      ioffx, ioffy, ioffz, smallnum,

     1      base_grid_nx, base_grid_ny, base_grid_nz,
     1      base_grid_firsthalf_nx, base_grid_firsthalf_ny,
     1      base_grid_firsthalf_nz,

     1      base_grid_xmin, base_grid_xmax,
     1      base_grid_ymin, base_grid_ymax,
     1      base_grid_zmin, base_grid_zmax,
     1      base_grid_xratio, base_grid_yratio, base_grid_zratio,

     1      my_level_ixmin, my_level_ixmax,
     1      my_level_iymin, my_level_iymax,
     1      my_level_izmin, my_level_izmax,

     1      my_xmin, my_xmax, my_dx0, my_xratio,
     1      my_cut_ix0, my_cut_xb0,

     1      my_ymin, my_ymax, my_dy0, my_yratio,
     1      my_cut_iy0, my_cut_yb0,

     1      my_zmin, my_zmax, my_dz0, my_zratio,
     1      my_cut_iz0, my_cut_zb0)

! ------------------------------
! Now build the finer/coarser level mesh.

! ------------------------------
! Build mesh in the x-direction.

      IF ( my_cut_ix0 .LE. 0) THEN

! -----

        my_xb_save ( 0) = my_xmin

! -----

        DO ixx = 1, my_level_ixmax + ioffx

        my_dx_save ( ixx) = my_dx0 * my_xratio**ixx / my_xratio
        my_xb_save ( ixx) = my_xb_save ( ixx - 1)
     1                    + my_dx_save ( ixx)
        my_xc_save ( ixx) = 0.5 * ( my_xb_save ( ixx)
     1                            + my_xb_save ( ixx - 1) )

        END DO

! -----

        DO ixx = 0, my_level_ixmin - ioffx, -1

        my_dx_save ( ixx) = my_dx_save ( my_level_ixmin - ixx)
        my_xb_save ( ixx - 1) = my_xb_save ( ixx)
     1                        - my_dx_save ( ixx)
        my_xc_save ( ixx) = 0.5 * ( my_xb_save ( ixx)
     1                            + my_xb_save ( ixx - 1) )

        END DO

! -----

      ELSE

! -----

        my_xb_save ( my_cut_ix0 - 1) = my_cut_xb0

! -----

        DO ixx = my_cut_ix0, my_level_ixmax + ioffx

        my_dx_save ( ixx) = my_dx0
     1          * my_xratio**( ixx - my_cut_ix0 + 1) / my_xratio
        my_xb_save ( ixx) = my_xb_save ( ixx - 1)
     1                    + my_dx_save ( ixx)
        my_xc_save ( ixx) = 0.5 * ( my_xb_save ( ixx)
     1                            + my_xb_save ( ixx - 1) )

        END DO

! -----

        DO ixx = my_cut_ix0 - 1, my_level_ixmin - ioffx, -1

        my_dx_save ( ixx) = my_dx0
     1          * my_xratio**( my_cut_ix0 - ixx) / my_xratio
        my_xb_save ( ixx - 1) = my_xb_save ( ixx)
     1                        - my_dx_save ( ixx)
        my_xc_save ( ixx) = 0.5 * ( my_xb_save ( ixx)
     1                            + my_xb_save ( ixx - 1) )

        END DO

! -----

      END IF

! ------------------------------
! Build mesh in the y-direction.

      IF ( my_cut_iy0 .LE. 0) THEN

! -----

        my_yb_save ( 0) = my_ymin

! -----

        DO iyy = 1, my_level_iymax + ioffy

        my_dy_save ( iyy) = my_dy0 * my_yratio**iyy / my_yratio
        my_yb_save ( iyy) = my_yb_save ( iyy - 1)
     1                    + my_dy_save ( iyy)
        my_yc_save ( iyy) = 0.5 * ( my_yb_save ( iyy)
     1                            + my_yb_save ( iyy - 1) )

        END DO

! -----

        DO iyy = 0, my_level_iymin - ioffy, -1

        my_dy_save ( iyy) = my_dy_save ( my_level_iymin - iyy)
        my_yb_save ( iyy - 1) = my_yb_save ( iyy)
     1                        - my_dy_save ( iyy)
        my_yc_save ( iyy) = 0.5 * ( my_yb_save ( iyy)
     1                            + my_yb_save ( iyy - 1) )

        END DO

! -----

      ELSE

! -----

        my_yb_save ( my_cut_iy0 - 1) = my_cut_yb0

! -----

        DO iyy = my_cut_iy0, my_level_iymax + ioffy

        my_dy_save ( iyy) = my_dy0
     1          * my_yratio**( iyy - my_cut_iy0 + 1) / my_yratio
        my_yb_save ( iyy) = my_yb_save ( iyy - 1)
     1                    + my_dy_save ( iyy)
        my_yc_save ( iyy) = 0.5 * ( my_yb_save ( iyy)
     1                            + my_yb_save ( iyy - 1) )

        END DO

! -----

        DO iyy = my_cut_iy0 - 1, my_level_iymin - ioffy, -1

        my_dy_save ( iyy) = my_dy0
     1          * my_yratio**( my_cut_iy0 - iyy) / my_yratio
        my_yb_save ( iyy - 1) = my_yb_save ( iyy)
     1                        - my_dy_save ( iyy)
        my_yc_save ( iyy) = 0.5 * ( my_yb_save ( iyy)
     1                            + my_yb_save ( iyy - 1) )

        END DO

! -----

      END IF

! ------------------------------
! Build mesh in the z-direction.

      IF ( my_cut_iz0 .LE. 0) THEN

! -----

        my_zb_save ( 0) = my_zmin

! -----

        DO izz = 1, my_level_izmax + ioffz

        my_dz_save ( izz) = my_dz0 * my_zratio**izz / my_zratio
        my_zb_save ( izz) = my_zb_save ( izz - 1)
     1                    + my_dz_save ( izz)
        my_zc_save ( izz) = 0.5 * ( my_zb_save ( izz)
     1                            + my_zb_save ( izz - 1) )

        END DO

! -----

        DO izz = 0, my_level_izmin - ioffz, -1

        my_dz_save ( izz) = my_dz_save ( my_level_izmin - izz)
        my_zb_save ( izz - 1) = my_zb_save ( izz)
     1                        - my_dz_save ( izz)
        my_zc_save ( izz) = 0.5 * ( my_zb_save ( izz)
     1                            + my_zb_save ( izz - 1) )

        END DO

! -----

      ELSE

! -----

        my_zb_save ( my_cut_iz0 - 1) = my_cut_zb0

! -----

        DO izz = my_cut_iz0, my_level_izmax + ioffz

        my_dz_save ( izz) = my_dz0
     1          * my_zratio**( izz - my_cut_iz0 + 1) / my_zratio
        my_zb_save ( izz) = my_zb_save ( izz - 1)
     1                    + my_dz_save ( izz)
        my_zc_save ( izz) = 0.5 * ( my_zb_save ( izz)
     1                            + my_zb_save ( izz - 1) )

        END DO

! -----

        DO izz = my_cut_iz0 - 1, my_level_izmin - ioffz, -1

        my_dz_save ( izz) = my_dz0
     1          * my_zratio**( my_cut_iz0 - izz) / my_zratio
        my_zb_save ( izz - 1) = my_zb_save ( izz)
     1                        - my_dz_save ( izz)
        my_zc_save ( izz) = 0.5 * ( my_zb_save ( izz)
     1                            + my_zb_save ( izz - 1) )

        END DO

! -----

      END IF

! ------------------------------
! Write out the information about the finer/coarser level mesh in the
! x-direction.

      WRITE ( 6, *)" ---------------------------------"
      WRITE ( 6, *)" finer/coarser level mesh in the x-direction"
      WRITE ( 6, *)" ---------------------------------"
      ixx = my_level_ixmin - ioffx - 1
      WRITE ( 6, *)" ixx = ", ixx
      WRITE ( 6, *)"my_xb_save = ", my_xb_save ( ixx)
      WRITE ( 6, *)" "

      DO ixx = my_level_ixmin - ioffx, my_level_ixmax + ioffx

      WRITE ( 6, *)" ixx = ", ixx
      WRITE ( 6, *)"my_dx_save = ", my_dx_save ( ixx),
     1             " my_xc_save = ", my_xc_save ( ixx)
      WRITE ( 6, *)"my_xb_save = ", my_xb_save ( ixx)
      WRITE ( 6, *)" "

      END DO

! ------------------------------
! Write out the information about the finer/coarser level mesh in the
! y-direction.

      WRITE ( 6, *)" ---------------------------------"
      WRITE ( 6, *)" finer/coarser level mesh in the y-direction"
      WRITE ( 6, *)" ---------------------------------"
      iyy = my_level_iymin - ioffy - 1
      WRITE ( 6, *)" iyy = ", iyy
      WRITE ( 6, *)"my_yb_save = ", my_yb_save ( iyy)
      WRITE ( 6, *)" "

      DO iyy = my_level_iymin - ioffy, my_level_iymax + ioffy

      WRITE ( 6, *)" iyy = ", iyy
      WRITE ( 6, *)"my_dy_save = ", my_dy_save ( iyy),
     1             " my_yc_save = ", my_yc_save ( iyy)
      WRITE ( 6, *)"my_yb_save = ", my_yb_save ( iyy)
      WRITE ( 6, *)" "

      END DO

! ------------------------------
! Write out the information about the finer/coarser level mesh in the
! z-direction.

      WRITE ( 6, *)" ---------------------------------"
      WRITE ( 6, *)" finer/coarser level mesh in the z-direction"
      WRITE ( 6, *)" ---------------------------------"
      izz = my_level_izmin - ioffz - 1
      WRITE ( 6, *)" izz = ", izz
      WRITE ( 6, *)"my_zb_save = ", my_zb_save ( izz)
      WRITE ( 6, *)" "

      DO izz = my_level_izmin - ioffz, my_level_izmax + ioffz

      WRITE ( 6, *)" izz = ", izz
      WRITE ( 6, *)"my_dz_save = ", my_dz_save ( izz),
     1             " my_zc_save = ", my_zc_save ( izz)
      WRITE ( 6, *)"my_zb_save = ", my_zb_save ( izz)
      WRITE ( 6, *)" "

      END DO

! ----------------------------------------------------------------------

      END PROGRAM TEST_MESHGEN
#endif
! End of debugging test code for MESHGEN.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_PRIMITIVES_MIN_MAX

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save,

     1        scan_ghost_zones,
     1        rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1        vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1        bfymin, bfymax, bfzmin, bfzmax, div_min, div_max,
     1        speciemin, speciemax)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Obtains primitive variables from conserved variables. Also obtains their
! min and max values as well as the divergence.

! Useful in producing on-processor diagnostics for AMR simulations or
! smaller tasks like flagging for refinement.

! Also produces "divvelgr, msonicgr, flattengr" which can be useful during
! prolongation.

! "scan_ghost_zones = 1/0" depending on whether we do or don't want ghost
! zones scanned for minima and maxima.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe, scan_ghost_zones

      REAL rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax, div_min, div_max

      REAL speciemin ( 0: NFLUID), speciemax ( 0: NFLUID)

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy


      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr, workgr

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
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1

      REAL ob_8pi

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b, c,
     1      dx, dy, dz, ob_dx, ob_dy, ob_dz

! ----------------------------------------------------------------------

! Variables that go through the EOS subroutine boundary.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      rho, eng, prs, bprs, temperature, entropy,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)

      IF ( n_cc_components .NE. NFIELD + 1 + NFLUID) THEN
      WRITE ( 6, *)"Stopping in MAKE_PRIMITIVES_MIN_MAX :",
     1             " n_cc_components = ", n_cc_components
      STOP
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in MAKE_PRIMITIVES_MIN_MAX : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in MAKE_PRIMITIVES_MIN_MAX : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in MAKE_PRIMITIVES_MIN_MAX : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in MAKE_PRIMITIVES_MIN_MAX : NFIE","LD"
      WRITE ( 6, *) NFIELD, ".NE. 7"
      STOP
      END IF

      IF ( ( scan_ghost_zones .NE. 0) .AND.
     1     ( scan_ghost_zones .NE. 1) ) THEN
      WRITE ( 6, *) "Stopping in MAKE_PRIMITIVES_MIN_MAX :",
     1              " scan_ghost_zones .NE. 0 or 1 ", scan_ghost_zones
      STOP
      END IF

! ----------------------------------------------------------------------

! Make zone centered primitive variables. The 6, 7, 8 components of the
! conserved variables contain the zone-averaged magnetic fields.

! Zone centered magnetic fields are always evaluated with second order
! accuracy here.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, rho, eng, bprs, specie, prs,         &
!$OMP&   temperature, entropy, gamma_eos, gamma_soundspeed, mean_mol_wt)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ob_8pi, gasconst, gamma, smallnum, prsfloor, rhofloor, molewt,&
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, entropygr,                            &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr,                            &
!$OMP&    mean_mol_wt_gr, msonicgr)                                     &
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

          u ( ixx, iyy, izz, :, 5) = 0.0 ! Zero out ADER's time rate of update.

          u ( ixx, iyy, izz, 1, 1) = AMAX1 ( rhofloor,
     1                               u ( ixx, iyy, izz, 1, 1))

          rhogr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 1, 1)

          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)

          vlxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 2, 1) * a ( ixx)
          vlygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 3, 1) * a ( ixx)
          vlzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 4, 1) * a ( ixx)

          bfxgr ( ixx, iyy, izz) =0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                                  + bfxbdy ( ixxm1, iyy, izz, 1))
          bfygr ( ixx, iyy, izz) =0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                                  + bfybdy ( ixx, iyym1, izz, 1))
          bfzgr ( ixx, iyy, izz) =0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                                  + bfzbdy ( ixx, iyy, izzm1, 1))

          u ( ixx, iyy, izz, 6, 1) = bfxgr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 7, 1) = bfygr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 8, 1) = bfzgr ( ixx, iyy, izz)

          b ( ixx) = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          bprs ( ixx) = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          rho ( ixx) = rhogr ( ixx, iyy, izz)

          eng ( ixx) = ( u ( ixx, iyy, izz, 5, 1) - b ( ixx)
     1                                            - bprs ( ixx) )

#if ( NFLUID >= 1)
          specie ( 1:NFLUID, ixx) =
     1         u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) * a ( ixx)
#endif

          END DO ! End of DO ixx = loop.

! ----------

! Choose between a simple or more complicated EOS. Species fraction 
! can be strictly restricted to lie between 0 and 1 if it is copied back.


#if ( NFLUID == 0)

          DO ixx = 1 - ioffx, ix1 + ioffx

          gamma_eos ( ixx) = gamma
          gamma_soundspeed ( ixx) = gamma
          mean_mol_wt ( ixx) = molewt

          prs ( ixx) = AMAX1 ( eng ( ixx) * ( gamma - 1.0), prsfloor)
          temperature ( ixx) = molewt * prs ( ixx)
     1                       / ( gasconst * rho ( ixx))

          END DO ! End of DO ixx = loop.

#else

      CALL RHO_ENG_IN_PRS_T_GAMA_OUT
     1        ( 1, ix1, ioffx,

     1          gasconst, gamma, smallnum, prsfloor, rhofloor,
     1          cp_specie, molewt_specie, heat_of_formation,

     1          rho, eng, specie,

     1          prs, temperature, entropy,
     1          mean_mol_wt, gamma_eos, gamma_soundspeed)
          
#endif

! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)

          prsgr ( ixx, iyy, izz) = prs ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
! This is the only place in the code where the "entropygr" variable is
! rebuilt for the present timestep.
          entropygr ( ixx, iyy, izz, 1) = prs ( ixx)
     1                         / rho ( ixx)**( gamma_eos ( ixx) - 1.0)
          entropygr ( ixx, iyy, izz, 2: 5) = 0.0
#endif

          gamma_eos_gr ( ixx, iyy, izz) = gamma_eos ( ixx)
          gamma_soundspeed_gr ( ixx, iyy, izz) = gamma_soundspeed ( ixx)
          mean_mol_wt_gr ( ixx, iyy, izz) = mean_mol_wt ( ixx)

#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) =
     1    specie ( 1:NFLUID, ixx) * rhogr ( ixx, iyy, izz)
#endif

! Choose either the full magnetosonic speed or sound speed for the flattener.
          msonicgr ( ixx, iyy, izz) =
!     1    ( 2.0 * bprs ( ixx) + gamma_soundspeed ( ixx)
!     1                      * prs ( ixx) ) * a ( ixx)     ! full msonic speed
     1     gamma_soundspeed ( ixx) * prs ( ixx) * a ( ixx) ! sound speed
          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! Make the undivided divergence of the velocity, "divvelgr". Use it along
! with the minimal value of "msonicgr" in a halo of zones around the
! current one to make the flattener, "flattengr".

! This is a minimalist interpretation of the flattener without upwinding.
! It is adequate for prolongation problems where it might be beneficial
! to include some bit of the physics of flattening into the reconstruction
! that is done during the prolongation step.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, dx, dy, dz, ob_dx, ob_dy, ob_dz)     &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, flatten_coef, xb_save, dx_save, xc_save,               &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    vlxgr, vlygr, vlzgr, divvelgr, msonicgr, flattengr, workgr)   &
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

! ----------

          IF ( igeom .EQ. 1) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 2) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 3) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1              * dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

! ----------

          divvelgr ( ixx, iyy, izz) =
     1    ( vlxgr ( ixxp1, iyy, izz) - vlxgr ( ixxm1, iyy, izz) )
     1  * ob_dx ( ixx)
     1  + ( vlygr ( ixx, iyyp1, izz) - vlygr ( ixx, iyym1, izz) )
     1  * ob_dy ( ixx)
     1  + ( vlzgr ( ixx, iyy, izzp1) - vlzgr ( ixx, iyy, izzm1) )
     1  * ob_dz ( ixx)

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)
     1  * AMIN1 ( dx ( ixx), dy ( ixx), dz ( ixx) )

          workgr ( ixx, iyy, izz) = AMIN1 (

     1  msonicgr ( ixxm1, iyy, izz), msonicgr ( ixx, iyy, izz),
     1  msonicgr ( ixxp1, iyy, izz),
     1  msonicgr ( ixxm1, iyyp1, izz), msonicgr ( ixx, iyyp1, izz),
     1  msonicgr ( ixxp1, iyyp1, izz),
     1  msonicgr ( ixxm1, iyym1, izz), msonicgr ( ixx, iyym1, izz),
     1  msonicgr ( ixxp1, iyym1, izz),

     1  msonicgr ( ixxm1, iyy, izzp1), msonicgr ( ixx, iyy, izzp1),
     1  msonicgr ( ixxp1, iyy, izzp1),
     1  msonicgr ( ixxm1, iyyp1, izzp1), msonicgr ( ixx, iyyp1, izzp1),
     1  msonicgr ( ixxp1, iyyp1, izzp1),
     1  msonicgr ( ixxm1, iyym1, izzp1), msonicgr ( ixx, iyym1, izzp1),
     1  msonicgr ( ixxp1, iyym1, izzp1),

     1  msonicgr ( ixxm1, iyy, izzm1), msonicgr ( ixx, iyy, izzm1),
     1  msonicgr ( ixxp1, iyy, izzm1),
     1  msonicgr ( ixxm1, iyyp1, izzm1), msonicgr ( ixx, iyyp1, izzm1),
     1  msonicgr ( ixxp1, iyyp1, izzm1),
     1  msonicgr ( ixxm1, iyym1, izzm1), msonicgr ( ixx, iyym1, izzm1),
     1  msonicgr ( ixxp1, iyym1, izzm1)

     1                                    )

          a ( ixx) = ABS ( divvelgr ( ixx, iyy, izz))
     1      / ( flatten_coef * workgr ( ixx, iyy, izz))

          a ( ixx) = a ( ixx) - 1.0

          IF ( a ( ixx) .LE. 0.0) THEN
            flattengr ( ixx, iyy, izz) = 0.0
          ELSE
            flattengr ( ixx, iyy, izz) = AMIN1 ( 1.0, a ( ixx))
          END IF

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! Optionally evaluate and write out the scaled, undivided divergence here.
! In regions with very small magnetic field, it might not be wise to scale
! the undivided divergence by the magnitude of the magnetic field.
! In such situations, either don't scale the divergence at all or
! uncomment the klugey fix-up provided below.

! This is also an opportunity to evaluate and write out the min and
! max of the flow variables if that is desired.

      div_min = 1.0e30
      div_max = - 1.0e30

      rhomin = 1.0e30
      rhomax = - 1.0e30

      vlxmin = 1.0e30
      vlxmax = - 1.0e30
      vlymin = 1.0e30
      vlymax = - 1.0e30
      vlzmin = 1.0e30
      vlzmax = - 1.0e30

      prsmin = 1.0e30
      prsmax = - 1.0e30

      bfxmin = 1.0e30
      bfxmax = - 1.0e30
      bfymin = 1.0e30
      bfymax = - 1.0e30
      bfzmin = 1.0e30
      bfzmax = - 1.0e30

      speciemin ( :) = 1.0e30
      speciemax ( :) = - 1.0e30

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, dx, dy, dz,                  & 
!$OMP&    ob_dx, ob_dy, ob_dz, a, b, c)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    scan_ghost_zones,                                             &
!$OMP&    igeom, smallnum, xb_save, dx_save, xc_save,                   &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    & 
!$OMP&    bfxbdy, bfybdy, bfzbdy)                                       &
!$OMP&           REDUCTION ( MIN: div_min) REDUCTION ( MAX: div_max)    &
!$OMP&           REDUCTION ( MIN: rhomin) REDUCTION ( MAX: rhomax)      &
!$OMP&           REDUCTION ( MIN: vlxmin) REDUCTION ( MAX: vlxmax)      &
!$OMP&           REDUCTION ( MIN: vlymin) REDUCTION ( MAX: vlymax)      &
!$OMP&           REDUCTION ( MIN: vlzmin) REDUCTION ( MAX: vlzmax)      &
!$OMP&           REDUCTION ( MIN: prsmin) REDUCTION ( MAX: prsmax)      &
!$OMP&           REDUCTION ( MIN: bfxmin) REDUCTION ( MAX: bfxmax)      &
!$OMP&           REDUCTION ( MIN: bfymin) REDUCTION ( MAX: bfymax)      &
!$OMP&           REDUCTION ( MIN: bfzmin) REDUCTION ( MAX: bfzmax)      &
!$OMP&           REDUCTION ( MIN: speciemin) REDUCTION ( MAX: speciemax)&
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - scan_ghost_zones * ioffz,
     1                             iz1 + scan_ghost_zones * ioffz
        DO iyy = 1 - scan_ghost_zones * ioffy,
     1                             iy1 + scan_ghost_zones * ioffy
          DO ixx = 1 - scan_ghost_zones * ioffx,
     1                             ix1 + scan_ghost_zones * ioffx

! ----------

          IF ( igeom .EQ. 1) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 2) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

          IF ( igeom .EQ. 3) THEN

          dx ( ixx) = dx_save ( ixx)
          dy ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          dz ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1              * dz_save ( izz)

          ob_dx ( ixx) = 1.0 / dx ( ixx)
          ob_dy ( ixx) = 1.0 / dy ( ixx)
          ob_dz ( ixx) = 1.0 / dz ( ixx)

          END IF

! ----------

          a ( ixx) = ( bfxbdy ( ixx, iyy, izz, 1)
     1               - bfxbdy ( ixx - 1, iyy, izz, 1) ) * ob_dx ( ixx)
     1             + ( bfybdy ( ixx, iyy, izz, 1)
     1               - bfybdy ( ixx, iyy - 1, izz, 1) ) * ob_dy ( ixx)
     1             + ( bfzbdy ( ixx, iyy, izz, 1)
     1               - bfzbdy ( ixx, iyy, izz - 1, 1) ) * ob_dz ( ixx)

          a ( ixx) = a ( ixx) * AMIN1 ( dx ( ixx), dy ( ixx),
     1                                             dz ( ixx) )

          b ( ixx) = u ( ixx, iyy, izz, 6, 1)**2
     1             + u ( ixx, iyy, izz, 7, 1)**2
     1             + u ( ixx, iyy, izz, 8, 1)**2
          b ( ixx) = SQRT ( b ( ixx))
          b ( ixx) = AMAX1 ( b ( ixx), smallnum)

          a ( ixx) = a ( ixx) / b ( ixx)
!          a ( ixx) = a ( ixx) / AMAX1 ( b ( ixx),
!     1                     msonicgr ( ixx, iyy, izz)
!     1                     * SQRT ( rhogr ( ixx, iyy, izz) ) )

          div_min = AMIN1 ( div_min, a ( ixx))
          div_max = AMAX1 ( div_max, a ( ixx))

! ----------

          rhomin = AMIN1 ( rhomin, rhogr ( ixx, iyy, izz))
          rhomax = AMAX1 ( rhomax, rhogr ( ixx, iyy, izz))

          vlxmin = AMIN1 ( vlxmin, vlxgr ( ixx, iyy, izz))
          vlxmax = AMAX1 ( vlxmax, vlxgr ( ixx, iyy, izz))
          vlymin = AMIN1 ( vlymin, vlygr ( ixx, iyy, izz))
          vlymax = AMAX1 ( vlymax, vlygr ( ixx, iyy, izz))
          vlzmin = AMIN1 ( vlzmin, vlzgr ( ixx, iyy, izz))
          vlzmax = AMAX1 ( vlzmax, vlzgr ( ixx, iyy, izz))

          prsmin = AMIN1 ( prsmin, prsgr ( ixx, iyy, izz))
          prsmax = AMAX1 ( prsmax, prsgr ( ixx, iyy, izz))

          bfxmin = AMIN1 ( bfxmin, bfxgr ( ixx, iyy, izz))
          bfxmax = AMAX1 ( bfxmax, bfxgr ( ixx, iyy, izz))
          bfymin = AMIN1 ( bfymin, bfygr ( ixx, iyy, izz))
          bfymax = AMAX1 ( bfymax, bfygr ( ixx, iyy, izz))
          bfzmin = AMIN1 ( bfzmin, bfzgr ( ixx, iyy, izz))
          bfzmax = AMAX1 ( bfzmax, bfzgr ( ixx, iyy, izz))

#if ( NFLUID >= 1)
          c ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)

          speciemin ( 1: NFLUID) = AMIN1 ( speciemin ( 1: NFLUID),
     1    u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) * c ( ixx) )

          speciemax ( 1: NFLUID) = AMAX1 ( speciemax ( 1: NFLUID),
     1    u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) * c ( ixx) )
#endif

! ----------

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_PRIMITIVES_MIN_MAX

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------













