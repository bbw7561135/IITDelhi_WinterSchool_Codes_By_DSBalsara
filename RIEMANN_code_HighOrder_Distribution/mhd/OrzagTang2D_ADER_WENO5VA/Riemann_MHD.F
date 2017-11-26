! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! This flag is set to "1" only for testing the 1D RS capabilities of the
! "MAKE_FLUX_PTWISE_?" subroutines. When it is set to "1", we use pure
! 1D Riemann solver in its pointwise form to build facial fluxes.
! When this flag is set to "0", the fluxes also take on contributions
! from the 2D Riemann solver.
! Electric fields are always built using the 2D Riemann solver.
#define USE_ONLY_1D_RS 0

! The "TRUE_2D" is only set to 1 if the problem is truly 2d with no
! variation in the third direction. It enables the 2d RS to operate more
! efficiently in the dimensionality where it comes into its own.
#define TRUE_2D 1

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM RIEMANN_MHD

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This is a very fast, three dimensional, fully
! functional MHD code that does divergence-free MHD using a
! second, third and fourth order accurate in time ADER
! timestep. The spatial interpolation is either TVD or WENO.
! This code was put together to go along with Balsara's text book.
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
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,

     1        ixx, iyy, izz, i_comp, i_rad, ifluid, istep, ntstep,
     1        nimgstep, ndumpstep

! Set size of problem domain.

      PARAMETER ( ix1 = 400, iy1 = 400, iz1 = 1,
     1            ioffx = 6, ioffy = 6, ioffz = 0,
     1            n_cc_components = NFIELD + 1 + NFLUID,
     1            n_cc_modes = NUM_SPACETIME_BASES, 
     1            n_fc_modes = 1)

      REAL time, dtcur, dtcurp1, timestop, timeimage, dtimage,
     1     dt_shift,

     1     xmin, xmax, ymin, ymax, zmin, zmax,

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

! Save rates of change of the source terms, though not for all the modes.
! This may become a problem for PNPM schemes with N > 0.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: dudt

! "bfxbdy, bfybdy, bfzbdy" hold facial magnetic fields and their
! moments in the transverse directions. "bfxbdy" is collocated at the
! center of the x-face and so on...

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

! "entropygr" holds the entropy density and is useful when pressure
! positivity switch is turned on.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

! "flux_x, flux_y, flux_z" hold the fluxes in each of the three directions.

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_y, flux_z,
     1                                      flux_temp

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
     1      entropy_flux_x, entropy_flux_y, entropy_flux_z

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

      gamma = 5.0/3.0     ! Polytropic index for the gas.
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

! Courant number: "cfl_coef" is based on the entire signal speed.
! When multi-d Riemann solvers are used, they can be set in the ranges
! [0.3,0.45] in 3d; [0.5,0.65] in 2d and [0.6,0.8] in 1d.
! At fourth order this number needs to be reduced because we are not using
! fourth order Gaussian integration in time.
      cfl_coef = 0.4

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
      prs_flatten_coef = 1.3

      fracdiff = 0.1  !Set fractional change needed for entropy fix.(deprecated)

! If the divergence of the velocity at a zone boundary exceeds this coef times
! the sound speed then we have to blend in some fraction of HLL RS.
      divratio_einfeldt_rs = 0.4

! Set up variables for ADER timestep iteration.
      ader_tolerance = 1.0e-16 ! Convergence tolerance for ADER, set it v.low. 
      max_iter_ader_cg = NUMORDER-1 ! Number of iterations; for CG, set it low.
      max_iter_ader_dg = NUMORDER  ! Number of iterations for ADER-DG.

! Set a threshold value that triggers the Lapidus viscosity as well as
! its value.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05

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
      dtimage = 1.0 / 20.0

! Time at which you want the simulation to stop.
      timestop = 1.0

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
      WRITE ( 6, *)"cfl_coef, mc_coef = ", cfl_coef, mc_coef
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
      WRITE ( 6, *)"ader_tolerance = ", ader_tolerance
      WRITE ( 6, *)"max_iter_ader_cg, max_iter_ader_dg = ",
     1              max_iter_ader_cg, max_iter_ader_dg
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Algorithms chosen for this run : "
      WRITE ( 6, *)" "

      WRITE ( 6, *)"Order of accuracy = ", NUMORDER
      WRITE ( 6, *)" "

#if ( SLOPE_LIMITER == 1)
      WRITE ( 6, *)"TVD-style limiter of appropriate order being used."
#else
      WRITE ( 6, *)"WENO limiter being used."
#endif
      WRITE ( 6, *)" "

#if ( LIMIT_TWICE == 1)
      WRITE ( 6, *)"Limiting done on conserved & characteristic vars."
#else
      WRITE ( 6, *)"Limiting done only on conserved vars."
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

#if ( BUILD_SOURCE_TERMS == 1)
      WRITE ( 6, *)"Source terms will be built."
#else
      WRITE ( 6, *)"Source terms will not be built."
#endif
      WRITE ( 6, *)" "

#if ( PRESSURE_POSITIVITY == 1)
      WRITE ( 6, *)"Positivity of density and pressure will be",
     1             " enforced as much as possible."
#else
      WRITE ( 6, *)"Positivity of density and pressure will",
     1             " not be enforced"
#endif
      WRITE ( 6, *)" "

      WRITE ( 6, *)"--------------------------------------------------"

! ----------------------------------------------------------------------

! If dimensioning is inconsistent, then write out some info. and stop.

      IF ( ( ix1 .GT. 1) .AND.
     1     ( ioffx .LT. 4 + 2 * ( NUMORDER / 4) ) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : increase ioffx = ", ioffx
      STOP
      END IF

      IF ( ( ix1 .EQ. 1) .AND. ( ioffx .NE. 0) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : ioffx must be 0 ", ioffx
      STOP
      END IF

      IF ( ( iy1 .GT. 1) .AND.
     1     ( ioffy .LT. 4 + 2 * ( NUMORDER / 4) ) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : increase ioffy = ", ioffy
      STOP
      END IF

      IF ( ( iy1 .EQ. 1) .AND. ( ioffy .NE. 0) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : ioffy must be 0 ", ioffy
      STOP
      END IF

      IF ( ( iz1 .GT. 1) .AND.
     1     ( ioffz .LT. 4 + 2 * ( NUMORDER / 4) ) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : increase ioffz = ", ioffz
      STOP
      END IF

      IF ( ( iz1 .EQ. 1) .AND. ( ioffz .NE. 0) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : ioffz must be 0 ", ioffz
      STOP
      END IF

      IF ( NTRACERS > NFLUID) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_MHD : NTRACER","S > NFLUI","D",
     1              NTRACERS, NFLUID
      STOP
      END IF

      IF ( ( TRUE_2D == 1) .AND. ( ioffx * ioffy * ioffz > 0) ) THEN
      WRITE ( 6, *)
     1  "Stopping in RIEMANN_MHD : TRUE _2D == 1 for 3d problem"
      STOP
      END IF

! ----------------------------------------------------------------------

! Set the size of the patch being simulated.

      xmin = 0.0
      xmax = 2.0
      ymin = 0.0
      ymax = 2.0
      zmin = -0.5
      zmax = 0.5

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

! Popular ones: 1 for inflow; 2 for outflow/continuitive; 3 for reflective;
! 4 for periodic in serial setting.

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

      a = cfl_coef         ! Save CFL number.
      cfl_coef = 0.05      ! Start problem with reduced CFL number.

      CALL INIT_ORZAG_TANG_XY

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u, bfxbdy, bfybdy, bfzbdy,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

      cfl_coef = a           ! Copy back the original CFL number.

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

      READ ( 5, *) time, dtcur, i_comp

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

      WRITE ( 5, *) time, dtcur, n_cc_components

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

      WRITE ( 6, *) time, dtcur, n_cc_components

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

! Limit the slopes of all the flow variables. The time rate of update,
! i.e. the modes "u ( ixx, iyy, izz, :, NUM_ZONECEN_BASES + 1:)", 
! are also zeroed out here.
! The flattener is also built. The primitive variables are built.

      CALL LIMIT_ALL_MOMENTS ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Enforce positivity of density and pressure, if needed, for a v. stringent
! problem.

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
! -----

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
      WRITE ( 6, *)"Stopping in RIEMANN_MHD: time .GE. timestop"
      WRITE ( 6, *)"time, timestop = ", time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      STOP
      END IF

! ----------------------------------------------------------------------

! This is a predictor step and can be done either by using full-up
! Riemann solvers at zone boundaries or ( more inexpensively) by using
! the conserved variables and their modal variation within the zone.

! For stringent problems, the former is more stabilizing than the latter.
! Consequently, the code defaults to a situation where we use
! full-up Riemann solvers at zone boundaries (not implemented for ADER).

      CALL MAKE_ADER_UPDT_PTWISE ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, dudt, bfxbdy, bfybdy, bfzbdy,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Now that the time rate of update, i.e. the modes
! "u ( ixx, iyy, :, izz, NUM_ZONECEN_BASES + 1:)", are built
! within a zone, we are in a position to build space-time-averaged fluxes.
! These are built at face centers and electric fields at edge-centers.

! The fluxes in this corrector stage do need to be zeroed out before
! subroutine entry.

! ------------------------------

      dt_shift = 1.0

      flux_x ( :, :, :, :) = 0.0
      prs_flux_x ( :, :, :) = 0.0
      entropy_flux_x ( :, :, :) = 0.0

      flux_y ( :, :, :, :) = 0.0
      prs_flux_y ( :, :, :) = 0.0
      entropy_flux_y ( :, :, :) = 0.0

      flux_z ( :, :, :, :) = 0.0
      prs_flux_z ( :, :, :) = 0.0
      entropy_flux_z ( :, :, :) = 0.0

      efxgr ( :, :, :) = 0.0
      efygr ( :, :, :) = 0.0
      efzgr ( :, :, :) = 0.0

! ------------------------------

! Build "efxgr" using 2D RS and add contributions from the 2D RS
! to "flux_y, flux_z".

      CALL APPLY_MULTID_RS_X_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        flux_y, flux_z, flux_temp, efxgr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Build "efygr" using 2D RS and add contributions from the 2D RS 
! to "flux_x, flux_z".

      CALL APPLY_MULTID_RS_Y_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        flux_x, flux_z, flux_temp, efygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Build "efzgr" using 2D RS and add contributions from the 2D RS
! to "flux_x, flux_y".

      CALL APPLY_MULTID_RS_Z_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        flux_x, flux_y, flux_temp, efzgr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        divvelgr, msonicgr, flattengr,

     1        x_indx_limits, y_indx_limits, z_indx_limits,
     1        bcarr_x, bcarr_y, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Zero out "flux_x" and collect contributions to it from 1d RS.

      CALL MAKE_FLUX_PTWISE_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        expand_by_one_zone, dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_x, prs_flux_x,
     1        entropygr, entropy_flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Zero out "flux_y" and collect contributions to it from 1d RS.

      CALL MAKE_FLUX_PTWISE_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        expand_by_one_zone, dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_y, prs_flux_y,
     1        entropygr, entropy_flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Zero out "flux_z" and collect contributions to it from 1d RS.

      CALL MAKE_FLUX_PTWISE_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        expand_by_one_zone, dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_z, prs_flux_z,
     1        entropygr, entropy_flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Update the flow variables in this corrector step and get new "dtcurp1".

      CALL UPDATE_VARS_CORRECTOR

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, dudt, bfxbdy, bfybdy, bfzbdy, flux_x, flux_y, flux_z,
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

      END PROGRAM RIEMANN_MHD

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

      SUBROUTINE LIMIT_ALL_MOMENTS ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        u, bfxbdy, bfybdy, bfzbdy,
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
! stability. Set LIMIT_TWICE to 1 to invoke characteristic limiting.

! SLOPE_LIMITER == 2 yields WENO limiting at 2nd, 3rd and 4th order.
! SLOPE_LIMITER == 1 yields MC limiting at 2nd order and PPM at 3rd order. At
!                    4th order it yields a hybrid r=5 WENO + Van Albada limiter.

! "flattengr" applies a flattener in regions of strong shock.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

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

! Variables that are internal to this subroutine.

      INTEGER ixx, iyy, izz, i_comp, ixx0, iyy0, izz0,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,
     1        ixxp3, iyyp3, izzp3, ixxm3, iyym3, izzm3,
     1        limit_characteristics

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: workgr

      REAL ob_8pi, ob_6, ob_12, ob_24, ob_60, r2_3, r13_12, r13_3,
     1     r781_20,

     1     a, b, c,d, e, f, g, h,
     1     psi_limiter, var_zone_min, var_zone_max,
     1     var_neighbor_min,var_neighbor_max,
     1     MinMod, MC_beta,

     1     u_ip2, u_ip1, u_i, u_im1, u_im2, u_im3, central_zone_wt,
     1     R1_0, R2_0, R3_0, R1_x, R2_x, R3_x, R1_x2, R2_x2, R3_x2,
     1     R1_x3, R2_x3, R3_x3,
     1     R1_xy, R2_xy, R3_xy, R4_xy, R5_xy,
     1     R1_yz, R2_yz, R3_yz, R4_yz, R5_yz,
     1     R1_xz, R2_xz, R3_xz, R4_xz, R5_xz,
     1     R1_x2y, R2_x2y, R3_x2y, R4_x2y, R5_x2y,
     1     R1_xy2, R2_xy2, R3_xy2, R4_xy2, R5_xy2,
     1     R1_y2z, R2_y2z, R3_y2z, R4_y2z, R5_y2z,
     1     R1_yz2, R2_yz2, R3_yz2, R4_yz2, R5_yz2,
     1     R1_x2z, R2_x2z, R3_x2z, R4_x2z, R5_x2z,
     1     R1_xz2, R2_xz2, R3_xz2, R4_xz2, R5_xz2,
     1     R1_xyz, R2_xyz, R3_xyz, R4_xyz, R5_xyz, R6_xyz,
     1     R7_xyz, R8_xyz,

     1     Is1, Is2, Is3, Is4, Is5, Is6, Is7, Is8,
     1     Is1_4, Is2_4, Is3_4, Is4_4, Is5_4,
     1     Wt1, Wt2, Wt3, Wt4, Wt5, Wt6, Wt7, Wt8,

     1     div_min, div_max, dx, dy, dz, ob_dx, ob_dy, ob_dz,
     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, flatten_1d

      REAL, DIMENSION ( 0: NFLUID, ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID, NUMORDER) :: u_1d

! MinMod limiter
      MinMod ( a, b) = AMIN1 ( ABS ( a), ABS ( b))
     1                    * ( SIGN ( 0.5, a) + SIGN ( 0.5, b))

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( a, b, mc_coef) = AMIN1 ( 0.5 * ( ABS ( a + b)),
     1                         mc_coef * ABS ( a), mc_coef * ABS ( b) )
     1                    * ( SIGN ( 0.5, a) + SIGN ( 0.5, b))

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

      ob_12 = 1.0 / 12.0

      ob_24 = 1.0 / 24.0

      ob_60 = 1.0 / 60.0

      r2_3 = 2.0 / 3.0

      r13_3 = 13.0 / 3.0

      r13_12 = 13.0 / 12.0

      r781_20 = 781.0 / 20.0

      central_zone_wt = 100.0

      limit_characteristics = LIMIT_TWICE

! ----------

      IF ( n_cc_components .NE. NFIELD + 1 + NFLUID) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS : n_cc_components = ",
     1              n_cc_components
      STOP
      END IF

      IF ( n_cc_modes .NE. NUM_SPACETIME_BASES) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS : n_cc_modes = ",
     1              n_cc_modes, NUM_SPACETIME_BASES
      STOP
      END IF

      IF ( n_fc_modes .NE. 1) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS : n_fc_modes = ",
     1              n_fc_modes
      STOP
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_MOMENTS : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_MOMENTS : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_MOMENTS : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

! ----------

! Make sure the number of offset zones is adequate at 2nd and 3rd order.

      IF ( ( NUMORDER .LE. 3) .AND. ( ix1 .GT. 1) .AND.
     1                              ( ioffx .LT. 4) ) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS, ioffx .LT. 4",
     1              NUMORDER, ioffx
      STOP
      END IF

      IF ( ( NUMORDER .LE. 3) .AND. ( iy1 .GT. 1) .AND.
     1                              ( ioffy .LT. 4) ) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS, ioffy .LT. 4",
     1              NUMORDER, ioffy
      STOP
      END IF

      IF ( ( NUMORDER .LE. 3) .AND. ( iz1 .GT. 1) .AND.
     1                              ( ioffz .LT. 4) ) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS, ioffz .LT. 4",
     1              NUMORDER, ioffz
      STOP
      END IF

! ----------

! Make sure the number of offset zones is adequate at 4th order.

      IF ( ( NUMORDER .EQ. 4) .AND. ( ix1 .GT. 1) .AND.
     1                              ( ioffx .LT. 6) ) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS, ioffx .LT. 6",
     1              NUMORDER, ioffx
      STOP
      END IF

      IF ( ( NUMORDER .EQ. 4) .AND. ( iy1 .GT. 1) .AND.
     1                              ( ioffy .LT. 6) ) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS, ioffy .LT. 6",
     1              NUMORDER, ioffy
      STOP
      END IF

      IF ( ( NUMORDER .EQ. 4) .AND. ( iz1 .GT. 1) .AND.
     1                              ( ioffz .LT. 6) ) THEN
      WRITE ( 6, *)"Stopping in LIMIT_ALL_MOMENTS, ioffz .LT. 6",
     1              NUMORDER, ioffz
      STOP
      END IF

! ----------

! MHD with Eos has 7 characteristic fields.

      IF ( NFIELD .NE. 7) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_MOMENTS : NFIE","LD"
      WRITE ( 6, *) NFIELD, ".NE. 7"
      STOP
      END IF

! ----------------------------------------------------------------------

! Make zone centered primitive variables.The 6, 7, 8 components of the
! conserved variables contain the zone-averaged magnetic fields.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP&    izzp3, izzm3, iyy, iyyp1, iyym1, iyyp2, iyym2, iyyp3, iyym3,  &
!$OMP&    ixx, ixxp1, ixxm1, ixxp2, ixxm2, ixxp3, ixxm3,                &
!$OMP&    u_ip2, u_ip1, u_i, u_im1, u_im2, u_im3,                       &
!$OMP&    R1_0, R2_0, R3_0, R1_x, R2_x, R3_x, R1_x2, R2_x2, R3_x2,      &
!$OMP&    R1_x3, R2_x3, R3_x3, Is1, Is2, Is3, Wt1, Wt2, Wt3,            &
!$OMP&    a, b, rho, eng, bprs, specie, prs,                            &
!$OMP&   temperature, entropy, gamma_eos, gamma_soundspeed, mean_mol_wt)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    ob_8pi, gasconst, gamma, smallnum, prsfloor, rhofloor, molewt,&
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    ob_6, ob_12, ob_24, ob_60, r13_3, r13_12, r781_20,            &
!$OMP&    central_zone_wt,                                              &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, entropygr,                            &
!$OMP&    temperaturegr, gamma_eos_gr, gamma_soundspeed_gr,             &
!$OMP&    mean_mol_wt_gr, msonicgr)                                     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz
      izzp3 = MIN0 ( izz + 3, iz1 + ioffz)
      izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
      izzm1 = MAX0 ( izz - 1, 1 - ioffz)
      izzm2 = MAX0 ( izz - 2, 1 - ioffz)
      izzm3 = MAX0 ( izz - 3, 1 - ioffz)

        DO iyy = 1 - ioffy, iy1 + ioffy
        iyyp3 = MIN0 ( iyy + 3, iy1 + ioffy)
        iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
        iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
        iyym2 = MAX0 ( iyy - 2, 1 - ioffy)
        iyym3 = MAX0 ( iyy - 3, 1 - ioffy)

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp3 = MIN0 ( ixx + 3, ix1 + ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)
          ixxm3 = MAX0 ( ixx - 3, 1 - ioffx)

! Zero out ADER's time rate of update at all orders.
          u ( ixx, iyy, izz, :,
     1        NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES) = 0.0

          u ( ixx, iyy, izz, 1, 1) = AMAX1 ( rhofloor,
     1                               u ( ixx, iyy, izz, 1, 1))

          rhogr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 1, 1)

          a = 1.0 / rhogr ( ixx, iyy, izz)

          vlxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 2, 1) * a
          vlygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 3, 1) * a
          vlzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 4, 1) * a

! ----------
#if ( NUMORDER == 2)

          u ( ixx, iyy, izz, 6, 1) = 0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                                  + bfxbdy ( ixxm1, iyy, izz, 1))
          u ( ixx, iyy, izz, 6, 2) = bfxbdy ( ixx, iyy, izz, 1)
     1                             - bfxbdy ( ixxm1, iyy, izz, 1)

          u ( ixx, iyy, izz, 7, 1) = 0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                                  + bfybdy ( ixx, iyym1, izz, 1))
          u ( ixx, iyy, izz, 7, 3) = bfybdy ( ixx, iyy, izz, 1)
     1                             - bfybdy ( ixx, iyym1, izz, 1)

          u ( ixx, iyy, izz, 8, 1) = 0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                                  + bfzbdy ( ixx, iyy, izzm1, 1))
          u ( ixx, iyy, izz, 8, 4) = bfzbdy ( ixx, iyy, izz, 1)
     1                             - bfzbdy ( ixx, iyy, izzm1, 1)

          bfxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 6, 1)
          bfygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 7, 1)
          bfzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 8, 1)

#endif
! End of "NUMORDER == 2".
! ----------
#if ( NUMORDER == 3)

! -----

          u_ip1 = bfxbdy ( ixxp1, iyy, izz, 1)
          u_i   = bfxbdy ( ixx, iyy, izz, 1)
          u_im1 = bfxbdy ( ixxm1, iyy, izz, 1)
          u_im2 = bfxbdy ( ixxm2, iyy, izz, 1)

          R1_0 = ( 8.0 * u_im1 - u_im2 + 5.0 * u_i) * ob_12
          R1_x = u_i - u_im1
          R1_x2 = 0.5 * ( u_i - 2.0 * u_im1 + u_im2)
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_0 = ( 8.0 * u_i - u_ip1 + 5 * u_im1) * ob_12
          R2_x = u_i - u_im1
          R2_x2 = 0.5 * ( u_ip1 - 2.0 * u_i + u_im1)
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = 1.0 / ( Is2 + smallnum)**4
          a = 1.0 / ( Wt1 + Wt2)
          Wt1 = a * Wt1
          Wt2 = a * Wt2

          u ( ixx, iyy, izz, 6, 1) = Wt1 * R1_0 + Wt2 * R2_0
          u ( ixx, iyy, izz, 6, 2) = Wt1 * R1_x + Wt2 * R2_x
          u ( ixx, iyy, izz, 6, 5) = Wt1 * R1_x2 + Wt2 * R2_x2
          bfxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 6, 1)

! -----

          u_ip1 = bfybdy ( ixx, iyyp1, izz, 1)
          u_i   = bfybdy ( ixx, iyy, izz, 1)
          u_im1 = bfybdy ( ixx, iyym1, izz, 1)
          u_im2 = bfybdy ( ixx, iyym2, izz, 1)

          R1_0 = ( 8.0 * u_im1 - u_im2 + 5.0 * u_i) * ob_12
          R1_x = u_i - u_im1
          R1_x2 = 0.5 * ( u_i - 2.0 * u_im1 + u_im2)
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_0 = ( 8.0 * u_i - u_ip1 + 5 * u_im1) * ob_12
          R2_x = u_i - u_im1
          R2_x2 = 0.5 * ( u_ip1 - 2.0 * u_i + u_im1)
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = 1.0 / ( Is2 + smallnum)**4
          a = 1.0 / ( Wt1 + Wt2)
          Wt1 = a * Wt1
          Wt2 = a * Wt2

          u ( ixx, iyy, izz, 7, 1) = Wt1 * R1_0 + Wt2 * R2_0
          u ( ixx, iyy, izz, 7, 3) = Wt1 * R1_x + Wt2 * R2_x
          u ( ixx, iyy, izz, 7, 6) = Wt1 * R1_x2 + Wt2 * R2_x2
          bfygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 7, 1)

! -----

          u_ip1 = bfzbdy ( ixx, iyy, izzp1, 1)
          u_i   = bfzbdy ( ixx, iyy, izz, 1)
          u_im1 = bfzbdy ( ixx, iyy, izzm1, 1)
          u_im2 = bfzbdy ( ixx, iyy, izzm2, 1)

          R1_0 = ( 8.0 * u_im1 - u_im2 + 5.0 * u_i) * ob_12
          R1_x = u_i - u_im1
          R1_x2 = 0.5 * ( u_i - 2.0 * u_im1 + u_im2)
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_0 = ( 8.0 * u_i - u_ip1 + 5 * u_im1) * ob_12
          R2_x = u_i - u_im1
          R2_x2 = 0.5 * ( u_ip1 - 2.0 * u_i + u_im1)
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = 1.0 / ( Is2 + smallnum)**4
          a = 1.0 / ( Wt1 + Wt2)
          Wt1 = a * Wt1
          Wt2 = a * Wt2

          u ( ixx, iyy, izz, 8, 1) = Wt1 * R1_0 + Wt2 * R2_0
          u ( ixx, iyy, izz, 8, 4) = Wt1 * R1_x + Wt2 * R2_x
          u ( ixx, iyy, izz, 8, 7) = Wt1 * R1_x2 + Wt2 * R2_x2
          bfzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 8, 1)

! -----
#endif
! End of "NUMORDER == 3".
! ----------
#if ( NUMORDER == 4)

! -----

          u_ip2 = bfxbdy ( ixxp2, iyy, izz, 1)
          u_ip1 = bfxbdy ( ixxp1, iyy, izz, 1)
          u_i   = bfxbdy ( ixx, iyy, izz, 1)
          u_im1 = bfxbdy ( ixxm1, iyy, izz, 1)
          u_im2 = bfxbdy ( ixxm2, iyy, izz, 1)
          u_im3 = bfxbdy ( ixxm3, iyy, izz, 1)

          R1_0 = ( 19.0 * u_im1 - 5.0 * u_im2 + u_im3
     1           + 9.0 * u_i) * ob_24
          R1_x = ( -57.0 * u_im1 - 3.0 * u_im2 + u_im3
     1           + 59.0 * u_i) * ob_60
          R1_x2 = 0.25 * ( -7.0 * u_im1 + 5.0 * u_im2 - u_im3
     1           + 3.0 * u_i)
          R1_x3 = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i) * ob_6
          Is1 = ( R1_x + 0.1 * R1_x3)**2 + r13_3 * R1_x2**2
     1           + r781_20 * R1_x3**2

          R2_0 = ( 13.0 * u_im1 - u_im2 + 13.0 * u_i - u_ip1) * ob_24
          R2_x = ( -63.0 * u_im1 + u_im2 + 63.0 * u_i - u_ip1) * ob_60
          R2_x2 = 0.25 * ( - u_im1 + u_im2 - u_i + u_ip1)
          R2_x3 = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1) * ob_6
          Is2 = ( R2_x + 0.1 * R2_x3)**2 + r13_3 * R2_x2**2
     1           + r781_20 * R2_x3**2

          R3_0 = ( 9.0 * u_im1 + 19.0 * u_i - 5.0 * u_ip1
     1           + u_ip2) * ob_24
          R3_x = ( -59.0 * u_im1 + 57.0 * u_i + 3.0 * u_ip1
     1           - u_ip2) * ob_60
          R3_x2 = 0.25 * ( 3.0 * u_im1 - 7.0 * u_i + 5.0 * u_ip1
     1           - u_ip2)
          R3_x3 = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2) * ob_6
          Is3 = ( R3_x + 0.1 * R3_x3)**2 + r13_3 * R3_x2**2
     1           + r781_20 * R3_x3**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4
          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          u ( ixx, iyy, izz, 6, 1) = Wt1 * R1_0 + Wt2 * R2_0
     1                             + Wt3 * R3_0
          u ( ixx, iyy, izz, 6, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                             + Wt3 * R3_x
          u ( ixx, iyy, izz, 6, 5) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                             + Wt3 * R3_x2
          u ( ixx, iyy, izz, 6, 11) = Wt1 * R1_x3 + Wt2 * R2_x3
     1                             + Wt3 * R3_x3
          bfxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 6, 1)

! -----

          u_ip2 = bfybdy ( ixx, iyyp2, izz, 1)
          u_ip1 = bfybdy ( ixx, iyyp1, izz, 1)
          u_i   = bfybdy ( ixx, iyy, izz, 1)
          u_im1 = bfybdy ( ixx, iyym1, izz, 1)
          u_im2 = bfybdy ( ixx, iyym2, izz, 1)
          u_im3 = bfybdy ( ixx, iyym3, izz, 1)

          R1_0 = ( 19.0 * u_im1 - 5.0 * u_im2 + u_im3
     1           + 9.0 * u_i) * ob_24
          R1_x = ( -57.0 * u_im1 - 3.0 * u_im2 + u_im3
     1           + 59.0 * u_i) * ob_60
          R1_x2 = 0.25 * ( -7.0 * u_im1 + 5.0 * u_im2 - u_im3
     1           + 3.0 * u_i)
          R1_x3 = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i) * ob_6
          Is1 = ( R1_x + 0.1 * R1_x3)**2 + r13_3 * R1_x2**2
     1           + r781_20 * R1_x3**2

          R2_0 = ( 13.0 * u_im1 - u_im2 + 13.0 * u_i - u_ip1) * ob_24
          R2_x = ( -63.0 * u_im1 + u_im2 + 63.0 * u_i - u_ip1) * ob_60
          R2_x2 = 0.25 * ( - u_im1 + u_im2 - u_i + u_ip1)
          R2_x3 = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1) * ob_6
          Is2 = ( R2_x + 0.1 * R2_x3)**2 + r13_3 * R2_x2**2
     1           + r781_20 * R2_x3**2

          R3_0 = ( 9.0 * u_im1 + 19.0 * u_i - 5.0 * u_ip1
     1           + u_ip2) * ob_24
          R3_x = ( -59.0 * u_im1 + 57.0 * u_i + 3.0 * u_ip1
     1           - u_ip2) * ob_60
          R3_x2 = 0.25 * ( 3.0 * u_im1 - 7.0 * u_i + 5.0 * u_ip1
     1           - u_ip2)
          R3_x3 = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2) * ob_6
          Is3 = ( R3_x + 0.1 * R3_x3)**2 + r13_3 * R3_x2**2
     1           + r781_20 * R3_x3**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4
          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          u ( ixx, iyy, izz, 7, 1) = Wt1 * R1_0 + Wt2 * R2_0
     1                             + Wt3 * R3_0
          u ( ixx, iyy, izz, 7, 3) = Wt1 * R1_x + Wt2 * R2_x
     1                             + Wt3 * R3_x
          u ( ixx, iyy, izz, 7, 6) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                             + Wt3 * R3_x2
          u ( ixx, iyy, izz, 7, 12) = Wt1 * R1_x3 + Wt2 * R2_x3
     1                             + Wt3 * R3_x3
          bfygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 7, 1)

! -----

          u_ip2 = bfzbdy ( ixx, iyy, izzp2, 1)
          u_ip1 = bfzbdy ( ixx, iyy, izzp1, 1)
          u_i   = bfzbdy ( ixx, iyy, izz, 1)
          u_im1 = bfzbdy ( ixx, iyy, izzm1, 1)
          u_im2 = bfzbdy ( ixx, iyy, izzm2, 1)
          u_im3 = bfzbdy ( ixx, iyy, izzm3, 1)

          R1_0 = ( 19.0 * u_im1 - 5.0 * u_im2 + u_im3
     1           + 9.0 * u_i) * ob_24
          R1_x = ( -57.0 * u_im1 - 3.0 * u_im2 + u_im3
     1           + 59.0 * u_i) * ob_60
          R1_x2 = 0.25 * ( -7.0 * u_im1 + 5.0 * u_im2 - u_im3
     1           + 3.0 * u_i)
          R1_x3 = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i) * ob_6
          Is1 = ( R1_x + 0.1 * R1_x3)**2 + r13_3 * R1_x2**2
     1           + r781_20 * R1_x3**2

          R2_0 = ( 13.0 * u_im1 - u_im2 + 13.0 * u_i - u_ip1) * ob_24
          R2_x = ( -63.0 * u_im1 + u_im2 + 63.0 * u_i - u_ip1) * ob_60
          R2_x2 = 0.25 * ( - u_im1 + u_im2 - u_i + u_ip1)
          R2_x3 = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1) * ob_6
          Is2 = ( R2_x + 0.1 * R2_x3)**2 + r13_3 * R2_x2**2
     1           + r781_20 * R2_x3**2

          R3_0 = ( 9.0 * u_im1 + 19.0 * u_i - 5.0 * u_ip1
     1           + u_ip2) * ob_24
          R3_x = ( -59.0 * u_im1 + 57.0 * u_i + 3.0 * u_ip1
     1           - u_ip2) * ob_60
          R3_x2 = 0.25 * ( 3.0 * u_im1 - 7.0 * u_i + 5.0 * u_ip1
     1           - u_ip2)
          R3_x3 = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2) * ob_6
          Is3 = ( R3_x + 0.1 * R3_x3)**2 + r13_3 * R3_x2**2
     1           + r781_20 * R3_x3**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4
          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          u ( ixx, iyy, izz, 8, 1) = Wt1 * R1_0 + Wt2 * R2_0
     1                             + Wt3 * R3_0
          u ( ixx, iyy, izz, 8, 4) = Wt1 * R1_x + Wt2 * R2_x
     1                             + Wt3 * R3_x
          u ( ixx, iyy, izz, 8, 7) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                             + Wt3 * R3_x2
          u ( ixx, iyy, izz, 8, 13) = Wt1 * R1_x3 + Wt2 * R2_x3
     1                             + Wt3 * R3_x3
          bfzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 8, 1)

! -----
#endif
! End of "NUMORDER == 4".
! ----------

          b = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          bprs ( ixx) = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          rho ( ixx) = rhogr ( ixx, iyy, izz)

          eng ( ixx) = ( u ( ixx, iyy, izz, 5, 1) - b - bprs ( ixx) )

#if ( NFLUID >= 1)
          specie ( 1:NFLUID, ixx) =
     1         u ( ixx, iyy, izz, 9: NFIELD + 1 + NFLUID, 1) * a
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

          a = 1.0 / rhogr ( ixx, iyy, izz)

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
!     1                      * prs ( ixx) ) * a        ! full msonic speed
     1     gamma_soundspeed ( ixx) * prs ( ixx) * a    ! sound speed
          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! Make the undivided divergence of the velocity, "divvelgr". Use it along
! with the minimal value of "msonicgr" in a halo of zones around the
! current one to make the flattener, "flattengr".

! "flattengr" is unity at strong shocks; zero in smooth flow.

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

          dx = dx_save ( ixx)
          dy = dy_save ( iyy)
          dz = dz_save ( izz)

          ob_dx = 1.0 / dx
          ob_dy = 1.0 / dy
          ob_dz = 1.0 / dz

          END IF

          IF ( igeom .EQ. 2) THEN

          dx = dx_save ( ixx)
          dy = xc_save ( ixx) * dy_save ( iyy)
          dz = dz_save ( izz)

          ob_dx = 1.0 / dx
          ob_dy = 1.0 / dy
          ob_dz = 1.0 / dz

          END IF

          IF ( igeom .EQ. 3) THEN

          dx = dx_save ( ixx)
          dy = xc_save ( ixx) * dy_save ( iyy)
          dz = xc_save ( ixx) * SIN ( yc_save ( iyy)) * dz_save ( izz)

          ob_dx = 1.0 / dx
          ob_dy = 1.0 / dy
          ob_dz = 1.0 / dz

          END IF

! ----------

          divvelgr ( ixx, iyy, izz) =
     1    ( vlxgr ( ixxp1, iyy, izz) - vlxgr ( ixxm1, iyy, izz) )
     1  * ob_dx
     1  + ( vlygr ( ixx, iyyp1, izz) - vlygr ( ixx, iyym1, izz) )
     1  * ob_dy
     1  + ( vlzgr ( ixx, iyy, izzp1) - vlzgr ( ixx, iyy, izzm1) )
     1  * ob_dz

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)
     1  * AMIN1 ( dx, dy, dz )

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

          a = ABS ( divvelgr ( ixx, iyy, izz))
     1      / ( flatten_coef * workgr ( ixx, iyy, izz))

          a = a - 1.0

          IF ( a .LE. 0.0) THEN
            flattengr ( ixx, iyy, izz) = 0.0
          ELSE
            flattengr ( ixx, iyy, izz) = AMIN1 ( 1.0, a)
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
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

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
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

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
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

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
! Limiting in the 1d conserved and/or characteristic variables
! done in the below three loop nests. The slopes are fed in just in
! case we wish to make a P1PM scheme in the future.
! ----------------------------------------------------------------------

! Limiting in the x-direction.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,            &
!$OMP&    bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma, smallnum, prsfloor, rhofloor, pi, gasconst,            &
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    limit_characteristics,                                        &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, flattengr)                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! ----------
          IF ( ioffx .GT. 0) THEN
! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx

! -----

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

! -----

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

! -----

          rhobar ( ixx) = rhogr ( ixx, iyy, izz)
          prsbar ( ixx) = prsgr ( ixx, iyy, izz)

          vlxbar ( ixx) = vlxgr ( ixx, iyy, izz)
          vlybar ( ixx) = vlygr ( ixx, iyy, izz)
          vlzbar ( ixx) = vlzgr ( ixx, iyy, izz)

          bfxbar ( ixx) = bfxgr ( ixx, iyy, izz)
          bfybar ( ixx) = bfygr ( ixx, iyy, izz)
          bfzbar ( ixx) = bfzgr ( ixx, iyy, izz)

#if ( NFLUID >= 1)
          a = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( ixx, 1: NFLUID) =
     1                  u ( ixx, iyy, izz, 9: n_cc_components, 1) * a
#endif

          flatten_1d ( ixx) = flattengr ( ixx, iyy, izz)

! -----

          END DO ! End of DO ixx = loop.

! ----------

          bfxbdy_1d ( 0 - ioffx: ix1 + ioffx) =
     1    bfxbdy ( 0 - ioffx: ix1 + ioffx, iyy, izz, 1)

          CALL LIMIT_1D_MHD_CONS
     1   ( 1, ix1, ioffx, limit_characteristics,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx

! -----

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

! -----
#if ( NUMORDER >= 3)

          u ( ixx, iyy, izz, 1, 5) = u_1d ( ixx, 1, 3)
          u ( ixx, iyy, izz, 2, 5) = u_1d ( ixx, 2, 3)
          u ( ixx, iyy, izz, 3, 5) = u_1d ( ixx, 3, 3)
          u ( ixx, iyy, izz, 4, 5) = u_1d ( ixx, 4, 3)
          u ( ixx, iyy, izz, 5, 5) = u_1d ( ixx, 5, 3)
          u ( ixx, iyy, izz, 7, 5) = u_1d ( ixx, 6, 3)
          u ( ixx, iyy, izz, 8, 5) = u_1d ( ixx, 7, 3)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 5) =
     1                       u_1d ( ixx, 8: NFIELD + NFLUID, 3)
#endif

# endif
! End of "NUMORDER >= 3"
! -----
#if ( NUMORDER >= 4)

          u ( ixx, iyy, izz, 1, 11) = u_1d ( ixx, 1, 4)
          u ( ixx, iyy, izz, 2, 11) = u_1d ( ixx, 2, 4)
          u ( ixx, iyy, izz, 3, 11) = u_1d ( ixx, 3, 4)
          u ( ixx, iyy, izz, 4, 11) = u_1d ( ixx, 4, 4)
          u ( ixx, iyy, izz, 5, 11) = u_1d ( ixx, 5, 4)
          u ( ixx, iyy, izz, 7, 11) = u_1d ( ixx, 6, 4)
          u ( ixx, iyy, izz, 8, 11) = u_1d ( ixx, 7, 4)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 11) =
     1                       u_1d ( ixx, 8: NFIELD + NFLUID, 4)
#endif

# endif
! End of "NUMORDER >= 4"
! -----

          END DO ! End of DO ixx = loop.

! ----------
          ELSE
! ----------

          u ( :, iyy, izz, :, 2) = 0.0
#if ( NUMORDER >= 3)
          u ( :, iyy, izz, :, 5) = 0.0
# endif
#if ( NUMORDER >= 4)
          u ( :, iyy, izz, :, 11) = 0.0
# endif

! ----------
          END IF ! End of "ioffx .GT. 0".
! ----------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Limiting in the y-direction.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,            &
!$OMP&    bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma, smallnum, prsfloor, rhofloor, pi, gasconst,            &
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    limit_characteristics,                                        &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, flattengr)                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

        DO ixx = 1 - ioffx, ix1 + ioffx

! ----------
          IF ( ioffy .GT. 0) THEN
! ----------

          DO iyy = 1 - ioffy, iy1 + ioffy

! -----

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

! -----

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

! -----

          rhobar ( iyy) = rhogr ( ixx, iyy, izz)
          prsbar ( iyy) = prsgr ( ixx, iyy, izz)

          vlxbar ( iyy) = vlygr ( ixx, iyy, izz)
          vlybar ( iyy) = vlzgr ( ixx, iyy, izz)
          vlzbar ( iyy) = vlxgr ( ixx, iyy, izz)

          bfxbar ( iyy) = bfygr ( ixx, iyy, izz)
          bfybar ( iyy) = bfzgr ( ixx, iyy, izz)
          bfzbar ( iyy) = bfxgr ( ixx, iyy, izz)

#if ( NFLUID >= 1)
          a = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( iyy, 1: NFLUID) =
     1                  u ( ixx, iyy, izz, 9: n_cc_components, 1) * a
#endif

          flatten_1d ( iyy) = flattengr ( ixx, iyy, izz)

! -----

          END DO ! End of DO iyy = loop.

! ----------

          bfxbdy_1d ( 0 - ioffy: iy1 + ioffy) =
     1    bfybdy ( ixx, 0 - ioffy: iy1 + ioffy, izz, 1)

          CALL LIMIT_1D_MHD_CONS
     1   ( 1, iy1, ioffy, limit_characteristics,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------

          DO iyy = 1 - ioffy, iy1 + ioffy

! -----

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

! -----
#if ( NUMORDER >= 3)

          u ( ixx, iyy, izz, 1, 6) = u_1d ( iyy, 1, 3)
          u ( ixx, iyy, izz, 3, 6) = u_1d ( iyy, 2, 3)
          u ( ixx, iyy, izz, 4, 6) = u_1d ( iyy, 3, 3)
          u ( ixx, iyy, izz, 2, 6) = u_1d ( iyy, 4, 3)
          u ( ixx, iyy, izz, 5, 6) = u_1d ( iyy, 5, 3)
          u ( ixx, iyy, izz, 8, 6) = u_1d ( iyy, 6, 3)
          u ( ixx, iyy, izz, 6, 6) = u_1d ( iyy, 7, 3)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 6) =
     1                       u_1d ( iyy, 8: NFIELD + NFLUID, 3)
#endif

#endif
! End of "NUMORDER >= 3"
! -----
#if ( NUMORDER >= 4)

          u ( ixx, iyy, izz, 1, 12) = u_1d ( iyy, 1, 4)
          u ( ixx, iyy, izz, 3, 12) = u_1d ( iyy, 2, 4)
          u ( ixx, iyy, izz, 4, 12) = u_1d ( iyy, 3, 4)
          u ( ixx, iyy, izz, 2, 12) = u_1d ( iyy, 4, 4)
          u ( ixx, iyy, izz, 5, 12) = u_1d ( iyy, 5, 4)
          u ( ixx, iyy, izz, 8, 12) = u_1d ( iyy, 6, 4)
          u ( ixx, iyy, izz, 6, 12) = u_1d ( iyy, 7, 4)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 12) =
     1                       u_1d ( iyy, 8: NFIELD + NFLUID, 4)
#endif

#endif
! End of "NUMORDER >= 4"
! -----

          END DO ! End of DO iyy = loop.

! ----------
          ELSE
! ----------

          u ( ixx, :, izz, :, 3) = 0
#if ( NUMORDER >= 3)
          u ( ixx, :, izz, :, 6) = 0.0
# endif
#if ( NUMORDER >= 4)
          u ( ixx, :, izz, :, 12) = 0.0
# endif

! ----------
          END IF ! End of "ioffy .GT. 0".
! ----------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Limiting in the z-direction.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,            &
!$OMP&    bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)          &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma, smallnum, prsfloor, rhofloor, pi, gasconst,            &
!$OMP&    cp_specie, molewt_specie, heat_of_formation,                  &
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    limit_characteristics,                                        &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,    &
!$OMP&    bfxbdy, bfybdy, bfzbdy, flattengr)                            &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO iyy = 1 - ioffy, iy1 + ioffy

! ----------
          IF ( ioffz .GT. 0) THEN
! ----------

          DO izz = 1 - ioffz, iz1 + ioffz

! -----

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

! -----

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

! -----

          rhobar ( izz) = rhogr ( ixx, iyy, izz)
          prsbar ( izz) = prsgr ( ixx, iyy, izz)

          vlxbar ( izz) = vlzgr ( ixx, iyy, izz)
          vlybar ( izz) = vlxgr ( ixx, iyy, izz)
          vlzbar ( izz) = vlygr ( ixx, iyy, izz)

          bfxbar ( izz) = bfzgr ( ixx, iyy, izz)
          bfybar ( izz) = bfxgr ( ixx, iyy, izz)
          bfzbar ( izz) = bfygr ( ixx, iyy, izz)

#if ( NFLUID >= 1)
          a = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( izz, 1: NFLUID) =
     1                  u ( ixx, iyy, izz, 9: n_cc_components, 1) * a
#endif

          flatten_1d ( izz) = flattengr ( ixx, iyy, izz)

! -----

          END DO ! End of DO izz = loop.

! ----------

          bfxbdy_1d ( 0 - ioffz: iz1 + ioffz) =
     1    bfzbdy ( ixx, iyy, 0 - ioffz: iz1 + ioffz, 1)

          CALL LIMIT_1D_MHD_CONS
     1   ( 1, iz1, ioffz, limit_characteristics,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------

          DO izz = 1 - ioffz, iz1 + ioffz

! -----

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

! -----
#if ( NUMORDER >= 3)

          u ( ixx, iyy, izz, 1, 7) = u_1d ( izz, 1, 3)
          u ( ixx, iyy, izz, 4, 7) = u_1d ( izz, 2, 3)
          u ( ixx, iyy, izz, 2, 7) = u_1d ( izz, 3, 3)
          u ( ixx, iyy, izz, 3, 7) = u_1d ( izz, 4, 3)
          u ( ixx, iyy, izz, 5, 7) = u_1d ( izz, 5, 3)
          u ( ixx, iyy, izz, 6, 7) = u_1d ( izz, 6, 3)
          u ( ixx, iyy, izz, 7, 7) = u_1d ( izz, 7, 3)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 7) =
     1                       u_1d ( izz, 8: NFIELD + NFLUID, 3)
#endif

#endif
! End of "NUMORDER >= 3"
! -----
#if ( NUMORDER >= 4)

          u ( ixx, iyy, izz, 1, 13) = u_1d ( izz, 1, 4)
          u ( ixx, iyy, izz, 4, 13) = u_1d ( izz, 2, 4)
          u ( ixx, iyy, izz, 2, 13) = u_1d ( izz, 3, 4)
          u ( ixx, iyy, izz, 3, 13) = u_1d ( izz, 4, 4)
          u ( ixx, iyy, izz, 5, 13) = u_1d ( izz, 5, 4)
          u ( ixx, iyy, izz, 6, 13) = u_1d ( izz, 6, 4)
          u ( ixx, iyy, izz, 7, 13) = u_1d ( izz, 7, 4)
#if ( NFLUID >= 1)
          u ( ixx, iyy, izz, 9: n_cc_components, 13) =
     1                       u_1d ( izz, 8: NFIELD + NFLUID, 4)
#endif

#endif
! End of "NUMORDER >= 4"
! -----

          END DO ! End of DO izz = loop.

! ----------
          ELSE
! ----------

          u ( ixx, iyy, :, :, 4) = 0.0
#if ( NUMORDER >= 3)
          u ( ixx, iyy, :, :, 7) = 0.0
#endif
#if ( NUMORDER >= 4)
          u ( ixx, iyy, :, :, 13) = 0.0
#endif

! ----------
          END IF ! End of "ioffz .GT. 0".
! ----------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO ixx = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! End of limiting on the 1d conserved and/or characteristic variables.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Limit the cross-terms at 3rd order.
#if ( NUMORDER == 3)
! ----------------------------------------------------------------------

! Make zone-centered cross-terms here.

      DO i_comp = 1, n_cc_components

! ------------------------------

! Performing WENO reconstruction of the cross-terms in the xy-plane.
! This is the 8th moment with basis functions P1[x] * P1[y].

      IF ( ( ioffx .GT. 1) .AND. ( ioffy .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, R1_xy, R2_xy, R3_xy, R4_xy,          &
!$OMP&    Is1, Is2, Is3, Is4, Wt1, Wt2, Wt3, Wt4)                       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, u, flattengr)                                       &
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

! -----
! This is a smaller part of a smoothness measure for the xy cross term

            a = 4.0 * ( u ( ixx, iyy, izz, i_comp, 5)**2
     1                + u ( ixx, iyy, izz, i_comp, 6)**2)
     1        + smallnum

! -----

            R1_xy = u ( ixxp1, iyyp1, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 6)

            Is1 = R1_xy**2 + a

! -----

            R2_xy = - u ( ixxp1, iyym1, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 6)

            Is2 = R2_xy**2 + a

! -----

            R3_xy = - u ( ixxm1, iyyp1, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 6)

            Is3 = R3_xy**2 + a

! -----

            R4_xy = u ( ixxm1, iyym1, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 6)

            Is4 = R4_xy**2 + a

! -----

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4

            b = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)

            Wt1 = Wt1 * b
            Wt2 = Wt2 * b
            Wt3 = Wt3 * b
            Wt4 = Wt4 * b

            u ( ixx, iyy, izz, i_comp, 8) = ( Wt1 * R1_xy + Wt2 * R2_xy
     1                                      + Wt3 * R3_xy + Wt4 * R4_xy)
     1                           * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 8) = 0.0

! ----------
      END IF
! End of "( ioffx .GT. 1) .AND. ( ioffy .GT. 1)"

! ------------------------------

! Performing WENO reconstruction of the cross-terms in the yz-plane.
! This is the 9th moment with basis functions P1[y] * P1[z].

      IF ( ( ioffy .GT. 1) .AND. ( ioffz .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, R1_yz, R2_yz, R3_yz, R4_yz,          &
!$OMP&    Is1, Is2, Is3, Is4, Wt1, Wt2, Wt3, Wt4)                       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, u, flattengr)                                       &
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

! -----
! This is a smaller part of a smoothness measure for the yz cross term

            a = 4.0 * ( u ( ixx, iyy, izz, i_comp, 6)**2
     1                + u ( ixx, iyy, izz, i_comp, 7)**2)
     1        + smallnum

! -----

            R1_yz = u ( ixx, iyyp1, izzp1, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 6)
     1         - u ( ixx, iyy, izz, i_comp, 7)

            Is1 = R1_yz**2 + a

! -----

            R2_yz = - u ( ixx, iyyp1, izzm1, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 6)
     1         + u ( ixx, iyy, izz, i_comp, 7)

            Is2 = R2_yz**2 + a

! -----

            R3_yz = - u ( ixx, iyym1, izzp1, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 6)
     1         + u ( ixx, iyy, izz, i_comp, 7)

            Is3 = R3_yz**2 + a

! -----

            R4_yz = u ( ixx, iyym1, izzm1, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 6)
     1         - u ( ixx, iyy, izz, i_comp, 7)

            Is4 = R4_yz**2 + a

! -----

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4

            b = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)

            Wt1 = Wt1 * b
            Wt2 = Wt2 * b
            Wt3 = Wt3 * b
            Wt4 = Wt4 * b

            u ( ixx, iyy, izz, i_comp, 9) = ( Wt1 * R1_yz + Wt2 * R2_yz
     1                                      + Wt3 * R3_yz + Wt4 * R4_yz)
     1                           * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 9) = 0.0

! ----------
      END IF
! End of "( ioffy .GT. 1) .AND. ( ioffz .GT. 1)"

! ------------------------------

! Performing WENO reconstruction of the cross-terms in the xz-plane.
! This is the 10th moment with basis functions P1[x] * P1[z].

      IF ( ( ioffx .GT. 1) .AND. ( ioffz .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, R1_xz, R2_xz, R3_xz, R4_xz,          &
!$OMP&    Is1, Is2, Is3, Is4, Wt1, Wt2, Wt3, Wt4)                       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, u, flattengr)                                       &
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

! -----
! This is a smaller part of a smoothness measure for the xz cross term

            a = 4.0 * ( u ( ixx, iyy, izz, i_comp, 5)**2
     1                + u ( ixx, iyy, izz, i_comp, 7)**2)
     1        + smallnum

! -----

            R1_xz = u ( ixxp1, iyy, izzp1, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 7)

            Is1 = R1_xz**2 + a

! -----

            R2_xz = - u ( ixxp1, iyy, izzm1, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 7)

            Is2 = R2_xz**2 + a

! -----

            R3_xz = - u ( ixxm1, iyy, izzp1, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 7)

            Is3 = R3_xz**2 + a

! -----

            R4_xz = u ( ixxm1, iyy, izzm1, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 7)

            Is4 = R4_xz**2 + a

! -----

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4

            b = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)

            Wt1 = Wt1 * b
            Wt2 = Wt2 * b
            Wt3 = Wt3 * b
            Wt4 = Wt4 * b

            u ( ixx, iyy, izz, i_comp, 10) =
     1                                 ( Wt1 * R1_xz + Wt2 * R2_xz
     1                                 + Wt3 * R3_xz + Wt4 * R4_xz)
     1                           * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 10) = 0.0

! ----------
      END IF
! End of "( ioffx .GT. 1) .AND. ( ioffz .GT. 1)"

! ------------------------------

      END DO
! End of "DO i_comp ="

! ----------------------------------------------------------------------
#endif
! End of "NUMORDER == 3"
! End of limiting the cross-terms at 3rd order.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Limit the cross-terms at 4th order.
#if ( NUMORDER == 4)
! ----------------------------------------------------------------------

! Make zone-centered cross-terms here.

      DO i_comp = 1, n_cc_components

! ------------------------------

! Performing WENO reconstruction of the cross-terms in the xy-plane.
! This is the 8th, 14th and 16th moment with basis functions P1[x] * P1[y],
! P2[x] * P1[y] and P1[x] * P2[y].

      IF ( ( ioffx .GT. 1) .AND. ( ioffy .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, R1_xy, R2_xy, R3_xy, R4_xy, R5_xy,                   &
!$OMP&    R1_x2y, R2_x2y, R3_x2y, R4_x2y, R5_x2y,                       &
!$OMP&    R1_xy2, R2_xy2, R3_xy2, R4_xy2, R5_xy2,                       &
!$OMP&    Is1, Is2, Is3, Is4, Is5, Is1_4, Is2_4, Is3_4, Is4_4, Is5_4,   &
!$OMP&    Wt1, Wt2, Wt3, Wt4, Wt5)                                      &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, r2_3, central_zone_wt, u, flattengr)                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - ioffz, iz1 + ioffz

        izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
        izzm1 = MAX0 ( izz - 1, 1 - ioffz)
        izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
        izzm2 = MAX0 ( izz - 2, 1 - ioffz)

          DO iyy = 1 - ioffy, iy1 + ioffy

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

            DO ixx = 1 - ioffx, ix1 + ioffx

            ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
            ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
            ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
            ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

! -----
! This is a smaller part of a smoothness measure for the xy cross term

            a = 3.0 * ( u ( ixx, iyy, izz, i_comp, 11)**2
     1                + u ( ixx, iyy, izz, i_comp, 12)**2)
     1        + 4.0 * ( u ( ixx, iyy, izz, i_comp, 5)**2
     1                + u ( ixx, iyy, izz, i_comp, 6)**2)
     1        + smallnum

! -----
! This is a smaller part of a smoothness measure for the x^2y and y^2x cross terms

            b = 36.0 *( u ( ixx, iyy, izz, i_comp, 11)**2
     1                + u ( ixx, iyy, izz, i_comp, 12)**2 )
     1        + smallnum

! -----

            R1_xy = 3.0*u ( ixxp1, iyyp1, izz, i_comp, 1)
     1          - 0.5*( u ( ixxp1, iyyp2, izz, i_comp, 1)
     1                 +u ( ixxp2, iyyp1, izz, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 5)
     1                 +u ( ixx, iyy, izz, i_comp, 6))
     1          - 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1          + 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12))
     1          - 1.5*( u ( ixx, iyy, izz, i_comp, 2)
     1                 +u ( ixx, iyy, izz, i_comp, 3)) 


            R1_x2y = - u ( ixxp1, iyyp1, izz, i_comp, 1)
     1            + 0.5*( u ( ixxp2, iyyp1, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            - 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1            - u ( ixx, iyy, izz, i_comp, 5)


            R1_xy2 = - u ( ixxp1, iyyp1, izz, i_comp, 1)
     1            + 0.5*( u ( ixxp1, iyyp2, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 5))
     1            - 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1            - u ( ixx, iyy, izz, i_comp, 6)

            Is1 = R1_xy**2 + r2_3 * ( R1_x2y**2 + R1_xy2**2) + a

            Is1_4 = 4.0 * ( R1_x2y**2 + R1_xy2**2 ) + b

! -----

            R2_xy = 3.0*u ( ixxm1, iyym1, izz, i_comp, 1)
     1          - 0.5*( u ( ixxm1, iyym2, izz, i_comp, 1)
     1                 +u ( ixxm2, iyym1, izz, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 5)
     1                 +u ( ixx, iyy, izz, i_comp, 6))
     1          - 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1          - 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12))
     1          + 1.5*( u ( ixx, iyy, izz, i_comp, 2)
     1                 +u ( ixx, iyy, izz, i_comp, 3))


            R2_x2y = u ( ixxm1, iyym1, izz, i_comp, 1)
     1          - 0.5*( u ( ixxm2, iyym1, izz, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 1)
     1                 -u ( ixx, iyy, izz, i_comp, 3)
     1                 +u ( ixx, iyy, izz, i_comp, 6))
     1           - 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1           + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1           + u ( ixx, iyy, izz, i_comp, 5)


            R2_xy2 = u ( ixxm1, iyym1, izz, i_comp, 1)
     1          - 0.5*( u ( ixxm1, iyym2, izz, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 1)
     1                 -u ( ixx, iyy, izz, i_comp, 2)
     1                 +u ( ixx, iyy, izz, i_comp, 5))
     1           - 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1           + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1           + u ( ixx, iyy, izz, i_comp, 6)

            Is2 = R2_xy**2 + r2_3 * (R2_x2y**2 + R2_xy2**2) + a

            Is2_4 = 4.0 * ( R2_x2y**2 + R2_xy2**2 ) + b

! -----

            R3_xy = - 3.0*u ( ixxp1, iyym1, izz, i_comp, 1)
     1            + 0.5*( u ( ixxp1, iyym2, izz, i_comp, 1)
     1                   +u ( ixxp2, iyym1, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 5)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            + 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1            - 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                       -u ( ixx, iyy, izz, i_comp, 12))
     1            + 1.5*( u ( ixx, iyy, izz, i_comp, 2)
     1                   -u ( ixx, iyy, izz, i_comp, 3))


            R3_x2y =   u ( ixxp1, iyym1, izz, i_comp, 1)
     1            - 0.5*( u ( ixxp2, iyym1, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   -u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1            + u ( ixx, iyy, izz, i_comp, 5)


            R3_xy2 = - u ( ixxp1, iyym1, izz, i_comp, 1)
     1            + 0.5*( u ( ixxp1, iyym2, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 5))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1            - u ( ixx, iyy, izz, i_comp, 6)

            Is3 = R3_xy**2 + r2_3 * ( R3_x2y**2 + R3_xy2**2) + a

            Is3_4 = 4.0 * ( R3_x2y**2 + R3_xy2**2 ) + b

! -----

            R4_xy = - 3.0*u ( ixxm1, iyyp1, izz, i_comp, 1)
     1            + 0.5*( u ( ixxm1, iyyp2, izz, i_comp, 1)
     1                   +u ( ixxm2, iyyp1, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 5)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            + 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1            + 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                       -u ( ixx, iyy, izz, i_comp, 12))
     1            + 1.5*(-u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 3))


            R4_x2y = - u ( ixxm1, iyyp1, izz, i_comp, 1)
     1            + 0.5*( u ( ixxm2, iyyp1, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1            - u ( ixx, iyy, izz, i_comp, 5)


            R4_xy2 =   u ( ixxm1, iyyp1, izz, i_comp, 1)
     1            - 0.5*( u ( ixxm1, iyyp2, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   -u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 5))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1            + u ( ixx, iyy, izz, i_comp, 6)

            Is4 = R4_xy**2 + r2_3 * ( R4_x2y**2 + R4_xy2**2) + a

            Is4_4 = 4.0 * ( R4_x2y**2 + R4_xy2**2 ) + b

! -----

            R5_xy = 0.25 * ( u ( ixxm1, iyym1, izz, i_comp, 1)
     1                   - u ( ixxm1, iyyp1, izz, i_comp, 1)
     1                   - u ( ixxp1, iyym1, izz, i_comp, 1)
     1                   + u ( ixxp1, iyyp1, izz, i_comp, 1) )

            R5_x2y = 0.05*(- 5.0 * u ( ixxm1, iyym1, izz, i_comp, 1)
     1                  + 5.0 * u ( ixxm1, iyyp1, izz, i_comp, 1)
     1                  - 5.0 * u ( ixxp1, iyym1, izz, i_comp, 1)
     1                  + 5.0 * u ( ixxp1, iyyp1, izz, i_comp, 1)
     1                  - 22.0 * u ( ixx, iyy, izz, i_comp, 12)
     1                  - 20.0 * u ( ixx, iyy, izz, i_comp, 3 ) )

            R5_xy2 = 0.05*(- 5.0 * u ( ixxm1, iyym1, izz, i_comp, 1)
     1                  - 5.0 * u ( ixxm1, iyyp1, izz, i_comp, 1)
     1                  + 5.0 * u ( ixxp1, iyym1, izz, i_comp, 1)
     1                  + 5.0 * u ( ixxp1, iyyp1, izz, i_comp, 1)
     1                  - 22.0 * u ( ixx, iyy, izz, i_comp, 11)
     1                  - 20.0 * u ( ixx, iyy, izz, i_comp, 2) )

            Is5 = R5_xy**2 + r2_3 * ( R5_x2y**2 + R5_xy2**2) + a

            Is5_4 = 4.0 * ( R5_x2y**2 + R5_xy2**2 ) + b


! -----
! Weights for the xy cross term

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4
            Wt5 = central_zone_wt / Is5**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5)

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c

            u ( ixx, iyy, izz, i_comp, 8) = ( Wt1 * R1_xy + Wt2 * R2_xy
     1                                      + Wt3 * R3_xy + Wt4 * R4_xy
     1                                      + Wt5 * R5_xy)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----
! Weights for the x^2y and xy^2 cross terms

            Wt1 = 1.0 / Is1_4**4
            Wt2 = 1.0 / Is2_4**4
            Wt3 = 1.0 / Is3_4**4
            Wt4 = 1.0 / Is4_4**4
            Wt5 = central_zone_wt / Is5_4**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5)

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c

            u ( ixx, iyy, izz, i_comp, 14) =
     1        ( Wt1 * R1_x2y + Wt2 * R2_x2y + Wt3 * R3_x2y
     1        + Wt4 * R4_x2y + Wt5 * R5_x2y)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

            u ( ixx, iyy, izz, i_comp, 16) =
     1        ( Wt1 * R1_xy2 + Wt2 * R2_xy2 + Wt3 * R3_xy2
     1        + Wt4 * R4_xy2 + Wt5 * R5_xy2)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 8) = 0.0
        u ( :, :, :, i_comp, 14) = 0.0
        u ( :, :, :, i_comp, 16) = 0.0

! ----------
      END IF
! End of "( ioffx .GT. 1) .AND. ( ioffy .GT. 1)"

! ------------------------------

! Performing WENO reconstruction of the cross-terms in the yz-plane.
! This is the 9th, 17th and 19th moment with basis functions P1[y] * P1[z],
! P2[y] * P1[z] and P1[y] * P2[z].

      IF ( ( ioffy .GT. 1) .AND. ( ioffz .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, R1_yz, R2_yz, R3_yz, R4_yz, R5_yz,                   &
!$OMP&    R1_y2z, R2_y2z, R3_y2z, R4_y2z, R5_y2z,                       &
!$OMP&    R1_yz2, R2_yz2, R3_yz2, R4_yz2, R5_yz2,                       &
!$OMP&    Is1, Is2, Is3, Is4, Is5, Is1_4, Is2_4, Is3_4, Is4_4, Is5_4,   &
!$OMP&    Wt1, Wt2, Wt3, Wt4, Wt5)                                      &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, r2_3, central_zone_wt, u, flattengr)                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - ioffz, iz1 + ioffz

        izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
        izzm1 = MAX0 ( izz - 1, 1 - ioffz)
        izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
        izzm2 = MAX0 ( izz - 2, 1 - ioffz)

          DO iyy = 1 - ioffy, iy1 + ioffy

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

            DO ixx = 1 - ioffx, ix1 + ioffx

            ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
            ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
            ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
            ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

! -----
! This is a smaller part of a smoothness measure for the yz term

            a = 3.0 * ( u ( ixx, iyy, izz, i_comp, 12)**2
     1                + u ( ixx, iyy, izz, i_comp, 13)**2)
     1        + 4.0 * ( u ( ixx, iyy, izz, i_comp, 6)**2
     1                + u ( ixx, iyy, izz, i_comp, 7)**2)
     1        + smallnum

! -----
! This is a smaller part of a smoothness measure for the y^2z and yz^2 cross terms

            b = 36.0 * ( u ( ixx, iyy, izz, i_comp, 12)**2
     1                 + u ( ixx, iyy, izz, i_comp, 13)**2 )
     1        + smallnum

! -----

            R1_yz = 3.0*u ( ixx, iyyp1, izzp1, i_comp, 1)
     1          - 0.5*( u ( ixx, iyyp1, izzp2, i_comp, 1)
     1                 +u ( ixx, iyyp2, izzp1, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 6)
     1                 +u ( ixx, iyy, izz, i_comp, 7))
     1          - 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1          + 1.35 * ( u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1          - 1.5*( u ( ixx, iyy, izz, i_comp, 3)
     1                 +u ( ixx, iyy, izz, i_comp, 4))


            R1_y2z = - u ( ixx, iyyp1, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixx, iyyp2, izzp1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 4)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            - 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1            - u ( ixx, iyy, izz, i_comp, 6)


            R1_yz2 = - u ( ixx, iyyp1, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixx, iyyp1, izzp2, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            - 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1            - u ( ixx, iyy, izz, i_comp, 7)

            Is1 = R1_yz**2 + r2_3 * ( R1_y2z**2 + R1_yz2**2) + a

            Is1_4 = 4.0 * ( R1_y2z**2 + R1_yz2**2 ) + b

! -----

            R2_yz = 3.0*u ( ixx, iyym1, izzm1, i_comp, 1)
     1          - 0.5*( u ( ixx, iyym1, izzm2, i_comp, 1)
     1                 +u ( ixx, iyym2, izzm1, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 6)
     1                 +u ( ixx, iyy, izz, i_comp, 7))
     1          - 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1          - 1.35 * ( u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1          + 1.5*( u ( ixx, iyy, izz, i_comp, 3)
     1                 +u ( ixx, iyy, izz, i_comp, 4))


            R2_y2z = u ( ixx, iyym1, izzm1, i_comp, 1)
     1          - 0.5*( u ( ixx, iyym2, izzm1, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 1)
     1                 -u ( ixx, iyy, izz, i_comp, 4)
     1                 +u ( ixx, iyy, izz, i_comp, 7))
     1           - 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1           + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1           + u ( ixx, iyy, izz, i_comp, 6)


            R2_yz2 = u ( ixx, iyym1, izzm1, i_comp, 1)
     1          - 0.5*( u ( ixx, iyym1, izzm2, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 1)
     1                 -u ( ixx, iyy, izz, i_comp, 3)
     1                 +u ( ixx, iyy, izz, i_comp, 6))
     1           - 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1           + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1           + u ( ixx, iyy, izz, i_comp, 7)

            Is2 = R2_yz**2 + r2_3 * ( R2_y2z**2 + R2_yz2**2) + a

            Is2_4 = 4.0 * ( R2_y2z**2 + R2_yz2**2 ) + b

! -----

            R3_yz = - 3.0*u ( ixx, iyyp1, izzm1, i_comp, 1)
     1            + 0.5*( u ( ixx, iyyp1, izzm2, i_comp, 1)
     1                   +u ( ixx, iyyp2, izzm1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 6)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1            - 1.35 * ( u ( ixx, iyy, izz, i_comp, 12)
     1                       -u ( ixx, iyy, izz, i_comp, 13))
     1            + 1.5*( u ( ixx, iyy, izz, i_comp, 3)
     1                   -u ( ixx, iyy, izz, i_comp, 4))


            R3_y2z =   u ( ixx, iyyp1, izzm1, i_comp, 1)
     1            - 0.5*( u ( ixx, iyyp2, izzm1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   -u ( ixx, iyy, izz, i_comp, 4)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1            + u ( ixx, iyy, izz, i_comp, 6)


            R3_yz2 = - u ( ixx, iyyp1, izzm1, i_comp, 1)
     1            + 0.5*( u ( ixx, iyyp1, izzm2, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1            - u ( ixx, iyy, izz, i_comp, 7)

            Is3 = R3_yz**2 + r2_3 * ( R3_y2z**2 + R3_yz2**2) + a

            Is3_4 = 4.0 * ( R3_y2z**2 + R3_yz2**2 ) + b

! -----

            R4_yz = - 3.0*u ( ixx, iyym1, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixx, iyym1, izzp2, i_comp, 1)
     1                   +u ( ixx, iyym2, izzp1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 6)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1            + 1.35 * ( u ( ixx, iyy, izz, i_comp, 12)
     1                       -u ( ixx, iyy, izz, i_comp, 13))
     1            + 1.5*(-u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 4))


            R4_y2z = - u ( ixx, iyym1, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixx, iyym2, izzp1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 4)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 12)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1            - u ( ixx, iyy, izz, i_comp, 6)


            R4_yz2 =   u ( ixx, iyym1, izzp1, i_comp, 1)
     1            - 0.5*( u ( ixx, iyym1, izzp2, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   -u ( ixx, iyy, izz, i_comp, 3)
     1                   +u ( ixx, iyy, izz, i_comp, 6))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 12)
     1            + u ( ixx, iyy, izz, i_comp, 7)

            Is4 = R4_yz**2 + r2_3 * ( R4_y2z**2 + R4_yz2**2) + a

            Is4_4 = 4.0 * ( R4_y2z**2 + R4_yz2**2 ) + b

! -----

           R5_yz = 0.25 * ( u ( ixx, iyym1, izzm1, i_comp, 1)
     1                  - u ( ixx, iyym1, izzp1, i_comp, 1)
     1                  - u ( ixx, iyyp1, izzm1, i_comp, 1)
     1                  + u ( ixx, iyyp1, izzp1, i_comp, 1) )

           R5_y2z = 0.05 *(- 5.0 * u ( ixx, iyym1, izzm1, i_comp, 1)
     1                  + 5.0 * u ( ixx, iyym1, izzp1, i_comp, 1)
     1                  - 5.0 * u ( ixx, iyyp1, izzm1, i_comp, 1)
     1                  + 5.0 * u ( ixx, iyyp1, izzp1, i_comp, 1)
     1                  - 22.0 * u ( ixx, iyy, izz, i_comp, 13)
     1                  - 20.0 * u ( ixx, iyy, izz, i_comp, 4 ) )

           R5_yz2 = 0.05 *(- 5.0 * u ( ixx, iyym1, izzm1, i_comp, 1)
     1                  - 5.0 * u ( ixx, iyym1, izzp1, i_comp, 1)
     1                  + 5.0 * u ( ixx, iyyp1, izzm1, i_comp, 1)
     1                  + 5.0 * u ( ixx, iyyp1, izzp1, i_comp, 1)
     1                  - 22.0 * u ( ixx, iyy, izz, i_comp, 12)
     1                  - 20.0 * u ( ixx, iyy, izz, i_comp, 3) )

            Is5 = R5_yz**2 + r2_3 * ( R5_y2z**2 + R5_yz2**2) + a

            Is5_4 = 4.0 * ( R5_y2z**2 + R5_yz2**2 ) + b

! -----
! Weights for the yz cross term

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4
            Wt5 = central_zone_wt / Is5**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5 )

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c

            u ( ixx, iyy, izz, i_comp, 9) = ( Wt1 * R1_yz + Wt2 * R2_yz
     1                                      + Wt3 * R3_yz + Wt4 * R4_yz
     1                                      + Wt5 * R5_yz)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----
! Weights for the y^2z and yz^2 cross terms

            Wt1 = 1.0 / Is1_4**4
            Wt2 = 1.0 / Is2_4**4
            Wt3 = 1.0 / Is3_4**4
            Wt4 = 1.0 / Is4_4**4
            Wt5 = central_zone_wt / Is5_4**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5)

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c

            u ( ixx, iyy, izz, i_comp, 17) =
     1        ( Wt1 * R1_y2z + Wt2 * R2_y2z + Wt3 * R3_y2z
     1        + Wt4 * R4_y2z + Wt5 * R5_y2z)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

            u ( ixx, iyy, izz, i_comp, 19) =
     1        ( Wt1 * R1_yz2 + Wt2 * R2_yz2 + Wt3 * R3_yz2
     1        + Wt4 * R4_yz2 + Wt5 * R5_yz2)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 9) = 0.0
        u ( :, :, :, i_comp, 17) = 0.0
        u ( :, :, :, i_comp, 19) = 0.0

! ----------
      END IF
! End of "( ioffy .GT. 1) .AND. ( ioffz .GT. 1)"

! ------------------------------

! Performing WENO reconstruction of the cross-terms in the xz-plane.
! This is the 10th, 15th and 18th moment with basis functions P1[x] * P1[z],
! P2[x] * P1[z] and P1[x] * P2[z].

      IF ( ( ioffx .GT. 1) .AND. ( ioffz .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, R1_xz, R2_xz, R3_xz, R4_xz, R5_xz,                   &
!$OMP&    R1_x2z, R2_x2z, R3_x2z, R4_x2z, R5_x2z,                       &
!$OMP&    R1_xz2, R2_xz2, R3_xz2, R4_xz2, R5_xz2,                       &
!$OMP&    Is1, Is2, Is3, Is4, Is5, Is1_4, Is2_4, Is3_4, Is4_4, Is5_4,   &
!$OMP&    Wt1, Wt2, Wt3, Wt4, Wt5)                                      &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, r2_3, central_zone_wt, u, flattengr)                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - ioffz, iz1 + ioffz

        izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
        izzm1 = MAX0 ( izz - 1, 1 - ioffz)
        izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
        izzm2 = MAX0 ( izz - 2, 1 - ioffz)

          DO iyy = 1 - ioffy, iy1 + ioffy

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

            DO ixx = 1 - ioffx, ix1 + ioffx

            ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
            ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
            ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
            ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

! -----
! This is a smaller part of a smoothness measure for the xz cross term

            a = 3.0 * ( u ( ixx, iyy, izz, i_comp, 11)**2
     1                     + u ( ixx, iyy, izz, i_comp, 13)**2)
     1               + 4.0 * (u ( ixx, iyy, izz, i_comp, 5)**2
     1                     + u ( ixx, iyy, izz, i_comp, 7)**2)
     1             + smallnum

! -----
! This is a smaller part of a smoothness measure for the x^2z and xz^2 cross terms

            b = 36.0 *( u ( ixx, iyy, izz, i_comp, 11)**2
     1                     + u ( ixx, iyy, izz, i_comp, 13)**2 )
     1                     + smallnum

! -----

            R1_xz = 3.0*u ( ixxp1, iyy, izzp1, i_comp, 1)
     1          - 0.5*( u ( ixxp1, iyy, izzp2, i_comp, 1)
     1                 +u ( ixxp2, iyy, izzp1, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 5)
     1                 +u ( ixx, iyy, izz, i_comp, 7))
     1          - 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1          + 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1          - 1.5*( u ( ixx, iyy, izz, i_comp, 2)
     1                 +u ( ixx, iyy, izz, i_comp, 4))


            R1_x2z = - u ( ixxp1, iyy, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixxp2, iyy, izzp1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 4)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            - 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1            - u ( ixx, iyy, izz, i_comp, 5)


            R1_xz2 = - u ( ixxp1, iyy, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixxp1, iyy, izzp2, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 5))
     1            - 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1            - u ( ixx, iyy, izz, i_comp, 7)

            Is1 = R1_xz**2 + r2_3 * ( R1_x2z**2 + R1_xz2**2) + a

            Is1_4 = 4.0 * ( R1_x2z**2 + R1_xz2**2 ) + b

! -----

            R2_xz = 3.0*u ( ixxm1, iyy, izzm1, i_comp, 1)
     1          - 0.5*( u ( ixxm1, iyy, izzm2, i_comp, 1)
     1                 +u ( ixxm2, iyy, izzm1, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 5)
     1                 +u ( ixx, iyy, izz, i_comp, 7))
     1          - 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1          - 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1          + 1.5*( u ( ixx, iyy, izz, i_comp, 2)
     1                 +u ( ixx, iyy, izz, i_comp, 4))


            R2_x2z = u ( ixxm1, iyy, izzm1, i_comp, 1)
     1          - 0.5*( u ( ixxm2, iyy, izzm1, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 1)
     1                 -u ( ixx, iyy, izz, i_comp, 4)
     1                 +u ( ixx, iyy, izz, i_comp, 7))
     1           - 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1           + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1           + u ( ixx, iyy, izz, i_comp, 5)


            R2_xz2 = u ( ixxm1, iyy, izzm1, i_comp, 1)
     1          - 0.5*( u ( ixxm1, iyy, izzm2, i_comp, 1)
     1                 +u ( ixx, iyy, izz, i_comp, 1)
     1                 -u ( ixx, iyy, izz, i_comp, 2)
     1                 +u ( ixx, iyy, izz, i_comp, 5))
     1           - 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1           + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1           + u ( ixx, iyy, izz, i_comp, 7)

            Is2 = R2_xz**2 + r2_3 * ( R2_x2z**2 + R2_xz2**2) + a

            Is2_4 = 4.0 * ( R2_x2z**2 + R2_xz2**2 ) + b

! -----

            R3_xz = - 3.0*u ( ixxp1, iyy, izzm1, i_comp, 1)
     1            + 0.5*( u ( ixxp1, iyy, izzm2, i_comp, 1)
     1                   +u ( ixxp2, iyy, izzm1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 5)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1            - 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                       -u ( ixx, iyy, izz, i_comp, 13))
     1            + 1.5*( u ( ixx, iyy, izz, i_comp, 2)
     1                   -u ( ixx, iyy, izz, i_comp, 4))


            R3_x2z =   u ( ixxp1, iyy, izzm1, i_comp, 1)
     1            - 0.5*( u ( ixxp2, iyy, izzm1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   -u ( ixx, iyy, izz, i_comp, 4)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1            + u ( ixx, iyy, izz, i_comp, 5)


            R3_xz2 = - u ( ixxp1, iyy, izzm1, i_comp, 1)
     1            + 0.5*( u ( ixxp1, iyy, izzm2, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 5))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1            - u ( ixx, iyy, izz, i_comp, 7)

            Is3 = R3_xz**2 + r2_3 * ( R3_x2z**2 + R3_xz2**2) + a

            Is3_4 = 4.0 * ( R3_x2z**2 + R3_xz2**2 ) + b

! -----

            R4_xz = - 3.0*u ( ixxm1, iyy, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixxm1, iyy, izzp2, i_comp, 1)
     1                   +u ( ixxm2, iyy, izzp1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 5)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 2.0*u ( ixx, iyy, izz, i_comp, 1)
     1            + 1.35 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                       -u ( ixx, iyy, izz, i_comp, 13))
     1            + 1.5*(-u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 4))


            R4_x2z = - u ( ixxm1, iyy, izzp1, i_comp, 1)
     1            + 0.5*( u ( ixxm2, iyy, izzp1, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 4)
     1                   +u ( ixx, iyy, izz, i_comp, 7))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 11)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 13)
     1            - u ( ixx, iyy, izz, i_comp, 5)


            R4_xz2 =   u ( ixxm1, iyy, izzp1, i_comp, 1)
     1            - 0.5*( u ( ixxm1, iyy, izzp2, i_comp, 1)
     1                   +u ( ixx, iyy, izz, i_comp, 1)
     1                   -u ( ixx, iyy, izz, i_comp, 2)
     1                   +u ( ixx, iyy, izz, i_comp, 5))
     1            + 3.0*u ( ixx, iyy, izz, i_comp, 13)
     1            + 0.55 * u ( ixx, iyy, izz, i_comp, 11)
     1            + u ( ixx, iyy, izz, i_comp, 7)

            Is4 = R4_xz**2 + r2_3 * ( R4_x2z**2 + R4_xz2**2) + a

            Is4_4 = 4.0 * ( R4_x2z**2 + R4_xz2**2 ) + b

! -----

           R5_xz = 0.25 * ( u ( ixxm1, iyy, izzm1, i_comp, 1)
     1                  - u ( ixxm1, iyy, izzp1, i_comp, 1)
     1                  - u ( ixxp1, iyy, izzm1, i_comp, 1)
     1                  + u ( ixxp1, iyy, izzp1, i_comp, 1) )

           R5_x2z = 0.05*(- 5.0 * u ( ixxm1, iyy, izzm1, i_comp, 1)
     1                 + 5.0 * u ( ixxm1, iyy, izzp1, i_comp, 1)
     1                 - 5.0 * u ( ixxp1, iyy, izzm1, i_comp, 1)
     1                 + 5.0 * u ( ixxp1, iyy, izzp1, i_comp, 1)
     1                 - 22.0 * u ( ixx, iyy, izz, i_comp, 13)
     1                 - 20.0 * u ( ixx, iyy, izz, i_comp, 4 ) )

           R5_xz2 = 0.05*(- 5.0 * u ( ixxm1, iyy, izzm1, i_comp, 1)
     1                 - 5.0 * u ( ixxm1, iyy, izzp1, i_comp, 1)
     1                 + 5.0 * u ( ixxp1, iyy, izzm1, i_comp, 1)
     1                 + 5.0 * u ( ixxp1, iyy, izzp1, i_comp, 1)
     1                 - 22.0 * u ( ixx, iyy, izz, i_comp, 11)
     1                 - 20.0 * u ( ixx, iyy, izz, i_comp, 2) )

            Is5 = R5_xz**2 + r2_3 * ( R5_x2z**2 + R5_xz2**2) + a

            Is5_4 = 4.0 * ( R5_x2z**2 + R5_xz2**2 ) + b

! -----
! Weights for the xz cross term

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4
            Wt5 = central_zone_wt / Is5**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5 )

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c

            u ( ixx, iyy, izz, i_comp, 10) =
     1                                 ( Wt1 * R1_xz + Wt2 * R2_xz
     1                                 + Wt3 * R3_xz + Wt4 * R4_xz
     1                                 + Wt5 * R5_xz)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----
! Weights for the x^2z and xz^2 cross terms

            Wt1 = 1.0 / Is1_4**4
            Wt2 = 1.0 / Is2_4**4
            Wt3 = 1.0 / Is3_4**4
            Wt4 = 1.0 / Is4_4**4
            Wt5 = central_zone_wt / Is5_4**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5)

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c

            u ( ixx, iyy, izz, i_comp, 15) =
     1        ( Wt1 * R1_x2z + Wt2 * R2_x2z + Wt3 * R3_x2z
     1        + Wt4 * R4_x2z + Wt5 * R5_x2z)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

            u ( ixx, iyy, izz, i_comp, 18) =
     1        ( Wt1 * R1_xz2 + Wt2 * R2_xz2 + Wt3 * R3_xz2
     1        + Wt4 * R4_xz2 + Wt5 * R5_xz2)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 10) = 0.0
        u ( :, :, :, i_comp, 15) = 0.0
        u ( :, :, :, i_comp, 18) = 0.0

! ----------
      END IF
! End of "( ioffx .GT. 1) .AND. ( ioffz .GT. 1)"

! ------------------------------

! Performing WENO reconstruction of the cross-terms in xyz-space.
! This is the 20th moment with basis functions P1[x]* P1[y] * P1[z].

      IF ( ( ioffx .GT. 1) .AND. ( ioffy .GT. 1) .AND.
     1     ( ioffz .GT. 1) ) THEN

! ----------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c,                                                      &
!$OMP&   R1_xyz, R2_xyz, R3_xyz, R4_xyz, R5_xyz, R6_xyz, R7_xyz, R8_xyz,&
!$OMP&    Is1, Is2, Is3, Is4, Is5, Is6, Is7, Is8,                       &
!$OMP&    Wt1, Wt2, Wt3, Wt4, Wt5, Wt6, Wt7, Wt8)                       &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, i_comp,  &
!$OMP&    smallnum, u, flattengr)                                       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - ioffz, iz1 + ioffz

        izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
        izzm1 = MAX0 ( izz - 1, 1 - ioffz)
        izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
        izzm2 = MAX0 ( izz - 2, 1 - ioffz)

          DO iyy = 1 - ioffy, iy1 + ioffy

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

            DO ixx = 1 - ioffx, ix1 + ioffx

            ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
            ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
            ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
            ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

! -----
! This is a smaller part of a smoothness measure for the third order term

            b = 36.0 *( u ( ixx, iyy, izz, i_comp, 11)**2
     1                     + u ( ixx, iyy, izz, i_comp, 12)**2
     1                     + u ( ixx, iyy, izz, i_comp, 13)**2)
     1             + 4.0 * (u ( ixx, iyy, izz, i_comp, 14)**2
     1                     + u ( ixx, iyy, izz, i_comp, 15)**2
     1                     + u ( ixx, iyy, izz, i_comp, 16)**2
     1                     + u ( ixx, iyy, izz, i_comp, 17)**2
     1                     + u ( ixx, iyy, izz, i_comp, 18)**2
     1                     + u ( ixx, iyy, izz, i_comp, 19)**2)
     1                     + smallnum

! -----

            R1_xyz = u ( ixxp1, iyyp1, izzp1, i_comp, 1)
     1         - 1.1 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 6)
     1         - u ( ixx, iyy, izz, i_comp, 7)
     1         - u ( ixx, iyy, izz, i_comp, 8)
     1         - u ( ixx, iyy, izz, i_comp, 9)
     1         - u ( ixx, iyy, izz, i_comp, 10)
     1         - u ( ixx, iyy, izz, i_comp, 14)
     1         - u ( ixx, iyy, izz, i_comp, 15)
     1         - u ( ixx, iyy, izz, i_comp, 16)
     1         - u ( ixx, iyy, izz, i_comp, 17)
     1         - u ( ixx, iyy, izz, i_comp, 18)
     1         - u ( ixx, iyy, izz, i_comp, 19)

           Is1 = R1_xyz**2 + b

! -----

            R2_xyz = - u ( ixxm1, iyyp1, izzp1, i_comp, 1)
     1         + 1.1 * (-u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 6)
     1         + u ( ixx, iyy, izz, i_comp, 7)
     1         - u ( ixx, iyy, izz, i_comp, 8)
     1         + u ( ixx, iyy, izz, i_comp, 9)
     1         - u ( ixx, iyy, izz, i_comp, 10)
     1         + u ( ixx, iyy, izz, i_comp, 14)
     1         + u ( ixx, iyy, izz, i_comp, 15)
     1         - u ( ixx, iyy, izz, i_comp, 16)
     1         + u ( ixx, iyy, izz, i_comp, 17)
     1         - u ( ixx, iyy, izz, i_comp, 18)
     1         + u ( ixx, iyy, izz, i_comp, 19)

           Is2 = R2_xyz**2 + b

! -----

            R3_xyz = - u ( ixxp1, iyym1, izzp1, i_comp, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     -u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 6)
     1         + u ( ixx, iyy, izz, i_comp, 7)
     1         - u ( ixx, iyy, izz, i_comp, 8)
     1         - u ( ixx, iyy, izz, i_comp, 9)
     1         + u ( ixx, iyy, izz, i_comp, 10)
     1         - u ( ixx, iyy, izz, i_comp, 14)
     1         + u ( ixx, iyy, izz, i_comp, 15)
     1         + u ( ixx, iyy, izz, i_comp, 16)
     1         + u ( ixx, iyy, izz, i_comp, 17)
     1         + u ( ixx, iyy, izz, i_comp, 18)
     1         - u ( ixx, iyy, izz, i_comp, 19)

           Is3 = R3_xyz**2 + b

! -----

            R4_xyz = - u ( ixxp1, iyyp1, izzm1, i_comp, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12)
     1                     -u ( ixx, iyy, izz, i_comp, 13))
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 6)
     1         + u ( ixx, iyy, izz, i_comp, 7)
     1         + u ( ixx, iyy, izz, i_comp, 8)
     1         - u ( ixx, iyy, izz, i_comp, 9)
     1         - u ( ixx, iyy, izz, i_comp, 10)
     1         + u ( ixx, iyy, izz, i_comp, 14)
     1         - u ( ixx, iyy, izz, i_comp, 15)
     1         + u ( ixx, iyy, izz, i_comp, 16)
     1         - u ( ixx, iyy, izz, i_comp, 17)
     1         + u ( ixx, iyy, izz, i_comp, 18)
     1         + u ( ixx, iyy, izz, i_comp, 19)

           Is4 = R4_xyz**2 + b

! -----

            R5_xyz = u ( ixxm1, iyym1, izzp1, i_comp, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12)
     1                     -u ( ixx, iyy, izz, i_comp, 13))
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 6)
     1         - u ( ixx, iyy, izz, i_comp, 7)
     1         - u ( ixx, iyy, izz, i_comp, 8)
     1         + u ( ixx, iyy, izz, i_comp, 9)
     1         + u ( ixx, iyy, izz, i_comp, 10)
     1         + u ( ixx, iyy, izz, i_comp, 14)
     1         - u ( ixx, iyy, izz, i_comp, 15)
     1         + u ( ixx, iyy, izz, i_comp, 16)
     1         - u ( ixx, iyy, izz, i_comp, 17)
     1         + u ( ixx, iyy, izz, i_comp, 18)
     1         + u ( ixx, iyy, izz, i_comp, 19)

           Is5 = R5_xyz**2 + b

! -----

            R6_xyz = u ( ixxm1, iyyp1, izzm1, i_comp, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     -u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         + u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 6)
     1         - u ( ixx, iyy, izz, i_comp, 7)
     1         + u ( ixx, iyy, izz, i_comp, 8)
     1         + u ( ixx, iyy, izz, i_comp, 9)
     1         - u ( ixx, iyy, izz, i_comp, 10)
     1         - u ( ixx, iyy, izz, i_comp, 14)
     1         + u ( ixx, iyy, izz, i_comp, 15)
     1         + u ( ixx, iyy, izz, i_comp, 16)
     1         + u ( ixx, iyy, izz, i_comp, 17)
     1         + u ( ixx, iyy, izz, i_comp, 18)
     1         - u ( ixx, iyy, izz, i_comp, 19)

           Is6 = R6_xyz**2 + b

! -----

            R7_xyz = u ( ixxp1, iyym1, izzm1, i_comp, 1)
     1         + 1.1 * (-u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1         - u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         + u ( ixx, iyy, izz, i_comp, 3)
     1         + u ( ixx, iyy, izz, i_comp, 4)
     1         - u ( ixx, iyy, izz, i_comp, 5)
     1         - u ( ixx, iyy, izz, i_comp, 6)
     1         - u ( ixx, iyy, izz, i_comp, 7)
     1         + u ( ixx, iyy, izz, i_comp, 8)
     1         - u ( ixx, iyy, izz, i_comp, 9)
     1         + u ( ixx, iyy, izz, i_comp, 10)
     1         + u ( ixx, iyy, izz, i_comp, 14)
     1         + u ( ixx, iyy, izz, i_comp, 15)
     1         - u ( ixx, iyy, izz, i_comp, 16)
     1         + u ( ixx, iyy, izz, i_comp, 17)
     1         - u ( ixx, iyy, izz, i_comp, 18)
     1         + u ( ixx, iyy, izz, i_comp, 19)

           Is7 = R7_xyz**2 + b

! -----

            R8_xyz = - u ( ixxm1, iyym1, izzm1, i_comp, 1)
     1         - 1.1 * ( u ( ixx, iyy, izz, i_comp, 11)
     1                     +u ( ixx, iyy, izz, i_comp, 12)
     1                     +u ( ixx, iyy, izz, i_comp, 13))
     1         + u ( ixx, iyy, izz, i_comp, 1)
     1         - u ( ixx, iyy, izz, i_comp, 2)
     1         - u ( ixx, iyy, izz, i_comp, 3)
     1         - u ( ixx, iyy, izz, i_comp, 4)
     1         + u ( ixx, iyy, izz, i_comp, 5)
     1         + u ( ixx, iyy, izz, i_comp, 6)
     1         + u ( ixx, iyy, izz, i_comp, 7)
     1         + u ( ixx, iyy, izz, i_comp, 8)
     1         + u ( ixx, iyy, izz, i_comp, 9)
     1         + u ( ixx, iyy, izz, i_comp, 10)
     1         - u ( ixx, iyy, izz, i_comp, 14)
     1         - u ( ixx, iyy, izz, i_comp, 15)
     1         - u ( ixx, iyy, izz, i_comp, 16)
     1         - u ( ixx, iyy, izz, i_comp, 17)
     1         - u ( ixx, iyy, izz, i_comp, 18)
     1         - u ( ixx, iyy, izz, i_comp, 19)

           Is8 = R8_xyz**2 + b

! -----

            Wt1 = 1.0 / Is1**4
            Wt2 = 1.0 / Is2**4
            Wt3 = 1.0 / Is3**4
            Wt4 = 1.0 / Is4**4
            Wt5 = 1.0 / Is5**4
            Wt6 = 1.0 / Is6**4
            Wt7 = 1.0 / Is7**4
            Wt8 = 1.0 / Is8**4

            c = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4 + Wt5 + Wt6 + Wt7 + Wt8 )

            Wt1 = Wt1 * c
            Wt2 = Wt2 * c
            Wt3 = Wt3 * c
            Wt4 = Wt4 * c
            Wt5 = Wt5 * c
            Wt6 = Wt6 * c
            Wt7 = Wt7 * c
            Wt8 = Wt8 * c

            u ( ixx, iyy, izz, i_comp, 20) =
     1        ( Wt1 * R1_xyz + Wt2 * R2_xyz + Wt3 * R3_xyz
     1        + Wt4 * R4_xyz + Wt5 * R5_xyz + Wt6 * R6_xyz
     1        + Wt7 * R7_xyz + Wt8 * R8_xyz)
     1                            * ( 1.0 - flattengr ( ixx, iyy, izz))

! -----

            END DO
! End of "DO ixx ="

          END DO
! End of "DO iyy ="

        END DO
! End of "DO ixx ="

!$OMP END PARALLEL DO

! ----------
      ELSE
! ----------

        u ( :, :, :, i_comp, 20) = 0.0

! ----------
      END IF
! End of "( ioffx .GT. 1) .AND. ( ioffy .GT. 1) .AND. ( ioffz .GT. 1)"

! ------------------------------

      END DO
! End of "DO i_comp ="

! ----------------------------------------------------------------------
#endif
! End of "NUMORDER == 4"
! End of limiting the cross-terms at 4th order.
! ----------------------------------------------------------------------
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

          DO ixx = 1 - ioffx, ix1 + ioffx
          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)

! ----------

          c = 1.0 - flattengr ( ixx, iyy, izz)

! ----------

          a = entropygr ( ixxp1, iyy, izz, 1)
     1      - entropygr ( ixx, iyy, izz, 1)

          b = entropygr ( ixx, iyy, izz, 1)
     1      - entropygr ( ixxm1, iyy, izz, 1)

          entropygr ( ixx, iyy, izz, 2) = MC_beta ( a, b, mc_coef) * c

! ----------

          a = entropygr ( ixx, iyyp1, izz, 1)
     1      - entropygr ( ixx, iyy, izz, 1)

          b = entropygr ( ixx, iyy, izz, 1)
     1      - entropygr ( ixx, iyym1, izz, 1)

          entropygr ( ixx, iyy, izz, 3) = MC_beta ( a, b, mc_coef) * c

! ----------

          a = entropygr ( ixx, iyy, izzp1, 1)
     1      - entropygr ( ixx, iyy, izz, 1)

          b = entropygr ( ixx, iyy, izz, 1)
     1      - entropygr ( ixx, iyy, izzm1, 1)

          entropygr ( ixx, iyy, izz, 4) = MC_beta ( a, b, mc_coef) * c

! ----------

          a =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          b =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          c =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          d =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 4)

          e =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          f =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          g =  entropygr ( ixx, iyy, izz, 1)
     1 + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          h =  entropygr ( ixx, iyy, izz, 1)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 2)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 3)
     1 - 0.5 * entropygr ( ixx, iyy, izz, 4)

          var_zone_min = AMIN1 ( a, b, c, d, e, f, g, h)

          var_zone_max = AMAX1 ( a, b, c, d, e, f, g, h)

! ----------

          var_neighbor_min = AMIN1 ( entropygr ( ixx, iyy, izz, 1),
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

          var_neighbor_max = AMAX1 ( entropygr ( ixx, iyy, izz, 1),
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

          b = ( var_neighbor_max - entropygr ( ixx, iyy, izz, 1) )
     1      / ( var_zone_max - entropygr ( ixx, iyy, izz, 1)
     1        + smallnum)

          c = ( entropygr ( ixx, iyy, izz, 1) - var_neighbor_min )
     1      / ( entropygr ( ixx, iyy, izz, 1) - var_zone_min
     1        + smallnum)

          a = AMIN1 ( 1.0, psi_limiter * AMIN1 ( b, c) )

          entropygr ( ixx, iyy, izz, 2) =
     1    a * entropygr ( ixx, iyy, izz, 2)
          entropygr ( ixx, iyy, izz, 3) =
     1    a * entropygr ( ixx, iyy, izz, 3)
          entropygr ( ixx, iyy, izz, 4) =
     1    a * entropygr ( ixx, iyy, izz, 4)

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
#endif
! End of "PRESSURE_POSITIVITY == 1".
! ----------------------------------------------------------------------
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
!$OMP&    bfxbdy, bfybdy, bfzbdy, msonicgr)                             &
!$OMP&           REDUCTION ( MIN: div_min) REDUCTION ( MAX: div_max)    &
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

          dx = dx_save ( ixx)
          dy = dy_save ( iyy)
          dz = dz_save ( izz)

          ob_dx = 1.0 / dx
          ob_dy = 1.0 / dy
          ob_dz = 1.0 / dz

          END IF

          IF ( igeom .EQ. 2) THEN

          dx = dx_save ( ixx)
          dy = xc_save ( ixx) * dy_save ( iyy)
          dz = dz_save ( izz)

          ob_dx = 1.0 / dx
          ob_dy = 1.0 / dy
          ob_dz = 1.0 / dz

          END IF

          IF ( igeom .EQ. 3) THEN

          dx = dx_save ( ixx)
          dy = xc_save ( ixx) * dy_save ( iyy)
          dz = xc_save ( ixx) * SIN ( yc_save ( iyy)) * dz_save ( izz)

          ob_dx = 1.0 / dx
          ob_dy = 1.0 / dy
          ob_dz = 1.0 / dz

          END IF

! ----------

          a = ( bfxbdy ( ixx, iyy, izz, 1)
     1        - bfxbdy ( ixx - 1, iyy, izz, 1) ) * ob_dx
     1      + ( bfybdy ( ixx, iyy, izz, 1)
     1        - bfybdy ( ixx, iyy - 1, izz, 1) ) * ob_dy
     1      + ( bfzbdy ( ixx, iyy, izz, 1)
     1        - bfzbdy ( ixx, iyy, izz - 1, 1) ) * ob_dz

          a = a * AMIN1 ( dx, dy, dz )

          b = u ( ixx, iyy, izz, 6, 1)**2 + u ( ixx, iyy, izz, 7, 1)**2
     1      + u ( ixx, iyy, izz, 8, 1)**2
          b = SQRT ( b)
          b = AMAX1 ( b, smallnum)

          a = a / b
!          a = a / AMAX1 ( b, msonicgr ( ixx, iyy, izz)
!     1                     * SQRT ( rhogr ( ixx, iyy, izz) ) )

          div_min = AMIN1 ( div_min, a)
          div_max = AMAX1 ( div_max, a)

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

! ----------

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      WRITE ( 6, *)"From LIMIT_ALL_MOMENTS: div_min, div_max = ",
     1              div_min, div_max
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

      END SUBROUTINE LIMIT_ALL_MOMENTS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LIMIT_1D_MHD_CONS
     1   ( oned_begin, oned_end, oned_off, limit_characteristics,

     1     gamma, smallnum, prsfloor, rhofloor, pi, gasconst,
     1     cp_specie, molewt_specie, heat_of_formation,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, bfxbdy_1d,
     1     bfxbar, bfybar, bfzbar, speciebar, flatten_1d, u_1d)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Does 1d conserved and/or characteristic limiting for the MHD variables.
! Characteristic variables are limited or not depending on
! "limit_characteristics = 0/1".
! The variables are assumed to come in with variable and slope
! already supplied in "u_1d". The primitive variables are also
! supplied and are used for building eigenvectors. The limited moments are
! returned in "u_1d ( :, :, NUMORDER)".

! We also check whether the boundary extrapolated states produce positive
! densities and pressures. If not, we zero out the slope in that zone.

! SLOPE_LIMITER == 2 yields WENO limiting at 2nd, 3rd and 4th order.
! SLOPE_LIMITER == 1 yields MC limiting at 2nd order and PPM at 3rd order. At
!                    4th order it yields a hybrid r=5 WENO + Van Albada limiter.

! LIMIT_TWICE == 1 triggers characteristic limiting along with limiting
! of the conserved variables. Otherwise, we just limit the conserved variables.

! ----------------------------------------------------------------------

      IMPLICIT NONE

! These variables are passed across the subroutine boundary.

      INTEGER oned_begin, oned_end, oned_off, limit_characteristics

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
     1                  NFIELD + NFLUID, NUMORDER) :: u_1d

! ----------------------------------------------------------------------

      INTEGER nriembegin, nriemend, iprimitive, nsys, ix, j, k, ifield,
     1        ix_shrink, debug_stop

      REAL ob_8pi, ob_3, ob_6, ob_12, ob_20, ob_24, ob_56, ob_60,
     1     ob_120, r13_3, r781_20, r123_455, r1421461_2275,
     1     a, b, c, d, e, tau, MinMod, MC_beta,
     1     u_ip3, u_ip2, u_ip1, u_i, u_im1, u_im2, u_im3,
     1     R1_x, R2_x, R3_x, R4_x, R1_x2, R2_x2, R3_x2, R4_x2,
     1     R1_x3, R2_x3, R3_x3, R4_x3, R1_x4,
     1     Is1, Is2, Is3, Is4, Wt1, Wt2, Wt3, Wt4,
     1     dd_u_0, dd_u_l, dd_u_r, dd_u_c, dd_u_lim,
     1     alpha_rt, alpha_lf, u_6,
     1     rh, en, vx, vy, vz, bx, by, bz

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: xbar,
     1     u_rt, u_lf, du_lim

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                                NFIELD + NFLUID) :: eigenval

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD + NFLUID, NFIELD + NFLUID) ::
     1                                    lf_eigenvec, rt_eigenvec

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID, NUMORDER) :: del_eigwt, del_u

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + NFLUID, - NUMORDER: NUMORDER) :: eigwt

! MinMod limiter
      MinMod ( a, b) = AMIN1 ( ABS ( a), ABS ( b))
     1                    * ( SIGN ( 0.5, a) + SIGN ( 0.5, b))

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( a, b, mc_coef) = AMIN1 ( 0.5 * ( ABS ( a + b)),
     1                         mc_coef * ABS ( a), mc_coef * ABS ( b) )
     1                    * ( SIGN ( 0.5, a) + SIGN ( 0.5, b))

! ----------------------------------------------------------------------

      nriembegin = oned_begin - oned_off
      nriemend = oned_end + oned_off
      iprimitive = 0

      nsys = NFIELD + NFLUID

      ob_8pi = 1.0 / ( 8.0 * pi)

      ob_3 = 1.0 / 3.0

      ob_6 = 1.0 / 6.0

      ob_12 = 1.0 / 12.0

      ob_20 = 1.0 / 20.0

      ob_24 = 1.0 / 24.0

      ob_56 = 1.0 / 56.0

      ob_60 = 1.0 / 60.0

      ob_120 = 1.0 / 120.0

      r13_3 = 13.0 / 3.0

      r123_455 = 123.0 / 455.0

      r781_20 = 781.0 / 20.0

      r1421461_2275 = 1421461.0 / 2275.0

      debug_stop = 0

! Set the range by which the zones need to be shrunk. 2 for 2nd and 3rd order;
! 3 for 4th order. PPM, however, requires 3 in characteristic space.

      ix_shrink = 3
      IF ( NUMORDER <= 3) ix_shrink = 2

! It is good to leave a 2 zone margin on either side where reconstruction
! is allowed to act.

      IF ( ( oned_begin - oned_off + ix_shrink > oned_begin - 2) .OR.
     1     ( oned_end + oned_off - ix_shrink < oned_end + 2) ) THEN

      WRITE ( 6, *)
     1     "Stopping in LIMIT_1D_MHD_CONS: Do not have 2 zone margin"
      WRITE ( 6, *)"left margins: ", oned_begin - oned_off + ix_shrink,
     1                               oned_begin - 2
      WRITE ( 6, *)"right margins: ", oned_end + oned_off - ix_shrink,
     1                                oned_end + 2

      WRITE ( 6, *)"order, limiter, characteristics = ", NUMORDER,
     1              SLOPE_LIMITER, limit_characteristics

      STOP
      END IF

! ------------------------------

      IF ( limit_characteristics .EQ. 1) THEN

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

      END IF
! End of "limit_characteristics .EQ. 1".

! ----------------------------------------------------------------------

! Now project into the space of eigenweights.

      IF ( limit_characteristics .EQ. 1) THEN

      eigwt ( oned_begin - oned_off: oned_end + oned_off, :, :) = 0.0

      DO j = - ix_shrink, ix_shrink

      DO k = 1, nsys

        DO ifield = 1, nsys

          DO ix = oned_begin - oned_off + ix_shrink,
     1            oned_end + oned_off - ix_shrink

          eigwt ( ix, ifield, j) = eigwt ( ix, ifield, j)
     1  + lf_eigenvec ( ix, ifield, k) * u_1d ( ix + j, k, 1)

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

      END DO ! End of "DO k ="

      END DO ! End of "DO j ="

      END IF
! End of "limit_characteristics .EQ. 1".

! ----------------------------------------------------------------------
#if ( NUMORDER == 2)
#if ( SLOPE_LIMITER == 1)
! ------------------------------------------------------------

! Second order TVD limiter being done here.

! MC-beta limiter in physical and characteristic space.

      DO ifield = 1, nsys

! --------------------

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_1d ( ix, ifield, 2) =
     1    MC_beta ( u_1d ( ix + 1, ifield, 1) - u_1d ( ix, ifield, 1),
     1              u_1d ( ix, ifield, 1) - u_1d ( ix - 1, ifield, 1),
     1              mc_coef )

        END DO ! End of "DO ix ="

! ----------------------------------------

        IF ( limit_characteristics .EQ. 1) THEN

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          del_eigwt ( ix, ifield, 2) =
     1    MC_beta ( eigwt ( ix, ifield, 1) - eigwt ( ix, ifield, 0),
     1              eigwt ( ix, ifield, 0) - eigwt ( ix, ifield, -1),
     1              mc_coef )

        END DO ! End of "DO ix ="

        END IF
! End of "limit_characteristics .EQ. 1".

! --------------------

      END DO ! End of "DO ifield ="

! ------------------------------------------------------------
#else
! ------------------------------------------------------------

! First moment of r=3 WENO limiter being done here.

! r=3 WENO limiter in physical and characteristic space.

      DO ifield = 1, nsys

! --------------------

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip2 = u_1d ( ix + 2, ifield, 1)
          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i   = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)
          u_im2 = u_1d ( ix - 2, ifield, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_x2 = - u_im1 + 0.5 * u_im2 + 0.5 * u_i
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_x2 = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_x2 = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          u_1d ( ix, ifield, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                          + Wt3 * R3_x

        END DO ! End of "DO ix ="

! ----------------------------------------

        IF ( limit_characteristics .EQ. 1) THEN

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip2 = eigwt ( ix, ifield, 2)
          u_ip1 = eigwt ( ix, ifield, 1)
          u_i   = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)
          u_im2 = eigwt ( ix, ifield, -2)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_x2 = - u_im1 + 0.5 * u_im2 + 0.5 * u_i
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_x2 = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_x2 = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          del_eigwt ( ix, ifield, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                               + Wt3 * R3_x

        END DO ! End of "DO ix ="

        END IF
! End of "limit_characteristics .EQ. 1".

! --------------------

      END DO ! End of "DO ifield ="

! ------------------------------------------------------------
#endif
! End of "SLOPE_LIMITER == 1/2"
#endif
! End of "NUMORDER == 2".
! ----------------------------------------------------------------------
#if ( NUMORDER == 3)
#if ( SLOPE_LIMITER == 1)
! ------------------------------------------------------------

! Third order PPM limiter being done here.

! PPM limiter in physical and characteristic space.

      DO ifield = 1, nsys

! --------------------

! This is traditional TVD-preserving PPM in physical space.

! -----
! Build eq. (1.8) of CW, which is just an MC limiter for uniform zones.

        DO ix = oned_begin - oned_off + 1,
     1          oned_end + oned_off - 1

          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i   = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)

          a = u_ip1 - u_i
          b = u_i - u_im1

          du_lim ( ix) = MC_beta ( a, b, 2.0)

        END DO ! End of "DO ix ="

! -----
! Now build eq. (1.6) of CW.

        DO ix = oned_begin - oned_off + ix_shrink - 1,
     1          oned_end + oned_off - ix_shrink

! Eq. (1.9) of CW. We are building values at the cell edges.

          u_ip2 = u_1d ( ix + 2, ifield, 1)
          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i   = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)

          u_rt ( ix) = u_i + 0.5 * ( u_ip1 - u_i )
     1               - ob_6 * ( du_lim ( ix + 1) - du_lim ( ix) )

          u_lf ( ix + 1) = u_rt ( ix)

        END DO ! End of "DO ix ="

! -----
! Now limit the parabolae, eq. (1.10) of CW. PPM flattens profiles entirely when
! "alpha_rt * alpha_lf > 0.0"; which causes excessive flattening. My modification
! is to retain the linear profiles even when that condition tests true.

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)

! Note that this may not be the full MC limiter because of the use of "mc_coef".
          du_lim ( ix) = MC_beta ( u_ip1 - u_i, u_i - u_im1, mc_coef)

          alpha_rt = u_rt ( ix) - u_i
          alpha_lf = u_lf ( ix) - u_i

          IF ( alpha_rt * alpha_lf > 0.0 ) THEN

! The next two lines are the original PPM choice (see first case in the bottom 
! paragraph of pg. 178 of CW), and they are commented out.
!            u_rt ( ix) = u_i
!            u_lf ( ix) = u_i

            u_rt ( ix) = u_i + 0.5 * du_lim ( ix)
            u_lf ( ix) = u_i - 0.5 * du_lim ( ix)

          END IF

          a = u_rt ( ix) - u_lf ( ix)
          b = u_i - 0.5 * ( u_rt ( ix) + u_lf ( ix) )

          IF ( a * b > ob_6 * a**2 ) THEN
            u_lf ( ix) = 3.0 * u_i - 2.0 * u_rt ( ix)
          ELSE IF ( - ob_6 * a**2 > a * b) THEN
            u_rt ( ix) = 3.0 * u_i - 2.0 * u_lf ( ix)
          END IF

          u_1d ( ix, ifield, 2) = u_rt ( ix) - u_lf ( ix)
          u_1d ( ix, ifield, 3) = - ( 6.0 * u_i
     1                            - 3.0 * ( u_lf ( ix) + u_rt ( ix)) )

        END DO ! End of "DO ix ="

! -----

! End of traditional PPM in physical space.

! --------------------------------------------------

        IF ( limit_characteristics .EQ. 1) THEN

! --------------------

! This is traditional TVD-preserving PPM in characteristic space.

! -----
! Build eq. (1.8) of CW, which is just an MC limiter for uniform zones.

        DO ix = oned_begin - oned_off + 1,
     1          oned_end + oned_off - 1

          u_ip1 = eigwt ( ix, ifield, 1)
          u_i   = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)

          a = u_ip1 - u_i
          b = u_i - u_im1

          du_lim ( ix) = MC_beta ( a, b, 2.0)

        END DO ! End of "DO ix ="

! -----
! Now build eq. (1.6) of CW.

        DO ix = oned_begin - oned_off + ix_shrink - 1,
     1          oned_end + oned_off - ix_shrink

! Eq. (1.9) of CW. We are building values at the cell edges.

          u_ip2 = eigwt ( ix, ifield, 2)
          u_ip1 = eigwt ( ix, ifield, 1)
          u_i   = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)

          u_rt ( ix) = u_i + 0.5 * ( u_ip1 - u_i )
     1               - ob_6 * ( du_lim ( ix + 1) - du_lim ( ix) )

          u_lf ( ix + 1) = u_rt ( ix)

        END DO ! End of "DO ix ="

! -----
! Now limit the parabolae, eq. (1.10) of CW. PPM flattens profiles entirely when
! "alpha_rt * alpha_lf > 0.0"; which causes excessive flattening. My modification
! is to retain the linear profiles even when that condition tests true.

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip1 = eigwt ( ix, ifield, 1)
          u_i = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)

! Note that this may not be the full MC limiter because of the use of "mc_coef".
          du_lim ( ix) = MC_beta ( u_ip1 - u_i, u_i - u_im1, mc_coef)

          alpha_rt = u_rt ( ix) - u_i
          alpha_lf = u_lf ( ix) - u_i

          IF ( alpha_rt * alpha_lf > 0.0 ) THEN

! The next two lines are the original PPM choice (see first case in the bottom
! paragraph of pg. 178 of CW), and they are commented out.
!            u_rt ( ix) = u_i
!            u_lf ( ix) = u_i

            u_rt ( ix) = u_i + 0.5 * du_lim ( ix)
            u_lf ( ix) = u_i - 0.5 * du_lim ( ix)

          END IF

          a = u_rt ( ix) - u_lf ( ix)
          b = u_i - 0.5 * ( u_rt ( ix) + u_lf ( ix) )

          IF ( a * b > ob_6 * a**2 ) THEN
            u_lf ( ix) = 3.0 * u_i - 2.0 * u_rt ( ix)
          ELSE IF ( - ob_6 * a**2 > a * b) THEN
            u_rt ( ix) = 3.0 * u_i - 2.0 * u_lf ( ix)
          END IF

          del_eigwt ( ix, ifield, 2) = u_rt ( ix) - u_lf ( ix)
          del_eigwt ( ix, ifield, 3) = - ( 6.0 * u_i
     1                            - 3.0 * ( u_lf ( ix) + u_rt ( ix)) )

        END DO ! End of "DO ix ="

! -----

! End of traditional PPM in characteristic space.

! --------------------

        END IF
! End of "limit_characteristics .EQ. 1".

! --------------------

      END DO ! End of "DO ifield ="

! ------------------------------------------------------------
#else
! ------------------------------------------------------------

! Third order r=3 WENO being done here.

! r=3 WENO limiter in physical and characteristic space.

      DO ifield = 1, nsys

! ----------

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip2 = u_1d ( ix + 2, ifield, 1)
          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i   = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)
          u_im2 = u_1d ( ix - 2, ifield, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_x2 = - u_im1 + 0.5 * u_im2 + 0.5 * u_i
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_x2 = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_x2 = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          u_1d ( ix, ifield, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                          + Wt3 * R3_x

          u_1d ( ix, ifield, 3) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                          + Wt3 * R3_x2

        END DO ! End of "DO ix ="

! ----------------------------------------

        IF ( limit_characteristics .EQ. 1) THEN

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip2 = eigwt ( ix, ifield, 2)
          u_ip1 = eigwt ( ix, ifield, 1)
          u_i   = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)
          u_im2 = eigwt ( ix, ifield, -2)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_x2 = - u_im1 + 0.5 * u_im2 + 0.5 * u_i
          Is1 = R1_x**2 + r13_3 * R1_x2**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_x2 = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_x2**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_x2 = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_x2**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          del_eigwt ( ix, ifield, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                               + Wt3 * R3_x

          del_eigwt ( ix, ifield, 3) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                               + Wt3 * R3_x2

        END DO ! End of "DO ix ="

        END IF
! End of "limit_characteristics .EQ. 1".

! --------------------

      END DO ! End of "DO ifield ="

! ------------------------------------------------------------
#endif
! End of "SLOPE_LIMITER == 1/2"
#endif
! End of "NUMORDER == 3".
! ----------------------------------------------------------------------
#if ( NUMORDER == 4)
#if ( SLOPE_LIMITER == 1)
! ------------------------------------------------------------

! r=5 WENO limiter, which is hybridized with a Van Albada limiter being done here.

! Hybrid r=5 WENO + Van Albada limiter in physical and characteristic space.

      DO ifield = 1, nsys

! --------------------

! This is hybrid r=5 WENO + Van Albada limiter in physical space.

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip2 = u_1d ( ix + 2, ifield, 1)
          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i   = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)
          u_im2 = u_1d ( ix - 2, ifield, 1)

          R1_x = ( - 82.0 * u_im1 + 11.0 * u_im2 + 82.0 * u_ip1
     1             - 11.0 * u_ip2 ) * ob_120
          R1_x2 = ( 40.0 * u_im1 - 3.0 * u_im2 - 74.0 * u_i
     1            + 40.0 * u_ip1 - 3.0 * u_ip2 ) * ob_56
          R1_x3 = ( 2.0 * u_im1 - u_im2 - 2.0 * u_ip1 + u_ip2 ) * ob_12
          R1_x4 = ( - 4.0 * u_im1 + u_im2 + 6.0 * u_i
     1              - 4.0 * u_ip1 + u_ip2 ) * ob_24

          Is1 = ( R1_x + 0.1 * R1_x3)**2
     1        + r13_3 * ( R1_x2 + r123_455 * R1_x4)**2
     1        + r781_20 * R1_x3**2 + r1421461_2275 * R1_x4**2

          R2_x = u_i - u_im1
          Is2 = R2_x**2

          R3_x = u_ip1 - u_i
          Is3 = R3_x**2

          tau = ( 0.5 * ( ABS ( Is1 - Is2) + ABS ( Is1 - Is3) ) )**2

          Wt1 = 0.98 * ( 1.0 + tau / ( Is1 + smallnum) )
          Wt2 = 0.01 * ( 1.0 + tau / ( Is2 + smallnum) )
          Wt3 = 0.01 * ( 1.0 + tau / ( Is3 + smallnum) )

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1 / 0.98
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          u_1d ( ix, ifield, 2) =
     1          Wt1 * ( R1_x - 0.01 * R2_x - 0.01 * R3_x)
     1        + Wt2 * R2_x + Wt3 * R3_x

          u_1d ( ix, ifield, 3) = Wt1 * R1_x2

          u_1d ( ix, ifield, 4) = Wt1 * R1_x3

       END DO ! End of "DO ix ="

! ----------------------------------------

! This is hybrid r=5 WENO + Van Albada limiter in characteristic space.

        IF ( limit_characteristics .EQ. 1) THEN

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip2 = eigwt ( ix, ifield, 2)
          u_ip1 = eigwt ( ix, ifield, 1)
          u_i   = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)
          u_im2 = eigwt ( ix, ifield, -2)

          R1_x = ( - 82.0 * u_im1 + 11.0 * u_im2 + 82.0 * u_ip1
     1             - 11.0 * u_ip2 ) * ob_120
          R1_x2 = ( 40.0 * u_im1 - 3.0 * u_im2 - 74.0 * u_i
     1            + 40.0 * u_ip1 - 3.0 * u_ip2 ) * ob_56
          R1_x3 = ( 2.0 * u_im1 - u_im2 - 2.0 * u_ip1 + u_ip2 ) * ob_12
          R1_x4 = ( - 4.0 * u_im1 + u_im2 + 6.0 * u_i
     1              - 4.0 * u_ip1 + u_ip2 ) * ob_24

          Is1 = ( R1_x + 0.1 * R1_x3)**2
     1        + r13_3 * ( R1_x2 + r123_455 * R1_x4)**2
     1        + r781_20 * R1_x3**2 + r1421461_2275 * R1_x4**2

          R2_x = u_i - u_im1
          Is2 = R2_x**2

          R3_x = u_ip1 - u_i
          Is3 = R3_x**2

          tau = ( 0.5 * ( ABS ( Is1 - Is2) + ABS ( Is1 - Is3) ) )**2

          Wt1 = 0.98 * ( 1.0 + tau / ( Is1 + smallnum) )
          Wt2 = 0.01 * ( 1.0 + tau / ( Is2 + smallnum) )
          Wt3 = 0.01 * ( 1.0 + tau / ( Is3 + smallnum) )

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1 / 0.98
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          del_eigwt ( ix, ifield, 2) = 
     1          Wt1 * ( R1_x - 0.01 * R2_x - 0.01 * R3_x) 
     1        + Wt2 * R2_x + Wt3 * R3_x

          del_eigwt ( ix, ifield, 3) = Wt1 * R1_x2

          del_eigwt ( ix, ifield, 4) = Wt1 * R1_x3

        END DO ! End of "DO ix ="

        END IF
! End of "limit_characteristics .EQ. 1".

! --------------------

      END DO ! End of "DO ifield ="

! ------------------------------------------------------------
#else
! ------------------------------------------------------------

! r=4 WENO being done here.

! r=4 WENO limiter in physical and characteristic space.

      DO ifield = 1, nsys

! --------------------

! This is r=4 WENO in physical space.

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip3 = u_1d ( ix + 3, ifield, 1)
          u_ip2 = u_1d ( ix + 2, ifield, 1)
          u_ip1 = u_1d ( ix + 1, ifield, 1)
          u_i   = u_1d ( ix, ifield, 1)
          u_im1 = u_1d ( ix - 1, ifield, 1)
          u_im2 = u_1d ( ix - 2, ifield, 1)
          u_im3 = u_1d ( ix - 3, ifield, 1)

          R1_x = ( - 177.0 * u_im1 + 87.0 * u_im2
     1           - 19.0 * u_im3 + 109.0 * u_i ) * ob_60
          R1_x2 = - 2.5 * u_im1 + 2.0 * u_im2 - 0.5 * u_im3 + u_i
          R1_x3 = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i ) * ob_6
          Is1 = ( R1_x + 0.1 * R1_x3)**2 + 13.0 * R1_x2**2 * ob_3
     1          + 39.05 * R1_x3**2

          R2_x = ( - 63.0 * u_im1 + 11.0 * u_im2
     1           + 33.0 * u_i + 19.0 * u_ip1 ) * ob_60
          R2_x2 = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R2_x3 = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1 ) * ob_6
          Is2 = ( R2_x + 0.1 * R2_x3)**2 + 13.0 * R2_x2**2 * ob_3
     1          + 39.05 * R2_x3**2

          R3_x = ( - 19.0 * u_im1 - 33.0 * u_i
     1           + 63.0 * u_ip1 - 11.0 * u_ip2 ) * ob_60
          R3_x2 = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R3_x3 = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2 ) * ob_6
          Is3 = ( R3_x + 0.1 * R3_x3)**2 + 13.0 * R3_x2**2 * ob_3
     1          + 39.05 * R3_x3**2

          R4_x = ( - 109.0 * u_i + 177.0 * u_ip1
     1           - 87.0 * u_ip2 + 19.0 * u_ip3 ) * ob_60
          R4_x2 = u_i - 2.5 * u_ip1 + 2.0 * u_ip2 - 0.5 * u_ip3
          R4_x3 = ( - u_i + 3.0 * u_ip1 - 3.0 * u_ip2 + u_ip3 ) * ob_6
          Is4 = ( R4_x + 0.1 * R4_x3)**2 + 13.0 * R4_x2**2 * ob_3
     1          + 39.05 * R4_x3**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = central_zone_wt / ( Is3 + smallnum)**4
          Wt4 = 1.0 / ( Is4 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3
          Wt4 = a * Wt4

          u_1d ( ix, ifield, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                          + Wt3 * R3_x + Wt4 * R4_x

          u_1d ( ix, ifield, 3) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                          + Wt3 * R3_x2 + Wt4 * R4_x2

          u_1d ( ix, ifield, 4) = Wt1 * R1_x3 + Wt2 * R2_x3
     1                          + Wt3 * R3_x3 + Wt4 * R4_x3

        END DO ! End of "DO ix ="

! ----------------------------------------

! This is r=4 WENO in characteristic space.

        IF ( limit_characteristics .EQ. 1) THEN

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

          u_ip3 = eigwt ( ix, ifield, 3)
          u_ip2 = eigwt ( ix, ifield, 2)
          u_ip1 = eigwt ( ix, ifield, 1)
          u_i   = eigwt ( ix, ifield, 0)
          u_im1 = eigwt ( ix, ifield, -1)
          u_im2 = eigwt ( ix, ifield, -2)
          u_im3 = eigwt ( ix, ifield, -3)

          R1_x = ( - 177.0 * u_im1 + 87.0 * u_im2
     1           - 19.0 * u_im3 + 109.0 * u_i ) * ob_60
          R1_x2 = - 2.5 * u_im1 + 2.0 * u_im2 - 0.5 * u_im3 + u_i
          R1_x3 = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i ) * ob_6
          Is1 = ( R1_x + 0.1 * R1_x3)**2 + 13.0 * R1_x2**2 * ob_3
     1          + 39.05 * R1_x3**2

          R2_x = ( - 63.0 * u_im1 + 11.0 * u_im2
     1           + 33.0 * u_i + 19.0 * u_ip1 ) * ob_60
          R2_x2 = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R2_x3 = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1 ) * ob_6
          Is2 = ( R2_x + 0.1 * R2_x3)**2 + 13.0 * R2_x2**2 * ob_3
     1          + 39.05 * R2_x3**2

          R3_x = ( - 19.0 * u_im1 - 33.0 * u_i
     1           + 63.0 * u_ip1 - 11.0 * u_ip2 ) * ob_60
          R3_x2 = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R3_x3 = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2 ) * ob_6
          Is3 = ( R3_x + 0.1 * R3_x3)**2 + 13.0 * R3_x2**2 * ob_3
     1          + 39.05 * R3_x3**2

          R4_x = ( - 109.0 * u_i + 177.0 * u_ip1
     1           - 87.0 * u_ip2 + 19.0 * u_ip3 ) * ob_60
          R4_x2 = u_i - 2.5 * u_ip1 + 2.0 * u_ip2 - 0.5 * u_ip3
          R4_x3 = ( - u_i + 3.0 * u_ip1 - 3.0 * u_ip2 + u_ip3 ) * ob_6
          Is4 = ( R4_x + 0.1 * R4_x3)**2 + 13.0 * R4_x2**2 * ob_3
     1          + 39.05 * R4_x3**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = central_zone_wt / ( Is3 + smallnum)**4
          Wt4 = 1.0 / ( Is4 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3
          Wt4 = a * Wt4

          del_eigwt ( ix, ifield, 2) = Wt1 * R1_x + Wt2 * R2_x
     1                               + Wt3 * R3_x + Wt4 * R4_x

          del_eigwt ( ix, ifield, 3) = Wt1 * R1_x2 + Wt2 * R2_x2
     1                               + Wt3 * R3_x2 + Wt4 * R4_x2

          del_eigwt ( ix, ifield, 4) = Wt1 * R1_x3 + Wt2 * R2_x3
     1                               + Wt3 * R3_x3 + Wt4 * R4_x3

        END DO ! End of "DO ix ="

        END IF
! End of "limit_characteristics .EQ. 1".

! --------------------

      END DO ! End of "DO ifield ="

! ------------------------------------------------------------
#endif
! End of "SLOPE_LIMITER == 1/2"
#endif
! End of "NUMORDER == 4".
! ----------------------------------------------------------------------

! Now project back into physical space.

      IF ( limit_characteristics .EQ. 1) THEN

      del_u ( oned_begin - oned_off: oned_end + oned_off, :, :) = 0.0

      DO k = 1, nsys

        DO ifield = 1, nsys

          DO ix = oned_begin - oned_off + ix_shrink,
     1            oned_end + oned_off - ix_shrink

! ----------

          del_u ( ix, k, 2) = del_u ( ix, k, 2)
     1                      + del_eigwt ( ix, ifield, 2)
     1                      * rt_eigenvec ( ix, k, ifield)

! ----------

#if ( NUMORDER >= 3)
          del_u ( ix, k, 3) = del_u ( ix, k, 3)
     1                      + del_eigwt ( ix, ifield, 3)
     1                      * rt_eigenvec ( ix, k, ifield)
#endif

! ----------

#if ( NUMORDER >= 4)
          del_u ( ix, k, 4) = del_u ( ix, k, 4)
     1                      + del_eigwt ( ix, ifield, 4)
     1                      * rt_eigenvec ( ix, k, ifield)
#endif

! ----------

          END DO ! End of "DO ix ="

        END DO ! End of "DO ifield ="

      END DO ! End of "DO k ="

      END IF
! End of "limit_characteristics .EQ. 1".

! ----------------------------------------------------------------------

! And choose the minimum of the physical space and characteristic space
! limiters. Also introduce a flattener. One could also choose pure
! characteristic limiting.
      
      DO k = 1, nsys

! ----------

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

        a = 1.0 - flatten_1d ( ix)

        u_1d ( ix, k, 2) = u_1d ( ix, k, 2) * a
#if ( NUMORDER >= 3)
        u_1d ( ix, k, 3) = u_1d ( ix, k, 3) * a
#endif
#if ( NUMORDER >= 4)
        u_1d ( ix, k, 4) = u_1d ( ix, k, 4) * a
#endif

        END DO ! End of "DO ix ="

! ----------
! Here we limit twice.

        IF ( limit_characteristics .EQ. 1) THEN

        DO ix = oned_begin - oned_off + ix_shrink,
     1          oned_end + oned_off - ix_shrink

        u_1d ( ix, k, 2) = MinMod ( u_1d ( ix, k, 2),
     1                             del_u ( ix, k, 2) * a )
#if ( NUMORDER >= 3)
        u_1d ( ix, k, 3) = MinMod ( u_1d ( ix, k, 3),
     1                             del_u ( ix, k, 3) * a )
#endif
#if ( NUMORDER >= 4)
        u_1d ( ix, k, 4) = MinMod ( u_1d ( ix, k, 4),
     1                             del_u ( ix, k, 4) * a )
#endif

        END DO ! End of "DO ix ="

        END IF
! End of "limit_characteristics .EQ. 1".

! ----------

      END DO ! End of "DO k ="

! ------------------------------
! For the small number of zones on either side for which limiting can't
! be done, zero out the moments.

      u_1d( oned_begin - oned_off: oned_begin - oned_off + ix_shrink-1,
     1                                :, 2: NUMORDER) = 0.0

      u_1d ( oned_end + oned_off - ix_shrink + 1: oned_end + oned_off,
     1                                :, 2: NUMORDER) = 0.0

! ----------------------------------------------------------------------

! Check the boundary extrapolated states for positivity of density & pressure.
! If that is not achieved, zero the slopes in that zone.

! ------------------------------
#if ( NUMORDER == 2)
! ------------------------------

      DO ix = oned_begin - oned_off + ix_shrink,
     1        oned_end + oned_off - ix_shrink

        rh = u_1d ( ix, 1, 1) + 0.5 * u_1d ( ix, 1, 2)
        en = u_1d ( ix, 5, 1) + 0.5 * u_1d ( ix, 5, 2)

        vx = u_1d ( ix, 2, 1) + 0.5 * u_1d ( ix, 2, 2)
        vy = u_1d ( ix, 3, 1) + 0.5 * u_1d ( ix, 3, 2)
        vz = u_1d ( ix, 4, 1) + 0.5 * u_1d ( ix, 4, 2)

        bx = bfxbdy_1d ( ix)
        by = u_1d ( ix, 6, 1) + 0.5 * u_1d ( ix, 6, 2)
        bz = u_1d ( ix, 7, 1) + 0.5 * u_1d ( ix, 7, 2)

        IF ( ( rh <= rhofloor) .OR.
     1       ( en - 0.5 * ( vx**2 + vy**2 + vz**2) / rh
     1       - ob_8pi * ( bx**2 + by**2 + bz**2) <= prsfloor) ) THEN
        u_1d ( ix, :, 2: NUMORDER) = 0.0 
        END IF

      END DO

! ----------

      DO ix = oned_begin - oned_off + ix_shrink,
     1        oned_end + oned_off - ix_shrink

        rh = u_1d ( ix, 1, 1) - 0.5 * u_1d ( ix, 1, 2)
        en = u_1d ( ix, 5, 1) - 0.5 * u_1d ( ix, 5, 2)

        vx = u_1d ( ix, 2, 1) - 0.5 * u_1d ( ix, 2, 2)
        vy = u_1d ( ix, 3, 1) - 0.5 * u_1d ( ix, 3, 2)
        vz = u_1d ( ix, 4, 1) - 0.5 * u_1d ( ix, 4, 2)

        bx = bfxbdy_1d ( ix - 1)
        by = u_1d ( ix, 6, 1) - 0.5 * u_1d ( ix, 6, 2)
        bz = u_1d ( ix, 7, 1) - 0.5 * u_1d ( ix, 7, 2)

        IF ( ( rh <= rhofloor) .OR.
     1       ( en - 0.5 * ( vx**2 + vy**2 + vz**2) / rh
     1       - ob_8pi * ( bx**2 + by**2 + bz**2) <= prsfloor) ) THEN
        u_1d ( ix, :, 2: NUMORDER) = 0.0
        END IF

      END DO

! ------------------------------
#endif
! End of "NUMORDER == 2".
! ------------------------------
#if ( NUMORDER == 3)
! ------------------------------

      DO ix = oned_begin - oned_off + ix_shrink,
     1        oned_end + oned_off - ix_shrink

        rh = u_1d ( ix, 1, 1) + 0.5 * u_1d ( ix, 1, 2)
     1     + ob_6 * u_1d ( ix, 1, 3)
        en = u_1d ( ix, 5, 1) + 0.5 * u_1d ( ix, 5, 2)
     1     + ob_6 * u_1d ( ix, 5, 3)

        vx = u_1d ( ix, 2, 1) + 0.5 * u_1d ( ix, 2, 2)
     1     + ob_6 * u_1d ( ix, 2, 3)
        vy = u_1d ( ix, 3, 1) + 0.5 * u_1d ( ix, 3, 2)
     1     + ob_6 * u_1d ( ix, 3, 3)
        vz = u_1d ( ix, 4, 1) + 0.5 * u_1d ( ix, 4, 2)
     1     + ob_6 * u_1d ( ix, 4, 3)

        bx = bfxbdy_1d ( ix)
        by = u_1d ( ix, 6, 1) + 0.5 * u_1d ( ix, 6, 2)
     1     + ob_6 * u_1d ( ix, 6, 3)
        bz = u_1d ( ix, 7, 1) + 0.5 * u_1d ( ix, 7, 2)
     1     + ob_6 * u_1d ( ix, 7, 3)

        IF ( ( rh <= rhofloor) .OR.
     1       ( en - 0.5 * ( vx**2 + vy**2 + vz**2) / rh
     1       - ob_8pi * ( bx**2 + by**2 + bz**2) <= prsfloor) ) THEN
        u_1d ( ix, :, 2: NUMORDER) = 0.0 
        END IF

      END DO

! ----------

      DO ix = oned_begin - oned_off + ix_shrink,
     1        oned_end + oned_off - ix_shrink

        rh = u_1d ( ix, 1, 1) - 0.5 * u_1d ( ix, 1, 2)
     1     + ob_6 * u_1d ( ix, 1, 3)
        en = u_1d ( ix, 5, 1) - 0.5 * u_1d ( ix, 5, 2)
     1     + ob_6 * u_1d ( ix, 5, 3)

        vx = u_1d ( ix, 2, 1) - 0.5 * u_1d ( ix, 2, 2)
     1     + ob_6 * u_1d ( ix, 2, 3)
        vy = u_1d ( ix, 3, 1) - 0.5 * u_1d ( ix, 3, 2)
     1     + ob_6 * u_1d ( ix, 3, 3)
        vz = u_1d ( ix, 4, 1) - 0.5 * u_1d ( ix, 4, 2)
     1     + ob_6 * u_1d ( ix, 4, 3)

        bx = bfxbdy_1d ( ix - 1)
        by = u_1d ( ix, 6, 1) - 0.5 * u_1d ( ix, 6, 2)
     1     + ob_6 * u_1d ( ix, 6, 3)
        bz = u_1d ( ix, 7, 1) - 0.5 * u_1d ( ix, 7, 2)
     1     + ob_6 * u_1d ( ix, 7, 3)

        IF ( ( rh <= rhofloor) .OR.
     1       ( en - 0.5 * ( vx**2 + vy**2 + vz**2) / rh
     1       - ob_8pi * ( bx**2 + by**2 + bz**2) <= prsfloor) ) THEN
        u_1d ( ix, :, 2: NUMORDER) = 0.0
        END IF

      END DO

! ------------------------------
#endif
! End of "NUMORDER == 3".
! ------------------------------
#if ( NUMORDER == 4)
! ------------------------------

      DO ix = oned_begin - oned_off + ix_shrink,
     1        oned_end + oned_off - ix_shrink

        rh = u_1d ( ix, 1, 1) + 0.5 * u_1d ( ix, 1, 2)
     1     + ob_6 * u_1d ( ix, 1, 3) + ob_20 * u_1d ( ix, 1, 4)
        en = u_1d ( ix, 5, 1) + 0.5 * u_1d ( ix, 5, 2)
     1     + ob_6 * u_1d ( ix, 5, 3) + ob_20 * u_1d ( ix, 5, 4)

        vx = u_1d ( ix, 2, 1) + 0.5 * u_1d ( ix, 2, 2)
     1     + ob_6 * u_1d ( ix, 2, 3) + ob_20 * u_1d ( ix, 2, 4)
        vy = u_1d ( ix, 3, 1) + 0.5 * u_1d ( ix, 3, 2)
     1     + ob_6 * u_1d ( ix, 3, 3) + ob_20 * u_1d ( ix, 3, 4)
        vz = u_1d ( ix, 4, 1) + 0.5 * u_1d ( ix, 4, 2)
     1     + ob_6 * u_1d ( ix, 4, 3) + ob_20 * u_1d ( ix, 4, 4)

        bx = bfxbdy_1d ( ix)
        by = u_1d ( ix, 6, 1) + 0.5 * u_1d ( ix, 6, 2)
     1     + ob_6 * u_1d ( ix, 6, 3) + ob_20 * u_1d ( ix, 6, 4)
        bz = u_1d ( ix, 7, 1) + 0.5 * u_1d ( ix, 7, 2)
     1     + ob_6 * u_1d ( ix, 7, 3) + ob_20 * u_1d ( ix, 7, 4)

        IF ( ( rh <= rhofloor) .OR.
     1       ( en - 0.5 * ( vx**2 + vy**2 + vz**2) / rh
     1       - ob_8pi * ( bx**2 + by**2 + bz**2) <= prsfloor) ) THEN
        u_1d ( ix, :, 2: NUMORDER) = 0.0 
        END IF

      END DO

! ----------

      DO ix = oned_begin - oned_off + ix_shrink,
     1        oned_end + oned_off - ix_shrink

        rh = u_1d ( ix, 1, 1) - 0.5 * u_1d ( ix, 1, 2)
     1     + ob_6 * u_1d ( ix, 1, 3) - ob_20 * u_1d ( ix, 1, 4)
        en = u_1d ( ix, 5, 1) - 0.5 * u_1d ( ix, 5, 2)
     1     + ob_6 * u_1d ( ix, 5, 3) - ob_20 * u_1d ( ix, 5, 4)

        vx = u_1d ( ix, 2, 1) - 0.5 * u_1d ( ix, 2, 2)
     1     + ob_6 * u_1d ( ix, 2, 3) - ob_20 * u_1d ( ix, 2, 4)
        vy = u_1d ( ix, 3, 1) - 0.5 * u_1d ( ix, 3, 2)
     1     + ob_6 * u_1d ( ix, 3, 3) - ob_20 * u_1d ( ix, 3, 4)
        vz = u_1d ( ix, 4, 1) - 0.5 * u_1d ( ix, 4, 2)
     1     + ob_6 * u_1d ( ix, 4, 3) - ob_20 * u_1d ( ix, 4, 4)

        bx = bfxbdy_1d ( ix - 1)
        by = u_1d ( ix, 6, 1) - 0.5 * u_1d ( ix, 6, 2)
     1     + ob_6 * u_1d ( ix, 6, 3) - ob_20 * u_1d ( ix, 6, 4)
        bz = u_1d ( ix, 7, 1) - 0.5 * u_1d ( ix, 7, 2)
     1     + ob_6 * u_1d ( ix, 7, 3) - ob_20 * u_1d ( ix, 7, 4)

        IF ( ( rh <= rhofloor) .OR.
     1       ( en - 0.5 * ( vx**2 + vy**2 + vz**2) / rh
     1       - ob_8pi * ( bx**2 + by**2 + bz**2) <= prsfloor) ) THEN
        u_1d ( ix, :, 2: NUMORDER) = 0.0
        END IF

      END DO

! ------------------------------
#endif
! End of "NUMORDER == 4".
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE LIMIT_1D_MHD_CONS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_ADER_UPDT_PTWISE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, 

     1        u, dudt, bfxbdy, bfybdy, bfzbdy,
     1        entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy" and produces an ADER-based time
! rate of update. This is stored in
! "u ( ixx, iyy, izz, :, NUM_ZONECEN_BASES + 1: NUM_SPACETIME_BASES)" within
! each zone. In this version the fluxes are not obtained from the
! Riemann solvers but rather by direct flux evaluation at dynamically
! active quadrature points in the reference element, see Balsara et al. (2009).

! It is effectively Stage I, i.e. the predictor stage, of a two-stage
! update.

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
     1                  n_cc_components) :: dudt

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

! Variables past this point are not passed across this subroutine's interface.

      INTEGER ixx, iyy, izz, i_euler, ivar, ivar_lo, corrupted_zone

      REAL dtmove, xcuse, ycuse, zcuse, dxuse, dyuse, dzuse

      REAL, DIMENSION ( NVAR, NUM_ZONECEN_BASES) :: w_1d

      REAL, DIMENSION ( 5) :: entropy_1d

      REAL, DIMENSION ( NUM_SPACETIME_QUADPTS) ::
     1                gamma_eos_quadpt, molewt_quadpt

      REAL, DIMENSION ( NVAR, NUM_SPACETIME_BASES) ::
     1                  u_1d, f_1d, g_1d, h_1d, s_1d

      REAL, DIMENSION ( NVAR) :: s_int

! ----------------------------------------------------------------------

      dtmove = dtcur

! ----------------------------------------------------------------------

! Now we cycle through the whole mesh, choosing ADER_CG.

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, i_euler, ivar, ivar_lo,      &
!$OMP&    corrupted_zone, xcuse, dxuse, ycuse, dyuse, zcuse, dzuse,     &
!$OMP&    gamma_eos_quadpt, molewt_quadpt,                              &
!$OMP&    w_1d, entropy_1d, u_1d, f_1d, g_1d, h_1d, s_1d, s_int)        &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz, dtmove,  &
!$OMP&    igeom, gamma, gasconst,                                       &
!$OMP&    pi, rhofloor, prsfloor, smallnum,                             &
!$OMP&    ader_tolerance, max_iter_ader_cg, max_iter_ader_dg,           &
!$OMP&    xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,         &
!$OMP&    zb_save, dz_save, zc_save,                                    &
!$OMP&    gamma_eos_gr, mean_mol_wt_gr, u, entropygr,                   &
!$OMP&    dudt)                                                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

! ----------------------------------------------------------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

! Provide the coordinates of each zone.
          xcuse = xc_save ( ixx)
          dxuse = dx_save ( ixx)

          ycuse = yc_save ( iyy)
          dyuse = dy_save ( iyy)

          zcuse = zc_save ( izz)
          dzuse = dz_save ( izz)

! Provide the spatially reconstructed moments in each zone. Also provide
! EOS information and relativistic information (if appropriate).
          w_1d ( :, 1:NUM_ZONECEN_BASES) =
     1                      u ( ixx, iyy, izz, :, 1:NUM_ZONECEN_BASES)

#if ( PRESSURE_POSITIVITY == 1)
          entropy_1d ( :) = entropygr ( ixx, iyy, izz, :)
#else
          entropy_1d ( :) = 0.0
#endif

          gamma_eos_quadpt ( 1:NUM_SPACETIME_QUADPTS) =
     1                                   gamma_eos_gr ( ixx, iyy, izz)
          molewt_quadpt ( 1:NUM_SPACETIME_QUADPTS) =
     1                                 mean_mol_wt_gr ( ixx, iyy, izz)

! ----------------------------------------------------------------------

! Call "ADER_CG_MODAL_3D_TIME_EVOL_MHD" 
! which give us the space-time representation of "u".

          CALL ADER_CG_MODAL_3D_MHD_PTWISE

     1           ( dtmove, igeom, pi, rhofloor, prsfloor, smallnum,
     1             gamma, gasconst,
     1             ader_tolerance, max_iter_ader_cg,

     1             corrupted_zone,

     1             xcuse, ycuse, zcuse, dxuse, dyuse, dzuse,

     1             w_1d, entropy_1d, gamma_eos_quadpt, molewt_quadpt,

     1             u_1d, f_1d, g_1d, h_1d, s_1d, s_int)

! ----------------------------------------------------------------------

! Return the space-time modal information that has been returned by
! the ADER subroutine. Pathological zones will have their moments zeroed here.

          u ( ixx, iyy, izz, :, 1:NUM_SPACETIME_BASES) =
     1                          u_1d ( :, 1:NUM_SPACETIME_BASES)

! Return the space-time averaged source terms.
          dudt ( ixx, iyy, izz, :) = s_int ( :)

          END DO
! End "DO ixx ="

! ----------------------------------------------------------------------

        END DO
! End "DO iyy ="

      END DO
! End "DO izz ="

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_ADER_UPDT_PTWISE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_FLUX_PTWISE_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        expand_by_one_zone, dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_x, prs_flux_x,
     1        entropygr, entropy_flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy" and puts out "flux_x and prs_flux_x".
! The latter arrays hold the fluxes at the x-face.

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_x".

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! "dt_shift" is the fraction of the timestep over which we want the fluxes.
! Usually 1.0, it can be less than that for AMR fractional steps.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      REAL t1, t2, t3, ob3, ob6, ob12, ob20, ob24, ob_6_oned,
     1     a, b, ob_8pi

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

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: U_Remi, U_Repl,
     1     flux_re, F_L, F_R, U_HLL, F_HLL,
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

! Time coefficients in terms of the fractional time "dt_shift".
      t1 = dt_shift
      t2 = dt_shift**2
      t3 = dt_shift**3

! Some useful constants for evaluating states.
      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob20 = 1.0 / 20.0
      ob24 = 1.0 / 24.0

! ----------

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
!$OMP&    nriembegin, nriemend, U_Remi, U_Repl,                         &

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
!$OMP&    igeom, one_x, expand_by_one_zone,                             &
!$OMP&    dtcur, dt_shift, iHLLEM_choice, nspecies_min, nspecies_max,   &
!$OMP&    t1, t2, t3, ob3, ob6, ob12, ob20, ob24, ob_6_oned, ob_8pi,    &

!$OMP&    gamma, smallnum, rhofloor, prsfloor, pi,                      &
!$OMP&    use_approx_eos, split_prsflux_re,                             &
!$OMP&    lapidus_threshold, lapidus_coef, flatten_coef,                &
!$OMP&    gasconst, cp_specie, molewt_specie, heat_of_formation,        &

!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &

!$OMP&    x_indx_limits, bcarr_x, u,                                    &
!$OMP&    flux_x, prs_flux_x, entropygr, entropy_flux_x,                &
!$OMP&    mean_mol_wt_gr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr,  &
!$OMP&    divvelgr)                                                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

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
! Build "U_Remi" and "U_Repl" at various orders.
! ------------------------------
#if ( NUMORDER == 2)
! ----------

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 2)
     1      + 0.5 * u ( ixx, iyy, izz, :, 5) * t1

      U_Repl ( :) = u ( ixxp1, iyy, izz, :, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, :, 2)
     1      + 0.5 * u ( ixxp1, iyy, izz, :, 5) * t1

! ----------
#endif
! ------------------------------
#if ( NUMORDER == 3)
! ----------

#if ( USE_ONLY_1D_RS != 1)

! These are u[1/2,0,0,t] and u[-1/2,0,0,t] averaged over time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 2)
     1      + ob6 * u ( ixx, iyy, izz, :, 5)
     1     - ob12 * u ( ixx, iyy, izz, :, 6)
     1     - ob12 * u ( ixx, iyy, izz, :, 7)
     1      + 0.5 * u ( ixx, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 12) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 13) * t1

      U_Repl ( :) = u ( ixxp1, iyy, izz, :, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, :, 2)
     1      + ob6 * u ( ixxp1, iyy, izz, :, 5)
     1     - ob12 * u ( ixxp1, iyy, izz, :, 6)
     1     - ob12 * u ( ixxp1, iyy, izz, :, 7)
     1      + 0.5 * u ( ixxp1, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixxp1, iyy, izz, :, 12) * t2
     1     - 0.25 * u ( ixxp1, iyy, izz, :, 13) * t1

#else

! These are u[1/2,y,z,t] and u[-1/2,y,z,t] averaged over y, z and time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 2)
     1      + ob6 * u ( ixx, iyy, izz, :, 5)
     1      + 0.5 * u ( ixx, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 12) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 13) * t1


      U_Repl ( :) = u ( ixxp1, iyy, izz, :, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, :, 2)
     1      + ob6 * u ( ixxp1, iyy, izz, :, 5)
     1      + 0.5 * u ( ixxp1, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixxp1, iyy, izz, :, 12) * t2
     1     - 0.25 * u ( ixxp1, iyy, izz, :, 13) * t1

#endif

! ----------
#endif
! ------------------------------
#if ( NUMORDER == 4)
! ----------

#if ( USE_ONLY_1D_RS != 1)

! These are u[1/2,0,0,t] and u[-1/2,0,0,t] averaged over time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 2)
     1      + ob6 * u ( ixx, iyy, izz, :, 5)
     1     - ob12 * u ( ixx, iyy, izz, :, 6)
     1     - ob12 * u ( ixx, iyy, izz, :, 7)
     1     + ob20 * u ( ixx, iyy, izz, :, 11)
     1     - ob24 * u ( ixx, iyy, izz, :, 16)
     1     - ob24 * u ( ixx, iyy, izz, :, 18)
     1      + 0.5 * u ( ixx, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 23) * t3
     1     + 0.25 * u ( ixx, iyy, izz, :, 24) * t1
     1      + ob6 * u ( ixx, iyy, izz, :, 27) * t2
     1     + ob12 * u ( ixx, iyy, izz, :, 30) * t1
     1     - ob24 * u ( ixx, iyy, izz, :, 31) * t1
     1     - ob24 * u ( ixx, iyy, izz, :, 32) * t1

      U_Repl ( :) = u ( ixxp1, iyy, izz, :, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, :, 2)
     1      + ob6 * u ( ixxp1, iyy, izz, :, 5)
     1     - ob12 * u ( ixxp1, iyy, izz, :, 6)
     1     - ob12 * u ( ixxp1, iyy, izz, :, 7)
     1     - ob20 * u ( ixxp1, iyy, izz, :, 11)
     1     + ob24 * u ( ixxp1, iyy, izz, :, 16)
     1     + ob24 * u ( ixxp1, iyy, izz, :, 18)
     1      + 0.5 * u ( ixxp1, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixxp1, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixxp1, iyy, izz, :, 23) * t3
     1     - 0.25 * u ( ixxp1, iyy, izz, :, 24) * t1
     1      - ob6 * u ( ixxp1, iyy, izz, :, 27) * t2
     1     + ob12 * u ( ixxp1, iyy, izz, :, 30) * t1
     1     - ob24 * u ( ixxp1, iyy, izz, :, 31) * t1
     1     - ob24 * u ( ixxp1, iyy, izz, :, 32) * t1

#else

! These are u[1/2,y,z,t] and u[-1/2,y,z,t] averaged over y, z and time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 2)
     1      + ob6 * u ( ixx, iyy, izz, :, 5)
     1     + ob20 * u ( ixx, iyy, izz, :, 11)
     1      + 0.5 * u ( ixx, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 23) * t3
     1     + 0.25 * u ( ixx, iyy, izz, :, 24) * t1
     1      + ob6 * u ( ixx, iyy, izz, :, 27) * t2
     1     + ob12 * u ( ixx, iyy, izz, :, 30) * t1

      U_Repl ( :) = u ( ixxp1, iyy, izz, :, 1)
     1      - 0.5 * u ( ixxp1, iyy, izz, :, 2)
     1      + ob6 * u ( ixxp1, iyy, izz, :, 5)
     1     - ob20 * u ( ixxp1, iyy, izz, :, 11)
     1      + 0.5 * u ( ixxp1, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixxp1, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixxp1, iyy, izz, :, 23) * t3
     1     - 0.25 * u ( ixxp1, iyy, izz, :, 24) * t1
     1      - ob6 * u ( ixxp1, iyy, izz, :, 27) * t2
     1     + ob12 * u ( ixxp1, iyy, izz, :, 30) * t1

#endif

! ----------
#endif
! ------------------------------

      bfxre = 0.5 * ( U_Remi ( 6) + U_Repl ( 6) )

! --------------------

      rho_L = AMAX1 ( U_Remi ( 1), rhofloor)
      a = 1.0 / rho_L

      vlx_L = U_Remi ( 2) * a
      vly_L = U_Remi ( 3) * a
      vlz_L = U_Remi ( 4) * a

      bfy_L = U_Remi ( 7)
      bfz_L = U_Remi ( 8)

! ----------

      mean_mol_wt_L = mean_mol_wt_gr ( ixx, iyy, izz)
      gamma_eos_L = gamma_eos_gr ( ixx, iyy, izz)
      gamma_soundspeed_L = gamma_soundspeed_gr ( ixx, iyy, izz)

! ----------

      entropy_L = entropygr ( ixx, iyy, izz, 1)
     1    + 0.5 * entropygr ( ixx, iyy, izz, 2)
     1+ dt_shift * entropygr ( ixx, iyy, izz, 5)

      b = 0.5 * rho_L * ( vlx_L**2 + vly_L**2 + vlz_L**2)
     1  + ( bfxre**2 + bfy_L**2 + bfz_L**2) * ob_8pi

      prs_L = ( U_Remi ( 5) - b) * ( gamma_eos_L - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( prs_L .LE. prsfloor) prs_L = entropy_L
     1                              * rho_L**( gamma_eos_L - 1.0)
#endif
      prs_L = AMAX1 ( prsfloor, prs_L)

! ----------

#if ( NFLUID >= 1)
      specie_L ( 1: NFLUID) = U_Remi ( nspecies_min: nspecies_max)
#endif

! --------------------

      rho_R = AMAX1 ( U_Repl ( 1), rhofloor)
      a = 1.0 / rho_R

      vlx_R = U_Repl ( 2) * a
      vly_R = U_Repl ( 3) * a
      vlz_R = U_Repl ( 4) * a

      bfy_R = U_Repl ( 7)
      bfz_R = U_Repl ( 8)

! ----------

      mean_mol_wt_R = mean_mol_wt_gr ( ixxp1, iyy, izz)
      gamma_eos_R = gamma_eos_gr ( ixxp1, iyy, izz)
      gamma_soundspeed_R = gamma_soundspeed_gr ( ixxp1, iyy, izz)

! ----------

      entropy_R = entropygr ( ixxp1, iyy, izz, 1)
     1    - 0.5 * entropygr ( ixxp1, iyy, izz, 2)
     1+ dt_shift * entropygr ( ixxp1, iyy, izz, 5)

      b = 0.5 * rho_R * ( vlx_R**2 + vly_R**2 + vlz_R**2)
     1  + ( bfxre**2 + bfy_R**2 + bfz_R**2) * ob_8pi

      prs_R = ( U_Repl ( 5) - b) * ( gamma_eos_R - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( prs_R .LE. prsfloor) prs_R = entropy_R
     1                              * rho_R**( gamma_eos_R - 1.0)
#endif
      prs_R = AMAX1 ( prsfloor, prs_R)

! ----------

#if ( NFLUID >= 1)
      specie_R ( 1: NFLUID) = U_Repl ( nspecies_min: nspecies_max)
#endif

! ------------------------------

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

          IF ( USE_ONLY_1D_RS == 1) THEN

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

          ELSE IF ( ( ixx == nriembegin + one_x) .AND.
     1                ( bcarr_x ( 1, iyy, izz) == 3) ) THEN

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
     1                ( bcarr_x ( 2, iyy, izz) == 3) ) THEN

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

      END SUBROUTINE MAKE_FLUX_PTWISE_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_FLUX_PTWISE_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        expand_by_one_zone, dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_y, prs_flux_y,
     1        entropygr, entropy_flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy" and puts out "flux_y and prs_flux_y".
! The latter arrays hold the fluxes at the y-faces.

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_y".

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! "dt_shift" is the fraction of the timestep over which we want the fluxes.
! Usually 1.0, it can be less than that for AMR fractional steps.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      REAL t1, t2, t3, ob3, ob6, ob12, ob20, ob24, ob_6_oned,
     1     a, b, ob_8pi

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

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: U_Remi, U_Repl,
     1     flux_re, F_L, F_R, U_HLL, F_HLL,
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

! Time coefficients in terms of the fractional time "dt_shift".
      t1 = dt_shift
      t2 = dt_shift**2
      t3 = dt_shift**3

! Some useful constants for evaluating states.
      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob20 = 1.0 / 20.0
      ob24 = 1.0 / 24.0

! ----------

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
!$OMP&    nriembegin, nriemend, U_Remi, U_Repl,                         &

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
!$OMP&    igeom, one_y, expand_by_one_zone,                             &
!$OMP&    dtcur, dt_shift, iHLLEM_choice, nspecies_min, nspecies_max,   &
!$OMP&    t1, t2, t3, ob3, ob6, ob12, ob20, ob24, ob_6_oned, ob_8pi,    &

!$OMP&    gamma, smallnum, rhofloor, prsfloor, pi,                      &
!$OMP&    use_approx_eos, split_prsflux_re,                             &
!$OMP&    lapidus_threshold, lapidus_coef, flatten_coef,                &
!$OMP&    gasconst, cp_specie, molewt_specie, heat_of_formation,        &

!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &

!$OMP&    y_indx_limits, bcarr_y, u,                                    &
!$OMP&    flux_y, prs_flux_y, entropygr, entropy_flux_y,                &
!$OMP&    mean_mol_wt_gr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr,  &
!$OMP&    divvelgr)                                                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

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
! Build "U_Remi" and "U_Repl" at various orders.
! ------------------------------
#if ( NUMORDER == 2)
! ----------

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 3)
     1      + 0.5 * u ( ixx, iyy, izz, :, 5) * t1

      U_Repl ( :) = u ( ixx, iyyp1, izz, :, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, :, 3)
     1      + 0.5 * u ( ixx, iyyp1, izz, :, 5) * t1

! ----------
#endif
! ------------------------------
#if ( NUMORDER == 3)
! ----------

#if ( USE_ONLY_1D_RS != 1)

! These are u[0,1/2,0,t] and u[0,-1/2,0,t] averaged over time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 3)
     1     - ob12 * u ( ixx, iyy, izz, :, 5)
     1      + ob6 * u ( ixx, iyy, izz, :, 6)
     1     - ob12 * u ( ixx, iyy, izz, :, 7)
     1      + 0.5 * u ( ixx, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 12) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 14) * t1

      U_Repl ( :) = u ( ixx, iyyp1, izz, :, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, :, 3)
     1     - ob12 * u ( ixx, iyyp1, izz, :, 5)
     1      + ob6 * u ( ixx, iyyp1, izz, :, 6)
     1     - ob12 * u ( ixx, iyyp1, izz, :, 7)
     1      + 0.5 * u ( ixx, iyyp1, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyyp1, izz, :, 12) * t2
     1     - 0.25 * u ( ixx, iyyp1, izz, :, 14) * t1

#else

! These are u[x,1/2,z,t] and u[x,-1/2,z,t] averaged over x, z and time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 3)
     1      + ob6 * u ( ixx, iyy, izz, :, 6)
     1      + 0.5 * u ( ixx, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 12) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 14) * t1

      U_Repl ( :) = u ( ixx, iyyp1, izz, :, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, :, 3)
     1      + ob6 * u ( ixx, iyyp1, izz, :, 6)
     1      + 0.5 * u ( ixx, iyyp1, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyyp1, izz, :, 12) * t2
     1     - 0.25 * u ( ixx, iyyp1, izz, :, 14) * t1

#endif

! ----------
#endif
! ------------------------------
#if ( NUMORDER == 4)
! ----------

#if ( USE_ONLY_1D_RS != 1)

! These are u[0,1/2,0,t] and u[0,-1/2,0,t] averaged over time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 3)
     1     - ob12 * u ( ixx, iyy, izz, :, 5)
     1      + ob6 * u ( ixx, iyy, izz, :, 6)
     1     - ob12 * u ( ixx, iyy, izz, :, 7)
     1     + ob20 * u ( ixx, iyy, izz, :, 12)
     1     - ob24 * u ( ixx, iyy, izz, :, 14)
     1     - ob24 * u ( ixx, iyy, izz, :, 19)
     1      + 0.5 * u ( ixx, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 23) * t3
     1     + 0.25 * u ( ixx, iyy, izz, :, 25) * t1
     1      + ob6 * u ( ixx, iyy, izz, :, 28) * t2
     1     - ob24 * u ( ixx, iyy, izz, :, 30) * t1
     1     + ob12 * u ( ixx, iyy, izz, :, 31) * t1
     1     - ob24 * u ( ixx, iyy, izz, :, 32) * t1

      U_Repl ( :) = u ( ixx, iyyp1, izz, :, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, :, 3)
     1     - ob12 * u ( ixx, iyyp1, izz, :, 5)
     1      + ob6 * u ( ixx, iyyp1, izz, :, 6)
     1     - ob12 * u ( ixx, iyyp1, izz, :, 7)
     1     - ob20 * u ( ixx, iyyp1, izz, :, 12)
     1     + ob24 * u ( ixx, iyyp1, izz, :, 14)
     1     + ob24 * u ( ixx, iyyp1, izz, :, 19)
     1      + 0.5 * u ( ixx, iyyp1, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyyp1, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyyp1, izz, :, 23) * t3
     1     - 0.25 * u ( ixx, iyyp1, izz, :, 25) * t1
     1      - ob6 * u ( ixx, iyyp1, izz, :, 28) * t2
     1     - ob24 * u ( ixx, iyyp1, izz, :, 30) * t1
     1     + ob12 * u ( ixx, iyyp1, izz, :, 31) * t1
     1     - ob24 * u ( ixx, iyyp1, izz, :, 32) * t1

#else

! These are u[x,1/2,z,t] and u[x,-1/2,z,t] averaged over x, z and time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 3)
     1      + ob6 * u ( ixx, iyy, izz, :, 6)
     1     + ob20 * u ( ixx, iyy, izz, :, 12)
     1      + 0.5 * u ( ixx, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 23) * t3
     1     + 0.25 * u ( ixx, iyy, izz, :, 25) * t1
     1      + ob6 * u ( ixx, iyy, izz, :, 28) * t2
     1     + ob12 * u ( ixx, iyy, izz, :, 31) * t1

      U_Repl ( :) = u ( ixx, iyyp1, izz, :, 1)
     1      - 0.5 * u ( ixx, iyyp1, izz, :, 3)
     1      + ob6 * u ( ixx, iyyp1, izz, :, 6)
     1     - ob20 * u ( ixx, iyyp1, izz, :, 12)
     1      + 0.5 * u ( ixx, iyyp1, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyyp1, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyyp1, izz, :, 23) * t3
     1     - 0.25 * u ( ixx, iyyp1, izz, :, 25) * t1
     1      - ob6 * u ( ixx, iyyp1, izz, :, 28) * t2
     1     + ob12 * u ( ixx, iyyp1, izz, :, 31) * t1

#endif

! ----------
#endif
! ------------------------------

      bfyre = 0.5 * ( U_Remi ( 7) + U_Repl ( 7) )

! --------------------

      rho_L = AMAX1 ( U_Remi ( 1), rhofloor)
      a = 1.0 / rho_L

      vlx_L = U_Remi ( 2) * a
      vly_L = U_Remi ( 3) * a
      vlz_L = U_Remi ( 4) * a

      bfx_L = U_Remi ( 6)
      bfz_L = U_Remi ( 8)

! ----------

      mean_mol_wt_L = mean_mol_wt_gr ( ixx, iyy, izz)
      gamma_eos_L = gamma_eos_gr ( ixx, iyy, izz)
      gamma_soundspeed_L = gamma_soundspeed_gr ( ixx, iyy, izz)

! ----------

      entropy_L = entropygr ( ixx, iyy, izz, 1)
     1    + 0.5 * entropygr ( ixx, iyy, izz, 3)
     1+ dt_shift * entropygr ( ixx, iyy, izz, 5)

      b = 0.5 * rho_L * ( vlx_L**2 + vly_L**2 + vlz_L**2)
     1  + ( bfx_L**2 + bfyre**2 + bfz_L**2) * ob_8pi

      prs_L = ( U_Remi ( 5) - b) * ( gamma_eos_L - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( prs_L .LE. prsfloor) prs_L = entropy_L
     1                              * rho_L**( gamma_eos_L - 1.0)
#endif
      prs_L = AMAX1 ( prsfloor, prs_L)

! ----------

#if ( NFLUID >= 1)
      specie_L ( 1: NFLUID) = U_Remi ( nspecies_min: nspecies_max)
#endif

! --------------------

      rho_R = AMAX1 ( U_Repl ( 1), rhofloor)
      a = 1.0 / rho_R

      vlx_R = U_Repl ( 2) * a
      vly_R = U_Repl ( 3) * a
      vlz_R = U_Repl ( 4) * a

      bfx_R = U_Repl ( 6)
      bfz_R = U_Repl ( 8)

! ----------

      mean_mol_wt_R = mean_mol_wt_gr ( ixx, iyyp1, izz)
      gamma_eos_R = gamma_eos_gr ( ixx, iyyp1, izz)
      gamma_soundspeed_R = gamma_soundspeed_gr ( ixx, iyyp1, izz)

! ----------

      entropy_R = entropygr ( ixx, iyyp1, izz, 1)
     1    - 0.5 * entropygr ( ixx, iyyp1, izz, 3)
     1+ dt_shift * entropygr ( ixx, iyyp1, izz, 5)

      b = 0.5 * rho_R * ( vlx_R**2 + vly_R**2 + vlz_R**2)
     1  + ( bfx_R**2 + bfyre**2 + bfz_R**2) * ob_8pi

      prs_R = ( U_Repl ( 5) - b) * ( gamma_eos_R - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( prs_R .LE. prsfloor) prs_R = entropy_R
     1                              * rho_R**( gamma_eos_R - 1.0)
#endif
      prs_R = AMAX1 ( prsfloor, prs_R)

! ----------

#if ( NFLUID >= 1)
      specie_R ( 1: NFLUID) = U_Repl ( nspecies_min: nspecies_max)
#endif

! ------------------------------

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

          IF ( USE_ONLY_1D_RS == 1) THEN

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

          ELSE IF ( ( iyy == nriembegin + one_y) .AND.
     1                ( bcarr_y ( ixx, 1, izz) == 3) ) THEN

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
     1                ( bcarr_y ( ixx, 2, izz) == 3) ) THEN

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

      END SUBROUTINE MAKE_FLUX_PTWISE_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_FLUX_PTWISE_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,

     1        expand_by_one_zone, dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy,
     1        flux_z, prs_flux_z,
     1        entropygr, entropy_flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr, bfxgr, bfygr, bfzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u, bfxbdy, bfybdy, bfzbdy" and puts out "flux_z and prs_flux_z".
! The latter arrays hold the fluxes at the z-faces.

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_z".

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! "dt_shift" is the fraction of the timestep over which we want the fluxes.
! Usually 1.0, it can be less than that for AMR fractional steps.

! This subroutine is effectively Stage II, i.e. the corrector stage,
! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone

      REAL dtcur, dt_shift

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_fc_modes) :: bfxbdy, bfybdy, bfzbdy

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

      REAL t1, t2, t3, ob3, ob6, ob12, ob20, ob24, ob_6_oned,
     1     a, b, ob_8pi

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

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

      REAL, DIMENSION ( NFIELD + NFLUID + 1) :: U_Remi, U_Repl,
     1     flux_re, F_L, F_R, U_HLL, F_HLL,
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

! Time coefficients in terms of the fractional time "dt_shift".
      t1 = dt_shift
      t2 = dt_shift**2
      t3 = dt_shift**3

! Some useful constants for evaluating states.
      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob20 = 1.0 / 20.0
      ob24 = 1.0 / 24.0

! ----------

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
!$OMP&    nriembegin, nriemend, U_Remi, U_Repl,                         &

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
!$OMP&    igeom, one_z, expand_by_one_zone,                             &
!$OMP&    dtcur, dt_shift, iHLLEM_choice, nspecies_min, nspecies_max,   &
!$OMP&    t1, t2, t3, ob3, ob6, ob12, ob20, ob24, ob_6_oned, ob_8pi,    &

!$OMP&    gamma, smallnum, rhofloor, prsfloor, pi,                      &
!$OMP&    use_approx_eos, split_prsflux_re,                             &
!$OMP&    lapidus_threshold, lapidus_coef, flatten_coef,                &
!$OMP&    gasconst, cp_specie, molewt_specie, heat_of_formation,        &

!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &

!$OMP&    z_indx_limits, bcarr_z, u,                                    &
!$OMP&    flux_z, prs_flux_z, entropygr, entropy_flux_z,                &
!$OMP&    mean_mol_wt_gr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr,  &
!$OMP&    divvelgr)                                                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

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
! Build "U_Remi" and "U_Repl" at various orders.
! ------------------------------
#if ( NUMORDER == 2)
! ----------

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 4)
     1      + 0.5 * u ( ixx, iyy, izz, :, 5) * t1

      U_Repl ( :) = u ( ixx, iyy, izzp1, :, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, :, 4)
     1      + 0.5 * u ( ixx, iyy, izzp1, :, 5) * t1

! ----------
#endif
! ------------------------------
#if ( NUMORDER == 3)
! ----------

#if ( USE_ONLY_1D_RS != 1)

! These are u[0,0,1/2,t] and u[0,0,-1/2,t] averaged over time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 4)
     1     - ob12 * u ( ixx, iyy, izz, :, 5)
     1     - ob12 * u ( ixx, iyy, izz, :, 6)
     1      + ob6 * u ( ixx, iyy, izz, :, 7)
     1      + 0.5 * u ( ixx, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 12) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 15) * t1

      U_Repl ( :) = u ( ixx, iyy, izzp1, :, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, :, 4)
     1     - ob12 * u ( ixx, iyy, izzp1, :, 5)
     1     - ob12 * u ( ixx, iyy, izzp1, :, 6)
     1      + ob6 * u ( ixx, iyy, izzp1, :, 7)
     1      + 0.5 * u ( ixx, iyy, izzp1, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izzp1, :, 12) * t2
     1     - 0.25 * u ( ixx, iyy, izzp1, :, 15) * t1

#else

! These are u[x,y,1/2,t] and u[x,y,-1/2,t] averaged over x, y and time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 4)
     1      + ob6 * u ( ixx, iyy, izz, :, 7)
     1      + 0.5 * u ( ixx, iyy, izz, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 12) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 15) * t1

      U_Repl ( :) = u ( ixx, iyy, izzp1, :, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, :, 4)
     1      + ob6 * u ( ixx, iyy, izzp1, :, 7)
     1      + 0.5 * u ( ixx, iyy, izzp1, :, 11) * t1
     1      + ob3 * u ( ixx, iyy, izzp1, :, 12) * t2
     1     - 0.25 * u ( ixx, iyy, izzp1, :, 15) * t1

#endif

! ----------
#endif
! ------------------------------
#if ( NUMORDER == 4)
! ----------

#if ( USE_ONLY_1D_RS != 1)

! These are u[0,0,1/2,t] and u[0,0,-1/2,t] averaged over time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 4)
     1     - ob12 * u ( ixx, iyy, izz, :, 5)
     1     - ob12 * u ( ixx, iyy, izz, :, 6)
     1      + ob6 * u ( ixx, iyy, izz, :, 7)
     1     + ob20 * u ( ixx, iyy, izz, :, 13)
     1     - ob24 * u ( ixx, iyy, izz, :, 15)
     1     - ob24 * u ( ixx, iyy, izz, :, 17)
     1      + 0.5 * u ( ixx, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 23) * t3
     1     + 0.25 * u ( ixx, iyy, izz, :, 26) * t1
     1      + ob6 * u ( ixx, iyy, izz, :, 29) * t2
     1     - ob24 * u ( ixx, iyy, izz, :, 30) * t1
     1     - ob24 * u ( ixx, iyy, izz, :, 31) * t1
     1     + ob12 * u ( ixx, iyy, izz, :, 32) * t1

      U_Repl ( :) = u ( ixx, iyy, izzp1, :, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, :, 4)
     1     - ob12 * u ( ixx, iyy, izzp1, :, 5)
     1     - ob12 * u ( ixx, iyy, izzp1, :, 6)
     1      + ob6 * u ( ixx, iyy, izzp1, :, 7)
     1     - ob20 * u ( ixx, iyy, izzp1, :, 13)
     1     + ob24 * u ( ixx, iyy, izzp1, :, 15)
     1     + ob24 * u ( ixx, iyy, izzp1, :, 17)
     1      + 0.5 * u ( ixx, iyy, izzp1, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izzp1, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izzp1, :, 23) * t3
     1     - 0.25 * u ( ixx, iyy, izzp1, :, 26) * t1
     1      - ob6 * u ( ixx, iyy, izzp1, :, 29) * t2
     1     - ob24 * u ( ixx, iyy, izzp1, :, 30) * t1
     1     - ob24 * u ( ixx, iyy, izzp1, :, 31) * t1
     1     + ob12 * u ( ixx, iyy, izzp1, :, 32) * t1

#else

! These are u[x,y,1/2,t] and u[x,y,-1/2,t] averaged over x, y and time.

      U_Remi ( :) = u ( ixx, iyy, izz, :, 1)
     1      + 0.5 * u ( ixx, iyy, izz, :, 4)
     1      + ob6 * u ( ixx, iyy, izz, :, 7)
     1     + ob20 * u ( ixx, iyy, izz, :, 13)
     1      + 0.5 * u ( ixx, iyy, izz, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izz, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izz, :, 23) * t3
     1     + 0.25 * u ( ixx, iyy, izz, :, 26) * t1
     1      + ob6 * u ( ixx, iyy, izz, :, 29) * t2
     1     + ob12 * u ( ixx, iyy, izz, :, 32) * t1

      U_Repl ( :) = u ( ixx, iyy, izzp1, :, 1)
     1      - 0.5 * u ( ixx, iyy, izzp1, :, 4)
     1      + ob6 * u ( ixx, iyy, izzp1, :, 7)
     1     - ob20 * u ( ixx, iyy, izzp1, :, 13)
     1      + 0.5 * u ( ixx, iyy, izzp1, :, 21) * t1
     1      + ob3 * u ( ixx, iyy, izzp1, :, 22) * t2
     1     + 0.25 * u ( ixx, iyy, izzp1, :, 23) * t3
     1     - 0.25 * u ( ixx, iyy, izzp1, :, 26) * t1
     1      - ob6 * u ( ixx, iyy, izzp1, :, 29) * t2
     1     + ob12 * u ( ixx, iyy, izzp1, :, 32) * t1

#endif

! ----------
#endif
! ------------------------------

      bfzre = 0.5 * ( U_Remi ( 8) + U_Repl ( 8) )

! --------------------

      rho_L = AMAX1 ( U_Remi ( 1), rhofloor)
      a = 1.0 / rho_L

      vlx_L = U_Remi ( 2) * a
      vly_L = U_Remi ( 3) * a
      vlz_L = U_Remi ( 4) * a

      bfx_L = U_Remi ( 6)
      bfy_L = U_Remi ( 7)

! ----------

      mean_mol_wt_L = mean_mol_wt_gr ( ixx, iyy, izz)
      gamma_eos_L = gamma_eos_gr ( ixx, iyy, izz)
      gamma_soundspeed_L = gamma_soundspeed_gr ( ixx, iyy, izz)

! ----------

      entropy_L = entropygr ( ixx, iyy, izz, 1)
     1    + 0.5 * entropygr ( ixx, iyy, izz, 4)
     1+ dt_shift * entropygr ( ixx, iyy, izz, 5)

      b = 0.5 * rho_L * ( vlx_L**2 + vly_L**2 + vlz_L**2)
     1  + ( bfx_L**2 + bfy_L**2 + bfzre**2) * ob_8pi

      prs_L = ( U_Remi ( 5) - b) * ( gamma_eos_L - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( prs_L .LE. prsfloor) prs_L = entropy_L
     1                              * rho_L**( gamma_eos_L - 1.0)
#endif
      prs_L = AMAX1 ( prsfloor, prs_L)

! ----------

#if ( NFLUID >= 1)
      specie_L ( 1: NFLUID) = U_Remi ( nspecies_min: nspecies_max)
#endif

! --------------------

      rho_R = AMAX1 ( U_Repl ( 1), rhofloor)
      a = 1.0 / rho_R

      vlx_R = U_Repl ( 2) * a
      vly_R = U_Repl ( 3) * a
      vlz_R = U_Repl ( 4) * a

      bfx_R = U_Repl ( 6)
      bfy_R = U_Repl ( 7)

! ----------

      mean_mol_wt_R = mean_mol_wt_gr ( ixx, iyy, izzp1)
      gamma_eos_R = gamma_eos_gr ( ixx, iyy, izzp1)
      gamma_soundspeed_R = gamma_soundspeed_gr ( ixx, iyy, izzp1)

! ----------

      entropy_R = entropygr ( ixx, iyy, izzp1, 1)
     1    - 0.5 * entropygr ( ixx, iyy, izzp1, 4)
     1+ dt_shift * entropygr ( ixx, iyy, izzp1, 5)

      b = 0.5 * rho_R * ( vlx_R**2 + vly_R**2 + vlz_R**2)
     1  + ( bfx_R**2 + bfy_R**2 + bfzre**2) * ob_8pi

      prs_R = ( U_Repl ( 5) - b) * ( gamma_eos_R - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( prs_R .LE. prsfloor) prs_R = entropy_R
     1                              * rho_R**( gamma_eos_R - 1.0)
#endif
      prs_R = AMAX1 ( prsfloor, prs_R)

! ----------

#if ( NFLUID >= 1)
      specie_R ( 1: NFLUID) = U_Repl ( nspecies_min: nspecies_max)
#endif

! ------------------------------

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

          IF ( USE_ONLY_1D_RS == 1) THEN

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

          ELSE IF ( ( izz == nriembegin + one_z) .AND.
     1                ( bcarr_z ( ixx, iyy, 1) == 3) ) THEN

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
     1                ( bcarr_z ( ixx, iyy, 2) == 3) ) THEN

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

      END SUBROUTINE MAKE_FLUX_PTWISE_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE APPLY_MULTID_RS_X_EDGE

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        flux_y, flux_z, flux_temp, efxgr,

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

! Since this is an ADER-style update, it is suggested that the update
! be done in two stages:

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
     1        n_cc_components, n_cc_modes, n_fc_modes,
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

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_y, flux_z, flux_temp

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) :: efxgr

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

      REAL ob_6_multid, ob_6_oned,
     1     t1, t2, t3, ob3, ob6, ob12, ob20, ob24

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + 1 + NFLUID) ::

     1      U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn,

     1      U_Rt_Up_Pt, U_Lf_Up_Pt, U_Lf_Dn_Pt, U_Rt_Dn_Pt

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

! In true 2d, only the x-edges contribute to a problem whose variation
! is entirely in the yz-plane.

      IF ( ( TRUE_2D == 1) .AND. ( ioffx /= 0) ) RETURN

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

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

! ----------

! Time coefficients in terms of the fractional time "dt_shift".
      t1 = dt_shift
      t2 = dt_shift**2
      t3 = dt_shift**3

! Some useful constants for evaluating states.
      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob20 = 1.0 / 20.0
      ob24 = 1.0 / 24.0

! ----------------------------------------------------------------------

        multid_wavemodel_num = 1  ! 1 includes supersonic; 2 is fully subsonic.

! "4 * ob_6_multid + 4 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux

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
!$OMP&    U_Rt_Up_Pt, U_Lf_Up_Pt, U_Lf_Dn_Pt, U_Rt_Dn_Pt,               &
!$OMP&    entropy_Rt_Up, entropy_Lf_Up, entropy_Lf_Dn, entropy_Rt_Dn,   &

!$OMP&    gamma_eos_Rt_Up, gamma_eos_Lf_Up, gamma_eos_Lf_Dn,            &
!$OMP&    gamma_eos_Rt_Dn,                                              &
!$OMP&    gamma_sdsp_Rt_Up, gamma_sdsp_Lf_Up, gamma_sdsp_Lf_Dn,         &
!$OMP&    gamma_sdsp_Rt_Dn,                                             &

!$OMP&    blend_HLL, my_Prs_Gas_Star,                                   &
!$OMP&    my_blend_HLL, my_blend_MDHLL, my_Prs_Tot_Star,                &
!$OMP&    U_vrtx, F_vrtx, G_vrtx, gamma_eos, gamma_sdsp,                &
!$OMP&    U_Star, F_Star, G_Star,                                       &
!$OMP&    U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,     &
!$OMP&    entropy_vrtx, F_entropy_flux_x, G_entropy_flux_y,             &

!$OMP&    temp_01, temp_02, temp_03, temp_04, temp_05)                  &

!$OMP&            SHARED ( ix1, iy1, iz1, ioffx, ioffy, ioffz,          &
!$OMP&    expand_by_one_zone, ioffx_expand, ioffy_expand, ioffz_expand, &
!$OMP&    ob_6_multid, ob_6_oned, blend_mdhll_to_mdhllc,                &
!$OMP&    dt_shift, smallnum, prsfloor, rhofloor, pi, idirection,       &
!$OMP&    t1, t2, t3, ob3, ob6, ob12, ob20, ob24,                       &

!$OMP&    nzone_vrtx, Vlx_ALE, Vly_ALE, nx, ny,                         &
!$OMP&    multid_wavemodel_num, split_prsflux_re,                       &

!$OMP&    x_indx_limits, y_indx_limits, z_indx_limits,                  &
!$OMP&    bcarr_x, bcarr_y, bcarr_z, dx_save, dy_save, dz_save,         &
!$OMP&    u, entropygr, gamma_eos_gr, gamma_soundspeed_gr, flattengr,   &
!$OMP&    efxgr, flux_y, flux_z, flux_temp)                             &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz_expand, iz1 + ioffz_expand

      izzm1 = MAX0 ( izz - 1, 0 - ioffz)
      izz0 = MAX0 ( izz, 1 - ioffz)
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)

        DO iyy = 0 - ioffy_expand, iy1 + ioffy_expand

        iyym1 = MAX0 ( iyy - 1, 0 - ioffy)
        iyy0 = MAX0 ( iyy, 1 - ioffy)
        iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)

! ----------------------------------------------------------------------
! Load in the four states that surround the x-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------
! Indices ( iyy0, izz0) always map to an interior point on the mesh including
! its ghost zones. This is needed because ( iyy, izz) map to the edges.

          nriembegin = x_indx_limits ( 1, iyy0, izz0) - ioffx_expand
          nriemend = x_indx_limits ( 2, iyy0, izz0) + ioffx_expand

          DO ixx = nriembegin, nriemend

! ------------------------------
! First, we build the four states at various orders.
! ------------------------------
#if ( NUMORDER == 2)
! --------------------

          U_Lf_Dn ( ixx, :) = u ( ixx, iyy0, izz0, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 4)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 5) * t1

          U_Lf_Dn_Pt ( ixx, :) = U_Lf_Dn ( ixx, :)

! ----------

          U_Rt_Dn ( ixx, :) = u ( ixx, iyyp1, izz0, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 4)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 5) * t1

          U_Rt_Dn_Pt ( ixx, :) = U_Rt_Dn ( ixx, :)

! ----------

          U_Lf_Up ( ixx, :) = u ( ixx, iyy0, izzp1, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyy0, izzp1, :, 4)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 5) * t1

          U_Lf_Up_Pt ( ixx, :) = U_Lf_Up ( ixx, :)

! ----------

          U_Rt_Up ( ixx, :) = u ( ixx, iyyp1, izzp1, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 4)
     1                + 0.5 * u ( ixx, iyyp1, izzp1, :, 5) * t1

          U_Rt_Up_Pt ( ixx, :) = U_Rt_Up ( ixx, :)

! --------------------
#endif
! ------------------------------
#if ( NUMORDER == 3)
! --------------------

          U_Lf_Dn ( ixx, :) = u ( ixx, iyy0, izz0, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 4)
     1                + ob6 * u ( ixx, iyy0, izz0, :, 6)
     1                + ob6 * u ( ixx, iyy0, izz0, :, 7)
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 9)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 11)
     1           * t1 + ob3 * u ( ixx, iyy0, izz0, :, 12) * t2
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 14) * t1
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 15) * t1

          U_Lf_Dn_Pt ( ixx, :) = u ( ixx, iyy0, izz0, :, 1)
     1                   + 0.5 * u ( ixx, iyy0, izz0, :, 3)
     1                   + 0.5 * u ( ixx, iyy0, izz0, :, 4)
     1                  - ob12 * u ( ixx, iyy0, izz0, :, 5)
     1                   + ob6 * u ( ixx, iyy0, izz0, :, 6)
     1                   + ob6 * u ( ixx, iyy0, izz0, :, 7)
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 9)
     1                   + 0.5 * u ( ixx, iyy0, izz0, :, 11) * t1
     1                   + ob3 * u ( ixx, iyy0, izz0, :, 12) * t2
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 14) * t1
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 15) * t1

! ----------

          U_Rt_Dn ( ixx, :) = u ( ixx, iyyp1, izz0, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 4)
     1                + ob6 * u ( ixx, iyyp1, izz0, :, 6)
     1                + ob6 * u ( ixx, iyyp1, izz0, :, 7)
     1               - 0.25 * u ( ixx, iyyp1, izz0, :, 9)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 11) * t1
     1                + ob3 * u ( ixx, iyyp1, izz0, :, 12) * t2
     1               - 0.25 * u ( ixx, iyyp1, izz0, :, 14) * t1
     1               + 0.25 * u ( ixx, iyyp1, izz0, :, 15) * t1

          U_Rt_Dn_Pt ( ixx, :) = u ( ixx, iyyp1, izz0, :, 1)
     1                   - 0.5 * u ( ixx, iyyp1, izz0, :, 3)
     1                   + 0.5 * u ( ixx, iyyp1, izz0, :, 4)
     1                  - ob12 * u ( ixx, iyyp1, izz0, :, 5)
     1                   + ob6 * u ( ixx, iyyp1, izz0, :, 6)
     1                   + ob6 * u ( ixx, iyyp1, izz0, :, 7)
     1                  - 0.25 * u ( ixx, iyyp1, izz0, :, 9)
     1                   + 0.5 * u ( ixx, iyyp1, izz0, :, 11) * t1
     1                   + ob3 * u ( ixx, iyyp1, izz0, :, 12) * t2
     1                  - 0.25 * u ( ixx, iyyp1, izz0, :, 14) * t1
     1                  + 0.25 * u ( ixx, iyyp1, izz0, :, 15) * t1

! ----------

          U_Lf_Up ( ixx, :) = u ( ixx, iyy0, izzp1, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyy0, izzp1, :, 4)
     1                + ob6 * u ( ixx, iyy0, izzp1, :, 6)
     1                + ob6 * u ( ixx, iyy0, izzp1, :, 7)
     1               - 0.25 * u ( ixx, iyy0, izzp1, :, 9)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 11) * t1
     1                + ob3 * u ( ixx, iyy0, izzp1, :, 12) * t2
     1               + 0.25 * u ( ixx, iyy0, izzp1, :, 14) * t1
     1               - 0.25 * u ( ixx, iyy0, izzp1, :, 15) * t1

          U_Lf_Up_Pt ( ixx, :) = u ( ixx, iyy0, izzp1, :, 1)
     1                   + 0.5 * u ( ixx, iyy0, izzp1, :, 3)
     1                   - 0.5 * u ( ixx, iyy0, izzp1, :, 4)
     1                  - ob12 * u ( ixx, iyy0, izzp1, :, 5)
     1                   + ob6 * u ( ixx, iyy0, izzp1, :, 6)
     1                   + ob6 * u ( ixx, iyy0, izzp1, :, 7)
     1                  - 0.25 * u ( ixx, iyy0, izzp1, :, 9)
     1                   + 0.5 * u ( ixx, iyy0, izzp1, :, 11) * t1
     1                   + ob3 * u ( ixx, iyy0, izzp1, :, 12) * t2
     1                  + 0.25 * u ( ixx, iyy0, izzp1, :, 14) * t1
     1                  - 0.25 * u ( ixx, iyy0, izzp1, :, 15) * t1

! ----------

          U_Rt_Up ( ixx, :) = u ( ixx, iyyp1, izzp1, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 4)
     1                + ob6 * u ( ixx, iyyp1, izzp1, :, 6)
     1                + ob6 * u ( ixx, iyyp1, izzp1, :, 7)
     1               + 0.25 * u ( ixx, iyyp1, izzp1, :, 9)
     1                + 0.5 * u ( ixx, iyyp1, izzp1, :, 11) * t1
     1                + ob3 * u ( ixx, iyyp1, izzp1, :, 12) * t2
     1               - 0.25 * u ( ixx, iyyp1, izzp1, :, 14) * t1
     1               - 0.25 * u ( ixx, iyyp1, izzp1, :, 15) * t1

          U_Rt_Up_Pt ( ixx, :) = u ( ixx, iyyp1, izzp1, :, 1)
     1                   - 0.5 * u ( ixx, iyyp1, izzp1, :, 3)
     1                   - 0.5 * u ( ixx, iyyp1, izzp1, :, 4)
     1                  - ob12 * u ( ixx, iyyp1, izzp1, :, 5)
     1                   + ob6 * u ( ixx, iyyp1, izzp1, :, 6)
     1                   + ob6 * u ( ixx, iyyp1, izzp1, :, 7)
     1                  + 0.25 * u ( ixx, iyyp1, izzp1, :, 9)
     1                   + 0.5 * u ( ixx, iyyp1, izzp1, :, 11) * t1
     1                   + ob3 * u ( ixx, iyyp1, izzp1, :, 12) * t2
     1                  - 0.25 * u ( ixx, iyyp1, izzp1, :, 14) * t1
     1                  - 0.25 * u ( ixx, iyyp1, izzp1, :, 15) * t1

! --------------------
#endif
! ------------------------------
#if ( NUMORDER == 4)
! --------------------

          U_Lf_Dn ( ixx, :) = u ( ixx, iyy0, izz0, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 4)
     1                + ob6 * u ( ixx, iyy0, izz0, :, 6)
     1                + ob6 * u ( ixx, iyy0, izz0, :, 7)
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 9)
     1               + ob20 * u ( ixx, iyy0, izz0, :, 12)
     1               + ob20 * u ( ixx, iyy0, izz0, :, 13)
     1               + ob12 * u ( ixx, iyy0, izz0, :, 17)
     1               + ob12 * u ( ixx, iyy0, izz0, :, 19)
     1                + 0.5 * u ( ixx, iyy0, izz0, :, 21) * t1
     1                + ob3 * u ( ixx, iyy0, izz0, :, 22) * t2
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 23) * t3
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 25) * t1
     1               + 0.25 * u ( ixx, iyy0, izz0, :, 26) * t1
     1                + ob6 * u ( ixx, iyy0, izz0, :, 28) * t2
     1                + ob6 * u ( ixx, iyy0, izz0, :, 29) * t2
     1               + ob12 * u ( ixx, iyy0, izz0, :, 31) * t1
     1               + ob12 * u ( ixx, iyy0, izz0, :, 32) * t1
     1              + 0.125 * u ( ixx, iyy0, izz0, :, 34) * t1

          U_Lf_Dn_Pt ( ixx, :) = u ( ixx, iyy0, izz0, :, 1)
     1                   + 0.5 * u ( ixx, iyy0, izz0, :, 3)
     1                   + 0.5 * u ( ixx, iyy0, izz0, :, 4)
     1                  - ob12 * u ( ixx, iyy0, izz0, :, 5)
     1                   + ob6 * u ( ixx, iyy0, izz0, :, 6)
     1                   + ob6 * u ( ixx, iyy0, izz0, :, 7)
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 9)
     1                  + ob20 * u ( ixx, iyy0, izz0, :, 12)
     1                  + ob20 * u ( ixx, iyy0, izz0, :, 13)
     1                  - ob24 * u ( ixx, iyy0, izz0, :, 14)
     1                  - ob24 * u ( ixx, iyy0, izz0, :, 15)
     1                  + ob12 * u ( ixx, iyy0, izz0, :, 17)
     1                  + ob12 * u ( ixx, iyy0, izz0, :, 19)
     1                   + 0.5 * u ( ixx, iyy0, izz0, :, 21) * t1
     1                   + ob3 * u ( ixx, iyy0, izz0, :, 22) * t2
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 23) * t3
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 25) * t1
     1                  + 0.25 * u ( ixx, iyy0, izz0, :, 26) * t1
     1                   + ob6 * u ( ixx, iyy0, izz0, :, 28) * t2
     1                   + ob6 * u ( ixx, iyy0, izz0, :, 29) * t2
     1                  - ob24 * u ( ixx, iyy0, izz0, :, 30) * t1
     1                  + ob12 * u ( ixx, iyy0, izz0, :, 31) * t1
     1                  + ob12 * u ( ixx, iyy0, izz0, :, 32) * t1
     1                 + 0.125 * u ( ixx, iyy0, izz0, :, 34) * t1

! ----------

          U_Rt_Dn ( ixx, :) = u ( ixx, iyyp1, izz0, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izz0, :, 3)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 4)
     1                + ob6 * u ( ixx, iyyp1, izz0, :, 6)
     1                + ob6 * u ( ixx, iyyp1, izz0, :, 7)
     1               - 0.25 * u ( ixx, iyyp1, izz0, :, 9)
     1               - ob20 * u ( ixx, iyyp1, izz0, :, 12)
     1               + ob20 * u ( ixx, iyyp1, izz0, :, 13)
     1               + ob12 * u ( ixx, iyyp1, izz0, :, 17)
     1               - ob12 * u ( ixx, iyyp1, izz0, :, 19)
     1                + 0.5 * u ( ixx, iyyp1, izz0, :, 21) * t1
     1                + ob3 * u ( ixx, iyyp1, izz0, :, 22) * t2
     1               + 0.25 * u ( ixx, iyyp1, izz0, :, 23) * t3
     1               - 0.25 * u ( ixx, iyyp1, izz0, :, 25) * t1
     1               + 0.25 * u ( ixx, iyyp1, izz0, :, 26) * t1
     1                - ob6 * u ( ixx, iyyp1, izz0, :, 28) * t2
     1                + ob6 * u ( ixx, iyyp1, izz0, :, 29) * t2
     1               + ob12 * u ( ixx, iyyp1, izz0, :, 31) * t1
     1               + ob12 * u ( ixx, iyyp1, izz0, :, 32) * t1
     1              - 0.125 * u ( ixx, iyyp1, izz0, :, 34) * t1

          U_Rt_Dn_Pt ( ixx, :) = u ( ixx, iyyp1, izz0, :, 1)
     1                   - 0.5 * u ( ixx, iyyp1, izz0, :, 3)
     1                   + 0.5 * u ( ixx, iyyp1, izz0, :, 4)
     1                  - ob12 * u ( ixx, iyyp1, izz0, :, 5)
     1                   + ob6 * u ( ixx, iyyp1, izz0, :, 6)
     1                   + ob6 * u ( ixx, iyyp1, izz0, :, 7)
     1                  - 0.25 * u ( ixx, iyyp1, izz0, :, 9)
     1                  - ob20 * u ( ixx, iyyp1, izz0, :, 12)
     1                  + ob20 * u ( ixx, iyyp1, izz0, :, 13)
     1                  + ob24 * u ( ixx, iyyp1, izz0, :, 14)
     1                  - ob24 * u ( ixx, iyyp1, izz0, :, 15)
     1                  + ob12 * u ( ixx, iyyp1, izz0, :, 17)
     1                  - ob12 * u ( ixx, iyyp1, izz0, :, 19)
     1                   + 0.5 * u ( ixx, iyyp1, izz0, :, 21) * t1
     1                   + ob3 * u ( ixx, iyyp1, izz0, :, 22) * t2
     1                  + 0.25 * u ( ixx, iyyp1, izz0, :, 23) * t3
     1                  - 0.25 * u ( ixx, iyyp1, izz0, :, 25) * t1
     1                  + 0.25 * u ( ixx, iyyp1, izz0, :, 26) * t1
     1                   - ob6 * u ( ixx, iyyp1, izz0, :, 28) * t2
     1                   + ob6 * u ( ixx, iyyp1, izz0, :, 29) * t2
     1                  - ob24 * u ( ixx, iyyp1, izz0, :, 30) * t1
     1                  + ob12 * u ( ixx, iyyp1, izz0, :, 31) * t1
     1                  + ob12 * u ( ixx, iyyp1, izz0, :, 32) * t1
     1                 - 0.125 * u ( ixx, iyyp1, izz0, :, 34) * t1

! ----------

          U_Lf_Up ( ixx, :) = u ( ixx, iyy0, izzp1, :, 1)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyy0, izzp1, :, 4)
     1                + ob6 * u ( ixx, iyy0, izzp1, :, 6)
     1                + ob6 * u ( ixx, iyy0, izzp1, :, 7)
     1               - 0.25 * u ( ixx, iyy0, izzp1, :, 9)
     1               + ob20 * u ( ixx, iyy0, izzp1, :, 12)
     1               - ob20 * u ( ixx, iyy0, izzp1, :, 13)
     1               - ob12 * u ( ixx, iyy0, izzp1, :, 17)
     1               + ob12 * u ( ixx, iyy0, izzp1, :, 19)
     1                + 0.5 * u ( ixx, iyy0, izzp1, :, 21) * t1
     1                + ob3 * u ( ixx, iyy0, izzp1, :, 22) * t2
     1               + 0.25 * u ( ixx, iyy0, izzp1, :, 23) * t3
     1               + 0.25 * u ( ixx, iyy0, izzp1, :, 25) * t1
     1               - 0.25 * u ( ixx, iyy0, izzp1, :, 26) * t1
     1                + ob6 * u ( ixx, iyy0, izzp1, :, 28) * t2
     1                - ob6 * u ( ixx, iyy0, izzp1, :, 29) * t2
     1               + ob12 * u ( ixx, iyy0, izzp1, :, 31) * t1
     1               + ob12 * u ( ixx, iyy0, izzp1, :, 32) * t1
     1              - 0.125 * u ( ixx, iyy0, izzp1, :, 34) * t1

          U_Lf_Up_Pt ( ixx, :) = u ( ixx, iyy0, izzp1, :, 1)
     1                   + 0.5 * u ( ixx, iyy0, izzp1, :, 3)
     1                   - 0.5 * u ( ixx, iyy0, izzp1, :, 4)
     1                  - ob12 * u ( ixx, iyy0, izzp1, :, 5)
     1                   + ob6 * u ( ixx, iyy0, izzp1, :, 6)
     1                   + ob6 * u ( ixx, iyy0, izzp1, :, 7)
     1                  - 0.25 * u ( ixx, iyy0, izzp1, :, 9)
     1                  + ob20 * u ( ixx, iyy0, izzp1, :, 12)
     1                  - ob20 * u ( ixx, iyy0, izzp1, :, 13)
     1                  - ob24 * u ( ixx, iyy0, izzp1, :, 14)
     1                  + ob24 * u ( ixx, iyy0, izzp1, :, 15)
     1                  - ob12 * u ( ixx, iyy0, izzp1, :, 17)
     1                  + ob12 * u ( ixx, iyy0, izzp1, :, 19)
     1                   + 0.5 * u ( ixx, iyy0, izzp1, :, 21) * t1
     1                   + ob3 * u ( ixx, iyy0, izzp1, :, 22) * t2
     1                  + 0.25 * u ( ixx, iyy0, izzp1, :, 23) * t3
     1                  + 0.25 * u ( ixx, iyy0, izzp1, :, 25) * t1
     1                  - 0.25 * u ( ixx, iyy0, izzp1, :, 26) * t1
     1                   + ob6 * u ( ixx, iyy0, izzp1, :, 28) * t2
     1                   - ob6 * u ( ixx, iyy0, izzp1, :, 29) * t2
     1                  - ob24 * u ( ixx, iyy0, izzp1, :, 30) * t1
     1                  + ob12 * u ( ixx, iyy0, izzp1, :, 31) * t1
     1                  + ob12 * u ( ixx, iyy0, izzp1, :, 32) * t1
     1                 - 0.125 * u ( ixx, iyy0, izzp1, :, 34) * t1

! ----------

          U_Rt_Up ( ixx, :) = u ( ixx, iyyp1, izzp1, :, 1)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 3)
     1                - 0.5 * u ( ixx, iyyp1, izzp1, :, 4)
     1                + ob6 * u ( ixx, iyyp1, izzp1, :, 6)
     1                + ob6 * u ( ixx, iyyp1, izzp1, :, 7)
     1               + 0.25 * u ( ixx, iyyp1, izzp1, :, 9)
     1               - ob20 * u ( ixx, iyyp1, izzp1, :, 12)
     1               - ob20 * u ( ixx, iyyp1, izzp1, :, 13)
     1               - ob12 * u ( ixx, iyyp1, izzp1, :, 17)
     1               - ob12 * u ( ixx, iyyp1, izzp1, :, 19)
     1                + 0.5 * u ( ixx, iyyp1, izzp1, :, 21) * t1
     1                + ob3 * u ( ixx, iyyp1, izzp1, :, 22) * t2
     1               + 0.25 * u ( ixx, iyyp1, izzp1, :, 23) * t3
     1               - 0.25 * u ( ixx, iyyp1, izzp1, :, 25) * t1
     1               - 0.25 * u ( ixx, iyyp1, izzp1, :, 26) * t1
     1                - ob6 * u ( ixx, iyyp1, izzp1, :, 28) * t2
     1                - ob6 * u ( ixx, iyyp1, izzp1, :, 29) * t2
     1               + ob12 * u ( ixx, iyyp1, izzp1, :, 31) * t1
     1               + ob12 * u ( ixx, iyyp1, izzp1, :, 32) * t1
     1              + 0.125 * u ( ixx, iyyp1, izzp1, :, 34) * t1

          U_Rt_Up_Pt ( ixx, :) = u ( ixx, iyyp1, izzp1, :, 1)
     1                   - 0.5 * u ( ixx, iyyp1, izzp1, :, 3)
     1                   - 0.5 * u ( ixx, iyyp1, izzp1, :, 4)
     1                  - ob12 * u ( ixx, iyyp1, izzp1, :, 5)
     1                   + ob6 * u ( ixx, iyyp1, izzp1, :, 6)
     1                   + ob6 * u ( ixx, iyyp1, izzp1, :, 7)
     1                  + 0.25 * u ( ixx, iyyp1, izzp1, :, 9)
     1                  - ob20 * u ( ixx, iyyp1, izzp1, :, 12)
     1                  - ob20 * u ( ixx, iyyp1, izzp1, :, 13)
     1                  + ob24 * u ( ixx, iyyp1, izzp1, :, 14)
     1                  + ob24 * u ( ixx, iyyp1, izzp1, :, 15)
     1                  - ob12 * u ( ixx, iyyp1, izzp1, :, 17)
     1                  - ob12 * u ( ixx, iyyp1, izzp1, :, 19)
     1                   + 0.5 * u ( ixx, iyyp1, izzp1, :, 21) * t1
     1                   + ob3 * u ( ixx, iyyp1, izzp1, :, 22) * t2
     1                  + 0.25 * u ( ixx, iyyp1, izzp1, :, 23) * t3
     1                  - 0.25 * u ( ixx, iyyp1, izzp1, :, 25) * t1
     1                  - 0.25 * u ( ixx, iyyp1, izzp1, :, 26) * t1
     1                   - ob6 * u ( ixx, iyyp1, izzp1, :, 28) * t2
     1                   - ob6 * u ( ixx, iyyp1, izzp1, :, 29) * t2
     1                  - ob24 * u ( ixx, iyyp1, izzp1, :, 30) * t1
     1                  + ob12 * u ( ixx, iyyp1, izzp1, :, 31) * t1
     1                  + ob12 * u ( ixx, iyyp1, izzp1, :, 32) * t1
     1                 + 0.125 * u ( ixx, iyyp1, izzp1, :, 34) * t1

! --------------------
#endif
! ------------------------------
! Build any entropy-related data for the four states.

          entropy_Lf_Dn ( ixx) = entropygr ( ixx, iyy0, izz0, 1)
     1                   + 0.5 * entropygr ( ixx, iyy0, izz0, 3)
     1                   + 0.5 * entropygr ( ixx, iyy0, izz0, 4)
     1              + dt_shift * entropygr ( ixx, iyy0, izz0, 5)

          entropy_Rt_Dn ( ixx) = entropygr ( ixx, iyyp1, izz0, 1)
     1                   - 0.5 * entropygr ( ixx, iyyp1, izz0, 3)
     1                   + 0.5 * entropygr ( ixx, iyyp1, izz0, 4)
     1              + dt_shift * entropygr ( ixx, iyyp1, izz0, 5)

          entropy_Lf_Up ( ixx) = entropygr ( ixx, iyy0, izzp1, 1)
     1                   + 0.5 * entropygr ( ixx, iyy0, izzp1, 3)
     1                   - 0.5 * entropygr ( ixx, iyy0, izzp1, 4)
     1              + dt_shift * entropygr ( ixx, iyy0, izzp1, 5)

          entropy_Rt_Up ( ixx) = entropygr ( ixx, iyyp1, izzp1, 1)
     1                   - 0.5 * entropygr ( ixx, iyyp1, izzp1, 3)
     1                   - 0.5 * entropygr ( ixx, iyyp1, izzp1, 4)
     1              + dt_shift * entropygr ( ixx, iyyp1, izzp1, 5)

! ------------------------------
! Build any EOS-related data for the four states.

          gamma_eos_Lf_Dn ( ixx) = gamma_eos_gr ( ixx, iyy0, izz0)
          gamma_sdsp_Lf_Dn ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyy0, izz0)

          gamma_eos_Rt_Dn ( ixx) = gamma_eos_gr ( ixx, iyyp1, izz0)
          gamma_sdsp_Rt_Dn ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyyp1, izz0)

          gamma_eos_Lf_Up ( ixx) = gamma_eos_gr ( ixx, iyy0, izzp1)
          gamma_sdsp_Lf_Up ( ixx) =
     1                      gamma_soundspeed_gr ( ixx, iyy0, izzp1)

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

              U_Rt_Dn_Pt ( ixx, :) = U_Lf_Dn_Pt ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 2, izz0) == 3) THEN

              U_Rt_Dn ( ixx, :) = U_Lf_Dn ( ixx, :)
              U_Rt_Dn ( ixx, 3) = - U_Rt_Dn ( ixx, 3)
              U_Rt_Dn ( ixx, 6) = - U_Rt_Dn ( ixx, 6)
              U_Rt_Dn ( ixx, 8) = - U_Rt_Dn ( ixx, 8)

              U_Rt_Dn_Pt ( ixx, :) = U_Lf_Dn_Pt ( ixx, :)
              U_Rt_Dn_Pt ( ixx, 3) = - U_Rt_Dn_Pt ( ixx, 3)
              U_Rt_Dn_Pt ( ixx, 6) = - U_Rt_Dn_Pt ( ixx, 6)
              U_Rt_Dn_Pt ( ixx, 8) = - U_Rt_Dn_Pt ( ixx, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixx, 2, izzp1) == 2) THEN

              U_Rt_Up ( ixx, :) = U_Lf_Up ( ixx, :)

              U_Rt_Up_Pt ( ixx, :) = U_Lf_Up_Pt ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 2, izzp1) == 3) THEN

              U_Rt_Up ( ixx, :) = U_Lf_Up ( ixx, :)
              U_Rt_Up ( ixx, 3) = - U_Rt_Up ( ixx, 3)
              U_Rt_Up ( ixx, 6) = - U_Rt_Up ( ixx, 6)
              U_Rt_Up ( ixx, 8) = - U_Rt_Up ( ixx, 8)

              U_Rt_Up_Pt ( ixx, :) = U_Lf_Up_Pt ( ixx, :)
              U_Rt_Up_Pt ( ixx, 3) = - U_Rt_Up_Pt ( ixx, 3)
              U_Rt_Up_Pt ( ixx, 6) = - U_Rt_Up_Pt ( ixx, 6)
              U_Rt_Up_Pt ( ixx, 8) = - U_Rt_Up_Pt ( ixx, 8)

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

              U_Lf_Dn_Pt ( ixx, :) = U_Rt_Dn_Pt ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 1, izz0) == 3) THEN

              U_Lf_Dn ( ixx, :) = U_Rt_Dn ( ixx, :)
              U_Lf_Dn ( ixx, 3) = - U_Lf_Dn ( ixx, 3)
              U_Lf_Dn ( ixx, 6) = - U_Lf_Dn ( ixx, 6)
              U_Lf_Dn ( ixx, 8) = - U_Lf_Dn ( ixx, 8)

              U_Lf_Dn_Pt ( ixx, :) = U_Rt_Dn_Pt ( ixx, :)
              U_Lf_Dn_Pt ( ixx, 3) = - U_Lf_Dn_Pt ( ixx, 3)
              U_Lf_Dn_Pt ( ixx, 6) = - U_Lf_Dn_Pt ( ixx, 6)
              U_Lf_Dn_Pt ( ixx, 8) = - U_Lf_Dn_Pt ( ixx, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixx, 1, izzp1) == 2) THEN

              U_Lf_Up ( ixx, :) = U_Rt_Up ( ixx, :)

              U_Lf_Up_Pt ( ixx, :) = U_Rt_Up_Pt ( ixx, :)

            ELSE IF ( bcarr_y ( ixx, 1, izzp1) == 3) THEN

              U_Lf_Up ( ixx, :) = U_Rt_Up ( ixx, :)
              U_Lf_Up ( ixx, 3) = - U_Lf_Up ( ixx, 3)
              U_Lf_Up ( ixx, 6) = - U_Lf_Up ( ixx, 6)
              U_Lf_Up ( ixx, 8) = - U_Lf_Up ( ixx, 8)

              U_Lf_Up_Pt ( ixx, :) = U_Rt_Up_Pt ( ixx, :)
              U_Lf_Up_Pt ( ixx, 3) = - U_Lf_Up_Pt ( ixx, 3)
              U_Lf_Up_Pt ( ixx, 6) = - U_Lf_Up_Pt ( ixx, 6)
              U_Lf_Up_Pt ( ixx, 8) = - U_Lf_Up_Pt ( ixx, 8)

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

              U_Lf_Up_Pt ( ixx, :) = U_Lf_Dn_Pt ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyy0, 2) == 3) THEN

              U_Lf_Up ( ixx, :) = U_Lf_Dn ( ixx, :)
              U_Lf_Up ( ixx, 4) = - U_Lf_Up ( ixx, 4)
              U_Lf_Up ( ixx, 6) = - U_Lf_Up ( ixx, 6)
              U_Lf_Up ( ixx, 7) = - U_Lf_Up ( ixx, 7)

              U_Lf_Up_Pt ( ixx, :) = U_Lf_Dn_Pt ( ixx, :)
              U_Lf_Up_Pt ( ixx, 4) = - U_Lf_Up_Pt ( ixx, 4)
              U_Lf_Up_Pt ( ixx, 6) = - U_Lf_Up_Pt ( ixx, 6)
              U_Lf_Up_Pt ( ixx, 7) = - U_Lf_Up_Pt ( ixx, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixx, iyyp1, 2) == 2) THEN

              U_Rt_Up ( ixx, :) = U_Rt_Dn ( ixx, :)

              U_Rt_Up_Pt ( ixx, :) = U_Rt_Dn_Pt ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyyp1, 2) == 3) THEN

              U_Rt_Up ( ixx, :) = U_Rt_Dn ( ixx, :)
              U_Rt_Up ( ixx, 4) = - U_Rt_Up ( ixx, 4)
              U_Rt_Up ( ixx, 6) = - U_Rt_Up ( ixx, 6)
              U_Rt_Up ( ixx, 7) = - U_Rt_Up ( ixx, 7)

              U_Rt_Up_Pt ( ixx, :) = U_Rt_Dn_Pt ( ixx, :)
              U_Rt_Up_Pt ( ixx, 4) = - U_Rt_Up_Pt ( ixx, 4)
              U_Rt_Up_Pt ( ixx, 6) = - U_Rt_Up_Pt ( ixx, 6)
              U_Rt_Up_Pt ( ixx, 7) = - U_Rt_Up_Pt ( ixx, 7)

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

              U_Lf_Dn_Pt ( ixx, :) = U_Lf_Up_Pt ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyy0, 1) == 3) THEN

              U_Lf_Dn ( ixx, :) = U_Lf_Up ( ixx, :)
              U_Lf_Dn ( ixx, 4) = - U_Lf_Dn ( ixx, 4)
              U_Lf_Dn ( ixx, 6) = - U_Lf_Dn ( ixx, 6)
              U_Lf_Dn ( ixx, 7) = - U_Lf_Dn ( ixx, 7)

              U_Lf_Dn_Pt ( ixx, :) = U_Lf_Up_Pt ( ixx, :)
              U_Lf_Dn_Pt ( ixx, 4) = - U_Lf_Dn_Pt ( ixx, 4)
              U_Lf_Dn_Pt ( ixx, 6) = - U_Lf_Dn_Pt ( ixx, 6)
              U_Lf_Dn_Pt ( ixx, 7) = - U_Lf_Dn_Pt ( ixx, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixx, iyyp1, 1) == 2) THEN

              U_Rt_Dn ( ixx, :) = U_Rt_Up ( ixx, :)

              U_Rt_Dn_Pt ( ixx, :) = U_Rt_Up_Pt ( ixx, :)

            ELSE IF ( bcarr_z ( ixx, iyyp1, 1) == 3) THEN

              U_Rt_Dn ( ixx, :) = U_Rt_Up ( ixx, :)
              U_Rt_Dn ( ixx, 4) = - U_Rt_Dn ( ixx, 4)
              U_Rt_Dn ( ixx, 6) = - U_Rt_Dn ( ixx, 6)
              U_Rt_Dn ( ixx, 7) = - U_Rt_Dn ( ixx, 7)

              U_Rt_Dn_Pt ( ixx, :) = U_Rt_Up_Pt ( ixx, :)
              U_Rt_Dn_Pt ( ixx, 4) = - U_Rt_Dn_Pt ( ixx, 4)
              U_Rt_Dn_Pt ( ixx, 6) = - U_Rt_Dn_Pt ( ixx, 6)
              U_Rt_Dn_Pt ( ixx, 7) = - U_Rt_Dn_Pt ( ixx, 7)

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
#if ( USE_ONLY_1D_RS != 1)
! ------------------------------

          my_blend_HLL = blend_HLL ( ixx)
          my_blend_MDHLL = blend_mdhll_to_mdhllc

          U_vrtx ( 1, :) = U_Rt_Up_Pt ( ixx, :)
          U_vrtx ( 2, :) = U_Lf_Up_Pt ( ixx, :)
          U_vrtx ( 3, :) = U_Lf_Dn_Pt ( ixx, :)
          U_vrtx ( 4, :) = U_Rt_Dn_Pt ( ixx, :)

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

            flux_temp ( ixx, iyy, izz + 1, :)=
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
     1                                 flux_z ( ixx, iyy + 1, izz, :)
     1                               + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_z ( ixx, iyy + 1, izz, 4) =
     1    flux_z ( ixx, iyy + 1, izz, 4) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
#endif
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
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        flux_x, flux_z, flux_temp, efygr,

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
     1        n_cc_components, n_cc_modes, n_fc_modes,
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

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_z, flux_temp

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) :: efygr

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

      REAL ob_6_multid, ob_6_oned,
     1     t1, t2, t3, ob3, ob6, ob12, ob20, ob24

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + 1 + NFLUID) ::

     1      U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn,

     1      U_Rt_Up_Pt, U_Lf_Up_Pt, U_Lf_Dn_Pt, U_Rt_Dn_Pt

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

! In true 2d, only the y-edges contribute to a problem whose variation
! is entirely in the xz-plane.

      IF ( ( TRUE_2D == 1) .AND. ( ioffy /= 0) ) RETURN

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

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

! ----------

! Time coefficients in terms of the fractional time "dt_shift".
      t1 = dt_shift
      t2 = dt_shift**2
      t3 = dt_shift**3

! Some useful constants for evaluating states.
      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob20 = 1.0 / 20.0
      ob24 = 1.0 / 24.0

! ----------------------------------------------------------------------

        multid_wavemodel_num = 1  ! 1 includes supersonic; 2 is fully subsonic.

! "4 * ob_6_multid + 4 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux

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
!$OMP&    U_Rt_Up_Pt, U_Lf_Up_Pt, U_Lf_Dn_Pt, U_Rt_Dn_Pt,               &
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
!$OMP&    t1, t2, t3, ob3, ob6, ob12, ob20, ob24,                       &

!$OMP&    nzone_vrtx, Vlx_ALE, Vly_ALE, nx, ny,                         &
!$OMP&    multid_wavemodel_num, split_prsflux_re,                       &

!$OMP&    x_indx_limits, y_indx_limits, z_indx_limits,                  &
!$OMP&    bcarr_x, bcarr_y, bcarr_z, dx_save, dy_save, dz_save,         &
!$OMP&    u, entropygr, gamma_eos_gr, gamma_soundspeed_gr, flattengr,   &
!$OMP&    efygr, flux_x, flux_z, flux_temp)                             &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 0 - ioffz_expand, iz1 + ioffz_expand

      izzm1 = MAX0 ( izz - 1, 0 - ioffz)
      izz0 = MAX0 ( izz, 1 - ioffz)
      izzp1 = MIN0 ( izz + 1, iz1 + ioffz)

        DO ixx = 0 - ioffx_expand, ix1 + ioffx_expand

        ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)
        ixx0 = MAX0 ( ixx, 1 - ioffx)
        ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)

! ----------------------------------------------------------------------
! Load in the four states that surround the y-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------
! Indices ( ixx0, izz0) always map to an interior point on the mesh including
! its ghost zones. This is needed because ( ixx, izz) map to the edges.

          nriembegin = y_indx_limits ( ixx0, 1, izz0) - ioffy_expand
          nriemend = y_indx_limits ( ixx0, 2, izz0) + ioffy_expand

          DO iyy = nriembegin, nriemend

! ------------------------------
! First, we build the four states at various orders.
! ------------------------------
#if ( NUMORDER == 2)
! --------------------

          U_Lf_Dn ( iyy, :) = u ( ixx0, iyy, izz0, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 4)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 5) * t1

          U_Lf_Dn_Pt ( iyy, :) = U_Lf_Dn ( iyy, :)

! ----------

          U_Lf_Up ( iyy, :) = u ( ixxp1, iyy, izz0, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 4)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 5) * t1

          U_Lf_Up_Pt ( iyy, :) = U_Lf_Up ( iyy, :)

! ----------

          U_Rt_Dn ( iyy, :) = u ( ixx0, iyy, izzp1, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixx0, iyy, izzp1, :, 4)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 5) * t1

          U_Rt_Dn_Pt ( iyy, :) = U_Rt_Dn ( iyy, :)

! ----------

          U_Rt_Up ( iyy, :) = u ( ixxp1, iyy, izzp1, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 4)
     1                + 0.5 * u ( ixxp1, iyy, izzp1, :, 5) * t1

          U_Rt_Up_Pt ( iyy, :) = U_Rt_Up ( iyy, :)

! --------------------
#endif
! ------------------------------
#if ( NUMORDER == 3)
! --------------------

          U_Lf_Dn ( iyy, :) = u ( ixx0, iyy, izz0, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 4)
     1                + ob6 * u ( ixx0, iyy, izz0, :, 5)
     1                + ob6 * u ( ixx0, iyy, izz0, :, 7)
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 10)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 11) * t1
     1                + ob3 * u ( ixx0, iyy, izz0, :, 12) * t2
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 13) * t1
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 15) * t1

          U_Lf_Dn_Pt ( iyy, :) = u ( ixx0, iyy, izz0, :, 1)
     1                   + 0.5 * u ( ixx0, iyy, izz0, :, 2)
     1                   + 0.5 * u ( ixx0, iyy, izz0, :, 4)
     1                   + ob6 * u ( ixx0, iyy, izz0, :, 5)
     1                  - ob12 * u ( ixx0, iyy, izz0, :, 6)
     1                   + ob6 * u ( ixx0, iyy, izz0, :, 7)
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 10)
     1                   + 0.5 * u ( ixx0, iyy, izz0, :, 11) * t1
     1                   + ob3 * u ( ixx0, iyy, izz0, :, 12) * t2
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 13) * t1
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 15) * t1

! ----------

          U_Lf_Up ( iyy, :) = u ( ixxp1, iyy, izz0, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 4)
     1                + ob6 * u ( ixxp1, iyy, izz0, :, 5)
     1                + ob6 * u ( ixxp1, iyy, izz0, :, 7)
     1               - 0.25 * u ( ixxp1, iyy, izz0, :, 10)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 11) * t1
     1                + ob3 * u ( ixxp1, iyy, izz0, :, 12) * t2
     1               - 0.25 * u ( ixxp1, iyy, izz0, :, 13) * t1
     1               + 0.25 * u ( ixxp1, iyy, izz0, :, 15) * t1

          U_Lf_Up_Pt ( iyy, :) = u ( ixxp1, iyy, izz0, :, 1)
     1                   - 0.5 * u ( ixxp1, iyy, izz0, :, 2)
     1                   + 0.5 * u ( ixxp1, iyy, izz0, :, 4)
     1                   + ob6 * u ( ixxp1, iyy, izz0, :, 5)
     1                  - ob12 * u ( ixxp1, iyy, izz0, :, 6)
     1                   + ob6 * u ( ixxp1, iyy, izz0, :, 7)
     1                  - 0.25 * u ( ixxp1, iyy, izz0, :, 10)
     1                   + 0.5 * u ( ixxp1, iyy, izz0, :, 11) * t1
     1                   + ob3 * u ( ixxp1, iyy, izz0, :, 12) * t2
     1                  - 0.25 * u ( ixxp1, iyy, izz0, :, 13) * t1
     1                  + 0.25 * u ( ixxp1, iyy, izz0, :, 15) * t1

! ----------

          U_Rt_Dn ( iyy, :) = u ( ixx0, iyy, izzp1, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixx0, iyy, izzp1, :, 4)
     1                + ob6 * u ( ixx0, iyy, izzp1, :, 5)
     1                + ob6 * u ( ixx0, iyy, izzp1, :, 7)
     1               - 0.25 * u ( ixx0, iyy, izzp1, :, 10)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 11) * t1
     1                + ob3 * u ( ixx0, iyy, izzp1, :, 12) * t2
     1               + 0.25 * u ( ixx0, iyy, izzp1, :, 13) * t1
     1               - 0.25 * u ( ixx0, iyy, izzp1, :, 15) * t1

          U_Rt_Dn_Pt ( iyy, :) = u ( ixx0, iyy, izzp1, :, 1)
     1                   + 0.5 * u ( ixx0, iyy, izzp1, :, 2)
     1                   - 0.5 * u ( ixx0, iyy, izzp1, :, 4)
     1                   + ob6 * u ( ixx0, iyy, izzp1, :, 5)
     1                  - ob12 * u ( ixx0, iyy, izzp1, :, 6)
     1                   + ob6 * u ( ixx0, iyy, izzp1, :, 7)
     1                  - 0.25 * u ( ixx0, iyy, izzp1, :, 10)
     1                   + 0.5 * u ( ixx0, iyy, izzp1, :, 11) * t1
     1                   + ob3 * u ( ixx0, iyy, izzp1, :, 12) * t2
     1                  + 0.25 * u ( ixx0, iyy, izzp1, :, 13) * t1
     1                  - 0.25 * u ( ixx0, iyy, izzp1, :, 15) * t1

! ----------

          U_Rt_Up ( iyy, :) = u ( ixxp1, iyy, izzp1, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 4)
     1                + ob6 * u ( ixxp1, iyy, izzp1, :, 5)
     1                + ob6 * u ( ixxp1, iyy, izzp1, :, 7)
     1               + 0.25 * u ( ixxp1, iyy, izzp1, :, 10)
     1                + 0.5 * u ( ixxp1, iyy, izzp1, :, 11) * t1
     1                + ob3 * u ( ixxp1, iyy, izzp1, :, 12) * t2
     1               - 0.25 * u ( ixxp1, iyy, izzp1, :, 13) * t1
     1               - 0.25 * u ( ixxp1, iyy, izzp1, :, 15) * t1

          U_Rt_Up_Pt ( iyy, :) = u ( ixxp1, iyy, izzp1, :, 1)
     1                   - 0.5 * u ( ixxp1, iyy, izzp1, :, 2)
     1                   - 0.5 * u ( ixxp1, iyy, izzp1, :, 4)
     1                   + ob6 * u ( ixxp1, iyy, izzp1, :, 5)
     1                  - ob12 * u ( ixxp1, iyy, izzp1, :, 6)
     1                   + ob6 * u ( ixxp1, iyy, izzp1, :, 7)
     1                  + 0.25 * u ( ixxp1, iyy, izzp1, :, 10)
     1                   + 0.5 * u ( ixxp1, iyy, izzp1, :, 11) * t1
     1                   + ob3 * u ( ixxp1, iyy, izzp1, :, 12) * t2
     1                  - 0.25 * u ( ixxp1, iyy, izzp1, :, 13) * t1
     1                  - 0.25 * u ( ixxp1, iyy, izzp1, :, 15) * t1

! --------------------
#endif
! ------------------------------
#if ( NUMORDER == 4)
! --------------------

          U_Lf_Dn ( iyy, :) = u ( ixx0, iyy, izz0, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 4)
     1                + ob6 * u ( ixx0, iyy, izz0, :, 5)
     1                + ob6 * u ( ixx0, iyy, izz0, :, 7)
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 10)
     1               + ob20 * u ( ixx0, iyy, izz0, :, 11)
     1               + ob20 * u ( ixx0, iyy, izz0, :, 13)
     1               + ob12 * u ( ixx0, iyy, izz0, :, 15)
     1               + ob12 * u ( ixx0, iyy, izz0, :, 18)
     1                + 0.5 * u ( ixx0, iyy, izz0, :, 21) * t1
     1                + ob3 * u ( ixx0, iyy, izz0, :, 22) * t2
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 23) * t3
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 24) * t1
     1               + 0.25 * u ( ixx0, iyy, izz0, :, 26) * t1
     1                + ob6 * u ( ixx0, iyy, izz0, :, 27) * t2
     1                + ob6 * u ( ixx0, iyy, izz0, :, 29) * t2
     1               + ob12 * u ( ixx0, iyy, izz0, :, 30) * t1
     1               + ob12 * u ( ixx0, iyy, izz0, :, 32) * t1
     1              + 0.125 * u ( ixx0, iyy, izz0, :, 35) * t1

          U_Lf_Dn_Pt ( iyy, :) = u ( ixx0, iyy, izz0, :, 1)
     1                   + 0.5 * u ( ixx0, iyy, izz0, :, 2)
     1                   + 0.5 * u ( ixx0, iyy, izz0, :, 4)
     1                   + ob6 * u ( ixx0, iyy, izz0, :, 5)
     1                  - ob12 * u ( ixx0, iyy, izz0, :, 6)
     1                   + ob6 * u ( ixx0, iyy, izz0, :, 7)
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 10)
     1                  + ob20 * u ( ixx0, iyy, izz0, :, 11)
     1                  + ob20 * u ( ixx0, iyy, izz0, :, 13)
     1                  + ob12 * u ( ixx0, iyy, izz0, :, 15)
     1                  - ob24 * u ( ixx0, iyy, izz0, :, 16)
     1                  - ob24 * u ( ixx0, iyy, izz0, :, 17)
     1                  + ob12 * u ( ixx0, iyy, izz0, :, 18)
     1                   + 0.5 * u ( ixx0, iyy, izz0, :, 21) * t1
     1                   + ob3 * u ( ixx0, iyy, izz0, :, 22) * t2
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 23) * t3
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 24) * t1
     1                  + 0.25 * u ( ixx0, iyy, izz0, :, 26) * t1
     1                   + ob6 * u ( ixx0, iyy, izz0, :, 27) * t2
     1                   + ob6 * u ( ixx0, iyy, izz0, :, 29) * t2
     1                  + ob12 * u ( ixx0, iyy, izz0, :, 30) * t1
     1                  - ob24 * u ( ixx0, iyy, izz0, :, 31) * t1
     1                  + ob12 * u ( ixx0, iyy, izz0, :, 32) * t1
     1                 + 0.125 * u ( ixx0, iyy, izz0, :, 35) * t1

! ----------

          U_Lf_Up ( iyy, :) = u ( ixxp1, iyy, izz0, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izz0, :, 2)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 4)
     1                + ob6 * u ( ixxp1, iyy, izz0, :, 5)
     1                + ob6 * u ( ixxp1, iyy, izz0, :, 7)
     1               - 0.25 * u ( ixxp1, iyy, izz0, :, 10)
     1               - ob20 * u ( ixxp1, iyy, izz0, :, 11)
     1               + ob20 * u ( ixxp1, iyy, izz0, :, 13)
     1               + ob12 * u ( ixxp1, iyy, izz0, :, 15)
     1               - ob12 * u ( ixxp1, iyy, izz0, :, 18)
     1                + 0.5 * u ( ixxp1, iyy, izz0, :, 21) * t1
     1                + ob3 * u ( ixxp1, iyy, izz0, :, 22) * t2
     1               + 0.25 * u ( ixxp1, iyy, izz0, :, 23) * t3
     1               - 0.25 * u ( ixxp1, iyy, izz0, :, 24) * t1
     1               + 0.25 * u ( ixxp1, iyy, izz0, :, 26) * t1
     1                - ob6 * u ( ixxp1, iyy, izz0, :, 27) * t2
     1                + ob6 * u ( ixxp1, iyy, izz0, :, 29) * t2
     1               + ob12 * u ( ixxp1, iyy, izz0, :, 30) * t1
     1               + ob12 * u ( ixxp1, iyy, izz0, :, 32) * t1
     1              - 0.125 * u ( ixxp1, iyy, izz0, :, 35) * t1

          U_Lf_Up_Pt ( iyy, :) = u ( ixxp1, iyy, izz0, :, 1)
     1                   - 0.5 * u ( ixxp1, iyy, izz0, :, 2)
     1                   + 0.5 * u ( ixxp1, iyy, izz0, :, 4)
     1                   + ob6 * u ( ixxp1, iyy, izz0, :, 5)
     1                  - ob12 * u ( ixxp1, iyy, izz0, :, 6)
     1                   + ob6 * u ( ixxp1, iyy, izz0, :, 7)
     1                  - 0.25 * u ( ixxp1, iyy, izz0, :, 10)
     1                  - ob20 * u ( ixxp1, iyy, izz0, :, 11)
     1                  + ob20 * u ( ixxp1, iyy, izz0, :, 13)
     1                  + ob12 * u ( ixxp1, iyy, izz0, :, 15)
     1                  + ob24 * u ( ixxp1, iyy, izz0, :, 16)
     1                  - ob24 * u ( ixxp1, iyy, izz0, :, 17)
     1                  - ob12 * u ( ixxp1, iyy, izz0, :, 18)
     1                   + 0.5 * u ( ixxp1, iyy, izz0, :, 21) * t1
     1                   + ob3 * u ( ixxp1, iyy, izz0, :, 22) * t2
     1                  + 0.25 * u ( ixxp1, iyy, izz0, :, 23) * t3
     1                  - 0.25 * u ( ixxp1, iyy, izz0, :, 24) * t1
     1                  + 0.25 * u ( ixxp1, iyy, izz0, :, 26) * t1
     1                   - ob6 * u ( ixxp1, iyy, izz0, :, 27) * t2
     1                   + ob6 * u ( ixxp1, iyy, izz0, :, 29) * t2
     1                  + ob12 * u ( ixxp1, iyy, izz0, :, 30) * t1
     1                  - ob24 * u ( ixxp1, iyy, izz0, :, 31) * t1
     1                  + ob12 * u ( ixxp1, iyy, izz0, :, 32) * t1
     1                 - 0.125 * u ( ixxp1, iyy, izz0, :, 35) * t1

! ----------

          U_Rt_Dn ( iyy, :) = u ( ixx0, iyy, izzp1, :, 1)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixx0, iyy, izzp1, :, 4)
     1                + ob6 * u ( ixx0, iyy, izzp1, :, 5)
     1                + ob6 * u ( ixx0, iyy, izzp1, :, 7)
     1               - 0.25 * u ( ixx0, iyy, izzp1, :, 10)
     1               + ob20 * u ( ixx0, iyy, izzp1, :, 11)
     1               - ob20 * u ( ixx0, iyy, izzp1, :, 13)
     1               - ob12 * u ( ixx0, iyy, izzp1, :, 15)
     1               + ob12 * u ( ixx0, iyy, izzp1, :, 18)
     1                + 0.5 * u ( ixx0, iyy, izzp1, :, 21) * t1
     1                + ob3 * u ( ixx0, iyy, izzp1, :, 22) * t2
     1               + 0.25 * u ( ixx0, iyy, izzp1, :, 23) * t3
     1               + 0.25 * u ( ixx0, iyy, izzp1, :, 24) * t1
     1               - 0.25 * u ( ixx0, iyy, izzp1, :, 26) * t1
     1                + ob6 * u ( ixx0, iyy, izzp1, :, 27) * t2
     1                - ob6 * u ( ixx0, iyy, izzp1, :, 29) * t2
     1               + ob12 * u ( ixx0, iyy, izzp1, :, 30) * t1
     1               + ob12 * u ( ixx0, iyy, izzp1, :, 32) * t1
     1              - 0.125 * u ( ixx0, iyy, izzp1, :, 35) * t1

          U_Rt_Dn_Pt ( iyy, :) = u ( ixx0, iyy, izzp1, :, 1)
     1                   + 0.5 * u ( ixx0, iyy, izzp1, :, 2)
     1                   - 0.5 * u ( ixx0, iyy, izzp1, :, 4)
     1                   + ob6 * u ( ixx0, iyy, izzp1, :, 5)
     1                  - ob12 * u ( ixx0, iyy, izzp1, :, 6)
     1                   + ob6 * u ( ixx0, iyy, izzp1, :, 7)
     1                  - 0.25 * u ( ixx0, iyy, izzp1, :, 10)
     1                  + ob20 * u ( ixx0, iyy, izzp1, :, 11)
     1                  - ob20 * u ( ixx0, iyy, izzp1, :, 13)
     1                  - ob12 * u ( ixx0, iyy, izzp1, :, 15)
     1                  - ob24 * u ( ixx0, iyy, izzp1, :, 16)
     1                  + ob24 * u ( ixx0, iyy, izzp1, :, 17)
     1                  + ob12 * u ( ixx0, iyy, izzp1, :, 18)
     1                   + 0.5 * u ( ixx0, iyy, izzp1, :, 21) * t1
     1                   + ob3 * u ( ixx0, iyy, izzp1, :, 22) * t2
     1                  + 0.25 * u ( ixx0, iyy, izzp1, :, 23) * t3
     1                  + 0.25 * u ( ixx0, iyy, izzp1, :, 24) * t1
     1                  - 0.25 * u ( ixx0, iyy, izzp1, :, 26) * t1
     1                   + ob6 * u ( ixx0, iyy, izzp1, :, 27) * t2
     1                   - ob6 * u ( ixx0, iyy, izzp1, :, 29) * t2
     1                  + ob12 * u ( ixx0, iyy, izzp1, :, 30) * t1
     1                  - ob24 * u ( ixx0, iyy, izzp1, :, 31) * t1
     1                  + ob12 * u ( ixx0, iyy, izzp1, :, 32) * t1
     1                 - 0.125 * u ( ixx0, iyy, izzp1, :, 35) * t1

! ----------

          U_Rt_Up ( iyy, :) = u ( ixxp1, iyy, izzp1, :, 1)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 2)
     1                - 0.5 * u ( ixxp1, iyy, izzp1, :, 4)
     1                + ob6 * u ( ixxp1, iyy, izzp1, :, 5)
     1                + ob6 * u ( ixxp1, iyy, izzp1, :, 7)
     1               + 0.25 * u ( ixxp1, iyy, izzp1, :, 10)
     1               - ob20 * u ( ixxp1, iyy, izzp1, :, 11)
     1               - ob20 * u ( ixxp1, iyy, izzp1, :, 13)
     1               - ob12 * u ( ixxp1, iyy, izzp1, :, 15)
     1               - ob12 * u ( ixxp1, iyy, izzp1, :, 18)
     1                + 0.5 * u ( ixxp1, iyy, izzp1, :, 21) * t1
     1                + ob3 * u ( ixxp1, iyy, izzp1, :, 22) * t2
     1               + 0.25 * u ( ixxp1, iyy, izzp1, :, 23) * t3
     1               - 0.25 * u ( ixxp1, iyy, izzp1, :, 24) * t1
     1               - 0.25 * u ( ixxp1, iyy, izzp1, :, 26) * t1
     1                - ob6 * u ( ixxp1, iyy, izzp1, :, 27) * t2
     1                - ob6 * u ( ixxp1, iyy, izzp1, :, 29) * t2
     1               + ob12 * u ( ixxp1, iyy, izzp1, :, 30) * t1
     1               + ob12 * u ( ixxp1, iyy, izzp1, :, 32) * t1
     1              + 0.125 * u ( ixxp1, iyy, izzp1, :, 35) * t1

          U_Rt_Up_Pt ( iyy, :) = u ( ixxp1, iyy, izzp1, :, 1)
     1                   - 0.5 * u ( ixxp1, iyy, izzp1, :, 2)
     1                   - 0.5 * u ( ixxp1, iyy, izzp1, :, 4)
     1                   + ob6 * u ( ixxp1, iyy, izzp1, :, 5)
     1                  - ob12 * u ( ixxp1, iyy, izzp1, :, 6)
     1                   + ob6 * u ( ixxp1, iyy, izzp1, :, 7)
     1                  + 0.25 * u ( ixxp1, iyy, izzp1, :, 10)
     1                  - ob20 * u ( ixxp1, iyy, izzp1, :, 11)
     1                  - ob20 * u ( ixxp1, iyy, izzp1, :, 13)
     1                  - ob12 * u ( ixxp1, iyy, izzp1, :, 15)
     1                  + ob24 * u ( ixxp1, iyy, izzp1, :, 16)
     1                  + ob24 * u ( ixxp1, iyy, izzp1, :, 17)
     1                  - ob12 * u ( ixxp1, iyy, izzp1, :, 18)
     1                   + 0.5 * u ( ixxp1, iyy, izzp1, :, 21) * t1
     1                   + ob3 * u ( ixxp1, iyy, izzp1, :, 22) * t2
     1                  + 0.25 * u ( ixxp1, iyy, izzp1, :, 23) * t3
     1                  - 0.25 * u ( ixxp1, iyy, izzp1, :, 24) * t1
     1                  - 0.25 * u ( ixxp1, iyy, izzp1, :, 26) * t1
     1                   - ob6 * u ( ixxp1, iyy, izzp1, :, 27) * t2
     1                   - ob6 * u ( ixxp1, iyy, izzp1, :, 29) * t2
     1                  + ob12 * u ( ixxp1, iyy, izzp1, :, 30) * t1
     1                  - ob24 * u ( ixxp1, iyy, izzp1, :, 31) * t1
     1                  + ob12 * u ( ixxp1, iyy, izzp1, :, 32) * t1
     1                 + 0.125 * u ( ixxp1, iyy, izzp1, :, 35) * t1

! --------------------
#endif
! ------------------------------
! Build any entropy-related data for the four states.

          entropy_Lf_Dn ( iyy) = entropygr ( ixx0, iyy, izz0, 1)
     1                   + 0.5 * entropygr ( ixx0, iyy, izz0, 2)
     1                   + 0.5 * entropygr ( ixx0, iyy, izz0, 4)
     1              + dt_shift * entropygr ( ixx0, iyy, izz0, 5)

          entropy_Lf_Up ( iyy) = entropygr ( ixxp1, iyy, izz0, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyy, izz0, 2)
     1                   + 0.5 * entropygr ( ixxp1, iyy, izz0, 4)
     1              + dt_shift * entropygr ( ixxp1, iyy, izz0, 5)

          entropy_Rt_Dn ( iyy) = entropygr ( ixx0, iyy, izzp1, 1)
     1                   + 0.5 * entropygr ( ixx0, iyy, izzp1, 2)
     1                   - 0.5 * entropygr ( ixx0, iyy, izzp1, 4)
     1              + dt_shift * entropygr ( ixx0, iyy, izzp1, 5)

          entropy_Rt_Up ( iyy) = entropygr ( ixxp1, iyy, izzp1, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyy, izzp1, 2)
     1                   - 0.5 * entropygr ( ixxp1, iyy, izzp1, 4)
     1              + dt_shift * entropygr ( ixxp1, iyy, izzp1, 5)

! ------------------------------
! Build any EOS-related data for the four states.

          gamma_eos_Lf_Dn ( iyy) = gamma_eos_gr ( ixx0, iyy, izz0)
          gamma_sdsp_Lf_Dn ( iyy) =
     1                      gamma_soundspeed_gr ( ixx0, iyy, izz0)

          gamma_eos_Lf_Up ( iyy) = gamma_eos_gr ( ixxp1, iyy, izz0)
          gamma_sdsp_Lf_Up ( iyy) =
     1                      gamma_soundspeed_gr ( ixxp1, iyy, izz0)

          gamma_eos_Rt_Dn ( iyy) = gamma_eos_gr ( ixx0, iyy, izzp1)
          gamma_sdsp_Rt_Dn ( iyy) =
     1                      gamma_soundspeed_gr ( ixx0, iyy, izzp1)

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

              U_Lf_Up_Pt ( iyy, :) = U_Lf_Dn_Pt ( iyy, :)

            ELSE IF ( bcarr_x ( 2, iyy, izz0) == 3) THEN

              U_Lf_Up ( iyy, :) = U_Lf_Dn ( iyy, :)
              U_Lf_Up ( iyy, 2) = - U_Lf_Up ( iyy, 2)
              U_Lf_Up ( iyy, 7) = - U_Lf_Up ( iyy, 7)
              U_Lf_Up ( iyy, 8) = - U_Lf_Up ( iyy, 8)

              U_Lf_Up_Pt ( iyy, :) = U_Lf_Dn_Pt ( iyy, :)
              U_Lf_Up_Pt ( iyy, 2) = - U_Lf_Up_Pt ( iyy, 2)
              U_Lf_Up_Pt ( iyy, 7) = - U_Lf_Up_Pt ( iyy, 7)
              U_Lf_Up_Pt ( iyy, 8) = - U_Lf_Up_Pt ( iyy, 8)

            END IF

! ----------

            IF ( bcarr_x ( 2, iyy, izzp1) == 2) THEN

              U_Rt_Up ( iyy, :) = U_Rt_Dn ( iyy, :)

              U_Rt_Up_Pt ( iyy, :) = U_Rt_Dn_Pt ( iyy, :)

            ELSE IF ( bcarr_x ( 2, iyy, izzp1) == 3) THEN

              U_Rt_Up ( iyy, :) = U_Rt_Dn ( iyy, :)
              U_Rt_Up ( iyy, 2) = - U_Rt_Up ( iyy, 2)
              U_Rt_Up ( iyy, 7) = - U_Rt_Up ( iyy, 7)
              U_Rt_Up ( iyy, 8) = - U_Rt_Up ( iyy, 8)

              U_Rt_Up_Pt ( iyy, :) = U_Rt_Dn_Pt ( iyy, :)
              U_Rt_Up_Pt ( iyy, 2) = - U_Rt_Up_Pt ( iyy, 2)
              U_Rt_Up_Pt ( iyy, 7) = - U_Rt_Up_Pt ( iyy, 7)
              U_Rt_Up_Pt ( iyy, 8) = - U_Rt_Up_Pt ( iyy, 8)

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

              U_Lf_Dn_Pt ( iyy, :) = U_Lf_Up_Pt ( iyy, :)

            ELSE IF ( bcarr_x ( 1, iyy, izz0) == 3) THEN

              U_Lf_Dn ( iyy, :) = U_Lf_Up ( iyy, :)
              U_Lf_Dn ( iyy, 2) = - U_Lf_Dn ( iyy, 2)
              U_Lf_Dn ( iyy, 7) = - U_Lf_Dn ( iyy, 7)
              U_Lf_Dn ( iyy, 8) = - U_Lf_Dn ( iyy, 8)

              U_Lf_Dn_Pt ( iyy, :) = U_Lf_Up_Pt ( iyy, :)
              U_Lf_Dn_Pt ( iyy, 2) = - U_Lf_Dn_Pt ( iyy, 2)
              U_Lf_Dn_Pt ( iyy, 7) = - U_Lf_Dn_Pt ( iyy, 7)
              U_Lf_Dn_Pt ( iyy, 8) = - U_Lf_Dn_Pt ( iyy, 8)

            END IF

! ----------

            IF ( bcarr_x ( 1, iyy, izzp1) == 2) THEN

              U_Rt_Dn ( iyy, :) = U_Rt_Up ( iyy, :)

              U_Rt_Dn_Pt ( iyy, :) = U_Rt_Up_Pt ( iyy, :)

            ELSE IF ( bcarr_x ( 1, iyy, izzp1) == 3) THEN

              U_Rt_Dn ( iyy, :) = U_Rt_Up ( iyy, :)
              U_Rt_Dn ( iyy, 2) = - U_Rt_Dn ( iyy, 2)
              U_Rt_Dn ( iyy, 7) = - U_Rt_Dn ( iyy, 7)
              U_Rt_Dn ( iyy, 8) = - U_Rt_Dn ( iyy, 8)

              U_Rt_Dn_Pt ( iyy, :) = U_Rt_Up_Pt ( iyy, :)
              U_Rt_Dn_Pt ( iyy, 2) = - U_Rt_Dn_Pt ( iyy, 2)
              U_Rt_Dn_Pt ( iyy, 7) = - U_Rt_Dn_Pt ( iyy, 7)
              U_Rt_Dn_Pt ( iyy, 8) = - U_Rt_Dn_Pt ( iyy, 8)

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

              U_Rt_Dn_Pt ( iyy, :) = U_Lf_Dn_Pt ( iyy, :)

            ELSE IF ( bcarr_z ( ixx0, iyy, 2) == 3) THEN

              U_Rt_Dn ( iyy, :) = U_Lf_Dn ( iyy, :)
              U_Rt_Dn ( iyy, 4) = - U_Rt_Dn ( iyy, 4)
              U_Rt_Dn ( iyy, 6) = - U_Rt_Dn ( iyy, 6)
              U_Rt_Dn ( iyy, 7) = - U_Rt_Dn ( iyy, 7)

              U_Rt_Dn_Pt ( iyy, :) = U_Lf_Dn_Pt ( iyy, :)
              U_Rt_Dn_Pt ( iyy, 4) = - U_Rt_Dn_Pt ( iyy, 4)
              U_Rt_Dn_Pt ( iyy, 6) = - U_Rt_Dn_Pt ( iyy, 6)
              U_Rt_Dn_Pt ( iyy, 7) = - U_Rt_Dn_Pt ( iyy, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixxp1, iyy, 2) == 2) THEN

              U_Rt_Up ( iyy, :) = U_Lf_Up ( iyy, :)

              U_Rt_Up_Pt ( iyy, :) = U_Lf_Up_Pt ( iyy, :)

            ELSE IF ( bcarr_z ( ixxp1, iyy, 2) == 3) THEN

              U_Rt_Up ( iyy, :) = U_Lf_Up ( iyy, :)
              U_Rt_Up ( iyy, 4) = - U_Rt_Up ( iyy, 4)
              U_Rt_Up ( iyy, 6) = - U_Rt_Up ( iyy, 6)
              U_Rt_Up ( iyy, 7) = - U_Rt_Up ( iyy, 7)

              U_Rt_Up_Pt ( iyy, :) = U_Lf_Up_Pt ( iyy, :)
              U_Rt_Up_Pt ( iyy, 4) = - U_Rt_Up_Pt ( iyy, 4)
              U_Rt_Up_Pt ( iyy, 6) = - U_Rt_Up_Pt ( iyy, 6)
              U_Rt_Up_Pt ( iyy, 7) = - U_Rt_Up_Pt ( iyy, 7)

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

              U_Lf_Dn_Pt ( iyy, :) = U_Rt_Dn_Pt ( iyy, :)

            ELSE IF ( bcarr_z ( ixx0, iyy, 1) == 3) THEN

              U_Lf_Dn ( iyy, :) = U_Rt_Dn ( iyy, :)
              U_Lf_Dn ( iyy, 4) = - U_Lf_Dn ( iyy, 4)
              U_Lf_Dn ( iyy, 6) = - U_Lf_Dn ( iyy, 6)
              U_Lf_Dn ( iyy, 7) = - U_Lf_Dn ( iyy, 7)

              U_Lf_Dn_Pt ( iyy, :) = U_Rt_Dn_Pt ( iyy, :)
              U_Lf_Dn_Pt ( iyy, 4) = - U_Lf_Dn_Pt ( iyy, 4)
              U_Lf_Dn_Pt ( iyy, 6) = - U_Lf_Dn_Pt ( iyy, 6)
              U_Lf_Dn_Pt ( iyy, 7) = - U_Lf_Dn_Pt ( iyy, 7)

            END IF

! ----------

            IF ( bcarr_z ( ixxp1, iyy, 1) == 2) THEN

              U_Lf_Up ( iyy, :) = U_Rt_Up ( iyy, :)

              U_Lf_Up_Pt ( iyy, :) = U_Rt_Up_Pt ( iyy, :)

            ELSE IF ( bcarr_z ( ixxp1, iyy, 1) == 3) THEN

              U_Lf_Up ( iyy, :) = U_Rt_Up ( iyy, :)
              U_Lf_Up ( iyy, 4) = - U_Lf_Up ( iyy, 4)
              U_Lf_Up ( iyy, 6) = - U_Lf_Up ( iyy, 6)
              U_Lf_Up ( iyy, 7) = - U_Lf_Up ( iyy, 7)

              U_Lf_Up_Pt ( iyy, :) = U_Rt_Up_Pt ( iyy, :)
              U_Lf_Up_Pt ( iyy, 4) = - U_Lf_Up_Pt ( iyy, 4)
              U_Lf_Up_Pt ( iyy, 6) = - U_Lf_Up_Pt ( iyy, 6)
              U_Lf_Up_Pt ( iyy, 7) = - U_Lf_Up_Pt ( iyy, 7)

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
#if ( USE_ONLY_1D_RS != 1)
! ------------------------------

          my_blend_HLL = blend_HLL ( iyy)
          my_blend_MDHLL = blend_mdhll_to_mdhllc

          U_vrtx ( 1, :) = U_Rt_Up_Pt ( iyy, :)
          U_vrtx ( 2, :) = U_Lf_Up_Pt ( iyy, :)
          U_vrtx ( 3, :) = U_Lf_Dn_Pt ( iyy, :)
          U_vrtx ( 4, :) = U_Rt_Dn_Pt ( iyy, :)

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

            flux_temp ( ixx, iyy, izz + 1, :)=
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

            flux_z ( ixx + 1, iyy, izz, :)=
     1                                 flux_z ( ixx + 1, iyy, izz, :)
     1                               + ob_6_multid * F_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_z ( ixx + 1, iyy, izz, 4) =
     1    flux_z ( ixx + 1, iyy, izz, 4) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
#endif
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
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, bfxbdy, bfybdy, bfzbdy, entropygr,
     1        flux_x, flux_y, flux_temp, efzgr,

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
     1        n_cc_components, n_cc_modes, n_fc_modes,
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

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: flux_x, flux_y, flux_temp

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) :: efzgr

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

      REAL ob_6_multid, ob_6_oned,
     1     t1, t2, t3, ob3, ob6, ob12, ob20, ob24

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD + 1 + NFLUID) ::

     1      U_Rt_Up, U_Lf_Up, U_Lf_Dn, U_Rt_Dn,

     1      U_Rt_Up_Pt, U_Lf_Up_Pt, U_Lf_Dn_Pt, U_Rt_Dn_Pt

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

! In true 2d, only the z-edges contribute to a problem whose variation
! is entirely in the xy-plane.

      IF ( ( TRUE_2D == 1) .AND. ( ioffz /= 0) ) RETURN

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

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

! ----------

! Time coefficients in terms of the fractional time "dt_shift".
      t1 = dt_shift
      t2 = dt_shift**2
      t3 = dt_shift**3

! Some useful constants for evaluating states.
      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob20 = 1.0 / 20.0
      ob24 = 1.0 / 24.0

! ----------------------------------------------------------------------

        multid_wavemodel_num = 1  ! 1 includes supersonic; 2 is fully subsonic.

! "4 * ob_6_multid + 4 * ob_6_oned == 1"
        ob_6_multid = blend_multid_flux

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
!$OMP&    U_Rt_Up_Pt, U_Lf_Up_Pt, U_Lf_Dn_Pt, U_Rt_Dn_Pt,               &
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
!$OMP&    t1, t2, t3, ob3, ob6, ob12, ob20, ob24,                       &

!$OMP&    nzone_vrtx, Vlx_ALE, Vly_ALE, nx, ny,                         &
!$OMP&    multid_wavemodel_num, split_prsflux_re,                       &

!$OMP&    x_indx_limits, y_indx_limits, z_indx_limits,                  &
!$OMP&    bcarr_x, bcarr_y, bcarr_z, dx_save, dy_save, dz_save,         &
!$OMP&    u, entropygr, gamma_eos_gr, gamma_soundspeed_gr, flattengr,   &
!$OMP&    efzgr, flux_x, flux_y, flux_temp)                             &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 0 - ioffy_expand, iy1 + ioffy_expand

      iyym1 = MAX0 ( iyy - 1, 0 - ioffy)
      iyy0 = MAX0 ( iyy, 1 - ioffy)
      iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)

        DO ixx = 0 - ioffx_expand, ix1 + ioffx_expand

        ixxm1 = MAX0 ( ixx - 1, 0 - ioffx)
        ixx0 = MAX0 ( ixx, 1 - ioffx)
        ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)

! ----------------------------------------------------------------------
! Load in the four states that surround the z-edge for the
! multidimensional Riemann solver.
! ----------------------------------------------------------------------
! Indices ( ixx0, iyy0) always map to an interior point on the mesh including
! its ghost zones. This is needed because ( ixx, iyy) map to the edges.

          nriembegin = z_indx_limits ( ixx0, iyy0, 1) - ioffz_expand
          nriemend = z_indx_limits ( ixx0, iyy0, 2) + ioffz_expand

          DO izz = nriembegin, nriemend

! ------------------------------
! First, we build the four states at various orders.
! ------------------------------
#if ( NUMORDER == 2)
! --------------------

          U_Lf_Dn ( izz, :) = u ( ixx0, iyy0, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 3)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 5) * t1

          U_Lf_Dn_Pt ( izz, :) = U_Lf_Dn ( izz, :)

! ----------

          U_Rt_Dn ( izz, :) = u ( ixxp1, iyy0, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 3)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 5) * t1

          U_Rt_Dn_Pt ( izz, :) = U_Rt_Dn ( izz, :)

! ----------

          U_Lf_Up ( izz, :) = u ( ixx0, iyyp1, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixx0, iyyp1, izz, :, 3)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 5) * t1

          U_Lf_Up_Pt ( izz, :) = U_Lf_Up ( izz, :)

! ----------

          U_Rt_Up ( izz, :) = u ( ixxp1, iyyp1, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 3)
     1                + 0.5 * u ( ixxp1, iyyp1, izz, :, 5) * t1

          U_Rt_Up_Pt ( izz, :) = U_Rt_Up ( izz, :)

! --------------------
#endif
! ------------------------------
#if ( NUMORDER == 3)
! --------------------

          U_Lf_Dn ( izz, :) = u ( ixx0, iyy0, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 3)
     1                + ob6 * u ( ixx0, iyy0, izz, :, 5)
     1                + ob6 * u ( ixx0, iyy0, izz, :, 6)
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 8)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 11) * t1
     1                + ob3 * u ( ixx0, iyy0, izz, :, 12) * t2
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 13) * t1
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 14) * t1

          U_Lf_Dn_Pt ( izz, :) = u ( ixx0, iyy0, izz, :, 1)
     1                   + 0.5 * u ( ixx0, iyy0, izz, :, 2)
     1                   + 0.5 * u ( ixx0, iyy0, izz, :, 3)
     1                   + ob6 * u ( ixx0, iyy0, izz, :, 5)
     1                   + ob6 * u ( ixx0, iyy0, izz, :, 6)
     1                  - ob12 * u ( ixx0, iyy0, izz, :, 7)
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 8)
     1                   + 0.5 * u ( ixx0, iyy0, izz, :, 11) * t1
     1                   + ob3 * u ( ixx0, iyy0, izz, :, 12) * t2
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 13) * t1
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 14) * t1

! ----------

          U_Rt_Dn ( izz, :) = u ( ixxp1, iyy0, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 3)
     1                + ob6 * u ( ixxp1, iyy0, izz, :, 5)
     1                + ob6 * u ( ixxp1, iyy0, izz, :, 6)
     1               - 0.25 * u ( ixxp1, iyy0, izz, :, 8)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 11) * t1
     1                + ob3 * u ( ixxp1, iyy0, izz, :, 12) * t2
     1               - 0.25 * u ( ixxp1, iyy0, izz, :, 13) * t1
     1               + 0.25 * u ( ixxp1, iyy0, izz, :, 14) * t1

          U_Rt_Dn_Pt ( izz, :) = u ( ixxp1, iyy0, izz, :, 1)
     1                   - 0.5 * u ( ixxp1, iyy0, izz, :, 2)
     1                   + 0.5 * u ( ixxp1, iyy0, izz, :, 3)
     1                   + ob6 * u ( ixxp1, iyy0, izz, :, 5)
     1                   + ob6 * u ( ixxp1, iyy0, izz, :, 6)
     1                  - ob12 * u ( ixxp1, iyy0, izz, :, 7)
     1                  - 0.25 * u ( ixxp1, iyy0, izz, :, 8)
     1                   + 0.5 * u ( ixxp1, iyy0, izz, :, 11) * t1
     1                   + ob3 * u ( ixxp1, iyy0, izz, :, 12) * t2
     1                  - 0.25 * u ( ixxp1, iyy0, izz, :, 13) * t1
     1                  + 0.25 * u ( ixxp1, iyy0, izz, :, 14) * t1

! ----------

          U_Lf_Up ( izz, :) = u ( ixx0, iyyp1, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixx0, iyyp1, izz, :, 3)
     1                + ob6 * u ( ixx0, iyyp1, izz, :, 5)
     1                + ob6 * u ( ixx0, iyyp1, izz, :, 6)
     1               - 0.25 * u ( ixx0, iyyp1, izz, :, 8)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 11) * t1
     1                + ob3 * u ( ixx0, iyyp1, izz, :, 12) * t2
     1               + 0.25 * u ( ixx0, iyyp1, izz, :, 13) * t1
     1               - 0.25 * u ( ixx0, iyyp1, izz, :, 14) * t1

          U_Lf_Up_Pt ( izz, :) = u ( ixx0, iyyp1, izz, :, 1)
     1                   + 0.5 * u ( ixx0, iyyp1, izz, :, 2)
     1                   - 0.5 * u ( ixx0, iyyp1, izz, :, 3)
     1                   + ob6 * u ( ixx0, iyyp1, izz, :, 5)
     1                   + ob6 * u ( ixx0, iyyp1, izz, :, 6)
     1                  - ob12 * u ( ixx0, iyyp1, izz, :, 7)
     1                  - 0.25 * u ( ixx0, iyyp1, izz, :, 8)
     1                   + 0.5 * u ( ixx0, iyyp1, izz, :, 11) * t1
     1                   + ob3 * u ( ixx0, iyyp1, izz, :, 12) * t2
     1                  + 0.25 * u ( ixx0, iyyp1, izz, :, 13) * t1
     1                  - 0.25 * u ( ixx0, iyyp1, izz, :, 14) * t1

! ----------

          U_Rt_Up ( izz, :) = u ( ixxp1, iyyp1, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 3)
     1                + ob6 * u ( ixxp1, iyyp1, izz, :, 5)
     1                + ob6 * u ( ixxp1, iyyp1, izz, :, 6)
     1               + 0.25 * u ( ixxp1, iyyp1, izz, :, 8)
     1                + 0.5 * u ( ixxp1, iyyp1, izz, :, 11) * t1
     1                + ob3 * u ( ixxp1, iyyp1, izz, :, 12) * t2
     1               - 0.25 * u ( ixxp1, iyyp1, izz, :, 13) * t1
     1               - 0.25 * u ( ixxp1, iyyp1, izz, :, 14) * t1

          U_Rt_Up_Pt ( izz, :) = u ( ixxp1, iyyp1, izz, :, 1)
     1                   - 0.5 * u ( ixxp1, iyyp1, izz, :, 2)
     1                   - 0.5 * u ( ixxp1, iyyp1, izz, :, 3)
     1                   + ob6 * u ( ixxp1, iyyp1, izz, :, 5)
     1                   + ob6 * u ( ixxp1, iyyp1, izz, :, 6)
     1                  - ob12 * u ( ixxp1, iyyp1, izz, :, 7)
     1                  + 0.25 * u ( ixxp1, iyyp1, izz, :, 8)
     1                   + 0.5 * u ( ixxp1, iyyp1, izz, :, 11) * t1
     1                   + ob3 * u ( ixxp1, iyyp1, izz, :, 12) * t2
     1                  - 0.25 * u ( ixxp1, iyyp1, izz, :, 13) * t1
     1                  - 0.25 * u ( ixxp1, iyyp1, izz, :, 14) * t1

! --------------------
#endif
! ------------------------------
#if ( NUMORDER == 4)
! --------------------

          U_Lf_Dn ( izz, :) = u ( ixx0, iyy0, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 3)
     1                + ob6 * u ( ixx0, iyy0, izz, :, 5)
     1                + ob6 * u ( ixx0, iyy0, izz, :, 6)
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 8)
     1               + ob20 * u ( ixx0, iyy0, izz, :, 11)
     1               + ob20 * u ( ixx0, iyy0, izz, :, 12)
     1               + ob12 * u ( ixx0, iyy0, izz, :, 14)
     1               + ob12 * u ( ixx0, iyy0, izz, :, 16)
     1                + 0.5 * u ( ixx0, iyy0, izz, :, 21) * t1
     1                + ob3 * u ( ixx0, iyy0, izz, :, 22) * t2
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 23) * t3
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 24) * t1
     1               + 0.25 * u ( ixx0, iyy0, izz, :, 25) * t1
     1                + ob6 * u ( ixx0, iyy0, izz, :, 27) * t2
     1                + ob6 * u ( ixx0, iyy0, izz, :, 28) * t2
     1               + ob12 * u ( ixx0, iyy0, izz, :, 30) * t1
     1               + ob12 * u ( ixx0, iyy0, izz, :, 31) * t1
     1              + 0.125 * u ( ixx0, iyy0, izz, :, 33) * t1

          U_Lf_Dn_Pt ( izz, :) = u ( ixx0, iyy0, izz, :, 1)
     1                   + 0.5 * u ( ixx0, iyy0, izz, :, 2)
     1                   + 0.5 * u ( ixx0, iyy0, izz, :, 3)
     1                   + ob6 * u ( ixx0, iyy0, izz, :, 5)
     1                   + ob6 * u ( ixx0, iyy0, izz, :, 6)
     1                  - ob12 * u ( ixx0, iyy0, izz, :, 7)
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 8)
     1                  + ob20 * u ( ixx0, iyy0, izz, :, 11)
     1                  + ob20 * u ( ixx0, iyy0, izz, :, 12)
     1                  + ob12 * u ( ixx0, iyy0, izz, :, 14)
     1                  + ob12 * u ( ixx0, iyy0, izz, :, 16)
     1                  - ob24 * u ( ixx0, iyy0, izz, :, 18)
     1                  - ob24 * u ( ixx0, iyy0, izz, :, 19)
     1                   + 0.5 * u ( ixx0, iyy0, izz, :, 21) * t1
     1                   + ob3 * u ( ixx0, iyy0, izz, :, 22) * t2
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 23) * t3
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 24) * t1
     1                  + 0.25 * u ( ixx0, iyy0, izz, :, 25) * t1
     1                   + ob6 * u ( ixx0, iyy0, izz, :, 27) * t2
     1                   + ob6 * u ( ixx0, iyy0, izz, :, 28) * t2
     1                  + ob12 * u ( ixx0, iyy0, izz, :, 30) * t1
     1                  + ob12 * u ( ixx0, iyy0, izz, :, 31) * t1
     1                  - ob24 * u ( ixx0, iyy0, izz, :, 32) * t1
     1                 + 0.125 * u ( ixx0, iyy0, izz, :, 33) * t1

! ----------

          U_Rt_Dn ( izz, :) = u ( ixxp1, iyy0, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyy0, izz, :, 2)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 3)
     1                + ob6 * u ( ixxp1, iyy0, izz, :, 5)
     1                + ob6 * u ( ixxp1, iyy0, izz, :, 6)
     1               - 0.25 * u ( ixxp1, iyy0, izz, :, 8)
     1               - ob20 * u ( ixxp1, iyy0, izz, :, 11)
     1               + ob20 * u ( ixxp1, iyy0, izz, :, 12)
     1               + ob12 * u ( ixxp1, iyy0, izz, :, 14)
     1               - ob12 * u ( ixxp1, iyy0, izz, :, 16)
     1                + 0.5 * u ( ixxp1, iyy0, izz, :, 21) * t1
     1                + ob3 * u ( ixxp1, iyy0, izz, :, 22) * t2
     1               + 0.25 * u ( ixxp1, iyy0, izz, :, 23) * t3
     1               - 0.25 * u ( ixxp1, iyy0, izz, :, 24) * t1
     1               + 0.25 * u ( ixxp1, iyy0, izz, :, 25) * t1
     1                - ob6 * u ( ixxp1, iyy0, izz, :, 27) * t2
     1                + ob6 * u ( ixxp1, iyy0, izz, :, 28) * t2
     1               + ob12 * u ( ixxp1, iyy0, izz, :, 30) * t1
     1               + ob12 * u ( ixxp1, iyy0, izz, :, 31) * t1
     1              - 0.125 * u ( ixxp1, iyy0, izz, :, 33) * t1

          U_Rt_Dn_Pt ( izz, :) = u ( ixxp1, iyy0, izz, :, 1)
     1                   - 0.5 * u ( ixxp1, iyy0, izz, :, 2)
     1                   + 0.5 * u ( ixxp1, iyy0, izz, :, 3)
     1                   + ob6 * u ( ixxp1, iyy0, izz, :, 5)
     1                   + ob6 * u ( ixxp1, iyy0, izz, :, 6)
     1                  - ob12 * u ( ixxp1, iyy0, izz, :, 7)
     1                  - 0.25 * u ( ixxp1, iyy0, izz, :, 8)
     1                  - ob20 * u ( ixxp1, iyy0, izz, :, 11)
     1                  + ob20 * u ( ixxp1, iyy0, izz, :, 12)
     1                  + ob12 * u ( ixxp1, iyy0, izz, :, 14)
     1                  - ob12 * u ( ixxp1, iyy0, izz, :, 16)
     1                  + ob24 * u ( ixxp1, iyy0, izz, :, 18)
     1                  - ob24 * u ( ixxp1, iyy0, izz, :, 19)
     1                   + 0.5 * u ( ixxp1, iyy0, izz, :, 21) * t1
     1                   + ob3 * u ( ixxp1, iyy0, izz, :, 22) * t2
     1                  + 0.25 * u ( ixxp1, iyy0, izz, :, 23) * t3
     1                  - 0.25 * u ( ixxp1, iyy0, izz, :, 24) * t1
     1                  + 0.25 * u ( ixxp1, iyy0, izz, :, 25) * t1
     1                   - ob6 * u ( ixxp1, iyy0, izz, :, 27) * t2
     1                   + ob6 * u ( ixxp1, iyy0, izz, :, 28) * t2
     1                  + ob12 * u ( ixxp1, iyy0, izz, :, 30) * t1
     1                  + ob12 * u ( ixxp1, iyy0, izz, :, 31) * t1
     1                  - ob24 * u ( ixxp1, iyy0, izz, :, 32) * t1
     1                 - 0.125 * u ( ixxp1, iyy0, izz, :, 33) * t1

! ----------

          U_Lf_Up ( izz, :) = u ( ixx0, iyyp1, izz, :, 1)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixx0, iyyp1, izz, :, 3)
     1                + ob6 * u ( ixx0, iyyp1, izz, :, 5)
     1                + ob6 * u ( ixx0, iyyp1, izz, :, 6)
     1               - 0.25 * u ( ixx0, iyyp1, izz, :, 8)
     1               + ob20 * u ( ixx0, iyyp1, izz, :, 11)
     1               - ob20 * u ( ixx0, iyyp1, izz, :, 12)
     1               - ob12 * u ( ixx0, iyyp1, izz, :, 14)
     1               + ob12 * u ( ixx0, iyyp1, izz, :, 16)
     1                + 0.5 * u ( ixx0, iyyp1, izz, :, 21) * t1
     1                + ob3 * u ( ixx0, iyyp1, izz, :, 22) * t2
     1               + 0.25 * u ( ixx0, iyyp1, izz, :, 23) * t3
     1               + 0.25 * u ( ixx0, iyyp1, izz, :, 24) * t1
     1               - 0.25 * u ( ixx0, iyyp1, izz, :, 25) * t1
     1                + ob6 * u ( ixx0, iyyp1, izz, :, 27) * t2
     1                - ob6 * u ( ixx0, iyyp1, izz, :, 28) * t2
     1               + ob12 * u ( ixx0, iyyp1, izz, :, 30) * t1
     1               + ob12 * u ( ixx0, iyyp1, izz, :, 31) * t1
     1              - 0.125 * u ( ixx0, iyyp1, izz, :, 33) * t1

          U_Lf_Up_Pt ( izz, :) = u ( ixx0, iyyp1, izz, :, 1)
     1                   + 0.5 * u ( ixx0, iyyp1, izz, :, 2)
     1                   - 0.5 * u ( ixx0, iyyp1, izz, :, 3)
     1                   + ob6 * u ( ixx0, iyyp1, izz, :, 5)
     1                   + ob6 * u ( ixx0, iyyp1, izz, :, 6)
     1                  - ob12 * u ( ixx0, iyyp1, izz, :, 7)
     1                  - 0.25 * u ( ixx0, iyyp1, izz, :, 8)
     1                  + ob20 * u ( ixx0, iyyp1, izz, :, 11)
     1                  - ob20 * u ( ixx0, iyyp1, izz, :, 12)
     1                  - ob12 * u ( ixx0, iyyp1, izz, :, 14)
     1                  + ob12 * u ( ixx0, iyyp1, izz, :, 16)
     1                  - ob24 * u ( ixx0, iyyp1, izz, :, 18)
     1                  + ob24 * u ( ixx0, iyyp1, izz, :, 19)
     1                   + 0.5 * u ( ixx0, iyyp1, izz, :, 21) * t1
     1                   + ob3 * u ( ixx0, iyyp1, izz, :, 22) * t2
     1                  + 0.25 * u ( ixx0, iyyp1, izz, :, 23) * t3
     1                  + 0.25 * u ( ixx0, iyyp1, izz, :, 24) * t1
     1                  - 0.25 * u ( ixx0, iyyp1, izz, :, 25) * t1
     1                   + ob6 * u ( ixx0, iyyp1, izz, :, 27) * t2
     1                   - ob6 * u ( ixx0, iyyp1, izz, :, 28) * t2
     1                  + ob12 * u ( ixx0, iyyp1, izz, :, 30) * t1
     1                  + ob12 * u ( ixx0, iyyp1, izz, :, 31) * t1
     1                  - ob24 * u ( ixx0, iyyp1, izz, :, 32) * t1
     1                 - 0.125 * u ( ixx0, iyyp1, izz, :, 33) * t1

! ----------

          U_Rt_Up ( izz, :) = u ( ixxp1, iyyp1, izz, :, 1)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 2)
     1                - 0.5 * u ( ixxp1, iyyp1, izz, :, 3)
     1                + ob6 * u ( ixxp1, iyyp1, izz, :, 5)
     1                + ob6 * u ( ixxp1, iyyp1, izz, :, 6)
     1               + 0.25 * u ( ixxp1, iyyp1, izz, :, 8)
     1               - ob20 * u ( ixxp1, iyyp1, izz, :, 11)
     1               - ob20 * u ( ixxp1, iyyp1, izz, :, 12)
     1               - ob12 * u ( ixxp1, iyyp1, izz, :, 14)
     1               - ob12 * u ( ixxp1, iyyp1, izz, :, 16)
     1                + 0.5 * u ( ixxp1, iyyp1, izz, :, 21) * t1
     1                + ob3 * u ( ixxp1, iyyp1, izz, :, 22) * t2
     1               + 0.25 * u ( ixxp1, iyyp1, izz, :, 23) * t3
     1               - 0.25 * u ( ixxp1, iyyp1, izz, :, 24) * t1
     1               - 0.25 * u ( ixxp1, iyyp1, izz, :, 25) * t1
     1                - ob6 * u ( ixxp1, iyyp1, izz, :, 27) * t2
     1                - ob6 * u ( ixxp1, iyyp1, izz, :, 28) * t2
     1               + ob12 * u ( ixxp1, iyyp1, izz, :, 30) * t1
     1               + ob12 * u ( ixxp1, iyyp1, izz, :, 31) * t1
     1              + 0.125 * u ( ixxp1, iyyp1, izz, :, 33) * t1

          U_Rt_Up_Pt ( izz, :) = u ( ixxp1, iyyp1, izz, :, 1)
     1                   - 0.5 * u ( ixxp1, iyyp1, izz, :, 2)
     1                   - 0.5 * u ( ixxp1, iyyp1, izz, :, 3)
     1                   + ob6 * u ( ixxp1, iyyp1, izz, :, 5)
     1                   + ob6 * u ( ixxp1, iyyp1, izz, :, 6)
     1                  - ob12 * u ( ixxp1, iyyp1, izz, :, 7)
     1                  + 0.25 * u ( ixxp1, iyyp1, izz, :, 8)
     1                  - ob20 * u ( ixxp1, iyyp1, izz, :, 11)
     1                  - ob20 * u ( ixxp1, iyyp1, izz, :, 12)
     1                  - ob12 * u ( ixxp1, iyyp1, izz, :, 14)
     1                  - ob12 * u ( ixxp1, iyyp1, izz, :, 16)
     1                  + ob24 * u ( ixxp1, iyyp1, izz, :, 18)
     1                  + ob24 * u ( ixxp1, iyyp1, izz, :, 19)
     1                   + 0.5 * u ( ixxp1, iyyp1, izz, :, 21) * t1
     1                   + ob3 * u ( ixxp1, iyyp1, izz, :, 22) * t2
     1                  + 0.25 * u ( ixxp1, iyyp1, izz, :, 23) * t3
     1                  - 0.25 * u ( ixxp1, iyyp1, izz, :, 24) * t1
     1                  - 0.25 * u ( ixxp1, iyyp1, izz, :, 25) * t1
     1                   - ob6 * u ( ixxp1, iyyp1, izz, :, 27) * t2
     1                   - ob6 * u ( ixxp1, iyyp1, izz, :, 28) * t2
     1                  + ob12 * u ( ixxp1, iyyp1, izz, :, 30) * t1
     1                  + ob12 * u ( ixxp1, iyyp1, izz, :, 31) * t1
     1                  - ob24 * u ( ixxp1, iyyp1, izz, :, 32) * t1
     1                 + 0.125 * u ( ixxp1, iyyp1, izz, :, 33) * t1

! --------------------
#endif
! ------------------------------
! Build any entropy-related data for the four states.

          entropy_Lf_Dn ( izz) = entropygr ( ixx0, iyy0, izz, 1)
     1                   + 0.5 * entropygr ( ixx0, iyy0, izz, 2)
     1                   + 0.5 * entropygr ( ixx0, iyy0, izz, 3)
     1              + dt_shift * entropygr ( ixx0, iyy0, izz, 5)

          entropy_Rt_Dn ( izz) = entropygr ( ixxp1, iyy0, izz, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyy0, izz, 2)
     1                   + 0.5 * entropygr ( ixxp1, iyy0, izz, 3)
     1              + dt_shift * entropygr ( ixxp1, iyy0, izz, 5)

          entropy_Lf_Up ( izz) = entropygr ( ixx0, iyyp1, izz, 1)
     1                   + 0.5 * entropygr ( ixx0, iyyp1, izz, 2)
     1                   - 0.5 * entropygr ( ixx0, iyyp1, izz, 3)
     1              + dt_shift * entropygr ( ixx0, iyyp1, izz, 5)

          entropy_Rt_Up ( izz) = entropygr ( ixxp1, iyyp1, izz, 1)
     1                   - 0.5 * entropygr ( ixxp1, iyyp1, izz, 2)
     1                   - 0.5 * entropygr ( ixxp1, iyyp1, izz, 3)
     1              + dt_shift * entropygr ( ixxp1, iyyp1, izz, 5)

! ------------------------------
! Build any EOS-related data for the four states.

          gamma_eos_Lf_Dn ( izz) = gamma_eos_gr ( ixx0, iyy0, izz)
          gamma_sdsp_Lf_Dn ( izz) =
     1                      gamma_soundspeed_gr ( ixx0, iyy0, izz)

          gamma_eos_Rt_Dn ( izz) = gamma_eos_gr ( ixxp1, iyy0, izz)
          gamma_sdsp_Rt_Dn ( izz) =
     1                      gamma_soundspeed_gr ( ixxp1, iyy0, izz)

          gamma_eos_Lf_Up ( izz) = gamma_eos_gr ( ixx0, iyyp1, izz)
          gamma_sdsp_Lf_Up ( izz) =
     1                      gamma_soundspeed_gr ( ixx0, iyyp1, izz)

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

              U_Rt_Dn_Pt ( izz, :) = U_Lf_Dn_Pt ( izz, :)

            ELSE IF ( bcarr_x ( 2, iyy0, izz) == 3) THEN

              U_Rt_Dn ( izz, :) = U_Lf_Dn ( izz, :)
              U_Rt_Dn ( izz, 2) = - U_Rt_Dn ( izz, 2)
              U_Rt_Dn ( izz, 7) = - U_Rt_Dn ( izz, 7)
              U_Rt_Dn ( izz, 8) = - U_Rt_Dn ( izz, 8)

              U_Rt_Dn_Pt ( izz, :) = U_Lf_Dn_Pt ( izz, :)
              U_Rt_Dn_Pt ( izz, 2) = - U_Rt_Dn_Pt ( izz, 2)
              U_Rt_Dn_Pt ( izz, 7) = - U_Rt_Dn_Pt ( izz, 7)
              U_Rt_Dn_Pt ( izz, 8) = - U_Rt_Dn_Pt ( izz, 8)

            END IF

! ----------

            IF ( bcarr_x ( 2, iyyp1, izz) == 2) THEN

              U_Rt_Up ( izz, :) = U_Lf_Up ( izz, :)

              U_Rt_Up_Pt ( izz, :) = U_Lf_Up_Pt ( izz, :)

            ELSE IF ( bcarr_x ( 2, iyyp1, izz) == 3) THEN

              U_Rt_Up ( izz, :) = U_Lf_Up ( izz, :)
              U_Rt_Up ( izz, 2) = - U_Rt_Up ( izz, 2)
              U_Rt_Up ( izz, 7) = - U_Rt_Up ( izz, 7)
              U_Rt_Up ( izz, 8) = - U_Rt_Up ( izz, 8)

              U_Rt_Up_Pt ( izz, :) = U_Lf_Up_Pt ( izz, :)
              U_Rt_Up_Pt ( izz, 2) = - U_Rt_Up_Pt ( izz, 2)
              U_Rt_Up_Pt ( izz, 7) = - U_Rt_Up_Pt ( izz, 7)
              U_Rt_Up_Pt ( izz, 8) = - U_Rt_Up_Pt ( izz, 8)

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

              U_Lf_Dn_Pt ( izz, :) = U_Rt_Dn_Pt ( izz, :)

            ELSE IF ( bcarr_x ( 1, iyy0, izz) == 3) THEN

              U_Lf_Dn ( izz, :) = U_Rt_Dn ( izz, :)
              U_Lf_Dn ( izz, 2) = - U_Lf_Dn ( izz, 2)
              U_Lf_Dn ( izz, 7) = - U_Lf_Dn ( izz, 7)
              U_Lf_Dn ( izz, 8) = - U_Lf_Dn ( izz, 8)

              U_Lf_Dn_Pt ( izz, :) = U_Rt_Dn_Pt ( izz, :)
              U_Lf_Dn_Pt ( izz, 2) = - U_Lf_Dn_Pt ( izz, 2)
              U_Lf_Dn_Pt ( izz, 7) = - U_Lf_Dn_Pt ( izz, 7)
              U_Lf_Dn_Pt ( izz, 8) = - U_Lf_Dn_Pt ( izz, 8)

            END IF

! ----------

            IF ( bcarr_x ( 1, iyyp1, izz) == 2) THEN

              U_Lf_Up ( izz, :) = U_Rt_Up ( izz, :)

              U_Lf_Up_Pt ( izz, :) = U_Rt_Up_Pt ( izz, :)

            ELSE IF ( bcarr_x ( 1, iyyp1, izz) == 3) THEN

              U_Lf_Up ( izz, :) = U_Rt_Up ( izz, :)
              U_Lf_Up ( izz, 2) = - U_Lf_Up ( izz, 2)
              U_Lf_Up ( izz, 7) = - U_Lf_Up ( izz, 7)
              U_Lf_Up ( izz, 8) = - U_Lf_Up ( izz, 8)

              U_Lf_Up_Pt ( izz, :) = U_Rt_Up_Pt ( izz, :)
              U_Lf_Up_Pt ( izz, 2) = - U_Lf_Up_Pt ( izz, 2)
              U_Lf_Up_Pt ( izz, 7) = - U_Lf_Up_Pt ( izz, 7)
              U_Lf_Up_Pt ( izz, 8) = - U_Lf_Up_Pt ( izz, 8)

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

              U_Lf_Up_Pt ( izz, :) = U_Lf_Dn_Pt ( izz, :)

            ELSE IF ( bcarr_y ( ixx0, 2, izz) == 3) THEN

              U_Lf_Up ( izz, :) = U_Lf_Dn ( izz, :)
              U_Lf_Up ( izz, 3) = - U_Lf_Up ( izz, 3)
              U_Lf_Up ( izz, 6) = - U_Lf_Up ( izz, 6)
              U_Lf_Up ( izz, 8) = - U_Lf_Up ( izz, 8)

              U_Lf_Up_Pt ( izz, :) = U_Lf_Dn_Pt ( izz, :)
              U_Lf_Up_Pt ( izz, 3) = - U_Lf_Up_Pt ( izz, 3)
              U_Lf_Up_Pt ( izz, 6) = - U_Lf_Up_Pt ( izz, 6)
              U_Lf_Up_Pt ( izz, 8) = - U_Lf_Up_Pt ( izz, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixxp1, 2, izz) == 2) THEN

              U_Rt_Up ( izz, :) = U_Rt_Dn  ( izz, :)

              U_Rt_Up_Pt ( izz, :) = U_Rt_Dn_Pt  ( izz, :)

            ELSE IF ( bcarr_y ( ixxp1, 2, izz) == 3) THEN

              U_Rt_Up ( izz, :) = U_Rt_Dn  ( izz, :)
              U_Rt_Up ( izz, 3) = - U_Rt_Up ( izz, 3)
              U_Rt_Up ( izz, 6) = - U_Rt_Up ( izz, 6)
              U_Rt_Up ( izz, 8) = - U_Rt_Up ( izz, 8)

              U_Rt_Up_Pt ( izz, :) = U_Rt_Dn_Pt  ( izz, :)
              U_Rt_Up_Pt ( izz, 3) = - U_Rt_Up_Pt ( izz, 3)
              U_Rt_Up_Pt ( izz, 6) = - U_Rt_Up_Pt ( izz, 6)
              U_Rt_Up_Pt ( izz, 8) = - U_Rt_Up_Pt ( izz, 8)

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

              U_Lf_Dn_Pt ( izz, :) = U_Lf_Up_Pt ( izz, :)

            ELSE IF ( bcarr_y ( ixx0, 1, izz) == 3) THEN

              U_Lf_Dn ( izz, :) = U_Lf_Up ( izz, :)
              U_Lf_Dn ( izz, 3) = - U_Lf_Dn ( izz, 3)
              U_Lf_Dn ( izz, 6) = - U_Lf_Dn ( izz, 6)
              U_Lf_Dn ( izz, 8) = - U_Lf_Dn ( izz, 8)

              U_Lf_Dn_Pt ( izz, :) = U_Lf_Up_Pt ( izz, :)
              U_Lf_Dn_Pt ( izz, 3) = - U_Lf_Dn_Pt ( izz, 3)
              U_Lf_Dn_Pt ( izz, 6) = - U_Lf_Dn_Pt ( izz, 6)
              U_Lf_Dn_Pt ( izz, 8) = - U_Lf_Dn_Pt ( izz, 8)

            END IF

! ----------

            IF ( bcarr_y ( ixxp1, 1, izz) == 2) THEN

              U_Rt_Dn  ( izz, :) = U_Rt_Up ( izz, :)

              U_Rt_Dn_Pt  ( izz, :) = U_Rt_Up_Pt ( izz, :)

            ELSE IF ( bcarr_y ( ixxp1, 1, izz) == 3) THEN

              U_Rt_Dn ( izz, :) = U_Rt_Up ( izz, :)
              U_Rt_Dn ( izz, 3) = - U_Rt_Dn ( izz, 3)
              U_Rt_Dn ( izz, 6) = - U_Rt_Dn ( izz, 6)
              U_Rt_Dn ( izz, 8) = - U_Rt_Dn ( izz, 8)

              U_Rt_Dn_Pt ( izz, :) = U_Rt_Up_Pt ( izz, :)
              U_Rt_Dn_Pt ( izz, 3) = - U_Rt_Dn_Pt ( izz, 3)
              U_Rt_Dn_Pt ( izz, 6) = - U_Rt_Dn_Pt ( izz, 6)
              U_Rt_Dn_Pt ( izz, 8) = - U_Rt_Dn_Pt ( izz, 8)

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
#if ( USE_ONLY_1D_RS != 1)
! ------------------------------

          my_blend_HLL = blend_HLL ( izz)
          my_blend_MDHLL = blend_mdhll_to_mdhllc

          U_vrtx ( 1, :) = U_Rt_Up_Pt ( izz, :)
          U_vrtx ( 2, :) = U_Lf_Up_Pt ( izz, :)
          U_vrtx ( 3, :) = U_Lf_Dn_Pt ( izz, :)
          U_vrtx ( 4, :) = U_Rt_Dn_Pt ( izz, :)

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
     1                                 flux_y ( ixx + 1, iyy, izz, :)
     1                               + ob_6_multid * G_Star ( :)

! 100% of the pressure flux will be provided by the 1D RS. The other parts
! of the flux are still multid.

            IF ( split_prsflux_re == 1) flux_y ( ixx + 1, iyy, izz, 3) =
     1    flux_y ( ixx + 1, iyy, izz, 3) - ob_6_multid * my_Prs_Tot_Star

          END IF

! ------------------------------
#endif
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

      SUBROUTINE UPDATE_VARS_CORRECTOR

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes, n_fc_modes,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, dudt, bfxbdy, bfybdy, bfzbdy, flux_x, flux_y, flux_z,
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

! 1) Makes a one-step update for "u" using the
! "flux_x, flux_y, flux_z" variables. If pressure positivity
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
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components) :: dudt

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

      REAL a, b, c, d, e, f, ob_8pi, ob_4pi,
     1     vlx_accn, vly_accn, vlz_accn,

     1     rho_temp, vlx_temp, vly_temp, vlz_temp,
     1     bfx_temp, bfy_temp, bfz_temp,

     1     delx1, delx2, dely1, dely2, delz1, delz2, areaint,

     1     div_min, div_max,
     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax

! ----------------------------------------------------------------------

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      nspecies_min = NFIELD + 2
      nspecies_max = NFIELD + 1 + NFLUID

      IF ( ( expand_by_one_zone .LT. 0) .OR.
     1     ( expand_by_one_zone .GT. 1) ) THEN
      WRITE ( 6, *)
     1    "Stopping in UPDATE_VARS_CORRECTOR : expand_by_one_zone = ",
     1     expand_by_one_zone
      STOP
      END IF

! ----------------------------------------------------------------------
! Update the zone centered variables.
! ----------------------------------------------------------------------

! The "dudt" term has already taken on contributions from the source terms.
! So there is no need to zero it out here.
! The time rate of change of entropy does have to be zeroed out here.

#if ( PRESSURE_POSITIVITY == 1)
!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, nspecies_max, entropygr)      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

          entropygr ( ixx, iyy, izz, 5) = 0.0

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.
!$OMP END PARALLEL DO
#endif
! End of "PRESSURE_POSITIVITY == 1".
! ----------------------------------------------------------------------
! Cartesian

      IF ( igeom .EQ. 1) THEN

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c)                     &
!$OMP&            SHARED ( ix1, iy1, iz1,                               &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, dudt, flux_x, flux_y, flux_z,                              &
!$OMP&    prs_flux_x, prs_flux_y, prs_flux_z,                           &
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! -----

          a = dtcur / dx_save ( ixx)

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) ) * a
#endif

          dudt ( ixx, iyy, izz, 2) = dudt ( ixx, iyy, izz, 2)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a

! -----

          a = dtcur / dy_save ( iyy)

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a
#endif

          dudt ( ixx, iyy, izz, 3) = dudt ( ixx, iyy, izz, 3)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a

! -----

          a = dtcur / dz_save ( izz)

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a
#endif

          dudt ( ixx, iyy, izz, 4) = dudt ( ixx, iyy, izz, 4)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a

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

      STOP'not done yet'

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
!$OMP&    u, dudt, flux_x, flux_y, flux_z,                              &
!$OMP&    prs_flux_x, prs_flux_y, prs_flux_z,                           &
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! -----

          a = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2 )
          a = dtcur / a

          b = xb_save ( ixx)

          c = xb_save ( ixx - 1)

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ) * a
#endif

          a = dtcur / dx_save ( ixx)

          dudt ( ixx, iyy, izz, 2) = dudt ( ixx, iyy, izz, 2)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a

! -----

          a = xc_save ( ixx) * dy_save ( iyy)
          a = dtcur / a

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a
#endif

          a = dtcur / dy_save ( iyy)

          dudt ( ixx, iyy, izz, 3) = dudt ( ixx, iyy, izz, 3)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a

! -----

          a = dtcur / dz_save ( izz)

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a
#endif

          dudt ( ixx, iyy, izz, 4) = dudt ( ixx, iyy, izz, 4)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a

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
!$OMP&            SHARED ( ix1, iy1, iz1, dudt, xc_save)                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          dudt ( ixx, iyy, izz, 3) = dudt ( ixx, iyy, izz, 3)
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

      STOP'not done yet'

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
!$OMP&    u, dudt, flux_x, flux_y, flux_z,                              &
!$OMP&    prs_flux_x, prs_flux_y, prs_flux_z,                           &
!$OMP&    entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z)    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! -----

          a = ( xb_save ( ixx)**3 - xb_save ( ixx - 1)**3 ) / 3.0
          a = dtcur / a

          b = xb_save ( ixx)**2

          c = xb_save ( ixx - 1)**2

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ) * a
#endif

          a = dtcur / dx_save ( ixx)

          dudt ( ixx, iyy, izz, 2) = dudt ( ixx, iyy, izz, 2)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a

! -----

          a = - xc_save ( ixx) * ( COS ( yb_save ( iyy))
     1                           - COS ( yb_save ( iyy - 1)) )
          a = dtcur / a

          b = SIN ( yb_save ( iyy))

          c = SIN ( yb_save ( iyy - 1))

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) * b
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) * c ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz) * b
     1    - entropy_flux_y ( ixx, iyy - 1, izz) * c ) * a
#endif

          a = xc_save ( ixx) * dy_save ( iyy)
          a = dtcur / a

          dudt ( ixx, iyy, izz, 3) = dudt ( ixx, iyy, izz, 3)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a

! -----

          a = xc_save ( ixx) * SIN ( yc_save ( iyy)) * dz_save ( izz)
          a = dtcur / a

          dudt ( ixx, iyy, izz, 1: nspecies_max) =
     1    dudt ( ixx, iyy, izz, 1: nspecies_max)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a
#endif

          a = dtcur / dz_save ( izz)

          dudt ( ixx, iyy, izz, 4) = dudt ( ixx, iyy, izz, 4)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Transform the z-angular momentum density back to a momentum density.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, dudt, xc_save, yc_save)       &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          dudt ( ixx, iyy, izz, 4) = dudt ( ixx, iyy, izz, 4)
     1                    / ( xc_save ( ixx) * SIN ( yc_save ( iyy)) )
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------

! Now, actually make the update. Do it in the same way for all possible
! geometries.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, nspecies_max, u, dudt,        &
!$OMP&    entropygr)                                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1          u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1        + dudt ( ixx, iyy, izz, 1: nspecies_max)

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

          dely1 = xb_save ( ixx) * dy_save ( iyy)
          dely2 = xb_save ( ixx) * dy_save ( iyy)
          delz1 = dz_save ( izz)
          delz2 = dz_save ( izz)

          areaint = xb_save ( ixx) * dy_save ( iyy) * dz_save ( izz)

          IF ( areaint .GT. smallnum) THEN

          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - dely1 * efygr ( ixx, iyy, izz)
     1              + dely2 * efygr ( ixx, iyy, izz - 1)
     1              + delz1 * efzgr ( ixx, iyy, izz)
     1              - delz2 * efzgr ( ixx, iyy - 1, izz)
     1             ) / areaint

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

          delx1 = dx_save ( ixx)
          delx2 = dx_save ( ixx)
          dely1 = xb_save ( ixx) * dy_save ( iyy)
          dely2 = xb_save ( ixx - 1) * dy_save ( iyy)

          areaint = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
     1            * dy_save ( iyy)

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - delx1 * efxgr ( ixx, iyy, izz)
     1              + delx2 * efxgr ( ixx, iyy - 1, izz)
     1              + dely1 * efygr ( ixx, iyy, izz)
     1              - dely2 * efygr ( ixx - 1, iyy, izz)
     1             ) / areaint

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

          dely1 = xb_save ( ixx) * dy_save ( iyy)
          dely2 = xb_save ( ixx) * dy_save ( iyy)
          delz1 = xb_save ( ixx) * SIN ( yb_save ( iyy))
     1          * dz_save ( izz)
          delz2 = xb_save ( ixx) * SIN ( yb_save ( iyy - 1))
     1          * dz_save ( izz)

          areaint = - xb_save ( ixx)**2 * dz_save ( izz)
     1         * ( COS ( yb_save ( iyy)) - COS ( yb_save ( iyy - 1)) )

          IF ( areaint .GT. smallnum) THEN

          bfxbdy ( ixx, iyy, izz, 1) = bfxbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - dely1 * efygr ( ixx, iyy, izz)
     1              + dely2 * efygr ( ixx, iyy, izz - 1)
     1              + delz1 * efzgr ( ixx, iyy, izz)
     1              - delz2 * efzgr ( ixx, iyy - 1, izz)
     1             ) / areaint

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

          delx1 = dx_save ( ixx)
          delx2 = dx_save ( ixx)
          delz1 = xb_save ( ixx) * SIN ( yb_save ( iyy))
     1          * dz_save ( izz)
          delz2 = xb_save ( ixx - 1) * SIN ( yb_save ( iyy))
     1          * dz_save ( izz)

          areaint = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
     1            * SIN ( yb_save ( iyy)) * dz_save ( izz)

          IF ( areaint .GT. smallnum) THEN

          bfybdy ( ixx, iyy, izz, 1) = bfybdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( delx1 * efxgr ( ixx, iyy, izz)
     1            - delx2 * efxgr ( ixx, iyy, izz - 1)
     1            - delz1 * efzgr ( ixx, iyy, izz)
     1            + delz2 * efzgr ( ixx - 1, iyy, izz)
     1             ) / areaint

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

          delx1 = dx_save ( ixx)
          delx2 = dx_save ( ixx)
          dely1 = xb_save ( ixx) * dy_save ( iyy)
          dely2 = xb_save ( ixx - 1) * dy_save ( iyy)

          areaint = 0.5 * ( xb_save ( ixx)**2 - xb_save ( ixx - 1)**2)
     1            * dy_save ( iyy)

          bfzbdy ( ixx, iyy, izz, 1) = bfzbdy ( ixx, iyy, izz, 1)
     1  - dtcur * ( - delx1 * efxgr ( ixx, iyy, izz)
     1              + delx2 * efxgr ( ixx, iyy - 1, izz)
     1              + dely1 * efygr ( ixx, iyy, izz)
     1              - dely2 * efygr ( ixx - 1, iyy, izz)
     1             ) / areaint

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

          a = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a
     1                     + vlxgr ( ixx, iyy, izz) )
          vly_temp = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a
     1                     + vlygr ( ixx, iyy, izz) )
          vlz_temp = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a
     1                     + vlzgr ( ixx, iyy, izz) )

          vlx_accn = 0.0  ! Load in with x-acceleration. Similarly for others.
          vly_accn = 0.0
          vlz_accn = 0.0

          vlx_accn = dtcur * rho_temp * vlx_accn
          vly_accn = dtcur * rho_temp * vly_accn
          vlz_accn = dtcur * rho_temp * vlz_accn

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1) + vlx_accn
          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1) + vly_accn
          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1) + vlz_accn

          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + vlx_accn
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5) + vly_accn
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5) + vlz_accn

          a = vlx_accn * vlx_temp
     1      + vly_accn * vly_temp
     1      + vlz_accn * vlz_temp

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1) + a
          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + a

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

          a = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a
     1                     + vlxgr ( ixx, iyy, izz) )
          vly_temp = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a
     1                     + vlygr ( ixx, iyy, izz) )
          vlz_temp = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a
     1                     + vlzgr ( ixx, iyy, izz) )

          bfx_temp = 0.25 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                      + bfxbdy ( ixxm1, iyy, izz, 1))
     1             + 0.5 * bfxgr ( ixx, iyy, izz)
          bfy_temp = 0.25 * ( bfybdy ( ixx, iyy, izz, 1)
     1                      + bfybdy ( ixx, iyym1, izz, 1))
     1             + 0.5 * bfygr ( ixx, iyy, izz)
          bfz_temp = 0.25 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                      + bfzbdy ( ixx, iyy, izzm1, 1))
     1             + 0.5 * bfzgr ( ixx, iyy, izz)

          vlx_accn = 0.0  ! Load in with x-acceleration. Similarly for others.
          vly_accn = 0.0
          vlz_accn = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn = dtcur * rho_temp * vlx_accn
!          vly_accn = dtcur * rho_temp * vly_accn
!          vlz_accn = dtcur * rho_temp * vlz_accn

            b =  vlx_accn
     1      + ( rho_temp * vly_temp**2 - ob_4pi * bfy_temp**2 )
     1      * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1) + b
          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + b

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1) + vly_accn
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5) + vly_accn

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1) + vlz_accn
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5) + vlz_accn

          c =  vlx_accn * vlx_temp
     1       + vly_accn * vly_temp
     1       + vlz_accn * vlz_temp

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1) + c
          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + c

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

          a = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                     + rhogr ( ixx, iyy, izz) )

          vlx_temp = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a
     1                     + vlxgr ( ixx, iyy, izz) )
          vly_temp = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a
     1                     + vlygr ( ixx, iyy, izz) )
          vlz_temp = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a
     1                     + vlzgr ( ixx, iyy, izz) )

          bfx_temp = 0.25 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                      + bfxbdy ( ixxm1, iyy, izz, 1))
     1             + 0.5 * bfxgr ( ixx, iyy, izz)
          bfy_temp = 0.25 * ( bfybdy ( ixx, iyy, izz, 1)
     1                      + bfybdy ( ixx, iyym1, izz, 1))
     1             + 0.5 * bfygr ( ixx, iyy, izz)
          bfz_temp = 0.25 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                      + bfzbdy ( ixx, iyy, izzm1, 1))
     1             + 0.5 * bfzgr ( ixx, iyy, izz)

          c = COS ( yc_save ( iyy)) / SIN ( yc_save ( iyy))

          vlx_accn = 0.0  ! Load in with x-acceleration. Similarly for others.
          vly_accn = 0.0
          vlz_accn = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn = dtcur * rho_temp * vlx_accn
!          vly_accn = dtcur * rho_temp * vly_accn
!          vlz_accn = dtcur * rho_temp * vlz_accn

          d = + vlx_accn
     1    + ( rho_temp * ( vly_temp**2 + vlz_temp**2 )
     1        - ob_4pi * ( bfy_temp**2 + bfz_temp**2 ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1) + d
          u ( ixx, iyy, izz, 2, 5) = u ( ixx, iyy, izz, 2, 5) + d

          e = vly_accn
     1    - ( rho_temp * ( vlx_temp * vly_temp - vlz_temp**2 * c )
     1        - ob_4pi * ( bfx_temp * bfy_temp - bfz_temp**2 * c ) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1) + e
          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5) + e

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1) + vlz_accn
          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5) + vlz_accn

          f = vlx_accn * vlx_temp
     1      + vly_accn * vly_temp
     1      + vlz_accn * vlz_temp

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1) + f
          u ( ixx, iyy, izz, 5, 5) = u ( ixx, iyy, izz, 5, 5) + f

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE UPDATE_VARS_CORRECTOR

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
     1                  1 - ioffz : iz1 + ioffz, 5) :: entropygr

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

      REAL a, b, c, d, e, ob_8pi, ob_4pi,

     1     delx1, delx2, dely1, dely2, delz1, delz2, areaint

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
!$OMP&    igeom, ob_8pi, cfl_coef, rhofloor, prsfloor,                  &
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

          a = 1.0 / rhogr ( ixx, iyy, izz)

          vlxgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 2, 1) * a
          vlygr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 3, 1) * a
          vlzgr ( ixx, iyy, izz) = u ( ixx, iyy, izz, 4, 1) * a

          bfxgr ( ixx, iyy, izz) =0.5 * ( bfxbdy ( ixx, iyy, izz, 1)
     1                                  + bfxbdy ( ixxm1, iyy, izz, 1))
          bfygr ( ixx, iyy, izz) =0.5 * ( bfybdy ( ixx, iyy, izz, 1)
     1                                  + bfybdy ( ixx, iyym1, izz, 1))
          bfzgr ( ixx, iyy, izz) =0.5 * ( bfzbdy ( ixx, iyy, izz, 1)
     1                                  + bfzbdy ( ixx, iyy, izzm1, 1))

          u ( ixx, iyy, izz, 6, 1) = bfxgr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 7, 1) = bfygr ( ixx, iyy, izz)
          u ( ixx, iyy, izz, 8, 1) = bfzgr ( ixx, iyy, izz)

          e = ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          b = 0.5 * rhogr ( ixx, iyy, izz) * e

          c = ob_8pi
     1      * ( bfxgr ( ixx, iyy, izz)**2 + bfygr ( ixx, iyy, izz)**2
     1        + bfzgr ( ixx, iyy, izz)**2)

          prsgr ( ixx, iyy, izz) = ( u ( ixx, iyy, izz, 5, 1) - b - c )
     1                         * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

! ----------

#if ( PRESSURE_POSITIVITY == 1)
          IF ( prsgr ( ixx, iyy, izz) .LE. prsfloor) THEN

! If the pressure that is evaluated from conserved variables is negative then,
! and then only, replace it with a pressure that is evaluated from the
! entropy density. This is one of ONLY two places where the energy density is
! "touched up". Thus, this is one of ONLY two places where we trade
! conservation in return for physical consistency. The other place is
! in the "RECONST_POS_RHO_PRS_MHD" routines.

          d = prsgr ( ixx, iyy, izz)

          prsgr ( ixx, iyy, izz) = entropygr ( ixx, iyy, izz, 1)
     1 * rhogr ( ixx, iyy, izz)**( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1)
     1    + ( prsgr ( ixx, iyy, izz) - d )
     1    / ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

          END IF
#endif

! ----------

! Even without pressure positivity, we minimally elevate the energy density
! so that the pressure has its floor value.

          IF ( prsgr ( ixx, iyy, izz) .LT. prsfloor) THEN
          prsgr ( ixx, iyy, izz) = prsfloor
          u ( ixx, iyy, izz, 5, 1) = b + c + prsfloor
     1                         / ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
          END IF

! ----------

          msonicgr ( ixx, iyy, izz) =
     1    ( 2.0 * c + gamma_soundspeed_gr ( ixx, iyy, izz)
     1              * prsgr ( ixx, iyy, izz) ) * a
          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

          e = SQRT ( e)

! ----------

          IF ( igeom .EQ. 1) THEN

            d = AMIN1 ( dx_save ( ixx), dy_save ( iyy), dz_save ( izz))

          ELSE IF ( igeom .EQ. 2) THEN

            d = AMIN1 ( dx_save ( ixx), xc_save ( ixx) * dy_save ( iyy),
     1                dz_save ( izz) )

          ELSE IF ( igeom .EQ. 3) THEN

            d = AMIN1 ( dx_save ( ixx), xc_save ( ixx) * dy_save ( iyy),
     1        xc_save ( ixx) * SIN ( yc_save ( iyy)) * dz_save ( izz) )

          END IF

! ----------

          d = cfl_coef * d / ( msonicgr ( ixx, iyy, izz) + e)

          dtcurp1 = AMIN1 ( dtcurp1, d)

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
!$OMP&    bfxbdy, bfybdy, bfzbdy, gamma_eos_gr, gamma_soundspeed_gr,    &
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















