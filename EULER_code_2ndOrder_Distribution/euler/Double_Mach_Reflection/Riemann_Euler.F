! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM RIEMANN_EULER

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This is a very fast, three dimensional, second order accurate, fully
! functional Hydrodynamics code. The spatial interpolation
! is either TVD or WENO. This code was put together to go along with
! Balsara's text book. It can be compiled and run with the
! lines below. Several sample Euler test problems are provided with the code. 

! The Trustees of the Univ. of Notre Dame may assert some rights to
! ownership.

! Code operates in 1, 2 or 3 dimensions. A dimension is switched off
! by putting "ioff? = 0" and "i?1 = 1" in that direction. Thus to switch
! off the z-direction, set "ioffz = 0" and "iz1 = 1". If a direction
! is dynamically active, it is a good idea to set "ioff? = 4" in that
! direction.

! More diagnostic information can be written out per timestep by
! uncommenting lines at the end of "LIMIT_ALL_SLOPES_EULER".

! ----------------------------------------------------------------------

! Compilation lines (Best to run code in 64 bit mode, even on 32 bit CPUs):

! cp Riemann_Euler_src Riemann_Euler.F
! cp Euler_RS_src Euler_RS.F
! cp Eos_Euler_src Eos_Euler.F
! cp Applications_Euler_src Apps.F

! Intel Fortran Compiler:
! Use "-openmp" to enable OpenMP parallelization.
! (To do once-in-a-while OpenMP consistency checks use
! "-openmp -openmp-report2 -par-report3" flags.)
! ( For large problems use "-shared-intel -mcmodel=large" in compile line
! to overcome the 32 bit pointers that are default with ifort.)
! ifort -autodouble -w -O2 -fpp -c Riemann_Euler.F
! ifort -autodouble -w -O2 -fpp -c Euler_RS.F
! ifort -autodouble -w -O2 -fpp -c Eos_Euler.F
! ifort -autodouble -w -O2 -fpp -c Apps.F
! ifort -autodouble -w -O2 -o xeuler Riemann_Euler.o Euler_RS.o Eos_Euler.o \
!        Apps.o

! PathScale Compiler (use -m32 for 32 bit processors):
! Use "-openmp" to enable OpenMP parallelization.
! pathf90 -i8 -r8 -O2 -ftpp -m64 -c Riemann_Euler.F
! pathf90 -i8 -r8 -O2 -ftpp -m64 -c Euler_RS.F
! pathf90 -i8 -r8 -O2 -ftpp -m64 -c Eos_Euler.F
! pathf90 -i8 -r8 -O2 -ftpp -m64 -c Apps.F
! pathf90 -i8 -r8 -O2 -ftpp -m64 -o xeuler Riemann_Euler.o Euler_RS.o \
!        Eos_Euler.o Apps.o

! PGI Compiler:
! Use "-mp" to enable OpenMP parallelization.
! pgf90 -fastsse -i8 -Mr8 -Mr8intrinsics -silent -c Riemann_Euler.F
! pgf90 -fastsse -i8 -Mr8 -Mr8intrinsics -silent -c Euler_RS.F
! pgf90 -fastsse -i8 -Mr8 -Mr8intrinsics -silent -c Eos_Euler.F
! pgf90 -fastsse -i8 -Mr8 -Mr8intrinsics -silent -c Apps.F
! pgf90 -fastsse -i8 -Mr8 -Mr8intrinsics -silent -o xeuler Riemann_Euler.o \
!        Euler_RS.o Eos_Euler.o Apps.o

! GNU fortran Compiler:
! Use "-fopenmp" to enable OpenMP parallelization.
! gfortran -fdefault-integer-8 -fdefault-real-8 -O2 -c Riemann_Euler.F
! gfortran -fdefault-integer-8 -fdefault-real-8 -O2 -c Euler_RS.F 
! gfortran -fdefault-integer-8 -fdefault-real-8 -O2 -c Eos_Euler.F
! gfortran -fdefault-integer-8 -fdefault-real-8 -O2 -c Apps.F
! gfortran -fdefault-integer-8 -fdefault-real-8 -O2 -o xeuler Riemann_Euler.o \
!        Euler_RS.o Eos_Euler.o Apps.o

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        num_level, num_grid, num_pe,
     1        expand_by_one_zone,

     1        n_rad_components, n_rad_modes, n_super_timestep_fld,

     1        ixx, iyy, izz, i_comp, i_rad, ifluid, istep, ntstep,
     1        nimgstep, ndumpstep

! Set size of problem domain.

      PARAMETER ( ix1 = 480, iy1 = 1, iz1 = 120,
     1            ioffx = 4, ioffy = 0, ioffz = 4,
     1            n_cc_components_euler = NFIELD_EULER + NFLUID_EULER,
     1            n_cc_modes_euler = 5)

      PARAMETER ( n_rad_components = NUM_FLD_BINS, n_rad_modes = 5)

      REAL time, dtcur, dtcurp1, timestop, timeimage, dtimage,

     1     xmin, xmax, ymin, ymax, zmin, zmax,

     1     dt_shift,

     1     timing_info1, timing_info2, a

      CHARACTER*8 rhofil_euler, vlxfil_euler, vlyfil_euler,
     1            vlzfil_euler, prsfil_euler,
     1            dumpfil,
     1            speciefil_euler, temp_speciefil_euler,
     1            radfil, temp_radfil

! ------------------------------

! Big static arrays declared past this point.

! "u_euler" holds zone-centered solution as well as all its space-time moments.
! The order of variables in "u" is "density, x,y,z-momentum densities,
! energy density, species density".

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler, n_cc_modes_euler) ::

     1                  u_euler

! "entropy_euler_gr" holds the entropy density and is useful when pressure
! positivity switch is turned on.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes_euler) :: entropy_euler_gr

! "flux_euler_x, flux_euler_y, flux_euler_z" hold the fluxes in each of
! the three directions.

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_cc_components_euler) ::

     1      flux_euler_x, flux_euler_y, flux_euler_z

! "prs_flux_euler_x, prs_flux_euler_y, prs_flux_euler_z" hold the fluxes
! in each of the three directions.

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::

     1      prs_flux_euler_x, prs_flux_euler_y, prs_flux_euler_z,
     1      entropy_flux_euler_x, entropy_flux_euler_y, 
     1      entropy_flux_euler_z

! Zone centered primitive variables and some information for divergence of
! velocity and flattener. These are derived variables.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, msonic_euler_gr, divvel_euler_gr,
     1          flatten_euler_gr

! Real*4 variable used for imaging variables.

      REAL*4, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: r4gr

! ----------
#if ( NUM_FLD_BINS >= 1)
! Declare variables for FLD-based radiative transfer.

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_rad_components, n_rad_modes) :: u_rad

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz,
     1                  n_rad_components) ::
     1        rad_flux_x, rad_flux_y, rad_flux_z

#endif
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

! Courant number: It is provably stable up to 0.33 in 3d; 0.5 in 2d and
! 1.0 in 1d. When Riemann solvers are used in the predictor step
! (which may not be the default state in which the code wakes up.),
! it most likely stable up to 0.45 in 3d; 0.6 in 2d and 1.0 in 1d.
      cfl_coef = 0.5

! If the MC_beta limiter is being used, decide on how much of a compression
! factor one wants to put in. Ranges between 1.0 and 2.0.
! 1.0 makes it MinMod; 2.0 makes it full van Leer MC.
      mc_coef = 2.0

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
! its value.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05

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
      dtimage = 0.2 / 20.0

! Time at which you want the simulation to stop.
      timestop = 0.2

! ------------------------------

! For serial version, set the following to defaults.
      num_level = 0    ! level number is zero.
      num_grid = 0     ! grid number is zero.
      num_pe = 0       ! processor number is zero.

! ------------------------------

! Set default multispecies values, just to make sure they are not NANS.
! For an actual application these will have to be reset to relevant values.

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
      WRITE ( 6, *)"gamma_euler, gasconst_euler, molewt_euler = ",
     1              gamma_euler, gasconst_euler, molewt_euler
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
      WRITE ( 6, *)"HLLC Riemann solver being used."
#else
      WRITE ( 6, *)"HLL Riemann solver being used."
#endif
      WRITE ( 6, *)" "

#if ( PREDICT_WITH_RS == 1)
      WRITE ( 6, *)"Predictor step in ADER done with Riemann solver",
     1             ": slower but more stable."
#else
      WRITE ( 6, *)"Predictor step in ADER done w/o Riemann solver",
     1             ": faster."
#endif
      WRITE ( 6, *)" "

#if ( MULTID_LIMITING == 1)
      WRITE ( 6, *)"Multidimensional TVD limiting will be applied",
     1        " to all variables."
#else
      WRITE ( 6, *)"Multidimensional TVD limiting not being used."
#endif
      WRITE ( 6, *)" "

#if ( PRESSURE_POSITIVITY == 1)
      WRITE ( 6, *)"Positivity of pressure will be enforced as much",
     1          " as possible."
#else
      WRITE ( 6, *)"Positivity of pressure will not be enforced."
#endif
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

! If dimensioning is inconsistent, then write out some info. and stop.

      IF ( ( ix1 .GT. 1) .AND. ( ioffx .LT. 4) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER : increase ioffx= ",ioffx
      STOP
      END IF

      IF ( ( ix1 .EQ. 1) .AND. ( ioffx .NE. 0) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER : ioffx must be 0 ",ioffx
      STOP
      END IF

      IF ( ( iy1 .GT. 1) .AND. ( ioffy .LT. 4) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER : increase ioffy= ",ioffy
      STOP
      END IF

      IF ( ( iy1 .EQ. 1) .AND. ( ioffy .NE. 0) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER : ioffy must be 0 ",ioffy
      STOP
      END IF

      IF ( ( iz1 .GT. 1) .AND. ( ioffz .LT. 4) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER : increase ioffz= ",ioffz
      STOP
      END IF

      IF ( ( iz1 .EQ. 1) .AND. ( ioffz .NE. 0) ) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER : ioffz must be 0 ", ioffz
      STOP
      END IF

      IF ( NTRACERS_EULER > NFLUID_EULER) THEN
      WRITE ( 6, *)"Stopping in RIEMANN_EULER: NTRACER","S > NFLUI","D",
     1              NTRACERS_EULER, NFLUID_EULER
      STOP
      END IF

! ----------------------------------------------------------------------

! Set the size of the patch being simulated.

      xmin = 0.0
      xmax = 4.0
      ymin = -1.0
      ymax = 1.0
      zmin = 0.0
      zmax = 1.0

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

      rhofil_euler = "rhoa0001"
      vlxfil_euler = "vlxa0001"
      vlyfil_euler = "vlya0001"
      vlzfil_euler = "vlza0001"
      prsfil_euler = "prsa0001"

      speciefil_euler = "sp010001"

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
! See "PAD_EULER_BOUNDARY_X,Y,Z" for meanings of boundary condition arrays.
! Boundary condition number 17 is just reserved for "do nothing", i.e.
! for situations when an external driver sets the boundary conditions.

! ------------------------------

! Bottom and top yz planes.

      DO iyy = 1 - ioffy, iy1 + ioffy

        DO izz = 1 - ioffz, iz1 + ioffz

        x_indx_limits ( 1, iyy, izz) = 1
        x_indx_limits ( 2, iyy, izz) = ix1

        IF ( ABS ( xmin - base_grid_xmin) .LT. smallnum) THEN
          bcarr_x ( 1, iyy, izz) = 1
        ELSE
          bcarr_x ( 1, iyy, izz) = 17
        END IF

        IF ( ABS ( xmax - base_grid_xmax) .LT. smallnum) THEN
          bcarr_x ( 2, iyy, izz) = 2
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
          bcarr_y ( ixx, 1, izz) = 2
        ELSE
          bcarr_y ( ixx, 1, izz) = 17
        END IF

        IF ( ABS ( ymax - base_grid_ymax) .LT. smallnum) THEN
          bcarr_y ( ixx, 2, izz) = 2
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

        IF ( xc_save ( ixx) .LE. 0.16666666666667) THEN
          bcarr_z ( ixx, iyy, 1) = 1
        ELSE
          bcarr_z ( ixx, iyy, 1) = 3
        END IF

        ELSE
          bcarr_z ( ixx, iyy, 1) = 17
        END IF

        IF ( ABS ( zmax - base_grid_zmax) .LT. smallnum) THEN
          bcarr_z ( ixx, iyy, 2) = 17
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

      CALL INIT_DOUBLE_MACH_2D_XZ_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        dtcur,

     1        u_euler,

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

      READ ( 5, *) time, dtcur, i_comp, i_rad

      READ ( 5, *) dumpfil

      READ ( 5, *) rhofil_euler
      READ ( 5, *) vlxfil_euler
      READ ( 5, *) vlyfil_euler
      READ ( 5, *) vlzfil_euler
      READ ( 5, *) prsfil_euler
      READ ( 5, *) speciefil_euler
      READ ( 5, *) radfil

      CLOSE ( 5)

      IF ( i_comp .NE. n_cc_components_euler) THEN
      WRITE ( 6, *)
     1      "Stopping because : i_comp .NE. n_cc_components_euler ",
     1              i_comp, n_cc_components_euler
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


      READ ( 98) ( ( ( ( u_euler ( ixx, iyy, izz, i_comp, 1),
     1    ixx = 1, ix1), iyy = 1, iy1), izz = 1, iz1),
     1    i_comp = 1, n_cc_components_euler)

#if ( NUM_FLD_BINS >= 1)
      READ ( 98) ( ( ( ( u_rad ( ixx, iyy, izz, i_rad, 1), ixx=1, ix1),
     1    iyy = 1, iy1), izz = 1, iz1), i_rad = 1, n_rad_components)
#endif

      CLOSE ( 98)

! ----------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *) dumpfil,
     1             " read in with time, dtcur = ", time, dtcur
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Image file names = ", rhofil_euler," ",
     1            vlxfil_euler," ", vlyfil_euler," ", vlzfil_euler," ",
     1            prsfil_euler," ", speciefil_euler," ", radfil
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

      WRITE ( 98) ( ( ( ( u_euler ( ixx, iyy, izz, i_comp, 1),
     1    ixx = 1, ix1), iyy = 1, iy1), izz = 1, iz1),
     1    i_comp = 1, n_cc_components_euler)

#if ( NUM_FLD_BINS >= 1)
      WRITE ( 98) ( ( ( ( u_rad ( ixx, iyy, izz, i_rad, 1),ixx=1, ix1),
     1    iyy = 1, iy1), izz = 1, iz1), i_rad = 1, n_rad_components)
#endif

      CLOSE ( 98)

! ----------

      WRITE ( 6, *) dumpfil,
     1      " written out with time, dtcur, n_cc_components_euler = ",
     1      time, dtcur, n_cc_components_euler
      WRITE ( 6, *)" ix1, iy1, iz1 = ", ix1, iy1, iz1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Image file names = ", rhofil_euler," ",
     1            vlxfil_euler," ", vlyfil_euler," ", vlzfil_euler," ",
     1            prsfil_euler," ", speciefil_euler," ", radfil
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

! ----------

      OPEN ( unit = 5, file = "riemann.in")

      WRITE ( 5, *) time, dtcur, n_cc_components_euler, n_rad_components

      WRITE ( 5, *) dumpfil

      WRITE ( 5, *) rhofil_euler
      WRITE ( 5, *) vlxfil_euler
      WRITE ( 5, *) vlyfil_euler
      WRITE ( 5, *) vlzfil_euler
      WRITE ( 5, *) prsfil_euler
      WRITE ( 5, *) speciefil_euler
      WRITE ( 5, *) radfil

      CLOSE ( 5)

! ----------

      WRITE ( 6, *) time, dtcur, n_cc_components_euler, n_rad_components

      WRITE ( 6, *) dumpfil

      WRITE ( 6, *) rhofil_euler
      WRITE ( 6, *) vlxfil_euler
      WRITE ( 6, *) vlyfil_euler
      WRITE ( 6, *) vlzfil_euler
      WRITE ( 6, *) prsfil_euler
      WRITE ( 6, *) speciefil_euler
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

      CALL PAD_EULER_BOUNDARY_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        time,

     1        u_euler, bcarr_x, x_indx_limits,

     1        xc_save, yc_save, zc_save)

      CALL PAD_EULER_BOUNDARY_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        time,

     1        u_euler, bcarr_y, y_indx_limits,

     1        xc_save, yc_save, zc_save)

      CALL PAD_EULER_BOUNDARY_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        time,

     1        u_euler, bcarr_z, z_indx_limits,

     1        xc_save, yc_save, zc_save)

      CALL ADD_DOUBLE_MACH_2D_XZ_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        xmin, xmax, ymin, ymax, zmin, zmax,

     1        time,

     1        u_euler,

     1        xb_save, dx_save, xc_save,
     1        yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Limit the slopes of all the flow variables. The time rate of update
! "u ( ixx, iyy, izz, :, 5)" is also zeroed out here. The flattener
! is also built. The primitive variables are built.

      CALL LIMIT_ALL_SLOPES_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,

     1        u_euler, entropy_euler_gr,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

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

      WRITE ( 6, *)"Writing image files = ", rhofil_euler," ",
     1            vlxfil_euler," ", vlyfil_euler," ", vlzfil_euler," ",
     1            prsfil_euler," ", speciefil_euler
      WRITE ( 6, *)" "

! -----

      OPEN ( unit = 99, file = rhofil_euler)
      r4gr ( :, :, :) = rho_euler_gr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlxfil_euler)
      r4gr ( :, :, :) = vlx_euler_gr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlyfil_euler)
      r4gr ( :, :, :) = vly_euler_gr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlzfil_euler)
      r4gr ( :, :, :) = vlz_euler_gr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = prsfil_euler)
      r4gr ( :, :, :) = prs_euler_gr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

! -----

#if ( NFLUID_EULER >= 1)

      temp_speciefil_euler = speciefil_euler

      DO ifluid = 1, NFLUID_EULER

      OPEN ( unit = 99, file = temp_speciefil_euler)
      r4gr ( :, :, :) = u_euler ( :, :, :, NFIELD_EULER + ifluid, 1)
     1                / rho_euler_gr ( :, :, :)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

      CALL INCNAME_SP ( temp_speciefil_euler, 1)

      END DO

#endif
! -----
#if ( NUM_FLD_BINS >= 1)

      temp_radfil = radfil

      DO i_rad = 1, n_rad_components

      OPEN ( unit = 99, file = temp_radfil)
      r4gr ( :, :, :) = u_rad ( :, :, :, i_rad, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ix1),
     1                     iyy = 1, iy1), izz = 1, iz1)
      CLOSE ( 99)

      CALL INCNAME_SP ( temp_radfil, 1)

      END DO

#endif
! ----------

      CALL INCNAME_8 ( rhofil_euler, 1)
      CALL INCNAME_8 ( vlxfil_euler, 1)
      CALL INCNAME_8 ( vlyfil_euler, 1)
      CALL INCNAME_8 ( vlzfil_euler, 1)
      CALL INCNAME_8 ( prsfil_euler, 1)

      CALL INCNAME_8 ( speciefil_euler, 1)

      CALL INCNAME_8 ( radfil, 1)

! ----------------------------------------------------------------------
      END IF ! End of "( time .GT. timestop)..." condition.
! ----------------------------------------------------------------------

! This is the principal stopping point in this code.

      IF ( time .GE. timestop) THEN
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Stopping in RIEMANN_EULER: time .GE. timestop"
      WRITE ( 6, *)"time, timestop = ", time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      STOP
      END IF

! ----------------------------------------------------------------------

! Build fluxes in all directions at "dt_shift = 0.0" and use that to
! build the first, i.e. predictor, stage. The fluxes don't need to be zeroed
! and are only used to build the time rate of update
! "u ( ixx, iyy, izz, :, 5)" within a zone.

! This is a predictor step and can be done either by using full-up
! Riemann solvers at zone boundaries or ( more inexpensively) by using
! the conserved variables and their linear variation within the zone.

! For stringent problems, the former is more stabilizing than the latter.
! Consequently, the code defaults to a situation where we use
! full-up Riemann solvers at zone boundaries.

! ------------------------------
#if ( PREDICT_WITH_RS == 1)
! ------------------------------

      dt_shift = 0.0

      CALL MAKE_EULER_FLUX_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u_euler, flux_euler_x, prs_flux_euler_x,
     1        entropy_euler_gr, entropy_flux_euler_x,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      dt_shift = 0.0

      CALL MAKE_EULER_FLUX_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u_euler, flux_euler_y, prs_flux_euler_y,
     1        entropy_euler_gr, entropy_flux_euler_y,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

      dt_shift = 0.0

      CALL MAKE_EULER_FLUX_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u_euler, flux_euler_z, prs_flux_euler_z,
     1        entropy_euler_gr, entropy_flux_euler_z,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------
#else
! ------------------------------

      CALL MAKE_ADER_UPDT_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur,

     1        u_euler, entropy_euler_gr,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------
#endif
! End of "PREDICT_WITH_RS"
! ----------------------------------------------------------------------

! Add in the source terms to the predictor.

      CALL ADD_PREDICTOR_SOURCE_TERMS_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur,

     1        u_euler,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

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

! Zero out "flux_euler_x"
      flux_euler_x ( :, :, :, :) = 0.0
      prs_flux_euler_x ( :, :, :) = 0.0
      entropy_flux_euler_x ( :, :, :) = 0.0

      dt_shift = 0.5

      CALL MAKE_EULER_FLUX_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u_euler, flux_euler_x, prs_flux_euler_x,
     1        entropy_euler_gr, entropy_flux_euler_x,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Zero out "flux_euler_y"
      flux_euler_y ( :, :, :, :) = 0.0
      prs_flux_euler_y ( :, :, :) = 0.0
      entropy_flux_euler_y ( :, :, :) = 0.0

      dt_shift = 0.5

      CALL MAKE_EULER_FLUX_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u_euler, flux_euler_y, prs_flux_euler_y,
     1        entropy_euler_gr, entropy_flux_euler_y,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ------------------------------

! Zero out "flux_euler_z"
      flux_euler_z ( :, :, :, :) = 0.0
      prs_flux_euler_z ( :, :, :) = 0.0
      entropy_flux_euler_z ( :, :, :) = 0.0

      dt_shift = 0.5

      CALL MAKE_EULER_FLUX_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u_euler, flux_euler_z, prs_flux_euler_z,
     1        entropy_euler_gr, entropy_flux_euler_z,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Update the flow variables in this corrector step and get new "dtcurp1".

      CALL UPDATE_VARS_CORRECTOR_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        expand_by_one_zone,
     1        dtcur,

     1        u_euler, flux_euler_x, flux_euler_y, flux_euler_z,
     1        prs_flux_euler_x, prs_flux_euler_y, prs_flux_euler_z,
     1        entropy_euler_gr, entropy_flux_euler_x, 
     1        entropy_flux_euler_y, entropy_flux_euler_z,


     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------

! Make the next timestep that will be minimized in the ensuing subroutine.
      dtcurp1 = 1.25 * dtcur

      CALL EVALUATE_TIMESTEP_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components_euler, n_cc_modes_euler,
     1        dtcur, dtcurp1,

     1        u_euler, entropy_euler_gr,

     1        rho_euler_gr, prs_euler_gr, vlx_euler_gr, vly_euler_gr,
     1          vlz_euler_gr,
     1        gamma_eos_euler_gr, gamma_soundspeed_euler_gr,
     1          mean_mol_wt_euler_gr,
     1        temperature_euler_gr, divvel_euler_gr, msonic_euler_gr,
     1          flatten_euler_gr,

     1        x_indx_limits,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------

! Update "time" and write out a small amount of timestep info.
! Gather timing information.
! If the "time" is about to exceed "timestop" in the next time step,
! curtail the next timestep so as to stop exactly at "timestop".

      CALL CPU_TIME ( timing_info2)

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

      END PROGRAM RIEMANN_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PAD_EULER_BOUNDARY_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        time,

     1        u, bcarr_x, x_indx_limits,

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
     1        n_cc_components, n_cc_modes

      REAL time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

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
     1     species_tmp ( 0: NFLUID_EULER)

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, species_tmp)                                         &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    x_indx_limits, bcarr_x, u)                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - ioffz, iz1 + ioffz

      DO iyy = 1 - ioffy, iy1 + ioffy

      oned_begin = x_indx_limits ( 1, iyy, izz)
      oned_end = x_indx_limits ( 2, iyy, izz)

! ----------------------------------------------------------------------

! Inflow at lower x-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 8.0
      momx_tmp = 8.0 * 4.125 * SQRT ( 3.0)
      momy_tmp = 0.0
      momz_tmp = - 8.0 * 4.125
      eng_tmp = 116.5 / ( gamma_euler - 1.0) + 0.5 * 8.0 * 8.25 ** 2
      species_tmp ( :) = 0.0

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 1) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx

      u ( ilhs, iyy, izz, 1, 1) = rho_tmp
      u ( ilhs, iyy, izz, 2, 1) = momx_tmp
      u ( ilhs, iyy, izz, 3, 1) = momy_tmp
      u ( ilhs, iyy, izz, 4, 1) = momz_tmp
      u ( ilhs, iyy, izz, 5, 1) = eng_tmp
#if ( NFLUID_EULER >= 1)
      u ( ilhs, iyy, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID_EULER)
#endif

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at lower x-boundary.

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 2) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx
      irhs = oned_begin

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

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

      END DO

      END IF

! ------------------------------

! Periodic at lower x-boundary.

      IF ( bcarr_x ( 1, iyy, izz) .EQ. 4) THEN

      DO ixx = 1, ioffx

      ilhs = oned_begin - ixx
      irhs = oned_end + 1 - ixx

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, species_tmp)                                         &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    x_indx_limits, bcarr_x, u)                                    &
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
      species_tmp ( :) = 0.0

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 1) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx

      u ( ilhs, iyy, izz, 1, 1) = rho_tmp
      u ( ilhs, iyy, izz, 2, 1) = momx_tmp
      u ( ilhs, iyy, izz, 3, 1) = momy_tmp
      u ( ilhs, iyy, izz, 4, 1) = momz_tmp
      u ( ilhs, iyy, izz, 5, 1) = eng_tmp
#if ( NFLUID_EULER >= 1)
      u ( ilhs, iyy, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID_EULER)
#endif

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at upper x-boundary.

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 2) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx
      irhs = oned_end

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

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

      END DO

      END IF

! ------------------------------

! Periodic at upper x-boundary.

      IF ( bcarr_x ( 2, iyy, izz) .EQ. 4) THEN

      DO ixx = 1, ioffx

      ilhs = oned_end + ixx
      irhs = oned_begin + ixx - 1

      u ( ilhs, iyy, izz, :, 1) = u ( irhs, iyy, izz, :, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE PAD_EULER_BOUNDARY_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PAD_EULER_BOUNDARY_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        time,

     1        u, bcarr_y, y_indx_limits,

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
     1        n_cc_components, n_cc_modes

      REAL time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

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
     1     species_tmp ( 0: NFLUID_EULER)

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, species_tmp)                                         &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    y_indx_limits, bcarr_y, u)                                    &
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
      species_tmp ( :) = 0.0

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 1) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy

      u ( ixx, ilhs, izz, 1, 1) = rho_tmp
      u ( ixx, ilhs, izz, 2, 1) = momx_tmp
      u ( ixx, ilhs, izz, 3, 1) = momy_tmp
      u ( ixx, ilhs, izz, 4, 1) = momz_tmp
      u ( ixx, ilhs, izz, 5, 1) = eng_tmp
#if ( NFLUID_EULER >= 1)
      u ( ixx, ilhs, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID_EULER)
#endif

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at lower y-boundary.

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 2) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy
      irhs = oned_begin

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

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

      END DO

      END IF

! ------------------------------

! Periodic at lower y-boundary.

      IF ( bcarr_y ( ixx, 1, izz) .EQ. 4) THEN

      DO iyy = 1, ioffy

      ilhs = oned_begin - iyy
      irhs = oned_end + 1 - iyy

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, species_tmp)                                         &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    y_indx_limits, bcarr_y, u)                                    &
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
      species_tmp ( :) = 0.0

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 1) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy

      u ( ixx, ilhs, izz, 1, 1) = rho_tmp
      u ( ixx, ilhs, izz, 2, 1) = momx_tmp
      u ( ixx, ilhs, izz, 3, 1) = momy_tmp
      u ( ixx, ilhs, izz, 4, 1) = momz_tmp
      u ( ixx, ilhs, izz, 5, 1) = eng_tmp
#if ( NFLUID_EULER >= 1)
      u ( ixx, ilhs, izz, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID_EULER)
#endif

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at upper y-boundary.

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 2) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy
      irhs = oned_end

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

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

      END DO

      END IF

! ------------------------------

! Periodic at upper y-boundary.

      IF ( bcarr_y ( ixx, 2, izz) .EQ. 4) THEN

      DO iyy = 1, ioffy

      ilhs = oned_end + iyy
      irhs = oned_begin + iyy - 1

      u ( ixx, ilhs, izz, :, 1) = u ( ixx, irhs, izz, :, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE PAD_EULER_BOUNDARY_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE PAD_EULER_BOUNDARY_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        time,

     1        u, bcarr_z, z_indx_limits,

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
     1        n_cc_components, n_cc_modes

      REAL time

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

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
     1     species_tmp ( 0: NFLUID_EULER)

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, species_tmp)                                         &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    z_indx_limits, bcarr_z, u)                                    &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 1 - ioffy, iy1 + ioffy

      DO ixx = 1 - ioffx, ix1 + ioffx

      oned_begin = z_indx_limits ( ixx, iyy, 1)
      oned_end = z_indx_limits ( ixx, iyy, 2)

! ----------------------------------------------------------------------

! Inflow at lower z-boundary. Fill in below variables with desired (space
! and time-dependent) values.

      rho_tmp = 8.0
      momx_tmp = 8.0 * 4.125 * SQRT ( 3.0)
      momy_tmp = 0.0
      momz_tmp = - 8.0 * 4.125
      eng_tmp = 116.5 / ( gamma_euler - 1.0) + 0.5 * 8.0 * 8.25 ** 2
      species_tmp ( :) = 0.0

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 1) THEN

      DO izz = 1, ioffz

      ilhs = oned_begin - izz

      u ( ixx, iyy, ilhs, 1, 1) = rho_tmp
      u ( ixx, iyy, ilhs, 2, 1) = momx_tmp
      u ( ixx, iyy, ilhs, 3, 1) = momy_tmp
      u ( ixx, iyy, ilhs, 4, 1) = momz_tmp
      u ( ixx, iyy, ilhs, 5, 1) = eng_tmp
#if ( NFLUID_EULER >= 1)
      u ( ixx, iyy, ilhs, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID_EULER)
#endif

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at lower z-boundary.

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 2) THEN

      DO izz = 1, ioffz

      ilhs = oned_begin - izz
      irhs = oned_begin

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

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

      END DO

      END IF

! ------------------------------

! Periodic at lower z-boundary.

      IF ( bcarr_z ( ixx, iyy, 1) .EQ. 4) THEN

      DO izz = 1, ioffz

      ilhs = oned_begin - izz
      irhs = oned_end + 1 - izz

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, ilhs, irhs,                  &
!$OMP&    oned_begin, oned_end, rho_tmp, momx_tmp, momy_tmp, momz_tmp,  &
!$OMP&    eng_tmp, species_tmp)                                         &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, nspecies_min, nspecies_max,                            &
!$OMP&    z_indx_limits, bcarr_z, u)                                    &
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
      species_tmp ( :) = 0.0

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 1) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz

      u ( ixx, iyy, ilhs, 1, 1) = rho_tmp
      u ( ixx, iyy, ilhs, 2, 1) = momx_tmp
      u ( ixx, iyy, ilhs, 3, 1) = momy_tmp
      u ( ixx, iyy, ilhs, 4, 1) = momz_tmp
      u ( ixx, iyy, ilhs, 5, 1) = eng_tmp
#if ( NFLUID_EULER >= 1)
      u ( ixx, iyy, ilhs, nspecies_min: nspecies_max, 1) =
     1                     species_tmp ( 1: NFLUID_EULER)
#endif

      END DO

      END IF

! ------------------------------

! Outflow/Continuitive at upper z-boundary.

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 2) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz
      irhs = oned_end

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

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

      END DO

      END IF

! ------------------------------

! Periodic at upper z-boundary.

      IF ( bcarr_z ( ixx, iyy, 2) .EQ. 4) THEN

      DO izz = 1, ioffz

      ilhs = oned_end + izz
      irhs = oned_begin + izz - 1

      u ( ixx, iyy, ilhs, :, 1) = u ( ixx, iyy, irhs, :, 1)

      END DO

      END IF

! ----------------------------------------------------------------------

      END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE PAD_EULER_BOUNDARY_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LIMIT_ALL_SLOPES_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,

     1        u, entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Obtains primitive variables from conserved variables.

! Limits the conserved variables "u" in all three directions.
! Makes their undivided differences (slopes).
! The conserved variables in "u" are "density, x,y,z-momentum densities,
! energy density, species density".

! Characteristic limiting can also be done in addition to limiting on
! the conserved variables. It does increase the cost but provides more
! stability. Set LIMIT_TWICE to 1 to invoke additional limiting.

! If pressure positivity switch is on, it also builds up an "entropygr"
! variable which contains the entropy density. The slopes of this variable
! are set conservatively and it is evolved in time with the other variables.
! The caveat is that it is evolved using physical pressures unless the
! pressures go negative. If that happens, the pressures are supplemented
! by the entropy variable.

! SLOPE_LIMITER if set to 1, gives TVD limiting.
! Otherwise, gives r=3 centered WENO limiting.

! "flattengr" applies a flattener in regions of strong shock.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      INTEGER ixx, iyy, izz, i_comp,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: workgr

      REAL psi_limiter, central_zone_wt, ob_6, r13_12,
     1     MinMod, MC_beta, ta, tb,
     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, flatten_1d,

     1     a, b, c, d, e, f, g, h,
     1     var_zone_min, var_zone_max,
     1     var_neighbor_min, var_neighbor_max,

     1     u_ip2, u_ip1, u_i, u_im1, u_im2,
     1     R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3,

     1     dx, dy, dz, ob_dx, ob_dy, ob_dz

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD_EULER + NFLUID_EULER, 2) :: u_1d

! MinMod limiter
!      MinMod ( ta, tb) = AMIN1 ( ABS ( ta), ABS ( tb))
!     1                    * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( ta, tb, mc_coef) = AMIN1 ( 0.5 * ( ABS ( ta + tb)),
     1                      mc_coef * ABS ( ta), mc_coef * ABS ( tb) )
     1                   * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! ----------------------------------------------------------------------

! Variables that go through the EOS subroutine boundary.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      rho, eng, prs, temperature, entropy,
     1      mean_mol_wt, gamma_eos, gamma_soundspeed

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie

! ----------------------------------------------------------------------

      ob_6 = 1.0 / 6.0

      r13_12 = 13.0 / 12.0

      central_zone_wt = 100.0

! ----------

      IF ( n_cc_components .NE. NFIELD_EULER + NFLUID_EULER) THEN
      WRITE ( 6, *)
     1      "Stopping in LIMIT_ALL_SLOPES_EULER : n_cc_components = ",
     1              n_cc_components
      STOP
      END IF

      IF ( n_cc_modes .NE. 5) THEN
      WRITE ( 6, *)
     1      "Stopping in LIMIT_ALL_SLOPES_EULER : n_cc_modes = ",
     1              n_cc_modes
      STOP
      END IF

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES_EULER : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES_EULER : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES_EULER : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

      IF ( NFIELD_EULER .NE. 5) THEN
      WRITE ( 6, *)
     1  "Stopping in LIMIT_ALL_SLOPES_EULER : NFIELD_","EULER"
      WRITE ( 6, *) NFIELD_EULER, ".NE. 5"
      STOP
      END IF

! ----------------------------------------------------------------------

! Make zone centered primitive variables.

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, rho, eng, specie, prs,               &
!$OMP&   temperature, entropy, gamma_eos, gamma_soundspeed, mean_mol_wt)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gasconst_euler, gamma_euler, smallnum, prsfloor, rhofloor,    &
!$OMP&    molewt_euler,                                                 &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, entropygr,              &
!$OMP&    temperaturegr, gamma_eos_gr, gamma_soundspeed_gr,             &
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

          b ( ixx) = 0.5 * rhogr ( ixx, iyy, izz)
     1      * ( vlxgr ( ixx, iyy, izz)**2 + vlygr ( ixx, iyy, izz)**2
     1        + vlzgr ( ixx, iyy, izz)**2)

          rho ( ixx) = rhogr ( ixx, iyy, izz)

          eng ( ixx) = ( u ( ixx, iyy, izz, 5, 1) - b ( ixx) )

#if ( NFLUID_EULER >= 1)
          specie ( 1:NFLUID_EULER, ixx) =
     1         u ( ixx, iyy, izz, 6: NFIELD_EULER + NFLUID_EULER, 1)
     1       * a ( ixx)
#endif

          END DO ! End of DO ixx = loop.

! ----------

! Choose between a simple or more complicated EOS. Species fraction
! can be strictly restricted to lie between 0 and 1 if it is copied back.

#if ( NFLUID_EULER == 0)

          DO ixx = 1 - ioffx, ix1 + ioffx

          gamma_eos ( ixx) = gamma_euler
          gamma_soundspeed ( ixx) = gamma_euler
          mean_mol_wt ( ixx) = molewt_euler

          prs ( ixx) = AMAX1 ( eng ( ixx) * ( gamma_euler - 1.0),
     1                                                  prsfloor)
          temperature ( ixx) = molewt_euler * prs ( ixx)
     1                       / ( gasconst_euler * rho ( ixx))

          END DO ! End of DO ixx = loop.

#else

      CALL RHO_ENG_IN_PRS_T_GAMA_OUT_EULER
     1        ( 1, ix1, ioffx,

     1       gasconst_euler, gamma_euler, smallnum, prsfloor, rhofloor,
     1   cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,

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

#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz, 6: NFIELD_EULER + NFLUID_EULER, 1) =
     1    specie ( 1:NFLUID_EULER, ixx) * rhogr ( ixx, iyy, izz)
#endif

          msonicgr ( ixx, iyy, izz) =
     1     gamma_soundspeed ( ixx) * prs ( ixx) * a ( ixx)   ! sound speed
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
! TVD/r=3 WENO limiters at second order start here.
! ----------------------------------------------------------------------

! Now make the slopes (with flattening) for the zone-centered
! conserved variables. We do not need to do this for the magnetic field
! components because they will be done later.

! ------------------------------
#if ( SLOPE_LIMITER == 1)
! ------------------------------
! MC_beta Limiter.

      DO i_comp = 1, n_cc_components

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c)                                   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    i_comp, mc_coef, u, flattengr)                                &
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

          a ( ixx) = u ( ixxp1, iyy, izz, i_comp, 1)
     1             - u ( ixx, iyy, izz, i_comp, 1)

          b ( ixx) = u ( ixx, iyy, izz, i_comp, 1)
     1             - u ( ixxm1, iyy, izz, i_comp, 1)

          END DO

! --------------------

          DO ixx = 1 - ioffx, ix1 + ioffx

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

      END DO ! End of DO i_comp = loop.

! ------------------------------
#else
! ------------------------------
! WENO Limiter.

      DO i_comp = 1, n_cc_components

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, izzp2, izzm2,            &
!$OMP& iyy, iyyp1, iyym1, iyyp2, iyym2, ixx, ixxp1, ixxm1, ixxp2, ixxm2,&
!$OMP&    a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                     &
!$OMP&    R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3)               &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    i_comp, central_zone_wt, smallnum, r13_12, u, flattengr)      &
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

          c ( ixx) = 1.0 - flattengr ( ixx, iyy, izz)

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

! --------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

      END DO ! End of DO i_comp = loop.

! ------------------------------
#endif
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

      psi_limiter = 0.5    ! For tough problems, use this value for all vars.
!      IF ( i_comp .EQ. 1) psi_limiter = 1.0

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, izzp1, izzm1, iyy, iyyp1, iyym1,       &
!$OMP&    ixx, ixxp1, ixxm1, a, b, c, d, e, f, g, h,                    &
!$OMP&    var_zone_min, var_zone_max, var_neighbor_min,var_neighbor_max)&
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    i_comp, smallnum, psi_limiter, u)                             &
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

! ------------------------------

      END DO ! End of DO i_comp = loop.

! ----------------------------------------------------------------------
#endif
! End of "MULTID_LIMITING == 1"
! ----------------------------------------------------------------------
! If we additionally want limiting on the characteristic variables then we
! do the below three loop nests.
#if ( LIMIT_TWICE == 1)
! ----------------------------------------------------------------------

! Characteristic limiting in the x-direction.

      IF ( ioffx .GT. 0) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a,                           &
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar,                       &
!$OMP&    speciebar, flatten_1d, u_1d)                                  &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,    &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, flattengr)              &
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
#if ( NFLUID_EULER >= 1)
          u_1d ( ixx, 6: NFIELD_EULER + NFLUID_EULER, 1) =
     1                       u ( ixx, iyy, izz, 6: n_cc_components, 1)
#endif

          u_1d ( ixx, 1, 2) = u ( ixx, iyy, izz, 1, 2)
          u_1d ( ixx, 2, 2) = u ( ixx, iyy, izz, 2, 2)
          u_1d ( ixx, 3, 2) = u ( ixx, iyy, izz, 3, 2)
          u_1d ( ixx, 4, 2) = u ( ixx, iyy, izz, 4, 2)
          u_1d ( ixx, 5, 2) = u ( ixx, iyy, izz, 5, 2)
#if ( NFLUID_EULER >= 1)
          u_1d ( ixx, 6: NFIELD_EULER + NFLUID_EULER, 2) = 
     1                       u ( ixx, iyy, izz, 6: n_cc_components, 2)
#endif

          rhobar ( ixx) = rhogr ( ixx, iyy, izz)
          prsbar ( ixx) = prsgr ( ixx, iyy, izz)

          vlxbar ( ixx) = vlxgr ( ixx, iyy, izz)
          vlybar ( ixx) = vlygr ( ixx, iyy, izz)
          vlzbar ( ixx) = vlzgr ( ixx, iyy, izz)

#if ( NFLUID_EULER >= 1)
          a ( ixx) = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( 1: NFLUID_EULER, ixx) =
     1           u ( ixx, iyy, izz, 6: n_cc_components, 1) * a ( ixx)
#endif

          flatten_1d ( ixx) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO ixx = loop.

! ----------

          CALL LIMIT_1D_EULER_CONS ( 1, ix1, ioffx,

     1     gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,
     1     cp_specie_euler, molewt_specie_euler,
     1       heat_of_formation_euler,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     speciebar, flatten_1d, u_1d)

! ----------

          DO ixx = 1 - ioffx, ix1 + ioffx

          u ( ixx, iyy, izz, 1, 2) = u_1d ( ixx, 1, 2)
          u ( ixx, iyy, izz, 2, 2) = u_1d ( ixx, 2, 2)
          u ( ixx, iyy, izz, 3, 2) = u_1d ( ixx, 3, 2)
          u ( ixx, iyy, izz, 4, 2) = u_1d ( ixx, 4, 2)
          u ( ixx, iyy, izz, 5, 2) = u_1d ( ixx, 5, 2)
#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz, 6: n_cc_components, 2) =
     1             u_1d ( ixx, 6: NFIELD_EULER + NFLUID_EULER, 2)
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
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar,                       &
!$OMP&    speciebar, flatten_1d, u_1d)                                  &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,    &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, flattengr)              &
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
#if ( NFLUID_EULER >= 1)
          u_1d ( iyy, 6: NFIELD_EULER + NFLUID_EULER, 1) =
     1                       u ( ixx, iyy, izz, 6: n_cc_components, 1)
#endif

          u_1d ( iyy, 1, 2) = u ( ixx, iyy, izz, 1, 3)
          u_1d ( iyy, 2, 2) = u ( ixx, iyy, izz, 3, 3)
          u_1d ( iyy, 3, 2) = u ( ixx, iyy, izz, 4, 3)
          u_1d ( iyy, 4, 2) = u ( ixx, iyy, izz, 2, 3)
          u_1d ( iyy, 5, 2) = u ( ixx, iyy, izz, 5, 3)
#if ( NFLUID_EULER >= 1)
          u_1d ( iyy, 6: NFIELD_EULER + NFLUID_EULER, 2) =
     1                       u ( ixx, iyy, izz, 6: n_cc_components, 3)
#endif

          rhobar ( iyy) = rhogr ( ixx, iyy, izz)
          prsbar ( iyy) = prsgr ( ixx, iyy, izz)

          vlxbar ( iyy) = vlygr ( ixx, iyy, izz)
          vlybar ( iyy) = vlzgr ( ixx, iyy, izz)
          vlzbar ( iyy) = vlxgr ( ixx, iyy, izz)

#if ( NFLUID_EULER >= 1)
          a ( iyy) = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( 1: NFLUID_EULER, iyy) =
     1           u ( ixx, iyy, izz, 6: n_cc_components, 1) * a ( iyy)
#endif

          flatten_1d ( iyy) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO iyy = loop.

! ----------

          CALL LIMIT_1D_EULER_CONS ( 1, iy1, ioffy,

     1     gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,
     1     cp_specie_euler, molewt_specie_euler,
     1       heat_of_formation_euler,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     speciebar, flatten_1d, u_1d)

! ----------

          DO iyy = 1 - ioffy, iy1 + ioffy

          u ( ixx, iyy, izz, 1, 3) = u_1d ( iyy, 1, 2)
          u ( ixx, iyy, izz, 3, 3) = u_1d ( iyy, 2, 2)
          u ( ixx, iyy, izz, 4, 3) = u_1d ( iyy, 3, 2)
          u ( ixx, iyy, izz, 2, 3) = u_1d ( iyy, 4, 2)
          u ( ixx, iyy, izz, 5, 3) = u_1d ( iyy, 5, 2)
#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz, 6: n_cc_components, 3) =
     1             u_1d ( iyy, 6: NFIELD_EULER + NFLUID_EULER, 2)
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
!$OMP&    rhobar, prsbar, vlxbar, vlybar, vlzbar,                       &
!$OMP&    speciebar, flatten_1d, u_1d)                                  &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,    &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    mc_coef, central_zone_wt, n_cc_components,                    &
!$OMP&    u, rhogr, prsgr, vlxgr, vlygr, vlzgr, flattengr)              &
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
#if ( NFLUID_EULER >= 1)
          u_1d ( izz, 6: NFIELD_EULER + NFLUID_EULER, 1) =
     1                       u ( ixx, iyy, izz, 6: n_cc_components, 1)
#endif

          u_1d ( izz, 1, 2) = u ( ixx, iyy, izz, 1, 4)
          u_1d ( izz, 2, 2) = u ( ixx, iyy, izz, 4, 4)
          u_1d ( izz, 3, 2) = u ( ixx, iyy, izz, 2, 4)
          u_1d ( izz, 4, 2) = u ( ixx, iyy, izz, 3, 4)
          u_1d ( izz, 5, 2) = u ( ixx, iyy, izz, 5, 4)
#if ( NFLUID_EULER >= 1)
          u_1d ( izz, 6: NFIELD_EULER + NFLUID_EULER, 2) =
     1                       u ( ixx, iyy, izz, 6: n_cc_components, 4)
#endif

          rhobar ( izz) = rhogr ( ixx, iyy, izz)
          prsbar ( izz) = prsgr ( ixx, iyy, izz)

          vlxbar ( izz) = vlzgr ( ixx, iyy, izz)
          vlybar ( izz) = vlxgr ( ixx, iyy, izz)
          vlzbar ( izz) = vlygr ( ixx, iyy, izz)

#if ( NFLUID_EULER >= 1)
          a ( izz) = 1.0 / rhogr ( ixx, iyy, izz)
          speciebar ( 1: NFLUID_EULER, izz) =
     1           u ( ixx, iyy, izz, 6: n_cc_components, 1) * a ( izz)
#endif

          flatten_1d ( izz) = flattengr ( ixx, iyy, izz)

          END DO ! End of DO izz = loop.

! ----------

          CALL LIMIT_1D_EULER_CONS ( 1, iz1, ioffz,

     1     gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,
     1     cp_specie_euler, molewt_specie_euler,
     1       heat_of_formation_euler,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     speciebar, flatten_1d, u_1d)

! ----------

          DO izz = 1 - ioffz, iz1 + ioffz

          u ( ixx, iyy, izz, 1, 4) = u_1d ( izz, 1, 2)
          u ( ixx, iyy, izz, 4, 4) = u_1d ( izz, 2, 2)
          u ( ixx, iyy, izz, 2, 4) = u_1d ( izz, 3, 2)
          u ( ixx, iyy, izz, 3, 4) = u_1d ( izz, 4, 2)
          u ( ixx, iyy, izz, 5, 4) = u_1d ( izz, 5, 2)
#if ( NFLUID_EULER >= 1)
          u ( ixx, iyy, izz, 6: n_cc_components, 4) =
     1             u_1d ( izz, 6: NFIELD_EULER + NFLUID_EULER, 2)
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

! For problems with very strong shocks, one
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

! This is an opportunity to evaluate and write out the min and
! max of the flow variables if that is desired.

#if ( WRITE_DIAGNOSTICS == 1)

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

! ------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, dx, dy, dz,                  &
!$OMP&    ob_dx, ob_dy, ob_dz, a, b, c)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    rhogr, prsgr, vlxgr, vlygr, vlzgr)                            &
!$OMP&           REDUCTION ( MIN: rhomin) REDUCTION ( MAX: rhomax)      &
!$OMP&           REDUCTION ( MIN: vlxmin) REDUCTION ( MAX: vlxmax)      &
!$OMP&           REDUCTION ( MIN: vlymin) REDUCTION ( MAX: vlymax)      &
!$OMP&           REDUCTION ( MIN: vlzmin) REDUCTION ( MAX: vlzmax)      &
!$OMP&           REDUCTION ( MIN: prsmin) REDUCTION ( MAX: prsmax)      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1

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

! ----------

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! ------------------------------

      WRITE ( 6, *)"From LIMIT_ALL_SLOPES_EULER: = "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"rhomin, rhomax, prsmin, prsmax = ",
     1              rhomin, rhomax, prsmin, prsmax
      WRITE ( 6, *)" "
      WRITE ( 6, *)"vlxmin, vlxmax, vlymin, vlymax, vlzmin, vlzmax = ",
     1              vlxmin, vlxmax, vlymin, vlymax, vlzmin, vlzmax
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE LIMIT_ALL_SLOPES_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE LIMIT_1D_EULER_CONS ( oned_begin, oned_end, oned_off,

     1     gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,
     1     cp_specie_euler, molewt_specie_euler,
     1       heat_of_formation_euler,

     1     mc_coef, central_zone_wt,

     1     rhobar, prsbar, vlxbar, vlybar, vlzbar,
     1     speciebar, flatten_1d, u_1d)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Does characteristic limiting for the Euler variables. The variables are
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

      REAL gamma_euler, smallnum, prsfloor, rhofloor, gasconst_euler,
     1     mc_coef, central_zone_wt

      REAL, DIMENSION ( 0 : NFLUID_EULER) ::
     1     cp_specie_euler, molewt_specie_euler, heat_of_formation_euler

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhobar, prsbar, vlxbar, vlybar, vlzbar, flatten_1d

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     speciebar

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD_EULER + NFLUID_EULER, 2) :: u_1d

! ----------------------------------------------------------------------

      INTEGER nriembegin, nriemend, iprimitive, nsys, ix, j, k, ifield,
     1        debug_stop

      REAL r13_12, ta, tb, MinMod, MC_beta

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b,
     1     u_ip2, u_ip1, u_i, u_im1, u_im2,
     1     R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3,
     1     rh, en, vx, vy, vz

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD_EULER + NFLUID_EULER) ::
     1     eigenval, del_eigwt, del_u

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD_EULER + NFLUID_EULER, NFIELD_EULER + NFLUID_EULER) ::
     1     lf_eigenvec, rt_eigenvec

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  NFIELD_EULER + NFLUID_EULER, -2 : 2) :: eigwt

! MinMod limiter
      MinMod ( ta, tb) = AMIN1 ( ABS ( ta), ABS ( tb))
     1                    * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( ta, tb, mc_coef) = AMIN1 ( 0.5 * ( ABS ( ta + tb)),
     1                      mc_coef * ABS ( ta), mc_coef * ABS ( tb) )
     1                    * ( SIGN ( 0.5, ta) + SIGN ( 0.5, tb))

! ----------------------------------------------------------------------

      COMMON /  LIMIT_1D_EULER_CONS_COM / a, b,

     1     u_ip2, u_ip1, u_i, u_im1, u_im2,
     1     R1_x, R2_x, R3_x, Is1, Is2, Is3, Wt1, Wt2, Wt3,
     1     rh, en, vx, vy, vz,

     1     eigenval, del_eigwt, del_u,

     1     lf_eigenvec, rt_eigenvec,

     1     eigwt,

     1     nriembegin, nriemend, iprimitive, nsys, ix, j, k, ifield,
     1     debug_stop,
     1     r13_12

!$OMP THREADPRIVATE (/LIMIT_1D_EULER_CONS_COM/)

! ----------------------------------------------------------------------

      nriembegin = oned_begin - oned_off
      nriemend = oned_end + oned_off
      iprimitive = 0

      nsys = NFIELD_EULER + NFLUID_EULER

      r13_12 = 13.0 / 12.0

      debug_stop = 0

! ------------------------------

      CALL EIGENVEC_EULER

     1             ( nriembegin, nriemend, iprimitive,

     1               gamma_euler, smallnum, rhofloor, prsfloor,
     1               debug_stop,

     1               gasconst_euler,
     1               cp_specie_euler, molewt_specie_euler,
     1                  heat_of_formation_euler,

     1               rhobar, prsbar,
     1               vlxbar, vlybar, vlzbar,
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
          Is1 ( ix) = r13_12 * ( - 2.0 * u_im1 ( ix) + u_im2 ( ix)
     1                           + u_i ( ix))**2
     1              + 0.25 * ( - 4.0 * u_im1 ( ix) + u_im2 ( ix)
     1                         + 3.0 * u_i ( ix))**2

          R2_x ( ix) = 0.5 * ( u_ip1 ( ix) - u_im1 ( ix) )
          Is2 ( ix) = r13_12 * ( u_im1 ( ix) - 2.0 * u_i ( ix)
     1                         + u_ip1 ( ix) )**2
     1              + 0.25 * ( u_im1 ( ix) - u_ip1 ( ix) )**2

          R3_x ( ix) = - 1.5 * u_i ( ix) + 2.0 * u_ip1 ( ix)
     1                 - 0.5 * u_ip2 ( ix)
          Is3 ( ix) = r13_12 * ( u_i ( ix) - 2.0 * u_ip1 ( ix)
     1                         + u_ip2 ( ix))**2
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
     1               + Wt2 ( ix) * R2_x ( ix) + Wt3 ( ix) * R3_x ( ix)

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

        IF ( ( rh ( ix) .LE. 0.0) .OR.
     1       ( en ( ix) - 0.5 * ( vx ( ix)**2 + vy ( ix)**2
     1                  + vz ( ix)**2) / rh ( ix) .LE. 0.0) ) THEN
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

        IF ( ( rh ( ix) .LE. 0.0) .OR.
     1       ( en ( ix) - 0.5 * ( vx ( ix)**2 + vy ( ix)**2
     1                  + vz ( ix)**2) / rh ( ix) .LE. 0.0) ) THEN
        u_1d ( ix, :, 2) = 0.0
        END IF

      END DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE LIMIT_1D_EULER_CONS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_ADER_UPDT_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone,
     1        dtcur, 

     1        u, entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and produces an ADER-based time
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

! This subroutine is effectively Stage I, i.e. the predictor stage, ! of a two-stage update.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
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

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b, c

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho_pl, rho_mi, prs_pl, prs_mi,
     1     eng_pl, eng_mi,

     1     vlx_pl, vlx_mi, vly_pl, vly_mi,
     1     vlz_pl, vlz_mi,

     1     prs_flux_pl, prs_flux_mi, entropy_pl, entropy_mi

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD_EULER + NFLUID_EULER) ::
     1     flux_pl, flux_mi

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_ADER_UPDT_EULER : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_ADER_UPDT_EULER : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_ADER_UPDT_EULER : expand 1d arrays"
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
!$OMP&    rho_pl, prs_pl, eng_pl, vlx_pl, vly_pl, vlz_pl,               &
!$OMP&    flux_pl, prs_flux_pl, entropy_pl,                             &
!$OMP&    rho_mi, prs_mi, eng_mi, vlx_mi, vly_mi, vlz_mi,               &
!$OMP&    flux_mi, prs_flux_mi, entropy_mi)                             &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, smallnum, rhofloor, prsfloor,                          &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, entropygr, gamma_eos_gr)                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

! ----------------------------------------------------------------------

! Obtain primitive variables; evaluate fluxes and add them in at x-faces.

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

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

! ----------

          rho_pl ( ixx) = AMAX1 ( rhofloor, rho_pl ( ixx))
          a ( ixx) = 1.0 / rho_pl ( ixx)

          vlx_pl ( ixx) = vlx_pl ( ixx) * a ( ixx)
          vly_pl ( ixx) = vly_pl ( ixx) * a ( ixx)
          vlz_pl ( ixx) = vlz_pl ( ixx) * a ( ixx)

          b ( ixx) = 0.5 * rho_pl ( ixx) * ( vlx_pl ( ixx)**2
     1      + vly_pl ( ixx)**2 + vlz_pl ( ixx)**2)

          prs_pl ( ixx) = ( eng_pl ( ixx) - b ( ixx) )
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

          prs_mi ( ixx) = ( eng_mi ( ixx) - b ( ixx) )
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

          prs_flux_pl ( ixx) = prs_pl ( ixx)

          flux_pl ( ixx, 3) = rho_pl ( ixx) * vlx_pl ( ixx)
     1                                      * vly_pl ( ixx)

          flux_pl ( ixx, 4) = rho_pl ( ixx) * vlx_pl ( ixx)
     1                                      * vlz_pl ( ixx)

          flux_pl ( ixx, 5) = ( eng_pl ( ixx) + prs_pl ( ixx) )
     1                      * vlx_pl ( ixx)

#if ( NFLUID_EULER >= 1)
          flux_pl ( ixx, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 2) )
     1          * vlx_pl ( ixx)
#endif

! ----------

          flux_mi ( ixx, 1) = rho_mi ( ixx) * vlx_mi ( ixx)

          flux_mi ( ixx, 2) = rho_mi ( ixx) * vlx_mi ( ixx)**2

          prs_flux_mi ( ixx) = prs_mi ( ixx)

          flux_mi ( ixx, 3) = rho_mi ( ixx) * vlx_mi ( ixx)
     1                                      * vly_mi ( ixx)

          flux_mi ( ixx, 4) = rho_mi ( ixx) * vlx_mi ( ixx)
     1                                      * vlz_mi ( ixx)

          flux_mi ( ixx, 5) = ( eng_mi ( ixx) + prs_mi ( ixx) )
     1                      * vlx_mi ( ixx)

#if ( NFLUID_EULER >= 1)
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

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

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

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

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

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

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
!$OMP&    rho_pl, prs_pl, eng_pl, vlx_pl, vly_pl, vlz_pl,               &
!$OMP&    flux_pl, prs_flux_pl, entropy_pl,                             &
!$OMP&    rho_mi, prs_mi, eng_mi, vlx_mi, vly_mi, vlz_mi,               &
!$OMP&    flux_mi, prs_flux_mi, entropy_mi)                             &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, smallnum, rhofloor, prsfloor,                          &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, entropygr, gamma_eos_gr)                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

! ----------------------------------------------------------------------

! Obtain primitive variables; evaluate fluxes and add them in at y-faces.

          DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

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

! ----------

          rho_pl ( iyy) = AMAX1 ( rhofloor, rho_pl ( iyy))
          a ( iyy) = 1.0 / rho_pl ( iyy)

          vlx_pl ( iyy) = vlx_pl ( iyy) * a ( iyy)
          vly_pl ( iyy) = vly_pl ( iyy) * a ( iyy)
          vlz_pl ( iyy) = vlz_pl ( iyy) * a ( iyy)

          b ( iyy) = 0.5 * rho_pl ( iyy) * ( vlx_pl ( iyy)**2
     1      + vly_pl ( iyy)**2 + vlz_pl ( iyy)**2)

          prs_pl ( iyy) = ( eng_pl ( iyy) - b ( iyy) )
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

          prs_mi ( iyy) = ( eng_mi ( iyy) - b ( iyy) )
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

          flux_pl ( iyy, 3) = rho_pl ( iyy) * vly_pl ( iyy)**2

          prs_flux_pl ( iyy) = prs_pl ( iyy)

          flux_pl ( iyy, 4) = rho_pl ( iyy) * vly_pl ( iyy)
     1                      * vlz_pl ( iyy)

          flux_pl ( iyy, 5) = ( eng_pl ( iyy) + prs_pl ( iyy) )
     1                      * vly_pl ( iyy)

#if ( NFLUID_EULER >= 1)
          flux_pl ( iyy, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 3) )
     1          * vly_pl ( iyy)
#endif

! ----------

          flux_mi ( iyy, 1) = rho_mi ( iyy) * vly_mi ( iyy)

          flux_mi ( iyy, 2) = rho_mi ( iyy) * vly_mi ( iyy)
     1                      * vlx_mi ( iyy)

          flux_mi ( iyy, 3) = rho_mi ( iyy) * vly_mi ( iyy)**2

          prs_flux_mi ( iyy) = prs_mi ( iyy)

          flux_mi ( iyy, 4) = rho_mi ( iyy) * vly_mi ( iyy)
     1                      * vlz_mi ( iyy)

          flux_mi ( iyy, 5) = ( eng_mi ( iyy) + prs_mi ( iyy) )
     1                      * vly_mi ( iyy)

#if ( NFLUID_EULER >= 1)
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

          DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

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

          DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

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

          DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

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
!$OMP&    rho_pl, prs_pl, eng_pl, vlx_pl, vly_pl, vlz_pl,               &
!$OMP&    flux_pl, prs_flux_pl, entropy_pl,                             &
!$OMP&    rho_mi, prs_mi, eng_mi, vlx_mi, vly_mi, vlz_mi,               &
!$OMP&    flux_mi, prs_flux_mi, entropy_mi)                             &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, smallnum, rhofloor, prsfloor,                          &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, entropygr, gamma_eos_gr)                                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)


      DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

        DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

! ----------------------------------------------------------------------

! Obtain primitive variables; evaluate fluxes and add them in at z-faces.

          DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

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

! ----------

          rho_pl ( izz) = AMAX1 ( rhofloor, rho_pl ( izz))
          a ( izz) = 1.0 / rho_pl ( izz)

          vlx_pl ( izz) = vlx_pl ( izz) * a ( izz)
          vly_pl ( izz) = vly_pl ( izz) * a ( izz)
          vlz_pl ( izz) = vlz_pl ( izz) * a ( izz)

          b ( izz) = 0.5 * rho_pl ( izz) * ( vlx_pl ( izz)**2
     1      + vly_pl ( izz)**2 + vlz_pl ( izz)**2)

          prs_pl ( izz) = ( eng_pl ( izz) - b ( izz) )
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

          prs_mi ( izz) = ( eng_mi ( izz) - b ( izz) )
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

          flux_pl ( izz, 3) = rho_pl ( izz) * vlz_pl ( izz)
     1                      * vly_pl ( izz)

          flux_pl ( izz, 4) = rho_pl ( izz) * vlz_pl ( izz)**2

          prs_flux_pl ( izz) = prs_pl ( izz)

          flux_pl ( izz, 5) = ( eng_pl ( izz) + prs_pl ( izz) )
     1                      * vlz_pl ( izz)

#if ( NFLUID_EULER >= 1)
          flux_pl ( izz, nspecies_min: nspecies_max) =
     1            ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1      + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 4) )
     1          * vlz_pl ( izz)
#endif

! ----------

          flux_mi ( izz, 1) = rho_mi ( izz) * vlz_mi ( izz)

          flux_mi ( izz, 2) = rho_mi ( izz) * vlz_mi ( izz)
     1                      * vlx_mi ( izz)

          flux_mi ( izz, 3) = rho_mi ( izz) * vlz_mi ( izz)
     1                      * vly_mi ( izz)

          flux_mi ( izz, 4) = rho_mi ( izz) * vlz_mi ( izz)**2

          prs_flux_mi ( izz) = prs_mi ( izz)

          flux_mi ( izz, 5) = ( eng_mi ( izz) + prs_mi ( izz) )
     1                      * vlz_mi ( izz)

#if ( NFLUID_EULER >= 1)
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

          DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

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

          DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

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

      END SUBROUTINE MAKE_ADER_UPDT_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_EULER_FLUX_X

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, flux_x, prs_flux_x,
     1        entropygr, entropy_flux_x,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        x_indx_limits, bcarr_x,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and puts out "flux_x and prs_flux_x".

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_x".

! Since this is an ADER-style update, it is suggested that the update
! be done in two stages:

! Stage I : With the slopes evaluated, enter with
! "dt_shift = 0.0". This gives
! fluxes that are used to obtain the time rate of update
! "u ( ixx, iyy, izz, :, 5)" within a zone.

! Stage II : With the slopes and the time rate of update
! evaluated, enter with "dt_shift = 0.5".
! This reloads the flux with time-centered
! values, suitable for a spatially and temporally second order accurate update.

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
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
     1                  n_cc_components) :: flux_x

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1     prs_flux_x, entropy_flux_x

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      INTEGER ixx, iyy, izz, nriembegin, nriemend, one_x,
     1        iprimitive, use_approx_eos, split_prsflux_re, debug_stop,
     1        nspecies_min, nspecies_max

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b, c

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,
     1     entropyremi, entropyrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr, lapidus_arr

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD_EULER + NFLUID_EULER, NFIELD_EULER + NFLUID_EULER) ::
     1    lf_eigenvec_re, rt_eigenvec_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD_EULER + NFLUID_EULER) ::
     1     eigenval_re, eigenwt_re, conserved_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD_EULER + NFLUID_EULER + 1) ::
     1     flux_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD_EULER + NFLUID_EULER) ::
     1      zonebdy_flux, zonebdy_conserved

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------

! The computation is not dynamically active in the x-direction
! ( as indicated by setting "ioffx = 0" on the outside) so return.

      IF ( ioffx .EQ. 0) RETURN

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffx .LT. ONEDFLUIDBEGIN) .OR.
     1     ( ix1 + ioffx .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_EULER_FLUX_X : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffx, present_base = ",
     1              0 - ioffx, ONEDFLUIDBEGIN
      WRITE ( 6, *)"ix1 + ioffx, present_bound = ",
     1              ix1 + ioffx, ONEDFLUIDEND
      STOP
      END IF

      one_x = MIN0 ( expand_by_one_zone, ioffx)

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

! ----------

      iprimitive = 0     ! Eigenstructure constructed in conserved variables.
      use_approx_eos = 1 ! Just use simple EOS with constant polytropic index.

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c, debug_stop,         &
!$OMP&    nriembegin, nriemend,                                         &
!$OMP&    rhoremi, prsremi, vlxremi, vlyremi, vlzremi, specieremi,      &
!$OMP&    mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,      &
!$OMP&    entropyremi, msonic_arr, divvel_arr, lapidus_arr,             &
!$OMP&    rhorepl, prsrepl, vlxrepl, vlyrepl, vlzrepl, specierepl,      &
!$OMP&    mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,      &
!$OMP&    entropyrepl, zonebdy_conserved, zonebdy_flux,                 &
!$OMP&    eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,      &
!$OMP&    conserved_re, flux_re, entropy_flux_re,                       &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, specie0, do_ein)                &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, one_x, expand_by_one_zone,                             &
!$OMP&    dtcur, dt_shift, nspecies_min, nspecies_max,                  &
!$OMP&    gamma_euler, gasconst_euler, smallnum, rhofloor, prsfloor,    &
!$OMP&    ientropyfix, iprimitive, use_approx_eos, fracdiff,            &
!$OMP&    split_prsflux_re, divratio_einfeldt_rs,                       &
!$OMP&    lapidus_threshold, lapidus_coef,                              &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    x_indx_limits, bcarr_x, u, flux_x, prs_flux_x, entropygr,     &
!$OMP&    entropy_flux_x, mean_mol_wt_gr, gamma_eos_gr,                 &
!$OMP&    gamma_soundspeed_gr, msonicgr, divvelgr)                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( expand_by_one_zone, ioffz),
     1       iz1 + MIN0 ( expand_by_one_zone, ioffz)

        DO iyy = 1 - MIN0 ( expand_by_one_zone, ioffy),
     1         iy1 + MIN0 ( expand_by_one_zone, ioffy)

        debug_stop = 0

! ----------------------------------------------------------------------

! Load in the left and right states for the Riemann solver.

          nriembegin = x_indx_limits ( 1, iyy, izz) - 1 - one_x
          nriemend = x_indx_limits ( 2, iyy, izz) + one_x

          DO ixx = nriembegin, nriemend

! ----------

          rhoremi ( ixx) = u ( ixx, iyy, izz, 1, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 1, 2)
     1        + dt_shift * u ( ixx, iyy, izz, 1, 5)

          vlxremi ( ixx) = u ( ixx, iyy, izz, 2, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 2, 2)
     1        + dt_shift * u ( ixx, iyy, izz, 2, 5)

          vlyremi ( ixx) = u ( ixx, iyy, izz, 3, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 3, 2)
     1        + dt_shift * u ( ixx, iyy, izz, 3, 5)

          vlzremi ( ixx) = u ( ixx, iyy, izz, 4, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 4, 2)
     1        + dt_shift * u ( ixx, iyy, izz, 4, 5)

          prsremi ( ixx) = u ( ixx, iyy, izz, 5, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 5, 2)
     1        + dt_shift * u ( ixx, iyy, izz, 5, 5)

          mean_mol_wt_remi ( ixx) = mean_mol_wt_gr ( ixx, iyy, izz)
          gamma_eos_remi ( ixx) = gamma_eos_gr ( ixx, iyy, izz)
          gamma_soundspeed_remi ( ixx) =
     1                         gamma_soundspeed_gr ( ixx, iyy, izz)

! ----------

          rhorepl ( ixx) = u ( ixx + 1, iyy, izz, 1, 1)
     1             - 0.5 * u ( ixx + 1, iyy, izz, 1, 2)
     1        + dt_shift * u ( ixx + 1, iyy, izz, 1, 5)

          vlxrepl ( ixx) = u ( ixx + 1, iyy, izz, 2, 1)
     1             - 0.5 * u ( ixx + 1, iyy, izz, 2, 2)
     1        + dt_shift * u ( ixx + 1, iyy, izz, 2, 5)

          vlyrepl ( ixx) = u ( ixx + 1, iyy, izz, 3, 1)
     1             - 0.5 * u ( ixx + 1, iyy, izz, 3, 2)
     1        + dt_shift * u ( ixx + 1, iyy, izz, 3, 5)

          vlzrepl ( ixx) = u ( ixx + 1, iyy, izz, 4, 1)
     1             - 0.5 * u ( ixx + 1, iyy, izz, 4, 2)
     1        + dt_shift * u ( ixx + 1, iyy, izz, 4, 5)

          prsrepl ( ixx) = u ( ixx + 1, iyy, izz, 5, 1)
     1             - 0.5 * u ( ixx + 1, iyy, izz, 5, 2)
     1        + dt_shift * u ( ixx + 1, iyy, izz, 5, 5)

          mean_mol_wt_repl ( ixx) = mean_mol_wt_gr ( ixx + 1, iyy, izz)
          gamma_eos_repl ( ixx) = gamma_eos_gr ( ixx + 1, iyy, izz)
          gamma_soundspeed_repl ( ixx) =
     1                         gamma_soundspeed_gr ( ixx + 1, iyy, izz)

! ----------

          rhoremi ( ixx) = AMAX1 ( rhofloor, rhoremi ( ixx))
          a ( ixx) = 1.0 / rhoremi ( ixx)

          vlxremi ( ixx) = vlxremi ( ixx) * a ( ixx)
          vlyremi ( ixx) = vlyremi ( ixx) * a ( ixx)
          vlzremi ( ixx) = vlzremi ( ixx) * a ( ixx)

          b ( ixx) = 0.5 * rhoremi ( ixx) * ( vlxremi ( ixx)**2
     1      + vlyremi ( ixx)**2 + vlzremi ( ixx)**2)

          prsremi ( ixx) = ( prsremi ( ixx) - b ( ixx))
     1                   * ( gamma_eos_remi ( ixx) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropyremi ( ixx) = entropygr ( ixx, iyy, izz, 1)
     1                 + 0.5 * entropygr ( ixx, iyy, izz, 2)
!     1            + dy_shift * entropygr ( ixx, iyy, izz, 3)
!     1            + dz_shift * entropygr ( ixx, iyy, izz, 4)
     1            + dt_shift * entropygr ( ixx, iyy, izz, 5)
          IF ( prsremi ( ixx) .LE. prsfloor) prsremi ( ixx) =
     1    entropyremi ( ixx)
     1  * rhoremi ( ixx)**( gamma_eos_remi ( ixx) - 1.0)
#endif
          prsremi ( ixx) = AMAX1 ( prsfloor, prsremi ( ixx))

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, ixx) =
     1             u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1     + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 2)
     1+ dt_shift * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 5)

          specieremi ( 1: NFLUID_EULER, ixx) =
     1    specieremi ( 1: NFLUID_EULER, ixx) * a ( ixx)
#endif

! ----------

          rhorepl ( ixx) = AMAX1 ( rhofloor, rhorepl ( ixx))
          a ( ixx) = 1.0 / rhorepl ( ixx)

          vlxrepl ( ixx) = vlxrepl ( ixx) * a ( ixx)
          vlyrepl ( ixx) = vlyrepl ( ixx) * a ( ixx)
          vlzrepl ( ixx) = vlzrepl ( ixx) * a ( ixx)

          b ( ixx) = 0.5 * rhorepl ( ixx) * ( vlxrepl ( ixx)**2
     1      + vlyrepl ( ixx)**2 + vlzrepl ( ixx)**2)

          prsrepl ( ixx) = ( prsrepl ( ixx) - b ( ixx))
     1                   * ( gamma_eos_repl ( ixx) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropyrepl ( ixx) = entropygr ( ixx + 1, iyy, izz, 1)
     1                 - 0.5 * entropygr ( ixx + 1, iyy, izz, 2)
!     1            + dy_shift * entropygr ( ixx + 1, iyy, izz, 3)
!     1            + dz_shift * entropygr ( ixx + 1, iyy, izz, 4)
     1            + dt_shift * entropygr ( ixx + 1, iyy, izz, 5)
          IF ( prsrepl ( ixx) .LE. prsfloor) prsrepl ( ixx) =
     1    entropyrepl ( ixx)
     1  * rhorepl ( ixx)**( gamma_eos_repl ( ixx) - 1.0)
#endif
          prsrepl ( ixx) = AMAX1 ( prsfloor, prsrepl ( ixx))

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, ixx) =
     1             u ( ixx + 1, iyy, izz, nspecies_min: nspecies_max, 1)
     1     - 0.5 * u ( ixx + 1, iyy, izz, nspecies_min: nspecies_max, 2)
     1+ dt_shift * u ( ixx + 1, iyy, izz, nspecies_min: nspecies_max, 5)

          specierepl ( 1: NFLUID_EULER, ixx) =
     1    specierepl ( 1: NFLUID_EULER, ixx) * a ( ixx)
#endif

! ----------

          msonic_arr ( ixx) = 0.5 * ( msonicgr ( ixx, iyy, izz)
     1                              + msonicgr ( ixx + 1, iyy, izz))
          divvel_arr ( ixx) = 0.5 * ( divvelgr ( ixx, iyy, izz)
     1                              + divvelgr ( ixx + 1, iyy, izz))

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

! Zero out all other variables that go through RS interface. This is not
! needed but is done to overcome deficiencies in certain compilers.
! They stop allocating memory for automatic arrays if they are not initialized.

          lf_eigenvec_re ( ixx, :, :) = 0.0
          rt_eigenvec_re ( ixx, :, :) = 0.0
          eigenval_re ( ixx, :) = 0.0
          eigenwt_re ( ixx, :) = 0.0

          conserved_re ( ixx, :) = 0.0
          flux_re ( ixx, :) = 0.0

          zonebdy_flux ( ixx, :, :) = 0.0
          zonebdy_conserved ( ixx, :, :) = 0.0

          rho0 ( ixx) = 0.0
          prs0 ( ixx) = 0.0
          vlx0 ( ixx) = 0.0
          vly0 ( ixx) = 0.0
          vlz0 ( ixx) = 0.0
          entropy_flux_re ( ixx) = 0.0

          specie0 ( :, ixx) = 0.0

          do_ein ( ixx) = 0

! ----------

          END DO ! End of DO ixx = loop.

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Continuitive boundary conditions in the x-direction shown below.

          IF ( bcarr_x ( 1, iyy, izz) .EQ. 2) THEN

          rhoremi ( nriembegin + one_x) = rhorepl ( nriembegin + one_x)
          prsremi ( nriembegin + one_x) = prsrepl ( nriembegin + one_x)

          vlxremi ( nriembegin + one_x) = vlxrepl ( nriembegin + one_x)
          vlyremi ( nriembegin + one_x) = vlyrepl ( nriembegin + one_x)
          vlzremi ( nriembegin + one_x) = vlzrepl ( nriembegin + one_x)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_x) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_x)
#endif

          END IF

          IF ( bcarr_x ( 2, iyy, izz) .EQ. 2) THEN

          rhorepl ( nriemend - one_x) = rhoremi ( nriemend - one_x)
          prsrepl ( nriemend - one_x) = prsremi ( nriemend - one_x)

          vlxrepl ( nriemend - one_x) = vlxremi ( nriemend - one_x)
          vlyrepl ( nriemend - one_x) = vlyremi ( nriemend - one_x)
          vlzrepl ( nriemend - one_x) = vlzremi ( nriemend - one_x)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_x) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_x)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Reflection boundary conditions in the x-direction shown below.

          IF ( bcarr_x ( 1, iyy, izz) .EQ. 3) THEN

          rhoremi ( nriembegin + one_x) = rhorepl ( nriembegin + one_x)
          prsremi ( nriembegin + one_x) = prsrepl ( nriembegin + one_x)

          vlxremi ( nriembegin + one_x)= -vlxrepl ( nriembegin + one_x)
          vlyremi ( nriembegin + one_x) = vlyrepl ( nriembegin + one_x)
          vlzremi ( nriembegin + one_x) = vlzrepl ( nriembegin + one_x)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_x) = 
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_x)
#endif

          END IF

          IF ( bcarr_x ( 2, iyy, izz) .EQ. 3) THEN

          rhorepl ( nriemend - one_x) = rhoremi ( nriemend - one_x)
          prsrepl ( nriemend - one_x) = prsremi ( nriemend - one_x)

          vlxrepl ( nriemend - one_x)= -vlxremi ( nriemend - one_x)
          vlyrepl ( nriemend - one_x) = vlyremi ( nriemend - one_x)
          vlzrepl ( nriemend - one_x) = vlzremi ( nriemend - one_x)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_x) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_x)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Periodic boundary conditions in the x-direction shown below.Works for serial.

          IF ( bcarr_x ( 1, iyy, izz) .EQ. 4) THEN

          rhoremi ( nriembegin + one_x) = rhoremi ( nriemend - one_x)
          prsremi ( nriembegin + one_x) = prsremi ( nriemend - one_x)

          vlxremi ( nriembegin + one_x) = vlxremi ( nriemend - one_x)
          vlyremi ( nriembegin + one_x) = vlyremi ( nriemend - one_x)
          vlzremi ( nriembegin + one_x) = vlzremi ( nriemend - one_x)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_x) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_x)
#endif

          END IF

          IF ( bcarr_x ( 2, iyy, izz) .EQ. 4) THEN

          rhorepl ( nriemend - one_x) = rhorepl ( nriembegin + one_x)
          prsrepl ( nriemend - one_x) = prsrepl ( nriembegin + one_x)

          vlxrepl ( nriemend - one_x) = vlxrepl ( nriembegin + one_x)
          vlyrepl ( nriemend - one_x) = vlyrepl ( nriembegin + one_x)
          vlzrepl ( nriemend - one_x) = vlzrepl ( nriembegin + one_x)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_x) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_x)
#endif

          END IF

! ----------------------------------------------------------------------

! Call the Riemann solver. HLL is the default.

#if ( USE_RIEMANN_SOLVER == 1)
      CALL EULER_RIEM_ROE
#elif ( USE_RIEMANN_SOLVER == 2)
      CALL EULER_RIEM_HLLC
#elif ( USE_RIEMANN_SOLVER == 3)
      CALL EULER_RIEM_HLL
#else
      CALL EULER_RIEM_HLLEM_DB
#endif

     1   ( nriembegin, nriemend,

     1     gamma_euler, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst_euler,
     1  cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlx0, vly0, vlz0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

! Add in the contributions from the lapidus viscosity to the 1d fluxes.

          IF ( lapidus_coef .GT. 0.0) THEN

          DO ixx = nriembegin, nriemend

          flux_re ( ixx, 1) = flux_re ( ixx, 1) + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, 1, 1) - u ( ixx + 1, iyy, izz, 1, 1) )

          flux_re ( ixx, 2) = flux_re ( ixx, 2) + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, 2, 1) - u ( ixx + 1, iyy, izz, 2, 1) )
          flux_re ( ixx, 3) = flux_re ( ixx, 3) + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, 3, 1) - u ( ixx + 1, iyy, izz, 3, 1) )
          flux_re ( ixx, 4) = flux_re ( ixx, 4) + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, 4, 1) - u ( ixx + 1, iyy, izz, 4, 1) )

          flux_re ( ixx, 5) = flux_re ( ixx, 5) + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, 5, 1) - u ( ixx + 1, iyy, izz, 5, 1) )

#if ( NFLUID_EULER >= 1)
          flux_re ( ixx, nspecies_min + 1: nspecies_max + 1) =
     1    flux_re ( ixx, nspecies_min + 1: nspecies_max + 1)
     1  + lapidus_arr ( ixx)
     1  * ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1    - u ( ixx + 1, iyy, izz, nspecies_min: nspecies_max, 1) )
#endif

          END DO

          END IF

! ----------------------------------------------------------------------

! Hand back fluxes. Fluxes have been evaluated at face centers.

          DO ixx = nriembegin, nriemend

          flux_x ( ixx, iyy, izz, 1) = flux_re ( ixx, 1)
          flux_x ( ixx, iyy, izz, 2) = flux_re ( ixx, 2)
          flux_x ( ixx, iyy, izz, 3) = flux_re ( ixx, 3)
          flux_x ( ixx, iyy, izz, 4) = flux_re ( ixx, 4)
          flux_x ( ixx, iyy, izz, 5) = flux_re ( ixx, 5)

          prs_flux_x ( ixx, iyy, izz) = flux_re ( ixx, 6)

#if ( NFLUID_EULER >= 1)
          flux_x ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1            flux_re ( ixx, nspecies_min + 1: nspecies_max + 1)
#endif

#if ( PRESSURE_POSITIVITY == 1)
          entropy_flux_x ( ixx, iyy, izz) = entropy_flux_re ( ixx)
#endif

          END DO ! End of DO ixx = loop.

! ----------

! If we are doing the predictor stage then hand back the scaled time
! rate of update to the zone-centered "u ( ixx, iyy, izz, :, 5)".
! Note shift to zone-centered index.

          IF ( ABS ( dt_shift) .LT. smallnum) THEN

! ----------
! Cartesian

          IF ( igeom .EQ. 1) THEN

          DO ixx = nriembegin + 1, nriemend

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

          END DO ! End of DO ixx = loop.

          END IF

! -----
! Cylindrical

          IF ( igeom .EQ. 2) THEN

          DO ixx = nriembegin + 1, nriemend

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

          END DO ! End of DO ixx = loop.

          END IF

! -----
! Spherical

          IF ( igeom .EQ. 3) THEN

          DO ixx = nriembegin + 1, nriemend

          a ( ixx) = ( xb_save ( ixx)**3
     1               - xb_save ( ixx - 1)**3 ) / 3.0
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

          END DO ! End of DO ixx = loop.

          END IF

! ----------

          END IF ! End of "dt_shift .EQ. 0.0"

! ----------------------------------------------------------------------

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_EULER_FLUX_X

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_EULER_FLUX_Y

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, flux_y, prs_flux_y,
     1        entropygr, entropy_flux_y,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        y_indx_limits, bcarr_y,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and puts out "flux_y and prs_flux_y". 

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_y".

! Since this is an ADER-style update, it is suggested that the update 
! be done in two stages:

! Stage I : With the slopes evaluated, enter with
! "dt_shift = 0.0". This gives 
! fluxes that are used to obtain the time rate of update 
! "u ( ixx, iyy, izz, :, 5)" within a zone.

! Stage II : With the slopes and the time rate of update
! evaluated, enter with "dt_shift = 0.5".
! This reloads the flux with time-centered
! values, suitable for a spatially and temporally second order accurate update.

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
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
     1                  n_cc_components) :: flux_y

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1     prs_flux_y, entropy_flux_y

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      INTEGER ixx, iyy, izz, nriembegin, nriemend, one_y,
     1        iprimitive, use_approx_eos, split_prsflux_re, debug_stop,
     1        nspecies_min, nspecies_max

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b, c

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,
     1     entropyremi, entropyrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr, lapidus_arr

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD_EULER + NFLUID_EULER, NFIELD_EULER + NFLUID_EULER) ::
     1    lf_eigenvec_re, rt_eigenvec_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD_EULER + NFLUID_EULER) ::
     1     eigenval_re, eigenwt_re, conserved_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD_EULER + NFLUID_EULER + 1) ::
     1     flux_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD_EULER + NFLUID_EULER) ::
     1      zonebdy_flux, zonebdy_conserved

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------

! The computation is not dynamically active in the y-direction
! ( as indicated by setting "ioffy = 0" on the outside) so return.

      IF ( ioffy .EQ. 0) RETURN

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffy .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iy1 + ioffy .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_EULER_FLUX_Y : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffy, present_base = ",
     1              0 - ioffy, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iy1 + ioffy, present_bound = ",
     1              iy1 + ioffy, ONEDFLUIDEND
      STOP
      END IF

      one_y = MIN0 ( expand_by_one_zone, ioffy)

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

! ----------

      iprimitive = 0     ! Eigenstructure constructed in conserved variables.
      use_approx_eos = 1 ! Just use simple EOS with constant polytropic index.

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c, debug_stop,         &
!$OMP&    nriembegin, nriemend,                                         &
!$OMP&    rhoremi, prsremi, vlxremi, vlyremi, vlzremi, specieremi,      &
!$OMP&    mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,      &
!$OMP&    entropyremi, msonic_arr, divvel_arr, lapidus_arr,             &
!$OMP&    rhorepl, prsrepl, vlxrepl, vlyrepl, vlzrepl, specierepl,      &
!$OMP&    mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,      &
!$OMP&    entropyrepl, zonebdy_conserved, zonebdy_flux,                 &
!$OMP&    eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,      &
!$OMP&    conserved_re, flux_re, entropy_flux_re,                       &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, specie0, do_ein)                &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, one_y, expand_by_one_zone,                             &
!$OMP&    dtcur, dt_shift, nspecies_min, nspecies_max,                  &
!$OMP&    gamma_euler, gasconst_euler, smallnum, rhofloor, prsfloor,    &
!$OMP&    ientropyfix, iprimitive, use_approx_eos, fracdiff,            &
!$OMP&    split_prsflux_re, divratio_einfeldt_rs,                       &
!$OMP&    lapidus_threshold, lapidus_coef,                              &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    y_indx_limits, bcarr_y, u, flux_y, prs_flux_y, entropygr,     &
!$OMP&    entropy_flux_y, mean_mol_wt_gr, gamma_eos_gr,                 &
!$OMP&    gamma_soundspeed_gr, msonicgr, divvelgr)                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( expand_by_one_zone, ioffz),
     1       iz1 + MIN0 ( expand_by_one_zone, ioffz)

        DO ixx = 1 - MIN0 ( expand_by_one_zone, ioffx),
     1         ix1 + MIN0 ( expand_by_one_zone, ioffx)

        debug_stop = 0

! ----------------------------------------------------------------------

! Load in the left and right states for the Riemann solver.

          nriembegin = y_indx_limits ( ixx, 1, izz) - 1 - one_y
          nriemend = y_indx_limits ( ixx, 2, izz) + one_y

          DO iyy = nriembegin, nriemend

! ----------

          rhoremi ( iyy) = u ( ixx, iyy, izz, 1, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 1, 3)
     1        + dt_shift * u ( ixx, iyy, izz, 1, 5)

          vlxremi ( iyy) = u ( ixx, iyy, izz, 2, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 2, 3)
     1        + dt_shift * u ( ixx, iyy, izz, 2, 5)

          vlyremi ( iyy) = u ( ixx, iyy, izz, 3, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 3, 3)
     1        + dt_shift * u ( ixx, iyy, izz, 3, 5)

          vlzremi ( iyy) = u ( ixx, iyy, izz, 4, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 4, 3)
     1        + dt_shift * u ( ixx, iyy, izz, 4, 5)

          prsremi ( iyy) = u ( ixx, iyy, izz, 5, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 5, 3)
     1        + dt_shift * u ( ixx, iyy, izz, 5, 5)

          mean_mol_wt_remi ( iyy) = mean_mol_wt_gr ( ixx, iyy, izz)
          gamma_eos_remi ( iyy) = gamma_eos_gr ( ixx, iyy, izz)
          gamma_soundspeed_remi ( iyy) =
     1                         gamma_soundspeed_gr ( ixx, iyy, izz) 

! ----------

          rhorepl ( iyy) = u ( ixx, iyy + 1, izz, 1, 1)
     1             - 0.5 * u ( ixx, iyy + 1, izz, 1, 3)
     1        + dt_shift * u ( ixx, iyy + 1, izz, 1, 5)

          vlxrepl ( iyy) = u ( ixx, iyy + 1, izz, 2, 1)
     1             - 0.5 * u ( ixx, iyy + 1, izz, 2, 3)
     1        + dt_shift * u ( ixx, iyy + 1, izz, 2, 5)

          vlyrepl ( iyy) = u ( ixx, iyy + 1, izz, 3, 1)
     1             - 0.5 * u ( ixx, iyy + 1, izz, 3, 3)
     1        + dt_shift * u ( ixx, iyy + 1, izz, 3, 5)

          vlzrepl ( iyy) = u ( ixx, iyy + 1, izz, 4, 1)
     1             - 0.5 * u ( ixx, iyy + 1, izz, 4, 3)
     1        + dt_shift * u ( ixx, iyy + 1, izz, 4, 5)

          prsrepl ( iyy) = u ( ixx, iyy + 1, izz, 5, 1)
     1             - 0.5 * u ( ixx, iyy + 1, izz, 5, 3)
     1        + dt_shift * u ( ixx, iyy + 1, izz, 5, 5)

          mean_mol_wt_repl ( iyy) = mean_mol_wt_gr ( ixx, iyy + 1, izz)
          gamma_eos_repl ( iyy) = gamma_eos_gr ( ixx, iyy + 1, izz)
          gamma_soundspeed_repl ( iyy) =
     1                         gamma_soundspeed_gr ( ixx, iyy + 1, izz)

! ----------

          rhoremi ( iyy) = AMAX1 ( rhofloor, rhoremi ( iyy))
          a ( iyy) = 1.0 / rhoremi ( iyy)

          vlxremi ( iyy) = vlxremi ( iyy) * a ( iyy)
          vlyremi ( iyy) = vlyremi ( iyy) * a ( iyy)
          vlzremi ( iyy) = vlzremi ( iyy) * a ( iyy)

          b ( iyy) = 0.5 * rhoremi ( iyy) * ( vlxremi ( iyy)**2
     1      + vlyremi ( iyy)**2 + vlzremi ( iyy)**2)

          prsremi ( iyy) = ( prsremi ( iyy) - b ( iyy))
     1                   * ( gamma_eos_remi ( iyy) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropyrepl ( iyy) = entropygr ( ixx, iyy + 1, izz, 1)
!     1            + dx_shift * entropygr ( ixx, iyy + 1, izz, 2)
     1                 - 0.5 * entropygr ( ixx, iyy + 1, izz, 3)
!     1            + dz_shift * entropygr ( ixx, iyy + 1, izz, 4)
     1            + dt_shift * entropygr ( ixx, iyy + 1, izz, 5)
          IF ( prsrepl ( iyy) .LE. prsfloor) prsrepl ( iyy) =
     1    entropyrepl ( iyy)
     1  * rhorepl ( iyy)**( gamma_eos_repl ( iyy) - 1.0)
#endif
          prsremi ( iyy) = AMAX1 ( prsfloor, prsremi ( iyy))

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, iyy) =
     1             u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1     + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 3)
     1+ dt_shift * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 5)

          specieremi ( 1: NFLUID_EULER, iyy) =
     1    specieremi ( 1: NFLUID_EULER, iyy) * a ( iyy)
#endif

! ----------

          rhorepl ( iyy) = AMAX1 ( rhofloor, rhorepl ( iyy))
          a ( iyy) = 1.0 / rhorepl ( iyy)

          vlxrepl ( iyy) = vlxrepl ( iyy) * a ( iyy)
          vlyrepl ( iyy) = vlyrepl ( iyy) * a ( iyy)
          vlzrepl ( iyy) = vlzrepl ( iyy) * a ( iyy)

          b ( iyy) = 0.5 * rhorepl ( iyy) * ( vlxrepl ( iyy)**2
     1      + vlyrepl ( iyy)**2 + vlzrepl ( iyy)**2)

          prsrepl ( iyy) = ( prsrepl ( iyy) - b ( iyy))
     1                   * ( gamma_eos_repl ( iyy) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropyrepl ( iyy) = entropygr ( ixx, iyy + 1, izz, 1)
!     1            + dx_shift * entropygr ( ixx, iyy + 1, izz, 2)
     1                 - 0.5 * entropygr ( ixx, iyy + 1, izz, 3)
!     1            + dz_shift * entropygr ( ixx, iyy + 1, izz, 4)
     1            + dt_shift * entropygr ( ixx, iyy + 1, izz, 5)
          IF ( prsrepl ( iyy) .LE. prsfloor) prsrepl ( iyy) =
     1    entropyrepl ( iyy)
     1  * rhorepl ( iyy)**( gamma_eos_repl ( iyy) - 1.0)
#endif
          prsrepl ( iyy) = AMAX1 ( prsfloor, prsrepl ( iyy))

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, iyy) =
     1             u ( ixx, iyy + 1, izz, nspecies_min: nspecies_max, 1)
     1     - 0.5 * u ( ixx, iyy + 1, izz, nspecies_min: nspecies_max, 3)
     1+ dt_shift * u ( ixx, iyy + 1, izz, nspecies_min: nspecies_max, 5)

          specierepl ( 1: NFLUID_EULER, iyy) =
     1    specierepl ( 1: NFLUID_EULER, iyy) * a ( iyy)
#endif

! ----------

          msonic_arr ( iyy) = 0.5 * ( msonicgr ( ixx, iyy, izz)
     1                              + msonicgr ( ixx, iyy + 1, izz))
          divvel_arr ( iyy) = 0.5 * ( divvelgr ( ixx, iyy, izz)
     1                              + divvelgr ( ixx, iyy + 1, izz))

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

! Zero out all other variables that go through RS interface. This is not
! needed but is done to overcome deficiencies in certain compilers.
! They stop allocating memory for automatic arrays if they are not initialized.

          lf_eigenvec_re ( iyy, :, :) = 0.0
          rt_eigenvec_re ( iyy, :, :) = 0.0
          eigenval_re ( iyy, :) = 0.0
          eigenwt_re ( iyy, :) = 0.0

          conserved_re ( iyy, :) = 0.0
          flux_re ( iyy, :) = 0.0

          zonebdy_flux ( iyy, :, :) = 0.0
          zonebdy_conserved ( iyy, :, :) = 0.0

          rho0 ( iyy) = 0.0
          prs0 ( iyy) = 0.0
          vlx0 ( iyy) = 0.0
          vly0 ( iyy) = 0.0
          vlz0 ( iyy) = 0.0
          entropy_flux_re ( iyy) = 0.0

          specie0 ( :, iyy) = 0.0

          do_ein ( iyy) = 0

! ----------

          END DO ! End of DO iyy = loop.

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Continuitive boundary conditions in the y-direction shown below.

          IF ( bcarr_y ( ixx, 1, izz) .EQ. 2) THEN

          rhoremi ( nriembegin + one_y) = rhorepl ( nriembegin + one_y)
          prsremi ( nriembegin + one_y) = prsrepl ( nriembegin + one_y)

          vlxremi ( nriembegin + one_y) = vlxrepl ( nriembegin + one_y)
          vlyremi ( nriembegin + one_y) = vlyrepl ( nriembegin + one_y)
          vlzremi ( nriembegin + one_y) = vlzrepl ( nriembegin + one_y)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_y) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_y)
#endif

          END IF

          IF ( bcarr_y ( ixx, 2, izz) .EQ. 2) THEN

          rhorepl ( nriemend - one_y) = rhoremi ( nriemend - one_y)
          prsrepl ( nriemend - one_y) = prsremi ( nriemend - one_y)

          vlxrepl ( nriemend - one_y) = vlxremi ( nriemend - one_y)
          vlyrepl ( nriemend - one_y) = vlyremi ( nriemend - one_y)
          vlzrepl ( nriemend - one_y) = vlzremi ( nriemend - one_y)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_y) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_y)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Reflection boundary conditions in the y-direction shown below.

          IF ( bcarr_y ( ixx, 1, izz) .EQ. 3) THEN

          rhoremi ( nriembegin + one_y) = rhorepl ( nriembegin + one_y)
          prsremi ( nriembegin + one_y) = prsrepl ( nriembegin + one_y)

          vlxremi ( nriembegin + one_y) = vlxrepl ( nriembegin + one_y)
          vlyremi ( nriembegin + one_y)= -vlyrepl ( nriembegin + one_y)
          vlzremi ( nriembegin + one_y) = vlzrepl ( nriembegin + one_y)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_y) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_y)
#endif

          END IF

          IF ( bcarr_y ( ixx, 2, izz) .EQ. 3) THEN

          rhorepl ( nriemend - one_y) = rhoremi ( nriemend - one_y)
          prsrepl ( nriemend - one_y) = prsremi ( nriemend - one_y)

          vlxrepl ( nriemend - one_y) = vlxremi ( nriemend - one_y)
          vlyrepl ( nriemend - one_y)= -vlyremi ( nriemend - one_y)
          vlzrepl ( nriemend - one_y) = vlzremi ( nriemend - one_y)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_y) = 
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_y)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Periodic boundary conditions in the y-direction shown below. Works in serial.

          IF ( bcarr_y ( ixx, 1, izz) .EQ. 4) THEN

          rhoremi ( nriembegin + one_y) = rhoremi ( nriemend - one_y)
          prsremi ( nriembegin + one_y) = prsremi ( nriemend - one_y)

          vlxremi ( nriembegin + one_y) = vlxremi ( nriemend - one_y)
          vlyremi ( nriembegin + one_y) = vlyremi ( nriemend - one_y)
          vlzremi ( nriembegin + one_y) = vlzremi ( nriemend - one_y)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_y) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_y)
#endif

          END IF

          IF ( bcarr_y ( ixx, 2, izz) .EQ. 4) THEN

          rhorepl ( nriemend - one_y) = rhorepl ( nriembegin + one_y)
          prsrepl ( nriemend - one_y) = prsrepl ( nriembegin + one_y)

          vlxrepl ( nriemend - one_y) = vlxrepl ( nriembegin + one_y)
          vlyrepl ( nriemend - one_y) = vlyrepl ( nriembegin + one_y)
          vlzrepl ( nriemend - one_y) = vlzrepl ( nriembegin + one_y)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_y) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_y)
#endif

          END IF

! ----------------------------------------------------------------------

! Call the Riemann solver with cyclic rotation of input variables.
! HLL is the default.

#if ( USE_RIEMANN_SOLVER == 1)
      CALL EULER_RIEM_ROE
#elif ( USE_RIEMANN_SOLVER == 2)
      CALL EULER_RIEM_HLLC
#elif ( USE_RIEMANN_SOLVER == 3)
      CALL EULER_RIEM_HLL
#else
      CALL EULER_RIEM_HLLEM_DB
#endif

     1   ( nriembegin, nriemend,

     1     gamma_euler, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst_euler,
     1  cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlyremi, vlyrepl,
     1        vlzremi, vlzrepl,
     1        vlxremi, vlxrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vly0, vlz0, vlx0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

! Add in the contributions from the lapidus viscosity to the 1d fluxes.

          IF ( lapidus_coef .GT. 0.0) THEN

          DO iyy = nriembegin, nriemend

          flux_re ( iyy, 1) = flux_re ( iyy, 1) + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, 1, 1) - u ( ixx, iyy + 1, izz, 1, 1) )

          flux_re ( iyy, 4) = flux_re ( iyy, 4) + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, 2, 1) - u ( ixx, iyy + 1, izz, 2, 1) )
          flux_re ( iyy, 2) = flux_re ( iyy, 2) + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, 3, 1) - u ( ixx, iyy + 1, izz, 3, 1) )
          flux_re ( iyy, 3) = flux_re ( iyy, 3) + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, 4, 1) - u ( ixx, iyy + 1, izz, 4, 1) )

          flux_re ( iyy, 5) = flux_re ( iyy, 5) + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, 5, 1) - u ( ixx, iyy + 1, izz, 5, 1) )

#if ( NFLUID_EULER >= 1)
          flux_re ( iyy, nspecies_min + 1: nspecies_max + 1) =
     1    flux_re ( iyy, nspecies_min + 1: nspecies_max + 1)
     1  + lapidus_arr ( iyy)
     1  * ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1    - u ( ixx, iyy + 1, izz, nspecies_min: nspecies_max, 1) )
#endif

          END DO

          END IF

! ----------------------------------------------------------------------

! Hand back fluxes. Use cyclic rotation for fluxes.
! Fluxes have been evaluated at face centers.

          DO iyy = nriembegin, nriemend

          flux_y ( ixx, iyy, izz, 1) = flux_re ( iyy, 1)
          flux_y ( ixx, iyy, izz, 2) = flux_re ( iyy, 4)
          flux_y ( ixx, iyy, izz, 3) = flux_re ( iyy, 2)
          flux_y ( ixx, iyy, izz, 4) = flux_re ( iyy, 3)
          flux_y ( ixx, iyy, izz, 5) = flux_re ( iyy, 5)

          prs_flux_y ( ixx, iyy, izz) = flux_re ( iyy, 6)

#if ( NFLUID_EULER >= 1)
          flux_y ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1            flux_re ( iyy, nspecies_min + 1: nspecies_max + 1)
#endif

#if ( PRESSURE_POSITIVITY == 1)
          entropy_flux_y ( ixx, iyy, izz) = entropy_flux_re ( iyy)
#endif

          END DO ! End of DO iyy = loop.

! ----------

! If we are doing the predictor stage then hand back the scaled time
! rate of update to the zone-centered "u ( ixx, iyy, izz, :, 5)".
! Note shift to zone-centered index.

          IF ( ABS ( dt_shift) .LT. smallnum) THEN

! ----------
! Cartesian

          IF ( igeom .EQ. 1) THEN

          DO iyy = nriembegin + 1, nriemend

          a ( iyy) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( iyy)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( iyy)
#endif

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( iyy)

          END DO ! End of DO iyy = loop.

          END IF

! -----
! Cylindrical

          IF ( igeom .EQ. 2) THEN

          DO iyy = nriembegin + 1, nriemend

          a ( iyy) = xc_save ( ixx) * dy_save ( iyy)
          a ( iyy) = dtcur / a ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( iyy)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( iyy)
#endif

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( iyy)

          END DO ! End of DO iyy = loop.

          END IF

! -----
! Spherical

          IF ( igeom .EQ. 3) THEN

          DO iyy = nriembegin + 1, nriemend

          a ( iyy) = - xc_save ( ixx) * ( COS ( yb_save ( iyy))
     1                           - COS ( yb_save ( iyy - 1)) )
          a ( iyy) = dtcur / a ( iyy)

          b ( iyy) = SIN ( yb_save ( iyy))

          c ( iyy) = SIN ( yb_save ( iyy - 1))

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) * b ( iyy)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) * c ( iyy) )
     1  * a ( iyy)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_y ( ixx, iyy, izz) * b ( iyy)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) * c ( iyy) ) * a ( iyy)
#endif

          a ( iyy) = xc_save ( ixx) * dy_save ( iyy)
          a ( iyy) = dtcur / a ( iyy)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( iyy)

          END DO ! End of DO iyy = loop.

          END IF

! ----------

          END IF ! End of "dt_shift .EQ. 0.0"

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_EULER_FLUX_Y

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MAKE_EULER_FLUX_Z

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone,
     1        dtcur, dt_shift,

     1        u, flux_z, prs_flux_z,
     1        entropygr, entropy_flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        z_indx_limits, bcarr_z,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and puts out "flux_z and prs_flux_z".

! If pressure positivity is called for, it also takes in "entropygr" and
! puts out "entropy_flux_z".

! Since this is an ADER-style update, it is suggested that the update 
! be done in two stages:

! Stage I : With the slopes evaluated, enter with
! "dt_shift = 0.0". This gives 
! fluxes that are used to obtain the time rate of update 
! "u ( ixx, iyy, izz, :, 5)" within a zone.

! Stage II : With the slopes and the time rate of update
! evaluated, enter with "dt_shift = 0.5".
! This reloads the flux with time-centered
! values, suitable for a spatially and temporally second order accurate update.

! "expand_by_one_zone = 1" allows one to expand the dynamics by one
! zone outside the physical domain. This is useful for periodic
! domains in serial setting as well as parallel computations.

! "dtcur" is the timestep.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
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
     1                  n_cc_components) :: flux_z

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1      prs_flux_z, entropy_flux_z

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      INTEGER ixx, iyy, izz, nriembegin, nriemend, one_z,
     1        iprimitive, use_approx_eos, split_prsflux_re, debug_stop,
     1        nspecies_min, nspecies_max

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: a, b

! ----------------------------------------------------------------------

! These are the arrays that go through the Riemann solver interface.

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rhoremi, rhorepl, prsremi, prsrepl,
     1     entropyremi, entropyrepl,

     1     vlxremi, vlxrepl, vlyremi, vlyrepl,
     1     vlzremi, vlzrepl,

     1     mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,
     1     mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,

     1     msonic_arr, divvel_arr, lapidus_arr

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     specieremi, specierepl

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1    NFIELD_EULER + NFLUID_EULER, NFIELD_EULER + NFLUID_EULER) ::
     1    lf_eigenvec_re, rt_eigenvec_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1      NFIELD_EULER + NFLUID_EULER) ::
     1     eigenval_re, eigenwt_re, conserved_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                 NFIELD_EULER + NFLUID_EULER + 1) ::
     1     flux_re

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND,
     1                  0 : 1, NFIELD_EULER + NFLUID_EULER) ::
     1      zonebdy_flux, zonebdy_conserved

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1     rho0, prs0, vlx0, vly0, vlz0,
     1     entropy_flux_re

      REAL, DIMENSION ( 0: NFLUID_EULER,
     1                  ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1      specie0

      INTEGER, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) :: do_ein

! ----------------------------------------------------------------------

! The computation is not dynamically active in the z-direction
! ( as indicated by setting "ioffz = 0" on the outside) so return.

      IF ( ioffz .EQ. 0) RETURN

! ----------

! Check limits of the domain against the length of 1d static arrays.

      IF ( ( 0 - ioffz .LT. ONEDFLUIDBEGIN) .OR.
     1     ( iz1 + ioffz .GT. ONEDFLUIDEND) ) THEN
      WRITE ( 6, *)"Stopping in MAKE_EULER_FLUX_Z : expand 1d arrays"
      WRITE ( 6, *)"0 - ioffz, present_base = ",
     1              0 - ioffz, ONEDFLUIDBEGIN
      WRITE ( 6, *)"iz1 + ioffz, present_bound = ",
     1              iz1 + ioffz, ONEDFLUIDEND
      STOP
      END IF

      one_z = MIN0 ( expand_by_one_zone, ioffz)

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

! ----------

      iprimitive = 0     ! Eigenstructure constructed in conserved variables.
      use_approx_eos = 1 ! Just use simple EOS with constant polytropic index.

      IF ( igeom .EQ. 1) THEN
      split_prsflux_re = 0  ! Don't want to split the momentum flux into prs+..
      ELSE
      split_prsflux_re = 1
      END IF

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, debug_stop,            &
!$OMP&    nriembegin, nriemend,                                         &
!$OMP&    rhoremi, prsremi, vlxremi, vlyremi, vlzremi, specieremi,      &
!$OMP&    mean_mol_wt_remi, gamma_eos_remi, gamma_soundspeed_remi,      &
!$OMP&    entropyremi, msonic_arr, divvel_arr, lapidus_arr,             &
!$OMP&    rhorepl, prsrepl, vlxrepl, vlyrepl, vlzrepl, specierepl,      &
!$OMP&    mean_mol_wt_repl, gamma_eos_repl, gamma_soundspeed_repl,      &
!$OMP&    entropyrepl, zonebdy_conserved, zonebdy_flux,                 &
!$OMP&    eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,      &
!$OMP&    conserved_re, flux_re, entropy_flux_re,                       &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, specie0, do_ein)                &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    igeom, one_z, expand_by_one_zone,                             &
!$OMP&    dtcur, dt_shift, nspecies_min, nspecies_max,                  &
!$OMP&    gamma_euler, gasconst_euler, smallnum, rhofloor, prsfloor,    &
!$OMP&    ientropyfix, iprimitive, use_approx_eos, fracdiff,            &
!$OMP&    split_prsflux_re, divratio_einfeldt_rs,                       &
!$OMP&    lapidus_threshold, lapidus_coef,                              &
!$OMP&    cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,&
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    z_indx_limits, bcarr_z, u, flux_z, prs_flux_z, entropygr,     &
!$OMP&    entropy_flux_z, mean_mol_wt_gr, gamma_eos_gr,                 &
!$OMP&    gamma_soundspeed_gr, msonicgr, divvelgr)                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 1 - MIN0 ( expand_by_one_zone, ioffy),
     1       iy1 + MIN0 ( expand_by_one_zone, ioffy)

        DO ixx = 1 - MIN0 ( expand_by_one_zone, ioffx),
     1         ix1 + MIN0 ( expand_by_one_zone, ioffx)

        debug_stop = 0

! ----------------------------------------------------------------------

! Load in the left and right states for the Riemann solver.

          nriembegin = z_indx_limits ( ixx, iyy, 1) - 1 - one_z
          nriemend = z_indx_limits ( ixx, iyy, 2) + one_z

          DO izz = nriembegin, nriemend

! ----------

          rhoremi ( izz) = u ( ixx, iyy, izz, 1, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 1, 4)
     1        + dt_shift * u ( ixx, iyy, izz, 1, 5)

          vlxremi ( izz) = u ( ixx, iyy, izz, 2, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 2, 4)
     1        + dt_shift * u ( ixx, iyy, izz, 2, 5)

          vlyremi ( izz) = u ( ixx, iyy, izz, 3, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 3, 4)
     1        + dt_shift * u ( ixx, iyy, izz, 3, 5)

          vlzremi ( izz) = u ( ixx, iyy, izz, 4, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 4, 4)
     1        + dt_shift * u ( ixx, iyy, izz, 4, 5)

          prsremi ( izz) = u ( ixx, iyy, izz, 5, 1)
     1             + 0.5 * u ( ixx, iyy, izz, 5, 4)
     1        + dt_shift * u ( ixx, iyy, izz, 5, 5)

          mean_mol_wt_remi ( izz) = mean_mol_wt_gr ( ixx, iyy, izz)
          gamma_eos_remi ( izz) = gamma_eos_gr ( ixx, iyy, izz)
          gamma_soundspeed_remi ( izz) =
     1                         gamma_soundspeed_gr ( ixx, iyy, izz) 

! ----------

          rhorepl ( izz) = u ( ixx, iyy, izz + 1, 1, 1)
     1             - 0.5 * u ( ixx, iyy, izz + 1, 1, 4)
     1        + dt_shift * u ( ixx, iyy, izz + 1, 1, 5)

          vlxrepl ( izz) = u ( ixx, iyy, izz + 1, 2, 1)
     1             - 0.5 * u ( ixx, iyy, izz + 1, 2, 4)
     1        + dt_shift * u ( ixx, iyy, izz + 1, 2, 5)

          vlyrepl ( izz) = u ( ixx, iyy, izz + 1, 3, 1)
     1             - 0.5 * u ( ixx, iyy, izz + 1, 3, 4)
     1        + dt_shift * u ( ixx, iyy, izz + 1, 3, 5)

          vlzrepl ( izz) = u ( ixx, iyy, izz + 1, 4, 1)
     1             - 0.5 * u ( ixx, iyy, izz + 1, 4, 4)
     1        + dt_shift * u ( ixx, iyy, izz + 1, 4, 5)

          prsrepl ( izz) = u ( ixx, iyy, izz + 1, 5, 1)
     1             - 0.5 * u ( ixx, iyy, izz + 1, 5, 4)
     1        + dt_shift * u ( ixx, iyy, izz + 1, 5, 5)

          mean_mol_wt_repl ( izz) = mean_mol_wt_gr ( ixx, iyy, izz + 1)
          gamma_eos_repl ( izz) = gamma_eos_gr ( ixx, iyy, izz + 1)
          gamma_soundspeed_repl ( izz) =
     1                         gamma_soundspeed_gr ( ixx, iyy, izz + 1) 

! ----------

          rhoremi ( izz) = AMAX1 ( rhofloor, rhoremi ( izz))
          a ( izz) = 1.0 / rhoremi ( izz)

          vlxremi ( izz) = vlxremi ( izz) * a ( izz)
          vlyremi ( izz) = vlyremi ( izz) * a ( izz)
          vlzremi ( izz) = vlzremi ( izz) * a ( izz)

          b ( izz) = 0.5 * rhoremi ( izz) * ( vlxremi ( izz)**2
     1      + vlyremi ( izz)**2 + vlzremi ( izz)**2)

          prsremi ( izz) = ( prsremi ( izz) - b ( izz))
     1                   * ( gamma_eos_remi ( izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropyremi ( izz) = entropygr ( ixx, iyy, izz, 1)
!     1            + dx_shift * entropygr ( ixx, iyy, izz, 2)
!     1            + dy_shift * entropygr ( ixx, iyy, izz, 3)
     1                 + 0.5 * entropygr ( ixx, iyy, izz, 4)
     1            + dt_shift * entropygr ( ixx, iyy, izz, 5)
          IF ( prsremi ( izz) .LE. prsfloor) prsremi ( izz) =
     1    entropyremi ( izz)
     1  * rhoremi ( izz)**( gamma_eos_remi ( izz) - 1.0)
#endif
          prsremi ( izz) = AMAX1 ( prsfloor, prsremi ( izz))

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, izz) =
     1             u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1     + 0.5 * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 4)
     1+ dt_shift * u ( ixx, iyy, izz, nspecies_min: nspecies_max, 5)

          specieremi ( 1: NFLUID_EULER, izz) =
     1    specieremi ( 1: NFLUID_EULER, izz) * a ( izz)
#endif

! ----------

          rhorepl ( izz) = AMAX1 ( rhofloor, rhorepl ( izz))
          a ( izz) = 1.0 / rhorepl ( izz)

          vlxrepl ( izz) = vlxrepl ( izz) * a ( izz)
          vlyrepl ( izz) = vlyrepl ( izz) * a ( izz)
          vlzrepl ( izz) = vlzrepl ( izz) * a ( izz)

          b ( izz) = 0.5 * rhorepl ( izz) * ( vlxrepl ( izz)**2
     1      + vlyrepl ( izz)**2 + vlzrepl ( izz)**2)

          prsrepl ( izz) = ( prsrepl ( izz) - b ( izz))
     1                   * ( gamma_eos_repl ( izz) - 1.0)
#if ( PRESSURE_POSITIVITY == 1)
          entropyrepl ( izz) = entropygr ( ixx, iyy, izz + 1, 1)
!     1            + dx_shift * entropygr ( ixx, iyy, izz + 1, 2)
!     1            + dy_shift * entropygr ( ixx, iyy, izz + 1, 3)
     1                 - 0.5 * entropygr ( ixx, iyy, izz + 1, 4)
     1            + dt_shift * entropygr ( ixx, iyy, izz + 1, 5)
          IF ( prsrepl ( izz) .LE. prsfloor) prsrepl ( izz) =
     1    entropyrepl ( izz)
     1  * rhorepl ( izz)**( gamma_eos_repl ( izz) - 1.0)
#endif
          prsrepl ( izz) = AMAX1 ( prsfloor, prsrepl ( izz))

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, izz) =
     1             u ( ixx, iyy, izz + 1, nspecies_min: nspecies_max, 1)
     1     - 0.5 * u ( ixx, iyy, izz + 1, nspecies_min: nspecies_max, 4)
     1+ dt_shift * u ( ixx, iyy, izz + 1, nspecies_min: nspecies_max, 5)

          specierepl ( 1: NFLUID_EULER, izz) =
     1    specierepl ( 1: NFLUID_EULER, izz) * a ( izz)
#endif

! ----------

          msonic_arr ( izz) = 0.5 * ( msonicgr ( ixx, iyy, izz)
     1                              + msonicgr ( ixx, iyy, izz + 1))
          divvel_arr ( izz) = 0.5 * ( divvelgr ( ixx, iyy, izz)
     1                              + divvelgr ( ixx, iyy, izz + 1))

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

! Zero out all other variables that go through RS interface. This is not
! needed but is done to overcome deficiencies in certain compilers.
! They stop allocating memory for automatic arrays if they are not initialized.

          lf_eigenvec_re ( izz, :, :) = 0.0
          rt_eigenvec_re ( izz, :, :) = 0.0
          eigenval_re ( izz, :) = 0.0
          eigenwt_re ( izz, :) = 0.0

          conserved_re ( izz, :) = 0.0
          flux_re ( izz, :) = 0.0

          zonebdy_flux ( izz, :, :) = 0.0
          zonebdy_conserved ( izz, :, :) = 0.0

          rho0 ( izz) = 0.0
          prs0 ( izz) = 0.0
          vlx0 ( izz) = 0.0
          vly0 ( izz) = 0.0
          vlz0 ( izz) = 0.0
          entropy_flux_re ( izz) = 0.0

          specie0 ( :, izz) = 0.0

          do_ein ( izz) = 0

! ----------

          END DO ! End of DO izz = loop.

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Continuitive boundary conditions in the z-direction shown below.

          IF ( bcarr_z ( ixx, iyy, 1) .EQ. 2) THEN

          rhoremi ( nriembegin + one_z) = rhorepl ( nriembegin + one_z)
          prsremi ( nriembegin + one_z) = prsrepl ( nriembegin + one_z)

          vlxremi ( nriembegin + one_z) = vlxrepl ( nriembegin + one_z)
          vlyremi ( nriembegin + one_z) = vlyrepl ( nriembegin + one_z)
          vlzremi ( nriembegin + one_z) = vlzrepl ( nriembegin + one_z)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_z) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_z)
#endif

          END IF

          IF ( bcarr_z ( ixx, iyy, 2) .EQ. 2) THEN

          rhorepl ( nriemend - one_z) = rhoremi ( nriemend - one_z)
          prsrepl ( nriemend - one_z) = prsremi ( nriemend - one_z)

          vlxrepl ( nriemend - one_z) = vlxremi ( nriemend - one_z)
          vlyrepl ( nriemend - one_z) = vlyremi ( nriemend - one_z)
          vlzrepl ( nriemend - one_z) = vlzremi ( nriemend - one_z)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_z) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_z)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Reflection boundary conditions in the z-direction shown below.

          IF ( bcarr_z ( ixx, iyy, 1) .EQ. 3) THEN

          rhoremi ( nriembegin + one_z) = rhorepl ( nriembegin + one_z)
          prsremi ( nriembegin + one_z) = prsrepl ( nriembegin + one_z)

          vlxremi ( nriembegin + one_z) = vlxrepl ( nriembegin + one_z)
          vlyremi ( nriembegin + one_z) = vlyrepl ( nriembegin + one_z)
          vlzremi ( nriembegin + one_z)= -vlzrepl ( nriembegin + one_z)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_z) = 
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_z)
#endif

          END IF

          IF ( bcarr_z ( ixx, iyy, 2) .EQ. 3) THEN

          rhorepl ( nriemend - one_z) = rhoremi ( nriemend - one_z)
          prsrepl ( nriemend - one_z) = prsremi ( nriemend - one_z)

          vlxrepl ( nriemend - one_z) = vlxremi ( nriemend - one_z)
          vlyrepl ( nriemend - one_z) = vlyremi ( nriemend - one_z)
          vlzrepl ( nriemend - one_z)= -vlzremi ( nriemend - one_z)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_z) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_z)
#endif

          END IF

! ----------------------------------------------------------------------

! Enforce some boundary conditions by hand if needed.
! Periodic boundary conditions in the z-direction shown below. Works in serial.

          IF ( bcarr_z ( ixx, iyy, 1) .EQ. 4) THEN

          rhoremi ( nriembegin + one_z) = rhoremi ( nriemend - one_z)
          prsremi ( nriembegin + one_z) = prsremi ( nriemend - one_z)

          vlxremi ( nriembegin + one_z) = vlxremi ( nriemend - one_z)
          vlyremi ( nriembegin + one_z) = vlyremi ( nriemend - one_z)
          vlzremi ( nriembegin + one_z) = vlzremi ( nriemend - one_z)

#if ( NFLUID_EULER >= 1)
          specieremi ( 1: NFLUID_EULER, nriembegin + one_z) =
     1    specieremi ( 1: NFLUID_EULER, nriemend - one_z)
#endif

          END IF

          IF ( bcarr_z ( ixx, iyy, 2) .EQ. 4) THEN

          rhorepl ( nriemend - one_z) = rhorepl ( nriembegin + one_z)
          prsrepl ( nriemend - one_z) = prsrepl ( nriembegin + one_z)

          vlxrepl ( nriemend - one_z) = vlxrepl ( nriembegin + one_z)
          vlyrepl ( nriemend - one_z) = vlyrepl ( nriembegin + one_z)
          vlzrepl ( nriemend - one_z) = vlzrepl ( nriembegin + one_z)

#if ( NFLUID_EULER >= 1)
          specierepl ( 1: NFLUID_EULER, nriemend - one_z) =
     1    specierepl ( 1: NFLUID_EULER, nriembegin + one_z)
#endif

          END IF

! ----------------------------------------------------------------------

! Call the Riemann solver with cyclic rotation of input variables.
! HLL is the default.

#if ( USE_RIEMANN_SOLVER == 1)
      CALL EULER_RIEM_ROE
#elif ( USE_RIEMANN_SOLVER == 2)
      CALL EULER_RIEM_HLLC
#elif ( USE_RIEMANN_SOLVER == 3)
      CALL EULER_RIEM_HLL
#else
      CALL EULER_RIEM_HLLEM_DB
#endif

     1   ( nriembegin, nriemend,

     1     gamma_euler, smallnum, prsfloor, rhofloor, fracdiff,
     1     ientropyfix, iprimitive, use_approx_eos,
     1     split_prsflux_re, divratio_einfeldt_rs, debug_stop,

     1     gasconst_euler,
     1  cp_specie_euler, molewt_specie_euler, heat_of_formation_euler,

     1        rhoremi, rhorepl,
     1        prsremi, prsrepl,
     1        vlzremi, vlzrepl,
     1        vlxremi, vlxrepl,
     1        vlyremi, vlyrepl,
     1        specieremi, specierepl,

     1        gamma_eos_remi, gamma_eos_repl,
     1        gamma_soundspeed_remi, gamma_soundspeed_repl,
     1        mean_mol_wt_remi, mean_mol_wt_repl,

     1        msonic_arr, divvel_arr,

     1        zonebdy_conserved, zonebdy_flux,

     1     eigenval_re, lf_eigenvec_re, rt_eigenvec_re, eigenwt_re,
     1     conserved_re, flux_re, entropy_flux_re,

     1     rho0, prs0, vlz0, vlx0, vly0, specie0,
     1     do_ein )

! ----------------------------------------------------------------------

! Add in the contributions from the lapidus viscosity to the 1d fluxes.

          IF ( lapidus_coef .GT. 0.0) THEN

          DO izz = nriembegin, nriemend

          flux_re ( izz, 1) = flux_re ( izz, 1) + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, 1, 1) - u ( ixx, iyy, izz + 1, 1, 1) )

          flux_re ( izz, 3) = flux_re ( izz, 3) + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, 2, 1) - u ( ixx, iyy, izz + 1, 2, 1) )
          flux_re ( izz, 4) = flux_re ( izz, 4) + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, 3, 1) - u ( ixx, iyy, izz + 1, 3, 1) )
          flux_re ( izz, 2) = flux_re ( izz, 2) + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, 4, 1) - u ( ixx, iyy, izz + 1, 4, 1) )

          flux_re ( izz, 5) = flux_re ( izz, 5) + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, 5, 1) - u ( ixx, iyy, izz + 1, 5, 1) )

#if ( NFLUID_EULER >= 1)
          flux_re ( izz, nspecies_min + 1: nspecies_max + 1) =
     1    flux_re ( izz, nspecies_min + 1: nspecies_max + 1)
     1  + lapidus_arr ( izz)
     1  * ( u ( ixx, iyy, izz, nspecies_min: nspecies_max, 1)
     1    - u ( ixx, iyy, izz + 1, nspecies_min: nspecies_max, 1) )
#endif

          END DO

          END IF

! ----------------------------------------------------------------------

! Hand back fluxes. Use cyclic rotation for fluxes.
! Fluxes have been evaluated at face centers.

          DO izz = nriembegin, nriemend

          flux_z ( ixx, iyy, izz, 1) = flux_re ( izz, 1)
          flux_z ( ixx, iyy, izz, 2) = flux_re ( izz, 3)
          flux_z ( ixx, iyy, izz, 3) = flux_re ( izz, 4)
          flux_z ( ixx, iyy, izz, 4) = flux_re ( izz, 2)
          flux_z ( ixx, iyy, izz, 5) = flux_re ( izz, 5)

          prs_flux_z ( ixx, iyy, izz) = flux_re ( izz, 6)

#if ( NFLUID_EULER >= 1)
          flux_z ( ixx, iyy, izz, nspecies_min: nspecies_max) =
     1              flux_re ( izz, nspecies_min + 1: nspecies_max + 1)
#endif

#if ( PRESSURE_POSITIVITY == 1)
          entropy_flux_z ( ixx, iyy, izz) = entropy_flux_re ( izz)
#endif

          END DO ! End of DO izz = loop.

! ----------

! If we are doing the predictor stage then hand back the scaled time
! rate of update to the zone-centered "u ( ixx, iyy, izz, :, 5)".
! Note shift to zone-centered index.

          IF ( ABS ( dt_shift) .LT. smallnum) THEN

! ----------
! Cartesian & Cylindrical

          IF ( ( igeom .EQ. 1) .OR. ( igeom .EQ. 2)) THEN

          DO izz = nriembegin + 1, nriemend

          a ( izz) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( izz)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( izz)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( izz)

          END DO ! End of DO izz = loop.

          END IF

! -----
! Spherical

          IF ( igeom .EQ. 3) THEN

          DO izz = nriembegin + 1, nriemend

          a ( izz) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1             * dz_save ( izz)
          a ( izz) = dtcur / a ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 5) = 
     1    u ( ixx, iyy, izz, 1: nspecies_max, 5)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( izz)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 5) = entropygr ( ixx, iyy, izz, 5)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( izz)
#endif

          u ( ixx, iyy, izz, 4, 5) = u ( ixx, iyy, izz, 4, 5)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( izz)

          END DO ! End of DO izz = loop.

          END IF

! ----------

          END IF ! End of "dt_shift .EQ. 0.0"

! ----------------------------------------------------------------------

        END DO ! End of DO ixx = loop.

      END DO ! End of DO iyy = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE MAKE_EULER_FLUX_Z

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ADD_PREDICTOR_SOURCE_TERMS_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone,
     1        dtcur, 

     1        u,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
     1        gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1        temperaturegr, divvelgr, msonicgr, flattengr,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in "u" and adds in source terms to the time
! rate of update. This is stored in "u ( ixx, iyy, izz, :, 5)" within
! each zone.

! It is effectively Stage I, i.e. the predictor stage, of a two-stage
! update.

! External body forces have to be added in by hand.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone

      REAL dtcur

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::  a, b, c,
     1        vlx_accn, vly_accn, vlz_accn

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

! ----------------------------------------------------------------------

! Only external body forces need to be included in cartesian geometry.
! This part is just an instantiation of how acceleration terms contribute
! to the momentum and energy densities in the predictor part of the timestep.

#if ( 0 == 1)
      IF ( igeom .EQ. 1) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c,                     &
!$OMP&    vlx_accn, vly_accn, vlz_accn)                                 &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, nspecies_min, nspecies_max,                            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr)                         &
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
!$OMP&    dtcur, nspecies_min, nspecies_max,                            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr)                         &
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
     1      + ( rhogr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)**2 )
     1      * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1      + vly_accn ( ixx)
     1      - ( rhogr ( ixx, iyy, izz) * vlxgr ( ixx, iyy, izz)
     1                                 * vlygr ( ixx, iyy, izz) )
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
!$OMP&    dtcur, nspecies_min, nspecies_max,                            &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr)                         &
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
     1    + rhogr ( ixx, iyy, izz) * ( vlygr ( ixx, iyy, izz)**2
     1                               + vlzgr ( ixx, iyy, izz)**2 )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 5) = u ( ixx, iyy, izz, 3, 5)
     1    + vly_accn ( ixx) - rhogr ( ixx, iyy, izz)
     1    * ( vlxgr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)
     1      - vlzgr ( ixx, iyy, izz)**2 * c ( ixx) )
     1    * dtcur / xc_save ( ixx)

          a ( ixx) = rhogr ( ixx, iyy, izz) * vlxgr ( ixx, iyy, izz)
     1                               * vlzgr ( ixx, iyy, izz)

          b ( ixx) = rhogr ( ixx, iyy, izz) * vlygr ( ixx, iyy, izz)
     1                               * vlzgr ( ixx, iyy, izz)

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

      END SUBROUTINE ADD_PREDICTOR_SOURCE_TERMS_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE UPDATE_VARS_CORRECTOR_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        expand_by_one_zone,
     1        dtcur,

     1        u, flux_x, flux_y, flux_z,
     1        prs_flux_x, prs_flux_y, prs_flux_z,
     1       entropygr, entropy_flux_x, entropy_flux_y, entropy_flux_z,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

! 2) Puts in the time-centered source terms in various geometries.
! External body forces have to be added in by hand.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "riemann.com"

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
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
     1                  n_cc_components) :: flux_x, flux_y, flux_z

      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::

     1      prs_flux_x, prs_flux_y, prs_flux_z,
     1      entropy_flux_x, entropy_flux_y, entropy_flux_z

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::

     1     a, b, c, d, e, vlx_accn, vly_accn, vlz_accn,

     1     rho_temp, vlx_temp, vly_temp, vlz_temp,

     1     delx1, delx2, dely1, dely2, delz1, delz2, areaint

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

      IF ( ( expand_by_one_zone .LT. 0) .OR.
     1     ( expand_by_one_zone .GT. 1) ) THEN
      WRITE ( 6, *)
     1"Stopping in UPDATE_VARS_CORRECTOR_EULER : expand_by_one_zone = ",
     1       expand_by_one_zone
      STOP
      END IF

! ----------------------------------------------------------------------
! Update the zone centered variables.
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

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_x ( ixx, iyy, izz)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
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

! Transform the y-momentum density to an angular momentum density.
! Make corresponding changes for the fluxes for that zone-centered component.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, u, xc_save)                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1                             * xc_save ( ixx)
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

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

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_y ( ixx, iyy, izz)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dy_save ( iyy)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1  - ( prs_flux_z ( ixx, iyy, izz)
     1    - prs_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)

! -----

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ------------------------------

! Transform the y-angular momentum density back to a momentum density.
! As long as the fluxes are not reused (no AMR for non-cartesian) this is ok.
! As long as the coordinate system does not touch the origin, it is easy
! to go from angular momentum fluxes to momentum fluxes.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, u, xc_save)                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
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

! Transform the z-momentum density to an angular momentum density.
! Make corresponding changes for the fluxes for that zone-centered component.

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx)                              &
!$OMP&            SHARED ( ix1, iy1, iz1, u, xc_save, yc_save)          &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1, iz1
        DO iyy = 1, iy1
          DO ixx = 1, ix1
          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1                       * xc_save ( ixx) * SIN ( yc_save ( iyy))
          END DO
        END DO
      END DO

!$OMP END PARALLEL DO

! -----

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

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_x ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_x ( ixx - 1, iyy, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_x ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_x ( ixx - 1, iyy, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dx_save ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1)
     1  - ( prs_flux_x ( ixx, iyy, izz)
     1    - prs_flux_x ( ixx - 1, iyy, izz) ) * a ( ixx)

! -----

          a ( ixx) = - xc_save ( ixx) * ( COS ( yb_save ( iyy))
     1                                  - COS ( yb_save ( iyy - 1)) )
          a ( ixx) = dtcur / a ( ixx)

          b ( ixx) = SIN ( yb_save ( iyy))

          c ( ixx) = SIN ( yb_save ( iyy - 1))

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_y ( ixx, iyy, izz, 1: nspecies_max) * b ( ixx)
     1    - flux_y ( ixx, iyy - 1, izz, 1: nspecies_max) * c ( ixx) )
     1  * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_y ( ixx, iyy, izz) * b ( ixx)
     1    - entropy_flux_y ( ixx, iyy - 1, izz) * c ( ixx) ) * a ( ixx)
#endif

          a ( ixx) = xc_save ( ixx) * dy_save ( iyy)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1  - ( prs_flux_y ( ixx, iyy, izz)
     1    - prs_flux_y ( ixx, iyy - 1, izz) ) * a ( ixx)

! -----

          a ( ixx) = xc_save ( ixx) * SIN ( yc_save ( iyy))
     1             * dz_save ( izz)
          a ( ixx) = dtcur / a ( ixx)

          u ( ixx, iyy, izz, 1: nspecies_max, 1) =
     1    u ( ixx, iyy, izz, 1: nspecies_max, 1)
     1  - ( flux_z ( ixx, iyy, izz, 1: nspecies_max) 
     1    - flux_z ( ixx, iyy, izz - 1, 1: nspecies_max) ) * a ( ixx)

#if ( PRESSURE_POSITIVITY == 1)
          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)
     1  - ( entropy_flux_z ( ixx, iyy, izz)
     1    - entropy_flux_z ( ixx, iyy, izz - 1) ) * a ( ixx)
#endif

          a ( ixx) = dtcur / dz_save ( izz)

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
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
          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
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
!$OMP&    vlx_temp, vly_temp, vlz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr)                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                    + rhogr ( ixx, iyy, izz) )
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

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1)
     1                      + vlx_accn ( ixx) * vlx_temp ( ixx)
     1                      + vly_accn ( ixx) * vly_temp ( ixx)
     1                      + vlz_accn ( ixx) * vlz_temp ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 1".
#endif

! ----------------------------------------------------------------------
! Cylindrical

      IF ( igeom .EQ. 2) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c, rho_temp,           &
!$OMP&    vlx_temp, vly_temp, vlz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr)                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                    + rhogr ( ixx, iyy, izz) )
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

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1)
     1      + vlx_accn ( ixx) + rho_temp ( ixx) * vly_temp ( ixx)**2
     1                        * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1                             + vly_accn ( ixx)

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1                             + vlz_accn ( ixx)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1)
     1                             + vlx_accn ( ixx) * vlx_temp ( ixx)
     1                             + vly_accn ( ixx) * vly_temp ( ixx)
     1                             + vlz_accn ( ixx) * vlz_temp ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 2".

! ----------------------------------------------------------------------
! Spherical

      IF ( igeom .EQ. 3) THEN

!$OMP PARALLEL DO PRIVATE ( izz, iyy, ixx, a, b, c, rho_temp,           &
!$OMP&    vlx_temp, vly_temp, vlz_temp, vlx_accn, vly_accn, vlz_accn)   &
!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&    dtcur, nspecies_min, nspecies_max, xb_save, dx_save, xc_save, &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    u, rhogr, vlxgr, vlygr, vlzgr, prsgr)                         &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO izz = 1 - MIN0 ( 1, ioffz), iz1 + MIN0 ( 1, ioffz)

        DO iyy = 1 - MIN0 ( 1, ioffy), iy1 + MIN0 ( 1, ioffy)

          DO ixx = 1 - MIN0 ( 1, ioffx), ix1 + MIN0 ( 1, ioffx)

          a ( ixx) = 1.0 / u ( ixx, iyy, izz, 1, 1)

          rho_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 1, 1)
     1                    + rhogr ( ixx, iyy, izz) )
          vlx_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 2, 1) * a ( ixx)
     1                    + vlxgr ( ixx, iyy, izz) )
          vly_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 3, 1) * a ( ixx)
     1                    + vlygr ( ixx, iyy, izz) )
          vlz_temp ( ixx) = 0.5 * ( u ( ixx, iyy, izz, 4, 1) * a ( ixx)
     1                    + vlzgr ( ixx, iyy, izz) )

          c ( ixx) = COS ( yc_save ( iyy)) / SIN ( yc_save ( iyy))

! Load in with x-acceleration. Similarly for others.
          vlx_accn ( ixx) = 0.0
          vly_accn ( ixx) = 0.0
          vlz_accn ( ixx) = 0.0

! Uncomment if body forces are actually present.
!          vlx_accn ( ixx) = dtcur * rho_temp ( ixx) * vlx_accn ( ixx)
!          vly_accn ( ixx) = dtcur * rho_temp ( ixx) * vly_accn ( ixx)
!          vlz_accn ( ixx) = dtcur * rho_temp ( ixx) * vlz_accn ( ixx)

          u ( ixx, iyy, izz, 2, 1) = u ( ixx, iyy, izz, 2, 1)
     1  + vlx_accn ( ixx)
     1  + rho_temp ( ixx) * ( vly_temp ( ixx)**2 + vlz_temp ( ixx)**2 )
     1  * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 3, 1) = u ( ixx, iyy, izz, 3, 1)
     1    + vly_accn ( ixx)
     1    - rho_temp ( ixx) * ( vlx_temp ( ixx) * vly_temp ( ixx)
     1                        - vlz_temp ( ixx)**2 * c ( ixx) )
     1    * dtcur / xc_save ( ixx)

          u ( ixx, iyy, izz, 4, 1) = u ( ixx, iyy, izz, 4, 1)
     1                             + vlz_accn ( ixx)

          u ( ixx, iyy, izz, 5, 1) = u ( ixx, iyy, izz, 5, 1)
     1                             + vlx_accn ( ixx) * vlx_temp ( ixx)
     1                             + vly_accn ( ixx) * vly_temp ( ixx)
     1                             + vlz_accn ( ixx) * vlz_temp ( ixx)

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

      END IF ! End of "igeom .EQ. 3".

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE UPDATE_VARS_CORRECTOR_EULER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE EVALUATE_TIMESTEP_EULER

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_components, n_cc_modes,
     1        dtcur, dtcurp1,

     1        u, entropygr,

     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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
     1        n_cc_components, n_cc_modes

      REAL dtcur, dtcurp1

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) :: entropygr

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1        rhogr, prsgr, vlxgr, vlygr, vlzgr,
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

      REAL, DIMENSION ( ONEDFLUIDBEGIN : ONEDFLUIDEND) ::
     1        a, b, c, d, e

! ----------------------------------------------------------------------

      nspecies_min = NFIELD_EULER + 1
      nspecies_max = NFIELD_EULER + NFLUID_EULER

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
!$OMP&    igeom, cfl_coef, rhofloor, prsfloor,                          &
!$OMP&    xb_save, dx_save, xc_save,                                    &
!$OMP&    yb_save, dy_save, yc_save, zb_save, dz_save, zc_save,         &
!$OMP&    x_indx_limits, u, entropygr, rhogr, prsgr, vlxgr, vlygr,      &
!$OMP&    vlzgr, gamma_eos_gr, gamma_soundspeed_gr, msonicgr)           &
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

          e ( ixx) = ( vlxgr ( ixx, iyy, izz)**2
     1      + vlygr ( ixx, iyy, izz)**2 + vlzgr ( ixx, iyy, izz)**2)

          b ( ixx) = 0.5 * rhogr ( ixx, iyy, izz) * e ( ixx)

          prsgr ( ixx, iyy, izz) = ( u ( ixx, iyy, izz, 5, 1)
     1                             - b ( ixx) )
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
          u ( ixx, iyy, izz, 5, 1) = b ( ixx) + prsfloor
     1                         / ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)
          END IF

! ----------

          msonicgr ( ixx, iyy, izz) =
     1             gamma_soundspeed_gr ( ixx, iyy, izz)
     1           * prsgr ( ixx, iyy, izz) * a ( ixx)
          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

          e ( ixx) = SQRT ( e ( ixx))

! ----------

          IF ( igeom .EQ. 1) THEN

            d ( ixx) = AMIN1 ( dx_save ( ixx), dy_save ( iyy),
     1                         dz_save ( izz))

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

          d ( ixx) = cfl_coef * d ( ixx)
     1             / ( msonicgr ( ixx, iyy, izz) + e ( ixx))

          dtcurp1 = AMIN1 ( dtcurp1, d ( ixx))

! ----------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE EVALUATE_TIMESTEP_EULER

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








