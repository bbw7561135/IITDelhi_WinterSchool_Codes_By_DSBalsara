! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#include "directives"

! Set the maximum number of zones that can come together at a vertex.
!#define MAX_ZONES_PER_VRTX 4

! This version only uses HLLEM and only works with HLL. I.e. it is
! the equivalent of "USE_HLLEM 1" and "NWAVES_1DRS = 2".

! The only trustable outputs from this code are "U_Star, F_Star, G_Star".
! I.e. for 1D RS, call the 1D HLLEM RS separately.

#define NWAVES_1DRS 2

! "BUILD_SUBSTRUCTURE = 1" builds sub-structure in the strongly
! interacting state. It can improve the solution in the subsonic case.
! However, we frequently use only a subset of the eigenspace. Because the projection
! into that eigenspace is imperfect, we could generate relic velocity in the third
! dimension. If that is undesirable, turn off the substructure. It only increases
! dissipation by a modest amount.
#define BUILD_SUBSTRUCTURE 1

! For v. stringent problems, one can ask for the wave speeds in the 1D subsonic
! case to be extremized. In that event, set below parameter to
! 0 (no expansion) exact shock-capturing speeds from HLL/C RS used;
! 1 (expand waves by a little);
! or 2 (more extreme fix; rebuild wave model; 1D RS all-subsonic).
#define EXTREMIZE_WAVESPEEDS 0

! There are lots of unnecessary warnings. Decide whether to report them.
#define REPORT_WARNING 0

! ADER schemes might benefit from providing input states and also fluxes.
! "SAVE_FLUX_MODES 0" when only state is input; "1" when fluxes also input.
! Also watch out for setting of "PRESSURE_POSITIVITY".
#define SAVE_FLUX_MODES 0

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MuSIC_MHD_2DRS_1_STRUCT

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, prsfloor, rhofloor, pi,

     1     blend_HLL, blend_MDHLL,
     1     Prs_Tot_Star, Prs_Gas_Star,

     1     Vlx_ALE, Vly_ALE,

     1     U, F, G, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     U_Star, F_Star, G_Star,
     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Implements the genuinely two dimensional MuSIC Riemann solver in
! 2d for MHD. The self-similarity formulation from Balsara (2014) JCP
! is used. It is referred to as B14. The multidimensional MuSIC RS
! based on 1D HLLEM is from Balsara et al. (2016) JCP. It is referred
! to as Betal16. Dumbser and Balsara (2016) JCP is referred to as DB16.

! The equations in this subroutine are keyed to B14 and Betal16 and DB16.
! The equations in the other subroutines might be keyed to BDA14
! (Balsara, Dumbser & Abgrall 2014 JCP).

! This version is specific to structured meshes.

! On input:

! "nzone_vrtx" number of zones/states that come together at a vertex.

! "multid_wavemodel_num" Make a choice for the multidimensional wave model.
! See notes in "MULTID_WAVE_MODEL" for more info. In this version, it
! is suppressed when "MULTID_WAVE_MODEL_MINIMALIST" is invoked to build
! the multidimensional wave model.

! We can handle up to "MAX_ZONES_PER_VRTX" states coming together at the
! vertex in question. Of these, we only use the first "nzone_vrtx" of them.

! "NFIELD" is always set to seven. "density, x,y,z-momentum density,
! total energy density, bfx, bfy, bfz" is the order of the variables.

! "NFLUID" is zero when there are no multifluids. If multiple fluid
! species are present, it adds that many components of species mass density
! to the vector of conserved variables.

! "smallnum" is a small number used for wavespeed math. Set to 1.0e-12.
! "prsfloor" is the floor value for the pressure.
! "rhofloor" is the floor value of the density.

! "U, gamma_eos, gamma_sdsp" conserved states in an anti-clockwise, cyclic
! arrangement around the vertex of interest. These are provided as input. 

! The use of fluxes from "F, G", instead of building them from
! scratch, is only in the case of ADER schemes where such information can
! be provided. Even then, it is indicated by setting "SAVE_FLUX_MODES" to 1.

! "nx, ny" unit normals along the zone boundaries that separate the states.
! The first unit normal separates the state 1 and 2. The arrangement of
! unit normals is similarly anti-clockwise and cyclic around the vertex.

! "entropy" entropy density arranged in anti-clockwise, cyclic fashion.
! If the pressure hits floor value, it gives us a second chance to recover
! from that kind of bad situation.

! "blend_HLL, blend_MDHLL" amount of 1d and multid. HLL flux to be
! blended in with the multid. HLLC flux. A small fraction may be 
! stabilizing for some problems.

! "Vlx_ALE, Vly_ALE" Velocity of the ALE mesh, specified as input.
! The state and fluxes are evaluated along the characteristic lines that
! move with this velocity.


! On output:

! "U_Star, F_Star, G_Star" multid. HLLC state and fluxes in the global frame.
! This is the thing we work hard to find.

! "Prs_Tot_Star" Total gas + magnetic pressure in the multid. HLLC state.
! "Prs_Gas_Star" gas pressure in the multid. HLLC state.
! This is a DEPRECATED feature and could be removed from the subroutine
! interface.

! "U_HLLC, F_HLLC, G_HLLC" the state and fluxes from the 1d Riemann
! problems that are solved as a part of the algorithm in obtaining the
! multid. HLLC state and fluxes. "F_HLLC" is especially useful because
! it is the 1d flux in the direction that is orthogonal to the zone
! boundaries, i.e. in the direction "(-ny, nx)". In other words, it is
! the flux one would traditionally use in a 1d Riemann solver approach.

! "F_HLLC_Global_x, G_HLLC_Global_y" The previous 1d fluxes may also be
! needed in a global coordinate system and are provided here.

! "F_entropy_flux_x, G_entropy_flux_y" analogous to above for the flux
! of entropy density. Evolving entropy density gives us a second chance
! to recover from a situation where pressure has hit floor value.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx, multid_wavemodel_num

      REAL smallnum, prsfloor, rhofloor, pi,

! Amount of HLL we will have to blend in for stabilization.
     1     blend_HLL, blend_MDHLL,

! Velocity of vertex; total and gas pressure in resolved state.
     1     Vlx_ALE, Vly_ALE,
     1     Prs_Tot_Star, Prs_Gas_Star

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables and fluxes that come together at a vertex.
     1      U, F, G,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

! entropy variables that go through RS interface.
     1      entropy, F_entropy_flux_x, G_entropy_flux_y,

! Components of the unit normal.
     1      nx, ny

! Final, strongly-interacting state and Global x- and y-fluxes.
! These are the states that straddle the time axis.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: U_Star, F_Star, G_Star

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      INTEGER iHLLEM_choice, i, j, k, ip1, im1, itempaa,

! Is the wave model supersonic? Do we have a good wave model? Can we linearize
! around the state that we found?
     1        Is_SuperSonic, Is_SuperSonic_State, Ifail_WavModl,
     1        Ifail, Linearize_Around_U_Star,

     1        iprimitive, idirection, debug_stop

      REAL fourpi, ob_4pi, ob_8pi, expand_wavespeeds,
     1     tempaa, tempab, tempac, tempad, tempae, tempaf, tempag,
     1     temp_0, temp_x, temp_y, l0, l1, l2,
     1     small_frac, small_speed, small_rho,
     1     Mass_Matrix_Lin_Var_x, Mass_Matrix_Lin_Var_y,

! Variables used for multid. limiting.
     1     psi_limiter, min_limiter, largest_min_limiter,

! Supersonic speeds in the global frame; always built even if flow isn't
! supersonic. Also provide the maximal speed in this RP. (output)
     1     Vlx_Centroid, Vly_Centroid, S_Radial_Max,

! "yF_Lo, yF_Hi" are points of intersection where x = 0 intersects wave model.
     1     yF_Lo, yF_Hi,

! "xG_Lo, xG_Hi" are points of intersection where y = 0 intersects wave model.
     1     xG_Lo, xG_Hi,

! First 10 moments of x and y over areas contributing to "U", "F" and "G".
     1     Int_Momnt_U ( 10), Int_Momnt_F ( 10), Int_Momnt_G ( 10),

! Density, velocities and magnetic fields in the resolved subsonic Multid
! state. This is the state around which we linearize the wave model.
     1     Rho_Star, Vlx_Star, Vly_Star, Vlz_Star, Eng_Tot_Star,
     1     Bfx_Star, Bfy_Star, Bfz_Star, V2_Star, B2_Star,
     1     V_Dot_B_Star, Gamma_Star, ob_Rho, gasconst, Xbar,

! Density and pressure min and max from all the incoming states.
! Also, density and pressure min and max at all the vertices of the wave model.
     1     Rho_Neighbor_Min, Rho_Neighbor_Max, Prs_Neighbor_Min,
     1     Prs_Neighbor_Max, Rho_Vertex_Min, Rho_Vertex_Max,
     1     Prs_Vertex_Min, Prs_Vertex_Max,

! Speeds relative to the centroid of the wave model.
     1     Sx_h, Sy_h,

! Speeds that bound the wave model.
     1     S_L, S_R, S_D, S_U

      REAL, DIMENSION ( 0: NFLUID) :: Specie_Star


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! We also retain the primitive variables.
     1      V,

! The four resolved states and their fluxes in the x- and y-directions.
     1     U_Minus, U_Star_Minus, U_DblStar_Minus,
     1                            U_DblStar_Plus, U_Star_Plus, U_Plus,
     1     F_Minus, F_Star_Minus, F_DblStar_Minus,
     1                            F_DblStar_Plus, F_Star_Plus, F_Plus,
     1     G_Minus, G_Star_Minus, G_DblStar_Minus,
     1                            G_DblStar_Plus, G_Star_Plus, G_Plus,

! Averages of the two resolved states and two resolved y-fluxes. From 1D HLL RS.
     1      U_Star_Avg, F_Star_Avg, G_Star_Avg


      INTEGER, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) :: IF_Mid, IG_Mid


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! The left, left-Alfven, middle, right-Alfven, and right speeds of the
! 1D HLLD RS in the Tau-direction,
! "S_Tau_Min, S_Tau_Alf_Minus, S_Tau_Mid, S_Tau_Alf_Plus, S_Tau_Max".
! W/o Alfven, it is also used for 1D HLLC RS.
     1      S_Tau_Min, S_Tau_Alf_Minus, S_Tau_Mid,
     1                 S_Tau_Alf_Plus, S_Tau_Max,

! The extremal speeds in the Neta-direction, "S_Neta_Min, S_Neta_Max".
     1      S_Neta_Min, S_Neta_Max,

! Fastest speed along the principal direction. Extents of the local
! wavefront in the Tau-direction. Intersection points of the wave fronts.
     1      S_Neta_Multid, L_Tau_Min, L_Tau_Max, x, y,

! Unit vector along the edge of a wave model. Dot products of that unit vector
! with "n" and "t" unit vectors.
     1     neta_x, neta_y, neta_dot_n, neta_dot_t,

! "LF_Tau_Mid, LG_Tau_Mid" are mid-points of the facial integration.
     1     LF_Tau_Mid, LG_Tau_Mid


      REAL, DIMENSION ( MAX_ZONES_PER_VRTX, NWAVES_1DRS) ::

! Wave speeds in the 1D Riemann problem; unparametrized and parametrized.
     1     S_1DRS, L_1DRS


! Conserved states and appropriate fluxes for the 1d Riemann fans.
! Constant states assumed in the Riemann fans -- i.e. we are using
! approximate RS.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID,
     1                  MAX_ZONES_PER_VRTX, 0: NWAVES_1DRS) ::

     1      U_1DRS, F_1DRS, G_1DRS


      REAL, DIMENSION ( NFIELD + 1 + NFLUID) ::

! Panel integrals in the faces of the wave model.
     1      PanelInt_U0, PanelInt_Ux, PanelInt_Uy,

! Variables that bear the sub-structure.
     1      U_Lin_Var, U_Lin_Var_x, U_Lin_Var_y, U_temp,
     1      U_Lim_Lin_Var_x, U_Lim_Lin_Var_y,

! These are used for the multidimensional limiting.
     1      U_Neighbor_Min, U_Neighbor_Max, U_Vertex_Min, U_Vertex_Max

! ----------------------------------------------------------------------

! Panel integrals for HLL fluxes in the faces of the wave model.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) ::

     1      HLL_PanelInt_U0, HLL_PanelInt_Ux, HLL_PanelInt_Uy,

     1      HLL_U_Bar, HLL_F_Bar, HLL_G_Bar

! ----------------------------------------------------------------------

! Matrices used for the least squares process.

      REAL, DIMENSION ( (NWAVES_1DRS + 1) * MAX_ZONES_PER_VRTX, 2) ::
     1                A_LS

      REAL, DIMENSION ( NFIELD + 1 + NFLUID,
     1                  (NWAVES_1DRS + 1) * MAX_ZONES_PER_VRTX) ::
     1                Rhs_LS

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, 2) :: R_22

      REAL A22 ( 2, 2), Inv_A22 ( 2, 2)

! ----------------------------------------------------------------------

! Eigenvalues and eigenweights in x and y-directions.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: eigenval_x, eigenval_y,
     1      eigenwt_x, eigenwt_y, delta_x, delta_y

      REAL, DIMENSION ( NFIELD + 1 + NFLUID, NFIELD + 1 + NFLUID) ::

! These are the characteristic matrices or portions of the characteristic
! matrices that pick out certain inner waves in the wave foliation.
     1     Abar_InnerWaves, Bbar_InnerWaves,

! Eigenvectors in the x and y-directions.
     1      lf_eigenvec_x, rt_eigenvec_x, lf_eigenvec_y, rt_eigenvec_y

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 / fourpi

      small_frac = 1.0e-6

! "psi_limiter" Should lie between [0.5, 1.0].

      psi_limiter = 1.0


! If too much limiting is needed, it may be safest not to have sub-structure.
! Setting "largest_min_limiter = 0.0" means any and all sub-structure is allowed
! as long as it is consistent with the multid. limiter. Setting 
! "largest_min_limiter = 1.0" (or larger than 1.0) means no allowed sub-structure.
! We can make it self-adjusting by keying it to "blend_HLL".
! Since use of "largest_min_limiter" is dicey, we set it to zero and 
! actually don't use it.

      largest_min_limiter = AMAX1 ( 0.2, blend_HLL)
!      largest_min_limiter = 0.0


! Factor by which subsonic wavespeeds are expanded just to be on the safe side.
! "expand_wavespeeds = 1.01" means 1% expansion of wave model.
! "expand_wavespeeds = 1.05" means 5% expansion of wave model.
      expand_wavespeeds = 1.01


! "iHLLEM_choice = 1" gives original Dumbser-Balsara HLLEM; it has
! slighly less dissipation for fast-moving waves. "iHLLEM_choice = 2" gives
! Balsara et al HLLEM; it is more stable for strong shocks.
! Both will capture stationary intermediate waves.
! "iHLLEM_choice = 2" is safer because it builds slightly smaller gradients.

! "iHLLEM_choice = 1 or 2 for HLLEM; 0 for HLL". For this code 2 is best
! because the 2D RS exactly retrieves the 1D HLLEM RS.
      iHLLEM_choice = 2

! ----------------------------------------------------------------------

! Check to make sure that number of zones around the vertex is within allowed
! limits.

      IF ( MAX_ZONES_PER_VRTX < nzone_vrtx) THEN
      WRITE ( 6, *)"STOPPING in MuSIC_MHD_2DRS_1_STRUCT because:"
      WRITE ( 6, *)"NZONES_","PER_VRTX < nzone_vrtx",
     1              MAX_ZONES_PER_VRTX, nzone_vrtx
      STOP
      END IF

! Check to make sure that we are solving the MHD system.

      IF ( NFIELD /= 7) THEN
      WRITE ( 6, *)"STOPPING in MuSIC_MHD_2DRS_1_STRUCT because:"
      WRITE ( 6, *)"NFIELD /= 7", NFIELD
      STOP
      END IF

! Check to make sure that suitable multidimensional wave model is being used.

      IF ( ( multid_wavemodel_num < 0) .OR.
     1     ( multid_wavemodel_num > 3) ) THEN
      WRITE ( 6, *)"STOPPING in MuSIC_MHD_2DRS_1_STRUCT because:"
      WRITE ( 6, *)"multid_wavemodel_num = ", multid_wavemodel_num
!      WRITE ( 6, *)"Present version only deals with subsonic cases"
      STOP
      END IF

! Check to make sure that our 1D RS is appropriate.

      IF ( NWAVES_1DRS /= 2) THEN
      WRITE ( 6, *)"STOPPING in MuSIC_MHD_2DRS_1_STRUCT because:"
      WRITE ( 6, *)"NWAVES_","1DRS = ", NWAVES_1DRS
      STOP
      END IF

! ----------------------------------------------------------------------
! STEP I) Call "ONE_D_MHD_PRIM_F_G_WAVSPD" to get the wave speeds and build
! primitives and fluxes. This also builds "S_Tau_Min, S_Tau_Max, S_Neta_Min,
! S_Neta_Max" that are needed in the multidimensional, subsonic wave model.
! ----------------------------------------------------------------------

      CALL ONE_D_MHD_PRIM_F_G_WAVSPD ( nzone_vrtx,

     1     smallnum, prsfloor, rhofloor, pi, expand_wavespeeds,

     1     U, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     V, F, G,
     1     S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max)

! ----------------------------------------------------------------------
! STEP II) Call "ONE_D_MHD_RS_HLL_TAU" to solve the 1d HLL Riemann
! problem in the "Tau" direction. The multidimensional, wave model can
! be expanded here.

! ----------------------------------------------------------------------

! 1D HLL RS is called. While this subroutine has all the same
! interfaces as "ONE_D_MHD_RS_HLLC_TAU", it is actually a cheater.
! It just builds the HLL RS and puts appropriate placeholders for HLLC
! variables.
! I.e., it is based on the premise that an HLLC without a jump at the
! contact is an HLL.

      CALL ONE_D_MHD_RS_HLL_TAU ( nzone_vrtx,

     1     smallnum, prsfloor, rhofloor, pi, blend_HLL,
     1     expand_wavespeeds,

     1     U, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     V, F, G,

     1     S_Tau_Min, S_Tau_Mid, S_Tau_Max, S_Neta_Min, S_Neta_Max,

     1     U_Minus, U_Star_Minus, U_Star_Plus, U_Plus,
     1     F_Minus, F_Star_Minus, F_Star_Plus, F_Plus,
     1     G_Minus, G_Star_Minus, G_Star_Plus, G_Plus,

     1     U_Star_Avg, F_Star_Avg, G_Star_Avg,

     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)

! ----------------------------------------------------------------------
! STEP III) Store wave speeds (and states) in appropriate fashion for
! later use.
! ----------------------------------------------------------------------

! Wave foliation for 1D HLL RS.

      DO i = 1, nzone_vrtx

      S_1DRS ( i, 1) = S_Tau_Min ( i)
      S_1DRS ( i, 2) = S_Tau_Max ( i)

      U_1DRS ( :, i, 0) = U_Minus ( i, :)
      U_1DRS ( :, i, 1) = U_Star_Avg ( i, :)
      U_1DRS ( :, i, 2) = U_Plus ( i, :)

      F_1DRS ( :, i, 0) = F_Minus ( i, :)
      F_1DRS ( :, i, 1) = F_Star_Avg ( i, :)
      F_1DRS ( :, i, 2) = F_Plus ( i, :)

      G_1DRS ( :, i, 0) = G_Minus ( i, :)
      G_1DRS ( :, i, 1) = G_Star_Avg ( i, :)
      G_1DRS ( :, i, 2) = G_Plus ( i, :)

      END DO

! ----------

#if ( 0 == 1)

      WRITE ( 6, *)" "
      WRITE ( 6, *)"----------------------------------------"
      DO i = 1, nzone_vrtx
      WRITE ( 6, *)"---------- Face i = ", i," ----------"
      WRITE ( 6, *)"S_1DRS = ", S_1DRS ( i, :)
      WRITE ( 6, *)"----------"
        DO j = 0, NWAVES_1DRS
          WRITE ( 6, *)"Face i = ", i," wave j = ", j
          WRITE ( 6, *)"U_1DRS = ", i," ",U_1DRS ( :, i, j)
          WRITE ( 6, *)" "
          WRITE ( 6, *)"F_1DRS = ", i," ",F_1DRS ( :, i, j)
          WRITE ( 6, *)" "
          WRITE ( 6, *)"G_1DRS = ", i," ",G_1DRS ( :, i, j)
          WRITE ( 6, *)" "
          WRITE ( 6, *)"----------"
        END DO
      END DO
      WRITE ( 6, *)" "
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------
! STEP IV) Now form the multidimensional wave model.
! ----------------------------------------------------------------------

      Ifail = 0

      CALL MULTID_WAVE_MODEL

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, blend_HLL, expand_wavespeeds, Vlx_ALE, Vly_ALE,

     1     U, gamma_eos, gamma_sdsp, S_1DRS,

     1     S_Neta_Min, S_Neta_Max, S_Tau_Min, S_Tau_Max, nx, ny,

     1     Vlx_Centroid, Vly_Centroid, S_Radial_Max,
     1     Is_SuperSonic, Is_SuperSonic_State,

     1     S_Neta_Multid, L_Tau_Min, L_Tau_Max, x, y,
     1     neta_x, neta_y, neta_dot_n, neta_dot_t, L_1DRS,

     1     IF_Mid, LF_Tau_Mid, yF_Lo, yF_Hi,
     1     IG_Mid, LG_Tau_Mid, xG_Lo, xG_Hi,
     1     Int_Momnt_U, Int_Momnt_F, Int_Momnt_G, Ifail)

      IF ( Ifail == 1) THEN
      WRITE ( 6, *)"STOPPING in MuSIC_MHD_2DRS_1_STRUCT because:"
     1      ," after MULTID_WAVE_MODEL, Ifail = ", Ifail
      STOP
      END IF


#if ( 0 == 1)

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)"S_Neta_Multid = ", S_Neta_Multid ( 1: nzone_vrtx)

#endif

! ----------------------------------------------------------------------
! STEP V) The states and fluxes in all the supersonic situations can now
! be returned if a supersonic state has been found.

! "Prs_Tot_Star, Prs_Gas_Star" will be filled with good values in this case.
! We just return "U_Star, F_Star, G_Star" from appropriate parts of the
! 1D Riemann fans.

! In principle, for HLLEM, we could have made a finer interpolation in the
! linear profile. But returning a lumped flux, as is done here, is an acceptable
! compromise.

! ----------------------------------------------------------------------

      IF ( Is_SuperSonic /= 0) THEN

      i = Is_SuperSonic
      j = Is_SuperSonic_State

      U_Star ( :) = U_1DRS ( :, i, j)

      F_Star ( :) = G_1DRS ( :, i, j) * nx ( i)
     1            - F_1DRS ( :, i, j) * ny ( i)

      G_Star ( :) = G_1DRS ( :, i, j) * ny ( i)
     1            + F_1DRS ( :, i, j) * nx ( i)


      Gamma_Star = 0.0

      DO i = 1, nzone_vrtx
      Gamma_Star = Gamma_Star + gamma_eos ( i)
      END DO

      Gamma_Star = Gamma_Star / nzone_vrtx


      Rho_Star = AMAX1 ( U_Star ( 1), rhofloor)

      ob_Rho = 1.0 / Rho_Star

      Vlx_Star = ob_Rho * U_Star ( 2)
      Vly_Star = ob_Rho * U_Star ( 3)
      Vlz_Star = ob_Rho * U_Star ( 4)

      Eng_Tot_Star = U_Star ( 5)

      Bfx_Star = U_Star ( 6)
      Bfy_Star = U_Star ( 7)
      Bfz_Star = U_Star ( 8)


      V2_Star = Vlx_Star**2 + Vly_Star**2 + Vlz_Star**2

      B2_Star = Bfx_Star**2 + Bfy_Star**2 + Bfz_Star**2

      V_Dot_B_Star = Vlx_Star * Bfx_Star + Vly_Star * Bfy_Star
     1             + Vlz_Star * Bfz_Star


      Prs_Gas_Star = Eng_Tot_Star - 0.5 * Rho_Star * V2_Star
     1             - ob_8pi * B2_Star
      Prs_Gas_Star = Prs_Gas_Star * ( Gamma_Star - 1.0)
      Prs_Gas_Star = AMAX1 ( Prs_Gas_Star, prsfloor)

      Prs_Tot_Star = Prs_Gas_Star + ob_8pi * B2_Star


      RETURN

      END IF

! ----------------------------------------------------------------------
! STEP VI) Build just the panel integrals needed for HLL. They are
! "HLL_PanelInt_U0, HLL_PanelInt_Ux, HLL_PanelInt_Uy".

! Use the lowest of these modes to build "HLL_U_Bar ( :)", "HLL_F_Bar ( :)"
! and "HLL_G_Bar ( :)".
! ----------------------------------------------------------------------

      HLL_PanelInt_U0 ( :) = 0.0
      HLL_PanelInt_Ux ( :) = 0.0
      HLL_PanelInt_Uy ( :) = 0.0

! ------------------------------------------------------------
! Looping over the side panels of the wave model. Because this is HLL,
! there are only three waves in each side panel.

      DO i = 1, nzone_vrtx

! ------------------------------
! First wave in this face. Integrate "U_Minus" etc. over
! [ L_Tau_Min ( i), L_1DRS ( i, 1)].

      U_temp ( :) = F_Minus ( i, :) * neta_dot_t ( i)
     1            + G_Minus ( i, :) * neta_dot_n ( i)
     1            - U_Minus ( i, :) * S_Neta_Multid ( i)

! ----------

      l1 = L_Tau_Min ( i)
      l2 = L_1DRS ( i, 1)

! ----------

      temp_0 = l2 - l1

      HLL_PanelInt_U0 ( :) = HLL_PanelInt_U0 ( :) + temp_0 * U_temp ( :)

! ----------

      temp_x = x ( i) * ( l2 - l1)
     1       - neta_y ( i) * 0.5 * ( l2**2 - l1**2)

      HLL_PanelInt_Ux ( :) = HLL_PanelInt_Ux ( :) + temp_x * U_temp ( :)

! ----------

      temp_y = y ( i) * ( l2 - l1)
     1       + neta_x ( i) * 0.5 * ( l2**2 - l1**2)

      HLL_PanelInt_Uy ( :) = HLL_PanelInt_Uy ( :) + temp_y * U_temp ( :)

! ------------------------------
! Last wave in this face. Integrate "U_Plus" etc. over
! [ L_1DRS ( i, NWAVES_1DRS), L_Tau_Max ( i)].

      U_temp ( :) = F_Plus ( i, :) * neta_dot_t ( i)
     1            + G_Plus ( i, :) * neta_dot_n ( i)
     1            - U_Plus ( i, :) * S_Neta_Multid ( i)

! ----------

      l1 = L_1DRS ( i, NWAVES_1DRS)
      l2 = L_Tau_Max ( i)

! ----------

      temp_0 = l2 - l1

      HLL_PanelInt_U0 ( :) = HLL_PanelInt_U0 ( :) + temp_0 * U_temp ( :)

! ----------

      temp_x = x ( i) * ( l2 - l1)
     1       - neta_y ( i) * 0.5 * ( l2**2 - l1**2)

      HLL_PanelInt_Ux ( :) = HLL_PanelInt_Ux ( :) + temp_x * U_temp ( :)

! ----------

      temp_y = y ( i) * ( l2 - l1)
     1       + neta_x ( i) * 0.5 * ( l2**2 - l1**2)

      HLL_PanelInt_Uy ( :) = HLL_PanelInt_Uy ( :) + temp_y * U_temp ( :)

! ------------------------------
! Intermediate HLL wave in this face. Integrate "U_Star_Avg" etc. over
! [ L_1DRS ( i, 1), L_1DRS ( i, NWAVES_1DRS)].

      U_temp ( :) = F_Star_Avg ( i, :) * neta_dot_t ( i)
     1            + G_Star_Avg ( i, :) * neta_dot_n ( i)
     1            - U_Star_Avg ( i, :) * S_Neta_Multid ( i)

! ----------

      l1 = L_1DRS ( i, 1)
      l2 = L_1DRS ( i, NWAVES_1DRS)

! ----------

      temp_0 = l2 - l1

      HLL_PanelInt_U0 ( :) = HLL_PanelInt_U0 ( :) + temp_0 * U_temp ( :)

! ----------
! "temp_x, temp_y" are integrals of "x" and "y" along this segment.

      temp_x = x ( i) * ( l2 - l1)
     1       - neta_y ( i) * 0.5 * ( l2**2 - l1**2)

      HLL_PanelInt_Ux ( :) = HLL_PanelInt_Ux ( :) + temp_x * U_temp ( :)

! ----------

      temp_y = y ( i) * ( l2 - l1)
     1       + neta_x ( i) * 0.5 * ( l2**2 - l1**2)

      HLL_PanelInt_Uy ( :) = HLL_PanelInt_Uy ( :) + temp_y * U_temp ( :)

! ------------------------------

      END DO  ! End of "DO i =" loop.Looping over the side panels of wave model.

! ------------------------------------------------------------

! This is eqn. (28) of B14. It gives us "HLL_U_Bar ( :)", the lowest order
! state.

      tempaa = 0.5 / Int_Momnt_U ( 1)

      HLL_U_Bar ( :) = - tempaa * HLL_PanelInt_U0 ( :)

! ------------------------------

! This is eqn. (33) of B14. It gives us "HLL_F_Bar ( :), the lowest order
! x-flux.

      tempaa = 1.0 / Int_Momnt_U ( 1)

      HLL_F_Bar ( :) = Vlx_Centroid * HLL_U_Bar ( :)
     1               + tempaa * HLL_PanelInt_Ux ( :)

! ------------------------------

! This is eqn. (34) of B14. It gives us "HLL_G_Bar ( :)", the lowest order
! y-flux.

      tempaa = 1.0 / Int_Momnt_U ( 1)

      HLL_G_Bar ( :) = Vly_Centroid * HLL_U_Bar ( :)
     1               + tempaa * HLL_PanelInt_Uy ( :)

! ------------------------------

      IF ( ( BUILD_SUBSTRUCTURE /= 1) .AND. ( NWAVES_1DRS == 2) ) THEN
! If we get here, we are truly doing an HLL RS without substructure. So we
! are done.

! ----------

      U_Star ( :) = HLL_U_Bar ( :)

      F_Star ( :) = HLL_F_Bar ( :)

      G_Star ( :) = HLL_G_Bar ( :)

      U_Lin_Var_x ( :) = 0.0
      U_Lin_Var_y ( :) = 0.0

      U_Lim_Lin_Var_x ( :) = 0.0
      U_Lim_Lin_Var_y ( :) = 0.0

      Linearize_Around_U_Star = 0

! ----------

      Gamma_Star = 0.0

      DO i = 1, nzone_vrtx
      Gamma_Star = Gamma_Star + gamma_eos ( i)
      END DO

      Gamma_Star = Gamma_Star / nzone_vrtx


      Rho_Star = AMAX1 ( U_Star ( 1), rhofloor)

      ob_Rho = 1.0 / Rho_Star

      Vlx_Star = ob_Rho * U_Star ( 2)
      Vly_Star = ob_Rho * U_Star ( 3)
      Vlz_Star = ob_Rho * U_Star ( 4)

      Eng_Tot_Star = U_Star ( 5)

      Bfx_Star = U_Star ( 6)
      Bfy_Star = U_Star ( 7)
      Bfz_Star = U_Star ( 8)


      V2_Star = Vlx_Star**2 + Vly_Star**2 + Vlz_Star**2

      B2_Star = Bfx_Star**2 + Bfy_Star**2 + Bfz_Star**2

      V_Dot_B_Star = Vlx_Star * Bfx_Star + Vly_Star * Bfy_Star
     1             + Vlz_Star * Bfz_Star


      Prs_Gas_Star = Eng_Tot_Star - 0.5 * Rho_Star * V2_Star
     1             - ob_8pi * B2_Star
      Prs_Gas_Star = Prs_Gas_Star * ( Gamma_Star - 1.0)
      Prs_Gas_Star = AMAX1 ( Prs_Gas_Star, prsfloor)

      Prs_Tot_Star = Prs_Gas_Star + ob_8pi * B2_Star

! ----------

      RETURN  ! For multid. HLL RS w/o sub-structure, we are ready to return.

! ----------

      END IF
! End of "IF ( ( BUILD_SUBSTRUCTURE /= 1) .AND. ( NWAVES_1DRS == 2) )".

! ------------------------------

#if ( 0 == 1)
      WRITE ( 6, *)"----------"
      WRITE ( 6, *)"Writing HLL_U_Bar, HLL_F_Bar and HLL_G_Bar:"
      WRITE ( 6, *)"HLL_U_Bar = ", HLL_U_Bar ( :)
      WRITE ( 6, *)"HLL_F_Bar = ", HLL_F_Bar ( :)
      WRITE ( 6, *)"HLL_G_Bar = ", HLL_G_Bar ( :)
#endif

! ----------------------------------------------------------------------
! STEP VII) Build the panel integral for "PanelInt_U0", i.e. zeroth
! moment over the wave model. Use it to build "U_Star". Also build
! "PanelInt_Ux" and "PanelInt_Uy", the x and y-moments. Use them
! to build "F_Star" and "G_Star" respectively.

! Then use "U_Star" and the states at the boundary of the wave model to build
! the full gradients "U_Lin_Var_x, U_Lin_Var_y". This is done as a
! Galerkin projection.

! This is the most important step because it builds the strongly-interacting
! state as well as its linear profile. The linear profile coming out of this
! step already has some reasonable limiting done to it if HLLI is used.

! ----------------------------------------------------------------------

! First build "U_Star", "F_Star" and "G_Star".
! For multidimensional HLLEM, we just start with the multidimensional HLL.
! Later we will add additional sub-structure to it.

      U_Star ( :) = HLL_U_Bar ( :)

      F_Star ( :) = HLL_F_Bar ( :)

      G_Star ( :) = HLL_G_Bar ( :)


#if ( 0 == 1)

      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)"After Panel Integral:"
      WRITE ( 6, *)"U_Star ( :) = ", U_Star ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"F_Star ( :) = ", F_Star ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"G_Star ( :) = ", G_Star ( :)
      WRITE ( 6, *)" "

#endif

! ------------------------------------------------------------
! Build the full gradients for "U_Lin_Var_x, U_Lin_Var_y" using Galerkin
! projection. See eqns. (5.8) and (5.9) of Betal16. For Cartesian
! meshes, these become eqns. (3.14) and (3.15) of Betal16.

! Looping over the side panels of the wave model. Because this is HLL,
! there are only three constant states in each side panel.
! ------------------------------------------------------------

      U_Lin_Var_x ( :) = 0.0
      U_Lin_Var_y ( :) = 0.0

      Mass_Matrix_Lin_Var_x = 0.0
      Mass_Matrix_Lin_Var_y = 0.0

! ------------------------------------------------------------

      DO i = 1, nzone_vrtx

! ------------------------------
! First state in this face. Integrate "U_Minus" etc. over
! [ L_Tau_Min ( i), L_1DRS ( i, 1)].

        U_temp ( :) = U_Minus ( i, :) - U_Star ( :)

! ----------

        l1 = L_Tau_Min ( i)
        l2 = L_1DRS ( i, 1)

! ----------

        Mass_Matrix_Lin_Var_x = Mass_Matrix_Lin_Var_x
     1 + ABS ( neta_x ( i))
     1                  * ( ( ( l2**3 - l1**3 ) * neta_y ( i)**2 ) / 3.0
     1                    - ( l2**2 - l1**2 ) * neta_y ( i) * x ( i)
     1                    + ( l2 - l1 ) * x ( i)**2 )

        Mass_Matrix_Lin_Var_y = Mass_Matrix_Lin_Var_y
     1 + ABS ( neta_y ( i))
     1                  * ( ( ( l2**3 - l1**3 ) * neta_x ( i)**2 ) / 3.0
     1                    - ( l2**2 - l1**2 ) * neta_x ( i) * y ( i)
     1                    + ( l2 - l1 ) * y ( i)**2 )

! ----------

        temp_x = x ( i) * ( l2 - l1)
     1         - neta_y ( i) * 0.5 * ( l2**2 - l1**2)

        temp_x = temp_x * ABS ( neta_x ( i))

        U_Lin_Var_x ( :) = U_Lin_Var_x ( :) + temp_x * U_temp ( :)

! ----------

        temp_y = y ( i) * ( l2 - l1)
     1         + neta_x ( i) * 0.5 * ( l2**2 - l1**2)

        temp_y = temp_y * ABS ( neta_y ( i))

        U_Lin_Var_y ( :) = U_Lin_Var_y ( :) + temp_y * U_temp ( :)

! ------------------------------
! Last state in this face. Integrate "U_Plus" etc. over
! [ L_1DRS ( i, NWAVES_1DRS), L_Tau_Max ( i)].

        U_temp ( :) = U_Plus ( i, :) - U_Star ( :)

! ----------

        l1 = L_1DRS ( i, NWAVES_1DRS)
        l2 = L_Tau_Max ( i)

! ----------

        Mass_Matrix_Lin_Var_x = Mass_Matrix_Lin_Var_x
     1 + ABS ( neta_x ( i))
     1                  * ( ( ( l2**3 - l1**3 ) * neta_y ( i)**2 ) / 3.0
     1                    - ( l2**2 - l1**2 ) * neta_y ( i) * x ( i)
     1                    + ( l2 - l1 ) * x ( i)**2 )

        Mass_Matrix_Lin_Var_y = Mass_Matrix_Lin_Var_y
     1 + ABS ( neta_y ( i))
     1                  * ( ( ( l2**3 - l1**3 ) * neta_x ( i)**2 ) / 3.0
     1                    - ( l2**2 - l1**2 ) * neta_x ( i) * y ( i)
     1                    + ( l2 - l1 ) * y ( i)**2 )

! ----------

        temp_x = x ( i) * ( l2 - l1)
     1         - neta_y ( i) * 0.5 * ( l2**2 - l1**2)

        temp_x = temp_x * ABS ( neta_x ( i))

        U_Lin_Var_x ( :) = U_Lin_Var_x ( :) + temp_x * U_temp ( :)

! ----------

        temp_y = y ( i) * ( l2 - l1)
     1         + neta_x ( i) * 0.5 * ( l2**2 - l1**2)

        temp_y = temp_y * ABS ( neta_y ( i))

        U_Lin_Var_y ( :) = U_Lin_Var_y ( :) + temp_y * U_temp ( :)

! ------------------------------
! Intermediate HLL wave in this face. Integrate "U_Star_Avg" etc. over
! [ L_1DRS ( i, 1), L_1DRS ( i, NWAVES_1DRS)].

        U_temp ( :) = U_Star_Avg ( i, :) - U_Star ( :)

! ----------

        l1 = L_1DRS ( i, 1)
        l2 = L_1DRS ( i, NWAVES_1DRS)

! ----------

        Mass_Matrix_Lin_Var_x = Mass_Matrix_Lin_Var_x
     1 + ABS ( neta_x ( i))
     1                  * ( ( ( l2**3 - l1**3 ) * neta_y ( i)**2 ) / 3.0
     1                    - ( l2**2 - l1**2 ) * neta_y ( i) * x ( i)
     1                    + ( l2 - l1 ) * x ( i)**2 )

        Mass_Matrix_Lin_Var_y = Mass_Matrix_Lin_Var_y
     1 + ABS ( neta_y ( i))
     1                  * ( ( ( l2**3 - l1**3 ) * neta_x ( i)**2 ) / 3.0
     1                    - ( l2**2 - l1**2 ) * neta_x ( i) * y ( i)
     1                    + ( l2 - l1 ) * y ( i)**2 )

! ----------

        temp_x = x ( i) * ( l2 - l1)
     1         - neta_y ( i) * 0.5 * ( l2**2 - l1**2)

        temp_x = temp_x * ABS ( neta_x ( i))

        U_Lin_Var_x ( :) = U_Lin_Var_x ( :) + temp_x * U_temp ( :)

! ----------

        temp_y = y ( i) * ( l2 - l1)
     1         + neta_x ( i) * 0.5 * ( l2**2 - l1**2)

        temp_y = temp_y * ABS ( neta_y ( i))

        U_Lin_Var_y ( :) = U_Lin_Var_y ( :) + temp_y * U_temp ( :)

! ------------------------------

      END DO  ! End of "DO i =" loop.Looping over the side panels of wave model.

! ------------------------------------------------------------
! This completes the two sets of gradients for "U_Lin_Var_x, U_Lin_Var_y" using
! Galerkin projection. Then apply the limiter to the two choices.

      tempaa = 1.0 / Mass_Matrix_Lin_Var_x
      U_Lin_Var_x ( :) = U_Lin_Var_x ( :) * tempaa

      tempab = 1.0 / Mass_Matrix_Lin_Var_y
      U_Lin_Var_y ( :) = U_Lin_Var_y ( :) * tempab

! ------------------------------------------------------------

! Evaluate "[ S_L, S_R] x [ S_D, S_U]" :-
! A general approach consists of looking at the moments
! given in "Int_Momnt_U ( :)".
! This is done below using eqns. (4.6) and (4.7) of Betal16.
! For Cartesian meshes, the two approaches are identical.

! Extent of wave model in the x-direction.
      tempaa = SQRT ( 12.0 * Int_Momnt_U ( 4) / Int_Momnt_U ( 1) )

! Extent of wave model in the y-direction.
      tempab = SQRT ( 12.0 * Int_Momnt_U ( 6) / Int_Momnt_U ( 1) )

      S_L = - 0.5 * tempaa
      S_R = 0.5 * tempaa

      S_D = - 0.5 * tempab
      S_U = 0.5 * tempab

! ----------
! Now that the extent of the wave model has been found in both directions,
! shift it using the centroid of the wave model, stored in
! "(Vlx_Centroid, Vly_Centroid)".

      S_L = S_L + Vlx_Centroid
      S_R = S_R + Vlx_Centroid

      S_D = S_D + Vly_Centroid
      S_U = S_U + Vly_Centroid

! Since this is truly the subsonic case, we enforce coverage of the origin.
! While this will always be ensured in the Cartesian case, there are some
! corner cases associated with polygonal wave models where this needs to 
! be enforced.

! Later on, these extremal speeds will also be expanded to contain all the
! eigenvalues in either direction.

      S_L = AMIN1 ( S_L, - smallnum)
      S_R = AMAX1 ( S_R, smallnum)

      S_D = AMIN1 ( S_D, - smallnum)
      S_U = AMAX1 ( S_U, smallnum)

! ------------------------------------------------------------

#if ( 0 == 1)

      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)"After constructing S_L, S_R, S_D, S_U and slopes "
      WRITE ( 6, *)"S_L, S_R = ", S_L, S_R
      WRITE ( 6, *)"S_D, S_U = ", S_D, S_U
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_Lin_Var_x ( :) = ", U_Lin_Var_x ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_Lin_Var_y ( :) = ", U_Lin_Var_y ( :)
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------
! STEP VIII) From "U_Star" build primitive variables for this state.
! "(Rho_Star, Vlx,y,z_Star, Prs_Gas_Star, Bfx,y,z_Star)".
! We can linearize around this state if it is physical later on.

! ----------------------------------------------------------------------

      Gamma_Star = 0.0

      DO i = 1, nzone_vrtx
      Gamma_Star = Gamma_Star + gamma_eos ( i)
      END DO

      Gamma_Star = Gamma_Star / nzone_vrtx


      Rho_Star = AMAX1 ( U_Star ( 1), rhofloor)

      ob_Rho = 1.0 / Rho_Star

      Vlx_Star = ob_Rho * U_Star ( 2)
      Vly_Star = ob_Rho * U_Star ( 3)
      Vlz_Star = ob_Rho * U_Star ( 4)

      Eng_Tot_Star = U_Star ( 5)

      Bfx_Star = U_Star ( 6)
      Bfy_Star = U_Star ( 7)
      Bfz_Star = U_Star ( 8)


      V2_Star = Vlx_Star**2 + Vly_Star**2 + Vlz_Star**2

      B2_Star = Bfx_Star**2 + Bfy_Star**2 + Bfz_Star**2

      V_Dot_B_Star = Vlx_Star * Bfx_Star + Vly_Star * Bfy_Star
     1             + Vlz_Star * Bfz_Star


      Prs_Gas_Star = Eng_Tot_Star - 0.5 * Rho_Star * V2_Star
     1             - ob_8pi * B2_Star
      Prs_Gas_Star = Prs_Gas_Star * ( Gamma_Star - 1.0)
      Prs_Gas_Star = AMAX1 ( Prs_Gas_Star, prsfloor)

      Prs_Tot_Star = Prs_Gas_Star + ob_8pi * B2_Star

! ------------------------------
! First, find the min and max density and pressure from the neighboring states.
! If the density and pressure lie within this range, we can linearize
! around a physical state. Otherwise, we cannot.

      Rho_Neighbor_Min = V ( 1, 1)
      Rho_Neighbor_Max = V ( 1, 1)

      Prs_Neighbor_Min = V ( 1, 5)
      Prs_Neighbor_Max = V ( 1, 5)

      DO i = 2, nzone_vrtx

      Rho_Neighbor_Min = AMIN1 ( Rho_Neighbor_Min, V ( i, 1) )
      Rho_Neighbor_Max = AMAX1 ( Rho_Neighbor_Max, V ( i, 1) )

      Prs_Neighbor_Min = AMIN1 ( Prs_Neighbor_Min, V ( i, 5) )
      Prs_Neighbor_Max = AMAX1 ( Prs_Neighbor_Max, V ( i, 5) )

      END DO

! ------------------------------
! Insisting that "Linearize_Around_U_Star = 0" all the time would switch
! off the parts of code that evaluate sub-structure.

! We allow a little latitude, in the density and pressure variation
! when we decide when to allow linearization and when not to allow it.

      Linearize_Around_U_Star = 0

      IF ( ( BUILD_SUBSTRUCTURE == 1) .AND.

     1     ( Rho_Neighbor_Min * 0.9 < Rho_Star) .AND.
     1     ( Rho_Star < Rho_Neighbor_Max * 1.1) .AND.

     1     ( Prs_Neighbor_Min * 0.9 < Prs_Gas_Star) .AND.
     1     ( Prs_Gas_Star < Prs_Neighbor_Max * 1.1) .AND.

     1     ( Rho_Star > 2.0 * rhofloor ) .AND.
     1     ( Prs_Gas_Star > 2.0 * prsfloor))

     1                             Linearize_Around_U_Star = 1

! ----------

#if ( 0 == 1)

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)"Linearize_Around_U_Star, Pass 1 = ",
     1              Linearize_Around_U_Star
      WRITE ( 6, *)"Rho_Star, Prs_Gas_Star = ", Rho_Star, Prs_Gas_Star
      WRITE ( 6, *)" "

#endif

! ------------------------------

      IF ( Linearize_Around_U_Star == 0) THEN

! ----------

        IF ( ( Rho_Star <= 2.0 * rhofloor ) .OR.
     1       ( Prs_Gas_Star <= 2.0 * prsfloor) ) THEN

! Density and pressure in the resolved "U_Star ( :)" state are really messed
! up if we get here. If we get here, the safest solution is to use HLL and exit.

! ----------

          U_Star ( :) = HLL_U_Bar ( :)

          F_Star ( :) = HLL_F_Bar ( :)

          G_Star ( :) = HLL_G_Bar ( :)

! -----

          Gamma_Star = 0.0

          DO i = 1, nzone_vrtx
          Gamma_Star = Gamma_Star + gamma_eos ( i)
          END DO

          Gamma_Star = Gamma_Star / nzone_vrtx


          Rho_Star = AMAX1 ( U_Star ( 1), rhofloor)

          ob_Rho = 1.0 / Rho_Star

          Vlx_Star = ob_Rho * U_Star ( 2)
          Vly_Star = ob_Rho * U_Star ( 3)
          Vlz_Star = ob_Rho * U_Star ( 4)

          Eng_Tot_Star = U_Star ( 5)

          Bfx_Star = U_Star ( 6)
          Bfy_Star = U_Star ( 7)
          Bfz_Star = U_Star ( 8)


          V2_Star = Vlx_Star**2 + Vly_Star**2 + Vlz_Star**2

          B2_Star = Bfx_Star**2 + Bfy_Star**2 + Bfz_Star**2

          V_Dot_B_Star = Vlx_Star * Bfx_Star + Vly_Star * Bfy_Star
     1                 + Vlz_Star * Bfz_Star


          Prs_Gas_Star = Eng_Tot_Star - 0.5 * Rho_Star * V2_Star
     1                 - ob_8pi * B2_Star
          Prs_Gas_Star = Prs_Gas_Star * ( Gamma_Star - 1.0)
          Prs_Gas_Star = AMAX1 ( Prs_Gas_Star, prsfloor)

          Prs_Tot_Star = Prs_Gas_Star + ob_8pi * B2_Star

! ----------

        END IF  ! End of "IF ( ( Rho_Star <= 2.0 * rhofloor )".

! ----------

        U_Lin_Var_x ( :) = 0.0
        U_Lin_Var_y ( :) = 0.0

        U_Lim_Lin_Var_x ( :) = 0.0
        U_Lim_Lin_Var_y ( :) = 0.0

        RETURN  ! Returning with multid. HLL flux.

      END IF  ! End of "IF ( Linearize_Around_U_Star == 0)".

! ----------------------------------------------------------------------
! STEP IX) Further restrict the linear variation in "U_Lin_Var_x, U_Lin_Var_y".

! Find the min and max variables "U_Neighbor_Min/Max" by using the states at
! neighboring nodes. Also build "U_Vertex_Min/Max" by using the linearization
! within the wave model.

! Use Barth & Frederickson limiter to limit the slopes. Needs non-zero range
! in order to work. So the ranges are expanded by "smallnum" to make it work.
! ----------------------------------------------------------------------

      IF ( Linearize_Around_U_Star == 1) THEN

! ------------------------------
! First, find min and max from the neighboring (incoming) states; also find
! min and max at the vertices of the wave model. Ensure non-zero range.
! The positiity of density and pressure will serve as an additional check later.

      tempaa = 100.0 * smallnum

      DO j = 1, NFIELD + 1 + NFLUID

      tempab = AMAX1 ( tempaa * ABS ( U_Star ( j)), smallnum)

      U_Neighbor_Min ( j) = U_Star ( j) - tempab
      U_Neighbor_Max ( j) = U_Star ( j) + tempab

      U_Vertex_Min ( j) = U_Star ( j) - tempab
      U_Vertex_Max ( j) = U_Star ( j) + tempab

      END DO


      Rho_Vertex_Min = 1.0e30
      Rho_Vertex_Max = - 1.0e30

      Prs_Vertex_Min = 1.0e30
      Prs_Vertex_Max = - 1.0e30

! ----------

      DO i = 1, nzone_vrtx

! -----

      U_Neighbor_Min ( :) = AMIN1 ( U_Neighbor_Min ( :), U ( i, :))
      U_Neighbor_Max ( :) = AMAX1 ( U_Neighbor_Max ( :), U ( i, :))

! -----

      U_temp ( :) = U_Star ( :) + U_Lin_Var_x ( :) * x ( i)
     1                          + U_Lin_Var_y ( :) * y ( i)

      U_Vertex_Min ( :) = AMIN1 ( U_Vertex_Min ( :), U_temp ( :))
      U_Vertex_Max ( :) = AMAX1 ( U_Vertex_Max ( :), U_temp ( :))

! -----

      Rho_Vertex_Min = AMIN1 ( Rho_Vertex_Min, U_temp ( 1))
      Rho_Vertex_Max = AMAX1 ( Rho_Vertex_Max, U_temp ( 1))

      tempaa = U_temp ( 5)
     1       - 0.5 * ( U_temp ( 2)**2 + U_temp ( 3)**2 + U_temp ( 4)**2)
     1             / AMAX1 ( U_temp ( 1), rhofloor)
     1       - ob_8pi * ( U_temp ( 6)**2 + U_temp ( 7)**2
     1                  + U_temp ( 8)**2)
      tempaa = ( Gamma_Star - 1.0) * tempaa

      Prs_Vertex_Min = AMIN1 ( Prs_Vertex_Min, tempaa)
      Prs_Vertex_Max = AMAX1 ( Prs_Vertex_Max, tempaa)

! -----

      END DO

! ----------
! Second, limit the slopes in component-by-component fashion.
! Keep track of the most extreme limiting that we do. If we have done
! a lot of limiting, it is safest to go back to the original fluxes
! without any sub-structure.

      min_limiter = 1.0

      DO j = 1, NFIELD + 1 + NFLUID

      tempab = ( U_Neighbor_Max ( j) - U_Star ( j) )
     1       / ( U_Vertex_Max ( j) - U_Star ( j))

      tempac = ( U_Star ( j) - U_Neighbor_Min ( j) )
     1       / ( U_Star ( j) - U_Vertex_Min ( j))

      tempaa = AMAX1 ( 0.0, AMIN1 ( 1.0,
     1                 psi_limiter * AMIN1 ( tempab, tempac) ) )

      min_limiter = AMIN1 ( min_limiter, tempaa)

! This is the only place where the profile is being reduced in this code
! (though it can be zeroed elsewhere). For idealized problems
! (like ideal contacts or shear/Alfven waves), it may be ok to comment
! out the next two lines.
      U_Lin_Var_x ( j) = tempaa * U_Lin_Var_x ( j)
      U_Lin_Var_y ( j) = tempaa * U_Lin_Var_y ( j)

      END DO

! ----------
! Third, (if too much limiting is needed, or) if the densities and pressures
! at the vertices are too low, then it is best to bail.
! In that case, the HLL-type solution is the safest one, so we revert to it.

! For constant states, i.e. for "Rho_Neighbor_Max == Rho_Neighbor_Min" situation
! or for "Prs_Neighbor_Max == Prs_Neighbor_Min" situation, the logic still allows
! a profile to develop.

! Since use of "largest_min_limiter" is dicey, we actually don't use it.

! ----------

! More conservative choices adopted here.
      tempaa = 1.1 * Rho_Neighbor_Max
      tempab = 0.9 * Rho_Neighbor_Min

      tempac = 1.1 * Prs_Neighbor_Max
      tempad = 0.9 * Prs_Neighbor_Min

! ----------

#if ( 0 == 1)

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)"Making sure profile does not hit Rho, Prs floor"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Rho_Vertex_Min, tempab = ", Rho_Vertex_Min, tempab
      WRITE ( 6, *)"Prs_Vertex_Min, tempad = ", Prs_Vertex_Min, tempad
      WRITE ( 6, *)" "

#endif

! ----------

      IF (

!     1     ( min_limiter < largest_min_limiter) .OR.

     1     ( Rho_Vertex_Min <= tempab) .OR.  ! Ok to comment out for idealized tests.
     1     ( Prs_Vertex_Min <= tempad) .OR.  ! Ok to comment out for idealized tests.

     1     ( Rho_Vertex_Min <= rhofloor) .OR.
     1     ( Prs_Vertex_Min <= prsfloor) ) THEN

        U_Star ( :) = HLL_U_Bar ( :)
        F_Star ( :) = HLL_F_Bar ( :)
        G_Star ( :) = HLL_G_Bar ( :)

        Linearize_Around_U_Star = 0

        U_Lin_Var_x ( :) = 0.0
        U_Lin_Var_y ( :) = 0.0

        U_Lim_Lin_Var_x ( :) = 0.0
        U_Lim_Lin_Var_y ( :) = 0.0

        RETURN  ! Return with multid. HLL fluxes.

      END IF

! ------------------------------

      END IF  ! End of "IF ( Linearize_Around_U_Star == 1)".

! ------------------------------

#if ( 0 == 1)

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)"Linearize_Around_U_Star, Pass 2 after multid Lim= ",
     1              Linearize_Around_U_Star
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_Lin_Var_x ( :) = ", U_Lin_Var_x ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_Lin_Var_y ( :) = ", U_Lin_Var_y ( :)
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------
! STEP X) From the "U_Star" state, we build eigenvalues,
! "eigenval_x, eigenval_y" and eigenvectors,
! "lf_eigenvec_x, rt_eigenvec_x, lf_eigenvec_y, rt_eigenvec_y".
! These eigenvectors only contain the linearly degenerate waves.

! Use these eigenvalues and eigenvectors to project the full gradients
! from "U_Lin_Var_x, U_Lin_Var_y" so that only the linearly degenerate part
! is retained in "U_Lim_Lin_Var_x, U_Lim_Lin_Var_y".

! ----------------------------------------------------------------------

      IF ( Linearize_Around_U_Star == 1) THEN

! ------------------------------------------------------------

      iprimitive = 0  ! Asking for eigenvectors in conserved variables.

      gasconst = 1.0  ! Gas constant is irrelevant for these eigenvectors.

      Xbar = 0.0  ! These are pure eigenvectors, so "Xbar = 0.0".

! ------------------------------

      idirection = 1  ! Asking for eigenvectors in the x-direction.

      debug_stop = 0

#if ( COMPLETE_EIGENSYSTEM == 0)
      CALL EIGENVECMHD_PTWISE_LINDEG
#else
      CALL EIGENVECMHD_PTWISE_8
#endif

     1             ( iprimitive, idirection,

     1               Gamma_Star, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,

     1               Rho_Star, Prs_Gas_Star,
     1               Vlx_Star, Vly_Star, Vlz_Star,
     1               Bfx_Star, Bfy_Star, Bfz_Star, Xbar,
     1               Specie_Star,

     1               eigenval_x,
     1               lf_eigenvec_x, rt_eigenvec_x )

! --------------------
! Expand "S_L, S_R" based on "eigenval_x ( :)". We assume that all
! x-directional eigenvalues are provided by the previous subroutine.

      S_R = AMAX1 ( S_R, eigenval_x ( 7))
      S_L = AMIN1 ( S_L, eigenval_x ( 1))

! ------------------------------------------------------------

      idirection = 2  ! Asking for eigenvectors in the y-direction.

      debug_stop = 0

#if ( COMPLETE_EIGENSYSTEM == 0)
      CALL EIGENVECMHD_PTWISE_LINDEG
#else
      CALL EIGENVECMHD_PTWISE_8
#endif

     1             ( iprimitive, idirection,

     1               Gamma_Star, smallnum, prsfloor, rhofloor, pi,
     1               debug_stop,

     1               gasconst,

     1               Rho_Star, Prs_Gas_Star,
     1               Vlx_Star, Vly_Star, Vlz_Star,
     1               Bfx_Star, Bfy_Star, Bfz_Star, Xbar,
     1               Specie_Star,

     1               eigenval_y,
     1               lf_eigenvec_y, rt_eigenvec_y )

! --------------------
! Expand "S_D, S_U" based on "eigenval_y ( :)". We assume that all
! y-directional eigenvalues are provided by the previous subroutine.

      S_U = AMAX1 ( S_U, eigenval_y ( 7))
      S_D = AMIN1 ( S_D, eigenval_y ( 1))

! ------------------------------------------------------------

! Zero out eigenweights and "delta" for all the waves in both directions.

      eigenwt_x ( :) = 0.0
      eigenwt_y ( :) = 0.0

      delta_x ( :) = 0.0
      delta_y ( :) = 0.0

! Characteristically projected linear variations. They contain only the
! contribution from the linearly degenerate waves.

      U_Lim_Lin_Var_x ( :) = 0.0
      U_Lim_Lin_Var_y ( :) = 0.0

! ----------------------------------------

      DO i = 1, 7  ! Looping over waves.

! ----------------------------------------

#if ( COMPLETE_EIGENSYSTEM == 0)

! Entropy wave and Alfven waves only.
      IF ( ( i == 1) .OR. ( i == 3) .OR. ( i == 5) .OR. ( i == 7)) CYCLE

#elif ( COMPLETE_EIGENSYSTEM == 1)

! Entropy wave, Alfven waves and slow magnetosonic waves.
! The Alfven + slow msonics carry most of the shear in the flow.
      IF ( ( i == 1) .OR. ( i == 7)) CYCLE

#else

! Fast magnetosonic waves also included. The fast msonics carry most of the
! pressure jump.

#endif

! ----------------------------------------
! Build eigenweights of linearly degenerate fields.

        DO j = 1, NFIELD + NFLUID + 1

          eigenwt_x ( i) = eigenwt_x ( i)
     1                   + lf_eigenvec_x ( i, j) * U_Lin_Var_x ( j)

          eigenwt_y ( i) = eigenwt_y ( i)
     1                   + lf_eigenvec_y ( i, j) * U_Lin_Var_y ( j)

        END DO

! ----------------------------------------
! Improve eigenweights of linearly degenerate fields. Also make sure that
! eigenvalues do not exceed the boundaries of the wave model.

        eigenval_x ( i) = AMIN1 ( AMAX1 ( eigenval_x ( i), S_L), S_R)

        eigenval_y ( i) = AMIN1 ( AMAX1 ( eigenval_y ( i), S_D), S_U)

! ------------------------------

        IF ( iHLLEM_choice == 1) THEN

! ------------------------------

! Choose "delta_x, delta_y" from eqn. (28) of DB16.

        delta_x ( i) = 1.0
     1          - AMIN1 ( eigenval_x ( i), 0.0) / ( S_L - smallnum)
     1          - AMAX1 ( eigenval_x ( i), 0.0) / ( S_R + smallnum)

        delta_y ( i) = 1.0
     1          - AMIN1 ( eigenval_y ( i), 0.0) / ( S_D - smallnum)
     1          - AMAX1 ( eigenval_y ( i), 0.0) / ( S_U + smallnum)

! ------------------------------

        ELSE

! ------------------------------

! Choose "delta_x, delta_y" eqns. (3.20) and (3.21) from Betal16.

        delta_x ( i) = - 2.0 * S_R * S_L / ( S_R - S_L)**2  ! Baseline value.

        tempaa = 0.5 * ( S_R - S_L)**2 - eigenval_x ( i) * ( S_R + S_L)

        tempab = S_R * AMIN1 ( eigenval_x ( i), 0.0)
     1         + S_L * AMAX1 ( eigenval_x ( i), 0.0) - S_R * S_L

        IF ( tempaa > 1.0e-6 * ( S_R - S_L)**2 ) delta_x ( i) =
     1       AMAX1 ( 0.0, AMIN1 ( delta_x ( i), tempab / tempaa) )

! -----

        delta_y ( i) = - 2.0 * S_U * S_D / ( S_U - S_D)**2  ! Baseline value.

        tempaa = 0.5 * ( S_U - S_D)**2 - eigenval_y ( i) * ( S_U + S_D)

        tempab = S_U * AMIN1 ( eigenval_y ( i), 0.0)
     1         + S_D * AMAX1 ( eigenval_y ( i), 0.0) - S_U * S_D

        IF ( tempaa > 1.0e-6 * ( S_U - S_D)**2 ) delta_y ( i) =
     1       AMAX1 ( 0.0, AMIN1 ( delta_y ( i), tempab / tempaa) )

! ------------------------------

        END IF  ! End of "IF ( iHLLEM_choice == 1)".

! ------------------------------

        tempac = 1.0 - AMAX1 ( blend_HLL, blend_MDHLL)

        delta_x ( i) = tempac * delta_x ( i)
        delta_y ( i) = tempac * delta_y ( i)

! ----------------------------------------
! Build gradient that only has contributions from linearly degenerate fields.
! These are eqns. (3.16) and (3.17) of Betal16.
! These could also be eqns. (5.10) and (5.11) of Betal16 because those
! equations are more suitable for use with divided differences.

          tempaa = 2.0 * delta_x ( i) * eigenwt_x ( i)

          tempab = 2.0 * delta_y ( i) * eigenwt_y ( i)

! ----------

          U_Lim_Lin_Var_x ( :) = U_Lim_Lin_Var_x ( :)
     1                         + tempaa * rt_eigenvec_x ( :, i)

          U_Lim_Lin_Var_y ( :) = U_Lim_Lin_Var_y ( :)
     1                         + tempab * rt_eigenvec_y ( :, i)

! ----------------------------------------

      END DO  ! End of "DO i = 1, 7"; I.e., Looping over waves.

! ------------------------------------------------------------

      END IF  ! End of "IF ( Linearize_Around_U_Star == 1)".

! ------------------------------

#if ( 0 == 1)

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)"Linearize_Around_U_Star, Pass 3, after char proj= ",
     1              Linearize_Around_U_Star
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_Lim_Lin_Var_x ( :) = ", U_Lim_Lin_Var_x ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U_Lim_Lin_Var_y ( :) = ", U_Lim_Lin_Var_y ( :)
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------
! STEP XI) "U_Star" is fully built here. The constant parts of
! "F_Star" and "G_Star" also take on contributions because of the
! substructure and those contributions are added in here.

! It would be prudent to test whether "U_Star ( :)" corresponds to
! a physical state. If it does not have a physical state, just use
! "HLL_U_Bar ( :), HLL_F_Bar ( :), HLL_G_Bar ( :)".

! At the end of this step and the next, we return "U_Star, F_Star and G_Star".

! ----------------------------------------------------------------------

! If we are working with a moving mesh with nodal velocity (Vlx_ALE, Vly_ALE)
! then all that matters is the relative velocity of the node w.r.t.
! ( Vlx_Centroid, Vly_Centroid).
! See eqns. (5.12), (5.13) and (5.14) of Betal16.
! Also see eqns. (5.5) and (5.6) of Betal16.

      Sx_h = Vlx_ALE - Vlx_Centroid
      Sy_h = Vly_ALE - Vly_Centroid

! ------------------------------

      IF ( Linearize_Around_U_Star == 1) THEN

! ------------------------------

! Build "U_Lin_Var ( :)" which is the linear variation  that is added to
! "U_Star ( :)". See eqn. (5.17) of Betal16.

      U_Lin_Var ( :) = U_Lim_Lin_Var_x ( :) * Sx_h
     1               + U_Lim_Lin_Var_y ( :) * Sy_h

! This is our final "U_Star" with sub-structure added in.

      U_Star ( :) = U_Star ( :) + U_Lin_Var ( :)

! ----------

! Also build "U_Lin_Var_x ( :), U_Lin_Var_y ( :)" which carry the partial
! contributions to the x- and y-fluxes.

      U_Lin_Var_x ( :) = U_Lim_Lin_Var_x ( :) * Sx_h

      U_Lin_Var_y ( :) = U_Lim_Lin_Var_y ( :) * Sy_h

! ------------------------------

! Eqns. (5.5) and (5.6) of Betal16 also dictate that the constant parts of
! "F_Star ( :), G_Star ( :)" should take on extra terms. Provide
! those extra terms here.

      tempaa = 3.0 * Int_Momnt_U ( 4) / Int_Momnt_U ( 1)

      tempab = 3.0 * Int_Momnt_U ( 6) / Int_Momnt_U ( 1)

      F_Star ( :) = F_Star ( :) + tempaa * U_Lim_Lin_Var_x ( :)

      G_Star ( :) = G_Star ( :) + tempab * U_Lim_Lin_Var_y ( :)

! ------------------------------
! Build primitive variables from our final expression for "U_Star ( :)".

      Gamma_Star = 0.0

      DO i = 1, nzone_vrtx
      Gamma_Star = Gamma_Star + gamma_eos ( i)
      END DO

      Gamma_Star = Gamma_Star / nzone_vrtx


      Rho_Star = AMAX1 ( U_Star ( 1), rhofloor)

      ob_Rho = 1.0 / Rho_Star

      Vlx_Star = ob_Rho * U_Star ( 2)
      Vly_Star = ob_Rho * U_Star ( 3)
      Vlz_Star = ob_Rho * U_Star ( 4)

      Eng_Tot_Star = U_Star ( 5)

      Bfx_Star = U_Star ( 6)
      Bfy_Star = U_Star ( 7)
      Bfz_Star = U_Star ( 8)


      V2_Star = Vlx_Star**2 + Vly_Star**2 + Vlz_Star**2

      B2_Star = Bfx_Star**2 + Bfy_Star**2 + Bfz_Star**2

      V_Dot_B_Star = Vlx_Star * Bfx_Star + Vly_Star * Bfy_Star
     1             + Vlz_Star * Bfz_Star


      Prs_Gas_Star = Eng_Tot_Star - 0.5 * Rho_Star * V2_Star
     1             - ob_8pi * B2_Star
      Prs_Gas_Star = Prs_Gas_Star * ( Gamma_Star - 1.0)
      Prs_Gas_Star = AMAX1 ( Prs_Gas_Star, prsfloor)

      Prs_Tot_Star = Prs_Gas_Star + ob_8pi * B2_Star

      Specie_Star ( 0) = 0.0

#if ( NFLUID >= 1)
      DO i = 1, NFLUID
      Specie_Star ( i) = ob_Rho * U_Star ( NFIELD + 1 + i)
      END DO
#endif

! --------------------
! If the resolved state density or pressure are below floor values, it still
! makes sense to use HLL states and fluxes and exit.

      IF ( ( Rho_Star <= rhofloor) .OR.
     1     ( Prs_Gas_Star <= prsfloor) ) THEN

        U_Star ( :) = HLL_U_Bar ( :)

        F_Star ( :) = HLL_F_Bar ( :)
        G_Star ( :) = HLL_G_Bar ( :)

        Linearize_Around_U_Star = 0

        U_Lin_Var_x ( :) = 0.0
        U_Lin_Var_y ( :) = 0.0

        U_Lim_Lin_Var_x ( :) = 0.0
        U_Lim_Lin_Var_y ( :) = 0.0

        RETURN

      END IF

! ------------------------------

      END IF  ! End of "IF ( Linearize_Around_U_Star == 1)".

! ------------------------------

#if ( 0 == 1)

      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)"Before multiplication with Abar, Bbar:"
      WRITE ( 6, *)"U_Star ( :) = ", U_Star ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"F_Star ( :) = ", F_Star ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"G_Star ( :) = ", G_Star ( :)
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------
! STEP XII) Build "Abar_InnerWaves, Bbar_InnerWaves". These are the only
! waves we want to have propagating within the strongly interacting state.
! They are obtained by setting the remaining eigenvalues and eigenvectors to
! zero.

! These two matrices act like projection operators. They project the linear
! variation of the fluxes on to the part that carries fluctuations in the
! linearly degenerate waves. The nonlinear, i.e. magnetosonic, waves
! are projected out.

! See eqns. (5.13) and (5.14) of Betal16. Also see eqn. (3.22) of Betal16
! for information about constructing the characteristic matrices.

! ----------------------------------------------------------------------

      IF ( Linearize_Around_U_Star == 1) THEN

! ------------------------------------------------------------
! Now use "eigenval_x, lf_eigenvec_x, rt_eigenvec_x" to build
! "Abar_InnerWaves". This is eqn. (3.22) of Betal16.

      Abar_InnerWaves ( :, :) = 0.0

! --------------------
! Project the waves to build a (syncopated) characteristic matrix in the x-direction.

      DO k = 1, 7

! ----------

#if ( COMPLETE_EIGENSYSTEM == 0)

! Entropy wave and Alfven waves only.
      IF ( ( k == 1) .OR. ( k == 3) .OR. ( k == 5) .OR. ( k == 7)) CYCLE

#elif ( COMPLETE_EIGENSYSTEM == 1)

! Entropy wave, Alfven waves and slow magnetosonic waves. The Alfven + slow msonics
! cary most of the shear in the flow.
      IF ( ( k == 1) .OR. ( k == 7)) CYCLE

#else

! Fast magnetosonic waves also included. The fast msonics carry most of the
! pressure jump.

#endif

! ----------

        DO i = 1, NFIELD + NFLUID + 1

          DO j = 1, NFIELD + NFLUID + 1

          Abar_InnerWaves ( i, j) = Abar_InnerWaves ( i, j)
     1 + rt_eigenvec_x ( i, k) * eigenval_x ( k) * lf_eigenvec_x ( k, j)

          END DO

        END DO

! ----------

      END DO

! ------------------------------------------------------------
! Now use "eigenval_y, lf_eigenvec_y, rt_eigenvec_y" to build
! "Bbar_InnerWaves". This is eqn. (3.22) of Betal16.

      Bbar_InnerWaves ( :, :) = 0.0

! --------------------
! Project the waves to build a (syncopated) characteristic matrix in the y-direction.

      DO k = 1, 7

! ----------

#if ( COMPLETE_EIGENSYSTEM == 0)

! Entropy wave and Alfven waves only.
      IF ( ( k == 1) .OR. ( k == 3) .OR. ( k == 5) .OR. ( k == 7)) CYCLE

#elif ( COMPLETE_EIGENSYSTEM == 1)

! Entropy wave, Alfven waves and slow magnetosonic waves. The Alfven + slow msonics
! cary most of the shear in the flow.
      IF ( ( k == 1) .OR. ( k == 7)) CYCLE

#else

! Fast magnetosonic waves also included. The fast msonics carry most of the
! pressure jump.

#endif

! ----------

        DO i = 1, NFIELD + NFLUID + 1

          DO j = 1, NFIELD + NFLUID + 1

          Bbar_InnerWaves ( i, j) = Bbar_InnerWaves ( i, j)
     1 + rt_eigenvec_y ( i, k) * eigenval_y ( k) * lf_eigenvec_y ( k, j)

          END DO

        END DO

! ----------

      END DO

! ------------------------------------------------------------
! Now add to fluxes "F_Star ( :), G_Star ( :)" by using the characteristic
! matrices. This entails adding the contribution from "U_Lin_Var_x ( :)" and
! "U_Lin_Var_y ( :)" with the help of the characteristic matrices.

! Notice that "U_Star ( :)" has already been set in the previous step.
! Here we just improve the HLL fluxes to include the additional structure.

! Use eqns. (3.17) and (3.18) of Betal16.
! Alternatively, see eqns. (5.18) and (5.19) of Betal16.

        DO i = 1, NFIELD + 1 + NFLUID

          DO j = 1, NFIELD + 1 + NFLUID

          F_Star ( i) = F_Star ( i)
     1                + Abar_InnerWaves ( i, j) * U_Lin_Var_x ( j)

          G_Star ( i) = G_Star ( i)
     1                + Bbar_InnerWaves ( i, j) * U_Lin_Var_y ( j)

          END DO

        END DO

! ------------------------------------------------------------
! Blend in some HLL flux if some of it is needed for stability's sake.

        tempaa = AMAX1 ( blend_HLL, blend_MDHLL)
        tempab = 1.0 - tempaa

        U_Star ( :) = tempaa * HLL_U_Bar ( :) + tempab * U_Star ( :)

        F_Star ( :) = tempaa * HLL_F_Bar ( :) + tempab * F_Star ( :)

        G_Star ( :) = tempaa * HLL_G_Bar ( :) + tempab * G_Star ( :)

! ----------------------------------------------------------------------

      END IF
! End of "IF ( Linearize_Around_U_Star == 1)"

! ------------------------------

#if ( 0 == 1)

      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)"After multiplication with Abar, Bbar:"
      WRITE ( 6, *)"U_Star ( :) = ", U_Star ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"F_Star ( :) = ", F_Star ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"G_Star ( :) = ", G_Star ( :)
      WRITE ( 6, *)" "

#endif

! ----------------------------------------------------------------------

      END SUBROUTINE MuSIC_MHD_2DRS_1_STRUCT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FLIP_V2U ( nzone_vrtx, idirection, pi, V, U,
     1                      gamma_eos)

! ----------------------------------------------------------------------

! The 2d RS is instituted in the xy-direction with the z-direction
! being the direction of the edge. When the edge is along the
! x-direction (idirection = 1) or y-direction ((idirection = 2) we
! have to rotate the velocity and magnetic field components.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx, idirection, i

      REAL pi, ob_8pi

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) :: V, U

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) :: gamma_eos

      REAL temp_01, temp_02, temp_03

! ----------------------------------------------------------------------
! Build conserved variables "U" from primitive variables "V".

      ob_8pi = 1.0 / ( 8.0 * pi)

      DO i = 1, nzone_vrtx

      U ( i, 1) = V ( i, 1)

      U ( i, 2) = V ( i, 1) * V ( i, 2)
      U ( i, 3) = V ( i, 1) * V ( i, 3)
      U ( i, 4) = V ( i, 1) * V ( i, 4)

      U ( i, 5) = V ( i, 5) / ( gamma_eos ( i) - 1.0)
     1        + 0.5 * V ( i, 1) * ( V ( i, 2)**2 + V ( i, 3)**2
     1                            + V ( i, 4)**2 )
     1        + ob_8pi * ( V ( i, 6)**2 + V ( i, 7)**2 + V ( i, 8)**2)

      U ( i, 6) = V ( i, 6)
      U ( i, 7) = V ( i, 7)
      U ( i, 8) = V ( i, 8)

#if ( NFLUID >= 1)
      U ( i, 9: NFIELD + 1 + NFLUID) = V ( i, 9: NFIELD + 1 + NFLUID)
#endif

      END DO

! ----------------------------------------------------------------------

      IF ( idirection == 1) THEN

! ----------

      DO i = 1, nzone_vrtx

      temp_01 = U ( i, 2)
      temp_02 = U ( i, 3)
      temp_03 = U ( i, 4)

      U ( i, 2) = temp_02
      U ( i, 3) = temp_03
      U ( i, 4) = temp_01

      temp_01 = U ( i, 6)
      temp_02 = U ( i, 7)
      temp_03 = U ( i, 8)

      U ( i, 6) = temp_02
      U ( i, 7) = temp_03
      U ( i, 8) = temp_01

      END DO

! ----------

      ELSE IF ( idirection == 2) THEN

! ----------

      DO i = 1, nzone_vrtx

      temp_01 = U ( i, 2)
      temp_02 = U ( i, 3)
      temp_03 = U ( i, 4)

      U ( i, 2) = temp_03
      U ( i, 3) = temp_01
      U ( i, 4) = temp_02

      temp_01 = U ( i, 6)
      temp_02 = U ( i, 7)
      temp_03 = U ( i, 8)

      U ( i, 6) = temp_03
      U ( i, 7) = temp_01
      U ( i, 8) = temp_02

      END DO

! ----------

      END IF

! ----------------------------------------------------------------------

      END SUBROUTINE FLIP_V2U

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FLIP_U ( nzone_vrtx, idirection, U)

! ----------------------------------------------------------------------

! The 2d RS is instituted in the xy-direction with the z-direction
! being the direction of the edge. When the edge is along the
! x-direction (idirection = 1) or y-direction ((idirection = 2) we
! have to rotate the velocity and magnetic field components.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx, idirection, i

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) :: U

      REAL temp_01, temp_02, temp_03

! ----------------------------------------------------------------------

      IF ( idirection == 1) THEN

! ----------

      DO i = 1, nzone_vrtx

      temp_01 = U ( i, 2)
      temp_02 = U ( i, 3)
      temp_03 = U ( i, 4)

      U ( i, 2) = temp_02
      U ( i, 3) = temp_03
      U ( i, 4) = temp_01

      temp_01 = U ( i, 6)
      temp_02 = U ( i, 7)
      temp_03 = U ( i, 8)

      U ( i, 6) = temp_02
      U ( i, 7) = temp_03
      U ( i, 8) = temp_01

      END DO

! ----------

      ELSE IF ( idirection == 2) THEN

! ----------

      DO i = 1, nzone_vrtx

      temp_01 = U ( i, 2)
      temp_02 = U ( i, 3)
      temp_03 = U ( i, 4)

      U ( i, 2) = temp_03
      U ( i, 3) = temp_01
      U ( i, 4) = temp_02

      temp_01 = U ( i, 6)
      temp_02 = U ( i, 7)
      temp_03 = U ( i, 8)

      U ( i, 6) = temp_03
      U ( i, 7) = temp_01
      U ( i, 8) = temp_02

      END DO

! ----------

      END IF

! ----------------------------------------------------------------------

      END SUBROUTINE FLIP_U

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE UNFLIP_FG ( nzone_vrtx, idirection, F_Star, G_Star,
     1                       F_HLLC_Global_x, G_HLLC_Global_y)

! ----------------------------------------------------------------------

! The 2d RS is instituted in the xy-direction with the z-direction 
! being the direction of the edge. When the edge is along the 
! x-direction (idirection = 1) or y-direction ((idirection = 2) we 
! have to back-rotate the velocity and magnetic field components in
! the fluxes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx, idirection, i

      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: F_Star, G_Star

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) :: 
     1                  F_HLLC_Global_x, G_HLLC_Global_y

      REAL temp_01, temp_02, temp_03

! ----------------------------------------------------------------------

      IF ( idirection == 1) THEN

! ----------

      temp_01 = F_Star ( 2)
      temp_02 = F_Star ( 3)
      temp_03 = F_Star ( 4)

      F_Star ( 2) = temp_03
      F_Star ( 3) = temp_01
      F_Star ( 4) = temp_02

      temp_01 = F_Star ( 6)
      temp_02 = F_Star ( 7)
      temp_03 = F_Star ( 8)

      F_Star ( 6) = temp_03
      F_Star ( 7) = temp_01
      F_Star ( 8) = temp_02

! -----

      temp_01 = G_Star ( 2)
      temp_02 = G_Star ( 3)
      temp_03 = G_Star ( 4)

      G_Star ( 2) = temp_03
      G_Star ( 3) = temp_01
      G_Star ( 4) = temp_02

      temp_01 = G_Star ( 6)
      temp_02 = G_Star ( 7)
      temp_03 = G_Star ( 8)

      G_Star ( 6) = temp_03
      G_Star ( 7) = temp_01
      G_Star ( 8) = temp_02

! -----

      DO i = 1, nzone_vrtx

! -----

      temp_01 = F_HLLC_Global_x ( i, 2)
      temp_02 = F_HLLC_Global_x ( i, 3)
      temp_03 = F_HLLC_Global_x ( i, 4)

      F_HLLC_Global_x ( i, 2) = temp_03
      F_HLLC_Global_x ( i, 3) = temp_01
      F_HLLC_Global_x ( i, 4) = temp_02

      temp_01 = F_HLLC_Global_x ( i, 6)
      temp_02 = F_HLLC_Global_x ( i, 7)
      temp_03 = F_HLLC_Global_x ( i, 8)

      F_HLLC_Global_x ( i, 6) = temp_03
      F_HLLC_Global_x ( i, 7) = temp_01
      F_HLLC_Global_x ( i, 8) = temp_02

! -----

      temp_01 = G_HLLC_Global_y ( i, 2)
      temp_02 = G_HLLC_Global_y ( i, 3)
      temp_03 = G_HLLC_Global_y ( i, 4)

      G_HLLC_Global_y ( i, 2) = temp_03
      G_HLLC_Global_y ( i, 3) = temp_01
      G_HLLC_Global_y ( i, 4) = temp_02

      temp_01 = G_HLLC_Global_y ( i, 6)
      temp_02 = G_HLLC_Global_y ( i, 7)
      temp_03 = G_HLLC_Global_y ( i, 8)

      G_HLLC_Global_y ( i, 6) = temp_03
      G_HLLC_Global_y ( i, 7) = temp_01
      G_HLLC_Global_y ( i, 8) = temp_02

! -----

      END DO

! ----------

      ELSE IF ( idirection == 2) THEN

! ----------

      temp_01 = F_Star ( 2)
      temp_02 = F_Star ( 3)
      temp_03 = F_Star ( 4)

      F_Star ( 2) = temp_02
      F_Star ( 3) = temp_03
      F_Star ( 4) = temp_01

      temp_01 = F_Star ( 6)
      temp_02 = F_Star ( 7)
      temp_03 = F_Star ( 8)

      F_Star ( 6) = temp_02
      F_Star ( 7) = temp_03
      F_Star ( 8) = temp_01

! -----

      temp_01 = G_Star ( 2)
      temp_02 = G_Star ( 3)
      temp_03 = G_Star ( 4)

      G_Star ( 2) = temp_02
      G_Star ( 3) = temp_03
      G_Star ( 4) = temp_01

      temp_01 = G_Star ( 6)
      temp_02 = G_Star ( 7)
      temp_03 = G_Star ( 8)

      G_Star ( 6) = temp_02
      G_Star ( 7) = temp_03
      G_Star ( 8) = temp_01

! -----

      DO i = 1, nzone_vrtx

! -----

      temp_01 = F_HLLC_Global_x ( i, 2)
      temp_02 = F_HLLC_Global_x ( i, 3)
      temp_03 = F_HLLC_Global_x ( i, 4)

      F_HLLC_Global_x ( i, 2) = temp_02
      F_HLLC_Global_x ( i, 3) = temp_03
      F_HLLC_Global_x ( i, 4) = temp_01

      temp_01 = F_HLLC_Global_x ( i, 6)
      temp_02 = F_HLLC_Global_x ( i, 7)
      temp_03 = F_HLLC_Global_x ( i, 8)

      F_HLLC_Global_x ( i, 6) = temp_02
      F_HLLC_Global_x ( i, 7) = temp_03
      F_HLLC_Global_x ( i, 8) = temp_01

! -----

      temp_01 = G_HLLC_Global_y ( i, 2)
      temp_02 = G_HLLC_Global_y ( i, 3)
      temp_03 = G_HLLC_Global_y ( i, 4)

      G_HLLC_Global_y ( i, 2) = temp_02
      G_HLLC_Global_y ( i, 3) = temp_03
      G_HLLC_Global_y ( i, 4) = temp_01

      temp_01 = G_HLLC_Global_y ( i, 6)
      temp_02 = G_HLLC_Global_y ( i, 7)
      temp_03 = G_HLLC_Global_y ( i, 8)

      G_HLLC_Global_y ( i, 6) = temp_02
      G_HLLC_Global_y ( i, 7) = temp_03
      G_HLLC_Global_y ( i, 8) = temp_01

! -----

      END DO

! ----------

      END IF

! ----------------------------------------------------------------------

      END SUBROUTINE UNFLIP_FG

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
 
      SUBROUTINE ONE_D_MHD_PRIM_F_G_WAVSPD ( nzone_vrtx,

     1     smallnum, prsfloor, rhofloor, pi, expand_wavespeeds,

     1     U, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     V, F, G,
     1     S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max)

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Takes in "U, gamma_eos, gamma_sdsp, nx, ny" and uses them to build
! the corresponding primitives "V" and x- and y-fluxes "F, G".

! Here, "nx, ny" are components of the normal vector along the principal
! directions of wave propagation. The vector ( nx ( 1), ny ( 1)) is the
! mesh-aligned unit normal between states U ( 1, :) and U ( 2, :). The
! convention continues cyclically.

! Builds the extremal speeds in both the Tau- and Neta-directions
! and stores them in "S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max".

! Inside the subroutine, we rotate into the local frames to find the
! extremal speeds. The variables and fluxes are always returned in the
! global space.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx

      REAL smallnum, prsfloor, rhofloor, pi, expand_wavespeeds


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Two state variables that come together at a zone edge. We build their
! primitive variables and associated x- and y-fluxes.
     1      U, V, F, G

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices; normal vectors along the direction of wave propagation.
     1      gamma_eos, gamma_sdsp, nx, ny,

! Entropy density.
     1      entropy,

! The extremal speeds in both the Tau- and Neta-directions
     1      S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      INTEGER i, j, ip1

      REAL fourpi, ob_4pi, ob_8pi

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

     1      Rho, Prs, Etot, Vlx, Vly, Vlz,
     1      Bfx, Bfy, Bfz, B2, V_dot_B,

     1      Rho_use, Prs_use, gamma_sdsp_use, gamma_eos_use,
     1      Total_Enthalpy_remi, Total_Enthalpy_repl,
     1      Total_Enthalpy_use,
     1      Vl_Neta_use, Vl_Tau_use, Vlz_use,
     1      Bf_Neta_use, Bf_Tau_use, Bfz_use, B2_use, X_use,

     1      temp_01, temp_02, temp_03, temp_04

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 / fourpi

      IF ( MAX_ZONES_PER_VRTX < nzone_vrtx) THEN
      WRITE ( 6, *)"STOPPING in ONE_D_MHD_PRIM_F_G_WAVSPD because : "
      WRITE ( 6, *)"NZONES_","PER_VRTX < nzone_vrtx",
     1              MAX_ZONES_PER_VRTX, nzone_vrtx
      STOP
      END IF

      IF ( NFIELD /= 7) THEN
      WRITE ( 6, *)"STOPPING in ONE_D_MHD_PRIM_F_G_WAVSPD because : "
      WRITE ( 6, *)"NFIELD /= 7", NFIELD
      STOP
      END IF

! ----------------------------------------------------------------------
! STEP I) Build "V" from "U". Also build the named primitive variables.
! Use them to build the flux vectors "F, G".
! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

! ------------------------------

      Rho ( i) = AMAX1 ( U ( i, 1), rhofloor)
      temp_01 ( i) = 1.0 / Rho ( i)

      Vlx ( i) = U ( i, 2) * temp_01 ( i)
      Vly ( i) = U ( i, 3) * temp_01 ( i)
      Vlz ( i) = U ( i, 4) * temp_01 ( i)

      Etot ( i) = U ( i, 5)

      Bfx ( i) = U ( i, 6)
      Bfy ( i) = U ( i, 7)
      Bfz ( i) = U ( i, 8)

      Prs ( i) = U ( i, 5)
     1            - 0.5 * Rho ( i) * ( Vlx ( i)**2 + Vly ( i)**2
     1                               + Vlz ( i)**2 )
     1            - ob_8pi * ( Bfx ( i)**2 + Bfy ( i)**2
     1                       + Bfz ( i)**2 )

      Prs ( i) = ( gamma_eos ( i) - 1.0 ) * Prs ( i)
#if ( PRESSURE_POSITIVITY == 1)
      IF ( Prs ( i) <= prsfloor) Prs ( i) = entropy ( i)
     1                           * Rho ( i)**( gamma_eos ( i) - 1.0 )
#endif
      Prs ( i) = AMAX1 ( Prs ( i), prsfloor)

      B2 ( i) = Bfx ( i)**2 + Bfy ( i)**2 + Bfz ( i)**2

      V_dot_B ( i) = Vlx ( i) * Bfx ( i) + Vly ( i) * Bfy ( i)
     1             + Vlz ( i) * Bfz ( i)

! -----

      V ( i, 1) = Rho ( i)

      V ( i, 2) = Vlx ( i)
      V ( i, 3) = Vly ( i)
      V ( i, 4) = Vlz ( i)

      V ( i, 5) = Prs ( i)

      V ( i, 6) = Bfx ( i)
      V ( i, 7) = Bfy ( i)
      V ( i, 8) = Bfz ( i)

#if ( NFLUID >= 1)
      V ( i, 9: NFIELD + 1 + NFLUID) = temp_01 ( i)
     1                               * U ( i, 9: NFIELD + 1 + NFLUID)
#endif

! ------------------------------

      END DO

! ------------------------------------------------------------
#if ( SAVE_FLUX_MODES != 1)
! ------------------------------------------------------------

      DO j = 1, NFIELD + 1 + NFLUID

      DO i = 1, nzone_vrtx

! ------------------------------

      F ( i, j) = Vlx ( i) * U ( i, j)

! -----

      G ( i, j) = Vly ( i) * U ( i, j)

! ------------------------------

      END DO

      END DO

! ------------------------------------------------------------

      DO i = 1, nzone_vrtx

! ------------------------------

      F ( i, 2) = F ( i, 2)
     1             + Prs ( i) + ob_8pi * B2 ( i)
     1             - ob_4pi * Bfx ( i)**2

      F ( i, 3) = F ( i, 3)
     1             - ob_4pi * Bfx ( i) * Bfy ( i)

      F ( i, 4) = F ( i, 4)
     1             - ob_4pi * Bfx ( i) * Bfz ( i)

      F ( i, 5) = F ( i, 5)
     1         + ( Prs ( i) + ob_8pi * B2 ( i) ) * Vlx ( i)
     1         - Bfx ( i) * V_dot_B ( i) * ob_4pi

      F ( i, 6) = F ( i, 6) - Vlx ( i) * Bfx ( i)

      F ( i, 7) = F ( i, 7) - Vly ( i) * Bfx ( i)

      F ( i, 8) = F ( i, 8) - Vlz ( i) * Bfx ( i)

! -----

      G ( i, 2) = G ( i, 2)
     1             - ob_4pi * Bfx ( i) * Bfy ( i)

      G ( i, 3) = G ( i, 3)
     1             + Prs ( i) + ob_8pi * B2 ( i)
     1             - ob_4pi * Bfy ( i)**2

      G ( i, 4) = G ( i, 4)
     1             - ob_4pi * Bfy ( i) * Bfz ( i)

      G ( i, 5) = G ( i, 5)
     1         + ( Prs ( i) + ob_8pi * B2 ( i) ) * Vly ( i)
     1         - Bfy ( i) * V_dot_B ( i) * ob_4pi

      G ( i, 6) = G ( i, 6) - Vlx ( i) * Bfy ( i)

      G ( i, 7) = G ( i, 7) - Vly ( i) * Bfy ( i)

      G ( i, 8) = G ( i, 8) - Vlz ( i) * Bfy ( i)

! ------------------------------

      END DO

! ----------------------------------------------------------------------
#endif
! ----------------------------------------------------------------------
! STEP II) Build "S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max". This
! requires viewing each state in each of the two local frames of the
! two zone boundaries on either side of it. We also view the roe-averaged
! ( or arithmetic-averaged state) in the local frame of the associated
! zone boundary.
! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1 ! cyclic mapping of states to boundaries.

! ------------------------------
! Speeds from i^th state assigned to i^th boundary.

      Vl_Neta_use ( i) = Vlx ( i) * nx ( i) + Vly ( i) * ny ( i)
      Vl_Tau_use ( i) = - Vlx ( i) * ny ( i) + Vly ( i) * nx ( i)

      Bf_Neta_use ( i) = Bfx ( i) * nx ( i) + Bfy ( i) * ny ( i)
      Bf_Tau_use ( i) = - Bfx ( i) * ny ( i) + Bfy ( i) * nx ( i)

      Rho_use ( i) = Rho ( i)
      Prs_use ( i) = Prs ( i)

      gamma_sdsp_use ( i) = gamma_sdsp ( i)

      B2_use ( i) = B2 ( i)

! ------------------------------
! For a pure HLL/C/E RS, the i^th state only contributes to "S_Tau_Min".

      temp_01 ( i) = 1.0 / Rho_use ( i)

      temp_02 ( i) = ( gamma_sdsp_use ( i) * Prs_use ( i)
     1             + ob_4pi * B2_use ( i) ) * temp_01 ( i)

      temp_03 ( i) = ob_4pi * gamma_sdsp_use ( i) * Prs_use ( i)
     1             * Bf_Tau_use ( i)**2 * temp_01 ( i)**2

      temp_04 ( i) = AMAX1 ( temp_02 ( i)**2 - 4.0 * temp_03 ( i), 0.0)

      temp_04 ( i) = 0.5 * ( temp_02 ( i) + SQRT ( temp_04 ( i) ) )

      temp_04 ( i) = SQRT ( temp_04 ( i) )

      S_Tau_Min ( i) = Vl_Tau_use ( i) - temp_04 ( i)

! Commented just to get an exact representation of the speeds for HLL/C RS.
! Feel free to uncomment just to increase extremal speeds for the sake of
! stability.
#if ( EXTREMIZE_WAVESPEEDS >= 1)
      S_Tau_Max ( i) = Vl_Tau_use ( i) + temp_04 ( i)
#endif

! ------------------------------

      temp_03 ( i) = ob_4pi * gamma_sdsp_use ( i) * Prs_use ( i)
     1             * Bf_Neta_use ( i)**2 * temp_01 ( i)**2

      temp_04 ( i) = AMAX1 ( temp_02 ( i)**2 - 4.0 * temp_03 ( i), 0.0)

      temp_04 ( i) = 0.5 * ( temp_02 ( i) + SQRT ( temp_04 ( i) ) )

      temp_04 ( i) = SQRT ( temp_04 ( i) )

      S_Neta_Min ( i) = Vl_Neta_use ( i) - temp_04 ( i)

      S_Neta_Max ( i) = Vl_Neta_use ( i) + temp_04 ( i)

! ------------------------------
! Speeds from (i+1)^th state assigned to i^th boundary.

      Vl_Neta_use ( i) = Vlx ( ip1) * nx ( i) + Vly ( ip1) * ny ( i)
      Vl_Tau_use ( i) = - Vlx ( ip1) * ny ( i) + Vly ( ip1) * nx ( i)

      Bf_Neta_use ( i) = Bfx ( ip1) * nx ( i) + Bfy ( ip1) * ny ( i)
      Bf_Tau_use ( i) = - Bfx ( ip1) * ny ( i) + Bfy ( ip1) * nx ( i)

      Rho_use ( i) = Rho ( ip1)
      Prs_use ( i) = Prs ( ip1)

      gamma_sdsp_use ( i) = gamma_sdsp ( ip1)

      B2_use ( i) = B2 ( ip1)

! ------------------------------
! For a pure HLL/C/E RS, the (i+1)^th state only contributes to "S_Tau_Max".

      temp_01 ( i) = 1.0 / Rho_use ( i)

      temp_02 ( i) = ( gamma_sdsp_use ( i) * Prs_use ( i)
     1             + ob_4pi * B2_use ( i) ) * temp_01 ( i)

      temp_03 ( i) = ob_4pi * gamma_sdsp_use ( i) * Prs_use ( i)
     1             * Bf_Tau_use ( i)**2 * temp_01 ( i)**2

      temp_04 ( i) = AMAX1 ( temp_02 ( i)**2 - 4.0 * temp_03 ( i), 0.0)

      temp_04 ( i) = 0.5 * ( temp_02 ( i) + SQRT ( temp_04 ( i) ) )

      temp_04 ( i) = SQRT ( temp_04 ( i) )

! Commented just to get an exact representation of the speeds for HLL/C RS.
! Feel free to uncomment just to increase extremal speeds for the sake of
! stability.
#if ( EXTREMIZE_WAVESPEEDS >= 1)
      S_Tau_Min ( i) = AMIN1 ( S_Tau_Min ( i),
     1                        Vl_Tau_use ( i) - temp_04 ( i))
#endif

! Commented just to get an exact representation of the speeds for HLL/C RS.
! Feel free to uncomment just to increase extremal speeds for the sake of
! stability.
#if ( EXTREMIZE_WAVESPEEDS >= 1)
      S_Tau_Max ( i) = AMAX1 ( S_Tau_Max ( i),
     1                        Vl_Tau_use ( i) + temp_04 ( i) )
#else
      S_Tau_Max ( i) = Vl_Tau_use ( i) + temp_04 ( i)
#endif

! ------------------------------

      temp_03 ( i) = ob_4pi * gamma_sdsp_use ( i) * Prs_use ( i)
     1             * Bf_Neta_use ( i)**2 * temp_01 ( i)**2

      temp_04 ( i) = AMAX1 ( temp_02 ( i)**2 - 4.0 * temp_03 ( i), 0.0)

      temp_04 ( i) = 0.5 * ( temp_02 ( i) + SQRT ( temp_04 ( i) ) )

      temp_04 ( i) = SQRT ( temp_04 ( i) )

      S_Neta_Min ( i) = AMIN1 ( S_Neta_Min ( i),
     1                         Vl_Neta_use ( i) - temp_04 ( i))

      S_Neta_Max ( i) = AMAX1 ( S_Neta_Max ( i),
     1                         Vl_Neta_use ( i) + temp_04 ( i))

! ------------------------------
! Speeds from averaging i^th and (i+1)^th states assigned to i^th boundary.
! Roe-averages built here. But default to arithmetic averages if the former
! are unphysical.

! -----
! Build "Rho_use" and "gamma_*_use".

      temp_01 ( i) = SQRT ( Rho ( ip1) / Rho ( i) )
      temp_02 ( i) = 1.0 / ( 1.0 + temp_01 ( i) )

      gamma_sdsp_use ( i) =0.5 * ( gamma_sdsp ( i) + gamma_sdsp ( ip1))
      gamma_eos_use ( i) =0.5 * ( gamma_eos ( i) + gamma_eos ( ip1))

      Rho_use ( i) = SQRT ( Rho ( i) * Rho ( ip1) )

! -----
! Build "Vl_Neta_use, Vl_Tau_use, Vlz_use".

      Vl_Neta_use ( i) = ( Vlx ( i) * nx ( i) + Vly ( i)  * ny ( i) )
     1 + ( Vlx ( ip1) * nx ( i) + Vly ( ip1)  * ny ( i) ) * temp_01 ( i)
      Vl_Neta_use ( i) = Vl_Neta_use ( i) * temp_02 ( i)

      Vl_Tau_use ( i) = ( -Vlx ( i) * ny ( i) + Vly ( i) * nx ( i) )
     1 + ( -Vlx ( ip1) * ny ( i) + Vly ( ip1) * nx ( i) ) * temp_01 ( i)
      Vl_Tau_use ( i) = Vl_Tau_use ( i) * temp_02 ( i)

      Vlz_use ( i) = ( Vlz ( i) + Vlz ( ip1) * temp_01 ( i) )
     1             * temp_02 ( i)

! -----
! Build "Bf_Neta_use, Bf_Tau_use, Bfz_use, X_use, B2_use".

      Bf_Neta_use ( i) =
     1   ( Bfx ( i) * nx ( i) + Bfy ( i) * ny ( i) ) * temp_01 ( i)
     1 + ( Bfx ( ip1) * nx ( i) + Bfy ( ip1) * ny ( i) )
      Bf_Neta_use ( i) = Bf_Neta_use ( i) * temp_02 ( i)

      Bf_Tau_use ( i) =
     1   ( -Bfx ( i) * ny ( i) + Bfy ( i) * nx ( i) ) * temp_01 ( i)
     1 + ( -Bfx ( ip1) * ny ( i) + Bfy ( ip1) * nx ( i) )
      Bf_Tau_use ( i) = Bf_Tau_use ( i) * temp_02 ( i)

      Bfz_use ( i) = ( Bfz ( i) * temp_01 ( i) + Bfz ( ip1) )
     1             * temp_02 ( i)

      X_use ( i) = ( ( Bfx ( ip1) * nx ( i) + Bfy ( ip1) * ny ( i) )
     1            - ( Bfx ( i) * nx ( i) + Bfy ( i) * ny ( i) ) )**2
     1          + ( ( -Bfx ( ip1) * ny ( i) + Bfy ( ip1) * nx ( i) )
     1            - ( -Bfx ( i) * ny ( i) + Bfy ( i) * nx ( i) ) )**2
     1          + ( Bfz ( ip1) - Bfz ( i) )**2
      X_use ( i) = 0.5 * X_use ( i)
     1          / ( SQRT ( Rho ( i)) + SQRT ( Rho ( ip1)) )**2

      B2_use ( i) = Bf_Neta_use ( i)**2 + Bf_Tau_use ( i)**2
     1            + Bfz_use ( i)**2

! -----
! Build "Total_Enthalpy_*" and therefrom "Prs_use".

      Total_Enthalpy_remi ( i) =
     1     ( gamma_eos ( i) * Prs ( i) / ( gamma_eos ( i) - 1.0)
     1     + 0.5 * Rho ( i)
     1     * ( Vlx ( i)**2 + Vly ( i)**2 + Vlz ( i)**2)
     1     + ob_4pi * ( Bfx ( i)**2 + Bfy ( i)**2 + Bfz ( i)**2)
     1                           ) / Rho ( i)

      Total_Enthalpy_repl ( i) =
     1     ( gamma_eos ( ip1) * Prs ( ip1) / ( gamma_eos ( ip1) - 1.0)
     1     + 0.5 * Rho ( ip1) 
     1     * ( Vlx ( ip1)**2 + Vly ( ip1)**2 + Vlz ( ip1)**2)
     1     + ob_4pi * ( Bfx ( ip1)**2 + Bfy ( ip1)**2 + Bfz ( ip1)**2)
     1                            ) / Rho ( ip1)

      Total_Enthalpy_use ( i) = ( Total_Enthalpy_remi ( i)
     1     + Total_Enthalpy_repl ( i) * temp_01 ( i) ) * temp_02 ( i)

      Prs_use ( i) = Total_Enthalpy_use ( i)
     1     - 0.5 * ( Vl_Neta_use ( i)**2 + Vl_Tau_use ( i)**2
     1             + Vlz_use ( i)**2 )
     1     - ob_4pi * B2_use ( i) / Rho_use ( i)

      Prs_use ( i) = ( gamma_eos_use ( i) - 1.0) * Rho_use ( i)
     1             * Prs_use ( i) / gamma_eos_use ( i)

      Prs_use ( i) = Prs_use ( i)
     1   - ( gamma_eos_use ( i) - 2.0) * Rho_use ( i) * X_use ( i)
     1   / ( fourpi * gamma_eos_use ( i) )

      Prs_use ( i) = AMAX1 ( Prs_use ( i), prsfloor)

! ----------

      IF ( ( Rho_use ( i) <= rhofloor) .OR.
     1     ( Prs_use ( i) <= prsfloor) ) THEN

      Vl_Neta_use ( i) = 0.5 * ( ( Vlx ( i) + Vlx ( ip1)) * nx ( i)
     1                         + ( Vly ( i) + Vly ( ip1)) * ny ( i) )
      Vl_Tau_use ( i) = 0.5 * ( - ( Vlx ( i) + Vlx ( ip1)) * ny ( i)
     1                          + ( Vly ( i) + Vly ( ip1)) * nx ( i) )

      Bf_Neta_use ( i) = 0.5 * ( ( Bfx ( i) + Bfx ( ip1)) * nx ( i)
     1                         + ( Bfy ( i) + Bfy ( ip1)) * ny ( i) )
      Bf_Tau_use ( i) = 0.5 * ( - ( Bfx ( i) + Bfx ( ip1)) * ny ( i)
     1                          + ( Bfy ( i) + Bfy ( ip1)) * nx ( i) )

      Rho_use ( i) = 0.5 * ( Rho ( i) + Rho ( ip1))
      Prs_use ( i) = 0.5 * ( Prs ( i) + Prs ( ip1))

      gamma_sdsp_use ( i) =0.5 * ( gamma_sdsp ( i) + gamma_sdsp ( ip1))

      B2_use ( i) = Bf_Neta_use ( i)**2 + Bf_Tau_use ( i)**2
     1            + 0.25 * ( Bfz ( i) + Bfz ( ip1))**2

      X_use ( i) = 0.0

      END IF

! ------------------------------
! For a pure HLL/C/E RS, the Roe-averaged state contributes to
! "S_Tau_Min" and "S_Tau_Max".

      temp_01 ( i) = 1.0 / Rho_use ( i)

      temp_02 ( i) = ( gamma_sdsp_use ( i) * Prs_use ( i)
     1             + ob_4pi * B2_use ( i) ) * temp_01 ( i)

      temp_03 ( i) = ob_4pi * gamma_sdsp_use ( i) * Prs_use ( i)
     1             * Bf_Tau_use ( i)**2 * temp_01 ( i)**2

      temp_04 ( i) = AMAX1 ( temp_02 ( i)**2 - 4.0 * temp_03 ( i), 0.0)

      temp_04 ( i) = 0.5 * ( temp_02 ( i) + SQRT ( temp_04 ( i) ) )

      temp_04 ( i) = SQRT ( temp_04 ( i) )

      S_Tau_Min ( i) = AMIN1 ( S_Tau_Min ( i),
     1                        Vl_Tau_use ( i) - temp_04 ( i))

      S_Tau_Max ( i) = AMAX1 ( S_Tau_Max ( i),
     1                        Vl_Tau_use ( i) + temp_04 ( i))

! ------------------------------

      temp_03 ( i) = ob_4pi * gamma_sdsp_use ( i) * Prs_use ( i)
     1             * Bf_Neta_use ( i)**2 * temp_01 ( i)**2

      temp_04 ( i) = AMAX1 ( temp_02 ( i)**2 - 4.0 * temp_03 ( i), 0.0)

      temp_04 ( i) = 0.5 * ( temp_02 ( i) + SQRT ( temp_04 ( i) ) )

      temp_04 ( i) = SQRT ( temp_04 ( i) )

      S_Neta_Min ( i) = AMIN1 ( S_Neta_Min ( i),
     1                         Vl_Neta_use ( i) - temp_04 ( i))

      S_Neta_Max ( i) = AMAX1 ( S_Neta_Max ( i),
     1                         Vl_Neta_use ( i) + temp_04 ( i))

! ------------------------------

      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE ONE_D_MHD_PRIM_F_G_WAVSPD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ONE_D_MHD_RS_HLL_TAU ( nzone_vrtx,

     1     smallnum, prsfloor, rhofloor, pi, blend_HLL,
     1     expand_wavespeeds,

     1     U, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     V, F, G,

     1     S_Tau_Min, S_Tau_Mid, S_Tau_Max, S_Neta_Min, S_Neta_Max,

     1     U_Minus, U_Star_Minus, U_Star_Plus, U_Plus,
     1     F_Minus, F_Star_Minus, F_Star_Plus, F_Plus,
     1     G_Minus, G_Star_Minus, G_Star_Plus, G_Plus,

     1     U_Star_Avg, F_Star_Avg, G_Star_Avg,

     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! Builds a 1d Riemann solver in the Tau-direction, returning the very
! variables that are needed in a 2d RS.

! Here the HLLC RS is forced to default to an HLL RS. Thus the intermediate
! steps in formulating an HLLC RS are avoided.

! As input:-

! We come in with conserved variables "U" arranged cyclically around a
! vertex in an anti-clockwise direction. We also come in with polytropic
! indices. Also come in with some sense of how much HLL we will have to
! blend in for stabilization.

! We also come in with vectors of primitive variables "V" and
! x- and y-directional fluxes, "F, G". These again have a cyclic and
! anti-clockwise arrangement around the vertex in question.

! The normal vector components "nx, ny" are also provided.
! Here "Neta = ( nx, ny)" is the unit normal along faces between triangles; 
! "Tau = ( -ny, nx)" is orthogonal to "Neta" and arranged to have an
! anti-clockwise direction since this is also the ordering of the incoming
! states.

! The above variables are then used to find the left, middle and right
! speeds of the HLLC RS in Tau-direction, "S_Tau_Min, S_Tau_Mid, S_Tau_Max".
! "S_Tau_Min, S_Tau_Max" are input but can be reset here. "S_Tau_Mid"
! is output.

! The extremal speeds in the Neta-direction are "S_Neta_Min, S_Neta_Max".
! They can also be reset here if needed.

! As output:-

! The HLL state and flux, "U_HLL and F_HLL" are built. Using them,
! we also build the total gas + magnetic pressure in the resolved
! state "Prs_Tot_Star". The HLL RS also gives other "Rho_Star" etc. states.
! Other "F_Add and G_Add" terms are built, which simplify the process of
! obtaining the HLLC RS.

! The HLLC formulation then yields the four conserved vectors (two bounding
! & two intermediate states), "U_Minus, U_Star_Minus, U_Star_Plus, U_Plus".

! These can be used to build the Tau-directional fluxes (two bounding &
! two intermediate fluxes), "F_Minus, F_Star_Minus, F_Star_Plus, F_Plus". 

! We also build the Neta-directional fluxes (two bounding &
! two intermediate fluxes), "G_Minus, G_Star_Minus, G_Star_Plus, G_Plus".

! For calculational convenience, we also return
! "U_Star_Avg, G_Star_Avg", which are the averages of the two
! resolved states and two resolved neta-fluxes.

! Lastly, the state and Tau-flux from the HLLC RS is also returned in
! "U_HLLC, F_HLLC, G_HLLC". These are the resolved numerical fluxes in
! the Tau- and Neta-directions.

! "F_HLLC_Global_x, G_HLLC_Global_y" are just rotated versions of the same
! in the global coordinate system. They are especially useful in building
! the supersonic fluxes.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx

      REAL smallnum, prsfloor, rhofloor, pi, expand_wavespeeds,

! Amount of HLL we will have to blend in for stabilization.
     1     blend_HLL


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables that come together at a zone edge.
! We also come in with their primitive variables and associated
! x- and y-fluxes.
     1      U, V, F, G,

! The two resolved states and their fluxes in the x- and y-directions.
     1      U_Minus, U_Star_Minus, U_Star_Plus, U_Plus,
     1      F_Minus, F_Star_Minus, F_Star_Plus, F_Plus,
     1      G_Minus, G_Star_Minus, G_Star_Plus, G_Plus,

! Averages of the two resolved states and two resolved y-fluxes.
     1      U_Star_Avg, F_Star_Avg, G_Star_Avg,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

! Components of the unit normal.
     1      nx, ny,

! The left, middle and right speeds of the HLLC RS in the Tau-direction,
! "S_Tau_Min, S_Tau_Mid, S_Tau_Max".
     1      S_Tau_Min, S_Tau_Mid, S_Tau_Max,

! The extremal speeds in the Neta-direction, "S_Neta_Min, S_Neta_Max".
     1      S_Neta_Min, S_Neta_Max,

     1      entropy, F_entropy_flux_x, G_entropy_flux_y

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      INTEGER i, j, ip1, iredo_HLL

      REAL fourpi, ob_4pi, ob_8pi, use_blend_HLL, tempaa


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Primitives obtained from the resolved HLL state. Useful for building
! "F_Add and G_Add".
! The total gas + magnetic pressure in the resolved state "Prs_Tot_Star".
     1      Rho_Star, Prs_Tot_Star, Vlx_Star, Vly_Star, Vlz_Star,
     1      Bfx_Star, Bfy_Star, Bfz_Star, V_dot_B_Star,
     1      B2_Star, Prs_Star, gamma_sdsp_Star,

! HLLC really requires deft handling of the pressures to work well.
! We build the supplementary variables here.
     &      Rho_Minus, Rho_Plus, Prs_Tot_Minus, Prs_Tot_Plus,           
     &      Vtau_Minus, Vtau_Plus, Btau_Minus, Btau_Plus,               

! Velocity and magnetic field components in the "Neta, Tau" directions.
     1      Vl_Tau, Bf_Tau, Vl_Neta, Bf_Neta,

! Polytropic indices on either side.
     1      gamma_eos_Lf, gamma_eos_Rt,
     1      gamma_sdsp_Lf, gamma_sdsp_Rt,

! Redone versions of wave speeds in case we want to expand the wave model.
     1      S_Tau_Min_Prev, S_Tau_Max_Prev,
     1      S_Neta_Min_Prev, S_Neta_Max_Prev,

     1      entropy_flux_Minus, entropy_flux_Plus, entropy_HLL,

     1      temp_01, temp_02, temp_03, temp_04, temp_05, temp_06


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Vectors that carry the HLL solution, needed for building the "_Star" states.
     1      U_HLL, F_HLL, G_HLL,

! Vectors that are added to speed * conserved vector to get fluxes.
     1      F_Add, G_Add

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 / fourpi

      IF ( MAX_ZONES_PER_VRTX < nzone_vrtx) THEN
      WRITE ( 6, *)"STOPPING in ONE_D_MHD_RS_HLL_TAU because : "
      WRITE ( 6, *)"NZONES_","PER_VRTX < nzone_vrtx",
     1              MAX_ZONES_PER_VRTX, nzone_vrtx
      STOP
      END IF

      IF ( NFIELD /= 7) THEN
      WRITE ( 6, *)"STOPPING in ONE_D_MHD_RS_HLL_TAU because : "
      WRITE ( 6, *)"NFIELD /= 7", NFIELD
      STOP
      END IF

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Since the resolved flux is being evaluated in two stages, reinitialize the
! conserved variables. Also do the same for the fluxes.

! Note though that the fluxes are in the facially rotated frame. Thus
! "F_Minus and F_Plus" are along the "Tau" direction and
! "G_Minus and G_Plus" are in the "Neta" direction.

! Here "Neta = ( nx, ny)" is the unit normal along faces between triangles;
! "Tau = ( -ny, nx)" is orthogonal to "Neta" and arranged to have an
! anti-clockwise direction since this is also the ordering of the incoming
! states.

! Care should be taken to ensure that all the variables needed later on
! in this subroutine are reinitialized here.

! See eqn. (B.2) and eqns. after it in BDA14 to see how
! "F_Minus, F_Plus" and "G_Minus, G_Plus" are defined.

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1  ! cyclic mapping of states to boundaries.

! -----

      U_Minus ( i, :) = U ( i, :)
      F_Minus ( i, :) = - ny ( i) * F ( i, :) + nx ( i) * G ( i, :)
      G_Minus ( i, :) = nx ( i) * F ( i, :) + ny ( i) * G ( i, :)

      Rho_Minus ( i) = V ( i, 1)
      Prs_Tot_Minus ( i) = V ( i, 5) + ob_8pi * ( V ( i, 6)**2          
     &                           + V ( i, 7)**2 + V ( i, 8)**2 )
      Vtau_Minus ( i) = - ny ( i) * V ( i, 2) + nx ( i) * V ( i, 3)
      Btau_Minus ( i) = - ny ( i) * V ( i, 6) + nx ( i) * V ( i, 7)

      gamma_sdsp_Lf ( i) = gamma_sdsp ( i)
      gamma_eos_Lf ( i) = gamma_eos ( i)

! -----

      U_Plus ( i, :) = U ( ip1, :)
      F_Plus ( i, :) = - ny ( i) * F ( ip1, :) + nx ( i) * G ( ip1, :)
      G_Plus ( i, :) = nx ( i) * F ( ip1, :) + ny ( i) * G ( ip1, :)

      Rho_Plus ( i) = V ( ip1, 1)
      Prs_Tot_Plus ( i) = V ( ip1, 5) + ob_8pi * ( V ( ip1, 6)**2       
     &                            + V ( ip1, 7)**2 + V ( ip1, 8)**2 )
      Vtau_Plus ( i) = - ny ( i) * V ( ip1, 2) + nx ( i) * V ( ip1, 3)
      Btau_Plus ( i) = - ny ( i) * V ( ip1, 6) + nx ( i) * V ( ip1, 7)

      gamma_sdsp_Rt ( i) = gamma_sdsp ( ip1)
      gamma_eos_Rt ( i) = gamma_eos ( ip1)

! -----

      END DO

! ----------------------------------------------------------------------
! STEP_1d_HLLC_TAU I) First build the HLL state and flux "U_HLL and F_HLL".
! Use it to build "S_Tau_Mid, Prs_Tot_Star, Rho_Star, Vlx_Star, Vly_Star,
! Vlz_Star, Bfx_Star, Bfy_Star, Bfz_Star" from the intermediate state
! of the HLL RS.
! ----------------------------------------------------------------------

! We also build "Vl_Tau, Bf_Tau, Vl_Neta, Bf_Neta" from the above states.

! Use the above to build "F_Add, G_Add", which facilitate HLLC construction.
! Then use the above to build "G_HLL".

! This is eqns. (A.8) and (A.9) of BDA14.

      DO i = 1, nzone_vrtx

! ------------------------------

      temp_01 ( i) = 1.0 / ( S_Tau_Max ( i) - S_Tau_Min ( i) )
      temp_02 ( i) = S_Tau_Max ( i) * temp_01 ( i)
      temp_03 ( i) = S_Tau_Min ( i) * temp_01 ( i)

! ------------------------------

      END DO

! ------------------------------------------------------------

      DO j = 1, NFIELD + 1 + NFLUID

      DO i = 1, nzone_vrtx

! ------------------------------

      U_HLL ( i, j) = temp_02 ( i) * U_Plus ( i, j)
     1              - temp_03 ( i) * U_Minus ( i, j)
     1              - temp_01 ( i)
     1              * ( F_Plus ( i, j) - F_Minus ( i, j) )

      F_HLL ( i, j) = temp_02 ( i) * F_Minus ( i, j)
     1              - temp_03 ( i) * F_Plus ( i, j)
     1             + S_Tau_Max ( i) * S_Tau_Min ( i) * temp_01 ( i)
     1              * ( U_Plus ( i, j) - U_Minus ( i, j) )

! ------------------------------

      END DO

      END DO

! ------------------------------------------------------------

      DO i = 1, nzone_vrtx

! ------------------------------
! Find Intermediate states for HLL RS. This is eqn. (B.4) of BDA14.
! Recall that these states are in terms of x,y,z-coordinates so that
! we get x,y,z-velocities and magnetic fields.

! In the off-chance that this is a situation where the longitudinal field
! is discontinuous, the RS includes compensating terms related to the
! jump in the longitudinal field.

      Rho_Star ( i) = AMAX1 ( U_HLL ( i, 1), rhofloor)
      temp_04 ( i) = 1.0 / Rho_Star ( i)

      Vlx_Star ( i) = U_HLL ( i, 2) * temp_04 ( i)
      Vly_Star ( i) = U_HLL ( i, 3) * temp_04 ( i)
      Vlz_Star ( i) = U_HLL ( i, 4) * temp_04 ( i)

      Bfx_Star ( i) = U_HLL ( i, 6)
      Bfy_Star ( i) = U_HLL ( i, 7)
      Bfz_Star ( i) = U_HLL ( i, 8)

      B2_Star ( i) = Bfx_Star ( i)**2 + Bfy_Star ( i)**2
     1             + Bfz_Star ( i)**2

      V_dot_B_Star ( i) = Vlx_Star ( i) * Bfx_Star ( i)
     1                  + Vly_Star ( i) * Bfy_Star ( i)
     1                  + Vlz_Star ( i) * Bfz_Star ( i)

! Velocity and magnetic field in the rotated frame.
      Vl_Tau ( i) = - ny ( i) * Vlx_Star ( i) + nx ( i) * Vly_Star ( i)
      Bf_Tau ( i) = - ny ( i) * Bfx_Star ( i) + nx ( i) * Bfy_Star ( i)

      Vl_Neta ( i) = nx ( i) * Vlx_Star ( i) + ny ( i) * Vly_Star ( i)
      Bf_Neta ( i) = nx ( i) * Bfx_Star ( i) + ny ( i) * Bfy_Star ( i)

! This important speed mediates the entire HLLC Riemann solver.

       S_Tau_Mid ( i) = ( Prs_Tot_Plus ( i) - Prs_Tot_Minus ( i)
     1   + Rho_Minus ( i) * Vtau_Minus ( i)
     1                    * ( S_Tau_Min ( i) - Vtau_Minus ( i) )
     1   - Rho_Plus ( i) * Vtau_Plus ( i)
     1                   * ( S_Tau_Max ( i) - Vtau_Plus ( i) )
     1   - ob_4pi * ( Btau_Plus ( i)**2 - Btau_Minus ( i)**2 ) )
     1 / ( Rho_Minus ( i) * ( S_Tau_Min ( i) - Vtau_Minus ( i) )
     1   - Rho_Plus ( i) * ( S_Tau_Max ( i) - Vtau_Plus ( i) ) )

! Restore consistency in the longitudinal velocity.
        Vl_Tau ( i) = S_Tau_Mid ( i)

! Total pressure from HLL resolved state. Note, fluxes are rotated into the
! "tau" direction.

      Prs_Tot_Star ( i) = 0.5 * ( Prs_Tot_Minus ( i) + Prs_Tot_Plus ( i)
     1   + Rho_Minus ( i) * ( S_Tau_Min ( i) - Vtau_Minus ( i) )
     1                    * ( S_Tau_Mid ( i) - Vtau_Minus ( i) )
     1   + Rho_Plus ( i) * ( S_Tau_Max ( i) - Vtau_Plus ( i) )
     1                   * ( S_Tau_Mid ( i) - Vtau_Plus ( i) )
     1   + ob_4pi * ( 2.0 * Bf_Tau ( i)**2
     1              - Btau_Plus ( i)**2 - Btau_Minus ( i)**2 ) )

      Prs_Tot_Star ( i) = AMAX1 ( Prs_Tot_Star ( i), prsfloor)

      Prs_Star ( i) = Prs_Tot_Star ( i) - ob_8pi * B2_Star ( i)

      Prs_Star ( i) = AMAX1 ( Prs_Star ( i), prsfloor)

! ------------------------------

      END DO

! ----------------------------------------------------------------------
! Build the longitudinal and transverse entropy flux. It is always an HLL flux. 

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1  ! cyclic mapping of states to boundaries.

! ------------------------------

      temp_01 ( i) = AMIN1 ( S_Tau_Min ( i), - smallnum)  ! Minimum speed
      temp_02 ( i) = AMAX1 ( S_Tau_Max ( i), smallnum)  ! Maximum speed

      temp_03 ( i) = 1.0 / ( temp_02 ( i) - temp_01 ( i) )
      temp_04 ( i) = temp_02 ( i) * temp_03 ( i)
      temp_05 ( i) = temp_01 ( i) * temp_03 ( i)
      temp_06 ( i) = temp_02 ( i) * temp_01 ( i) * temp_03 ( i)

      entropy_flux_Minus ( i) = entropy ( i) * Vtau_Minus ( i)
      entropy_flux_Plus ( i) = entropy ( ip1) * Vtau_Plus ( i)

      entropy_HLL ( i) = temp_04 ( i) * entropy ( ip1)
     1              - temp_05 ( i) * entropy ( i)
     1              - temp_03 ( i) * ( entropy_flux_Plus ( i)
     1                               - entropy_flux_Minus ( i) )

      F_entropy_flux_x ( i) = temp_04 ( i) * entropy_flux_Minus ( i)
     1              - temp_05 ( i) * entropy_flux_Plus ( i)
     1              + temp_06 ( i) * ( entropy ( ip1) - entropy ( i) )

      G_entropy_flux_y ( i) = entropy_HLL ( i) * Vl_Tau ( i)

! ------------------------------

      END DO

! ----------------------------------------------------------------------
! STEP_1d_HLLC_TAU III) Evaluate "F_Add" and "G_Add", then build "G_HLL".
! This completes all the steps for the HLL RS.
! ----------------------------------------------------------------------

! Only the first five components of "F_Add" and "G_Add" matter.
! This is eqn. (B.5) of BDA14.

      DO i = 1, nzone_vrtx

! ------------------------------

      F_Add ( i, 1) = 0.0

      F_Add ( i, 2) = - ny ( i) * Prs_Tot_Star ( i)
     1              - ob_4pi * Bfx_Star ( i) * Bf_Tau ( i)

      F_Add ( i, 3) = nx ( i) * Prs_Tot_Star ( i)
     1              - ob_4pi * Bfy_Star ( i) * Bf_Tau ( i)

      F_Add ( i, 4) = - ob_4pi * Bfz_Star ( i) * Bf_Tau ( i)

      F_Add ( i, 5) = Prs_Tot_Star ( i) * Vl_Tau ( i)
     1              - ob_4pi * Bf_Tau ( i) * V_dot_B_Star ( i)

      F_Add ( i, 6) = 0.0

      F_Add ( i, 7) = 0.0

      F_Add ( i, 8) = 0.0

#if ( NFLUID >= 1)
      F_Add ( i, 9: NFIELD + 1 + NFLUID) = 0.0
#endif

! -----

      G_Add ( i, 1) = 0.0

      G_Add ( i, 2) = nx ( i) * Prs_Tot_Star ( i) 
     1              - ob_4pi * Bfx_Star ( i) * Bf_Neta ( i)

      G_Add ( i, 3) = ny ( i) * Prs_Tot_Star ( i)
     1              - ob_4pi * Bfy_Star ( i) * Bf_Neta ( i)

      G_Add ( i, 4) = - ob_4pi * Bfz_Star ( i) * Bf_Neta ( i)

      G_Add ( i, 5) = Prs_Tot_Star ( i) * Vl_Neta ( i)
     1              - ob_4pi * Bf_Neta ( i) * V_dot_B_Star ( i)

      G_Add ( i, 6) = 0.0

      G_Add ( i, 7) = 0.0

      G_Add ( i, 8) = 0.0

#if ( NFLUID >= 1)
      G_Add ( i, 9: NFIELD + 1 + NFLUID) = 0.0
#endif

! ------------------------------

      END DO

! ------------------------------------------------------------

! This evaluation will be wrong for components 6, 7 and 8.
! This is eqn. (B.5) of BDA14.

      DO j = 1, NFIELD + 1 + NFLUID

      DO i = 1, nzone_vrtx

! -----

      G_HLL ( i, j) = Vl_Neta ( i) * U_HLL ( i, j) + G_Add ( i, j)

! -----

      END DO

      END DO

! ------------------------------

! We fix components 6, 7 and 8 here.

      DO i = 1, nzone_vrtx

! -----

      temp_04 ( i) = Vlx_Star ( i) * Bfy_Star ( i)
     1             - Vly_Star ( i) * Bfx_Star ( i)

      G_HLL ( i, 6) = - ny ( i) * temp_04 ( i)
      G_HLL ( i, 7) = nx ( i) * temp_04 ( i)
      G_HLL ( i, 8) = Bfz_Star ( i) * Vl_Neta ( i)
     1              - Vlz_Star ( i) * Bf_Neta ( i)

! -----

      END DO

! ----------------------------------------------------------------------
! Here we default the 1d HLLC into 1d HLL.

      DO i = 1, nzone_vrtx

      U_Star_Minus ( i, :) = U_HLL ( i, :)
      F_Star_Minus ( i, :) = F_HLL ( i, :)
      G_Star_Minus ( i, :) = G_HLL ( i, :)

      U_Star_Plus ( i, :) = U_HLL ( i, :)
      F_Star_Plus ( i, :) = F_HLL ( i, :)
      G_Star_Plus ( i, :) = G_HLL ( i, :)

      END DO

! ----------------------------------------------------------------------
! "U_Star_Avg, G_Star_Avg" just carry the information that will be used
! to build the strongly interacting state and fluxes in the
! multidimensional HLL RS.

      DO j = 1, NFIELD + 1 + NFLUID

      DO i = 1, nzone_vrtx

! ------------------------------

      U_Star_Avg ( i, j) = U_HLL ( i, j)

      F_Star_Avg ( i, j) = F_HLL ( i, j)

      G_Star_Avg ( i, j) = G_HLL ( i, j)

! ------------------------------

      END DO

      END DO

! ----------------------------------------------------------------------
! STEP_1d_HLLC_TAU IV) Take care of exception conditions. If they are
! encountered, we return an HLL state instead of an HLLC state.
! ----------------------------------------------------------------------

      use_blend_HLL = 1.0

        S_Tau_Mid ( i) = AMAX1 ( AMIN1 ( S_Tau_Mid ( i),
     1                           S_Tau_Max ( i) ), S_Tau_Min ( i) )

! ----------------------------------------------------------------------
! STEP_1d_HLLC_TAU V) Obtain the Tau-direction flux from the HLLC RS
! and return it in "F_HLLC". Similarly, return the Neta-direction flux
! in "G_HLLC".
! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

! ------------------------------

      IF ( S_Tau_Min ( i) >= 0.0 ) THEN

        U_HLLC ( i, :) = U_Minus ( i, :)
        F_HLLC ( i, :) = F_Minus ( i, :)
        G_HLLC ( i, :) = G_Minus ( i, :)

      ELSE IF ( S_Tau_Max ( i) <= 0.0 ) THEN

        U_HLLC ( i, :) = U_Plus ( i, :)
        F_HLLC ( i, :) = F_Plus ( i, :)
        G_HLLC ( i, :) = G_Plus ( i, :)

      ELSE
! This is defaulted to an HLL RS because we set 
! "U_Star_Minus = U_Star_Plus = U_HLL" and so on in a previous step.

        U_HLLC ( i, :) = U_Star_Minus ( i, :)
        F_HLLC ( i, :) = F_Star_Minus ( i, :)
        G_HLLC ( i, :) = G_Star_Minus ( i, :)

      END IF

! ----------

        F_HLLC_Global_x ( i, :) = G_HLLC ( i, :) * nx ( i)
     1                          - F_HLLC ( i, :) * ny ( i) 
        G_HLLC_Global_y ( i, :) = G_HLLC ( i, :) * ny ( i)
     1                          + F_HLLC ( i, :) * nx ( i)

! ------------------------------

      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE ONE_D_MHD_RS_HLL_TAU

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE MULTID_WAVE_MODEL                                      

     &   ( nzone_vrtx, multid_wavemodel_num,                            

     &     smallnum, blend_HLL, expand_wavespeeds, Vlx_ALE, Vly_ALE,    

     &     U, gamma_eos, gamma_sdsp, S_1DRS,                            

     &     S_Neta_Min, S_Neta_Max, S_Tau_Min, S_Tau_Max, nx, ny,        

     &     Vlx_Centroid, Vly_Centroid, S_Radial_Max,                    
     &     Is_SuperSonic, Is_SuperSonic_State,                          

     &     S_Neta_Multid, L_Tau_Min, L_Tau_Max, x, y,                   
     &     neta_x, neta_y, neta_dot_n, neta_dot_t, L_1DRS,              

     &     IF_Mid, LF_Tau_Mid, yF_Lo, yF_Hi,                            
     &     IG_Mid, LG_Tau_Mid, xG_Lo, xG_Hi,                            
     &     Int_Momnt_U, Int_Momnt_F, Int_Momnt_G, ifail)

! ----------------------------------------------------------------------

! Written by : Dinshaw S. Balsara

! On input:-

! "multid_wavemodel_num" choice of multid wave model. 0 gives original
! wave model; 1 gives Variant # 1; 2 gives Variant # 2; 3 gives Variant # 3.

! Takes in "S_Neta_Min, S_Neta_Max, S_Tau_Min, S_Tau_Max" wave speeds
! and unit normal "nx, ny" to yield "S_Neta_Multid", which is the speed
! of fastest wave propagation along "(nx, ny)". It will also build a lot
! of other geometric information associated with the wave speeds.

! "S_1DRS ( 1: nzone_vrtx, NWAVES_1DRS)" comes in with the wave speeds
! of all the 1d RP speeds. It is left untouched.

! "blend_HLL" tells us how much HLL we need to blend in. If it is large and
! if the problem is subsonic (and extreme), one may use "EXTREMIZE_WAVESPEEDS"
! to reset "S_Tau_Min, S_Tau_Max" to increase dissipation. 

! "expand_wavespeeds" Factor by which subsonic wavespeeds are expanded
! just to be on the safe side.

! "Vlx_ALE, Vly_ALE" Velocity with which the ALE mesh moves. Usually, the
! fluxes and state relative to that direction is desired. It is used to
! discern whether the wave model is supersonic relative to this velocity.

! ------------------------------

! On output:-

! "S_Neta_Multid" is the largest speed in the "Neta" direction.
! For supersonic flows it can be negative in certain directions.

! "L_Tau_Min, L_Tau_Max" are the minimum and maximum speeds within each
! face. They are formed by the intersection of the sides of the polygon
! along its principal directions. They are measured relative to the
! starting vertex of each face.

! "(x(i), y(i))" is the coordinate point in global space for the
! starting point of the each i^th face. "(x(i+1), y(i+1))" is the
! ending point of the same face.

! "Is_SuperSonic" is set to 0 if wave model is subsonic. Otherwise, it is
! set to the side number of the wave model that contains the supersonic
! flow. This facilitates the search for the fluxes in the supersonic case.

! "Is_SuperSonic_State" gives the state number [0: NWAVES_1DRS] in the
! face picked out by "Is_SuperSonic". If the wave model is supersonic
! relative to "Vlx_ALE, Vly_ALE" then "Is_SuperSonic, Is_SuperSonic_State"
! fully specify the information that is to be returned by the multid RS.

! Our best estimate of the velocity of the centroid of the wave model is 
! given as "Vlx_Centroid, Vly_Centroid".

! The largest, positive speed in this Riemann problem is "S_Radial_Max".

! "neta_x, neta_y" is the unit vector joining i^th and (i+1)^th vertex of
! the wave model.

! "neta_dot_n, neta_dot_t" are dot products between unit vector "neta"
! and unit vectors "n" and "t". They are very useful for flux evaluation
! on the boundary of the wave model.

! "L_1DRS ( 1: nzone_vrtx, NWAVES_1DRS)" parametrizes "S_1DRS" within
! each face so that it goes from "L_Tau_Min" to "L_Tau_Max".

! "IF_Mid ( i)" is 1 if the entire i^th face has x <= 0 and will contribute
! to the evaluation of the x-flux. It is 0 if the entire i^th face has x > 0
! and will not contribute at all to the evaluation of the x-flux.
! If x ( ip1) <= 0 and x ( i) > 0 then only [ LF_Tau_Mid ( i), L_Tau_Max ( i)]
! will contribute to the x-flux. In that case "IF_Mid ( i) = - 1".
! If x ( ip1) > 0 and x ( i) <= 0 then only [ L_Tau_Min ( i), LF_Tau_Mid ( i)]
! will contribute to the x-flux. In that case "IF_Mid ( i) = - 2".
! If "IF_Mid ( i)" is negative for the i^th face, then "yF_Lo, yF_Hi" will
! contain the y-coordinates of the two faces that have x = 0.

! "IG_Mid ( i)" is 1 if the entire i^th face has y <= 0 and will contribute
! to the evaluation of the y-flux. It is 0 if the entire i^th face has y > 0
! and will not contribute at all to the evaluation of the y-flux.
! If y ( ip1) <= 0 and y ( i) > 0 then only [ LG_Tau_Mid ( i), L_Tau_Max ( i)]
! will contribute to the y-flux. In that case "IG_Mid ( i) = - 1".
! If y ( ip1) > 0 and y ( i) <= 0 then only [ L_Tau_Min ( i), LG_Tau_Mid ( i)]
! will contribute to the y-flux. In that case "IG_Mid ( i) = - 2".
! If "IG_Mid ( i)" is negative for the i^th face, then "xG_Lo, xG_Hi" will
! contain the x-coordinates of the two faces that have y = 0.

! "Int_Momnt_U ( 10), Int_Momnt_F ( 10), Int_Momnt_G ( 10)" Integration of
! the first 10 moments of x and y over the areas of the wave model that
! contribute to the evaluation of "U", "F" and "G".

! "ifail", it is set to 0 to begin with. There are extreme circumstances
! where this wave model may produce a non-convex shape. I.e. we cannot fully
! guarantee a convex shape for the wave model. If a non-convex shape is found,
! we return with "ifail = 1". In that case, use the larger wave models
! from BDA14.

! --------------------------------------------------

! Three variants of the wave model can be produced. They are listed below
! in order of increasing dissipation:

! VARIANT # 1: In addition to "S_Neta_Max_Ex" to form the largest wave
! speed, we also use "S_Neta_Min_Ex" to form "S_Neta_Multid".
! I.e., the maximization of "S_Neta_Multid" does not only depend on
! "S_Neta_Max_Ex" but rather it also depends on "S_Neta_Min_Ex".
! Uncomment the commented out lines in
! "S_Neta_Multid ( i) = AMAX1 ( S_Neta_Multid ( i)".

! VARIANT # 2: We use the subsonic version of the RS to get all the
! fluxes so that "S_Neta_Multid" is always positive. Uncomment the line,
! "S_Neta_Multid ( 1: nzone_vrtx) = ( expand_wavespeeds - 1.0) * S_Radial_Max".

! VARIANT # 3: We maximize "S_Neta_Multid" by using "S_Radial_Max".
! This gives a multid analogue of LLF, though one can put a contact in like in
! HLLC. Uncomment the line, "S_Neta_Multid ( 1: nzone_vrtx) = S_Radial_Max".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx, multid_wavemodel_num,                         

! Is the wave model supersonic? (output)
     &        Is_SuperSonic, Is_SuperSonic_State, ifail

      REAL smallnum, blend_HLL, expand_wavespeeds,                      

! ALE speed with which the vertex of the mesh moves (input).
     &     Vlx_ALE, Vly_ALE,                                            

! Supersonic speeds in the global frame; always built even if flow isn't
! supersonic. Also provide the maximal speed in this RP. (output)
     &     Vlx_Centroid, Vly_Centroid, S_Radial_Max,                    

! "yF_Lo, yF_Hi" are points of intersection where x = 0 intersects wave model.
     &     yF_Lo, yF_Hi,                                                

! "xG_Lo, xG_Hi" are points of intersection where y = 0 intersects wave model.
     &     xG_Lo, xG_Hi,                                                

! First 10 moments of x and y over areas contributing to "U", "F" and "G".
     &     Int_Momnt_U ( 10), Int_Momnt_F ( 10), Int_Momnt_G ( 10)


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX, NFIELD + 1 + NFLUID) :: 

! Cyclically arranged state variables that come together at a vertex.
     &      U


      INTEGER, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) :: IF_Mid, IG_Mid


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::                      

! Extremal speeds in the Neta- and Tau-directions. Unit vectors in
! principal directions. (inputs)
     &     S_Neta_Min, S_Neta_Max, S_Tau_Min, S_Tau_Max, nx, ny,        

! Polytropic indices on either side.
     &      gamma_eos, gamma_sdsp,                                      

! Fastest speed along the principal direction. Extents of the local
! wavefront. Intersection points of the wave fronts. (outputs)
     &     S_Neta_Multid, L_Tau_Min, L_Tau_Max, x, y,                   

! Unit vector along the edge of a wave model. Dot products of that unit vector
! with "n" and "t" unit vectors.
     &     neta_x, neta_y, neta_dot_n, neta_dot_t,                      

! "LF_Tau_Mid, LG_Tau_Mid" are mid-points of the facial integration.
     &     LF_Tau_Mid, LG_Tau_Mid


      REAL, DIMENSION ( MAX_ZONES_PER_VRTX, NWAVES_1DRS) ::             

! Wave speeds in the 1D Riemann problem; unparametrized and parametrized.
     &     S_1DRS, L_1DRS

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      INTEGER i, j, im1, ip1, iter, ifound,                             
     &        itempaa, itempab, Is_Subsonic

      REAL S_Min, S_Max,                                                
     &     XX_Start, XX_End, YY_Start, YY_End,                          
     &     tempaa, tempab, tempac, tempad,                              
     &     temp_x, temp_y, x1, y1, x2, y2,                              
     &     xcen, ycen, moments_arr ( 10)

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::                      

     &      neta, tau, S_Neta_Multid_Sav,                               

! Slightly expanded version of above-mentioned wave speeds.
     &      S_Neta_Min_Ex, S_Neta_Max_Ex, S_Tau_Min_Ex, S_Tau_Max_Ex

! ----------------------------------------------------------------------
! STEP I) Slightly expand the wavespeeds relative to the centroids of each
! individual Monge cone from each individual one-dimensional Riemann solver.
! The alternative is to keep those speeds unchanged.

! This is eqns. (B.1) and (B.2) of BDA14.

      ifail = 0  ! For now, the wave model is good.

      S_Radial_Max = -1.0e30

      DO i = 1, nzone_vrtx

      tempaa = 0.5 * ( S_Neta_Min ( i) + S_Neta_Max ( i) )
      tempab = 0.5 * ( S_Tau_Min ( i) + S_Tau_Max ( i) )

      S_Neta_Max_Ex ( i) = tempaa                                       
     &            + expand_wavespeeds * ( S_Neta_Max ( i) - tempaa)
      S_Neta_Min_Ex ( i) = tempaa                                       
     &            + expand_wavespeeds * ( S_Neta_Min ( i) - tempaa)

      S_Tau_Max_Ex ( i) = tempab                                        
     &            + expand_wavespeeds * ( S_Tau_Max ( i) - tempab)
      S_Tau_Min_Ex ( i) = tempab                                        
     &            + expand_wavespeeds * ( S_Tau_Min ( i) - tempab)

!      S_Neta_Max_Ex ( i) = S_Neta_Max ( i)
!      S_Neta_Min_Ex ( i) = S_Neta_Min ( i)

!      S_Tau_Max_Ex ( i) = S_Tau_Max ( i)
!      S_Tau_Min_Ex ( i) = S_Tau_Min ( i)

      S_Radial_Max = AMAX1 ( S_Radial_Max,                              
     &               SQRT ( S_Neta_Max ( i)**2 + S_Tau_Max ( i)**2),    
     &               SQRT ( S_Neta_Max ( i)**2 + S_Tau_Min ( i)**2),    
     &               SQRT ( S_Neta_Min ( i)**2 + S_Tau_Max ( i)**2),    
     &               SQRT ( S_Neta_Min ( i)**2 + S_Tau_Min ( i)**2) )

      END DO  ! End of "DO i =" loop.

! ----------------------------------------------------------------------
! STEP II) Evaluate "S_Neta_Multid ( i)" along the "nx ( i), ny ( i)"
! direction. This consists of visiting all directions "nx ( j), ny ( j)"
! and finding the maximum propagation speed along the former direction.

! The process is made possible by the fact that we know the shape and
! orientation of the ellipse in its own local coordinate system.
! We know the shape of the rectangle that circumscribes the ellipse.
! By transforming into that system, we can find the distance we need to
! travel along the unit normal "neta, tau" to hit the boundary of
! that rectangle.
! ----------------------------------------------------------------------

! n(j) = ( n(i) dot Neta (j)) Neta (j) + ( n(i) dot Tau (j)) Tau (j)
! is a projection that gives us the normal along the i^th principal
! direction in the ( Neta (j), Tau (j)) frame.

! Also make min and max signal speeds that one can use later on to
! ensure that the wave model is good.

! ------------------------------------------------------------
! First Sub-Step: Form "S_Neta_Multid".

! This is eqn. (B.3) of BDA14.

      S_Neta_Multid ( 1: nzone_vrtx) = -1.0e30

      IF ( multid_wavemodel_num == 2) THEN

! VARIANT # 2: Use line below so that all zones are treated subsonically.
        S_Neta_Multid ( 1: nzone_vrtx) = ( expand_wavespeeds - 1.0)     
     &                                 * S_Radial_Max

      ELSE IF ( multid_wavemodel_num == 3) THEN

! VARIANT # 3: Use line below to get something like a multid LLF RS.
        S_Neta_Multid ( 1: nzone_vrtx) = S_Radial_Max

      END IF

      S_Min = 1.0e30
      S_Max = -1.0e30

! --------------------

! This is eqns. (B.4) and (B.5) of BDA14.

      DO i = 1, nzone_vrtx

! ----------

        IF ( multid_wavemodel_num == 0) THEN

! ----------

          DO j = 1, nzone_vrtx

          neta ( j) = nx ( i) * nx ( j) + ny ( i) * ny ( j)
          tau ( j) = - nx ( i) * ny ( j) + ny ( i) * nx ( j)

          S_Neta_Multid ( i) = AMAX1 ( S_Neta_Multid ( i)               
     &  , S_Neta_Max_Ex ( j) * neta ( j) + S_Tau_Max_Ex ( j) * tau ( j) 
     &  , S_Neta_Max_Ex ( j) * neta ( j) + S_Tau_Min_Ex ( j) * tau ( j) 
     &                                                                 )

          END DO ! End of "DO j =".

! ----------

        ELSE

! ----------
! VARIANT # 1: Strictly speaking, the "S_Neta_Min_Ex ( j)" does not need
! to contribute to the expansion of "S_Neta_Multid ( i)". Use if needed.

          DO j = 1, nzone_vrtx

          neta ( j) = nx ( i) * nx ( j) + ny ( i) * ny ( j)
          tau ( j) = - nx ( i) * ny ( j) + ny ( i) * nx ( j)

          S_Neta_Multid ( i) = AMAX1 ( S_Neta_Multid ( i)               
     &  , S_Neta_Max_Ex ( j) * neta ( j) + S_Tau_Max_Ex ( j) * tau ( j) 
     &  , S_Neta_Max_Ex ( j) * neta ( j) + S_Tau_Min_Ex ( j) * tau ( j) 
     &  , S_Neta_Min_Ex ( j) * neta ( j) + S_Tau_Max_Ex ( j) * tau ( j) 
     &  , S_Neta_Min_Ex ( j) * neta ( j) + S_Tau_Min_Ex ( j) * tau ( j) 
     &                                                                 )

          END DO ! End of "DO j =".

! ----------

        END IF  ! End of "IF ( multid_wavemodel_num == 0)".

! ----------

      S_Min = AMIN1 ( S_Min, ABS ( S_Neta_Multid ( i) ) )
      S_Max = AMAX1 ( S_Max, ABS ( S_Neta_Multid ( i) ) )

      END DO ! End of "DO i =".

! ------------------------------------------------------------
! Second Sub-Step: Iteratively make "S_Neta_Multid" converge to its largest
! value.

! This is eqn. (B.7) of the BDA14.

      DO iter = 1, 3 * nzone_vrtx

! ------------------------------------------------------------

      S_Neta_Multid_Sav ( 1: nzone_vrtx)=S_Neta_Multid ( 1: nzone_vrtx)

! --------------------------------------------------

      DO i = 1, nzone_vrtx

        DO j = 1, nzone_vrtx

        neta ( j) = nx ( i) * nx ( j) + ny ( i) * ny ( j)
        tau ( j) = - nx ( i) * ny ( j) + ny ( i) * nx ( j)

        S_Neta_Multid ( i) = AMAX1 ( S_Neta_Multid ( i),                
     &   S_Neta_Multid ( j) * neta ( j) + S_Tau_Max_Ex ( j) * tau ( j), 
     &   S_Neta_Multid ( j) * neta ( j) + S_Tau_Min_Ex ( j) * tau ( j) )

        END DO ! End of "DO j =".

      S_Radial_Max = AMAX1 ( S_Radial_Max, S_Neta_Multid ( i),          
     &              SQRT ( S_Neta_Multid ( i)**2 + S_Tau_Max ( i)**2),  
     &              SQRT ( S_Neta_Multid ( i)**2 + S_Tau_Min ( i)**2) )

      END DO ! End of "DO i =".

! --------------------------------------------------
! Now check for convergence. If we have converged, we can exit the
! cycle of iterations.

!      itempaa = 0

!      DO i = 1, nzone_vrtx
!      IF ( ABS ( S_Neta_Multid ( i) - S_Neta_Multid_Sav ( i) ) >=       
!     &                         smallnum * S_Radial_Max ) itempaa = 1
!      END DO

!      IF ( itempaa == 0) EXIT

! --------------------------------------------------
! Tentatively build ( x, y) for each "S_Neta_Multid".

      DO i = 1, nzone_vrtx

      im1 = i - 1
      IF ( im1 < 1) im1 = nzone_vrtx

      tempaa = 1.0 / ( nx ( im1) * ny ( i) - ny ( im1) * nx ( i))

      x ( i) = tempaa * ( S_Neta_Multid ( im1) * ny ( i)                
     &                  - S_Neta_Multid ( i) * ny ( im1) )

      y ( i) = - tempaa * ( S_Neta_Multid ( im1) * nx ( i)              
     &                    - S_Neta_Multid ( i) * nx ( im1) )

      END DO

! ------------------------------
! Tentatively build "L_Tau_Min, L_Tau_Max" so that we can ensure that they
! are within bounds. If they are, we can exit the cycle of iterations.

      itempaa = 0

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

      L_Tau_Min ( i) = - x ( i) * ny ( i) + y ( i) * nx ( i)

      L_Tau_Max ( i) = - x ( ip1) * ny ( i) + y ( ip1) * nx ( i)

      tempaa = smallnum * S_Radial_Max

      IF ( ( L_Tau_Min ( i) > S_Tau_Min ( i) + tempaa) .OR.             
     &     ( S_Tau_Min ( i) > S_Tau_Max ( i) + tempaa) .OR.             
     &     ( S_Tau_Max ( i) > L_Tau_Max ( i) + tempaa) ) itempaa = 1

      END DO

      IF ( itempaa == 0) EXIT

! ------------------------------------------------------------

      END DO ! End of "DO iter =".

! ------------------------------------------------------------
! Stop if convergence is not obtained. Actually, good to do this test only
! if subsonic. For this reason, it is deferred to the end of the next step.

!      IF ( itempaa == 1) THEN
!      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :",           
!     &             " S_Neta_Multid did not converge"
!      STOP
!      END IF

! ----------------------------------------------------------------------
! STEP III) Decide whether the Riemann Problem is subsonic or supersonic
! in a multidimensional sense.
! Also set up "Vlx_Centroid, Vly_Centroid".

! If it is subsonic then we can give ourselves the option of inflating the
! wave speeds by some amount so that the waves don't touch the origin. This
! can even be done for subsonic and supersonic cases to yield the
! multidimensional analogue of the LLF Riemann solver.
! ----------------------------------------------------------------------

      Is_Subsonic = 1

      Vlx_Centroid = 0.0
      Vly_Centroid = 0.0

      DO i = 1, nzone_vrtx

      neta ( i) = 0.5 * ( S_Neta_Max ( i) + S_Neta_Min ( i) )
      tau ( i) = 0.5 * ( S_Tau_Max ( i) + S_Tau_Min ( i) )

      IF ( S_Neta_Multid ( i) <= 0.0) Is_Subsonic = 0

      Vlx_Centroid = Vlx_Centroid                                       
     &               + neta ( i) * nx ( i) - tau ( i) * ny ( i)
      Vly_Centroid = Vly_Centroid                                       
     &               + neta ( i) * ny ( i) + tau ( i) * nx ( i)

      END DO

      Vlx_Centroid = Vlx_Centroid / REAL ( nzone_vrtx)
      Vly_Centroid = Vly_Centroid / REAL ( nzone_vrtx)

! -----

! Stop if convergence was not obtained in the subsonic cases.

      IF ( ( Is_Subsonic == 0) .AND. ( itempaa == 1) ) THEN
      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :",           
     &             " S_Neta_Multid did not converge"
#if ( 1 == 1)
      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"---------- i = ", i," -------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"nx, ny = ", nx ( i), ny ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U = ", U ( i, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma = ", gamma_eos ( i), gamma_sdsp ( i)
      WRITE ( 6, *)" "

      END DO
#endif
      STOP
      END IF

! -----

! Stop if wave model is supersonic but has a vanishingly small supersonic
! flow speed.

      IF ( ( Is_Subsonic == 0) .AND.                                    
     &     ( SQRT ( Vlx_Centroid**2 + Vly_Centroid**2)                  
     &                                         < smallnum) ) THEN
      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :"
      WRITE ( 6, *)"supersonic speeds too small ", Vlx_Centroid,        
     &                             Vly_Centroid, S_Radial_Max
      STOP
      END IF

! ------------------------------
! Use "S_Min, S_Max, S_Radial_Max" to increase "S_Neta_Multid ( i)"
! if needed. This ensures that the areas will be positive in the
! subsonic case. Expand all the waves by 0.01% as a safety measure.

! For 1d HLL, one can often treat the supersonic case as marginally subsonic.
! The same can be accomplished for multid-HLL by commenting out
! "IF ( Is_Subsonic == 1) THEN" and "END IF" below.

      IF ( Is_Subsonic == 1) THEN

      S_Radial_Max = expand_wavespeeds * S_Radial_Max

! The coefficient in front of "S_Radial_Max" determines how much extra
! dissipation is introduced; with a coefficient of unity being the
! maximum amount of dissipation.
      tempaa = AMAX1 ( (expand_wavespeeds - 1.0) * ABS ( S_Radial_Max), 
     &                 100.0 * smallnum)

      DO i = 1, nzone_vrtx

      S_Neta_Multid ( i) =                                              
     &       AMAX1 ( expand_wavespeeds * S_Neta_Multid ( i), tempaa)

      END DO

      END IF

! ------------------------------
! USEFUL_CHECK) If one wants the multidimensional analogue of an LLF
! flux, then use the below code. This doesn't work because it introduces
! downstream coupling.

#if ( 0 == 1)
      Is_Subsonic = 1
      S_Neta_Multid ( 1: nzone_vrtx) = ABS ( S_Radial_Max)
#endif

! ----------------------------------------------------------------------
! STEP IV) Find "x, y", i.e. the locations where the wave fronts in each
! of the principal directions intersect. This is eqn. (2.4) of the paper.

! At least for subsonic cases, we do insist on convexity of the
! wave model. Do all the tests imaginable to prove that it is convex.
! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

      im1 = i - 1
      IF ( im1 < 1) im1 = nzone_vrtx

      tempaa = 1.0 / ( nx ( im1) * ny ( i) - ny ( im1) * nx ( i))

      x ( i) = tempaa * ( S_Neta_Multid ( im1) * ny ( i)                
     &                  - S_Neta_Multid ( i) * ny ( im1) )

      y ( i) = - tempaa * ( S_Neta_Multid ( im1) * nx ( i)              
     &                    - S_Neta_Multid ( i) * nx ( im1) )

      END DO

! ------------------------------
! If "blend_HLL" is large and if the problem is subsonic, one may reset
! "S_Tau_Min, S_Tau_Max" to increase dissipation.
! May only be needed for extreme problems and even there it is iffy.

#if ( EXTREMIZE_WAVESPEEDS >= 1 && 0 == 1)
      IF ( ( Is_Subsonic == 1) .AND. ( blend_HLL >= 0.5) ) THEN

        DO i = 1, nzone_vrtx

        ip1 = i + 1
        IF ( ip1 > nzone_vrtx) ip1 = 1

        S_Tau_Max ( i) = AMAX1 ( S_Tau_Max ( i),                        
     &                     - x ( ip1) * ny ( i) + y ( ip1) * nx ( i) )

        S_Tau_Min ( i) = AMIN1 ( S_Tau_Min ( i),                        
     &                     - x ( i) * ny ( i) + y ( i) * nx ( i) )

        END DO

      END IF
#endif

! ----------------------------------------------------------------------
! STEP V) Now that we have endpoints of each wavefront "(x, y)" and
! "S_Neta_Multid", we can find "L_Tau_Min, L_Tau_Max".
! We just project "(x, y)" on to "Tau = ( -ny, nx)".
! ----------------------------------------------------------------------

! USEFUL_CHECK) Also check to make sure that we have the ordered inequality:
! "L_Tau_Min <= S_Tau_Min <= S_Tau_Max <= L_Tau_Max"; with some margins.
! Since "S_Neta_Multid" was marginally inflated at the end of STEP III, the
! margins are not needed.

      itempaa = 0

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

!      L_Tau_Min ( i) =                                                  
!     &        - ( x ( i) - S_Neta_Multid ( i) * nx ( i) ) * ny ( i)     
!     &        + ( y ( i) - S_Neta_Multid ( i) * ny ( i) ) * nx ( i)

!      L_Tau_Max ( i) =                                                  
!     &        - ( x ( ip1) - S_Neta_Multid ( i) * nx ( i) ) * ny ( i)   
!     &        + ( y ( ip1) - S_Neta_Multid ( i) * ny ( i) ) * nx ( i)

      L_Tau_Min ( i) = - x ( i) * ny ( i) + y ( i) * nx ( i)

      L_Tau_Max ( i) = - x ( ip1) * ny ( i) + y ( ip1) * nx ( i)

      tempaa = smallnum * S_Radial_Max

      IF ( ( L_Tau_Min ( i) > S_Tau_Min ( i) + tempaa) .OR.             
     &     ( S_Tau_Min ( i) > S_Tau_Max ( i) + tempaa) .OR.             
     &     ( S_Tau_Max ( i) > L_Tau_Max ( i) + tempaa) ) itempaa = 1

!      IF ( ( L_Tau_Min ( i) > S_Tau_Min ( i) ) .OR.                     
!     &     ( S_Tau_Min ( i) > S_Tau_Max ( i) ) .OR.                     
!     &     ( S_Tau_Max ( i) > L_Tau_Max ( i) ) ) itempaa = 1

      END DO

! -----

      IF ( ( Is_Subsonic == 0) .AND. ( itempaa == 1) ) THEN
      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :"
      WRITE ( 6, *)" L_Tau_Min, L_Tau_Max don't bound all the",         
     &             " wave speeds in each face."
#if ( 1 == 1)
      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"---------- i = ", i," -------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"nx, ny = ", nx ( i), ny ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U = ", U ( i, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma = ", gamma_eos ( i), gamma_sdsp ( i)
      WRITE ( 6, *)" "

      END DO
#endif
      STOP
      END IF

! ----------------------------------------------------------------------
! STEP VI) Set "ifail" if the wave model is non-convex.
! ----------------------------------------------------------------------

! TEST I) The following test just makes sure that the intersection
! points of all the wavefronts are not tangled. The assumption is that
! the original unit normals "Neta = ( nx, ny)" are cyclically arranged
! with an anti-clockwise orientation. As a result, the intersection
! points between the wavefronts associated with those normals should
! also follow that same cyclic arrangement.

! We check that the cross product (r (i) - r ( j)) X (r (ip1) - r ( j))
! is positive. If it has become negative, the wave model has become non-convex.

      tempaa = 1.0e30

      DO j = 1, nzone_vrtx

        DO i = 1, nzone_vrtx

        ip1 = i + 1
        IF ( ip1 > nzone_vrtx) ip1 = 1

        IF ( ( j == i) .OR. ( j == ip1) ) CYCLE

        tempaa = AMIN1 ( tempaa,                                        
     &                   ( x ( i) - x ( j)) * ( y ( ip1) - y ( j))      
     &                 - ( x ( ip1) - x ( j)) * ( y ( i) - y ( j)) )

        END DO ! End of "DO i =" loop.

      END DO

      IF ( tempaa <= 0.0) THEN

      ifail = 1

      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :",           
     &        " intersection points in wave model have become tangled."
#if ( 1 == 1)
      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"---------- i = ", i," -------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"nx, ny = ", nx ( i), ny ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U = ", U ( i, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma = ", gamma_eos ( i), gamma_sdsp ( i)
      WRITE ( 6, *)" "

      END DO
#endif
      STOP
      END IF

! ------------------------------

! TEST II) Now check the cross product (r (ip1) - r ( i)) X (r (im1) - r ( i)).
! If it is at all negative, or zero, it means that the wave model has
! become non-convex.

      tempaa = 1.0e30

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

      im1 = i - 1
      IF ( im1 < 1) im1 = nzone_vrtx

      tempaa = AMIN1 ( tempaa,                                          
     &                 ( x ( ip1) - x ( i)) * ( y ( im1) - y ( i))      
     &               - ( y ( ip1) - y ( i)) * ( x ( im1) - x ( i)) )

      END DO ! End of "DO i =" loop.

      IF ( tempaa <= 0.0) THEN

      ifail = 1

      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :",           
     &        " wave model has become non-convex."
#if ( 1 == 1)
      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"---------- i = ", i," -------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"nx, ny = ", nx ( i), ny ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U = ", U ( i, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"gamma = ", gamma_eos ( i), gamma_sdsp ( i)
      WRITE ( 6, *)" "

      END DO
#endif
      STOP
      END IF

! ----------------------------------------------------------------------
! STEP VII) Build the extra ( neta_x ( i), neta_y ( i)) which is the
! outward-pointing unit vector that connects vertex ( x ( i), y ( i))
! to ( x ( ip1), y ( ip1)).

! For this type of wave model, the dot products between
! ( neta_x ( i), neta_y ( i)) and unit vector n = ( nx ( i), ny ( i))
! or unit vector t = ( - ny ( i), nx ( i)) is trivial.
! Can only differ by a sign.

! In this code, "L_Tau_Min ( i) = 0.0" and "L_Tau_Max ( i)" is just the
! length of that side of the wave model. I.e. it just parametrizes the length
! of that side.
! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

! ----------

      neta_x ( i) = y ( ip1) - y ( i)
      neta_y ( i) = x ( i) - x ( ip1)

      tempaa = 1.0 / SQRT ( neta_x ( i)**2 + neta_y ( i)**2)

      neta_x ( i) = neta_x ( i) * tempaa
      neta_y ( i) = neta_y ( i) * tempaa

! ----------

      S_Neta_Multid ( i) = 0.5 * ( x ( i) + x ( ip1) ) * neta_x ( i)    
     &                   + 0.5 * ( y ( i) + y ( ip1) ) * neta_y ( i)

! ----------

      neta_dot_n ( i) = neta_x ( i) * nx ( i) + neta_y ( i) * ny ( i)
      neta_dot_t ( i) = - neta_x ( i) * ny ( i) + neta_y ( i) * nx ( i)

! ----------

      L_Tau_Min ( i) = 0.0
      L_Tau_Max ( i) = SQRT ( (y ( ip1) - y ( i))**2                    
     &                      + (x ( ip1) - x ( i))**2)

! ----------

      END DO

! ----------------------------------------------------------------------
! STEP VIII) In each side "i" we seek the intersection of the ith wavefront:

! neta_x ( i) * x + neta_y ( i) * y = neta_x ( i) * x ( i) + neta_y ( i) * y ( i)

! with the wavefront of the jth wave in the ith 1D RP.

! - ny ( i) * x + nx ( i) * y = S_1DRS ( i, j)

! We then use the intersection point to find "L_1DRS ( i, j)".

! For this plan to suceed, the angle between ( neta_x ( i), neta_y ( i)) and
! ( nx ( i), ny ( i)) must lie between -90 and 90 degrees.
! Because we set ( neta_x ( i), neta_y ( i)) = ( nx ( i), ny ( i)) in the
! previous step, it is assured for this choice of wave model.

! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

      DO j = 1, NWAVES_1DRS

! ----------

      tempaa = neta_x ( i) * x ( i) + neta_y ( i) * y ( i)

      tempab = neta_x ( i) * nx ( i) + neta_y ( i) * ny ( i)

      tempab = 1.0 / tempab

      temp_x = tempab * ( tempaa * nx ( i)                              
     &                  - S_1DRS ( i, j) * neta_y ( i) )

      temp_y = tempab * ( S_1DRS ( i, j) * neta_x ( i)                  
     &                  + tempaa * ny ( i) )

      L_1DRS ( i, j) = SQRT ( ( temp_x - x ( i))**2                     
     &                      + ( temp_y - y ( i))**2 )

! ----------

      END DO  ! End of "DO j =" loop.

      END DO  ! End of "DO i =" loop.

! ----------------------------------------------------------------------
! STEP IX) Find "Vlx_Centroid, Vly_Centroid". Eventually we will shift all
! "( x ( i), y ( i))" to be centered around that centroid.
! ----------------------------------------------------------------------

! First find arithmetic averages of the vertices, "temp_x, temp_y".
! This point is sure to lie within the wave model and serves as a temporary
! origin relative to which the centroid calculations are done.

      temp_x = 0.0
      temp_y = 0.0

      DO i = 1, nzone_vrtx

      temp_x = temp_x + x ( i)
      temp_y = temp_y + y ( i)

      END DO

      temp_x = temp_x / nzone_vrtx
      temp_y = temp_y / nzone_vrtx

! ----------

      tempaa = 0.0

      Vlx_Centroid = 0.0
      Vly_Centroid = 0.0

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

      x1 = x ( i) - temp_x
      y1 = y ( i) - temp_y
      x2 = x ( ip1) - temp_x
      y2 = y ( ip1) - temp_y

      CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,                    
     &                               xcen, ycen, moments_arr )

      tempaa = tempaa + moments_arr ( 1)
      Vlx_Centroid = Vlx_Centroid + moments_arr ( 2)
      Vly_Centroid = Vly_Centroid + moments_arr ( 3)

      END DO

      Vlx_Centroid = temp_x + Vlx_Centroid / tempaa
      Vly_Centroid = temp_y + Vly_Centroid / tempaa

! ----------
! Make sure that the outward pointing normals are truly pointing outward of
! the wave model

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

      temp_x = 0.5 * ( x ( i) + x ( ip1)) - Vlx_Centroid
      temp_y = 0.5 * ( y ( i) + y ( ip1)) - Vly_Centroid

      IF ( neta_x ( i) * temp_x + neta_y ( i) * temp_y <= 0.0) ifail = 1

      END DO

      IF ( ifail == 1) THEN
      WRITE ( 6, *)"STOPPING in MULTID_WAVE_MODEL because :",           
     &        " Normals are not outward pointing w.r.t. centroid."
      STOP
      END IF

! ----------------------------------------------------------------------
! STEP X) The centroid is still at ( Vlx_Centroid, Vly_Centroid) 
! system. The vertex of the mesh is at ( Vlx_ALE, Vly_ALE).
! If these two points lie on opposite sides of any of the edges of the
! multidimensional wave model then we are supersonic. Otherwise, we are
! subsonic.

! Because of the "<" in the IF conditional, the boundary of the wave model
! is excluded from being a supersonic case. Might as well give it a
! stabilizing subsonic solution if it is on the edge.

! The boundary which has supersonic flow will have its number returned
! in "Is_SuperSonic". Otherwise, it will return 0 for subsonic.
! Also do the computational geometry analysis to find the exact state
! "Is_SuperSonic_State" within that boundary; range [ 0, NWAVES_1DRS].

! At the end of this step, we can shift ( x ( i), y ( i)) so that the
! wave model is centered on the centroid "( Vlx_Centroid, Vly_Centroid)".
! ----------------------------------------------------------------------

      Is_SuperSonic = 0
      Is_SuperSonic_State = - 1000000  ! Dummy value.

      ifound = 0

! ------------------------------

      SIDE: DO i = 1, nzone_vrtx  ! Named do loop with name "SIDE".

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

      im1 = i - 1
      IF ( im1 < 1) im1 = nzone_vrtx

! ------------------------------

      tempaa = neta_x ( i) * ( Vlx_Centroid - x ( i) )                  
     &       + neta_y ( i) * ( Vly_Centroid - y ( i) )

      tempab = neta_x ( i) * ( Vlx_ALE - x ( i) )                       
     &       + neta_y ( i) * ( Vly_ALE - y ( i) )

! --------------------

      IF ( tempaa * tempab < 0.0) THEN

      ifound = 1

! --------------------

        tempac = neta_x ( im1) * ( Vly_ALE - y ( i))                    
     &         - neta_y ( im1) * ( Vlx_ALE - x ( i))

        tempad = neta_y ( i) * ( Vlx_ALE - x ( i))                      
     &         - neta_x ( i) * ( Vly_ALE - y ( i))

        IF ( ( tempac >= 0.0) .AND. ( tempad >= 0.0) ) THEN

          Is_SuperSonic = i
          Is_SuperSonic_State = 0

          EXIT SIDE

        END IF

! ----------

        tempac = - ny ( i) * ( Vlx_ALE - x ( i))                        
     &           + nx ( i) * ( Vly_ALE - y ( i))

        tempad = - ny ( i) * Vlx_ALE + nx ( i) * Vly_ALE                
     &         - S_1DRS ( i, 1)

        IF ( tempac * tempad <= 0.0) THEN

          Is_SuperSonic = i
          Is_SuperSonic_State = 0

          EXIT SIDE

        END IF

! ----------

        tempac = - ny ( i) * ( Vlx_ALE - x ( ip1))                      
     &           + nx ( i) * ( Vly_ALE - y ( ip1))

        tempad = - ny ( i) * Vlx_ALE + nx ( i) * Vly_ALE                
     &         - S_1DRS ( i, NWAVES_1DRS)

        IF ( tempac * tempad <= 0.0) THEN

          Is_SuperSonic = i
          Is_SuperSonic_State = NWAVES_1DRS

          EXIT SIDE

        END IF

! ----------

        tempac = neta_x ( i) * ( Vly_ALE - y ( ip1))                    
     &         - neta_y ( i) * ( Vlx_ALE - x ( ip1))

        tempad = neta_y ( ip1) * ( Vlx_ALE - x ( ip1))                  
     &         - neta_x ( ip1) * ( Vly_ALE - y ( ip1))

        IF ( ( tempac >= 0.0) .AND. ( tempad >= 0.0) ) THEN

          Is_SuperSonic = i
          Is_SuperSonic_State = NWAVES_1DRS

          EXIT SIDE

        END IF

! ----------
! Scan the intermediate waves in the 1D RP.

        DO j = 1, NWAVES_1DRS - 1

        tempac = - ny ( i) * Vlx_ALE + nx ( i) * Vly_ALE                
     &         - S_1DRS ( i, j)

        tempad = - ny ( i) * Vlx_ALE + nx ( i) * Vly_ALE                
     &         - S_1DRS ( i, j + 1)

        IF ( tempac * tempad <= 0.0) THEN

          Is_SuperSonic = i
          Is_SuperSonic_State = j

          EXIT SIDE

        END IF

        END DO

! --------------------

      END IF  ! End of "IF ( tempaa * tempab < 0.0)".

! ------------------------------

      END DO SIDE  ! End of "DO i =" loop.

! ------------------------------
! Should not detect a situation where wave model is supersonic but more
! detailed search did not find a state.

      IF ( ( ifound == 1) .AND. ( Is_SuperSonic == 0) ) THEN
      WRITE ( 6, *)"Stopping in MULTID_WAVE_MODEL because: "
      WRITE ( 6, *)"inconsisten values : ifound, Is_SuperSonic = ",     
     &              ifound, Is_SuperSonic
      STOP
      END IF

! ------------------------------------------------------------

! Now shift all ( x ( i), y ( i)) so that ( Vlx_Centroid, Vly_Centroid) is in
! the center of the wave model. We can always shift back as and when needed.

      DO i = 1, nzone_vrtx

      x ( i) = x ( i) - Vlx_Centroid
      y ( i) = y ( i) - Vly_Centroid

      END DO

! ----------------------------------------------------------------------
! STEP XI) Build "IF_Mid, LF_Tau_Mid, yF_Mid, IG_Mid, LG_Tau_Mid, xG_Mid".
! These extra variables are very helpful in carrying out the facial integrals
! in the faces of the multidimensional wave model.

! The wave model must be shifted to the centroid before this step can be done.
! ----------------------------------------------------------------------

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

! ----------

      IF ( ( x ( ip1) <= 0.0) .AND. ( x ( i) <= 0.0) ) THEN
! In the i^th face of wave model, integrate "[ L_Tau_Min ( i), L_Tau_Max ( i)]"

        IF_Mid ( i) = 1
        LF_Tau_Mid ( i) = 0.0  ! Dummy value; immaterial here.

      ELSE IF ( ( x ( ip1) <= 0.0) .AND. ( x ( i) > 0.0) ) THEN
! In the i^th face of wave model, integrate "[ LF_Tau_Mid ( i), L_Tau_Max ( i)]"

        IF_Mid ( i) = - 1
        LF_Tau_Mid ( i) = x ( i) / neta_y ( i)
        yF_Hi = y ( i) + neta_x ( i) * LF_Tau_Mid ( i)

      ELSE IF ( ( x ( ip1) > 0.0) .AND. ( x ( i) <= 0.0) ) THEN
! In the i^th face of wave model, integrate "[ L_Tau_Min ( i), LF_Tau_Mid ( i)]"

        IF_Mid ( i) = - 2
        LF_Tau_Mid ( i) = x ( i) / neta_y ( i)
        yF_Lo = y ( i) + neta_x ( i) * LF_Tau_Mid ( i)

      ELSE
! The i^th face of wave model does not contribute to the x-flux.

        IF_Mid ( i) = 0
        LF_Tau_Mid ( i) = 0.0  ! Dummy value; immaterial here.

      END IF

! ----------

      IF ( ( y ( ip1) <= 0.0) .AND. ( y ( i) <= 0.0) ) THEN
! In the i^th face of wave model, integrate "[ L_Tau_Min ( i), L_Tau_Max ( i)]"

        IG_Mid ( i) = 1
        LG_Tau_Mid ( i) = 0.0  ! Dummy value; immaterial here.

      ELSE IF ( ( y ( ip1) <= 0.0) .AND. ( y ( i) > 0.0) ) THEN
! In the i^th face of wave model, integrate "[ LG_Tau_Mid ( i), L_Tau_Max ( i)]"

        IG_Mid ( i) = - 1
        LG_Tau_Mid ( i) = - y ( i) / neta_x ( i)
        xG_Lo = x ( i) - neta_y ( i) * LG_Tau_Mid ( i)

      ELSE IF ( ( y ( ip1) > 0.0) .AND. ( y ( i) <= 0.0) ) THEN
! In the i^th face of wave model, integrate "[ L_Tau_Min ( i), LG_Tau_Mid ( i)]"

        IG_Mid ( i) = - 2
        LG_Tau_Mid ( i) = - y ( i) / neta_x ( i)
        xG_Hi = x ( i) - neta_y ( i) * LG_Tau_Mid ( i)

      ELSE
! The i^th face of wave model does not contribute to the y-flux.

        IG_Mid ( i) = 0
        LG_Tau_Mid ( i) = 0.0  ! Dummy value; immaterial here.

      END IF

! ----------

      END DO

! ----------------------------------------------------------------------
! STEP XII) Evaluate "Int_Momnt_U ( 10), Int_Momnt_F ( 10), Int_Momnt_G ( 10)"
! These are the integrals of the first 10 moments of x and y over the areas
! of the wave model that contribute to the evaluation of "U", "F" and "G".

! The wave model must be shifted to the centroid before this step can be done.
! ----------------------------------------------------------------------

      Int_Momnt_U ( :) = 0.0
      Int_Momnt_F ( :) = 0.0
      Int_Momnt_G ( :) = 0.0

      DO i = 1, nzone_vrtx

      ip1 = i + 1
      IF ( ip1 > nzone_vrtx) ip1 = 1

! ----------
! Build moments over the entire wave model.
! The contributions that are common with the x <= 0 part or y <= 0 part are
! added in.

        x1 = x ( i)
        y1 = y ( i)
        x2 = x ( ip1)
        y2 = y ( ip1)

        CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,                  
     &                                 xcen, ycen, moments_arr )

        Int_Momnt_U ( :) = Int_Momnt_U ( :) + moments_arr ( :)

        IF ( IF_Mid ( i) == 1)                                          
     &  Int_Momnt_F ( :) = Int_Momnt_F ( :) + moments_arr ( :)

        IF ( IG_Mid ( i) == 1)                                          
     &  Int_Momnt_G ( :) = Int_Momnt_G ( :) + moments_arr ( :)

! ----------
! Build moments over the x <= 0 part of the wave model. The common parts are
! added in.

      IF ( IF_Mid ( i) == - 1) THEN

        x1 = 0.0
        y1 = yF_Hi
        x2 = x ( ip1)
        y2 = y ( ip1)

        CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,                  
     &                                 xcen, ycen, moments_arr )

        Int_Momnt_F ( :) = Int_Momnt_F ( :) + moments_arr ( :)

      ELSE IF ( IF_Mid ( i) == - 2) THEN

        x1 = x ( i)
        y1 = y ( i)
        x2 = 0.0
        y2 = yF_Lo

        CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,                  
     &                                 xcen, ycen, moments_arr )

        Int_Momnt_F ( :) = Int_Momnt_F ( :) + moments_arr ( :)

      END IF

! ----------
! Build moments over the y <= 0 part of the wave model. The common parts are
! added in.

      IF ( IG_Mid ( i) == - 1) THEN

        x1 = xG_Lo
        y1 = 0.0
        x2 = x ( ip1)
        y2 = y ( ip1)

        CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,                  
     &                                 xcen, ycen, moments_arr )

        Int_Momnt_G ( :) = Int_Momnt_G ( :) + moments_arr ( :)

      ELSE IF ( IG_Mid ( i) == - 2) THEN

        x1 = x ( i)
        y1 = y ( i)
        x2 = xG_Hi
        y2 = 0.0

        CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,                  
     &                                 xcen, ycen, moments_arr )

        Int_Momnt_G ( :) = Int_Momnt_G ( :) + moments_arr ( :)

      END IF

! ----------

      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE MULTID_WAVE_MODEL

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,
     1           xcen, ycen, moments_arr )

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Three points in the xy-plane are specified cyclincally. They are:
! O the origin, P1 ( x1, y1) and P2 ( x2, y2). These are the inputs.

! The outputs are : (xcen, ycen) the centroid and the vector of moments
! given by "moments_arr ( 10)". The elements of the array are the following:
! moments_arr ( 1) -- zeroth moment
! moments_arr ( 2: 3) -- First moment in x and y.
! moments_arr ( 4: 6) -- x^2 moment; x * y moment; y^2 moment.
! moments_arr ( 7: 10) -- x^3 moment; x^2 * y moment; x * y^2 moment;
!                         y^3 moment

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL x1, y1, x2, y2, xcen, ycen, moments_arr ( 10)

! ----------------------------------------------------------------------

      REAL ct, st, x1p, x2p, y2p, p2, p12,
     1     ob3, ob6, ob12, ob15, ob20,
     1     tempaa, tempab, tempac, tempad

! ----------------------------------------------------------------------

      ob3 = 1.0 / 3.0
      ob6 = 1.0 / 6.0
      ob12 = 1.0 / 12.0
      ob15 = 1.0 / 15.0
      ob20 = 1.0 / 20.0

! STEP I) Form "ct, st" and rotate the axes so that we have segment OP1
! along the new x-axes. Then find coordinates of P1 and P2 in the rotated
! frame as (x1p, 0.0) and (x2p, y2p).
! Also find centroid.
! Also build reciprocals of the slopes in the rotated frame.

      tempaa = x1**2 + y1**2

      IF ( tempaa <= 1.0e-30) THEN
      WRITE ( 6, *)"Stopping in XY_MOMENTS_OF_TRIANGLES : ",
     1      " Could not rotate the axes; tempaa = ", tempaa
      STOP
      END IF

      tempaa = 1.0 / SQRT ( tempaa)
      ct = x1 * tempaa
      st = y1 * tempaa

! Find the coordinates of P1 and P2 in the rotated frame.

      x1p = ct * x1 + st * y1

      x2p = ct * x2 + st * y2
      y2p = - st * x2 + ct * y2

! If "y2p" is zero or negative, the points were not cyclic or the area of
! the triangle is zero. Either way, we are in trouble.

      IF ( y2p <= 1.0e-30) THEN
      WRITE ( 6, *)"Stopping in XY_MOMENTS_OF_TRIANGLES : ",
     1      " Points not cyclic or area zero; y2p = ", y2p
      STOP
      END IF

! ------------------------------

      xcen = ob3 * ( x1 + x2)
      ycen = ob3 * ( y1 + y2)

! ------------------------------

      p2 = x2p / y2p

      p12 = ( x2p - x1p) / y2p

! ----------------------------------------------------------------------
! STEP II) Find the zeroth and first moments.

      moments_arr ( 1) = x1p * y2p + 0.5 * p12 * y2p**2
     1                 - 0.5 * p2 * y2p**2

! ------------------------------

! x-moment in rotated frame.

      tempaa = 0.5 * x1p**2 * y2p + 0.5 * p12 * x1p * y2p**2
     1       + ob6 * p12**2 * y2p**3 - ob6 * p2**2 * y2p**3

! y-moment in rotated frame.

      tempab = 0.5 * x1p * y2p**2 + ob3 * p12 * y2p**3
     1       - ob3 * p2 * y2p**3

! ------------------------------

! First x-moment in non-rotated frame.

      moments_arr ( 2) = ct * tempaa - st * tempab

! First y-moment in non-rotated frame.

      moments_arr ( 3) = st * tempaa + ct * tempab

! ----------------------------------------------------------------------
! STEP III) Find the second moments.

! x^2-moment in rotated frame.

      tempaa = ob3 * x1p**3 * y2p + 0.5 * p12 * x1p**2 * y2p**2
     1       + ob3 * p12**2 * x1p * y2p**3 + ob12 * p12**3 * y2p**4
     1       - ob12 * p2**3 * y2p**4

! x*y-moment in rotated frame.

      tempab = 0.25 * x1p**2 * y2p**2 + ob3 * p12 * x1p * y2p**3
     1       + 0.125 * p12**2 * y2p**4 - 0.125 * p2**2 * y2p**4

! y^2-moment in rotated frame.

      tempac = ob3 * x1p * y2p**3 + 0.25 * p12 * y2p**4
     1       - 0.25 * p2 * y2p**4

! ------------------------------

! Second x^2-moment in non-rotated frame.

      moments_arr ( 4) = ct**2 * tempaa - 2.0 * ct * st * tempab
     1                 + st**2 * tempac

! Second x*y-moment in non-rotated frame.

      moments_arr ( 5) = ct * st * tempaa + ( ct**2 - st**2) * tempab
     1                 - ct * st * tempac

! Second y^2-moment in non-rotated frame.

      moments_arr ( 6) = st**2 * tempaa + 2.0 * ct * st * tempab
     1                 + ct**2 * tempac

! ----------------------------------------------------------------------
! STEP IV) Find the third moments.

! x^3-moment in rotated frame.

      tempaa = 0.25 * x1p**4 * y2p + 0.5 * p12 * x1p**3 * y2p**2
     1       + 0.5 * p12**2 * x1p**2 * y2p**3
     1       + 0.25 * p12**3 * x1p * y2p**4
     1       + ob20 * p12**4 * y2p**5 - ob20 * p2**4 * y2p**5

! x^2 * y -moment in rotated frame.

      tempab = ob6 * x1p**3 * y2p**2 + ob3 * p12 * x1p**2 * y2p**3
     1       + 0.25 * p12**2 * x1p * y2p**4 + ob15 * p12**3 * y2p**5
     1       - ob15 * p2**3 * y2p**5

! x * y^2 -moment in rotated frame.
      
      tempac = ob6 * x1p**2 * y2p**3 + 0.25 * p12 * x1p * y2p**4
     1       + 0.1 * p12**2 * y2p**5 - 0.1 * p2**2 * y2p**5

! y^3 -moment in rotated frame.

      tempad = 0.25 * x1p * y2p**4 + 0.2 * p12 * y2p**5
     1       - 0.2 * p2 * y2p**5

! ------------------------------

! Third x^3-moment in non-rotated frame.

      moments_arr ( 7) = ct**3 * tempaa - 3.0 * ct**2 * st * tempab
     1                 + 3.0 * ct * st**2 * tempac - st**3 * tempad

! Third x^2 * y -moment in non-rotated frame.

      moments_arr ( 8) = ct**2 * st * tempaa
     1                 + ( ct**3 - 2.0 * ct * st**2 ) * tempab
     1                 + ( st**3 - 2.0 * ct**2 * st ) * tempac
     1                 + ct * st**2 * tempad

! Third x * y^2 -moment in non-rotated frame.

      moments_arr ( 9) = ct * st**2 * tempaa
     1                 + ( 2 * ct**2 * st - st**3 ) * tempab
     1                 + ( ct**3 - 2 * ct * st**2 ) * tempac
     1                 - ct**2 * st * tempad

! Third y^3 -moment in non-rotated frame.

      moments_arr ( 10) = st**3 * tempaa + 3.0 * ct * st**2 * tempab
     1                  + 3.0 * ct**2 * st * tempac + ct**3 * tempad

! ----------------------------------------------------------------------

      END SUBROUTINE XY_MOMENTS_OF_TRIANGLES

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)
      PROGRAM TEST_MULTID_RS

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER i, j, nzone_vrtx, multid_wavemodel_num,

! Is the wave model subsonic?
     1        Is_SuperSonic, Is_SuperSonic_State, ifail

      REAL smallnum, prsfloor, rhofloor, pi, tempaa, expand_wavespeeds,

! Amount of HLL we will have to blend in for stabilization.
     1     blend_HLL, blend_MDHLL,

! Supersonic speeds in the global frame; always built even if flow isn't 
! supersonic.
     1     Vlx_Centroid, Vly_Centroid, S_Radial_Max,

! "yF_Lo, yF_Hi" are points of intersection where x = 0 intersects wave model.
     1     yF_Lo, yF_Hi,

! "xG_Lo, xG_Hi" are points of intersection where y = 0 intersects wave model.
     1     xG_Lo, xG_Hi,

! First 10 moments of x and y over areas contributing to "U", "F" and "G".
     1     Int_Momnt_U ( 10), Int_Momnt_F ( 10), Int_Momnt_G ( 10),


! Unit normal to the contact discontinuity; total and gas pressure in
! resolved state.
     1     Vlx_ALE, Vly_ALE,
     1     Prs_Tot_Star, Prs_Gas_Star


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables that come together at a zone edge.
! We also come in with their primitive variables and associated
! x- and y-fluxes.
     1      U, V, F, G,

! The two resolved states and their fluxes in the x- and y-directions.
     1      U_Minus, U_Star_Minus, U_Star_Plus, U_Plus,
     1      F_Minus, F_Star_Minus, F_Star_Plus, F_Plus,
     1      G_Minus, G_Star_Minus, G_Star_Plus, G_Plus,

! Averages of the two resolved states and two resolved y-fluxes.
     1      U_Star_Avg, G_Star_Avg,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y


      INTEGER, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) :: IF_Mid, IG_Mid


      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

     1      entropy, F_entropy_flux_x, G_entropy_flux_y,

! Components of the unit normal.
     1      nx, ny,

! The left, middle and right speeds of the HLLC RS in the Tau-direction,
! "S_Tau_Min, S_Tau_Mid, S_Tau_Max".
     1      S_Tau_Min, S_Tau_Mid, S_Tau_Max,

! The extremal speeds in the Neta-direction, "S_Neta_Min, S_Neta_Max".
     1      S_Neta_Min, S_Neta_Max,

! Fastest speed along the principal direction. Extents of the local
! wavefront. Intersection points of the wave fronts.
     1      S_Neta_Multid, L_Tau_Min, L_Tau_Max, x, y,

! Unit vector along the edge of a wave model. Dot products of that unit vector
! with "n" and "t" unit vectors.
     1     neta_x, neta_y, neta_dot_n, neta_dot_t,

! "LF_Tau_Mid, LG_Tau_Mid" are mid-points of the facial integration.
     1     LF_Tau_Mid, LG_Tau_Mid


      REAL, DIMENSION ( MAX_ZONES_PER_VRTX, NWAVES_1DRS) ::

! Wave speeds in the 1D Riemann problem; unparametrized and parametrized.
     1     S_1DRS, L_1DRS


! Final, strongly-interacting state and Global x- and y-fluxes.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: U_Star, F_Star, G_Star

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "test_multid_rs.hllem.out")

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      multid_wavemodel_num = 0  ! Choice of multidimensional wave model (0-3).

! Set up some constants.
      smallnum = 1.0e-12
      prsfloor = 1.0e-12
      rhofloor = 1.0e-8

      pi = 4.0 * ATAN ( 1.0)

      expand_wavespeeds = 1.05

      blend_HLL = 0.0
      blend_MDHLL = 0.0

      Vlx_ALE = 0.0 
      Vly_ALE = 0.0

! ------------------------------

! Build vector of primitives for all "nzone_vrtx" zones.
      DO i = 1, nzone_vrtx

      gamma_eos ( i) = 1.4
      gamma_sdsp ( i) = gamma_eos ( i)

      V ( i, 1) = gamma_eos ( i)

      V ( i, 2) = 0.2 + i * 0.1
      V ( i, 3) = 0.3 - i * 0.15
      V ( i, 4) = 0.3

      V ( i, 5) = 1.0

      V ( i, 6) = 1.0 + i * 0.3
      V ( i, 7) = -2.0 + i * 0.4
      V ( i, 8) = 1.5

      END DO

! Build vector of conserved variables for all "nzone_vrtx" zones.
      DO i = 1, nzone_vrtx

      entropy ( i) = V ( i, 5) / V ( i, 1)**( gamma_eos ( i) - 1.0)

      U ( i, 1) = V ( i, 1)

      U ( i, 2) = V ( i, 1) * V ( i, 2)
      U ( i, 3) = V ( i, 1) * V ( i, 3)
      U ( i, 4) = V ( i, 1) * V ( i, 4)

      U ( i, 5) = V ( i, 5) / ( gamma_eos ( i) - 1.0)
     1          + 0.5 * V ( i, 1) * ( V ( i, 2)**2 + V ( i, 3)**2
     1                              + V ( i, 4)**2 )
     1          + ( V ( i, 6)**2 + V ( i, 7)**2 + V ( i, 8)**2)
     1          / ( 8.0 * pi)

      U ( i, 6) = V ( i, 6)
      U ( i, 7) = V ( i, 7)
      U ( i, 8) = V ( i, 8)

      END DO

! ------------------------------

! Make unit normals along the principal directions. Default is cartesian
! mesh test with four unit normals.

#if ( 0 == 1)
      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0
#endif

      DO i = 1, nzone_vrtx

      tempaa = 2.0 * pi * ( i - 1.0) / REAL ( nzone_vrtx) + 0.5 * pi

      nx ( i) = COS ( tempaa)
      ny ( i) = SIN ( tempaa)

      END DO

! ----------------------------------------------------------------------
! Call "ONE_D_MHD_PRIM_F_G_WAVSPD" and write out the wave speeds.

      CALL ONE_D_MHD_PRIM_F_G_WAVSPD ( nzone_vrtx,

     1     smallnum, prsfloor, rhofloor, pi, expand_wavespeeds,

     1     U, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     V, F, G,
     1     S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max)

! ------------------------------

      DO i = 1, nzone_vrtx

! Make an approximate fake foliation of the intermediate waves.

      DO j = 1, NWAVES_1DRS

      S_1DRS ( i, j) = S_Tau_Min ( i)
     1               + j * ( S_Tau_Max ( i) - S_Tau_Min ( i))
     1               / ( NWAVES_1DRS + 1.0)

      END DO

      END DO  ! End of "DO i =" loop.


      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"After ONE_D_MHD_PRIM_F_G_WAVSPD"
      WRITE ( 6, *)" "


      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)"i = ", i,
     1             " S_Tau_Min, S_Tau_Max, S_Neta_Min, S_Neta_Max = "
      WRITE ( 6, *) S_Tau_Min ( i), S_Tau_Max ( i), S_Neta_Min ( i),
     1              S_Neta_Max ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"V = "
      WRITE ( 6, *) ( V ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"U = "
      WRITE ( 6, *) ( U ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"F = "
      WRITE ( 6, *) ( F ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"G = "
      WRITE ( 6, *) ( G ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "

      END DO

! ----------------------------------------------------------------------
! Call either "MULTID_WAVE_MODEL" or "MULTID_WAVE_MODEL_MINIMALIST"
! and write out relevant variables.

      CALL MULTID_WAVE_MODEL
!      CALL MULTID_WAVE_MODEL_MINIMALIST

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, blend_HLL, expand_wavespeeds, Vlx_ALE, Vly_ALE,

     1     U, gamma_eos, gamma_sdsp, S_1DRS,

     1     S_Neta_Min, S_Neta_Max, S_Tau_Min, S_Tau_Max, nx, ny,

     1     Vlx_Centroid, Vly_Centroid, S_Radial_Max,
     1     Is_SuperSonic, Is_SuperSonic_State,

     1     S_Neta_Multid, L_Tau_Min, L_Tau_Max, x, y,
     1     neta_x, neta_y, neta_dot_n, neta_dot_t, L_1DRS,

     1     IF_Mid, LF_Tau_Mid, yF_Lo, yF_Hi,
     1     IG_Mid, LG_Tau_Mid, xG_Lo, xG_Hi,
     1     Int_Momnt_U, Int_Momnt_F, Int_Momnt_G, ifail)

! ------------------------------

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"After MULTID_WAVE_MODEL"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ifail = ", ifail
      WRITE ( 6, *)" "
      WRITE ( 6, *)"First the scalars:"
      WRITE ( 6, *)" "

      WRITE ( 6, *)"Is_SuperSonic = ", Is_SuperSonic
      WRITE ( 6, *)"Is_SuperSonic_State = ", Is_SuperSonic_State
      WRITE ( 6, *)"Vlx_ALE, Vly_ALE = ", Vlx_ALE, Vly_ALE
      WRITE ( 6, *)" "

      WRITE ( 6, *)"Vlx_Centroid, Vly_Centroid, S_Radial_Max = ",
     1              Vlx_Centroid, Vly_Centroid, S_Radial_Max
      WRITE ( 6, *)" "

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yF_Lo, yF_Hi = ", yF_Lo, yF_Hi
      WRITE ( 6, *)"xG_Lo, xG_Hi = ", xG_Lo, xG_Hi
      WRITE ( 6, *)" "

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Int_Momnt_U ( :) = ", Int_Momnt_U ( :)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Int_Momnt_F ( :) = ", Int_Momnt_F ( :)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Int_Momnt_G ( :) = ", Int_Momnt_G ( :)
      WRITE ( 6, *)" "


      WRITE ( 6, *)"Next the vectors:"
      WRITE ( 6, *)" "

      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"i = ", i," S_Neta_Multid = ", S_Neta_Multid ( i)
      WRITE ( 6, *)" "
      DO j = 1, nzone_vrtx
      WRITE ( 6, *)"j = ", j," vertex = ( ", x ( j)," , ", y ( j)," ) "
      END DO
      WRITE ( 6, *)" "

      WRITE ( 6, *)"L_Tau_Min, L_Tau_Max = ",
     1              L_Tau_Min ( i), L_Tau_Max ( i)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"nx, ny, x, y = ", nx ( i), ny ( i), x ( i), y ( i)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"neta_x, neta_y = ", neta_x ( i), neta_y ( i)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"neta_dot_n, neta_dot_t = ",
     1              neta_dot_n ( i), neta_dot_t ( i)
      WRITE ( 6, *)" "
      
      WRITE ( 6, *)"----------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"S_1DRS ( i, :) = ", S_1DRS ( i, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"L_1DRS ( i, :) = ", L_1DRS ( i, :)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"IF_Mid, IG_Mid = ", IF_Mid ( i), IG_Mid ( i)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"LF_Tau_Mid, LG_Tau_Mid = ",
     1              LF_Tau_Mid ( i), LG_Tau_Mid ( i)
      WRITE ( 6, *)" "

      END DO

! ----------------------------------------------------------------------
! Call "MuSIC_MHD_2DRS_1_STRUCT".

      CALL MuSIC_MHD_2DRS_1_STRUCT

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, prsfloor, rhofloor, pi,

     1     blend_HLL, blend_MDHLL,
     1     Prs_Tot_Star, Prs_Gas_Star,

     1     Vlx_ALE, Vly_ALE,

     1     U, F, G, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     U_Star, F_Star, G_Star,
     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)

! ------------------------------

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"After MuSIC_MHD_2DRS_1_STRUCT : "
      WRITE ( 6, *)"----------"
      WRITE ( 6, *)"U_Star = "
      WRITE ( 6, *) ( U_Star ( j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"F_Star = "
      WRITE ( 6, *) ( F_Star ( j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"G_Star = "
      WRITE ( 6, *) ( G_Star ( j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Prs_Tot_Star, Prs_Gas_Star = ",
     1              Prs_Tot_Star, Prs_Gas_Star
      WRITE ( 6, *)" "

      DO i = 1, nzone_vrtx

      WRITE ( 6, *)"----------"
      WRITE ( 6, *)"U_HLLC = "
      WRITE ( 6, *) ( U_HLLC ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"F_HLLC = "
      WRITE ( 6, *) ( F_HLLC ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"G_HLLC = "
      WRITE ( 6, *) ( G_HLLC ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"F_HLLC_Global_x = "
      WRITE ( 6, *)
     1  ( F_HLLC_Global_x ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"G_HLLC_Global_y = "
      WRITE ( 6, *)
     1  ( G_HLLC_Global_y ( i, j), j = 1, NFIELD + 1 + NFLUID)
      WRITE ( 6, *)" "

      END DO

! ----------------------------------------------------------------------

      END PROGRAM TEST_MULTID_RS
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)

      PROGRAM PLOT_TWO_D_RP

! For this code, we don't want a dissipative Riemann solver. Instead,
! we want the solution to the true multidimensional Riemann problem. The
! settings are, therefore, a little different from the ones needed for
! obtaining numerical fluxes. They are:
! "NWAVES_1DRS = 2, 3" .
! "EXTREMIZE_WAVESPEEDS = 0".
! "multid_wavemodel_num = 0, 1", "expand_wavespeeds = 1.05, 1.1".
! Also uncomment the lines which use "Vlx_ALE, Vly_ALE" in the multid RS.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nzone_vrtx, multid_wavemodel_num

      REAL smallnum, prsfloor, rhofloor, pi, gamma,

! Amount of HLL we will have to blend in for stabilization.
     1     blend_HLL, blend_MDHLL,

! Velocity of vertex; total and gas pressure in resolved state.
     1     Vlx_ALE, Vly_ALE,
     1     Prs_Tot_Star, Prs_Gas_Star

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX,
     1                  NFIELD + 1 + NFLUID) ::

! Cyclically arranged state variables that come together at a vertex.
     1      U, V, F, G,

! Resolved numerical HLLC state and flux in the Tau- and Neta-directions.
! Also provide the same fluxes in the global frame of reference.
     1      U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y

      REAL, DIMENSION ( 1 : MAX_ZONES_PER_VRTX) ::

! Polytropic indices on either side.
     1      gamma_eos, gamma_sdsp,

! entropy variables that go through RS interface.
     1      entropy, F_entropy_flux_x, G_entropy_flux_y,

! Components of the unit normal.
     1      nx, ny

! Final, strongly-interacting state and Global x- and y-fluxes.
! These are the states that straddle the time axis.
      REAL, DIMENSION ( NFIELD + 1 + NFLUID) :: U_Star, F_Star, G_Star

! ----------------------------------------------------------------------

      INTEGER ix1, iy1, i, j, k, iproblem 

      PARAMETER ( ix1 = 400, iy1 = 400)

      REAL xmin, xmax, ymin, ymax, dx, dy, expand_wavespeeds

      REAL, DIMENSION ( ix1, iy1, NFIELD + 1 + NFLUID) :: U_Plot

      REAL*4 r4gr ( ix1, iy1)

      REAL, DIMENSION ( ix1) :: xcsave

      REAL, DIMENSION ( iy1) :: ycsave

! ----------------------------------------------------------------------

! Define the domain and zone-centered coordinates.

      xmin = -2.5
      xmax = 2.5
      ymin = -2.5
      ymax = 2.5

      dx = ( xmax - xmin) / ix1
      dy = ( ymax - ymin) / iy1

      DO i = 1, ix1
      xcsave ( i) = xmin + ( i - 0.5) * dx
      END DO

      DO j = 1, iy1
      ycsave ( j) = ymin + ( j - 0.5) * dy
      END DO

! ----------------------------------------------------------------------

      nzone_vrtx = 4  ! Number of zones coming together at the vertex.

      multid_wavemodel_num = 0  ! Choice of multidimensional wave model (0-3).

! Set up some constants.
      smallnum = 1.0e-12
      prsfloor = 1.0e-12
      rhofloor = 1.0e-8

      pi = 4.0 * ATAN ( 1.0)

      expand_wavespeeds = 1.05

      blend_HLL = 0.0
      blend_MDHLL = 0.0

      gamma = 1.4
      gamma_eos ( :) = gamma
      gamma_sdsp ( :) = gamma_eos ( :)

! ------------------------------

! Make unit normals along the principal directions. Default is cartesian
! mesh test with four unit normals.

      nx ( 1) = 0.0
      ny ( 1) = 1.0

      nx ( 2) = -1.0
      ny ( 2) = 0.0

      nx ( 3) = 0.0
      ny ( 3) = -1.0

      nx ( 4) = 1.0
      ny ( 4) = 0.0

! ------------------------------

      iproblem = 1

      V ( :, :) = 0.0

      IF ( iproblem == 1) THEN

! First 2d Riemann problem.
! S_R = 1.91; S_L = -1.32; S_U = 1.91; S_D = -1.32 for this problem.
! Use that information to set the domain boundaries.
! "[-2.5, 2.5]X[-2.5, 2.5]" is a good set of values for
! "[xmin, xmax]X[ymin, ymax]".
! When plotting the density data, use "datamin = 0.5, datamax = 1.8".

      V ( 1, 1) = 0.5313
      V ( 1, 5) = 0.4
      V ( 1, 2) = 0.0
      V ( 1, 3) = 0.0

      V ( 2, 1) = 1.0
      V ( 2, 5) = 1.0
      V ( 2, 2) = 0.7276
      V ( 2, 3) = 0.0

      V ( 3, 1) = 0.8
      V ( 3, 5) = 1.0
      V ( 3, 2) = 0.0
      V ( 3, 3) = 0.0

      V ( 4, 1) = 1.0
      V ( 4, 5) = 1.0
      V ( 4, 2) = 0.0
      V ( 4, 3) = 0.7276

      ELSE IF ( iproblem == 2) THEN

! Second 2d Riemann problem.
! S_R = 2.09; S_L = -1.18; S_U = 2.09; S_D = -1.18 for this problem.
! Use that information to set the domain boundaries.
! "[-2.0, 2.0]X[-2.0, 2.0]" is a good set of values for
! "[xmin, xmax]X[ymin, ymax]".
! When plotting the density data, use "datamin = 0.5, datamax = 1.8".

      V ( 1, 1) = 1.5
      V ( 1, 5) = 1.5
      V ( 1, 2) = 0.0
      V ( 1, 3) = 0.0

      V ( 2, 1) = 0.5323
      V ( 2, 5) = 0.3
      V ( 2, 2) = 1.206
      V ( 2, 3) = 0.0

      V ( 3, 1) = 0.1379
      V ( 3, 5) = 0.029
      V ( 3, 2) = 1.206
      V ( 3, 3) = 1.206

      V ( 4, 1) = 0.5323
      V ( 4, 5) = 0.3
      V ( 4, 2) = 0.0
      V ( 4, 3) = 1.206

      ELSE

! "[-2.5, 2.5]X[-2.5, 2.5]" is a good set of values for
! "[xmin, xmax]X[ymin, ymax]".

      V ( 1, 1) = 1.1
      V ( 1, 5) = 1.1
      V ( 1, 2) = 0.0
      V ( 1, 3) = 0.0

      V ( 2, 1) = 0.506474820143885
      V ( 2, 5) = 0.35
      V ( 2, 2) = 0.893871177417743
      V ( 2, 3) = 0.0

      V ( 3, 1) = 1.1
      V ( 3, 5) = 1.1
      V ( 3, 2) = 0.893871177417743
      V ( 3, 3) = 0.893871177417743

      V ( 4, 1) = 0.506474820143885
      V ( 4, 5) = 0.35
      V ( 4, 2) = 0.0
      V ( 4, 3) = 0.893871177417743

      END IF

! ------------------------------

! Build vector of conserved variables for all "nzone_vrtx" zones.
      DO i = 1, nzone_vrtx

      entropy ( i) = V ( i, 5) / V ( i, 1)**( gamma_eos ( i) - 1.0)

      U ( i, 1) = V ( i, 1)

      U ( i, 2) = V ( i, 1) * V ( i, 2)
      U ( i, 3) = V ( i, 1) * V ( i, 3)
      U ( i, 4) = V ( i, 1) * V ( i, 4)

      U ( i, 5) = V ( i, 5) / ( gamma_eos ( i) - 1.0)
     1          + 0.5 * V ( i, 1) * ( V ( i, 2)**2 + V ( i, 3)**2
     1                              + V ( i, 4)**2 )
     1          + ( V ( i, 6)**2 + V ( i, 7)**2 + V ( i, 8)**2)
     1          / ( 8.0 * pi)

      U ( i, 6) = V ( i, 6)
      U ( i, 7) = V ( i, 7)
      U ( i, 8) = V ( i, 8)

      END DO

! ----------------------------------------------------------------------

      DO j = 1, iy1

        DO i = 1, ix1

        Vlx_ALE = xcsave ( i)
        Vly_ALE = ycsave ( j)

! Call "MuSIC_MHD_2DRS_1_STRUCT".

      CALL MuSIC_MHD_2DRS_1_STRUCT

     1   ( nzone_vrtx, multid_wavemodel_num,

     1     smallnum, prsfloor, rhofloor, pi,

     1     blend_HLL, blend_MDHLL,
     1     Prs_Tot_Star, Prs_Gas_Star,

     1     Vlx_ALE, Vly_ALE,

     1     U, F, G, gamma_eos, gamma_sdsp, nx, ny,
     1     entropy,

     1     U_Star, F_Star, G_Star,
     1     U_HLLC, F_HLLC, G_HLLC, F_HLLC_Global_x, G_HLLC_Global_y,
     1     F_entropy_flux_x, G_entropy_flux_y)


        U_Plot ( i, j, :) = U_Star ( :)


        END DO

      END DO 

! ----------------------------------------------------------------------
! Now write out the density variables that have been stored in "U_Plot".

      OPEN ( UNIT = 99, FILE = "rh030001")

      r4gr ( :, :) = U_Plot ( :, :, 1)

      WRITE ( 99, *) ( ( r4gr ( i, j), i = 1, ix1), j = 1, iy1)

      CLOSE ( 99)

! ----------

      OPEN ( UNIT = 99, FILE = "vx030001")

      r4gr ( :, :) = U_Plot ( :, :, 2) / U_Plot ( :, :, 1)

      WRITE ( 99, *) ( ( r4gr ( i, j), i = 1, ix1), j = 1, iy1)

      CLOSE ( 99)

! ----------

      OPEN ( UNIT = 99, FILE = "vy030001")

      r4gr ( :, :) = U_Plot ( :, :, 3) / U_Plot ( :, :, 1)

      WRITE ( 99, *) ( ( r4gr ( i, j), i = 1, ix1), j = 1, iy1)

      CLOSE ( 99)

! ----------

      OPEN ( UNIT = 99, FILE = "pr030001")

      r4gr ( :, :) = U_Plot ( :, :, 5) - 0.5 * ( U_Plot ( :, :, 2)**2
     1             + U_Plot ( :, :, 3)**2 ) / U_Plot ( :, :, 1)
      r4gr ( :, :) = ( gamma - 1.0) * r4gr ( :, :)

      WRITE ( 99, *) ( ( r4gr ( i, j), i = 1, ix1), j = 1, iy1)

      CLOSE ( 99)

! ----------------------------------------------------------------------

      END PROGRAM PLOT_TWO_D_RP
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( 0 == 1)
      PROGRAM TEST_XY_MOMENTS_OF_TRIANGLES

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER i, ip1, npolygon

      PARAMETER ( npolygon = 5)

      REAL x1, y1, x2, y2, xcen, ycen, moments_arr ( 10),
     1     pi, theta, rpoly, xpoly ( npolygon), ypoly ( npolygon),
     1     xshift, yshift, a, b,
     1     add_moments ( 10)

! ----------------------------------------------------------------------

      pi = 4.0 * ATAN ( 1.0)

      theta = 2.0 * pi / npolygon

      rpoly = 3.0

! Shift must be small enough to ensure that origin stays inside the polygon.
      xshift = 0.2 * rpoly
      yshift = -0.3 * rpoly

! ----------------------------------------------------------------------
! First and simplest test. Hand in any triangle with vertices O, P1, P2
! counterclockwise. The first moments divided by the zeroth moment should
! return the centroid.

! For certain right triangles, the second moments are easy to evaluate.
! For example, with "x1 = 1.0, y1 = 0.0, x2 = 1.0, y2 = 1.0" we get
! x^2 moment 0.25; x*y moment to be 0.125 and y^2 moment to be 1/12.
! Rotate/reflect/expand/contract this triangle to test other variations.

! The third moments are also easy to evaluate for the same triangles.
! For example, with "x1 = 1.0, y1 = 0.0, x2 = 1.0, y2 = 1.0" we get
! x^3 moment 0.2 ; x**2 * y moment to be 0.1 ; x * y**2 moment to be 1/15 ;
! and y^3 moment to be 0.05 .

      x1 = 1.0
      y1 = 0.5

      x2 = 0.7
      y2 = 2.0

      CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,
     1           xcen, ycen, moments_arr )

      WRITE ( 6, *)" "
      WRITE ( 6, *)"First Test: "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xcen, ycen = ", xcen, ycen
      WRITE ( 6, *)"Match with : ", moments_arr ( 2) / moments_arr ( 1),
     1              moments_arr ( 3) / moments_arr ( 1)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Second moments 4, 5, 6 = ", moments_arr ( 4),
     1              moments_arr ( 5), moments_arr ( 6)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Third moments 7, 8, 9, 10 = ", moments_arr ( 7),
     1              moments_arr ( 8), moments_arr ( 9),
     1              moments_arr ( 10)

! ----------------------------------------------------------------------
! Second test. Set up vertices of a regular polygon centered on the origin.
! The center of mass must be zero.

! For a square/rectangle, the second moments are easy to evaluate.
! Uncomment the lines below to see how. For the rectangle [-a,a]x[-b,b]
! the x^2 moment is "4.0 * a**3 * b / 3.0". The x*y moment is zero.
! The y^2 moment is "4.0 * a * b**3 / 3.0".

      DO i = 1, npolygon
      xpoly ( i) = rpoly * COS ( i * theta + 0.1796 * pi)
      ypoly ( i) = rpoly * SIN ( i * theta + 0.1796 * pi)
      END DO

! Vertices of a rectangle. Must set npolygon = 4.
#if ( 0 == 1)
      a = 2.0
      b = 1.0
      xpoly ( 1) = a
      ypoly ( 1) = b
      xpoly ( 2) = -a
      ypoly ( 2) = b
      xpoly ( 3) = -a
      ypoly ( 3) = -b
      xpoly ( 4) = a
      ypoly ( 4) = -b
#endif

      add_moments ( :) = 0.0

      DO i = 1, npolygon

      ip1 = i + 1
      IF ( ip1 > npolygon) ip1 = 1

      x1 = xpoly ( i)
      y1 = ypoly ( i)
      x2 = xpoly ( ip1)
      y2 = ypoly ( ip1)

      CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,
     1           xcen, ycen, moments_arr )

      add_moments ( :) = add_moments ( :) + moments_arr ( :)

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Second Test: Center of mass of regular polygon"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"must be origin for this test : ",
     1             add_moments ( 2) / add_moments ( 1),
     1             add_moments ( 3) / add_moments ( 1)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Second moments, x * y moment must be zero"
      WRITE ( 6, *)"add_moments 4, 5, 6 = ", add_moments ( 4),
     1              add_moments ( 5), add_moments ( 6)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Third moments 7, 8, 9, 10 = ", moments_arr ( 7),
     1              moments_arr ( 8), moments_arr ( 9), 
     1              moments_arr ( 10)

! ----------------------------------------------------------------------
! Third test. Set up vertices of a shifted regular polygon.
! The shift must match up.

      DO i = 1, npolygon
      xpoly ( i) = xshift + rpoly * COS ( i * theta + 0.1796 * pi)
      ypoly ( i) = yshift + rpoly * SIN ( i * theta + 0.1796 * pi)
      END DO


      add_moments ( :) = 0.0

      DO i = 1, npolygon

      ip1 = i + 1
      IF ( ip1 > npolygon) ip1 = 1

      x1 = xpoly ( i)
      y1 = ypoly ( i)
      x2 = xpoly ( ip1)
      y2 = ypoly ( ip1)

      CALL XY_MOMENTS_OF_TRIANGLES ( x1, y1, x2, y2,
     1           xcen, ycen, moments_arr )

      add_moments ( :) = add_moments ( :) + moments_arr ( :)

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Third Test : Shifted regular polygon"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Original center of mass : ", xshift, yshift
      WRITE ( 6, *)"Must match with : ",
     1             add_moments ( 2) / add_moments ( 1),
     1             add_moments ( 3) / add_moments ( 1)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Second moments 4, 5, 6 = ", add_moments ( 4),
     1              add_moments ( 5), add_moments ( 6)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Third moments 7, 8, 9, 10 = ", moments_arr ( 7),
     1              moments_arr ( 8), moments_arr ( 9), 
     1              moments_arr ( 10)

! ----------------------------------------------------------------------

      END PROGRAM TEST_XY_MOMENTS_OF_TRIANGLES
#endif

! ----------------------------------------------------------------------






























