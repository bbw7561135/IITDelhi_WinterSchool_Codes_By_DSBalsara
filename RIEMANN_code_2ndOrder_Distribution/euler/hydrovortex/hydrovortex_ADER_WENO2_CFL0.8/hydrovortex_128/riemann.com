! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      INTEGER igeom, ientropyfix

      REAL pi, gamma, molewt, gasconst, rhofloor, prsfloor, smallnum,
     1     gamma_euler, molewt_euler, gasconst_euler,
     1     isospeed, barotropic_stiffness, lightspd, stefan,

     1     cfl_coef_sdsp, cfl_coef_advct,
     1     mc_coef, flatten_coef, prs_flatten_coef,
     1     fracdiff, divratio_einfeldt_rs,
     1     lapidus_threshold, lapidus_coef,
     1     blend_multid_flux, blend_mdhll_to_mdhllc,

     1     base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1     base_grid_ymax, base_grid_zmin, base_grid_zmax,

     1     gamma_specie ( 0 : NFLUID), cp_specie ( 0 : NFLUID),
     1     cv_specie ( 0 : NFLUID), molewt_specie ( 0 : NFLUID),
     1     heat_of_formation ( 0 : NFLUID),

     1     gamma_specie_euler ( 0 : NFLUID_EULER),
     1     cp_specie_euler ( 0 : NFLUID_EULER),
     1     cv_specie_euler ( 0 : NFLUID_EULER),
     1     molewt_specie_euler ( 0 : NFLUID_EULER),
     1     heat_of_formation_euler ( 0 : NFLUID_EULER)

      COMMON/RIEMANN_COM/ igeom, ientropyfix,
     1     pi, gamma, molewt, gasconst, rhofloor, prsfloor, smallnum,
     1     gamma_euler, molewt_euler, gasconst_euler,
     1     isospeed, barotropic_stiffness, lightspd, stefan,

     1     cfl_coef_sdsp, cfl_coef_advct,
     1     mc_coef, flatten_coef, prs_flatten_coef,
     1     fracdiff, divratio_einfeldt_rs,
     1     lapidus_threshold, lapidus_coef,
     1     blend_multid_flux, blend_mdhll_to_mdhllc,

     1     base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1     base_grid_ymax, base_grid_zmin, base_grid_zmax,

     1     gamma_specie, cp_specie,
     1     cv_specie, molewt_specie,
     1     heat_of_formation,

     1     gamma_specie_euler, cp_specie_euler,
     1     cv_specie_euler, molewt_specie_euler,
     1     heat_of_formation_euler

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------







