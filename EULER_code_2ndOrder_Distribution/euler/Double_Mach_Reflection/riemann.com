! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      INTEGER igeom, ientropyfix

      REAL pi, gamma_euler, molewt_euler, gasconst_euler,
     1     rhofloor, prsfloor, smallnum,
     1     isospeed, barotropic_stiffness, lightspd, stefan,

     1     cfl_coef, mc_coef, flatten_coef, prs_flatten_coef,
     1     fracdiff, divratio_einfeldt_rs,
     1     lapidus_threshold, lapidus_coef,

     1     base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1     base_grid_ymax, base_grid_zmin, base_grid_zmax,

     1     gamma_specie_euler ( 0 : NFLUID_EULER),
     1     cp_specie_euler ( 0 : NFLUID_EULER),
     1     cv_specie_euler ( 0 : NFLUID_EULER),
     1     molewt_specie_euler ( 0 : NFLUID_EULER),
     1     heat_of_formation_euler ( 0 : NFLUID_EULER)

      COMMON/RIEMANN_COM/ igeom, ientropyfix,
     1     pi, gamma_euler, molewt_euler, gasconst_euler,
     1     rhofloor, prsfloor, smallnum,
     1     isospeed, barotropic_stiffness, lightspd, stefan,

     1     cfl_coef, mc_coef, flatten_coef, prs_flatten_coef,
     1     fracdiff, divratio_einfeldt_rs,
     1     lapidus_threshold, lapidus_coef,

     1     base_grid_xmin, base_grid_xmax, base_grid_ymin,
     1     base_grid_ymax, base_grid_zmin, base_grid_zmax,

     1     gamma_specie_euler, cp_specie_euler,
     1     cv_specie_euler, molewt_specie_euler,
     1     heat_of_formation_euler

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------







