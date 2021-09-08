MODULE cable_implicit_main_mod
  
CONTAINS

SUBROUTINE cable_implicit_main( cycleno, & ! num_cycles
              row_length, rows, land_pts, ntiles, sm_levels,                  &
              Fland,                                                          &
              !forcing & increments: temp and humidity 
              tl_1, qw_1, dtl1_1, dqw1_1, ctctq1,                             &
              !prog: UM soi , per land_pt
              ftl_1, ftl_surft, fqw_1, fqw_surft,                             &
              tstar_surft, surf_ht_flux_land, surf_htf_surft,                 &
              ecan_surft, esoil_surft, ei_surft, radnet_surft,                &
              t1p5m_surft, q1p5m_surft, melt_surft,                           &
              dtstar_surft,                                                   &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soil_cbl )
  !subrs called 
USE cable_implicit_driv_mod, ONLY: cable_implicit_driver
USE cable_implicit_unpack_mod, ONLY: implicit_unpack

!cable progs are set here
USE cable_prognostic_info_mod, ONLY:                                          &
  SoilTemp        =>  SoilTemp_CABLE,                                         &
  SoilMoisture    =>  SoilMoisture_CABLE,                                     &
  FrozenSoilFrac  => FrozenSoilFrac_CABLE,                                    &
  SnowDepth   => SnowDepth_CABLE,                                             &
  SnowMass    => SnowMass_CABLE,                                              &
  SnowTemp    => SnowTemp_CABLE,                                              &
  SnowDensity => SnowDensity_CABLE,                                           &
  SnowAge     => SnowAge_CABLE,                                               &
  ThreeLayerSnowFlag  => ThreeLayerSnowFlag_CABLE,                            &
  OneLyrSnowDensity   => OneLyrSnowDensity_CABLE


!data
USE cable_types_mod, ONLY: mp
USE cable_types_mod, ONLY: L_tile_pts
 
  ! processor number, timestep number & width
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, cable_runtime
  
USE jules_surface_types_mod, ONLY: npft
 
!cable progs are set here
USE cable_prognostic_info_mod, ONLY:                                          &
  soil_temp_cable   =>  SoilTemp_CABLE,                                       &
  soil_moist_cable  =>  SoilMoisture_CABLE,                                   &
  soil_froz_frac_cable => FrozenSoilFrac_CABLE,                               &
  snow_dpth_cable => SnowDepth_CABLE,                                         &
  snow_mass_cable => SnowMass_CABLE,                                          &
  snow_temp_cable => SnowTemp_CABLE,                                          &
  snow_rho_cable  => SnowDensity_CABLE,                                       &
  snow_flg_cable  => ThreeLayerSnowFlag_CABLE,                                &
  snow_age_cable  => SnowAge_CABLE,                                           &
  snow_avg_rho_cable => OneLyrSnowDensity_CABLE

USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soil_parameter_type
USE cable_params_mod,         ONLY: soilin
USE ancil_info, ONLY: dim_cs1, dim_cs2

USE p_s_parms, ONLY: sthu_soilt, sthf_soilt
USE prognostics, ONLY: t_soil_soilt, smcl_soilt, snowdepth_surft 
USE prognostics, ONLY: snow_depth => snowdepth_surft 
USE prognostics, ONLY: canopy_gb 
USE prognostics, ONLY: canopy_tile => canopy_surft 
!Imports for driving and flux variables
USE forcing, ONLY: ls_rain_cable => ls_rain_ij,                               &
                    ls_snow_cable => ls_snow_ij,                              &
                    conv_rain_cable => con_rain_ij,                           &
                    conv_snow_cable => con_snow_ij


USE p_s_parms,      ONLY: smvcst_cable => smvcst_soilt
USE cable_other_constants_mod, ONLY: nrb
IMPLICIT NONE
  
!___ re-decl input args
INTEGER :: cycleno
INTEGER :: row_length,rows, land_pts, ntiles, sm_levels
  
REAL,  DIMENSION(land_pts) ::                                                 &
    fland
   
REAL, DIMENSION(row_length,rows) ::                                           &
 tl_1, qw_1,       &  !forcing: temp and humidity 
 dtl1_1, dqw1_1       !forcing: increments to temp and humidity 

!Ticket #132 needs ctctq1
REAL, DIMENSION(row_length,rows) ::                                           &
 ctctq1               !forcing: information for temp and humidity increment
  
!prog: UM soil quantities, non-tiled  aggregate over CABLE tiles per land_pt
REAL, DIMENSION(land_pts,sm_levels) ::                                        &
  t_soil,                                                                     &
  smcl,                                                                       &
  sthf,                                                                       &
  sthu
INTEGER :: i_day_number 
REAL, DIMENSION( land_pts, ntiles ) ::                                        &
  snow_surft, &     ! snow ammount on tile nee:snow_tile 
  ftl_surft,                                                                  &
  fqw_surft,                                                                  &
  tstar_surft,                                                                &
  surf_htf_surft,                                                             &
  ecan_surft, esoil_surft,                                                    &
  ei_surft, radnet_surft,                                                     &
  gs_surft,                                                                   &
  t1p5m_surft,                                                                &
  q1p5m_surft,                                                                &
  melt_surft,                                                                 &
  dtstar_surft
    
REAL, DIMENSION( land_pts ) ::                                                &
  gs               ! conductance for use in dust scheme
REAL :: dtrad(mp)
    
REAL, DIMENSION(row_length,rows) ::                                           &
  ftl_1,                                                                      &
  fqw_1,                                                                      &
  surf_ht_flux_land

REAL, DIMENSION( land_pts ) ::                                                &
  gpp,   & ! IN Gross primary productivity (kg C/m2/s).
  npp,   & ! IN Net primary productivity   (kg C/m2/s).
  resp_p   ! IN Plant respiration (kg C/m2/s).

REAL, DIMENSION( land_pts,ntiles ) ::                                         &
  gpp_ft,   & ! IN Gross primary productivity  on PFTs (kg C/m2/s).
  npp_ft,   & ! IN Net primary productivity  on PFTs (kg C/m2/s).
  g_leaf,   & ! IN Leaf turnover rate (/360days).
  resp_p_ft   ! IN Plant respiration on PFTs  (kg C/m2/s).
   
REAL :: resp_s(land_pts,dim_cs1) ! IN Soil respiration (kg C/m2/s).
REAL :: resp_s_tot(dim_cs2)      ! IN Total soil resp'n (kg C/m2/s).
 
!___ local vars

! UM type vars but no feedback 
REAL, DIMENSION(land_pts,ntiles,sm_levels) :: STHU_surft
REAL, DIMENSION(land_pts,ntiles) ::  tot_alb, TRANSP_surft
REAL, DIMENSION(land_pts,ntiles,3) :: snow_cond        

LOGICAL, SAVE :: first_call = .TRUE.
INTEGER,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable

TYPE(air_type),       INTENT(INOUT)  :: air_cbl
TYPE(met_type),       INTENT(INOUT)  :: met_cbl
TYPE(radiation_type),       INTENT(INOUT)  :: rad_cbl
TYPE(roughness_type),     INTENT(INOUT)  :: rough_cbl
TYPE(canopy_type),    INTENT(INOUT)  :: canopy_cbl
TYPE(soil_snow_type),     INTENT(INOUT)  :: ssnow_cbl
TYPE(bgc_pool_type),       INTENT(INOUT)  :: bgc_cbl
TYPE(balances_type),       INTENT(INOUT)  :: bal_cbl
TYPE(sum_flux_type),  INTENT(INOUT)  :: sum_flux_cbl
TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(soil_parameter_type),  INTENT(INOUT) ::  soil_cbl
  ! std template args  
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_implicit_main"
REAL :: dummy(1)

t_soil(:,:) = t_soil_soilt(:,1,:) 
smcl(:,:) = smcl_soilt(:,1,:) 
sthu(:,:) = sthu_soilt(:,1,:)  
sthf(:,:) = sthf_soilt(:,1,:) 
  
! FLAGS def. specific call to CABLE from UM
cable_runtime%um          = .FALSE.
cable_runtime%um_implicit = .TRUE.
cable_runtime%um_explicit = .FALSE.
  
isnow_flg_cable = INT(snow_flg_cable)

CALL cable_implicit_driver( i_day_number, cycleno, &! num_cycles  
      row_length,rows, land_pts, ntiles, npft, sm_levels,                     &
      dim_cs1, dim_cs2, Fland,                                                &
      LS_RAIN_cable, CONV_RAIN_cable, LS_SNOW_cable, CONV_SNOW_cable,         &
      dtl1_1, dqw1_1, ctctq1, t_soil, soil_temp_cable,                        &
      smcl, soil_moist_cable,                                                 &
      !aquifer_moist_cable, 
      REAL(kwidth_gl),                                                        &
      SMVCST_cable, sthf, soil_froz_frac_cable, sthu,                         &
      snow_surft, snow_avg_rho_cable, isnow_flg_cable, snow_dpth_cable,       &
      snow_mass_cable, snow_rho_cable, snow_temp_cable,                       &
      ftl_1, FTL_surft, fqw_1, FQW_surft,                                     &
      TSTAR_surft, surf_ht_flux_land, ECAN_surft, ESOIL_surft,                &
      EI_surft, RADNET_surft, SNOW_AGE_cable, canopy_tile, gs, gs_surft,           &
      T1P5M_surft, Q1P5M_surft, canopy_gb, MELT_surft,                        &
      !NPP, NPP_FT, GPP, GPP_FT, RESP_S, RESP_S_TOT,              &
      !RESP_P, RESP_P_FT, G_LEAF, 
      tl_1, qw_1, SURF_HTF_surft,                                             &
      !C_pool_casa, N_pool_casa, P_pool_casa, LAI_casa, PHENPHASE_casa,       &
      !NPP_PFT_ACC, RSP_W_PFT_ACC,
      dtrad,                                                                  &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soil_cbl )
 
CALL implicit_unpack( cycleno, row_length,rows, land_pts, ntiles,             &
                      L_tile_pts, npft,                                       &
                     sm_levels, dim_cs1, dim_cs2, t_soil, soil_temp_cable,    &
                     smcl, soil_moist_cable, &!aquifer_moist_cable,           &
                     SMVCST_cable, sthf,                                      &
                     soil_froz_frac_cable, sthu, STHU_surft, snow_surft,      &
                     snow_avg_rho_cable, isnow_flg_cable, snow_dpth_cable,    &
                     snow_mass_cable, snow_rho_cable, snow_temp_cable,        &
                     snow_cond, ftl_1, FTL_surft, fqw_1, FQW_surft,           &
                     TSTAR_surft, surf_ht_flux_land, ECAN_surft,              &
                     ESOIL_surft, EI_surft, RADNET_surft, tot_alb,            &
                     SNOW_AGE_cable, canopy_tile, gs,GS_surft,                     &
                     T1P5M_surft, Q1P5M_surft,                                &
                     canopy_gb, fland, MELT_surft,                            &
                     !NPP, NPP_FT, GPP, GPP_FT, RESP_S, RESP_S_TOT,          &
                     !RESP_P, RESP_P_FT, G_LEAF, TRANSP_surft,  &
                     !NPP_PFT_ACC, RSP_W_PFT_ACC,
                     SURF_HTF_surft,                                          &
                     dtrad, dtstar_surft,                                     &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soil_cbl )

snow_flg_cable = REAL(isnow_flg_cable)
t_soil_soilt(:,1,:) = t_soil(:,:) 
smcl_soilt(:,1,:)   = smcl(:,:)   
sthu_soilt(:,1,:)   = sthu(:,:)   
sthf_soilt(:,1,:)   = sthf(:,:)   
    
first_call = .FALSE.        

cable_runtime%um_implicit = .FALSE.

!-------- End Unique subroutine body -----------

dummy(1) =  SUM(snow_surft(:,:) + ftl_surft(:,:) + fqw_surft(:,:) + tstar_surft(:,:) + &
surf_htf_surft(:,:) + radnet_surft(:,:) + dtstar_surft(:,:) ) 
 
RETURN

END SUBROUTINE cable_implicit_main
  
END MODULE cable_implicit_main_mod


