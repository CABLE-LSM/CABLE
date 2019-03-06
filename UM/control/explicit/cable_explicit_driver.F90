!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
! 
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located 
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: UM code sf_exch
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

module cable_explicit_driv_mod
  
contains

SUBROUTINE cable_explicit_driver( row_length, rows, land_pts, ntiles,npft,     &
                                  sm_levels, timestep, latitude, longitude,    &
                                  land_index, tile_frac,  tile_pts, tile_index,&
                                  bexp, hcon, satcon, sathh, smvcst,           &
                                  smvcwt, smvccl, albsoil, slope_avg,slope_std,&
                                  dz_gw,perm_gw,drain_gw,snow_tile,            &
                                  snow_rho1l, snow_age, snow_flg3l, snow_rho3l,&
                                  snow_depth3l, snow_tmp3l, snow_mass3l,       &
                                  lw_down, cos_zenith_angle, surf_down_sw,     &
                                  ls_rain, ls_snow, tl_1, qw_1, vshr_land,     &
                                  pstar, z1_tq, z1_uv, visc_sublayer_depth,    &
                                  canopy_tile, Fland, CO2_MMR,                 &
                                  ! + Extra atmospheric co2 variables
                                  !CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,
                                  !L_CO2_INTERACTIVE, 
                                  smcl_tile, smgw_tile,sthf_tile, sthu,        &
                                  tsoil_tile, canht_ft, lai_ft,                &
                                  sin_theta_latitude, dzsoil,                  &
                                  FTL_TILE, FQW_TILE, TSTAR_TILE,              &
                                  U_S, U_S_STD_TILE, CD_TILE, CH_TILE,         &
                                  RADNET_TILE, FRACA, RESFS, RESFT,            &
                                  Z0H_TILE, Z0M_TILE,                          &
                                  RECIP_L_MO_TILE, EPOT_TILE,                  &
                                  CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,          &
                                  SOIL_ORDER, NIDEP, NIFIX, PWEA, PDUST,       &
                                  GLAI, PHENPHASE, NPP_FT_ACC, RESP_W_FT_ACC,  &
                                  endstep, timestep_number, mype )    
  !subrs called 
  USE cable_um_init_mod, ONLY : interface_UM_data
  USE cable_cbm_module, ONLY : cbm

  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir
  
  USE cable_diag_module
  
  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
  !--- vars common to CABLE declared 
  USE cable_common_module!, ONLY : cable_runtime, cable_user, ktau_gl,         &
                         !         knode_gl, kwidth_gl, kend_gl,               &
                         !         report_version_no,                          &
                         !         l_vcmaxFeedbk, l_laiFeedbk
   
  USE cable_data_module, ONLY : cable

  !--- reads runtime and user switches and reports
  USE cable_um_tech_mod, ONLY : cable_um_runtime_vars, air, bgc, canopy,      &
                                met, bal, rad, rough, soil, ssnow, sum_flux,  &
                                veg, basic_diag 
  
  USE cable_def_types_mod, ONLY : mp, ms

  USE cable_decs_mod, only : L_tile_pts!, rho_water
  USE cable_data_module, only : PHYS

  USE casavariable
  USE casa_types_mod
  USE cable_climate_mod
  !made this a module so can build in the UM re:dependencies etc
  !did the same dor sli_main. promote everything to modules
  USE casa_cable
   
  implicit none
  
  !___ re-decl input args
   
  !___IN: UM dimensions, array indexes, flags
  INTEGER ::                                                      & 
    row_length, rows, & ! UM grid resolution
    land_pts,         & ! # of land points being processed
    ntiles,           & ! # of tiles 
    npft,             & ! # of plant functional types
    sm_levels           ! # of soil layers 

  ! index of land points being processed
  INTEGER, DIMENSION(land_pts) :: land_index 

  ! # of land points on each tile
  INTEGER,  DIMENSION(ntiles) :: tile_pts 
  
  INTEGER,  DIMENSION(land_pts, ntiles) ::                         & 
    tile_index ,& ! index of tile points being processed
    snow_flg3l   ! 3 layer snow flag

  !___UM parameters: water density, soil layer thicknesses 
  REAL,  DIMENSION(sm_levels) :: dzsoil

  !___UM soil/snow/radiation/met vars
  REAL,  DIMENSION(land_pts) :: & 
    bexp,    & ! => parameter b in Campbell equation 
    hcon,    & ! Soil thermal conductivity (W/m/K).
    satcon,  & ! hydraulic conductivity @ saturation [mm/s]
    sathh,   &
    smvcst,  &
    smvcwt,  &
    smvccl,  &
    albsoil, &
    fland 
  
  REAL,  DIMENSION(land_pts) :: & 
    slope_avg,&
    slope_std,&
      dz_gw,sy_gw,perm_gw,drain_gw

  
  REAL,  DIMENSION(row_length,rows) :: &
    cos_zenith_angle, & 
    latitude,   &
    longitude,  &
    sw_down,    & ! NOT SW forcing. surf_down_sw IS 
    lw_down,    &
    ls_rain,    &
    ls_snow,    &
    tl_1,       &
    qw_1,       &  
    vshr_land,  &
    pstar,      &
    z1_tq,      &
    z1_uv
  
  REAL, DIMENSION(row_length, rows, 4) ::                         &
    surf_down_sw 
   
  REAL,  DIMENSION(land_pts, ntiles) ::                         &
    snow_tile,    &     
    tile_frac,    &    
    snow_rho1l,   &
    snow_age,     &
    canopy_tile,  &
    visc_sublayer_depth 

  REAL, DIMENSION(land_pts, npft) ::                              &
    canht_ft, lai_ft 
   
  REAL,  DIMENSION(land_pts, ntiles,3) ::                       &
    snow_cond,     &
    snow_rho3l,    &
    snow_depth3l,  &
    snow_mass3l,   &
    snow_tmp3l
   
  REAL, DIMENSION(land_pts, sm_levels) ::                         &
    sthu 
   
  REAL, DIMENSION(land_pts, ntiles, sm_levels) :: & 
    sthu_tile, &
    sthf_tile, &
    smcl_tile, &
    tsoil_tile
   
  REAL :: co2_mmr
   
  ! rml 2/7/13 Extra atmospheric co2 variables
  !LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
  !INTEGER, INTENT(IN) ::                              &
  !  CO2_DIM_LEN                                      &
  !  ,CO2_DIM_ROW
  !REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio
  
  REAL :: sin_theta_latitude(row_length,rows) 
    
  !___return fluxes AND miscelaneous 
  REAL, DIMENSION(land_pts,ntiles) :: &
    FTL_TILE,   &  ! Surface FTL for land tiles     
    FQW_TILE,   &  ! Surface FQW for land tiles     
    TSTAR_TILE,       & 
    Z0H_TILE,         &
    Z0M_TILE,   &
    CD_TILE,    &     ! Drag coefficient
    CH_TILE,    &     ! Transfer coefficient for heat & moisture
    U_S_STD_TILE, &      ! Surface friction velocity
    RADNET_TILE,   & ! Surface net radiation
    RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                     ! factor for fraction (1-FRACA) of snow-free land tiles
    RESFT,         & ! Total resistance factor.
                     ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                     ! 1 for snow.    
    FRACA,         & ! Fraction of surface moisture
    RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
    EPOT_TILE,    &
    smgw_tile

  REAL, DIMENSION(row_length,rows)  :: &
    U_S               ! Surface friction velocity (m/s)
   
  ! end step of experiment, this step, step width, processor num
  INTEGER :: endstep, timestep_number, mype
  REAL ::  timestep     
   
  !CASA vars 
  REAL, DIMENSION(land_pts,ntiles,10) :: &
    CPOOL_TILE,    & ! Carbon Pools
    NPOOL_TILE       ! Nitrogen Pools

  REAL, DIMENSION(land_pts,ntiles,12) :: &
    PPOOL_TILE       ! Phosphorus Pools

  REAL, DIMENSION(land_pts) :: &
    SOIL_ORDER,    & ! Soil Order (1 to 12)
    NIDEP,         & ! Nitrogen Deposition
    NIFIX,         & ! Nitrogen Fixation
    PWEA,          & ! Phosphorus from Weathering
    PDUST            ! Phosphorus from Dust

  REAL, DIMENSION(land_pts,ntiles) :: &
    GLAI,         &  ! Leaf Area Index for Prognostics LAI
    PHENPHASE,    &  ! Phenology Phase for Casa-CNP
    NPP_FT_ACC,   &
    RESP_W_FT_ACC
  
  !___ local vars
  !jhan: this can be moved and USEd ?
  TYPE (climate_type) :: climate     ! climate variables
   
  CHARACTER(LEN=200), PARAMETER ::                                            &
      runtime_vars_file = 'cable.nml'! path/name namelist def runtime vars

  !___ 1st call in RUN (!=ktau_gl -see below) 
  LOGICAL, SAVE :: first_cable_call = .TRUE.

  real  :: rho_water, rho_ice
  
  ! std template args 
  character(len=*), parameter :: subr_name = "cable_explicit_driver"
   
# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
  
  !--- user FLAGS, variables etc def. in cable.nml is read on 
  !--- first time step of each run. these variables are read at 
  !--- runtime and for the most part do not require a model rebuild.
  IF(first_cable_call) THEN
     CALL cable_um_runtime_vars(runtime_vars_file) 
     first_cable_call = .FALSE.
  ENDIF      

  rho_water = PHYS%density_liq

  IF (cable_user%GW_model) then
     rho_ice = PHYS%density_ice
  ELSE
     rho_ice = PHYS%density_liq
  ENDIF      

  !---------------------------------------------------------------------!
  !--- initialize CABLE using UM forcings etc. these args are passed ---!
  !--- down from UM.                                                 ---! 
  !---------------------------------------------------------------------!
  CALL interface_UM_data( row_length, rows, land_pts, ntiles, npft,            &
                          sm_levels, ktau_gl, latitude, longitude,             &
                          land_index, tile_frac, tile_pts, tile_index,         &
                          bexp, hcon, satcon, sathh, smvcst, smvcwt,           &
                          smvccl, albsoil, slope_avg,slope_std,dz_gw,          &
                          perm_gw,drain_gw,snow_tile, snow_rho1l,              &
                          snow_age, snow_flg3l, snow_rho3l, snow_cond,         &
                          snow_depth3l, snow_tmp3l, snow_mass3l, sw_down,      &
                          lw_down, cos_zenith_angle, surf_down_sw, ls_rain,    &
                          ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,        &
                          z1_uv, rho_water, rho_ice,L_tile_pts,                &
                          visc_sublayer_depth, canopy_tile, Fland,             &
                          CO2_MMR, & ! pass 3d co2 through to cable if required
                          !CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                          sthu_tile, smcl_tile, smgw_tile,sthf_tile,           &
                          sthu, tsoil_tile, canht_ft, lai_ft,                  &
                          sin_theta_latitude, dzsoil,                          &
                          CPOOL_TILE, NPOOL_TILE, PPOOL_TILE, SOIL_ORDER,      &
                          NIDEP, NIFIX, PWEA, PDUST, GLAI, PHENPHASE,          &
                          NPP_FT_ACC,RESP_W_FT_ACC )


  !---------------------------------------------------------------------!
  !--- Feedback prognostic vcmax and daily LAI from casaCNP to CABLE ---!
  !---------------------------------------------------------------------!
  IF(l_vcmaxFeedbk) call casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)
  IF(l_laiFeedbk) veg%vlai(:) = casamet%glai(:)

  canopy%oldcansto=canopy%cansto
  rad%otrad = rad%trad

  !---------------------------------------------------------------------!
  !--- cbm "mainly" controls the calling of model components         ---!  
  !---------------------------------------------------------------------!
  CALL cbm( ktau_gl,timestep, air, bgc, canopy, met, bal,                      &
            rad, rough, soil, ssnow, sum_flux, veg, climate )
  
  !-------- End Unique subroutine body -----------

  fprintf_dir=trim(fprintf_dir_root)//trim(unique_subdir)//"/"
  if(L_cable_fprint) then 
    !basics to std output stream
    if (knode_gl == 0 .and. ktau_gl == 1)  call cable_fprintf(subr_name, .true.) 
    !more detailed output
    vname=trim(subr_name//'_')
    call cable_fprintf( cDiag00, vname, knode_gl, ktau_gl, .true. )
  endif

  if(L_cable_Pyfprint) then 
  endif

  return

End subroutine cable_explicit_driver

End module cable_explicit_driv_mod
