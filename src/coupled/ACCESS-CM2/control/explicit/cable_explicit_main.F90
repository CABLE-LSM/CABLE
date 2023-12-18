!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose:
!
! Called from: JULES: surf_couple_
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_explicit_main_mod
  
contains

SUBROUTINE cable_explicit_main(                                                &
# if defined(UM_JULES)
            ! grid, model, dimensions. PFT frac per landpoint etc   
            mype, timestep_width, timestep_number, endstep,                    &
            cycleno, numcycles,                                                &
            row_length, rows, land_pts, ntiles,                                &
            npft, sm_levels,                                                   &
            latitude, longitude,                                               &
            land_index, tile_frac, tile_pts, tile_index,                       &
            ! Soil parameters **jhan: could be undef.@some point  issue here
            bexp, hcon, satcon, sathh,                                         &
            smvcst, smvcwt, smvccl, albsoil,                                   &
            ! packs/ unpacked ssnow% snowd 
            snow_tile, lw_down, cosine_zenith_angle,                           &
            surf_down_sw, ls_rain, ls_snow,                                    &
            tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,           &
            Fland, CO2_MMR, sthu, canht_ft, lai_ft ,                           &
            sin_theta_latitude, dzsoil, FTL_TILE, FQW_TILE,                    &
            TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE,                   &
            RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,              &
            RECIP_L_MO_TILE, EPOT_TILE &
# else  
            !jhan:this looks like I was testing standalone
            ,timestep_number                                                   &
# endif
            )
  !subrs called 
  USE cable_explicit_driv_mod, ONLY : cable_explicit_driver
  USE cable_expl_unpack_mod, ONLY : cable_expl_unpack

  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
!jhan:this looks like I was testing standalone
# if !defined(UM_JULES)
  USE model_grid_mod, ONLY: latitude, longitude
# endif
# if defined(UM_JULES)
  USE atm_fields_real_mod, ONLY : soil_temp_cable, soil_moist_cable,           &
                                  soil_froz_frac_cable, snow_dpth_cable,       &
                                  snow_mass_cable, snow_temp_cable,            &
                                  snow_rho_cable, snow_avg_rho_cable,          &
                                  snow_age_cable, snow_flg_cable,              &
                                  C_pool_casa, N_pool_casa, P_pool_casa,       &
                                  SOIL_ORDER_casa, N_DEP_casa, N_FIX_casa,     &
                                  P_DUST_casa, P_weath_casa, LAI_casa,         &
                                  PHENPHASE_casa, NPP_PFT_ACC, RSP_W_PFT_ACC,  &
                                  aquifer_moist_cable,aquifer_thickness_cable, &
                                  slope_avg_cable,slope_std_cable,&
                                  visc_sublayer_depth,aquifer_perm_cable,&
                                  aquifer_draindens_cable
#endif
  !make availble these vars and CABLE % types - to be tidied
  USE cable_decs_mod, only : L_tile_pts, rho_water
  USE cable_um_tech_mod, ONLY : air, bgc, canopy, met, bal, rad, rough, soil,  &
                                ssnow, sum_flux, veg
  implicit none
 
  !___ re-decl input args
  integer :: mype, timestep_number, endstep, cycleno, numcycles
  real :: timestep_width
  INTEGER ::                                                      & 
    row_length, rows, & ! UM grid resolution
    land_pts,         & ! # of land points being processed
    ntiles,           & ! # of tiles 
    npft,             & ! # of plant functional types
    sm_levels          ! # of soil layers 

# if defined(UM_JULES)
  REAL,  DIMENSION(row_length,rows) :: latitude,  longitude
# endif   
  INTEGER, DIMENSION(land_pts) :: land_index  ! index of land points processed

  INTEGER,  DIMENSION(ntiles) :: tile_pts ! # of land points on each tile

  INTEGER,  DIMENSION(land_pts, ntiles) :: tile_index ! index of tile points 

  REAL, DIMENSION(land_pts, ntiles) :: tile_frac
   
  !___UM soil/snow/radiation/met vars
  REAL,  DIMENSION(land_pts) :: & 
    bexp,    & ! => parameter b in Campbell equation 
    hcon,    & ! Soil thermal conductivity (W/m/K).
    satcon,  & ! hydraulic conductivity @ saturation [mm/s]
    sathh,   &
    smvcst,  &
    smvcwt,  &
    smvccl,  &
    albsoil

  REAL,  DIMENSION(land_pts, ntiles) :: snow_tile
   
  REAL,  DIMENSION(row_length,rows) ::                                         &
   cosine_zenith_angle,   &
   lw_down,               &
   ls_rain,               &
   ls_snow

  REAL, DIMENSION(row_length, rows, 4) ::                         &
    surf_down_sw 
  
  REAL,  DIMENSION(row_length,rows) ::                             &
    tl_1,       &
    qw_1,       &  
    vshr_land,  &
    pstar,      &
    z1_tq,      &
    z1_uv

  REAL, DIMENSION(land_pts, ntiles) :: canopy_tile
  
  REAL, DIMENSION(land_pts) :: fland
  
  REAL :: co2_mmr

  REAL, DIMENSION(land_pts, sm_levels) :: sthu 
   
  REAL, DIMENSION(land_pts, npft) :: canht_ft, lai_ft 
  
  REAL :: sin_theta_latitude(row_length,rows) 
    
  REAL,  DIMENSION(sm_levels) :: dzsoil

  !___return miscelaneous 
  REAL, DIMENSION(land_pts,ntiles) :: &
    FTL_TILE,    & ! Surface FTL for land tiles     
    FQW_TILE,    & ! Surface FQW for land tiles     
    TSTAR_TILE,  & ! radiative temperature of surface 
    U_S,         & ! Surface friction velocity (m/s)
    U_S_STD_TILE,& ! Surface friction velocity
    CD_TILE,     & ! Drag coefficient
    CH_TILE,     &  ! Transfer coefficient for heat & moisture
    RADNET_TILE, & ! Surface net radiation
    FRACA,       & ! Fraction of surface moisture
    RESFS,       & ! Combined soil, stomatal & aerodynamic resistance
                   ! factor for fraction (1-FRACA) of snow-free land tiles
    RESFT,       & ! Total resistance factor.
                   ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                   ! 1 for snow.    
    Z0H_TILE,    &
    Z0M_TILE,    &
    RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
    EPOT_TILE

!jhan:this looks like I was testing standalone
# if !defined(UM_JULES)
  integer :: timestep_number  
# endif

  !___ local vars
  integer,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable
  logical, save :: first_call = .true.
  real :: radians_degrees
  REAL,  DIMENSION(row_length,rows) :: latitude_deg, longitude_deg

  ! std template args 
  character(len=*), parameter :: subr_name = "cable_explicit_main"

  !-------- Unique subroutine body -----------
  
  !--- initialize cable_runtime% switches 
  cable_runtime%um =          .TRUE.
  cable_runtime%um_explicit = .TRUE.
   
  ! initialize processor number, timestep width & number, endstep 
  ! UM overwrites these defaults. Not yet done in StandAlone 
  if( first_call ) then
    knode_gl =0; kwidth_gl = 1200.; kend_gl=-1
# if defined(UM_JULES)
    knode_gl  = mype 
    kwidth_gl = int(timestep_width)
    kend_gl   = endstep   
# endif

  !--- Convert lat/long to degrees
    radians_degrees = 180.0 / ( 4.0*atan(1.0) ) ! 180 / PI
  latitude_deg  = latitude * radians_degrees
  longitude_deg  = longitude * radians_degrees
  endif

  ktau_gl   = timestep_number
  
  if( .NOT. allocated(L_tile_pts) ) allocate( L_tile_pts(land_pts, ntiles) ) 

  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  isnow_flg_cable = int(snow_flg_cable)

  call cable_explicit_driver( row_length, rows, land_pts, ntiles,npft,         &
                              sm_levels, timestep_width,                       &
                              latitude_deg, longitude_deg,                     &
                              land_index, tile_frac,  tile_pts, tile_index,    &
                              bexp, hcon, satcon, sathh, smvcst,               &
                              smvcwt,  smvccl, albsoil,                        &
                              slope_avg_cable,slope_std_cable,                 &
                              aquifer_thickness_cable,                         &
                              aquifer_perm_cable,aquifer_draindens_cable,      &
                              snow_tile,    &
                              snow_avg_rho_cable, snow_age_cable,              &
                              isnow_flg_cable, snow_rho_cable, snow_dpth_cable,&
                              snow_temp_cable, snow_mass_cable,                &
                              lw_down, cosine_zenith_angle, surf_down_sw,      &
                              ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar,  &
                              z1_tq, z1_uv,visc_sublayer_depth,  canopy_tile,  &
                              Fland, CO2_MMR,                                  &
                              soil_moist_cable, aquifer_moist_cable,           &
                              soil_froz_frac_cable, sthu,                      &
                              soil_temp_cable, canht_ft, lai_ft,               &
                              sin_theta_latitude, dzsoil, FTL_TILE, FQW_TILE,  &
                              TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE, &
                              RADNET_TILE, FRACA, RESFS, RESFT,                &
                              Z0H_TILE, Z0M_TILE,                              &
                              RECIP_L_MO_TILE, EPOT_TILE,                      &
                              C_pool_casa, N_pool_casa, P_pool_casa,           &
                              SOIL_ORDER_casa, N_DEP_casa, N_FIX_casa,         &
                              P_weath_casa, P_DUST_casa, LAI_casa,             &
                              PHENPHASE_casa, NPP_PFT_ACC, RSP_W_PFT_ACC,      &
                              endstep, timestep_number, mype )    

  
  !----------------------------------------------------------------------------
  !--- CALL _unpack to unpack variables from CABLE back to UM format to return
  !----------------------------------------------------------------------------
   call cable_expl_unpack( latitude, longitude, FTL_TILE, FQW_TILE, TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq, &
                           canopy%fes, canopy%fev )

  !-------- End Unique subroutine body -----------
  
  cable_runtime%um_explicit = .FALSE.
  first_call = .false.        

return

End subroutine cable_explicit_main
  
End module cable_explicit_main_mod

