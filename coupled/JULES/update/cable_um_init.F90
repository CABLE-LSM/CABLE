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
! Purpose: Initialize and update CABLE variables from UM forcing, calls to 
!          memory allocation and initialization subroutines
!
! Called from: cable_explicit_driver
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: No significant change from v1.8
! ==============================================================================

MODULE cable_um_init_mod
IMPLICIT NONE   
PUBLIC :: interface_UM_data
PRIVATE  

CONTAINS

SUBROUTINE interface_UM_data( mp, row_length, rows, land_pts, ntiles,         &
                              npft, sm_levels, itimestep, latitude, longitude,&
                              land_index, tile_frac, tile_pts, tile_index,    &
                              bexp, hcon, satcon, sathh, smvcst, smvcwt,      &
                              smvccl, albsoil, slope_avg,slope_std,dz_gw,     &
                              perm_gw,drain_gw,snow_tile, snow_rho1l,         &
                              snow_age, isnow_flg3l, snow_rho3l, snow_cond,   &
                              snow_depth3l, snow_tmp3l, snow_mass3l, sw_down, &
                              lw_down, cos_zenith_angle, surf_down_sw, ls_rain, &
                              ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,   &
                              z1_uv, rho_water, rho_ice,L_tile_pts,           &
                              visc_sublayer_depth, canopy_tile, Fland,        &
                              co2_mmr, & !r935 CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE
                              sthu_tile, smcl_tile, smgw_tile,sthf_tile, sthu,&
                              tsoil_tile, canht_ft, lai_ft, sin_theta_latitude, &
                              dzsoil, cpool_tile, npool_tile, ppool_tile,     &
                              soil_order, nidep, nifix, pwea, pdust, glai,    &
                              phenphase, npp_ft_acc, resp_w_ft_acc,           &
                              air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl, &
                              ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl,      &
                              veg_cbl, soilin, soil_cbl )
!subrs:
USE cable_um_init_subrs_mod, ONLY: initialize_veg
USE cable_um_init_subrs_mod, ONLY: initialize_roughness
USE cable_um_init_subrs_mod, ONLY: initialize_soilsnow
USE cable_um_init_subrs_mod, ONLY: initialize_canopy
USE cable_um_init_subrs_mod, ONLY: initialize_radiation
USE cable_um_init_subrs_mod, ONLY: init_bgc_vars
USE cable_um_init_subrs_mod, ONLY: init_sumflux_zero
!H!USE cable_um_init_subrs_mod, ONLY: init_respiration
USE cable_um_init_soil_mod,  ONLY: initialize_soil
USE cable_um_tech_mod,       ONLY: alloc_um_interface_types ! mem. allocation subr (um1%) 
USE cbl_soil_snow_init_special_module, ONLY: spec_init_soil_snow
!H!  USE casa_um_inout_mod
!data:   
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
USE cable_params_mod,         ONLY: soilin_type
USE cable_params_mod,         ONLY: soil_parameter_type
USE cable_um_tech_mod,        ONLY: um1                      ! um1% type UM basics 4 convenience
USE cable_common_module,      ONLY: cable_user               ! namelist user  configs
USE cable_common_module,      ONLY: l_casacnp
USE cable_common_module,      ONLY: knode_gl
USE cable_common_module,      ONLY: kwidth_gl
USE cable_other_constants_mod, ONLY : CLAI_THRESH => LAI_THRESH

!-------------------------------------------------------------------------- 
!--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
!-------------------------------------------------------------------------- 
!___UM dimensions, array indexes, flags
INTEGER ::                                                                    &
   mp,               & ! active pts CABLE 
   row_length, rows, & ! UM resolution
   land_pts,         & ! number of land_pts
   ntiles,           & ! number of tiles
   npft,             & ! number of Plant Functional Types
   sm_levels           ! number of soil layers

INTEGER, DIMENSION(land_pts) ::                                               &
   land_index  ! index of land point 

INTEGER, DIMENSION(ntiles) ::                                                 &
   tile_pts    ! number of land_pts per tile type

INTEGER, DIMENSION(land_pts, ntiles) ::                                       &
   tile_index, &  ! index of tile 
   isnow_flg3l    ! flag for 3-layer snow 

!___UM parameters 
INTEGER :: itimestep
REAL :: rho_water ,rho_ice
REAL, DIMENSION(sm_levels) ::                                                 &
   dzsoil

!___UM soil/snow/radiation/met vars
REAL, DIMENSION(land_pts) ::                                                  &
   bexp, &     !
   hcon, &     !  
   satcon, &   !
   sathh,  &   !
   smvcst, &   !
   smvcwt, &   !
   smvccl, &   !
   albsoil,&   !
   fland       !
    
REAL, DIMENSION(land_pts) ::                                                  &
   slope_avg,                                                                 &
   slope_std,                                                                 &
   dz_gw,                                                                     &
   sy_gw,                                                                     &
   perm_gw,                                                                   &
   drain_gw

REAL, DIMENSION(row_length,rows) ::                                           &
   sw_down, &        !
   cos_zenith_angle  !

REAL, DIMENSION(row_length,rows) ::                                           &
   latitude, longitude,                                                       &
   lw_down, &  !
   ls_rain, &  !
   ls_snow, &  !   
   tl_1,    &  !
   qw_1,    &  !
   vshr_land,& !
   pstar,   &  !
   z1_tq,   &  !
   z1_uv       !

REAL, DIMENSION(land_pts, ntiles) ::                                          &
   snow_tile   !

REAL, DIMENSION(land_pts, ntiles) ::                                          &
   tile_frac, &   !   
   snow_rho1l,&   !
   snow_age     !

REAL,DIMENSION(land_pts, ntiles) ::                                           &
   surf_down_sw 

REAL, DIMENSION(land_pts, npft) ::                                            &
   canht_ft, & !
   lai_ft      !
REAL,DIMENSION(land_pts, ntiles) ::                                           &
   canopy_tile !

REAL,DIMENSION(land_pts, ntiles) ::                                           &
   visc_sublayer_depth !
REAL, DIMENSION(land_pts, ntiles,3) ::                                        &
   snow_cond   !

REAL, DIMENSION(land_pts, ntiles,3) ::                                        &
   snow_rho3l, &     !
   snow_depth3l, &   ! 
   snow_mass3l,  &   ! 
   snow_tmp3l

REAL, DIMENSION(land_pts, sm_levels) ::                                       &
   sthu  !

REAL, DIMENSION(land_pts, ntiles, sm_levels) ::                               &
   sthu_tile, &   !
   sthf_tile, &   !
   smcl_tile, &   !
   tsoil_tile     !

REAL, DIMENSION(land_pts, ntiles) ::                                          &
   smgw_tile

REAL :: co2_mmr
LOGICAL  :: l_co2_interactive
INTEGER ::                                                                    &
   co2_dim_len                                                                &
  ,co2_dim_row 
REAL :: co2_3d(1,1)  ! co2 mass mixing ratio

LOGICAL,DIMENSION(land_pts, ntiles) ::                                        &
   L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0

REAL, DIMENSION(row_length,rows) ::                                           &
   sin_theta_latitude

REAL, DIMENSION(land_pts,ntiles,10) ::                                        &
   cpool_tile,                                                                &
   npool_tile

REAL, DIMENSION(land_pts,ntiles,12) ::                                        &
   ppool_tile

REAL, DIMENSION(land_pts) ::                                                  &
   soil_order,                                                                &
   nidep,                                                                     &
   nifix,                                                                     &
   pwea,                                                                      &
   pdust

REAL, DIMENSION(land_pts,ntiles) ::                                           &
   glai,                                                                      &
   phenphase

REAL, DIMENSION(land_pts,ntiles) ::                                           &
   npp_ft_acc,                                                                &
   resp_w_ft_acc

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
TYPE(soilin_type),  INTENT(INOUT) ::  soilin  
TYPE(soil_parameter_type),  INTENT(INOUT) ::  soil_cbl
!------------------------------------------------------------------------- 
!--- end INPUT ARGS FROM cable_explicit_driver() -------------------------
!------------------------------------------------------------------------- 

!___ local decs
!___defs 1st call to CABLE in this run. necesary in UM & coupled.
LOGICAL, SAVE :: first_call = .TRUE. 
INTEGER :: i,j

INTEGER, DIMENSION(land_pts, ntiles) ::                                       &
  tile_index_mp   !tile # on each land point, ntile for packing to CABLE

REAL, DIMENSION(land_pts, ntiles) ::  clobbered_htveg

!--- logn, vegparmnew can be set thru cable.nml
INTEGER :: logn = 6       ! 6=write to std out
LOGICAL :: vegparmnew = .TRUE.   ! true=read std veg params false=CASA file 

  !---------------------------------------------------------------------!
  !--- code to create type um1% conaining UM basic vars describing    --! 
  !--- dimensions etc, which are used frequently throughout interfaces. !
  !--- and then assign every time explicit is called (shouldn't need to)!
  !---------------------------------------------------------------------!
CALL alloc_um_interface_types( row_length, rows, land_pts,                    &
                                  ntiles, sm_levels )
um1%l_tile_pts = l_tile_pts 

CALL assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,             &
                              npft, sm_levels, itimestep, latitude,           &
                              longitude, land_index, tile_frac,               &
                              tile_pts, tile_index, l_tile_pts,               &
                              rho_water,rho_ice  )

!---------------------------------------------------------------------!
!--- CABLE vars are initialized/updated from passed UM vars       ----!
!---------------------------------------------------------------------!

CALL initialize_veg( clobbered_htveg, land_pts, npft, ntiles, sm_levels, mp,     &
                     canht_ft, lai_ft, dzsoil, veg_cbl,         &
                    tile_pts, tile_index, tile_frac, L_tile_pts,                 &
                    CLAI_thresh )

CALL initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,              &
                      smvccl, albsoil, tsoil_tile, sthu, sthu_tile,           &
                      dzsoil, slope_avg, slope_std, dz_gw,                    &
                      perm_gw,drain_gw, soilin, veg_cbl, soil_cbl,            &
                      ssnow_cbl ) 

CALL initialize_roughness( z1_tq, z1_uv, clobbered_htveg, rough_cbl,          &
                           veg_cbl ) 
   
CALL initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile,           &
                          smgw_tile, snow_tile, snow_rho1l, snow_age,         &
                          isnow_flg3l, snow_rho3l, snow_cond,                 &
                          snow_depth3l, snow_mass3l, snow_tmp3l,              &
                          fland, sin_theta_latitude, soil_cbl, ssnow_cbl,     &
                          met_cbl, bal_cbl, veg_cbl ) 

CALL initialize_canopy(canopy_tile,visc_sublayer_depth, canopy_cbl)

CALL initialize_radiation( sw_down, lw_down, cos_zenith_angle,                &
                           surf_down_sw, sin_theta_latitude, ls_rain,         &
                           ls_snow, tl_1, qw_1, vshr_land, pstar,             &
                           co2_mmr,co2_3d,co2_dim_len,co2_dim_row,            &
                           l_co2_interactive, rad_cbl, met_cbl, soil_cbl )   

IF ( first_call ) THEN
call spec_init_soil_snow( real(kwidth_gl), soil_cbl, ssnow_cbl, canopy_cbl, met_cbl, bal_cbl, veg_cbl)
  CALL init_bgc_vars( bgc_cbl, veg_cbl ) 
  CALL init_sumflux_zero(sum_flux_cbl) 

  !--- initialize respiration for CASA-CNP
  !IF (l_casacnp) THEN ?
  !H!CALL init_respiration(npp_ft_acc,resp_w_ft_acc, canopy_cbl)

  !H!  ! Lestevens 28 Sept 2012 - Initialize CASA-CNP here
  !H!     if (l_casacnp) then
  !H!       if (knode_gl==0) then
  !H!         print *, '  '; print *, 'CASA_log:'
  !H!         print *, '  Calling CasaCNP - Initialise '
  !H!         print *, '  l_casacnp = ',l_casacnp
  !H!         print *, 'End CASA_log:'; print *, '  '
  !H!       endif
  !H!       call init_casacnp(sin_theta_latitude,cpool_tile,npool_tile,&
  !H!                         ppool_tile,soil_order,nidep,nifix,pwea,pdust,&
  !H!                         GLAI,PHENPHASE)
  !H!     endif

  first_call = .FALSE. 
END IF      
RETURN
END SUBROUTINE interface_UM_data
                                   
!============================================================================

SUBROUTINE assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,       &
                                    npft, sm_levels, timestep, latitude,      &
                                    longitude, land_index, tile_frac,         &
                                    tile_pts, tile_index, l_tile_pts,         &
                                    rho_water,rho_ice  )
USE cable_um_tech_mod,   ONLY: um1
USE cable_common_module, ONLY: cable_user

INTEGER, INTENT(IN) :: row_length, rows, land_pts, ntiles, npft, sm_levels
INTEGER, INTENT(IN) :: timestep 
REAL, INTENT(IN), DIMENSION(row_length,rows) :: latitude, longitude 
REAL,INTENT(IN):: rho_water,rho_ice
INTEGER, INTENT(IN), DIMENSION(land_pts)  :: land_index 
INTEGER, INTENT(IN), DIMENSION(ntiles)  :: tile_pts 
INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles)  :: tile_index
REAL, INTENT(IN), DIMENSION(land_pts, ntiles)  :: tile_frac 
LOGICAL, INTENT(IN), DIMENSION(land_pts,ntiles)  :: l_tile_pts 
     
um1%row_length = row_length
um1%rows       = rows
um1%land_pts   = land_pts
um1%ntiles     = ntiles   
um1%npft       = npft   
um1%sm_levels  = sm_levels
um1%timestep   = timestep
um1%latitude   = latitude
um1%longitude  = longitude
um1%land_index = land_index
um1%tile_frac  = tile_frac
um1%tile_pts   = tile_pts
um1%tile_index = tile_index
um1%rho_water  = rho_water 
um1%rho_ice    = rho_ice

END SUBROUTINE assign_um_basics_to_um1

END MODULE cable_um_init_mod
