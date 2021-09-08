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
! Purpose: Routines to pass UM variables into appropriate CABLE variables and 
!          to map parameters for each surface type to CABLE arrays
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Rewrites of code in v1.8 (ACCESS1.3)
! ==============================================================================

MODULE cable_um_init_subrs_mod
   
IMPLICIT NONE

CONTAINS

!========================================================================
          
SUBROUTINE initialize_veg( clobbered_htveg, land_pts, npft, ntiles, ms, mp,      &
                           canht_ft, lai_ft, soil_zse, veg_cbl,                  &
                    tile_pts, tile_index, tile_frac, L_tile_pts,                 &
                    CLAI_thresh )

USE cable_params_mod,  ONLY: veg_parameter_type
USE cbl_LAI_canopy_height_mod,  ONLY: limit_HGT_LAI

INTEGER ::   mp                ! active pts CABLE
REAL :: CLAI_thresh  
INTEGER, DIMENSION(ntiles) ::                                                 &
   tile_pts    ! number of land_pts per tile type

REAL, DIMENSION(ntiles) ::                                                 &
   tile_frac
INTEGER, DIMENSION(land_pts, ntiles) ::                                       &
   tile_index     ! index of tile 
LOGICAL,DIMENSION(land_pts, ntiles) ::                                        &
   L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0

INTEGER :: land_pts
INTEGER :: ntiles
REAL :: clobbered_htveg(land_pts, ntiles)
INTEGER :: npft
INTEGER :: ms ! soil levels
TYPE(veg_parameter_type), INTENT(INOUT) :: veg_cbl
REAL, INTENT(IN) :: canht_ft(land_pts, npft)
REAL, INTENT(IN) :: lai_ft(land_pts, npft) 
REAL, DIMENSION(ms) :: soil_zse 

! limit IN height, LAI  and initialize existing cable % types
CALL limit_HGT_LAI( clobbered_htveg, veg_cbl%vlai, veg_cbl%hc, mp, land_pts, ntiles,           &
                    tile_pts, tile_index, tile_frac, L_tile_pts,              &
                    LAI_ft, canht_ft, CLAI_thresh )
     
END SUBROUTINE initialize_veg

       
SUBROUTINE initialize_radiation( sw_down, lw_down, cos_zenith_angle,          &
                                 surf_down_sw, sin_theta_latitude, ls_rain,   &
                                 ls_snow, tl_1, qw_1, vshr_land, pstar,       &
! rml 2/7/13 pass 3d co2 through to cable if required
                   co2_mmr,co2_3d,co2_dim_len,co2_dim_row,l_co2_interactive, rad_cbl, met_cbl , soil_cbl )   

USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_met_type_mod, ONLY: met_type
USE cable_params_mod, ONLY: soil_parameter_type
USE cable_other_constants_mod, ONLY : crad_thresh  => rad_thresh
USE cable_def_types_mod, ONLY: mp, mstype
USE cable_um_tech_mod,   ONLY: um1, conv_rain_prevstep, conv_snow_prevstep
USE cable_common_module, ONLY: cable_runtime, cable_user, ktau_gl

TYPE(radiation_type),     INTENT(INOUT)  :: rad_cbl  
TYPE(met_type),           INTENT(INOUT)  :: met_cbl  
TYPE(soil_parameter_type),     INTENT(INOUT)  :: soil_cbl 

REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: sw_down
   
REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::                      &
   lw_down,                                                                   &
   sin_theta_latitude
   
REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: cos_zenith_angle

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) :: surf_down_sw 
   
REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::                      &
   ls_rain,                                                                   &
   ls_snow,                                                                   &
   tl_1,                                                                      &
   qw_1,                                                                      &
   vshr_land,                                                                 &
   pstar
   
REAL, INTENT(IN) :: co2_mmr
! rml 2/7/13 Extra atmospheric co2 variables
LOGICAL, INTENT(IN) :: l_co2_interactive
INTEGER, INTENT(IN) ::                                                        &
   co2_dim_len                                                                &
  ,co2_dim_row
REAL, INTENT(IN) :: co2_3d(:,:)  ! co2 mass mixing ratio
             
!___defs 1st call to CABLE in this run. OK in UM & coupled
LOGICAL, SAVE :: first_call= .TRUE.
     
LOGICAL :: skip =.TRUE. 
REAL, DIMENSION(mstype) :: dummy 
INTEGER :: i,j
INTEGER, SAVE :: k = 0
      
dummy = 0.0 
      
IF ( first_call ) THEN
  rad_cbl%albedo_T = soil_cbl%albsoil(:,1)
  ALLOCATE( conv_rain_prevstep(mp), conv_snow_prevstep(mp) )
  conv_rain_prevstep = 0.0 
  conv_snow_prevstep = 0.0
END IF   
      
CALL um2cable_rr( cos_zenith_angle, met_cbl%coszen)
CALL um2cable_lp( 0.5 * surf_down_sw, dummy , met_cbl%fsd(:,1), soil_cbl%isoilm, skip )
CALL um2cable_lp( 0.5 * surf_down_sw, dummy , met_cbl%fsd(:,2), soil_cbl%isoilm, skip )
!H!DO i = 1, mp
!H!  IF ( met_cbl%fsd(i,1) < ( .1 * crad_thresh) ) met_cbl%fsd(i,1) = 0.0
!H!  IF ( met_cbl%fsd(i,2) < ( .1 * crad_thresh) ) met_cbl%fsd(i,2) = 0.0
!H!END DO
CALL um2cable_rr( lw_down, met_cbl%fld)
CALL um2cable_rr( (ls_rain * um1%timestep), met_cbl%precip)
CALL um2cable_rr( (ls_snow * um1%timestep), met_cbl%precip_sn)
CALL um2cable_rr( tl_1, met_cbl%tk)
CALL um2cable_rr( qw_1, met_cbl%qv)
CALL um2cable_rr( vshr_land, met_cbl%ua)
CALL um2cable_rr( pstar * 0.01, met_cbl%pmb)
      
!---re-set some of CABLE's forcing variables
met_cbl%precip   =  met_cbl%precip + met_cbl%precip_sn 
!met%precip   =  (met%precip + conv_rain_prevstep) &
!               + (met%precip_sn +  conv_snow_prevstep)
!               + (met%precip_sn +  conv_rain_prevstep)
met_cbl%tvair =     met_cbl%tk
met_cbl%tvrad =     met_cbl%tk
!H!met_cbl%coszen =    MAX(met_cbl%coszen,1e-8)

!initialise rad%trad on first call only
IF (first_call) THEN
  rad_cbl%trad = met_cbl%tk
  first_call = .FALSE.
END IF 

!---this is necessary clobrring at present 
!H!WHERE (met_cbl%ua < 0.001 ) met_cbl%ua = 0.001
      
!H!! rml 24/2/11 Set atmospheric CO2 seen by cable to CO2_MMR (value seen 
!H!! by radiation scheme).  Option in future to have cable see interactive 
!H!! (3d) CO2 field Convert CO2 from kg/kg to mol/mol ( m_air, 
!H!! 28.966 taken from include/constant/ccarbon.h file )
!H!! r935 rml 2/7/13 Add in co2_interactive option
!H!!IF (L_CO2_INTERACTIVE) THEN
!H!!  CALL um2cable_rr(CO2_3D, met%ca)
!H!!ELSE
!H!met_cbl%ca = co2_mmr
!H!!ENDIF
!H!met_cbl%ca = met_cbl%ca * 28.966 / 44.0

!H!WHERE (met_cbl%coszen < crad_thresh ) 
!H!  rad_cbl%fbeam(:,1) = REAL(0) 
!H!  rad_cbl%fbeam(:,2) = REAL(0) 
!H!  rad_cbl%fbeam(:,3) = REAL(0) 
!H!END WHERE

!--- CABLE radiation type forcings, not set by um2cable_met_rad(
!--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
CALL um2cable_rr( um1%longitude,rad_cbl%longitude )

END SUBROUTINE initialize_radiation

!========================================================================
!========================================================================
!========================================================================
          
SUBROUTINE initialize_canopy(canopy_tile,visc_sublayer_dz, canopy_cbl)
USE cable_um_tech_mod,   ONLY: um1 
USE cable_common_module, ONLY: cable_runtime, cable_user
USE cable_canopy_type_mod, ONLY: canopy_type
TYPE(canopy_type),     INTENT(INOUT)  :: canopy_cbl
   
REAL, INTENT(IN),DIMENSION(um1%land_pts, um1%ntiles) :: canopy_tile
REAL, INTENT(IN),DIMENSION(um1%land_pts, um1%ntiles) :: visc_sublayer_dz
   
! defs 1st call to CABLE in this run. OK in UM & coupled
LOGICAL, SAVE :: first_call= .TRUE.
   
   !--- %ga is computed (on LHS only) in define_canopy and 
   !--- then used in soilsnow() in implicit call, then unpacked
IF ( first_call ) THEN
  canopy_cbl%ga = 0.0
  canopy_cbl%us = 0.1 !to match Loobos
  canopy_cbl%fes_cor = 0.0
  canopy_cbl%fhs_cor = 0.0
  canopy_cbl%fwsoil = 1.0
  first_call = .FALSE.
END IF

!---set canopy storage (already in dim(land_pts,ntiles) ) 
canopy_cbl%cansto = PACK(canopy_tile, um1%l_tile_pts)
canopy_cbl%oldcansto = canopy_cbl%cansto
canopy_cbl%sublayer_dz(:) = 0.0  !junk
canopy_cbl%sublayer_dz(:) = PACK(visc_sublayer_dz(:,:),um1%l_tile_pts) 
!H!this is not used in HAC, not unpacked anyway AND why it didnt show up
!aserror until after mods to outpput.nml I'll never know
!H!where (canopy%sublayer_dz .lt. 1.0e-8) canopy%sublayer_dz = 1.0e-8
!H!where (canopy%sublayer_dz .gt. 1.0) canopy%sublayer_dz = 1.0

IF (first_call ) THEN

  WRITE(6,*) 'maxval canopy%sublayer_dz',MAXVAL(canopy_cbl%sublayer_dz,dim = 1)
  WRITE(6,*) 'minval canopy%sublayer_dz',MINVAL(canopy_cbl%sublayer_dz,dim = 1)

END IF

END SUBROUTINE initialize_canopy

!========================================================================
!========================================================================
!========================================================================
 
SUBROUTINE initialize_soilsnow( smvcst, tsoil_tile, sthf_tile,smcl_tile,smgw_tile, &
                                snow_tile, snow_rho1l, snow_age, isnow_flg3l, &
                                snow_rho3l, snow_cond, snow_depth3l,          &
                                snow_mass3l, snow_tmp3l, fland,               &
                                sin_theta_latitude, soil_cbl, ssnow_cbl, met_cbl, bal_cbl, veg_cbl )

USE cable_def_types_mod,  ONLY: mp, msn, ms, r_2,mstype
USE cable_um_tech_mod,   ONLY: um1
USE cable_common_module, ONLY: cable_runtime, cable_user, ktau_gl
   
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_met_type_mod, ONLY: met_type
USE cable_balances_type_mod, ONLY: balances_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soil_parameter_type

TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(soil_parameter_type),   INTENT(INOUT) :: soil_cbl
TYPE(soil_snow_type),     INTENT(INOUT)  :: ssnow_cbl
TYPE(met_type),     INTENT(INOUT)  :: met_cbl
TYPE(balances_type),     INTENT(INOUT)  :: bal_cbl


REAL, INTENT(IN), DIMENSION(um1%land_pts) :: smvcst
   
REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles, um1%sm_levels) ::       &
   sthf_tile, &   !
   smcl_tile, &   !
   tsoil_tile     !

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) ::                      &
   smgw_tile

INTEGER, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) :: isnow_flg3l 

REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles) :: snow_tile

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) ::                      &
   snow_rho1l, &  !
   snow_age     !

REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles,3) :: snow_cond

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles,3) ::                    &
   snow_rho3l,    & !
   snow_depth3l,  & !
   snow_mass3l,   & !
   snow_tmp3l       !
   
REAL, INTENT(IN), DIMENSION(um1%land_pts) :: fland 
   
REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: sin_theta_latitude
   
INTEGER :: i,j,k,l,n
REAL  :: zsetot, max_snow_depth = 50000.0
REAL, ALLOCATABLE:: fwork(:,:,:), sfact(:), fvar(:), rtemp(:),                &
                    tot_mass_tmp(:,:,:), ice_vol_tmp(:,:,:)
LOGICAL :: skip =.TRUE. 
LOGICAL, SAVE :: first_call = .TRUE.
REAL, DIMENSION(mstype) :: dummy 
      
dummy = 0.0 

!     not sure if this is in restart file hence repeated again
!H!IF ( first_call) THEN 
!H!  ssnow_cbl%pudsto = 0.0 
!H!END IF
!H!ssnow_cbl%pudsmx = 0.0
!H!ssnow_cbl%wbtot1 = 0
!H!ssnow_cbl%wbtot2 = 0
!H!ssnow_cbl%wb_lake = 0.0

 !why is snow updated from um values every timestep
 !but soil moisture not?
      
!line removed Jun 2018 to assist with water conservation
!in coupled runs alongside the ice-berg scheme
!snow_tile = MIN(max_snow_depth, snow_tile)

ssnow_cbl%snowd  = PACK(snow_tile,um1%l_tile_pts)
ssnow_cbl%ssdnn  = PACK(snow_rho1l,um1%l_tile_pts)  
ssnow_cbl%isflag = PACK(INT(isnow_flg3l),um1%l_tile_pts)  
      
DO j = 1, msn
         
  ssnow_cbl%sdepth(:,j)= PACK(snow_depth3l(:,:,j),um1%l_tile_pts)
  ssnow_cbl%smass(:,j) = PACK(snow_mass3l(:,:,j),um1%l_tile_pts)  
  ssnow_cbl%ssdn(:,j)  = PACK(snow_rho3l(:,:,j),um1%l_tile_pts)  
  ssnow_cbl%tggsn(:,j) = PACK(snow_tmp3l(:,:,j),um1%l_tile_pts)  
  ssnow_cbl%sconds(:,j)= PACK(snow_cond(:,:,j),um1%l_tile_pts)  
         
END DO 
!ssnow%wb_lake = MAX( ssnow%wbtot2 - ssnow%wbtot1, 0.)
       
DO j = 1,um1%sm_levels
  ssnow_cbl%tgg(:,j) = PACK(tsoil_tile(:,:,j),um1%l_tile_pts)
END DO 
ssnow_cbl%snage = PACK(snow_age, um1%l_tile_pts)

ssnow_cbl%GWwb(:) = PACK(smgw_tile(:,:),um1%l_tile_pts)
!H!this is not used in HAC, not unpacked anyway AND why it didnt show up
 !H!where (ssnow%GWwb .gt. soil%GWssat_vec) ssnow%GWwb = soil%GWssat_vec

IF ( first_call) THEN 
        
  ssnow_cbl%wbtot = 0.0
  ssnow_cbl%wb_lake = 0.0
  ssnow_cbl%totwblake = 0.0  ! wb_lake integrated over river timestep
  ssnow_cbl%tggav = 0.0
  ssnow_cbl%rtsoil = 50.0
  ssnow_cbl%rtsoil = 100.0 !to match offline for Loobos
  ssnow_cbl%t_snwlr = 0.05

  ! snow depth from prev timestep 
  ssnow_cbl%osnowd  = PACK(snow_tile,um1%l_tile_pts)  

  zsetot = SUM(soil_cbl%zse)
  DO k = 1, um1%sm_levels
    ssnow_cbl%tggav = ssnow_cbl%tggav  + soil_cbl%zse(k) * ssnow_cbl%tgg(:,k) / zsetot
  END DO
     
  ! not updated 
  ALLOCATE( sfact( mp ) )
  sfact = 0.68
  WHERE (soil_cbl%albsoil(:,1) <= 0.14) 
    sfact = 0.5
  ELSEWHERE (soil_cbl%albsoil(:,1) > 0.14 .AND. soil_cbl%albsoil(:,1) <= 0.20)
    sfact = 0.62
  END WHERE
  ssnow_cbl%albsoilsn(:,2) = 2.0 * soil_cbl%albsoil(:,1) / (1.0 + sfact)
  ssnow_cbl%albsoilsn(:,1) = sfact * ssnow_cbl%albsoilsn(:,2)
  DEALLOCATE( sfact )

  ALLOCATE( fvar(um1%land_pts ) )
  DO n = 1,um1%ntiles
    DO k = 1,um1%tile_pts(n)
      l = um1%tile_index(k,n)
      fvar(l) = REAL(l)
    END DO
  END DO
  CALL um2cable_lp( fland,dummy , ssnow_cbl%fland, soil_cbl%isoilm, skip )
  CALL um2cable_lp( fvar, dummy, ssnow_cbl%ifland, soil_cbl%isoilm, skip )
  DEALLOCATE( fvar )
         
  !--- updated via smcl,sthf etc 
  !previous code required rho_water == rho_ice
  !GWmodel correctly uses rho_ice ~= 0.92*rho_water
  !code now handles when densities are ice==water and ice!=water
  IF (ALLOCATED(tot_mass_tmp)) DEALLOCATE(tot_mass_tmp)
  IF (ALLOCATED(ice_vol_tmp)) DEALLOCATE(ice_vol_tmp)

  ALLOCATE( ice_vol_tmp(um1%land_pts,um1%ntiles,um1%sm_levels) )
  ALLOCATE( tot_mass_tmp(um1%land_pts,um1%ntiles,um1%sm_levels) )

  ice_vol_tmp(:,:,:) = 0.0
  tot_mass_tmp(:,:,:) = 0.0

  DO n = 1,um1%ntiles                                                   
    DO k = 1,um1%tile_pts(n)                                           
      i = um1%tile_index(k,n)                                      
      DO j = 1,um1%sm_levels
        tot_mass_tmp(i,n,j) = smcl_tile(i,n,j)
        ice_vol_tmp(i,n,j) = sthf_tile(i,n,j) !matching Loobos * smvcst(i)
      END DO ! J
    END DO
  END DO
   
  DO j = 1,um1%sm_levels
     !ice volume
    ssnow_cbl%wb(:,j) =     pack(tot_mass_tmp(:,:,J),um1%l_tile_pts) !match offline 
    ssnow_cbl%wbice(:,j) = PACK(ice_vol_tmp(:,:,j),um1%l_tile_pts)
    !ssnow_cbl%wbice(:,j) = MAX(0.0,ssnow_cbl%wbice(:,j))  !should not be needed -- mrd561
    !liq volume  from (tot_mass - ice_mass) / (dz*rho_liq)
    ssnow_cbl%wbliq(:,j)= ( ssnow_cbl%wb(:,j) -  ssnow_cbl%wbice(:,j) )  / (soil_cbl%zse(j) * um1%rho_water)   
!jhan: both already per layer WHY is per layer needed??    
    ssnow_cbl%wbliq(:,j)= ( ssnow_cbl%wb(:,j) -  ssnow_cbl%wbice(:,j) )  * um1%rho_water   
    !ssnow_cbl%wbliq(:,j)= (PACK(tot_mass_tmp(:,:,j),um1%l_tile_pts) - & !total mass
    !                   ssnow_cbl%wbice(:,j) * soil_cbl%zse(j) * um1%rho_ice & !subtract ice mass
    !                   ) / (soil_cbl%zse(j) * um1%rho_water)            !convert units
    !ssnow_cbl%wb(:,j)    = ssnow_cbl%wbice(:,j) + ssnow_cbl%wbliq(:,j) 
!!match offline ssnow_cbl%wb(:,j)= pack(tot_mass_tmp(:,:,J),um1%l_tile_pts) !match offline 
    !WHERE( veg_cbliveg == 16 ) ssnow%wb(:,J) = 0.95*soil%ssat
    !WHERE( veg_cbliveg == 16 ) ssnow%wb(:,J) = soil%sfc
  END DO
!match offline ssnow_cbl%wbliq = ssnow_cbl%wb - ssnow_cbl%wbice !match offline 
         
  DEALLOCATE( tot_mass_tmp )
  DEALLOCATE( ice_vol_tmp )

         
  ssnow_cbl%owetfac = MAX( 0.0, MIN( 1.0,                                     &
                  ( ssnow_cbl%wb(:,1) - soil_cbl%swilt ) /                    &
                  ( MAX(0.083, (soil_cbl%sfc - soil_cbl%swilt) ) )            &
                  ) )

  ! Temporay fix for accounting for reduction of soil evaporation 
  ! due to freezing
  WHERE ( ssnow_cbl%wbice(:,1) > 0.0 )                                        &
     ! Prevents divide by zero at glaciated points where both 
     ! wb and wbice=0.
     ssnow_cbl%owetfac = ssnow_cbl%owetfac * ( 1.0 - ssnow_cbl%wbice(:,1) /   &
                     ssnow_cbl%wb(:,1) )**2
      
  !jhan: do we want to do this before %owetfac is set 
  DO j = 1, um1%sm_levels
     !should be removed!!!!!!!! This cannot conserve if there are any
     !dynamics 
    WHERE ( soil_cbl%isoilm == 9 ) ! permanent ice: remove hard-wired no. in future
      ssnow_cbl%wb(:,j) = 0.95 * soil_cbl%ssat
      ssnow_cbl%wbice(:,j) = 0.85 * ssnow_cbl%wb(:,j)
    END WHERE
    !no not force rho_water==rho_ice==1000.0
    ssnow_cbl%wbtot = ssnow_cbl%wbtot + soil_cbl%zse(j) *                     &
                               (ssnow_cbl%wbliq(:,j) * um1%rho_water+         &
                                ssnow_cbl%wbice(:,j) * um1%rho_ice )                     
  END DO
  IF (cable_user%gw_model) THEN
    ssnow_cbl%wbtot = ssnow_cbl%wbtot + ssnow_cbl%GWwb(:) * soil_cbl%GWdz * um1%rho_water
  END IF
     
  bal_cbl%wbtot0 = ssnow_cbl%wbtot

! !H!   !---set antartic flag using  sin_theta_latitude(row_length,rows)
! !H!   ALLOCATE( fwork(1,um1%land_pts,um1%ntiles) )
! !H!   fwork = 0.0
! !H!   DO n = 1,um1%ntiles                     
! !H!     DO k = 1,um1%tile_pts(n)
! !H!       l = um1%tile_index(k,n)
! !H!       j=(um1%land_index(l) - 1) / um1%row_length + 1
! !H!       i = um1%land_index(l) - (j-1) * um1%row_length
! !H!       IF ( sin_theta_latitude(i,j)  <   -0.91 ) fwork(1,l,n) = 1.0
! !H!     END DO
! !H!   END DO
! !H!   ssnow_cbl%iantrct = PACK(fwork(1,:,:),um1%l_tile_pts)
! !H!         
! !H!   DEALLOCATE( fwork )

  !! SLI specific initialisations:
  !IF(cable_user%SOIL_STRUC=='sli') THEN
  !   ssnow%h0(:)        = 0.0
  !   ssnow%S(:,:)       = ssnow%wb(:,:)/SPREAD(soil%ssat,2,ms)
  !   ssnow%snowliq(:,:) = 0.0
  !   ssnow%Tsurface     = 25.0
  !   ssnow%nsnow        = 0
  !   !ssnow%Tsoil        = ssnow%tgg - 273.16
  !   ssnow%kth          = 0.3
  !   ssnow%lE           = 0.
  !   ! vh ! should be calculated from soil moisture or be in restart file
  !   !ssnow%sconds(:,:)  = 0.06_r_2    ! vh snow thermal cond (W m-2 K-1),
  !   ! should be in restart file
  !END IF

  first_call = .FALSE.

END IF ! END: if (first_call)       

!mrd561
!should be initialized but not sure if here
ssnow_cbl%qrecharge    = 0.0
ssnow_cbl%wtd          = 1.0
ssnow_cbl%rtevap_sat   = 0.0
ssnow_cbl%rtevap_unsat = 0.0
ssnow_cbl%satfrac      = 0.5
ssnow_cbl%qhz          = 0.0
ssnow_cbl%qhlev        = 0.0
ssnow_cbl%wbliq = ssnow_cbl%wb - ssnow_cbl%wbice

!H!! SLI specific initialisations:
!H!IF (cable_user%soil_struc == 'sli') THEN
!H!  ssnow_cbl%h0(:)        = 0.0
!H!  ssnow_cbl%s(:,:)       = ssnow_cbl%wb(:,:) /soil_cbl%ssat_vec !SPREAD(soil%ssat,2,ms)
!H!  ssnow_cbl%snowliq(:,:) = 0.0
!H!  ssnow_cbl%Tsurface     = 25.0  !why not ssnow%tgg(:,1) - 273.16
!H!  ssnow_cbl%nsnow        = 0
!H!  ssnow_cbl%Tsoil        = ssnow_cbl%tgg - 273.16
!H!  ssnow_cbl%kth          = 0.3
!H!  ssnow_cbl%lE           = 0.0
!H!  ! vh ! should be calculated from soil moisture or be in restart file
!H!  ssnow_cbl%sconds(:,:)  = 0.06_r_2    ! vh snow thermal cond (W m-2 K-1),
!H!  ! should be in restart file
!H!END IF
!H!
!H!
!H!IF (cable_user%gw_model) THEN
!H!  ssnow_cbl%wb_lake(:) = 0.0  !already prevent drainage unless fully
!H!                          !saturated and can take from GWwb
!H!  ssnow_cbl%wbtot1 = 0.0
!H!  ssnow_cbl%wbtot2 = 0.0
!H!  ssnow_cbl%wbtot1 = MAX(0.0,                                                 &
!H!                     (soil_cbl%sfc_vec(:,1) - ssnow_cbl%wb(:,1)) * soil_cbl%zse(1) / soil_cbl%GWdz(:) )
!H!
!H!  WHERE ( veg_cbl%iveg == 16 .AND.                                            &
!H!         ssnow_cbl%wb(:,1) < soil_cbl%sfc_vec(:,1) .AND.                      &
!H!         ssnow_cbl%GWwb(:)  >   ssnow_cbl%wbtot1)
!H!
!H!    ssnow_cbl%wb(:,1) = soil_cbl%sfc_vec(:,1)
!H!    ssnow_cbl%GWwb(:) = ssnow_cbl%GWwb(:) -  ssnow_cbl%wbtot1(:)
!H!
!H!  END WHERE
!H!
!H!  ssnow_cbl%wbtot1 = 0.0
!H!          
!H!ELSE
!H!  !     DO J=1, msn
!H!  ssnow_cbl%wbtot1 = 0.0
!H!  ssnow_cbl%wbtot2 = 0.0
!H!  DO j = 1, 1
!H!
!H!    WHERE ( veg_cbl%iveg == 16 .AND. ssnow_cbl%wb(:,j) < soil_cbl%sfc ) 
!H!        ! lakes: remove hard-wired number in future version
!H!      ssnow_cbl%wbtot1 = ssnow_cbl%wbtot1 + REAL( ssnow_cbl%wb(:,j) ) * 1000.0 * &
!H!                     soil_cbl%zse(j)
!H!      ssnow_cbl%wb(:,j) = soil_cbl%sfc
!H!      ssnow_cbl%wbtot2 = ssnow_cbl%wbtot2 + REAL( ssnow_cbl%wb(:,j) ) * 1000.0 * &
!H!                     soil_cbl%zse(j)
!H!    END WHERE
!H!
!H!  END DO
!H!  ssnow_cbl%wb_lake = MAX( ssnow_cbl%wbtot2 - ssnow_cbl%wbtot1, 0.0)
!H!END IF 

END SUBROUTINE initialize_soilsnow
 
!========================================================================
!========================================================================
!========================================================================
          
SUBROUTINE initialize_roughness( z1_tq, z1_uv, htveg, rough_cbl, veg_cbl )  
USE cable_um_tech_mod,   ONLY: um1
USE cable_common_module, ONLY: ktau_gl
USE cable_def_types_mod, ONLY: mp
USE cable_common_module, ONLY: cable_runtime, cable_user
USE cable_roughness_type_mod, ONLY: roughness_type

USE cable_params_mod,         ONLY: veg_parameter_type

TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(roughness_type),     INTENT(INOUT)  :: rough_cbl

   
REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::  z1_tq, z1_uv
REAL, DIMENSION(um1%land_pts, um1%ntiles) :: htveg
INTEGER :: i,j,k,l,n
REAL, ALLOCATABLE, DIMENSION(:,:) :: jhruff, jhwork

   !--- CABLE roughness type forcings
CALL um2cable_rr( z1_tq, rough_cbl%za_tq)
CALL um2cable_rr( z1_uv, rough_cbl%za_uv)

END SUBROUTINE initialize_roughness

!--- UM met forcing vars needed by CABLE commonly have UM dimensions
!---(row_length,rows), which is no good to CABLE. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE um2cable_rr(umvar,cablevar)
USE cable_def_types_mod, ONLY: mp
USE cable_um_tech_mod,   ONLY:um1
 
REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: umvar   
REAL, INTENT(INOUT), DIMENSION(mp) :: cablevar
REAL, DIMENSION(um1%land_pts,um1%ntiles) :: fvar   
INTEGER :: n,k,l,j,i

fvar = 0.0
DO n = 1,um1%ntiles                     
  DO k = 1,um1%tile_pts(n)
    l = um1%tile_index(k,n)
    j=(um1%land_index(l) - 1) / um1%row_length + 1
    i = um1%land_index(l) - (j-1) * um1%row_length
    fvar(l,n) = umvar(i,j)
  END DO
END DO
cablevar =  PACK(fvar,um1%l_tile_pts)

END SUBROUTINE um2cable_rr

!========================================================================

SUBROUTINE um2cable_irr(umvar,cablevar)
USE cable_def_types_mod, ONLY: mp
USE cable_um_tech_mod,   ONLY:um1
 
INTEGER, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: umvar   
INTEGER, INTENT(INOUT), DIMENSION(mp) :: cablevar
INTEGER, DIMENSION(um1%land_pts,um1%ntiles) :: fvar   
INTEGER :: n,k,l,j,i

fvar = 0.0
DO n = 1,um1%ntiles                     
  DO k = 1,um1%tile_pts(n)
    l = um1%tile_index(k,n)
    j=(um1%land_index(l) - 1) / um1%row_length + 1
    i = um1%land_index(l) - (j-1) * um1%row_length
    fvar(l,n) = umvar(i,j)
  END DO
END DO
cablevar =  PACK(fvar,um1%l_tile_pts)

END SUBROUTINE um2cable_irr


!========================================================================
!========================================================================

!--- UM met forcing vars needed by CABLE which have UM dimensions
!---(land_points)[_lp], which is no good to cable. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE um2cable_lp(umvar, defaultin, cablevar, soiltype, skip )
USE cable_def_types_mod, ONLY: mp, mstype
USE cable_um_tech_mod,   ONLY:um1
  
REAL, INTENT(IN), DIMENSION(um1%land_pts) :: umvar
REAL, INTENT(IN), DIMENSION(mstype) :: defaultin    
REAL, INTENT(INOUT), DIMENSION(mp) :: cablevar
INTEGER, INTENT(INOUT), DIMENSION(mp) :: soiltype
REAL, DIMENSION(:,:), ALLOCATABLE:: fvar   
LOGICAL, OPTIONAL :: skip
INTEGER :: n,k,l,i

         
ALLOCATE( fvar(um1%land_pts,um1%ntiles) )
fvar = 0.0

! loop over Ntiles
DO n = 1,um1%ntiles
   ! loop over number of points per tile
  DO k = 1,um1%tile_pts(n)
     ! index of each point per tile in an array of dim=(land_pts,ntiles)
    l = um1%tile_index(k,n)
    ! at this point fvar=umvar, ELSE=0.0 
    fvar(l,n) = umvar(l)
    ! unless explicitly SKIPPED by including arg in subr call
    IF ( .NOT. PRESENT(skip) ) THEN
       ! on perma frost tile, set fvar=defaultin
      IF ( n == um1%ntiles ) THEN
        fvar(l,n) =  defaultin(9)
      END IF
    END IF
  END DO
END DO
     
cablevar     =  PACK(fvar,um1%l_tile_pts)
  
! unless explicitly SKIPPED by including arg in subr call
IF ( .NOT. PRESENT(skip) ) THEN
  DO i = 1,mp
     ! soiltype=9 for perma-frost tiles 
    IF (soiltype(i) == 9) cablevar(i) =  defaultin(9)         
  END DO        
END IF
   
DEALLOCATE(fvar)

END SUBROUTINE um2cable_lp
 

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_bgc_vars(bgc_cbl, veg_cbl) 
USE cable_def_types_mod, ONLY: ncs, ncp 
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
IMPLICIT NONE

TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(bgc_pool_type),       INTENT(INOUT)  :: bgc_cbl

   
INTEGER :: k

! note that ratecp and ratecs are the same for all veg at the moment. (BP)
DO k = 1,ncp
  bgc_cbl%cplant(:,k) = veg_cbl%cplant(:,k)
  bgc_cbl%ratecp(k) = veg_cbl%ratecp(1,k)
END DO
DO k = 1,ncs
  bgc_cbl%csoil(:,k) = veg_cbl%csoil(:,k)
  bgc_cbl%ratecs(k) = veg_cbl%ratecs(1,k)
END DO
END SUBROUTINE init_bgc_vars

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_sumflux_zero(sum_flux_cbl) 
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
TYPE(sum_flux_type),  INTENT(INOUT)  :: sum_flux_cbl
sum_flux_cbl%sumpn = 0.0; sum_flux_cbl%sumrp = 0.0; sum_flux_cbl%sumrpw = 0.0
sum_flux_cbl%sumrpr = 0.0; sum_flux_cbl%sumrs = 0.0; sum_flux_cbl%sumrd = 0.0
sum_flux_cbl%dsumpn = 0.0; sum_flux_cbl%dsumrp = 0.0; sum_flux_cbl%dsumrs = 0.0
sum_flux_cbl%dsumrd = 0.0; sum_flux_cbl%sumxrp = 0.0;  sum_flux_cbl%sumxrs = 0.0
END SUBROUTINE init_sumflux_zero 

END MODULE cable_um_init_subrs_mod





