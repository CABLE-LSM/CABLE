MODULE cable_hyd_main_mod
  
CONTAINS

SUBROUTINE cable_hyd_main( land_pts, ntiles, lying_snow, SNOW_surft, surf_roff, &
                           sub_surf_roff, tot_tfall,                          &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soil_cbl )
  
  !subrs called 
USE cable_hyd_driv_mod, ONLY: cable_hyd_driver
  
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
  !processor number, timestep number / width, endstep
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY: cable_runtime

!data
USE cable_types_mod, ONLY: L_tile_pts
 
IMPLICIT NONE
 
!___ re-decl input args

INTEGER :: land_pts, ntiles
  
REAL :: snow_surft(land_pts,ntiles)
  
REAL, DIMENSION(land_pts) ::                                                  &
  lying_snow,    & ! OUT Gridbox snowmass (kg/m2)        
  sub_surf_roff,  & ! OUT Sub-surface runoff (kg/m2/s).
  surf_roff,      & ! OUT Surface runoff (kg/m2/s).
  tot_tfall         ! OUT Total throughfall (kg/m2/s).

!___ local vars
LOGICAL, SAVE :: first_call = .TRUE.
  
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
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_hyd_main"
 
!H!# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
  
  !--- initialize cable_runtime% switches 
cable_runtime%um =          .false.
cable_runtime%um_hydrology =.TRUE.
  
CALL cable_hyd_driver( land_pts, ntiles, L_tile_pts, lying_snow, SNOW_surft,  &
                       surf_roff, sub_surf_roff, tot_tfall, ssnow_cbl, canopy_cbl, veg_cbl )
  
cable_runtime%um_hydrology =.FALSE.
  
!-------- End Unique subroutine body -----------
  
first_call = .FALSE.        

RETURN

END SUBROUTINE cable_hyd_main
  
END MODULE cable_hyd_main_mod

