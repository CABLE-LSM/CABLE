!This is ONLY here as CABLE sciencecode has not yet been upgraded and is USEing cable_define_types everywhere
!little does it know that it is actually USE these new *_cble versions
MODULE cable_def_types_mod

USE cable_types_mod!!,          ONLY: mp, l_tile_pts
USE cable_air_type_mod,       ONLY: air_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_climate_type_mod,   ONLY: climate_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soil_parameter_type

IMPLICIT NONE

END MODULE cable_def_types_mod

