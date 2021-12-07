MODULE cbl_ssnow_data_mod

USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
       veg_parameter_type, canopy_type, met_type,        &
       balances_type, r_2, ms, mp
!distribute these per sbr
USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ
USE cable_phys_constants_mod, ONLY : CHL => HL
USE cable_phys_constants_mod, ONLY : CHLF => HLF
USE cable_phys_constants_mod, ONLY : Cdensity_liq => density_liq
USE cable_phys_constants_mod, ONLY : Ccgsnow => cgsnow
USE cable_phys_constants_mod, ONLY : Ccswat => cswat
USE cable_phys_constants_mod, ONLY : Ccsice => csice
USE cable_phys_constants_mod, ONLY : Ccs_rho_wat => cs_rho_wat
USE cable_phys_constants_mod, ONLY : Ccs_rho_ice => cs_rho_ice

USE cable_common_module, ONLY: cable_user,snow_ccnsw,snmin,&
       max_ssdn,max_sconds,frozen_limit,&
       max_glacier_snowd

  IMPLICIT NONE

END MODULE cbl_ssnow_data_mod
