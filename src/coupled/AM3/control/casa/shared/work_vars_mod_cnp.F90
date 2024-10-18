MODULE work_vars_mod_cnp

!!USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
!!USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_data_type
!!USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
!!USE cable_sum_flux_type_mod,  ONLY: sum_flux_data_type

USE casa_biome_type_mod,      ONLY: casa_biome_data_type
USE casa_biome_type_mod,      ONLY: casa_biome_type
USE casa_pool_type_mod,       ONLY: casa_pool_data_type
USE casa_pool_type_mod,       ONLY: casa_pool_type
USE casa_flux_type_mod,       ONLY: casa_flux_data_type
USE casa_flux_type_mod,       ONLY: casa_flux_type
USE casa_met_type_mod,        ONLY: casa_met_data_type
USE casa_met_type_mod,        ONLY: casa_met_type
USE casa_balance_type_mod,    ONLY: casa_bal_data_type        
USE casa_balance_type_mod,    ONLY: casa_bal_type             
USE phenology_type_mod,       ONLY: phenology_data_type
USE phenology_type_mod,       ONLY: phenology_type

IMPLICIT NONE

PUBLIC :: work_vars_type
PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WORK_VARS_MOD_CNP'

TYPE :: work_vars_type

  TYPE( casa_biome_type ) :: casabiome
  TYPE( casa_pool_type  ) :: casapool
  TYPE( casa_pool_type  ) :: sum_casapool
  TYPE( casa_flux_type  ) :: casaflux
  TYPE( casa_flux_type  ) :: sum_casaflux
  TYPE( casa_met_type   ) :: casamet
  TYPE( casa_bal_type   ) :: casabal
  TYPE( phenology_type )  :: phen 

END TYPE work_vars_type

! data arrays need to be declared outside of the work% TYPE
TYPE( casa_biome_data_type ), PUBLIC, TARGET :: casabiome_data              
TYPE( casa_pool_data_type  ), PUBLIC, TARGET :: casapool_data
TYPE( casa_pool_data_type  ), PUBLIC, TARGET :: sum_casapool_data
TYPE( casa_flux_data_type  ), PUBLIC, TARGET :: casaflux_data
TYPE( casa_flux_data_type  ), PUBLIC, TARGET :: sum_casaflux_data
TYPE( casa_met_data_type   ), PUBLIC, TARGET :: casamet_data
TYPE( casa_bal_data_type   ), PUBLIC, TARGET :: casabal_data
TYPE( phenology_data_type ),  PUBLIC, TARGET :: phen_data 

END MODULE work_vars_mod_cnp

