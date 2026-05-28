MODULE casavariable

!legacy - USE all the other MODULES here
! data declarations
USE casa_biome_type_mod,    ONLY: casa_biome_data_type,      &
                                  casa_biome => casa_biome_type
USE casa_pool_type_mod,     ONLY: casa_pool_data_type,       &
                                  casa_pool  => casa_pool_type
USE casa_flux_type_mod,     ONLY: casa_flux_data_type,       &
                                  casa_flux  => casa_flux_type
USE casa_met_type_mod,      ONLY: casa_met_data_type,        &
                                  casa_met   =>  casa_met_type
USE casa_balance_type_mod,  ONLY: casa_bal_data_type,        &
                                  casa_bal   => casa_bal_type             
USE phenology_type_mod,     ONLY: phenology_data_type,       &
                                  phenology_type
IMPLICIT NONE

PUBLIC :: zero_sum_casa
PUBLIC :: update_sum_casa
PUBLIC :: alloc_casavariable

CONTAINS

SUBROUTINE alloc_casavariable( casabiome, casabiome_data,                      &
                               casapool,  casapool_data,                       &
                               casaflux,  casaflux_data,                       &
                               casamet,   casamet_data,                        &
                               casabal,   casabal_data,                        &
                               sum_casapool, sum_casapool_data,                &
                               sum_casaflux, sum_casaflux_data,                &
                               phen, phen_data, mp )

! subrs                               
USE casa_biome_type_mod, ONLY: zero_casa_biome_data_type
USE casa_biome_type_mod, ONLY: assoc_casa_biome_type

USE casa_pool_type_mod, ONLY: alloc_casa_pool_data_type
USE casa_pool_type_mod, ONLY: assoc_casa_pool_type
USE casa_pool_type_mod, ONLY: zero_casa_pool_data_type

USE casa_flux_type_mod, ONLY: alloc_casa_flux_data_type
USE casa_flux_type_mod, ONLY: zero_casa_flux_data_type
USE casa_flux_type_mod, ONLY: assoc_casa_flux_type

USE casa_met_type_mod, ONLY: alloc_casa_met_data_type
USE casa_met_type_mod, ONLY: assoc_casa_met_type
USE casa_met_type_mod, ONLY: zero_casa_met_data_type

USE sum_casa_pool_type_mod, ONLY: alloc_sum_casa_pool_data_type
USE sum_casa_pool_type_mod, ONLY: assoc_sum_casa_pool_type

USE sum_casa_flux_type_mod, ONLY: alloc_sum_casa_flux_data_type
USE sum_casa_flux_type_mod, ONLY: assoc_sum_casa_flux_type

USE casa_balance_type_mod, ONLY: alloc_casa_bal_data_type  
USE casa_balance_type_mod, ONLY: zero_casa_bal_data_type   
USE casa_balance_type_mod, ONLY: assoc_casa_bal_type       

USE phenology_type_mod,    ONLY: alloc_phenology_data_type
USE phenology_type_mod,    ONLY: zero_phenology_data_type
USE phenology_type_mod,    ONLY: assoc_phenology_type

IMPLICIT NONE

TYPE (casa_biome),           INTENT(INOUT) :: casabiome
TYPE (casa_biome_data_type), INTENT(INOUT) :: casabiome_data
TYPE (casa_pool),            INTENT(INOUT) :: casapool
TYPE (casa_pool_data_type),  INTENT(INOUT) :: casapool_data
TYPE (casa_flux),            INTENT(INOUT) :: casaflux
TYPE (casa_flux_data_type),  INTENT(INOUT) :: casaflux_data
TYPE (casa_pool),            INTENT(INOUT) :: sum_casapool
TYPE (casa_pool_data_type),  INTENT(INOUT) :: sum_casapool_data
TYPE (casa_flux),            INTENT(INOUT) :: sum_casaflux
TYPE (casa_flux_data_type),  INTENT(INOUT) :: sum_casaflux_data
TYPE (casa_bal),             INTENT(INOUT) :: casabal 
TYPE (casa_bal_data_type),   INTENT(INOUT) :: casabal_data 
TYPE (casa_met),             INTENT(INOUT) :: casamet 
TYPE (casa_met_data_type),   INTENT(INOUT) :: casamet_data 
TYPE (phenology_type),       INTENT(INOUT) :: phen 
TYPE (phenology_data_type),  INTENT(INOUT) :: phen_data 

INTEGER,                     INTENT(INOUT) :: mp

CALL zero_casa_biome_data_type( casabiome_data )   
! allocating as 1s can also be done here following cable types, given lsm_id=casa
CALL alloc_casa_pool_data_type( casapool_data, mp )   
CALL zero_casa_pool_data_type( casapool_data )   

CALL alloc_casa_flux_data_type( casaflux_data, mp )   
CALL zero_casa_flux_data_type( casaflux_data ) 

CALL alloc_casa_met_data_type( casamet_data, mp )   
CALL zero_casa_met_data_type( casamet_data )   

CALL alloc_sum_casa_pool_data_type( sum_casapool_data, mp )   

CALL alloc_sum_casa_flux_data_type( sum_casaflux_data, mp )   

CALL alloc_casa_bal_data_type( casabal_data, mp )   
CALL zero_casa_bal_data_type( casabal_data )      

CALL alloc_phenology_data_type( phen_data, mp )   
CALL zero_phenology_data_type( phen_data )      

!IF ( cable_runtime%um ) THEN 
  CALL assoc_casa_biome_type( casabiome, casabiome_data )      
  CALL assoc_casa_pool_type  ( casapool,  casapool_data  )      
  CALL assoc_casa_flux_type ( casaflux,  casaflux_data  )      
  CALL assoc_casa_met_type  ( casamet,   casamet_data   )      
  CALL assoc_casa_bal_type  ( casabal,   casabal_data   )             
  CALL assoc_sum_casa_pool_type( sum_casapool, sum_casapool_data )
  CALL assoc_sum_casa_flux_type( sum_casaflux, sum_casaflux_data )
  CALL assoc_phenology_type( phen, phen_data )
!END IF

RETURN                               
END SUBROUTINE alloc_casavariable

SUBROUTINE zero_sum_casa( sum_casapool, sum_casaflux )

USE sum_casa_pool_type_mod, ONLY: zero_sum_casa_pool_data
USE sum_casa_flux_type_mod, ONLY: zero_sum_casa_flux_data

TYPE (casa_pool_data_type), INTENT(INOUT) :: sum_casapool
TYPE (casa_flux_data_type), INTENT(INOUT) :: sum_casaflux

CALL zero_sum_casa_pool_data( sum_casapool ) 
CALL zero_sum_casa_flux_data( sum_casaflux ) 

RETURN
END SUBROUTINE zero_sum_casa

SUBROUTINE update_sum_casa( sum_casapool, sum_casaflux, casapool, casaflux,    &
                            sum_now, average_now, nsteps )

USE sum_casa_pool_type_mod, ONLY: update_sum_casa_pool
USE sum_casa_flux_type_mod, ONLY: update_sum_casa_flux

IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: sum_casapool
TYPE (casa_flux_data_type), INTENT(INOUT) :: sum_casaflux
TYPE (casa_pool_data_type), INTENT(IN)    :: casapool
TYPE (casa_flux_data_type), INTENT(IN)    :: casaflux
LOGICAL,               INTENT(IN)    :: sum_now, average_now
INTEGER,               INTENT(IN)    :: nsteps

CALL update_sum_casa_pool( sum_casapool, casapool, sum_now, average_now, nsteps)
                            
CALL update_sum_casa_flux( sum_casapool, sum_casaflux, casapool, casaflux,     &
                           sum_now, average_now, nsteps )

RETURN
END SUBROUTINE update_sum_casa

END MODULE casavariable
