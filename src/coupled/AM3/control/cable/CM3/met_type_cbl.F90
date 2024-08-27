MODULE cable_met_type_mod

IMPLICIT NONE

PUBLIC :: met_type
PUBLIC :: met_data_type
PUBLIC :: alloc_met_type
PUBLIC :: dealloc_met_type
PUBLIC :: assoc_met_type
PUBLIC :: nullify_met_cbl

! Meterological data:
TYPE met_data_type

  INTEGER, ALLOCATABLE :: year      (:)   ! local time year AD
  INTEGER, ALLOCATABLE :: moy       (:)   ! local time month of year 
  REAL, ALLOCATABLE    :: ca        (:)   ! CO2 concentration (mol/mol)                 
  REAL, ALLOCATABLE    :: doy       (:)   ! local time day of year = days since
  REAL, ALLOCATABLE    :: hod       (:)   ! local hour of day                            
  REAL, ALLOCATABLE    :: ofsd      (:)   ! downward SW radiation (W/m2)         
  REAL, ALLOCATABLE    :: fld       (:)   ! downward LW radiation (W/m2)          
  REAL, ALLOCATABLE    :: precip    (:)   ! rainfall (liquid+solid)(mm/dels)             
  REAL, ALLOCATABLE    :: precip_sn (:)   ! solid preipitation only (mm/dels)            
  REAL, ALLOCATABLE    :: tk        (:)   ! surface air temperature (oK)                 
  REAL, ALLOCATABLE    :: tvair     (:)   ! within canopy air temperature (oK)           
  REAL, ALLOCATABLE    :: tvrad     (:)   ! radiative veg. temperature (K)         
  REAL, ALLOCATABLE    :: pmb       (:)   ! surface air pressure (mbar)                  
  REAL, ALLOCATABLE    :: ua        (:)   ! surface wind speed (m/s)                     
  REAL, ALLOCATABLE    :: qv        (:)   ! surface specific humidity (g/g)              
  REAL, ALLOCATABLE    :: coszen    (:)   ! cos(zenith angle of sun)                   
  REAL, ALLOCATABLE    :: Ndep      (:)   ! nitrogen deposition (gN m-2 d-1)              
  REAL, ALLOCATABLE    :: qvair     (:)   ! in canopy specific humidity (g/g)        
  REAL, ALLOCATABLE    :: da        (:)   ! H2O vap pres deficit at ref height (Pa)
  REAL, ALLOCATABLE    :: dva       (:)   ! H2O vap pres deficit in canopy 
  REAL, ALLOCATABLE    :: fsd       (:,:) ! downward SW radiation (W/m2)
                
END TYPE met_data_type

TYPE met_type

  INTEGER, POINTER :: year      (:)   ! local time year AD
  INTEGER, POINTER :: moy       (:)   ! local time month of year 
  REAL, POINTER    :: ca        (:)   ! CO2 concentration (mol/mol)                 
  REAL, POINTER    :: doy       (:)   ! local time day of year = days since
  REAL, POINTER    :: hod       (:)   ! local hour of day                            
  REAL, POINTER    :: ofsd      (:)   ! downward SW radiation (W/m2)         
  REAL, POINTER    :: fld       (:)   ! downward LW radiation (W/m2)          
  REAL, POINTER    :: precip    (:)   ! rainfall (liquid+solid)(mm/dels)             
  REAL, POINTER    :: precip_sn (:)   ! solid preipitation only (mm/dels)            
  REAL, POINTER    :: tk        (:)   ! surface air temperature (oK)                 
  REAL, POINTER    :: tvair     (:)   ! within canopy air temperature (oK)           
  REAL, POINTER    :: tvrad     (:)   ! radiative veg. temperature (K)         
  REAL, POINTER    :: pmb       (:)   ! surface air pressure (mbar)                  
  REAL, POINTER    :: ua        (:)   ! surface wind speed (m/s)                     
  REAL, POINTER    :: qv        (:)   ! surface specific humidity (g/g)              
  REAL, POINTER    :: coszen    (:)   ! cos(zenith angle of sun)                   
  REAL, POINTER    :: Ndep      (:)   ! nitrogen deposition (gN m-2 d-1)              
  REAL, POINTER    :: qvair     (:)   ! in canopy specific humidity (g/g)        
  REAL, POINTER    :: da        (:)   ! H2O vap pres deficit at ref height (Pa)
  REAL, POINTER    :: dva       (:)   ! H2O vap pres deficit in canopy 
  REAL, POINTER    :: fsd       (:,:) ! downward SW radiation (W/m2)
                
END TYPE met_type

CONTAINS

SUBROUTINE alloc_met_type(met, mp)
USE grid_constants_mod_cbl,   ONLY: swb  ! # SW bands 
IMPLICIT NONE

TYPE(met_data_type), INTENT(INOUT) :: met
INTEGER, INTENT(IN) :: mp

ALLOCATE( met% year      (mp) )
ALLOCATE( met% moy       (mp) )
ALLOCATE( met% ca        (mp) )
ALLOCATE( met% doy       (mp) )
ALLOCATE( met% hod       (mp) )
ALLOCATE( met% ofsd      (mp) )
ALLOCATE( met% fld       (mp) )
ALLOCATE( met% precip    (mp) )
ALLOCATE( met% precip_sn (mp) )
ALLOCATE( met% tk        (mp) )
ALLOCATE( met% tvair     (mp) )
ALLOCATE( met% tvrad     (mp) )
ALLOCATE( met% pmb       (mp) )
ALLOCATE( met% ua        (mp) )
ALLOCATE( met% qv        (mp) )
ALLOCATE( met% coszen    (mp) )
ALLOCATE( met% Ndep      (mp) )
ALLOCATE( met% qvair     (mp) )
ALLOCATE( met% da        (mp) )
ALLOCATE( met% dva       (mp) )
ALLOCATE( met% fsd       (mp,swb) )
               
met % year      (:)   = 0.0      
met % moy       (:)   = 0.0      
met % ca        (:)   = 0.0      
met % doy       (:)   = 0.0      
met % hod       (:)   = 0.0      
met % ofsd      (:)   = 0.0      
met % fld       (:)   = 0.0      
met % precip    (:)   = 0.0      
met % precip_sn (:)   = 0.0      
met % tk        (:)   = 0.0      
met % tvair     (:)   = 0.0      
met % tvrad     (:)   = 0.0      
met % pmb       (:)   = 0.0      
met % ua        (:)   = 0.0      
met % qv        (:)   = 0.0      
met % coszen    (:)   = 0.0      
met % Ndep      (:)   = 0.0      
met % qvair     (:)   = 0.0      
met % da        (:)   = 0.0      
met % dva       (:)   = 0.0      
met % fsd       (:,:) = 0.0      

RETURN
END SUBROUTINE alloc_met_type

SUBROUTINE dealloc_met_type(met)
IMPLICIT NONE

TYPE(met_type), INTENT(inout) :: met

DEALLOCATE ( met % year      )
DEALLOCATE ( met % moy       )
DEALLOCATE ( met % ca        )
DEALLOCATE ( met % doy       )
DEALLOCATE ( met % hod       )
DEALLOCATE ( met % ofsd      )
DEALLOCATE ( met % fld       )
DEALLOCATE ( met % precip    )
DEALLOCATE ( met % precip_sn )
DEALLOCATE ( met % tk        )
DEALLOCATE ( met % tvair     )
DEALLOCATE ( met % tvrad     )
DEALLOCATE ( met % pmb       )
DEALLOCATE ( met % ua        )
DEALLOCATE ( met % qv        )
DEALLOCATE ( met % coszen    )
DEALLOCATE ( met % Ndep      )
DEALLOCATE ( met % qvair     )
DEALLOCATE ( met % da        )
DEALLOCATE ( met % dva       )
DEALLOCATE ( met % fsd       )

RETURN
END SUBROUTINE dealloc_met_type

SUBROUTINE assoc_met_type(met, met_data )
! Description:
!   Associate the CABLE work pointers in the derived type structure
IMPLICIT NONE

!Arguments
TYPE(met_type),      INTENT(IN OUT)         :: met
TYPE(met_data_type), INTENT(IN OUT), TARGET :: met_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_met_cbl(met)

met% year      => met_data% year                  
met% moy       => met_data% moy                   
met% ca        => met_data% ca                    
met% doy       => met_data% doy                   
met% hod       => met_data% hod                   
met% ofsd      => met_data% ofsd                  
met% fld       => met_data% fld                   
met% precip    => met_data% precip                
met% precip_sn => met_data% precip_sn             
met% tk        => met_data% tk                    
met% tvair     => met_data% tvair                 
met% tvrad     => met_data% tvrad                 
met% pmb       => met_data% pmb                   
met% ua        => met_data% ua                    
met% qv        => met_data% qv                    
met% coszen    => met_data% coszen                
met% Ndep      => met_data% Ndep                  
met% qvair     => met_data% qvair                 
met% da        => met_data% da                    
met% dva       => met_data% dva                   
met% fsd       => met_data% fsd                   

RETURN
END SUBROUTINE assoc_met_type

SUBROUTINE nullify_met_cbl( met )
! Description:
!   Nullify the CABLE work pointers in the derived type structure
IMPLICIT NONE

!Arguments
TYPE(met_type), INTENT(IN OUT) :: met 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( met % year      )
NULLIFY( met % moy       )
NULLIFY( met % ca        )
NULLIFY( met % doy       )
NULLIFY( met % hod       )
NULLIFY( met % ofsd      )
NULLIFY( met % fld       )
NULLIFY( met % precip    )
NULLIFY( met % precip_sn )
NULLIFY( met % tk        )
NULLIFY( met % tvair     )
NULLIFY( met % tvrad     )
NULLIFY( met % pmb       )
NULLIFY( met % ua        )
NULLIFY( met % qv        )
NULLIFY( met % coszen    )
NULLIFY( met % Ndep      )
NULLIFY( met % qvair     )
NULLIFY( met % da        )
NULLIFY( met % dva       )
NULLIFY( met % fsd       )

RETURN

END SUBROUTINE nullify_met_cbl

END MODULE cable_met_type_mod
