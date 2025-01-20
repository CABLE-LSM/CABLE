MODULE cable_climate_type_mod
!* Description:
! This module defines the climate and climate_data types and allocates 
! arrays and pointers for these types.

IMPLICIT NONE

PUBLIC :: climate_type
PUBLIC :: climate_data_type
PUBLIC :: alloc_climate_type
PUBLIC :: dealloc_climate_type
PUBLIC :: assoc_climate_type
PUBLIC :: nullify_climate_cbl

! Climate data:
TYPE climate_data_type

  INTEGER :: nyear_average = 20
  INTEGER :: nday_average  = 31
  INTEGER :: nyears    
    !! number of years in climate record
  INTEGER :: doy 
    !! day of year

  INTEGER, ALLOCATABLE :: chilldays (:)          
    !! length of chilling period (period with T<5deg) 
  INTEGER, ALLOCATABLE :: iveg      (:)          
    !! potential vegetation type based on climatic constraints 
  INTEGER, ALLOCATABLE :: biome     (:)  

  REAL, ALLOCATABLE :: dtemp               (:)   
    !! daily temperature                                            
  REAL, ALLOCATABLE :: dmoist              (:)   
    !! daily moisture availability                                  
  REAL, ALLOCATABLE :: mtemp               (:)   
    !! mean temperature over the last 31 days                       
  REAL, ALLOCATABLE :: qtemp               (:)   
    !! mean temperature over the last 91 days                       
  REAL, ALLOCATABLE :: mmoist              (:)   
    !! monthly moisture availability                                
  REAL, ALLOCATABLE :: mtemp_min           (:)   
    !! minimum monthly temperature                                  
  REAL, ALLOCATABLE :: mtemp_max           (:)   
    !! maximum monthly temperature                                  
  REAL, ALLOCATABLE :: qtemp_max           (:)   
    !! mean temperature of the warmest quarter (so far this year)   
  REAL, ALLOCATABLE :: mtemp_min20         (:)   
    !! minimum monthly temperature, averaged over 20 y              
  REAL, ALLOCATABLE :: mtemp_max20         (:)   
    !! maximum monthly temperature, averaged over 20 y              
  REAL, ALLOCATABLE :: atemp_mean          (:)   
    !! annual average temperature                                   
  REAL, ALLOCATABLE :: AGDD5               (:)                                                                  
  REAL, ALLOCATABLE :: GDD5                (:)   
    !! growing degree day sum relative to 5deg base temperature     
  REAL, ALLOCATABLE :: AGDD0               (:)                                                                  
  REAL, ALLOCATABLE :: GDD0                (:)   
    !! growing degree day sum relative to 0deg base temperature     
  REAL, ALLOCATABLE :: alpha_PT            (:)   
    !! ratio of annual evap to annual PT evap                       
  REAL, ALLOCATABLE :: evap_PT             (:)   
    !! annual PT evap [mm]                                          
  REAL, ALLOCATABLE :: aevap               (:)   
    !! annual evap [mm]                                             
  REAL, ALLOCATABLE :: alpha_PT20          (:)                                                                  
  REAL, ALLOCATABLE :: qtemp_max_last_year (:)   
    !! mean temperature of the warmest quarter (last calendar year) 

  REAL, ALLOCATABLE :: mtemp_min_20        (:,:) 
    !! minimum monthly temperatures for the last 20 y   
  REAL, ALLOCATABLE :: mtemp_max_20        (:,:) 
    !! maximum monthly temperatures for the last 20 y   
  REAL, ALLOCATABLE :: dtemp_31            (:,:) 
    !! daily temperature for the last 31 days           
  REAL, ALLOCATABLE :: dmoist_31           (:,:) 
    !! daily moisture availability for the last 31 days 
  REAL, ALLOCATABLE :: alpha_PT_20         (:,:) 
    !! priestley Taylor Coefft for last 20 y            
  REAL, ALLOCATABLE :: dtemp_91            (:,:) 
    !! daily temperature for the last 91 days           

END TYPE climate_data_type

TYPE climate_type

  INTEGER, POINTER :: nyear_average
  INTEGER, POINTER :: nday_average
  INTEGER, POINTER :: nyears                 
    !! number of years in climate record
  INTEGER, POINTER :: doy                    
    !! day of year

  INTEGER, POINTER :: chilldays (:)          
    !! length of chilling period (period with T<5deg) 
  INTEGER, POINTER :: iveg      (:)          
    !! potential vegetation type based on climatic constraints 
  INTEGER, POINTER :: biome     (:)  

  REAL, POINTER :: dtemp               (:)   
    !! daily temperature                                            
  REAL, POINTER :: dmoist              (:)   
    !! daily moisture availability                                  
  REAL, POINTER :: mtemp               (:)   
    !! mean temperature over the last 31 days                       
  REAL, POINTER :: qtemp               (:)   
    !! mean temperature over the last 91 days                       
  REAL, POINTER :: mmoist              (:)   
    !! monthly moisture availability                                
  REAL, POINTER :: mtemp_min           (:)   
    !! minimum monthly temperature                                  
  REAL, POINTER :: mtemp_max           (:)   
    !! maximum monthly temperature                                  
  REAL, POINTER :: qtemp_max           (:)   
    !! mean temperature of the warmest quarter (so far this year)   
  REAL, POINTER :: mtemp_min20         (:)   
    !! minimum monthly temperature, averaged over 20 y              
  REAL, POINTER :: mtemp_max20         (:)   
    !! maximum monthly temperature, averaged over 20 y              
  REAL, POINTER :: atemp_mean          (:)   
    !! annual average temperature                                   
  REAL, POINTER :: AGDD5               (:)                                                                  
  REAL, POINTER :: GDD5                (:)   
    !! growing degree day sum relative to 5deg base temperature     
  REAL, POINTER :: AGDD0               (:)                                                                  
  REAL, POINTER :: GDD0                (:)   
    !! growing degree day sum relative to 0deg base temperature     
  REAL, POINTER :: alpha_PT            (:)   
    !! ratio of annual evap to annual PT evap                       
  REAL, POINTER :: evap_PT             (:)   
    !! annual PT evap [mm]                                          
  REAL, POINTER :: aevap               (:)   
    !! annual evap [mm]                                             
  REAL, POINTER :: alpha_PT20          (:)                                                                  
  REAL, POINTER :: qtemp_max_last_year (:)   
    !! mean temperature of the warmest quarter (last calendar year) 

  REAL, POINTER :: mtemp_min_20        (:,:) 
    !! minimum monthly temperatures for the last 20 y   
  REAL, POINTER :: mtemp_max_20        (:,:) 
    !! maximum monthly temperatures for the last 20 y   
  REAL, POINTER :: dtemp_31            (:,:) 
    !! daily temperature for the last 31 days           
  REAL, POINTER :: dmoist_31           (:,:) 
    !! daily moisture availability for the last 31 days 
  REAL, POINTER :: alpha_PT_20         (:,:) 
    !! priestley Taylor Coefft for last 20 y            
  REAL, POINTER :: dtemp_91            (:,:) 
    !! daily temperature for the last 91 days           

END TYPE climate_type

CONTAINS

SUBROUTINE alloc_climate_type(climate, mp)
!* Description:
! Allocate and initialise arrays in the climate_data type

USE grid_constants_mod_cbl,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: niter            ! number of iterations for za/L
USE cable_common_module,      ONLY: cable_runtime

IMPLICIT NONE

TYPE(climate_data_type), INTENT(INOUT) :: climate
INTEGER, INTENT(IN) :: mp
!! Number of tiles/patches in the CABLE simulation

INTEGER :: ny
INTEGER :: nd
INTEGER, PARAMETER :: ns = 91 

ny = climate%nyear_average
nd = climate%nday_average

IF ( cable_runtime%um ) THEN
  ALLOCATE( climate% chilldays   (1) )
  ALLOCATE( climate% iveg        (1) )
  ALLOCATE( climate% biome       (1) )
  ALLOCATE( climate% dtemp       (1) )
  ALLOCATE( climate% dmoist      (1) )
  ALLOCATE( climate% mtemp       (1) )
  ALLOCATE( climate% qtemp       (1) )
  ALLOCATE( climate% mmoist      (1) )
  ALLOCATE( climate% mtemp_min   (1) )
  ALLOCATE( climate% mtemp_max   (1) )
  ALLOCATE( climate% qtemp_max   (1) )
  ALLOCATE( climate% mtemp_min20 (1) )
  ALLOCATE( climate% mtemp_max20 (1) )
  ALLOCATE( climate% atemp_mean  (1) )
  ALLOCATE( climate% AGDD5       (1) )
  ALLOCATE( climate% GDD5        (1) )
  ALLOCATE( climate% AGDD0       (1) )
  ALLOCATE( climate% GDD0        (1) )
  ALLOCATE( climate% alpha_PT    (1) )
  ALLOCATE( climate% evap_PT     (1) )
  ALLOCATE( climate% aevap       (1) )
  ALLOCATE( climate% alpha_PT20  (1) )
  ALLOCATE( climate% qtemp_max_last_year (1) )
  ALLOCATE( climate% mtemp_min_20 (1,1) )
  ALLOCATE( climate% mtemp_max_20 (1,1) )
  ALLOCATE( climate% dtemp_31     (1,1) )
  ALLOCATE( climate% dmoist_31    (1,1) )
  ALLOCATE( climate% alpha_PT_20  (1,1) )
  ALLOCATE( climate% dtemp_91     (1,1) )

ELSE

  ALLOCATE( climate% chilldays   (mp) )
  ALLOCATE( climate% iveg        (mp) )
  ALLOCATE( climate% biome       (mp) )
  ALLOCATE( climate% dtemp       (mp) )
  ALLOCATE( climate% dmoist      (mp) )
  ALLOCATE( climate% mtemp       (mp) )
  ALLOCATE( climate% qtemp       (mp) )
  ALLOCATE( climate% mmoist      (mp) )
  ALLOCATE( climate% mtemp_min   (mp) )
  ALLOCATE( climate% mtemp_max   (mp) )
  ALLOCATE( climate% qtemp_max   (mp) )
  ALLOCATE( climate% mtemp_min20 (mp) )
  ALLOCATE( climate% mtemp_max20 (mp) )
  ALLOCATE( climate% atemp_mean  (mp) )
  ALLOCATE( climate% AGDD5       (mp) )
  ALLOCATE( climate% GDD5        (mp) )
  ALLOCATE( climate% AGDD0       (mp) )
  ALLOCATE( climate% GDD0        (mp) )
  ALLOCATE( climate% alpha_PT    (mp) )
  ALLOCATE( climate% evap_PT     (mp) )
  ALLOCATE( climate% aevap       (mp) )
  ALLOCATE( climate% alpha_PT20  (mp) )
  ALLOCATE( climate% qtemp_max_last_year (mp) )
  ALLOCATE( climate% mtemp_min_20 (mp,ny) )
  ALLOCATE( climate% mtemp_max_20 (mp,ny) )
  ALLOCATE( climate% dtemp_31     (mp,nd) )
  ALLOCATE( climate% dmoist_31    (mp,nd) )
  ALLOCATE( climate% alpha_PT_20  (mp,ny) )
  ALLOCATE( climate% dtemp_91     (mp,ns) )

END  IF

climate% chilldays    (:)   = 0.0      
climate% iveg         (:)   = 0.0      
climate% biome        (:)   = 0.0      
climate% dtemp        (:)   = 0.0      
climate% dmoist       (:)   = 0.0      
climate% mtemp        (:)   = 0.0      
climate% qtemp        (:)   = 0.0      
climate% mmoist       (:)   = 0.0      
climate% mtemp_min    (:)   = 0.0      
climate% mtemp_max    (:)   = 0.0      
climate% qtemp_max    (:)   = 0.0      
climate% mtemp_min20  (:)   = 0.0      
climate% mtemp_max20  (:)   = 0.0      
climate% atemp_mean   (:)   = 0.0      
climate% AGDD5        (:)   = 0.0      
climate% GDD5         (:)   = 0.0      
climate% AGDD0        (:)   = 0.0      
climate% GDD0         (:)   = 0.0      
climate% alpha_PT     (:)   = 0.0      
climate% evap_PT      (:)   = 0.0      
climate% aevap        (:)   = 0.0      
climate% alpha_PT20   (:)   = 0.0      
climate% qtemp_max_last_year   (:)   = 0.0      
climate% mtemp_min_20 (:,:) = 0.0      
climate% mtemp_max_20 (:,:) = 0.0      
climate% dtemp_31     (:,:) = 0.0      
climate% dmoist_31    (:,:) = 0.0      
climate% alpha_PT_20  (:,:) = 0.0      
climate% dtemp_91     (:,:) = 0.0      

RETURN
END SUBROUTINE alloc_climate_type

SUBROUTINE dealloc_climate_type(climate)
!* Description:
! Deallocate arrays in the climate_data type

TYPE(climate_type), INTENT(inout) :: climate

DEALLOCATE ( climate% chilldays      )
DEALLOCATE ( climate% iveg           )
DEALLOCATE ( climate% biome          )
DEALLOCATE ( climate% dtemp          )
DEALLOCATE ( climate% dmoist         )
DEALLOCATE ( climate% mtemp          )
DEALLOCATE ( climate% qtemp          )
DEALLOCATE ( climate% mmoist         )
DEALLOCATE ( climate% mtemp_min      )
DEALLOCATE ( climate% mtemp_max      )
DEALLOCATE ( climate% qtemp_max      )
DEALLOCATE ( climate% mtemp_min20    )
DEALLOCATE ( climate% mtemp_max20    )
DEALLOCATE ( climate% atemp_mean     )
DEALLOCATE ( climate% AGDD5          )
DEALLOCATE ( climate% GDD5           )
DEALLOCATE ( climate% AGDD0          )
DEALLOCATE ( climate% GDD0           )
DEALLOCATE ( climate% alpha_PT       )
DEALLOCATE ( climate% evap_PT        )
DEALLOCATE ( climate% aevap          )
DEALLOCATE ( climate% alpha_PT20     )
DEALLOCATE ( climate % qtemp_max_last_year )
DEALLOCATE ( climate % mtemp_min_20  )
DEALLOCATE ( climate % mtemp_max_20  )
DEALLOCATE ( climate % dtemp_31      )
DEALLOCATE ( climate % dmoist_31     )
DEALLOCATE ( climate % alpha_PT_20   )
DEALLOCATE ( climate % dtemp_91      )

RETURN
END SUBROUTINE dealloc_climate_type

SUBROUTINE assoc_climate_type(climate, climate_data )
!* Description:
!   Associate the climate type pointers to the climate_data arrays.

IMPLICIT NONE

!Arguments
TYPE(climate_type),      INTENT(IN OUT)         :: climate
TYPE(climate_data_type), INTENT(IN OUT), TARGET :: climate_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_climate_cbl(climate)

climate% chilldays                 => climate_data% chilldays  
climate% iveg                      => climate_data% iveg       
climate% biome                     => climate_data% biome      
climate% dtemp                     => climate_data% dtemp      
climate% dmoist                    => climate_data% dmoist     
climate% mtemp                     => climate_data% mtemp      
climate% qtemp                     => climate_data% qtemp      
climate% mmoist                    => climate_data% mmoist     
climate% mtemp_min                 => climate_data% mtemp_min  
climate% mtemp_max                 => climate_data% mtemp_max  
climate% qtemp_max                 => climate_data% qtemp_max  
climate% mtemp_min20               => climate_data% mtemp_min20
climate% mtemp_max20               => climate_data% mtemp_max20
climate% atemp_mean                => climate_data% atemp_mean 
climate% AGDD5                     => climate_data% AGDD5      
climate% GDD5                      => climate_data% GDD5       
climate% AGDD0                     => climate_data% AGDD0      
climate% GDD0                      => climate_data% GDD0       
climate% alpha_PT                  => climate_data% alpha_PT   
climate% evap_PT                   => climate_data% evap_PT    
climate% aevap                     => climate_data% aevap      
climate% alpha_PT20                => climate_data% alpha_PT20 
climate% qtemp_max_last_year       => climate_data% qtemp_max_last_year 
climate% mtemp_min_20              => climate_data% mtemp_min_20        
climate% mtemp_max_20              => climate_data% mtemp_max_20        
climate% dtemp_31                  => climate_data% dtemp_31            
climate% dmoist_31                 => climate_data% dmoist_31           
climate% alpha_PT_20               => climate_data% alpha_PT_20         
climate% dtemp_91                  => climate_data% dtemp_91            
           
RETURN
END SUBROUTINE assoc_climate_type

SUBROUTINE nullify_climate_cbl( climate )
!* Description:
!   Nullify the climate type pointers.

IMPLICIT NONE

!Arguments
TYPE(climate_type), INTENT(IN OUT) :: climate 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( climate% chilldays     )
NULLIFY( climate% iveg          )
NULLIFY( climate% biome         )
NULLIFY( climate% dtemp         )
NULLIFY( climate% dmoist        )
NULLIFY( climate% mtemp         )
NULLIFY( climate% qtemp         )
NULLIFY( climate% mmoist        )
NULLIFY( climate% mtemp_min     )
NULLIFY( climate% mtemp_max     )
NULLIFY( climate% qtemp_max     )
NULLIFY( climate% mtemp_min20   )
NULLIFY( climate% mtemp_max20   )
NULLIFY( climate% atemp_mean    )
NULLIFY( climate% AGDD5         )
NULLIFY( climate% GDD5          )
NULLIFY( climate% AGDD0         )
NULLIFY( climate% GDD0          )
NULLIFY( climate% alpha_PT      )
NULLIFY( climate% evap_PT       )
NULLIFY( climate% aevap         )
NULLIFY( climate% alpha_PT20    )
NULLIFY( climate % qtemp_max_last_year )
NULLIFY( climate % mtemp_min_20 )
NULLIFY( climate % mtemp_max_20 )
NULLIFY( climate % dtemp_31     )
NULLIFY( climate % dmoist_31    )
NULLIFY( climate % alpha_PT_20  )
NULLIFY( climate % dtemp_91     )

RETURN

END SUBROUTINE nullify_climate_cbl

END MODULE cable_climate_type_mod
