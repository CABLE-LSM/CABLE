MODULE cable_radiation_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: radiation_type
PUBLIC :: radiation_data_type
PUBLIC :: alloc_radiation_type
PUBLIC :: dealloc_radiation_type
PUBLIC :: assoc_radiation_type
PUBLIC :: nullify_radiation_cbl

! Radiation variables:
TYPE radiation_data_type

  REAL, ALLOCATABLE :: transb    (:) ! fraction SW beam tranmitted through canopy          
  REAL, ALLOCATABLE :: albedo_T  (:) ! canopy+soil albedo for VIS+NIR                      
  REAL, ALLOCATABLE :: longitude (:) ! longitude                                           
  REAL, ALLOCATABLE :: workp1    (:) ! absorbed short-wave radiation for soil              
  REAL, ALLOCATABLE :: workp2    (:) ! absorbed short-wave radiation for soil              
  REAL, ALLOCATABLE :: workp3    (:) ! absorbed short-wave radiation for soil              
  REAL, ALLOCATABLE :: extkb     (:) ! beam radiation extinction coeff                     
  REAL, ALLOCATABLE :: extkd2    (:) ! diffuse 2D radiation extinction coeff               
  REAL, ALLOCATABLE :: extkd     (:) ! diffuse radiation extinction coeff (-)              
  REAL, ALLOCATABLE :: flws      (:) ! soil long-wave radiation                            
  REAL, ALLOCATABLE :: latitude  (:) ! latitude                                            
  REAL, ALLOCATABLE :: lwabv     (:) ! long wave absorbed by vegetation                    
  REAL, ALLOCATABLE :: qssabs    (:) ! absorbed short-wave radiation for soil              
  REAL, ALLOCATABLE :: transd    (:) ! frac SW diffuse transmitted through canopy          
  REAL, ALLOCATABLE :: trad      (:) !  radiative temperature (soil and veg)               
  REAL, ALLOCATABLE :: otrad     (:) ! radiative temperature on previous timestep (ACCESS) 

  REAL, ALLOCATABLE :: fvlai     (:,:) ! leaf area index of big leaf                    
  REAL, ALLOCATABLE :: rhocdf    (:,:) ! canopy diffuse reflectance (-)                 
  REAL, ALLOCATABLE :: rniso     (:,:) ! sum(rad%qcan, 3) total abs by canopy (W/m2)    
  REAL, ALLOCATABLE :: scalex    (:,:) ! scaling PARAMETER for big leaf                 
  REAL, ALLOCATABLE :: albedo    (:,:) ! canopy+soil albedo                             
  REAL, ALLOCATABLE :: reffdf    (:,:) ! effective conopy diffuse reflectance           
  REAL, ALLOCATABLE :: reffbm    (:,:) ! effective conopy beam reflectance              
  REAL, ALLOCATABLE :: extkbm    (:,:) ! modified k beam(6.20)(for leaf scattering)     
  REAL, ALLOCATABLE :: extkdm    (:,:) ! modified k diffuse(6.20)(for leaf scattering)  
  REAL, ALLOCATABLE :: fbeam     (:,:) ! beam fraction                                  
  REAL, ALLOCATABLE :: cexpkbm   (:,:) ! canopy beam transmittance                      
  REAL, ALLOCATABLE :: cexpkdm   (:,:) ! canopy diffuse transmittance                   
  REAL, ALLOCATABLE :: rhocbm    (:,:) ! modified canopy beam reflectance(6.21)         
  REAL, ALLOCATABLE :: gradis    (:,:) ! radiative conductance                          
  
  REAL, ALLOCATABLE :: qcan      (:,:,:) ! absorbed radiation for canopy (W/m^2) 

END TYPE radiation_data_type

TYPE radiation_type

  REAL, POINTER :: transb    (:) ! fraction SW beam tranmitted through canopy          
  REAL, POINTER :: albedo_T  (:) ! canopy+soil albedo for VIS+NIR                      
  REAL, POINTER :: longitude (:) ! longitude                                           
  REAL, POINTER :: workp1    (:) ! absorbed short-wave radiation for soil              
  REAL, POINTER :: workp2    (:) ! absorbed short-wave radiation for soil              
  REAL, POINTER :: workp3    (:) ! absorbed short-wave radiation for soil              
  REAL, POINTER :: extkb     (:) ! beam radiation extinction coeff                     
  REAL, POINTER :: extkd2    (:) ! diffuse 2D radiation extinction coeff               
  REAL, POINTER :: extkd     (:) ! diffuse radiation extinction coeff (-)              
  REAL, POINTER :: flws      (:) ! soil long-wave radiation                            
  REAL, POINTER :: latitude  (:) ! latitude                                            
  REAL, POINTER :: lwabv     (:) ! long wave absorbed by vegetation                    
  REAL, POINTER :: qssabs    (:) ! absorbed short-wave radiation for soil              
  REAL, POINTER :: transd    (:) ! frac SW diffuse transmitted through canopy          
  REAL, POINTER :: trad      (:) !  radiative temperature (soil and veg)               
  REAL, POINTER :: otrad     (:) ! radiative temperature on previous timestep (ACCESS) 

  REAL, POINTER :: fvlai     (:,:) ! leaf area index of big leaf                    
  REAL, POINTER :: rhocdf    (:,:) ! canopy diffuse reflectance (-)                 
  REAL, POINTER :: rniso     (:,:) ! sum(rad%qcan, 3) total abs by canopy (W/m2)    
  REAL, POINTER :: scalex    (:,:) ! scaling PARAMETER for big leaf                 
  REAL, POINTER :: albedo    (:,:) ! canopy+soil albedo                             
  REAL, POINTER :: reffdf    (:,:) ! effective conopy diffuse reflectance           
  REAL, POINTER :: reffbm    (:,:) ! effective conopy beam reflectance              
  REAL, POINTER :: extkbm    (:,:) ! modified k beam(6.20)(for leaf scattering)     
  REAL, POINTER :: extkdm    (:,:) ! modified k diffuse(6.20)(for leaf scattering)  
  REAL, POINTER :: fbeam     (:,:) ! beam fraction                                  
  REAL, POINTER :: cexpkbm   (:,:) ! canopy beam transmittance                      
  REAL, POINTER :: cexpkdm   (:,:) ! canopy diffuse transmittance                   
  REAL, POINTER :: rhocbm    (:,:) ! modified canopy beam reflectance(6.21)         
  REAL, POINTER :: gradis    (:,:) ! radiative conductance                          
  
  REAL, POINTER :: qcan      (:,:,:) ! absorbed radiation for canopy (W/m^2) 

END TYPE radiation_type

CONTAINS

SUBROUTINE alloc_radiation_type(radiation, mp)

USE grid_constants_mod_cbl,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: nrb, swb         ! # Radiation/SW bands 

IMPLICIT NONE

TYPE(radiation_data_type), INTENT(INOUT) :: radiation
INTEGER, INTENT(IN) :: mp

ALLOCATE( radiation% transb    (mp) )
ALLOCATE( radiation% albedo_T  (mp) )
ALLOCATE( radiation% longitude (mp) )
ALLOCATE( radiation% workp1    (mp) )
ALLOCATE( radiation% workp2    (mp) )
ALLOCATE( radiation% workp3    (mp) )
ALLOCATE( radiation% extkb     (mp) )
ALLOCATE( radiation% extkd2    (mp) )
ALLOCATE( radiation% extkd     (mp) )
ALLOCATE( radiation% flws      (mp) )
ALLOCATE( radiation% latitude  (mp) )
ALLOCATE( radiation% lwabv     (mp) )
ALLOCATE( radiation% qssabs    (mp) )
ALLOCATE( radiation% transd    (mp) )
ALLOCATE( radiation% trad      (mp) )
ALLOCATE( radiation% otrad     (mp) )

ALLOCATE( radiation% fvlai     (mp,mf ) )
ALLOCATE( radiation% rhocdf    (mp,nrb) )
ALLOCATE( radiation% rniso     (mp,mf ) )
ALLOCATE( radiation% scalex    (mp,mf ) )
ALLOCATE( radiation% albedo    (mp,nrb) )
ALLOCATE( radiation% reffdf    (mp,nrb) )
ALLOCATE( radiation% reffbm    (mp,nrb) )
ALLOCATE( radiation% extkbm    (mp,nrb) )
ALLOCATE( radiation% extkdm    (mp,nrb) )
ALLOCATE( radiation% fbeam     (mp,nrb) )
ALLOCATE( radiation% cexpkbm   (mp,swb) )
ALLOCATE( radiation% cexpkdm   (mp,swb) )
ALLOCATE( radiation% rhocbm    (mp,nrb) )
ALLOCATE( radiation% gradis    (mp,mf ) )
                     
ALLOCATE( radiation% qcan      (mp,mf,nrb) )

radiation % transb    (:)     = 0.0      
radiation % albedo_T  (:)     = 0.0      
radiation % longitude (:)     = 0.0      
radiation % workp1    (:)     = 0.0      
radiation % workp2    (:)     = 0.0      
radiation % workp3    (:)     = 0.0      
radiation % extkb     (:)     = 0.0      
radiation % extkd2    (:)     = 0.0      
radiation % extkd     (:)     = 0.0      
radiation % flws      (:)     = 0.0      
radiation % latitude  (:)     = 0.0      
radiation % lwabv     (:)     = 0.0      
radiation % qssabs    (:)     = 0.0      
radiation % transd    (:)     = 0.0      
radiation % trad      (:)     = 0.0      
radiation % otrad     (:)     = 0.0      
radiation % fvlai     (:,:)   = 0.0      
radiation % rhocdf    (:,:)   = 0.0      
radiation % rniso     (:,:)   = 0.0      
radiation % scalex    (:,:)   = 0.0      
radiation % albedo    (:,:)   = 0.0      
radiation % reffdf    (:,:)   = 0.0      
radiation % reffbm    (:,:)   = 0.0      
radiation % extkbm    (:,:)   = 0.0      
radiation % extkdm    (:,:)   = 0.0      
radiation % fbeam     (:,:)   = 0.0      
radiation % cexpkbm   (:,:)   = 0.0      
radiation % cexpkdm   (:,:)   = 0.0      
radiation % rhocbm    (:,:)   = 0.0      
radiation % gradis    (:,:)   = 0.0      
radiation % qcan      (:,:,:) = 0.0      

RETURN
END SUBROUTINE alloc_radiation_type

SUBROUTINE dealloc_radiation_type(radiation)

TYPE(radiation_type), INTENT(inout) :: radiation

DEALLOCATE ( radiation % transb    )
DEALLOCATE ( radiation % albedo_T  )
DEALLOCATE ( radiation % longitude )
DEALLOCATE ( radiation % workp1    )
DEALLOCATE ( radiation % workp2    )
DEALLOCATE ( radiation % workp3    )
DEALLOCATE ( radiation % extkb     )
DEALLOCATE ( radiation % extkd2    )
DEALLOCATE ( radiation % extkd     )
DEALLOCATE ( radiation % flws      )
DEALLOCATE ( radiation % latitude  )
DEALLOCATE ( radiation % lwabv     )
DEALLOCATE ( radiation % qssabs    )
DEALLOCATE ( radiation % transd    )
DEALLOCATE ( radiation % trad      )
DEALLOCATE ( radiation % otrad     )
DEALLOCATE ( radiation % fvlai     )
DEALLOCATE ( radiation % rhocdf    )
DEALLOCATE ( radiation % rniso     )
DEALLOCATE ( radiation % scalex    )
DEALLOCATE ( radiation % albedo    )
DEALLOCATE ( radiation % reffdf    )
DEALLOCATE ( radiation % reffbm    )
DEALLOCATE ( radiation % extkbm    )
DEALLOCATE ( radiation % extkdm    )
DEALLOCATE ( radiation % fbeam     )
DEALLOCATE ( radiation % cexpkbm   )
DEALLOCATE ( radiation % cexpkdm   )
DEALLOCATE ( radiation % rhocbm    )
DEALLOCATE ( radiation % gradis    )
DEALLOCATE ( radiation % qcan      )

RETURN
END SUBROUTINE dealloc_radiation_type

SUBROUTINE assoc_radiation_type(radiation, radiation_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(radiation_type),      INTENT(IN OUT)         :: radiation
TYPE(radiation_data_type), INTENT(IN OUT), TARGET :: radiation_data
CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_radiation_cbl(radiation)

radiation% transb    => radiation_data% transb   
radiation% albedo_T  => radiation_data% albedo_T 
radiation% longitude => radiation_data% longitude
radiation% workp1    => radiation_data% workp1   
radiation% workp2    => radiation_data% workp2   
radiation% workp3    => radiation_data% workp3   
radiation% extkb     => radiation_data% extkb    
radiation% extkd2    => radiation_data% extkd2   
radiation% extkd     => radiation_data% extkd    
radiation% flws      => radiation_data% flws     
radiation% latitude  => radiation_data% latitude 
radiation% lwabv     => radiation_data% lwabv    
radiation% qssabs    => radiation_data% qssabs   
radiation% transd    => radiation_data% transd   
radiation% trad      => radiation_data% trad     
radiation% otrad     => radiation_data% otrad    
radiation% fvlai     => radiation_data% fvlai    
radiation% rhocdf    => radiation_data% rhocdf   
radiation% rniso     => radiation_data% rniso    
radiation% scalex    => radiation_data% scalex   
radiation% albedo    => radiation_data% albedo   
radiation% reffdf    => radiation_data% reffdf   
radiation% reffbm    => radiation_data% reffbm   
radiation% extkbm    => radiation_data% extkbm   
radiation% extkdm    => radiation_data% extkdm   
radiation% fbeam     => radiation_data% fbeam    
radiation% cexpkbm   => radiation_data% cexpkbm  
radiation% cexpkdm   => radiation_data% cexpkdm  
radiation% rhocbm    => radiation_data% rhocbm   
radiation% gradis    => radiation_data% gradis   
radiation% qcan      => radiation_data% qcan     

RETURN
END SUBROUTINE assoc_radiation_type

SUBROUTINE nullify_radiation_cbl( radiation )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(radiation_type), INTENT(IN OUT) :: radiation 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( radiation % transb    )
NULLIFY( radiation % albedo_T  )
NULLIFY( radiation % longitude )
NULLIFY( radiation % workp1    )
NULLIFY( radiation % workp2    )
NULLIFY( radiation % workp3    )
NULLIFY( radiation % extkb     )
NULLIFY( radiation % extkd2    )
NULLIFY( radiation % extkd     )
NULLIFY( radiation % flws      )
NULLIFY( radiation % latitude  )
NULLIFY( radiation % lwabv     )
NULLIFY( radiation % qssabs    )
NULLIFY( radiation % transd    )
NULLIFY( radiation % trad      )
NULLIFY( radiation % otrad     )
NULLIFY( radiation % fvlai     )
NULLIFY( radiation % rhocdf    )
NULLIFY( radiation % rniso     )
NULLIFY( radiation % scalex    )
NULLIFY( radiation % albedo    )
NULLIFY( radiation % reffdf    )
NULLIFY( radiation % reffbm    )
NULLIFY( radiation % extkbm    )
NULLIFY( radiation % extkdm    )
NULLIFY( radiation % fbeam     )
NULLIFY( radiation % cexpkbm   )
NULLIFY( radiation % cexpkdm   )
NULLIFY( radiation % rhocbm    )
NULLIFY( radiation % gradis    )
NULLIFY( radiation % qcan      )

RETURN

END SUBROUTINE nullify_radiation_cbl

END MODULE cable_radiation_type_mod
