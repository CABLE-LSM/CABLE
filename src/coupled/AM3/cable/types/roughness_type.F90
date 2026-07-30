MODULE cable_roughness_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: roughness_type
PUBLIC :: roughness_data_type
PUBLIC :: alloc_roughness_type
PUBLIC :: dealloc_roughness_type
PUBLIC :: assoc_roughness_type
PUBLIC :: nullify_roughness_cbl

! Roughness variables:
TYPE roughness_data_type

  REAL, ALLOCATABLE :: disp       (:) ! zero-plane displacement                         
  REAL, ALLOCATABLE :: hruff      (:) ! canopy height above snow level                  
  REAL, ALLOCATABLE :: hruff_grmx (:) ! max ht of canopy from tiles on same grid        
  REAL, ALLOCATABLE :: rt0us      (:) ! eq. 3.54, SCAM manual (CSIRO tech report 132)   
  REAL, ALLOCATABLE :: rt1usa     (:) ! resistance from disp to hruf                    
  REAL, ALLOCATABLE :: rt1usb     (:) ! resist fr hruf to zruffs (zref if zref<zruffs)  
  REAL, ALLOCATABLE :: rt1        (:) ! 1/aerodynamic conductance                       
  REAL, ALLOCATABLE :: za_uv      (:) ! level of lowest atmospheric model layer         
  REAL, ALLOCATABLE :: za_tq      (:) ! level of lowest atmospheric model layer         
  REAL, ALLOCATABLE :: z0m        (:) ! roughness length                                
  REAL, ALLOCATABLE :: zref_uv    (:) ! Reference height for met forcing                
  REAL, ALLOCATABLE :: zref_tq    (:) ! Reference height for met forcing                
  REAL, ALLOCATABLE :: zruffs     (:) ! SCALAR Roughness sublayer depth (ground=origin) 
  REAL, ALLOCATABLE :: z0soilsn   (:) ! roughness length of bare soil surface           
  REAL, ALLOCATABLE :: z0soil     (:) ! roughness length of bare soil surface           
  ! "coexp": coefficient in exponential in-canopy wind profile
  ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
  ! canopy and roughness-sublayer U(z) at z=h
  REAL, ALLOCATABLE :: coexp      (:) ! Extinction coef for wind profile in canopy 
  ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
  REAL, ALLOCATABLE :: usuh       (:) ! Friction velocity/windspeed at canopy height 
  REAL, ALLOCATABLE :: term2      (:)  
  REAL, ALLOCATABLE :: term3      (:)  
  REAL, ALLOCATABLE :: term5      (:)  
  REAL, ALLOCATABLE :: term6      (:)  
  REAL, ALLOCATABLE :: term6a     (:) ! for aerodyn resist. calc. 

END TYPE roughness_data_type

TYPE roughness_type

  REAL, POINTER :: disp       (:) ! zero-plane displacement                         
  REAL, POINTER :: hruff      (:) ! canopy height above snow level                  
  REAL, POINTER :: hruff_grmx (:) ! max ht of canopy from tiles on same grid        
  REAL, POINTER :: rt0us      (:) ! eq. 3.54, SCAM manual (CSIRO tech report 132)   
  REAL, POINTER :: rt1usa     (:) ! resistance from disp to hruf                    
  REAL, POINTER :: rt1usb     (:) ! resist fr hruf to zruffs (zref if zref<zruffs)  
  REAL, POINTER :: rt1        (:) ! 1/aerodynamic conductance                       
  REAL, POINTER :: za_uv      (:) ! level of lowest atmospheric model layer         
  REAL, POINTER :: za_tq      (:) ! level of lowest atmospheric model layer         
  REAL, POINTER :: z0m        (:) ! roughness length                                
  REAL, POINTER :: zref_uv    (:) ! Reference height for met forcing                
  REAL, POINTER :: zref_tq    (:) ! Reference height for met forcing                
  REAL, POINTER :: zruffs     (:) ! SCALAR Roughness sublayer depth (ground=origin) 
  REAL, POINTER :: z0soilsn   (:) ! roughness length of bare soil surface           
  REAL, POINTER :: z0soil     (:) ! roughness length of bare soil surface           
  REAL, POINTER :: coexp      (:) ! Extinction coef for wind profile in canopy 
  REAL, POINTER :: usuh       (:) ! Friction velocity/windspeed at canopy height 
  REAL, POINTER :: term2      (:)  
  REAL, POINTER :: term3      (:)  
  REAL, POINTER :: term5      (:)  
  REAL, POINTER :: term6      (:)  
  REAL, POINTER :: term6a     (:) ! for aerodyn resist. calc. 

END TYPE roughness_type

CONTAINS

SUBROUTINE alloc_roughness_type(roughness, mp)

USE grid_constants_mod_cbl,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: niter            ! number of iterations for za/L

IMPLICIT NONE

TYPE(roughness_data_type), INTENT(INOUT) :: roughness
INTEGER, INTENT(IN) :: mp

ALLOCATE( roughness% disp       (mp) )
ALLOCATE( roughness% hruff      (mp) )
ALLOCATE( roughness% hruff_grmx (mp) )
ALLOCATE( roughness% rt0us      (mp) )
ALLOCATE( roughness% rt1usa     (mp) )
ALLOCATE( roughness% rt1usb     (mp) )
ALLOCATE( roughness% rt1        (mp) )
ALLOCATE( roughness% za_uv      (mp) )
ALLOCATE( roughness% za_tq      (mp) )
ALLOCATE( roughness% z0m        (mp) )
ALLOCATE( roughness% zref_uv    (mp) )
ALLOCATE( roughness% zref_tq    (mp) )
ALLOCATE( roughness% zruffs     (mp) )
ALLOCATE( roughness% z0soilsn   (mp) )
ALLOCATE( roughness% z0soil     (mp) )
ALLOCATE( roughness% coexp      (mp) )
ALLOCATE( roughness% usuh       (mp) )
ALLOCATE( roughness% term2      (mp) )
ALLOCATE( roughness% term3      (mp) )
ALLOCATE( roughness% term5      (mp) )
ALLOCATE( roughness% term6      (mp) )
ALLOCATE( roughness% term6a     (mp) )

roughness % disp       (:)      = 0.0      
roughness % hruff      (:)      = 0.0      
roughness % hruff_grmx (:)      = 0.0      
roughness % rt0us      (:)      = 0.0      
roughness % rt1usa     (:)      = 0.0      
roughness % rt1usb     (:)      = 0.0      
roughness % rt1        (:)      = 0.0      
roughness % za_uv      (:)      = 0.0      
roughness % za_tq      (:)      = 0.0      
roughness % z0m        (:)      = 0.0      
roughness % zref_uv    (:)      = 0.0      
roughness % zref_tq    (:)      = 0.0      
roughness % zruffs     (:)      = 0.0      
roughness % z0soilsn   (:)      = 0.0      
roughness % z0soil     (:)      = 0.0      
roughness % coexp      (:)      = 0.0      
roughness % usuh       (:)      = 0.0      
roughness % term2      (:)      = 0.0      
roughness % term3      (:)      = 0.0      
roughness % term5      (:)      = 0.0      
roughness % term6      (:)      = 0.0      
roughness % term6a     (:)      = 0.0      

RETURN
END SUBROUTINE alloc_roughness_type

SUBROUTINE dealloc_roughness_type(roughness)

TYPE(roughness_type), INTENT(inout) :: roughness

DEALLOCATE ( roughness % disp       )
DEALLOCATE ( roughness % hruff      )
DEALLOCATE ( roughness % hruff_grmx )
DEALLOCATE ( roughness % rt0us      )
DEALLOCATE ( roughness % rt1usa     )
DEALLOCATE ( roughness % rt1usb     )
DEALLOCATE ( roughness % rt1        )
DEALLOCATE ( roughness % za_uv      )
DEALLOCATE ( roughness % za_tq      )
DEALLOCATE ( roughness % z0m        )
DEALLOCATE ( roughness % zref_uv    )
DEALLOCATE ( roughness % zref_tq    )
DEALLOCATE ( roughness % zruffs     )
DEALLOCATE ( roughness % z0soilsn   )
DEALLOCATE ( roughness % z0soil     )
DEALLOCATE ( roughness % coexp      )
DEALLOCATE ( roughness % usuh       )
DEALLOCATE ( roughness % term2      )
DEALLOCATE ( roughness % term3      )
DEALLOCATE ( roughness % term5      )
DEALLOCATE ( roughness % term6      )
DEALLOCATE ( roughness % term6a     )

RETURN
END SUBROUTINE dealloc_roughness_type

SUBROUTINE assoc_roughness_type(roughness, roughness_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(roughness_type),      INTENT(IN OUT)         :: roughness
TYPE(roughness_data_type), INTENT(IN OUT), TARGET :: roughness_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_roughness_cbl(roughness)

roughness% disp             => roughness_data% disp       
roughness% hruff            => roughness_data% hruff      
roughness% hruff_grmx       => roughness_data% hruff_grmx 
roughness% rt0us            => roughness_data% rt0us      
roughness% rt1usa           => roughness_data% rt1usa     
roughness% rt1usb           => roughness_data% rt1usb     
roughness% rt1              => roughness_data% rt1        
roughness% za_uv            => roughness_data% za_uv      
roughness% za_tq            => roughness_data% za_tq      
roughness% z0m              => roughness_data% z0m        
roughness% zref_uv          => roughness_data% zref_uv    
roughness% zref_tq          => roughness_data% zref_tq    
roughness% zruffs           => roughness_data% zruffs     
roughness% z0soilsn         => roughness_data% z0soilsn   
roughness% z0soil           => roughness_data% z0soil     
roughness% coexp            => roughness_data% coexp      
roughness% usuh             => roughness_data% usuh       
roughness% term2            => roughness_data% term2      
roughness% term3            => roughness_data% term3      
roughness% term5            => roughness_data% term5      
roughness% term6            => roughness_data% term6      
roughness% term6a           => roughness_data% term6a     

RETURN
END SUBROUTINE assoc_roughness_type

SUBROUTINE nullify_roughness_cbl( roughness )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(roughness_type), INTENT(IN OUT) :: roughness 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( roughness % disp       )
NULLIFY( roughness % hruff      )
NULLIFY( roughness % hruff_grmx )
NULLIFY( roughness % rt0us      )
NULLIFY( roughness % rt1usa     )
NULLIFY( roughness % rt1usb     )
NULLIFY( roughness % rt1        )
NULLIFY( roughness % za_uv      )
NULLIFY( roughness % za_tq      )
NULLIFY( roughness % z0m        )
NULLIFY( roughness % zref_uv    )
NULLIFY( roughness % zref_tq    )
NULLIFY( roughness % zruffs     )
NULLIFY( roughness % z0soilsn   )
NULLIFY( roughness % z0soil     )
NULLIFY( roughness % coexp      )
NULLIFY( roughness % usuh       )
NULLIFY( roughness % term2      )
NULLIFY( roughness % term3      )
NULLIFY( roughness % term5      )
NULLIFY( roughness % term6      )
NULLIFY( roughness % term6a     )

RETURN

END SUBROUTINE nullify_roughness_cbl

END MODULE cable_roughness_type_mod
