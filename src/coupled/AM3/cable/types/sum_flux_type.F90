MODULE cable_sum_flux_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: sum_flux_type
PUBLIC :: sum_flux_data_type
PUBLIC :: alloc_sum_flux_type
PUBLIC :: dealloc_sum_flux_type
PUBLIC :: assoc_sum_flux_type
PUBLIC :: nullify_sum_flux_cbl

! Cumulative flux variables:
TYPE sum_flux_data_type

  REAL, ALLOCATABLE :: sumpn  (:)   ! sum of canopy photosynthesis (g C m-2)
  REAL, ALLOCATABLE :: sumrp  (:)   ! sum of plant respiration (g C m-2)    
  REAL, ALLOCATABLE :: sumrpw (:)   ! sum of plant respiration (g C m-2)     
  REAL, ALLOCATABLE :: sumrpr (:)   ! sum of plant respiration (g C m-2)     
  REAL, ALLOCATABLE :: sumrs  (:)   ! sum of soil respiration (g C m-2)      
  REAL, ALLOCATABLE :: sumrd  (:)   ! sum of daytime respiration (g C m-2)   
  REAL, ALLOCATABLE :: dsumpn (:)   ! daily sumpn                           
  REAL, ALLOCATABLE :: dsumrp (:)   ! daily sumrp                           
  REAL, ALLOCATABLE :: dsumrs (:)   ! daily sumrs                           
  REAL, ALLOCATABLE :: dsumrd (:)   ! daily sumrd                           
  REAL, ALLOCATABLE :: sumxrp (:)   ! sum plant resp. modifier                 
  REAL, ALLOCATABLE :: sumxrs (:)   ! sum soil resp. modifier                  

END TYPE sum_flux_data_type

TYPE sum_flux_type

  REAL, POINTER :: sumpn  (:)   ! sum of canopy photosynthesis (g C m-2)
  REAL, POINTER :: sumrp  (:)   ! sum of plant respiration (g C m-2)    
  REAL, POINTER :: sumrpw (:)   ! sum of plant respiration (g C m-2)     
  REAL, POINTER :: sumrpr (:)   ! sum of plant respiration (g C m-2)     
  REAL, POINTER :: sumrs  (:)   ! sum of soil respiration (g C m-2)      
  REAL, POINTER :: sumrd  (:)   ! sum of daytime respiration (g C m-2)   
  REAL, POINTER :: dsumpn (:)   ! daily sumpn                           
  REAL, POINTER :: dsumrp (:)   ! daily sumrp                           
  REAL, POINTER :: dsumrs (:)   ! daily sumrs                           
  REAL, POINTER :: dsumrd (:)   ! daily sumrd                           
  REAL, POINTER :: sumxrp (:)   ! sum plant resp. modifier                 
  REAL, POINTER :: sumxrs (:)   ! sum soil resp. modifier                  

END TYPE sum_flux_type

CONTAINS

SUBROUTINE alloc_sum_flux_type(sum_flux, mp)

USE grid_constants_mod_cbl,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: niter            ! number of iterations for za/L

IMPLICIT NONE

TYPE(sum_flux_data_type), INTENT(INOUT) :: sum_flux
INTEGER, INTENT(IN) :: mp

ALLOCATE( sum_flux% sumpn  (mp) )
ALLOCATE( sum_flux% sumrp  (mp) )
ALLOCATE( sum_flux% sumrpw (mp) )
ALLOCATE( sum_flux% sumrpr (mp) )
ALLOCATE( sum_flux% sumrs  (mp) )
ALLOCATE( sum_flux% sumrd  (mp) )
ALLOCATE( sum_flux% dsumpn (mp) )
ALLOCATE( sum_flux% dsumrp (mp) )
ALLOCATE( sum_flux% dsumrs (mp) )
ALLOCATE( sum_flux% dsumrd (mp) )
ALLOCATE( sum_flux% sumxrp (mp) )
ALLOCATE( sum_flux% sumxrs (mp) )

sum_flux % sumpn  (:) = 0.0      
sum_flux % sumrp  (:) = 0.0      
sum_flux % sumrpw (:) = 0.0      
sum_flux % sumrpr (:) = 0.0      
sum_flux % sumrs  (:) = 0.0      
sum_flux % sumrd  (:) = 0.0      
sum_flux % dsumpn (:) = 0.0      
sum_flux % dsumrp (:) = 0.0      
sum_flux % dsumrs (:) = 0.0      
sum_flux % dsumrd (:) = 0.0      
sum_flux % sumxrp (:) = 0.0      
sum_flux % sumxrs (:) = 0.0      

RETURN
END SUBROUTINE alloc_sum_flux_type

SUBROUTINE dealloc_sum_flux_type(sum_flux)

TYPE(sum_flux_type), INTENT(inout) :: sum_flux

DEALLOCATE ( sum_flux % sumpn  )
DEALLOCATE ( sum_flux % sumrp  )
DEALLOCATE ( sum_flux % sumrpw )
DEALLOCATE ( sum_flux % sumrpr )
DEALLOCATE ( sum_flux % sumrs  )
DEALLOCATE ( sum_flux % sumrd  )
DEALLOCATE ( sum_flux % dsumpn )
DEALLOCATE ( sum_flux % dsumrp )
DEALLOCATE ( sum_flux % dsumrs )
DEALLOCATE ( sum_flux % dsumrd )
DEALLOCATE ( sum_flux % sumxrp )
DEALLOCATE ( sum_flux % sumxrs )

RETURN
END SUBROUTINE dealloc_sum_flux_type

SUBROUTINE assoc_sum_flux_type(sum_flux, sum_flux_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(sum_flux_type),      INTENT(IN OUT)         :: sum_flux
TYPE(sum_flux_data_type), INTENT(IN OUT), TARGET :: sum_flux_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_sum_flux_cbl(sum_flux)

sum_flux% sumpn        => sum_flux_data% sumpn  
sum_flux% sumrp        => sum_flux_data% sumrp  
sum_flux% sumrpw       => sum_flux_data% sumrpw 
sum_flux% sumrpr       => sum_flux_data% sumrpr 
sum_flux% sumrs        => sum_flux_data% sumrs  
sum_flux% sumrd        => sum_flux_data% sumrd  
sum_flux% dsumpn       => sum_flux_data% dsumpn 
sum_flux% dsumrp       => sum_flux_data% dsumrp 
sum_flux% dsumrs       => sum_flux_data% dsumrs 
sum_flux% dsumrd       => sum_flux_data% dsumrd 
sum_flux% sumxrp       => sum_flux_data% sumxrp 
sum_flux% sumxrs       => sum_flux_data% sumxrs 

RETURN
END SUBROUTINE assoc_sum_flux_type

SUBROUTINE nullify_sum_flux_cbl( sum_flux )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(sum_flux_type), INTENT(IN OUT) :: sum_flux 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( sum_flux % sumpn  )
NULLIFY( sum_flux % sumrp  )
NULLIFY( sum_flux % sumrpw )
NULLIFY( sum_flux % sumrpr )
NULLIFY( sum_flux % sumrs  )
NULLIFY( sum_flux % sumrd  )
NULLIFY( sum_flux % dsumpn )
NULLIFY( sum_flux % dsumrp )
NULLIFY( sum_flux % dsumrs )
NULLIFY( sum_flux % dsumrd )
NULLIFY( sum_flux % sumxrp )
NULLIFY( sum_flux % sumxrs )

RETURN

END SUBROUTINE nullify_sum_flux_cbl

END MODULE cable_sum_flux_type_mod
