MODULE cable_bgc_pool_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: bgc_pool_type
PUBLIC :: bgc_pool_data_type
PUBLIC :: alloc_bgc_pool_type
PUBLIC :: dealloc_bgc_pool_type
PUBLIC :: assoc_bgc_pool_type
PUBLIC :: nullify_bgc_pool_cbl

TYPE bgc_pool_data_type

  REAL, ALLOCATABLE :: ratecp (:)   ! plant carbon rate constant (1/year)
  REAL, ALLOCATABLE :: ratecs (:)   ! soil carbon rate constant (1/year)
  REAL, ALLOCATABLE :: cplant (:,:) ! plant carbon (g C/m2)) 
  REAL, ALLOCATABLE :: csoil  (:,:) ! soil carbon (g C/m2)

END TYPE bgc_pool_data_type

TYPE bgc_pool_type

  REAL, POINTER :: ratecp (:)   ! plant carbon rate constant (1/year)
  REAL, POINTER :: ratecs (:)   ! soil carbon rate constant (1/year)
  REAL, POINTER :: cplant (:,:) ! plant carbon (g C/m2)) 
  REAL, POINTER :: csoil  (:,:) ! soil carbon (g C/m2)

END TYPE bgc_pool_type

CONTAINS

SUBROUTINE alloc_bgc_pool_type(bgc_pool, mp)

USE grid_constants_mod_cbl,   ONLY: nsCs  ! # soil carbon stores
USE grid_constants_mod_cbl,   ONLY: nvCs  ! # vegetation carbon stores

IMPLICIT NONE

TYPE(bgc_pool_data_type), INTENT(INOUT) :: bgc_pool
INTEGER, INTENT(IN) :: mp

ALLOCATE( bgc_pool% ratecp (nvCs) )
ALLOCATE( bgc_pool% ratecs (nsCs) )
ALLOCATE( bgc_pool% cplant (mp,nvCs) )
ALLOCATE( bgc_pool% csoil  (mp,nsCs) )

bgc_pool %  ratecp(:)      = 0.0      
bgc_pool %  ratecs(:)      = 0.0      
bgc_pool %  cplant(:,:)    = 0.0      
bgc_pool %  csoil (:,:)    = 0.0      

RETURN
END SUBROUTINE alloc_bgc_pool_type

SUBROUTINE dealloc_bgc_pool_type(bgc_pool)

TYPE(bgc_pool_type), INTENT(inout) :: bgc_pool

DEALLOCATE ( bgc_pool % ratecp)
DEALLOCATE ( bgc_pool % ratecs)
DEALLOCATE ( bgc_pool % cplant)
DEALLOCATE ( bgc_pool % csoil )

RETURN
END SUBROUTINE dealloc_bgc_pool_type

SUBROUTINE assoc_bgc_pool_type(bgc_pool, bgc_pool_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(bgc_pool_type),      INTENT(IN OUT)         :: bgc_pool
TYPE(bgc_pool_data_type), INTENT(IN OUT), TARGET :: bgc_pool_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_bgc_pool_cbl(bgc_pool)

bgc_pool% ratecp      => bgc_pool_data% ratecp
bgc_pool% ratecs      => bgc_pool_data% ratecs
bgc_pool% cplant      => bgc_pool_data% cplant
bgc_pool% csoil       => bgc_pool_data% csoil 

RETURN
END SUBROUTINE assoc_bgc_pool_type

SUBROUTINE nullify_bgc_pool_cbl( bgc_pool )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(bgc_pool_type), INTENT(IN OUT) :: bgc_pool 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( bgc_pool % ratecp )
NULLIFY( bgc_pool % ratecs )
NULLIFY( bgc_pool % cplant )
NULLIFY( bgc_pool % csoil  )

RETURN

END SUBROUTINE nullify_bgc_pool_cbl

END MODULE cable_bgc_pool_type_mod
