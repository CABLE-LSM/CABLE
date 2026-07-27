MODULE cable_air_type_mod

IMPLICIT NONE

PUBLIC :: air_type
PUBLIC :: air_data_type
PUBLIC :: alloc_air_type
PUBLIC :: dealloc_air_type
PUBLIC :: assoc_air_type
PUBLIC :: nullify_air_cbl

TYPE air_data_type

  REAL, ALLOCATABLE :: rho (:)   ! dry air density (kg m-3)
  REAL, ALLOCATABLE :: volm(:)   ! molar volume (m3 mol-1)
  REAL, ALLOCATABLE :: rlam(:)   ! latent heat for water (j/kg)
  REAL, ALLOCATABLE :: qsat(:)   ! saturation specific humidity
  REAL, ALLOCATABLE :: epsi(:)   ! d(qsat)/dt ((kg/kg)/k)
  REAL, ALLOCATABLE :: visc(:)   ! air kinematic viscosity (m2/s)
  REAL, ALLOCATABLE :: psyc(:)   ! psychrometric constant
  REAL, ALLOCATABLE :: dsatdk(:) ! d(es)/dt (mb/k)
  REAL, ALLOCATABLE :: cmolar(:) ! conv. from m/s to mol/m2/s
      
END TYPE air_data_type

TYPE air_type

  REAL, POINTER :: rho (:)   ! dry air density (kg m-3)
  REAL, POINTER :: volm(:)   ! molar volume (m3 mol-1)
  REAL, POINTER :: rlam(:)   ! latent heat for water (j/kg)
  REAL, POINTER :: qsat(:)   ! saturation specific humidity
  REAL, POINTER :: epsi(:)   ! d(qsat)/dt ((kg/kg)/k)
  REAL, POINTER :: visc(:)   ! air kinematic viscosity (m2/s)
  REAL, POINTER :: psyc(:)   ! psychrometric constant
  REAL, POINTER :: dsatdk(:) ! d(es)/dt (mb/k)
  REAL, POINTER :: cmolar(:) ! conv. from m/s to mol/m2/s
      
END TYPE air_type

CONTAINS

SUBROUTINE alloc_air_type(air, mp)
IMPLICIT NONE

TYPE(air_data_type), INTENT(INOUT) :: air
INTEGER, INTENT(IN) :: mp

ALLOCATE( air% rho   (mp) )
ALLOCATE( air% volm  (mp) )
ALLOCATE( air% rlam  (mp) )
ALLOCATE( air% qsat  (mp) )
ALLOCATE( air% epsi  (mp) )
ALLOCATE( air% visc  (mp) )
ALLOCATE( air% psyc  (mp) )
ALLOCATE( air% dsatdk(mp) )
ALLOCATE( air% cmolar(mp) )

air% rho (:)   = 0.0        
air% volm(:)   = 0.0        
air% rlam(:)   = 0.0        
air% qsat(:)   = 0.0        
air% epsi(:)   = 0.0        
air% visc(:)   = 0.0        
air% psyc(:)   = 0.0        
air% dsatdk(:) = 0.0      
air% cmolar(:) = 0.0      

RETURN
END SUBROUTINE alloc_air_type

SUBROUTINE dealloc_air_type(air)
IMPLICIT NONE

TYPE(air_type), INTENT(inout) :: air

DEALLOCATE( air% rho    )
DEALLOCATE( air% volm   )
DEALLOCATE( air% rlam   )
DEALLOCATE( air% qsat   )
DEALLOCATE( air% epsi   )
DEALLOCATE( air% visc   )
DEALLOCATE( air% psyc   )
DEALLOCATE( air% dsatdk )
DEALLOCATE( air% cmolar )

RETURN
END SUBROUTINE dealloc_air_type

SUBROUTINE assoc_air_type(air, air_data )
! Description:
!   Associate the CABLE work pointers in the derived type structure
IMPLICIT NONE

!Arguments
TYPE(air_type),      INTENT(IN OUT)         :: air
TYPE(air_data_type), INTENT(IN OUT), TARGET :: air_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''

!End of header

CALL nullify_air_cbl(air)

air% rho    => air_data% rho 
air% volm   => air_data% volm
air% rlam   => air_data% rlam
air% qsat   => air_data% qsat
air% epsi   => air_data% epsi
air% visc   => air_data% visc
air% psyc   => air_data% psyc
air% dsatdk => air_data% dsatdk
air% cmolar => air_data% cmolar

RETURN
END SUBROUTINE assoc_air_type

SUBROUTINE nullify_air_cbl( air )
! Description:
!   Nullify the CABLE work pointers in the derived type structure
IMPLICIT NONE

!Arguments
TYPE(air_type), INTENT(IN OUT) :: air 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_AIR_VARS_CBL'
!End of header

NULLIFY( air % rho    )
NULLIFY( air % volm   )
NULLIFY( air % rlam   )
NULLIFY( air % qsat   )
NULLIFY( air % epsi   )
NULLIFY( air % visc   )
NULLIFY( air % psyc   )
NULLIFY( air % dsatdk )
NULLIFY( air % cmolar )

RETURN

END SUBROUTINE nullify_air_cbl

END MODULE cable_air_type_mod
