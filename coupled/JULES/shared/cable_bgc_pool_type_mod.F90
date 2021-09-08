MODULE cable_bgc_pool_type_mod

USE cable_types_mod, ONLY: ncp, ncs

IMPLICIT NONE
  
TYPE bgc_pool_type

  REAL, DIMENSION(:,:), POINTER ::                                            &
       cplant,   & ! plant carbon (g C/m2))
       csoil      ! soil carbon (g C/m2)
     
  !H!REAL, DIMENSION(:,:), POINTER ::                                         &
  REAL, DIMENSION(:), POINTER ::                                              &
       ratecp,                                                                &
       ratecs

END TYPE bgc_pool_type

!Instantiation:
TYPE(bgc_pool_type) :: bgc_cbl

END MODULE cable_bgc_pool_type_mod
