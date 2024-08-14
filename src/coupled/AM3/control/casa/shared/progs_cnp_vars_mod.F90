MODULE progs_cnp_vars_mod

IMPLICIT NONE

PUBLIC :: progs_cnp_vars_alloc
PUBLIC :: progs_cnp_vars_assoc
PUBLIC :: progs_cnp_vars_data_type
PUBLIC :: progs_cnp_vars_type
PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROGS_CNP_VARS_MOD'
! Prognostic Fields for CASA-CNP to be initialized from IO
TYPE :: progs_cnp_vars_data_type

  REAL, ALLOCATABLE :: C_pool_casa(:,:,:)
  REAL, ALLOCATABLE :: N_pool_casa(:,:,:)
  REAL, ALLOCATABLE :: P_pool_casa(:,:,:)
  REAL, ALLOCATABLE :: soil_order_casa(:)
  REAL, ALLOCATABLE :: N_dep_casa(:)
  REAL, ALLOCATABLE :: N_fix_casa(:)
  REAL, ALLOCATABLE :: P_dust_casa(:)
  REAL, ALLOCATABLE :: P_weath_casa(:)
  REAL, ALLOCATABLE :: LAI_casa(:,:)
  REAL, ALLOCATABLE :: phenphase_casa(:,:)
  REAL, ALLOCATABLE :: wood_hvest_C(:,:,:)
  REAL, ALLOCATABLE :: wood_hvest_N(:,:,:)
  REAL, ALLOCATABLE :: wood_hvest_P(:,:,:)
  REAL, ALLOCATABLE :: thinning(:,:)
  REAL, ALLOCATABLE :: prev_yr_sfrac(:,:)

END TYPE progs_cnp_vars_data_type

TYPE :: progs_cnp_vars_type

  REAL, POINTER, PUBLIC :: C_pool_casa(:,:,:)
  REAL, POINTER, PUBLIC :: N_pool_casa(:,:,:)
  REAL, POINTER, PUBLIC :: P_pool_casa(:,:,:)
  REAL, POINTER, PUBLIC :: soil_order_casa(:)
  REAL, POINTER, PUBLIC :: N_dep_casa(:)         
  REAL, POINTER, PUBLIC :: N_fix_casa(:)         
  REAL, POINTER, PUBLIC :: P_dust_casa(:)         
  REAL, POINTER, PUBLIC :: P_weath_casa(:)         
  REAL, POINTER, PUBLIC :: LAI_casa(:,:)         
  REAL, POINTER, PUBLIC :: phenphase_casa(:,:)         
  REAL, POINTER, PUBLIC :: wood_hvest_C(:,:,:)         
  REAL, POINTER, PUBLIC :: wood_hvest_N(:,:,:)         
  REAL, POINTER, PUBLIC :: wood_hvest_P(:,:,:)         
  REAL, POINTER, PUBLIC :: thinning(:,:)                  
  REAL, POINTER, PUBLIC :: prev_yr_sfrac(:,:)

END TYPE progs_cnp_vars_type

CONTAINS

SUBROUTINE progs_cnp_vars_alloc(land_pts, nsurft, progs_cnp_vars_data )
                                
USE grid_constants_mod_cbl, ONLY : nCpool_casa, nNpool_casa, nPpool_casa
USE casadimension,          ONLY: mwood

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts, nsurft

TYPE(progs_cnp_vars_data_type), INTENT(OUT) :: progs_cnp_vars_data

ALLOCATE ( progs_cnp_vars_data % C_pool_casa     ( land_pts, nsurft, nCpool_casa ) )    
ALLOCATE ( progs_cnp_vars_data % N_pool_casa     ( land_pts, nsurft, nNpool_casa ) )         
ALLOCATE ( progs_cnp_vars_data % P_pool_casa     ( land_pts, nsurft, nPpool_casa ) )         
ALLOCATE ( progs_cnp_vars_data % soil_order_casa ( land_pts                      ) )
ALLOCATE ( progs_cnp_vars_data % N_dep_casa      ( land_pts                      ) )
ALLOCATE ( progs_cnp_vars_data % N_fix_casa      ( land_pts                      ) )
ALLOCATE ( progs_cnp_vars_data % P_dust_casa     ( land_pts                      ) )
ALLOCATE ( progs_cnp_vars_data % P_weath_casa    ( land_pts                      ) )
ALLOCATE ( progs_cnp_vars_data % LAI_casa        ( land_pts, nsurft              ) )
ALLOCATE ( progs_cnp_vars_data % phenphase_casa  ( land_pts, nsurft              ) )
ALLOCATE ( progs_cnp_vars_data % wood_hvest_C    ( land_pts, nsurft, mwood       ) )
ALLOCATE ( progs_cnp_vars_data % wood_hvest_N    ( land_pts, nsurft, mwood       ) )      
ALLOCATE ( progs_cnp_vars_data % wood_hvest_P    ( land_pts, nsurft, mwood       ) )      
ALLOCATE ( progs_cnp_vars_data % thinning        ( land_pts, nsurft              ) )
ALLOCATE ( progs_cnp_vars_data % prev_yr_sfrac   ( land_pts, nsurft              ) )

progs_cnp_vars_data % C_pool_casa(:,:,:)  = 0.0
progs_cnp_vars_data % N_pool_casa(:,:,:)  = 0.0
progs_cnp_vars_data % P_pool_casa(:,:,:)  = 0.0
progs_cnp_vars_data % soil_order_casa(:)  = 0.0
progs_cnp_vars_data % N_dep_casa(:)       = 0.0
progs_cnp_vars_data % N_fix_casa(:)       = 0.0  
progs_cnp_vars_data % P_dust_casa(:)      = 0.0  
progs_cnp_vars_data % P_weath_casa(:)     = 0.0  
progs_cnp_vars_data % LAI_casa(:,:)       = 0.0 
progs_cnp_vars_data % phenphase_casa(:,:) = 0.0
progs_cnp_vars_data % wood_hvest_C(:,:,:) = 0.0
progs_cnp_vars_data % wood_hvest_N(:,:,:) = 0.0
progs_cnp_vars_data % wood_hvest_P(:,:,:) = 0.0
progs_cnp_vars_data % thinning(:,:)       = 0.0
progs_cnp_vars_data % prev_yr_sfrac(:,:)  = 0.0

RETURN
END SUBROUTINE progs_cnp_vars_alloc

SUBROUTINE progs_cnp_vars_dealloc(progs_cnp_vars_data )

IMPLICIT NONE

!Arguments
TYPE(progs_cnp_vars_data_type) :: progs_cnp_vars_data

DEALLOCATE ( progs_cnp_vars_data % C_pool_casa     ) 
DEALLOCATE ( progs_cnp_vars_data % N_pool_casa     )  
DEALLOCATE ( progs_cnp_vars_data % P_pool_casa     )  
DEALLOCATE ( progs_cnp_vars_data % soil_order_casa )  
DEALLOCATE ( progs_cnp_vars_data % N_dep_casa      )  
DEALLOCATE ( progs_cnp_vars_data % N_fix_casa      )  
DEALLOCATE ( progs_cnp_vars_data % P_dust_casa     )  
DEALLOCATE ( progs_cnp_vars_data % P_weath_casa    )  
DEALLOCATE ( progs_cnp_vars_data % LAI_casa        )  
DEALLOCATE ( progs_cnp_vars_data % phenphase_casa  )  
DEALLOCATE ( progs_cnp_vars_data % wood_hvest_C    )  
DEALLOCATE ( progs_cnp_vars_data % wood_hvest_N    )  
DEALLOCATE ( progs_cnp_vars_data % wood_hvest_P    )  
DEALLOCATE ( progs_cnp_vars_data % thinning        )  
DEALLOCATE ( progs_cnp_vars_data % prev_yr_sfrac   )  

RETURN
END SUBROUTINE progs_cnp_vars_dealloc

SUBROUTINE progs_cnp_vars_assoc(progs_cnp_vars, progs_cnp_vars_data )

IMPLICIT NONE

!Arguments
TYPE(progs_cnp_vars_type), INTENT(IN OUT) :: progs_cnp_vars
TYPE(progs_cnp_vars_data_type), INTENT(IN OUT), TARGET :: progs_cnp_vars_data

progs_cnp_vars % C_pool_casa      => progs_cnp_vars_data % C_pool_casa
progs_cnp_vars % N_pool_casa      => progs_cnp_vars_data % N_pool_casa
progs_cnp_vars % P_pool_casa      => progs_cnp_vars_data % P_pool_casa
progs_cnp_vars % soil_order_casa  => progs_cnp_vars_data % soil_order_casa
progs_cnp_vars % N_dep_casa       => progs_cnp_vars_data % N_dep_casa
progs_cnp_vars % N_fix_casa       => progs_cnp_vars_data % N_fix_casa
progs_cnp_vars % P_dust_casa      => progs_cnp_vars_data % P_dust_casa
progs_cnp_vars % P_weath_casa     => progs_cnp_vars_data % P_weath_casa
progs_cnp_vars % LAI_casa         => progs_cnp_vars_data % LAI_casa
progs_cnp_vars % phenphase_casa   => progs_cnp_vars_data % phenphase_casa
progs_cnp_vars % wood_hvest_C     => progs_cnp_vars_data % wood_hvest_C
progs_cnp_vars % wood_hvest_N     => progs_cnp_vars_data % wood_hvest_N
progs_cnp_vars % wood_hvest_P     => progs_cnp_vars_data % wood_hvest_P
progs_cnp_vars % thinning         => progs_cnp_vars_data % thinning
progs_cnp_vars % prev_yr_sfrac    => progs_cnp_vars_data % prev_yr_sfrac

RETURN
END SUBROUTINE progs_cnp_vars_assoc

SUBROUTINE progs_cnp_vars_nullify(progs_cnp_vars)

IMPLICIT NONE

!Arguments
TYPE(progs_cnp_vars_type) :: progs_cnp_vars

NULLIFY ( progs_cnp_vars % C_pool_casa     ) 
NULLIFY ( progs_cnp_vars % N_pool_casa     )  
NULLIFY ( progs_cnp_vars % P_pool_casa     )  
NULLIFY ( progs_cnp_vars % soil_order_casa )  
NULLIFY ( progs_cnp_vars % N_dep_casa      )  
NULLIFY ( progs_cnp_vars % N_fix_casa      )  
NULLIFY ( progs_cnp_vars % P_dust_casa     )  
NULLIFY ( progs_cnp_vars % P_weath_casa    )  
NULLIFY ( progs_cnp_vars % LAI_casa        )  
NULLIFY ( progs_cnp_vars % phenphase_casa  )  
NULLIFY ( progs_cnp_vars % wood_hvest_C    )  
NULLIFY ( progs_cnp_vars % wood_hvest_N    )  
NULLIFY ( progs_cnp_vars % wood_hvest_P    )  
NULLIFY ( progs_cnp_vars % thinning        )  
NULLIFY ( progs_cnp_vars % prev_yr_sfrac   )  

RETURN
 
END SUBROUTINE progs_cnp_vars_nullify
 
END MODULE progs_cnp_vars_mod
