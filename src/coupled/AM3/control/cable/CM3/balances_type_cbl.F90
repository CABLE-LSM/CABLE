MODULE cable_balances_type_mod

IMPLICIT NONE

PUBLIC :: balances_type
PUBLIC :: balances_data_type
PUBLIC :: alloc_balances_type
PUBLIC :: dealloc_balances_type
PUBLIC :: assoc_balances_type
PUBLIC :: nullify_balances_cbl

! Energy and water balance variables:
TYPE balances_data_type

  REAL, ALLOCATABLE :: drybal(:)            ! energy balance for dry canopy               
  REAL, ALLOCATABLE :: ebal(:)              ! energy balance per time step (W/m^2)      (:)  
  REAL, ALLOCATABLE :: ebal_tot(:)          ! cumulative energy balance (W/m^2)         (:)  
  REAL, ALLOCATABLE :: ebal_cncheck(:)      ! energy balance consistency check (W/m^2)  (:)  
  REAL, ALLOCATABLE :: ebal_tot_cncheck(:)  ! cumulative energy balance (W/m^2)         (:)  
  REAL, ALLOCATABLE :: ebaltr(:)            ! energy balance per time step (W/m^2)      (:)  
  REAL, ALLOCATABLE :: ebal_tottr(:)        ! cumulative energy balance (W/m^2)         (:)  
  REAL, ALLOCATABLE :: evap_tot(:)          ! cumulative evapotranspiration (mm/dels)   (:)  
  REAL, ALLOCATABLE :: osnowd0(:)           ! snow depth, first time step               (:)  
  REAL, ALLOCATABLE :: precip_tot(:)        ! cumulative precipitation (mm/dels)        (:)  
  REAL, ALLOCATABLE :: rnoff_tot(:)         ! cumulative runoff (mm/dels)               (:)  
  REAL, ALLOCATABLE :: wbal(:)              ! water balance per time step (mm/dels)     (:)  
  REAL, ALLOCATABLE :: wbal_tot(:)          ! cumulative water balance (mm/dels)        (:)  
  REAL, ALLOCATABLE :: wbtot0(:)            ! total soil water (mm), first time step    (:)  
  REAL, ALLOCATABLE :: wetbal(:)            ! energy balance for wet canopy             (:)  
  REAL, ALLOCATABLE :: cansto0(:)           ! canopy water storage (mm)                 (:)  
  REAL, ALLOCATABLE :: owbtot(:)            ! total soil water (mm), first time step    (:)  
  REAL, ALLOCATABLE :: evapc_tot(:)         ! cumulative evapotranspiration (mm/dels)   (:)  
  REAL, ALLOCATABLE :: evaps_tot(:)         ! cumulative evapotranspiration (mm/dels)   (:)  
  REAL, ALLOCATABLE :: rnof1_tot(:)         ! cumulative runoff (mm/dels)               (:)  
  REAL, ALLOCATABLE :: rnof2_tot(:)         ! cumulative runoff (mm/dels)               (:)  
  REAL, ALLOCATABLE :: snowdc_tot(:)        ! cumulative runoff (mm/dels)               (:)  
  REAL, ALLOCATABLE :: wbal_tot1(:)         ! cumulative water balance (mm/dels)        (:)  
  REAL, ALLOCATABLE :: delwc_tot(:)         ! energy balance for wet canopy             (:)  
  REAL, ALLOCATABLE :: qasrf_tot(:)         ! heat advected to the snow by precip.      (:)  
  REAL, ALLOCATABLE :: qfsrf_tot(:)         ! energy of snowpack phase changes          (:)  
  REAL, ALLOCATABLE :: qssrf_tot(:)         ! energy of snowpack phase changes          (:)  
  REAL, ALLOCATABLE :: Radbal    (:)  
  REAL, ALLOCATABLE :: EbalSoil  (:)  
  REAL, ALLOCATABLE :: Ebalveg   (:)  
  REAL, ALLOCATABLE :: Radbalsum (:)  

END TYPE balances_data_type

TYPE balances_type

  REAL, POINTER :: drybal(:)            ! energy balance for dry canopy               
  REAL, POINTER :: ebal(:)              ! energy balance per time step (W/m^2)      (:)  
  REAL, POINTER :: ebal_tot(:)          ! cumulative energy balance (W/m^2)         (:)  
  REAL, POINTER :: ebal_cncheck(:)      ! energy balance consistency check (W/m^2)  (:)  
  REAL, POINTER :: ebal_tot_cncheck(:)  ! cumulative energy balance (W/m^2)         (:)  
  REAL, POINTER :: ebaltr(:)            ! energy balance per time step (W/m^2)      (:)  
  REAL, POINTER :: ebal_tottr(:)        ! cumulative energy balance (W/m^2)         (:)  
  REAL, POINTER :: evap_tot(:)          ! cumulative evapotranspiration (mm/dels)   (:)  
  REAL, POINTER :: osnowd0(:)           ! snow depth, first time step               (:)  
  REAL, POINTER :: precip_tot(:)        ! cumulative precipitation (mm/dels)        (:)  
  REAL, POINTER :: rnoff_tot(:)         ! cumulative runoff (mm/dels)               (:)  
  REAL, POINTER :: wbal(:)              ! water balance per time step (mm/dels)     (:)  
  REAL, POINTER :: wbal_tot(:)          ! cumulative water balance (mm/dels)        (:)  
  REAL, POINTER :: wbtot0(:)            ! total soil water (mm), first time step    (:)  
  REAL, POINTER :: wetbal(:)            ! energy balance for wet canopy             (:)  
  REAL, POINTER :: cansto0(:)           ! canopy water storage (mm)                 (:)  
  REAL, POINTER :: owbtot(:)            ! total soil water (mm), first time step    (:)  
  REAL, POINTER :: evapc_tot(:)         ! cumulative evapotranspiration (mm/dels)   (:)  
  REAL, POINTER :: evaps_tot(:)         ! cumulative evapotranspiration (mm/dels)   (:)  
  REAL, POINTER :: rnof1_tot(:)         ! cumulative runoff (mm/dels)               (:)  
  REAL, POINTER :: rnof2_tot(:)         ! cumulative runoff (mm/dels)               (:)  
  REAL, POINTER :: snowdc_tot(:)        ! cumulative runoff (mm/dels)               (:)  
  REAL, POINTER :: wbal_tot1(:)         ! cumulative water balance (mm/dels)        (:)  
  REAL, POINTER :: delwc_tot(:)         ! energy balance for wet canopy             (:)  
  REAL, POINTER :: qasrf_tot(:)         ! heat advected to the snow by precip.      (:)  
  REAL, POINTER :: qfsrf_tot(:)         ! energy of snowpack phase changes          (:)  
  REAL, POINTER :: qssrf_tot(:)         ! energy of snowpack phase changes          (:)  
  REAL, POINTER :: Radbal    (:)  
  REAL, POINTER :: EbalSoil  (:)  
  REAL, POINTER :: Ebalveg   (:)  
  REAL, POINTER :: Radbalsum (:)  

END TYPE balances_type

CONTAINS

SUBROUTINE alloc_balances_type(balances, mp)
IMPLICIT NONE

TYPE(balances_data_type), INTENT(INOUT) :: balances
INTEGER, INTENT(IN) :: mp

ALLOCATE( balances% drybal           (mp) )
ALLOCATE( balances% ebal             (mp) )
ALLOCATE( balances% ebal_tot         (mp) )
ALLOCATE( balances% ebal_cncheck     (mp) )
ALLOCATE( balances% ebal_tot_cncheck (mp) )
ALLOCATE( balances% ebaltr           (mp) )
ALLOCATE( balances% ebal_tottr       (mp) )
ALLOCATE( balances% evap_tot         (mp) )
ALLOCATE( balances% osnowd0          (mp) )
ALLOCATE( balances% precip_tot       (mp) )
ALLOCATE( balances% rnoff_tot        (mp) )
ALLOCATE( balances% wbal             (mp) )
ALLOCATE( balances% wbal_tot         (mp) )
ALLOCATE( balances% wbtot0           (mp) )
ALLOCATE( balances% wetbal           (mp) )
ALLOCATE( balances% cansto0          (mp) )
ALLOCATE( balances% owbtot           (mp) )
ALLOCATE( balances% evapc_tot        (mp) )
ALLOCATE( balances% evaps_tot        (mp) )
ALLOCATE( balances% rnof1_tot        (mp) )
ALLOCATE( balances% rnof2_tot        (mp) )
ALLOCATE( balances% snowdc_tot       (mp) )
ALLOCATE( balances% wbal_tot1        (mp) )
ALLOCATE( balances% delwc_tot        (mp) )
ALLOCATE( balances% qasrf_tot        (mp) )
ALLOCATE( balances% qfsrf_tot        (mp) )
ALLOCATE( balances% qssrf_tot        (mp) )
ALLOCATE( balances% Radbal           (mp) )
ALLOCATE( balances% EbalSoil         (mp) )
ALLOCATE( balances% Ebalveg          (mp) )
ALLOCATE( balances% Radbalsum        (mp) )

balances % drybal           (:)      = 0.0      
balances % ebal             (:)      = 0.0      
balances % ebal_tot         (:)      = 0.0      
balances % ebal_cncheck     (:)      = 0.0      
balances % ebal_tot_cncheck (:)      = 0.0      
balances % ebaltr           (:)      = 0.0      
balances % ebal_tottr       (:)      = 0.0      
balances % evap_tot         (:)      = 0.0      
balances % osnowd0          (:)      = 0.0      
balances % precip_tot       (:)      = 0.0      
balances % rnoff_tot        (:)      = 0.0      
balances % wbal             (:)      = 0.0      
balances % wbal_tot         (:)      = 0.0      
balances % wbtot0           (:)      = 0.0      
balances % wetbal           (:)      = 0.0      
balances % cansto0          (:)      = 0.0      
balances % owbtot           (:)      = 0.0      
balances % evapc_tot        (:)      = 0.0      
balances % evaps_tot        (:)      = 0.0      
balances % rnof1_tot        (:)      = 0.0      
balances % rnof2_tot        (:)      = 0.0      
balances % snowdc_tot       (:)      = 0.0      
balances % wbal_tot1        (:)      = 0.0      
balances % delwc_tot        (:)      = 0.0      
balances % qasrf_tot        (:)      = 0.0      
balances % qfsrf_tot        (:)      = 0.0      
balances % qssrf_tot        (:)      = 0.0      
balances % Radbal           (:)      = 0.0      
balances % EbalSoil         (:)      = 0.0      
balances % Ebalveg          (:)      = 0.0      
balances % Radbalsum        (:)      = 0.0      

RETURN
END SUBROUTINE alloc_balances_type

SUBROUTINE dealloc_balances_type(balances)
IMPLICIT NONE

TYPE(balances_type), INTENT(inout) :: balances

DEALLOCATE( balances% drybal           )
DEALLOCATE( balances% ebal             )
DEALLOCATE( balances% ebal_tot         )
DEALLOCATE( balances% ebal_cncheck     )
DEALLOCATE( balances% ebal_tot_cncheck )
DEALLOCATE( balances% ebaltr           )
DEALLOCATE( balances% ebal_tottr       )
DEALLOCATE( balances% evap_tot         )
DEALLOCATE( balances% osnowd0          )
DEALLOCATE( balances% precip_tot       )
DEALLOCATE( balances% rnoff_tot        )
DEALLOCATE( balances% wbal             )
DEALLOCATE( balances% wbal_tot         )
DEALLOCATE( balances% wbtot0           )
DEALLOCATE( balances% wetbal           )
DEALLOCATE( balances% cansto0          )
DEALLOCATE( balances% owbtot           )
DEALLOCATE( balances% evapc_tot        )
DEALLOCATE( balances% evaps_tot        )
DEALLOCATE( balances% rnof1_tot        )
DEALLOCATE( balances% rnof2_tot        )
DEALLOCATE( balances% snowdc_tot       )
DEALLOCATE( balances% wbal_tot1        )
DEALLOCATE( balances% delwc_tot        )
DEALLOCATE( balances% qasrf_tot        )
DEALLOCATE( balances% qfsrf_tot        )
DEALLOCATE( balances% qssrf_tot        )
DEALLOCATE( balances% Radbal           )
DEALLOCATE( balances% EbalSoil         )
DEALLOCATE( balances% Ebalveg          )
DEALLOCATE( balances% Radbalsum        )

RETURN
END SUBROUTINE dealloc_balances_type

SUBROUTINE assoc_balances_type(balances, balances_data )
! Description:
!   Associate the CABLE work pointers in the derived type structure
IMPLICIT NONE

!Arguments
TYPE(balances_type),      INTENT(IN OUT)         :: balances
TYPE(balances_data_type), INTENT(IN OUT), TARGET :: balances_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_balances_cbl(balances)

balances% drybal                 => balances_data% drybal           
balances% ebal                   => balances_data% ebal             
balances% ebal_tot               => balances_data% ebal_tot         
balances% ebal_cncheck           => balances_data% ebal_cncheck     
balances% ebal_tot_cncheck       => balances_data% ebal_tot_cncheck 
balances% ebaltr                 => balances_data% ebaltr           
balances% ebal_tottr             => balances_data% ebal_tottr       
balances% evap_tot               => balances_data% evap_tot         
balances% osnowd0                => balances_data% osnowd0          
balances% precip_tot             => balances_data% precip_tot       
balances% rnoff_tot              => balances_data% rnoff_tot        
balances% wbal                   => balances_data% wbal             
balances% wbal_tot               => balances_data% wbal_tot         
balances% wbtot0                 => balances_data% wbtot0           
balances% wetbal                 => balances_data% wetbal           
balances% cansto0                => balances_data% cansto0          
balances% owbtot                 => balances_data% owbtot           
balances% evapc_tot              => balances_data% evapc_tot        
balances% evaps_tot              => balances_data% evaps_tot        
balances% rnof1_tot              => balances_data% rnof1_tot        
balances% rnof2_tot              => balances_data% rnof2_tot        
balances% snowdc_tot             => balances_data% snowdc_tot       
balances% wbal_tot1              => balances_data% wbal_tot1        
balances% delwc_tot              => balances_data% delwc_tot        
balances% qasrf_tot              => balances_data% qasrf_tot        
balances% qfsrf_tot              => balances_data% qfsrf_tot        
balances% qssrf_tot              => balances_data% qssrf_tot        
balances% Radbal                 => balances_data% Radbal           
balances% EbalSoil               => balances_data% EbalSoil         
balances% Ebalveg                => balances_data% Ebalveg          
balances% Radbalsum              => balances_data% Radbalsum        

RETURN
END SUBROUTINE assoc_balances_type

SUBROUTINE nullify_balances_cbl( balances )
! Description:
!   Nullify the CABLE work pointers in the derived type structure
IMPLICIT NONE

!Arguments
TYPE(balances_type), INTENT(IN OUT) :: balances 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( balances % drybal          )
NULLIFY( balances % ebal            )
NULLIFY( balances % ebal_tot        )
NULLIFY( balances % ebal_cncheck    )
NULLIFY( balances % ebal_tot_cncheck)
NULLIFY( balances % ebaltr          )
NULLIFY( balances % ebal_tottr      )
NULLIFY( balances % evap_tot        )
NULLIFY( balances % osnowd0         )
NULLIFY( balances % precip_tot      )
NULLIFY( balances % rnoff_tot       )
NULLIFY( balances % wbal            )
NULLIFY( balances % wbal_tot        )
NULLIFY( balances % wbtot0          )
NULLIFY( balances % wetbal          )
NULLIFY( balances % cansto0         )
NULLIFY( balances % owbtot          )
NULLIFY( balances % evapc_tot       )
NULLIFY( balances % evaps_tot       )
NULLIFY( balances % rnof1_tot       )
NULLIFY( balances % rnof2_tot       )
NULLIFY( balances % snowdc_tot      )
NULLIFY( balances % wbal_tot1       )
NULLIFY( balances % delwc_tot       )
NULLIFY( balances % qasrf_tot       )
NULLIFY( balances % qfsrf_tot       )
NULLIFY( balances % qssrf_tot       )
NULLIFY( balances % Radbal          )
NULLIFY( balances % EbalSoil        )
NULLIFY( balances % Ebalveg         )
NULLIFY( balances % Radbalsum       )
 
RETURN
 
END SUBROUTINE nullify_balances_cbl
 
END MODULE cable_balances_type_mod
 
 
 
 
 
 
 
 
