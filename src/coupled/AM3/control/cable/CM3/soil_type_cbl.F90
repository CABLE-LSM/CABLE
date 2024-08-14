MODULE cable_soil_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: soil_type
PUBLIC :: soil_data_type
PUBLIC :: alloc_soil_type
PUBLIC :: dealloc_soil_type
PUBLIC :: assoc_soil_type
PUBLIC :: nullify_soil_cbl

! Soil parameters:
TYPE soil_data_type

  INTEGER, ALLOCATABLE :: isoilm (:)  ! integer soil type
  REAL, ALLOCATABLE :: bch       (:)  ! parameter b in Campbell equation                  
  REAL, ALLOCATABLE :: c3        (:)  ! c3 drainage coeff (fraction)                      
  REAL, ALLOCATABLE :: clay      (:)  ! fraction of soil which is clay                    
  REAL, ALLOCATABLE :: css       (:)  ! soil specific heat capacity [kJ/kg/K]             
  REAL, ALLOCATABLE :: hsbh      (:)  ! difsat * etasat (=hyds*abs(sucs)*bch)             
  REAL, ALLOCATABLE :: hyds      (:)  ! hydraulic cond@ saturation [m/s], Ksat   
  REAL, ALLOCATABLE :: i2bp3     (:)  ! par. one in K vis suction (=nint(bch)+2)
  REAL, ALLOCATABLE :: ibp2      (:)  ! par. two in K vis suction (fn of pbch)            
  REAL, ALLOCATABLE :: rhosoil   (:)  ! soil density [kg/m3]                              
  REAL, ALLOCATABLE :: sand      (:)  ! fraction of soil which is sand                    
  REAL, ALLOCATABLE :: sfc       (:)  ! vol H2O @ field capacity                          
  REAL, ALLOCATABLE :: silt      (:)  ! fraction of soil which is silt                    
  REAL, ALLOCATABLE :: ssat      (:)  ! vol H2O @ saturation                              
  REAL, ALLOCATABLE :: sucs      (:)  ! suction at saturation (m)                         
  REAL, ALLOCATABLE :: swilt     (:)  ! vol H2O @ wilting                                 
  REAL, ALLOCATABLE :: zse       (:)  ! soil layer thickness (1=top) [m]         
  REAL, ALLOCATABLE :: zshh      (:)  ! dist b/n consecutive layer midpoints (m)
  REAL, ALLOCATABLE :: soilcol   (:)  ! color per patches/tiles Ticket #27
  REAL, ALLOCATABLE :: albsoilf  (:)  ! soil reflectance Ticket #27
  REAL, ALLOCATABLE :: albsoil   (:,:)! soil reflectance

  REAL(r_2), ALLOCATABLE :: heat_cap_lower_limit (:,:)  
  REAL(r_2), ALLOCATABLE :: zse_vec         (:,:)  
  REAL(r_2), ALLOCATABLE :: css_vec         (:,:)  
  REAL(r_2), ALLOCATABLE :: cnsd_vec        (:,:)  
  REAL(r_2), ALLOCATABLE :: cnsd            (:) ! thermal cond dry soil [W/m/K] 
  REAL(r_2), ALLOCATABLE :: pwb_min         (:) ! working var (swilt/ssat)**ibp2      

  REAL(r_2), ALLOCATABLE :: drain_dens        (:) ! mean dist to rivers/streams
  REAL(r_2), ALLOCATABLE :: elev              (:) ! elevation above sea level                        
  REAL(r_2), ALLOCATABLE :: elev_std          (:) ! elevation above sea level                        
  REAL(r_2), ALLOCATABLE :: slope             (:) ! mean slope of grid cell                          
  REAL(r_2), ALLOCATABLE :: slope_std         (:) ! stddev of grid cell slope                       

  ! Parameters for GW module that vary with soil layer
  REAL(r_2), ALLOCATABLE :: sucs_vec    (:,:) ! psi at saturation in [mm]                        
  REAL(r_2), ALLOCATABLE :: hyds_vec    (:,:) ! sat hydraulic cond [mm/s]         
  REAL(r_2), ALLOCATABLE :: bch_vec     (:,:) ! C and H B [none]                                 
  REAL(r_2), ALLOCATABLE :: clay_vec    (:,:) ! fraction of soil that is clay              
  REAL(r_2), ALLOCATABLE :: sand_vec    (:,:) ! fraction of soil that is sand              
  REAL(r_2), ALLOCATABLE :: silt_vec    (:,:) ! fraction of soil that is silt              
  REAL(r_2), ALLOCATABLE :: org_vec     (:,:) ! frac soil made of organic soils
  REAL(r_2), ALLOCATABLE :: rhosoil_vec (:,:) ! soil density  [kg/m3]
  REAL(r_2), ALLOCATABLE :: ssat_vec    (:,:) ! vol H2O content at sat [mm3/mm3] 
  REAL(r_2), ALLOCATABLE :: watr        (:,:) ! resid soil H2O content [mm3/mm3]     
  REAL(r_2), ALLOCATABLE :: sfc_vec     (:,:) ! field capcacity (hk = 1 mm/day)                  
  REAL(r_2), ALLOCATABLE :: swilt_vec   (:,:) ! wilting point (hk = 0.02 mm/day)                 

  ! Parameters for GW module for the aquifer
  REAL(r_2), ALLOCATABLE :: GWsucs_vec   (:) ! head in the aquifer [mm]                               
  REAL(r_2), ALLOCATABLE :: GWhyds_vec   (:) ! satur hydraulic cond [mm/s] 
  REAL(r_2), ALLOCATABLE :: GWbch_vec    (:) ! clapp and horn b [none]               
  REAL(r_2), ALLOCATABLE :: GWssat_vec   (:) ! saturated water content [mm3/mm3]       
  REAL(r_2), ALLOCATABLE :: GWwatr       (:) ! residual water content [mm3/mm3]        
  REAL(r_2), ALLOCATABLE :: GWz          (:) ! node depth of the aquifer    [m]                       
  REAL(r_2), ALLOCATABLE :: GWdz         (:) ! thickness of the aquifer   [m]                         
  REAL(r_2), ALLOCATABLE :: GWrhosoil_vec(:) ! density of substrate [kg/m3]               
  
  ! Additional SLI parameters
  INTEGER(r_2), ALLOCATABLE :: nhorizons (:)  ! number of soil horizons
  INTEGER(r_2), ALLOCATABLE :: ishorizon (:,:)  ! horizon number 1:nhorizons
  REAL(r_2),    ALLOCATABLE :: clitt     (:)  ! litter (tC/ha)
  REAL(r_2),    ALLOCATABLE :: zeta      (:)  ! macropore parameter
  REAL(r_2),    ALLOCATABLE :: fsatmax   (:)  ! variably saturated area parameter
  !REAL(r_2), DIMENSION(:,:), POINTER :: swilt_vec ! vol H2O @ wilting
  !REAL(r_2), DIMENSION(:,:), POINTER :: ssat_vec  ! vol H2O @ sat
  !REAL(r_2), DIMENSION(:,:), POINTER :: sfc_vec   ! vol H2O @ fc

END TYPE soil_data_type

TYPE soil_type

  INTEGER, POINTER :: isoilm (:)   ! integer soil type
  REAL, POINTER :: bch       (:)   ! parameter b in Campbell equation                  
  REAL, POINTER :: c3        (:)   ! c3 drainage coeff (fraction)                      
  REAL, POINTER :: clay      (:)   ! fraction of soil which is clay                    
  REAL, POINTER :: css       (:)   ! soil specific heat capacity [kJ/kg/K]             
  REAL, POINTER :: hsbh      (:)   ! difsat * etasat (=hyds*abs(sucs)*bch)             
  REAL, POINTER :: hyds      (:)   ! hydraulic cond@ saturation [m/s], Ksat   
  REAL, POINTER :: i2bp3     (:)   ! par. one in K vis suction (=nint(bch)+2)          
  REAL, POINTER :: ibp2      (:)   ! par. two in K vis suction (fn of pbch)            
  REAL, POINTER :: rhosoil   (:)   ! soil density [kg/m3]                              
  REAL, POINTER :: sand      (:)   ! fraction of soil which is sand                    
  REAL, POINTER :: sfc       (:)   ! vol H2O @ field capacity                          
  REAL, POINTER :: silt      (:)   ! fraction of soil which is silt                    
  REAL, POINTER :: ssat      (:)   ! vol H2O @ saturation                              
  REAL, POINTER :: sucs      (:)   ! suction at saturation (m)                         
  REAL, POINTER :: swilt     (:)   ! vol H2O @ wilting                                 
  REAL, POINTER :: zse       (:)   ! thickness of each soil layer (1=top) [m]         
  REAL, POINTER :: zshh      (:)   ! dist b/n consecutive layer midpoints (m)
  REAL, POINTER :: soilcol   (:)   ! color per patches/tiles Ticket #27
  REAL, POINTER :: albsoilf  (:)   ! soil reflectance Ticket #27
  REAL, POINTER :: albsoil   (:,:) ! soil reflectance 

  REAL(r_2), POINTER :: heat_cap_lower_limit (:,:)  
  REAL(r_2), POINTER :: zse_vec        (:,:)  
  REAL(r_2), POINTER :: css_vec        (:,:)  
  REAL(r_2), POINTER :: cnsd_vec       (:,:)  
  REAL(r_2), POINTER :: cnsd           (:) ! thermal cond dry soil [W/m/K] 
  REAL(r_2), POINTER :: pwb_min        (:) ! working var (swilt/ssat)**ibp2      

  REAL(r_2), POINTER :: drain_dens        (:) ! mean dist to rivers/streams) 
  REAL(r_2), POINTER :: elev              (:) ! elevation above sea level                        
  REAL(r_2), POINTER :: elev_std          (:) ! elevation above sea level                        
  REAL(r_2), POINTER :: slope             (:) ! mean slope of grid cell                          
  REAL(r_2), POINTER :: slope_std         (:) ! stddev of grid cell slope                       

  ! Parameters for GW module that vary with soil layer
  REAL(r_2), POINTER :: sucs_vec    (:,:) ! psi at saturation in [mm]                        
  REAL(r_2), POINTER :: hyds_vec    (:,:) ! sat hydraulic cond [mm/s]         
  REAL(r_2), POINTER :: bch_vec     (:,:) ! C and H B [none]                                 
  REAL(r_2), POINTER :: clay_vec    (:,:) ! fraction of soil that is clay              
  REAL(r_2), POINTER :: sand_vec    (:,:) ! fraction of soil that is sand              
  REAL(r_2), POINTER :: silt_vec    (:,:) ! fraction of soil that is silt              
  REAL(r_2), POINTER :: org_vec     (:,:) ! frac soil made of organic soils
  REAL(r_2), POINTER :: rhosoil_vec (:,:) ! soil density  [kg/m3]
  REAL(r_2), POINTER :: ssat_vec    (:,:) ! vol H2O content at sat [mm3/mm3] 
  REAL(r_2), POINTER :: watr        (:,:) ! resid soil H2O content [mm3/mm3]     
  REAL(r_2), POINTER :: sfc_vec     (:,:) ! field capcacity (hk = 1 mm/day)                  
  REAL(r_2), POINTER :: swilt_vec   (:,:) ! wilting point (hk = 0.02 mm/day)                 

  ! Parameters for GW module for the aquifer
  REAL(r_2), POINTER :: GWsucs_vec   (:) ! head in the aquifer [mm]                               
  REAL(r_2), POINTER :: GWhyds_vec   (:) ! satur hydraulic cond [mm/s] 
  REAL(r_2), POINTER :: GWbch_vec    (:) ! clapp and horn b [none]               
  REAL(r_2), POINTER :: GWssat_vec   (:) ! saturated water content [mm3/mm3]       
  REAL(r_2), POINTER :: GWwatr       (:) ! residual water content [mm3/mm3]        
  REAL(r_2), POINTER :: GWz          (:) ! node depth of the aquifer    [m]                       
  REAL(r_2), POINTER :: GWdz         (:) ! thickness of the aquifer   [m]                         
  REAL(r_2), POINTER :: GWrhosoil_vec(:) ! density of substrate [kg/m3]               
  
  ! Additional SLI parameters
  INTEGER(r_2), POINTER :: nhorizons (:)  ! number of soil horizons
  INTEGER(r_2), POINTER :: ishorizon (:,:)  ! horizon number 1:nhorizons
  REAL(r_2),    POINTER :: clitt     (:)  ! litter (tC/ha)
  REAL(r_2),    POINTER :: zeta      (:)  ! macropore parameter
  REAL(r_2),    POINTER :: fsatmax   (:)  ! variably saturated area parameter
  !REAL(r_2), DIMENSION(:,:), POINTER :: swilt_vec ! vol H2O @ wilting
  !REAL(r_2), DIMENSION(:,:), POINTER :: ssat_vec  ! vol H2O @ sat
  !REAL(r_2), DIMENSION(:,:), POINTER :: sfc_vec   ! vol H2O @ fc

END TYPE soil_type

CONTAINS

SUBROUTINE alloc_soil_type(soil, mp)

USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: swb              ! # SW bands 

IMPLICIT NONE

TYPE(soil_data_type), INTENT(INOUT) :: soil
INTEGER, INTENT(IN) :: mp

ALLOCATE( soil% isoilm   (mp) )
ALLOCATE( soil% bch      (mp) )
ALLOCATE( soil% c3       (mp) )
ALLOCATE( soil% clay     (mp) )
ALLOCATE( soil% css      (mp) )
ALLOCATE( soil% hsbh     (mp) )
ALLOCATE( soil% hyds     (mp) )
ALLOCATE( soil% i2bp3    (mp) )
ALLOCATE( soil% ibp2     (mp) )
ALLOCATE( soil% rhosoil  (mp) )
ALLOCATE( soil% sand     (mp) )
ALLOCATE( soil% sfc      (mp) )
ALLOCATE( soil% silt     (mp) )
ALLOCATE( soil% ssat     (mp) )
ALLOCATE( soil% sucs     (mp) )
ALLOCATE( soil% swilt    (mp) )
ALLOCATE( soil% zse      (nsl) )
ALLOCATE( soil% zshh     (nsl+1) )
ALLOCATE( soil% soilcol  (mp) )
ALLOCATE( soil% albsoilf (mp) )
ALLOCATE( soil% albsoil  (mp,swb) )

ALLOCATE( soil% heat_cap_lower_limit (mp,nsl) )
ALLOCATE( soil% zse_vec  (mp,nsl) )
ALLOCATE( soil% css_vec  (mp,nsl) )
ALLOCATE( soil% cnsd_vec (mp,nsl) )
ALLOCATE( soil% cnsd     (mp) )
ALLOCATE( soil% pwb_min  (mp) )

ALLOCATE( soil% drain_dens (mp) ) ! mean dist to rivers/streams
ALLOCATE( soil% elev       (mp) ) ! elevation above sea level  
ALLOCATE( soil% elev_std   (mp) ) ! elevation above sea level  
ALLOCATE( soil% slope      (mp) ) ! mean slope of grid cell    
ALLOCATE( soil% slope_std  (mp) ) ! stddev of grid cell slope  

! Parameters for GW module that vary with soil layer
ALLOCATE( soil% sucs_vec    (mp,nsl) )
ALLOCATE( soil% hyds_vec    (mp,nsl) )
ALLOCATE( soil% bch_vec     (mp,nsl) )
ALLOCATE( soil% clay_vec    (mp,nsl) )
ALLOCATE( soil% sand_vec    (mp,nsl) )
ALLOCATE( soil% silt_vec    (mp,nsl) )
ALLOCATE( soil% org_vec     (mp,nsl) )
ALLOCATE( soil% rhosoil_vec (mp,nsl) )
ALLOCATE( soil% ssat_vec    (mp,nsl) )
ALLOCATE( soil% watr        (mp,nsl) )
ALLOCATE( soil% sfc_vec     (mp,nsl) )
ALLOCATE( soil% swilt_vec   (mp,nsl) )

! Parameters for GW module for the aquifer
ALLOCATE( soil% GWhyds_vec   (mp) )
ALLOCATE( soil% GWsucs_vec   (mp) )
ALLOCATE( soil% GWbch_vec    (mp) )
ALLOCATE( soil% GWssat_vec   (mp) )
ALLOCATE( soil% GWwatr       (mp) )
ALLOCATE( soil% GWz          (mp) )
ALLOCATE( soil% GWdz         (mp) )
ALLOCATE( soil% GWrhosoil_vec(mp) )
  
! Additional SLI parameters
ALLOCATE( soil% nhorizons(mp) )
ALLOCATE( soil% ishorizon(mp,nsl) )
ALLOCATE( soil% clitt    (mp) )
ALLOCATE( soil% zeta     (mp) )
ALLOCATE( soil% fsatmax  (mp) )

soil % isoilm         (:) = 0.0      
soil % bch            (:) = 0.0      
soil % c3             (:) = 0.0      
soil % clay           (:) = 0.0      
soil % css            (:) = 0.0      
soil % hsbh           (:) = 0.0      
soil % hyds           (:) = 0.0      
soil % i2bp3          (:) = 0.0      
soil % ibp2           (:) = 0.0      
soil % rhosoil        (:) = 0.0      
soil % sand           (:) = 0.0      
soil % sfc            (:) = 0.0      
soil % silt           (:) = 0.0      
soil % ssat           (:) = 0.0      
soil % sucs           (:) = 0.0      
soil % swilt          (:) = 0.0      
soil % zse            (:) = 0.0      
soil % zshh           (:) = 0.0      
soil % soilcol        (:) = 0.0      
soil % albsoilf       (:) = 0.0      
soil % albsoil        (:,:) = 0.0      
soil % zse_vec        (:,:) = 0.0      
soil % css_vec        (:,:) = 0.0      
soil % cnsd_vec       (:,:) = 0.0      
soil % cnsd           (:) = 0.0      
soil % pwb_min        (:) = 0.0      
soil % drain_dens     (:) = 0.0      
soil % elev           (:) = 0.0      
soil % elev_std       (:) = 0.0      
soil % slope          (:) = 0.0      
soil % slope_std      (:) = 0.0      
soil % sucs_vec       (:,:) = 0.0      
soil % hyds_vec       (:,:) = 0.0      
soil % bch_vec        (:,:) = 0.0      
soil % clay_vec       (:,:) = 0.0      
soil % sand_vec       (:,:) = 0.0      
soil % silt_vec       (:,:) = 0.0      
soil % org_vec        (:,:) = 0.0      
soil % rhosoil_vec    (:,:) = 0.0      
soil % ssat_vec       (:,:) = 0.0      
soil % watr           (:,:) = 0.0      
soil % sfc_vec        (:,:) = 0.0      
soil % swilt_vec      (:,:) = 0.0      
soil % GWhyds_vec     (:) = 0.0      
soil % GWsucs_vec     (:) = 0.0      
soil % GWbch_vec      (:) = 0.0      
soil % GWssat_vec     (:) = 0.0      
soil % GWwatr         (:) = 0.0      
soil % GWz            (:) = 0.0      
soil % GWdz           (:) = 0.0      
soil % GWrhosoil_vec  (:) = 0.0      
soil % nhorizons      (:) = 0.0      
soil % ishorizon      (:,:) = 0.0      
soil % clitt          (:) = 0.0      
soil % zeta           (:) = 0.0      
soil % fsatmax        (:) = 0.0      
soil % heat_cap_lower_limit(:,:) = 0.0      

RETURN
END SUBROUTINE alloc_soil_type

SUBROUTINE dealloc_soil_type(soil)

TYPE(soil_type), INTENT(inout) :: soil

DEALLOCATE ( soil % isoilm        )
DEALLOCATE ( soil % bch           )
DEALLOCATE ( soil % c3            )
DEALLOCATE ( soil % clay          )
DEALLOCATE ( soil % css           )
DEALLOCATE ( soil % hsbh          )
DEALLOCATE ( soil % hyds          )
DEALLOCATE ( soil % i2bp3         )
DEALLOCATE ( soil % ibp2          )
DEALLOCATE ( soil % rhosoil       )
DEALLOCATE ( soil % sand          )
DEALLOCATE ( soil % sfc           )
DEALLOCATE ( soil % silt          )
DEALLOCATE ( soil % ssat          )
DEALLOCATE ( soil % sucs          )
DEALLOCATE ( soil % swilt         )
DEALLOCATE ( soil % zse           )
DEALLOCATE ( soil % zshh          )
DEALLOCATE ( soil % soilcol       )
DEALLOCATE ( soil % albsoilf      )
DEALLOCATE ( soil % albsoil       )
DEALLOCATE ( soil % heat_cap_lower_limit )
DEALLOCATE ( soil % zse_vec       )
DEALLOCATE ( soil % css_vec       )
DEALLOCATE ( soil % cnsd_vec      )
DEALLOCATE ( soil % cnsd          )
DEALLOCATE ( soil % pwb_min       )
DEALLOCATE ( soil % drain_dens    )
DEALLOCATE ( soil % elev          )
DEALLOCATE ( soil % elev_std      )
DEALLOCATE ( soil % slope         )
DEALLOCATE ( soil % slope_std     )
DEALLOCATE ( soil % sucs_vec      )
DEALLOCATE ( soil % hyds_vec      )
DEALLOCATE ( soil % bch_vec       )
DEALLOCATE ( soil % clay_vec      )
DEALLOCATE ( soil % sand_vec      )
DEALLOCATE ( soil % silt_vec      )
DEALLOCATE ( soil % org_vec       )
DEALLOCATE ( soil % rhosoil_vec   )
DEALLOCATE ( soil % ssat_vec      )
DEALLOCATE ( soil % watr          )
DEALLOCATE ( soil % sfc_vec       )
DEALLOCATE ( soil % swilt_vec     )
DEALLOCATE ( soil % GWhyds_vec    )
DEALLOCATE ( soil % GWsucs_vec    )
DEALLOCATE ( soil % GWbch_vec     )
DEALLOCATE ( soil % GWssat_vec    )
DEALLOCATE ( soil % GWwatr        )
DEALLOCATE ( soil % GWz           )
DEALLOCATE ( soil % GWdz          )
DEALLOCATE ( soil % GWrhosoil_vec )
DEALLOCATE ( soil % nhorizons     )
DEALLOCATE ( soil % ishorizon     )
DEALLOCATE ( soil % clitt         )
DEALLOCATE ( soil % zeta          )
DEALLOCATE ( soil % fsatmax       )

RETURN
END SUBROUTINE dealloc_soil_type

SUBROUTINE assoc_soil_type(soil, soil_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(soil_type),      INTENT(IN OUT)         :: soil
TYPE(soil_data_type), INTENT(IN OUT), TARGET :: soil_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_soil_cbl(soil)

soil% isoilm        => soil_data% isoilm         
soil% bch           => soil_data% bch            
soil% c3            => soil_data% c3             
soil% clay          => soil_data% clay           
soil% css           => soil_data% css            
soil% hsbh          => soil_data% hsbh           
soil% hyds          => soil_data% hyds           
soil% i2bp3         => soil_data% i2bp3          
soil% ibp2          => soil_data% ibp2           
soil% rhosoil       => soil_data% rhosoil        
soil% sand          => soil_data% sand           
soil% sfc           => soil_data% sfc            
soil% silt          => soil_data% silt           
soil% ssat          => soil_data% ssat           
soil% sucs          => soil_data% sucs           
soil% swilt         => soil_data% swilt          
soil% zse           => soil_data% zse            
soil% zshh          => soil_data% zshh           
soil% soilcol       => soil_data% soilcol        
soil% albsoilf      => soil_data% albsoilf       
soil% albsoil       => soil_data% albsoil        
soil% heat_cap_lower_limit => soil_data% heat_cap_lower_limit  
soil% zse_vec       => soil_data% zse_vec        
soil% css_vec       => soil_data% css_vec        
soil% cnsd_vec      => soil_data% cnsd_vec       
soil% cnsd          => soil_data% cnsd           
soil% pwb_min       => soil_data% pwb_min        
soil% drain_dens    => soil_data% drain_dens     
soil% elev          => soil_data% elev           
soil% elev_std      => soil_data% elev_std       
soil% slope         => soil_data% slope          
soil% slope_std     => soil_data% slope_std      
soil% sucs_vec      => soil_data% sucs_vec       
soil% hyds_vec      => soil_data% hyds_vec       
soil% bch_vec       => soil_data% bch_vec        
soil% clay_vec      => soil_data% clay_vec       
soil% sand_vec      => soil_data% sand_vec       
soil% silt_vec      => soil_data% silt_vec       
soil% org_vec       => soil_data% org_vec        
soil% rhosoil_vec   => soil_data% rhosoil_vec    
soil% ssat_vec      => soil_data% ssat_vec       
soil% watr          => soil_data% watr           
soil% sfc_vec       => soil_data% sfc_vec        
soil% swilt_vec     => soil_data% swilt_vec      
soil% GWhyds_vec    => soil_data% GWhyds_vec     
soil% GWsucs_vec    => soil_data% GWsucs_vec     
soil% GWbch_vec     => soil_data% GWbch_vec      
soil% GWssat_vec    => soil_data% GWssat_vec     
soil% GWwatr        => soil_data% GWwatr         
soil% GWz           => soil_data% GWz            
soil% GWdz          => soil_data% GWdz           
soil% GWrhosoil_vec => soil_data% GWrhosoil_vec  
soil% nhorizons     => soil_data% nhorizons      
soil% ishorizon     => soil_data% ishorizon      
soil% clitt         => soil_data% clitt          
soil% zeta          => soil_data% zeta           
soil% fsatmax       => soil_data% fsatmax        

RETURN
END SUBROUTINE assoc_soil_type

SUBROUTINE nullify_soil_cbl( soil )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(soil_type), INTENT(IN OUT) :: soil 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( soil %  isoilm        )
NULLIFY( soil %  bch           )
NULLIFY( soil %  c3            )
NULLIFY( soil %  clay          )
NULLIFY( soil %  css           )
NULLIFY( soil %  hsbh          )
NULLIFY( soil %  hyds          )
NULLIFY( soil %  i2bp3         )
NULLIFY( soil %  ibp2          )
NULLIFY( soil %  rhosoil       )
NULLIFY( soil %  sand          )
NULLIFY( soil %  sfc           )
NULLIFY( soil %  silt          )
NULLIFY( soil %  ssat          )
NULLIFY( soil %  sucs          )
NULLIFY( soil %  swilt         )
NULLIFY( soil %  zse           )
NULLIFY( soil %  zshh          )
NULLIFY( soil %  soilcol       )
NULLIFY( soil %  albsoilf      )
NULLIFY( soil %  albsoil       )
NULLIFY( soil %  heat_cap_lower_limit )
NULLIFY( soil %  zse_vec       )
NULLIFY( soil %  css_vec       )
NULLIFY( soil %  cnsd_vec      )
NULLIFY( soil %  cnsd          )
NULLIFY( soil %  pwb_min       )
NULLIFY( soil %  drain_dens    )
NULLIFY( soil %  elev          )
NULLIFY( soil %  elev_std      )
NULLIFY( soil %  slope         )
NULLIFY( soil %  slope_std     )
NULLIFY( soil %  sucs_vec      )
NULLIFY( soil %  hyds_vec      )
NULLIFY( soil %  bch_vec       )
NULLIFY( soil %  clay_vec      )
NULLIFY( soil %  sand_vec      )
NULLIFY( soil %  silt_vec      )
NULLIFY( soil %  org_vec       )
NULLIFY( soil %  rhosoil_vec   )
NULLIFY( soil %  ssat_vec      )
NULLIFY( soil %  watr          )
NULLIFY( soil %  sfc_vec       )
NULLIFY( soil %  swilt_vec     )
NULLIFY( soil %  GWhyds_vec    )
NULLIFY( soil %  GWsucs_vec    )
NULLIFY( soil %  GWbch_vec     )
NULLIFY( soil %  GWssat_vec    )
NULLIFY( soil %  GWwatr        )
NULLIFY( soil %  GWz           )
NULLIFY( soil %  GWdz          )
NULLIFY( soil %  GWrhosoil_vec )
NULLIFY( soil %  nhorizons     )
NULLIFY( soil %  ishorizon     )
NULLIFY( soil %  clitt         )
NULLIFY( soil %  zeta          )
NULLIFY( soil %  fsatmax       )

RETURN

END SUBROUTINE nullify_soil_cbl

END MODULE cable_soil_type_mod
