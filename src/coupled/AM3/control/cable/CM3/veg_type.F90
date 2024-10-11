MODULE cable_veg_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: veg_type
PUBLIC :: veg_data_type
PUBLIC :: alloc_veg_type
PUBLIC :: dealloc_veg_type
PUBLIC :: assoc_veg_type
PUBLIC :: nullify_veg_cbl

! Vegetation parameters:
TYPE veg_data_type

  INTEGER, ALLOCATABLE :: iveg (:)      ! vegetation(+nveg) type 
  INTEGER, ALLOCATABLE :: iLU  (:)      ! land use type
	LOGICAL, ALLOCATABLE :: deciduous (:) ! flag used for phenology fix 

	REAL, ALLOCATABLE :: canst1  (:)    ! max intercepted water by canopy (mm/LAI)              
	REAL, ALLOCATABLE :: dleaf   (:)    ! chararacteristc legnth of leaf (m)                
	REAL, ALLOCATABLE :: ejmax   (:)    ! max pot. electron transp rate top leaf(mol/m2/s)      
	REAL, ALLOCATABLE :: meth    (:)    ! method for calculation of canopy fluxes and temp.     
	REAL, ALLOCATABLE :: frac4   (:)    ! fraction of c4 plants                                
	REAL, ALLOCATABLE :: hc      (:)    ! roughness height of canopy (veg - snow)               
	REAL, ALLOCATABLE :: vlai    (:)    ! leaf area index                                       
	REAL, ALLOCATABLE :: xalbnir (:)                                                            
	REAL, ALLOCATABLE :: rp20    (:)    ! plant respiration coefficient at 20 C                 
	REAL, ALLOCATABLE :: rpcoef  (:)    ! temperature coef nonleaf plant respiration (1/C)      
	REAL, ALLOCATABLE :: rs20    (:)    ! soil respiration at 20 C [mol m-2 s-1]                
	REAL, ALLOCATABLE :: shelrb  (:)    ! sheltering factor (dimensionless)                     
	REAL, ALLOCATABLE :: vegcf   (:)    ! kdcorbin, 08/10                                       
	REAL, ALLOCATABLE :: tminvj  (:)    ! min temperature of the start of photosynthesis        
	REAL, ALLOCATABLE :: toptvj  (:)    ! opt temperature of the start of photosynthesis        
	REAL, ALLOCATABLE :: tmaxvj  (:)    ! max temperature of the start of photosynthesis        
	REAL, ALLOCATABLE :: vbeta   (:)    !                                                       
	REAL, ALLOCATABLE :: vcmax   (:)    ! max RuBP carboxylation rate top leaf (mol/m2/s)       
	REAL, ALLOCATABLE :: xfang   (:)    ! leaf angle PARAMETER                                
	REAL, ALLOCATABLE :: extkn   (:)    ! extinction coef for vertical                          
	REAL, ALLOCATABLE :: vlaimax (:)    ! extinction coef for vertical                          
	REAL, ALLOCATABLE :: wai     (:)    ! wood area index (stem+branches+twigs)                
	REAL, ALLOCATABLE :: a1gs    (:)    ! a1 parameter in stomatal conductance model            
	REAL, ALLOCATABLE :: d0gs    (:)    ! d0 in stomatal conductance model                
	REAL, ALLOCATABLE :: alpha   (:)    ! initial slope of J-Q response curve                
	REAL, ALLOCATABLE :: convex  (:)    ! convexity of J-Q response curve                
	REAL, ALLOCATABLE :: cfrd    (:)    ! ratio of day respiration to vcmax                
	REAL, ALLOCATABLE :: gswmin  (:)    ! minimal stomatal conductance                          
	REAL, ALLOCATABLE :: conkc0  (:)    ! Michaelis-menton constant for carboxylase             
	REAL, ALLOCATABLE :: conko0  (:)    ! Michaelis-menton constant for oxygenase               
	REAL, ALLOCATABLE :: ekc     (:)    ! activation energy for caroxylagse                
	REAL, ALLOCATABLE :: eko     (:)    ! acvtivation enegery for oxygenase                
	REAL, ALLOCATABLE :: g0      (:)    ! Belinda's stomatal model intercept, Ticket #56.
	REAL, ALLOCATABLE :: g1      (:)    ! Belinda's stomatal model slope, Ticket #56.   

	REAL, ALLOCATABLE :: refl    (:,:)  
	REAL, ALLOCATABLE :: taul    (:,:)  
	REAL, ALLOCATABLE :: froot   (:,:)  ! fraction of root in each soil layer

  ! Additional  veg parameters:
	REAL(r_2), ALLOCATABLE :: rootbeta (:) ! parameter for estimating vertical root mass distribution (froot) 
	REAL(r_2), ALLOCATABLE :: gamma    (:) ! parameter in root efficiency function (Lai and Katul 2000)       
	REAL(r_2), ALLOCATABLE :: ZR       (:) ! maximum rooting depth (cm)                                       
	REAL(r_2), ALLOCATABLE :: F10      (:) ! fraction of roots in top 10 cm                                   
	REAL(r_2), ALLOCATABLE :: clitt    (:) !                                                                  
  
  ! Additional POP veg param
  INTEGER,   ALLOCATABLE :: disturbance_interval  (:,:)      
	REAL(r_2), ALLOCATABLE :: disturbance_intensity (:,:) !                                                                  

END TYPE veg_data_type

TYPE veg_type

  INTEGER, POINTER :: iveg (:)      ! vegetation(+nveg) type 
  INTEGER, POINTER :: iLU  (:)      ! land use type
	LOGICAL, POINTER :: deciduous (:) ! flag used for phenology fix 

	REAL, POINTER :: canst1  (:)    ! max intercepted water by canopy (mm/LAI)              
	REAL, POINTER :: dleaf   (:)    ! chararacteristc legnth of leaf (m)                
	REAL, POINTER :: ejmax   (:)    ! max pot. electron transp rate top leaf(mol/m2/s)      
	REAL, POINTER :: meth    (:)    ! method for calculation of canopy fluxes and temp.     
	REAL, POINTER :: frac4   (:)    ! fraction of c4 plants                                
	REAL, POINTER :: hc      (:)    ! roughness height of canopy (veg - snow)               
	REAL, POINTER :: vlai    (:)    ! leaf area index                                       
	REAL, POINTER :: xalbnir (:)                                                            
	REAL, POINTER :: rp20    (:)    ! plant respiration coefficient at 20 C                 
	REAL, POINTER :: rpcoef  (:)    ! temperature coef nonleaf plant respiration (1/C)      
	REAL, POINTER :: rs20    (:)    ! soil respiration at 20 C [mol m-2 s-1]                
	REAL, POINTER :: shelrb  (:)    ! sheltering factor (dimensionless)                     
	REAL, POINTER :: vegcf   (:)    ! kdcorbin, 08/10                                       
	REAL, POINTER :: tminvj  (:)    ! min temperature of the start of photosynthesis        
	REAL, POINTER :: toptvj  (:)    ! opt temperature of the start of photosynthesis        
	REAL, POINTER :: tmaxvj  (:)    ! max temperature of the start of photosynthesis        
	REAL, POINTER :: vbeta   (:)    !                                                       
	REAL, POINTER :: vcmax   (:)    ! max RuBP carboxylation rate top leaf (mol/m2/s)       
	REAL, POINTER :: xfang   (:)    ! leaf angle PARAMETER                                
	REAL, POINTER :: extkn   (:)    ! extinction coef for vertical                          
	REAL, POINTER :: vlaimax (:)    ! extinction coef for vertical                          
	REAL, POINTER :: wai     (:)    ! wood area index (stem+branches+twigs)                
	REAL, POINTER :: a1gs    (:)    ! a1 parameter in stomatal conductance model            
	REAL, POINTER :: d0gs    (:)    ! d0 in stomatal conductance model                
	REAL, POINTER :: alpha   (:)    ! initial slope of J-Q response curve                
	REAL, POINTER :: convex  (:)    ! convexity of J-Q response curve                
	REAL, POINTER :: cfrd    (:)    ! ratio of day respiration to vcmax                
	REAL, POINTER :: gswmin  (:)    ! minimal stomatal conductance                          
	REAL, POINTER :: conkc0  (:)    ! Michaelis-menton constant for carboxylase             
	REAL, POINTER :: conko0  (:)    ! Michaelis-menton constant for oxygenase               
	REAL, POINTER :: ekc     (:)    ! activation energy for caroxylagse                
	REAL, POINTER :: eko     (:)    ! acvtivation enegery for oxygenase                
	REAL, POINTER :: g0      (:)    ! Belinda's stomatal model intercept, Ticket #56.
	REAL, POINTER :: g1      (:)    ! Belinda's stomatal model slope, Ticket #56.   

	REAL, POINTER :: refl    (:,:)  
	REAL, POINTER :: taul    (:,:)  
	REAL, POINTER :: froot   (:,:)  ! fraction of root in each soil layer

  ! Additional  veg parameters:
	REAL(r_2), POINTER :: rootbeta (:) ! parameter for estimating vertical root mass distribution (froot) 
	REAL(r_2), POINTER :: gamma    (:) ! parameter in root efficiency function (Lai and Katul 2000)       
	REAL(r_2), POINTER :: ZR       (:) ! maximum rooting depth (cm)                                       
	REAL(r_2), POINTER :: F10      (:) ! fraction of roots in top 10 cm                                   
	REAL(r_2), POINTER :: clitt    (:) !                                                                  
  
  ! Additional POP veg param
  INTEGER,   POINTER :: disturbance_interval  (:,:)      
	REAL(r_2), POINTER :: disturbance_intensity (:,:) !                                                                  

END TYPE veg_type

CONTAINS

SUBROUTINE alloc_veg_type(veg, mp)

USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: swb              ! # Radiation SW bands 

IMPLICIT NONE

TYPE(veg_data_type), INTENT(INOUT) :: veg
INTEGER, INTENT(IN) :: mp

ALLOCATE( veg% iveg       (mp) )
ALLOCATE( veg% iLU        (mp) )
ALLOCATE( veg% deciduous  (mp) )
ALLOCATE( veg% canst1     (mp) )
ALLOCATE( veg% dleaf      (mp) )
ALLOCATE( veg% ejmax      (mp) )
ALLOCATE( veg% meth       (mp) )
ALLOCATE( veg% frac4      (mp) )
ALLOCATE( veg% hc         (mp) )
ALLOCATE( veg% vlai       (mp) )
ALLOCATE( veg% xalbnir    (mp) )
ALLOCATE( veg% rp20       (mp) )
ALLOCATE( veg% rpcoef     (mp) )
ALLOCATE( veg% rs20       (mp) )
ALLOCATE( veg% shelrb     (mp) )
ALLOCATE( veg% vegcf      (mp) )
ALLOCATE( veg% tminvj     (mp) )
ALLOCATE( veg% toptvj     (mp) )
ALLOCATE( veg% tmaxvj     (mp) )
ALLOCATE( veg% vbeta      (mp) )
ALLOCATE( veg% vcmax      (mp) )
ALLOCATE( veg% xfang      (mp) )
ALLOCATE( veg% extkn      (mp) )
ALLOCATE( veg% vlaimax    (mp) )
ALLOCATE( veg% wai        (mp) )
ALLOCATE( veg% a1gs       (mp) )
ALLOCATE( veg% d0gs       (mp) )
ALLOCATE( veg% alpha      (mp) )
ALLOCATE( veg% convex     (mp) )
ALLOCATE( veg% cfrd       (mp) )
ALLOCATE( veg% gswmin     (mp) )
ALLOCATE( veg% conkc0     (mp) )
ALLOCATE( veg% conko0     (mp) )
ALLOCATE( veg% ekc        (mp) )
ALLOCATE( veg% eko        (mp) )
ALLOCATE( veg% g0         (mp) )
ALLOCATE( veg% g1         (mp) )
ALLOCATE( veg% refl       (mp,swb) )
ALLOCATE( veg% taul       (mp,swb) )
ALLOCATE( veg% froot      (mp,nsl) )
ALLOCATE( veg% rootbeta   (mp) )
ALLOCATE( veg% gamma      (mp) )
ALLOCATE( veg% ZR         (mp) )
ALLOCATE( veg% F10        (mp) )
ALLOCATE( veg% clitt      (mp) )
ALLOCATE( veg% disturbance_interval  (mp,2) ) !jhan:2??
ALLOCATE( veg% disturbance_intensity (mp,2) ) !jhan:2??

veg % iveg      (:)      = 0.0      
veg % iLU       (:)      = 0.0      
veg % deciduous (:)      = .FALSE. 
veg % canst1    (:)      = 0.0      
veg % dleaf     (:)      = 0.0      
veg % ejmax     (:)      = 0.0      
veg % meth      (:)      = 0.0      
veg % frac4     (:)      = 0.0      
veg % hc        (:)      = 0.0      
veg % vlai      (:)      = 0.0      
veg % xalbnir   (:)      = 0.0      
veg % rp20      (:)      = 0.0      
veg % rpcoef    (:)      = 0.0      
veg % rs20      (:)      = 0.0      
veg % shelrb    (:)      = 0.0      
veg % vegcf     (:)      = 0.0      
veg % tminvj    (:)      = 0.0      
veg % toptvj    (:)      = 0.0      
veg % tmaxvj    (:)      = 0.0      
veg % vbeta     (:)      = 0.0      
veg % vcmax     (:)      = 0.0      
veg % xfang     (:)      = 0.0      
veg % extkn     (:)      = 0.0      
veg % vlaimax   (:)      = 0.0      
veg % wai       (:)      = 0.0      
veg % a1gs      (:)      = 0.0      
veg % d0gs      (:)      = 0.0      
veg % alpha     (:)      = 0.0      
veg % convex    (:)      = 0.0      
veg % cfrd      (:)      = 0.0      
veg % gswmin    (:)      = 0.0      
veg % conkc0    (:)      = 0.0      
veg % conko0    (:)      = 0.0      
veg % ekc       (:)      = 0.0      
veg % eko       (:)      = 0.0      
veg % g0        (:)      = 0.0      
veg % g1        (:)      = 0.0      
veg % refl      (:,:)      = 0.0      
veg % taul      (:,:)      = 0.0      
veg % froot     (:,:)      = 0.0      
veg % rootbeta  (:)      = 0.0      
veg % gamma     (:)      = 0.0      
veg % ZR        (:)      = 0.0      
veg % F10       (:)      = 0.0      
veg % clitt     (:)      = 0.0      
veg % disturbance_interval  (:,:) = 0.0      
veg % disturbance_intensity (:,:) = 0.0      

RETURN
END SUBROUTINE alloc_veg_type

SUBROUTINE dealloc_veg_type(veg)

TYPE(veg_type), INTENT(inout) :: veg

DEALLOCATE ( veg % iveg      )
DEALLOCATE ( veg % iLU       )
DEALLOCATE ( veg % deciduous )
DEALLOCATE ( veg % canst1    )
DEALLOCATE ( veg % dleaf     )
DEALLOCATE ( veg % ejmax     )
DEALLOCATE ( veg % meth      )
DEALLOCATE ( veg % frac4     )
DEALLOCATE ( veg % hc        )
DEALLOCATE ( veg % vlai      )
DEALLOCATE ( veg % xalbnir   )
DEALLOCATE ( veg % rp20      )
DEALLOCATE ( veg % rpcoef    )
DEALLOCATE ( veg % rs20      )
DEALLOCATE ( veg % shelrb    )
DEALLOCATE ( veg % vegcf     )
DEALLOCATE ( veg % tminvj    )
DEALLOCATE ( veg % toptvj    )
DEALLOCATE ( veg % tmaxvj    )
DEALLOCATE ( veg % vbeta     )
DEALLOCATE ( veg % vcmax     )
DEALLOCATE ( veg % xfang     )
DEALLOCATE ( veg % extkn     )
DEALLOCATE ( veg % vlaimax   )
DEALLOCATE ( veg % wai       )
DEALLOCATE ( veg % a1gs      )
DEALLOCATE ( veg % d0gs      )
DEALLOCATE ( veg % alpha     )
DEALLOCATE ( veg % convex    )
DEALLOCATE ( veg % cfrd      )
DEALLOCATE ( veg % gswmin    )
DEALLOCATE ( veg % conkc0    )
DEALLOCATE ( veg % conko0    )
DEALLOCATE ( veg % ekc       )
DEALLOCATE ( veg % eko       )
DEALLOCATE ( veg % g0        )
DEALLOCATE ( veg % g1        )
DEALLOCATE ( veg % refl      )
DEALLOCATE ( veg % taul      )
DEALLOCATE ( veg % froot     )
DEALLOCATE ( veg % rootbeta  )
DEALLOCATE ( veg % gamma     )
DEALLOCATE ( veg % ZR        )
DEALLOCATE ( veg % F10       )
DEALLOCATE ( veg % clitt     )
DEALLOCATE ( veg % disturbance_interval )
DEALLOCATE ( veg % disturbance_intensity)

RETURN
END SUBROUTINE dealloc_veg_type

SUBROUTINE assoc_veg_type(veg, veg_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(veg_type),      INTENT(IN OUT)         :: veg
TYPE(veg_data_type), INTENT(IN OUT), TARGET :: veg_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''
!End of header

CALL nullify_veg_cbl(veg)

veg% iveg            => veg_data% iveg      
veg% iLU             => veg_data% iLU       
veg% deciduous       => veg_data% deciduous 
veg% canst1          => veg_data% canst1    
veg% dleaf           => veg_data% dleaf     
veg% ejmax           => veg_data% ejmax     
veg% meth            => veg_data% meth      
veg% frac4           => veg_data% frac4     
veg% hc              => veg_data% hc        
veg% vlai            => veg_data% vlai      
veg% xalbnir         => veg_data% xalbnir   
veg% rp20            => veg_data% rp20      
veg% rpcoef          => veg_data% rpcoef    
veg% rs20            => veg_data% rs20      
veg% shelrb          => veg_data% shelrb    
veg% vegcf           => veg_data% vegcf     
veg% tminvj          => veg_data% tminvj    
veg% toptvj          => veg_data% toptvj    
veg% tmaxvj          => veg_data% tmaxvj    
veg% vbeta           => veg_data% vbeta     
veg% vcmax           => veg_data% vcmax     
veg% xfang           => veg_data% xfang     
veg% extkn           => veg_data% extkn     
veg% vlaimax         => veg_data% vlaimax   
veg% wai             => veg_data% wai       
veg% a1gs            => veg_data% a1gs      
veg% d0gs            => veg_data% d0gs      
veg% alpha           => veg_data% alpha     
veg% convex          => veg_data% convex    
veg% cfrd            => veg_data% cfrd      
veg% gswmin          => veg_data% gswmin    
veg% conkc0          => veg_data% conkc0    
veg% conko0          => veg_data% conko0    
veg% ekc             => veg_data% ekc       
veg% eko             => veg_data% eko       
veg% g0              => veg_data% g0        
veg% g1              => veg_data% g1        
veg% refl            => veg_data% refl      
veg% taul            => veg_data% taul      
veg% froot           => veg_data% froot     
veg% rootbeta        => veg_data% rootbeta  
veg% gamma           => veg_data% gamma     
veg% ZR              => veg_data% ZR        
veg% F10             => veg_data% F10       
veg% clitt           => veg_data% clitt     
veg% disturbance_interval  => veg_data% disturbance_interval 
veg% disturbance_intensity => veg_data% disturbance_intensity

RETURN
END SUBROUTINE assoc_veg_type

SUBROUTINE nullify_veg_cbl( veg )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(veg_type), INTENT(IN OUT) :: veg 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_CBL_TYPES'
!End of header

NULLIFY( veg % iveg      )
NULLIFY( veg % iLU       )
NULLIFY( veg % deciduous )
NULLIFY( veg % canst1    )
NULLIFY( veg % dleaf     )
NULLIFY( veg % ejmax     )
NULLIFY( veg % meth      )
NULLIFY( veg % frac4     )
NULLIFY( veg % hc        )
NULLIFY( veg % vlai      )
NULLIFY( veg % xalbnir   )
NULLIFY( veg % rp20      )
NULLIFY( veg % rpcoef    )
NULLIFY( veg % rs20      )
NULLIFY( veg % shelrb    )
NULLIFY( veg % vegcf     )
NULLIFY( veg % tminvj    )
NULLIFY( veg % toptvj    )
NULLIFY( veg % tmaxvj    )
NULLIFY( veg % vbeta     )
NULLIFY( veg % vcmax     )
NULLIFY( veg % xfang     )
NULLIFY( veg % extkn     )
NULLIFY( veg % vlaimax   )
NULLIFY( veg % wai       )
NULLIFY( veg % a1gs      )
NULLIFY( veg % d0gs      )
NULLIFY( veg % alpha     )
NULLIFY( veg % convex    )
NULLIFY( veg % cfrd      )
NULLIFY( veg % gswmin    )
NULLIFY( veg % conkc0    )
NULLIFY( veg % conko0    )
NULLIFY( veg % ekc       )
NULLIFY( veg % eko       )
NULLIFY( veg % g0        )
NULLIFY( veg % g1        )
NULLIFY( veg % refl      )
NULLIFY( veg % taul      )
NULLIFY( veg % froot     )
NULLIFY( veg % rootbeta  )
NULLIFY( veg % gamma     )
NULLIFY( veg % ZR        )
NULLIFY( veg % F10       )
NULLIFY( veg % clitt     )
NULLIFY( veg % disturbance_interval  )
NULLIFY( veg % disturbance_intensity )

RETURN

END SUBROUTINE nullify_veg_cbl

END MODULE cable_veg_type_mod
