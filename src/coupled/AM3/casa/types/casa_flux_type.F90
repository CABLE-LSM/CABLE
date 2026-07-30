MODULE casa_flux_type_mod

USE cable_other_constants_mod, ONLY: r_2 ! currently DOUBLE precision

IMPLICIT NONE

PUBLIC :: casa_flux_data_type
PUBLIC :: casa_flux_type
PUBLIC :: alloc_casa_flux_data_type
PUBLIC :: assoc_casa_flux_type

TYPE casa_flux_data_type

  REAL(r_2), ALLOCATABLE :: Cgpp          (:)
  REAL(r_2), ALLOCATABLE :: Cnpp          (:)
  REAL(r_2), ALLOCATABLE :: Cnbp          (:)
  REAL(r_2), ALLOCATABLE :: Crp           (:)
  REAL(r_2), ALLOCATABLE :: Crgplant      (:)
  REAL(r_2), ALLOCATABLE :: Nminfix       (:)
  REAL(r_2), ALLOCATABLE :: Nminuptake    (:)
  REAL(r_2), ALLOCATABLE :: Plabuptake    (:)
  REAL(r_2), ALLOCATABLE :: Clabloss      (:) 
  REAL(r_2), ALLOCATABLE :: fracClabile   (:)
  REAL(r_2), ALLOCATABLE :: stemnpp       (:) 
  REAL(r_2), ALLOCATABLE :: frac_sapwood  (:) 
  REAL(r_2), ALLOCATABLE :: Cplant_turnover_tot  (:) 
  REAL(r_2), ALLOCATABLE :: sapwood_area  (:)   
  REAL(r_2), ALLOCATABLE :: Crsoil        (:)
  REAL(r_2), ALLOCATABLE :: Nmindep       (:)
  REAL(r_2), ALLOCATABLE :: Nminloss      (:)
  REAL(r_2), ALLOCATABLE :: Nminleach     (:)
  REAL(r_2), ALLOCATABLE :: Nupland       (:)
  REAL(r_2), ALLOCATABLE :: Nlittermin    (:)
  REAL(r_2), ALLOCATABLE :: Nsmin         (:)
  REAL(r_2), ALLOCATABLE :: Nsimm         (:)
  REAL(r_2), ALLOCATABLE :: Nsnet         (:)
  REAL(r_2), ALLOCATABLE :: fNminloss     (:)
  REAL(r_2), ALLOCATABLE :: fNminleach    (:)
  REAL(r_2), ALLOCATABLE :: Pdep          (:)
  REAL(r_2), ALLOCATABLE :: Pwea          (:)
  REAL(r_2), ALLOCATABLE :: Pleach        (:)
  REAL(r_2), ALLOCATABLE :: Ploss         (:)
  REAL(r_2), ALLOCATABLE :: Pupland       (:)
  REAL(r_2), ALLOCATABLE :: Plittermin    (:)
  REAL(r_2), ALLOCATABLE :: Psmin         (:)
  REAL(r_2), ALLOCATABLE :: Psimm         (:)
  REAL(r_2), ALLOCATABLE :: Psnet         (:)
  REAL(r_2), ALLOCATABLE :: fPleach       (:)
  REAL(r_2), ALLOCATABLE :: kplab         (:)
  REAL(r_2), ALLOCATABLE :: kpsorb        (:)
  REAL(r_2), ALLOCATABLE :: kpocc         (:)
  REAL(r_2), ALLOCATABLE :: kmlabp        (:)
  REAL(r_2), ALLOCATABLE :: Psorbmax      (:)
  REAL(r_2), ALLOCATABLE :: Cnep          (:)   
  REAL(r_2), ALLOCATABLE :: FluxCtoCO2    (:)  
  REAL(r_2), ALLOCATABLE :: FluxCtohwp    (:)  
  REAL(r_2), ALLOCATABLE :: FluxNtohwp    (:)  
  REAL(r_2), ALLOCATABLE :: FluxPtohwp    (:)  
  REAL(r_2), ALLOCATABLE :: FluxCtoclear  (:)    
  REAL(r_2), ALLOCATABLE :: FluxNtoclear  (:)    
  REAL(r_2), ALLOCATABLE :: FluxPtoclear  (:)    
  REAL(r_2), ALLOCATABLE :: CtransferLUC  (:)    
  REAL(r_2), ALLOCATABLE :: meangpp       (:)     
  REAL(r_2), ALLOCATABLE :: meanrleaf     (:)         
  REAL(r_2), ALLOCATABLE :: fracCalloc    (:,:)     
  REAL(r_2), ALLOCATABLE :: fracNalloc    (:,:)     
  REAL(r_2), ALLOCATABLE :: fracPalloc    (:,:)     
  REAL(r_2), ALLOCATABLE :: Crmplant      (:,:)  
  REAL(r_2), ALLOCATABLE :: kplant        (:,:)
  REAL(r_2), ALLOCATABLE :: klitter       (:,:)
  REAL(r_2), ALLOCATABLE :: ksoil         (:,:)
  REAL(r_2), ALLOCATABLE :: fromLtoCO2    (:,:)
  REAL(r_2), ALLOCATABLE :: fromStoCO2    (:,:)
  REAL(r_2), ALLOCATABLE :: FluxCtolitter (:,:)
  REAL(r_2), ALLOCATABLE :: FluxNtolitter (:,:)
  REAL(r_2), ALLOCATABLE :: FluxPtolitter (:,:)
  REAL(r_2), ALLOCATABLE :: FluxCtosoil   (:,:)
  REAL(r_2), ALLOCATABLE :: FluxNtosoil   (:,:)
  REAL(r_2), ALLOCATABLE :: FluxPtosoil   (:,:)
  REAL(r_2), ALLOCATABLE :: fromPtoL      (:,:,:)
  REAL(r_2), ALLOCATABLE :: fromLtoS      (:,:,:)
  REAL(r_2), ALLOCATABLE :: fromStoS      (:,:,:)       
  REAL(r_2), ALLOCATABLE :: Cplant_turnover                     (:,:)    
  REAL(r_2), ALLOCATABLE :: Cplant_turnover_disturbance         (:)    
  REAL(r_2), ALLOCATABLE :: Cplant_turnover_crowding            (:)   
  REAL(r_2), ALLOCATABLE :: Cplant_turnover_resource_limitation (:) 

END TYPE casa_flux_data_type

TYPE casa_flux_type

  REAL(r_2), POINTER, PUBLIC :: Cgpp          (:)
  REAL(r_2), POINTER, PUBLIC :: Cnpp          (:)
  REAL(r_2), POINTER, PUBLIC :: Cnbp          (:)
  REAL(r_2), POINTER, PUBLIC :: Crp           (:)
  REAL(r_2), POINTER, PUBLIC :: Crgplant      (:)
  REAL(r_2), POINTER, PUBLIC :: Nminfix       (:)
  REAL(r_2), POINTER, PUBLIC :: Nminuptake    (:)
  REAL(r_2), POINTER, PUBLIC :: Plabuptake    (:)
  REAL(r_2), POINTER, PUBLIC :: Clabloss      (:) 
  REAL(r_2), POINTER, PUBLIC :: fracClabile   (:)
  REAL(r_2), POINTER, PUBLIC :: stemnpp       (:) 
  REAL(r_2), POINTER, PUBLIC :: frac_sapwood  (:) 
  REAL(r_2), POINTER, PUBLIC :: Cplant_turnover_tot  (:) 
  REAL(r_2), POINTER, PUBLIC :: sapwood_area  (:)   
  REAL(r_2), POINTER, PUBLIC :: Crsoil        (:)
  REAL(r_2), POINTER, PUBLIC :: Nmindep       (:)
  REAL(r_2), POINTER, PUBLIC :: Nminloss      (:)
  REAL(r_2), POINTER, PUBLIC :: Nminleach     (:)
  REAL(r_2), POINTER, PUBLIC :: Nupland       (:)
  REAL(r_2), POINTER, PUBLIC :: Nlittermin    (:)
  REAL(r_2), POINTER, PUBLIC :: Nsmin         (:)
  REAL(r_2), POINTER, PUBLIC :: Nsimm         (:)
  REAL(r_2), POINTER, PUBLIC :: Nsnet         (:)
  REAL(r_2), POINTER, PUBLIC :: fNminloss     (:)
  REAL(r_2), POINTER, PUBLIC :: fNminleach    (:)
  REAL(r_2), POINTER, PUBLIC :: Pdep          (:)
  REAL(r_2), POINTER, PUBLIC :: Pwea          (:)
  REAL(r_2), POINTER, PUBLIC :: Pleach        (:)
  REAL(r_2), POINTER, PUBLIC :: Ploss         (:)
  REAL(r_2), POINTER, PUBLIC :: Pupland       (:)
  REAL(r_2), POINTER, PUBLIC :: Plittermin    (:)
  REAL(r_2), POINTER, PUBLIC :: Psmin         (:)
  REAL(r_2), POINTER, PUBLIC :: Psimm         (:)
  REAL(r_2), POINTER, PUBLIC :: Psnet         (:)
  REAL(r_2), POINTER, PUBLIC :: fPleach       (:)
  REAL(r_2), POINTER, PUBLIC :: kplab         (:)
  REAL(r_2), POINTER, PUBLIC :: kpsorb        (:)
  REAL(r_2), POINTER, PUBLIC :: kpocc         (:)
  REAL(r_2), POINTER, PUBLIC :: kmlabp        (:)
  REAL(r_2), POINTER, PUBLIC :: Psorbmax      (:)
  REAL(r_2), POINTER, PUBLIC :: Cnep          (:)   
  REAL(r_2), POINTER, PUBLIC :: FluxCtoCO2    (:)  
  REAL(r_2), POINTER, PUBLIC :: FluxCtohwp    (:)  
  REAL(r_2), POINTER, PUBLIC :: FluxNtohwp    (:)  
  REAL(r_2), POINTER, PUBLIC :: FluxPtohwp    (:)  
  REAL(r_2), POINTER, PUBLIC :: FluxCtoclear  (:)    
  REAL(r_2), POINTER, PUBLIC :: FluxNtoclear  (:)    
  REAL(r_2), POINTER, PUBLIC :: FluxPtoclear  (:)    
  REAL(r_2), POINTER, PUBLIC :: CtransferLUC  (:)    
  REAL(r_2), POINTER, PUBLIC :: meangpp       (:)     
  REAL(r_2), POINTER, PUBLIC :: meanrleaf     (:)         
  REAL(r_2), POINTER, PUBLIC :: fracCalloc    (:,:)     
  REAL(r_2), POINTER, PUBLIC :: fracNalloc    (:,:)     
  REAL(r_2), POINTER, PUBLIC :: fracPalloc    (:,:)     
  REAL(r_2), POINTER, PUBLIC :: Crmplant      (:,:)  
  REAL(r_2), POINTER, PUBLIC :: kplant        (:,:)
  REAL(r_2), POINTER, PUBLIC :: klitter       (:,:)
  REAL(r_2), POINTER, PUBLIC :: ksoil         (:,:)
  REAL(r_2), POINTER, PUBLIC :: fromLtoCO2    (:,:)
  REAL(r_2), POINTER, PUBLIC :: fromStoCO2    (:,:)
  REAL(r_2), POINTER, PUBLIC :: FluxCtolitter (:,:)
  REAL(r_2), POINTER, PUBLIC :: FluxNtolitter (:,:)
  REAL(r_2), POINTER, PUBLIC :: FluxPtolitter (:,:)
  REAL(r_2), POINTER, PUBLIC :: FluxCtosoil   (:,:)
  REAL(r_2), POINTER, PUBLIC :: FluxNtosoil   (:,:)
  REAL(r_2), POINTER, PUBLIC :: FluxPtosoil   (:,:)
  REAL(r_2), POINTER, PUBLIC :: fromPtoL      (:,:,:)
  REAL(r_2), POINTER, PUBLIC :: fromLtoS      (:,:,:)
  REAL(r_2), POINTER, PUBLIC :: fromStoS      (:,:,:)       
  REAL(r_2), POINTER, PUBLIC :: Cplant_turnover                     (:,:)    
  REAL(r_2), POINTER, PUBLIC :: Cplant_turnover_disturbance         (:)    
  REAL(r_2), POINTER, PUBLIC :: Cplant_turnover_crowding            (:)   
  REAL(r_2), POINTER, PUBLIC :: Cplant_turnover_resource_limitation (:) 

END TYPE casa_flux_type

CONTAINS

SUBROUTINE alloc_casa_flux_data_type( casaflux_data, arraysize ) 

USE casadimension, ONLY: mplant, mlitter, msoil

IMPLICIT NONE

TYPE (casa_flux_data_type), INTENT(INOUT) :: casaflux_data 
INTEGER,                    INTENT(IN)    :: arraysize

ALLOCATE ( casaflux_data % Cgpp         (  arraysize ) ) 
ALLOCATE ( casaflux_data % Cnpp         (  arraysize ) ) 
ALLOCATE ( casaflux_data % Cnbp         (  arraysize ) ) 
ALLOCATE ( casaflux_data % Crp          (  arraysize ) ) 
ALLOCATE ( casaflux_data % Crgplant     (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nminfix      (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nminuptake   (  arraysize ) ) 
ALLOCATE ( casaflux_data % Plabuptake   (  arraysize ) ) 
ALLOCATE ( casaflux_data % Clabloss     (  arraysize ) ) 
ALLOCATE ( casaflux_data % fracClabile  (  arraysize ) ) 
ALLOCATE ( casaflux_data % stemnpp      (  arraysize ) ) 
ALLOCATE ( casaflux_data % frac_sapwood (  arraysize ) ) 
ALLOCATE ( casaflux_data % Cplant_turnover_tot (  arraysize ) ) 
ALLOCATE ( casaflux_data % sapwood_area (  arraysize ) ) 
ALLOCATE ( casaflux_data % Crsoil       (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nmindep      (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nminloss     (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nminleach    (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nupland      (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nlittermin   (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nsmin        (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nsimm        (  arraysize ) ) 
ALLOCATE ( casaflux_data % Nsnet        (  arraysize ) ) 
ALLOCATE ( casaflux_data % fNminloss    (  arraysize ) ) 
ALLOCATE ( casaflux_data % fNminleach   (  arraysize ) ) 
ALLOCATE ( casaflux_data % Pdep         (  arraysize ) ) 
ALLOCATE ( casaflux_data % Pwea         (  arraysize ) ) 
ALLOCATE ( casaflux_data % Pleach       (  arraysize ) ) 
ALLOCATE ( casaflux_data % Ploss        (  arraysize ) ) 
ALLOCATE ( casaflux_data % Pupland      (  arraysize ) ) 
ALLOCATE ( casaflux_data % Plittermin   (  arraysize ) ) 
ALLOCATE ( casaflux_data % Psmin        (  arraysize ) ) 
ALLOCATE ( casaflux_data % Psimm        (  arraysize ) ) 
ALLOCATE ( casaflux_data % Psnet        (  arraysize ) ) 
ALLOCATE ( casaflux_data % fPleach      (  arraysize ) ) 
ALLOCATE ( casaflux_data % kplab        (  arraysize ) ) 
ALLOCATE ( casaflux_data % kpsorb       (  arraysize ) ) 
ALLOCATE ( casaflux_data % kpocc        (  arraysize ) ) 
ALLOCATE ( casaflux_data % kmlabp       (  arraysize ) ) 
ALLOCATE ( casaflux_data % Psorbmax     (  arraysize ) ) 
ALLOCATE ( casaflux_data % Cnep         (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxCtoCO2   (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxCtohwp   (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxNtohwp   (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxPtohwp   (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxCtoclear (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxNtoclear (  arraysize ) ) 
ALLOCATE ( casaflux_data % FluxPtoclear (  arraysize ) ) 
ALLOCATE ( casaflux_data % CtransferLUC (  arraysize ) ) 
ALLOCATE ( casaflux_data % meangpp      (  arraysize ) ) 
ALLOCATE ( casaflux_data % meanrleaf    (  arraysize ) ) 
  
ALLOCATE ( casaflux_data % fracCalloc    ( arraysize, mplant          ) )
ALLOCATE ( casaflux_data % fracNalloc    ( arraysize, mplant          ) )
ALLOCATE ( casaflux_data % fracPalloc    ( arraysize, mplant          ) )
ALLOCATE ( casaflux_data % Crmplant      ( arraysize, mplant          ) )
ALLOCATE ( casaflux_data % kplant        ( arraysize, mplant          ) )
ALLOCATE ( casaflux_data % klitter       ( arraysize, mlitter         ) ) 
ALLOCATE ( casaflux_data % ksoil         ( arraysize, msoil           ) )
ALLOCATE ( casaflux_data % fromLtoCO2    ( arraysize, mlitter         ) ) 
ALLOCATE ( casaflux_data % fromStoCO2    ( arraysize, msoil           ) ) 
ALLOCATE ( casaflux_data % FluxCtolitter ( arraysize, mlitter         ) ) 
ALLOCATE ( casaflux_data % FluxNtolitter ( arraysize, mlitter         ) ) 
ALLOCATE ( casaflux_data % FluxPtolitter ( arraysize, mlitter         ) ) 
ALLOCATE ( casaflux_data % FluxCtosoil   ( arraysize, msoil           ) )       
ALLOCATE ( casaflux_data % FluxNtosoil   ( arraysize, msoil           ) )       
ALLOCATE ( casaflux_data % FluxPtosoil   ( arraysize, msoil           ) )       
ALLOCATE ( casaflux_data % fromPtoL      ( arraysize, mlitter, mplant ) )         
ALLOCATE ( casaflux_data % fromLtoS      ( arraysize, msoil, mlitter  ) )         
ALLOCATE ( casaflux_data % fromStoS      ( arraysize, msoil, msoil    ) )         
  
ALLOCATE ( casaflux_data % Cplant_turnover                    ( arraysize, mplant ) )
ALLOCATE ( casaflux_data % Cplant_turnover_disturbance        ( arraysize ) )     
ALLOCATE ( casaflux_data % Cplant_turnover_crowding           ( arraysize ) )    
ALLOCATE ( casaflux_data % Cplant_turnover_resource_limitation( arraysize ) )   

END SUBROUTINE alloc_casa_flux_data_type

SUBROUTINE zero_casa_flux_data_type( casaflux_data ) 

IMPLICIT NONE

TYPE (casa_flux_data_type), INTENT(INOUT) :: casaflux_data 

casaflux_data % Cgpp          (:)     = 0.0 
casaflux_data % Cnpp          (:)     = 0.0
casaflux_data % Cnbp          (:)     = 0.0
casaflux_data % Crp           (:)     = 0.0
casaflux_data % Crgplant      (:)     = 0.0
casaflux_data % Nminfix       (:)     = 0.0
casaflux_data % Nminuptake    (:)     = 0.0
casaflux_data % Plabuptake    (:)     = 0.0
casaflux_data % Clabloss      (:)     = 0.0 
casaflux_data % fracClabile   (:)     = 0.0
casaflux_data % stemnpp       (:)     = 0.0 
casaflux_data % frac_sapwood  (:)     = 0.0 
casaflux_data % Cplant_turnover_tot  (:)     = 0.0 
casaflux_data % sapwood_area  (:)     = 0.0   
casaflux_data % Crsoil        (:)     = 0.0
casaflux_data % Nmindep       (:)     = 0.0
casaflux_data % Nminloss      (:)     = 0.0
casaflux_data % Nminleach     (:)     = 0.0 
casaflux_data % Nupland       (:)     = 0.0
casaflux_data % Nlittermin    (:)     = 0.0
casaflux_data % Nsmin         (:)     = 0.0
casaflux_data % Nsimm         (:)     = 0.0
casaflux_data % Nsnet         (:)     = 0.0
casaflux_data % fNminloss     (:)     = 0.0
casaflux_data % fNminleach    (:)     = 0.0 
casaflux_data % Pdep          (:)     = 0.0
casaflux_data % Pwea          (:)     = 0.0 
casaflux_data % Pleach        (:)     = 0.0 
casaflux_data % Ploss         (:)     = 0.0 
casaflux_data % Pupland       (:)     = 0.0
casaflux_data % Plittermin    (:)     = 0.0
casaflux_data % Psmin         (:)     = 0.0
casaflux_data % Psimm         (:)     = 0.0 
casaflux_data % Psnet         (:)     = 0.0
casaflux_data % fPleach       (:)     = 0.0
casaflux_data % kplab         (:)     = 0.0
casaflux_data % kpsorb        (:)     = 0.0
casaflux_data % kpocc         (:)     = 0.0
casaflux_data % kmlabp        (:)     = 0.0
casaflux_data % Psorbmax      (:)     = 0.0 
casaflux_data % Cnep          (:)     = 0.0    
casaflux_data % FluxCtoCO2    (:)     = 0.0   
casaflux_data % FluxCtohwp    (:)     = 0.0   
casaflux_data % FluxNtohwp    (:)     = 0.0   
casaflux_data % FluxPtohwp    (:)     = 0.0   
casaflux_data % FluxCtoclear  (:)     = 0.0     
casaflux_data % FluxNtoclear  (:)     = 0.0     
casaflux_data % FluxPtoclear  (:)     = 0.0    
casaflux_data % CtransferLUC  (:)     = 0.0    
casaflux_data % meangpp       (:)     = 0.0     
casaflux_data % meanrleaf     (:)     = 0.0             
casaflux_data % fracCalloc    (:,:)   = 0.0     
casaflux_data % fracNalloc    (:,:)   = 0.0     
casaflux_data % fracPalloc    (:,:)   = 0.0     
casaflux_data % Crmplant      (:,:)   = 0.0  
casaflux_data % kplant        (:,:)   = 0.0
casaflux_data % klitter       (:,:)   = 0.0
casaflux_data % ksoil         (:,:)   = 0.0
casaflux_data % fromLtoCO2    (:,:)   = 0.0 
casaflux_data % fromStoCO2    (:,:)   = 0.0
casaflux_data % FluxCtolitter (:,:)   = 0.0 
casaflux_data % FluxNtolitter (:,:)   = 0.0 
casaflux_data % FluxPtolitter (:,:)   = 0.0 
casaflux_data % FluxCtosoil   (:,:)   = 0.0
casaflux_data % FluxNtosoil   (:,:)   = 0.0
casaflux_data % FluxPtosoil   (:,:)   = 0.0
casaflux_data % fromPtoL      (:,:,:) = 0.0
casaflux_data % fromLtoS      (:,:,:) = 0.0
casaflux_data % fromStoS      (:,:,:) = 0.0        
casaflux_data % Cplant_turnover                     (:,:) = 0.0
casaflux_data % Cplant_turnover_disturbance         (:)   = 0.0
casaflux_data % Cplant_turnover_crowding            (:)   = 0.0
casaflux_data % Cplant_turnover_resource_limitation (:)   = 0.0

RETURN
END SUBROUTINE zero_casa_flux_data_type

SUBROUTINE assoc_casa_flux_type( casaflux, casaflux_data ) 

IMPLICIT NONE

TYPE (casa_flux_type),      INTENT(INOUT) :: casaflux
TYPE (casa_flux_data_type), INTENT(INOUT), TARGET :: casaflux_data

casaflux % Cgpp          => casaflux_data % Cgpp            
casaflux % Cnpp          => casaflux_data % Cnpp           
casaflux % Cnpp          => casaflux_data % Cnbp           
casaflux % Crp           => casaflux_data % Crp            
casaflux % Crgplant      => casaflux_data % Crgplant       
casaflux % Nminfix       => casaflux_data % Nminfix        
casaflux % Nminuptake    => casaflux_data % Nminuptake     
casaflux % Plabuptake    => casaflux_data % Plabuptake     
casaflux % Clabloss      => casaflux_data % Clabloss       
casaflux % fracClabile   => casaflux_data % fracClabile    
casaflux % stemnpp       => casaflux_data % stemnpp        
casaflux % frac_sapwood  => casaflux_data % frac_sapwood   
casaflux % Cplant_turnover_tot  => casaflux_data % Cplant_turnover_tot
casaflux % sapwood_area  => casaflux_data % sapwood_area   
casaflux % Crsoil        => casaflux_data % Crsoil         
casaflux % Nmindep       => casaflux_data % Nmindep        
casaflux % Nminloss      => casaflux_data % Nminloss       
casaflux % Nminleach     => casaflux_data % Nminleach      
casaflux % Nupland       => casaflux_data % Nupland        
casaflux % Nlittermin    => casaflux_data % Nlittermin     
casaflux % Nsmin         => casaflux_data % Nsmin          
casaflux % Nsimm         => casaflux_data % Nsimm          
casaflux % Nsnet         => casaflux_data % Nsnet          
casaflux % fNminloss     => casaflux_data % fNminloss      
casaflux % fNminleach    => casaflux_data % fNminleach     
casaflux % Pdep          => casaflux_data % Pdep           
casaflux % Pwea          => casaflux_data % Pwea           
casaflux % Pleach        => casaflux_data % Pleach         
casaflux % Ploss         => casaflux_data % Ploss          
casaflux % Pupland       => casaflux_data % Pupland        
casaflux % Plittermin    => casaflux_data % Plittermin     
casaflux % Psmin         => casaflux_data % Psmin          
casaflux % Psimm         => casaflux_data % Psimm          
casaflux % Psnet         => casaflux_data % Psnet          
casaflux % fPleach       => casaflux_data % fPleach        
casaflux % kplab         => casaflux_data % kplab          
casaflux % kpsorb        => casaflux_data % kpsorb         
casaflux % kpocc         => casaflux_data % kpocc          
casaflux % kmlabp        => casaflux_data % kmlabp         
casaflux % Psorbmax      => casaflux_data % Psorbmax       
casaflux % Cnep          => casaflux_data % Cnep           
casaflux % FluxCtoCO2    => casaflux_data % FluxCtoCO2     
casaflux % FluxCtohwp    => casaflux_data % FluxCtohwp     
casaflux % FluxNtohwp    => casaflux_data % FluxNtohwp     
casaflux % FluxPtohwp    => casaflux_data % FluxPtohwp     
casaflux % FluxCtoclear  => casaflux_data % FluxCtoclear   
casaflux % FluxNtoclear  => casaflux_data % FluxNtoclear   
casaflux % FluxPtoclear  => casaflux_data % FluxPtoclear   
casaflux % CtransferLUC  => casaflux_data % CtransferLUC   
casaflux % meangpp       => casaflux_data % meangpp        
casaflux % meanrleaf     => casaflux_data % meanrleaf      

casaflux % fracCalloc    => casaflux_data % fracCalloc     
casaflux % fracNalloc    => casaflux_data % fracNalloc     
casaflux % fracPalloc    => casaflux_data % fracPalloc     
casaflux % Crmplant      => casaflux_data % Crmplant       
casaflux % kplant        => casaflux_data % kplant         
casaflux % klitter       => casaflux_data % klitter         
casaflux % ksoil         => casaflux_data % ksoil          
casaflux % fromLtoCO2    => casaflux_data % fromLtoCO2      
casaflux % fromStoCO2    => casaflux_data % fromStoCO2      
casaflux % FluxCtolitter => casaflux_data % FluxCtolitter   
casaflux % FluxNtolitter => casaflux_data % FluxNtolitter   
casaflux % FluxPtolitter => casaflux_data % FluxPtolitter   
casaflux % FluxCtosoil   => casaflux_data % FluxCtosoil           
casaflux % FluxNtosoil   => casaflux_data % FluxNtosoil           
casaflux % FluxPtosoil   => casaflux_data % FluxPtosoil           
casaflux % fromPtoL      => casaflux_data % fromPtoL           
casaflux % fromLtoS      => casaflux_data % fromLtoS           
casaflux % fromStoS      => casaflux_data % fromStoS           

casaflux % Cplant_turnover     => casaflux_data % Cplant_turnover                           
casaflux % Cplant_turnover_disturbance                                         &
                               => casaflux_data % Cplant_turnover_disturbance               
casaflux % Cplant_turnover_crowding                                            &
                               => casaflux_data % Cplant_turnover_crowding                
casaflux % Cplant_turnover_resource_limitation                                 &
                         => casaflux_data % Cplant_turnover_resource_limitation

END SUBROUTINE assoc_casa_flux_type

END MODULE casa_flux_type_mod
