MODULE sum_casa_flux_type_mod

IMPLICIT NONE

PUBLIC :: alloc_sum_casa_flux_data_type
PUBLIC :: assoc_sum_casa_flux_type
PUBLIC :: zero_sum_casa_flux_data
PUBLIC :: update_sum_casa_flux

CONTAINS

SUBROUTINE alloc_sum_casa_flux_data_type( sum_casaflux_data, arraysize)
       
USE casadimension, ONLY: mplant, mlitter, msoil
USE casa_flux_type_mod, ONLY: casa_flux_data_type

IMPLICIT NONE

TYPE (casa_flux_data_type), INTENT(INOUT) :: sum_casaflux_data
INTEGER,                    INTENT(IN)    :: arraysize

ALLOCATE ( sum_casaflux_data % Cgpp         (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Cnpp         (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Crp          (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Crgplant     (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nminfix      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nminuptake   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Plabuptake   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Clabloss     (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % fracClabile  (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % stemnpp      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % frac_sapwood (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % sapwood_area (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Crsoil       (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nmindep      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nminloss     (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nminleach    (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nupland      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nlittermin   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nsmin        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nsimm        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Nsnet        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % fNminloss    (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % fNminleach   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Pdep         (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Pwea         (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Pleach       (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Ploss        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Pupland      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Plittermin   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Psmin        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Psimm        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Psnet        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % fPleach      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % kplab        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % kpsorb       (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % kpocc        (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % kmlabp       (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Psorbmax     (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % Cnep         (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxCtoCO2   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxCtohwp   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxNtohwp   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxPtohwp   (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxCtoclear (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxNtoclear (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % FluxPtoclear (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % CtransferLUC (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % meangpp      (  arraysize ) ) 
ALLOCATE ( sum_casaflux_data % meanrleaf    (  arraysize ) ) 
  
ALLOCATE ( sum_casaflux_data % fracCalloc    ( arraysize, mplant          ) )
ALLOCATE ( sum_casaflux_data % fracNalloc    ( arraysize, mplant          ) )
ALLOCATE ( sum_casaflux_data % fracPalloc    ( arraysize, mplant          ) )
ALLOCATE ( sum_casaflux_data % Crmplant      ( arraysize, mplant          ) )
ALLOCATE ( sum_casaflux_data % kplant        ( arraysize, mplant          ) )
ALLOCATE ( sum_casaflux_data % klitter       ( arraysize, mlitter         ) ) 
ALLOCATE ( sum_casaflux_data % ksoil         ( arraysize, msoil           ) )
ALLOCATE ( sum_casaflux_data % fromLtoCO2    ( arraysize, mlitter         ) ) 
ALLOCATE ( sum_casaflux_data % fromStoCO2    ( arraysize, msoil           ) ) 
ALLOCATE ( sum_casaflux_data % FluxCtolitter ( arraysize, mlitter         ) ) 
ALLOCATE ( sum_casaflux_data % FluxNtolitter ( arraysize, mlitter         ) ) 
ALLOCATE ( sum_casaflux_data % FluxPtolitter ( arraysize, mlitter         ) ) 
ALLOCATE ( sum_casaflux_data % FluxCtosoil   ( arraysize, msoil           ) )       
ALLOCATE ( sum_casaflux_data % FluxNtosoil   ( arraysize, msoil           ) )       
ALLOCATE ( sum_casaflux_data % FluxPtosoil   ( arraysize, msoil           ) )       
ALLOCATE ( sum_casaflux_data % fromPtoL      ( arraysize, mlitter, mplant ) )         
ALLOCATE ( sum_casaflux_data % fromLtoS      ( arraysize, msoil, mlitter  ) )         
ALLOCATE ( sum_casaflux_data % fromStoS      ( arraysize, msoil, msoil    ) )         
  
ALLOCATE ( sum_casaflux_data % Cplant_turnover                    ( arraysize, mplant ) )
ALLOCATE ( sum_casaflux_data % Cplant_turnover_disturbance        ( arraysize ) )     
ALLOCATE ( sum_casaflux_data % Cplant_turnover_crowding           ( arraysize ) )    
ALLOCATE ( sum_casaflux_data % Cplant_turnover_resource_limitation( arraysize ) )   

RETURN
END SUBROUTINE alloc_sum_casa_flux_data_type

SUBROUTINE assoc_sum_casa_flux_type( sum_casaflux, sum_casaflux_data ) 

USE casa_flux_type_mod, ONLY: casa_flux_data_type
USE casa_flux_type_mod, ONLY: casa_flux_type

 IMPLICIT NONE

TYPE (casa_flux_type),      INTENT(INOUT) :: sum_casaflux
TYPE (casa_flux_data_type), INTENT(INOUT), TARGET :: sum_casaflux_data

sum_casaflux % Cgpp          => sum_casaflux_data % Cgpp           
sum_casaflux % Cnpp          => sum_casaflux_data % Cnpp           
sum_casaflux % Crp           => sum_casaflux_data % Crp            
sum_casaflux % Crgplant      => sum_casaflux_data % Crgplant       
sum_casaflux % Nminfix       => sum_casaflux_data % Nminfix        
sum_casaflux % Nminuptake    => sum_casaflux_data % Nminuptake     
sum_casaflux % Plabuptake    => sum_casaflux_data % Plabuptake     
sum_casaflux % Clabloss      => sum_casaflux_data % Clabloss       
sum_casaflux % fracClabile   => sum_casaflux_data % fracClabile    
sum_casaflux % stemnpp       => sum_casaflux_data % stemnpp        
sum_casaflux % frac_sapwood  => sum_casaflux_data % frac_sapwood   
sum_casaflux % sapwood_area  => sum_casaflux_data % sapwood_area   
sum_casaflux % Crsoil        => sum_casaflux_data % Crsoil         
sum_casaflux % Nmindep       => sum_casaflux_data % Nmindep        
sum_casaflux % Nminloss      => sum_casaflux_data % Nminloss       
sum_casaflux % Nminleach     => sum_casaflux_data % Nminleach      
sum_casaflux % Nupland       => sum_casaflux_data % Nupland        
sum_casaflux % Nlittermin    => sum_casaflux_data % Nlittermin     
sum_casaflux % Nsmin         => sum_casaflux_data % Nsmin          
sum_casaflux % Nsimm         => sum_casaflux_data % Nsimm          
sum_casaflux % Nsnet         => sum_casaflux_data % Nsnet          
sum_casaflux % fNminloss     => sum_casaflux_data % fNminloss      
sum_casaflux % fNminleach    => sum_casaflux_data % fNminleach     
sum_casaflux % Pdep          => sum_casaflux_data % Pdep           
sum_casaflux % Pwea          => sum_casaflux_data % Pwea           
sum_casaflux % Pleach        => sum_casaflux_data % Pleach         
sum_casaflux % Ploss         => sum_casaflux_data % Ploss          
sum_casaflux % Pupland       => sum_casaflux_data % Pupland        
sum_casaflux % Plittermin    => sum_casaflux_data % Plittermin     
sum_casaflux % Psmin         => sum_casaflux_data % Psmin          
sum_casaflux % Psimm         => sum_casaflux_data % Psimm          
sum_casaflux % Psnet         => sum_casaflux_data % Psnet          
sum_casaflux % fPleach       => sum_casaflux_data % fPleach        
sum_casaflux % kplab         => sum_casaflux_data % kplab          
sum_casaflux % kpsorb        => sum_casaflux_data % kpsorb         
sum_casaflux % kpocc         => sum_casaflux_data % kpocc          
sum_casaflux % kmlabp        => sum_casaflux_data % kmlabp         
sum_casaflux % Psorbmax      => sum_casaflux_data % Psorbmax       
sum_casaflux % Cnep          => sum_casaflux_data % Cnep           
sum_casaflux % FluxCtoCO2    => sum_casaflux_data % FluxCtoCO2     
sum_casaflux % FluxCtohwp    => sum_casaflux_data % FluxCtohwp     
sum_casaflux % FluxNtohwp    => sum_casaflux_data % FluxNtohwp     
sum_casaflux % FluxPtohwp    => sum_casaflux_data % FluxPtohwp     
sum_casaflux % FluxCtoclear  => sum_casaflux_data % FluxCtoclear   
sum_casaflux % FluxNtoclear  => sum_casaflux_data % FluxNtoclear   
sum_casaflux % FluxPtoclear  => sum_casaflux_data % FluxPtoclear   
sum_casaflux % CtransferLUC  => sum_casaflux_data % CtransferLUC   
sum_casaflux % meangpp       => sum_casaflux_data % meangpp        
sum_casaflux % meanrleaf     => sum_casaflux_data % meanrleaf      
           
sum_casaflux % fracCalloc    => sum_casaflux_data % fracCalloc     
sum_casaflux % fracNalloc    => sum_casaflux_data % fracNalloc     
sum_casaflux % fracPalloc    => sum_casaflux_data % fracPalloc     
sum_casaflux % Crmplant      => sum_casaflux_data % Crmplant       
sum_casaflux % kplant        => sum_casaflux_data % kplant         
sum_casaflux % klitter       => sum_casaflux_data % klitter        
sum_casaflux % ksoil         => sum_casaflux_data % ksoil          
sum_casaflux % fromLtoCO2    => sum_casaflux_data % fromLtoCO2     
sum_casaflux % fromStoCO2    => sum_casaflux_data % fromStoCO2     
sum_casaflux % FluxCtolitter => sum_casaflux_data % FluxCtolitter 
sum_casaflux % FluxNtolitter => sum_casaflux_data % FluxNtolitter 
sum_casaflux % FluxPtolitter => sum_casaflux_data % FluxPtolitter 
sum_casaflux % FluxCtosoil   => sum_casaflux_data % FluxCtosoil    
sum_casaflux % FluxNtosoil   => sum_casaflux_data % FluxNtosoil    
sum_casaflux % FluxPtosoil   => sum_casaflux_data % FluxPtosoil    
sum_casaflux % fromPtoL      => sum_casaflux_data % fromPtoL       
sum_casaflux % fromLtoS      => sum_casaflux_data % fromLtoS       
sum_casaflux % fromStoS      => sum_casaflux_data % fromStoS       

sum_casaflux % Cplant_turnover                => sum_casaflux % Cplant_turnover                     
sum_casaflux % Cplant_turnover_disturbance    =>                               &
                                 sum_casaflux_data % Cplant_turnover_disturbance         
sum_casaflux % Cplant_turnover_crowding       =>                               &
                                    sum_casaflux_data % Cplant_turnover_crowding            
sum_casaflux % Cplant_turnover_resource_limitation =>                          &
                         sum_casaflux_data % Cplant_turnover_resource_limitation 

RETURN
END SUBROUTINE assoc_sum_casa_flux_type

SUBROUTINE zero_sum_casa_flux_data( sum_casaflux )

USE casa_flux_type_mod, ONLY: casa_flux_data_type

IMPLICIT NONE

TYPE (casa_flux_data_type), INTENT(INOUT) :: sum_casaflux

sum_casaflux % Cgpp          = 0.0 
sum_casaflux % Cnpp          = 0.0 
sum_casaflux % Crp           = 0.0 
sum_casaflux % Crgplant      = 0.0 
sum_casaflux % Nminfix       = 0.0 
sum_casaflux % Nminuptake    = 0.0 
sum_casaflux % Plabuptake    = 0.0 
sum_casaflux % Clabloss      = 0.0 
sum_casaflux % fracClabile   = 0.0 
sum_casaflux % stemnpp       = 0.0 
sum_casaflux % frac_sapwood  = 0.0 
sum_casaflux % sapwood_area  = 0.0 
sum_casaflux % Crsoil        = 0.0 
sum_casaflux % Nmindep       = 0.0 
sum_casaflux % Nminloss      = 0.0 
sum_casaflux % Nminleach     = 0.0 
sum_casaflux % Nupland       = 0.0 
sum_casaflux % Nlittermin    = 0.0 
sum_casaflux % Nsmin         = 0.0 
sum_casaflux % Nsimm         = 0.0 
sum_casaflux % Nsnet         = 0.0 
sum_casaflux % fNminloss     = 0.0 
sum_casaflux % fNminleach    = 0.0 
sum_casaflux % Pdep          = 0.0 
sum_casaflux % Pwea          = 0.0 
sum_casaflux % Pleach        = 0.0 
sum_casaflux % Ploss         = 0.0 
sum_casaflux % Pupland       = 0.0 
sum_casaflux % Plittermin    = 0.0 
sum_casaflux % Psmin         = 0.0 
sum_casaflux % Psimm         = 0.0 
sum_casaflux % Psnet         = 0.0 
sum_casaflux % fPleach       = 0.0 
sum_casaflux % kplab         = 0.0 
sum_casaflux % kpsorb        = 0.0 
sum_casaflux % kpocc         = 0.0 
sum_casaflux % kmlabp        = 0.0 
sum_casaflux % Psorbmax      = 0.0 
sum_casaflux % Cnep          = 0.0 
sum_casaflux % FluxCtoCO2    = 0.0 
sum_casaflux % FluxCtohwp    = 0.0 
sum_casaflux % FluxNtohwp    = 0.0 
sum_casaflux % FluxPtohwp    = 0.0 
sum_casaflux % FluxCtoclear  = 0.0 
sum_casaflux % FluxNtoclear  = 0.0 
sum_casaflux % FluxPtoclear  = 0.0 
sum_casaflux % CtransferLUC  = 0.0 
sum_casaflux % meangpp       = 0.0 
sum_casaflux % meanrleaf     = 0.0 
            
sum_casaflux % fracCalloc    = 0.0 
sum_casaflux % fracNalloc    = 0.0 
sum_casaflux % fracPalloc    = 0.0 
sum_casaflux % Crmplant      = 0.0 
sum_casaflux % kplant        = 0.0 
sum_casaflux % klitter       = 0.0 
sum_casaflux % ksoil         = 0.0 
sum_casaflux % fromLtoCO2    = 0.0 
sum_casaflux % fromStoCO2    = 0.0 
sum_casaflux % FluxCtolitter = 0.0
sum_casaflux % FluxNtolitter = 0.0
sum_casaflux % FluxPtolitter = 0.0
sum_casaflux % FluxCtosoil   = 0.0 
sum_casaflux % FluxNtosoil   = 0.0 
sum_casaflux % FluxPtosoil   = 0.0 
sum_casaflux % fromPtoL      = 0.0 
sum_casaflux % fromLtoS      = 0.0 
sum_casaflux % fromStoS      = 0.0 

sum_casaflux % Cplant_turnover                     = 0.0 
sum_casaflux % Cplant_turnover_disturbance         = 0.0 
sum_casaflux % Cplant_turnover_crowding            = 0.0 
sum_casaflux % Cplant_turnover_resource_limitation = 0.0 

END SUBROUTINE zero_sum_casa_flux_data

!jhan:pass individual elements of casapool
!do the calcshere with the pointers to the _data fields
SUBROUTINE update_sum_casa_flux( sum_casapool, sum_casaflux,         &
                              casapool, casaflux, sum_now, average_now, nsteps )

USE casa_flux_type_mod, ONLY: casa_flux_data_type
USE casa_pool_type_mod, ONLY: casa_pool_data_type

IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: sum_casapool
TYPE (casa_flux_data_type), INTENT(INOUT) :: sum_casaflux
TYPE (casa_pool_data_type), INTENT(IN)    :: casapool
TYPE (casa_flux_data_type), INTENT(IN)    :: casaflux
LOGICAL,               INTENT(IN)    :: sum_now, average_now
INTEGER,               INTENT(IN)    :: nsteps

!local vars
INTEGER             :: rsteps

rsteps =  REAL(nsteps)

IF (sum_now) THEN

  sum_casaflux%Cgpp        = sum_casaflux%Cgpp        + casaflux%Cgpp
  sum_casaflux%Cnpp        = sum_casaflux%Cnpp        + casaflux%Cnpp
  sum_casaflux%Crp         = sum_casaflux%Crp         + casaflux%Crp
  sum_casaflux%Crgplant    = sum_casaflux%Crgplant    + casaflux%Crgplant
  sum_casaflux%Nminfix     = sum_casaflux%Nminfix     + casaflux%Nminfix
  sum_casaflux%Nminuptake  = sum_casaflux%Nminuptake  + casaflux%Nminuptake
  sum_casaflux%Plabuptake  = sum_casaflux%Plabuptake  + casaflux%Plabuptake
  sum_casaflux%Clabloss    = sum_casaflux%Clabloss    + casaflux%Clabloss
  sum_casaflux%fracClabile = sum_casaflux%fracClabile + casaflux%fracClabile
  
  sum_casaflux%fracCalloc  = sum_casaflux%fracCalloc  +                        &
                             ( casaflux%fracCalloc * casapool%cplant )
  sum_casaflux%fracNalloc  = sum_casaflux%fracNalloc + casaflux%fracNalloc
  sum_casaflux%fracPalloc  = sum_casaflux%fracPalloc + casaflux%fracPalloc
  
  sum_casaflux%kplant      = sum_casaflux%kplant     +                         &
                             ( casaflux%kplant * casapool%cplant )
  
  sum_casaflux%Crmplant    = sum_casaflux%Crmplant   + casaflux%Crmplant
  sum_casaflux%fromPtoL    = sum_casaflux%fromPtoL   + casaflux%fromPtoL
  sum_casaflux%Cnep        = sum_casaflux%Cnep       + casaflux%Cnep
  sum_casaflux%Crsoil      = sum_casaflux%Crsoil     + casaflux%Crsoil
  sum_casaflux%Nmindep     = sum_casaflux%Nmindep    + casaflux%Nmindep
  sum_casaflux%Nminloss    = sum_casaflux%Nminloss   + casaflux%Nminloss
  sum_casaflux%Nminleach   = sum_casaflux%Nminleach  + casaflux%Nminleach
  sum_casaflux%Nupland     = sum_casaflux%Nupland    + casaflux%Nupland
  sum_casaflux%Nlittermin  = sum_casaflux%Nlittermin + casaflux%Nlittermin
  sum_casaflux%Nsmin       = sum_casaflux%Nsmin      + casaflux%Nsmin
  sum_casaflux%Nsimm       = sum_casaflux%Nsimm      + casaflux%Nsimm
  sum_casaflux%Nsnet       = sum_casaflux%Nsnet      + casaflux%Nsnet
  sum_casaflux%fNminloss   = sum_casaflux%fNminloss  + casaflux%fNminloss
  sum_casaflux%fNminleach  = sum_casaflux%fNminleach + casaflux%fNminleach
  sum_casaflux%Pdep        = sum_casaflux%Pdep       + casaflux%Pdep
  sum_casaflux%Pwea        = sum_casaflux%Pwea       + casaflux%Pwea
  sum_casaflux%Pleach      = sum_casaflux%Pleach     + casaflux%Pleach
  sum_casaflux%Ploss       = sum_casaflux%Ploss      + casaflux%Ploss
  sum_casaflux%Pupland     = sum_casaflux%Pupland    + casaflux%Pupland
  sum_casaflux%Plittermin  = sum_casaflux%Plittermin + casaflux%Plittermin
  sum_casaflux%Psmin       = sum_casaflux%Psmin      + casaflux%Psmin
  sum_casaflux%Psimm       = sum_casaflux%Psimm      + casaflux%Psimm
  sum_casaflux%Psnet       = sum_casaflux%Psnet      + casaflux%Psnet
  sum_casaflux%fPleach     = sum_casaflux%fPleach    + casaflux%fPleach
  sum_casaflux%kplab       = sum_casaflux%kplab      + casaflux%kplab
  sum_casaflux%kpsorb      = sum_casaflux%kpsorb      + casaflux%kpsorb
  sum_casaflux%kpocc       = sum_casaflux%kpocc      + casaflux%kpocc
  sum_casaflux%kmlabP      = sum_casaflux%kmlabP     + casaflux%kmlabP
  sum_casaflux%Psorbmax    = sum_casaflux%Psorbmax   + casaflux%Psorbmax
  sum_casaflux%klitter     = sum_casaflux%klitter    + casaflux%klitter
  sum_casaflux%ksoil       = sum_casaflux%ksoil      + casaflux%ksoil
  sum_casaflux%fromLtoS    = sum_casaflux%fromLtoS   + casaflux%fromLtoS
  sum_casaflux%fromStoS    = sum_casaflux%fromStoS   + casaflux%fromStoS
  sum_casaflux%fromLtoCO2  = sum_casaflux%fromLtoCO2 + casaflux%fromLtoCO2
  sum_casaflux%fromStoCO2  = sum_casaflux%fromStoCO2 + casaflux%fromStoCO2
  
  sum_casaflux%stemnpp      = sum_casaflux%stemnpp      + casaflux%stemnpp
  sum_casaflux%frac_sapwood = sum_casaflux%frac_sapwood + casaflux%frac_sapwood
  sum_casaflux%sapwood_area = sum_casaflux%sapwood_area + casaflux%sapwood_area
  sum_casaflux%FluxCtoco2   = sum_casaflux%FluxCtoco2   + casaflux%FluxCtoco2
  sum_casaflux%FluxCtosoil  = sum_casaflux%FluxCtosoil  + casaflux%FluxCtosoil
  sum_casaflux%FluxNtosoil  = sum_casaflux%FluxNtosoil  + casaflux%FluxNtosoil
  sum_casaflux%FluxPtosoil  = sum_casaflux%FluxPtosoil  + casaflux%FluxPtosoil
  
  sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter  +                   &
                               casaflux%FluxCtolitter
  
  sum_casaflux%FluxNtolitter = sum_casaflux%FluxNtolitter +                    &
                               casaflux%FluxNtolitter
  
  sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter  +                   &
                               casaflux%FluxPtolitter

  sum_casaflux%Cplant_turnover = sum_casaflux%Cplant_turnover +                &
                                 casaflux%Cplant_turnover

  sum_casaflux%Cplant_turnover_disturbance =                                   &
                                    sum_casaflux%Cplant_turnover_disturbance + &
                                    casaflux%Cplant_turnover_disturbance

  sum_casaflux%Cplant_turnover_crowding =                                      &
                                       sum_casaflux%Cplant_turnover_crowding + &
                                       casaflux%Cplant_turnover_crowding

  sum_casaflux%Cplant_turnover_resource_limitation =                           &  
                           sum_casaflux%Cplant_turnover_resource_limitation +  &
                           casaflux%Cplant_turnover_resource_limitation
ENDIF

IF (average_now) THEN

  WHERE ( sum_casapool%Cplant .GT. 1.e-12 )
    sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc / sum_casapool%Cplant
    sum_casaflux%kplant     =  sum_casaflux%kplant / sum_casapool%Cplant
  ELSEWHERE
    sum_casaflux%fracCalloc = 0.0
    sum_casaflux%kplant     = 0.0
  ENDWHERE

  sum_casaflux%Cgpp          = sum_casaflux%Cgpp          / rsteps
  sum_casaflux%Cnpp          = sum_casaflux%Cnpp          / rsteps
  sum_casaflux%Crp           = sum_casaflux%Crp           / rsteps
  sum_casaflux%Crgplant      = sum_casaflux%Crgplant      / rsteps
  sum_casaflux%Nminfix       = sum_casaflux%Nminfix       / rsteps
  sum_casaflux%Nminuptake    = sum_casaflux%Nminuptake    / rsteps
  sum_casaflux%Plabuptake    = sum_casaflux%Plabuptake    / rsteps
  sum_casaflux%Clabloss      = sum_casaflux%Clabloss      / rsteps
  sum_casaflux%fracClabile   = sum_casaflux%fracClabile   / rsteps
  sum_casaflux%fracNalloc    = sum_casaflux%fracNalloc    / rsteps
  sum_casaflux%fracPalloc    = sum_casaflux%fracPalloc    / rsteps
  sum_casaflux%Crmplant      = sum_casaflux%Crmplant      / rsteps
  sum_casaflux%fromPtoL      = sum_casaflux%fromPtoL      / rsteps
  sum_casaflux%Cnep          = sum_casaflux%Cnep          / rsteps
  sum_casaflux%Crsoil        = sum_casaflux%Crsoil        / rsteps
  sum_casaflux%Nmindep       = sum_casaflux%Nmindep       / rsteps
  sum_casaflux%Nminloss      = sum_casaflux%Nminloss      / rsteps
  sum_casaflux%Nminleach     = sum_casaflux%Nminleach     / rsteps
  sum_casaflux%Nupland       = sum_casaflux%Nupland       / rsteps
  sum_casaflux%Nlittermin    = sum_casaflux%Nlittermin    / rsteps
  sum_casaflux%Nsmin         = sum_casaflux%Nsmin         / rsteps
  sum_casaflux%Nsimm         = sum_casaflux%Nsimm         / rsteps
  sum_casaflux%Nsnet         = sum_casaflux%Nsnet         / rsteps
  sum_casaflux%fNminloss     = sum_casaflux%fNminloss     / rsteps
  sum_casaflux%fNminleach    = sum_casaflux%fNminleach    / rsteps
  sum_casaflux%Pdep          = sum_casaflux%Pdep          / rsteps
  sum_casaflux%Pwea          = sum_casaflux%Pwea          / rsteps
  sum_casaflux%Pleach        = sum_casaflux%Pleach        / rsteps
  sum_casaflux%Ploss         = sum_casaflux%Ploss         / rsteps
  sum_casaflux%Pupland       = sum_casaflux%Pupland       / rsteps
  sum_casaflux%Plittermin    = sum_casaflux%Plittermin    / rsteps
  sum_casaflux%Psmin         = sum_casaflux%Psmin         / rsteps
  sum_casaflux%Psimm         = sum_casaflux%Psimm         / rsteps
  sum_casaflux%Psnet         = sum_casaflux%Psnet         / rsteps
  sum_casaflux%fPleach       = sum_casaflux%fPleach       / rsteps
  sum_casaflux%kplab         = sum_casaflux%kplab         / rsteps
  sum_casaflux%kpsorb        = sum_casaflux%kpsorb        / rsteps
  sum_casaflux%kpocc         = sum_casaflux%kpocc         / rsteps
  sum_casaflux%kmlabP        = sum_casaflux%kmlabP        / rsteps
  sum_casaflux%Psorbmax      = sum_casaflux%Psorbmax      / rsteps
  sum_casaflux%klitter       = sum_casaflux%klitter       / rsteps
  sum_casaflux%ksoil         = sum_casaflux%ksoil         / rsteps
  sum_casaflux%fromLtoS      = sum_casaflux%fromLtoS      / rsteps
  sum_casaflux%fromStoS      = sum_casaflux%fromStoS      / rsteps
  sum_casaflux%fromLtoCO2    = sum_casaflux%fromLtoCO2    / rsteps
  sum_casaflux%fromStoCO2    = sum_casaflux%fromStoCO2    / rsteps
  sum_casaflux%stemnpp       = sum_casaflux%stemnpp       / rsteps
  sum_casaflux%frac_sapwood  = sum_casaflux%frac_sapwood  / rsteps
  sum_casaflux%sapwood_area  = sum_casaflux%sapwood_area  / rsteps
  sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter / rsteps
  sum_casaflux%FluxNtolitter = sum_casaflux%FluxNtolitter / rsteps
  sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter / rsteps
  sum_casaflux%FluxCtosoil   = sum_casaflux%FluxCtosoil   / rsteps
  sum_casaflux%FluxNtosoil   = sum_casaflux%FluxNtosoil   / rsteps
  sum_casaflux%FluxPtosoil   = sum_casaflux%FluxPtosoil   / rsteps
  sum_casaflux%FluxCtoco2    = sum_casaflux%FluxCtoco2    / rsteps

  sum_casaflux%Cplant_turnover = sum_casaflux%Cplant_turnover / rsteps

  sum_casaflux%Cplant_turnover_disturbance =                                   &
                                   casaflux%Cplant_turnover_disturbance / rsteps

  sum_casaflux% Cplant_turnover_crowding   =                                   &
                                  sum_casaflux%Cplant_turnover_crowding / rsteps

  sum_casaflux%Cplant_turnover_resource_limitation =                           &
                       sum_casaflux%Cplant_turnover_resource_limitation / rsteps

ENDIF

END SUBROUTINE update_sum_casa_flux

END MODULE sum_casa_flux_type_mod
