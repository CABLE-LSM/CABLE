MODULE sum_casa_pool_type_mod

USE cable_other_constants_mod, ONLY: r_2 ! currently DOUBLE precision

IMPLICIT NONE
 
PUBLIC :: alloc_sum_casa_pool_data_type
PUBLIC :: assoc_sum_casa_pool_type
PUBLIC :: zero_sum_casa_pool_data
PUBLIC :: update_sum_casa_pool

CONTAINS

SUBROUTINE alloc_sum_casa_pool_data_type( sum_casapool_data, arraysize ) 

USE casa_pool_type_mod, ONLY: casa_pool_data_type
USE casadimension, ONLY: mplant, mlitter, msoil, mwood

IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: sum_casapool_data 
INTEGER,                    INTENT(IN)    :: arraysize

ALLOCATE ( sum_casapool_data % Clabile        ( arraysize ) )
ALLOCATE ( sum_casapool_data % dClabiledt     ( arraysize ) )
ALLOCATE ( sum_casapool_data % Ctot           ( arraysize ) )
ALLOCATE ( sum_casapool_data % Ctot_0         ( arraysize ) )
ALLOCATE ( sum_casapool_data % Nsoilmin       ( arraysize ) )
ALLOCATE ( sum_casapool_data % Psoillab       ( arraysize ) )
ALLOCATE ( sum_casapool_data % Psoilsorb      ( arraysize ) )
ALLOCATE ( sum_casapool_data % Psoilocc       ( arraysize ) )
ALLOCATE ( sum_casapool_data % dNsoilmindt    ( arraysize ) )
ALLOCATE ( sum_casapool_data % dPsoillabdt    ( arraysize ) )
ALLOCATE ( sum_casapool_data % dPsoilsorbdt   ( arraysize ) )
ALLOCATE ( sum_casapool_data % dPsoiloccdt    ( arraysize ) )

ALLOCATE ( sum_casapool_data % Cplant         ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % Nplant         ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % Pplant         ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % dCplantdt      ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % dNplantdt      ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % dPplantdt      ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % ratioNCplant   ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % ratioNPplant   ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % Clitter        ( arraysize, mlitter ) ) 
ALLOCATE ( sum_casapool_data % Nlitter        ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % Plitter        ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % dClitterdt     ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % dNlitterdt     ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % dPlitterdt     ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % ratioNClitter  ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % ratioNPlitter  ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % Csoil          ( arraysize, msoil ) )   
ALLOCATE ( sum_casapool_data % Nsoil          ( arraysize, msoil ) )    
ALLOCATE ( sum_casapool_data % Psoil          ( arraysize, msoil ) )    
ALLOCATE ( sum_casapool_data % dCsoildt       ( arraysize, msoil ) ) 
ALLOCATE ( sum_casapool_data % dNsoildt       ( arraysize, msoil ) ) 
ALLOCATE ( sum_casapool_data % dPsoildt       ( arraysize, msoil ) ) 
ALLOCATE ( sum_casapool_data % ratioNCsoil    ( arraysize, msoil ) ) 
ALLOCATE ( sum_casapool_data % ratioNCsoilnew ( arraysize, msoil ) ) 
ALLOCATE ( sum_casapool_data % ratioNPsoil    ( arraysize, msoil ) ) 
ALLOCATE ( sum_casapool_data % ratioNCsoilmin ( arraysize, msoil ) )
ALLOCATE ( sum_casapool_data % ratioNCsoilmax ( arraysize, msoil ) )
ALLOCATE ( sum_casapool_data % ratioPcsoil    ( arraysize, msoil ) )
ALLOCATE ( sum_casapool_data % ratioPcplant   ( arraysize, mplant ) ) 
ALLOCATE ( sum_casapool_data % ratioPclitter  ( arraysize, mlitter ) )  
ALLOCATE ( sum_casapool_data % cwoodprod      ( arraysize, mwood ) ) 
ALLOCATE ( sum_casapool_data % nwoodprod      ( arraysize, mwood ) ) 
ALLOCATE ( sum_casapool_data % pwoodprod      ( arraysize, mwood ) ) 

END SUBROUTINE alloc_sum_casa_pool_data_type

SUBROUTINE zero_sum_casa_pool_data( sum_casapool ) 

USE casa_pool_type_mod, ONLY: casa_pool_data_type

IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: sum_casapool

sum_casapool % Clabile        = 0.0 
sum_casapool % dClabiledt     = 0.0 
sum_casapool % Ctot           = 0.0 
sum_casapool % Ctot_0         = 0.0 
sum_casapool % Nsoilmin       = 0.0 
sum_casapool % Psoillab       = 0.0 
sum_casapool % Psoilsorb      = 0.0 
sum_casapool % Psoilocc       = 0.0 
sum_casapool % dNsoilmindt    = 0.0 
sum_casapool % dPsoillabdt    = 0.0 
sum_casapool % dPsoilsorbdt   = 0.0 
sum_casapool % dPsoiloccdt    = 0.0 

sum_casapool % Cplant         = 0.0 
sum_casapool % Nplant         = 0.0 
sum_casapool % Pplant         = 0.0 
sum_casapool % dCplantdt      = 0.0 
sum_casapool % dNplantdt      = 0.0 
sum_casapool % dPplantdt      = 0.0 
sum_casapool % ratioNCplant   = 0.0 
sum_casapool % ratioNPplant   = 0.0 
sum_casapool % Clitter        = 0.0 
sum_casapool % Nlitter        = 0.0 
sum_casapool % Plitter        = 0.0 
sum_casapool % dClitterdt     = 0.0 
sum_casapool % dNlitterdt     = 0.0 
sum_casapool % dPlitterdt     = 0.0 
sum_casapool % ratioNClitter  = 0.0 
sum_casapool % ratioNPlitter  = 0.0 
sum_casapool % Csoil          = 0.0 
sum_casapool % Nsoil          = 0.0 
sum_casapool % Psoil          = 0.0 
sum_casapool % dCsoildt       = 0.0 
sum_casapool % dNsoildt       = 0.0 
sum_casapool % dPsoildt       = 0.0 
sum_casapool % ratioNCsoil    = 0.0 
sum_casapool % ratioNCsoilnew = 0.0 
sum_casapool % ratioNPsoil    = 0.0 
sum_casapool % ratioNCsoilmin = 0.0 
sum_casapool % ratioNCsoilmax = 0.0 
sum_casapool % ratioPcsoil    = 0.0 
sum_casapool % ratioPcplant   = 0.0 
sum_casapool % ratioPclitter  = 0.0 
sum_casapool % cwoodprod      = 0.0 
sum_casapool % nwoodprod      = 0.0 
sum_casapool % pwoodprod      = 0.0 

RETURN
END SUBROUTINE zero_sum_casa_pool_data

SUBROUTINE assoc_sum_casa_pool_type( sum_casapool, sum_casapool_data ) 

USE casa_pool_type_mod, ONLY: casa_pool_data_type
USE casa_pool_type_mod, ONLY: casa_pool_type

IMPLICIT NONE

TYPE (casa_pool_type),      INTENT(INOUT) :: sum_casapool
TYPE (casa_pool_data_type), INTENT(INOUT), TARGET :: sum_casapool_data

sum_casapool % Clabile        => sum_casapool_data % Clabile         
sum_casapool % dClabiledt     => sum_casapool_data % dClabiledt      
sum_casapool % Ctot           => sum_casapool_data % Ctot            
sum_casapool % Ctot_0         => sum_casapool_data % Ctot_0          
sum_casapool % Nsoilmin       => sum_casapool_data % Nsoilmin        
sum_casapool % Psoillab       => sum_casapool_data % Psoillab        
sum_casapool % Psoilsorb      => sum_casapool_data % Psoilsorb       
sum_casapool % Psoilocc       => sum_casapool_data % Psoilocc        
sum_casapool % dNsoilmindt    => sum_casapool_data % dNsoilmindt     
sum_casapool % dPsoillabdt    => sum_casapool_data % dPsoillabdt     
sum_casapool % dPsoilsorbdt   => sum_casapool_data % dPsoilsorbdt    
sum_casapool % dPsoiloccdt    => sum_casapool_data % dPsoiloccdt     

sum_casapool % Cplant         => sum_casapool_data % Cplant          
sum_casapool % Nplant         => sum_casapool_data % Nplant          
sum_casapool % Pplant         => sum_casapool_data % Pplant          
sum_casapool % dCplantdt      => sum_casapool_data % dCplantdt       
sum_casapool % dNplantdt      => sum_casapool_data % dNplantdt       
sum_casapool % dPplantdt      => sum_casapool_data % dPplantdt       
sum_casapool % ratioNCplant   => sum_casapool_data % ratioNCplant    
sum_casapool % ratioNPplant   => sum_casapool_data % ratioNPplant    
sum_casapool % Clitter        => sum_casapool_data % Clitter         
sum_casapool % Nlitter        => sum_casapool_data % Nlitter         
sum_casapool % Plitter        => sum_casapool_data % Plitter         
sum_casapool % dClitterdt     => sum_casapool_data % dClitterdt      
sum_casapool % dNlitterdt     => sum_casapool_data % dNlitterdt      
sum_casapool % dPlitterdt     => sum_casapool_data % dPlitterdt      
sum_casapool % ratioNClitter  => sum_casapool_data % ratioNClitter   
sum_casapool % ratioNPlitter  => sum_casapool_data % ratioNPlitter   
sum_casapool % Csoil          => sum_casapool_data % Csoil           
sum_casapool % Nsoil          => sum_casapool_data % Nsoil           
sum_casapool % Psoil          => sum_casapool_data % Psoil           
sum_casapool % dCsoildt       => sum_casapool_data % dCsoildt        
sum_casapool % dNsoildt       => sum_casapool_data % dNsoildt        
sum_casapool % dPsoildt       => sum_casapool_data % dPsoildt        
sum_casapool % ratioNCsoil    => sum_casapool_data % ratioNCsoil     
sum_casapool % ratioNCsoilnew => sum_casapool_data % ratioNCsoilnew  
sum_casapool % ratioNPsoil    => sum_casapool_data % ratioNPsoil     
sum_casapool % ratioNCsoilmin => sum_casapool_data % ratioNCsoilmin  
sum_casapool % ratioNCsoilmax => sum_casapool_data % ratioNCsoilmax  
sum_casapool % ratioPcsoil    => sum_casapool_data % ratioPcsoil     
sum_casapool % ratioPcplant   => sum_casapool_data % ratioPcplant    
sum_casapool % ratioPclitter  => sum_casapool_data % ratioPclitter   
sum_casapool % cwoodprod      => sum_casapool_data % cwoodprod       
sum_casapool % nwoodprod      => sum_casapool_data % nwoodprod       
sum_casapool % pwoodprod      => sum_casapool_data % pwoodprod       

RETURN
END SUBROUTINE assoc_sum_casa_pool_type

SUBROUTINE update_sum_casa_pool( sum_casapool, casapool, sum_now, average_now,    &
                            nsteps )

USE casa_pool_type_mod, ONLY: casa_pool_data_type
                                      
IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: sum_casapool
TYPE (casa_pool_data_type), INTENT(IN) :: casapool
LOGICAL, INTENT(IN) :: sum_now, average_now
INTEGER, INTENT(IN) :: nsteps

!local vars
INTEGER             :: rsteps

rsteps =   REAL( nsteps )

IF (sum_now) THEN

  sum_casapool%Clabile        = sum_casapool%Clabile   + casapool%Clabile
  sum_casapool%dClabiledt     = sum_casapool%Clabile   + casapool%Clabile
  sum_casapool%Cplant         = sum_casapool%Cplant    + casapool%Cplant
  sum_casapool%Nplant         = sum_casapool%Nplant    + casapool%Nplant
  sum_casapool%Pplant         = sum_casapool%Pplant    + casapool%Pplant
  sum_casapool%dCplantdt      = sum_casapool%dCplantdt + casapool%dCplantdt
  sum_casapool%dNplantdt      = sum_casapool%dNplantdt + casapool%dNplantdt
  sum_casapool%dPplantdt      = sum_casapool%dPplantdt + casapool%dPplantdt
  sum_casapool%ratioNCplant   = sum_casapool%ratioNCplant                      &
                                + casapool%ratioNCplant
  sum_casapool%ratioNPplant   = sum_casapool%ratioNPplant                      &
                                + casapool%ratioNPplant
  sum_casapool%Nsoilmin       = sum_casapool%Nsoilmin  + casapool%Nsoilmin
  sum_casapool%Psoillab       = sum_casapool%Psoillab  + casapool%Psoillab
  sum_casapool%Psoilsorb      = sum_casapool%Psoilsorb + casapool%Psoilsorb
                                
  sum_casapool%Psoilocc       = sum_casapool%Psoilocc     + casapool%Psoilocc
  sum_casapool%dNsoilmindt    = sum_casapool%dNsoilmindt  + casapool%dNsoilmindt
  sum_casapool%dPsoillabdt    = sum_casapool%dPsoillabdt  + casapool%dPsoillabdt
  sum_casapool%dPsoilsorbdt   = sum_casapool%dPsoilsorbdt                      &
                                + casapool%dPsoilsorbdt
  sum_casapool%dPsoiloccdt    = sum_casapool%dPsoiloccdt  + casapool%dPsoiloccdt
  sum_casapool%Clitter        = sum_casapool%Clitter      + casapool%Clitter
  sum_casapool%Nlitter        = sum_casapool%Nlitter      + casapool%Nlitter
  sum_casapool%Plitter        = sum_casapool%Plitter      + casapool%Plitter
  sum_casapool%dClitterdt     = sum_casapool%dClitterdt   + casapool%dClitterdt
  sum_casapool%dNlitterdt     = sum_casapool%dNlitterdt   + casapool%dNlitterdt
  sum_casapool%dPlitterdt     = sum_casapool%dPlitterdt   + casapool%dPlitterdt
  sum_casapool%ratioNClitter  = sum_casapool%ratioNClitter                     &
                                + casapool%ratioNClitter
  sum_casapool%ratioNPlitter  = sum_casapool%ratioNPlitter                     &
                                + casapool%ratioNPlitter
  sum_casapool%Csoil          = sum_casapool%Csoil         + casapool%Csoil
  sum_casapool%Nsoil          = sum_casapool%Nsoil         + casapool%Nsoil
  sum_casapool%Psoil          = sum_casapool%Psoil         + casapool%Psoil
  sum_casapool%dCsoildt       = sum_casapool%dCsoildt      + casapool%dCsoildt
  sum_casapool%dNsoildt       = sum_casapool%dNsoildt      + casapool%dNsoildt
  sum_casapool%dPsoildt       = sum_casapool%dPsoildt      + casapool%dPsoildt
  sum_casapool%ratioNCsoil    = sum_casapool%ratioNCsoil                       &
                                + casapool%ratioNCsoil
  sum_casapool%ratioNPsoil    = sum_casapool%ratioNPsoil                       &
                                + casapool%ratioNPsoil
  sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew                    &
                                + casapool%ratioNCsoilnew
  sum_casapool%ratioNCsoilmin = sum_casapool%ratioNCsoilmin                    &
                                + casapool%ratioNCsoilmin
  sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax                    &
                                + casapool%ratioNCsoilmax
  sum_casapool%ratioPcsoil    = sum_casapool%ratioPcsoil  + casapool%ratioPcsoil
  sum_casapool%ratioPcplant   = sum_casapool%ratioPcplant                      &
                                + casapool%ratioPcplant
  sum_casapool%ratioPclitter  = sum_casapool%ratioPclitter                     &
                                + casapool%ratioPclitter

ENDIF

IF (average_now) THEN
  
  sum_casapool%Clabile        = sum_casapool%Clabile        / rsteps
  sum_casapool%dClabiledt     = sum_casapool%Clabile        / rsteps
  sum_casapool%Nplant         = sum_casapool%Nplant         / rsteps
  sum_casapool%Pplant         = sum_casapool%Pplant         / rsteps
  sum_casapool%dCplantdt      = sum_casapool%dCplantdt      / rsteps
  sum_casapool%dNplantdt      = sum_casapool%dNplantdt      / rsteps
  sum_casapool%dPplantdt      = sum_casapool%dPplantdt      / rsteps
  sum_casapool%ratioNCplant   = sum_casapool%ratioNCplant   / rsteps
  sum_casapool%ratioNPplant   = sum_casapool%ratioNPplant   / rsteps
  sum_casapool%Nsoilmin       = sum_casapool%Nsoilmin       / rsteps
  sum_casapool%Psoillab       = sum_casapool%Psoillab       / rsteps
  sum_casapool%Psoilsorb      = sum_casapool%Psoilsorb      / rsteps
  sum_casapool%Psoilocc       = sum_casapool%Psoilocc       / rsteps
  sum_casapool%dNsoilmindt    = sum_casapool%dNsoilmindt    / rsteps
  sum_casapool%dPsoillabdt    = sum_casapool%dPsoillabdt    / rsteps
  sum_casapool%dPsoilsorbdt   = sum_casapool%dPsoilsorbdt   / rsteps
  sum_casapool%dPsoiloccdt    = sum_casapool%dPsoiloccdt    / rsteps
  sum_casapool%Clitter        = sum_casapool%Clitter        / rsteps
  sum_casapool%Nlitter        = sum_casapool%Nlitter        / rsteps
  sum_casapool%Plitter        = sum_casapool%Plitter        / rsteps
  sum_casapool%dClitterdt     = sum_casapool%dClitterdt     / rsteps
  sum_casapool%dNlitterdt     = sum_casapool%dNlitterdt     / rsteps
  sum_casapool%dPlitterdt     = sum_casapool%dPlitterdt     / rsteps
  sum_casapool%ratioNClitter  = sum_casapool%ratioNClitter  / rsteps
  sum_casapool%ratioNPlitter  = sum_casapool%ratioNPlitter  / rsteps
  sum_casapool%Csoil          = sum_casapool%Csoil          / rsteps
  sum_casapool%Nsoil          = sum_casapool%Nsoil          / rsteps
  sum_casapool%Psoil          = sum_casapool%Psoil          / rsteps
  sum_casapool%dCsoildt       = sum_casapool%dCsoildt       / rsteps
  sum_casapool%dNsoildt       = sum_casapool%dNsoildt       / rsteps
  sum_casapool%dPsoildt       = sum_casapool%dPsoildt       / rsteps
  sum_casapool%ratioNCsoil    = sum_casapool%ratioNCsoil    / rsteps
  sum_casapool%ratioNPsoil    = sum_casapool%ratioNPsoil    / rsteps
  sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew / rsteps
  sum_casapool%ratioNCsoilmin = sum_casapool%ratioNCsoilmin / rsteps
  sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax / rsteps
  sum_casapool%ratioPcsoil    = sum_casapool%ratioPcsoil    / rsteps
  sum_casapool%ratioPcplant   = sum_casapool%ratioPcplant   / rsteps
  sum_casapool%ratioPclitter  = sum_casapool%ratioPclitter  / rsteps

ENDIF

RETURN
END SUBROUTINE update_sum_casa_pool

END MODULE sum_casa_pool_type_mod
