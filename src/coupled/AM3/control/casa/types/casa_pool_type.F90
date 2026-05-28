MODULE casa_pool_type_mod

USE cable_other_constants_mod, ONLY: r_2 ! currently DOUBLE precision

IMPLICIT NONE
 
PUBLIC :: casa_pool_data_type
PUBLIC :: casa_pool_type
PUBLIC :: alloc_casa_pool_data_type
PUBLIC :: assoc_casa_pool_type
PUBLIC :: zero_casa_pool_data_type

TYPE casa_pool_data_type

 REAL(r_2), ALLOCATABLE :: Clabile        (:)
 REAL(r_2), ALLOCATABLE :: dClabiledt     (:)
 REAL(r_2), ALLOCATABLE :: Ctot           (:)         
 REAL(r_2), ALLOCATABLE :: Ctot_0         (:)         
 REAL(r_2), ALLOCATABLE :: Nsoilmin       (:)
 REAL(r_2), ALLOCATABLE :: Psoillab       (:)
 REAL(r_2), ALLOCATABLE :: Psoilsorb      (:)
 REAL(r_2), ALLOCATABLE :: Psoilocc       (:)
 REAL(r_2), ALLOCATABLE :: dNsoilmindt    (:)
 REAL(r_2), ALLOCATABLE :: dPsoillabdt    (:)
 REAL(r_2), ALLOCATABLE :: dPsoilsorbdt   (:)
 REAL(r_2), ALLOCATABLE :: dPsoiloccdt    (:)      

 REAL(r_2), ALLOCATABLE :: Cplant         (:,:) 
 REAL(r_2), ALLOCATABLE :: Nplant         (:,:) 
 REAL(r_2), ALLOCATABLE :: Pplant         (:,:) 
 REAL(r_2), ALLOCATABLE :: dCplantdt      (:,:)    
 REAL(r_2), ALLOCATABLE :: dNplantdt      (:,:)  
 REAL(r_2), ALLOCATABLE :: dPplantdt      (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNCplant   (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNPplant   (:,:)  
 REAL(r_2), ALLOCATABLE :: Clitter        (:,:)  
 REAL(r_2), ALLOCATABLE :: Nlitter        (:,:)       
 REAL(r_2), ALLOCATABLE :: Plitter        (:,:)  
 REAL(r_2), ALLOCATABLE :: dClitterdt     (:,:)  
 REAL(r_2), ALLOCATABLE :: dNlitterdt     (:,:)  
 REAL(r_2), ALLOCATABLE :: dPlitterdt     (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNClitter  (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNPlitter  (:,:)  
 REAL(r_2), ALLOCATABLE :: Csoil          (:,:)  
 REAL(r_2), ALLOCATABLE :: Nsoil          (:,:)  
 REAL(r_2), ALLOCATABLE :: Psoil          (:,:)  
 REAL(r_2), ALLOCATABLE :: dCsoildt       (:,:)  
 REAL(r_2), ALLOCATABLE :: dNsoildt       (:,:)  
 REAL(r_2), ALLOCATABLE :: dPsoildt       (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNCsoil    (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNCsoilnew (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNPsoil    (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNCsoilmin (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioNCsoilmax (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioPcsoil    (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioPcplant   (:,:)  
 REAL(r_2), ALLOCATABLE :: ratioPclitter  (:,:)  
 REAL(r_2), ALLOCATABLE :: cwoodprod      (:,:)  
 REAL(r_2), ALLOCATABLE :: nwoodprod      (:,:)  
 REAL(r_2), ALLOCATABLE :: pwoodprod      (:,:)  

END TYPE casa_pool_data_type

TYPE casa_pool_type

 REAL(r_2), POINTER, PUBLIC :: Clabile        (:)
 REAL(r_2), POINTER, PUBLIC :: dClabiledt     (:)
 REAL(r_2), POINTER, PUBLIC :: Ctot           (:)         
 REAL(r_2), POINTER, PUBLIC :: Ctot_0         (:)         
 REAL(r_2), POINTER, PUBLIC :: Nsoilmin       (:)
 REAL(r_2), POINTER, PUBLIC :: Psoillab       (:)
 REAL(r_2), POINTER, PUBLIC :: Psoilsorb      (:)
 REAL(r_2), POINTER, PUBLIC :: Psoilocc       (:)
 REAL(r_2), POINTER, PUBLIC :: dNsoilmindt    (:)
 REAL(r_2), POINTER, PUBLIC :: dPsoillabdt    (:)
 REAL(r_2), POINTER, PUBLIC :: dPsoilsorbdt   (:)
 REAL(r_2), POINTER, PUBLIC :: dPsoiloccdt    (:)      

 REAL(r_2), POINTER, PUBLIC :: Cplant         (:,:) 
 REAL(r_2), POINTER, PUBLIC :: Nplant         (:,:) 
 REAL(r_2), POINTER, PUBLIC :: Pplant         (:,:) 
 REAL(r_2), POINTER, PUBLIC :: dCplantdt      (:,:)    
 REAL(r_2), POINTER, PUBLIC :: dNplantdt      (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dPplantdt      (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNCplant   (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNPplant   (:,:)  
 REAL(r_2), POINTER, PUBLIC :: Clitter        (:,:)  
 REAL(r_2), POINTER, PUBLIC :: Nlitter        (:,:)       
 REAL(r_2), POINTER, PUBLIC :: Plitter        (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dClitterdt     (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dNlitterdt     (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dPlitterdt     (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNClitter  (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNPlitter  (:,:)  
 REAL(r_2), POINTER, PUBLIC :: Csoil          (:,:)  
 REAL(r_2), POINTER, PUBLIC :: Nsoil          (:,:)  
 REAL(r_2), POINTER, PUBLIC :: Psoil          (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dCsoildt       (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dNsoildt       (:,:)  
 REAL(r_2), POINTER, PUBLIC :: dPsoildt       (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNCsoil    (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNCsoilnew (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNPsoil    (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNCsoilmin (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioNCsoilmax (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioPcsoil    (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioPcplant   (:,:)  
 REAL(r_2), POINTER, PUBLIC :: ratioPclitter  (:,:)  
 REAL(r_2), POINTER, PUBLIC :: cwoodprod      (:,:)  
 REAL(r_2), POINTER, PUBLIC :: nwoodprod      (:,:)  
 REAL(r_2), POINTER, PUBLIC :: pwoodprod      (:,:)  

END TYPE casa_pool_type

CONTAINS

SUBROUTINE alloc_casa_pool_data_type( casapool_data, arraysize ) 

USE casadimension, ONLY: mplant, mlitter, msoil, mwood

IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: casapool_data 
INTEGER,                    INTENT(IN)    :: arraysize

ALLOCATE ( casapool_data % Clabile        ( arraysize ) )
ALLOCATE ( casapool_data % dClabiledt     ( arraysize ) )
ALLOCATE ( casapool_data % Ctot           ( arraysize ) )
ALLOCATE ( casapool_data % Ctot_0         ( arraysize ) )
ALLOCATE ( casapool_data % Nsoilmin       ( arraysize ) )
ALLOCATE ( casapool_data % Psoillab       ( arraysize ) )
ALLOCATE ( casapool_data % Psoilsorb      ( arraysize ) )
ALLOCATE ( casapool_data % Psoilocc       ( arraysize ) )
ALLOCATE ( casapool_data % dNsoilmindt    ( arraysize ) )
ALLOCATE ( casapool_data % dPsoillabdt    ( arraysize ) )
ALLOCATE ( casapool_data % dPsoilsorbdt   ( arraysize ) )
ALLOCATE ( casapool_data % dPsoiloccdt    ( arraysize ) )

ALLOCATE ( casapool_data % Cplant         ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % Nplant         ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % Pplant         ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % dCplantdt      ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % dNplantdt      ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % dPplantdt      ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % ratioNCplant   ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % ratioNPplant   ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % Clitter        ( arraysize, mlitter ) ) 
ALLOCATE ( casapool_data % Nlitter        ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % Plitter        ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % dClitterdt     ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % dNlitterdt     ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % dPlitterdt     ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % ratioNClitter  ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % ratioNPlitter  ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % Csoil          ( arraysize, msoil ) )   
ALLOCATE ( casapool_data % Nsoil          ( arraysize, msoil ) )    
ALLOCATE ( casapool_data % Psoil          ( arraysize, msoil ) )    
ALLOCATE ( casapool_data % dCsoildt       ( arraysize, msoil ) ) 
ALLOCATE ( casapool_data % dNsoildt       ( arraysize, msoil ) ) 
ALLOCATE ( casapool_data % dPsoildt       ( arraysize, msoil ) ) 
ALLOCATE ( casapool_data % ratioNCsoil    ( arraysize, msoil ) ) 
ALLOCATE ( casapool_data % ratioNCsoilnew ( arraysize, msoil ) ) 
ALLOCATE ( casapool_data % ratioNPsoil    ( arraysize, msoil ) ) 
ALLOCATE ( casapool_data % ratioNCsoilmin ( arraysize, msoil ) )
ALLOCATE ( casapool_data % ratioNCsoilmax ( arraysize, msoil ) )
ALLOCATE ( casapool_data % ratioPcsoil    ( arraysize, msoil ) )
ALLOCATE ( casapool_data % ratioPcplant   ( arraysize, mplant ) ) 
ALLOCATE ( casapool_data % ratioPclitter  ( arraysize, mlitter ) )  
ALLOCATE ( casapool_data % cwoodprod      ( arraysize, mwood ) ) 
ALLOCATE ( casapool_data % nwoodprod      ( arraysize, mwood ) ) 
ALLOCATE ( casapool_data % pwoodprod      ( arraysize, mwood ) ) 

END SUBROUTINE alloc_casa_pool_data_type

SUBROUTINE zero_casa_pool_data_type( casapool_data ) 

IMPLICIT NONE

TYPE (casa_pool_data_type), INTENT(INOUT) :: casapool_data 

casapool_data % Clabile        = 0.0 
casapool_data % dClabiledt     = 0.0 
casapool_data % Ctot           = 0.0 
casapool_data % Ctot_0         = 0.0 
casapool_data % Nsoilmin       = 0.0 
casapool_data % Psoillab       = 0.0 
casapool_data % Psoilsorb      = 0.0 
casapool_data % Psoilocc       = 0.0 
casapool_data % dNsoilmindt    = 0.0 
casapool_data % dPsoillabdt    = 0.0 
casapool_data % dPsoilsorbdt   = 0.0 
casapool_data % dPsoiloccdt    = 0.0 
                              
casapool_data % Cplant         = 0.0 
casapool_data % Nplant         = 0.0 
casapool_data % Pplant         = 0.0 
casapool_data % dCplantdt      = 0.0 
casapool_data % dNplantdt      = 0.0 
casapool_data % dPplantdt      = 0.0 
casapool_data % ratioNCplant   = 0.0 
casapool_data % ratioNPplant   = 0.0 
casapool_data % Clitter        = 0.0 
casapool_data % Nlitter        = 0.0 
casapool_data % Plitter        = 0.0 
casapool_data % dClitterdt     = 0.0 
casapool_data % dNlitterdt     = 0.0 
casapool_data % dPlitterdt     = 0.0 
casapool_data % ratioNClitter  = 0.0 
casapool_data % ratioNPlitter  = 0.0 
casapool_data % Csoil          = 0.0 
casapool_data % Nsoil          = 0.0 
casapool_data % Psoil          = 0.0 
casapool_data % dCsoildt       = 0.0 
casapool_data % dNsoildt       = 0.0 
casapool_data % dPsoildt       = 0.0 
casapool_data % ratioNCsoil    = 0.0 
casapool_data % ratioNCsoilnew = 0.0 
casapool_data % ratioNPsoil    = 0.0 
casapool_data % ratioNCsoilmin = 0.0 
casapool_data % ratioNCsoilmax = 0.0 
casapool_data % ratioPcsoil    = 0.0 
casapool_data % ratioPcplant   = 0.0 
casapool_data % ratioPclitter  = 0.0 
casapool_data % cwoodprod      = 0.0 
casapool_data % nwoodprod      = 0.0 
casapool_data % pwoodprod      = 0.0 

RETURN
END SUBROUTINE zero_casa_pool_data_type

SUBROUTINE assoc_casa_pool_type( casapool, casapool_data ) 

IMPLICIT NONE

TYPE (casa_pool_type),      INTENT(INOUT) :: casapool
TYPE (casa_pool_data_type), INTENT(INOUT), TARGET :: casapool_data

casapool % Clabile        => casapool_data % Clabile         
casapool % dClabiledt     => casapool_data % dClabiledt      
casapool % Ctot           => casapool_data % Ctot            
casapool % Ctot_0         => casapool_data % Ctot_0          
casapool % Nsoilmin       => casapool_data % Nsoilmin        
casapool % Psoillab       => casapool_data % Psoillab        
casapool % Psoilsorb      => casapool_data % Psoilsorb       
casapool % Psoilocc       => casapool_data % Psoilocc        
casapool % dNsoilmindt    => casapool_data % dNsoilmindt     
casapool % dPsoillabdt    => casapool_data % dPsoillabdt     
casapool % dPsoilsorbdt   => casapool_data % dPsoilsorbdt    
casapool % dPsoiloccdt    => casapool_data % dPsoiloccdt     

casapool % Cplant         => casapool_data % Cplant          
casapool % Nplant         => casapool_data % Nplant          
casapool % Pplant         => casapool_data % Pplant          
casapool % dCplantdt      => casapool_data % dCplantdt       
casapool % dNplantdt      => casapool_data % dNplantdt       
casapool % dPplantdt      => casapool_data % dPplantdt       
casapool % ratioNCplant   => casapool_data % ratioNCplant    
casapool % ratioNPplant   => casapool_data % ratioNPplant    
casapool % Clitter        => casapool_data % Clitter         
casapool % Nlitter        => casapool_data % Nlitter         
casapool % Plitter        => casapool_data % Plitter         
casapool % dClitterdt     => casapool_data % dClitterdt      
casapool % dNlitterdt     => casapool_data % dNlitterdt      
casapool % dPlitterdt     => casapool_data % dPlitterdt      
casapool % ratioNClitter  => casapool_data % ratioNClitter   
casapool % ratioNPlitter  => casapool_data % ratioNPlitter   
casapool % Csoil          => casapool_data % Csoil           
casapool % Nsoil          => casapool_data % Nsoil           
casapool % Psoil          => casapool_data % Psoil           
casapool % dCsoildt       => casapool_data % dCsoildt        
casapool % dNsoildt       => casapool_data % dNsoildt        
casapool % dPsoildt       => casapool_data % dPsoildt        
casapool % ratioNCsoil    => casapool_data % ratioNCsoil     
casapool % ratioNCsoilnew => casapool_data % ratioNCsoilnew  
casapool % ratioNPsoil    => casapool_data % ratioNPsoil     
casapool % ratioNCsoilmin => casapool_data % ratioNCsoilmin  
casapool % ratioNCsoilmax => casapool_data % ratioNCsoilmax  
casapool % ratioPcsoil    => casapool_data % ratioPcsoil     
casapool % ratioPcplant   => casapool_data % ratioPcplant    
casapool % ratioPclitter  => casapool_data % ratioPclitter   
casapool % cwoodprod      => casapool_data % cwoodprod       
casapool % nwoodprod      => casapool_data % nwoodprod       
casapool % pwoodprod      => casapool_data % pwoodprod       

RETURN
END SUBROUTINE assoc_casa_pool_type

END MODULE casa_pool_type_mod
