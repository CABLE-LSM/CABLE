MODULE casa_balance_type_mod

USE cable_other_constants_mod, ONLY: r_2 ! currently DOUBLE precision

IMPLICIT NONE

PUBLIC :: casa_bal_data_type
PUBLIC :: casa_bal_type
PUBLIC :: alloc_casa_bal_data_type
PUBLIC :: zero_casa_bal_data_type
PUBLIC :: assoc_casa_bal_type

TYPE casa_bal_data_type

  REAL(r_2), ALLOCATABLE :: FCgppyear     (:)
  REAL(r_2), ALLOCATABLE :: FCnppyear     (:)
  REAL(r_2), ALLOCATABLE :: FCrmleafyear  (:)
  REAL(r_2), ALLOCATABLE :: FCrmwoodyear  (:)
  REAL(r_2), ALLOCATABLE :: FCrmrootyear  (:)
  REAL(r_2), ALLOCATABLE :: FCrgrowyear   (:)
  REAL(r_2), ALLOCATABLE :: FCrpyear      (:)
  REAL(r_2), ALLOCATABLE :: FCrsyear      (:)
  REAL(r_2), ALLOCATABLE :: FCneeyear     (:)
  REAL(r_2), ALLOCATABLE :: dCdtyear      (:)
  REAL(r_2), ALLOCATABLE :: LAImax        (:)
  REAL(r_2), ALLOCATABLE :: Cleafmean     (:)
  REAL(r_2), ALLOCATABLE :: Crootmean     (:)
  REAL(r_2), ALLOCATABLE :: FNdepyear     (:)
  REAL(r_2), ALLOCATABLE :: FNfixyear     (:)
  REAL(r_2), ALLOCATABLE :: FNsnetyear    (:)
  REAL(r_2), ALLOCATABLE :: FNupyear      (:)
  REAL(r_2), ALLOCATABLE :: FNleachyear   (:)
  REAL(r_2), ALLOCATABLE :: FNlossyear    (:)
  REAL(r_2), ALLOCATABLE :: FPweayear     (:)
  REAL(r_2), ALLOCATABLE :: FPdustyear    (:)
  REAL(r_2), ALLOCATABLE :: FPsnetyear    (:)
  REAL(r_2), ALLOCATABLE :: FPupyear      (:)
  REAL(r_2), ALLOCATABLE :: FPleachyear   (:)
  REAL(r_2), ALLOCATABLE :: FPlossyear    (:)
  REAL(r_2), ALLOCATABLE :: nsoilminlast  (:)
  REAL(r_2), ALLOCATABLE :: psoillablast  (:)
  REAL(r_2), ALLOCATABLE :: psoilsorblast (:)
  REAL(r_2), ALLOCATABLE :: psoilocclast  (:)
  REAL(r_2), ALLOCATABLE :: cbalance      (:)
  REAL(r_2), ALLOCATABLE :: nbalance      (:)
  REAL(r_2), ALLOCATABLE :: pbalance      (:)
  REAL(r_2), ALLOCATABLE :: sumcbal       (:)
  REAL(r_2), ALLOCATABLE :: sumnbal       (:)
  REAL(r_2), ALLOCATABLE :: sumpbal       (:)
  REAL(r_2), ALLOCATABLE :: clabilelast   (:)
  REAL(r_2), ALLOCATABLE :: glaimon       (:,:)
  REAL(r_2), ALLOCATABLE :: glaimonx      (:,:)
  REAL(r_2), ALLOCATABLE :: cplantlast    (:,:)
  REAL(r_2), ALLOCATABLE :: nplantlast    (:,:)
  REAL(r_2), ALLOCATABLE :: pplantlast    (:,:)
  REAL(r_2), ALLOCATABLE :: clitterlast   (:,:)
  REAL(r_2), ALLOCATABLE :: nlitterlast   (:,:)
  REAL(r_2), ALLOCATABLE :: plitterlast   (:,:)
  REAL(r_2), ALLOCATABLE :: csoillast     (:,:)
  REAL(r_2), ALLOCATABLE :: nsoillast     (:,:)
  REAL(r_2), ALLOCATABLE :: psoillast     (:,:)
 
END TYPE casa_bal_data_type

TYPE casa_bal_type

  REAL(r_2), POINTER, PUBLIC :: FCgppyear     (:)
  REAL(r_2), POINTER, PUBLIC :: FCnppyear     (:)
  REAL(r_2), POINTER, PUBLIC :: FCrmleafyear  (:)
  REAL(r_2), POINTER, PUBLIC :: FCrmwoodyear  (:)
  REAL(r_2), POINTER, PUBLIC :: FCrmrootyear  (:)
  REAL(r_2), POINTER, PUBLIC :: FCrgrowyear   (:)
  REAL(r_2), POINTER, PUBLIC :: FCrpyear      (:)
  REAL(r_2), POINTER, PUBLIC :: FCrsyear      (:)
  REAL(r_2), POINTER, PUBLIC :: FCneeyear     (:)
  REAL(r_2), POINTER, PUBLIC :: dCdtyear      (:)
  REAL(r_2), POINTER, PUBLIC :: LAImax        (:)
  REAL(r_2), POINTER, PUBLIC :: Cleafmean     (:)
  REAL(r_2), POINTER, PUBLIC :: Crootmean     (:)
  REAL(r_2), POINTER, PUBLIC :: FNdepyear     (:)
  REAL(r_2), POINTER, PUBLIC :: FNfixyear     (:)
  REAL(r_2), POINTER, PUBLIC :: FNsnetyear    (:)
  REAL(r_2), POINTER, PUBLIC :: FNupyear      (:)
  REAL(r_2), POINTER, PUBLIC :: FNleachyear   (:)
  REAL(r_2), POINTER, PUBLIC :: FNlossyear    (:)
  REAL(r_2), POINTER, PUBLIC :: FPweayear     (:)
  REAL(r_2), POINTER, PUBLIC :: FPdustyear    (:)
  REAL(r_2), POINTER, PUBLIC :: FPsnetyear    (:)
  REAL(r_2), POINTER, PUBLIC :: FPupyear      (:)
  REAL(r_2), POINTER, PUBLIC :: FPleachyear   (:)
  REAL(r_2), POINTER, PUBLIC :: FPlossyear    (:)
  REAL(r_2), POINTER, PUBLIC :: nsoilminlast  (:)
  REAL(r_2), POINTER, PUBLIC :: psoillablast  (:)
  REAL(r_2), POINTER, PUBLIC :: psoilsorblast (:)
  REAL(r_2), POINTER, PUBLIC :: psoilocclast  (:)
  REAL(r_2), POINTER, PUBLIC :: cbalance      (:)
  REAL(r_2), POINTER, PUBLIC :: nbalance      (:)
  REAL(r_2), POINTER, PUBLIC :: pbalance      (:)
  REAL(r_2), POINTER, PUBLIC :: sumcbal       (:)
  REAL(r_2), POINTER, PUBLIC :: sumnbal       (:)
  REAL(r_2), POINTER, PUBLIC :: sumpbal       (:)
  REAL(r_2), POINTER, PUBLIC :: clabilelast   (:)
  REAL(r_2), POINTER, PUBLIC :: glaimon       (:,:)
  REAL(r_2), POINTER, PUBLIC :: glaimonx      (:,:)
  REAL(r_2), POINTER, PUBLIC :: cplantlast    (:,:)
  REAL(r_2), POINTER, PUBLIC :: nplantlast    (:,:)
  REAL(r_2), POINTER, PUBLIC :: pplantlast    (:,:)
  REAL(r_2), POINTER, PUBLIC :: clitterlast   (:,:)
  REAL(r_2), POINTER, PUBLIC :: nlitterlast   (:,:)
  REAL(r_2), POINTER, PUBLIC :: plitterlast   (:,:)
  REAL(r_2), POINTER, PUBLIC :: csoillast     (:,:)
  REAL(r_2), POINTER, PUBLIC :: nsoillast     (:,:)
  REAL(r_2), POINTER, PUBLIC :: psoillast     (:,:)
 
END TYPE casa_bal_type

CONTAINS

SUBROUTINE alloc_casa_bal_data_type( casabal_data, arraysize ) 

USE casadimension, ONLY: mplant, mlitter, msoil

IMPLICIT NONE

TYPE (casa_bal_data_type), INTENT(INOUT) :: casabal_data 
INTEGER,                   INTENT(IN)    :: arraysize

ALLOCATE ( casabal_data % FCgppyear     ( arraysize ) )
ALLOCATE ( casabal_data % FCnppyear     ( arraysize ) )
ALLOCATE ( casabal_data % FCrmleafyear  ( arraysize ) )
ALLOCATE ( casabal_data % FCrmwoodyear  ( arraysize ) )
ALLOCATE ( casabal_data % FCrmrootyear  ( arraysize ) )
ALLOCATE ( casabal_data % FCrgrowyear   ( arraysize ) )
ALLOCATE ( casabal_data % FCrpyear      ( arraysize ) )
ALLOCATE ( casabal_data % FCrsyear      ( arraysize ) )
ALLOCATE ( casabal_data % FCneeyear     ( arraysize ) )
ALLOCATE ( casabal_data % dCdtyear      ( arraysize ) )
ALLOCATE ( casabal_data % LAImax        ( arraysize ) )
ALLOCATE ( casabal_data % Cleafmean     ( arraysize ) )
ALLOCATE ( casabal_data % Crootmean     ( arraysize ) )
ALLOCATE ( casabal_data % FNdepyear     ( arraysize ) )
ALLOCATE ( casabal_data % FNfixyear     ( arraysize ) )
ALLOCATE ( casabal_data % FNsnetyear    ( arraysize ) )
ALLOCATE ( casabal_data % FNupyear      ( arraysize ) )
ALLOCATE ( casabal_data % FNleachyear   ( arraysize ) )
ALLOCATE ( casabal_data % FNlossyear    ( arraysize ) )
ALLOCATE ( casabal_data % FPweayear     ( arraysize ) )
ALLOCATE ( casabal_data % FPdustyear    ( arraysize ) )
ALLOCATE ( casabal_data % FPsnetyear    ( arraysize ) )
ALLOCATE ( casabal_data % FPupyear      ( arraysize ) )
ALLOCATE ( casabal_data % FPleachyear   ( arraysize ) )
ALLOCATE ( casabal_data % FPlossyear    ( arraysize ) )
ALLOCATE ( casabal_data % nsoilminlast  ( arraysize ) )
ALLOCATE ( casabal_data % psoillablast  ( arraysize ) )
ALLOCATE ( casabal_data % psoilsorblast ( arraysize ) )
ALLOCATE ( casabal_data % psoilocclast  ( arraysize ) )
ALLOCATE ( casabal_data % cbalance      ( arraysize ) )
ALLOCATE ( casabal_data % nbalance      ( arraysize ) )
ALLOCATE ( casabal_data % pbalance      ( arraysize ) )
ALLOCATE ( casabal_data % sumcbal       ( arraysize ) )
ALLOCATE ( casabal_data % sumnbal       ( arraysize ) )
ALLOCATE ( casabal_data % sumpbal       ( arraysize ) )
ALLOCATE ( casabal_data % clabilelast   ( arraysize ) )
ALLOCATE ( casabal_data % glaimon       ( arraysize, 12 ) )       ! 12 months
ALLOCATE ( casabal_data % glaimonx      ( arraysize, 12 ) )       ! 12 months
ALLOCATE ( casabal_data % cplantlast    ( arraysize, mplant  ) )
ALLOCATE ( casabal_data % nplantlast    ( arraysize, mplant  ) )
ALLOCATE ( casabal_data % pplantlast    ( arraysize, mplant  ) )
ALLOCATE ( casabal_data % clitterlast   ( arraysize, mlitter ) )
ALLOCATE ( casabal_data % nlitterlast   ( arraysize, mlitter ) )
ALLOCATE ( casabal_data % plitterlast   ( arraysize, mlitter ) )
ALLOCATE ( casabal_data % csoillast     ( arraysize, msoil   ) )
ALLOCATE ( casabal_data % nsoillast     ( arraysize, msoil   ) )
ALLOCATE ( casabal_data % psoillast     ( arraysize, msoil   ) )

END SUBROUTINE alloc_casa_bal_data_type

SUBROUTINE zero_casa_bal_data_type( casabal_data ) 

IMPLICIT NONE

TYPE (casa_bal_data_type), INTENT(INOUT) :: casabal_data 

casabal_data % FCgppyear     (:)     = 0.0
casabal_data % FCnppyear     (:)     = 0.0
casabal_data % FCrmleafyear  (:)     = 0.0
casabal_data % FCrmwoodyear  (:)     = 0.0
casabal_data % FCrmrootyear  (:)     = 0.0
casabal_data % FCrgrowyear   (:)     = 0.0
casabal_data % FCrpyear      (:)     = 0.0
casabal_data % FCrsyear      (:)     = 0.0
casabal_data % FCneeyear     (:)     = 0.0
casabal_data % dCdtyear      (:)     = 0.0
casabal_data % LAImax        (:)     = 0.0
casabal_data % Cleafmean     (:)     = 0.0
casabal_data % Crootmean     (:)     = 0.0
casabal_data % FNdepyear     (:)     = 0.0
casabal_data % FNfixyear     (:)     = 0.0
casabal_data % FNsnetyear    (:)     = 0.0
casabal_data % FNupyear      (:)     = 0.0
casabal_data % FNleachyear   (:)     = 0.0
casabal_data % FNlossyear    (:)     = 0.0
casabal_data % FPweayear     (:)     = 0.0
casabal_data % FPdustyear    (:)     = 0.0
casabal_data % FPsnetyear    (:)     = 0.0
casabal_data % FPupyear      (:)     = 0.0
casabal_data % FPleachyear   (:)     = 0.0
casabal_data % FPlossyear    (:)     = 0.0
casabal_data % nsoilminlast  (:)     = 0.0
casabal_data % psoillablast  (:)     = 0.0
casabal_data % psoilsorblast (:)     = 0.0
casabal_data % psoilocclast  (:)     = 0.0
casabal_data % cbalance      (:)     = 0.0
casabal_data % nbalance      (:)     = 0.0
casabal_data % pbalance      (:)     = 0.0
casabal_data % sumcbal       (:)     = 0.0
casabal_data % sumnbal       (:)     = 0.0
casabal_data % sumpbal       (:)     = 0.0
casabal_data % clabilelast   (:)     = 0.0
casabal_data % glaimon       (:,:)   = 0.0
casabal_data % glaimonx      (:,:)   = 0.0
casabal_data % cplantlast    (:,:)   = 0.0
casabal_data % nplantlast    (:,:)   = 0.0
casabal_data % pplantlast    (:,:)   = 0.0
casabal_data % clitterlast   (:,:)   = 0.0
casabal_data % nlitterlast   (:,:)   = 0.0
casabal_data % plitterlast   (:,:)   = 0.0
casabal_data % csoillast     (:,:)   = 0.0
casabal_data % nsoillast     (:,:)   = 0.0
casabal_data % psoillast     (:,:)   = 0.0

RETURN
END SUBROUTINE zero_casa_bal_data_type

SUBROUTINE assoc_casa_bal_type( casabal, casabal_data ) 

IMPLICIT NONE

TYPE (casa_bal_type),      INTENT(INOUT) :: casabal
TYPE (casa_bal_data_type), INTENT(INOUT), TARGET :: casabal_data

casabal % FCgppyear     => casabal_data % FCgppyear     
casabal % FCnppyear     => casabal_data % FCnppyear     
casabal % FCrmleafyear  => casabal_data % FCrmleafyear  
casabal % FCrmwoodyear  => casabal_data % FCrmwoodyear  
casabal % FCrmrootyear  => casabal_data % FCrmrootyear  
casabal % FCrgrowyear   => casabal_data % FCrgrowyear   
casabal % FCrpyear      => casabal_data % FCrpyear      
casabal % FCrsyear      => casabal_data % FCrsyear      
casabal % FCneeyear     => casabal_data % FCneeyear     
casabal % dCdtyear      => casabal_data % dCdtyear      
casabal % LAImax        => casabal_data % LAImax        
casabal % Cleafmean     => casabal_data % Cleafmean     
casabal % Crootmean     => casabal_data % Crootmean     
casabal % FNdepyear     => casabal_data % FNdepyear     
casabal % FNfixyear     => casabal_data % FNfixyear     
casabal % FNsnetyear    => casabal_data % FNsnetyear    
casabal % FNupyear      => casabal_data % FNupyear      
casabal % FNleachyear   => casabal_data % FNleachyear   
casabal % FNlossyear    => casabal_data % FNlossyear    
casabal % FPweayear     => casabal_data % FPweayear     
casabal % FPdustyear    => casabal_data % FPdustyear    
casabal % FPsnetyear    => casabal_data % FPsnetyear    
casabal % FPupyear      => casabal_data % FPupyear      
casabal % FPleachyear   => casabal_data % FPleachyear   
casabal % FPlossyear    => casabal_data % FPlossyear    
casabal % nsoilminlast  => casabal_data % nsoilminlast  
casabal % psoillablast  => casabal_data % psoillablast  
casabal % psoilsorblast => casabal_data % psoilsorblast 
casabal % psoilocclast  => casabal_data % psoilocclast  
casabal % cbalance      => casabal_data % cbalance      
casabal % nbalance      => casabal_data % nbalance      
casabal % pbalance      => casabal_data % pbalance      
casabal % sumcbal       => casabal_data % sumcbal       
casabal % sumnbal       => casabal_data % sumnbal       
casabal % sumpbal       => casabal_data % sumpbal       
casabal % clabilelast   => casabal_data % clabilelast   
casabal % glaimon       => casabal_data % glaimon       
casabal % glaimonx      => casabal_data % glaimonx      
casabal % cplantlast    => casabal_data % cplantlast    
casabal % nplantlast    => casabal_data % nplantlast    
casabal % pplantlast    => casabal_data % pplantlast    
casabal % clitterlast   => casabal_data % clitterlast   
casabal % nlitterlast   => casabal_data % nlitterlast   
casabal % plitterlast   => casabal_data % plitterlast   
casabal % csoillast     => casabal_data % csoillast     
casabal % nsoillast     => casabal_data % nsoillast     
casabal % psoillast     => casabal_data % psoillast     

RETURN
END SUBROUTINE assoc_casa_bal_type

END MODULE casa_balance_type_mod
