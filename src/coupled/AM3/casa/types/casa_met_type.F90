MODULE casa_met_type_mod

USE cable_other_constants_mod, ONLY: r_2 ! currently DOUBLE precision

IMPLICIT NONE

PUBLIC :: casa_met_data_type
PUBLIC :: casa_met_type
PUBLIC :: alloc_casa_met_data_type
PUBLIC :: assoc_casa_met_type

TYPE casa_met_data_type

  INTEGER, ALLOCATABLE   :: lnonwood       (:)
  INTEGER, ALLOCATABLE   :: iveg2          (:)
  INTEGER, ALLOCATABLE   :: ijgcm          (:)
  INTEGER, ALLOCATABLE   :: isorder        (:)

  REAL(r_2), ALLOCATABLE :: glai           (:)
  REAL(r_2), ALLOCATABLE :: Tairk          (:)
  REAL(r_2), ALLOCATABLE :: precip         (:)
  REAL(r_2), ALLOCATABLE :: tsoilavg       (:)
  REAL(r_2), ALLOCATABLE :: moistavg       (:)
  REAL(r_2), ALLOCATABLE :: btran          (:)
  REAL(r_2), ALLOCATABLE :: lat            (:)
  REAL(r_2), ALLOCATABLE :: lon            (:)
  REAL(r_2), ALLOCATABLE :: areacell       (:)
  
  REAL(r_2), ALLOCATABLE :: Tsoil          (:,:)
  REAL(r_2), ALLOCATABLE :: moist          (:,:)
  REAL(r_2), ALLOCATABLE :: Tairkspin      (:,:)
  REAL(r_2), ALLOCATABLE :: cgppspin       (:,:)
  REAL(r_2), ALLOCATABLE :: crmplantspin_1 (:,:) 
  REAL(r_2), ALLOCATABLE :: crmplantspin_2 (:,:) 
  REAL(r_2), ALLOCATABLE :: crmplantspin_3 (:,:) 
  REAL(r_2), ALLOCATABLE :: Tsoilspin_1    (:,:) 
  REAL(r_2), ALLOCATABLE :: Tsoilspin_2    (:,:) 
  REAL(r_2), ALLOCATABLE :: Tsoilspin_3    (:,:) 
  REAL(r_2), ALLOCATABLE :: Tsoilspin_4    (:,:) 
  REAL(r_2), ALLOCATABLE :: Tsoilspin_5    (:,:) 
  REAL(r_2), ALLOCATABLE :: Tsoilspin_6    (:,:) 
  REAL(r_2), ALLOCATABLE :: moistspin_1    (:,:) 
  REAL(r_2), ALLOCATABLE :: moistspin_2    (:,:) 
  REAL(r_2), ALLOCATABLE :: moistspin_3    (:,:) 
  REAL(r_2), ALLOCATABLE :: moistspin_4    (:,:) 
  REAL(r_2), ALLOCATABLE :: moistspin_5    (:,:) 
  REAL(r_2), ALLOCATABLE :: moistspin_6    (:,:) 
  REAL(r_2), ALLOCATABLE :: mtempspin      (:,:)

END TYPE casa_met_data_type

TYPE casa_met_type

  INTEGER, PUBLIC, POINTER   :: lnonwood       (:)
  INTEGER, PUBLIC, POINTER   :: iveg2          (:)
  INTEGER, PUBLIC, POINTER   :: ijgcm          (:)
  INTEGER, PUBLIC, POINTER   :: isorder        (:)

  REAL(r_2), PUBLIC, POINTER :: glai           (:)
  REAL(r_2), PUBLIC, POINTER :: Tairk          (:)
  REAL(r_2), PUBLIC, POINTER :: precip         (:)
  REAL(r_2), PUBLIC, POINTER :: tsoilavg       (:)
  REAL(r_2), PUBLIC, POINTER :: moistavg       (:)
  REAL(r_2), PUBLIC, POINTER :: btran          (:)
  REAL(r_2), PUBLIC, POINTER :: lat            (:)
  REAL(r_2), PUBLIC, POINTER :: lon            (:)
  REAL(r_2), PUBLIC, POINTER :: areacell       (:)
  
  REAL(r_2), PUBLIC, POINTER :: Tsoil          (:,:)
  REAL(r_2), PUBLIC, POINTER :: moist          (:,:)
  REAL(r_2), PUBLIC, POINTER :: Tairkspin      (:,:)
  REAL(r_2), PUBLIC, POINTER :: cgppspin       (:,:)
  REAL(r_2), PUBLIC, POINTER :: crmplantspin_1 (:,:) 
  REAL(r_2), PUBLIC, POINTER :: crmplantspin_2 (:,:) 
  REAL(r_2), PUBLIC, POINTER :: crmplantspin_3 (:,:) 
  REAL(r_2), PUBLIC, POINTER :: Tsoilspin_1    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: Tsoilspin_2    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: Tsoilspin_3    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: Tsoilspin_4    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: Tsoilspin_5    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: Tsoilspin_6    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: moistspin_1    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: moistspin_2    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: moistspin_3    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: moistspin_4    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: moistspin_5    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: moistspin_6    (:,:) 
  REAL(r_2), PUBLIC, POINTER :: mtempspin      (:,:)

END TYPE casa_met_type

CONTAINS

SUBROUTINE alloc_casa_met_data_type( casamet_data, arraysize ) 

USE casadimension,          ONLY: mdyear
USE grid_constants_mod_cbl, ONLY: nsl ! # of soil layers [6]

IMPLICIT NONE

TYPE (casa_met_data_type), INTENT(INOUT) :: casamet_data 
INTEGER,                   INTENT(IN)    :: arraysize

ALLOCATE( casamet_data % lnonwood       ( arraysize ) ) 
ALLOCATE( casamet_data % iveg2          ( arraysize ) ) 
ALLOCATE( casamet_data % ijgcm          ( arraysize ) ) 
ALLOCATE( casamet_data % isorder        ( arraysize ) )  
ALLOCATE( casamet_data % glai           ( arraysize ) ) 
ALLOCATE( casamet_data % Tairk          ( arraysize ) ) 
ALLOCATE( casamet_data % precip         ( arraysize ) ) 
ALLOCATE( casamet_data % tsoilavg       ( arraysize ) ) 
ALLOCATE( casamet_data % moistavg       ( arraysize ) ) 
ALLOCATE( casamet_data % btran          ( arraysize ) ) 
ALLOCATE( casamet_data % lat            ( arraysize ) )  
ALLOCATE( casamet_data % lon            ( arraysize ) )  
ALLOCATE( casamet_data % areacell       ( arraysize ) ) 

ALLOCATE( casamet_data % Tsoil          ( arraysize, nsl ) )
ALLOCATE( casamet_data % moist          ( arraysize, nsl ) )
ALLOCATE( casamet_data % Tairkspin      ( arraysize, mdyear ) )
ALLOCATE( casamet_data % cgppspin       ( arraysize, mdyear ) )    
ALLOCATE( casamet_data % crmplantspin_1 ( arraysize, mdyear ) )
ALLOCATE( casamet_data % crmplantspin_2 ( arraysize, mdyear ) )
ALLOCATE( casamet_data % crmplantspin_3 ( arraysize, mdyear ) )
ALLOCATE( casamet_data % Tsoilspin_1    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % Tsoilspin_2    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % Tsoilspin_3    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % Tsoilspin_4    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % Tsoilspin_5    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % Tsoilspin_6    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % moistspin_1    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % moistspin_2    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % moistspin_3    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % moistspin_4    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % moistspin_5    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % moistspin_6    ( arraysize, mdyear ) )
ALLOCATE( casamet_data % mtempspin      ( arraysize, mdyear ) )

RETURN
END SUBROUTINE alloc_casa_met_data_type


SUBROUTINE zero_casa_met_data_type( casamet_data ) 

IMPLICIT NONE

TYPE (casa_met_data_type), INTENT(INOUT) :: casamet_data 

casamet_data % lnonwood       (:)   = 0.0         
casamet_data % iveg2          (:)   = 0.0               
casamet_data % ijgcm          (:)   = 0.0                 
casamet_data % isorder        (:)   = 0.0              
casamet_data % glai           (:)   = 0.0   
casamet_data % Tairk          (:)   = 0.0   
casamet_data % precip         (:)   = 0.0   
casamet_data % tsoilavg       (:)   = 0.0  
casamet_data % moistavg       (:)   = 0.0  
casamet_data % btran          (:)   = 0.0  
casamet_data % lat            (:)   = 0.0  
casamet_data % lon            (:)   = 0.0  
casamet_data % areacell       (:)   = 0.0  
                              
casamet_data % Tsoil          (:,:) = 0.0 
casamet_data % moist          (:,:) = 0.0
casamet_data % Tairkspin      (:,:) = 0.0
casamet_data % cgppspin       (:,:) = 0.0
casamet_data % crmplantspin_1 (:,:) = 0.0          
casamet_data % crmplantspin_2 (:,:) = 0.0         
casamet_data % crmplantspin_3 (:,:) = 0.0          
casamet_data % Tsoilspin_1    (:,:) = 0.0
casamet_data % Tsoilspin_2    (:,:) = 0.0
casamet_data % Tsoilspin_3    (:,:) = 0.0
casamet_data % Tsoilspin_4    (:,:) = 0.0
casamet_data % Tsoilspin_5    (:,:) = 0.0
casamet_data % Tsoilspin_6    (:,:) = 0.0
casamet_data % moistspin_1    (:,:) = 0.0
casamet_data % moistspin_2    (:,:) = 0.0
casamet_data % moistspin_3    (:,:) = 0.0
casamet_data % moistspin_4    (:,:) = 0.0
casamet_data % moistspin_5    (:,:) = 0.0
casamet_data % moistspin_6    (:,:) = 0.0
casamet_data % mtempspin      (:,:) = 0.0

RETURN
END SUBROUTINE zero_casa_met_data_type


SUBROUTINE assoc_casa_met_type( casamet, casamet_data ) 

IMPLICIT NONE

TYPE (casa_met_type),      INTENT(INOUT) :: casamet
TYPE (casa_met_data_type), INTENT(INOUT), TARGET :: casamet_data

casamet % lnonwood       =>  casamet_data % lnonwood              
casamet % iveg2          =>  casamet_data % iveg2                       
casamet % ijgcm          =>  casamet_data % ijgcm                         
casamet % isorder        =>  casamet_data % isorder                    
casamet % glai           =>  casamet_data % glai            
casamet % Tairk          =>  casamet_data % Tairk           
casamet % precip         =>  casamet_data % precip          
casamet % tsoilavg       =>  casamet_data % tsoilavg       
casamet % moistavg       =>  casamet_data % moistavg       
casamet % btran          =>  casamet_data % btran          
casamet % lat            =>  casamet_data % lat            
casamet % lon            =>  casamet_data % lon            
casamet % areacell       =>  casamet_data % areacell       
                        
casamet % Tsoil          =>  casamet_data % Tsoil          
casamet % moist          =>  casamet_data % moist          
casamet % Tairkspin      =>  casamet_data % Tairkspin      
casamet % cgppspin       =>  casamet_data % cgppspin       
casamet % crmplantspin_1 =>  casamet_data % crmplantspin_1         
casamet % crmplantspin_2 =>  casamet_data % crmplantspin_2        
casamet % crmplantspin_3 =>  casamet_data % crmplantspin_3         
casamet % Tsoilspin_1    =>  casamet_data % Tsoilspin_1    
casamet % Tsoilspin_2    =>  casamet_data % Tsoilspin_2    
casamet % Tsoilspin_3    =>  casamet_data % Tsoilspin_3    
casamet % Tsoilspin_4    =>  casamet_data % Tsoilspin_4    
casamet % Tsoilspin_5    =>  casamet_data % Tsoilspin_5    
casamet % Tsoilspin_6    =>  casamet_data % Tsoilspin_6    
casamet % moistspin_1    =>  casamet_data % moistspin_1    
casamet % moistspin_2    =>  casamet_data % moistspin_2    
casamet % moistspin_3    =>  casamet_data % moistspin_3    
casamet % moistspin_4    =>  casamet_data % moistspin_4    
casamet % moistspin_5    =>  casamet_data % moistspin_5    
casamet % moistspin_6    =>  casamet_data % moistspin_6    
casamet % mtempspin      =>  casamet_data % mtempspin      

RETURN
END SUBROUTINE assoc_casa_met_type

END MODULE casa_met_type_mod

