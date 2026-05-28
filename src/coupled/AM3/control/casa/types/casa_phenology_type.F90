MODULE phenology_type_mod

USE cable_other_constants_mod, ONLY: r_2          ! currently DOUBLE precision
                                                  ! ntype_max=17
USE casadimension,             ONLY: mvtype, mphase, mdyear

IMPLICIT NONE

PUBLIC :: phenology_data_type
PUBLIC :: phenology_type
PUBLIC :: alloc_phenology_data_type
PUBLIC :: zero_phenology_data_type
PUBLIC :: assoc_phenology_type
                                     
TYPE phenology_data_type

  INTEGER,   ALLOCATABLE :: phase          (:)
  REAL(r_2), ALLOCATABLE :: TKshed         (:)
  INTEGER,   ALLOCATABLE :: doyphase       (:,:)
  REAL,      ALLOCATABLE :: phen           (:)    ! fraction of max LAI
  REAL,      ALLOCATABLE :: aphen          (:)    ! annual leaf on sum
  INTEGER,   ALLOCATABLE :: phasespin      (:,:)
  INTEGER,   ALLOCATABLE :: doyphasespin_1 (:,:)
  INTEGER,   ALLOCATABLE :: doyphasespin_2 (:,:)
  INTEGER,   ALLOCATABLE :: doyphasespin_3 (:,:)
  INTEGER,   ALLOCATABLE :: doyphasespin_4 (:,:)
             
END TYPE phenology_data_type
 
TYPE phenology_type

  INTEGER,   PUBLIC,  POINTER :: phase          (:)
  REAL(r_2), PUBLIC,  POINTER :: TKshed         (:)
  INTEGER,   PUBLIC,  POINTER :: doyphase       (:,:)
  REAL,      PUBLIC,  POINTER :: phen           (:)  
  REAL,      PUBLIC,  POINTER :: aphen          (:)  
  INTEGER,   PUBLIC,  POINTER :: phasespin      (:,:)
  INTEGER,   PUBLIC,  POINTER :: doyphasespin_1 (:,:)
  INTEGER,   PUBLIC,  POINTER :: doyphasespin_2 (:,:)
  INTEGER,   PUBLIC,  POINTER :: doyphasespin_3 (:,:)
  INTEGER,   PUBLIC,  POINTER :: doyphasespin_4 (:,:)
  
END TYPE phenology_type
 
CONTAINS

SUBROUTINE alloc_phenology_data_type( phen_data, arraysize )

IMPLICIT NONE

TYPE( phenology_data_type ), INTENT(INOUT) :: phen_data
INTEGER,                          INTENT(IN)    :: arraysize

ALLOCATE( phen_data % TKshed         ( mvtype            ) )
ALLOCATE( phen_data % phase          ( arraysize         ) )
ALLOCATE( phen_data % doyphase       ( arraysize, mphase ) )
ALLOCATE( phen_data % phen           ( arraysize         ) )
ALLOCATE( phen_data % aphen          ( arraysize         ) )
ALLOCATE( phen_data % phasespin      ( arraysize, mdyear ) )
ALLOCATE( phen_data % doyphasespin_1 ( arraysize, mdyear ) )
ALLOCATE( phen_data % doyphasespin_2 ( arraysize, mdyear ) )
ALLOCATE( phen_data % doyphasespin_3 ( arraysize, mdyear ) )
ALLOCATE( phen_data % doyphasespin_4 ( arraysize, mdyear ) )

RETURN
END SUBROUTINE alloc_phenology_data_type

SUBROUTINE zero_phenology_data_type( phen_data ) 

IMPLICIT NONE

TYPE (phenology_data_type), INTENT(INOUT) :: phen_data 

phen_data % phase          (:)   = 0
phen_data % doyphase       (:,:) = 0
phen_data % TKshed         (:)   = 0.0 
phen_data % phen           (:)   = 0.0   
phen_data % aphen          (:)   = 0.0   
phen_data % phasespin      (:,:) = 0.0   
phen_data % doyphasespin_1 (:,:) = 0.0   
phen_data % doyphasespin_2 (:,:) = 0.0   
phen_data % doyphasespin_3 (:,:) = 0.0   
phen_data % doyphasespin_4 (:,:) = 0.0   

RETURN
END SUBROUTINE zero_phenology_data_type

SUBROUTINE assoc_phenology_type( phen, phen_data ) 

IMPLICIT NONE

TYPE (phenology_type),      INTENT(INOUT)         :: phen
TYPE (phenology_data_type), INTENT(INOUT), TARGET :: phen_data

phen % phase          => phen_data % phase                 
phen % doyphase       => phen_data % doyphase              
phen % TKshed         => phen_data % TKshed                 
phen % phen           => phen_data % phen                     
phen % aphen          => phen_data % aphen                    
phen % phasespin      => phen_data % phasespin                
phen % doyphasespin_1 => phen_data % doyphasespin_1           
phen % doyphasespin_2 => phen_data % doyphasespin_2           
phen % doyphasespin_3 => phen_data % doyphasespin_3           
phen % doyphasespin_4 => phen_data % doyphasespin_4           

RETURN
END SUBROUTINE assoc_phenology_type

END MODULE phenology_type_mod

