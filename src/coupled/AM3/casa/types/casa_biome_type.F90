MODULE casa_biome_type_mod

USE cable_other_constants_mod, ONLY: r_2          ! currently DOUBLE precision
USE grid_constants_mod_cbl,    ONLY: ntype_max    ! # veg types [13],non-veg=4,
                                                  ! ntype_max=17
USE casadimension,             ONLY: mplant, mlitter, mso, msoil

IMPLICIT NONE

PUBLIC :: casa_biome_data_type
PUBLIC :: casa_biome_type
PUBLIC :: zero_casa_biome_data_type
PUBLIC :: assoc_casa_biome_type
                                     
TYPE casa_biome_data_type

  INTEGER   :: ivt2            ( ntype_max )
  REAL(r_2) :: xkleafcoldmax   ( ntype_max )
  REAL(r_2) :: xkleafcoldexp   ( ntype_max )
  REAL(r_2) :: xkleafdrymax    ( ntype_max )
  REAL(r_2) :: xkleafdryexp    ( ntype_max )
  REAL(r_2) :: glaimax         ( ntype_max )
  REAL(r_2) :: glaimin         ( ntype_max )
  REAL(r_2) :: sla             ( ntype_max )
  REAL(r_2) :: ratiofrootleaf  ( ntype_max )
  REAL(r_2) :: kroot           ( ntype_max )
  REAL(r_2) :: krootlen        ( ntype_max )
  REAL(r_2) :: rootdepth       ( ntype_max )
  REAL(r_2) :: kuptake         ( ntype_max )
  REAL(r_2) :: kminN           ( ntype_max )
  REAL(r_2) :: kuplabP         ( ntype_max )
  REAL(r_2) :: kclabrate       ( ntype_max )
  REAL(r_2) :: xnpmax          ( ntype_max )
  REAL(r_2) :: q10soil         ( ntype_max )
  REAL(r_2) :: xkoptlitter     ( ntype_max )
  REAL(r_2) :: xkoptsoil       ( ntype_max )
  REAL(r_2) :: xkplab          ( mso       )
  REAL(r_2) :: xkpsorb         ( mso       )
  REAL(r_2) :: xkpocc          ( mso       )
  REAL(r_2) :: prodptase       ( ntype_max )
  REAL(r_2) :: costnpup        ( ntype_max )
  REAL(r_2) :: maxfinelitter   ( ntype_max )
  REAL(r_2) :: maxcwd          ( ntype_max )
  REAL(r_2) :: nintercept      ( ntype_max )
  REAL(r_2) :: nslope          ( ntype_max )
  REAL(r_2) :: plantrate       ( ntype_max, mplant )
  REAL(r_2) :: rmplant         ( ntype_max, mplant )
  REAL(r_2) :: fracnpptoP      ( ntype_max, mplant )
  REAL(r_2) :: fraclignin      ( ntype_max, mplant )
  REAL(r_2) :: fraclabile      ( ntype_max, mplant )
  REAL(r_2) :: ratioNCplantmin ( ntype_max, mplant )
  REAL(r_2) :: ratioNCplantmax ( ntype_max, mplant )
  REAL(r_2) :: ratioNPplantmin ( ntype_max, mplant )
  REAL(r_2) :: ratioNPplantmax ( ntype_max, mplant )
  REAL(r_2) :: fracLigninplant ( ntype_max, mplant )
  REAL(r_2) :: ftransNPtoL     ( ntype_max, mplant )
  REAL(r_2) :: ftransPPtoL     ( ntype_max, mplant )
  REAL(r_2) :: litterrate      ( ntype_max, mlitter)
  REAL(r_2) :: ratioPcplantmin ( ntype_max, mplant )
  REAL(r_2) :: ratioPcplantmax ( ntype_max, mplant )
  REAL(r_2) :: soilrate        ( ntype_max, msoil  )

END TYPE casa_biome_data_type
 
TYPE casa_biome_type

  INTEGER, POINTER, PUBLIC   :: ivt2            (:)
  REAL(r_2), POINTER, PUBLIC :: xkleafcoldmax   (:)
  REAL(r_2), POINTER, PUBLIC :: xkleafcoldexp   (:)
  REAL(r_2), POINTER, PUBLIC :: xkleafdrymax    (:)
  REAL(r_2), POINTER, PUBLIC :: xkleafdryexp    (:)
  REAL(r_2), POINTER, PUBLIC :: glaimax         (:)
  REAL(r_2), POINTER, PUBLIC :: glaimin         (:)
  REAL(r_2), POINTER, PUBLIC :: sla             (:)
  REAL(r_2), POINTER, PUBLIC :: ratiofrootleaf  (:)
  REAL(r_2), POINTER, PUBLIC :: kroot           (:)
  REAL(r_2), POINTER, PUBLIC :: krootlen        (:)
  REAL(r_2), POINTER, PUBLIC :: rootdepth       (:)
  REAL(r_2), POINTER, PUBLIC :: kuptake         (:)
  REAL(r_2), POINTER, PUBLIC :: kminN           (:)
  REAL(r_2), POINTER, PUBLIC :: kuplabP         (:)
  REAL(r_2), POINTER, PUBLIC :: kclabrate       (:)
  REAL(r_2), POINTER, PUBLIC :: xnpmax          (:)
  REAL(r_2), POINTER, PUBLIC :: q10soil         (:)
  REAL(r_2), POINTER, PUBLIC :: xkoptlitter     (:)
  REAL(r_2), POINTER, PUBLIC :: xkoptsoil       (:)
  REAL(r_2), POINTER, PUBLIC :: xkplab          (:)
  REAL(r_2), POINTER, PUBLIC :: xkpsorb         (:)
  REAL(r_2), POINTER, PUBLIC :: xkpocc          (:)
  REAL(r_2), POINTER, PUBLIC :: prodptase       (:)
  REAL(r_2), POINTER, PUBLIC :: costnpup        (:)
  REAL(r_2), POINTER, PUBLIC :: maxfinelitter   (:)
  REAL(r_2), POINTER, PUBLIC :: maxcwd          (:)
  REAL(r_2), POINTER, PUBLIC :: nintercept      (:)
  REAL(r_2), POINTER, PUBLIC :: nslope          (:)
  REAL(r_2), POINTER, PUBLIC :: plantrate       (:,:)
  REAL(r_2), POINTER, PUBLIC :: rmplant         (:,:)
  REAL(r_2), POINTER, PUBLIC :: fracnpptoP      (:,:)
  REAL(r_2), POINTER, PUBLIC :: fraclignin      (:,:)
  REAL(r_2), POINTER, PUBLIC :: fraclabile      (:,:)
  REAL(r_2), POINTER, PUBLIC :: ratioNCplantmin (:,:)
  REAL(r_2), POINTER, PUBLIC :: ratioNCplantmax (:,:)
  REAL(r_2), POINTER, PUBLIC :: ratioNPplantmin (:,:)
  REAL(r_2), POINTER, PUBLIC :: ratioNPplantmax (:,:)
  REAL(r_2), POINTER, PUBLIC :: fracLigninplant (:,:)
  REAL(r_2), POINTER, PUBLIC :: ftransNPtoL     (:,:)
  REAL(r_2), POINTER, PUBLIC :: ftransPPtoL     (:,:)
  REAL(r_2), POINTER, PUBLIC :: litterrate      (:,:)
  REAL(r_2), POINTER, PUBLIC :: ratioPcplantmin (:,:)
  REAL(r_2), POINTER, PUBLIC :: ratioPcplantmax (:,:)

  REAL(r_2), POINTER, PUBLIC :: soilrate        (:,:)

END TYPE casa_biome_type
 
CONTAINS

SUBROUTINE zero_casa_biome_data_type( casabiome_data ) 

IMPLICIT NONE

TYPE (casa_biome_data_type), INTENT(INOUT) :: casabiome_data 

casabiome_data % ivt2            (:)   = 0
casabiome_data % xkleafcoldmax   (:)   = 0.0
casabiome_data % xkleafcoldexp   (:)   = 0.0
casabiome_data % xkleafdrymax    (:)   = 0.0
casabiome_data % xkleafdryexp    (:)   = 0.0
casabiome_data % glaimax         (:)   = 0.0
casabiome_data % glaimin         (:)   = 0.0
casabiome_data % sla             (:)   = 0.0
casabiome_data % ratiofrootleaf  (:)   = 0.0
casabiome_data % kroot           (:)   = 0.0
casabiome_data % krootlen        (:)   = 0.0
casabiome_data % rootdepth       (:)   = 0.0
casabiome_data % kuptake         (:)   = 0.0
casabiome_data % kminN           (:)   = 0.0
casabiome_data % kuplabP         (:)   = 0.0
casabiome_data % kclabrate       (:)   = 0.0
casabiome_data % xnpmax          (:)   = 0.0
casabiome_data % q10soil         (:)   = 0.0
casabiome_data % xkoptlitter     (:)   = 0.0
casabiome_data % xkoptsoil       (:)   = 0.0
casabiome_data % xkplab          (:)   = 0.0
casabiome_data % xkpsorb         (:)   = 0.0
casabiome_data % xkpocc          (:)   = 0.0
casabiome_data % prodptase       (:)   = 0.0
casabiome_data % costnpup        (:)   = 0.0
casabiome_data % maxfinelitter   (:)   = 0.0
casabiome_data % maxcwd          (:)   = 0.0
casabiome_data % nintercept      (:)   = 0.0
casabiome_data % nslope          (:)   = 0.0
casabiome_data % plantrate       (:,:) = 0.0
casabiome_data % rmplant         (:,:) = 0.0
casabiome_data % fracnpptoP      (:,:) = 0.0
casabiome_data % fraclignin      (:,:) = 0.0
casabiome_data % fraclabile      (:,:) = 0.0
casabiome_data % ratioNCplantmin (:,:) = 0.0
casabiome_data % ratioNCplantmax (:,:) = 0.0
casabiome_data % ratioNPplantmin (:,:) = 0.0
casabiome_data % ratioNPplantmax (:,:) = 0.0
casabiome_data % fracLigninplant (:,:) = 0.0
casabiome_data % ftransNPtoL     (:,:) = 0.0
casabiome_data % ftransPPtoL     (:,:) = 0.0
casabiome_data % litterrate      (:,:) = 0.0
casabiome_data % ratioPcplantmin (:,:) = 0.0
casabiome_data % ratioPcplantmax (:,:) = 0.0
casabiome_data % soilrate        (:,:) = 0.0

RETURN
END SUBROUTINE zero_casa_biome_data_type


SUBROUTINE assoc_casa_biome_type( casabiome, casabiome_data ) 

IMPLICIT NONE

TYPE (casa_biome_type),      INTENT(INOUT) :: casabiome
TYPE (casa_biome_data_type), INTENT(INOUT), TARGET :: casabiome_data

casabiome % ivt2            =>  casabiome_data % ivt2               
casabiome % xkleafcoldmax   =>  casabiome_data % xkleafcoldmax   
casabiome % xkleafcoldexp   =>  casabiome_data % xkleafcoldexp   
casabiome % xkleafdrymax    =>  casabiome_data % xkleafdrymax    
casabiome % xkleafdryexp    =>  casabiome_data % xkleafdryexp    
casabiome % glaimax         =>  casabiome_data % glaimax         
casabiome % glaimin         =>  casabiome_data % glaimin         
casabiome % sla             =>  casabiome_data % sla             
casabiome % ratiofrootleaf  =>  casabiome_data % ratiofrootleaf  
casabiome % kroot           =>  casabiome_data % kroot           
casabiome % krootlen        =>  casabiome_data % krootlen        
casabiome % rootdepth       =>  casabiome_data % rootdepth       
casabiome % kuptake         =>  casabiome_data % kuptake         
casabiome % kminN           =>  casabiome_data % kminN           
casabiome % kuplabP         =>  casabiome_data % kuplabP         
casabiome % kclabrate       =>  casabiome_data % kclabrate       
casabiome % xnpmax          =>  casabiome_data % xnpmax          
casabiome % q10soil         =>  casabiome_data % q10soil         
casabiome % xkoptlitter     =>  casabiome_data % xkoptlitter     
casabiome % xkoptsoil       =>  casabiome_data % xkoptsoil       
casabiome % xkplab          =>  casabiome_data % xkplab          
casabiome % xkpsorb         =>  casabiome_data % xkpsorb         
casabiome % xkpocc          =>  casabiome_data % xkpocc          
casabiome % prodptase       =>  casabiome_data % prodptase       
casabiome % costnpup        =>  casabiome_data % costnpup        
casabiome % maxfinelitter   =>  casabiome_data % maxfinelitter   
casabiome % maxcwd          =>  casabiome_data % maxcwd          
casabiome % nintercept      =>  casabiome_data % nintercept      
casabiome % nslope          =>  casabiome_data % nslope          
casabiome % plantrate       =>  casabiome_data % plantrate       
casabiome % rmplant         =>  casabiome_data % rmplant         
casabiome % fracnpptoP      =>  casabiome_data % fracnpptoP      
casabiome % fraclignin      =>  casabiome_data % fraclignin      
casabiome % fraclabile      =>  casabiome_data % fraclabile      
casabiome % ratioNCplantmin =>  casabiome_data % ratioNCplantmin 
casabiome % ratioNCplantmax =>  casabiome_data % ratioNCplantmax 
casabiome % ratioNPplantmin =>  casabiome_data % ratioNPplantmin 
casabiome % ratioNPplantmax =>  casabiome_data % ratioNPplantmax 
casabiome % fracLigninplant =>  casabiome_data % fracLigninplant 
casabiome % ftransNPtoL     =>  casabiome_data % ftransNPtoL     
casabiome % ftransPPtoL     =>  casabiome_data % ftransPPtoL     
casabiome % litterrate      =>  casabiome_data % litterrate      
casabiome % ratioPcplantmin =>  casabiome_data % ratioPcplantmin 
casabiome % ratioPcplantmax =>  casabiome_data % ratioPcplantmax 
casabiome % soilrate        =>  casabiome_data % soilrate

RETURN
END SUBROUTINE assoc_casa_biome_type


END MODULE casa_biome_type_mod
