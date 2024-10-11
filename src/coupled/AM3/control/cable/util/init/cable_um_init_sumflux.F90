MODULE cable_um_init_sumflux_mod
   
IMPLICIT NONE
PUBLIC init_sumflux_zero

CONTAINS

SUBROUTINE init_sumflux_zero( sum_flux ) 

USE cable_def_types_mod, ONLY: sum_flux_type
                               
IMPLICIT NONE

TYPE(sum_flux_type), INTENT(OUT) :: sum_flux

sum_flux%sumpn  = 0.0 
sum_flux%sumrp  = 0.0 
sum_flux%sumrpw = 0.0
sum_flux%sumrpr = 0.0 
sum_flux%sumrs  = 0.0 
sum_flux%sumrd  = 0.0
sum_flux%dsumpn = 0.0 
sum_flux%dsumrp = 0.0 
sum_flux%dsumrs = 0.0
sum_flux%dsumrd = 0.0 
sum_flux%sumxrp = 0.0  
sum_flux%sumxrs = 0.0
RETURN
END SUBROUTINE init_sumflux_zero 

END MODULE cable_um_init_sumflux_mod




