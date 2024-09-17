MODULE cable_um_init_bgc_mod
   
IMPLICIT NONE
PUBLIC init_bgc_vars

CONTAINS
 
SUBROUTINE init_bgc_vars( pars, bgc, veg ) 
USE cable_def_types_mod, ONLY : ncs, ncp 
USE cable_def_types_mod, ONLY: bgc_pool_type, veg_parameter_type
USE params_io_mod_cbl,   ONLY: params_io_data_type
IMPLICIT NONE
TYPE(params_io_data_type), INTENT(IN)  :: pars
TYPE(veg_parameter_type),  INTENT(IN)  :: veg        ! vegetation parameters
TYPE(bgc_pool_type),       INTENT(OUT) :: bgc

INTEGER :: k

! note that ratecp and ratecs are the same for all veg at the moment. (BP)
DO k=1,ncp
  bgc%cplant(:,k) = pars%vegin_cplant(k,veg%iveg)
  bgc%ratecp(k)   = pars%vegin_ratecp(k,1)
ENDDO
DO k=1,ncs
  bgc%csoil(:,k) = pars%vegin_csoil(k,veg%iveg)
  bgc%ratecs(k)  = pars%vegin_ratecs(k,1)
ENDDO

END SUBROUTINE init_bgc_vars

END MODULE cable_um_init_bgc_mod




