MODULE init_cable_runtime_opts_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: init_cable_runtime_options

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='def_cable_grid_mod'

CONTAINS

SUBROUTINE init_cable_runtime_options()

USE cable_runtime_opts_mod ,ONLY: cable_user
USE cable_runtime_opts_mod ,ONLY: satuparam
USE cable_runtime_opts_mod ,ONLY: wiltparam

IMPLICIT NONE

!call def_cable_runtime_opts()
!USE declared in here cable_runtime_opts_mod ,ONLY : cable_user
!already declared with a size no need to allocate
!call alloc_cable_state 
!grab these from the read values as with the params
!in the first instance HW here 

!H!cable_user%diag_soil_resp='ON'
!H-test!cable_user%fwsoil_switch='standard'
!H-test!cable_user%gs_switch='leuning'
cable_user%fwsoil_switch='standard'
cable_user%gs_switch='medlyn'
!H!match Loobos 
!H!cable_user%l_rev_corr = .TRUE.
!H!cable_user%l_revised_coupling = .TRUE.
cable_user%soil_thermal_fix = .TRUE.
cable_user%ssnow_potev='HDM'
!jh:these only matter if calling hydraulic_resistribution
!satuparam=0.8
!wiltparam=0.5

END SUBROUTINE init_cable_runtime_options

END MODULE init_cable_runtime_opts_mod
