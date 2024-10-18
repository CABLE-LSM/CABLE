MODULE prognostic_bank_mod_cbl

TYPE ProgBank
   
  REAL, DIMENSION(:,:), ALLOCATABLE ::                                       &
    tsoil, smcl, sthf,                                                       & 
    snow_depth,                                                              &
    snow_mass, snow_tmp, snow_rho ,                                          &
    sli_s, sli_tsoil, sli_sconds, sli_snowliq,                               &
    cplant, csoil !carbon variables
  
  REAL, DIMENSION(:), ALLOCATABLE ::                                         &
    snow_rho1l, snow_age, snow_tile,                                         &
    ocanopy,                                                                 &
    fes_cor,fhs_cor, osnowd,owetfac,otss,gwwb,tss0,                        &
    puddle, rtsoil, wblake, gwaq,                                            &
    sli_h0, sli_tsurf
  
  INTEGER, DIMENSION(:), ALLOCATABLE ::                                      &
    snow_flg3l, sli_nsnow
  
END TYPE ProgBank

CONTAINS

SUBROUTINE cable_store_prognostics( pb, mp, nsl, nsCs,nvCs, ssnow, canopy, bgc )

USE cable_common_module,   ONLY: cable_user
USE cable_def_types_mod,   ONLY: soil_snow_type, bgc_pool_type 
USE cable_canopy_type_mod, ONLY: canopy_type

IMPLICIT NONE
INTEGER, INTENT(IN) :: mp        ! # active tiles 
INTEGER, INTENT(IN) :: nsl       ! # soil layers
INTEGER, INTENT(IN) :: nsCs      ! # soil carbon stores
INTEGER, INTENT(IN) :: nvCs      ! # vegetation carbon stores

TYPE (ProgBank),      INTENT(INOUT):: pb
TYPE(soil_snow_type), INTENT(INOUT) :: ssnow 
TYPE(canopy_type),    INTENT(INOUT) :: canopy
TYPE(bgc_pool_type),  INTENT(INOUT) :: bgc

IF (.NOT. ALLOCATED(pb%tsoil) ) then

  ALLOCATE( pb%tsoil(mp,nsl) )
  ALLOCATE( pb%smcl(mp,nsl) )
  ALLOCATE( pb%sthf(mp,nsl) )
  ALLOCATE( pb%snow_depth(mp,3) )
  ALLOCATE( pb%snow_mass(mp,3) )
  ALLOCATE( pb%snow_tmp(mp,3) )
  ALLOCATE( pb%snow_rho(mp,3) )
  ALLOCATE( pb%snow_rho1l(mp) )
  ALLOCATE( pb%snow_age(mp) )
  ALLOCATE( pb%snow_flg3l(mp) )
  ALLOCATE( pb%snow_tile(mp) )
  ALLOCATE( pb%puddle(mp) )
  ALLOCATE( pb%owetfac(mp) )
  ALLOCATE( pb%rtsoil(mp) )
  ALLOCATE( pb%wblake(mp) )
  ALLOCATE(pb%tss0(mp) )
  ALLOCATE( pb%ocanopy(mp) )
  ALLOCATE( pb%fes_cor(mp) )
  !carbon variables - may not be needed unless CASA
  ALLOCATE( pb%cplant(mp,nvcs) )
  ALLOCATE( pb%csoil(mp,nscs) )
  !GW model variables no need to test always has a value
  !so do not introduce issues with restarting a GW run 
  ALLOCATE (pb%GWaq(mp) )
  
  !SLI variables - Jhan please check the second dimension
  IF (cable_user%soil_struc=='sli') THEN
    ALLOCATE(pb%sli_nsnow(mp) )
    ALLOCATE(pb%sli_s(mp,nsl) )
    ALLOCATE(pb%sli_tsoil(mp,nsl) )
    ALLOCATE(pb%sli_sconds(mp,3) )
    ALLOCATE(pb%sli_h0(mp) )
    ALLOCATE(pb%sli_tsurf(mp) )
    ALLOCATE(pb%sli_snowliq(mp,3) )
  ENDIF
    
ENDIF !.NOT. allocated   

pb%tsoil      = ssnow%tgg
pb%smcl       = ssnow%wb
pb%sthf       = ssnow%wbice
pb%snow_depth = ssnow%sdepth
pb%snow_mass  = ssnow%smass
pb%snow_tmp   = ssnow%tggsn
pb%snow_rho   = ssnow%ssdn
pb%snow_rho1l = ssnow%ssdnn
pb%snow_age   = ssnow%snage
pb%snow_flg3l = ssnow%isflag
pb%snow_tile  = ssnow%snowd
pb%puddle     = ssnow%pudsto
pb%rtsoil     = ssnow%rtsoil !?needed in restart?
pb%owetfac    = ssnow%owetfac
pb%wblake     = ssnow%wb_lake 
pb%tss0       = ssnow%tss
pb%ocanopy    = canopy%oldcansto
pb%fes_cor    = canopy%fes_cor
pb%cplant     = bgc%cplant  ! may not be needed unless CASA
pb%csoil      = bgc%csoil   ! may not be needed unless CASA
pb%gwaq       = ssnow%gwwb  !GW model ?issues with restarting? 
    
!SLI variables
IF (cable_user%soil_struc=='sli') THEN
  pb%sli_nsnow   = ssnow%nsnow
  pb%sli_s       = ssnow%s
  pb%tsoil       = ssnow%tsoil
  pb%sli_sconds  = ssnow%sconds
  pb%sli_h0      = ssnow%h0
  pb%sli_Tsurf   = ssnow%tsurface
  pb%sli_snowliq = ssnow%snowliq
ENDIF

RETURN
END SUBROUTINE cable_store_prognostics

SUBROUTINE cable_reinstate_prognostics( pb, ssnow, canopy, bgc )

USE cable_def_types_mod,   ONLY: soil_snow_type, bgc_pool_type 
USE cable_canopy_type_mod, ONLY: canopy_type
USE cable_common_module,   ONLY: cable_user

IMPLICIT NONE
TYPE (ProgBank),      INTENT(INOUT) :: pb
TYPE(soil_snow_type), INTENT(INOUT) :: ssnow 
TYPE(canopy_type),    INTENT(INOUT) :: canopy
TYPE(bgc_pool_type),  INTENT(INOUT) :: bgc

ssnow%tgg        = pb%tsoil
ssnow%wb         = pb%smcl
ssnow%wbice      = pb%sthf
ssnow%sdepth     = pb%snow_depth
ssnow%smass      = pb%snow_mass
ssnow%tggsn      = pb%snow_tmp
ssnow%ssdn       = pb%snow_rho
ssnow%ssdnn      = pb%snow_rho1l
ssnow%snage      = pb%snow_age
ssnow%isflag     = pb%snow_flg3l
ssnow%snowd      = pb%snow_tile
ssnow%pudsto     = pb%puddle
ssnow%rtsoil     = pb%rtsoil  !?needed
ssnow%owetfac    = pb%owetfac
ssnow%wb_lake    = pb%wblake  
ssnow%tss        = pb%tss0
canopy%oldcansto = pb%ocanopy
canopy%fes_cor   = pb%fes_cor
bgc%cplant       = pb%cplant ! may not be needed unless CASA
bgc%csoil        = pb%csoil   ! may not be needed unless CASA
ssnow%GWwb       = pb%GWaq! GW model variables

!SLI variables
IF (cable_user%soil_struc=='sli') THEN
  ssnow%nsnow    = pb%sli_nsnow
  ssnow%S        = pb%sli_s
  ssnow%Tsoil    = pb%sli_tsoil
  ssnow%sconds   = pb%sli_sconds
  ssnow%h0       = pb%sli_h0
  ssnow%Tsurface = pb%sli_tsurf
  ssnow%snowliq  = pb%sli_snowliq
ENDIF

RETURN
END SUBROUTINE cable_reinstate_prognostics
 
END MODULE prognostic_bank_mod_cbl
