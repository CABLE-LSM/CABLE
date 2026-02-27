MODULE cable_canopy_type_mod

USE cable_other_constants_mod, ONLY: r_2

IMPLICIT NONE

PUBLIC :: canopy_type
PUBLIC :: canopy_data_type
PUBLIC :: alloc_canopy_type
PUBLIC :: dealloc_canopy_type
PUBLIC :: assoc_canopy_type
PUBLIC :: nullify_canopy_cbl

! Canopy/vegetation variables:
TYPE canopy_data_type

  REAL, ALLOCATABLE, PUBLIC :: cansto(:)            ! canopy water storage (mm)
  REAL, ALLOCATABLE, PUBLIC :: cduv(:)              ! drag coefficient for momentum
  REAL, ALLOCATABLE, PUBLIC :: delwc(:)             ! change in canopy water store (mm/dels)
  REAL, ALLOCATABLE, PUBLIC :: dewmm(:)             ! dewfall (mm)
  REAL, ALLOCATABLE, PUBLIC :: fe(:)                ! total latent heat (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: fh(:)                ! total sensible heat (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: fpn(:)               ! plant photosynthesis (g C m-2 s-1)
  REAL, ALLOCATABLE, PUBLIC :: frp(:)               ! plant respiration (g C m-2 s-1)
  REAL, ALLOCATABLE, PUBLIC :: frpw(:)              ! plant respiration (woody component) (g C m-2 s-1)
  REAL, ALLOCATABLE, PUBLIC :: frpr(:)              ! plant respiration (root component) (g C m-2 s-1)
  REAL, ALLOCATABLE, PUBLIC :: frs(:)               ! soil respiration (g C m-2 s-1)
  REAL, ALLOCATABLE, PUBLIC :: fnee(:)              ! net carbon flux (g C m-2 s-1)
  REAL, ALLOCATABLE, PUBLIC :: frday(:)             ! daytime leaf resp
  REAL, ALLOCATABLE, PUBLIC :: fnv(:)               ! net rad. avail. to canopy (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: fev(:)               ! latent hf from canopy (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: epot(:)              ! total potential evaporation
  REAL, ALLOCATABLE, PUBLIC :: fnpp(:)              ! npp flux
  REAL, ALLOCATABLE, PUBLIC :: fevw_pot(:)          ! potential lat heat from canopy
  REAL, ALLOCATABLE, PUBLIC :: gswx_T(:)            ! ! stom cond for water
  REAL, ALLOCATABLE, PUBLIC :: cdtq(:)              ! drag coefficient for momentum
  REAL, ALLOCATABLE, PUBLIC :: wetfac_cs(:)         !
  REAL, ALLOCATABLE, PUBLIC :: fevw(:)              ! lat heat fl wet canopy (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: fhvw(:)              ! sens heatfl from wet canopy (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: oldcansto(:)         ! canopy water storage (mm)
  REAL, ALLOCATABLE, PUBLIC :: fhv(:)               ! sens heatfl from canopy (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: fns(:)               ! net rad avail to soil (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: fhs(:)               ! sensible heat flux from soil
  REAL, ALLOCATABLE, PUBLIC :: fhs_cor(:)           ! 
  REAL, ALLOCATABLE, PUBLIC :: ga(:)                ! ground heat flux (W/m2) ???
  REAL, ALLOCATABLE, PUBLIC :: ghflux(:)            ! ground heat flux (W/m2) ???
  REAL, ALLOCATABLE, PUBLIC :: precis(:)            ! throughfall to soil, after snow (mm)
  REAL, ALLOCATABLE, PUBLIC :: qscrn(:)             ! specific humudity at screen height (g/g)
  REAL, ALLOCATABLE, PUBLIC :: rnet(:)              ! net radiation absorbed by surface (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: rniso(:)             !isothermal net radiation absorbed by surface (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: segg(:)              ! latent heatfl from soil mm
  REAL, ALLOCATABLE, PUBLIC :: sghflux(:)           ! ground heat flux (W/m2) ???
  REAL, ALLOCATABLE, PUBLIC :: through(:)           ! canopy throughfall (mm)
  REAL, ALLOCATABLE, PUBLIC :: through_sn(:)        ! canopy snow throughfall (equal to precip_sn) (mm)
  REAL, ALLOCATABLE, PUBLIC :: spill(:)             ! can.storage excess after dewfall (mm)
  REAL, ALLOCATABLE, PUBLIC :: tscrn(:)             ! air temperature at screen height (oC)
  REAL, ALLOCATABLE, PUBLIC :: wcint(:)             ! canopy rainfall interception (mm)
  REAL, ALLOCATABLE, PUBLIC :: tv(:)                ! vegetation temp (K)
  REAL, ALLOCATABLE, PUBLIC :: us(:)                ! friction velocity
  REAL, ALLOCATABLE, PUBLIC :: uscrn(:)             ! wind speed at screen height (m/s)
  REAL, ALLOCATABLE, PUBLIC :: vlaiw(:)             ! lai adj for snow depth for calc of resistances
  REAL, ALLOCATABLE, PUBLIC :: rghlai(:)            ! lai adj for snow depth for calc of resistances
  REAL, ALLOCATABLE, PUBLIC :: fwet(:)              ! fraction of canopy wet
  REAL, ALLOCATABLE, PUBLIC :: fns_cor(:)           ! correction to net rad avail to soil (W/m2)
  REAL, ALLOCATABLE, PUBLIC :: ga_cor(:)            ! correction to ground heat flux (W/m2)    
  REAL, ALLOCATABLE, PUBLIC :: gswx(:,:)            ! stom cond for water
  REAL, ALLOCATABLE, PUBLIC :: zetar(:,:)           ! stability parameter (ref height)
  REAL, ALLOCATABLE, PUBLIC :: zetash(:,:)          ! stability parameter (shear height)
  REAL(r_2), ALLOCATABLE, PUBLIC :: fess(:)         ! latent heatfl from soil (w/m2)
  REAL(r_2), ALLOCATABLE, PUBLIC :: fesp(:)         ! latent heatfl from soil (w/m2)
  REAL(r_2), ALLOCATABLE, PUBLIC :: dgdtg(:)        ! derivative of gflux wrt soil temp
  REAL(r_2), ALLOCATABLE, PUBLIC :: fes(:)          ! latent heatfl from soil (w/m2)
  REAL(r_2), ALLOCATABLE, PUBLIC :: fes_cor(:)      ! latent heatfl from soil (w/m2)
  REAL(r_2), ALLOCATABLE, PUBLIC :: fevc(:)         ! dry canopy transpiration (w/m2)
  REAL(r_2), ALLOCATABLE, PUBLIC :: ofes(:)         ! latent heatfl from soil (w/m2)
  REAL(r_2), ALLOCATABLE, PUBLIC :: sublayer_dz(:)  !
  REAL(r_2), ALLOCATABLE, PUBLIC :: gw(:,:)         ! dry canopy conductance (ms-1) edit vh 6/7/09
  REAL(r_2), ALLOCATABLE, PUBLIC :: ancj(:,:,:)     ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
  REAL(r_2), ALLOCATABLE, PUBLIC :: tlfy(:,:)       ! sunlit and shaded leaf temperatures
  REAL(r_2), ALLOCATABLE, PUBLIC :: ecy(:,:)        ! sunlit and shaded leaf transpiration (dry canopy)
  REAL(r_2), ALLOCATABLE, PUBLIC :: ecx(:,:)        ! sunlit and shaded leaf latent heat flux
  REAL(r_2), ALLOCATABLE, PUBLIC :: ci(:,:,:)       ! intra-cellular CO2 vh 6/7/09
  REAL(r_2), ALLOCATABLE, PUBLIC :: fwsoil(:)       !
  ! vh_js ! !litter thermal conductivity (Wm-2K-1) and vapour diffusivity (m2s-1)
  REAL(r_2), ALLOCATABLE, PUBLIC :: kthLitt(:)      !
  REAL(r_2), ALLOCATABLE, PUBLIC :: DvLitt(:)       !
  !SSEB - new variables limits on correction terms - for future use
  !REAL(r_2), DIMENSION(:), POINTER ::                                     &
  !  fescor_upp,& ! upper limit on the correction term fes_cor (W/m2)
  !  fescor_low   ! lower limit on the correction term fes_cor (W/m2)

END TYPE canopy_data_type

! Canopy/vegetation variables:
TYPE canopy_type

  REAL, POINTER, PUBLIC :: cansto(:)            ! canopy water storage (mm)
  REAL, POINTER, PUBLIC :: cduv(:)              ! drag coefficient for momentum
  REAL, POINTER, PUBLIC :: delwc(:)             ! change in canopy water store (mm/dels)
  REAL, POINTER, PUBLIC :: dewmm(:)             ! dewfall (mm)
  REAL, POINTER, PUBLIC :: fe(:)                ! total latent heat (W/m2)
  REAL, POINTER, PUBLIC :: fh(:)                ! total sensible heat (W/m2)
  REAL, POINTER, PUBLIC :: fpn(:)               ! plant photosynthesis (g C m-2 s-1)
  REAL, POINTER, PUBLIC :: frp(:)               ! plant respiration (g C m-2 s-1)
  REAL, POINTER, PUBLIC :: frpw(:)              ! plant respiration (woody component) (g C m-2 s-1)
  REAL, POINTER, PUBLIC :: frpr(:)              ! plant respiration (root component) (g C m-2 s-1)
  REAL, POINTER, PUBLIC :: frs(:)               ! soil respiration (g C m-2 s-1)
  REAL, POINTER, PUBLIC :: fnee(:)              ! net carbon flux (g C m-2 s-1)
  REAL, POINTER, PUBLIC :: frday(:)             ! daytime leaf resp
  REAL, POINTER, PUBLIC :: fnv(:)               ! net rad. avail. to canopy (W/m2)
  REAL, POINTER, PUBLIC :: fev(:)               ! latent hf from canopy (W/m2)
  REAL, POINTER, PUBLIC :: epot(:)              ! total potential evaporation
  REAL, POINTER, PUBLIC :: fnpp(:)              ! npp flux
  REAL, POINTER, PUBLIC :: fevw_pot(:)          ! potential lat heat from canopy
  REAL, POINTER, PUBLIC :: gswx_T(:)            ! ! stom cond for water
  REAL, POINTER, PUBLIC :: cdtq(:)              ! drag coefficient for momentum
  REAL, POINTER, PUBLIC :: wetfac_cs(:)         !
  REAL, POINTER, PUBLIC :: fevw(:)              ! lat heat fl wet canopy (W/m2)
  REAL, POINTER, PUBLIC :: fhvw(:)              ! sens heatfl from wet canopy (W/m2)
  REAL, POINTER, PUBLIC :: oldcansto(:)         ! canopy water storage (mm)
  REAL, POINTER, PUBLIC :: fhv(:)               ! sens heatfl from canopy (W/m2)
  REAL, POINTER, PUBLIC :: fns(:)               ! net rad avail to soil (W/m2)
  REAL, POINTER, PUBLIC :: fhs(:)               ! sensible heat flux from soil
  REAL, POINTER, PUBLIC :: fhs_cor(:)           ! 
  REAL, POINTER, PUBLIC :: ga(:)                ! ground heat flux (W/m2) ???
  REAL, POINTER, PUBLIC :: ghflux(:)            ! ground heat flux (W/m2) ???
  REAL, POINTER, PUBLIC :: precis(:)            ! throughfall to soil, after snow (mm)
  REAL, POINTER, PUBLIC :: qscrn(:)             ! specific humudity at screen height (g/g)
  REAL, POINTER, PUBLIC :: rnet(:)              ! net radiation absorbed by surface (W/m2)
  REAL, POINTER, PUBLIC :: rniso(:)             !isothermal net radiation absorbed by surface (W/m2)
  REAL, POINTER, PUBLIC :: segg(:)              ! latent heatfl from soil mm
  REAL, POINTER, PUBLIC :: sghflux(:)           ! ground heat flux (W/m2) ???
  REAL, POINTER, PUBLIC :: through(:)           ! canopy throughfall (mm)
  REAL, POINTER, PUBLIC :: through_sn(:)        ! canopy snow throughfall (equal to precip_sn) (mm)
  REAL, POINTER, PUBLIC :: spill(:)             ! can.storage excess after dewfall (mm)
  REAL, POINTER, PUBLIC :: tscrn(:)             ! air temperature at screen height (oC)
  REAL, POINTER, PUBLIC :: wcint(:)             ! canopy rainfall interception (mm)
  REAL, POINTER, PUBLIC :: tv(:)                ! vegetation temp (K)
  REAL, POINTER, PUBLIC :: us(:)                ! friction velocity
  REAL, POINTER, PUBLIC :: uscrn(:)             ! wind speed at screen height (m/s)
  REAL, POINTER, PUBLIC :: vlaiw(:)             ! lai adj for snow depth for calc of resistances
  REAL, POINTER, PUBLIC :: rghlai(:)            ! lai adj for snow depth for calc of resistances
  REAL, POINTER, PUBLIC :: fwet(:)              ! fraction of canopy wet
  REAL, POINTER, PUBLIC :: fns_cor(:)           ! correction to net rad avail to soil (W/m2)
  REAL, POINTER, PUBLIC :: ga_cor(:)            ! correction to ground heat flux (W/m2)    
  REAL, POINTER, PUBLIC :: gswx(:,:)            ! stom cond for water
  REAL, POINTER, PUBLIC :: zetar(:,:)           ! stability parameter (ref height)
  REAL, POINTER, PUBLIC :: zetash(:,:)          ! stability parameter (shear height)
  REAL(r_2), POINTER, PUBLIC :: fess(:)         ! latent heatfl from soil (w/m2)
  REAL(r_2), POINTER, PUBLIC :: fesp(:)         ! latent heatfl from soil (w/m2)
  REAL(r_2), POINTER, PUBLIC :: dgdtg(:)        ! derivative of gflux wrt soil temp
  REAL(r_2), POINTER, PUBLIC :: fes(:)          ! latent heatfl from soil (w/m2)
  REAL(r_2), POINTER, PUBLIC :: fes_cor(:)      ! latent heatfl from soil (w/m2)
  REAL(r_2), POINTER, PUBLIC :: fevc(:)         ! dry canopy transpiration (w/m2)
  REAL(r_2), POINTER, PUBLIC :: ofes(:)         ! latent heatfl from soil (w/m2)
  REAL(r_2), POINTER, PUBLIC :: sublayer_dz(:)  !
  REAL(r_2), POINTER, PUBLIC :: gw(:,:)         ! dry canopy conductance (ms-1) edit vh 6/7/09
  REAL(r_2), POINTER, PUBLIC :: ancj(:,:,:)     ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
  REAL(r_2), POINTER, PUBLIC :: tlfy(:,:)       ! sunlit and shaded leaf temperatures
  REAL(r_2), POINTER, PUBLIC :: ecy(:,:)        ! sunlit and shaded leaf transpiration (dry canopy)
  REAL(r_2), POINTER, PUBLIC :: ecx(:,:)        ! sunlit and shaded leaf latent heat flux
  REAL(r_2), POINTER, PUBLIC :: ci(:,:,:)       ! intra-cellular CO2 vh 6/7/09
  REAL(r_2), POINTER, PUBLIC :: fwsoil(:)       !
  REAL(r_2), POINTER, PUBLIC :: kthLitt(:)      !!litter thermal conductivity (Wm-2K-1) and vapour diffusivity (m2s-1)
  REAL(r_2), POINTER, PUBLIC :: DvLitt(:)       !
  !SSEB - new variables limits on correction terms - for future use
  !REAL(r_2), DIMENSION(:), POINTER ::                                     &
  !  fescor_upp,& ! upper limit on the correction term fes_cor (W/m2)
  !  fescor_low   ! lower limit on the correction term fes_cor (W/m2)

END TYPE canopy_type


CONTAINS

SUBROUTINE alloc_canopy_type(var, mp)

USE grid_constants_mod_cbl,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_mod_cbl,   ONLY: nsl              ! # soil layers                
USE grid_constants_mod_cbl,   ONLY: niter            ! number of iterations for za/L

IMPLICIT NONE

TYPE(canopy_data_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE( var% fess(mp) )
ALLOCATE( var% fesp(mp) )
ALLOCATE( var% cansto(mp) )
ALLOCATE( var% cduv(mp) )
ALLOCATE( var% delwc(mp) )
ALLOCATE( var% dewmm(mp) )
ALLOCATE( var% dgdtg(mp) )
ALLOCATE( var% fe(mp) )
ALLOCATE( var% fh(mp) )
ALLOCATE( var% fpn(mp) )
ALLOCATE( var% frp(mp) )
ALLOCATE( var% frpw(mp) )
ALLOCATE( var% frpr(mp) )
ALLOCATE( var% frs(mp) )
ALLOCATE( var% fnee(mp) )
ALLOCATE( var% frday(mp) )
ALLOCATE( var% fnv(mp) )
ALLOCATE( var% fev(mp) )
ALLOCATE( var% fevc(mp) )
ALLOCATE( var% fhv(mp) )
ALLOCATE( var% fns(mp) )
ALLOCATE( var% fhs(mp) )
ALLOCATE( var% fhs_cor(mp) )
ALLOCATE( var% ga(mp) )
ALLOCATE( var% ghflux(mp) )
ALLOCATE( var% precis(mp) )
ALLOCATE( var% qscrn(mp) )
ALLOCATE( var% rnet(mp) )
ALLOCATE( var% rniso(mp) )
ALLOCATE( var% segg(mp) )
ALLOCATE( var% sghflux(mp) )
ALLOCATE( var% through(mp) )
ALLOCATE( var% through_sn(mp) )
ALLOCATE( var% spill(mp) )
ALLOCATE( var% tscrn(mp) )
ALLOCATE( var% wcint(mp) )
ALLOCATE( var% tv(mp) )
ALLOCATE( var% us(mp) )
ALLOCATE( var% uscrn(mp) )
ALLOCATE( var% rghlai(mp) )
ALLOCATE( var% vlaiw(mp) )
ALLOCATE( var% fwet(mp) )
ALLOCATE( var% fns_cor(mp) )    !REV_CORR variable
ALLOCATE( var% ga_cor(mp) )     !REV_CORR variable
ALLOCATE( var% epot(mp) )
ALLOCATE( var% fnpp(mp) )
ALLOCATE( var% fevw_pot(mp) )
ALLOCATE( var% gswx_T(mp) )
ALLOCATE( var% cdtq(mp) )
ALLOCATE( var% wetfac_cs(mp) )
ALLOCATE( var% fevw(mp) )
ALLOCATE( var% fhvw(mp) )
ALLOCATE( var% fes(mp) )
ALLOCATE( var% fes_cor(mp) )
ALLOCATE( var% gswx(mp,mf) )
ALLOCATE( var% oldcansto(mp) )
ALLOCATE( var% zetar(mp,niter) )
ALLOCATE( var% zetash(mp,niter) )
ALLOCATE( var % fwsoil(mp) )
ALLOCATE( var % ofes(mp) )
ALLOCATE( var%sublayer_dz(mp) )
ALLOCATE( var % gw(mp,mf) )     ! dry canopy conductance (ms-1) edit vh 6/7/09
ALLOCATE( var % ancj(mp,mf,3) ) ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
ALLOCATE( var % tlfy(mp,mf) )   ! sunlit and shaded leaf temperatures
ALLOCATE( var % ecy(mp,mf) )    ! sunlit and shaded leaf transpiration (dry canopy)
ALLOCATE( var % ecx(mp,mf) )    ! sunlit and shaded leaf latent heat flux
ALLOCATE( var % ci(mp,mf,3) )   ! intra-cellular CO2 vh 6/7/09
ALLOCATE(var % kthLitt(mp))
ALLOCATE(var % DvLitt(mp))
!ALLOCATE( var% fescor_upp(mp) )  !SSEB variable
!ALLOCATE( var% fescor_low(mp) )  !SSEB variable

var % cansto(:)      = 0.0      
var % cduv(:)        = 0.0      
var % delwc(:)       = 0.0      
var % dewmm(:)       = 0.0      
var % fe(:)          = 0.0      
var % fh(:)          = 0.0      
var % fpn(:)         = 0.0      
var % frp(:)         = 0.0      
var % frpw(:)        = 0.0      
var % frpr(:)        = 0.0      
var % frs(:)         = 0.0      
var % fnee(:)        = 0.0      
var % frday(:)       = 0.0      
var % fnv(:)         = 0.0      
var % fev(:)         = 0.0      
var % epot(:)        = 0.0      
var % fnpp(:)        = 0.0      
var % fevw_pot(:)    = 0.0      
var % gswx_T(:)      = 0.0      
var % cdtq(:)        = 0.0      
var % wetfac_cs(:)   = 0.0      
var % fevw(:)        = 0.0      
var % fhvw(:)        = 0.0      
var % oldcansto(:)   = 0.0      
var % fhv(:)         = 0.0      
var % fns(:)         = 0.0      
var % fhs(:)         = 0.0      
var % fhs_cor(:)     = 0.0      
var % ga(:)          = 0.0      
var % ghflux(:)      = 0.0      
var % precis(:)      = 0.0      
var % qscrn(:)       = 0.0      
var % rnet(:)        = 0.0      
var % rniso(:)       = 0.0      
var % segg(:)        = 0.0      
var % sghflux(:)     = 0.0      
var % through(:)     = 0.0      
var % through_sn(:)  = 0.0      
var % spill(:)       = 0.0      
var % tscrn(:)       = 0.0      
var % wcint(:)       = 0.0      
var % tv(:)          = 0.0      
var % us(:)          = 0.0      
var % uscrn(:)       = 0.0      
var % vlaiw(:)       = 0.0      
var % rghlai(:)      = 0.0      
var % fwet(:)        = 0.0      
var % fns_cor(:)     = 0.0      
var % ga_cor(:)      = 0.0      
var % gswx(:,:)      = 0.0      
var % zetar(:,:)     = 0.0      
var % zetash(:,:)    = 0.0      
var % fess(:)        = 0.0_r_2 
var % fesp(:)        = 0.0_r_2 
var % dgdtg(:)       = 0.0_r_2 
var % fes(:)         = 0.0_r_2 
var % fes_cor(:)     = 0.0_r_2 
var % fevc(:)        = 0.0_r_2 
var % ofes(:)        = 0.0_r_2 
var % sublayer_dz(:) = 0.0_r_2 
var % gw(:,:)        = 0.0_r_2 
var % ancj(:,:,:)    = 0.0_r_2 
var % tlfy(:,:)      = 0.0_r_2 
var % ecy(:,:)       = 0.0_r_2 
var % ecx(:,:)       = 0.0_r_2 
var % ci(:,:,:)      = 0.0_r_2 
var % fwsoil(:)      = 0.0_r_2 
var % kthLitt(:)     = 0.0_r_2 
var % DvLitt(:)      = 0.0_r_2 


RETURN
END SUBROUTINE alloc_canopy_type

SUBROUTINE dealloc_canopy_type(var)

   TYPE(canopy_type), INTENT(inout) :: var

   DEALLOCATE ( var % fess )
   DEALLOCATE ( var % fesp )
   DEALLOCATE( var% cansto )
   DEALLOCATE( var% cduv )
   DEALLOCATE( var% delwc )
   DEALLOCATE( var% dewmm )
   DEALLOCATE( var% dgdtg )
   DEALLOCATE( var% fe )
   DEALLOCATE( var% fh )
   DEALLOCATE( var% fpn )
   DEALLOCATE( var% frp )
   DEALLOCATE( var% frpw )
   DEALLOCATE( var% frpr )
   DEALLOCATE( var% frs )
   DEALLOCATE( var% fnee )
   DEALLOCATE( var% frday )
   DEALLOCATE( var% fnv )
   DEALLOCATE( var% fev )
   DEALLOCATE( var% fevc )
   DEALLOCATE( var% fhv )
   DEALLOCATE( var% fns )
   DEALLOCATE( var% fhs )
   DEALLOCATE( var% fhs_cor )
   DEALLOCATE( var% ga )
   DEALLOCATE( var% ghflux )
   DEALLOCATE( var% precis )
   DEALLOCATE( var% qscrn )
   DEALLOCATE( var% rnet )
   DEALLOCATE( var% rniso )
   DEALLOCATE( var% segg )
   DEALLOCATE( var% sghflux )
   DEALLOCATE( var% through )
   DEALLOCATE( var% through_sn )
   DEALLOCATE( var% spill )
   DEALLOCATE( var% tscrn )
   DEALLOCATE( var% wcint )
   DEALLOCATE( var% tv )
   DEALLOCATE( var% us )
   DEALLOCATE( var% uscrn )
   DEALLOCATE( var% rghlai )
   DEALLOCATE( var% vlaiw )
   DEALLOCATE( var% fwet )
   DEALLOCATE( var% fns_cor )   !REV_CORR variable
   DEALLOCATE( var% ga_cor )    !REV_CORR variable
   DEALLOCATE( var% epot )
   DEALLOCATE( var% fnpp )
   DEALLOCATE( var% fevw_pot )
   DEALLOCATE( var% gswx_T )
   DEALLOCATE( var% cdtq )
   DEALLOCATE( var% wetfac_cs )
   DEALLOCATE( var% fevw )
   DEALLOCATE( var% fhvw )
   DEALLOCATE( var% fes )
   DEALLOCATE( var% fes_cor )
   DEALLOCATE( var% gswx )
   DEALLOCATE( var% oldcansto )
   DEALLOCATE( var% zetar )
   DEALLOCATE( var% zetash )
   DEALLOCATE ( var % fwsoil )
   DEALLOCATE ( var % ofes )
   DEALLOCATE( var% sublayer_dz )
   DEALLOCATE (var % kthLitt)
   DEALLOCATE (var % DvLitt)
   !DEALLOCATE( var% fescor_upp ) !SSEB variable
   !DEALLOCATE( var% fescor_low ) !SSEB variable

END SUBROUTINE dealloc_canopy_type

SUBROUTINE assoc_canopy_type(canopy, canopy_data )

! Description:
!   Associate the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(canopy_type),      INTENT(IN OUT)         :: canopy
TYPE(canopy_data_type), INTENT(IN OUT), TARGET :: canopy_data

CHARACTER(LEN=*), PARAMETER :: RoutineName=''

!End of header

CALL nullify_canopy_cbl(canopy)

canopy% cansto      => canopy_data% cansto
canopy% cduv        => canopy_data% cduv
canopy% delwc       => canopy_data% delwc
canopy% dewmm       => canopy_data% dewmm
canopy% fe          => canopy_data% fe
canopy% fh          => canopy_data% fh
canopy% fpn         => canopy_data% fpn        
canopy% frp         => canopy_data% frp        
canopy% frpw        => canopy_data% frpw
canopy% frpr        => canopy_data% frpr
canopy% frs         => canopy_data% frs
canopy% fnee        => canopy_data% fnee
canopy% frday       => canopy_data% frday
canopy% fnv         => canopy_data% fnv
canopy% fev         => canopy_data% fev
canopy% epot        => canopy_data% epot
canopy% fnpp        => canopy_data% fnpp
canopy% fevw_pot    => canopy_data% fevw_pot
canopy% gswx_T      => canopy_data% gswx_T
canopy% cdtq        => canopy_data% cdtq
canopy% wetfac_cs   => canopy_data% wetfac_cs
canopy% fevw        => canopy_data% fevw
canopy% fhvw        => canopy_data% fhvw
canopy% oldcansto   => canopy_data% oldcansto
canopy% fhv         => canopy_data% fhv        
canopy% fns         => canopy_data% fns        
canopy% fhs         => canopy_data% fhs        
canopy% fhs_cor     => canopy_data% fhs_cor
canopy% ga          => canopy_data% ga
canopy% ghflux      => canopy_data% ghflux
canopy% precis      => canopy_data% precis
canopy% qscrn       => canopy_data% qscrn
canopy% rnet        => canopy_data% rnet
canopy% rniso       => canopy_data% rniso
canopy% segg        => canopy_data% segg
canopy% sghflux     => canopy_data% sghflux
canopy% through     => canopy_data% through
canopy% through_sn  => canopy_data% through_sn
canopy% spill       => canopy_data% spill
canopy% tscrn       => canopy_data% tscrn
canopy% wcint       => canopy_data% wcint
canopy% tv          => canopy_data% tv
canopy% us          => canopy_data% us
canopy% uscrn       => canopy_data% uscrn
canopy% vlaiw       => canopy_data% vlaiw
canopy% rghlai      => canopy_data% rghlai
canopy% fwet        => canopy_data% fwet
canopy% fns_cor     => canopy_data% fns_cor
canopy% ga_cor      => canopy_data% ga_cor
canopy% gswx        => canopy_data% gswx
canopy% zetar       => canopy_data% zetar
canopy% zetash      => canopy_data% zetash
canopy% fess        => canopy_data% fess
canopy% fesp        => canopy_data% fesp
canopy% dgdtg       => canopy_data% dgdtg
canopy% fes         => canopy_data% fes
canopy% fes_cor     => canopy_data% fes_cor
canopy% fevc        => canopy_data% fevc
canopy% ofes        => canopy_data% ofes
canopy% sublayer_dz => canopy_data% sublayer_dz
canopy% gw          => canopy_data% gw
canopy% ancj        => canopy_data% ancj
canopy% tlfy        => canopy_data% tlfy
canopy% ecy         => canopy_data% ecy
canopy% ecx         => canopy_data% ecx
canopy% ci          => canopy_data% ci
canopy% fwsoil      => canopy_data% fwsoil
canopy% kthLitt     => canopy_data% kthLitt
canopy% DvLitt      => canopy_data% DvLitt
 
RETURN
END SUBROUTINE assoc_canopy_type

!===============================================================================
SUBROUTINE nullify_canopy_cbl( var )

! Description:
!   Nullify the CABLE work pointers in the derived type structure

IMPLICIT NONE

!Arguments
TYPE(canopy_type), INTENT(IN OUT) :: var 

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_WORK_VARS_CBL'

!End of header

   NULLIFY( var % fess )
   NULLIFY( var % fesp )
   NULLIFY( var% cansto )
   NULLIFY( var% cduv )
   NULLIFY( var% delwc )
   NULLIFY( var% dewmm )
   NULLIFY( var% dgdtg )
   NULLIFY( var% fe )
   NULLIFY( var% fh )
   NULLIFY( var% fpn )
   NULLIFY( var% frp )
   NULLIFY( var% frpw )
   NULLIFY( var% frpr )
   NULLIFY( var% frs )
   NULLIFY( var% fnee )
   NULLIFY( var% frday )
   NULLIFY( var% fnv )
   NULLIFY( var% fev )
   NULLIFY( var% fevc )
   NULLIFY( var% fhv )
   NULLIFY( var% fns )
   NULLIFY( var% fhs )
   NULLIFY( var% fhs_cor )
   NULLIFY( var% ga )
   NULLIFY( var% ghflux )
   NULLIFY( var% precis )
   NULLIFY( var% qscrn )
   NULLIFY( var% rnet )
   NULLIFY( var% rniso )
   NULLIFY( var% segg )
   NULLIFY( var% sghflux )
   NULLIFY( var% through )
   NULLIFY( var% through_sn )
   NULLIFY( var% spill )
   NULLIFY( var% tscrn )
   NULLIFY( var% wcint )
   NULLIFY( var% tv )
   NULLIFY( var% us )
   NULLIFY( var% uscrn )
   NULLIFY( var% rghlai )
   NULLIFY( var% vlaiw )
   NULLIFY( var% fwet )
   NULLIFY( var% fns_cor )   !REV_CORR variable
   NULLIFY( var% ga_cor )    !REV_CORR variable
   NULLIFY( var% epot )
   NULLIFY( var% fnpp )
   NULLIFY( var% fevw_pot )
   NULLIFY( var% gswx_T )
   NULLIFY( var% cdtq )
   NULLIFY( var% wetfac_cs )
   NULLIFY( var% fevw )
   NULLIFY( var% fhvw )
   NULLIFY( var% fes )
   NULLIFY( var% fes_cor )
   NULLIFY( var% gswx )
   NULLIFY( var% oldcansto )
   NULLIFY( var% zetar )
   NULLIFY( var% zetash )
   NULLIFY( var % fwsoil )
   NULLIFY( var % ofes )
   NULLIFY( var% sublayer_dz )
   NULLIFY( var % kthLitt)! liiter resistances to heat and vapour transfer
   NULLIFY( var % DvLitt)
   !DEALLOCATE( var% fescor_upp ) !SSEB variable
   !DEALLOCATE( var% fescor_low ) !SSEB variable

RETURN

END SUBROUTINE nullify_canopy_cbl

END MODULE cable_canopy_type_mod
