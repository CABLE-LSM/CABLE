MODULE cable_canopy_type_mod

USE cable_types_mod, ONLY: r_2

IMPLICIT NONE

  ! Canopy/vegetation variables:
TYPE canopy_type

  REAL, DIMENSION(:), POINTER ::                                              &
       cansto,  & ! canopy water storage (mm)
       cduv,    & ! drag coefficient for momentum
       delwc,   & ! change in canopy water store (mm/dels)
       dewmm,   & ! dewfall (mm)
       fe,      & ! total latent heat (W/m2)
       fh,      & ! total sensible heat (W/m2)
       fpn,     & ! plant photosynthesis (g C m-2 s-1)
       frp,     & ! plant respiration (g C m-2 s-1)
       frpw,    & ! plant respiration (woody component) (g C m-2 s-1)
       frpr,    & ! plant respiration (root component) (g C m-2 s-1)
       frs,     & ! soil respiration (g C m-2 s-1)
       fnee,    & ! net carbon flux (g C m-2 s-1)
       frday,   & ! daytime leaf resp
       fnv,     & ! net rad. avail. to canopy (W/m2)
       fev,     & ! latent hf from canopy (W/m2)
       epot,    & ! total potential evaporation
       fnpp,    & ! npp flux
       fevw_pot,& ! potential lat heat from canopy
       gswx_T,  & ! ! stom cond for water
       cdtq,    & ! drag coefficient for momentum
       wetfac_cs,&!
       fevw,    & ! lat heat fl wet canopy (W/m2)
       fhvw,    & ! sens heatfl from wet canopy (W/m2)
       oldcansto,&! canopy water storage (mm)
       fhv,     & ! sens heatfl from canopy (W/m2)
       fns,     & ! net rad avail to soil (W/m2)
       fhs,     & ! sensible heat flux from soil
       fhs_cor,                                                               &
       ga,      & ! ground heat flux (W/m2) ???
       ghflux,  & ! ground heat flux (W/m2) ???
       precis,  & ! throughfall to soil, after snow (mm)
       qscrn,   & ! specific humudity at screen height (g/g)
       rnet,    & ! net radiation absorbed by surface (W/m2)
       rniso,    & !isothermal net radiation absorbed by surface (W/m2)
       segg,    & ! latent heatfl from soil mm
       sghflux, & ! ground heat flux (W/m2) ???
       through, & ! canopy throughfall (mm)
       through_sn, & ! canopy snow throughfall (equal to precip_sn) (mm)
       spill,   & ! can.storage excess after dewfall (mm)
       tscrn,   & ! air temperature at screen height (oC)
       wcint,   & ! canopy rainfall interception (mm)
       tv,      & ! vegetation temp (K)
       us,      & ! friction velocity
       uscrn,   & ! wind speed at screen height (m/s)
       vlaiw,   & ! lai adj for snow depth for calc of resistances
       rghlai,  & ! lai adj for snow depth for calc of resistances
       fwet       ! fraction of canopy wet

  !INH - new REV_CORR coupling variables
  REAL, DIMENSION(:), POINTER ::                                              &
       fns_cor, & ! correction to net rad avail to soil (W/m2)
       ga_cor  ! correction to ground heat flux (W/m2)

  REAL, DIMENSION(:,:), POINTER ::                                            &
       evapfbl,                                                               &
       gswx,    & ! stom cond for water
       zetar, &   ! stability parameter (ref height)
                             ! vh_js !
       zetash      ! stability parameter (shear height)

  REAL(r_2), DIMENSION(:), POINTER ::                                         &
       fess,    & ! latent heatfl from soil (W/m2)
       fesp,    & ! latent heatfl from soil (W/m2)
       dgdtg,   & ! derivative of gflux wrt soil temp
       fes,     & ! latent heatfl from soil (W/m2)
       fes_cor, & ! latent heatfl from soil (W/m2)
       fevc,     &  ! dry canopy transpiration (W/m2)
       ofes     ! latent heatfl from soil (W/m2)

  !SSEB - new variables limits on correction terms - for future use
  !REAL(r_2), DIMENSION(:), POINTER ::                                     &
  !  fescor_upp,& ! upper limit on the correction term fes_cor (W/m2)
  !  fescor_low   ! lower limit on the correction term fes_cor (W/m2)

  REAL(r_2), DIMENSION(:), POINTER ::                                         &
       sublayer_dz

  ! Additional variables:
  REAL(r_2), DIMENSION(:,:),   POINTER :: gw     ! dry canopy conductance (ms-1) edit vh 6/7/09
  REAL(r_2), DIMENSION(:,:,:), POINTER :: ancj   ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
  REAL(r_2), DIMENSION(:,:),   POINTER :: tlfy   ! sunlit and shaded leaf temperatures
  REAL(r_2), DIMENSION(:,:),   POINTER :: ecy    ! sunlit and shaded leaf transpiration (dry canopy)
  REAL(r_2), DIMENSION(:,:),   POINTER :: ecx    ! sunlit and shaded leaf latent heat flux
  REAL(r_2), DIMENSION(:,:,:), POINTER :: ci     ! intra-cellular CO2 vh 6/7/09
  REAL(r_2), DIMENSION(:),     POINTER :: fwsoil !

  ! vh_js ! !litter thermal conductivity (Wm-2K-1) and vapour diffusivity (m2s-1)
  REAL(r_2), DIMENSION(:), POINTER :: kthLitt, DvLitt


END TYPE canopy_type

!Instantiation:
TYPE(canopy_type) :: canopy_cbl

CONTAINS

SUBROUTINE alloc_canopy_type(var, mp)

USE grid_constants_cbl_mod,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_cbl_mod,   ONLY: niter            ! number of iterations for za/L

IMPLICIT NONE

TYPE(canopy_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp
!INTEGER, INTENT(IN) :: ms
INTEGER, PARAMETER  :: ms = 6

ALLOCATE ( var % fess(mp) )
ALLOCATE ( var % fesp(mp) )
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
ALLOCATE ( var % evapfbl(mp,ms) )
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
!ALLOCATE( var% fescor_upp(mp) )  !SSEB variable
!ALLOCATE( var% fescor_low(mp) )  !SSEB variable
ALLOCATE( var% gswx(mp,mf) )
ALLOCATE( var% oldcansto(mp) )
ALLOCATE( var% zetar(mp,niter) )
ALLOCATE( var% zetash(mp,niter) )
ALLOCATE ( var % fwsoil(mp) )
ALLOCATE ( var % ofes(mp) )

ALLOCATE( var%sublayer_dz(mp) )

ALLOCATE ( var % gw(mp,mf) )     ! dry canopy conductance (ms-1) edit vh 6/7/09
ALLOCATE ( var % ancj(mp,mf,3) ) ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
ALLOCATE ( var % tlfy(mp,mf) )   ! sunlit and shaded leaf temperatures
ALLOCATE ( var % ecy(mp,mf) )    ! sunlit and shaded leaf transpiration (dry canopy)
ALLOCATE ( var % ecx(mp,mf) )    ! sunlit and shaded leaf latent heat flux
ALLOCATE ( var % ci(mp,mf,3) )   ! intra-cellular CO2 vh 6/7/09
ALLOCATE ( var % fwsoil (mp) )

! vh_js ! liiter resistances to heat and vapour transfer
ALLOCATE (var % kthLitt(mp))
ALLOCATE (var % DvLitt(mp))

END SUBROUTINE alloc_canopy_type


END MODULE cable_canopy_type_mod
