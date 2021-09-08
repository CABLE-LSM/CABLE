MODULE allocate_cable_state_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: alloc_cable_state

!Functions allocating cable% types overloaded 
!jhan: Alloc routines could all initialise to NaN or zero for debugging?
PUBLIC :: alloc_cbm_var

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_CABLE_ARRAYS_MOD'

INTERFACE alloc_cbm_var
MODULE PROCEDURE alloc_balances_type,                                         &
  alloc_soil_snow_type,                                                       &
  alloc_canopy_type,                                                          &
  alloc_radiation_type,                                                       &
  alloc_roughness_type,                                                       &
  alloc_air_type,                                                             &
  alloc_met_type,                                                             &
  alloc_sum_flux_type,                                                        &
  alloc_bgc_pool_type !,                                                        &
  !alloc_climate_type
END INTERFACE

CONTAINS

SUBROUTINE alloc_cable_state ( mp, air_cbl, met_cbl, rad_cbl, rough_cbl,      &
                         canopy_cbl, ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl )

USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type

IMPLICIT NONE

INTEGER :: mp

TYPE(air_type)        :: air_cbl
TYPE(met_type)        :: met_cbl
TYPE(radiation_type)  :: rad_cbl
TYPE(roughness_type)  :: rough_cbl
TYPE(canopy_type)     :: canopy_cbl
TYPE(soil_snow_type)  :: ssnow_cbl
TYPE(bgc_pool_type)   :: bgc_cbl
TYPE(balances_type)   :: bal_cbl
TYPE(sum_flux_type)   :: sum_flux_cbl

LOGICAL, SAVE :: first_call =.TRUE.

IF ( first_call  ) THEN
  CALL alloc_cbm_var(air_cbl, mp)
  CALL alloc_cbm_var(met_cbl, mp)
  CALL alloc_cbm_var(rad_cbl, mp)
  CALL alloc_cbm_var(rough_cbl, mp)
  CALL alloc_cbm_var(canopy_cbl, mp)
  CALL alloc_cbm_var(ssnow_cbl, mp)
  CALL alloc_cbm_var(bgc_cbl, mp)
  CALL alloc_cbm_var(bal_cbl, mp)
  CALL alloc_cbm_var(sum_flux_cbl, mp)
END IF

first_call =.FALSE.
END SUBROUTINE alloc_cable_state

SUBROUTINE alloc_balances_type(var, mp)

USE cable_balances_type_mod,  ONLY: balances_type

IMPLICIT NONE

TYPE(balances_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE( var% drybal(mp) )
ALLOCATE( var% ebal(mp) )
ALLOCATE( var% ebal_tot(mp) )
ALLOCATE( var% ebaltr(mp) )
ALLOCATE( var% ebal_tottr(mp) )
ALLOCATE( var% ebal_cncheck(mp) )
ALLOCATE( var% ebal_tot_cncheck(mp) )
ALLOCATE( var% evap_tot(mp) )
ALLOCATE( var% osnowd0(mp) )
ALLOCATE( var% precip_tot(mp) )
ALLOCATE( var% rnoff_tot(mp) )
ALLOCATE( var% wbal(mp) )
ALLOCATE( var% wbal_tot(mp) )
ALLOCATE( var% wbtot0(mp) )
ALLOCATE( var% wetbal(mp) )
ALLOCATE( var% cansto0(mp) )
ALLOCATE( var% evapc_tot(mp) )
ALLOCATE( var% evaps_tot(mp) )
ALLOCATE( var% rnof1_tot(mp) )
ALLOCATE( var% rnof2_tot(mp) )
ALLOCATE( var% snowdc_tot(mp) )
ALLOCATE( var% wbal_tot1(mp) )
ALLOCATE( var% owbtot(mp) )
ALLOCATE( var% delwc_tot(mp) )
ALLOCATE( var% qasrf_tot(mp) )
ALLOCATE( var% qfsrf_tot(mp) )
ALLOCATE( var% qssrf_tot(mp) )

ALLOCATE( var% Radbal(mp) )
ALLOCATE( var% EbalSoil(mp) )
ALLOCATE( var% Ebalveg(mp) )
ALLOCATE( var% Radbalsum(mp) )

END SUBROUTINE alloc_balances_type

SUBROUTINE alloc_soil_snow_type(var, mp)

!USE jules_surface_types_mod,  ONLY: n_tiles => ntype !this was not available
!module in ESM1.5 *AND* convolutes where we are getting dims from
USE grid_constants_cbl_mod,   ONLY: trb !total # rad "bands"
USE grid_constants_cbl_mod,   ONLY: tsl !total # snow layers
!USE jules_soil_mod,           ONLY: ms => sm_levels   ! number of soil levels
!!Achtung! - how many places are these defined!
USE cable_types_mod,  ONLY: n_tiles
USE cable_types_mod,  ONLY : ms
USE cable_soil_snow_type_mod, ONLY: soil_snow_type

IMPLICIT NONE

TYPE(soil_snow_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE ( var % iantrct(mp) )
ALLOCATE ( var % pudsto(mp) )
ALLOCATE ( var % pudsmx(mp) )
ALLOCATE ( var % dtmlt(mp,3) )
ALLOCATE( var% albsoilsn(mp,trb) )
ALLOCATE( var% cls(mp) )
ALLOCATE( var% dfn_dtg(mp) )
ALLOCATE( var% dfh_dtg(mp) )
ALLOCATE( var% dfe_ddq(mp) )
ALLOCATE( var% ddq_dtg(mp) )
ALLOCATE( var% dfe_dtg(mp) )    !REV_CORR variable
ALLOCATE( var% evapsn(mp) )
ALLOCATE( var% fwtop(mp) )
ALLOCATE( var% fwtop1(mp) )
ALLOCATE( var% fwtop2(mp) )
ALLOCATE( var% fwtop3(mp) )
ALLOCATE( var% gammzz(mp,ms) )
ALLOCATE( var% isflag(mp) )
ALLOCATE( var% osnowd(mp) )
ALLOCATE( var% potev(mp) )
ALLOCATE( var% runoff(mp) )
ALLOCATE( var% rnof1(mp) )
ALLOCATE( var% rnof2(mp) )
ALLOCATE( var% rtsoil(mp) )
ALLOCATE( var% rtsoil_expl(mp) )
ALLOCATE( var% sconds(mp,tsl) )
ALLOCATE( var% sdepth(mp,tsl) )
ALLOCATE( var% smass(mp,tsl) )
ALLOCATE( var% snage(mp) )
ALLOCATE( var% snowd(mp) )
ALLOCATE( var% smelt(mp) )
ALLOCATE( var% ssdn(mp,tsl) )
ALLOCATE( var% ssdnn(mp) )
ALLOCATE( var% tgg(mp,ms) )
ALLOCATE( var% tggsn(mp,tsl) )
ALLOCATE( var% tss(mp) )
ALLOCATE( var% tss_p(mp) )
ALLOCATE( var% deltss(mp) )
ALLOCATE( var% owb1(mp) )
ALLOCATE( var% wb(mp,ms) )
ALLOCATE( var% wbice(mp,ms) )
ALLOCATE( var% wblf(mp,ms) )
ALLOCATE( var%wbtot(mp) )
ALLOCATE( var%wbtot1(mp) )
ALLOCATE( var%wbtot2(mp) )
ALLOCATE( var%wb_lake(mp) )
ALLOCATE( var%totwblake(mp) )
ALLOCATE( var%sinfil(mp) )
ALLOCATE( var%evapfbl(mp,ms) )
ALLOCATE( var%qstss(mp) )
ALLOCATE( var%wetfac(mp) )
ALLOCATE( var%owetfac(mp) )
ALLOCATE( var%t_snwlr(mp) )
ALLOCATE( var%wbfice(mp,ms) )
ALLOCATE( var%tggav(mp) )
ALLOCATE( var%otgg(mp) )
ALLOCATE( var%otss(mp) )
ALLOCATE( var%otss_0(mp) )
ALLOCATE( var%tprecip(mp) )
ALLOCATE( var%tevap(mp) )
ALLOCATE( var%trnoff(mp) )
ALLOCATE( var%totenbal(mp) )
ALLOCATE( var%totenbal2(mp) )
ALLOCATE( var%fland(mp) )
ALLOCATE( var%ifland(mp) )
ALLOCATE( var%tilefrac(mp,n_tiles) )
ALLOCATE( var%qasrf(mp) )
ALLOCATE( var%qfsrf(mp) )
ALLOCATE( var%qssrf(mp) )

!mrd561
!MD
!Aquifer variables
ALLOCATE( var%GWwb(mp) )
ALLOCATE( var%GWhk(mp) )
ALLOCATE( var%GWdhkdw(mp) )
ALLOCATE( var%GWdsmpdw(mp) )
ALLOCATE( var%wtd(mp) )
ALLOCATE( var%GWsmp(mp) )
ALLOCATE( var%GWwbeq(mp) )
ALLOCATE( var%GWzq(mp) )
ALLOCATE( var%qhz(mp) )
ALLOCATE( var%qhlev(mp,ms+1) )
ALLOCATE( var%satfrac(mp) )
ALLOCATE( var%Qrecharge(mp) )
ALLOCATE( var%rh_srf(mp) )
ALLOCATE( var%rtevap_unsat(mp) )
ALLOCATE( var%rtevap_sat(mp) )
ALLOCATE( var%rt_qh_sublayer(mp) )
!soil moisture variables
ALLOCATE( var%wbeq(mp,ms) )
ALLOCATE( var%zq(mp,ms) )
ALLOCATE( var%icefrac(mp,ms) )
ALLOCATE( var%fracice(mp,ms) )
ALLOCATE( var%hk(mp,ms) )
ALLOCATE( var%smp(mp,ms) )
ALLOCATE( var%dhkdw(mp,ms) )
ALLOCATE( var%dsmpdw(mp,ms) )
ALLOCATE( var%wbliq(mp,ms) )
ALLOCATE( var%wmliq(mp,ms) )
ALLOCATE( var%wmice(mp,ms) )
ALLOCATE( var%wmtot(mp,ms) )

! Allocate variables for SLI soil model:
!IF(cable_user%SOIL_STRUC=='sli') THEN
ALLOCATE ( var % s(mp,ms) )
ALLOCATE ( var % Tsoil(mp,ms) )
ALLOCATE ( var % sl(mp) )
ALLOCATE ( var % tl(mp) )
ALLOCATE ( var % h0(mp) )
ALLOCATE ( var % rex(mp,ms) )
ALLOCATE ( var % wflux(mp,0:ms) )
ALLOCATE ( var % delwcol(mp) )
ALLOCATE ( var % zdelta(mp) )
ALLOCATE ( var % kth(mp,ms) )
ALLOCATE ( var % Tsurface(mp) )
ALLOCATE ( var % lE(mp) )
ALLOCATE ( var % evap(mp) )
ALLOCATE ( var % ciso(mp,ms+1) )
ALLOCATE ( var % cisoL(mp) )
ALLOCATE ( var % rlitt(mp) )
ALLOCATE ( var % thetai(mp,ms) )
ALLOCATE ( var % snowliq(mp,3) )
ALLOCATE ( var % nsteps(mp) )
ALLOCATE ( var % nsnow(mp) )
ALLOCATE ( var % TsurfaceFR(mp) )
ALLOCATE ( var % Ta_daily(mp,100))
ALLOCATE ( var % Qadv_daily(mp) )
ALLOCATE ( var % G0_daily(mp) )
ALLOCATE ( var % Qevap_daily(mp) )
ALLOCATE ( var % Qprec_daily(mp) )
ALLOCATE ( var % Qprec_snow_daily(mp) )

!END IF

END SUBROUTINE alloc_soil_snow_type

SUBROUTINE alloc_canopy_type(var, mp)

!USE jules_soil_mod,           ONLY: ms => sm_levels  ! number of soil levels
!!Achtung! - how many places are these defined!
USE cable_types_mod,  ONLY : ms
USE grid_constants_cbl_mod,   ONLY: mf               ! # leaves (sunlit/shaded)
USE grid_constants_cbl_mod,   ONLY: niter            ! number of iterations for za/L
USE cable_canopy_type_mod,    ONLY: canopy_type

IMPLICIT NONE

TYPE(canopy_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

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

!! vh_js !! liiter resistances to heat and vapour transfer
ALLOCATE (var % kthLitt(mp))
ALLOCATE (var % DvLitt(mp))

END SUBROUTINE alloc_canopy_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_radiation_type(var, mp)

USE grid_constants_cbl_mod,   ONLY: trb
USE grid_constants_cbl_mod,   ONLY: mf !# leaves (sunlit/shaded)
USE grid_constants_cbl_mod,   ONLY: swb => n_rad_sw_bands ! number of spectral bands VISual/NIR
USE cable_radiation_type_mod, ONLY: radiation_type

IMPLICIT NONE

TYPE(radiation_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE( var% albedo(mp,trb) )
ALLOCATE( var% extkb(mp) )
ALLOCATE( var% extkd2(mp) )
ALLOCATE( var% extkd(mp) )
ALLOCATE( var% flws(mp) )
ALLOCATE( var% fvlai(mp,mf) )
ALLOCATE( var% latitude(mp) )
ALLOCATE( var% lwabv(mp) )
ALLOCATE( var% qcan(mp,mf,trb) )
ALLOCATE( var% qssabs(mp) )
ALLOCATE( var% rhocdf(mp,trb) )
ALLOCATE( var% rniso(mp,mf) )
ALLOCATE( var% scalex(mp,mf) )
ALLOCATE( var% transd(mp) )
ALLOCATE( var% trad(mp) )
ALLOCATE( var% otrad(mp) )
ALLOCATE( var% reffdf(mp,trb) )
ALLOCATE( var% reffbm(mp,trb) )
ALLOCATE( var% extkbm(mp,trb) )
ALLOCATE( var% extkdm(mp,trb) )
ALLOCATE( var% cexpkbm(mp,swb) )
ALLOCATE( var% cexpkdm(mp,swb) )
ALLOCATE( var% fbeam(mp,trb) )
ALLOCATE( var% rhocbm(mp,trb) )
ALLOCATE( var% transb(mp) )
ALLOCATE( var% albedo_T(mp) )
ALLOCATE( var% gradis(mp,mf) )
ALLOCATE( var% longitude(mp) )
ALLOCATE( var% workp1(mp) )
ALLOCATE( var% workp2(mp) )
ALLOCATE( var% workp3(mp) )

END SUBROUTINE alloc_radiation_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_roughness_type(var, mp)

USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type

IMPLICIT NONE

TYPE(roughness_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE ( var % coexp(mp) )
ALLOCATE ( var % disp(mp) )
ALLOCATE ( var % hruff(mp) )
ALLOCATE ( var % hruff_grmx(mp) )
ALLOCATE ( var % rt0us(mp) )
ALLOCATE ( var % rt1usa(mp) )
ALLOCATE ( var % rt1usb(mp) )
ALLOCATE ( var % rt1(mp) )
ALLOCATE ( var % term2(mp) )
ALLOCATE ( var % term3(mp) )
ALLOCATE ( var % term5(mp) )
ALLOCATE ( var % term6(mp) )
ALLOCATE ( var % term6a(mp) )
ALLOCATE ( var % usuh(mp) )
ALLOCATE ( var % za_uv(mp) )
ALLOCATE ( var % za_tq(mp) )
ALLOCATE ( var % z0m(mp) )
ALLOCATE ( var % zref_uv(mp) )
ALLOCATE ( var % zref_tq(mp) )
ALLOCATE ( var % zruffs(mp) )
ALLOCATE ( var % z0soilsn(mp) )
ALLOCATE ( var % z0soil(mp) )



END SUBROUTINE alloc_roughness_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_air_type(var, mp)

USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type

IMPLICIT NONE

TYPE(air_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE ( var % rho(mp) )
ALLOCATE ( var % volm(mp) )
ALLOCATE ( var % rlam(mp) )
ALLOCATE ( var % qsat(mp) )
ALLOCATE ( var % epsi(mp) )
ALLOCATE ( var % visc(mp) )
ALLOCATE ( var % psyc(mp) )
ALLOCATE ( var % dsatdk(mp) )
ALLOCATE ( var % cmolar(mp) )

END SUBROUTINE alloc_air_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_met_type(var, mp)

USE grid_constants_cbl_mod,   ONLY: swb => n_rad_sw_bands ! number of spectral bands VISual/NIR
USE cable_met_type_mod,       ONLY: met_type

IMPLICIT NONE

TYPE(met_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE ( var % ca(mp) )
ALLOCATE ( var % year(mp) )
ALLOCATE ( var % moy(mp) )
ALLOCATE ( var % doy(mp) )
ALLOCATE ( var % hod(mp) )
ALLOCATE ( var % fsd(mp,swb) )
ALLOCATE ( var % ofsd(mp) )
ALLOCATE ( var % fld(mp) )
ALLOCATE ( var % precip(mp) )
ALLOCATE ( var % precip_sn(mp) )
ALLOCATE ( var % tk(mp) )
ALLOCATE ( var % tvair(mp) )
ALLOCATE ( var % tvrad(mp) )
ALLOCATE ( var % pmb(mp) )
ALLOCATE ( var % ua(mp) )
ALLOCATE ( var % qv(mp) )
ALLOCATE ( var % qvair(mp) )
ALLOCATE ( var % da(mp) )
ALLOCATE ( var % dva(mp) )
ALLOCATE ( var % coszen(mp) )
ALLOCATE ( var % Ndep(mp) )
ALLOCATE ( var % Pdep(mp) )

END SUBROUTINE alloc_met_type

! ------------------------------------------------------------------------------

!H!  SUBROUTINE alloc_climate_type(var, mp)
!H!
!H!USE cable_air_type_mod,       ONLY : air_type
!H!USE cable_met_type_mod,       ONLY : met_type
!H!USE cable_radiation_type_mod, ONLY : radiation_type
!H!USE cable_roughness_type_mod, ONLY : roughness_type
!H!USE cable_canopy_type_mod,    ONLY : canopy_type
!H!USE cable_soil_snow_type_mod, ONLY : soil_snow_type
!H!USE cable_bgc_pool_type_mod,  ONLY : bgc_pool_type
!H!USE cable_balances_type_mod,  ONLY : balances_type
!H!USE cable_sum_flux_type_mod,  ONLY : sum_flux_type
!H!
!H!implicit none
!H!
!H!    TYPE(climate_type), INTENT(inout) :: var
!H!    INTEGER, INTENT(in) :: mp
!H!    INTEGER :: ny, nd
!H!    ny = var%nyear_average
!H!    nd = var%nday_average
!H!    PRINT*, 'ny', ny
!H!    PRINT*, 'nd', nd
!H!    !   ALLOCATE ( var %  nyears )
!H!    !   ALLOCATE ( var %  doy )
!H!    ALLOCATE ( var %  dtemp(mp) )
!H!    ALLOCATE ( var %  dmoist(mp) )
!H!    ALLOCATE ( var % mtemp(mp) )
!H!    ALLOCATE ( var % qtemp(mp) )
!H!    ALLOCATE ( var % mmoist(mp) )
!H!    ALLOCATE ( var % mtemp_min(mp) )
!H!    ALLOCATE ( var %  mtemp_max20(mp) )
!H!    ALLOCATE ( var % mtemp_min20(mp) )
!H!    ALLOCATE ( var %  mtemp_max(mp) )
!H!    ALLOCATE ( var %  qtemp_max(mp) )
!H!    ALLOCATE ( var %  qtemp_max_last_year(mp) )
!H!    ALLOCATE ( var % atemp_mean(mp) )
!H!    ALLOCATE ( var % AGDD5(mp) )
!H!    ALLOCATE ( var % GDD5(mp) )
!H!    ALLOCATE ( var % AGDD0(mp) )
!H!    ALLOCATE ( var % GDD0(mp) )
!H!    ALLOCATE ( var % chilldays(mp) )
!H!    ALLOCATE ( var % iveg(mp) )
!H!    ALLOCATE ( var % biome(mp) )
!H!    ALLOCATE ( var % alpha_PT(mp) )
!H!    ALLOCATE ( var % alpha_PT20(mp) )
!H!    ALLOCATE ( var % evap_PT(mp) )
!H!    ALLOCATE ( var % aevap(mp) )
!H!
!H!    ALLOCATE ( var % mtemp_min_20(mp,ny) )
!H!    ALLOCATE ( var %     mtemp_max_20(mp,ny) )
!H!    ALLOCATE ( var %     dtemp_31(mp,nd) )
!H!    ALLOCATE ( var %     dmoist_31(mp,nd) )
!H!    ALLOCATE ( var %     dtemp_91(mp,91) )
!H!    ALLOCATE ( var %     alpha_PT_20(mp,ny) )
!H!
!H!  END SUBROUTINE alloc_climate_type

  ! ------------------------------------------------------------------------------

SUBROUTINE alloc_sum_flux_type(var, mp)

USE cable_sum_flux_type_mod,  ONLY: sum_flux_type

IMPLICIT NONE

TYPE(sum_flux_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE ( var % sumpn(mp) )
ALLOCATE ( var % sumrp(mp) )
ALLOCATE ( var % sumrpw(mp) )
ALLOCATE ( var % sumrpr(mp) )
ALLOCATE ( var % sumrs(mp) )
ALLOCATE ( var % sumrd(mp) )
ALLOCATE ( var % dsumpn(mp) )
ALLOCATE ( var % dsumrp(mp) )
ALLOCATE ( var % dsumrs(mp) )
ALLOCATE ( var % dsumrd(mp) )
ALLOCATE ( var % sumxrp(mp) )
ALLOCATE ( var % sumxrs(mp) )

END SUBROUTINE alloc_sum_flux_type

SUBROUTINE alloc_bgc_pool_type(var, mp)

USE grid_constants_cbl_mod,   ONLY: ncp             ! # vegetation carbon stores
USE grid_constants_cbl_mod,   ONLY: ncs             ! # soil carbon stores

USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type

IMPLICIT NONE

TYPE(bgc_pool_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

ALLOCATE ( var % cplant(mp,ncp) )
ALLOCATE ( var % csoil(mp,ncs) )
ALLOCATE ( var % ratecp(ncp) )
ALLOCATE ( var % ratecs(ncs) )
!H!ALLOCATE ( var % ratecp(mp,ncp) )
!H!ALLOCATE ( var % ratecs(mp,ncs) )

END SUBROUTINE alloc_bgc_pool_type




END MODULE allocate_cable_state_mod
