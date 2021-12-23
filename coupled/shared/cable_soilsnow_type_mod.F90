MODULE cable_soil_snow_type_mod

USE cable_types_mod, ONLY: r_2

IMPLICIT NONE

  ! Soil and snow variables:
TYPE soil_snow_type

  INTEGER, DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow

  REAL, DIMENSION(:), POINTER ::                                              &
       iantrct, & ! pointer to Antarctic land points
       pudsto,  & ! puddle storage
       pudsmx,  & ! puddle storage
       cls,     & ! factor for latent heat
       dfn_dtg, & ! d(canopy%fns)/d(ssnow%tgg)
       dfh_dtg, & ! d(canopy%fhs)/d(ssnow%tgg)
       dfe_ddq, & ! d(canopy%fes)/d(dq)        - REV_CORR: no longer necessary
       ddq_dtg, & ! d(dq)/d(ssnow%tgg)         - REV_CORR: no longer necessary
       dfe_dtg, & ! d(canopy%fes)/d(ssnow%tgg) - REV_CORR: covers above vars
       evapsn,  & ! snow evaporation
       fwtop,   & ! water flux to the soil
       fwtop1,  & ! water flux to the soil
       fwtop2,  & ! water flux to the soil
       fwtop3,  & ! water flux to the soil
       osnowd,  & ! snow depth from previous time step
       potev,   & ! potential evapotranspiration
       runoff,  & ! total runoff (mm/dels)
       rnof1,   & ! surface runoff (mm/dels)
       rnof2,   & ! deep drainage (mm/dels)
       rtsoil,  & ! turbulent resistance for soil
       rtsoil_expl,  & ! turbulent resistance for soil (explicit path)
       wbtot1,  & ! total soil water (mm)
       wbtot2,  & ! total soil water (mm)
       wb_lake,                                                               &
       totwblake, & !daily integrated wb_lake: used in ACCESS
       sinfil,                                                                &
       qstss,                                                                 &
       wetfac,  & ! surface wetness fact. at current time step
       owetfac, & ! surface wetness fact. at previous time step
       t_snwlr, & ! top snow layer depth in 3 layer snowpack
       tggav,   & ! mean soil temperature in K
       otgg,    & ! soil temperature in K
       otss,    & ! surface temperature (weighted soil, snow)
       otss_0,  & ! surface temperature (weighted soil, snow)
       tprecip,                                                               &
       tevap,                                                                 &
       trnoff,                                                                &
       totenbal,&!
       totenbal2,                                                             &
       fland,   & ! factor for latent heat
       ifland,  & ! integer soil type
       qasrf,   & ! heat advected to the snow by precip.
       qfsrf,   & ! energy of snowpack phase changes
       qssrf,   & ! sublimation
       snage,   & ! snow age
       snowd,   & ! snow depth (liquid water)
       smelt,   & ! snow melt
       ssdnn,   & ! average snow density
       tss,     & ! surface temperature (weighted soil, snow)
       tss_p,   & ! surface temperature (weighted soil, snow)
       deltss,  & ! surface temperature (weighted soil, snow)
       owb1       ! surface temperature (weighted soil, snow)

  REAL, DIMENSION(:,:), POINTER ::                                            &
       sconds,     & !
       sdepth,     & ! snow depth
       smass,      & ! snow mass
       ssdn,       & ! snow densities
       tgg,        & ! soil temperature in K
       tggsn,      & ! snow temperature in K
       dtmlt,      & ! water flux to the soil
       albsoilsn,  & ! soil + snow reflectance
       evapfbl,    & !
       tilefrac      ! factor for latent heat


  REAL(r_2), DIMENSION(:), POINTER ::                                         &
       wbtot   ! total soil water (mm)

  REAL(r_2), DIMENSION(:,:), POINTER ::                                       &
       gammzz,  & ! heat capacity for each soil layer
       wb,      & ! volumetric soil moisture (solid+liq)
       wbice,   & ! soil ice
       wblf,    & !
       wbfice     !

  !mrd561
  !MD variables for the revised soil moisture + GW scheme
  REAL(r_2), DIMENSION(:), POINTER   ::                                       &
       GWwb,    &  ! water content in aquifer [mm3/mm3]
       GWhk,    &  ! aquifer hydraulic conductivity  [mm/s]
       GWdhkdw, &  ! aquifer d(hk) over d(water content) [(mm/s)/(mm3/mm3)]
       GWdsmpdw,&  ! aquifer d(smp) / dw   [(mm)/(mm3/mm3)]
       wtd,     &  ! water table depth   [mm]
       GWsmp,   &  ! aquifer soil matric potential [mm]
       GWwbeq,  &  ! equilibrium aquifer water content [mm3/mm3]
       GWzq,    &  ! equilibrium aquifer smp   [mm]
       qhz, & ! horizontal hydraulic conductivity in 1D gw model for soil layers  [mm/s]
       satfrac,                                                               &
       Qrecharge,                                                             &
       rh_srf,                                                                &
       rtevap_sat,                                                            &
       rtevap_unsat,                                                          &
       rt_qh_sublayer

  REAL(r_2), DIMENSION(:,:), POINTER  ::                                      &
       wbeq,    &    ! equilibrium water content [mm3/mm3]
       zq,      &    ! equilibrium smp       [mm]
       icefrac, &    ! ice fraction  [none]  -> ice mass / total mass
       fracice, &    ! alternate ice fraction  [none] - parameterized
       hk,      &    ! hydraulic conductivity for soil layers [mm/s]
       smp,     &    ! soil matric potential for soil layers         [mm]
       dhkdw, & ! d(hydraulic conductivity ) d(water) for soil layers [(mm/s)/(mm3/mm3)]
       dsmpdw,  &    ! d(smp)/ d(water) for soil layers   [(mm)/(mm3/mm3)]
       wbliq,   &    ! volumetric liquid water content  [mm3/mm3]
       wmliq,   &    !water mass [mm] liq
       wmice,   &    !water mass [mm] ice
       wmtot,   &    !water mass [mm] liq+ice ->total
       qhlev
  ! Additional SLI variables:
  REAL(r_2), DIMENSION(:,:), POINTER :: s         ! moisture content relative to sat value    (edit vh 23/01/08)
  REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil         !     Tsoil (deg C)
  REAL(r_2), DIMENSION(:),   POINTER :: sl        ! litter moisture content relative to sat value (edit vh 23/01/08)
  REAL(r_2), DIMENSION(:),   POINTER :: tl        ! litter temperature in K     (edit vh 23/01/08)
  REAL(r_2), DIMENSION(:),   POINTER :: h0        ! pond height in m            (edit vh 23/01/08)
  REAL(r_2), DIMENSION(:,:), POINTER :: rex       ! root extraction from each layer (mm/dels)
  REAL(r_2), DIMENSION(:,:), POINTER :: wflux     ! water flux at layer boundaries (mm s-1)
  REAL(r_2), DIMENSION(:),   POINTER :: delwcol   ! change in water column (mm / dels)
  REAL(r_2), DIMENSION(:),   POINTER :: zdelta    ! water table depth           (edit vh 23/06/08)
  REAL(r_2), DIMENSION(:,:), POINTER :: kth       ! thermal conductivity           (edit vh 29/07/08)
  REAL(r_2), DIMENSION(:),   POINTER :: Tsurface  !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
  REAL(r_2), DIMENSION(:),   POINTER :: lE        ! soil latent heat flux
  REAL(r_2), DIMENSION(:),   POINTER :: evap      ! soil evaporation (mm / dels)
  REAL(r_2), DIMENSION(:,:), POINTER :: ciso      ! concentration of minor isotopologue in soil water (kg m-3 water)
  REAL(r_2), DIMENSION(:),   POINTER :: cisoL     ! concentration of minor isotopologue in litter water (kg m-3 water)
  REAL(r_2), DIMENSION(:),   POINTER :: rlitt     ! resistance to heat/moisture transfer through litter (m-1 s)
  REAL(r_2), DIMENSION(:,:), POINTER :: thetai    ! volumetric ice content (MC)
  REAL(r_2), DIMENSION(:,:), POINTER :: snowliq   ! liquid snow content (mm H2O)
  REAL(r_2), DIMENSION(:),   POINTER :: nsteps    ! number of iterations at each timestep
  REAL(r_2), DIMENSION(:),   POINTER :: TsurfaceFR  !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
  REAL(r_2), DIMENSION(:,:), POINTER :: Ta_daily        ! air temp averaged over last 24h
  INTEGER, DIMENSION(:),     POINTER :: nsnow ! number of layers in snow-pack (0-nsnow_max)
  REAL(r_2), DIMENSION(:),   POINTER :: Qadv_daily  ! advective heat flux into surface , daily average (W m-2)
  REAL(r_2), DIMENSION(:),   POINTER :: G0_daily  ! conductive heat flux into surface , daily average (W m-2)
  REAL(r_2), DIMENSION(:),   POINTER :: Qevap_daily ! evaporative flux at surface, daily average (m s-1)
  REAL(r_2), DIMENSION(:),   POINTER :: Qprec_daily ! liquid precip, daily average (m s-1)
  REAL(r_2), DIMENSION(:),   POINTER :: Qprec_snow_daily ! solid precip, daily average (m s-1)



END TYPE soil_snow_type

!Instantiation:
TYPE(soil_snow_type) :: ssnow_cbl

CONTAINS

SUBROUTINE alloc_soil_snow_type(var, mp)

USE grid_constants_cbl_mod,   ONLY: ntiles  ! # land-cover types
USE grid_constants_cbl_mod,   ONLY: nsl     ! # soil layers
USE grid_constants_cbl_mod,   ONLY: nsnl    ! # snow layers
USE grid_constants_cbl_mod,   ONLY: nrb     ! # rad "bands"

IMPLICIT NONE

TYPE(soil_snow_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp
INTEGER, PARAMETER :: SLI_nsnow_max = 2 !see nsnow_max in SLI_numbers() 

ALLOCATE( var % iantrct(mp) )
ALLOCATE( var % pudsto(mp) )
ALLOCATE( var % pudsmx(mp) )
ALLOCATE( var % dtmlt(mp,nsnl) )
ALLOCATE( var% albsoilsn(mp,nrb) )
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
ALLOCATE( var% gammzz(mp,nsl) )
ALLOCATE( var% isflag(mp) )
ALLOCATE( var% osnowd(mp) )
ALLOCATE( var% potev(mp) )
ALLOCATE( var% runoff(mp) )
ALLOCATE( var% rnof1(mp) )
ALLOCATE( var% rnof2(mp) )
ALLOCATE( var% rtsoil(mp) )
ALLOCATE( var% rtsoil_expl(mp) )
ALLOCATE( var% sconds(mp,nsnl) )
ALLOCATE( var% sdepth(mp,nsnl) )
ALLOCATE( var% smass(mp,nsnl) )
ALLOCATE( var% snage(mp) )
ALLOCATE( var% snowd(mp) )
ALLOCATE( var% smelt(mp) )
ALLOCATE( var% ssdn(mp,nsnl) )
ALLOCATE( var% ssdnn(mp) )
ALLOCATE( var% tgg(mp,nsl) )
ALLOCATE( var% tggsn(mp,nsnl) )
ALLOCATE( var% tss(mp) )
ALLOCATE( var% tss_p(mp) )
ALLOCATE( var% deltss(mp) )
ALLOCATE( var% owb1(mp) )
ALLOCATE( var% wb(mp,nsl) )
ALLOCATE( var% wbice(mp,nsl) )
ALLOCATE( var% wblf(mp,nsl) )
ALLOCATE( var%wbtot(mp) )
ALLOCATE( var%wbtot1(mp) )
ALLOCATE( var%wbtot2(mp) )
ALLOCATE( var%wb_lake(mp) )
ALLOCATE( var%totwblake(mp) )
ALLOCATE( var%sinfil(mp) )
ALLOCATE( var%evapfbl(mp,nsl) )
ALLOCATE( var%qstss(mp) )
ALLOCATE( var%wetfac(mp) )
ALLOCATE( var%owetfac(mp) )
ALLOCATE( var%t_snwlr(mp) )
ALLOCATE( var%wbfice(mp,nsl) )
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
ALLOCATE( var%tilefrac(mp,ntiles) )
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
ALLOCATE( var%qhlev(mp,nsl+1) )
ALLOCATE( var%satfrac(mp) )
ALLOCATE( var%Qrecharge(mp) )
ALLOCATE( var%rh_srf(mp) )
ALLOCATE( var%rtevap_unsat(mp) )
ALLOCATE( var%rtevap_sat(mp) )
ALLOCATE( var%rt_qh_sublayer(mp) )
!soil moisture variables
ALLOCATE( var%wbeq(mp,nsl) )
ALLOCATE( var%zq(mp,nsl) )
ALLOCATE( var%icefrac(mp,nsl) )
ALLOCATE( var%fracice(mp,nsl) )
ALLOCATE( var%hk(mp,nsl) )
ALLOCATE( var%smp(mp,nsl) )
ALLOCATE( var%dhkdw(mp,nsl) )
ALLOCATE( var%dsmpdw(mp,nsl) )
ALLOCATE( var%wbliq(mp,nsl) )
ALLOCATE( var%wmliq(mp,nsl) )
ALLOCATE( var%wmice(mp,nsl) )
ALLOCATE( var%wmtot(mp,nsl) )

! Allocate variables for SLI soil model:
ALLOCATE( var % s(mp,nsl) )
ALLOCATE( var % Tsoil(mp,nsl) )
ALLOCATE( var % sl(mp) )
ALLOCATE( var % tl(mp) )
ALLOCATE( var % h0(mp) )
ALLOCATE( var % rex(mp,nsl) )
ALLOCATE( var % delwcol(mp) )
ALLOCATE( var % zdelta(mp) )
ALLOCATE( var % kth(mp,nsl) )
ALLOCATE( var % Tsurface(mp) )
ALLOCATE( var % lE(mp) )
ALLOCATE( var % evap(mp) )
ALLOCATE( var % cisoL(mp) )
ALLOCATE( var % rlitt(mp) )
ALLOCATE( var % thetai(mp,nsl) )
ALLOCATE( var % snowliq(mp,SLI_nsnow_max) )
ALLOCATE( var % nsteps(mp) )
ALLOCATE( var % nsnow(mp) )
ALLOCATE( var % TsurfaceFR(mp) )
ALLOCATE( var % Qadv_daily(mp) )
ALLOCATE( var % G0_daily(mp) )
ALLOCATE( var % Qevap_daily(mp) )
ALLOCATE( var % Qprec_daily(mp) )
ALLOCATE( var % Qprec_snow_daily(mp) )

!END IF

END SUBROUTINE alloc_soil_snow_type


END MODULE cable_soil_snow_type_mod
