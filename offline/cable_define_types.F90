!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: defines parameters, variables and derived types, allocation and
!          deallocation of these derived types
!
! Contact: Bernard.Pak@csiro.au
!
! History: Brings together define_dimensions and define_types from v1.4b
!          rs20 now in veg% instead of soil%
!          fes split into fess and fesp (though fes still defined)
!
! Jan 2016: Now includes climate% for use in climate variables required for
! prognostic phenology and potential veg type
! ==============================================================================
!#define UM_BUILD yes
MODULE cable_def_types_mod

  ! Contains all variables which are not subroutine-internal

  IMPLICIT NONE

  SAVE

  PUBLIC

  !---CABLE default KINDs for representing INTEGER/REAL values
  !---at least 10-digit precision

  INTEGER :: mp,    & ! # total no of patches/tiles
       mvtype,& ! total # vegetation types,   from input
#ifdef UM_BUILD 
       mstype=9,& ! total # soil types, needs to de defined atCompile TimeForNow
#else       
       mstype,& ! total # soil types,         from input
#endif
       mland                           ! # land grid cells

  INTEGER, PARAMETER ::                                                        &
       i_d  = KIND(9), &
#ifdef UM_BUILD 
       r_2  = KIND(1.0),&!SELECTED_REAL_KIND(12, 50), &
#else       
       r_2  = KIND(1.d0),&!SELECTED_REAL_KIND(12, 50), &
#endif
       n_tiles = 17,  & ! # possible no of different
       ncp = 3,       & ! # vegetation carbon stores
       ncs = 2,       & ! # soil carbon stores
       mf = 2,        & ! # leaves (sunlit, shaded)
       nrb = 3,       & ! # radiation bands
       msn = 3,       & ! max # snow layers
       swb = 2,       & ! # shortwave bands
       niter = 4,     & ! number of iterations for za/L
                                !      ms = 12          ! # soil layers
       ms = 6         ! # soil layers - standard
  !       ms = 13          ! for Loetschental experiment
  INTEGER, PARAMETER :: n_ktherm = 3
  !   PRIVATE :: r_2, ms, msn, mf, nrb, ncp, ncs

  ! .............................................................................

  ! Energy and water balance variables:
  TYPE balances_type

     REAL, DIMENSION(:), POINTER ::                                           &
          drybal,           & ! energy balance for dry canopy
          ebal,             & ! energy balance per time step (W/m^2)
          ebal_tot,         & ! cumulative energy balance (W/m^2)
          ebal_cncheck,     & ! energy balance consistency check (W/m^2)
          ebal_tot_cncheck, & ! cumulative energy balance (W/m^2)
          ebaltr,           & ! energy balance per time step (W/m^2)
          ebal_tottr,       & ! cumulative energy balance (W/m^2)
          evap_tot,         & ! cumulative evapotranspiration (mm/dels)
          osnowd0,          & ! snow depth, first time step
          precip_tot,       & ! cumulative precipitation (mm/dels)
          rnoff_tot,        & ! cumulative runoff (mm/dels)
          wbal,             & ! water balance per time step (mm/dels)
          wbal_tot,         & ! cumulative water balance (mm/dels)
          wbtot0,           & ! total soil water (mm), first time step
          wetbal,           & ! energy balance for wet canopy
          cansto0,          & ! canopy water storage (mm)
          owbtot,           & ! total soil water (mm), first time step
          evapc_tot,        & ! cumulative evapotranspiration (mm/dels)
          evaps_tot,        & ! cumulative evapotranspiration (mm/dels)
          rnof1_tot,        & ! cumulative runoff (mm/dels)
          rnof2_tot,        & ! cumulative runoff (mm/dels)
          snowdc_tot,       & ! cumulative runoff (mm/dels)
          wbal_tot1,        & ! cumulative water balance (mm/dels)
          delwc_tot,        & ! energy balance for wet canopy
          qasrf_tot,        & ! heat advected to the snow by precip.
          qfsrf_tot,        & ! energy of snowpack phase changes
          qssrf_tot, &        ! energy of snowpack phase changes
          Radbal, &
          EbalSoil, &
          Ebalveg, &
          Radbalsum

  END TYPE balances_type

  ! .............................................................................

  ! Soil parameters:
  TYPE soil_parameter_type

     INTEGER, DIMENSION(:), POINTER ::                                        &
          isoilm     ! integer soil type

     REAL, DIMENSION(:), POINTER ::                                           &
          bch,     & ! parameter b in Campbell equation
          c3,      & ! c3 drainage coeff (fraction)
          clay,    & ! fraction of soil which is clay
          css,     & ! soil specific heat capacity [kJ/kg/K]
          hsbh,    & ! difsat * etasat (=hyds*abs(sucs)*bch)
          hyds,    & ! hydraulic conductivity @ saturation [m/s], Ksat
          i2bp3,   & ! par. one in K vis suction (=nint(bch)+2)
          ibp2,    & ! par. two in K vis suction (fn of pbch)
          rhosoil, & ! soil density [kg/m3]
          sand,    & ! fraction of soil which is sand
          sfc,     & ! vol H2O @ field capacity
          silt,    & ! fraction of soil which is silt
          ssat,    & ! vol H2O @ saturation
          sucs,    & ! suction at saturation (m)
          swilt,   & ! vol H2O @ wilting
          zse,     & ! thickness of each soil layer (1=top) in m
          zshh,    & ! distance between consecutive layer midpoints (m)
                                ! vars intro for Ticket #27
          soilcol, & ! keep color for all patches/tiles
          albsoilf   ! soil reflectance

     REAL(r_2), DIMENSION(:,:), POINTER :: &
          heat_cap_lower_limit

     REAL(r_2), DIMENSION(:,:), POINTER :: &
          zse_vec,css_vec,cnsd_vec

     REAL(r_2), DIMENSION(:), POINTER ::                                      &
          cnsd,    & ! thermal conductivity of dry soil [W/m/K]
          pwb_min    ! working variable (swilt/ssat)**ibp2

     REAL, DIMENSION(:,:), POINTER ::                                         &
          albsoil    ! soil reflectance (2nd dim. BP 21Oct2009)
     !mrd561
     !MD parameters for GW module that vary with soil layer
     REAL(r_2), DIMENSION(:,:), POINTER ::                                    &
          sucs_vec, & !psi at saturation in [mm]
          hyds_vec,  & !saturated hydraulic conductivity  [mm/s]
          bch_vec, & !C and H B [none]
          clay_vec,  & !fraction of soil that is clay [frac]
          sand_vec,  & !fraction of soil that is sand [frac]
          silt_vec,  & !fraction of soil that is silt [frac]
          org_vec,   & !fration of soil made of organic soils [frac]
          rhosoil_vec,& !soil density  [kg/m3]
          ssat_vec, & !volumetric water content at saturation [mm3/mm3]
          watr,   & !residual water content of the soil [mm3/mm3]
          sfc_vec, & !field capcacity (hk = 1 mm/day)
          swilt_vec     ! wilting point (hk = 0.02 mm/day)

     REAL(r_2), DIMENSION(:), POINTER ::                                      &
          drain_dens,&!  drainage density ( mean dist to rivers/streams )
          elev, &  !elevation above sea level
          elev_std, &  !elevation above sea level
          slope,  &  !mean slope of grid cell
          slope_std  !stddev of grid cell slope

     !MD parameters for GW module for the aquifer
     REAL(r_2), DIMENSION(:), POINTER ::                                       &
          GWsucs_vec,  &  !head in the aquifer [mm]
          GWhyds_vec,   &  !saturated hydraulic conductivity of the aquifer [mm/s]
          GWbch_vec,  & !clapp and horn b of the aquifer   [none]
          GWssat_vec,  & !saturated water content of the aquifer [mm3/mm3]
          GWwatr,    & !residual water content of the aquifer [mm3/mm3]
          GWz,       & !node depth of the aquifer    [m]
          GWdz,      & !thickness of the aquifer   [m]
          GWrhosoil_vec    !density of the aquifer substrate [kg/m3]

     ! Additional SLI parameters
     INTEGER,   DIMENSION(:),   POINTER :: nhorizons ! number of soil horizons
     INTEGER,   DIMENSION(:,:), POINTER :: ishorizon ! horizon number 1:nhorizons
     REAL(r_2), DIMENSION(:),   POINTER :: clitt     ! litter (tC/ha)
     REAL(r_2), DIMENSION(:),   POINTER :: zeta      ! macropore parameter
     REAL(r_2), DIMENSION(:),   POINTER :: fsatmax   ! variably saturated area parameter
     !REAL(r_2), DIMENSION(:,:), POINTER :: swilt_vec ! vol H2O @ wilting
     !REAL(r_2), DIMENSION(:,:), POINTER :: ssat_vec  ! vol H2O @ sat
     !REAL(r_2), DIMENSION(:,:), POINTER :: sfc_vec   ! vol H2O @ fc

  END TYPE soil_parameter_type

  ! .............................................................................

  ! Soil and snow variables:
  TYPE soil_snow_type

     INTEGER, DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow

     REAL, DIMENSION(:), POINTER ::                                           &
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
          wbtot1,  & ! total soil water (mm)
          wbtot2,  & ! total soil water (mm)
          wb_lake, &
          totwblake, & !daily integrated wb_lake: used in ACCESS
          sinfil,  &
          qstss,   &
          wetfac,  & ! surface wetness fact. at current time step
          owetfac, & ! surface wetness fact. at previous time step
          t_snwlr, & ! top snow layer depth in 3 layer snowpack
          tggav,   & ! mean soil temperature in K
          otgg,    & ! soil temperature in K
          otss,    & ! surface temperature (weighted soil, snow)
          otss_0,  & ! surface temperature (weighted soil, snow)
          tprecip, &
          tevap,   &
          trnoff,  &
          totenbal,&!
          totenbal2,&
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

     REAL, DIMENSION(:,:), POINTER ::                                         &
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


     REAL(r_2), DIMENSION(:), POINTER ::                                      &
          wbtot   ! total soil water (mm)

     REAL(r_2), DIMENSION(:,:), POINTER ::                                    &
          gammzz,  & ! heat capacity for each soil layer
          wb,      & ! volumetric soil moisture (solid+liq)
          wbice,   & ! soil ice
          wblf,    & !
          wbfice     !

     !mrd561
     !MD variables for the revised soil moisture + GW scheme
     REAL(r_2), DIMENSION(:), POINTER   ::                                     &
          GWwb,    &  ! water content in aquifer [mm3/mm3]
          GWhk,    &  ! aquifer hydraulic conductivity  [mm/s]
          GWdhkdw, &  ! aquifer d(hk) over d(water content) [(mm/s)/(mm3/mm3)]
          GWdsmpdw,&  ! aquifer d(smp) / dw   [(mm)/(mm3/mm3)]
          wtd,     &  ! water table depth   [mm]
          GWsmp,   &  ! aquifer soil matric potential [mm]
          GWwbeq,  &  ! equilibrium aquifer water content [mm3/mm3]
          GWzq,    &  ! equilibrium aquifer smp   [mm]
          qhz,     &  ! horizontal hydraulic conductivity in 1D gw model for soil layers  [mm/s]
          satfrac, &
          Qrecharge,&
          rh_srf,   &
          rtevap_sat,&
          rtevap_unsat,&
          rt_qh_sublayer

     REAL(r_2), DIMENSION(:,:), POINTER  ::                                     &
          wbeq,    &    ! equilibrium water content [mm3/mm3]
          zq,      &    ! equilibrium smp       [mm]
          icefrac, &    ! ice fraction  [none]  -> ice mass / total mass
          fracice, &    ! alternate ice fraction  [none] - parameterized
          hk,      &    ! hydraulic conductivity for soil layers [mm/s]
          smp,     &    ! soil matric potential for soil layers         [mm]
          dhkdw,   &    ! d(hydraulic conductivity ) d(water) for soil layers [(mm/s)/(mm3/mm3)]
          dsmpdw,  &    ! d(smp)/ d(water) for soil layers   [(mm)/(mm3/mm3)]
          wbliq,   &    ! volumetric liquid water content  [mm3/mm3]
          wmliq,   &    !water mass [mm] liq
          wmice,   &    !water mass [mm] ice
          wmtot,   &    !water mass [mm] liq+ice ->total
          qhlev
     ! Additional SLI variables:
     REAL(r_2), DIMENSION(:,:), POINTER :: S         ! moisture content relative to sat value    (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil         !     Tsoil (deg C)
     REAL(r_2), DIMENSION(:),   POINTER :: SL        ! litter moisture content relative to sat value (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:),   POINTER :: TL        ! litter temperature in K     (edit vh 23/01/08)
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

  ! .............................................................................

  ! Vegetation parameters:
  TYPE veg_parameter_type

     INTEGER, DIMENSION(:), POINTER ::                                        &
          iveg , &      ! vegetation type
          iLU ! land use type
     REAL, DIMENSION(:), POINTER ::                                           &
          canst1,  & ! max intercepted water by canopy (mm/LAI)
          dleaf,   & ! chararacteristc legnth of leaf (m)
          ejmax,   & ! max pot. electron transp rate top leaf(mol/m2/s)
          meth,    & ! method for calculation of canopy fluxes and temp.
          frac4,   & ! fraction of c4 plants
          hc,      & ! roughness height of canopy (veg - snow)
          vlai,    & ! leaf area index
          xalbnir, &
          rp20,    & ! plant respiration coefficient at 20 C
          rpcoef,  & ! temperature coef nonleaf plant respiration (1/C)
          rs20,    & ! soil respiration at 20 C [mol m-2 s-1]
          shelrb,  & ! sheltering factor (dimensionless)
          vegcf,   & ! kdcorbin, 08/10
          tminvj,  & ! min temperature of the start of photosynthesis
          toptvj,  & ! opt temperature of the start of photosynthesis
          tmaxvj,  & ! max temperature of the start of photosynthesis
          vbeta,   & !
          vcmax,   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
          xfang,   & ! leaf angle PARAMETER
          extkn,   & ! extinction coef for vertical
          vlaimax, & ! extinction coef for vertical
          wai,     & ! wood area index (stem+branches+twigs)
          a1gs,    & ! a1 parameter in stomatal conductance model
          d0gs,    & ! d0 in stomatal conductance model
          alpha,   & ! initial slope of J-Q response curve
          convex,  & ! convexity of J-Q response curve
          cfrd,    & ! ratio of day respiration to vcmax
          gswmin,  & ! minimal stomatal conductance
          conkc0,  &  ! Michaelis-menton constant for carboxylase
          conko0,  &  ! Michaelis-menton constant for oxygenase
          ekc,     &  ! activation energy for caroxylagse
          eko,     &  ! acvtivation enegery for oxygenase
          g0,      & ! Belinda's stomatal model intercept, Ticket #56.
          g1         ! Belinda's stomatal model slope, Ticket #56.

     LOGICAL, DIMENSION(:), POINTER ::                                        &
          deciduous ! flag used for phenology fix

     REAL, DIMENSION(:,:), POINTER ::                                         &
          refl,    &
          taul,    &
          froot      ! fraction of root in each soil layer

     ! Additional  veg parameters:
     REAL(r_2), DIMENSION(:), POINTER :: rootbeta ! parameter for estimating vertical root mass distribution (froot)
     REAL(r_2), DIMENSION(:), POINTER :: gamma    ! parameter in root efficiency function (Lai and Katul 2000)
     REAL(r_2), DIMENSION(:), POINTER :: ZR       ! maximum rooting depth (cm)
     REAL(r_2), DIMENSION(:), POINTER :: F10      ! fraction of roots in top 10 cm

     REAL(r_2), DIMENSION(:), POINTER :: clitt     !

     ! Additional POP veg param
     INTEGER, DIMENSION(:,:), POINTER ::  disturbance_interval
     REAL(r_2), DIMENSION(:,:), POINTER ::  disturbance_intensity

  END TYPE veg_parameter_type

  ! .............................................................................

  ! Canopy/vegetation variables:
  TYPE canopy_type


     REAL, DIMENSION(:), POINTER ::                                           &
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
          fhs_cor, &
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
     REAL, DIMENSION(:), POINTER ::                                           &
          fns_cor, & ! correction to net rad avail to soil (W/m2)
          ga_cor  ! correction to ground heat flux (W/m2)

     REAL, DIMENSION(:,:), POINTER ::                                         &
          evapfbl, &
          gswx,    & ! stom cond for water
          zetar, &   ! stability parameter (ref height)
                                !! vh_js !!
          zetash      ! stability parameter (shear height)

     REAL(r_2), DIMENSION(:), POINTER ::                                      &
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

     REAL(r_2), DIMENSION(:), POINTER :: &
          sublayer_dz

     ! Additional variables:
     REAL(r_2), DIMENSION(:,:),   POINTER :: gw     ! dry canopy conductance (ms-1) edit vh 6/7/09
     REAL(r_2), DIMENSION(:,:,:), POINTER :: ancj   ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
     REAL(r_2), DIMENSION(:,:),   POINTER :: tlfy   ! sunlit and shaded leaf temperatures
     REAL(r_2), DIMENSION(:,:),   POINTER :: ecy    ! sunlit and shaded leaf transpiration (dry canopy)
     REAL(r_2), DIMENSION(:,:),   POINTER :: ecx    ! sunlit and shaded leaf latent heat flux
     REAL(r_2), DIMENSION(:,:,:), POINTER :: ci     ! intra-cellular CO2 vh 6/7/09
     REAL(r_2), DIMENSION(:),     POINTER :: fwsoil !

     !! vh_js !! !litter thermal conductivity (Wm-2K-1) and vapour diffusivity (m2s-1)
     REAL(r_2), DIMENSION(:), POINTER :: kthLitt, DvLitt


  END TYPE canopy_type

  ! .............................................................................

  ! Radiation variables:
  TYPE radiation_type

     REAL, DIMENSION(:), POINTER   ::                                         &
          transb,  & ! fraction SW beam tranmitted through canopy
          albedo_T,& ! canopy+soil albedo for VIS+NIR
          longitude,&! longitude
          workp1,  & ! absorbed short-wave radiation for soil
          workp2,  & ! absorbed short-wave radiation for soil
          workp3,  & ! absorbed short-wave radiation for soil
          extkb,   & ! beam radiation extinction coeff
          extkd2,  & ! diffuse 2D radiation extinction coeff
          extkd,   & ! diffuse radiation extinction coeff (-)
          flws,    & ! soil long-wave radiation
          latitude,& ! latitude
          lwabv,   & ! long wave absorbed by vegetation
          qssabs,  & ! absorbed short-wave radiation for soil
          transd,  & ! frac SW diffuse transmitted through canopy
          trad,    & !  radiative temperature (soil and veg)
          otrad      ! radiative temperature on previous timestep (ACCESS)

     REAL, DIMENSION(:,:), POINTER  ::                                        &
          fvlai,   & ! leaf area index of big leaf
          rhocdf,  & ! canopy diffuse reflectance (-)
          rniso,   & ! sum(rad%qcan, 3) total abs by canopy (W/m2)
          scalex,  & ! scaling PARAMETER for big leaf
          albedo,  & ! canopy+soil albedo
          reffdf,  & ! effective conopy diffuse reflectance
          reffbm,  & ! effective conopy beam reflectance
          extkbm,  & ! modified k beam(6.20)(for leaf scattering)
          extkdm,  & ! modified k diffuse(6.20)(for leaf scattering)
          fbeam,   & ! beam fraction
          cexpkbm, & ! canopy beam transmittance
          cexpkdm, & ! canopy diffuse transmittance
          rhocbm,  & ! modified canopy beam reflectance(6.21)
          gradis     ! radiative conductance

     REAL, DIMENSION(:,:,:), POINTER ::                                       &
          qcan ! absorbed radiation for canopy (W/m^2)


  END TYPE radiation_type

  ! .............................................................................

  ! Roughness variables:
  TYPE roughness_type

     REAL, DIMENSION(:), POINTER ::                                           &
          disp,    & ! zero-plane displacement
          hruff,   & ! canopy height above snow level
          hruff_grmx,&! max ht of canopy from tiles on same grid
          rt0us,   & ! eq. 3.54, SCAM manual (CSIRO tech report 132)
          rt1usa,  & ! resistance from disp to hruf
          rt1usb,  & ! resist fr hruf to zruffs (zref if zref<zruffs)
          rt1,     & ! 1/aerodynamic conductance
          za_uv,   & ! level of lowest atmospheric model layer
          za_tq,   & ! level of lowest atmospheric model layer
          z0m,     & ! roughness length
          zref_uv, & ! Reference height for met forcing
          zref_tq, & ! Reference height for met forcing
          zruffs,  & ! SCALAR Roughness sublayer depth (ground=origin)
          z0soilsn,& ! roughness length of bare soil surface
          z0soil     ! roughness length of bare soil surface

     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL, DIMENSION(:), POINTER ::                                           &
          coexp ! Extinction coef for wind profile in canopy

     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL, DIMENSION(:), POINTER ::                                           &
          usuh ! Friction velocity/windspeed at canopy height

     REAL, DIMENSION(:), POINTER ::                                           &
          term2, term3, term5, term6, term6a ! for aerodyn resist. calc.



  END TYPE roughness_type

  ! .............................................................................

  ! Air variables:
  TYPE air_type

     REAL, DIMENSION(:), POINTER ::                                           &
          rho,     & ! dry air density (kg m-3)
          volm,    & ! molar volume (m3 mol-1)
          rlam,    & ! latent heat for water (j/kg)
          qsat,    & ! saturation specific humidity
          epsi,    & ! d(qsat)/dT ((kg/kg)/K)
          visc,    & ! air kinematic viscosity (m2/s)
          psyc,    & ! psychrometric constant
          dsatdk,  & ! d(es)/dT (mb/K)
          cmolar     ! conv. from m/s to mol/m2/s

  END TYPE air_type

  ! .............................................................................

  ! Meterological data:
  TYPE met_type

     INTEGER, DIMENSION(:), POINTER ::                                        &
          year,    & ! local time year AD
          moy        ! local time month of year

     REAL, DIMENSION(:), POINTER ::                                           &
          ca,      & ! CO2 concentration (mol/mol)
          doy,     & ! local time day of year = days since 0 hr 1st Jan
          hod,     & ! local hour of day
          ofsd,    & ! downward short-wave radiation (W/m2)
          fld,     & ! downward long-wave radiation (W/m2)
          precip,  & ! rainfall (liquid+solid)(mm/dels)
          precip_sn,&! solid preipitation only (mm/dels)
          tk,      & ! surface air temperature (oK)
          tvair,   & ! within canopy air temperature (oK)
          tvrad,   & ! radiative vegetation temperature (K)
          pmb,     & ! surface air pressure (mbar)
          ua,      & ! surface wind speed (m/s)
          qv,      & ! surface specific humidity (g/g)
          qvair,   & ! within canopy specific humidity (g/g)
          da,      & ! water vap pressure deficit at ref height (Pa)
          dva,     & ! in canopy water vap pressure deficit (Pa)
          coszen,   &  ! cos(zenith angle of sun)
          Ndep,     &   ! nitrogen deposition (gN m-2 d-1)
          Pdep ! P deposition (gP m-2 d-1)

     REAL, DIMENSION(:,:), POINTER ::                                         &
          fsd  ! downward short-wave radiation (W/m2)

  END TYPE met_type

  ! .............................................................................

  ! Climate data:
  TYPE climate_type

     INTEGER :: nyear_average = 20
     INTEGER :: nday_average  = 31
     !      INTEGER, POINTER ::                                                  &
     INTEGER ::                                                  &
          nyears, & ! number of years in climate record
          doy ! day of year

     INTEGER, DIMENSION(:), POINTER ::                                   &
          chilldays, &   ! length of chilling period (period with T<5deg)
          iveg, &        ! potential vegetation type based on climatic constraints
          biome

     REAL, DIMENSION(:), POINTER ::                                           &
          dtemp,        & ! daily temperature
          dmoist,        & ! daily moisture availability
          mtemp,       & ! mean temperature over the last 31 days
          qtemp,       & ! mean temperature over the last 91 days
          mmoist,        & ! monthly moisture availability
          mtemp_min,   & ! minimum monthly temperature
          mtemp_max,   & ! maximum monhtly temperature
          qtemp_max,   & ! mean temperature of the warmest quarter (so far this year)
          qtemp_max_last_year,   & ! mean temperature of the warmest quarter (last calendar year)
          mtemp_min20,   & ! minimum monthly temperature, averaged over 20 y
          mtemp_max20,   & ! maximum monhtly temperature, averaged over 20 y
          atemp_mean,  & ! annual average temperature
          AGDD5,       &
          GDD5,        & ! growing degree day sum relative to 5deg base temperature
          AGDD0,        & !
          GDD0,        & ! growing degree day sum relative to 0deg base temperature
          alpha_PT,    & ! ratio of annual evap to annual PT evap
          evap_PT,    & ! annual PT evap [mm]
          aevap , &       ! annual evap [mm]
          alpha_PT20

     REAL, DIMENSION(:,:), POINTER ::                                   &
          mtemp_min_20, & ! mimimum monthly temperatures for the last 20 y
          mtemp_max_20, & ! maximum monthly temperatures for the last 20 y
          dtemp_31 , &    ! daily temperature for the last 31 days
          dmoist_31 , &    ! daily moisture availability for the last 31 days
          alpha_PT_20, &      ! priestley Taylor Coefft for last 20 y
          dtemp_91     ! daily temperature for the last 91 days

  END TYPE climate_type

  ! .............................................................................

  ! Cumulative flux variables:
  TYPE sum_flux_type

     REAL, DIMENSION(:), POINTER ::                                           &
          sumpn,   & ! sum of canopy photosynthesis (g C m-2)
          sumrp,   & ! sum of plant respiration (g C m-2)
          sumrpw,  & ! sum of plant respiration (g C m-2)
          sumrpr,  & ! sum of plant respiration (g C m-2)
          sumrs,   & ! sum of soil respiration (g C m-2)
          sumrd,   & ! sum of daytime respiration (g C m-2)
          dsumpn,  & ! daily sumpn
          dsumrp,  & ! daily sumrp
          dsumrs,  & ! daily sumrs
          dsumrd,  & ! daily sumrd
          sumxrp,  & ! sum plant resp. modifier
          sumxrs     ! sum soil resp. modifier

  END TYPE sum_flux_type

  ! .............................................................................

  TYPE bgc_pool_type

     REAL, DIMENSION(:,:), POINTER ::                                         &
          cplant,  & ! plant carbon (g C/m2))
          csoil   ! soil carbon (g C/m2)


     REAL, DIMENSION(ncp)  :: ratecp ! plant carbon rate constant (1/year)

     REAL, DIMENSION(ncs)  :: ratecs ! soil carbon rate constant (1/year)

  END TYPE bgc_pool_type

  ! .............................................................................

  ! Functions for allocating these types
  ! All overloaded so code only needs to call alloc_cbm_var
  ! Alloc routines could all initialise to NaN or zero for debugging?
  ! Don't need the mp argument here as it's a module variable.
  PUBLIC :: alloc_cbm_var
  PRIVATE :: alloc_bgc_pool_type, dealloc_bgc_pool_type

  INTERFACE alloc_cbm_var
     MODULE PROCEDURE alloc_balances_type,                                    &
          alloc_soil_parameter_type,                                            &
          alloc_soil_snow_type,                                                 &
          alloc_veg_parameter_type,                                             &
          alloc_canopy_type,                                                    &
          alloc_radiation_type,                                                 &
          alloc_roughness_type,                                                 &
          alloc_air_type,                                                       &
          alloc_met_type,                                                       &
          alloc_sum_flux_type,                                                  &
          alloc_bgc_pool_type ,                                                 &
          alloc_climate_type
  END INTERFACE

  INTERFACE dealloc_cbm_var
     MODULE PROCEDURE dealloc_balances_type,                                  &
          dealloc_soil_parameter_type,                                          &
          dealloc_soil_snow_type,                                               &
          dealloc_veg_parameter_type,                                           &
          dealloc_canopy_type,                                                  &
          dealloc_radiation_type,                                               &
          dealloc_roughness_type,                                               &
          dealloc_air_type,                                                     &
          dealloc_met_type,                                                     &
          dealloc_sum_flux_type,                                                &
          dealloc_bgc_pool_type
  END INTERFACE


CONTAINS

  SUBROUTINE alloc_balances_type(var, mp)

    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

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

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_soil_parameter_type(var, mp)

    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    ALLOCATE( var% bch(mp) )
    ALLOCATE( var% c3(mp) )
    ALLOCATE( var% clay(mp) )
    ALLOCATE( var% css(mp) )
    ALLOCATE( var% hsbh(mp) )
    ALLOCATE( var% hyds(mp) )
    ALLOCATE( var% i2bp3(mp) )
    ALLOCATE( var% ibp2(mp) )
    ALLOCATE( var% isoilm(mp) )
    ALLOCATE( var% rhosoil(mp) )
    ALLOCATE( var% sand(mp) )
    ALLOCATE( var% sfc(mp) )
    ALLOCATE( var% silt(mp) )
    ALLOCATE( var% ssat(mp) )
    ALLOCATE( var% sucs(mp) )
    ALLOCATE( var% swilt(mp) )
    ALLOCATE( var% zse(ms) )
    ALLOCATE( var% zshh(ms+1) )
    ALLOCATE( var% cnsd(mp) )
    ALLOCATE( var% albsoil(mp, nrb) )
    ALLOCATE( var% pwb_min(mp) )
    ALLOCATE( var% albsoilf(mp) )
    ALLOCATE( var% soilcol(mp) )
    !mrd561
    !MD
    !Aquifer properties
    ALLOCATE( var%GWhyds_vec(mp) )
    ALLOCATE( var%GWsucs_vec(mp) )
    ALLOCATE( var%GWbch_vec(mp) )
    ALLOCATE( var%GWssat_vec(mp) )
    ALLOCATE( var%GWwatr(mp) )
    var%GWwatr(:) = 0.05
    ALLOCATE( var%GWz(mp) )
    ALLOCATE( var%GWdz(mp) )
    ALLOCATE( var%GWrhosoil_vec(mp) )
    !soil properties (vary by layer)
    ALLOCATE( var% zse_vec(mp,ms) )
    ALLOCATE( var% heat_cap_lower_limit(mp,ms) )
    ALLOCATE( var% css_vec(mp,ms) )
    ALLOCATE( var% cnsd_vec(mp,ms) )
    ALLOCATE( var%hyds_vec(mp,ms) )
    ALLOCATE( var%sucs_vec(mp,ms) )
    ALLOCATE( var%bch_vec(mp,ms) )
    ALLOCATE( var%ssat_vec(mp,ms) )
    ALLOCATE( var%watr(mp,ms) )
    var%watr(:,:) = 0.05
    ALLOCATE( var%sfc_vec(mp,ms) )
    ALLOCATE( var%swilt_vec(mp,ms) )
    ALLOCATE( var%sand_vec(mp,ms) )
    ALLOCATE( var%clay_vec(mp,ms) )
    ALLOCATE( var%silt_vec(mp,ms) )
    ALLOCATE( var%org_vec(mp,ms) )
    ALLOCATE( var%rhosoil_vec(mp,ms) )

    ALLOCATE( var%drain_dens(mp) )
    ALLOCATE( var%elev(mp) )
    ALLOCATE( var%elev_std(mp) )
    ALLOCATE( var%slope(mp) )
    ALLOCATE( var%slope_std(mp) )

    ! Allocate variables for SLI soil model:
    ALLOCATE ( var % nhorizons(mp) )
    ALLOCATE ( var % ishorizon(mp,ms) )
    ALLOCATE ( var % clitt(mp) )
    ALLOCATE ( var % zeta(mp) )
    ALLOCATE ( var % fsatmax(mp) )
    !ALLOCATE ( var % swilt_vec(mp,ms) )
    !ALLOCATE ( var % ssat_vec(mp,ms) )
    !ALLOCATE ( var % sfc_vec(mp,ms) )
    IF(.NOT.(ASSOCIATED(var % swilt_vec))) ALLOCATE ( var % swilt_vec(mp,ms) )
    IF(.NOT.(ASSOCIATED(var % ssat_vec))) ALLOCATE ( var % ssat_vec(mp,ms) )
    IF(.NOT.(ASSOCIATED(var % sfc_vec))) ALLOCATE ( var % sfc_vec(mp,ms) )


  END SUBROUTINE alloc_soil_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_soil_snow_type(var, mp)

    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    ALLOCATE ( var % iantrct(mp) )
    ALLOCATE ( var % pudsto(mp) )
    ALLOCATE ( var % pudsmx(mp) )
    ALLOCATE ( var % dtmlt(mp,3) )
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
    ALLOCATE( var% gammzz(mp,ms) )
    ALLOCATE( var% isflag(mp) )
    ALLOCATE( var% osnowd(mp) )
    ALLOCATE( var% potev(mp) )
    ALLOCATE( var% runoff(mp) )
    ALLOCATE( var% rnof1(mp) )
    ALLOCATE( var% rnof2(mp) )
    ALLOCATE( var% rtsoil(mp) )
    ALLOCATE( var% sconds(mp,msn) )
    ALLOCATE( var% sdepth(mp,msn) )
    ALLOCATE( var% smass(mp,msn) )
    ALLOCATE( var% snage(mp) )
    ALLOCATE( var% snowd(mp) )
    ALLOCATE( var% smelt(mp) )
    ALLOCATE( var% ssdn(mp,msn) )
    ALLOCATE( var% ssdnn(mp) )
    ALLOCATE( var% tgg(mp,ms) )
    ALLOCATE( var% tggsn(mp,msn) )
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
    ALLOCATE ( var % S(mp,ms) )
    ALLOCATE ( var % Tsoil(mp,ms) )
    ALLOCATE ( var % SL(mp) )
    ALLOCATE ( var % TL(mp) )
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

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_veg_parameter_type(var, mp)

    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    ALLOCATE( var% canst1(mp) )
    ALLOCATE( var% dleaf(mp) )
    ALLOCATE( var% ejmax(mp) )
    ALLOCATE( var% iveg(mp) )
    ALLOCATE( var% iLU(mp) )
    ALLOCATE( var% meth(mp) )
    ALLOCATE( var% frac4(mp) )
    ALLOCATE( var% hc(mp) )
    ALLOCATE( var% vlai(mp) )
    ALLOCATE( var% xalbnir(mp) )
    ALLOCATE( var% rp20(mp) )
    ALLOCATE( var% rpcoef(mp) )
    ALLOCATE( var% rs20(mp) )
    ALLOCATE( var% shelrb(mp) )
    ALLOCATE( var% vegcf(mp) )
    ALLOCATE( var% tminvj(mp) )
    ALLOCATE( var% toptvj(mp) )
    ALLOCATE( var% tmaxvj(mp) )
    ALLOCATE( var% vbeta(mp) )
    ALLOCATE( var% vcmax(mp) )
    ALLOCATE( var% xfang(mp) )
    ALLOCATE( var%extkn(mp) )
    ALLOCATE( var%wai(mp) )
    ALLOCATE( var%deciduous(mp) )
    ALLOCATE( var%froot(mp,ms) )
    !was nrb(=3), but never uses (:,3) in model
#ifdef UM_BUILD 
    ALLOCATE( var%refl(mp,nrb) ) !jhan:swb?
    ALLOCATE( var%taul(mp,nrb) )
#else
    ALLOCATE( var%refl(mp,2) ) !jhan:swb?
    ALLOCATE( var%taul(mp,2) )
#endif
    ALLOCATE( var%vlaimax(mp) )
    ALLOCATE( var%a1gs(mp) )
    ALLOCATE( var%d0gs(mp) )
    ALLOCATE( var%alpha(mp) )
    ALLOCATE( var%convex(mp) )
    ALLOCATE( var%cfrd(mp) )
    ALLOCATE( var%gswmin(mp) )
    ALLOCATE( var%conkc0(mp) )
    ALLOCATE( var%conko0(mp) )
    ALLOCATE( var%ekc(mp) )
    ALLOCATE( var%eko(mp) )
    ALLOCATE( var% g0(mp) )   ! Ticket #56.
    ALLOCATE( var% g1(mp) )   ! Ticket #56.


    ALLOCATE ( var % rootbeta(mp) )
    ALLOCATE ( var % gamma(mp) )
    ALLOCATE ( var % F10(mp) )
    ALLOCATE ( var % ZR(mp) )
    ALLOCATE ( var % clitt(mp) )

    ALLOCATE ( var % disturbance_interval(mp,2) )
    ALLOCATE ( var % disturbance_intensity(mp,2) )



  END SUBROUTINE alloc_veg_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_canopy_type(var, mp)

    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

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
    ALLOCATE( var% zetar(mp,NITER) )
    ALLOCATE( var% zetash(mp,NITER) )
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

    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    ALLOCATE( var% albedo(mp,nrb) )
    ALLOCATE( var% extkb(mp) )
    ALLOCATE( var% extkd2(mp) )
    ALLOCATE( var% extkd(mp) )
    ALLOCATE( var% flws(mp) )
    ALLOCATE( var% fvlai(mp,mf) )
    ALLOCATE( var% latitude(mp) )
    ALLOCATE( var% lwabv(mp) )
    ALLOCATE( var% qcan(mp,mf,nrb) )
    ALLOCATE( var% qssabs(mp) )
    ALLOCATE( var% rhocdf(mp,nrb) )
    ALLOCATE( var% rniso(mp,mf) )
    ALLOCATE( var% scalex(mp,mf) )
    ALLOCATE( var% transd(mp) )
    ALLOCATE( var% trad(mp) )
    ALLOCATE( var% otrad(mp) )
    ALLOCATE( var% reffdf(mp,nrb) )
    ALLOCATE( var% reffbm(mp,nrb) )
    ALLOCATE( var% extkbm(mp,nrb) )
    ALLOCATE( var% extkdm(mp,nrb) )
    ALLOCATE( var% cexpkbm(mp,swb) )
    ALLOCATE( var% cexpkdm(mp,swb) )
    ALLOCATE( var% fbeam(mp,nrb) )
    ALLOCATE( var% rhocbm(mp,nrb) )
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

    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

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

    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

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

    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

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

  SUBROUTINE alloc_climate_type(var, mp)

    IMPLICIT NONE

    TYPE(climate_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    INTEGER :: ny, nd
    ny = var%nyear_average
    nd = var%nday_average
    PRINT*, 'ny', ny
    PRINT*, 'nd', nd
    !   ALLOCATE ( var %  nyears )
    !   ALLOCATE ( var %  doy )
    ALLOCATE ( var %  dtemp(mp) )
    ALLOCATE ( var %  dmoist(mp) )
    ALLOCATE ( var % mtemp(mp) )
    ALLOCATE ( var % qtemp(mp) )
    ALLOCATE ( var % mmoist(mp) )
    ALLOCATE ( var % mtemp_min(mp) )
    ALLOCATE ( var %  mtemp_max20(mp) )
    ALLOCATE ( var % mtemp_min20(mp) )
    ALLOCATE ( var %  mtemp_max(mp) )
    ALLOCATE ( var %  qtemp_max(mp) )
    ALLOCATE ( var %  qtemp_max_last_year(mp) )
    ALLOCATE ( var % atemp_mean(mp) )
    ALLOCATE ( var % AGDD5(mp) )
    ALLOCATE ( var % GDD5(mp) )
    ALLOCATE ( var % AGDD0(mp) )
    ALLOCATE ( var % GDD0(mp) )
    ALLOCATE ( var % chilldays(mp) )
    ALLOCATE ( var % iveg(mp) )
    ALLOCATE ( var % biome(mp) )
    ALLOCATE ( var % alpha_PT(mp) )
    ALLOCATE ( var % alpha_PT20(mp) )
    ALLOCATE ( var % evap_PT(mp) )
    ALLOCATE ( var % aevap(mp) )

    ALLOCATE ( var % mtemp_min_20(mp,ny) )
    ALLOCATE ( var %     mtemp_max_20(mp,ny) )
    ALLOCATE ( var %     dtemp_31(mp,nd) )
    ALLOCATE ( var %     dmoist_31(mp,nd) )
    ALLOCATE ( var %     dtemp_91(mp,91) )
    ALLOCATE ( var %     alpha_PT_20(mp,ny) )

  END SUBROUTINE alloc_climate_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_sum_flux_type(var, mp)

    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

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

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_bgc_pool_type(var, mp)

    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    ALLOCATE ( var % cplant(mp,ncp) )
    ALLOCATE ( var % csoil(mp,ncs) )

  END SUBROUTINE alloc_bgc_pool_type

  ! ------------------------------------------------------------------------------

  ! Begin deallocation routines:
  SUBROUTINE dealloc_balances_type(var)

    TYPE(balances_type), INTENT(inout) :: var

    DEALLOCATE( var% drybal )
    DEALLOCATE( var% ebal )
    DEALLOCATE( var% ebal_tot )
    DEALLOCATE( var% ebaltr )
    DEALLOCATE( var% ebal_tottr )
    DEALLOCATE( var% ebal_cncheck )
    DEALLOCATE( var% ebal_tot_cncheck )
    DEALLOCATE( var% evap_tot)
    DEALLOCATE( var% osnowd0 )
    DEALLOCATE( var% precip_tot )
    DEALLOCATE( var% rnoff_tot )
    DEALLOCATE( var% wbal )
    DEALLOCATE( var% wbal_tot )
    DEALLOCATE( var% wbtot0 )
    DEALLOCATE( var% wetbal )
    DEALLOCATE( var% cansto0 )
    DEALLOCATE( var% evapc_tot )
    DEALLOCATE( var% evaps_tot )
    DEALLOCATE( var% rnof1_tot )
    DEALLOCATE( var% rnof2_tot )
    DEALLOCATE( var% snowdc_tot )
    DEALLOCATE( var% wbal_tot1 )
    DEALLOCATE( var% owbtot )
    DEALLOCATE( var% delwc_tot )
    DEALLOCATE( var% qasrf_tot )
    DEALLOCATE( var% qfsrf_tot )
    DEALLOCATE( var% qssrf_tot )

    DEALLOCATE( var% Radbal )
    DEALLOCATE( var% Ebalsoil )
    DEALLOCATE( var% Ebalveg )
    DEALLOCATE( var% Radbalsum )

  END SUBROUTINE dealloc_balances_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_soil_parameter_type(var)

    TYPE(soil_parameter_type), INTENT(inout) :: var

    DEALLOCATE( var% bch )
    DEALLOCATE( var% c3 )
    DEALLOCATE( var% clay )
    DEALLOCATE( var% css )
    DEALLOCATE( var% hsbh )
    DEALLOCATE( var% hyds )
    DEALLOCATE( var% i2bp3 )
    DEALLOCATE( var% ibp2 )
    DEALLOCATE( var% isoilm )
    DEALLOCATE( var% rhosoil )
    DEALLOCATE( var% sand )
    DEALLOCATE( var% sfc )
    DEALLOCATE( var% silt )
    DEALLOCATE( var% ssat )
    DEALLOCATE( var% sucs )
    DEALLOCATE( var% swilt )
    DEALLOCATE( var% zse )
    DEALLOCATE( var% zshh )
    DEALLOCATE( var% cnsd )
    DEALLOCATE( var% albsoil )
    DEALLOCATE( var% cnsd )
    DEALLOCATE( var% pwb_min)
    DEALLOCATE( var% albsoilf )
    DEALLOCATE( var% soilcol )
    !mrd561
    !MD
    !Aquifer properties
    DEALLOCATE( var%GWhyds_vec )
    DEALLOCATE( var%GWsucs_vec )
    DEALLOCATE( var%GWbch_vec )
    DEALLOCATE( var%GWssat_vec )
    DEALLOCATE( var%GWwatr )
    DEALLOCATE( var%GWz )
    DEALLOCATE( var%GWdz )
    DEALLOCATE( var%GWrhosoil_vec )
    !soil properties (vary by layer)
    DEALLOCATE( var% zse_vec )
    DEALLOCATE( var% heat_cap_lower_limit )
    DEALLOCATE( var% css_vec )
    DEALLOCATE( var% cnsd_vec )
    DEALLOCATE( var%hyds_vec )
    DEALLOCATE( var%sucs_vec )
    DEALLOCATE( var%bch_vec )
    DEALLOCATE( var%ssat_vec )
    DEALLOCATE( var%watr )
    DEALLOCATE( var%sfc_vec )
    DEALLOCATE( var%swilt_vec )
    DEALLOCATE( var%sand_vec )
    DEALLOCATE( var%clay_vec )
    DEALLOCATE( var%silt_vec )
    DEALLOCATE( var%org_vec  )
    DEALLOCATE( var%rhosoil_vec )
    DEALLOCATE( var%drain_dens )
    DEALLOCATE( var%elev )
    DEALLOCATE( var%elev_std )
    DEALLOCATE( var%slope )
    DEALLOCATE( var%slope_std )
    ! Deallocate variables for SLI soil model:
    !IF(cable_user%SOIL_STRUC=='sli') THEN
    DEALLOCATE ( var % nhorizons)
    DEALLOCATE ( var % ishorizon)
    DEALLOCATE ( var % clitt )
    DEALLOCATE ( var % zeta )
    DEALLOCATE ( var % fsatmax )
    !DEALLOCATE ( var % swilt_vec )
    !DEALLOCATE ( var % ssat_vec )
    !DEALLOCATE ( var % sfc_vec )
    IF(ASSOCIATED(var % swilt_vec)) DEALLOCATE ( var % swilt_vec )
    IF(ASSOCIATED(var % ssat_vec)) DEALLOCATE ( var % ssat_vec )
    IF(ASSOCIATED(var % sfc_vec)) DEALLOCATE ( var % sfc_vec )
    !END IF


  END SUBROUTINE dealloc_soil_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_soil_snow_type(var)

    TYPE(soil_snow_type), INTENT(inout) :: var

    DEALLOCATE ( var % iantrct )
    DEALLOCATE ( var % pudsto )
    DEALLOCATE ( var % pudsmx )
    DEALLOCATE ( var % dtmlt )
    DEALLOCATE( var% albsoilsn )
    DEALLOCATE( var% cls )
    DEALLOCATE( var% dfn_dtg )
    DEALLOCATE( var% dfh_dtg )
    DEALLOCATE( var% dfe_ddq )
    DEALLOCATE( var% ddq_dtg )
    DEALLOCATE( var% dfe_dtg )  !REV_CORR variable
    DEALLOCATE( var% evapsn )
    DEALLOCATE( var% fwtop )
    DEALLOCATE( var% fwtop1 )
    DEALLOCATE( var% fwtop2 )
    DEALLOCATE( var% fwtop3 )
    DEALLOCATE( var% gammzz )
    DEALLOCATE( var% isflag )
    DEALLOCATE( var% osnowd )
    DEALLOCATE( var% potev )
    DEALLOCATE( var% runoff )
    DEALLOCATE( var% rnof1 )
    DEALLOCATE( var% rnof2 )
    DEALLOCATE( var% rtsoil )
    DEALLOCATE( var% sconds )
    DEALLOCATE( var% sdepth )
    DEALLOCATE( var% smass )
    DEALLOCATE( var% snage )
    DEALLOCATE( var% snowd )
    DEALLOCATE( var% smelt )
    DEALLOCATE( var% ssdn )
    DEALLOCATE( var% ssdnn )
    DEALLOCATE( var% tgg )
    DEALLOCATE( var% tggsn )
    DEALLOCATE( var% tss )
    DEALLOCATE( var% tss_p )
    DEALLOCATE( var% deltss )
    DEALLOCATE( var% owb1 )
    DEALLOCATE( var% wb )
    DEALLOCATE( var% wbice )
    DEALLOCATE( var% wblf )
    DEALLOCATE( var%wbtot )
    DEALLOCATE( var%wbtot1 )
    DEALLOCATE( var%wbtot2 )
    DEALLOCATE( var%wb_lake )
    DEALLOCATE( var%totwblake )
    DEALLOCATE( var%sinfil )
    DEALLOCATE( var%evapfbl)
    DEALLOCATE( var%qstss)
    DEALLOCATE( var%wetfac )
    DEALLOCATE( var%owetfac )
    DEALLOCATE( var%t_snwlr )
    DEALLOCATE( var%wbfice )
    DEALLOCATE( var%tggav )
    DEALLOCATE( var%otgg )
    DEALLOCATE( var%otss )
    DEALLOCATE( var%otss_0 )
    DEALLOCATE( var%tprecip )
    DEALLOCATE( var%tevap )
    DEALLOCATE( var%trnoff )
    DEALLOCATE( var%totenbal )
    DEALLOCATE( var%totenbal2 )
    DEALLOCATE( var%fland )
    DEALLOCATE( var%ifland )
    DEALLOCATE( var%tilefrac )
    DEALLOCATE( var%qasrf )
    DEALLOCATE( var%qfsrf )
    DEALLOCATE( var%qssrf )
    !MD
    !Aquifer variables
    DEALLOCATE( var%GWwb )
    DEALLOCATE( var%GWhk )
    DEALLOCATE( var%GWdhkdw )
    DEALLOCATE( var%GWdsmpdw )
    DEALLOCATE( var%wtd )
    DEALLOCATE( var%GWsmp )
    DEALLOCATE( var%GWwbeq )
    DEALLOCATE( var%GWzq )
    DEALLOCATE( var%qhz )
    DEALLOCATE( var%qhlev )
    DEALLOCATE( var%satfrac )
    DEALLOCATE( var%Qrecharge )
    DEALLOCATE( var%rh_srf )
    DEALLOCATE( var%rtevap_unsat )
    DEALLOCATE( var%rtevap_sat )
    DEALLOCATE( var%rt_qh_sublayer )
    !soil moisture variables
    DEALLOCATE( var%wbeq )
    DEALLOCATE( var%zq )
    DEALLOCATE( var%icefrac )
    DEALLOCATE( var%fracice )
    DEALLOCATE( var%hk )
    DEALLOCATE( var%smp )
    DEALLOCATE( var%dhkdw )
    DEALLOCATE( var%dsmpdw )
    DEALLOCATE( var%wbliq )
    DEALLOCATE( var%wmliq )
    DEALLOCATE( var%wmice )
    DEALLOCATE( var%wmtot )

    !IF(cable_user%SOIL_STRUC=='sli') THEN
    DEALLOCATE ( var % S )
    DEALLOCATE ( var % Tsoil )
    DEALLOCATE ( var % SL )
    DEALLOCATE ( var % TL )
    DEALLOCATE ( var % h0)
    DEALLOCATE ( var % rex )
    DEALLOCATE ( var % wflux )
    DEALLOCATE ( var % delwcol )
    DEALLOCATE ( var % zdelta )
    DEALLOCATE ( var % kth )
    DEALLOCATE ( var % Tsurface )
    DEALLOCATE ( var % lE )
    DEALLOCATE ( var % evap )
    DEALLOCATE ( var % ciso )
    DEALLOCATE ( var % cisoL )
    DEALLOCATE ( var % rlitt )
    DEALLOCATE ( var % thetai )
    DEALLOCATE (var % snowliq)
    DEALLOCATE (var % nsteps)
    DEALLOCATE (var % nsnow)
    DEALLOCATE ( var % TsurfaceFR )
    DEALLOCATE ( var % Ta_daily )
    DEALLOCATE ( var % G0_daily )
    DEALLOCATE ( var % Qadv_daily )
    DEALLOCATE ( var % Qevap_daily )
    DEALLOCATE ( var % Qprec_daily )
    DEALLOCATE ( var % Qprec_snow_daily )

    ! END IF

  END SUBROUTINE dealloc_soil_snow_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_veg_parameter_type(var)

    TYPE(veg_parameter_type), INTENT(inout) :: var

    DEALLOCATE( var% canst1 )
    DEALLOCATE( var% dleaf )
    DEALLOCATE( var% ejmax )
    DEALLOCATE( var% iveg )
    DEALLOCATE( var% iLU )
    DEALLOCATE( var% meth )
    DEALLOCATE( var% frac4 )
    DEALLOCATE( var% hc )
    DEALLOCATE( var% vlai )
    DEALLOCATE( var% xalbnir )
    DEALLOCATE( var% rp20 )
    DEALLOCATE( var% rpcoef )
    DEALLOCATE( var% rs20 )
    DEALLOCATE( var% shelrb )
    DEALLOCATE( var% vegcf )
    DEALLOCATE( var% tminvj )
    DEALLOCATE( var% toptvj )
    DEALLOCATE( var% tmaxvj )
    DEALLOCATE( var% vbeta)
    DEALLOCATE( var% vcmax )
    DEALLOCATE( var% xfang )
    DEALLOCATE( var%extkn )
    DEALLOCATE( var%wai )
    DEALLOCATE( var%deciduous )
    DEALLOCATE( var%froot)
    DEALLOCATE( var%refl )
    DEALLOCATE( var%taul )
    DEALLOCATE( var%a1gs )
    DEALLOCATE( var%d0gs )
    DEALLOCATE( var%alpha )
    DEALLOCATE( var%convex )
    DEALLOCATE( var%cfrd )
    DEALLOCATE( var%gswmin )
    DEALLOCATE( var%conkc0 )
    DEALLOCATE( var%conko0 )
    DEALLOCATE( var%ekc )
    DEALLOCATE( var%eko )
    DEALLOCATE( var%g0 ) ! Ticket #56.
    DEALLOCATE( var%g1 ) ! Ticket #56.

    ! Deallocate variables for SLI soil model:
    !IF(cable_user%SOIL_STRUC=='sli') THEN
    DEALLOCATE ( var % rootbeta )
    DEALLOCATE ( var % gamma ) ! vh 20/07/09
    DEALLOCATE ( var % F10 )
    DEALLOCATE ( var % ZR )
    DEALLOCATE ( var % CLitt )
    DEALLOCATE ( var % disturbance_interval )
    DEALLOCATE ( var % disturbance_intensity )
    IF(ASSOCIATED(var % gamma)) DEALLOCATE ( var % gamma )
    ! END IF

  END SUBROUTINE dealloc_veg_parameter_type

  ! ------------------------------------------------------------------------------

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
    DEALLOCATE ( var % evapfbl )
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
    !DEALLOCATE( var% fescor_upp ) !SSEB variable
    !DEALLOCATE( var% fescor_low ) !SSEB variable
    DEALLOCATE( var% gswx )
    DEALLOCATE( var% oldcansto )
    DEALLOCATE( var% zetar )
    DEALLOCATE( var% zetash )
    DEALLOCATE ( var % fwsoil )
    DEALLOCATE ( var % ofes )
    DEALLOCATE( var% sublayer_dz )

    !! vh_js !! liiter resistances to heat and vapour transfer
    DEALLOCATE (var % kthLitt)
    DEALLOCATE (var % DvLitt)

  END SUBROUTINE dealloc_canopy_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_radiation_type(var)

    TYPE(radiation_type), INTENT(inout) :: var

    DEALLOCATE( var% albedo )
    DEALLOCATE( var% extkb )
    DEALLOCATE( var% extkd2 )
    DEALLOCATE( var% extkd )
    DEALLOCATE( var% flws )
    DEALLOCATE( var% fvlai )
    DEALLOCATE( var% latitude )
    DEALLOCATE( var% lwabv )
    DEALLOCATE( var% qcan )
    DEALLOCATE( var% qssabs )
    DEALLOCATE( var% rhocdf )
    DEALLOCATE( var% rniso )
    DEALLOCATE( var% scalex )
    DEALLOCATE( var% transd )
    DEALLOCATE( var% trad )
    DEALLOCATE( var% otrad )
    DEALLOCATE( var% reffdf )
    DEALLOCATE( var% reffbm )
    DEALLOCATE( var% extkbm )
    DEALLOCATE( var% extkdm )
    DEALLOCATE( var% fbeam )
    DEALLOCATE( var% cexpkbm )
    DEALLOCATE( var% cexpkdm )
    DEALLOCATE( var% rhocbm )
    DEALLOCATE( var% transb )
    DEALLOCATE( var% albedo_T )
    DEALLOCATE( var% gradis )
    DEALLOCATE( var% longitude )
    DEALLOCATE( var% workp1 )
    DEALLOCATE( var% workp2 )
    DEALLOCATE( var% workp3 )

  END SUBROUTINE dealloc_radiation_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_roughness_type(var)

    TYPE(roughness_type), INTENT(inout) :: var

    DEALLOCATE ( var % coexp )
    DEALLOCATE ( var % disp )
    DEALLOCATE ( var % hruff )
    DEALLOCATE ( var % hruff_grmx )
    DEALLOCATE ( var % rt0us )
    DEALLOCATE ( var % rt1usa )
    DEALLOCATE ( var % rt1usb )
    DEALLOCATE ( var % rt1 )
    DEALLOCATE ( var % term2 )
    DEALLOCATE ( var % term3 )
    DEALLOCATE ( var % term5 )
    DEALLOCATE ( var % term6 )
    DEALLOCATE ( var % term6a )
    DEALLOCATE ( var % usuh )
    DEALLOCATE ( var % za_uv )
    DEALLOCATE ( var % za_tq )
    DEALLOCATE ( var % z0m )
    DEALLOCATE ( var % zref_uv )
    DEALLOCATE ( var % zref_tq )
    DEALLOCATE ( var % zruffs )
    DEALLOCATE ( var % z0soilsn )
    DEALLOCATE ( var % z0soil )

  END SUBROUTINE dealloc_roughness_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_air_type(var)

    TYPE(air_type), INTENT(inout) :: var

    DEALLOCATE ( var % rho )
    DEALLOCATE ( var % volm )
    DEALLOCATE ( var % rlam )
    DEALLOCATE ( var % qsat )
    DEALLOCATE ( var % epsi )
    DEALLOCATE ( var % visc )
    DEALLOCATE ( var % psyc )
    DEALLOCATE ( var % dsatdk )
    DEALLOCATE ( var % cmolar )

  END SUBROUTINE dealloc_air_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_met_type(var)

    TYPE(met_type), INTENT(inout) :: var

    DEALLOCATE ( var % ca )
    DEALLOCATE ( var % year )
    DEALLOCATE ( var % moy )
    DEALLOCATE ( var % doy )
    DEALLOCATE ( var % hod )
    DEALLOCATE ( var % fsd )
    DEALLOCATE ( var % ofsd )
    DEALLOCATE ( var % fld )
    DEALLOCATE ( var % precip )
    DEALLOCATE ( var % precip_sn )
    DEALLOCATE ( var % tk )
    DEALLOCATE ( var % tvair )
    DEALLOCATE ( var % tvrad )
    DEALLOCATE ( var % pmb )
    DEALLOCATE ( var % ua )
    DEALLOCATE ( var % qv )
    DEALLOCATE ( var % qvair )
    DEALLOCATE ( var % da )
    DEALLOCATE ( var % dva )
    DEALLOCATE ( var % coszen )
    DEALLOCATE ( var % Ndep )
    DEALLOCATE ( var % Pdep )

  END SUBROUTINE dealloc_met_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_sum_flux_type(var)

    TYPE(sum_flux_type), INTENT(inout) :: var

    DEALLOCATE ( var % sumpn )
    DEALLOCATE ( var % sumrp )
    DEALLOCATE ( var % sumrpw )
    DEALLOCATE ( var % sumrpr )
    DEALLOCATE ( var % sumrs )
    DEALLOCATE ( var % sumrd )
    DEALLOCATE ( var % dsumpn )
    DEALLOCATE ( var % dsumrp )
    DEALLOCATE ( var % dsumrs )
    DEALLOCATE ( var % dsumrd )
    DEALLOCATE ( var % sumxrp )
    DEALLOCATE ( var % sumxrs )

  END SUBROUTINE dealloc_sum_flux_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_bgc_pool_type(var)

    TYPE(bgc_pool_type), INTENT(inout) :: var

    DEALLOCATE ( var % cplant )
    DEALLOCATE ( var % csoil )

  END SUBROUTINE dealloc_bgc_pool_type


END MODULE cable_def_types_mod
