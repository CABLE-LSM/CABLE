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

MODULE cable_def_types_mod

  ! Contains all variables which are not subroutine-internal

  IMPLICIT NONE

  SAVE

  PUBLIC

  !---CABLE default KINDs for representing INTEGER/REAL values
  !---at least 10-digit precision

  INTEGER :: mp, & ! # total no of patches/tiles
       mvtype,   & ! total # vegetation types,   from input
       mstype,   & ! total # soil types,         from input
       mland       ! # land grid cells

  INTEGER, PARAMETER :: &
       ! i_d  = KIND(9), &                  ! probably unintended
       i_d  = SELECTED_INT_KIND(9), &       ! as in pop.f90
       ! r_2  = SELECTED_REAL_KIND(12, 50), & ! old
       r_2  = KIND(1.0d0), &                ! as in pop.f90
       n_tiles = 17,  & ! # possible no of different
       ncp = 3,       & ! # vegetation carbon stores
       ncs = 2,       & ! # soil carbon stores
       mf = 2,        & ! # leaves (sunlit, shaded)
       nrb = 3,       & ! # radiation bands
       msn = 3,       & ! max # snow layers
       swb = 2,       & ! # shortwave bands
       niter = 4,     & ! number of iterations for za/L
                                ! ms = 12          ! # soil layers
       ms = 6           ! # soil layers - standard
  ! ms = 13          ! for Loetschental experiment
  ! ms = 10

  !   PRIVATE :: r_2, ms, msn, mf, nrb, ncp, ncs

  ! .............................................................................

  ! Energy and water balance variables:
  TYPE balances_type

     REAL, DIMENSION(:), POINTER :: &
          drybal => null(),           & ! energy balance for dry canopy
          ebal => null(),             & ! energy balance per time step (W/m^2)
          ebal_tot => null(),         & ! cumulative energy balance (W/m^2)
          ebal_cncheck => null(),     & ! energy balance consistency check (W/m^2)
          ebal_tot_cncheck => null(), & ! cumulative energy balance (W/m^2)
          ebaltr => null(),           & ! energy balance per time step (W/m^2)
          ebal_tottr => null(),       & ! cumulative energy balance (W/m^2)
          evap_tot => null(),         & ! cumulative evapotranspiration (mm/dels)
          osnowd0 => null(),          & ! snow depth, first time step
          precip_tot => null(),       & ! cumulative precipitation (mm/dels)
          rnoff_tot => null(),        & ! cumulative runoff (mm/dels)
          wbal => null(),             & ! water balance per time step (mm/dels)
          wbal_tot => null(),         & ! cumulative water balance (mm/dels)
          wbtot0 => null(),           & ! total soil water (mm), first time step
          wetbal => null(),           & ! energy balance for wet canopy
          cansto0 => null(),          & ! canopy water storage (mm)
          owbtot => null(),           & ! total soil water (mm), first time step
          evapc_tot => null(),        & ! cumulative evapotranspiration (mm/dels)
          evaps_tot => null(),        & ! cumulative evapotranspiration (mm/dels)
          rnof1_tot => null(),        & ! cumulative runoff (mm/dels)
          rnof2_tot => null(),        & ! cumulative runoff (mm/dels)
          snowdc_tot => null(),       & ! cumulative runoff (mm/dels)
          wbal_tot1 => null(),        & ! cumulative water balance (mm/dels)
          delwc_tot => null(),        & ! energy balance for wet canopy
          qasrf_tot => null(),        & ! heat advected to the snow by precip.
          qfsrf_tot => null(),        & ! energy of snowpack phase changes
          qssrf_tot => null(), &        ! energy of snowpack phase changes
          Radbal => null(), &
          EbalSoil => null(), &
          Ebalveg => null(), &
          Radbalsum => null()

  END TYPE balances_type

  ! .............................................................................

  ! Soil parameters:
  TYPE soil_parameter_type

     INTEGER, DIMENSION(:), POINTER :: &
          isoilm     ! integer soil type

     REAL, DIMENSION(:), POINTER :: &
          bch => null(),     & ! parameter b in Campbell equation
          c3 => null(),      & ! c3 drainage coeff (fraction)
          clay => null(),    & ! fraction of soil which is clay
          css => null(),     & ! soil specific heat capacity [kJ/kg/K]
          hsbh => null(),    & ! difsat * etasat (=hyds*abs(sucs)*bch)
          hyds => null(),    & ! hydraulic conductivity @ saturation [m/s], Ksat
          i2bp3 => null(),   & ! par. one in K vis suction (=nint(bch)+2)
          ibp2 => null(),    & ! par. two in K vis suction (fn of pbch)
          rhosoil => null(), & ! soil density [kg/m3]
          sand => null(),    & ! fraction of soil which is sand
          sfc => null(),     & ! vol H2O @ field capacity
          silt => null(),    & ! fraction of soil which is silt
          ssat => null(),    & ! vol H2O @ saturation
          sucs => null(),    & ! suction at saturation (m)
          swilt => null(),   & ! vol H2O @ wilting
          zse => null(),     & ! thickness of each soil layer (1=top) in m
          zshh => null(),    & ! distance between consecutive layer midpoints (m)
                                ! vars intro for Ticket #27
          soilcol => null(), & ! keep color for all patches/tiles
          albsoilf => null()   ! soil reflectance

     REAL(r_2), DIMENSION(:), POINTER :: &
          cnsd => null(),    & ! thermal conductivity of dry soil [W/m/K]
          pwb_min => null()    ! working variable (swilt/ssat)**ibp2

     REAL, DIMENSION(:,:), POINTER :: &
          albsoil    ! soil reflectance (2nd dim. BP 21Oct2009)

     ! Additional SLI parameters
     INTEGER,   DIMENSION(:),   POINTER :: nhorizons => null() ! number of soil horizons
     INTEGER,   DIMENSION(:,:), POINTER :: ishorizon => null() ! horizon number 1:nhorizons
     REAL(r_2), DIMENSION(:),   POINTER :: clitt => null()     ! litter (tC/ha)
     REAL(r_2), DIMENSION(:),   POINTER :: zeta => null()      ! macropore parameter
     REAL(r_2), DIMENSION(:),   POINTER :: fsatmax => null()   ! variably saturated area parameter
     REAL(r_2), DIMENSION(:,:), POINTER :: swilt_vec => null() ! vol H2O @ wilting
     REAL(r_2), DIMENSION(:,:), POINTER :: ssat_vec => null()  ! vol H2O @ sat
     REAL(r_2), DIMENSION(:,:), POINTER :: sfc_vec => null()   ! vol H2O @ fc

  END TYPE soil_parameter_type

  ! .............................................................................

  ! Soil and snow variables:
  TYPE soil_snow_type

     INTEGER, DIMENSION(:), POINTER :: isflag => null() ! 0 => no snow 1 => snow

     REAL, DIMENSION(:), POINTER :: &
          iantrct => null(), & ! pointer to Antarctic land points
          pudsto => null(),  & ! puddle storage
          pudsmx => null(),  & ! puddle storage
          cls => null(),     & ! factor for latent heat
          dfn_dtg => null(), & ! d(canopy%fns)/d(ssnow%tgg)
          dfh_dtg => null(), & ! d(canopy%fhs)/d(ssnow%tgg)
          dfe_ddq => null(), & ! d(canopy%fes)/d(dq)
          ddq_dtg => null(), & ! d(dq)/d(ssnow%tgg)
          evapsn => null(),  & ! snow evaporation
          fwtop => null(),   & ! water flux to the soil
          fwtop1 => null(),  & ! water flux to the soil
          fwtop2 => null(),  & ! water flux to the soil
          fwtop3 => null(),  & ! water flux to the soil
          osnowd => null(),  & ! snow depth from previous time step
          potev => null(),   & ! potential evapotranspiration
          runoff => null(),  & ! total runoff (mm/dels)
          rnof1 => null(),   & ! surface runoff (mm/dels)
          rnof2 => null(),   & ! deep drainage (mm/dels)
          rtsoil => null(),  & ! turbulent resistance for soil
          wbtot1 => null(),  & ! total soil water (mm)
          wbtot2 => null(),  & ! total soil water (mm)
          wb_lake => null(), &
          sinfil => null(),  &
          qstss => null(),   &
          wetfac => null(),  & ! surface wetness fact. at current time step
          owetfac => null(), & ! surface wetness fact. at previous time step
          t_snwlr => null(), & ! top snow layer depth in 3 layer snowpack
          tggav => null(),   & ! mean soil temperature in K
          otgg => null(),    & ! soil temperature in K
          otss => null(),    & ! surface temperature (weighted soil, snow)
          otss_0 => null(),  & ! surface temperature (weighted soil, snow)
          tprecip => null(), &
          tevap => null(),   &
          trnoff => null(),  &
          totenbal => null(),  & !
          totenbal2 => null(), &
          fland => null(),   & ! factor for latent heat
          ifland => null(),  & ! integer soil type
          qasrf => null(),   & ! heat advected to the snow by precip.
          qfsrf => null(),   & ! energy of snowpack phase changes
          qssrf => null(),   & ! sublimation
          snage => null(),   & ! snow age
          snowd => null(),   & ! snow depth (liquid water)
          smelt => null(),   & ! snow melt
          ssdnn => null(),   & ! average snow density
          tss => null(),     & ! surface temperature (weighted soil, snow)
          tss_p => null(),   & ! surface temperature (weighted soil, snow)
          deltss => null(),  & ! surface temperature (weighted soil, snow)
          owb1 => null()       ! surface temperature (weighted soil, snow)

     REAL, DIMENSION(:,:), POINTER :: &
          sconds => null(),     & !
          sdepth => null(),     & ! snow depth
          smass => null(),      & ! snow mass
          ssdn => null(),       & ! snow densities
          tgg => null(),        & ! soil temperature in K
          tggsn => null(),      & ! snow temperature in K
          dtmlt => null(),      & ! water flux to the soil
          albsoilsn => null(),  & ! soil + snow reflectance
          evapfbl => null(),    & !
          tilefrac => null()      ! factor for latent heat

     REAL(r_2), DIMENSION(:), POINTER :: &
          wbtot => null()   ! total soil water (mm)

     REAL(r_2), DIMENSION(:,:), POINTER :: &
          gammzz => null(),  & ! heat capacity for each soil layer
          wb => null(),      & ! volumetric soil moisture (solid+liq)
          wbice => null(),   & ! soil ice
          wblf => null(),    & !
          wbfice => null()     !

     ! Additional SLI variables:
     REAL(r_2), DIMENSION(:,:), POINTER :: S => null()         ! moisture content relative to sat value    (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil => null()     ! Tsoil (deg C)
     REAL(r_2), DIMENSION(:),   POINTER :: SL => null()        ! litter moisture content relative to sat value (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:),   POINTER :: TL => null()        ! litter temperature in K     (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:),   POINTER :: h0 => null()        ! pond height in m            (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: rex => null()       ! root extraction from each layer (mm/dels)
     REAL(r_2), DIMENSION(:,:), POINTER :: wflux => null()     ! water flux at layer boundaries (mm s-1)
     REAL(r_2), DIMENSION(:),   POINTER :: delwcol => null()   ! change in water column (mm / dels)
     REAL(r_2), DIMENSION(:),   POINTER :: zdelta => null()    ! water table depth           (edit vh 23/06/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: kth => null()       ! thermal conductivity           (edit vh 29/07/08)
     REAL(r_2), DIMENSION(:),   POINTER :: Tsurface => null()  !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:),   POINTER :: lE => null()      ! soil latent heat flux
     REAL(r_2), DIMENSION(:),   POINTER :: evap => null()    ! soil evaporation (mm / dels)
     REAL(r_2), DIMENSION(:,:), POINTER :: ciso => null()    ! concentration of minor isotopologue in soil water (kg m-3 water)
     REAL(r_2), DIMENSION(:),   POINTER :: cisoL => null()   ! concentration of minor isotopologue in litter water (kg m-3 water)
     REAL(r_2), DIMENSION(:),   POINTER :: rlitt => null()   ! resistance to heat/moisture transfer through litter (m-1 s)
     REAL(r_2), DIMENSION(:,:), POINTER :: thetai => null()  ! volumetric ice content (MC)
     REAL(r_2), DIMENSION(:,:), POINTER :: snowliq => null() ! liquid snow content (mm H2O)
     REAL(r_2), DIMENSION(:),   POINTER :: nsteps => null()  ! number of iterations at each timestep
     REAL(r_2), DIMENSION(:),   POINTER :: TsurfaceFR => null() ! temperature at surface (soil, pond or litter) (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: Ta_daily => null() ! air temp averaged over last 24h
     INTEGER, DIMENSION(:),     POINTER :: nsnow => null() ! number of layers in snow-pack (0-nsnow_max)
     REAL(r_2), DIMENSION(:),   POINTER :: Qadv_daily => null() ! advective heat flux into surface , daily average (W m-2)
     REAL(r_2), DIMENSION(:),   POINTER :: G0_daily => null()  ! conductive heat flux into surface , daily average (W m-2)
     REAL(r_2), DIMENSION(:),   POINTER :: Qevap_daily => null() ! evaporative flux at surface, daily average (m s-1)
     REAL(r_2), DIMENSION(:),   POINTER :: Qprec_daily => null() ! liquid precip, daily average (m s-1)
     REAL(r_2), DIMENSION(:),   POINTER :: Qprec_snow_daily => null() ! solid precip, daily average (m s-1)
     real(r_2), dimension(:),   pointer :: E_fusion_sn => null()
     real(r_2), dimension(:),   pointer :: E_sublimation_sn => null()
     real(r_2), dimension(:),   pointer :: latent_heat_sn => null()
     real(r_2), dimension(:),   pointer :: evap_liq_sn => null()
     real(r_2), dimension(:),   pointer :: surface_melt => null()
     real(r_2), dimension(:),   pointer :: Qadv_rain_sn => null()

  END TYPE soil_snow_type

  ! .............................................................................

  ! Vegetation parameters:
  TYPE veg_parameter_type

     INTEGER, DIMENSION(:), POINTER :: &
          iveg => null(),   & ! vegetation type
          ivegp => null(),  & ! dominant potential vegetation type
          iLU => null()        ! land use type

     REAL, DIMENSION(:), POINTER :: &
          canst1 => null(),  & ! max intercepted water by canopy (mm/LAI)
          dleaf => null(),   & ! chararacteristc legnth of leaf (m)
          ejmax => null(),   & ! max pot. electron transp rate top leaf(mol/m2/s)
          ejmax_shade => null(), & ! max pot. electron transp rate top leaf(mol/m2/s)
          ejmax_sun => null(),   & ! max pot. electron transp rate top leaf(mol/m2/s)
          meth => null(),    & ! method for calculation of canopy fluxes and temp.
          frac4 => null(),   & ! fraction of c4 plants
          hc => null(),      & ! roughness height of canopy (veg - snow)
          vlai => null(),    & ! leaf area index
          xalbnir => null(), &
          rp20 => null(),    & ! plant respiration coefficient at 20 C
          rpcoef => null(),  & ! temperature coef nonleaf plant respiration (1/C)
          rs20 => null(),    & ! soil respiration at 20 C [mol m-2 s-1]
          shelrb => null(),  & ! sheltering factor (dimensionless)
          vegcf => null(),   & ! kdcorbin, 08/10
          tminvj => null(),  & ! min temperature of the start of photosynthesis
          toptvj => null(),  & ! opt temperature of the start of photosynthesis
          tmaxvj => null(),  & ! max temperature of the start of photosynthesis
          vbeta => null(),   & !
          vcmax => null(),   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
          vcmax_shade => null(), & ! max RuBP carboxylation rate top leaf (mol/m2/s)
          vcmax_sun => null(),   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
          xfang => null(),   & ! leaf angle PARAMETER
          extkn => null(),   & ! extinction coef for vertical
          vlaimax => null(), & ! extinction coef for vertical
          wai => null(),     & ! wood area index (stem+branches+twigs)
          a1gs => null(),    & ! a1 parameter in stomatal conductance model
          d0gs => null(),    & ! d0 in stomatal conductance model
          alpha => null(),   & ! initial slope of J-Q response curve
          convex => null(),  & ! convexity of J-Q response curve
          cfrd => null(),    & ! ratio of day respiration to vcmax
          gswmin => null(),  & ! minimal stomatal conductance
          conkc0 => null(),  & ! Michaelis-menton constant for carboxylase
          conko0 => null(),  & ! Michaelis-menton constant for oxygenase
          ekc => null(),     & ! activation energy for caroxylagse
          eko => null(),     & ! acvtivation enegery for oxygenase
          g0 => null(),      & ! Belinda's stomatal model intercept, Ticket #56.
          g1 => null(),      & ! Belinda's stomatal model slope, Ticket #56.
          vcmaxcc => null(), & ! max Cc-based carboxylation rate top leaf (mol/m2/s)
          ejmaxcc => null(), & ! max Cc-based RuBP regeneration rate top leaf (mol/m2/s)
          gmmax => null(),   & ! max. mesophyll conductance at 25degC top leaf (mol/2/s)
          gm => null(),      & ! mesophyll conductance adjusted for N content (mol/2/s)
          c4kci => null(),   & ! C4 plants: initial slope of An-Ci response curve (Ci-based)
          c4kcc => null(),   & ! C4 plants: initial slope of An-Ci response curve (Cc-based)
          bjv => null()        ! Jmax-Vcmax ratio at 25degC

     LOGICAL, DIMENSION(:), POINTER :: &
          deciduous => null() ! flag used for phenology fix

     REAL, DIMENSION(:,:), POINTER :: &
          refl => null(),    &
          taul => null(),    &
          froot => null()      ! fraction of root in each soil layer

     ! Additional  veg parameters:
     REAL(r_2), DIMENSION(:), POINTER :: rootbeta => null() ! parameter for estimating vertical root mass distribution (froot)
     REAL(r_2), DIMENSION(:), POINTER :: gamma => null()    ! parameter in root efficiency function (Lai and Katul 2000)
     REAL(r_2), DIMENSION(:), POINTER :: ZR => null()       ! maximum rooting depth (cm)
     REAL(r_2), DIMENSION(:), POINTER :: F10 => null()      ! fraction of roots in top 10 cm

     REAL(r_2), DIMENSION(:), POINTER :: clitt => null()     !

     ! Additional POP veg param
     INTEGER,   DIMENSION(:,:), POINTER :: disturbance_interval => null()
     REAL(r_2), DIMENSION(:,:), POINTER :: disturbance_intensity => null()

  END TYPE veg_parameter_type

  ! .............................................................................

  ! Canopy/vegetation variables:
  TYPE canopy_type

     REAL, DIMENSION(:), POINTER :: &
          cansto => null(),  & ! canopy water storage (mm)
          cduv => null(),    & ! drag coefficient for momentum
          delwc => null(),   & ! change in canopy water store (mm/dels)
          dewmm => null(),   & ! dewfall (mm)
          fe => null(),      & ! total latent heat (W/m2)
          fh => null(),      & ! total sensible heat (W/m2)
          fpn => null(),     & ! plant photosynthesis (g C m-2 s-1)
          frp => null(),     & ! plant respiration (g C m-2 s-1)
          frpw => null(),    & ! plant respiration (woody component) (g C m-2 s-1)
          frpr => null(),    & ! plant respiration (root component) (g C m-2 s-1)
          frs => null(),     & ! soil respiration (g C m-2 s-1)
          fnee => null(),    & ! net carbon flux (g C m-2 s-1)
          frday => null(),   & ! daytime leaf resp
          fnv => null(),     & ! net rad. avail. to canopy (W/m2)
          fev => null(),     & ! latent hf from canopy (W/m2)
          epot => null(),    & ! total potential evaporation
          fnpp => null(),    & ! npp flux
          fevw_pot => null(),& ! potential lat heat from canopy
          gswx_T => null(),  & ! stom cond for water
          cdtq => null(),    & ! drag coefficient for momentum
          wetfac_cs => null(),&!
          fevw => null(),    & ! lat heat fl wet canopy (W/m2)
          fhvw => null(),    & ! sens heatfl from wet canopy (W/m2)
          oldcansto => null(),&! canopy water storage (mm)
          fhv => null(),     & ! sens heatfl from canopy (W/m2)
          fns => null(),     & ! net rad avail to soil (W/m2)
          fhs => null(),     & ! sensible heat flux from soil
          fhs_cor => null(), &
          ga => null(),      & ! ground heat flux (W/m2) ???
          ghflux => null(),  & ! ground heat flux (W/m2) ???
          precis => null(),  & ! throughfall to soil, after snow (mm)
          qscrn => null(),   & ! specific humudity at screen height (g/g)
          rnet => null(),    & ! net radiation absorbed by surface (W/m2)
          rniso => null(),   & ! isothermal net radiation absorbed by surface (W/m2)
          segg => null(),    & ! latent heatfl from soil mm
          sghflux => null(), & ! ground heat flux (W/m2) ???
          through => null(), & ! canopy throughfall (mm)
          spill => null(),   & ! can.storage excess after dewfall (mm)
          tscrn => null(),   & ! air temperature at screen height (oC)
          wcint => null(),   & ! canopy rainfall interception (mm)
          tv => null(),      & ! vegetation temp (K)
          us => null(),      & ! friction velocity
          uscrn => null(),   & ! wind speed at screen height (m/s)
          vlaiw => null(),   & ! lai adj for snow depth for calc of resistances
          rghlai => null(),  & ! lai adj for snow depth for calc of resistances
          fwet => null()       ! fraction of canopy wet

     REAL, DIMENSION(:,:), POINTER :: &
          evapfbl => null(), &
          gswx => null(),    & ! stom cond for water
          zetar => null(),   & ! stability parameter (ref height)
                                !! vh_js !!
          zetash => null()     ! stability parameter (shear height)

     REAL(r_2), DIMENSION(:), POINTER :: &
          fess => null(),    & ! latent heatfl from soil (W/m2)
          fesp => null(),    & ! latent heatfl from soil (W/m2)
          dgdtg => null(),   & ! derivative of gflux wrt soil temp
          fes => null(),     & ! latent heatfl from soil (W/m2)
          fes_cor => null(), & ! latent heatfl from soil (W/m2)
          fevc => null(),    & ! dry canopy transpiration (W/m2)
          ofes => null(),    & ! latent heatfl from soil (W/m2)
          A_sh => null(),    & ! net photosynthesis from shaded leaves
          A_sl => null(),    & ! net photosynthesis from sunlit leaves
          A_slC => null(),   & ! net photosynthesis from sunlit leaves (rubisco limited)
          A_shC => null(),   & ! net photosynthesis from shaded leaves  (rubisco limited)
          A_slJ => null(),   & ! net photosynthesis from sunlit leaves (rubp limited)
          A_shJ => null(),   & ! net photosynthesis from shaded leaves  (rubp limited)
          GPP_sh => null(), &  ! gross photosynthesis from shaded leaves
          GPP_sl => null(), &  ! gross photosynthesis from sunlit leaves
          fevc_sl => null(), &  ! dry canopy transpiration sunlit leaves (W/m2)
          fevc_sh => null(), &  ! dry canopy transpiration shaded leaves (W/m2)
          eta_A_cs => null(),& ! elasticity of net photosynthesis wrt cs, mulitplied by net  photosythesis
          dAdcs => null(),   & ! sensitivity of net photosynthesis wrt cs, mulitplied by net photosythesis
          eta_GPP_cs => null(), &      ! elasticity of gross photosynthesis wrt cs, mulitplied by net photosythesis
          eta_A_cs_sl => null(), &      ! elasticity of net photosynthesis wrt cs, sl leaves
          eta_A_cs_sh => null(), &      ! elasticity of net photosynthesis wrt cs, sh leaves
          eta_fevc_cs_sl => null(), &      ! elasticity of net photosynthesis wrt cs, sl leaves
          eta_fevc_cs_sh => null(), &      ! elasticity of net photosynthesis wrt cs, sh leaves
          eta_fevc_cs => null(), & ! elasticity of transpiration wrt cs, mulitplied by transpiration
          cs => null(),      & ! leaf surface CO2 (ppm), mulitplied by net photosythesis
          cs_sl => null(),   & ! leaf surface CO2 (ppm) (sunlit)
          cs_sh => null(),   & ! leaf surface CO2 (ppm) (shaded)
                                ! ci_sl => null(), &    !  leaf internal CO2 (ppm) (sunlit)
                                ! ci_sh => null(), &    !  leaf internal CO2 (ppm) (shaded)
          tlf => null(),     & ! dry leaf temperature
          dlf => null()        ! dryleaf vp minus in-canopy vp (Pa)

     ! Additional variables:
     REAL(r_2), DIMENSION(:,:),   POINTER :: gw => null()     ! dry canopy conductance (ms-1) edit vh 6/7/09
     REAL(r_2), DIMENSION(:,:,:), POINTER :: ancj => null()   ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
     REAL(r_2), DIMENSION(:,:),   POINTER :: tlfy => null()   ! sunlit and shaded leaf temperatures
     REAL(r_2), DIMENSION(:,:),   POINTER :: ecy => null()    ! sunlit and shaded leaf transpiration (dry canopy)
     REAL(r_2), DIMENSION(:,:),   POINTER :: ecx => null()    ! sunlit and shaded leaf latent heat flux
     ! REAL(r_2), DIMENSION(:,:,:), POINTER :: ci => null()     ! intra-cellular CO2 vh 6/7/09
     REAL(r_2), DIMENSION(:),     POINTER :: fwsoil => null() !

     ! vh_js - litter thermal conductivity (Wm-2K-1) and vapour diffusivity (m2s-1)
     REAL(r_2), DIMENSION(:), POINTER :: kthLitt => null()
     REAL(r_2), DIMENSION(:), POINTER :: DvLitt => null()

     ! 13C
     real(r_2), dimension(:,:), pointer :: An => null()        ! sunlit and shaded net assimilation [mol(CO2)/m2/s]
     real(r_2), dimension(:,:), pointer :: Rd => null()        ! sunlit and shaded leaf respiration [mol(CO2)/m2/s]
     logical,   dimension(:),   pointer :: isc3 => null()      ! C3 / C4 mask
     real(r_2), dimension(:,:), pointer :: vcmax => null()     ! max RuBP carboxylation rate [mol(CO2)/m2/s]
     ! CO2 compensation point in absence of dark respiration [mol(CO2)/mol(air)]
     real(r_2), dimension(:,:), pointer :: gammastar => null()
     real(r_2), dimension(:,:), pointer :: gsc => null()       ! stomatal conductance for CO2 [mol(CO2)/m^2/s]
     real(r_2), dimension(:,:), pointer :: gbc => null()       ! leaf boundary layer conductance for CO2 [mol(CO2)/m^2/s]
     real(r_2), dimension(:,:), pointer :: gac => null()       ! aerodynamic conductance for CO2 [mol(CO2)/m^2/s]
     real(r_2), dimension(:,:), pointer :: ci => null()        ! stomatal CO2 concentration [mol(CO2)/mol(air)]

  END TYPE canopy_type

  ! .............................................................................

  ! Radiation variables:
  TYPE radiation_type

     REAL, DIMENSION(:), POINTER :: &
          transb => null(),  & ! fraction SW beam tranmitted through canopy
          albedo_T => null(),& ! canopy+soil albedo for VIS+NIR
          longitude => null(),&! longitude
          workp1 => null(),  & ! absorbed short-wave radiation for soil
          workp2 => null(),  & ! absorbed short-wave radiation for soil
          workp3 => null(),  & ! absorbed short-wave radiation for soil
          extkb => null(),   & ! beam radiation extinction coeff
          extkd2 => null(),  & ! diffuse 2D radiation extinction coeff
          extkd => null(),   & ! diffuse radiation extinction coeff (-)
          flws => null(),    & ! soil long-wave radiation
          latitude => null(),& ! latitude
          lwabv => null(),   & ! long wave absorbed by vegetation
          qssabs => null(),  & ! absorbed short-wave radiation for soil
          transd => null(),  & ! frac SW diffuse transmitted through canopy
          trad => null()       !  radiative temperature (soil and veg)

     REAL, DIMENSION(:,:), POINTER :: &
          fvlai => null(),   & ! leaf area index of big leaf
          rhocdf => null(),  & ! canopy diffuse reflectance (-)
          rniso => null(),   & ! sum(rad%qcan, 3) total abs by canopy (W/m2)
          scalex => null(),  & ! scaling PARAMETER for big leaf
          albedo => null(),  & ! canopy+soil albedo
          reffdf => null(),  & ! effective conopy diffuse reflectance
          reffbm => null(),  & ! effective conopy beam reflectance
          extkbm => null(),  & ! modified k beam(6.20)(for leaf scattering)
          extkdm => null(),  & ! modified k diffuse(6.20)(for leaf scattering)
          fbeam => null(),   & ! beam fraction
          cexpkbm => null(), & ! canopy beam transmittance
          cexpkdm => null(), & ! canopy diffuse transmittance
          rhocbm => null(),  & ! modified canopy beam reflectance(6.21)
          gradis => null()     ! radiative conductance

     REAL, DIMENSION(:,:,:), POINTER :: &
          qcan => null() ! absorbed radiation for canopy (W/m^2)

  END TYPE radiation_type

  ! .............................................................................

  ! Roughness variables:
  TYPE roughness_type

     REAL, DIMENSION(:), POINTER :: &
          disp => null(),    & ! zero-plane displacement
          hruff => null(),   & ! canopy height above snow level
          hruff_grmx => null(),&! max ht of canopy from tiles on same grid
          rt0us => null(),   & ! eq. 3.54, SCAM manual (CSIRO tech report 132)
          rt1usa => null(),  & ! resistance from disp to hruf
          rt1usb => null(),  & ! resist fr hruf to zruffs (zref if zref<zruffs)
          rt1 => null(),     & ! 1/aerodynamic conductance
          za_uv => null(),   & ! level of lowest atmospheric model layer
          za_tq => null(),   & ! level of lowest atmospheric model layer
          z0m => null(),     & ! roughness length
          zref_uv => null(), & ! Reference height for met forcing
          zref_tq => null(), & ! Reference height for met forcing
          zruffs => null(),  & ! SCALAR Roughness sublayer depth (ground=origin)
          z0soilsn => null(),& ! roughness length of bare soil surface
          z0soil => null()     ! roughness length of bare soil surface

     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL, DIMENSION(:), POINTER :: &
          coexp => null() ! Extinction coef for wind profile in canopy

     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL, DIMENSION(:), POINTER :: &
          usuh => null() ! Friction velocity/windspeed at canopy height

     REAL, DIMENSION(:), POINTER :: & ! for aerodyn resist. calc.
          term2 => null(), term3 => null(), term5 => null(), term6 => null(), term6a => null()

  END TYPE roughness_type

  ! .............................................................................

  ! Air variables:
  TYPE air_type

     REAL, DIMENSION(:), POINTER :: &
          rho => null(),     & ! dry air density (kg m-3)
          volm => null(),    & ! molar volume (m3 mol-1)
          rlam => null(),    & ! latent heat for water (j/kg)
          qsat => null(),    & ! saturation specific humidity
          epsi => null(),    & ! d(qsat)/dT ((kg/kg)/K)
          visc => null(),    & ! air kinematic viscosity (m2/s)
          psyc => null(),    & ! psychrometric constant
          dsatdk => null(),  & ! d(es)/dT (mb/K)
          cmolar => null()     ! conv. from m/s to mol/m2/s

  END TYPE air_type

  ! .............................................................................

  ! Meterological data:
  TYPE met_type

     INTEGER, DIMENSION(:), POINTER :: &
          year => null(),    & ! local time year AD
          moy => null()        ! local time month of year

     REAL, DIMENSION(:), POINTER :: &
          ca => null(),      & ! CO2 concentration (mol/mol)
          doy => null(),     & ! local time day of year = days since 0 hr 1st Jan
          hod => null(),     & ! local hour of day
          ofsd => null(),    & ! downward short-wave radiation (W/m2)
          fld => null(),     & ! downward long-wave radiation (W/m2)
          precip => null(),  & ! rainfall (liquid+solid)(mm/dels)
          precip_sn => null(),&! solid preipitation only (mm/dels)
          tk => null(),      & ! surface air temperature (oK)
          tvair => null(),   & ! within canopy air temperature (oK)
          tvrad => null(),   & ! radiative vegetation temperature (K)
          pmb => null(),     & ! surface air pressure (mbar)
          ua => null(),      & ! surface wind speed (m/s)
          qv => null(),      & ! surface specific humidity (g/g)
          qvair => null(),   & ! within canopy specific humidity (g/g)
          da => null(),      & ! water vap pressure deficit at ref height (Pa)
          dva => null(),     & ! in canopy water vap pressure deficit (Pa)
          coszen => null(),  & ! cos(zenith angle of sun)
          Ndep => null(),    & ! nitrogen deposition (gN m-2 d-1)
          Pdep => null(),    & ! P deposition (gP m-2 d-1)
          u10 => null(),     & ! 10 m horizontal wind (m/s)
          rhum => null()       ! relative humidity (%)
     REAL, DIMENSION(:,:), POINTER :: &
          fsd => null()  ! downward short-wave radiation (W/m2)

  END TYPE met_type

  ! .............................................................................

  ! Climate data:
  TYPE climate_type

     INTEGER :: nyear_average = 20
     INTEGER :: nday_average  = 31
     !      INTEGER, POINTER ::                                                  &
     INTEGER :: &
          nyears, & ! number of years in climate record
          doy ! day of year

     INTEGER, DIMENSION(:), POINTER :: &
          chilldays => null(), &  ! length of chilling period (period with T<5deg)
          iveg => null(), &       ! potential vegetation type based on climatic constraints
          biome => null(), &
          GMD => null(), &        ! growing moisture days (== number days since min moisture threshold)
          modis_igbp => null(), & ! IGBP biome classification
          DSLR => null(), &       ! days since last rain
          NDAY_Nesterov => null()

     REAL, DIMENSION(:), POINTER :: &
          dtemp => null(),        &                ! daily mean temperature
          dmoist => null(),        &               ! daily moisture availability
          dmoist_min => null(),        &           ! minimum daily moisture availability over the year
          dmoist_min20 => null(),        &         ! min daily moisture avail over the year, averaged over 20 y
          dmoist_max => null(),        &           ! maximum daily moisture availability over the year
          dmoist_max20 => null(),        &         ! max daily moisture avail over the year, averaged over 20 y
          mtemp => null(),       &                 ! mean temperature over the last 31 days
          qtemp => null(),       &                 ! mean temperature over the last 91 days
          mmoist => null(),        &               ! monthly moisture availability
          mtemp_min => null(),   &                 ! minimum monthly temperature
          mtemp_max => null(),   &                 ! maximum monhtly temperature
          qtemp_max => null(),   &                 ! mean temperature of the warmest quarter (so far this year)
          qtemp_max_last_year => null(),   &       ! mean temperature of the warmest quarter (last calendar year)
          mtemp_min20 => null(),   &               ! minimum monthly temperature, averaged over 20 y
          mtemp_max20 => null(),   &               ! maximum monhtly temperature, averaged over 20 y
          atemp_mean => null(),  &                 ! annual average temperature
          AGDD5 => null(),       &
          GDD5 => null(),        &                 ! growing degree day sum relative to 5deg base temperature
          AGDD0 => null(),        &                !
          GDD0 => null(),        &                 ! growing degree day sum relative to 0deg base temperature
          alpha_PT => null(),    &                 ! ratio of annual evap to annual PT evap
          evap_PT => null(),    &                  ! annual PT evap [mm]
          aevap => null(), &                       ! annual evap [mm]
          alpha_PT20 => null(), &
          GDD0_rec => null(), &                    ! growing degree day sum related to spring photosynthetic recovery
          frec => null(), &                        ! fractional photosynthetic recovery
          dtemp_min => null(), &                   ! daily minimum temperature
          fdorm => null(), &                       ! dormancy fraction (1 prior to first autumn frost; 0 after 10 severe frosts)
          fapar_ann_max => null(), &               ! maximum midday fpar so far this year
          fapar_ann_max_last_year => null(), &     ! maximum midday fpar last year
          AvgAnnMaxFAPAR => null(), &              ! average annual maximum FAPAR
          dtemp_max => null(), &                   ! daily maximum temperature
          drhum => null(),      &                  ! daily average relative humidity
          du10_max => null(),   &                  ! daily max wind speed at 10 m
          dprecip => null(),   &                   ! daily total precip (mm)
          aprecip => null(),   &                   ! total precip accumulated over the current year (mm)
          aprecip_av20 => null(), &                ! annual precip averaged over the last 20 y
          last_precip => null(), &                 ! rainfall accumulated since last day without rain
          KBDI => null(),  &                       ! Keetch-Byram-Drought-Index (Keetch, 1968)
          FFDI => null(), &                        ! Forest Fire Danger Index
          D_MacArthur => null(), &                 ! MacArthur Drought Factor
          Nesterov_current => null(), &            ! current nesterov index
          Nesterov_ann_max => null(), &            ! annual maximum nesterov index (current year)
          Nesterov_ann_max_last_year  => null(), & ! annual maximum nesterov index (last year)
          Nesterov_ann_running_max => null()

     REAL, DIMENSION(:,:), POINTER :: &
          mtemp_min_20 => null(), &    ! mimimum monthly temperatures for the last 20 y
          mtemp_max_20 => null(), &    ! maximum monthly temperatures for the last 20 y
          dmoist_min_20 => null(), &   ! min daily moisture for the last 20 y
          dmoist_max_20 => null(), &   ! max daily moisture for the last 20 y
          dtemp_31  => null(), &       ! daily temperature for the last 31 days
          dmoist_31  => null(), &      ! daily moisture availability for the last 31 days
          alpha_PT_20 => null(), &     ! priestley Taylor Coefft for last 20 y
          dtemp_91 => null(), &        ! daily temperature for the last 91 days
          APAR_leaf_sun => null(), &   ! flux of PAR absrobed by sunlit leaves (subdiurnal time-step)
          APAR_leaf_shade => null(), & ! flux of PAR absrobed by shaded leaves (subdiurnal time-step)
          Dleaf_sun => null(), &       ! leaf-air vapour pressure difference (sun leaves, subdiurnal time-step)
          Dleaf_shade => null(),  &    ! leaf-air vapour pressure difference (shade leaves, subdiurnal time-step)
          Tleaf_sun => null(), &       ! leaf T (sun leaves, subdiurnal time-step)
          Tleaf_shade => null(), &     ! leaf T  (shade leaves, subdiurnal time-step)
          cs_sun => null(), &          ! sun leaf cs (ppm CO2)
          cs_shade => null(), &        ! shade leaf cs (ppm CO2)
          scalex_sun => null(), &      ! canopy depth scaling factor on vcmax and jmax (sun leaves)
          scalex_shade => null(), &    ! canopy depth scaling factor on vcmax and jmax (shade leaves)
          fwsoil => null(), &          ! soil-moisture modifier to stomatal conductance
          aprecip_20 => null(), &      ! annual average rainfall for the last 20 years
          Rd_sun => null(), &
          Rd_shade => null()
  END TYPE climate_type

  ! .............................................................................

  ! Cumulative flux variables:
  TYPE sum_flux_type

     REAL, DIMENSION(:), POINTER :: &
          sumpn => null(),   & ! sum of canopy photosynthesis (g C m-2)
          sumrp => null(),   & ! sum of plant respiration (g C m-2)
          sumrpw => null(),  & ! sum of plant respiration (g C m-2)
          sumrpr => null(),  & ! sum of plant respiration (g C m-2)
          sumrs => null(),   & ! sum of soil respiration (g C m-2)
          sumrd => null(),   & ! sum of daytime respiration (g C m-2)
          dsumpn => null(),  & ! daily sumpn
          dsumrp => null(),  & ! daily sumrp
          dsumrs => null(),  & ! daily sumrs
          dsumrd => null(),  & ! daily sumrd
          sumxrp => null(),  & ! sum plant resp. modifier
          sumxrs => null()     ! sum soil resp. modifier

  END TYPE sum_flux_type

  ! .............................................................................

  TYPE bgc_pool_type

     REAL, DIMENSION(:,:), POINTER :: &
          cplant => null(),  & ! plant carbon (g C/m2))
          csoil => null()      ! soil carbon (g C/m2)

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
     MODULE PROCEDURE alloc_balances_type, &
          alloc_soil_parameter_type,         &
          alloc_soil_snow_type,              &
          alloc_veg_parameter_type,          &
          alloc_canopy_type,                 &
          alloc_radiation_type,              &
          alloc_roughness_type,              &
          alloc_air_type,                    &
          alloc_met_type,                    &
          alloc_sum_flux_type,               &
          alloc_bgc_pool_type,               &
          alloc_climate_type
  END INTERFACE alloc_cbm_var

  INTERFACE dealloc_cbm_var
     MODULE PROCEDURE dealloc_balances_type, &
          dealloc_soil_parameter_type,         &
          dealloc_soil_snow_type,              &
          dealloc_veg_parameter_type,          &
          dealloc_canopy_type,                 &
          dealloc_radiation_type,              &
          dealloc_roughness_type,              &
          dealloc_air_type,                    &
          dealloc_met_type,                    &
          dealloc_sum_flux_type,               &
          dealloc_bgc_pool_type
  END INTERFACE dealloc_cbm_var

  public :: zero_cbm_var
  interface zero_cbm_var
     module procedure zero_balances_type, &
          zero_soil_parameter_type,         &
          zero_soil_snow_type,              &
          zero_veg_parameter_type,          &
          zero_canopy_type,                 &
          zero_radiation_type,              &
          zero_roughness_type,              &
          zero_air_type,                    &
          zero_met_type,                    &
          zero_sum_flux_type,               &
          zero_bgc_pool_type,               &
          zero_climate_type
  end interface zero_cbm_var

  public :: print_cbm_var
  interface print_cbm_var
     module procedure print_balances_type, &
          print_soil_parameter_type,         &
          print_soil_snow_type,              &
          print_veg_parameter_type,          &
          print_canopy_type,                 &
          print_radiation_type,              &
          print_roughness_type,              &
          print_air_type,                    &
          print_met_type,                    &
          print_sum_flux_type,               &
          print_bgc_pool_type,               &
          print_climate_type
  end interface print_cbm_var

CONTAINS

  SUBROUTINE alloc_balances_type(var, mp)

    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%drybal(mp))
    allocate(var%ebal(mp))
    allocate(var%ebal_tot(mp))
    allocate(var%ebaltr(mp))
    allocate(var%ebal_tottr(mp))
    allocate(var%ebal_cncheck(mp))
    allocate(var%ebal_tot_cncheck(mp))
    allocate(var%evap_tot(mp))
    allocate(var%osnowd0(mp))
    allocate(var%precip_tot(mp))
    allocate(var%rnoff_tot(mp))
    allocate(var%wbal(mp))
    allocate(var%wbal_tot(mp))
    allocate(var%wbtot0(mp))
    allocate(var%wetbal(mp))
    allocate(var%cansto0(mp))
    allocate(var%evapc_tot(mp))
    allocate(var%evaps_tot(mp))
    allocate(var%rnof1_tot(mp))
    allocate(var%rnof2_tot(mp))
    allocate(var%snowdc_tot(mp))
    allocate(var%wbal_tot1(mp))
    allocate(var%owbtot(mp))
    allocate(var%delwc_tot(mp))
    allocate(var%qasrf_tot(mp))
    allocate(var%qfsrf_tot(mp))
    allocate(var%qssrf_tot(mp))

    allocate(var%Radbal(mp))
    allocate(var%EbalSoil(mp))
    allocate(var%Ebalveg(mp))
    allocate(var%Radbalsum(mp))

  END SUBROUTINE alloc_balances_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_soil_parameter_type(var, mp)

    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%bch(mp))
    allocate(var%c3(mp))
    allocate(var%clay(mp))
    allocate(var%css(mp))
    allocate(var%hsbh(mp))
    allocate(var%hyds(mp))
    allocate(var%i2bp3(mp))
    allocate(var%ibp2(mp))
    allocate(var%isoilm(mp))
    allocate(var%rhosoil(mp))
    allocate(var%sand(mp))
    allocate(var%sfc(mp))
    allocate(var%silt(mp))
    allocate(var%ssat(mp))
    allocate(var%sucs(mp))
    allocate(var%swilt(mp))
    allocate(var%zse(ms))
    allocate(var%zshh(ms+1))
    allocate(var%cnsd(mp))
    allocate(var%albsoil(mp, nrb))
    allocate(var%pwb_min(mp))
    allocate(var%albsoilf(mp))
    allocate(var%soilcol(mp))

    ! Allocate variables for SLI soil model:
    allocate(var%nhorizons(mp))
    allocate(var%ishorizon(mp,ms))
    allocate(var%clitt(mp))
    allocate(var%zeta(mp))
    allocate(var%fsatmax(mp))
    allocate(var%swilt_vec(mp,ms))
    allocate(var%ssat_vec(mp,ms))
    allocate(var%sfc_vec(mp,ms))
    if (.not.(associated(var%swilt_vec))) allocate(var%swilt_vec(mp,ms))
    if (.not.(associated(var%ssat_vec)))  allocate(var%ssat_vec(mp,ms))
    if (.not.(associated(var%sfc_vec)))   allocate(var%sfc_vec(mp,ms))

  END SUBROUTINE alloc_soil_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_soil_snow_type(var, mp)

    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%iantrct(mp))
    allocate(var%pudsto(mp))
    allocate(var%pudsmx(mp))
    allocate(var%dtmlt(mp,3))
    allocate(var%albsoilsn(mp,nrb))
    allocate(var%cls(mp))
    allocate(var%dfn_dtg(mp))
    allocate(var%dfh_dtg(mp))
    allocate(var%dfe_ddq(mp))
    allocate(var%ddq_dtg(mp))
    allocate(var%evapsn(mp))
    allocate(var%fwtop(mp))
    allocate(var%fwtop1(mp))
    allocate(var%fwtop2(mp))
    allocate(var%fwtop3(mp))
    allocate(var%gammzz(mp,ms))
    allocate(var%isflag(mp))
    allocate(var%osnowd(mp))
    allocate(var%potev(mp))
    allocate(var%runoff(mp))
    allocate(var%rnof1(mp))
    allocate(var%rnof2(mp))
    allocate(var%rtsoil(mp))
    allocate(var%sconds(mp,msn))
    allocate(var%sdepth(mp,msn))
    allocate(var%smass(mp,msn))
    allocate(var%snage(mp))
    allocate(var%snowd(mp))
    allocate(var%smelt(mp))
    allocate(var%ssdn(mp,msn))
    allocate(var%ssdnn(mp))
    allocate(var%tgg(mp,ms))
    allocate(var%tggsn(mp,msn))
    allocate(var%tss(mp))
    allocate(var%tss_p(mp))
    allocate(var%deltss(mp))
    allocate(var%owb1(mp))
    allocate(var%wb(mp,ms))
    allocate(var%wbice(mp,ms))
    allocate(var%wblf(mp,ms))
    allocate(var%wbtot(mp))
    allocate(var%wbtot1(mp))
    allocate(var%wbtot2(mp))
    allocate(var%wb_lake(mp))
    allocate(var%sinfil(mp))
    allocate(var%evapfbl(mp,ms))
    allocate(var%qstss(mp))
    allocate(var%wetfac(mp))
    allocate(var%owetfac(mp))
    allocate(var%t_snwlr(mp))
    allocate(var%wbfice(mp,ms))
    allocate(var%tggav(mp))
    allocate(var%otgg(mp))
    allocate(var%otss(mp))
    allocate(var%otss_0(mp))
    allocate(var%tprecip(mp))
    allocate(var%tevap(mp))
    allocate(var%trnoff(mp))
    allocate(var%totenbal(mp))
    allocate(var%totenbal2(mp))
    allocate(var%fland(mp))
    allocate(var%ifland(mp))
    allocate(var%tilefrac(mp,n_tiles))
    allocate(var%qasrf(mp))
    allocate(var%qfsrf(mp))
    allocate(var%qssrf(mp))

    ! Allocate variables for SLI soil model:
    ! IF(cable_user%SOIL_STRUC=='sli') THEN
    allocate(var%S(mp,ms))
    allocate(var%Tsoil(mp,ms))
    allocate(var%SL(mp))
    allocate(var%TL(mp))
    allocate(var%h0(mp))
    allocate(var%rex(mp,ms))
    allocate(var%wflux(mp,0:ms))
    allocate(var%delwcol(mp))
    allocate(var%zdelta(mp))
    allocate(var%kth(mp,ms))
    allocate(var%Tsurface(mp))
    allocate(var%lE(mp))
    allocate(var%evap(mp))
    allocate(var%ciso(mp,ms+1))
    allocate(var%cisoL(mp))
    allocate(var%rlitt(mp))
    allocate(var%thetai(mp,ms))
    allocate(var%snowliq(mp,msn))
    allocate(var%nsteps(mp))
    allocate(var%nsnow(mp))
    allocate(var%TsurfaceFR(mp))
    allocate(var%Ta_daily(mp,100))
    allocate(var%Qadv_daily(mp))
    allocate(var%G0_daily(mp))
    allocate(var%Qevap_daily(mp))
    allocate(var%Qprec_daily(mp))
    allocate(var%Qprec_snow_daily(mp))
    allocate(var%E_fusion_sn(mp))
    allocate(var%E_sublimation_sn(mp))
    allocate(var%latent_heat_sn(mp))
    allocate(var%evap_liq_sn(mp))
    allocate(var%surface_melt(mp))
    allocate(var%Qadv_rain_sn(mp))
    ! END IF

  END SUBROUTINE alloc_soil_snow_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_veg_parameter_type(var, mp)

    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%canst1(mp))
    allocate(var%dleaf(mp))
    allocate(var%ejmax(mp))
    allocate(var%ejmax_shade(mp))
    allocate(var%ejmax_sun(mp))
    allocate(var%iveg(mp))
    allocate(var%ivegp(mp))
    allocate(var%iLU(mp))
    allocate(var%meth(mp))
    allocate(var%frac4(mp))
    allocate(var%hc(mp))
    allocate(var%vlai(mp))
    allocate(var%xalbnir(mp))
    allocate(var%rp20(mp))
    allocate(var%rpcoef(mp))
    allocate(var%rs20(mp))
    allocate(var%shelrb(mp))
    allocate(var%vegcf(mp))
    allocate(var%tminvj(mp))
    allocate(var%toptvj(mp))
    allocate(var%tmaxvj(mp))
    allocate(var%vbeta(mp))
    allocate(var%vcmax(mp))
    allocate(var%vcmax_shade(mp))
    allocate(var%vcmax_sun(mp))
    allocate(var%xfang(mp))
    allocate(var%extkn(mp))
    allocate(var%wai(mp))
    allocate(var%deciduous(mp))
    allocate(var%froot(mp,ms))
    !was nrb(=3), but never uses (:,3) in model
    allocate(var%refl(mp,2)) !jhan:swb?
    allocate(var%taul(mp,2))
    allocate(var%vlaimax(mp))
    allocate(var%a1gs(mp))
    allocate(var%d0gs(mp))
    allocate(var%alpha(mp))
    allocate(var%convex(mp))
    allocate(var%cfrd(mp))
    allocate(var%gswmin(mp))
    allocate(var%conkc0(mp))
    allocate(var%conko0(mp))
    allocate(var%ekc(mp))
    allocate(var%eko(mp))
    allocate(var%g0(mp))   ! Ticket #56.
    allocate(var%g1(mp))   ! Ticket #56.
    allocate(var%vcmaxcc(mp))
    allocate(var%ejmaxcc(mp))
    allocate(var%gmmax(mp))
    allocate(var%gm(mp))
    allocate(var%c4kci(mp))
    allocate(var%c4kcc(mp))
    allocate(var%bjv(mp))

    allocate(var%rootbeta(mp))
    allocate(var%gamma(mp))
    allocate(var%F10(mp))
    allocate(var%ZR(mp))
    allocate(var%clitt(mp))

    allocate(var%disturbance_interval(mp,2))
    allocate(var%disturbance_intensity(mp,2))

  END SUBROUTINE alloc_veg_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_canopy_type(var, mp)

    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%fess(mp))
    allocate(var%fesp(mp))
    allocate(var%cansto(mp))
    allocate(var%cduv(mp))
    allocate(var%delwc(mp))
    allocate(var%dewmm(mp))
    allocate(var%dgdtg(mp))
    allocate(var%fe(mp))
    allocate(var%fh(mp))
    allocate(var%fpn(mp))
    allocate(var%frp(mp))
    allocate(var%frpw(mp))
    allocate(var%frpr(mp))
    allocate(var%frs(mp))
    allocate(var%fnee(mp))
    allocate(var%frday(mp))
    allocate(var%fnv(mp))
    allocate(var%fev(mp))
    allocate(var%fevc(mp))
    allocate(var%fhv(mp))
    allocate(var%fns(mp))
    allocate(var%fhs(mp))
    allocate(var%fhs_cor(mp))
    allocate(var%ga(mp))
    allocate(var%ghflux(mp))
    allocate(var%precis(mp))
    allocate(var%qscrn(mp))
    allocate(var%rnet(mp))
    allocate(var%rniso(mp))
    allocate(var%segg(mp))
    allocate(var%sghflux(mp))
    allocate(var%through(mp))
    allocate(var%spill(mp))
    allocate(var%tscrn(mp))
    allocate(var%wcint(mp))
    allocate(var%tv(mp))
    allocate(var%us(mp))
    allocate(var%uscrn(mp))
    allocate(var%rghlai(mp))
    allocate(var%vlaiw(mp))
    allocate(var%fwet(mp))
    allocate(var%A_sh(mp))
    allocate(var%A_sl(mp))
    allocate(var%A_slC(mp))
    allocate(var%A_shC(mp))
    allocate(var%A_slJ(mp))
    allocate(var%A_shJ(mp))
    allocate(var%GPP_sh(mp))
    allocate(var%GPP_sl(mp))
    allocate(var%fevc_sh(mp))
    allocate(var%fevc_sl(mp))
    allocate(var%eta_GPP_cs(mp))
    allocate(var%eta_fevc_cs(mp))
    allocate(var%eta_A_cs(mp))
    allocate(var%eta_A_cs_sh(mp))
    allocate(var%eta_A_cs_sl(mp))
    allocate(var%eta_fevc_cs_sh(mp))
    allocate(var%eta_fevc_cs_sl(mp))
    allocate(var%cs(mp))
    allocate(var%dAdcs(mp))
    allocate(var%cs_sl(mp))
    allocate(var%cs_sh(mp))
    ! allocate(var%ci_sl(mp))
    ! allocate(var%ci_sh(mp))
    allocate(var%tlf(mp))
    allocate(var%dlf(mp))

    allocate(var%evapfbl(mp,ms))
    allocate(var%epot(mp))
    allocate(var%fnpp(mp))
    allocate(var%fevw_pot(mp))
    allocate(var%gswx_T(mp))
    allocate(var%cdtq(mp))
    allocate(var%wetfac_cs(mp))
    allocate(var%fevw(mp))
    allocate(var%fhvw(mp))
    allocate(var%fes(mp))
    allocate(var%fes_cor(mp))
    allocate(var%gswx(mp,mf))
    allocate(var%oldcansto(mp))
    allocate(var%zetar(mp,NITER))
    allocate(var%zetash(mp,NITER))
    allocate(var%fwsoil(mp))
    allocate(var%ofes(mp))

    allocate(var%gw(mp,mf))     ! dry canopy conductance (ms-1) edit vh 6/7/09
    allocate(var%ancj(mp,mf,3)) ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
    allocate(var%tlfy(mp,mf))   ! sunlit and shaded leaf temperatures
    allocate(var%ecy(mp,mf))    ! sunlit and shaded leaf transpiration (dry canopy)
    allocate(var%ecx(mp,mf))    ! sunlit and shaded leaf latent heat flux
    ! allocate(var%ci(mp,mf,3))   ! intra-cellular CO2 vh 6/7/09
    allocate(var%fwsoil(mp))

    ! vh_js - litter resistances to heat and vapour transfer
    allocate(var%kthLitt(mp))
    allocate(var%DvLitt(mp))

    ! 13C
    allocate(var%An(mp,mf))
    allocate(var%Rd(mp,mf))
    allocate(var%isc3(mp))
    allocate(var%vcmax(mp,mf))
    allocate(var%gammastar(mp,mf))
    allocate(var%gsc(mp,mf))
    allocate(var%gbc(mp,mf))
    allocate(var%gac(mp,mf))
    allocate(var%ci(mp,mf))

  END SUBROUTINE alloc_canopy_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_radiation_type(var, mp)

    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%albedo(mp,nrb))
    allocate(var%extkb(mp))
    allocate(var%extkd2(mp))
    allocate(var%extkd(mp))
    allocate(var%flws(mp))
    allocate(var%fvlai(mp,mf))
    allocate(var%latitude(mp))
    allocate(var%lwabv(mp))
    allocate(var%qcan(mp,mf,nrb))
    allocate(var%qssabs(mp))
    allocate(var%rhocdf(mp,nrb))
    allocate(var%rniso(mp,mf))
    allocate(var%scalex(mp,mf))
    allocate(var%transd(mp))
    allocate(var%trad(mp))
    allocate(var%reffdf(mp,nrb))
    allocate(var%reffbm(mp,nrb))
    allocate(var%extkbm(mp,nrb))
    allocate(var%extkdm(mp,nrb))
    allocate(var%cexpkbm(mp,swb))
    allocate(var%cexpkdm(mp,swb))
    allocate(var%fbeam(mp,nrb))
    allocate(var%rhocbm(mp,nrb))
    allocate(var%transb(mp))
    allocate(var%albedo_T(mp))
    allocate(var%gradis(mp,mf))
    allocate(var%longitude(mp))
    allocate(var%workp1(mp))
    allocate(var%workp2(mp))
    allocate(var%workp3(mp))

  END SUBROUTINE alloc_radiation_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_roughness_type(var, mp)

    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%coexp(mp))
    allocate(var%disp(mp))
    allocate(var%hruff(mp))
    allocate(var%hruff_grmx(mp))
    allocate(var%rt0us(mp))
    allocate(var%rt1usa(mp))
    allocate(var%rt1usb(mp))
    allocate(var%rt1(mp))
    allocate(var%term2(mp))
    allocate(var%term3(mp))
    allocate(var%term5(mp))
    allocate(var%term6(mp))
    allocate(var%term6a(mp))
    allocate(var%usuh(mp))
    allocate(var%za_uv(mp))
    allocate(var%za_tq(mp))
    allocate(var%z0m(mp))
    allocate(var%zref_uv(mp))
    allocate(var%zref_tq(mp))
    allocate(var%zruffs(mp))
    allocate(var%z0soilsn(mp))
    allocate(var%z0soil(mp))

  END SUBROUTINE alloc_roughness_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_air_type(var, mp)

    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%rho(mp))
    allocate(var%volm(mp))
    allocate(var%rlam(mp))
    allocate(var%qsat(mp))
    allocate(var%epsi(mp))
    allocate(var%visc(mp))
    allocate(var%psyc(mp))
    allocate(var%dsatdk(mp))
    allocate(var%cmolar(mp))

  END SUBROUTINE alloc_air_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_met_type(var, mp)

    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%ca(mp))
    allocate(var%year(mp))
    allocate(var%moy(mp))
    allocate(var%doy(mp))
    allocate(var%hod(mp))
    allocate(var%fsd(mp,swb))
    allocate(var%ofsd(mp))
    allocate(var%fld(mp))
    allocate(var%precip(mp))
    allocate(var%precip_sn(mp))
    allocate(var%tk(mp))
    allocate(var%tvair(mp))
    allocate(var%tvrad(mp))
    allocate(var%pmb(mp))
    allocate(var%ua(mp))
    allocate(var%qv(mp))
    allocate(var%qvair(mp))
    allocate(var%da(mp))
    allocate(var%dva(mp))
    allocate(var%coszen(mp))
    allocate(var%Ndep(mp))
    allocate(var%Pdep(mp))
    allocate(var%rhum(mp))
    allocate(var%u10(mp))

  END SUBROUTINE alloc_met_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_climate_type(var, mp, ktauday)

    IMPLICIT NONE

    TYPE(climate_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp, ktauday
    INTEGER :: ny, nd
    ny = var%nyear_average
    nd = var%nday_average

    allocate(var%chilldays(mp))
    allocate(var%iveg(mp))
    allocate(var%biome(mp))
    allocate(var%GMD(mp))
    allocate(var%modis_igbp(mp))
    allocate(var%DSLR(mp))
    allocate(var%NDAY_Nesterov(mp))

    allocate(var%dtemp(mp))
    allocate(var%dmoist(mp))
    allocate(var%dmoist_min(mp))
    allocate(var%dmoist_min20(mp))
    allocate(var%dmoist_max(mp))
    allocate(var%dmoist_max20(mp))
    allocate(var%mtemp(mp))
    allocate(var%qtemp(mp))
    allocate(var%mmoist(mp))
    allocate(var%mtemp_min(mp))
    allocate(var%mtemp_max(mp))
    allocate(var%qtemp_max(mp))
    allocate(var%qtemp_max_last_year(mp))
    allocate(var%mtemp_min20(mp))
    allocate(var%mtemp_max20(mp))
    allocate(var%atemp_mean(mp))
    allocate(var%AGDD5(mp))
    allocate(var%GDD5(mp))
    allocate(var%AGDD0(mp))
    allocate(var%GDD0(mp))
    allocate(var%alpha_PT(mp))
    allocate(var%evap_PT(mp))
    allocate(var%aevap(mp))
    allocate(var%alpha_PT20(mp))
    allocate(var%GDD0_rec(mp))
    allocate(var%frec(mp))
    allocate(var%dtemp_min(mp))
    allocate(var%fdorm(mp))
    allocate(var%fapar_ann_max(mp))
    allocate(var%fapar_ann_max_last_year(mp))
    allocate(var%AvgAnnMaxFAPAR(mp))
    allocate(var%dtemp_max(mp))
    allocate(var%drhum(mp))
    allocate(var%du10_max(mp))
    allocate(var%dprecip(mp))
    allocate(var%aprecip(mp))
    allocate(var%aprecip_av20(mp))
    allocate(var%last_precip(mp))
    allocate(var%KBDI(mp))
    allocate(var%FFDI(mp))
    allocate(var%D_MacArthur(mp))
    allocate(var%Nesterov_Current(mp))
    allocate(var%Nesterov_ann_max(mp))
    allocate(var%Nesterov_ann_max_last_year(mp))
    allocate(var%Nesterov_ann_running_max(mp))

    allocate(var%mtemp_min_20(mp,ny))
    allocate(var%mtemp_max_20(mp,ny))
    allocate(var%dmoist_min_20(mp,ny))
    allocate(var%dmoist_max_20(mp,ny))
    allocate(var%dtemp_31(mp,nd))
    allocate(var%dmoist_31(mp,nd))
    allocate(var%alpha_PT_20(mp,ny))
    allocate(var%dtemp_91(mp,91))
    allocate(var%APAR_leaf_sun(mp,ktauday*5))
    allocate(var%APAR_leaf_shade(mp,ktauday*5))
    allocate(var%Dleaf_sun(mp,ktauday*5))
    allocate(var%Dleaf_shade(mp,ktauday*5))
    allocate(var%Tleaf_sun(mp,ktauday*5))
    allocate(var%Tleaf_shade(mp,ktauday*5))
    allocate(var%cs_sun(mp,ktauday*5))
    allocate(var%cs_shade(mp,ktauday*5))
    allocate(var%scalex_sun(mp,ktauday*5))
    allocate(var%scalex_shade(mp,ktauday*5))
    allocate(var%fwsoil(mp,ktauday*5))
    allocate(var%aprecip_20(mp,ny))
    allocate(var%Rd_sun(mp,ktauday*5))
    allocate(var%Rd_shade(mp,ktauday*5))

  END SUBROUTINE alloc_climate_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_sum_flux_type(var, mp)

    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%sumpn(mp))
    allocate(var%sumrp(mp))
    allocate(var%sumrpw(mp))
    allocate(var%sumrpr(mp))
    allocate(var%sumrs(mp))
    allocate(var%sumrd(mp))
    allocate(var%dsumpn(mp))
    allocate(var%dsumrp(mp))
    allocate(var%dsumrs(mp))
    allocate(var%dsumrd(mp))
    allocate(var%sumxrp(mp))
    allocate(var%sumxrs(mp))

  END SUBROUTINE alloc_sum_flux_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE alloc_bgc_pool_type(var, mp)

    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    allocate(var%cplant(mp,ncp))
    allocate(var%csoil(mp,ncs))

  END SUBROUTINE alloc_bgc_pool_type

  ! ------------------------------------------------------------------------------

  ! Begin deallocation routines:
  SUBROUTINE dealloc_balances_type(var)

    TYPE(balances_type), INTENT(inout) :: var

    deallocate(var%drybal)
    deallocate(var%ebal)
    deallocate(var%ebal_tot)
    deallocate(var%ebaltr)
    deallocate(var%ebal_tottr)
    deallocate(var%ebal_cncheck)
    deallocate(var%ebal_tot_cncheck)
    deallocate(var%evap_tot)
    deallocate(var%osnowd0)
    deallocate(var%precip_tot)
    deallocate(var%rnoff_tot)
    deallocate(var%wbal)
    deallocate(var%wbal_tot)
    deallocate(var%wbtot0)
    deallocate(var%wetbal)
    deallocate(var%cansto0)
    deallocate(var%evapc_tot)
    deallocate(var%evaps_tot)
    deallocate(var%rnof1_tot)
    deallocate(var%rnof2_tot)
    deallocate(var%snowdc_tot)
    deallocate(var%wbal_tot1)
    deallocate(var%owbtot)
    deallocate(var%delwc_tot)
    deallocate(var%qasrf_tot)
    deallocate(var%qfsrf_tot)
    deallocate(var%qssrf_tot)

    deallocate(var%Radbal)
    deallocate(var%Ebalsoil)
    deallocate(var%Ebalveg)
    deallocate(var%Radbalsum)

  END SUBROUTINE dealloc_balances_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_soil_parameter_type(var)

    TYPE(soil_parameter_type), INTENT(inout) :: var

    deallocate(var%bch)
    deallocate(var%c3)
    deallocate(var%clay)
    deallocate(var%css)
    deallocate(var%hsbh)
    deallocate(var%hyds)
    deallocate(var%i2bp3)
    deallocate(var%ibp2)
    deallocate(var%isoilm)
    deallocate(var%rhosoil)
    deallocate(var%sand)
    deallocate(var%sfc)
    deallocate(var%silt)
    deallocate(var%ssat)
    deallocate(var%sucs)
    deallocate(var%swilt)
    deallocate(var%zse)
    deallocate(var%zshh)
    deallocate(var%cnsd)
    deallocate(var%albsoil)
    deallocate(var%cnsd)
    deallocate(var%pwb_min)
    deallocate(var%albsoilf)
    deallocate(var%soilcol)
    ! Deallocate variables for SLI soil model:
    ! IF(cable_user%SOIL_STRUC=='sli') THEN
    deallocate(var%nhorizons)
    deallocate(var%ishorizon)
    deallocate(var%clitt)
    deallocate(var%zeta)
    deallocate(var%fsatmax)
    deallocate(var%swilt_vec)
    deallocate(var%ssat_vec)
    deallocate(var%sfc_vec)
    if (associated(var%swilt_vec)) deallocate(var%swilt_vec)
    if (associated(var%ssat_vec))  deallocate(var%ssat_vec)
    if (associated(var%sfc_vec))   deallocate(var%sfc_vec)
    ! END IF

  END SUBROUTINE dealloc_soil_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_soil_snow_type(var)

    TYPE(soil_snow_type), INTENT(inout) :: var

    deallocate(var%iantrct)
    deallocate(var%pudsto)
    deallocate(var%pudsmx)
    deallocate(var%dtmlt)
    deallocate(var%albsoilsn)
    deallocate(var%cls)
    deallocate(var%dfn_dtg)
    deallocate(var%dfh_dtg)
    deallocate(var%dfe_ddq)
    deallocate(var%ddq_dtg)
    deallocate(var%evapsn)
    deallocate(var%fwtop)
    deallocate(var%fwtop1)
    deallocate(var%fwtop2)
    deallocate(var%fwtop3)
    deallocate(var%gammzz)
    deallocate(var%isflag)
    deallocate(var%osnowd)
    deallocate(var%potev)
    deallocate(var%runoff)
    deallocate(var%rnof1)
    deallocate(var%rnof2)
    deallocate(var%rtsoil)
    deallocate(var%sconds)
    deallocate(var%sdepth)
    deallocate(var%smass)
    deallocate(var%snage)
    deallocate(var%snowd)
    deallocate(var%smelt)
    deallocate(var%ssdn)
    deallocate(var%ssdnn)
    deallocate(var%tgg)
    deallocate(var%tggsn)
    deallocate(var%tss)
    deallocate(var%tss_p)
    deallocate(var%deltss)
    deallocate(var%owb1)
    deallocate(var%wb)
    deallocate(var%wbice)
    deallocate(var%wblf)
    deallocate(var%wbtot)
    deallocate(var%wbtot1)
    deallocate(var%wbtot2)
    deallocate(var%wb_lake)
    deallocate(var%sinfil)
    deallocate(var%evapfbl)
    deallocate(var%qstss)
    deallocate(var%wetfac)
    deallocate(var%owetfac)
    deallocate(var%t_snwlr)
    deallocate(var%wbfice)
    deallocate(var%tggav)
    deallocate(var%otgg)
    deallocate(var%otss)
    deallocate(var%otss_0)
    deallocate(var%tprecip)
    deallocate(var%tevap)
    deallocate(var%trnoff)
    deallocate(var%totenbal)
    deallocate(var%totenbal2)
    deallocate(var%fland)
    deallocate(var%ifland)
    deallocate(var%tilefrac)
    deallocate(var%qasrf)
    deallocate(var%qfsrf)
    deallocate(var%qssrf)

    ! IF(cable_user%SOIL_STRUC=='sli') THEN
    deallocate(var%S)
    deallocate(var%Tsoil)
    deallocate(var%SL)
    deallocate(var%TL)
    deallocate(var%h0)
    deallocate(var%rex)
    deallocate(var%wflux)
    deallocate(var%delwcol)
    deallocate(var%zdelta)
    deallocate(var%kth)
    deallocate(var%Tsurface)
    deallocate(var%lE)
    deallocate(var%evap)
    deallocate(var%ciso)
    deallocate(var%cisoL)
    deallocate(var%rlitt)
    deallocate(var%thetai)
    deallocate(var%snowliq)
    deallocate(var%nsteps)
    deallocate(var%nsnow)
    deallocate(var%TsurfaceFR)
    deallocate(var%Ta_daily)
    deallocate(var%G0_daily)
    deallocate(var%Qadv_daily)
    deallocate(var%Qevap_daily)
    deallocate(var%Qprec_daily)
    deallocate(var%Qprec_snow_daily)

    deallocate(var%E_fusion_sn)
    deallocate(var%E_sublimation_sn)
    deallocate(var%latent_heat_sn)
    deallocate(var%evap_liq_sn)
    deallocate(var%surface_melt)
    deallocate(var%Qadv_rain_sn)
    ! END IF

  END SUBROUTINE dealloc_soil_snow_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_veg_parameter_type(var)

    TYPE(veg_parameter_type), INTENT(inout) :: var

    deallocate(var%canst1)
    deallocate(var%dleaf)
    deallocate(var%ejmax)
    deallocate(var%ejmax_shade)
    deallocate(var%ejmax_sun)
    deallocate(var%iveg)
    deallocate(var%ivegp)
    deallocate(var%iLU)
    deallocate(var%meth)
    deallocate(var%frac4)
    deallocate(var%hc)
    deallocate(var%vlai)
    deallocate(var%xalbnir)
    deallocate(var%rp20)
    deallocate(var%rpcoef)
    deallocate(var%rs20)
    deallocate(var%shelrb)
    deallocate(var%vegcf)
    deallocate(var%tminvj)
    deallocate(var%toptvj)
    deallocate(var%tmaxvj)
    deallocate(var%vbeta)
    deallocate(var%vcmax)
    deallocate(var%vcmax_shade)
    deallocate(var%vcmax_sun)
    deallocate(var%xfang)
    deallocate(var%extkn)
    deallocate(var%wai)
    deallocate(var%deciduous)
    deallocate(var%froot)
    deallocate(var%refl)
    deallocate(var%taul)
    deallocate(var%a1gs)
    deallocate(var%d0gs)
    deallocate(var%alpha)
    deallocate(var%convex)
    deallocate(var%cfrd)
    deallocate(var%gswmin)
    deallocate(var%conkc0)
    deallocate(var%conko0)
    deallocate(var%ekc)
    deallocate(var%eko)
    deallocate(var%g0) ! Ticket #56.
    deallocate(var%g1) ! Ticket #56.

    ! Deallocate variables for SLI soil model:
    ! IF(cable_user%SOIL_STRUC=='sli') THEN
    deallocate(var%rootbeta)
    deallocate(var%gamma) ! vh 20/07/09
    deallocate(var%F10)
    deallocate(var%ZR)
    deallocate(var%CLitt)
    deallocate(var%disturbance_interval)
    deallocate(var%disturbance_intensity)
    ! END IF

  END SUBROUTINE dealloc_veg_parameter_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_canopy_type(var)

    TYPE(canopy_type), INTENT(inout) :: var

    deallocate(var%fess)
    deallocate(var%fesp)
    deallocate(var%cansto)
    deallocate(var%cduv)
    deallocate(var%delwc)
    deallocate(var%dewmm)
    deallocate(var%dgdtg)
    deallocate(var%fe)
    deallocate(var%fh)
    deallocate(var%fpn)
    deallocate(var%frp)
    deallocate(var%frpw)
    deallocate(var%frpr)
    deallocate(var%frs)
    deallocate(var%fnee)
    deallocate(var%frday)
    deallocate(var%fnv)
    deallocate(var%fev)
    deallocate(var%fevc)
    deallocate(var%fhv)
    deallocate(var%fns)
    deallocate(var%fhs)
    deallocate(var%fhs_cor)
    deallocate(var%ga)
    deallocate(var%ghflux)
    deallocate(var%precis)
    deallocate(var%qscrn)
    deallocate(var%rnet)
    deallocate(var%rniso)
    deallocate(var%segg)
    deallocate(var%sghflux)
    deallocate(var%through)
    deallocate(var%spill)
    deallocate(var%tscrn)
    deallocate(var%wcint)
    deallocate(var%tv)
    deallocate(var%us)
    deallocate(var%uscrn)
    deallocate(var%rghlai)
    deallocate(var%vlaiw)
    deallocate(var%fwet)
    deallocate(var%A_sh)
    deallocate(var%A_sl)
    deallocate(var%A_slC)
    deallocate(var%A_shC)
    deallocate(var%A_slJ)
    deallocate(var%A_shJ)
    deallocate(var%eta_A_cs)
    deallocate(var%cs)
    deallocate(var%dAdcs)
    deallocate(var%cs_sl)
    deallocate(var%cs_sh)
    deallocate(var%tlf)
    deallocate(var%dlf)

    deallocate(var%evapfbl)
    deallocate(var%epot)
    deallocate(var%fnpp)
    deallocate(var%fevw_pot)
    deallocate(var%gswx_T)
    deallocate(var%cdtq)
    deallocate(var%wetfac_cs)
    deallocate(var%fevw)
    deallocate(var%fhvw)
    deallocate(var%fes)
    deallocate(var%fes_cor)
    deallocate(var%gswx)
    deallocate(var%oldcansto)
    deallocate(var%zetar)
    deallocate(var%zetash)
    deallocate(var%fwsoil)
    deallocate(var%ofes)

    !! vh_js !! litter resistances to heat and vapour transfer
    deallocate(var%kthLitt)
    deallocate(var%DvLitt)

  END SUBROUTINE dealloc_canopy_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_radiation_type(var)

    TYPE(radiation_type), INTENT(inout) :: var

    deallocate(var%albedo)
    deallocate(var%extkb)
    deallocate(var%extkd2)
    deallocate(var%extkd)
    deallocate(var%flws)
    deallocate(var%fvlai)
    deallocate(var%latitude)
    deallocate(var%lwabv)
    deallocate(var%qcan)
    deallocate(var%qssabs)
    deallocate(var%rhocdf)
    deallocate(var%rniso)
    deallocate(var%scalex)
    deallocate(var%transd)
    deallocate(var%trad)
    deallocate(var%reffdf)
    deallocate(var%reffbm)
    deallocate(var%extkbm)
    deallocate(var%extkdm)
    deallocate(var%fbeam)
    deallocate(var%cexpkbm)
    deallocate(var%cexpkdm)
    deallocate(var%rhocbm)
    deallocate(var%transb)
    deallocate(var%albedo_T)
    deallocate(var%gradis)
    deallocate(var%longitude)
    deallocate(var%workp1)
    deallocate(var%workp2)
    deallocate(var%workp3)

  END SUBROUTINE dealloc_radiation_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_roughness_type(var)

    TYPE(roughness_type), INTENT(inout) :: var

    deallocate(var%coexp)
    deallocate(var%disp)
    deallocate(var%hruff)
    deallocate(var%hruff_grmx)
    deallocate(var%rt0us)
    deallocate(var%rt1usa)
    deallocate(var%rt1usb)
    deallocate(var%rt1)
    deallocate(var%term2)
    deallocate(var%term3)
    deallocate(var%term5)
    deallocate(var%term6)
    deallocate(var%term6a)
    deallocate(var%usuh)
    deallocate(var%za_uv)
    deallocate(var%za_tq)
    deallocate(var%z0m)
    deallocate(var%zref_uv)
    deallocate(var%zref_tq)
    deallocate(var%zruffs)
    deallocate(var%z0soilsn)
    deallocate(var%z0soil)

  END SUBROUTINE dealloc_roughness_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_air_type(var)

    TYPE(air_type), INTENT(inout) :: var

    deallocate(var%rho)
    deallocate(var%volm)
    deallocate(var%rlam)
    deallocate(var%qsat)
    deallocate(var%epsi)
    deallocate(var%visc)
    deallocate(var%psyc)
    deallocate(var%dsatdk)
    deallocate(var%cmolar)

  END SUBROUTINE dealloc_air_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_met_type(var)

    TYPE(met_type), INTENT(inout) :: var

    deallocate(var%ca)
    deallocate(var%year)
    deallocate(var%moy)
    deallocate(var%doy)
    deallocate(var%hod)
    deallocate(var%fsd)
    deallocate(var%ofsd)
    deallocate(var%fld)
    deallocate(var%precip)
    deallocate(var%precip_sn)
    deallocate(var%tk)
    deallocate(var%tvair)
    deallocate(var%tvrad)
    deallocate(var%pmb)
    deallocate(var%ua)
    deallocate(var%qv)
    deallocate(var%qvair)
    deallocate(var%da)
    deallocate(var%dva)
    deallocate(var%coszen)
    deallocate(var%Ndep)
    deallocate(var%Pdep)
    deallocate(var%rhum)
    deallocate(var%u10)

  END SUBROUTINE dealloc_met_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_sum_flux_type(var)

    TYPE(sum_flux_type), INTENT(inout) :: var

    deallocate(var%sumpn)
    deallocate(var%sumrp)
    deallocate(var%sumrpw)
    deallocate(var%sumrpr)
    deallocate(var%sumrs)
    deallocate(var%sumrd)
    deallocate(var%dsumpn)
    deallocate(var%dsumrp)
    deallocate(var%dsumrs)
    deallocate(var%dsumrd)
    deallocate(var%sumxrp)
    deallocate(var%sumxrs)

  END SUBROUTINE dealloc_sum_flux_type

  ! ------------------------------------------------------------------------------

  SUBROUTINE dealloc_bgc_pool_type(var)

    TYPE(bgc_pool_type), INTENT(inout) :: var

    deallocate(var%cplant)
    deallocate(var%csoil)

  END SUBROUTINE dealloc_bgc_pool_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_balances_type(var)

    TYPE(balances_type), INTENT(inout) :: var

    var%drybal           = 0
    var%ebal             = 0
    var%ebal_tot         = 0
    var%ebaltr           = 0
    var%ebal_tottr       = 0
    var%ebal_cncheck     = 0
    var%ebal_tot_cncheck = 0
    var%evap_tot         = 0
    var%osnowd0          = 0
    var%precip_tot       = 0
    var%rnoff_tot        = 0
    var%wbal             = 0
    var%wbal_tot         = 0
    var%wbtot0           = 0
    var%wetbal           = 0
    var%cansto0          = 0
    var%evapc_tot        = 0
    var%evaps_tot        = 0
    var%rnof1_tot        = 0
    var%rnof2_tot        = 0
    var%snowdc_tot       = 0
    var%wbal_tot1        = 0
    var%owbtot           = 0
    var%delwc_tot        = 0
    var%qasrf_tot        = 0
    var%qfsrf_tot        = 0
    var%qssrf_tot        = 0

    var%Radbal    = 0
    var%EbalSoil  = 0
    var%Ebalveg   = 0
    var%Radbalsum = 0

  END SUBROUTINE zero_balances_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_soil_parameter_type(var)

    TYPE(soil_parameter_type), INTENT(inout) :: var

    var%bch      = 0
    var%c3       = 0
    var%clay     = 0
    var%css      = 0
    var%hsbh     = 0
    var%hyds     = 0
    var%i2bp3    = 0
    var%ibp2     = 0
    var%isoilm   = 0
    var%rhosoil  = 0
    var%sand     = 0
    var%sfc      = 0
    var%silt     = 0
    var%ssat     = 0
    var%sucs     = 0
    var%swilt    = 0
    var%zse      = 0
    var%zshh     = 0
    var%cnsd     = 0
    var%albsoil  = 0
    var%pwb_min  = 0
    var%albsoilf = 0
    var%soilcol  = 0

    var%nhorizons = 0
    var%ishorizon = 0
    var%clitt     = 0
    var%zeta      = 0
    var%fsatmax   = 0
    var%swilt_vec = 0
    var%ssat_vec  = 0
    var%sfc_vec   = 0
    var%swilt_vec = 0
    var%ssat_vec  = 0
    var%sfc_vec   = 0

  END SUBROUTINE zero_soil_parameter_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_soil_snow_type(var)

    TYPE(soil_snow_type), INTENT(inout) :: var

    var%iantrct   = 0
    var%pudsto    = 0
    var%pudsmx    = 0
    var%dtmlt     = 0
    var%albsoilsn = 0
    var%cls       = 0
    var%dfn_dtg   = 0
    var%dfh_dtg   = 0
    var%dfe_ddq   = 0
    var%ddq_dtg   = 0
    var%evapsn    = 0
    var%fwtop     = 0
    var%fwtop1    = 0
    var%fwtop2    = 0
    var%fwtop3    = 0
    var%gammzz    = 0
    var%isflag    = 0
    var%osnowd    = 0
    var%potev     = 0
    var%runoff    = 0
    var%rnof1     = 0
    var%rnof2     = 0
    var%rtsoil    = 0
    var%sconds    = 0
    var%sdepth    = 0
    var%smass     = 0
    var%snage     = 0
    var%snowd     = 0
    var%smelt     = 0
    var%ssdn      = 0
    var%ssdnn     = 0
    var%tgg       = 0
    var%tggsn     = 0
    var%tss       = 0
    var%tss_p     = 0
    var%deltss    = 0
    var%owb1      = 0
    var%wb        = 0
    var%wbice     = 0
    var%wblf      = 0
    var%wbtot     = 0
    var%wbtot1    = 0
    var%wbtot2    = 0
    var%wb_lake   = 0
    var%sinfil    = 0
    var%evapfbl   = 0
    var%qstss     = 0
    var%wetfac    = 0
    var%owetfac   = 0
    var%t_snwlr   = 0
    var%wbfice    = 0
    var%tggav     = 0
    var%otgg      = 0
    var%otss      = 0
    var%otss_0    = 0
    var%tprecip   = 0
    var%tevap     = 0
    var%trnoff    = 0
    var%totenbal  = 0
    var%totenbal2 = 0
    var%fland     = 0
    var%ifland    = 0
    var%tilefrac  = 0
    var%qasrf     = 0
    var%qfsrf     = 0
    var%qssrf     = 0

    var%S                = 0
    var%Tsoil            = 0
    var%SL               = 0
    var%TL               = 0
    var%h0               = 0
    var%rex              = 0
    var%wflux            = 0
    var%delwcol          = 0
    var%zdelta           = 0
    var%kth              = 0
    var%Tsurface         = 0
    var%lE               = 0
    var%evap             = 0
    var%ciso             = 0
    var%cisoL            = 0
    var%rlitt            = 0
    var%thetai           = 0
    var%snowliq          = 0
    var%nsteps           = 0
    var%nsnow            = 0
    var%TsurfaceFR       = 0
    var%Ta_daily         = 0
    var%Qadv_daily       = 0
    var%G0_daily         = 0
    var%Qevap_daily      = 0
    var%Qprec_daily      = 0
    var%Qprec_snow_daily = 0
    var%E_fusion_sn      = 0
    var%E_sublimation_sn = 0
    var%latent_heat_sn   = 0
    var%evap_liq_sn      = 0
    var%surface_melt     = 0
    var%Qadv_rain_sn     = 0

  END SUBROUTINE zero_soil_snow_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_veg_parameter_type(var)

    TYPE(veg_parameter_type), INTENT(inout) :: var

    var%canst1      = 0
    var%dleaf       = 0
    var%ejmax       = 0
    var%ejmax_shade = 0
    var%ejmax_sun   = 0
    var%iveg        = 0
    var%ivegp       = 0
    var%iLU         = 0
    var%meth        = 0
    var%frac4       = 0
    var%hc          = 0
    var%vlai        = 0
    var%xalbnir     = 0
    var%rp20        = 0
    var%rpcoef      = 0
    var%rs20        = 0
    var%shelrb      = 0
    var%vegcf       = 0
    var%tminvj      = 0
    var%toptvj      = 0
    var%tmaxvj      = 0
    var%vbeta       = 0
    var%vcmax       = 0
    var%vcmax_shade = 0
    var%vcmax_sun   = 0
    var%xfang       = 0
    var%extkn       = 0
    var%wai         = 0
    var%deciduous   = .false.
    var%froot       = 0
    var%refl        = 0
    var%taul        = 0
    var%vlaimax     = 0
    var%a1gs        = 0
    var%d0gs        = 0
    var%alpha       = 0
    var%convex      = 0
    var%cfrd        = 0
    var%gswmin      = 0
    var%conkc0      = 0
    var%conko0      = 0
    var%ekc         = 0
    var%eko         = 0
    var%g0          = 0
    var%g1          = 0
    var%vcmaxcc     = 0
    var%ejmaxcc     = 0
    var%gmmax       = 0
    var%gm          = 0
    var%c4kci       = 0
    var%c4kcc       = 0
    var%bjv         = 0

    var%rootbeta = 0
    var%gamma    = 0
    var%F10      = 0
    var%ZR       = 0
    var%clitt    = 0

    var%disturbance_interval  = 0
    var%disturbance_intensity = 0

  END SUBROUTINE zero_veg_parameter_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_canopy_type(var)

    TYPE(canopy_type), INTENT(inout) :: var

    var%fess           = 0
    var%fesp           = 0
    var%cansto         = 0
    var%cduv           = 0
    var%delwc          = 0
    var%dewmm          = 0
    var%dgdtg          = 0
    var%fe             = 0
    var%fh             = 0
    var%fpn            = 0
    var%frp            = 0
    var%frpw           = 0
    var%frpr           = 0
    var%frs            = 0
    var%fnee           = 0
    var%frday          = 0
    var%fnv            = 0
    var%fev            = 0
    var%fevc           = 0
    var%fhv            = 0
    var%fns            = 0
    var%fhs            = 0
    var%fhs_cor        = 0
    var%ga             = 0
    var%ghflux         = 0
    var%precis         = 0
    var%qscrn          = 0
    var%rnet           = 0
    var%rniso          = 0
    var%segg           = 0
    var%sghflux        = 0
    var%through        = 0
    var%spill          = 0
    var%tscrn          = 0
    var%wcint          = 0
    var%tv             = 0
    var%us             = 0
    var%uscrn          = 0
    var%rghlai         = 0
    var%vlaiw          = 0
    var%fwet           = 0
    var%A_sh           = 0
    var%A_sl           = 0
    var%A_slC          = 0
    var%A_shC          = 0
    var%A_slJ          = 0
    var%A_shJ          = 0
    var%GPP_sh         = 0
    var%GPP_sl         = 0
    var%fevc_sh        = 0
    var%fevc_sl        = 0
    var%eta_GPP_cs     = 0
    var%eta_fevc_cs    = 0
    var%eta_A_cs       = 0
    var%eta_A_cs_sh    = 0
    var%eta_A_cs_sl    = 0
    var%eta_fevc_cs_sh = 0
    var%eta_fevc_cs_sl = 0
    var%cs             = 0
    var%dAdcs          = 0
    var%cs_sl          = 0
    var%cs_sh          = 0
    var%tlf            = 0
    var%dlf            = 0

    var%evapfbl   = 0
    var%epot      = 0
    var%fnpp      = 0
    var%fevw_pot  = 0
    var%gswx_T    = 0
    var%cdtq      = 0
    var%wetfac_cs = 0
    var%fevw      = 0
    var%fhvw      = 0
    var%fes       = 0
    var%fes_cor   = 0
    var%gswx      = 0
    var%oldcansto = 0
    var%zetar     = 0
    var%zetash    = 0
    var%fwsoil    = 0
    var%ofes      = 0

    var%gw     = 0
    var%ancj   = 0
    var%tlfy   = 0
    var%ecy    = 0
    var%ecx    = 0
    var%fwsoil = 0

    var%kthLitt = 0
    var%DvLitt  = 0

    var%An        = 0
    var%Rd        = 0
    var%isc3      = .false.
    var%vcmax     = 0
    var%gammastar = 0
    var%gsc       = 0
    var%gbc       = 0
    var%gac       = 0
    var%ci        = 0

  END SUBROUTINE zero_canopy_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_radiation_type(var)

    TYPE(radiation_type), INTENT(inout) :: var

    var%albedo    = 0
    var%extkb     = 0
    var%extkd2    = 0
    var%extkd     = 0
    var%flws      = 0
    var%fvlai     = 0
    var%latitude  = 0
    var%lwabv     = 0
    var%qcan      = 0
    var%qssabs    = 0
    var%rhocdf    = 0
    var%rniso     = 0
    var%scalex    = 0
    var%transd    = 0
    var%trad      = 0
    var%reffdf    = 0
    var%reffbm    = 0
    var%extkbm    = 0
    var%extkdm    = 0
    var%cexpkbm   = 0
    var%cexpkdm   = 0
    var%fbeam     = 0
    var%rhocbm    = 0
    var%transb    = 0
    var%albedo_T  = 0
    var%gradis    = 0
    var%longitude = 0
    var%workp1    = 0
    var%workp2    = 0
    var%workp3    = 0

  END SUBROUTINE zero_radiation_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_roughness_type(var)

    TYPE(roughness_type), INTENT(inout) :: var

    var%coexp      = 0
    var%disp       = 0
    var%hruff      = 0
    var%hruff_grmx = 0
    var%rt0us      = 0
    var%rt1usa     = 0
    var%rt1usb     = 0
    var%rt1        = 0
    var%term2      = 0
    var%term3      = 0
    var%term5      = 0
    var%term6      = 0
    var%term6a     = 0
    var%usuh       = 0
    var%za_uv      = 0
    var%za_tq      = 0
    var%z0m        = 0
    var%zref_uv    = 0
    var%zref_tq    = 0
    var%zruffs     = 0
    var%z0soilsn   = 0
    var%z0soil     = 0

  END SUBROUTINE zero_roughness_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_air_type(var)

    TYPE(air_type), INTENT(inout) :: var

    var%rho    = 0
    var%volm   = 0
    var%rlam   = 0
    var%qsat   = 0
    var%epsi   = 0
    var%visc   = 0
    var%psyc   = 0
    var%dsatdk = 0
    var%cmolar = 0

  END SUBROUTINE zero_air_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_met_type(var)

    TYPE(met_type), INTENT(inout) :: var

    var%ca        = 0
    var%year      = 0
    var%moy       = 0
    var%doy       = 0
    var%hod       = 0
    var%fsd       = 0
    var%ofsd      = 0
    var%fld       = 0
    var%precip    = 0
    var%precip_sn = 0
    var%tk        = 0
    var%tvair     = 0
    var%tvrad     = 0
    var%pmb       = 0
    var%ua        = 0
    var%qv        = 0
    var%qvair     = 0
    var%da        = 0
    var%dva       = 0
    var%coszen    = 0
    var%Ndep      = 0
    var%Pdep      = 0
    var%rhum      = 0
    var%u10       = 0

  END SUBROUTINE zero_met_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_climate_type(var)

    IMPLICIT NONE

    TYPE(climate_type), INTENT(inout) :: var

    var%chilldays     = 0
    var%iveg          = 0
    var%biome         = 0
    var%GMD           = 0
    var%modis_igbp    = 0
    var%DSLR          = 0
    var%NDAY_Nesterov = 0

    var%dtemp                      = 0
    var%dmoist                     = 0
    var%dmoist_min                 = 0
    var%dmoist_min20               = 0
    var%dmoist_max                 = 0
    var%dmoist_max20               = 0
    var%mtemp                      = 0
    var%qtemp                      = 0
    var%mmoist                     = 0
    var%mtemp_min                  = 0
    var%mtemp_max                  = 0
    var%qtemp_max                  = 0
    var%qtemp_max_last_year        = 0
    var%mtemp_min20                = 0
    var%mtemp_max20                = 0
    var%atemp_mean                 = 0
    var%AGDD5                      = 0
    var%GDD5                       = 0
    var%AGDD0                      = 0
    var%GDD0                       = 0
    var%alpha_PT                   = 0
    var%evap_PT                    = 0
    var%aevap                      = 0
    var%alpha_PT20                 = 0
    var%GDD0_rec                   = 0
    var%frec                       = 0
    var%dtemp_min                  = 0
    var%fdorm                      = 0
    var%fapar_ann_max              = 0
    var%fapar_ann_max_last_year    = 0
    var%AvgAnnMaxFAPAR             = 0
    var%dtemp_max                  = 0
    var%drhum                      = 0
    var%du10_max                   = 0
    var%dprecip                    = 0
    var%aprecip                    = 0
    var%aprecip_av20               = 0
    var%last_precip                = 0
    var%KBDI                       = 0
    var%FFDI                       = 0
    var%D_MacArthur                = 0
    var%Nesterov_Current           = 0
    var%Nesterov_ann_max           = 0
    var%Nesterov_ann_max_last_year = 0
    var%Nesterov_ann_running_max   = 0

    var%mtemp_min_20    = 0
    var%mtemp_max_20    = 0
    var%dmoist_min_20   = 0
    var%dmoist_max_20   = 0
    var%dtemp_31        = 0
    var%dmoist_31       = 0
    var%alpha_PT_20     = 0
    var%dtemp_91        = 0
    var%APAR_leaf_sun   = 0
    var%APAR_leaf_shade = 0
    var%Dleaf_sun       = 0
    var%Dleaf_shade     = 0
    var%Tleaf_sun       = 0
    var%Tleaf_shade     = 0
    var%cs_sun          = 0
    var%cs_shade        = 0
    var%scalex_sun      = 0
    var%scalex_shade    = 0
    var%fwsoil          = 0
    var%aprecip_20      = 0
    var%Rd_sun          = 0
    var%Rd_shade        = 0

  END SUBROUTINE zero_climate_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_sum_flux_type(var)

    TYPE(sum_flux_type), INTENT(inout) :: var

    var%sumpn  = 0
    var%sumrp  = 0
    var%sumrpw = 0
    var%sumrpr = 0
    var%sumrs  = 0
    var%sumrd  = 0
    var%dsumpn = 0
    var%dsumrp = 0
    var%dsumrs = 0
    var%dsumrd = 0
    var%sumxrp = 0
    var%sumxrs = 0

  END SUBROUTINE zero_sum_flux_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE zero_bgc_pool_type(var)

    TYPE(bgc_pool_type), INTENT(inout) :: var

    var%cplant = 0
    var%csoil  = 0

  END SUBROUTINE zero_bgc_pool_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_balances_type(var)

    TYPE(balances_type), INTENT(in) :: var

    write(*,*) 'bal%drybal ', var%drybal
    write(*,*) 'bal%ebal ', var%ebal
    write(*,*) 'bal%ebal_tot ', var%ebal_tot
    write(*,*) 'bal%ebaltr ', var%ebaltr
    write(*,*) 'bal%ebal_tottr ', var%ebal_tottr
    write(*,*) 'bal%ebal_cncheck ', var%ebal_cncheck
    write(*,*) 'bal%ebal_tot_cncheck ', var%ebal_tot_cncheck
    write(*,*) 'bal%evap_tot ', var%evap_tot
    write(*,*) 'bal%osnowd0 ', var%osnowd0
    write(*,*) 'bal%precip_tot ', var%precip_tot
    write(*,*) 'bal%rnoff_tot ', var%rnoff_tot
    write(*,*) 'bal%wbal ', var%wbal
    write(*,*) 'bal%wbal_tot ', var%wbal_tot
    write(*,*) 'bal%wbtot0 ', var%wbtot0
    write(*,*) 'bal%wetbal ', var%wetbal
    write(*,*) 'bal%cansto0 ', var%cansto0
    write(*,*) 'bal%evapc_tot ', var%evapc_tot
    write(*,*) 'bal%evaps_tot ', var%evaps_tot
    write(*,*) 'bal%rnof1_tot ', var%rnof1_tot
    write(*,*) 'bal%rnof2_tot ', var%rnof2_tot
    write(*,*) 'bal%snowdc_tot ', var%snowdc_tot
    write(*,*) 'bal%wbal_tot1 ', var%wbal_tot1
    write(*,*) 'bal%owbtot ', var%owbtot
    write(*,*) 'bal%delwc_tot ', var%delwc_tot
    write(*,*) 'bal%qasrf_tot ', var%qasrf_tot
    write(*,*) 'bal%qfsrf_tot ', var%qfsrf_tot
    write(*,*) 'bal%qssrf_tot ', var%qssrf_tot

    write(*,*) 'bal%Radbal ', var%Radbal
    write(*,*) 'bal%EbalSoil ', var%EbalSoil
    write(*,*) 'bal%Ebalveg ', var%Ebalveg
    write(*,*) 'bal%Radbalsum ', var%Radbalsum

  END SUBROUTINE print_balances_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_soil_parameter_type(var)

    TYPE(soil_parameter_type), INTENT(in) :: var

    write(*,*) 'soil%bch ', var%bch
    write(*,*) 'soil%c3 ', var%c3
    write(*,*) 'soil%clay ', var%clay
    write(*,*) 'soil%css ', var%css
    write(*,*) 'soil%hsbh ', var%hsbh
    write(*,*) 'soil%hyds ', var%hyds
    write(*,*) 'soil%i2bp3 ', var%i2bp3
    write(*,*) 'soil%ibp2 ', var%ibp2
    write(*,*) 'soil%isoilm ', var%isoilm
    write(*,*) 'soil%rhosoil ', var%rhosoil
    write(*,*) 'soil%sand ', var%sand
    write(*,*) 'soil%sfc ', var%sfc
    write(*,*) 'soil%silt ', var%silt
    write(*,*) 'soil%ssat ', var%ssat
    write(*,*) 'soil%sucs ', var%sucs
    write(*,*) 'soil%swilt ', var%swilt
    write(*,*) 'soil%zse ', var%zse
    write(*,*) 'soil%zshh ', var%zshh
    write(*,*) 'soil%cnsd ', var%cnsd
    write(*,*) 'soil%albsoil ', var%albsoil
    write(*,*) 'soil%pwb_min ', var%pwb_min
    write(*,*) 'soil%albsoilf ', var%albsoilf
    write(*,*) 'soil%soilcol ', var%soilcol

    write(*,*) 'soil%nhorizons ', var%nhorizons
    write(*,*) 'soil%ishorizon ', var%ishorizon
    write(*,*) 'soil%clitt ', var%clitt
    write(*,*) 'soil%zeta ', var%zeta
    write(*,*) 'soil%fsatmax ', var%fsatmax
    write(*,*) 'soil%swilt_vec ', var%swilt_vec
    write(*,*) 'soil%ssat_vec ', var%ssat_vec
    write(*,*) 'soil%sfc_vec ', var%sfc_vec
    write(*,*) 'soil%swilt_vec ', var%swilt_vec
    write(*,*) 'soil%ssat_vec ', var%ssat_vec
    write(*,*) 'soil%sfc_vec ', var%sfc_vec

  END SUBROUTINE print_soil_parameter_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_soil_snow_type(var)

    TYPE(soil_snow_type), INTENT(in) :: var

    write(*,*) 'ssnow%iantrct ', var%iantrct
    write(*,*) 'ssnow%pudsto ', var%pudsto
    write(*,*) 'ssnow%pudsmx ', var%pudsmx
    write(*,*) 'ssnow%dtmlt ', var%dtmlt
    write(*,*) 'ssnow%albsoilsn ', var%albsoilsn
    write(*,*) 'ssnow%cls ', var%cls
    write(*,*) 'ssnow%dfn_dtg ', var%dfn_dtg
    write(*,*) 'ssnow%dfh_dtg ', var%dfh_dtg
    write(*,*) 'ssnow%dfe_ddq ', var%dfe_ddq
    write(*,*) 'ssnow%ddq_dtg ', var%ddq_dtg
    write(*,*) 'ssnow%evapsn ', var%evapsn
    write(*,*) 'ssnow%fwtop ', var%fwtop
    write(*,*) 'ssnow%fwtop1 ', var%fwtop1
    write(*,*) 'ssnow%fwtop2 ', var%fwtop2
    write(*,*) 'ssnow%fwtop3 ', var%fwtop3
    write(*,*) 'ssnow%gammzz ', var%gammzz
    write(*,*) 'ssnow%isflag ', var%isflag
    write(*,*) 'ssnow%osnowd ', var%osnowd
    write(*,*) 'ssnow%potev ', var%potev
    write(*,*) 'ssnow%runoff ', var%runoff
    write(*,*) 'ssnow%rnof1 ', var%rnof1
    write(*,*) 'ssnow%rnof2 ', var%rnof2
    write(*,*) 'ssnow%rtsoil ', var%rtsoil
    write(*,*) 'ssnow%sconds ', var%sconds
    write(*,*) 'ssnow%sdepth ', var%sdepth
    write(*,*) 'ssnow%smass ', var%smass
    write(*,*) 'ssnow%snage ', var%snage
    write(*,*) 'ssnow%snowd ', var%snowd
    write(*,*) 'ssnow%smelt ', var%smelt
    write(*,*) 'ssnow%ssdn ', var%ssdn
    write(*,*) 'ssnow%ssdnn ', var%ssdnn
    write(*,*) 'ssnow%tgg ', var%tgg
    write(*,*) 'ssnow%tggsn ', var%tggsn
    write(*,*) 'ssnow%tss ', var%tss
    write(*,*) 'ssnow%tss_p ', var%tss_p
    write(*,*) 'ssnow%deltss ', var%deltss
    write(*,*) 'ssnow%owb1 ', var%owb1
    write(*,*) 'ssnow%wb ', var%wb
    write(*,*) 'ssnow%wbice ', var%wbice
    write(*,*) 'ssnow%wblf ', var%wblf
    write(*,*) 'ssnow%wbtot ', var%wbtot
    write(*,*) 'ssnow%wbtot1 ', var%wbtot1
    write(*,*) 'ssnow%wbtot2 ', var%wbtot2
    write(*,*) 'ssnow%wb_lake ', var%wb_lake
    write(*,*) 'ssnow%sinfil ', var%sinfil
    write(*,*) 'ssnow%evapfbl ', var%evapfbl
    write(*,*) 'ssnow%qstss ', var%qstss
    write(*,*) 'ssnow%wetfac ', var%wetfac
    write(*,*) 'ssnow%owetfac ', var%owetfac
    write(*,*) 'ssnow%t_snwlr ', var%t_snwlr
    write(*,*) 'ssnow%wbfice ', var%wbfice
    write(*,*) 'ssnow%tggav ', var%tggav
    write(*,*) 'ssnow%otgg ', var%otgg
    write(*,*) 'ssnow%otss ', var%otss
    write(*,*) 'ssnow%otss_0 ', var%otss_0
    write(*,*) 'ssnow%tprecip ', var%tprecip
    write(*,*) 'ssnow%tevap ', var%tevap
    write(*,*) 'ssnow%trnoff ', var%trnoff
    write(*,*) 'ssnow%totenbal ', var%totenbal
    write(*,*) 'ssnow%totenbal2 ', var%totenbal2
    write(*,*) 'ssnow%fland ', var%fland
    write(*,*) 'ssnow%ifland ', var%ifland
    write(*,*) 'ssnow%tilefrac ', var%tilefrac
    write(*,*) 'ssnow%qasrf ', var%qasrf
    write(*,*) 'ssnow%qfsrf ', var%qfsrf
    write(*,*) 'ssnow%qssrf ', var%qssrf

    write(*,*) 'ssnow%S ', var%S
    write(*,*) 'ssnow%Tsoil ', var%Tsoil
    write(*,*) 'ssnow%SL ', var%SL
    write(*,*) 'ssnow%TL ', var%TL
    write(*,*) 'ssnow%h0 ', var%h0
    write(*,*) 'ssnow%rex ', var%rex
    write(*,*) 'ssnow%wflux ', var%wflux
    write(*,*) 'ssnow%delwcol ', var%delwcol
    write(*,*) 'ssnow%zdelta ', var%zdelta
    write(*,*) 'ssnow%kth ', var%kth
    write(*,*) 'ssnow%Tsurface ', var%Tsurface
    write(*,*) 'ssnow%lE ', var%lE
    write(*,*) 'ssnow%evap ', var%evap
    write(*,*) 'ssnow%ciso ', var%ciso
    write(*,*) 'ssnow%cisoL ', var%cisoL
    write(*,*) 'ssnow%rlitt ', var%rlitt
    write(*,*) 'ssnow%thetai ', var%thetai
    write(*,*) 'ssnow%snowliq ', var%snowliq
    write(*,*) 'ssnow%nsteps ', var%nsteps
    write(*,*) 'ssnow%nsnow ', var%nsnow
    write(*,*) 'ssnow%TsurfaceFR ', var%TsurfaceFR
    write(*,*) 'ssnow%Ta_daily ', var%Ta_daily
    write(*,*) 'ssnow%Qadv_daily ', var%Qadv_daily
    write(*,*) 'ssnow%G0_daily ', var%G0_daily
    write(*,*) 'ssnow%Qevap_daily ', var%Qevap_daily
    write(*,*) 'ssnow%Qprec_daily ', var%Qprec_daily
    write(*,*) 'ssnow%Qprec_snow_daily ', var%Qprec_snow_daily
    write(*,*) 'ssnow%E_fusion_sn ', var%E_fusion_sn
    write(*,*) 'ssnow%E_sublimation_sn ', var%E_sublimation_sn
    write(*,*) 'ssnow%latent_heat_sn ', var%latent_heat_sn
    write(*,*) 'ssnow%evap_liq_sn ', var%evap_liq_sn
    write(*,*) 'ssnow%surface_melt ', var%surface_melt
    write(*,*) 'ssnow%Qadv_rain_sn ', var%Qadv_rain_sn

  END SUBROUTINE print_soil_snow_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_veg_parameter_type(var)

    TYPE(veg_parameter_type), INTENT(in) :: var

    write(*,*) 'veg%canst1 ', var%canst1
    write(*,*) 'veg%dleaf ', var%dleaf
    write(*,*) 'veg%ejmax ', var%ejmax
    write(*,*) 'veg%ejmax_shade ', var%ejmax_shade
    write(*,*) 'veg%ejmax_sun ', var%ejmax_sun
    write(*,*) 'veg%iveg ', var%iveg
    write(*,*) 'veg%ivegp ', var%ivegp
    write(*,*) 'veg%iLU ', var%iLU
    write(*,*) 'veg%meth ', var%meth
    write(*,*) 'veg%frac4 ', var%frac4
    write(*,*) 'veg%hc ', var%hc
    write(*,*) 'veg%vlai ', var%vlai
    write(*,*) 'veg%xalbnir ', var%xalbnir
    write(*,*) 'veg%rp20 ', var%rp20
    write(*,*) 'veg%rpcoef ', var%rpcoef
    write(*,*) 'veg%rs20 ', var%rs20
    write(*,*) 'veg%shelrb ', var%shelrb
    write(*,*) 'veg%vegcf ', var%vegcf
    write(*,*) 'veg%tminvj ', var%tminvj
    write(*,*) 'veg%toptvj ', var%toptvj
    write(*,*) 'veg%tmaxvj ', var%tmaxvj
    write(*,*) 'veg%vbeta ', var%vbeta
    write(*,*) 'veg%vcmax ', var%vcmax
    write(*,*) 'veg%vcmax_shade ', var%vcmax_shade
    write(*,*) 'veg%vcmax_sun ', var%vcmax_sun
    write(*,*) 'veg%xfang ', var%xfang
    write(*,*) 'veg%extkn ', var%extkn
    write(*,*) 'veg%wai ', var%wai
    write(*,*) 'veg%deciduous ', var%deciduous
    write(*,*) 'veg%froot ', var%froot
    write(*,*) 'veg%refl ', var%refl
    write(*,*) 'veg%taul ', var%taul
    write(*,*) 'veg%vlaimax ', var%vlaimax
    write(*,*) 'veg%a1gs ', var%a1gs
    write(*,*) 'veg%d0gs ', var%d0gs
    write(*,*) 'veg%alpha ', var%alpha
    write(*,*) 'veg%convex ', var%convex
    write(*,*) 'veg%cfrd ', var%cfrd
    write(*,*) 'veg%gswmin ', var%gswmin
    write(*,*) 'veg%conkc0 ', var%conkc0
    write(*,*) 'veg%conko0 ', var%conko0
    write(*,*) 'veg%ekc ', var%ekc
    write(*,*) 'veg%eko ', var%eko
    write(*,*) 'veg%g0 ', var%g0
    write(*,*) 'veg%g1 ', var%g1
    write(*,*) 'veg%vcmaxcc ', var%vcmaxcc
    write(*,*) 'veg%ejmaxcc ', var%ejmaxcc
    write(*,*) 'veg%gmmax ', var%gmmax
    write(*,*) 'veg%gm ', var%gm
    write(*,*) 'veg%c4kci ', var%c4kci
    write(*,*) 'veg%c4kcc ', var%c4kcc
    write(*,*) 'veg%bjv ', var%bjv

    write(*,*) 'veg%rootbeta ', var%rootbeta
    write(*,*) 'veg%gamma ', var%gamma
    write(*,*) 'veg%F10 ', var%F10
    write(*,*) 'veg%ZR ', var%ZR
    write(*,*) 'veg%clitt ', var%clitt

    write(*,*) 'veg%disturbance_interval ', var%disturbance_interval
    write(*,*) 'veg%disturbance_intensity ', var%disturbance_intensity

  END SUBROUTINE print_veg_parameter_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_canopy_type(var)

    TYPE(canopy_type), INTENT(in) :: var

    write(*,*) 'canopy%fess ', var%fess
    write(*,*) 'canopy%fesp ', var%fesp
    write(*,*) 'canopy%cansto ', var%cansto
    write(*,*) 'canopy%cduv ', var%cduv
    write(*,*) 'canopy%delwc ', var%delwc
    write(*,*) 'canopy%dewmm ', var%dewmm
    write(*,*) 'canopy%dgdtg ', var%dgdtg
    write(*,*) 'canopy%fe ', var%fe
    write(*,*) 'canopy%fh ', var%fh
    write(*,*) 'canopy%fpn ', var%fpn
    write(*,*) 'canopy%frp ', var%frp
    write(*,*) 'canopy%frpw ', var%frpw
    write(*,*) 'canopy%frpr ', var%frpr
    write(*,*) 'canopy%frs ', var%frs
    write(*,*) 'canopy%fnee ', var%fnee
    write(*,*) 'canopy%frday ', var%frday
    write(*,*) 'canopy%fnv ', var%fnv
    write(*,*) 'canopy%fev ', var%fev
    write(*,*) 'canopy%fevc ', var%fevc
    write(*,*) 'canopy%fhv ', var%fhv
    write(*,*) 'canopy%fns ', var%fns
    write(*,*) 'canopy%fhs ', var%fhs
    write(*,*) 'canopy%fhs_cor ', var%fhs_cor
    write(*,*) 'canopy%ga ', var%ga
    write(*,*) 'canopy%ghflux ', var%ghflux
    write(*,*) 'canopy%precis ', var%precis
    write(*,*) 'canopy%qscrn ', var%qscrn
    write(*,*) 'canopy%rnet ', var%rnet
    write(*,*) 'canopy%rniso ', var%rniso
    write(*,*) 'canopy%segg ', var%segg
    write(*,*) 'canopy%sghflux ', var%sghflux
    write(*,*) 'canopy%through ', var%through
    write(*,*) 'canopy%spill ', var%spill
    write(*,*) 'canopy%tscrn ', var%tscrn
    write(*,*) 'canopy%wcint ', var%wcint
    write(*,*) 'canopy%tv ', var%tv
    write(*,*) 'canopy%us ', var%us
    write(*,*) 'canopy%uscrn ', var%uscrn
    write(*,*) 'canopy%rghlai ', var%rghlai
    write(*,*) 'canopy%vlaiw ', var%vlaiw
    write(*,*) 'canopy%fwet ', var%fwet
    write(*,*) 'canopy%A_sh ', var%A_sh
    write(*,*) 'canopy%A_sl ', var%A_sl
    write(*,*) 'canopy%A_slC ', var%A_slC
    write(*,*) 'canopy%A_shC ', var%A_shC
    write(*,*) 'canopy%A_slJ ', var%A_slJ
    write(*,*) 'canopy%A_shJ ', var%A_shJ
    write(*,*) 'canopy%GPP_sh ', var%GPP_sh
    write(*,*) 'canopy%GPP_sl ', var%GPP_sl
    write(*,*) 'canopy%fevc_sh ', var%fevc_sh
    write(*,*) 'canopy%fevc_sl ', var%fevc_sl
    write(*,*) 'canopy%eta_GPP_cs ', var%eta_GPP_cs
    write(*,*) 'canopy%eta_fevc_cs ', var%eta_fevc_cs
    write(*,*) 'canopy%eta_A_cs ', var%eta_A_cs
    write(*,*) 'canopy%eta_A_cs_sh ', var%eta_A_cs_sh
    write(*,*) 'canopy%eta_A_cs_sl ', var%eta_A_cs_sl
    write(*,*) 'canopy%eta_fevc_cs_sh ', var%eta_fevc_cs_sh
    write(*,*) 'canopy%eta_fevc_cs_sl ', var%eta_fevc_cs_sl
    write(*,*) 'canopy%cs ', var%cs
    write(*,*) 'canopy%dAdcs ', var%dAdcs
    write(*,*) 'canopy%cs_sl ', var%cs_sl
    write(*,*) 'canopy%cs_sh ', var%cs_sh
    write(*,*) 'canopy%tlf ', var%tlf
    write(*,*) 'canopy%dlf ', var%dlf

    write(*,*) 'canopy%evapfbl ', var%evapfbl
    write(*,*) 'canopy%epot ', var%epot
    write(*,*) 'canopy%fnpp ', var%fnpp
    write(*,*) 'canopy%fevw_pot ', var%fevw_pot
    write(*,*) 'canopy%gswx_T ', var%gswx_T
    write(*,*) 'canopy%cdtq ', var%cdtq
    write(*,*) 'canopy%wetfac_cs ', var%wetfac_cs
    write(*,*) 'canopy%fevw ', var%fevw
    write(*,*) 'canopy%fhvw ', var%fhvw
    write(*,*) 'canopy%fes ', var%fes
    write(*,*) 'canopy%fes_cor ', var%fes_cor
    write(*,*) 'canopy%gswx ', var%gswx
    write(*,*) 'canopy%oldcansto ', var%oldcansto
    write(*,*) 'canopy%zetar ', var%zetar
    write(*,*) 'canopy%zetash ', var%zetash
    write(*,*) 'canopy%fwsoil ', var%fwsoil
    write(*,*) 'canopy%ofes ', var%ofes

    write(*,*) 'canopy%gw ', var%gw
    write(*,*) 'canopy%ancj ', var%ancj
    write(*,*) 'canopy%tlfy ', var%tlfy
    write(*,*) 'canopy%ecy ', var%ecy
    write(*,*) 'canopy%ecx ', var%ecx
    write(*,*) 'canopy%fwsoil ', var%fwsoil

    write(*,*) 'canopy%kthLitt ', var%kthLitt
    write(*,*) 'canopy%DvLitt ', var%DvLitt

    write(*,*) 'canopy%An ', var%An
    write(*,*) 'canopy%Rd ', var%Rd
    write(*,*) 'canopy%isc3 ', var%isc3
    write(*,*) 'canopy%vcmax ', var%vcmax
    write(*,*) 'canopy%gammastar ', var%gammastar
    write(*,*) 'canopy%gsc ', var%gsc
    write(*,*) 'canopy%gbc ', var%gbc
    write(*,*) 'canopy%gac ', var%gac
    write(*,*) 'canopy%ci ', var%ci

  END SUBROUTINE print_canopy_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_radiation_type(var)

    TYPE(radiation_type), INTENT(in) :: var

    write(*,*) 'rad%albedo ', var%albedo
    write(*,*) 'rad%extkb ', var%extkb
    write(*,*) 'rad%extkd2 ', var%extkd2
    write(*,*) 'rad%extkd ', var%extkd
    write(*,*) 'rad%flws ', var%flws
    write(*,*) 'rad%fvlai ', var%fvlai
    write(*,*) 'rad%latitude ', var%latitude
    write(*,*) 'rad%lwabv ', var%lwabv
    write(*,*) 'rad%qcan ', var%qcan
    write(*,*) 'rad%qssabs ', var%qssabs
    write(*,*) 'rad%rhocdf ', var%rhocdf
    write(*,*) 'rad%rniso ', var%rniso
    write(*,*) 'rad%scalex ', var%scalex
    write(*,*) 'rad%transd ', var%transd
    write(*,*) 'rad%trad ', var%trad
    write(*,*) 'rad%reffdf ', var%reffdf
    write(*,*) 'rad%reffbm ', var%reffbm
    write(*,*) 'rad%extkbm ', var%extkbm
    write(*,*) 'rad%extkdm ', var%extkdm
    write(*,*) 'rad%cexpkbm ', var%cexpkbm
    write(*,*) 'rad%cexpkdm ', var%cexpkdm
    write(*,*) 'rad%fbeam ', var%fbeam
    write(*,*) 'rad%rhocbm ', var%rhocbm
    write(*,*) 'rad%transb ', var%transb
    write(*,*) 'rad%albedo_T ', var%albedo_T
    write(*,*) 'rad%gradis ', var%gradis
    write(*,*) 'rad%longitude ', var%longitude
    write(*,*) 'rad%workp1 ', var%workp1
    write(*,*) 'rad%workp2 ', var%workp2
    write(*,*) 'rad%workp3 ', var%workp3

  END SUBROUTINE print_radiation_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_roughness_type(var)

    TYPE(roughness_type), INTENT(in) :: var

    write(*,*) 'rough%coexp ', var%coexp
    write(*,*) 'rough%disp ', var%disp
    write(*,*) 'rough%hruff ', var%hruff
    write(*,*) 'rough%hruff_grmx ', var%hruff_grmx
    write(*,*) 'rough%rt0us ', var%rt0us
    write(*,*) 'rough%rt1usa ', var%rt1usa
    write(*,*) 'rough%rt1usb ', var%rt1usb
    write(*,*) 'rough%rt1 ', var%rt1
    write(*,*) 'rough%term2 ', var%term2
    write(*,*) 'rough%term3 ', var%term3
    write(*,*) 'rough%term5 ', var%term5
    write(*,*) 'rough%term6 ', var%term6
    write(*,*) 'rough%term6a ', var%term6a
    write(*,*) 'rough%usuh ', var%usuh
    write(*,*) 'rough%za_uv ', var%za_uv
    write(*,*) 'rough%za_tq ', var%za_tq
    write(*,*) 'rough%z0m ', var%z0m
    write(*,*) 'rough%zref_uv ', var%zref_uv
    write(*,*) 'rough%zref_tq ', var%zref_tq
    write(*,*) 'rough%zruffs ', var%zruffs
    write(*,*) 'rough%z0soilsn ', var%z0soilsn
    write(*,*) 'rough%z0soil ', var%z0soil

  END SUBROUTINE print_roughness_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_air_type(var)

    TYPE(air_type), INTENT(in) :: var

    write(*,*) 'air%rho ', var%rho
    write(*,*) 'air%volm ', var%volm
    write(*,*) 'air%rlam ', var%rlam
    write(*,*) 'air%qsat ', var%qsat
    write(*,*) 'air%epsi ', var%epsi
    write(*,*) 'air%visc ', var%visc
    write(*,*) 'air%psyc ', var%psyc
    write(*,*) 'air%dsatdk ', var%dsatdk
    write(*,*) 'air%cmolar ', var%cmolar

  END SUBROUTINE print_air_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_met_type(var)

    TYPE(met_type), INTENT(in) :: var

    write(*,*) 'met%ca ', var%ca
    write(*,*) 'met%year ', var%year
    write(*,*) 'met%moy ', var%moy
    write(*,*) 'met%doy ', var%doy
    write(*,*) 'met%hod ', var%hod
    write(*,*) 'met%fsd ', var%fsd
    write(*,*) 'met%ofsd ', var%ofsd
    write(*,*) 'met%fld ', var%fld
    write(*,*) 'met%precip ', var%precip
    write(*,*) 'met%precip_sn ', var%precip_sn
    write(*,*) 'met%tk ', var%tk
    write(*,*) 'met%tvair ', var%tvair
    write(*,*) 'met%tvrad ', var%tvrad
    write(*,*) 'met%pmb ', var%pmb
    write(*,*) 'met%ua ', var%ua
    write(*,*) 'met%qv ', var%qv
    write(*,*) 'met%qvair ', var%qvair
    write(*,*) 'met%da ', var%da
    write(*,*) 'met%dva ', var%dva
    write(*,*) 'met%coszen ', var%coszen
    write(*,*) 'met%Ndep ', var%Ndep
    write(*,*) 'met%Pdep ', var%Pdep
    write(*,*) 'met%rhum ', var%rhum
    write(*,*) 'met%u10 ', var%u10

  END SUBROUTINE print_met_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_climate_type(var)

    IMPLICIT NONE

    TYPE(climate_type), INTENT(in) :: var

    write(*,*) 'climate%chilldays ', var%chilldays
    write(*,*) 'climate%iveg ', var%iveg
    write(*,*) 'climate%biome ', var%biome
    write(*,*) 'climate%GMD ', var%GMD
    write(*,*) 'climate%modis_igbp ', var%modis_igbp
    write(*,*) 'climate%DSLR ', var%DSLR
    write(*,*) 'climate%NDAY_Nesterov ', var%NDAY_Nesterov

    write(*,*) 'climate%dtemp ', var%dtemp
    write(*,*) 'climate%dmoist ', var%dmoist
    write(*,*) 'climate%dmoist_min ', var%dmoist_min
    write(*,*) 'climate%dmoist_min20 ', var%dmoist_min20
    write(*,*) 'climate%dmoist_max ', var%dmoist_max
    write(*,*) 'climate%dmoist_max20 ', var%dmoist_max20
    write(*,*) 'climate%mtemp ', var%mtemp
    write(*,*) 'climate%qtemp ', var%qtemp
    write(*,*) 'climate%mmoist ', var%mmoist
    write(*,*) 'climate%mtemp_min ', var%mtemp_min
    write(*,*) 'climate%mtemp_max ', var%mtemp_max
    write(*,*) 'climate%qtemp_max ', var%qtemp_max
    write(*,*) 'climate%qtemp_max_last_year ', var%qtemp_max_last_year
    write(*,*) 'climate%mtemp_min20 ', var%mtemp_min20
    write(*,*) 'climate%mtemp_max20 ', var%mtemp_max20
    write(*,*) 'climate%atemp_mean ', var%atemp_mean
    write(*,*) 'climate%AGDD5 ', var%AGDD5
    write(*,*) 'climate%GDD5 ', var%GDD5
    write(*,*) 'climate%AGDD0 ', var%AGDD0
    write(*,*) 'climate%GDD0 ', var%GDD0
    write(*,*) 'climate%alpha_PT ', var%alpha_PT
    write(*,*) 'climate%evap_PT ', var%evap_PT
    write(*,*) 'climate%aevap ', var%aevap
    write(*,*) 'climate%alpha_PT20 ', var%alpha_PT20
    write(*,*) 'climate%GDD0_rec ', var%GDD0_rec
    write(*,*) 'climate%frec ', var%frec
    write(*,*) 'climate%dtemp_min ', var%dtemp_min
    write(*,*) 'climate%fdorm ', var%fdorm
    write(*,*) 'climate%fapar_ann_max ', var%fapar_ann_max
    write(*,*) 'climate%fapar_ann_max_last_year ', var%fapar_ann_max_last_year
    write(*,*) 'climate%AvgAnnMaxFAPAR ', var%AvgAnnMaxFAPAR
    write(*,*) 'climate%dtemp_max ', var%dtemp_max
    write(*,*) 'climate%drhum ', var%drhum
    write(*,*) 'climate%du10_max ', var%du10_max
    write(*,*) 'climate%dprecip ', var%dprecip
    write(*,*) 'climate%aprecip ', var%aprecip
    write(*,*) 'climate%aprecip_av20 ', var%aprecip_av20
    write(*,*) 'climate%last_precip ', var%last_precip
    write(*,*) 'climate%KBDI ', var%KBDI
    write(*,*) 'climate%FFDI ', var%FFDI
    write(*,*) 'climate%D_MacArthur ', var%D_MacArthur
    write(*,*) 'climate%Nesterov_Current ', var%Nesterov_Current
    write(*,*) 'climate%Nesterov_ann_max ', var%Nesterov_ann_max
    write(*,*) 'climate%Nesterov_ann_max_last_year ', var%Nesterov_ann_max_last_year
    write(*,*) 'climate%Nesterov_ann_running_max ', var%Nesterov_ann_running_max

    write(*,*) 'climate%mtemp_min_20 ', var%mtemp_min_20
    write(*,*) 'climate%mtemp_max_20 ', var%mtemp_max_20
    write(*,*) 'climate%dmoist_min_20 ', var%dmoist_min_20
    write(*,*) 'climate%dmoist_max_20 ', var%dmoist_max_20
    write(*,*) 'climate%dtemp_31 ', var%dtemp_31
    write(*,*) 'climate%dmoist_31 ', var%dmoist_31
    write(*,*) 'climate%alpha_PT_20 ', var%alpha_PT_20
    write(*,*) 'climate%dtemp_91 ', var%dtemp_91
    write(*,*) 'climate%APAR_leaf_sun ', var%APAR_leaf_sun
    write(*,*) 'climate%APAR_leaf_shade ', var%APAR_leaf_shade
    write(*,*) 'climate%Dleaf_sun ', var%Dleaf_sun
    write(*,*) 'climate%Dleaf_shade ', var%Dleaf_shade
    write(*,*) 'climate%Tleaf_sun ', var%Tleaf_sun
    write(*,*) 'climate%Tleaf_shade ', var%Tleaf_shade
    write(*,*) 'climate%cs_sun ', var%cs_sun
    write(*,*) 'climate%cs_shade ', var%cs_shade
    write(*,*) 'climate%scalex_sun ', var%scalex_sun
    write(*,*) 'climate%scalex_shade ', var%scalex_shade
    write(*,*) 'climate%fwsoil ', var%fwsoil
    write(*,*) 'climate%aprecip_20 ', var%aprecip_20
    write(*,*) 'climate%Rd_sun ', var%Rd_sun
    write(*,*) 'climate%Rd_shade ', var%Rd_shade

  END SUBROUTINE print_climate_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_sum_flux_type(var)

    TYPE(sum_flux_type), INTENT(in) :: var

    write(*,*) 'sumflux%sumpn ', var%sumpn
    write(*,*) 'sumflux%sumrp ', var%sumrp
    write(*,*) 'sumflux%sumrpw ', var%sumrpw
    write(*,*) 'sumflux%sumrpr ', var%sumrpr
    write(*,*) 'sumflux%sumrs ', var%sumrs
    write(*,*) 'sumflux%sumrd ', var%sumrd
    write(*,*) 'sumflux%dsumpn ', var%dsumpn
    write(*,*) 'sumflux%dsumrp ', var%dsumrp
    write(*,*) 'sumflux%dsumrs ', var%dsumrs
    write(*,*) 'sumflux%dsumrd ', var%dsumrd
    write(*,*) 'sumflux%sumxrp ', var%sumxrp
    write(*,*) 'sumflux%sumxrs ', var%sumxrs

  END SUBROUTINE print_sum_flux_type


  ! ------------------------------------------------------------------------------


  SUBROUTINE print_bgc_pool_type(var)

    TYPE(bgc_pool_type), INTENT(in) :: var

    write(*,*) 'bgc%cplant ', var%cplant
    write(*,*) 'bgc%csoil ', var%csoil

  END SUBROUTINE print_bgc_pool_type


  ! ------------------------------------------------------------------------------

END MODULE cable_def_types_mod
