! ==============================================================================
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
module cable_def_types_mod

  ! Contains all variables which are not subroutine-internal

  implicit none

  private

  !--- CABLE default KINDs for representing INTEGER/REAL values
  !--- at least 10-digit precision

  integer, public :: &
       mp,     & ! # total no of patches/tiles
       mvtype, & ! total # vegetation types, from input
       mstype, & ! total # soil types,       from input
       mland     ! # land grid cells

  integer, public, parameter :: &
       ! i_d  = KIND(9), &                  ! probably unintended
       i_d  = SELECTED_INT_kind(9), &       ! as in pop.f90
       ! r_2  = SELECTED_REAL_KIND(12, 50), & ! old
       r_2  = kind(1.0d0), &                ! as in pop.f90
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

  ! defined types
  public :: balances_type
  public :: soil_parameter_type
  public :: soil_snow_type
  public :: veg_parameter_type
  public :: canopy_type
  public :: radiation_type
  public :: roughness_type
  public :: air_type
  public :: met_type
  public :: climate_type
  public :: sum_flux_type
  public :: bgc_pool_type

  ! functions on the defined types
  public :: alloc_cbm_var         ! allocate all variables in a type
  public :: dealloc_cbm_var       ! deallocate all variables in a type
  public :: print_cbm_var         ! print all variables of a type to screen
  public :: read_netcdf_cbm_var   ! read all variables of a type from a netcdf file
  public :: write_netcdf_cbm_var  ! write all variables of a type in a netcdf file
  public :: zero_cbm_var          ! set to zero all variables in a type

  ! general functions
  public :: nc_err

  ! .............................................................................

  ! Energy and water balance variables:
  type balances_type

     real, dimension(:), pointer :: &
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

  end type balances_type

  ! .............................................................................

  ! Soil parameters:
  type soil_parameter_type

     integer, dimension(:), pointer :: &
          isoilm     ! integer soil type

     real, dimension(:), pointer :: &
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

     real(r_2), dimension(:), pointer :: &
          cnsd => null(),    & ! thermal conductivity of dry soil [W/m/K]
          pwb_min => null()    ! working variable (swilt/ssat)**ibp2

     real, dimension(:,:), pointer :: &
          albsoil    ! soil reflectance (2nd dim. BP 21Oct2009)

     ! Additional SLI parameters
     integer,   dimension(:),   pointer :: nhorizons => null() ! number of soil horizons
     integer,   dimension(:,:), pointer :: ishorizon => null() ! horizon number 1:nhorizons
     real(r_2), dimension(:),   pointer :: clitt => null()     ! litter (tC/ha)
     real(r_2), dimension(:),   pointer :: zeta => null()      ! macropore parameter
     real(r_2), dimension(:),   pointer :: fsatmax => null()   ! variably saturated area parameter
     real(r_2), dimension(:,:), pointer :: swilt_vec => null() ! vol H2O @ wilting
     real(r_2), dimension(:,:), pointer :: ssat_vec => null()  ! vol H2O @ sat
     real(r_2), dimension(:,:), pointer :: sfc_vec => null()   ! vol H2O @ fc

  end type soil_parameter_type

  ! .............................................................................

  ! Soil and snow variables:
  type soil_snow_type

     integer, dimension(:), pointer :: isflag => null() ! 0 => no snow 1 => snow

     real, dimension(:), pointer :: &
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

     real, dimension(:,:), pointer :: &
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

     real(r_2), dimension(:), pointer :: &
          wbtot => null()   ! total soil water (mm)

     real(r_2), dimension(:,:), pointer :: &
          gammzz => null(),  & ! heat capacity for each soil layer
          wb => null(),      & ! volumetric soil moisture (solid+liq)
          wbice => null(),   & ! soil ice
          wblf => null(),    & !
          wbfice => null()     !

     ! Additional SLI variables:
     real(r_2), dimension(:,:), pointer :: S => null()         ! moisture content relative to sat value    (edit vh 23/01/08)
     real(r_2), dimension(:,:), pointer :: Tsoil => null()     ! Tsoil (deg C)
     real(r_2), dimension(:),   pointer :: SL => null()        ! litter moisture content relative to sat value (edit vh 23/01/08)
     real(r_2), dimension(:),   pointer :: TL => null()        ! litter temperature in K     (edit vh 23/01/08)
     real(r_2), dimension(:),   pointer :: h0 => null()        ! pond height in m            (edit vh 23/01/08)
     real(r_2), dimension(:,:), pointer :: rex => null()       ! root extraction from each layer (mm/dels)
     real(r_2), dimension(:,:), pointer :: wflux => null()     ! water flux at layer boundaries (mm s-1)
     real(r_2), dimension(:),   pointer :: delwcol => null()   ! change in water column (mm / dels)
     real(r_2), dimension(:),   pointer :: zdelta => null()    ! water table depth           (edit vh 23/06/08)
     real(r_2), dimension(:,:), pointer :: kth => null()       ! thermal conductivity           (edit vh 29/07/08)
     real(r_2), dimension(:),   pointer :: Tsurface => null()  !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
     real(r_2), dimension(:),   pointer :: lE => null()      ! soil latent heat flux
     real(r_2), dimension(:),   pointer :: evap => null()    ! soil evaporation (mm / dels)
     real(r_2), dimension(:,:), pointer :: ciso => null()    ! concentration of minor isotopologue in soil water (kg m-3 water)
     real(r_2), dimension(:),   pointer :: cisoL => null()   ! concentration of minor isotopologue in litter water (kg m-3 water)
     real(r_2), dimension(:),   pointer :: rlitt => null()   ! resistance to heat/moisture transfer through litter (m-1 s)
     real(r_2), dimension(:,:), pointer :: thetai => null()  ! volumetric ice content (MC)
     real(r_2), dimension(:,:), pointer :: snowliq => null() ! liquid snow content (mm H2O)
     real(r_2), dimension(:),   pointer :: nsteps => null()  ! number of iterations at each timestep
     real(r_2), dimension(:),   pointer :: TsurfaceFR => null() ! temperature at surface (soil, pond or litter) (edit vh 22/10/08)
     real(r_2), dimension(:,:), pointer :: Ta_daily => null() ! air temp averaged over last 24h
     integer, dimension(:),     pointer :: nsnow => null() ! number of layers in snow-pack (0-nsnow_max)
     real(r_2), dimension(:),   pointer :: Qadv_daily => null() ! advective heat flux into surface , daily average (W m-2)
     real(r_2), dimension(:),   pointer :: G0_daily => null()  ! conductive heat flux into surface , daily average (W m-2)
     real(r_2), dimension(:),   pointer :: Qevap_daily => null() ! evaporative flux at surface, daily average (m s-1)
     real(r_2), dimension(:),   pointer :: Qprec_daily => null() ! liquid precip, daily average (m s-1)
     real(r_2), dimension(:),   pointer :: Qprec_snow_daily => null() ! solid precip, daily average (m s-1)
     real(r_2), dimension(:),   pointer :: E_fusion_sn => null()
     real(r_2), dimension(:),   pointer :: E_sublimation_sn => null()
     real(r_2), dimension(:),   pointer :: latent_heat_sn => null()
     real(r_2), dimension(:),   pointer :: evap_liq_sn => null()
     real(r_2), dimension(:),   pointer :: surface_melt => null()
     real(r_2), dimension(:),   pointer :: Qadv_rain_sn => null()

  end type soil_snow_type

  ! .............................................................................

  ! Vegetation parameters:
  type veg_parameter_type

     integer, dimension(:), pointer :: &
          iveg => null(),   & ! vegetation type
          ivegp => null(),  & ! dominant potential vegetation type
          iLU => null()        ! land use type

     real, dimension(:), pointer :: &
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

     logical, dimension(:), pointer :: &
          deciduous => null() ! flag used for phenology fix

     real, dimension(:,:), pointer :: &
          refl => null(),    &
          taul => null(),    &
          froot => null()      ! fraction of root in each soil layer

     ! Additional  veg parameters:
     real(r_2), dimension(:), pointer :: rootbeta => null() ! parameter for estimating vertical root mass distribution (froot)
     real(r_2), dimension(:), pointer :: gamma => null()    ! parameter in root efficiency function (Lai and Katul 2000)
     real(r_2), dimension(:), pointer :: ZR => null()       ! maximum rooting depth (cm)
     real(r_2), dimension(:), pointer :: F10 => null()      ! fraction of roots in top 10 cm

     real(r_2), dimension(:), pointer :: clitt => null()     !

     ! Additional POP veg param
     integer,   dimension(:,:), pointer :: disturbance_interval => null()
     real(r_2), dimension(:,:), pointer :: disturbance_intensity => null()

  end type veg_parameter_type

  ! .............................................................................

  ! Canopy/vegetation variables:
  type canopy_type

     real, dimension(:), pointer :: &
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

     real, dimension(:,:), pointer :: &
          evapfbl => null(), &
          gswx => null(),    & ! stom cond for water
          zetar => null(),   & ! stability parameter (ref height)
                                !! vh_js !!
          zetash => null()     ! stability parameter (shear height)

     real(r_2), dimension(:), pointer :: &
          fess => null(),    & ! latent heatfl from soil (W/m2)
          fesp => null(),    & ! latent heatfl from soil (W/m2)
          dgdtg => null(),   & ! derivative of gflux wrt soil temp
          fes => null(),     & ! latent heatfl from soil (W/m2)
          fes_cor => null(), & ! latent heatfl from soil (W/m2)
          fevc => null(),    & ! dry canopy transpiration (W/m2)
          ofes => null(),    & ! latent heatfl from soil (W/m2)
          A_sl => null(),    & ! net photosynthesis from sunlit leaves
          A_sh => null(),    & ! net photosynthesis from shaded leaves
          A_slC => null(),   & ! net photosynthesis from sunlit leaves (rubisco limited)
          A_shC => null(),   & ! net photosynthesis from shaded leaves  (rubisco limited)
          A_slJ => null(),   & ! net photosynthesis from sunlit leaves (rubp limited)
          A_shJ => null(),   & ! net photosynthesis from shaded leaves  (rubp limited)
          GPP_sl => null(), &  ! gross photosynthesis from sunlit leaves
          GPP_sh => null(), &  ! gross photosynthesis from shaded leaves
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
     real(r_2), dimension(:,:),   pointer :: gw => null()     ! dry canopy conductance (ms-1) edit vh 6/7/09
     real(r_2), dimension(:,:,:), pointer :: ancj => null()   ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
     real(r_2), dimension(:,:),   pointer :: tlfy => null()   ! sunlit and shaded leaf temperatures
     real(r_2), dimension(:,:),   pointer :: ecy => null()    ! sunlit and shaded leaf transpiration (dry canopy)
     real(r_2), dimension(:,:),   pointer :: ecx => null()    ! sunlit and shaded leaf latent heat flux
     ! REAL(r_2), DIMENSION(:,:,:), POINTER :: ci => null()     ! intra-cellular CO2 vh 6/7/09
     real(r_2), dimension(:),     pointer :: fwsoil => null() !

     ! vh_js - litter thermal conductivity (Wm-2K-1) and vapour diffusivity (m2s-1)
     real(r_2), dimension(:), pointer :: kthLitt => null()
     real(r_2), dimension(:), pointer :: DvLitt => null()

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

  end type canopy_type

  ! .............................................................................

  ! Radiation variables:
  type radiation_type

     real, dimension(:), pointer :: &
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

     real, dimension(:,:), pointer :: &
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

     real, dimension(:,:,:), pointer :: &
          qcan => null() ! absorbed radiation for canopy (W/m^2)

  end type radiation_type

  ! .............................................................................

  ! Roughness variables:
  type roughness_type

     real, dimension(:), pointer :: &
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
     real, dimension(:), pointer :: &
          coexp => null() ! Extinction coef for wind profile in canopy

     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     real, dimension(:), pointer :: &
          usuh => null() ! Friction velocity/windspeed at canopy height

     real, dimension(:), pointer :: & ! for aerodyn resist. calc.
          term2 => null(), &
          term3 => null(), &
          term5 => null(), &
          term6 => null(), &
          term6a => null()

  end type roughness_type

  ! .............................................................................

  ! Air variables:
  type air_type

     real, dimension(:), pointer :: &
          rho => null(),     & ! dry air density (kg m-3)
          volm => null(),    & ! molar volume (m3 mol-1)
          rlam => null(),    & ! latent heat for water (j/kg)
          qsat => null(),    & ! saturation specific humidity
          epsi => null(),    & ! d(qsat)/dT ((kg/kg)/K)
          visc => null(),    & ! air kinematic viscosity (m2/s)
          psyc => null(),    & ! psychrometric constant
          dsatdk => null(),  & ! d(es)/dT (mb/K)
          cmolar => null()     ! conv. from m/s to mol/m2/s

  end type air_type

  ! .............................................................................

  ! Meterological data:
  type met_type

     integer, dimension(:), pointer :: &
          year => null(),    & ! local time year AD
          moy => null()        ! local time month of year

     real, dimension(:), pointer :: &
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
          rhum => null(),    & ! relative humidity (%)
          fdiff => null()      ! fraction of diffuse radiation
     real, dimension(:,:), pointer :: &
          fsd => null()  ! downward short-wave radiation (W/m2)

  end type met_type

  ! .............................................................................

  ! Climate data:
  type climate_type

     integer :: nyear_average = 20
     integer :: nday_average  = 31
     !      INTEGER, POINTER ::                                                  &
     integer :: &
          nyears, & ! number of years in climate record
          doy ! day of year

     integer, dimension(:), pointer :: &
          chilldays => null(), &  ! length of chilling period (period with T<5deg)
          iveg => null(), &       ! potential vegetation type based on climatic constraints
          biome => null(), &
          GMD => null(), &        ! growing moisture days (== number days since min moisture threshold)
          modis_igbp => null(), & ! IGBP biome classification
          DSLR => null(), &       ! days since last rain
          NDAY_Nesterov => null()

     real, dimension(:), pointer :: &
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

     real, dimension(:,:), pointer :: &
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
  end type climate_type

  ! .............................................................................

  ! Cumulative flux variables:
  type sum_flux_type

     real, dimension(:), pointer :: &
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

  end type sum_flux_type

  ! .............................................................................

  type bgc_pool_type

     real, dimension(:,:), pointer :: &
          cplant => null(),  & ! plant carbon (g C/m2))
          csoil => null()      ! soil carbon (g C/m2)

     real, dimension(ncp)  :: ratecp ! plant carbon rate constant (1/year)
     real, dimension(ncs)  :: ratecs ! soil carbon rate constant (1/year)

  end type bgc_pool_type

  ! .............................................................................

  ! Functions for allocating and deallocating these types
  ! Functions are overloaded, i.e. all called by calling alloc_cbm_var.
  interface alloc_cbm_var
     module procedure alloc_balances_type, &
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
  end interface alloc_cbm_var

  interface dealloc_cbm_var
     module procedure dealloc_balances_type, &
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
  end interface dealloc_cbm_var

  ! Functions for setting to zero all variables in types
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

  ! Functions to print to screen all variables in types
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

  ! Functions to read from netcdf file all variables in types
  interface read_netcdf_cbm_var
     module procedure read_netcdf_climate_type
  end interface read_netcdf_cbm_var

  ! Functions to write to netcdf file all variables in types
  interface write_netcdf_cbm_var
     module procedure write_netcdf_climate_type
  end interface write_netcdf_cbm_var

contains

  subroutine alloc_balances_type(bal, mp)

    type(balances_type), intent(inout) :: bal
    integer, intent(in) :: mp

    allocate(bal%drybal(mp))
    allocate(bal%ebal(mp))
    allocate(bal%ebal_tot(mp))
    allocate(bal%ebal_cncheck(mp))
    allocate(bal%ebal_tot_cncheck(mp))
    allocate(bal%ebaltr(mp))
    allocate(bal%ebal_tottr(mp))
    allocate(bal%evap_tot(mp))
    allocate(bal%osnowd0(mp))
    allocate(bal%precip_tot(mp))
    allocate(bal%rnoff_tot(mp))
    allocate(bal%wbal(mp))
    allocate(bal%wbal_tot(mp))
    allocate(bal%wbtot0(mp))
    allocate(bal%wetbal(mp))
    allocate(bal%cansto0(mp))
    allocate(bal%owbtot(mp))
    allocate(bal%evapc_tot(mp))
    allocate(bal%evaps_tot(mp))
    allocate(bal%rnof1_tot(mp))
    allocate(bal%rnof2_tot(mp))
    allocate(bal%snowdc_tot(mp))
    allocate(bal%wbal_tot1(mp))
    allocate(bal%delwc_tot(mp))
    allocate(bal%qasrf_tot(mp))
    allocate(bal%qfsrf_tot(mp))
    allocate(bal%qssrf_tot(mp))
    allocate(bal%Radbal(mp))
    allocate(bal%EbalSoil(mp))
    allocate(bal%Ebalveg(mp))
    allocate(bal%Radbalsum(mp))

  end subroutine alloc_balances_type

  ! ------------------------------------------------------------------

  subroutine alloc_soil_parameter_type(soil, mp)

    type(soil_parameter_type), intent(inout) :: soil
    integer, intent(in) :: mp

    allocate(soil%isoilm(mp))
    allocate(soil%bch(mp))
    allocate(soil%c3(mp))
    allocate(soil%clay(mp))
    allocate(soil%css(mp))
    allocate(soil%hsbh(mp))
    allocate(soil%hyds(mp))
    allocate(soil%i2bp3(mp))
    allocate(soil%ibp2(mp))
    allocate(soil%rhosoil(mp))
    allocate(soil%sand(mp))
    allocate(soil%sfc(mp))
    allocate(soil%silt(mp))
    allocate(soil%ssat(mp))
    allocate(soil%sucs(mp))
    allocate(soil%swilt(mp))
    allocate(soil%zse(ms))
    allocate(soil%zshh(ms+1))
    allocate(soil%soilcol(mp))
    allocate(soil%albsoilf(mp))
    allocate(soil%cnsd(mp))
    allocate(soil%pwb_min(mp))
    allocate(soil%albsoil(mp, nrb))
    ! Allocate variables for SLI soil model:
    allocate(soil%nhorizons(mp))
    allocate(soil%ishorizon(mp,ms))
    allocate(soil%clitt(mp))
    allocate(soil%zeta(mp))
    allocate(soil%fsatmax(mp))
    allocate(soil%swilt_vec(mp,ms))
    allocate(soil%ssat_vec(mp,ms))
    allocate(soil%sfc_vec(mp,ms))

  end subroutine alloc_soil_parameter_type

  ! ------------------------------------------------------------------

  subroutine alloc_soil_snow_type(ssnow, mp)

    type(soil_snow_type), intent(inout) :: ssnow
    integer, intent(in) :: mp

    allocate(ssnow%isflag(mp))
    allocate(ssnow%iantrct(mp))
    allocate(ssnow%pudsto(mp))
    allocate(ssnow%pudsmx(mp))
    allocate(ssnow%cls(mp))
    allocate(ssnow%dfn_dtg(mp))
    allocate(ssnow%dfh_dtg(mp))
    allocate(ssnow%dfe_ddq(mp))
    allocate(ssnow%ddq_dtg(mp))
    allocate(ssnow%evapsn(mp))
    allocate(ssnow%fwtop(mp))
    allocate(ssnow%fwtop1(mp))
    allocate(ssnow%fwtop2(mp))
    allocate(ssnow%fwtop3(mp))
    allocate(ssnow%osnowd(mp))
    allocate(ssnow%potev(mp))
    allocate(ssnow%runoff(mp))
    allocate(ssnow%rnof1(mp))
    allocate(ssnow%rnof2(mp))
    allocate(ssnow%rtsoil(mp))
    allocate(ssnow%wbtot1(mp))
    allocate(ssnow%wbtot2(mp))
    allocate(ssnow%wb_lake(mp))
    allocate(ssnow%sinfil(mp))
    allocate(ssnow%qstss(mp))
    allocate(ssnow%wetfac(mp))
    allocate(ssnow%owetfac(mp))
    allocate(ssnow%t_snwlr(mp))
    allocate(ssnow%tggav(mp))
    allocate(ssnow%otgg(mp))
    allocate(ssnow%otss(mp))
    allocate(ssnow%otss_0(mp))
    allocate(ssnow%tprecip(mp))
    allocate(ssnow%tevap(mp))
    allocate(ssnow%trnoff(mp))
    allocate(ssnow%totenbal(mp))
    allocate(ssnow%totenbal2(mp))
    allocate(ssnow%fland(mp))
    allocate(ssnow%ifland(mp))
    allocate(ssnow%qasrf(mp))
    allocate(ssnow%qfsrf(mp))
    allocate(ssnow%qssrf(mp))
    allocate(ssnow%snage(mp))
    allocate(ssnow%snowd(mp))
    allocate(ssnow%smelt(mp))
    allocate(ssnow%ssdnn(mp))
    allocate(ssnow%tss(mp))
    allocate(ssnow%tss_p(mp))
    allocate(ssnow%deltss(mp))
    allocate(ssnow%owb1(mp))
    allocate(ssnow%sconds(mp,msn))
    allocate(ssnow%sdepth(mp,msn))
    allocate(ssnow%smass(mp,msn))
    allocate(ssnow%ssdn(mp,msn))
    allocate(ssnow%tgg(mp,ms))
    allocate(ssnow%tggsn(mp,msn))
    allocate(ssnow%dtmlt(mp,3))
    allocate(ssnow%albsoilsn(mp,nrb))
    allocate(ssnow%evapfbl(mp,ms))
    allocate(ssnow%tilefrac(mp,n_tiles))
    allocate(ssnow%wbtot(mp))
    allocate(ssnow%gammzz(mp,ms))
    allocate(ssnow%wb(mp,ms))
    allocate(ssnow%wbice(mp,ms))
    allocate(ssnow%wblf(mp,ms))
    allocate(ssnow%wbfice(mp,ms))
    ! Allocate variables for SLI soil model:
    allocate(ssnow%S(mp,ms))
    allocate(ssnow%Tsoil(mp,ms))
    allocate(ssnow%SL(mp))
    allocate(ssnow%TL(mp))
    allocate(ssnow%h0(mp))
    allocate(ssnow%rex(mp,ms))
    allocate(ssnow%wflux(mp,0:ms))
    allocate(ssnow%delwcol(mp))
    allocate(ssnow%zdelta(mp))
    allocate(ssnow%kth(mp,ms))
    allocate(ssnow%Tsurface(mp))
    allocate(ssnow%lE(mp))
    allocate(ssnow%evap(mp))
    allocate(ssnow%ciso(mp,ms+1))
    allocate(ssnow%cisoL(mp))
    allocate(ssnow%rlitt(mp))
    allocate(ssnow%thetai(mp,ms))
    allocate(ssnow%snowliq(mp,msn))
    allocate(ssnow%nsteps(mp))
    allocate(ssnow%TsurfaceFR(mp))
    allocate(ssnow%Ta_daily(mp,100))
    allocate(ssnow%nsnow(mp))
    allocate(ssnow%Qadv_daily(mp))
    allocate(ssnow%G0_daily(mp))
    allocate(ssnow%Qevap_daily(mp))
    allocate(ssnow%Qprec_daily(mp))
    allocate(ssnow%Qprec_snow_daily(mp))
    allocate(ssnow%E_fusion_sn(mp))
    allocate(ssnow%E_sublimation_sn(mp))
    allocate(ssnow%latent_heat_sn(mp))
    allocate(ssnow%evap_liq_sn(mp))
    allocate(ssnow%surface_melt(mp))
    allocate(ssnow%Qadv_rain_sn(mp))

  end subroutine alloc_soil_snow_type

  ! ------------------------------------------------------------------

  subroutine alloc_veg_parameter_type(veg, mp)

    type(veg_parameter_type), intent(inout) :: veg
    integer, intent(in) :: mp

    allocate(veg%iveg(mp))
    allocate(veg%ivegp(mp))
    allocate(veg%iLU(mp))
    allocate(veg%canst1(mp))
    allocate(veg%dleaf(mp))
    allocate(veg%ejmax(mp))
    allocate(veg%ejmax_shade(mp))
    allocate(veg%ejmax_sun(mp))
    allocate(veg%meth(mp))
    allocate(veg%frac4(mp))
    allocate(veg%hc(mp))
    allocate(veg%vlai(mp))
    allocate(veg%xalbnir(mp))
    allocate(veg%rp20(mp))
    allocate(veg%rpcoef(mp))
    allocate(veg%rs20(mp))
    allocate(veg%shelrb(mp))
    allocate(veg%vegcf(mp))
    allocate(veg%tminvj(mp))
    allocate(veg%toptvj(mp))
    allocate(veg%tmaxvj(mp))
    allocate(veg%vbeta(mp))
    allocate(veg%vcmax(mp))
    allocate(veg%vcmax_shade(mp))
    allocate(veg%vcmax_sun(mp))
    allocate(veg%xfang(mp))
    allocate(veg%extkn(mp))
    allocate(veg%vlaimax(mp))
    allocate(veg%wai(mp))
    allocate(veg%a1gs(mp))
    allocate(veg%d0gs(mp))
    allocate(veg%alpha(mp))
    allocate(veg%convex(mp))
    allocate(veg%cfrd(mp))
    allocate(veg%gswmin(mp))
    allocate(veg%conkc0(mp))
    allocate(veg%conko0(mp))
    allocate(veg%ekc(mp))
    allocate(veg%eko(mp))
    allocate(veg%g0(mp))   ! Ticket #56.
    allocate(veg%g1(mp))   ! Ticket #56.
    allocate(veg%vcmaxcc(mp))
    allocate(veg%ejmaxcc(mp))
    allocate(veg%gmmax(mp))
    allocate(veg%gm(mp))
    allocate(veg%c4kci(mp))
    allocate(veg%c4kcc(mp))
    allocate(veg%bjv(mp))
    allocate(veg%deciduous(mp))
    ! was nrb(=3), but never uses (:,3) in model
    allocate(veg%refl(mp,2)) ! jhan:swb?
    allocate(veg%taul(mp,2))
    allocate(veg%froot(mp,ms))
    allocate(veg%rootbeta(mp))
    allocate(veg%gamma(mp))
    allocate(veg%ZR(mp))
    allocate(veg%F10(mp))
    allocate(veg%clitt(mp))
    allocate(veg%disturbance_interval(mp,2))
    allocate(veg%disturbance_intensity(mp,2))

  end subroutine alloc_veg_parameter_type

  ! ------------------------------------------------------------------

  subroutine alloc_canopy_type(canopy, mp)

    type(canopy_type), intent(inout) :: canopy
    integer, intent(in) :: mp

    allocate(canopy%cansto(mp))
    allocate(canopy%cduv(mp))
    allocate(canopy%delwc(mp))
    allocate(canopy%dewmm(mp))
    allocate(canopy%fe(mp))
    allocate(canopy%fh(mp))
    allocate(canopy%fpn(mp))
    allocate(canopy%frp(mp))
    allocate(canopy%frpw(mp))
    allocate(canopy%frpr(mp))
    allocate(canopy%frs(mp))
    allocate(canopy%fnee(mp))
    allocate(canopy%frday(mp))
    allocate(canopy%fnv(mp))
    allocate(canopy%fev(mp))
    allocate(canopy%epot(mp))
    allocate(canopy%fnpp(mp))
    allocate(canopy%fevw_pot(mp))
    allocate(canopy%gswx_T(mp))
    allocate(canopy%cdtq(mp))
    allocate(canopy%wetfac_cs(mp))
    allocate(canopy%fevw(mp))
    allocate(canopy%fhvw(mp))
    allocate(canopy%oldcansto(mp))
    allocate(canopy%fhv(mp))
    allocate(canopy%fns(mp))
    allocate(canopy%fhs(mp))
    allocate(canopy%fhs_cor(mp))
    allocate(canopy%ga(mp))
    allocate(canopy%ghflux(mp))
    allocate(canopy%precis(mp))
    allocate(canopy%qscrn(mp))
    allocate(canopy%rnet(mp))
    allocate(canopy%rniso(mp))
    allocate(canopy%segg(mp))
    allocate(canopy%sghflux(mp))
    allocate(canopy%through(mp))
    allocate(canopy%spill(mp))
    allocate(canopy%tscrn(mp))
    allocate(canopy%wcint(mp))
    allocate(canopy%tv(mp))
    allocate(canopy%us(mp))
    allocate(canopy%uscrn(mp))
    allocate(canopy%vlaiw(mp))
    allocate(canopy%rghlai(mp))
    allocate(canopy%fwet(mp))
    allocate(canopy%evapfbl(mp,ms))
    allocate(canopy%gswx(mp,mf))
    allocate(canopy%zetar(mp,NITER))
    allocate(canopy%zetash(mp,NITER))
    allocate(canopy%fess(mp))
    allocate(canopy%fesp(mp))
    allocate(canopy%dgdtg(mp))
    allocate(canopy%fes(mp))
    allocate(canopy%fes_cor(mp))
    allocate(canopy%fevc(mp))
    allocate(canopy%ofes(mp))
    allocate(canopy%A_sl(mp))
    allocate(canopy%A_sh(mp))
    allocate(canopy%A_slC(mp))
    allocate(canopy%A_shC(mp))
    allocate(canopy%A_slJ(mp))
    allocate(canopy%A_shJ(mp))
    allocate(canopy%GPP_sl(mp))
    allocate(canopy%GPP_sh(mp))
    allocate(canopy%fevc_sl(mp))
    allocate(canopy%fevc_sh(mp))
    allocate(canopy%eta_A_cs(mp))
    allocate(canopy%dAdcs(mp))
    allocate(canopy%eta_GPP_cs(mp))
    allocate(canopy%eta_A_cs_sl(mp))
    allocate(canopy%eta_A_cs_sh(mp))
    allocate(canopy%eta_fevc_cs_sl(mp))
    allocate(canopy%eta_fevc_cs_sh(mp))
    allocate(canopy%eta_fevc_cs(mp))
    allocate(canopy%cs(mp))
    allocate(canopy%cs_sl(mp))
    allocate(canopy%cs_sh(mp))
    ! allocate(canopy%ci_sl(mp))
    ! allocate(canopy%ci_sh(mp))
    allocate(canopy%tlf(mp))
    allocate(canopy%dlf(mp))
    allocate(canopy%gw(mp,mf))     ! dry canopy conductance (ms-1) edit vh 6/7/09
    allocate(canopy%ancj(mp,mf,3)) ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
    allocate(canopy%tlfy(mp,mf))   ! sunlit and shaded leaf temperatures
    allocate(canopy%ecy(mp,mf))    ! sunlit and shaded leaf transpiration (dry canopy)
    allocate(canopy%ecx(mp,mf))    ! sunlit and shaded leaf latent heat flux
    ! allocate(canopy%ci(mp,mf,3))   ! intra-cellular CO2 vh 6/7/09
    allocate(canopy%fwsoil(mp))
    ! vh_js - litter resistances to heat and vapour transfer
    allocate(canopy%kthLitt(mp))
    allocate(canopy%DvLitt(mp))
    ! 13C
    allocate(canopy%An(mp,mf))
    allocate(canopy%Rd(mp,mf))
    allocate(canopy%isc3(mp))
    allocate(canopy%vcmax(mp,mf))
    allocate(canopy%gammastar(mp,mf))
    allocate(canopy%gsc(mp,mf))
    allocate(canopy%gbc(mp,mf))
    allocate(canopy%gac(mp,mf))
    allocate(canopy%ci(mp,mf))

  end subroutine alloc_canopy_type

  ! ------------------------------------------------------------------

  subroutine alloc_radiation_type(rad, mp)

    type(radiation_type), intent(inout) :: rad
    integer, intent(in) :: mp

    allocate(rad%transb(mp))
    allocate(rad%albedo_T(mp))
    allocate(rad%longitude(mp))
    allocate(rad%workp1(mp))
    allocate(rad%workp2(mp))
    allocate(rad%workp3(mp))
    allocate(rad%extkb(mp))
    allocate(rad%extkd2(mp))
    allocate(rad%extkd(mp))
    allocate(rad%flws(mp))
    allocate(rad%latitude(mp))
    allocate(rad%lwabv(mp))
    allocate(rad%qssabs(mp))
    allocate(rad%transd(mp))
    allocate(rad%trad(mp))
    allocate(rad%fvlai(mp,mf))
    allocate(rad%rhocdf(mp,nrb))
    allocate(rad%rniso(mp,mf))
    allocate(rad%scalex(mp,mf))
    allocate(rad%albedo(mp,nrb))
    allocate(rad%reffdf(mp,nrb))
    allocate(rad%reffbm(mp,nrb))
    allocate(rad%extkbm(mp,nrb))
    allocate(rad%extkdm(mp,nrb))
    allocate(rad%fbeam(mp,nrb))
    allocate(rad%cexpkbm(mp,swb))
    allocate(rad%cexpkdm(mp,swb))
    allocate(rad%rhocbm(mp,nrb))
    allocate(rad%gradis(mp,mf))
    allocate(rad%qcan(mp,mf,nrb))

  end subroutine alloc_radiation_type

  ! ------------------------------------------------------------------

  subroutine alloc_roughness_type(rough, mp)

    type(roughness_type), intent(inout) :: rough
    integer, intent(in) :: mp

    allocate(rough%disp(mp))
    allocate(rough%hruff(mp))
    allocate(rough%hruff_grmx(mp))
    allocate(rough%rt0us(mp))
    allocate(rough%rt1usa(mp))
    allocate(rough%rt1usb(mp))
    allocate(rough%rt1(mp))
    allocate(rough%za_uv(mp))
    allocate(rough%za_tq(mp))
    allocate(rough%z0m(mp))
    allocate(rough%zref_uv(mp))
    allocate(rough%zref_tq(mp))
    allocate(rough%zruffs(mp))
    allocate(rough%z0soilsn(mp))
    allocate(rough%z0soil(mp))
    allocate(rough%coexp(mp))
    allocate(rough%usuh(mp))
    allocate(rough%term2(mp))
    allocate(rough%term3(mp))
    allocate(rough%term5(mp))
    allocate(rough%term6(mp))
    allocate(rough%term6a(mp))

  end subroutine alloc_roughness_type

  ! ------------------------------------------------------------------

  subroutine alloc_air_type(air, mp)

    type(air_type), intent(inout) :: air
    integer, intent(in) :: mp

    allocate(air%rho(mp))
    allocate(air%volm(mp))
    allocate(air%rlam(mp))
    allocate(air%qsat(mp))
    allocate(air%epsi(mp))
    allocate(air%visc(mp))
    allocate(air%psyc(mp))
    allocate(air%dsatdk(mp))
    allocate(air%cmolar(mp))

  end subroutine alloc_air_type

  ! ------------------------------------------------------------------

  subroutine alloc_met_type(met, mp)

    type(met_type), intent(inout) :: met
    integer, intent(in) :: mp

    allocate(met%year(mp))
    allocate(met%moy(mp))
    allocate(met%ca(mp))
    allocate(met%doy(mp))
    allocate(met%hod(mp))
    allocate(met%ofsd(mp))
    allocate(met%fld(mp))
    allocate(met%precip(mp))
    allocate(met%precip_sn(mp))
    allocate(met%tk(mp))
    allocate(met%tvair(mp))
    allocate(met%tvrad(mp))
    allocate(met%pmb(mp))
    allocate(met%ua(mp))
    allocate(met%qv(mp))
    allocate(met%qvair(mp))
    allocate(met%da(mp))
    allocate(met%dva(mp))
    allocate(met%coszen(mp))
    allocate(met%Ndep(mp))
    allocate(met%Pdep(mp))
    allocate(met%u10(mp))
    allocate(met%rhum(mp))
    allocate(met%fdiff(mp))
    allocate(met%fsd(mp,swb))

  end subroutine alloc_met_type

  ! ------------------------------------------------------------------

  subroutine alloc_climate_type(climate, mp, ktauday)

    implicit none

    type(climate_type), intent(inout) :: climate
    integer,            intent(in)    :: mp
    integer,            intent(in)    :: ktauday

    integer :: ny, nd

    ny = climate%nyear_average
    nd = climate%nday_average

    allocate(climate%chilldays(mp))
    allocate(climate%iveg(mp))
    allocate(climate%biome(mp))
    allocate(climate%GMD(mp))
    allocate(climate%modis_igbp(mp))
    allocate(climate%DSLR(mp))
    allocate(climate%NDAY_Nesterov(mp))
    allocate(climate%dtemp(mp))
    allocate(climate%dmoist(mp))
    allocate(climate%dmoist_min(mp))
    allocate(climate%dmoist_min20(mp))
    allocate(climate%dmoist_max(mp))
    allocate(climate%dmoist_max20(mp))
    allocate(climate%mtemp(mp))
    allocate(climate%qtemp(mp))
    allocate(climate%mmoist(mp))
    allocate(climate%mtemp_min(mp))
    allocate(climate%mtemp_max(mp))
    allocate(climate%qtemp_max(mp))
    allocate(climate%qtemp_max_last_year(mp))
    allocate(climate%mtemp_min20(mp))
    allocate(climate%mtemp_max20(mp))
    allocate(climate%atemp_mean(mp))
    allocate(climate%AGDD5(mp))
    allocate(climate%GDD5(mp))
    allocate(climate%AGDD0(mp))
    allocate(climate%GDD0(mp))
    allocate(climate%alpha_PT(mp))
    allocate(climate%evap_PT(mp))
    allocate(climate%aevap(mp))
    allocate(climate%alpha_PT20(mp))
    allocate(climate%GDD0_rec(mp))
    allocate(climate%frec(mp))
    allocate(climate%dtemp_min(mp))
    allocate(climate%fdorm(mp))
    allocate(climate%fapar_ann_max(mp))
    allocate(climate%fapar_ann_max_last_year(mp))
    allocate(climate%AvgAnnMaxFAPAR(mp))
    allocate(climate%dtemp_max(mp))
    allocate(climate%drhum(mp))
    allocate(climate%du10_max(mp))
    allocate(climate%dprecip(mp))
    allocate(climate%aprecip(mp))
    allocate(climate%aprecip_av20(mp))
    allocate(climate%last_precip(mp))
    allocate(climate%KBDI(mp))
    allocate(climate%FFDI(mp))
    allocate(climate%D_MacArthur(mp))
    allocate(climate%Nesterov_Current(mp))
    allocate(climate%Nesterov_ann_max(mp))
    allocate(climate%Nesterov_ann_max_last_year(mp))
    allocate(climate%Nesterov_ann_running_max(mp))
    allocate(climate%mtemp_min_20(mp, ny))
    allocate(climate%mtemp_max_20(mp, ny))
    allocate(climate%dmoist_min_20(mp, ny))
    allocate(climate%dmoist_max_20(mp, ny))
    allocate(climate%dtemp_31(mp, nd))
    allocate(climate%dmoist_31(mp, nd))
    allocate(climate%alpha_PT_20(mp, ny))
    allocate(climate%dtemp_91(mp, 91))
    allocate(climate%APAR_leaf_sun(mp, ktauday*5))
    allocate(climate%APAR_leaf_shade(mp, ktauday*5))
    allocate(climate%Dleaf_sun(mp, ktauday*5))
    allocate(climate%Dleaf_shade(mp, ktauday*5))
    allocate(climate%Tleaf_sun(mp, ktauday*5))
    allocate(climate%Tleaf_shade(mp, ktauday*5))
    allocate(climate%cs_sun(mp, ktauday*5))
    allocate(climate%cs_shade(mp, ktauday*5))
    allocate(climate%scalex_sun(mp, ktauday*5))
    allocate(climate%scalex_shade(mp, ktauday*5))
    allocate(climate%fwsoil(mp, ktauday*5))
    allocate(climate%aprecip_20(mp, ny))
    allocate(climate%Rd_sun(mp, ktauday*5))
    allocate(climate%Rd_shade(mp, ktauday*5))

  end subroutine alloc_climate_type

  ! ------------------------------------------------------------------

  subroutine alloc_sum_flux_type(sum_flux, mp)

    type(sum_flux_type), intent(inout) :: sum_flux
    integer, intent(in) :: mp

    allocate(sum_flux%sumpn(mp))
    allocate(sum_flux%sumrp(mp))
    allocate(sum_flux%sumrpw(mp))
    allocate(sum_flux%sumrpr(mp))
    allocate(sum_flux%sumrs(mp))
    allocate(sum_flux%sumrd(mp))
    allocate(sum_flux%dsumpn(mp))
    allocate(sum_flux%dsumrp(mp))
    allocate(sum_flux%dsumrs(mp))
    allocate(sum_flux%dsumrd(mp))
    allocate(sum_flux%sumxrp(mp))
    allocate(sum_flux%sumxrs(mp))

  end subroutine alloc_sum_flux_type

  ! ------------------------------------------------------------------

  subroutine alloc_bgc_pool_type(bgc, mp)

    type(bgc_pool_type), intent(inout) :: bgc
    integer, intent(in) :: mp

    allocate(bgc%cplant(mp,ncp))
    allocate(bgc%csoil(mp,ncs))

  end subroutine alloc_bgc_pool_type

  ! ------------------------------------------------------------------

  ! Begin deallocation routines:
  subroutine dealloc_balances_type(bal)

    type(balances_type), intent(inout) :: bal

    deallocate(bal%drybal)
    deallocate(bal%ebal)
    deallocate(bal%ebal_tot)
    deallocate(bal%ebaltr)
    deallocate(bal%ebal_tottr)
    deallocate(bal%ebal_cncheck)
    deallocate(bal%ebal_tot_cncheck)
    deallocate(bal%evap_tot)
    deallocate(bal%osnowd0)
    deallocate(bal%precip_tot)
    deallocate(bal%rnoff_tot)
    deallocate(bal%wbal)
    deallocate(bal%wbal_tot)
    deallocate(bal%wbtot0)
    deallocate(bal%wetbal)
    deallocate(bal%cansto0)
    deallocate(bal%evapc_tot)
    deallocate(bal%evaps_tot)
    deallocate(bal%rnof1_tot)
    deallocate(bal%rnof2_tot)
    deallocate(bal%snowdc_tot)
    deallocate(bal%wbal_tot1)
    deallocate(bal%owbtot)
    deallocate(bal%delwc_tot)
    deallocate(bal%qasrf_tot)
    deallocate(bal%qfsrf_tot)
    deallocate(bal%qssrf_tot)
    deallocate(bal%Radbal)
    deallocate(bal%Ebalsoil)
    deallocate(bal%Ebalveg)
    deallocate(bal%Radbalsum)

  end subroutine dealloc_balances_type

  ! ------------------------------------------------------------------

  subroutine dealloc_soil_parameter_type(soil)

    type(soil_parameter_type), intent(inout) :: soil

    deallocate(soil%bch)
    deallocate(soil%c3)
    deallocate(soil%clay)
    deallocate(soil%css)
    deallocate(soil%hsbh)
    deallocate(soil%hyds)
    deallocate(soil%i2bp3)
    deallocate(soil%ibp2)
    deallocate(soil%isoilm)
    deallocate(soil%rhosoil)
    deallocate(soil%sand)
    deallocate(soil%sfc)
    deallocate(soil%silt)
    deallocate(soil%ssat)
    deallocate(soil%sucs)
    deallocate(soil%swilt)
    deallocate(soil%zse)
    deallocate(soil%zshh)
    deallocate(soil%cnsd)
    deallocate(soil%albsoil)
    deallocate(soil%cnsd)
    deallocate(soil%pwb_min)
    deallocate(soil%albsoilf)
    deallocate(soil%soilcol)
    ! Deallocate variables for SLI soil model:
    deallocate(soil%nhorizons)
    deallocate(soil%ishorizon)
    deallocate(soil%clitt)
    deallocate(soil%zeta)
    deallocate(soil%fsatmax)
    deallocate(soil%swilt_vec)
    deallocate(soil%ssat_vec)
    deallocate(soil%sfc_vec)
    if (associated(soil%swilt_vec)) deallocate(soil%swilt_vec)
    if (associated(soil%ssat_vec))  deallocate(soil%ssat_vec)
    if (associated(soil%sfc_vec))   deallocate(soil%sfc_vec)

  end subroutine dealloc_soil_parameter_type

  ! ------------------------------------------------------------------

  subroutine dealloc_soil_snow_type(ssnow)

    type(soil_snow_type), intent(inout) :: ssnow

    deallocate(ssnow%iantrct)
    deallocate(ssnow%pudsto)
    deallocate(ssnow%pudsmx)
    deallocate(ssnow%dtmlt)
    deallocate(ssnow%albsoilsn)
    deallocate(ssnow%cls)
    deallocate(ssnow%dfn_dtg)
    deallocate(ssnow%dfh_dtg)
    deallocate(ssnow%dfe_ddq)
    deallocate(ssnow%ddq_dtg)
    deallocate(ssnow%evapsn)
    deallocate(ssnow%fwtop)
    deallocate(ssnow%fwtop1)
    deallocate(ssnow%fwtop2)
    deallocate(ssnow%fwtop3)
    deallocate(ssnow%gammzz)
    deallocate(ssnow%isflag)
    deallocate(ssnow%osnowd)
    deallocate(ssnow%potev)
    deallocate(ssnow%runoff)
    deallocate(ssnow%rnof1)
    deallocate(ssnow%rnof2)
    deallocate(ssnow%rtsoil)
    deallocate(ssnow%sconds)
    deallocate(ssnow%sdepth)
    deallocate(ssnow%smass)
    deallocate(ssnow%snage)
    deallocate(ssnow%snowd)
    deallocate(ssnow%smelt)
    deallocate(ssnow%ssdn)
    deallocate(ssnow%ssdnn)
    deallocate(ssnow%tgg)
    deallocate(ssnow%tggsn)
    deallocate(ssnow%tss)
    deallocate(ssnow%tss_p)
    deallocate(ssnow%deltss)
    deallocate(ssnow%owb1)
    deallocate(ssnow%wb)
    deallocate(ssnow%wbice)
    deallocate(ssnow%wblf)
    deallocate(ssnow%wbtot)
    deallocate(ssnow%wbtot1)
    deallocate(ssnow%wbtot2)
    deallocate(ssnow%wb_lake)
    deallocate(ssnow%sinfil)
    deallocate(ssnow%evapfbl)
    deallocate(ssnow%qstss)
    deallocate(ssnow%wetfac)
    deallocate(ssnow%owetfac)
    deallocate(ssnow%t_snwlr)
    deallocate(ssnow%wbfice)
    deallocate(ssnow%tggav)
    deallocate(ssnow%otgg)
    deallocate(ssnow%otss)
    deallocate(ssnow%otss_0)
    deallocate(ssnow%tprecip)
    deallocate(ssnow%tevap)
    deallocate(ssnow%trnoff)
    deallocate(ssnow%totenbal)
    deallocate(ssnow%totenbal2)
    deallocate(ssnow%fland)
    deallocate(ssnow%ifland)
    deallocate(ssnow%tilefrac)
    deallocate(ssnow%qasrf)
    deallocate(ssnow%qfsrf)
    deallocate(ssnow%qssrf)
    deallocate(ssnow%S)
    deallocate(ssnow%Tsoil)
    deallocate(ssnow%SL)
    deallocate(ssnow%TL)
    deallocate(ssnow%h0)
    deallocate(ssnow%rex)
    deallocate(ssnow%wflux)
    deallocate(ssnow%delwcol)
    deallocate(ssnow%zdelta)
    deallocate(ssnow%kth)
    deallocate(ssnow%Tsurface)
    deallocate(ssnow%lE)
    deallocate(ssnow%evap)
    deallocate(ssnow%ciso)
    deallocate(ssnow%cisoL)
    deallocate(ssnow%rlitt)
    deallocate(ssnow%thetai)
    deallocate(ssnow%snowliq)
    deallocate(ssnow%nsteps)
    deallocate(ssnow%nsnow)
    deallocate(ssnow%TsurfaceFR)
    deallocate(ssnow%Ta_daily)
    deallocate(ssnow%G0_daily)
    deallocate(ssnow%Qadv_daily)
    deallocate(ssnow%Qevap_daily)
    deallocate(ssnow%Qprec_daily)
    deallocate(ssnow%Qprec_snow_daily)
    deallocate(ssnow%E_fusion_sn)
    deallocate(ssnow%E_sublimation_sn)
    deallocate(ssnow%latent_heat_sn)
    deallocate(ssnow%evap_liq_sn)
    deallocate(ssnow%surface_melt)
    deallocate(ssnow%Qadv_rain_sn)

  end subroutine dealloc_soil_snow_type

  ! ------------------------------------------------------------------

  subroutine dealloc_veg_parameter_type(veg)

    type(veg_parameter_type), intent(inout) :: veg

    deallocate(veg%canst1)
    deallocate(veg%dleaf)
    deallocate(veg%ejmax)
    deallocate(veg%ejmax_shade)
    deallocate(veg%ejmax_sun)
    deallocate(veg%iveg)
    deallocate(veg%ivegp)
    deallocate(veg%iLU)
    deallocate(veg%meth)
    deallocate(veg%frac4)
    deallocate(veg%hc)
    deallocate(veg%vlai)
    deallocate(veg%xalbnir)
    deallocate(veg%rp20)
    deallocate(veg%rpcoef)
    deallocate(veg%rs20)
    deallocate(veg%shelrb)
    deallocate(veg%vegcf)
    deallocate(veg%tminvj)
    deallocate(veg%toptvj)
    deallocate(veg%tmaxvj)
    deallocate(veg%vbeta)
    deallocate(veg%vcmax)
    deallocate(veg%vcmax_shade)
    deallocate(veg%vcmax_sun)
    deallocate(veg%xfang)
    deallocate(veg%extkn)
    deallocate(veg%wai)
    deallocate(veg%deciduous)
    deallocate(veg%froot)
    deallocate(veg%refl)
    deallocate(veg%taul)
    deallocate(veg%a1gs)
    deallocate(veg%d0gs)
    deallocate(veg%alpha)
    deallocate(veg%convex)
    deallocate(veg%cfrd)
    deallocate(veg%gswmin)
    deallocate(veg%conkc0)
    deallocate(veg%conko0)
    deallocate(veg%ekc)
    deallocate(veg%eko)
    deallocate(veg%g0) ! Ticket #56.
    deallocate(veg%g1) ! Ticket #56.
    ! Deallocate variables for SLI soil model:
    deallocate(veg%rootbeta)
    deallocate(veg%gamma) ! vh 20/07/09
    deallocate(veg%F10)
    deallocate(veg%ZR)
    deallocate(veg%CLitt)
    deallocate(veg%disturbance_interval)
    deallocate(veg%disturbance_intensity)

  end subroutine dealloc_veg_parameter_type

  ! ------------------------------------------------------------------

  subroutine dealloc_canopy_type(canopy)

    type(canopy_type), intent(inout) :: canopy

    deallocate(canopy%fess)
    deallocate(canopy%fesp)
    deallocate(canopy%cansto)
    deallocate(canopy%cduv)
    deallocate(canopy%delwc)
    deallocate(canopy%dewmm)
    deallocate(canopy%dgdtg)
    deallocate(canopy%fe)
    deallocate(canopy%fh)
    deallocate(canopy%fpn)
    deallocate(canopy%frp)
    deallocate(canopy%frpw)
    deallocate(canopy%frpr)
    deallocate(canopy%frs)
    deallocate(canopy%fnee)
    deallocate(canopy%frday)
    deallocate(canopy%fnv)
    deallocate(canopy%fev)
    deallocate(canopy%fevc)
    deallocate(canopy%fhv)
    deallocate(canopy%fns)
    deallocate(canopy%fhs)
    deallocate(canopy%fhs_cor)
    deallocate(canopy%ga)
    deallocate(canopy%ghflux)
    deallocate(canopy%precis)
    deallocate(canopy%qscrn)
    deallocate(canopy%rnet)
    deallocate(canopy%rniso)
    deallocate(canopy%segg)
    deallocate(canopy%sghflux)
    deallocate(canopy%through)
    deallocate(canopy%spill)
    deallocate(canopy%tscrn)
    deallocate(canopy%wcint)
    deallocate(canopy%tv)
    deallocate(canopy%us)
    deallocate(canopy%uscrn)
    deallocate(canopy%rghlai)
    deallocate(canopy%vlaiw)
    deallocate(canopy%fwet)
    deallocate(canopy%A_sh)
    deallocate(canopy%A_sl)
    deallocate(canopy%A_slC)
    deallocate(canopy%A_shC)
    deallocate(canopy%A_slJ)
    deallocate(canopy%A_shJ)
    deallocate(canopy%eta_A_cs)
    deallocate(canopy%cs)
    deallocate(canopy%dAdcs)
    deallocate(canopy%cs_sl)
    deallocate(canopy%cs_sh)
    deallocate(canopy%tlf)
    deallocate(canopy%dlf)
    deallocate(canopy%evapfbl)
    deallocate(canopy%epot)
    deallocate(canopy%fnpp)
    deallocate(canopy%fevw_pot)
    deallocate(canopy%gswx_T)
    deallocate(canopy%cdtq)
    deallocate(canopy%wetfac_cs)
    deallocate(canopy%fevw)
    deallocate(canopy%fhvw)
    deallocate(canopy%fes)
    deallocate(canopy%fes_cor)
    deallocate(canopy%gswx)
    deallocate(canopy%oldcansto)
    deallocate(canopy%zetar)
    deallocate(canopy%zetash)
    deallocate(canopy%fwsoil)
    deallocate(canopy%ofes)
    !! vh_js !! litter resistances to heat and vapour transfer
    deallocate(canopy%kthLitt)
    deallocate(canopy%DvLitt)

  end subroutine dealloc_canopy_type

  ! ------------------------------------------------------------------

  subroutine dealloc_radiation_type(rad)

    type(radiation_type), intent(inout) :: rad

    deallocate(rad%albedo)
    deallocate(rad%extkb)
    deallocate(rad%extkd2)
    deallocate(rad%extkd)
    deallocate(rad%flws)
    deallocate(rad%fvlai)
    deallocate(rad%latitude)
    deallocate(rad%lwabv)
    deallocate(rad%qcan)
    deallocate(rad%qssabs)
    deallocate(rad%rhocdf)
    deallocate(rad%rniso)
    deallocate(rad%scalex)
    deallocate(rad%transd)
    deallocate(rad%trad)
    deallocate(rad%reffdf)
    deallocate(rad%reffbm)
    deallocate(rad%extkbm)
    deallocate(rad%extkdm)
    deallocate(rad%fbeam)
    deallocate(rad%cexpkbm)
    deallocate(rad%cexpkdm)
    deallocate(rad%rhocbm)
    deallocate(rad%transb)
    deallocate(rad%albedo_T)
    deallocate(rad%gradis)
    deallocate(rad%longitude)
    deallocate(rad%workp1)
    deallocate(rad%workp2)
    deallocate(rad%workp3)

  end subroutine dealloc_radiation_type

  ! ------------------------------------------------------------------

  subroutine dealloc_roughness_type(rough)

    type(roughness_type), intent(inout) :: rough

    deallocate(rough%coexp)
    deallocate(rough%disp)
    deallocate(rough%hruff)
    deallocate(rough%hruff_grmx)
    deallocate(rough%rt0us)
    deallocate(rough%rt1usa)
    deallocate(rough%rt1usb)
    deallocate(rough%rt1)
    deallocate(rough%term2)
    deallocate(rough%term3)
    deallocate(rough%term5)
    deallocate(rough%term6)
    deallocate(rough%term6a)
    deallocate(rough%usuh)
    deallocate(rough%za_uv)
    deallocate(rough%za_tq)
    deallocate(rough%z0m)
    deallocate(rough%zref_uv)
    deallocate(rough%zref_tq)
    deallocate(rough%zruffs)
    deallocate(rough%z0soilsn)
    deallocate(rough%z0soil)

  end subroutine dealloc_roughness_type

  ! ------------------------------------------------------------------

  subroutine dealloc_air_type(air)

    type(air_type), intent(inout) :: air

    deallocate(air%rho)
    deallocate(air%volm)
    deallocate(air%rlam)
    deallocate(air%qsat)
    deallocate(air%epsi)
    deallocate(air%visc)
    deallocate(air%psyc)
    deallocate(air%dsatdk)
    deallocate(air%cmolar)

  end subroutine dealloc_air_type

  ! ------------------------------------------------------------------

  subroutine dealloc_met_type(met)

    type(met_type), intent(inout) :: met

    deallocate(met%ca)
    deallocate(met%year)
    deallocate(met%moy)
    deallocate(met%doy)
    deallocate(met%hod)
    deallocate(met%fsd)
    deallocate(met%ofsd)
    deallocate(met%fld)
    deallocate(met%precip)
    deallocate(met%precip_sn)
    deallocate(met%tk)
    deallocate(met%tvair)
    deallocate(met%tvrad)
    deallocate(met%pmb)
    deallocate(met%ua)
    deallocate(met%qv)
    deallocate(met%qvair)
    deallocate(met%da)
    deallocate(met%dva)
    deallocate(met%coszen)
    deallocate(met%Ndep)
    deallocate(met%Pdep)
    deallocate(met%rhum)
    deallocate(met%fdiff)
    deallocate(met%u10)

  end subroutine dealloc_met_type

  ! ------------------------------------------------------------------

  subroutine dealloc_sum_flux_type(sum_flux)

    type(sum_flux_type), intent(inout) :: sum_flux

    deallocate(sum_flux%sumpn)
    deallocate(sum_flux%sumrp)
    deallocate(sum_flux%sumrpw)
    deallocate(sum_flux%sumrpr)
    deallocate(sum_flux%sumrs)
    deallocate(sum_flux%sumrd)
    deallocate(sum_flux%dsumpn)
    deallocate(sum_flux%dsumrp)
    deallocate(sum_flux%dsumrs)
    deallocate(sum_flux%dsumrd)
    deallocate(sum_flux%sumxrp)
    deallocate(sum_flux%sumxrs)

  end subroutine dealloc_sum_flux_type

  ! ------------------------------------------------------------------

  subroutine dealloc_bgc_pool_type(bgc)

    type(bgc_pool_type), intent(inout) :: bgc

    deallocate(bgc%cplant)
    deallocate(bgc%csoil)

  end subroutine dealloc_bgc_pool_type


  ! ------------------------------------------------------------------


  subroutine zero_balances_type(bal)

    type(balances_type), intent(inout) :: bal

    bal%drybal           = 0
    bal%ebal             = 0
    bal%ebal_tot         = 0
    bal%ebaltr           = 0
    bal%ebal_tottr       = 0
    bal%ebal_cncheck     = 0
    bal%ebal_tot_cncheck = 0
    bal%evap_tot         = 0
    bal%osnowd0          = 0
    bal%precip_tot       = 0
    bal%rnoff_tot        = 0
    bal%wbal             = 0
    bal%wbal_tot         = 0
    bal%wbtot0           = 0
    bal%wetbal           = 0
    bal%cansto0          = 0
    bal%evapc_tot        = 0
    bal%evaps_tot        = 0
    bal%rnof1_tot        = 0
    bal%rnof2_tot        = 0
    bal%snowdc_tot       = 0
    bal%wbal_tot1        = 0
    bal%owbtot           = 0
    bal%delwc_tot        = 0
    bal%qasrf_tot        = 0
    bal%qfsrf_tot        = 0
    bal%qssrf_tot        = 0
    bal%Radbal    = 0
    bal%EbalSoil  = 0
    bal%Ebalveg   = 0
    bal%Radbalsum = 0

  end subroutine zero_balances_type


  ! ------------------------------------------------------------------


  subroutine zero_soil_parameter_type(soil)

    type(soil_parameter_type), intent(inout) :: soil

    soil%bch      = 0
    soil%c3       = 0
    soil%clay     = 0
    soil%css      = 0
    soil%hsbh     = 0
    soil%hyds     = 0
    soil%i2bp3    = 0
    soil%ibp2     = 0
    soil%isoilm   = 0
    soil%rhosoil  = 0
    soil%sand     = 0
    soil%sfc      = 0
    soil%silt     = 0
    soil%ssat     = 0
    soil%sucs     = 0
    soil%swilt    = 0
    soil%zse      = 0
    soil%zshh     = 0
    soil%cnsd     = 0
    soil%albsoil  = 0
    soil%pwb_min  = 0
    soil%albsoilf = 0
    soil%soilcol  = 0
    soil%nhorizons = 0
    soil%ishorizon = 0
    soil%clitt     = 0
    soil%zeta      = 0
    soil%fsatmax   = 0
    soil%swilt_vec = 0
    soil%ssat_vec  = 0
    soil%sfc_vec   = 0
    soil%swilt_vec = 0
    soil%ssat_vec  = 0
    soil%sfc_vec   = 0

  end subroutine zero_soil_parameter_type


  ! ------------------------------------------------------------------


  subroutine zero_soil_snow_type(ssnow)

    type(soil_snow_type), intent(inout) :: ssnow

    ssnow%iantrct   = 0
    ssnow%pudsto    = 0
    ssnow%pudsmx    = 0
    ssnow%dtmlt     = 0
    ssnow%albsoilsn = 0
    ssnow%cls       = 0
    ssnow%dfn_dtg   = 0
    ssnow%dfh_dtg   = 0
    ssnow%dfe_ddq   = 0
    ssnow%ddq_dtg   = 0
    ssnow%evapsn    = 0
    ssnow%fwtop     = 0
    ssnow%fwtop1    = 0
    ssnow%fwtop2    = 0
    ssnow%fwtop3    = 0
    ssnow%gammzz    = 0
    ssnow%isflag    = 0
    ssnow%osnowd    = 0
    ssnow%potev     = 0
    ssnow%runoff    = 0
    ssnow%rnof1     = 0
    ssnow%rnof2     = 0
    ssnow%rtsoil    = 0
    ssnow%sconds    = 0
    ssnow%sdepth    = 0
    ssnow%smass     = 0
    ssnow%snage     = 0
    ssnow%snowd     = 0
    ssnow%smelt     = 0
    ssnow%ssdn      = 0
    ssnow%ssdnn     = 0
    ssnow%tgg       = 0
    ssnow%tggsn     = 0
    ssnow%tss       = 0
    ssnow%tss_p     = 0
    ssnow%deltss    = 0
    ssnow%owb1      = 0
    ssnow%wb        = 0
    ssnow%wbice     = 0
    ssnow%wblf      = 0
    ssnow%wbtot     = 0
    ssnow%wbtot1    = 0
    ssnow%wbtot2    = 0
    ssnow%wb_lake   = 0
    ssnow%sinfil    = 0
    ssnow%evapfbl   = 0
    ssnow%qstss     = 0
    ssnow%wetfac    = 0
    ssnow%owetfac   = 0
    ssnow%t_snwlr   = 0
    ssnow%wbfice    = 0
    ssnow%tggav     = 0
    ssnow%otgg      = 0
    ssnow%otss      = 0
    ssnow%otss_0    = 0
    ssnow%tprecip   = 0
    ssnow%tevap     = 0
    ssnow%trnoff    = 0
    ssnow%totenbal  = 0
    ssnow%totenbal2 = 0
    ssnow%fland     = 0
    ssnow%ifland    = 0
    ssnow%tilefrac  = 0
    ssnow%qasrf     = 0
    ssnow%qfsrf     = 0
    ssnow%qssrf     = 0
    ssnow%S                = 0
    ssnow%Tsoil            = 0
    ssnow%SL               = 0
    ssnow%TL               = 0
    ssnow%h0               = 0
    ssnow%rex              = 0
    ssnow%wflux            = 0
    ssnow%delwcol          = 0
    ssnow%zdelta           = 0
    ssnow%kth              = 0
    ssnow%Tsurface         = 0
    ssnow%lE               = 0
    ssnow%evap             = 0
    ssnow%ciso             = 0
    ssnow%cisoL            = 0
    ssnow%rlitt            = 0
    ssnow%thetai           = 0
    ssnow%snowliq          = 0
    ssnow%nsteps           = 0
    ssnow%nsnow            = 0
    ssnow%TsurfaceFR       = 0
    ssnow%Ta_daily         = 0
    ssnow%Qadv_daily       = 0
    ssnow%G0_daily         = 0
    ssnow%Qevap_daily      = 0
    ssnow%Qprec_daily      = 0
    ssnow%Qprec_snow_daily = 0
    ssnow%E_fusion_sn      = 0
    ssnow%E_sublimation_sn = 0
    ssnow%latent_heat_sn   = 0
    ssnow%evap_liq_sn      = 0
    ssnow%surface_melt     = 0
    ssnow%Qadv_rain_sn     = 0

  end subroutine zero_soil_snow_type


  ! ------------------------------------------------------------------


  subroutine zero_veg_parameter_type(veg)

    type(veg_parameter_type), intent(inout) :: veg

    veg%canst1      = 0
    veg%dleaf       = 0
    veg%ejmax       = 0
    veg%ejmax_shade = 0
    veg%ejmax_sun   = 0
    veg%iveg        = 0
    veg%ivegp       = 0
    veg%iLU         = 0
    veg%meth        = 0
    veg%frac4       = 0
    veg%hc          = 0
    veg%vlai        = 0
    veg%xalbnir     = 0
    veg%rp20        = 0
    veg%rpcoef      = 0
    veg%rs20        = 0
    veg%shelrb      = 0
    veg%vegcf       = 0
    veg%tminvj      = 0
    veg%toptvj      = 0
    veg%tmaxvj      = 0
    veg%vbeta       = 0
    veg%vcmax       = 0
    veg%vcmax_shade = 0
    veg%vcmax_sun   = 0
    veg%xfang       = 0
    veg%extkn       = 0
    veg%wai         = 0
    veg%deciduous   = .false.
    veg%froot       = 0
    veg%refl        = 0
    veg%taul        = 0
    veg%vlaimax     = 0
    veg%a1gs        = 0
    veg%d0gs        = 0
    veg%alpha       = 0
    veg%convex      = 0
    veg%cfrd        = 0
    veg%gswmin      = 0
    veg%conkc0      = 0
    veg%conko0      = 0
    veg%ekc         = 0
    veg%eko         = 0
    veg%g0          = 0
    veg%g1          = 0
    veg%vcmaxcc     = 0
    veg%ejmaxcc     = 0
    veg%gmmax       = 0
    veg%gm          = 0
    veg%c4kci       = 0
    veg%c4kcc       = 0
    veg%bjv         = 0
    veg%rootbeta = 0
    veg%gamma    = 0
    veg%F10      = 0
    veg%ZR       = 0
    veg%clitt    = 0
    veg%disturbance_interval  = 0
    veg%disturbance_intensity = 0

  end subroutine zero_veg_parameter_type


  ! ------------------------------------------------------------------


  subroutine zero_canopy_type(canopy)

    type(canopy_type), intent(inout) :: canopy

    canopy%fess           = 0
    canopy%fesp           = 0
    canopy%cansto         = 0
    canopy%cduv           = 0
    canopy%delwc          = 0
    canopy%dewmm          = 0
    canopy%dgdtg          = 0
    canopy%fe             = 0
    canopy%fh             = 0
    canopy%fpn            = 0
    canopy%frp            = 0
    canopy%frpw           = 0
    canopy%frpr           = 0
    canopy%frs            = 0
    canopy%fnee           = 0
    canopy%frday          = 0
    canopy%fnv            = 0
    canopy%fev            = 0
    canopy%fevc           = 0
    canopy%fhv            = 0
    canopy%fns            = 0
    canopy%fhs            = 0
    canopy%fhs_cor        = 0
    canopy%ga             = 0
    canopy%ghflux         = 0
    canopy%precis         = 0
    canopy%qscrn          = 0
    canopy%rnet           = 0
    canopy%rniso          = 0
    canopy%segg           = 0
    canopy%sghflux        = 0
    canopy%through        = 0
    canopy%spill          = 0
    canopy%tscrn          = 0
    canopy%wcint          = 0
    canopy%tv             = 0
    canopy%us             = 0
    canopy%uscrn          = 0
    canopy%rghlai         = 0
    canopy%vlaiw          = 0
    canopy%fwet           = 0
    canopy%A_sh           = 0
    canopy%A_sl           = 0
    canopy%A_slC          = 0
    canopy%A_shC          = 0
    canopy%A_slJ          = 0
    canopy%A_shJ          = 0
    canopy%GPP_sh         = 0
    canopy%GPP_sl         = 0
    canopy%fevc_sh        = 0
    canopy%fevc_sl        = 0
    canopy%eta_GPP_cs     = 0
    canopy%eta_fevc_cs    = 0
    canopy%eta_A_cs       = 0
    canopy%eta_A_cs_sh    = 0
    canopy%eta_A_cs_sl    = 0
    canopy%eta_fevc_cs_sh = 0
    canopy%eta_fevc_cs_sl = 0
    canopy%cs             = 0
    canopy%dAdcs          = 0
    canopy%cs_sl          = 0
    canopy%cs_sh          = 0
    canopy%tlf            = 0
    canopy%dlf            = 0
    canopy%evapfbl   = 0
    canopy%epot      = 0
    canopy%fnpp      = 0
    canopy%fevw_pot  = 0
    canopy%gswx_T    = 0
    canopy%cdtq      = 0
    canopy%wetfac_cs = 0
    canopy%fevw      = 0
    canopy%fhvw      = 0
    canopy%fes       = 0
    canopy%fes_cor   = 0
    canopy%gswx      = 0
    canopy%oldcansto = 0
    canopy%zetar     = 0
    canopy%zetash    = 0
    canopy%fwsoil    = 0
    canopy%ofes      = 0
    canopy%gw     = 0
    canopy%ancj   = 0
    canopy%tlfy   = 0
    canopy%ecy    = 0
    canopy%ecx    = 0
    canopy%fwsoil = 0
    canopy%kthLitt = 0
    canopy%DvLitt  = 0
    canopy%An        = 0
    canopy%Rd        = 0
    canopy%isc3      = .false.
    canopy%vcmax     = 0
    canopy%gammastar = 0
    canopy%gsc       = 0
    canopy%gbc       = 0
    canopy%gac       = 0
    canopy%ci        = 0

  end subroutine zero_canopy_type


  ! ------------------------------------------------------------------


  subroutine zero_radiation_type(rad)

    type(radiation_type), intent(inout) :: rad

    rad%albedo    = 0
    rad%extkb     = 0
    rad%extkd2    = 0
    rad%extkd     = 0
    rad%flws      = 0
    rad%fvlai     = 0
    rad%latitude  = 0
    rad%lwabv     = 0
    rad%qcan      = 0
    rad%qssabs    = 0
    rad%rhocdf    = 0
    rad%rniso     = 0
    rad%scalex    = 0
    rad%transd    = 0
    rad%trad      = 0
    rad%reffdf    = 0
    rad%reffbm    = 0
    rad%extkbm    = 0
    rad%extkdm    = 0
    rad%cexpkbm   = 0
    rad%cexpkdm   = 0
    rad%fbeam     = 0
    rad%rhocbm    = 0
    rad%transb    = 0
    rad%albedo_T  = 0
    rad%gradis    = 0
    rad%longitude = 0
    rad%workp1    = 0
    rad%workp2    = 0
    rad%workp3    = 0

  end subroutine zero_radiation_type


  ! ------------------------------------------------------------------


  subroutine zero_roughness_type(rought)

    type(roughness_type), intent(inout) :: rought

    rought%coexp      = 0
    rought%disp       = 0
    rought%hruff      = 0
    rought%hruff_grmx = 0
    rought%rt0us      = 0
    rought%rt1usa     = 0
    rought%rt1usb     = 0
    rought%rt1        = 0
    rought%term2      = 0
    rought%term3      = 0
    rought%term5      = 0
    rought%term6      = 0
    rought%term6a     = 0
    rought%usuh       = 0
    rought%za_uv      = 0
    rought%za_tq      = 0
    rought%z0m        = 0
    rought%zref_uv    = 0
    rought%zref_tq    = 0
    rought%zruffs     = 0
    rought%z0soilsn   = 0
    rought%z0soil     = 0

  end subroutine zero_roughness_type


  ! ------------------------------------------------------------------


  subroutine zero_air_type(air)

    type(air_type), intent(inout) :: air

    air%rho    = 0
    air%volm   = 0
    air%rlam   = 0
    air%qsat   = 0
    air%epsi   = 0
    air%visc   = 0
    air%psyc   = 0
    air%dsatdk = 0
    air%cmolar = 0

  end subroutine zero_air_type


  ! ------------------------------------------------------------------


  subroutine zero_met_type(met)

    type(met_type), intent(inout) :: met

    met%ca        = 0
    met%year      = 0
    met%moy       = 0
    met%doy       = 0
    met%hod       = 0
    met%fsd       = 0
    met%ofsd      = 0
    met%fld       = 0
    met%precip    = 0
    met%precip_sn = 0
    met%tk        = 0
    met%tvair     = 0
    met%tvrad     = 0
    met%pmb       = 0
    met%ua        = 0
    met%qv        = 0
    met%qvair     = 0
    met%da        = 0
    met%dva       = 0
    met%coszen    = 0
    met%Ndep      = 0
    met%Pdep      = 0
    met%rhum      = 0
    met%fdiff     = 0
    met%u10       = 0

  end subroutine zero_met_type


  ! ------------------------------------------------------------------


  subroutine zero_climate_type(climate)

    implicit none

    type(climate_type), intent(inout) :: climate

    climate%chilldays     = 0
    climate%iveg          = 0
    climate%biome         = 0
    climate%GMD           = 0
    climate%modis_igbp    = 0
    climate%DSLR          = 0
    climate%NDAY_Nesterov = 0
    climate%dtemp                      = 0
    climate%dmoist                     = 0
    climate%dmoist_min                 = 0
    climate%dmoist_min20               = 0
    climate%dmoist_max                 = 0
    climate%dmoist_max20               = 0
    climate%mtemp                      = 0
    climate%qtemp                      = 0
    climate%mmoist                     = 0
    climate%mtemp_min                  = 0
    climate%mtemp_max                  = 0
    climate%qtemp_max                  = 0
    climate%qtemp_max_last_year        = 0
    climate%mtemp_min20                = 0
    climate%mtemp_max20                = 0
    climate%atemp_mean                 = 0
    climate%AGDD5                      = 0
    climate%GDD5                       = 0
    climate%AGDD0                      = 0
    climate%GDD0                       = 0
    climate%alpha_PT                   = 0
    climate%evap_PT                    = 0
    climate%aevap                      = 0
    climate%alpha_PT20                 = 0
    climate%GDD0_rec                   = 0
    climate%frec                       = 0
    climate%dtemp_min                  = 0
    climate%fdorm                      = 0
    climate%fapar_ann_max              = 0
    climate%fapar_ann_max_last_year    = 0
    climate%AvgAnnMaxFAPAR             = 0
    climate%dtemp_max                  = 0
    climate%drhum                      = 0
    climate%du10_max                   = 0
    climate%dprecip                    = 0
    climate%aprecip                    = 0
    climate%aprecip_av20               = 0
    climate%last_precip                = 0
    climate%KBDI                       = 0
    climate%FFDI                       = 0
    climate%D_MacArthur                = 0
    climate%Nesterov_Current           = 0
    climate%Nesterov_ann_max           = 0
    climate%Nesterov_ann_max_last_year = 0
    climate%Nesterov_ann_running_max   = 0
    climate%mtemp_min_20    = 0
    climate%mtemp_max_20    = 0
    climate%dmoist_min_20   = 0
    climate%dmoist_max_20   = 0
    climate%dtemp_31        = 0
    climate%dmoist_31       = 0
    climate%alpha_PT_20     = 0
    climate%dtemp_91        = 0
    climate%APAR_leaf_sun   = 0
    climate%APAR_leaf_shade = 0
    climate%Dleaf_sun       = 0
    climate%Dleaf_shade     = 0
    climate%Tleaf_sun       = 0
    climate%Tleaf_shade     = 0
    climate%cs_sun          = 0
    climate%cs_shade        = 0
    climate%scalex_sun      = 0
    climate%scalex_shade    = 0
    climate%fwsoil          = 0
    climate%aprecip_20      = 0
    climate%Rd_sun          = 0
    climate%Rd_shade        = 0

  end subroutine zero_climate_type


  ! ------------------------------------------------------------------


  subroutine zero_sum_flux_type(sum_flux)

    type(sum_flux_type), intent(inout) :: sum_flux

    sum_flux%sumpn  = 0
    sum_flux%sumrp  = 0
    sum_flux%sumrpw = 0
    sum_flux%sumrpr = 0
    sum_flux%sumrs  = 0
    sum_flux%sumrd  = 0
    sum_flux%dsumpn = 0
    sum_flux%dsumrp = 0
    sum_flux%dsumrs = 0
    sum_flux%dsumrd = 0
    sum_flux%sumxrp = 0
    sum_flux%sumxrs = 0

  end subroutine zero_sum_flux_type


  ! ------------------------------------------------------------------


  subroutine zero_bgc_pool_type(bgc)

    type(bgc_pool_type), intent(inout) :: bgc

    bgc%cplant = 0
    bgc%csoil  = 0

  end subroutine zero_bgc_pool_type


  ! ------------------------------------------------------------------


  subroutine print_balances_type(bal)

    type(balances_type), intent(in) :: bal

    write(*,*) 'bal%drybal ', bal%drybal
    write(*,*) 'bal%ebal ', bal%ebal
    write(*,*) 'bal%ebal_tot ', bal%ebal_tot
    write(*,*) 'bal%ebaltr ', bal%ebaltr
    write(*,*) 'bal%ebal_tottr ', bal%ebal_tottr
    write(*,*) 'bal%ebal_cncheck ', bal%ebal_cncheck
    write(*,*) 'bal%ebal_tot_cncheck ', bal%ebal_tot_cncheck
    write(*,*) 'bal%evap_tot ', bal%evap_tot
    write(*,*) 'bal%osnowd0 ', bal%osnowd0
    write(*,*) 'bal%precip_tot ', bal%precip_tot
    write(*,*) 'bal%rnoff_tot ', bal%rnoff_tot
    write(*,*) 'bal%wbal ', bal%wbal
    write(*,*) 'bal%wbal_tot ', bal%wbal_tot
    write(*,*) 'bal%wbtot0 ', bal%wbtot0
    write(*,*) 'bal%wetbal ', bal%wetbal
    write(*,*) 'bal%cansto0 ', bal%cansto0
    write(*,*) 'bal%evapc_tot ', bal%evapc_tot
    write(*,*) 'bal%evaps_tot ', bal%evaps_tot
    write(*,*) 'bal%rnof1_tot ', bal%rnof1_tot
    write(*,*) 'bal%rnof2_tot ', bal%rnof2_tot
    write(*,*) 'bal%snowdc_tot ', bal%snowdc_tot
    write(*,*) 'bal%wbal_tot1 ', bal%wbal_tot1
    write(*,*) 'bal%owbtot ', bal%owbtot
    write(*,*) 'bal%delwc_tot ', bal%delwc_tot
    write(*,*) 'bal%qasrf_tot ', bal%qasrf_tot
    write(*,*) 'bal%qfsrf_tot ', bal%qfsrf_tot
    write(*,*) 'bal%qssrf_tot ', bal%qssrf_tot
    write(*,*) 'bal%Radbal ', bal%Radbal
    write(*,*) 'bal%EbalSoil ', bal%EbalSoil
    write(*,*) 'bal%Ebalveg ', bal%Ebalveg
    write(*,*) 'bal%Radbalsum ', bal%Radbalsum

  end subroutine print_balances_type


  ! ------------------------------------------------------------------


  subroutine print_soil_parameter_type(soil)

    type(soil_parameter_type), intent(in) :: soil

    write(*,*) 'soil%bch ', soil%bch
    write(*,*) 'soil%c3 ', soil%c3
    write(*,*) 'soil%clay ', soil%clay
    write(*,*) 'soil%css ', soil%css
    write(*,*) 'soil%hsbh ', soil%hsbh
    write(*,*) 'soil%hyds ', soil%hyds
    write(*,*) 'soil%i2bp3 ', soil%i2bp3
    write(*,*) 'soil%ibp2 ', soil%ibp2
    write(*,*) 'soil%isoilm ', soil%isoilm
    write(*,*) 'soil%rhosoil ', soil%rhosoil
    write(*,*) 'soil%sand ', soil%sand
    write(*,*) 'soil%sfc ', soil%sfc
    write(*,*) 'soil%silt ', soil%silt
    write(*,*) 'soil%ssat ', soil%ssat
    write(*,*) 'soil%sucs ', soil%sucs
    write(*,*) 'soil%swilt ', soil%swilt
    write(*,*) 'soil%zse ', soil%zse
    write(*,*) 'soil%zshh ', soil%zshh
    write(*,*) 'soil%cnsd ', soil%cnsd
    write(*,*) 'soil%albsoil ', soil%albsoil
    write(*,*) 'soil%pwb_min ', soil%pwb_min
    write(*,*) 'soil%albsoilf ', soil%albsoilf
    write(*,*) 'soil%soilcol ', soil%soilcol
    write(*,*) 'soil%nhorizons ', soil%nhorizons
    write(*,*) 'soil%ishorizon ', soil%ishorizon
    write(*,*) 'soil%clitt ', soil%clitt
    write(*,*) 'soil%zeta ', soil%zeta
    write(*,*) 'soil%fsatmax ', soil%fsatmax
    write(*,*) 'soil%swilt_vec ', soil%swilt_vec
    write(*,*) 'soil%ssat_vec ', soil%ssat_vec
    write(*,*) 'soil%sfc_vec ', soil%sfc_vec
    write(*,*) 'soil%swilt_vec ', soil%swilt_vec
    write(*,*) 'soil%ssat_vec ', soil%ssat_vec
    write(*,*) 'soil%sfc_vec ', soil%sfc_vec

  end subroutine print_soil_parameter_type


  ! ------------------------------------------------------------------


  subroutine print_soil_snow_type(ssnow)

    type(soil_snow_type), intent(in) :: ssnow

    write(*,*) 'ssnow%iantrct ', ssnow%iantrct
    write(*,*) 'ssnow%pudsto ', ssnow%pudsto
    write(*,*) 'ssnow%pudsmx ', ssnow%pudsmx
    write(*,*) 'ssnow%dtmlt ', ssnow%dtmlt
    write(*,*) 'ssnow%albsoilsn ', ssnow%albsoilsn
    write(*,*) 'ssnow%cls ', ssnow%cls
    write(*,*) 'ssnow%dfn_dtg ', ssnow%dfn_dtg
    write(*,*) 'ssnow%dfh_dtg ', ssnow%dfh_dtg
    write(*,*) 'ssnow%dfe_ddq ', ssnow%dfe_ddq
    write(*,*) 'ssnow%ddq_dtg ', ssnow%ddq_dtg
    write(*,*) 'ssnow%evapsn ', ssnow%evapsn
    write(*,*) 'ssnow%fwtop ', ssnow%fwtop
    write(*,*) 'ssnow%fwtop1 ', ssnow%fwtop1
    write(*,*) 'ssnow%fwtop2 ', ssnow%fwtop2
    write(*,*) 'ssnow%fwtop3 ', ssnow%fwtop3
    write(*,*) 'ssnow%gammzz ', ssnow%gammzz
    write(*,*) 'ssnow%isflag ', ssnow%isflag
    write(*,*) 'ssnow%osnowd ', ssnow%osnowd
    write(*,*) 'ssnow%potev ', ssnow%potev
    write(*,*) 'ssnow%runoff ', ssnow%runoff
    write(*,*) 'ssnow%rnof1 ', ssnow%rnof1
    write(*,*) 'ssnow%rnof2 ', ssnow%rnof2
    write(*,*) 'ssnow%rtsoil ', ssnow%rtsoil
    write(*,*) 'ssnow%sconds ', ssnow%sconds
    write(*,*) 'ssnow%sdepth ', ssnow%sdepth
    write(*,*) 'ssnow%smass ', ssnow%smass
    write(*,*) 'ssnow%snage ', ssnow%snage
    write(*,*) 'ssnow%snowd ', ssnow%snowd
    write(*,*) 'ssnow%smelt ', ssnow%smelt
    write(*,*) 'ssnow%ssdn ', ssnow%ssdn
    write(*,*) 'ssnow%ssdnn ', ssnow%ssdnn
    write(*,*) 'ssnow%tgg ', ssnow%tgg
    write(*,*) 'ssnow%tggsn ', ssnow%tggsn
    write(*,*) 'ssnow%tss ', ssnow%tss
    write(*,*) 'ssnow%tss_p ', ssnow%tss_p
    write(*,*) 'ssnow%deltss ', ssnow%deltss
    write(*,*) 'ssnow%owb1 ', ssnow%owb1
    write(*,*) 'ssnow%wb ', ssnow%wb
    write(*,*) 'ssnow%wbice ', ssnow%wbice
    write(*,*) 'ssnow%wblf ', ssnow%wblf
    write(*,*) 'ssnow%wbtot ', ssnow%wbtot
    write(*,*) 'ssnow%wbtot1 ', ssnow%wbtot1
    write(*,*) 'ssnow%wbtot2 ', ssnow%wbtot2
    write(*,*) 'ssnow%wb_lake ', ssnow%wb_lake
    write(*,*) 'ssnow%sinfil ', ssnow%sinfil
    write(*,*) 'ssnow%evapfbl ', ssnow%evapfbl
    write(*,*) 'ssnow%qstss ', ssnow%qstss
    write(*,*) 'ssnow%wetfac ', ssnow%wetfac
    write(*,*) 'ssnow%owetfac ', ssnow%owetfac
    write(*,*) 'ssnow%t_snwlr ', ssnow%t_snwlr
    write(*,*) 'ssnow%wbfice ', ssnow%wbfice
    write(*,*) 'ssnow%tggav ', ssnow%tggav
    write(*,*) 'ssnow%otgg ', ssnow%otgg
    write(*,*) 'ssnow%otss ', ssnow%otss
    write(*,*) 'ssnow%otss_0 ', ssnow%otss_0
    write(*,*) 'ssnow%tprecip ', ssnow%tprecip
    write(*,*) 'ssnow%tevap ', ssnow%tevap
    write(*,*) 'ssnow%trnoff ', ssnow%trnoff
    write(*,*) 'ssnow%totenbal ', ssnow%totenbal
    write(*,*) 'ssnow%totenbal2 ', ssnow%totenbal2
    write(*,*) 'ssnow%fland ', ssnow%fland
    write(*,*) 'ssnow%ifland ', ssnow%ifland
    write(*,*) 'ssnow%tilefrac ', ssnow%tilefrac
    write(*,*) 'ssnow%qasrf ', ssnow%qasrf
    write(*,*) 'ssnow%qfsrf ', ssnow%qfsrf
    write(*,*) 'ssnow%qssrf ', ssnow%qssrf
    write(*,*) 'ssnow%S ', ssnow%S
    write(*,*) 'ssnow%Tsoil ', ssnow%Tsoil
    write(*,*) 'ssnow%SL ', ssnow%SL
    write(*,*) 'ssnow%TL ', ssnow%TL
    write(*,*) 'ssnow%h0 ', ssnow%h0
    write(*,*) 'ssnow%rex ', ssnow%rex
    write(*,*) 'ssnow%wflux ', ssnow%wflux
    write(*,*) 'ssnow%delwcol ', ssnow%delwcol
    write(*,*) 'ssnow%zdelta ', ssnow%zdelta
    write(*,*) 'ssnow%kth ', ssnow%kth
    write(*,*) 'ssnow%Tsurface ', ssnow%Tsurface
    write(*,*) 'ssnow%lE ', ssnow%lE
    write(*,*) 'ssnow%evap ', ssnow%evap
    write(*,*) 'ssnow%ciso ', ssnow%ciso
    write(*,*) 'ssnow%cisoL ', ssnow%cisoL
    write(*,*) 'ssnow%rlitt ', ssnow%rlitt
    write(*,*) 'ssnow%thetai ', ssnow%thetai
    write(*,*) 'ssnow%snowliq ', ssnow%snowliq
    write(*,*) 'ssnow%nsteps ', ssnow%nsteps
    write(*,*) 'ssnow%nsnow ', ssnow%nsnow
    write(*,*) 'ssnow%TsurfaceFR ', ssnow%TsurfaceFR
    write(*,*) 'ssnow%Ta_daily ', ssnow%Ta_daily
    write(*,*) 'ssnow%Qadv_daily ', ssnow%Qadv_daily
    write(*,*) 'ssnow%G0_daily ', ssnow%G0_daily
    write(*,*) 'ssnow%Qevap_daily ', ssnow%Qevap_daily
    write(*,*) 'ssnow%Qprec_daily ', ssnow%Qprec_daily
    write(*,*) 'ssnow%Qprec_snow_daily ', ssnow%Qprec_snow_daily
    write(*,*) 'ssnow%E_fusion_sn ', ssnow%E_fusion_sn
    write(*,*) 'ssnow%E_sublimation_sn ', ssnow%E_sublimation_sn
    write(*,*) 'ssnow%latent_heat_sn ', ssnow%latent_heat_sn
    write(*,*) 'ssnow%evap_liq_sn ', ssnow%evap_liq_sn
    write(*,*) 'ssnow%surface_melt ', ssnow%surface_melt
    write(*,*) 'ssnow%Qadv_rain_sn ', ssnow%Qadv_rain_sn

  end subroutine print_soil_snow_type


  ! ------------------------------------------------------------------


  subroutine print_veg_parameter_type(veg)

    type(veg_parameter_type), intent(in) :: veg

    write(*,*) 'veg%canst1 ', veg%canst1
    write(*,*) 'veg%dleaf ', veg%dleaf
    write(*,*) 'veg%ejmax ', veg%ejmax
    write(*,*) 'veg%ejmax_shade ', veg%ejmax_shade
    write(*,*) 'veg%ejmax_sun ', veg%ejmax_sun
    write(*,*) 'veg%iveg ', veg%iveg
    write(*,*) 'veg%ivegp ', veg%ivegp
    write(*,*) 'veg%iLU ', veg%iLU
    write(*,*) 'veg%meth ', veg%meth
    write(*,*) 'veg%frac4 ', veg%frac4
    write(*,*) 'veg%hc ', veg%hc
    write(*,*) 'veg%vlai ', veg%vlai
    write(*,*) 'veg%xalbnir ', veg%xalbnir
    write(*,*) 'veg%rp20 ', veg%rp20
    write(*,*) 'veg%rpcoef ', veg%rpcoef
    write(*,*) 'veg%rs20 ', veg%rs20
    write(*,*) 'veg%shelrb ', veg%shelrb
    write(*,*) 'veg%vegcf ', veg%vegcf
    write(*,*) 'veg%tminvj ', veg%tminvj
    write(*,*) 'veg%toptvj ', veg%toptvj
    write(*,*) 'veg%tmaxvj ', veg%tmaxvj
    write(*,*) 'veg%vbeta ', veg%vbeta
    write(*,*) 'veg%vcmax ', veg%vcmax
    write(*,*) 'veg%vcmax_shade ', veg%vcmax_shade
    write(*,*) 'veg%vcmax_sun ', veg%vcmax_sun
    write(*,*) 'veg%xfang ', veg%xfang
    write(*,*) 'veg%extkn ', veg%extkn
    write(*,*) 'veg%wai ', veg%wai
    write(*,*) 'veg%deciduous ', veg%deciduous
    write(*,*) 'veg%froot ', veg%froot
    write(*,*) 'veg%refl ', veg%refl
    write(*,*) 'veg%taul ', veg%taul
    write(*,*) 'veg%vlaimax ', veg%vlaimax
    write(*,*) 'veg%a1gs ', veg%a1gs
    write(*,*) 'veg%d0gs ', veg%d0gs
    write(*,*) 'veg%alpha ', veg%alpha
    write(*,*) 'veg%convex ', veg%convex
    write(*,*) 'veg%cfrd ', veg%cfrd
    write(*,*) 'veg%gswmin ', veg%gswmin
    write(*,*) 'veg%conkc0 ', veg%conkc0
    write(*,*) 'veg%conko0 ', veg%conko0
    write(*,*) 'veg%ekc ', veg%ekc
    write(*,*) 'veg%eko ', veg%eko
    write(*,*) 'veg%g0 ', veg%g0
    write(*,*) 'veg%g1 ', veg%g1
    write(*,*) 'veg%vcmaxcc ', veg%vcmaxcc
    write(*,*) 'veg%ejmaxcc ', veg%ejmaxcc
    write(*,*) 'veg%gmmax ', veg%gmmax
    write(*,*) 'veg%gm ', veg%gm
    write(*,*) 'veg%c4kci ', veg%c4kci
    write(*,*) 'veg%c4kcc ', veg%c4kcc
    write(*,*) 'veg%bjv ', veg%bjv
    write(*,*) 'veg%rootbeta ', veg%rootbeta
    write(*,*) 'veg%gamma ', veg%gamma
    write(*,*) 'veg%F10 ', veg%F10
    write(*,*) 'veg%ZR ', veg%ZR
    write(*,*) 'veg%clitt ', veg%clitt
    write(*,*) 'veg%disturbance_interval ', veg%disturbance_interval
    write(*,*) 'veg%disturbance_intensity ', veg%disturbance_intensity

  end subroutine print_veg_parameter_type


  ! ------------------------------------------------------------------


  subroutine print_canopy_type(canopy)

    type(canopy_type), intent(in) :: canopy

    write(*,*) 'canopy%fess ', canopy%fess
    write(*,*) 'canopy%fesp ', canopy%fesp
    write(*,*) 'canopy%cansto ', canopy%cansto
    write(*,*) 'canopy%cduv ', canopy%cduv
    write(*,*) 'canopy%delwc ', canopy%delwc
    write(*,*) 'canopy%dewmm ', canopy%dewmm
    write(*,*) 'canopy%dgdtg ', canopy%dgdtg
    write(*,*) 'canopy%fe ', canopy%fe
    write(*,*) 'canopy%fh ', canopy%fh
    write(*,*) 'canopy%fpn ', canopy%fpn
    write(*,*) 'canopy%frp ', canopy%frp
    write(*,*) 'canopy%frpw ', canopy%frpw
    write(*,*) 'canopy%frpr ', canopy%frpr
    write(*,*) 'canopy%frs ', canopy%frs
    write(*,*) 'canopy%fnee ', canopy%fnee
    write(*,*) 'canopy%frday ', canopy%frday
    write(*,*) 'canopy%fnv ', canopy%fnv
    write(*,*) 'canopy%fev ', canopy%fev
    write(*,*) 'canopy%fevc ', canopy%fevc
    write(*,*) 'canopy%fhv ', canopy%fhv
    write(*,*) 'canopy%fns ', canopy%fns
    write(*,*) 'canopy%fhs ', canopy%fhs
    write(*,*) 'canopy%fhs_cor ', canopy%fhs_cor
    write(*,*) 'canopy%ga ', canopy%ga
    write(*,*) 'canopy%ghflux ', canopy%ghflux
    write(*,*) 'canopy%precis ', canopy%precis
    write(*,*) 'canopy%qscrn ', canopy%qscrn
    write(*,*) 'canopy%rnet ', canopy%rnet
    write(*,*) 'canopy%rniso ', canopy%rniso
    write(*,*) 'canopy%segg ', canopy%segg
    write(*,*) 'canopy%sghflux ', canopy%sghflux
    write(*,*) 'canopy%through ', canopy%through
    write(*,*) 'canopy%spill ', canopy%spill
    write(*,*) 'canopy%tscrn ', canopy%tscrn
    write(*,*) 'canopy%wcint ', canopy%wcint
    write(*,*) 'canopy%tv ', canopy%tv
    write(*,*) 'canopy%us ', canopy%us
    write(*,*) 'canopy%uscrn ', canopy%uscrn
    write(*,*) 'canopy%rghlai ', canopy%rghlai
    write(*,*) 'canopy%vlaiw ', canopy%vlaiw
    write(*,*) 'canopy%fwet ', canopy%fwet
    write(*,*) 'canopy%A_sh ', canopy%A_sh
    write(*,*) 'canopy%A_sl ', canopy%A_sl
    write(*,*) 'canopy%A_slC ', canopy%A_slC
    write(*,*) 'canopy%A_shC ', canopy%A_shC
    write(*,*) 'canopy%A_slJ ', canopy%A_slJ
    write(*,*) 'canopy%A_shJ ', canopy%A_shJ
    write(*,*) 'canopy%GPP_sh ', canopy%GPP_sh
    write(*,*) 'canopy%GPP_sl ', canopy%GPP_sl
    write(*,*) 'canopy%fevc_sh ', canopy%fevc_sh
    write(*,*) 'canopy%fevc_sl ', canopy%fevc_sl
    write(*,*) 'canopy%eta_GPP_cs ', canopy%eta_GPP_cs
    write(*,*) 'canopy%eta_fevc_cs ', canopy%eta_fevc_cs
    write(*,*) 'canopy%eta_A_cs ', canopy%eta_A_cs
    write(*,*) 'canopy%eta_A_cs_sh ', canopy%eta_A_cs_sh
    write(*,*) 'canopy%eta_A_cs_sl ', canopy%eta_A_cs_sl
    write(*,*) 'canopy%eta_fevc_cs_sh ', canopy%eta_fevc_cs_sh
    write(*,*) 'canopy%eta_fevc_cs_sl ', canopy%eta_fevc_cs_sl
    write(*,*) 'canopy%cs ', canopy%cs
    write(*,*) 'canopy%dAdcs ', canopy%dAdcs
    write(*,*) 'canopy%cs_sl ', canopy%cs_sl
    write(*,*) 'canopy%cs_sh ', canopy%cs_sh
    write(*,*) 'canopy%tlf ', canopy%tlf
    write(*,*) 'canopy%dlf ', canopy%dlf
    write(*,*) 'canopy%evapfbl ', canopy%evapfbl
    write(*,*) 'canopy%epot ', canopy%epot
    write(*,*) 'canopy%fnpp ', canopy%fnpp
    write(*,*) 'canopy%fevw_pot ', canopy%fevw_pot
    write(*,*) 'canopy%gswx_T ', canopy%gswx_T
    write(*,*) 'canopy%cdtq ', canopy%cdtq
    write(*,*) 'canopy%wetfac_cs ', canopy%wetfac_cs
    write(*,*) 'canopy%fevw ', canopy%fevw
    write(*,*) 'canopy%fhvw ', canopy%fhvw
    write(*,*) 'canopy%fes ', canopy%fes
    write(*,*) 'canopy%fes_cor ', canopy%fes_cor
    write(*,*) 'canopy%gswx ', canopy%gswx
    write(*,*) 'canopy%oldcansto ', canopy%oldcansto
    write(*,*) 'canopy%zetar ', canopy%zetar
    write(*,*) 'canopy%zetash ', canopy%zetash
    write(*,*) 'canopy%fwsoil ', canopy%fwsoil
    write(*,*) 'canopy%ofes ', canopy%ofes
    write(*,*) 'canopy%gw ', canopy%gw
    write(*,*) 'canopy%ancj ', canopy%ancj
    write(*,*) 'canopy%tlfy ', canopy%tlfy
    write(*,*) 'canopy%ecy ', canopy%ecy
    write(*,*) 'canopy%ecx ', canopy%ecx
    write(*,*) 'canopy%fwsoil ', canopy%fwsoil
    write(*,*) 'canopy%kthLitt ', canopy%kthLitt
    write(*,*) 'canopy%DvLitt ', canopy%DvLitt
    write(*,*) 'canopy%An ', canopy%An
    write(*,*) 'canopy%Rd ', canopy%Rd
    write(*,*) 'canopy%isc3 ', canopy%isc3
    write(*,*) 'canopy%vcmax ', canopy%vcmax
    write(*,*) 'canopy%gammastar ', canopy%gammastar
    write(*,*) 'canopy%gsc ', canopy%gsc
    write(*,*) 'canopy%gbc ', canopy%gbc
    write(*,*) 'canopy%gac ', canopy%gac
    write(*,*) 'canopy%ci ', canopy%ci

  end subroutine print_canopy_type


  ! ------------------------------------------------------------------


  subroutine print_radiation_type(rad)

    type(radiation_type), intent(in) :: rad

    write(*,*) 'rad%albedo ', rad%albedo
    write(*,*) 'rad%extkb ', rad%extkb
    write(*,*) 'rad%extkd2 ', rad%extkd2
    write(*,*) 'rad%extkd ', rad%extkd
    write(*,*) 'rad%flws ', rad%flws
    write(*,*) 'rad%fvlai ', rad%fvlai
    write(*,*) 'rad%latitude ', rad%latitude
    write(*,*) 'rad%lwabv ', rad%lwabv
    write(*,*) 'rad%qcan ', rad%qcan
    write(*,*) 'rad%qssabs ', rad%qssabs
    write(*,*) 'rad%rhocdf ', rad%rhocdf
    write(*,*) 'rad%rniso ', rad%rniso
    write(*,*) 'rad%scalex ', rad%scalex
    write(*,*) 'rad%transd ', rad%transd
    write(*,*) 'rad%trad ', rad%trad
    write(*,*) 'rad%reffdf ', rad%reffdf
    write(*,*) 'rad%reffbm ', rad%reffbm
    write(*,*) 'rad%extkbm ', rad%extkbm
    write(*,*) 'rad%extkdm ', rad%extkdm
    write(*,*) 'rad%cexpkbm ', rad%cexpkbm
    write(*,*) 'rad%cexpkdm ', rad%cexpkdm
    write(*,*) 'rad%fbeam ', rad%fbeam
    write(*,*) 'rad%rhocbm ', rad%rhocbm
    write(*,*) 'rad%transb ', rad%transb
    write(*,*) 'rad%albedo_T ', rad%albedo_T
    write(*,*) 'rad%gradis ', rad%gradis
    write(*,*) 'rad%longitude ', rad%longitude
    write(*,*) 'rad%workp1 ', rad%workp1
    write(*,*) 'rad%workp2 ', rad%workp2
    write(*,*) 'rad%workp3 ', rad%workp3

  end subroutine print_radiation_type


  ! ------------------------------------------------------------------


  subroutine print_roughness_type(rought)

    type(roughness_type), intent(in) :: rought

    write(*,*) 'rough%coexp ', rought%coexp
    write(*,*) 'rough%disp ', rought%disp
    write(*,*) 'rough%hruff ', rought%hruff
    write(*,*) 'rough%hruff_grmx ', rought%hruff_grmx
    write(*,*) 'rough%rt0us ', rought%rt0us
    write(*,*) 'rough%rt1usa ', rought%rt1usa
    write(*,*) 'rough%rt1usb ', rought%rt1usb
    write(*,*) 'rough%rt1 ', rought%rt1
    write(*,*) 'rough%term2 ', rought%term2
    write(*,*) 'rough%term3 ', rought%term3
    write(*,*) 'rough%term5 ', rought%term5
    write(*,*) 'rough%term6 ', rought%term6
    write(*,*) 'rough%term6a ', rought%term6a
    write(*,*) 'rough%usuh ', rought%usuh
    write(*,*) 'rough%za_uv ', rought%za_uv
    write(*,*) 'rough%za_tq ', rought%za_tq
    write(*,*) 'rough%z0m ', rought%z0m
    write(*,*) 'rough%zref_uv ', rought%zref_uv
    write(*,*) 'rough%zref_tq ', rought%zref_tq
    write(*,*) 'rough%zruffs ', rought%zruffs
    write(*,*) 'rough%z0soilsn ', rought%z0soilsn
    write(*,*) 'rough%z0soil ', rought%z0soil

  end subroutine print_roughness_type


  ! ------------------------------------------------------------------


  subroutine print_air_type(air)

    type(air_type), intent(in) :: air

    write(*,*) 'air%rho ', air%rho
    write(*,*) 'air%volm ', air%volm
    write(*,*) 'air%rlam ', air%rlam
    write(*,*) 'air%qsat ', air%qsat
    write(*,*) 'air%epsi ', air%epsi
    write(*,*) 'air%visc ', air%visc
    write(*,*) 'air%psyc ', air%psyc
    write(*,*) 'air%dsatdk ', air%dsatdk
    write(*,*) 'air%cmolar ', air%cmolar

  end subroutine print_air_type


  ! ------------------------------------------------------------------


  subroutine print_met_type(met)

    type(met_type), intent(in) :: met

    write(*,*) 'met%ca ', met%ca
    write(*,*) 'met%year ', met%year
    write(*,*) 'met%moy ', met%moy
    write(*,*) 'met%doy ', met%doy
    write(*,*) 'met%hod ', met%hod
    write(*,*) 'met%fsd ', met%fsd
    write(*,*) 'met%ofsd ', met%ofsd
    write(*,*) 'met%fld ', met%fld
    write(*,*) 'met%precip ', met%precip
    write(*,*) 'met%precip_sn ', met%precip_sn
    write(*,*) 'met%tk ', met%tk
    write(*,*) 'met%tvair ', met%tvair
    write(*,*) 'met%tvrad ', met%tvrad
    write(*,*) 'met%pmb ', met%pmb
    write(*,*) 'met%ua ', met%ua
    write(*,*) 'met%qv ', met%qv
    write(*,*) 'met%qvair ', met%qvair
    write(*,*) 'met%da ', met%da
    write(*,*) 'met%dva ', met%dva
    write(*,*) 'met%coszen ', met%coszen
    write(*,*) 'met%Ndep ', met%Ndep
    write(*,*) 'met%Pdep ', met%Pdep
    write(*,*) 'met%rhum ', met%rhum
    write(*,*) 'met%fdiff ', met%fdiff
    write(*,*) 'met%u10 ', met%u10

  end subroutine print_met_type


  ! ------------------------------------------------------------------


  subroutine print_climate_type(met)

    implicit none

    type(climate_type), intent(in) :: met

    write(*,*) 'climate%chilldays ', met%chilldays
    write(*,*) 'climate%iveg ', met%iveg
    write(*,*) 'climate%biome ', met%biome
    write(*,*) 'climate%GMD ', met%GMD
    write(*,*) 'climate%modis_igbp ', met%modis_igbp
    write(*,*) 'climate%DSLR ', met%DSLR
    write(*,*) 'climate%NDAY_Nesterov ', met%NDAY_Nesterov
    write(*,*) 'climate%dtemp ', met%dtemp
    write(*,*) 'climate%dmoist ', met%dmoist
    write(*,*) 'climate%dmoist_min ', met%dmoist_min
    write(*,*) 'climate%dmoist_min20 ', met%dmoist_min20
    write(*,*) 'climate%dmoist_max ', met%dmoist_max
    write(*,*) 'climate%dmoist_max20 ', met%dmoist_max20
    write(*,*) 'climate%mtemp ', met%mtemp
    write(*,*) 'climate%qtemp ', met%qtemp
    write(*,*) 'climate%mmoist ', met%mmoist
    write(*,*) 'climate%mtemp_min ', met%mtemp_min
    write(*,*) 'climate%mtemp_max ', met%mtemp_max
    write(*,*) 'climate%qtemp_max ', met%qtemp_max
    write(*,*) 'climate%qtemp_max_last_year ', met%qtemp_max_last_year
    write(*,*) 'climate%mtemp_min20 ', met%mtemp_min20
    write(*,*) 'climate%mtemp_max20 ', met%mtemp_max20
    write(*,*) 'climate%atemp_mean ', met%atemp_mean
    write(*,*) 'climate%AGDD5 ', met%AGDD5
    write(*,*) 'climate%GDD5 ', met%GDD5
    write(*,*) 'climate%AGDD0 ', met%AGDD0
    write(*,*) 'climate%GDD0 ', met%GDD0
    write(*,*) 'climate%alpha_PT ', met%alpha_PT
    write(*,*) 'climate%evap_PT ', met%evap_PT
    write(*,*) 'climate%aevap ', met%aevap
    write(*,*) 'climate%alpha_PT20 ', met%alpha_PT20
    write(*,*) 'climate%GDD0_rec ', met%GDD0_rec
    write(*,*) 'climate%frec ', met%frec
    write(*,*) 'climate%dtemp_min ', met%dtemp_min
    write(*,*) 'climate%fdorm ', met%fdorm
    write(*,*) 'climate%fapar_ann_max ', met%fapar_ann_max
    write(*,*) 'climate%fapar_ann_max_last_year ', met%fapar_ann_max_last_year
    write(*,*) 'climate%AvgAnnMaxFAPAR ', met%AvgAnnMaxFAPAR
    write(*,*) 'climate%dtemp_max ', met%dtemp_max
    write(*,*) 'climate%drhum ', met%drhum
    write(*,*) 'climate%du10_max ', met%du10_max
    write(*,*) 'climate%dprecip ', met%dprecip
    write(*,*) 'climate%aprecip ', met%aprecip
    write(*,*) 'climate%aprecip_av20 ', met%aprecip_av20
    write(*,*) 'climate%last_precip ', met%last_precip
    write(*,*) 'climate%KBDI ', met%KBDI
    write(*,*) 'climate%FFDI ', met%FFDI
    write(*,*) 'climate%D_MacArthur ', met%D_MacArthur
    write(*,*) 'climate%Nesterov_Current ', met%Nesterov_Current
    write(*,*) 'climate%Nesterov_ann_max ', met%Nesterov_ann_max
    write(*,*) 'climate%Nesterov_ann_max_last_year ', met%Nesterov_ann_max_last_year
    write(*,*) 'climate%Nesterov_ann_running_max ', met%Nesterov_ann_running_max
    ! commented because of large output: uncomment if needed
    ! write(*,*) 'climate%mtemp_min_20 ', met%mtemp_min_20
    ! write(*,*) 'climate%mtemp_max_20 ', met%mtemp_max_20
    ! write(*,*) 'climate%dmoist_min_20 ', met%dmoist_min_20
    ! write(*,*) 'climate%dmoist_max_20 ', met%dmoist_max_20
    ! write(*,*) 'climate%dtemp_31 ', met%dtemp_31
    ! write(*,*) 'climate%dmoist_31 ', met%dmoist_31
    ! write(*,*) 'climate%alpha_PT_20 ', met%alpha_PT_20
    ! write(*,*) 'climate%dtemp_91 ', met%dtemp_91
    ! write(*,*) 'climate%APAR_leaf_sun ', met%APAR_leaf_sun
    ! write(*,*) 'climate%APAR_leaf_shade ', met%APAR_leaf_shade
    ! write(*,*) 'climate%Dleaf_sun ', met%Dleaf_sun
    ! write(*,*) 'climate%Dleaf_shade ', met%Dleaf_shade
    ! write(*,*) 'climate%Tleaf_sun ', met%Tleaf_sun
    ! write(*,*) 'climate%Tleaf_shade ', met%Tleaf_shade
    ! write(*,*) 'climate%cs_sun ', met%cs_sun
    ! write(*,*) 'climate%cs_shade ', met%cs_shade
    ! write(*,*) 'climate%scalex_sun ', met%scalex_sun
    ! write(*,*) 'climate%scalex_shade ', met%scalex_shade
    ! write(*,*) 'climate%fwsoil ', met%fwsoil
    ! write(*,*) 'climate%aprecip_20 ', met%aprecip_20
    ! write(*,*) 'climate%Rd_sun ', met%Rd_sun
    ! write(*,*) 'climate%Rd_shade ', met%Rd_shade

  end subroutine print_climate_type


  ! ------------------------------------------------------------------


  subroutine print_sum_flux_type(sum_flux)

    type(sum_flux_type), intent(in) :: sum_flux

    write(*,*) 'sumflux%sumpn ', sum_flux%sumpn
    write(*,*) 'sumflux%sumrp ', sum_flux%sumrp
    write(*,*) 'sumflux%sumrpw ', sum_flux%sumrpw
    write(*,*) 'sumflux%sumrpr ', sum_flux%sumrpr
    write(*,*) 'sumflux%sumrs ', sum_flux%sumrs
    write(*,*) 'sumflux%sumrd ', sum_flux%sumrd
    write(*,*) 'sumflux%dsumpn ', sum_flux%dsumpn
    write(*,*) 'sumflux%dsumrp ', sum_flux%dsumrp
    write(*,*) 'sumflux%dsumrs ', sum_flux%dsumrs
    write(*,*) 'sumflux%dsumrd ', sum_flux%dsumrd
    write(*,*) 'sumflux%sumxrp ', sum_flux%sumxrp
    write(*,*) 'sumflux%sumxrs ', sum_flux%sumxrs

  end subroutine print_sum_flux_type


  ! ------------------------------------------------------------------


  subroutine print_bgc_pool_type(bgc)

    type(bgc_pool_type), intent(in) :: bgc

    write(*,*) 'bgc%cplant ', bgc%cplant
    write(*,*) 'bgc%csoil ', bgc%csoil

  end subroutine print_bgc_pool_type


  ! ------------------------------------------------------------------


  subroutine nc_err(status, ivar)

    use netcdf, only: nf90_noerr, nf90_strerror
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif

    integer, intent(in)              :: status
    integer, intent(inout), optional :: ivar

#ifdef __MPI__
    integer :: ierr
#endif

    if (status /= nf90_noerr) then
       write(*,*) "netCDF error"
       write(*,*) trim(nf90_strerror(status))
#ifdef __MPI__
       call MPI_Abort(0, 171, ierr)
#else
       stop 171
#endif
    else
       if (present(ivar)) ivar = ivar + 1
    end if

  end subroutine nc_err


  ! ------------------------------------------------------------------


  subroutine read_netcdf_climate_type(filename, climate)

    use netcdf, only: nf90_open, nf90_nowrite, &
         nf90_inq_varid, nf90_get_var, nf90_close
#ifdef __MPI__
    use mpi, only: MPI_Abort
#endif

    implicit none

    character(len=*),   intent(in)  :: filename
    type(climate_type), intent(inout) :: climate

    logical :: existfile
    integer :: fid, vid
#ifdef __MPI__
    integer :: ierr
#endif

    ! open netCDF file
    inquire(file=trim(filename), exist=existfile)
    if (.not. existfile) then
       write(*,*) filename, ' does not exist!'
#ifdef __MPI__
       call MPI_Abort(0, 172, ierr)
#else
       stop 172
#endif
    endif

    ! open netCDF file
    call nc_err(nf90_open(trim(filename), nf90_nowrite, fid))

    ! read variables
    ! integer scalars
    call nc_err(nf90_inq_varid(fid, 'nyear_average', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nyear_average))
    call nc_err(nf90_inq_varid(fid, 'nday_average', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nday_average))
    call nc_err(nf90_inq_varid(fid, 'nyears', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nyears))
    call nc_err(nf90_inq_varid(fid, 'doy', vid))
    call nc_err(nf90_get_var(fid, vid, climate%doy))
    ! integer vectors
    call nc_err(nf90_inq_varid(fid, 'chilldays', vid))
    call nc_err(nf90_get_var(fid, vid, climate%chilldays))
    call nc_err(nf90_inq_varid(fid, 'iveg', vid))
    call nc_err(nf90_get_var(fid, vid, climate%iveg))
    call nc_err(nf90_inq_varid(fid, 'biome', vid))
    call nc_err(nf90_get_var(fid, vid, climate%biome))
    call nc_err(nf90_inq_varid(fid, 'gmd', vid))
    call nc_err(nf90_get_var(fid, vid, climate%gmd))
    call nc_err(nf90_inq_varid(fid, 'modis_igbp', vid))
    call nc_err(nf90_get_var(fid, vid, climate%modis_igbp))
    call nc_err(nf90_inq_varid(fid, 'dslr', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dslr))
    call nc_err(nf90_inq_varid(fid, 'nday_nesterov', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nday_nesterov))
    ! real vectors
    call nc_err(nf90_inq_varid(fid, 'dtemp', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dtemp))
    call nc_err(nf90_inq_varid(fid, 'dmoist', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist))
    call nc_err(nf90_inq_varid(fid, 'dmoist_min', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_min))
    call nc_err(nf90_inq_varid(fid, 'dmoist_min20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_min20))
    call nc_err(nf90_inq_varid(fid, 'dmoist_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_max))
    call nc_err(nf90_inq_varid(fid, 'dmoist_max20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_max20))
    call nc_err(nf90_inq_varid(fid, 'mtemp', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp))
    call nc_err(nf90_inq_varid(fid, 'qtemp', vid))
    call nc_err(nf90_get_var(fid, vid, climate%qtemp))
    call nc_err(nf90_inq_varid(fid, 'mmoist', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mmoist))
    call nc_err(nf90_inq_varid(fid, 'mtemp_min', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp_min))
    call nc_err(nf90_inq_varid(fid, 'mtemp_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp_max))
    call nc_err(nf90_inq_varid(fid, 'qtemp_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%qtemp_max))
    call nc_err(nf90_inq_varid(fid, 'qtemp_max_last_year', vid))
    call nc_err(nf90_get_var(fid, vid, climate%qtemp_max_last_year))
    call nc_err(nf90_inq_varid(fid, 'mtemp_min20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp_min20))
    call nc_err(nf90_inq_varid(fid, 'mtemp_max20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp_max20))
    call nc_err(nf90_inq_varid(fid, 'atemp_mean', vid))
    call nc_err(nf90_get_var(fid, vid, climate%atemp_mean))
    call nc_err(nf90_inq_varid(fid, 'agdd5', vid))
    call nc_err(nf90_get_var(fid, vid, climate%agdd5))
    call nc_err(nf90_inq_varid(fid, 'gdd5', vid))
    call nc_err(nf90_get_var(fid, vid, climate%gdd5))
    call nc_err(nf90_inq_varid(fid, 'agdd0', vid))
    call nc_err(nf90_get_var(fid, vid, climate%agdd0))
    call nc_err(nf90_inq_varid(fid, 'gdd0', vid))
    call nc_err(nf90_get_var(fid, vid, climate%gdd0))
    call nc_err(nf90_inq_varid(fid, 'alpha_pt', vid))
    call nc_err(nf90_get_var(fid, vid, climate%alpha_pt))
    call nc_err(nf90_inq_varid(fid, 'evap_pt', vid))
    call nc_err(nf90_get_var(fid, vid, climate%evap_pt))
    call nc_err(nf90_inq_varid(fid, 'aevap', vid))
    call nc_err(nf90_get_var(fid, vid, climate%aevap))
    call nc_err(nf90_inq_varid(fid, 'alpha_pt20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%alpha_pt20))
    call nc_err(nf90_inq_varid(fid, 'gdd0_rec', vid))
    call nc_err(nf90_get_var(fid, vid, climate%gdd0_rec))
    call nc_err(nf90_inq_varid(fid, 'frec', vid))
    call nc_err(nf90_get_var(fid, vid, climate%frec))
    call nc_err(nf90_inq_varid(fid, 'dtemp_min', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dtemp_min))
    call nc_err(nf90_inq_varid(fid, 'fdorm', vid))
    call nc_err(nf90_get_var(fid, vid, climate%fdorm))
    call nc_err(nf90_inq_varid(fid, 'fapar_ann_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%fapar_ann_max))
    call nc_err(nf90_inq_varid(fid, 'fapar_ann_max_last_year', vid))
    call nc_err(nf90_get_var(fid, vid, climate%fapar_ann_max_last_year))
    call nc_err(nf90_inq_varid(fid, 'avgannmaxfapar', vid))
    call nc_err(nf90_get_var(fid, vid, climate%avgannmaxfapar))
    call nc_err(nf90_inq_varid(fid, 'dtemp_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dtemp_max))
    call nc_err(nf90_inq_varid(fid, 'drhum', vid))
    call nc_err(nf90_get_var(fid, vid, climate%drhum))
    call nc_err(nf90_inq_varid(fid, 'du10_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%du10_max))
    call nc_err(nf90_inq_varid(fid, 'dprecip', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dprecip))
    call nc_err(nf90_inq_varid(fid, 'aprecip', vid))
    call nc_err(nf90_get_var(fid, vid, climate%aprecip))
    call nc_err(nf90_inq_varid(fid, 'aprecip_av20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%aprecip_av20))
    call nc_err(nf90_inq_varid(fid, 'last_precip', vid))
    call nc_err(nf90_get_var(fid, vid, climate%last_precip))
    call nc_err(nf90_inq_varid(fid, 'kbdi', vid))
    call nc_err(nf90_get_var(fid, vid, climate%kbdi))
    call nc_err(nf90_inq_varid(fid, 'ffdi', vid))
    call nc_err(nf90_get_var(fid, vid, climate%ffdi))
    call nc_err(nf90_inq_varid(fid, 'd_macarthur', vid))
    call nc_err(nf90_get_var(fid, vid, climate%d_macarthur))
    call nc_err(nf90_inq_varid(fid, 'nesterov_current', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nesterov_current))
    call nc_err(nf90_inq_varid(fid, 'nesterov_ann_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nesterov_ann_max))
    call nc_err(nf90_inq_varid(fid, 'nesterov_ann_max_last_year', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nesterov_ann_max_last_year))
    call nc_err(nf90_inq_varid(fid, 'nesterov_ann_running_max', vid))
    call nc_err(nf90_get_var(fid, vid, climate%nesterov_ann_running_max))
    ! real arrays, [dim1, dim2]
    call nc_err(nf90_inq_varid(fid, 'mtemp_min_20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp_min_20))
    call nc_err(nf90_inq_varid(fid, 'mtemp_max_20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%mtemp_max_20))
    call nc_err(nf90_inq_varid(fid, 'dmoist_min_20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_min_20))
    call nc_err(nf90_inq_varid(fid, 'dmoist_max_20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_max_20))
    call nc_err(nf90_inq_varid(fid, 'alpha_pt_20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%alpha_pt_20))
    call nc_err(nf90_inq_varid(fid, 'aprecip_20', vid))
    call nc_err(nf90_get_var(fid, vid, climate%aprecip_20))
    ! real arrays, [dim1, dim3]
    call nc_err(nf90_inq_varid(fid, 'dtemp_31', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dtemp_31))
    call nc_err(nf90_inq_varid(fid, 'dmoist_31', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dmoist_31))
    ! real arrays, [dim1, dim4]
    call nc_err(nf90_inq_varid(fid, 'dtemp_91', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dtemp_91))
    ! real arrays, [dim1, dim5]
    call nc_err(nf90_inq_varid(fid, 'apar_leaf_sun', vid))
    call nc_err(nf90_get_var(fid, vid, climate%apar_leaf_sun))
    call nc_err(nf90_inq_varid(fid, 'apar_leaf_shade', vid))
    call nc_err(nf90_get_var(fid, vid, climate%apar_leaf_shade))
    call nc_err(nf90_inq_varid(fid, 'dleaf_sun', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dleaf_sun))
    call nc_err(nf90_inq_varid(fid, 'dleaf_shade', vid))
    call nc_err(nf90_get_var(fid, vid, climate%dleaf_shade))
    call nc_err(nf90_inq_varid(fid, 'tleaf_sun', vid))
    call nc_err(nf90_get_var(fid, vid, climate%tleaf_sun))
    call nc_err(nf90_inq_varid(fid, 'tleaf_shade', vid))
    call nc_err(nf90_get_var(fid, vid, climate%tleaf_shade))
    call nc_err(nf90_inq_varid(fid, 'cs_sun', vid))
    call nc_err(nf90_get_var(fid, vid, climate%cs_sun))
    call nc_err(nf90_inq_varid(fid, 'cs_shade', vid))
    call nc_err(nf90_get_var(fid, vid, climate%cs_shade))
    call nc_err(nf90_inq_varid(fid, 'scalex_sun', vid))
    call nc_err(nf90_get_var(fid, vid, climate%scalex_sun))
    call nc_err(nf90_inq_varid(fid, 'scalex_shade', vid))
    call nc_err(nf90_get_var(fid, vid, climate%scalex_shade))
    call nc_err(nf90_inq_varid(fid, 'fwsoil', vid))
    call nc_err(nf90_get_var(fid, vid, climate%fwsoil))
    call nc_err(nf90_inq_varid(fid, 'rd_sun', vid))
    call nc_err(nf90_get_var(fid, vid, climate%rd_sun))
    call nc_err(nf90_inq_varid(fid, 'rd_shade', vid))
    call nc_err(nf90_get_var(fid, vid, climate%rd_shade))

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine read_netcdf_climate_type


  ! ------------------------------------------------------------------


  subroutine write_netcdf_climate_type(filename, climate)

    use netcdf, only: nf90_create, nf90_clobber, nf90_64bit_offset, &
         nf90_def_dim, nf90_def_var, nf90_int, nf90_float, nf90_enddef, &
         nf90_put_var, nf90_close

    implicit none

    character(len=*),   intent(in) :: filename
    type(climate_type), intent(in) :: climate

    integer :: fid
    integer :: dimid1, dimid2, dimid3, dimid4, dimid5
    integer :: i
    integer, dimension(78) :: vid

    ! create netCDF file
    call nc_err(nf90_create(trim(filename), ior(nf90_clobber, nf90_64bit_offset), fid))

    ! define dimensions
    ! land
    call nc_err(nf90_def_dim(fid, 'dim1', size(climate%dtemp, 1), dimid1))
    ! number of years (stored for 20 yr running means)
    call nc_err(nf90_def_dim(fid, 'dim2', size(climate%mtemp_min_20, 2), dimid2))
    ! number of days (stored for 31 day monthly means)
    call nc_err(nf90_def_dim(fid, 'dim3', size(climate%dtemp_31, 2), dimid3))
    ! number of days (stored for 91 day quarterly means)
    call nc_err(nf90_def_dim(fid, 'dim4', size(climate%dtemp_91, 2), dimid4))
    ! number of 5 days of sub-diurnal time-steps (stored for leaf photosynthesis drivers)
    call nc_err(nf90_def_dim(fid, 'dim5', size(climate%fwsoil, 2), dimid5))

    ! define variables
    i = 1
    ! define integer scalar variables
    call nc_err(nf90_def_var(fid, 'nyear_average', nf90_int, vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nday_average', nf90_int, vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nyears', nf90_int, vid(i)), i)
    call nc_err(nf90_def_var(fid, 'doy', nf90_int, vid(i)), i)

    ! define integer vector variables
    call nc_err(nf90_def_var(fid, 'chilldays', nf90_int, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'iveg', nf90_int, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'biome', nf90_int, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'gmd', nf90_int, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'modis_igbp', nf90_int, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dslr', nf90_int, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nday_nesterov', nf90_int, [dimid1], vid(i)), i)

    ! define real vector variables
    call nc_err(nf90_def_var(fid, 'dtemp', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_min', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_min20', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_max20', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtemp', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'qtemp', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mmoist', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtemp_min', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtemp_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'qtemp_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'qtemp_max_last_year', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtemp_min20', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtemp_max20', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'atemp_mean', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'agdd5', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'gdd5', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'agdd0', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'gdd0', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'alpha_pt', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'evap_pt', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'aevap', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'alpha_pt20', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'gdd0_rec', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'frec', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dtemp_min', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fdorm', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fapar_ann_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fapar_ann_max_last_year', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'avgannmaxfapar', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dtemp_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'drhum', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'du10_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dprecip', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'aprecip', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'aprecip_av20', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'last_precip', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'kbdi', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'ffdi', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'd_macarthur', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nesterov_current', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nesterov_ann_max', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nesterov_ann_max_last_year', nf90_float, [dimid1], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'nesterov_ann_running_max', nf90_float, [dimid1], vid(i)), i)

    ! define real array variables, [dim1, dim2]
    call nc_err(nf90_def_var(fid, 'mtemp_min_20', nf90_float, [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'mtemp_max_20', nf90_float, [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_min_20', nf90_float, [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_max_20', nf90_float, [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'alpha_pt_20', nf90_float, [dimid1, dimid2], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'aprecip_20', nf90_float, [dimid1, dimid2], vid(i)), i)

    ! define real array variables, [dim1, dim3]
    call nc_err(nf90_def_var(fid, 'dtemp_31', nf90_float, [dimid1, dimid3], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dmoist_31', nf90_float, [dimid1, dimid3], vid(i)), i)

    ! define real array variables, [dim1, dim4]
    call nc_err(nf90_def_var(fid, 'dtemp_91', nf90_float, [dimid1, dimid4], vid(i)), i)

    ! define real array variables, [dim1, dim5]
    call nc_err(nf90_def_var(fid, 'apar_leaf_sun', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'apar_leaf_shade', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dleaf_sun', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'dleaf_shade', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tleaf_sun', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'tleaf_shade', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cs_sun', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'cs_shade', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'scalex_sun', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'scalex_shade', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'fwsoil', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rd_sun', nf90_float, [dimid1, dimid5], vid(i)), i)
    call nc_err(nf90_def_var(fid, 'rd_shade', nf90_float, [dimid1, dimid5], vid(i)), i)

    ! end define mode
    call nc_err(nf90_enddef(fid))

    ! put variables
    i = 1
    ! integer scalars
    call nc_err(nf90_put_var(fid, vid(i), climate%nyear_average), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nday_average), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nyears), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%doy), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%chilldays), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%iveg), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%biome), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%gmd), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%modis_igbp), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dslr), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nday_nesterov), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dtemp), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_min), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_min20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_max20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%qtemp), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mmoist), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp_min), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%qtemp_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%qtemp_max_last_year), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp_min20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp_max20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%atemp_mean), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%agdd5), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%gdd5), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%agdd0), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%gdd0), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%alpha_pt), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%evap_pt), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%aevap), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%alpha_pt20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%gdd0_rec), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%frec), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dtemp_min), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%fdorm), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%fapar_ann_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%fapar_ann_max_last_year), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%avgannmaxfapar), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dtemp_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%drhum), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%du10_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dprecip), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%aprecip), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%aprecip_av20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%last_precip), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%kbdi), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%ffdi), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%d_macarthur), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nesterov_current), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nesterov_ann_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nesterov_ann_max_last_year), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%nesterov_ann_running_max), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp_min_20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%mtemp_max_20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_min_20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_max_20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%alpha_pt_20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%aprecip_20), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dtemp_31), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dmoist_31), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dtemp_91), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%apar_leaf_sun), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%apar_leaf_shade), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dleaf_sun), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%dleaf_shade), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%tleaf_sun), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%tleaf_shade), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%cs_sun), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%cs_shade), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%scalex_sun), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%scalex_shade), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%fwsoil), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%rd_sun), i)
    call nc_err(nf90_put_var(fid, vid(i), climate%rd_shade), i)

    ! close NetCDF file
    call nc_err(nf90_close(fid))

  end subroutine write_netcdf_climate_type

end module cable_def_types_mod
