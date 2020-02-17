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

   INTEGER :: mp,    & ! # total no of patches/tiles
              mvtype,& ! total # vegetation types,   from input
              mstype,& ! total # soil types,         from input
              mland                           ! # land grid cells

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

      REAL, DIMENSION(:), POINTER ::                                           &
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

      INTEGER, DIMENSION(:), POINTER ::                                        &
         isoilm     ! integer soil type

      REAL, DIMENSION(:), POINTER ::                                           &
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

      REAL(r_2), DIMENSION(:), POINTER ::                                      &
         cnsd => null(),    & ! thermal conductivity of dry soil [W/m/K]
         pwb_min => null()    ! working variable (swilt/ssat)**ibp2

      REAL, DIMENSION(:,:), POINTER ::                                         &
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

      REAL, DIMENSION(:), POINTER ::                                           &
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
         totenbal => null(),&!
         totenbal2 => null(),&
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

      REAL, DIMENSION(:,:), POINTER ::                                         &
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

      REAL(r_2), DIMENSION(:), POINTER ::                                      &
         wbtot => null()   ! total soil water (mm)

      REAL(r_2), DIMENSION(:,:), POINTER ::                                    &
         gammzz => null(),  & ! heat capacity for each soil layer
         wb => null(),      & ! volumetric soil moisture (solid+liq)
         wbice => null(),   & ! soil ice
         wblf => null(),    & !
         wbfice => null()     !

     ! Additional SLI variables:
     REAL(r_2), DIMENSION(:,:), POINTER :: S => null()         ! moisture content relative to sat value    (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil => null()         !     Tsoil (deg C)
     REAL(r_2), DIMENSION(:),   POINTER :: SL => null()        ! litter moisture content relative to sat value (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:),   POINTER :: TL => null()        ! litter temperature in K     (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:),   POINTER :: h0 => null()        ! pond height in m            (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: rex => null()       ! root extraction from each layer (mm/dels)
     REAL(r_2), DIMENSION(:,:), POINTER :: wflux => null()     ! water flux at layer boundaries (mm s-1)
     REAL(r_2), DIMENSION(:),   POINTER :: delwcol => null()   ! change in water column (mm / dels)
     REAL(r_2), DIMENSION(:),   POINTER :: zdelta => null()    ! water table depth           (edit vh 23/06/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: kth => null()       ! thermal conductivity           (edit vh 29/07/08)
     REAL(r_2), DIMENSION(:),   POINTER :: Tsurface => null()  !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:),   POINTER :: lE => null()        ! soil latent heat flux
     REAL(r_2), DIMENSION(:),   POINTER :: evap => null()      ! soil evaporation (mm / dels)
     REAL(r_2), DIMENSION(:,:), POINTER :: ciso => null()      ! concentration of minor isotopologue in soil water (kg m-3 water)
     REAL(r_2), DIMENSION(:),   POINTER :: cisoL => null()     ! concentration of minor isotopologue in litter water (kg m-3 water)
     REAL(r_2), DIMENSION(:),   POINTER :: rlitt => null()     ! resistance to heat/moisture transfer through litter (m-1 s)
     REAL(r_2), DIMENSION(:,:), POINTER :: thetai => null()    ! volumetric ice content (MC)
     REAL(r_2), DIMENSION(:,:), POINTER :: snowliq => null()   ! liquid snow content (mm H2O)
     REAL(r_2), DIMENSION(:),   POINTER :: nsteps => null()    ! number of iterations at each timestep
     REAL(r_2), DIMENSION(:),   POINTER :: TsurfaceFR => null()  !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: Ta_daily => null()        ! air temp averaged over last 24h
     INTEGER, DIMENSION(:),     POINTER :: nsnow => null() ! number of layers in snow-pack (0-nsnow_max)
     REAL(r_2), DIMENSION(:),   POINTER :: Qadv_daily => null()  ! advective heat flux into surface , daily average (W m-2)
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
         gmmax => null()      ! max. mesophyll conductance at 25degC top leaf

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

      REAL, DIMENSION(:), POINTER ::                                           &
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
               
      REAL, DIMENSION(:,:), POINTER ::                                         &
         evapfbl => null(), &
         gswx => null(),    & ! stom cond for water
         zetar => null(),   & ! stability parameter (ref height)
         !! vh_js !!
         zetash => null()     ! stability parameter (shear height)

      REAL(r_2), DIMENSION(:), POINTER ::                                      &
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
     REAL(r_2), DIMENSION(:), POINTER :: kthLitt => null(), DvLitt => null()

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

      REAL, DIMENSION(:), POINTER   ::                                         &
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

      REAL, DIMENSION(:,:), POINTER  ::                                        &
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

      REAL, DIMENSION(:,:,:), POINTER ::                                       &
         qcan => null() ! absorbed radiation for canopy (W/m^2)

  END TYPE radiation_type

  ! .............................................................................

   ! Roughness variables:
   TYPE roughness_type

      REAL, DIMENSION(:), POINTER ::                                           &
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
      REAL, DIMENSION(:), POINTER ::                                           &
         coexp => null() ! Extinction coef for wind profile in canopy

      ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
      REAL, DIMENSION(:), POINTER ::                                           &
         usuh => null() ! Friction velocity/windspeed at canopy height

      REAL, DIMENSION(:), POINTER ::                                           &
         term2 => null(), term3 => null(), term5 => null(), term6 => null(), term6a => null() ! for aerodyn resist. calc.

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

      INTEGER, DIMENSION(:), POINTER ::                                        &
         year => null(),    & ! local time year AD
         moy => null()        ! local time month of year

      REAL, DIMENSION(:), POINTER ::                                           &
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
      REAL, DIMENSION(:,:), POINTER ::                                         &
         fsd => null()  ! downward short-wave radiation (W/m2)

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
       chilldays => null(), &   ! length of chilling period (period with T<5deg)
       iveg => null(), &        ! potential vegetation type based on climatic constraints
       biome => null(), &
       GMD => null(), &           ! growing moisture days (== number days since min moisture threshold)
       modis_igbp => null(), &    ! IGBP biome classification
       DSLR => null(), &           ! days since last rain
       NDAY_Nesterov => null()
 

      REAL, DIMENSION(:), POINTER ::                                           &
      dtemp => null(),        & ! daily mean temperature
      dmoist => null(),        & ! daily moisture availability
      dmoist_min => null(),        & ! minimum daily moisture availability over the year
      dmoist_min20 => null(),        & ! min daily moisture avail over the year, averaged over 20 y
      dmoist_max => null(),        & ! maximum daily moisture availability over the year
      dmoist_max20 => null(),        & ! max daily moisture avail over the year, averaged over 20 y
      mtemp => null(),       & ! mean temperature over the last 31 days
      qtemp => null(),       & ! mean temperature over the last 91 days
      mmoist => null(),        & ! monthly moisture availability
      mtemp_min => null(),   & ! minimum monthly temperature
      mtemp_max => null(),   & ! maximum monhtly temperature
      qtemp_max => null(),   & ! mean temperature of the warmest quarter (so far this year)
      qtemp_max_last_year => null(),   & ! mean temperature of the warmest quarter (last calendar year)
      mtemp_min20 => null(),   & ! minimum monthly temperature, averaged over 20 y
      mtemp_max20 => null(),   & ! maximum monhtly temperature, averaged over 20 y
      atemp_mean => null(),  & ! annual average temperature
      AGDD5 => null(),       &
      GDD5 => null(),        & ! growing degree day sum relative to 5deg base temperature
      AGDD0 => null(),        & ! 
      GDD0 => null(),        & ! growing degree day sum relative to 0deg base temperature
      alpha_PT => null(),    & ! ratio of annual evap to annual PT evap
      evap_PT => null(),    & ! annual PT evap [mm]
      aevap => null(), &       ! annual evap [mm]
      alpha_PT20 => null(), &   
      GDD0_rec => null(), &    ! growing degree day sum related to spring photosynthetic recovery
      frec => null(), &           ! fractional photosynthetic recovery
      dtemp_min => null(), &      ! daily minimum temperature
      fdorm => null(), & ! dormancy fraction (1 prior to first autumn frost; 0 after 10 severe frosts)
      fapar_ann_max => null(), & ! maximum midday fpar so far this year
      fapar_ann_max_last_year => null(), & ! maximum midday fpar last year
      AvgAnnMaxFAPAR => null(), &  ! average annual maximum FAPAR
      dtemp_max => null(), &       ! daily maximum temperature
      drhum => null(),      &       ! daily average relative humidity
      du10_max => null(),   &       ! daily max wind speed at 10 m
      dprecip => null(),   &            ! daily total precip (mm)
      aprecip => null(),   &            ! total precip accumulated over the current year (mm)
      aprecip_av20 => null(), &       ! annual precip averaged over the last 20 y
      last_precip => null(), &        ! rainfall accumulated since last day without rain
      KBDI => null(),  &         ! Keetch-Byram-Drought-Index (Keetch, 1968)
      FFDI => null(), &          ! Forest Fire Danger Index
      D_MacArthur => null(), &           ! MacArthur Drought Factor
      Nesterov_current => null(), &         ! current nesterov index
      Nesterov_ann_max => null(), &            ! annual maximum nesterov index (current year)
      Nesterov_ann_max_last_year  => null(), &            ! annual maximum nesterov index (last year)
      Nesterov_ann_running_max => null()

      REAL, DIMENSION(:,:), POINTER ::                                   &
      mtemp_min_20 => null(), & ! mimimum monthly temperatures for the last 20 y
      mtemp_max_20 => null(), & ! maximum monthly temperatures for the last 20 y
      dmoist_min_20 => null(), & ! min daily moisture for the last 20 y
      dmoist_max_20 => null(), & ! max daily moisture for the last 20 y
      dtemp_31  => null(), &    ! daily temperature for the last 31 days
      dmoist_31  => null(), &    ! daily moisture availability for the last 31 days
      alpha_PT_20 => null(), &      ! priestley Taylor Coefft for last 20 y
      dtemp_91 => null(), &     ! daily temperature for the last 91 days
      APAR_leaf_sun => null(), &  ! flux of PAR absrobed by sunlit leaves (subdiurnal time-step)
      APAR_leaf_shade => null(), & ! flux of PAR absrobed by shaded leaves (subdiurnal time-step)
      Dleaf_sun => null(), & ! leaf-air vapour pressure difference (sun leaves, subdiurnal time-step)
      Dleaf_shade => null(),  & ! leaf-air vapour pressure difference (shade leaves, subdiurnal time-step)
      Tleaf_sun => null(), & ! leaf T (sun leaves, subdiurnal time-step)
      Tleaf_shade => null(), &  ! leaf T  (shade leaves, subdiurnal time-step)
      cs_sun => null(), &     ! sun leaf cs (ppm CO2)
      cs_shade => null(), &          ! shade leaf cs (ppm CO2)
      scalex_sun => null(), & ! canopy depth scaling factor on vcmax and jmax (sun leaves)
      scalex_shade => null(), & ! canopy depth scaling factor on vcmax and jmax (shade leaves)
      fwsoil => null(), &         ! soil-moisture modifier to stomatal conductance
      aprecip_20 => null(), &    ! annual average rainfall for the last 20 years
      Rd_sun => null(), &
      Rd_shade => null()
   END TYPE climate_type

   ! .............................................................................

   ! Cumulative flux variables:
   TYPE sum_flux_type

      REAL, DIMENSION(:), POINTER ::                                           &
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

      REAL, DIMENSION(:,:), POINTER ::                                         &
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

   allocate( var% drybal(mp) )
   allocate( var% ebal(mp) )
   allocate( var% ebal_tot(mp) )
   allocate( var% ebaltr(mp) )
   allocate( var% ebal_tottr(mp) )
   allocate( var% ebal_cncheck(mp) )
   allocate( var% ebal_tot_cncheck(mp) )
   allocate( var% evap_tot(mp) )
   allocate( var% osnowd0(mp) )
   allocate( var% precip_tot(mp) )
   allocate( var% rnoff_tot(mp) )
   allocate( var% wbal(mp) )
   allocate( var% wbal_tot(mp) )
   allocate( var% wbtot0(mp) )
   allocate( var% wetbal(mp) )
   allocate( var% cansto0(mp) )
   allocate( var% evapc_tot(mp) )
   allocate( var% evaps_tot(mp) )
   allocate( var% rnof1_tot(mp) )
   allocate( var% rnof2_tot(mp) )
   allocate( var% snowdc_tot(mp) )
   allocate( var% wbal_tot1(mp) )
   allocate( var% owbtot(mp) )
   allocate( var% delwc_tot(mp) )
   allocate( var% qasrf_tot(mp) )
   allocate( var% qfsrf_tot(mp) )
   allocate( var% qssrf_tot(mp) )

    allocate( var% Radbal(mp) )
    allocate( var% EbalSoil(mp) )
    allocate( var% Ebalveg(mp) )
    allocate( var% Radbalsum(mp) )

END SUBROUTINE alloc_balances_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_soil_parameter_type(var, mp)

   TYPE(soil_parameter_type), INTENT(inout) :: var
   INTEGER, INTENT(in) :: mp

   allocate( var% bch(mp) )
   allocate( var% c3(mp) )
   allocate( var% clay(mp) )
   allocate( var% css(mp) )
   allocate( var% hsbh(mp) )
   allocate( var% hyds(mp) )
   allocate( var% i2bp3(mp) )
   allocate( var% ibp2(mp) )
   allocate( var% isoilm(mp) )
   allocate( var% rhosoil(mp) )
   allocate( var% sand(mp) )
   allocate( var% sfc(mp) )
   allocate( var% silt(mp) )
   allocate( var% ssat(mp) )
   allocate( var% sucs(mp) )
   allocate( var% swilt(mp) )
   allocate( var% zse(ms) )
   allocate( var% zshh(ms+1) )
   allocate( var% cnsd(mp) )
   allocate( var% albsoil(mp, nrb) )
   allocate( var% pwb_min(mp) )
   allocate( var% albsoilf(mp) )
   allocate( var% soilcol(mp) )

   ! Allocate variables for SLI soil model:
   ALLOCATE ( var % nhorizons(mp) )
   ALLOCATE ( var % ishorizon(mp,ms) )
   ALLOCATE ( var % clitt(mp) )
   ALLOCATE ( var % zeta(mp) )
   ALLOCATE ( var % fsatmax(mp) )
   ALLOCATE ( var % swilt_vec(mp,ms) )
   ALLOCATE ( var % ssat_vec(mp,ms) )
   ALLOCATE ( var % sfc_vec(mp,ms) )
   IF(.NOT.(ASSOCIATED(var % swilt_vec))) ALLOCATE ( var % swilt_vec(mp,ms) )
   IF(.NOT.(ASSOCIATED(var % ssat_vec))) ALLOCATE ( var % ssat_vec(mp,ms) )
   IF(.NOT.(ASSOCIATED(var % sfc_vec))) ALLOCATE ( var % sfc_vec(mp,ms) )


END SUBROUTINE alloc_soil_parameter_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_soil_snow_type(var, mp)

   TYPE(soil_snow_type), INTENT(inout) :: var
   INTEGER, INTENT(in) :: mp

   ALLOCATE( var%iantrct(mp) )
   ALLOCATE( var%pudsto(mp) )
   ALLOCATE( var%pudsmx(mp) )
   ALLOCATE( var%dtmlt(mp,3) )
   ALLOCATE( var%albsoilsn(mp,nrb) )
   ALLOCATE( var%cls(mp) )
   ALLOCATE( var%dfn_dtg(mp) )
   ALLOCATE( var%dfh_dtg(mp) )
   ALLOCATE( var%dfe_ddq(mp) )
   ALLOCATE( var%ddq_dtg(mp) )
   ALLOCATE( var%evapsn(mp) )
   ALLOCATE( var%fwtop(mp) )
   ALLOCATE( var%fwtop1(mp) )
   ALLOCATE( var%fwtop2(mp) )
   ALLOCATE( var%fwtop3(mp) )
   ALLOCATE( var%gammzz(mp,ms) )
   ALLOCATE( var%isflag(mp) )
   ALLOCATE( var%osnowd(mp) )
   ALLOCATE( var%potev(mp) )
   ALLOCATE( var%runoff(mp) )
   ALLOCATE( var%rnof1(mp) )
   ALLOCATE( var%rnof2(mp) )
   ALLOCATE( var%rtsoil(mp) )
   ALLOCATE( var%sconds(mp,msn) )
   ALLOCATE( var%sdepth(mp,msn) )
   ALLOCATE( var%smass(mp,msn) )
   ALLOCATE( var%snage(mp) )
   ALLOCATE( var%snowd(mp) )
   ALLOCATE( var%smelt(mp) )
   ALLOCATE( var%ssdn(mp,msn) )
   ALLOCATE( var%ssdnn(mp) )
   ALLOCATE( var%tgg(mp,ms) )
   ALLOCATE( var%tggsn(mp,msn) )
   ALLOCATE( var%tss(mp) )
   ALLOCATE( var%tss_p(mp) )
   ALLOCATE( var%deltss(mp) )
   ALLOCATE( var%owb1(mp) )
   ALLOCATE( var%wb(mp,ms) )
   ALLOCATE( var%wbice(mp,ms) )
   ALLOCATE( var%wblf(mp,ms) )
   ALLOCATE( var%wbtot(mp) )
   ALLOCATE( var%wbtot1(mp) )
   ALLOCATE( var%wbtot2(mp) )
   ALLOCATE( var%wb_lake(mp) )
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

    ! Allocate variables for SLI soil model:
    !IF(cable_user%SOIL_STRUC=='sli') THEN
    ALLOCATE( var%S(mp,ms) )
    ALLOCATE( var%Tsoil(mp,ms) )
    ALLOCATE( var%SL(mp) )
    ALLOCATE( var%TL(mp) )
    ALLOCATE( var%h0(mp) )
    ALLOCATE( var%rex(mp,ms) )
    ALLOCATE( var%wflux(mp,0:ms) )
    ALLOCATE( var%delwcol(mp) )
    ALLOCATE( var%zdelta(mp) )
    ALLOCATE( var%kth(mp,ms) )
    ALLOCATE( var%Tsurface(mp) )
    ALLOCATE( var%lE(mp) )
    ALLOCATE( var%evap(mp) )
    ALLOCATE( var%ciso(mp,ms+1) )
    ALLOCATE( var%cisoL(mp) )
    ALLOCATE( var%rlitt(mp) )
    ALLOCATE( var%thetai(mp,ms) )
    ALLOCATE( var%snowliq(mp,3) )
    ALLOCATE( var%nsteps(mp) )
    ALLOCATE( var%nsnow(mp) )
    ALLOCATE( var%TsurfaceFR(mp) )
    ALLOCATE( var%Ta_daily(mp,100))
    ALLOCATE( var%Qadv_daily(mp) )
    ALLOCATE( var%G0_daily(mp) )
    ALLOCATE( var%Qevap_daily(mp) )
    ALLOCATE( var%Qprec_daily(mp) )
    ALLOCATE( var%Qprec_snow_daily(mp) )
    ALLOCATE( var%E_fusion_sn(mp) )
    ALLOCATE( var%E_sublimation_sn(mp) )
    ALLOCATE( var%latent_heat_sn(mp) )
    ALLOCATE( var%evap_liq_sn(mp) )
    ALLOCATE( var%surface_melt(mp) )
    ALLOCATE( var%Qadv_rain_sn(mp))

    !END IF

END SUBROUTINE alloc_soil_snow_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_veg_parameter_type(var, mp)

   TYPE(veg_parameter_type), INTENT(inout) :: var
   INTEGER, INTENT(in) :: mp

   ALLOCATE( var% canst1(mp) )
   ALLOCATE( var% dleaf(mp) )
   ALLOCATE( var% ejmax(mp) )
   ALLOCATE( var% ejmax_shade(mp) )
   ALLOCATE( var% ejmax_sun(mp) )
   ALLOCATE( var% iveg(mp) )
   ALLOCATE( var% ivegp(mp) )
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
   ALLOCATE( var% vcmax_shade(mp) )
   ALLOCATE( var% vcmax_sun(mp) )
   ALLOCATE( var% xfang(mp) )
   ALLOCATE( var%extkn(mp) )
   ALLOCATE( var%wai(mp) )
   ALLOCATE( var%deciduous(mp) )
   ALLOCATE( var%froot(mp,ms) )
   !was nrb(=3), but never uses (:,3) in model
   ALLOCATE( var%refl(mp,2) ) !jhan:swb?
   ALLOCATE( var%taul(mp,2) )
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
   ALLOCATE( var%vcmaxcc(mp) )
   ALLOCATE( var%ejmaxcc(mp) )
   ALLOCATE( var%gmmax(mp) )


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
   ALLOCATE( var% spill(mp) )
   ALLOCATE( var% tscrn(mp) )
   ALLOCATE( var% wcint(mp) )
   ALLOCATE( var% tv(mp) )
   ALLOCATE( var% us(mp) )
   ALLOCATE( var% uscrn(mp) )
   ALLOCATE( var% rghlai(mp) )
   ALLOCATE( var% vlaiw(mp) )
   ALLOCATE( var% fwet(mp) )
   ALLOCATE( var% A_sh(mp) )
   ALLOCATE( var% A_sl(mp) )
   ALLOCATE( var% A_slC(mp) )
   ALLOCATE( var% A_shC(mp) )
   ALLOCATE( var% A_slJ(mp) )
   ALLOCATE( var% A_shJ(mp) )
   ALLOCATE( var% GPP_sh(mp))
   ALLOCATE( var% GPP_sl(mp))
   ALLOCATE( var% fevc_sh(mp))
   ALLOCATE( var% fevc_sl(mp))
   ALLOCATE( var% eta_GPP_cs(mp) )
   ALLOCATE( var% eta_fevc_cs(mp) )
   ALLOCATE( var% eta_A_cs(mp) )
   ALLOCATE( var% eta_A_cs_sh(mp) )
   ALLOCATE( var% eta_A_cs_sl(mp) )
   ALLOCATE( var% eta_fevc_cs_sh(mp) )
   ALLOCATE( var% eta_fevc_cs_sl(mp) )
   ALLOCATE( var% cs(mp) )
   ALLOCATE( var% dAdcs(mp) )
   ALLOCATE( var% cs_sl(mp) )
   ALLOCATE( var% cs_sh(mp) )
   ! ALLOCATE( var% ci_sl(mp) )
   ! ALLOCATE( var% ci_sh(mp) )
   ALLOCATE( var% tlf(mp) )
   ALLOCATE( var% dlf(mp) )

   ALLOCATE( var% evapfbl(mp,ms) )
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
   ALLOCATE( var% zetar(mp,NITER) )
   ALLOCATE( var% zetash(mp,NITER) )
   ALLOCATE( var % fwsoil(mp) )
   ALLOCATE( var % ofes(mp) )

   ALLOCATE( var % gw(mp,mf) )     ! dry canopy conductance (ms-1) edit vh 6/7/09
   ALLOCATE( var % ancj(mp,mf,3) ) ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
   ALLOCATE( var % tlfy(mp,mf) )   ! sunlit and shaded leaf temperatures
   ALLOCATE( var % ecy(mp,mf) )    ! sunlit and shaded leaf transpiration (dry canopy)
   ALLOCATE( var % ecx(mp,mf) )    ! sunlit and shaded leaf latent heat flux
   ! ALLOCATE( var % ci(mp,mf,3) )   ! intra-cellular CO2 vh 6/7/09
   ALLOCATE( var % fwsoil (mp) )

   ! vh_js - litter resistances to heat and vapour transfer
   ALLOCATE(var % kthLitt(mp))
   ALLOCATE(var % DvLitt(mp))

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
   ALLOCATE ( var % rhum(mp) )
   ALLOCATE ( var % u10(mp) )
END SUBROUTINE alloc_met_type

! ------------------------------------------------------------------------------

SUBROUTINE alloc_climate_type(var, mp, ktauday)

   IMPLICIT NONE

   TYPE(climate_type), INTENT(inout) :: var
   INTEGER, INTENT(in) :: mp, ktauday
   INTEGER :: ny, nd
   ny = var%nyear_average
   nd = var%nday_average

!   ALLOCATE ( var %  nyears )
!   ALLOCATE ( var %  doy )
   ALLOCATE ( var %  dtemp(mp) )
   ALLOCATE ( var %  dmoist(mp) )
   ALLOCATE ( var %  dmoist_min(mp) )
   ALLOCATE ( var %  dmoist_min20(mp) )
   ALLOCATE ( var %  dmoist_max(mp) )
   ALLOCATE ( var %  dmoist_max20(mp) )
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
   ALLOCATE ( var % GMD(mp) )
   ALLOCATE ( var % alpha_PT(mp) )
   ALLOCATE ( var % alpha_PT20(mp) )
   ALLOCATE ( var % evap_PT(mp) )
   ALLOCATE ( var % aevap(mp) )
   ALLOCATE ( var % GDD0_rec(mp) )
   ALLOCATE ( var % frec(mp) )
   ALLOCATE ( var % dtemp_min(mp) )
   ALLOCATE ( var % fdorm(mp) )
   ALLOCATE ( var % fapar_ann_max(mp) )
   ALLOCATE ( var % fapar_ann_max_last_year(mp) )
   ALLOCATE ( var % modis_igbp(mp) )
   ALLOCATE ( var % DSLR(mp) )
   ALLOCATE ( var % NDAY_Nesterov(mp) )
   ALLOCATE ( var % AvgAnnMaxFAPAR(mp) )
   ALLOCATE ( var % dtemp_max(mp) )
   ALLOCATE ( var % drhum(mp) )
   ALLOCATE ( var % du10_max(mp) )
   ALLOCATE ( var % dprecip(mp) )
   ALLOCATE ( var % aprecip(mp) )
   ALLOCATE ( var % aprecip_av20(mp) )
   ALLOCATE ( var % last_precip(mp) )
   ALLOCATE ( var % KBDI(mp) )
   ALLOCATE ( var % FFDI(mp) )
   ALLOCATE ( var % D_MacArthur(mp) )
   ALLOCATE ( var % Nesterov_Current(mp) )
   ALLOCATE ( var % Nesterov_ann_max(mp) )
   ALLOCATE ( var % Nesterov_ann_max_last_year(mp) )
   ALLOCATE ( var % Nesterov_ann_running_max(mp) )

   ALLOCATE ( var % mtemp_min_20(mp,ny) )
   ALLOCATE ( var %     mtemp_max_20(mp,ny) )
   ALLOCATE ( var % dmoist_min_20(mp,ny) )
   ALLOCATE ( var % dmoist_max_20(mp,ny) )
   ALLOCATE ( var %     dtemp_31(mp,nd) )
   ALLOCATE ( var %     dmoist_31(mp,nd) )
   ALLOCATE ( var %     dtemp_91(mp,91) )
   ALLOCATE ( var %     alpha_PT_20(mp,ny) )
   ALLOCATE ( var %   APAR_leaf_sun(mp,ktauday*5) )
   ALLOCATE ( var %   APAR_leaf_shade(mp,ktauday*5) )
   ALLOCATE ( var %   Dleaf_sun(mp,ktauday*5) )
   ALLOCATE ( var %   Dleaf_shade(mp,ktauday*5) )
   ALLOCATE ( var %   Tleaf_sun(mp,ktauday*5) )
   ALLOCATE ( var %   Tleaf_shade(mp,ktauday*5) )
   ALLOCATE ( var %   cs_sun(mp,ktauday*5) )
   ALLOCATE ( var %   cs_shade(mp,ktauday*5) )
   ALLOCATE ( var %   scalex_sun(mp,ktauday*5) )
   ALLOCATE ( var %   scalex_shade(mp,ktauday*5) )
   ALLOCATE ( var %   fwsoil(mp,ktauday*5) )
   ALLOCATE ( var % aprecip_20(mp,ny) )
   ALLOCATE ( var %   Rd_sun(mp,ktauday*5) )
   ALLOCATE ( var %   Rd_shade(mp,ktauday*5) )

   
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
    ! Deallocate variables for SLI soil model:
    !IF(cable_user%SOIL_STRUC=='sli') THEN
    DEALLOCATE ( var % nhorizons)
    DEALLOCATE ( var % ishorizon)
    DEALLOCATE ( var % clitt )
    DEALLOCATE ( var % zeta )
    DEALLOCATE ( var % fsatmax )
    DEALLOCATE ( var % swilt_vec )
    DEALLOCATE ( var % ssat_vec )
    DEALLOCATE ( var % sfc_vec )
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

    DEALLOCATE ( var % E_fusion_sn)
    DEALLOCATE ( var % E_sublimation_sn)
    DEALLOCATE ( var % latent_heat_sn)
    DEALLOCATE ( var % evap_liq_sn)
    DEALLOCATE ( var % surface_melt)
    DEALLOCATE ( var % Qadv_rain_sn)

    ! END IF

END SUBROUTINE dealloc_soil_snow_type

! ------------------------------------------------------------------------------

SUBROUTINE dealloc_veg_parameter_type(var)

   TYPE(veg_parameter_type), INTENT(inout) :: var

   DEALLOCATE( var% canst1 )
   DEALLOCATE( var% dleaf )
   DEALLOCATE( var% ejmax )
   DEALLOCATE( var% ejmax_shade )
   DEALLOCATE( var% ejmax_sun )
   DEALLOCATE( var% iveg )
   DEALLOCATE( var% ivegp )
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
   DEALLOCATE( var% vcmax_shade )
   DEALLOCATE( var% vcmax_sun )
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
   DEALLOCATE( var% spill )
   DEALLOCATE( var% tscrn )
   DEALLOCATE( var% wcint )
   DEALLOCATE( var% tv )
   DEALLOCATE( var% us )
   DEALLOCATE( var% uscrn )
   DEALLOCATE( var% rghlai )
   DEALLOCATE( var% vlaiw )
   DEALLOCATE( var% fwet )
   DEALLOCATE( var% A_sh )
   DEALLOCATE( var% A_sl )
   DEALLOCATE( var% A_slC )
   DEALLOCATE( var% A_shC)
   DEALLOCATE( var% A_slJ )
   DEALLOCATE( var% A_shJ )
   DEALLOCATE( var% eta_A_cs )
   DEALLOCATE( var% cs )
   DEALLOCATE( var% dAdcs )
   DEALLOCATE( var% cs_sl)
   DEALLOCATE( var% cs_sh )
   DEALLOCATE( var% tlf )
   DEALLOCATE( var% dlf )

   
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
   DEALLOCATE( var% gswx )
   DEALLOCATE( var% oldcansto )
   DEALLOCATE( var% zetar )
   DEALLOCATE( var% zetash )
   DEALLOCATE ( var % fwsoil )
   DEALLOCATE ( var % ofes )

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
   DEALLOCATE ( var % rhum )
   DEALLOCATE ( var % u10 )
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
