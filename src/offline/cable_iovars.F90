!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory cTYPE(casa_flux_type), INTENT(IN) :: casaflux ! casa fluxesontaining CABLE code.
!
! ==============================================================================
! Purpose: Defines input/output related variables for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Development by Gab Abramowitz
!          Additional code to use multiple vegetation types per grid-cell (patches)
!
! ==============================================================================
MODULE cable_IO_vars_module

  USE cable_def_types_mod, ONLY : r_2, mvtype, mstype

  IMPLICIT NONE

  PUBLIC
  PRIVATE r_2, mvtype, mstype

  ! ============ Timing variables =====================
  REAL :: shod ! start time hour-of-day

  INTEGER :: sdoy,smoy,syear ! start time day-of-year month and year

  CHARACTER(LEN=200) :: timeunits ! timing info read from nc file

  CHARACTER(LEN=10) :: calendar ! 'standard' if using leap years (set by 
                                ! leaps namelist option), else 'noleap'

  CHARACTER(LEN=3) :: time_coord ! GMT or LOCal time variables

  REAL(r_2),POINTER,DIMENSION(:) :: timevar ! time variable from file

  INTEGER,DIMENSION(12) ::                                                    &
       daysm = (/31,28,31,30,31,30,31,31,30,31,30,31/),                         &
       daysml = (/31,29,31,30,31,30,31,31,30,31,30,31/),                        &
       lastday = (/31,59,90,120,151,181,212,243,273,304,334,365/),              &
       lastdayl = (/31,60,91,121,152,182,213,244,274,305,335,366/)

  LOGICAL :: leaps   ! use leap year timing?

  ! ============ Structure variables ===================
  REAL, POINTER,DIMENSION(:) :: latitude, longitude

  REAL,POINTER, DIMENSION(:,:) :: lat_all, lon_all ! lat and lon

  CHARACTER(LEN=4) :: metGrid ! Either 'land' or 'mask'

  INTEGER,POINTER,DIMENSION(:,:) :: mask ! land/sea mask from met file

  INTEGER,POINTER,DIMENSION(:) :: land_x,land_y ! indicies of land in mask

  INTEGER ::                                                                  &
       xdimsize,ydimsize,   & ! sizes of x and y dimensions
       ngridcells             ! number of gridcells in simulation

  ! For vegetated surface type
  TYPE patch_type

     REAL ::                                                                  &
          frac,    &  ! fractional cover of each veg patch
          latitude,&
          longitude

  END TYPE patch_type


  TYPE land_type

     INTEGER ::                                                               &
          nap,     & ! number of active (>0%) patches (<=max_vegpatches)
          cstart,  & ! pos of 1st gridcell veg patch in main arrays
          cend,    & ! pos of last gridcell veg patch in main arrays
          ilat,    & ! replacing land_y  ! ??
          ilon       ! replacing land_x  ! ??

  END TYPE land_type


  TYPE(land_type),DIMENSION(:),POINTER :: landpt
  TYPE(patch_type), DIMENSION(:), POINTER :: patch

  INTEGER ::                                                                  &
       max_vegpatches,   & ! The maximum # of patches in any grid cell
       nmetpatches         ! size of patch dimension in met file, if exists

  ! =============== File details ==========================
   TYPE globalMet_type
     LOGICAL           ::                                                     &
       l_gpcc,&! = .FALSE., &         ! ypwang following Chris Lu (30/oct/2012)
       l_gswp,&!= .FALSE. , &         ! BP May 2013
       l_ncar,&! = .FALSE., &         ! BP Dec 2013
       l_access ! = .FALSE.          ! BP May 2013

      CHARACTER(LEN=99) ::                                                     &
         rainf, &
         snowf, &
         LWdown, &
         SWdown, &
         PSurf, &
         Qair, &
         Tair, &
         wind

   END TYPE globalMet_type
   
   TYPE(globalMet_type) :: globalMetfile
 
  TYPE gswp_type

     CHARACTER(LEN=200) ::                                                     &
          rainf, &
          snowf, &
          LWdown, &
          SWdown, &
          PSurf, &
          Qair, &
          Tair, &
          wind, &
          mask

  END TYPE gswp_type

  TYPE(gswp_type)      :: gswpfile


  INTEGER ::                                                                  &
       ncciy,      & ! year number (& switch) for gswp run
       ncid_rin,   & ! input netcdf restart file ID
       logn          ! log file unit number

  LOGICAL ::                                                                  &
!       verbose,    & ! print init and param details of all grid cells?
!       soilparmnew   ! read IGBP new soil map. Q.Zhang @ 12/20/2010

! above replaced by below - rk4417 - phase2

       verbose! print init and param details of all grid cells? ! MMY,    & 
!   soilparmnew   ! read IGBP new soil map. Q.Zhang @ 12/20/2010
! MMY @Oct2022 change to use soilparmnew by default


  ! ================ Veg and soil type variables ============================
  INTEGER, POINTER ::                                                         &
       soiltype_metfile(:,:),  & ! user defined soil type (from met file)
       vegtype_metfile(:,:)      ! user-def veg type (from met file)

   REAL, POINTER :: vegpatch_metfile(:,:) ! Anna: patchfrac for user-def vegtype


  TYPE parID_type ! model parameter IDs in netcdf file

     INTEGER :: bch,latitude,clay,css,rhosoil,hyds,rs20,sand,sfc,silt,        &
          ssat,sucs,swilt,froot,zse,canst1,dleaf,meth,za_tq,za_uv,             &
          ejmax,frac4,hc,lai,rp20,rpcoef,shelrb, vbeta, xalbnir,               &
          vcmax,xfang,ratecp,ratecs,refsbare,isoil,iveg,albsoil,               &
          taul,refl,tauw,refw,wai,vegcf,extkn,tminvj,tmaxvj,                   &
          veg_class,soil_class,mvtype,mstype,patchfrac,                        &
          WatSat,GWWatSat,SoilMatPotSat,GWSoilMatPotSat,                       &
          HkSat,GWHkSat,FrcSand,FrcClay,Clappb,Watr,GWWatr,sfc_vec,forg,swilt_vec, &
          slope,slope_std,GWdz,SatFracmax,Qhmax,QhmaxEfold,HKefold,HKdepth,   &
          sand_vec,clay_vec,bch_vec,org_vec,elev,elev_std  ! inserted by rk4417 - phase2

     INTEGER :: ishorizon,nhorizons,clitt, &
          zeta,fsatmax, &
          gamma,ZR,F10

     INTEGER :: g0,g1 ! Ticket #56

  END TYPE parID_type

  ! =============== Logical  variables ============================
  TYPE input_details_type

     LOGICAL ::                                                               &
          Wind,    & ! T => 'Wind' is present; F => use vector component wind
          LWdown,  & ! T=> downward longwave is present in met file
          CO2air,  & ! T=> air CO2 concentration is present in met file
          PSurf,   & ! T=> surface air pressure is present in met file
          Snowf,   & ! T=> snowfall variable is present in met file
          avPrecip,& ! T=> ave rainfall present in met file (use for spinup)
          LAI,     & ! T=> LAI is present in the met file
          LAI_T,   & ! T=> LAI is time dependent, for each time step
          LAI_M,   & ! T=> LAI is time dependent, for each month
          LAI_P,   & ! T=> LAI is patch dependent
          parameters,&! TRUE if non-default parameters are found
          initial, & ! switched to TRUE when initialisation data are loaded
          patch,   & ! T=> met file have a subgrid veg/soil patch dimension
          laiPatch   ! T=> LAI file have a subgrid veg patch dimension

  END TYPE input_details_type

  TYPE(input_details_type) :: exists

  TYPE output_inclusion_type

     ! Which variables to include in output file, values initialised here
     ! and can be reset by namelist file read in driver:
     ! Groups of output variables:

     LOGICAL ::                                                               &
          met = .FALSE.,       & ! input met data
          flux = .FALSE.,      &  ! convective, runoff, NEE
          radiation = .FALSE., & ! net rad, albedo
          carbon = .FALSE.,    & ! NEE, GPP, NPP, stores
          soil = .FALSE.,      &  ! soil states
          snow = .FALSE.,      &  ! snow states
          veg = .FALSE.,       & ! vegetation states
          params = .FALSE.,    & ! input parameters used to produce run
          balances = .FALSE.,  & ! energy and water balances
          restart = .FALSE.,   & ! create restart file?
          ensemble = .FALSE.,  & ! are we creating an ensemble run?
          patch = .FALSE. , &   ! should patch-specific info be written
                                ! to output file?
                                ! vh_js !
          casa = .FALSE.       ! additional casa outputs (C stores and plant turnover)

     ! Should output grid follow met file 'default'; force with 'land' or 'mask':
     CHARACTER(LEN=7) ::                                                      &
          grid = 'default', &
          averaging = 'all' ! 'all', 'daily', 'monthly', 'user6'(6hrly)

     INTEGER ::                                                               &
          interval ! in case of 'user6' above, interval will be 6

     ! variables specified individually:
     LOGICAL ::                                                               &
          SWdown = .FALSE.,    & ! 6 downward short-wave radiation [W/m2]
          LWdown = .FALSE.,    & ! 7 downward long-wave radiation [W/m2]
          Rainf = .FALSE.,     & ! 8 rainfall [kg/m2/s]
          Snowf = .FALSE.,     & ! 9 snowfall [kg/m2/s]
          PSurf = .FALSE.,     & ! 10 surface pressure [Pa]
          Tair = .FALSE.,      & ! 11 surface air temperature [K]
          Qair = .FALSE.,      & ! 12 specific humidity [kg/kg]
          Tscrn = .FALSE.,     & !    screen level air temperature [oC]
          Tex = .FALSE.,       & !    extremes in screen level temperature [oC]
          Qscrn = .FALSE.,     & !    screen level specific humidity [kg/kg]
          CO2air = .FALSE.,    & ! 13 CO2 concentration [ppmv]
          Wind = .FALSE.,      & ! 14 windspeed [m/s]
          Wind_N = .FALSE.,    & ! 15 surface wind speed, N component [m/s]
          Wind_E = .FALSE.,    & ! 16 surface wind speed, E component [m/s]
          LAI = .FALSE.,       & !
          Qmom = .FALSE.,      & !    momentum flux [kg/m/s2]
          Qh = .FALSE.,        & ! 17 sensible heat flux [W/m2]
          Qle = .FALSE.,       & ! 18 latent heat flux [W/m2]
          Qg = .FALSE.,        & ! 19 ground heat flux [W/m2]
          SWnet = .FALSE.,     & ! 20 net shortwave [W/m2]
          LWnet = .FALSE.,     & ! 21 net longwave [W/m2]
          Evap = .FALSE.,      & ! 22 total evapotranspiration [kg/m2/s]
          Ewater = .FALSE.,    & ! 23 evap. from surface water storage [kg/m2/s]
          ESoil = .FALSE.,     & ! 24 bare soil evaporation [kg/m2/s]
          TVeg = .FALSE.,      & ! 25 vegetation transpiration [kg/m2/s]
          ECanop = .FALSE.,    & ! 26 interception evaporation [kg/m2/s]
          PotEvap = .FALSE.,   & ! 27 potential evapotranspiration [kg/m2/s]
          ACond = .FALSE.,     & ! 28 aerodynamic conductance [m/s]
          SoilWet = .FALSE.,   & ! 29 total soil wetness [-]
          Albedo = .FALSE.,    & ! 30 albedo [-]
          visAlbedo = .FALSE., & ! vars intro for Ticket #27
          nirAlbedo = .FALSE., & ! vars intro for Ticket #27
          VegT = .FALSE.,      & ! 31 vegetation temperature [K]
          SoilTemp = .FALSE.,  & ! 32 av.layer soil temperature [K]
          SoilMoist = .FALSE., & ! 33 av.layer soil moisture [kg/m2]
          SoilMoistIce = .FALSE., & ! 33 av.layer soil frozen moisture [kg/m2]
          Qs = .FALSE.,        & ! 34 surface runoff [kg/m2/s]
          Qsb = .FALSE.,       &! 35 subsurface runoff [kg/m2/s]
          DelSoilMoist = .FALSE., & ! 36 change in soilmoisture
                                ! (sum layers) [kg/m2]
          DelSWE = .FALSE.,    & ! 37 change in snow water equivalent [kg/m2]
          DelIntercept = .FALSE.,& ! 38 change in interception storage [kg/m2]
          SnowT = .FALSE.,     & ! 39 snow surface temp [K]
          BaresoilT = .FALSE., & ! 40 surface bare soil temp [K]
          AvgSurfT = .FALSE.,  & ! 41 Average surface temperature [K]
          RadT = .FALSE.,      & ! 42 Radiative surface temperature [K]
          SWE = .FALSE.,       & ! 43 snow water equivalent [kg/m2]
          SnowMelt = .FALSE.,       & ! 43 snow melt [kg/m2/s] !vh!
          RootMoist = .FALSE., & ! 44 root zone soil moisture [kg/m2]
          CanopInt = .FALSE.,  & ! 45 total canopy water storage [kg/m2]
          NEE  = .FALSE.,      & ! 46 net ecosystem exchange [umol/m2/s]
          NPP  = .FALSE.,      & ! 47 net primary production of C
                                ! by veg [umol/m2/s]
          GPP = .FALSE.,       & ! 48 gross primary production C
                                ! by veg [umol/m2/s]
          AutoResp = .FALSE.,  & ! 49 autotrophic respiration [umol/m2/s]
          LeafResp = .FALSE.,  & ! 51 autotrophic respiration [umol/m2/s]
          HeteroResp = .FALSE.,& ! 50 heterotrophic respiration [umol/m2/s]
          SnowDepth = .FALSE., & ! actual depth of snow in [m]
                                !variables
          Rnet = .FALSE.,      & ! net absorbed radiation [W/m2]
          HVeg = .FALSE.,      & ! sensible heat from vegetation [W/m2]
          HSoil = .FALSE.,     & ! sensible heat from soil [W/m2]
          RnetSoil = .FALSE.,     & ! sensible heat from soil [W/m2] !vh!
          Ebal = .FALSE.,      & ! cumulative energy balance [W/m2]
          Wbal = .FALSE.,      & ! cumulative water balance [W/m2]
                                ! vh_js ! added CanT and fwsoil to the list
          CanT = .FALSE.,      & ! within-canopy temperature [K]
          Fwsoil = .FALSE.,      & ! soil moisture modifier to stomatal conductance
          Area = .FALSE., & ! patch area in km2
                                !mrd561
                                !MD GW
          GWMoist = .FALSE.,   & ! water balance of aquifer [mm3/mm3]
          WatTable = .FALSE.,  & ! water table depth [m]
          Qrecharge=.FALSE.,   &  !recharge to /from auqifer
          SatFrac=.FALSE.,       & ! Saturated Fraction of Gridcell (tile)

                                ! vh_js ! additional casa variables
          NBP = .FALSE., &
          dCdt = .FALSE., &
          TotSoilCarb = .FALSE.,   &
          TotLivBiomass = .FALSE., &
          TotLittCarb = .FALSE., &
          SoilCarbFast = .FALSE., &
          SoilCarbSlow = .FALSE., &
          SoilCarbPassive = .FALSE., &
          LittCarbMetabolic = .FALSE., &
          LittCarbStructural = .FALSE., &
          LittCarbCWD = .FALSE., &
          PlantCarbLeaf = .FALSE., &
          PlantCarbFineRoot = .FALSE., &
          PlantCarbWood = .FALSE., &
          PlantTurnover = .FALSE., &
          PlantTurnoverLeaf = .FALSE., &
          PlantTurnoverFineRoot = .FALSE., &
          PlantTurnoverWood = .FALSE., &
          PlantTurnoverWoodDist = .FALSE., &
          PlantTurnoverWoodCrowding = .FALSE., &
          PlantTurnoverWoodResourceLim = .FALSE., &
          LandUseFlux = .FALSE., &
                                !parameters
          bch = .FALSE.,       & ! parameter b in Campbell equation 1985
          latitude = .FALSE.,  & ! site latitude
          clay = .FALSE.,      & ! fraction of clay in soil
          css = .FALSE.,       & ! heat capacity of soil minerals [J/kg/C]
          rhosoil = .FALSE.,   & ! soil density [kg/m3]
          hyds = .FALSE.,      & ! hydraulic conductivity @ saturation [m/s], Ksat
          rs20 = .FALSE.,      & ! soil respiration at 20 C [dimensionless],
                                ! (0.1 - 10), prop to om
          sand  = .FALSE.,     & ! fraction of sand in soil
          sfc = .FALSE.,       & ! vol H2O @ field capacity
          silt  = .FALSE.,     & ! fraction of silt in soil
          ssat = .FALSE.,      & ! vol H2O @ saturation
          sucs = .FALSE.,      & ! suction at saturation [m]
          swilt = .FALSE.,     & ! vol H2O @ wilting
          froot = .FALSE.,     & ! fraction of roots in each soil layer
          zse = .FALSE.,       & ! thickness of each soil layer (1=top) (m)
          canst1 = .FALSE.,    & ! max intercepted water by canopy [mm/LAI]
                                ! (0.08 - 0.12) {avoid}
          dleaf = .FALSE.,     & ! chararacteristic length of leaf [m],
                                ! (0.005 - 0.2) pine -> tropical
          ejmax  = .FALSE.,    & ! max pot. electron transport rate
                                ! top leaf[mol/m2/s](1e-5 - 3e-4) {use}
          frac4  = .FALSE.,    & ! fraction of c4 plants [-]
          hc = .FALSE.,        & ! height of canopy [m]
          rp20  = .FALSE.,     & ! plant respiration coefficient at
                                ! 20 C [-] 0.1 - 10 (frp 0 - 15e-6 mol/m2/s)
          g0   = .FALSE.,      & ! Ticket #56
          g1   = .FALSE.,      & ! Ticket #56
          rpcoef  = .FALSE.,   & ! temperature coef nonleaf plant
                                ! respiration [1/C] (0.8 - 1.5)
          shelrb  = .FALSE.,   & ! sheltering factor [-] {avoid - insensitive?}
          vcmax  = .FALSE.,    & ! maximum RuBP carboxylation rate
                                ! top leaf [mol/m2/s](5e-6 - 1.5e-4){use}
          xfang  = .FALSE.,    & ! leaf angle PARAMETER (dimensionless)
                                ! (v leaf -1.0 horiz 1.0 sphere 0 (-1 - 1))
          wai    = .FALSE.,    & ! wood area index
          vegcf  = .FALSE.,    & !
          extkn  = .FALSE.,    & !
          ratecp = .FALSE.,    & ! plant carbon pool rate constant (1/year)
          ratecs = .FALSE.,    & ! soil carbon pool rate constant (1/year)
          albsoil = .FALSE.,   & ! soil reflectance [-]
          taul = .FALSE.,      & ! leaf transmissivity [-](V:0.07 - 0.15
                                ! NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
          refl = .FALSE.,      & ! leaf reflectance [-](V:0.07 - 0.15 \
                                ! NIR: 0.3 - 0.6 IR: 0.0 - 0.05)
          tminvj = .FALSE.,    & ! min temperature of the start of
                                ! photosynthesis(leaf phenology)[-] (-10 - 10)
          tmaxvj  = .FALSE.,   & ! max temperature of the start of
                                ! photosynthesis(leaf phenology)[-] (-5 - 15)
          vbeta = .FALSE.,     & ! stomatal sensitivity to soil water
          xalbnir = .FALSE.,   & ! modifier for albedo in near ir band
          iveg  = .FALSE.,     & ! vegetation type from global index
          patchfrac  = .FALSE.,& ! fractional cover of each veg/soil patch
          isoil  = .FALSE.,    & ! soil type from global index
          meth  = .FALSE.,     & ! method for solving turbulence in canopy scheme
          za  = .FALSE.,       & ! something to do with roughness ????
          elev = .false.,&       !mean subgrid elev    ! inserted 2 lines - rk4417 - phase2
          elev_std=.false.,&     !stddev of subgrid elev
          slope = .FALSE.,&      !mean subgrid slope
          slope_std=.FALSE.,&    !stddev of subgrid slope
          GWdz=.FALSE.,&         !aquifer thickness
          SatFracmax=.FALSE.,&
          Qhmax=.FALSE.,&
          QhmaxEfold=.FALSE.,&
          HKefold=.FALSE.,&
!          HKdepth ! commented out by rk4417 - phase2
          HKdepth=.false.,&  ! added this block - rk4417 - phase2
          SMP=.false.,&
          SMP_hys=.false.,&
          WB_hys=.false.,&
          SSAT_hys=.false.,&
          WATR_hys=.false.,&
          hys_fac=.false.

  END TYPE output_inclusion_type


  TYPE(output_inclusion_type),SAVE :: output
  TYPE(output_inclusion_type),SAVE :: patchout ! do we want patch-specific info

  ENUM, BIND(C)
    ENUMERATOR :: NO_CHECK = 0
    ENUMERATOR :: ON_TIMESTEP = 1
    ENUMERATOR :: ON_WRITE = 2
    ENUMERATOR :: RANGE_CHECK
  END ENUM
  TYPE checks_type
    LOGICAL :: energy_bal, mass_bal
    INTEGER(KIND(RANGE_CHECK)) :: ranges  ! 0 = NO , 1 = TIMESTEP , 2 = WRITE
    LOGICAL :: exit
  END TYPE checks_type

  TYPE(checks_type) :: check ! what types of checks to perform

  ! ============== Proxy input variables ================================
  REAL,POINTER,DIMENSION(:)  :: PrecipScale! precip scaling per site for spinup
  REAL,POINTER,DIMENSION(:,:)  :: defaultLAI ! in case met file/host model
  ! has no LAI
  REAL :: fixedCO2 ! CO2 level if CO2air not in met file

  ! For threading:
  !$OMP THREADPRIVATE(landpt,patch)
CONTAINS
  SUBROUTINE set_group_output_values

    !*#Purpose:
    ! Set individual variables to output according to the values of the group options from the namelist entries in `output%`.
    IF (output%params) THEN
      output%iveg = .TRUE.
      output%patchfrac = .TRUE.
      output%isoil = .TRUE.
      output%bch = .TRUE.
      output%clay = .TRUE.
      output%sand = .TRUE.
      output%silt = .TRUE.
      output%css = .TRUE.
      output%rhosoil = .TRUE.
      output%hyds = .TRUE.
      output%sucs = .TRUE.
      output%rs20 = .TRUE.
      output%ssat = .TRUE.
      output%sfc = .TRUE.
      output%swilt = .TRUE.
      output%albsoil = .TRUE.
      output%canst1 = .TRUE.
      output%dleaf = .TRUE.
      output%ejmax = .TRUE.
      output%vcmax = .TRUE.
      output%frac4 = .TRUE.
      output%hc = .TRUE.
      output%rp20 = .TRUE.
      output%g0 = .TRUE.
      output%g1 = .TRUE.
      output%rpcoef = .TRUE.
      output%shelrb = .TRUE.
      output%xfang = .TRUE.
      output%wai = .TRUE.
      output%vegcf = .TRUE.
      output%extkn = .TRUE.
      output%tminvj = .TRUE.
      output%tmaxvj = .TRUE.
      output%vbeta = .TRUE.
      output%xalbnir = .TRUE.
      output%meth = .TRUE.
      output%za = .TRUE.
      output%ratecp = .TRUE.
      output%ratecs = .TRUE.
      output%froot = .TRUE.
      output%zse = .TRUE.
      output%slope = .TRUE.
      output%slope_std = .TRUE.
      output%GWdz = .TRUE.
    END IF

    IF (output%met) THEN
      output%Swdown = .TRUE.
      output%Lwdown = .TRUE.
      output%Rainf = .TRUE.
      output%Snowf = .TRUE.
      output%PSurf = .TRUE.
      output%Tair = .TRUE.
      output%Qair = .TRUE.
      output%Wind = .TRUE.
      output%CO2air = .TRUE.
    END IF

    IF (output%flux) THEN
      output%Qmom = .TRUE.
      output%Qh = .TRUE.
      output%Qle = .TRUE.
      output%Qg = .TRUE.
      output%Qs = .TRUE.
      output%Qsb = .TRUE.
      output%Evap = .TRUE.
      output%ECanop = .TRUE.
      output%PotEvap = .TRUE.
      output%TVeg = .TRUE.
      output%ESoil = .TRUE.
      output%HVeg = .TRUE.
      output%HSoil = .TRUE.
      output%RNetSoil = .TRUE.
      output%NEE = .TRUE.
    END IF

    IF (output%soil) THEN
      output%SoilMoist = .TRUE.
      output%SoilTemp = .TRUE.
      output%BaresoilT = .TRUE.
      output%WatTable = .TRUE.
      output%GWMoist = .TRUE.
      output%SatFrac = .TRUE.
      output%Qrecharge = .TRUE.
    END IF

    IF (output%snow) THEN
      output%SWE = .TRUE.
      output%SnowT = .TRUE.
      output%SnowDepth = .TRUE.
    END IF

    IF (output%radiation) THEN
      output%Swnet = .TRUE.
      output%Lwnet = .TRUE.
      output%Rnet = .TRUE.
      output%Albedo = .TRUE.
      output%RadT = .TRUE.
    END IF

    IF (output%veg) THEN
      output%Tscrn = .TRUE.
      output%Tex = .TRUE.
      output%Qscrn = .TRUE.
      output%VegT = .TRUE.
      output%CanT = .TRUE.
      output%Fwsoil = .TRUE.
      output%CanopInt = .TRUE.
      output%LAI = .TRUE.
    END IF

    IF (output%balances) THEN
      output%Ebal = .TRUE.
      output%Wbal = .TRUE.
    END IF

    IF (output%carbon) THEN
      output%GPP = .TRUE.
      output%NPP = .TRUE.
      output%NBP = .TRUE.
      output%NEE = .TRUE.
      output%AutoResp = .TRUE.
      output%LeafResp = .TRUE.
      output%HeteroResp = .TRUE.
    END IF

  END SUBROUTINE set_group_output_values

END MODULE cable_IO_vars_module
