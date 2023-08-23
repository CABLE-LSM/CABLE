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
! Purpose: Output module for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Developed by Gab Abramowitz
!          Output of additional variables and parameters relative to v1.4b
!
!
! ==============================================================================
! CALLed from:    cable_driver.F90
!
! MODULEs used:   cable_abort_module
!                 cable_common_module
!                 cable_checks_module
!                 cable_def_types_mod
!                 cable_IO_vars_module
!                 cable_write_module
!                 netcdf
!
! CALLs:          open_output_file
!                 write_output
!                 close_output_file
!                 create_restart
!
MODULE cable_output_module

  USE cable_abort_module,  ONLY: cable_abort, nc_abort
  USE cable_def_types_mod
  USE casavariable,        ONLY: casa_pool, casa_flux, casa_met
  USE cable_IO_vars_module
  USE cable_checks_module, ONLY: mass_balance, energy_balance, ranges
  USE cable_write_module
  USE netcdf
  USE cable_common_module, ONLY: filename, CurYear,IS_LEAPYEAR, cable_user ! , calcsoilalbedo
  ! 13C
  USE cable_c13o2_def,     only: c13o2_pool,  c13o2_flux

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: open_output_file, write_output, close_output_file, create_restart

  INTEGER :: ncid_out ! output data netcdf file ID
  TYPE out_varID_type ! output variable IDs in netcdf file
     INTEGER :: SWdown, FracDiff, LWdown, Wind, Wind_E, PSurf, &
          Tair, Qair, Rainf, Snowf, CO2air, &
          Qle, Qh, Qg, NEE, fbeam, SWnet, &
          LWnet, SoilMoist, SoilMoistPFT, SoilTemp, Albedo, &
          visAlbedo, nirAlbedo, SoilMoistIce, &
          Qs, Qsb, Evap, BaresoilT, SWE, SnowT, &
          RadT, VegT, Ebal, Wbal, AutoResp, RootResp, &
          StemResp,LeafResp, HeteroResp, GPP, NPP, LAI, &
          ECanop, TVeg, ESoil, CanopInt, SnowDepth, &
          HVeg, HSoil, Rnet, tvar, CanT,Fwsoil, RnetSoil, SnowMelt, &
                                ! vh_mc ! additional variables for ESM-SnowMIP
          hfds, hfdsn, hfls, hfmlt, hfrs, hfsbl, hfss, rlus, rsus, &
          esn, evspsbl, evspsblsoi, evspsblveg, mrrob, mrros, sbl, &
          snm, snmsl, tran, albs, albsn, cw, lqsn, lwsnl, mrfsofr, &
          mrlqso, mrlsl, snc, snd, snw, snwc, tcs, tgs, ts, tsl, &
          tsn, tsns, &
          NBP, TotSoilCarb, TotLivBiomass, &
          TotLittCarb, SoilCarbFast, SoilCarbSlow, SoilCarbPassive, &
          LittCarbMetabolic, LittCarbStructural, LittCarbCWD, &
          PlantCarbLeaf, PlantCarbFineRoot, PlantCarbWood, &
          PlantTurnover, PlantTurnoverLeaf, PlantTurnoverFineRoot, &
          PlantTurnoverWood, PlantTurnoverWoodDist, PlantTurnoverWoodCrowding, &
          PlantTurnoverWoodResourceLim, dCdt, Area, LandUseFlux, patchfrac, &
          vcmax, ejmax, hc, &
          GPP_sl, GPP_sh, GPP_slC, GPP_shC, GPP_slJ, GPP_shJ, &
          eta_GPP_cs,  eta_TVeg_cs, dGPPdcs, CO2s, gsw_sl, gsw_sh, gsw_TVeg, &
          An_sl, An_sh, ci_sl, ci_sh, scalex_sl, scalex_sh, dlf, &
          vcmax_ts, jmax_ts, qcan_sl, qcan_sh, &
          An, Rd, cplant, clitter, csoil, clabile, &
          A13n, aDisc13, c13plant, c13litter, c13soil, c13labile
  END TYPE out_varID_type
  TYPE(out_varID_type) :: ovid ! netcdf variable IDs for output variables

  TYPE(parID_type) :: opid ! netcdf variable IDs for output variables

  TYPE output_temporary_type
     REAL(KIND=4), POINTER, DIMENSION(:) :: SWdown => null()   ! 5 downward short-wave radiation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: FracDiff => null() ! 6 Fraction of diffuse radiation
     REAL(KIND=4), POINTER, DIMENSION(:) :: LWdown => null()   ! 7 downward long-wave radiation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Rainf => null()    ! 8 rainfall [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Snowf => null()    ! 9 snowfall [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: PSurf => null()    ! 10 surface pressure [Pa]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tair => null()     ! 11 surface air temperature
     ! [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qair => null()   ! 12 specific humidity [kg/kg]
     REAL(KIND=4), POINTER, DIMENSION(:) :: CO2air => null() ! 13 CO2 concentration [ppmv]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wind => null()   ! 14 windspeed [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wind_N => null() ! 15 surface wind speed, N component [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wind_E => null() ! 16 surface wind speed, E component [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LAI => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qh => null()     ! 17 sensible heat flux [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qle => null()    ! 18 latent heat flux [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qg => null()     ! 19 ground heat flux [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: fbeam => null()  ! 20 fracion of direct radiation
     REAL(KIND=4), POINTER, DIMENSION(:) :: SWnet => null()  ! 20 net shortwave [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LWnet => null()  ! 21 net longwave [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Evap => null()   ! 22 total evapotranspiration [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Ewater => null() ! 23 evap. from surface water storage [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: ESoil => null()  ! 24 bare soil evaporation [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: TVeg => null()   ! 25 vegetation transpiration [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: ECanop => null() ! 26 interception evaporation [kg/m2/s]
     ! 27 potential evapotranspiration [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: PotEvap => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: ACond => null()   ! 28 aerodynamic conductance [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilWet => null() ! 29 total soil wetness [-]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Albedo => null()  ! 30 albedo [-]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: visAlbedo => null()  ! vars intro for Ticket #27
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: nirAlbedo => null()  ! vars intro for Ticket #27
     REAL(KIND=4), POINTER, DIMENSION(:) :: VegT => null()    ! 31 vegetation temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilTemp => null()  ! 32 av.layer soil temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilMoist => null() ! 33 av.layer soil moisture [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:)   :: SoilMoistPFT => null() ! 33 soil moisture per PFT [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilMoistIce => null() ! 33 av.layer soil frozen moisture [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qs => null()  ! 34 surface runoff [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qsb => null() ! 35 subsurface runoff [kg/m2/s]
     ! 36 change in soilmoisture (sum layers) [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: DelSoilMoist => null()
     ! 37 change in snow water equivalent [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: DelSWE => null()
     ! 38 change in interception storage [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: DelIntercept => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: SnowT => null()     ! 39 snow surface temp [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: BaresoilT => null() ! 40 surface bare soil temp [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: AvgSurfT => null()  ! 41 Average surface temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: RadT => null()      ! 42 Radiative surface temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SWE => null()       ! 43 snow water equivalent [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: RootMoist => null() ! 44 root zone soil moisture [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: CanopInt => null()  ! 45 total canopy water storage [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: NEE => null()       ! 46 net ecosystem exchange [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: NPP => null()       ! 47 net primary production of C by veg [umol/m2/s]
     ! 48 gross primary production C by veg [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: AutoResp => null()   ! 49 autotrophic respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LeafResp => null()   ! 51 autotrophic respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: HeteroResp => null() ! 50 heterotrophic respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SnowDepth => null()  ! actual depth of snow in [m]
     ! vh_mc ! additional variables for ESM-SnowMIP
     real(kind=4), pointer, dimension(:)   :: hfds => null()  ! downward heat flux at ground surface [W/m2]
     real(kind=4), pointer, dimension(:)   :: hfdsn => null() ! downward heat flux into snowpack [W/m2]
     real(kind=4), pointer, dimension(:)   :: hfls => null()  ! surface upward latent heat flux [W/m2]
     real(kind=4), pointer, dimension(:)   :: hfmlt => null() ! energy of fusion [W/m2]
     real(kind=4), pointer, dimension(:)   :: hfrs => null()  ! heat transferred to snowpack by rain [W/m2]
     real(kind=4), pointer, dimension(:)   :: hfsbl => null() ! energy of sublimation [W/m2]
     real(kind=4), pointer, dimension(:)   :: hfss => null()  ! surface upward sensible heat flux [W/m2]
     real(kind=4), pointer, dimension(:)   :: rlus => null()  ! surface upwelling longwave radiation [W/m2]
     real(kind=4), pointer, dimension(:)   :: rsus => null()  ! surface upwelling shortwave radiation [W/m2]
     real(kind=4), pointer, dimension(:)   :: esn => null()   ! liquid water evaporation from snowpack [kg/m2/s]
     ! total water vapour flux from the surface to the atmosphere [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: evspsbl => null()
     real(kind=4), pointer, dimension(:)   :: evspsblsoi => null() ! evaporation and sublimation from soil [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: evspsblveg => null() ! evaporation and sublimation from canopy [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: mrrob => null()   ! subsurface runoff [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: mrros => null()   ! surface runoff [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: sbl => null()     ! sublimation of snow [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: snm => null()     ! surface snow melt [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: snmsl => null()   ! water flowing out of snowpack [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: tran => null()    ! transpiration [kg/m2/s]
     real(kind=4), pointer, dimension(:)   :: albs => null()    ! surface albedo [-]
     real(kind=4), pointer, dimension(:)   :: albsn => null()   ! snow albedo [-]
     real(kind=4), pointer, dimension(:)   :: cw => null()      ! total canopy water storage [kg/m2]
     real(kind=4), pointer, dimension(:)   :: lqsn => null()    ! mass fraction of liquid water in snowpack [-]
     real(kind=4), pointer, dimension(:)   :: lwsnl => null()   ! liquid water content of snowpack [kg/m2]
     real(kind=4), pointer, dimension(:,:) :: mrfsofr => null() ! mass fractions of frozen water in soil layers [-]
     real(kind=4), pointer, dimension(:,:) :: mrlqso => null()  ! mass fractions of unfrozen water in soil layers [-]
     real(kind=4), pointer, dimension(:,:) :: mrlsl => null()   ! masses of frozen and unfrozen moisture in soil layers [kg/m2]
     real(kind=4), pointer, dimension(:)   :: snc => null()     ! snow area fraction [-]
     real(kind=4), pointer, dimension(:)   :: snd => null()     ! snowdepth [m]
     real(kind=4), pointer, dimension(:)   :: snw => null()     ! mass of snowpack [kg/m2]
     real(kind=4), pointer, dimension(:)   :: snwc => null()    ! mass of snow intercepted by vegetation [kg/m2]
     real(kind=4), pointer, dimension(:)   :: tcs => null()     ! vegetation canopy temperature [K]
     real(kind=4), pointer, dimension(:)   :: tgs => null()     ! temperature of bare soil [K]
     real(kind=4), pointer, dimension(:)   :: ts => null()      ! surface temperature [K]
     real(kind=4), pointer, dimension(:,:) :: tsl => null()     ! temperatures of soil layers [K]
     real(kind=4), pointer, dimension(:)   :: tsn => null()     ! snow internal temperature [K]
     real(kind=4), pointer, dimension(:)   :: tsns => null()    ! snow surface temperature [K]
     ! Non-Alma variables
     REAL(KIND=4), POINTER, DIMENSION(:) :: Rnet => null()  ! net absorbed radiation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: HVeg => null()  ! sensible heat from vegetation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: HSoil => null() ! sensible heat from soil [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: RnetSoil => null() ! latent heat from soil [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SnowMelt => null() ! snow melt [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Ebal => null()  ! cumulative energy balance [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wbal => null()  ! cumulative water balance [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: CanT => null()  ! within-canopy temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Fwsoil => null()  ! soil-moisture modfier to stomatal conductance [-]

     ! [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: NBP => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: dCdt => null()
     ! [kg C /m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: TotSoilCarb => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: TotLivBiomass => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: TotLittCarb => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilCarbFast => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilCarbSlow => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilCarbPassive => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: LittCarbMetabolic => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: LittCarbStructural => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: LittCarbCWD => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantCarbLeaf => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantCarbFineRoot => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantCarbWood => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnover => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverLeaf => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverFineRoot => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWood => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWoodDist => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWoodCrowding => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWoodResourceLim => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: Area => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: LandUseFlux => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: vcmax => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: ejmax => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: patchfrac => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: hc => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP_sl => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP_sh => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP_slC => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP_shC => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP_slJ => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP_shJ => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: An_sl => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: An_sh=> null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: ci_sl => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: ci_sh => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: scalex_sl => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: scalex_sh => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: dlf => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: eta_GPP_cs => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: eta_TVeg_cs => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: dGPPdcs => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: CO2s => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: gsw_TVeg => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: vcmax_ts => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: jmax_ts => null()
     REAL(KIND=4), POINTER, DIMENSION(:) :: gsw_sl => null()   ! stomatal conductance (sunlit leaves)
     REAL(KIND=4), POINTER, DIMENSION(:) :: gsw_sh => null()   ! stomatal conductance (shaded leaves)
     REAL(KIND=4), POINTER, DIMENSION(:) :: RootResp => null()   ! autotrophic root respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: StemResp => null()   ! autotrophic stem respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: qcan_sl => null()   ! absorbed radiation by canopy (sunlit leaves) [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: qcan_sh => null()   ! absorbed radiation by canopy (shaded leaves) [W/m2]
     REAL(KIND=r_2), POINTER, DIMENSION(:)   :: An => null()        ! total net assimilation
     REAL(KIND=r_2), POINTER, DIMENSION(:)   :: Rd => null()        ! total leaf respiration
     REAL(KIND=r_2), POINTER, DIMENSION(:,:) :: cplant => null()    ! plant carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:,:) :: clitter => null()   ! litter carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:,:) :: csoil => null()     ! soil carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:)   :: clabile => null()   ! excess carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:)   :: A13n => null()      ! total 13C net assimilation
     REAL(KIND=r_2), POINTER, DIMENSION(:)   :: aDisc13 => null()   ! 13C discrimination times gross assimilation
     REAL(KIND=r_2), POINTER, DIMENSION(:,:) :: c13plant => null()  ! 13C plant carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:,:) :: c13litter => null() ! 13C litter carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:,:) :: c13soil => null()   ! 13C soil carbon pools
     REAL(KIND=r_2), POINTER, DIMENSION(:)   :: c13labile => null() ! 13C excess carbon pools
  END TYPE output_temporary_type
  TYPE(output_temporary_type), SAVE :: out

  INTEGER :: ok   ! netcdf error status

  interface toreal4
     module procedure dp2sp, i2sp, sp2sp
  end interface toreal4

CONTAINS

  SUBROUTINE open_output_file(dels, soil, veg, bgc, rough)

    use casadimension, only: mplant, mlitter, msoil

    ! Creates netcdf output file, defines all variables
    ! and writes parameters to it if requested by user.
    REAL,                      INTENT(IN) :: dels ! time step size
    TYPE(soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE(veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
    TYPE(bgc_pool_type),       INTENT(IN) :: bgc
    TYPE(roughness_type),      INTENT(IN) :: rough
    ! REAL, POINTER,DIMENSION(:,:) :: surffrac ! fraction of each surf type

    INTEGER :: xID, yID, zID, radID, soilID, soilcarbID, &
         plantcarbID, tID, landID, patchID ! dimension IDs
    INTEGER :: latID, lonID, llatvID, llonvID ! time,lat,lon variable ID
    INTEGER :: xvID, yvID   ! coordinate variable IDs for GrADS readability
    !    INTEGER :: surffracID         ! surface fraction varaible ID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
    integer :: nplantid, nlitterid, nsoilid

    real(kind=4), parameter :: zero4 = real(0.0,4)

    ! Create output file:
#ifdef __NETCDF3__
    ok = NF90_CREATE(trim(filename%out), ior(nf90_clobber,nf90_64bit_offset), ncid_out)
#else
    ok = NF90_CREATE(trim(filename%out), ior(nf90_clobber,ior(nf90_netcdf4,nf90_classic_model)), ncid_out)
#endif
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating output file ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_out, 'x', xdimsize, xID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining x dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'y', ydimsize, yID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining y dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Define patch dimension, whether it's used or not:
    ok = NF90_DEF_DIM(ncid_out, 'patch', max_vegpatches, patchID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error defining patch dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    !   ! Define surftype dimension (currently only used for surffrac variable):
    !    ok = NF90_DEF_DIM(ncid_out,'surftype',4,surftypeID)
    !    IF (ok /= NF90_NOERR) CALL nc_abort &
    !         (ok,'Error defining syrftype dimension in output file. '// &
    !         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining vertical soil dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'rad', nrb, radID) ! number of radiation bands
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining radiation dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'soil_carbon_pools', ncs, soilcarbID) ! # pools
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining soil carbon pool dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out,'plant_carbon_pools',ncp,plantcarbID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error defining plant carbon pool dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'time', NF90_UNLIMITED, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error defining time dimension in output file. '// &
         '(SUBROUTINE open_output_file)')

    ! 13C
    if (cable_user%c13o2) then
       ok = NF90_DEF_DIM(ncid_out,'nplant',mplant,nplantid)
       if (ok /= NF90_NOERR) &
            call nc_abort(ok,'Error defining nplant dimension in output file. (SUBROUTINE open_output_file)')
       ok = NF90_DEF_DIM(ncid_out,'nlitter',mlitter,nlitterid)
       if (ok /= NF90_NOERR) &
            call nc_abort(ok,'Error defining nlitter dimension in output file. (SUBROUTINE open_output_file)')
       ok = NF90_DEF_DIM(ncid_out,'nsoil',msoil,nsoilid)
       if (ok /= NF90_NOERR) &
            call nc_abort(ok,'Error defining nsoil dimension in output file. (SUBROUTINE open_output_file)')
    endif

    IF(output%grid == 'mask' .OR. output%grid == 'ALMA' .OR. &
         (metGrid == 'mask' .AND. output%grid == 'default')) THEN
       ! for land/sea mask type grid:
       ! Atmospheric 'z' dim of size 1 to comply with ALMA grid type:
       ok = NF90_DEF_DIM(ncid_out, 'z', 1, zID)
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error defining z dimension in output file. '// &
            '(SUBROUTINE open_output_file)')
    ELSE IF(output%grid == 'land' .OR. &
         (metGrid == 'land' .AND. output%grid == 'default')) THEN
       ! For land only compression grid:
       ok = NF90_DEF_DIM(ncid_out, 'land', mland, landID) ! number of land
       ! points
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error defining land dimension in output file. '// &
            '(SUBROUTINE open_output_file)')

       ok = NF90_DEF_VAR(ncid_out, 'local_lat', NF90_FLOAT, (/landID/), llatvID)
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error defining land lat variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       ok = NF90_PUT_ATT(ncid_out, llatvID, 'units', "degrees_north")
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error defining local lat variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ok = NF90_DEF_VAR(ncid_out, 'local_lon', NF90_FLOAT, (/landID/), llonvID)
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error defining land lon variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       ok = NF90_PUT_ATT(ncid_out, llonvID, 'units', "degrees_east")
       IF (ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error defining local lon variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')

    END IF
    ! Define "time" variable and its attributes:
    ok = NF90_DEF_VAR(ncid_out, 'time', NF90_DOUBLE, (/tID/), ovid%tvar)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'calendar', calendar)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')

    ! Define latitude and longitude variable (ALMA):
    ok = NF90_DEF_VAR(ncid_out, 'latitude', NF90_FLOAT, (/xID, yID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining latitude variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, latID, 'units', 'degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining latitude variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_VAR(ncid_out, 'longitude', NF90_FLOAT, (/xID, yID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining longitude variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining longitude variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Write "cordinate variables" to enable reading by GrADS:
    ok = NF90_DEF_VAR(ncid_out, 'x', NF90_FLOAT, (/xID/), xvID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining "x" variable (for GrADS) in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, xvID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error writing x coordinate variable (GrADS) units in output '// &
         'file. (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, xvID, 'comment', &
         'x coordinate variable for GrADS compatibility')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error writing x variables comment in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_VAR(ncid_out, 'y', NF90_FLOAT, (/yID/), yvID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining "y" variable (for GrADS) in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, yvID, 'units', 'degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error writing y coordinate variable (GrADS) units in output '// &
         'file. (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, yvID, 'comment', &
         'y coordinate variable for GrADS compatibility')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error writing y variables comment in output file. '// &
         '(SUBROUTINE open_output_file)')
    !   ! Define fraction of each surface type:
    !   CALL define_ovar(ncid_out,surffracID,'surffrac','-', &
    !       'Fraction of each surface type: vegetated; urban; lake; land ice', &
    !       .FALSE.,surftypeID,'surftype',xID,yID,zID,landID,patchID)

    !=============DEFINE OUTPUT VARIABLES=======================================
    ! Define met forcing variables in output file and allocate temp output vars:
    IF(output%met .OR. output%SWdown) THEN
       CALL define_ovar(ncid_out, &
            ovid%SWdown, 'SWdown', 'W/m^2', 'Downward shortwave radiation', &
            patchout%SWdown, 'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SWdown(mp))
       out%SWdown = zero4 ! initialise
    END IF
    IF(output%met .OR. output%FracDiff) THEN
     CALL define_ovar(ncid_out, &
          ovid%FracDiff, 'FracDiff', '-', 'Fraction of diffuse radiation', &
          patchout%FracDiff, 'dummy', xID, yID, zID, landID, patchID, tID)
     ALLOCATE(out%FracDiff(mp))
     out%FracDiff = zero4 ! initialise
    END IF
    IF(output%met .OR. output%LWdown) THEN
       CALL define_ovar(ncid_out, ovid%LWdown, 'LWdown', 'W/m^2', &
            'Downward longwave radiation', patchout%LWdown, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LWdown(mp))
       out%LWdown = zero4 ! initialise
    END IF
    IF(output%met .OR. output%Tair) THEN
       CALL define_ovar(ncid_out, ovid%Tair, &
            'Tair', 'K', 'Surface air temperature', patchout%Tair, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Tair(mp))
       out%Tair = zero4 ! initialise
    END IF
    IF(output%met .OR. output%Rainf) THEN
       CALL define_ovar(ncid_out, ovid%Rainf, 'Rainf', &
            'kg/m^2/s', 'Rainfall+snowfall', patchout%Rainf, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Rainf(mp))
       out%Rainf = zero4 ! initialise
    END IF
    IF(output%met .OR. output%Snowf) THEN
       CALL define_ovar(ncid_out, ovid%Snowf, 'Snowf', &
            'kg/m^2/s', 'Snowfall', patchout%Snowf, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Snowf(mp))
       out%Snowf = zero4 ! initialise
    END IF
    IF(output%met .OR. output%Qair) THEN
       CALL define_ovar(ncid_out, ovid%Qair, 'Qair', &
            'kg/kg', 'Surface specific humidity', patchout%Qair, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qair(mp))
       out%Qair = zero4 ! initialise
    END IF
    IF(output%met .OR. output%Wind) THEN
       CALL define_ovar(ncid_out, ovid%Wind, 'Wind', &
            'm/s', 'Scalar surface wind speed', patchout%Wind, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Wind(mp))
       out%Wind = zero4 ! initialise
    END IF
    IF(output%met .OR. output%PSurf) THEN
       CALL define_ovar(ncid_out, ovid%PSurf, 'PSurf', &
            'hPa', 'Surface air pressure', patchout%PSurf, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PSurf(mp))
       out%PSurf = zero4 ! initialise
    END IF
    IF(output%met .OR. output%CO2air) THEN
       CALL define_ovar(ncid_out, ovid%CO2air, 'CO2air', 'ppmv', &
            'Surface air CO2 concentration', patchout%CO2air, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CO2air(mp))
       out%CO2air = zero4 ! initialise
    END IF
    ! Define surface flux variables in output file and allocate temp output
    ! vars:
    IF(output%flux .OR. output%Qle) THEN
       CALL define_ovar(ncid_out, ovid%Qle, 'Qle', 'W/m^2', &
            'Surface latent heat flux',patchout%Qle,'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qle(mp))
       out%Qle = zero4 ! initialise
    END IF
    IF(output%flux .OR. output%Qh) THEN
       CALL define_ovar(ncid_out,ovid%Qh,'Qh', 'W/m^2', &
            'Surface sensible heat flux',patchout%Qh,'dummy', &
            xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%Qh(mp))
       out%Qh = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%Qg) THEN
       CALL define_ovar(ncid_out, ovid%Qg, 'Qg', 'W/m^2', &
            'Surface ground heat flux', patchout%Qg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qg(mp))
       out%Qg = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%Qs) THEN
       CALL define_ovar(ncid_out, ovid%Qs, 'Qs', &
            'kg/m^2/s', 'Surface runoff', patchout%Qs, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qs(mp))
       out%Qs = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%Qsb) THEN
       CALL define_ovar(ncid_out, ovid%Qsb, 'Qsb', 'kg/m^2/s', &
            'Subsurface runoff', patchout%Qsb, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qsb(mp))
       out%Qsb = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%Evap) THEN
       CALL define_ovar(ncid_out, ovid%Evap,'Evap', 'kg/m^2/s', &
            'Total evapotranspiration', patchout%Evap, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Evap(mp))
       out%Evap = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%ECanop) THEN
       CALL define_ovar(ncid_out, ovid%Ecanop, 'ECanop', 'kg/m^2/s', &
            'Wet canopy evaporation', patchout%ECanop, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ECanop(mp))
       out%ECanop = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%TVeg) THEN
       CALL define_ovar(ncid_out, ovid%TVeg, 'TVeg', 'kg/m^2/s', &
            'Vegetation transpiration', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TVeg(mp))
       out%TVeg = zero4 ! initialise
    END IF

    IF(output%flux) THEN
       CALL define_ovar(ncid_out, ovid%gsw_sl, 'gsw_sl', 'mol/m^2/s', &
            'stomatal conductance sl leaves', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%gsw_sl(mp))
       out%gsw_sl = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%gsw_sh, 'gsw_sh', 'mol/m^2/s', &
            'stomatal conductance sh leaves', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%gsw_sh(mp))
       out%gsw_sh = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%scalex_sl, 'scale_factor_sunlit', '[ ]', &
            'canopy scaling factor sl leaves', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%scalex_sl(mp))
       out%scalex_sl = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%scalex_sh, 'scale_factor_shaded', '[ ]', &
            'canopy scaling factor sh leaves', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%scalex_sh(mp))
       out%scalex_sh = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%dlf, 'leaf_to_air_vpd', 'kPa', &
            'leaf to air vapour pressure difference', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%dlf(mp))
       out%dlf = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%qcan_sl, 'qcan_sl', 'W/m^2', &
            'absorbed radiation by canopy, sunlit leaves', patchout%Qcan, 'radiation', &
            xID, yID, zID, landID, patchID, radID, tID)
       ALLOCATE(out%qcan_sl(mp,nrb))
       out%qcan_sl = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%qcan_sh, 'qcan_sh', 'W/m^2', &
            'absorbed radiation by canopy, shaded leaves', patchout%Qcan, 'radiation', &
            xID, yID, zID, landID, patchID, radID, tID)
       ALLOCATE(out%qcan_sh(mp,nrb))
       out%qcan_sh = zero4 ! initialise
    END IF

    IF(output%flux .OR. output%ESoil) THEN
       CALL define_ovar(ncid_out, ovid%ESoil, 'ESoil', 'kg/m^2/s', &
            'Evaporation from soil', patchout%ESoil, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ESoil(mp))
       out%ESoil = zero4 ! initialise
    END IF
    IF(output%flux .OR. output%HVeg) THEN
       CALL define_ovar(ncid_out, ovid%HVeg, 'HVeg', 'W/m^2', &
            'Sensible heat from vegetation', patchout%HVeg, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HVeg(mp))
       out%HVeg = zero4 ! initialise
    END IF
    IF(output%flux .OR. output%HSoil) THEN
       CALL define_ovar(ncid_out, ovid%HSoil, 'HSoil', 'W/m^2', &
            'Sensible heat from soil', patchout%HSoil, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HSoil(mp))
       out%HSoil = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%RnetSoil, 'RnetSoil', 'W/m^2', &
            'Net radiation absorbed by ground', patchout%RnetSoil, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RnetSoil(mp))
       out%RnetSoil = zero4 ! initialise
    END IF
    IF(output%flux .OR. output%carbon .OR. output%NEE) THEN
       CALL define_ovar(ncid_out, ovid%NEE, 'NEE', 'umol/m^2/s', &
            'Net ecosystem exchange of CO2', patchout%NEE, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NEE(mp))
       out%NEE = zero4 ! initialise
    END IF


    ! Define soil state variables in output file and allocate temp output vars:
    IF(output%soil .OR. output%SoilMoist) THEN
       CALL define_ovar(ncid_out, ovid%SoilMoist, 'SoilMoist', 'm^3/m^3', &
            'Average layer soil moisture', patchout%SoilMoist, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       CALL define_ovar(ncid_out, ovid%SoilMoistPFT, 'SoilMoistPFT', 'kg m-2', &
            'Average soil moisture per PFT (whole soil column)', patchout%SoilMoistPFT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%SoilMoistIce, 'SoilMoistIce', 'm^3/m^3', &
            'Average layer frozen soil moisture', patchout%SoilMoistIce, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilMoist(mp,ms))
       ALLOCATE(out%SoilMoistPFT(mp))
       ALLOCATE(out%SoilMoistIce(mp,ms))
       out%SoilMoist = zero4 ! initialise
       out%SoilMoistPFT = zero4 ! initialise
       out%SoilMoistIce = zero4 ! initialise
    END IF
    IF(output%soil .OR. output%SoilTemp) THEN
       CALL define_ovar(ncid_out, ovid%SoilTemp, 'SoilTemp', 'K', &
            'Average layer soil temperature', patchout%SoilTemp, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilTemp(mp,ms))
       out%SoilTemp = zero4 ! initialise
    END IF
    IF(output%soil .OR. output%BaresoilT) THEN
       CALL define_ovar(ncid_out, ovid%BaresoilT, 'BaresoilT', &
            'K', 'Bare soil temperature', patchout%BaresoilT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%BaresoilT(mp))
       out%BaresoilT = zero4 ! initialise
    END IF
    ! Define snow state variables in output file and allocate temp output vars:
    IF(output%snow .OR. output%SWE) THEN
       CALL define_ovar(ncid_out, ovid%SWE, 'SWE', 'kg/m^2', &
            'Snow water equivalent', patchout%SWE, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SWE(mp))
       out%SWE = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%SnowMelt, 'SnowMelt', 'kg/m^2/s', &
            'Snow Melt Rate', patchout%SnowMelt, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowMelt(mp))
       out%SnowMelt = zero4 ! initialise
    END IF
    IF(output%snow .OR. output%SnowT) THEN
       CALL define_ovar(ncid_out, ovid%SnowT, 'SnowT', 'K', &
            'Snow surface temperature', patchout%SnowT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowT(mp))
       out%SnowT = zero4 ! initialise
    END IF
    IF(output%snow .OR. output%SnowDepth) THEN
       CALL define_ovar(ncid_out, ovid%SnowDepth, 'SnowDepth', &
            'm', 'Snow depth', patchout%SnowDepth, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowDepth(mp))
       out%SnowDepth = zero4 ! initialise
    END IF
    ! Define radiative variables in output file and allocate temp output vars:
    IF(output%radiation) THEN
       CALL define_ovar(ncid_out, ovid%fbeam, 'fbeam', '-', &
            'fraction of direct radiation (visible)', &
            patchout%fbeam, 'dummy', xID, yID, zID, landID, &
            patchID, tID)
       ALLOCATE(out%fbeam(mp))
       out%fbeam = zero4 ! initialise
    END IF
    IF(output%radiation .OR. output%SWnet) THEN
       CALL define_ovar(ncid_out, ovid%SWnet, 'SWnet', 'W/m^2', &
            'Net shortwave radiation absorbed by surface', &
            patchout%SWnet, 'dummy', xID, yID, zID, landID, &
            patchID, tID)
       ALLOCATE(out%SWnet(mp))
       out%SWnet = zero4 ! initialise
    END IF
    IF(output%radiation .OR. output%LWnet) THEN
       CALL define_ovar(ncid_out, ovid%LWnet, 'LWnet', 'W/m^2', &
            'Net longwave radiation absorbed by surface', &
            patchout%LWnet, 'dummy', xID, yID, zID, landID, &
            patchID, tID)
       ALLOCATE(out%LWnet(mp))
       out%LWnet = zero4 ! initialise
    END IF
    IF(output%radiation .OR. output%Rnet) THEN
       CALL define_ovar(ncid_out, ovid%Rnet, 'Rnet', 'W/m^2', &
            'Net radiation absorbed by surface', patchout%Rnet, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Rnet(mp))
       out%Rnet = zero4 ! initialise
    END IF
    IF(output%radiation .OR. output%Albedo) THEN
       CALL define_ovar(ncid_out, ovid%Albedo, 'Albedo', '-', &
            'Surface albedo', patchout%Albedo, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Albedo(mp))
       out%Albedo = zero4 ! initialise
    END IF

    ! output calc of soil albedo based on colour? - Ticket #27
    !IF (calcsoilalbedo) THEN
    IF(output%radiation .OR. output%visAlbedo) THEN
       CALL define_ovar(ncid_out, ovid%visAlbedo, 'visAlbedo', '-', &
            'Surface vis albedo:total, beam, diffuse', patchout%visAlbedo, &
            'radiation',xID, yID, zID, landID, patchID,radID, tID)
       ALLOCATE(out%visAlbedo(mp,nrb))
       out%visAlbedo = zero4 ! initialise
    END IF
    IF(output%radiation .OR. output%nirAlbedo) THEN
       CALL define_ovar(ncid_out, ovid%nirAlbedo, 'nirAlbedo', '-', &
            'Surface nir albedo:total, beam, diffuse', patchout%nirAlbedo, &
            'radiation',xID, yID, zID, landID, patchID,radID, tID)
       ALLOCATE(out%nirAlbedo(mp,nrb))
       out%nirAlbedo = zero4 ! initialise
    END IF

    !END IF

    IF(output%radiation .OR. output%RadT) THEN
       CALL define_ovar(ncid_out, ovid%RadT, 'RadT', 'K', &
            'Radiative surface temperature', patchout%RadT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RadT(mp))
       out%RadT = zero4 ! initialise
    END IF
    ! Define vegetation variables in output file and allocate temp output vars:
    IF(output%veg .OR. output%VegT) THEN
       CALL define_ovar(ncid_out, ovid%VegT, 'VegT', 'K', &
            'Average vegetation temperature', patchout%VegT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%VegT(mp))
       out%VegT = zero4 ! initialise
    END IF
    IF(output%veg .OR. output%CanT) THEN
       CALL define_ovar(ncid_out, ovid%CanT, 'CanT', 'K', &
            'Within-canopy temperature', patchout%CanT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CanT(mp))
       out%CanT = zero4 ! initialise
    END IF
    IF(output%veg .OR. output%Fwsoil) THEN
       CALL define_ovar(ncid_out, ovid%Fwsoil, 'Fwsoil', '[-]', &
            'soil moisture modifier to stomatal conductance', patchout%Fwsoil, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Fwsoil(mp))
       out%Fwsoil = zero4 ! initialise
    END IF
    IF(output%veg .OR. output%CanopInt) THEN
       CALL define_ovar(ncid_out, ovid%CanopInt, 'CanopInt', 'kg/m^2', &
            'Canopy intercepted water storage', patchout%CanopInt, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CanopInt(mp))
       out%CanopInt = zero4 ! initialise
    END IF
    IF(output%veg .OR. output%LAI) THEN
       CALL define_ovar(ncid_out, ovid%LAI, 'LAI', '-', &
            'Leaf area index', patchout%LAI, 'dummy', xID, &
            yID, zID, landID, patchID, tID)
       ALLOCATE(out%LAI(mp))
       out%LAI = zero4 ! initialise
    END IF

    ! Alexis
    IF(output%veg) THEN
       CALL define_ovar(ncid_out, ovid%vcmax, 'vcmax', '-', &
            'Vcmax at 25 degC [molC/m^2/s]', patchout%LAI, 'dummy', xID, &
            yID, zID, landID, patchID, tID)
       ALLOCATE(out%vcmax(mp))
       out%vcmax = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%ejmax, 'jmax', '-', &
            'Jmax at 25 degC [mol(e)/m^2/s]', patchout%LAI, 'dummy', xID, &
            yID, zID, landID, patchID, tID)
       ALLOCATE(out%ejmax(mp))
       out%ejmax = zero4 ! initialise
    END IF
    ! Alexis

    ! Define balance variables in output file and allocate temp output vars:
    IF(output%balances .OR. output%Ebal) THEN
       CALL define_ovar(ncid_out, ovid%Ebal, 'Ebal', 'W/m^2', &
            'Cumulative energy imbalance', patchout%Ebal, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Ebal(mp))
       out%Ebal = zero4 ! initialise
    END IF
    IF(output%balances .OR. output%Wbal) THEN
       CALL define_ovar(ncid_out, ovid%Wbal, 'Wbal', 'kg/m^2', &
            'Cumulative water imbalance', patchout%Wbal, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Wbal(mp))
       out%Wbal = zero4 ! initialise
    END IF
    ! Define carbon variables in output file and allocate temp output vars:
    IF(output%carbon .OR. output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%AutoResp, 'AutoResp', 'umol/m^2/s', &
            'Autotrophic respiration', patchout%AutoResp, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%AutoResp(mp))
       out%AutoResp = zero4 ! initialise
    END IF
    IF(output%casa .OR. output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%RootResp, 'RootResp', 'umol/m^2/s', &
            'Fine Root Autotrophic respiration', patchout%AutoResp, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RootResp(mp))
       out%RootResp = zero4 ! initialise
    END IF

    IF(output%casa .OR. output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%StemResp, 'StemResp', 'umol/m^2/s', &
            'StemWood Autotrophic respiration', patchout%AutoResp, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%StemResp(mp))
       out%StemResp = zero4 ! initialise
    END IF

    IF(output%carbon .OR. output%LeafResp) THEN
       CALL define_ovar(ncid_out, ovid%LeafResp, 'LeafResp', 'umol/m^2/s', &
            'Leaf respiration', patchout%LeafResp, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LeafResp(mp))
       out%LeafResp = zero4 ! initialise
    END IF
    IF(output%carbon .OR. output%HeteroResp) THEN
       CALL define_ovar(ncid_out, ovid%HeteroResp, 'HeteroResp', 'umol/m^2/s', &
            'Heterotrophic respiration', patchout%HeteroResp, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HeteroResp(mp))
       out%HeteroResp = zero4 ! initialise
    END IF
    IF(output%carbon.OR.output%GPP) THEN
       CALL define_ovar(ncid_out, ovid%GPP, 'GPP', 'umol/m^2/s', &
            'Gross primary production', patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP(mp))
       out%GPP = zero4 ! initialise
    END IF
    IF(output%GPP_components) THEN
       CALL define_ovar(ncid_out, ovid%An_sl, 'Anet_sunlit', 'umol/m^2/s', &
            'Anet from sunlit leaves', patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%An_sl(mp))
       out%An_sl = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%An_sh, 'Anet_shaded', 'umol/m^2/s', &
            'Anet from shaded leaves', patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%An_sh(mp))
       out%An_sh = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%ci_sl, 'ci_sl', 'ppm', &
            'ci, sunlit leaves', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ci_sl(mp))
       out%ci_sl = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%ci_sh, 'ci_sh', 'ppm', &
            'ci, sunlit leaves', patchout%TVeg, 'dummy', &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ci_sh(mp))
       out%ci_sh = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%GPP_sl, 'GPP_sunlit', 'umol/m^2/s', &
            'Gross primary production from sunlit leaves', patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP_sl(mp))
       out%GPP_sl = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%GPP_sh, 'GPP_shaded', 'umol/m^2/s', &
            'Gross primary production from shaded leaves', patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP_sh(mp))
       out%GPP_sh = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%GPP_slC, 'GPP_sunlit_C', 'umol/m^2/s', &
            'Gross primary production from carboxylation-rate-limited sunlit leaves', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP_slC(mp))
       out%GPP_slC = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%GPP_shC, 'GPP_shaded_C', 'umol/m^2/s', &
            'Gross primary production from carboxylation-rate-limited shaded leaves', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP_shC(mp))
       out%GPP_shC = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%GPP_slJ, 'GPP_sunlit_J', 'umol/m^2/s', &
            'Gross primary production from electron-transport-rate-limited sunlit leaves', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP_slJ(mp))
       out%GPP_slJ = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%GPP_shJ, 'GPP_shaded_J', 'umol/m^2/s', &
            'Gross primary production from electron-transport-rate-limited shaded leaves', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP_shJ(mp))
       out%GPP_shJ = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%eta_GPP_cs, 'eta_GPP_cs', 'umol/m^2/s', &
            'elasticity of Gross primary production wrt cs, multiplied by GPP', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%eta_GPP_cs(mp))
       out%eta_GPP_cs = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%eta_TVeg_cs, 'eta_TVeg_cs', 'kg/m^2/s', &
            'elasticity of Transpiration wrt cs, multiplied by Transpiration', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%eta_TVeg_cs(mp))
       out%eta_TVeg_cs = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%gsw_TVeg, 'gsw_TVeg', 'mol/m^s/s * kg/m^2/s', &
            'stomatal conductance, multiplied by Transpiration', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%gsw_TVeg(mp))
       out%gsw_TVeg = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%dGPPdcs, 'dGPPdcs', '(umol/m^2/s)^2', &
            'sensitivity of Gross primary production wrt cs, multiplied by GPP', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%dGPPdcs(mp))
       out%dGPPdcs = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%CO2s, 'CO2s', 'ppm umol/m^2/s', &
            'CO2 concentration at leaf surface , multiplied by GPP', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CO2s(mp))
       out%CO2s = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%vcmax_ts, 'vcmax_weighted', 'mol/m^2/s', &
            'vcmax weighted by sunlit/shaded LAI', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%vcmax_ts(mp))
       out%vcmax_ts = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%jmax_ts, 'jmax_weighted', 'mol/m^2/s', &
            'jmax weighted by sunlit/shaded LAI', &
            patchout%GPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%jmax_ts(mp))
       out%jmax_ts = zero4 ! initialise
    END IF

    ! vh_mc ! additional variables for ESM-SnowMIP
    IF (output%snowmip) THEN
       ! define
       CALL define_ovar(ncid_out, ovid%hfds, 'hfds', 'W/m2', &
            'downward heat flux at ground surface', patchout%hfds, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%hfdsn, 'hfdsn', 'W/m2', &
            'downward heat flux into snowpack', patchout%hfdsn, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%hfls, 'hfls', 'W/m2', &
            'surface upward latent heat flux', patchout%hfls, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%hfmlt, 'hfmlt', 'W/m2', &
            'energy of fusion', patchout%hfmlt, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%hfrs, 'hfrs', 'W/m2', &
            'heat transferred to snowpack by rain', patchout%hfrs, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%hfsbl, 'hfsbl', 'W/m2', &
            'energy of sublimation', patchout%hfsbl, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%hfss, 'hfss', 'W/m2', &
            'surface upward sensible heat flux', patchout%hfss, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%rlus, 'rlus', 'W/m2', &
            'surface upwelling longwave radiation', patchout%rlus, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%rsus, 'rsus', 'W/m2', &
            'surface upwelling shortwave radiation', patchout%rsus, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%esn, 'esn', 'kg/m2/s', &
            'liquid water evaporation from snowpack', patchout%esn, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%evspsbl, 'evspsbl', 'kg/m2/s', &
            'total water vapour flux from the surface to the atmosphere', patchout%evspsbl, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%evspsblsoi, 'evspsblsoi', 'kg/m2/s', &
            'evaporation and sublimation from soil', patchout%evspsblsoi, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%evspsblveg, 'evspsblveg', 'kg/m2/s', &
            'evaporation and sublimation from canopy', patchout%evspsblveg, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%mrrob, 'mrrob', 'kg/m2/s', &
            'subsurface runoff', patchout%mrrob, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%mrros, 'mrros', 'kg/m2/s', &
            'surface runoff', patchout%mrros, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%sbl, 'sbl', 'kg/m2/s', &
            'sublimation of snow', patchout%sbl, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%snm, 'snm', 'kg/m2/s', &
            'surface snow melt', patchout%snm, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%snmsl, 'snmsl', 'kg/m2/s', &
            'water flowing out of snowpack', patchout%snmsl, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%tran, 'tran', 'kg/m2/s', &
            'transpiration', patchout%tran, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%albs, 'albs', '-', &
            'surface albedo', patchout%albs, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%albsn, 'albsn', '-', &
            'snow albedo', patchout%albsn, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%cw, 'cw', 'kg/m2', &
            'total canopy water storage', patchout%cw, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%lqsn, 'lqsn', '-', &
            'mass fraction of liquid water in snowpack', patchout%lqsn, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%lwsnl, 'lwsnl', 'kg/m2', &
            'liquid water content of snowpack', patchout%lwsnl, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%mrfsofr, 'mrfsofr', '-', &
            'mass fractions of frozen water in soil layers', patchout%mrfsofr, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       CALL define_ovar(ncid_out, ovid%mrlqso, 'mrlqso', '-', &
            'mass fractions of unfrozen water in soil layers', patchout%mrlqso, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       CALL define_ovar(ncid_out, ovid%mrlsl, 'mrlsl', 'kg/m2', &
            'masses of frozen and unfrozen moisture in soil layers', patchout%mrlsl, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       CALL define_ovar(ncid_out, ovid%snc, 'snc', '-', &
            'snow area fraction', patchout%snc, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%snd, 'snd', 'm', &
            'snowdepth', patchout%snd, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%snw, 'snw', 'kg/m2', &
            'mass of snowpack', patchout%snw, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%snwc, 'snwc', 'kg/m2', &
            'mass of snow intercepted by vegetation', patchout%snwc, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%tcs, 'tcs', 'K', &
            'vegetation canopy temperature', patchout%tcs, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%tgs, 'tgs', 'K', &
            'temperature of bare soil', patchout%tgs, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%ts, 'ts', 'K', &
            'surface temperature', patchout%ts, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%tsl, 'tsl', 'K', &
            'temperatures of soil layers', patchout%tsl, &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       CALL define_ovar(ncid_out, ovid%tsn, 'tsn', 'K', &
            'snow internal temperature', patchout%tsn, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       CALL define_ovar(ncid_out, ovid%tsns, 'tsns', 'K', &
            'snow surface temperature', patchout%tsns, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ! allocate
       allocate(out%hfds(mp))
       allocate(out%hfdsn(mp))
       allocate(out%hfls(mp))
       allocate(out%hfmlt(mp))
       allocate(out%hfrs(mp))
       allocate(out%hfsbl(mp))
       allocate(out%hfss(mp))
       allocate(out%rlus(mp))
       allocate(out%rsus(mp))
       allocate(out%esn(mp))
       allocate(out%evspsbl(mp))
       allocate(out%evspsblsoi(mp))
       allocate(out%evspsblveg(mp))
       allocate(out%mrrob(mp))
       allocate(out%mrros(mp))
       allocate(out%sbl(mp))
       allocate(out%snm(mp))
       allocate(out%snmsl(mp))
       allocate(out%tran(mp))
       allocate(out%albs(mp))
       allocate(out%albsn(mp))
       allocate(out%cw(mp))
       allocate(out%lqsn(mp))
       allocate(out%lwsnl(mp))
       allocate(out%mrfsofr(mp,ms))
       allocate(out%mrlqso(mp,ms))
       allocate(out%mrlsl(mp,ms))
       allocate(out%snc(mp))
       allocate(out%snd(mp))
       allocate(out%snw(mp))
       allocate(out%snwc(mp))
       allocate(out%tcs(mp))
       allocate(out%tgs(mp))
       allocate(out%ts(mp))
       allocate(out%tsl(mp,ms))
       allocate(out%tsn(mp))
       allocate(out%tsns(mp))
       ! initialise
       out%hfds       = zero4
       out%hfdsn      = zero4
       out%hfls       = zero4
       out%hfmlt      = zero4
       out%hfrs       = zero4
       out%hfsbl      = zero4
       out%hfss       = zero4
       out%rlus       = zero4
       out%rsus       = zero4
       out%esn        = zero4
       out%evspsbl    = zero4
       out%evspsblsoi = zero4
       out%evspsblveg = zero4
       out%mrrob      = zero4
       out%mrros      = zero4
       out%sbl        = zero4
       out%snm        = zero4
       out%snmsl      = zero4
       out%tran       = zero4
       out%albs       = zero4
       out%albsn      = zero4
       out%cw         = zero4
       out%lqsn       = zero4
       out%lwsnl      = zero4
       out%mrfsofr    = zero4
       out%mrlqso     = zero4
       out%mrlsl      = zero4
       out%snc        = zero4
       out%snd        = zero4
       out%snw        = zero4
       out%snwc       = zero4
       out%tcs        = zero4
       out%tgs        = zero4
       out%ts         = zero4
       out%tsl        = zero4
       out%tsn        = zero4
       out%tsns       = zero4
    END IF

    IF(output%carbon .OR. output%NPP) THEN
       CALL define_ovar(ncid_out, ovid%NPP, 'NPP', 'umol/m^2/s', &
            'Net primary production', patchout%NPP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NPP(mp))
       out%NPP = zero4 ! initialise
    END IF

    IF(output%carbon .OR. output%NBP) THEN
       CALL define_ovar(ncid_out, ovid%NBP, 'NBP', 'umol/m^2/s', &
            'Net Biosphere Production (uptake +ve)' &
            , patchout%NBP, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NBP(mp))
       out%NBP = zero4 ! initialise
    ENDIF

    IF (output%casa) THEN

       CALL define_ovar(ncid_out, ovid%dCdt, 'dCdt', 'umol/m^2/s', &
            'Carbon accumulation rate (uptake +ve)', patchout%dCdt, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%dCdt(mp))
       out%dCdt = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%TotSoilCarb, 'TotSoilCarb', 'kg C/m^2', &
            'Total Soil and Litter Carbon', patchout%TotSoilCarb, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TotSoilCarb(mp))
       out%TotSoilCarb = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%TotLittCarb, 'TotLittCarb', 'kg C/m^2', &
            'Total Litter Carbon', patchout%TotLittCarb, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TotLittCarb(mp))
       out%TotLittCarb = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%SoilCarbFast, 'SoilCarbFast', 'kg C/m^2', &
            'Soil Carbon: Fast Turnover', patchout%SoilCarbFast, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SoilCarbFast(mp))
       out%SoilCarbFast = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%SoilCarbSlow, 'SoilCarbSlow', 'kg C/m^2', &
            'Soil Carbon: Slow Turnover', patchout%SoilCarbSlow, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SoilCarbSlow(mp))
       out%SoilCarbSlow = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%SoilCarbPassive, 'SoilCarbPassive', 'kg C/m^2', &
            'Soil Carbon: Passive', patchout%SoilCarbPassive, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SoilCarbPassive(mp))
       out%SoilCarbPassive = zero4 ! initialise


       CALL define_ovar(ncid_out, ovid%LittCarbMetabolic, 'LittCarbMetabolic', 'kg C/m^2', &
            'Litter Carbon: metabolic', patchout%LittCarbMetabolic, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LittCarbMetabolic(mp))
       out%LittCarbMetabolic = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%LittCarbStructural, 'LittCarbStructural', 'kg C/m^2', &
            'Litter Carbon: structural', patchout%LittCarbStructural, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LittCarbStructural(mp))
       out%LittCarbStructural = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%LittCarbCWD, 'LittCarbCWD', 'kg C/m^2', &
            'Litter Carbon: CWD', patchout%LittCarbCWD, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LittCarbCWD(mp))
       out%LittCarbCWD = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%PlantCarbLeaf, 'PlantCarbLeaf', 'kg C/m^2', &
            'Plant Carbon: leaf', patchout%PlantCarbLeaf, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantCarbLeaf(mp))
       out%PlantCarbLeaf = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%PlantCarbFineRoot, 'PlantCarbFineRoot', 'kg C/m^2', &
            'Plant Carbon: Fine roots', patchout%PlantCarbFineRoot, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantCarbFineRoot(mp))
       out%PlantCarbFineRoot = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%PlantCarbWood, 'PlantCarbWood', 'kg C/m^2', &
            'Plant Carbon: wood (above- and below-ground', patchout%PlantCarbWood, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantCarbWood(mp))
       out%PlantCarbWood = zero4 ! initialise

       CALL define_ovar(ncid_out, ovid%TotLivBiomass, 'TotLivBiomass', 'kg C/m^2', &
            'Total Biomass', patchout%TotLivBiomass, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TotLivBiomass(mp))
       out%TotLivBiomass = zero4 ! initialise


       CALL define_ovar(ncid_out, ovid%PlantTurnover, 'PlantTurnover', 'umol/m^2/s', &
            'Total Biomass Turnover', patchout%PlantTurnover, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnover(mp))
       out%PlantTurnover = zero4

       CALL define_ovar(ncid_out, ovid%PlantTurnoverLeaf, 'PlantTurnoverLeaf ', &
            'umol/m^2/s', &
            'Leaf Biomass Turnover', patchout%PlantTurnoverLeaf, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverLeaf(mp))
       out%PlantTurnoverLeaf = zero4

       CALL define_ovar(ncid_out, ovid%PlantTurnoverFineRoot, 'PlantTurnoverFineRoot ', &
            'umol/m^2/s', &
            'FineRoot Biomass Turnover', patchout%PlantTurnoverFineRoot, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverFineRoot(mp))
       out%PlantTurnoverFineRoot = zero4

       CALL define_ovar(ncid_out, ovid%PlantTurnoverWood, 'PlantTurnoverWood ', &
            'umol/m^2/s', &
            'Woody Biomass Turnover', patchout%PlantTurnoverWood, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWood(mp))
       out%PlantTurnoverWood = zero4

       CALL define_ovar(ncid_out, ovid%PlantTurnoverWoodDist, 'PlantTurnoverWoodDist ', &
            'umol/m^2/s', &
            'Woody Biomass Turnover (disturbance)', patchout%PlantTurnoverWoodDist, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWoodDist(mp))
       out%PlantTurnoverWoodDist = zero4

       CALL define_ovar(ncid_out, ovid%PlantTurnoverWoodCrowding, 'PlantTurnoverWoodCrowding ', &
            'umol/m^2/s', &
            'Woody Biomass Turnover (crowding)', patchout%PlantTurnoverWoodCrowding, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWoodCrowding(mp))
       out%PlantTurnoverWoodCrowding = zero4

       CALL define_ovar(ncid_out, ovid%PlantTurnoverWoodResourceLim, 'PlantTurnoverWoodResourceLim ', &
            'umol/m^2/s', &
            'Woody Biomass Turnover (Resource Limitation)', patchout%PlantTurnoverWoodResourceLim, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWoodResourceLim(mp))
       out%PlantTurnoverWoodResourceLim = zero4

       IF (cable_user%POPLUC) THEN
          CALL define_ovar(ncid_out, ovid%LandUseFlux, 'LandUseFlux ', &
               'umol/m^2/s', &
               'Sum of wood harvest and clearing fluxes', patchout%LandUseFlux, &
               'dummy', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%LandUseFlux(mp))
          out%LandUseFlux = zero4
       ENDIF

    ENDIF

    ! 13C
    if (cable_user%c13o2 .and. output%c13o2) then
       ! 12C
       allocate(out%An(mp))
       allocate(out%Rd(mp))
       allocate(out%cplant(mp,mplant))
       allocate(out%clitter(mp,mlitter))
       allocate(out%csoil(mp,msoil))
       allocate(out%clabile(mp))
       ! 13C
       allocate(out%A13n(mp))
       allocate(out%aDisc13(mp))
       allocate(out%c13plant(mp,mplant))
       allocate(out%c13litter(mp,mlitter))
       allocate(out%c13soil(mp,msoil))
       allocate(out%c13labile(mp))
       ! 12C
       out%An      = 0.0_r_2
       out%Rd      = 0.0_r_2
       out%cplant  = 0.0_r_2
       out%clitter = 0.0_r_2
       out%csoil   = 0.0_r_2
       out%clabile = 0.0_r_2
       ! 13C
       out%A13n      = 0.0_r_2
       out%aDisc13   = 0.0_r_2
       out%c13plant  = 0.0_r_2
       out%c13litter = 0.0_r_2
       out%c13soil   = 0.0_r_2
       out%c13labile = 0.0_r_2
       ! all defined as float
       ! 12C
       call define_ovar(ncid_out, ovid%An, 'An', 'mol/m^2/s', &
            'Net assimilation', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, tID)
       call define_ovar(ncid_out, ovid%Rd, 'Rd', 'mol/m^2/s', &
            'Leaf respiration', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, tID)
       call define_ovar(ncid_out, ovid%cplant, 'Cplant', 'kgC/m^2', &
            'Plant carbon pools', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, nplantid, tID)
       call define_ovar(ncid_out, ovid%clitter, 'Clitter', 'kgC/m^2', &
            'Litter carbon pools', patchout%c13o2, &
                                !'littercarbon', xID, yID, zID, landID, patchID, nlitterid, tID)
            'generic', xID, yID, zID, landID, patchID, nlitterid, tID)
       call define_ovar(ncid_out, ovid%csoil, 'Csoil', 'kgC/m^2', &
            'Soil carbon pools', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, nsoilid, tID)
       call define_ovar(ncid_out, ovid%clabile, 'Clabile', 'kgC/m^2', &
            'Excess carbon pool', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, tID)
       ! 13C
       call define_ovar(ncid_out, ovid%A13n, 'An13', 'mol/m^2/s', &
            'Net assimilation of 13CO2', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, tID)
       call define_ovar(ncid_out, ovid%aDisc13, 'aDisc13', '1*mol/m^2/s', &
            'Leaf discrimination times gross assimilation', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, tID)
       call define_ovar(ncid_out, ovid%c13plant, 'C13plant', 'kg13C/m^2', &
            'Plant carbon pools', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, nplantid, tID)
       call define_ovar(ncid_out, ovid%c13litter, 'C13litter', 'kg13C/m^2', &
            'Litter carbon pools', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, nlitterid, tID)
       call define_ovar(ncid_out, ovid%c13soil, 'C13soil', 'kg13C/m^2', &
            'Soil carbon pools', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, nsoilid, tID)
       call define_ovar(ncid_out, ovid%c13labile, 'C13labile', 'kg13C/m^2', &
            'Excess carbon pool', patchout%c13o2, &
            'generic', xID, yID, zID, landID, patchID, tID)
    endif

    !! vh_js !!
    CALL define_ovar(ncid_out, ovid%Area, 'Area', 'km2', &
         'Patch Area', patchout%Area, &
         'dummy', xID, yID, zID, landID, patchID, tID)
    ALLOCATE(out%Area(mp))
    out%Area = zero4 ! initialise

    ! Define CABLE parameters in output file:
    IF(output%params .OR. output%iveg) CALL define_ovar(ncid_out, opid%iveg, &
         'iveg', '-', 'Vegetation type', patchout%iveg, 'integer', &
         xID, yID, zID, landID, patchID)

    IF (cable_user%POPLUC) THEN

       CALL define_ovar(ncid_out, opid%patchfrac, 'patchfrac', '-', &
            'Fractional cover of vegetation patches', patchout%patchfrac, 'real', &
            xID, yID, zID, landID, patchID, tID)

    ELSE

       IF((output%params .OR. output%patchfrac) &
            .AND. (patchout%patchfrac .OR. output%patch)) &
            CALL define_ovar(ncid_out, opid%patchfrac, 'patchfrac', '-', &
            'Fractional cover of vegetation patches', patchout%patchfrac, 'real', &
            xID, yID, zID, landID, patchID)

    ENDIF


    IF(output%params .OR. output%isoil) CALL define_ovar(ncid_out, opid%isoil, &
         'isoil', '-', 'Soil type', patchout%isoil, 'integer', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%bch) CALL define_ovar(ncid_out, opid%bch, &
         'bch', '-', 'Parameter b, Campbell eqn 1985', patchout%bch, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%clay) CALL define_ovar(ncid_out, opid%clay, &
         'clay', '-', 'Fraction of soil which is clay', patchout%clay, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%sand) CALL define_ovar(ncid_out, opid%sand, &
         'sand', '-', 'Fraction of soil which is sand', patchout%sand, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%silt) CALL define_ovar(ncid_out, opid%silt, &
         'silt', '-', 'Fraction of soil which is silt', patchout%silt, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%ssat) CALL define_ovar(ncid_out, opid%ssat, &
         'ssat', '-', 'Fraction of soil volume which is water @ saturation', &
         patchout%ssat, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%sfc) CALL define_ovar(ncid_out, opid%sfc, &
         'sfc', '-', 'Fraction of soil volume which is water @ field capacity', &
         patchout%sfc, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%swilt) CALL define_ovar(ncid_out, opid%swilt, &
         'swilt', '-', 'Fraction of soil volume which is water @ wilting point', &
         patchout%swilt, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%hyds) CALL define_ovar(ncid_out, opid%hyds, &
         'hyds', 'm/s', 'Hydraulic conductivity @ saturation', &
         patchout%hyds, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%sucs) CALL define_ovar(ncid_out, opid%sucs, &
         'sucs', 'm', 'Suction @ saturation', &
         patchout%sucs, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%css) CALL define_ovar(ncid_out, opid%css, &
         'css', 'J/kg/C', 'Heat capacity of soil minerals', &
         patchout%css, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%rhosoil) CALL define_ovar(ncid_out, &
         opid%rhosoil, 'rhosoil', 'kg/m^3', 'Density of soil minerals', &
         patchout%rhosoil, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%rs20) CALL define_ovar(ncid_out, opid%rs20, &
         'rs20', '-', 'Soil respiration coefficient at 20C', &
         patchout%rs20, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .AND. output%albsoil) CALL define_ovar(ncid_out, &
         opid%albsoil, 'albsoil', '-', &
         'Snow free shortwave soil reflectance fraction', &
         patchout%albsoil, radID, 'radiation', xID, yID, zID, landID, patchID)
    !! vh_js !!
    IF (cable_user%CALL_POP) THEN
       IF(output%params .OR. output%hc) CALL define_ovar(ncid_out, opid%hc, &
            'hc', 'm', 'Height of canopy', patchout%hc, &
            'real', xID, yID, zID, landID, patchID,tID)
    ELSE
       IF(output%params .OR. output%hc) CALL define_ovar(ncid_out, opid%hc, &
            'hc', 'm', 'Height of canopy', patchout%hc, &
            'real', xID, yID, zID, landID, patchID)
    ENDIF

    IF(output%params .OR. output%canst1) CALL define_ovar(ncid_out, &
         opid%canst1, 'canst1', 'mm/LAI', 'Max water intercepted by canopy', &
         patchout%canst1, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%dleaf) CALL define_ovar(ncid_out, opid%dleaf, &
         'dleaf', 'm', 'Chararacteristic length of leaf', &
         patchout%dleaf, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%frac4) CALL define_ovar(ncid_out, opid%frac4, &
         'frac4', '-', 'Fraction of plants which are C4', &
         patchout%frac4, 'real', xID, yID, zID, landID, patchID)
    !IF (output%params .OR. output%ejmax) CALL define_ovar(ncid_out, opid%ejmax, &
    !    'ejmax', 'mol/m^2/s', 'Max potential electron transport rate top leaf', &
    !    patchout%ejmax, 'real', xID, yID, zID, landID, patchID)
    ! Alexis
    !IF (output%params .OR. output%vcmax) CALL define_ovar(ncid_out, opid%vcmax, &
    !    'vcmax', 'mol/m^2/s', 'Maximum RuBP carboxylation rate top leaf', &
    !    patchout%vcmax, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%rp20) CALL define_ovar(ncid_out, opid%rp20, &
         'rp20', '-', 'Plant respiration coefficient at 20C', &
         patchout%rp20, 'real', xID, yID, zID, landID, patchID)
    ! Ticket #56
    IF(output%params .OR. output%g0) CALL define_ovar(ncid_out, opid%g0, &
         'g0', '-', 'g0 term in Medlyn Stom Cond. Param', &
         patchout%g0, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%g1) CALL define_ovar(ncid_out, opid%g1, &
         'g1', '-', 'g1 term in Medlyn Stom Cond. Param', &
         patchout%g1, 'real', xID, yID, zID, landID, patchID)
    ! end Ticket #56

    IF(output%params .OR. output%rpcoef) CALL define_ovar(ncid_out, &
         opid%rpcoef, 'rpcoef', '1/C', &
         'Temperature coef nonleaf plant respiration', &
         patchout%rpcoef, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%shelrb) CALL define_ovar(ncid_out, &
         opid%shelrb, 'shelrb', '-', 'Sheltering factor', patchout%shelrb, &
         'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%xfang) CALL define_ovar(ncid_out, opid%xfang, &
         'xfang', '-', 'Leaf angle parameter',patchout%xfang, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%wai) CALL define_ovar(ncid_out, opid%wai, &
         'wai', '-', 'Wood area index', patchout%wai, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%vegcf) CALL define_ovar(ncid_out, opid%vegcf, &
         'vegcf', '-', 'vegcf', patchout%vegcf, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%extkn) CALL define_ovar(ncid_out, opid%extkn, &
         'extkn', '-', 'Nitrogen extinction coef for vert. canopy profile', &
         patchout%extkn, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%tminvj) CALL define_ovar(ncid_out, &
         opid%tminvj, 'tminvj', 'C', &
         'Min temperature for the start of photosynthesis', &
         patchout%tminvj, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%tmaxvj) CALL define_ovar(ncid_out, &
         opid%tmaxvj, 'tmaxvj', 'C', 'Max temperature for photosynthesis', &
         patchout%tmaxvj, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%vbeta) CALL define_ovar(ncid_out, opid%vbeta, &
         'vbeta', '-', 'Stomatal sensitivity to soil water', &
         patchout%vbeta, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%xalbnir) CALL define_ovar(ncid_out, &
         opid%xalbnir, 'xalbnir', '-', 'Modifier for albedo in near ir band', &
         patchout%xalbnir, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%meth) CALL define_ovar(ncid_out, opid%meth, &
         'meth', '-', 'Canopy turbulence parameterisation choice', &
         patchout%meth, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%za) THEN
       CALL define_ovar(ncid_out, opid%za_uv, 'za_uv', 'm', &
            'Reference height (lowest atm. model layer) for momentum', &
            patchout%za, 'real', xID, yID, zID, landID, patchID)
       CALL define_ovar(ncid_out, opid%za_tq, 'za_tq', 'm', &
            'Reference height (lowest atm. model layer) for scalars', &
            patchout%za, 'real', xID, yID, zID, landID, patchID)
    ENDIF
    IF(output%params .AND. output%ratecp) CALL define_ovar(ncid_out, &
         opid%ratecp, 'ratecp', '1/year', 'Plant carbon rate constant', &
         patchout%ratecp, plantcarbID, 'plantcarbon', xID, yID, zID, &
         landID, patchID)
    IF(output%params .AND. output%ratecs) CALL define_ovar(ncid_out, &
         opid%ratecs, 'ratecs', '1/year', 'Soil carbon rate constant', &
         patchout%ratecs, soilcarbID, 'soilcarbon', xID, yID, zID, &
         landID, patchID)
    IF(output%params .AND. output%zse) CALL define_ovar(ncid_out, opid%zse, &
         'zse', 'm', 'Depth of each soil layer', &
         patchout%zse, soilID, 'soil', xID, yID, zID, landID, patchID)
    IF(output%params .AND. output%froot) CALL define_ovar(ncid_out, opid%froot, &
         'froot', '-', 'Fraction of roots in each soil layer', &
         patchout%froot, soilID, 'soil', xID, yID, zID, landID, patchID)

    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Production", &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out,NF90_GLOBAL,"Source", &
         'CABLE LSM output file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "CABLE_input_file", &
         TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Determine output aggregation details:
    IF(output%averaging(1:4) == 'user') THEN
       ! User-specified aggregation interval for output:
       READ(output%averaging(5:7), *) output%interval
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging", &
            TRIM(output%averaging(5:7))//'-hourly output')
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Convert interval value from hours to time steps (for use in output
       ! write):
       output%interval = output%interval * 3600/INT(dels)
    ELSE IF(output%averaging(1:3) == 'all') THEN ! output all timesteps
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging", &
            TRIM(output%averaging)//' timesteps recorded')
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Set output interval to be one time step
       output%interval = 1
    ELSE IF(output%averaging(1:2) == 'mo') THEN ! monthly output
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging", &
            TRIM(output%averaging))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Output interval will be determined dynamically for monthly output
    ELSE IF(output%averaging(1:2) == 'da') THEN ! daily output
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging", &
            TRIM(output%averaging))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Set output interval to be # time steps in 24 hours:
       output%interval = 3600*24/INT(dels)
    ELSE
       CALL cable_abort ('Unknown output averaging interval specified '// &
            'in namelist file. (SUBROUTINE open_output_file)')
    END IF

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating output file ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Write latitude and longitude variables:
    ok = NF90_PUT_VAR(ncid_out, latID, lat_all)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error writing latitude variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_VAR(ncid_out, lonID, lon_all)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error writing longitude variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    IF (output%grid == 'land' .OR. &
         (metGrid == 'land' .AND. output%grid == 'default'))  THEN

       ok = NF90_PUT_VAR(ncid_out, llatvID, toreal4(latitude) )
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error writing loc lat variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

       ok = NF90_PUT_VAR(ncid_out, llonvID, toreal4(longitude) )
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok, 'Error writing loc lon variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ENDIF

    ! Write GrADS coordinate variables
    ok = NF90_PUT_VAR(ncid_out, xvID, lon_all(:, 1))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing GrADS x coordinate variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_VAR(ncid_out, yvID, lat_all(1, :))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing GrADS y coordinate variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Write model parameters if requested:
    IF(output%params .OR. output%iveg) CALL write_ovar(ncid_out, opid%iveg, &
         'iveg', toreal4(veg%iveg), ranges%iveg, patchout%iveg, 'integer')

    IF (.not.cable_user%POPLUC) THEN
       IF ((output%params .OR. output%patchfrac) .AND. (patchout%patchfrac .OR. output%patch)) &
            CALL write_ovar(ncid_out, opid%patchfrac, 'patchfrac', &
            toreal4(patch(:)%frac), (/0.0, 1.0/), patchout%patchfrac, 'real')
    ENDIF
    IF(output%params .OR. output%isoil) CALL write_ovar(ncid_out, opid%isoil, &
         'isoil', toreal4(soil%isoilm), ranges%isoil, patchout%isoil,'integer')
    IF(output%params .OR. output%bch) CALL write_ovar(ncid_out, opid%bch, &
         'bch', toreal4(soil%bch), ranges%bch, patchout%bch, 'real')
    IF(output%params .OR. output%clay) CALL write_ovar(ncid_out, opid%clay, &
         'clay', toreal4(soil%clay), ranges%clay, patchout%clay, 'real')
    IF(output%params .OR. output%sand) CALL write_ovar(ncid_out, opid%sand, &
         'sand', toreal4(soil%sand), ranges%sand, patchout%sand, 'real')
    IF(output%params .OR. output%silt) CALL write_ovar(ncid_out, opid%silt, &
         'silt', toreal4(soil%silt), ranges%silt, patchout%silt, 'real')
    IF(output%params .OR. output%css) CALL write_ovar(ncid_out, opid%css, &
         'css', toreal4(soil%css), ranges%css, patchout%css, 'real')
    IF(output%params .OR. output%rhosoil) CALL write_ovar(ncid_out, &
         opid%rhosoil, 'rhosoil',toreal4(soil%rhosoil), &
         ranges%rhosoil, patchout%rhosoil, 'real')
    IF(output%params .OR. output%hyds) CALL write_ovar(ncid_out, opid%hyds, &
         'hyds', toreal4(soil%hyds), ranges%hyds, patchout%hyds, 'real')
    IF(output%params .OR. output%sucs) CALL write_ovar(ncid_out, opid%sucs, &
         'sucs', toreal4(soil%sucs), ranges%sucs, patchout%sucs, 'real')
    IF(output%params .OR. output%rs20) CALL write_ovar(ncid_out, opid%rs20, &
         'rs20', toreal4(veg%rs20), ranges%rs20, patchout%rs20, 'real')
    !         'rs20',toreal4(soil%rs20),ranges%rs20,patchout%rs20,'real')
    IF(output%params .OR. output%ssat) CALL write_ovar(ncid_out, opid%ssat, &
         'ssat', toreal4(soil%ssat), ranges%ssat, patchout%ssat, 'real')
    IF(output%params .OR. output%sfc) CALL write_ovar(ncid_out, opid%sfc, &
         'sfc', toreal4(soil%sfc), ranges%sfc, patchout%sfc, 'real')
    IF(output%params .OR. output%swilt) CALL write_ovar(ncid_out, opid%swilt, &
         'swilt', toreal4(soil%swilt), ranges%swilt, patchout%swilt, 'real')
    IF(output%params .AND. output%albsoil) CALL write_ovar(ncid_out, &
         opid%albsoil, 'albsoil', toreal4(soil%albsoil), &
         ranges%albsoil, patchout%albsoil, 'radiation')
    IF(output%params .OR. output%canst1) CALL write_ovar(ncid_out, &
         opid%canst1, 'canst1', toreal4(veg%canst1), &
         ranges%canst1, patchout%canst1, 'real')
    IF(output%params .OR. output%dleaf) CALL write_ovar(ncid_out, opid%dleaf, &
         'dleaf', toreal4(veg%dleaf), ranges%dleaf, patchout%dleaf, 'real')
    !    IF(output%params .OR. output%ejmax) CALL write_ovar(ncid_out, opid%ejmax, &
    !              'ejmax', toreal4(veg%ejmax), ranges%ejmax, patchout%ejmax, 'real')
    !Alexis
    !IF(output%params .OR. output%vcmax) CALL write_ovar(ncid_out, opid%vcmax, &
    !          'vcmax', toreal4(veg%vcmax), ranges%vcmax, patchout%vcmax, 'real')
    IF(output%params .OR. output%frac4) CALL write_ovar(ncid_out, opid%frac4, &
         'frac4', toreal4(veg%frac4), ranges%frac4, patchout%frac4, 'real')

    IF (.not.cable_user%CALL_POP) THEN
       IF(output%params .OR. output%hc) CALL write_ovar(ncid_out, opid%hc, &
            'hc', toreal4(veg%hc), ranges%hc, patchout%hc, 'real')
    ENDIF
    IF(output%params .OR. output%rp20) CALL write_ovar(ncid_out, opid%rp20, &
         'rp20', toreal4(veg%rp20),ranges%rp20, patchout%rp20, 'real')
    ! Ticket #56
    IF(output%params .OR. output%g0) CALL write_ovar(ncid_out, opid%g0, &
         'g0', toreal4(veg%g0),ranges%g0, patchout%g0, 'real')
    IF(output%params .OR. output%g1) CALL write_ovar(ncid_out, opid%g1, &
         'g1', toreal4(veg%g1),ranges%g1, patchout%g1, 'real')
    ! End Ticket #56
    IF(output%params .OR. output%rpcoef) CALL write_ovar(ncid_out, &
         opid%rpcoef, 'rpcoef', toreal4(veg%rpcoef), &
         ranges%rpcoef, patchout%rpcoef, 'real')
    IF(output%params .OR. output%shelrb) CALL write_ovar(ncid_out, &
         opid%shelrb, 'shelrb', toreal4(veg%shelrb), &
         ranges%shelrb, patchout%shelrb, 'real')
    IF(output%params .OR. output%xfang) CALL write_ovar(ncid_out, opid%xfang, &
         'xfang', toreal4(veg%xfang), ranges%xfang, patchout%xfang, 'real')
    IF(output%params .OR. output%wai) CALL write_ovar(ncid_out, opid%wai, &
         'wai', toreal4(veg%wai), ranges%wai, patchout%wai, 'real')
    IF(output%params .OR. output%vegcf) CALL write_ovar(ncid_out, opid%vegcf, &
         'vegcf', toreal4(veg%vegcf), ranges%vegcf, patchout%vegcf, 'real')
    IF(output%params .OR. output%extkn) CALL write_ovar(ncid_out, opid%extkn, &
         'extkn', toreal4(veg%extkn), ranges%extkn, patchout%extkn, 'real')
    IF(output%params .OR. output%tminvj) CALL write_ovar(ncid_out, &
         opid%tminvj, 'tminvj', toreal4(veg%tminvj), &
         ranges%tminvj, patchout%tminvj, 'real')
    IF(output%params .OR. output%tmaxvj) CALL write_ovar(ncid_out, &
         opid%tmaxvj, 'tmaxvj', toreal4(veg%tmaxvj), &
         ranges%tmaxvj, patchout%tmaxvj, 'real')
    IF(output%params .OR. output%vbeta) CALL write_ovar(ncid_out, opid%vbeta, &
         'vbeta', toreal4(veg%vbeta), ranges%vbeta, patchout%vbeta, 'real')
    IF(output%params .OR. output%xalbnir) CALL write_ovar(ncid_out, &
         opid%xalbnir, 'xalbnir', toreal4(veg%xalbnir), &
         ranges%xalbnir, patchout%xalbnir, 'real')
    IF(output%params .OR. output%meth) CALL write_ovar(ncid_out, opid%meth, &
         'meth', toreal4(veg%meth), ranges%meth, patchout%meth, 'integer')
    IF(output%params .OR. output%za) THEN
       CALL write_ovar(ncid_out, opid%za_uv, &
            'za_uv', toreal4(rough%za_uv), ranges%za, patchout%za, 'real')
       CALL write_ovar(ncid_out, opid%za_tq, &
            'za_tq', toreal4(rough%za_tq), ranges%za, patchout%za, 'real')
    ENDIF
    IF(output%params .AND. output%ratecp) CALL write_ovar(ncid_out, &
         opid%ratecp, 'ratecp',SPREAD(toreal4(bgc%ratecp),1,mp), ranges%ratecp, &
         patchout%ratecp,'plantcarbon')! no spatial dim at present
    IF(output%params .AND. output%ratecs) CALL write_ovar(ncid_out, &
         opid%ratecs, 'ratecs', SPREAD(toreal4(bgc%ratecs), 1, mp), ranges%ratecs, &
         patchout%ratecs, 'soilcarbon')! no spatial dim at present
    IF(output%params .AND. output%froot) CALL write_ovar (ncid_out, opid%froot, &
         'froot', toreal4(veg%froot), ranges%froot, patchout%froot, 'soil')
    IF(output%params .AND. output%zse) CALL write_ovar(ncid_out, opid%zse, &
         'zse', SPREAD(toreal4(soil%zse), 1, mp),ranges%zse, &
         patchout%zse, 'soil')! no spatial dim at present

  END SUBROUTINE open_output_file

  !=============================================================================

  SUBROUTINE write_output(dels, ktau, met, canopy, casaflux, casapool, casamet, ssnow, &
       rad, bal, air, soil, veg, SBOLTZ, EMLEAF, EMSOIL, c13o2pools, c13o2flux)
    ! Writes model output variables and, if requested, calls
    ! energy and mass balance routines. This subroutine is called
    ! each timestep, but may only write to the output file periodically,
    ! depending on whether the user has specified that output should be
    ! aggregated, e.g. to monthly or 6-hourly averages.
    REAL, INTENT(IN)              :: dels ! time step size
    INTEGER, INTENT(IN)           :: ktau ! timestep number in loop which include spinup
    REAL, INTENT(IN) :: SBOLTZ, EMLEAF, EMSOIL
    TYPE(met_type), INTENT(IN)         :: met  ! met data
    TYPE(canopy_type), INTENT(IN)      :: canopy ! canopy variable data
    TYPE(soil_snow_type), INTENT(IN)   :: ssnow ! soil data
    TYPE(soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE(radiation_type), INTENT(IN)  :: rad   ! radiation data
    TYPE(air_type), INTENT(IN)        :: air
    TYPE(veg_parameter_type), INTENT(IN) :: veg ! vegetation parameters
    TYPE(casa_flux), INTENT(IN) :: casaflux ! casa fluxes
    TYPE(casa_pool), INTENT(IN) :: casapool ! casa fluxes
    TYPE(balances_type), INTENT(INOUT) :: bal
    TYPE (casa_met), INTENT(IN) :: casamet
    ! 13C
    TYPE(c13o2_pool), INTENT(IN) :: c13o2pools ! 13CO2 pools
    TYPE(c13o2_flux), INTENT(IN) :: c13o2flux  ! 13CO2 fluxes

    REAL(r_2), DIMENSION(1) :: timetemp ! temporary variable for storing time
    ! value
    LOGICAL :: writenow ! write to output file during this time step?
    INTEGER, SAVE :: out_timestep ! counter for output time steps
    INTEGER, SAVE :: out_month ! counter for output month
    ! INTEGER, DIMENSION(mp) :: realyear ! fix problem for yr b4 leap yr
    INTEGER :: backtrack  ! modify timetemp for averaged output

    INTEGER :: dday ! number of past-years days for monthly output LN
    INTEGER :: iy   ! Counter
    INTEGER, SAVE   :: YearStart
    INTEGER :: ok ! for netcdf sync
    real(kind=4) :: rinterval
    real(r_2)    :: r2interval, gd2umols
    real(kind=4), parameter :: zero4 = real(0.0,4)
    real, dimension(mp) :: totlai

    ! logical :: opened
    ! integer :: varid
    gd2umols = 1.0_r_2 / (86400.0_r_2 * 1.201e-5_r_2)

    ! IF asked to check mass/water balance:
    IF (check%mass_bal) CALL mass_balance(dels, ktau, ssnow, soil, canopy, met, air, bal)

    ! IF asked to check energy balance:
    IF (check%energy_bal) CALL energy_balance(met, rad, canopy, bal, ssnow, SBOLTZ, EMLEAF, EMSOIL)

    ! Initialise output time step counter and month counter:
    IF (ktau == 1) THEN
       out_timestep = 0
       out_month = 0
       ! use met%year(1) instead of CABLE_USER%YearStart for non-GSWP forcing and leap years
       IF ( TRIM(cable_user%MetType) .EQ. '' ) then
          YearStart = met%year(1)
       ELSE
          YearStart = CABLE_USER%YearStart
       ENDIF
    END IF
    ! Decide on output averaging regime:
    IF(output%averaging(1:3) == 'all') THEN ! write every time step to file
       ! Set flag to write data for current time step:
       writenow = .TRUE.
       ! Set output time step to be current model time step:
       out_timestep = ktau
       backtrack = 0
    ELSE IF (output%averaging(1:4) == 'user' .OR. output%averaging(1:2)=='da') THEN
       ! user defined output interval or daily output
       IF (MOD(ktau, output%interval) == 0) THEN ! i.e.ktau divisible by
          ! interval
          ! write to output file this time step
          writenow = .TRUE.
          ! increment output time step counter:
          out_timestep = out_timestep + 1
          backtrack = output%interval / 2
       ELSE
          writenow = .FALSE.
       END IF
    ELSE IF(output%averaging(1:2) == 'mo') THEN ! write monthly averages to file
       ! ! realyear = met%year
       ! realyear = CurYear
       ! IF(ktau >= 365*24*3600/INT(dels)) THEN
       !    WHERE(met%doy == 1) realyear = realyear - 1   ! last timestep of year
       ! END IF

       ! LN Inserted for multiyear output
       dday = 0
       DO iy=YearStart, CurYear-1
          IF (IS_LEAPYEAR(iy) .AND. leaps) THEN
             dday = dday + 366
          ELSE
             dday = dday + 365
          ENDIF
       END DO
       ! LN Inserted for multiyear output

       ! Are we using leap year calendar?
       IF (leaps) THEN
          ! If currently a leap year:
          if (is_leapyear(CurYear)) then
             !! vh_js !!
             IF(ANY(NINT(real(lastdayl+dday) * 24. * 3600. / dels) == ktau)) THEN
                out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
                ! write to output file this time step
                writenow = .TRUE.
                ! increment output time step counter:
                out_timestep = out_timestep + 1
                ! set numbr of time steps in output period
                output%interval = daysml(out_month) * 24 * 3600 / INT(dels)
             ELSE
                writenow = .FALSE.
             END IF
          ELSE ! not currently a leap year
             ! last time step of month
             !! vh_js !!
             IF(ANY(NINT(real(lastday+dday) * 24. * 3600. / dels) == ktau)) THEN
                ! increment output month counter
                out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
                ! write to output file this time step
                writenow = .TRUE.
                ! increment output time step counter:
                out_timestep = out_timestep + 1
                ! set numbr of time steps in output period
                output%interval = daysm(out_month) * 24 * 3600 / INT(dels)
             ELSE
                writenow = .FALSE.
             END IF
          END IF
       ELSE ! not using leap year timing in this run
          !! vh_js !!
          IF(ANY(NINT(real((lastday+dday))*24.*3600./real(INT(dels)))==ktau)) THEN ! last time step of month
             ! increment output month counter
             out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
             ! write to output file this time step
             writenow = .TRUE.
             ! increment output time step counter:
             out_timestep = out_timestep + 1
             ! set number of time steps in output period
             output%interval = daysm(out_month) * 24 * 3600 / INT(dels)
          ELSE
             writenow = .FALSE.
          END IF
       END IF ! using leap year timing or not
       backtrack = output%interval / 2

    ELSE ! type of output aggregation
       CALL cable_abort('Unknown output averaging request in namelist file.'// &
            '(SUBROUTINE write_output)')
    END IF

    ! Note that size of averaging interval, output%interval, is set when opening
    ! output file unless output is monthly (in which case it's set above)

    ! If this time step is an output time step:
    IF (writenow) THEN
       ! Write to temporary time variable:
       timetemp(1) = REAL(REAL(ktau-backtrack)*dels,r_2)
       ! inquire(unit=ncid_out, opened=opened)
       ! if (.not. opened) ok = NF90_OPEN(filename%out, NF90_WRITE, ncid_out)
       ! ok = NF90_INQ_VARID(ncid_out, 'time', varid)
       ! Write time variable for this output time step:
       ok = NF90_PUT_VAR(ncid_out, ovid%tvar, timetemp, &
            start = (/out_timestep/), count = (/1/))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
            'Error writing time variable to ' &
            //TRIM(filename%out)// '(SUBROUTINE write_output)')
       rinterval  = toreal4(1) / toreal4(output%interval)
       r2interval = 1.0_r_2 / real(output%interval,r_2)

    END IF

    ! Arguments to write_ovar: current time step; output file netcdf file ID;
    ! netcdf variable ID; variable name; variable data; variable ranges;
    ! non-land fill value; include patch info for this var; any specific
    ! formatting info; met variables for reporting in case of abort.



    ! SWdown:  downward short-wave radiation [W/m^2]
    IF(output%met .OR. output%SWdown) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SWdown = out%SWdown + toreal4(met%fsd(:, 1) + met%fsd(:, 2))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SWdown = out%SWdown * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SWdown, 'SWdown', &
               out%SWdown, ranges%SWdown, patchout%SWdown, 'default', met)
          ! Reset temporary output variable:
          out%SWdown = zero4
       END IF
    END IF
    ! FracDiff: Diffuse Radiation Fraction 
    IF(output%met .OR. output%FracDiff) THEN
      ! Add current timestep's value to total of temporary output variable:
      out%FracDiff = out%FracDiff + toreal4(met%fdiff)
      IF(writenow) THEN
        ! Divide accumulated variable by number of accumulated time steps:
        out%FracDiff = out%FracDiff * rinterval
        ! Write value to file:
        CALL write_ovar(out_timestep, ncid_out, ovid%FracDiff, 'FracDiff', &
             out%FracDiff, ranges%FracDiff, patchout%FracDiff, 'default', met)
        ! Reset temporary output variable:
        out%FracDiff = zero4
      END IF
    END IF
    ! LWdown: downward long-wave radiation [W/m^2]
    IF(output%met .OR. output%LWdown) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LWdown = out%LWdown + toreal4(met%fld)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LWdown = out%LWdown * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LWdown, 'LWdown', &
               out%LWdown, ranges%LWdown, patchout%LWdown, 'default', met)
          ! Reset temporary output variable:
          out%LWdown = zero4
       END IF
    END IF
    ! Tair: surface air temperature [K]
    IF(output%met .OR. output%Tair) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Tair = out%Tair + toreal4(met%tk)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Tair = out%Tair * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Tair, 'Tair', out%Tair, &
               ranges%Tair, patchout%Tair, 'ALMA', met)
          ! Reset temporary output variable:
          out%Tair = zero4
       END IF
    END IF
    ! Rainf: rainfall [kg/m^2/s]
    IF(output%met .OR. output%Rainf) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Rainf = out%Rainf + toreal4(met%precip / dels)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Rainf = out%Rainf * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Rainf, 'Rainf', &
               out%Rainf, ranges%Rainf, patchout%Rainf, 'default', met)
          ! Reset temporary output variable:
          out%Rainf = zero4
       END IF
    END IF
    ! Snowf: snowfall [kg/m^2/s]
    IF(output%met .OR. output%Snowf) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Snowf = out%Snowf + toreal4(met%precip_sn / dels)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Snowf = out%Snowf * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Snowf, 'Snowf', &
               out%Snowf, ranges%Snowf, patchout%Snowf, 'default', met)
          ! Reset temporary output variable:
          out%Snowf = zero4
       END IF
    END IF
    ! PSurf: surface pressure [Pa]
    IF(output%met .OR. output%PSurf) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%PSurf = out%PSurf + toreal4(met%pmb)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PSurf = out%PSurf * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PSurf, 'PSurf', &
               out%PSurf, ranges%PSurf, patchout%PSurf, 'default', met)
          ! Reset temporary output variable:
          out%PSurf = zero4
       END IF
    END IF
    ! Qair: specific humidity [kg/kg]
    IF(output%met .OR. output%Qair) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qair = out%Qair + toreal4(met%qv)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qair = out%Qair * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qair, 'Qair', out%Qair, &
               ranges%Qair, patchout%Qair, 'ALMA', met)
          ! Reset temporary output variable:
          out%Qair = zero4
       END IF
    END IF
    ! Wind: windspeed [m/s]
    IF(output%met .OR. output%Wind) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Wind = out%Wind + toreal4(met%ua)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Wind = out%Wind * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Wind, 'Wind', out%Wind, &
               ranges%Wind, patchout%Wind, 'ALMA', met)
          ! Reset temporary output variable:
          out%Wind = zero4
       END IF
    END IF
    ! CO2air: CO2 concentration [ppmv]
    IF(output%met .OR. output%CO2air) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%CO2air = out%CO2air + toreal4(met%ca * 1000000.0)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%CO2air = out%CO2air * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%CO2air, 'CO2air', &
               OUT%CO2air, ranges%CO2air, patchout%CO2air, 'ALMA', met)
          ! Reset temporary output variable:
          out%CO2air = zero4
       END IF
    END IF
    !-----------------------WRITE FLUX DATA-------------------------------------
    ! Qle: latent heat flux [W/m^2]
    IF(output%flux .OR. output%Qle) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qle = out%Qle + toreal4(canopy%fe)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qle = out%Qle * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qle, 'Qle', out%Qle, &
               ranges%Qle, patchout%Qle, 'default', met)
          ! Reset temporary output variable:
          out%Qle = zero4
       END IF
    END IF
    ! Qh: sensible heat flux [W/m^2]
    IF(output%flux .OR. output%Qh) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qh = out%Qh + toreal4(canopy%fh)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qh = out%Qh * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qh, 'Qh', out%Qh, &
               ranges%Qh, patchout%Qh, 'default', met)
          ! Reset temporary output variable:
          out%Qh = zero4
       END IF
    END IF
    ! Qg: ground heat flux [W/m^2]
    IF(output%flux .OR. output%Qg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qg = out%Qg + toreal4(canopy%ga)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qg = out%Qg * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qg, 'Qg', out%Qg, &
               ranges%Qg, patchout%Qg, 'default', met)
          ! Reset temporary output variable:
          out%Qg = zero4
       END IF
    END IF
    ! Qs: surface runoff [kg/m^2/s]
    IF(output%flux .OR. output%Qs) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qs = out%Qs + ssnow%rnof1 / dels
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qs = out%Qs * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qs, 'Qs', out%Qs, &
               ranges%Qs, patchout%Qs, 'default', met)
          ! Reset temporary output variable:
          out%Qs = zero4
       END IF
    END IF
    ! Qsb: subsurface runoff [kg/m^2/s]
    IF(output%flux .OR. output%Qsb) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qsb = out%Qsb + toreal4(ssnow%rnof2 / dels)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qsb = out%Qsb * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qsb, 'Qsb', out%Qsb, &
               ranges%Qsb, patchout%Qsb, 'default', met)
          ! Reset temporary output variable:
          out%Qsb = zero4
       END IF
    END IF
    ! Evap: total evapotranspiration [kg/m^2/s]
    IF(output%flux .OR. output%Evap) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Evap = out%Evap + toreal4(canopy%fe / air%rlam)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Evap = out%Evap * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Evap, 'Evap', out%Evap, &
               ranges%Evap, patchout%Evap, 'default', met)
          ! Reset temporary output variable:
          out%Evap = zero4
       END IF
    END IF
    ! ECanop: interception evaporation [kg/m^2/s]
    IF(output%flux .OR. output%ECanop) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%ECanop = out%ECanop + toreal4(canopy%fevw / air%rlam)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%ECanop = out%ECanop * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%ECanop, 'ECanop', &
               out%ECanop, ranges%ECanop, patchout%ECanop, 'default', met)
          ! Reset temporary output variable:
          out%ECanop = zero4
       END IF
    END IF
    IF(output%flux) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%gsw_sl = out%gsw_sl + toreal4(canopy%gswx(:,1))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%gsw_sl = out%gsw_sl * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%gsw_sl, 'gsw_sl', out%gsw_sl, &
               ranges%gsw_sl, patchout%Tveg, 'default', met)
          ! Reset temporary output variable:
          out%gsw_sl = zero4
       END IF
       ! Add current timestep's value to total of temporary output variable:
       out%gsw_sh = out%gsw_sh + toreal4(canopy%gswx(:,2))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%gsw_sh = out%gsw_sh * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%gsw_sh, 'gsw_sh', out%gsw_sh, &
               ranges%gsw_sh, patchout%Tveg, 'default', met)
          ! Reset temporary output variable:
          out%gsw_sh = zero4
       END IF

    END IF


    ! TVeg: vegetation transpiration [kg/m^2/s]
    IF(output%flux .OR. output%TVeg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%TVeg = out%TVeg + toreal4(canopy%fevc / air%rlam)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%TVeg = out%TVeg * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%TVeg, 'TVeg', out%TVeg, &
               ranges%TVeg, patchout%TVeg, 'default', met)
          ! Reset temporary output variable:
          out%TVeg = zero4
       END IF
    END IF

    ! ESoil: bare soil evaporation [kg/m^2/s]
    IF(output%flux .OR. output%ESoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       IF(cable_user%SOIL_STRUC=='sli') THEN
          out%ESoil = out%ESoil + toreal4(ssnow%evap/dels) !vh!
       ELSE
          out%ESoil = out%ESoil + toreal4(canopy%fes / air%rlam)
       ENDIF
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%ESoil = out%ESoil * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%ESoil, 'ESoil', &
               out%ESoil, ranges%ESoil, patchout%ESoil, 'default', met)
          ! Reset temporary output variable:
          out%ESoil = zero4
       END IF
    END IF
    ! HVeg: sensible heat from vegetation [W/m^2]
    IF(output%flux .OR. output%HVeg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%HVeg = out%HVeg + toreal4(canopy%fhv)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%HVeg = out%HVeg * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%HVeg, 'HVeg', out%HVeg, &
               ranges%HVeg, patchout%HVeg, 'default', met)
          ! Reset temporary output variable:
          out%HVeg = zero4
       END IF
    END IF
    ! HSoil: sensible heat from soil [W/m^2]
    IF(output%flux .OR. output%HSoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%HSoil = out%HSoil + toreal4(canopy%fhs)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%HSoil = out%HSoil * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%HSoil, 'HSoil', &
               out%HSoil, ranges%HSoil, patchout%HSoil, 'default', met)
          ! Reset temporary output variable:
          out%HSoil = zero4
       END IF

       out%RnetSoil = out%RnetSoil + toreal4(canopy%fns)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%RnetSoil = out%RnetSoil * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%RnetSoil, 'RnetSoil', &
               out%RnetSoil, ranges%HSoil, patchout%HSoil, 'default', met)
          ! Reset temporary output variable:
          out%RnetSoil = zero4
       END IF
    END IF
    ! NEE: net ecosystem exchange [umol/m^2/s]
    IF(output%flux .OR. output%carbon .OR. output%NEE) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%NEE = out%NEE + toreal4(canopy%fnee / 1.201E-5)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%NEE = out%NEE * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%NEE, 'NEE', out%NEE, &
               ranges%NEE, patchout%NEE, 'default', met)
          ! Reset temporary output variable:
          out%NEE = zero4
       END IF
    END IF

    ! Qcan: absorbed radiation by sunlit canopy [W/m^2]
    IF(output%flux) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%qcan_sl = out%qcan_sl + toreal4(rad%qcan(:,1,:))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%qcan_sl = out%qcan_sl * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%qcan_sl, 'qcan_sl', out%qcan_sl, &
               ranges%Qcan, patchout%Qcan, 'radiation', met)
          ! Reset temporary output variable:
          out%qcan_sl = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%qcan_sh = out%qcan_sh + toreal4(rad%qcan(:,2,:))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%qcan_sh = out%qcan_sh * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%qcan_sh, 'qcan_sh', out%qcan_sh, &
               ranges%Qcan, patchout%Qcan, 'radiation', met)
          ! Reset temporary output variable:
          out%qcan_sh = zero4
       END IF
    END IF

    !-----------------------WRITE SOIL STATE DATA-------------------------------
    ! SoilMoist: av.layer soil moisture [kg/m^2]
    IF(output%soil .OR. output%SoilMoist) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SoilMoist = out%SoilMoist + toreal4(ssnow%wb)
       out%SoilMoistPFT = out%SoilMoistPFT + toreal4(sum(ssnow%wb * 1000.0_r_2 * real(spread(soil%zse,1,mp),r_2),2))
       out%SoilMoistIce = out%SoilMoistIce + toreal4(ssnow%wbice)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilMoist = out%SoilMoist * rinterval
          out%SoilMoistPFT = out%SoilMoistPFT * rinterval
          out%SoilMoistIce = out%SoilMoistIce * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilMoist, 'SoilMoist', &
               out%SoilMoist, ranges%SoilMoist, patchout%SoilMoist, 'soil', met)
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilMoistPFT, 'SoilMoistPFT', &
               out%SoilMoistPFT, ranges%SoilMoistPFT, patchout%SoilMoist, 'default', met)
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilMoistIce, 'SoilMoistIce', &
               out%SoilMoistIce, ranges%SoilMoist, patchout%SoilMoistIce, 'soil', met)
          ! Reset temporary output variable:
          out%SoilMoist = zero4
          out%SoilMoistPFT = zero4
          out%SoilMoistIce = zero4
       END IF
    END IF
    ! SoilTemp: av.layer soil temperature [K]
    IF(output%soil .OR. output%SoilTemp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SoilTemp = out%SoilTemp + toreal4(ssnow%tgg)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilTemp = out%SoilTemp * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilTemp, 'SoilTemp', &
               out%SoilTemp, ranges%SoilTemp, patchout%SoilTemp, 'soil', met)
          ! Reset temporary output variable:
          out%SoilTemp = zero4
       END IF
    END IF
    ! BaresoilT: surface bare soil temp [K]
    IF (output%soil .OR. output%BaresoilT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%BaresoilT = out%BaresoilT + toreal4(ssnow%tgg(:, 1))
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%BaresoilT = out%BaresoilT * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%BaresoilT, 'BaresoilT', &
               out%BaresoilT, ranges%BaresoilT, patchout%BaresoilT, 'default', met)
          ! Reset temporary output variable:
          out%BaresoilT = zero4
       END IF
    END IF
    !----------------------WRITE SNOW STATE DATA--------------------------------
    ! SWE: snow water equivalent [kg/m^2]
    IF(output%snow .OR. output%SWE) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SWE = out%SWE + toreal4(ssnow%snowd)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SWE = out%SWE * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SWE, 'SWE', out%SWE, &
               ranges%SWE, patchout%SWE, 'default', met)
          ! Reset temporary output variable:
          out%SWE = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%SnowMelt = out%SnowMelt + toreal4(ssnow%smelt)/dels
       ! temp test vh !
       !out%SnowMelt = out%SnowMelt + toreal4(ssnow%nsteps)/dels
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SnowMelt = out%SnowMelt * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SnowMelt, 'SnowMelt', out%SnowMelt, &
               (/-99999.0, 9999999.0/), patchout%SnowMelt, 'default', met)
          ! Reset temporary output variable:
          out%SnowMelt = zero4
       END IF

    END IF
    ! SnowT: snow surface temp [K]
    IF (output%snow .OR. output%SnowT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SnowT = out%SnowT + toreal4(ssnow%tggsn(:, 1))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SnowT = out%SnowT * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SnowT, 'SnowT', &
               out%SnowT, ranges%SnowT, patchout%SnowT, 'default', met)
          ! Reset temporary output variable:
          out%SnowT = zero4
       END IF
    END IF
    ! SnowDepth: actual depth of snow in [m]
    IF (output%snow .OR. output%SnowDepth) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SnowDepth = out%SnowDepth + toreal4(SUM(ssnow%sdepth, 2))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SnowDepth = out%SnowDepth * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SnowDepth, 'SnowDepth', &
               out%SnowDepth, ranges%SnowDepth, patchout%SnowDepth, 'default', met)
          ! Reset temporary output variable:
          out%SnowDepth = zero4
       END IF
    END IF
    !-------------------------WRITE RADIATION DATA------------------------------
    ! fraction of direct radiation [0-1]
    IF (output%radiation) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%fbeam = out%fbeam + toreal4(rad%fbeam(:,1))
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%fbeam = out%fbeam * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%fbeam, 'fbeam', &
               out%fbeam, ranges%Albedo, patchout%fbeam, 'default', met)
          ! Reset temporary output variable:
          out%fbeam = zero4
       END IF
    END IF

    ! SWnet: net shortwave [W/m^2]
    IF (output%radiation .OR. output%SWnet) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SWnet = out%SWnet + toreal4(SUM(rad%qcan(:, :, 1), 2) + &
            SUM(rad%qcan(:, :, 2), 2) + rad%qssabs)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SWnet = out%SWnet * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SWnet, 'SWnet', &
               out%SWnet, ranges%SWnet, patchout%SWnet, 'default', met)
          ! Reset temporary output variable:
          out%SWnet = zero4
       END IF
    END IF
    ! LWnet: net longwave [W/m^2]
    IF (output%radiation .OR. output%LWnet) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LWnet = out%LWnet + &
            toreal4(met%fld - sboltz * emleaf * canopy%tv ** 4 * (1 - rad%transd) - rad%flws * rad%transd)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LWnet = out%LWnet * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LWnet, 'LWnet', &
               out%LWnet, ranges%LWnet, patchout%LWnet, 'default', met)
          ! Reset temporary output variable:
          out%LWnet = zero4
       END IF
    END IF
    ! Rnet: net absorbed radiation [W/m^2]
    IF (output%radiation .OR. output%Rnet) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Rnet = out%Rnet + toreal4(met%fld - sboltz * emleaf * canopy%tv ** 4 * &
            (1 - rad%transd) -rad%flws * rad%transd + SUM(rad%qcan(:, :, 1), 2) + &
            SUM(rad%qcan(:, :, 2), 2) + rad%qssabs)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Rnet = out%Rnet * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Rnet, 'Rnet', out%Rnet, &
               ranges%Rnet, patchout%Rnet, 'default', met)
          ! Reset temporary output variable:
          out%Rnet = zero4
       END IF
    END IF
    ! Albedo:
    IF (output%radiation .OR. output%Albedo) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Albedo = out%Albedo + toreal4((rad%albedo(:, 1) + rad%albedo(:, 2)) * 0.5)
       !! output calc of soil albedo based on colour? - Ticket #27
       !IF (calcsoilalbedo) THEN

       IF ( mod(ktau,INT(24.0*3600.0/dels)) == INT(24.0*3600.0/dels)/2 &
            .OR. output%averaging(1:3) == 'all') THEN

          out%visAlbedo(:,1) = out%visAlbedo(:,1) + toreal4(rad%albedo(:, 1) )
          out%nirAlbedo(:,1) = out%nirAlbedo(:,1) + toreal4(rad%albedo(:, 2) )

          out%visAlbedo(:,2) = out%visAlbedo(:,2) + toreal4(rad%reffbm(:, 1) )
          out%nirAlbedo(:,2) = out%nirAlbedo(:,2) + toreal4(rad%reffbm(:, 2) )

          out%visAlbedo(:,3) = out%visAlbedo(:,3) + toreal4(rad%reffdf(:, 1) )
          out%nirAlbedo(:,3) = out%nirAlbedo(:,3) + toreal4(rad%reffdf(:, 2) )
       ENDIF
       !END IF

       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Albedo = out%Albedo * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Albedo, 'Albedo', &
               out%Albedo, ranges%Albedo, patchout%Albedo, 'default', met)
          ! Reset temporary output variable:
          out%Albedo = zero4

          ! output calc of soil albedo based on colour? - Ticket #27
          !IF (calcsoilalbedo) THEN
          IF (output%averaging(1:3) == 'all') THEN
             out%visAlbedo = out%visAlbedo * rinterval
          ELSE
             out%visAlbedo = out%visAlbedo * rinterval * &
                  INT(24.0*3600.0/dels)
          ENDIF
          CALL write_ovar(out_timestep, ncid_out, ovid%visAlbedo, 'visAlbedo',&
               out%visAlbedo, ranges%visAlbedo, patchout%visAlbedo, 'radiation', met)
          out%visAlbedo = zero4
          IF (output%averaging(1:3) == 'all') THEN
             out%nirAlbedo = out%nirAlbedo * rinterval
          ELSE
             out%nirAlbedo = out%nirAlbedo * rinterval* &
                  INT(24.0*3600.0/dels)
          ENDIF
          CALL write_ovar(out_timestep, ncid_out, ovid%nirAlbedo, 'nirAlbedo',&
               out%nirAlbedo, ranges%nirAlbedo, patchout%nirAlbedo, 'radiation', met)
          out%nirAlbedo = zero4
          !END IF
       END IF
    END IF
    ! RadT: Radiative surface temperature [K]
    IF (output%radiation .OR. output%RadT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%RadT = out%RadT + toreal4((((1.0 - rad%transd) * emleaf * sboltz * &
            canopy%tv ** 4 + rad%transd * emsoil * sboltz * (ssnow%tss) ** 4) / sboltz)**0.25)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%RadT = out%RadT * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%RadT, 'RadT', out%RadT, &
               ranges%RadT, patchout%RadT, 'default', met)
          ! Reset temporary output variable:
          out%RadT = zero4
       END IF
    END IF
    !------------------------WRITE VEGETATION DATA------------------------------
    ! VegT: vegetation temperature [K]
    IF(output%veg .OR. output%VegT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%VegT = out%VegT + toreal4(canopy%tv)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%VegT = out%VegT * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%VegT, 'VegT', out%VegT, &
               ranges%VegT, patchout%VegT, 'default', met)
          ! Reset temporary output variable:
          out%VegT = zero4
       END IF
    END IF
    ! CanT: within-canopy temperature [K]
    IF (output%veg .OR. output%CanT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%CanT = out%CanT + toreal4(met%tvair)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%CanT = out%CanT * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%CanT, 'CanT', out%CanT, &
               ranges%CanT, patchout%CanT, 'default', met)
          ! Reset temporary output variable:
          out%CanT = zero4
       END IF
    END IF
    IF (output%veg .OR. output%Fwsoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Fwsoil = out%Fwsoil + toreal4(canopy%fwsoil)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Fwsoil = out%Fwsoil * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Fwsoil, 'Fwsoil', out%Fwsoil, &
               ranges%Fwsoil, patchout%Fwsoil, 'default', met)
          ! Reset temporary output variable:
          out%Fwsoil = zero4
       END IF
    END IF
    ! CanopInt: total canopy water storage [kg/m^2]
    IF (output%veg .OR. output%CanopInt) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%CanopInt = out%CanopInt + toreal4(canopy%cansto)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%CanopInt = out%CanopInt * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%CanopInt, 'CanopInt', &
               out%CanopInt, ranges%CanopInt, patchout%CanopInt, 'default', met)
          ! Reset temporary output variable:
          out%CanopInt = zero4
       END IF
    END IF
    ! LAI:
    IF (output%veg .OR. output%LAI) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LAI = out%LAI + toreal4(veg%vlai)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LAI = out%LAI * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LAI, 'LAI', out%LAI, &
               ranges%LAI, patchout%LAI, 'default', met)
          ! Reset temporary output variable:
          out%LAI = zero4
       END IF
    END IF

    !Alexis
    IF(output%veg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%vcmax = out%vcmax + toreal4(veg%vcmax)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%vcmax = out%vcmax * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%vcmax, 'vcmax', out%vcmax, &
               ranges%vcmax, patchout%LAI, 'default', met)
          ! Reset temporary output variable:
          out%vcmax = zero4
       END IF

       out%ejmax = out%ejmax + toreal4(veg%ejmax)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%ejmax = out%ejmax * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%ejmax, 'jmax', out%ejmax, &
               ranges%ejmax, patchout%LAI, 'default', met)
          ! Reset temporary output variable:
          out%ejmax = zero4
       END IF

    END IF

    !------------------------WRITE BALANCES DATA--------------------------------
    ! Ebal: cumulative energy balance [W/m^2]
    IF (output%balances .OR. output%Ebal) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Ebal = out%Ebal + toreal4(bal%ebal_tot)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Ebal = out%Ebal * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Ebal, 'Ebal', out%Ebal, &
               ranges%Ebal, patchout%Ebal, 'default', met)
          ! Reset temporary output variable:
          out%Ebal = zero4
       END IF
    END IF
    ! Wbal: cumulative water balance  [kg/m^2/s]
    IF (output%balances .OR. output%Wbal) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Wbal = out%Wbal + toreal4(bal%wbal_tot)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Wbal = out%Wbal * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Wbal, 'Wbal', &
               out%Wbal, ranges%Wbal, patchout%Wbal, 'default', met)
          ! Reset temporary output variable:
          out%Wbal = zero4
       END IF
    END IF
    !------------------------WRITE CARBON DATA----------------------------------
    ! GPP: gross primary production C by veg [umol/m^2/s]
    !      added frday in the calculation of GPP (BP may08)
    IF (output%carbon .OR. output%GPP) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%GPP = out%GPP + toreal4((-1.0 * canopy%fpn + canopy%frday) / 1.201E-5)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%GPP = out%GPP * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%GPP, 'GPP', out%GPP, &
               ranges%GPP, patchout%GPP, 'default', met)
          ! Reset temporary output variable:
          out%GPP = zero4
       END IF
    END IF

    ! components of GPP
    IF (output%GPP_components) THEN
       out%scalex_sl   = out%scalex_sl   + toreal4(rad%scalex(:,1))
       out%scalex_sh   = out%scalex_sh   + toreal4(rad%scalex(:,2))
       out%dlf         = out%dlf         + toreal4(canopy%dlf)
       out%An_sl       = out%An_sl       + toreal4(canopy%A_sl*1.0e6_r_2)
       out%An_sh       = out%An_sh       + toreal4(canopy%A_sh*1.0e6_r_2)
       out%ci_sl       = out%ci_sl       + toreal4(canopy%ci(:,1)*1.0e6_r_2)
       out%ci_sh       = out%ci_sh       + toreal4(canopy%ci(:,2)*1.0e6_r_2)
       out%GPP_sl      = out%GPP_sl      + toreal4(canopy%GPP_sl*1.0e6_r_2)
       out%GPP_sh      = out%GPP_sh      + toreal4(canopy%GPP_sh*1.0e6_r_2)
       out%GPP_slC     = out%GPP_slC     + toreal4(canopy%A_slC*1.0e6_r_2)
       out%GPP_shC     = out%GPP_shC     + toreal4(canopy%A_shC*1.0e6_r_2)
       out%GPP_slJ     = out%GPP_slJ     + toreal4(canopy%A_slJ*1.0e6_r_2)
       out%GPP_shJ     = out%GPP_shJ     + toreal4(canopy%A_shJ*1.0e6_r_2)
       out%eta_GPP_cs  = out%eta_GPP_cs  + toreal4(canopy%eta_A_cs*1.0e6_r_2)
       out%eta_TVeg_cs = out%eta_TVeg_cs + toreal4(canopy%eta_fevc_cs/air%rlam)
       out%dGPPdcs     = out%dGPPdcs     + toreal4(canopy%dAdcs*1.0e6_r_2)
       out%CO2s        = out%CO2s        + toreal4(canopy%cs*1.0e6_r_2)
       totlai(:) = rad%fvlai(:,1) + rad%fvlai(:,2)
       where (totlai(:) .gt. 0.01)
          out%gsw_TVeg =  out%gsw_TVeg + &
               toreal4((canopy%gswx(:,1) * rad%fvlai(:,1)/totlai(:) + &
               canopy%gswx(:,2) * rad%fvlai(:,2)/totlai(:)) * &
               real(canopy%fevc) / air%rlam)
          out%vcmax_ts =  out%vcmax_ts + &
               toreal4(veg%vcmax_sun * rad%fvlai(:,1)/totlai(:) + &
               veg%vcmax_shade * rad%fvlai(:,2)/totlai(:))
          out%jmax_ts =  out%jmax_ts + &
               toreal4(veg%ejmax_sun * rad%fvlai(:,1)/totlai(:) + &
               veg%ejmax_shade * rad%fvlai(:,2)/totlai(:))
       elsewhere
          out%gsw_TVeg = out%gsw_TVeg
          out%vcmax_ts = out%vcmax_ts + toreal4(veg%vcmax_shade)
          out%jmax_ts  = out%jmax_ts  + toreal4(veg%ejmax_shade)
       endwhere

       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%scalex_sl   = out%scalex_sl   * rinterval
          out%scalex_sh   = out%scalex_sh   * rinterval
          out%dlf         = out%dlf         * rinterval
          out%An_sl       = out%An_sl       * rinterval
          out%An_sh       = out%An_sh       * rinterval
          out%ci_sl       = out%ci_sl       * rinterval
          out%ci_sh       = out%ci_sh       * rinterval
          out%GPP_sl      = out%GPP_sl      * rinterval
          out%GPP_sh      = out%GPP_sh      * rinterval
          out%GPP_slC     = out%GPP_slC     * rinterval
          out%GPP_shC     = out%GPP_shC     * rinterval
          out%GPP_slJ     = out%GPP_slJ     * rinterval
          out%GPP_shJ     = out%GPP_shJ     * rinterval
          out%eta_GPP_cs  = out%eta_GPP_cs  * rinterval
          out%eta_TVeg_cs = out%eta_TVeg_cs * rinterval
          out%gsw_TVeg    = out%gsw_TVeg    * rinterval
          out%dGPPdcs     = out%dGPPdcs     * rinterval
          out%CO2s        = out%CO2s        * rinterval
          out%vcmax_ts    = out%vcmax_ts    * rinterval
          out%jmax_ts     = out%jmax_ts     * rinterval

          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%scalex_sl, 'scalex_sl', out%scalex_sl, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%scalex_sh, 'scalex_sh', out%scalex_sh, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%dlf, 'leaf to air vpd', out%dlf, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%An_sl, 'Anet_sl', out%An_sl, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%An_sh, 'Anet_sh', out%An_sh, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%ci_sl, 'ci_sl', out%ci_sl, &
               ranges%CO2air, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%ci_sh, 'ci_sh', out%ci_sh, &
               ranges%CO2air, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%GPP_sl, 'GPP_sl', out%GPP_sl, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%GPP_sh, 'GPP_sh', out%GPP_sh, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%GPP_slC, 'GPP_slC', out%GPP_slC, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%GPP_shC, 'GPP_shC', out%GPP_shC, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%GPP_slJ, 'GPP_slJ', out%GPP_slJ, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%GPP_shJ, 'GPP_shJ', out%GPP_shJ, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%eta_GPP_cs, 'eta_GPP_cs', &
               out%eta_GPP_cs, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%eta_TVeg_cs, 'eta_TVeg_cs', &
               out%eta_TVeg_cs, &
               ranges%TVeg, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%gsw_TVeg, 'gsw_TVeg', &
               out%gsw_TVeg, &
               ranges%TVeg, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%dGPPdcs, 'dGPPdcs', &
               out%dGPPdcs, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%CO2s, 'CO2s', &
               out%CO2s, &
               ranges%GPP, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%vcmax_ts, 'LAI-weighted vcmax', &
               out%vcmax_ts, &
               ranges%vcmax, patchout%GPP, 'default', met)

          CALL write_ovar(out_timestep, ncid_out, ovid%jmax_ts, 'LAI-weighted jmax', &
               out%jmax_ts, &
               ranges%ejmax, patchout%GPP, 'default', met)

          ! Reset temporary output variable:
          out%scalex_sl   = zero4
          out%scalex_sh   = zero4
          out%dlf         = zero4
          out%An_sl       = zero4
          out%An_sh       = zero4
          out%ci_sl       = zero4
          out%ci_sh       = zero4
          out%GPP_sl      = zero4
          out%GPP_sh      = zero4
          out%GPP_slC     = zero4
          out%GPP_shC     = zero4
          out%GPP_slJ     = zero4
          out%GPP_shJ     = zero4
          out%eta_GPP_cs  = zero4
          out%eta_TVeg_cs = zero4
          out%gsw_TVeg    = zero4
          out%dGPPdcs     = zero4
          out%CO2s        = zero4
          out%vcmax_ts    = zero4
          out%jmax_ts     = zero4

       ENDIF

    ENDIF

    ! NPP: net primary production of C by veg [umol/m^2/s]
    IF (output%carbon .OR. output%NPP) THEN
       ! Add current timestep's value to total of temporary output variable:
       !out%NPP = out%NPP + REAL((-1.0 * canopy%fpn - canopy%frp &
       !     - casaflux%clabloss/86400.0) / 1.201E-5, 4)
       ! vh ! expression below can be slightly different form that above in cases where
       ! leaf maintenance respiration is reduced in CASA
       ! (relative to its original value calculated in cable_canopy)
       ! in order to avoid negative carbon stores.
       IF(output%casa) THEN
          out%NPP = out%NPP + toreal4(casaflux%cnpp * gd2umols)
       ELSE
          out%NPP = out%NPP + toreal4((-1.0 * canopy%fpn - canopy%frp &
               - casaflux%clabloss/86400.0) / 1.201E-5)
       ENDIF
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%NPP = out%NPP * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%NPP, 'NPP', out%NPP, &
               ranges%NPP, patchout%NPP, 'default', met)
          ! Reset temporary output variable:
          out%NPP = zero4
       END IF
    END IF
    ! AutoResp: autotrophic respiration [umol/m^2/s]
    IF(output%carbon .OR. output%AutoResp) THEN
       ! Add current timestep's value to total of temporary output variable:
       !out%AutoResp = out%AutoResp + REAL((canopy%frp + canopy%frday + casaflux%clabloss/86400.0) &
       !                                    / 1.201E-5, 4)
       ! vh ! expression below can be slightly different from that above in cases where
       ! leaf maintenance respiration is reduced in CASA
       ! (relative to its original value calculated in cable_canopy)
       ! in order to avoid negative carbon stores.
       IF (output%casa) THEN
          out%AutoResp = out%AutoResp + toreal4(canopy%frday / 1.201E-5) + &
               toreal4((casaflux%crmplant(:,2) + casaflux%crmplant(:,3) + &
               casaflux%crgplant + casaflux%clabloss) * gd2umols)
       ELSE
          out%AutoResp = out%AutoResp + toreal4((canopy%frp + canopy%frday) / 1.201e-5)
       ENDIF

       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%AutoResp = out%AutoResp * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%AutoResp, 'AutoResp', &
               out%AutoResp, ranges%AutoResp, patchout%AutoResp, 'default', met)
          ! Reset temporary output variable:
          out%AutoResp = zero4
       END IF

       IF (output%casa) THEN
          out%RootResp = out%RootResp + toreal4(casaflux%crmplant(:,3) * gd2umols) !+ &
          ! toreal4(0.3*casaflux%crmplant(:,2)/86400.0/ 1.201E-5)
          IF(writenow) THEN
             ! Divide accumulated variable by number of accumulated time steps:
             out%RootResp = out%RootResp * rinterval
             ! Write value to file:
             CALL write_ovar(out_timestep, ncid_out, ovid%RootResp, 'RootResp', &
                  out%RootResp, ranges%AutoResp, patchout%AutoResp, 'default', met)
             ! Reset temporary output variable:
             out%RootResp = zero4
          END IF
       END IF

       IF (output%casa) THEN
          out%StemResp = out%StemResp + toreal4(casaflux%crmplant(:,2) * gd2umols)
          IF(writenow) THEN
             ! Divide accumulated variable by number of accumulated time steps:
             out%StemResp = out%StemResp * rinterval
             ! Write value to file:
             CALL write_ovar(out_timestep, ncid_out, ovid%StemResp, 'StemResp', &
                  out%StemResp, ranges%AutoResp, patchout%AutoResp, 'default', met)
             ! Reset temporary output variable:
             out%StemResp = zero4
          END IF
       END IF

    END IF

    ! LeafResp: Leaf respiration [umol/m^2/s]
    IF(output%carbon .OR. output%LeafResp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LeafResp = out%LeafResp + toreal4(canopy%frday / 1.201e-5)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LeafResp = out%LeafResp * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LeafResp, 'LeafResp', &
               out%LeafResp, ranges%LeafResp, patchout%LeafResp, 'default', met)
          ! Reset temporary output variable:
          out%LeafResp = zero4
       END IF
    END IF
    ! HeteroResp: heterotrophic respiration [umol/m^2/s]
    IF(output%carbon .OR. output%HeteroResp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%HeteroResp = out%HeteroResp + toreal4(canopy%frs / 1.201E-5)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%HeteroResp = out%HeteroResp * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%HeteroResp, &
               'HeteroResp', out%HeteroResp, ranges%HeteroResp, &
               patchout%HeteroResp, 'default', met)
          ! Reset temporary output variable:
          out%HeteroResp = zero4
       END IF
    END IF

    ! output patch area
    IF(output%casa) THEN
       out%Area = toreal4(casamet%areacell / 1.0e6_r_2) ! km2
       IF(writenow) THEN
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Area, 'Area', out%Area, &
               ranges%Area, patchout%Area, 'default', met)
       END IF
       IF (cable_user%POPLUC) THEN
          ! output patch fraction
          IF(writenow) THEN
             CALL write_ovar(out_timestep, ncid_out, opid%patchfrac, 'patchfrac', &
                  toreal4(patch(:)%frac), (/0.0, 1.0/), patchout%patchfrac, 'default',met)

          END IF
       ENDIF
    ENDIF
    IF (cable_user%CALL_POP) THEN
       IF(writenow) THEN
          IF(output%params .OR. output%hc) CALL write_ovar(out_timestep,ncid_out, opid%hc, &
               'hc', toreal4(veg%hc), ranges%hc, patchout%hc, 'default', met)
       ENDIF
    ENDIF

    ! vh_mc ! additional variables for ESM-SnowMIP
    IF (output%snowmip) THEN
       ! Add current timestep's value to total of temporary output variable
       out%hfds       = out%hfds        + toreal4(canopy%ga)
       out%hfdsn      = out%hfdsn       + toreal4(merge(canopy%ga,0.,sum(ssnow%sdepth,2)>0.))
       out%hfls       = out%hfls        + toreal4(canopy%fe)
       out%hfmlt      = out%hfmlt       + toreal4(ssnow%E_fusion_sn) ! energy of fusion [W/m2]
       out%hfrs       = out%hfrs        + toreal4(ssnow%Qadv_rain_sn)
       out%hfsbl      = out%hfsbl       + toreal4(ssnow%E_SUBLIMATION_SN)
       out%hfss       = out%hfss        + toreal4(canopy%fh)
       out%rlus       = out%rlus        + toreal4(emsoil * sboltz * &
            (rad%transd*ssnow%otss**4 + (1.0-rad%transd)*canopy%tv**4))
       out%rsus       = out%rsus        + toreal4(rad%albedo(:,1)*met%fsd(:,1) + rad%albedo(:,2)*met%fsd(:,2))
       out%esn        = out%esn         + toreal4(ssnow%evap_liq_sn)
       out%evspsbl    = out%evspsbl     + toreal4(ssnow%evap/dels + max(real(canopy%fevc),0.0)/air%rlam)
       IF (cable_user%soil_struc=='sli') THEN
          out%evspsblsoi = out%evspsblsoi  + toreal4(ssnow%evap / dels)
       ELSE
          out%evspsblsoi = out%evspsblsoi  + toreal4(canopy%fes / real(air%rlam,r_2))
       ENDIF
       out%evspsblveg = out%evspsblveg  + toreal4(canopy%fevw / air%rlam)
       out%mrrob      = out%mrrob       + toreal4(ssnow%rnof2 / dels)
       out%mrros      = out%mrros       + toreal4(ssnow%rnof1 / dels)
       out%sbl        = out%sbl         + toreal4(ssnow%E_sublimation_sn)
       out%snm        = out%snm         + toreal4(ssnow%surface_melt)
       out%snmsl      = out%snmsl       + toreal4(ssnow%smelt)
       out%tran       = out%tran        + toreal4(canopy%fevc / real(air%rlam,r_2))
       out%albs       = out%albs        + toreal4((rad%albedo(:,1) + rad%albedo(:,2)) * 0.5)
       out%albsn      = out%albsn       + toreal4((ssnow%albsoilsn(:,1) + ssnow%albsoilsn(:,2)) * 0.5)
       out%cw         = out%cw          + toreal4(canopy%cansto)
       where (sum(ssnow%sdepth,2)>0.)
          out%lqsn    = out%lqsn        + toreal4(sum(ssnow%snowliq,2)/ssnow%snowd)
       elsewhere
          out%lqsn    = out%lqsn        + toreal4(0.)
       endwhere
       out%lwsnl      = out%lwsnl       + toreal4(sum(ssnow%snowliq,2))
       out%mrfsofr    = out%mrfsofr     + toreal4(ssnow%wbice/ssnow%wb)
       out%mrlqso     = out%mrlqso      + toreal4((ssnow%wb-ssnow%wbice)/ssnow%wb)
       out%mrlsl      = out%mrlsl       + toreal4(ssnow%wb*1000.0_r_2*real(spread(soil%zse,1,mp),r_2))
       out%snc        = out%snc         + toreal4(merge(1.,0.,sum(ssnow%sdepth,2)>0.))
       out%snd        = out%snd         + toreal4(sum(ssnow%sdepth, 2))
       out%snw        = out%snw         + toreal4(ssnow%snowd)
       out%snwc       = out%snwc        + toreal4(0.)
       out%tcs        = out%tcs         + toreal4(canopy%tv)
       out%tgs        = out%tgs         + toreal4(ssnow%tgg(:,1))
       out%ts         = out%ts          + toreal4(ssnow%tss)
       ! wrong dimension! Suggest using existing SoilTemp output variable
       out%tsl        = out%tsl         + toreal4(ssnow%tgg)
       ! out%tsl        = out%tsl         + toreal4(0)
       out%tsn        = out%tsn         + toreal4(ssnow%tggsn(:, 1))
       out%tsns       = out%tsns        + toreal4(ssnow%tss)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps
          out%hfds       = out%hfds        * rinterval
          out%hfdsn      = out%hfdsn       * rinterval
          out%hfls       = out%hfls        * rinterval
          out%hfmlt      = out%hfmlt       * rinterval
          out%hfrs       = out%hfrs        * rinterval
          out%hfsbl      = out%hfsbl       * rinterval
          out%hfss       = out%hfss        * rinterval
          out%rlus       = out%rlus        * rinterval
          out%rsus       = out%rsus        * rinterval
          out%esn        = out%esn         * rinterval
          out%evspsbl    = out%evspsbl     * rinterval
          out%evspsblsoi = out%evspsblsoi  * rinterval
          out%evspsblveg = out%evspsblveg  * rinterval
          out%mrrob      = out%mrrob       * rinterval
          out%mrros      = out%mrros       * rinterval
          out%sbl        = out%sbl         * rinterval
          out%snm        = out%snm         * rinterval
          out%snmsl      = out%snmsl       * rinterval
          out%tran       = out%tran        * rinterval
          out%albs       = out%albs        * rinterval
          out%albsn      = out%albsn       * rinterval
          out%cw         = out%cw          * rinterval
          out%lqsn       = out%lqsn        * rinterval
          out%lwsnl      = out%lwsnl       * rinterval
          out%mrfsofr    = out%mrfsofr     * rinterval
          out%mrlqso     = out%mrlqso      * rinterval
          out%mrlsl      = out%mrlsl       * rinterval
          out%snc        = out%snc         * rinterval
          out%snd        = out%snd         * rinterval
          out%snw        = out%snw         * rinterval
          out%snwc       = out%snwc        * rinterval
          out%tcs        = out%tcs         * rinterval
          out%tgs        = out%tgs         * rinterval
          out%ts         = out%ts          * rinterval
          out%tsl        = out%tsl         * rinterval
          out%tsn        = out%tsn         * rinterval
          out%tsns       = out%tsns        * rinterval
          ! Write value to file. Use ranges%Ebal for -999999 to 999999
          call write_ovar(out_timestep, ncid_out, ovid%hfds, 'hfds', &
               out%hfds, ranges%Ebal, patchout%hfds, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%hfdsn, 'hfdsn', &
               out%hfdsn, ranges%Ebal, patchout%hfdsn, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%hfls, 'hfls', &
               out%hfls, ranges%Ebal, patchout%hfls, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%hfmlt, 'hfmlt', &
               out%hfmlt, ranges%Ebal, patchout%hfmlt, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%hfrs, 'hfrs', &
               out%hfrs, ranges%Ebal, patchout%hfrs, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%hfsbl, 'hfsbl', &
               out%hfsbl, ranges%Ebal, patchout%hfsbl, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%hfss, 'hfss', &
               out%hfss, ranges%Ebal, patchout%hfss, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%rlus, 'rlus', &
               out%rlus, ranges%Ebal, patchout%rlus, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%rsus, 'rsus', &
               out%rsus, ranges%Ebal, patchout%rsus, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%esn, 'esn', &
               out%esn, ranges%Ebal, patchout%esn, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%evspsbl, 'evspsbl', &
               out%evspsbl, ranges%Ebal, patchout%evspsbl, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%evspsblsoi, 'evspsblsoi', &
               out%evspsblsoi, ranges%Ebal, patchout%evspsblsoi, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%evspsblveg, 'evspsblveg', &
               out%evspsblveg, ranges%Ebal, patchout%evspsblveg, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%mrrob, 'mrrob', &
               out%mrrob, ranges%Ebal, patchout%mrrob, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%mrros, 'mrros', &
               out%mrros, ranges%Ebal, patchout%mrros, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%sbl, 'sbl', &
               out%sbl, ranges%Ebal, patchout%sbl, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%snm, 'snm', &
               out%snm, ranges%Ebal, patchout%snm, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%snmsl, 'snmsl', &
               out%snmsl, ranges%Ebal, patchout%snmsl, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%tran, 'tran', &
               out%tran, ranges%Ebal, patchout%tran, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%albs, 'albs', &
               out%albs, ranges%Ebal, patchout%albs, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%albsn, 'albsn', &
               out%albsn, ranges%Ebal, patchout%albsn, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%cw, 'cw', &
               out%cw, ranges%Ebal, patchout%cw, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%lqsn, 'lqsn', &
               out%lqsn, ranges%Ebal, patchout%lqsn, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%lwsnl, 'lwsnl', &
               out%lwsnl, ranges%Ebal, patchout%lwsnl, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%mrfsofr, 'mrfsofr', &
               out%mrfsofr, ranges%Ebal, patchout%mrfsofr, 'soil', met)
          call write_ovar(out_timestep, ncid_out, ovid%mrlqso, 'mrlqso', &
               out%mrlqso, ranges%Ebal, patchout%mrlqso, 'soil', met)
          call write_ovar(out_timestep, ncid_out, ovid%mrlsl, 'mrlsl', &
               out%mrlsl, ranges%Ebal, patchout%mrlsl, 'soil', met)
          call write_ovar(out_timestep, ncid_out, ovid%snc, 'snc', &
               out%snc, ranges%Ebal, patchout%snc, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%snd, 'snd', &
               out%snd, ranges%Ebal, patchout%snd, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%snw, 'snw', &
               out%snw, ranges%Ebal, patchout%snw, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%snwc, 'snwc', &
               out%snwc, ranges%Ebal, patchout%snwc, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%tcs, 'tcs', &
               out%tcs, ranges%Ebal, patchout%tcs, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%tgs, 'tgs', &
               out%tgs, ranges%Ebal, patchout%tgs, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%ts, 'ts', &
               out%ts, ranges%Ebal, patchout%ts, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%tsl, 'tsl', &
               out%tsl, ranges%Ebal, patchout%tsl, 'soil', met)
          call write_ovar(out_timestep, ncid_out, ovid%tsn, 'tsn', &
               out%tsn, ranges%Ebal, patchout%tsn, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%tsns, 'tsns', &
               out%tsns, ranges%Ebal, patchout%tsns, 'default', met)
          ! Reset temporary output variable
          out%hfds       = zero4
          out%hfdsn      = zero4
          out%hfls       = zero4
          out%hfmlt      = zero4
          out%hfrs       = zero4
          out%hfsbl      = zero4
          out%hfss       = zero4
          out%rlus       = zero4
          out%rsus       = zero4
          out%esn        = zero4
          out%evspsbl    = zero4
          out%evspsblsoi = zero4
          out%evspsblveg = zero4
          out%mrrob      = zero4
          out%mrros      = zero4
          out%sbl        = zero4
          out%snm        = zero4
          out%snmsl      = zero4
          out%tran       = zero4
          out%albs       = zero4
          out%albsn      = zero4
          out%cw         = zero4
          out%lqsn       = zero4
          out%lwsnl      = zero4
          out%mrfsofr    = zero4
          out%mrlqso     = zero4
          out%mrlsl      = zero4
          out%snc        = zero4
          out%snd        = zero4
          out%snw        = zero4
          out%snwc       = zero4
          out%tcs        = zero4
          out%tgs        = zero4
          out%ts         = zero4
          out%tsl        = zero4
          out%tsn        = zero4
          out%tsns       = zero4
       END IF
    END IF

    ! NBP and turnover fluxes [umol/m^2/s]
    IF (output%carbon .OR. output%NBP) THEN
       ! Add current timestep's value to total of temporary output variable:
       IF (cable_user%POPLUC) THEN
          out%NBP = out%NBP + (-toreal4((casaflux%Crsoil-casaflux%cnpp &
               - casapool%dClabiledt) * gd2umols)) !- &
          !REAL((casaflux%FluxCtohwp + casaflux%FluxCtoclear  )/86400.0_r_2 &
          !/ 1.201E-5_r_2, 4)
       ELSE
          out%NBP = out%NBP + (-toreal4((casaflux%Crsoil-casaflux%cnpp &
               - casapool%dClabiledt) * gd2umols))
       ENDIF

       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%NBP = out%NBP * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%NBP, 'NBP', out%NBP, &
               ranges%NEE, patchout%NBP, 'default', met)
          ! Reset temporary output variable:
          out%NBP = zero4
       END IF
    ENDIF
    IF (output%casa) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%dCdt = out%dCdt + toreal4((casapool%ctot-casapool%ctot_0) * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%dCdt = out%dCdt * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%dCdt, 'dCdt', out%dCdt, &
               ranges%NEE, patchout%dCdt, 'default', met)
          ! Reset temporary output variable:
          out%dCdt = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnover = out%PlantTurnover + &
            toreal4(sum(casaflux%Cplant_turnover,2) * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnover = out%PlantTurnover * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnover, 'PlantTurnover', out%PlantTurnover, &
               ranges%NEE, patchout%PlantTurnover, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnover = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnoverLeaf = out%PlantTurnoverLeaf + &
            toreal4(casaflux%Cplant_turnover(:,1) * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnoverLeaf = out%PlantTurnoverLeaf * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnoverLeaf, 'PlantTurnoverLeaf', out%PlantTurnoverLeaf, &
               ranges%NEE, patchout%PlantTurnoverLeaf, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnoverLeaf = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnoverFineRoot = out%PlantTurnoverFineRoot + &
            toreal4(casaflux%Cplant_turnover(:,3) * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnoverFineRoot = out%PlantTurnoverFineRoot * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnoverFineRoot, 'PlantTurnoverFineRoot', &
               out%PlantTurnoverFineRoot, &
               ranges%NEE, patchout%PlantTurnoverFineRoot, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnoverFineRoot = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnoverWood = out%PlantTurnoverWood + &
            toreal4((casaflux%Cplant_turnover(:,2)) * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnoverWood = out%PlantTurnoverWood * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnoverWood, 'PlantTurnoverWood', out%PlantTurnoverWood, &
               ranges%NEE, patchout%PlantTurnoverWood, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnoverWood = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnoverWoodDist = out%PlantTurnoverWoodDist + &
            toreal4(casaflux%Cplant_turnover_disturbance * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnoverWoodDist = out%PlantTurnoverWoodDist * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnoverWoodDist, 'PlantTurnoverWoodDist', &
               out%PlantTurnoverWoodDist, &
               ranges%NEE, patchout%PlantTurnoverWoodDist, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnoverWoodDist = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnoverWoodCrowding = out%PlantTurnoverWoodCrowding + &
            toreal4(casaflux%Cplant_turnover_crowding * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnoverWoodCrowding = out%PlantTurnoverWoodCrowding * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnoverWoodCrowding, 'PlantTurnoverWoodCrowding', &
               out%PlantTurnoverWoodCrowding, &
               ranges%NEE, patchout%PlantTurnoverWoodCrowding, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnoverWoodCrowding = zero4
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%PlantTurnoverWoodResourceLim = out%PlantTurnoverWoodResourceLim + &
            toreal4(casaflux%Cplant_turnover_resource_limitation * gd2umols)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantTurnoverWoodResourceLim = out%PlantTurnoverWoodResourceLim * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantTurnoverWoodResourceLim, &
               'PlantTurnoverWoodResourceLim', &
               out%PlantTurnoverWoodResourceLim, &
               ranges%NEE, patchout%PlantTurnoverWoodResourceLim, 'default', met)
          ! Reset temporary output variable:
          out%PlantTurnoverWoodResourceLim = zero4
       END IF

       IF (cable_user%POPLUC) THEN
          ! Add current timestep's value to total of temporary output variable:
          out%LandUseFlux = out%LandUseFlux + &
               toreal4((casaflux%FluxCtohwp + casaflux%FluxCtoclear) * gd2umols)
          IF (writenow) THEN
             ! Divide accumulated variable by number of accumulated time steps:
             out%LandUseFlux = out%LandUseFlux * rinterval
             ! Write value to file:
             CALL write_ovar(out_timestep, ncid_out, ovid%LandUseFlux, 'LandUseFlux', &
                  out%LandUseFlux, &
                  ranges%NEE, patchout%LandUseFlux, 'default', met)
             ! Reset temporary output variable:
             out%LandUseFlux = zero4
          END IF
       ENDIF

    END IF

    ! plant carbon [kg C m-2]
    IF (output%casa) THEN
       out%TotSoilCarb = out%TotSoilCarb + toreal4((SUM(casapool%csoil,2)+SUM(casapool%clitter,2)) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%TotSoilCarb = out%TotSoilCarb * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%TotSoilCarb, 'TotSoilCarb', out%TotSoilCarb, &
               ranges%TotSoilCarb, patchout%TotSoilCarb, 'default', met)
          ! Reset temporary output variable:
          out%TotSoilCarb = zero4
       END IF

       out%TotLittCarb = out%TotLittCarb + toreal4(SUM(casapool%clitter,2) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%TotLittCarb = out%TotLittCarb * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%TotLittCarb, 'TotLittCarb', out%TotLittCarb, &
               ranges%TotLittCarb, patchout%TotLittCarb, 'default', met)
          ! Reset temporary output variable:
          out%TotLittCarb = zero4
       END IF

       out%SoilCarbFast = out%SoilCarbFast + toreal4(casapool%csoil(:,1) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilCarbFast = out%SoilCarbFast * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilCarbFast, 'SoilCarbFast' &
               , out%SoilCarbFast, &
               ranges%TotSoilCarb, patchout%SoilCarbFast, 'default', met)
          ! Reset temporary output variable:
          out%SoilCarbFast = zero4
       END IF

       out%SoilCarbSlow = out%SoilCarbSlow + toreal4(casapool%csoil(:,2)/ 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilCarbSlow = out%SoilCarbSlow * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilCarbSlow, 'SoilCarbSlow' &
               , out%SoilCarbSlow, &
               ranges%TotSoilCarb, patchout%SoilCarbSlow, 'default', met)
          ! Reset temporary output variable:
          out%SoilCarbSlow = zero4
       END IF

       out%SoilCarbPassive = out%SoilCarbPassive + toreal4(casapool%csoil(:,3) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilCarbPassive = out%SoilCarbPassive * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilCarbPassive, 'SoilCarbPassive' &
               , out%SoilCarbPassive, &
               ranges%TotSoilCarb, patchout%SoilCarbPassive, 'default', met)
          ! Reset temporary output variable:
          out%SoilCarbPassive = zero4
       END IF

       out%LittCarbMetabolic = out%LittCarbMetabolic + toreal4(casapool%clitter(:,1) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LittCarbMetabolic = out%LittCarbMetabolic * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LittCarbMetabolic, 'LittCarbMetabolic', out%LittCarbMetabolic, &
               ranges%TotLittCarb, patchout%LittCarbMetabolic, 'default', met)
          ! Reset temporary output variable:
          out%LittCarbMetabolic = zero4
       END IF

       out%LittCarbStructural = out%LittCarbStructural + toreal4(casapool%clitter(:,2) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LittCarbStructural = out%LittCarbStructural * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LittCarbStructural, 'LittCarbStructural', out%LittCarbStructural, &
               ranges%TotLittCarb, patchout%LittCarbStructural, 'default', met)
          ! Reset temporary output variable:
          out%LittCarbStructural = zero4
       END IF

       out%LittCarbCWD = out%LittCarbCWD + toreal4(casapool%clitter(:,3) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LittCarbCWD = out%LittCarbCWD * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LittCarbCWD, 'LittCarbCWD', out%LittCarbCWD, &
               ranges%TotLittCarb, patchout%LittCarbCWD, 'default', met)
          ! Reset temporary output variable:
          out%LittCarbCWD = zero4
       END IF

       out%PlantCarbLeaf = out%PlantCarbLeaf + toreal4(casapool%cplant(:,1) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantCarbLeaf = out%PlantCarbLeaf * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantCarbLeaf, 'PlantCarbLeaf', out%PlantCarbLeaf, &
               ranges%TotLittCarb, patchout%PlantCarbLeaf, 'default', met)
          ! Reset temporary output variable:
          out%PlantCarbLeaf = zero4
       END IF

       out%PlantCarbFineRoot = out%PlantCarbFineRoot + toreal4(casapool%cplant(:,3) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantCarbFineRoot = out%PlantCarbFineRoot * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantCarbFineRoot, 'PlantCarbFineRoot', &
               out%PlantCarbFineRoot, &
               ranges%TotLittCarb, patchout%PlantCarbFineRoot, 'default', met)
          ! Reset temporary output variable:
          out%PlantCarbFineRoot = zero4
       END IF

       out%PlantCarbWood = out%PlantCarbWood + toreal4(casapool%cplant(:,2) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PlantCarbWood = out%PlantCarbWood * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PlantCarbWood, 'PlantCarbWood', out%PlantCarbWood, &
               ranges%TotLittCarb, patchout%PlantCarbWood, 'default', met)
          ! Reset temporary output variable:
          out%PlantCarbWood = zero4
       END IF

       out%TotLivBiomass = out%TotLivBiomass + toreal4((SUM(casapool%cplant,2)) / 1000.0_r_2)
       IF (writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%TotLivBiomass = out%TotLivBiomass * rinterval
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%TotLivBiomass, 'TotLivBiomass', out%TotLivBiomass, &
               ranges%TotLivBiomass, patchout%TotLivBiomass, 'default', met)
          ! Reset temporary output variable:
          out%TotLivBiomass = zero4
       END IF

    END IF

    ! 13C
    if (cable_user%c13o2 .and. output%c13o2) then
       ! Add current timestep's value to total of temporary output variable
       ! 12C
       out%An      = out%An      + sum(canopy%An,2)
       out%Rd      = out%Rd      + sum(canopy%Rd,2)
       out%cplant  = out%cplant  + casapool%cplant
       out%clitter = out%clitter + casapool%clitter
       out%csoil   = out%csoil   + casapool%csoil
       out%clabile = out%clabile + casapool%clabile
       ! 13C
       out%A13n      = out%A13n      + sum(c13o2flux%An,2)
       !MC - Check with Vanessa that no leaf area
       out%aDisc13   = out%aDisc13   + sum((canopy%An+canopy%Rd)*c13o2flux%Disc,2)
       out%c13plant  = out%c13plant  + c13o2pools%cplant
       out%c13litter = out%c13litter + c13o2pools%clitter
       out%c13soil   = out%c13soil   + c13o2pools%csoil
       out%c13labile = out%c13labile + c13o2pools%clabile
       if (writenow) then
          ! Divide accumulated variable by number of accumulated time steps
          ! 12C
          out%An      = out%An      * r2interval
          out%Rd      = out%Rd      * r2interval
          out%cplant  = out%cplant  * r2interval
          out%clitter = out%clitter * r2interval
          out%csoil   = out%csoil   * r2interval
          out%clabile = out%clabile * r2interval
          ! 13C
          out%A13n      = out%A13n      * r2interval
          out%aDisc13   = out%aDisc13   * r2interval
          out%c13plant  = out%c13plant  * r2interval
          out%c13litter = out%c13litter * r2interval
          out%c13soil   = out%c13soil   * r2interval
          out%c13labile = out%c13labile * r2interval
          ! Write value to file
          ! 12C
          call write_ovar(out_timestep, ncid_out, ovid%An, 'An', &
               toreal4(out%An), ranges%Wbal, patchout%c13o2, 'default', met) ! Wbal ranges -999999 to 999999
          call write_ovar(out_timestep, ncid_out, ovid%Rd, 'Rd', &
               toreal4(out%Rd), ranges%Wbal, patchout%c13o2, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%cplant, 'Cplant', &
               toreal4(out%cplant), ranges%TotLivBiomass, patchout%c13o2, 'generic', met)
          call write_ovar(out_timestep, ncid_out, ovid%clitter, 'Clitter', &
               toreal4(out%clitter), ranges%TotLittCarb, patchout%c13o2, 'generic', met)
          call write_ovar(out_timestep, ncid_out, ovid%csoil, 'Csoil', &
               toreal4(out%csoil), ranges%TotSoilCarb, patchout%c13o2, 'generic', met)
          call write_ovar(out_timestep, ncid_out, ovid%clabile, 'Clabile', &
               toreal4(out%clabile), ranges%TotLivBiomass, patchout%c13o2, 'default', met)
          ! 13C
          call write_ovar(out_timestep, ncid_out, ovid%A13n, 'An13', &
               toreal4(out%A13n), ranges%Wbal, patchout%c13o2, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%aDisc13, 'aDisc13', &
               toreal4(out%aDisc13), ranges%Wbal, patchout%c13o2, 'default', met)
          call write_ovar(out_timestep, ncid_out, ovid%c13plant, 'C13plant', &
               toreal4(out%c13plant), ranges%TotLivBiomass, patchout%c13o2, 'generic', met)
          call write_ovar(out_timestep, ncid_out, ovid%c13litter, 'C13litter', &
               toreal4(out%c13litter), ranges%TotLittCarb, patchout%c13o2, 'generic', met)
          call write_ovar(out_timestep, ncid_out, ovid%c13soil, 'C13soil', &
               toreal4(out%c13soil), ranges%TotSoilCarb, patchout%c13o2, 'generic', met)
          call write_ovar(out_timestep, ncid_out, ovid%c13labile, 'C13labile', &
               toreal4(out%c13labile), ranges%TotLivBiomass, patchout%c13o2, 'default', met)
          ! Reset temporary output variable
          ! 12C
          out%An      = 0.0_r_2
          out%Rd      = 0.0_r_2
          out%cplant  = 0.0_r_2
          out%clitter = 0.0_r_2
          out%csoil   = 0.0_r_2
          out%clabile = 0.0_r_2
          ! 13C
          out%A13n      = 0.0_r_2
          out%aDisc13   = 0.0_r_2
          out%c13plant  = 0.0_r_2
          out%c13litter = 0.0_r_2
          out%c13soil   = 0.0_r_2
          out%c13labile = 0.0_r_2
       endif
    endif

    !MC - Do we need this?
    !     From the netcdf documentation:
    !         The function NF90 SYNC offers a way to synchronize the disk copy of a netCDF dataset
    !         with in-memory buffers. There are two reasons you might want to synchronize after writes:
    !         - To minimize data loss in case of abnormal termination, or
    !         - To make data available to other processes for reading immediately after it is written.
    !     Both is not the case.
    ok = nf90_sync(ncid_out)

  END SUBROUTINE write_output

  !=============================================================================

  SUBROUTINE close_output_file(bal)
    ! Closes output file, reports cumulative mass and energy
    ! balances, and deallocates variables.

    implicit none

    TYPE(balances_type), INTENT(IN) :: bal

    INTEGER :: i ! do loop counter

    ! Close file
    ok = NF90_CLOSE(ncid_out)
    ncid_out = -1
    IF (ok /= NF90_NOERR) &
         CALL nc_abort(ok, 'Error closing output file '//TRIM(filename%out)//'(SUBROUTINE close_output_file)')

    ! Report balance info to log file if verbose writing is requested:
    IF (output%balances .AND. verbose) THEN
       WRITE(logn,*) ''
       DO i = 1, mland
          WRITE(logn,'(A51,I7,1X,A11,E12.4,A6)') &
               ' Cumulative energy balance for each patch in site #', &
               i,'is (W/m^2):'
          WRITE(logn,*) &
               bal%ebal_tot(landpt(i)%cstart:landpt(i)%cstart + landpt(i)%nap - 1)
          WRITE(logn,'(A50,I7,1X,A8,E12.4,A3)') &
               ' Cumulative water balance for each patch in site #', &
               i,'is (mm):'
          WRITE(logn,*) &
               bal%wbal_tot(landpt(i)%cstart:landpt(i)%cstart + landpt(i)%nap - 1)
          WRITE(logn,*) ''
       END DO
    END IF

    ! Successful run!
    WRITE(logn,*) ''
    WRITE(logn,*) 'Run finished and output file closed.'

  END SUBROUTINE close_output_file

  !=============================================================================

  SUBROUTINE create_restart(logn, dels, ktau, soil, veg, ssnow, &
       canopy, rad, bgc, bal)

    ! Creates a restart file for CABLE using a land only grid with mland
    ! land points and max_vegpatches veg/soil patches (some of which may
    ! not be active). It uses CABLE's internal variable names.

    implicit none

    integer,                   intent(in) :: logn   ! log file number
    real,                      intent(in) :: dels   ! time step size
    integer,                   intent(in) :: ktau   ! timestep number in loop which include spinup
    type(soil_parameter_type), intent(in) :: soil   ! soil parameters
    type(veg_parameter_type),  intent(in) :: veg    ! vegetation parameters
    type(soil_snow_type),      intent(in) :: ssnow  ! soil and snow variables
    type(bgc_pool_type),       intent(in) :: bgc    ! carbon pool variables
    type(canopy_type),         intent(in) :: canopy ! vegetation variables
    type(radiation_type),      intent(in) :: rad    ! radiation variables
    type(balances_type),       intent(in) :: bal    ! energy and water balance variables
    ! integer, intent(in) :: mvtype
    ! integer, intent(in) :: mstype

    type(parID_type) :: rpid ! parameter IDs for restart nc file
    integer :: ncid_restart ! netcdf restart file id
    ! real, pointer,dimension(:,:) :: surffrac ! fraction of each surf type
    integer :: dummy ! dummy argument in subroutine call
    integer :: mlandID, mpID, radID, soilID, napID, &
         soilcarbID, plantcarbID, tID, snowID ! dimension IDs
    !    integer :: mlandID, surftypeID, patchID, radID, soilID, &
    !         soilcarbID, plantcarbID, tID, snowID ! dimension IDs
    integer :: tvarID, latID, lonID !,surffracID ! time,lat,lon variable ID
    integer :: tggID, wbID, wbiceID, tssID, ssdnnID, ssdnID, osnowdID, &
         smassID, sdepthID, snageID, snowdID, rtsoilID, isflagID, &
         canstoID, albsoilsnID, gammzzID, tggsnID, sghfluxID, &
         ghfluxID, runoffID, rnof1ID, rnof2ID, gaID, dgdtgID, &
         fevID, fesID, fhsID, wbtot0ID, osnowd0ID, cplantID, &
         csoilID, tradID, cfrdID, vlaiID, albedoID, otssID, &
         wetfacID, fwsoilID, hcID, usID
    integer :: h0ID, snowliqID, SID, TsurfaceID, scondsID, nsnowID, TsoilID
    character(len=10) :: todaydate, nowtime ! used to timestamp netcdf file
    ! character         :: FRST_OUT*100, CYEAR*4
    character :: frst_out*200, cyear*4

    dummy = 0 ! initialise

    WRITE(logn, '(A24)') ' Writing restart file...'
    IF ( TRIM(filename%path) .EQ. '' ) filename%path = './'
    frst_out = TRIM(filename%path)//'/'//TRIM(filename%restart_out)
    ! Look for explicit restart file (netCDF). If not, asssume input is path
    IF ( INDEX(TRIM(frst_out),'.nc',BACK=.TRUE.) .NE. LEN_TRIM(frst_out)-2 ) THEN
       WRITE( CYEAR,FMT="(I4)" ) CurYear + 1
       frst_out = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//&
            '_'//CYEAR//'_cable_rst.nc'
    ENDIF

    ! Create output file:
#ifdef __NETCDF3__
    ok = NF90_CREATE(trim(frst_out), ior(nf90_clobber,nf90_64bit_offset), ncid_restart)
#else
    ok = NF90_CREATE(trim(frst_out), ior(nf90_clobber,ior(nf90_netcdf4,nf90_classic_model)), ncid_restart)
#endif
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_restart, 'mland', mland, mlandID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining mland dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'mp', mp, mpID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining mp dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining vertical soil dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'snow', 3, snowID) ! number of snow layers
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining vertical snow dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'rad', nrb, radID) ! number of rad. bands
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining radiation dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'soil_carbon_pools', ncs, soilcarbID)
    ! number of soil carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining soil carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'plant_carbon_pools', ncp, plantcarbID)
    ! number of plant carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining plant carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'time', 1, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time dimension in restart file. '// &
         '(SUBROUTINE create_restart)')

    ! Define "time" variable and its attributes:
    ok=NF90_DEF_VAR(ncid_restart,'time',NF90_DOUBLE,(/tID/),tvarID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'calendar', calendar)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time variable attribute calendar in restart file. '// &
         '(SUBROUTINE create_restart)')

    ! Define latitude and longitude variable:
    ok=NF90_DEF_VAR(ncid_restart, 'latitude', NF90_FLOAT, (/mlandID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining latitude variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart,latID,'units','degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining latitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok=NF90_DEF_VAR(ncid_restart, 'longitude', NF90_FLOAT, (/mlandID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining longitude variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining longitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')

    ! Define number of active patches variable:
    ok = NF90_DEF_VAR(ncid_restart, 'nap', NF90_INT, (/mlandID/), napID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining nap variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, napID, 'long_name', &
         'Number of active patches')
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error defining nap variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define patch fraction variable:
    ok=nf90_def_var(ncid_restart, 'patchfrac', nf90_double, (/mpid/), &
         rpid%patchfrac)
    if (ok /= nf90_noerr) call nc_abort &
         (ok, 'Error defining patchfrac variable in restart file. '// &
         '(Subroutine create_restart)')
    ok = nf90_put_att(ncid_restart, rpid%patchfrac, 'long_name', &
         'Fraction of vegetated grid cell area occupied by a '// &
         'vegetation/soil patch')
    if (ok /= nf90_noerr) call nc_abort &
         (ok, 'Error defining patchfrac variable attributes in restart file. '// &
         '(Subroutine create_restart)')

    ! mvtype (Number of vegetation types):
    ok = NF90_DEF_VAR(ncid_restart, 'mvtype', NF90_INT, rpid%mvtype)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining mvtype variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%mvtype, "long_name", &
         "Number of vegetation types")
    ! mstype (Number of soil types):
    ok = NF90_DEF_VAR(ncid_restart, 'mstype', NF90_INT, rpid%mstype)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining mstype variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%mstype, "long_name", &
         "Number of soil types")

    !======begin defining state variables=======================================
    ! Interface arguments: netcdf file ID, variableID, variable name, variable
    ! units, variable long name, YES to write patch info (as this is a restart
    ! file), OPTIONAL extra dimension ID (e.g. for soil dimensioned variables),
    ! dimension switch to indicate what extra dimension is real or integer for
    ! single dim variables, xdimID,ydimID, zdimID (all three not used here),
    ! land dim ID, patch dim ID, YES we're writing a restart file.
    !------------------define soil states---------------------------------------
    CALL define_ovar(ncid_restart, tggID, 'tgg', 'K', &
         'Average layer soil temperature', &
         .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wbID, 'wb', 'vol/vol', &
         'Average layer volumetric soil moisture', &
         .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wbiceID, 'wbice', 'vol/vol', &
         'Average layer volumetric soil ice', &
         .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, tssID, 'tss', 'K', &
         'Combined soil/snow temperature', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, albsoilsnID, 'albsoilsn', '-', &
         'Combined soil/snow albedo', &
         .TRUE., radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rtsoilID, 'rtsoil', '??', &
         'Turbulent resistance for soil', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gammzzID, 'gammzz', 'J/kg/C', &
         'Heat capacity for each soil layer', &
         .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, runoffID, 'runoff', 'mm/timestep', &
         'Total runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rnof1ID, 'rnof1', 'mm/timestep', &
         'Surface runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rnof2ID, 'rnof2', 'mm/timestep', &
         'Subsurface runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! MC - 20220303 - add 2 ssnow parameter to restart
    CALL define_ovar(ncid_restart, otssID, 'otss', 'K', &
         'Combined soil/snow temperature at previous timestep', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wetfacID, 'wetfac', 'K', &
         'Surface wetness factor', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !---------------define snow states------------------------------------------
    CALL define_ovar(ncid_restart, tggsnID, 'tggsn', 'K', &
         'Average layer snow temperature', &
         .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ssdnnID, 'ssdnn', 'kg/m^3', &
         'Average snow density', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ssdnID, 'ssdn', 'kg/m^3', &
         'Average layer snow density', &
         .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, snowdID, 'snowd', 'mm', &
         'Liquid water eqivalent snow depth', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, snageID, 'snage', '??', &
         'Snow age', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, smassID, 'smass', 'kg/m^2', &
         'Average layer snow mass', &
         .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, sdepthID, 'sdepth', 'm', &
         'Snow layer depth', .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, osnowdID, 'osnowd', 'mm', &
         'Previous time step snow depth in water equivalent', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, isflagID, 'isflag', '-', &
         'Snow layer scheme flag', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    !----------------define canopy states----------------------------------
    CALL define_ovar(ncid_restart, canstoID, 'cansto', 'mm', &
         'Canopy surface water storage', .TRUE., 'real', 0, 0, 0, &
         mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ghfluxID, 'ghflux', 'W/m^2?', &
         '????', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, sghfluxID, 'sghflux', 'W/m^2?', &
         '????', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gaID, 'ga', 'W/m^2', &
         'Ground heat flux', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, dgdtgID, 'dgdtg', 'W/m^2/K', &
         'Derivative of ground heat flux wrt soil temperature', .TRUE., &
         'r2', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fevID, 'fev', 'W/m^2', &
         'Latent heat flux from vegetation', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fesID, 'fes', 'W/m^2', &
         'Latent heat flux from soil', &
         .TRUE., 'r2', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fhsID, 'fhs', 'W/m^2', &
         'Sensible heat flux from soil', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! MC - 20220303 - add 2 canopy parameters to restart
    CALL define_ovar(ncid_restart, fwsoilID, 'fwsoil', '-', &
         'Water limitation factor', .TRUE., &
         'r2', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, usID, 'us', 'm/s', &
         'Friction velocity', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !--------------biogeochemical variables------------------------
    CALL define_ovar(ncid_restart, cplantID, 'cplant', 'gC/m^2', &
         'Plant carbon stores', &
         .TRUE., plantcarbID, 'plantcarbon', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, csoilID, 'csoil', 'gC/m^2', &
         'Soil carbon stores', &
         .TRUE., soilcarbID, 'soilcarbon', 0, 0, 0, mpID, dummy, .TRUE.)
    !-------------------others---------------------------------
    CALL define_ovar(ncid_restart, wbtot0ID, 'wbtot0', 'mm', &
         'Initial time step soil water total', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, osnowd0ID, 'osnowd0', 'mm', &
         'Initial time step snow water total', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, albedoID, 'albedo', '-', &
         'Albedo for shortwave and NIR radiation', &
         .TRUE., radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, tradID, 'trad', 'K', &
         'Surface radiative temperature (soil/snow/veg inclusive)', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! MC - 20220303 - add 3 veg parameters to restart
    CALL define_ovar(ncid_restart, cfrdID, 'cfrd', '-', &
         'Ratio of day respiration to vcmax', &
         .true., 'real', 0, 0, 0, mpID, dummy, .true.)
    CALL define_ovar(ncid_restart, vlaiID, 'vlai', '-', &
         'Leaf area index', &
         .true., 'real', 0, 0, 0, mpID, dummy, .true.)
    CALL define_ovar(ncid_restart, hcID, 'hc', 'm', &
         'Roughness height of canopy', &
         .true., 'real', 0, 0, 0, mpID, dummy, .true.)
    !---------------------MODEL PARAMETERS---------------------------------
    WRITE(logn,'(A43)') '   Writing model parameters to restart file'
    CALL define_ovar(ncid_restart, rpid%iveg, 'iveg', '-', &
         'Vegetation type', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%isoil, 'isoil', '-', &
         'Soil type', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%clay, 'clay', '-', &
    !                  'Fraction of soil which is clay', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%sand, 'sand', '-', &
    !                  'Fraction of soil which is sand', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%silt, 'silt', '-', &
    !                  'Fraction of soil which is silt', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%ssat, 'ssat', '-', &
    !                  'Fraction of soil volume which is water @ saturation', &
    !                 .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%sfc, 'sfc', '-', &
    !                 'Fraction of soil volume which is water @ field capacity', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%swilt, 'swilt', '-', &
    !                  'Fraction of soil volume which is water @ wilting point', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! zse (depth of each soil layer):
    ok = NF90_DEF_VAR(ncid_restart, 'zse', NF90_FLOAT, (/soilID/), rpid%zse)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining zse variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%zse, "long_name", &
         "Depth of each soil layer")
    ok = NF90_PUT_ATT(ncid_restart, rpid%zse, "units", "m")
    ! CALL define_ovar(ncid_restart, rpid%froot, 'froot', '-', &
    !                  'Fraction of roots in each soil layer', &
    !                   .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%bch, 'bch', '-', &
    !                  'Parameter b, Campbell eqn 1985', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%hyds, 'hyds', 'm/s', &
    !                  'Hydraulic conductivity @ saturation', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%sucs, 'sucs', 'm', &
    !                  'Suction @ saturation', .TRUE., &
    !                  'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%css, 'css', 'J/kg/C', &
    !                  'Heat capacity of soil minerals', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%rhosoil, 'rhosoil', 'kg/m^3', &
    !                  'Density of soil minerals', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%rs20, 'rs20', '-', &
    !                  'Soil respiration coefficient at 20C', &
    !                   .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%albsoil, 'albsoil', '-', &
    !                  'Soil reflectance', .TRUE., &
    !                  radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%hc, 'hc', 'm', &
    !                  'Height of canopy', .TRUE., &
    !                  'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%canst1, 'canst1', 'mm/LAI', &
    !                  'Max water intercepted by canopy', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%dleaf, 'dleaf', 'm', &
    !                  'Chararacteristic length of leaf', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%frac4, 'frac4', '-', &
    !                  'Fraction of plants which are C4', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%ejmax, 'ejmax', 'mol/m^2/s', &
    !                  'Max potential electron transport rate top leaf', .TRUE., &
    !                  'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%vcmax, 'vcmax', 'mol/m^2/s', &
    !                  'Maximum RuBP carboxylation rate top leaf', .TRUE., &
    !                  'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%rp20, 'rp20', '-', &
    !                  'Plant respiration coefficient at 20C', .TRUE., 'real', &
    !                  0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%g0, 'g0', '-', &
    !                  'g0 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
    !                  0, 0, 0, mpID, dummy, .TRUE.) ! Ticket #56
    ! CALL define_ovar(ncid_restart, rpid%g1, 'g1', '-', &
    !                  'g1 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
    !                  0, 0, 0, mpID, dummy, .TRUE.)  ! Ticket #56
    ! CALL define_ovar(ncid_restart, rpid%rpcoef, 'rpcoef', '1/C', &
    !                  'Temperature coef nonleaf plant respiration', .TRUE., &
    !                  'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%shelrb, 'shelrb', '-', &
    !           'Sheltering factor', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%xfang, 'xfang', '-', &
    !        'Leaf angle parameter', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%wai, 'wai', '-', &
    !             'Wood area index', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%vegcf, 'vegcf', '-', &
    !                  'vegcf', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%extkn, 'extkn', '-', &
    !                  'Extinction coef for vertical nitrogen profile', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%tminvj, 'tminvj', 'C', &
    !                  'Min temperature for the start of photosynthesis', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%tmaxvj, 'tmaxvj', 'C', &
    !                  'Max temperature for the start of photosynthesis', &
    !                   .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%vbeta, 'vbeta', '-', &
    !                  'Stomatal sensitivity to soil water', &
    !                   .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%xalbnir, 'xalbnir', '-', &
    !                  'modifier for albedo in near ir band', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! ! ratecp (Plant carbon rate constant):
    ! ok = NF90_DEF_VAR(ncid_restart, 'ratecp', NF90_FLOAT, (/plantcarbID/), &
    !                   rpid%ratecp)
    ! IF (ok /= NF90_NOERR) CALL nc_abort &
    !                  (ok, 'Error defining ratecp variable in restart file. '// &
    !                   '(SUBROUTINE create_restart)')
    ! ok = NF90_PUT_ATT(ncid_restart, rpid%ratecp, "long_name", &
    !                   "Plant carbon rate constant")
    ! ok = NF90_PUT_ATT(ncid_restart, rpid%ratecp, "units", "1/year")
    ! ! ratecs (Soil carbon rate constant):
    ! ok = NF90_DEF_VAR(ncid_restart, 'ratecs', NF90_FLOAT, (/soilcarbID/), &
    !                   rpid%ratecs)
    ! IF (ok /= NF90_NOERR) CALL nc_abort &
    !                  (ok, 'Error defining ratecs variable in restart file. '// &
    !                   '(SUBROUTINE create_restart)')
    ! ok = NF90_PUT_ATT(ncid_restart, rpid%ratecs, "long_name", &
    !                   "Soil carbon rate constant")
    ! ok = NF90_PUT_ATT(ncid_restart, rpid%ratecs, "units", "1/year")
    ! CALL define_ovar(ncid_restart, rpid%meth, 'meth', '-', &
    !                  'Canopy turbulence parameterisation switch', &
    !                  .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%za_uv, 'za_uv', 'm', &
    !                 'Reference height (lowest atm. model layer) for momentum', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! CALL define_ovar(ncid_restart, rpid%za_tq, 'za_tq', 'm', &
    !                  'Reference height (lowest atm. model layer) for scalars', &
    !                  .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !
    ! IF(cable_user%SOIL_STRUC=='sli'.OR.cable_user%FWSOIL_SWITCH=='Haverd2013') THEN
    !   CALL define_ovar(ncid_restart,rpid%gamma,'gamma','-', &
    !         'Parameter in root efficiency function (Lai and Katul 2000)', &
    !         .TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
    ! ENDIF
    ! Soil-Litter-Iso soil model
    IF (cable_user%SOIL_STRUC=='sli') THEN
       ! Parameters for SLI:
       ! CALL define_ovar(ncid_restart,rpid%nhorizons,'nhorizons','-', &
       !      'Number of soil horizons',.TRUE.,'integer',0,0,0,mpID,dummy,.TRUE.)
       ! CALL define_ovar(ncid_restart,rpid%zeta,'zeta','[ ]', &
       !      'exponent factor in Topmodel eq',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       ! CALL define_ovar(ncid_restart,rpid%fsatmax,'fsatmax','[ ]', &
       !      'param in Topmodel eq',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       ! CALL define_ovar(ncid_restart,rpid%ishorizon,'ishorizon','-', &
       !      'Horizon number',.TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
       ! CALL define_ovar(ncid_restart,rpid%clitt,'clitt','tC/ha', &
       !      'Litter layer carbon content',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       ! CALL define_ovar(ncid_restart,rpid%ZR,'ZR','cm', &
       !      'Maximum rooting depth',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       ! CALL define_ovar(ncid_restart,rpid%F10,'F10','-', &
       !      'Fraction of roots in top 10 cm', &
       !      .TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       ! Variables for SLI:
       CALL define_ovar(ncid_restart, SID, 'S', '-', &
            'Fractional soil moisture content relative to saturated value', &
            .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
       CALL define_ovar(ncid_restart, TsoilID, 'Tsoil', 'degC', &
            'Tsoil', &
            .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
       !MC - Should be r_2 but specific subroutine not coded yet in cable_write
       CALL define_ovar(ncid_restart, snowliqID, 'snowliq', 'mm', &
            'liquid water content of snowpack', &
            .TRUE., snowID, 'r2snow', 0, 0, 0, mpID, dummy, .TRUE.)
       CALL define_ovar(ncid_restart,scondsID,'sconds','Wm-1K-1', &
            'thermal cond of snowpack', &
            .TRUE.,snowID,'snow',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart, h0ID, 'h0', 'm', &
            'Pond height above soil', &
            .TRUE., 'r2', 0, 0, 0, mpID, dummy, .TRUE.)
       CALL define_ovar(ncid_restart,nsnowID,'nsnow','-', &
            'number of snow layers', &
            .TRUE.,'integer',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart, TsurfaceID, 'Tsurface', 'degC', &
            'soil or snow surface T', &
            .TRUE., 'r2', 0, 0, 0, mpID, dummy, .TRUE.)
    END IF ! SLI soil model

    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Production", &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to ' &
         //TRIM(frst_out)// ' (SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Source", &
         'CABLE LSM restart file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to ' &
         //TRIM(frst_out)// ' (SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "CABLE_input_file", &
         TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to ' &
         //TRIM(frst_out)// ' (SUBROUTINE create_restart)')

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write time variable:
    ok = NF90_PUT_VAR(ncid_restart, tvarID, real(ktau,r_2)*real(dels,r_2))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error time variable to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write latitude and longitude variables:
    ok = NF90_PUT_VAR(ncid_restart, latID, latitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
         'Error writing latitude variable to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, lonID, longitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
         'Error writing longitude variable to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write number of active patches for each land grid cell:
    ok = NF90_PUT_VAR(ncid_restart, napID, landpt(:)%nap)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
         'Error writing nap variable to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write vegetated patch fractions
    ok = NF90_PUT_VAR(ncid_restart, rpid%patchfrac, &
         patch(:)%frac, start = (/1/), count = (/mp/))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing patchfrac to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write number of veg and soil types
    ok = NF90_PUT_VAR(ncid_restart, rpid%mvtype, mvtype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
         'Error writing mvtype parameter to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rpid%mstype, mstype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
         'Error writing mstype parameter to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write parameters:
    CALL write_ovar(ncid_restart, rpid%iveg, 'iveg', toreal4(veg%iveg), &
         ranges%iveg, .TRUE., 'integer', .TRUE.)
    CALL write_ovar(ncid_restart, rpid%isoil, 'isoil', toreal4(soil%isoilm), &
         ranges%isoil, .TRUE., 'integer', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%bch, 'bch', toreal4(soil%bch), &
    !                  ranges%bch, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%clay, 'clay', toreal4(soil%clay), &
    !                  ranges%clay, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%sand, 'sand', toreal4(soil%sand), &
    !                  ranges%sand, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%silt, 'silt', toreal4(soil%silt), &
    !                  ranges%silt, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%css, 'css', toreal4(soil%css), &
    !                  ranges%css, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%rhosoil, 'rhosoil', &
    !                  toreal4(soil%rhosoil), ranges%rhosoil, .TRUE., 'real', &
    !                  .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%hyds, 'hyds', toreal4(soil%hyds), &
    !                  ranges%hyds, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%sucs, 'sucs', toreal4(soil%sucs), &
    !                  ranges%sucs, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%rs20, 'rs20', toreal4(veg%rs20), &
    !                  ranges%rs20, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%ssat, 'ssat', toreal4(soil%ssat), &
    !                  ranges%ssat, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%sfc, 'sfc', toreal4(soil%sfc), &
    !                  ranges%sfc, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%swilt, 'swilt', toreal4(soil%swilt), &
    !                  ranges%swilt, .TRUE., 'real', .TRUE.)
    ! ! Soil dimensioned variables/parameters:
    ! CALL write_ovar (ncid_restart, rpid%froot, 'froot', toreal4(veg%froot), &
    !                  ranges%froot, .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, tggID, 'tgg', toreal4(ssnow%tgg), &
         ranges%SoilTemp, .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, wbID, 'wb', ssnow%wb, ranges%SoilMoist, &
         .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, wbiceID, 'wbice', ssnow%wbice, &
         ranges%SoilMoist, .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, gammzzID, 'gammzz', ssnow%gammzz, &
         (/-99999.0, 9999999.0/), .TRUE., 'soil', .TRUE.)
    ! Snow dimensioned variables/parameters:
    CALL write_ovar (ncid_restart, ssdnID, 'ssdn', toreal4(ssnow%ssdn), &
         (/0.0, 9999.0/), .TRUE., 'snow', .TRUE.)
    CALL write_ovar (ncid_restart, smassID, 'smass', toreal4(ssnow%smass), &
         (/0.0, 9999.0/), .TRUE., 'snow', .TRUE.)
    CALL write_ovar (ncid_restart, sdepthID, 'sdepth', toreal4(ssnow%sdepth), &
         (/0.0, 9999.0/), .TRUE., 'snow', .TRUE.)
    CALL write_ovar (ncid_restart, tggsnID, 'tggsn', toreal4(ssnow%tggsn), &
         (/100.0, 300.0/), .TRUE., 'snow', .TRUE.)
    ! Other dims
    CALL write_ovar (ncid_restart, albsoilsnID, 'albsoilsn', &
         toreal4(ssnow%albsoilsn), (/0.0, 1.0/), .TRUE., 'radiation', .TRUE.)
    CALL write_ovar (ncid_restart, cplantID, 'cplant', toreal4(bgc%cplant), &
         (/-99999.0, 9999999.0/), .TRUE., 'plantcarbon', .TRUE.)
    CALL write_ovar (ncid_restart, csoilID, 'csoil', toreal4(bgc%csoil), &
         (/-99999.0, 9999999.0/), .TRUE., 'soilcarbon', .TRUE.)
    ok = NF90_PUT_VAR(ncid_restart, rpid%zse, toreal4(soil%zse))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing zse parameter to ' &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ! ! Single dim:
    ! CALL write_ovar (ncid_restart, rpid%albsoil, 'albsoil', &
    !                  toreal4(soil%albsoil), ranges%albsoil, .TRUE., &
    !                  'radiation', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%canst1, 'canst1', toreal4(veg%canst1), &
    !                  ranges%canst1, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%dleaf, 'dleaf', toreal4(veg%dleaf), &
    !                  ranges%dleaf, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%ejmax, 'ejmax', toreal4(veg%ejmax), &
    !                  ranges%ejmax, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%vcmax, 'vcmax', toreal4(veg%vcmax), &
    !                  ranges%vcmax, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%frac4, 'frac4', toreal4(veg%frac4), &
    !                  ranges%frac4, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%hc, 'hc', toreal4(veg%hc), &
    !                  ranges%hc, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%rp20, 'rp20', toreal4(veg%rp20), &
    !                  ranges%rp20, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%g0, 'g0', toreal4(veg%g0), &
    !                  ranges%g0, .TRUE., 'real', .TRUE.) ! Ticket #56
    ! CALL write_ovar (ncid_restart, rpid%g1, 'g1', toreal4(veg%g1), &
    !                  ranges%g1, .TRUE., 'real', .TRUE.) ! Ticket #56
    ! CALL write_ovar (ncid_restart, rpid%rpcoef, 'rpcoef', toreal4(veg%rpcoef), &
    !                  ranges%rpcoef, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%shelrb, 'shelrb', toreal4(veg%shelrb), &
    !                  ranges%shelrb, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%xfang, 'xfang', toreal4(veg%xfang), &
    !                  ranges%xfang, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%wai, 'wai', toreal4(veg%wai), &
    !                  ranges%wai, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%vegcf, 'vegcf', toreal4(veg%vegcf), &
    !                  ranges%vegcf, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%extkn, 'extkn', toreal4(veg%extkn), &
    !                  ranges%extkn, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%tminvj, 'tminvj', toreal4(veg%tminvj), &
    !                  ranges%tminvj, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%tmaxvj, 'tmaxvj', toreal4(veg%tmaxvj), &
    !                  ranges%tmaxvj, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%vbeta, 'vbeta', toreal4(veg%vbeta), &
    !                  ranges%vbeta, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%xalbnir, 'xalbnir', &
    !                  toreal4(veg%xalbnir), ranges%xalbnir, .TRUE., &
    !                  'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%tmaxvj, 'tmaxvj', toreal4(veg%tmaxvj), &
    !                  ranges%tmaxvj, .TRUE., 'real', .TRUE.)
    ! ok = NF90_PUT_VAR(ncid_restart, rpid%ratecp, toreal4(bgc%ratecp))
    ! IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
    !                                    'Error writing ratecp parameter to ' &
    !      //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ! ok = NF90_PUT_VAR(ncid_restart, rpid%ratecs, toreal4(bgc%ratecs))
    ! IF(ok /= NF90_NOERR) CALL nc_abort(ok, &
    !                                    'Error writing ratecs parameter to ' &
    !      //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ! CALL write_ovar (ncid_restart, rpid%meth, 'meth', toreal4(veg%meth), &
    !                  ranges%meth, .TRUE., 'integer', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%za_uv, 'za_uv', toreal4(rough%za_uv), &
    !                  ranges%za, .TRUE., 'real', .TRUE.)
    ! CALL write_ovar (ncid_restart, rpid%za_tq, 'za_tq', toreal4(rough%za_tq), &
    !                  ranges%za, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, tssID, 'tss', toreal4(ssnow%tss), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, ssdnnID, 'ssdnn', toreal4(ssnow%ssdnn), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, osnowdID, 'osnowd', toreal4(ssnow%osnowd), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, snageID, 'snage', toreal4(ssnow%snage), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, snowdID, 'snowd', toreal4(ssnow%snowd), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rtsoilID, 'rtsoil', toreal4(ssnow%rtsoil), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, isflagID, 'isflag', toreal4(ssnow%isflag), &
         (/-99999.0, 9999999.0/), .TRUE., 'integer', .TRUE.)
    CALL write_ovar (ncid_restart, canstoID, 'cansto', toreal4(canopy%cansto), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, sghfluxID, 'sghflux', &
         toreal4(canopy%sghflux), (/-99999.0, 9999999.0/), &
         .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, ghfluxID, 'ghflux', toreal4(canopy%ghflux), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, runoffID, 'runoff', toreal4(ssnow%runoff), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rnof1ID, 'rnof1', toreal4(ssnow%rnof1), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rnof2ID, 'rnof2', toreal4(ssnow%rnof2), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    ! MC - 20220303 - add 2 ssnow parameter to restart
    CALL write_ovar(ncid_restart, otssID, 'otss', toreal4(ssnow%otss), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar(ncid_restart, wetfacID, 'wetfac', toreal4(ssnow%wetfac), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar(ncid_restart, gaID, 'ga', toreal4(canopy%ga), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar(ncid_restart, dgdtgID, 'dgdtg', canopy%dgdtg, &
         (/-99999.0, 9999999.0/), .TRUE., 'r2', .TRUE.)
    CALL write_ovar (ncid_restart, fevID, 'fev', toreal4(canopy%fev), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, fesID, 'fes', canopy%fes, &
         (/-99999.0, 9999999.0/), .TRUE., 'r2', .TRUE.)
    CALL write_ovar (ncid_restart, fhsID, 'fhs', toreal4(canopy%fhs), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    ! MC - 20220303 - add 2 canopy parameters to restart
    CALL write_ovar(ncid_restart, fwsoilID, 'fwsoil', canopy%fwsoil, &
         (/-99999.0, 9999999.0/), .TRUE., 'r2', .TRUE.)
    CALL write_ovar(ncid_restart, usID, 'us', toreal4(canopy%us), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, wbtot0ID, 'wbtot0', toreal4(bal%wbtot0), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, osnowd0ID, 'osnowd0', toreal4(bal%osnowd0), &
         (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, albedoID, 'albedo', toreal4(rad%albedo), &
         ranges%Albedo, .TRUE., 'radiation', .TRUE.)
    CALL write_ovar (ncid_restart, tradID, 'trad', &
         toreal4(rad%trad), ranges%RadT, .TRUE., 'real', .TRUE.)
    ! MC - 20220303 - add 3 veg parameters to restart
    CALL write_ovar(ncid_restart, cfrdID, 'cfrd', &
         toreal4(veg%cfrd), (/-99999.0, 9999999.0/), .true., 'real', .true.)
    CALL write_ovar(ncid_restart, vlaiID, 'vlai', &
         toreal4(veg%vlai), (/-99999.0, 9999999.0/), .true., 'real', .true.)
    CALL write_ovar(ncid_restart, hcID, 'hc', &
         toreal4(veg%hc), (/-99999.0, 9999999.0/), .true., 'real', .true.)
    ! IF(cable_user%SOIL_STRUC=='sli'.OR.cable_user%FWSOIL_SWITCH=='Haverd2013') THEN
    !    CALL write_ovar (ncid_restart,rpid%gamma,'gamma', &
    !         toreal4(veg%gamma),(/-99999.0,99999.0/),.TRUE.,'real',.TRUE.)
    ! ENDIF
    !
    IF (cable_user%SOIL_STRUC=='sli') THEN
       ! Write SLI parameters:
       ! CALL write_ovar (ncid_restart,rpid%nhorizons,'nhorizons', &
       !      toreal4(soil%nhorizons),(/-99999.0,99999.0/),.TRUE.,'integer',.TRUE.)
       ! CALL write_ovar (ncid_restart,rpid%ishorizon,'ishorizon', &
       !      toreal4(soil%ishorizon),(/-99999.0,99999.0/),.TRUE.,'soil',.TRUE.)
       ! CALL write_ovar (ncid_restart,rpid%clitt,'clitt', &
       !      toreal4(veg%clitt),(/-99999.0,99999.0/),.TRUE.,'real',.TRUE.)
       ! CALL write_ovar (ncid_restart,rpid%ZR,'ZR', &
       !      toreal4(veg%ZR),(/-99999.0,99999.0/),.TRUE.,'real',.TRUE.)
       ! CALL write_ovar (ncid_restart,rpid%F10,'F10', &
       !      toreal4(veg%F10),(/-99999.0,99999.0/),.TRUE.,'real',.TRUE.)
       ! Write SLI variables:
       CALL write_ovar(ncid_restart,SID, 'S', ssnow%S, &
            (/0.0,1.5/), .TRUE., 'soil', .TRUE.)
       CALL write_ovar(ncid_restart, TsoilID, 'Tsoil', ssnow%Tsoil, &
            (/-100.0,100.0/), .TRUE., 'soil', .TRUE.)
       ! CALL write_ovar(ncid_restart, snowliqID, 'snowliq', ssnow%snowliq, &
       !      (/-99999.0,99999.0/), .TRUE., 'snow', .TRUE.) ! vh ! needs fixing
       CALL write_ovar(ncid_restart, snowliqID, 'snowliq', ssnow%snowliq, &
            (/-99999.0,99999.0/), .TRUE., 'snow', .TRUE.)
       CALL write_ovar(ncid_restart, scondsID, 'sconds', ssnow%sconds, &
            (/-99999.0,99999.0/), .TRUE., 'snow', .TRUE.)
       CALL write_ovar(ncid_restart, h0ID, 'h0', ssnow%h0, &
            (/-99999.0,99999.0/), .TRUE., 'real', .TRUE.)
       CALL write_ovar(ncid_restart, nsnowID, 'nsnow', toreal4(ssnow%nsnow), &
            (/-99999.0,99999.0/), .TRUE., 'integer',.TRUE.)
       CALL write_ovar(ncid_restart, TsurfaceID, 'Tsurface', ssnow%Tsurface, &
            (/-99999.0,99999.0/), .TRUE., 'real', .TRUE.)
    END IF

    ! Close restart file
    ok = NF90_CLOSE(ncid_restart)
    ncid_restart = -1

    WRITE(logn, '(A36)') '   Restart file complete and closed.'

  END SUBROUTINE create_restart

  ! ------------------------------------------------------------------

  elemental pure function dp2sp(var)

    use cable_def_types_mod, only: r_2

    implicit none

    real(r_2),   intent(in) :: var
    real(kind=4)            :: dp2sp

    real(r_2),    parameter :: tini = real(tiny(1.0),r_2)
    real(kind=4), parameter :: zero = real(0.0,4)

    if (abs(var) > tini) then
       dp2sp = real(var,4)
    else
       dp2sp = zero
    endif

    return

  end function dp2sp


  elemental pure function i2sp(var)

    implicit none

    integer  ,   intent(in) :: var
    real(kind=4)            :: i2sp

    i2sp = real(var,4)

    return

  end function i2sp


  elemental pure function sp2sp(var)

    implicit none

    real,        intent(in) :: var
    real(kind=4)            :: sp2sp

    real,         parameter :: tini = tiny(1.0)
    real(kind=4), parameter :: zero = real(0.0,4)

    if (abs(var) > tini) then
       sp2sp = real(var,4)
    else
       sp2sp = zero
    endif

    return

  end function sp2sp

  ! ------------------------------------------------------------------

END MODULE cable_output_module
