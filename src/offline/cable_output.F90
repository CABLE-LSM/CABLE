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


  USE cable_abort_module, ONLY: abort, nc_abort
  USE cable_def_types_mod
  USE casavariable, ONLY: casa_pool, casa_flux, casa_met
  USE cable_IO_vars_module
  USE cable_checks_module, ONLY: mass_balance, energy_balance, ranges, check_range
  USE cable_write_module
  USE netcdf
  USE cable_common_module, ONLY: filename, calcsoilalbedo, CurYear,IS_LEAPYEAR, cable_user,&
       gw_params
  IMPLICIT NONE
  PRIVATE
  PUBLIC open_output_file, write_output, close_output_file, create_restart, check_and_write, output_par_settings_type
  INTEGER :: ncid_out, ncid_restart ! output/restart data netcdf file ID
  REAL :: missing_value = -999999.0 ! for netcdf output
  TYPE out_varID_type ! output variable IDs in netcdf file
     INTEGER ::      SWdown, LWdown, Wind, Wind_E, PSurf,                       &
          Tair, Qair, Tscrn, Qscrn, Rainf, Snowf, CO2air,            &
          Tmx, Tmn, Txx, Tnn,                                        &
          Qmom, Qle, Qh, Qg, NEE, SWnet,                             &
          LWnet, SoilMoist, SoilTemp, Albedo,                        &
          visAlbedo, nirAlbedo, SoilMoistIce,                        &
          Qs, Qsb, Evap, PotEvap, BaresoilT, SWE, SnowT,             &
          RadT, VegT, Ebal, Wbal, AutoResp, RootResp,                &
          StemResp, LeafResp, HeteroResp, GPP, NPP, LAI,             &
          ECanop, TVeg, ESoil, CanopInt, SnowDepth,                  &
          HVeg, HSoil, Rnet, tvar, CanT,Fwsoil, RnetSoil, SnowMelt, &
          NBP, TotSoilCarb, TotLivBiomass, &
          TotLittCarb, SoilCarbFast, SoilCarbSlow, SoilCarbPassive, &
          LittCarbMetabolic, LittCarbStructural, LittCarbCWD, &
          PlantCarbLeaf, PlantCarbFineRoot, PlantCarbWood, &
          PlantTurnover, PlantTurnoverLeaf, PlantTurnoverFineRoot, &
          PlantTurnoverWood, PlantTurnoverWoodDist, PlantTurnoverWoodCrowding, &
          PlantTurnoverWoodResourceLim, dCdt, Area, LandUseFlux, patchfrac, &
          vcmax,hc,WatTable,GWMoist,SatFrac,Qrecharge
  END TYPE out_varID_type
  TYPE(out_varID_type) :: ovid ! netcdf variable IDs for output variables
  TYPE(parID_type) :: opid ! netcdf variable IDs for output variables
  TYPE output_temporary_type
     REAL(KIND=4), POINTER, DIMENSION(:) :: SWdown ! 6 downward short-wave
     ! radiation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LWdown ! 7 downward long-wave
     ! radiation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Rainf  ! 8 rainfall [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Snowf  ! 9 snowfall [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: PSurf  ! 10 surface pressure [Pa]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tair   ! 11 surface air temperature
     ! [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qair   ! 12 specific humidity [kg/kg]
     !INH new output variables
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tscrn  ! -- screen-level air
     ! temperature [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qscrn  ! -- screen level specific
     ! humidity [kg/kg]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tmx    ! -- averaged daily maximum
     ! screen level temp [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Txx    ! -- max screen level temp
     ! in averaging period [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tmn    ! -- averaged daily minimum
     ! screen level temp [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tnn    ! -- min screen level temp
     ! in averaging period [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tdaymx ! -- daily maximum
     ! screen level temp [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Tdaymn ! -- daily maximum
     ! screen level temp [oC]
     REAL(KIND=4), POINTER, DIMENSION(:) :: CO2air ! 13 CO2 concentration [ppmv]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wind   ! 14 windspeed [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wind_N ! 15 surface wind speed, N
     ! component [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wind_E ! 16 surface wind speed, E
     ! component [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LAI
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qmom   ! -- momentum flux [kg/m/s2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qh     ! 17 sensible heat flux [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qle    ! 18 latent heat flux [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qg     ! 19 ground heat flux [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SWnet  ! 20 net shortwave [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LWnet  ! 21 net longwave [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Evap   ! 22 total evapotranspiration
     ! [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Ewater ! 23 evap. from surface water
     ! storage [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: ESoil  ! 24 bare soil evaporation
     ! [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: TVeg   ! 25 vegetation transpiration
     ! [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: ECanop ! 26 interception evaporation
     ! [kg/m2/s]
     ! 27 potential evapotranspiration [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: PotEvap
     REAL(KIND=4), POINTER, DIMENSION(:) :: ACond   ! 28 aerodynamic conductance
     ! [m/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilWet ! 29 total soil wetness [-]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Albedo  ! 30 albedo [-]
     REAL(KIND=4), POINTER, DIMENSION(:) :: visAlbedo  ! vars intro for Ticket #27
     REAL(KIND=4), POINTER, DIMENSION(:) :: nirAlbedo  ! vars intro for Ticket #27
     REAL(KIND=4), POINTER, DIMENSION(:) :: VegT    ! 31 vegetation temperature
     ! [K]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilTemp  ! 32 av.layer soil
     ! temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilMoist ! 33 av.layer soil
     ! moisture [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilMoistIce ! 33 av.layer soil
     ! frozen moisture [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qs  ! 34 surface runoff [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Qsb ! 35 subsurface runoff [kg/m2/s]
     ! 36 change in soilmoisture (sum layers) [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: DelSoilMoist
     ! 37 change in snow water equivalent [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: DelSWE
     ! 38 change in interception storage [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: DelIntercept
     REAL(KIND=4), POINTER, DIMENSION(:) :: SnowT     ! 39 snow surface temp [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: BaresoilT ! 40 surface bare soil
     ! temp [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: AvgSurfT  ! 41 Average surface
     ! temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: RadT      ! 42 Radiative surface
     ! temperature [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SWE       ! 43 snow water equivalent
     ! [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: RootMoist ! 44 root zone soil
     ! moisture [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: CanopInt  ! 45 total canopy water
     ! storage [kg/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: NEE       ! 46 net ecosystem exchange
     ! [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: NPP       ! 47 net primary production
     ! of C by veg [umol/m2/s]
     ! 48 gross primary production C by veg [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: GPP
     REAL(KIND=4), POINTER, DIMENSION(:) :: AutoResp   ! 49 autotrophic
     ! respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: LeafResp   ! 51 autotrophic
     ! respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: HeteroResp ! 50 heterotrophic
     ! respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SnowDepth  ! actual depth of snow in
     ! [m]
     ! Non-Alma variables
     REAL(KIND=4), POINTER, DIMENSION(:) :: Rnet  ! net absorbed radiation [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: HVeg  ! sensible heat from vegetation
     ! [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: HSoil ! sensible heat from soil
     ! [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: RnetSoil ! latent heat from soil
     ! [kg/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: SnowMelt ! snow melt
     ! [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Ebal  ! cumulative energy balance
     ! [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Wbal  ! cumulative water balance
     ! [W/m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: CanT  ! within-canopy temperature
     ! [K]
     REAL(KIND=4), POINTER, DIMENSION(:) :: Fwsoil  ! soil-moisture modfier to stomatal conductance
     ! [-]

     ![umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: NBP
     REAL(KIND=4), POINTER, DIMENSION(:) :: dCdt
     ! [kg C /m2]
     REAL(KIND=4), POINTER, DIMENSION(:) :: TotSoilCarb
     REAL(KIND=4), POINTER, DIMENSION(:) :: TotLivBiomass
     REAL(KIND=4), POINTER, DIMENSION(:) :: TotLittCarb
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilCarbFast
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilCarbSlow
     REAL(KIND=4), POINTER, DIMENSION(:) :: SoilCarbPassive
     REAL(KIND=4), POINTER, DIMENSION(:) :: LittCarbMetabolic
     REAL(KIND=4), POINTER, DIMENSION(:) :: LittCarbStructural
     REAL(KIND=4), POINTER, DIMENSION(:) :: LittCarbCWD
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantCarbLeaf
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantCarbFineRoot
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantCarbWood
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnover
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverLeaf
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverFineRoot
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWood
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWoodDist
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWoodCrowding
     REAL(KIND=4), POINTER, DIMENSION(:) :: PlantTurnoverWoodResourceLim
     REAL(KIND=4), POINTER, DIMENSION(:) :: Area
     REAL(KIND=4), POINTER, DIMENSION(:) :: LandUseFlux
     REAL(KIND=4), POINTER, DIMENSION(:) :: vcmax
     REAL(KIND=4), POINTER, DIMENSION(:) :: patchfrac
     REAL(KIND=4), POINTER, DIMENSION(:) :: hc
     REAL(KIND=4), POINTER, DIMENSION(:)   :: SatFrac         !Saturated Fraction of Grid Cell
     REAL(KIND=4), POINTER, DIMENSION(:)   :: Qrecharge         !recharge rate Grid Cell
     REAL(KIND=4), POINTER, DIMENSION(:)   :: GWMoist       ! water balance of aquifer [mm3/mm3]
     REAL(KIND=4), POINTER, DIMENSION(:)   :: WatTable      ! water table depth [m]

     REAL(KIND=4), POINTER, DIMENSION(:) :: RootResp   !  autotrophic root respiration [umol/m2/s]
     REAL(KIND=4), POINTER, DIMENSION(:) :: StemResp   !  autotrophic stem respiration [umol/m2/s]
  END TYPE output_temporary_type

  TYPE output_var_settings_type
     TYPE(met_type), POINTER :: met
     LOGICAL :: writenow
     ! Optional
     CHARACTER(LEN=15) :: dimswitch = 'default'
  END TYPE output_var_settings_type
  TYPE output_par_settings_type
     TYPE(met_type), POINTER :: met
     ! Optional
     LOGICAL :: restart = .FALSE.
     CHARACTER(LEN=15) :: dimswitch = 'default'
  END TYPE output_par_settings_type
  TYPE(output_temporary_type), SAVE :: out
  INTEGER :: ok   ! netcdf error status


  ! Some Additional internal variables
  INTEGER :: out_timestep ! counter for output time steps

  INTERFACE check_and_write
      MODULE PROCEDURE :: check_and_write_d1
      MODULE PROCEDURE :: check_and_write_d2
      MODULE PROCEDURE :: check_and_write_d1_p
      MODULE PROCEDURE :: check_and_write_d2_p
  END INTERFACE check_and_write
  
  INTERFACE generate_out_write_acc
      MODULE PROCEDURE :: generate_out_write_acc_d1
      MODULE PROCEDURE :: generate_out_write_acc_d2
  END INTERFACE generate_out_write_acc
CONTAINS

  SUBROUTINE open_output_file(dels, soil, veg, bgc, rough, met)
    ! Creates netcdf output file, defines all variables
    ! and writes parameters to it if requested by user.
    REAL, INTENT(IN) :: dels ! time step size
    TYPE (soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)  :: veg  ! vegetation parameters
    TYPE (bgc_pool_type), INTENT(IN)       :: bgc
    TYPE (roughness_type), INTENT(IN)      :: rough
    TYPE (met_type), TARGET, INTENT(IN)       :: met
    ! REAL, POINTER :: surffrac(:, :) ! fraction of each surf type

    INTEGER :: xID, yID, zID, radID, soilID, soilcarbID,                  &
         plantcarbID, tID, landID, patchID ! dimension IDs
    INTEGER :: latID, lonID, llatvID, llonvID ! time,lat,lon variable ID
    INTEGER :: xvID, yvID   ! coordinate variable IDs for GrADS readability
    !    INTEGER :: surffracID         ! surface fraction varaible ID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file

    TYPE(output_par_settings_type) :: out_settings

    out_settings = output_par_settings_type(met=met, restart=.FALSE.)

    ! Create output file:
    ok = NF90_CREATE(filename%out, NF90_CLOBBER, ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating output file '       &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ! Put the file in define mode:
    ok = NF90_REDEF(ncid_out)
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_out, 'x', xdimsize, xID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining x dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'y', ydimsize, yID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining y dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Define patch dimension, whether it's used or not:
    ok = NF90_DEF_DIM(ncid_out, 'patch', max_vegpatches, patchID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok,'Error defining patch dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    !   ! Define surftype dimension (currently only used for surffrac variable):
    !    ok = NF90_DEF_DIM(ncid_out,'surftype',4,surftypeID)
    !    IF (ok /= NF90_NOERR) CALL nc_abort &
    !         (ok,'Error defining syrftype dimension in output file. '// &
    !         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining vertical soil dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'rad', nrb, radID) ! number of radiation bands
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining radiation dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'soil_carbon_pools', ncs, soilcarbID) ! # pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining soil carbon pool dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out,'plant_carbon_pools',ncp,plantcarbID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok,'Error defining plant carbon pool dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'time', NF90_UNLIMITED, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok,'Error defining time dimension in output file. '// &
         '(SUBROUTINE open_output_file)')
    IF(output%grid == 'mask' .OR. output%grid == 'ALMA' .OR.                   &
         (metGrid == 'mask' .AND. output%grid == 'default')) THEN
       ! for land/sea mask type grid:
       ! Atmospheric 'z' dim of size 1 to comply with ALMA grid type:
       ok = NF90_DEF_DIM(ncid_out, 'z', 1, zID)
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining z dimension in output file. '// &
            '(SUBROUTINE open_output_file)')
    ELSE IF(output%grid == 'land' .OR.                                         &
         (metGrid == 'land' .AND. output%grid == 'default')) THEN
       ! For land only compression grid:
       ok = NF90_DEF_DIM(ncid_out, 'land', mland, landID) ! number of land
       ! points
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining land dimension in output file. '// &
            '(SUBROUTINE open_output_file)')

       ok = NF90_DEF_VAR(ncid_out, 'local_lat', NF90_FLOAT, (/landID/), llatvID)
       IF (ok /= NF90_NOERR) CALL nc_abort                                        &
            (ok, 'Error defining land lat variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       ok = NF90_PUT_ATT(ncid_out, llatvID, 'units', "degrees_north")
       IF (ok /= NF90_NOERR) CALL nc_abort                                        &
            (ok, 'Error defining local lat variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ok = NF90_DEF_VAR(ncid_out, 'local_lon', NF90_FLOAT, (/landID/), llonvID)
       IF (ok /= NF90_NOERR) CALL nc_abort                                        &
            (ok, 'Error defining land lon variable in output file. '// &
            '(SUBROUTINE open_output_file)')
       ok = NF90_PUT_ATT(ncid_out, llonvID, 'units', "degrees_east")
       IF (ok /= NF90_NOERR) CALL nc_abort                                        &
            (ok, 'Error defining local lon variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')

    END IF
    ! Define "time" variable and its attributes:
    ok = NF90_DEF_VAR(ncid_out, 'time', NF90_DOUBLE, (/tID/), ovid%tvar)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'calendar', calendar)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Define latitude and longitude variable (ALMA):
    ok = NF90_DEF_VAR(ncid_out, 'latitude', NF90_FLOAT, (/xID, yID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, latID, 'units', 'degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_VAR(ncid_out, 'longitude', NF90_FLOAT, (/xID, yID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining longitude variable in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining longitude variable attributes in output file. '// &
         '(SUBROUTINE open_output_file)')
    ! Write "cordinate variables" to enable reading by GrADS:
    ok = NF90_DEF_VAR(ncid_out, 'x', NF90_FLOAT, (/xID/), xvID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining "x" variable (for GrADS) in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, xvID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error writing x coordinate variable (GrADS) units in output '// &
         'file. (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, xvID, 'comment',                               &
         'x coordinate variable for GrADS compatibility')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error writing x variables comment in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_VAR(ncid_out, 'y', NF90_FLOAT, (/yID/), yvID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining "y" variable (for GrADS) in output file. '// &
         '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, yvID, 'units', 'degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error writing y coordinate variable (GrADS) units in output '//  &
         'file. (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, yvID, 'comment',                               &
         'y coordinate variable for GrADS compatibility')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error writing y variables comment in output file. '// &
         '(SUBROUTINE open_output_file)')
    !   ! Define fraction of each surface type:
    !   CALL define_ovar(ncid_out,surffracID,'surffrac','-', &
    !       'Fraction of each surface type: vegetated; urban; lake; land ice', &
    !       .FALSE.,surftypeID,'surftype',xID,yID,zID,landID,patchID)

    !=============DEFINE OUTPUT VARIABLES=======================================
    ! Define met forcing variables in output file and allocate temp output vars:
    
    IF(output%SWdown) THEN
       CALL define_ovar(ncid_out,                                              &
            ovid%SWdown, 'SWdown', 'W/m^2', 'Downward shortwave radiation',    &
            patchout%SWdown, 'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SWdown(mp))
       out%SWdown = 0.0 ! initialise
    END IF
    IF(output%LWdown) THEN
       CALL define_ovar(ncid_out, ovid%LWdown, 'LWdown', 'W/m^2',              &
            'Downward longwave radiation', patchout%LWdown,        &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LWdown(mp))
       out%LWdown = 0.0 ! initialise
    END IF
    IF(output%Tair) THEN
       CALL define_ovar(ncid_out, ovid%Tair,                                   &
            'Tair', 'K', 'Surface air temperature', patchout%Tair, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Tair(mp))
       out%Tair = 0.0 ! initialise
    END IF
    IF(output%Rainf) THEN
       CALL define_ovar(ncid_out, ovid%Rainf, 'Rainf',                         &
            'kg/m^2/s', 'Rainfall+snowfall', patchout%Rainf,       &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Rainf(mp))
       out%Rainf = 0.0 ! initialise
    END IF
    IF(output%Snowf) THEN
       CALL define_ovar(ncid_out, ovid%Snowf, 'Snowf',                         &
            'kg/m^2/s', 'Snowfall', patchout%Snowf,                &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Snowf(mp))
       out%Snowf = 0.0 ! initialise
    END IF
    IF(output%Qair) THEN
       CALL define_ovar(ncid_out, ovid%Qair, 'Qair',                           &
            'kg/kg', 'Surface specific humidity', patchout%Qair,   &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qair(mp))
       out%Qair = 0.0 ! initialise
    END IF
    IF(output%Wind) THEN
       CALL define_ovar(ncid_out, ovid%Wind, 'Wind',                           &
            'm/s', 'Scalar surface wind speed', patchout%Wind,     &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Wind(mp))
       out%Wind = 0.0 ! initialise
    END IF
    IF(output%PSurf) THEN
       CALL define_ovar(ncid_out, ovid%PSurf, 'PSurf',                         &
            'hPa', 'Surface air pressure', patchout%PSurf,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PSurf(mp))
       out%PSurf = 0.0 ! initialise
    END IF
    IF(output%CO2air) THEN
       CALL define_ovar(ncid_out, ovid%CO2air, 'CO2air', 'ppmv',               &
            'Surface air CO2 concentration', patchout%CO2air,      &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CO2air(mp))
       out%CO2air = 0.0 ! initialise
    END IF
    ! Define surface flux variables in output file and allocate temp output
    ! vars:
    IF(output%Qmom) THEN
       CALL define_ovar(ncid_out, ovid%Qmom, 'Qmom', 'kg/m/s2',                &
            'Surface momentum flux',patchout%Qmom,'dummy',       &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qmom(mp))
       out%Qmom = 0.0 ! initialise
    END IF
    IF(output%Qle) THEN
       CALL define_ovar(ncid_out, ovid%Qle, 'Qle', 'W/m^2',                    &
            'Surface latent heat flux',patchout%Qle,'dummy',       &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qle(mp))
       out%Qle = 0.0 ! initialise
    END IF
    IF(output%Qh) THEN
       CALL define_ovar(ncid_out,ovid%Qh,'Qh', 'W/m^2',                        &
            'Surface sensible heat flux',patchout%Qh,'dummy',      &
            xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%Qh(mp))
       out%Qh = 0.0 ! initialise
    END IF

    IF(output%Qg) THEN
       CALL define_ovar(ncid_out, ovid%Qg, 'Qg', 'W/m^2',                      &
            'Surface ground heat flux', patchout%Qg, 'dummy',      &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qg(mp))
       out%Qg = 0.0 ! initialise
    END IF
    IF(output%Qs) THEN
       CALL define_ovar(ncid_out, ovid%Qs, 'Qs',                               &
            'kg/m^2/s', 'Surface runoff', patchout%Qs, 'dummy',    &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qs(mp))
       out%Qs = 0.0 ! initialise
    END IF
    IF(output%Qsb) THEN
       CALL define_ovar(ncid_out, ovid%Qsb, 'Qsb', 'kg/m^2/s',                 &
            'Subsurface runoff', patchout%Qsb, 'dummy',            &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qsb(mp))
       out%Qsb = 0.0 ! initialise
    END IF
    IF(output%Evap) THEN
       CALL define_ovar(ncid_out, ovid%Evap,'Evap', 'kg/m^2/s',                &
            'Total evapotranspiration', patchout%Evap, 'dummy',    &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Evap(mp))
       out%Evap = 0.0 ! initialise
    END IF
    IF(output%PotEvap) THEN
       CALL define_ovar(ncid_out, ovid%PotEvap,'PotEvap', 'kg/m^2/s',          &
            'Potential evaporation', patchout%PotEvap, 'dummy',     &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PotEvap(mp))
       out%PotEvap = 0.0 ! initialise
    END IF
    IF(output%ECanop) THEN
       CALL define_ovar(ncid_out, ovid%Ecanop, 'ECanop', 'kg/m^2/s',           &
            'Wet canopy evaporation', patchout%ECanop, 'dummy',    &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ECanop(mp))
       out%ECanop = 0.0 ! initialise
    END IF
    IF(output%TVeg) THEN
       CALL define_ovar(ncid_out, ovid%TVeg, 'TVeg', 'kg/m^2/s',               &
            'Vegetation transpiration', patchout%TVeg, 'dummy',    &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TVeg(mp))
       out%TVeg = 0.0 ! initialise
    END IF
    IF(output%ESoil) THEN
       CALL define_ovar(ncid_out, ovid%ESoil, 'ESoil', 'kg/m^2/s',             &
            'Evaporation from soil', patchout%ESoil, 'dummy',      &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ESoil(mp))
       out%ESoil = 0.0 ! initialise
    END IF
    IF(output%HVeg) THEN
       CALL define_ovar(ncid_out, ovid%HVeg, 'HVeg', 'W/m^2',                  &
            'Sensible heat from vegetation', patchout%HVeg,        &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HVeg(mp))
       out%HVeg = 0.0 ! initialise
    END IF
    IF(output%HSoil) THEN
       CALL define_ovar(ncid_out, ovid%HSoil, 'HSoil', 'W/m^2',                &
            'Sensible heat from soil', patchout%HSoil, 'dummy',    &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HSoil(mp))
       out%HSoil = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%RnetSoil, 'RnetSoil', 'W/m^2',                &
            'Net radiation absorbed by ground', patchout%RnetSoil, 'dummy',    &
            xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RnetSoil(mp))
       out%RnetSoil = 0.0 ! initialise
    END IF
    IF(output%NEE) THEN
       CALL define_ovar(ncid_out, ovid%NEE, 'NEE', 'umol/m^2/s',               &
            'Net ecosystem exchange of CO2', patchout%NEE,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NEE(mp))
       out%NEE = 0.0 ! initialise
    END IF



    ! Define soil state variables in output file and allocate temp output vars:
    IF(output%SoilMoist) THEN
       CALL define_ovar(ncid_out, ovid%SoilMoist, 'SoilMoist', 'm^3/m^3',      &
            'Average layer soil moisture', patchout%SoilMoist,     &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       CALL define_ovar(ncid_out, ovid%SoilMoistIce, 'SoilMoistIce', 'm^3/m^3',      &
            'Average layer frozen soil moisture', patchout%SoilMoistIce,     &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilMoist(mp,ms))
       ALLOCATE(out%SoilMoistIce(mp,ms))
       out%SoilMoist = 0.0 ! initialise
       out%SoilMoistIce = 0.0 ! initialise
    END IF
    IF(output%SoilTemp) THEN
       CALL define_ovar(ncid_out, ovid%SoilTemp, 'SoilTemp', 'K',              &
            'Average layer soil temperature', patchout%SoilTemp,   &
            'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilTemp(mp,ms))
       out%SoilTemp = 0.0 ! initialise
    END IF
    IF(output%BaresoilT) THEN
       CALL define_ovar(ncid_out, ovid%BaresoilT, 'BaresoilT',                 &
            'K', 'Bare soil temperature', patchout%BaresoilT,      &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%BaresoilT(mp))
       out%BaresoilT = 0.0 ! initialise
    END IF
    ! Define snow state variables in output file and allocate temp output vars:
    IF(output%SWE) THEN
       CALL define_ovar(ncid_out, ovid%SWE, 'SWE', 'kg/m^2',                   &
            'Snow water equivalent', patchout%SWE,                 &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SWE(mp))
       out%SWE = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%SnowMelt, 'SnowMelt', 'kg/m^2/s',                   &
            'Snow Melt Rate', patchout%SnowMelt,                 &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowMelt(mp))
       out%SnowMelt = 0.0 ! initialise
    END IF
    IF(output%SnowT) THEN
       CALL define_ovar(ncid_out, ovid%SnowT, 'SnowT', 'K',                    &
            'Snow surface temperature', patchout%SnowT,            &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowT(mp))
       out%SnowT = 0.0 ! initialise
    END IF
    IF(output%SnowDepth) THEN
       CALL define_ovar(ncid_out, ovid%SnowDepth, 'SnowDepth',                 &
            'm', 'Snow depth', patchout%SnowDepth,                 &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowDepth(mp))
       out%SnowDepth = 0.0 ! initialise
    END IF
    ! Define radiative variables in output file and allocate temp output vars:
    IF(output%SWnet) THEN
       CALL define_ovar(ncid_out, ovid%SWnet, 'SWnet', 'W/m^2',                &
            'Net shortwave radiation absorbed by surface',         &
            patchout%SWnet, 'dummy', xID, yID, zID, landID,       &
            patchID, tID)
       ALLOCATE(out%SWnet(mp))
       out%SWnet = 0.0 ! initialise
    END IF
    IF(output%LWnet) THEN
       CALL define_ovar(ncid_out, ovid%LWnet, 'LWnet', 'W/m^2',                &
            'Net longwave radiation absorbed by surface',          &
            patchout%LWnet, 'dummy', xID, yID, zID, landID,        &
            patchID, tID)
       ALLOCATE(out%LWnet(mp))
       out%LWnet = 0.0 ! initialise
    END IF
    IF(output%Rnet) THEN
       CALL define_ovar(ncid_out, ovid%Rnet, 'Rnet', 'W/m^2',                  &
            'Net radiation absorbed by surface', patchout%Rnet,    &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Rnet(mp))
       out%Rnet = 0.0 ! initialise
    END IF
    IF(output%Albedo) THEN
       CALL define_ovar(ncid_out, ovid%Albedo, 'Albedo', '-',                  &
            'Surface albedo', patchout%Albedo,                     &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Albedo(mp))
       out%Albedo = 0.0 ! initialise
    END IF

    ! output calc of soil albedo based on colour? - Ticket #27
    IF (calcsoilalbedo) THEN
       IF(output%visAlbedo) THEN
          CALL define_ovar(ncid_out, ovid%visAlbedo, 'visAlbedo', '-',          &
               'Surface vis albedo', patchout%visAlbedo,              &
               'dummy', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%visAlbedo(mp))
          out%visAlbedo = 0.0 ! initialise
       END IF
       IF(output%nirAlbedo) THEN
          CALL define_ovar(ncid_out, ovid%nirAlbedo, 'nirAlbedo', '-',          &
               'Surface nir albedo', patchout%nirAlbedo,              &
               'dummy', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%nirAlbedo(mp))
          out%nirAlbedo = 0.0 ! initialise
       END IF
    END IF

    IF(output%RadT) THEN
       CALL define_ovar(ncid_out, ovid%RadT, 'RadT', 'K',                      &
            'Radiative surface temperature', patchout%RadT,        &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RadT(mp))
       out%RadT = 0.0 ! initialise
    END IF
    ! Define vegetation variables in output file and allocate temp output vars:
    ! REV_CORR - new output variables.
    IF(output%Tscrn) THEN
       CALL define_ovar(ncid_out, ovid%Tscrn,                                  &
            'Tscrn', 'oC', 'screen level air temperature', &
            patchout%Tscrn, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Tscrn(mp))
       out%Tscrn = 0.0 ! initialise
    END IF
    IF (output%Tex) THEN
       IF((output%averaging(1:2) == 'da').OR.(output%averaging(1:2)=='mo')) THEN
          CALL define_ovar(ncid_out, ovid%Txx,                                 &
               'Txx', 'oC', 'max screen-level T in reporting period',&
               patchout%Tex,                                         &
               'ALMA', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%Txx(mp))
          out%Txx = -1.0E6 !initialise extremes at unreasonable value
          CALL define_ovar(ncid_out, ovid%Tnn,                                &
               'Tnn', 'oC', 'min screen-level T in reporting period',&
               patchout%Tex,                                         &
               'ALMA', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%Tnn(mp))
          out%Tnn = 1.0E6 !initialise extremes at unreasonable value
       ENDIF
       IF (output%averaging(1:2)=='mo') THEN
          !%Tdaymx is the current day max T - this is a working variable)
          !%Tmx is average of those values to be output
          CALL define_ovar(ncid_out, ovid%Tmx,                                &
               'Tmx', 'oC', 'averaged daily maximum screen-level T', &
               patchout%Tex,                                         &
               'ALMA', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%Tmx(mp), out%Tdaymx(mp))
          out%Tmx = 0.0 !initialise average
          out%Tdaymx = -1.0E6 !initialise extremes at unreasonable value
          CALL define_ovar(ncid_out, ovid%Tmn,                                &
               'Tmn', 'oC', 'averaged daily minimum screen-level T', &
               patchout%Tex,                                         &
               'ALMA', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%Tmn(mp),out%Tdaymn(mp))
          out%Tmn = 0.0
          out%Tdaymn = 1.0E6
       ENDIF
    ENDIF
    IF(output%Qscrn) THEN
       CALL define_ovar(ncid_out, ovid%Qscrn,                                  &
            'Qscrn', 'kg/kg', 'screen level specific humdity',     &
            patchout%Qscrn, &
            'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qscrn(mp))
       out%Qscrn = 0.0 ! initialise
    END IF
    IF(output%VegT) THEN
       CALL define_ovar(ncid_out, ovid%VegT, 'VegT', 'K',                      &
            'Average vegetation temperature', patchout%VegT,       &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%VegT(mp))
       out%VegT = 0.0 ! initialise
    END IF
    IF(output%CanT) THEN
       CALL define_ovar(ncid_out, ovid%CanT, 'CanT', 'K', &
            'Within-canopy temperature', patchout%CanT, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CanT(mp))
       out%CanT = 0.0 ! initialise
    END IF
    IF(output%Fwsoil) THEN
       CALL define_ovar(ncid_out, ovid%Fwsoil, 'Fwsoil', '[-]', &
            'soil moisture modifier to stomatal conductance', patchout%Fwsoil, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Fwsoil(mp))
       out%Fwsoil = 0.0 ! initialise
    END IF
    IF(output%CanopInt) THEN
       CALL define_ovar(ncid_out, ovid%CanopInt, 'CanopInt', 'kg/m^2',         &
            'Canopy intercepted water storage', patchout%CanopInt, &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CanopInt(mp))
       out%CanopInt = 0.0 ! initialise
    END IF
    IF(output%LAI) THEN
       CALL define_ovar(ncid_out, ovid%LAI, 'LAI', '-',                        &
            'Leaf area index', patchout%LAI, 'dummy', xID,         &
            yID, zID, landID, patchID, tID)
       ALLOCATE(out%LAI(mp))
       out%LAI = 0.0 ! initialise
    END IF
    ! Define balance variables in output file and allocate temp output vars:
    IF(output%Ebal) THEN
       CALL define_ovar(ncid_out, ovid%Ebal, 'Ebal', 'W/m^2',                  &
            'Cumulative energy imbalance', patchout%Ebal,          &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Ebal(mp))
       out%Ebal = 0.0 ! initialise
    END IF
    IF(output%Wbal) THEN
       CALL define_ovar(ncid_out, ovid%Wbal, 'Wbal', 'kg/m^2',                 &
            'Cumulative water imbalance', patchout%Wbal,           &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Wbal(mp))
       out%Wbal = 0.0 ! initialise
    END IF
    ! Define carbon variables in output file and allocate temp output vars:
    IF(output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%AutoResp, 'AutoResp', 'umol/m^2/s',     &
            'Autotrophic respiration', patchout%AutoResp,          &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%AutoResp(mp))
       out%AutoResp = 0.0 ! initialise
    END IF
    IF(output%casa .AND. output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%RootResp, 'RootResp', 'umol/m^2/s',     &
            'Fine Root Autotrophic respiration', patchout%AutoResp,          &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RootResp(mp))
       out%RootResp = 0.0 ! initialise
    END IF

    IF(output%casa .AND. output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%StemResp, 'StemResp', 'umol/m^2/s',     &
            'StemWood Autotrophic respiration', patchout%AutoResp,          &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%StemResp(mp))
       out%StemResp = 0.0 ! initialise
    END IF

    IF(output%LeafResp) THEN
       CALL define_ovar(ncid_out, ovid%LeafResp, 'LeafResp', 'umol/m^2/s',     &
            'Leaf respiration', patchout%LeafResp,                 &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LeafResp(mp))
       out%LeafResp = 0.0 ! initialise
    END IF
    IF(output%HeteroResp) THEN
       CALL define_ovar(ncid_out, ovid%HeteroResp, 'HeteroResp', 'umol/m^2/s', &
            'Heterotrophic respiration', patchout%HeteroResp,      &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HeteroResp(mp))
       out%HeteroResp = 0.0 ! initialise
    END IF
    IF(output%GPP) THEN
       CALL define_ovar(ncid_out, ovid%GPP, 'GPP', 'umol/m^2/s',               &
            'Gross primary production', patchout%GPP,              &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP(mp))
       out%GPP = 0.0 ! initialise


    END IF





    IF(output%NPP) THEN
       CALL define_ovar(ncid_out, ovid%NPP, 'NPP', 'umol/m^2/s',               &
            'Net primary production', patchout%NPP,                &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NPP(mp))
       out%NPP = 0.0 ! initialise
    END IF

    !MD groundwater related variables
    IF(output%WatTable) THEN
       CALL define_ovar(ncid_out, ovid%WatTable, 'WatTable', 'm',      &
            'Water Table Depth', patchout%WatTable,     &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%WatTable(mp))
       out%WatTable = 0.0 ! initialise
    END IF
    IF(output%GWMoist) THEN
       CALL define_ovar(ncid_out, ovid%GWMoist, 'GWMoist', 'mm3/mm3',      &
            'Aquifer mositure content', patchout%GWMoist,     &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GWMoist(mp))
       out%GWMoist = 0.0 ! initialise
    END IF

    IF(output%SatFrac) THEN
       CALL define_ovar(ncid_out, ovid%SatFrac, 'SatFrac', 'unitless',      &
            'Saturated Fraction of Gridcell', patchout%SatFrac,     &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SatFrac(mp))
       out%SatFrac = 0.0 ! initialise
    END IF

    IF(output%Qrecharge) THEN
       CALL define_ovar(ncid_out, ovid%Qrecharge, 'Qrecharge', 'mm/s',      &
            'Recharge to or from Aquifer', patchout%Qrecharge,     &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qrecharge(mp))
       out%Qrecharge = 0.0 ! initialise
    END IF

    IF(output%casa) THEN
       CALL define_ovar(ncid_out, ovid%NBP, 'NBP', 'umol/m^2/s',               &
            'Net Biosphere Production (uptake +ve)', patchout%NBP,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NBP(mp))
       out%NBP = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%dCdt, 'dCdt', 'umol/m^2/s',               &
            'Carbon accumulation rate (uptake +ve)', patchout%dCdt,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%dCdt(mp))
       out%dCdt = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%TotSoilCarb, 'TotSoilCarb', 'kg C/m^2',               &
            'Total Soil and Litter Carbon', patchout%TotSoilCarb,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TotSoilCarb(mp))
       out%TotSoilCarb = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%TotLittCarb, 'TotLittCarb', 'kg C/m^2',               &
            'Total Litter Carbon', patchout%TotLittCarb,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TotLittCarb(mp))
       out%TotLittCarb = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%SoilCarbFast, 'SoilCarbFast', 'kg C/m^2',               &
            'Soil Carbon: Fast Turnover', patchout%SoilCarbFast,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SoilCarbFast(mp))
       out%SoilCarbFast = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%SoilCarbSlow, 'SoilCarbSlow', 'kg C/m^2',               &
            'Soil Carbon: Slow Turnover', patchout%SoilCarbSlow,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SoilCarbSlow(mp))
       out%SoilCarbSlow = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%SoilCarbPassive, 'SoilCarbPassive', 'kg C/m^2',               &
            'Soil Carbon: Passive', patchout%SoilCarbPassive,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SoilCarbPassive(mp))
       out%SoilCarbPassive = 0.0 ! initialise


       CALL define_ovar(ncid_out, ovid%LittCarbMetabolic, 'LittCarbMetabolic', 'kg C/m^2',               &
            'Litter Carbon: metabolic', patchout%LittCarbMetabolic,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LittCarbMetabolic(mp))
       out%LittCarbMetabolic = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%LittCarbStructural, 'LittCarbStructural', 'kg C/m^2',               &
            'Litter Carbon: structural', patchout%LittCarbStructural,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LittCarbStructural(mp))
       out%LittCarbStructural = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%LittCarbCWD, 'LittCarbCWD', 'kg C/m^2',               &
            'Litter Carbon: CWD', patchout%LittCarbCWD,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LittCarbCWD(mp))
       out%LittCarbCWD = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%PlantCarbLeaf, 'PlantCarbLeaf', 'kg C/m^2',               &
            'Plant Carbon: leaf', patchout%PlantCarbLeaf,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantCarbLeaf(mp))
       out%PlantCarbLeaf = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%PlantCarbFineRoot, 'PlantCarbFineRoot', 'kg C/m^2',               &
            'Plant Carbon: Fine roots', patchout%PlantCarbFineRoot,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantCarbFineRoot(mp))
       out%PlantCarbFineRoot = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%PlantCarbWood, 'PlantCarbWood', 'kg C/m^2',               &
            'Plant Carbon: wood (above- and below-ground', patchout%PlantCarbWood,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantCarbWood(mp))
       out%PlantCarbWood = 0.0 ! initialise

       CALL define_ovar(ncid_out, ovid%TotLivBiomass, 'TotLivBiomass', 'kg C/m^2',               &
            'Total Biomass', patchout%TotLivBiomass,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TotLivBiomass(mp))
       out%TotLivBiomass = 0.0 ! initialise


       CALL define_ovar(ncid_out, ovid%PlantTurnover, 'PlantTurnover', 'umol/m^2/s',               &
            'Total Biomass Turnover', patchout%PlantTurnover,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnover(mp))
       out%PlantTurnover = 0.0

       CALL define_ovar(ncid_out, ovid%PlantTurnoverLeaf, 'PlantTurnoverLeaf ', &
            'umol/m^2/s',               &
            'Leaf Biomass Turnover', patchout%PlantTurnoverLeaf,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverLeaf(mp))
       out%PlantTurnoverLeaf = 0.0

       CALL define_ovar(ncid_out, ovid%PlantTurnoverFineRoot, 'PlantTurnoverFineRoot ', &
            'umol/m^2/s',               &
            'FineRoot Biomass Turnover', patchout%PlantTurnoverFineRoot,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverFineRoot(mp))
       out%PlantTurnoverFineRoot = 0.0


       CALL define_ovar(ncid_out, ovid%PlantTurnoverWood, 'PlantTurnoverWood ', &
            'umol/m^2/s',               &
            'Woody Biomass Turnover', patchout%PlantTurnoverWood,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWood(mp))
       out%PlantTurnoverWood = 0.0

       CALL define_ovar(ncid_out, ovid%PlantTurnoverWoodDist, 'PlantTurnoverWoodDist ', &
            'umol/m^2/s',               &
            'Woody Biomass Turnover (disturbance)', patchout%PlantTurnoverWoodDist,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWoodDist(mp))
       out%PlantTurnoverWoodDist = 0.0


       CALL define_ovar(ncid_out, ovid%PlantTurnoverWoodCrowding, 'PlantTurnoverWoodCrowding ', &
            'umol/m^2/s',               &
            'Woody Biomass Turnover (crowding)', patchout%PlantTurnoverWoodCrowding,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWoodCrowding(mp))
       out%PlantTurnoverWoodCrowding = 0.0

       CALL define_ovar(ncid_out, ovid%PlantTurnoverWoodResourceLim, 'PlantTurnoverWoodResourceLim ', &
            'umol/m^2/s',               &
            'Woody Biomass Turnover (Resource Limitation)', patchout%PlantTurnoverWoodResourceLim,         &
            'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PlantTurnoverWoodResourceLim(mp))
       out%PlantTurnoverWoodResourceLim = 0.0

       IF (cable_user%POPLUC) THEN

          CALL define_ovar(ncid_out, ovid%LandUseFlux, 'LandUseFlux ', &
               'umol/m^2/s',               &
               'Sum of wood harvest and clearing fluxes', patchout%LandUseFlux,         &
               'dummy', xID, yID, zID, landID, patchID, tID)
          ALLOCATE(out%LandUseFlux(mp))
          out%LandUseFlux = 0.0
       ENDIF


    END IF

    ! vh_js !
    CALL define_ovar(ncid_out, ovid%Area, 'Area', 'km2',               &
         'Patch Area', patchout%Area,         &
         'dummy', xID, yID, zID, landID, patchID, tID)
    ALLOCATE(out%Area(mp))
    out%Area = 0.0 ! initialise


    ! Define CABLE parameters in output file:
    IF(output%iveg) CALL define_ovar(ncid_out, opid%iveg,   &
         'iveg', '-', 'Vegetation type', patchout%iveg, 'integer', &
         xID, yID, zID, landID, patchID)

    IF (cable_user%POPLUC) THEN

       CALL define_ovar(ncid_out, opid%patchfrac, 'patchfrac', '-',          &
            'Fractional cover of vegetation patches', patchout%patchfrac, 'real', &
            xID, yID, zID, landID, patchID, tID)

    ELSE

       IF((output%patchfrac)                                   &
            .AND. (patchout%patchfrac .OR. output%patch))                         &
            CALL define_ovar(ncid_out, opid%patchfrac, 'patchfrac', '-',          &
            'Fractional cover of vegetation patches', patchout%patchfrac, 'real', &
            xID, yID, zID, landID, patchID)

    ENDIF


    IF(output%isoil) CALL define_ovar(ncid_out, opid%isoil, &
         'isoil', '-', 'Soil type', patchout%isoil, 'integer', &
         xID, yID, zID, landID, patchID)
    IF(output%bch) CALL define_ovar(ncid_out, opid%bch,     &
         'bch', '-', 'Parameter b, Campbell eqn 1985', patchout%bch, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%clay) CALL define_ovar(ncid_out, opid%clay,   &
         'clay', '-', 'Fraction of soil which is clay', patchout%clay, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%sand) CALL define_ovar(ncid_out, opid%sand,   &
         'sand', '-', 'Fraction of soil which is sand', patchout%sand, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%silt) CALL define_ovar(ncid_out, opid%silt,   &
         'silt', '-', 'Fraction of soil which is silt', patchout%silt, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%ssat) CALL define_ovar(ncid_out, opid%ssat,   &
         'ssat', '-', 'Fraction of soil volume which is water @ saturation', &
         patchout%ssat, 'real', xID, yID, zID, landID, patchID)
    IF(output%sfc) CALL define_ovar(ncid_out, opid%sfc,     &
         'sfc', '-', 'Fraction of soil volume which is water @ field capacity', &
         patchout%sfc, 'real', xID, yID, zID, landID, patchID)
    IF(output%swilt) CALL define_ovar(ncid_out, opid%swilt, &
         'swilt', '-', 'Fraction of soil volume which is water @ wilting point', &
         patchout%swilt, 'real', xID, yID, zID, landID, patchID)
    IF(output%hyds) CALL define_ovar(ncid_out, opid%hyds,   &
         'hyds', 'm/s', 'Hydraulic conductivity @ saturation', &
         patchout%hyds, 'real', xID, yID, zID, landID, patchID)
    IF(output%sucs) CALL define_ovar(ncid_out, opid%sucs,   &
         'sucs', 'm', 'Suction @ saturation', &
         patchout%sucs, 'real', xID, yID, zID, landID, patchID)
    IF(output%css) CALL define_ovar(ncid_out, opid%css,     &
         'css', 'J/kg/C', 'Heat capacity of soil minerals', &
         patchout%css, 'real', xID, yID, zID, landID, patchID)
    IF(output%rhosoil) CALL define_ovar(ncid_out,           &
         opid%rhosoil, 'rhosoil', 'kg/m^3', 'Density of soil minerals', &
         patchout%rhosoil, 'real', xID, yID, zID, landID, patchID)
    IF(output%rs20) CALL define_ovar(ncid_out, opid%rs20,   &
         'rs20', '-', 'Soil respiration coefficient at 20C', &
         patchout%rs20, 'real', xID, yID, zID, landID, patchID)
    IF(output%albsoil) CALL define_ovar(ncid_out,           &
         opid%albsoil, 'albsoil', '-', &
         'Snow free shortwave soil reflectance fraction', &
         patchout%albsoil, radID, 'radiation', xID, yID, zID, landID, patchID)
    ! vh_js !
    IF (cable_user%CALL_POP) THEN
       IF(output%hc) CALL define_ovar(ncid_out, opid%hc,    &
            'hc', 'm', 'Height of canopy', patchout%hc, &
            'real', xID, yID, zID, landID, patchID,tID)
    ELSE
       IF(output%hc) CALL define_ovar(ncid_out, opid%hc,    &
            'hc', 'm', 'Height of canopy', patchout%hc, &
            'real', xID, yID, zID, landID, patchID)
    ENDIF

    IF(output%canst1) CALL define_ovar(ncid_out,            &
         opid%canst1, 'canst1', 'mm/LAI', 'Max water intercepted by canopy', &
         patchout%canst1, 'real', xID, yID, zID, landID, patchID)
    IF(output%dleaf) CALL define_ovar(ncid_out, opid%dleaf, &
         'dleaf', 'm', 'Chararacteristic length of leaf', &
         patchout%dleaf, 'real', xID, yID, zID, landID, patchID)
    IF(output%frac4) CALL define_ovar(ncid_out, opid%frac4, &
         'frac4', '-', 'Fraction of plants which are C4', &
         patchout%frac4, 'real', xID, yID, zID, landID, patchID)
    IF(output%ejmax) CALL define_ovar(ncid_out, opid%ejmax, &
         'ejmax', 'mol/m^2/s', 'Max potential electron transport rate top leaf', &
         patchout%ejmax, 'real', xID, yID, zID, landID, patchID)
    IF(output%vcmax) CALL define_ovar(ncid_out, opid%vcmax, &
         'vcmax', 'mol/m^2/s', 'Maximum RuBP carboxylation rate top leaf', &
         patchout%vcmax, 'real', xID, yID, zID, landID, patchID)
    IF(output%rp20) CALL define_ovar(ncid_out, opid%rp20,   &
         'rp20', '-', 'Plant respiration coefficient at 20C', &
         patchout%rp20, 'real', xID, yID, zID, landID, patchID)
    ! Ticket #56
    IF(output%g0) CALL define_ovar(ncid_out, opid%g0,   &
         'g0', '-', 'g0 term in Medlyn Stom Cond. Param', &
         patchout%g0, 'real', xID, yID, zID, landID, patchID)
    IF(output%g1) CALL define_ovar(ncid_out, opid%g1,   &
         'g1', '-', 'g1 term in Medlyn Stom Cond. Param', &
         patchout%g1, 'real', xID, yID, zID, landID, patchID)
    ! end Ticket #56

    IF(output%rpcoef) CALL define_ovar(ncid_out,            &
         opid%rpcoef, 'rpcoef', '1/C', &
         'Temperature coef nonleaf plant respiration', &
         patchout%rpcoef, 'real', xID, yID, zID, landID, patchID)
    IF(output%shelrb) CALL define_ovar(ncid_out,            &
         opid%shelrb, 'shelrb', '-', 'Sheltering factor', patchout%shelrb, &
         'real', xID, yID, zID, landID, patchID)
    IF(output%xfang) CALL define_ovar(ncid_out, opid%xfang, &
         'xfang', '-', 'Leaf angle parameter',patchout%xfang, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%wai) CALL define_ovar(ncid_out, opid%wai,     &
         'wai', '-', 'Wood area index', patchout%wai, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%vegcf) CALL define_ovar(ncid_out, opid%vegcf, &
         'vegcf', '-', 'vegcf', patchout%vegcf, 'real', &
         xID, yID, zID, landID, patchID)
    IF(output%extkn) CALL define_ovar(ncid_out, opid%extkn, &
         'extkn', '-', 'Nitrogen extinction coef for vert. canopy profile', &
         patchout%extkn, 'real', xID, yID, zID, landID, patchID)
    IF(output%tminvj) CALL define_ovar(ncid_out,            &
         opid%tminvj, 'tminvj', 'C', &
         'Min temperature for the start of photosynthesis', &
         patchout%tminvj, 'real', xID, yID, zID, landID, patchID)
    IF(output%tmaxvj) CALL define_ovar(ncid_out,            &
         opid%tmaxvj, 'tmaxvj', 'C', 'Max temperature for photosynthesis', &
         patchout%tmaxvj, 'real', xID, yID, zID, landID, patchID)
    IF(output%vbeta) CALL define_ovar(ncid_out, opid%vbeta, &
         'vbeta', '-', 'Stomatal sensitivity to soil water', &
         patchout%vbeta, 'real', xID, yID, zID, landID, patchID)
    IF(output%xalbnir) CALL define_ovar(ncid_out,           &
         opid%xalbnir, 'xalbnir', '-', 'Modifier for albedo in near ir band', &
         patchout%xalbnir, 'real', xID, yID, zID, landID, patchID)
    IF(output%meth) CALL define_ovar(ncid_out, opid%meth,   &
         'meth', '-', 'Canopy turbulence parameterisation choice', &
         patchout%meth, 'real', xID, yID, zID, landID, patchID)
    IF(output%za) THEN
       CALL define_ovar(ncid_out, opid%za_uv, 'za_uv', 'm',                     &
            'Reference height (lowest atm. model layer) for momentum', &
            patchout%za, 'real', xID, yID, zID, landID, patchID)
       CALL define_ovar(ncid_out, opid%za_tq, 'za_tq', 'm',                     &
            'Reference height (lowest atm. model layer) for scalars', &
            patchout%za, 'real', xID, yID, zID, landID, patchID)
    ENDIF
    IF(output%ratecp) CALL define_ovar(ncid_out,            &
         opid%ratecp, 'ratecp', '1/year', 'Plant carbon rate constant', &
         patchout%ratecp, plantcarbID, 'plantcarbon', xID, yID, zID, &
         landID, patchID)
    IF(output%ratecs) CALL define_ovar(ncid_out,            &
         opid%ratecs, 'ratecs', '1/year', 'Soil carbon rate constant', &
         patchout%ratecs, soilcarbID, 'soilcarbon', xID, yID, zID, &
         landID, patchID)
    IF(output%zse) CALL define_ovar(ncid_out, opid%zse,     &
         'zse', 'm', 'Depth of each soil layer', &
         patchout%zse, soilID, 'soil', xID, yID, zID, landID, patchID)
    IF(output%froot) CALL define_ovar(ncid_out, opid%froot, &
         'froot', '-', 'Fraction of roots in each soil layer', &
         patchout%froot, soilID, 'soil', xID, yID, zID, landID, patchID)

    !    IF(output%slope) CALL define_ovar(ncid_out, opid%slope,   &
    !           'slope', '-', 'Mean subgrid topographic slope', &
    !                          patchout%slope, 'real', xID, yID, zID, landID, patchID)
    !
    !    IF(output%slope_std) CALL define_ovar(ncid_out, opid%slope_std,   &
    !           'slope_std', '-', 'Mean subgrid topographic slope_std', &
    !                          patchout%slope_std, 'real', xID, yID, zID, landID, patchID)
    !
    !    IF(output%GWdz) CALL define_ovar(ncid_out, opid%GWdz,   &
    !           'GWdz', '-', 'Mean aquifer layer thickness ', &
    !                          patchout%GWdz, 'real', xID, yID, zID, landID, patchID)
    !
    IF(output%params .AND. cable_user%gw_model) THEN
       CALL define_ovar(ncid_out, opid%Qhmax,   &
            'Qhmax', 'mm/s', 'Maximum subsurface drainage ', &
            patchout%Qhmax, 'real', xID, yID, zID, landID, patchID)
       CALL define_ovar(ncid_out, opid%QhmaxEfold,   &
            'QhmaxEfold', 'm', 'Maximum subsurface drainage decay rate', &
            patchout%QhmaxEfold, 'real', xID, yID, zID, landID, patchID)
       CALL define_ovar(ncid_out, opid%SatFracmax,   &
            'SatFracmax', '-', 'Controls max saturated fraction ', &
            patchout%SatFracmax, 'real', xID, yID, zID, landID, patchID)
       CALL define_ovar(ncid_out, opid%HKefold,   &
            'HKefold', '1/m', 'Rate HK decays with depth ', &
            patchout%HKefold, 'real', xID, yID, zID, landID, patchID)
       CALL define_ovar(ncid_out, opid%HKdepth,   &
            'HKdepth', 'm', 'Depth at which HKsat(z) is HKsat(0) ', &
            patchout%HKdepth, 'real', xID, yID, zID, landID, patchID)
    END IF


    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Production",                     &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out,NF90_GLOBAL,"Source", &
         'CABLE LSM output file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "CABLE_input_file",               &
         TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Determine output aggregation details:
    IF(output%averaging(1:4) == 'user') THEN
       ! User-specified aggregation interval for output:
       READ(output%averaging(5:7), *) output%interval
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
            TRIM(output%averaging(5:7))//'-hourly output')
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Convert interval value from hours to time steps (for use in output
       ! write):
       output%interval = output%interval * 3600/INT(dels)
    ELSE IF(output%averaging(1:3) == 'all') THEN ! output all timesteps
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
            TRIM(output%averaging)//' timesteps recorded')
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Set output interval to be one time step
       output%interval = 1
    ELSE IF(output%averaging(1:2) == 'mo') THEN ! monthly output
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
            TRIM(output%averaging))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Output interval will be determined dynamically for monthly output
    ELSE IF(output%averaging(1:2) == 'da') THEN ! daily output
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
            TRIM(output%averaging))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
            'Error writing global detail to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Set output interval to be # time steps in 24 hours:
       output%interval = 3600*24/INT(dels)
    ELSE
       CALL abort ('Unknown output averaging interval specified '//            &
            'in namelist file. (SUBROUTINE open_output_file)')
    END IF

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating output file '       &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Write latitude and longitude variables:

    ok = NF90_PUT_VAR(ncid_out, latID, REAL(lat_all, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
         (ok, 'Error writing latitude variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_VAR(ncid_out, lonID, REAL(lon_all, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
         (ok, 'Error writing longitude variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    IF (output%grid == 'land' .OR. &
         (metGrid == 'land' .AND. output%grid == 'default'))  THEN

       ok = NF90_PUT_VAR(ncid_out, llatvID, REAL(latitude) )
       IF(ok /= NF90_NOERR) CALL nc_abort                                         &
            (ok, 'Error writing loc lat variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

       ok = NF90_PUT_VAR(ncid_out, llonvID, REAL(longitude) )
       IF(ok /= NF90_NOERR) CALL nc_abort                                         &
            (ok, 'Error writing loc lon variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ENDIF

    ! Write GrADS coordinate variables
    ok = NF90_PUT_VAR(ncid_out, xvID, REAL(lon_all(:, 1), 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
         (ok, 'Error writing GrADS x coordinate variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_VAR(ncid_out, yvID, REAL(lat_all(1, :), 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
         (ok, 'Error writing GrADS y coordinate variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Write model parameters if requested:

    !~ Patch
    out_settings%dimswitch = "real"
    IF (output%patchfrac .AND. (patchout%patchfrac .OR. output%patch)) THEN
      CALL check_and_write(opid%patchfrac, 'patchfrac',                  &
           REAL(patch(:)%frac, 4), ranges%patchfrac, patchout%patchfrac, out_settings)
    END IF

    !~ Soil
    out_settings%dimswitch = "integer"
    IF (output%isoil) THEN
      CALL check_and_write(opid%isoil,  &
           'isoil', REAL(soil%isoilm, 4), ranges%isoil, patchout%isoil, out_settings)
    END IF
    out_settings%dimswitch = "real"
    IF (output%bch) THEN
      CALL check_and_write(opid%bch,      &
           'bch', REAL(soil%bch, 4), ranges%bch, patchout%bch, out_settings)
    END IF
    IF (output%clay) THEN
      CALL check_and_write(opid%clay,    &
           'clay', REAL(soil%clay, 4), ranges%clay, patchout%clay, out_settings)
    END IF
    IF (output%sand) THEN
      CALL check_and_write(opid%sand,    &
           'sand', REAL(soil%sand, 4), ranges%sand, patchout%sand, out_settings)
    END IF
    IF (output%silt) THEN
      CALL check_and_write(opid%silt,    &
           'silt', REAL(soil%silt, 4), ranges%silt, patchout%silt, out_settings)
    END IF
    IF (output%css) THEN
      CALL check_and_write(opid%css,      &
           'css', REAL(soil%css, 4), ranges%css, patchout%css, out_settings)
    END IF
    IF (output%rhosoil) THEN
      CALL check_and_write(opid%rhosoil, 'rhosoil',REAL(soil%rhosoil,4), &
           ranges%rhosoil, patchout%rhosoil, out_settings)
    END IF
    IF (output%hyds) THEN
      CALL check_and_write(opid%hyds,    &
           'hyds', REAL(soil%hyds, 4), ranges%hyds, patchout%hyds, out_settings)
    END IF
    IF (output%sucs) THEN
      CALL check_and_write(opid%sucs,    &
           'sucs', REAL(soil%sucs, 4), ranges%sucs, patchout%sucs, out_settings)
    END IF
    IF (output%rs20) THEN
      CALL check_and_write(opid%rs20,    &
           'rs20', REAL(veg%rs20, 4), ranges%rs20, patchout%rs20, out_settings)
    !           'rs20',REAL(soil%rs20,4),ranges%rs20,patchout%rs20,out_settings)
    END IF
    IF (output%ssat) THEN
      CALL check_and_write(opid%ssat,    &
           'ssat', REAL(soil%ssat, 4), ranges%ssat, patchout%ssat, out_settings)
    END IF
    IF (output%sfc) THEN
      CALL check_and_write(opid%sfc,      &
           'sfc', REAL(soil%sfc, 4), ranges%sfc, patchout%sfc, out_settings)
    END IF
    IF (output%swilt) THEN
      CALL check_and_write(opid%swilt,  &
           'swilt', REAL(soil%swilt, 4), ranges%swilt, patchout%swilt, out_settings)
    END IF

    !    IF (output%slope) THEN
    !      CALL check_and_write(ncid_out, opid%slope,    &
    !                   'slope', REAL(soil%slope, 4), ranges%slope, patchout%slope, out_settings)
    !    END IF
    !    IF (output%slope_std) THEN
    !      CALL check_and_write(opid%slope_std,    &
    !                   'slope_std', REAL(soil%slope_std, 4), ranges%slope_std, patchout%slope_std, out_settings)
    !    END IF
    !    IF (output%GWdz) THEN
    !      CALL check_and_write(opid%GWdz,    &
    !                   'GWdz', REAL(soil%GWdz, 4), ranges%GWdz, patchout%GWdz, out_settings)
    !    END IF

    IF (output%albsoil) THEN
      out_settings%dimswitch = "radiation"
      CALL  check_and_write(opid%albsoil, 'albsoil', REAL(soil%albsoil, 4), &
            ranges%albsoil, patchout%albsoil, out_settings)
    END IF

    IF (output%zse) THEN
      out_settings%dimswitch = "soil"
      CALL  check_and_write(opid%zse,      &
            'zse', SPREAD(REAL(soil%zse, 4), 1, mp),ranges%zse, &
            patchout%zse, out_settings)! no spatial dim at present
    END IF

    !~ Veg
    out_settings%dimswitch = "integer"
    IF (output%iveg) THEN
      CALL check_and_write(opid%iveg,    &
           'iveg', REAL(veg%iveg, 4), ranges%iveg, patchout%iveg, out_settings)
    END IF
    IF (output%meth) THEN
      CALL check_and_write(opid%meth,    &
           'meth', REAL(veg%meth, 4), ranges%meth, patchout%meth, out_settings)
    END IF

    out_settings%dimswitch = "real"
    IF (output%canst1) THEN
      CALL check_and_write(opid%canst1, 'canst1', REAL(veg%canst1, 4), &
           ranges%canst1, patchout%canst1, out_settings)
    END IF
    IF (output%dleaf) THEN
      CALL  check_and_write(opid%dleaf,  &
           'dleaf', REAL(veg%dleaf, 4), ranges%dleaf, patchout%dleaf, out_settings)
    END IF
    IF (output%ejmax) THEN
      CALL check_and_write(opid%ejmax,  &
           'ejmax', REAL(veg%ejmax, 4), ranges%ejmax, patchout%ejmax, out_settings)
    END IF
    IF (output%vcmax) THEN
    CALL check_and_write(opid%vcmax,  &
         'vcmax', REAL(veg%vcmax, 4), ranges%vcmax, patchout%vcmax, out_settings)
    END IF
    IF (output%frac4) THEN
      CALL check_and_write(opid%frac4,  &
           'frac4', REAL(veg%frac4, 4), ranges%frac4, patchout%frac4, out_settings)
    END IF
    IF (.NOT.cable_user%CALL_POP .and. output%hc) THEN
      CALL check_and_write(opid%hc,        &
            'hc', REAL(veg%hc, 4), ranges%hc, patchout%hc, out_settings)
    END IF
    IF (output%rp20) THEN
      CALL check_and_write(opid%rp20,    &
           'rp20', REAL(veg%rp20, 4),ranges%rp20, patchout%rp20, out_settings)
    END IF

    ! Ticket #56
    IF (output%g0) THEN
      CALL check_and_write(opid%g0,    &
           'g0', REAL(veg%g0, 4),ranges%g0, patchout%g0, out_settings)
    END IF
    IF (output%g1) THEN
      CALL check_and_write(opid%g1,    &
           'g1', REAL(veg%g1, 4),ranges%g1, patchout%g1, out_settings)
    END IF

    ! End Ticket #56
    IF (output%rpcoef) THEN
      CALL check_and_write(opid%rpcoef, 'rpcoef', REAL(veg%rpcoef, 4), &
           ranges%rpcoef, patchout%rpcoef, out_settings)
    END IF
    IF (output%shelrb) THEN
      CALL check_and_write(opid%shelrb, 'shelrb', REAL(veg%shelrb, 4), &
           ranges%shelrb, patchout%shelrb, out_settings)
    END IF
    IF (output%xfang) THEN
      CALL check_and_write(opid%xfang,  &
           'xfang', REAL(veg%xfang, 4), ranges%xfang, patchout%xfang, out_settings)
    END IF
    IF (output%wai) THEN
      CALL  check_and_write(opid%wai,      &
           'wai', REAL(veg%wai, 4), ranges%wai, patchout%wai, out_settings)
    END IF
    IF (output%vegcf) THEN
      CALL check_and_write(opid%vegcf,  &
           'vegcf', REAL(veg%vegcf, 4), ranges%vegcf, patchout%vegcf, out_settings)
    END IF
    IF (output%extkn) THEN
      CALL check_and_write(opid%extkn,  &
           'extkn', REAL(veg%extkn, 4), ranges%extkn, patchout%extkn, out_settings)
    END IF
    IF (output%tminvj) THEN
      CALL check_and_write(opid%tminvj, 'tminvj', REAL(veg%tminvj, 4), &
           ranges%tminvj, patchout%tminvj, out_settings)
    END IF
    IF (output%tmaxvj) THEN
      CALL check_and_write(opid%tmaxvj, 'tmaxvj', REAL(veg%tmaxvj, 4), &
           ranges%tmaxvj, patchout%tmaxvj, out_settings)
    END IF
    IF (output%vbeta) THEN
    CALL check_and_write(opid%vbeta,  &
        'vbeta', REAL(veg%vbeta, 4), ranges%vbeta, patchout%vbeta, out_settings)
    END IF
    IF (output%xalbnir) THEN
    CALL check_and_write(opid%xalbnir, 'xalbnir', REAL(veg%xalbnir, 4), &
         ranges%xalbnir, patchout%xalbnir, out_settings)
    END IF
    IF (output%froot) THEN
      out_settings%dimswitch = "soil"
      CALL check_and_write(opid%froot, &
           'froot', REAL(veg%froot, 4), ranges%froot, patchout%froot, out_settings)
    END IF

    !~ Rough
    out_settings%dimswitch = "real"
    IF (output%za) THEN
      CALL check_and_write(opid%za_uv,                                    &
           'za_uv', REAL(rough%za_uv, 4), ranges%za, patchout%za, out_settings)
    END IF
    IF (output%za) THEN
      CALL check_and_write(opid%za_tq,                                    &
           'za_tq', REAL(rough%za_tq, 4), ranges%za, patchout%za, out_settings)
    END IF

    !~ bgc
    IF (output%ratecp) THEN
      out_settings%dimswitch = "plantcarbon"
      CALL check_and_write(opid%ratecp, 'ratecp',SPREAD(REAL(bgc%ratecp,4),1,mp), ranges%ratecp, &
           patchout%ratecp, out_settings)! no spatial dim at present
    END IF
    IF (output%ratecs) THEN
      out_settings%dimswitch = "soilcarbon"
      CALL check_and_write(opid%ratecs, 'ratecs', SPREAD(REAL(bgc%ratecs, 4), 1, mp), ranges%ratecs,  &
           patchout%ratecs, out_settings)! no spatial dim at present
    END IF

    !~ gwmodel
    out_settings%dimswitch = "real"
    IF (output%params .AND. cable_user%gw_model) THEN
      CALL check_and_write(opid%SatFracmax,    &
           'SatFracmax', SPREAD(REAL(gw_params%MaxSatFraction,4),1,mp), &
           ranges%gw_default, patchout%SatFracmax, out_settings)

      CALL check_and_write(opid%Qhmax,    &
           'Qhmax', SPREAD(REAL(gw_params%MaxHorzDrainRate, 4),1,mp), &
           ranges%gw_default, patchout%Qhmax, out_settings)

      CALL check_and_write(opid%QhmaxEfold,    &
           'QhmaxEfold', SPREAD(REAL(gw_params%EfoldHorzDrainRate, 4),1,mp), &
           ranges%gw_default, patchout%QhmaxEfold, out_settings)

      CALL check_and_write(opid%HKefold,    &
           'HKefold', SPREAD(REAL(gw_params%hkrz, 4),1,mp), &
           ranges%gw_default, patchout%HKefold, out_settings)

      CALL check_and_write(opid%HKdepth,    &
           'HKdepth', SPREAD(REAL(gw_params%zdepth, 4),1,mp), &
            ranges%gw_default, patchout%HKdepth, out_settings)
    END IF

  END SUBROUTINE open_output_file

  !=============================================================================
  SUBROUTINE write_output(dels, ktau, met, canopy, casaflux, casapool, casamet, ssnow, &
                          rad, bal, air, soil, veg, SBOLTZ, EMLEAF, EMSOIL)
    ! Writes model output variables and, if requested, calls
    ! energy and mass balance routines. This subroutine is called
    ! each timestep, but may only write to the output file periodically,
    ! depending on whether the user has specified that output should be
    ! aggregated, e.g. to monthly or 6-hourly averages.
    REAL, INTENT(IN)              :: dels ! time step size
    INTEGER, INTENT(IN)           :: ktau ! timestep number in loop which include spinup
    REAL, INTENT(IN) :: SBOLTZ, EMLEAF, EMSOIL
    TYPE(met_type), TARGET, INTENT(IN)         :: met  ! met data
    TYPE(canopy_type), INTENT(IN)      :: canopy ! canopy variable data
    TYPE(soil_snow_type), INTENT(IN)   :: ssnow ! soil data
    TYPE(soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE(radiation_type), INTENT(IN)  :: rad   ! radiation data
    TYPE(air_type), INTENT(IN)        :: air
    TYPE(veg_parameter_type), INTENT(IN) :: veg ! vegetation parameters
    TYPE(casa_flux), INTENT(IN) :: casaflux ! casa fluxes
    TYPE(casa_pool), INTENT(IN) :: casapool ! casa fluxes
    TYPE(balances_type), INTENT(INOUT) :: bal
    TYPE(casa_met), INTENT(IN) :: casamet

    REAL(r_2) :: timetemp(1) ! temporary variable for storing time
    ! value
    INTEGER, SAVE :: out_month ! counter for output month
    INTEGER :: realyear(mp) ! fix problem for yr b4 leap yr
    INTEGER :: backtrack  ! modify timetemp for averaged output

    INTEGER :: dday ! number of past-years days for monthly output LN
    INTEGER :: iy   ! Counter
    !MC - use met%year(1) instead of CABLE_USER%YearStart for non-GSWP forcing and leap years
    INTEGER, SAVE :: YearStart

    INTEGER :: ok

    TYPE(output_var_settings_type) :: out_settings

    ! Temporary accumulation variable to be passed (we expect implicit-type conversion on assignment)
    ! Assumption: All variables have size mp, including CASA
    REAL(4) :: temp_acc(mp)

    out_settings = output_var_settings_type(met=met, writenow=.FALSE., dimswitch='default')

    ! IF asked to check mass/water balance:
    IF(check%mass_bal) CALL mass_balance(dels, ktau, ssnow, soil, canopy,            &
         met,air,bal)

    ! IF asked to check energy balance:
    IF(check%energy_bal) CALL energy_balance(dels,ktau,met,rad,                     &
         canopy,bal,ssnow,                 &
         SBOLTZ, EMLEAF, EMSOIL )
     
    ! Initialise output time step counter and month counter:
    IF(ktau == 1) THEN
       out_timestep = 0
       out_month = 0
       !MC - use met%year(1) instead of CABLE_USER%YearStart for non-GSWP forcing and leap years
       IF ( TRIM(cable_user%MetType) .EQ. '' ) THEN
          YearStart = met%year(1)
       ELSE
          YearStart = CABLE_USER%YearStart
       ENDIF
    END IF
    ! Decide on output averaging regime:
    IF(output%averaging(1:3) == 'all') THEN ! write every time step to file
       ! Set flag to write data for current time step:
       out_settings%writenow = .TRUE.
       ! Set output time step to be current model time step:
       out_timestep = ktau
       backtrack = 0
    ELSE IF(output%averaging(1:4) == 'user' .OR. output%averaging(1:2)=='da')  &
         THEN
       ! user defined output interval or daily output
       IF(MOD(ktau, output%interval) == 0) THEN ! i.e.ktau divisible by
          ! interval
          ! write to output file this time step
          out_settings%writenow = .TRUE.
          ! increment output time step counter:
          out_timestep = out_timestep + 1
          backtrack = output%interval / 2
       ELSE
          out_settings%writenow = .FALSE.
       END IF
    ELSE IF(output%averaging(1:2) == 'mo') THEN ! write monthly averages to file
       !realyear = met%year
       realyear = REAL(CurYear)
       IF(ktau >= 365*24*3600/INT(dels)) THEN
          WHERE(met%doy == 1) realyear = realyear - 1   ! last timestep of year
       END IF

       ! LN Inserted for multiyear output
       dday = 0
       !MC - use met%year(1) instead of CABLE_USER%YearStart for non-GSWP forcing and leap years
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
          IF (is_leapyear(CurYear)) THEN
             ! vh_js !
             IF(ANY(INT(REAL(lastdayl+dday) * 24. * 3600. / dels) == ktau)) THEN
                out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
                ! write to output file this time step
                out_settings%writenow = .TRUE.
                ! increment output time step counter:
                out_timestep = out_timestep + 1
                ! set numbr of time steps in output period
                output%interval = daysml(out_month) * 24 * 3600 / INT(dels)
             ELSE
                out_settings%writenow = .FALSE.
             END IF
          ELSE ! not currently a leap year
             ! last time step of month
             ! vh_js !
             IF(ANY(INT(REAL(lastday+dday) * 24. * 3600. / dels) == ktau)) THEN
                ! increment output month counter
                out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
                ! write to output file this time step
                out_settings%writenow = .TRUE.
                ! increment output time step counter:
                out_timestep = out_timestep + 1
                ! set numbr of time steps in output period
                output%interval = daysm(out_month) * 24 * 3600 / INT(dels)
             ELSE
                out_settings%writenow = .FALSE.
             END IF
          END IF
       ELSE ! not using leap year timing in this run

          ! vh_js !
          IF(ANY(INT((REAL((lastday+dday))*24.*3600./REAL(INT(dels))))==ktau)) THEN ! last time step of month
             ! IF(ANY(((lastday+dday)*24*3600/INT(dels))==ktau)) THEN ! last time step of month
             ! increment output month counter
             out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
             ! write to output file this time step
             out_settings%writenow = .TRUE.
             ! increment output time step counter:
             out_timestep = out_timestep + 1
             ! set numbr of time steps in output period
             output%interval = daysm(out_month) * 24 * 3600 / INT(dels)
          ELSE
             out_settings%writenow = .FALSE.
          END IF
       END IF ! using leap year timing or not
       backtrack = output%interval / 2


    ELSE ! type of output aggregation
       CALL abort('Unknown output averaging request in namelist file.'//       &
            '(SUBROUTINE write_output)')
    END IF

    ! Note that size of averaging interval, output%interval, is set when opening
    ! output file unless output is monthly (in which case it's set above)

    ! If this time step is an output time step:
    IF(out_settings%writenow) THEN
       ! Write to temporary time variable:
       timetemp(1) = DBLE(REAL(ktau-backtrack)*dels)
       ! Write time variable for this output time step:
       ok = NF90_PUT_VAR(ncid_out, ovid%tvar, timetemp,                        &
            start = (/out_timestep/), count = (/1/))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
            'Error writing time variable to ' &
            //TRIM(filename%out)// '(SUBROUTINE write_output)')
    END IF

    ! Arguments to generate_out_write_acc: current time step; output file netcdf file ID;
    ! netcdf variable ID; variable name; variable data; variable ranges;
    ! non-land fill value; include patch info for this var; any specific
    ! formatting info; met variables for reporting in case of abort.

    !-----------------------WRITE MET DATA-------------------------------------
    out_settings%dimswitch = 'default'
    IF (output%SWdown) THEN
      ! SWdown:  downward short-wave radiation [W/m^2]
      CALL generate_out_write_acc(ovid%SWdown, 'SWdown', out%SWdown, REAL(met%fsd(:, 1) + met%fsd(:, 2)), ranges%SWdown, patchout%SWdown, out_settings)
    END IF
    IF (output%LWdown) THEN
      ! LWdown: downward long-wave radiation [W/m^2]
      CALL generate_out_write_acc(ovid%LWdown, 'LWdown', out%LWdown, REAL(met%fld, 4), ranges%LWdown, patchout%LWdown, out_settings)
    END IF
    IF (output%Rainf) THEN
      ! Rainf: rainfall [kg/m^2/s]
      CALL generate_out_write_acc(ovid%Rainf, 'Rainf', out%Rainf, REAL(met%precip/dels, 4), ranges%Rainf, patchout%Rainf, out_settings)
    END IF
    IF (output%Snowf) THEN
      ! Snowf: snowfall [kg/m^2/s]
      CALL generate_out_write_acc(ovid%Snowf, 'Snowf', out%Snowf, REAL(met%precip_sn/dels, 4), ranges%Snowf, patchout%Snowf, out_settings)
    END IF
    IF (output%PSurf) THEN
      ! PSurf: surface pressure [Pa]
      CALL generate_out_write_acc(ovid%PSurf, 'PSurf', out%PSurf, REAL(met%pmb, 4), ranges%PSurf, patchout%PSurf, out_settings)
    END IF

    out_settings%dimswitch = 'ALMA'
    IF (output%Tair) THEN
      ! Tair: surface air temperature [K]
      CALL generate_out_write_acc(ovid%Tair, 'Tair', out%Tair, REAL(met%tk, 4), ranges%Tair, patchout%Tair, out_settings)
    END IF
    IF (output%Qair) THEN
      ! Qair: specific humidity [kg/kg]
      CALL generate_out_write_acc(ovid%Qair, 'Qair', out%Qair, REAL(met%qv, 4), ranges%Qair, patchout%Qair, out_settings)
    END IF
    IF (output%Wind) THEN
      ! Wind: windspeed [m/s]
      CALL generate_out_write_acc(ovid%Wind, 'Wind', out%Wind, REAL(met%ua, 4), ranges%Wind, patchout%Wind, out_settings)
    END IF
    IF (output%CO2air) THEN
      ! CO2air: CO2 concentration [ppmv]
      CALL generate_out_write_acc(ovid%CO2air, 'CO2air', out%CO2air, REAL(met%ca*1000000.0, 4), ranges%CO2air, patchout%CO2air, out_settings)
    END IF

    !-----------------------WRITE FLUX DATA-------------------------------------
    out_settings%dimswitch = 'default'
    ! Qmom: momentum flux [kg/m/s2] INH
    IF (output%Qmom) THEN
    CALL generate_out_write_acc(ovid%Qmom, 'Qmom', out%Qmom, REAL(air%rho, 4)*(REAL(canopy%us, 4)**2.), ranges%Qmom, patchout%Qmom, out_settings)
    END IF
    IF (output%Qmom) THEN
      ! Qle: latent heat flux [W/m^2]
      CALL generate_out_write_acc(ovid%Qle, 'Qle', out%Qle, REAL(canopy%fe, 4), ranges%Qle, patchout%Qle, out_settings)
    END IF
    IF (output%Qh) THEN
      ! Qh: sensible heat flux [W/m^2]
      CALL generate_out_write_acc(ovid%Qh, 'Qh', out%Qh, REAL(canopy%fh, 4), ranges%Qh, patchout%Qh, out_settings)
    END IF
    IF (output%Qg) THEN
      ! Qg: ground heat flux [W/m^2]
      CALL generate_out_write_acc(ovid%Qg, 'Qg', out%Qg, REAL(canopy%ga, 4), ranges%Qg, patchout%Qg, out_settings)
    END IF
    IF (output%Qs) THEN
      ! Qs: surface runoff [kg/m^2/s]
      CALL generate_out_write_acc(ovid%Qs, 'Qs', out%Qs, REAL(ssnow%rnof1/dels, 4), ranges%Qs, patchout%Qs, out_settings)
    END IF
    IF (output%Qsb) THEN
      ! Qsb: subsurface runoff [kg/m^2/s]
      CALL generate_out_write_acc(ovid%Qsb, 'Qsb', out%Qsb, REAL(ssnow%rnof2/dels, 4), ranges%Qsb, patchout%Qsb, out_settings)
    END IF
    IF (output%Evap) THEN
      ! Evap: total evapotranspiration [kg/m^2/s]
      CALL generate_out_write_acc(ovid%Evap, 'Evap', out%Evap, REAL(canopy%fe/air%rlam, 4), ranges%Evap, patchout%Evap, out_settings)
    END IF
    IF (output%PotEvap) THEN
      ! PotEVap: potential evapotranspiration [kg/m^2/s]
      CALL generate_out_write_acc(ovid%PotEvap, 'PotEvap', out%PotEvap, REAL(canopy%epot/dels, 4), ranges%PotEvap, patchout%PotEvap, out_settings)
    END IF
    IF (output%ECanop) THEN
      ! ECanop: interception evaporation [kg/m^2/s]
      CALL generate_out_write_acc(ovid%ECanop, 'ECanop', out%ECanop, REAL(canopy%fevw/air%rlam, 4), ranges%ECanop, patchout%ECanop, out_settings)
    END IF
    IF (output%TVeg) THEN
      ! TVeg: vegetation transpiration [kg/m^2/s]
      CALL generate_out_write_acc(ovid%TVeg, 'TVeg', out%TVeg, REAL(canopy%fevc/air%rlam, 4), ranges%TVeg, patchout%TVeg, out_settings)
    END IF

    IF (output%Esoil) THEN
      ! ESoil: bare soil evaporation [kg/m^2/s]
      IF (cable_user%SOIL_STRUC == 'sli') THEN
        temp_acc = ssnow%evap/dels !vh!
      ELSE
        temp_acc = canopy%fes/air%rlam
      END IF
      CALL generate_out_write_acc(ovid%Esoil, 'Esoil', out%Esoil, temp_acc, ranges%Esoil, patchout%Esoil, out_settings)
    END IF

    IF (output%HVeg) THEN
      ! HVeg: sensible heat from vegetation [W/m^2]
      CALL generate_out_write_acc(ovid%HVeg, 'HVeg', out%HVeg, REAL(canopy%fhv, 4), ranges%HVeg, patchout%HVeg, out_settings)
    END IF
    IF (output%HSoil) THEN
      ! HSoil: sensible heat from soil [W/m^2]
    CALL generate_out_write_acc(ovid%HSoil, 'HSoil', out%HSoil, REAL(canopy%fhs, 4), ranges%HSoil, patchout%HSoil, out_settings)
    END IF
    IF (output%RNetSoil) THEN
      CALL generate_out_write_acc(ovid%RNetSoil, 'RNetSoil', out%RNetSoil, REAL(canopy%fns, 4), ranges%HSoil, patchout%HSoil, out_settings)
    END IF
    IF (output%NEE) THEN
      ! NEE: net ecosystem exchange [umol/m^2/s]
      CALL generate_out_write_acc(ovid%NEE, 'NEE', out%NEE, REAL(canopy%fnee/1.201E-5, 4), ranges%NEE, patchout%NEE, out_settings)
    END IF

    !-----------------------WRITE SOIL STATE DATA-------------------------------

    out_settings%dimswitch = 'soil'
    IF (output%SoilMoist) THEN
      ! SoilMoist: av.layer soil moisture [kg/m^2]
      CALL generate_out_write_acc(ovid%SoilMoist, 'SoilMoist', out%SoilMoist, REAL(ssnow%wb, 4), ranges%SoilMoist, patchout%SoilMoistIce, out_settings)
      CALL generate_out_write_acc(ovid%SoilMoistIce, 'SoilMoistIce', out%SoilMoistIce, REAL(ssnow%wbice, 4), ranges%SoilMoist, patchout%SoilMoistIce, out_settings)
    END IF
    IF (output%SoilTemp) THEN
      ! SoilTemp: av.layer soil temperature [K]
      CALL generate_out_write_acc(ovid%SoilTemp, 'SoilTemp', out%SoilTemp, REAL(ssnow%tgg, 4), ranges%SoilTemp, patchout%SoilTemp, out_settings)
    END IF

    out_settings%dimswitch = 'default'
    IF (output%BaresoilT) THEN
      ! BaresoilT: surface bare soil temp [K]
      CALL generate_out_write_acc(ovid%BaresoilT, 'BaresoilT', out%BaresoilT, REAL(ssnow%tgg(:, 1), 4), ranges%BaresoilT, patchout%BaresoilT, out_settings)
    END IF
    !MD Write the hydrology output data from the groundwater module calculations
    IF (cable_user%GW_MODEL) THEN
      IF (output%WatTable) THEN
        !water table depth
        CALL generate_out_write_acc(ovid%WatTable, 'WatTable', out%WatTable, REAL(ssnow%wtd/1000.0, 4), ranges%WatTable, patchout%WatTable, out_settings)
      END IF
      IF (output%GWMoist) THEN
        !aquifer water content
        CALL generate_out_write_acc(ovid%GWMoist, 'GWMoist', out%GWMoist, REAL(ssnow%GWwb, 4), ranges%GWwb, patchout%GWMoist, out_settings)
      END IF
      IF (output%SatFrac) THEN
        !write(*,*) 'Qinfl'    !MDeck
        CALL generate_out_write_acc(ovid%SatFrac, 'SatFrac', out%SatFrac, REAL(ssnow%satfrac, 4), ranges%SatFrac, patchout%SatFrac, out_settings)
      END IF
    END IF

    IF (output%Qrecharge) THEN
      ! recharge rate
      CALL generate_out_write_acc(ovid%Qrecharge, 'Qrecharge', out%Qrecharge, REAL(ssnow%Qrecharge, 4), ranges%Qrecharge, patchout%Qrecharge, out_settings)
    END IF

    !----------------------WRITE SNOW STATE DATA--------------------------------
    IF (output%SWE) THEN
      ! SWE: snow water equivalent [kg/m^2]
      CALL generate_out_write_acc(ovid%SWE, 'SWE', out%SWE, REAL(ssnow%snowd, 4), ranges%SWE, patchout%SWE, out_settings)
      CALL generate_out_write_acc(ovid%SnowMelt, 'SnowMelt', out%SnowMelt, REAL(ssnow%smelt/dels, 4), ranges%SnowMelt, patchout%SnowMelt, out_settings)
    END IF

    IF (output%SnowT) THEN
      ! SnowT: snow surface temp [K]
      CALL generate_out_write_acc(ovid%SnowT, 'SnowT', out%SnowT, REAL(ssnow%tggsn(:, 1), 4), ranges%SnowT, patchout%SnowT, out_settings)
    END IF
    IF (output%SnowDepth) THEN
      ! SnowDepth: actual depth of snow in [m]
      CALL generate_out_write_acc(ovid%SnowDepth, 'SnowDepth', out%SnowDepth, REAL(SUM(ssnow%sdepth, 2), 4), ranges%SnowDepth, patchout%SnowDepth, out_settings)
    END IF

    !-------------------------WRITE RADIATION DATA------------------------------
    IF (output%Swnet) THEN
      ! SWnet: net shortwave [W/m^2]
      temp_acc = SUM(rad%qcan(:, :, 1), 2) + SUM(rad%qcan(:, :, 2), 2) + rad%qssabs
      CALL generate_out_write_acc(ovid%Swnet, 'Swnet', out%Swnet, temp_acc, ranges%Swnet, patchout%Swnet, out_settings)
    END IF
    IF (output%Lwnet) THEN
      ! LWnet: net longwave [W/m^2]
      temp_acc = met%fld - sboltz*emleaf*canopy%tv**4*(1 - rad%transd) - &
                 rad%flws*rad%transd
      CALL generate_out_write_acc(ovid%Lwnet, 'Lwnet', out%Lwnet, temp_acc, ranges%Lwnet, patchout%Lwnet, out_settings)
    END IF
    IF (output%Rnet) THEN
      ! Rnet: net absorbed radiation [W/m^2]
      temp_acc = met%fld - sboltz*emleaf*canopy%tv**4* &
                 (1 - rad%transd) - rad%flws*rad%transd + &
                 SUM(rad%qcan(:, :, 1), 2) + &
                 SUM(rad%qcan(:, :, 2), 2) + rad%qssabs
      CALL generate_out_write_acc(ovid%Rnet, 'Rnet', out%Rnet, temp_acc, ranges%Rnet, patchout%Rnet, out_settings)
    END IF

    IF (output%Albedo) THEN
      ! Albedo:
      CALL generate_out_write_acc(ovid%Albedo, 'Albedo', out%Albedo, REAL((rad%albedo(:, 1) + rad%albedo(:, 2)) &
                                                                                       *0.5, 4), ranges%Albedo, patchout%Albedo, out_settings)
      IF (calcsoilalbedo) THEN
        CALL generate_out_write_acc(ovid%visAlbedo, 'visAlbedo', out%visAlbedo, REAL(rad%albedo(:, 1), 4), ranges%visAlbedo, patchout%visAlbedo, out_settings)
        CALL generate_out_write_acc(ovid%nirAlbedo, 'nirAlbedo', out%nirAlbedo, REAL(rad%albedo(:, 2), 4), ranges%nirAlbedo, patchout%nirAlbedo, out_settings)
      END IF
    END IF


    ! RadT: Radiative surface temperature [K]
    IF (output%RadT) THEN
      temp_acc = (((1.0 - rad%transd)*emleaf*sboltz* &
                   canopy%tv**4 + rad%transd*emsoil*sboltz*(ssnow%tss)**4)/sboltz)**0.25
      CALL generate_out_write_acc(ovid%RadT, 'RadT', out%RadT, temp_acc, ranges%RadT, patchout%RadT, out_settings)
    END IF

    !------------------------WRITE VEGETATION DATA------------------------------

    out_settings%dimswitch = 'ALMA'
    IF (output%Tscrn) THEN
      ! Tscrn: screen level air temperature [oC]
      CALL generate_out_write_acc(ovid%Tscrn, 'Tscrn', out%Tscrn, REAL(canopy%tscrn, 4), ranges%Tscrn, patchout%Tscrn, out_settings)
    END IF

    !INH - extremes in screen level air temperature [oC]
    IF (output%Tex) THEN
      !if 'daily' then only daily values - using variables Txx and Tnn
      IF (output%averaging(1:2) == 'da') THEN
        DO iy = 1, mp
          out%Txx(iy) = MAX(out%Txx(iy), REAL(canopy%tscrn(iy), 4))
          out%Tnn(iy) = MIN(out%Tnn(iy), REAL(canopy%tscrn(iy), 4))
        END DO
        IF (out_settings%writenow) THEN
          CALL check_and_write(ovid%Txx, 'Txx', &
                               out%Txx, out%Txx, ranges%Tscrn, patchout%Tex, out_settings)
          CALL check_and_write(ovid%Tnn, 'Tnn', &
                               out%Tnn, out%Tnn, ranges%Tscrn, patchout%Tex, out_settings)
          !Reset temporary output variables:
          out%Txx = -1.0E6
          out%Tnn = 1.0E6
        END IF
      END IF

      IF (output%averaging(1:2) == 'mo') THEN
        !if monthly then both full extremes and averaged extremes
        DO iy = 1, mp
          out%Txx(iy) = MAX(out%Txx(iy), REAL(canopy%tscrn(iy), 4))
          out%Tnn(iy) = MIN(out%Tnn(iy), REAL(canopy%tscrn(iy), 4))
          out%Tdaymx(iy) = MAX(out%Tdaymx(iy), REAL(canopy%tscrn(iy), 4))
          out%Tdaymn(iy) = MIN(out%Tdaymn(iy), REAL(canopy%tscrn(iy), 4))
        END DO
        !take copy of day's max/min for averaged output - reset Tdaymx/mn
        IF (MOD(ktau, 24*3600/INT(dels)) == 0) THEN
          out%Tmx = out%Tmx + out%Tdaymx
          out%Tmn = out%Tmn + out%Tdaymn
          out%Tdaymx = -1.0E6
          out%Tdaymn = 1.0E6
        END IF
        IF (out_settings%writenow) THEN
          !divide by number of records in average (dels*%interval/24/3600)
          out%Tmx = REAL(86400, 4)*out%Tmx/REAL(output%interval*INT(dels), 4)
          out%Tmn = REAL(86400, 4)*out%Tmn/REAL(output%interval*INT(dels), 4)
          !write to file
          CALL check_and_write(ovid%Txx, 'Txx', &
                               out%Txx, out%Txx, ranges%Tscrn, patchout%Tex, out_settings)
          CALL check_and_write(ovid%Tnn, 'Tnn', &
                               out%Tnn, out%Tnn, ranges%Tscrn, patchout%Tex, out_settings)
          CALL check_and_write(ovid%Tmx, 'Tmx', &
                               out%Tmx, out%Tmx, ranges%Tscrn, patchout%Tex, out_settings)
          CALL check_and_write(ovid%Tmn, 'Tmn', &
                               out%Tmn, out%Tmn, ranges%Tscrn, patchout%Tex, out_settings)
          !Reset temporary output variables:
          out%Txx = -1.0E6
          out%Tnn = 1.0E6
          out%Tmx = 0.0
          out%Tmn = 0.0
        END IF
      END IF
    END IF

    IF (output%Qscrn) THEN
      ! Qscrn: screen level specific humdity [kg/kg]
      CALL generate_out_write_acc(ovid%qscrn, 'Qscrn', out%qscrn, REAL(canopy%qscrn, 4), ranges%Qscrn, patchout%Qscrn, out_settings)
    END IF

    out_settings%dimswitch = 'default'
    IF (output%VegT) THEN
      ! VegT: vegetation temperature [K]
      CALL generate_out_write_acc(ovid%VegT, 'VegT', out%VegT, REAL(canopy%tv, 4), ranges%VegT, patchout%VegT, out_settings)
    END IF
    IF (output%CanT) THEN
      ! CanT: within-canopy temperature [K]
      CALL generate_out_write_acc(ovid%CanT, 'CanT', out%CanT, REAL(met%tvair, 4), ranges%CanT, patchout%CanT, out_settings)
    END IF
    IF (output%Fwsoil) THEN
      ! Fwsoil
      CALL generate_out_write_acc(ovid%Fwsoil, 'Fwsoil', out%Fwsoil, REAL(canopy%fwsoil, 4), ranges%Fwsoil, patchout%Fwsoil, out_settings)
    END IF
    IF (output%CanopInt) THEN
      ! CanopInt: total canopy water storage [kg/m^2]
      CALL generate_out_write_acc(ovid%CanopInt, 'CanopInt', out%CanopInt, REAL(canopy%cansto, 4), ranges%CanopInt, patchout%CanopInt, out_settings)
    END IF
    IF (output%LAI) THEN
      ! LAI:
      CALL generate_out_write_acc(ovid%LAI, 'LAI', out%LAI, REAL(veg%vlai, 4), ranges%LAI, patchout%LAI, out_settings)
    END IF
    !------------------------WRITE BALANCES DATA--------------------------------
    IF (output%Ebal) THEN
      ! Ebal: cumulative energy balance [W/m^2]
      CALL generate_out_write_acc(ovid%Ebal, 'Ebal', out%Ebal, REAL(bal%ebal_tot, 4), ranges%Ebal, patchout%Ebal, out_settings)
    END IF
    IF (output%Wbal) THEN
      ! Wbal: cumulative water balance  [kg/m^2/s]
      CALL generate_out_write_acc(ovid%Wbal, 'Wbal', out%Wbal, REAL(bal%wbal_tot, 4), ranges%Wbal, patchout%Wbal, out_settings)
    END IF
    !------------------------WRITE CARBON DATA----------------------------------
    ! GPP: gross primary production C by veg [umol/m^2/s]
    !      added frday in the calculation of GPP (BP may08)

    ! temp_acc = REAL((-1.0 * canopy%fpn) / 1.201E-5, 4)
    IF (output%GPP) THEN
      CALL generate_out_write_acc(ovid%GPP, 'GPP', out%GPP, REAL((-1.0*canopy%fpn + canopy%frday) &
        /1.201E-5, 4), ranges%GPP, patchout%GPP, out_settings)
    END IF

    ! NPP: net primary production of C by veg [umol/m^2/s]
    IF (output%NPP) THEN
      ! Add current timestep's value to total of temporary output variable:
      !out%NPP = out%NPP + REAL((-1.0 * canopy%fpn - canopy%frp &
      !     - casaflux%clabloss/86400.0) / 1.201E-5, 4)
      ! vh ! expression below can be slightly different form that above in cases where
      ! leaf maintenance respiration is reduced in CASA
      ! (relative to its original value calculated in cable_canopy)
      ! in order to avoid negative carbon stores.
      IF (output%casa) THEN
        temp_acc = casaflux%cnpp/86400.0/1.201E-5
      ELSE
        temp_acc = (-1.0*canopy%fpn - canopy%frp)/1.201E-5 ! &
      END IF
      CALL generate_out_write_acc(ovid%NPP, 'NPP', out%NPP, temp_acc, ranges%NPP, patchout%NPP, out_settings)
    END IF

    ! AutoResp: autotrophic respiration [umol/m^2/s]
    IF (output%AutoResp) THEN
      ! Add current timestep's value to total of temporary output variable:
      !out%AutoResp = out%AutoResp + REAL((canopy%frp + canopy%frday + casaflux%clabloss/86400.0)          &
      !                                    / 1.201E-5, 4)
      ! vh ! expression below can be slightly different form that above in cases where
      ! leaf maintenance respiration is reduced in CASA
      ! (relative to its original value calculated in cable_canopy)
      ! in order to avoid negative carbon stores.

      IF (output%casa) THEN
        temp_acc = canopy%frday/1.201E-5 + &
                   (casaflux%crmplant(:, 2)/86400.0 + casaflux%crmplant(:, 3)/86400.0 + &
                    casaflux%crgplant/86400.0 + casaflux%clabloss/86400.)/1.201E-5
      ELSE
        temp_acc = (canopy%frp + canopy%frday)/1.201E-5
      END IF
      CALL generate_out_write_acc(ovid%AutoResp, 'AutoResp', out%AutoResp, temp_acc, ranges%AutoResp, patchout%AutoResp, out_settings)
      IF (output%casa) THEN
        ! rootresp alt: REAL(0.3*casaflux%crmplant(:,2)/86400.0/ 1.201E-5, 4)
        CALL generate_out_write_acc(ovid%RootResp, 'RootResp', out%RootResp, &
                                    REAL(casaflux%crmplant(:, 3)/86400.0/1.201E-5, 4), ranges%AutoResp, patchout%AutoResp, out_settings)
        CALL generate_out_write_acc(ovid%StemResp, 'StemResp', out%StemResp, &
                                    REAL(casaflux%crmplant(:, 2)/86400.0/1.201E-5, 4), ranges%AutoResp, patchout%AutoResp, out_settings)
      END IF
    END IF

    IF (output%LeafResp) THEN
      ! LeafResp: Leaf respiration [umol/m^2/s]
      CALL generate_out_write_acc(ovid%LeafResp, 'LeafResp', out%LeafResp, REAL(canopy%frday/1.201E-5, 4), ranges%LeafResp, patchout%LeafResp, out_settings)
    END IF
    IF (output%HeteroResp) THEN
      ! HeteroResp: heterotrophic respiration [umol/m^2/s]
      CALL generate_out_write_acc(ovid%HeteroResp, 'HeteroResp', out%HeteroResp, REAL(canopy%frs/1.201E-5, 4), ranges%HeteroResp, patchout%HeteroResp, out_settings)
    END IF

    ! output patch area
    IF (output%casa) THEN
      out%Area = casamet%areacell/1e6 ! km2

      CALL check_and_write(ovid%Area, 'Area', out%Area, out%Area, ranges%Area, patchout%Area, out_settings)
      IF (cable_user%POPLUC) THEN
        CALL check_and_write(ovid%patchfrac, 'patchfrac', REAL(patch(:)%frac, 4), REAL(patch(:)%frac, 4), ranges%Area, patchout%Area, out_settings)
      END IF
      IF (cable_user%CALL_POP) THEN
        CALL check_and_write(ovid%hc, 'hc', REAL(veg%hc, 4), REAL(veg%hc, 4), ranges%hc, patchout%hc, out_settings)
      END IF

      ! NBP and turnover fluxes [umol/m^2/s]
      IF (output%NBP) THEN
        IF (cable_user%POPLUC) THEN
          temp_acc = -(casaflux%Crsoil - casaflux%cnpp &
                       - casapool%dClabiledt)/86400.0 &
                     /1.201E-5 !-  &
          !REAL((casaflux%FluxCtohwp + casaflux%FluxCtoclear  )/86400.0 &
          !/ 1.201E-5, 4)
        ELSE
          temp_acc = -(casaflux%Crsoil - casaflux%cnpp &
                       - casapool%dClabiledt)/86400.0 &
                     /1.201E-5
        END IF
        CALL generate_out_write_acc(ovid%NBP, 'NBP', out%NBP, temp_acc, ranges%NEE, patchout%NBP, out_settings)
      END IF

      !------------------------WRITE REMAINING CASA DATA----------------------------------
      CALL generate_out_write_acc(ovid%dCdt, 'dCdt', out%dCdt, &
                                  REAL((casapool%ctot - casapool%ctot_0)/86400.0/1.201E-5, 4), ranges%NEE, patchout%dCdt, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnover, 'PlantTurnover', out%PlantTurnover, &
                                  REAL((SUM(casaflux%Cplant_turnover, 2))/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnover, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnoverLeaf, 'PlantTurnoverLeaf', out%PlantTurnoverLeaf, &
                                  REAL((casaflux%Cplant_turnover(:, 1))/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnoverLeaf, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnoverFineRoot, 'PlantTurnoverFineRoot', out%PlantTurnoverFineRoot, &
                                  REAL((casaflux%Cplant_turnover(:, 3))/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnoverFineRoot, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnoverWood, 'PlantTurnoverWood', out%PlantTurnoverWood, &
                                  REAL((casaflux%Cplant_turnover(:, 2))/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnoverWood, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnoverWoodDist, 'PlantTurnoverWoodDist', out%PlantTurnoverWoodDist, &
                                  REAL(casaflux%Cplant_turnover_disturbance/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnoverWoodDist, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnoverWoodCrowding, 'PlantTurnoverWoodCrowding', out%PlantTurnoverWoodCrowding, &
                                  REAL(casaflux%Cplant_turnover_crowding/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnoverWoodCrowding, out_settings)
      CALL generate_out_write_acc(ovid%PlantTurnoverWoodResourceLim, 'PlantTurnoverWoodResourceLim', out%PlantTurnoverWoodResourceLim, &
                                  REAL((casaflux%Cplant_turnover_resource_limitation)/86400.0/1.201E-5, 4), ranges%NEE, patchout%PlantTurnoverWoodResourceLim, out_settings)
      IF (cable_user%POPLUC) THEN
        CALL generate_out_write_acc(ovid%LandUseFlux, 'LandUseFlux', out%LandUseFlux, &
                                    REAL((casaflux%FluxCtohwp + casaflux%FluxCtoclear)/86400.0/1.201E-5, 4), ranges%NEE, patchout%LandUseFlux, out_settings)
      END IF

      ! plant carbon [kg C m-2]
      CALL generate_out_write_acc(ovid%TotSoilCarb, 'TotSoilCarb', out%TotSoilCarb, REAL((SUM(casapool%csoil, 2) + SUM(casapool%clitter, 2)) &
                                                                                                    /1000.0, 4), ranges%TotSoilCarb, patchout%TotSoilCarb, out_settings)
      CALL generate_out_write_acc(ovid%TotLittCarb, 'TotLittCarb', out%TotLittCarb, REAL(SUM(casapool%clitter, 2)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%TotLittCarb, out_settings)

      ! csoil
      CALL generate_out_write_acc(ovid%SoilCarbFast, 'SoilCarbFast', out%SoilCarbFast, REAL(casapool%csoil(:, 1)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%SoilCarbFast, out_settings)
      CALL generate_out_write_acc(ovid%SoilCarbSlow, 'SoilCarbSlow', out%SoilCarbSlow, REAL(casapool%csoil(:, 2)/1000.0, 4), &
                                  ranges%TotSoilCarb, patchout%SoilCarbSlow, out_settings)
      CALL generate_out_write_acc(ovid%SoilCarbPassive, 'SoilCarbPassive', out%SoilCarbPassive, REAL(casapool%csoil(:, 3)/1000.0, 4), &
                                  ranges%TotSoilCarb, patchout%SoilCarbPassive, out_settings)

      ! clitter
      CALL generate_out_write_acc(ovid%LittCarbMetabolic, 'LittCarbMetabolic', out%LittCarbMetabolic, REAL(casapool%clitter(:, 1)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%LittCarbMetabolic, out_settings)
      CALL generate_out_write_acc(ovid%LittCarbStructural, 'LittCarbStructural', out%LittCarbStructural, REAL(casapool%clitter(:, 2)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%LittCarbStructural, out_settings)
      CALL generate_out_write_acc(ovid%LittCarbCWD, 'LittCarbCWD', out%LittCarbCWD, REAL(casapool%clitter(:, 3)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%LittCarbCWD, out_settings)

      ! cplant
      CALL generate_out_write_acc(ovid%PlantCarbLeaf, 'PlantCarbLeaf', out%PlantCarbLeaf, REAL(casapool%cplant(:, 1)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%PlantCarbLeaf, out_settings)
      CALL generate_out_write_acc(ovid%PlantCarbWood, 'PlantCarbWood', out%PlantCarbWood, REAL(casapool%cplant(:, 2)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%PlantCarbWood, out_settings)
      CALL generate_out_write_acc(ovid%PlantCarbFineRoot, 'PlantCarbFineRoot', out%PlantCarbFineRoot, REAL(casapool%cplant(:, 3)/1000.0, 4), &
                                  ranges%TotLittCarb, patchout%PlantCarbFineRoot, out_settings)

      CALL generate_out_write_acc(ovid%TotLivBiomass, 'TotLivBiomass', out%TotLivBiomass, REAL((SUM(casapool%cplant, 2)) &
                                                                                                /1000.0, 4), ranges%TotLivBiomass, patchout%TotLivBiomass, out_settings)
    END IF

    IF (cable_user%sync_nc_file) &
      ok = NF90_SYNC(ncid_out)

  END SUBROUTINE write_output

  SUBROUTINE check_and_write_d1(varID, vname, out_var, acc_val, vrange, writepatch, out_settings)
    INTEGER, INTENT(IN) :: varID ! variable's netcdf ID
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    REAL(4), INTENT(IN) :: out_var(:)
    REAL(4), INTENT(IN) :: acc_val(:)
    REAL, INTENT(IN) :: vrange(2) ! max and min for variable
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    TYPE(output_var_settings_type), INTENT(IN) :: out_settings ! output specific settings

    IF ((check%ranges .EQ. ON_TIMESTEP) .OR. (out_settings%writenow .AND. (check%ranges .EQ. ON_WRITE))) THEN
      CALL check_range(vname, acc_val, vrange, out_timestep, out_settings%met)
    END IF

    IF (out_settings%writenow) THEN
      ! Write value to file:
      CALL write_ovar(out_timestep, ncid_out, varID, vname, &
                      out_var, writepatch, out_settings%dimswitch, out_settings%met)
    END IF
  END SUBROUTINE check_and_write_d1

  SUBROUTINE check_and_write_d2(varID, vname, out_var, acc_val, vrange, writepatch, out_settings)
    INTEGER, INTENT(IN) :: varID ! variable's netcdf ID
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    REAL(4), INTENT(IN) :: out_var(:, :)
    REAL(4), INTENT(IN) :: acc_val(:, :)
    REAL, INTENT(IN) :: vrange(2) ! max and min for variable
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    TYPE(output_var_settings_type), INTENT(IN) :: out_settings  ! output specific settings

    IF ((check%ranges .EQ. ON_TIMESTEP) .OR. (out_settings%writenow .AND. (check%ranges .EQ. ON_WRITE))) THEN
      CALL check_range(vname, acc_val, vrange, out_timestep, out_settings%met)
    END IF

    IF (out_settings%writenow) THEN
      ! Write value to file:
      CALL write_ovar(out_timestep, ncid_out, varID, vname, &
                      out_var, writepatch, out_settings%dimswitch, out_settings%met)
    END IF
  END SUBROUTINE check_and_write_d2

  SUBROUTINE check_and_write_d1_p(parID, pname, out_par, prange, writepatch, out_settings)
    INTEGER, INTENT(IN) :: parID ! parameter netcdf ID
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of parameter
    REAL(4), INTENT(IN) :: out_par(:)
    REAL, INTENT(IN) :: prange(2) ! max and min for parameter
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this par?
    TYPE(output_par_settings_type), INTENT(IN) :: out_settings  ! output specific settings

    INTEGER :: ncid_file

    CALL check_range(pname, out_par, prange, out_timestep, out_settings%met)

    IF (out_settings%restart) THEN
      ncid_file = ncid_restart
    ELSE
      ncid_file = ncid_out
    END IF
    ! Write value to file:
    CALL write_ovar(ncid_file, parID, pname, out_par, &
                    writepatch, out_settings%dimswitch, out_settings%restart)

  END SUBROUTINE check_and_write_d1_p

  SUBROUTINE check_and_write_d2_p(parID, pname, out_par, prange, writepatch, out_settings)
    INTEGER, INTENT(IN) :: parID ! parameter netcdf ID
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of parameter
    REAL(4), INTENT(IN) :: out_par(:, :)
    REAL, INTENT(IN) :: prange(2) ! max and min for parameter
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this par?
    TYPE(output_par_settings_type), INTENT(IN) :: out_settings ! output specific settings

    INTEGER :: ncid_file

    CALL check_range(pname, out_par, prange, out_timestep, out_settings%met)

    IF (out_settings%restart) THEN
      ncid_file = ncid_restart
    ELSE
      ncid_file = ncid_out
    END IF
    ! Write value to file:
    CALL write_ovar(ncid_file, parID, pname, out_par, &
                      writepatch, out_settings%dimswitch, out_settings%restart)

  END SUBROUTINE check_and_write_d2_p

  SUBROUTINE generate_out_write_acc_d1(varID, vname, out_var, acc_val, vrange, writepatch, out_settings)
    INTEGER, INTENT(IN) :: varID ! variable's netcdf ID
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    REAL(4), INTENT(INOUT) :: out_var(:)
    REAL(4), INTENT(IN) :: acc_val(:)
    REAL, INTENT(IN) :: vrange(2) ! max and min for variable
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    TYPE(output_var_settings_type), INTENT(IN) :: out_settings ! met data

    ! Accumulate out_var until interval timesteps
    out_var = out_var + acc_val
    IF (out_settings%writenow) THEN
      out_var = out_var/REAL(output%interval, 4)
    END IF

    CALL check_and_write(varID, vname, out_var, acc_val, vrange, writepatch, out_settings)

    ! Reset the value if it has been written to file
    IF (out_settings%writenow) THEN
      out_var = 0.0
    END IF

  END SUBROUTINE generate_out_write_acc_d1

  SUBROUTINE generate_out_write_acc_d2(varID, vname, out_var, acc_val, vrange, writepatch, out_settings)
    INTEGER, INTENT(IN) :: varID ! variable's netcdf ID
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    REAL(4), INTENT(INOUT) :: out_var(:, :)
    REAL(4), INTENT(IN) :: acc_val(:, :)
    REAL, INTENT(IN) :: vrange(2) ! max and min for variable
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    TYPE(output_var_settings_type), INTENT(IN) :: out_settings ! met data

    ! Accumulate out_var until interval timesteps
    out_var = out_var + acc_val
    IF (out_settings%writenow) THEN
      out_var = out_var/REAL(output%interval, 4)
    END IF

    CALL check_and_write(varID, vname, out_var, acc_val, vrange, writepatch, out_settings)

    ! Reset the value if it has been written to file
    IF (out_settings%writenow) THEN
      out_var = 0.0
    END IF

  END SUBROUTINE generate_out_write_acc_d2

  !=============================================================================
  SUBROUTINE close_output_file(bal, air, bgc, canopy, met,                     &
       rad, rough, soil, ssnow, sum_flux, veg)
    ! Closes output file, reports cumulative mass and energy
    ! balances, and deallocates variables.
    TYPE (met_type), INTENT(INOUT)       :: met
    TYPE (air_type), INTENT(INOUT)       :: air
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (radiation_type),INTENT(INOUT)  :: rad
    TYPE (sum_flux_type), INTENT(INOUT)  :: sum_flux
    TYPE(balances_type),INTENT(INOUT)    :: bal

    INTEGER :: i ! do loop counter

    ! Close file
    ok = NF90_CLOSE(ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing output file '        &
         //TRIM(filename%out)// '(SUBROUTINE close_output_file)')

    ! Report balance info to log file if verbose writing is requested:
    IF(output%balances .AND. verbose) THEN
       WRITE(logn, *)
       DO i = 1, mland
          WRITE(logn, '(A51,I7,1X,A11,E12.4,A6)')                              &
               ' Cumulative energy balance for each patch in site #',         &
               i,'is (W/m^2):'
          WRITE(logn, *)                                                       &
               bal%ebal_tot(landpt(i)%cstart:landpt(i)%cstart +               &
               landpt(i)%nap - 1)
          WRITE(logn,'(A50,I7,1X,A8,E12.4,A3)')                                &
               ' Cumulative water balance for each patch in site #',          &
               i,'is (mm):'
          WRITE(logn, *)                                                       &
               bal%wbal_tot(landpt(i)%cstart:landpt(i)%cstart +               &
               landpt(i)%nap - 1)
          WRITE(logn, *)
       END DO
    END IF

    ! Successful run!
    WRITE(logn, *)
    WRITE(logn, *) 'Run finished and output file closed.'

  END SUBROUTINE close_output_file
  !=============================================================================
  SUBROUTINE create_restart(logn, dels, ktau, soil, veg, ssnow,                      &
       canopy, rough, rad, bgc, bal, met)
    ! Creates a restart file for CABLE using a land only grid with mland
    ! land points and max_vegpatches veg/soil patches (some of which may
    ! not be active). It uses CABLE's internal variable names.
    INTEGER, INTENT(IN) :: logn ! log file number
    REAL, INTENT(IN) :: dels ! time step size
    INTEGER, INTENT(IN)           :: ktau ! timestep number in loop which include spinup
    TYPE (met_type), TARGET, INTENT(IN)             :: met ! meteorological data
    TYPE (soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)  :: veg  ! vegetation parameters
    TYPE (soil_snow_type), INTENT(IN)      :: ssnow  ! soil and snow variables
    TYPE (bgc_pool_type), INTENT(IN)       :: bgc    ! carbon pool variables
    TYPE (canopy_type), INTENT(IN)         :: canopy ! vegetation variables
    TYPE (roughness_type), INTENT(IN)      :: rough  ! roughness varibles
    TYPE (radiation_type), INTENT(IN)  :: rad ! radiation variables
    TYPE (balances_type), INTENT(IN) :: bal ! energy and water balance variables
    !    INTEGER, INTENT(IN) :: mvtype
    !    INTEGER, INTENT(IN) :: mstype

    TYPE(parID_type) :: rpid ! parameter IDs for restart nc file


    TYPE(output_par_settings_type) :: out_settings

    LOGICAL, PARAMETER :: patchout_var = .TRUE.

    ! REAL, POINTER, :: surffrac(:, :) ! fraction of each surf type
    INTEGER :: dummy ! dummy argument in subroutine call
    INTEGER :: mlandID, mpID, radID, soilID, napID,                       &
         soilcarbID, plantcarbID, tID, snowID ! dimension IDs
    !    INTEGER :: mlandID, surftypeID, patchID, radID, soilID, &
    !         soilcarbID, plantcarbID, tID, snowID ! dimension IDs
    INTEGER :: tvarID, latID, lonID !,surffracID ! time,lat,lon variable ID
    INTEGER :: tggID, wbID, wbiceID, tssID, ssdnnID, ssdnID, osnowdID,    &
         smassID, sdepthID, snageID, snowdID, rtsoilID, isflagID,   &
         canstoID, albsoilsnID, gammzzID, tggsnID, sghfluxID,       &
         ghfluxID, runoffID, rnof1ID, rnof2ID, gaID, dgdtgID,       &
         fevID, fesID, fhsID, wbtot0ID, osnowd0ID, cplantID,        &
         csoilID, tradID, albedoID, gwID
    INTEGER :: h0ID, snowliqID, SID, TsurfaceID, scondsID, nsnowID, TsoilID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
    ! CHARACTER         :: FRST_OUT*100, CYEAR*4
    CHARACTER         :: FRST_OUT*200, CYEAR*4

    out_settings = output_par_settings_type(met=met, restart=.TRUE.)

    dummy = 0 ! initialise

    WRITE(logn, '(A24)') ' Writing restart file...'
    frst_out = TRIM(filename%restart_out)
    ! Look for explicit restart file (netCDF). If not, asssume input is path
    IF ( INDEX(TRIM(frst_out),'.nc',BACK=.TRUE.) .NE. LEN_TRIM(frst_out)-2 ) THEN
       WRITE( CYEAR,FMT="(I4)" ) CurYear + 1
       frst_out = TRIM(cable_user%RunIden)//'_'//CYEAR//'_cable_rst.nc'
    ENDIF

    ! Create output file:
    ok = NF90_CREATE(frst_out, NF90_CLOBBER, ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file '      &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ! Put the file in define mode:
    ok = NF90_REDEF(ncid_restart)
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_restart, 'mland', mland, mlandID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mland dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'mp', mp, mpID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mp dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining vertical soil dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'snow', 3, snowID) ! number of snow layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining vertical snow dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'rad', nrb, radID) ! number of rad. bands
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining radiation dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'soil_carbon_pools', ncs, soilcarbID)
    ! number of soil carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining soil carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'plant_carbon_pools', ncp, plantcarbID)
    ! number of plant carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining plant carbon pool dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'time', 1, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok, 'Error defining time dimension in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define "time" variable and its attributes:
    ok=NF90_DEF_VAR(ncid_restart,'time',NF90_DOUBLE,(/tID/),tvarID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'calendar', calendar)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining time variable attribute calendar in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define latitude and longitude variable:
    ok=NF90_DEF_VAR(ncid_restart, 'latitude', NF90_FLOAT, (/mlandID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart,latID,'units','degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok=NF90_DEF_VAR(ncid_restart, 'longitude', NF90_FLOAT, (/mlandID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining longitude variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining longitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define number of active patches variable:
    ok = NF90_DEF_VAR(ncid_restart, 'nap', NF90_FLOAT, (/mlandID/), napID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining nap variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, napID, 'long_name',                        &
         'Number of active patches')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok,'Error defining nap variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! Define patch fraction variable:
    ok=NF90_DEF_VAR(ncid_restart, 'patchfrac', NF90_FLOAT, (/mpID/),           &
         rpid%patchfrac)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining patchfrac variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%patchfrac, 'long_name',               &
         'Fraction of vegetated grid cell area occupied by a '//  &
         'vegetation/soil patch')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining patchfrac variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ! mvtype (Number of vegetation types):
    ok = NF90_DEF_VAR(ncid_restart, 'mvtype', NF90_INT, rpid%mvtype)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mvtype variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%mvtype, "long_name",                  &
         "Number of vegetation types")
    ! mstype (Number of soil types):
    ok = NF90_DEF_VAR(ncid_restart, 'mstype', NF90_INT, rpid%mstype)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining mstype variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%mstype, "long_name",                  &
         "Number of soil types")

    !======begin defining state variables=======================================
    ! Interface arguments: netcdf file ID, variableID, variable name, variable
    ! units, variable long name, YES to write patch info (as this is a restart
    ! file), OPTIONAL extra dimension ID (e.g. for soil dimensioned variables),
    ! dimension switch to indicate what extra dimension is real or integer for
    ! single dim variables, xdimID,ydimID, zdimID (all three not used here),
    ! land dim ID, patch dim ID, YES we're writing a restart file.
    !------------------define soil states---------------------------------------
    CALL define_ovar(ncid_restart, tggID, 'tgg', 'K',                          &
         'Average layer soil temperature',                         &
         .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wbID, 'wb', 'vol/vol',                      &
         'Average layer volumetric soil moisture',                 &
         .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wbiceID, 'wbice', 'vol/vol',                &
         'Average layer volumetric soil ice',                      &
         .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, tssID, 'tss', 'K',                          &
         'Combined soil/snow temperature',                         &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, albsoilsnID, 'albsoilsn', '-',              &
         'Combined soil/snow albedo',                              &
         .TRUE., radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rtsoilID, 'rtsoil', '??',                   &
         'Turbulent resistance for soil', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gammzzID, 'gammzz', 'J/kg/C',               &
         'Heat capacity for each soil layer',                      &
         .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, runoffID, 'runoff', 'mm/timestep',          &
         'Total runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rnof1ID, 'rnof1', 'mm/timestep',            &
         'Surface runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rnof2ID, 'rnof2', 'mm/timestep',            &
         'Subsurface runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !---------------define snow states------------------------------------------
    CALL define_ovar(ncid_restart, tggsnID, 'tggsn', 'K',                      &
         'Average layer snow temperature',                         &
         .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ssdnnID, 'ssdnn', 'kg/m^3',                 &
         'Average snow density', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ssdnID, 'ssdn', 'kg/m^3',                   &
         'Average layer snow density',                             &
         .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, snowdID, 'snowd', 'mm',                     &
         'Liquid water eqivalent snow depth',                      &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, snageID, 'snage', '??',                     &
         'Snow age', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, smassID, 'smass', 'kg/m^2',                 &
         'Average layer snow mass',                                &
         .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, sdepthID, 'sdepth', 'm',                    &
         'Snow layer depth', .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, osnowdID, 'osnowd', 'mm',                   &
         'Previous time step snow depth in water equivalent',      &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, isflagID, 'isflag', '-',                    &
         'Snow layer scheme flag', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    !----------------define canopy states----------------------------------
    CALL define_ovar(ncid_restart, canstoID, 'cansto', 'mm',                   &
         'Canopy surface water storage', .TRUE., 'real', 0, 0, 0,  &
         mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ghfluxID, 'ghflux', 'W/m^2?',               &
         '????', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, sghfluxID, 'sghflux', 'W/m^2?',             &
         '????', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gaID, 'ga', 'W/m^2',                        &
         'Ground heat flux', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, dgdtgID, 'dgdtg', 'W/m^2/K',                &
         'Derivative of ground heat flux wrt soil temperature', .TRUE., &
         'r2', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fevID, 'fev', 'W/m^2',                      &
         'Latent heat flux from vegetation', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fesID, 'fes', 'W/m^2',                      &
         'Latent heat flux from soil',                             &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fhsID, 'fhs', 'W/m^2',                      &
         'Sensible heat flux from soil',                           &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !--------------biogeochemical variables------------------------
    CALL define_ovar(ncid_restart, cplantID, 'cplant', 'gC/m^2',               &
         'Plant carbon stores',                                    &
         .TRUE., plantcarbID, 'plantcarbon', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, csoilID, 'csoil', 'gC/m^2',                 &
         'Soil carbon stores',                                     &
         .TRUE., soilcarbID, 'soilcarbon', 0, 0, 0, mpID, dummy, .TRUE.)
    !-------------------others---------------------------------
    CALL define_ovar(ncid_restart, wbtot0ID, 'wbtot0', 'mm',                   &
         'Initial time step soil water total',                     &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, osnowd0ID, 'osnowd0', 'mm',                 &
         'Initial time step snow water total',                     &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, albedoID, 'albedo', '-',                    &
         'Albedo for shortwave and NIR radiation',                 &
         .TRUE., radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, tradID, 'trad', 'K',                        &
         'Surface radiative temperature (soil/snow/veg inclusive)', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !---------------------MODEL PARAMETERS---------------------------------
    WRITE(logn,'(A43)') '   Writing model parameters to restart file'
    CALL define_ovar(ncid_restart, rpid%iveg, 'iveg', '-',                     &
         'Vegetation type', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%isoil, 'isoil', '-',                   &
         'Soil type', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%clay, 'clay', '-',                     &
!$                     'Fraction of soil which is clay',                         &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%sand, 'sand', '-',                     &
!$                     'Fraction of soil which is sand',                         &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%silt, 'silt', '-',                     &
!$                     'Fraction of soil which is silt',                         &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%ssat, 'ssat', '-',                     &
!$                     'Fraction of soil volume which is water @ saturation',    &
!$                    .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%sfc, 'sfc', '-',                       &
!$                    'Fraction of soil volume which is water @ field capacity', &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%swilt, 'swilt', '-',                   &
!$                     'Fraction of soil volume which is water @ wilting point', &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! zse (depth of each soil layer):
    ok = NF90_DEF_VAR(ncid_restart, 'zse', NF90_FLOAT, (/soilID/), rpid%zse)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining zse variable in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%zse, "long_name",                     &
         "Depth of each soil layer")
    ok = NF90_PUT_ATT(ncid_restart, rpid%zse, "units", "m")
!$    CALL define_ovar(ncid_restart, rpid%froot, 'froot', '-',                   &
!$                     'Fraction of roots in each soil layer',                   &
!$                      .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%bch, 'bch', '-',                       &
!$                     'Parameter b, Campbell eqn 1985',                         &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%hyds, 'hyds', 'm/s',                   &
!$                     'Hydraulic conductivity @ saturation',                    &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%sucs, 'sucs', 'm',                     &
!$                     'Suction @ saturation', .TRUE.,                           &
!$                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%css, 'css', 'J/kg/C',                  &
!$                     'Heat capacity of soil minerals',                         &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%rhosoil, 'rhosoil', 'kg/m^3',          &
!$                     'Density of soil minerals',                               &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%rs20, 'rs20', '-',                     &
!$                     'Soil respiration coefficient at 20C',                    &
!$                      .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%albsoil, 'albsoil', '-',               &
         'Soil reflectance', .TRUE.,                               &
         radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%hc, 'hc', 'm',                         &
!$                     'Height of canopy', .TRUE.,                               &
!$                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%canst1, 'canst1', 'mm/LAI',            &
!$                     'Max water intercepted by canopy',                        &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%dleaf, 'dleaf', 'm',                   &
!$                     'Chararacteristic length of leaf',                        &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%frac4, 'frac4', '-',                   &
!$                     'Fraction of plants which are C4',                        &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%ejmax, 'ejmax', 'mol/m^2/s',           &
!$                     'Max potential electron transport rate top leaf', .TRUE., &
!$                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%vcmax, 'vcmax', 'mol/m^2/s',           &
!$                     'Maximum RuBP carboxylation rate top leaf', .TRUE.,       &
!$                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%rp20, 'rp20', '-',                     &
!$                     'Plant respiration coefficient at 20C', .TRUE., 'real',   &
!$                     0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%g0, 'g0', '-',                     &
!$                     'g0 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
!$                     0, 0, 0, mpID, dummy, .TRUE.) ! Ticket #56
!$    CALL define_ovar(ncid_restart, rpid%g1, 'g1', '-',                     &
!$                     'g1 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
!$                     0, 0, 0, mpID, dummy, .TRUE.)  ! Ticket #56
!$    CALL define_ovar(ncid_restart, rpid%rpcoef, 'rpcoef', '1/C',               &
!$                     'Temperature coef nonleaf plant respiration', .TRUE.,     &
!$                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%shelrb, 'shelrb', '-',                 &
!$              'Sheltering factor', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%xfang, 'xfang', '-',                   &
!$           'Leaf angle parameter', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%wai, 'wai', '-',                       &
!$                'Wood area index', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%vegcf, 'vegcf', '-',                   &
!$                     'vegcf', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%extkn, 'extkn', '-',                   &
!$                     'Extinction coef for vertical nitrogen profile',          &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%tminvj, 'tminvj', 'C',                 &
!$                     'Min temperature for the start of photosynthesis',        &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%tmaxvj, 'tmaxvj', 'C',                 &
!$                     'Max temperature for the start of photosynthesis',        &
!$                      .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%vbeta, 'vbeta', '-',                   &
!$                     'Stomatal sensitivity to soil water',                     &
!$                      .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%xalbnir, 'xalbnir', '-',               &
!$                     'modifier for albedo in near ir band',                    &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    ! ratecp (Plant carbon rate constant):
!$    ok = NF90_DEF_VAR(ncid_restart, 'ratecp', NF90_FLOAT, (/plantcarbID/),     &
!$                      rpid%ratecp)
!$    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
!$                     (ok, 'Error defining ratecp variable in restart file. '// &
!$                      '(SUBROUTINE create_restart)')
!$    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecp, "long_name",                  &
!$                      "Plant carbon rate constant")
!$    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecp, "units", "1/year")
!$    ! ratecs (Soil carbon rate constant):
!$    ok = NF90_DEF_VAR(ncid_restart, 'ratecs', NF90_FLOAT, (/soilcarbID/),      &
!$                      rpid%ratecs)
!$    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
!$                     (ok, 'Error defining ratecs variable in restart file. '// &
!$                      '(SUBROUTINE create_restart)')
!$    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecs, "long_name",                  &
!$                      "Soil carbon rate constant")
!$    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecs, "units", "1/year")
!$    CALL define_ovar(ncid_restart, rpid%meth, 'meth', '-',                     &
!$                     'Canopy turbulence parameterisation switch',              &
!$                     .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%za_uv, 'za_uv', 'm',                   &
!$                    'Reference height (lowest atm. model layer) for momentum', &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
!$    CALL define_ovar(ncid_restart, rpid%za_tq, 'za_tq', 'm',                   &
!$                     'Reference height (lowest atm. model layer) for scalars', &
!$                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gwID, 'GWwb', 'mm3/mm3','GW water content', &
         .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)

!$    IF(cable_user%SOIL_STRUC=='sli'.OR.cable_user%FWSOIL_SWITCH=='Haverd2013') THEN
!$      CALL define_ovar(ncid_restart,rpid%gamma,'gamma','-', &
!$            'Parameter in root efficiency function (Lai and Katul 2000)', &
!$            .TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
!$    ENDIF
    ! Soil-Litter-Iso soil model
    IF(cable_user%SOIL_STRUC=='sli') THEN
       ! Parameters for SLI:
!$       CALL define_ovar(ncid_restart,rpid%nhorizons,'nhorizons','-', &
!$            'Number of soil horizons',.TRUE.,'integer',0,0,0,mpID,dummy,.TRUE.)
!$       CALL define_ovar(ncid_restart,rpid%zeta,'zeta','[ ]', &
!$            'exponent factor in Topmodel eq',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
!$       CALL define_ovar(ncid_restart,rpid%fsatmax,'fsatmax','[ ]', &
!$            'param in Topmodel eq',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
!$       CALL define_ovar(ncid_restart,rpid%ishorizon,'ishorizon','-', &
!$            'Horizon number',.TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
!$       CALL define_ovar(ncid_restart,rpid%clitt,'clitt','tC/ha', &
!$            'Litter layer carbon content',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
!$       CALL define_ovar(ncid_restart,rpid%ZR,'ZR','cm', &
!$            'Maximum rooting depth',.TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
!$       CALL define_ovar(ncid_restart,rpid%F10,'F10','-', &
!$            'Fraction of roots in top 10 cm', &
!$            .TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       ! Variables for SLI:
       CALL define_ovar(ncid_restart,SID,'S','-',&
            'Fractional soil moisture content relative to saturated value', &
            .TRUE.,soilID,'soil',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart,TsoilID,'Tsoil','degC',&
            'Tsoil', &
            .TRUE.,soilID,'soil',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart,snowliqID,'snowliq','mm',&
            'liquid water content of snowpack', &
            .TRUE.,snowID,'snow',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart,scondsID,'sconds','Wm-1K-1',&
            'thermal cond of snowpack', &
            .TRUE.,snowID,'snow',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart,h0ID,'h0','m',&
            'Pond height above soil', &
            .TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart,nsnowID,'nsnow','-',&
            'number of snow layers', &
            .TRUE.,'integer',0,0,0,mpID,dummy,.TRUE.)
       CALL define_ovar(ncid_restart,TsurfaceID,'Tsurface','degC',&
            'soil or snow surface T', &
            .TRUE.,'real',0,0,0,mpID,dummy,.TRUE.)
    END IF ! SLI soil model

    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Production",                 &
         TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(frst_out)// ' (SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Source",                     &
         'CABLE LSM restart file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(frst_out)// ' (SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "CABLE_input_file",           &
         TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(frst_out)// ' (SUBROUTINE create_restart)')

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file '      &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write time variable:
    ok = NF90_PUT_VAR(ncid_restart, tvarID, REAL(REAL(ktau) * dels, r_2))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error time variable to '           &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write latitude and longitude variables:
    ok = NF90_PUT_VAR(ncid_restart, latID, latitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing latitude variable to '   &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, lonID, longitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing longitude variable to '  &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write number of active patches for each land grid cell:
    ok = NF90_PUT_VAR(ncid_restart, napID, landpt(:)%nap)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing nap variable to '        &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write vegetated patch fractions
    ok = NF90_PUT_VAR(ncid_restart, rpid%patchfrac,                            &
         patch(:)%frac, start = (/1/), count = (/mp/))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing patchfrac to '       &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write number of veg and soil types
    ok = NF90_PUT_VAR(ncid_restart, rpid%mvtype,mvtype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing mvtype parameter to '    &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rpid%mstype,mstype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
         'Error writing mstype parameter to '    &
         //TRIM(frst_out)// '(SUBROUTINE create_restart)')

    ! Write parameters:
    !~ veg and soil
    out_settings%dimswitch = "integer"
    CALL check_and_write(rpid%iveg, &
                         'iveg', REAL(veg%iveg, 4), ranges%iveg, patchout_var, out_settings)
    CALL check_and_write(rpid%isoil, 'isoil', REAL(soil%isoilm, 4), &
                         ranges%isoil, patchout_var, out_settings)
    out_settings%dimswitch = "real"
!$  CALL check_and_write(rpid%bch, 'bch', REAL(soil%bch, 4), &
!$                       ranges%bch, patchout_var, out_settings)
!$  CALL check_and_write(rpid%bch, 'bch', REAL(soil%bch, 4), &
!$                       ranges%bch, patchout_var, out_settings)
!$  CALL check_and_write(rpid%clay, 'clay', REAL(soil%clay, 4), &
!$                       ranges%clay, patchout_var, out_settings)
!$  CALL check_and_write(rpid%sand, 'sand', REAL(soil%sand, 4), &
!$                       ranges%sand, patchout_var, out_settings)
!$  CALL check_and_write(rpid%silt, 'silt', REAL(soil%silt, 4), &
!$                       ranges%silt, patchout_var, out_settings)
!$  CALL check_and_write(rpid%css, 'css', REAL(soil%css, 4), &
!$                       ranges%css, patchout_var, out_settings)
!$  CALL check_and_write(rpid%rhosoil, 'rhosoil', &
!$                       REAL(soil%rhosoil, 4), ranges%rhosoil, patchout_var, &
!$                       out_settings)
!$  CALL check_and_write(rpid%hyds, 'hyds', REAL(soil%hyds, 4), &
!$                       ranges%hyds, patchout_var, out_settings)
!$  CALL check_and_write(rpid%sucs, 'sucs', REAL(soil%sucs, 4), &
!$                       ranges%sucs, patchout_var, out_settings)
!$  CALL check_and_write(rpid%rs20, 'rs20', REAL(veg%rs20, 4), &
!$                       ranges%rs20, patchout_var, out_settings)
!$  CALL check_and_write(rpid%ssat, 'ssat', REAL(soil%ssat, 4), &
!$                       ranges%ssat, patchout_var, out_settings)
!$  CALL check_and_write(rpid%sfc, 'sfc', REAL(soil%sfc, 4), &
!$                       ranges%sfc, patchout_var, out_settings)
!$  CALL check_and_write(rpid%swilt, 'swilt', REAL(soil%swilt, 4), &
!$                       ranges%swilt, patchout_var, out_settings)
    ! Soil dimensioned variables/parameters:
    out_settings%dimswitch = "soil"
!$  CALL check_and_write(rpid%froot, 'froot', REAL(veg%froot, 4), &
!$                       ranges%froot, patchout_var, out_settings)

    !~ ssnow
    !~~ Soil dimensioned variables/parameters:
    out_settings%dimswitch = "soil"
    CALL check_and_write(tggID, 'tgg', REAL(ssnow%tgg, 4), &
                         ranges%SoilTemp, patchout_var, out_settings)
    CALL check_and_write(wbID, 'wb', REAL(ssnow%wb, 4), ranges%SoilMoist, &
                         patchout_var, out_settings)
    CALL check_and_write(wbiceID, 'wbice', REAL(ssnow%wbice, 4), &
                         ranges%SoilMoist, patchout_var, out_settings)
    CALL check_and_write(gammzzID, 'gammzz', REAL(ssnow%gammzz, 4), &
                         ranges%default_l, patchout_var, out_settings)
    !~~ Snow dimensioned variables/parameters:
    out_settings%dimswitch = "snow"
    CALL check_and_write(ssdnID, 'ssdn', REAL(ssnow%ssdn, 4), &
                         ranges%ssdn, patchout_var, out_settings)
    CALL check_and_write(smassID, 'smass', REAL(ssnow%smass, 4), &
                         ranges%smass, patchout_var, out_settings)
    CALL check_and_write(sdepthID, 'sdepth', REAL(ssnow%sdepth, 4), &
                         ranges%sdepth, patchout_var, out_settings)
    CALL check_and_write(tggsnID, 'tggsn', REAL(ssnow%tggsn, 4), &
                         ranges%tggsn, patchout_var, out_settings)
    !~~ Other dims
    out_settings%dimswitch = "radiation"
    CALL check_and_write(albsoilsnID, 'albsoilsn', &
                         REAL(ssnow%albsoilsn, 4), ranges%albsoiln, patchout_var, out_settings)
    out_settings%dimswitch = "plantcarbon"
    CALL check_and_write(cplantID, 'cplant', REAL(bgc%cplant, 4), &
                         ranges%default_l, patchout_var, out_settings)
    out_settings%dimswitch = "soilcarbon"
    CALL check_and_write(csoilID, 'csoil', REAL(bgc%csoil, 4), &
                         ranges%default_l, patchout_var, out_settings)
    ok = NF90_PUT_VAR(ncid_restart, rpid%zse, REAL(soil%zse, 4))
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing zse parameter to ' &
                                        //TRIM(frst_out)//'(SUBROUTINE create_restart)')
    ! Single dim:
    out_settings%dimswitch = "radiation"
    CALL check_and_write(rpid%albsoil, 'albsoil', &
                         REAL(soil%albsoil, 4), ranges%albsoil, patchout_var, out_settings)
!$  out_settings%dimswitch = "real"
!$  CALL check_and_write(rpid%canst1, 'canst1', REAL(veg%canst1, 4), &
!$                       ranges%canst1, patchout_var, out_settings)
!$  CALL check_and_write(rpid%dleaf, 'dleaf', REAL(veg%dleaf, 4), &
!$                       ranges%dleaf, patchout_var, out_settings)
!$  CALL check_and_write(rpid%ejmax, 'ejmax', REAL(veg%ejmax, 4), &
!$                       ranges%ejmax, patchout_var, out_settings)
!$  CALL check_and_write(rpid%vcmax, 'vcmax', REAL(veg%vcmax, 4), &
!$                       ranges%vcmax, patchout_var, out_settings)
!$  CALL check_and_write(rpid%frac4, 'frac4', REAL(veg%frac4, 4), &
!$                       ranges%frac4, patchout_var, out_settings)
!$  CALL check_and_write(rpid%hc, 'hc', REAL(veg%hc, 4), &
!$                       ranges%hc, patchout_var, out_settings)
!$  CALL check_and_write(rpid%rp20, 'rp20', REAL(veg%rp20, 4), &
!$                       ranges%rp20, patchout_var, out_settings)
!$  CALL check_and_write(rpid%g0, 'g0', REAL(veg%g0, 4), &
!$                       ranges%g0, patchout_var, out_settings) ! Ticket #56
!$  CALL check_and_write(rpid%g1, 'g1', REAL(veg%g1, 4), &
!$                       ranges%g1, patchout_var, out_settings) ! Ticket #56
!$  CALL check_and_write(rpid%rpcoef, 'rpcoef', REAL(veg%rpcoef, 4), &
!$                       ranges%rpcoef, patchout_var, out_settings)
!$  CALL check_and_write(rpid%shelrb, 'shelrb', REAL(veg%shelrb, 4), &
!$                       ranges%shelrb, patchout_var, out_settings)
!$  CALL check_and_write(rpid%xfang, 'xfang', REAL(veg%xfang, 4), &
!$                       ranges%xfang, patchout_var, out_settings)
!$  CALL check_and_write(rpid%wai, 'wai', REAL(veg%wai, 4), &
!$                       ranges%wai, patchout_var, out_settings)
!$  CALL check_and_write(rpid%vegcf, 'vegcf', REAL(veg%vegcf, 4), &
!$                       ranges%vegcf, patchout_var, out_settings)
!$  CALL check_and_write(rpid%extkn, 'extkn', REAL(veg%extkn, 4), &
!$                       ranges%extkn, patchout_var, out_settings)
!$  CALL check_and_write(rpid%tminvj, 'tminvj', REAL(veg%tminvj, 4), &
!$                       ranges%tminvj, patchout_var, out_settings)
!$  CALL check_and_write(rpid%tmaxvj, 'tmaxvj', REAL(veg%tmaxvj, 4), &
!$                       ranges%tmaxvj, patchout_var, out_settings)
!$  CALL check_and_write(rpid%vbeta, 'vbeta', REAL(veg%vbeta, 4), &
!$                       ranges%vbeta, patchout_var, out_settings)
!$  CALL check_and_write(rpid%xalbnir, 'xalbnir', &
!$                       REAL(veg%xalbnir, 4), ranges%xalbnir, patchout_var, &
!$                       out_settings)
!$  CALL check_and_write(rpid%tmaxvj, 'tmaxvj', REAL(veg%tmaxvj, 4), &
!$                       ranges%tmaxvj, patchout_var, out_settings)
!$  ok = NF90_PUT_VAR(rpid%ratecp, REAL(bgc%ratecp, 4))
!$  IF (ok /= NF90_NOERR) CALL nc_abort(ok, &
!$                                      'Error writing ratecp parameter to ' &
!$                                      //TRIM(frst_out)//'(SUBROUTINE create_restart)')
!$  ok = NF90_PUT_VAR(rpid%ratecs, REAL(bgc%ratecs, 4))
!$  IF (ok /= NF90_NOERR) CALL nc_abort(ok, &
!$                                      'Error writing ratecs parameter to ' &
!$                                      //TRIM(frst_out)//'(SUBROUTINE create_restart)')
!$  out_settings%dimswitch = "integer"
!$  CALL check_and_write(rpid%meth, 'meth', REAL(veg%meth, 4), &
!$                       ranges%meth, patchout_var, out_settings)
!$  out_settings%dimswitch = "real"
!$  CALL check_and_write(rpid%za_uv, 'za_uv', REAL(rough%za_uv, 4), &
!$                       ranges%za, patchout_var, out_settings)
!$  CALL check_and_write(rpid%za_tq, 'za_tq', REAL(rough%za_tq, 4), &
!$                       ranges%za, patchout_var, out_settings)
    out_settings%dimswitch = "r2"
    CALL check_and_write(dgdtgID, 'dgdtg', REAL(canopy%dgdtg, 4), &
                         ranges%default_l, patchout_var, out_settings)
    out_settings%dimswitch = "integer"
    CALL check_and_write(isflagID, 'isflag', REAL(ssnow%isflag, 4), &
                         ranges%default_l, patchout_var, out_settings)
    out_settings%dimswitch = "real"
    CALL check_and_write(gwID, 'GWwb', REAL(ssnow%GWwb, 4), &
                         ranges%GWwb, patchout_var, out_settings)
    CALL check_and_write(tssID, 'tss', REAL(ssnow%tss, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(ssdnnID, 'ssdnn', REAL(ssnow%ssdnn, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(osnowdID, 'osnowd', REAL(ssnow%osnowd, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(snageID, 'snage', REAL(ssnow%snage, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(snowdID, 'snowd', REAL(ssnow%snowd, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(rtsoilID, 'rtsoil', REAL(ssnow%rtsoil, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(canstoID, 'cansto', REAL(canopy%cansto, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(sghfluxID, 'sghflux', &
                         REAL(canopy%sghflux, 4), ranges%default_l, &
                         patchout_var, out_settings)
    CALL check_and_write(ghfluxID, 'ghflux', REAL(canopy%ghflux, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(runoffID, 'runoff', REAL(ssnow%runoff, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(rnof1ID, 'rnof1', REAL(ssnow%rnof1, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(rnof2ID, 'rnof2', REAL(ssnow%rnof2, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(gaID, 'ga', REAL(canopy%ga, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(fevID, 'fev', REAL(canopy%fev, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(fesID, 'fes', REAL(canopy%fes, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(fhsID, 'fhs', REAL(canopy%fhs, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(wbtot0ID, 'wbtot0', REAL(bal%wbtot0, 4), &
                         ranges%default_l, patchout_var, out_settings)
    CALL check_and_write(osnowd0ID, 'osnowd0', REAL(bal%osnowd0, 4), &
                         ranges%default_l, patchout_var, out_settings)

    !~ Radiation
    out_settings%dimswitch = "radiation"
    CALL check_and_write(albedoID, 'albedo', REAL(rad%albedo, 4), &
                         ranges%Albedo, patchout_var, out_settings)
    out_settings%dimswitch = "real"
    CALL check_and_write(tradID, 'trad', &
                         REAL(rad%trad, 4), ranges%RadT, patchout_var, out_settings)

!$  IF (cable_user%SOIL_STRUC == 'sli' .OR. cable_user%FWSOIL_SWITCH == 'Haverd2013') THEN
!$    CALL check_and_write(rpid%gamma, 'gamma', &
!$                         REAL(veg%gamma, 4), ranges%default_s, patchout_var, out_settings)
!$  END IF
!$

    IF (cable_user%SOIL_STRUC == 'sli') THEN
      ! Write SLI parameters:
      out_settings%dimswitch = "integer"
!$    CALL check_and_write(rpid%nhorizons, 'nhorizons', &
!$                         REAL(soil%nhorizons, 4), ranges%default_s, patchout_var, out_settings)
      CALL check_and_write(nsnowID, 'nsnow', REAL(ssnow%nsnow, 4), &
                           ranges%default_s, patchout_var, out_settings)
      out_settings%dimswitch = "soil"
!$    CALL check_and_write(rpid%ishorizon, 'ishorizon', &
!$                         REAL(soil%ishorizon, 4), ranges%default_s, patchout_var, out_settings)
      CALL check_and_write(SID, 'S', REAL(ssnow%S, 4), &
                         ranges%S, patchout_var, out_settings)
      CALL check_and_write(TsoilID, 'Tsoil', REAL(ssnow%Tsoil, 4), &
                         ranges%Tsoil, patchout_var, out_settings)
      out_settings%dimswitch = "real"
!$    CALL check_and_write(rpid%clitt, 'clitt', &
!$                         REAL(veg%clitt, 4), ranges%default_s, patchout_var, out_settings)
!$    CALL check_and_write(rpid%ZR, 'ZR', &
!$                         REAL(veg%ZR, 4), ranges%default_s, patchout_var, out_settings)
!$    CALL check_and_write(rpid%F10, 'F10', &
!$                         REAL(veg%F10, 4), ranges%default_s, patchout_var, out_settings)
      CALL check_and_write(TsurfaceID, 'Tsurface', REAL(ssnow%Tsurface, 4), &
                           ranges%default_s, patchout_var, out_settings)
      CALL check_and_write(h0ID, 'h0', REAL(ssnow%h0, 4), &
                           ranges%default_s, patchout_var, out_settings)
      out_settings%dimswitch = "snow"
      CALL check_and_write(snowliqID, 'snowliq', REAL(ssnow%snowliq, 4), &
                           ranges%default_s, patchout_var, out_settings)
      CALL check_and_write(scondsID, 'sconds', REAL(ssnow%sconds, 4), &
                           ranges%default_s, patchout_var, out_settings)
    END IF

    ! Close restart file
    ok = NF90_CLOSE(ncid_restart)

    WRITE (logn, '(A36)') '   Restart file complete and closed.'

  END SUBROUTINE create_restart


END MODULE cable_output_module
