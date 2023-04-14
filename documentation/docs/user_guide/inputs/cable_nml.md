
# Namelist file cable.nml #

The namelist file for CABLE is cable.nml. For offline applications it
should be located in the directory in which you are running CABLE. For
ACCESS cases, cable.nml should be located in ~/CABLE-AUX/UM on the
machine that CABLE is being executed on.

The cable.nml file includes some settings that are common across all CABLE
applications and others that are required only for offline
applications. The following are annotated examples of cable.nml:

!!! Note

    If your namelist file has entries that do not appear here please check the
    [obsolete and deprecated user features](../obsolete_and_deprecated_features/obsolete_and_deprecated_features.md) page.


## For offline applications ##

``` fortran-free-form
&cable
   filename%met = '/projects/access/CABLE-AUX/offline/TumbaFluxnet.1.3_met.nc'
   filename%out = 'out_cable.nc'
   filename%log = 'log_cable.txt'
   filename%restart_in  = ' '
   filename%restart_out = './restart_out.nc'
   filename%type    = '/projects/access/CABLE-AUX/offline/gridinfo_CSIRO_1x1.nc'
   vegparmnew = .TRUE.  ! using new format when true
   soilparmnew = .TRUE.  ! using new format when true
   spinup = .TRUE.  ! do we spin up the model?
   delsoilM = 0.001   ! allowed variation in soil moisture for spin up
   delsoilT = 0.01    ! allowed variation in soil temperature for spin up
   output%restart = .TRUE.  ! should a restart file be created?
   output%met = .TRUE.  ! input met data
   output%flux = .TRUE.  ! convective, runoff, NEE
   output%soil = .TRUE.  ! soil states
   output%snow = .TRUE.  ! snow states
   output%radiation = .TRUE.  ! net rad, albedo
   output%carbon    = .TRUE.  ! NEE, GPP, NPP, stores
   output%veg       = .TRUE.  ! vegetation states
   output%params    = .TRUE.  ! input parameters used to produce run
   output%balances  = .TRUE.  ! energy and water balances
   check%ranges     = .FALSE.  ! variable ranges, input and output
   check%energy_bal = .TRUE.  ! energy balance
   check%mass_bal   = .TRUE.  ! water/mass balance
   verbose = .TRUE. ! write details of every grid cell init and params to log?
   leaps = .TRUE. ! calculate timing with leap years?
   logn = 88      ! log file number - declared in input module
   fixedCO2 = 350.0   ! if not found in met file, in ppmv
   spincasa      = .FALSE.     ! spin casa before running the model if TRUE, and should be set to FALSE if spincasainput = .TRUE.
   l_casacnp     = .FALSE.  ! using casaCNP with CABLE
   l_laiFeedbk   = .FALSE.  ! using prognostic LAI
   l_vcmaxFeedbk = .FALSE.  ! using prognostic Vcmax
   icycle = 0   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P
   casafile%cnpipool='/projects/access/CABLE-AUX/core/biogeochem/poolcnpInTumbarumba.csv'  !
   casafile%cnpbiome='/projects/access/CABLE-AUX/core/biogeochem/pftlookup_csiro_v16_17tiles.csv'  ! biome specific BGC parameters
   casafile%cnpepool='poolcnpOut.csv'    ! end of run pool size
   casafile%cnpmetout='casamet.nc'                ! output daily met forcing for spinning casacnp
   casafile%cnpmetin=''          ! list of daily met files for spinning casacnp
   casafile%phen='/projects/access/CABLE-AUX/core/biogeochem/modis_phenology_csiro.txt'        ! modis phenology
   casafile%cnpflux='cnpfluxOut.csv'
   ncciy = 0 ! 0 for not using gswp; 4-digit year input for year of gswp met
   gswpfile%rainf = 'gswp/Rainf_gswp1987.nc'
   gswpfile%snowf = 'gswp/Snowf_gswp1987.nc'
   gswpfile%LWdown= 'gswp/LWdown_srb1987.nc'
   gswpfile%SWdown= 'gswp/SWdown_srb1987.nc'
   gswpfile%PSurf = 'gswp/PSurf_ecor1987.nc'
   gswpfile%Qair  = 'gswp/Qair_cru1987.nc'
   gswpfile%Tair  = 'gswp/Tair_cru1987.nc'
   gswpfile%wind  = 'gswp/Wind_ncep1987.nc'
   redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
   wiltParam = 0.5
   satuParam = 0.8
   cable_user%FWSOIL_SWITCH = 'standard'        ! choices are: 'standard', 'non-linear extrapolation' or 'Lai and Ktaul 2000'
   cable_user%DIAG_SOIL_RESP = 'ON '
   cable_user%LEAF_RESPIRATION = 'ON '
   cable_user%RUN_DIAG_LEVEL= 'BASIC'        ! choices are: 'BASIC', 'NONE'
   cable_user%CONSISTENCY_CHECK= .TRUE.      ! TRUE outputs combined fluxes at each timestep for comparison to a control run
   cable_user%CASA_DUMP_READ = .FALSE.      ! TRUE reads CASA forcing from netcdf format
   cable_user%CASA_DUMP_WRITE = .FALSE.      ! TRUE outputs CASA forcing in netcdf format
   cable_user%SSNOW_POTEV= 'HDM'      ! Humidity Deficit Method
&end
```

| Namelist variable            | Type      | Available values                    | Default value                                                                | Description                                                                                              |
|------------------------------|-----------|-------------------------------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| filename%met                 | character(len=500) | any string of max. 500 characters | '/projects/access/CABLE-AUX/offline/TumbaFluxnet.1.3_met.nc'          | Meteorological forcing file                                                                              |
| filename%out                 | character(len=500) | any string of max. 500 characters | 'out_cable.nc'                                                        | Output file                                                                                              |
| filename%log                 | character(len=500) | any string of max. 500 characters | 'log_cable.txt'                                                       | Log file                                                                                                 |
| filename%restart_in          | character(len=500) | any string of max. 500 characters | ' '                                                                   | Input restart file                                                                                       |
| filename%restart_out         | character(len=500) | any string of max. 500 characters | './restart_out.nc'                                                    | Output restart file                                                                                      |
| filename%type                | character(len=500) | any string of max. 500 characters | '/projects/access/CABLE-AUX/offline/gridinfo_CSIRO_1x1.nc'            | Static information file about grid, soil and vegetation                                                  |
| vegparmnew                   | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Use new format for vegetation parameter files when true                                                  |
| soilparmnew                  | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Use new format for soil parameter files when true                                                        |
| spinup                       | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Spin up the model when .TRUE.                                                                            |
| delsoilM                     | real      |  any real number                    | 0.001                                                                        | Allowed variation in soil moisture for spin up                                                           |
| delsoilT                     | real      |  any real number                    | 0.01                                                                         | Allowed variation in soil temperature for spin up                                                        |
| output%restart               | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Create a restart file when .TRUE.                                                                        |
| output%met                   | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Use input meteorological data when .TRUE.                                                                |
| output%flux                  | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Convective, runoff, NEE                                                                                  |
| output%soil                  | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Use soil states when .TRUE.                                                                              |
| output%snow                  | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Use snow states when .TRUE.                                                                              |
| output%radiation             | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Net rad, albedo                                                                                          |
| output%carbon                | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | NEE, GPP, NPP, stores                                                                                    |
| output%veg                   | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Vegetation states                                                                                        |
| output%params                | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Input parameters used to produce run                                                                     |
| output%balances              | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Energy and water balances                                                                                |
| check%ranges                 | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | Variable ranges, input and output                                                                        |
| check%energy_bal             | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Energy balance                                                                                           |
| check%mass_bal               | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Water/mass balance                                                                                       |
| verbose                      | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Write details of every grid cell initialisation and parameters to log.                                   |
| leaps                        | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | Calculate timing with leap years?                                                                        |
| logn                         | integer   | any integer number                  | 88                                                                           | Log file number - declared in input module                                                               |
| fixedCO2                     | real      | any real number                     | 350                                                                          | CO2 concentration if not found in met file, in ppmv                                                      |
| spincasa                     | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | Spin CASA before running the model if TRUE, and should be set to FALSE if spincasainput                  |
| l_casacnp                    | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | Use CASA-CNP with CABLE                                                                                  |
| l_laiFeedbk                  | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | Use prognostic LAI                                                                                       |
| l_vcmaxFeedbk                | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | Use prognostic Vcmax                                                                                     |
| icycle                       | integer   | 0, 1, 2, 3, >10                     | 0                                                                            | 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P, >10 excludes call to CBM (main call to cable)  |
| casafile%cnpipool            | character(len=99) | any string of max. 99 characters | '/projects/access/CABLE-AUX/core/biogeochem/poolcnpInTumbarumba.csv'    | Initial CNP pools file.                                                                                  |
| casafile%cnpbiome            | character(len=99) | any string of max. 99 characters | '/projects/access/CABLE-AUX/core/biogeochem/pftlookup_csiro_v16_17tiles.csv' | Biome specific biogeochemical parameters                                                            |
| casafile%cnpepool            | character(len=99) | any string of max. 99 characters | 'poolcnpOut.csv'                                                        | End of run pool size                                                                                     |
| casafile%cnpmetout           | character(len=99) | any string of max. 99 characters | 'casamet.nc'                                                            | Output daily met forcing for spinning CASA-CNP                                                           |
| casafile%cnpmetin            | character(len=99) | any string of max. 99 characters | ''                                                                      | List of daily met files for spinning CASA-CNP                                                            |
| casafile%phen                | character(len=99) | any string of max. 99 characters | '/projects/access/CABLE-AUX/core/biogeochem/modis_phenology_csiro.txt'  | Modis phenology                                                                                          |
| casafile%cnpflux             | character(len=99) | any string of max. 99 characters | 'cnpfluxOut.csv'                                                        | Output file for CNP fluxes                                                                               |
| ncciy                        | integer   |                                     | 0                                                                            | GSWP year. 0 for not using GSWP; 4-digit year input for year of GSWP meteorological forcing.             |
| gswpfile%rainf               | character(len=200) | any string of max. 200 characters | 'gswp/Rainf_gswp1987.nc'                                              | Input file for rain meteorological forcing.                                                              |
| gswpfile%snowf               | character(len=200) | any string of max. 200 characters | 'gswp/Snowf_gswp1987.nc'                                              | Input file for snow meteorological forcing.                                                              |
| gswpfile%LWdown              | character(len=200) | any string of max. 200 characters | 'gswp/LWdown_srb1987.nc'                                              | Input file for long wave radiation meteorological forcing.                                               |
| gswpfile%SWdown              | character(len=200) | any string of max. 200 characters | 'gswp/SWdown_srb1987.nc'                                              | Input file for short wave radiation meteorological forcing.                                              |
| gswpfile%PSurf               | character(len=200) | any string of max. 200 characters | 'gswp/PSurf_ecor1987.nc'                                              | Input file for surface pressure meteorological forcing.                                                  |
| gswpfile%Qair                | character(len=200) | any string of max. 200 characters | 'gswp/Qair_cru1987.nc'                                                | Input file for humidity meteorological forcing.                                                          |
| gswpfile%Tair                | character(len=200) | any string of max. 200 characters | 'gswp/Tair_cru1987.nc'                                                | Input file for temperature meteorological forcing.                                                       |
| gswpfile%wind                | character(len=200) | any string of max. 200 characters | 'gswp/Wind_ncep1987.nc'                                               | Input file for wind meteorological forcing.                                                              |
| redistrb                     | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | Turn on/off the hydraulic redistribution.                                                                |
| wiltParam                    | real      |                                     | 0.5                                                                          | Wilt parameter in hydraulic redistribution module.                                                       |
| satuParam                    | real      |                                     | 0.8                                                                          | Saturation parameter in hydraulic redistribution module.                                                 |
| cable_user%FWSOIL_SWITCH     | character | 'standard' 'non-linear extrapolation' 'Lai and Ktaul 2000' | 'standard'                                            | Controls root water uptake function.                                                                     |
| cable_user%DIAG_SOIL_RESP    | character(len=3) | 'ON ' 'OFF'                  | 'ON '                                                                        | Controls soil respiration scheme when CASA-CNP not used.                                                 |
| cable_user%LEAF_RESPIRATION  | character(len=3) | 'ON ' 'OFF'                  | 'ON '                                                                        | Controls what is output into the GPP stash variable.                                                     |
| cable_user%RUN_DIAG_LEVEL    | character | 'BASIC' 'NONE'                      | 'BASIC'                                                                      | Controls output from CABLE to standard out.                                                              |
| cable_user%CONSISTENCY_CHECK | logical   | .TRUE. .FALSE.                      | .TRUE.                                                                       | TRUE outputs combined fluxes at each timestep for comparison to a control run                            |
| cable_user%CASA_DUMP_READ    | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | TRUE reads CASA forcing from netcdf format                                                               |
| cable_user%CASA_DUMP_WRITE   | logical   | .TRUE. .FALSE.                      | .FALSE.                                                                      | TRUE outputs CASA forcing in netcdf format                                                               |
| cable_user%SSNOW_POTEV       | character(len=3) | 'HDM' 'P-M'                  | 'HDM'                                                                        | Use the "Humidity Deficit Method" or  "Penman-Monteith" method                                           |


## For ACCESS applications ##


``` fortran-free-form
&cable
!!
!! cable namelist variables applicable to ACCESS/UM runs
!! default parameter files are available in the following directory on raijin
!! vegetation parameters in def_veg_params.txt or for vbeta=1 veg_params_vbeta1.txt
!! vbeta=1 recommended but should be tested for your applications
!
  filename%veg     = '/projects/access/CABLE-AUX/core/biogeophys/veg_params_vbeta1.txt'
!
!! or checkout to your own area
! filename%veg     = '~/CABLE-AUX/core/biogeophys/def_veg_params.txt' !relative to your home dir
!!
!! soil parameters in def_soil_params.txt. Not used except for isoil=9 (permanent ice) and
!! rhosoil and css values (isoil=2 used for all non-ice tiles).
!! Other soil parameters from spatially explicit ancillary files.
!
  filename%soil    = '/projects/access/CABLE-AUX/core/biogeophys/def_soil_params.txt'
!
! filename%soil    = '~/CABLE-AUX/core/biogeophys/def_soil_params.txt' !relative to your home dir
!
!! cable_user flags
!
  cable_user%CABLE_RUNTIME_COUPLED = .TRUE. ! controls whether snow initialised if
                                            ! ktau_gl=1
                                            ! must be set true for coupled run
                                            ! use true for atmosphere only runs if
                                            ! want to take snow from start dump
                                            ! default is .false.
  cable_user%FWSOIL_SWITCH = 'standard' ! Controls root water uptake function
                                        ! choices are:
                                        ! 'standard'
                                        ! 'non-linear extrapolation'
                                        ! 'Lai and Ktaul 2000' NB TYPO IS IN CODE (should be Katul)
  cable_user%DIAG_SOIL_RESP = 'ON '  ! Controls soil respiration scheme when CASA-CNP not used
                                     ! 'ON' uses scheme tuned with parameter veg%vegcf
                                     ! 'OFF' uses scheme tuned with parameter veg%rs20
  cable_user%LEAF_RESPIRATION = 'OFF' ! Controls what is output into the GPP stash variable
                                      ! 'OFF' is true GPP (photosynthesis + leaf respiration)
                                      ! 'ON' is only photosynthesis, canopy%fpn
  cable_user%RUN_DIAG_LEVEL= 'BASIC'  ! Controls output from CABLE to standard out
                                      ! choices are:
                                      ! 'BASIC'
                                      ! 'NONE'
                                      ! 'zero' (in cable_explicit_driver - what does this do?)
  cable_user%ssnow_POTEV = ''    ! Controls method used for soil evaporation
                                 ! Choices are:
                                 ! 'P-M' Penman_Monteith, otherwise defaults to
                                 ! humidity deficit method
!!--------------Test Switches--------------
  cable_user%l_new_roughness_soil  = .FALSE.   ! for new soil roughness, E.Kowalczyk Mar14
  cable_user%l_new_runoff_speed    = .FALSE.   ! for new increase in runoff speed, E.Kowalczyk Mar14
  cable_user%l_new_reduce_soilevp  = .FALSE.   ! to reduce soil evaporation, E.Kowalczyk Mar14
!!-----------------------------------------
!
  redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
  wiltParam = 0.5     ! only used if hydraulic redistribution .true.
  satuParam = 0.8     ! only used if hydraulic redistribution .true.
!
! CASA-CNP flags
  l_casacnp = .true. ! redundant? actually controlled by icycle
  icycle = 3 ! Controls whether CASA-CNP is run and for which species
             ! Choices are:
             ! 0 CASA-CNP not run
             ! 1 Carbon only
             ! 2 Carbon and nitrogen
             ! 3 Carbon, nitrogen and phosphorus
  l_laiFeedbk   = .TRUE.  ! using prognostic LAI
  l_vcmaxFeedbk = .TRUE.  ! using prognostic Vcmax
! filenames for CASA-CNP input files
  casafile%cnpbiome='~/CABLE-AUX/core/biogeochem/pftlookup_csiro_v16_17tiles.csv' ! biome specific BGC parameters
  casafile%phen='~/CABLE-AUX/core/biogeochem/modis_phenology_csiro.txt' ! phenology by latitude (modis derived)
&end
```
