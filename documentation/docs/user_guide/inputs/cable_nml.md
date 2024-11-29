
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
    [obsolete and deprecated user features][obsolete] page.


## List of namelist variables

| Namelist variable                | Type               | Available values                                           | Default value                    | Description                                                                                             |
|----------------------------------|--------------------|------------------------------------------------------------|----------------------------------|---------------------------------------------------------------------------------------------------------|
| filename%met                     | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Meteorological forcing file name.                                                                       |
| filename%path                    | character(len=500) | any string of max. 500 characters                          | './'                             | Path for output and restart files for CABLE and CASA.                                                   |
| filename%out                     | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Output file name.                                                                                       |
| filename%log                     | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Execution log file name.                                                                                |
| filename%restart_in              | character(len=500) | any string of max. 500 characters                          | ' '                              | Input restart file name. `''` will first search for `./<RunIden>_<year>_cable_rst.nc` and if not found, it uses default parameters. |
| filename%restart_out             | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Output restart file name.                                                                               |
| filename%type                    | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Static information file about grid, soil and vegetation.                                                |
| filename%lai                     | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Default leaf area index file name.                                                                      |
| filename%soilcolor               | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Soil color file name. E.g. "CABLE-AUX/offline/soilcolor_global_1x1.nc"             |
| filename%gw_elev                 | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Elevation data filename. This data is used by the groundwater option only.                              |
| filename%fxpft                   | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Plant functional type fraction, wood harvest and secondary harvest file.                                |
| filename%fxluh2cable             | character(len=500) | any string of max. 500 characters                          | uninitialised                    | 12 land-use states into 17 CABLE plant functional types mapping file name.                              |
| filename%gridnew                 | character(len=500) | any string of max. 500 characters                          | uninitialised                    | Updated gridinfo file name.                                                                             |
| filename%trunk_sumbal            | character(len=500) | any string of max. 500 characters                          | '.trunk_sumbal'                  | Input filename to read combined energy and water balance at each timestep (control run). Used when `consistency_check` is TRUE |
| filename%new_sumbal              | character(len=500) | any string of max. 500 characters                          | 'new_sumbal'                     | Output filename to write combined energy and water balance at each timestep (current run). Used when `consistency_check` is TRUE |
| vegparmnew                       | logical            | .TRUE. .FALSE.                                             | .FALSE. but was .TRUE. in ESM1.5 | Use new format for vegetation parameter files when .TRUE.                                               |
| soilparmnew                      | logical            | .TRUE. .FALSE.                                             | uninitialised                    | Use new format for soil parameter files when .TRUE.                                                     |
| spinup                           | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Spin up the model when .TRUE.                                                                           |
| delsoilM                         | real               | any real number                                            | uninitialised                    | Allowed variation in soil moisture for spin up.                                                         |
| delsoilT                         | real               | any real number                                            | uninitialised                    | Allowed variation in soil temperature for spin up.                                                      |
| output%restart                   | logical            | .TRUE. .FALSE.                                             | uninitialised                    | Create a restart file when .TRUE..                                                                      |
| output%met                       | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output input meteorological data when .TRUE..                                                           |
| output%flux                      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output convective, runoff, NEE when TRUE                                                                |
| output%soil                      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output soil states when .TRUE..                                                                         |
| output%snow                      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output snow states when .TRUE..                                                                         |
| output%radiation                 | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output Net radiation and albedo when .TRUE..                                                            |
| output%carbon                    | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output NEE, GPP, NPP, stores when .TRUE..                                                               |
| output%veg                       | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output vegetation states when .TRUE..                                                                   |
| output%params                    | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output the input parameters used to produce run then .TRUE..                                            |
| output%balances                  | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Output energy and water balances when .TRUE..                                                           |
| output%grid                      | character(len=7)   | 'default' 'land' 'mask' 'ALMA'                             | 'default'                        | Output grid convention. `'land'` outputs land-only points, `'mask'` outputs masked spatial grids. `'ALMA'` uses the Assistance for Land-surface Modelling Activities convention. `'default'` uses the whichever convention the input meteorological file is using. |
| output%averaging                 | character(len=7)   | 'all' 'daily' 'monthly' 'user6'                            | 'all'                            | Output averaging.                                                                                       |
| check%ranges                     | integer            |  0 (`NO_CHECK`), 1 (`ON_TIMESTEP`), 2 (`ON_WRITE`)         | uninitialised                    | Check input and output variables at certain timesteps against valid ranges.                             |
| check%exit                       | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Behaviour on failed range checks. If true , write to console and exit the program, else write to logfile as warning |
| check%energy_bal                 | logical            | .TRUE. .FALSE.                                             | uninitialised                    | Check the energy balance.                                                                               |
| check%mass_bal                   | logical            | .TRUE. .FALSE.                                             | uninttialised                    | Check the water/mass balance.                                                                           |
| verbose                          | logical            | .TRUE. .FALSE.                                             | uninitialised                    | Write details of every grid cell initialisation and parameters to log.                                  |
| leaps                            | logical            | .TRUE. .FALSE.                                             | depends on configuration         | Calculate timing with leap years.                                                                       |
| logn                             | integer            | any integer number                                         | uninitialised                    | Log file number - declared in input module.                                                             |
| fixedCO2                         | real               | any real number                                            | 350                              | CO$_2$ concentration if not found in met file, in units ppmv.                                           |
| spincasa                         | logical            | .TRUE. .FALSE.                                             | uninitialised                    | Spin CASA before running the model if TRUE, and should be set to FALSE if spincasainput is .TRUE..      |
| l_casacnp                        | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use CASA-CNP with CABLE.                                                                                |
| l_landuse                        | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use land-use change.                                                                                    |
| l_laiFeedbk                      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use prognostic leaf area index.                                                                         |
| l_vcmaxFeedbk                    | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use prognostic Vcmax.                                                                                   |
| CASAonly                         | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | ONLY run CASA-CNP.                                                                                      |
| icycle                           | integer            | 0, 1, 2, 3, >10                                            | 0                                | 0 for not using CASA-CNP, 1 for C, 2 for C+N, 3 for C+N+P, >10 excludes call to CBM (main call to cable).|
| casafile%cnpipool                | character(len=99)  | any string of max. 99 characters                           | ''                               | Initial CNP pool size file name.                                                                        |
| casafile%cnpbiome                | character(len=99)  | any string of max. 99 characters                           | uninitialised                    | Biome specific biogeochemical parameters file name.                                                     |
| casafile%cnpepool                | character(len=99)  | any string of max. 99 characters                           | uninitialised                    | End of run pool size file name.                                                                         |
| casafile%cnpmetout               | character(len=99)  | any string of max. 99 characters                           | uninitialised                    | Output daily met forcing for spinning CASA-CNP file name.                                               |
| casafile%cnpmetin                | character(len=99)  | any string of max. 99 characters                           | uninitialised                    | List of daily met files for spinning CASA-CNP file name.                                                |
| casafile%phen                    | character(len=99)  | any string of max. 99 characters                           | uninitialised                    | Modis phenology file name.                                                                              |
| casafile%cnpflux                 | character(len=99)  | any string of max. 99 characters                           | uninitialised                    | Output file name for CNP fluxes.                                                                        |
| casafile%c2cdumppath             | character(len=99)  | any string of max. 99 characters                           | ''                               | CABLE to CASA dump for CASA spin up.                                                                    |
| ncciy                            | integer            | any integer                                                | uninitialised                    | GSWP year. 0 for not using GSWP; 4-digit year input for year of GSWP meteorological forcing.            |
| gswpfile%rainf                   | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for rain meteorological forcing.                                                             |
| gswpfile%snowf                   | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for snow meteorological forcing.                                                             |
| gswpfile%LWdown                  | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for long wave radiation meteorological forcing.                                              |
| gswpfile%SWdown                  | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for short wave radiation meteorological forcing.                                             |
| gswpfile%PSurf                   | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for surface pressure meteorological forcing.                                                 |
| gswpfile%Qair                    | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for humidity meteorological forcing.                                                         |
| gswpfile%Tair                    | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for temperature meteorological forcing.                                                      |
| gswpfile%wind                    | character(len=200) | any string of max. 200 characters                          | uninitialised                    | Input file for wind meteorological forcing.                                                             |
| redistrb                         | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Turn on/off the hydraulic redistribution.                                                               |
| wiltParam                        | real               | any real                                                   | 0.0                              | Wilt parameter in hydraulic redistribution module.                                                      |
| satuParam                        | real               | any real                                                   | 0.0                              | Saturation parameter in hydraulic redistribution module.                                                |
| cable_user%FWSOIL_SWITCH         | character          | 'standard' 'non-linear extrapolation' 'Lai and Ktaul 2000' | ''                               | Controls root water uptake function.                                                                    |
| cable_user%DIAG_SOIL_RESP        | character(len=3)   | 'ON' 'OFF'                                                 | ''                               | Controls soil respiration scheme when CASA-CNP not used.                                                |
| cable_user%LEAF_RESPIRATION      | character(len=3)   | 'ON' 'OFF'                                                 | uninitialised                    | Controls what is output into the gross primary productivity stash variable.                             |
| cable_user%RUN_DIAG_LEVEL        | character(len=5)   | 'BASIC' 'NONE'                                             | uninitialised                    | Controls output from CABLE to standard out.                                                             |
| cable_user%CONSISTENCY_CHECK     | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | TRUE outputs combined fluxes at each timestep for comparison to a control run.                          |
| cable_user%CASA_DUMP_READ        | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | TRUE reads CASA forcing from netcdf format.                                                             |
| cable_user%CASA_DUMP_WRITE       | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | TRUE outputs CASA forcing in netcdf format.                                                             |
| cable_user%SSNOW_POTEV           | character(len=3)   | 'HDM' 'P-M'                                                | ''                               | Use the “Humidity Deficit Method” or “Penman-Monteith” method.                                          |
| cable_user%GS_SWITCH             | character(len=20)  | 'medlyn' 'leuning'                                         | 'medlyn'                         | Stomatal conductance model.                                                                             |
| cable_user%l_revised_coupling    | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Revised coupling on implicit step of ACCESS-CM2.                                                        |
| cable_user%l_rev_corr            | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Apply revised sensitvity/correction terms to soilsnow energy balance.                                   |
| cable_user%soil_thermal_fix      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use alternative soil conductivity implementation.                                                       |
| cable_user%phenology_switch      | character(len=20)  | 'MODIS' 'climate'                                          | 'MODIS'                          | Use prescribed MODIS phenology or climate dependant phenology.                                          |
| cable_user%RunIden               | character(len=10)  | any string of max. 10 characters                           | 'STANDARD'                       | Run identifier string for input/output files.                                                           |
| cable_user%MetType               | character(len=6)   | '' 'gswp' 'gswp3' 'plum' 'cru' 'site' 'bios'               | ' '                              | Type of input meteorological data.                                                                      |
| cable_user%soil_struc            | character(len=20)  | 'default' 'sli'                                            | 'default'                        | Use default or soil-litter-iso soil model.                                                              |
| cable_user%POP_out               | character(len=3)   | 'epi' 'rst' 'ini'                                          | 'rst'                            | POP restart file type. `'epi'` is end of year state `'rst'` is a standard restart file, `'ini'` is an initialisation restart file. |
| cable_user%POP_rst               | character(len=50)  | any string of max. 50 characters                           | ' '                              | POP restart file directory.                                                                             |
| cable_user%casa_out_freq         | character(len=8)   | 'daily' 'monthly' 'annually'                               | 'annually'                       | CASA output frequency.                                                                                  |
| cable_user%vcmax                 | character(len=10)  | 'standard' 'Walker2014'                                    | 'standard'                       | Alternative functional form of $V_{cmax}$ (maximum RuBP carboxylation rate). For `'Walker2014'` see [Walker et al., 2014][Walker] |
| cable_user%POPLUC_RunType        | character(len=10)  | 'static' 'init' 'restart'                                  | 'static'                         | Type of POP-LUC run. `'init'` is initial POP-LUC run, `'restart'` is a restart run, `'static'` disables dynamic vegetation and sets `cable_user%POPLUC = .FALSE.`.|
| cable_user%call_pop              | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Enable POP.                                                                                             |
| cable_user%POP_fromZero          | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Run POP initialisation.                                                                                 |
| cable_user%call_climate          | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Track climate variables for use in phenology and potential plant functional type modules.               |
| cable_user%climate_fromzero      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Initialise climate variables with 0.                                                                    |
| cable_user%CASA_fromzero         | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Initialise CASA from 0 without restart files.                                                           |
| cable_user%POPLUC                | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Enable POP-LUC.                                                                                         |
| cable_user%casa_spin_startyear   | integer            | any integer                                                | 1950                             | CASA spinup start year.                                                                                 |
| cable_user%casa_spin_endyear     | integer            | any integer                                                | 1960                             | CASA spinup end year.                                                                                   |
| cable_user%yearstart             | integer            | any integer                                                | 0                                | Start year of run.                                                                                      |
| cable_user%yearend               | integer            | any integer                                                | 0                                | End year of run.                                                                                        |
| cable_user%casa_nrep             | integer            | any integer                                                | 1                                | Number of repetitions of CASA spinup.                                                                   |
| cable_user%smrf_name             | character(len=10)  | 'CASA-CNP' 'SOILN' 'TRIFFID'                               | uninitialised                    | Soil moisture respiration function.                                                                     |
| cable_user%strf_name             | character(len=10)  | 'CASA-CNP' 'SOILN' 'TRIFFID'                               | uninitialised                    | Soil temperature respiration function.                                                                  |
| cable_user%LogWorker             | logical            | .TRUE. .FALSE.                                             | .TRUE.                           | Write output of each worker.                                                                            |
| cable_user%l_new_roughness_soil  | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use alternative roughness and flux resistance scheme.                                                   |
| cable_user%l_new_runoff_speed    | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use alternative soil water drainage scheme.                                                             |
| cable_user%l_new_reduce_soilevp  | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use alternative soil evaporation scheme.                                                                |
| cable_user%srf                   | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Switch for customized soil respiration.                                                                 |
| cable_user%litter                | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Add litter layer on top of the soil surface.                                                            |
| cable_user%gw_model              | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Enable alternate soil model that accounts for groundwater.                                              |
| cable_user%gswp3                 | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use GSWP3 forcing.                                                                                      |
| cable_user%or_evap               | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Enable Dani Or's evaporation scheme.                                                                    |
| cable_user%sync_nc_file          | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Write netCDF output to file without closing during run.                                                 |
| cable_user%max_spins             | integer            | any integer                                                | 1                                | Maximum number of iterations for a spinup. CABLE will act as if convergence has been reached after max_spins iterations even if it was not reached. |
| cable_user%fix_access_roots      | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Use plant functional type dependent roots in ACCESS                                                     |
| cable_user%access13roots         | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Switch to use ACCESS1.3 %froot.                                                                         |
| cable_user%l_limit_labile        | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Limit labile in spinup.                                                                                 |
| cable_user%NtilesThruMetFile     | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | Specify Ntiles through met file.                                                                        |
| cable_user%l_ice_consistency     | logical            | .TRUE. .FALSE.                                             | .FALSE.                          | If true, ensures consistency between soil and vegetation tiles with permanent ice. All tiles with permanent ice for soil will have ice for vegetation and vice-versa. All the parameters for these new ice tiles are updated to the ice parameters.                       |

## For offline applications ##

```fortran
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
   check%ranges     = 1  ! Range-checks on every timestep
   check%exit = .TRUE.  ! Exit the program if range checks fail
   check%energy_bal = .TRUE.  ! energy balance
   check%mass_bal   = .TRUE.  ! water/mass balance
   verbose = .TRUE. ! write details of every grid cell init and params to log?
   leaps = .TRUE. ! calculate timing with leap years
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

## For ACCESS applications ##


```fortran
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

<!-- markdown-link-check-disable-line --> [Walker]: https://doi.org/10.1002/ece3.1173
[obsolete]: ../other_resources/obsolete_and_deprecated_features/obsolete_and_deprecated_features.md