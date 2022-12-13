
# Namelist file cable.nml #

The namelist file for CABLE is cable.nml. For offline applications it
should be located in the directory in which you are running CABLE. For
ACCESS cases, cable.nml should be located in ~/CABLE-AUX/UM on the
machine that CABLE is being executed on.

The cable.nml file includes some settings that are common across all CABLE
applications and others that are required only for offline
applications. The following are annotated examples of cable.nml:

## For offline applications ##

``` fortran-free-form
&cable
   filename%met = '/projects/access/CABLE-AUX/offline/TumbaFluxnet.1.3_met.nc'
   filename%out = 'out_cable.nc'
   filename%log = 'log_cable.txt'
   filename%restart_in  = ' '
   filename%restart_out = './restart_out.nc'
   filename%type    = '/projects/access/CABLE-AUX/offline/gridinfo_CSIRO_1x1.nc'
   filename%veg    = '/projects/access/CABLE-AUX/core/biogeophys/veg_params_vbeta1.txt'
   filename%soil    = '/projects/access/CABLE-AUX/core/biogeophys/def_soil_params.txt'
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
   cable_user%FWSOIL_SWITCH = 'standard'        ! choices are:
                                                 ! 1. standard
                                                 ! 2. non-linear extrapolation
                                                 ! 3. Lai and Ktaul 2000
   cable_user%DIAG_SOIL_RESP = 'ON '
   cable_user%LEAF_RESPIRATION = 'ON '
   cable_user%RUN_DIAG_LEVEL= 'BASIC'        ! choices are:
                                                 ! 1. BASIC
                                                 ! 1. NONE
   cable_user%CONSISTENCY_CHECK= .TRUE.      ! TRUE outputs combined fluxes at each timestep for comparison to a control run
   cable_user%CASA_DUMP_READ = .FALSE.      ! TRUE reads CASA forcing from netcdf format
   cable_user%CASA_DUMP_WRITE = .FALSE.      ! TRUE outputs CASA forcing in netcdf format
   cable_user%SSNOW_POTEV= 'HDM'      ! Humidity Deficit Method
&end
```


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
