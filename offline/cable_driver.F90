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
! Purpose: Offline driver for CABLE
! Contact: Bernard.Pak@csiro.au
!
! History: Since 1.4b, capability to run global offline (ncciy = YEAR),
!          inclusion of call to CASA-CNP (icycle>0)
!          exclusion of call to cbm (icycle>10)
!          soil_snow_type now ssnow (instead of ssoil)
!
!
! ========================================untitled======================================
! Uses:           cable_def_types_mod
!                 cable_IO_vars_module
!                 cable_common_module
!                 cable_input_module
!                 cable_output_module
!                 cable_cbm_module
!                 casadimension
!                 casavariable
!
! CALLs:       open_met_file
!              load_parameters
!              open_output_file
!              get_met_data
!              casa_feedback
!              cbm
!              bgcdriver
!              sumcflux
!              write_output
!              casa_poolout
!              casa_fluxout
!              create_restart
!              close_met_file
!              close_output_file
!              prepareFiles
!
!
! input  file: [SiteName].nc
!              poolcnpIn[SiteName].csv -- for CASA-CNP only
!              gridinfo_CSIRO_1x1.nc
!              def_veg_params.txt
!              def_soil_params.txt -- nearly redundant, can be switched on
!              restart_in.nc -- not strictly required
!
! output file: log_cable.txt
!              out_cable.nc
!              restart_out.nc
!              poolcnpOut.csv -- from CASA-CNP
!==============================================================================

PROGRAM cable_offline_driver

  USE cable_def_types_mod
  USE cable_io_vars_module, ONLY: logn, gswpfile, ncciy, leaps, &
       verbose, fixedCO2, output, check, patchout, soilparmnew, &
       timeunits, exists, calendar, landpt
  USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user, &
       cable_runtime, filename, &
       redistrb, wiltParam, satuParam, CurYear, &
       IS_LEAPYEAR, IS_CASA_TIME, calcsoilalbedo, get_unit, &
       report_version_no, kwidth_gl
  use cable_data_module,    only: driver_type, point2constants
  use cable_input_module,   only: open_met_file, load_parameters, get_met_data, close_met_file, &
       ncid_rain, ncid_snow, ncid_lw, ncid_sw, ncid_ps, ncid_qa, ncid_ta, ncid_wd
  use cable_output_module,  only: create_restart, open_output_file, write_output, close_output_file
  use cable_write_module,   only: nullify_write
  USE cable_cbm_module
  USE cable_diag_module
  !mpidiff
  USE cable_climate_mod

  ! modules related to CASA-CNP
  USE casadimension,        ONLY: icycle
  USE casavariable,         ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
       !mpidiff
       casa_met, casa_balance, zero_sum_casa, update_sum_casa
  USE phenvariable,         ONLY: phen_variable

  !! vh_js !!
  ! modules related to POP
  USE POP_Types,     ONLY: POP_TYPE
  USE POPLUC_Types,  ONLY: POPLUC_Type
  USE POPLUC_Module, ONLY: WRITE_LUC_OUTPUT_NC, WRITE_LUC_OUTPUT_GRID_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, POPLUC_set_patchfrac
  USE POP_Constants, ONLY: HEIGHT_BINS, NCOHORT_MAX, shootfrac

  ! Fire Model BLAZE
  USE BLAZE_MOD,     ONLY: TYPE_BLAZE, INI_BLAZE, BLAZE_ACCOUNTING,  WRITE_BLAZE_OUTPUT_NC
  USE SIMFIRE_MOD,   ONLY: TYPE_SIMFIRE, INI_SIMFIRE

  ! 13C
  use cable_c13o2_def,         only: c13o2_delta_atm, c13o2_flux, c13o2_pool, c13o2_luc, &
       c13o2_update_sum_pools, c13o2_zero_sum_pools
  use cable_c13o2,             only: c13o2_save_luc, c13o2_update_luc, &
       c13o2_write_restart_flux, c13o2_write_restart_pools, c13o2_write_restart_luc, &
       c13o2_create_output, c13o2_write_output, c13o2_close_output, c13o2_nvars_output
  use cable_c13o2,             only: c13o2_print_delta_flux, c13o2_print_delta_pools, c13o2_print_delta_luc
  use mo_isotope,              only: isoratio ! vpdbc13
  use mo_c13o2_photosynthesis, only: c13o2_discrimination_simple, c13o2_discrimination
  use cable_data_module,       only: icanopy_type

  ! PLUME-MIP only
  USE CABLE_PLUME_MIP,      ONLY: PLUME_MIP_TYPE, PLUME_MIP_GET_MET,&
       PLUME_MIP_INIT

  USE CABLE_CRU,            ONLY: CRU_TYPE, CRU_GET_SUBDIURNAL_MET, CRU_INIT, cru_close
  USE CABLE_site,           ONLY: site_TYPE, site_INIT, site_GET_CO2_Ndep

  ! BIOS only
  USE cable_bios_met_obs_params,   ONLY:  cable_bios_read_met, cable_bios_init, &
                                          cable_bios_load_params, &
                                          cable_bios_load_climate_params

  ! LUC_EXPT only
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE, LUC_EXPT_INIT, close_luh2

#ifdef __NAG__
  USE F90_UNIX
#endif

  IMPLICIT NONE

  ! CABLE namelist: model configuration, runtime/user switches
  CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml'

  ! timing variables
  INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
  INTEGER, PARAMETER ::  mloop  = 4 !MCTEST 30  ! CASA-CNP PreSpinup loops
  INTEGER :: LALLOC ! allocation coefficient for passing to spincasa

  INTEGER        ::                                                           &
       ktau,       &  ! increment equates to timestep, resets if spinning up
       ktau_tot,   &  ! NO reset when spinning up, total timesteps by model
       kend,       &  ! no. of time steps in run
       !CLN kstart = 1, &  ! timestep to start at
       koffset = 0, &  ! timestep to start at
       koffset_met = 0, &  !offfset for site met data ('site' only)
       ktauday,    &  ! day counter for CASA-CNP
       idoy,       &  ! day of year (1:365) counter for CASA-CNP
       nyear,      &  ! year counter for CASA-CNP
       casa_it,    &  ! number of calls to CASA-CNP
       YYYY,       &  !
       RYEAR,      &  !
       RRRR,       &  !
       NRRRR,      &  !
       ctime,      &  ! day count for casacnp
       LOY, &         ! days in year
       count_sum_casa, & ! number of time steps over which casa pools &
  !and fluxes are aggregated (for output)
       wlogn = 10001

  REAL :: dels                        ! time step size in seconds

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: GSWP_MID
  CHARACTER :: dum*9, str1*9, str2*9, str3*9

  ! CABLE variables
  TYPE(met_type)       :: met     ! met input variables
  TYPE(air_type)       :: air     ! air property variables
  TYPE(canopy_type)    :: canopy  ! vegetation variables
  TYPE(radiation_type) :: rad     ! radiation variables
  TYPE(roughness_type) :: rough   ! roughness varibles
  TYPE(balances_type)  :: bal     ! energy and water balance variables
  TYPE(soil_snow_type) :: ssnow   ! soil and snow variables
  !mpidiff
  TYPE(climate_type)   :: climate     ! climate variables

  ! CABLE parameters
  ! JK: C is used for both driver_type and i_canopy_type constants throughout the code,
  !     should be avoided ?
  TYPE(soil_parameter_type) :: soil ! soil parameters
  TYPE(veg_parameter_type)  :: veg  ! vegetation parameters
  TYPE(driver_type)    :: C         ! constants used locally
  TYPE(icanopy_type)   :: PHOTO     ! photosynthesis constants

  TYPE(sum_flux_type)  :: sum_flux ! cumulative flux variables
  TYPE(bgc_pool_type)  :: bgc  ! carbon pool variables

  ! CASA-CNP variables
  TYPE(casa_biome)     :: casabiome
  TYPE(casa_pool)      :: casapool
  TYPE(casa_flux)      :: casaflux
  TYPE(casa_pool)      :: sum_casapool
  TYPE(casa_flux)      :: sum_casaflux
  TYPE(casa_met)       :: casamet
  TYPE(casa_balance)   :: casabal
  TYPE(phen_variable)  :: phen
  !! vh_js !!
  TYPE(POP_TYPE)        :: POP
  TYPE(POPLUC_TYPE)     :: POPLUC
  TYPE(PLUME_MIP_TYPE)  :: PLUME
  TYPE(CRU_TYPE)        :: CRU
  TYPE(site_TYPE)       :: site
  TYPE(LUC_EXPT_TYPE)   :: LUC_EXPT
  CHARACTER             :: cyear*4
  CHARACTER             :: ncfile*99

  ! BLAZE variables
  TYPE(TYPE_BLAZE)    :: BLAZE
  TYPE(TYPE_SIMFIRE)  :: SIMFIRE

  ! 13C
  type(c13o2_flux)  :: c13o2flux
  type(c13o2_pool)  :: c13o2pools, sum_c13o2pools
  type(c13o2_luc)   :: c13o2luc
  real(r_2), dimension(:,:), allocatable :: casasave
  real(r_2), dimension(:,:), allocatable :: lucsave
  ! I/O
  integer :: c13o2_outfile_id
  character(len=40), dimension(c13o2_nvars_output) :: c13o2_vars
  integer,           dimension(c13o2_nvars_output) :: c13o2_var_ids
  ! discrimination
  integer :: ileaf
  real(r_2), dimension(:,:), allocatable :: gpp, diff
  real(r_2), dimension(:),   allocatable :: Ra 
  ! delta-13C of atmospheric CO2
  integer            :: iunit, ios
  real               :: iyear
  integer            :: c13o2_atm_syear, c13o2_atm_eyear
  character(len=100) :: header

  ! declare vars for switches (default .FALSE.) etc declared thru namelist
  LOGICAL, SAVE           :: &
       vegparmnew = .FALSE.,       & ! using new format input file (BP dec 2007)
       spinup = .FALSE.,           & ! model spinup to soil state equilibrium?
       spinConv = .FALSE.,         & ! has spinup converged?
       spincasainput = .FALSE.,    & ! TRUE: SAVE input req'd to spin CASA-CNP;
                                ! FALSE: READ input to spin CASA-CNP
       spincasa = .FALSE.,         & ! TRUE: CASA-CNP Will spin mloop times,
                                ! FALSE: no spin up
       l_casacnp = .FALSE.,        & ! using CASA-CNP with CABLE
       l_laiFeedbk = .FALSE.,      & ! using prognostic LAI
       l_vcmaxFeedbk = .FALSE.,    & ! using prognostic Vcmax
       CASAONLY      = .FALSE.,    & ! ONLY Run CASA-CNP
       CALL1 = .TRUE.,             &
       SPINon= .TRUE.

  REAL              :: &
       delsoilM,         & ! allowed variation in soil moisture for spin up
       delsoilT            ! allowed variation in soil temperature for spin up

 INTEGER :: Metyear, Y, LOYtmp

  ! temporary storage for soil moisture/temp. in spin up mode
  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
       soilMtemp, &
       soilTtemp

  ! timing
  REAL:: etime ! Declare the type of etime(), For receiving user and system time, total time

  !___ unique unit/file identifiers for cable_diag: arbitrarily 5 here
  INTEGER, SAVE :: iDiagZero=0, iDiag1=0, iDiag2=0, iDiag3=0, iDiag4=0

  ! switches etc defined thru namelist (by default cable.nml)
  NAMELIST/CABLE/ &
       filename,         & ! TYPE, containing input filenames
       vegparmnew,       & ! use new soil param. method
       soilparmnew,      & ! use new soil param. method
       calcsoilalbedo,   & ! albedo considers soil color Ticket #27
       spinup,           & ! spinup model (soil) to steady state
       delsoilM,delsoilT,& !
       output,           &
       patchout,         &
       check,            &
       verbose,          &
       leaps,            &
       logn,             &
       fixedCO2,         &
       spincasainput,    &
       spincasa,         &
       l_casacnp,        &
       l_laiFeedbk,      &
       l_vcmaxFeedbk,    &
       icycle,           &
       casafile,         &
       ncciy,            &
       gswpfile,         &
       redistrb,         &
       wiltParam,        &
       satuParam,        &
       cable_user           ! additional USER switches
  
  !mpidiff
  INTEGER :: i,x,kk

  ! Vars for standard for quasi-bitwise reproducability b/n runs
  ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
  CHARACTER(len=30), PARAMETER :: &
       Ftrunk_sumbal  = ".trunk_sumbal", &
       Fnew_sumbal    = "new_sumbal"

  REAL(r_2) :: &
       trunk_sumbal = 0.0, & !
       new_sumbal = 0.0,   &
       new_sumfpn = 0.0,   &
       new_sumfe = 0.0

  INTEGER :: nkend=0
  INTEGER :: ioerror
  INTEGER :: count_bal = 0
  ! END header

  ! Open, read and close the namelist file.
  OPEN( 10, FILE = CABLE_NAMELIST )
  READ( 10, NML=CABLE )   !where NML=CABLE defined above
  CLOSE(10)

  ! Open, read and close the consistency check file.
  ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
  IF (cable_user%consistency_check) THEN
     OPEN( 11, FILE = Ftrunk_sumbal,STATUS='old',ACTION='READ',IOSTAT=ioerror )
     IF (ioerror==0) THEN
        READ( 11, * ) trunk_sumbal  ! written by previous trunk version
     ENDIF
     CLOSE(11)
  ENDIF

  ! Open log file:
  OPEN(logn,FILE=filename%log)

  CALL report_version_no( logn )

  IF( IARGC() > 0 ) THEN
     CALL GETARG(1, filename%met)
     CALL GETARG(2, casafile%cnpipool)
  ENDIF

  ! INITIALISATION depending on nml settings
  IF (TRIM(cable_user%MetType) .EQ. 'gswp') THEN
     IF ( CABLE_USER%YearStart.eq.0 .and. ncciy.gt.0) THEN
        CABLE_USER%YearStart = ncciy
        CABLE_USER%YearEnd = ncciy
     ELSEIF  ( CABLE_USER%YearStart.eq.0 .and. ncciy.eq.0) THEN
        write(*,*) 'undefined start year for gswp met: '
        write(*,*) 'enter value for ncciy or'
        write(*,*) '(CABLE_USER%YearStart and  CABLE_USER%YearEnd) in cable.nml'

        write(logn,*) 'undefined start year for gswp met: '
        write(logn,*) 'enter value for ncciy or'
        write(logn,*) '(CABLE_USER%YearStart and  CABLE_USER%YearEnd) in cable.nml'

        stop
     ENDIF
  ENDIF

  CurYear = CABLE_USER%YearStart

  IF ( icycle .GE. 11 ) THEN
     icycle                     = icycle - 10
     CASAONLY                   = .TRUE.
     CABLE_USER%CASA_DUMP_READ  = .TRUE.
     CABLE_USER%CASA_DUMP_WRITE = .FALSE.
  ELSEIF ( icycle .EQ. 0 ) THEN
     CABLE_USER%CASA_DUMP_READ  = .FALSE.
     spincasa                   = .FALSE.
     !! vh_js !!
     CABLE_USER%CALL_POP        = .FALSE.
     CABLE_USER%CALL_BLAZE      = .FALSE.
  ENDIF

  !! vh_js !!
  IF (icycle.gt.0) THEN
     l_casacnp = .TRUE.
  ELSE
     l_casacnp = .FALSE.
  ENDIF

  !! vh_js !! suggest LALLOC should ulitmately be a switch in the .nml file
  IF (CABLE_USER%CALL_POP) THEN
     LALLOC = 3 ! for use with POP: makes use of pipe model to partition between stem and leaf
  ELSE
     LALLOC = 0 ! default
  ENDIF

  ! IF ( .NOT. spinup ) THEN
  !    IF ( spincasa ) THEN
  !       spincasa = .FALSE.
  !       WRITE(*,*)    "spinup == .FALSE. -> spincasa set to .F."
  !       WRITE(logn,*) "spinup == .FALSE. -> spincasa set to .F."
  !    ENDIF
  ! ENDIF

  IF ( TRIM(cable_user%MetType) .EQ. 'gpgs' ) THEN
     leaps = .TRUE.
     calendar = "standard"
     cable_user%MetType = 'gswp'
  ENDIF

  cable_runtime%offline = .TRUE.

  ! associate pointers used locally with global definitions
  CALL point2constants(C)

  IF( l_casacnp  .AND. ( icycle == 0 .OR. icycle > 3 ) )                   &
       STOP 'icycle must be 1 to 3 when using casaCNP'
  !IF( ( l_laiFeedbk .OR. l_vcmaxFeedbk ) )       &
  !   STOP 'casaCNP required to get prognostic LAI or Vcmax'
  IF( l_vcmaxFeedbk .AND. icycle < 1 )                                     &
       STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
  IF( icycle > 0 .AND. ( .NOT. soilparmnew ) )                             &
       STOP 'casaCNP must use new soil parameters'

  NRRRR = MERGE(MAX(CABLE_USER%CASA_NREP,1), 1, CASAONLY)
  ! casa time count
  ctime = 0

  ! INISTUFF

  ! Open met data and get site information from netcdf file. (NON-GSWP ONLY!)
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites.

  ! IF ( TRIM(cable_user%MetType) .NE. "gswp" .AND. &
  !      TRIM(cable_user%MetType) .NE. "gpgs" .AND. &
  !      TRIM(cable_user%MetType) .NE. "plum" .AND. &
  !      TRIM(cable_user%MetType) .NE. "cru") THEN
  IF (TRIM(cable_user%MetType) .EQ. 'site' .OR. &
       TRIM(cable_user%MetType) .EQ. '') THEN
     CALL open_met_file( dels, koffset, kend, spinup, C%TFRZ )
     IF ( koffset .NE. 0 .AND. CABLE_USER%CALL_POP ) THEN
        WRITE(*,*) "When using POP, episode must start at Jan 1st!"
        STOP 991
     ENDIF
  ELSE IF ( NRRRR .GT. 1 ) THEN
     IF(.NOT.ALLOCATED(GSWP_MID)) ALLOCATE( GSWP_MID( 8, CABLE_USER%YearStart:CABLE_USER%YearEnd ) )
  ENDIF

  ! 13C
  ! Read atmospheric delta-13C values
  if (cable_user%c13o2) then
     ! get start and end year
     call get_unit(iunit)
     open(iunit, file=trim(cable_user%c13o2_delta_atm_file), status="old", action="read")
     read(iunit, fmt=*, iostat=ios) header
     read(iunit, fmt=*, iostat=ios) iyear
     c13o2_atm_syear = floor(iyear)
     do while (ios==0)
        read(iunit, fmt=*, iostat=ios) iyear
     end do
     c13o2_atm_eyear = floor(iyear)
     ! allocate d13ca(start:end)
     allocate(c13o2_delta_atm(c13o2_atm_syear:c13o2_atm_eyear))
     ! read annual global mean values
     rewind(iunit)
     read(iunit, fmt=*, iostat=ios) header
     do while (ios==0)
        read(iunit, fmt=*, iostat=ios) iyear, c13o2_delta_atm(floor(iyear))
     end do
     close(iunit)
     c13o2_delta_atm = c13o2_delta_atm/1000._r_2
  endif

  ! outer loop - spinup loop no. ktau_tot :
  RYEAR = 0
  ktau_tot = 0
  SPINon   = .TRUE.
  ! YearStart = CABLE_USER%YearStart
  ! YearEnd = CABLE_USER%YearEnd
  ! cable_user%CASA_SPIN_ENDYEAR

  SPINLOOP:DO WHILE ( SPINon )

     NREP: DO RRRR = 1, NRRRR
        IF (TRIM(cable_user%MetType) .EQ. "bios") THEN
           CALL CPU_TIME(etime)

           CALL cable_bios_init(dels,curyear,met,kend,ktauday)

           koffset   = 0
           leaps = .true.

           write(str1,'(i4)') curyear
           str1 = adjustl(str1)
           write(str2,'(i02)') 1
           str2 = adjustl(str2)
           write(str3,'(i02)') 1
           str3 = adjustl(str3)
           timeunits="seconds since "//trim(str1)//"-"//trim(str2)//"-"//trim(str3)//" 00:00:00"
           calendar = 'standard'
        ENDIF

        print*, "CABLE_USER%YearStart,  CABLE_USER%YearEnd", CABLE_USER%YearStart,  CABLE_USER%YearEnd

        YEAR: DO YYYY= CABLE_USER%YearStart,  CABLE_USER%YearEnd
           CurYear = YYYY
           IF ( leaps .AND. IS_LEAPYEAR( YYYY ) ) THEN
              LOY = 366
           ELSE
              LOY = 365
           ENDIF
           ! Check for gswp run
           IF ( TRIM(cable_user%MetType) .EQ. 'gswp' ) THEN
              ncciy = CurYear

              CALL prepareFiles(ncciy)
              IF ( RRRR .EQ. 1 ) THEN
                 CALL open_met_file( dels, koffset, kend, spinup, C%TFRZ )
                 IF (leaps.and.is_leapyear(YYYY).and.kend.eq.2920) THEN
                    STOP 'LEAP YEAR INCOMPATIBILITY WITH INPUT MET !!!'
                 ENDIF
                 IF ( NRRRR .GT. 1 ) THEN
                    GSWP_MID(1,YYYY) = ncid_rain
                    GSWP_MID(2,YYYY) = ncid_snow
                    GSWP_MID(3,YYYY) = ncid_lw
                    GSWP_MID(4,YYYY) = ncid_sw
                    GSWP_MID(5,YYYY) = ncid_ps
                    GSWP_MID(6,YYYY) = ncid_qa
                    GSWP_MID(7,YYYY) = ncid_ta
                    GSWP_MID(8,YYYY) = ncid_wd
                 ENDIF
              ELSE
                 ncid_rain = GSWP_MID(1,YYYY)
                 ncid_snow = GSWP_MID(2,YYYY)
                 ncid_lw   = GSWP_MID(3,YYYY)
                 ncid_sw   = GSWP_MID(4,YYYY)
                 ncid_ps   = GSWP_MID(5,YYYY)
                 ncid_qa   = GSWP_MID(6,YYYY)
                 ncid_ta   = GSWP_MID(7,YYYY)
                 ncid_wd   = GSWP_MID(8,YYYY)
                 kend      = ktauday * LOY
              ENDIF
           ELSE IF ( TRIM(cable_user%MetType) .EQ. 'plum' ) THEN
              ! PLUME experiment setup using WATCH
              IF ( CALL1 ) THEN
                 CALL CPU_TIME(etime)
                 CALL PLUME_MIP_INIT( PLUME )

                 dels      = PLUME%dt
                 koffset   = 0
                 leaps = PLUME%LeapYears

                 write(str1,'(i4)') CurYear
                 str1 = adjustl(str1)
                 write(str2,'(i02)') 1
                 str2 = adjustl(str2)
                 write(str3,'(i02)') 1
                 str3 = adjustl(str3)
                 timeunits="seconds since "//trim(str1)//"-"//trim(str2)//"-"//trim(str3)//" 00:00"
                 if (leaps) then
                    calendar = "standard"
                 else
                    calendar = "noleap"
                 endif
              ENDIF
              IF ( .NOT. PLUME%LeapYears ) LOY = 365
              kend = NINT(24.0*3600.0/dels) * LOY

           ELSE IF ( TRIM(cable_user%MetType) .EQ. 'bios' ) THEN

              kend = NINT(24.0*3600.0/dels) * LOY

           ELSE IF ( TRIM(cable_user%MetType) .EQ. 'cru' ) THEN
              ! TRENDY experiment using CRU-NCEP
              IF ( CALL1 ) THEN
                 CALL CPU_TIME(etime)
                 CALL CRU_INIT( CRU )

                 dels      = CRU%dtsecs
                 koffset   = 0
                 leaps = .false.         ! No leap years in CRU-NCEP
                 exists%Snowf = .false.  ! No snow in CRU-NCEP, so ensure it will
                                         ! be determined from temperature in CABLE
                 write(str1,'(i4)') CurYear
                 str1 = adjustl(str1)
                 write(str2,'(a)') '01'
                 str2 = adjustl(str2)
                 write(str3,'(a)') '01'
                 str3 = adjustl(str3)
                 timeunits="seconds since "//trim(str1)//"-"//trim(str2)//"-"//trim(str3)//" 00:00:00"
                 calendar = "noleap"
              ENDIF
              LOY = 365
              kend = NINT(24.0*3600.0/dels) * LOY
           ELSE IF ( TRIM(cable_user%MetType) .EQ. 'site' ) THEN
              ! site experiment eg AmazonFace (spinup or  transient run type)
              IF ( CALL1 ) THEN
                 CALL CPU_TIME(etime)
                 CALL site_INIT( site )
                 write(str1,'(i4)') CurYear
                 str1 = adjustl(str1)
                 write(str2,'(i02)') 1
                 str2 = adjustl(str2)
                 write(str3,'(i02)') 1
                 str3 = adjustl(str3)
                 timeunits="seconds since "//trim(str1)//"-"//trim(str2)//"-"//trim(str3)//" 00:00:00"
                 calendar = 'standard'
              ENDIF

              ! get koffset to add to time-step of sitemet
              IF (TRIM(site%RunType)=='historical') THEN
                 MetYear = CurYear
                 LOY = 365
                 IF (IS_LEAPYEAR(MetYear)) LOY = 366
                 kend = NINT(24.0*3600.0/dels) * LOY
              ELSEIF (TRIM(site%RunType)=='spinup' .OR. TRIM(site%RunType)=='transient') THEN
                 ! setting met year so we end the spin-up at the end of the site data-years.
                 MetYear = site%spinstartyear + &
                      MOD(CurYear- &
                      (site%spinstartyear-(site%spinendyear-site%spinstartyear +1)*100), &
                      (site%spinendyear-site%spinstartyear +1))
                 LOY = 365
                 kend = NINT(24.0*3600.0/dels) * LOY

              ENDIF
              write(*,*) 'MetYear: ', MetYear
              write(*,*) 'Simulation Year: ', CurYear
              koffset_met = 0
              if (MetYear .gt. site%spinstartyear) then
                 DO Y = site%spinstartyear, MetYear-1
                    LOYtmp = 365
                    IF (IS_LEAPYEAR(Y)) LOYtmp = 366
                    koffset_met = koffset_met + INT( REAL(LOYtmp) * 86400./REAL(dels) )
                 ENDDO
              endif
              ! LOY = 365
              ! IF (IS_LEAPYEAR(MetYear)) LOY = 366
              ! kend = NINT(24.0*3600.0/dels) * LOY
              ! ! get koffset to add to time-step of sitemet
              ! IF (TRIM(site%RunType)=='historical') THEN
              !    MetYear = CurYear
              ! ELSEIF (TRIM(site%RunType)=='spinup' .OR. TRIM(site%RunType)=='transient') THEN
              !    ! setting met year so we end the spin-up at the end of the site data-years.
              !    MetYear = site%spinstartyear + &
              !         MOD(CurYear- &
              !         (site%spinstartyear-(site%spinendyear-site%spinstartyear +1)*100), &
              !         (site%spinendyear-site%spinstartyear +1))
              ! ENDIF
              ! write(*,*) 'MetYear: ', MetYear
              ! write(*,*) 'Simulation Year: ', CurYear
              ! koffset_met = 0
              ! if (MetYear .gt. site%spinstartyear) then
              !    DO Y = site%spinstartyear, MetYear-1
              !       LOYtmp = 365
              !       IF (IS_LEAPYEAR(Y)) LOYtmp = 366
              !       koffset_met = koffset_met + INT( REAL(LOYtmp) * 86400./REAL(dels) )
              !    ENDDO
              ! endif
           ENDIF ! cable_user%MetType

           ! somethings (e.g. CASA-CNP) only need to be done once per day
           ktauday=INT(24.0*3600.0/dels)

           ! Checks where parameters and initialisations should be loaded from.
           ! If they can be found in either the met file or restart file, they will
           ! load from there, with the met file taking precedence. Otherwise, they'll
           ! be chosen from a coarse global grid of veg and soil types, based on
           ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
           IF ( CALL1 ) THEN
              IF (cable_user%POPLUC) CALL LUC_EXPT_INIT(LUC_EXPT)

              !! vh_js !!
              CALL load_parameters( met, air, ssnow, veg,climate,bgc,          &
                   soil, canopy, rough, rad, sum_flux,                  &
                   bal, logn, vegparmnew, casabiome, casapool,          &
                   casaflux, sum_casapool, sum_casaflux, &
                   casamet, casabal, phen, POP, spinup,               &
                   C%EMSOIL, C%TFRZ, LUC_EXPT, POPLUC, BLAZE, SIMFIRE, &
                   c13o2flux, c13o2pools, sum_c13o2pools, c13o2luc)
              ! 13C
              if (cable_user%c13o2) then
                 allocate(gpp(size(canopy%An,1),size(canopy%An,2)))
                 allocate(Ra(size(canopy%An,1)))
                 allocate(diff(size(canopy%An,1),size(canopy%An,2)))
                 allocate(casasave(c13o2pools%ntile,c13o2pools%npools))
                 if (cable_user%popluc) allocate(lucsave(c13o2luc%nland,c13o2luc%npools))
              endif

              if (cable_user%POPLUC .and. trim(cable_user%POPLUC_RunType) .eq. 'static') cable_user%POPLUC= .FALSE.

              ! Having read the default parameters, if this is a bios run we will now
              ! overwrite the subset of them required for bios.
              if ( trim(cable_user%MetType) .eq. 'bios' ) call cable_bios_load_params(soil)

              ! Open output file:
              IF (.NOT.CASAONLY) THEN
                 IF ( TRIM(filename%out) .EQ. '' ) THEN
                    IF ( CABLE_USER%YEARSTART .GT. 0 ) THEN
                       WRITE(dum, FMT="(I4,'_',I4)") CABLE_USER%YEARSTART, CABLE_USER%YEAREND
                       filename%out = TRIM(filename%path)//'/'//&
                            TRIM(cable_user%RunIden)//'_'//&
                            TRIM(dum)//'_cable_out.nc'
                    ELSE
                       filename%out = TRIM(filename%path)//'/'//&
                            TRIM(cable_user%RunIden)//'_cable_out.nc'
                    ENDIF
                 ENDIF
                 IF (RRRR.EQ.1) THEN
                    CALL nullify_write() ! nullify pointers
                    CALL open_output_file( dels, soil, veg, bgc, rough )
                 ENDIF
              ENDIF

              ssnow%otss_0 = ssnow%tgg(:,1)
              ssnow%otss = ssnow%tgg(:,1)
              ssnow%tss = ssnow%tgg(:,1)
              canopy%fes_cor = 0.
              canopy%fhs_cor = 0.
              met%ofsd = 0.1

              CALL zero_sum_casa(sum_casapool, sum_casaflux)
              ! 13C
              if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
              count_sum_casa = 0

              spinConv = .FALSE. ! initialise spinup convergence variable
              IF (.NOT.spinup) spinConv=.TRUE.

              if (cable_user%call_climate) then
                 call alloc_cbm_var(climate,mp,ktauday)
                 call climate_init(climate, mp, ktauday)
                 if (.not.cable_user%climate_fromzero) call read_climate_restart_nc(climate, ktauday)
              endif

              if (trim(cable_user%MetType) .eq. 'cru') &
                   casamet%glai = 1.0  ! initialise glai for use in cable_roughness

              ! additional params needed for BLAZE
              if ( trim(cable_user%MetType) .eq. 'bios' ) call cable_bios_load_climate_params(climate)


              IF ((icycle>0) .AND. spincasa) THEN
                 write(*,*) 'EXT spincasacnp enabled with mloop=', mloop
                 CALL spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                      casaflux,casamet,casabal,phen,POP,climate,LALLOC, c13o2flux, c13o2pools)
                 SPINon   = .FALSE.
                 SPINconv = .FALSE.
              ELSEIF ( casaonly .AND. (.NOT. spincasa)) THEN
                 CALL CASAONLY_LUC(dels,kstart,kend,veg,soil,casabiome,casapool, &
                      casaflux,casamet,casabal,phen,POP,climate,LALLOC, LUC_EXPT, POPLUC, &
                      sum_casapool, sum_casaflux, c13o2flux, c13o2pools, sum_c13o2pools, c13o2luc)
                 SPINon   = .FALSE.
                 SPINconv = .FALSE.
                 ktau     = kend
              ENDIF

              IF (cable_user%CALL_BLAZE) THEN
                 !CLN ?VH is rad%lat/lon below correct?
                 CALL INI_BLAZE( cable_user%CALL_POP, cable_user%BURNT_AREA, &
                      cable_user%BLAZE_TSTEP, mland, rad%latitude(landpt(:)%cstart), &
                      rad%longitude(landpt(:)%cstart), BLAZE )

                 !CLN IF ( .NOT. spinup) CALL READ_BLAZE_RESTART(...)

                 IF ( TRIM(cable_user%BURNT_AREA) == "SIMFIRE" ) THEN
                    CALL INI_SIMFIRE(mland,cable_user%SIMFIRE_REGION,SIMFIRE, &
                         climate%modis_igbp(landpt(:)%cstart) ) !CLN here we need to check for the SIMFIRE biome setting
                 ENDIF
              ENDIF

           ENDIF ! CALL1

           ! globally (WRT code) accessible kend through USE cable_common_module
           kwidth_gl = int(dels)
           kend_gl   = kend
           knode_gl  = 0

           IF (casaonly) THEN
              ! CALL1 = .FALSE.
              EXIT
           ENDIF

           ! time step loop over ktau
           DO ktau=kstart, kend

              ! increment total timstep counter
              ktau_tot = ktau_tot + 1

              ! globally (WRT code) accessible kend through USE cable_common_module
              ktau_gl = ktau_tot

              !idoy =INT( MOD(REAL(CEILING(REAL((ktau+koffset)/ktauday))),REAL(LOY)))

              idoy =INT(MOD(REAL(CEILING((real(ktau+koffset))/real(ktauday))),REAL(LOY)))
              !  write(*,*) 'idoy: ', ktauday, ktau, koffset, idoy, (REAL((real(ktau+koffset))/real(ktauday)))

              IF ( idoy .EQ. 0 ) idoy = LOY

              ! needed for CASA-CNP
              nyear =INT((kend+koffset)/(LOY*ktauday))

              ! Get met data and LAI, set time variables.
              ! Rainfall input may be augmented for spinup purposes:
              IF ( TRIM(cable_user%MetType) .EQ. 'plum' ) THEN
                 IF (( .NOT. CASAONLY ) .OR. (CASAONLY.and.CALL1))  THEN
                    CALL PLUME_MIP_GET_MET(PLUME, MET, YYYY, ktau, kend, &
                         (YYYY.EQ.CABLE_USER%YearEnd .AND. ktau.EQ.kend))

                 ENDIF
              ELSE IF ( TRIM(cable_user%MetType) .EQ. 'bios' ) THEN
                 IF (( .NOT. CASAONLY ).OR. (CASAONLY.and.CALL1))  THEN
                    CALL  cable_bios_read_met(MET, CurYear, ktau, kend, &
                         (YYYY.EQ.CABLE_USER%YearEnd .AND. ktau.EQ.kend), dels )
                 END IF
              ELSE IF ( TRIM(cable_user%MetType) .EQ. 'cru' ) THEN
                 IF (( .NOT. CASAONLY ).OR. (CASAONLY.and.CALL1))  THEN
                    CALL CRU_GET_SUBDIURNAL_MET(CRU, met, &
                         YYYY, ktau, kend, &
                         YYYY.EQ.CABLE_USER%YearEnd)
                 ENDIF
              ELSE
                 IF (TRIM(cable_user%MetType) .EQ. 'site') &
                      CALL get_met_data( spinup, spinConv, met, soil, &
                      rad, veg, kend, dels, C%TFRZ, ktau+koffset_met, &
                      kstart+koffset_met )
                 IF (TRIM(cable_user%MetType) .EQ. '') &
                      CALL get_met_data( spinup, spinConv, met, soil, &
                      rad, veg, kend, dels, C%TFRZ, ktau+koffset, &
                      kstart+koffset )
                 IF (TRIM(cable_user%MetType) .EQ. 'site' ) THEN
                    CALL site_get_CO2_Ndep(site)
                    WHERE (met%ca .eq. fixedCO2 /1000000.0)  ! not read in metfile
                       met%ca = site%CO2 / 1.e+6
                    ENDWHERE
                    met%Ndep = site%Ndep  *1000./10000./365. ! kg ha-1 y-1 > g m-2 d-1
                    met%Pdep = site%Pdep  *1000./10000./365. ! kg ha-1 y-1 > g m-2 d-1
                    met%fsd = max(met%fsd,0.0)
                 ENDIF
              ENDIF ! cable_user%MetType

              IF (TRIM(cable_user%MetType) .EQ. '') THEN
                 CurYear = met%year(1)
                 IF ( leaps .AND. IS_LEAPYEAR(CurYear) ) THEN
                    LOY = 366
                 ELSE
                    LOY = 365
                 ENDIF
              ENDIF
              met%ofsd = met%fsd(:,1) + met%fsd(:,2)
              canopy%oldcansto = canopy%cansto
              ! Zero out lai where there is no vegetation acc. to veg. index
              WHERE (veg%iveg(:) .GE. 14) veg%vlai = 0.

              ! 13C
              ! if (cable_user%c13o2) c13o2flux%ca = 0.992_r_2 * real(met%ca,r_2) ! * vpdbc13 / vpdbc13 ! -8 permil
              if (cable_user%c13o2) then
                 if ((CurYear < c13o2_atm_syear) .or. (CurYear > c13o2_atm_eyear)) then
                    write(*,*) 'Current year ', CurYear, 'not in atmospheric delta-13C (min/max): ', &
                         c13o2_atm_syear, c13o2_atm_eyear
                    stop
                 endif
                 c13o2flux%ca = (c13o2_delta_atm(CurYear) + 1.0_r_2) * real(met%ca,r_2) ! * vpdbc13 / vpdbc13
              endif

              ! At first time step of year, set tile area according to updated LU areas
              ! and zero casa fluxes

              IF (ktau == 1) THEN
                 if (icycle>1) CALL casa_cnpflux(casaflux,casapool,casabal,.TRUE.)
                 if (CABLE_USER%POPLUC) CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
              ENDIF

              IF ( .NOT. CASAONLY ) THEN

                 ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
                 IF (l_vcmaxFeedbk) then
                    IF (MOD(ktau,ktauday) == 1) THEN
                       CALL casa_feedback( ktau, veg, casabiome, casapool, casamet, climate, ktauday )
                    ENDIF
                 ELSE !JK: finite gm only effective if l_vcmaxFeedbk = .TRUE.
                    veg%vcmax_shade = veg%vcmax
                    veg%ejmax_shade = veg%ejmax

                    veg%vcmax_sun = veg%vcmax
                    veg%ejmax_sun = veg%ejmax
                 ENDIF

                 IF (l_laiFeedbk .and. icycle>0) veg%vlai(:) = casamet%glai(:)
                 ! Call land surface scheme for this timestep, all grid points:
                 
                 CALL cbm(ktau, dels, air, bgc, canopy, met, &
                          bal, rad, rough, soil, ssnow, &
                          sum_flux, veg, climate)
                 ! 13C
                 if (cable_user%c13o2) then
                    gpp  = canopy%An + canopy%Rd
                    Ra   = isoratio(c13o2flux%ca, real(met%ca,r_2), 1.0_r_2)
                    diff = canopy%An - (spread(real(met%ca,r_2),2,mf)-canopy%ci) * &
                         (1.0_r_2/(1.0_r_2/canopy%gac+1.0_r_2/canopy%gbc+1.0_r_2/canopy%gsc))
                    ! if (any(abs(canopy%An*diff) > epsilon(1._r_2))) then
                    !    print*, 'VV01 ', veg%frac4
                    !    print*, 'VV02 ', canopy%vcmax
                    !    print*, 'VV03 ', canopy%gammastar
                    !    print*, 'VV04 ', canopy%gsc
                    !    print*, 'VV05 ', canopy%gbc
                    !    print*, 'VV06 ', canopy%gac
                    !    print*, 'VV07 ', met%ca
                    !    print*, 'VV08 ', canopy%ci
                    !    print*, 'VV09 ', canopy%An
                    !    print*, 'VV10 ', (spread(real(met%ca,r_2),2,mf)-canopy%ci) * &
                    !         (1.0_r_2/(1.0_r_2/canopy%gac+1.0_r_2/canopy%gbc+1.0_r_2/canopy%gsc))
                    !    print*, 'VV11 ', diff
                    ! endif
                    do ileaf=1, mf
                       if (cable_user%c13o2_simple_disc) then
                          call c13o2_discrimination_simple( &
                               ! -- Input
                               ! isc3
                               real(dels,r_2), canopy%isc3, &
                               ! GPP and Leaf respiration
                               gpp(:,ileaf), canopy%Rd(:,ileaf), &
                               ! Ambient and stomatal CO2 concentration
                               real(met%ca,r_2), canopy%ci(:,ileaf), &
                               ! leaf temperature
                               real(canopy%tlf,r_2), &
                               ! Ambient isotope ratio
                               Ra, &
                               ! -- Inout
                               ! Starch pool and its isotope ratio
                               c13o2flux%Vstarch, c13o2flux%Rstarch, &
                               ! -- Output
                               ! discrimination
                               c13o2flux%Disc(:,ileaf), &
                               ! 13CO2 flux
                               c13o2flux%An(:,ileaf) )
                       else
                          call c13o2_discrimination( &
                               ! -- Input
                               real(dels,r_2), canopy%isc3, &
                               ! Photosynthesis variables
                               ! Vcmax<< because of temperature dependence of Vcmax
                               canopy%vcmax(:,ileaf), gpp(:,ileaf), canopy%Rd(:,ileaf), canopy%gammastar(:,ileaf), &
                               ! Could use Vcmax25?
                               ! real(veg%vcmax,r_2), gpp(:,ileaf), canopy%Rd(:,ileaf), canopy%gammastar(:,ileaf), &
                               ! CO2 concentrations
                               real(met%ca,r_2), canopy%ci(:,ileaf), &
                               ! Conductances
                               canopy%gac(:,ileaf), canopy%gbc(:,ileaf), canopy%gsc(:,ileaf), &
                               ! leaf temperature
                               real(canopy%tlf,r_2), &
                               ! Ambient isotope ratio
                               Ra, &
                               ! -- Inout
                               ! Starch pool and isotope ratios of pools for respiration
                               c13o2flux%Vstarch, c13o2flux%Rstarch, &
                               c13o2flux%Rsucrose(:,ileaf), c13o2flux%Rphoto(:,ileaf), &
                               ! -- Output
                               ! discrimination
                               c13o2flux%Disc(:,ileaf), &
                               ! 13CO2 flux
                               c13o2flux%An(:,ileaf) )
                       endif ! cable_user%c13o2_simple_disc
                    end do ! ileaf=1:mf
                    ! c13o2flux%An = 1.005_r_2 * canopy%An !  * vpdbc13 / vpdbc13 ! Test 5 permil
                 endif ! cable_user%c13o2

                 if (cable_user%CALL_climate) then
                    call cable_climate(ktau_tot, kstart, kend, ktauday, idoy, LOY, met, &
                                       climate, canopy, veg, ssnow, air, rad, dels, mp)
                 endif

                 ssnow%smelt = ssnow%smelt*dels
                 ssnow%rnof1 = ssnow%rnof1*dels
                 ssnow%rnof2 = ssnow%rnof2*dels
                 ssnow%runoff = ssnow%runoff*dels

              ELSE IF ( IS_CASA_TIME("dread", yyyy, ktau, kstart, koffset, kend, ktauday, logn) ) THEN
                 ! CLN Read from file instead
                 WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart+koffset)/(LOY*ktauday))
                 ncfile       = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
                 casa_it = NINT( REAL(ktau / ktauday) )
                 CALL read_casa_dump( ncfile, casamet, casaflux, phen, climate, c13o2flux, casa_it, kend, .FALSE. )
              ENDIF

              !jhan this is insufficient testing. condition for
              !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
              IF(icycle >0 .OR.       CABLE_USER%CASA_DUMP_WRITE ) THEN
                 !! vh_js !!
                 CALL bgcdriver( ktau, kstart, kend, dels, met,                     &
                      ssnow, canopy, veg, soil, climate, casabiome,                  &
                      casapool, casaflux, casamet, casabal,                 &
                      phen, pop, spinConv, spinup, ktauday, idoy, loy,              &
                      CABLE_USER%CASA_DUMP_READ, CABLE_USER%CASA_DUMP_WRITE,   &
                      LALLOC, c13o2flux, c13o2pools )

                 IF(MOD((ktau-kstart+1),ktauday)==0) THEN ! end of day

                    IF ( cable_user%CALL_BLAZE ) THEN
                       CALL BLAZE_ACCOUNTING(BLAZE, met, climate, ktau, dels, YYYY, idoy)

                       call blaze_driver(blaze%ncells,blaze, simfire, casapool, casaflux, &
                            casamet, climate, real(shootfrac), idoy, YYYY, 1)

                       call write_blaze_output_nc( BLAZE, ktau.EQ.kend .AND. YYYY.EQ.cable_user%YearEnd)
                    ENDIF

                 ENDIF

                 IF(MOD((ktau-kstart+1),ktauday)==0) THEN ! end of day

                    !mpidiff
                    ! update time-aggregates of casa pools and fluxes
                    CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                         & .TRUE. , .FALSE., 1)
                    ! 13C
                    if (cable_user%c13o2) call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, .true., .false., 1)
                    count_sum_casa = count_sum_casa + 1
                 ENDIF


                 IF( ((MOD((ktau-kstart+1),ktauday)==0) .AND.  &
                      MOD((ktau-kstart+1)/ktauday,LOY)==0) )THEN ! end of year
                    IF (CABLE_USER%POPLUC) THEN
                       ! Dynamic LUC
                       CALL LUCdriver(casabiome, casapool, casaflux, POP, LUC_EXPT, POPLUC, veg, c13o2pools)
                    ENDIF

                    ! one annual time-step of POP
                    CALL POPdriver(casaflux,casabal,veg, POP)

                    IF (CABLE_USER%POPLUC) THEN
                       ! Dynamic LUC: update casa pools according to LUC transitions
                       ! 13C
                       if (cable_user%c13o2) call c13o2_save_luc(casapool, popluc, casasave, lucsave)
                       CALL POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
                       ! 13C
                       if (cable_user%c13o2) then
#ifdef __C13DEBUG__
                          call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc, casapool)
#else
                          call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc)
#endif
                       endif
                       ! Dynamic LUC: write output
                       CALL WRITE_LUC_OUTPUT_NC( POPLUC, YYYY, ( YYYY.EQ.cable_user%YearEnd ))
                    ENDIF

                    !!CALL BLAZE_TURNOVER()

                    !! CLN BLAZE TURNOVER

                 ENDIF

              ENDIF
              ! WRITE CASA OUTPUT
              IF (icycle >0) THEN

                 IF ( IS_CASA_TIME("write", yyyy, ktau, kstart, koffset, kend, ktauday, logn) ) THEN
                    ctime = ctime + 1
                    !mpidiff
                    CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                         .FALSE., .TRUE., count_sum_casa)
                    ! 13C
                    if (cable_user%c13o2) &
                         call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, .false., .true., count_sum_casa)
                    CALL WRITE_CASA_OUTPUT_NC (veg, casamet, sum_casapool, casabal, sum_casaflux, &
                         CASAONLY, ctime, ( ktau.EQ.kend .AND. YYYY .EQ. cable_user%YearEnd.AND. RRRR .EQ.NRRRR ) )
                    ! 13C
                    if (cable_user%c13o2) then
                       if (ctime == 1) then
                          call c13o2_create_output(casamet, sum_c13o2pools, c13o2_outfile_id, c13o2_vars, c13o2_var_ids)
                       endif
                       call c13o2_write_output(c13o2_outfile_id, c13o2_vars, c13o2_var_ids, ctime, sum_c13o2pools)
                       if ( (ktau == kend) .and. (YYYY == cable_user%YearEnd) .and. (RRRR == NRRRR) ) &
                            call c13o2_close_output(c13o2_outfile_id)
                    end if
                    !mpidiff
                    count_sum_casa = 0
                    CALL zero_sum_casa(sum_casapool, sum_casaflux)
                    ! 13C
                    if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
                 ENDIF


                 IF (((.NOT.spinup).OR.(spinup.AND.spinConv)).and. &
                      MOD((ktau-kstart+1),ktauday)==0) THEN
                    IF ( CABLE_USER%CASA_DUMP_WRITE )  THEN

                       IF (TRIM(cable_user%MetType).EQ.'' .OR. &
                            TRIM(cable_user%MetType).EQ.'site' ) THEN

                          WRITE(CYEAR,FMT="(I4)") CurYear
                       ELSE
                          !CLN CHECK FOR LEAP YEAR
                          WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart)/(LOY*ktauday))
                       ENDIF
                       ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'

                       IF (TRIM(cable_user%MetType).EQ.'' ) THEN
                          CALL write_casa_dump(ncfile, casamet , casaflux, phen, climate, c13o2flux, &
                               int(climate%doy), LOY)
                       ELSE
                          CALL write_casa_dump(ncfile, casamet , casaflux, phen, climate, c13o2flux, &
                               idoy, kend/ktauday)
                       ENDIF

                    ENDIF
                 ENDIF

              ENDIF


              IF ( .NOT. CASAONLY ) THEN
                 ! sumcflux is pulled out of subroutine cbm
                 ! so that casaCNP can be called before adding the fluxes
                 ! (Feb 2008, YP)
                 CALL sumcflux( ktau, kstart, kend, dels, bgc,              &
                      canopy, soil, ssnow, sum_flux, veg,                   &
                      met, casaflux, l_vcmaxFeedbk )
              ENDIF

              ! Write timestep's output to file if either: we're not spinning up
              ! or we're spinning up and the spinup has converged:
              IF ( (.NOT. CASAONLY) .AND. spinConv ) THEN
                 !mpidiff
                 IF ( TRIM(cable_user%MetType) .EQ. 'plum'  .OR.  &
                      TRIM(cable_user%MetType) .EQ. 'cru'   .OR.  &
                      TRIM(cable_user%MetType) .EQ. 'bios'  .OR.  &
                      TRIM(cable_user%MetType) .EQ. 'gswp'  .OR.  &
                      TRIM(cable_user%MetType) .EQ. 'site' ) then
                    CALL write_output( dels, ktau_tot, met, canopy, casaflux, casapool, casamet, &
                         ssnow, rad, bal, air, soil, veg, C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux )
                 ELSE
                    CALL write_output( dels, ktau,     met, canopy, casaflux, casapool, casamet, &
                         ssnow, rad, bal, air, soil, veg, C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux )
                 ENDIF
              ENDIF


              ! dump bitwise reproducible testing data
              IF( cable_user%RUN_DIAG_LEVEL == 'zero') THEN
                 IF (.NOT.CASAONLY) THEN
                    IF((.NOT.spinup).OR.(spinup.AND.spinConv))                        &
                         CALL cable_diag( iDiagZero, "FLUXES", mp, kend, ktau,                        &
                         knode_gl, "FLUXES",                         &
                         canopy%fe + canopy%fh )
                 ENDIF
              ENDIF

              ! Check this run against standard for quasi-bitwise reproducability
              ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
              IF(cable_user%consistency_check) THEN

                 count_bal = count_bal +1;
                 new_sumbal = new_sumbal + SUM(bal%wbal)/mp +  SUM(bal%ebal)/mp
                 new_sumfpn = new_sumfpn + SUM(canopy%fpn)/mp
                 new_sumfe = new_sumfe + SUM(canopy%fe)/mp
                 if (ktau == kend) write(*,*) ''
                 if (ktau == kend) write(*,*) "time-space-averaged energy & water balances"
                 if (ktau == kend) write(*,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                      sum(bal%ebal_tot)/mp/count_bal, sum(bal%wbal_tot)/mp/count_bal
                 if (ktau == kend) write(*,*) "time-space-averaged latent heat and net photosynthesis"
                 if (ktau == kend) write(*,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s]",  &
                      new_sumfe/count_bal, new_sumfpn/count_bal
                 if (ktau == kend) write(logn,*) ''
                 if (ktau == kend) write(logn,*) "time-space-averaged energy & water balances"
                 if (ktau == kend) write(logn,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                      sum(bal%ebal_tot)/mp/count_bal, sum(bal%wbal_tot)/mp/count_bal
                 if (ktau == kend) write(logn,*) "time-space-averaged latent heat and net photosynthesis"
                 if (ktau == kend) write(logn,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s]",  &
                      new_sumfe/count_bal, new_sumfpn/count_bal


                 ! vh ! commented code below detects Nans in evaporation flux and stops if there are any.
                 do kk=1,mp
                    if( canopy%fe(kk).NE.( canopy%fe(kk))) THEN
                       write(*,*) 'fe nan', kk, ktau,met%qv(kk), met%precip(kk),met%precip_sn(kk), &
                            met%fld(kk), met%fsd(kk,:), met%tk(kk), met%ua(kk), ssnow%potev(kk), met%pmb(kk), &
                            canopy%ga(kk), ssnow%tgg(kk,:), canopy%fwsoil(kk)


                       stop
                    endif
                    if ( casaflux%cnpp(kk).NE. casaflux%cnpp(kk)) then
                       write(*,*) 'npp nan', kk, ktau,  casaflux%cnpp(kk)
                       !stop

                    endif


                    !if (canopy%fwsoil(kk).eq.0.0) then
                    !   write(*,*) 'zero fwsoil', ktau, canopy%fpn(kk)
                    !endif


                 enddo

                 IF( ktau == kend ) THEN
                    nkend = nkend+1

                    IF( ABS(new_sumbal-trunk_sumbal) < 1.e-7) THEN

                       PRINT *, ""
                       PRINT *, &
                            "NB. Offline-serial runs spinup cycles:", nkend
                       PRINT *, &
                            "Internal check shows this version reproduces the trunk sumbal"

                    ELSE

                       PRINT *, ""
                       PRINT *, &
                            "NB. Offline-serial runs spinup cycles:", nkend
                       PRINT *, &
                            "Internal check shows in this version new_sumbal != trunk sumbal"
                       PRINT *, &
                            "Writing new_sumbal to the file:", TRIM(Fnew_sumbal)

                       OPEN( 12, FILE = Fnew_sumbal )
                       WRITE( 12, '(F20.7)' ) new_sumbal  ! written by previous trunk version
                       CLOSE(12)

                    ENDIF
                 ENDIF

              ENDIF

              CALL1 = .FALSE.


           END DO ! END Do loop over timestep ktau

           CALL1 = .FALSE.

           !jhan this is insufficient testing. condition for
           !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
           ! see if spinup (if conducting one) has converged:
           !IF(spinup.AND..NOT.spinConv) THEN
           IF(spinup.AND.(.NOT.spinConv).AND.(.NOT.CASAONLY) ) THEN

              ! Write to screen and log file:
              WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),&
                   ' of data set complete...'
              WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',                &
                   INT(ktau_tot/kend), ' of data set complete...'


              ! IF not 1st run through whole dataset:
!!$           IF( MOD( ktau_tot, kend ) .EQ. 0 .AND. ktau_Tot .GT. kend .AND. &
!!$                YYYY.EQ. CABLE_USER%YearEnd .OR. ( NRRRR .GT. 1 .AND. &
!!$                RRRR.EQ. NRRRR) ) THEN

              IF( MOD( ktau_tot, kend ) .EQ. 0 .AND. ktau_Tot .GT. kend .AND. &
                   YYYY.EQ. CABLE_USER%YearEnd ) THEN

                 ! evaluate spinup
                 IF( ANY( ABS(ssnow%wb-soilMtemp)>delsoilM).OR.               &
                      ANY( ABS(ssnow%tgg-soilTtemp)>delsoilT) ) THEN

                    ! No complete convergence yet
                    PRINT *, 'ssnow%wb : ', ssnow%wb
                    PRINT *, 'soilMtemp: ', soilMtemp
                    PRINT *, 'ssnow%tgg: ', ssnow%tgg
                    PRINT *, 'soilTtemp: ', soilTtemp

                 ELSE ! spinup has converged

                    spinConv = .TRUE.
                    ! Write to screen and log file:
                    WRITE(*,'(A33)') ' Spinup has converged - final run'
                    WRITE(logn,'(A52)')                                       &
                         ' Spinup has converged - final run - writing all data'
                    WRITE(logn,'(A37,F8.5,A28)')                              &
                         ' Criteria: Change in soil moisture < ',             &
                         delsoilM, ' in any layer over whole run'
                    WRITE(logn,'(A40,F8.5,A28)' )                             &
                         '           Change in soil temperature < ',           &
                         delsoilT, ' in any layer over whole run'
                 END IF

              ELSE ! allocate variables for storage

                 IF (.NOT.ALLOCATED(soilMtemp)) ALLOCATE(  soilMtemp(mp,ms) )
                 IF (.NOT.ALLOCATED(soilTtemp)) ALLOCATE(  soilTtemp(mp,ms) )

              END IF

              ! store soil moisture and temperature
              IF ( YYYY.EQ. CABLE_USER%YearEnd ) THEN
                 soilTtemp = ssnow%tgg
                 soilMtemp = REAL(ssnow%wb)
              ENDIF

           ELSE

              ! if not spinning up, or spin up has converged, exit:
              IF ( SpinOn ) THEN
                 write(*,*) "setting SPINON -> FALSE", YYYY, RRRR
                 SPINon = .FALSE.
              END IF

           END IF


           IF((.NOT.(spinup.OR.casaonly)).OR.(spinup.AND.spinConv)) THEN
              IF (icycle > 0) THEN
                 CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
              END IF

           ENDIF

           IF ( .NOT. (spinup.OR.casaonly) .OR. spinconv ) THEN
              if ( NRRRR .GT. 1 ) THEN
                 RYEAR = YYYY + ( CABLE_USER%YearEnd - CABLE_USER%YearStart + 1 ) &
                      * ( RRRR - 1 )
              ELSE
                 RYEAR = YYYY
              END if
              IF ( cable_user%CALL_POP.and.POP%np.gt.0 ) then

                 IF (TRIM(cable_user%POP_out).eq.'epi') THEN
                    write(*,*) 'writing episodic POP output'
                    CALL POP_IO( pop, casamet, RYEAR, 'WRITE_EPI', &
                         (YYYY.EQ.CABLE_USER%YearEnd .AND. RRRR.EQ.NRRRR) )
                 ENDIF
              ENDIF

              !! CLN WRT BLAZE_OUT here

           ENDIF

           ! Close met data input file:
           IF ( TRIM(cable_user%MetType) .EQ. "gswp" .AND. &
                RRRR .EQ. NRRRR ) THEN
              CALL close_met_file
              IF ( YYYY .EQ. CABLE_USER%YearEnd .AND. &
                   NRRRR .GT. 1 ) DEALLOCATE ( GSWP_MID )
           ENDIF

           ! set tile area according to updated LU areas
           IF (CABLE_USER%POPLUC) THEN
              CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
           ENDIF

           IF ((icycle.gt.0).AND.(.NOT.casaonly)) THEN
              ! re-initalise annual flux sums
              casabal%FCgppyear=0.0
              casabal%FCrpyear=0.0
              casabal%FCnppyear=0
              casabal%FCrsyear=0.0
              casabal%FCneeyear=0.0
           ENDIF
           CALL CPU_TIME(etime)
           PRINT *, 'Finished. ', etime, ' seconds needed for year'

           ! 13C - While testing
           if (cable_user%c13o2) then
              call c13o2_print_delta_flux(c13o2flux)
              call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
              if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
           endif

        END DO YEAR

     END DO NREP

  END DO SPINLOOP

  IF ( SpinConv .AND. .NOT. CASAONLY ) THEN
     ! Close output file and deallocate main variables:
     CALL close_output_file( bal, air, bgc, canopy, met,                      &
          rad, rough, soil, ssnow,                                            &
          sum_flux, veg )
  ENDIF

  IF ( cable_user%CALL_POP.and.POP%np.gt.0 ) THEN
     !mpidiff
     IF ( CASAONLY .or. cable_user%pop_fromzero &
          .or.TRIM(cable_user%POP_out).eq.'ini' ) THEN
        CALL POP_IO( pop, casamet, RYEAR+1, 'WRITE_INI', .TRUE.)
     ELSE
        CALL POP_IO( pop, casamet, RYEAR+1, 'WRITE_RST', .TRUE.)
     ENDIF
  ENDIF

  !!CLN BLAZE WRITE RST

  IF (icycle > 0) THEN
     !CALL casa_poolout( ktau, veg, soil, casabiome,              &
     !     casapool, casaflux, casamet, casabal, phen )
     CALL write_casa_restart_nc( casamet, casapool,casaflux,phen, CASAONLY )
     ! 13C
     if (cable_user%c13o2) then
        call c13o2_write_restart_pools(c13o2pools)
        if (cable_user%POPLUC) call c13o2_write_restart_luc(c13o2luc)
        ! While testing
        call c13o2_print_delta_flux(c13o2flux)
        call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
        if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
     endif
  END IF

  IF (cable_user%POPLUC .AND. .NOT. CASAONLY ) THEN
     CALL WRITE_LUC_RESTART_NC ( POPLUC, YYYY )
  ENDIF

  IF ( .NOT. CASAONLY ) THEN
     ! Write restart file if requested:
     IF (output%restart) &
          CALL create_restart(logn, dels, ktau, soil, veg, ssnow, canopy, rough, rad, bgc, bal, met)
     !mpidiff
     if (cable_user%CALL_climate) &
          CALL WRITE_CLIMATE_RESTART_NC(climate, ktauday)
     ! 13C
     if (cable_user%c13o2) call c13o2_write_restart_flux(c13o2flux)
     !--- LN ------------------------------------------[
  ENDIF

  IF ( TRIM(cable_user%MetType) .NE. "gswp" .AND. &
       TRIM(cable_user%MetType) .NE. "bios" .AND. &
       TRIM(cable_user%MetType) .NE. "plum" .AND. &
       TRIM(cable_user%MetType) .NE. "cru" ) CALL close_met_file

  if (trim(cable_user%MetType) == 'cru') call cru_close(CRU)
  
  if (cable_user%POPLUC) call close_luh2(LUC_EXPT)
  
  call cpu_time(etime)
  write(logn,*) 'Finished. ', etime, ' seconds needed for ', kend,' hours'
  close(logn) ! Close log file
  write(*,*) 'Finished. ', etime, ' seconds needed for ', kend,' hours'

END PROGRAM cable_offline_driver


! ***************************************************************************************


SUBROUTINE prepareFiles(ncciy)

  USE cable_IO_vars_module, ONLY: logn, gswpfile

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ncciy

  WRITE(logn,*) 'CABLE offline global run using gswp forcing for ', ncciy
  PRINT *,      'CABLE offline global run using gswp forcing for ', ncciy

  CALL renameFiles(logn,gswpfile%rainf,ncciy,'rainf')
  CALL renameFiles(logn,gswpfile%snowf,ncciy,'snowf')
  CALL renameFiles(logn,gswpfile%LWdown,ncciy,'LWdown')
  CALL renameFiles(logn,gswpfile%SWdown,ncciy,'SWdown')
  CALL renameFiles(logn,gswpfile%PSurf,ncciy,'PSurf')
  CALL renameFiles(logn,gswpfile%Qair,ncciy,'Qair')
  CALL renameFiles(logn,gswpfile%Tair,ncciy,'Tair')
  CALL renameFiles(logn,gswpfile%wind,ncciy,'wind')

END SUBROUTINE prepareFiles


SUBROUTINE renameFiles(logn, inFile, ncciy, inName)

  IMPLICIT NONE

  INTEGER,            INTENT(IN)    :: logn
  CHARACTER(LEN=200), INTENT(INOUT) :: inFile
  INTEGER,            INTENT(IN)    :: ncciy
  CHARACTER(LEN=*),   INTENT(IN)    :: inName

  INTEGER :: nn
  INTEGER :: idummy

  nn = INDEX(inFile,'19')
  READ(inFile(nn:nn+3),'(i4)') idummy
  WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
  WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)

END SUBROUTINE renameFiles


! ***************************************************************************************
! subroutine for reading LU input data, zeroing biomass in empty secondary forest tiles
! and tranferring LUC-based age weights for secondary forest to POP structure
SUBROUTINE LUCdriver( casabiome, casapool, casaflux, POP, LUC_EXPT, POPLUC, veg, c13o2pools )

  USE cable_def_types_mod, ONLY: r_2, veg_parameter_type, mland
  USE cable_carbon_module
  USE cable_common_module, ONLY: CABLE_USER, is_casa_time, CurYear
  USE cable_IO_vars_module, ONLY: logn, landpt, patch
  USE casadimension
  USE casaparm
  USE casavariable
  USE POP_Types,  Only: POP_TYPE
  USE POPMODULE,            ONLY: POPStep, POP_init_single
  USE TypeDef,              ONLY: i4b, dp
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE, read_LUH2, &
       ptos,ptog,stog,gtos,grassfrac, pharv, smharv, syharv
  USE POPLUC_Types
  USE POPLUC_Module, ONLY: POPLUCStep, POPLUC_weights_Transfer, WRITE_LUC_OUTPUT_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, READ_LUC_RESTART_NC
  ! 13C
  use cable_c13o2_def, only: c13o2_pool

  IMPLICIT NONE

  TYPE(casa_biome),         INTENT(INOUT) :: casabiome
  TYPE(casa_pool),          INTENT(INOUT) :: casapool
  TYPE(casa_flux),          INTENT(INOUT) :: casaflux
  TYPE(POP_TYPE),           INTENT(INOUT) :: POP
  TYPE(LUC_EXPT_TYPE),      INTENT(INOUT) :: LUC_EXPT
  TYPE(POPLUC_TYPE),        INTENT(INOUT) :: POPLUC
  TYPE(veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
  ! 13C
  type(c13o2_pool),         intent(inout) :: c13o2pools

  integer ::  k, j, l, yyyy

  yyyy = CurYear

  LUC_EXPT%CTSTEP = yyyy -  LUC_EXPT%FirstYear + 1

  CALL READ_LUH2(LUC_EXPT)

  DO k=1,mland
     POPLUC%ptos(k)   = LUC_EXPT%INPUT(ptos)%VAL(k)
     POPLUC%ptog(k)   = LUC_EXPT%INPUT(ptog)%VAL(k)
     POPLUC%stog(k)   = LUC_EXPT%INPUT(stog)%VAL(k)
     POPLUC%gtop(k)   = 0.0_r_2
     POPLUC%gtos(k)   = LUC_EXPT%INPUT(gtos)%VAL(k)
     POPLUC%pharv(k)  = LUC_EXPT%INPUT(pharv)%VAL(k)
     POPLUC%smharv(k) = LUC_EXPT%INPUT(smharv)%VAL(k)
     POPLUC%syharv(k) = LUC_EXPT%INPUT(syharv)%VAL(k)
     POPLUC%thisyear  = yyyy
  ENDDO

  ! zero secondary forest tiles in POP where secondary forest area is zero
  DO k=1,mland
     if ((POPLUC%frac_primf(k)-POPLUC%frac_forest(k))==0.0_r_2 &
          .and. (.not.LUC_EXPT%prim_only(k))) then
        j = landpt(k)%cstart+1
        do l=1,size(POP%Iwood)
           if( POP%Iwood(l) == j) then
              CALL POP_init_single(POP,veg%disturbance_interval,l)
              exit
           endif
        enddo

        casapool%cplant(j,leaf) = 0.01_r_2
        casapool%nplant(j,leaf)= casabiome%ratioNCplantmin(veg%iveg(j),leaf)* casapool%cplant(j,leaf)
        casapool%pplant(j,leaf)= casabiome%ratioPCplantmin(veg%iveg(j),leaf)* casapool%cplant(j,leaf)

        casapool%cplant(j,froot) = 0.01_r_2
        casapool%nplant(j,froot)= casabiome%ratioNCplantmin(veg%iveg(j),froot)* casapool%cplant(j,froot)
        casapool%pplant(j,froot)= casabiome%ratioPCplantmin(veg%iveg(j),froot)* casapool%cplant(j,froot)

        casapool%cplant(j,wood) = 0.01_r_2
        casapool%nplant(j,wood)= casabiome%ratioNCplantmin(veg%iveg(j),wood)* casapool%cplant(j,wood)
        casapool%pplant(j,wood)= casabiome%ratioPCplantmin(veg%iveg(j),wood)* casapool%cplant(j,wood)
        casaflux%frac_sapwood(j) = 1.0_r_2

        ! 13C
        if (cable_user%c13o2) then
           c13o2pools%cplant(j,leaf)  = 0.01_r_2 ! * vpdbc13 / vpdbc13 ! Divide by 13C
           c13o2pools%cplant(j,wood)  = 0.01_r_2 ! * vpdbc13 / vpdbc13 ! so that about same numerical precision as 12C
           c13o2pools%cplant(j,froot) = 0.01_r_2 ! * vpdbc13 / vpdbc13 !
        endif
     endif
  ENDDO

  CALL POPLUCStep(POPLUC,yyyy)

  CALL POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

END SUBROUTINE LUCdriver
