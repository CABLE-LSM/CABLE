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

  use cable_def_types_mod
  use cable_io_vars_module, only: logn, gswpfile, ncciy, leaps, &
       verbose, fixedCO2, output, check, patchout, soilparmnew, &
       timeunits, exists, calendar, landpt
  use cable_common_module,  only: ktau_gl, kend_gl, knode_gl, cable_user, &
       cable_runtime, filename, &
       redistrb, wiltParam, satuParam, CurYear, &
       IS_LEAPYEAR, IS_CASA_TIME, calcsoilalbedo, get_unit, &
       report_version_no, kwidth_gl
  use cable_data_module,    only: driver_type, point2constants
  use cable_input_module,   only: open_met_file, load_parameters, get_met_data, close_met_file, &
       ncid_rain, ncid_snow, ncid_lw, ncid_sw, ncid_ps, ncid_qa, ncid_ta, ncid_wd
  use cable_output_module,  only: create_restart, open_output_file, write_output, close_output_file
  use cable_write_module,   only: nullify_write
  use cable_cbm_module
  use cable_diag_module
  !mpidiff
  use cable_climate_mod

  ! modules related to CASA-CNP
  use casadimension,        only: icycle
  use casavariable,         only: casafile, casa_biome, casa_pool, casa_flux, casa_timeunits, &
       !mpidiff
       casa_met, casa_balance, zero_sum_casa, update_sum_casa, print_casa_var
  use phenvariable,         only: phen_variable
  use casa_cable,           only: bgcdriver, POPdriver, read_casa_dump, write_casa_dump, casa_feedback, sumcflux
  use cable_spincasacnp,    only: spincasacnp
  use cable_casaonly_luc,   only: casaonly_luc
  use casa_inout,           only: casa_cnpflux, casa_fluxout, write_casa_restart_nc, write_casa_output_nc

  !! vh_js !!
  ! modules related to POP
  use POP_Types,     only: POP_TYPE
  use POPLUC_Types,  only: POPLUC_Type
  use POPLUC_Module, only: WRITE_LUC_OUTPUT_NC, WRITE_LUC_OUTPUT_GRID_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, POPLUC_set_patchfrac
  use POP_Constants, only: rshootfrac
  use cable_pop_io,  only: pop_io

  ! Fire Model BLAZE
  use blaze_drv,     only: blaze_driver
  use BLAZE_MOD,     only: TYPE_BLAZE, INI_BLAZE, BLAZE_ACCOUNTING, WRITE_BLAZE_OUTPUT_NC
  use SIMFIRE_MOD,   only: TYPE_SIMFIRE, INI_SIMFIRE

  ! gm
  use cable_adjust_JV_gm_module, only: read_gm_LUT, LUT_VcmaxJmax, LUT_gm, LUT_Vcmax, LUT_Rd

  ! 13C
  use cable_c13o2_def,         only: c13o2_delta_atm, c13o2_flux, c13o2_pool, c13o2_luc, &
       c13o2_update_sum_pools, c13o2_zero_sum_pools
  use cable_c13o2,             only: c13o2_save_luc, c13o2_update_luc, &
       c13o2_write_restart_flux, c13o2_write_restart_pools, c13o2_write_restart_luc, &
       c13o2_create_output, c13o2_write_output, c13o2_close_output, c13o2_nvars_output, &
       c13o2_sanity_pools, c13o2_sanity_luc
  use cable_c13o2,             only: c13o2_print_delta_flux, c13o2_print_delta_pools, c13o2_print_delta_luc
  use mo_isotope,              only: isoratio ! vpdbc13
  use mo_c13o2_photosynthesis, only: c13o2_discrimination_simple, c13o2_discrimination
  use mo_utils,                  only: eq, ne
  ! use mo_isotope,              only: delta1000

  ! PLUME-MIP only
  use CABLE_PLUME_MIP,      only: PLUME_MIP_TYPE, PLUME_MIP_GET_MET, &
       PLUME_MIP_INIT

  use CABLE_CRU,            only: CRU_TYPE, CRU_GET_SUBDIURNAL_MET, CRU_INIT, cru_close
  use CABLE_site,           only: site_TYPE, site_INIT, site_GET_CO2_Ndep

  ! BIOS only
  use cable_bios_met_obs_params,   only:  cable_bios_read_met, cable_bios_init, &
                                          cable_bios_load_params, &
                                          cable_bios_load_climate_params

  ! LUC_EXPT only
  use CABLE_LUC_EXPT, only: LUC_EXPT_TYPE, LUC_EXPT_INIT, close_luh2

#ifdef __NAG__
  use F90_UNIX
#endif

  implicit none

  ! CABLE namelist: model configuration, runtime/user switches
  character(len=200), parameter :: cable_namelist='cable.nml'

  ! timing variables
  integer, parameter :: kstart = 1   ! start of simulation
  integer, parameter :: mloop  = 30  ! CASA-CNP PreSpinup loops
  integer :: LALLOC ! allocation coefficient for passing to spincasa

  integer :: &
       ktau, &            ! increment equates to timestep, resets if spinning up
       ktau_tot, &        ! NO reset when spinning up, total timesteps by model
       kend, &            ! no. of time steps in run
                          !CLN kstart = 1, &  ! timestep to start at
       koffset = 0, &     ! timestep to start at
       koffset_met = 0, & !offfset for site met data ('site' only)
       ktauday, &         ! day counter for CASA-CNP
       idoy, &            ! day of year (1:365) counter for CASA-CNP
       nyear, &           ! year counter for CASA-CNP
       casa_it, &         ! number of calls to CASA-CNP
       YYYY, &            !
       RYEAR, &           !
       RRRR, &            !
       NRRRR, &           !
       LOY, &             ! days in year
       count_sum_casa     ! number of time steps over which casa pools &
  integer :: ctime ! time counter for casacnp
  !and fluxes are aggregated (for output)

  real :: dels ! time step size in seconds

  integer, dimension(:,:), allocatable :: GSWP_MID
  character(len=9) :: dum
  character(len=4) :: str1

  ! CABLE variables
  type(met_type)       :: met     ! met input variables
  type(air_type)       :: air     ! air property variables
  type(canopy_type)    :: canopy  ! vegetation variables
  type(radiation_type) :: rad     ! radiation variables
  type(roughness_type) :: rough   ! roughness varibles
  type(balances_type)  :: bal     ! energy and water balance variables
  type(soil_snow_type) :: ssnow   ! soil and snow variables
  !mpidiff
  type(climate_type)   :: climate ! climate variables

  ! CABLE parameters
  ! JK: C is used for both driver_type and i_canopy_type constants throughout the code,
  !     should be avoided ?
  type(soil_parameter_type) :: soil ! soil parameters
  type(veg_parameter_type)  :: veg  ! vegetation parameters
  type(driver_type)         :: C    ! constants used locally

  type(sum_flux_type)  :: sum_flux ! cumulative flux variables
  type(bgc_pool_type)  :: bgc  ! carbon pool variables

  ! CASA-CNP variables
  type(casa_biome)     :: casabiome
  type(casa_pool)      :: casapool
  type(casa_flux)      :: casaflux
  type(casa_pool)      :: sum_casapool
  type(casa_flux)      :: sum_casaflux
  type(casa_met)       :: casamet
  type(casa_balance)   :: casabal
  type(phen_variable)  :: phen
  !! vh_js !!
  type(POP_TYPE)        :: POP
  type(POPLUC_TYPE)     :: POPLUC
  type(PLUME_MIP_TYPE)  :: PLUME
  type(CRU_TYPE)        :: CRU
  type(site_TYPE)       :: site
  type(LUC_EXPT_TYPE)   :: LUC_EXPT
  character             :: cyear*4
  character             :: ncfile*99

  ! BLAZE variables
  type(TYPE_BLAZE)    :: BLAZE
  type(TYPE_SIMFIRE)  :: SIMFIRE

  ! 13C
  type(c13o2_flux) :: c13o2flux
  type(c13o2_pool) :: c13o2pools, sum_c13o2pools
  type(c13o2_luc)  :: c13o2luc
  real(r_2), dimension(:,:), allocatable :: casasave
  real(r_2), dimension(:,:), allocatable :: lucsave
  ! I/O
  integer :: c13o2_outfile_id
  character(len=40), dimension(c13o2_nvars_output) :: c13o2_vars
  integer,           dimension(c13o2_nvars_output) :: c13o2_var_ids
  ! discrimination
  integer :: ileaf
  real(r_2), dimension(:,:), allocatable :: gpp ! , diff
  real(r_2), dimension(:),   allocatable :: Ra
  ! delta-13C of atmospheric CO2
  integer            :: iunit, ios
  real               :: iyear
  integer            :: c13o2_atm_syear, c13o2_atm_eyear
  character(len=100) :: header

  ! declare vars for switches (default .FALSE.) etc declared thru namelist
  logical :: &
       vegparmnew    = .false., & ! using new format input file (BP dec 2007)
       spinup        = .false., & ! model spinup to soil state equilibrium?
       spinConv      = .false., & ! has spinup converged?
       spincasainput = .false., & ! TRUE: SAVE input req'd to spin CASA-CNP
                                  ! FALSE: READ input to spin CASA-CNP
       spincasa      = .false., & ! TRUE: CASA-CNP Will spin mloop times,
                                  ! FALSE: no spin up
       l_casacnp     = .false., & ! using CASA-CNP with CABLE
       l_laiFeedbk   = .false., & ! using prognostic LAI
       l_vcmaxFeedbk = .false., & ! using prognostic Vcmax
       CASAONLY      = .false., & ! ONLY Run CASA-CNP
       CALL1         = .true., &
       SPINon        = .true.

  logical :: CASA_TIME
  logical :: first_casa_write
  logical :: liseod, liseoy ! is end of day, is end of year

  real :: &
       delsoilM, & ! allowed variation in soil moisture for spin up
       delsoilT    ! allowed variation in soil temperature for spin up

 integer :: Metyear, Y, LOYtmp

  ! temporary storage for soil moisture/temp. in spin up mode
  real, allocatable, dimension(:,:) :: &
       soilMtemp, &
       soilTtemp

  ! timing
  real:: etime ! Declare the type of etime(), For receiving user and system time, total time

  !___ unique unit/file identifiers for cable_diag: arbitrarily 5 here
  integer :: iDiagZero=0

  ! switches etc defined thru namelist (by default cable.nml)
  namelist /cablenml/ &
       filename, &       ! TYPE, containing input filenames
       vegparmnew, &     ! use new soil param. method
       soilparmnew, &    ! use new soil param. method
       calcsoilalbedo, & ! albedo considers soil color Ticket #27
       spinup, &         ! spinup model (soil) to steady state
       delsoilM, &
       delsoilT, &
       output, &
       patchout, &
       check, &
       verbose, &
       leaps, &
       logn, &
       fixedCO2, &
       spincasainput, &
       spincasa, &
       l_casacnp, &
       l_laiFeedbk, &
       l_vcmaxFeedbk, &
       icycle, &
       casafile, &
       ncciy, &
       gswpfile, &
       redistrb, &
       wiltParam, &
       satuParam, &
       cable_user  ! additional USER switches

  !mpidiff
  integer :: kk

  ! Vars for standard for quasi-bitwise reproducability b/n runs
  ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
  character(len=30), parameter :: &
       Ftrunk_sumbal = ".trunk_sumbal", &
       Fnew_sumbal   = "new_sumbal"

  real(r_2) :: &
       trunk_sumbal = 0.0_r_2, & !
       new_sumbal   = 0.0_r_2, &
       new_sumfpn   = 0.0_r_2, &
       new_sumfe    = 0.0_r_2

  integer :: nkend=0
  integer :: ioerror
  integer :: count_bal = 0

  ! command line arguments
  integer :: narg, len1, len2
  character(len=500) :: arg1
  character(len=200) :: arg2

  ! END header

  ! Open, read and close the namelist file.
  open(10, file=cable_namelist)
  read(10, nml=cablenml)   !where nml=cable defined above
  close(10)

  ! Open, read and close the consistency check file.
  ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
  if (cable_user%consistency_check) then
     open(11, file=Ftrunk_sumbal, status='old', action='READ', iostat=ioerror)
     if (ioerror==0) then
        read(11,*) trunk_sumbal  ! written by previous trunk version
     end if
     close(11)
  end if

  ! Open log file:
  open(logn, file=filename%log)

  call report_version_no(logn)

  narg = command_argument_count()
  if (narg > 0) then
     call get_command_argument(1, arg1, len1)
     filename%met = arg1(1:len1)
     call get_command_argument(2, arg2, len2)
     casafile%cnpipool = arg2(1:len2)
  end if

  ! INITIALISATION depending on nml settings
  if (trim(cable_user%MetType) == 'gswp') then
     if ((CABLE_USER%YearStart == 0) .and. (ncciy > 0)) then
        CABLE_USER%YearStart = ncciy
        CABLE_USER%YearEnd = ncciy
     else if ((CABLE_USER%YearStart == 0) .and. (ncciy == 0)) then
        write(*,*) 'undefined start year for gswp met: '
        write(*,*) 'enter value for ncciy or'
        write(*,*) '(CABLE_USER%YearStart and  CABLE_USER%YearEnd) in cable.nml'

        write(logn,*) 'undefined start year for gswp met: '
        write(logn,*) 'enter value for ncciy or'
        write(logn,*) '(CABLE_USER%YearStart and  CABLE_USER%YearEnd) in cable.nml'
        stop 201
     end if
  end if

  CurYear = CABLE_USER%YearStart

  if (icycle >= 11) then
     icycle                     = icycle - 10
     CASAONLY                   = .true.
     CABLE_USER%CASA_DUMP_READ  = .true.
     CABLE_USER%CASA_DUMP_WRITE = .false.
  else if (icycle == 0) then
     CABLE_USER%CASA_DUMP_READ  = .false.
     spincasa                   = .false.
     !! vh_js !!
     CABLE_USER%CALL_POP        = .false.
     CABLE_USER%CALL_BLAZE      = .false.
  end if

  !! vh_js !!
  if (icycle > 0) then
     l_casacnp = .true.
  else
     l_casacnp = .false.
  end if

  !! vh_js !! suggest LALLOC should ulitmately be a switch in the .nml file
  if (CABLE_USER%CALL_POP) then
     LALLOC = 3 ! for use with POP: makes use of pipe model to partition between stem and leaf
  else
     LALLOC = 0 ! default
  end if

  ! IF ( .NOT. spinup ) THEN
  !    IF ( spincasa ) THEN
  !       spincasa = .FALSE.
  !       WRITE(*,*)    "spinup == .FALSE. -> spincasa set to .F."
  !       WRITE(logn,*) "spinup == .FALSE. -> spincasa set to .F."
  !    END IF
  ! END IF

  if (trim(cable_user%MetType) == 'gpgs') then
     leaps = .true.
     calendar = "standard"
     cable_user%MetType = 'gswp'
  end if

  cable_runtime%offline = .true.

  ! associate pointers used locally with global definitions
  call point2constants(C)

  if (l_casacnp  .and. ((icycle == 0) .or. (icycle > 3))) then
     stop 'icycle must be 1 to 3 when using casaCNP'
  end if
  ! IF( ( l_laiFeedbk .OR. l_vcmaxFeedbk ) ) &
  !    STOP 'casaCNP required to get prognostic LAI or Vcmax'
  if ((l_laifeedbk .or. l_vcmaxfeedbk) .and. (.not. l_casacnp)) then
     stop 'casaCNP required to get prognostic LAI or Vcmax'
  end if
  if (l_vcmaxFeedbk .and. (icycle < 2)) then
     stop 'icycle must be 2 to 3 to get prognostic Vcmax'
  end if
  if ((icycle > 0) .and. (.not. soilparmnew)) then
     stop 'casaCNP must use new soil parameters'
  end if

  NRRRR = merge(max(CABLE_USER%CASA_NREP, 1), 1, CASAONLY)

  ! casa time count
  ctime = 0
  first_casa_write = .true.

  ! INISTUFF

  ! Open met data and get site information from netcdf file. (NON-GSWP ONLY!)
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites.

  ! IF ( TRIM(cable_user%MetType) .NE. "gswp" .AND. &
  !      TRIM(cable_user%MetType) .NE. "gpgs" .AND. &
  !      TRIM(cable_user%MetType) .NE. "plume" .AND. &
  !      TRIM(cable_user%MetType) .NE. "cru") THEN
  if ((trim(cable_user%MetType) .eq. 'site') .or. &
       (trim(cable_user%MetType) .eq. '')) then
     call open_met_file( dels, koffset, kend, spinup, C%TFRZ )
     if ((koffset /= 0) .and. CABLE_USER%CALL_POP) then
        write(*,*) "When using POP, episode must start at Jan 1st!"
        stop 202
     end if
  else if (NRRRR > 1) then
     if (.not. allocated(GSWP_MID)) &
          allocate(GSWP_MID(8, CABLE_USER%YearStart:CABLE_USER%YearEnd))
  end if

  ! gm lookup table
  if (cable_user%explicit_gm .and. (len_trim(cable_user%gm_LUT_file) > 1)) then
     write(*,*) 'Reading gm LUT file'
     call read_gm_LUT(trim(cable_user%gm_LUT_file), LUT_VcmaxJmax, LUT_gm, &
          LUT_Vcmax, LUT_Rd)
  end if

  ! 13C
  ! Read atmospheric delta-13C values
  if (cable_user%c13o2) then
     ! get start and end year
     call get_unit(iunit)
     open(iunit, file=trim(cable_user%c13o2_delta_atm_file), status="old", &
          action="read")
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
     c13o2_delta_atm = c13o2_delta_atm / 1000._r_2
  end if

  ! outer loop - spinup loop no. ktau_tot :
  RYEAR    = 0
  ktau_tot = 0
  SPINon   = .true.
  ! YearStart = CABLE_USER%YearStart
  ! YearEnd   = CABLE_USER%YearEnd
  ! cable_user%CASA_SPIN_ENDYEAR

  SPINLOOP: do while (SPINon)

     NREP: do RRRR=1, NRRRR

        if (trim(cable_user%MetType) == "bios") then
           call CPU_time(etime)
           call cable_bios_init(dels,curyear,met,kend,ktauday)
           koffset   = 0
           leaps = .true.
           write(str1,'(i4)') curyear
           str1 = adjustl(str1)
           timeunits = "seconds since " // trim(str1) // "-01-01 00:00:00"
           calendar = "standard"
           casa_timeunits = "days since " // trim(str1) // "-01-01 00:00:00"
        end if

        write(*,*) "CABLE_USER%YearStart,  CABLE_USER%YearEnd", CABLE_USER%YearStart,  CABLE_USER%YearEnd

        YEAR: do YYYY=CABLE_USER%YearStart,  CABLE_USER%YearEnd

           CurYear = YYYY
           if (leaps .and. IS_LEAPYEAR(YYYY)) then
              LOY = 366
           else
              LOY = 365
           end if

           if (trim(cable_user%MetType) == 'gswp') then
              ! GSWP run
              ncciy = CurYear
              write(*,*) 'Looking for global offline run info.'
              call preparefiles(ncciy)
              if (RRRR == 1) then
                 call open_met_file(dels, koffset, kend, spinup, C%TFRZ)
                 if (leaps .and. is_leapyear(YYYY) .and. (kend == 2920)) then
                    stop 'LEAP YEAR INCOMPATIBILITY WITH INPUT MET !!!'
                 end if
                 if (NRRRR > 1) then
                    GSWP_MID(1, YYYY) = ncid_rain
                    GSWP_MID(2, YYYY) = ncid_snow
                    GSWP_MID(3, YYYY) = ncid_lw
                    GSWP_MID(4, YYYY) = ncid_sw
                    GSWP_MID(5, YYYY) = ncid_ps
                    GSWP_MID(6, YYYY) = ncid_qa
                    GSWP_MID(7, YYYY) = ncid_ta
                    GSWP_MID(8, YYYY) = ncid_wd
                 end if
              else
                 ncid_rain = GSWP_MID(1, YYYY)
                 ncid_snow = GSWP_MID(2, YYYY)
                 ncid_lw   = GSWP_MID(3, YYYY)
                 ncid_sw   = GSWP_MID(4, YYYY)
                 ncid_ps   = GSWP_MID(5, YYYY)
                 ncid_qa   = GSWP_MID(6, YYYY)
                 ncid_ta   = GSWP_MID(7, YYYY)
                 ncid_wd   = GSWP_MID(8, YYYY)
                 kend      = ktauday * LOY
              end if
           else if (trim(cable_user%MetType) == 'plume') then
              ! PLUME experiment setup using WATCH
              if (CALL1) then
                 call cpu_time(etime)
                 call plume_mip_init( plume )
                 dels      = PLUME%dt
                 koffset   = 0
                 leaps = PLUME%LeapYears
                 write(str1,'(i4)') CurYear
                 str1 = adjustl(str1)
                 timeunits="seconds since " // trim(str1) // "-01-01 00:00:00"
                 if (leaps) then
                    calendar = "standard"
                 else
                    calendar = "noleap"
                 end if
                 casa_timeunits="days since " // trim(str1) // "-01-01 00:00:00"
              end if
              if (.not. PLUME%LeapYears) LOY = 365
              kend = nint(24.0 * 3600.0 / dels) * LOY
           else if (trim(cable_user%MetType) == 'bios') then
              ! BIOS run
              kend = nint(24.0 * 3600.0 / dels) * LOY
           else if (trim(cable_user%MetType) == 'cru') then
              ! TRENDY experiment using CRU-NCEP
              if (CALL1) then
                 call cpu_time(etime)
                 call cru_init(cru)
                 dels         = cru%dtsecs
                 koffset      = 0
                 if (CRU%MetVersion == "VERIFY_2021") then
                    leaps = .true.
                    calendar = "standard"
                    if (IS_LEAPYEAR(CurYear)) then
                       LOY = 366
                    else
                       LOY = 365
                    end if
                 else
                    leaps    = .false. ! No leap years in CRU-NCEP
                    calendar = "noleap"
                    LOY      = 365
                 end if
                 exists%Snowf = .false. ! No snow in CRU-NCEP, so ensure it will
                                        ! be determined from temperature in CABLE
                 write(str1, '(i4)') CurYear
                 str1 = adjustl(str1)
                 timeunits = "seconds since " // trim(str1) // "-01-01 00:00:00"
                 casa_timeunits = "days since " // trim(str1) // "-01-01 00:00:00"
              end if
              kend = nint(24.0 * 3600.0 / dels) * LOY
           else if (trim(cable_user%MetType) == 'site') then
              ! site experiment eg AmazonFace (spinup or transient run type)
              if (CALL1) then
                 call cpu_time(etime)
                 call site_init(site)
                 write(str1,'(i4)') CurYear
                 str1 = adjustl(str1)
                 timeunits = "seconds since " // trim(str1) // "-01-01 00:00:00"
                 calendar  = "standard"
                 casa_timeunits="days since " // trim(str1) // "-01-01 00:00:00"
              end if
              ! get koffset to add to time-step of sitemet
              if (trim(site%RunType) == 'historical') then
                 MetYear = CurYear
                 leaps = .true.
                 LOY = 365
                 if (IS_LEAPYEAR(MetYear)) LOY = 366
                 kend = nint(24.0 * 3600.0 / dels) * LOY
              else if ((trim(site%RunType) == 'spinup') .or. &
                   trim(site%RunType) == 'transient') then
                 ! setting met year so we end the spin-up at the end
                 ! of the site data-years.
                 MetYear = site%spinstartyear + &
                      mod(CurYear - &
                      (site%spinstartyear-(site%spinendyear-site%spinstartyear +1)*100), &
                      (site%spinendyear-site%spinstartyear+1))
                 leaps = .false.
                 LOY = 365
                 kend = nint(24.0 * 3600.0 / dels) * LOY
              end if
              write(*,*) 'MetYear: ', MetYear
              write(*,*) 'Simulation Year: ', CurYear
              koffset_met = 0
              if (MetYear > site%spinstartyear) then
                 do Y = site%spinstartyear, MetYear-1
                    LOYtmp = 365
                    if (IS_LEAPYEAR(Y)) LOYtmp = 366
                    koffset_met = koffset_met + int(real(LOYtmp) * 86400./real(dels))
                 end do
              end if
           end if ! cable_user%MetType

           ! somethings (e.g. CASA-CNP) only need to be done once per day
           ktauday = int(24.0 * 3600.0 / dels)

           ! Checks where parameters and initialisations should be loaded from.
           ! If they can be found in either the met file or restart file, they will
           ! load from there, with the met file taking precedence. Otherwise, they'll
           ! be chosen from a coarse global grid of veg and soil types, based on
           ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
           if (CALL1) then
              if (cable_user%POPLUC) call LUC_EXPT_INIT(LUC_EXPT)

              ! 13C
              call load_parameters(met, air, ssnow, veg, bgc, &
                   soil, canopy, rough, rad, sum_flux, &
                   bal, logn, vegparmnew, casabiome, casapool, &
                   casaflux, sum_casapool, sum_casaflux, &
                   casamet, casabal, phen, POP, spinup, &
                   C%EMSOIL, C%TFRZ, LUC_EXPT, POPLUC, &
                   c13o2flux, c13o2pools, sum_c13o2pools, c13o2luc)

              ! 13C
              if (cable_user%c13o2) then
                 allocate(gpp(size(canopy%An, 1), size(canopy%An, 2)))
                 allocate(Ra(size(canopy%An, 1)))
                 ! allocate(diff(size(canopy%An, 1), size(canopy%An, 2)))
                 allocate(casasave(c13o2pools%ntile, c13o2pools%npools))
                 if (cable_user%popluc) &
                      allocate(lucsave(c13o2luc%nland, c13o2luc%npools))
              end if

              if (cable_user%POPLUC .and. &
                   (trim(cable_user%POPLUC_RunType) == 'static')) &
                   cable_user%POPLUC = .false.

              ! Having read the default parameters, if this is a bios run we
              ! will now overwrite the subset of them required for bios.
              if (trim(cable_user%MetType) .eq. 'bios') &
                   call cable_bios_load_params(soil)

              ! Open output file
              if (.not. CASAONLY) then
                 if (trim(filename%out) == '' ) then
                    if (CABLE_USER%YEARSTART > 0) then
                       write(dum, FMT="(I4,'_',I4)") CABLE_USER%YEARSTART, &
                            CABLE_USER%YEAREND
                       filename%out = trim(filename%path) // '/' // &
                            trim(cable_user%RunIden) // '_' // &
                            trim(dum) // '_cable_out.nc'
                    else
                       filename%out = trim(filename%path) // '/' // &
                            trim(cable_user%RunIden) // '_cable_out.nc'
                    end if
                 end if
                 if (RRRR == 1) then
                    call nullify_write() ! nullify pointers
                    call open_output_file(dels, soil, veg, bgc, rough)
                 end if
              end if

              canopy%fes_cor = 0.0_r_2
              canopy%fhs_cor = 0.0
              met%ofsd       = 0.1 ! not used

              call zero_sum_casa(sum_casapool, sum_casaflux)
              ! 13C
              if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
              count_sum_casa = 0

              spinConv = .false. ! initialise spinup convergence variable
              if (.not.spinup) spinConv=.true.

              if (cable_user%call_climate) then
                 call alloc_cbm_var(climate, mp, ktauday)
                 call zero_cbm_var(climate)
                 call climate_init(climate)
                 if (.not. cable_user%climate_fromzero) &
                      call read_climate_restart_nc(climate)
              end if

              !MC - do we need that? casamet also set only if icycle>0
              ! if (trim(cable_user%MetType) .eq. 'cru') then
              !    casamet%glai = 1.0_r_2 ! initialise glai for use in cable_roughness
              !    where (veg%iveg(:) .ge. 14) casamet%glai = 0.0_r_2
              ! end if
              if ((trim(cable_user%MetType) == 'cru') .and. (icycle > 0)) then
                 if (all(real(casamet%glai(:)) == 0.)) then
                    write(*,*) 'W A R N I N G : deleted setting casamet%glai for CRU.'
                    stop 205
                 end if
              end if

              ! additional params needed for BLAZE
              if (trim(cable_user%MetType) == 'bios') &
                   call cable_bios_load_climate_params(climate)

              if (cable_user%CALL_BLAZE) then
                 print*, "CLN BLAZE INIT"
                 call INI_BLAZE(mland, rad%latitude(landpt(:)%cstart), &
                      rad%longitude(landpt(:)%cstart), BLAZE)


                 if (trim(BLAZE%BURNT_AREA_SRC) == "SIMFIRE") then
                    print*, "CLN SIMFIRE INIT"
                    !CLN here we need to check for the SIMFIRE biome setting
                    call INI_SIMFIRE(mland ,SIMFIRE, &
                         climate%modis_igbp(landpt(:)%cstart))
                 end if
              end if

              if ((icycle>0) .and. spincasa) then
                 write(*,*) 'EXT spincasacnp enabled with mloop=', mloop
                 call spincasacnp(dels, kstart, kend, mloop, veg, soil, &
                      casabiome, casapool, casaflux, casamet, casabal, phen, &
                      POP, climate, LALLOC, c13o2flux, c13o2pools, &
                      BLAZE, SIMFIRE)
                 if (cable_user%c13o2) &
                      call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                 ! Set 13C and 12C LUC pools to 0 if any < epsilon(1.0_dp)
                 SPINon   = .false.
                 SPINconv = .false.
                 !MC - MPI sets CASAONLY = .true. here
              ! ELSE IF ( casaonly .AND. (.NOT. spincasa)) THEN
              else if (casaonly .and. (.not. spincasa) .and. &
                   cable_user%popluc) then
                 write(*,*) 'EXT CASAONLY_LUC'
                 call CASAONLY_LUC(dels, kstart, kend, veg, soil, &
                      casabiome, casapool, casaflux, casamet, casabal, phen, &
                      POP, climate, LALLOC, LUC_EXPT, POPLUC, &
                      sum_casapool, sum_casaflux, c13o2flux, c13o2pools, &
                      sum_c13o2pools, c13o2luc)
                 if (cable_user%c13o2) then
                    call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                    call c13o2_sanity_luc(popluc, c13o2luc)
                 end if
                 SPINon   = .false.
                 SPINconv = .false.
              end if

           end if ! CALL1

           ! globally (WRT code) accessible kend through USE cable_common_module
           kwidth_gl = int(dels)
           kend_gl   = kend
           knode_gl  = 0

           !MC IF (casaonly) THEN
           if (spincasa .or. casaonly) then
              ! CALL1 = .false.
              exit
           end if

           ! time step loop over ktau
           do ktau=kstart, kend

              ! increment total timstep counter
              ktau_tot = ktau_tot + 1

              ! globally (WRT code) accessible kend through USE cable_common_module
              ktau_gl = ktau_tot

              idoy = int( mod(real(ceiling((real(ktau+koffset)) / &
                   real(ktauday))), real(loy)) )
              if (idoy .eq. 0) idoy = LOY

              ! needed for CASA-CNP
              !MC - diff to MPI
              nyear = int((kend+koffset) / (LOY*ktauday))

              ! Get met data and LAI, set time variables.
              ! Rainfall input may be augmented for spinup purposes:
              if (trim(cable_user%MetType) == 'plume') then
                 if ((.not. CASAONLY) .or. (CASAONLY .and. CALL1)) then
                    call PLUME_MIP_GET_MET(PLUME, MET, YYYY, ktau, kend, &
                         ((YYYY == CABLE_USER%YearEnd) .and. (ktau == kend)))
                 end if
              else if (trim(cable_user%MetType) == 'bios') then
                 if ((.not. CASAONLY) .or. (CASAONLY .and. CALL1)) then
                    call cable_bios_read_met(MET, CurYear, ktau, dels)
                 end if
              else if (trim(cable_user%MetType) == 'cru') then
                 if ((.not. CASAONLY) .or. (CASAONLY .and. CALL1)) then
                    call CRU_GET_SUBDIURNAL_MET(CRU, met, YYYY, ktau, kend)
                 end if
              else
                 if (trim(cable_user%MetType) == 'site') &
                      call get_met_data(spinup, spinConv, met, &
                      rad, veg, dels, C%TFRZ, ktau+koffset_met, &
                      kstart+koffset_met)
                 if (cable_user%perturb_Ta) &
                      met%tk = met%tk + cable_user%Ta_perturbation
                 if (trim(cable_user%MetType) == '') &
                      call get_met_data(spinup, spinConv, met, &
                      rad, veg, dels, C%TFRZ, ktau+koffset, &
                      kstart+koffset)
                 if (trim(cable_user%MetType) == 'site') then
                    call site_get_CO2_Ndep(site)
                    where (eq(met%ca, fixedCO2 / 1000000.0))
                       ! not read in metfile
                       met%ca = site%CO2 / 1.e+6
                    END WHERE
                    ! kg ha-1 y-1 > g m-2 d-1
                    met%Ndep = site%Ndep * 1000. / 10000. / 365.
                    ! kg ha-1 y-1 > g m-2 d-1
                    met%Pdep = site%Pdep * 1000. / 10000. / 365.
                    met%fsd = max(met%fsd,0.0)
                 end if
              end if ! cable_user%MetType

              if (trim(cable_user%MetType) == '') then
                 CurYear = met%year(1)
                 if (leaps .and. IS_LEAPYEAR(CurYear)) then
                    LOY = 366
                 else
                    LOY = 365
                 end if
              end if
              met%ofsd = met%fsd(:,1) + met%fsd(:,2)
              canopy%oldcansto = canopy%cansto
              ! Zero out lai where there is no vegetation acc. to veg. index
              where (veg%iveg(:) >= 14) veg%vlai = 0.

              ! 13C
              ! if (cable_user%c13o2) c13o2flux%ca = 0.992_r_2 * real(met%ca,r_2) ! * vpdbc13 / vpdbc13 ! -8 permil
              if (cable_user%c13o2) then
                 if ((CurYear < c13o2_atm_syear) .or. (CurYear > c13o2_atm_eyear)) then
                    write(*,*) 'Current year ', CurYear, &
                         'not in atmospheric delta-13C (min/max): ', &
                         c13o2_atm_syear, c13o2_atm_eyear
                    stop 203
                 end if
                 c13o2flux%ca = (c13o2_delta_atm(CurYear) + 1.0_r_2) * &
                      real(met%ca, r_2) ! * vpdbc13 / vpdbc13
              end if

              ! At first time step of year, set tile area according to
              ! updated LU areas and zero casa fluxes
              if (ktau == 1) then
                 if (icycle > 1) &
                      call casa_cnpflux(casaflux, casapool, casabal, .true.)
                 if (CABLE_USER%POPLUC) then
                    call POPLUC_set_patchfrac(POPLUC, LUC_EXPT)
                 end if
              end if

              if ((trim(cable_user%MetType) == '') .or. &
                   (trim(cable_user%MetType) == 'site')) then
                 CASA_TIME = IS_CASA_TIME("write", metyear, ktau, kstart, &
                      koffset, ktauday, logn)
              else
                 CASA_TIME = IS_CASA_TIME("write", yyyy, ktau, kstart, &
                      koffset, ktauday, logn)
              end if
              liseod = mod((ktau-kstart+1), ktauday) == 0
              liseoy = mod((ktau-kstart+1) / ktauday, LOY) == 0

              if (.not. CASAONLY) then
                 ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
                 if (l_vcmaxFeedbk) then
                    if (mod(ktau, ktauday) == 1) then
                       call casa_feedback(ktau, veg, casabiome, casapool, &
                            casamet, climate, ktauday)
                    end if
                 else ! JK: finite gm only effective if l_vcmaxFeedbk = .TRUE.
                    veg%vcmax_shade = veg%vcmax
                    veg%ejmax_shade = veg%ejmax
                    veg%vcmax_sun = veg%vcmax
                    veg%ejmax_sun = veg%ejmax
                 end if

                 if (l_laiFeedbk .and. (icycle > 0)) &
                      veg%vlai(:) = real(casamet%glai(:))

                 ! if ((ktau == kstart) .or. (ktau == kend)) then
                 ! if ((ktau == 9) .or. (ktau == kend)) then
                 !    print*, 'CC00 ', YYYY, ktau, kstart, kend
                 !    call print_cbm_var(canopy)
                 ! end if

                 call cbm(ktau, dels, air, bgc, canopy, met, &
                      bal, rad, rough, soil, ssnow, &
                      veg, climate)

                 ! if ((ktau == kstart) .or. (ktau == kend)) then
                 ! if ((ktau == 9) .or. (ktau == kend)) then
                 !    print*, 'CC01 ', YYYY, ktau, kstart, kend
                 !    call print_cbm_var(canopy)
                 ! end if

                 ! 13C
                 if (cable_user%c13o2) then
                    gpp  = canopy%An + canopy%Rd
                    Ra   = isoratio(c13o2flux%ca, real(met%ca,r_2), 1.0_r_2)
                    ! diff = canopy%An - (spread(real(met%ca,r_2),2,mf)-canopy%ci) * &
                    !      (1.0_r_2/(1.0_r_2/canopy%gac+1.0_r_2/canopy%gbc+1.0_r_2/canopy%gsc))
                    !MCTest
                    c13o2flux%An       = canopy%An
                    ! c13o2flux%An       = 1.005_r_2 * canopy%An !  * vpdbc13 / vpdbc13 ! Test 5 permil
                    c13o2flux%Disc     = 0.0_r_2
                    c13o2flux%Vstarch  = c13o2flux%Vstarch + 1.0e-6_r_2
                    c13o2flux%Rstarch  = c13o2flux%Rstarch
                    c13o2flux%Rsucrose = c13o2flux%Rsucrose
                    c13o2flux%Rphoto   = c13o2flux%Rphoto
                    ! do ileaf=1, mf
                    !    if (cable_user%c13o2_simple_disc) then
                    !       call c13o2_discrimination_simple( &
                    !            ! -- Input
                    !            ! isc3
                    !            real(dels,r_2), canopy%isc3, &
                    !            ! GPP and Leaf respiration
                    !            gpp(:,ileaf), canopy%Rd(:,ileaf), &
                    !            ! Ambient and stomatal CO2 concentration
                    !            real(met%ca,r_2), canopy%ci(:,ileaf), &
                    !            ! leaf temperature
                    !            real(canopy%tlf,r_2), &
                    !            ! Ambient isotope ratio
                    !            Ra, &
                    !            ! -- Inout
                    !            ! Starch pool and its isotope ratio
                    !            c13o2flux%Vstarch, c13o2flux%Rstarch, &
                    !            ! -- Output
                    !            ! discrimination
                    !            c13o2flux%Disc(:,ileaf), &
                    !            ! 13CO2 flux
                    !            c13o2flux%An(:,ileaf) )
                    !    else
                    !       call c13o2_discrimination( &
                    !            ! -- Input
                    !            real(dels,r_2), canopy%isc3, &
                    !            ! Photosynthesis variables
                    !            ! Vcmax<< because of temperature dependence of Vcmax
                    !            canopy%vcmax(:,ileaf), gpp(:,ileaf), canopy%Rd(:,ileaf), canopy%gammastar(:,ileaf), &
                    !            ! Could use Vcmax25?
                    !            ! real(veg%vcmax,r_2), gpp(:,ileaf), canopy%Rd(:,ileaf), canopy%gammastar(:,ileaf), &
                    !            ! CO2 concentrations
                    !            real(met%ca,r_2), canopy%ci(:,ileaf), &
                    !            ! Conductances
                    !            canopy%gac(:,ileaf), canopy%gbc(:,ileaf), canopy%gsc(:,ileaf), &
                    !            ! leaf temperature
                    !            real(canopy%tlf,r_2), &
                    !            ! Ambient isotope ratio
                    !            Ra, &
                    !            ! -- Inout
                    !            ! Starch pool and isotope ratios of pools for respiration
                    !            c13o2flux%Vstarch, c13o2flux%Rstarch, &
                    !            c13o2flux%Rsucrose(:,ileaf), c13o2flux%Rphoto(:,ileaf), &
                    !            ! -- Output
                    !            ! discrimination
                    !            c13o2flux%Disc(:,ileaf), &
                    !            ! 13CO2 flux
                    !            c13o2flux%An(:,ileaf) )
                    !    end if ! cable_user%c13o2_simple_disc
                    ! end do ! ileaf=1:mf
                    !MCTest
                 end if ! cable_user%c13o2

                 if (cable_user%CALL_climate) &
                      call cable_climate(ktau_tot, kstart, ktauday, idoy, LOY, met, &
                      climate, canopy, veg, ssnow, rad, dels, mp)

                 ssnow%smelt  = ssnow%smelt  * dels
                 ssnow%rnof1  = ssnow%rnof1  * dels
                 ssnow%rnof2  = ssnow%rnof2  * dels
                 ssnow%runoff = ssnow%runoff * dels
                 ! if ((ktau == kstart) .or. (ktau == kend)) then
                 ! if ((ktau == 9) .or. (ktau == kend)) then
                 !    ! CALL load_parameters( met, air, ssnow, veg, bgc, &
                 !    !      soil, canopy, rough, rad, sum_flux, &
                 !    !      bal, logn, vegparmnew, casabiome, casapool, &
                 !    !      casaflux, sum_casapool, sum_casaflux, &
                 !    !      casamet, casabal, phen, POP, spinup, &
                 !    !      C%EMSOIL, C%TFRZ, LUC_EXPT, POPLUC, &
                 !    !      c13o2flux, c13o2pools, sum_c13o2pools, c13o2luc)
                 !    print*, 'MM01 ', YYYY, ktau, kstart, kend
                 !    call print_cbm_var(canopy)
                 ! end if

              else if (IS_CASA_TIME("dread", yyyy, ktau, kstart, koffset, ktauday, logn)) then
                 if (cable_user%casa_dump_read) then
                    ! CLN Read from file instead
                    write(CYEAR, fmt="(I4)") CurYear + &
                         int((ktau-kstart+koffset) / (LOY*ktauday))
                    ncfile       = trim(casafile%c2cdumppath) // 'c2c_' // CYEAR // '_dump.nc'
                    casa_it = nint(real(ktau / ktauday))
                    call read_casa_dump(ncfile, casamet, casaflux, phen, climate, c13o2flux, casa_it, kend, .false.)
                 end if
              end if ! .not. casaonly

              if ((icycle > 0) .or. CABLE_USER%CASA_DUMP_WRITE) then

                 call bgcdriver(ktau, kstart, dels, met, &
                      ssnow, canopy, veg, soil, climate, casabiome, &
                      casapool, casaflux, casamet, casabal, &
                      phen, pop, ktauday, idoy, loy, &
                      CABLE_USER%CASA_DUMP_READ, &
                      LALLOC, c13o2flux, c13o2pools)
                 if (cable_user%c13o2) &
                      call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                 ! if (any(delta1000(c13o2pools%cplant, casapool%cplant, 1.0_r_2, 0.0_r_2, tiny(1.0_r_2)) > 0.0_r_2)) then
                 !    print*, 'CC04.01 ', casapool%cplant
                 !    print*, 'CC04.02 ', c13o2pools%cplant
                 !    print*, 'CC04.03 ', delta1000(c13o2pools%cplant, casapool%cplant, 1.0_r_2, 0.0_r_2, tiny(1.0_r_2))
                 ! end if
                 ! if ((ktau == kstart) .or. (ktau == kend)) then
                 ! if ((ktau == 9) .or. (ktau == kend)) then
                 !    print*, 'MM02 ', YYYY, ktau, kstart, kend
                 !    call print_cbm_var(canopy)
                 ! end if

                 if (liseod) then ! end of day
                    if (cable_user%CALL_BLAZE) then
                       print*, "CLN CAlling BLAZE"
                       call BLAZE_ACCOUNTING(BLAZE, climate, ktau, dels, YYYY, idoy)
                       call blaze_driver(blaze%ncells, blaze, simfire, casapool, casaflux, &
                            casamet, climate, rshootfrac, idoy, YYYY, 1, POP, veg)

                       call write_blaze_output_nc(BLAZE, &
                            (ktau == kend) .and. (YYYY == cable_user%YearEnd))
                    end if
                 end if

                 if (liseod .and. liseoy) then ! end of year

                    if (CABLE_USER%POPLUC) then
                       ! Dynamic LUC
                       call LUCdriver(casabiome, casapool, casaflux, POP, LUC_EXPT, POPLUC, veg, c13o2pools)
                       if (cable_user%c13o2) &
                            call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                    end if

                    ! one annual time-step of POP
                    call POPdriver(casaflux, casabal, veg, POP)

                    if (CABLE_USER%POPLUC) then
                       ! Dynamic LUC: update casa pools according to LUC transitions
                       ! 13C
                       if (cable_user%c13o2) &
                            call c13o2_save_luc(casapool, popluc, casasave, lucsave)
                       call POP_LUC_CASA_transfer(POPLUC, POP, LUC_EXPT, &
                            casapool, casabal, casaflux, ktauday)
                       ! 13C
                       if (cable_user%c13o2) then
#ifdef __C13DEBUG__
                          call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc, casapool)
#else
                          call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc)
#endif
                          if (cable_user%c13o2) then
                             call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                             call c13o2_sanity_luc(popluc, c13o2luc)
                          end if
                       end if
                       ! Dynamic LUC: write output
                       if (output%grid(1:3) == 'lan') then
                          call WRITE_LUC_OUTPUT_NC(POPLUC, YYYY, (YYYY == cable_user%YearEnd))
                       else
                          call WRITE_LUC_OUTPUT_GRID_NC(POPLUC, YYYY, (YYYY == cable_user%YearEnd))
                       end if
                    end if ! POPLUC

                 end if ! end of day and end of year

                 if (liseod) then ! end of day
                    ctime = ctime + 1  ! update casa time
                    ! update time-aggregates of casa pools and fluxes
                    count_sum_casa = count_sum_casa + 1
                    call update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                         .true., CASA_TIME, count_sum_casa) ! end of year
                    ! 13C
                    if (cable_user%c13o2) &
                         call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, &
                         .true., CASA_TIME, count_sum_casa)
                 end if

              end if ! icycle > 0  or  casa_dump_write

              ! WRITE CASA OUTPUT
              if (icycle > 0) then

                 if (CASA_TIME) then
                    call WRITE_CASA_OUTPUT_NC(veg, casamet, sum_casapool, &
                         casabal, sum_casaflux, CASAONLY, ctime, &
                         (ktau == kend) .and. (YYYY == cable_user%YearEnd) .and. &
                         (RRRR == NRRRR))

                    ! 13C
                    if (cable_user%c13o2) then
                       if (first_casa_write) then
                          call c13o2_create_output(casamet, sum_c13o2pools, &
                               c13o2_outfile_id, c13o2_vars, c13o2_var_ids)
                          first_casa_write = .false.
                       end if
                       call c13o2_write_output(c13o2_outfile_id, c13o2_vars, &
                            c13o2_var_ids, ctime, sum_c13o2pools)
                       if ((ktau == kend) .and. (YYYY == cable_user%YearEnd) .and. (RRRR == NRRRR)) &
                            call c13o2_close_output(c13o2_outfile_id)
                    end if
                    count_sum_casa = 0
                    call zero_sum_casa(sum_casapool, sum_casaflux)
                    ! 13C
                    if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
                 end if

                 if ( ( (.not.spinup) .or. (spinup.and.spinConv) ) .and. liseod ) then
                    if ( CABLE_USER%CASA_DUMP_WRITE )  then

                       if (trim(cable_user%MetType).eq.'' .or. &
                            trim(cable_user%MetType).eq.'site' ) then
                          write(CYEAR,FMT="(I4)") CurYear
                       else
                          !CLN CHECK FOR LEAP YEAR
                          write(CYEAR,FMT="(I4)") CurYear + int((ktau-kstart)/(LOY*ktauday))
                       end if
                       ncfile = trim(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'

                       if (trim(cable_user%MetType).eq.'' ) then
                          call write_casa_dump(ncfile, casamet , casaflux, phen, climate, c13o2flux, &
                               int(climate%doy), LOY)
                       else
                          call write_casa_dump(ncfile, casamet , casaflux, phen, climate, c13o2flux, &
                               idoy, kend/ktauday)
                       end if

                    end if
                 end if

              end if ! icycle > 0


              if (.not. CASAONLY) then
                 ! sumcflux is pulled out of subroutine cbm
                 ! so that casaCNP can be called before adding the fluxes
                 ! (Feb 2008, YP)
                 call sumcflux(ktau, kstart, dels, &
                      canopy, sum_flux, &
                      casaflux, l_vcmaxFeedbk)
              end if

              ! Write timestep's output to file if either: we're not spinning up
              ! or we're spinning up and the spinup has converged:
              if ((.not. CASAONLY) .and. spinConv) then
                 if ( trim(cable_user%MetType) == 'plume' .or. &
                      trim(cable_user%MetType) == 'cru' .or. &
                      trim(cable_user%MetType) == 'bios' .or. &
                      trim(cable_user%MetType) == 'gswp' .or. &
                      trim(cable_user%MetType) == 'site' ) then
                    call write_output(dels, ktau_tot, met, canopy, casaflux, &
                         casapool, casamet, ssnow, rad, bal, air, soil, veg, &
                         C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux)
                 else
                    call write_output(dels, ktau, met, canopy, casaflux, &
                         casapool, casamet, ssnow, rad, bal, air, soil, veg, &
                         C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux )
                 end if
              end if

              ! dump bitwise reproducible testing data
              if (cable_user%RUN_DIAG_LEVEL == 'zero') then
                 if (.not. CASAONLY) then
                    if((.not. spinup) .or. (spinup .and. spinConv)) &
                         call cable_diag(iDiagZero, "FLUXES", mp, kend, ktau, &
                         knode_gl, "FLUXES", canopy%fe + canopy%fh)
                 end if
              end if

              ! Check this run against standard for quasi-bitwise reproducability
              ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
              if (cable_user%consistency_check) then
                 count_bal = count_bal + 1
                 new_sumbal = new_sumbal + sum(bal%wbal) / mp + &
                      sum(bal%ebal)/mp
                 new_sumfpn = new_sumfpn + sum(canopy%fpn) / mp
                 new_sumfe = new_sumfe + sum(canopy%fe) / mp
                 if (ktau == kend) write(*,*) ''
                 if (ktau == kend) write(*,*) "time-space-averaged energy & water balances"
                 if (ktau == kend) write(*,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                      sum(bal%ebal_tot) / mp / count_bal, sum(bal%wbal_tot) / mp / count_bal
                 if (ktau == kend) write(*,*) "time-space-averaged latent heat and net photosynthesis"
                 if (ktau == kend) write(*,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s]", &
                      new_sumfe / count_bal, new_sumfpn / count_bal
                 if (ktau == kend) write(logn,*) ''
                 if (ktau == kend) write(logn,*) "time-space-averaged energy & water balances"
                 if (ktau == kend) write(logn,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                      sum(bal%ebal_tot) / mp / count_bal, sum(bal%wbal_tot) / mp / count_bal
                 if (ktau == kend) write(logn,*) "time-space-averaged latent heat and net photosynthesis"
                 if (ktau == kend) write(logn,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s]", &
                      new_sumfe / count_bal, new_sumfpn / count_bal

                 ! vh ! commented code below detects Nans in evaporation flux and stops if there are any.
                 do kk=1, mp
                    if (ne(canopy%fe(kk), canopy%fe(kk))) then
                       write(*,*) 'fe nan', kk, ktau,met%qv(kk), met%precip(kk),met%precip_sn(kk), &
                            met%fld(kk), met%fsd(kk,:), met%tk(kk), met%ua(kk), ssnow%potev(kk), met%pmb(kk), &
                            canopy%ga(kk), ssnow%tgg(kk,:), canopy%fwsoil(kk)
                       stop 204
                    end if
                    if (ne(casaflux%cnpp(kk), casaflux%cnpp(kk))) then
                       write(*,*) 'npp nan', kk, ktau,  casaflux%cnpp(kk)
                       !stop
                    end if
                    ! if (canopy%fwsoil(kk).eq.0.0) then
                    !    write(*,*) 'zero fwsoil', ktau, canopy%fpn(kk)
                    ! end if
                 end do

                 if (ktau == kend) then
                    nkend = nkend + 1
                    if ( abs(new_sumbal-trunk_sumbal) < 1.e-7) then
                       write(*,*) ""
                       write(*,*) "NB. Offline-serial runs spinup cycles:", nkend
                       write(*,*) "Internal check shows this version reproduces the trunk sumbal"
                    else
                       write(*,*) ""
                       write(*,*) "NB. Offline-serial runs spinup cycles:", nkend
                       write(*,*) "Internal check shows in this version new_sumbal != trunk sumbal"
                       write(*,*) "Writing new_sumbal to the file:", trim(Fnew_sumbal)
                       open(12, file=Fnew_sumbal)
                       write(12, '(F20.7)') new_sumbal  ! written by previous trunk version
                       close(12)
                    end if
                 end if

              end if ! consistency_check

              CALL1 = .false.

           end do ! END Do loop over timestep ktau

           CALL1 = .false.

           ! see if spinup (if conducting one) has converged
           if (spinup .and. (.not.spinConv) .and. (.not.CASAONLY)) then

              ! Write to screen and log file:
              write(*,'(A18,I3,A24)') ' Spinning up: run ', int(ktau_tot / kend), &
                   ' of data set complete...'
              write(logn,'(A18,I3,A24)') ' Spinning up: run ', &
                   int(ktau_tot / kend), ' of data set complete...'

              ! IF not 1st run through whole dataset:
              if( (mod(ktau_tot, kend) == 0) .and. (ktau_tot > kend) .and. &
                   (YYYY == CABLE_USER%YearEnd) ) then

                 ! evaluate spinup
                 if( any(abs(ssnow%wb-soilMtemp) > delsoilM) .or. &
                      any(abs(ssnow%tgg-soilTtemp) > delsoilT) ) then
                    ! No complete convergence yet
                    write(*,*) 'ssnow%wb : ', ssnow%wb
                    write(*,*) 'soilMtemp: ', soilMtemp
                    write(*,*) 'ssnow%tgg: ', ssnow%tgg
                    write(*,*) 'soilTtemp: ', soilTtemp
                 else ! spinup has converged
                    spinConv = .true.
                    ! Write to screen and log file:
                    write(*,'(A33)') ' Spinup has converged - final run'
                    write(logn,'(A52)') &
                         ' Spinup has converged - final run - writing all data'
                    write(logn,'(A37,F8.5,A28)') &
                         ' Criteria: Change in soil moisture < ', &
                         delsoilM, ' in any layer over whole run'
                    write(logn,'(A40,F8.5,A28)' ) &
                         '           Change in soil temperature < ', &
                         delsoilT, ' in any layer over whole run'
                 end if

              else ! allocate variables for storage

                 if (.not. allocated(soilMtemp)) allocate(soilMtemp(mp,ms))
                 if (.not. allocated(soilTtemp)) allocate(soilTtemp(mp,ms))

              end if

              ! store soil moisture and temperature
              if (YYYY == CABLE_USER%YearEnd) then
                 soilTtemp = ssnow%tgg
                 soilMtemp = real(ssnow%wb)
              end if

           else

              ! if not spinning up, or spin up has converged, exit
              if (SpinOn) then
                 write(*,*) "setting SPINON -> FALSE", YYYY, RRRR
                 SPINon = .false.
              end if

           end if ! spinup .and. ...

           if ((.not. (spinup .or. casaonly)) .or. (spinup .and. spinConv)) then
              if (icycle > 0) then
                 call casa_fluxout(nyear, veg, soil, casabal, casamet)
              end if
           end if

           if (.not. (spinup.or.casaonly) .or. spinconv) then
              if (NRRRR > 1) then
                 RYEAR = YYYY &
                      + (CABLE_USER%YearEnd - CABLE_USER%YearStart + 1) &
                      * (RRRR - 1)
              else
                 RYEAR = YYYY
              end if
              if (cable_user%CALL_POP .and. (POP%np > 0)) then
                 if (trim(cable_user%POP_out) == 'epi') then
                    write(*,*) 'writing episodic POP output'
                    call POP_IO(pop, casamet, RYEAR, 'WRITE_EPI', &
                         (YYYY == CABLE_USER%YearEnd .and. (RRRR == NRRRR)))
                 end if
              end if
              !! CLN WRT BLAZE_OUT here
           end if

           ! Close met data input file:
           if ( (trim(cable_user%MetType) == "gswp") .and. &
                (RRRR .eq. NRRRR) ) then
              call close_met_file()
              if ( (YYYY == CABLE_USER%YearEnd) .and. &
                   (NRRRR == 1) ) deallocate( GSWP_MID )
           end if

           ! set tile area according to updated LU areas
           if (CABLE_USER%POPLUC) then
              call POPLUC_set_patchfrac(POPLUC, LUC_EXPT)
           end if

           if ((icycle > 0) .and. (.not. casaonly)) then
              ! re-initalise annual flux sums
              casabal%FCgppyear = 0.0_r_2
              casabal%FCrpyear  = 0.0_r_2
              casabal%FCnppyear = 0.0_r_2
              casabal%FCrsyear  = 0.0_r_2
              casabal%FCneeyear = 0.0_r_2
           end if

           call CPU_time(etime)
           write(*,*) 'Finished. ', etime, ' seconds needed for year'

           ! 13C - While testing
           if (cable_user%c13o2) then
              call c13o2_print_delta_flux(c13o2flux)
              call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
              if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
           end if

        end do YEAR

     end do NREP

  end do SPINLOOP

  if (SpinConv .and. (.not. CASAONLY)) then
     ! Close output file and deallocate main variables:
     call close_output_file(bal)
  end if

  if (cable_user%CALL_POP .and. (POP%np > 0)) then
     if ( CASAONLY .or. cable_user%pop_fromzero .or. &
          (trim(cable_user%POP_out) == 'ini') ) then
        call POP_IO(pop, casamet, RYEAR+1, 'WRITE_INI', .true.)
     else
        call POP_IO(pop, casamet, RYEAR+1, 'WRITE_RST', .true.)
     end if
  end if

  !!CLN BLAZE WRITE RST

  IF (icycle > 0) THEN
     !CALL casa_poolout( ktau, veg, soil, casabiome, &
     !     casapool, casaflux, casamet, casabal, phen )
     CALL write_casa_restart_nc(casabiome, casamet, casapool, casaflux, casabal, phen)

     ! 13C
     if (cable_user%c13o2) then
        call c13o2_write_restart_pools(casamet, c13o2pools)
        if (cable_user%POPLUC) call c13o2_write_restart_luc(popluc, c13o2luc)
        ! While testing
        call c13o2_print_delta_flux(c13o2flux)
        call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
        if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
     end if
  END IF

  IF (cable_user%POPLUC .AND. .NOT. CASAONLY ) THEN
     CALL WRITE_LUC_RESTART_NC(POPLUC)
  END IF

  IF ( .NOT. CASAONLY ) THEN
     ! Write restart file if requested:
     IF (output%restart) &
          CALL create_restart(logn, dels, ktau, soil, veg, ssnow, canopy, rad, bgc, bal)
     !mpidiff
     if (cable_user%CALL_climate) &
          CALL WRITE_CLIMATE_RESTART_NC(climate)
     ! 13C
     if (cable_user%c13o2) call c13o2_write_restart_flux(casamet, c13o2flux)
  END IF

  if ( trim(cable_user%MetType) /= "gswp" .and. &
       trim(cable_user%MetType) /= "bios" .and. &
       trim(cable_user%MetType) /= "plume" .and. &
       trim(cable_user%MetType) /= "cru" ) call close_met_file

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
  write(*,*)    'CABLE offline global run using gswp forcing for ', ncciy

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

  USE cable_def_types_mod,  ONLY: r_2, veg_parameter_type, mland
  USE cable_carbon_module
  USE cable_common_module,  ONLY: CABLE_USER, CurYear
  USE cable_IO_vars_module, ONLY: landpt
  USE casadimension
  USE casaparm
  USE casavariable
  USE POP_Types,            Only: POP_TYPE
  USE POPMODULE,            ONLY: POP_init_single
  USE CABLE_LUC_EXPT,       ONLY: LUC_EXPT_TYPE, read_LUH2, &
       ptos, ptog, stog, gtos, pharv, smharv, syharv
  USE POPLUC_Types
  USE POPLUC_Module,        ONLY: POPLUCStep, POPLUC_weights_Transfer
  ! 13C
  use cable_c13o2_def,      only: c13o2_pool
  use mo_utils,             only: eq

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
     POPLUC%ptos(k)   = real(LUC_EXPT%INPUT(ptos)%VAL(k), r_2)
     POPLUC%ptog(k)   = real(LUC_EXPT%INPUT(ptog)%VAL(k), r_2)
     POPLUC%stog(k)   = real(LUC_EXPT%INPUT(stog)%VAL(k), r_2)
     POPLUC%gtop(k)   = 0.0_r_2
     POPLUC%gtos(k)   = real(LUC_EXPT%INPUT(gtos)%VAL(k),   r_2)
     POPLUC%pharv(k)  = real(LUC_EXPT%INPUT(pharv)%VAL(k),  r_2)
     POPLUC%smharv(k) = real(LUC_EXPT%INPUT(smharv)%VAL(k), r_2)
     POPLUC%syharv(k) = real(LUC_EXPT%INPUT(syharv)%VAL(k), r_2)

     !MC - from MPI code
     ! POPLUC%ptoc(k) = real(LUC_EXPT%INPUT(ptoc)%VAL(k), r_2)
     ! POPLUC%ptoq(k) = real(LUC_EXPT%INPUT(ptoq)%VAL(k), r_2)
     ! POPLUC%stoc(k) = real(LUC_EXPT%INPUT(stoc)%VAL(k), r_2)
     ! POPLUC%stoq(k) = real(LUC_EXPT%INPUT(stoq)%VAL(k), r_2)
     ! POPLUC%ctos(k) = real(LUC_EXPT%INPUT(ctos)%VAL(k), r_2)
     ! POPLUC%qtos(k) = real(LUC_EXPT%INPUT(qtos)%VAL(k), r_2)

     POPLUC%thisyear  = yyyy

  END DO

  ! zero secondary forest tiles in POP where secondary forest area is zero
  DO k=1,mland
     if ( eq((POPLUC%frac_primf(k)-POPLUC%frac_forest(k)), 0.0_r_2) &
          .and. (.not. LUC_EXPT%prim_only(k)) ) then
        j = landpt(k)%cstart+1
        do l=1,size(POP%Iwood)
           if( POP%Iwood(l) == j) then
              CALL POP_init_single(POP,veg%disturbance_interval,l)
              exit
           end if
        end do

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
        end if
     end if
  END DO

  CALL POPLUCStep(POPLUC, yyyy)

  CALL POPLUC_weights_transfer(POPLUC, POP, LUC_EXPT)

END SUBROUTINE LUCdriver
