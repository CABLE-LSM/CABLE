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
! Purpose: Offline driver for mpi master in CABLE global run
!
! Contact: Bernard.Pak@csiro.au
!
! History: Since 1.4b, capability to run global offline (ncciy = YEAR),
!          inclusion of call to CASA-CNP (icycle>0)
!          soil_snow_type now ssnow (instead of ssoil)
!
!          MPI wrapper developed by Maciej Golebiewski (2012)
!          Modified from cable_driver.F90 in CABLE-2.0_beta r171 by B Pak
!
! ==============================================================================
! Uses:        mpi
!              cable_mpicommon
!              cable_def_types_mod
!              cable_IO_vars_module
!              cable_common_module
!              cable_data_module
!              cable_input_module
!              cable_output_module
!              cable_cbm_module
!              casadimension
!              casavariable
!              phenvariable
!
! CALLs:       point2constants
!              open_met_file
!              load_parameters
!              open_output_file
!              get_met_data
!              write_output
!              casa_poolout
!              casa_fluxout
!              create_restart
!              close_met_file
!              close_output_file
!              prepareFiles
!              find_extents
!              master_decomp
!              master_cable_params
!              master_casa_params
!              master_intypes
!              master_outtypes
!              master_casa_types
!              master_restart_types
!              master_send_input
!              master_receive
!              master_end
!              master_casa_dump_types
!              master_casa_LUC_types
!              ! 13C
!              master_c13o2_flux_params
!              master_c13o2_pool_params
!              master_c13o2_luc_params
!              master_c13o2_flux_types
!              master_c13o2_pool_types
!              master_c13o2_luc_types
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
MODULE cable_mpimaster

  use cable_mpicommon

  implicit none

  save

  private

  ! number of workers; set in master_decomp
  integer :: wnp

  ! TODO: m3d_t mat_t and vec_t to be factored out from here and from master_outtypes
  ! MPI: slices of 3D arrays
  integer, allocatable, dimension(:,:) :: m3d_t
  ! MPI: slices of matrices (2D)
  integer, allocatable, dimension(:,:) :: mat_t

  ! MPI: parts of vectors (1D)
  ! MPI: vec_t dimension is wnp; as each worker gets a single indexed with nvec blocks
  integer, allocatable, dimension(:) :: vec_t

  ! MPI derived datatype handles for sending input data to the workers
  integer, allocatable, dimension(:) :: inp_ts

  ! MPI derived datatype handles for receiving output from the workers
  integer, allocatable, dimension(:) :: recv_ts

  ! master's struct for receiving restart data from the workers
  integer, allocatable, dimension(:) :: restart_ts

  ! CASA related derived types

  ! MPI derived datatype handles for receiving casa results from the workers
  ! and restart values
  integer, allocatable, dimension(:) :: casa_ts

  ! MPI derived datatype handles for send/receiving casa dump values from the workers
  integer, allocatable, dimension(:) :: casa_dump_ts

  ! MPI derived datatype handles for send/receiving casa pool values (needed for LUC)
  !  from the workers
  integer, allocatable, dimension(:) :: casa_LUC_ts

  !CLN  ! MPI derived datatype handles for receiving casa restart values from the workers
  !CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_restart_ts

  ! climate derived type
  integer, allocatable, dimension(:) :: climate_ts

  ! MPI derived datatype handles for Sending/receiving vals results for BLAZE
  integer, allocatable, dimension(:) :: blaze_in_ts
  integer, allocatable, dimension(:) :: blaze_out_ts
  integer, allocatable, dimension(:) :: blaze_restart_ts

  ! POP related derived types

  ! MPI derived datatype handles for receiving POP results from the workers
  integer :: pop_ts

  ! 13C
  ! MPI derived datatype handles for receiving c13o2 results from the workers
  ! and restart values
  integer, allocatable, dimension(:) :: c13o2_flux_ts
  integer, allocatable, dimension(:) :: c13o2_pool_ts
  integer, allocatable, dimension(:) :: c13o2_luc_ts

  ! MPI: isend request array for scattering input data to the workers
  integer, allocatable, dimension(:) :: inp_req
  ! MPI: isend status array for scattering input data to the workers
  integer, allocatable, dimension(:,:) :: inp_stats

  ! MPI: irecv request array for gathering results from the workers
  integer, allocatable, dimension(:) :: recv_req
  ! MPI: irecv status array for gathering results from the workers
  integer, allocatable, dimension(:,:) :: recv_stats

  ! MPI: landpoints decomposition; global info used by the master process
  type(lpdecomp_t), allocatable, dimension(:) :: wland

  public :: mpidrv_master

contains

  subroutine mpidrv_master(comm)

    use mpi

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
    use cable_input_module,   only: open_met_file, load_parameters, &
         get_met_data,close_met_file
    use cable_output_module,  only: create_restart, open_output_file, &
         write_output, close_output_file
    use cable_write_module,   only: nullify_write
    ! USE cable_cbm_module
    use cable_climate_mod

    ! modules related to CASA-CNP
    use casadimension,        only: icycle
    use casavariable,         only: casafile, casa_biome, casa_pool, casa_flux, casa_timeunits, &
         casa_met, casa_balance, zero_sum_casa, update_sum_casa
    use phenvariable,         only: phen_variable
    use casa_cable,           only: write_casa_dump
    use casa_inout,           only: casa_fluxout, write_casa_restart_nc, write_casa_output_nc
    use casa_inout,           only: casa_cnpflux

    !CLN added
    ! modules related to POP
    use POP_Types,            only: POP_TYPE
    use cable_pop_io,         only: pop_io
    use POPLUC_Types,         only: POPLUC_Type
    use POPLUC_Module,        only: WRITE_LUC_OUTPUT_NC, WRITE_LUC_OUTPUT_GRID_NC, &
         POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, POPLUC_set_patchfrac, &
         READ_LUC_RESTART_NC, alloc_popluc
    ! LUC_EXPT only
    use CABLE_LUC_EXPT,       only: LUC_EXPT_TYPE, LUC_EXPT_INIT, close_luh2

    ! modules related to fire
    use BLAZE_MOD,            only: TYPE_BLAZE, INI_BLAZE, WRITE_BLAZE_OUTPUT_NC
    use BLAZE_MPI,            only: MASTER_BLAZE_TYPES ! , MASTER_SIMFIRE_TYPES
    use SIMFIRE_MOD,          only: TYPE_SIMFIRE, INI_SIMFIRE

    !MCJK - check if need in mpimaster and mpiworker
    ! gm
    use cable_adjust_JV_gm_module, only: read_gm_LUT, LUT_VcmaxJmax, LUT_gm, LUT_Vcmax, LUT_Rd
    !MCJK - check if need in mpimaster and mpiworker

    ! 13C
    use cable_c13o2_def,      only: c13o2_delta_atm, c13o2_flux, c13o2_pool, c13o2_luc, &
         c13o2_update_sum_pools, c13o2_zero_sum_pools
    use cable_c13o2,          only: c13o2_save_luc, c13o2_update_luc, &
         c13o2_write_restart_flux, c13o2_write_restart_pools, c13o2_write_restart_luc, &
         c13o2_create_output, c13o2_write_output, c13o2_close_output, c13o2_nvars_output, &
         c13o2_sanity_pools, c13o2_sanity_luc
    use cable_c13o2,          only: c13o2_print_delta_flux, c13o2_print_delta_pools, c13o2_print_delta_luc
    use mo_utils,             only: ne

    ! PLUME-MIP only
    use CABLE_PLUME_MIP,      only: PLUME_MIP_TYPE, PLUME_MIP_GET_MET, &
         PLUME_MIP_INIT

    use CABLE_CRU,            only: CRU_TYPE, CRU_GET_SUBDIURNAL_MET, CRU_INIT ! , cru_close

    ! BIOS only
    use cable_bios_met_obs_params, only: cable_bios_read_met, cable_bios_init, &
         cable_bios_load_params, cable_bios_load_climate_params

    implicit none

    ! MPI:
    integer               :: comm ! MPI communicator for comms with the workers

    ! CABLE namelist: model configuration, runtime/user switches
    character(len=200), parameter :: cable_namelist='cable.nml'

    ! timing variables
    integer, parameter ::  kstart = 1   ! start of simulation

    integer :: &
         ktau, &  ! increment equates to timestep, resets if spinning up
         ktau_tot, &  ! NO reset when spinning up, total timesteps by model
         kend, &  ! no. of time steps in run
         !CLN kstart = 1, &  ! timestep to start at
         koffset = 0, &  ! timestep to start at
         ktauday, &  ! day counter for CASA-CNP
         idoy, &  ! day of year (1:365) counter for CASA-CNP
         nyear, &  ! year counter for CASA-CNP
         YYYY, &  !
         LOY, &  ! Length of Year
         maxdiff(2), &  ! location of maximum in convergence test
         count_sum_casa  ! number of time steps over which casa pools &
    integer :: ctime ! time for casacnp

    real :: dels                    ! time step size in seconds
    character(len=9) :: dum         ! dummy char for filename generation
    character(len=4) :: str1        ! dummy char for filename generation

    ! CABLE variables
    type(met_type)       :: met     ! met input variables: see below for imet in MPI variables
    type(air_type)       :: air     ! air property variables
    type(canopy_type)    :: canopy  ! vegetation variables
    type(radiation_type) :: rad     ! radiation variables
    type(roughness_type) :: rough   ! roughness varibles
    type(balances_type)  :: bal     ! energy and water balance variables
    type(soil_snow_type) :: ssnow   ! soil and snow variables
    type(climate_type)   :: climate ! climate variables

    ! CABLE parameters
    type(soil_parameter_type) :: soil     ! soil parameters
    type(veg_parameter_type)  :: veg      ! vegetation parameters: see below for iveg in MPI variables
    type(driver_type)         :: C        ! constants used locally
    type(sum_flux_type)       :: sum_flux ! cumulative flux variables
    type(bgc_pool_type)       :: bgc      ! carbon pool variables

    ! CASA-CNP variables
    type(casa_biome)     :: casabiome
    type(casa_pool)      :: casapool
    type(casa_flux)      :: casaflux
    type(casa_pool)      :: sum_casapool
    type(casa_flux)      :: sum_casaflux
    type(casa_met)       :: casamet
    type(casa_balance)   :: casabal
    type(phen_variable)  :: phen
    type(POP_TYPE)       :: POP
    type(POPLUC_TYPE)    :: POPLUC
    type(LUC_EXPT_TYPE)  :: LUC_EXPT
    type(PLUME_MIP_TYPE) :: PLUME
    type(CRU_TYPE)       :: CRU
    character            :: cyear*4
    character            :: ncfile*99

    ! BLAZE variables
    type(TYPE_BLAZE)    :: BLAZE
    type(TYPE_SIMFIRE)  :: SIMFIRE

    ! 13C
    type(c13o2_flux)  :: c13o2flux
    type(c13o2_pool)  :: c13o2pools, sum_c13o2pools
    type(c13o2_luc)   :: c13o2luc
    real(r_2), dimension(:,:), allocatable :: casasave
    real(r_2), dimension(:,:), allocatable :: lucsave
    ! I/O
    logical :: first_casa_write
    integer :: c13o2_outfile_id
    character(len=40), dimension(c13o2_nvars_output) :: c13o2_vars
    integer,           dimension(c13o2_nvars_output) :: c13o2_var_ids
    ! delta-13C of atmospheric CO2
    integer            :: iunit, ios
    real               :: iyear
    integer            :: c13o2_atm_syear, c13o2_atm_eyear
    character(len=100) :: header

    ! declare vars for switches (default .FALSE.) etc declared thru namelist
    logical, save           :: &
         vegparmnew    = .false., & ! using new format input file (BP dec 2007)
         spinup        = .false., & ! model spinup to soil state equilibrium?
         spinConv      = .false., & ! has spinup converged?
         spincasainput = .false., & ! TRUE: SAVE input req'd to spin CASA-CNP;
                                    ! FALSE: READ input to spin CASA-CNP
         spincasa      = .false., & ! TRUE: CASA-CNP Will spin mloop,
                                    ! FALSE: no spin up
         l_casacnp     = .false., & ! using CASA-CNP with CABLE
         l_laiFeedbk   = .false., & ! using prognostic LAI
         l_vcmaxFeedbk = .false., & ! using prognostic Vcmax
         CASAONLY      = .false., & ! ONLY Run CASA-CNP
         CALL1         = .true.

         real              :: &
         delsoilM, & ! allowed variation in soil moisture for spin up
         delsoilT            ! allowed variation in soil temperature for spin up

    ! temporary storage for soil moisture/temp. in spin up mode
    real, allocatable, dimension(:,:) :: &
         soilMtemp, &
         soilTtemp

    ! MPI:
    type(met_type)           :: imet  ! read ahead met input variables
    type(veg_parameter_type) :: iveg  ! MPI read ahead vegetation parameters
    logical :: loop_exit     ! MPI: exit flag for bcast to workers
    integer :: iktau    ! read ahead index of time step = 1 ..  kend
    integer :: oktau    ! ktau = 1 ..  kend for output
    integer :: icomm ! separate dupes of MPI communicator for send and recv
    integer :: ocomm ! separate dupes of MPI communicator for send and recv
    integer :: ierr
    integer :: rank, off, cnt

    ! Vars for standard for quasi-bitwise reproducability b/n runs
    ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
    character(len=30), parameter :: &
         Ftrunk_sumbal  = ".trunk_sumbal", &
         Fnew_sumbal    = "new_sumbal"

    real(r_2), save :: &
         trunk_sumbal = 0.0, & !
         new_sumbal   = 0.0, &
         new_sumfpn   = 0.0, &
         new_sumfe    = 0.0

    integer :: count_bal = 0
    integer :: nkend=0
    integer :: ioerror=0

    logical :: casa_time
    logical :: liseod, liseoy ! is end of day, is end of year

    ! switches etc defined thru namelist (by default cable.nml)
    namelist /cablenml/ &
         filename, & ! TYPE, containing input filenames
         vegparmnew, & ! use new soil param. method
         soilparmnew, & ! use new soil param. method
         calcsoilalbedo, & ! vars intro for Ticket #27
         spinup, & ! spinup model (soil) to steady state
         delsoilM,delsoilT, & !
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
         cable_user           ! additional USER switches

    integer :: kk
    integer :: lalloc
    integer, parameter :: mloop = 30 ! CASA-CNP PreSpinup loops
    real :: etime, etimelast

    ! command line arguments
    integer :: narg, len1, len2
    character(len=500) :: arg1
    character(len=200) :: arg2

    ! end header

    etimelast = 0.0

    ! Open, read and close the namelist file.
    open(10, file=cable_namelist, status="old", action="read")
    read(10, nml=cablenml)   !where nml=cable defined above
    close(10)

    ! Open, read and close the consistency check file.
    ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
    if (cable_user%consistency_check) then
       open(11, file=ftrunk_sumbal, status='old', action='read', iostat=ioerror)
       if (ioerror==0) then
          read(11,*) trunk_sumbal  ! written by previous trunk version
       end if
       close(11)
    end if

    ! Open log file:
    open(logn, file=filename%log)

    call report_version_no(logn)

    ! if (iargc() > 0) then
    !    call getarg(1, filename%met)
    !    call getarg(2, casafile%cnpipool)
    ! end if
    narg = command_argument_count()
    if (narg > 0) then
       call get_command_argument(1, arg1, len1)
       filename%met = arg1(1:len1)
       call get_command_argument(2, arg2, len2)
       casafile%cnpipool = arg2(1:len2)
    end if

    ! INITIALISATION depending on nml settings
    if (trim(cable_user%MetType) == 'gswp') THEN
       if ((cable_user%YearStart == 0) .and. (ncciy > 0)) then
          cable_user%YearStart = ncciy
          cable_user%YearEnd   = ncciy
       else if ((cable_user%YearStart == 0) .and. (ncciy == 0)) then
          write(*,*) 'undefined start year for gswp met: '
          write(*,*) 'enter value for ncciy or'
          write(*,*) '(cable_user%YearStart and  cable_user%YearEnd) in cable.nml'
          write(logn,*) 'undefined start year for gswp met:'
          write(logn,*) 'enter value for ncciy or'
          write(logn,*) '(cable_user%YearStart and  cable_user%YearEnd) in cable.nml'
          call MPI_Abort(comm, 8, ierr)
       end if
    end if

    CurYear = cable_user%YearStart

    if (icycle >= 11) then
       icycle                     = icycle - 10
       casaonly                   = .true.
       cable_user%casa_dump_read  = .true.
       cable_user%casa_dump_write = .false.
    else if (icycle == 0) then
       cable_user%casa_dump_read  = .false.
       cable_user%call_pop        = .false.
       cable_user%call_blaze      = .false.
    end if

    !! vh_js !!
    if (icycle > 0) then
       l_casacnp = .true.
    else
       l_casacnp = .false.
    end if

    !! vh_js !! suggest LALLOC should ulitmately be a switch in the .nml file
    if (cable_user%call_pop) then
       lalloc = 3 ! for use with POP: makes use of pipe model to partition between stem and leaf
    else
       lalloc = 0 ! default
    end if

    if (trim(cable_user%MetType) == 'gpgs') then
       leaps = .true.
       cable_user%MetType = 'gswp'
    else if (trim(cable_user%MetType) == 'bios') then
       leaps = .true.
    end if

    cable_runtime%offline = .true.

    ! associate pointers used locally with global definitions
    call point2constants(C)

    if (l_casacnp  .and. ((icycle == 0) .or. (icycle > 3))) then
       write(*,*) 'icycle must be 1 to 3 when using casaCNP'
       call MPI_Abort(comm, 9, ierr)
    end if
    if ((l_laifeedbk .or. l_vcmaxfeedbk) .and. (.not. l_casacnp)) then
       write(*,*) 'casaCNP required to get prognostic LAI or Vcmax'
       call MPI_Abort(comm, 10, ierr)
    end if
    if (l_vcmaxfeedbk .and. (icycle < 2)) then
       write(*,*) 'icycle must be 2 to 3 to get prognostic Vcmax'
       call MPI_Abort(comm, 11, ierr)
    end if
    if ((icycle > 0) .and. (.not. soilparmnew) ) then
       write(*,*) 'casaCNP must use new soil parameters'
       call MPI_Abort(comm, 12, ierr)
    end if

    ! casa time count
    ctime = 0
    first_casa_write = .true.

    ! Initialise settings depending on met dataset

    ! Open met data and get site information from netcdf file. (NON-GSWP ONLY!)
    ! This retrieves time step size, number of timesteps, starting date,
    ! latitudes, longitudes, number of sites.
    if (trim(cable_user%MetType) .ne. "gswp" .and. &
        trim(cable_user%MetType) .ne. "gpgs" .and. &
        trim(cable_user%MetType) .ne. "plume" .and. &
        trim(cable_user%MetType) .ne. "bios" .and. &
        trim(cable_user%MetType) .ne. "cru") then
       call open_met_file( dels, koffset, kend, spinup, C%TFRZ )
       if ((koffset /= 0) .and. cable_user%CALL_POP) then
          write(*,*) "When using POP, episode must start on Jan 1st!"
          call MPI_Abort(comm, 13, ierr)
       end if
    end if

    !MCJK - check if need in mpimaster and mpiworker
    ! Read gm lookup table
    if (cable_user%explicit_gm .and. (len_trim(cable_user%gm_LUT_file) > 1)) then
        write(*,*) 'Reading gm LUT file'
        call read_gm_LUT(trim(cable_user%gm_LUT_file), LUT_VcmaxJmax, &
             LUT_gm, LUT_Vcmax, LUT_Rd)
    end if
    !MCJK - check if need in mpimaster and mpiworker

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

    ! Tell the workers if we're leaping
    call MPI_Bcast(leaps, 1, MPI_LOGICAL, 0, comm, ierr)

    ! outer loop - spinup loop no. ktau_tot
    ktau_tot = 0
    ktau     = 0
    SPINLOOP: do

       ! Bios initialisation
       if (trim(cable_user%MetType) == "bios") then
          call cpu_time(etime)
          call cable_bios_init(dels, curyear, met, kend, ktauday)
          koffset = 0
          leaps   = .true.
          write(str1,'(i4)') curyear
          str1 = adjustl(str1)
          timeunits = "seconds since " // trim(str1) // "-01-01 00:00:00"
          calendar = "standard"
          casa_timeunits = "days since " // trim(str1) // "-01-01 00:00:00"
       end if

       write(*,*) "CABLE_USER%YearStart,  CABLE_USER%YearEnd", CABLE_USER%YearStart,  CABLE_USER%YearEnd

       ! Loop through simulation years
       YEARLOOP: do YYYY=cable_user%YearStart, cable_user%YearEnd

          CurYear = YYYY
          if (leaps .and. is_leapyear(YYYY)) then
             LOY = 366
          else
             LOY = 365
          end if

          ! CLN From here extract CALL1 and put LOY computation
          !     to end of CALL1 block
          if (trim(cable_user%MetType) == 'gswp') then
             ncciy = CurYear
             write(*,*) 'Looking for global offline run info.'
             call preparefiles(ncciy)
             call open_met_file( dels, koffset, kend, spinup, c%tfrz )
          else if (trim(cable_user%MetType) == 'plume') then
             if (CALL1) then
                call cpu_time(etime)
                call plume_mip_init( plume )
                dels    = plume%dt
                koffset = 0
                leaps   = plume%LeapYears
                write(str1,'(i4)') CurYear
                str1 = adjustl(str1)
                timeunits = "seconds since "//trim(str1)//"-01-01 00:00:00"
                if (leaps) then
                   calendar = "standard"
                else
                   calendar = "noleap"
                end if
                casa_timeunits = "days since "//trim(str1)//"-01-01 00:00:00"
             end if
             if (.not. PLUME%LeapYears) LOY = 365
             kend = nint(24.0 * 3600.0 / dels) * LOY
          else if (trim(cable_user%MetType) == 'bios') then
             kend = nint(24.0 * 3600.0 / dels) * LOY
          else if (trim(cable_user%MetType) == 'cru') then
             if (CALL1) then
                call cpu_time(etime)
                call cru_init( cru )
                dels         = cru%dtsecs
                koffset      = 0
                if (CRU%MetVersion == "VERIFY_2021") then
                   leaps = .true.
                   calendar = "standard"
                   if ( IS_LEAPYEAR(CurYear) ) then
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
                write(str1,'(i4)') CurYear
                str1 = adjustl(str1)
                timeunits = "seconds since " // trim(str1) // "-01-01 00:00:00"
                calendar  = "standard"
                casa_timeunits = "days since " // trim(str1) // "-01-01 00:00:00"
             end if
             kend = nint(24.0 * 3600.0 / dels) * LOY
             ! no ( trim(cable_user%MetType) .eq. 'site' ) in MPI
          end if ! cable_user%MetType
          ! CLN To here extract CALL1 and put LOY computation to end of CALL1 block

          ! some things (e.g. CASA-CNP) only need to be done once per day
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
                allocate(casasave(c13o2pools%ntile, c13o2pools%npools))
                if (cable_user%popluc) &
                     allocate(lucsave(c13o2luc%nland, c13o2luc%npools))
             end if

             !par disabled as blaze init moved below
             ! ! Abort, if an error occurred during BLAZE/SIMFIRE init !CLN check again
             ! IF (BLAZE%ERR) CALL MPI_Abort(comm,0,ierr)

             if (cable_user%POPLUC .and. &
                  (trim(cable_user%POPLUC_RunType) == 'static')) &
                  cable_user%POPLUC = .false.

             ! Having read the default parameters, if this is a bios run we will now
             ! overwrite the subset of them required for bios.
             if (trim(cable_user%MetType) .eq. 'bios') &
                  call cable_bios_load_params(soil)

             ! Open output file
             if (.not. CASAONLY) then
                if (trim(filename%out) == '') then
                   if (cable_user%YEARSTART > 0) then
                      write(dum, FMT="(I4,'_',I4)") cable_user%YEARSTART, &
                           cable_user%YEAREND
                      filename%out = trim(filename%path) // '/' // &
                           trim(cable_user%RunIden) // '_' // &
                           trim(dum) // '_cable_out.nc'
                   else
                      filename%out = trim(filename%path) // '/' // &
                           trim(cable_user%RunIden) // '_cable_out.nc'
                   end if
                end if
                if (YYYY.eq.cable_user%YEARSTART) then
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

             spinConv = .FALSE. ! initialise spinup convergence variable
             if (.not. spinup) spinConv = .true.

             !TRUNK not in trunk
             if (cable_user%call_climate) then
                call alloc_cbm_var(climate, mp, ktauday)
                call zero_cbm_var(climate)
                call climate_init(climate)
                if (.not. cable_user%climate_fromzero) &
                     call READ_CLIMATE_RESTART_NC(climate)
             end if

             !MC - do we need that? casamet also set only if icycle>0
             ! if (trim(cable_user%MetType) .eq. 'cru') then
             !    casamet%glai = 1.0_r_2 ! initialise glai for use in cable_roughness
             !    where (veg%iveg(:) .ge. 14) casamet%glai = 0.0_r_2
             ! end if
             if ((trim(cable_user%MetType) == 'cru') .and. (icycle > 0)) then
                if (all(real(casamet%glai(:)) == 0.)) then
                   write(*,*) 'W A R N I N G : deleted setting casamet%glai for CRU.'
                   call MPI_Abort(comm, 50, ierr)
                end if
             end if

             if (trim(cable_user%MetType) == 'bios') &
                  call cable_bios_load_climate_params(climate)

             ! MPI: above was standard serial code
             !      now it's time to initialize the workers

             ! MPI: bcast to workers so that they do not need to open
             !      the met file themselves
             call MPI_Bcast(dels, 1, MPI_REAL, 0, comm, ierr)

          end if ! Call1

          call MPI_Bcast(kend, 1, MPI_INTEGER, 0, comm, ierr)

          if (CALL1) then
             ! MPI: need to know extents before creating datatypes
             call find_extents()

             ! MPI: calculate and broadcast landpoint decomposition to the workers
             call master_decomp(comm, mland)

             ! MPI: set up stuff for new irecv isend code that separates completion
             ! from posting of requests
             ! wnp is set in master_decomp above
             allocate(inp_req(wnp))
             allocate(inp_stats(MPI_STATUS_SIZE, wnp))
             allocate(recv_req(wnp))
             allocate(recv_stats(MPI_STATUS_SIZE, wnp))
             call MPI_Comm_dup(comm, icomm, ierr)
             call MPI_Comm_dup(comm, ocomm, ierr)

             ! MPI: data set in load_parameter is now scattered out to the workers
             CALL master_cable_params(comm, met, air, ssnow, veg, bgc, &
                  soil, canopy, rough, rad, sum_flux, bal)

             if (cable_user%call_climate) &
                  call master_climate_types(comm, climate)

             ! MPI: mvtype and mstype send out here instead of inside
             !      master_casa_params so that old CABLE carbon module
             !      can use them. (BP May 2013)
             CALL MPI_Bcast(mvtype, 1, MPI_INTEGER, 0, comm, ierr)
             CALL MPI_Bcast(mstype, 1, MPI_INTEGER, 0, comm, ierr)

             ! MPI: casa parameters scattered only if cnp module is active
             IF (icycle > 0) THEN
                ! MPI:
                CALL master_casa_params(comm, casabiome, casapool, casaflux, &
                     casamet, casabal, phen)
                ! 13C
                if (cable_user%c13o2) then
                   CALL master_c13o2_flux_params(comm, c13o2flux)
                   CALL master_c13o2_pool_params(comm, c13o2pools)
                   if (cable_user%popluc) then
                      call master_c13o2_luc_params(comm, c13o2luc)
                   end if
                end if

                ! Create and Send POP ini/restart data
                if (cable_user%call_pop) then
                   call master_pop_types(comm, pop)
                end if

                ! Fire init
                if (cable_user%call_blaze) then
                   call ini_blaze(mland, rad%latitude(landpt(:)%cstart), &
                        rad%longitude(landpt(:)%cstart), blaze)

                   !par blaze restart not required uses climate data
                   !create handles for restart-data
                   allocate(blaze_restart_ts(wnp))
                   allocate(blaze_out_ts(wnp))
                   allocate(blaze_in_ts(wnp))
                   call master_blaze_types(comm, wland, wnp, mland, blaze, &
                        blaze_restart_ts, blaze_in_ts, blaze_out_ts)
                   if (trim(blaze%burnt_area_src) == "SIMFIRE" ) then
                      !MCfire - mpiworker receives ktau_gl as tag
                      call master_send_input(comm, blaze_in_ts, ktau)
                      !CLN here we need to check for the SIMFIRE biome &
                      !    setting
                      call INI_SIMFIRE(mland ,SIMFIRE, &
                         climate%modis_igbp(landpt(:)%cstart))
                      WRITE(logn,*)"After ini_simf"
                      CALL FLUSH(logn)
                   end if ! simfire
                end if ! call blaze
             end if ! icycle

             ! MPI: allocate read ahead buffers for input met and veg data
             call alloc_cbm_var(imet, mp)
             call alloc_cbm_var(iveg, mp)
             call zero_cbm_var(imet)
             call zero_cbm_var(iveg)

             ! MPI: create inp_t types to scatter input data to the workers
             ! at the start of every timestep
             !   CALL master_intypes (comm,met,veg)
             ! for read ahead use the new variables
             call master_intypes(comm, imet, iveg)

             ! MPI: create recv_t types to receive results from the workers
             ! at the end of every timestep
             call master_outtypes(comm, met, canopy, ssnow, rad, bal, air, soil, veg)
             if (cable_user%c13o2) call master_c13o2_flux_types(comm, c13o2flux)

             ! MPI: create type for receiving casa results
             ! only if cnp module is active
             if (icycle > 0) then
                call master_casa_types(comm, casabiome, casapool, casaflux, &
                     casamet, casabal, phen)
                ! 13C
                if (cable_user%c13o2) call master_c13o2_pool_types(comm, c13o2pools)

                if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                   ! 13C
                   call master_casa_dump_types(comm, casamet, casaflux, phen, climate, c13o2flux)
                end if
                if (cable_user%popluc) then
                   CALL master_casa_LUC_types(comm, casapool, casabal, casaflux)
                   ! 13C
                   if (cable_user%c13o2) then
                      call master_c13o2_luc_types(comm, c13o2luc)
                   end if
                end if
             end if ! icycle > 0

             ! MPI: create type to send restart data back to the master
             ! only if restart file is to be created
             if (output%restart) then
                call master_restart_types(comm, canopy, air, veg, ssnow)
             end if

             ! CALL master_sumcasa_types(comm, sum_casapool, sum_casaflux)
             if ((icycle > 0) .and. spincasa) then
                write(*,*) 'EXT spincasacnp enabled with mloop= ', mloop, dels, kstart, kend
                ! 13C
                call master_spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                     casaflux,casamet,casabal,phen,POP,climate, &
                     c13o2flux, c13o2pools, icomm, ocomm)
                if (cable_user%c13o2) call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                SPINconv = .false.
                CASAONLY = .true.
             else if (casaonly .and. (.not. spincasa) .and. cable_user%popluc) then
                ! 13C
                write(*,*) 'EXT CASAONLY_LUC'
                call master_CASAONLY_LUC(dels,kstart,kend,veg,casabiome,casapool, &
                     casaflux,casamet,casabal,phen,POP,climate, LUC_EXPT, POPLUC, &
                     c13o2flux, c13o2pools, c13o2luc, icomm, ocomm)
                if (cable_user%c13o2) then
                   call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                   if (cable_user%POPLUC) call c13o2_sanity_luc(popluc, c13o2luc)
                end if
                SPINconv = .false.
             end if

             ! MPI: mostly original serial code follows...
          end if ! CALL1

          ! globally (WRT code) accessible kend through USE cable_common_module
          ktau_gl   = 0
          kwidth_gl = int(dels)
          kend_gl   = kend
          knode_gl  = 0

          ! MPI: separate time step counters for reading and writing
          ! (ugly, I know)
          iktau = ktau_gl
          oktau = ktau_gl

          ! MPI: read ahead
          iktau = iktau + 1

          ! MPI: flip ktau_gl
          !tmp_kgl = ktau_gl
          ktau_gl = iktau

          if (spincasa .or. casaonly) then
             exit
          end if

          ! Read ahead: send input before workes start ktau loop
          if (.not. CASAONLY) then
             if (trim(cable_user%MetType) == 'plume') then
                call PLUME_MIP_GET_MET(PLUME, iMET, YYYY, iktau, kend, &
                     ((YYYY == cable_user%YearEnd) .and. (iktau == kend)))
             else if (trim(cable_user%MetType) == 'bios') then
                call  cable_bios_read_met(iMET, CurYear, iktau, dels)
             else if (trim(cable_user%MetType) == 'cru') then
                call CRU_GET_SUBDIURNAL_MET(CRU, imet, YYYY, iktau, kend)
                ! MC - Added two lines defining iveg and lai in input veg: iveg
                !      because this is done in get_met_data but not in cru_get_subdiurnal_met
                !MCcheck
                iveg%iveg = veg%iveg
                iveg%vlai = veg%vlai
             else
                call get_met_data(spinup, spinConv, imet, &
                     rad, iveg, dels, C%TFRZ, iktau+koffset, &
                     kstart+koffset)
             end if
          end if

          if (trim(cable_user%MetType) == '') then
             CurYear = met%year(1)
             if (leaps .and. IS_LEAPYEAR(CurYear)) then
                LOY = 366
             else
                LOY = 365
             end if
          end if
          imet%ofsd = imet%fsd(:,1) + imet%fsd(:,2)
          canopy%oldcansto = canopy%cansto
          ! Zero out lai where there is no vegetation acc. to veg. index
          where (iveg%iveg(:) >= 14) iveg%vlai = 0.

          ! Read ahead: send input before workes start ktau loop
          if (.not. CASAONLY) then
             ! MPI: scatter input data to the workers
             call master_send_input(icomm, inp_ts, iktau)
             !MC - The else statement that sends casa_dump data is not commented out in the trunk.
             !     And in mpiworker there is:
             !         ELSE IF ( IS_CASA_TIME("dread", yyyy, ktau, kstart, koffset, ktauday, wlogn) ) then
          ! ELSE
             ! CALL master_send_input(icomm, casa_dump_ts, iktau)
             ! CALL MPI_Waitall(wnp, inp_req, inp_stats, ierr)
          end if

          ! 13C
          if (cable_user%c13o2) then
             if ((CurYear < c13o2_atm_syear) .or. (CurYear > c13o2_atm_eyear)) then
                write(*,*) 'Current year ', CurYear, &
                     'not in atmospheric delta-13C (min/max): ', &
                     c13o2_atm_syear, c13o2_atm_eyear
                call MPI_Abort(comm, 14, ierr)
             end if
             c13o2flux%ca = (c13o2_delta_atm(CurYear) + 1.0_r_2) * &
                  real(imet%ca, r_2) ! * vpdbc13 / vpdbc13
             ! broadcast
             do rank=1, wnp
                off = wland(rank)%patch0
                cnt = wland(rank)%npatch
                call MPI_Send(c13o2flux%ca(off), cnt, MPI_Double_precision, &
                     rank, 0, icomm, ierr)
             end do
          end if

          ! IF (.NOT.spincasa) THEN
          ! time step loop over ktau
          KTAULOOP: do ktau=kstart, kend-1
             ! globally (WRT code) accessible kend through USE cable_common_module
             iktau = iktau + 1
             oktau = oktau + 1
             ktau_tot = ktau_tot + 1
             ktau_gl  = iktau

             met%year = imet%year
             met%doy  = imet%doy

             ! some things (e.g. CASA-CNP) only need to be done once per day
             idoy = int( mod(real(ceiling((real(ktau+koffset)) / &
                  real(ktauday))), real(loy)) )
             if (idoy .eq. 0) idoy = LOY

             ! needed for CASA-CNP
             !MC - diff to serial
             nyear = int((kend-kstart+1) / (loy*ktauday))

             ! Get met data and LAI, set time variables.
             ! Rainfall input may be augmented for spinup purposes:
             !          met%ofsd = met%fsd(:,1) + met%fsd(:,2)
             if (trim(cable_user%MetType) == 'plume') then
                call plume_mip_get_met(plume, imet, yyyy, iktau, kend, &
                     (yyyy == cable_user%YearEnd) .and. (iktau == kend))
             else if (trim(cable_user%MetType) == 'bios') then
                if ((.not. CASAONLY) .or. (CASAONLY .and. CALL1)) then
                   call cable_bios_read_met(imet, CurYear, iktau, dels)
                end if
             else if (trim(cable_user%MetType) == 'cru') then
                call cru_get_subdiurnal_met(cru, imet, YYYY, iktau, kend)
                ! MC - Added two lines defining iveg and lai in input veg: iveg
                !      because this is done in get_met_data but not in cru_get_subdiurnal_met
                iveg%iveg = veg%iveg ! LAI and veg class not in CRU
                iveg%vlai = veg%vlai
             else
                call get_met_data(spinup, spinconv, imet, &
                     rad, iveg, dels, c%tfrz, iktau+koffset, &
                     kstart+koffset)
             end if

             if (trim(cable_user%MetType) == '') then
                CurYear = imet%year(1)
                if (leaps .and. IS_LEAPYEAR(CurYear)) then
                   LOY = 366
                else
                   LOY = 365
                end if
             end if
             imet%ofsd = imet%fsd(:,1) + imet%fsd(:,2)
             canopy%oldcansto = canopy%cansto
             ! Zero out lai where there is no vegetation acc. to veg. index
             where (iveg%iveg(:) >= 14) iveg%vlai = 0.

             ! At first time step of year, set tile area according to updated LU areas
             if (ktau == 1) then
                if (icycle > 1) &
                     call casa_cnpflux(casaflux, casapool, casabal, .true.)
                if (cable_user%POPLUC) then
                   call POPLUC_set_patchfrac(POPLUC, LUC_EXPT)
                end if
             end if

             casa_time = IS_CASA_TIME("write", yyyy, oktau, kstart, koffset, &
                  ktauday, logn)
             liseod = mod((oktau-kstart+1+koffset), ktauday) == 0
             liseoy = mod((oktau-kstart+1+koffset) / ktauday, LOY) == 0

             if (.not. CASAONLY) then

                if (icycle > 0) then
                   ! receive casa update from worker
                   if (liseod) then
                      call master_receive(ocomm, oktau, casa_ts)
                      ! 13C
                      if (cable_user%c13o2) &
                           call master_receive(ocomm, oktau, c13o2_pool_ts)

                      ctime = ctime + 1
                      ! update time-aggregates of casa pools and fluxes
                      count_sum_casa = count_sum_casa + 1
                      call update_sum_casa(sum_casapool, sum_casaflux, &
                           casapool, casaflux, .true., casa_time, &
                           count_sum_casa)
                      ! 13C
                      if (cable_user%c13o2) &
                           call c13o2_update_sum_pools(sum_c13o2pools, &
                           c13o2pools, .true., casa_time, count_sum_casa)

                      if (cable_user%call_blaze) then
                         call master_receive(ocomm, oktau, blaze_out_ts)
                         BLAZE%time =  BLAZE%time + 86400
                         call write_blaze_output_nc(BLAZE, &
                              (ktau == kend) .AND. (YYYY == cable_user%YearEnd))
                      end if
                   end if ! end of day

                   ! receive casa dump requirements from worker
                   if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                      if (((.not.spinup) .or. (spinup .and. spinConv)) .and. &
                           is_casa_time("dwrit", yyyy, oktau, kstart, koffset, ktauday, logn)) then
                         call master_receive(ocomm, oktau, casa_dump_ts)
                      end if
                   end if
                end if ! icycle>0
                ! MPI: receive this time step's results from the workers
                call master_receive(ocomm, oktau, recv_ts)
                if (cable_user%c13o2) call master_receive(ocomm, oktau, c13o2_flux_ts)
                ! MPI: scatter input data to the workers
                call master_send_input(icomm, inp_ts, iktau)
                ! 13C
                if (cable_user%c13o2) then
                   ! already checked that CurYear is in input file
                   c13o2flux%ca = (c13o2_delta_atm(CurYear) + 1.0_r_2) * real(imet%ca,r_2)
                   do rank=1, wnp
                      off = wland(rank)%patch0
                      cnt = wland(rank)%npatch
                      call MPI_Send(c13o2flux%ca(off), cnt, MPI_DOUBLE_PRECISION, rank, 0, icomm, ierr)
                   end do
                end if
                if (((.not. spinup) .or. (spinup .and. spinConv)) .and. liseod) then
                   if (cable_user%casa_dump_write) then
                      !cln check for leap year
                      write(cyear, fmt="(I4)") curyear + &
                           int((oktau-kstart+koffset) / (loy*ktauday))
                      ncfile = trim(casafile%c2cdumppath) // 'c2c_' // cyear // '_dump.nc'
                      ! 13C
                      call write_casa_dump(ncfile, casamet, casaflux, phen, &
                           climate, c13o2flux, idoy, kend/ktauday)
                   end if
                end if
             else if (is_casa_time("dread", yyyy, iktau, kstart, koffset, &
                  ktauday, logn)) then ! .not. casaonly
                if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                   call master_send_input(icomm, casa_dump_ts, iktau)
                end if

             end if ! .not. CASAONLY

             met%ofsd = met%fsd(:,1) + met%fsd(:,2)
             canopy%oldcansto = canopy%cansto

             ! Zero out lai where there is no vegetation acc. to veg. index
             where (iveg%iveg(:) .ge. 14) iveg%vlai = 0.
             ! Write time step's output to file if either: we're not spinning up
             ! or we're spinning up and the spinup has converged:
             ! MPI: TODO: pull mass and energy balance calculation from write_output
             ! and refactor into worker code

             ktau_gl = oktau
             if ((.not. spinup) .or. (spinup .and. spinConv)) then
                if (icycle > 0) then
                   if (casa_time) then
                      !ctime = ctime + 1
                      ! count_sum_casa = count_sum_casa + 1
                      ! call update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                      !      .true., casa_time, count_sum_casa)
                      ! ! 13C
                      ! if (cable_user%c13o2) &
                      !      call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, .true., casa_time, count_sum_casa)
                      ! count_sum_casa = 0

                      call write_casa_output_nc(veg, casamet, sum_casapool, &
                           casabal, sum_casaflux, CASAONLY, ctime, &
                           (ktau == kend) .and. (yyyy == cable_user%YearEnd))
                      ! 13C
                      if (cable_user%c13o2) then
                         if (first_casa_write) then
                            call c13o2_create_output(casamet, sum_c13o2pools, &
                                 c13o2_outfile_id, c13o2_vars, c13o2_var_ids)
                            first_casa_write = .false.
                         end if
                         call c13o2_write_output(c13o2_outfile_id, c13o2_vars, &
                              c13o2_var_ids, ctime, sum_c13o2pools)
                      end if
                      if (casa_time) then
                         count_sum_casa = 0
                         call zero_sum_casa(sum_casapool, sum_casaflux)
                         ! 13C
                         if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
                      end if
                   end if
                end if

                if ((.not. CASAONLY) .and. spinConv) then
                   if (trim(cable_user%mettype) == 'plume' &
                       .or. trim(cable_user%mettype) == 'cru' &
                       .or. trim(cable_user%mettype) == 'bios' &
                       .or. trim(cable_user%mettype) == 'gswp') then
                      ! 13C
                      call write_output(dels, ktau_tot, met, canopy, casaflux, &
                           casapool, casamet, ssnow, rad, bal, air, soil, veg, &
                           c%sboltz, c%emleaf, c%emsoil, c13o2pools, c13o2flux)
                   else
                      ! 13C
                      call write_output(dels, ktau, met, canopy, casaflux, &
                           casapool, casamet, ssnow, rad, bal, air, soil, veg, &
                           c%sboltz, c%emleaf, c%emsoil, c13o2pools, c13o2flux)
                   end if
                end if

             end if ! .not. spinup .or. (spinup.and.spinConv)

             !---------------------------------------------------------------------!
             ! Check this run against standard for quasi-bitwise reproducability   !
             ! Check triggered by cable_user%consistency_check=.TRUE. in cable.nml !
             !---------------------------------------------------------------------!
             if (cable_user%consistency_check) then
                count_bal = count_bal + 1
                new_sumbal = new_sumbal + sum(bal%wbal) / mp + &
                     sum(bal%ebal) / mp
                new_sumfpn = new_sumfpn + sum(canopy%fpn) / mp
                new_sumfe  = new_sumfe  + sum(canopy%fe)  / mp
                ! if (ktau == kend-1) PRINT*, "time-space-averaged energy & water balances"
                ! if (ktau == kend-1) PRINT*,"Ebal_tot[Wm-2], Wbal_tot[mm]", &
                !      sum(bal%ebal_tot)/mp/count_bal, sum(bal%wbal_tot)/mp/count_bal
                ! if (ktau == kend-1) PRINT*, "time-space-averaged latent heat and net photosynthesis"
                ! if (ktau == kend-1) PRINT*, "sum_fe[Wm-2], sum_fpn[umol/m2/s]", &
                !      new_sumfe/count_bal, new_sumfpn/count_bal

                ! check for Nans in biophysical outputs and abort if there are any
                if (any(ne(canopy%fe, canopy%fe))) then
                   do kk=1, mp
                      if (ne(canopy%fe(kk), canopy%fe(kk))) then
                         write(*,*) 'Nan in evap flux,', kk, patch(kk)%latitude, patch(kk)%longitude
                         write(*,*) 'fe nan', kk, ktau,met%qv(kk), met%precip(kk),met%precip_sn(kk), &
                              met%fld(kk), met%fsd(kk,:), met%tk(kk), met%ua(kk), &
                              ssnow%potev(kk), met%pmb(kk), &
                              canopy%ga(kk), ssnow%tgg(kk,:), canopy%fwsoil(kk), &
                              rad%fvlai(kk,:) ,  rad%fvlai(kk,1), &
                              rad%fvlai(kk,2), canopy%vlaiw(kk)
                         call MPI_Abort(comm, 15, ierr)
                      end if
                   end do
                end if

                if (ktau==(kend-1)) then
                   nkend = nkend+1
                   if( abs(new_sumbal-trunk_sumbal) < 1.e-7) then
                      write(*,*) ""
                      write(*,*) "NB. Offline-parallel runs spinup cycles:", nkend
                      write(*,*) "Internal check shows this version reproduces the trunk sumbal"
                   else
                      write(*,*) ""
                      write(*,*) "NB. Offline-parallel runs spinup cycles:", nkend
                      write(*,*) "Internal check shows in this version new_sumbal != trunk sumbal"
                      write(*,*) "The difference is: ", new_sumbal - trunk_sumbal
                      write(*,*) "Writing new_sumbal to the file:", trim(Fnew_sumbal)
                   end if
                end if
             end if ! cable_user%consistency_check

             CALL1 = .false.

          end do KTAULOOP ! END Do loop over timestep ktau

          CALL1 = .false.

          ! MPI: read ahead tail to receive (last step and write)
          met%year = imet%year
          met%doy  = imet%doy
          oktau    = oktau + 1
          ktau_tot = ktau_tot + 1
          ktau_gl  = oktau

          casa_time = IS_CASA_TIME("write", yyyy, ktau, kstart, koffset, ktauday, logn)
          if (.not. CASAONLY) then
             if (icycle > 0) then
                call master_receive(ocomm, oktau, casa_ts)
                ! 13C
                if (cable_user%c13o2) call master_receive(ocomm, oktau, c13o2_pool_ts)
                if (cable_user%call_blaze) then
                  !par recv blaze_out_ts
                  call master_receive(ocomm, oktau, blaze_out_ts)
                  BLAZE%time =  BLAZE%time + 86400
                  call write_blaze_output_nc(BLAZE, &
                       (ktau == kend) .and. (YYYY == cable_user%YearEnd))
                end if
                if ( ((.not. spinup) .or. (spinup .and. spinConv)) .and. &
                     IS_CASA_TIME("dwrit", yyyy, oktau, kstart, &
                     koffset, ktauday, logn) ) then
                   if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                      call master_receive(ocomm, oktau, casa_dump_ts)
                   end if
                end if
             end if ! icycle > 0
             call master_receive(ocomm, oktau, recv_ts)
             if (cable_user%c13o2) call master_receive(ocomm, oktau, c13o2_flux_ts)
          end if ! .not. CASAONLY

          met%ofsd = met%fsd(:,1) + met%fsd(:,2)
          canopy%oldcansto = canopy%cansto

          if ( trim(cable_user%MetType) .eq. "gswp" ) &
               call close_met_file()

          if ((icycle > 0) .and. cable_user%CALL_POP) then

             if (cable_user%POPLUC) then
                ! master receives casa updates required for LUC calculations here
                call master_receive(ocomm, 0, casa_LUC_ts)
                ! 13C
                if (cable_user%c13o2) then
                   call master_receive(ocomm, 0, c13o2_flux_ts)
                   call master_receive(ocomm, 0, c13o2_pool_ts)
                   call master_receive(ocomm, 0, c13o2_luc_ts)
                end if
                ! Dynamic LUC
                ! 13C
                call LUCdriver(casabiome, casapool, casaflux, &
                     POP, LUC_EXPT, POPLUC, veg, c13o2pools)
                ! transfer POP updates to workers
                off = 1
                do rank = 1, wnp
                   if ( rank .gt. 1 ) off = off + wland(rank-1)%npop_iwood
                   cnt = wland(rank)%npop_iwood
                   call MPI_Send( POP%pop_grid(off), cnt, pop_ts, rank, 0, icomm, ierr )
                end do
             end if ! POPLUC

             ! one annual time-step of POP (worker calls POP here)
             call master_receive_pop(POP, ocomm)

             if (cable_user%POPLUC) then
                ! Dynamic LUC: update casa pools according to LUC transitions
                ! 13C
                if (cable_user%c13o2) call c13o2_save_luc(casapool, popluc, &
                     casasave, lucsave)
                call POP_LUC_CASA_transfer(POPLUC, POP, LUC_EXPT, casapool, &
                     casabal, casaflux, ktauday)
                ! 13C
                if (cable_user%c13o2) then
                   call c13o2_update_luc(casasave, lucsave, popluc, &
                        luc_expt%prim_only, c13o2pools, c13o2luc)
                   call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                   call c13o2_sanity_luc(popluc, c13o2luc)
                end if

                ! Dynamic LUC: write output
                if (output%grid(1:3) == 'lan') then
                   call WRITE_LUC_OUTPUT_NC(POPLUC, YYYY, &
                        (YYYY == cable_user%YearEnd))
                else
                   call WRITE_LUC_OUTPUT_GRID_NC(POPLUC, YYYY, &
                        (YYYY == cable_user%YearEnd))
                end if

                ! send updates for CASA pools, resulting from LUC
                call master_send_input(icomm, casa_LUC_ts, nyear)
                ! 13C
                if (cable_user%c13o2) then
                   call master_send_input(icomm, c13o2_flux_ts, nyear)
                   call master_send_input(icomm, c13o2_pool_ts, nyear)
                   call master_send_input(icomm, c13o2_luc_ts, nyear)
                end if
             end if ! POPLUC

          end if ! icycle>0 .and. cable_user%CALL_POP

          if (icycle > 0) then
             ctime = ctime + 1
             ! update time-aggregates of casa pools and fluxes
             count_sum_casa = count_sum_casa + 1
             call update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                  .true., casa_time, count_sum_casa)
             ! 13C
             if (cable_user%c13o2) &
                  call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, .true., casa_time, count_sum_casa)
          end if

          ! WRITE OUTPUT
          if ((.not. spinup).or.(spinup .and. spinConv)) then
             if (icycle > 0) then
                ! ctime = ctime + 1
                !TRUNK no if but call write_casa in any case
                ! if ( is_casa_time("write", yyyy, ktau, kstart, koffset, ktauday, logn) ) then
                if (casa_time) then
                   ! count_sum_casa = count_sum_casa + 1
                   ! call update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
                   !      .true., .true., count_sum_casa)
                   ! ! 13C
                   ! if (cable_user%c13o2) &
                   !      call c13o2_update_sum_pools(sum_c13o2pools, c13o2pools, .true., .true., count_sum_casa)
                   ! call write_casa_output_nc(veg, casamet, casapool, casabal, casaflux, &
                   !      CASAONLY, ctime, (ktau==kend) .and. (YYYY==cable_user%YearEnd))
                   call write_casa_output_nc(veg, casamet, sum_casapool, &
                        casabal, sum_casaflux, CASAONLY, ctime, &
                        (ktau == kend) .and. (YYYY == cable_user%YearEnd))
                   ! 13C
                   if (cable_user%c13o2) then
                      call c13o2_write_output(c13o2_outfile_id, c13o2_vars, &
                           c13o2_var_ids, ctime, sum_c13o2pools)
                      if (YYYY == cable_user%YearEnd) &
                           call c13o2_close_output(c13o2_outfile_id)
                   end if
                   count_sum_casa = 0
                   call zero_sum_casa(sum_casapool, sum_casaflux)
                   ! 13C
                   if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
                end if
                if ( cable_user%CALL_POP ) then
                   ! CALL master_receive_pop(POP, ocomm)
                   ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
                   if (trim(cable_user%POP_out) == 'epi') then
                      call POP_IO(pop, casamet, CurYear, 'WRITE_EPI', &
                           (CurYear == cable_user%YearEnd))
                   end if
                end if
             end if ! icycle > 0

             if ( ((.not. spinup) .or. (spinup .and. spinConv)) .and. &
                  (mod((ktau-kstart+1), ktauday) == 0) ) then
                ! ctime = ctime + 86400  ! update casa time
                if (cable_user%CASA_DUMP_WRITE)  then
                   !CLN CHECK FOR LEAP YEAR
                   write(CYEAR,FMT="(I4)") CurYear + &
                        int((ktau-kstart) / (LOY*ktauday))
                   ncfile = trim(casafile%c2cdumppath) // 'c2c_' // CYEAR // '_dump.nc'
                   ! 13C
                   call write_casa_dump(ncfile, casamet , casaflux, phen, &
                        climate, c13o2flux, LOY, kend/ktauday)
                end if
             end if

             if ((.not. CASAONLY) .and. spinConv) then
                if (trim(cable_user%MetType) == 'plume' &
                     .or. trim(cable_user%MetType) == 'cru' &
                     .or. trim(cable_user%MetType) == 'bios' &
                     .or. trim(cable_user%MetType) == 'gswp') then
                   ! 13C
                   call write_output(dels, ktau_tot, met, canopy, casaflux, &
                        casapool, casamet, ssnow, rad, bal, air, soil, veg, &
                        C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux)
                else
                   ! 13C
                   call write_output(dels, ktau, met, canopy, casaflux, &
                        casapool, casamet, ssnow, rad, bal, air, soil, veg, &
                        C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux )
                end if
             end if

             if (cable_user%consistency_check) then
                count_bal = count_bal + 1
                new_sumbal = new_sumbal + sum(bal%wbal) / mp + &
                     sum(bal%ebal) / mp
                new_sumfpn = new_sumfpn + sum(canopy%fpn) / mp
                new_sumfe = new_sumfe + sum(canopy%fe) / mp
                if (ktau == kend) write(*,*) ""
                if (ktau == kend) write(*,*) "time-space-averaged energy & water balances"
                if (ktau == kend) write(*,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                     sum(bal%ebal_tot) / mp / count_bal, sum(bal%wbal_tot) / mp / count_bal
                if (ktau == kend) write(*,*) "time-space-averaged latent heat and net photosynthesis"
                if (ktau == kend) write(*,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s]", &
                     new_sumfe / count_bal, new_sumfpn / count_bal
                if (ktau == kend) write(logn,*) ""
                if (ktau == kend) write(logn,*) "time-space-averaged energy & water balances"
                if (ktau == kend) write(logn,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep] ", &
                     sum(bal%ebal_tot) / mp / count_bal, sum(bal%wbal_tot) / mp / count_bal
                if (ktau == kend) write(logn,*) "time-space-averaged latent heat and net photosynthesis"
                if (ktau == kend) write(logn,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s] ", &
                     new_sumfe / count_bal, new_sumfpn / count_bal
             end if

          end if ! (.NOT.spinup).OR.(spinup.AND.spinConv)

          ! set tile area according to updated LU areas
          if (cable_user%POPLUC) &
               call POPLUC_set_patchfrac(POPLUC, LUC_EXPT)

       end do YEARLOOP

       if (spincasa.or.casaonly) then
          exit
       end if

       ! see if spinup (if conducting one) has converged
       if (spinup .and. (.not. spinConv) .and. (.not. CASAONLY)) then
          ! Write to screen and log file:
          write(*,'(A18,I3,A24)') ' Spinning up: run ', int(ktau_tot / kend), &
               ' of data set complete...'
          write(logn,'(A18,I3,A24)') ' Spinning up: run ', int(ktau_tot / kend), &
               ' of data set complete...'

          ! IF not 1st run through whole dataset:
          if (int(ktau_tot / kend) > 1) then

             ! evaluate spinup
             if ( any(abs(ssnow%wb-soilMtemp) > delsoilM) .or. &
                  any(abs(ssnow%tgg-soilTtemp) > delsoilT) ) then
                ! No complete convergence yet
                maxdiff = maxloc(abs(ssnow%wb-soilMtemp))
                write(*,*) 'Example location of moisture non-convergence: ',maxdiff
                write(*,*) 'ssnow%wb : ', ssnow%wb(maxdiff(1), maxdiff(2))
                write(*,*) 'soilMtemp: ', soilMtemp(maxdiff(1), maxdiff(2))
                maxdiff = maxloc(abs(ssnow%tgg-soilTtemp))
                write(*,*) 'Example location of temperature non-convergence: ',maxdiff
                write(*,*) 'ssnow%tgg: ', ssnow%tgg(maxdiff(1), maxdiff(2))
                write(*,*) 'soilTtemp: ', soilTtemp(maxdiff(1), maxdiff(2))
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

          end if ! INT( ktau_tot/kend ) > 1

          if (YYYY == cable_user%YearEnd) then
             ! store soil moisture and temperature
             soilTtemp = ssnow%tgg
             soilMtemp = real(ssnow%wb)
          end if

          ! MPI:
          loop_exit = .false.

       else

          ! if not spinning up, or spin up has converged, exit:
          ! EXIT
          ! MPI:
          loop_exit = .true.

       end if ! spinup.AND..NOT.spinConv.AND. .NOT. CASAONLY

       ! MPI: let the workers know whether it's time to quit
       call MPI_Bcast(loop_exit, 1, MPI_LOGICAL, 0, comm, ierr)

       if (loop_exit) then
          exit
       end if

    end do SPINLOOP

    IF ((icycle > 0) .and. (.not. spincasa) .and. (.not. casaonly)) THEN
       ! MPI: gather casa results from all the workers
       CALL master_receive(ocomm, ktau_gl, casa_ts)
       ! 13C
       if (cable_user%c13o2) then
          call master_receive(ocomm, ktau_gl, c13o2_flux_ts)
          call master_receive(ocomm, ktau_gl, c13o2_pool_ts)
       end if
       if (cable_user%call_blaze) then
         !par recv blaze_out_ts
         call master_receive(ocomm, ktau_gl, blaze_out_ts)
       end if

       ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
       ! CALL casa_poolout( ktau, veg, soil, casabiome, &
       !     casapool, casaflux, casamet, casabal, phen )
       CALL casa_fluxout(nyear, veg, soil, casabal, casamet)
       ! CALL write_casa_restart_nc ( casamet, casapool, met, CASAONLY )
       ! CALL write_casa_restart_nc(casamet, casapool, casaflux, phen, CASAONLY)
       CALL write_casa_restart_nc(casabiome, casamet, casapool, casaflux, casabal, phen)
       ! 13C
       if (cable_user%c13o2) then
          call c13o2_write_restart_flux(casamet, c13o2flux)
          call c13o2_write_restart_pools(casamet, c13o2pools)
          if (cable_user%POPLUC) call c13o2_write_restart_luc(popluc, c13o2luc)
          ! While testing
          print*, 'not spincasa and not casaonly'
          call c13o2_print_delta_flux(c13o2flux)
          call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
          if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
       end if

       if (cable_user%CALL_POP .and. (POP%np > 0)) then
          if ( CASAONLY .or. cable_user%POP_fromZero &
               .or. (trim(cable_user%POP_out) == 'ini') ) then
             call POP_IO(pop, casamet, CurYear+1, 'WRITE_INI', .true.)
          else
             call POP_IO(pop, casamet, CurYear+1, 'WRITE_RST', .true.)
          end if
       end if
       if (cable_user%POPLUC .and. .not. CASAONLY ) then
          call WRITE_LUC_RESTART_NC(POPLUC)
       end if
    else if (icycle > 0) then
       ! 13C
       ! While testing
       if (cable_user%c13o2) then
          print*, 'spincasa or casaonly'
          call c13o2_print_delta_flux(c13o2flux)
          call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
          if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
       end if
    end if

    ! Write restart file if requested:
    if (output%restart .and. (.not. CASAONLY)) then
       call master_receive(comm, ktau_gl, restart_ts)
       call create_restart(logn, dels, ktau, soil, veg, ssnow, &
            canopy, rad, bgc, bal)

       if (cable_user%CALL_climate) then
          call master_receive(comm, ktau_gl, climate_ts)
          call WRITE_CLIMATE_RESTART_NC(climate)
       end if
    end if

    ! Close met data input file:
    if ( trim(cable_user%MetType) /= "gswp" .and. &
         trim(cable_user%MetType) /= "bios" .and. &
         trim(cable_user%MetType) /= "plume" .and. &
         trim(cable_user%MetType) /= "cru") call close_met_file()
    if (.not. CASAONLY) then
       ! Close output file and deallocate main variables:
       call close_output_file(bal)
      ! WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck
    end if

    if (cable_user%POPLUC) call close_luh2(LUC_EXPT)

    ! Close log file
    close(logn)

    call CPU_time(etime)
    write(*,*) 'Master End in ', etime, ' seconds.'
    ! MPI: cleanup
    call master_end(icycle, output%restart)

    return

  end subroutine mpidrv_master


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


  SUBROUTINE renameFiles(logn,inFile,ncciy,inName)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: logn
    CHARACTER(LEN=200), INTENT(INOUT) :: inFile
    INTEGER, INTENT(IN) :: ncciy
    CHARACTER(LEN=*),  INTENT(IN)    :: inName

    INTEGER :: nn
    INTEGER :: idummy

    nn = INDEX(inFile,'19')
    READ(inFile(nn:nn+3),'(i4)') idummy
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)

  END SUBROUTINE renameFiles

  !CLNSUBROUTINE prepareFiles(ncciy)
  !CLN  USE cable_IO_vars_module, ONLY: logn,gswpfile
  !CLN  IMPLICIT NONE
  !CLN  INTEGER, INTENT(IN) :: ncciy
  !CLN
  !CLN  WRITE(logn,*) 'CABLE offline global run using gswp forcing for ', ncciy
  !CLN  PRINT *,      'CABLE offline global run using gswp forcing for ', ncciy
  !CLN
  !CLN  CALL renameFiles(logn,gswpfile%rainf,16,ncciy,'rainf')
  !CLN  CALL renameFiles(logn,gswpfile%snowf,16,ncciy,'snowf')
  !CLN  CALL renameFiles(logn,gswpfile%LWdown,16,ncciy,'LWdown')
  !CLN  CALL renameFiles(logn,gswpfile%SWdown,16,ncciy,'SWdown')
  !CLN  CALL renameFiles(logn,gswpfile%PSurf,16,ncciy,'PSurf')
  !CLN  CALL renameFiles(logn,gswpfile%Qair,14,ncciy,'Qair')
  !CLN  CALL renameFiles(logn,gswpfile%Tair,14,ncciy,'Tair')
  !CLN  CALL renameFiles(logn,gswpfile%wind,15,ncciy,'wind')
  !CLN
  !CLNEND SUBROUTINE prepareFiles
  !CLN
  !CLN
  !CLNSUBROUTINE renameFiles(logn,inFile,nn,ncciy,inName)
  !CLN  IMPLICIT NONE
  !CLN  INTEGER, INTENT(IN) :: logn
  !CLN  INTEGER, INTENT(IN) :: nn
  !CLN  INTEGER, INTENT(IN) :: ncciy
  !CLN  CHARACTER(LEN=99), INTENT(INOUT) :: inFile
  !CLN  CHARACTER(LEN=*),  INTENT(IN)    :: inName
  !CLN  INTEGER :: idummy
  !CLN
  !CLN  READ(inFile(nn:nn+3),'(i4)') idummy
  !CLN  IF (idummy < 1983 .OR. idummy > 1995) THEN
  !CLN    PRINT *, 'Check position of the year number in input gswp file', inFile
  !CLN    STOP
  !CLN  ELSE
  !CLN    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
  !CLN    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
  !CLN  END IF
  !CLN
  !CLNEND SUBROUTINE renameFiles


! ============== PRIVATE SUBROUTINES USED ONLY BY THE MPI MASTER ===============


! MPI: calculates and sends grid decomposition info to the workers
SUBROUTINE master_decomp(comm, mland)

  use mpi

  USE cable_IO_vars_module, ONLY : landpt

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: comm  ! MPI communicator to talk to the workers
  INTEGER, INTENT(IN) :: mland ! total number of landpoints in the global grid
  ! INTEGER, INTENT(IN) :: mp    ! total number of land patches in the global grid

  INTEGER :: lpw  ! average number of landpoints per worker
  INTEGER :: rank, rest, nxt, pcnt, ierr, i, tmp
  INTEGER :: patchcnt  ! sum of patches for a range landpoints

  ! how many workers do we have?
  CALL MPI_Comm_size(comm, wnp, ierr)
  wnp = wnp - 1

  ALLOCATE(wland(wnp), STAT=ierr)
  IF (ierr /= 0) THEN
     ! TODO: print an error message
     write(*,*) 'master-decomp MPI_ABORT'
     CALL MPI_Abort(comm, 16, ierr)
  END IF

  ! MPI: calculate landpoint distribution among the workers
  ! this version distributes landpoints rather than active patches,
  ! but perhaps this will be easy to change in the future?
  lpw = mland / wnp
  IF (lpw == 0) THEN
     write(*,*) 'number of workers > than mland',wnp,mland
     CALL MPI_Abort(comm, 17, ierr)
  END IF
  rest = MOD(mland, wnp)
  nxt = 1

  DO rank = 1, wnp
     wland(rank)%landp0 = nxt
     pcnt = lpw
     ! MPI: poor man's load balancing:
     ! - difference between number of landpoints assigned to
     ! different workers is 1 at most
     ! - inherent load balance because calculations are done on "patches"
     ! whose number differ between landpoints
     ! TODO: - the above to be addressed in the next version
     IF (rest > 0) THEN
        pcnt = pcnt + 1
        rest = rest - 1
     END IF
     wland(rank)%nland = pcnt

     ! MPI: let each worker know their assignement
     ! in this version of cable workers care only about the number of points
     ! CALL MPI_Send(nxt, 1, MPI_INTEGER, rank, 0, comm, ierr)
     CALL MPI_Send(pcnt, 1, MPI_INTEGER, rank, 0, comm, ierr)

     ! MPI: should be the same as landpt(nxt)%cstart
     wland(rank)%patch0 = landpt(nxt)%cstart
     ! MPI: calculate no of patches for pcnt landpoint starting from nxt
     ! MPI: TODO: workers only need to know the number of their patches
     ! or maybe not (if they need patch displacement in the input)

     ! MPI: find number of patches in all landpoints assigned to this
     ! worker (difference between last patch of last landpoint and first of
     ! first)
     patchcnt = landpt(nxt+pcnt-1)%cend - landpt(nxt)%cstart + 1
     wland(rank)%npatch = patchcnt

     ! MPI: development check
     tmp = 0
     DO i = 1, pcnt
        tmp = tmp + landpt(nxt+i-1)%nap
     END DO
     IF (patchcnt /= tmp) THEN
        WRITE(*,*) 'invalid patch number for worker ', patchcnt, tmp, rank
        CALL MPI_Abort(comm, 18, ierr)
     END IF

     ! MPI: at this point workers can't determine patches on their own
     ! so we need to send it explicitely
     ! in this version of cable workers care only about the number of patches
     ! CALL MPI_Send (wland(rank)%patch0, 1, MPI_INTEGER, rank, 0, comm, ierr)
     CALL MPI_Send(wland(rank)%npatch, 1, MPI_INTEGER, rank, 0, comm, ierr)

     nxt = nxt + pcnt
  END DO

  RETURN

END SUBROUTINE master_decomp


! MPI: creates param_t type for the master to scatter the default parameters
! to the workers
! then sends the parameters
! and finally frees the MPI type
SUBROUTINE master_cable_params(comm, met, air, ssnow, veg, bgc, soil, canopy, rough, rad, sum_flux, bal)

  use mpi

  USE cable_def_types_mod
  USE cable_IO_vars_module
  use cable_mpicommon, only: nparam, add_address_hvector

  IMPLICIT NONE

  ! subroutine arguments

  INTEGER,                   INTENT(IN)    :: comm ! MPI communicator
  TYPE(met_type),            INTENT(INOUT) :: met
  TYPE(air_type),            INTENT(INOUT) :: air
  TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow
  TYPE(veg_parameter_type),  INTENT(INOUT) :: veg
  TYPE(bgc_pool_type),       INTENT(INOUT) :: bgc
  TYPE(soil_parameter_type), INTENT(INOUT) :: soil
  TYPE(canopy_type),         INTENT(INOUT) :: canopy
  TYPE(roughness_type),      INTENT(INOUT) :: rough
  TYPE(radiation_type),      INTENT(INOUT) :: rad
  TYPE(sum_flux_type),       INTENT(INOUT) :: sum_flux
  TYPE(balances_type),       INTENT(INOUT) :: bal

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: ierr
  INTEGER, ALLOCATABLE, DIMENSION(:) :: param_t

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride
  INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks
  INTEGER :: rank

  INTEGER :: landpt_t, patch_t

  INTEGER :: nxt, pcnt, off, cnt


  ! create MPI types for exchanging slices of landpt and patch arrays
  CALL decomp_types(landpt_t, patch_t)

  ! MPI: TODO: replace sends with isends
  DO rank = 1, wnp

     ! MPI: send worker's landpts
     nxt = wland(rank)%landp0
     pcnt = wland(rank)%nland
     CALL MPI_Send(landpt(nxt), pcnt, landpt_t, rank, 0, comm, ierr)

     ! MPI: send worker's patch
     nxt = wland(rank)%patch0
     pcnt = wland(rank)%npatch
     CALL MPI_Send(patch(nxt), pcnt, patch_t, rank, 0, comm, ierr)

  END DO

  ! MPI: TODO: free landp_t and patch_t types?

  ntyp = nparam

  ALLOCATE(param_t(wnp))

  ALLOCATE(blen(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mp * extr1
  r2stride = mp * extr2
  istride  = mp * extid

  ! default type is byte, to be overriden for multi-D types
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  DO rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     I1LEN  = cnt * extid
     llen  = cnt * extl

     bidx = 0

     ! the order of variables follows argument list
     ! the order of fields within follows alloc_*_type subroutines

     ! ----------- met --------------

     call add_address_hvector(met%year, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%moy, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%ca, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%doy, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%hod, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%ofsd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%fld, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%precip, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%precip_sn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%tk, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%tvair, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%tvrad, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%pmb, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%ua, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%qv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%qvair, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%da, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%dva, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%coszen, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%Ndep, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%Pdep, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%u10, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%rhum, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(met%fsd, off, cnt, mp, displs, blen, types, bidx)

     ! ----------- air --------------

     call add_address_hvector(air%rho, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%volm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%rlam, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%qsat, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%epsi, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%visc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%psyc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%dsatdk, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(air%cmolar, off, cnt, mp, displs, blen, types, bidx)

     ! ----------- ssnow --------------

     call add_address_hvector(ssnow%isflag, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%iantrct, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%pudsto, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%pudsmx, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%cls, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%dfn_dtg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%dfh_dtg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%dfe_ddq, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%ddq_dtg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%evapsn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%fwtop, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%fwtop1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%fwtop2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%fwtop3, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%osnowd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%potev, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%runoff, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%rnof1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%rnof2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%rtsoil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wbtot1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wbtot2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wb_lake, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%sinfil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%qstss, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wetfac, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%owetfac, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%t_snwlr, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tggav, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%otgg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%otss, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%otss_0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tprecip, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tevap, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%trnoff, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%totenbal, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%totenbal2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%fland, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%ifland, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%qasrf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%qfsrf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%qssrf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%snage, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%snowd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%smelt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%ssdnn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tss, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tss_p, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%deltss, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%owb1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%sconds, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%sdepth, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%smass, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%ssdn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tgg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tggsn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%dtmlt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%albsoilsn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%evapfbl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%tilefrac, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wbtot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%gammzz, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wb, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wbice, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wblf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wbfice, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%S, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Tsoil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%SL, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%TL, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%h0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%rex, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%wflux, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%delwcol, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%zdelta, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%kth, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Tsurface, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%lE, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%evap, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%ciso, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%cisoL, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%rlitt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%thetai, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%snowliq, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%nsteps, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%TsurfaceFR, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Ta_daily, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%nsnow, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Qadv_daily, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%G0_daily, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Qevap_daily, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Qprec_daily, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Qprec_snow_daily, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%E_fusion_sn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%E_sublimation_sn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%latent_heat_sn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%evap_liq_sn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%surface_melt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(ssnow%Qadv_rain_sn, off, cnt, mp, displs, blen, types, bidx)

     ! ----------- veg --------------

     call add_address_hvector(veg%iveg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ivegp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%iLU, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%canst1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%dleaf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ejmax, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ejmax_shade, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ejmax_sun, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%meth, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%frac4, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%hc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vlai, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%xalbnir, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%rp20, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%rpcoef, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%rs20, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%shelrb, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vegcf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%tminvj, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%toptvj, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%tmaxvj, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vbeta, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vcmax, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vcmax_shade, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vcmax_sun, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%xfang, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%extkn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vlaimax, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%wai, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%a1gs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%d0gs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%alpha, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%convex, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%cfrd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%gswmin, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%conkc0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%conko0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ekc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%eko, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%g0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%g1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%vcmaxcc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ejmaxcc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%gmmax, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%gm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%c4kci, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%c4kcc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%bjv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%deciduous, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%refl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%taul, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%froot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%rootbeta, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%gamma, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%ZR, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%F10, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%clitt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%disturbance_interval, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(veg%disturbance_intensity, off, cnt, mp, displs, blen, types, bidx)

     ! ----------- bgc --------------

     call add_address_hvector(bgc%cplant, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bgc%csoil, off, cnt, mp, displs, blen, types, bidx)
     ! ncp, each worker gets the same copy of whole array
     call add_address_hvector(bgc%ratecp, 1, ncp, 1, displs, blen, types, bidx)
     ! ncs, each worker gets the same copy of whole array
     call add_address_hvector(bgc%ratecs, 1, ncs, 1, displs, blen, types, bidx)

     ! ----------- soil --------------

     call add_address_hvector(soil%isoilm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%bch, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%c3, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%clay, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%css, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%hsbh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%hyds, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%i2bp3, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%ibp2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%rhosoil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%sand, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%sfc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%silt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%ssat, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%sucs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%swilt, off, cnt, mp, displs, blen, types, bidx)
     ! ms, each worker gets the same copy of whole array
     call add_address_hvector(soil%zse, 1, ms, 1, displs, blen, types, bidx)
     ! ms+1, each worker gets the same copy of whole array
     call add_address_hvector(soil%zshh, 1, ms+1, 1, displs, blen, types, bidx)
     call add_address_hvector(soil%soilcol, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%albsoilf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%cnsd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%pwb_min, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%albsoil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%nhorizons, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%ishorizon, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%clitt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%zeta, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%fsatmax, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%swilt_vec, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%ssat_vec, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(soil%sfc_vec, off, cnt, mp, displs, blen, types, bidx)

     ! ----------- canopy --------------

     call add_address_hvector(canopy%cansto, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%cduv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%delwc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%dewmm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fe, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fpn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%frp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%frpw, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%frpr, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%frs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fnee, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%frday, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fnv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fev, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%epot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fnpp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fevw_pot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gswx_T, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%cdtq, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%wetfac_cs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fevw, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fhvw, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%oldcansto, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fhv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fns, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fhs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fhs_cor, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ga, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ghflux, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%precis, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%qscrn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%rnet, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%rniso, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%segg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%sghflux, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%through, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%spill, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%tscrn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%wcint, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%tv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%us, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%uscrn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%vlaiw, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%rghlai, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fwet, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%evapfbl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gswx, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%zetar, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%zetash, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fess, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fesp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%dgdtg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fes, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fes_cor, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fevc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ofes, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%A_sl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%A_sh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%A_slC, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%A_shC, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%A_slJ, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%A_shJ, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%GPP_sl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%GPP_sh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fevc_sl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fevc_sh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_A_cs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%dAdcs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_GPP_cs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_A_cs_sl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_A_cs_sh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_fevc_cs_sl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_fevc_cs_sh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%eta_fevc_cs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%cs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%cs_sl, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%cs_sh, off, cnt, mp, displs, blen, types, bidx)
     ! call add_address_hvector(canopy%ci_sl, off, cnt, mp, displs, blen, types, bidx)
     ! call add_address_hvector(canopy%ci_sh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%tlf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%dlf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gw, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ancj, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%tlfy, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ecy, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ecx, off, cnt, mp, displs, blen, types, bidx)
     ! call add_address_hvector(canopy%ci, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%fwsoil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%kthLitt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%DvLitt, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%An, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%Rd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%isc3, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%vcmax, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gammastar, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gsc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gbc, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%gac, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(canopy%ci, off, cnt, mp, displs, blen, types, bidx)

     ! ------- rough -------

     call add_address_hvector(rough%disp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%hruff, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%hruff_grmx, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%rt0us, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%rt1usa, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%rt1usb, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%rt1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%za_uv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%za_tq, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%z0m, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%zref_uv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%zref_tq, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%zruffs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%z0soilsn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%z0soil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%coexp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%usuh, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%term2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%term3, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%term5, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%term6, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rough%term6a, off, cnt, mp, displs, blen, types, bidx)

     ! --------rad --------

     call add_address_hvector(rad%transb, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%albedo_T, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%longitude, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%workp1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%workp2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%workp3, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%extkb, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%extkd2, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%extkd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%flws, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%latitude, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%lwabv, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%qssabs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%transd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%trad, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%fvlai, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%rhocdf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%rniso, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%scalex, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%albedo, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%reffdf, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%reffbm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%extkbm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%extkdm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%fbeam, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%cexpkbm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%cexpkdm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%rhocbm, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%gradis, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(rad%qcan, off, cnt, mp, displs, blen, types, bidx)

     ! ------- sum_flux -----

     call add_address_hvector(sum_flux%sumpn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumrp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumrpw, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumrpr, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumrs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumrd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%dsumpn, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%dsumrp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%dsumrs, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%dsumrd, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumxrp, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(sum_flux%sumxrs, off, cnt, mp, displs, blen, types, bidx)

     ! ------- bal ----

     call add_address_hvector(bal%drybal, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%ebal, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%ebal_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%ebal_cncheck, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%ebal_tot_cncheck, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%ebaltr, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%ebal_tottr, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%evap_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%osnowd0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%precip_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%rnoff_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%wbal, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%wbal_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%wbtot0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%wetbal, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%cansto0, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%owbtot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%evapc_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%evaps_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%rnof1_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%rnof2_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%snowdc_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%wbal_tot1, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%delwc_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%qasrf_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%qfsrf_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%qssrf_tot, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%Radbal, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%EbalSoil, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%Ebalveg, off, cnt, mp, displs, blen, types, bidx)
     call add_address_hvector(bal%Radbalsum, off, cnt, mp, displs, blen, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of param_t fields ', bidx, ', fix it (01)!'
        call MPI_Abort(comm, 19, ierr)
     end if

     call MPI_Type_create_struct(bidx, blen, displs, types, param_t(rank), ierr)
     call MPI_Type_commit(param_t(rank), ierr)

     call MPI_Type_size(param_t(rank), tsize, ierr)
     call MPI_Type_get_extent(param_t(rank), tmplb, text, ierr)

     WRITE(*,*) 'master to rank param_t blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  END DO ! rank

  WRITE (*,*) 'total cable params size sent to all workers: ', localtotal

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blen)

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  CALL MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total cable params size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
     WRITE(*,*) 'error: total length of cable params sent and received differ'
     CALL MPI_Abort (comm, 20, ierr)
  END IF

  ! so, now send all the parameters
  CALL master_send_input(comm, param_t, 0)
  !  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank=1, wnp
     CALL MPI_Type_Free(param_t(rank), ierr)
  END DO

  DEALLOCATE(param_t)

  ! all CABLE parameters have been transferred to the workers by now
  RETURN

END SUBROUTINE master_cable_params


! MPI: creates casa_ts types to broadcast/scatter the default casa parameters
! to all the workers
! then sends them
! and finally frees the MPI type
SUBROUTINE master_casa_params(comm, casabiome, casapool, casaflux, casamet, casabal, phen)

  use mpi

  USE cable_def_types_mod
  use casadimension, only: mlitter, mplant, msoil, mphase, mdyear, mso
  USE casavariable
  USE phenvariable
  use cable_mpicommon, only: ncasaparam, add_address_1block, add_address_hvector

  IMPLICIT NONE

  ! sub arguments
  INTEGER, INTENT(IN) :: comm  ! MPI communicator
  TYPE(casa_biome),    INTENT(INOUT) :: casabiome
  TYPE(casa_pool),     INTENT(INOUT) :: casapool
  TYPE(casa_flux),     INTENT(INOUT) :: casaflux
  TYPE(casa_met),      INTENT(INOUT) :: casamet
  TYPE(casa_balance),  INTENT(INOUT) :: casabal
  TYPE(phen_variable), INTENT(INOUT) :: phen

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: ierr
  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_t

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride
  INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks

  INTEGER :: rank, off, cnt

  !  moved to calling before this subroutine (BP May 2013)
  !  CALL MPI_Bcast (mvtype, 1, MPI_INTEGER, 0, comm, ierr)
  !  CALL MPI_Bcast (mstype, 1, MPI_INTEGER, 0, comm, ierr)

  ntyp = ncasaparam

  ALLOCATE(casa_t(wnp))

  ALLOCATE(blocks(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mp * extr1
  r2stride = mp * extr2
  istride = mp * extid

  ! default type is byte, to be overriden for multi-D types
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  DO rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     I1LEN = cnt * extid
     llen = cnt * extl

     bidx = 0

     ! ------- casabiome -----

     ! all workes get the same casabiome
     call add_address_1block(casabiome%ivt2, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkleafcoldmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkleafcoldexp, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkleafdrymax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkleafdryexp, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%glaimax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%glaimin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%sla, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratiofrootleaf, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%kroot, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%krootlen, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%rootdepth, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%kuptake, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%kminN, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%KuplabP, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%kclabrate, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xnpmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%q10soil, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkoptlitter, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkoptsoil, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%maxfinelitter, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%maxcwd, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%prodptase, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%costnpup, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkplab, 1, mso, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkpsorb, 1, mso, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkpocc, 1, mso, displs, blocks, types, bidx)
     call add_address_1block(casabiome%nintercept, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%nslope, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%plantrate, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%rmplant, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fracnpptoP, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fraclignin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fraclabile, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNCplantmin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNCplantmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioPCplantmin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioPCplantmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fracLigninplant, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ftransNPtoL, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ftransPPtoL, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%litterrate, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%soilrate, 1, mvtype, displs, blocks, types, bidx)
     ! added by ln
     call add_address_1block(casabiome%ratioNPplantmin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNPplantmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%la_to_sa, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%vcmax_scalar, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%disturbance_interval, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_EnzPool, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_KMO2, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_KMcp, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_Ea, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_alpha, 1, mvtype, displs, blocks, types, bidx)

     ! ------ casapool ----

     call add_address_hvector(casapool%Clabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dClabiledt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Cplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Pplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dCplantdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNplantdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPplantdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioPCplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nsoilmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoillab, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoilsorb, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoilocc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNsoilmindt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoillabdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoilsorbdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoiloccdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Clitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nlitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Plitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dClitterdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNlitterdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPlitterdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNClitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioPClitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Csoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dCsoildt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNsoildt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoildt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioPCsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoilnew, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoilmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoilmax, off, cnt, mp, displs, blocks, types, bidx)
     ! added by LN
     call add_address_hvector(casapool%ratioNPplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNPlitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNPsoil, off, cnt, mp, displs, blocks, types, bidx)

     ! ------- casaflux ----

     call add_address_hvector(casaflux%Cgpp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cnpp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crgplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminfix, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminuptake, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Plabuptake, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Clabloss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracClabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracCalloc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracNalloc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracPalloc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crmplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplant_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplant_tot, off, cnt, mp, displs, blocks, types, bidx)
     ! gol124: temp
     ! 3D
     call add_address_hvector(casaflux%fromPtoL, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromPtoL_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cnep, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nmindep, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminloss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nupland, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nlittermin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nsmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nsimm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nsnet, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fNminloss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fNminleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pdep, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pwea, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Ploss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pupland, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Plittermin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psimm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psnet, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fPleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplab, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kpsorb, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kpocc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kmlabP, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psorbmax, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%frac_sapwood, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%sapwood_area, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fHarvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fCrop, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%klitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%klitter_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%klitter_tot, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%ksoil, off, cnt, mp, displs, blocks, types, bidx)
     ! gol124: temp only
     call add_address_hvector(casaflux%fromLtoS, off, cnt, mp, displs, blocks, types, bidx)
     ! gol124: temp only
     call add_address_hvector(casaflux%fromStoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromLtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromStoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxCtolitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxNtolitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxPtolitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxCtosoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxNtosoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxPtosoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxCtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     ! 13C
     call add_address_hvector(casaflux%FluxFromPtoL, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromLtoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromStoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromPtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromLtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromStoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromPtoHarvest, off, cnt, mp, displs, blocks, types, bidx)

     ! ------- casamet ----

     call add_address_hvector(casamet%glai, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tairk, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%precip, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%tsoilavg, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistavg, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%btran, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%lnonwood, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moist, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%iveg2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%ijgcm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%isorder, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%lat, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%lon, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%areacell, off, cnt, mp, displs, blocks, types, bidx)

     ! ------- casabal ----

     call add_address_hvector(casabal%FCgppyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCnppyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrmleafyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrmwoodyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrmrootyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrgrowyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrpyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrsyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCneeyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNdepyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNfixyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNsnetyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNupyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNleachyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNlossyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPweayear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPdustyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPsnetyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPupyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPleachyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPlossyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%glaimon, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%glaimonx, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%cplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%pplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%clitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nlitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%plitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%csoillast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nsoillast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoillast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nsoilminlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoillablast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoilsorblast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoilocclast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%cbalance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nbalance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%pbalance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%sumcbal, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%sumnbal, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%sumpbal, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%clabilelast, off, cnt, mp, displs, blocks, types, bidx)

     ! ------- phen -------

     call add_address_hvector(phen%phase, off, cnt, mp, displs, blocks, types, bidx)
     ! all workers get the same TKshed
     call add_address_hvector(phen%TKshed, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphase, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%phen, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%aphen, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%phasespin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_1, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_3, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_4, off, cnt, mp, displs, blocks, types, bidx)

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid number of casa_t param fields ',bidx,', fix it (02)!'
        CALL MPI_Abort(comm, 21, ierr)
     END IF

     CALL MPI_Type_create_struct(bidx, blocks, displs, types, casa_t(rank), ierr)
     CALL MPI_Type_commit(casa_t(rank), ierr)

     CALL MPI_Type_size(casa_t(rank), tsize, ierr)
     CALL MPI_Type_get_extent(casa_t(rank), tmplb, text, ierr)

     WRITE(*,*) 'master to rank casa_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  END DO ! rank

  WRITE(*,*) 'total casa params size sent to all workers: ', localtotal
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  remotetotal = 0
  CALL MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total casa params size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
     WRITE(*,*) 'error: total length of casa params sent and received differ'
     CALL MPI_Abort (comm, 22, ierr)
  END IF

  CALL MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  CALL master_send_input(comm, casa_t, 0)
  !  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank=1, wnp
     CALL MPI_Type_Free(casa_t(rank), ierr)
  END DO

  DEALLOCATE(casa_t)

  ! all casa parameters have been sent to the workers by now

  RETURN

END SUBROUTINE master_casa_params

! MPI: creates inp_t types to send input data to the workers
! input data means arrays read by get_met_data; each worker receives
! only its own slice of the arrays
SUBROUTINE master_intypes(comm, met, veg)

  use mpi

  USE cable_def_types_mod
  use cable_mpicommon, only: ninput, add_address_hvector

  IMPLICIT NONE

  ! Arguments
  integer,                  intent(in) :: comm
  type(met_type),           intent(in) :: met  ! meteorological data
  type(veg_parameter_type), intent(in) :: veg  ! LAI retrieved from file

  ! Local variables

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
  integer :: tsize, localtotal, remotetotal

  integer(kind=MPI_ADDRESS_KIND) :: r1stride
  integer :: r1len, r2len, i1len, llen ! block lengths
  integer :: bidx ! block index
  integer :: ntyp ! total number of blocks
  integer :: rank ! worker rank
  integer :: off  ! first patch index for a worker
  integer :: cnt  ! mp for a worker
  integer :: ierr

  ALLOCATe(inp_ts(wnp))

  ! max total number of fields to receive: met + veg fields
  !ntyp = 10 + 1
  ntyp = ninput

  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! chunks of all 1D vectors are contiguous blocks of memory so just send them
  ! as blocks of bytes
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  do rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     I1LEN = cnt * extid
     llen = cnt * extl

     r1stride = mp * extr1

     bidx = 0

     ! met fields

     call add_address_hvector(met%fsd, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%tk, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%pmb, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%qv, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%ua, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%precip, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%precip_sn, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%fld, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%ca, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%coszen, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%Ndep, off, cnt, mp, displs, blocks, types, bidx)

     ! veg fields

     call add_address_hvector(veg%vlai, off, cnt, mp, displs, blocks, types, bidx)

     ! additional field missing from previous versions;
     ! added when trying to fix a bug in the new mpi code
     ! the order of these new fields follows the order of their
     ! declaration in cable_define_types.F90

     call add_address_hvector(met%year, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%moy, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%doy, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%hod, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%u10, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(met%rhum, off, cnt, mp, displs, blocks, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid intype nmat, nvec or n3d constant, fix it (03)!'
        call MPI_Abort(comm, 23, ierr)
     end if

     ! marshall all fields into a single MPI derived datatype for worker rank
     CALL MPI_Type_create_struct(bidx, blocks, displs, types, inp_ts(rank), ierr)
     CALL MPI_Type_commit(inp_ts(rank), ierr)
     CALL MPI_Type_size(inp_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent(inp_ts(rank), tmplb, text, ierr)

     write(*,*) 'master to ', rank, ': intype (01) struct blocks, size, extent and lb: ', &
                 bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  end do ! rank

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  write(*,*) 'total input data size sent to all workers: ', localtotal

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  call MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total input data size received by all workers: ', remotetotal

  if (localtotal /= remotetotal) then
     write(*,*) 'error: total length of input data sent and received differ (01)'
     call MPI_Abort(comm, 24, ierr)
  end if

  return

END SUBROUTINE master_intypes

! MPI: creates out_t types to receive the results from the workers
SUBROUTINE master_outtypes(comm,met,canopy,ssnow,rad,bal,air,soil,veg)

  use mpi

  USE cable_def_types_mod
  use cable_mpicommon, only: n3d, nmat, nvec, add_address_hvector

  IMPLICIT NONE

  INTEGER :: comm ! MPI communicator to talk to the workers

  TYPE(met_type),            INTENT(IN)    :: met
  TYPE(canopy_type),         INTENT(IN)    :: canopy
  TYPE(soil_snow_type),      INTENT(IN)    :: ssnow
  TYPE(radiation_type),      INTENT(IN)    :: rad
  TYPE(balances_type),       INTENT(INOUT) :: bal
  TYPE(air_type),            INTENT(IN)    :: air
  TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
  TYPE(veg_parameter_type),  INTENT(IN)    :: veg  ! vegetation parameters

  ! MPI: displacement (address) vector for vectors (1D arrays)
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: vaddr
  ! MPI: displacement (address) vector for matrices (2D arrays)
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: maddr
  ! MPI: displacement (address) vector for 3D arrays
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: m3daddr

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: midx, vidx, ierr

  ! base index to make types indexing easier
  INTEGER :: istart

  INTEGER :: i

  ! MPI: calculate the sizes/extents of Fortran types used by
  ! CABLE
  ! gol124: commented out because called earlier in mpidrv_master
  ! before it calls master_intypes
  ! CALL find_extents

  ! MPI: allocate matrices to hold datatype handles
  ALLOCATE(m3d_t(n3d, wnp))
  ALLOCATE(mat_t(nmat, wnp))
  ALLOCATE(vec_t(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = n3d + nmat + nvec
  ALLOCATE(blocks(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  ! MPI: allocate vector to hold handles for the combined type
  ALLOCATE(recv_ts(wnp))

  ! MPI: type_struct that combines all data for a single worker

  ! MPI: allocate address vectors for 3D matrices
  ! currently only 1 is 3D but to make future expansion easier
  ! let's handle it the same way as 3D and 2D
  ALLOCATE(m3daddr(n3d))

  ! MPI: allocate address vectors for matrices
  ALLOCATE(maddr(nmat))
  ! MPI: allocate address vectors for vectors
  ALLOCATE(vaddr(nvec))
  ! MPI: allocate vector block lengths
  ALLOCATE(blen(nvec))

  ! MPI: TODO: global var that holds the total number of patches/landunits?
  r1stride = mp * extr1
  r2stride = mp * extr2

  totalrecv = 0
  ! create a separate datatype for each worker
  ! TODO: replace all hvectors with mpi_type_create_subarray!
  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2

     ! ------------- 3D arrays -------------

     ! rad 3D
     ! midx = 0
     ! call add_address_hvector(rad%qcan, off, cnt, displs, blocks, types, midx)
     ! m3d_t(midx, rank) = types(midx)
     ! m3daddr(midx) = displs(midx)
     ! call MPI_Type_commit(m3d_t(midx, rank), ierr)
     CALL MPI_Get_address (rad%qcan(off,1,1), m3daddr(1), ierr)
     CALL MPI_Type_create_hvector (mf*nrb, r1len, r1stride, MPI_BYTE, &
 &                        m3d_t(1, rank), ierr)
     CALL MPI_Type_commit (m3d_t(1, rank), ierr)

     ! ------------- 2D arrays -------------

     ! MPI: an hvector type for each vector, maddr contains displacements
     ! for bundling these hvectors into the struct later
     ! block length is number of patches/worker * type extent
     ! stride is global number of patches * type extent
     ! repeat/no of blocks is the 2nd rank

     ! met 2D
     midx = 0

     midx = midx + 1
     CALL MPI_Get_address (met%fsd(off,1), maddr(midx), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     ! canopy 2D
     !     ! REAL(r_1)
     !     CALL MPI_Get_address (canopy%rwater(off,1), maddr(midx), ierr) ! 1
     !     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
     ! &                        mat_t(midx, rank), ierr)
     !     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     !     midx = midx + 1
     ! REAL(r_2)
     ! MPI: gol124: backport to r1134 changes r_2 to r_1
     ! MPI: gol124: in newest CABLE-cnp it's r_2 again
     midx = midx + 1
     CALL MPI_Get_address (canopy%evapfbl(off,1), maddr(midx), ierr) ! 2
     ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address (canopy%gswx(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address (canopy%zetar(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector (niter, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     ! 13C
     ! done as canopy%gswx
     midx = midx + 1
     CALL MPI_Get_address(canopy%An(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%Rd(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%vcmax(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%gammastar(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%gsc(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%gbc(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%gac(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(canopy%ci(off,1), maddr(midx), ierr) ! 2
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     ! ssnow 2D
     midx = midx + 1
     CALL MPI_Get_address(ssnow%dtmlt(off,1), maddr(midx), ierr)
     CALL MPI_Type_create_hvector(3, r1len, r1stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%albsoilsn(off,1), maddr(midx), ierr) ! 3
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (ssnow%gammzz(off,1), maddr(midx), ierr) ! 4
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%sconds(off,1), maddr(midx), ierr) ! 5
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%sdepth(off,1), maddr(midx), ierr) ! 6
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%smass(off,1), maddr(midx), ierr) ! 7
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     ! MPI: r1134 does not know about this field, comment out
     !midx = midx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (ssnow%dtmlt(off,1), maddr(midx), ierr) ! 8
     !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
     ! &                        mat_t(midx, rank), ierr)
     !CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%ssdn(off,1), maddr(midx), ierr) ! 9
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%tgg(off,1), maddr(midx), ierr) ! 10
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%tggsn(off,1), maddr(midx), ierr) ! 11
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (ssnow%wb(off,1), maddr(midx), ierr) ! 12
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%evapfbl(off,1), maddr(midx), ierr) ! 12
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%wbfice(off,1), maddr(midx), ierr) ! 13
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (ssnow%wbice(off,1), maddr(midx), ierr) ! 14
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (ssnow%wblf(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     ! additional for sli
     midx = midx + 1
     CALL MPI_Get_address(ssnow%S(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector(ms, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(ssnow%Tsoil(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector(ms, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(ssnow%thetai(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector(ms, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(ssnow%snowliq(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector(msn, r2len, r2stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address(ssnow%sconds(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector(msn, r1len, r1stride, MPI_BYTE, &
          mat_t(midx, rank), ierr)
     CALL MPI_Type_commit(mat_t(midx, rank), ierr)
     ! end additional for sli

     ! rad 2D
     midx = midx + 1
     CALL MPI_Get_address (rad%fbeam(off,1), maddr(midx), ierr) ! 102
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%albedo(off,1), maddr(midx), ierr) ! 16
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%fvlai(off,1), maddr(midx), ierr) ! 17
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (rad%gradis(off,1), maddr(midx), ierr) ! 18
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%rhocdf(off,1), maddr(midx), ierr) ! 19
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%rniso(off,1), maddr(midx), ierr) ! 20
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%scalex(off,1), maddr(midx), ierr) ! 21
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%reffdf(off,1), maddr(midx), ierr) ! 22
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%reffbm(off,1), maddr(midx), ierr) ! 23
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%extkbm(off,1), maddr(midx), ierr) ! 24
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%extkdm(off,1), maddr(midx), ierr) ! 25
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%cexpkbm(off,1), maddr(midx), ierr) ! 26
     ! Maciej: cexpkbm is mp*swb
     !     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
     ! &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%cexpkdm(off,1), maddr(midx), ierr) ! 27
     ! Maciej: cexpkdm is mp*swb
     !     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
     ! &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%rhocbm(off,1), maddr(midx), ierr) ! 27
     ! Maciej: rhocbm is mp*nrb
     !     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
     ! &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     ! air 2D - all fields 1D - skipped

     ! soil 2D
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%albsoil(off,1), maddr(midx), ierr) ! 28
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     ! veg 2D
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%refl(off,1), maddr(midx), ierr) ! 29
     CALL MPI_Type_create_hvector (2, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%taul(off,1), maddr(midx), ierr) ! 29
     CALL MPI_Type_create_hvector (2, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%froot(off,1), maddr(midx), ierr) ! 29
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
 &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     ! MPI: sanity check
     IF (midx /= nmat) THEN
        WRITE(*,*) 'master: outtype invalid nmat ',midx,' constant, fix it (04)!'
        CALL MPI_Abort (comm, 25, ierr)
     END IF

     ! ------------- 1D arrays -------------

     vidx = 0

     ! met

     call add_address_hvector(met%ca, off, cnt, mp, vaddr, blen, types, vidx)
     ! MPI: gol124: changed to 2D and moved up when Bernard
     ! ported to CABLE_r491
     call add_address_hvector(met%fld, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%precip, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%precip_sn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%tk, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%tvair, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%tvrad, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%pmb, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%ua, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%qv, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%qvair, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%da, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%dva, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(met%coszen, off, cnt, mp, vaddr, blen, types, vidx)

     ! canopy

     call add_address_hvector(canopy%fess, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fesp, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%cansto, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%cduv, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%delwc, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%dewmm, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%dgdtg, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fe, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fh, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fpn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%A_sh, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%A_sl, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%A_slC, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%A_shC, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%A_slJ, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%A_shJ, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%GPP_sh, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%GPP_sl, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%eta_A_cs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%dAdcs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%cs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%eta_GPP_cs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%eta_fevc_cs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%dlf, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%vcmax_shade, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%vcmax_sun, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%ejmax_shade, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%ejmax_sun, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%frp, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%frpw, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%frpr, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%frs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fnee, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%frday, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fnv, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fev, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fevc, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fevw, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fhv, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fhvw, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fns, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fes, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fes_cor, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fhs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fhs_cor, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fwet, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%epot, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fnpp, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fevw_pot, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%gswx_T, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%cdtq, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%wetfac_cs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%ga, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%ghflux, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%precis, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%qscrn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%rnet, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%segg, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%sghflux, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%spill, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%through, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%tscrn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%tv, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%us, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%uscrn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%vlaiw, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%rghlai, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%wcint, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(canopy%fwsoil, off, cnt, mp, vaddr, blen, types, vidx)
     ! 13C
     call add_address_hvector(canopy%isc3, off, cnt, mp, vaddr, blen, types, vidx)

     ! ssnow

     call add_address_hvector(ssnow%pudsto, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%pudsmx, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%cls, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%dfn_dtg, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%dfh_dtg, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%dfe_ddq, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%ddq_dtg, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%evapsn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%fwtop, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%fwtop1, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%fwtop2, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%fwtop3, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%isflag, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%osnowd, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%potev, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%pwb_min, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%runoff, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%rnof1, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%rnof2, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%rtsoil, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%snage, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%snowd, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%smelt, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%ssdnn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%tss, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%wbtot, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%wb_lake, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%sinfil, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%qstss, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%wetfac, off, cnt, mp, vaddr, blen, types, vidx)
     ! MPI: TODO: maybe not needed for transfer to master?
     call add_address_hvector(ssnow%owetfac, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%t_snwlr, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%tggav, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%otss, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%otss_0, off, cnt, mp, vaddr, blen, types, vidx)

     ! rad

     call add_address_hvector(rad%extkb, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%extkd2, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%extkd, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%flws, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%latitude, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%lwabv, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%qssabs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%transd, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%trad, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(rad%transb, off, cnt, mp, vaddr, blen, types, vidx)

     ! bal

     ! MPI: gol124: changed to 2D and moved up when Bernard
     call add_address_hvector(bal%drybal, off, cnt, mp, vaddr, blen, types, vidx)
     ! MPI: remove ebal from exchanged data, calculate temp val on the master
     call add_address_hvector(bal%osnowd0, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(bal%wbtot0, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(bal%wetbal, off, cnt, mp, vaddr, blen, types, vidx)

     ! air

     call add_address_hvector(air%rho, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%volm, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%rlam, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%qsat, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%epsi, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%visc, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%psyc, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%dsatdk, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(air%cmolar, off, cnt, mp, vaddr, blen, types, vidx)

     ! soil

     call add_address_hvector(soil%bch, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%c3, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%clay, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%cnsd, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%css, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%hsbh, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%hyds, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%i2bp3, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%ibp2, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%isoilm, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%rhosoil, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%rs20, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%sand, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%sfc, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%silt, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%ssat, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%sucs, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(soil%swilt, off, cnt, mp, vaddr, blen, types, vidx)

     ! veg

     call add_address_hvector(veg%iveg, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%meth, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%vlai, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%canst1, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%ejmax, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%frac4, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%wai, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%vegcf, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%tminvj, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%tmaxvj, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%vbeta, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%xalbnir, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%hc, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%shelrb, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%vcmax, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%xfang, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%dleaf, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%rp20, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%rpcoef, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%extkn, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(veg%deciduous, off, cnt, mp, vaddr, blen, types, vidx)

     ! additional for SLI
     call add_address_hvector(ssnow%Tsurface, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%h0, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%delwcol, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%evap, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%nsnow, off, cnt, mp, vaddr, blen, types, vidx)
     call add_address_hvector(ssnow%nsteps, off, cnt, mp, vaddr, blen, types, vidx)
     ! end additional for SLI

     ! MPI: sanity check
     if (vidx /= nvec) then
        write(*,*) 'master: outtype invalid nvec ', vidx, ' constant, fix it (05)!'
        call MPI_Abort(comm, 26, ierr)
     end if

     ! MPI: all vectors into a single hindexed type for each worker
     call MPI_Type_create_hindexed(vidx, blen, vaddr, MPI_BYTE, vec_t(rank), ierr)
     call MPI_Type_commit(vec_t(rank), ierr)

     ! MPI: TODO: each worker's all types into a single struct
     ! 1 slice of 3D array (hvector), displacement offset
     ! TODO: displ(1) not set!
     istart = 0
     blocks(istart + 1) = 1 ! 1 copy of hvector
     displs(istart + 1) = m3daddr(istart + 1)
     types(istart + 1) = m3d_t(istart + 1,rank)

     ! 29 slices of matrices (hvector), displacements maddr array
     istart = istart + n3d
     DO i = 1, midx
        !blocks(istart + i) = extent or size of mat_t?
        blocks(istart + i) = 1 ! single copy of mat_t
        displs(istart + i) = maddr(i)
        types(istart + i) = mat_t(i, rank)
     END DO
     ! 163 hindexed, displacements vaddr
     istart = istart + midx
     DO i = 1, vidx
        blocks(istart + i) = blen(i)
        displs(istart + i) = vaddr(i)
        types(istart + i) = MPI_BYTE
     END DO
     ! total: 184

     CALL MPI_Type_create_struct (vidx+midx+1, blocks, displs, types, recv_ts(rank), ierr)
     CALL MPI_Type_commit (recv_ts(rank), ierr)

     CALL MPI_Type_size (recv_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (recv_ts(rank), tmplb, text, ierr)

     WRITE(*,*) 'master: data recv from ', rank, ': size, extent, lb: ', tsize, text, tmplb

     totalrecv = totalrecv + tsize

  END DO

  WRITE(*,*) 'total data size received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total data size sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE(*,*) 'error master: totalsend and totalrecv differ',totalsend,totalrecv
          CALL MPI_Abort (comm, 27, ierr)
  END IF

  DEALLOCATE(blen)
  DEALLOCATE(vaddr)
  DEALLOCATE(maddr)
  DEALLOCATE(m3daddr)
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN

END SUBROUTINE master_outtypes


! MPI: creates handles for receiving casa final results from the workers
subroutine master_casa_types(comm, casabiome, casapool, casaflux, casamet, &
     casabal, phen)

#ifndef __INTEL__
  use mpi,            only: MPI_Address_kind, MPI_Abort, &
       MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
       MPI_Type_get_extent, MPI_Type_free, MPI_Reduce, &
       MPI_In_place, MPI_Integer, MPI_Sum, MPI_Byte
#else
  use mpi
#endif
  use cable_def_types_mod, only: mp, mvtype
  use casadimension,   only: mso
  use casavariable,    only: casa_biome, casa_pool, casa_flux, casa_met, &
       casa_balance, ncasa_biome, ncasa_pool, ncasa_flux, ncasa_met, &
       ncasa_bal
  use phenvariable,    only: phen_variable, ncasa_phen
  use cable_mpicommon, only: add_address_hvector

  implicit none

  integer :: comm ! MPI communicator to talk to the workers

  type(casa_biome),    intent(INOUT) :: casabiome
  type(casa_pool),     intent(INOUT) :: casapool
  type(casa_flux),     intent(INOUT) :: casaflux
  type(casa_met),      intent(INOUT) :: casamet
  type(casa_balance),  intent(INOUT) :: casabal
  type(phen_variable), intent(INOUT) :: phen

  ! MPI: temp arrays for marshalling all types into a struct
  integer, allocatable, dimension(:) :: blocks
  integer(MPI_Address_kind), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types
  integer :: ntyp ! number of worker's types

  integer :: i

  integer :: tsize, totalrecv, totalsend
  integer(MPI_Address_kind) :: text, tmplb

  integer :: rank, off, cnt
  integer :: bidx, ierr

  allocate(casa_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = ncasa_biome + ncasa_pool + ncasa_flux + ncasa_met + &
       ncasa_bal + ncasa_phen
  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  do rank=1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     bidx = 0

     ! casabiome

     call add_address_hvector(casabiome%ivt2, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkleafcoldmax, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkleafcoldexp, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkleafdrymax, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkleafdryexp, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%glaimax, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%glaimin, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%sla, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratiofrootleaf, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%kroot, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%krootlen, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%rootdepth, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%kuptake, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%kminN, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%KuplabP, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%kclabrate, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xnpmax, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%q10soil, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkoptlitter, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkoptsoil, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkplab, 1, mso, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkpsorb, 1, mso, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%xkpocc, 1, mso, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%prodptase, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%costnpup, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%maxfinelitter, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%maxcwd, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%nintercept, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%nslope, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%la_to_sa, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%vcmax_scalar, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%disturbance_interval, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%DAMM_EnzPool, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%DAMM_KMO2, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%DAMM_KMcp, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%DAMM_Ea, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%DAMM_alpha, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%plantrate, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%rmplant, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%fracnpptoP, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%fraclignin, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%fraclabile, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratioNCplantmin, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratioNCplantmax, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratioNPplantmin, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratioNPplantmax, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%fracLigninplant, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ftransNPtoL, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ftransPPtoL, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%litterrate, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratioPcplantmax, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%ratioPcplantmin, 1, mvtype, mvtype, displs, blocks, types, bidx)
     call add_address_hvector(casabiome%soilrate, 1, mvtype, mvtype, displs, blocks, types, bidx)

     ! casapool

     call add_address_hvector(casapool%Clabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dClabiledt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Ctot, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Ctot_0, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Cplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Pplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dCplantdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNplantdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPplantdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNPplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nsoilmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoillab, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoilsorb, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoilocc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNsoilmindt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoillabdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoilsorbdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoiloccdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Clitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nlitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Plitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dClitterdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNlitterdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPlitterdt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNClitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNPlitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Csoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Nsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%Psoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dCsoildt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dNsoildt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%dPsoildt, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNPsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoilnew, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoilmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioNCsoilmax, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioPCsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioPCplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ratioPClitter, off, cnt, mp, displs, blocks, types, bidx)

     ! casaflux

     call add_address_hvector(casaflux%Cgpp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cnpp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crgplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminfix, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminuptake, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Plabuptake, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Clabloss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracClabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%stemnpp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%frac_sapwood, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%sapwood_area, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Charvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nharvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pharvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fharvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fcrop, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracCalloc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracNalloc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fracPalloc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crmplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cplant_turnover, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromPtoL, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cnep, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Crsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nmindep, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminloss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nminleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nupland, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nlittermin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nsmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nsimm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Nsnet, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fNminloss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fNminleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pdep, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pwea, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Ploss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Pupland, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Plittermin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psimm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psnet, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fPleach, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplab, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kpsorb, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kpocc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kmlabP, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Psorbmax, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cplant_turnover_disturbance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cplant_turnover_crowding, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%Cplant_turnover_resource_limitation, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%klitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%ksoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromLtoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromStoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromLtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromStoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtolitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxNtolitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxPtolitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtosoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxNtosoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxPtosoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtohwp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxNtohwp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxPtohwp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtoclear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxNtoclear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxPtoclear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%CtransferLUC, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fromPtoL_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%klitter_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%klitter_tot, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplant_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%kplant_tot, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtoCO2_plant_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtoCO2_litter_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxfromPtoCO2_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fluxfromLtoCO2_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxNtoAtm_fire, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromPtoL, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromLtoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromStoS, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromPtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromLtoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromStoCO2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxFromPtoHarvest, off, cnt, mp, displs, blocks, types, bidx)

     ! casamet

     call add_address_hvector(casamet%glai, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tairk, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%precip, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%tsoilavg, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistavg, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%btran, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%lnonwood, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moist, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%iveg2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%ijgcm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%isorder, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%lat, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%lon, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%areacell, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tairkspin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%cgppspin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%crmplantspin_1, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%crmplantspin_2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%crmplantspin_3, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoilspin_1, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoilspin_2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoilspin_3, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoilspin_4, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoilspin_5, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%Tsoilspin_6, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistspin_1, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistspin_2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistspin_3, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistspin_4, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistspin_5, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moistspin_6, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%mtempspin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%frecspin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%cAn12spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%cAn13spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%dprecip_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%aprecip_av20_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%du10_max_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%drhum_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%dtemp_max_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%dtemp_min_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%KBDI_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%D_MacArthur_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%FFDI_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%last_precip_spin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%DSLR_spin, off, cnt, mp, displs, blocks, types, bidx)

     ! casabal

     call add_address_hvector(casabal%FCgppyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCnppyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrmleafyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrmwoodyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrmrootyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrgrowyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrpyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCrsyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FCneeyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%dCdtyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%LAImax, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%Cleafmean, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%Crootmean, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNdepyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNfixyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNsnetyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNupyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNleachyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FNlossyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPweayear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPdustyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPsnetyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPupyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPleachyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%FPlossyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%glaimon, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%glaimonx, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%cplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%pplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%clitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nlitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%plitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%csoillast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nsoillast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoillast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nsoilminlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoillablast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoilsorblast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%psoilocclast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%cbalance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%nbalance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%pbalance, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%sumcbal, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%sumnbal, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%sumpbal, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%clabilelast, off, cnt, mp, displs, blocks, types, bidx)

     ! phen

     call add_address_hvector(phen%Tkshed, 1, mvtype, 1, displs, blocks, types, bidx)
     call add_address_hvector(phen%phase, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphase, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%phen, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%aphen, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%phasespin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_1, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_2, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_3, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphasespin_4, off, cnt, mp, displs, blocks, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of casa fields, fix it (06)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 28, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, &
          casa_ts(rank), ierr)
     call MPI_Type_commit(casa_ts(rank), ierr)

     call MPI_Type_size(casa_ts(rank), tsize, ierr)
     call MPI_Type_get_extent(casa_ts(rank), tmplb, text, ierr)

     write(*,*) 'casa results recv from worker, size, extent, lb: ', &
          rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     ! TODO: also free partial types for intypes, outtypes etc.
     do i=1, ntyp
        if (types(i) /= MPI_Byte) then
           call MPI_Type_free(types(i), ierr)
        end if
     end do

  end do

  write(*,*) 'total size of casa results received from all workers: ', &
       totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  call MPI_Reduce(MPI_In_place, totalsend, 1, MPI_Integer, MPI_Sum, &
       0, comm, ierr)

  write(*,*) 'total size of casa results sent by all workers: ', totalsend

  if (totalrecv /= totalsend) then
     write(*,*) 'error: casa results totalsend and totalrecv differ'
     call MPI_Abort(comm, 29, ierr)
  end if

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  return

end subroutine master_casa_types


subroutine master_climate_types(comm, climate)

#ifndef __INTEL__
  use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Abort, &
       MPI_Byte, &
       MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
       MPI_Type_get_extent, MPI_Reduce, MPI_In_place, MPI_Integer, &
       MPI_Sum, MPI_Abort, MPI_Isend, MPI_Bottom, MPI_Waitall
#else
  use mpi
#endif
  use cable_def_types_mod, only: climate_type, mp

  integer,            intent(in)    :: comm
  type(climate_type), intent(inout) :: climate

  ! MPI: temp arrays for marshalling all types into a struct
  integer, allocatable, dimension(:) :: blocks
  integer(MPI_Address_kind), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types
  integer :: ntyp ! number of worker's types

  integer :: tsize, totalrecv, totalsend
  integer(MPI_Address_kind) :: text, tmplb

  integer :: rank, off, cnt
  integer :: bidx, ierr

  allocate(climate_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nclimate

  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

 do rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     bidx = 0

     ! ------------- scalars  -------------

     bidx = bidx + 1
     call MPI_Get_address(climate%nyear_average, displs(bidx), ierr)
     blocks(bidx) = extid
     types(bidx)  = MPI_Byte

     bidx = bidx + 1
     call MPI_Get_address(climate%nday_average, displs(bidx), ierr)
     blocks(bidx) = extid
     types(bidx)  = MPI_Byte

     bidx = bidx + 1
     call MPI_Get_address(climate%nyears, displs(bidx), ierr)
     blocks(bidx) = extid
     types(bidx)  = MPI_Byte

     bidx = bidx + 1
     call MPI_Get_address(climate%doy, displs(bidx), ierr)
     blocks(bidx) = extid
     types(bidx)  = MPI_Byte

     ! ------------- arrays -------------

     call add_address_hvector(climate%chilldays, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%iveg, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%biome, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%GMD, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%modis_igbp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%DSLR, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%NDAY_Nesterov, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dtemp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_min, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_min20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_max20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%qtemp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mmoist, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp_min, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%qtemp_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%qtemp_max_last_year, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp_min20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp_max20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%atemp_mean, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%AGDD5, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%GDD5, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%AGDD0, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%GDD0, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%alpha_PT, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%evap_PT, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%aevap, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%alpha_PT20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%GDD0_rec, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%frec, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dtemp_min, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%fdorm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%fapar_ann_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%fapar_ann_max_last_year, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%AvgAnnMaxFAPAR, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dtemp_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%drhum, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%du10_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dprecip, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%aprecip, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%aprecip_av20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%last_precip, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%KBDI, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%FFDI, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%D_MacArthur, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Nesterov_Current, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Nesterov_ann_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Nesterov_ann_max_last_year, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Nesterov_ann_running_max, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp_min_20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%mtemp_max_20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_min_20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_max_20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dtemp_31, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dmoist_31, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%alpha_PT_20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%dtemp_91, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%APAR_leaf_sun, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%APAR_leaf_shade, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Dleaf_sun, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Dleaf_shade, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Tleaf_sun, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Tleaf_shade, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%cs_sun, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%cs_shade, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%scalex_sun, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%scalex_shade, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%fwsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%aprecip_20, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Rd_sun, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%Rd_shade, off, cnt, mp, displs, blocks, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of climate fields, fix it (07)!'
        call MPI_Abort(comm, 30, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, &
          climate_ts(rank), ierr)
     call MPI_Type_commit(climate_ts(rank), ierr)

     call MPI_Type_size(climate_ts(rank), tsize, ierr)
     call MPI_Type_get_extent(climate_ts(rank), tmplb, text, ierr)

     write(*,*) 'climate results recv from worker, size, extent, lb: ', &
          rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize
  end do

  write(*,*) 'total size of climate results received from all workers: ', &
       totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  call MPI_Reduce(MPI_In_place, totalsend, 1, MPI_Integer, MPI_Sum, &
       0, comm, ierr)

  write(*,*) 'total size of climate results sent by all workers: ', &
       totalsend

  if (totalrecv /= totalsend) then
     write(*,*) 'error: climate results totalsend and totalrecv differ'
     call MPI_Abort(comm, 31, ierr)
  end if

  do rank=1, wnp
     call MPI_Isend(MPI_Bottom, 1, climate_ts(rank), rank, 0, comm, &
          inp_req(rank), ierr)
  end do

  call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  return

end subroutine master_climate_types


! MPI: creates datatype handles to receive restart data from workers
subroutine master_restart_types(comm, canopy, air, veg, ssnow)

#ifndef __INTEL__
  use mpi,                 only: MPI_Address_kind, MPI_Abort, &
       MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
       MPI_Type_get_extent, MPI_Type_free, MPI_Reduce, &
       MPI_In_place, MPI_Integer, MPI_Sum
#else
  use mpi
#endif
  use cable_def_types_mod, only: mp, ms, canopy_type, air_type, &
       veg_parameter_type, soil_snow_type
  use cable_mpicommon,     only: nrestart, add_address_hvector

  implicit none

  integer :: comm ! MPI communicator to talk to the workers

  type(canopy_type),        intent(in) :: canopy
  type(air_type),           intent(in) :: air
  type(veg_parameter_type), intent(in) :: veg
  type(soil_snow_type),     intent(in) :: ssnow

  ! MPI: temp arrays for marshalling all types into a struct
  integer, allocatable, dimension(:) :: blocks
  integer(MPI_Address_kind), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types
  integer :: ntyp ! number of worker's types

  integer :: last2d, i

  integer :: tsize, totalrecv, totalsend
  integer(MPI_Address_kind) :: text, tmplb

  integer :: rank, off, cnt
  integer :: bidx, ierr

  allocate(restart_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nrestart
  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  do rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     bidx = 0

     ! ------------- 2D arrays -------------

     call add_address_hvector(canopy%evapfbl, off, cnt, mp, displs, blocks, types, bidx)

     last2d = bidx

     ! ------------- 1D vectors -------------

     call add_address_hvector(canopy%cduv, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%dewmm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%dgdtg, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%frpw, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%frpr, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%fnv, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%rho, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%volm, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%qsat, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%epsi, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%visc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%psyc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%dsatdk, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(air%cmolar, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(ssnow%otss, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(ssnow%wetfac, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%fwsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(canopy%us, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(veg%cfrd, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(veg%vlai, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(veg%hc, off, cnt, mp, displs, blocks, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of restart fields, fix it (09)!'
        write(*,*) 'bidx: ', bidx
        write(*,*) 'ntyp: ', ntyp
        call MPI_Abort(comm, 32, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, &
          restart_ts(rank), ierr)
     call MPI_Type_commit(restart_ts(rank), ierr)

     call MPI_Type_size(restart_ts(rank), tsize, ierr)
     call MPI_Type_get_extent(restart_ts(rank), tmplb, text, ierr)

     write(*,*) 'restart results recv from worker, size, extent, lb: ', &
          rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     do i=1, last2d
        call MPI_Type_free(types(i), ierr)
     end do
  end do

  write(*,*) 'total size of restart fields received from all workers: ', &
       totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  call MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
       0, comm, ierr)

  write(*,*) 'total size of restart fields sent by all workers: ', &
       totalsend

  if (totalrecv /= totalsend) then
     write(*,*) 'error: restart fields totalsend and totalrecv differ'
     call MPI_Abort(comm, 33, ierr)
  end if

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  return

end subroutine master_restart_types


! MPI: Casa - dump read and write
SUBROUTINE master_casa_dump_types(comm, casamet, casaflux, phen, climate, c13o2flux )

  use mpi

  use casadimension,       only: icycle, mphase, mplant
  use casavariable,        only: casa_met, casa_flux
  use cable_def_types_mod, only: mp, ms, climate_type
  use phenvariable
  use cable_mpicommon,     only: ncdumprw, add_address_hvector
  ! 13C
  use cable_common_module, only: cable_user
  use cable_c13o2_def,     only: c13o2_flux

  IMPLICIT NONE

  INTEGER,             INTENT(IN)    :: comm
  TYPE(casa_met),      INTENT(INOUT) :: casamet
  TYPE(casa_flux),     INTENT(INOUT) :: casaflux
  TYPE(phen_variable), INTENT(INOUT) :: phen
  TYPE(climate_type),  INTENT(INOUT) :: climate
  ! 13C
  type(c13o2_flux),    intent(inout) :: c13o2flux

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, Istride
  INTEGER :: r1len, r2len, I1LEN ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks
  INTEGER :: rank ! worker rank
  INTEGER :: off  ! first patch index for a worker
  INTEGER :: cnt  ! mp for a worker
  INTEGER :: ierr

  ALLOCATE(casa_dump_ts(wnp))

  ntyp = ncdumprw + icycle - 1
  ! 13C
  if (cable_user%c13o2) ntyp = ntyp + 2
  if (cable_user%call_blaze) ntyp = ntyp + 11

  ALLOCATE(blocks(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  ! chunks of all 1D vectors are contiguous blocks of memory so just send them
  ! as blocks of bytes
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  DO rank=1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1

     r1stride = mp * extr1

     I1LEN = cnt * extid
     istride = mp * extid

     r2len = cnt * extr2
     r2stride = mp * extr2

     ! casamet fields

     bidx = 0

     call add_address_hvector(casamet%tairk, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%tsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casamet%moist, off, cnt, mp, displs, blocks, types, bidx)

     ! casaflux fields

     call add_address_hvector(casaflux%cgpp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%crmplant, off, cnt, mp, displs, blocks, types, bidx)

     ! phen fields

     call add_address_hvector(phen%phase, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(phen%doyphase, off, cnt, mp, displs, blocks, types, bidx)

     ! climate fields (for acclimation of autotrophic resp)

     call add_address_hvector(climate%qtemp_max_last_year, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(climate%frec, off, cnt, mp, displs, blocks, types, bidx)
     ! climate fields (for BLAZE)

     if (cable_user%call_blaze) then
        call add_address_hvector(climate%dprecip, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%aprecip_av20, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%du10_max, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%drhum, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%dtemp_max, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%dtemp_min, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%KBDI, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%D_MacArthur, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%FFDI, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%DSLR, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(climate%last_precip, off, cnt, mp, displs, blocks, types, bidx)
     end if

     ! N and P deposition

     if (icycle > 1) then
        call add_address_hvector(casaflux%Nmindep, off, cnt, mp, displs, blocks, types, bidx)
     end if

     if (icycle > 2) then
        call add_address_hvector(casaflux%Pdep, off, cnt, mp, displs, blocks, types, bidx)
     end if

     ! 13C

     if (cable_user%c13o2) then
        call add_address_hvector(c13o2flux%cAn12, off, cnt, mp, displs, blocks, types, bidx)
        call add_address_hvector(c13o2flux%cAn, off, cnt, mp, displs, blocks, types, bidx)
     end if

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid intype in master_casa_dump, fix it (10)!'
        CALL MPI_Abort(comm, 34, ierr)
     END IF

     ! marshall all fields into a single MPI derived datatype for worker rank
     CALL MPI_Type_create_struct(bidx, blocks, displs, types, casa_dump_ts(rank), ierr)
     CALL MPI_Type_commit(casa_dump_ts(rank), ierr)

     CALL MPI_Type_size(casa_dump_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent(casa_dump_ts(rank), tmplb, text, ierr)

     WRITE(*,*) 'master to ', rank, ': intype (03) struct blocks, size, extent and lb: ', &
          bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  END DO ! ranks

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  WRITE(*,*) 'total casa_dump size sent to all workers: ', localtotal

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  ! write(*,*) 'b3.1 reduce wk', MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr
  ! call flush(6)
  CALL MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
  ! write(*,*) 'b3.2 reduce wk', MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr
  ! call flush(6)

  WRITE(*,*) 'total input data size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
     WRITE(*,*) 'error: total length of input data sent and received differ (03)'
     CALL MPI_Abort(comm, 35, ierr)
  END IF

  !  DO rank = 1, wnp
  !
  !     CALL MPI_ISend (MPI_BOTTOM, 1, casa_dump_ts(rank), rank, 0, comm, &
  ! &               inp_req(rank), ierr)
  !
  !  END DO
  !
  !  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

END SUBROUTINE master_casa_dump_types

! #############################################################################################################

! MPI: Casa-LUC: exchanging casapools between master and worker, as required for LUC updates
SUBROUTINE master_casa_LUC_types(comm, casapool, casabal, casaflux)

  use mpi

  use cable_def_types_mod, only: mp
  USE casadimension,       ONLY: mplant, mlitter, msoil
  USE casavariable,        ONLY: casa_pool, casa_balance, casa_flux
  use cable_mpicommon,     only: nLUCrw, add_address_hvector

  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: comm
  TYPE(casa_pool),    INTENT(IN) :: casapool
  TYPE(casa_balance), INTENT(IN) :: casabal
  TYPE(casa_flux),    INTENT(IN) :: casaflux

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, Istride
  INTEGER :: r1len, r2len, i1len ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks
  INTEGER :: rank ! worker rank
  INTEGER :: off  ! first patch index for a worker
  INTEGER :: cnt  ! mp for a worker
  INTEGER :: ierr
  integer :: last2d, i

  ALLOCATE(casa_LUC_ts(wnp))

  ntyp = nLUCrw

  ALLOCATE(blocks(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  ! chunks of all 1D vectors are contiguous blocks of memory so just send them
  ! as blocks of bytes
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  DO rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len    = cnt * extr1
     r1stride = mp  * extr1

     i1len    = cnt * extid
     istride  = mp  * extid

     r2len    = cnt * extr2
     r2stride = mp  * extr2

     bidx = 0

     ! casapool fields 2D
     call add_address_hvector(casapool%cplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%clitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%csoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%nplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%nlitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%nsoil, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%pplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%plitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%psoil, off, cnt, mp, displs, blocks, types, bidx)
     ! casabal fields 2D
     call add_address_hvector(casabal%cplantlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%clitterlast, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%csoillast, off, cnt, mp, displs, blocks, types, bidx)

     last2d = bidx

     ! casapool fields 1D
     call add_address_hvector(casapool%Nsoilmin, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%clabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casapool%ctot, off, cnt, mp, displs, blocks, types, bidx)
     ! casabal fields 1D
     call add_address_hvector(casabal%FCneeyear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casabal%clabilelast, off, cnt, mp, displs, blocks, types, bidx)
     ! casaflux fields 1D
     call add_address_hvector(casaflux%fHarvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%cHarvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%nHarvest, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%fcrop, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtohwp, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%FluxCtoclear, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(casaflux%CtransferLUC, off, cnt, mp, displs, blocks, types, bidx)

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid intype in master_casa_LUC, fix it (11)!'
        CALL MPI_Abort(comm, 36, ierr)
     END IF

     ! marshall all fields into a single MPI derived datatype for worker rank

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, casa_LUC_ts(rank), ierr)
     CALL MPI_Type_commit (casa_LUC_ts(rank), ierr)

     CALL MPI_Type_size (casa_LUC_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (casa_LUC_ts(rank), tmplb, text, ierr)

     WRITE(*,*) 'master to ',rank,': intype (04) struct blocks, size, extent and lb: ', &
          bidx,tsize,text,tmplb

     localtotal = localtotal + tsize

     ! free the partial types used for matrices
     do i=1, last2d
        call MPI_Type_free(types(i), ierr)
     end do

  END DO ! ranks

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  WRITE(*,*) 'total casa_LUC size sent to all workers: ', localtotal

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  CALL MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total input data size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
     WRITE(*,*) 'error: total length of input data sent and received differ (04)'
     CALL MPI_Abort(comm, 37, ierr)
  END IF

END SUBROUTINE master_casa_LUC_types

! *******************************************************************************************

! MPI
! Creates pop_ts types to broadcast/cscatter the default POP parameters
! to all workers

SUBROUTINE master_pop_types(comm, pop)

  use mpi
  USE POP_mpi
  USE POP_types,          ONLY: pop_type
  USE cable_common_module,ONLY: cable_user
  USE casavariable,       ONLY: casa_met

  IMPLICIT NONE

  INTEGER,        INTENT(IN)    :: comm
  TYPE(pop_type), INTENT(INOUT) :: pop

  INTEGER :: rank ! worker rank
  INTEGER :: off  ! first patch index for a worker
  INTEGER :: cnt  ! mp for a worker
  INTEGER :: ierr
  INTEGER :: prev_mp
  INTEGER, ALLOCATABLE :: w_iwood(:)

  ! Also send Pop relevant info to workers.

  prev_mp = 0
  DO rank = 1, wnp

     off = wland(rank)%patch0

     ! number of forest or shrub patches per worker
     wland(rank)%npop_iwood = &
          COUNT(POP%Iwood .GE.off .AND. POP%Iwood .LT. (off+wland(rank)%npatch))

     CALL MPI_Send (wland(rank)%npop_iwood, 1, MPI_INTEGER, rank, 0, comm, ierr)

     ! Index mapping for pop-pixels for each worker
     ALLOCATE( wland(rank)%iwood(wland(rank)%npop_iwood) )
     ALLOCATE( w_iwood(wland(rank)%npop_iwood) )

     ! Array of indeces in master domain
     wland(rank)%iwood(:) = PACK(POP%Iwood, (POP%Iwood.GE.off .AND. POP%Iwood .LT. (off+wland(rank)%npatch)))

     ! Index of patch in worker domain
     w_iwood(:) = wland(rank)%iwood(:) - prev_mp
     CALL MPI_Send (w_iwood, wland(rank)%npop_iwood, MPI_INTEGER, rank, 0, comm, ierr )

     DEALLOCATE(w_iwood)

     ! accounted patches so far
     prev_mp = prev_mp + wland(rank)%npatch
  END DO

  ! once generate a structure fitting the pop derived types
  CALL create_pop_gridcell_type(pop_ts, comm)

  ! Send restart stuff once
  IF ( .NOT. cable_user%POP_fromZero ) THEN
     off = 1
     DO rank = 1, wnp

        IF ( rank .GT. 1 ) off = off + wland(rank-1)%npop_iwood

        cnt = wland(rank)%npop_iwood
        ! Maciej: worker_pop_types always call MPI_Recv, even if cnt is zero
        !        IF ( cnt .EQ. 0 ) CYCLE
        CALL MPI_Send( POP%pop_grid(off), cnt, pop_ts, rank, 0, comm, ierr )
        CALL MPI_Send( POP%it_pop(off), cnt, MPI_INTEGER, rank, 0, comm, ierr )

     END DO

  END IF

END SUBROUTINE master_pop_types

SUBROUTINE master_receive_pop(POP, comm)

  USE MPI
  USE POP_mpi
  USE POP_Types, ONLY: pop_type

  IMPLICIT NONE

  INTEGER,        INTENT(IN)    :: comm
  TYPE(pop_type), INTENT(INOUT) :: pop

  INTEGER :: ierr, rank, off, cnt, stat(MPI_STATUS_SIZE)

  off = 1
  DO rank=1, wnp

     IF (rank .GT. 1) off = off + wland(rank-1)%npop_iwood

     cnt = wland(rank)%npop_iwood
     IF (cnt .EQ. 0) CYCLE

     CALL MPI_Recv(POP%pop_grid(off), cnt, pop_ts, rank, 0, comm, stat, ierr)
     CALL MPI_Recv(POP%it_pop(off), cnt, MPI_INTEGER, rank, 0, comm, stat, ierr)

  END DO

END SUBROUTINE master_receive_pop

! 13C
! MPI: creates c13o2_flux_ts types to broadcast/scatter the default c13o2_flux parameters to all the workers
! then sends them
! and finally frees the MPI type
subroutine master_c13o2_flux_params(comm, c13o2flux)

  use mpi
  ! use mpi,                 only: &
  !      MPI_ADDRESS_KIND, MPI_STATUS_SIZE, MPI_BYTE, &
  !      MPI_Get_address, MPI_Type_create_hvector, &
  !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
  !      MPI_Type_get_extent, MPI_Reduce, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, &
  !      MPI_Barrier
  use cable_def_types_mod, only: mp, mf
  use cable_c13o2_def,     only: c13o2_flux
  use cable_mpicommon,     only: nc13o2_flux, add_address_hvector

  implicit none

  integer,          intent(in)    :: comm  ! mpi communicator
  type(c13o2_flux), intent(inout) :: c13o2flux

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types

  ! temp vars for verifying block number and tot1al length of inp_t
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
  integer :: tsize, localtotal, remotetotal

  integer :: ierr
  integer, allocatable, dimension(:) :: c13o2_flux_t

  integer(kind=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride
  integer :: r1len, r2len, i1len, llen ! block lengths
  integer :: bidx ! block index
  integer :: ntyp ! total number of blocks

  integer :: rank, off, cnt

  ntyp = nc13o2_flux

  allocate(c13o2_flux_t(wnp))

  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mp * extr1
  r2stride = mp * extr2
  istride  = mp * extid

  ! default type is byte, to be overriden for multi-D types
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  do rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0 ! global module variable
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     i1len = cnt * extid
     llen  = cnt * extl

     bidx = 0

     ! 1D
     call add_address_hvector(c13o2flux%ca, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%cAn12, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%cAn, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%RAn, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Vstarch, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Rstarch, off, cnt, mp, displs, blocks, types, bidx)
     ! 2D
     call add_address_hvector(c13o2flux%An, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Disc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Rsucrose, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Rphoto, off, cnt, mp, displs, blocks, types, bidx)
     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_flux_t param fields ', bidx, ', fix it (13)!'
        call MPI_Abort(comm, 38, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_flux_t(rank), ierr)
     call MPI_Type_commit(c13o2_flux_t(rank), ierr)

     call MPI_Type_size(c13o2_flux_t(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_flux_t(rank), tmplb, text, ierr)

     write(*,*) 'master to rank c13o2_flux_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  end do ! rank

  write(*,*) 'total c13o2_flux params size sent to all workers: ', localtotal
  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  ! MPI: check whether total size of received data equals total data sent by all the workers
  remotetotal = 0
  call MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total c13o2_flux params size received by all workers: ', remotetotal

  if (localtotal /= remotetotal) then
     write(*,*) 'error: total length of c13o2_flux params sent and received differ'
     call MPI_Abort(comm, 39, ierr)
  end if

  call MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  call master_send_input(comm, c13o2_flux_t, 0)
  ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  do rank=1, wnp
     call MPI_Type_Free(c13o2_flux_t(rank), ierr)
  end do

  deallocate(c13o2_flux_t)

  ! all c13o2_flux parameters have been sent to the workers by now
  return

end subroutine master_c13o2_flux_params


! MPI: creates c13o2_pool_ts types to broadcast/scatter the default c13o2_pool parameters to all the workers
! then sends them
! and finally frees the MPI type
subroutine master_c13o2_pool_params(comm, c13o2pools)

  use mpi
  ! use mpi,                 only: &
  !      MPI_ADDRESS_KIND, MPI_STATUS_SIZE, MPI_BYTE, &
  !      MPI_Get_address, MPI_Type_create_hvector, &
  !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
  !      MPI_Type_get_extent, MPI_Reduce, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, &
  !      MPI_Barrier
  use cable_def_types_mod, only: mp
  use casadimension,       only: mplant, mlitter, msoil
  use cable_c13o2_def,     only: c13o2_pool
  use cable_mpicommon,     only: nc13o2_pool, add_address_hvector

  implicit none

  integer,          intent(in)    :: comm  ! mpi communicator
  type(c13o2_pool), intent(inout) :: c13o2pools

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types

  ! temp vars for verifying block number and tot1al length of inp_t
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
  integer :: tsize, localtotal, remotetotal

  integer :: ierr
  integer, allocatable, dimension(:) :: c13o2_pool_t

  integer(kind=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride
  integer :: r1len, r2len, i1len, llen ! block lengths
  integer :: bidx ! block index
  integer :: ntyp ! total number of blocks

  integer :: rank, off, cnt

  ntyp = nc13o2_pool

  allocate(c13o2_pool_t(wnp))

  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mp * extr1
  r2stride = mp * extr2
  istride  = mp * extid

  ! default type is byte, to be overriden for multi-D types
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  do rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%patch0 ! global module variable
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     i1len = cnt * extid
     llen  = cnt * extl

     bidx = 0

     ! 1D
     call add_address_hvector(c13o2pools%clabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2pools%charvest, off, cnt, mp, displs, blocks, types, bidx)
     ! 2D
     call add_address_hvector(c13o2pools%cplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2pools%clitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2pools%csoil, off, cnt, mp, displs, blocks, types, bidx)
     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_pool_t param fields ', bidx, ', fix it (14)!'
        call MPI_Abort(comm, 40, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_pool_t(rank), ierr)
     call MPI_Type_commit(c13o2_pool_t(rank), ierr)

     call MPI_Type_size(c13o2_pool_t(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_pool_t(rank), tmplb, text, ierr)

     write(*,*) 'master to rank c13o2_pool_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  end do ! rank

  write(*,*) 'total c13o2_pool params size sent to all workers: ', localtotal

  ! MPI: check whether total size of received data equals total data sent by all the workers
  remotetotal = 0
  call MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total c13o2_pool params size received by all workers: ', remotetotal

  if (localtotal /= remotetotal) then
     write(*,*) 'error: total length of c13o2_pool params sent and received differ'
     call MPI_Abort(comm, 41, ierr)
  end if

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  call MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  call master_send_input(comm, c13o2_pool_t, 0)
  ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  do rank=1, wnp
     call MPI_Type_Free(c13o2_pool_t(rank), ierr)
  end do

  deallocate(c13o2_pool_t)

  ! all c13o2_pool parameters have been sent to the workers by now
  return

end subroutine master_c13o2_pool_params


! MPI: creates c13o2_luc_ts types to broadcast/scatter the default c13o2_luc parameters to all the workers
! then sends them
! and finally frees the MPI type
subroutine master_c13o2_luc_params(comm, c13o2luc)

  use mpi
  ! use mpi,                 only: &
  !      MPI_ADDRESS_KIND, MPI_STATUS_SIZE, MPI_BYTE, &
  !      MPI_Get_address, MPI_Type_create_hvector, &
  !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
  !      MPI_Type_get_extent, MPI_Reduce, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, &
  !      MPI_Barrier
  use cable_def_types_mod, only: mland
  ! use cable_def_types_mod, only: mp
  use cable_c13o2_def,     only: c13o2_luc
  use cable_mpicommon,     only: nc13o2_luc, add_address_hvector

  implicit none

  integer,         intent(in)    :: comm  ! mpi communicator
  type(c13o2_luc), intent(inout) :: c13o2luc

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types

  ! temp vars for verifying block number and tot1al length of inp_t
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
  integer :: tsize, localtotal, remotetotal

  integer :: ierr
  integer, allocatable, dimension(:) :: c13o2_luc_t

  integer(kind=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride
  integer :: r1len, r2len, i1len, llen ! block lengths
  integer :: bidx ! block index
  integer :: ntyp ! total number of blocks

  integer :: rank, off, cnt

  allocate(c13o2_luc_t(wnp))

  ntyp = nc13o2_luc
  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mland * extr1
  r2stride = mland * extr2
  istride  = mland * extid
  ! r1stride = mp * extr1
  ! r2stride = mp * extr2
  ! istride  = mp * extid

  ! default type is byte, to be overriden for multi-D types
  types = MPI_BYTE

  ! total size of input data sent to all workers
  localtotal = 0

  ! create a separate MPI derived datatype for each worker
  do rank = 1, wnp

     ! starting patch and number for each worker rank
     off = wland(rank)%landp0 ! global module variable
     cnt = wland(rank)%nland

     r1len = cnt * extr1
     r2len = cnt * extr2
     i1len = cnt * extid
     llen  = cnt * extl

     bidx = 0

     ! 1D
     call add_address_hvector(c13o2luc%cagric, off, cnt, mland, displs, blocks, types, bidx)
     ! 2D
     call add_address_hvector(c13o2luc%charvest, off, cnt, mland, displs, blocks, types, bidx)
     call add_address_hvector(c13o2luc%cclearance, off, cnt, mland, displs, blocks, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_luc_t param fields ', bidx, ', fix it (15)!'
        call MPI_Abort(comm, 42, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_luc_t(rank), ierr)
     call MPI_Type_commit(c13o2_luc_t(rank), ierr)

     call MPI_Type_size(c13o2_luc_t(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_luc_t(rank), tmplb, text, ierr)

     write(*,*) 'master to rank c13o2_luc_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  end do ! rank

  write(*,*) 'total c13o2_luc params size sent to all workers: ', localtotal
  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  ! MPI: check whether total size of received data equals total data sent by all the workers
  remotetotal = 0
  call MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total c13o2_luc params size received by all workers: ', remotetotal

  if (localtotal /= remotetotal) then
     write(*,*) 'error: total length of c13o2_luc params sent and received differ'
     call MPI_Abort(comm, 43, ierr)
  end if

  call MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  call master_send_input(comm, c13o2_luc_t, 0)
  ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  do rank=1, wnp
     call MPI_Type_Free(c13o2_luc_t(rank), ierr)
  end do

  deallocate(c13o2_luc_t)

  ! all c13o2_luc parameters have been sent to the workers by now
  return

end subroutine master_c13o2_luc_params


! MPI: creates handles for receiving c13o2_flux final results from the workers
subroutine master_c13o2_flux_types(comm, c13o2flux)

  use mpi
  ! use mpi,                 only: &
  !      MPI_BYTE, MPI_ADDRESS_KIND, MPI_Get_address, MPI_Type_create_hvector, &
  !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
  !      MPI_Type_get_extent, MPI_Type_free, MPI_Reduce, MPI_IN_PLACE, &
  !      MPI_INTEGER, MPI_SUM
  use cable_def_types_mod, only: mp, mf
  use cable_c13o2_def,     only: c13o2_flux
  use cable_mpicommon,     only: nc13o2_flux, add_address_hvector

  implicit none

  integer,          intent(in)    :: comm ! mpi communicator to talk to the workers
  type(c13o2_flux), intent(inout) :: c13o2flux

  ! MPI: temp arrays for marshalling all types into a struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types
  integer :: ntyp ! number of worker's types

  integer :: last2d, i

  ! mpi: block lengths and strides for hvector representing matrices
  integer :: r1len, r2len, i1len
  integer(kind=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  integer :: tsize, totalrecv, totalsend
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb

  integer :: rank, off, cnt
  integer :: bidx, ierr

  allocate(c13o2_flux_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nc13o2_flux
  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  r1stride = mp * extr1
  r2stride = mp * extr2
  istride  = mp * extid
  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  do rank=1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     i1len = cnt * extid

     bidx = 0

     ! ------------- 2D arrays -------------

     call add_address_hvector(c13o2flux%An, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Disc, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Rsucrose, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Rphoto, off, cnt, mp, displs, blocks, types, bidx)
     ! ------------- 1D vectors -------------

     last2d = bidx

     call add_address_hvector(c13o2flux%ca, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%cAn12, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%cAn, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%RAn, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Vstarch, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2flux%Rstarch, off, cnt, mp, displs, blocks, types, bidx)
     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_flux fields, fix it (16)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 44, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_flux_ts(rank), ierr)
     call MPI_Type_commit(c13o2_flux_ts(rank), ierr)

     call MPI_Type_size(c13o2_flux_ts(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_flux_ts(rank), tmplb, text, ierr)

     write(*,*) 'c13o2_flux results recv from worker, size, extent, lb: ', &
          rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     ! TODO: also free partial types for intypes, outtypes etc.
     do i=1, last2d
        call MPI_Type_free(types(i), ierr)
     end do

  end do ! rank

  write(*,*) 'total size of c13o2_flux results received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  call MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
       0, comm, ierr)

  write(*,*) 'total size of c13o2_flux results sent by all workers: ', totalsend

  if (totalrecv /= totalsend) then
     write(*,*) 'error: c13o2_flux results totalsend and totalrecv differ'
     call MPI_Abort(comm, 45, ierr)
  end if

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  return

end subroutine master_c13o2_flux_types


! MPI: creates handles for receiving c13o2_pool final results from the workers
subroutine master_c13o2_pool_types(comm, c13o2pools)

  use mpi
  ! use mpi,                 only: &
  !      MPI_BYTE, MPI_ADDRESS_KIND, MPI_Get_address, MPI_Type_create_hvector, &
  !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
  !      MPI_Type_get_extent, MPI_Type_free, MPI_Reduce, MPI_IN_PLACE, &
  !      MPI_INTEGER, MPI_SUM
  use cable_def_types_mod, only: mp
  use casadimension,       only: mplant, mlitter, msoil
  use cable_c13o2_def,     only: c13o2_pool
  use cable_mpicommon,     only: nc13o2_pool, add_address_hvector

  implicit none

  integer,          intent(in)    :: comm ! mpi communicator to talk to the workers
  type(c13o2_pool), intent(inout) :: c13o2pools

  ! MPI: temp arrays for marshalling all types into a struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types
  integer :: ntyp ! number of worker's types

  integer :: last2d, i

  ! mpi: block lengths and strides for hvector representing matrices
  integer :: r1len, r2len, i1len
  integer(kind=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  integer :: tsize, totalrecv, totalsend
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb

  integer :: rank, off, cnt
  integer :: bidx, ierr

  allocate(c13o2_pool_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nc13o2_pool
  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  r1stride = mp * extr1
  r2stride = mp * extr2
  istride  = mp * extid
  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  do rank=1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     i1len = cnt * extid

     bidx = 0

     ! ------------- 2D arrays -------------

     call add_address_hvector(c13o2pools%cplant, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2pools%clitter, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2pools%csoil, off, cnt, mp, displs, blocks, types, bidx)
     ! ------------- 1D vectors -------------

     last2d = bidx

     call add_address_hvector(c13o2pools%clabile, off, cnt, mp, displs, blocks, types, bidx)
     call add_address_hvector(c13o2pools%charvest, off, cnt, mp, displs, blocks, types, bidx)
     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_pool fields, fix it (17)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 46, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_pool_ts(rank), ierr)
     call MPI_Type_commit(c13o2_pool_ts(rank), ierr)

     call MPI_Type_size(c13o2_pool_ts(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_pool_ts(rank), tmplb, text, ierr)

     write(*,*) 'c13o2_pool results recv from worker, size, extent, lb: ', &
          rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     ! TODO: also free partial types for intypes, outtypes etc.
     do i=1, last2d
        call MPI_Type_free(types(i), ierr)
     end do

  end do ! rank

  write(*,*) 'total size of c13o2_pool results received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  call MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
       0, comm, ierr)

  write(*,*) 'total size of c13o2_pool results sent by all workers: ', totalsend

  if (totalrecv /= totalsend) then
     write(*,*) 'error: c13o2_pool results totalsend and totalrecv differ'
     call MPI_Abort(comm, 47, ierr)
  end if

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  return

end subroutine master_c13o2_pool_types


! MPI: creates handles for receiving c13o2_luc final results from the workers
subroutine master_c13o2_luc_types(comm, c13o2luc)

  use mpi
  ! use mpi,                 only: &
  !      MPI_BYTE, MPI_ADDRESS_KIND, MPI_Get_address, MPI_Type_create_hvector, &
  !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
  !      MPI_Type_get_extent, MPI_Type_free, MPI_Reduce, MPI_IN_PLACE, &
  !      MPI_INTEGER, MPI_SUM
  use cable_def_types_mod, only: mland
  ! use cable_def_types_mod, only: mp
  use cable_c13o2_def,     only: c13o2_luc
  use cable_mpicommon,     only: nc13o2_luc, add_address_hvector

  implicit none

  integer,         intent(in)    :: comm ! mpi communicator to talk to the workers
  type(c13o2_luc), intent(inout) :: c13o2luc

  ! MPI: temp arrays for marshalling all types into a struct
  integer, allocatable, dimension(:) :: blocks
  integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
  integer, allocatable, dimension(:) :: types
  integer :: ntyp ! number of worker's types

  integer :: last2d, i

  ! mpi: block lengths and strides for hvector representing matrices
  integer :: r1len, r2len, i1len
  integer(kind=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  integer :: tsize, totalrecv, totalsend
  integer(kind=MPI_ADDRESS_KIND) :: text, tmplb

  integer :: rank, off, cnt
  integer :: bidx, ierr

  allocate(c13o2_luc_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nc13o2_luc
  allocate(blocks(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! chunks of all 1D vectors are contiguous blocks of memory so just send them
  ! as blocks of bytes
  types = MPI_BYTE ! initialize

  r1stride = mland * extr1
  r2stride = mland * extr2
  istride  = mland * extid
  ! r1stride = mp * extr1
  ! r2stride = mp * extr2
  ! istride  = mp * extid

  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  do rank=1, wnp
     off = wland(rank)%landp0
     cnt = wland(rank)%nland

     r1len = cnt * extr1
     r2len = cnt * extr2
     i1len = cnt * extid

     bidx = 0

     ! ------------- 2D arrays -------------

     call add_address_hvector(c13o2luc%charvest, off, cnt, mland, displs, blocks, types, bidx)
     call add_address_hvector(c13o2luc%cclearance, off, cnt, mland, displs, blocks, types, bidx)

     ! ------------- 1D vectors -------------

     last2d = bidx

     call add_address_hvector(c13o2luc%cagric, off, cnt, mland, displs, blocks, types, bidx)

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_luc fields, fix it (18)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 48, ierr)
     end if

     call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_luc_ts(rank), ierr)
     call MPI_Type_commit(c13o2_luc_ts(rank), ierr)

     call MPI_Type_size(c13o2_luc_ts(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_luc_ts(rank), tmplb, text, ierr)

     write(*,*) 'c13o2_luc results recv from worker, size, extent, lb: ', &
          rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     do i=1, last2d
        call MPI_Type_free(types(i), ierr)
     end do

  end do ! rank

  deallocate(types)
  deallocate(displs)
  deallocate(blocks)

  write(*,*) 'total size of c13o2_luc results received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  call MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total size of c13o2_luc results sent by all workers: ', totalsend

  if (totalrecv /= totalsend) then
     write(*,*) 'error: c13o2_luc results totalsend and totalrecv differ'
     call MPI_Abort(comm, 49, ierr)
  end if

  return

end subroutine master_c13o2_luc_types
! 13C


! MPI: scatters input data for timestep ktau to all workers
SUBROUTINE master_send_input(comm, dtypes, ktau)

  use mpi

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: comm
  INTEGER, DIMENSION(:), INTENT(IN) :: dtypes
  INTEGER, INTENT(IN) :: ktau    ! timestep

  INTEGER :: rank, ierr

!  IF (.NOT. ALLOCATED(inp_req)) THEN
!     ALLOCATE (inp_req(wnp))
!  END IF

  DO rank=1, wnp
     CALL MPI_ISend(MPI_BOTTOM, 1, dtypes(rank), rank, ktau, comm, inp_req(rank), ierr)
  END DO

  !IF (.NOT. ALLOCATED(inp_stats)) THEN
  !   ALLOCATE (inp_stats(MPI_STATUS_SIZE, wnp))
  !END IF
  CALL MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  RETURN

END SUBROUTINE master_send_input

! receives model output variables from the workers for a single timestep
! uses the timestep value as the message tag
SUBROUTINE master_receive(comm, ktau, types)

  use mpi

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: comm
  INTEGER, INTENT(IN) :: ktau    ! timestep
  ! array with MPI datatype handles for receiving
  INTEGER, DIMENSION(:), INTENT(IN) :: types

  ! local vars
  INTEGER :: rank, ierr

  !IF (.NOT. ALLOCATED(recv_req)) THEN
  !   ALLOCATE (recv_req(wnp))
  !END IF
  !IF (.NOT. ALLOCATED(recv_stats)) THEN
  !   ALLOCATE (recv_stats(MPI_STATUS_SIZE, wnp))
  !END IF

  DO rank=1, wnp
      CALL MPI_Irecv(MPI_BOTTOM, 1, types(rank), rank, ktau, comm, recv_req(rank), ierr)
  END DO

  ! MPI: we need all timesteps data before processing/saving
  !write(*,*) 'waitall 3'
  CALL MPI_Waitall(wnp, recv_req, recv_stats, ierr)

  RETURN

END SUBROUTINE master_receive

! TODO: receives variables that are required by create_restart
! but not write_output
! gol124: how about call master_receive (comm, ktau, restart_ts)
! instead of a separate receive_restart sub?
!SUBROUTINE receive_restart (comm,ktau,dels,soil,veg,ssnow, &
! &              canopy,rough,rad,bgc,bal)
!
!
!  RETURN
!
!END SUBROUTINE receive_restart


! frees memory used for data structures specific to the master
SUBROUTINE master_end(icycle, restart)

  use mpi
  use cable_common_module, only: cable_user

  IMPLICIT NONE

  INTEGER :: icycle ! casa flag
  LOGICAL :: restart ! restart flag

  INTEGER :: rank, i, ierr

  ! MPI: free MPI types
  DO rank=1, wnp

     CALL MPI_Type_free(inp_ts(rank), ierr)

     CALL MPI_Type_free(recv_ts(rank), ierr)

     IF (restart) THEN
        CALL MPI_Type_free(restart_ts(rank), ierr)
     END IF

     IF (icycle>0) THEN
        CALL MPI_Type_free(casa_ts(rank), ierr)
        ! 13C
        if (cable_user%c13o2) then
           call MPI_Type_free(c13o2_flux_ts(rank), ierr)
           call MPI_Type_free(c13o2_pool_ts(rank), ierr)
        end if
     END IF
     !MC - LUC is not freed. Does freeing objects close files as well?

     ! gol124: TODO: m3d_t, mat_t and vec_t can
     ! be freed at the end of master_outtypes
     DO i=1, n3d
        CALL MPI_Type_free(m3d_t(i, rank), ierr)
     END DO

     DO i=1, nmat
        CALL MPI_Type_free(mat_t(i, rank), ierr)
     END DO

     CALL MPI_Type_free(vec_t(rank), ierr)

  END DO

  ! MPI: free arrays for receive requests and statuses
  IF (ALLOCATED(recv_stats)) THEN
     DEALLOCATE(recv_stats)
  END IF
  IF (ALLOCATED(recv_req)) THEN
     DEALLOCATE(recv_req)
  END IF

  ! MPI: free arrays for send requests and statuses
  IF (ALLOCATED(inp_stats)) THEN
     DEALLOCATE(inp_stats)
  END IF
  IF (ALLOCATED(inp_req)) THEN
     DEALLOCATE(inp_req)
  END IF

  ! MPI: free derived datatype handle array
  DEALLOCATE(inp_ts)

  DEALLOCATE(recv_ts)

  IF (restart) THEN
     DEALLOCATE(restart_ts)
  END IF

  IF (icycle>0) THEN
     DEALLOCATE(casa_ts)
     if (cable_user%c13o2) then
        deallocate(c13o2_flux_ts)
        deallocate(c13o2_pool_ts)
     end if
  END IF

  ! MPI: free partial derived datatype handle arrays
  DEALLOCATE(vec_t)
  DEALLOCATE(mat_t)
  DEALLOCATE(m3d_t)

  ! MPI: free landpoint decomposition info
  DEALLOCATE(wland)

  RETURN

END SUBROUTINE master_end

subroutine copy_save_variables(metin, phenin, metout, phenout)
  ! Copy variables from dump files from in to out to save them
  ! because they will come back as zero from the MPI workers.
  ! Use the same function to copy them back after receiving
  ! (the zeros) from the workers.
  !
  ! Used in spincasa and casaonly_luc.

  use casavariable,        only: casa_met
  use phenvariable,        only: phen_variable
  use cable_common_module, only: cable_user

  implicit none

  type(casa_met),      intent(in)    :: metin
  type(phen_variable), intent(in)    :: phenin
  type(casa_met),      intent(inout) :: metout
  type(phen_variable), intent(inout) :: phenout

  metout%Tairkspin      = metin%Tairkspin
  metout%Tsoilspin_1    = metin%Tsoilspin_1
  metout%Tsoilspin_2    = metin%Tsoilspin_2
  metout%Tsoilspin_3    = metin%Tsoilspin_3
  metout%Tsoilspin_4    = metin%Tsoilspin_4
  metout%Tsoilspin_5    = metin%Tsoilspin_5
  metout%Tsoilspin_6    = metin%Tsoilspin_6
  metout%moistspin_1    = metin%moistspin_1
  metout%moistspin_2    = metin%moistspin_2
  metout%moistspin_3    = metin%moistspin_3
  metout%moistspin_4    = metin%moistspin_4
  metout%moistspin_5    = metin%moistspin_5
  metout%moistspin_6    = metin%moistspin_6
  metout%cgppspin       = metin%cgppspin
  metout%crmplantspin_1 = metin%crmplantspin_1
  metout%crmplantspin_2 = metin%crmplantspin_2
  metout%crmplantspin_3 = metin%crmplantspin_3
  phenout%phasespin      = phenin%phasespin
  phenout%doyphasespin_1 = phenin%doyphasespin_1
  phenout%doyphasespin_2 = phenin%doyphasespin_2
  phenout%doyphasespin_3 = phenin%doyphasespin_3
  phenout%doyphasespin_4 = phenin%doyphasespin_4
  metout%mtempspin = metin%mtempspin
  metout%frecspin  = metin%frecspin
  if (cable_user%c13o2) then
     metout%cAn12spin = metin%cAn12spin
     metout%cAn13spin = metin%cAn13spin
  end if
  if (cable_user%call_blaze) then
     metout%dprecip_spin      = metin%dprecip_spin
     metout%aprecip_av20_spin = metin%aprecip_av20_spin
     metout%du10_max_spin     = metin%du10_max_spin
     metout%drhum_spin        = metin%drhum_spin
     metout%dtemp_max_spin    = metin%dtemp_max_spin
     metout%dtemp_max_spin    = metin%dtemp_max_spin
     metout%KBDI_spin         = metin%KBDI_spin
     metout%D_MacArthur_spin  = metin%D_MacArthur_spin
     metout%FFDI_spin         = metin%FFDI_spin
     metout%DSLR_spin         = metin%DSLR_spin
     metout%last_precip_spin  = metin%last_precip_spin
  end if

end subroutine copy_save_variables

! 13C
SUBROUTINE master_spincasacnp(dels, kstart, kend, mloop, veg, soil, casabiome, casapool, &
     casaflux, casamet, casabal, phen, POP, climate, c13o2flux, c13o2pools, icomm, ocomm)

  use cable_def_types_mod
  use cable_carbon_module
  use cable_common_module, only: cable_user
  use casadimension
  use casaparm
  use casavariable
  use phenvariable
  use casa_cable,          only: read_casa_dump
  use casa_inout,          only: casa_fluxout, write_casa_restart_nc ! , casa_poolout
  use POP_types,           only: POP_type
  use cable_pop_io,        only: pop_io
  ! 13C
  use cable_c13o2_def,     only: c13o2_flux, c13o2_pool
  use cable_c13o2,         only: c13o2_write_restart_pools, c13o2_sanity_pools

  implicit none

  !!cln  character(len=99), intent(in)  :: fcnpspin
  real,    intent(in)    :: dels
  integer, intent(in)    :: kstart
  integer, intent(in)    :: kend
  integer, intent(in)    :: mloop
  type(veg_parameter_type),  intent(inout) :: veg  ! vegetation parameters
  type(soil_parameter_type), intent(inout) :: soil ! soil parameters
  type(casa_biome),          intent(inout) :: casabiome
  type(casa_pool),           intent(inout) :: casapool
  type(casa_flux),           intent(inout) :: casaflux
  type(casa_met),            intent(inout) :: casamet
  type(casa_balance),        intent(inout) :: casabal
  type(phen_variable),       intent(inout) :: phen
  type(POP_type),            intent(inout) :: POP
  type(climate_type),        intent(inout) :: climate
  ! 13C
  type(c13o2_flux),          intent(inout) :: c13o2flux
  type(c13o2_pool),          intent(inout) :: c13o2pools
  ! communicator for error-messages
  integer, intent(in)  :: icomm, ocomm

  ! local variables
  integer                  :: myearspin, nyear, nloop1
  character(len=99)        :: ncfile
  character(len=4)         :: cyear
  integer                  :: ktau, ktauday, nday, idoy, ktauy, nloop

  type(casa_met)      :: casametsave
  type(phen_variable) :: phensave

  call alloc_casa_var(casametsave, size(casamet%tairk, 1))
  call alloc_phenvariable(phensave, size(phen%phase, 1))
  call zero_casa_var(casametsave)
  call zero_phenvariable(phensave)

  ktauday = nint(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday
  ktau    = 0

  myearspin = cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1

  ! compute the mean fluxes and residence time of each carbon pool

  do nyear=1, myearspin
     write(cyear,fmt="(I4)") cable_user%casa_spin_startyear + nyear - 1
     ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
     ! 13C
     call read_casa_dump(ncfile, casamet, casaflux, phen, climate, c13o2flux, ktau, kend, .true.)

     do idoy=1, mdyear
        ktau = (idoy-1)*ktauday + 1
        casamet%tairk(:)   = casamet%Tairkspin(:,idoy)
        casamet%tsoil(:,1) = casamet%Tsoilspin_1(:,idoy)
        casamet%tsoil(:,2) = casamet%Tsoilspin_2(:,idoy)
        casamet%tsoil(:,3) = casamet%Tsoilspin_3(:,idoy)
        casamet%tsoil(:,4) = casamet%Tsoilspin_4(:,idoy)
        casamet%tsoil(:,5) = casamet%Tsoilspin_5(:,idoy)
        casamet%tsoil(:,6) = casamet%Tsoilspin_6(:,idoy)
        casamet%moist(:,1) = casamet%moistspin_1(:,idoy)
        casamet%moist(:,2) = casamet%moistspin_2(:,idoy)
        casamet%moist(:,3) = casamet%moistspin_3(:,idoy)
        casamet%moist(:,4) = casamet%moistspin_4(:,idoy)
        casamet%moist(:,5) = casamet%moistspin_5(:,idoy)
        casamet%moist(:,6) = casamet%moistspin_6(:,idoy)
        casaflux%cgpp(:)       = casamet%cgppspin(:,idoy)
        casaflux%crmplant(:,1) = casamet%crmplantspin_1(:,idoy)
        casaflux%crmplant(:,2) = casamet%crmplantspin_2(:,idoy)
        casaflux%crmplant(:,3) = casamet%crmplantspin_3(:,idoy)
        phen%phase(:)      = phen%phasespin(:,idoy)
        phen%doyphase(:,1) = phen%doyphasespin_1(:,idoy)
        phen%doyphase(:,2) = phen%doyphasespin_2(:,idoy)
        phen%doyphase(:,3) = phen%doyphasespin_3(:,idoy)
        phen%doyphase(:,4) = phen%doyphasespin_4(:,idoy)
        climate%qtemp_max_last_year(:) = real(casamet%mtempspin(:,idoy))
        climate%frec(:)                = real(casamet%frecspin(:,idoy))
        ! casaflux%Nmindep and casaflux%Pdep set in read_casa_dump
        ! 13C
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
        end if
        ! BLAZE
        if (cable_user%call_blaze) then
           climate%dprecip(:)      = real(casamet%dprecip_spin(:,idoy))
           climate%aprecip_av20(:) = real(casamet%aprecip_av20_spin(:,idoy))
           climate%du10_max(:)     = real(casamet%du10_max_spin(:,idoy))
           climate%drhum(:)        = real(casamet%drhum_spin(:,idoy))
           climate%dtemp_max(:)    = real(casamet%dtemp_max_spin(:,idoy))
           climate%dtemp_min(:)    = real(casamet%dtemp_max_spin(:,idoy))
           climate%KBDI(:)         = real(casamet%KBDI_spin(:,idoy))
           climate%D_MacArthur(:)  = real(casamet%D_MacArthur_spin(:,idoy))
           climate%FFDI(:)         = real(casamet%FFDI_spin(:,idoy))
           climate%DSLR(:)         = casamet%DSLR_spin(:,idoy)
           climate%last_precip(:)  = real(casamet%last_precip_spin(:,idoy))
        end if

        call master_send_input(icomm, casa_dump_ts, idoy)
      end do

  end do

  nloop1= max(1,mloop-3)

  do nloop=1, mloop
     write(*,*) 'nloop =', nloop
     !!CLN  OPEN(91,file=fcnpspin)
     !!CLN  read(91,*)
     do nyear=1, myearspin

        write(cyear,FMT="(I4)") cable_user%casa_spin_startyear + nyear - 1
        ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
        ! 13C
        call read_casa_dump(ncfile, casamet, casaflux, phen, climate, c13o2flux, ktau, kend, .true.)

        do idoy=1, mdyear
           ktauy = idoy*ktauday
           casamet%tairk(:)   = casamet%Tairkspin(:,idoy)
           casamet%tsoil(:,1) = casamet%Tsoilspin_1(:,idoy)
           casamet%tsoil(:,2) = casamet%Tsoilspin_2(:,idoy)
           casamet%tsoil(:,3) = casamet%Tsoilspin_3(:,idoy)
           casamet%tsoil(:,4) = casamet%Tsoilspin_4(:,idoy)
           casamet%tsoil(:,5) = casamet%Tsoilspin_5(:,idoy)
           casamet%tsoil(:,6) = casamet%Tsoilspin_6(:,idoy)
           casamet%moist(:,1) = casamet%moistspin_1(:,idoy)
           casamet%moist(:,2) = casamet%moistspin_2(:,idoy)
           casamet%moist(:,3) = casamet%moistspin_3(:,idoy)
           casamet%moist(:,4) = casamet%moistspin_4(:,idoy)
           casamet%moist(:,5) = casamet%moistspin_5(:,idoy)
           casamet%moist(:,6) = casamet%moistspin_6(:,idoy)
           casaflux%cgpp(:)       = casamet%cgppspin(:,idoy)
           casaflux%crmplant(:,1) = casamet%crmplantspin_1(:,idoy)
           casaflux%crmplant(:,2) = casamet%crmplantspin_2(:,idoy)
           casaflux%crmplant(:,3) = casamet%crmplantspin_3(:,idoy)
           phen%phase(:)      = phen%phasespin(:,idoy)
           phen%doyphase(:,1) = phen%doyphasespin_1(:,idoy)
           phen%doyphase(:,2) = phen%doyphasespin_2(:,idoy)
           phen%doyphase(:,3) = phen%doyphasespin_3(:,idoy)
           phen%doyphase(:,4) = phen%doyphasespin_4(:,idoy)
           climate%qtemp_max_last_year(:) = real(casamet%mtempspin(:,idoy))
           climate%frec(:)                = real(casamet%frecspin(:,idoy))
           ! casaflux%Nmindep and casaflux%Pdep set in read_casa_dump
           ! 13C
           if (cable_user%c13o2) then
              c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
              c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
           end if
           ! BLAZE
           if (cable_user%call_blaze) then
              climate%dprecip(:)      = real(casamet%dprecip_spin(:,idoy))
              climate%aprecip_av20(:) = real(casamet%aprecip_av20_spin(:,idoy))
              climate%du10_max(:)     = real(casamet%du10_max_spin(:,idoy))
              climate%drhum(:)        = real(casamet%drhum_spin(:,idoy))
              climate%dtemp_max(:)    = real(casamet%dtemp_max_spin(:,idoy))
              climate%dtemp_min(:)    = real(casamet%dtemp_max_spin(:,idoy))
              climate%KBDI(:)         = real(casamet%KBDI_spin(:,idoy))
              climate%D_MacArthur(:)  = real(casamet%D_MacArthur_spin(:,idoy))
              climate%FFDI(:)         = real(casamet%FFDI_spin(:,idoy))
              climate%DSLR(:)         = casamet%DSLR_spin(:,idoy)
              climate%last_precip(:)  = real(casamet%last_precip_spin(:,idoy))
           end if

           call master_send_input(icomm, casa_dump_ts, idoy)
        end do ! end doy

     end do   ! end of nyear

  end do     ! end of nloop

  ! save spin variables (for output) because returned as zero from workers
  call copy_save_variables(casamet, phen, casametsave, phensave)

  ! write(*,*) 'b4 master receive casa'
  call master_receive(ocomm, 0, casa_ts)
  ! write(*,*) 'after master receive casa'

  ! rewrite saved variables
  call copy_save_variables(casametsave, phensave, casamet, phen)

  ! call casa_poolout(ktau, veg, soil, casabiome, &
  !      casapool, casaflux, casamet, casabal, phen)
  call casa_fluxout(CABLE_USER%CASA_SPIN_STARTYEAR+myearspin-1, &
       veg, soil, casabal, casamet)

  ! call write_casa_restart_nc(casamet, casapool, casaflux, phen, .true.)
  call write_casa_restart_nc(casabiome, casamet, casapool, casaflux, &
       casabal, phen)
  ! 13C
  if (cable_user%c13o2) then
     call master_receive(ocomm, 0, c13o2_pool_ts)
     call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
     call c13o2_write_restart_pools(casamet, c13o2pools)
  end if

  if ( cable_user%call_POP .and. (POP%np.gt.0) ) then
     ! write(*,*) 'b4 master receive pop'
     call master_receive_pop(POP, ocomm)
     ! write(*,*) 'after master receive pop'
     call POP_io(pop, casamet, myearspin, 'WRITE_INI', .true.)
  end if

END SUBROUTINE master_spincasacnp

!*********************************************************************************************


SUBROUTINE master_CASAONLY_LUC(dels, kstart, kend, veg, casabiome, casapool, &
     casaflux, casamet, casabal, phen, POP, climate, LUC_EXPT, POPLUC, &
     ! 13C
     c13o2flux, c13o2pools, c13o2luc, &
     icomm, ocomm)

#ifndef __INTEL__
  use mpi, only: MPI_Send
#else
  use mpi
#endif
  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module,  ONLY: cable_user
  USE cable_IO_vars_module, ONLY: landpt, output
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  use casa_inout,           only: write_casa_restart_nc
  use casa_cable,           only: read_casa_dump
  USE POP_Types,            only: POP_TYPE
  USE POPMODULE,            ONLY: POP_init_single
  use cable_pop_io,         only: pop_io
  USE TypeDef,              ONLY: dp
  USE CABLE_LUC_EXPT,       ONLY: LUC_EXPT_TYPE, read_LUH2, &
       ptos, ptog, stog, gtos, pharv, smharv, syharv, &
       ptoc, ptoq, stoc, stoq, ctos, qtos
  USE POPLUC_Types
  USE POPLUC_Module,        ONLY: POPLUCStep, POPLUC_weights_Transfer, WRITE_LUC_OUTPUT_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, &
       POPLUC_set_patchfrac,  WRITE_LUC_OUTPUT_GRID_NC
  ! 13C
  use cable_c13o2_def,      only: c13o2_flux, c13o2_pool, c13o2_luc
  use cable_c13o2,          only: c13o2_save_luc, c13o2_update_luc, &
       c13o2_write_restart_pools, c13o2_write_restart_luc, &
       c13o2_sanity_pools, c13o2_sanity_luc
  use mo_utils,             only: eq

  IMPLICIT NONE

  real,                      intent(in)    :: dels
  integer,                   intent(in)    :: kstart
  integer,                   intent(in)    :: kend
  type(veg_parameter_type),  intent(inout) :: veg  ! vegetation parameters
  type(casa_biome),          intent(inout) :: casabiome
  type(casa_pool),           intent(inout) :: casapool
  type(casa_flux),           intent(inout) :: casaflux
  type(casa_met),            intent(inout) :: casamet
  type(casa_balance),        intent(inout) :: casabal
  type(phen_variable),       intent(inout) :: phen
  type(pop_type),            intent(inout) :: pop
  type(climate_type),        intent(inout) :: climate
  type(luc_expt_type),       intent(inout) :: luc_expt
  type(popluc_type),         intent(inout) :: popluc
  ! 13C
  type(c13o2_flux),          intent(inout) :: c13o2flux
  type(c13o2_pool),          intent(inout) :: c13o2pools
  type(c13o2_luc),           intent(inout) :: c13o2luc
  !type (casa_pool)   , intent(inout) :: sum_casapool
  !type (casa_flux)   , intent(inout) :: sum_casaflux
  ! communicator for error-messages
  integer,                   intent(in)    :: icomm, ocomm

  ! local variables
  integer           :: myearspin, nyear, yyyy, nyear_dump
  character(len=99) :: ncfile
  character(len=4)  :: cyear
  integer           :: ktau, ktauday, nday, idoy

  type(casa_met)      :: casametsave
  type(phen_variable) :: phensave

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  integer :: k, j, l
  integer :: rank, off, cnt, ierr
  ! 13C
  real(dp), dimension(:,:), allocatable :: casasave
  real(dp), dimension(:,:), allocatable :: lucsave


  !  if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
  !  if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
  !  if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
  !  if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
  !  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))
  !
  !  IF (cable_user%CALL_POP) THEN
  !     Iw = POP%Iwood
  !  END IF

  ktauday = int(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday
  !  ctime = 0
  !  CALL zero_sum_casa(sum_casapool, sum_casaflux)
  !       count_sum_casa = 0

  call alloc_casa_var(casametsave, size(casamet%tairk, 1))
  call alloc_phenvariable(phensave, size(phen%phase, 1))
  call zero_casa_var(casametsave)
  call zero_phenvariable(phensave)

  ! 13C
  if (cable_user%c13o2) then
     allocate(casasave(c13o2pools%ntile,c13o2pools%npools))
     allocate(lucsave(c13o2luc%nland,c13o2luc%npools))
  end if

  myearspin = cable_user%yearend - cable_user%yearstart + 1
  yyyy      = cable_user%yearstart - 1

  do nyear=1, myearspin
     yyyy = yyyy+1


     nyear_dump = mod(nyear, cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1)
     if (nyear_dump == 0) &
           nyear_dump = cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1

     write(cyear,fmt="(I4)") cable_user%casa_spin_startyear + nyear_dump - 1
     ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
     write(*,'(a,i04,a)') 'casaonly_LUC ', YYYY, ' '//trim(ncfile)
     ! 13C
     call read_casa_dump(trim(ncfile), casamet, casaflux, phen, climate, c13o2flux, 1, 1, .true.)

     !!CLN901  format(A99)
     do idoy=1, mdyear
        ktau = (idoy-1)*ktauday + ktauday

        casamet%tairk(:)   = casamet%Tairkspin(:,idoy)
        casamet%tsoil(:,1) = casamet%Tsoilspin_1(:,idoy)
        casamet%tsoil(:,2) = casamet%Tsoilspin_2(:,idoy)
        casamet%tsoil(:,3) = casamet%Tsoilspin_3(:,idoy)
        casamet%tsoil(:,4) = casamet%Tsoilspin_4(:,idoy)
        casamet%tsoil(:,5) = casamet%Tsoilspin_5(:,idoy)
        casamet%tsoil(:,6) = casamet%Tsoilspin_6(:,idoy)
        casamet%moist(:,1) = casamet%moistspin_1(:,idoy)
        casamet%moist(:,2) = casamet%moistspin_2(:,idoy)
        casamet%moist(:,3) = casamet%moistspin_3(:,idoy)
        casamet%moist(:,4) = casamet%moistspin_4(:,idoy)
        casamet%moist(:,5) = casamet%moistspin_5(:,idoy)
        casamet%moist(:,6) = casamet%moistspin_6(:,idoy)
        casaflux%cgpp(:)       = casamet%cgppspin(:,idoy)
        casaflux%crmplant(:,1) = casamet%crmplantspin_1(:,idoy)
        casaflux%crmplant(:,2) = casamet%crmplantspin_2(:,idoy)
        casaflux%crmplant(:,3) = casamet%crmplantspin_3(:,idoy)
        phen%phase(:)      = phen%phasespin(:,idoy)
        phen%doyphase(:,1) = phen%doyphasespin_1(:,idoy)
        phen%doyphase(:,2) = phen%doyphasespin_2(:,idoy)
        phen%doyphase(:,3) = phen%doyphasespin_3(:,idoy)
        phen%doyphase(:,4) = phen%doyphasespin_4(:,idoy)
        climate%qtemp_max_last_year(:) = real(casamet%mtempspin(:,idoy))
        climate%frec(:)                = real(casamet%frecspin(:,idoy))
        ! casaflux%Nmindep and casaflux%Pdep set in read_casa_dump
        ! 13C
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
        end if
        ! BLAZE
        if (cable_user%call_blaze) then
           climate%dprecip(:)      = real(casamet%dprecip_spin(:,idoy))
           climate%aprecip_av20(:) = real(casamet%aprecip_av20_spin(:,idoy))
           climate%du10_max(:)     = real(casamet%du10_max_spin(:,idoy))
           climate%drhum(:)        = real(casamet%drhum_spin(:,idoy))
           climate%dtemp_max(:)    = real(casamet%dtemp_max_spin(:,idoy))
           climate%dtemp_min(:)    = real(casamet%dtemp_max_spin(:,idoy))
           climate%KBDI(:)         = real(casamet%KBDI_spin(:,idoy))
           climate%D_MacArthur(:)  = real(casamet%D_MacArthur_spin(:,idoy))
           climate%FFDI(:)         = real(casamet%FFDI_spin(:,idoy))
           climate%DSLR(:)         = casamet%DSLR_spin(:,idoy)
           climate%last_precip(:)  = real(casamet%last_precip_spin(:,idoy))
        end if
        call master_send_input(icomm, casa_dump_ts, idoy)

        IF (idoy==mdyear) THEN ! end of year

           ! save spin variables (for output) because returned as zero from workers
           call copy_save_variables(casamet, phen, casametsave, phensave)

           ! get casa update from workers
           call master_receive(ocomm, 0, casa_ts)
           CALL master_receive(ocomm, 0, casa_LUC_ts)

           ! rewrite saved variables
           call copy_save_variables(casametsave, phensave, casamet, phen)

           ! 13C
           if (cable_user%c13o2) then
              call master_receive(ocomm, 0, c13o2_flux_ts)
              call master_receive(ocomm, 0, c13o2_pool_ts)
              call master_receive(ocomm, 0, c13o2_luc_ts)
           end if

           LUC_EXPT%CTSTEP = yyyy -  LUC_EXPT%FirstYear + 1

           CALL READ_LUH2(LUC_EXPT)

           DO k=1,mland
              POPLUC%ptos(k)   = real(LUC_EXPT%INPUT(ptos)%VAL(k), dp)
              POPLUC%ptog(k)   = real(LUC_EXPT%INPUT(ptog)%VAL(k), dp)
              POPLUC%stog(k)   = real(LUC_EXPT%INPUT(stog)%VAL(k), dp)
              POPLUC%gtop(k)   = 0.0_dp
              POPLUC%gtos(k)   = real(LUC_EXPT%INPUT(gtos)%VAL(k), dp)
              POPLUC%pharv(k)  = real(LUC_EXPT%INPUT(pharv)%VAL(k), dp)
              POPLUC%smharv(k) = real(LUC_EXPT%INPUT(smharv)%VAL(k), dp)
              POPLUC%syharv(k) = real(LUC_EXPT%INPUT(syharv)%VAL(k), dp)

              POPLUC%ptoc(k) = real(LUC_EXPT%INPUT(ptoc)%VAL(k), dp)
              POPLUC%ptoq(k) = real(LUC_EXPT%INPUT(ptoq)%VAL(k), dp)
              POPLUC%stoc(k) = real(LUC_EXPT%INPUT(stoc)%VAL(k), dp)
              POPLUC%stoq(k) = real(LUC_EXPT%INPUT(stoq)%VAL(k), dp)
              POPLUC%ctos(k) = real(LUC_EXPT%INPUT(ctos)%VAL(k), dp)
              POPLUC%qtos(k) = real(LUC_EXPT%INPUT(qtos)%VAL(k), dp)

              POPLUC%thisyear = yyyy
           END DO

           ! set landuse index for secondary forest POP landscapes
           DO k=1, POP%np
              IF (yyyy.eq.LUC_EXPT%YearStart) THEN
                 if (veg%iLU(POP%Iwood(k)).eq.2) then
                    !POP%LU(k) = 2
                    POP%pop_grid(k)%LU = 2
                 end if
              end if
           END DO

           ! zero secondary forest tiles in POP where secondary forest area is zero
           DO k=1,mland
              if ( eq(POPLUC%primf(k)-POPLUC%frac_forest(k), 0.0_dp) &
                   .and. (.not. LUC_EXPT%prim_only(k)) ) then

                 j = landpt(k)%cstart+1
                 do l=1,size(POP%Iwood)
                    if( POP%Iwood(l) == j) then
                       CALL POP_init_single(POP,veg%disturbance_interval,l)
                       exit
                    end if
                 end do

                 casapool%cplant(j,leaf) = 0.01_dp
                 casapool%nplant(j,leaf)= casabiome%ratioNCplantmin(veg%iveg(j),leaf) * casapool%cplant(j,leaf)
                 casapool%pplant(j,leaf)= casabiome%ratioPCplantmin(veg%iveg(j),leaf) * casapool%cplant(j,leaf)

                 casapool%cplant(j,froot) = 0.01_dp
                 casapool%nplant(j,froot)= casabiome%ratioNCplantmin(veg%iveg(j),froot) * casapool%cplant(j,froot)
                 casapool%pplant(j,froot)= casabiome%ratioPCplantmin(veg%iveg(j),froot) * casapool%cplant(j,froot)

                 casapool%cplant(j,wood) = 0.01_dp
                 casapool%nplant(j,wood)= casabiome%ratioNCplantmin(veg%iveg(j),wood) * casapool%cplant(j,wood)
                 casapool%pplant(j,wood)= casabiome%ratioPCplantmin(veg%iveg(j),wood) * casapool%cplant(j,wood)
                 casaflux%frac_sapwood(j) = 1.0_dp

                 ! 13C
                 if (cable_user%c13o2) then
                    c13o2pools%cplant(j,leaf)  = 0.01_dp ! * vpdbc13 / vpdbc13 ! Divide by 13C
                    c13o2pools%cplant(j,wood)  = 0.01_dp ! * vpdbc13 / vpdbc13 ! so that about same numerical precision as 12C
                    c13o2pools%cplant(j,froot) = 0.01_dp ! * vpdbc13 / vpdbc13 !
                 end if
              end if
           END DO

           CALL POPLUCStep(POPLUC,yyyy)

           CALL POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

           ! transfer POP updates to workers
           off = 1
           DO rank=1, wnp
              IF (rank .GT. 1) off = off + wland(rank-1)%npop_iwood
              cnt = wland(rank)%npop_iwood
              CALL MPI_Send(POP%pop_grid(off), cnt, pop_ts, rank, 0, icomm, ierr)
           END DO

           ! workers call POP here

           ! CALL POPdriver(casaflux,casabal,veg, POP)
           CALL master_receive_pop(POP, ocomm)

           ! 13C
           if (cable_user%c13o2) call c13o2_save_luc(casapool, popluc, casasave, lucsave)
           CALL POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
           ! 13C
           if (cable_user%c13o2) then
              call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc)
              call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
              call c13o2_sanity_luc(popluc, c13o2luc)
           end if

           IF (output%grid(1:3) == 'lan') THEN
              CALL WRITE_LUC_OUTPUT_NC(POPLUC, YYYY, (YYYY.EQ.cable_user%YearEnd))
           ELSE
              CALL WRITE_LUC_OUTPUT_GRID_NC(POPLUC, YYYY, (YYYY.EQ.cable_user%YearEnd))
           END IF

           CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)

        END IF  ! end of year

     end do
     ! send updates for CASA pools, resulting from LUC
     !NEW CALL master_send_input(icomm, casa_ts, nyear)
     CALL master_send_input(icomm, casa_LUC_ts, nyear)
     ! 13C
     if (cable_user%c13o2) then
        call master_send_input(icomm, c13o2_flux_ts, nyear)
        call master_send_input(icomm, c13o2_pool_ts, nyear)
        call master_send_input(icomm, c13o2_luc_ts, nyear)
     end if
  end do ! year=1,myearspin

  ! 13C
  if (cable_user%c13o2) then
     call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
     if (cable_user%POPLUC) call c13o2_sanity_luc(popluc, c13o2luc)
  end if

  CALL WRITE_LUC_RESTART_NC(POPLUC)
  ! CALL write_casa_restart_nc(casamet, casapool, casaflux, phen, .TRUE.)
  CALL write_casa_restart_nc(casabiome, casamet, casapool, casaflux, casabal, phen)
  ! 13C
  if (cable_user%c13o2) then
     call c13o2_write_restart_pools(casamet, c13o2pools)
     if (cable_user%POPLUC) call c13o2_write_restart_luc(popluc, c13o2luc)
  end if

  CALL POP_IO(pop, casamet, myearspin, 'WRITE_INI', .TRUE.)

END SUBROUTINE master_CASAONLY_LUC


!********************************************************************************************************
! subroutine for reading LU input data, zeroing biomass in empty secondary forest tiles
! and tranferring LUC-based age weights for secondary forest to POP structure


SUBROUTINE LUCdriver(casabiome, casapool, casaflux, POP, LUC_EXPT, POPLUC, veg, &
     ! 13C
     c13o2pools)

  USE cable_def_types_mod , ONLY: r_2, veg_parameter_type, mland
  USE cable_carbon_module
  USE cable_common_module,  ONLY: cable_user, CurYear
  USE cable_IO_vars_module, ONLY: landpt
  USE casadimension
  USE casaparm
  USE casavariable
  USE POP_Types,            ONLY: POP_TYPE
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

  write(*,*) 'cablecasa_LUC', CurYear
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

     ! MC - not in serial code
     ! POPLUC%ptoc(k) = real(LUC_EXPT%INPUT(ptoc)%VAL(k), r_2)
     ! POPLUC%ptoq(k) = real(LUC_EXPT%INPUT(ptoq)%VAL(k), r_2)
     ! POPLUC%stoc(k) = real(LUC_EXPT%INPUT(stoc)%VAL(k), r_2)
     ! POPLUC%stoq(k) = real(LUC_EXPT%INPUT(stoq)%VAL(k), r_2)
     ! POPLUC%ctos(k) = real(LUC_EXPT%INPUT(ctos)%VAL(k), r_2)
     ! POPLUC%qtos(k) = real(LUC_EXPT%INPUT(qtos)%VAL(k), r_2)

     POPLUC%thisyear = yyyy
  END DO
  ! zero secondary forest tiles in POP where secondary forest area is zero
  DO k=1,mland
     if ( eq(POPLUC%frac_primf(k)-POPLUC%frac_forest(k), 0.0_r_2) &
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

  CALL POPLUCStep(POPLUC,yyyy)

  CALL POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

END SUBROUTINE LUCdriver

!*********************************************************************************************************

END MODULE cable_mpimaster
