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
! Uses:           mpi
!                 cable_mpicommon
!                 cable_def_types_mod
!                 cable_IO_vars_module
!                 cable_common_module
!                 cable_data_module
!                 cable_input_module
!                 cable_output_module
!                 cable_cbm_module
!                 casadimension
!                 casavariable
!                 phenvariable
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

  USE cable_mpicommon

  IMPLICIT NONE

  SAVE

  PRIVATE

  ! number of workers; set in master_decomp
  INTEGER :: wnp

  ! TODO: m3d_t mat_t and vec_t to be factored out from here and from master_outtypes
  ! MPI: slices of 3D arrays
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: m3d_t
  ! MPI: slices of matrices (2D)
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mat_t

  ! MPI: parts of vectors (1D)
  ! MPI: vec_t dimension is wnp; as each worker gets a single indexed with nvec blocks
  INTEGER, ALLOCATABLE, DIMENSION(:) :: vec_t

  ! MPI derived datatype handles for sending input data to the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: inp_ts

  ! MPI derived datatype handles for receiving output from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_ts

  ! master's struct for receiving restart data from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: restart_ts

  ! CASA related derived types

  ! MPI derived datatype handles for receiving casa results from the workers
  ! and restart values
  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_ts

  ! MPI derived datatype handles for send/receiving casa dump values from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_dump_ts

  ! MPI derived datatype handles for send/receiving casa pool values (needed for LUC)
  !  from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_LUC_ts

  !CLN  ! MPI derived datatype handles for receiving casa restart values from the workers
  !CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_restart_ts

  ! climate derived type
  INTEGER, ALLOCATABLE, DIMENSION(:) :: climate_ts

  ! MPI derived datatype handles for Sending/receiving vals results for BLAZE
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blaze_out_ts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blaze_restart_ts

  ! POP related derived types

  ! MPI derived datatype handles for receiving POP results from the workers
  INTEGER :: pop_ts

  ! 13C
  ! MPI derived datatype handles for receiving c13o2 results from the workers
  ! and restart values
  INTEGER, ALLOCATABLE, DIMENSION(:) :: c13o2_flux_ts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: c13o2_pool_ts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: c13o2_luc_ts

  ! MPI: isend request array for scattering input data to the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: inp_req
  ! MPI: isend status array for scattering input data to the workers
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: inp_stats

  ! MPI: irecv request array for gathering results from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_req
  ! MPI: irecv status array for gathering results from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: recv_stats

  ! MPI: landpoints decomposition; global info used by the master process
  TYPE(lpdecomp_t), ALLOCATABLE, DIMENSION(:) :: wland

  PUBLIC :: mpidrv_master

CONTAINS

  SUBROUTINE mpidrv_master(comm)

    use mpi

    USE cable_def_types_mod
    USE cable_io_vars_module, ONLY: logn, gswpfile, ncciy, leaps, &
         verbose, fixedCO2, output, check, patchout, soilparmnew, &
         timeunits, exists, calendar, landpt, latitude, longitude
    USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user, &
         cable_runtime, filename, &
         redistrb, wiltParam, satuParam, CurYear, &
         IS_LEAPYEAR, IS_CASA_TIME, calcsoilalbedo, get_unit, &
         report_version_no, kwidth_gl
    USE cable_data_module,    ONLY: driver_type, icanopy_type, point2constants
    USE cable_input_module,   ONLY: open_met_file, load_parameters, &
         get_met_data,close_met_file
    USE cable_output_module,  ONLY: create_restart, open_output_file, &
         write_output, close_output_file
    USE cable_write_module,   ONLY: nullify_write
    USE cable_cbm_module
    USE cable_climate_mod

    ! modules related to CASA-CNP
    USE casadimension,        ONLY: icycle
    USE casavariable,         ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
         casa_met, casa_balance, zero_sum_casa, update_sum_casa
    USE phenvariable,         ONLY: phen_variable

    !CLN added
    ! modules related to POP
    USE POP_Types,            ONLY: POP_TYPE
    USE POPLUC_Types,         ONLY: POPLUC_Type
    USE POPLUC_Module,        ONLY: WRITE_LUC_OUTPUT_NC, WRITE_LUC_OUTPUT_GRID_NC, &
         POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, POPLUC_set_patchfrac, &
         READ_LUC_RESTART_NC, alloc_popluc
    ! LUC_EXPT only
    USE CABLE_LUC_EXPT,       ONLY: LUC_EXPT_TYPE, LUC_EXPT_INIT, close_luh2

    ! modules related to fire
    USE BLAZE_MOD,            ONLY: TYPE_BLAZE, INI_BLAZE, WRITE_BLAZE_OUTPUT_NC
    USE BLAZE_MPI,            ONLY: MASTER_BLAZE_TYPES, MASTER_SIMFIRE_TYPES
    USE SIMFIRE_MOD,          ONLY: TYPE_SIMFIRE, INI_SIMFIRE

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

    ! BIOS only
    USE cable_bios_met_obs_params,   ONLY:  cable_bios_read_met, cable_bios_init, &
         cable_bios_load_params, cable_bios_load_climate_params

    IMPLICIT NONE

    ! MPI:
    INTEGER               :: comm ! MPI communicator for comms with the workers

    ! CABLE namelist: model configuration, runtime/user switches
    CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml'

    ! timing variables
    INTEGER, PARAMETER ::  kstart = 1   ! start of simulation

    INTEGER ::        &
         ktau,        &  ! increment equates to timestep, resets if spinning up
         ktau_tot,    &  ! NO reset when spinning up, total timesteps by model
         kend,        &  ! no. of time steps in run
                                !CLN      kstart = 1, &  ! timestep to start at
         koffset = 0, &  ! timestep to start at
         ktauday,     &  ! day counter for CASA-CNP
         idoy,        &  ! day of year (1:365) counter for CASA-CNP
         nyear,       &  ! year counter for CASA-CNP
         ctime,       &  ! time for casacnp
         YYYY,        &  !
         LOY,         &  ! Length of Year
         maxdiff(2)      ! location of maximum in convergence test

    REAL :: dels                    ! time step size in seconds
    CHARACTER(len=9) :: dum         ! dummy char for fileName generation
    CHARACTER(len=4) :: str1        ! dummy char for fileName generation

    ! CABLE variables
    TYPE(met_type)       :: met     ! met input variables: see below for imet in MPI variables
    TYPE(air_type)       :: air     ! air property variables
    TYPE(canopy_type)    :: canopy  ! vegetation variables
    TYPE(radiation_type) :: rad     ! radiation variables
    TYPE(roughness_type) :: rough   ! roughness varibles
    TYPE(balances_type)  :: bal     ! energy and water balance variables
    TYPE(soil_snow_type) :: ssnow   ! soil and snow variables
    TYPE(climate_type)   :: climate ! climate variables

    ! CABLE parameters
    TYPE(soil_parameter_type) :: soil     ! soil parameters
    TYPE(veg_parameter_type)  :: veg      ! vegetation parameters: see below for iveg in MPI variables
    TYPE(driver_type)         :: C        ! constants used locally
    TYPE(sum_flux_type)       :: sum_flux ! cumulative flux variables
    TYPE(bgc_pool_type)       :: bgc      ! carbon pool variables

    ! CASA-CNP variables
    TYPE(casa_biome)     :: casabiome
    TYPE(casa_pool)      :: casapool
    TYPE(casa_flux)      :: casaflux
    TYPE(casa_pool)      :: sum_casapool
    TYPE(casa_flux)      :: sum_casaflux
    TYPE(casa_met)       :: casamet
    TYPE(casa_balance)   :: casabal
    TYPE(phen_variable)  :: phen
    TYPE(POP_TYPE)       :: POP
    TYPE(POPLUC_TYPE)    :: POPLUC
    TYPE(LUC_EXPT_TYPE)  :: LUC_EXPT
    TYPE(PLUME_MIP_TYPE) :: PLUME
    TYPE(CRU_TYPE)       :: CRU
    CHARACTER            :: cyear*4
    CHARACTER            :: ncfile*99

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
    ! delta-13C of atmospheric CO2
    integer            :: iunit, ios
    real               :: iyear
    integer            :: c13o2_atm_syear, c13o2_atm_eyear
    character(len=100) :: header

    ! declare vars for switches (default .FALSE.) etc declared thru namelist
    LOGICAL, SAVE           :: &
         vegparmnew    = .FALSE., & ! using new format input file (BP dec 2007)
         spinup        = .FALSE., & ! model spinup to soil state equilibrium?
         spinConv      = .FALSE., & ! has spinup converged?
         spincasainput = .FALSE., & ! TRUE: SAVE input req'd to spin CASA-CNP;
                                    ! FALSE: READ input to spin CASA-CNP
         spincasa      = .FALSE., & ! TRUE: CASA-CNP Will spin mloop,
                                    ! FALSE: no spin up
         l_casacnp     = .FALSE., & ! using CASA-CNP with CABLE
         l_laiFeedbk   = .FALSE., & ! using prognostic LAI
         l_vcmaxFeedbk = .FALSE., & ! using prognostic Vcmax
         CASAONLY      = .FALSE., & ! ONLY Run CASA-CNP
         CALL1         = .TRUE.

         REAL              :: &
         delsoilM,         & ! allowed variation in soil moisture for spin up
         delsoilT            ! allowed variation in soil temperature for spin up

    ! temporary storage for soil moisture/temp. in spin up mode
    REAL, ALLOCATABLE, DIMENSION(:,:) :: &
         soilMtemp, &
         soilTtemp

    ! MPI:
    TYPE(met_type)           :: imet  ! read ahead met input variables
    TYPE(veg_parameter_type) :: iveg  ! MPI read ahead vegetation parameters
    LOGICAL :: loop_exit     ! MPI: exit flag for bcast to workers
    INTEGER :: iktau    ! read ahead index of time step = 1 ..  kend
    INTEGER :: oktau    ! ktau = 1 ..  kend for output
    INTEGER :: icomm ! separate dupes of MPI communicator for send and recv
    INTEGER :: ocomm ! separate dupes of MPI communicator for send and recv
    INTEGER :: ierr
    INTEGER :: rank, off, cnt

    ! Vars for standard for quasi-bitwise reproducability b/n runs
    ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
    CHARACTER(len=30), PARAMETER :: &
         Ftrunk_sumbal  = ".trunk_sumbal", &
         Fnew_sumbal    = "new_sumbal"

    REAL(r_2), SAVE :: &
         trunk_sumbal = 0.0, & !
         new_sumbal   = 0.0, &
         new_sumfpn   = 0.0, &
         new_sumfe    = 0.0

    INTEGER :: count_bal = 0
    INTEGER :: nkend=0
    INTEGER :: ioerror=0


    ! switches etc defined thru namelist (by default cable.nml)
    NAMELIST/CABLE/                  &
         filename,         & ! TYPE, containing input filenames
         vegparmnew,       & ! use new soil param. method
         soilparmnew,      & ! use new soil param. method
         calcsoilalbedo,   & ! vars intro for Ticket #27
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
    
    integer :: kk
    integer :: lalloc
    integer, parameter :: mloop = 30 !MCTEST 30   ! CASA-CNP PreSpinup loops
    real :: etime, etimelast

    ! command line arguments
    integer :: narg, len1, len2
    character(len=500) :: arg1
    character(len=200) :: arg2

    ! end header
    
    
    etimelast = 0.0
    
    ! Open, read and close the namelist file.
    open(10, file = cable_namelist, status="old", action="read")
    read(10, nml=cable )   !where nml=cable defined above
    close(10)
    
    ! Open, read and close the consistency check file.
    ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
    if (cable_user%consistency_check) then
       open(11, file=ftrunk_sumbal, status='old', action='read', iostat=ioerror)
       if (ioerror==0) then
          read(11,*) trunk_sumbal  ! written by previous trunk version
       endif
       close(11)
    endif

    ! Open log file:
    open(logn,file=filename%log)

    call report_version_no(logn)

    ! if (iargc() > 0) then
    !    call getarg(1, filename%met)
    !    call getarg(2, casafile%cnpipool)
    ! endif
    narg = command_argument_count()
    if (narg > 0) then
       call get_command_argument(1, arg1, len1)
       filename%met = arg1(1:len1)
       call get_command_argument(2, arg2, len2)
       casafile%cnpipool = arg2(1:len2)
    endif

    ! INITIALISATION depending on nml settings
    if (trim(cable_user%MetType) .EQ. 'gswp') THEN
       if (cable_user%YearStart.eq.0 .and. ncciy.gt.0) then
          cable_user%YearStart = ncciy
          cable_user%YearEnd   = ncciy
       elseif (cable_user%YearStart.eq.0 .and. ncciy.eq.0) then
          write(*,*) 'undefined start year for gswp met: '
          write(*,*) 'enter value for ncciy or'
          write(*,*) '(cable_user%YearStart and  cable_user%YearEnd) in cable.nml'
          write(logn,*) 'undefined start year for gswp met:'
          write(logn,*) 'enter value for ncciy or'
          write(logn,*) '(cable_user%YearStart and  cable_user%YearEnd) in cable.nml'
          stop
       endif
    endif

    CurYear = cable_user%YearStart

    if (icycle .ge. 11) then
       icycle                     = icycle - 10
       casaonly                   = .true.
       cable_user%casa_dump_read  = .true.
       cable_user%casa_dump_write = .false.
    elseif (icycle .eq. 0) then
       cable_user%casa_dump_read  = .false.
       cable_user%call_pop        = .false.
       cable_user%call_blaze      = .false.
    endif

    !! vh_js !!
    if (icycle.gt.0) then
       l_casacnp = .true.
    else
       l_casacnp = .false.
    endif

    !! vh_js !! suggest LALLOC should ulitmately be a switch in the .nml file
    if (cable_user%call_pop) then
       lalloc = 3 ! for use with POP: makes use of pipe model to partition between stem and leaf
    else
       lalloc = 0 ! default
    endif

    if ( trim(cable_user%MetType) .eq. 'gpgs' ) then
       leaps = .true.
       cable_user%MetType = 'gswp'
    elseif ( trim(cable_user%MetType) .eq. 'bios' ) then
       leaps = .true.
    endif

    cable_runtime%offline = .true.

    ! associate pointers used locally with global definitions
    call point2constants( C )

    if ( l_casacnp  .and. (icycle == 0 .or. icycle > 3) ) &
         stop 'icycle must be 1 to 3 when using casaCNP'
    if ( (l_laifeedbk .or. l_vcmaxfeedbk) .and. (.not. l_casacnp) ) &
         stop 'casaCNP required to get prognostic LAI or Vcmax'
    if (l_vcmaxfeedbk .and. icycle < 2) &
         stop 'icycle must be 2 to 3 to get prognostic Vcmax'
    if ( icycle > 0 .and. (.not. soilparmnew) ) &
         stop 'casaCNP must use new soil parameters'

    ! casa time count
    ctime = 0

    ! Iinitialise settings depending on met dataset

    ! Open met data and get site information from netcdf file. (NON-GSWP ONLY!)
    ! This retrieves time step size, number of timesteps, starting date,
    ! latitudes, longitudes, number of sites.
    IF (TRIM(cable_user%MetType) .NE. "gswp" .AND. &
        TRIM(cable_user%MetType) .NE. "gpgs" .AND. &
        TRIM(cable_user%MetType) .NE. "plum" .AND. &
        TRIM(cable_user%MetType) .NE. "bios" .AND. &
        TRIM(cable_user%MetType) .NE. "cru") THEN
       CALL open_met_file( dels, koffset, kend, spinup, C%TFRZ )
       IF ( koffset .NE. 0 .AND. cable_user%CALL_POP ) THEN
          WRITE(*,*) "When using POP, episode must start on Jan 1st!"
          STOP 991
       ENDIF
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

    ! Tell the workers if we're leaping
    ! print*, 'MASTER Send 01 ', leaps
    call MPI_Bcast(leaps, 1, MPI_LOGICAL, 0, comm, ierr)

    ! outer loop - spinup loop no. ktau_tot :
    ktau_tot = 0
    ktau     = 0
    SPINLOOP: do

       ! Bios initialisation
       if (trim(cable_user%MetType) .eq. "bios") then
          call cpu_time(etime)
          call cable_bios_init(dels, curyear, met, kend, ktauday)
          call MPI_Bcast(mland, 1, MPI_INTEGER, 0, comm, ierr)
          call MPI_Bcast(latitude, mland, MPI_REAL, 0, comm, ierr)
          call MPI_Bcast(longitude, mland, MPI_REAL, 0, comm, ierr)
          koffset = 0
          leaps   = .true.
          write(str1,'(i4)') curyear
          str1 = adjustl(str1)
          timeunits="seconds since "//trim(str1)//"-01-01 00:00:00"
          calendar = 'standard'
       endif

       ! Loop through simulation years
       YEARLOOP: do YYYY=cable_user%YearStart, cable_user%YearEnd
          CurYear = YYYY
          if (leaps .and. is_leapyear(YYYY)) then
             LOY = 366
          else
             LOY = 365
          endif

          ! CLN From here extract CALL1 and put LOY computation to end of CALL1 block
          if (trim(cable_user%MetType) .eq. 'plum') then
             if ( CALL1 ) THEN
                call plume_mip_init( plume )
                dels    = plume%dt
                koffset = 0
                leaps   = plume%LeapYears
                write(str1,'(i4)') CurYear
                str1 = adjustl(str1)
                timeunits="seconds since "//trim(str1)//"-01-01 00:00:00"
             endif
             kend = nint(24.0*3600.0/dels) * LOY
          else if (trim(cable_user%MetType) .eq. 'bios') then
             kend = nint(24.0*3600.0/dels) * LOY
          else if (trim(cable_user%MetType) .eq. 'cru') then
             if (CALL1) then
                call cpu_time(etime)
                call cru_init( cru )
                dels         = cru%dtsecs
                koffset      = 0
                leaps        = .false.  ! No leap years in CRU-NCEP
                exists%Snowf = .false.  ! No snow in CRU-NCEP, so ensure it will
                                        ! be determined from temperature in CABLE
                write(str1,'(i4)') CurYear
                str1 = adjustl(str1)
                timeunits="seconds since "//trim(str1)//"-01-01 00:00:00"
                calendar = "standard"
             endif
             LOY = 365
             kend = nint(24.0*3600.0/dels) * LOY
          else if (trim(cable_user%MetType) .eq. 'gswp') then
             ncciy = CurYear
             write(*,*) 'Looking for global offline run info.'
             call preparefiles(ncciy)
             call open_met_file( dels, koffset, kend, spinup, c%tfrz )
          endif ! cable_user%MetType
          ! CLN To here extract CALL1 and put LOY computation to end of CALL1 block

          ! some things (e.g. CASA-CNP) only need to be done once per day
          ktauday = int(24.0*3600.0/dels)

          ! Checks where parameters and initialisations should be loaded from.
          ! If they can be found in either the met file or restart file, they will
          ! load from there, with the met file taking precedence. Otherwise, they'll
          ! be chosen from a coarse global grid of veg and soil types, based on
          ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
          IF ( CALL1 ) THEN
             IF (cable_user%POPLUC) THEN
                CALL LUC_EXPT_INIT(LUC_EXPT)
                call alloc_POPLUC(POPLUC, mland)
                POPLUC%np = mland
                !MC - Why in MPI but not in cable_driver?
                !     READ_LUC_RESTART_NC is called in POPLUC_init, which is called in load_parameters
                IF (trim(cable_user%POPLUC_RunType) .EQ. 'restart') THEN
                   CALL READ_LUC_RESTART_NC(POPLUC)
                   LUC_EXPT%grass    = real(POPLUC%grass)
                   LUC_EXPT%primaryf = real(min(POPLUC%primf, 1.0_r_2-POPLUC%grass))
                   LUC_EXPT%secdf    = max((1.0 -  LUC_EXPT%grass - LUC_EXPT%primaryf), 0.0)
                   LUC_EXPT%crop     = real(max(min(POPLUC%crop, POPLUC%grass), 0.0_r_2))
                   LUC_EXPT%past     = max(min(LUC_EXPT%grass - LUC_EXPT%crop, real(POPLUC%past)), 0.0)
                ENDIF
             ENDIF

             !! vh_js !!
             CALL load_parameters(met, air, ssnow, veg, climate, bgc, &
                  soil, canopy, rough, rad, sum_flux, &
                  bal, logn, vegparmnew, casabiome, casapool, &
                  casaflux, sum_casapool, sum_casaflux, &
                  casamet, casabal, phen, POP, spinup, &
                  C%EMSOIL, C%TFRZ, LUC_EXPT, POPLUC, BLAZE, SIMFIRE, &
                  c13o2flux, c13o2pools, sum_c13o2pools, c13o2luc)
             ! print*, 'CD01   ', casamet%glai
             ! 13C
             if (cable_user%c13o2) then
                allocate(casasave(c13o2pools%ntile,c13o2pools%npools))
                if (cable_user%popluc) allocate(lucsave(c13o2luc%nland,c13o2luc%npools))
             endif

             !par disabled as blaze init moved below
             ! Abort, if an error occurred during BLAZE/SIMFIRE init
             !IF (BLAZE%ERR) CALL MPI_Abort(comm,0,ierr)

             IF (cable_user%POPLUC .AND. TRIM(cable_user%POPLUC_RunType) .EQ. 'static') &
                  cable_user%POPLUC= .FALSE.

             !TRUNK not in trunk
             if (cable_user%call_climate) then
                CALL alloc_cbm_var(climate, mp, ktauday)
                CALL climate_init(climate, mp, ktauday)
                if (.NOT. cable_user%climate_fromzero) &
                     CALL READ_CLIMATE_RESTART_NC(climate, ktauday)
             endif

             ! Having read the default parameters, if this is a bios run we will now
             ! overwrite the subset of them required for bios.
             IF ( TRIM(cable_user%MetType) .EQ. 'bios' ) THEN
                CALL cable_bios_load_params(soil)
             ENDIF

             ssnow%otss_0   = ssnow%tgg(:,1)
             ssnow%otss     = ssnow%tgg(:,1)
             ssnow%tss      = ssnow%tgg(:,1)
             canopy%fes_cor = 0.
             canopy%fhs_cor = 0.
             met%ofsd       = 0.1

             !! CLN BLAZE
             !! IF ( cable_user%CALL_BLAZE ) &
              ! additional params needed for BLAZE
             if ( trim(cable_user%MetType) .eq. 'bios' ) call cable_bios_load_climate_params(climate)


             IF (.NOT. spinup) spinConv = .TRUE.

             ! MPI: above was standard serial code
             ! now it's time to initialize the workers

             ! MPI: bcast to workers so that they don't need to open the met file themselves
             ! print*, 'MASTER Send 02 ', dels
             CALL MPI_Bcast(dels, 1, MPI_REAL, 0, comm, ierr)
          ENDIF ! Call1

          ! print*, 'MASTER Send 03 ', kend
          CALL MPI_Bcast(kend, 1, MPI_INTEGER, 0, comm, ierr)
          ! print*, 'MASTER Sent 03 '

          IF ( CALL1 ) THEN
             ! MPI: need to know extents before creating datatypes
             CALL find_extents()

             ! MPI: calculate and broadcast landpoint decomposition to the workers
             ! print*, 'MASTER Send 04 comm'
             CALL master_decomp(comm, mland, mp)

             ! MPI: set up stuff for new irecv isend code that separates completion
             ! from posting of requests
             ! wnp is set in master_decomp above
             ALLOCATE(inp_req(wnp))
             ALLOCATE(inp_stats(MPI_STATUS_SIZE, wnp))
             ALLOCATE(recv_req(wnp))
             ALLOCATE(recv_stats(MPI_STATUS_SIZE, wnp))
             ! print*, 'MASTER Send 05 icomm'
             CALL MPI_Comm_dup(comm, icomm, ierr)
             ! print*, 'MASTER Send 06 ocomm'
             CALL MPI_Comm_dup(comm, ocomm, ierr)

             ! MPI: data set in load_parameter is now scattered out to the workers
             ! print*, 'MASTER Send 07 params'
             CALL master_cable_params(comm, met, air, ssnow, veg, bgc, soil, canopy, &
                  rough, rad, sum_flux, bal)

             ! print*, 'MASTER Send 08 climate'
             !TRUNK IF (cable_user%call_climate) &
             CALL master_climate_types(comm, climate, ktauday)

             ! MPI: mvtype and mstype send out here instead of inside master_casa_params
             !      so that old CABLE carbon module can use them. (BP May 2013)
             ! print*, 'MASTER Send 09 ', mvtype
             CALL MPI_Bcast(mvtype, 1, MPI_INTEGER, 0, comm, ierr)
             ! print*, 'MASTER Send 10 ', mstype
             CALL MPI_Bcast(mstype, 1, MPI_INTEGER, 0, comm, ierr)

             ! MPI: casa parameters scattered only if cnp module is active
             IF (icycle>0) THEN
                ! MPI:
                ! print*, 'MASTER Send 11 casa'
                CALL master_casa_params(comm, casabiome, casapool, casaflux, casamet, casabal, phen)
                ! 13C
                if (cable_user%c13o2) then
                   ! print*, 'MASTER Send 12 c13o2_flux'
                   CALL master_c13o2_flux_params(comm, c13o2flux)
                   ! print*, 'MASTER Send 13 c13o2_pool'
                   CALL master_c13o2_pool_params(comm, c13o2pools)
                   if (cable_user%popluc) then
                      ! print*, 'MASTER Send 13.1 c13o2_luc'
                      call master_c13o2_luc_params(comm, c13o2luc)
                   endif
                endif

                ! Create and Send POP ini/restart data
                if (cable_user%call_pop) then
                   ! print*, 'MASTER Send 14 pop'
                   call master_pop_types(comm, casamet, pop)
                endif

                ! Fire init and
                if (cable_user%call_blaze) then
                   call ini_blaze( mland, rad%latitude(landpt(:)%cstart), &
                                   rad%longitude(landpt(:)%cstart), blaze )

                   !par blaze restart not required uses climate data
                   !create handles for restart-data
                   ! print*, 'MASTER Send 15 blaze'
                   allocate(blaze_restart_ts(wnp))
                   allocate(blaze_out_ts(wnp))
                   call master_blaze_types(comm, wland, wnp, mland, blaze, blaze_restart_ts, blaze_out_ts)
                   !if (.not. spinup) then
                   !   ! print*, 'MASTER Send 16 blaze restart'
                   !   call master_send_input(comm, blaze_restart_ts, ktau)
                   !endif
                   !deallocate(blaze_restart_ts)
                   !deallocate(blaze_out_ts)
                   if (trim(blaze%burnt_area_src) == "SIMFIRE" ) then
                      call INI_SIMFIRE(mland ,SIMFIRE, &
                         climate%modis_igbp(landpt(:)%cstart) ) !CLN here we need to check for the SIMFIRE biome setting

                      !par blaze restart not required uses climate data
                      ! print*, 'MASTER Send 17 simfire'
                      !allocate(simfire_restart_ts(wnp))
                      !allocate(simfire_inp_ts(wnp))
                      !allocate(simfire_out_ts(wnp))
                      !call master_simfire_types(comm, wland, wnp, mland, simfire, &
                      !     simfire_restart_ts, simfire_inp_ts, simfire_out_ts)
                      !if (.not. spinup) then
                      !   ! print*, 'MASTER Send 18 simfire restart'
                      !   call master_send_input(comm, simfire_restart_ts, ktau)
                      !endif
                      !deallocate(simfire_restart_ts)
                      !deallocate(simfire_inp_ts)
                      !deallocate(simfire_out_ts)
                   end if
                end if
             END IF ! icycle

             ! MPI: allocate read ahead buffers for input met and veg data
             CALL alloc_cbm_var(imet, mp)
             CALL alloc_cbm_var(iveg, mp)

             ! MPI: create inp_t types to scatter input data to the workers
             ! at the start of every timestep
             !   CALL master_intypes (comm,met,veg)
             ! for read ahead use the new variables
             ! print*, 'MASTER Send 19 intypes'
             CALL master_intypes(comm, imet, iveg)

             ! MPI: create recv_t types to receive results from the workers
             ! at the end of every timestep
             ! print*, 'MASTER Send 20 outtypes'
             CALL master_outtypes(comm, met, canopy, ssnow, rad, bal, air, soil, veg)

             ! MPI: create type for receiving casa results
             ! only if cnp module is active
             if (icycle>0) then
                ! print*, 'MASTER Send 21 casa'
                call master_casa_types(comm, casapool, casaflux, casamet, casabal, phen)
                ! 13C
                if (cable_user%c13o2) then
                   ! print*, 'MASTER Send 22 c13o2_flux'
                   call master_c13o2_flux_types(comm, c13o2flux)
                   ! print*, 'MASTER Send 23 c13o2_pool'
                   call master_c13o2_pool_types(comm, c13o2pools)
                endif

                if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                   ! print*, 'MASTER Send 24 casa_dump'
                   ! 13C
                   call master_casa_dump_types(comm, casamet, casaflux, phen, climate, c13o2flux)
                endif
                if (cable_user%popluc) then
                   ! print*, 'MASTER Send 25 casa_LUC'
                   CALL master_casa_LUC_types(comm, casapool, casabal, casaflux)
                   ! 13C
                   if (cable_user%c13o2) then
                      ! print*, 'MASTER Send 26 c13o2_luc_type'
                      call master_c13o2_luc_types(comm, c13o2luc)
                   endif
                endif
             end if

             ! MPI: create type to send restart data back to the master
             ! only if restart file is to be created
             if (output%restart) then
                ! print*, 'MASTER Send 27 restart'
                call master_restart_types(comm, canopy, air)
             end if

             ! CALL zero_sum_casa(sum_casapool, sum_casaflux)
             ! ! 13C
             ! if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
             ! count_sum_casa = 0

             ! CALL master_sumcasa_types(comm, sum_casapool, sum_casaflux)
             IF (icycle>0 .AND. spincasa) THEN
                PRINT *, 'EXT spincasacnp enabled with mloop= ', mloop, dels, kstart, kend
                ! 13C
                CALL master_spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                     casaflux,casamet,casabal,phen,POP,climate,LALLOC, &
                     c13o2flux, c13o2pools, c13o2luc, icomm, ocomm)
                ! print*, 'CD02   ', casamet%glai
                SPINconv = .FALSE.
                CASAONLY = .TRUE.
                ktau_gl  = 0
                ktau     = 0
             ELSEIF ( casaonly .AND. (.NOT. spincasa) .AND. cable_user%popluc) THEN
                ! 13C
                PRINT *, 'EXT CASAONLY_LUC'
                CALL master_CASAONLY_LUC(dels,kstart,kend,veg,soil,casabiome,casapool, &
                     casaflux,casamet,casabal,phen,POP,climate,LALLOC, LUC_EXPT, POPLUC, &
                     c13o2flux, c13o2pools, c13o2luc, icomm, ocomm)
                ! print*, 'CD03   ', casamet%glai
                SPINconv = .FALSE.
                ktau_gl  = 0
                ktau     = 0
             ENDIF

             ! MPI: mostly original serial code follows...
          ENDIF ! CALL1

          ! Open output file:
          IF (.NOT.CASAONLY) THEN
             IF ( TRIM(filename%out) .EQ. '' ) THEN
                IF ( cable_user%YEARSTART .GT. 0 ) THEN
                   WRITE(dum, FMT="(I4,'_',I4)") cable_user%YEARSTART, &
                         cable_user%YEAREND
                   filename%out = TRIM(filename%path)//'/'//&
                        TRIM(cable_user%RunIden)//'_'//&
                        TRIM(dum)//'_cable_out.nc'
                ELSE
                   filename%out = TRIM(filename%path)//'/'//&
                        TRIM(cable_user%RunIden)//'_cable_out.nc'
                ENDIF
             ENDIF
             IF (YYYY.EQ.cable_user%YEARSTART) THEN
                CALL nullify_write() ! nullify pointers
                CALL open_output_file(dels, soil, veg, bgc, rough)
             ENDIF
          ENDIF

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

          IF (spincasa .OR. casaonly) THEN
             EXIT
          ENDIF

          IF (.NOT.casaonly) THEN
             IF ( TRIM(cable_user%MetType) .EQ. 'plum' ) THEN
                CALL PLUME_MIP_GET_MET(PLUME, iMET, YYYY, 1, kend, &
                     (YYYY.EQ.cable_user%YearEnd .AND. 1.EQ.kend))
                !CALL PLUME_MIP_GET_MET(PLUME, iMET, YYYY, 1, kend, &
                !     (YYYY.EQ.cable_user%YearEnd))
             ELSE IF ( TRIM(cable_user%MetType) .EQ. 'bios' ) THEN
                CALL  cable_bios_read_met(iMET, CurYear, iktau, kend, &
                       (YYYY.EQ.cable_user%YearEnd .AND. 1.EQ.kend), dels )
             ELSE IF ( TRIM(cable_user%MetType) .EQ. 'cru' ) THEN
                CALL CRU_GET_SUBDIURNAL_MET(CRU, imet, YYYY, 1, kend, &
                     (YYYY.EQ.cable_user%YearEnd))
                !MC - Added two lines defining iveg and lai in input veg: iveg
                !     because this is done in get_met_data but not in cru_get_subdiurnal_met
                iveg%iveg = veg%iveg
                iveg%vlai = veg%vlai
             ELSE
                CALL get_met_data( spinup, spinConv, imet, soil,   &
                     rad, iveg, kend, dels, C%TFRZ, iktau+koffset, &
                     kstart+koffset )
             ENDIF
          ENDIF

          ! IF ( CASAONLY .AND. IS_CASA_TIME("dread", yyyy, iktau, kstart, koffset, &
          !      kend, ktauday, logn) ) THEN
          !    WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart+koffset)/(LOY*ktauday))
          !    ncfile  = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
          !    casa_it = NINT( REAL(iktau / ktauday) )
          !    CALL read_casa_dump( ncfile, casamet, casaflux,phen, casa_it, kend, .FALSE. )
          ! ENDIF

          canopy%oldcansto=canopy%cansto

          ! Zero out lai where there is no vegetation acc. to veg. index
          WHERE ( iveg%iveg(:) .GE. 14 ) iveg%vlai = 0.

          ! Read ahead: send input before workes start ktau loop
          IF ( .NOT. CASAONLY ) THEN
             ! MPI: scatter input data to the workers
             ! print*, 'MASTER Send 28 input'
             CALL master_send_input(icomm, inp_ts, iktau)
             ! CALL MPI_Waitall(wnp, inp_req, inp_stats, ierr)
             !MC - The else statement that sends casa_dump data is not commented out in the trunk.
             !     And in mpiworker there is:
             !         ELSEIF ( IS_CASA_TIME("dread", yyyy, ktau, kstart, koffset, kend, ktauday, wlogn) ) then
          ! ELSE
             ! CALL master_send_input(icomm, casa_dump_ts, iktau)
             ! CALL MPI_Waitall(wnp, inp_req, inp_stats, ierr)
          ENDIF

          ! 13C
          if (cable_user%c13o2) then
             if ((CurYear < c13o2_atm_syear) .or. (CurYear > c13o2_atm_eyear)) then
                write(*,*) 'Current year ', CurYear, 'not in atmospheric delta-13C (min/max): ', &
                     c13o2_atm_syear, c13o2_atm_eyear
                call MPI_Abort(comm, 0, ierr)
             endif
             c13o2flux%ca = (c13o2_delta_atm(CurYear) + 1.0_r_2) * real(imet%ca,r_2) ! * vpdbc13 / vpdbc13
             ! broadcast
             ! print*, 'MASTER Send 29 Ca '
             do rank=1, wnp
                off = wland(rank)%patch0
                cnt = wland(rank)%npatch
                call MPI_Send(c13o2flux%ca(off), cnt, MPI_DOUBLE_PRECISION, rank, 0, icomm, ierr)
             end do
          endif

          ! IF (.NOT.spincasa) THEN
          ! time step loop over ktau
          ! print*, 'II01 ', YYYY
          KTAULOOP: do ktau=kstart, kend-1
             ! print*, 'II01.01 ', ktau
             ! globally (WRT code) accessible kend through USE cable_common_module
             iktau = iktau + 1
             oktau = oktau + 1
             ktau_tot = ktau_tot + 1
             ktau_gl  = iktau

             met%year = imet%year
             met%doy  = imet%doy

             ! some things (e.g. CASA-CNP) only need to be done once per day
             idoy = int(mod(real(ktau+koffset)/real(ktauday),real(loy)))
             if (idoy .eq. 0) idoy = loy

             ! needed for CASA-CNP
             nyear = int((kend-kstart+1)/(loy*ktauday))

             ! Get met data and LAI, set time variables.
             ! Rainfall input may be augmented for spinup purposes:
             !          met%ofsd = met%fsd(:,1) + met%fsd(:,2)
             ! print*, 'II01.02 '
             if (trim(cable_user%MetType) .eq. 'plum') then
                call plume_mip_get_met(plume, imet, yyyy, iktau, kend, &
                     yyyy.eq.cable_user%YearEnd .AND. iktau.EQ.kend)
             else if (trim(cable_user%MetType) .eq. 'bios') then
                if ((.not. CASAONLY) .or. (CASAONLY.and.CALL1)) then
                   call cable_bios_read_met(imet, CurYear, iktau, kend, &
                        yyyy.eq.cable_user%yearend .and. iktau.eq.kend, dels)
                end if
             else if (trim(cable_user%MetType) .eq. 'cru') then
                ! print*, 'II01.03 ', ktau
                call cru_get_subdiurnal_met(cru, imet, YYYY, iktau, kend, &
                     yyyy.eq.cable_user%YearEnd)
                ! print*, 'II01.04 '
                !MC - Added two lines defining iveg and lai in input veg: iveg
                !     because this is done in get_met_data but not in cru_get_subdiurnal_met
                iveg%iveg = veg%iveg ! LAI and veg class not in CRU
                iveg%vlai = veg%vlai
                if (CALL1)  casamet%glai = 1.0  ! initialise glai for use in cable_roughness
             else
                call get_met_data(spinup, spinconv, imet, soil, &
                     rad, iveg, kend, dels, c%tfrz, iktau+koffset, &
                     kstart+koffset)
             endif
             if (trim(cable_user%mettype) .ne. 'gswp') CurYear = met%year(1)
             ! print*, 'II01.05 '

             !  IF ( CASAONLY .AND. IS_CASA_TIME("dread", yyyy, iktau, kstart, koffset, &
             !       kend, ktauday, logn) )  THEN
             !     ! CLN READ FROM FILE INSTEAD !
             !     WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart+koffset)/(LOY*ktauday))
             !     ncfile  = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
             !     casa_it = NINT( REAL(iktau / ktauday) )
             !     CALL read_casa_dump( ncfile, casamet, casaflux, casa_it, kend, .FALSE. )
             !  ENDIF

             !  ! At first time step of year, set tile area according to updated LU areas
             !  IF (ktau == 1 .and. cable_user%POPLUC) THEN
             !    CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
             ! ENDIF

             if (.not. CASAONLY) then

                if (icycle > 0) then
                   ! print*, 'II10'
                   ! receive casa update from worker
                   if ( mod((oktau-kstart+1+koffset),ktauday).eq.0 ) then
                      ! print*, 'MASTER Receive 30 casa'
                      call master_receive(ocomm, oktau, casa_ts)
                      ! 13C
                      if (cable_user%c13o2) then
                         ! print*, 'MASTER Receive 31 c13o2_flux'
                         call master_receive(ocomm, oktau, c13o2_flux_ts)
                         ! print*, 'MASTER Receive 32 c13o2_flux'
                         call master_receive(ocomm, oktau, c13o2_pool_ts)
                      endif
                      ! write(*,*) 'after master_receive casa_ts'
                      ! print*, 'MASTER Receive waitall'
                      !call MPI_Waitall(wnp, recv_req, recv_stats, ierr)
                   endif

                   if (cable_user%call_blaze) then
                     if ( mod((oktau-kstart+1+koffset),ktauday).eq.0 ) then
                       !par recv blaze_out_ts
                       call master_receive(ocomm, oktau, blaze_out_ts)
                       BLAZE%time =  BLAZE%time + 86400
                       call write_blaze_output_nc( BLAZE, ktau.EQ.kend .AND. YYYY.EQ.cable_user%YearEnd)
                     endif
                   endif

                   ! write(*,*) 'after master_receive casa_ts waitall'
                   ! receive casa dump requirements from worker
                   if ( ((.not.spinup) .or. (spinup.and.spinConv)) .and. &
                        is_casa_time("dwrit", yyyy, oktau, kstart, &
                        koffset, kend, ktauday, logn) ) then
                      ! print*, 'MASTER Receive 33 dump'
                      call master_receive( ocomm, oktau, casa_dump_ts )
                      ! call MPI_Waitall(wnp, recv_req, recv_stats, ierr)
                   endif
                endif ! icycle>0

                ! MPI: receive this time step's results from the workers
                ! print*, 'MASTER Receive 34 recv'
                !MC veg%iveg is changing during receive for GNU compiler on Explor
                call master_receive(ocomm, oktau, recv_ts)
                ! print*, 'MASTER Received 34'
                ! call MPI_Waitall(wnp, recv_req, recv_stats, ierr)

                ! MPI: scatter input data to the workers
                ! print*, 'MASTER Send 35 input'
                call master_send_input(icomm, inp_ts, iktau)
                ! call MPI_Waitall (wnp, inp_req, inp_stats, ierr)

                ! 13C
                if (cable_user%c13o2) then
                   ! already checked that CurYear in input file
                   c13o2flux%ca = (c13o2_delta_atm(CurYear) + 1.0_r_2) * real(imet%ca,r_2)
                   ! print*, 'MASTER Send 36 Ca'
                   do rank=1, wnp
                      off = wland(rank)%patch0
                      cnt = wland(rank)%npatch
                      call MPI_Send(c13o2flux%ca(off), cnt, MPI_DOUBLE_PRECISION, rank, 0, icomm, ierr)
                   end do
                   ! print*, 'II13'
                endif

                ! IF ( ((.NOT.spinup).OR.(spinup.AND.spinConv)) .AND.   &
                !      ( IS_CASA_TIME("dwrit", yyyy, oktau, kstart, &
                !         koffset, kend, ktauday, logn) ) ) THEN
                !    WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart)/(LOY*ktauday))
                !    ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
                !    CALL write_casa_dump( ncfile, casamet , casaflux, idoy, &
                !         kend/ktauday )
                ! ENDIF

                IF ( mod((oktau-kstart+1+koffset),ktauday)==0 ) then  ! end of day
                   ctime = ctime + 86400
                ENDIF

                if ( ((.not.spinup) .or. (spinup.and.spinConv)) .and. &
                      mod((oktau-kstart+1+koffset),ktauday)==0 ) then
                   if ( cable_user%casa_dump_write )  then
                      !cln check for leap year
                      ! print*, 'II15'
                      write(cyear,fmt="(I4)") curyear + int((oktau-kstart+koffset)/(loy*ktauday))
                      ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
                      ! 13C
                      call write_casa_dump(ncfile, casamet, casaflux, phen, climate, &
                           c13o2flux, idoy, kend/ktauday)
                      ! print*, 'II16'
                   endif
                endif

             ! else if ( mod((iktau-kstart+1+koffset),ktauday)==0 ) then
             else if (is_casa_time("dread", yyyy, iktau, kstart, koffset, &
                  kend, ktauday, logn)) then

                ! print*, 'MASTER Send 35 dump'
                call master_send_input(icomm, casa_dump_ts, iktau)
                ! call MPI_Waitall (wnp, inp_req, inp_stats, ierr)

             endif ! .not. CASAONLY

             ! call MPI_Waitall(wnp, recv_req, recv_stats, ierr)
             ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

             ! print*, 'II17'
             met%ofsd = met%fsd(:,1) + met%fsd(:,2)
             canopy%oldcansto = canopy%cansto

             ! Zero out lai where there is no vegetation acc. to veg. index
             where (iveg%iveg(:) .ge. 14) iveg%vlai = 0.

             ! Write time step's output to file if either: we're not spinning up
             ! or we're spinning up and the spinup has converged:
             ! MPI: TODO: pull mass and energy balance calculation from write_output
             ! and refactor into worker code

             ktau_gl = oktau

             ! print*, 'II18'
             if ((.not.spinup) .or. (spinup.and.spinConv)) then
                if (icycle > 0) then
                   if ( is_casa_time("write", yyyy, oktau, kstart, &
                        koffset, kend, ktauday, logn) ) then
                      !ctime = ctime + 1
                      ! print*, 'II18.1'
                      call write_casa_output_nc(veg, casamet, casapool, casabal, casaflux, &
                           CASAONLY, ctime, ktau.eq.kend .and. yyyy.eq.cable_user%YearEnd)
                      ! 13C
                      if (cable_user%c13o2) then
                         if (ctime == 86400) then
                            ! call c13o2_create_output(casamet, sum_c13o2pools, c13o2_outfile_id, c13o2_vars, c13o2_var_ids)
                            ! print*, 'II18.2'
                            call c13o2_create_output(casamet, c13o2pools, c13o2_outfile_id, c13o2_vars, c13o2_var_ids)
                         endif
                         ! print*, 'II18.3'
                         call c13o2_write_output(c13o2_outfile_id, c13o2_vars, c13o2_var_ids, ctime, c13o2pools)
                      end if
                   endif
                endif

                ! print*, 'II18.4'
                if ((.not. CASAONLY) .and. spinConv) then
                   if (trim(cable_user%mettype) .eq. 'plum' &
                       .or. trim(cable_user%mettype) .eq. 'cru' &
                       .or. trim(cable_user%mettype) .eq. 'bios' &
                       .or. trim(cable_user%mettype) .eq. 'gswp') then
                      ! print*, 'II18.5'
                      ! 13C
                      call write_output(dels, ktau_tot, met, canopy, casaflux, casapool, &
                           casamet, ssnow, &
                           rad, bal, air, soil, veg, c%sboltz, &
                           c%emleaf, c%emsoil, c13o2pools, c13o2flux)
                   else
                      ! print*, 'II18.6'
                      ! 13C
                      call write_output(dels, ktau, met, canopy, casaflux, casapool, &
                           casamet, ssnow, &
                           rad, bal, air, soil, veg, c%sboltz, c%emleaf, c%emsoil, c13o2pools, c13o2flux)
                   endif
                end if
             endif

             !---------------------------------------------------------------------!
             ! Check this run against standard for quasi-bitwise reproducability   !
             ! Check triggered by cable_user%consistency_check=.TRUE. in cable.nml !
             !---------------------------------------------------------------------!
             ! print*, 'II19'
             IF (cable_user%consistency_check) THEN
                ! print*, 'II19.01'
                count_bal = count_bal +1;
                ! print*, 'II19.02'
                new_sumbal = new_sumbal + SUM(bal%wbal)/mp +  SUM(bal%ebal)/mp
                ! print*, 'II19.03'
                new_sumfpn = new_sumfpn + SUM(canopy%fpn)/mp
                ! print*, 'II19.04'
                new_sumfe  = new_sumfe + SUM(canopy%fe)/mp
                ! if (ktau == kend-1) PRINT*, "time-space-averaged energy & water balances"
                ! if (ktau == kend-1) PRINT*,"Ebal_tot[Wm-2], Wbal_tot[mm]", &
                !      sum(bal%ebal_tot)/mp/count_bal, sum(bal%wbal_tot)/mp/count_bal
                ! if (ktau == kend-1) PRINT*, "time-space-averaged latent heat and net photosynthesis"
                ! if (ktau == kend-1) PRINT*, "sum_fe[Wm-2], sum_fpn[umol/m2/s]",  &
                !      new_sumfe/count_bal, new_sumfpn/count_bal
                
                ! check for Nans in biophysical outputs and abort if there are any
                ! print*, 'II19.05'
                IF (any(canopy%fe .NE. canopy%fe)) THEN
                   DO kk=1,mp
                      IF (canopy%fe(kk) .NE. canopy%fe(kk)) THEN
                         WRITE(*,*) 'Nan in evap flux,', kk, patch(kk)%latitude, patch(kk)%longitude
                         write(*,*) 'fe nan', kk, ktau,met%qv(kk), met%precip(kk),met%precip_sn(kk), &
                              met%fld(kk), met%fsd(kk,:), met%tk(kk), met%ua(kk), &
                              ssnow%potev(kk), met%pmb(kk), &
                              canopy%ga(kk), ssnow%tgg(kk,:), canopy%fwsoil(kk), &
                              rad%fvlai(kk,:) ,  rad%fvlai(kk,1), &
                              rad%fvlai(kk,2), canopy%vlaiw(kk)
                         CALL MPI_Abort(comm, 0, ierr)
                      ENDIF
                   ENDDO
                ENDIF

                ! print*, 'II19.06'
                IF (ktau==(kend-1)) THEN
                   ! print*, 'II19.07'
                   nkend = nkend+1
                   IF( ABS(new_sumbal-trunk_sumbal) < 1.e-7) THEN
                      write(*,*) ""
                      write(*,*) "NB. Offline-parallel runs spinup cycles:", nkend
                      write(*,*) "Internal check shows this version reproduces the trunk sumbal"
                   ELSE
                      write(*,*) ""
                      write(*,*) "NB. Offline-parallel runs spinup cycles:", nkend
                      write(*,*) "Internal check shows in this version new_sumbal != trunk sumbal"
                      write(*,*) "The difference is: ", new_sumbal - trunk_sumbal
                      write(*,*) "Writing new_sumbal to the file:", TRIM(Fnew_sumbal)
                      !CLN OPEN( 12, FILE = Fnew_sumbal )
                      !CLN WRITE( 12, '(F20.7)' ) new_sumbal  ! written by previous trunk version
                      !CLN CLOSE(12)
                   ENDIF
                ENDIF
             ENDIF ! cable_user%consistency_check

             CALL1 = .false.

             !WRITE(*,*) " ktauloop end ", ktau, CurYear
             !if (met%doy(1).eq.10) then
             !CALL CPU_TIME(etime)
             !PRINT *, 'Master ktau ',ktau,  etime, etime-etimelast, ' seconds'
             !etimelast = etime
             ! endif

          end do KTAULOOP ! END Do loop over timestep ktau
          ! print*, 'II20'

          CALL1 = .false.

          ! MPI: read ahead tail to receive (last step and write)
          met%year = imet%year
          met%doy  = imet%doy
          oktau    = oktau + 1
          ktau_tot = ktau_tot + 1
          ktau_gl  = oktau

          IF ( .NOT. CASAONLY ) THEN
             IF ( icycle >0 ) THEN
                ! print*, 'MASTER Receive 37 casa'
                CALL master_receive(ocomm, oktau, casa_ts)
                ! 13C
                if (cable_user%c13o2) then
                   ! print*, 'MASTER Receive 38 c13o2_flux'
                   call master_receive(ocomm, oktau, c13o2_flux_ts)
                   ! print*, 'MASTER Receive 39 c13o2_pool'
                   call master_receive(ocomm, oktau, c13o2_pool_ts)
                endif
                if (cable_user%call_blaze) then
                  !par recv blaze_out_ts
                  call master_receive(ocomm, oktau, blaze_out_ts)
                  BLAZE%time =  BLAZE%time + 86400
                  call write_blaze_output_nc( BLAZE, ktau.EQ.kend .AND. YYYY.EQ.cable_user%YearEnd)
                end if
                ! CALL MPI_Waitall(wnp, recv_req, recv_stats, ierr)
                IF ( ((.NOT.spinup) .OR. (spinup.AND.spinConv)) .AND. &
                     IS_CASA_TIME("dwrit", yyyy, oktau, kstart, &
                     koffset, kend, ktauday, logn) ) THEN
                   ! print*, 'MASTER Receive 40 dump'
                   CALL master_receive( ocomm, oktau, casa_dump_ts )
                   ! CALL MPI_Waitall(wnp, recv_req, recv_stats, ierr)
                ENDIF
             ENDIF
             ! write(*,*) 'master_receive 6'
             ! print*, 'MASTER Receive 41 recv'
             CALL master_receive(ocomm, oktau, recv_ts)
             ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
          ENDIF

          met%ofsd = met%fsd(:,1) + met%fsd(:,2)
          canopy%oldcansto=canopy%cansto

          IF ( TRIM(cable_user%MetType) .EQ. "gswp" ) &
               CALL close_met_file()
          
          IF (icycle>0 .and. cable_user%CALL_POP)  THEN
             ! write(*,*) 'b4 annual calcs'
             IF (cable_user%POPLUC) THEN
                ! master receives casa updates required for LUC calculations here
                ! print*, 'MASTER Receive 42 luc'
                CALL master_receive(ocomm, 0, casa_LUC_ts)
                ! 13C
                if (cable_user%c13o2) then
                   ! print*, 'MASTER Receive 43 c13o2_luc'
                   call master_receive(ocomm, 0, c13o2_flux_ts)
                   call master_receive(ocomm, 0, c13o2_pool_ts)
                   call master_receive(ocomm, 0, c13o2_luc_ts)
                endif
                ! Dynamic LUC
                ! 13C
                CALL LUCdriver( casabiome,casapool,casaflux,POP,LUC_EXPT, POPLUC, veg, c13o2pools )
                ! transfer POP updates to workers
                ! print*, 'MASTER Send 44 pop'
                off = 1
                DO rank = 1, wnp
                   IF ( rank .GT. 1 ) off = off + wland(rank-1)%npop_iwood
                   cnt = wland(rank)%npop_iwood
                   CALL MPI_Send( POP%pop_grid(off), cnt, pop_ts, rank, 0, icomm, ierr )
                END DO
             ENDIF

             ! one annual time-step of POP (worker calls POP here)
             !CALL POPdriver(casaflux,casabal,veg, POP)
             ! print*, 'MASTER Receive 45 pop'
             CALL master_receive_pop(POP, ocomm)

             IF (cable_user%POPLUC) THEN
                ! Dynamic LUC: update casa pools according to LUC transitions
                ! 13C
                if (cable_user%c13o2) call c13o2_save_luc(casapool, popluc, casasave, lucsave)
                CALL POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
                ! 13C
                if (cable_user%c13o2) then
                   call c13o2_update_luc(casasave, lucsave, popluc, &
                        luc_expt%prim_only, c13o2pools, c13o2luc)
                endif
                
                ! Dynamic LUC: write output
                IF (output%grid(1:3) == 'lan') THEN
                   CALL WRITE_LUC_OUTPUT_NC( POPLUC, YYYY, ( YYYY.EQ.cable_user%YearEnd ))
                ELSE
                   CALL WRITE_LUC_OUTPUT_GRID_NC( POPLUC, YYYY, ( YYYY.EQ.cable_user%YearEnd ))
                   !CALL WRITE_LUC_OUTPUT_NC( POPLUC, YYYY, ( YYYY.EQ.cable_user%YearEnd ))
                ENDIF
             ENDIF

             if (cable_user%POPLUC) then
                ! send updates for CASA pools, resulting from LUC
                ! print*, 'MASTER Send 46 luc'
                CALL master_send_input(icomm, casa_LUC_ts, nyear)
                ! 13C
                if (cable_user%c13o2) then
                   ! print*, 'MASTER Send 47 c13o2_luc'
                   call master_send_input(icomm, c13o2_flux_ts, nyear)
                   call master_send_input(icomm, c13o2_pool_ts, nyear)
                   call master_send_input(icomm, c13o2_luc_ts, nyear)
                endif
             endif
          ENDIF ! icycle>0 .and. cable_user%CALL_POP

          ! WRITE OUTPUT
          IF ((.NOT.spinup).OR.(spinup.AND.spinConv)) THEN
             IF (icycle>0) THEN
                ctime = ctime + 86400
                !TRUNK no if but call write_casa in any case
                if ( is_casa_time("write", yyyy, ktau, kstart, koffset, kend, ktauday, logn) ) then
                   call write_casa_output_nc(veg, casamet, casapool, casabal, casaflux, &
                        CASAONLY, ctime, (ktau==kend) .and. (YYYY==cable_user%YearEnd))
                   ! 13C
                   if (cable_user%c13o2) then
                      call c13o2_write_output(c13o2_outfile_id, c13o2_vars, c13o2_var_ids, ctime, c13o2pools)
                      if (YYYY == cable_user%YearEnd) then
                         call c13o2_close_output(c13o2_outfile_id)
                      endif
                   end if
                endif
                IF ( cable_user%CALL_POP ) THEN
                   ! CALL master_receive_pop(POP, ocomm)
                   ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
                   IF ( TRIM(cable_user%POP_out).eq.'epi' ) THEN
                      CALL POP_IO( pop, casamet, CurYear, 'WRITE_EPI', &
                           (CurYear.EQ.cable_user%YearEnd) )
                   ENDIF
                ENDIF
             END IF

             ! IF ( cable_user%CASA_DUMP_WRITE )  THEN
             !    !CLN CHECK FOR LEAP YEAR
             !    WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart)/(LOY*ktauday))
             !    ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
             !    CALL write_casa_dump( ncfile, casamet , casaflux, idoy, &
             !         kend/ktauday )
             ! ENDIF

             IF (((.NOT.spinup).OR.(spinup.AND.spinConv)).and. &
                  MOD((ktau-kstart+1),ktauday)==0) THEN
                !ctime = ctime + 86400  ! update casa time
                IF ( cable_user%CASA_DUMP_WRITE )  THEN
                   !CLN CHECK FOR LEAP YEAR
                   WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart)/(LOY*ktauday))
                   ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
                   ! 13C
                   CALL write_casa_dump( ncfile, casamet , casaflux, phen, climate, c13o2flux, &
                   LOY,  kend/ktauday )
                ENDIF
             ENDIF

             IF ( (.NOT. CASAONLY) .AND. spinConv ) THEN
                IF ( TRIM(cable_user%MetType) .EQ. 'plum' &
                     .OR. TRIM(cable_user%MetType) .EQ. 'cru'   &
                     .OR. TRIM(cable_user%MetType) .EQ. 'bios'   &
                     .OR. TRIM(cable_user%MetType) .EQ. 'gswp') then
                   ! 13C
                   CALL write_output( dels, ktau_tot, met, canopy, casaflux, casapool, &
                        casamet, ssnow,         &
                        rad, bal, air, soil, veg, C%SBOLTZ,     &
                        C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux )
                ELSE
                   ! 13C
                   CALL write_output( dels, ktau, met, canopy, casaflux, casapool, casamet, &
                        ssnow,   &
                        rad, bal, air, soil, veg, C%SBOLTZ, C%EMLEAF, C%EMSOIL, c13o2pools, c13o2flux )
                ENDIF
             END IF

             IF (cable_user%consistency_check) THEN
                count_bal = count_bal +1;
                new_sumbal = new_sumbal + SUM(bal%wbal)/mp +  SUM(bal%ebal)/mp
                new_sumfpn = new_sumfpn + SUM(canopy%fpn)/mp
                new_sumfe = new_sumfe + SUM(canopy%fe)/mp
                if (ktau == kend) write(*,*) ""
                if (ktau == kend) write(*,*) "time-space-averaged energy & water balances"
                if (ktau == kend) write(*,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep]", &
                     sum(bal%ebal_tot)/mp/count_bal, sum(bal%wbal_tot)/mp/count_bal
                if (ktau == kend) write(*,*) "time-space-averaged latent heat and net photosynthesis"
                if (ktau == kend) write(*,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s]",  &
                     new_sumfe/count_bal, new_sumfpn/count_bal
                if (ktau == kend) write(logn,*) ""
                if (ktau == kend) write(logn,*) "time-space-averaged energy & water balances"
                if (ktau == kend) write(logn,*) "Ebal_tot[Wm-2], Wbal_tot[mm per timestep] ", &
                     sum(bal%ebal_tot)/mp/count_bal, sum(bal%wbal_tot)/mp/count_bal
                if (ktau == kend) write(logn,*) "time-space-averaged latent heat and net photosynthesis"
                if (ktau == kend) write(logn,*) "sum_fe[Wm-2], sum_fpn[umol/m2/s] ",  &
                     new_sumfe/count_bal, new_sumfpn/count_bal
             ENDIF

          END IF ! (.NOT.spinup).OR.(spinup.AND.spinConv)

          ! set tile area according to updated LU areas
          IF (cable_user%POPLUC) THEN
             CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
          ENDIF

       END DO YEARLOOP

       IF (spincasa.or.casaonly) THEN
          EXIT
       ENDIF

       ! jhan this is insufficient testing. condition for
       ! spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
       ! see if spinup (if conducting one) has converged:
       IF (spinup.AND..NOT.spinConv.AND. .NOT. CASAONLY) THEN
          ! Write to screen and log file:
          WRITE(*,'(A18,I3,A24)') ' Spinning up: run ', INT(ktau_tot/kend), &
               ' of data set complete...'
          WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ', INT(ktau_tot/kend), &
               ' of data set complete...'

          ! IF not 1st run through whole dataset:
          IF ( INT( ktau_tot/kend ) > 1 ) THEN

             ! evaluate spinup
             IF ( ANY( ABS(ssnow%wb-soilMtemp)>delsoilM).OR.                     &
                  ANY(ABS(ssnow%tgg-soilTtemp)>delsoilT) ) THEN
                ! No complete convergence yet
                !               PRINT *, 'ssnow%wb : ', ssnow%wb
                !               PRINT *, 'soilMtemp: ', soilMtemp
                !               PRINT *, 'ssnow%tgg: ', ssnow%tgg
                !               PRINT *, 'soilTtemp: ', soilTtemp
                maxdiff = MAXLOC(ABS(ssnow%wb-soilMtemp))
                write(*,*) 'Example location of moisture non-convergence: ',maxdiff
                write(*,*) 'ssnow%wb : ', ssnow%wb(maxdiff(1),maxdiff(2))
                write(*,*) 'soilMtemp: ', soilMtemp(maxdiff(1),maxdiff(2))
                maxdiff = MAXLOC(ABS(ssnow%tgg-soilTtemp))
                write(*,*) 'Example location of temperature non-convergence: ',maxdiff
                write(*,*) 'ssnow%tgg: ', ssnow%tgg(maxdiff(1),maxdiff(2))
                write(*,*) 'soilTtemp: ', soilTtemp(maxdiff(1),maxdiff(2))
             ELSE ! spinup has converged
                spinConv = .TRUE.
                ! Write to screen and log file:
                WRITE(*,'(A33)') ' Spinup has converged - final run'
                WRITE(logn,'(A52)')                                             &
                     ' Spinup has converged - final run - writing all data'
                WRITE(logn,'(A37,F8.5,A28)')                                    &
                     ' Criteria: Change in soil moisture < ',             &
                     delsoilM, ' in any layer over whole run'
                WRITE(logn,'(A40,F8.5,A28)' )                                   &
                     '           Change in soil temperature < ',          &
                     delsoilT, ' in any layer over whole run'
             END IF

          ELSE ! allocate variables for storage

             IF (.NOT.ALLOCATED(soilMtemp)) ALLOCATE(  soilMtemp(mp,ms) )
             IF (.NOT.ALLOCATED(soilTtemp)) ALLOCATE(  soilTtemp(mp,ms) )

          END IF ! INT( ktau_tot/kend ) > 1

          IF ( YYYY.EQ. cable_user%YearEnd ) THEN
             ! store soil moisture and temperature
             soilTtemp = ssnow%tgg
             soilMtemp = REAL(ssnow%wb)
          END IF
          ! MPI:
          loop_exit = .FALSE.

       ELSE

          ! if not spinning up, or spin up has converged, exit:
          ! EXIT
          ! MPI:
          loop_exit = .TRUE.

       END IF ! spinup.AND..NOT.spinConv.AND. .NOT. CASAONLY

       ! MPI: let the workers know whether it's time to quit
       ! print*, 'MASTER Send 48 loop_exit'
       CALL MPI_Bcast(loop_exit, 1, MPI_LOGICAL, 0, comm, ierr)

       IF (loop_exit) THEN
          EXIT
       END IF

    END DO SPINLOOP

    IF (icycle > 0 .and. (.not.spincasa).and. (.not.casaonly)) THEN
       ! MPI: gather casa results from all the workers
       ! print*, 'MASTER Receive 49 casa'
       CALL master_receive(ocomm, ktau_gl, casa_ts)
       ! 13C
       if (cable_user%c13o2) then
          ! print*, 'MASTER Receive 50 c13o2_flux'
          call master_receive(ocomm, ktau_gl, c13o2_flux_ts)
          ! print*, 'MASTER Receive 51 c13o2_pool'
          call master_receive(ocomm, ktau_gl, c13o2_pool_ts)
       endif
       if (cable_user%call_blaze) then
         !par recv blaze_out_ts
         call master_receive(ocomm, ktau_gl, blaze_out_ts)
       end if

       ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
       ! CALL casa_poolout( ktau, veg, soil, casabiome, &
       !     casapool, casaflux, casamet, casabal, phen )
       CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
       CALL write_casa_restart_nc( casamet, casapool,casaflux,phen,CASAONLY )
       ! CALL write_casa_restart_nc ( casamet, casapool, met, CASAONLY )
       ! 13C
       if (cable_user%c13o2) then
          call c13o2_write_restart_flux(c13o2flux)
          call c13o2_write_restart_pools(c13o2pools)
          if (cable_user%POPLUC) call c13o2_write_restart_luc(c13o2luc)
          ! While testing
          call c13o2_print_delta_flux(c13o2flux)
          call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
          if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
       endif

       IF ( cable_user%CALL_POP .and.POP%np.gt.0 ) THEN
          IF ( CASAONLY .OR. cable_user%POP_fromZero &
               .or.TRIM(cable_user%POP_out).eq.'ini' ) THEN
             CALL POP_IO( pop, casamet, CurYear+1, 'WRITE_INI', .TRUE.)
          ELSE
             CALL POP_IO( pop, casamet, CurYear+1, 'WRITE_RST', .TRUE.)
          ENDIF
       END IF
       IF (cable_user%POPLUC .and. .NOT. CASAONLY ) THEN
          ! print*, 'YYYY ', CurYear, YYYY, cable_user%YearEnd
          CALL WRITE_LUC_RESTART_NC ( POPLUC, YYYY )
       ENDIF
    else
       ! 13C
       ! While testing
       if (cable_user%c13o2) then
          call c13o2_print_delta_flux(c13o2flux)
          call c13o2_print_delta_pools(casapool, casaflux, c13o2pools)
          if (cable_user%POPLUC) call c13o2_print_delta_luc(popluc, c13o2luc)
       endif
    END IF

    ! Write restart file if requested:
    IF (output%restart .AND. (.NOT. CASAONLY)) THEN
       ! MPI: TODO: receive variables that are required by create_restart
       ! but not write_output
       ! CALL receive_restart(comm,ktau,dels,soil,veg,ssnow, &
       !                      canopy,rough,rad,bgc,bal)
       ! gol124: how about call master_receive (comm, ktau, restart_ts)
       ! instead of a separate receive_restart sub?
       ! print*, 'MASTER Receive 51 restart'
       CALL master_receive(comm, ktau_gl, restart_ts)
       ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)

       CALL create_restart(logn, dels, ktau, soil, veg, ssnow, &
            canopy, rough, rad, bgc, bal, met)

       if (cable_user%CALL_climate) then
          ! print*, 'MASTER Receive 52 climate'
          CALL master_receive(comm, ktau_gl, climate_ts)
          ! CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
          CALL WRITE_CLIMATE_RESTART_NC( climate, ktauday )
       END IF
    END IF

    ! Close met data input file:
    IF ( TRIM(cable_user%MetType) .NE. "gswp" .AND. &
         TRIM(cable_user%MetType) .NE. "bios" .AND. &
         TRIM(cable_user%MetType) .NE. "plum" .AND. &
         TRIM(cable_user%MetType) .NE. "cru") CALL close_met_file
    IF (.NOT. CASAONLY) THEN
       ! Close output file and deallocate main variables:
       CALL close_output_file(bal)
      ! WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck
    ENDIF

    ! if (trim(cable_user%MetType) == 'cru') call cru_close(CRU)
  
    if (cable_user%POPLUC) call close_luh2(LUC_EXPT)

    ! Close log file
    CLOSE(logn)

    CALL CPU_TIME(etime)
    write(*,*) 'Master End in ', etime, ' seconds.'
    ! MPI: cleanup
    CALL master_end(icycle, output%restart)

    RETURN

  END SUBROUTINE mpidrv_master


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
  !CLN  ENDIF
  !CLN
  !CLNEND SUBROUTINE renameFiles


! ============== PRIVATE SUBROUTINES USED ONLY BY THE MPI MASTER ===============


! MPI: calculates and sends grid decomposition info to the workers
SUBROUTINE master_decomp(comm, mland, mp)

  use mpi

  USE cable_IO_vars_module, ONLY : landpt

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: comm  ! MPI communicator to talk to the workers
  INTEGER, INTENT(IN) :: mland ! total number of landpoints in the global grid
  INTEGER, INTENT(IN) :: mp    ! total number of land patches in the global grid

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
     CALL MPI_Abort(comm, 0, ierr)
  END IF

  ! MPI: calculate landpoint distribution among the workers
  ! this version distributes landpoints rather than active patches,
  ! but perhaps this will be easy to change in the future?
  lpw = mland / wnp
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
        CALL MPI_Abort(comm, 0, ierr)
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
  ! switch soil colour albedo calc - Ticket #27
  USE cable_common_module, ONLY: calcsoilalbedo

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

  ! vars intro for Ticket #27
  IF (calcsoilalbedo) THEN
     ntyp = nparam + 1
  END IF

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

     bidx = bidx + 1
     CALL MPI_Get_address (met%ca(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%year(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (met%moy(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (met%doy(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%hod(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%fsd(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%fld(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%precip(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%precip_sn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%tk(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%tvair(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%tvrad(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%pmb(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%ua(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%qv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%qvair(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%da(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%dva(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%coszen(off), displs(bidx), ierr)
     blen(bidx) = r1len


     ! ----------- air --------------

     bidx = bidx + 1
     CALL MPI_Get_address (air%rho(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%volm(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%rlam(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%qsat(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%epsi(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%visc(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%psyc(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%dsatdk(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%cmolar(off), displs(bidx), ierr)
     blen(bidx) = r1len


     ! ----------- ssnow --------------

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%dtmlt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (3, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%pudsto(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%pudsmx(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%albsoilsn(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%cls(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%dfn_dtg(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%dfh_dtg(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%dfe_ddq(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%ddq_dtg(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%evapsn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%fwtop(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%fwtop1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%fwtop2(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%fwtop3(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%gammzz(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%isflag(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%osnowd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%potev(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%pwb_min(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%runoff(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%rnof1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%rnof2(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%rtsoil(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%sconds(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msn * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%sdepth(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msn * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%smass(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msn * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%snage(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%snowd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%smelt(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%ssdn(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msn * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%ssdnn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tgg(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tggsn(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msn * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tss(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wb(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wbfice(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wbice(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wblf(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     ! additional  for sli
     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%S(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1


     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%Tsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
!!$
     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%thetai(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1


     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%snowliq(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (3, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%Tsurface(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%h0(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%nsnow(off), displs(bidx), ierr)
     blen(bidx) = I1len
     ! end additional for sli


     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wbtot(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wb_lake(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%sinfil(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%evapfbl(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%qstss(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wetfac(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%owetfac(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%t_snwlr(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tggav(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%otss(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%otss_0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! ----------- veg --------------

     bidx = bidx + 1
     CALL MPI_Get_address (veg%canst1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%dleaf(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%ejmax(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%frac4(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%froot(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%hc(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%iveg(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (veg%ivegp(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (veg%ilu(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (veg%meth(off), displs(bidx), ierr)
! Maciej: veg%meth is REAL
!     blen(bidx) = I1LEN
     blen(bidx) = R1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (veg%rp20(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%rpcoef(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%shelrb(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%wai(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%vegcf(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%tminvj(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%tmaxvj(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%vbeta(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%xalbnir(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%vcmax(off), displs(bidx), ierr)
     blen(bidx) = r1len

    ! bidx = bidx + 1
    ! CALL MPI_Get_address (veg%vlai(off), displs(bidx), ierr)
    ! blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%xfang(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%extkn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%deciduous(off), displs(bidx), ierr)
! Maciej: deciduous is logical
!     blen(bidx) = r1len
     blen(bidx) = llen

     bidx = bidx + 1
     CALL MPI_Get_address (veg%a1gs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%d0gs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%alpha(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%convex(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%cfrd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%gswmin(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%conkc0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%conko0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%ekc(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%eko(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%clitt(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%zr(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%gamma(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%refl(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (2, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (veg%taul(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (2, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1


     bidx = bidx + 1
     CALL MPI_Get_address (veg%disturbance_interval(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (2, i1len, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (veg%disturbance_intensity(off,1), displs(bidx), ierr)
! Maciej: disturbance_intensity is REAL(r_2)
!     CALL MPI_Type_create_hvector (2, r1len, r1stride, MPI_BYTE, &
!          &                             types(bidx), ierr)
     CALL MPI_Type_create_hvector (2, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     ! Ticket #56, adding veg parms for Medlyn model
     bidx = bidx + 1
     CALL MPI_Get_address (veg%g0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%g1(off), displs(bidx), ierr)
     blen(bidx) = r1len
     ! Ticket #56, finish adding new veg parms

     bidx = bidx + 1
     CALL MPI_Get_address (veg%gmmax(off), displs(bidx), ierr)
     blen(bidx) = r1len


  ! ----------- bgc --------------

     bidx = bidx + 1
     CALL MPI_Get_address (bgc%cplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ncp * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bgc%csoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ncs, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ncs * r1len

     ! constant * ncp, each worker gets the same copy of whole array
     bidx = bidx + 1
     CALL MPI_Get_address (bgc%ratecp, displs(bidx), ierr)
     blen(bidx) = ncp * extr1

     ! constant * ncs, each worker gets the same copy of whole array
     bidx = bidx + 1
     CALL MPI_Get_address (bgc%ratecs, displs(bidx), ierr)
     blen(bidx) = ncs * extr1

     ! ----------- soil --------------

     bidx = bidx + 1
     CALL MPI_Get_address (soil%albsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%bch(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%c3(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%clay(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%cnsd(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%css(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%hsbh(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%hyds(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%i2bp3(off), displs(bidx), ierr)
     ! Maciej: i2bp3 is REAL
     !     blen(bidx) = I1LEN
     blen(bidx) = R1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (soil%ibp2(off), displs(bidx), ierr)
     ! Maciej: ibp2 is REAL
     !     blen(bidx) = I1LEN
     blen(bidx) = R1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (soil%isoilm(off), displs(bidx), ierr)
     ! Maciej isoilm is INTEGER
     !     blen(bidx) = r1len
     blen(bidx) = i1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%rhosoil(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%rs20(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%sand(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%sfc(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%silt(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%ssat(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%sucs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%swilt(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! the next two are extra for sli

     bidx = bidx + 1
     CALL MPI_Get_address (soil%zeta(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%fsatmax(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (soil%ishorizon(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, i1len, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     ! end extra sli

     ! constant * ms, each worker gets the same copy of whole array
     bidx = bidx + 1
     CALL MPI_Get_address (soil%zse, displs(bidx), ierr)
     blen(bidx) = ms * extr1

     ! constant * (ms+1), each worker gets the same copy of whole array
     bidx = bidx + 1
     CALL MPI_Get_address (soil%zshh, displs(bidx), ierr)
     blen(bidx) = (ms + 1) * extr1

     ! vars intro for Ticket #27
     IF (calcsoilalbedo) THEN
        bidx = bidx + 1
        CALL MPI_Get_address (soil%soilcol(off), displs(bidx), ierr)
        blen(bidx) = r1len
     END IF

     ! ----------- canopy --------------

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fess(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fesp(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%cansto(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%oldcansto(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%cduv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%delwc(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%dewmm(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%dgdtg(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fe(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fh(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fpn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frpw(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frpr(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fnee(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frday(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fnv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fev(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fevc(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fevw(off), displs(bidx), ierr)
     blen(bidx) = r1len

     !  bidx = bidx + 1
     !  CALL MPI_Get_address (canopy%potev_c(off), displs(bidx), ierr)
     !  blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fhv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fhvw(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fns(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fes(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fhs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fwet(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%ga(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%ghflux(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%precis(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%qscrn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%rnet(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%segg(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%sghflux(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%spill(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%through(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%tscrn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%tv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%us(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%uscrn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%vlaiw(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%rghlai(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%wcint(off), displs(bidx), ierr)
     blen(bidx) = r1len

     !  bidx = bidx + 1
     !  CALL MPI_Get_address (canopy%rwater(off,1), displs(bidx), ierr)
     !  CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
     !  &                             types(bidx), ierr)
     !  blen(bidx) = 1
     !  !blen(bidx) = ms * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%evapfbl(off,1), displs(bidx), ierr)
     ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%epot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fnpp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fevw_pot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%gswx_T(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%cdtq(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%wetfac_cs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fwsoil(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%gswx(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%zetar(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (niter, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     ! 13C
     bidx = bidx + 1
     CALL MPI_Get_address(canopy%An(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%Rd(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%isc3(off), displs(bidx), ierr)
     blen(bidx) = llen

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%vcmax(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%gammastar(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%gsc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%gbc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%gac(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%ci(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1


     ! ------- rough -------

     bidx = bidx + 1
     CALL MPI_Get_address (rough%coexp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%disp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%hruff(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%hruff_grmx(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%rt0us(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%rt1usa(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%rt1usb(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%rt1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%term2(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%term3(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%term5(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%term6(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%usuh(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%za_uv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%za_tq(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%z0m(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%zref_uv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%zref_tq(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%zruffs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%z0soilsn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rough%z0soil(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! --------rad --------

     bidx = bidx + 1
     CALL MPI_Get_address (rad%albedo(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%extkb(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%extkd2(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%extkd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%flws(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%fvlai(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mf * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%gradis(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mf * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%latitude(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%lwabv(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%qcan(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mf*nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     ! blen(bidx) = mf * nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%qssabs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%rhocdf(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%rniso(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mf * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%scalex(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mf * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%transd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%trad(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%reffdf(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%reffbm(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%extkbm(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%extkdm(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%fbeam(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%cexpkbm(off,1), displs(bidx), ierr)
     ! Maciej: cexpkbm is mp*swb
     !     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
     !          &                             types(bidx), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%cexpkdm(off,1), displs(bidx), ierr)
     ! Maciej: cexpkdm is mp*swb
     !     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
     !          &                             types(bidx), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = nrb * r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%rhocbm(off,1), displs(bidx), ierr)
     ! Maciej: rhocbm is mp*nrb
     !     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
     !          &                             types(bidx), ierr)
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (rad%transb(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! ------- sum_flux -----

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumpn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumrp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumrpw(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumrpr(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumrs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumrd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%dsumpn(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%dsumrp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%dsumrs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%dsumrd(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumxrp(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (sum_flux%sumxrs(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! ------- bal ----

     bidx = bidx + 1
     CALL MPI_Get_address (bal%drybal(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%ebal(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%ebal_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%ebal_cncheck(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%ebal_tot_cncheck(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%evap_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%osnowd0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%precip_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%rnoff_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%wbal(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%wbal_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%wbtot0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%wetbal(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%owbtot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%evapc_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%evaps_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%rnof1_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%rnof2_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%snowdc_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%wbal_tot1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%delwc_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%qasrf_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%qfsrf_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%qssrf_tot(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! additional field missing from previous versions;
     ! added when trying to fix a bug in the new mpi code
     ! the order of these new fields follows the order of their
     ! declaration in cable_define_types.F90

     bidx = bidx + 1
     CALL MPI_Get_address (bal%ebaltr(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%ebal_tottr(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (bal%cansto0(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%iantrct(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tss_p(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%deltss(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%owb1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%wbtot1(off), displs(bidx), ierr)
     blen(bidx) = r1len

     ! Maciej: duplicate!
     !     bidx = bidx + 1
     !     CALL MPI_Get_address (ssnow%wbtot1(off), displs(bidx), ierr)
     !     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tprecip(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tevap(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%trnoff(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%totenbal(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%totenbal2(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%fland(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%ifland(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%tilefrac(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (n_tiles, r1len, r1stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%qasrf(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%qfsrf(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (ssnow%qssrf(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (veg%vlaimax(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%albedo_T(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (rad%longitude(off), displs(bidx), ierr)
     blen(bidx) = r1len


     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of param_t fields ', bidx, ', fix it (01)!'
        call MPI_Abort(comm, 1, ierr)
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
     CALL MPI_Abort (comm, 0, ierr)
  END IF

  ! so, now send all the parameters
  CALL master_send_input(comm, param_t, 0)
  !  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank = 1, wnp
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

  USE casavariable
  USE phenvariable

  IMPLICIT NONE

  ! sub arguments
  INTEGER, INTENT(IN) :: comm  ! MPI communicator
  ! TODO: have these variables been already allocated?
  TYPE(casa_biome),    INTENT(INOUT) :: casabiome
  TYPE(casa_pool),     INTENT(INOUT) :: casapool
  TYPE(casa_flux),     INTENT(INOUT) :: casaflux
  TYPE(casa_met),      INTENT(INOUT) :: casamet
  TYPE(casa_balance),  INTENT(INOUT) :: casabal
  TYPE(phen_variable), INTENT(INOUT) :: phen

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: ierr
  ! INTEGER :: landp_t, patch_t, param_t
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

  ALLOCATE (casa_t(wnp))

  ALLOCATE (blen(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

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

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ivt2, displs(bidx), ierr)
     blen(bidx) = mvtype * extid

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkleafcoldmax, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkleafcoldexp, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkleafdrymax, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkleafdryexp, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%glaimax, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%glaimin, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%sla, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratiofrootleaf, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%kroot, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%krootlen, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%rootdepth, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%kuptake, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%kminN, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%KuplabP, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%kclabrate, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xnpmax, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%q10soil, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkoptlitter, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkoptsoil, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%maxfinelitter, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%maxcwd, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%prodptase, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%costnpup, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkplab, displs(bidx), ierr)
     blen(bidx) = mso * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkpsorb, displs(bidx), ierr)
     blen(bidx) = mso * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%xkpocc, displs(bidx), ierr)
     blen(bidx) = mso * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%nintercept, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%nslope, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%plantrate, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%rmplant, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%fracnpptoP, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%fraclignin, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%fraclabile, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratioNCplantmin, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratioNCplantmax, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratioPCplantmin, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratioPCplantmax, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%fracLigninplant, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ftransNPtoL, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ftransPPtoL, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%litterrate, displs(bidx), ierr)
     blen(bidx) = mvtype * mlitter * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%soilrate, displs(bidx), ierr)
     blen(bidx) = mvtype * msoil * extr2

     ! added by ln
     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratioNPplantmin, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%ratioNPplantmax, displs(bidx), ierr)
     blen(bidx) = mvtype * mplant * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%la_to_sa, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%vcmax_scalar, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%disturbance_interval, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

       bidx = bidx + 1
     CALL MPI_Get_address (casabiome%DAMM_EnzPool, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%DAMM_KMO2, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%DAMM_KMcp, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%DAMM_Ea, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (casabiome%DAMM_alpha, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     ! ------ casapool ----

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Clabile(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dClabiledt(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Cplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Nplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Pplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dCplantdt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dNplantdt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dPplantdt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioPCplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Nsoilmin(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Psoillab(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Psoilsorb(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Psoilocc(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dNsoilmindt(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dPsoillabdt(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dPsoilsorbdt(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dPsoiloccdt(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Clitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Nlitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Plitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dClitterdt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dNlitterdt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dPlitterdt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNClitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioPClitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Csoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Nsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Psoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dCsoildt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dNsoildt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%dPsoildt(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioPCsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCsoilnew(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCsoilmin(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCsoilmax(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     ! added by LN
     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNPplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNPlitter(off,1), displs(bidx), ierr)
     !  CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     ! Maciej
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNPsoil(off,1), displs(bidx), ierr)
     !  CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     ! Maciej
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     ! ------- casaflux ----

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cgpp(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cnpp(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crp(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crgplant(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminfix(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminuptake(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Plabuptake(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Clabloss(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracClabile(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracCalloc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracNalloc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracPalloc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crmplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kplant_fire(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     !blen(bidx) = mplant * r2len

     ! gol124: temp
     ! 3D
     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fromPtoL(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant * mlitter, r2len, r2stride, MPI_BYTE, &
                                  types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * mlitter * r2len

     ! bidx = bidx + 1
     ! CALL MPI_Get_address(casaflux%fromPtoL_fire(off,1,1), displs(bidx), ierr)
     ! CALL MPI_Type_create_hvector(mplant * mlitter, r2len, r2stride, MPI_BYTE, &
     !                              types(bidx), ierr)
     ! blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cnep(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crsoil(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nmindep(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminloss(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminleach(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nupland(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nlittermin(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nsmin(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nsimm(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nsnet(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fNminloss(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fNminleach(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Pdep(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Pwea(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Pleach(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Ploss(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Pupland(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Plittermin(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Psmin(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Psimm(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Psnet(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fPleach(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kplab(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kpsorb(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kpocc(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kmlabP(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Psorbmax(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%frac_sapwood(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%sapwood_area(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fHarvest(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fCrop(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%klitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%klitter_fire(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%ksoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     ! gol124: temp only
     ! 3D
     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromLtoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil * mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * mlitter * r2len

     ! gol124: temp only
     ! 3D
     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromStoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil * msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromLtoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromStoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxCtolitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxNtolitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxPtolitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxCtosoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxNtosoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxPtosoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxCtoCO2(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! 13C
     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromPtoL(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant*mlitter, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromLtoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter*msoil, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromStoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(msoil*msoil, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromPtoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromLtoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromStoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(msoil, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromPtoHarvest(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! ------- casamet ----

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%glai(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%Tairk(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%precip(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%tsoilavg(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%moistavg(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%btran(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%lnonwood(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%Tsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%moist(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = ms * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%iveg2(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%ijgcm(off), displs(bidx), ierr)
     ! Maciej: ijgcm is INTEGER
     blen(bidx) = i1len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%isorder(off), displs(bidx), ierr)
     ! Maciej: isorder is INTEGER
     blen(bidx) = i1len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%lat(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%lon(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%areacell(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! ------- casabal ----

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCgppyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCnppyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrmleafyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrmwoodyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrmrootyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrgrowyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrpyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrsyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCneeyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNdepyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNfixyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNsnetyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNupyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNleachyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNlossyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPweayear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPdustyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPsnetyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPupyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPleachyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPlossyear(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%glaimon(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (12, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = 12 * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%glaimonx(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (12, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = 12 * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%cplantlast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%nplantlast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%pplantlast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mplant * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%clitterlast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%nlitterlast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%plitterlast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%csoillast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%nsoillast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%psoillast(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1
     !blen(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%nsoilminlast(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%psoillablast(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%psoilsorblast(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%psoilocclast(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%cbalance(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%nbalance(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%pbalance(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%sumcbal(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%sumnbal(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%sumpbal(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%clabilelast(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! ------- phen -------

     bidx = bidx + 1
     CALL MPI_Get_address (phen%phase(off), displs(bidx), ierr)
     blen(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (phen%TKshed, displs(bidx), ierr)
     blen(bidx) = mvtype * extr2

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphase(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mphase, I1LEN, istride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%phen(off), displs(bidx), ierr)
     blen(bidx) =  r1len

     bidx = bidx + 1
     CALL MPI_Get_address (phen%aphen(off), displs(bidx), ierr)
     blen(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (phen%phasespin(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_1(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_3(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_4(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blen(bidx) = 1

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid number of casa_t param fields ',bidx,', fix it (02)!'
        CALL MPI_Abort(comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct(bidx, blen, displs, types, casa_t(rank), ierr)
     CALL MPI_Type_commit(casa_t(rank), ierr)

     CALL MPI_Type_size(casa_t(rank), tsize, ierr)
     CALL MPI_Type_get_extent(casa_t(rank), tmplb, text, ierr)

     WRITE(*,*) 'master to rank casa_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  END DO ! rank

  WRITE(*,*) 'total casa params size sent to all workers: ', localtotal
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blen)

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  remotetotal = 0
  CALL MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total casa params size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
     WRITE(*,*) 'error: total length of casa params sent and received differ'
     CALL MPI_Abort (comm, 0, ierr)
  END IF

  CALL MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  CALL master_send_input(comm, casa_t, 0)
  !  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank = 1, wnp
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

     bidx = bidx + 1
     call MPI_Get_address(met%fsd(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(swb, r1len, r1stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%tk(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%pmb(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%qv(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%ua(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%precip(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%precip_sn(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%fld(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%ca(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%coszen(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%Ndep(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! veg fields

     bidx = bidx + 1
     call MPI_Get_address(veg%vlai(off), displs(bidx), ierr)
     blocks(bidx) = r1len

    ! additional field missing from previous versions;
    ! added when trying to fix a bug in the new mpi code
    ! the order of these new fields follows the order of their
    ! declaration in cable_define_types.F90

     bidx = bidx + 1
     call MPI_Get_address(met%year(off), displs(bidx), ierr)
     blocks(bidx) = I1LEN

     bidx = bidx + 1
     call MPI_Get_address(met%moy(off), displs(bidx), ierr)
     blocks(bidx) = I1LEN

     bidx = bidx + 1
     call MPI_Get_address(met%doy(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%hod(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%u10(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     call MPI_Get_address(met%rhum(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid intype nmat, nvec or n3d constant, fix it (03)!'
        call MPI_Abort(comm, 1, ierr)
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
     call MPI_Abort(comm, 0, ierr)
  end if

  return

END SUBROUTINE master_intypes

! MPI: creates out_t types to receive the results from the workers
SUBROUTINE master_outtypes(comm,met,canopy,ssnow,rad,bal,air,soil,veg)

  use mpi

  USE cable_def_types_mod

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
     ! TODO: REAL(r_1) : rad%qcan(landunits,mf,nrb)
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
     !          &                        mat_t(midx, rank), ierr)
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
     !     &                        mat_t(midx, rank), ierr)
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

     ! additional  for sli
     midx = midx + 1
     CALL MPI_Get_address (ssnow%S(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address (ssnow%Tsoil(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address (ssnow%thetai(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     CALL MPI_Get_address (ssnow%snowliq(off,1), maddr(midx), ierr) ! 15
     CALL MPI_Type_create_hvector (3, r2len, r2stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

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
     !          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%cexpkdm(off,1), maddr(midx), ierr) ! 27
     ! Maciej: cexpkdm is mp*swb
     !     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
     !          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%rhocbm(off,1), maddr(midx), ierr) ! 27
     ! Maciej: rhocbm is mp*nrb
     !     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
     !          &                        mat_t(midx, rank), ierr)
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
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     ! ------------- 1D arrays -------------

     ! met
     vidx = 0

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%ca(off), vaddr(vidx), ierr) ! 1
     blen(vidx) = cnt * extr1

     ! gol124: not output, removed
     !vidx = vidx + 1
     ! INTEGER(i_d)
     !CALL MPI_Get_address (met%year(off), vaddr(vidx), ierr) ! 2
     !blen(vidx) = cnt * extid

     ! gol124: not output, removed
     !vidx = vidx + 1
     ! INTEGER(i_d)
     ! CALL MPI_Get_address (met%moy(off), vaddr(vidx), ierr) ! 3
     ! blen(vidx) = cnt * extid

     ! gol124: not output, removed
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (met%doy(off), vaddr(vidx), ierr) ! 4
     !blen(vidx) = cnt * extr1

     ! gol124: not output, removed
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (met%hod(off), vaddr(vidx), ierr) ! 5
     !blen(vidx) = cnt * extr1

     !vidx = vidx + 1
     ! REAL(r_1)
     ! MPI: gol124: changed to 2D and moved up when Bernard
     ! ported to CABLE_r491
     !CALL MPI_Get_address (met%fsd(off,1), vaddr(vidx), ierr)
     !blen(vidx) = cnt * extr1
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%fld(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%precip(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%precip_sn(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%tk(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%tvair(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%tvrad(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%pmb(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%ua(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%qv(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%qvair(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%da(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%dva(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (met%coszen(off), vaddr(vidx), ierr) ! 19
     blen(vidx) = cnt * extr1

     ! canopy
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fess(off), vaddr(vidx), ierr) ! 20
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fesp(off), vaddr(vidx), ierr) ! 20
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%cansto(off), vaddr(vidx), ierr) ! 20
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%cduv(off), vaddr(vidx), ierr) ! 21
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%delwc(off), vaddr(vidx), ierr) ! 22
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%dewmm(off), vaddr(vidx), ierr) ! 23
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (canopy%dgdtg(off), vaddr(vidx), ierr) ! 24
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fe(off), vaddr(vidx), ierr) ! 25
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fh(off), vaddr(vidx), ierr) ! 26
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fpn(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%A_sh(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%A_sl(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%A_slC(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%A_shC(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%A_slJ(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%A_shJ(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%GPP_sh(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%GPP_sl(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%eta_A_cs(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%dAdcs(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%cs(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%eta_GPP_cs(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2

     vidx = vidx + 1
     CALL MPI_Get_address (canopy%eta_fevc_cs(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     CALL MPI_Get_address (canopy%dlf(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr2

     vidx = vidx + 1
     CALL MPI_Get_address (veg%vcmax_shade(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     CALL MPI_Get_address (veg%vcmax_sun(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     CALL MPI_Get_address (veg%ejmax_shade(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     CALL MPI_Get_address (veg%ejmax_sun(off), vaddr(vidx), ierr) ! 27
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%frp(off), vaddr(vidx), ierr) ! 28
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%frpw(off), vaddr(vidx), ierr) ! 29
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%frpr(off), vaddr(vidx), ierr) ! 30
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%frs(off), vaddr(vidx), ierr) ! 31
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fnee(off), vaddr(vidx), ierr) ! 32
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%frday(off), vaddr(vidx), ierr) ! 33
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fnv(off), vaddr(vidx), ierr) ! 34
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fev(off), vaddr(vidx), ierr) ! 35
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (canopy%fevc(off), vaddr(vidx), ierr) ! 36
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (canopy%fevw(off), vaddr(vidx), ierr) ! 37
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! ! REAL(r_2)
     ! CALL MPI_Get_address (canopy%potev_c(off), vaddr(vidx), ierr) ! 38
     ! blen(vidx) = cnt * extr2
     ! vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fhv(off), vaddr(vidx), ierr) ! 39
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (canopy%fhvw(off), vaddr(vidx), ierr) ! 40
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fns(off), vaddr(vidx), ierr) ! 41
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fes(off), vaddr(vidx), ierr) ! 42
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fes_cor(off), vaddr(vidx), ierr) ! 42
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fhs(off), vaddr(vidx), ierr) ! 43
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fhs_cor(off), vaddr(vidx), ierr) ! 43
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fwet(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%epot(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fnpp(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%fevw_pot(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%gswx_T(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%cdtq(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%wetfac_cs(off), vaddr(vidx), ierr) ! 44
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%ga(off), vaddr(vidx), ierr) ! 45
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%ghflux(off), vaddr(vidx), ierr) ! 46
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%precis(off), vaddr(vidx), ierr) ! 47
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%qscrn(off), vaddr(vidx), ierr) ! 48
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%rnet(off), vaddr(vidx), ierr) ! 49
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%segg(off), vaddr(vidx), ierr) ! 50
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%sghflux(off), vaddr(vidx), ierr) ! 51
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%spill(off), vaddr(vidx), ierr) ! 52
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%through(off), vaddr(vidx), ierr) ! 53
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%tscrn(off), vaddr(vidx), ierr) ! 54
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%tv(off), vaddr(vidx), ierr) ! 55
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%us(off), vaddr(vidx), ierr) ! 56
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%uscrn(off), vaddr(vidx), ierr) ! 57
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%vlaiw(off), vaddr(vidx), ierr) ! 58
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%rghlai(off), vaddr(vidx), ierr) ! 58
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (canopy%wcint(off), vaddr(vidx), ierr) ! 59
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (canopy%fwsoil(off), vaddr(vidx), ierr) ! 59
     blen(vidx) = cnt * extr2
     
     ! 13C
     ! LOGICAL
     ! done as veg%deciduous
     vidx = vidx + 1
     CALL MPI_Get_address(canopy%isc3(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extl

     ! MPI: 2D vars moved above
     ! rwater
     ! evapfbl

     ! ssnow
     ! MPI: 2D vars moved above
     ! albsoilsn

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%pudsto(off), vaddr(vidx), ierr)
     blen(vidx) = r1len

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%pudsmx(off), vaddr(vidx), ierr)
     blen(vidx) = r1len

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%cls(off), vaddr(vidx), ierr) ! 60
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%dfn_dtg(off), vaddr(vidx), ierr) ! 61
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%dfh_dtg(off), vaddr(vidx), ierr) ! 62
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%dfe_ddq(off), vaddr(vidx), ierr) ! +1
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%ddq_dtg(off), vaddr(vidx), ierr) ! 63
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%evapsn(off), vaddr(vidx), ierr) ! 64
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%fwtop(off), vaddr(vidx), ierr) ! 65
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%fwtop1(off), vaddr(vidx), ierr)
     blen(vidx) = r1len

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%fwtop2(off), vaddr(vidx), ierr)
     blen(vidx) = r1len

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%fwtop3(off), vaddr(vidx), ierr)
     blen(vidx) = r1len

     ! MPI: 2D vars moved above
     ! gammzz
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (ssnow%isflag(off), vaddr(vidx), ierr) ! 66
     blen(vidx) = cnt * extid
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%osnowd(off), vaddr(vidx), ierr) ! 67
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%potev(off), vaddr(vidx), ierr) ! 68
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (soil%pwb_min(off), vaddr(vidx), ierr) ! 69
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%runoff(off), vaddr(vidx), ierr) ! 70
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%rnof1(off), vaddr(vidx), ierr) ! 71
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%rnof2(off), vaddr(vidx), ierr) ! 72
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%rtsoil(off), vaddr(vidx), ierr) ! 73
     blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! sconds
     ! sdepth
     ! smass
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%snage(off), vaddr(vidx), ierr) ! 74
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%snowd(off), vaddr(vidx), ierr) ! 75
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%smelt(off), vaddr(vidx), ierr) ! 76
     blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! dtmlt
     ! ssdn
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%ssdnn(off), vaddr(vidx), ierr) ! 77
     blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! tgg
     ! tggsn
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%tss(off), vaddr(vidx), ierr) ! 78
     blen(vidx) = cnt * extr1
     
     ! MPI: r1134 does not know about this field, comment out
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (ssnow%otss(off), vaddr(vidx), ierr) ! 79
     !blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! wb
     ! wbfice
     ! wbice
     ! wblf
     
     vidx = vidx + 1
     ! REAL(r_2)
     CALL MPI_Get_address (ssnow%wbtot(off), vaddr(vidx), ierr) ! 90
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%wb_lake(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%sinfil(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%qstss(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%wetfac(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     ! MPI: TODO: maybe not needed for transfer to master?
     CALL MPI_Get_address (ssnow%owetfac(off), vaddr(vidx), ierr) ! 92
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%t_snwlr(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%tggav(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%otss(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (ssnow%otss_0(off), vaddr(vidx), ierr) ! 91
     blen(vidx) = cnt * extr1

     ! rad
     ! MPI: 2D vars moved above
     ! albedo
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%extkb(off), vaddr(vidx), ierr) ! 93
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%extkd2(off), vaddr(vidx), ierr) ! 94
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%extkd(off), vaddr(vidx), ierr) ! 95
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%flws(off), vaddr(vidx), ierr) ! 96
     blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! fvlai
     ! gradis
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%latitude(off), vaddr(vidx), ierr) ! 97
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%lwabv(off), vaddr(vidx), ierr) !98
     blen(vidx) = cnt * extr1
     
     ! MPI: 3D vars moved above
     ! qcan
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%qssabs(off), vaddr(vidx), ierr) !99
     blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! rhocdf
     ! rniso
     ! scalex
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%transd(off), vaddr(vidx), ierr) ! 100
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%trad(off), vaddr(vidx), ierr) ! 101
     blen(vidx) = cnt * extr1

     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%transb(off), vaddr(vidx), ierr) ! 101
     blen(vidx) = cnt * extr1

     ! MPI: 2D vars moved above
     ! reffdf
     ! reffbm
     ! extkbm
     ! extkdm
     
     ! MPI: gol124: changed to 2D and moved up when Bernard
     ! ported to CABLE_r491
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (rad%fbeam(off,1), vaddr(vidx), ierr) ! 102
     !blen(vidx) = cnt * extr1
     ! MPI: 2D vars moved above
     ! cexpkbm
     ! cexpkdm

     ! bal
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (bal%drybal(off), vaddr(vidx), ierr) ! 103
     blen(vidx) = cnt * extr1
     
     ! MPI: remove ebal from exchanged data, calculate temp val on the master
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (bal%ebal(off), vaddr(vidx), ierr) ! 104
     !blen(vidx) = cnt * extr1
     ! MPI: remove ebal_tot from exchanged data, calculate val on the master
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (bal%ebal_tot(off), vaddr(vidx), ierr) ! 105
     !blen(vidx) = cnt * extr1
     ! MPI: remove seb from exchanged data, calculate temp val on the master
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (bal%seb(off), vaddr(vidx), ierr) ! 106
     !blen(vidx) = cnt * extr1
     ! MPI: remove seb_tot from exchanged data, calculate val on the master
     !vidx = vidx + 1
     ! REAL(r_1)
     !CALL MPI_Get_address (bal%seb_tot(off), vaddr(vidx), ierr) ! 107
     !blen(vidx) = cnt * extr1
     !vidx = vidx + 1
     ! REAL(r_1)
     ! MPI: remove evap_tot from exchanged data
     !CALL MPI_Get_address (bal%evap_tot(off), vaddr(vidx), ierr) ! 108
     !blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (bal%osnowd0(off), vaddr(vidx), ierr) ! 109
     blen(vidx) = cnt * extr1
     
     !vidx = vidx + 1
     ! REAL(r_1)
     ! MPI: remove precip_tot from exchanged data
     !CALL MPI_Get_address (bal%precip_tot(off), vaddr(vidx), ierr) ! 110
     !blen(vidx) = cnt * extr1
     !vidx = vidx + 1
     ! REAL(r_1)
     ! MPI: remove rnoff_tot from exchanged data
     !CALL MPI_Get_address (bal%rnoff_tot(off), vaddr(vidx), ierr) ! 111
     !blen(vidx) = cnt * extr1
     ! vidx = vidx + 1
     ! MPI: remove wbal from exchanged data
     ! REAL(r_1)
     ! CALL MPI_Get_address (bal%wbal(off), vaddr(vidx), ierr) ! 112
     !blen(vidx) = cnt * extr1
     !vidx = vidx + 1
     ! MPI: remove wbal_tot from exchanged data
     ! REAL(r_1)
     !CALL MPI_Get_address (bal%wbal_tot(off), vaddr(vidx), ierr) ! 113
     !blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (bal%wbtot0(off), vaddr(vidx), ierr) ! 114
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (bal%wetbal(off), vaddr(vidx), ierr) ! 115
     blen(vidx) = cnt * extr1

     ! air
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%rho(off), vaddr(vidx), ierr) ! 116
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%volm(off), vaddr(vidx), ierr) ! 117
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%rlam(off), vaddr(vidx), ierr) ! 118
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%qsat(off), vaddr(vidx), ierr) ! 119
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%epsi(off), vaddr(vidx), ierr) ! 120
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%visc(off), vaddr(vidx), ierr) ! 121
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%psyc(off), vaddr(vidx), ierr) ! 122
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%dsatdk(off), vaddr(vidx), ierr) ! 123
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (air%cmolar(off), vaddr(vidx), ierr) ! 124
     blen(vidx) = cnt * extr1

     ! soil
     ! MPI: 2D vars moved above
     ! albsoil
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%bch(off), vaddr(vidx), ierr) ! 125
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%c3(off), vaddr(vidx), ierr) ! 126
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%clay(off), vaddr(vidx), ierr) ! 127
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%cnsd(off), vaddr(vidx), ierr) ! 128
     blen(vidx) = cnt * extr2
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%css(off), vaddr(vidx), ierr) ! 129
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%hsbh(off), vaddr(vidx), ierr) ! 130
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%hyds(off), vaddr(vidx), ierr) ! 131
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (soil%i2bp3(off), vaddr(vidx), ierr) ! 132
     ! Maciej: i2bp3 is REAL
     !     blen(vidx) = cnt * extid
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (soil%ibp2(off), vaddr(vidx), ierr) ! 133
     ! Maciej: ibp2 is REAL
     !     blen(vidx) = cnt * extid
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (soil%isoilm(off), vaddr(vidx), ierr) ! 134
     blen(vidx) = cnt * extid
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%rhosoil(off), vaddr(vidx), ierr) ! 135
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%rs20(off), vaddr(vidx), ierr) ! 136
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%sand(off), vaddr(vidx), ierr) ! 137
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%sfc(off), vaddr(vidx), ierr) ! 138
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%silt(off), vaddr(vidx), ierr) ! 139
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%ssat(off), vaddr(vidx), ierr) ! 140
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%sucs(off), vaddr(vidx), ierr) ! 141
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (soil%swilt(off), vaddr(vidx), ierr) ! 142
     blen(vidx) = cnt * extr1

     ! veg
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (veg%iveg(off), vaddr(vidx), ierr) ! 143
     blen(vidx) = cnt * extid
     
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (veg%meth(off), vaddr(vidx), ierr) ! 144
     ! Maciej: veg%meth is REAL
     !     blen(vidx) = cnt * extid
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%vlai(off), vaddr(vidx), ierr) ! 145
     blen(vidx) = cnt * extr1
     
     ! MPI: 2D vars moved above
     ! froot
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%canst1(off), vaddr(vidx), ierr) ! 146
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%ejmax(off), vaddr(vidx), ierr) ! 147
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%frac4(off), vaddr(vidx), ierr) ! 148
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%wai(off), vaddr(vidx), ierr) ! 149
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%vegcf(off), vaddr(vidx), ierr) ! 150
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%tminvj(off), vaddr(vidx), ierr) ! 151
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%tmaxvj(off), vaddr(vidx), ierr) ! 152
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%vbeta(off), vaddr(vidx), ierr) ! 153
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%xalbnir(off), vaddr(vidx), ierr) ! 154
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%hc(off), vaddr(vidx), ierr) ! 155
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%shelrb(off), vaddr(vidx), ierr) ! 156
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%vcmax(off), vaddr(vidx), ierr) ! 157
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%xfang(off), vaddr(vidx), ierr) ! 158
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%dleaf(off), vaddr(vidx), ierr) ! 159
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%rp20(off), vaddr(vidx), ierr) ! 160
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%rpcoef(off), vaddr(vidx), ierr) ! 161
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (veg%extkn(off), vaddr(vidx), ierr) ! 162
     blen(vidx) = cnt * extr1
     
     vidx = vidx + 1
     ! LOGICAL
     CALL MPI_Get_address (veg%deciduous(off), vaddr(vidx), ierr) ! 163
     blen(vidx) = cnt * extl

      ! additional for SLI
     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%Tsurface(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr2

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%h0(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr2

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%delwcol(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr2

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%evap(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr2

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%nsnow(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extid

     vidx = vidx + 1
     CALL MPI_Get_address (ssnow%nsteps(off), vaddr(vidx), ierr)
     blen(vidx) = cnt * extr2

     ! end additional for SLI

     ! MPI: sanity check
     if (vidx /= nvec) then
        write(*,*) 'master: outtype invalid nvec ', vidx, ' constant, fix it (05)!'
        call MPI_Abort(comm, 1, ierr)
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
          CALL MPI_Abort (comm, 0, ierr)
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
SUBROUTINE master_casa_types(comm, casapool, casaflux, casamet, casabal, phen)

  use mpi

  USE cable_def_types_mod
  USE casadimension
  USE casavariable
  USE phenvariable

  IMPLICIT NONE

  INTEGER :: comm ! MPI communicator to talk to the workers

  TYPE(casa_pool),     INTENT(INOUT) :: casapool
  TYPE(casa_flux),     INTENT(INOUT) :: casaflux
  TYPE(casa_met),      INTENT(INOUT) :: casamet
  TYPE(casa_balance),  INTENT(INOUT) :: casabal
  TYPE(phen_variable), INTENT(INOUT) :: phen

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, I1LEN
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, ierr

  ALLOCATE (casa_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = ncasa_mat + ncasa_vec + (nphen - 1)
  ALLOCATE(blocks(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  r1stride = mp * extr1
  r2stride = mp * extr2
  istride = mp * extid
  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     I1LEN = cnt * extid

     bidx = 0

     ! ------------- 2D arrays -------------

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%cplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%clitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%csoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nlitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%pplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%plitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%psoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioPCplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNClitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioPClitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioNCsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%ratioPCsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphase(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mphase, I1LEN, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%phasespin(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_1(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_3(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (phen%doyphasespin_4(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mdyear, I1LEN, istride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cplant_turnover(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracCalloc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracNalloc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracPalloc(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crmplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%kplant_fire(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromPtoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromLtoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromStoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(msoil, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%kplant_tot(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%klitter_tot(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fluxfromPtoCO2_fire(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fluxfromLtoCO2_fire(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! ------------- 3D vectors -------------
     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromPtoL(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant * mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%klitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%klitter_fire(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%ksoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * r2len

     ! gol124: temp only
     ! 3D
     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromLtoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil * mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * mlitter * r2len

     ! gol124: temp only
     ! 3D
     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromStoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil * msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromLtoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fromStoCO2(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxCtolitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxNtolitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxPtolitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = mlitter * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxCtosoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxNtosoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fluxPtosoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, &
          &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = msoil * r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromPtoL(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mplant*mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromLtoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter*msoil, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromStoS(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(msoil*msoil, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fromPtoL_fire(off,1,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mlitter*mplant, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! ------------- 1D vectors -------------

     last2d = bidx

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%glai(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (phen%phase(off), displs(bidx), ierr)
     blocks(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address (phen%phen(off), displs(bidx), ierr)
     blocks(bidx) =  r1len

     bidx = bidx + 1
     CALL MPI_Get_address (phen%aphen(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%clabile(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Ctot(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Ctot_0(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nsoilmin(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%psoillab(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%psoilsorb(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%psoilocc(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%psorbmax(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%sumcbal(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%sumnbal(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%sumpbal(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCgppyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCnppyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrmleafyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrmwoodyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrmrootyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrgrowyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrpyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCrsyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCneeyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNdepyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNfixyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNsnetyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNupyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNleachyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FNlossyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPweayear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPdustyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPsnetyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPupyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPleachyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FPlossyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cgpp(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cnpp(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crp(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crgplant(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminfix(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminuptake(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Plabuptake(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Clabloss(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fracClabile(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cnep(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Crsoil(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nmindep(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminloss(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nminleach(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nupland(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nlittermin(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nsmin(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nsimm(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Nsnet(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%frac_sapwood(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%sapwood_area(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cplant_turnover_disturbance(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cplant_turnover_crowding(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%Cplant_turnover_resource_limitation(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fHarvest(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fCrop(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxFromPtoHarvest(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fluxCtoCO2_plant_fire(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fluxCtoCO2_litter_fire(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fluxNtoAtm_fire(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%stemnpp(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%fNminloss(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%FluxCtoCO2(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid number of casa fields, fix it (06)!'
        WRITE(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, casa_ts(rank), ierr)
     CALL MPI_Type_commit (casa_ts(rank), ierr)

     CALL MPI_Type_size (casa_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (casa_ts(rank), tmplb, text, ierr)

     WRITE(*,*) 'casa results recv from worker, size, extent, lb: ', rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     ! TODO: also free partial types for intypes, outtypes etc.
     DO i = 1, last2d
        CALL MPI_Type_free(types(i), ierr)
     END DO

  END DO

  WRITE(*,*) 'total size of casa results received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total size of casa results sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
     WRITE(*,*) 'error: casa results totalsend and totalrecv differ'
     CALL MPI_Abort(comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN

END SUBROUTINE master_casa_types

SUBROUTINE master_climate_types (comm, climate, ktauday)

  use mpi

  USE cable_def_types_mod, ONLY: climate_type, mp
  USE cable_climate_mod,   ONLY: climate_init, READ_CLIMATE_RESTART_NC

  TYPE(climate_type):: climate
  INTEGER :: comm, ktauday
  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, I1LEN
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, ierr, ny, nd, ndq, nsd

!!$  CALL climate_init (climate, mp, ktauday)
!!$  if (cable_user%call_climate .AND.(.NOT.cable_user%climate_fromzero)) &
!!$       CALL READ_CLIMATE_RESTART_NC (climate, ktauday)
  ALLOCATE (climate_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nclimate

  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  r1stride = mp * extr1
  r2stride = mp * extr2

  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

 DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2
     I1LEN = cnt * extid

     bidx = 0

     ! ------------- 2D arrays -------------
     
     ny = climate%nyear_average
     nd = climate%nday_average
     ndq = 91
     nsd = 5 * ktauday

     ! write(*,*) 'master, nd ny mp nsd', nd, ny,mp, nsd

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp_min_20(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ny, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     ! write(*,*) 'master', blocks, bidx, ny, r1len, displs, ierr, climate%mtemp_min_20(off,:)

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp_max_20(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ny, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%alpha_PT_20(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ny, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1


     bidx = bidx + 1
     CALL MPI_Get_address (climate%dtemp_31(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dtemp_91(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ndq, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1


     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_31(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1


     bidx = bidx + 1
     CALL MPI_Get_address (climate%APAR_leaf_sun(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%APAR_leaf_shade(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%Dleaf_sun(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%Dleaf_shade(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%Tleaf_sun(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%Tleaf_shade(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%cs_sun(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%cs_shade(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%scalex_sun(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%scalex_shade(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%fwsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (nsd, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_min_20(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ny, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_max_20(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ny, r1len, r1stride, MPI_BYTE, &
     &                                types(bidx), ierr)
     blocks(bidx) = 1
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%aprecip_20(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (ny, r1len, r1stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! ------------- 1D vectors -------------

     bidx = bidx + 1
     CALL MPI_Get_address (climate%chilldays(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%iveg(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%biome(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%GMD(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dtemp(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%qtemp(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mmoist(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp_min(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%qtemp_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%qtemp_max_last_year(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp_min20(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%mtemp_max20(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%atemp_mean(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%agdd5(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%gdd5(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%agdd0(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%gdd0(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%alpha_PT(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%evap_PT(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%aevap(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%frec(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%gdd0_rec(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%fdorm(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_min20(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_max20(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%fapar_ann_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%fapar_ann_max_last_year(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%modis_igbp(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%AvgAnnMaxFAPAR(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%DSLR(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%aprecip_av20(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%NDAY_Nesterov(off), displs(bidx), ierr)
     blocks(bidx) = i1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_min(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%dmoist_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%alpha_PT20(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%dtemp_min(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%dtemp_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%drhum(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%du10_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%dprecip(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%aprecip(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%last_precip(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%KBDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%FFDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%D_MacArthur(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%Nesterov_Current(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%Nesterov_ann_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%Nesterov_ann_max_last_year(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE
     
     bidx = bidx + 1
     CALL MPI_Get_address (climate%Nesterov_ann_running_max(off), displs(bidx), ierr)
     blocks(bidx) = r1len
     types(bidx)  = MPI_BYTE

     ! ------------- scalars  -------------

     bidx = bidx + 1
     CALL MPI_Get_address (climate%nyears, displs(bidx), ierr)
     blocks(bidx) = extid
     types(bidx)  = MPI_BYTE

     bidx = bidx + 1
     CALL MPI_Get_address (climate%doy, displs(bidx), ierr)
     blocks(bidx) = extid
     types(bidx)  = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid number of climate fields, fix it (07)!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, climate_ts(rank), ierr)
     CALL MPI_Type_commit (climate_ts(rank), ierr)

     CALL MPI_Type_size (climate_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (climate_ts(rank), tmplb, text, ierr)

     WRITE(*,*) 'climate results recv from worker, size, extent, lb: ', rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

  END DO

  WRITE(*,*) 'total size of climate results received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE(*,*) 'total size of climate results sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
     WRITE(*,*) 'error: climate results totalsend and totalrecv differ'
     CALL MPI_Abort(comm, 0, ierr)
  END IF

  DO rank = 1, wnp
     CALL MPI_ISend (MPI_BOTTOM, 1, climate_ts(rank), rank, 0, comm, inp_req(rank), ierr)
  END DO

  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN

END SUBROUTINE master_climate_types

!MPI

!CLNSUBROUTINE master_casa_restart_types( comm, casamet, casapool )
!CLN
!CLN    use mpi
!CLN
!CLN  USE casavariable, ONLY: casa_met, casa_pool
!CLN
!CLN  IMPLICIT NONE
!CLN
!CLN  INTEGER,INTENT(IN) :: comm
!CLN  TYPE (casa_met)    , INTENT(OUT) :: casamet
!CLN  TYPE (casa_pool)   , INTENT(OUT) :: casapool
!CLN
!CLN  ! local vars
!CLN
!CLN  ! temp arrays for marshalling all fields into a single struct
!CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
!CLN  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
!CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
!CLN
!CLN  ! temp vars for verifying block number and total length of inp_t
!CLN  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
!CLN  INTEGER :: tsize, localtotal, remotetotal
!CLN
!CLN  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride
!CLN  INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
!CLN  INTEGER :: bidx ! block index
!CLN  INTEGER :: ntyp ! total number of blocks
!CLN  INTEGER :: rank ! worker rank
!CLN  INTEGER :: off  ! first patch index for a worker
!CLN  INTEGER :: cnt  ! mp for a worker
!CLN  INTEGER :: ierr
!CLN
!CLN  ALLOCATE (casa_dump_ts(wnp))
!CLN
!CLN  ntyp = ncdumprw
!CLN
!CLN  ALLOCATE (blocks(ntyp))
!CLN  ALLOCATE (displs(ntyp))
!CLN  ALLOCATE (types(ntyp))
!CLN
!CLN  ! chunks of all 1D vectors are contiguous blocks of memory so just send them
!CLN  ! as blocks of bytes
!CLN  types = MPI_BYTE
!CLN
!CLN  ! total size of input data sent to all workers
!CLN  localtotal = 0
!CLN
!CLN  ! create a separate MPI derived datatype for each worker
!CLN  DO rank = 1, wnp
!CLN
!CLN     ! starting patch and number for each worker rank
!CLN     off = wland(rank)%patch0
!CLN     cnt = wland(rank)%npatch
!CLN
!CLN     r1len = cnt * extr1
!CLN
!CLN     r1stride = mp * extr1
!CLN
!CLN     ! casamet fields
!CLN
!CLN     bidx = 1
!CLN     CALL MPI_Get_address (casamet%glai(off), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len
!CLN
!CLN     ! casapool fields
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%clabile(off), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%psoillab(off), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%psoilsorb(off), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%psoilocc(off), displs(bidx), ierr)
!CLN     blocks(bidx) = r1len
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%cplant(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%nplant(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%pplant(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%cplant(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%clitter(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%nlitter(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%csoil(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%nsoil(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     bidx = bidx + 1
!CLN     CALL MPI_Get_address (casapool%psoil(off,1), displs(bidx), ierr)
!CLN     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
!CLN          types(bidx), ierr)
!CLN     blocks(bidx) = 1
!CLN
!CLN     ! MPI: sanity check
!CLN     IF (bidx /= ntyp) THEN
!CLN        WRITE (*,*) 'master: invalid intype in master_casa_dump, fix it (08)!'
!CLN        CALL MPI_Abort (comm, 1, ierr)
!CLN     END IF
!CLN
!CLN     ! marshall all fields into a single MPI derived datatype for worker rank
!CLN
!CLN     CALL MPI_Type_create_struct (bidx, blocks, displs, types, casa_restart_ts(rank), ierr)
!CLN     CALL MPI_Type_commit (casa_restart_ts(rank), ierr)
!CLN
!CLN     CALL MPI_Type_size (casa_restart_ts(rank), tsize, ierr)
!CLN     CALL MPI_Type_get_extent (casa_restart_ts(rank), tmplb, text, ierr)
!CLN
!CLN     WRITE (*,*) 'master to ',rank,': intype struct blocks, size, extent and lb: ', &
!CLN                 bidx,tsize,text,tmplb
!CLN
!CLN     localtotal = localtotal + tsize
!CLN
!CLN  END DO ! rank
!CLN
!CLN  DEALLOCATE(types)
!CLN  DEALLOCATE(displs)
!CLN  DEALLOCATE(blocks)
!CLN
!CLN  WRITE (*,*) 'total casa_dump size sent to all workers: ', localtotal
!CLN
!CLN  ! MPI: check whether total size of send input data equals total
!CLN  ! data received by all the workers
!CLN  remotetotal = 0
!CLN  CALL MPI_Reduce (MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
!CLN
!CLN  WRITE (*,*) 'total input data size received by all workers: ', remotetotal
!CLN
!CLN  IF (localtotal /= remotetotal) THEN
!CLN          WRITE (*,*) 'error: total length of input data sent and received differ (02)'
!CLN          CALL MPI_Abort (comm, 0, ierr)
!CLN  END IF
!CLN
!CLNEND SUBROUTINE master_casa_restart_types

! MPI: creates datatype handles to receive restart data from workers
SUBROUTINE master_restart_types(comm, canopy, air)

  use mpi

  USE cable_def_types_mod
  USE casadimension
  USE casavariable

  IMPLICIT NONE

  INTEGER :: comm ! MPI communicator to talk to the workers

  TYPE(canopy_type), INTENT(IN) :: canopy
  TYPE (air_type),INTENT(IN)        :: air
!  TYPE (casa_pool),           INTENT(INOUT) :: casapool
!  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
!  TYPE (casa_met),            INTENT(INOUT) :: casamet
!  TYPE (casa_balance),        INTENT(INOUT) :: casabal

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, ierr

  ALLOCATE(restart_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nrestart
  ALLOCATE(blocks(ntyp))
  ALLOCATE(displs(ntyp))
  ALLOCATE(types(ntyp))

  r1stride = mp * extr1
  r2stride = mp * extr2

  ! counter to sum total number of bytes receives from all workers
  totalrecv = 0

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx = 0

     ! ------------- 2D arrays -------------

     ! bidx = bidx + 1
     ! CALL MPI_Get_address (canopy%rwater(off,1), displs(bidx), ierr) ! 1
     ! CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
     ! &                             types(bidx), ierr)
     ! blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(canopy%evapfbl(off,1), displs(bidx), ierr) ! 2
     ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
     CALL MPI_Type_create_hvector(ms, r1len, r1stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     last2d = bidx

     ! ------------- 1D vectors -------------

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%cduv(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%dewmm(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%dgdtg(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frpw(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%frpr(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%fnv(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%rho(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%volm(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%qsat(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%epsi(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%visc(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%psyc(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%dsatdk(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (air%cmolar(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     types(last2d+1:bidx) = MPI_BYTE


     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid number of restart fields, fix it (09)!'
        WRITE(*,*) 'bidx: ', bidx
        WRITE(*,*) 'ntyp: ', ntyp
        CALL MPI_Abort(comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct(bidx, blocks, displs, types, restart_ts(rank), ierr)
     CALL MPI_Type_commit(restart_ts(rank), ierr)

     CALL MPI_Type_size(restart_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent(restart_ts(rank), tmplb, text, ierr)

     WRITE(*,*) 'restart results recv from worker, size, extent, lb: ', rank, tsize, text, tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free(types(i), ierr)
     END DO

  END DO

  WRITE(*,*) 'total size of restart fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  ! write(*,*) 'b5 reduce wk', MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr
  ! call flush(6)
  CALL MPI_Reduce(MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE(*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
     WRITE(*,*) 'error: restart fields totalsend and totalrecv differ'
     CALL MPI_Abort(comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN
  
END SUBROUTINE master_restart_types


! MPI: Casa - dump read and write
SUBROUTINE master_casa_dump_types(comm, casamet, casaflux, phen, climate, c13o2flux )

  use mpi

  use casadimension,       only: icycle
  use casavariable,        only: casa_met, casa_flux
  use cable_def_types_mod, only: climate_type
  use phenvariable
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

     bidx = 1
     CALL MPI_Get_address(casamet%tairk(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casamet%tsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(ms, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address(casamet%moist(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(ms, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! casaflux fields

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%cgpp(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address(casaflux%crmplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(ncp, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! phen fields

     bidx = bidx + 1
     CALL MPI_Get_address(phen%phase(off), displs(bidx), ierr)
     blocks(bidx) = I1LEN

     bidx = bidx + 1
     CALL MPI_Get_address(phen%doyphase(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector(mphase, I1LEN, istride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! climate fields

     bidx = bidx + 1
     CALL MPI_Get_address(climate%qtemp_max_last_year(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! N and P deposition
     if (icycle>1) then
        bidx = bidx + 1
        CALL MPI_Get_address(casaflux%Nmindep(off), displs(bidx), ierr)
        blocks(bidx) = r2len
     endif
     
     if (icycle>2) then
        bidx = bidx + 1
        CALL MPI_Get_address(casaflux%Pdep(off), displs(bidx), ierr)
        blocks(bidx) = r2len
     endif

     ! 13C
     ! c13o2 fields
     
     if (cable_user%c13o2) then
        bidx = bidx + 1
        CALL MPI_Get_address(c13o2flux%cAn12(off), displs(bidx), ierr)
        blocks(bidx) = r2len

        bidx = bidx + 1
        CALL MPI_Get_address(c13o2flux%cAn(off), displs(bidx), ierr)
        blocks(bidx) = r2len
     endif
     
     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid intype in master_casa_dump, fix it (10)!'
        CALL MPI_Abort(comm, 1, ierr)
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
     CALL MPI_Abort(comm, 0, ierr)
  END IF

  !  DO rank = 1, wnp
  !
  !     CALL MPI_ISend (MPI_BOTTOM, 1, casa_dump_ts(rank), rank, 0, comm, &
  !          &               inp_req(rank), ierr)
  !
  !  END DO
  !
  !  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

END SUBROUTINE master_casa_dump_types

! #############################################################################################################

! MPI: Casa-LUC: exchanging casapools between master and worker, as required for LUC updates
SUBROUTINE master_casa_LUC_types(comm, casapool, casabal, casaflux)

  use mpi

  USE casavariable, ONLY: casa_pool, casa_balance, casa_flux, mplant, mlitter, msoil
  
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

     ! casapool fields
     bidx = bidx + 1
     CALL MPI_Get_address (casapool%cplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%clitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%csoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nlitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%nsoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%pplant(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mplant, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%plitter(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%psoil(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (msoil, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     last2d = bidx

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%Nsoilmin(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casapool%clabile(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     ! casabal fields
     bidx = bidx + 1
     CALL MPI_Get_address (casabal%FCneeyear(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     ! casaflux fields
     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fHarvest(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%cHarvest(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%nHarvest(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     CALL MPI_Get_address (casaflux%fcrop(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     types(last2d+1:bidx) = MPI_BYTE


     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE(*,*) 'master: invalid intype in master_casa_LUC, fix it (11)!'
        CALL MPI_Abort(comm, 1, ierr)
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
     CALL MPI_Abort(comm, 0, ierr)
  END IF

END SUBROUTINE master_casa_LUC_types

! *******************************************************************************************

! MPI
! Creates pop_ts types to broadcast/cscatter the default POP parameters
! to all workers

SUBROUTINE master_pop_types(comm, casamet, pop)

  use mpi
  USE POP_mpi
  USE POP_types,          ONLY: pop_type
  USE cable_common_module,ONLY: cable_user
  USE casavariable,       ONLY: casa_met

  IMPLICIT NONE

  INTEGER,        INTENT(IN)    :: comm
  TYPE(casa_met), INTENT(IN)    :: casamet
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

!     CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ENDIF

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
  ! NO Waitall here as some workers might not be involved!!!

END SUBROUTINE master_receive_pop

!!CLNSUBROUTINE master_blaze_types (comm, BLAZE)
!!CLN
!!CLN  ! Send blaze restart data to workers
!!CLN
!!CLN  use mpi
!!CLN
!!CLN  USE blaze, ONLY: TYPE_BLAZE
!!CLN
!!CLN  IMPLICIT NONE
!!CLN
!!CLN  INTEGER :: comm ! MPI communicator to talk to the workers
!!CLN
!!CLN  TYPE(TYPE_BLAZE), INTENT(IN) :: BLAZE
!!CLN
!!CLN
!!CLN
!!CLN  ! MPI: temp arrays for marshalling all types into a struct
!!CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
!!CLN  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
!!CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
!!CLN  INTEGER :: ntyp ! number of worker's types
!!CLN
!!CLN  INTEGER :: last2d, i
!!CLN
!!CLN  ! MPI: block lenghts for hindexed representing all vectors
!!CLN  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
!!CLN
!!CLN  ! MPI: block lengths and strides for hvector representing matrices
!!CLN  INTEGER :: r1len, r2len
!!CLN  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride
!!CLN
!!CLN  INTEGER :: tsize, totalrecv, totalsend
!!CLN  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
!!CLN
!!CLN  INTEGER :: rank, off, cnt
!!CLN  INTEGER :: bidx, midx, vidx, ierr
!!CLN
!!CLN
!!CLN  ! Restart value handles (restart_blaze_ts())
!!CLN
!!CLN  ntyp = 9
!!CLN
!!CLN  !INTEGER,  DIMENSION(:),  ALLOCATABLE :: DSLR,
!!CLN  !REAL,     DIMENSION(:),  ALLOCATABLE :: RAINF, KBDI, LR
!!CLN  !REAL,     DIMENSION(:,:),ALLOCATABLE :: AnnRAINF, DEADWOOD, AGLB_w, AGLB_g, AGLit_w, AGLit_g
!!CLN
!!CLN  ALLOCATE (blocks(ntyp))
!!CLN  ALLOCATE (displs(ntyp))
!!CLN  ALLOCATE (types(ntyp))
!!CLN
!!CLN  istride  = mp * extid ! short integer
!!CLN  r1stride = mp * extr1 ! single precision
!!CLN  r2stride = mp * extr2 ! double precision
!!CLN
!!CLN  DO rank = 1, wnp
!!CLN     off = wland(rank)%patch0
!!CLN     cnt = wland(rank)%npatch
!!CLN
!!CLN     i1len = cnt * extid
!!CLN     r1len = cnt * extr1
!!CLN     r2len = cnt * extr2
!!CLN
!!CLN     bidx = 0
!!CLN
!!CLN     ! ------------- 2D arrays -------------
!!CLN
!!CLN     ! Annual (daily) rainfall (ncells,366)
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%AnnRainf(off,1), displs(bidx), ierr) ! 1
!!CLN     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
!!CLN     &                             types(bidx), ierr)
!!CLN     blocks(bidx) = 1
!!CLN
!!CLN     ! Above ground life woody biomass
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr) ! 2
!!CLN     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
!!CLN     &                             types(bidx), ierr)
!!CLN     blocks(bidx) = 1
!!CLN
!!CLN     ! Above ground life grassy biomass
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr) ! 3
!!CLN     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
!!CLN     &                             types(bidx), ierr)
!!CLN     blocks(bidx) = 1
!!CLN
!!CLN     ! Above ground woody litter
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr) ! 4
!!CLN     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
!!CLN     &                             types(bidx), ierr)
!!CLN     blocks(bidx) = 1
!!CLN
!!CLN     ! Above ground grassy litter
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr) ! 5
!!CLN     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
!!CLN     &                             types(bidx), ierr)
!!CLN     blocks(bidx) = 1
!!CLN
!!CLN     last2d = bidx
!!CLN
!!CLN     ! ------------- 1D vectors -------------
!!CLN
!!CLN     ! Integer days since last rainfall
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
!!CLN     blocks(bidx) = i1len
!!CLN
!!CLN     ! Real(sp) Last rainfall
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
!!CLN     blocks(bidx) = r1len
!!CLN
!!CLN     ! current KBDI
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
!!CLN     blocks(bidx) = r1len
!!CLN
!!CLN     ! DEADWOOD
!!CLN     bidx = bidx + 1
!!CLN     CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
!!CLN     blocks(bidx) = r1len
!!CLN
!!CLN     ! ------------- Wrap up -------------
!!CLN
!!CLN     types(last2d+1:bidx) = MPI_BYTE
!!CLN
!!CLN     ! MPI: sanity check
!!CLN     IF (bidx /= ntyp) THEN
!!CLN        WRITE (*,*) 'invalid blaze ntyp constant, fix it (12)!'
!!CLN        CALL MPI_Abort (comm, 1, ierr)
!!CLN     END IF
!!CLN
!!CLN     CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_restart_ts(rank), ierr)
!!CLN     CALL MPI_Type_commit (blaze_restart_ts(rank), ierr)
!!CLN
!!CLN     CALL MPI_Type_size (blaze_restart_ts(rank), tsize, ierr)
!!CLN     CALL MPI_Type_get_extent (blaze_restart_ts(rank), tmplb, text, ierr)
!!CLN
!!CLN     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
!!CLN   &       rank,tsize,text,tmplb
!!CLN
!!CLN     totalrecv = totalrecv + tsize
!!CLN
!!CLN     ! free the partial types used for matrices
!!CLN     DO i = 1, last2d
!!CLN        CALL MPI_Type_free (types(i), ierr)
!!CLN     END DO
!!CLN
!!CLN  END DO
!!CLN
!!CLN  WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv
!!CLN
!!CLN  ! MPI: check whether total size of received data equals total
!!CLN  ! data sent by all the workers
!!CLN  totalsend = 0
!!CLN  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
!!CLN    &     0, comm, ierr)
!!CLN
!!CLN  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend
!!CLN
!!CLN  IF (totalrecv /= totalsend) THEN
!!CLN          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
!!CLN          CALL MPI_Abort (comm, 0, ierr)
!!CLN  END IF
!!CLN
!!CLN  DEALLOCATE(types)
!!CLN  DEALLOCATE(displs)
!!CLN  DEALLOCATE(blocks)
!!CLN
!!CLN
!!CLN  ! Standard desired IO (restart_blaze_ts())
!!CLN
!!CLN
!!CLN
!!CLN
!!CLN
!!CLN  RETURN
!!CLN
!!CLNEND SUBROUTINE master_blaze_types
!!CLN

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
  use cable_mpicommon,     only: nc13o2_flux
  
  implicit none

  integer,          intent(in)    :: comm  ! mpi communicator
  type(c13o2_flux), intent(inout) :: c13o2flux

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blen
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

  allocate(blen(ntyp))
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
     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%ca(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%cAn12(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%cAn(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%RAn(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Vstarch(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Rstarch(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! 2D
     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%An(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Disc(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Rsucrose(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Rphoto(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_flux_t param fields ', bidx, ', fix it (13)!'
        call MPI_Abort(comm, 1, ierr)
     end if

     call MPI_Type_create_struct(bidx, blen, displs, types, c13o2_flux_t(rank), ierr)
     call MPI_Type_commit(c13o2_flux_t(rank), ierr)

     call MPI_Type_size(c13o2_flux_t(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_flux_t(rank), tmplb, text, ierr)

     write(*,*) 'master to rank c13o2_flux_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  end do ! rank

  write(*,*) 'total c13o2_flux params size sent to all workers: ', localtotal
  deallocate(types)
  deallocate(displs)
  deallocate(blen)

  ! MPI: check whether total size of received data equals total data sent by all the workers
  remotetotal = 0
  call MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total c13o2_flux params size received by all workers: ', remotetotal

  if (localtotal /= remotetotal) then
     write(*,*) 'error: total length of c13o2_flux params sent and received differ'
     call MPI_Abort(comm, 0, ierr)
  end if

  call MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  call master_send_input(comm, c13o2_flux_t, 0)
  ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  do rank = 1, wnp
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
  use cable_mpicommon,     only: nc13o2_pool
  
  implicit none

  integer,          intent(in)    :: comm  ! mpi communicator
  type(c13o2_pool), intent(inout) :: c13o2pools

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blen
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

  allocate(blen(ntyp))
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
     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%clabile(off), displs(bidx), ierr)
     blen(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%charvest(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! 2D
     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%cplant(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mplant, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%clitter(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%csoil(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(msoil, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_pool_t param fields ', bidx, ', fix it (14)!'
        call MPI_Abort(comm, 1, ierr)
     end if

     call MPI_Type_create_struct(bidx, blen, displs, types, c13o2_pool_t(rank), ierr)
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
     call MPI_Abort(comm, 0, ierr)
  end if

  deallocate(types)
  deallocate(displs)
  deallocate(blen)

  call MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  call master_send_input(comm, c13o2_pool_t, 0)
  ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  do rank = 1, wnp
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
  use cable_c13o2_def,     only: c13o2_luc
  use cable_mpicommon,     only: nc13o2_luc
  
  implicit none

  integer,         intent(in)    :: comm  ! mpi communicator
  type(c13o2_luc), intent(inout) :: c13o2luc

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  integer, allocatable, dimension(:) :: blen
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
  allocate(blen(ntyp))
  allocate(displs(ntyp))
  allocate(types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mland * extr1
  r2stride = mland * extr2
  istride  = mland * extid

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
     bidx = bidx + 1
     call MPI_Get_address(c13o2luc%cagric(off), displs(bidx), ierr)
     blen(bidx) = r2len

     ! 2D
     bidx = bidx + 1
     call MPI_Get_address(c13o2luc%charvest(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(c13o2luc%nharvest, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2luc%cclearance(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(c13o2luc%nclearance, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blen(bidx) = 1

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_luc_t param fields ', bidx, ', fix it (15)!'
        call MPI_Abort(comm, 1, ierr)
     end if

     call MPI_Type_create_struct(bidx, blen, displs, types, c13o2_luc_t(rank), ierr)
     call MPI_Type_commit(c13o2_luc_t(rank), ierr)

     call MPI_Type_size(c13o2_luc_t(rank), tsize, ierr)
     call MPI_Type_get_extent(c13o2_luc_t(rank), tmplb, text, ierr)

     write(*,*) 'master to rank c13o2_luc_t param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

     localtotal = localtotal + tsize

  end do ! rank

  write(*,*) 'total c13o2_luc params size sent to all workers: ', localtotal
  deallocate(types)
  deallocate(displs)
  deallocate(blen)

  ! MPI: check whether total size of received data equals total data sent by all the workers
  remotetotal = 0
  call MPI_Reduce(MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  write(*,*) 'total c13o2_luc params size received by all workers: ', remotetotal

  if (localtotal /= remotetotal) then
     write(*,*) 'error: total length of c13o2_luc params sent and received differ'
     call MPI_Abort(comm, 0, ierr)
  end if

  call MPI_Barrier(comm, ierr)

  ! so, now send all the parameters
  call master_send_input(comm, c13o2_luc_t, 0)
  ! call MPI_Waitall(wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  do rank = 1, wnp
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
  use cable_mpicommon,     only: nc13o2_flux

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

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%An(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Disc(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Rsucrose(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Rphoto(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mf, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! ------------- 1D vectors -------------

     last2d = bidx

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%ca(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%cAn12(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%cAn(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%RAn(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Vstarch(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2flux%Rstarch(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_flux fields, fix it (16)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 1, ierr)
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
     call MPI_Abort(comm, 0, ierr)
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
  use cable_mpicommon,     only: nc13o2_pool

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

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%cplant(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mplant, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%clitter(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(mlitter, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%csoil(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(msoil, r2len, r2stride, MPI_BYTE, &
          types(bidx), ierr)
     blocks(bidx) = 1

     ! ------------- 1D vectors -------------

     last2d = bidx

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%clabile(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     bidx = bidx + 1
     call MPI_Get_address(c13o2pools%charvest(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_pool fields, fix it (17)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 1, ierr)
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
     call MPI_Abort(comm, 0, ierr)
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
  use cable_c13o2_def,     only: c13o2_luc
  use cable_mpicommon,     only: nc13o2_luc

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

     bidx = bidx + 1
     call MPI_Get_address(c13o2luc%charvest(off,1), displs(bidx), ierr)
     ! get new types
     call MPI_Type_create_hvector(c13o2luc%nharvest, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     bidx = bidx + 1
     call MPI_Get_address(c13o2luc%cclearance(off,1), displs(bidx), ierr)
     call MPI_Type_create_hvector(c13o2luc%nclearance, r2len, r2stride, MPI_BYTE, types(bidx), ierr)
     blocks(bidx) = 1

     ! ------------- 1D vectors -------------

     last2d = bidx

     bidx = bidx + 1
     call MPI_Get_address(c13o2luc%cagric(off), displs(bidx), ierr)
     blocks(bidx) = r2len

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     if (bidx /= ntyp) then
        write(*,*) 'master: invalid number of c13o2_luc fields, fix it (18)!'
        write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
        call MPI_Abort(comm, 1, ierr)
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
     call MPI_Abort(comm, 0, ierr)
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
!            &              canopy,rough,rad,bgc,bal)
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
        endif
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
     endif
  END IF

  ! MPI: free partial derived datatype handle arrays
  DEALLOCATE(vec_t)
  DEALLOCATE(mat_t)
  DEALLOCATE(m3d_t)

  ! MPI: free landpoint decomposition info
  DEALLOCATE(wland)

  RETURN

END SUBROUTINE master_end

! 13C
SUBROUTINE master_spincasacnp(dels, kstart, kend, mloop, veg, soil, casabiome, casapool, &
     casaflux, casamet, casabal, phen, POP, climate, LALLOC, c13o2flux, c13o2pools, c13o2luc, icomm, ocomm)

  use cable_def_types_mod
  use cable_carbon_module
  use cable_common_module, only: cable_user
  use casadimension
  use casaparm
  use casavariable
  use phenvariable
  use POP_types,           only: POP_type
  use POPmodule,           only: POPstep
  ! 13C
  use cable_c13o2_def,     only: c13o2_flux, c13o2_pool, c13o2_luc
  use cable_c13o2,         only: c13o2_write_restart_pools, c13o2_write_restart_luc

  implicit none
  
  !!cln  character(len=99), intent(in)  :: fcnpspin
  real,    intent(in)    :: dels
  integer, intent(in)    :: kstart
  integer, intent(in)    :: kend
  integer, intent(in)    :: mloop
  integer, intent(in)    :: lalloc
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
  type(c13o2_luc),           intent(inout) :: c13o2luc
  ! communicator for error-messages
  integer, intent(in)  :: icomm, ocomm

  ! local variables
  integer                  :: myearspin,nyear, nloop1
  character(len=99)        :: ncfile
  character(len=4)         :: cyear
  integer                  :: ktau,ktauday,nday,idoy,ktauy,nloop

  ktauday = nint(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday
  ktau    = 0

  myearspin = cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1
  ! compute the mean fluxes and residence time of each carbon pool
  
  do nyear=1, myearspin
     write(cyear,fmt="(I4)") cable_user%casa_spin_startyear + nyear - 1
     ncfile = trim(casafile%c2cdumppath)//'c2c_'//cyear//'_dump.nc'
     print*, 'dumpfile', ncfile
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
        ! casaflux%Nmindep and casaflux%Pdep set in read_casa_dump
        ! 13C
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
        endif

        call master_send_input(icomm, casa_dump_ts, idoy)
      enddo

  enddo

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
           ! casaflux%Nmindep and casaflux%Pdep set in read_casa_dump
           ! 13C
           if (cable_user%c13o2) then
              c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
              c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
           endif

           call master_send_input(icomm, casa_dump_ts, idoy)
        enddo ! end doy

     enddo   ! end of nyear

  enddo     ! end of nloop
  
  ! write(*,*) 'b4 master receive casa'
  call master_receive(ocomm, 0, casa_ts)
  ! write(*,*) 'after master receive casa'
  call casa_poolout(ktau, veg, soil, casabiome, &
       casapool, casaflux, casamet, casabal, phen)

  call write_casa_restart_nc(casamet, casapool, casaflux, phen, .true.)
  ! 13C
  if (cable_user%c13o2) then
     call master_receive(ocomm, 0, c13o2_pool_ts)
     call c13o2_write_restart_pools(c13o2pools)
  endif

  if ( cable_user%call_POP .and. (POP%np.gt.0) ) then
     ! write(*,*) 'b4 master receive pop'
     call master_receive_pop(POP, ocomm)
     ! write(*,*) 'after master receive pop'
     call POP_io(pop, casamet, myearspin, 'WRITE_INI', .true.)
  endif

END SUBROUTINE master_spincasacnp

!*********************************************************************************************

SUBROUTINE master_CASAONLY_LUC(dels, kstart, kend, veg, soil, casabiome, casapool, &
     casaflux, casamet, casabal, phen, POP, climate, LALLOC, LUC_EXPT, POPLUC, &
     ! 13C
     c13o2flux, c13o2pools, c13o2luc, &
     icomm, ocomm)

  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module,  ONLY: cable_user, is_casa_time
  USE cable_IO_vars_module, ONLY: landpt, output
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE POP_Types,            only: POP_TYPE
  USE POPMODULE,            ONLY: POPStep, POP_init_single
  USE TypeDef,              ONLY: dp
  USE CABLE_LUC_EXPT,       ONLY: LUC_EXPT_TYPE, read_LUH2, &
       ptos, ptog, stog, gtos, pharv, smharv, syharv, &
       ptoc, ptoq, stoc, stoq, ctos, qtos
  USE POPLUC_Types
  USE POPLUC_Module,        ONLY: POPLUCStep, POPLUC_weights_Transfer, WRITE_LUC_OUTPUT_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, READ_LUC_RESTART_NC, &
       POPLUC_set_patchfrac,  WRITE_LUC_OUTPUT_GRID_NC
  ! 13C
  use cable_c13o2_def,      only: c13o2_flux, c13o2_pool, c13o2_luc
  use cable_c13o2,          only: c13o2_save_luc, c13o2_update_luc, &
       c13o2_write_restart_pools, c13o2_write_restart_luc

  IMPLICIT NONE
  
  real,                      intent(in)    :: dels
  integer,                   intent(in)    :: kstart
  integer,                   intent(in)    :: kend
  integer,                   intent(in)    :: lalloc
  type(veg_parameter_type),  intent(inout) :: veg  ! vegetation parameters
  type(soil_parameter_type), intent(inout) :: soil ! soil parameters
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
  integer           :: myearspin,nyear, yyyy, nyear_dump
  character(len=99) :: ncfile
  character(len=4)  :: cyear
  integer           :: ktau,ktauday,nday,idoy

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  integer :: k, j, l
  integer :: rank, off, cnt, ierr
  ! 13C
  real(dp), dimension(c13o2pools%ntile,c13o2pools%npools) :: casasave
  real(dp), dimension(c13o2luc%nland,c13o2luc%npools)     :: lucsave

  
  !  if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
  !  if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
  !  if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
  !  if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
  !  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))
  !
  !  IF (cable_user%CALL_POP) THEN
  !     Iw = POP%Iwood
  !  ENDIF

  ktauday = int(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday
  !  ctime = 0
  !  CALL zero_sum_casa(sum_casapool, sum_casaflux)
  !       count_sum_casa = 0

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
        ktau = (idoy-1)*ktauday +ktauday

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
        ! 13C
        if (cable_user%c13o2) then
           c13o2flux%cAn12(:) = casamet%cAn12spin(:,idoy)
           c13o2flux%cAn(:)   = casamet%cAn13spin(:,idoy)
        endif
        ! print*, 'MASTER Send CO01'
        call master_send_input(icomm, casa_dump_ts, idoy)

        ! ! zero annual sums
        ! if (idoy==1) CALL casa_cnpflux(casaflux,casapool,casabal,.TRUE.)
        
        ! CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
        !      casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
        !      xksoil,xkleaf,xkleafcold,xkleafdry,&
        !      cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
        !      nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
        !      pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
        !
        ! ! update time-aggregates of casa pools and fluxes
        ! CALL update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
        !                     & .TRUE. , .FALSE., 1)
        ! count_sum_casa = count_sum_casa + 1
        
        
        
        !    ! accumulate annual variables for use in POP
        !    IF(idoy==1 ) THEN
        !       casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * &
        !          0.7 ! (assumes 70% of wood NPP is allocated above ground)
        !       LAImax = casamet%glai
        !       Cleafmean = casapool%cplant(:,1)/real(mdyear)/1000.
        !       Crootmean = casapool%cplant(:,3)/real(mdyear)/1000.
        !    ELSE
        !       casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
        !       LAImax = max(casamet%glai, LAImax)
        !       Cleafmean = Cleafmean + casapool%cplant(:,1)/real(mdyear)/1000.
        !       Crootmean = Crootmean +casapool%cplant(:,3)/real(mdyear)/1000.
        !    ENDIF

        IF (idoy==mdyear) THEN ! end of year

           ! get casa update from workers
           ! print*, 'MASTER Receive CO02'
           call master_receive(ocomm, 0, casa_ts)
           CALL master_receive(ocomm, 0, casa_LUC_ts)

           ! 13C
           if (cable_user%c13o2) then
              ! print*, 'MASTER Receive CO02.1'
              call master_receive(ocomm, 0, c13o2_flux_ts)
              call master_receive(ocomm, 0, c13o2_pool_ts)
              call master_receive(ocomm, 0, c13o2_luc_ts)
           endif

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
           ENDDO

           ! set landuse index for secondary forest POP landscapes
           DO k=1, POP%np
              IF (yyyy.eq.LUC_EXPT%YearStart) THEN
                 if (veg%iLU(POP%Iwood(k)).eq.2) then
                    !POP%LU(k) = 2
                    POP%pop_grid(k)%LU = 2
                 endif
              endif
           ENDDO

           ! zero secondary forest tiles in POP where secondary forest area is zero
           DO k=1,mland
              if ((POPLUC%primf(k)-POPLUC%frac_forest(k))==0.0_dp &
                   .and. (.not.LUC_EXPT%prim_only(k))) then

                 j = landpt(k)%cstart+1
                 do l=1,size(POP%Iwood)
                    if( POP%Iwood(l) == j) then
                       CALL POP_init_single(POP,veg%disturbance_interval,l)
                       exit
                    endif
                 enddo

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
                 endif
              endif
           ENDDO

           CALL POPLUCStep(POPLUC,yyyy)

           CALL POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

           ! transfer POP updates to workers
           ! print*, 'MASTER Send CO03'
           off = 1
           DO rank=1, wnp
              IF (rank .GT. 1) off = off + wland(rank-1)%npop_iwood
              cnt = wland(rank)%npop_iwood
              CALL MPI_Send(POP%pop_grid(off), cnt, pop_ts, rank, 0, icomm, ierr)
           END DO

           ! workers call POP here

           ! CALL POPdriver(casaflux,casabal,veg, POP)
           ! print*, 'MASTER Receive CO04'
           CALL master_receive_pop(POP, ocomm)

           ! 13C
           if (cable_user%c13o2) call c13o2_save_luc(casapool, popluc, casasave, lucsave)
           CALL POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
           ! 13C
           if (cable_user%c13o2) call c13o2_update_luc(casasave, lucsave, popluc, luc_expt%prim_only, c13o2pools, c13o2luc)

           IF (output%grid(1:3) == 'lan') THEN
              CALL WRITE_LUC_OUTPUT_NC(POPLUC, YYYY, (YYYY.EQ.cable_user%YearEnd))
           ELSE
              CALL WRITE_LUC_OUTPUT_GRID_NC(POPLUC, YYYY, (YYYY.EQ.cable_user%YearEnd))
           ENDIF

           CALL POPLUC_set_patchfrac(POPLUC,LUC_EXPT)

        ENDIF  ! end of year

     enddo
     ! send updates for CASA pools, resulting from LUC
     ! print*, 'MASTER Send CO05'
     CALL master_send_input(icomm, casa_LUC_ts, nyear)
     ! 13C
     if (cable_user%c13o2) then
        ! print*, 'MASTER Send CO06'
        call master_send_input(icomm, c13o2_flux_ts, nyear)
        call master_send_input(icomm, c13o2_pool_ts, nyear)
        call master_send_input(icomm, c13o2_luc_ts, nyear)
     endif
  enddo ! year=1,myearspin
  
  CALL WRITE_LUC_RESTART_NC( POPLUC, YYYY )
  CALL write_casa_restart_nc(casamet, casapool, casaflux, phen, .TRUE.)
  ! 13C
  if (cable_user%c13o2) then
     call c13o2_write_restart_pools(c13o2pools)
     if (cable_user%POPLUC) call c13o2_write_restart_luc(c13o2luc)
  endif

  CALL POP_IO(pop, casamet, myearspin, 'WRITE_INI', .TRUE.)

END SUBROUTINE master_CASAONLY_LUC


!********************************************************************************************************
! subroutine for reading LU input data, zeroing biomass in empty secondary forest tiles
! and tranferring LUC-based age weights for secondary forest to POP structure


SUBROUTINE LUCdriver(casabiome,casapool, casaflux, POP, LUC_EXPT, POPLUC, veg, &
     ! 13C
     c13o2pools)

  USE cable_def_types_mod , ONLY: r_2, veg_parameter_type, mland
  USE cable_carbon_module
  USE cable_common_module,  ONLY: cable_user, is_casa_time, CurYear
  USE cable_IO_vars_module, ONLY: landpt
  USE casadimension
  USE casaparm
  USE casavariable
  USE POP_Types,            ONLY: POP_TYPE
  USE POPMODULE,            ONLY: POPStep, POP_init_single
  USE CABLE_LUC_EXPT,       ONLY: LUC_EXPT_TYPE, read_LUH2, &
       ptos, ptog, stog, gtos, pharv, smharv, syharv, &
       ptoc, ptoq, stoc, stoq, ctos, qtos
  USE POPLUC_Types
  USE POPLUC_Module,        ONLY: POPLUCStep, POPLUC_weights_Transfer, WRITE_LUC_OUTPUT_NC, &
       POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, READ_LUC_RESTART_NC
  ! 13C
  use cable_c13o2_def,      only: c13o2_pool

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
     POPLUC%gtos(k)   = real(LUC_EXPT%INPUT(gtos)%VAL(k), r_2)
     POPLUC%pharv(k)  = real(LUC_EXPT%INPUT(pharv)%VAL(k), r_2)
     POPLUC%smharv(k) = real(LUC_EXPT%INPUT(smharv)%VAL(k), r_2)
     POPLUC%syharv(k) = real(LUC_EXPT%INPUT(syharv)%VAL(k), r_2)

     POPLUC%ptoc(k) = real(LUC_EXPT%INPUT(ptoc)%VAL(k), r_2)
     POPLUC%ptoq(k) = real(LUC_EXPT%INPUT(ptoq)%VAL(k), r_2)
     POPLUC%stoc(k) = real(LUC_EXPT%INPUT(stoc)%VAL(k), r_2)
     POPLUC%stoq(k) = real(LUC_EXPT%INPUT(stoq)%VAL(k), r_2)
     POPLUC%ctos(k) = real(LUC_EXPT%INPUT(ctos)%VAL(k), r_2)
     POPLUC%qtos(k) = real(LUC_EXPT%INPUT(qtos)%VAL(k), r_2)

     POPLUC%thisyear = yyyy
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

!*********************************************************************************************************

END MODULE cable_mpimaster
