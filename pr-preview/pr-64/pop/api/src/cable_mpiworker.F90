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
! Purpose: Offline driver for mpi worker in CABLE global run
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
!              casa_feedback
!              cbm
!              bgcdriver
!              sumcflux
!              find_extents
!              worker_decomp
!              worker_cable_params
!              worker_casa_params
!              worker_intype
!              worker_outtype
!              worker_casa_type
!              worker_restart_type
!              worker_end
!              ! 13C
!              worker_c13o2_flux_params
!              worker_c13o2_pool_params
!              worker_c13o2_luc_params
!              worker_c13o2_flux_types
!              worker_c13o2_pool_types
!              worker_c13o2_luc_types
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
MODULE cable_mpiworker

  use cable_mpicommon
  use cable_IO_vars_module, only: wlogn

  implicit none

  save

  private

  ! MPI: MPI derived datatype for receiving input from the master
  integer :: inp_t

  ! MPI: MPI derived datatype for sending results back to the master
  integer :: send_t

  ! worker's struct for sending final casa results to the master
  integer :: casa_t

  ! worker's struct for rec'ing/sending final casa results to/from the master
  integer :: casa_dump_t

  ! worker's struct for rec'ing/sending casa pools to/from the master (for LUC calcs)
  integer :: casa_LUC_t

  ! worker's struct for rec'ing/sending final casa results to/from the master
  integer :: climate_t

  ! worker's struct for rec'ing/sending pop io to/from the master
  integer :: pop_t

  ! worker's struct for rec'ing/sending blaze restart to/from the master
  integer :: blaze_in_t
  integer :: blaze_out_t
  integer :: blaze_restart_t

  ! 13C
  ! MPI derived datatype handles for receiving c13o2 results from the workers and restart values
  integer :: c13o2_flux_t
  integer :: c13o2_pool_t
  integer :: c13o2_luc_t

  ! worker's struct for restart data to the master
  integer :: restart_t

  ! worker's logfile unit
  !INTEGER :: wlogn

  public :: mpidrv_worker

contains

  subroutine mpidrv_worker(comm)

    use mpi

    use cable_def_types_mod
    use cable_io_vars_module, only: logn, gswpfile, ncciy, leaps, &
         verbose, fixedCO2, output, check, patchout, soilparmnew, &
         landpt, latitude, longitude
    use cable_common_module,  only: ktau_gl, kend_gl, knode_gl, cable_user, &
         cable_runtime, filename, &
         redistrb, wiltParam, satuParam, CurYear, &
         IS_LEAPYEAR, IS_CASA_TIME, calcsoilalbedo, &
         kwidth_gl
    use cable_data_module,    only: driver_type, point2constants
    use cable_cbm_module,     only: cbm
    use cable_climate_mod

    ! modules related to CASA-CNP
    use casadimension,        only: icycle
    use casavariable,         only: casafile, casa_biome, casa_pool, casa_flux, &
         casa_met, casa_balance, print_casa_var
    use phenvariable,         only: phen_variable
    use casa_cable,           only: bgcdriver, POPdriver, casa_feedback, sumcflux
    use casa_inout,           only: casa_cnpflux

    !CLN added
    ! modules related to POP
    use POP_Types,            only: POP_TYPE
    use POP_Constants,        only: rshootfrac

    ! modules related to fire
    use blaze_drv,            only: blaze_driver
    use BLAZE_MOD,            only: TYPE_BLAZE, BLAZE_ACCOUNTING, INI_BLAZE
    use BLAZE_MPI,            only: WORKER_BLAZE_TYPES ! , WORKER_SIMFIRE_TYPES
    use SIMFIRE_MOD,          only: TYPE_SIMFIRE, INI_SIMFIRE

    ! gm
    use cable_adjust_JV_gm_module, only: read_gm_LUT, LUT_VcmaxJmax, LUT_gm, LUT_Vcmax, LUT_Rd

    ! 13C
    use cable_c13o2_def,         only: c13o2_flux, c13o2_pool, c13o2_luc ! , &
    ! c13o2_update_sum_pools, c13o2_zero_sum_pools
    use cable_c13o2,             only: c13o2_sanity_pools
    use mo_isotope,              only: isoratio ! vpdbc13
    use mo_c13o2_photosynthesis, only: c13o2_discrimination_simple, c13o2_discrimination

    implicit none

    ! MPI:
    integer :: comm ! MPI communicator for comms with the workers

    ! CABLE namelist: model configuration, runtime/user switches
    character(len=200), parameter :: cable_namelist='cable.nml'

    ! timing variables
    integer, parameter ::  kstart = 1  ! start of simulation
    integer, parameter ::  mloop  = 30 ! CASA-CNP PreSpinup loops

    integer :: &
         ktau, &  ! increment equates to timestep, resets if spinning up
         ktau_tot, &  ! NO reset when spinning up, total timesteps by model
         kend, &  ! no. of time steps in run
         !CLN      kstart = 1, &  ! timestep to start at
         koffset = 0, &  ! timestep to start at
         ktauday, &  ! day counter for CASA-CNP
         idoy, &  ! day of year (1:365) counter for CASA-CNP
         nyear, &  ! year counter for CASA-CNP
         YYYY, &  !
         LOY, &  ! Length of Year
         rank            ! Rank of this worker

    real      :: dels    ! time step size in seconds
    character :: cRank*4 ! for worker-logfiles

    ! CABLE variables
    type(met_type)       :: met     ! met input variables
    type(air_type)       :: air     ! air property variables
    type(canopy_type)    :: canopy  ! vegetation variables
    type(radiation_type) :: rad     ! radiation variables
    type(roughness_type) :: rough   ! roughness varibles
    type(balances_type)  :: bal     ! energy and water balance variables
    type(soil_snow_type) :: ssnow   ! soil and snow variables
    type(climate_type)   :: climate     ! climate variables

    ! CABLE parameters
    type(soil_parameter_type) :: soil     ! soil parameters
    type(veg_parameter_type)  :: veg      ! vegetation parameters
    type(driver_type)         :: C        ! constants used locally
    type(sum_flux_type)       :: sum_flux ! cumulative flux variables
    type(bgc_pool_type)       :: bgc      ! carbon pool variables

    ! CASA-CNP variables
    type(casa_biome)     :: casabiome
    type(casa_pool)      :: casapool
    type(casa_flux)      :: casaflux
    type(casa_met)       :: casamet
    type(casa_balance)   :: casabal
    type(phen_variable)  :: phen
    type(POP_TYPE)       :: POP

    ! BLAZE variables
    type(TYPE_BLAZE)     :: BLAZE
    type(TYPE_SIMFIRE)   :: SIMFIRE
    ! REAL(r_2), DIMENSION(:), ALLOCATABLE :: POP_TO, POP_CWD, POP_STR

    ! 13C
    type(c13o2_flux) :: c13o2flux
    type(c13o2_pool) :: c13o2pools ! , sum_c13o2pools
    type(c13o2_luc)  :: c13o2luc
    ! I/O
    ! discrimination
    ! integer :: ileaf
    real(r_2), dimension(:,:), allocatable :: gpp
    real(r_2), dimension(:),   allocatable :: Ra

    ! declare vars for switches (default .FALSE.) etc declared thru namelist
    logical, save :: &
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
         CALL1         = .true.

    logical :: casa_time
    logical :: liseod, liseoy ! is end of day, is end of year

    real :: &
         delsoilM, & ! allowed variation in soil moisture for spin up
         delsoilT    ! allowed variation in soil temperature for spin up

    ! MPI:
    logical :: loop_exit     ! MPI: exit flag for bcast to workers
    integer :: stat(MPI_STATUS_SIZE)
    integer :: icomm ! separate dupes of MPI communicator for send and recv
    integer :: ocomm ! separate dupes of MPI communicator for send and recv
    integer :: ierr
    real    :: etime, etimelast

    ! switches etc defined thru namelist (by default cable.nml)
    namelist /cablenml/ &
         filename, & ! TYPE, containing input filenames
         vegparmnew, & ! use new soil param. method
         soilparmnew, & ! use new soil param. method
         calcsoilalbedo, & ! switch: soil colour albedo - Ticket #27
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

    integer :: LALLOC

    ! command line arguments
    integer :: narg, len1, len2
    character(len=500) :: arg1
    character(len=200) :: arg2

    ! END header

    ! Maciej: make sure the variable does not go out of scope
    mp = 0

    ! Open, read and close the namelist file.
    open(10, file=cable_namelist)
    read(10, nml=cablenml)   !where nml=cable defined above
    close(10)

    narg = command_argument_count()
    if (narg > 0) then
       call get_command_argument(1, arg1, len1)
       filename%met = arg1(1:len1)
       call get_command_argument(2, arg2, len2)
       casafile%cnpipool = arg2(1:len2)
    end if

    if (CABLE_USER%POPLUC .and. (trim(CABLE_USER%POPLUC_RunType) == 'static')) &
         CABLE_USER%POPLUC= .false.

    ! Get worker's rank and determine logfile-unit

    ! MPI: TODO: find a way to preserve workers log messages somewhere
    ! (either separate files or collated by the master to a single file
    ! or perhaps use MPI-IO - but probably not gonna work with random length
    ! text strings)
    ! LN: Done!
    if (CABLE_USER%LogWorker) then
       call MPI_Comm_rank(comm, rank, ierr)
       write(cRank, fmt='(I4.4)') rank
       wlogn = 1000+rank
       open(wlogn, file="cable_log_" // cRank, status="REPLACE")
    else
       wlogn = 1000
       open(wlogn, file="/dev/null")
    end if

    ! INITIALISATION depending on nml settings

    CurYear = CABLE_USER%YearStart

    if (icycle >= 11) then
       icycle                     = icycle - 10
       CASAONLY                   = .true.
       CABLE_USER%CASA_DUMP_READ  = .true.
       CABLE_USER%CASA_DUMP_WRITE = .false.
    else if (icycle == 0) then
       CABLE_USER%CASA_DUMP_READ  = .false.
       spincasa                   = .false.
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
       ! for use with POP: makes use of pipe model to partition between stem and leaf
       LALLOC = 3
    else
       LALLOC = 0 ! default
    end if

    if (trim(cable_user%MetType) == 'gpgs') then
       cable_user%MetType = 'gswp'
    end if

    cable_runtime%offline = .true.

    ! associate pointers used locally with global definitions
    call point2constants(C)

    if (l_casacnp  .and. ((icycle == 0) .or. (icycle > 3))) then
       write(*,*) 'icycle must be 1 to 3 when using casaCNP'
       call MPI_Abort(comm, 50, ierr)
    end if
    if ((l_laiFeedbk .or. l_vcmaxFeedbk) .and. (.not. l_casacnp)) then
       write(*,*) 'casaCNP required to get prognostic LAI or Vcmax'
       call MPI_Abort(comm, 51, ierr)
    end if
    if (l_vcmaxFeedbk .and. (icycle < 2)) then
       write(*,*) 'icycle must be 2 to 3 to get prognostic Vcmax'
       call MPI_Abort(comm, 52, ierr)
    end if
    if ((icycle > 0) .and. (.not. soilparmnew)) then
       write(*,*) 'casaCNP must use new soil parameters'
       call MPI_Abort(comm, 53, ierr)
    end if

    ! gm lookup table
    if (cable_user%explicit_gm .and. (len_trim(cable_user%gm_LUT_file) > 1)) then
      write(*,*) 'Reading gm LUT file'
      call read_gm_LUT(cable_user%gm_LUT_file, LUT_VcmaxJmax, LUT_gm, &
           LUT_Vcmax, LUT_Rd)
    end if

    ! MPI: master only; necessary info will be received by MPI below

    ! Check for leap-year settings
    call MPI_Bcast(leaps, 1, MPI_LOGICAL, 0, comm, ierr)

    ktau_tot = 0
    SPINLOOP: do
       YEARLOOP: do YYYY=CABLE_USER%YearStart, CABLE_USER%YearEnd
          CurYear = YYYY
          if (leaps .and. IS_LEAPYEAR(YYYY)) then
             LOY = 366
          else
             LOY = 365
          end if

          if (CALL1) then
             if (.not. spinup) spinConv=.true.
             ! MPI: bcast to workers so that they don't need to open the met file themselves
             call MPI_Bcast(dels, 1, MPI_REAL, 0, comm, ierr)
          end if

          ! MPI: receive from master ending time fields
          call MPI_Bcast(kend, 1, MPI_INTEGER, 0, comm, ierr)

          if (CALL1) then
             ! MPI: need to know extents before creating datatypes
             call find_extents()

             ! MPI: receive decomposition info from the master
             call worker_decomp(comm)

             ! MPI: in overlap version sends and receives occur on separate comms
             call MPI_Comm_dup(comm, icomm, ierr)
             call MPI_Comm_dup(comm, ocomm, ierr)

             ! MPI: data set in load_parameter is now received from
             ! the master
             call worker_cable_params(comm, met, air, ssnow, veg, bgc, soil, &
                  canopy, rough, rad, sum_flux, bal)
             ! 13C
             if (cable_user%c13o2) then
                allocate(gpp(size(canopy%An, 1), size(canopy%An, 2)))
                allocate(Ra(size(canopy%An, 1)))
             end if

             ktauday = int(24.0 * 3600.0 / dels)
             if (cable_user%call_climate) &
                  call worker_climate_types(comm, climate, ktauday)

             ! MPI: mvtype and mstype send out here instead of inside worker_casa_params
             !      so that old CABLE carbon module can use them. (BP May 2013)
             call MPI_Bcast(mvtype, 1, MPI_INTEGER, 0, comm, ierr)
             call MPI_Bcast(mstype, 1, MPI_INTEGER, 0, comm, ierr)

             ! MPI: casa parameters received only if cnp module is active
             if (icycle > 0) then
                call worker_casa_params(comm, casabiome, casapool, casaflux, casamet, casabal, phen)
                ! 13C
                if (cable_user%c13o2) then
                   call worker_c13o2_flux_params(comm, c13o2flux)
                   call worker_c13o2_pool_params(comm, c13o2pools)
                   if (cable_user%popluc) then
                      call worker_c13o2_luc_params(comm, c13o2luc)
                   end if
                end if
                ! MPI: POP restart received only if pop module AND casa are active
                if (cable_user%call_pop) then
                   call worker_pop_types(comm, veg, pop)
                end if

                ! CLN:
                if (cable_user%call_blaze) then
                   call ini_blaze(mland, rad%latitude(landpt(:)%cstart), &
                        rad%longitude(landpt(:)%cstart), blaze)

                   !par blaze restart not required uses climate data
                   allocate(latitude(mland))
                   allocate(longitude(mland))
                   call worker_blaze_types(comm, mland, blaze, blaze_restart_t, &
                        blaze_in_t, blaze_out_t)
                   ! cln:  burnt_area
                   if ( trim(blaze%burnt_area_src) == "SIMFIRE" ) then
                      call MPI_recv(MPI_BOTTOM, 1, blaze_in_t, 0, ktau_gl, &
                           comm, stat, ierr)
                      !CLN here we need to check for the SIMFIRE biome setting
                      call INI_SIMFIRE(mland, SIMFIRE, climate%modis_igbp(landpt(:)%cstart))
                      !MC temporary
                      write(wlogn,*) "After ini_simf"
                      call flush(wlogn)
                      !par blaze restart not required uses climate data
                   end if
                end if
             end if ! icycle > 0

             ! MPI: create inp_t type to receive input data from the master
             ! at the start of every timestep
             call worker_intype(comm, met, veg)

             ! MPI: casa parameters received only if cnp module is active
             ! MPI: create send_t type to send the results to the master
             ! at the end of every timestep
             call worker_outtype(comm, met, canopy, ssnow, rad, bal, air, soil, veg)
             if (cable_user%c13o2) call worker_c13o2_flux_type(comm, c13o2flux)

             ! MPI: casa parameters received only if cnp module is active
             ! MPI: create type to send casa results back to the master
             ! only if cnp module is active
             if (icycle > 0) then
                call worker_casa_type(comm, casabiome, casapool, casaflux, &
                     casamet, casabal, phen)
                ! 13C
                if (cable_user%c13o2) call worker_c13o2_pool_type(comm, c13o2pools)

                if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                   call worker_casa_dump_types(comm, casamet, casaflux, phen, &
                        climate, c13o2flux)
                end if

                if (cable_user%popluc) then
                   call worker_casa_LUC_types(comm, casapool, casabal, casaflux)
                   ! MPI: casa parameters received only if cnp module is active
                   ! 13C
                   if (cable_user%c13o2) &
                        call worker_c13o2_luc_type(comm, c13o2luc)
                end if
             end if

             ! MPI: create type to send restart data back to the master
             ! only if restart file is to be created
             if (output%restart) then
                call worker_restart_type(comm, canopy, air, veg, ssnow)
             end if

             canopy%fes_cor = 0.0_r_2
             canopy%fhs_cor = 0.
             met%ofsd       = 0.1 ! not used
             !MC diff to serial code - pdep not set
             ! vh !
             met%pdep = 0.0
             !MC - do we need that? casamet also set only if icycle>0
             ! if (trim(cable_user%MetType) .eq. 'cru') then
             !    casamet%glai = 1.0_r_2
             !    where (veg%iveg(:) .ge. 14) casamet%glai = 0.0_r_2
             ! end if
             if ((trim(cable_user%MetType) == 'cru') .and. (icycle > 0)) then
                if (all(real(casamet%glai(:)) == 0.)) then
                   write(*,*) 'W A R N I N G : deleted setting casamet%glai for CRU.'
                   call MPI_Abort(comm, 70, ierr)
                end if
             end if

             ! CALL worker_sumcasa_types(comm, sum_casapool, sum_casaflux)
             ! ! 13C
             ! if (cable_user%c13o2) call c13o2_zero_sum_pools(sum_c13o2pools)
             ! count_sum_casa = 0

             if ((icycle > 0) .and. spincasa) then
                write(wlogn,*) 'EXT spincasacnp enabled with mloop=', mloop
                call worker_spincasacnp(dels, kstart, kend, mloop, veg, soil, &
                     casabiome, casapool, casaflux, casamet, casabal, phen, &
                     POP, climate, LALLOC, c13o2flux, c13o2pools, icomm, ocomm)
                if (cable_user%c13o2) &
                     call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                SPINconv = .false.
                CASAONLY = .true.
                ktau_gl  = 0
                ktau     = 0
             else if (casaonly .and. (.not. spincasa) .and. cable_user%popluc) then
                call worker_CASAONLY_LUC(dels, kstart, kend, veg, soil, &
                     casabiome, casapool, casaflux, casamet, casabal, phen, &
                     POP, climate, LALLOC, c13o2flux, c13o2pools, icomm, ocomm)
                if (cable_user%c13o2) &
                     call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                SPINconv = .false.
                ktau_gl  = 0
                ktau     = 0
             end if

          ! else ! CALL1

          !    if (icycle > 0) then
          !       ! re-initalise annual flux sums
          !       casabal%FCgppyear = 0.0
          !       casabal%FCrpyear  = 0.0
          !       casabal%FCnppyear = 0.0
          !       casabal%FCrsyear  = 0.0
          !       casabal%FCneeyear = 0.0
          !    end if

          end if ! CALL1

          ! globally (WRT code) accessible kend through USE cable_common_module
          ktau_gl   = 0
          kwidth_gl = int(dels)
          kend_gl   = kend
          knode_gl  = 0

          if (spincasa .or. casaonly) then
             exit
          end if

          ! IF (.NOT.spincasa) THEN
          ! time step loop over ktau
          KTAULOOP: do ktau=kstart, kend
             call CPU_time(etimelast)
             ! increment total timstep counter
             ktau_tot = ktau_tot + 1

             ! globally (WRT code) accessible kend through USE cable_common_module
             ktau_gl = ktau_gl + 1

             ! Check for today's day-of-year
             idoy = int( mod(real(ceiling((real(ktau+koffset)) / real(ktauday))), &
                  real(loy)) )
             if (idoy == 0) idoy = LOY

             ! needed for CASA-CNP
             !MC - diff to serial
             nyear = int((kend-kstart+1) / (LOY*ktauday))

             ! MPI: receive input data for this step from the master
             if (.not. CASAONLY) then
                call MPI_Recv(MPI_BOTTOM, 1, inp_t, 0, ktau_gl, icomm, stat, ierr)
                ! MPI: receive casa_dump_data for this step from the master
             else if (IS_CASA_TIME("dread", yyyy, ktau, kstart, koffset, &
                  ktauday, wlogn)) then
                if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                   call MPI_Recv(MPI_BOTTOM, 1, casa_dump_t, 0, ktau_gl, icomm, stat, ierr)
                end if
             end if

             ! Get met data and LAI, set time variables.
             ! Rainfall input may be augmented for spinup purposes:
             met%ofsd = met%fsd(:,1) + met%fsd(:,2)
             canopy%oldcansto = canopy%cansto
             ! Zero out lai where there is no vegetation acc. to veg. index
             where (veg%iveg(:) .ge. 14) veg%vlai = 0.

             ! MPI: some fields need explicit init, because we don't transfer
             ! them for better performance
             ! in the serial version this is done in get_met_data
             ! after input has been read from the file
             met%tvair = met%tk
             met%tvrad = met%tk

             ! 13C
             if (cable_user%c13o2) then
                call MPI_Recv(c13o2flux%ca(1), mp, MPI_DOUBLE_PRECISION, 0, 0, icomm, stat, ierr)
             end if

             if ((ktau==1) .and. (icycle > 1)) &
                  call casa_cnpflux(casaflux, casapool, casabal, .true.)

             ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
             casa_time = IS_CASA_TIME("write", yyyy, ktau, kstart, koffset, &
                  ktauday, logn)
             liseod = mod((ktau-kstart+1), ktauday) == 0
             liseoy = mod((ktau-kstart+1) / ktauday, LOY) == 0

             if (l_vcmaxFeedbk) then
                if (mod(ktau, ktauday) == 1) then
                   call casa_feedback(ktau, veg, casabiome, &
                        casapool, casamet, climate, ktauday)
                end if
             else
                veg%vcmax_shade = veg%vcmax
                veg%ejmax_shade = veg%ejmax
                veg%vcmax_sun = veg%vcmax
                veg%ejmax_sun = veg%ejmax
             end if

             !MC - if (l_laiFeedbk) veg%vlai(:) = real(casamet%glai(:))
             if (l_laiFeedbk .and. (icycle > 0)) &
                  veg%vlai(:) = real(casamet%glai(:))

             ! CALL land surface scheme for this timestep, all grid points:
             CALL cbm(ktau, dels, air, bgc, canopy, met, &
                  bal, rad, rough, soil, ssnow, &
                  veg, climate)

             ! 13C
             if (cable_user%c13o2) then
                gpp  = canopy%An + canopy%Rd
                Ra   = isoratio(c13o2flux%ca, real(met%ca,r_2), 1.0_r_2)
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

             !TRUNK - call of cable_climet before cbm
             if (cable_user%CALL_climate) &
                  CALL cable_climate(ktau, kstart, ktauday, idoy, LOY, met, &
                  climate, canopy, veg, ssnow, rad, dels, mp)

             ssnow%smelt  = ssnow%smelt  * dels
             ssnow%rnof1  = ssnow%rnof1  * dels
             ssnow%rnof2  = ssnow%rnof2  * dels
             ssnow%runoff = ssnow%runoff * dels

             call CPU_time(etime)
             write(wlogn,*) 'ktau, etime-etimelast ', ktau, etime-etimelast

             ! serial: IF ((icycle > 0) .OR. CABLE_USER%CASA_DUMP_WRITE) THEN
             if (icycle > 0) then

                call bgcdriver( ktau, kstart, dels, met, &
                     ssnow, canopy, veg, soil,climate, casabiome, &
                     casapool, casaflux, casamet, casabal, &
                     phen, pop, ktauday, idoy, loy, &
                     .false., LALLOC, c13o2flux, c13o2pools )
                if (cable_user%c13o2) &
                     call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
                write(wlogn,*) 'after bgcdriver ', MPI_BOTTOM, 1, casa_t, 0, ktau_gl, ocomm, ierr

                !TRUNK no if (liseod) then
                if (liseod) then
                   call MPI_Send(MPI_BOTTOM, 1, casa_t, 0, ktau_gl, ocomm, ierr)
                   ! 13C
                   if (cable_user%c13o2) &
                        call MPI_Send(MPI_BOTTOM, 1, c13o2_pool_t, 0, ktau_gl, ocomm, ierr)
                   write(wlogn,*) 'after casa mpi_send ', ktau
                end if

                if (liseod) then
                   if (cable_user%CALL_BLAZE) then
                      call BLAZE_ACCOUNTING(BLAZE, climate, ktau, dels, YYYY, idoy)
                      call blaze_driver(blaze%ncells, blaze, simfire, casapool, casaflux, &
                           casamet, climate, rshootfrac, idoy, YYYY, 1, POP, veg)
                      !par send blaze_out back
                      call MPI_Send(MPI_BOTTOM, 1, blaze_out_t, 0, ktau_gl, ocomm, ierr)
                   end if
                end if

                if (casa_time) then
                   write(wlogn,*) 'IN IS_CASA ', casapool%cplant(:,1)
                end if

                ! MPI: send the results back to the master
                if (((.not.spinup) .or. (spinup .and. spinConv)) .and. &
                     IS_CASA_TIME("dwrit", yyyy, ktau, kstart, &
                     koffset, ktauday, wlogn) ) then
                   if (cable_user%casa_dump_read .or. cable_user%casa_dump_write) then
                      call MPI_Send(MPI_BOTTOM, 1, casa_dump_t, 0, ktau_gl, ocomm, ierr)
                   end if
                end if

             end if ! icycle>0

             if (.not. casaonly) then
                ! sumcflux is pulled out of subroutine cbm
                ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
                call sumcflux( ktau, kstart, dels, &
                     canopy, sum_flux, &
                     casaflux, l_vcmaxFeedbk )
                write(wlogn, *) 'after sumcflux ', ktau
             end if
             ! MPI: send the results back to the master
             CALL MPI_Send(MPI_BOTTOM, 1, send_t, 0, ktau_gl, ocomm, ierr)
             if (cable_user%c13o2) call MPI_Send(MPI_BOTTOM, 1, c13o2_flux_t, 0, ktau_gl, ocomm, ierr)

             CALL1 = .false.
             write(wlogn, *) 'after send results back to master ', ktau

          end do KTAULOOP ! END Do loop over timestep ktau
          ! ELSE

          !    CALL1 = .FALSE.
          ! END IF

          flush(wlogn)

          if (icycle >0 .and. cable_user%CALL_POP) then

             if (CABLE_USER%POPLUC) then
                write(wlogn, *) 'before MPI_Send casa_LUC'
                ! worker sends casa updates required for LUC calculations here
                call MPI_Send(MPI_BOTTOM, 1, casa_LUC_t, 0, 0, ocomm, ierr)
                ! 13C
                if (cable_user%c13o2) then
                   call MPI_Send(MPI_BOTTOM, 1, c13o2_flux_t, 0, 0, ocomm, ierr)
                   call MPI_Send(MPI_BOTTOM, 1, c13o2_pool_t, 0, 0, ocomm, ierr)
                   call MPI_Send(MPI_BOTTOM, 1, c13o2_luc_t, 0, 0, ocomm, ierr)
                end if
                write(wlogn, *) 'after MPI_Send casa_LUC'
                ! master calls LUCDriver here
                ! worker receives casa and POP updates
                call MPI_Recv(POP%pop_grid(1), POP%np, pop_t, 0, 0, icomm, stat, ierr)
             end if
             ! one annual time-step of POP
             call POPdriver(casaflux, casabal, veg, POP)

             ! Call BLAZE again to compute turnovers depending on POP mortalities
             if ( cable_user%CALL_BLAZE ) then
                !MC - this is different to serial code
                call blaze_driver(blaze%ncells, blaze, simfire, casapool, casaflux, &
                     casamet, climate, rshootfrac, idoy, YYYY, 1, POP, veg)
             end if

             call worker_send_pop(POP, ocomm)

             if (CABLE_USER%POPLUC) then
                call MPI_Recv(MPI_BOTTOM, 1, casa_LUC_t, 0, nyear, icomm, stat, ierr)
                ! 13C
                if (cable_user%c13o2) then
                   call MPI_Recv(MPI_BOTTOM, 1, c13o2_flux_t, 0, nyear, icomm, stat, ierr)
                   call MPI_Recv(MPI_BOTTOM, 1, c13o2_pool_t, 0, nyear, icomm, stat, ierr)
                   call MPI_Recv(MPI_BOTTOM, 1, c13o2_luc_t, 0, nyear, icomm, stat, ierr)
                end if
             end if

          end if ! icycle >0 .and. cable_user%CALL_POP

          if ((icycle > 0) .and. (.not. casaonly)) then
             ! re-initalise annual flux sums
             casabal%FCgppyear = 0.0_r_2
             casabal%FCrpyear  = 0.0_r_2
             casabal%FCnppyear = 0.0_r_2
             casabal%FCrsyear  = 0.0_r_2
             casabal%FCneeyear = 0.0_r_2
          end if

       end do YEARLOOP

       if (spincasa .or. casaonly) then
          exit
       end if

       ! MPI: learn from the master whether it's time to quit
       call MPI_Bcast(loop_exit, 1, MPI_LOGICAL, 0, comm, ierr)

       if (loop_exit) then
          exit
       end if

    end do SPINLOOP

    if ((icycle > 0) .and. (.not. spincasa) .and. (.not. casaonly)) then
       ! MPI: send casa results back to the master
       call MPI_Send(MPI_BOTTOM, 1, casa_t, 0, ktau_gl, ocomm, ierr)
       ! 13C
       if (cable_user%c13o2) then
          call MPI_Send(MPI_BOTTOM, 1, c13o2_flux_t, 0, ktau_gl, ocomm, ierr)
          call MPI_Send(MPI_BOTTOM, 1, c13o2_pool_t, 0, ktau_gl, ocomm, ierr)
       end if
       if (cable_user%call_blaze) then
         !par send blaze_out back
         call MPI_Send(MPI_BOTTOM, 1, blaze_out_t, 0, ktau_gl, ocomm, ierr)
       end if
    end if

    ! Write restart file if requested:
    if (output%restart .and. (.not. CASAONLY)) then
       ! MPI: send variables that are required by create_restart
       call MPI_Send(MPI_BOTTOM, 1, restart_t, 0, ktau_gl, comm, ierr)
       ! MPI: output file written by master only
       if (cable_user%CALL_climate) then
          call MPI_Send(MPI_BOTTOM, 1, climate_t, 0, ktau_gl, comm, ierr)
       end if
    end if

    ! MPI: cleanup
    call worker_end(icycle, output%restart)

    ! Close log file
    ! MPI: closes handle to /dev/null in workers
    close(wlogn)

    return

  end subroutine mpidrv_worker



  ! ============== PRIVATE SUBROUTINES USED ONLY BY THE MPI WORKERS ===============


  ! MPI: receives grid decomposition info from the master
  SUBROUTINE worker_decomp(comm)

    use mpi
    USE cable_def_types_mod, ONLY: mland, mp

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: comm ! MPI communicator to talk to the workers

    INTEGER :: stat(MPI_STATUS_SIZE), ierr

    ! receive number of landpoints assigned to this worker
    CALL MPI_Recv(mland, 1, MPI_INTEGER, 0, 0, comm, stat, ierr)

    ! receive number of land patches assigned to this worker
    CALL MPI_Recv(mp, 1, MPI_INTEGER, 0, 0, comm, stat, ierr)

    RETURN

  END SUBROUTINE worker_decomp


  ! MPI: creates param_t type for the worker to receive the default parameters
  ! from the master process
  ! then receives the parameters
  ! and finally frees the MPI type
  SUBROUTINE worker_cable_params(comm, met, air, ssnow, veg, bgc, soil, canopy, &
       rough, rad, sum_flux, bal)

    use mpi

    USE cable_def_types_mod
    USE cable_IO_vars_module
    USE cable_input_module,  ONLY: allocate_cable_vars
    use cable_mpicommon,     only: nparam, add_address_1block

    IMPLICIT NONE

    ! subroutine arguments

    INTEGER, INTENT(IN) :: comm ! MPI communicator
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
    INTEGER :: tsize

    INTEGER :: stat(MPI_STATUS_SIZE), ierr
    INTEGER :: landp_t, patch_t, param_ts

    INTEGER :: r1len, r2len, i1len, llen ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank, ierr2, rcount, pos, offset

    CHARACTER, DIMENSION(:), ALLOCATABLE :: rbuf

    CALL MPI_Comm_rank(comm, rank, ierr)

    ! mp and mland should have been received previously by
    ! worker_decomp

    ! creates types to receive slices of landpt and patch arrays from the master
    CALL decomp_types(landp_t, patch_t)

    ! Allocate spatial heterogeneity variables:
    ALLOCATE(landpt(mland))

    ! and receive own slice from the master
    CALL MPI_Recv(landpt, mland, landp_t, 0, 0, comm, stat, ierr)

    !adjust cstart & cend for worker slice
    offset=landpt(1)%cstart
    landpt(:)%cstart=landpt(:)%cstart-offset+1
    landpt(:)%cend=landpt(:)%cend-offset+1

    CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssnow, &
         sum_flux,veg,mp)

    ! receive slice of patch array that was allocated above inside
    ! allocate_cable_vars
    CALL MPI_Recv(patch, mp, patch_t, 0, 0, comm, stat, ierr)

    ! MPI: TODO: probably not a bad idea to free landp_t and patch_t types
    ntyp = nparam

    ALLOCATE(blen(ntyp))
    ALLOCATE(displs(ntyp))
    ALLOCATE(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid
    llen = mp * extl

    bidx = 0

    ! the order of variables follows argument list
    ! the order of fields within follows alloc_*_type subroutines

    ! ----------- met --------------

    call add_address_1block(met%year, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%moy, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%ca, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%doy, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%hod, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%ofsd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%fld, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%precip, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%precip_sn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%tk, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%tvair, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%tvrad, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%pmb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%ua, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%qv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%qvair, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%da, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%dva, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%coszen, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%Ndep, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%Pdep, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%u10, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%rhum, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%fsd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(met%fdiff, 1, mp, displs, blen, types, bidx)

    ! ----------- air --------------

    call add_address_1block(air%rho, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%volm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%rlam, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%qsat, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%epsi, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%visc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%psyc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%dsatdk, 1, mp, displs, blen, types, bidx)
    call add_address_1block(air%cmolar, 1, mp, displs, blen, types, bidx)

    ! ----------- ssnow --------------

    call add_address_1block(ssnow%isflag, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%iantrct, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%pudsto, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%pudsmx, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%cls, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%dfn_dtg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%dfh_dtg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%dfe_ddq, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%ddq_dtg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%evapsn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%fwtop, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%fwtop1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%fwtop2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%fwtop3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%osnowd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%potev, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%runoff, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%rnof1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%rnof2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%rtsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wbtot1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wbtot2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wb_lake, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%sinfil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%qstss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wetfac, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%owetfac, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%t_snwlr, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tggav, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%otgg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%otss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%otss_0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tprecip, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tevap, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%trnoff, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%totenbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%totenbal2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%fland, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%ifland, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%qasrf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%qfsrf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%qssrf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%snage, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%snowd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%smelt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%ssdnn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tss_p, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%deltss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%owb1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%sconds, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%sdepth, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%smass, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%ssdn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tgg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tggsn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%dtmlt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%albsoilsn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%evapfbl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%tilefrac, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wbtot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%gammzz, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wbice, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wblf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wbfice, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%S, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Tsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%SL, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%TL, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%h0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%rex, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%wflux, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%delwcol, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%zdelta, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%kth, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Tsurface, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%lE, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%evap, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%ciso, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%cisoL, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%rlitt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%thetai, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%snowliq, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%nsteps, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%TsurfaceFR, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Ta_daily, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%nsnow, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Qadv_daily, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%G0_daily, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Qevap_daily, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Qprec_daily, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Qprec_snow_daily, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%E_fusion_sn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%E_sublimation_sn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%latent_heat_sn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%evap_liq_sn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%surface_melt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(ssnow%Qadv_rain_sn, 1, mp, displs, blen, types, bidx)

    ! ----------- veg --------------

    call add_address_1block(veg%iveg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ivegp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%iLU, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%canst1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%dleaf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ejmax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ejmax_shade, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ejmax_sun, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%meth, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%frac4, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%hc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vlai, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%xalbnir, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%rp20, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%rpcoef, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%rs20, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%shelrb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vegcf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%tminvj, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%toptvj, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%tmaxvj, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vbeta, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vcmax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vcmax_shade, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vcmax_sun, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%xfang, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%extkn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vlaimax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%wai, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%a1gs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%d0gs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%alpha, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%convex, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%cfrd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%gswmin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%conkc0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%conko0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ekc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%eko, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%g0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%g1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%vcmaxcc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ejmaxcc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%gmmax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%gm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%c4kci, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%c4kcc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%bjv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%deciduous, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%refl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%taul, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%froot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%rootbeta, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%gamma, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%ZR, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%F10, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%clitt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%disturbance_interval, 1, mp, displs, blen, types, bidx)
    call add_address_1block(veg%disturbance_intensity, 1, mp, displs, blen, types, bidx)

    ! ----------- bgc --------------

    call add_address_1block(bgc%cplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bgc%csoil, 1, mp, displs, blen, types, bidx)
    ! ncp, each worker gets the same copy of whole array
    call add_address_1block(bgc%ratecp, 1, ncp, displs, blen, types, bidx)
    ! ncs, each worker gets the same copy of whole array
    call add_address_1block(bgc%ratecs, 1, ncs, displs, blen, types, bidx)

    ! ----------- soil --------------

    call add_address_1block(soil%isoilm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%bch, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%c3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%clay, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%css, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%hsbh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%hyds, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%i2bp3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%ibp2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%rhosoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%sand, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%sfc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%silt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%ssat, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%sucs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%swilt, 1, mp, displs, blen, types, bidx)
    ! ms, each worker gets the same copy of whole array
    call add_address_1block(soil%zse, 1, ms, displs, blen, types, bidx)
    ! ms+1, each worker gets the same copy of whole array
    call add_address_1block(soil%zshh, 1, ms+1, displs, blen, types, bidx)
    call add_address_1block(soil%soilcol, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%albsoilf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%cnsd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%pwb_min, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%albsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%nhorizons, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%ishorizon, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%clitt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%zeta, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%fsatmax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%swilt_vec, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%ssat_vec, 1, mp, displs, blen, types, bidx)
    call add_address_1block(soil%sfc_vec, 1, mp, displs, blen, types, bidx)

    ! ----------- canopy --------------

    call add_address_1block(canopy%cansto, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%cduv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%delwc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%dewmm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fe, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fpn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%frp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%frpw, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%frpr, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%frs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fnee, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%frday, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fnv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fev, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%epot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fnpp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fevw_pot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gswx_T, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%cdtq, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%wetfac_cs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fevw, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fhvw, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%oldcansto, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fhv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fns, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fhs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fhs_cor, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ga, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ghflux, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%precis, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%qscrn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%rnet, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%rniso, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%segg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%sghflux, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%through, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%spill, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%tscrn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%wcint, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%tv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%us, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%uscrn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%vlaiw, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%rghlai, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fwet, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%evapfbl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gswx, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%zetar, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%zetash, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fess, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fesp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%dgdtg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fes, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fes_cor, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fevc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ofes, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%A_sl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%A_sh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%A_slC, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%A_shC, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%A_slJ, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%A_shJ, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%GPP_sl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%GPP_sh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fevc_sl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fevc_sh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_A_cs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%dAdcs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_GPP_cs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_A_cs_sl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_A_cs_sh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_fevc_cs_sl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_fevc_cs_sh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%eta_fevc_cs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%cs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%cs_sl, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%cs_sh, 1, mp, displs, blen, types, bidx)
    ! call add_address_1block(canopy%ci_sl, 1, mp, displs, blen, types, bidx)
    ! call add_address_1block(canopy%ci_sh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%tlf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%dlf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gw, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ancj, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%tlfy, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ecy, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ecx, 1, mp, displs, blen, types, bidx)
    ! call add_address_1block(canopy%ci, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%fwsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%kthLitt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%DvLitt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%An, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%Rd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%isc3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%vcmax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gammastar, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gsc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gbc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%gac, 1, mp, displs, blen, types, bidx)
    call add_address_1block(canopy%ci, 1, mp, displs, blen, types, bidx)

    ! ------- rough -------

    call add_address_1block(rough%disp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%hruff, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%hruff_grmx, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%rt0us, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%rt1usa, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%rt1usb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%rt1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%za_uv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%za_tq, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%z0m, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%zref_uv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%zref_tq, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%zruffs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%z0soilsn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%z0soil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%coexp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%usuh, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%term2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%term3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%term5, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%term6, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rough%term6a, 1, mp, displs, blen, types, bidx)

    ! --------rad --------

    call add_address_1block(rad%transb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%albedo_T, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%longitude, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%workp1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%workp2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%workp3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%extkb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%extkd2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%extkd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%flws, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%latitude, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%lwabv, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%qssabs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%transd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%trad, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%fvlai, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%rhocdf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%rniso, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%scalex, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%albedo, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%reffdf, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%reffbm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%extkbm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%extkdm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%fbeam, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%cexpkbm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%cexpkdm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%rhocbm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%gradis, 1, mp, displs, blen, types, bidx)
    call add_address_1block(rad%qcan, 1, mp, displs, blen, types, bidx)

    ! ------- sum_flux -----

    call add_address_1block(sum_flux%sumpn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumrp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumrpw, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumrpr, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumrs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumrd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%dsumpn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%dsumrp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%dsumrs, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%dsumrd, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumxrp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(sum_flux%sumxrs, 1, mp, displs, blen, types, bidx)

    ! ------- bal ----

    call add_address_1block(bal%drybal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%ebal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%ebal_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%ebal_cncheck, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%ebal_tot_cncheck, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%ebaltr, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%ebal_tottr, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%evap_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%osnowd0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%precip_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%rnoff_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%wbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%wbal_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%wbtot0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%wetbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%cansto0, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%owbtot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%evapc_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%evaps_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%rnof1_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%rnof2_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%snowdc_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%wbal_tot1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%delwc_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%qasrf_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%qfsrf_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%qssrf_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%Radbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%EbalSoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%Ebalveg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(bal%Radbalsum, 1, mp, displs, blen, types, bidx)

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,' invalid number of param_ts fields',bidx,', fix it (20)!'
       CALL MPI_Abort(comm, 54, ierr)
    END IF

    CALL MPI_Type_create_struct(bidx, blen, displs, types, param_ts, ierr)
    CALL MPI_Type_commit(param_ts, ierr)

    CALL MPI_Type_size(param_ts, tsize, ierr)
    CALL MPI_Type_get_extent(param_ts, tmplb, text, ierr)

    WRITE (*,*) 'worker param_ts blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...

    ! so, now receive all the parameters
    !    CALL MPI_Recv (MPI_BOTTOM, 1, param_ts, 0, 0, comm, stat, ierr)
    !   Maciej: buffered recv + unpac version
    ALLOCATE(rbuf(tsize))
    CALL MPI_Recv(rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    CALL MPI_Get_count(stat, param_ts, rcount, ierr2)

    IF (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       CALL MPI_Unpack(rbuf, tsize, pos, MPI_BOTTOM, rcount, param_ts, &
            comm, ierr)
       IF (ierr /= MPI_SUCCESS) write(*,*) 'cable param unpack error, rank: ', rank, ierr
    ELSE
       write(*,*) 'cable param recv rank err err2 rcount: ', rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)

    ! finally free the MPI type
    CALL MPI_Type_Free(param_ts, ierr)

    ! all CABLE parameters have been received from the master by now
    RETURN

  END SUBROUTINE worker_cable_params


  ! MPI: creates param_t type for the worker to receive the default casa
  ! parameters from the master process
  ! then receives them
  ! and finally frees the MPI type
  SUBROUTINE worker_casa_params(comm, casabiome, casapool, casaflux, &
       casamet, casabal, phen)

    use mpi
    use cable_def_types_mod
    use casadimension, only: mlitter, mplant, msoil, mphase, mdyear, mso
    use casavariable
    use phenvariable
    use cable_mpicommon, only: ncasaparam, add_address_1block

    implicit none

    ! sub arguments
    integer,             intent(in)  :: comm  ! MPI communicator
    type(casa_biome),    intent(inout) :: casabiome
    type(casa_pool),     intent(inout) :: casapool
    type(casa_flux),     intent(inout) :: casaflux
    type(casa_met),      intent(inout) :: casamet
    type(casa_balance),  intent(inout) :: casabal
    type(phen_variable), intent(inout) :: phen

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    integer, allocatable, dimension(:) :: blen
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
    integer :: tsize

    integer :: stat(MPI_STATUS_SIZE), ierr
    integer :: casa_ts

    integer :: r1len, r2len, i1len, llen ! block lengths
    integer :: bidx ! block index
    integer :: ntyp ! total number of blocks

    integer :: rank, off, ierr2, rcount, pos

    character, dimension(:), allocatable :: rbuf

    off = 1

    call MPI_Comm_rank(comm, rank, ierr)

    if (.not. associated(casabiome%ivt2)) then
       write (*,*) 'worker alloc casa and phen var with m patches: ', rank, mp
       CALL alloc_casa_var(casabiome)
       CALL alloc_casa_var(casapool, mp)
       CALL alloc_casa_var(casaflux, mp)
       CALL alloc_casa_var(casamet, mp)
       CALL alloc_casa_var(casabal, mp)
       call alloc_phenvariable(phen, mp)
       call zero_casa_var(casabiome)
       call zero_casa_var(casapool)
       call zero_casa_var(casaflux)
       call zero_casa_var(casamet)
       call zero_casa_var(casabal)
       call zero_phenvariable(phen)
    end if

    ntyp = ncasaparam

    allocate(blen(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid
    llen  = mp * extl

    bidx = 0

    ! ------- casabiome -----

    call add_address_1block(casabiome%ivt2, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkleafcoldmax, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkleafcoldexp, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkleafdrymax, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkleafdryexp, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%glaimax, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%glaimin, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%sla, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ratiofrootleaf, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%kroot, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%krootlen, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%rootdepth, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%kuptake, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%kminN, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%KuplabP, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%kclabrate, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xnpmax, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%q10soil, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkoptlitter, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkoptsoil, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%maxfinelitter, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%maxcwd, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%prodptase, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%costnpup, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkplab, 1, mso, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkpsorb, 1, mso, displs, blen, types, bidx)
    call add_address_1block(casabiome%xkpocc, 1, mso, displs, blen, types, bidx)
    call add_address_1block(casabiome%nintercept, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%nslope, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%plantrate, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%rmplant, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%fracnpptoP, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%fraclignin, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%fraclabile, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ratioNCplantmin, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ratioNCplantmax, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ratioPCplantmin, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ratioPCplantmax, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%fracLigninplant, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ftransNPtoL, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ftransPPtoL, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%litterrate, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%soilrate, 1, mvtype, displs, blen, types, bidx)
    ! added by ln
    call add_address_1block(casabiome%ratioNPplantmin, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%ratioNPplantmax, 1, mvtype, displs, blen, types, bidx)
    ! added by vh
    call add_address_1block(casabiome%la_to_sa, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%vcmax_scalar, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%disturbance_interval, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%DAMM_EnzPool, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%DAMM_KMO2, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%DAMM_KMcp, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%DAMM_Ea, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(casabiome%DAMM_alpha, 1, mvtype, displs, blen, types, bidx)

    ! ------ casapool ----

    call add_address_1block(casapool%Clabile, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dClabiledt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Cplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Nplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Pplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dCplantdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dNplantdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dPplantdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNCplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioPCplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Nsoilmin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Psoillab, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Psoilsorb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Psoilocc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dNsoilmindt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dPsoillabdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dPsoilsorbdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dPsoiloccdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Clitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Nlitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Plitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dClitterdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dNlitterdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dPlitterdt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNClitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioPClitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Csoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Nsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Psoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dCsoildt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dNsoildt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%dPsoildt, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNCsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioPCsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNCsoilnew, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNCsoilmin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNCsoilmax, 1, mp, displs, blen, types, bidx)
    ! added by LN
    call add_address_1block(casapool%ratioNPplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNPlitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ratioNPsoil, 1, mp, displs, blen, types, bidx)

    ! ------- casaflux ----

    call add_address_1block(casaflux%Cgpp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Cnpp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Crp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Crgplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nminfix, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nminuptake, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Plabuptake, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Clabloss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fracClabile, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fracCalloc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fracNalloc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fracPalloc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Crmplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kplant_fire, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kplant_tot, 1, mp, displs, blen, types, bidx)
    ! 3D
    call add_address_1block(casaflux%fromPtoL, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fromPtoL_fire, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Cnep, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Crsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nmindep, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nminloss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nminleach, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nupland, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nlittermin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nsmin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nsimm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Nsnet, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fNminloss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fNminleach, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Pdep, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Pwea, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Pleach, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Ploss, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Pupland, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Plittermin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Psmin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Psimm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Psnet, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fPleach, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kplab, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kpsorb, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kpocc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%kmlabP, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Psorbmax, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%frac_sapwood, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%sapwood_area, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fHarvest, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fCrop, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%klitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%klitter_fire, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%klitter_tot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%ksoil, 1, mp, displs, blen, types, bidx)
    ! 3D
    call add_address_1block(casaflux%fromLtoS, 1, mp, displs, blen, types, bidx)
    ! 3D
    call add_address_1block(casaflux%fromStoS, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fromLtoCO2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fromStoCO2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxCtolitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxNtolitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxPtolitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxCtosoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxNtosoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxPtosoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fluxCtoCO2, 1, mp, displs, blen, types, bidx)
    ! 13C
    call add_address_1block(casaflux%FluxFromPtoL, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxFromLtoS, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxFromStoS, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxFromPtoCO2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxFromLtoCO2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxFromStoCO2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxFromPtoHarvest, 1, mp, displs, blen, types, bidx)

    ! ------- casamet ----

    call add_address_1block(casamet%glai, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%Tairk, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%precip, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%tsoilavg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%moistavg, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%btran, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%lnonwood, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%Tsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%moist, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%iveg2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%ijgcm, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%isorder, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%lat, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%lon, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%areacell, 1, mp, displs, blen, types, bidx)

    ! ------- casabal ----

    call add_address_1block(casabal%FCgppyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCnppyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCrmleafyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCrmwoodyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCrmrootyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCrgrowyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCrpyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCrsyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCneeyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FNdepyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FNfixyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FNsnetyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FNupyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FNleachyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FNlossyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FPweayear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FPdustyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FPsnetyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FPupyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FPleachyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FPlossyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%glaimon, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%glaimonx, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%cplantlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%nplantlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%pplantlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%clitterlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%nlitterlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%plitterlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%csoillast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%nsoillast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%psoillast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%nsoilminlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%psoillablast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%psoilsorblast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%psoilocclast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%cbalance, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%nbalance, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%pbalance, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%sumcbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%sumnbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%sumpbal, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%clabilelast, 1, mp, displs, blen, types, bidx)

    ! ------- phen -------

    call add_address_1block(phen%phase, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%TKshed, 1, mvtype, displs, blen, types, bidx)
    call add_address_1block(phen%doyphase, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%phen, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%aphen, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%phasespin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%doyphasespin_1, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%doyphasespin_2, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%doyphasespin_3, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%doyphasespin_4, 1, mp, displs, blen, types, bidx)

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE(*,*) 'worker ', rank, ' invalid number of casa_ts param fields ', bidx, ', fix it (21)!'
       CALL MPI_Abort(comm, 55, ierr)
    END IF

    CALL MPI_Type_create_struct(bidx, blen, displs, types, casa_ts, ierr)
    CALL MPI_Type_commit(casa_ts, ierr)

    CALL MPI_Type_size(casa_ts, tsize, ierr)
    CALL MPI_Type_get_extent(casa_ts, tmplb, text, ierr)

    WRITE(*,*) 'worker casa_ts param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...

    CALL MPI_Barrier(comm, ierr)

    ! so, now receive all the parameters
    !    CALL MPI_Recv (MPI_BOTTOM, 1, casa_ts, 0, 0, comm, stat, ierr)

    !   Maciej: buffered recv + unpac version
    ALLOCATE(rbuf(tsize))
    CALL MPI_Recv(rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    CALL MPI_Get_count(stat, casa_ts, rcount, ierr2)

    IF (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       CALL MPI_Unpack(rbuf, tsize, pos, MPI_BOTTOM, rcount, casa_ts, &
            comm, ierr)
       IF (ierr /= MPI_SUCCESS) write(*,*) 'casa params unpack error, rank: ',rank,ierr
    ELSE
       write(*,*) 'casa params recv rank err err2 rcount: ', rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)

    ! finally free the MPI type
    CALL MPI_Type_Free(casa_ts, ierr)

    ! all casa parameters have been received from the master by now
    RETURN

  END SUBROUTINE worker_casa_params


  ! MPI: creates inp_t type to receive input data from the master
  SUBROUTINE worker_intype (comm,met,veg)

    use mpi

    USE cable_def_types_mod
    use cable_mpicommon, only: ninput, add_address_1block

    IMPLICIT NONE

    ! Arguments
    INTEGER,INTENT(IN) :: comm
    TYPE(met_type),INTENT(IN):: met ! meteorological data
    TYPE(veg_parameter_type),INTENT(IN) :: veg ! LAI retrieved from file

    ! Local variables

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks
    INTEGER :: ierr

    INTEGER :: rank

    CALL MPI_Comm_rank (comm, rank, ierr)

    r1len = mp * extr1
    r2len = mp * extr2
    I1LEN = mp * extid
    llen = mp * extl

    ! max total number of fields to receive: met + veg fields
    ! ntyp = 10 + 1
    ntyp = ninput

    ALLOCATE (blocks(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    bidx = 0

    ! met fields

    call add_address_1block(met%fsd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%fdiff, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%tk, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%pmb, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%qv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%ua, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%precip, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%precip_sn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%fld, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%ca, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%coszen, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%Ndep, 1, mp, displs, blocks, types, bidx)

    ! veg fields

    call add_address_1block(veg%vlai, 1, mp, displs, blocks, types, bidx)

    ! additional field missing from previous versions;
    ! added when trying to fix a bug in the new mpi code
    ! the order of these new fields follows the order of their
    ! declaration in cable_define_types.F90

    call add_address_1block(met%year, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%moy, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%doy, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%hod, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%u10, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%rhum, 1, mp, displs, blocks, types, bidx)

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,': invalid intype nmat, nvec or n3d constant, fix it (22)!'
       CALL MPI_Abort(comm, 56, ierr)
    END IF


    ! marshall all fields into a single MPI derived datatype

    ! all variables are contiguous blocks of memory so just send them
    ! as blocks of bytes
    types = MPI_BYTE

    CALL MPI_Type_create_struct (bidx, blocks, displs, types, inp_t, ierr)
    CALL MPI_Type_commit (inp_t, ierr)

    CALL MPI_Type_size (inp_t, tsize, ierr)
    CALL MPI_Type_get_extent (inp_t, tmplb, text, ierr)

    WRITE (*,*) 'worker ',rank,': intype struct blocks, size, extent and lb: ', &
         bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blocks)

    RETURN

  END SUBROUTINE worker_intype

  ! MPI: creates send_t type to send the results to the master
  !
  ! list of fields that master needs to receive for use in write_output:
  !
  ! air%          rlam
  !
  ! canopy%       delwc, fe, fev, fh, fhs, fhv, fevw, fevc, fes, fnee, fpn, frday,
  !               frp, frs, ga, through, spill, tv, cansto,
  !
  ! met%          precip, precip_sn, fld, fsd, tk, pmb, qv, ua, ca,
  !
  ! rad%          albedo, qcan, qssabs, transd, flws
  !
  ! soil%         zse,
  !
  ! ssnow%        wb, snowd, osnowd, runoff, cls, rnof1, rnof2, tgg, tggsn, sdepth, isflag
  !
  ! veg%          vlai
  !
  ! Total: 47
  SUBROUTINE worker_outtype(comm,met,canopy,ssnow,rad,bal,air,soil,veg)

    use mpi

    USE cable_def_types_mod
    use cable_mpicommon, only: n3d, nmat, nvec, add_address_1block

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

    ! MPI: temp arrays for marshalling all types into a struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types
    INTEGER :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    INTEGER :: r1len, r2len, I1LEN, llen

    INTEGER :: rank, off, cnt
    INTEGER :: bidx, ierr

    INTEGER :: tsize
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

    CALL MPI_Comm_rank (comm, rank, ierr)

    ! MPI: calculate the sizes/extents of Fortran types used by
    ! CABLE
    CALL find_extents

    ! MPI: allocate temp vectors used for marshalling
    ntyp = n3d + nmat + nvec
    ALLOCATE (blocks(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    ! MPI: should work because worker versions of CABLE variables
    ! are allocated with starting index set to patch0
    !off = wpatch%patch0
    !cnt = wpatch%npatch
    ! MPI: new version, all arrays indices run 1:mp
    off = 1
    cnt = mp

    r1len = cnt * extr1
    r2len = cnt * extr2
    I1LEN = cnt * extid
    llen  = cnt * extl

    bidx = 0

    ! ------------- 3D arrays -------------

    ! rad 3D
    call add_address_1block(rad%qcan, 1, mp, displs, blocks, types, bidx)

    ! ------------- 2D arrays -------------

    ! MPI: an hvector type for each vector, maddr contains displacements
    ! for bundling these hvectors into the struct later
    ! block length is number of patches/worker * type extent
    ! stride is global number of patches * type extent
    ! repeat/no of blocks is the 2nd rank

    ! met 2D
    call add_address_1block(met%fsd, 1, mp, displs, blocks, types, bidx)

    ! canopy 2D
    ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
    call add_address_1block(canopy%evapfbl, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%gswx, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%zetar, 1, mp, displs, blocks, types, bidx)
    ! 13C
    call add_address_1block(canopy%An, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%Rd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%vcmax, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%gammastar, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%gsc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%gbc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%gac, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%ci, 1, mp, displs, blocks, types, bidx)

    ! ssnow 2D
    call add_address_1block(ssnow%dtmlt, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%albsoilsn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%gammzz, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%sconds, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%sdepth, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%smass, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%ssdn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%tgg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%tggsn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wb, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%evapfbl, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wbfice, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wbice, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wblf, 1, mp, displs, blocks, types, bidx)
    ! additional  for sli
    call add_address_1block(ssnow%S, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%Tsoil, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%thetai, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%snowliq, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%sconds, 1, mp, displs, blocks, types, bidx)    ! end additional for sli

    ! rad 2D
    call add_address_1block(rad%fbeam, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%albedo, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%fvlai, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%gradis, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%rhocdf, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%rniso, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%scalex, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%reffdf, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%reffbm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%extkbm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%extkdm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%cexpkbm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%cexpkdm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%rhocbm, 1, mp, displs, blocks, types, bidx)

    ! air 2D - all fields 1D - skipped

    ! soil 2D
    call add_address_1block(soil%albsoil, 1, mp, displs, blocks, types, bidx)

    ! veg 2D
    call add_address_1block(veg%refl, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%taul, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%froot, 1, mp, displs, blocks, types, bidx)


    ! ------------- 1D arrays -------------

    ! met
    call add_address_1block(met%ca, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%fld, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%precip, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%precip_sn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%tk, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%tvair, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%tvrad, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%pmb, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%ua, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%qv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%qvair, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%da, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%dva, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%coszen, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(met%fdiff, 1, mp, displs, blocks, types, bidx)

    ! canopy
    call add_address_1block(canopy%fess, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fesp, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%cansto, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%cduv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%delwc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%dewmm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%dgdtg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fe, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fh, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fpn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%A_sh, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%A_sl, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%A_slC, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%A_shC, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%A_slJ, 1, mp, displs, blocks, types, bidx)

    call add_address_1block(canopy%A_shJ, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%GPP_sh, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%GPP_sl, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%eta_A_cs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%dAdcs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%cs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%eta_GPP_cs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%eta_fevc_cs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%dlf, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vcmax_shade, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vcmax_sun, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%ejmax_shade, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%ejmax_sun, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frp, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frpw, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frpr, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fnee, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frday, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fnv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fev, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fevc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fevw, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fhv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fhvw, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fns, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fes, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fes_cor, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fhs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fhs_cor, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fwet, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%epot, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fnpp, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fevw_pot, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%gswx_T, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%cdtq, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%wetfac_cs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%ga, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%ghflux, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%precis, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%qscrn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%rnet, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%segg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%sghflux, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%spill, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%through, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%tscrn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%tv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%us, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%uscrn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%vlaiw, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%rghlai, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%wcint, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fwsoil, 1, mp, displs, blocks, types, bidx)
    ! 13C
    call add_address_1block(canopy%isc3, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%pudsto, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%pudsmx, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%cls, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%dfn_dtg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%dfh_dtg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%dfe_ddq, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%ddq_dtg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%evapsn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%fwtop, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%fwtop1, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%fwtop2, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%fwtop3, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%isflag, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%osnowd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%potev, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%pwb_min, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%runoff, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%rnof1, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%rnof2, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%rtsoil, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%snage, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%snowd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%smelt, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%ssdnn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%tss, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wbtot, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wb_lake, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%sinfil, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%qstss, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wetfac, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%owetfac, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%t_snwlr, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%tggav, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%otss, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%otss_0, 1, mp, displs, blocks, types, bidx)

    ! rad
    call add_address_1block(rad%extkb, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%extkd2, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%extkd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%flws, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%latitude, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%lwabv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%qssabs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%transd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%trad, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(rad%transb, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(bal%drybal, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(bal%osnowd0, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(bal%wbtot0, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(bal%wetbal, 1, mp, displs, blocks, types, bidx)

    ! air
    call add_address_1block(air%rho, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%volm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%rlam, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%qsat, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%epsi, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%visc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%psyc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%dsatdk, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%cmolar, 1, mp, displs, blocks, types, bidx)

    ! soil
    call add_address_1block(soil%bch, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%c3, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%clay, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%cnsd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%css, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%hsbh, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%hyds, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%i2bp3, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%ibp2, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%isoilm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%rhosoil, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%rs20, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%sand, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%sfc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%silt, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%ssat, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%sucs, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(soil%swilt, 1, mp, displs, blocks, types, bidx)

    ! veg
    call add_address_1block(veg%iveg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%meth, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vlai, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%canst1, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%ejmax, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%frac4, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%wai, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vegcf, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%tminvj, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%tmaxvj, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vbeta, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%xalbnir, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%hc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%shelrb, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vcmax, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%xfang, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%dleaf, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%rp20, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%rpcoef, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%extkn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%deciduous, 1, mp, displs, blocks, types, bidx)
    ! additional for SLI
    call add_address_1block(ssnow%Tsurface, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%h0, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%delwcol, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%evap, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%nsnow, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%nsteps, 1, mp, displs, blocks, types, bidx)
    ! end additional for SLI

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,': invalid outtype nmat, nvec or n3d constant, fix it (23)!'
       CALL MPI_Abort(comm, 57, ierr)
    END IF

    types = MPI_BYTE

    CALL MPI_Type_create_struct (bidx, blocks, displs, types, send_t, ierr)
    CALL MPI_Type_commit (send_t, ierr)

    CALL MPI_Type_size (send_t, tsize, ierr)
    CALL MPI_Type_get_extent (send_t, tmplb, text, ierr)

    WRITE (*,*) 'worker ',rank,': struct blocks, size, extent and lb: ',bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    !mcd287  CALL MPI_Reduce (tsize, tsize, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blocks)

    RETURN

  END SUBROUTINE worker_outtype


  subroutine worker_casa_type(comm, casabiome, casapool, casaflux, casamet, &
       casabal, phen)

#ifndef __INTEL__
    use mpi,             only: MPI_Address_kind, MPI_Abort, &
         MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
         MPI_Type_get_extent, MPI_Reduce, MPI_Datatype_null, MPI_Integer, &
         MPI_Sum
#else
    use mpi
#endif
    use cable_def_types_mod, only: mp, mvtype
    use casadimension,   only: mso
    use casavariable,    only: casa_biome, casa_pool, casa_flux, casa_met, &
         casa_balance, ncasa_biome, ncasa_pool, ncasa_flux, ncasa_met, &
         ncasa_bal
    use phenvariable,    only: phen_variable, ncasa_phen
    use cable_mpicommon, only: add_address_1block

    implicit none

    ! MPI communicator to talk to the workers
    integer,             intent(in)    :: comm
    type(casa_biome),    intent(inout) :: casabiome
    type(casa_pool),     intent(inout) :: casapool
    type(casa_flux),     intent(inout) :: casaflux
    type(casa_met),      intent(inout) :: casamet
    type(casa_balance),  intent(inout) :: casabal
    type(phen_variable), intent(inout) :: phen

    ! local variables

    ! MPI: temp arrays for marshalling all types into a struct
    integer, allocatable, dimension(:) :: blocks
    integer(MPI_Address_kind), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types
    integer :: ntyp ! number of worker's types

    integer :: bidx, ierr

    integer :: tsize
    integer(MPI_Address_kind) :: text, tmplb

    ! MPI: allocate temp vectors used for marshalling
    ntyp = ncasa_biome + ncasa_pool + ncasa_flux + ncasa_met + &
         ncasa_bal + ncasa_phen
    allocate(blocks(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    bidx = 0

    ! casabiome

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
     call add_address_1block(casabiome%xkplab, 1, mso, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkpsorb, 1, mso, displs, blocks, types, bidx)
     call add_address_1block(casabiome%xkpocc, 1, mso, displs, blocks, types, bidx)
     call add_address_1block(casabiome%prodptase, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%costnpup, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%maxfinelitter, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%maxcwd, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%nintercept, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%nslope, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%la_to_sa, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%vcmax_scalar, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%disturbance_interval, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_EnzPool, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_KMO2, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_KMcp, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_Ea, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%DAMM_alpha, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%plantrate, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%rmplant, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fracnpptoP, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fraclignin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fraclabile, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNCplantmin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNCplantmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNPplantmin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioNPplantmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%fracLigninplant, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ftransNPtoL, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ftransPPtoL, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%litterrate, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioPcplantmax, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%ratioPcplantmin, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(casabiome%soilrate, 1, mvtype, displs, blocks, types, bidx)

     ! casapool

     call add_address_1block(casapool%Clabile, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dClabiledt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Ctot, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Ctot_0, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Cplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Nplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Pplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dCplantdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dNplantdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dPplantdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNCplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNPplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Nsoilmin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Psoillab, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Psoilsorb, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Psoilocc, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dNsoilmindt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dPsoillabdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dPsoilsorbdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dPsoiloccdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Clitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Nlitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Plitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dClitterdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dNlitterdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dPlitterdt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNClitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNPlitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Csoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Nsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%Psoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dCsoildt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dNsoildt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%dPsoildt, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNCsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNPsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNCsoilnew, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNCsoilmin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioNCsoilmax, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioPCsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioPCplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casapool%ratioPClitter, 1, mp, displs, blocks, types, bidx)

     ! casaflux

     call add_address_1block(casaflux%Cgpp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Cnpp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Crp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Crgplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nminfix, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nminuptake, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Plabuptake, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Clabloss, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fracClabile, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%stemnpp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%frac_sapwood, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%sapwood_area, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Charvest, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nharvest, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Pharvest, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fharvest, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fcrop, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fracCalloc, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fracNalloc, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fracPalloc, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Crmplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kplant, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Cplant_turnover, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fromPtoL, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Cnep, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Crsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nmindep, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nminloss, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nminleach, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nupland, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nlittermin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nsmin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nsimm, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Nsnet, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fNminloss, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fNminleach, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Pdep, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Pwea, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Pleach, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Ploss, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Pupland, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Plittermin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Psmin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Psimm, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Psnet, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fPleach, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kplab, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kpsorb, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kpocc, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kmlabP, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Psorbmax, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Cplant_turnover_disturbance, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Cplant_turnover_crowding, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%Cplant_turnover_resource_limitation, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%klitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%ksoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fromLtoS, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fromStoS, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fromLtoCO2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fromStoCO2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtolitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxNtolitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxPtolitter, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtosoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxNtosoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxPtosoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtoCO2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtohwp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxNtohwp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxPtohwp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtoclear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxNtoclear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxPtoclear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%CtransferLUC, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fromPtoL_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%klitter_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%klitter_tot, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kplant_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%kplant_tot, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtoCO2_plant_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxCtoCO2_litter_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fluxfromPtoCO2_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%fluxfromLtoCO2_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxNtoAtm_fire, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromPtoL, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromLtoS, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromStoS, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromPtoCO2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromLtoCO2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromStoCO2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casaflux%FluxFromPtoHarvest, 1, mp, displs, blocks, types, bidx)

     ! casamet

     call add_address_1block(casamet%glai, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tairk, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%precip, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%tsoilavg, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistavg, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%btran, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%lnonwood, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moist, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%iveg2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%ijgcm, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%isorder, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%lat, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%lon, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%areacell, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tairkspin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%cgppspin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%crmplantspin_1, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%crmplantspin_2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%crmplantspin_3, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoilspin_1, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoilspin_2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoilspin_3, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoilspin_4, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoilspin_5, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%Tsoilspin_6, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistspin_1, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistspin_2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistspin_3, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistspin_4, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistspin_5, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%moistspin_6, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%mtempspin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%frecspin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%cAn12spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%cAn13spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%dprecip_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%aprecip_av20_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%du10_max_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%drhum_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%dtemp_max_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%dtemp_min_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%KBDI_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%D_MacArthur_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%FFDI_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%last_precip_spin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casamet%DSLR_spin, 1, mp, displs, blocks, types, bidx)

     ! casabal

     call add_address_1block(casabal%FCgppyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCnppyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCrmleafyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCrmwoodyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCrmrootyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCrgrowyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCrpyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCrsyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FCneeyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%dCdtyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%LAImax, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%Cleafmean, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%Crootmean, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FNdepyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FNfixyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FNsnetyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FNupyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FNleachyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FNlossyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FPweayear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FPdustyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FPsnetyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FPupyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FPleachyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%FPlossyear, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%glaimon, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%glaimonx, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%cplantlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%nplantlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%pplantlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%clitterlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%nlitterlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%plitterlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%csoillast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%nsoillast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%psoillast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%nsoilminlast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%psoillablast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%psoilsorblast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%psoilocclast, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%cbalance, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%nbalance, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%pbalance, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%sumcbal, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%sumnbal, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%sumpbal, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(casabal%clabilelast, 1, mp, displs, blocks, types, bidx)

     ! phen

     call add_address_1block(phen%Tkshed, 1, mvtype, displs, blocks, types, bidx)
     call add_address_1block(phen%phase, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%doyphase, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%phen, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%aphen, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%phasespin, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%doyphasespin_1, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%doyphasespin_2, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%doyphasespin_3, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(phen%doyphasespin_4, 1, mp, displs, blocks, types, bidx)

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker: invalid number of casa fields, fix it (24)!'
       write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
       call MPI_Abort(comm, 58, ierr)
    end if

    call MPI_Type_create_struct(bidx, blocks, displs, types, casa_t, ierr)
    call MPI_Type_commit(casa_t, ierr)

    call MPI_Type_size(casa_t, tsize, ierr)
    call MPI_Type_get_extent(casa_t, tmplb, text, ierr)

    write(*,*) 'casa type struct blocks, size, extent and lb: ', &
         bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    call MPI_Reduce(tsize, MPI_Datatype_null, 1, MPI_Integer, MPI_Sum, &
         0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blocks)

    return

  end subroutine worker_casa_type


  subroutine worker_climate_types(comm, climate, ktauday)

#ifndef __INTEL__
    use mpi, only: MPI_Address_kind, MPI_Get_address, MPI_Status_size, &
         MPI_Abort, &
         MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
         MPI_Type_get_extent, MPI_Reduce, MPI_Datatype_null, MPI_Integer, &
         MPI_Sum, MPI_Recv, MPI_Byte, MPI_Get_count, MPI_SUCCESS, &
         MPI_Unpack, MPI_Bottom
#else
    use mpi
#endif
    use cable_def_types_mod, only: climate_type, alloc_cbm_var, mp, &
         zero_cbm_var
    use cable_mpicommon, only: nclimate, add_address_1block

    implicit none

    integer,            intent(in)   :: comm
    type(climate_type), intent(inout):: climate
    integer,            intent(in)   :: ktauday

    ! MPI: temp arrays for marshalling all types into a struct
    integer, allocatable, dimension(:) :: blocks
    integer(MPI_Address_kind), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types
    integer :: ntyp ! number of worker's types

    integer :: tsize, totalrecv
    integer(MPI_Address_kind) :: text, tmplb

    integer :: rank
    integer :: bidx, ierr
    integer :: stat(MPI_Status_size), ierr2, rcount, pos

    character, dimension(:), allocatable :: rbuf

    !TRUNK no call to alloc_cbm_var
    call alloc_cbm_var(climate, mp, ktauday)
    call zero_cbm_var(climate)

    !TRUNK but call to climate_init
    !CALL climate_init (climate)
    ! MPI: allocate temp vectors used for marshalling
    ntyp = nclimate

    allocate(blocks(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    ! counter to sum total number of bytes receives from all workers
    totalrecv = 0

    bidx = 0

    ! ------------- scalars  -------------

    bidx = bidx + 1
    call MPI_Get_address(climate%nyear_average, displs(bidx), ierr)
    blocks(bidx) = extid
    types(bidx)  = MPI_BYTE

    bidx = bidx + 1
    call MPI_Get_address(climate%nday_average, displs(bidx), ierr)
    blocks(bidx) = extid
    types(bidx)  = MPI_BYTE

    bidx = bidx + 1
    call MPI_Get_address(climate%nyears, displs(bidx), ierr)
    blocks(bidx) = extid
    types(bidx)  = MPI_BYTE

    bidx = bidx + 1
    call MPI_Get_address(climate%doy, displs(bidx), ierr)
    blocks(bidx) = extid
    types(bidx)  = MPI_BYTE

    ! ------------- arrays -------------

    call add_address_1block(climate%chilldays, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%iveg, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%biome, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%GMD, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%modis_igbp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%DSLR, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%NDAY_Nesterov, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dtemp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_min, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_min20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_max20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%qtemp, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mmoist, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp_min, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%qtemp_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%qtemp_max_last_year, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp_min20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp_max20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%atemp_mean, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%AGDD5, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%GDD5, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%AGDD0, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%GDD0, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%alpha_PT, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%evap_PT, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%aevap, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%alpha_PT20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%GDD0_rec, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%frec, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dtemp_min, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%fdorm, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%fapar_ann_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%fapar_ann_max_last_year, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%AvgAnnMaxFAPAR, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dtemp_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%drhum, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%du10_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dprecip, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%aprecip, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%aprecip_av20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%last_precip, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%KBDI, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%FFDI, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%D_MacArthur, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Nesterov_Current, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Nesterov_ann_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Nesterov_ann_max_last_year, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Nesterov_ann_running_max, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp_min_20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%mtemp_max_20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_min_20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_max_20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dtemp_31, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dmoist_31, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%alpha_PT_20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%dtemp_91, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%APAR_leaf_sun, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%APAR_leaf_shade, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Dleaf_sun, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Dleaf_shade, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Tleaf_sun, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Tleaf_shade, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%cs_sun, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%cs_shade, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%scalex_sun, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%scalex_shade, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%fwsoil, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%aprecip_20, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Rd_sun, 1, mp, displs, blocks, types, bidx)
     call add_address_1block(climate%Rd_shade, 1, mp, displs, blocks, types, bidx)

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write (*,*) 'worker: invalid number of climate fields, fix it (25)!'
       call MPI_Abort (comm, 59, ierr)
    end if

    call MPI_Type_create_struct (bidx, blocks, displs, types, climate_t, ierr)
    call MPI_Type_commit (climate_t, ierr)

    call MPI_Type_size (climate_t, tsize, ierr)
    call MPI_Type_get_extent (climate_t, tmplb, text, ierr)

    call MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    ! so, now receive all the parameters
    !    CALL MPI_Recv (MPI_BOTTOM, 1, climate_t, 0, 0, comm, stat, ierr)

    !   Maciej: buffered recv + unpac version
    allocate (rbuf(tsize))
    call MPI_Recv (rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    call MPI_Get_count (stat, climate_t, rcount, ierr2)

    if (ierr == MPI_SUCCESS .and. ierr2 == MPI_SUCCESS .and. rcount == 1) then
       pos = 0
       call MPI_Unpack (rbuf, tsize, pos, MPI_BOTTOM, rcount, climate_t, &
            comm, ierr)
       if (ierr /= MPI_SUCCESS) write(*,*) 'climate unpack error, rank: ', rank, ierr
    else
       write(*,*) 'climate recv rank err err2 rcount: ', rank, ierr, ierr2, rcount
    end if

    deallocate(rbuf)

    deallocate(types)
    deallocate(displs)
    deallocate(blocks)

    return

  end subroutine worker_climate_types


  ! MPI: creates restart_t type to send to the master the fields
  ! that are only required for the restart file but not included in the
  ! results sent at the end of each time step
  subroutine worker_restart_type(comm, canopy, air, veg, ssnow)

#ifndef __INTEL__
    use mpi,                 only: MPI_Address_kind, MPI_Comm_rank, &
         MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, &
         MPI_Type_size, MPI_Type_get_extent, MPI_Reduce, MPI_Datatype_null, &
         MPI_Integer, MPI_Sum
#else
    use mpi
#endif
    use cable_def_types_mod, only: mp, ms, canopy_type, air_type, &
         veg_parameter_type, soil_snow_type
    use cable_mpicommon,     only: nrestart, add_address_1block

    implicit none

    integer,                  intent(in)    :: comm
    type(canopy_type),        intent(inout) :: canopy
    type(air_type),           intent(inout) :: air
    type(veg_parameter_type), intent(inout) :: veg
    type(soil_snow_type),     intent(inout) :: ssnow

    ! MPI: temp arrays for marshalling all types into a struct
    integer, allocatable, dimension(:) :: blocks
    integer(MPI_Address_kind), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types
    integer :: ntyp ! number of worker's types

    integer :: rank
    integer :: bidx, ierr

    integer :: tsize
    integer(MPI_Address_kind) :: text, tmplb

    call MPI_Comm_rank(comm, rank, ierr)

    ! MPI: allocate temp vectors used for marshalling
    ntyp = nrestart
    allocate(blocks(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    bidx = 0

    ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
    call add_address_1block(canopy%evapfbl, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%cduv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%dewmm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%dgdtg, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frpw, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%frpr, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fnv, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%rho, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%volm, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%qsat, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%epsi, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%visc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%psyc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%dsatdk, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(air%cmolar, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%otss, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(ssnow%wetfac, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%fwsoil, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(canopy%us, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%cfrd, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%vlai, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(veg%hc, 1, mp, displs, blocks, types, bidx)

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write (*,*) 'invalid nrestart constant, fix it (26)!'
       call MPI_Abort(comm, 60, ierr)
    end if

    call MPI_Type_create_struct(bidx, blocks, displs, types, restart_t, ierr)
    call MPI_Type_commit(restart_t, ierr)

    call MPI_Type_size(restart_t, tsize, ierr)
    call MPI_Type_get_extent(restart_t, tmplb, text, ierr)

    write (*,*) 'restart struct blocks, size, extent and lb: ', &
         rank, bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    call MPI_Reduce(tsize, MPI_Datatype_null, 1, MPI_Integer, &
         MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blocks)

    return

  end subroutine worker_restart_type


  SUBROUTINE worker_casa_dump_types(comm, casamet, casaflux, phen, climate, c13o2flux)

    use mpi

    use casadimension,       only: icycle, mplant, mphase
    use casavariable,        only: casa_met, casa_flux
    use cable_def_types_mod, only: mp, ms, climate_type
    use phenvariable
    use cable_mpicommon,     only: ncdumprw, add_address_1block
    ! 13C
    use cable_common_module, only: cable_user
    use cable_c13o2_def,     only: c13o2_flux

    IMPLICIT NONE

    ! sub arguments
    INTEGER,             INTENT(IN)    :: comm  ! MPI communicator
    TYPE(casa_met),      INTENT(INOUT) :: casamet
    TYPE(casa_flux),     INTENT(INOUT) :: casaflux
    TYPE(phen_variable), INTENT(INOUT) :: phen
    TYPE(climate_type),  INTENT(INOUT) :: climate
    ! 13C
    type(c13o2_flux),    intent(inout) :: c13o2flux

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: ierr

    INTEGER :: r1len, r2len, I1LEN ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank

    CALL MPI_Comm_rank(comm, rank, ierr)

    ntyp = ncdumprw + icycle - 1
    ! 13C
    if (cable_user%c13o2) ntyp = ntyp + 2
    if (cable_user%call_blaze) ntyp = ntyp + 11

    ALLOCATE(blen(ntyp))
    ALLOCATE(displs(ntyp))
    ALLOCATE(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid

    bidx = 0

    ! ------- casamet ----

    call add_address_1block(casamet%Tairk, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%Tsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casamet%moist, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Cgpp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%Crmplant, 1, mp, displs, blen, types, bidx)
    ! ------- phen ----

    call add_address_1block(phen%phase, 1, mp, displs, blen, types, bidx)
    call add_address_1block(phen%doyphase, 1, mp, displs, blen, types, bidx)
    ! ------- climate (for acclimation of autotrophic resp) ----

    call add_address_1block(climate%qtemp_max_last_year, 1, mp, displs, blen, types, bidx)
    call add_address_1block(climate%frec, 1, mp, displs, blen, types, bidx)
    ! ------- climate (for BLAZE) ----

    if (cable_user%call_blaze) then
       call add_address_1block(climate%dprecip, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%aprecip_av20, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%du10_max, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%drhum, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%dtemp_max, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%dtemp_min, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%KBDI, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%D_MacArthur, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%FFDI, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%DSLR, 1, mp, displs, blen, types, bidx)
       call add_address_1block(climate%last_precip, 1, mp, displs, blen, types, bidx)
    end if

    ! ------- casaflux - N and P deposition ---

    if (icycle>1) then
       call add_address_1block(casaflux%Nmindep, 1, mp, displs, blen, types, bidx)
    end if

    if (icycle>2) then
       call add_address_1block(casaflux%Pdep, 1, mp, displs, blen, types, bidx)
    end if

    ! ------- c13o2 ----

    if (cable_user%c13o2) then
       call add_address_1block(c13o2flux%cAn12, 1, mp, displs, blen, types, bidx)
       call add_address_1block(c13o2flux%cAn, 1, mp, displs, blen, types, bidx)
    end if

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE(*,*) 'worker ', rank, ' invalid number of casa_dump_t param fields ', bidx, ', fix it (27)!'
       CALL MPI_Abort(comm, 61, ierr)
    END IF

    CALL MPI_Type_create_struct(bidx, blen, displs, types, casa_dump_t, ierr)
    CALL MPI_Type_commit(casa_dump_t, ierr)

    CALL MPI_Type_size(casa_dump_t, tsize, ierr)
    CALL MPI_Type_get_extent(casa_dump_t, tmplb, text, ierr)

    WRITE (*,*) 'worker casa_dump_t param blocks, size, extent and lb: ', rank, &
         bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    ! write(*,*) 'b27 reduce wk', tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr
    ! flush(6)
    CALL MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

    ! ! if anything went wrong the master will mpi_abort
    ! ! which mpi_recv below is going to catch...
    ! ! so, now receive all the parameters
    ! CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, 0, comm, stat, ierr)
    !
    ! ! finally free the MPI type
    ! CALL MPI_Type_Free (casa_dump_t, ierr)

    ! all casa parameters have been received from the master by now

  END SUBROUTINE worker_casa_dump_types


  SUBROUTINE worker_casa_LUC_types(comm, casapool, casabal, casaflux)

    use mpi
    use cable_def_types_mod, only: mp
    USE casadimension, ONLY: mplant, mlitter, msoil
    USE casavariable, ONLY: casa_pool, casa_balance, casa_flux
    use cable_mpicommon, only: nLUCrw, add_address_1block

    IMPLICIT NONE

    ! sub arguments
    INTEGER,            INTENT(IN) :: comm  ! MPI communicator
    TYPE(casa_pool),    INTENT(IN) :: casapool
    TYPE(casa_balance), INTENT(IN) :: casabal
    TYPE(casa_flux),    INTENT(IN) :: casaflux

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: ierr

    INTEGER :: r1len, r2len, I1LEN ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank, off

    CALL MPI_Comm_rank(comm, rank, ierr)

    ntyp = nLUCrw

    ALLOCATE(blen(ntyp))
    ALLOCATE(displs(ntyp))
    ALLOCATE(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid

    off  = 1

    bidx = 0

    call add_address_1block(casapool%cplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%clitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%csoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%nplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%nlitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%nsoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%pplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%plitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%psoil, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%cplantlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%clitterlast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%csoillast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%Nsoilmin, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%clabile, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casapool%ctot, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%FCneeyear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casabal%clabilelast, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fHarvest, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%cHarvest, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%nHarvest, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%fcrop, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxCtohwp, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%FluxCtoclear, 1, mp, displs, blen, types, bidx)
    call add_address_1block(casaflux%CtransferLUC, 1, mp, displs, blen, types, bidx)

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,' invalid number of casa_LUC_t param fields ',bidx,', fix it (28)!'
       CALL MPI_Abort (comm, 62, ierr)
    END IF

    CALL MPI_Type_create_struct (bidx, blen, displs, types, casa_LUC_t, ierr)
    CALL MPI_Type_commit (casa_LUC_t, ierr)

    CALL MPI_Type_size (casa_LUC_t, tsize, ierr)
    CALL MPI_Type_get_extent (casa_LUC_t, tmplb, text, ierr)

    WRITE (*,*) 'worker casa_LUC_t param blocks, size, extent and lb: ',rank, &
         bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

    ! ! if anything went wrong the master will mpi_abort
    ! ! which mpi_recv below is going to catch...
    ! ! so, now receive all the parameters
    ! CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, 0, comm, stat, ierr)

    ! ! finally free the MPI type
    ! Call MPI_Type_Free (casa_dump_t, ierr)

    ! all casa parameters have been received from the master by now

  END SUBROUTINE worker_casa_LUC_types


  SUBROUTINE worker_pop_types(comm, veg, pop)

    use mpi
    USE POP_mpi
    USE POPmodule,          ONLY: POP_INIT
    USE POP_types,          ONLY: pop_type
    USE cable_def_types_mod,ONLY: mp, veg_parameter_type
    USE cable_common_module,ONLY: cable_user

    IMPLICIT NONE

    INTEGER,         INTENT(IN)    :: comm
    TYPE (pop_type), INTENT(INOUT) :: pop
    TYPE (veg_parameter_type),INTENT(IN) :: veg

    INTEGER, DIMENSION(:), Allocatable :: Iwood, ainv
    INTEGER :: mp_pop, stat(MPI_STATUS_SIZE), ierr

    INTEGER :: rank, inv

    CALL MPI_Comm_rank (comm, rank, ierr)

    ! Get POP relevant info from Master
    CALL MPI_Recv ( mp_pop, 1, MPI_INTEGER, 0, 0, comm, stat, ierr )
    write(*,*) 'worker iwood to allocate', rank, mp_pop, mp
    !write(*,*) 'worker mppop', rank, mp_pop
    !ALLOCATE( POP%Iwood( mp_pop ) )
    ALLOCATE( Iwood( mp_pop ) )
    write(*,*) 'worker iwood allocated', rank, mp_pop
    !CALL MPI_Recv ( POP%Iwood, mp_pop, MPI_INTEGER, 0, 0, comm, stat, ierr )
    CALL MPI_Recv( Iwood, mp_pop, MPI_INTEGER, 0, 0, comm, stat, ierr )
    !write(*,*),'worker Iwood', rank, POP%Iwood
    ! Maciej
    IF (ANY (Iwood < 1) .OR. ANY (Iwood > mp)) THEN
       write(*,*) 'worker iwood values outside valid ranges', rank
       inv = COUNT(Iwood < 1)
       IF (inv .gt. 0) THEN
          write(*,*) 'no of values below 1: ', inv
          ALLOCATE (ainv(inv))
          ainv = PACK (Iwood, Iwood .lt. 1)
          write (*,*) 'values below 1: ', ainv
          DEALLOCATE(ainv)
       END IF
       inv = COUNT(Iwood > mp)
       IF (inv .gt. 0) THEN
          write(*,*) 'no of values above mp ', mp, inv
          ALLOCATE(ainv(inv))
          ainv = PACK(Iwood, Iwood .gt. mp)
          write (*,*) 'values above mp: ', ainv
          DEALLOCATE (ainv)
       END IF

    END IF
    ! maciej: unnecessary, pop_init starts with a call to alloc_pop
    ! CALL ALLOC_POP ( POP, mp_pop )
    !CALL POP_INIT( pop,veg%disturbance_interval,mp_pop,POP%Iwood )
    ! maciej: veg%disturbance_levels received earlier in worker_cable_params
    CALL POP_INIT( pop,veg%disturbance_interval,mp_pop,Iwood )

    DEALLOCATE (Iwood)

    ! Maciej: pop_t used also for sending results back to master, even if from zero
    CALL create_pop_gridcell_type (pop_t, comm)

    IF (.NOT. CABLE_USER%POP_fromZero ) THEN
       write(*,*) 'rank receiving pop_grid from master', rank
       CALL MPI_Recv( POP%pop_grid(1), mp_pop, pop_t, 0, 0, comm, stat, ierr )
       CALL MPI_Recv( POP%it_pop(1), mp_pop, MPI_INTEGER, 0, 0, comm, stat, ierr )
    END IF

  END SUBROUTINE worker_pop_types


  ! 13C
  ! MPI: creates param_t type for the worker to receive the default casa parameters from the master process
  ! then receives them
  ! and finally frees the MPI type
  subroutine worker_c13o2_flux_params(comm, c13o2flux)

    use mpi
    ! use mpi,                 only: &
    !      MPI_ADDRESS_KIND, MPI_STATUS_SIZE, MPI_Comm_rank, MPI_BYTE, MPI_Get_address, &
    !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
    !      MPI_Type_get_extent, MPI_Reduce, MPI_DATATYPE_NULL, MPI_INTEGER, MPI_Sum, &
    !      MPI_Barrier, MPI_Recv, MPI_Get_count, MPI_SUCCESS, MPI_Unpack, MPI_BOTTOM, MPI_Type_Free
    use cable_def_types_mod, only: mp, mf
    use cable_c13o2_def,     only: c13o2_flux, c13o2_alloc_flux, c13o2_zero_flux
    use cable_mpicommon, only: nc13o2_flux, add_address_1block

    implicit none

    integer,          intent(in)    :: comm  ! mpi communicator
    type(c13o2_flux), intent(inout) :: c13o2flux

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    integer, allocatable, dimension(:) :: blen
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
    integer :: tsize

    integer :: stat(MPI_STATUS_SIZE), ierr
    integer :: c13o2_flux_ts

    integer :: r1len, r2len, i1len, llen ! block lengths
    integer :: bidx ! block index
    integer :: ntyp ! total number of blocks

    integer :: rank, off, ierr2, rcount, pos

    character, dimension(:), allocatable :: rbuf

    off = 1

    call MPI_Comm_rank(comm, rank, ierr)

    if (.not. associated(c13o2flux%ca)) then
       write(*,*) 'worker alloc c13o2_flux with m patches: ', rank, mp
       call c13o2_alloc_flux(c13o2flux, mp)
       call c13o2_zero_flux(c13o2flux)
    end if

    ntyp = nc13o2_flux

    allocate(blen(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid
    llen  = mp * extl

    bidx = 0

    ! -------

    call add_address_1block(c13o2flux%ca, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%cAn12, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%cAn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%RAn, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%Vstarch, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%Rstarch, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%An, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%Disc, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%Rsucrose, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2flux%Rphoto, 1, mp, displs, blen, types, bidx)
    ! -------

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker ', rank, ' invalid number of c13o2_flux_ts param fields ', bidx, ', fix it (29)!'
       call MPI_Abort(comm, 63, ierr)
    end if

    call MPI_Type_create_struct(bidx, blen, displs, types, c13o2_flux_ts, ierr)
    call MPI_Type_commit(c13o2_flux_ts, ierr)

    call MPI_Type_size(c13o2_flux_ts, tsize, ierr)
    call MPI_Type_get_extent(c13o2_flux_ts, tmplb, text, ierr)

    write(*,*) 'worker c13o2_flux_ts param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total data sent by all the workers
    call MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...
    call MPI_Barrier(comm, ierr)

    ! so, now receive all the parameters
    ! CALL MPI_Recv (MPI_BOTTOM, 1, c13o2_flux_ts, 0, 0, comm, stat, ierr)

    ! Maciej: buffered recv + unpac version
    allocate(rbuf(tsize))
    call MPI_Recv(rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    call MPI_Get_count(stat, c13o2_flux_ts, rcount, ierr2)

    if (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       call MPI_Unpack(rbuf, tsize, pos, MPI_BOTTOM, rcount, c13o2_flux_ts, comm, ierr)
       if (ierr /= MPI_SUCCESS) write(*,*) 'c13o2_flux params unpack error, rank: ',rank,ierr
    ELSE
       write(*,*) 'c13o2_flux params recv rank err err2 rcount: ', rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)

    ! finally free the MPI type
    call MPI_Type_Free(c13o2_flux_ts, ierr)

    ! all c13o2_flux parameters have been received from the master by now
    return

  end subroutine worker_c13o2_flux_params


  subroutine worker_c13o2_pool_params(comm, c13o2pools)

    use mpi
    ! use mpi,                 only: &
    !      MPI_ADDRESS_KIND, MPI_STATUS_SIZE, MPI_Comm_rank, MPI_BYTE, MPI_Get_address, &
    !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
    !      MPI_Type_get_extent, MPI_Reduce, MPI_DATATYPE_NULL, MPI_INTEGER, MPI_Sum, &
    !      MPI_Barrier, MPI_Recv, MPI_Get_count, MPI_SUCCESS, MPI_Unpack, MPI_BOTTOM, MPI_Type_Free
    use cable_def_types_mod, only: mp
    use casadimension,       only: mplant, mlitter, msoil
    use cable_c13o2_def,     only: c13o2_pool, c13o2_alloc_pools, c13o2_zero_pools
    use cable_mpicommon,     only: nc13o2_pool, add_address_1block

    implicit none

    integer,          intent(in)    :: comm  ! mpi communicator
    type(c13o2_pool), intent(inout) :: c13o2pools

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    integer, allocatable, dimension(:) :: blen
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
    integer :: tsize

    integer :: stat(MPI_STATUS_SIZE), ierr
    integer :: c13o2_pool_ts

    integer :: r1len, r2len, i1len, llen ! block lengths
    integer :: bidx ! block index
    integer :: ntyp ! total number of blocks

    integer :: rank, off, ierr2, rcount, pos

    character, dimension(:), allocatable :: rbuf

    off = 1

    call MPI_Comm_rank(comm, rank, ierr)

    if (.not. associated(c13o2pools%cplant)) then
       write(*,*) 'worker alloc c13o2_pool with m patches: ', rank, mp
       call c13o2_alloc_pools(c13o2pools, mp)
       call c13o2_zero_pools(c13o2pools)
    end if

    ntyp = nc13o2_pool

    allocate(blen(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid
    llen  = mp * extl

    bidx = 0

    ! -------

    call add_address_1block(c13o2pools%clabile, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2pools%charvest, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2pools%cplant, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2pools%clitter, 1, mp, displs, blen, types, bidx)
    call add_address_1block(c13o2pools%csoil, 1, mp, displs, blen, types, bidx)
    ! -------

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker ', rank, ' invalid number of c13o2_pool_ts param fields ', bidx, ', fix it (30)!'
       call MPI_Abort(comm, 64, ierr)
    end if

    call MPI_Type_create_struct(bidx, blen, displs, types, c13o2_pool_ts, ierr)
    call MPI_Type_commit(c13o2_pool_ts, ierr)

    call MPI_Type_size(c13o2_pool_ts, tsize, ierr)
    call MPI_Type_get_extent(c13o2_pool_ts, tmplb, text, ierr)

    write(*,*) 'worker c13o2_pool_ts param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total data sent by all the workers
    call MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...
    call MPI_Barrier(comm, ierr)

    ! so, now receive all the parameters
    ! call MPI_Recv(MPI_BOTTOM, 1, c13o2_pool_ts, 0, 0, comm, stat, ierr)

    ! Maciej: buffered recv + unpac version
    allocate(rbuf(tsize))
    call MPI_Recv(rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    call MPI_Get_count(stat, c13o2_pool_ts, rcount, ierr2)

    if (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       call MPI_Unpack(rbuf, tsize, pos, MPI_BOTTOM, rcount, c13o2_pool_ts, comm, ierr)
       if (ierr /= MPI_SUCCESS) write(*,*) 'c13o2_pool params unpack error, rank: ',rank,ierr
    else
       write(*,*) 'c13o2_pool params recv rank err err2 rcount: ', rank, ierr, ierr2, rcount
    end if

    deallocate(rbuf)

    ! finally free the MPI type
    call MPI_Type_Free(c13o2_pool_ts, ierr)

    ! all c13o2_pool parameters have been received from the master by now
    return

  end subroutine worker_c13o2_pool_params


  subroutine worker_c13o2_luc_params(comm, c13o2luc)

    use mpi
    ! use mpi,                 only: &
    !      MPI_ADDRESS_KIND, MPI_STATUS_SIZE, MPI_Comm_rank, MPI_BYTE, MPI_Get_address, &
    !      MPI_Abort, MPI_Type_create_struct, MPI_Type_commit, MPI_Type_size, &
    !      MPI_Type_get_extent, MPI_Reduce, MPI_DATATYPE_NULL, MPI_INTEGER, MPI_Sum, &
    !      MPI_Barrier, MPI_Recv, MPI_Get_count, MPI_SUCCESS, MPI_Unpack, MPI_BOTTOM, MPI_Type_Free
    use cable_def_types_mod, only: mland
    ! use cable_def_types_mod, only: mp
    use cable_c13o2_def,     only: c13o2_luc, c13o2_alloc_luc, c13o2_zero_luc
    use cable_mpicommon,     only: nc13o2_luc, add_address_1block

    implicit none

    integer,         intent(in)    :: comm  ! mpi communicator
    type(c13o2_luc), intent(inout) :: c13o2luc

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    integer, allocatable, dimension(:) :: blen
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb
    integer :: tsize

    integer :: stat(MPI_STATUS_SIZE), ierr
    integer :: c13o2_luc_ts

    integer :: r1len, r2len, i1len, llen ! block lengths
    integer :: bidx ! block index
    integer :: ntyp ! total number of blocks

    integer :: rank, off, ierr2, rcount, pos

    character, dimension(:), allocatable :: rbuf

    off = 1

    call MPI_Comm_rank(comm, rank, ierr)

    if (.not. associated(c13o2luc%charvest)) then
       write(*,*) 'worker alloc c13o2_luc with m land points: ', rank, mland
       call c13o2_alloc_luc(c13o2luc, mland)
       ! write(*,*) 'worker alloc c13o2_luc with m land points: ', rank, mp
       ! call c13o2_alloc_luc(c13o2luc, mp)
       call c13o2_zero_luc(c13o2luc)
    end if

    ntyp = nc13o2_luc

    allocate(blen(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mland * extr1
    r2len = mland * extr2
    i1len = mland * extid
    llen  = mland * extl
    ! r1len = mp * extr1
    ! r2len = mp * extr2
    ! i1len = mp * extid
    ! llen  = mp * extl

    bidx = 0

    ! -------

    call add_address_1block(c13o2luc%cagric, 1, mland, displs, blen, types, bidx)
    call add_address_1block(c13o2luc%charvest, 1, mland, displs, blen, types, bidx)
    call add_address_1block(c13o2luc%cclearance, 1, mland, displs, blen, types, bidx)

    ! -------

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker ', rank, ' invalid number of c13o2_luc_ts param fields ', bidx, ', fix it (31)!'
       call MPI_Abort(comm, 65, ierr)
    end if

    call MPI_Type_create_struct(bidx, blen, displs, types, c13o2_luc_ts, ierr)
    call MPI_Type_commit(c13o2_luc_ts, ierr)

    call MPI_Type_size(c13o2_luc_ts, tsize, ierr)
    call MPI_Type_get_extent(c13o2_luc_ts, tmplb, text, ierr)

    write(*,*) 'worker c13o2_luc_ts param blocks, size, extent and lb: ', rank, bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total data sent by all the workers
    call MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...
    call MPI_Barrier(comm, ierr)

    ! so, now receive all the parameters
    ! CALL MPI_Recv (MPI_BOTTOM, 1, c13o2_luc_ts, 0, 0, comm, stat, ierr)

    ! Maciej: buffered recv + unpac version
    allocate(rbuf(tsize))
    call MPI_Recv(rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    call MPI_Get_count(stat, c13o2_luc_ts, rcount, ierr2)

    if (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       call MPI_Unpack(rbuf, tsize, pos, MPI_BOTTOM, rcount, c13o2_luc_ts, comm, ierr)
       if (ierr /= MPI_SUCCESS) write(*,*) 'c13o2_luc params unpack error, rank: ',rank,ierr
    ELSE
       write(*,*) 'c13o2_luc params recv rank err err2 rcount: ', rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)

    ! finally free the MPI type
    call MPI_Type_Free(c13o2_luc_ts, ierr)

    ! all c13o2_luc parameters have been received from the master by now
    return

  end subroutine worker_c13o2_luc_params


  ! creates MPI types for sending c13o2 results back to the master at the end of the simulation
  subroutine worker_c13o2_flux_type(comm, c13o2flux)

    use mpi
    ! use mpi,                 only: &
    !      MPI_BYTE, MPI_ADDRESS_KIND, MPI_Get_address, MPI_Abort, MPI_Type_create_struct, &
    !      MPI_Type_commit, MPI_Type_size, MPI_Type_get_extent, MPI_Reduce, &
    !      MPI_DATATYPE_NULL, MPI_INTEGER, MPI_Sum
    use cable_def_types_mod, only: mp, mf
    use cable_c13o2_def,     only: c13o2_flux
    use cable_mpicommon,     only: nc13o2_flux, add_address_1block

    implicit none

    ! subroutine arguments
    integer,          intent(in)    :: comm ! MPI communicator to talk to the workers
    type(c13o2_flux), intent(inout) :: c13o2flux

    ! local variables

    ! MPI: temp arrays for marshalling all types into a struct
    integer, allocatable, dimension(:) :: blocks
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types
    integer :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    integer :: r1len, r2len, i1len, llen

    integer :: off, cnt
    integer :: bidx, ierr

    integer :: tsize
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb

    ! MPI: allocate temp vectors used for marshalling
    ntyp = nc13o2_flux
    allocate(blocks(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    !off = wpatch%patch0
    !cnt = wpatch%npatch
    off = 1
    cnt = mp

    r1len = cnt * extr1
    r2len = cnt * extr2
    i1len = cnt * extid
    llen  = cnt * extl

    bidx = 0

    ! ------------- 2D arrays -------------

    call add_address_1block(c13o2flux%An, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%Disc, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%Rsucrose, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%Rphoto, 1, mp, displs, blocks, types, bidx)
    ! ------------- 1D vectors -------------

    call add_address_1block(c13o2flux%ca, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%cAn12, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%cAn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%RAn, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%Vstarch, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2flux%Rstarch, 1, mp, displs, blocks, types, bidx)
    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker: invalid number of c13o2_flux fields, fix it (32)!'
       write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
       call MPI_Abort(comm, 66, ierr)
    end if

    types = MPI_BYTE

    call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_flux_t, ierr)
    call MPI_Type_commit(c13o2_flux_t, ierr)

    call MPI_Type_size(c13o2_flux_t, tsize, ierr)
    call MPI_Type_get_extent(c13o2_flux_t, tmplb, text, ierr)

    write(*,*) 'c13o2_flux type struct blocks, size, extent and lb: ', bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    call MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blocks)

    return

  end subroutine worker_c13o2_flux_type


  subroutine worker_c13o2_pool_type(comm, c13o2pools)

    use mpi
    ! use mpi,                 only: &
    !      MPI_BYTE, MPI_ADDRESS_KIND, MPI_Get_address, MPI_Abort, MPI_Type_create_struct, &
    !      MPI_Type_commit, MPI_Type_size, MPI_Type_get_extent, MPI_Reduce, &
    !      MPI_DATATYPE_NULL, MPI_INTEGER, MPI_Sum
    use cable_def_types_mod, only: mp
    use casadimension,       only: mplant, mlitter, msoil
    use cable_c13o2_def,     only: c13o2_pool
    use cable_mpicommon,     only: nc13o2_pool, add_address_1block

    implicit none

    ! subroutine arguments
    integer,          intent(in)    :: comm ! MPI communicator to talk to the workers
    type(c13o2_pool), intent(inout) :: c13o2pools

    ! local variables

    ! MPI: temp arrays for marshalling all types into a struct
    integer, allocatable, dimension(:) :: blocks
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types
    integer :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    integer :: r1len, r2len, i1len, llen

    integer :: off, cnt
    integer :: bidx, ierr

    integer :: tsize
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb

    ! MPI: allocate temp vectors used for marshalling
    ntyp = nc13o2_pool
    allocate(blocks(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    !off = wpatch%patch0
    !cnt = wpatch%npatch
    off = 1
    cnt = mp

    r1len = cnt * extr1
    r2len = cnt * extr2
    i1len = cnt * extid
    llen  = cnt * extl

    bidx = 0

    ! ------------- 2D arrays -------------

    call add_address_1block(c13o2pools%cplant, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2pools%clitter, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2pools%csoil, 1, mp, displs, blocks, types, bidx)
    ! ------------- 1D vectors -------------

    call add_address_1block(c13o2pools%clabile, 1, mp, displs, blocks, types, bidx)
    call add_address_1block(c13o2pools%charvest, 1, mp, displs, blocks, types, bidx)
    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker: invalid number of c13o2_pool fields, fix it (33)!'
       write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
       call MPI_Abort(comm, 67, ierr)
    end if

    types = MPI_BYTE

    call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_pool_t, ierr)
    call MPI_Type_commit(c13o2_pool_t, ierr)

    call MPI_Type_size(c13o2_pool_t, tsize, ierr)
    call MPI_Type_get_extent(c13o2_pool_t, tmplb, text, ierr)

    write(*,*) 'c13o2_pool type struct blocks, size, extent and lb: ', bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    call MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blocks)

    return

  end subroutine worker_c13o2_pool_type


  subroutine worker_c13o2_luc_type(comm, c13o2luc)

    use mpi
    ! use mpi,                 only: &
    !      MPI_BYTE, MPI_ADDRESS_KIND, MPI_Get_address, MPI_Abort, MPI_Type_create_struct, &
    !      MPI_Type_commit, MPI_Type_size, MPI_Type_get_extent, MPI_Reduce, &
    !      MPI_DATATYPE_NULL, MPI_INTEGER, MPI_Sum
    use cable_def_types_mod, only: mland
    ! use cable_def_types_mod, only: mp
    use cable_c13o2_def,     only: c13o2_luc
    use cable_mpicommon,     only: nc13o2_luc, add_address_1block

    implicit none

    ! subroutine arguments
    integer,          intent(in)   :: comm ! MPI communicator to talk to the workers
    type(c13o2_luc), intent(inout) :: c13o2luc

    ! local variables

    ! MPI: temp arrays for marshalling all types into a struct
    integer, allocatable, dimension(:) :: blocks
    integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: displs
    integer, allocatable, dimension(:) :: types
    integer :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    integer :: r1len, r2len, i1len, llen

    integer :: off, cnt
    integer :: bidx, ierr

    integer :: tsize
    integer(kind=MPI_ADDRESS_KIND) :: text, tmplb

    ! MPI: allocate temp vectors used for marshalling
    ntyp = nc13o2_luc
    allocate(blocks(ntyp))
    allocate(displs(ntyp))
    allocate(types(ntyp))

    !off = wpatch%patch0
    !cnt = wpatch%npatch
    off = 1
    cnt = mland
    ! cnt = mp

    r1len = cnt * extr1
    r2len = cnt * extr2
    i1len = cnt * extid
    llen  = cnt * extl

    bidx = 0

    ! ------------- 2D arrays -------------

    call add_address_1block(c13o2luc%charvest, 1, mland, displs, blocks, types, bidx)
    call add_address_1block(c13o2luc%cclearance, 1, mland, displs, blocks, types, bidx)

    ! ------------- 1D vectors -------------

    call add_address_1block(c13o2luc%cagric, 1, mland, displs, blocks, types, bidx)

    ! MPI: sanity check
    if (bidx /= ntyp) then
       write(*,*) 'worker: invalid number of c13o2_luc fields, fix it (34)!'
       write(*,*) 'ntyp: ', ntyp, 'bidx: ', bidx
       call MPI_Abort(comm, 68, ierr)
    end if

    types = MPI_BYTE

    call MPI_Type_create_struct(bidx, blocks, displs, types, c13o2_luc_t, ierr)
    call MPI_Type_commit(c13o2_luc_t, ierr)

    call MPI_Type_size(c13o2_luc_t, tsize, ierr)
    call MPI_Type_get_extent(c13o2_luc_t, tmplb, text, ierr)

    write(*,*) 'c13o2_luc type struct blocks, size, extent and lb: ', bidx, tsize, text, tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    call MPI_Reduce(tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_Sum, 0, comm, ierr)

    deallocate(types)
    deallocate(displs)
    deallocate(blocks)

    return

  end subroutine worker_c13o2_luc_type



SUBROUTINE worker_send_pop(POP, comm)

    use mpi
    use POP_mpi
    use POP_Types, only: POP_type

    implicit none

    integer,        intent(in) :: comm
    type(POP_type), intent(in) :: POP

    integer :: ierr

    if (POP%np .eq. 0) return

    call MPI_Send(POP%pop_grid(1), POP%np, pop_t, 0, 0, comm, ierr)
    call MPI_Send(POP%it_pop(1), POP%np, MPI_INTEGER, 0, 0, comm, ierr)

END SUBROUTINE worker_send_pop


! frees memory used for worker's data structures
SUBROUTINE worker_end(icycle, restart)

    use mpi
    use cable_common_module, only: cable_user

    implicit none

    integer, intent(in) :: icycle ! casa flag
    logical, intent(in) :: restart

    integer :: ierr

    call MPI_Type_free(inp_t, ierr)

    call MPI_Type_free(send_t, ierr)

    if (icycle>0) then
       call MPI_Type_free(casa_t, ierr)
       ! 13C
       if (cable_user%c13o2) then
          call MPI_Type_free(c13o2_flux_t, ierr)
          call MPI_Type_free(c13o2_pool_t, ierr)
       end if
    end if
    !MC - LUC is not freed. Does freeing objects close files as well?

    if (restart) then
       call MPI_Type_free(restart_t, ierr)
    end if

    return

END SUBROUTINE worker_end


!*********************************************************************************************

SUBROUTINE worker_spincasacnp(dels, kstart, kend, mloop, &
     veg, soil, casabiome, casapool, casaflux, casamet, casabal, phen, POP, climate, &
     LALLOC, c13o2flux, c13o2pools, icomm, ocomm)

  use cable_def_types_mod
  use cable_carbon_module
  use cable_common_module, only: cable_user
  use casadimension
  use casaparm
  use casavariable
  use phenvariable
  use POP_types,            only: POP_type
  use casa_cable,           only: POPdriver, analyticpool
  use casa_inout,           only: biogeochem
  use TypeDef,              only: dp
  ! 13C
  use cable_c13o2_def,      only: c13o2_pool, c13o2_flux
  use cable_c13o2,          only: c13o2_save_casapool, c13o2_update_pools, c13o2_sanity_pools
  use mo_isotope,           only: isoratio

  use mpi
  ! modules related to fire
  use blaze_mod,            only: type_blaze
  use simfire_mod,          only: type_simfire
  use blaze_drv,            only: blaze_driver
  use POP_constants,        only: rshootfrac

  implicit none

  !!CLN  character(len=99), intent(in)  :: fcnpspin
  real,                      intent(in)    :: dels
  integer,                   intent(in)    :: kstart
  integer,                   intent(in)    :: kend
  integer,                   intent(in)    :: mloop
  integer,                   intent(in)    :: lalloc
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
  ! 13c
  type(c13o2_flux),          intent(inout) :: c13o2flux
  type(c13o2_pool),          intent(inout) :: c13o2pools
  ! communicator for error-messages
  integer,                   intent(in)  :: icomm, ocomm

  ! local variables
  real(r_2), dimension(:), allocatable, save  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
  real(r_2), dimension(:), allocatable, save  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
  real(r_2), dimension(:), allocatable, save  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
  real(r_2), dimension(:), allocatable, save  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
  real(r_2), dimension(:), allocatable, save  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
  ! chris 12/oct/2012 for spin up casa
  real(r_2), dimension(:), allocatable, save  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass
  real(r_2), dimension(:), allocatable, save  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

  ! local variables
  integer                  :: myearspin, nyear, nloop1, loy
  integer                  :: ktau,ktauday,nday,idoy,ktauy,nloop
  real(r_2), dimension(mp) :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  real(r_2), dimension(mp) :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  real(r_2), dimension(mp) :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  real(r_2), dimension(mp) :: xnplimit,  xknlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  integer, allocatable :: iw(:) ! array of indices corresponding to woody (shrub or forest) tiles
  real(dp) :: rday

  integer :: stat(MPI_STATUS_SIZE)
  integer :: ierr

  ! blaze variables
  type(type_blaze)   :: blaze
  type(type_simfire) :: simfire

  ! 13C
  real(r_2), dimension(:,:), allocatable :: casasave
  real(r_2), dimension(:), allocatable :: avg_c13leaf2met, avg_c13leaf2str, avg_c13root2met, &
       avg_c13root2str, avg_c13wood2cwd

  if (.not.allocated(iw)) allocate(iw(POP%np))

  !! vh_js !!
  if (cable_user%call_POP) Iw = POP%iwood

  ktauday = int(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday
  loy     = 365

  ! 13C
  if (cable_user%c13o2) allocate(casasave(c13o2pools%ntile,c13o2pools%npools))

  ! chris 12/oct/2012 for spin up casa
  if (.not.(allocated(avg_cleaf2met)))  allocate(avg_cleaf2met(mp), avg_cleaf2str(mp), avg_croot2met(mp), avg_croot2str(mp), &
       avg_cwood2cwd(mp), &
       avg_nleaf2met(mp), avg_nleaf2str(mp), avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp), &
       avg_pleaf2met(mp), avg_pleaf2str(mp), avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp), &
       avg_cgpp(mp),      avg_cnpp(mp),      avg_nuptake(mp),   avg_puptake(mp), &
       avg_xnplimit(mp),  avg_xkNlimiting(mp), avg_xklitter(mp), avg_xksoil(mp), &
       avg_rationcsoilmic(mp),avg_rationcsoilslow(mp),avg_rationcsoilpass(mp), &
       avg_nsoilmin(mp),  avg_psoillab(mp),    avg_psoilsorb(mp), avg_psoilocc(mp))
  ! 13C - allocate in any case even if cable_user%c13o2==.false. to pass to analytic soil and litter pools
  allocate(avg_c13leaf2met(mp))
  allocate(avg_c13leaf2str(mp))
  allocate(avg_c13root2met(mp))
  allocate(avg_c13root2str(mp))
  allocate(avg_c13wood2cwd(mp))

  myearspin = cable_user%casa_spin_endyear - cable_user%casa_spin_startyear + 1

  ! compute the mean fluxes and residence time of each carbon pool
  avg_cleaf2met       = 0.0_dp
  avg_cleaf2str       = 0.0_dp
  avg_croot2met       = 0.0_dp
  avg_croot2str       = 0.0_dp
  avg_cwood2cwd       = 0.0_dp
  avg_nleaf2met       = 0.0_dp
  avg_nleaf2str       = 0.0_dp
  avg_nroot2met       = 0.0_dp
  avg_nroot2str       = 0.0_dp
  avg_nwood2cwd       = 0.0_dp
  avg_pleaf2met       = 0.0_dp
  avg_pleaf2str       = 0.0_dp
  avg_proot2met       = 0.0_dp
  avg_proot2str       = 0.0_dp
  avg_pwood2cwd       = 0.0_dp
  avg_cgpp            = 0.0_dp
  avg_cnpp            = 0.0_dp
  avg_nuptake         = 0.0_dp
  avg_puptake         = 0.0_dp
  avg_xnplimit        = 0.0_dp
  avg_xkNlimiting     = 0.0_dp
  avg_xklitter        = 0.0_dp
  avg_xksoil          = 0.0_dp
  avg_nsoilmin        = 0.0_dp
  avg_psoillab        = 0.0_dp
  avg_psoilsorb       = 0.0_dp
  avg_psoilocc        = 0.0_dp
  avg_rationcsoilmic  = 0.0_dp
  avg_rationcsoilslow = 0.0_dp
  avg_rationcsoilpass = 0.0_dp
  ! 13C
  if (cable_user%c13o2) then
      avg_c13leaf2met = 0.0_dp
      avg_c13leaf2str = 0.0_dp
      avg_c13root2met = 0.0_dp
      avg_c13root2str = 0.0_dp
      avg_c13wood2cwd = 0.0_dp
   end if

  do nyear=1, myearspin
     ! WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
     ! ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
     ! call read_casa_dump( ncfile,casamet, casaflux, phen,climate, ktau ,kend,.TRUE. )
     do idoy=1, mdyear
        ktau = (idoy-1)*ktauday + 1
        call MPI_Recv(MPI_BOTTOM, 1, casa_dump_t, 0, idoy, icomm, stat, ierr)
        ! 13C
        if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
        call biogeochem(idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
             casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
             xksoil,xkleaf,xkleafcold,xkleafdry, &
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd, &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd, &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
        ! 13C
        if (cable_user%c13o2) then
           avg_c13leaf2met(:) = avg_c13leaf2met(:) + &
                cleaf2met(:) * isoratio(c13o2pools%cplant(:,leaf), casasave(:,leaf), 0.0_dp, tiny(1.0_dp)) ! 1.0_dp
           avg_c13leaf2str(:) = avg_c13leaf2str(:) + &
                cleaf2str(:) * isoratio(c13o2pools%cplant(:,leaf), casasave(:,leaf), 0.0_dp, tiny(1.0_dp))
           avg_c13root2met(:) = avg_c13root2met(:) + &
                croot2met(:) * isoratio(c13o2pools%cplant(:,froot), casasave(:,froot), 0.0_dp, tiny(1.0_dp))
           avg_c13root2str(:) = avg_c13root2str(:) + &
                croot2str(:) * isoratio(c13o2pools%cplant(:,froot), casasave(:,froot), 0.0_dp, tiny(1.0_dp))
           avg_c13wood2cwd(:) = avg_c13wood2cwd(:) + &
                cwood2cwd(:) * isoratio(c13o2pools%cplant(:,wood), casasave(:,wood), 0.0_dp, tiny(1.0_dp))
           call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
           if (cable_user%c13o2) call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
        end if

        if (cable_user%call_POP .and. POP%np.gt.0) then ! CALL_POP
           ! accumulate annual variables for use in POP
           if (mod(ktau/ktauday,loy)==1) then
              ! (assumes 70% of wood NPP is allocated above ground)
              casaflux%stemnpp  = casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
              casabal%LAImax    = casamet%glai
              casabal%Cleafmean = casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
              casabal%Crootmean = casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
           else
              casaflux%stemnpp  = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
              casabal%LAImax    = max(casamet%glai, casabal%LAImax)
              casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
              casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
           end if

           if (idoy==mdyear) then ! end of year
              call POPdriver(casaflux, casabal, veg, POP)
           end if  ! end of year
        else
           casaflux%stemnpp = 0.0_dp
        end if ! CALL_POP

        !CLN CALL BLAZE_DRIVER(...)

        ! WHERE(xkNlimiting .eq. 0)  !Chris Lu 4/June/2012
        !    xkNlimiting = 0.001
        ! END WHERE

        avg_cleaf2met       = avg_cleaf2met       + cleaf2met
        avg_cleaf2str       = avg_cleaf2str       + cleaf2str
        avg_croot2met       = avg_croot2met       + croot2met
        avg_croot2str       = avg_croot2str       + croot2str
        avg_cwood2cwd       = avg_cwood2cwd       + cwood2cwd

        avg_nleaf2met       = avg_nleaf2met       + nleaf2met
        avg_nleaf2str       = avg_nleaf2str       + nleaf2str
        avg_nroot2met       = avg_nroot2met       + nroot2met
        avg_nroot2str       = avg_nroot2str       + nroot2str
        avg_nwood2cwd       = avg_nwood2cwd       + nwood2cwd

        avg_pleaf2met       = avg_pleaf2met       + pleaf2met
        avg_pleaf2str       = avg_pleaf2str       + pleaf2str
        avg_proot2met       = avg_proot2met       + proot2met
        avg_proot2str       = avg_proot2str       + proot2str
        avg_pwood2cwd       = avg_pwood2cwd       + pwood2cwd

        avg_cgpp            = avg_cgpp            + casaflux%cgpp
        avg_cnpp            = avg_cnpp            + casaflux%cnpp
        avg_nuptake         = avg_nuptake         + casaflux%Nminuptake
        avg_puptake         = avg_puptake         + casaflux%Plabuptake

        avg_xnplimit        = avg_xnplimit        + xnplimit
        avg_xkNlimiting     = avg_xkNlimiting     + xkNlimiting
        avg_xklitter        = avg_xklitter        + xklitter
        avg_xksoil          = avg_xksoil          + xksoil

        avg_nsoilmin        = avg_nsoilmin        + casapool%nsoilmin
        avg_psoillab        = avg_psoillab        + casapool%psoillab
        avg_psoilsorb       = avg_psoilsorb       + casapool%psoilsorb
        avg_psoilocc        = avg_psoilocc        + casapool%psoilocc

        avg_rationcsoilmic  = avg_rationcsoilmic  + casapool%ratioNCsoilnew(:,mic)
        avg_rationcsoilslow = avg_rationcsoilslow + casapool%ratioNCsoilnew(:,slow)
        avg_rationcsoilpass = avg_rationcsoilpass + casapool%ratioNCsoilnew(:,pass)
     end do ! idoy=1, mdyear
  end do ! nyear=1, myearspin

  !!CLN    CLOSE(91)
  ! average
  rday = 1.0_dp / real(nday*myearspin, dp)

  avg_cleaf2met       = avg_cleaf2met       * rday
  avg_cleaf2str       = avg_cleaf2str       * rday
  avg_croot2met       = avg_croot2met       * rday
  avg_croot2str       = avg_croot2str       * rday
  avg_cwood2cwd       = avg_cwood2cwd       * rday

  avg_nleaf2met       = avg_nleaf2met       * rday
  avg_nleaf2str       = avg_nleaf2str       * rday
  avg_nroot2met       = avg_nroot2met       * rday
  avg_nroot2str       = avg_nroot2str       * rday
  avg_nwood2cwd       = avg_nwood2cwd       * rday

  avg_pleaf2met       = avg_pleaf2met       * rday
  avg_pleaf2str       = avg_pleaf2str       * rday
  avg_proot2met       = avg_proot2met       * rday
  avg_proot2str       = avg_proot2str       * rday
  avg_pwood2cwd       = avg_pwood2cwd       * rday

  avg_cgpp            = avg_cgpp            * rday
  avg_cnpp            = avg_cnpp            * rday
  avg_nuptake         = avg_nuptake         * rday
  avg_puptake         = avg_puptake         * rday

  avg_xnplimit        = avg_xnplimit        * rday
  avg_xkNlimiting     = avg_xkNlimiting     * rday
  avg_xklitter        = avg_xklitter        * rday
  avg_xksoil          = avg_xksoil          * rday

  avg_nsoilmin        = avg_nsoilmin        * rday
  avg_psoillab        = avg_psoillab        * rday
  avg_psoilsorb       = avg_psoilsorb       * rday
  avg_psoilocc        = avg_psoilocc        * rday

  avg_rationcsoilmic  = avg_rationcsoilmic  * rday
  avg_rationcsoilslow = avg_rationcsoilslow * rday
  avg_rationcsoilpass = avg_rationcsoilpass * rday

  ! 13C
  if (cable_user%c13o2) then
     avg_c13leaf2met = avg_c13leaf2met * rday
     avg_c13leaf2str = avg_c13leaf2str * rday
     avg_c13root2met = avg_c13root2met * rday
     avg_c13root2str = avg_c13root2str * rday
     avg_c13wood2cwd = avg_c13wood2cwd * rday
  end if

  call analyticpool(veg,soil,casabiome,casapool, &
       casaflux,casamet,casabal, &
       avg_cleaf2met,avg_cleaf2str,avg_croot2met,avg_croot2str,avg_cwood2cwd, &
       avg_nleaf2met,avg_nleaf2str,avg_nroot2met,avg_nroot2str,avg_nwood2cwd, &
       avg_pleaf2met,avg_pleaf2str,avg_proot2met,avg_proot2str,avg_pwood2cwd, &
       avg_cnpp, &
       avg_xkNlimiting,avg_xklitter,avg_xksoil, &
       avg_ratioNCsoilmic,avg_ratioNCsoilslow,avg_ratioNCsoilpass, &
       avg_nsoilmin,avg_psoillab,avg_psoilsorb,avg_psoilocc, &
       avg_c13leaf2met, avg_c13leaf2str, avg_c13root2met, &
       avg_c13root2str, avg_c13wood2cwd, c13o2pools)
  if (cable_user%c13o2) call c13o2_sanity_pools(casapool, casaflux, c13o2pools)

  nloop1= max(1,mloop-3)

  do nloop=1, mloop
     !!CLN  OPEN(91,file=fcnpspin)
     !!CLN  read(91,*)
     do nyear=1, myearspin
       ! WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
       ! ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
       ! call read_casa_dump( ncfile, casamet, casaflux, phen,climate, ktau, kend, .TRUE. )
        do idoy=1, mdyear
           ktauy = idoy*ktauday
           ktau  = (idoy-1)*ktauday + 1
           call MPI_Recv(MPI_BOTTOM, 1, casa_dump_t, 0, idoy, icomm, stat, ierr)

           ! 13C
           if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
           call biogeochem(idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf, &
                xkleafcold,xkleafdry, &
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd, &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd, &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
           ! 13C
           if (cable_user%c13o2) then
              call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
              call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
           end if

           !! CLN BLAZE here
           !!CLNblaze_spin_year = 1900 - myearspin + nyear
           !CLN CALL BLAZE_DRIVER(casapool,casaflux,shootfrac, met..., idoy,blaze_spin_year, DO_BLAZE_TO

           if (cable_user%call_POP .and. POP%np.gt.0) then ! CALL_POP
              ! accumulate annual variables for use in POP
              if (mod(ktau/ktauday,LOY) == 1) then
                 casaflux%stemnpp  =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
                 ! (assumes 70% of wood NPP is allocated above ground)
                 casabal%LAImax    = casamet%glai
                 casabal%Cleafmean = casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
                 casabal%Crootmean = casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
              else
                 casaflux%stemnpp  = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
                 casabal%LAImax    = max(casamet%glai, casabal%LAImax)
                 casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1) / real(LOY,dp) / 1000.0_dp
                 casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3) / real(LOY,dp) / 1000.0_dp
              end if

              if (idoy==mdyear) then ! end of year
                 call POPdriver(casaflux, casabal, veg, POP)
                 !CLN Check here accounting missing
                 if (cable_user%call_blaze) then
                    call blaze_driver(blaze%ncells, blaze, simfire, casapool, casaflux, &
                         casamet, climate, rshootfrac, idoy, 1900, 1, POP, veg)
                 end if
                 !! CLN BLAZE TURNOVER
              end if  ! end of year
           else
              casaflux%stemnpp = 0.0_dp
           end if ! CALL_POP
           write(wlogn,*) 'idoy ', idoy

        end do   ! end of idoy

     end do   ! end of nyear

  end do     ! end of nloop
  write(wlogn,*) 'b4 MPI_SEND'
  CALL MPI_Send(MPI_BOTTOM, 1, casa_t, 0, 0, ocomm, ierr)
  ! 13C
  if (cable_user%c13o2) then
     call MPI_Send(MPI_BOTTOM, 1, c13o2_pool_t, 0, 0, ocomm, ierr)
  end if
  write(wlogn,*) 'after MPI_SEND'
  IF(CABLE_USER%CALL_POP) CALL worker_send_pop(POP, ocomm)
  write(wlogn,*) 'cplant ', casapool%cplant(1,:)

END SUBROUTINE worker_spincasacnp

!*********************************************************************************************

SUBROUTINE worker_CASAONLY_LUC(dels, kstart, kend, veg, soil, casabiome, casapool, &
     casaflux, casamet, casabal, phen, POP, climate, LALLOC, &
     c13o2flux, c13o2pools, icomm, ocomm)

  USE cable_def_types_mod
  USE cable_carbon_module
  USE cable_common_module, ONLY: cable_user
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE POP_Types,           only: POP_TYPE
  use casa_cable,          only: POPdriver
  use casa_inout,          only: biogeochem
  USE TypeDef,             ONLY: dp
  use mpi
  ! 13C
  use cable_c13o2_def,     only: c13o2_flux, c13o2_pool
  use cable_c13o2,         only: c13o2_save_casapool, c13o2_update_pools, &
       c13o2_sanity_pools

  IMPLICIT NONE

  !!CLN  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
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
  ! 13c
  type(c13o2_flux),          intent(inout) :: c13o2flux
  type(c13o2_pool),          intent(inout) :: c13o2pools
  ! communicator for error-messages
  integer,                   intent(in)    :: icomm, ocomm

  ! local variables
  integer           :: myearspin, nyear
  integer           :: ktau, ktauday, nday, idoy
  real(r_2), dimension(mp) :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  real(r_2), dimension(mp) :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  real(r_2), dimension(mp) :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  real(r_2), dimension(mp) :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  real(dp) :: StemNPP(mp,2)

  integer :: stat(MPI_STATUS_SIZE)
  integer :: ierr, rank
  integer :: yyyy
  ! 13C
  real(dp), dimension(:,:), allocatable :: casasave

  ! 13C
  if (cable_user%c13o2) allocate(casasave(c13o2pools%ntile,c13o2pools%npools))

  ktauday = int(24.0*3600.0/dels)
  nday    = (kend-kstart+1)/ktauday

  myearspin = CABLE_USER%YEAREND - CABLE_USER%YEARSTART + 1
  yyyy      = CABLE_USER%YEARSTART - 1

  do nyear=1, myearspin
     do idoy=1, mdyear
        !CASAONLY_LUC ktau = (idoy-1)*ktauday + ktauday
        ktau = (idoy-1)*ktauday + 1
        CALL MPI_Recv(MPI_BOTTOM, 1, casa_dump_t, 0, idoy, icomm, stat, ierr)

        ! 13C
        if (cable_user%c13o2) call c13o2_save_casapool(casapool, casasave)
        CALL biogeochem(idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
             casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
             xksoil,xkleaf,xkleafcold,xkleafdry, &
             cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd, &
             nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd, &
             pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
        ! 13C
        if (cable_user%c13o2) then
           call c13o2_update_pools(casasave, casaflux, c13o2flux, c13o2pools)
           call c13o2_sanity_pools(casapool, casaflux, c13o2pools)
        end if

        !! CLN BLAZE here

        ! accumulate annual variables for use in POP
        IF (idoy==1) THEN
           ! (assumes 70% of wood NPP is allocated above ground)
           casaflux%stemnpp  =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
           casabal%LAImax    = casamet%glai
           casabal%Cleafmean = casapool%cplant(:,1) / real(mdyear,dp) / 1000.0_dp
           casabal%Crootmean = casapool%cplant(:,3) / real(mdyear,dp) / 1000.0_dp
        ELSE
           casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_dp
           casabal%LAImax    = max(casamet%glai, casabal%LAImax)
           casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1) / real(mdyear,dp) / 1000.0_dp
           casabal%Crootmean = casabal%Crootmean +casapool%cplant(:,3) / real(mdyear,dp) / 1000.0_dp
        END IF

        IF (idoy==mdyear) THEN ! end of year
           write(wlogn,*) 'b4 MPI_SEND, casa_LUC_t ', casapool%cplant(:,2)
           flush(wlogn)
           CALL MPI_Send(MPI_BOTTOM, 1, casa_t, 0, 0, ocomm, ierr)
           CALL MPI_Send(MPI_BOTTOM, 1, casa_LUC_t, 0, 0, ocomm, ierr)
           ! 13C
           if (cable_user%c13o2) then
              call MPI_Send(MPI_BOTTOM, 1, c13o2_flux_t, 0, 0, ocomm, ierr)
              call MPI_Send(MPI_BOTTOM, 1, c13o2_pool_t, 0, 0, ocomm, ierr)
              call MPI_Send(MPI_BOTTOM, 1, c13o2_luc_t, 0, 0, ocomm, ierr)
           end if
           write(wlogn,*) 'after MPI_SEND, casa_LUC_t ', casapool%cplant(:,2)
           flush(wlogn)
           StemNPP(:,1) = casaflux%stemnpp
           StemNPP(:,2) = 0.0_dp

           CALL MPI_Comm_rank(icomm, rank, ierr)
           write(wlogn,*)
           write(wlogn,*) 'rank receiving pop_grid from master ', rank
           ! write(wlogn,*) 'b4 MPI_Recv, pop_t cmass: ', POP%pop_grid%cmass_sum
           ! write(wlogn,*) 'b4 MPI_Recv, pop_t LU: ', POP%pop_grid%LU
           CALL MPI_Recv(POP%pop_grid(1), POP%np, pop_t, 0, 0, icomm, stat, ierr )
           ! write(wlogn,*)
           ! write(wlogn,*) 'after MPI_Recv, pop_t cmass: ', POP%pop_grid%cmass_sum
           write(wlogn,*) 'after MPI_Recv, pop_t'
           flush(wlogn)
           IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP
              write(wlogn,*) 'b4  POPdriver ', POP%pop_grid%cmass_sum
              CALL POPdriver(casaflux, casabal, veg, POP)
              !! CLN BLAZE TURNOVER
           END IF
           ! write(wlogn,*)
           ! write(wlogn,*) 'after POPstep cmass: ', POP%pop_grid%cmass_sum
           write(wlogn,*) 'after POPstep ',  POP%pop_grid%cmass_sum
           flush(wlogn)
           CALL worker_send_pop(POP, ocomm)
           write(wlogn,*) 'after worker_send_pop'
           flush(wlogn)
        END IF

     end do ! idoy=1,mdyear

     ! receive updates to CASA pools resulting from LUC
     write(wlogn,*)
     write(wlogn,*) 'b4 mpi_recv casa_LUC_t'
     !NEW CALL MPI_Recv(MPI_BOTTOM, 1, casa_t, 0, nyear, icomm, stat, ierr)
     CALL MPI_Recv(MPI_BOTTOM, 1, casa_LUC_t, 0, nyear, icomm, stat, ierr)
     ! 13C
     if (cable_user%c13o2) then
        call MPI_Recv(MPI_BOTTOM, 1, c13o2_flux_t, 0, nyear, icomm, stat, ierr)
        call MPI_Recv(MPI_BOTTOM, 1, c13o2_pool_t, 0, nyear, icomm, stat, ierr)
        call MPI_Recv(MPI_BOTTOM, 1, c13o2_luc_t, 0, nyear, icomm, stat, ierr)
     end if
     write(wlogn,*) 'after mpi_recv casa_LUC_t'
  end do

END SUBROUTINE WORKER_CASAONLY_LUC


END MODULE cable_mpiworker
