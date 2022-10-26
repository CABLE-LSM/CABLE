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
!                 casa_cable
!                 casa_inout_module
!
! CALLs:       
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

  USE cable_mpicommon
  USE cable_common_module,  ONLY: cable_user
  USE casa_inout_module
  USE casa_cable


  IMPLICIT NONE

  SAVE



  PRIVATE

  ! MPI: MPI derived datatype for receiving parameters from the master
  INTEGER :: param_t

  ! MPI: MPI derived datatype for receiving casa parameters from the master
  INTEGER :: casaparam_t

  ! MPI: MPI derived datatype for receiving input from the master
  INTEGER :: inp_t

  ! MPI: MPI derived datatype for sending results back to the master
  INTEGER :: send_t

  ! worker's struct for sending final casa results to the master
  INTEGER :: casa_t

  ! worker's struct for rec'ing/sending final casa results to/from the master
  INTEGER :: casa_dump_t

  ! worker's struct for rec'ing/sending casa pools to/from the master (for LUC calcs)
  INTEGER :: casa_LUC_t

  ! worker's struct for rec'ing/sending final casa results to/from the master
  INTEGER :: climate_t

  ! worker's struct for rec'ing/sending pop io to/from the master
  INTEGER :: pop_t

  ! worker's struct for restart data to the master
  INTEGER :: restart_t

  ! worker's logfile unit
  !INTEGER :: wlogn
  !debug moved to iovars -- easy to pass around

  PUBLIC :: mpidrv_worker
  REAL, allocatable  :: heat_cap_lower_limit(:,:)

CONTAINS

  SUBROUTINE mpidrv_worker (comm)

    USE mpi

    USE cable_def_types_mod
    USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps, globalMetfile,  &
         verbose, fixedCO2,output,check,patchout,    &
         patch_type,soilparmnew,&
         defaultLAI, wlogn
    USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user,     &
         cable_runtime, filename, myhome,            &
         redistrb, wiltParam, satuParam, CurYear,    &
         IS_LEAPYEAR, calcsoilalbedo,                &
         kwidth_gl, gw_params
  USE casa_ncdf_module, ONLY: is_casa_time
    USE cable_input_module,   ONLY: open_met_file,load_parameters,              &
         get_met_data,close_met_file
    USE cable_output_module,  ONLY: create_restart,open_output_file,            &
         write_output,close_output_file
    USE cable_cbm_module
    USE cable_climate_mod

    ! modules related to CASA-CNP
    USE casadimension,        ONLY: icycle
    USE casavariable,         ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
         casa_met, casa_balance
    USE phenvariable,         ONLY: phen_variable

    !CLN added
    ! modules related to POP
    USE POPmodule,            ONLY: POP_INIT
    USE POP_Types,            ONLY: POP_TYPE
    USE POP_Constants,        ONLY: HEIGHT_BINS, NCOHORT_MAX

    ! PLUME-MIP only
    USE CABLE_PLUME_MIP,      ONLY: PLUME_MIP_TYPE

    USE cable_namelist_util, ONLY : get_namelist_file_name,&
         CABLE_NAMELIST,arg_not_namelist
    
USE cbl_soil_snow_init_special_module
    IMPLICIT NONE

    ! MPI:
    INTEGER               :: comm ! MPI communicator for comms with the workers

    ! CABLE namelist: model configuration, runtime/user switches
    !CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml'

    ! timing variables
    INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
    INTEGER, PARAMETER ::  mloop  = 30   ! CASA-CNP PreSpinup loops

    INTEGER        ::                                                           &
         ktau,       &  ! increment equates to timestep, resets if spinning up
         ktau_tot,   &  ! NO reset when spinning up, total timesteps by model
         kend,       &  ! no. of time steps in run
                                !CLN      kstart = 1, &  ! timestep to start at
         koffset = 0, &  ! timestep to start at
         ktauday,    &  ! day counter for CASA-CNP
         idoy,       &  ! day of year (1:365) counter for CASA-CNP
         nyear,      &  ! year counter for CASA-CNP
         casa_it,    &  ! number of calls to CASA-CNP
         ctime,      &  ! day count for casacnp
         YYYY,       &  !
         LOY,        &  ! Length of Year
         count_sum_casa, & ! number of time steps over which casa pools &
                                !and fluxes are aggregated (for output)
         rank,       &  ! Rank of this worker
         maxdiff(2)     ! location of maximum in convergence test

    REAL      :: dels    ! time step size in seconds
    CHARACTER :: cRank*4 ! for worker-logfiles

    ! CABLE variables
    TYPE (met_type)       :: met     ! met input variables
    TYPE (air_type)       :: air     ! air property variables
    TYPE (canopy_type)    :: canopy  ! vegetation variables
    TYPE (radiation_type) :: rad     ! radiation variables
    TYPE (roughness_type) :: rough   ! roughness varibles
    TYPE (balances_type)  :: bal     ! energy and water balance variables
    TYPE (soil_snow_type) :: ssnow   ! soil and snow variables
    TYPE (climate_type)   :: climate     ! climate variables

    ! CABLE parameters
    TYPE (soil_parameter_type) :: soil ! soil parameters
    TYPE (veg_parameter_type)  :: veg  ! vegetation parameters

    TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
    TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables

    ! CASA-CNP variables
    TYPE (casa_biome)     :: casabiome
    TYPE (casa_pool)      :: casapool
    TYPE (casa_flux)      :: casaflux
    TYPE (casa_met)       :: casamet
    TYPE (casa_balance)   :: casabal
    TYPE (phen_variable)  :: phen
    TYPE (POP_TYPE)       :: POP
    TYPE (PLUME_MIP_TYPE) :: PLUME
    CHARACTER             :: cyear*4
    CHARACTER             :: ncfile*99

    ! declare vars for switches (default .FALSE.) etc declared thru namelist
    LOGICAL, SAVE           :: &
         vegparmnew    = .FALSE., & ! using new format input file (BP dec 2007)
         spinup        = .FALSE., & ! model spinup to soil state equilibrium?
         spinConv      = .FALSE., & ! has spinup converged?
         spincasa      = .FALSE., & ! TRUE: CASA-CNP Will spin mloop times,
         l_casacnp     = .FALSE., & ! using CASA-CNP with CABLE
         l_landuse     = .FALSE., & ! using LANDUSE                 
         l_laiFeedbk   = .FALSE., & ! using prognostic LAI
         l_vcmaxFeedbk = .FALSE., & ! using prognostic Vcmax
         CASAONLY      = .FALSE., & ! ONLY Run CASA-CNP
         CALL1         = .TRUE.

    REAL              :: &
         delsoilM,         & ! allowed variation in soil moisture for spin up
         delsoilT            ! allowed variation in soil temperature for spin up

    ! temporary storage for soil moisture/temp. in spin up mode
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: &
         soilMtemp,                         &
         soilTtemp

    ! MPI:
    LOGICAL :: loop_exit     ! MPI: exit flag for bcast to workers
    INTEGER :: stat(MPI_STATUS_SIZE)
    INTEGER :: icomm ! separate dupes of MPI communicator for send and recv
    INTEGER :: ocomm ! separate dupes of MPI communicator for send and recv
    INTEGER :: ierr
    CHARACTER(len=200):: Run

    ! switches etc defined thru namelist (by default cable.nml)
    NAMELIST/CABLE/                  &
         filename,         & ! TYPE, containing input filenames
         vegparmnew,       & ! use new soil param. method
         soilparmnew,      & ! use new soil param. method
         calcsoilalbedo,   & ! switch: soil colour albedo - Ticket #27
         spinup,           & ! spinup model (soil) to steady state
         delsoilM,delsoilT,& !
         output,           &
         patchout,         &
         check,            &
         verbose,          &
         leaps,            &
         logn,             &
         fixedCO2,         &
         spincasa,         &
         l_casacnp,        &
         l_landuse,        &
         l_laiFeedbk,      &
         l_vcmaxFeedbk,    &
         icycle,           &
         casafile,         &
         ncciy,            &
         gswpfile,         &
         globalMetfile,    &   
         redistrb,         &
         wiltParam,        &
         satuParam,        &
         cable_user,       &  ! additional USER switches
         gw_params

    INTEGER :: i,x,kk
    INTEGER :: LALLOC, iu
!For consistency w JAC
  REAL,ALLOCATABLE, SAVE :: c1(:,:)
  REAL,ALLOCATABLE, SAVE :: rhoch(:,:)
  REAL,ALLOCATABLE, SAVE :: xk(:,:)
    ! END header

    ! Maciej: make sure the variable does not go out of scope
    mp = 0

    ! Open, read and close the namelist file.
    OPEN( 10, FILE = CABLE_NAMELIST )
    READ( 10, NML=CABLE )   !where NML=CABLE defined above
    CLOSE(10)

    IF( IARGC() > 0 ) THEN
       CALL GETARG(1, filename%met)
       CALL GETARG(2, casafile%cnpipool)
    ENDIF


    IF (CABLE_USER%POPLUC .AND. TRIM(CABLE_USER%POPLUC_RunType) .EQ. 'static') &
         CABLE_USER%POPLUC= .FALSE.

    ! Get worker's rank and determine logfile-unit

    ! MPI: TODO: find a way to preserve workers log messages somewhere
    ! (either separate files or collated by the master to a single file
    ! or perhaps use MPI-IO - but probably not gonna work with random length
    ! text strings)
    ! LN: Done!
    IF ( CABLE_USER%LogWorker ) THEN
       CALL MPI_Comm_rank (comm, rank, ierr)
       WRITE(cRank,FMT='(I4.4)')rank
       wlogn = 1000+rank
       OPEN(wlogn,FILE="cable_log_"//cRank,STATUS="REPLACE")
    ELSE
       wlogn = 1000
       OPEN(wlogn, FILE="/dev/null")
    ENDIF

    ! INITIALISATION depending on nml settings

    CurYear = CABLE_USER%YearStart

    IF ( icycle .GE. 11 ) THEN
       icycle                     = icycle - 10
       CASAONLY                   = .TRUE.
       CABLE_USER%CASA_DUMP_READ  = .TRUE.
       CABLE_USER%CASA_DUMP_WRITE = .FALSE.
    ELSEIF ( icycle .EQ. 0 ) THEN
       CABLE_USER%CASA_DUMP_READ  = .FALSE.
       spincasa                   = .FALSE.
       CABLE_USER%CALL_POP        = .FALSE.
    ENDIF

    !! vh_js !!
    IF (icycle.GT.0) THEN
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

    IF ( TRIM(cable_user%MetType) .EQ. 'gpgs' ) THEN
       leaps = .TRUE.
       cable_user%MetType = 'gswp'
    ENDIF

    cable_runtime%offline = .TRUE.

    IF( l_casacnp  .AND. ( icycle == 0 .OR. icycle > 3 ) )                   &
         STOP 'icycle must be 1 to 3 when using casaCNP'
    IF( ( l_laiFeedbk .OR. l_vcmaxFeedbk ) .AND. ( .NOT. l_casacnp ) )       &
         STOP 'casaCNP required to get prognostic LAI or Vcmax'
    IF( l_vcmaxFeedbk .AND. icycle < 2 )                                     &
         STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
    IF( icycle > 0 .AND. ( .NOT. soilparmnew ) )                             &
         STOP 'casaCNP must use new soil parameters'

    ! Open log file:
    ! MPI: worker logs go to the black hole
    ! by opening the file we don't need to touch any of the code that writes
    ! to it and may be called somewhere by a worker
    ! OPEN(logn,FILE=filename%log)


    !! Check for gswp run
    ! MPI: done by the master only; if check fails then master MPI_Aborts
    ! everyone
    !IF (ncciy /= 0) THEN
    !
    !   PRINT *, 'Looking for global offline run info.'
    !
    !   IF (ncciy < 1986 .OR. ncciy > 1995) THEN
    !      PRINT *, 'Year ', ncciy, ' outside range of dataset!'
    !      STOP 'Please check input in namelist file.'
    !   ELSE
    !
    !      CALL prepareFiles(ncciy)
    !
    !   ENDIF
    !
    !ENDIF

    ! Open met data and get site information from netcdf file.
    ! This retrieves time step size, number of timesteps, starting date,
    ! latitudes, longitudes, number of sites.
    ! MPI: master only; necessary info will be received by MPI below
    !CALL open_met_file( dels, kend, spinup, C%TFRZ )

    ! Checks where parameters and initialisations should be loaded from.
    ! If they can be found in either the met file or restart file, they will
    ! load from there, with the met file taking precedence. Otherwise, they'll
    ! be chosen from a coarse global grid of veg and soil types, based on
    ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
    ! MPI: master only; necessary info will be received by MPI below
    !CALL load_parameters( met, air, ssnow, veg, bgc,                            &
    !                      soil, canopy, rough, rad, sum_flux,                   &
    !                      bal, logn, vegparmnew, casabiome, casapool,           &
    !                      casaflux, casamet, casabal, phen, C%EMSOIL,        &
    !                      C%TFRZ )

    ! Check for leap-year settings
    CALL MPI_Bcast (leaps, 1, MPI_LOGICAL, 0, comm, ierr)

    ktau_tot = 0
    SPINLOOP:DO
       YEAR: DO YYYY= CABLE_USER%YearStart,  CABLE_USER%YearEnd
          CurYear = YYYY

          IF ( leaps .AND. IS_LEAPYEAR( YYYY ) ) THEN
             LOY = 366
          ELSE
             LOY = 365
          ENDIF

          IF ( CALL1 ) THEN

             IF (.NOT.spinup)	spinConv=.TRUE.

             ! MPI: bcast to workers so that they don't need to open the met
             ! file themselves
             CALL MPI_Bcast (dels, 1, MPI_REAL, 0, comm, ierr)
          ENDIF

          ! MPI: receive from master ending time fields
          CALL MPI_Bcast (kend, 1, MPI_INTEGER, 0, comm, ierr)


          IF ( CALL1 ) THEN
             ! MPI: need to know extents before creating datatypes
             CALL find_extents

             ! MPI: receive decomposition info from the master
             CALL worker_decomp(comm)

             ! MPI: in overlap version sends and receives occur on separate comms
             CALL MPI_Comm_dup (comm, icomm, ierr)
             CALL MPI_Comm_dup (comm, ocomm, ierr)

             ! MPI: data set in load_parameter is now received from
             ! the master

             CALL worker_cable_params(comm, met,air,ssnow,veg,bgc,soil,canopy,&
                  &                        rough,rad,sum_flux,bal)

             !mrd561 debug
             WRITE(wlogn,*) ' ssat_vec min',MINVAL(soil%ssat_vec),MINLOC(soil%ssat_vec)
             WRITE(wlogn,*) ' sfc_vec min',MINVAL(soil%sfc_vec),MINLOC(soil%sfc_vec)
             WRITE(wlogn,*) ' wb min',MINVAL(ssnow%wb),MINLOC(ssnow%wb)
             CALL flush(wlogn)

             IF (cable_user%call_climate) THEN
                CALL worker_climate_types(comm, climate, ktauday )
             ENDIF

             ! MPI: mvtype and mstype send out here instead of inside worker_casa_params
             !      so that old CABLE carbon module can use them. (BP May 2013)
             CALL MPI_Bcast (mvtype, 1, MPI_INTEGER, 0, comm, ierr)
             CALL MPI_Bcast (mstype, 1, MPI_INTEGER, 0, comm, ierr)

             ! MPI: casa parameters received only if cnp module is active
             IF (icycle>0) THEN

                CALL worker_casa_params (comm,casabiome,casapool,casaflux,casamet,&
                     &                        casabal,phen)

                ! MPI: POP restart received only if pop module AND casa are active
                IF ( CABLE_USER%CALL_POP ) CALL worker_pop_types (comm,veg,casamet,pop)

             END IF

             ! MPI: create inp_t type to receive input data from the master
             ! at the start of every timestep
             CALL worker_intype (comm,met,veg)

             ! MPI: casa parameters received only if cnp module is active
             ! MPI: create send_t type to send the results to the master
             ! at the end of every timestep
             CALL worker_outtype (comm,met,canopy,ssnow,rad,bal,air,soil,veg)

             ! MPI: casa parameters received only if cnp module is active
             ! MPI: create type to send casa results back to the master
             ! only if cnp module is active
             IF (icycle>0) THEN
                CALL worker_casa_type (comm, casapool,casaflux, &
                     casamet,casabal, phen)

                IF ( CABLE_USER%CASA_DUMP_READ .OR. CABLE_USER%CASA_DUMP_WRITE ) &
                     CALL worker_casa_dump_types(comm, casamet, casaflux, phen, climate)
                WRITE(wlogn,*) 'cable_mpiworker, POPLUC: ',  CABLE_USER%POPLUC
                WRITE(*,*) 'cable_mpiworker, POPLUC: ',  CABLE_USER%POPLUC
                CALL flush(wlogn)
                IF ( CABLE_USER%POPLUC ) &
                     CALL worker_casa_LUC_types( comm, casapool, casabal)


                ! MPI: casa parameters received only if cnp module is active
             END IF

             ! MPI: create type to send restart data back to the master
             ! only if restart file is to be created
             IF(output%restart) THEN

                CALL worker_restart_type (comm, canopy, air)

             END IF

             ! Open output file:
             ! MPI: only the master writes to the files
             !CALL open_output_file( dels, soil, veg, bgc, rough )

             ssnow%otss_0   = ssnow%tgg(:,1)
             ssnow%otss     = ssnow%tgg(:,1)
             canopy%fes_cor = 0.
             canopy%fhs_cor = 0.
             met%ofsd = 0.1

             ! CALL worker_sumcasa_types(comm, sum_casapool, sum_casaflux)


             !count_sum_casa = 0


             IF( icycle>0 .AND. spincasa) THEN
                WRITE(wlogn,*) 'EXT spincasacnp enabled with mloop= ', mloop
                CALL worker_spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                     casaflux,casamet,casabal,phen,POP,climate,LALLOC, icomm, ocomm)
                SPINconv = .FALSE.
                CASAONLY                   = .TRUE.
                ktau_gl = 0
                ktau = 0

             ELSEIF ( casaonly .AND. (.NOT. spincasa) .AND. cable_user%popluc) THEN
                CALL worker_CASAONLY_LUC(dels,kstart,kend,veg,soil,casabiome,casapool, &
                     casaflux,casamet,casabal,phen,POP,climate,LALLOC, &
                     icomm, ocomm)
                SPINconv = .FALSE.
                ktau_gl = 0
                ktau = 0
             ENDIF

          ELSE
             IF (icycle.GT.0) THEN
                ! re-initalise annual flux sums
                casabal%FCgppyear =0.0
                casabal%FCrpyear  =0.0
                casabal%FCnppyear =0.0
                casabal%FCrsyear  =0.0
                casabal%FCneeyear =0.0
             ENDIF

          ENDIF !CALL1


          ! globally (WRT code) accessible kend through USE cable_common_module
          ktau_gl  = 0
          kwidth_gl = INT(dels)
          kend_gl  = kend
          knode_gl = 0

          IF (spincasa .OR. casaonly) THEN
             EXIT
          ENDIF

  if( .NOT. allocated(heat_cap_lower_limit) ) then
    allocate( heat_cap_lower_limit(mp,ms) ) 
    heat_cap_lower_limit = 0.01
  end if

  call spec_init_soil_snow(dels, soil, ssnow, canopy, met, bal, veg, heat_cap_lower_limit)

  
          ! IF (.NOT.spincasa) THEN
          ! time step loop over ktau
          KTAULOOP:DO ktau=kstart, kend

             ! increment total timstep counter
             ktau_tot = ktau_tot + 1

             WRITE(wlogn,*) 'ktau -',ktau_tot
             CALL flush(wlogn)

             ! globally (WRT code) accessible kend through USE cable_common_module
             ktau_gl = ktau_gl + 1

             ! somethings (e.g. CASA-CNP) only need to be done once per day
             ktauday=INT(24.0*3600.0/dels)
!!$             idoy = mod(ktau/ktauday,365)
!!$             IF(idoy==0) idoy=365
!!$
!!$             ! needed for CASA-CNP
!!$             nyear =INT((kend-kstart+1)/(365*ktauday))

             ! some things (e.g. CASA-CNP) only need to be done once per day
             idoy =INT( MOD((REAL(ktau+koffset)/REAL(ktauday)),REAL(LOY)))
             IF ( idoy .EQ. 0 ) idoy = LOY

             ! needed for CASA-CNP
             nyear =INT((kend-kstart+1)/(LOY*ktauday))




             canopy%oldcansto=canopy%cansto

             ! Get met data and LAI, set time variables.
             ! Rainfall input may be augmented for spinup purposes:
             met%ofsd = met%fsd(:,1) + met%fsd(:,2)
             ! MPI: input file read on the master only
             !CALL get_met_data( spinup, spinConv, met, soil,                    &
             !                   rad, veg, kend, dels, C%TFRZ, ktau )

             ! MPI: receive input data for this step from the master
             IF ( .NOT. CASAONLY ) THEN

                CALL MPI_Recv (MPI_BOTTOM, 1, inp_t, 0, ktau_gl, icomm, stat, ierr)

                ! MPI: receive casa_dump_data for this step from the master
             ELSEIF ( IS_CASA_TIME("dread", yyyy, ktau, kstart, koffset, &
                  kend, ktauday, wlogn) ) THEN
                CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, ktau_gl, icomm, stat, ierr)
             END IF

             ! MPI: some fields need explicit init, because we don't transfer
             ! them for better performance
             ! in the serial version this is done in get_met_data
             ! after input has been read from the file
             met%tvair = met%tk
             met%tvrad = met%tk

             ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
             IF (l_vcmaxFeedbk) CALL casa_feedback( ktau, veg, casabiome,    &
                  casapool, casamet )

             IF (l_laiFeedbk) veg%vlai(:) = REAL(casamet%glai(:))

             IF (cable_user%CALL_climate) &
                  CALL cable_climate(ktau_tot,kstart,kend,ktauday,idoy,LOY,met, &
                  climate, canopy, air, rad, dels, mp)


   IF (.NOT. allocated(c1)) ALLOCATE( c1(mp,nrb), rhoch(mp,nrb), xk(mp,nrb) )
             ! CALL land surface scheme for this timestep, all grid points:
   CALL cbm( ktau, dels, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg, climate, xk, c1, rhoch )

             ssnow%smelt  = ssnow%smelt*dels
             ssnow%rnof1  = ssnow%rnof1*dels
             ssnow%rnof2  = ssnow%rnof2*dels
             ssnow%runoff = ssnow%runoff*dels


             !jhan this is insufficient testing. condition for
             !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)

             IF(icycle >0) THEN

                CALL bgcdriver( ktau, kstart, kend, dels, met,          &
                     ssnow, canopy, veg, soil,climate, casabiome,       &
                     casapool, casaflux, casamet, casabal,              &
                     phen, pop, spinConv, spinup, ktauday, idoy, loy,   &
                     .FALSE., .FALSE., LALLOC )

                ! IF(MOD((ktau-kstart+1),ktauday)==0) THEN
                CALL MPI_Send (MPI_BOTTOM,1, casa_t,0,ktau_gl,ocomm,ierr)

                !  ENDIF

                IF ( IS_CASA_TIME("write", yyyy, ktau, kstart, &
                     koffset, kend, ktauday, wlogn) ) THEN
                   ! write(wlogn,*), 'IN IS_CASA', casapool%cplant(:,1)
                   !      CALL MPI_Send (MPI_BOTTOM,1, casa_t,0,ktau_gl,ocomm,ierr)
                ENDIF


                ! MPI: send the results back to the master
                IF( ((.NOT.spinup).OR.(spinup.AND.spinConv)) .AND. &
                     IS_CASA_TIME("dwrit", yyyy, ktau, kstart, &
                     koffset, kend, ktauday, logn))  &
                     CALL MPI_Send (MPI_BOTTOM, 1, casa_dump_t, 0, ktau_gl, ocomm, ierr)

             ENDIF

             ! sumcflux is pulled out of subroutine cbm
             ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
             CALL sumcflux( ktau, kstart, kend, dels, bgc,              &
                  canopy, soil, ssnow, sum_flux, veg,                   &
                  met, casaflux, l_vcmaxFeedbk )

             ! MPI: send the results back to the master
             CALL MPI_Send (MPI_BOTTOM, 1, send_t, 0, ktau_gl, ocomm, ierr)

             ! Write time step's output to file if either: we're not spinning up
             ! or we're spinning up and the spinup has converged:
             ! MPI: writing done only by the master
             !IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
             !   CALL write_output( dels, ktau, met, canopy, ssnow,                    &
             !                      rad, bal, air, soil, veg, C%SBOLTZ, &
             !                      C%EMLEAF, C%EMSOIL )


             CALL1 = .FALSE.

          END DO KTAULOOP ! END Do loop over timestep ktau
          ! ELSE

          CALL1 = .FALSE.
          ! ENDIF


          CALL flush(wlogn)
          IF (icycle >0 .AND.   cable_user%CALL_POP) THEN

             IF (CABLE_USER%POPLUC) THEN

                WRITE(wlogn,*), 'before MPI_Send casa_LUC'
                ! worker sends casa updates required for LUC calculations here
                CALL MPI_Send (MPI_BOTTOM, 1, casa_LUC_t, 0, 0, ocomm, ierr)
                WRITE(wlogn,*), 'after MPI_Send casa_LUC'
                ! master calls LUCDriver here
                ! worker receives casa and POP updates
                CALL MPI_Recv( POP%pop_grid(1), POP%np, pop_t, 0, 0, icomm, stat, ierr )

             ENDIF
             ! one annual time-step of POP
             CALL POPdriver(casaflux,casabal,veg, POP)

             CALL worker_send_pop (POP, ocomm)

             IF (CABLE_USER%POPLUC) &
                  CALL MPI_Recv (MPI_BOTTOM, 1, casa_LUC_t, 0, nyear, icomm, stat, ierr)

          ENDIF















          IF ( ((.NOT.spinup).OR.(spinup.AND.spinConv)).AND. &
               CABLE_USER%CALL_POP) THEN

             !CALL worker_send_pop (POP, ocomm)

          ENDIF


       END DO YEAR


       IF (spincasa .OR. casaonly) THEN
          EXIT
       ENDIF
       !!jhan this is insufficient testing. condition for
       !!spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
       !! see if spinup (if conducting one) has converged:
       !IF(spinup.AND..NOT.spinConv) THEN
       !
       !   ! Write to screen and log file:
       !   WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),      &
       !         ' of data set complete...'
       !   WRITE(wlogn,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),   &
       !         ' of data set complete...'
       !
       !   ! IF not 1st run through whole dataset:
       !   IF( INT( ktau_tot/kend ) > 1 ) THEN
       !
       !      ! evaluate spinup
       !      IF( ANY( ABS(ssnow%wb-soilMtemp)>delsoilM).OR.                     &
       !          ANY(ABS(ssnow%tgg-soilTtemp)>delsoilT) ) THEN
       !
       !         ! No complete convergence yet
       !         PRINT *, 'ssnow%wb : ', ssnow%wb
       !         PRINT *, 'soilMtemp: ', soilMtemp
       !         PRINT *, 'ssnow%tgg: ', ssnow%tgg
       !         PRINT *, 'soilTtemp: ', soilTtemp
       !
       !      ELSE ! spinup has converged
       !
       !         spinConv = .TRUE.
       !         ! Write to screen and log file:
       !         WRITE(*,'(A33)') ' Spinup has converged - final run'
       !         WRITE(logn,'(A52)')                                             &
       !                    ' Spinup has converged - final run - writing all data'
       !         WRITE(logn,'(A37,F8.5,A28)')                                    &
       !                    ' Criteria: Change in soil moisture < ',             &
       !                    delsoilM, ' in any layer over whole run'
       !         WRITE(logn,'(A40,F8.5,A28)' )                                   &
       !                    '           Change in soil temperature < ',          &
       !                    delsoilT, ' in any layer over whole run'
       !      END IF

       !   ELSE ! allocate variables for storage
       !
       !     ALLOCATE( soilMtemp(mp,ms), soilTtemp(mp,ms) )
       !
       !   END IF
       !
       !   ! store soil moisture and temperature
       !   soilTtemp = ssnow%tgg
       !   soilMtemp = REAL(ssnow%wb)

       !ELSE

       !   ! if not spinning up, or spin up has converged, exit:
       !   EXIT
       !
       !END IF

       ! MPI: learn from the master whether it's time to quit
       CALL MPI_Bcast (loop_exit, 1, MPI_LOGICAL, 0, comm, ierr)

       IF (loop_exit) THEN
          EXIT
       END IF

    END DO SPINLOOP

    IF (icycle > 0 .AND. (.NOT.spincasa).AND. (.NOT.casaonly)) THEN

       ! MPI: send casa results back to the master
       CALL MPI_Send (MPI_BOTTOM, 1, casa_t, 0, ktau_gl, ocomm, ierr)

       ! MPI: output file written by master only
       !CALL casa_poolout( ktau, veg, soil, casabiome,                           &
       !                   casapool, casaflux, casamet, casabal, phen )
       ! MPI: output file written by master only
       !CALL casa_fluxout( nyear, veg, soil, casabal, casamet)

    END IF

    ! Write restart file if requested:
    IF(output%restart .AND. (.NOT. CASAONLY)) THEN
       ! MPI: send variables that are required by create_restart
       CALL MPI_Send (MPI_BOTTOM, 1, restart_t, 0, ktau_gl, comm, ierr)
       ! MPI: output file written by master only

       IF (cable_user%CALL_climate) &
            CALL MPI_Send (MPI_BOTTOM, 1, climate_t, 0, ktau_gl, comm, ierr)
    END IF

    ! MPI: cleanup
    CALL worker_end(icycle, output%restart)
    ! MPI: open and close by master only
    ! Close met data input file:
    !CALL close_met_file
    ! MPI: open and close by master only
    ! Close output file and deallocate main variables:
    !CALL close_output_file( bal, air, bgc, canopy, met,                         &
    !                        rad, rough, soil, ssnow,                            &
    !                        sum_flux, veg )

    !WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck

    ! Close log file
    ! MPI: closes handle to /dev/null in workers
    CLOSE(wlogn)


    RETURN

  END SUBROUTINE mpidrv_worker



  ! ============== PRIVATE SUBROUTINES USED ONLY BY THE MPI WORKERS ===============


  ! MPI: receives grid decomposition info from the master
  SUBROUTINE worker_decomp (comm)

    USE mpi

    USE cable_def_types_mod, ONLY: mland, mp

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: comm ! MPI communicator to talk to the workers

    INTEGER :: stat(MPI_STATUS_SIZE), ierr

    ! receive number of landpoints assigned to this worker
    CALL MPI_Recv (mland, 1, MPI_INTEGER, 0, 0, comm, stat, ierr)

    ! receive number of land patches assigned to this worker
    CALL MPI_Recv (mp, 1, MPI_INTEGER, 0, 0, comm, stat, ierr)

    RETURN

  END SUBROUTINE worker_decomp


  ! MPI: creates param_t type for the worker to receive the default parameters
  ! from the master process
  ! then receives the parameters
  ! and finally frees the MPI type
  SUBROUTINE worker_cable_params (comm,met,air,ssnow,veg,bgc,soil,canopy,&
       rough,rad,sum_flux,bal)

    USE mpi

    USE cable_def_types_mod
    USE cable_IO_vars_module
    USE cable_input_module, ONLY: allocate_cable_vars
    USE cable_common_module,  ONLY: calcsoilalbedo

    IMPLICIT NONE

    ! subroutine arguments

    INTEGER, INTENT(IN) :: comm ! MPI communicator

    TYPE (met_type), INTENT(OUT) :: met
    TYPE (air_type), INTENT(OUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssnow
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal


    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: stat(MPI_STATUS_SIZE), ierr
    INTEGER :: landp_t, patch_t, param_t

    INTEGER :: r1len, r2len, i1len, llen ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank,  off, ierr2, rcount, pos

    CHARACTER, DIMENSION(:), ALLOCATABLE :: rbuf

    CALL MPI_Comm_rank (comm, rank, ierr)

    ! mp and mland should have been received previously by
    ! worker_decomp

    ! creates types to receive slices of landpt and patch arrays from the master
    CALL decomp_types (landp_t, patch_t)

    ! Allocate spatial heterogeneity variables:
    ALLOCATE(landpt(mland))

    ! and receive own slice from the master
    CALL MPI_Recv (landpt, mland, landp_t, 0, 0, comm, stat, ierr)

    CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssnow, &
         sum_flux,veg,mp)

    ! receive slice of patch array that was allocated above inside
    ! allocate_cable_vars
    CALL MPI_Recv (patch, mp, patch_t, 0, 0, comm, stat, ierr)

    ! MPI: TODO: probably not a bad idea to free landp_t and patch_t types

    ntyp = nparam

    ! ntyp increases if include ... Ticket #27
    IF (calcsoilalbedo) THEN
       ntyp = nparam + 1
    END IF

    ALLOCATE (blen(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

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

    bidx = bidx + 1
    CALL MPI_Get_address (met%ca, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%year, displs(bidx), ierr)
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%moy, displs(bidx), ierr)
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%doy, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%hod, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%fsd, displs(bidx), ierr)
    blen(bidx) = swb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%fld, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%precip, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%precip_sn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%tk, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%tvair, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%tvrad, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%pmb, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%ua, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%qv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%qvair, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%da, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%dva, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%coszen, displs(bidx), ierr)
    blen(bidx) = r1len


    ! ----------- air --------------

    bidx = bidx + 1
    CALL MPI_Get_address (air%rho, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%volm, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%rlam, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%qsat, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%epsi, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%visc, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%psyc, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%dsatdk, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (air%cmolar, displs(bidx), ierr)
    blen(bidx) = r1len


    ! ----------- ssnow --------------

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dtmlt, displs(bidx), ierr)
    blen(bidx) = 3 * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%pudsto, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%pudsmx, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%albsoilsn, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%cls, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfn_dtg, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfh_dtg, displs(bidx), ierr)
    blen(bidx) = r1len

    !INH - REV_CORR
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfe_dtg, displs(bidx), ierr)
    blen(bidx) = r1len

    !INH - REV_CORR - no longer needed
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfe_ddq, displs(bidx), ierr)
    blen(bidx) = r1len

    !INH - REV_CORR - no longer needed
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ddq_dtg, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%evapsn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop1, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop2, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop3, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%gammzz, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%isflag, displs(bidx), ierr)
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%osnowd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%potev, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%pwb_min, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%runoff, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%rnof1, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%rnof2, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%rtsoil, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%sconds, displs(bidx), ierr)
    blen(bidx) = msn * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%sdepth, displs(bidx), ierr)
    blen(bidx) = msn * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%smass, displs(bidx), ierr)
    blen(bidx) = msn * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%snage, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%snowd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%smelt, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ssdn, displs(bidx), ierr)
    blen(bidx) = msn * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ssdnn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tgg, displs(bidx), ierr)
    blen(bidx) = ms * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tggsn, displs(bidx), ierr)
    blen(bidx) = msn * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tss, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wb, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbfice, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbice, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wblf, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    ! additional  for sli
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%S, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%Tsoil, displs(bidx), ierr)
    blen(bidx) = ms * r2len
!!$
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%thetai, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%snowliq, displs(bidx), ierr)
    blen(bidx) = 3 * r2len


    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%Tsurface, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%h0, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%nsnow, displs(bidx), ierr)
    blen(bidx) = I1len
    ! end additional for sli



    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbtot, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wb_lake, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%sinfil, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%evapfbl, displs(bidx), ierr)
    blen(bidx) = ms * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%qstss, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wetfac, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%owetfac, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%t_snwlr, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tggav, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%otss, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%otss_0, displs(bidx), ierr)
    blen(bidx) = r1len

    ! ----------- veg --------------

    bidx = bidx + 1
    CALL MPI_Get_address (veg%canst1, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%dleaf, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%ejmax, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%frac4, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%froot, displs(bidx), ierr)
    blen(bidx) = ms * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%hc, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%iveg, displs(bidx), ierr)
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%meth, displs(bidx), ierr)
    ! Maciej: veg%meth is REAL
    !    blen(bidx) = i1len
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%rp20, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%rpcoef, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%shelrb, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%wai, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%vegcf, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%tminvj, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%tmaxvj, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%vbeta, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%xalbnir, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%vcmax, displs(bidx), ierr)
    blen(bidx) = r1len

    !  bidx = bidx + 1
    !  CALL MPI_Get_address (veg%vlai, displs(bidx), ierr)
    !  blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%xfang, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%extkn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%deciduous, displs(bidx), ierr)
    ! Maciej: deciduous is logical
    !    blen(bidx) = r1len
    blen(bidx) = llen

    bidx = bidx + 1
    CALL MPI_Get_address (veg%a1gs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%d0gs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%alpha, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%convex, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%cfrd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%gswmin, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%conkc0, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%conko0, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%ekc, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%eko, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%clitt, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%zr, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%gamma, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%refl, displs(bidx), ierr)
    blen(bidx) = 2 * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%taul, displs(bidx), ierr)
    blen(bidx) = 2 * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%disturbance_interval, displs(bidx), ierr)
    blen(bidx) = 2 * i1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%disturbance_intensity, displs(bidx), ierr)
    ! Maciej: disturbance_intensity is REAL(r_2)
    !    blen(bidx) = 2 * r1len
    blen(bidx) = 2 * r2len

    ! Ticket #56, adding new veg parms
    bidx = bidx + 1
    CALL MPI_Get_address (veg%g0, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%g1, displs(bidx), ierr)
    blen(bidx) = r1len
    ! Ticket #56, finish adding new veg parms

    ! ----------- bgc --------------

    bidx = bidx + 1
    CALL MPI_Get_address (bgc%cplant, displs(bidx), ierr)
    blen(bidx) = ncp * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bgc%csoil, displs(bidx), ierr)
    blen(bidx) = ncs * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bgc%ratecp, displs(bidx), ierr)
    blen(bidx) = ncp * extr1

    bidx = bidx + 1
    CALL MPI_Get_address (bgc%ratecs, displs(bidx), ierr)
    blen(bidx) = ncs * extr1

    ! ----------- soil --------------

    bidx = bidx + 1
    CALL MPI_Get_address (soil%albsoil, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%bch, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%c3, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%clay, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%cnsd, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%css, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%hsbh, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%hyds, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%i2bp3, displs(bidx), ierr)
    ! Maciej: i2bp3 is REAL
    !    blen(bidx) = i1len
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%ibp2, displs(bidx), ierr)
    ! Maciej: ibp2 is REAL
    !    blen(bidx) = i1len
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%isoilm, displs(bidx), ierr)
    ! Maciej isoilm is INTEGER
    !    blen(bidx) = r1len
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%rhosoil, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%rs20, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%sand, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%sfc, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%silt, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%ssat, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%sucs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%swilt, displs(bidx), ierr)
    blen(bidx) = r1len

    ! extra for sli
    bidx = bidx + 1
    CALL MPI_Get_address (soil%zeta, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%fsatmax, displs(bidx), ierr)
    blen(bidx) = r2len
    ! end extra for sil
    bidx = bidx + 1
    CALL MPI_Get_address (soil%zse, displs(bidx), ierr)
    blen(bidx) = ms * extr1

    bidx = bidx + 1
    CALL MPI_Get_address (soil%zshh, displs(bidx), ierr)
    blen(bidx) = (ms + 1) * extr1

    ! pass soilcolour albedo as well if including Ticket #27
    IF (calcsoilalbedo) THEN
       bidx = bidx + 1
       CALL MPI_Get_address (soil%soilcol, displs(bidx), ierr)
       blen(bidx) = r1len
    END IF

    ! ----------- canopy --------------

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fess, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fesp, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%cansto, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%oldcansto, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%cduv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%delwc, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%dewmm, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%dgdtg, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fe, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fh, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fpn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frpw, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frpr, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fnee, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frday, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fnv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fev, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fevc, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fevw, displs(bidx), ierr)
    blen(bidx) = r1len

    !  bidx = bidx + 1
    !  CALL MPI_Get_address (canopy%potev_c, displs(bidx), ierr)
    !  blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhvw, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fns, displs(bidx), ierr)
    blen(bidx) = r1len

    !INH - REV_CORR - temporary?
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fns_cor, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fes, displs(bidx), ierr)
    blen(bidx) = r2len

    !INH - REV_CORR - temporary?
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fes_cor, displs(bidx), ierr)
    blen(bidx) = r2len

    !INH - SSEB - temporary?
    !bidx = bidx + 1
    !CALL MPI_Get_address (canopy%fescor_upp, displs(bidx), ierr)
    !blen(bidx) = r2len

    !INH - SSEB - temporary?
    !bidx = bidx + 1
    !CALL MPI_Get_address (canopy%fescor_low, displs(bidx), ierr)
    !blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhs, displs(bidx), ierr)
    blen(bidx) = r1len

    !INH - REV_CORR - temporary?
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhs_cor, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fwet, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%ga, displs(bidx), ierr)
    blen(bidx) = r1len

    !INH - REV_CORR - temporary?
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%ga_cor, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%ghflux, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%precis, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%qscrn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%rnet, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%segg, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%sghflux, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%spill, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%through, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%tscrn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%tv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%us, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%uscrn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%vlaiw, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%rghlai, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%wcint, displs(bidx), ierr)
    blen(bidx) = r1len



    !  bidx = bidx + 1
    !  CALL MPI_Get_address (canopy%rwater, displs(bidx), ierr)
    !  blen(bidx) = ms * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%evapfbl, displs(bidx), ierr)
    ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
    blen(bidx) = ms * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%epot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fnpp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fevw_pot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%gswx_T, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%cdtq, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%wetfac_cs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fwsoil, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%gswx, displs(bidx), ierr)
    blen(bidx) = mf * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%zetar, displs(bidx), ierr)
    blen(bidx) = niter * r1len

    ! ------- rough -------

    bidx = bidx + 1
    CALL MPI_Get_address (rough%coexp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%disp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%hruff, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%hruff_grmx, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%rt0us, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%rt1usa, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%rt1usb, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%rt1, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%term2, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%term3, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%term5, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%term6, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%usuh, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%za_uv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%za_tq, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%z0m, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%zref_uv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%zref_tq, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%zruffs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%z0soilsn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rough%z0soil, displs(bidx), ierr)
    blen(bidx) = r1len

    ! --------rad --------

    bidx = bidx + 1
    CALL MPI_Get_address (rad%albedo, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkb, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkd2, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%flws, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%fvlai, displs(bidx), ierr)
    blen(bidx) = mf * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%gradis, displs(bidx), ierr)
    blen(bidx) = mf * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%latitude, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%lwabv, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%qcan, displs(bidx), ierr)
    blen(bidx) = mf * nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%qssabs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%rhocdf, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%rniso, displs(bidx), ierr)
    blen(bidx) = mf * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%scalex, displs(bidx), ierr)
    blen(bidx) = mf * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%transd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%trad, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%reffdf, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%reffbm, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkbm, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkdm, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%fbeam, displs(bidx), ierr)
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%cexpkbm, displs(bidx), ierr)
    ! Maciej: cexpkbm is mp*swb
    !    blen(bidx) = nrb * r1len
    blen(bidx) = swb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%cexpkdm, displs(bidx), ierr)
    ! Maciej: cexpkdm is mp*swb
    !    blen(bidx) = nrb * r1len
    blen(bidx) = swb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%rhocbm, displs(bidx), ierr)
    ! Maciej: rhocbm is mp*nrb
    !    blen(bidx) = swb * r1len
    blen(bidx) = nrb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%transb, displs(bidx), ierr)
    blen(bidx) = r1len

    ! ------- sum_flux -----

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumpn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumrp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumrpw, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumrpr, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumrs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumrd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%dsumpn, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%dsumrp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%dsumrs, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%dsumrd, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumxrp, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (sum_flux%sumxrs, displs(bidx), ierr)
    blen(bidx) = r1len

    ! ------- bal ----

    bidx = bidx + 1
    CALL MPI_Get_address (bal%drybal, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%ebal, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%ebal_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%ebal_cncheck, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%ebal_tot_cncheck, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%evap_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%osnowd0, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%precip_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%rnoff_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%wbal, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%wbal_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%wbtot0, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%wetbal, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%owbtot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%evapc_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%evaps_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%rnof1_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%rnof2_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%snowdc_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%wbal_tot1, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%delwc_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%qasrf_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%qfsrf_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%qssrf_tot, displs(bidx), ierr)
    blen(bidx) = r1len

    ! additional field missing from previous versions;
    ! added when trying to fix a bug in the new mpi code
    ! the order of these new fields follows the order of their
    ! declaration in cable_define_types.F90

    bidx = bidx + 1
    CALL MPI_Get_address (bal%ebaltr, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%ebal_tottr, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (bal%cansto0, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%iantrct, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tss_p, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%deltss, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%owb1, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbtot1, displs(bidx), ierr)
    blen(bidx) = r1len

    !    Maciej: duplicate!
    !    bidx = bidx + 1
    !    CALL MPI_Get_address (ssnow%wbtot1, displs(bidx), ierr)
    !    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tprecip, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tevap, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%trnoff, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%totenbal, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%totenbal2, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fland, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ifland, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tilefrac, displs(bidx), ierr)
    blen(bidx) = n_tiles * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%qasrf, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%qfsrf, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%qssrf, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (veg%vlaimax, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%albedo_T, displs(bidx), ierr)
    blen(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%longitude, displs(bidx), ierr)
    blen(bidx) = r1len

    !mrd add new GW parameters here
    !2D
    bidx = bidx + 1
    CALL MPI_Get_address (soil%ssat_vec, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%sucs_vec, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%hyds_vec, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%bch_vec, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%watr, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%swilt_vec, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%sfc_vec, displs(bidx), ierr)
    blen(bidx) = ms * r2len


    !1d
    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWssat_vec, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWsucs_vec, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWhyds_vec, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWbch_vec, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWwatr, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWz, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%GWdz, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%slope, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (soil%slope_std, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%GWwb, displs(bidx), ierr)
    blen(bidx) = r2len

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,' invalid number of param_t fields',bidx,', fix it!'
       CALL MPI_Abort (comm, 1, ierr)
    END IF

    CALL MPI_Type_create_struct (bidx, blen, displs, types, param_t, ierr)
    CALL MPI_Type_commit (param_t, ierr)

    CALL MPI_Type_size (param_t, tsize, ierr)
    CALL MPI_Type_get_extent (param_t, tmplb, text, ierr)

    WRITE (*,*) 'worker param_t blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...

    ! so, now receive all the parameters
    !    CALL MPI_Recv (MPI_BOTTOM, 1, param_t, 0, 0, comm, stat, ierr)
    !   Maciej: buffered recv + unpac version
    ALLOCATE (rbuf(tsize))
    CALL MPI_Recv (rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    CALL MPI_Get_count (stat, param_t, rcount, ierr2)

    IF (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       CALL MPI_Unpack (rbuf, tsize, pos, MPI_BOTTOM, rcount, param_t, &
            comm, ierr)
       IF (ierr /= MPI_SUCCESS) WRITE(*,*),'cable param unpack error, rank: ',rank,ierr
    ELSE
       WRITE(*,*),'cable param recv rank err err2 rcount: ',rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)


    ! finally free the MPI type
    CALL MPI_Type_Free (param_t, ierr)


    ! all CABLE parameters have been received from the master by now
    RETURN

  END SUBROUTINE worker_cable_params


  ! MPI: creates param_t type for the worker to receive the default casa
  ! parameters from the master process
  ! then receives them
  ! and finally frees the MPI type
  SUBROUTINE worker_casa_params (comm,casabiome,casapool,casaflux,casamet,&
       casabal,phen)

    USE mpi

    USE cable_def_types_mod

    USE casavariable
    USE phenvariable

    IMPLICIT NONE

    ! sub arguments
    INTEGER, INTENT(IN) :: comm  ! MPI communicator

    ! TODO: have these variables been already allocated?
    TYPE (casa_biome)   , INTENT(OUT) :: casabiome
    TYPE (casa_pool)   , INTENT(OUT) :: casapool
    TYPE (casa_flux)   , INTENT(OUT) :: casaflux
    TYPE (casa_met)    , INTENT(OUT) :: casamet
    TYPE (casa_balance), INTENT(OUT) :: casabal
    TYPE (phen_variable), INTENT(OUT)  :: phen

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: stat(MPI_STATUS_SIZE), ierr
    ! INTEGER :: landp_t, patch_t, param_t
    INTEGER :: casa_t

    INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank, off, ierr2, rcount, pos

    CHARACTER, DIMENSION(:), ALLOCATABLE :: rbuf

    off = 1

    CALL MPI_Comm_rank (comm, rank, ierr)

    IF (.NOT. ASSOCIATED (casabiome%ivt2)) THEN
       WRITE (*,*) 'worker alloc casa and phen var with m patches: ',rank,mp
       CALL alloc_casavariable (casabiome, casapool, &
            &      casaflux, casamet, casabal, mp)
       CALL alloc_phenvariable (phen, mp)
    END IF

    ntyp = ncasaparam

    ALLOCATE (blen(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    I1LEN = mp * extid
    llen = mp * extl

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

    !===================================================================
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

    !===================================================================
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

    ! ------ casapool ----

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Clabile, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dClabiledt, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Cplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Nplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Pplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dCplantdt, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dNplantdt, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dPplantdt, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioPCplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Nsoilmin, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Psoillab, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Psoilsorb, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Psoilocc, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dNsoilmindt, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dPsoillabdt, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dPsoilsorbdt, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dPsoiloccdt, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Clitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Nlitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Plitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dClitterdt, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dNlitterdt, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dPlitterdt, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNClitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioPClitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Csoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Nsoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Psoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dCsoildt, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dNsoildt, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%dPsoildt, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCsoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioPCsoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCsoilnew, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCsoilmin, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCsoilmax, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    ! added by LN
    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNPplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNPlitter, displs(bidx), ierr)
    !    blen(bidx) = mplant * r2len
    ! Maciej
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNPsoil, displs(bidx), ierr)
    !    blen(bidx) = mplant * r2len
    ! Maciej
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%cwoodprod, displs(bidx), ierr)
    blen(bidx) = mwood * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nwoodprod, displs(bidx), ierr)
    blen(bidx) = mwood * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%pwoodprod, displs(bidx), ierr)
    blen(bidx) = mwood * r2len

    ! ------- casaflux ----

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Cgpp, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Cnpp, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Crp, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Crgplant, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nminfix, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nminuptake, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Plabuptake, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Clabloss, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracClabile, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracCalloc, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracNalloc, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracPalloc, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Crmplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%kplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    ! 3D
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromPtoL, displs(bidx), ierr)
    blen(bidx) = mplant * mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Cnep, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Crsoil, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nmindep, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nminloss, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nminleach, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nupland, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nlittermin, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nsmin, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nsimm, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nsnet, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fNminloss, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fNminleach, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Pdep, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Pwea, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Pleach, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Ploss, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Pupland, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Plittermin, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Psmin, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Psimm, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Psnet, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fPleach, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%kplab, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%kpsorb, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%kpocc, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%kmlabP, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Psorbmax, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%frac_sapwood, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%sapwood_area, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%klitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%ksoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    ! 3D
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromLtoS, displs(bidx), ierr)
    blen(bidx) = msoil * mlitter * r2len

    ! 3D
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromStoS, displs(bidx), ierr)
    blen(bidx) = msoil * msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromLtoCO2, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromStoCO2, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxCtolitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxNtolitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxPtolitter, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxCtosoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxNtosoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxPtosoil, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxCtoCO2, displs(bidx), ierr)
    blen(bidx) = r2len

    ! ------- casamet ----

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%glai, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%Tairk, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%precip, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%tsoilavg, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%moistavg, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%btran, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%lnonwood, displs(bidx), ierr)
    blen(bidx) = I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%Tsoil, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%moist, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%iveg2, displs(bidx), ierr)
    blen(bidx) = I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%ijgcm, displs(bidx), ierr)
    ! Maciej: ijgcm is INTEGER
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%isorder, displs(bidx), ierr)
    ! Maciej: isorder is INTEGER
    blen(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%lat, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%lon, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%areacell, displs(bidx), ierr)
    blen(bidx) = r2len

    ! ------- casabal ----

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCgppyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCnppyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCrmleafyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCrmwoodyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCrmrootyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCrgrowyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCrpyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCrsyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCneeyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FNdepyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FNfixyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FNsnetyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FNupyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FNleachyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FNlossyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FPweayear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FPdustyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FPsnetyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FPupyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FPleachyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FPlossyear, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%glaimon, displs(bidx), ierr)
    blen(bidx) = 12 * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%glaimonx, displs(bidx), ierr)
    blen(bidx) = 12 * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%cplantlast, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%nplantlast, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%pplantlast, displs(bidx), ierr)
    blen(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%clitterlast, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%nlitterlast, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%plitterlast, displs(bidx), ierr)
    blen(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%csoillast, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%nsoillast, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%psoillast, displs(bidx), ierr)
    blen(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%nsoilminlast, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%psoillablast, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%psoilsorblast, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%psoilocclast, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%cbalance, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%nbalance, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%pbalance, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%sumcbal, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%sumnbal, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%sumpbal, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%clabilelast, displs(bidx), ierr)
    blen(bidx) = r2len

    ! ------- phen -------

    bidx = bidx + 1
    CALL MPI_Get_address (phen%phase, displs(bidx), ierr)
    blen(bidx) = I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (phen%TKshed, displs(bidx), ierr)
    blen(bidx) = mvtype * extr2

    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphase, displs(bidx), ierr)
    blen(bidx) = mphase * I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (phen%phen(off), displs(bidx), ierr)
    blen(bidx) =  r1len

    bidx = bidx + 1
    CALL MPI_Get_address (phen%aphen(off), displs(bidx), ierr)
    blen(bidx) = r1len


    bidx = bidx + 1
    CALL MPI_Get_address (phen%phasespin(off,1), displs(bidx), ierr)
    blen(bidx) = mdyear*I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_1(off,1), displs(bidx), ierr)
    blen(bidx) =  mdyear*I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_2(off,1), displs(bidx), ierr)
    blen(bidx) =  mdyear*I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_3(off,1), displs(bidx), ierr)
    blen(bidx) =  mdyear*I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_4(off,1), displs(bidx), ierr)
    blen(bidx) =  mdyear*I1LEN


    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,' invalid number of casa_t param fields ',bidx,', fix it!'
       CALL MPI_Abort (comm, 1, ierr)
    END IF

    CALL MPI_Type_create_struct (bidx, blen, displs, types, casa_t, ierr)
    CALL MPI_Type_commit (casa_t, ierr)

    CALL MPI_Type_size (casa_t, tsize, ierr)
    CALL MPI_Type_get_extent (casa_t, tmplb, text, ierr)

    WRITE (*,*) 'worker casa_t param blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

    ! if anything went wrong the master will mpi_abort
    ! which mpi_recv below is going to catch...

    CALL MPI_Barrier (comm, ierr)

    ! so, now receive all the parameters
    !    CALL MPI_Recv (MPI_BOTTOM, 1, casa_t, 0, 0, comm, stat, ierr)

    !   Maciej: buffered recv + unpac version
    ALLOCATE (rbuf(tsize))
    CALL MPI_Recv (rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    CALL MPI_Get_count (stat, casa_t, rcount, ierr2)

    IF (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       CALL MPI_Unpack (rbuf, tsize, pos, MPI_BOTTOM, rcount, casa_t, &
            comm, ierr)
       IF (ierr /= MPI_SUCCESS) WRITE(*,*),'casa params unpack error, rank: ',rank,ierr
    ELSE
       WRITE(*,*),'casa params recv rank err err2 rcount: ',rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)


    ! finally free the MPI type
    CALL MPI_Type_Free (casa_t, ierr)

    ! all casa parameters have been received from the master by now

    RETURN

  END SUBROUTINE worker_casa_params


  ! MPI: creates inp_t type to receive input data from the master
  SUBROUTINE worker_intype (comm,met,veg)

    USE mpi

    USE cable_def_types_mod

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

    ! met fields
    ! TODO: don't add optional field when not required
    ! (exists% flags)
    bidx = 0

    bidx = bidx + 1
    CALL MPI_Get_address (met%fsd, displs(bidx), ierr)
    blocks(bidx) = swb * r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%tk, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%pmb, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%qv, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%ua, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%precip, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%precip_sn, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%fld, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%ca, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%coszen, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%Ndep, displs(bidx), ierr)
    blocks(bidx) = r1len


    ! veg fields

    bidx = bidx + 1
    CALL MPI_Get_address (veg%vlai, displs(bidx), ierr)
    blocks(bidx) = r1len

    ! additional field missing from previous versions;
    ! added when trying to fix a bug in the new mpi code
    ! the order of these new fields follows the order of their
    ! declaration in cable_define_types.F90

    bidx = bidx + 1
    CALL MPI_Get_address (met%year, displs(bidx), ierr)
    blocks(bidx) = I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (met%moy, displs(bidx), ierr)
    blocks(bidx) = I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (met%doy, displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (met%hod, displs(bidx), ierr)
    blocks(bidx) = r1len


    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,': invalid intype nmat, nvec or n3d constant, fix it!'
       WRITE(*,*) 'bidx: ', bidx
       WRITE(*,*) 'nvec: ', nvec
       WRITE(*,*) 'n3d:', n3d
       WRITE(*,*) 'ntyp', ntyp
       CALL MPI_Abort (comm, 1, ierr)
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
  SUBROUTINE worker_outtype (comm,met,canopy,ssnow,rad,bal,air,soil,veg)

    USE mpi

    USE cable_def_types_mod

    IMPLICIT NONE

    INTEGER :: comm ! MPI communicator to talk to the workers

    TYPE(met_type), INTENT(IN) :: met
    TYPE(canopy_type), INTENT(IN) :: canopy
    TYPE(soil_snow_type), INTENT(IN) :: ssnow
    TYPE(radiation_type), INTENT(IN) :: rad
    TYPE (balances_type),INTENT(INOUT):: bal
    TYPE (air_type),INTENT(IN)     :: air
    TYPE (soil_parameter_type),INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type),INTENT(IN) :: veg ! vegetation parameters

    ! MPI: temp arrays for marshalling all types into a struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types
    INTEGER :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    INTEGER :: r1len, r2len, I1LEN, llen

    INTEGER :: rank, off, cnt
    INTEGER :: bidx, midx, vidx, ierr

    INTEGER :: tsize
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

    ! base index to make types indexing easier
    INTEGER :: istart

    INTEGER :: i

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
    llen = cnt * extl

    ! ------------- 3D arrays -------------

    ! rad 3D
    ! TODO: REAL(r_1) : rad%qcan(landunits,mf,nrb)
    ! CALL MPI_Type_create_hvector (mf*nrb, r1len, r1stride, MPI_BYTE, &
    !  &            m3d_t(1, rank), ierr)
    bidx = 1
    CALL MPI_Get_address (rad%qcan(off,1,1), displs(bidx), ierr)
    blocks(bidx) = r1len * mf * nrb

    ! ------------- 2D arrays -------------

    ! MPI: an hvector type for each vector, maddr contains displacements
    ! for bundling these hvectors into the struct later
    ! block length is number of patches/worker * type extent
    ! stride is global number of patches * type extent
    ! repeat/no of blocks is the 2nd rank

    ! met 2D
    bidx = bidx + 1
    CALL MPI_Get_address (met%fsd(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * swb

    ! canopy 2D
    !midx = 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%rwater(off,1), maddr(midx), ierr) ! 1
    !CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)

    ! TODO: skip, used for restart but not output
    !  bidx = bidx + 1
    !  CALL MPI_Get_address (canopy%rwater(off,1), displs(bidx), ierr)
    !  blocks(bidx) = r1len * ms

    ! midx = midx + 1
    ! REAL(r_2)
    ! CALL MPI_Get_address (canopy%evapfbl(off,1), maddr(midx), ierr) ! 2
    !CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%evapfbl(off,1), displs(bidx), ierr)
    ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
    blocks(bidx) = r1len * ms

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%gswx(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * mf

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%zetar(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * niter

    ! ssnow 2D
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dtmlt(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * 3

    !midx = midx + 1
    ! REAL(r_1)
    ! CALL MPI_Get_address (ssnow%albsoilsn(off,1), maddr(midx), ierr) ! 3
    ! CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    !midx = midx + 1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%albsoilsn(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    ! REAL(r_2)
    !CALL MPI_Get_address (ssnow%gammzz(off,1), maddr(midx), ierr) ! 4
    !CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%gammzz(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%sconds(off,1), maddr(midx), ierr) ! 5
    !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%sconds(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * msn

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%sdepth(off,1), maddr(midx), ierr) ! 6
    !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%sdepth(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * msn

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%smass(off,1), maddr(midx), ierr) ! 7
    !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%smass(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * msn

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%dtmlt(off,1), maddr(midx), ierr) ! 8
    !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    ! MPI: r1134 does not know about this field, comment out
    !bidx = bidx + 1
    !CALL MPI_Get_address (ssnow%dtmlt(off,1), displs(bidx), ierr)
    !blocks(bidx) = r1len * msn

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%ssdn(off,1), maddr(midx), ierr) ! 9
    !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ssdn(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * msn

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%tgg(off,1), maddr(midx), ierr) ! 10
    !CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tgg(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * ms

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%tggsn(off,1), maddr(midx), ierr) ! 11
    !CALL MPI_Type_create_hvector (msn, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tggsn(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * msn

    !midx = midx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (ssnow%wb(off,1), maddr(midx), ierr) ! 12
    !CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wb(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%evapfbl(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * ms

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%wbfice(off,1), maddr(midx), ierr) ! 13
    !CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbfice(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    !midx = midx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (ssnow%wbice(off,1), maddr(midx), ierr) ! 14
    !CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbice(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    !midx = midx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (ssnow%wblf(off,1), maddr(midx), ierr) ! 15
    !CALL MPI_Type_create_hvector (ms, r2len, r2stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wblf(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms
    ! additional  for sli

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%S(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%Tsoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%thetai(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * ms

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%snowliq(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * 3

    ! end additional for sli

    ! rad 2D
    bidx = bidx + 1
    CALL MPI_Get_address (rad%fbeam(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%albedo(off,1), maddr(midx), ierr) ! 16
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%albedo(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%fvlai(off,1), maddr(midx), ierr) ! 17
    !CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%fvlai(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * mf

    !midx = midx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (rad%gradis(off,1), maddr(midx), ierr) ! 18
    !CALL MPI_Type_create_hvector (mf, r2len, r2stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%gradis(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * mf

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%rhocdf(off,1), maddr(midx), ierr) ! 19
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%rhocdf(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%rniso(off,1), maddr(midx), ierr) ! 20
    !CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%rniso(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * mf

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%scalex(off,1), maddr(midx), ierr) ! 21
    !CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%scalex(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * mf

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%reffdf(off,1), maddr(midx), ierr) ! 22
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%reffdf(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%reffbm(off,1), maddr(midx), ierr) ! 23
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%reffbm(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%extkbm(off,1), maddr(midx), ierr) ! 24
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkbm(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%extkdm(off,1), maddr(midx), ierr) ! 25
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkdm(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%cexpkbm(off,1), maddr(midx), ierr) ! 26
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%cexpkbm(off,1), displs(bidx), ierr)
    ! Maciej: cexpkbm is mp*swb
    !    blocks(bidx) = r1len * nrb
    blocks(bidx) = r1len * swb

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%cexpkdm(off,1), maddr(midx), ierr) ! 27
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (rad%cexpkdm(off,1), displs(bidx), ierr)
    ! Maciej: cexpkdm is mp*swb
    !    blocks(bidx) = r1len * nrb
    blocks(bidx) = r1len * swb

    bidx = bidx + 1
    CALL MPI_Get_address (rad%rhocbm(off,1), displs(bidx), ierr)
    ! Maciej: rhocbm is mp*nrb
    !    blocks(bidx) = r1len * swb
    blocks(bidx) = r1len * nrb


    ! air 2D - all fields 1D - skipped

    ! soil 2D
    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%albsoil(off,1), maddr(midx), ierr) ! 28
    !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (soil%albsoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * nrb


    ! veg 2D
    bidx = bidx + 1
    CALL MPI_Get_address (veg%refl(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * 2

    bidx = bidx + 1
    CALL MPI_Get_address (veg%taul(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * 2

    !midx = midx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%froot(off,1), maddr(midx), ierr) ! 29
    !CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
    !  &            mat_t(midx, rank), ierr)
    bidx = bidx + 1
    CALL MPI_Get_address (veg%froot(off,1), displs(bidx), ierr)
    blocks(bidx) = r1len * ms



    ! ------------- 1D arrays -------------

    ! met
    !vidx = 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%ca(off), vaddr(vidx), ierr) ! 1
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%ca(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (met%year(off), vaddr(vidx), ierr) ! 2
    !blen(vidx) = cnt * extid
    ! gol124: not output, removed
    !bidx = bidx + 1
    !CALL MPI_Get_address (met%year(off), displs(bidx), ierr)
    !blocks(bidx) = I1LEN

    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (met%moy(off), vaddr(vidx), ierr) ! 3
    !blen(vidx) = cnt * extid
    ! gol124: not output, removed
    !bidx = bidx + 1
    !CALL MPI_Get_address (met%moy(off), displs(bidx), ierr)
    !blocks(bidx) = I1LEN

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%doy(off), vaddr(vidx), ierr) ! 4
    !blen(vidx) = cnt * extr1
    ! gol124: not output, removed
    !bidx = bidx + 1
    !CALL MPI_Get_address (met%doy(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%hod(off), vaddr(vidx), ierr) ! 5
    !blen(vidx) = cnt * extr1
    ! gol124: not output, removed
    !bidx = bidx + 1
    !CALL MPI_Get_address (met%hod(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%fsd(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    ! MPI: gol124: changed to 2D and move up when Bernard
    ! ported to CABLE_r491
    !bidx = bidx + 1
    !CALL MPI_Get_address (met%fsd(off,1), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%fld(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%fld(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%precip(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%precip(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%precip_sn(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%precip_sn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%tk(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%tk(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%tvair(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%tvair(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%tvrad(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%tvrad(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%pmb(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%pmb(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%ua(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%ua(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%qv(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%qv(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%qvair(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%qvair(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%da(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%da(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%dva(off), vaddr(vidx), ierr)
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%dva(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (met%coszen(off), vaddr(vidx), ierr) ! 19
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (met%coszen(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    ! canopy
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fess(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fesp(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%cansto(off), vaddr(vidx), ierr) ! 20
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%cansto(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%cduv(off), vaddr(vidx), ierr) ! 21
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%cduv(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%delwc(off), vaddr(vidx), ierr) ! 22
    !blen(vidx) = cnt * extr1

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%delwc(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%dewmm(off), vaddr(vidx), ierr) ! 23
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%dewmm(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (canopy%dgdtg(off), vaddr(vidx), ierr) ! 24
    !blen(vidx) = cnt * extr2

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%dgdtg(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fe(off), vaddr(vidx), ierr) ! 25
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fe(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fh(off), vaddr(vidx), ierr) ! 26
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fh(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fpn(off), vaddr(vidx), ierr) ! 27
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fpn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%frp(off), vaddr(vidx), ierr) ! 28
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frp(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%frpw(off), vaddr(vidx), ierr) ! 29
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frpw(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%frpr(off), vaddr(vidx), ierr) ! 30
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frpr(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%frs(off), vaddr(vidx), ierr) ! 31
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frs(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fnee(off), vaddr(vidx), ierr) ! 32
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fnee(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%frday(off), vaddr(vidx), ierr) ! 33
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%frday(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fnv(off), vaddr(vidx), ierr) ! 34
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fnv(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    ! gol124: MPI: DONE until here!!!

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fev(off), vaddr(vidx), ierr) ! 35
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fev(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (canopy%fevc(off), vaddr(vidx), ierr) ! 36
    !blen(vidx) = cnt * extr2
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fevc(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (canopy%fevw(off), vaddr(vidx), ierr) ! 37
    !blen(vidx) = cnt * extr2
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fevw(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (canopy%potev_c(off), vaddr(vidx), ierr) ! 38
    !blen(vidx) = cnt * extr2
    !  bidx = bidx + 1
    !  CALL MPI_Get_address (canopy%potev_c(off), displs(bidx), ierr)
    !  blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fhv(off), vaddr(vidx), ierr) ! 39
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhv(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (canopy%fhvw(off), vaddr(vidx), ierr) ! 40
    !blen(vidx) = cnt * extr2
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhvw(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fns(off), vaddr(vidx), ierr) ! 41
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fns(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fns_cor(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fes(off), vaddr(vidx), ierr) ! 42
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fes(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fes_cor(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fhs(off), vaddr(vidx), ierr) ! 43
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhs(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fhs_cor(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%fwet(off), vaddr(vidx), ierr) ! 44
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fwet(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%epot(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fnpp(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fevw_pot(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%gswx_T(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%cdtq(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%wetfac_cs(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%ga(off), vaddr(vidx), ierr) ! 45
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%ga(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%ga_cor(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%ghflux(off), vaddr(vidx), ierr) ! 46
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%ghflux(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%precis(off), vaddr(vidx), ierr) ! 47
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%precis(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%qscrn(off), vaddr(vidx), ierr) ! 48
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%qscrn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%rnet(off), vaddr(vidx), ierr) ! 49
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%rnet(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%segg(off), vaddr(vidx), ierr) ! 50
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%segg(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%sghflux(off), vaddr(vidx), ierr) ! 51
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%sghflux(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%spill(off), vaddr(vidx), ierr) ! 52
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%spill(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%through(off), vaddr(vidx), ierr) ! 53
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%through(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%tscrn(off), vaddr(vidx), ierr) ! 54
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%tscrn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%tv(off), vaddr(vidx), ierr) ! 55
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%tv(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%us(off), vaddr(vidx), ierr) ! 56
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%us(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%uscrn(off), vaddr(vidx), ierr) ! 57
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%uscrn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%vlaiw(off), vaddr(vidx), ierr) ! 58
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%vlaiw(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%rghlai(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (canopy%wcint(off), vaddr(vidx), ierr) ! 59
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (canopy%wcint(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%fwsoil(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    ! MPI: 2D vars moved above
    ! rwater
    ! evapfbl

    ! ssnow
    ! MPI: 2D vars moved above
    ! albsoilsn
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%pudsto(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%pudsmx(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%cls(off), vaddr(vidx), ierr) ! 60
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%cls(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%dfn_dtg(off), vaddr(vidx), ierr) ! 61
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfn_dtg(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%dfh_dtg(off), vaddr(vidx), ierr) ! 62
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfh_dtg(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfe_ddq(off), displs(bidx), ierr) ! +1
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%dfe_dtg(off), displs(bidx), ierr) ! +1
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%ddq_dtg(off), vaddr(vidx), ierr) ! 63
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ddq_dtg(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%evapsn(off), vaddr(vidx), ierr) ! 64
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%evapsn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%fwtop(off), vaddr(vidx), ierr) ! 65
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop1(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop2(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%fwtop3(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! gammzz
    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (ssnow%isflag(off), vaddr(vidx), ierr) ! 66
    !blen(vidx) = cnt * extid
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%isflag(off), displs(bidx), ierr)
    blocks(bidx) = I1LEN

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%osnowd(off), vaddr(vidx), ierr) ! 67
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%osnowd(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%potev(off), vaddr(vidx), ierr) ! 68
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%potev(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_2)
    !CALL MPI_Get_address (soil%pwb_min(off), vaddr(vidx), ierr) ! 69
    !blen(vidx) = cnt * extr2
    bidx = bidx + 1
    CALL MPI_Get_address (soil%pwb_min(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%runoff(off), vaddr(vidx), ierr) ! 70
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%runoff(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%rnof1(off), vaddr(vidx), ierr) ! 71
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%rnof1(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%rnof2(off), vaddr(vidx), ierr) ! 72
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%rnof2(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%rtsoil(off), vaddr(vidx), ierr) ! 73
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%rtsoil(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! sconds
    ! sdepth
    ! smass
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%snage(off), vaddr(vidx), ierr) ! 74
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%snage(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%snowd(off), vaddr(vidx), ierr) ! 75
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%snowd(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%smelt(off), vaddr(vidx), ierr) ! 76
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%smelt(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! dtmlt
    ! ssdn
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%ssdnn(off), vaddr(vidx), ierr) ! 77
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%ssdnn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! tgg
    ! tggsn
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%tss(off), vaddr(vidx), ierr) ! 78
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tss(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%otss(off), vaddr(vidx), ierr) ! 79
    !blen(vidx) = cnt * extr1
    ! MPI: r1134 does not know about this field, comment out
    !bidx = bidx + 1
    !CALL MPI_Get_address (ssnow%otss(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! wb
    ! wbfice
    ! wbice
    ! wblf
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%wbtot(off), vaddr(vidx), ierr) ! 90
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wbtot(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wb_lake(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%sinfil(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%qstss(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (ssnow%wetfac(off), vaddr(vidx), ierr) ! 91
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wetfac(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    ! MPI: TODO: maybe not needed for transfer to master?
    !CALL MPI_Get_address (ssnow%owetfac(off), vaddr(vidx), ierr) ! 92
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%owetfac(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%t_snwlr(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%tggav(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%otss(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%otss_0(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    ! rad
    ! MPI: 2D vars moved above
    ! albedo
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%extkb(off), vaddr(vidx), ierr) ! 93
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkb(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%extkd2(off), vaddr(vidx), ierr) ! 94
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkd2(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%extkd(off), vaddr(vidx), ierr) ! 95
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%extkd(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%flws(off), vaddr(vidx), ierr) ! 96
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%flws(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! fvlai
    ! gradis
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%latitude(off), vaddr(vidx), ierr) ! 97
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%latitude(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%lwabv(off), vaddr(vidx), ierr) !98
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%lwabv(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 3D vars moved above
    ! qcan
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%qssabs(off), vaddr(vidx), ierr) !99
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%qssabs(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! rhocdf
    ! rniso
    ! scalex
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%transd(off), vaddr(vidx), ierr) ! 100
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%transd(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%trad(off), vaddr(vidx), ierr) ! 101
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (rad%trad(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    bidx = bidx + 1
    CALL MPI_Get_address (rad%transb(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! reffdf
    ! reffbm
    ! extkbm
    ! extkdm
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (rad%fbeam(off), vaddr(vidx), ierr) ! 102
    !blen(vidx) = cnt * extr1
    ! MPI: gol124: changed to 2D and moved up when Bernard
    ! ported to CABLE_r491
    !bidx = bidx + 1
    !CALL MPI_Get_address (rad%fbeam(off,1), displs(bidx), ierr)
    !blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! cexpkbm
    ! cexpkdm

    ! bal
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%drybal(off), vaddr(vidx), ierr) ! 103
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (bal%drybal(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%ebal(off), vaddr(vidx), ierr) ! 104
    !blen(vidx) = cnt * extr1
    ! MPI: remove ebal from exchanged data, calculate temp val on the master
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%ebal(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%ebal_tot(off), vaddr(vidx), ierr) ! 105
    !blen(vidx) = cnt * extr1
    ! MPI: remove ebal_tot from exchanged data, calculate val on the master
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%ebal_tot(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%seb(off), vaddr(vidx), ierr) ! 106
    !blen(vidx) = cnt * extr1
    ! MPI: remove seb from exchanged data, calculate temp val on the master
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%seb(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%seb_tot(off), vaddr(vidx), ierr) ! 107
    !blen(vidx) = cnt * extr1
    ! MPI: remove seb_tot from exchanged data, calculate val on the master
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%seb_tot(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%evap_tot(off), vaddr(vidx), ierr) ! 108
    !blen(vidx) = cnt * extr1
    ! MPI: remove evap_tot from exchanged data
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%evap_tot(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%osnowd0(off), vaddr(vidx), ierr) ! 109
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (bal%osnowd0(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%precip_tot(off), vaddr(vidx), ierr) ! 110
    !blen(vidx) = cnt * extr1
    ! MPI: remove wbal from exchanged data
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%precip_tot(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%rnoff_tot(off), vaddr(vidx), ierr) ! 111
    !blen(vidx) = cnt * extr1
    ! MPI: remove wbal from exchanged data
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%rnoff_tot(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%wbal(off), vaddr(vidx), ierr) ! 112
    !blen(vidx) = cnt * extr1
    ! MPI: remove wbal from exchanged data
    ! bidx = bidx + 1
    ! CALL MPI_Get_address (bal%wbal(off), displs(bidx), ierr)
    ! blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%wbal_tot(off), vaddr(vidx), ierr) ! 113
    !blen(vidx) = cnt * extr1
    ! MPI: remove wbal_tot from exchanged data
    !bidx = bidx + 1
    !CALL MPI_Get_address (bal%wbal_tot(off), displs(bidx), ierr)
    !blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%wbtot0(off), vaddr(vidx), ierr) ! 114
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (bal%wbtot0(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (bal%wetbal(off), vaddr(vidx), ierr) ! 115
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (bal%wetbal(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    ! air
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%rho(off), vaddr(vidx), ierr) ! 116
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%rho(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%volm(off), vaddr(vidx), ierr) ! 117
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%volm(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%rlam(off), vaddr(vidx), ierr) ! 118
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (air%rlam(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%qsat(off), vaddr(vidx), ierr) ! 119
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%qsat(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%epsi(off), vaddr(vidx), ierr) ! 120
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%epsi(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%visc(off), vaddr(vidx), ierr) ! 121
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%visc(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%psyc(off), vaddr(vidx), ierr) ! 122
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%psyc(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%dsatdk(off), vaddr(vidx), ierr) ! 123
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%dsatdk(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (air%cmolar(off), vaddr(vidx), ierr) ! 124
    !blen(vidx) = cnt * extr1

    ! TODO: skip, used for restart but not output
    bidx = bidx + 1
    CALL MPI_Get_address (air%cmolar(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    ! soil
    ! MPI: 2D vars moved above
    ! albsoil
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%bch(off), vaddr(vidx), ierr) ! 125
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%bch(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%c3(off), vaddr(vidx), ierr) ! 126
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%c3(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%clay(off), vaddr(vidx), ierr) ! 127
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%clay(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%cnsd(off), vaddr(vidx), ierr) ! 128
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%cnsd(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%css(off), vaddr(vidx), ierr) ! 129
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%css(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%hsbh(off), vaddr(vidx), ierr) ! 130
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%hsbh(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%hyds(off), vaddr(vidx), ierr) ! 131
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%hyds(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (soil%i2bp3(off), vaddr(vidx), ierr) ! 132
    !blen(vidx) = cnt * extid
    bidx = bidx + 1
    CALL MPI_Get_address (soil%i2bp3(off), displs(bidx), ierr)
    ! Maciej: i2bp3 is REAL
    !    blocks(bidx) = I1LEN
    blocks(bidx) = R1LEN

    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (soil%ibp2(off), vaddr(vidx), ierr) ! 133
    !blen(vidx) = cnt * extid
    bidx = bidx + 1
    CALL MPI_Get_address (soil%ibp2(off), displs(bidx), ierr)
    ! Maciej: ibp2 is REAL
    !    blocks(bidx) = I1LEN
    blocks(bidx) = R1LEN

    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (soil%isoilm(off), vaddr(vidx), ierr) ! 134
    !blen(vidx) = cnt * extid
    bidx = bidx + 1
    CALL MPI_Get_address (soil%isoilm(off), displs(bidx), ierr)
    blocks(bidx) = I1LEN

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%rhosoil(off), vaddr(vidx), ierr) ! 135
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%rhosoil(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%rs20(off), vaddr(vidx), ierr) ! 136
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%rs20(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%sand(off), vaddr(vidx), ierr) ! 137
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%sand(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%sfc(off), vaddr(vidx), ierr) ! 138
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%sfc(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%silt(off), vaddr(vidx), ierr) ! 139
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%silt(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%ssat(off), vaddr(vidx), ierr) ! 140
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%ssat(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%sucs(off), vaddr(vidx), ierr) ! 141
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%sucs(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (soil%swilt(off), vaddr(vidx), ierr) ! 142
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (soil%swilt(off), displs(bidx), ierr)
    blocks(bidx) = r1len


    ! veg
    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (veg%iveg(off), vaddr(vidx), ierr) ! 143
    !blen(vidx) = cnt * extid
    bidx = bidx + 1
    CALL MPI_Get_address (veg%iveg(off), displs(bidx), ierr)
    blocks(bidx) = I1LEN

    !vidx = vidx + 1
    ! INTEGER(i_d)
    !CALL MPI_Get_address (veg%meth(off), vaddr(vidx), ierr) ! 144
    !blen(vidx) = cnt * extid
    bidx = bidx + 1
    CALL MPI_Get_address (veg%meth(off), displs(bidx), ierr)
    ! Maciej: veg%meth is REAL
    !    blocks(bidx) = I1LEN
    blocks(bidx) = R1LEN

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%vlai(off), vaddr(vidx), ierr) ! 145
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%vlai(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    ! MPI: 2D vars moved above
    ! froot
    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%canst1(off), vaddr(vidx), ierr) ! 146
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%canst1(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%ejmax(off), vaddr(vidx), ierr) ! 147
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%ejmax(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%frac4(off), vaddr(vidx), ierr) ! 148
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%frac4(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%wai(off), vaddr(vidx), ierr) ! 149
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%wai(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%vegcf(off), vaddr(vidx), ierr) ! 150
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%vegcf(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%tminvj(off), vaddr(vidx), ierr) ! 151
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%tminvj(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%tmaxvj(off), vaddr(vidx), ierr) ! 152
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%tmaxvj(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%vbeta(off), vaddr(vidx), ierr) ! 153
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%vbeta(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%xalbnir(off), vaddr(vidx), ierr) ! 154
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%xalbnir(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%hc(off), vaddr(vidx), ierr) ! 155
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%hc(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%shelrb(off), vaddr(vidx), ierr) ! 156
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%shelrb(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%vcmax(off), vaddr(vidx), ierr) ! 157
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%vcmax(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%xfang(off), vaddr(vidx), ierr) ! 158
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%xfang(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%dleaf(off), vaddr(vidx), ierr) ! 159
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%dleaf(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%rp20(off), vaddr(vidx), ierr) ! 160
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%rp20(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%rpcoef(off), vaddr(vidx), ierr) ! 161
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%rpcoef(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! REAL(r_1)
    !CALL MPI_Get_address (veg%extkn(off), vaddr(vidx), ierr) ! 162
    !blen(vidx) = cnt * extr1
    bidx = bidx + 1
    CALL MPI_Get_address (veg%extkn(off), displs(bidx), ierr)
    blocks(bidx) = r1len

    !vidx = vidx + 1
    ! LOGICAL
    !CALL MPI_Get_address (veg%deciduous(off), vaddr(vidx), ierr) ! 163
    !blen(vidx) = cnt * extl
    bidx = bidx + 1
    CALL MPI_Get_address (veg%deciduous(off), displs(bidx), ierr)
    blocks(bidx) = llen

    !mrd 1D GW
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%GWwb(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%wtd(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%satfrac(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%Qrecharge(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    ! additional for SLI
    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%Tsurface(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%h0(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%delwcol(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%evap(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%nsnow(off), displs(bidx), ierr)
    blocks(bidx) = i1len

    bidx = bidx + 1
    CALL MPI_Get_address (ssnow%nsteps(off), displs(bidx), ierr)
    blocks(bidx) = r2len
    ! end additional for SLI


    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,': invalid outtype nmat, nvec or n3d constant, fix it!'
       WRITE(*,*) 'bidx: ', bidx
       WRITE(*,*) 'nvec: ', nvec
       WRITE(*,*) 'n3d:', n3d
       WRITE(*,*) 'ntyp', ntyp
       CALL MPI_Abort (comm, 1, ierr)
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

  SUBROUTINE worker_time_update (met, kend, dels)

    USE cable_common_module, ONLY: ktau_gl
    USE cable_def_types_mod
    USE cable_IO_vars_module

    IMPLICIT NONE

    TYPE(met_type), INTENT(INOUT) :: met
    INTEGER, INTENT(IN) :: kend ! number of time steps in simulation
    REAL, INTENT(IN) :: dels ! time step size

    INTEGER :: i

    DO i=1,mland ! over all land points/grid cells
       ! First set timing variables:
       ! All timing details below are initially written to the first patch
       ! of each gridcell, then dumped to all patches for the gridcell.
       IF(ktau_gl==1) THEN ! initialise...
          SELECT CASE(time_coord)
          CASE('LOC')! i.e. use local time by default
             ! hour-of-day = starting hod
             met%hod(landpt(i)%cstart) = shod
             met%doy(landpt(i)%cstart) = sdoy
             met%moy(landpt(i)%cstart) = smoy
             met%year(landpt(i)%cstart) = syear
          CASE('GMT')! use GMT
             ! hour-of-day = starting hod + offset from GMT time:
             met%hod(landpt(i)%cstart) = shod + (longitude(i)/180.0)*12.0
             ! Note above that all met%* vars have dim mp,
             ! while longitude and latitude have dimension mland.
             met%doy(landpt(i)%cstart) = sdoy
             met%moy(landpt(i)%cstart) = smoy
             met%year(landpt(i)%cstart) = syear
          CASE DEFAULT
             CALL abort('Unknown time coordinate! ' &
                  //' (SUBROUTINE get_met_data)')
          END SELECT
       ELSE
          ! increment hour-of-day by time step size:
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) + dels/3600.0
       END IF
       !
       IF(met%hod(landpt(i)%cstart)<0.0) THEN ! may be -ve since longitude
          ! has range [-180,180]
          ! Reduce day-of-year by one and ammend hour-of-day:
          met%doy(landpt(i)%cstart) = met%doy(landpt(i)%cstart) - 1
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) + 24.0
          ! If a leap year AND we're using leap year timing:
          IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. &
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)).AND.leaps) THEN
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(0) ! ie Dec previous year
                met%moy(landpt(i)%cstart) = 12
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) - 1
                met%doy(landpt(i)%cstart) = 365 ! prev year not leap year as this is
             CASE(31) ! Jan
                met%moy(landpt(i)%cstart) = 1
             CASE(60) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(91) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(121)
                met%moy(landpt(i)%cstart) = 4
             CASE(152)
                met%moy(landpt(i)%cstart) = 5
             CASE(182)
                met%moy(landpt(i)%cstart) = 6
             CASE(213)
                met%moy(landpt(i)%cstart) = 7
             CASE(244)
                met%moy(landpt(i)%cstart) = 8
             CASE(274)
                met%moy(landpt(i)%cstart) = 9
             CASE(305)
                met%moy(landpt(i)%cstart) = 10
             CASE(335)
                met%moy(landpt(i)%cstart) = 11
             END SELECT
          ELSE ! not a leap year or not using leap year timing
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(0) ! ie Dec previous year
                met%moy(landpt(i)%cstart) = 12
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) - 1
                ! If previous year is a leap year
                IF((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. &
                     (MOD(syear,4)==0.AND.MOD(syear,400)==0)) THEN
                   met%doy(landpt(i)%cstart) = 366
                ELSE
                   met%doy(landpt(i)%cstart) = 365
                END IF
             CASE(31) ! Jan
                met%moy(landpt(i)%cstart) = 1
             CASE(59) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(90)
                met%moy(landpt(i)%cstart) = 3
             CASE(120)
                met%moy(landpt(i)%cstart) = 4
             CASE(151)
                met%moy(landpt(i)%cstart) = 5
             CASE(181)
                met%moy(landpt(i)%cstart) = 6
             CASE(212)
                met%moy(landpt(i)%cstart) = 7
             CASE(243)
                met%moy(landpt(i)%cstart) = 8
             CASE(273)
                met%moy(landpt(i)%cstart) = 9
             CASE(304)
                met%moy(landpt(i)%cstart) = 10
             CASE(334)
                met%moy(landpt(i)%cstart) = 11
             END SELECT
          END IF ! if leap year or not
       ELSE IF(met%hod(landpt(i)%cstart)>=24.0) THEN
          ! increment or GMT adj has shifted day
          ! Adjust day-of-year and hour-of-day:
          met%doy(landpt(i)%cstart) = met%doy(landpt(i)%cstart) + 1
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) - 24.0
          ! If a leap year AND we're using leap year timing:
          IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. &
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)).AND.leaps) THEN
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(32) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(61) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(92)
                met%moy(landpt(i)%cstart) = 4
             CASE(122)
                met%moy(landpt(i)%cstart) = 5
             CASE(153)
                met%moy(landpt(i)%cstart) = 6
             CASE(183)
                met%moy(landpt(i)%cstart) = 7
             CASE(214)
                met%moy(landpt(i)%cstart) = 8
             CASE(245)
                met%moy(landpt(i)%cstart) = 9
             CASE(275)
                met%moy(landpt(i)%cstart) = 10
             CASE(306)
                met%moy(landpt(i)%cstart) = 11
             CASE(336)
                met%moy(landpt(i)%cstart) = 12
             CASE(367)! end of year; increment
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) + 1
                met%moy(landpt(i)%cstart) = 1
                met%doy(landpt(i)%cstart) = 1
             END SELECT
             ! ELSE IF not leap year and Dec 31st, increment year
          ELSE
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(32) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(60) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(91)
                met%moy(landpt(i)%cstart) = 4
             CASE(121)
                met%moy(landpt(i)%cstart) = 5
             CASE(152)
                met%moy(landpt(i)%cstart) = 6
             CASE(182)
                met%moy(landpt(i)%cstart) = 7
             CASE(213)
                met%moy(landpt(i)%cstart) = 8
             CASE(244)
                met%moy(landpt(i)%cstart) = 9
             CASE(274)
                met%moy(landpt(i)%cstart) = 10
             CASE(305)
                met%moy(landpt(i)%cstart) = 11
             CASE(335)
                met%moy(landpt(i)%cstart) = 12
             CASE(366)! end of year; increment
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) + 1
                met%moy(landpt(i)%cstart) = 1
                met%doy(landpt(i)%cstart) = 1
             END SELECT
          END IF ! if leap year or not
       END IF ! if increment has pushed hod to a different day
       ! Now copy these values to all veg/soil patches in the current grid cell:
       met%hod(landpt(i)%cstart:landpt(i)%cend) = met%hod(landpt(i)%cstart)
       met%doy(landpt(i)%cstart:landpt(i)%cend) = met%doy(landpt(i)%cstart)
       met%moy(landpt(i)%cstart:landpt(i)%cend) = met%moy(landpt(i)%cstart)
       met%year(landpt(i)%cstart:landpt(i)%cend) = met%year(landpt(i)%cstart)
    ENDDO

    RETURN

  END SUBROUTINE worker_time_update

  ! creates MPI types for sending casa results back to the master at
  ! the end of the simulation
  !
  ! changes from old mpi version:
  !  phen: removed out because casa_poolout in this version
  !  is no longer writing phen%phase
  !

  !
  SUBROUTINE worker_casa_type (comm, casapool,casaflux, &
       casamet,casabal, phen)

    USE mpi
    ! USE cable_vars

    USE cable_def_types_mod
    USE casadimension
    USE casavariable
    !  gol124: commented out because casa_poolout in this version
    !  is no longer writing phen%phase
    USE phenvariable

    IMPLICIT NONE

    ! subroutine arguments
    INTEGER :: comm ! MPI communicator to talk to the workers
    TYPE (casa_pool),           INTENT(INOUT) :: casapool
    TYPE (casa_flux),           INTENT(INOUT) :: casaflux
    TYPE (casa_met),            INTENT(INOUT) :: casamet
    TYPE (casa_balance),        INTENT(INOUT) :: casabal
    TYPE (phen_variable),       INTENT(INOUT) :: phen

    ! local variables

    ! MPI: temp arrays for marshalling all types into a struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types
    INTEGER :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    INTEGER :: r1len, r2len, I1LEN, llen

    INTEGER :: rank, off, cnt
    INTEGER :: bidx, midx, vidx, ierr

    INTEGER :: tsize
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

    ! MPI: allocate temp vectors used for marshalling
    ntyp = ncasa_mat + ncasa_vec + (nphen -1)
    ALLOCATE (blocks(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    !off = wpatch%patch0
    !cnt = wpatch%npatch
    off = 1
    cnt = mp

    r1len = cnt * extr1
    r2len = cnt * extr2
    I1LEN = cnt * extid
    llen = cnt * extl

    bidx = 0

    ! ------------- 2D arrays -------------

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%cplant(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%clitter(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%csoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nplant(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nlitter(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nsoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%pplant(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%plitter(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%psoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCplant(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioPCplant(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNClitter(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioPClitter(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioNCsoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%ratioPCsoil(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%cwoodprod(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mwood

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nwoodprod(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mwood

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%pwoodprod(off,1), displs(bidx), ierr)
    blocks(bidx) = r2len * mwood

   ! 

    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphase(off,1), displs(bidx), ierr)
    blocks(bidx) = mphase*I1LEN



    bidx = bidx + 1
    CALL MPI_Get_address (phen%phasespin(off,1), displs(bidx), ierr)
    blocks(bidx) = mdyear*I1LEN


    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_1(off,1), displs(bidx), ierr)
    blocks(bidx) = mdyear*I1LEN


    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_2(off,1), displs(bidx), ierr)
    blocks(bidx) = mdyear*I1LEN


    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_3(off,1), displs(bidx), ierr)
    blocks(bidx) = mdyear*I1LEN


    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphasespin_4(off,1), displs(bidx), ierr)
    blocks(bidx) = mdyear*I1LEN

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Cplant_turnover(off,1), displs(bidx), ierr)
    blocks(bidx) = mplant*r2LEN


    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracCalloc, displs(bidx), ierr)
    blocks(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracNalloc, displs(bidx), ierr)
    blocks(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fracPalloc, displs(bidx), ierr)
    blocks(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Crmplant, displs(bidx), ierr)
    blocks(bidx) = mplant * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%kplant, displs(bidx), ierr)
    blocks(bidx) = mplant * r2len

    ! 3D
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromPtoL, displs(bidx), ierr)
    blocks(bidx) = mplant * mlitter * r2len


    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%klitter, displs(bidx), ierr)
    blocks(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%ksoil, displs(bidx), ierr)
    blocks(bidx) = msoil * r2len

    ! 3D
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromLtoS, displs(bidx), ierr)
    blocks(bidx) = msoil * mlitter * r2len

    ! 3D
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromStoS, displs(bidx), ierr)
    blocks(bidx) = msoil * msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromLtoCO2, displs(bidx), ierr)
    blocks(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fromStoCO2, displs(bidx), ierr)
    blocks(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxCtolitter, displs(bidx), ierr)
    blocks(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxNtolitter, displs(bidx), ierr)
    blocks(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxPtolitter, displs(bidx), ierr)
    blocks(bidx) = mlitter * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxCtosoil, displs(bidx), ierr)
    blocks(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxNtosoil, displs(bidx), ierr)
    blocks(bidx) = msoil * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%fluxPtosoil, displs(bidx), ierr)
    blocks(bidx) = msoil * r2len

















    ! ------------- 1D vectors -------------

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%glai(off), displs(bidx), ierr)
    blocks(bidx) = r2len

    !  gol124: commented out because casa_poolout in this version
    !  is no longer writing phen%phase
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
    !*****************************
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

    !*****************************
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

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker: invalid number of casa fields, fix it!'
       CALL MPI_Abort (comm, 1, ierr)
    END IF

    types = MPI_BYTE

    CALL MPI_Type_create_struct (bidx, blocks, displs, types, casa_t, ierr)
    CALL MPI_Type_commit (casa_t, ierr)

    CALL MPI_Type_size (casa_t, tsize, ierr)
    CALL MPI_Type_get_extent (casa_t, tmplb, text, ierr)

    WRITE (*,*) 'casa type struct blocks, size, extent and lb: ',bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blocks)

    RETURN

  END SUBROUTINE worker_casa_type



  SUBROUTINE worker_climate_types (comm, climate, ktauday )

    USE mpi

    USE cable_def_types_mod, ONLY: climate_type, alloc_cbm_var, mp
    USE cable_climate_mod, ONLY: climate_init

    IMPLICIT NONE

    TYPE(climate_type), INTENT(INOUT):: climate
    INTEGER, INTENT(IN) :: comm, ktauday

    ! MPI: temp arrays for marshalling all types into a struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types
    INTEGER :: ntyp ! number of worker's types

    INTEGER :: last2d, i

    ! MPI: block lenghts for hindexed representing all vectors
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

    ! MPI: block lengths and strides for hvector representing matrices
    INTEGER :: r1len, r2len, I1LEN
    INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

    INTEGER :: tsize, totalrecv, totalsend
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

    INTEGER :: rank, off, cnt
    INTEGER :: bidx, midx, vidx, ierr, ny, nd, ndq
    INTEGER :: stat(MPI_STATUS_SIZE), ierr2, rcount, pos


    CHARACTER, DIMENSION(:), ALLOCATABLE :: rbuf

    ! CALL alloc_cbm_var(climate,mp)

    IF (cable_user%call_climate) CALL climate_init ( climate, mp, ktauday )
    ! MPI: allocate temp vectors used for marshalling
    ntyp = nclimate

    ALLOCATE (blocks(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))
    off = 1

    ! counter to sum total number of bytes receives from all workers
    totalrecv = 0

    r1len = mp * extr1
    r2len = mp * extr2
    I1LEN = mp * extid

    bidx = 0

    ! ------------- 2D arrays -------------
    ny = climate%nyear_average
    nd = climate%nday_average
    ndq = 91

    bidx = bidx + 1
    CALL MPI_Get_address (climate%mtemp_min_20(off,1), displs(bidx), ierr)
    blocks(bidx) = ny*r1len
    types(bidx)  = MPI_BYTE


    bidx = bidx + 1
    CALL MPI_Get_address (climate%mtemp_max_20(off,1), displs(bidx), ierr)
    blocks(bidx) = ny*r1len
    types(bidx)  = MPI_BYTE

    bidx = bidx + 1
    CALL MPI_Get_address (climate%alpha_PT_20(off,1), displs(bidx), ierr)
    blocks(bidx) = ny*r1len
    types(bidx)  = MPI_BYTE


    bidx = bidx + 1
    CALL MPI_Get_address (climate%dtemp_31(off,1), displs(bidx), ierr)
    blocks(bidx) = nd*r1len
    types(bidx)  = MPI_BYTE

    bidx = bidx + 1
    CALL MPI_Get_address (climate%dtemp_91(off,1), displs(bidx), ierr)
    blocks(bidx) = ndq*r1len
    types(bidx)  = MPI_BYTE


    bidx = bidx + 1
    CALL MPI_Get_address (climate%dmoist_31(off,1), displs(bidx), ierr)
    blocks(bidx) = nd*r1len
    types(bidx)  = MPI_BYTE


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

    ! ------------- scalars  -------------

    bidx = bidx + 1
    CALL MPI_Get_address (climate%nyears, displs(bidx), ierr)
    blocks(bidx) = extid
    types(bidx)  = MPI_BYTE

    bidx = bidx + 1
    CALL MPI_Get_address (climate%doy, displs(bidx), ierr)
    blocks(bidx) = extid
    types(bidx)  = MPI_BYTE


!!$types = MPI_BYTE
    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker: invalid number of climate fields, fix it!'
       CALL MPI_Abort (comm, 1, ierr)
    END IF


    CALL MPI_Type_create_struct (bidx, blocks, displs, types, climate_t, ierr)
    CALL MPI_Type_commit (climate_t, ierr)

    CALL MPI_Type_size (climate_t, tsize, ierr)
    CALL MPI_Type_get_extent (climate_t, tmplb, text, ierr)

    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    ! so, now receive all the parameters
    !    CALL MPI_Recv (MPI_BOTTOM, 1, climate_t, 0, 0, comm, stat, ierr)

    !   Maciej: buffered recv + unpac version
    ALLOCATE (rbuf(tsize))
    CALL MPI_Recv (rbuf, tsize, MPI_BYTE, 0, 0, comm, stat, ierr)
    CALL MPI_Get_count (stat, climate_t, rcount, ierr2)

    IF (ierr == MPI_SUCCESS .AND. ierr2 == MPI_SUCCESS .AND. rcount == 1) THEN
       pos = 0
       CALL MPI_Unpack (rbuf, tsize, pos, MPI_BOTTOM, rcount, climate_t, &
            comm, ierr)
       IF (ierr /= MPI_SUCCESS) WRITE(*,*),'climate unpack error, rank: ',rank,ierr
    ELSE
       WRITE(*,*),'climate recv rank err err2 rcount: ',rank, ierr, ierr2, rcount
    END IF

    DEALLOCATE(rbuf)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blocks)

    RETURN

  END SUBROUTINE worker_climate_types
!!$
  ! MPI: creates restart_t type to send to the master the fields
  ! that are only required for the restart file but not included in the
  ! results sent at the end of each time step
  SUBROUTINE worker_restart_type (comm, canopy, air)

    USE mpi

    USE cable_def_types_mod

    IMPLICIT NONE

    INTEGER :: comm

    TYPE(canopy_type), INTENT(IN) :: canopy
    TYPE (air_type),INTENT(IN)     :: air

    ! MPI: temp arrays for marshalling all types into a struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types
    INTEGER :: ntyp ! number of worker's types

    ! MPI: block lengths and strides for hvector representing matrices
    INTEGER :: r1len, r2len, I1LEN, llen

    INTEGER :: rank, off, cnt
    INTEGER :: bidx, midx, vidx, ierr, nd, ny

    INTEGER :: tsize
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

    CALL MPI_Comm_rank (comm, rank, ierr)

    ! MPI: allocate temp vectors used for marshalling
    ntyp = nrestart
    ALLOCATE (blocks(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    off = 1
    cnt = mp
    bidx = 0

    r1len = cnt * extr1
    r2len = cnt * extr2
    I1LEN = cnt * extid
    llen = cnt * extl

    !  bidx = bidx + 1
    !  CALL MPI_Get_address (canopy%rwater(off,1), displs(bidx), ierr)
    !  blocks(bidx) = r1len * ms

    bidx = bidx + 1
    CALL MPI_Get_address (canopy%evapfbl(off,1), displs(bidx), ierr)
    ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
    blocks(bidx) = r1len * ms

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

    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'invalid nrestart constant, fix it!'
       CALL MPI_Abort (comm, 1, ierr)
    END IF

    types = MPI_BYTE

    CALL MPI_Type_create_struct (bidx, blocks, displs, types, restart_t, ierr)
    CALL MPI_Type_commit (restart_t, ierr)

    CALL MPI_Type_size (restart_t, tsize, ierr)
    CALL MPI_Type_get_extent (restart_t, tmplb, text, ierr)

    WRITE (*,*) 'restart struct blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    !mcd287  CALL MPI_Reduce (tsize, tsize, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
    WRITE(*,*) 'b4 reduce wk', tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr
    CALL flush(6)
    !call flush(wlogn)
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blocks)

    RETURN

  END SUBROUTINE worker_restart_type

  SUBROUTINE worker_casa_dump_types(comm, casamet, casaflux, phen, climate)

    USE mpi

    USE casavariable, ONLY: casa_met, casa_flux, mplant
    USE cable_def_types_mod, ONLY: climate_type
    USE phenvariable

    IMPLICIT NONE

    ! sub arguments
    INTEGER, INTENT(IN) :: comm  ! MPI communicator
    TYPE (casa_flux)   , INTENT(INOUT) :: casaflux
    TYPE (casa_met)    , INTENT(INOUT) :: casamet
    TYPE (phen_variable), INTENT(INOUT)  :: phen
    TYPE (climate_type):: climate

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: stat(MPI_STATUS_SIZE), ierr
    INTEGER :: landp_t, patch_t, param_t

    INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank

    CALL MPI_Comm_rank (comm, rank, ierr)

    ntyp = ncdumprw

    ALLOCATE (blen(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid

    bidx = 0

    ! ------- casamet ----

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%Tairk, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%Tsoil, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casamet%moist, displs(bidx), ierr)
    blen(bidx) = ms * r2len

    ! ------- casaflux ----

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Cgpp, displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Crmplant, displs(bidx), ierr)
    blen(bidx) = mplant * r2len



    !****************************************************************
    ! phen fields


    bidx = bidx + 1
    CALL MPI_Get_address (phen%phase, displs(bidx), ierr)
    blen(bidx) = I1LEN


    bidx = bidx + 1
    CALL MPI_Get_address (phen%doyphase, displs(bidx), ierr)
    blen(bidx) = mphase * i1len

    ! #294 - Avoid malformed var write for now 
    ! bidx = bidx + 1
    ! CALL MPI_Get_address (climate%mtemp_max, displs(bidx), ierr)
    ! blen(bidx) = r1len

    !****************************************************************
    ! Ndep
    bidx = bidx + 1
    CALL MPI_Get_address (casaflux%Nmindep, displs(bidx), ierr)
    blen(bidx) = r2len


    !******************************************************************


    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,' invalid number of casa_dump_t param fields ',bidx,', fix it!'
       CALL MPI_Abort (comm, 1, ierr)
    END IF

    CALL MPI_Type_create_struct (bidx, blen, displs, types, casa_dump_t, ierr)
    CALL MPI_Type_commit (casa_dump_t, ierr)

    CALL MPI_Type_size (casa_dump_t, tsize, ierr)
    CALL MPI_Type_get_extent (casa_dump_t, tmplb, text, ierr)

    WRITE (*,*) 'worker casa_dump_t param blocks, size, extent and lb: ',rank, &
         bidx,tsize,text,tmplb

    ! MPI: check whether total size of received data equals total
    ! data sent by all the workers
    CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    DEALLOCATE(types)
    DEALLOCATE(displs)
    DEALLOCATE(blen)

!!$ ! if anything went wrong the master will mpi_abort
!!$ ! which mpi_recv below is going to catch...
!!$ ! so, now receive all the parameters
!!$ CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, 0, comm, stat, ierr)
!!$
!!$ ! finally free the MPI type
!!$ CALL MPI_Type_Free (casa_dump_t, ierr)

    ! all casa parameters have been received from the master by now

  END SUBROUTINE worker_casa_dump_types

  SUBROUTINE worker_casa_LUC_types(comm, casapool, casabal)

    USE mpi

    USE casavariable, ONLY: casa_pool, mplant, mlitter, msoil, casa_balance

    IMPLICIT NONE

    ! sub arguments
    INTEGER, INTENT(IN) :: comm  ! MPI communicator
    TYPE (casa_pool)   , INTENT(IN) :: casapool
    TYPE (casa_balance),        INTENT(IN) :: casabal

    ! local vars

    ! temp arrays for marshalling all fields into a single struct
    INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
    INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: types

    ! temp vars for verifying block number and total length of inp_t
    INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
    INTEGER :: tsize

    INTEGER :: stat(MPI_STATUS_SIZE), ierr
    INTEGER :: landp_t, patch_t, param_t

    INTEGER :: r1len, r2len, I1LEN, llen ! block lengths
    INTEGER :: bidx ! block index
    INTEGER :: ntyp ! total number of blocks

    INTEGER :: rank, off

    CALL MPI_Comm_rank (comm, rank, ierr)

    ntyp = nLUCrw

    ALLOCATE (blen(ntyp))
    ALLOCATE (displs(ntyp))
    ALLOCATE (types(ntyp))

    ! default type is byte, to be overriden for multi-D types
    types = MPI_BYTE

    r1len = mp * extr1
    r2len = mp * extr2
    i1len = mp * extid
    off = 1
    bidx = 0

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%cplant(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%clitter(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%csoil(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nplant(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nlitter(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%nsoil(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%pplant(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * mplant

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%plitter(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * mlitter

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%psoil(off,1), displs(bidx), ierr)
    blen(bidx) = r2len * msoil

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%Nsoilmin(off), displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casapool%clabile(off), displs(bidx), ierr)
    blen(bidx) = r2len

    bidx = bidx + 1
    CALL MPI_Get_address (casabal%FCneeyear(off), displs(bidx), ierr)
    blen(bidx) = r2len



    ! MPI: sanity check
    IF (bidx /= ntyp) THEN
       WRITE (*,*) 'worker ',rank,' invalid number of casa_LUC_t param fields ',bidx,', fix it!'
       CALL MPI_Abort (comm, 1, ierr)
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

!!$ ! if anything went wrong the master will mpi_abort
!!$ ! which mpi_recv below is going to catch...
!!$ ! so, now receive all the parameters
!!$ CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, 0, comm, stat, ierr)
!!$
!!$ ! finally free the MPI type
!!$ CALL MPI_Type_Free (casa_dump_t, ierr)

    ! all casa parameters have been received from the master by now

  END SUBROUTINE worker_casa_LUC_types


  SUBROUTINE worker_pop_types(comm, veg, casamet, pop)

    USE mpi
    USE POP_mpi
    USE POPmodule,          ONLY: ALLOC_POP,POP_INIT
    USE POP_types,          ONLY: pop_type
    USE casavariable,       ONLY: casa_met
    USE cable_def_types_mod,ONLY: veg_parameter_type
    USE cable_common_module,ONLY: cable_user

    IMPLICIT NONE

    INTEGER,         INTENT(IN)    :: comm
    TYPE (casa_met), INTENT(IN)    :: casamet
    TYPE (pop_type), INTENT(INOUT) :: pop
    TYPE (veg_parameter_type),INTENT(IN) :: veg

    INTEGER, DIMENSION(:), ALLOCATABLE :: Iwood, ainv
    INTEGER :: mp_pop, stat(MPI_STATUS_SIZE), ierr

    INTEGER :: rank, inv

    CALL MPI_Comm_rank (comm, rank, ierr)

    ! Get POP relevant info from Master
    CALL MPI_Recv ( mp_pop, 1, MPI_INTEGER, 0, 0, comm, stat, ierr )
    WRITE(*,*),'worker iwood to allocate', rank, mp_pop, mp
    !write(*,*),'worker mppop', rank, mp_pop
    !ALLOCATE( POP%Iwood( mp_pop ) )
    ALLOCATE( Iwood( mp_pop ) )
    WRITE(*,*),'worker iwood allocated', rank, mp_pop
    !CALL MPI_Recv ( POP%Iwood, mp_pop, MPI_INTEGER, 0, 0, comm, stat, ierr )
    CALL MPI_Recv ( Iwood, mp_pop, MPI_INTEGER, 0, 0, comm, stat, ierr )
    !write(*,*),'worker Iwood', rank, POP%Iwood
    ! Maciej
    IF (ANY (Iwood < 1) .OR. ANY (Iwood > mp)) THEN
       WRITE(*,*),'worker iwood values outside valid ranges', rank
       inv = COUNT(Iwood < 1)
       IF (inv .GT. 0) THEN
          WRITE(*,*),'no of values below 1: ', inv
          ALLOCATE (ainv(inv))
          ainv = PACK (Iwood, Iwood .LT. 1)
          WRITE (*,*),'values below 1: ', ainv
          DEALLOCATE (ainv)
       END IF
       inv = COUNT(Iwood > mp)
       IF (inv .GT. 0) THEN
          WRITE(*,*),'no of values above mp ', mp, inv
          ALLOCATE (ainv(inv))
          ainv = PACK (Iwood, Iwood .GT. mp)
          WRITE (*,*),'values above mp: ', ainv
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
       WRITE(*,*),'rank receiving pop_grid from master', rank
       CALL MPI_Recv( POP%pop_grid(1), mp_pop, pop_t, 0, 0, comm, stat, ierr )
    END IF


  END SUBROUTINE worker_pop_types

  SUBROUTINE worker_send_pop (POP, comm)

    USE mpi
    USE POP_mpi
    USE POP_Types,           ONLY: pop_type

    IMPLICIT NONE

    INTEGER         ,INTENT(IN) :: comm
    TYPE (pop_type) ,INTENT(IN) :: pop

    INTEGER :: ierr

    IF ( POP%np .EQ. 0 ) RETURN

    CALL MPI_Send(POP%pop_grid(1), POP%np, pop_t, 0, 0, comm, ierr )


  END SUBROUTINE worker_send_pop


  ! frees memory used for worker's data structures
  SUBROUTINE worker_end(icycle, restart)

    USE mpi

    IMPLICIT NONE

    INTEGER :: icycle ! casa flag
    LOGICAL :: restart

    INTEGER :: ierr

    CALL MPI_Type_free (inp_t, ierr)

    CALL MPI_Type_free (send_t, ierr)

    IF (icycle>0) THEN
       CALL MPI_Type_free (casa_t, ierr)
    END IF

    IF (restart) THEN
       CALL MPI_Type_free (restart_t, ierr)
    END IF

    RETURN

  END SUBROUTINE worker_end


  !*********************************************************************************************

  SUBROUTINE worker_spincasacnp( dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
       casaflux,casamet,casabal,phen,POP,climate,LALLOC,  icomm, ocomm )

    ! USE cable_mpiworker
    USE cable_def_types_mod
    USE cable_carbon_module
    USE cable_common_module, ONLY: CABLE_USER
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    USE POP_Types,  ONLY: POP_TYPE
    USE POPMODULE,            ONLY: POPStep
    USE TypeDef,              ONLY: i4b, dp
    USE mpi

    !mrd561 debug
    USE cable_IO_vars_module, ONLY: wlogn

    IMPLICIT NONE
    !!CLN  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
    REAL,    INTENT(IN)    :: dels
    INTEGER, INTENT(IN)    :: kstart
    INTEGER, INTENT(IN)    :: kend
    INTEGER, INTENT(IN)    :: mloop
    INTEGER, INTENT(IN)    :: LALLOC
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal
    TYPE (phen_variable),         INTENT(INOUT) :: phen
    TYPE (POP_TYPE), INTENT(INOUT)     :: POP
    TYPE (climate_TYPE), INTENT(INOUT)     :: climate



    ! communicator for error-messages
    INTEGER, INTENT(IN)  :: icomm, ocomm
    TYPE (casa_met)  :: casaspin

    ! local variables
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
    !chris 12/oct/2012 for spin up casa
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass
    REAL(r_2), DIMENSION(:), ALLOCATABLE, SAVE  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_af
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_aw
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_ar
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_lf
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_lw
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_lr
  REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_annual_cnpp

    ! local variables
    INTEGER                  :: myearspin,nyear, nloop1, LOY
    CHARACTER(LEN=99)        :: ncfile
    CHARACTER(LEN=4)         :: cyear
    INTEGER                  :: ktau,ktauday,nday,idoy,ktaux,ktauy,nloop
    INTEGER, SAVE            :: ndays
    REAL,      DIMENSION(mp)      :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
    REAL,      DIMENSION(mp)      :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
    REAL,      DIMENSION(mp)      :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
    REAL,      DIMENSION(mp)      :: xcgpp,     xcnpp,     xnuptake,  xpuptake
    REAL,      DIMENSION(mp)      :: xnsoilmin, xpsoillab, xpsoilsorb,xpsoilocc
    REAL(r_2), DIMENSION(mp)      :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

    ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
    REAL,      DIMENSION(5,mvtype,mplant)  :: bmcplant,  bmnplant,  bmpplant
    REAL,      DIMENSION(5,mvtype,mlitter) :: bmclitter, bmnlitter, bmplitter
    REAL,      DIMENSION(5,mvtype,msoil)   :: bmcsoil,   bmnsoil,   bmpsoil
    REAL,      DIMENSION(5,mvtype)         :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
    REAL,      DIMENSION(mvtype)           :: bmarea
    INTEGER nptx,nvt,kloop

    REAL(dp)                               :: StemNPP(mp,2)
    INTEGER, ALLOCATABLE :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles


    INTEGER :: stat(MPI_STATUS_SIZE)
    INTEGER :: ierr



    IF (.NOT.ALLOCATED(Iw)) ALLOCATE(Iw(POP%np))


    !! vh_js !!
    IF (cable_user%CALL_POP) THEN

       Iw = POP%Iwood

    ENDIF

    ktauday=INT(24.0*3600.0/dels)
    nday=(kend-kstart+1)/ktauday
    LOY = 365
    !chris 12/oct/2012 for spin up casa
    IF (.NOT.(ALLOCATED(avg_cleaf2met)))  ALLOCATE(avg_cleaf2met(mp), avg_cleaf2str(mp), avg_croot2met(mp), avg_croot2str(mp), &
         avg_cwood2cwd(mp), &
         avg_nleaf2met(mp), avg_nleaf2str(mp), avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp), &
         avg_pleaf2met(mp), avg_pleaf2str(mp), avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp), &
         avg_cgpp(mp),      avg_cnpp(mp),      avg_nuptake(mp),   avg_puptake(mp),                     &
         avg_xnplimit(mp),  avg_xkNlimiting(mp), avg_xklitter(mp), avg_xksoil(mp),                      &
         avg_rationcsoilmic(mp),avg_rationcsoilslow(mp),avg_rationcsoilpass(mp),                        &
         avg_nsoilmin(mp),  avg_psoillab(mp),    avg_psoilsorb(mp), avg_psoilocc(mp))

  ALLOCATE(avg_af(mp))
  ALLOCATE(avg_aw(mp))
  ALLOCATE(avg_ar(mp))
  ALLOCATE(avg_lf(mp))
  ALLOCATE(avg_lw(mp))
  ALLOCATE(avg_lr(mp))
  ALLOCATE(avg_annual_cnpp(mp))
  avg_af = 0.0
  avg_aw = 0.0
  avg_ar = 0.0
  avg_lf = 0.0
  avg_lw = 0.0
  avg_lr = 0.0
  avg_annual_cnpp = 0.0

    myearspin = CABLE_USER%CASA_SPIN_ENDYEAR - CABLE_USER%CASA_SPIN_STARTYEAR + 1
    ! compute the mean fluxes and residence time of each carbon pool
    avg_cleaf2met=0.0; avg_cleaf2str=0.0; avg_croot2met=0.0; avg_croot2str=0.0; avg_cwood2cwd=0.0
    avg_nleaf2met=0.0; avg_nleaf2str=0.0; avg_nroot2met=0.0; avg_nroot2str=0.0; avg_nwood2cwd=0.0
    avg_pleaf2met=0.0; avg_pleaf2str=0.0; avg_proot2met=0.0; avg_proot2str=0.0; avg_pwood2cwd=0.0
    avg_cgpp=0.0;      avg_cnpp=0.0;      avg_nuptake=0.0;   avg_puptake=0.0
    avg_xnplimit=0.0;  avg_xkNlimiting=0.0; avg_xklitter=0.0; avg_xksoil=0.0
    avg_nsoilmin=0.0;  avg_psoillab=0.0;    avg_psoilsorb=0.0; avg_psoilocc=0.0
    avg_rationcsoilmic=0.0;avg_rationcsoilslow=0.0;avg_rationcsoilpass=0.0

    DO nyear=1,myearspin

       ! WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
       ! ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'

       !call read_casa_dump( ncfile,casamet, casaflux, phen,climate, ktau ,kend,.TRUE. )


       DO idoy=1,mdyear
          ktau=(idoy-1)*ktauday +1
          CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, idoy, icomm, stat, ierr)
          CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
               casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
               xksoil,xkleaf,xkleafcold,xkleafdry,&
               cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
               nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
               pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

          IF (cable_user%CALL_POP .AND. POP%np.GT.0) THEN ! CALL_POP
!!$           ! accumulate annual variables for use in POP
             IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
                casabal%LAImax = casamet%glai
                casabal%Cleafmean = casapool%cplant(:,1)/REAL(LOY)/1000.
                casabal%Crootmean = casapool%cplant(:,3)/REAL(LOY)/1000.
             ELSE
                casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                casabal%LAImax = MAX(casamet%glai, casabal%LAImax)
                casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/REAL(LOY)/1000.
                casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/REAL(LOY)/1000.
             ENDIF

             IF(idoy==mdyear) THEN ! end of year


                CALL POPdriver(casaflux,casabal,veg, POP)


             ENDIF  ! end of year
          ELSE
             casaflux%stemnpp = 0.
          ENDIF ! CALL_POP



          ! WHERE(xkNlimiting .eq. 0)  !Chris Lu 4/June/2012
          !    xkNlimiting = 0.001
          ! END WHERE

        ! Calculate average allocation fractions  (-) for the plant pools
        avg_af = avg_af + casaflux%fracCalloc(:,leaf)
        avg_aw = avg_aw + casaflux%fracCalloc(:,wood)
        avg_ar = avg_ar + casaflux%fracCalloc(:,froot)

        ! Calculate average turnover rates for the plant pools (yr-1)
        avg_lf = avg_lf + (casaflux%kplant(:,leaf) * REAL(LOY))
        avg_lw = avg_lw + (casaflux%kplant(:,wood) * REAL(LOY))
        avg_lr = avg_lr + (casaflux%kplant(:,froot) * REAL(LOY))

          avg_cleaf2met = avg_cleaf2met + cleaf2met
          avg_cleaf2str = avg_cleaf2str + cleaf2str
          avg_croot2met = avg_croot2met + croot2met
          avg_croot2str = avg_croot2str + croot2str
          avg_cwood2cwd = avg_cwood2cwd + cwood2cwd

          avg_nleaf2met = avg_nleaf2met + nleaf2met
          avg_nleaf2str = avg_nleaf2str + nleaf2str
          avg_nroot2met = avg_nroot2met + nroot2met
          avg_nroot2str = avg_nroot2str + nroot2str
          avg_nwood2cwd = avg_nwood2cwd + nwood2cwd

          avg_pleaf2met = avg_pleaf2met + pleaf2met
          avg_pleaf2str = avg_pleaf2str + pleaf2str
          avg_proot2met = avg_proot2met + proot2met
          avg_proot2str = avg_proot2str + proot2str
          avg_pwood2cwd = avg_pwood2cwd + pwood2cwd

          avg_cgpp      = avg_cgpp      + casaflux%cgpp
          avg_cnpp      = avg_cnpp      + casaflux%cnpp
          avg_nuptake   = avg_nuptake   + casaflux%Nminuptake
          avg_puptake   = avg_puptake   + casaflux%Plabuptake

          avg_xnplimit    = avg_xnplimit    + xnplimit
          avg_xkNlimiting = avg_xkNlimiting + xkNlimiting
          avg_xklitter    = avg_xklitter    + xklitter
          avg_xksoil      = avg_xksoil      + xksoil

          avg_nsoilmin    = avg_nsoilmin    + casapool%nsoilmin
          avg_psoillab    = avg_psoillab    + casapool%psoillab
          avg_psoilsorb   = avg_psoilsorb   + casapool%psoilsorb
          avg_psoilocc    = avg_psoilocc    + casapool%psoilocc

          avg_rationcsoilmic  = avg_rationcsoilmic  + casapool%ratioNCsoilnew(:,mic)
          avg_rationcsoilslow = avg_rationcsoilslow + casapool%ratioNCsoilnew(:,slow)
          avg_rationcsoilpass = avg_rationcsoilpass + casapool%ratioNCsoilnew(:,pass)
       ENDDO
    ENDDO

 ! Average the plant allocation fraction
  avg_af = avg_af / REAL(nday)
  avg_aw = avg_aw / REAL(nday)
  avg_ar = avg_ar / REAL(nday)

  ! Average the plant turnover fraction
  avg_lf = avg_lf / REAL(nday)
  avg_lw = avg_lw / REAL(nday)
  avg_lr = avg_lr / REAL(nday)

  ! Need the annual NPP to solve plant pools g C m-2 y-1
  avg_annual_cnpp = avg_cnpp / REAL(myearspin)

  avg_cleaf2met = avg_cleaf2met/REAL(nday)
  avg_cleaf2str = avg_cleaf2str/REAL(nday)
  avg_croot2met = avg_croot2met/REAL(nday)
  avg_croot2str = avg_croot2str/REAL(nday)
  avg_cwood2cwd = avg_cwood2cwd/REAL(nday)

  avg_nleaf2met = avg_nleaf2met/REAL(nday)
  avg_nleaf2str = avg_nleaf2str/REAL(nday)
  avg_nroot2met = avg_nroot2met/REAL(nday)
  avg_nroot2str = avg_nroot2str/REAL(nday)
  avg_nwood2cwd = avg_nwood2cwd/REAL(nday)

  avg_pleaf2met = avg_pleaf2met/REAL(nday)
  avg_pleaf2str = avg_pleaf2str/REAL(nday)
  avg_proot2met = avg_proot2met/REAL(nday)
  avg_proot2str = avg_proot2str/REAL(nday)
  avg_pwood2cwd = avg_pwood2cwd/REAL(nday)

  avg_cgpp      = avg_cgpp/REAL(nday)
  avg_cnpp      = avg_cnpp/REAL(nday)

  avg_nuptake   = avg_nuptake/REAL(nday)
  avg_puptake   = avg_puptake/REAL(nday)

  avg_xnplimit    = avg_xnplimit/REAL(nday)
  avg_xkNlimiting = avg_xkNlimiting/REAL(nday)
  avg_xklitter    = avg_xklitter/REAL(nday)

  avg_xksoil      = avg_xksoil/REAL(nday)

  avg_nsoilmin    = avg_nsoilmin/REAL(nday)
  avg_psoillab    = avg_psoillab/REAL(nday)
  avg_psoilsorb   = avg_psoilsorb/REAL(nday)
  avg_psoilocc    = avg_psoilocc/REAL(nday)

  avg_rationcsoilmic  = avg_rationcsoilmic  /REAL(nday)
  avg_rationcsoilslow = avg_rationcsoilslow /REAL(nday)
  avg_rationcsoilpass = avg_rationcsoilpass /REAL(nday)

    CALL analyticpool(kend,veg,soil,casabiome,casapool,                                          &
         casaflux,casamet,casabal,phen,                                         &
         avg_cleaf2met,avg_cleaf2str,avg_croot2met,avg_croot2str,avg_cwood2cwd, &
         avg_nleaf2met,avg_nleaf2str,avg_nroot2met,avg_nroot2str,avg_nwood2cwd, &
         avg_pleaf2met,avg_pleaf2str,avg_proot2met,avg_proot2str,avg_pwood2cwd, &
         avg_cgpp, avg_cnpp, avg_nuptake, avg_puptake,                          &
         avg_xnplimit,avg_xkNlimiting,avg_xklitter,avg_xksoil,                  &
         avg_ratioNCsoilmic,avg_ratioNCsoilslow,avg_ratioNCsoilpass,            &
       avg_nsoilmin,avg_psoillab,avg_psoilsorb,avg_psoilocc,                  &
       avg_af, avg_aw, avg_ar, avg_lf, avg_lw, avg_lr, avg_annual_cnpp)


    nloop1= MAX(1,mloop-3)

    DO nloop=1,mloop

       !!CLN  OPEN(91,file=fcnpspin)
       !!CLN  read(91,*)
       DO nyear=1,myearspin

          ! WRITE(CYEAR,FMT="(I4)") CABLE_USER%CASA_SPIN_STARTYEAR + nyear - 1
          ! ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
          ! call read_casa_dump( ncfile, casamet, casaflux, phen,climate, ktau, kend, .TRUE. )

          DO idoy=1,mdyear
             ktauy=idoy*ktauday
             ktau=(idoy-1)*ktauday +1
             CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, idoy, icomm, stat, ierr)

             CALL biogeochem(ktauy,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                  casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,&
                  xkleafcold,xkleafdry,&
                  cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                  nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                  pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)



             IF (cable_user%CALL_POP .AND. POP%np.GT.0) THEN ! CALL_POP

                ! accumulate annual variables for use in POP
                IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                   casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                   ! (assumes 70% of wood NPP is allocated above ground)
                   casabal%LAImax = casamet%glai
                   casabal%Cleafmean = casapool%cplant(:,1)/REAL(LOY)/1000.
                   casabal%Crootmean = casapool%cplant(:,3)/REAL(LOY)/1000.
                ELSE
                   casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * &
                        casaflux%fracCalloc(:,2) * 0.7
                   casabal%LAImax = MAX(casamet%glai, casabal%LAImax)
                   casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/REAL(LOY)/1000.
                   casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/REAL(LOY)/1000.
                ENDIF


                IF(idoy==mdyear) THEN ! end of year

                   CALL POPdriver(casaflux,casabal,veg, POP)

                ENDIF  ! end of year
             ELSE
                casaflux%stemnpp = 0.
             ENDIF ! CALL_POP


          ENDDO   ! end of idoy
       ENDDO   ! end of nyear

    ENDDO     ! end of nloop
    WRITE(wlogn,*) 'b4 MPI_SEND'
    CALL MPI_Send (MPI_BOTTOM, 1, casa_t, 0, 0, ocomm, ierr)
    WRITE(wlogn,*) 'after MPI_SEND'
    IF(CABLE_USER%CALL_POP) CALL worker_send_pop (POP, ocomm)
    WRITE(wlogn,*) 'cplant', casapool%cplant

  END SUBROUTINE worker_spincasacnp

  !*********************************************************************************************

  SUBROUTINE worker_CASAONLY_LUC( dels,kstart,kend,veg,soil,casabiome,casapool, &
       casaflux,casamet,casabal,phen,POP,climate,LALLOC,  icomm, ocomm )

    ! USE cable_mpiworker
    USE cable_def_types_mod
    USE cable_carbon_module
    USE cable_common_module, ONLY: CABLE_USER
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    USE POP_Types,  ONLY: POP_TYPE
    USE POPMODULE,            ONLY: POPStep
    USE TypeDef,              ONLY: i4b, dp
    USE mpi

    !mrd561 debug
    USE cable_IO_vars_module, ONLY: wlogn

    IMPLICIT NONE
    !!CLN  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
    REAL,    INTENT(IN)    :: dels
    INTEGER, INTENT(IN)    :: kstart
    INTEGER, INTENT(IN)    :: kend
    INTEGER, INTENT(IN)    :: LALLOC
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal
    TYPE (phen_variable),         INTENT(INOUT) :: phen
    TYPE (POP_TYPE), INTENT(INOUT)     :: POP
    TYPE (climate_TYPE), INTENT(INOUT)     :: climate



    ! communicator for error-messages
    INTEGER, INTENT(IN)  :: icomm, ocomm
    TYPE (casa_met)  :: casaspin

    ! local variables
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
    !chris 12/oct/2012 for spin up casa
    REAL,      DIMENSION(:), ALLOCATABLE, SAVE  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass
    REAL(r_2), DIMENSION(:), ALLOCATABLE, SAVE  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

    ! local variables
    INTEGER                  :: myearspin,nyear, nloop1
    CHARACTER(LEN=99)        :: ncfile
    CHARACTER(LEN=4)         :: cyear
    INTEGER                  :: ktau,ktauday,nday,idoy,ktaux,ktauy,nloop
    INTEGER, SAVE            :: ndays
    REAL,      DIMENSION(mp)      :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
    REAL,      DIMENSION(mp)      :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
    REAL,      DIMENSION(mp)      :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
    REAL,      DIMENSION(mp)      :: xcgpp,     xcnpp,     xnuptake,  xpuptake
    REAL,      DIMENSION(mp)      :: xnsoilmin, xpsoillab, xpsoilsorb,xpsoilocc
    REAL(r_2), DIMENSION(mp)      :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

    ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
    REAL,      DIMENSION(5,mvtype,mplant)  :: bmcplant,  bmnplant,  bmpplant
    REAL,      DIMENSION(5,mvtype,mlitter) :: bmclitter, bmnlitter, bmplitter
    REAL,      DIMENSION(5,mvtype,msoil)   :: bmcsoil,   bmnsoil,   bmpsoil
    REAL,      DIMENSION(5,mvtype)         :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
    REAL,      DIMENSION(mvtype)           :: bmarea
    INTEGER nptx,nvt,kloop

    REAL(dp)                               :: StemNPP(mp,2)


    INTEGER :: stat(MPI_STATUS_SIZE)
    INTEGER :: ierr, rank
    INTEGER :: yyyy



    ktauday=INT(24.0*3600.0/dels)
    nday=(kend-kstart+1)/ktauday

    myearspin = CABLE_USER%YEAREND - CABLE_USER%YEARSTART + 1
    yyyy = CABLE_USER%YEARSTART - 1

    DO nyear=1,myearspin
       DO idoy=1,mdyear
          ktau=(idoy-1)*ktauday +1
          CALL MPI_Recv (MPI_BOTTOM, 1, casa_dump_t, 0, idoy, icomm, stat, ierr)


          CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
               casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter, &
               xksoil,xkleaf,xkleafcold,xkleafdry,&
               cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
               nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
               pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)


          ! accumulate annual variables for use in POP
          IF(idoy==1 ) THEN
             casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
             casabal%LAImax = casamet%glai
             casabal%Cleafmean = casapool%cplant(:,1)/REAL(mdyear)/1000.
             casabal%Crootmean = casapool%cplant(:,3)/REAL(mdyear)/1000.
          ELSE
             casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
             casabal%LAImax = MAX(casamet%glai, casabal%LAImax)
             casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/REAL(mdyear)/1000.
             casabal%Crootmean = casabal%Crootmean +casapool%cplant(:,3)/REAL(mdyear)/1000.
          ENDIF

          IF(idoy==mdyear) THEN ! end of year
             WRITE(wlogn,*) 'b4 MPI_SEND,casa_LUC_t', casapool%cplant(:,2)
             CALL flush(wlogn)
             CALL MPI_Send (MPI_BOTTOM, 1, casa_LUC_t, 0, 0, ocomm, ierr)
             WRITE(wlogn,*) 'after MPI_SEND,casa_LUC_t', casapool%cplant(:,2)
             CALL flush(wlogn)
             StemNPP(:,1) = casaflux%stemnpp
             StemNPP(:,2) = 0.0

             CALL MPI_Comm_rank (icomm, rank, ierr)
             WRITE(wlogn,*)
             WRITE(wlogn,*),'rank receiving pop_grid from master', rank
!!$           write(wlogn,*) 'b4 MPI_Recv, pop_t cmass: ', POP%pop_grid%cmass_sum
!!$           write(wlogn,*) 'b4 MPI_Recv, pop_t LU: ', POP%pop_grid%LU
             CALL MPI_Recv( POP%pop_grid(1), POP%np, pop_t, 0, 0, icomm, stat, ierr )
!!$           write(wlogn,*)
!!$           write(wlogn,*) 'after MPI_Recv, pop_t cmass: ', POP%pop_grid%cmass_sum
             WRITE(wlogn,*) 'after MPI_Recv, pop_t '
             CALL flush(wlogn)
             IF (cable_user%CALL_POP .AND. POP%np.GT.0) THEN ! CALL_POP
                WRITE(wlogn,*), 'b4  POPdriver', POP%pop_grid%cmass_sum
                CALL POPdriver(casaflux,casabal,veg, POP)

             ENDIF
!!$           write(wlogn,*)
!!$           write(wlogn,*) 'after POPstep cmass: ', POP%pop_grid%cmass_sum
             WRITE(wlogn,*) 'after POPstep ',  POP%pop_grid%cmass_sum
             CALL flush(wlogn)
             CALL worker_send_pop (POP, ocomm)
             WRITE(wlogn,*) 'after worker_send_pop'
             CALL flush(wlogn)
          ENDIF


       ENDDO
       ! receive updates to CASA pools resulting from LUC
       WRITE(wlogn,*)
       WRITE(wlogn,*) 'b4 mpi_recv casa_LUC_t '
       CALL MPI_Recv (MPI_BOTTOM, 1, casa_LUC_t, 0, nyear, icomm, stat, ierr)
       WRITE(wlogn,*) 'after mpi_recv casa_LUC_t: '
    ENDDO

  END SUBROUTINE WORKER_CASAONLY_LUC


END MODULE cable_mpiworker
