!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
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

  ! worker's struct for restart data to the master
  INTEGER :: restart_t

  PUBLIC :: mpidrv_worker

CONTAINS

   SUBROUTINE mpidrv_worker (comm)

   USE mpi

   USE cable_def_types_mod
   USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps,                  &
                                   verbose, fixedCO2,output,check,patchout,    &
                                   patch_type,soilparmnew
   USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user,     &
                                   cable_runtime, filename,                    & 
                                   redistrb, wiltParam, satuParam
   USE cable_data_module,    ONLY: driver_type, point2constants
   USE cable_input_module,   ONLY: open_met_file,load_parameters,              &
                                   get_met_data,close_met_file
   USE cable_output_module,  ONLY: create_restart,open_output_file,            &
                                   write_output,close_output_file
   USE cable_cbm_module
   
   ! modules related to CASA-CNP
   USE casadimension,       ONLY: icycle 
   USE casavariable,        ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
                                  casa_met, casa_balance
   USE phenvariable,        ONLY: phen_variable

   IMPLICIT NONE

   ! MPI:
   INTEGER               :: comm ! MPI communicator for comms with the workers

   ! CABLE namelist: model configuration, runtime/user switches 
   CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml' 
   
   ! timing variables 
   INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
   
   INTEGER        ::                                                           &
      ktau,       &  ! increment equates to timestep, resets if spinning up
      ktau_tot,   &  ! NO reset when spinning up, total timesteps by model
      kend,       &  ! no. of time steps in run
      ktauday,    &  ! day counter for CASA-CNP
      idoy,       &  ! day of year (1:365) counter for CASA-CNP
      nyear          ! year counter for CASA-CNP

   REAL :: dels                        ! time step size in seconds
   
   ! CABLE variables
   TYPE (met_type)       :: met     ! met input variables
   TYPE (air_type)       :: air     ! air property variables
   TYPE (canopy_type)    :: canopy  ! vegetation variables
   TYPE (radiation_type) :: rad     ! radiation variables
   TYPE (roughness_type) :: rough   ! roughness varibles
   TYPE (balances_type)  :: bal     ! energy and water balance variables
   TYPE (soil_snow_type) :: ssnow   ! soil and snow variables
   
   ! CABLE parameters
   TYPE (soil_parameter_type) :: soil ! soil parameters	
   TYPE (veg_parameter_type)  :: veg  ! vegetation parameters	 
   TYPE (driver_type)    :: C         ! constants used locally  
   
   TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
   
   ! CASA-CNP variables 
   TYPE (casa_biome)     :: casabiome
   TYPE (casa_pool)      :: casapool
   TYPE (casa_flux)      :: casaflux
   TYPE (casa_met)       :: casamet
   TYPE (casa_balance)   :: casabal
   TYPE (phen_variable)  :: phen 
   
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
      l_vcmaxFeedbk = .FALSE.       ! using prognostic Vcmax
   
   
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

   ! switches etc defined thru namelist (by default cable.nml)
   NAMELIST/CABLE/                  &
                  filename,         & ! TYPE, containing input filenames 
                  vegparmnew,       & ! jhan: use new soil param. method
                  soilparmnew,      & ! jhan: use new soil param. method
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

   ! END header

   ! Open, read and close the namelist file.
   OPEN( 10, FILE = CABLE_NAMELIST )
      READ( 10, NML=CABLE )   !where NML=CABLE defined above
   CLOSE(10)

   IF( IARGC() > 0 ) THEN
      CALL GETARG(1, filename%met)
      CALL GETARG(2, casafile%cnpipool)
   ENDIF

    
   cable_runtime%offline = .TRUE.
   
   ! associate pointers used locally with global definitions
   CALL point2constants( C )
    
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

   ! MPI: TODO: find a way to preserve workers log messages somewhere
   ! (either separate files or collated by the master to a single file
   ! or perhaps use MPI-IO - but probably not gonna work with random length
   ! text strings)
   OPEN(logn,FILE="/dev/null")
 
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

   ! MPI: bcast to workers so that they don't need to open the met
   ! file themselves
   CALL MPI_Bcast (dels, 1, MPI_REAL, 0, comm, ierr)
   CALL MPI_Bcast (kend, 1, MPI_INTEGER, 0, comm, ierr)

   ! MPI: receive from master starting time fields
   !CALL bcast_start_time (comm)

   ! MPI: need to know extents before creating datatypes
   CALL find_extents

   ! MPI: receive decomposition info from the master
   call worker_decomp(comm)

   ! MPI: in overlap version sends and receives occur on separate comms
   CALL MPI_Comm_dup (comm, icomm, ierr)
   CALL MPI_Comm_dup (comm, ocomm, ierr)

   ! MPI: data set in load_parameter is now received from
   ! the master
   CALL worker_cable_params(comm, met,air,ssnow,veg,bgc,soil,canopy,&
   &                        rough,rad,sum_flux,bal)

   ! MPI: casa parameters received only if cnp module is active
   IF (icycle>0) THEN
     ! MPI:
     CALL worker_casa_params (comm,casabiome,casapool,casaflux,casamet,&
     &                        casabal,phen)
   END IF

   ! MPI: create inp_t type to receive input data from the master
   ! at the start of every timestep
   CALL worker_intype (comm,met,veg)

   ! MPI: create send_t type to send the results to the master
   ! at the end of every timestep
   CALL worker_outtype (comm,met,canopy,ssnow,rad,bal,air,soil,veg)

   ! MPI: create type to send casa results back to the master
   ! only if cnp module is active
   IF (icycle>0) THEN
      CALL worker_casa_type (comm, casapool,casaflux, &
                          casamet,casabal, phen)
   END IF

   ! MPI: create type to send restart data back to the master
   ! only if restart file is to be created
   IF(output%restart) THEN
      CALL worker_restart_type (comm, canopy, air)
   END IF

   ! Open output file:
   ! MPI: only the master writes to the files
   !CALL open_output_file( dels, soil, veg, bgc, rough )
 
   ssnow%otss_0 = ssnow%tgg(:,1)
   ssnow%otss = ssnow%tgg(:,1)
   canopy%fes_cor = 0.
   canopy%fhs_cor = 0.
   met%ofsd = 0.1
   
   ! outer loop - spinup loop no. ktau_tot :
   ktau_tot = 0 
   DO

      ! globally (WRT code) accessible kend through USE cable_common_module
      ktau_gl = 0
      kend_gl = kend
      knode_gl = 0
      
      ! time step loop over ktau
      DO ktau=kstart, kend 
         
         ! increment total timstep counter
         ktau_tot = ktau_tot + 1
         
         ! globally (WRT code) accessible kend through USE cable_common_module
         ktau_gl = ktau_gl + 1
         
         ! somethings (e.g. CASA-CNP) only need to be done once per day  
         ktauday=int(24.0*3600.0/dels)
         idoy = mod(ktau/ktauday,365)
         IF(idoy==0) idoy=365
         
         ! needed for CASA-CNP
         nyear =INT((kend-kstart+1)/(365*ktauday))
   
         canopy%oldcansto=canopy%cansto
   
         ! Get met data and LAI, set time variables.
         ! Rainfall input may be augmented for spinup purposes:
         met%ofsd = met%fsd(:,1) + met%fsd(:,2)
         ! MPI: input file read on the master only
         !CALL get_met_data( spinup, spinConv, met, soil,                    &
         !                   rad, veg, kend, dels, C%TFRZ, ktau ) 

         ! MPI: receive input data for this step from the master
         CALL MPI_Recv (MPI_BOTTOM, 1, inp_t, 0, ktau_gl, icomm, stat, ierr)

         ! MPI: some fields need explicit init, because we don't transfer
         ! them for better performance
         ! in the serial version this is done in get_met_data
         ! after input has been read from the file
         met%tvair = met%tk
         met%tvrad = met%tk

         ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
         IF (l_vcmaxFeedbk) CALL casa_feedback( ktau, veg, casabiome,    &
                                                casapool, casamet )
   
         IF (l_laiFeedbk) veg%vlai(:) = casamet%glai(:)
   
         ! CALL land surface scheme for this timestep, all grid points:
         CALL cbm( dels, air, bgc, canopy, met,                             &
                   bal, rad, rough, soil, ssnow,                            &
                   sum_flux, veg )
   
         ssnow%smelt = ssnow%smelt*dels
         ssnow%rnof1 = ssnow%rnof1*dels
         ssnow%rnof2 = ssnow%rnof2*dels
         ssnow%runoff = ssnow%runoff*dels
   
   
         !jhan this is insufficient testing. condition for 
         !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
         IF(icycle >0) THEN
            call bgcdriver( ktau, kstart, kend, dels, met,                  &
                            ssnow, canopy, veg, soil, casabiome,               &
                            casapool, casaflux, casamet, casabal,              &
                            phen, spinConv, spinup, ktauday, idoy,             &
                            .FALSE., .FALSE. )
         ENDIF 
   
         ! sumcflux is pulled out of subroutine cbm
         ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
         CALL sumcflux( ktau, kstart, kend, dels, bgc,                      &
                        canopy, soil, ssnow, sum_flux, veg,                    &
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
   
       END DO ! END Do loop over timestep ktau



   
      !!jhan this is insufficient testing. condition for 
      !!spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
      !! see if spinup (if conducting one) has converged:
      !IF(spinup.AND..NOT.spinConv) THEN
      !   
      !   ! Write to screen and log file:
      !   WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),      &
      !         ' of data set complete...'
      !   WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),   &
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

   END DO

   IF (icycle > 0) THEN

      ! MPI: send casa results back to the master
      CALL MPI_Send (MPI_BOTTOM, 1, casa_t, 0, ktau_gl, comm, ierr)
      ! MPI: output file written by master only
      !CALL casa_poolout( ktau, veg, soil, casabiome,                           &
      !                   casapool, casaflux, casamet, casabal, phen )
      ! MPI: output file written by master only
      !CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
  
   END IF

   ! Write restart file if requested:
   IF(output%restart) THEN
      ! MPI: send variables that are required by create_restart
      CALL MPI_Send (MPI_BOTTOM, 1, restart_t, 0, ktau_gl, comm, ierr)
      ! MPI: output file written by master only
      !CALL create_restart( logn, dels, ktau, soil, veg, ssnow,                 &
      !                     canopy, rough, rad, bgc, bal )
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
   CLOSE(logn)

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

  INTEGER :: r1len, r2len, ilen, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks

  INTEGER :: rank

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

  ALLOCATE (blen(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  ! default type is byte, to be overriden for multi-D types
  types = MPI_BYTE

  r1len = mp * extr1
  r2len = mp * extr2
  ilen = mp * extid
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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (met%moy, displs(bidx), ierr)
  blen(bidx) = ilen

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

  bidx = bidx + 1
  CALL MPI_Get_address (ssnow%dfe_ddq, displs(bidx), ierr)
  blen(bidx) = r1len

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
  blen(bidx) = ilen

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (veg%meth, displs(bidx), ierr)
  blen(bidx) = ilen

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

  bidx = bidx + 1
  CALL MPI_Get_address (veg%vlai, displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%xfang, displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%extkn, displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%deciduous, displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%refl, displs(bidx), ierr)
  blen(bidx) = 2 * r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%taul, displs(bidx), ierr)
  blen(bidx) = 2 * r1len

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (soil%ibp2, displs(bidx), ierr)
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (soil%isoilm, displs(bidx), ierr)
  blen(bidx) = r1len

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

  bidx = bidx + 1
  CALL MPI_Get_address (soil%zse, displs(bidx), ierr)
  blen(bidx) = ms * extr1

  bidx = bidx + 1
  CALL MPI_Get_address (soil%zshh, displs(bidx), ierr)
  blen(bidx) = (ms + 1) * extr1

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

  bidx = bidx + 1
  CALL MPI_Get_address (canopy%fes, displs(bidx), ierr)
  blen(bidx) = r2len

  bidx = bidx + 1
  CALL MPI_Get_address (canopy%fhs, displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (canopy%fwet, displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (canopy%ga, displs(bidx), ierr)
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
  blen(bidx) = nrb * r1len

  bidx = bidx + 1
  CALL MPI_Get_address (rad%cexpkdm, displs(bidx), ierr)
  blen(bidx) = nrb * r1len

  bidx = bidx + 1
  CALL MPI_Get_address (rad%rhocbm, displs(bidx), ierr)
  blen(bidx) = swb * r1len

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

  bidx = bidx + 1
  CALL MPI_Get_address (ssnow%wbtot1, displs(bidx), ierr)
  blen(bidx) = r1len

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
  CALL MPI_Recv (MPI_BOTTOM, 1, param_t, 0, 0, comm, stat, ierr)

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

  INTEGER :: r1len, r2len, ilen, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks

  INTEGER :: rank

  CALL MPI_Bcast (mvtype, 1, MPI_INTEGER, 0, comm, ierr)
  CALL MPI_Bcast (mstype, 1, MPI_INTEGER, 0, comm, ierr)

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
  ilen = mp * extid
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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%Tsoil, displs(bidx), ierr)
  blen(bidx) = ms * r2len

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%moist, displs(bidx), ierr)
  blen(bidx) = ms * r2len

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%iveg2, displs(bidx), ierr)
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%ijgcm, displs(bidx), ierr)
  blen(bidx) = r2len

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%isorder, displs(bidx), ierr)
  blen(bidx) = r2len

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (phen%TKshed, displs(bidx), ierr)
  blen(bidx) = mvtype * extr2

  bidx = bidx + 1
  CALL MPI_Get_address (phen%doyphase, displs(bidx), ierr)
  blen(bidx) = mphase * ilen

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
  CALL MPI_Recv (MPI_BOTTOM, 1, casa_t, 0, 0, comm, stat, ierr)

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

  INTEGER :: r1len, r2len, ilen, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks
  INTEGER :: ierr

  INTEGER :: rank

  CALL MPI_Comm_rank (comm, rank, ierr)

  r1len = mp * extr1
  r2len = mp * extr2
  ilen = mp * extid
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

  bidx = 1
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
   blocks(bidx) = ilen

   bidx = bidx + 1
   CALL MPI_Get_address (met%moy, displs(bidx), ierr)
   blocks(bidx) = ilen

   bidx = bidx + 1
   CALL MPI_Get_address (met%doy, displs(bidx), ierr)
   blocks(bidx) = r1len

   bidx = bidx + 1
   CALL MPI_Get_address (met%hod, displs(bidx), ierr)
   blocks(bidx) = r1len


  ! MPI: sanity check
  IF (bidx /= ntyp) THEN
     WRITE (*,*) 'worker ',rank,': invalid intype nmat, nvec or n3d constant, fix it!'
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
  INTEGER :: r1len, r2len, ilen, llen

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
  ilen = cnt * extid
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
  blocks(bidx) = r1len * nrb

  !midx = midx + 1
  ! REAL(r_1)
  !CALL MPI_Get_address (rad%cexpkdm(off,1), maddr(midx), ierr) ! 27
  !CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
  !  &            mat_t(midx, rank), ierr)
  bidx = bidx + 1
  CALL MPI_Get_address (rad%cexpkdm(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * nrb

  bidx = bidx + 1
  CALL MPI_Get_address (rad%rhocbm(off,1), displs(bidx), ierr)
  blocks(bidx) = r1len * swb


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
  !blocks(bidx) = ilen

  !vidx = vidx + 1
  ! INTEGER(i_d)
  !CALL MPI_Get_address (met%moy(off), vaddr(vidx), ierr) ! 3
  !blen(vidx) = cnt * extid
  ! gol124: not output, removed
  !bidx = bidx + 1
  !CALL MPI_Get_address (met%moy(off), displs(bidx), ierr)
  !blocks(bidx) = ilen

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
  blocks(bidx) = ilen

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
  blocks(bidx) = ilen

  !vidx = vidx + 1
  ! INTEGER(i_d)
  !CALL MPI_Get_address (soil%ibp2(off), vaddr(vidx), ierr) ! 133
  !blen(vidx) = cnt * extid
  bidx = bidx + 1
  CALL MPI_Get_address (soil%ibp2(off), displs(bidx), ierr)
  blocks(bidx) = ilen

  !vidx = vidx + 1
  ! INTEGER(i_d)
  !CALL MPI_Get_address (soil%isoilm(off), vaddr(vidx), ierr) ! 134
  !blen(vidx) = cnt * extid
  bidx = bidx + 1
  CALL MPI_Get_address (soil%isoilm(off), displs(bidx), ierr)
  blocks(bidx) = ilen

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
  blocks(bidx) = ilen

  !vidx = vidx + 1
  ! INTEGER(i_d)
  !CALL MPI_Get_address (veg%meth(off), vaddr(vidx), ierr) ! 144
  !blen(vidx) = cnt * extid
  bidx = bidx + 1
  CALL MPI_Get_address (veg%meth(off), displs(bidx), ierr)
  blocks(bidx) = ilen

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

  ! MPI: sanity check
  IF (bidx /= ntyp) THEN
     WRITE (*,*) 'worker ',rank,': invalid outtype nmat, nvec or n3d constant, fix it!'
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
  INTEGER :: r1len, r2len, ilen, llen

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr

  INTEGER :: tsize
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  ! MPI: allocate temp vectors used for marshalling
  ntyp = ncasa_mat + ncasa_vec
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  !off = wpatch%patch0
  !cnt = wpatch%npatch
  off = 1
  cnt = mp

  r1len = cnt * extr1
  r2len = cnt * extr2
  ilen = cnt * extid
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

  ! ------------- 1D vectors -------------

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%glai(off), displs(bidx), ierr)
  blocks(bidx) = r2len

!  gol124: commented out because casa_poolout in this version
!  is no longer writing phen%phase
  bidx = bidx + 1
  CALL MPI_Get_address (phen%phase(off), displs(bidx), ierr)
  blocks(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (casapool%clabile(off), displs(bidx), ierr)
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

  WRITE (*,*) 'struct blocks, size, extent and lb: ',bidx,tsize,text,tmplb

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN

END SUBROUTINE worker_casa_type

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
  INTEGER :: r1len, r2len, ilen, llen

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr

  INTEGER :: tsize
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

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
  ilen = cnt * extid
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

  WRITE (*,*) 'restart struct blocks, size, extent and lb: ',bidx,tsize,text,tmplb

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
!mcd287  CALL MPI_Reduce (tsize, tsize, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN

END SUBROUTINE worker_restart_type

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

END MODULE cable_mpiworker

