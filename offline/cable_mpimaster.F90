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

  ! TODO: m3d_t mat_t and vec_t to be factored out from here and from
  ! master_outtypes
  ! MPI: slices of 3D arrays
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: m3d_t
  ! MPI: slices of matrices (2D)
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mat_t

  ! MPI: parts of vectors (1D)
  ! MPI: vec_t dimension is wnp; as each worker gets a single hindexed
  ! with nvec blocks
  INTEGER, ALLOCATABLE, DIMENSION(:) :: vec_t

  ! MPI derived datatype handles for sending input data to the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: inp_ts

  ! MPI derived datatype handles for receiving output from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_ts

  ! MPI derived datatype handles for receiving casa results from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_ts

  ! master's struct for receiving restart data from the workers
  INTEGER, ALLOCATABLE, DIMENSION(:) :: restart_ts

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

SUBROUTINE mpidrv_master (comm)

   USE mpi

   USE cable_def_types_mod
   USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps,                  &
                                   verbose, fixedCO2,output,check,patchout,    &
                                   patch_type,soilparmnew
   USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user,     &
                                   cable_runtime, filename, redistrb,          & 
                                   report_version_no, wiltParam, satuParam
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
      nyear,      &  ! year counter for CASA-CNP
      maxdiff(2)     ! location of maximum in convergence test

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
   TYPE (met_type)       :: imet  ! read ahead met input variables
   TYPE (veg_parameter_type) :: iveg  ! MPI read ahead vegetation parameters
   LOGICAL :: loop_exit     ! MPI: exit flag for bcast to workers
   INTEGER :: iktau    ! read ahead index of time step = 1 ..  kend
   INTEGER :: oktau    ! ktau = 1 ..  kend for output
   INTEGER :: tmp_kgl    ! temp for ktau_gl
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

   ! Open log file:
   OPEN(logn,FILE=filename%log)
 
   CALL report_version_no( logn )
    
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

   ! Check for gswp run
   IF (ncciy /= 0) THEN
      
      PRINT *, 'Looking for global offline run info.'
      
      IF (ncciy < 1986 .OR. ncciy > 1995) THEN
         PRINT *, 'Year ', ncciy, ' outside range of dataset!'
         PRINT *, 'Please check input in namelist file.'

         ! MPI: stop all the workers, too!
         CALL MPI_Abort (comm, 1, ierr)

         STOP
      ELSE
         
         CALL prepareFiles(ncciy)
      
      ENDIF
   
   ENDIF
   

   ! Open met data and get site information from netcdf file.
   ! This retrieves time step size, number of timesteps, starting date,
   ! latitudes, longitudes, number of sites. 
   CALL open_met_file( dels, kend, spinup, C%TFRZ )
 
   ! Checks where parameters and initialisations should be loaded from.
   ! If they can be found in either the met file or restart file, they will 
   ! load from there, with the met file taking precedence. Otherwise, they'll
   ! be chosen from a coarse global grid of veg and soil types, based on 
   ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
   CALL load_parameters( met, air, ssnow, veg, bgc,                            &
                         soil, canopy, rough, rad, sum_flux,                   &
                         bal, logn, vegparmnew, casabiome, casapool,           &
                         casaflux, casamet, casabal, phen, C%EMSOIL,        &
                         C%TFRZ )

   ! MPI: above was standard serial code
   ! now it's time to initialize the workers

   ! MPI: bcast to workers so that they don't need to open the met
   ! file themselves
   CALL MPI_Bcast (dels, 1, MPI_REAL, 0, comm, ierr)
   CALL MPI_Bcast (kend, 1, MPI_INTEGER, 0, comm, ierr)

   ! MPI: need to know extents before creating datatypes
   CALL find_extents

   ! MPI: calculate and broadcast landpoint decomposition to the workers
   CALL master_decomp(comm, mland, mp)

   ! MPI: set up stuff for new irecv isend code that separates completion
   ! from posting of requests
   ALLOCATE (inp_req(wnp))
   ALLOCATE (inp_stats(MPI_STATUS_SIZE, wnp))
   ALLOCATE (recv_req(wnp))
   ALLOCATE (recv_stats(MPI_STATUS_SIZE, wnp))
   CALL MPI_Comm_dup (comm, icomm, ierr)
   CALL MPI_Comm_dup (comm, ocomm, ierr)

   ! MPI: data set in load_parameter is now scattered out to the
   ! workers
   CALL master_cable_params(comm, met,air,ssnow,veg,bgc,soil,canopy,&
   &                         rough,rad,sum_flux,bal)

   ! MPI: casa parameters scattered only if cnp module is active
   IF (icycle>0) THEN
     ! MPI:
     CALL master_casa_params (comm,casabiome,casapool,casaflux,casamet,&
     &                        casabal,phen)
   END IF

   ! MPI: allocate read ahead buffers for input met and veg data
   CALL alloc_cbm_var (imet, mp)
   CALL alloc_cbm_var (iveg, mp)

   ! MPI: create inp_t types to scatter input data to the workers
   ! at the start of every timestep
   !CALL master_intypes (comm,met,veg)
   ! for read ahead use the new variables
   CALL master_intypes (comm,imet,iveg)

   ! MPI: create out_t types to receive results from the workers
   ! at the end of every timestep
   CALL master_outtypes (comm,met,canopy,ssnow,rad,bal,air,soil,veg)

   ! MPI: create type for receiving casa results
   ! only if cnp module is active
   IF (icycle>0) THEN
      CALL master_casa_types (comm, casapool, casaflux, &
                              casamet, casabal, phen)
   END IF

   ! MPI: create type to send restart data back to the master
   ! only if restart file is to be created
   IF(output%restart) THEN
      CALL master_restart_types (comm, canopy, air)
   END IF

   ! MPI: mostly original serial code follows...

   ! Open output file:
   CALL open_output_file( dels, soil, veg, bgc, rough )
 
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

      ! MPI: separate time step counters for reading and writing
      ! (ugly, I know)
      iktau = ktau_gl
      oktau = ktau_gl

      ! MPI: read ahead
      iktau = iktau + 1

      canopy%oldcansto=canopy%cansto

      imet%ofsd = imet%fsd(:,1) + imet%fsd(:,2)

      ! MPI: flip ktau_gl
      !tmp_kgl = ktau_gl
      ktau_gl = iktau

      CALL get_met_data( spinup, spinConv, imet, soil,                    &
                         rad, iveg, kend, dels, C%TFRZ, iktau )

      ! MPI: scatter input data to the workers
      CALL master_send_input (icomm, inp_ts, iktau)
      CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)


      ! time step loop over ktau
      DO ktau=kstart, kend - 1

!         ! increment total timstep counter
!         ktau_tot = ktau_tot + 1
         iktau = iktau + 1
         oktau = oktau + 1

         met%year = imet%year
         met%doy = imet%doy

         ! globally (WRT code) accessible kend through USE cable_common_module
         ktau_tot= ktau_tot + 1
         ktau_gl = iktau
         
         ! somethings (e.g. CASA-CNP) only need to be done once per day  
         ktauday=int(24.0*3600.0/dels)
         idoy = mod(ktau/ktauday,365)
         IF(idoy==0) idoy=365
         
         ! needed for CASA-CNP
         nyear =INT((kend-kstart+1)/(365*ktauday))
   
!         canopy%oldcansto=canopy%cansto
   
         ! Get met data and LAI, set time variables.
         ! Rainfall input may be augmented for spinup purposes:
!          met%ofsd = met%fsd(:,1) + met%fsd(:,2)
         CALL get_met_data( spinup, spinConv, imet, soil,                    &
                            rad, iveg, kend, dels, C%TFRZ, iktau ) 

         ! MPI: receive this time step's results from the workers
         CALL master_receive (ocomm, oktau, recv_ts)

         ! MPI: scatter input data to the workers
         CALL master_send_input (icomm, inp_ts, iktau)

         CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
         CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

         met%ofsd = met%fsd(:,1) + met%fsd(:,2)
         canopy%oldcansto=canopy%cansto

         ! Write time step's output to file if either: we're not spinning up 
         ! or we're spinning up and the spinup has converged:
         ! MPI: TODO: pull mass and energy balance calculation from write_output
         ! and refactor into worker code
         ktau_gl = oktau
         IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
            CALL write_output( dels, ktau, met, canopy, ssnow,              &
                               rad, bal, air, soil, veg, C%SBOLTZ, &
                               C%EMLEAF, C%EMSOIL )
   
       END DO ! END Do loop over timestep ktau


    ! MPI: read ahead tail to receive (last step and write)
    met%year = imet%year
    met%doy = imet%doy
    oktau = oktau + 1
    ktau_tot= ktau_tot + 1
    ktau_gl = oktau
    CALL master_receive (ocomm, oktau, recv_ts)
    CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)
    met%ofsd = met%fsd(:,1) + met%fsd(:,2)
    canopy%oldcansto=canopy%cansto
    IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
            CALL write_output( dels, ktau, met, canopy, ssnow,         &
                               rad, bal, air, soil, veg, C%SBOLTZ,     &
                               C%EMLEAF, C%EMSOIL )

   
      !jhan this is insufficient testing. condition for 
      !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
      ! see if spinup (if conducting one) has converged:
      IF(spinup.AND..NOT.spinConv) THEN
         
         ! Write to screen and log file:
         WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),      &
               ' of data set complete...'
         WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),   &
               ' of data set complete...'
         
         ! IF not 1st run through whole dataset:
         IF( INT( ktau_tot/kend ) > 1 ) THEN 
            
            ! evaluate spinup
            IF( ANY( ABS(ssnow%wb-soilMtemp)>delsoilM).OR.                     &
                ANY(ABS(ssnow%tgg-soilTtemp)>delsoilT) ) THEN
               
               ! No complete convergence yet
!               PRINT *, 'ssnow%wb : ', ssnow%wb
!               PRINT *, 'soilMtemp: ', soilMtemp
!               PRINT *, 'ssnow%tgg: ', ssnow%tgg
!               PRINT *, 'soilTtemp: ', soilTtemp
               maxdiff = MAXLOC(ABS(ssnow%wb-soilMtemp))
               PRINT *, 'Example location of moisture non-convergence: ',maxdiff
               PRINT *, 'ssnow%wb : ', ssnow%wb(maxdiff(1),maxdiff(2))
               PRINT *, 'soilMtemp: ', soilMtemp(maxdiff(1),maxdiff(2))
               maxdiff = MAXLOC(ABS(ssnow%tgg-soilTtemp))
               PRINT *, 'Example location of temperature non-convergence: ',maxdiff
               PRINT *, 'ssnow%tgg: ', ssnow%tgg(maxdiff(1),maxdiff(2))
               PRINT *, 'soilTtemp: ', soilTtemp(maxdiff(1),maxdiff(2))

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
         
           ALLOCATE( soilMtemp(mp,ms), soilTtemp(mp,ms) )
         
         END IF
         
         ! store soil moisture and temperature
         soilTtemp = ssnow%tgg
         soilMtemp = REAL(ssnow%wb)

         ! MPI:
         loop_exit = .FALSE.

      ELSE

         ! if not spinning up, or spin up has converged, exit:
         ! EXIT
         ! MPI:
         loop_exit = .TRUE.

      END IF

      ! MPI: let the workers know whether it's time to quit
      CALL MPI_Bcast (loop_exit, 1, MPI_LOGICAL, 0, comm, ierr)

      IF (loop_exit) THEN
            EXIT
      END IF

   END DO

   IF (icycle > 0) THEN
      ! MPI: gather casa results from all the workers
      CALL master_receive (comm, ktau_gl, casa_ts)
      CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)

      CALL casa_poolout( ktau, veg, soil, casabiome,                           &
                         casapool, casaflux, casamet, casabal, phen )

      CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
  
   END IF

   ! Write restart file if requested:
   IF(output%restart) THEN
      ! MPI: TODO: receive variables that are required by create_restart
      ! but not write_output
      !CALL receive_restart (comm,ktau,dels,soil,veg,ssnow, &
      !       &              canopy,rough,rad,bgc,bal)
      ! gol124: how about call master_receive (comm, ktau, restart_ts)
      ! instead of a separate receive_restart sub?
      CALL master_receive (comm, ktau_gl, restart_ts)
      CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)

      CALL create_restart( logn, dels, ktau, soil, veg, ssnow,                 &
                           canopy, rough, rad, bgc, bal )
   END IF

   ! MPI: cleanup
   call master_end (icycle, output%restart)
      
   ! Close met data input file:
   CALL close_met_file
 
   ! Close output file and deallocate main variables:
   CALL close_output_file( bal, air, bgc, canopy, met,                         &
                           rad, rough, soil, ssnow,                            &
                           sum_flux, veg )

   WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck

   ! Close log file
   CLOSE(logn)

   RETURN

END SUBROUTINE mpidrv_master


SUBROUTINE prepareFiles(ncciy)
  USE cable_IO_vars_module, ONLY: logn,gswpfile
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncciy

  WRITE(logn,*) 'CABLE offline global run using gswp forcing for ', ncciy
  PRINT *,      'CABLE offline global run using gswp forcing for ', ncciy

  CALL renameFiles(logn,gswpfile%rainf,16,ncciy,'rainf')
  CALL renameFiles(logn,gswpfile%snowf,16,ncciy,'snowf')
  CALL renameFiles(logn,gswpfile%LWdown,16,ncciy,'LWdown')
  CALL renameFiles(logn,gswpfile%SWdown,16,ncciy,'SWdown')
  CALL renameFiles(logn,gswpfile%PSurf,16,ncciy,'PSurf')
  CALL renameFiles(logn,gswpfile%Qair,14,ncciy,'Qair')
  CALL renameFiles(logn,gswpfile%Tair,14,ncciy,'Tair')
  CALL renameFiles(logn,gswpfile%wind,15,ncciy,'wind')

END SUBROUTINE prepareFiles


SUBROUTINE renameFiles(logn,inFile,nn,ncciy,inName)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: logn
  INTEGER, INTENT(IN) :: nn
  INTEGER, INTENT(IN) :: ncciy
  CHARACTER(LEN=99), INTENT(INOUT) :: inFile
  CHARACTER(LEN=*),  INTENT(IN)    :: inName
  INTEGER :: idummy

  READ(inFile(nn:nn+3),'(i4)') idummy
  IF (idummy < 1983 .OR. idummy > 1995) THEN
    PRINT *, 'Check position of the year number in input gswp file', inFile
    STOP
  ELSE
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
  ENDIF

END SUBROUTINE renameFiles


! ============== PRIVATE SUBROUTINES USED ONLY BY THE MPI MASTER ===============


! MPI: calculates and sends grid decomposition info to the workers
SUBROUTINE master_decomp (comm, mland, mp)

  USE mpi

  USE cable_IO_vars_module, ONLY : landpt, patch

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: comm ! MPI communicator to talk to the workers
  INTEGER, INTENT(IN)   :: mland ! total number of landpoints in the global grid
  INTEGER, INTENT(IN)   :: mp ! total number of land patches in the global grid

  INTEGER :: lpw  ! average number of landpoints per worker
  INTEGER :: rank, rest, nxt, pcnt, ierr, i, tmp
  INTEGER :: patchcnt  ! sum of patches for a range landpoints

  ! how many workers do we have?
  CALL MPI_Comm_size (comm, wnp, ierr)
  wnp = wnp - 1

  ALLOCATE (wland(wnp), STAT=ierr)
  IF (ierr /= 0) THEN
          ! TODO: print an error message
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
     ! CALL MPI_Send (nxt, 1, MPI_INTEGER, rank, 0, comm, ierr)
     CALL MPI_Send (pcnt, 1, MPI_INTEGER, rank, 0, comm, ierr)

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
        WRITE (*,*) 'invalid patch number for worker ', &
        &           patchcnt,tmp,rank
        CALL MPI_Abort (comm, 0, ierr)
     END IF

     ! MPI: at this point workers can't determine patches on their own
     ! so we need to send it explicitely
     ! in this version of cable workers care only about the number of patches
     ! CALL MPI_Send (wland(rank)%patch0, 1, MPI_INTEGER, rank, 0, comm, ierr)
     CALL MPI_Send (wland(rank)%npatch, 1, MPI_INTEGER, rank, 0, comm, ierr)

     nxt = nxt + pcnt
  END DO

  RETURN

END SUBROUTINE master_decomp


! MPI: creates param_t type for the master to scatter the default parameters
! to the workers
! then sends the parameters
! and finally frees the MPI type
SUBROUTINE master_cable_params (comm,met,air,ssnow,veg,bgc,soil,canopy,&
                                rough,rad,sum_flux,bal)

  USE mpi

  USE cable_def_types_mod
  USE cable_IO_vars_module

  IMPLICIT NONE

  ! subroutine arguments

  INTEGER, INTENT(IN) :: comm ! MPI communicator

  TYPE (met_type), INTENT(INOUT) :: met
  TYPE (air_type), INTENT(INOUT) :: air
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
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: stat(MPI_STATUS_SIZE), ierr
  INTEGER, ALLOCATABLE, DIMENSION(:) :: param_ts

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride
  INTEGER :: r1len, r2len, ilen, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks
  INTEGER :: rank

  INTEGER :: landpt_t, patch_t

  INTEGER :: nxt, pcnt, off, cnt

  ! create MPI types for exchanging slices of landpt and patch arrays
  CALL decomp_types (landpt_t, patch_t)

  ! MPI: TODO: replace sends with isends
  DO rank = 1, wnp

     ! MPI: send worker's landpts
     nxt = wland(rank)%landp0
     pcnt = wland(rank)%nland
     CALL MPI_Send (landpt(nxt), pcnt, landpt_t, rank, 0, comm, ierr)

     ! MPI: send worker's patch
     nxt = wland(rank)%patch0
     pcnt = wland(rank)%npatch
     CALL MPI_Send (patch(nxt), pcnt, patch_t, rank, 0, comm, ierr)

  END DO

  ! MPI: TODO: free landp_t and patch_t types?

  ntyp = nparam

  ALLOCATE (param_ts(wnp))

  ALLOCATE (blen(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  ! MPI: array strides for multi-dimensional types
  r1stride = mp * extr1
  r2stride = mp * extr2

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
  ilen = cnt * extid
  llen = cnt * extl

  bidx = 0

  ! the order of variables follows argument list
  ! the order of fields within follows alloc_*_type subroutines

  ! ----------- met --------------

  bidx = bidx + 1
  CALL MPI_Get_address (met%ca(off), displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (met%year(off), displs(bidx), ierr)
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (met%moy(off), displs(bidx), ierr)
  blen(bidx) = ilen

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
  blen(bidx) = ilen

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (veg%meth(off), displs(bidx), ierr)
  blen(bidx) = ilen

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

  bidx = bidx + 1
  CALL MPI_Get_address (veg%vlai(off), displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%xfang(off), displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%extkn(off), displs(bidx), ierr)
  blen(bidx) = r1len

  bidx = bidx + 1
  CALL MPI_Get_address (veg%deciduous(off), displs(bidx), ierr)
  blen(bidx) = r1len

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (soil%ibp2(off), displs(bidx), ierr)
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (soil%isoilm(off), displs(bidx), ierr)
  blen(bidx) = r1len

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

  ! constant * ms, each worker gets the same copy of whole array
  bidx = bidx + 1
  CALL MPI_Get_address (soil%zse, displs(bidx), ierr)
  blen(bidx) = ms * extr1

  ! constant * (ms+1), each worker gets the same copy of whole array
  bidx = bidx + 1
  CALL MPI_Get_address (soil%zshh, displs(bidx), ierr)
  blen(bidx) = (ms + 1) * extr1

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
  CALL MPI_Get_address (canopy%gswx(off,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (mf, r1len, r1stride, MPI_BYTE, &
  &                             types(bidx), ierr)
  blen(bidx) = 1

  bidx = bidx + 1
  CALL MPI_Get_address (canopy%zetar(off,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (niter, r1len, r1stride, MPI_BYTE, &
  &                             types(bidx), ierr)
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
  CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
  &                             types(bidx), ierr)
  blen(bidx) = 1
  !blen(bidx) = nrb * r1len

  bidx = bidx + 1
  CALL MPI_Get_address (rad%cexpkdm(off,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
  &                             types(bidx), ierr)
  blen(bidx) = 1
  !blen(bidx) = nrb * r1len

  bidx = bidx + 1
  CALL MPI_Get_address (rad%rhocbm(off,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
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

  bidx = bidx + 1
  CALL MPI_Get_address (ssnow%wbtot1(off), displs(bidx), ierr)
  blen(bidx) = r1len

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
  IF (bidx /= ntyp) THEN
     WRITE (*,*) 'master: invalid number of param_t fields ',bidx,', fix it!'
     CALL MPI_Abort (comm, 1, ierr)
  END IF

  CALL MPI_Type_create_struct (bidx, blen, displs, types, param_ts(rank), ierr)
  CALL MPI_Type_commit (param_ts(rank), ierr)

  CALL MPI_Type_size (param_ts(rank), tsize, ierr)
  CALL MPI_Type_get_extent (param_ts(rank), tmplb, text, ierr)

  WRITE (*,*) 'master to rank param_t blocks, size, extent and lb: ',bidx,tsize,text,tmplb

  localtotal = localtotal + tsize

  END DO ! rank

  WRITE (*,*) 'total cable params size sent to all workers: ', localtotal

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blen)

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  CALL MPI_Reduce (MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE (*,*) 'total cable params size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
          WRITE (*,*) 'error: total length of cable params sent and received differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  ! so, now send all the parameters
  CALL master_send_input (comm, param_ts, 0)
  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank = 1, wnp
     CALL MPI_Type_Free (param_ts(rank), ierr)
  END DO

  DEALLOCATE (param_ts)

  ! all CABLE parameters have been transferred to the workers by now
  RETURN

END SUBROUTINE master_cable_params


! MPI: creates casa_ts types to broadcast/scatter the default casa parameters
! to all the workers
! then sends them
! and finally frees the MPI type
SUBROUTINE master_casa_params (comm,casabiome,casapool,casaflux,casamet,&
                               casabal,phen)

  USE mpi

  USE cable_def_types_mod

  USE casavariable
  USE phenvariable

  IMPLICIT NONE

  ! sub arguments
  INTEGER, INTENT(IN) :: comm  ! MPI communicator

  ! TODO: have these variables been already allocated?
  TYPE (casa_biome)   , INTENT(INOUT) :: casabiome
  TYPE (casa_pool)   , INTENT(INOUT) :: casapool
  TYPE (casa_flux)   , INTENT(INOUT) :: casaflux
  TYPE (casa_met)    , INTENT(INOUT) :: casamet
  TYPE (casa_balance), INTENT(INOUT) :: casabal
  TYPE (phen_variable), INTENT(INOUT)  :: phen 

  ! local vars

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: stat(MPI_STATUS_SIZE), ierr
  ! INTEGER :: landp_t, patch_t, param_t
  INTEGER, ALLOCATABLE, DIMENSION(:) :: casa_ts

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride, istride
  INTEGER :: r1len, r2len, ilen, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks

  INTEGER :: rank, off, cnt

  CALL MPI_Bcast (mvtype, 1, MPI_INTEGER, 0, comm, ierr)
  CALL MPI_Bcast (mstype, 1, MPI_INTEGER, 0, comm, ierr)

  ntyp = ncasaparam

  ALLOCATE (casa_ts(wnp))

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
  ilen = cnt * extid
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
  !blen(bidx) = mplant * r2len

  ! gol124: temp
  ! 3D
  bidx = bidx + 1
  CALL MPI_Get_address (casaflux%fromPtoL(off,1,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (mplant * mlitter, r2len, r2stride, MPI_BYTE, &
  &                             types(bidx), ierr)
  blen(bidx) = 1
  !blen(bidx) = mplant * mlitter * r2len

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
  CALL MPI_Get_address (casaflux%klitter(off,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (mlitter, r2len, r2stride, MPI_BYTE, &
  &                             types(bidx), ierr)
  blen(bidx) = 1
  !blen(bidx) = mlitter * r2len

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
  blen(bidx) = ilen

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%ijgcm(off), displs(bidx), ierr)
  blen(bidx) = r2len

  bidx = bidx + 1
  CALL MPI_Get_address (casamet%isorder(off), displs(bidx), ierr)
  blen(bidx) = r2len

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
  blen(bidx) = ilen

  bidx = bidx + 1
  CALL MPI_Get_address (phen%TKshed, displs(bidx), ierr)
  blen(bidx) = mvtype * extr2

  bidx = bidx + 1
  CALL MPI_Get_address (phen%doyphase(off,1), displs(bidx), ierr)
  CALL MPI_Type_create_hvector (mphase, ilen, istride, MPI_BYTE, &
  &                             types(bidx), ierr)
  blen(bidx) = 1
  !blen(bidx) = mphase * ilen

  ! MPI: sanity check
  IF (bidx /= ntyp) THEN
     WRITE (*,*) 'master: invalid number of casa_t param fields ',bidx,', fix it!'
     CALL MPI_Abort (comm, 1, ierr)
  END IF

  CALL MPI_Type_create_struct (bidx, blen, displs, types, casa_ts(rank), ierr)
  CALL MPI_Type_commit (casa_ts(rank), ierr)

  CALL MPI_Type_size (casa_ts(rank), tsize, ierr)
  CALL MPI_Type_get_extent (casa_ts(rank), tmplb, text, ierr)

  WRITE (*,*) 'master to rank casa_t param blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb

  localtotal = localtotal + tsize

  END DO ! rank

  WRITE (*,*) 'total casa params size sent to all workers: ', localtotal
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blen)

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  remotetotal = 0
  CALL MPI_Reduce (MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE (*,*) 'total casa params size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
          WRITE (*,*) 'error: total length of casa params sent and received differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  CALL MPI_Barrier (comm, ierr)

  ! so, now send all the parameters
  CALL master_send_input (comm, casa_ts, 0)
  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank = 1, wnp
     CALL MPI_Type_Free (casa_ts(rank), ierr)
  END DO

  DEALLOCATE (casa_ts)

  ! all casa parameters have been sent to the workers by now

  RETURN

END SUBROUTINE master_casa_params


! MPI: creates inp_t types to send input data to the workers
! input data means arrays read by get_met_data; each worker receives
! only its own slice of the arrays
SUBROUTINE master_intypes (comm,met,veg)

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
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride
  INTEGER :: r1len, r2len, ilen, llen ! block lengths
  INTEGER :: bidx ! block index
  INTEGER :: ntyp ! total number of blocks
  INTEGER :: rank ! worker rank
  INTEGER :: off  ! first patch index for a worker
  INTEGER :: cnt  ! mp for a worker
  INTEGER :: ierr

  ALLOCATE (inp_ts(wnp))

  ! max total number of fields to receive: met + veg fields
  !ntyp = 10 + 1
  ntyp = ninput

  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

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

     r1len = cnt * extr1
     r2len = cnt * extr2
     ilen = cnt * extid
     llen = cnt * extl

     r1stride = mp * extr1

     ! met fields

     bidx = 1
     CALL MPI_Get_address (met%fsd(off,1), displs(bidx), ierr)
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1
     !blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%tk(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%pmb(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%qv(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%ua(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%precip(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%precip_sn(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%fld(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%ca(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%coszen(off), displs(bidx), ierr)
     blocks(bidx) = r1len


     ! veg fields

     bidx = bidx + 1
     CALL MPI_Get_address (veg%vlai(off), displs(bidx), ierr)
     blocks(bidx) = r1len

    ! additional field missing from previous versions;
    ! added when trying to fix a bug in the new mpi code
    ! the order of these new fields follows the order of their
    ! declaration in cable_define_types.F90

     bidx = bidx + 1
     CALL MPI_Get_address (met%year(off), displs(bidx), ierr)
     blocks(bidx) = ilen

     bidx = bidx + 1
     CALL MPI_Get_address (met%moy(off), displs(bidx), ierr)
     blocks(bidx) = ilen

     bidx = bidx + 1
     CALL MPI_Get_address (met%doy(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     bidx = bidx + 1
     CALL MPI_Get_address (met%hod(off), displs(bidx), ierr)
     blocks(bidx) = r1len


     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'master: invalid intype nmat, nvec or n3d constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF


     ! marshall all fields into a single MPI derived datatype for worker rank
  
     CALL MPI_Type_create_struct (bidx, blocks, displs, types, inp_ts(rank), ierr)
     CALL MPI_Type_commit (inp_ts(rank), ierr)

     CALL MPI_Type_size (inp_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (inp_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'master to ',rank,': intype struct blocks, size, extent and lb: ', &
                 bidx,tsize,text,tmplb

     localtotal = localtotal + tsize

  END DO ! rank

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  WRITE (*,*) 'total input data size sent to all workers: ', localtotal

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  CALL MPI_Reduce (MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE (*,*) 'total input data size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
          WRITE (*,*) 'error: total length of input data sent and received differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  RETURN

END SUBROUTINE master_intypes

! MPI: creates out_t types to receive the results from the workers
SUBROUTINE master_outtypes (comm,met,canopy,ssnow,rad,bal,air,soil,veg)

  USE mpi

  USE cable_def_types_mod

  IMPLICIT NONE

  INTEGER :: comm ! MPI communicator to talk to the workers

  TYPE(met_type), INTENT(IN) :: met
  TYPE(canopy_type), INTENT(IN) :: canopy
  TYPE(soil_snow_type), INTENT(IN) :: ssnow
  TYPE(radiation_type), INTENT(IN) :: rad
  TYPE (balances_type),INTENT(INOUT):: bal 
  TYPE (air_type),INTENT(IN)        :: air
  TYPE (soil_parameter_type),INTENT(IN) :: soil ! soil parameters
  TYPE (veg_parameter_type),INTENT(IN) :: veg ! vegetation parameters

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
  ALLOCATE (m3d_t(n3d, wnp))
  ALLOCATE (mat_t(nmat, wnp))
  ALLOCATE (vec_t(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = n3d + nmat + nvec
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  ! MPI: allocate vector to hold handles for the combined type
  ALLOCATE (recv_ts(wnp))

  ! MPI: type_struct that combines all data for a single worker

  ! MPI: allocate address vectors for 3D matrices
  ! currently only 1 is 3D but to make future expansion easier
  ! let's handle it the same way as 3D and 2D
  ALLOCATE(m3daddr(n3d))

  ! MPI: allocate address vectors for matrices
  ALLOCATE (maddr(nmat))
  ! MPI: allocate address vectors for vectors
  ALLOCATE (vaddr(nvec))
  ! MPI: allocate vector block lengths 
  ALLOCATE (blen(nvec))

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

     ! ssnow 2D
     midx = midx + 1
     CALL MPI_Get_address (ssnow%dtmlt(off,1), maddr(midx), ierr)
     CALL MPI_Type_create_hvector (3, r1len, r1stride, MPI_BYTE, &
     &                             mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

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
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)
     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%cexpkdm(off,1), maddr(midx), ierr) ! 27
     CALL MPI_Type_create_hvector (nrb, r1len, r1stride, MPI_BYTE, &
          &                        mat_t(midx, rank), ierr)
     CALL MPI_Type_commit (mat_t(midx, rank), ierr)

     midx = midx + 1
     ! REAL(r_1)
     CALL MPI_Get_address (rad%rhocbm(off,1), maddr(midx), ierr) ! 27
     CALL MPI_Type_create_hvector (swb, r1len, r1stride, MPI_BYTE, &
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
        WRITE (*,*) 'master: outtype invalid nmat ',midx,' constant, fix it!'
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
!     ! REAL(r_2)
!     CALL MPI_Get_address (canopy%potev_c(off), vaddr(vidx), ierr) ! 38
!     blen(vidx) = cnt * extr2
!     vidx = vidx + 1
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
     blen(vidx) = cnt * extid
     vidx = vidx + 1
     ! INTEGER(i_d)
     CALL MPI_Get_address (soil%ibp2(off), vaddr(vidx), ierr) ! 133
     blen(vidx) = cnt * extid
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
     blen(vidx) = cnt * extid
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
     
     ! MPI: sanity check
     IF (vidx /= nvec) THEN
        WRITE (*,*) 'master: outtype invalid nvec ',vidx,' constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     ! MPI: all vectors into a single hindexed type for each worker
     CALL MPI_Type_create_hindexed (vidx, blen, vaddr, MPI_BYTE, vec_t(rank), &
          &                         ierr)
     CALL MPI_Type_commit (vec_t(rank), ierr)

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

     WRITE (*,*) 'master: data recv from ',rank,': size, extent, lb: ', &
   &       tsize,text,tmplb

     totalrecv = totalrecv + tsize

  END DO

  WRITE (*,*) 'total data size received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total data size sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error master: totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE (blen)
  DEALLOCATE (vaddr)
  DEALLOCATE (maddr)
  DEALLOCATE(m3daddr)
  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN
END SUBROUTINE master_outtypes

! MPI: creates handles for receiving casa final results from the workers
SUBROUTINE master_casa_types (comm, casapool, casaflux, &
                              casamet, casabal, phen)

  USE mpi

  USE cable_def_types_mod
  USE casadimension
  USE casavariable
  USE phenvariable

  IMPLICIT NONE

  INTEGER :: comm ! MPI communicator to talk to the workers

  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  TYPE (phen_variable), INTENT(INOUT)  :: phen 

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len, ilen
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr

  ALLOCATE (casa_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = ncasa_mat + ncasa_vec
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
     ilen = cnt * extid

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

     last2d = bidx

     ! ------------- 1D vectors -------------

     bidx = bidx + 1
     CALL MPI_Get_address (casamet%glai(off), displs(bidx), ierr)
     blocks(bidx) = r2len

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

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'master: invalid number of casa fields, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, casa_ts(rank), ierr)
     CALL MPI_Type_commit (casa_ts(rank), ierr)

     CALL MPI_Type_size (casa_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (casa_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'casa results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     ! TODO: also free partial types for intypes, outtypes etc.
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of casa results received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of casa results sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: casa results totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN

END SUBROUTINE master_casa_types

! MPI: creates datatype handles to receive restart data from workers
SUBROUTINE master_restart_types (comm, canopy, air)

  USE mpi

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

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

  INTEGER :: tsize, totalrecv, totalsend
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr

  ALLOCATE (restart_ts(wnp))

  ! MPI: allocate temp vectors used for marshalling
  ntyp = nrestart
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

     bidx = 0

     ! ------------- 2D arrays -------------

!     bidx = bidx + 1
!     CALL MPI_Get_address (canopy%rwater(off,1), displs(bidx), ierr) ! 1
!     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
!     &                             types(bidx), ierr)
!     blocks(bidx) = 1

     bidx = bidx + 1
     CALL MPI_Get_address (canopy%evapfbl(off,1), displs(bidx), ierr) ! 2
     ! MPI: gol124: changed to r1 when Bernard ported to CABLE_r491
     CALL MPI_Type_create_hvector (ms, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
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
        WRITE (*,*) 'master: invalid number of restart fields, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, restart_ts(rank), ierr)
     CALL MPI_Type_commit (restart_ts(rank), ierr)

     CALL MPI_Type_size (restart_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (restart_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  RETURN
END SUBROUTINE master_restart_types

! MPI: scatters input data for timestep ktau to all workers
SUBROUTINE master_send_input (comm, dtypes, ktau)

  USE mpi

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: comm
  INTEGER, DIMENSION(:), INTENT(IN) :: dtypes
  INTEGER, INTENT(IN) :: ktau    ! timestep

  INTEGER :: rank, ierr

!  IF (.NOT. ALLOCATED(inp_req)) THEN
!     ALLOCATE (inp_req(wnp))
!  END IF

  DO rank = 1, wnp
     CALL MPI_Isend (MPI_BOTTOM, 1, dtypes(rank), rank, ktau, comm, &
     &               inp_req(rank), ierr)
  END DO

  !IF (.NOT. ALLOCATED(inp_stats)) THEN
  !   ALLOCATE (inp_stats(MPI_STATUS_SIZE, wnp))
  !END IF

  !CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  RETURN

END SUBROUTINE master_send_input

! receives model output variables from the workers for a single timestep
! uses the timestep value as the message tag
SUBROUTINE master_receive(comm, ktau, types)

  USE mpi

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

  DO rank = 1, wnp
     CALL MPI_Irecv (MPI_BOTTOM, 1, types(rank), rank, ktau, comm, &
     &               recv_req(rank), ierr)
  END DO

  ! MPI: we need all timesteps data before processing/saving
  !CALL MPI_Waitall (wnp, recv_req, recv_stats, ierr)

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
SUBROUTINE master_end (icycle, restart)

  USE mpi

  IMPLICIT NONE

  INTEGER :: icycle ! casa flag
  LOGICAL :: restart ! restart flag

  INTEGER :: rank, i, ierr

  ! MPI: free MPI types
  DO rank = 1, wnp

     CALL MPI_Type_free (inp_ts(rank), ierr)

     CALL MPI_Type_free (recv_ts(rank), ierr)

     IF(restart) THEN
        CALL MPI_Type_free (restart_ts(rank), ierr)
     END IF

     IF (icycle>0) THEN
        CALL MPI_Type_free (casa_ts(rank), ierr)
     END IF

     ! gol124: TODO: m3d_t, mat_t and vec_t can
     ! be freed at the end of master_outtypes
     DO i = 1, n3d
        CALL MPI_Type_free (m3d_t(i, rank), ierr)
     END DO

     DO i = 1, nmat
        CALL MPI_Type_free (mat_t(i, rank), ierr)
     END DO

     CALL MPI_Type_free (vec_t(rank), ierr)

  END DO

  ! MPI: free arrays for receive requests and statuses
  IF (ALLOCATED(recv_stats)) THEN
     DEALLOCATE (recv_stats)
  END IF
  IF (ALLOCATED(recv_req)) THEN
     DEALLOCATE (recv_req)
  END IF

  ! MPI: free arrays for send requests and statuses
  IF (ALLOCATED(inp_stats)) THEN
     DEALLOCATE (inp_stats)
  END IF
  IF (ALLOCATED(inp_req)) THEN
     DEALLOCATE (inp_req)
  END IF


  ! MPI: free derived datatype handle array
  DEALLOCATE (inp_ts)

  DEALLOCATE (recv_ts)

  IF(restart) THEN
     DEALLOCATE (restart_ts)
  END IF

  IF (icycle>0) THEN
     DEALLOCATE (casa_ts)
  END IF

  ! MPI: free partial derived datatype handle arrays
  DEALLOCATE (vec_t)
  DEALLOCATE (mat_t)
  DEALLOCATE (m3d_t)

  ! MPI: free landpoint decomposition info
  DEALLOCATE (wland)

  RETURN

END SUBROUTINE master_end

END MODULE cable_mpimaster

