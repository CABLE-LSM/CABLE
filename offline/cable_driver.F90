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
!
! Contact: Bernard.Pak@csiro.au
!
! History: Since 1.4b, capability to run global offline (ncciy = YEAR),
!          inclusion of call to CASA-CNP (icycle>0)
!          soil_snow_type now ssnow (instead of ssoil)
!
!
! ==============================================================================
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
   USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps,                  &
                                   verbose, fixedCO2,output,check,patchout,    &
                                   patch_type,soilparmnew
   USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user,     &
                                   cable_runtime, filename, redistrb,          & 
                                   report_version_no, wiltParam, satuParam,    &
                                   calcsoilalbedo
   USE cable_data_module,    ONLY: driver_type, point2constants
   USE cable_input_module,   ONLY: open_met_file,load_parameters,              &
                                   get_met_data,close_met_file
   USE cable_output_module,  ONLY: create_restart,open_output_file,            &
                                   write_output,close_output_file
   USE cable_cbm_module
   
   USE cable_diag_module
   
   ! modules related to CASA-CNP
   USE casadimension,       ONLY: icycle 
   USE casavariable,        ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
                                  casa_met, casa_balance
   USE phenvariable,        ONLY: phen_variable
  USE casa_inout_mod

use bgcdriver_mod, ONLY : bgcdriver
   IMPLICIT NONE
   
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

   !___ unique unit/file identifiers for cable_diag: arbitrarily 5 here 
   INTEGER, SAVE :: iDiagZero=0, iDiag1=0, iDiag2=0, iDiag3=0, iDiag4=0

   ! switches etc defined thru namelist (by default cable.nml)
   NAMELIST/CABLE/                  &
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

   ! Vars for standard for quasi-bitwise reproducability b/n runs
   ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
   CHARACTER(len=30), PARAMETER ::                                             &
      Ftrunk_sumbal  = ".trunk_sumbal",                                        &
      Fnew_sumbal    = "new_sumbal"

   DOUBLE PRECISION ::                                                                     &
      trunk_sumbal = 0.0, & !
      new_sumbal = 0.0

   INTEGER :: nkend=0
   INTEGER :: ioerror

   ! END header

   ! Open, read and close the namelist file.
   OPEN( 10, FILE = CABLE_NAMELIST )
      READ( 10, NML=CABLE )   !where NML=CABLE defined above
   CLOSE(10)

   ! Open, read and close the consistency check file.
   ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
   IF(cable_user%consistency_check) THEN 
      OPEN( 11, FILE = Ftrunk_sumbal,STATUS='old',ACTION='READ',IOSTAT=ioerror )
         IF(ioerror==0) then
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
         STOP 'Please check input in namelist file.'
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
         CALL get_met_data( spinup, spinConv, met, soil,                    &
                            rad, veg, kend, dels, C%TFRZ, ktau ) 
   
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
            call bgcdriver( ktau, kstart, kend, dels, met,                     &
                            ssnow, canopy, veg, soil, casabiome,               &
                            casapool, casaflux, casamet, casabal,              &
                            phen, spinConv, spinup, ktauday, idoy,             &
                            .FALSE., .FALSE. )
         ENDIF 
   
         ! sumcflux is pulled out of subroutine cbm
         ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
         CALL sumcflux( ktau, kstart, kend, dels, bgc,                         &
                        canopy, soil, ssnow, sum_flux, veg,                    &
                        met, casaflux, l_vcmaxFeedbk )
   
         ! Write time step's output to file if either: we're not spinning up 
         ! or we're spinning up and the spinup has converged:
         IF((.NOT.spinup).OR.(spinup.AND.spinConv))                            &
            CALL write_output( dels, ktau, met, canopy, ssnow,                 &
                               rad, bal, air, soil, veg, C%SBOLTZ,             &
                               C%EMLEAF, C%EMSOIL )
   
         ! dump bitwise reproducible testing data
         IF( cable_user%RUN_DIAG_LEVEL == 'zero') THEN
            IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
               call cable_diag( iDiagZero, "FLUXES", mp, kend, ktau,                   &
                                knode_gl, "FLUXES",                            &
                          canopy%fe + canopy%fh )
         ENDIF
   
         ! Check this run against standard for quasi-bitwise reproducability
         ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
         IF(cable_user%consistency_check) THEN 
            
            new_sumbal = new_sumbal + SUM(bal%wbal_tot) + SUM(bal%ebal_tot)          &
                             + SUM(bal%ebal_tot_cncheck)
  
            IF( ktau == kend ) THEN
               nkend = nkend+1

               IF( abs(new_sumbal-trunk_sumbal) < 1.e-7) THEN

                  print *, ""
                  print *, &
                  "NB. Offline-serial runs spinup cycles:", nkend
                  print *, &
                  "Internal check shows this version reproduces the trunk sumbal"
               
               ELSE

                  print *, ""
                  print *, &
                  "NB. Offline-serial runs spinup cycles:", nkend
                  print *, &
                  "Internal check shows in this version new_sumbal != trunk sumbal"
                  print *, &
                  "Writing new_sumbal to the file:", TRIM(Fnew_sumbal)
                        
                  OPEN( 12, FILE = Fnew_sumbal )
                     WRITE( 12, '(F20.7)' ) new_sumbal  ! written by previous trunk version
                  CLOSE(12)
               
               ENDIF   
            ENDIF   
            
         ENDIF

         
      END DO ! END Do loop over timestep ktau



   
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
               maxdiff = MAXLOC(ABS(ssnow%wb-soilMtemp))
               PRINT *, 'Example location of moisture non-convergence: ', &
                        maxdiff
               PRINT *, 'ssnow%wb : ', ssnow%wb(maxdiff(1),maxdiff(2))
               PRINT *, 'soilMtemp: ', soilMtemp(maxdiff(1),maxdiff(2))
               maxdiff = MAXLOC(ABS(ssnow%tgg-soilTtemp))
               PRINT *, 'Example location of temperature non-convergence: ', &
                        maxdiff
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

      ELSE

         ! if not spinning up, or spin up has converged, exit:
         EXIT
       
      END IF

   END DO

   IF (icycle > 0) THEN
      
      CALL casa_poolout( ktau, veg, soil, casabiome,                           &
                         casapool, casaflux, casamet, casabal, phen )

      CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
  
   END IF

   ! Write restart file if requested:
   IF(output%restart)                                                          &
      CALL create_restart( logn, dels, ktau, soil, veg, ssnow,                 &
                           canopy, rough, rad, bgc, bal )
      
   ! Close met data input file:
   CALL close_met_file
 
   ! Close output file and deallocate main variables:
   CALL close_output_file( bal, air, bgc, canopy, met,                         &
                           rad, rough, soil, ssnow,                            &
                           sum_flux, veg )

   WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck
   
   ! Close log file
   CLOSE(logn)

END PROGRAM cable_offline_driver


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






