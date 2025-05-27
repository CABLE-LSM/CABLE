! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation 
! (CSIRO) ABN 41 687 119 230.

MODULE cable_driver_common_mod
  !! Module for CABLE offline driver common routines.
  USE cable_common_module, ONLY : &
    filename,                     &
    calcsoilalbedo,               &
    redistrb,                     &
    wiltParam,                    &
    satuParam,                    &
    snmin,                        &
    cable_user,                   &
    gw_params,                    &
    cable_runtime
  USE cable_IO_vars_module, ONLY : &
    soilparmnew,                   &
    output,                        &
    patchout,                      &
    check,                         &
    verbose,                       &
    leaps,                         &
    logn,                          &
    fixedCO2,                      &
    ncciy,                         &
    gswpfile,                      &
    globalMetfile,                 &
    set_group_output_values,       &
    timeunits,                     &
    exists,                        &
    calendar
  USE casadimension, ONLY : icycle
  USE casavariable, ONLY : casafile
  USE cable_namelist_util, ONLY : &
    get_namelist_file_name,       &
    CABLE_NAMELIST,               &
    arg_not_namelist
  USE cable_mpi_mod, ONLY : mpi_grp_t
  USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ
  USE cable_input_module, ONLY : open_met_file
  USE CABLE_PLUME_MIP, ONLY : PLUME_MIP_TYPE, PLUME_MIP_INIT
  USE CABLE_CRU, ONLY : CRU_TYPE, CRU_INIT
  USE CABLE_site, ONLY : site_TYPE, site_INIT
USE read_namelists_mod_cbl, ONLY: read_cable_namelist

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: CASAONLY_ICYCLE_MIN = 10
  INTEGER, PARAMETER :: N_MET_FORCING_VARIABLES_GSWP = 8
    !! Number of GSWP met forcing variables (rain, snow, lw, sw, ps, qa, ta, wd)
  REAL, PARAMETER :: CONSISTENCY_CHECK_TOLERANCE = 1.e-7
    !! Max tolerance value for quasi-bit reproducibility checks

  LOGICAL, SAVE, PUBLIC :: vegparmnew    = .FALSE. ! using new format input file (BP dec 2007)
  LOGICAL, SAVE, PUBLIC :: spinup        = .FALSE. ! model spinup to soil state equilibrium?
  LOGICAL, SAVE, PUBLIC :: spincasa      = .FALSE. ! TRUE: CASA-CNP Will spin mloop times, FALSE: no spin up
  LOGICAL, SAVE, PUBLIC :: CASAONLY      = .FALSE. ! ONLY Run CASA-CNP
  LOGICAL, SAVE, PUBLIC :: l_casacnp     = .FALSE. ! using CASA-CNP with CABLE
  LOGICAL, SAVE, PUBLIC :: l_landuse     = .FALSE. ! using CASA-CNP with CABLE
  LOGICAL, SAVE, PUBLIC :: l_laiFeedbk   = .FALSE. ! using prognostic LAI
  LOGICAL, SAVE, PUBLIC :: l_vcmaxFeedbk = .FALSE. ! using prognostic Vcmax

  REAL, SAVE, PUBLIC :: delsoilM ! allowed variation in soil moisture for spin up
  REAL, SAVE, PUBLIC :: delsoilT ! allowed variation in soil temperature for spin up
  REAL, SAVE, PUBLIC :: delgwM = 1e-4

INTEGER, SAVE, PUBLIC :: LALLOC = 0            ! alloc coeff passed to spincasa
INTEGER, PARAMETER    :: nmlunitnumber = 88813 ! this is to satisfy UM method 
                                               ! where a shared_unitnumber is used
  PUBLIC :: cable_driver_init
  PUBLIC :: cable_driver_init_gswp
  PUBLIC :: cable_driver_init_plume
  PUBLIC :: cable_driver_init_cru
  PUBLIC :: cable_driver_init_site
  PUBLIC :: cable_driver_init_default
  PUBLIC :: prepareFiles
  PUBLIC :: renameFiles
  PUBLIC :: prepareFiles_princeton
  PUBLIC :: LUCdriver
  PUBLIC :: compare_consistency_check_values

CONTAINS

  SUBROUTINE cable_driver_init(mpi_grp, NRRRR)
    !! Model initialisation routine for the CABLE offline driver.
    TYPE(mpi_grp_t), INTENT(IN) :: mpi_grp !! MPI group to use
    INTEGER, INTENT(OUT) :: NRRRR !! Number of repeated spin-up cycles

    INTEGER :: ioerror, unit
    CHARACTER(len=4) :: cRank ! for worker-logfiles

    !check to see if first argument passed to cable is
    !the name of the namelist file
    !if not use cable.nml
    CALL get_namelist_file_name()

    IF (mpi_grp%rank == 0) THEN
      WRITE(*,*) "THE NAME LIST IS ", CABLE_NAMELIST
    END IF

    CALL read_cable_namelist( nmlunitnumber, CABLE_NAMELIST, vegparmnew,      &
                                             spinup    ,      &
                                             spincasa  ,      &
                                             CASAONLY  ,      &
                                             l_casacnp ,      &
                                             l_landuse ,      &
                                             l_laiFeedbk ,    &
                                             l_vcmaxFeedbk    )

    cable_runtime%offline = .TRUE.

    ! Open log file:
    IF (mpi_grp%rank == 0) THEN
      OPEN(logn, FILE=filename%log)
    ELSE IF (cable_user%logworker) THEN
      WRITE(cRank, FMT='(I4.4)') mpi_grp%rank
      OPEN(NEWUNIT=logn, FILE="cable_log_"//cRank, STATUS="REPLACE")
    ELSE
      OPEN(NEWUNIT=logn, FILE="/dev/null")
    END IF

    IF (IARGC() > 0 .AND. arg_not_namelist) THEN
      CALL GETARG(1, filename%met)
      CALL GETARG(2, casafile%cnpipool)
    END IF

    ! Initialise flags to output individual variables according to group
    ! options from the namelist file
    IF (mpi_grp%rank == 0) THEN
      CALL set_group_output_values()
    END IF

    IF (TRIM(cable_user%POPLUC_RunType) == 'static') THEN
      cable_user%POPLUC= .FALSE.
    END IF

    ! TODO(Sean): we should not be setting namelist parameters in the following if
    ! block - all options are all configurable via the namelist file and is
    ! unclear that these options are being overwritten. A better approach would be
    ! to error for bad combinations of namelist parameters.
    IF (icycle > CASAONLY_ICYCLE_MIN) THEN
      icycle                     = icycle - CASAONLY_ICYCLE_MIN
      CASAONLY                   = .TRUE.
      CABLE_USER%CASA_DUMP_READ  = .TRUE.
      CABLE_USER%CASA_DUMP_WRITE = .FALSE.
    ELSE IF (icycle == 0) THEN
      CABLE_USER%CASA_DUMP_READ  = .FALSE.
      spincasa                   = .FALSE.
      CABLE_USER%CALL_POP        = .FALSE.
    END IF

    ! TODO(Sean): overwriting l_casacnp defeats the purpose of it being a namelist
    ! option - we should either deprecate the l_casacnp option or not overwrite
    ! its value.
    l_casacnp = icycle > 0

    IF (l_casacnp .AND. (icycle == 0 .OR. icycle > 3)) THEN
      STOP 'icycle must be 1 to 3 when using casaCNP'
    END IF
    IF ((l_laiFeedbk .OR. l_vcmaxFeedbk) .AND. (.NOT. l_casacnp)) THEN
      STOP 'casaCNP required to get prognostic LAI or Vcmax'
    END IF
    IF (l_vcmaxFeedbk .AND. (icycle < 2 .OR. icycle > 3)) THEN
      STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
    END IF
    IF (icycle > 0 .AND. (.NOT. soilparmnew)) THEN
      STOP 'casaCNP must use new soil parameters'
    END IF

    ! vh_js ! suggest LALLOC should ulitmately be a switch in the .nml file
    IF (cable_user%CALL_POP) THEN
      LALLOC = 3 ! for use with POP: makes use of pipe model to partition between stem and leaf
    END IF

    NRRRR = MERGE(MAX(cable_user%CASA_NREP,1), 1, CASAONLY)

  END SUBROUTINE cable_driver_init

  SUBROUTINE cable_driver_init_gswp(mpi_grp, GSWP_MID, NRRRR)
    !! Model initialisation routine (GSWP specific).
    TYPE(mpi_grp_t), INTENT(IN) :: mpi_grp !! MPI group to use
    INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL :: GSWP_MID(:,:) !! NetCDF file IDs for GSWP met forcing
    INTEGER, INTENT(IN), OPTIONAL :: NRRRR !! Number of repeated spin-up cycles

    IF (cable_user%YearStart == 0) THEN
      IF (ncciy == 0) THEN
        IF (mpi_grp%rank == 0) THEN
          PRINT*, 'undefined start year for gswp met: '
          PRINT*, 'enter value for ncciy or'
          PRINT*, '(CABLE_USER%YearStart and  CABLE_USER%YearEnd) in cable.nml'
        END IF
        WRITE(logn,*) 'undefined start year for gswp met: '
        WRITE(logn,*) 'enter value for ncciy or'
        WRITE(logn,*) '(CABLE_USER%YearStart and  CABLE_USER%YearEnd) in cable.nml'
        STOP
      END IF
      cable_user%YearStart = ncciy
      cable_user%YearEnd = ncciy
    END IF

    IF (.NOT. (PRESENT(GSWP_MID) .AND. PRESENT(NRRRR))) RETURN

    IF (NRRRR > 1 .AND. (.NOT. ALLOCATED(GSWP_MID))) THEN
      ALLOCATE(GSWP_MID(N_MET_FORCING_VARIABLES_GSWP, cable_user%YearStart:cable_user%YearEnd))
    END IF

  END SUBROUTINE cable_driver_init_gswp

  SUBROUTINE cable_driver_init_site(site)
    !* Model initialisation routine (site met specific).
    ! Site experiment, e.g. AmazonFace (spinup or transient run type).
    TYPE (site_TYPE), INTENT(OUT) :: site

    CHARACTER(len=9) :: str1, str2, str3

    IF (.NOT. l_casacnp) THEN
      WRITE(*,*) "MetType=site only works with CASA-CNP turned on"
      STOP 991
    END IF

    CALL site_INIT( site )
    WRITE(str1,'(i4)') cable_user%YearStart
    str1 = ADJUSTL(str1)
    WRITE(str2,'(i2)') 1
    str2 = ADJUSTL(str2)
    WRITE(str3,'(i2)') 1
    str3 = ADJUSTL(str3)
    timeunits="seconds since "//TRIM(str1)//"-"//TRIM(str2)//"-"//TRIM(str3)//"00:00"
    calendar = 'standard'

  END SUBROUTINE cable_driver_init_site

  SUBROUTINE cable_driver_init_default(dels, koffset, kend)
    !! Model initialisation routine (default met specific).
    REAL, INTENT(OUT) :: dels !! Time step size in seconds
    INTEGER, INTENT(OUT) :: koffset !! Timestep to start at
    INTEGER, INTENT(OUT) :: kend !! No. of time steps in run

    ! Open met data and get site information from netcdf file.
    ! This retrieves time step size, number of timesteps, starting date,
    ! latitudes, longitudes, number of sites.
    CALL open_met_file(dels, koffset, kend, spinup, CTFRZ)
    IF (koffset /= 0 .AND. cable_user%CALL_POP) THEN
      WRITE(*,*) "When using POP, episode must start at Jan 1st!"
      STOP 991
    END IF

  END SUBROUTINE cable_driver_init_default

  SUBROUTINE cable_driver_init_plume(dels, koffset, PLUME)
    !* Model initialisation routine (PLUME specific).
    ! PLUME experiment setup using WATCH.
    REAL, INTENT(OUT) :: dels !! Time step size in seconds
    INTEGER, INTENT(OUT) :: koffset !! Timestep to start at
    TYPE(PLUME_MIP_TYPE), INTENT(OUT) :: PLUME

    CHARACTER(len=9) :: str1, str2, str3

    CALL PLUME_MIP_INIT(PLUME)
    dels = PLUME%dt
    koffset = 0
    leaps = PLUME%LeapYears
    WRITE(str1,'(i4)') cable_user%YearStart
    str1 = ADJUSTL(str1)
    WRITE(str2,'(i2)') 1
    str2 = ADJUSTL(str2)
    WRITE(str3,'(i2)') 1
    str3 = ADJUSTL(str3)
    timeunits="seconds since "//TRIM(str1)//"-"//TRIM(str2)//"-"//TRIM(str3)//"00:00"

  END SUBROUTINE cable_driver_init_plume

  SUBROUTINE cable_driver_init_cru(dels, koffset, CRU)
    !* Model initialisation routine (CRU specific).
    ! TRENDY experiment using CRU-NCEP.
    REAL, INTENT(OUT) :: dels !! Time step size in seconds
    INTEGER, INTENT(OUT) :: koffset !! Timestep to start at
    TYPE(CRU_TYPE), INTENT(OUT) :: CRU

    CHARACTER(len=9) :: str1, str2, str3

    CALL CRU_INIT(CRU)
    dels = CRU%dtsecs
    koffset = 0
    leaps = .FALSE. ! No leap years in CRU-NCEP
    exists%Snowf = .FALSE.
      ! No snow in CRU-NCEP, so ensure it will be determined from temperature
      ! in CABLE
    WRITE(str1,'(i4)') cable_user%YearStart
    str1 = ADJUSTL(str1)
    WRITE(str2,'(i2)') 1
    str2 = ADJUSTL(str2)
    WRITE(str3,'(i2)') 1
    str3 = ADJUSTL(str3)
    timeunits="seconds since "//TRIM(str1)//"-"//TRIM(str2)//"-"//TRIM(str3)//"00:00"

  END SUBROUTINE cable_driver_init_cru

  SUBROUTINE prepareFiles(ncciy)
    !* Select the correct files given the year for filenames following 
    ! the gswp format
    
    INTEGER, INTENT(IN) :: ncciy !! Year to select met. forcing data.

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

  SUBROUTINE renameFiles(logn,inFile,ncciy,inName)
    !* Replace the year in the filename with the value of ncciy 
    ! for the gswp file format.

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logn !! Log file unit number
    INTEGER, INTENT(IN) :: ncciy !! Year to use in replacement in filenames
    INTEGER:: nn
    CHARACTER(LEN=200), INTENT(INOUT) :: inFile
    CHARACTER(LEN=*),  INTENT(IN)    :: inName
    INTEGER :: idummy

    nn = INDEX(inFile,'19')
    READ(inFile(nn:nn+3),'(i4)') idummy
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)

  END SUBROUTINE renameFiles

  SUBROUTINE prepareFiles_princeton(ncciy)
    !* Select the correct files given the year for filenames following the 
    ! princeton format
    INTEGER, INTENT(IN) :: ncciy
  
    WRITE(logn,*) 'CABLE offline global run using princeton forcing for ', ncciy
    PRINT *,     'CABLE offline global run using princeton forcing for ', ncciy
  
    CALL renameFiles_princeton(logn,gswpfile%rainf,ncciy,'rainf')
    CALL renameFiles_princeton(logn,gswpfile%LWdown,ncciy,'LWdown')
    CALL renameFiles_princeton(logn,gswpfile%SWdown,ncciy,'SWdown')
    CALL renameFiles_princeton(logn,gswpfile%PSurf,ncciy,'PSurf')
    CALL renameFiles_princeton(logn,gswpfile%Qair,ncciy,'Qair')
    CALL renameFiles_princeton(logn,gswpfile%Tair,ncciy,'Tair')
    CALL renameFiles_princeton(logn,gswpfile%wind,ncciy,'wind')
 
  END SUBROUTINE prepareFiles_princeton
 
  SUBROUTINE renameFiles_princeton(logn,inFile,ncciy,inName)
    !* Replace the year in the filename with the value of ncciy for 
    ! the princeton format
    INTEGER, INTENT(IN) :: logn,ncciy
    INTEGER:: nn
    CHARACTER(LEN=200), INTENT(INOUT) :: inFile
    CHARACTER(LEN=*),  INTENT(IN)        :: inName
    INTEGER :: idummy
  
    nn = INDEX(inFile,'19')
    READ(inFile(nn:nn+3),'(i4)') idummy
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    nn = INDEX(inFile,'19', BACK=.TRUE.)
    READ(inFile(nn:nn+3),'(i4)') idummy
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    READ(inFile(nn-5:nn-2),'(i4)') idummy
    WRITE(inFile(nn-5:nn-2),'(i4.4)') ncciy
    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
 
  END SUBROUTINE renameFiles_princeton

  !==============================================================================
  ! subroutine for 
  SUBROUTINE LUCdriver( casabiome,casapool, &
    casaflux,POP,LUC_EXPT, POPLUC, veg )
      !* Reading LU input data, zeroing biomass in empty secondary forest tiles
      ! and tranferring LUC-based age weights for secondary forest to POP structure
    USE cable_def_types_mod , ONLY: veg_parameter_type, mland
    USE cable_carbon_module
    USE cable_common_module, ONLY: CABLE_USER, CurYear
    USE casa_ncdf_module, ONLY: is_casa_time
    USE cable_IO_vars_module, ONLY: logn, landpt, patch
    USE casadimension
    USE casaparm
    USE casavariable
    USE POP_Types,  ONLY: POP_TYPE
    USE POPMODULE,            ONLY: POPStep, POP_init_single
    USE TypeDef,              ONLY: i4b, dp
    USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE, read_LUH2,&
      ptos,ptog,stog,gtos,grassfrac, pharv, smharv, syharv
    USE POPLUC_Types
    USE POPLUC_Module, ONLY: POPLUCStep, POPLUC_weights_Transfer, WRITE_LUC_OUTPUT_NC, &
      POP_LUC_CASA_transfer,  WRITE_LUC_RESTART_NC, READ_LUC_RESTART_NC

    IMPLICIT NONE

    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (POP_TYPE), INTENT(INOUT)     :: POP
    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    TYPE (veg_parameter_type),    INTENT(IN) :: veg  ! vegetation parameters

    INTEGER ::  k, j, l, yyyy

    WRITE(*,*) 'cablecasa_LUC', CurYear
    yyyy = CurYear

    LUC_EXPT%CTSTEP = yyyy -  LUC_EXPT%FirstYear + 1

    CALL READ_LUH2(LUC_EXPT)

    DO k=1,mland
      POPLUC%ptos(k) = LUC_EXPT%INPUT(ptos)%VAL(k)
      POPLUC%ptog(k) = LUC_EXPT%INPUT(ptog)%VAL(k)
      POPLUC%stop(k) = 0.0
      POPLUC%stog(k) = LUC_EXPT%INPUT(stog)%VAL(k)
      POPLUC%gtop(k) = 0.0
      POPLUC%gtos(k) = LUC_EXPT%INPUT(gtos)%VAL(k)
      POPLUC%pharv(k) = LUC_EXPT%INPUT(pharv)%VAL(k)
      POPLUC%smharv(k) = LUC_EXPT%INPUT(smharv)%VAL(k)
      POPLUC%syharv(k) = LUC_EXPT%INPUT(syharv)%VAL(k)
      POPLUC%thisyear = yyyy
    ENDDO

    ! zero secondary forest tiles in POP where secondary forest area is zero
    DO k=1,mland
      IF ((POPLUC%frac_primf(k)-POPLUC%frac_forest(k))==0.0 &
           .AND. (.NOT.LUC_EXPT%prim_only(k))) THEN
        j = landpt(k)%cstart+1
        DO l=1,SIZE(POP%Iwood)
          IF( POP%Iwood(l) == j) THEN
            CALL POP_init_single(POP,veg%disturbance_interval,l)
            EXIT
          ENDIF
        ENDDO

        casapool%cplant(j,leaf) = 0.01
        casapool%nplant(j,leaf)= casabiome%ratioNCplantmin(veg%iveg(j),leaf)* casapool%cplant(j,leaf)
        casapool%pplant(j,leaf)= casabiome%ratioPCplantmin(veg%iveg(j),leaf)* casapool%cplant(j,leaf)

        casapool%cplant(j,froot) = 0.01
        casapool%nplant(j,froot)= casabiome%ratioNCplantmin(veg%iveg(j),froot)* casapool%cplant(j,froot)
        casapool%pplant(j,froot)= casabiome%ratioPCplantmin(veg%iveg(j),froot)* casapool%cplant(j,froot)

        casapool%cplant(j,wood) = 0.01
        casapool%nplant(j,wood)= casabiome%ratioNCplantmin(veg%iveg(j),wood)* casapool%cplant(j,wood)
        casapool%pplant(j,wood)= casabiome%ratioPCplantmin(veg%iveg(j),wood)* casapool%cplant(j,wood)
        casaflux%frac_sapwood(j) = 1.0

      ENDIF
    ENDDO

    CALL POPLUCStep(POPLUC,yyyy)

    CALL POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

  END SUBROUTINE LUCdriver

  SUBROUTINE compare_consistency_check_values(new_sumbal)
    !* Compare reference values for quasi-bitwise reproducibility checks and write
    ! reference value on failed reproducibility.
    DOUBLE PRECISION, INTENT(IN) :: new_sumbal !! Reference value for current run.

    DOUBLE PRECISION :: trunk_sumbal
    INTEGER :: unit, ioerror

    OPEN(newunit=unit, file=filename%trunk_sumbal, status='OLD', action='READ', iostat=ioerror)
    IF(ioerror == 0) READ(unit, *) trunk_sumbal
    CLOSE(unit)

    IF(ioerror == 0 .AND. ABS(new_sumbal - trunk_sumbal) < CONSISTENCY_CHECK_TOLERANCE) THEN
      PRINT *, "Internal check shows this version reproduces the trunk sumbal"
      RETURN
    END IF

    PRINT *, "Internal check shows in this version new_sumbal != trunk sumbal"
    PRINT *, "Writing new_sumbal to the file:", TRIM(filename%new_sumbal)
    OPEN(newunit=unit, file=filename%new_sumbal)
    WRITE(unit, '(F20.7)') new_sumbal
    CLOSE(unit)

  END SUBROUTINE compare_consistency_check_values

END MODULE cable_driver_common_mod
