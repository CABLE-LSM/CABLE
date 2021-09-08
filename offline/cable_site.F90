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
! Purpose: site information to drive biophysical+biogeochemical site run
! Contact: vanessa.haverd@csiro.au
!
! History: July 2017
!
!
! ========================================untitled======================================
! Uses:		  CABLE_COMMON_MODULE
!		  cable_IO_vars_module
!
! CALLs:       site_init
!              site_GET_CO2_Ndep
! ==============================================================================

MODULE CABLE_site


  USE casa_ncdf_module, ONLY: HANDLE_ERR, GET_UNIT
  USE CABLE_COMMON_MODULE, ONLY: CurYear !  current year of multiannual run

  USE cable_IO_vars_module, ONLY: &  ! Selected cable_iovars.F90 variables:
       logn
  ! in CRU-NCEP. Setting this ensures snow will be determined in CABLE from temperature.

  IMPLICIT NONE

  TYPE site_TYPE
     CHARACTER(len=15)  :: RunType ! 'spinup', 'transient', 'AMB', 'ELE'
     REAL, DIMENSION(:) ,ALLOCATABLE :: CO2VALS  ! Global annual CO2 values (dim is the number of years of data, or 1 if time-invariant)
     REAL, DIMENSION(:) ,ALLOCATABLE :: NdepVALS  ! Global annual Ndep values (dim is the number of years of data, or 1 if time-invariant)
     REAL, DIMENSION(:) ,ALLOCATABLE :: PdepVALS  ! Global annual Pdep values (dim is the number of years of data, or 1 if time-invariant)
     INTEGER  :: mland                   ! Number of land cells
     CHARACTER(len=200) :: CO2NdepFile   ! CO2Ndepfile with path
     INTEGER :: spinstartyear
     INTEGER :: spinendyear
     REAL :: spinCO2 ! ppm in 1850
     REAL :: spinNdep  ! kgNha-1y-1 in 1850
     REAL :: spinPdep  ! kgPha-1y-1 in 1850
     REAL :: CO2   ! CO2 for current time step
     REAL :: Ndep  ! Ndep for current time step
     REAL :: Pdep  ! Pdep for current time step
  END TYPE site_TYPE

  TYPE (site_TYPE):: site  ! Define the variable CRU, of type CRU_TYPE

CONTAINS

  !**************************************************************************************************

  SUBROUTINE site_INIT( site )

    ! Initialise the contents of the site defined type collection, from the site namelist file

    !**************************************************************************************************


    IMPLICIT NONE

    TYPE (site_TYPE):: site

    INTEGER              :: nmlunit    ! Unit number for reading namelist file

    ! Temporary local names for site% variables as they are read from the namelist file.


    CHARACTER(len=15)    :: RunType
    CHARACTER(len=200)   :: CO2NdepFile
    INTEGER :: spinstartyear
    INTEGER :: spinendyear
    REAL :: spinCO2 ! ppm in 1850
    REAL :: spinNdep  ! kgNha-1y-1 in 1850
    REAL :: spinPdep  ! kgPha-1y-1 in 1850

    ! Flag for errors
    LOGICAL              :: ERR = .FALSE.

    NAMELIST /siteNML/ RunType, CO2NdepFile, spinstartyear, spinendyear, spinCO2, &
         spinNdep, spinPdep

    ! Read site namelist settings
    CALL GET_UNIT(nmlunit)  ! CABLE routine finds spare unit number
    OPEN (nmlunit,FILE="site.nml",STATUS='OLD',ACTION='READ')
    READ (nmlunit,NML=siteNML)
    CLOSE(nmlunit)

    ! Assign namelist settings to corresponding CRU defined-type elements
    site%RunType = RunType
    site%CO2NDepFile      = CO2NdepFile
    site%spinstartyear = spinstartyear
    site%spinendyear          = spinendyear
    site%spinCO2 = spinCO2
    site%spinNdep = spinNdep
    site%spinPdep = spinPdep
    ! Print settings
    WRITE(*   ,*)"========================================= SITE INFO  ============"
    WRITE(*   ,*)"site settings chosen:"
    WRITE(*   ,*)" RunType: ",TRIM( site%RunType)
    WRITE(*   ,*)" CO2NDepFile: ",TRIM(site%CO2NdepFile)
    WRITE(*   ,*)" spin start year               : ",site%spinstartyear
    WRITE(*   ,*)" spin end year : ",site%spinendyear
    WRITE(*   ,*)" CO2 value for spinup [ppm]  : ", site%spinCO2
    WRITE(*   ,*)" Ndep value for spinup [kg n ha-1 y-1] ", site%spinNdep
    WRITE(*   ,*)" Pdep value for spinup [kg P ha-1 y-1] ", site%spinPdep
    WRITE(logn,*)"========================================= SITE INFO ============"
    WRITE(logn,*)"site settings chosen:"
    WRITE(logn,*)" RunType: ",TRIM( site%RunType)
    WRITE(logn,*)" CO2NDepFile: ",TRIM(site%CO2NdepFile)
    WRITE(logn,*)" spin start year      : ", site%spinstartyear
    WRITE(logn,*)" spin end year : ",site%spinendyear
    WRITE(logn,*)" CO2 value for spinup [ppm]  : ", site%spinCO2
    WRITE(logn,*)" Ndep value for spinup [kg n ha-1 y-1] ", site%spinNdep
    WRITE(logn,*)" Pdep value for spinup [kg n ha-1 y-1] ", site%spinPdep


    WRITE(*   ,*)"========================================= site ============"
    WRITE(logn,*)"========================================= site ============"

    site%mland = 1



  END SUBROUTINE site_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE site_GET_CO2_Ndep( site)

    ! Get CO2 and N-dep values for use with a site run. Assign a static value if specified otherwise
    ! on the first call read all the annual values from a file into the site%CO2VALS and site%Ndep arrays

    IMPLICIT NONE

    TYPE(site_TYPE) :: site           ! site structure
    INTEGER              :: i, iunit, iyear, IOS = 0
    LOGICAL,        SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine

    ! For S0_TRENDY, use only static 1860 CO2 value and return immediately
    IF ( TRIM(site%RunType) .EQ. "spinup") THEN
       site%CO2 = site%spinCO2  ! CO2 in ppm for spinup
       site%Ndep = site%spinNdep
       site%Pdep = site%spinPdep
       RETURN

       ! If not spinup, varying CO2 and Ndep values will be used...
    ELSE

       ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
       ! values, open the (ascii) CO2 file and read the values into the array.
       IF (CALL1) THEN
          ALLOCATE( site%CO2VALS( 1850:2100 ) )
          ALLOCATE( site%NdepVALS( 1850:2100 ))
          ALLOCATE( site%PdepVALS( 1850:2100 ))
          CALL GET_UNIT(iunit)
          OPEN (iunit, FILE=TRIM(site%CO2NdepFILE), STATUS="OLD", ACTION="READ")
          ! get past header
          READ(iunit, *)
          DO WHILE( IOS .EQ. 0 )
             READ(iunit, FMT=* , IOSTAT=IOS) iyear, site%CO2VALS(iyear), site%NdepVALS(iyear), &
                  site%PdepVALS(iyear)
          END DO
          CLOSE(iunit)
          CALL1 = .FALSE.
       END IF
91     FORMAT (i4,',',f5.1,',',f4.2,',',f5.3)
       ! In all varying CO2 cases, return the element of the array for the current year
       ! as a single CO2 value.
       !



       site%CO2 = site%CO2VALS( CurYear )
       site%Ndep = site%NdepVALS( CurYear )
       site%Pdep = site%PdepVALS( CurYear )

    END IF

  END SUBROUTINE site_GET_CO2_Ndep

  !**************************************************************************************************

END MODULE CABLE_SITE
