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
! Purpose: Reads vegetation and soil parameter files, fills vegin, soilin
!          NB. Most soil parameters overwritten by spatially explicit datasets
!          input as ancillary file (for ACCESS) or surface data file (for offline)
!          Module enables accessibility of variables throughout CABLE
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: v2.0 vegin%dleaf now calculated from leaf length and width
!          Parameter files were read elsewhere in v1.8 (init_subrs)
!
!
! ==============================================================================

MODULE cable_common_module

  USE cable_pft_params_mod, ONLY : vegin
  USE cable_soil_params_mod, ONLY : soilin

  IMPLICIT NONE

  !---allows reference to "gl"obal timestep in run (from atm_step)
  !---total number of timesteps, and processing node
  INTEGER, SAVE :: ktau_gl, kend_gl, knode_gl, kwidth_gl

  LOGICAL :: L_fudge = .FALSE.

  INTEGER, SAVE :: CurYear  ! current year of multiannual run

  ! user switches turned on/off by the user thru namelists

  ! trunk modifications protected by these switches
  TYPE hide_switches
     LOGICAL ::                                                               &
                                ! L.Stevens - Test Switches
          L_NEW_ROUGHNESS_SOIL  = .FALSE., & ! from Ticket?
          L_NEW_RUNOFF_SPEED    = .FALSE., & ! from Ticket?
          L_NEW_REDUCE_SOILEVP  = .FALSE.    ! from Ticket?

  END TYPE hide_switches

  ! instantiate internal switches
  TYPE (hide_switches), SAVE :: hide


  ! set from environment variable $HOME
  CHARACTER(LEN=200) ::                                                       &
       myhome

  ! switch to calc sil albedo using soil colour - Ticket #27
  LOGICAL :: calcsoilalbedo = .FALSE.
  !---Lestevens Sept2012
  !---CASACNP switches and cycle index
  LOGICAL, SAVE :: l_casacnp,l_laiFeedbk,l_vcmaxFeedbk

  !---CABLE runtime switches def in this type
  TYPE kbl_internal_switches
     LOGICAL :: um = .FALSE., um_explicit = .FALSE., um_implicit = .FALSE.,   &
          um_radiation = .FALSE., um_hydrology = .FALSE.
     LOGICAL :: offline = .FALSE., mk3l = .FALSE.
  END TYPE kbl_internal_switches

  ! instantiate internal switches
  TYPE(kbl_internal_switches), SAVE :: cable_runtime

  ! user switches turned on/off by the user thru namelists
  ! CABLE-2.0 user switches all in single namelist file cable.nml
  ! clean these up for new namelist(s) format
  TYPE kbl_user_switches
     !jhan: this is redundant now we all use filename%veg?
     CHARACTER(LEN=200) ::                                                    &
          VEG_PARS_FILE  !

     CHARACTER(LEN=20) ::                                                     &
          FWSOIL_SWITCH, &     !
          PHENOLOGY_SWITCH = 'MODIS'   ! alternative is 'climate'
     !--- LN ------------------------------------------[

     ! Ticket #56
     CHARACTER(LEN=20) ::                                                     &
          GS_SWITCH='leuning'

     CHARACTER(LEN=10) :: RunIden       = 'STANDARD'  !
     CHARACTER(LEN=6)  :: MetType       = ' ' !
     CHARACTER(LEN=20) :: SOIL_STRUC    = "default" ! 'default' or 'sli'
     CHARACTER(LEN=3)  :: POP_out       = 'rst' ! POP output type ('epi' or 'rst')
     CHARACTER(LEN=50) :: POP_rst       = ' ' !
     CHARACTER(LEN=8)  :: CASA_OUT_FREQ = 'annually' ! 'daily', 'monthly', 'annually'
     CHARACTER(LEN=10)  :: vcmax = 'standard' ! "standard" or "Walker2014"
     CHARACTER(LEN=10)  :: POPLUC_RunType = 'static' ! 'static', 'init', 'restart'

     LOGICAL ::                                                               &
          CALL_POP               = .FALSE., & !
          POP_fromZero           = .FALSE., &
          CALL_Climate           = .FALSE., &
          Climate_fromZero       = .FALSE., &
          CASA_fromZero          = .FALSE., &
          POPLUC                 = .FALSE.

     INTEGER  :: &
          CASA_SPIN_STARTYEAR = 1950, &
          CASA_SPIN_ENDYEAR   = 1960, &
          YEARSTART           = 0, &
          YEAREND             = 0, &
          CASA_NREP           = 1
     !--- LN ------------------------------------------]

     CHARACTER(LEN=5) ::                                                      &
          RUN_DIAG_LEVEL  !

     CHARACTER(LEN=3) ::                                                      &
          SSNOW_POTEV,      & !
          DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical)
          LEAF_RESPIRATION    ! either ON or OFF (jhan:Make Logical)

     ! Custom soil respiration - see Ticket #42
     CHARACTER(LEN=10) ::                                                     &
          SMRF_NAME,   & ! Soil Moist Respiration Function
          STRF_NAME      ! Soil Temp Respiration Function

     LOGICAL ::                                                               &
          INITIALIZE_MAPPING    = .FALSE., & !
          CONSISTENCY_CHECK     = .FALSE., & !
          CASA_DUMP_READ        = .FALSE., & !
          CASA_DUMP_WRITE       = .FALSE., & !
          CABLE_RUNTIME_COUPLED = .TRUE. , & !
          LogWorker             = .TRUE. , & ! Write Output of each worker
                                ! L.Stevens - Test Switches
          L_NEW_ROUGHNESS_SOIL  = .FALSE., & !
          L_NEW_RUNOFF_SPEED    = .FALSE., & !
          L_NEW_REDUCE_SOILEVP  = .FALSE., & !

                                ! Switch for customized soil respiration - see Ticket #42
          SRF = .FALSE., &

                                !! vh_js !!
          litter = .FALSE.

     !INH - new switch for revised coupling on implicit step of ACCESS-CM2 Ticket #132
     LOGICAL :: l_revised_coupling = .FALSE.

     !INH -apply revised sensitvity/correction terms to soilsnow energy balance
     LOGICAL :: L_REV_CORR = .FALSE.     !switch to revert to unchanged code

     !MD
     LOGICAL :: GW_MODEL = .FALSE.
     LOGICAL :: alt_forcing = .FALSE.

     !using GSWP3 forcing?
     LOGICAL :: GSWP3 = .FALSE.
     LOGICAL :: or_evap = .FALSE.
     LOGICAL :: test_new_gw=.FALSE.
     LOGICAL :: sync_nc_file=.FALSE.
     INTEGER :: max_spins = -1
     LOGICAL :: fix_access_roots = .FALSE.  !use pft dependent roots in ACCESS
     !ticket#179
     LOGICAL :: soil_thermal_fix=.FALSE.
     !ACCESS roots
     LOGICAL :: access13roots = .FALSE.     !switch to use ACCESS1.3 %froot

  END TYPE kbl_user_switches

  ! instantiate internal switches
  TYPE(kbl_user_switches), SAVE :: cable_user

  ! external files read/written by CABLE
  TYPE filenames_type

     CHARACTER(LEN=500) ::                                                        &
          met,        & ! name of file for CABLE input
          path='./',       & ! path for output and restart files for CABLE and CASA
          out,        & ! name of file for CABLE output
          log,        & ! name of file for execution log
          restart_in = ' ', & ! name of restart file to read
          restart_out,& ! name of restart file to read
          LAI,        & ! name of file for default LAI
          TYPE,       & ! file for default veg/soil type
          veg,        & ! file for vegetation parameters
          soil,       & ! name of file for soil parameters
          soilcolor,  & ! file for soil color(soilcolor_global_1x1.nc)
          inits,      & ! name of file for initialisations
          soilIGBP,   & ! name of file for IGBP soil map
          gw_elev       !name of file for gw/elevation data

  END TYPE filenames_type

  TYPE(filenames_type) :: filename

  ! hydraulic_redistribution switch _soilsnow module
  LOGICAL ::                                                                  &
       redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution

  ! hydraulic_redistribution parameters _soilsnow module
  REAL :: wiltParam=0.5, satuParam=0.8

  TYPE organic_soil_params
     !Below are the soil properties for fully organic soil

     REAL ::    &
          hyds_vec_organic  = 1.0e-4,&
          sucs_vec_organic = 10.3,   &
          clappb_organic = 2.91,     &
          ssat_vec_organic = 0.9,    &
          watr_organic   = 0.1,     &
          sfc_vec_hk      = 1.157407e-06, &
          swilt_vec_hk      = 2.31481481e-8

  END TYPE organic_soil_params

  TYPE gw_parameters_type

     REAL ::                   &
          MaxHorzDrainRate=2e-4,  & !anisintropy * q_max [qsub]
          EfoldHorzDrainRate=2.0, & !e fold rate of q_horz
          MaxSatFraction=2500.0,     & !parameter controll max sat fraction
          hkrz=0.5,               & !hyds_vec variation with z
          zdepth=1.5,             & !level where hyds_vec(z) = hyds_vec(no z)
          frozen_frac=0.05,       & !ice fraction to determine first non-frozen layer for qsub
          SoilEvapAlpha = 1.0,    & !modify field capacity dependence of soil evap limit
          IceAlpha=3.0,           &
          IceBeta=1.0

     REAL :: ice_impedence=5.0

     TYPE(organic_soil_params) :: org

     INTEGER :: level_for_satfrac = 6
     LOGICAL :: ssgw_ice_switch = .FALSE.

     LOGICAL :: subsurface_sat_drainage = .TRUE.

  END TYPE gw_parameters_type

  TYPE(gw_parameters_type), SAVE :: gw_params

  REAL, SAVE ::        &!should be able to change parameters!!!
       max_glacier_snowd=1100.0,&
       snow_ccnsw = 2.0, &
                                !jh!an:clobber - effectively force single layer snow
                                !snmin = 100.0,      & ! for 1-layer;
       snmin = 1.,          & ! for 3-layer;
       max_ssdn = 750.0,    & !
       max_sconds = 2.51,   & !
       frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95)

  !CABLE_LSM: soil/veg params types & subr deleted here
  ! vn10.6-CABLE hacks-hardwires these
  !use these as the basis for namelist vars/files later in offline apps

  !CABLE_LSM: verify these are set if commented here
  !   !---parameters, tolerances, etc. could be set in _directives.h
  !jhan:cable.nml   real, parameter :: RAD_TOLS = 1.0e-2

  !jhan:temporary measure. improve hiding
  !   real, dimension(:,:), pointer,save :: c1, rhoch

  INTERFACE fudge_out
     MODULE PROCEDURE fudge_out_r2D, fudge_out_r1D, fudge_out_r3D, fudge_out_i2D
  END INTERFACE

CONTAINS

  !--- LN ------------------------------------------[
  SUBROUTINE HANDLE_ERR( status, msg )
    ! LN 06/2013
    USE netcdf
    INTEGER :: status
    CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: msg
    IF(status /= NF90_noerr) THEN
       WRITE(*,*)"netCDF error:"
       IF ( PRESENT( msg ) ) WRITE(*,*)msg
       !#define Vanessas_common
       !#ifdef Vanessas_common
       WRITE(*,*) TRIM(NF90_strerror(INT(status,4)))
       !#else
       !       WRITE(*,*) "UM builds with -i8. Therefore call to nf90_strerror is ", &
       !       " invalid. Quick fix to eliminate for now. Build NF90 with -i8, force -i4?"
       !#endif
       STOP -1
    END IF
  END SUBROUTINE HANDLE_ERR

  SUBROUTINE GET_UNIT (IUNIT)

    ! Find an unused unit for intermediate use
    ! PLEASE, use it ONLY when you OPEN AND CLOSE WITHIN THE SAME CALL
    ! or there could be interferences with other files!!!
    ! LN 05/2014

    IMPLICIT NONE

    INTEGER,INTENT(OUT) :: IUNIT
    INTEGER :: i
    LOGICAL :: is_open = .FALSE.

    DO i = 200, 10000
       INQUIRE ( UNIT=i, OPENED=is_open )
       IF ( .NOT. is_open ) EXIT
    END DO
    IUNIT = i

  END SUBROUTINE GET_UNIT

  ELEMENTAL FUNCTION IS_LEAPYEAR( YYYY )
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: YYYY
    LOGICAL :: IS_LEAPYEAR

    IS_LEAPYEAR = .FALSE.
    IF ( ( ( MOD( YYYY,  4 ) .EQ. 0 .AND. MOD( YYYY, 100 ) .NE. 0 ) .OR. &
         MOD( YYYY,400 ) .EQ. 0 ) ) IS_LEAPYEAR = .TRUE.

  END FUNCTION IS_LEAPYEAR

  FUNCTION LEAP_DAY( YYYY )
    IMPLICIT NONE
    INTEGER :: YYYY, LEAP_DAY

    IF ( IS_LEAPYEAR ( YYYY ) ) THEN
       LEAP_DAY = 1
    ELSE
       LEAP_DAY = 0
    END IF
  END FUNCTION LEAP_DAY

  SUBROUTINE YMDHMS2DOYSOD( YYYY,MM,DD,HOUR,MINUTE,SECOND,DOY,SOD )

    ! Compute Day-of-year and second-of-day from given date and time or

    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: YYYY,MM,DD,HOUR,MINUTE,SECOND
    INTEGER,INTENT(OUT) :: DOY,SOD

    !  LOGICAL :: IS_LEAPYEAR
    INTEGER, DIMENSION(12) :: MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    IF ( IS_LEAPYEAR( YYYY ) ) MONTH(2) = 29

    IF ( DD .GT. MONTH(MM) .OR. DD .LT. 1 .OR. &
         MM .GT. 12 .OR. MM .LT. 1 ) THEN
       WRITE(*,*)"Wrong date entered in YMDHMS2DOYSOD "
       WRITE(*,*)"DATE : ",YYYY,MM,DD
       STOP
    ENDIF
    DOY = DD
    IF ( MM .GT. 1 ) DOY = DOY + SUM( MONTH( 1:MM-1 ) )
    SOD = HOUR * 3600 + MINUTE * 60 + SECOND

  END SUBROUTINE YMDHMS2DOYSOD

  SUBROUTINE DOYSOD2YMDHMS( YYYY,DOY,SOD,MM,DD,HOUR,MINUTE,SECOND )

    ! Compute Day-of-year and second-of-day from given date and time or

    IMPLICIT NONE

    INTEGER,INTENT(IN)           :: YYYY,DOY,SOD
    INTEGER,INTENT(OUT)          :: MM,DD
    INTEGER,INTENT(OUT),OPTIONAL :: HOUR,MINUTE,SECOND

    !  LOGICAL :: IS_LEAPYEAR
    INTEGER :: MON, i
    INTEGER, DIMENSION(12) :: MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    IF ( IS_LEAPYEAR( YYYY ) ) MONTH(2) = 29

    IF ( SOD .GE. 86400 .OR. SOD .LT. 0 .OR. &
         DOY .GT. SUM(MONTH) .OR. DOY .LT. 1 ) THEN
       WRITE(*,*)"Wrong date entered in DOYSOD2YMDHMS "
       WRITE(*,*)"YYYY DOY SOD : ",YYYY,DOY,SOD
       STOP
    ENDIF

    MON = 0
    DO i = 1, 12
       IF ( MON + MONTH(i) .LT. DOY ) THEN
          MON = MON + MONTH(i)
       ELSE
          MM  = i
          DD  = DOY - MON
          EXIT
       ENDIF
    END DO
    IF ( PRESENT ( HOUR ) ) HOUR   = INT( REAL(SOD)/3600. )
    IF ( PRESENT (MINUTE) ) MINUTE = INT( ( REAL(SOD) - REAL(HOUR)*3600.) / 60. )
    IF ( PRESENT (SECOND) ) SECOND = SOD - HOUR*3600 - MINUTE*60

  END SUBROUTINE DOYSOD2YMDHMS

  SUBROUTINE LAND2XY( xdimsize, landgrid, x, y )

    ! Convert landgrid to x and y (indices for lat and lon) as
    ! used in CABLE
    ! LN 08/2015

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: xdimsize, landgrid
    INTEGER, INTENT(OUT) :: x, y

    y = INT(REAL((landGrid-1))/REAL(xdimsize)) + 1
    x = landGrid - (y-1) * xdimsize

  END SUBROUTINE LAND2XY

  SUBROUTINE XY2LAND( xdimsize, x, y, landgrid )

    ! Convert x and y (indices for lat and lon) to landgrid
    ! as used in CABLE
    ! LN 08/2015

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: xdimsize, x, y
    INTEGER, INTENT(OUT) :: landgrid

    landgrid = x + ( y - 1 ) * xdimsize

  END SUBROUTINE XY2LAND
  !--- LN ------------------------------------------]

  ! get svn revision number and status
  SUBROUTINE report_version_no( logn )

#ifdef NAG
    USE F90_UNIX_ENV, ONLY: getenv
#endif
    INTEGER, INTENT(IN) :: logn
    ! set from environment variable $HOME
    CHARACTER(LEN=200) ::                                                       &
         myhome,       & ! $HOME (POSIX) environment/shell variable
         fcablerev,    & ! recorded svn revision number at build time
         icable_status   ! recorded svn STATUS at build time (ONLY 200 chars of it)


    INTEGER :: icable_rev, ioerror

    CALL getenv("HOME", myhome)
    fcablerev = TRIM(myhome)//TRIM("/.cable_rev")

    OPEN(440,FILE=TRIM(fcablerev),STATUS='old',ACTION='READ',IOSTAT=ioerror)

    IF(ioerror==0) THEN
       ! get svn revision number (see WRITE comments)
       READ(440,*) icable_rev
    ELSE
       icable_rev=0 !default initialization
       PRINT *, "We'll keep running but the generated revision number "
       PRINT *, " in the log & file will be meaningless."
    ENDIF


    WRITE(logn,*) ''
    WRITE(logn,*) 'Revision nuber: ', icable_rev
    WRITE(logn,*) ''
    WRITE(logn,*)'This is the latest revision of you workin copy as sourced '
    WRITE(logn,*)'by the SVN INFO command at build time. Please note that the'
    WRITE(logn,*)'accuracy of this number is dependent on how recently you '
    WRITE(logn,*)'used SVN UPDATE.'

    ! get svn status (see WRITE comments)
    ! (jhan: make this output prettier & not limitted to 200 chars)
    WRITE(logn,*)'SVN STATUS indicates that you have (at least) the following'
    WRITE(logn,*)'local changes: '
    IF(ioerror==0) THEN
       READ(440,'(A)',IOSTAT=ioerror) icable_status
       WRITE(logn,*) TRIM(icable_status)
       WRITE(logn,*) ''
    ELSE
       WRITE(logn,*) '.cable_rev file does not exist,'
       WRITE(logn,*) 'suggesting you did not build libcable here'
       WRITE(logn,*) ''
    ENDIF

    CLOSE(440)

  END SUBROUTINE report_version_no

  SUBROUTINE init_veg_from_vegin(ifmp,fmp, veg, soil_zse )
    USE cable_def_types_mod, ONLY : veg_parameter_type, ms
    INTEGER ::  ifmp,  & ! start local mp, # landpoints (jhan:when is this not 1 )
         fmp     ! local mp, # landpoints
    REAL, DIMENSION(ms) :: soil_zse

    TYPE(veg_parameter_type) :: veg

    INTEGER :: is
    REAL :: totdepth
    INTEGER :: h

    ! Prescribe parameters for current gridcell based on veg/soil type (which
    ! may have loaded from default value file or met file):
    DO h = ifmp, fmp          ! over each patch in current grid
       veg%frac4(h)    = vegin%frac4(veg%iveg(h))
       veg%taul(h,1)    = vegin%taul(1,veg%iveg(h))
       veg%taul(h,2)    = vegin%taul(2,veg%iveg(h))
       veg%refl(h,1)    = vegin%refl(1,veg%iveg(h))
       veg%refl(h,2)    = vegin%refl(2,veg%iveg(h))
       veg%canst1(h)   = vegin%canst1(veg%iveg(h))
       veg%dleaf(h)    = vegin%dleaf(veg%iveg(h))
       veg%vcmax(h)    = vegin%vcmax(veg%iveg(h))
       veg%ejmax(h)    = vegin%ejmax(veg%iveg(h))
       veg%hc(h)       = vegin%hc(veg%iveg(h))
       veg%xfang(h)    = vegin%xfang(veg%iveg(h))
       veg%vbeta(h)    = vegin%vbeta(veg%iveg(h))
       veg%xalbnir(h)  = vegin%xalbnir(veg%iveg(h))
       veg%rp20(h)     = vegin%rp20(veg%iveg(h))
       veg%rpcoef(h)   = vegin%rpcoef(veg%iveg(h))
       veg%rs20(h)     = vegin%rs20(veg%iveg(h))
       veg%shelrb(h)   = vegin%shelrb(veg%iveg(h))
       veg%wai(h)      = vegin%wai(veg%iveg(h))
       veg%a1gs(h)     = vegin%a1gs(veg%iveg(h))
       veg%d0gs(h)     = vegin%d0gs(veg%iveg(h))
       veg%vegcf(h)    = vegin%vegcf(veg%iveg(h))
       veg%extkn(h)    = vegin%extkn(veg%iveg(h))
       veg%tminvj(h)   = vegin%tminvj(veg%iveg(h))
       veg%tmaxvj(h)   = vegin%tmaxvj(veg%iveg(h))
       veg%g0(h)       = vegin%g0(veg%iveg(h)) ! Ticket #56
       veg%g1(h)       = vegin%g1(veg%iveg(h)) ! Ticket #56
       veg%a1gs(h)   = vegin%a1gs(veg%iveg(h))
       veg%d0gs(h)   = vegin%d0gs(veg%iveg(h))
       veg%alpha(h)  = vegin%alpha(veg%iveg(h))
       veg%convex(h) = vegin%convex(veg%iveg(h))
       veg%cfrd(h)   = vegin%cfrd(veg%iveg(h))
       veg%gswmin(h) = vegin%gswmin(veg%iveg(h))
       veg%conkc0(h) = vegin%conkc0(veg%iveg(h))
       veg%conko0(h) = vegin%conko0(veg%iveg(h))
       veg%ekc(h)    = vegin%ekc(veg%iveg(h))
       veg%eko(h)    = vegin%eko(veg%iveg(h))
       veg%rootbeta(h)  = vegin%rootbeta(veg%iveg(h))
       veg%zr(h)       = vegin%zr(veg%iveg(h))
       veg%clitt(h)    = vegin%clitt(veg%iveg(h))

       ! mgk576, hydraulics stuff
       veg%sf(h)       = vegin%sf(veg%iveg(h))       ! mgk576
       veg%psi_f(h)    = vegin%psi_f(veg%iveg(h))    ! mgk576
       veg%X_hyd(h)    = vegin%X_hyd(veg%iveg(h))    ! mgk576
       veg%p50(h)      = vegin%p50(veg%iveg(h))      ! mgk576
       veg%s50(h)      = vegin%s50(veg%iveg(h))      ! mgk576
       veg%kp_sat(h)   = vegin%kp_sat(veg%iveg(h))      ! mgk576
       veg%Cl(h)       = vegin%Cl(veg%iveg(h))      ! mgk576
       veg%Cs(h)       = vegin%Cs(veg%iveg(h))      ! mgk576

    END DO ! over each veg patch in land point

    ! calculate vegin%froot from using rootbeta and soil depth
    ! (Jackson et al. 1996, Oceologica, 108:389-411)
    totdepth = 0.0
    DO is = 1, ms-1
       totdepth = totdepth + soil_zse(is) * 100.0  ! unit in centimetres
       print*, "wtf", is, veg%rootbeta, totdepth, soil_zse(is), SUM(soil_zse)
       veg%froot(:, is) = MIN( 1.0, 1.0-veg%rootbeta(:)**totdepth )
    END DO
    veg%froot(:, ms) = 1.0 - veg%froot(:, ms-1)
    DO is = ms-1, 2, -1
       veg%froot(:, is) = veg%froot(:, is)-veg%froot(:,is-1)
    END DO





  END SUBROUTINE init_veg_from_vegin


  FUNCTION IS_CASA_TIME(iotype, yyyy, ktau, kstart, koffset, kend, ktauday, logn)

    ! Correctly determine if it is time to dump-read or standard-write
    ! casa output from cable_driver.
    ! Writing casa-dump data is handled in casa_cable and therefore not \
    ! captured here
    !cable_common module was intended to be unequivocally common to all
    !applications. iovars is an offline module and so not appropriate to include
    !here. Suggested FIX is to move decs of vars needed (e.g. leaps) to here, and
    !then use common in iovars
#ifdef Vanessas_common
    USE cable_IO_vars_module, ONLY: leaps
#endif
    IMPLICIT NONE

    LOGICAL   :: IS_CASA_TIME
    INTEGER  ,INTENT(IN) :: yyyy, ktau, kstart, koffset, kend, ktauday, logn
    CHARACTER,INTENT(IN) :: iotype*5
    LOGICAL   :: is_eod, is_eom, is_eoy
    INTEGER   :: doy, m
    INTEGER, DIMENSION(12) :: MONTH

    is_eom       = .FALSE.
    is_eoy       = .FALSE.
    IS_CASA_TIME = .FALSE.

    MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    is_eod = ( MOD((ktau-kstart+1+koffset),ktauday).EQ.0 )
    IF ( .NOT. is_eod ) RETURN    ! NO if it is not end of day

#ifdef Vanessas_common
    IF ( IS_LEAPYEAR( YYYY ) .AND. leaps ) THEN
       MONTH(2) = 29
    ELSE
       MONTH(2) = 28
    ENDIF
#endif

    ! Check for reading from dump-file (hard-wired to daily casa-timestep)
    IF ( iotype .EQ. "dread" ) THEN
       IF ( CABLE_USER%CASA_DUMP_READ )  IS_CASA_TIME = .TRUE.
       ! Check for writing of casa dump output
    ELSE IF ( iotype .EQ. "dwrit" ) THEN
       IF ( CABLE_USER%CASA_DUMP_WRITE ) IS_CASA_TIME = .TRUE.
       ! Check for writing of casa standard output
    ELSE IF ( iotype .EQ. "write" ) THEN

       doy = NINT(REAL(ktau-kstart+1+koffset)/REAL(ktauday))
       DO m = 1, 12
          IF ( doy .EQ. SUM(MONTH(1:m)) ) THEN
             is_eom = .TRUE.
             IF ( m .EQ. 12 ) is_eoy = .TRUE.
             EXIT
          ENDIF
       END DO

       SELECT CASE ( TRIM(CABLE_USER%CASA_OUT_FREQ) )
       CASE ("daily"   ) ; IS_CASA_TIME = .TRUE.
       CASE ("monthly" ) ; IF ( is_eom ) IS_CASA_TIME = .TRUE.
       CASE ("annually") ; IF ( is_eoy ) IS_CASA_TIME = .TRUE.
       END SELECT
    ELSE
       WRITE(logn,*)"Wrong statement 'iotype'", iotype, "in call to IS_CASA_TIME"
       WRITE(*   ,*)"Wrong statement 'iotype'", iotype, "in call to IS_CASA_TIME"
       STOP -1
    ENDIF

  END FUNCTION IS_CASA_TIME


  SUBROUTINE fudge_out_i2D( i,j, var, varname, vzero, vval )
    ! interfaces on these
    INTEGER :: i,j
    INTEGER, DIMENSION(:,:) :: var
    ! ft changes with interface
    CHARACTER(len=*), PARAMETER :: &
         ft = '(  "fudge: ", A10, "(", I2.1, ",", I2.1, X, ") = ", I1.1 )'

    CHARACTER(len=*) :: varname
    LOGICAL :: vzero
    INTEGER :: vval

    ! content changes with interface
    var = var(i,j)
    IF( (vzero) ) var = vval
    WRITE (6, ft) varname,i, var(i,j)
  END SUBROUTINE fudge_out_i2D


  SUBROUTINE fudge_out_r1D( i, var, varname, vzero, vval )
    ! interfaces on these
    INTEGER :: i
    REAL, DIMENSION(:) :: var
    ! ft changes with interface
    CHARACTER(len=*), PARAMETER :: &
         ft = '(  "fudge: ", A10, "(", I2.1, X, ") = ", F15.3 )'

    CHARACTER(len=*) :: varname
    LOGICAL :: vzero
    REAL :: vval

    ! content changes with interface
    var = var(i)
    IF( (vzero) ) var = vval
    WRITE (6, ft) varname,i, var(i)
  END SUBROUTINE fudge_out_r1D

  SUBROUTINE fudge_out_r2D( i,j, var, varname, vzero, vval )
    ! interfaces on these
    INTEGER :: i,j
    REAL, DIMENSION(:,:) :: var
    ! ft changes with interface
    CHARACTER(len=*), PARAMETER :: &
         ft = '(  "fudge: ", A10, "(", I2.1, ",", I2.1, X, ") = ", F15.3 )'

    CHARACTER(len=*) :: varname
    LOGICAL :: vzero
    REAL :: vval

    ! content changes with interface
    var = var(i,j)
    IF( (vzero) ) var = vval
    WRITE (6, ft) varname,i,j, var(i,j)
  END SUBROUTINE fudge_out_r2D

  SUBROUTINE fudge_out_r3D( i,j,k, var, varname, vzero, vval )
    ! interfaces on these
    INTEGER :: i,j,k
    REAL, DIMENSION(:,:,:) :: var
    ! ft changes with interface
    CHARACTER(len=*), PARAMETER :: &
         ft = '(  "fudge: ", A10, "(",  I2.1, ",",I2.1, ",", I2.1, X, ") = ", F15.3 )'

    CHARACTER(len=*) :: varname
    LOGICAL :: vzero
    REAL :: vval

    ! content changes with interface
    var = var(i,j,k)
    IF( (vzero) ) var = vval
    WRITE (6, ft) varname,i,j,k, var(i,j,k)
  END SUBROUTINE fudge_out_r3D



END MODULE cable_common_module
