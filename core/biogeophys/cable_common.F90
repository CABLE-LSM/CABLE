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
  
  IMPLICIT NONE

  !---allows reference to "gl"obal timestep in run (from atm_step)
  !---total number of timesteps, and processing node
  INTEGER, SAVE :: ktau_gl, kend_gl, knode_gl, kwidth_gl
  INTEGER, SAVE :: CurYear  ! current year of multiannual run

  ! user switches turned on/off by the user thru namelists

  ! trunk modifications protected by these switches
  TYPE hide_switches
     LOGICAL ::                                                               &
          ! L.Stevens - Test Switches
          L_NEW_ROUGHNESS_SOIL  = .FALSE., & ! from Ticket?
          L_NEW_RUNOFF_SPEED    = .FALSE., & ! from Ticket?
          L_NEW_REDUCE_SOILEVP  = .FALSE., & ! from Ticket?
          Ticket46 = .FALSE.,              & !
          !jhan: default should be FALSE, bu set up nml etc
          Ticket49Bug1 = .false.,           & !
          Ticket49Bug2 = .false.,           & !
          Ticket49Bug3 = .false.,           & !
          Ticket49Bug4 = .false.,           & !
          Ticket49Bug5 = .false.,           & !
          Ticket49Bug6 = .false.              !

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
          um_radiation = .FALSE.
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
          PHENOLOGY_SWITCH = 'climate'   ! alternative is 'climate'
    !--- LN ------------------------------------------[

     ! Ticket #56
     CHARACTER(LEN=20) ::                                                     &
    ! GS_SWITCH='leuning'
     GS_SWITCH='medlyn'

     CHARACTER(LEN=20) :: g0_switch = 'default'  ! 'maximum',takes max of g0 and A*X, 'standard' is additive
     CHARACTER(LEN=10) :: RunIden       = 'STANDARD'  !
     CHARACTER(LEN=4)  :: MetType       = ' ' !
     CHARACTER(LEN=20) :: SOIL_STRUC    = "default" ! 'default' or 'sli'
     CHARACTER(LEN=3)  :: POP_out       = 'rst' ! POP output type ('epi' or 'rst' or 'ini')
     CHARACTER(LEN=50) :: POP_rst       = '' !
     CHARACTER(LEN=200) :: POP_restart_in = ''
     CHARACTER(LEN=200) :: POP_restart_out = ''
     CHARACTER(LEN=200) :: POP_outfile       = '' !
     CHARACTER(LEN=200) :: climate_restart_in = ''
     CHARACTER(LEN=200) :: climate_restart_out = ''
     CHARACTER(LEN=200) :: LUC_outfile       = '' !
     CHARACTER(LEN=200) :: LUC_restart_in = ''
     CHARACTER(LEN=200) :: LUC_restart_out = ''
     CHARACTER(LEN=8)  :: CASA_OUT_FREQ = 'annually' ! 'daily', 'monthly', 'annually'
     CHARACTER(LEN=10)  :: vcmax = 'standard' ! "standard" or "Walker2014"
     CHARACTER(LEN=10)  :: POPLUC_RunType = 'static' ! 'static', 'init', 'restart'
     CHARACTER(LEN=200) :: BLAZE_outfile       = '' !
     LOGICAL ::                                                               &
          CALL_POP               = .FALSE., & !
          POP_fromZero           = .FALSE., &
          CALL_Climate           = .FALSE., &
          CALL_BLAZE             = .FALSE., &
          Climate_fromZero       = .FALSE., &
          CASA_fromZero          = .FALSE., &
          POPLUC                 = .FALSE., &
          explicit_gm            = .FALSE., &     ! explicit (finite) mesophyll conductance
          acclim_autoresp            = .TRUE., &
          coordinate_photosyn    = .TRUE., &
          acclimate_photosyn      = .FALSE., &
          acclimate_autoresp_seasonal = .FALSE., &  ! acclimates to last 30 d, otherwise annual.
          limit_labile           = .FALSE., &
          Cumberland_soil        = .FALSE.   ! sets special CP soil params in calbe_sli_utils.F90

     INTEGER  ::  &
          CASA_SPIN_STARTYEAR = 1950, &
          CASA_SPIN_ENDYEAR   = 1960, &
          YEARSTART           = 0, &
          YEAREND             = 0, &
          CASA_NREP           = 1
     CHARACTER(len=7) :: &
          BURNT_AREA          = "SIMFIRE" ! either SIMFIRE or GFED31
     CHARACTER(len=8) :: &
          BLAZE_TSTEP         = "DAILY"   ! either DAILY, MONTHLY, ANNUALLY
     CHARACTER(len=6) :: &
          SIMFIRE_REGION      = "GLOBAL"  ! either GLOBAL, EUROPE, ANZ

    !--- LN ------------------------------------------]

     CHARACTER(LEN=5) ::                                                      &
          RUN_DIAG_LEVEL  !

     CHARACTER(LEN=3) ::                                                      &
          SSNOW_POTEV,      & !
          DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical)
          LEAF_RESPIRATION    ! either ON or OFF (jhan:Make Logical)

     ! Custom soil respiration - see Ticket #42
     CHARACTER(LEN=15) ::                                                     &
          !SMRF_NAME = 'Trudinger2016',   & ! Soil Moist Respiration Function
          !STRF_NAME = 'CASA-CNP'     ! Soil Temp Respiration Function
          !STRF_NAME = 'LT1994'    ! Soil Temp Respiration Function
          STRF_NAME = 'CASA-CNP',  &     ! DAMM Reverse M-M Enzyme Kinetics (Sihi et al, AFM 2018)
          SMRF_NAME = 'CASA-CNP'       ! ditto
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
          SRF = .TRUE., &

          !! vh_js !!
          litter = .FALSE.

     logical ::            c13o2 = .false.              ! switch 13CO2 calculations on
     logical ::            c13o2_simple_disc = .false.  ! simple or full leaf discrimination
     character(len=200) :: c13o2_delta_atm_file = ''    ! atmospheric delta 13CO2 values
     character(len=200) :: c13o2_outfile = ''           ! 13C Casa and LUC output file
     character(len=200) :: c13o2_restart_in_flux = ''   ! 13CO2 restart Canopy input file
     character(len=200) :: c13o2_restart_out_flux = ''  ! 13CO2 restart Canopy output file
     character(len=200) :: c13o2_restart_in_pools = ''  ! 13CO2 restart Casa input file
     character(len=200) :: c13o2_restart_out_pools = '' ! 13CO2 restart Casa output file
     character(len=200) :: c13o2_restart_in_luc = ''    ! 13CO2 restart LUC input file
     character(len=200) :: c13o2_restart_out_luc = ''   ! 13CO2 restart LUC output file

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
          type,       & ! file for default veg/soil type
          veg,        & ! file for vegetation parameters
          soil,       & ! name of file for soil parameters
          soilcolor,  & ! file for soil color(soilcolor_global_1x1.nc)
          inits,      & ! name of file for initialisations
          soilIGBP      ! name of file for IGBP soil map

  END TYPE filenames_type

  TYPE(filenames_type) :: filename

  ! hydraulic_redistribution switch _soilsnow module
  LOGICAL ::                                                                  &
       redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution

  ! hydraulic_redistribution parameters _soilsnow module
  REAL :: wiltParam=0.5, satuParam=0.8


  ! soil parameters read from file(filename%soil def. in cable.nml)
  ! & veg parameters read from file(filename%veg def. in cable.nml)
  TYPE soilin_type

     REAL, DIMENSION(:),ALLOCATABLE ::                                        &
          silt,    & !
          clay,    & !
          sand,    & !
          swilt,   & !
          sfc,     & !
          ssat,    & !
          bch,     & !
          hyds,    & !
          sucs,    & !
          rhosoil, & !
          css,     & !
          c3         !

  END TYPE soilin_type


  TYPE vegin_type

     REAL, DIMENSION(:),ALLOCATABLE ::                                        &
          canst1,     & !
          dleaf,      & !
          length,     & !
          width,      & !
          vcmax,      & !
          ejmax,      & !
          vcmaxcc,    & ! Cc-based Vcmax25
          ejmaxcc,    & ! Cc-based Jmax25
          gmmax,      & ! gmmax25
          hc,         & !
          xfang,      & !
          rp20,       & !
          rpcoef,     & !
          rs20,       & !
          wai,        & !
          rootbeta,   & !
          shelrb,     & !
          vegcf,      & !
          frac4,      & !
          xalbnir,    & !
          extkn,      & !
          tminvj,     & !
          tmaxvj,     & !
          vbeta,      &
          a1gs,       &
          d0gs,       &
          alpha,      &
          convex,     &
          cfrd,       &
          gswmin,     &
          conkc0,     &
          conko0,     &
          ekc,        &
          eko,        &
          g0,         & !  Ticket #56
          g1,         & !  Ticket #56
          zr,         &
          clitt,      &
          gamma

     REAL, DIMENSION(:,:),ALLOCATABLE ::                                      &
          froot,      & !
          cplant,     & !
          csoil,      & !
          ratecp,     & !
          ratecs,     & !
          refl,     & !
          taul        !

  END TYPE vegin_type

  CHARACTER(LEN=70), DIMENSION(:), POINTER ::                                 &
       veg_desc => null(),   & ! decriptions of veg type
       soil_desc => null()     ! decriptns of soil type

  TYPE(soilin_type), SAVE  :: soilin
  TYPE(vegin_type),  SAVE  :: vegin

  !   !---parameters, tolerances, etc. could be set in _directives.h
  !jhan:cable.nml   real, parameter :: RAD_TOLS = 1.0e-2

  !jhan:temporary measure. improve hiding
  !   real, dimension(:,:), pointer,save :: c1, rhoch

CONTAINS


  SUBROUTINE get_type_parameters(logn,vegparmnew, classification)

    ! Gets parameter values for each vegetation type and soil type.

    USE cable_def_types_mod, ONLY : mvtype, ms, ncs, ncp, mstype, nrb

    INTEGER,INTENT(IN) :: logn     ! log file unit number

    CHARACTER(LEN=4), INTENT(INOUT), OPTIONAL :: classification

    LOGICAL,INTENT(IN)      :: vegparmnew ! new format input file

    CHARACTER(LEN=80) :: comments
    CHARACTER(LEN=10) :: vegtypetmp
    CHARACTER(LEN=25) :: vegnametmp

    REAL    :: notused
    INTEGER :: ioerror ! input error integer
    INTEGER :: a, jveg ! do loop counter



    !================= Read in vegetation type specifications: ============
    OPEN(40,FILE=filename%veg,STATUS='old',ACTION='READ',IOSTAT=ioerror)

    IF(ioerror/=0) then
       STOP 'CABLE_log: Cannot open veg type definitions.'
    ENDIF

    IF (vegparmnew) THEN

       ! assume using IGBP/CSIRO vegetation types
       READ(40,*) comments
       READ(40,*) mvtype
       IF( present(classification) )                                         &
            WRITE(classification,'(a4)') comments(1:4)

    ELSE

       ! assume using CASA vegetation types
       !classification = 'CASA'
       READ(40,*)
       READ(40,*)
       READ(40,*) mvtype ! read # vegetation types
       READ(40,*)
       READ(40,*)
       comments = 'CASA'

    END IF

    WRITE(logn, '(A31,I3,1X,A10)') '  Number of vegetation types = ',        &
         mvtype,TRIM(comments)


    ! Allocate memory for type-specific vegetation parameters:
    ALLOCATE (                                                               &
         vegin%canst1( mvtype ), vegin%dleaf( mvtype ),                        &
         vegin%length( mvtype ), vegin%width( mvtype ),                        &
         vegin%vcmax( mvtype ),  vegin%ejmax( mvtype ),                        &
         vegin%vcmaxcc(mvtype), vegin%ejmaxcc(mvtype),                         &
         vegin%gmmax(mvtype),                                                  &
         vegin%hc( mvtype ), vegin%xfang( mvtype ),                            &
         vegin%rp20( mvtype ), vegin%rpcoef( mvtype ),                         &
         vegin%rs20( mvtype ), vegin%wai( mvtype ),                            &
         vegin%rootbeta( mvtype ), vegin%shelrb( mvtype ),                     &
         vegin%vegcf( mvtype ), vegin%frac4( mvtype ),                         &
         vegin%xalbnir( mvtype ), vegin%extkn( mvtype ),                       &
         vegin%tminvj( mvtype ), vegin%tmaxvj( mvtype ),                       &
         vegin%vbeta( mvtype ), vegin%froot( ms, mvtype ),                     &
         vegin%cplant( ncp, mvtype ), vegin%csoil( ncs, mvtype ),              &
         vegin%ratecp( ncp, mvtype ), vegin%ratecs( ncs, mvtype ),             &
         vegin%refl( nrb, mvtype ), vegin%taul( nrb, mvtype ),                 &
         veg_desc( mvtype ),                                                   &
         vegin%a1gs(mvtype), vegin%d0gs(mvtype),                               &
         vegin%alpha(mvtype),vegin%convex(mvtype),vegin%cfrd(mvtype),          &
         vegin%gswmin(mvtype),vegin%conkc0(mvtype), vegin%conko0(mvtype),      &
         vegin%ekc(mvtype), vegin%eko(mvtype),                                 &
         ! Ticket #56
         vegin%g0( mvtype ), vegin%g1( mvtype ),                               &
         !! vh_veg_params !!
         vegin%zr(mvtype), vegin%clitt(mvtype), vegin%gamma(mvtype))


    IF( vegparmnew ) THEN    ! added to read new format (BP dec 2007)

       ! Read in parameter values for each vegetation type:
       DO a = 1,mvtype

          READ(40,*) jveg, vegtypetmp, vegnametmp

          IF( jveg .GT. mvtype )                                             &
               STOP 'jveg out of range in parameter file'

          veg_desc(jveg) = vegnametmp

          READ(40,*) vegin%hc(jveg), vegin%xfang(jveg), vegin%width(jveg),   &
               vegin%length(jveg), vegin%frac4(jveg)
          ! only refl(1:2) and taul(1:2) used
          READ(40,*) vegin%refl(1:3,jveg) ! rhowood not used ! BP may2011
          READ(40,*) vegin%taul(1:3,jveg) ! tauwood not used ! BP may2011
          READ(40,*) notused, notused, notused, vegin%xalbnir(jveg)
          READ(40,*) notused, vegin%wai(jveg), vegin%canst1(jveg),           &
               vegin%shelrb(jveg), vegin%vegcf(jveg), vegin%extkn(jveg)

          READ(40,*) vegin%vcmax(jveg), vegin%rp20(jveg),                    &
               vegin%rpcoef(jveg),                                     &
               vegin%rs20(jveg)
          READ(40,*) vegin%tminvj(jveg), vegin%tmaxvj(jveg),                 &
               vegin%vbeta(jveg), vegin%rootbeta(jveg),                      &
               vegin%zr(jveg), vegin%clitt(jveg)
          READ(40,*) vegin%cplant(1:3,jveg), vegin%csoil(1:2,jveg)
          ! rates not currently set to vary with veg type
          READ(40,*) vegin%ratecp(1:3,jveg), vegin%ratecs(1:2,jveg)
          READ(40,*) vegin%a1gs(jveg), vegin%d0gs(jveg), vegin%alpha(jveg), vegin%convex(jveg), vegin%cfrd(jveg)
          READ(40,*) vegin%gswmin(jveg), vegin%conkc0(jveg), vegin%conko0(jveg), vegin%ekc(jveg), vegin%eko(jveg)
          READ(40,*) vegin%g0(jveg), vegin%g1(jveg)      ! Ticket #56
          READ(40,*) vegin%gamma(jveg), vegin%gmmax(jveg)
       END DO

    ELSE

       DO a = 1,mvtype
          READ(40,'(8X,A70)') veg_desc(a) ! Read description of each veg type
       END DO

       READ(40,*); READ(40,*)
       READ(40,*) vegin%canst1
       READ(40,*) vegin%width
       READ(40,*) vegin%length
       READ(40,*) vegin%vcmax
       READ(40,*) vegin%hc
       READ(40,*) vegin%xfang
       READ(40,*) vegin%rp20
       READ(40,*) vegin%rpcoef
       READ(40,*) vegin%rs20
       READ(40,*) vegin%shelrb
       READ(40,*) vegin%frac4
       READ(40,*) vegin%wai
       READ(40,*) vegin%vegcf
       READ(40,*) vegin%extkn
       READ(40,*) vegin%tminvj
       READ(40,*) vegin%tmaxvj
       READ(40,*) vegin%vbeta
       READ(40,*) vegin%xalbnir
       READ(40,*) vegin%rootbeta
       READ(40,*) vegin%cplant(1,:)
       READ(40,*) vegin%cplant(2,:)
       READ(40,*) vegin%cplant(3,:)
       READ(40,*) vegin%csoil(1,:)
       READ(40,*) vegin%csoil(2,:)
       READ(40,*)
       READ(40,*) vegin%ratecp(:,1)
       READ(40,*) vegin%a1gs
       READ(40,*) vegin%d0gs
       READ(40,*) vegin%alpha
       READ(40,*) vegin%convex
       READ(40,*) vegin%cfrd
       READ(40,*) vegin%gswmin
       READ(40,*) vegin%conkc0
       READ(40,*) vegin%conko0
       READ(40,*) vegin%ekc
       READ(40,*) vegin%eko

       ! Set ratecp to be the same for all veg types:
       vegin%ratecp(1,:)=vegin%ratecp(1,1)
       vegin%ratecp(2,:)=vegin%ratecp(2,1)
       vegin%ratecp(3,:)=vegin%ratecp(3,1)
       READ(40,*)
       READ(40,*) vegin%ratecs(:,1)
       vegin%ratecs(1,:)=vegin%ratecs(1,1)
       vegin%ratecs(2,:)=vegin%ratecs(2,1)

       ! old table does not have taul and refl ! BP may2011
       vegin%taul(1,:) = 0.07
       vegin%taul(2,:) = 0.425
       vegin%taul(3,:) = 0.0
       vegin%refl(1,:) = 0.07
       vegin%refl(2,:) = 0.425
       vegin%refl(3,:) = 0.0

       READ(40,*) vegin%g0 ! Ticket #56
       READ(40,*) vegin%g1 ! Ticket #56

    ENDIF

    WRITE(6,*)'CABLE_log:Closing veg params file: ',trim(filename%veg)

    CLOSE(40)

    ! new calculation dleaf since April 2012 (cable v1.8 did not use width)
    vegin%dleaf = SQRT(vegin%width * vegin%length)


    !================= Read in soil type specifications: ============
    OPEN(40,FILE=filename%soil,STATUS='old',ACTION='READ',IOSTAT=ioerror)

    IF(ioerror/=0) then
       STOP 'CABLE_log: Cannot open soil type definitions.'
    ENDIF

    READ(40,*); READ(40,*)
    READ(40,*) mstype ! Number of soil types
    READ(40,*); READ(40,*)

    ALLOCATE ( soil_desc(mstype) )
    ALLOCATE ( soilin%silt(mstype), soilin%clay(mstype), soilin%sand(mstype) )
    ALLOCATE ( soilin%swilt(mstype), soilin%sfc(mstype), soilin%ssat(mstype) )
    ALLOCATE ( soilin%bch(mstype), soilin%hyds(mstype), soilin%sucs(mstype) )
    ALLOCATE ( soilin%rhosoil(mstype), soilin%css(mstype) )

    DO a = 1,mstype
       READ(40,'(8X,A70)') soil_desc(a) ! Read description of each soil type
    END DO

    READ(40,*); READ(40,*)
    READ(40,*) soilin%silt
    READ(40,*) soilin%clay
    READ(40,*) soilin%sand
    READ(40,*) soilin%swilt
    READ(40,*) soilin%sfc
    READ(40,*) soilin%ssat
    READ(40,*) soilin%bch
    READ(40,*) soilin%hyds
    READ(40,*) soilin%sucs
    READ(40,*) soilin%rhosoil
    READ(40,*) soilin%css

    CLOSE(40)

  END SUBROUTINE get_type_parameters

    !--- LN ------------------------------------------[
  subroutine handle_err(status, msg) ! LN 06/2013

    use netcdf

    integer,          intent(in)           :: status
    character(len=*), intent(in), optional :: msg

    if (status /= NF90_noerr) then
       write(*,*) "netCDF error:"
       if (present(msg)) write(*,*) msg
!#define Vanessas_common
!#ifdef Vanessas_common
       write(*,*) trim(NF90_strerror(status))
!#else
!       WRITE(*,*) "UM builds with -i8. Therefore call to nf90_strerror is ", &
!       " invalid. Quick fix to eliminate for now. Build NF90 with -i8, force -i4?"
!#endif
       stop -1
    end if

  end subroutine handle_err


  subroutine get_unit(iunit)

    ! Find an unused unit for intermediate use
    ! PLEASE, use it ONLY when you OPEN AND CLOSE WITHIN THE SAME CALL
    ! or there could be interferences with other files!!!
    ! LN 05/2014

    implicit none

    integer, intent(out) :: iunit

    integer :: i
    logical :: is_open = .false.

    do i=200, 10000
       inquire(unit=i, opened=is_open)
       if (.not. is_open) exit
    end do
    iunit = i

  end subroutine get_unit


  elemental pure function is_leapyear(yyyy)

    implicit none

    integer, intent(in) :: yyyy
    logical             :: is_leapyear

    is_leapyear = .false.
    if ( (mod(yyyy,4).eq.0 .and. mod(yyyy,100).ne.0) .or. &
         mod(yyyy,400).eq.0 ) is_leapyear = .true.

  end function is_leapyear


  elemental pure function leap_day(yyyy)

    implicit none

    integer, intent(in) :: yyyy
    integer             :: leap_day

    if (is_leapyear(yyyy)) then
       leap_day = 1
    else
       leap_day = 0
    end if

  end function leap_day


  subroutine ymdhms2doysod(yyyy, mm, dd, hour, minute, second, doy, sod)

    ! Compute Day-of-year and second-of-day from given date and time or

    implicit none

    integer, intent(in)  :: yyyy, mm, dd, hour, minute, second
    integer, intent(out) :: doy, sod

    !  logical :: is_leapyear
    integer, dimension(12) :: month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    if (is_leapyear(yyyy)) month(2) = 29

    if (dd.gt.month(mm) .or. dd.lt.1 .or. &
         mm.gt.12 .or. mm.lt.1) then
       write(*,*) "Wrong date entered in subroutine ymdhms2doysod"
       write(*,*) "Date (ymdhms): ", yyyy, mm, dd, hour, minute, second
       stop 3
    endif

    doy = dd
    if (mm .gt. 1) doy = doy + sum(month(1:mm-1))

    sod = hour * 3600 + minute * 60 + second

  end subroutine ymdhms2doysod


  subroutine doysod2ymdhms(yyyy, doy, sod, mm, dd, hour, minute, second)

    ! Compute Day-of-year and second-of-day from given date and time or

    implicit none

    integer, intent(in)            :: yyyy, doy, sod
    integer, intent(out)           :: mm, dd
    integer, intent(out), optional :: hour, minute, second

    !  logical :: is_leapyear
    integer :: mon, i
    integer :: ihour, iminute, isecond
    integer, dimension(12) :: month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    if (is_leapyear(yyyy)) month(2) = 29

    if (sod.ge.86400 .or. sod.lt.0 .or. &
         doy.gt.sum(month) .or. doy.lt.1) then
       write(*,*) "Wrong date entered in subroutine doysod2ymdhms"
       write(*,*) "yyyy doy sod: ",yyyy, doy, sod
       stop 3
    endif

    mon = 0
    do i=1, 12
       if (mon+month(i) .lt. doy) then
          mon = mon + month(i)
       else
          mm = i
          dd = doy - mon
          exit
       endif
    end do

    ihour   = int(real(sod)/3600.)
    iminute = int((real(sod) - real(ihour)*3600.) / 60.)
    isecond = sod - ihour*3600 - iminute*60

    if (present(hour))   hour   = ihour
    if (present(minute)) minute = iminute
    if (present(second)) second = isecond

  end subroutine doysod2ymdhms


  subroutine land2xy(xdimsize, landgrid, x, y)

    ! Convert landgrid to x and y (indices for lat and lon) as
    ! used in CABLE
    ! LN 08/2015

    implicit none

    integer, intent(in)  :: xdimsize, landgrid
    integer, intent(out) :: x, y

    y = int(real(landgrid-1)/real(xdimsize)) + 1
    x = landgrid - (y-1) * xdimsize

  end subroutine land2xy


  subroutine xy2land(xdimsize, x, y, landgrid)

    ! Convert x and y (indices for lat and lon) to landgrid
    ! as used in CABLE
    ! LN 08/2015

    implicit none

    integer, intent(in)  :: xdimsize, x, y
    integer, intent(out) :: landgrid

    landgrid = x + (y-1) * xdimsize

  end subroutine xy2land
    !--- LN ------------------------------------------]


  ! get svn revision number and status
  SUBROUTINE report_version_no( logn )

#ifdef __NAG__
    USE F90_UNIX_ENV, only: getenv
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

    OPEN(440, FILE=TRIM(fcablerev), STATUS='old', ACTION='READ', IOSTAT=ioerror)

    IF(ioerror==0) then
       ! get svn revision number (see WRITE comments)
       READ(440,*) icable_rev
    ELSE
       icable_rev = 0 !default initialization
       write(*,'(a)') "We'll keep running but the generated revision number in the log & file will be meaningless."
    ENDIF


    WRITE(logn,*) ''
    WRITE(logn,*) 'Revision nuber: ', icable_rev
    WRITE(logn,*) ''
    WRITE(logn,*) 'This is the latest revision of you workin copy as sourced '
    WRITE(logn,*) 'by the SVN INFO command at build time. Please note that the'
    WRITE(logn,*) 'accuracy of this number is dependent on how recently you '
    WRITE(logn,*) 'used SVN UPDATE.'

    ! get svn status (see WRITE comments)
    ! (jhan: make this output prettier & not limitted to 200 chars)
    WRITE(logn,*) 'SVN STATUS indicates that you have (at least) the following'
    WRITE(logn,*) 'local changes: '
    IF(ioerror==0) then
       READ(440,'(A)',IOSTAT=ioerror) icable_status
       WRITE(logn,*) TRIM(icable_status)
       WRITE(logn,*) ''
    else
       WRITE(logn,*) '.cable_rev file does not exist, suggesting you did not build libcable here.'
       WRITE(logn,*) ''
    endif

    CLOSE(440)

  END SUBROUTINE report_version_no


  SUBROUTINE init_veg_from_vegin(ifmp,fmp, veg)
     use cable_def_types_mod, ONLY : veg_parameter_type
     integer ::  ifmp,  & ! start local mp, # landpoints (jhan:when is this not 1 )
                 fmp     ! local mp, # landpoints

     type(veg_parameter_type) :: veg

     integer :: h

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
            veg%vcmaxcc(h)  = vegin%vcmaxcc(veg%iveg(h))
            veg%ejmaxcc(h)  = vegin%ejmaxcc(veg%iveg(h))
            veg%gmmax(h)    = vegin%gmmax(veg%iveg(h))
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
            veg%froot(h,:)  = vegin%froot(:, veg%iveg(h))
            veg%zr(h)       = vegin%zr(veg%iveg(h))
            veg%clitt(h)    = vegin%clitt(veg%iveg(h))
         END DO ! over each veg patch in land point

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
    IF ( iotype .eq. "dread" ) THEN
       IF ( CABLE_USER%CASA_DUMP_READ )  IS_CASA_TIME = .TRUE.
    ! Check for writing of casa dump output
    ELSE IF ( iotype .eq. "dwrit" ) THEN
       IF ( CABLE_USER%CASA_DUMP_WRITE ) IS_CASA_TIME = .TRUE.
    ! Check for writing of casa standard output
    ELSE IF ( iotype .eq. "write" ) THEN

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


  FUNCTION Esatf(TC)
    !-------------------------------------------------------------------------------
    ! At temperature TC [deg C], return saturation water vapour pressure Esatf [mb]
    ! from Teten formula.
    ! MRR, xx/1987
    ! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
    !                 just like intrinsic functions.
    ! MRR, 12-mar-02: Convert Qsatf (specific humidity routine) to Esatf
    !-------------------------------------------------------------------------------

    IMPLICIT NONE
    REAL, INTENT(in):: TC           ! temp [deg C]
    REAL:: Esatf                    ! saturation vapour pressure [mb]
    REAL:: TCtmp                    ! local
    REAL,PARAMETER:: A = 6.106      ! Teten coefficients
    REAL,PARAMETER:: B = 17.27      ! Teten coefficients
    REAL,PARAMETER:: C = 237.3      ! Teten coefficients
!-------------------------------------------------------------------------------
    TCtmp = TC                          ! preserve TC
    IF (TCtmp.GT.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
    IF (TCtmp.LT.-40.0) TCtmp = -40.0
    Esatf = A*EXP(B*TCtmp/(C+TCtmp))    ! sat vapour pressure (mb)

  END FUNCTION Esatf


END MODULE cable_common_module
