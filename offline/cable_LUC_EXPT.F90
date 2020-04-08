MODULE CABLE_LUC_EXPT

  use cable_common_module,  only: is_leapyear, leap_day, handle_err, get_unit
  use cable_io_vars_module, only: logn, land_x, land_y, landpt, latitude, longitude
  use cable_def_types_mod,  only: mland

  implicit none

  type luc_input_type
     real, dimension(:), allocatable :: val
  end type luc_input_type

  type luc_expt_type
     character(len=400)   :: TransitionFilePath, ClimateFile, Run, NotPrimOnlyFile
     logical              :: DirectRead, READrst, WRITErst
     logical, allocatable :: prim_only(:)
     logical, allocatable :: ptos(:), ptog(:), stog(:), gtos(:)
     logical, allocatable :: ptoc(:), ptoq(:), stoc(:), stoq(:), ctos(:), qtos(:)
     integer, allocatable :: ivegp(:)
     integer, allocatable :: biome(:)
     integer :: YearStart, YearEnd, nfile
     integer :: ctstep
     real,    allocatable :: primaryf(:), grass(:), secdf(:), crop(:), past(:)
     real,    allocatable :: mtemp_min20(:)
     character(len=200),   dimension(17) :: TransFile
     character(len=12) ,   dimension(17) :: var_name
     integer,              dimension(17) :: f_id, v_id
     type(luc_input_type), dimension(17) :: input
     integer :: year, ydimsize, xdimsize, nrec, FirstYear
  end type luc_expt_type
  type(luc_expt_type), save :: luc_expt

  integer, parameter :: &
       ptos      = 1, &
       ptog      = 2, &
       stog      = 3, &
       gtos      = 4, &
       grassfrac = 5, &
       primffrac = 6, &
       pharv     = 7, &
       smharv    = 8, &
       syharv    = 9, &
       ptoc      = 10, &
       ptoq      = 11, &
       stoc      = 12, &
       stoq      = 13, &
       ctos      = 14, &
       qtos      = 15, &
       cropfrac  = 16, &
       pastfrac  = 17

CONTAINS

  ! ------------------------------------------------------------------
  ! CABLE_LUC_EXPT routines
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------

  SUBROUTINE LUC_EXPT_INIT(LUC_EXPT)

    use cable_common_module,       only: cable_user
    use cable_bios_met_obs_params, only: cable_bios_load_biome
    use netcdf,                    only: nf90_open, nf90_nowrite, nf90_inq_varid, nf90_inq_dimid, &
         nf90_inquire_dimension, nf90_inq_varid, nf90_get_att, nf90_get_var, nf90_close

    implicit none

    type(luc_expt_type), intent(inout) :: luc_expt

    REAL :: tmp
    REAL, ALLOCATABLE :: tmparr(:,:)
    INTEGER :: i, k, x
    INTEGER :: fID, timID,latID, lonID, tdimsize, xdimsize, ydimsize
    INTEGER :: xds, yds
    INTEGER :: STATUS,  iu
    CHARACTER(len=15)    :: Run
    CHARACTER(len=400)   :: TransitionFilePath, ClimateFile, NotPrimOnlyFile
    LOGICAL :: DirectRead
    INTEGER :: YearStart, YearEnd
    REAL, ALLOCATABLE :: tmpvec(:), tmparr3(:,:,:), CPC(:)
    INTEGER:: MVG(mland)
    INTEGER :: NotPrimOnly_fID, NotPrimOnly_vID
    INTEGER :: TimeVarID, Idash
    CHARACTER(len=100)    :: time_units
    CHARACTER(len=4) :: yearstr
    REAL :: projection_factor

    NAMELIST /LUCNML/  TransitionFilePath, ClimateFile, Run, DirectRead, YearStart, YearEnd, &
         NotPrimOnlyFile

    ALLOCATE( LUC_EXPT%prim_only(mland) )
    ALLOCATE( LUC_EXPT%ivegp(mland) )
    ALLOCATE( LUC_EXPT%biome(mland) )
    ALLOCATE( LUC_EXPT%ptos(mland) )
    ALLOCATE( LUC_EXPT%ptog(mland) )
    ALLOCATE( LUC_EXPT%stog(mland) )
    ALLOCATE( LUC_EXPT%gtos(mland) )
    ALLOCATE( LUC_EXPT%ptoc(mland) )
    ALLOCATE( LUC_EXPT%ptoq(mland) )
    ALLOCATE( LUC_EXPT%ctos(mland) )
    ALLOCATE( LUC_EXPT%qtos(mland) )
    ALLOCATE( LUC_EXPT%stoc(mland) )
    ALLOCATE( LUC_EXPT%stoq(mland) )
    ALLOCATE( LUC_EXPT%primaryf(mland) )
    ALLOCATE( LUC_EXPT%secdf(mland) )
    ALLOCATE( LUC_EXPT%grass(mland) )
    ALLOCATE( LUC_EXPT%crop(mland) )
    ALLOCATE( LUC_EXPT%past(mland) )
    ALLOCATE( LUC_EXPT%mtemp_min20(mland) )
    !MCINI
    call luc_expt_zero(LUC_EXPT)

    ALLOCATE( CPC(mland))
    LUC_EXPT%NotPrimOnlyFile = 'none'
    ! READ LUC_EXPT settings
    CALL GET_UNIT(iu)
    OPEN(iu,FILE="LUC.nml",STATUS='OLD',ACTION='READ')
    READ(iu,NML=LUCNML)
    CLOSE(iu)
    LUC_EXPT%TransitionFilePath = TransitionFilePath
    LUC_EXPT%ClimateFile        = ClimateFile
    LUC_EXPT%DirectRead         = DirectRead
    LUC_EXPT%YearStart          = YearStart
    LUC_EXPT%YearEnd            = YearEnd
    LUC_EXPT%NotPrimOnlyFile    = NotPrimOnlyFile

    write(*,'(a)') "================== LUC_EXPT  ============"
    write(*,'(a)') "LUC_EXPT settings chosen:"
    write(*,'(a)') " TransitionFilePath: "//trim(LUC_EXPT%TransitionFilePath)
    write(*,'(a)') " ClimateFile       : "//trim(LUC_EXPT%ClimateFile)

    ! Transition Filenames and variables
    LUC_EXPT%TransFile(1) = TRIM(LUC_EXPT%TransitionFilePath)//"/ptos.nc"
    LUC_EXPT%TransFile(2) = TRIM(LUC_EXPT%TransitionFilePath)//"/ptog.nc"
    LUC_EXPT%TransFile(3) = TRIM(LUC_EXPT%TransitionFilePath)//"/stog.nc"
    LUC_EXPT%TransFile(4) = TRIM(LUC_EXPT%TransitionFilePath)//"/gtos.nc"
    LUC_EXPT%TransFile(5) = TRIM(LUC_EXPT%TransitionFilePath)//"/grass.nc"
    LUC_EXPT%TransFile(6) = TRIM(LUC_EXPT%TransitionFilePath)//"/primaryf.nc"
    LUC_EXPT%TransFile(7) = TRIM(LUC_EXPT%TransitionFilePath)//"/pharv.nc"
    LUC_EXPT%TransFile(8) = TRIM(LUC_EXPT%TransitionFilePath)//"/smharv.nc"
    LUC_EXPT%TransFile(9) = TRIM(LUC_EXPT%TransitionFilePath)//"/syharv.nc"

    LUC_EXPT%TransFile(10) = TRIM(LUC_EXPT%TransitionFilePath)//"/ptoc.nc"
    LUC_EXPT%TransFile(11) = TRIM(LUC_EXPT%TransitionFilePath)//"/ptoq.nc"
    LUC_EXPT%TransFile(12) = TRIM(LUC_EXPT%TransitionFilePath)//"/stoc.nc"
    LUC_EXPT%TransFile(13) = TRIM(LUC_EXPT%TransitionFilePath)//"/stoq.nc"
    LUC_EXPT%TransFile(14) = TRIM(LUC_EXPT%TransitionFilePath)//"/ctos.nc"
    LUC_EXPT%TransFile(15) = TRIM(LUC_EXPT%TransitionFilePath)//"/qtos.nc"
    LUC_EXPT%TransFile(16) = TRIM(LUC_EXPT%TransitionFilePath)//"/crop.nc"
    LUC_EXPT%TransFile(17) = TRIM(LUC_EXPT%TransitionFilePath)//"/past.nc"

    LUC_EXPT%VAR_NAME(1) = 'ptos'
    LUC_EXPT%VAR_NAME(2) = 'ptog'
    LUC_EXPT%VAR_NAME(3) = 'stog'
    LUC_EXPT%VAR_NAME(4) = 'gtos'
    LUC_EXPT%VAR_NAME(5) = 'grass'
    LUC_EXPT%VAR_NAME(6) = 'primaryf'
    ! #ifdef __CRU2017__
    !     LUC_EXPT%VAR_NAME(7) = 'pharv'
    ! #else
    !     LUC_EXPT%VAR_NAME(7) = 'ptos'
    ! #endif
    LUC_EXPT%VAR_NAME(7) = 'pharv'
    LUC_EXPT%VAR_NAME(8) = 'smharv'
    LUC_EXPT%VAR_NAME(9) = 'syharv'

    LUC_EXPT%VAR_NAME(10) = 'ptoc'
    LUC_EXPT%VAR_NAME(11) = 'ptoq'
    LUC_EXPT%VAR_NAME(12) = 'stoc'
    LUC_EXPT%VAR_NAME(13) = 'stoq'
    LUC_EXPT%VAR_NAME(14) = 'ctos'
    LUC_EXPT%VAR_NAME(15) = 'qtos'
    LUC_EXPT%VAR_NAME(16) = 'crop'
    LUC_EXPT%VAR_NAME(17) = 'past'

    LUC_EXPT%nfile = 17
    DO x = 1, LUC_EXPT%nfile
       ALLOCATE( LUC_EXPT%INPUT(x)%VAL(mland) )
    END DO

    ! OPEN LUC INPUT FILES
    LUC_EXPT%f_id = -1
    DO i=1, LUC_EXPT%nfile

       write(*,'(a)') 'LUC input data file: '//trim(LUC_EXPT%TransFile(i))
       WRITE(logn,*)  'LUC input data file: ', LUC_EXPT%TransFile(i)

       STATUS = NF90_OPEN(TRIM(LUC_EXPT%TransFile(i)), NF90_NOWRITE, LUC_EXPT%F_ID(i))
       ! print*, 'OOpen30 ', i, LUC_EXPT%nfile, LUC_EXPT%F_ID(i), TRIM(LUC_EXPT%TransFile(i))
       CALL HANDLE_ERR(STATUS, "Opening LUH2 file "//LUC_EXPT%TransFile(i) )
       STATUS = NF90_INQ_VARID(LUC_EXPT%F_ID(i),TRIM(LUC_EXPT%VAR_NAME(i)), LUC_EXPT%V_ID(i))
       CALL HANDLE_ERR(STATUS, "Inquiring LUC_EXPT var "//TRIM(LUC_EXPT%VAR_NAME(i))// &
            " in "//LUC_EXPT%TransFile(i) )

       ! inquire dimensions
       IF (i.eq.1) THEN
          FID = LUC_EXPT%F_ID(i)
          STATUS = NF90_INQ_DIMID(FID,'lat',latID)
          STATUS = NF90_INQUIRE_DIMENSION(FID,latID,len=ydimsize)
          CALL HANDLE_ERR(STATUS, "Inquiring 'lat'"//TRIM(LUC_EXPT%TransFile(i)))
          LUC_EXPT%ydimsize = ydimsize

          STATUS = NF90_INQ_DIMID(FID,'lon',lonID)
          STATUS = NF90_INQUIRE_DIMENSION(FID,lonID,len=xdimsize)
          CALL HANDLE_ERR(STATUS, "Inquiring 'lon'"//TRIM(LUC_EXPT%TransFile(i)))
          LUC_EXPT%xdimsize = xdimsize

          STATUS = NF90_INQ_DIMID(FID,'time',timID)
          STATUS = NF90_INQUIRE_Dimension(FID,timID,len=tdimsize)
          CALL HANDLE_ERR(STATUS, "Inquiring 'time'"//TRIM(LUC_EXPT%TransFile(i)))
          LUC_EXPT%nrec = tdimsize

          ! STATUS = NF90_GET_VAR( Luc_expt%f_id(i), timID, tmp, &
          !      start=(/1,1,1/) )
          ! CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
          STATUS = NF90_INQ_VARID(FID,"time",TimeVarID)
          STATUS = nf90_get_att(FID, TimeVarID, "units", time_units)
          Idash = SCAN(time_units, '-')
          yearstr = time_units(Idash-4:Idash-1)
          read(yearstr,*)  LUC_EXPT%FirstYEAR
          ! write(*,*) 'LUH2 time units: ', TRIM(time_units), Idash, time_units(Idash-4:Idash-1)
          write(*,*) 'LUH2 first year', LUC_EXPT%FirstYEAR
          xds = LUC_EXPT%xdimsize
          yds = LUC_EXPT%ydimsize
       ENDIF
       ! print*, 'OOpened30'
       ! write(*,*) 'length LUH2 data: ', tdimsize
    ENDDO

    write(*,*) 'LUH2 first year', LUC_EXPT%FirstYEAR
    ! LUC_EXPT%FirstYEAR = 850
    ! Set internal counter
    LUC_EXPT%CTSTEP = 1 +  LUC_EXPT%YearStart - LUC_EXPT%FirstYEAR

    allocate(tmparr(xds,yds))

    ! READ initial states
    i = grassfrac
    IF ( LUC_EXPT%DirectRead ) THEN
       DO k = 1, mland
          STATUS = NF90_GET_VAR( LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmp, &
               start=(/land_x(k),land_y(k),LUC_EXPT%CTSTEP/) )
          CALL HANDLE_ERR(STATUS, "Reading direct from "//LUC_EXPT%TransFile(i) )
          LUC_EXPT%grass(k) = tmp
       END DO
    ELSE
       STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr, &
            start=(/1,1,LUC_EXPT%CTSTEP/),count=(/xds,yds,1/) )
       CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
       DO k = 1, mland
          LUC_EXPT%grass(k) = tmparr( land_x(k), land_y(k) )
       END DO
    ENDIF

    i = primffrac
    IF ( LUC_EXPT%DirectRead ) THEN
       DO k = 1, mland
          STATUS = NF90_GET_VAR( LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmp, &
               start=(/land_x(k),land_y(k),LUC_EXPT%CTSTEP/) )
          CALL HANDLE_ERR(STATUS, "Reading direct from "//LUC_EXPT%TransFile(i) )
          LUC_EXPT%primaryf(k) = tmp
       END DO
    ELSE
       STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr, &
            start=(/1,1,LUC_EXPT%CTSTEP/),count=(/xds,yds,1/) )
       CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
       DO k = 1, mland
          LUC_EXPT%primaryf(k) = tmparr( land_x(k), land_y(k) )
       END DO

    ENDIF

    i = cropfrac
    IF ( LUC_EXPT%DirectRead ) THEN
       DO k = 1, mland
          STATUS = NF90_GET_VAR( LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmp, &
               start=(/land_x(k),land_y(k),LUC_EXPT%CTSTEP/) )
          CALL HANDLE_ERR(STATUS, "Reading direct from "//LUC_EXPT%TransFile(i) )
          LUC_EXPT%crop(k) = tmp
       END DO
    ELSE
       STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr, &
            start=(/1,1,LUC_EXPT%CTSTEP/),count=(/xds,yds,1/) )
       CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
       DO k = 1, mland
          LUC_EXPT%crop(k) = tmparr( land_x(k), land_y(k) )
       END DO
    ENDIF

    i = pastfrac
    IF ( LUC_EXPT%DirectRead ) THEN
       DO k = 1, mland
          STATUS = NF90_GET_VAR( LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmp, &
               start=(/land_x(k),land_y(k),LUC_EXPT%CTSTEP/) )
          CALL HANDLE_ERR(STATUS, "Reading direct from "//LUC_EXPT%TransFile(i) )
          LUC_EXPT%past(k) = tmp
       END DO
    ELSE
       STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr, &
            start=(/1,1,LUC_EXPT%CTSTEP/),count=(/xds,yds,1/) )
       CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i))
       DO k = 1, mland
          LUC_EXPT%past(k) = tmparr( land_x(k), land_y(k) )
       END DO
    ENDIF

    LUC_EXPT%grass    = min(LUC_EXPT%grass, 1.0)
    LUC_EXPT%primaryf = min(LUC_EXPT%primaryf, 1.0-LUC_EXPT%grass)
    LUC_EXPT%secdf    = max(1.0-LUC_EXPT%grass-LUC_EXPT%primaryf, 0.0)
    LUC_EXPT%crop     = max(min(LUC_EXPT%crop, LUC_EXPT%grass), 0.0)
    LUC_EXPT%past     = max(min(LUC_EXPT%grass-LUC_EXPT%crop, LUC_EXPT%past), 0.0)

    IF (TRIM(cable_user%MetType) .EQ. "bios") THEN

       ! read bios parameter file to NVIS Major Vegetation Group "biomes"
       CALL cable_bios_load_biome(MVG)
       ! adjust fraction woody cover based on Major Vegetation Group
       LUC_EXPT%biome = MVG
       LUC_EXPT%ivegp = 2
       projection_factor = 0.65
       WHERE (LUC_EXPT%biome .eq. 1)
          CPC = 0.89
       ELSEWHERE (LUC_EXPT%biome .eq. 2)
          CPC = 0.81
       ELSEWHERE (LUC_EXPT%biome .eq. 3)
          CPC = 0.79
       ELSEWHERE (LUC_EXPT%biome .eq. 4)
          CPC = 0.50
       ELSEWHERE (LUC_EXPT%biome .eq. 5)
          CPC = 0.31
       ELSEWHERE (LUC_EXPT%biome .eq. 6)
          CPC = 0.15
       ELSEWHERE (LUC_EXPT%biome .eq. 7)
          CPC = 0.37
       ELSEWHERE (LUC_EXPT%biome .eq. 8)
          CPC = 0.27
       ELSEWHERE (LUC_EXPT%biome .eq. 9)
          CPC = 0.23
       ELSEWHERE (LUC_EXPT%biome .eq. 10)
          CPC = 0.24
       ELSEWHERE (LUC_EXPT%biome .eq. 11)
          CPC = 0.19
       ELSEWHERE (LUC_EXPT%biome .eq. 12)
          CPC = 0.25
       ELSEWHERE (LUC_EXPT%biome .eq. 13)
          CPC = 0.14
       ELSEWHERE (LUC_EXPT%biome .eq. 14)
          CPC = 0.33
       ELSEWHERE (LUC_EXPT%biome .eq. 15)
          CPC = 0.29
       ELSEWHERE (LUC_EXPT%biome .eq. 16)
          CPC = 0.13
       ELSEWHERE (LUC_EXPT%biome .eq. 17)
          CPC = 0.21
       ELSEWHERE (LUC_EXPT%biome .eq. 18)
          CPC = 0.34
       ELSEWHERE (LUC_EXPT%biome .eq. 19)
          CPC = 0.05
       ELSEWHERE (LUC_EXPT%biome .eq. 20)
          CPC = 0.16
       ELSEWHERE (LUC_EXPT%biome .eq. 21)
          CPC = 0.11
       ELSEWHERE (LUC_EXPT%biome .eq. 22)
          CPC = 0.06
       ELSEWHERE (LUC_EXPT%biome .eq. 23)
          CPC = 1.0
       ELSEWHERE (LUC_EXPT%biome .eq. 24)
          CPC = 0.04
       ELSEWHERE (LUC_EXPT%biome .eq. 25)
          CPC= 0.1
       ELSEWHERE (LUC_EXPT%biome .eq. 26)
          CPC= 0.1
       ELSEWHERE (LUC_EXPT%biome .eq. 27)
          CPC = 0.02
       ELSEWHERE (LUC_EXPT%biome .eq. 28)
          CPC = 0.1
       ELSEWHERE (LUC_EXPT%biome .eq. 29)
          CPC = 0.1
       ELSEWHERE (LUC_EXPT%biome .eq. 30)
          CPC = 0.5
       ELSEWHERE (LUC_EXPT%biome .eq. 31)
          CPC= 0.20
       ELSEWHERE (LUC_EXPT%biome .eq. 32)
          CPC = 0.24
       ELSEWHERE
          CPC= 0.1
       ENDWHERE

       CPC = min(CPC/projection_factor, 1.0)

       ! write(*,*)  LUC_EXPT%grass(93), LUC_EXPT%primaryf(93), LUC_EXPT%secdf(93), CPC

       LUC_EXPT%grass    = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf) * (1.0-CPC)
       LUC_EXPT%primaryf = LUC_EXPT%primaryf * CPC
       LUC_EXPT%secdf    = LUC_EXPT%secdf    * CPC
       ! write(*,*)  LUC_EXPT%grass(93), LUC_EXPT%primaryf(93), LUC_EXPT%secdf(93)
    ELSE
       CALL READ_ClimateFile(LUC_EXPT)
       ! hot desert
       WHERE (LUC_EXPT%biome .eq. 15)
          LUC_EXPT%ivegp = 14
       ENDWHERE

       WHERE (LUC_EXPT%biome .eq. 3 .or. LUC_EXPT%biome .eq. 11) ! savanna/ xerophytic woods
          LUC_EXPT%grass    = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf) * 0.6
          LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.4
          LUC_EXPT%secdf    = LUC_EXPT%secdf    * 0.4
       ELSEWHERE (LUC_EXPT%biome .eq. 12 .or. LUC_EXPT%biome .eq. 13 & ! shrub
            .or. LUC_EXPT%biome .eq. 15 .or. LUC_EXPT%biome .eq. 16  )
          LUC_EXPT%grass    = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf) * 0.8
          LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.2
          LUC_EXPT%secdf    = LUC_EXPT%secdf    * 0.2
       ELSEWHERE (LUC_EXPT%biome .eq. 7 .or. LUC_EXPT%biome .eq. 8 &  ! boreal
            .or. LUC_EXPT%biome .eq. 9 .or. LUC_EXPT%biome .eq. 10  )
          LUC_EXPT%grass    = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf) * 0.2
          LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.8
          LUC_EXPT%secdf    = LUC_EXPT%secdf    * 0.8
       ELSEWHERE (LUC_EXPT%biome .eq. 5 .or. LUC_EXPT%biome .eq. 6 ) ! DBL
          LUC_EXPT%grass    = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf) * 0.3
          LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.7
          LUC_EXPT%secdf    = LUC_EXPT%secdf    * 0.7
       END WHERE
    ENDIF

    ! write(59,*) TRIM(LUC_EXPT%NotPrimOnlyFile), (TRIM(LUC_EXPT%NotPrimOnlyFile).EQ.'none')
    ! READ transitions from primary to see if primary remains primary
    if (TRIM(LUC_EXPT%NotPrimOnlyFile).EQ.'none')   THEN
       LUC_EXPT%prim_only = .TRUE.
       IF(.NOT.ALLOCATED(tmpvec)) ALLOCATE(tmpvec(tdimsize))
       IF(.NOT.ALLOCATED(tmparr3)) ALLOCATE(tmparr3(xds,yds,tdimsize))
       DO i=1,2 ! ptos and ptog
          IF ( LUC_EXPT%DirectRead ) THEN
             DO k = 1, mland
                STATUS = NF90_GET_VAR( LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmpvec, &
                     start=(/land_x(k),land_y(k),1/), &
                     count=(/1,1,tdimsize/) )
                CALL HANDLE_ERR(STATUS, "Reading direct from "//LUC_EXPT%TransFile(i) )
                !IF (sum(tmpvec).gt.1e-3 .OR. LUC_EXPT%primaryf(k).lt.0.99) LUC_EXPT%prim_only(k) = .FALSE.
                IF (sum(tmpvec).gt.1e-3) LUC_EXPT%prim_only(k) = .FALSE.
             END DO
          ELSE
             STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr3, &
                  start=(/1,1,1/),count=(/xds,yds,tdimsize/) )
             CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
             DO k = 1, mland
                tmpvec = tmparr3(  land_x(k), land_y(k) , :)
                ! IF (sum(tmpvec).gt.1e-3.OR. LUC_EXPT%primaryf(k).lt.0.99) LUC_EXPT%prim_only(k) = .FALSE.
                IF (sum(tmpvec).gt.1e-3) LUC_EXPT%prim_only(k) = .FALSE.
             END DO
          ENDIF
       END DO
    ELSE
       tmparr = 0.0
       LUC_EXPT%prim_only = .TRUE.
       Status = NF90_OPEN(TRIM(NotPrimOnlyFile), NF90_NOWRITE, NotPrimOnly_fID)
       ! print*, 'OOpen31 ', NotPrimOnly_fID, TRIM(NotPrimOnlyFile)
       CALL HANDLE_ERR(STATUS, "Opening NotPrimOnlyFile"//TRIM(NotPrimOnlyFile ))
       Status = NF90_INQ_VARID( NotPrimOnly_fID,'cum_frac_prim_loss',  NotPrimOnly_vID)
       CALL HANDLE_ERR(STATUS, "Inquiring cum_frac_prim_loss in "//TRIM(NotPrimOnlyFile ) )
       STATUS = NF90_GET_VAR(NotPrimOnly_FID, NotPrimOnly_vID , tmparr, &
            start=(/1,1/),count=(/xds,yds/) )
       CALL HANDLE_ERR(STATUS, "Reading from "//TRIM(NotPrimOnlyFile ) )
       i = 0
       DO k = 1, mland
          if (tmparr( land_x(k), land_y(k)) .gt. 1e-3) then
             LUC_EXPT%prim_only(k) = .FALSE.
             i = i+1
          endif
       ENDDO
       PRINT*, "number of not prim_only grid-cells: ", i
       PRINT*, "number grid-cells: ", mland
       ! print*, 'OClose31 ', NotPrimOnly_fID
       STATUS = NF90_CLOSE(NotPrimOnly_fID)
       CALL HANDLE_ERR(STATUS, "Closing NotPrimOnly "//TRIM(NotPrimOnlyFile))
       NotPrimOnly_fID = -1
    ENDIF

    IF (TRIM(cable_user%MetType) .EQ. "bios") THEN
       WHERE (LUC_EXPT%prim_only .eqv. .TRUE.)
          LUC_EXPT%secdf    = 0.0
          LUC_EXPT%primaryf = 1.0
          LUC_EXPT%grass    = 0.0
          LUC_EXPT%grass    = LUC_EXPT%primaryf * (1.0-CPC)
          LUC_EXPT%primaryf = LUC_EXPT%primaryf * CPC
       ENDWHERE
       DEALLOCATE (CPC)
    ELSE
       ! set secondary vegetation area to be zero where land use transitions don't occur
       ! set grass component of primary vegetation cover
       WHERE (LUC_EXPT%prim_only .eqv. .TRUE.)
          LUC_EXPT%secdf    = 0.0
          LUC_EXPT%primaryf = 1.0
          LUC_EXPT%grass    = 0.0
          WHERE (LUC_EXPT%biome .eq. 3 .or. LUC_EXPT%biome .eq. 11) ! savanna/ xerophytic woods
             LUC_EXPT%grass    = LUC_EXPT%primaryf * 0.6
             LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.4
          ELSEWHERE (LUC_EXPT%biome .eq. 12 .or. LUC_EXPT%biome .eq. 13 &
               .or. LUC_EXPT%biome .eq. 15 .or. LUC_EXPT%biome .eq. 16  ) ! shrub
             LUC_EXPT%grass    = LUC_EXPT%primaryf * 0.8
             LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.2
          ELSEWHERE (LUC_EXPT%biome .eq. 7 .or. LUC_EXPT%biome .eq. 8 &
               .or. LUC_EXPT%biome .eq. 9 .or. LUC_EXPT%biome .eq. 10) ! boreal
             LUC_EXPT%grass    = LUC_EXPT%primaryf * 0.2
             LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.8
          ELSEWHERE (LUC_EXPT%biome .eq. 5 .or. LUC_EXPT%biome .eq. 6 ) ! DBL
             LUC_EXPT%grass    = LUC_EXPT%primaryf * 0.3
             LUC_EXPT%primaryf = LUC_EXPT%primaryf * 0.7
          END WHERE
       END WHERE
    ENDIF
    ! print*, 'ORead30 '

  END SUBROUTINE LUC_EXPT_INIT

  subroutine luc_expt_zero(luc_expt)

    implicit none

    type(luc_expt_type), intent(inout) :: luc_expt

    luc_expt%prim_only   = .false.
    luc_expt%ivegp       = 0
    luc_expt%biome       = 0
    luc_expt%ptos        = .false.
    luc_expt%ptog        = .false.
    luc_expt%stog        = .false.
    luc_expt%gtos        = .false.
    luc_expt%ptoc        = .false.
    luc_expt%ptoq        = .false.
    luc_expt%ctos        = .false.
    luc_expt%qtos        = .false.
    luc_expt%stoc        = .false.
    luc_expt%stoq        = .false.
    luc_expt%primaryf    = 0.
    luc_expt%secdf       = 0.
    luc_expt%grass       = 0.
    luc_expt%crop        = 0.
    luc_expt%past        = 0.
    luc_expt%mtemp_min20 = 0.

  end subroutine luc_expt_zero
  
  ! ------------------------------------------------------------------


  SUBROUTINE LUC_EXPT_SET_TILES(inVeg, inPfrac, LUC_EXPT)

    IMPLICIT NONE

    INTEGER,             INTENT(INOUT) :: inVeg(:,:,:)
    REAL,                INTENT(INOUT) :: inPFrac(:,:,:)
    TYPE(LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT
    
    INTEGER :: k, m, n

    DO k=1, mland
       m = landpt(k)%ilon
       n = landpt(k)%ilat

       if (inVeg(m,n,1).LT.11) THEN ! vegetated

          if (LUC_EXPT%prim_only(k) ) then

             inVeg(m,n,1)     = LUC_EXPT%ivegp(k)
             inVeg(m,n,2:3)   = 0
             inPFrac(m,n,2:3) = 0.0
             inPFrac(m,n,1)   = 1.0
             if ( LUC_EXPT%grass(k) .gt. 0.01) then
                IF (LUC_EXPT%mtemp_min20(k) .LE. 15.5) THEN
                   inVeg(m,n,2) = 6 ! C3 grass
                ELSE
                   inVeg(m,n,2) = 7 ! C4 grass
                ENDIF
                inPFrac(m,n,1) = min(LUC_EXPT%primaryf(k),1.0)
                inPFrac(m,n,2) = 1.0 - min(LUC_EXPT%primaryf(k),1.0)
             endif

          elseif ((.NOT.LUC_EXPT%prim_only(k)) ) then
             
             inVeg(m,n,1) = LUC_EXPT%ivegp(k)
             inVeg(m,n,2) = LUC_EXPT%ivegp(k)

             if (LUC_EXPT%mtemp_min20(k) .LE. 15.5) THEN
                inVeg(m,n,3) = 6 ! C3 grass
             ELSE
                inVeg(m,n,3) = 7 ! C4 grass
             ENDIF
             inPFrac(m,n,1) = min(LUC_EXPT%primaryf(k),1.0)
             inPFrac(m,n,2) = min(LUC_EXPT%secdf(k),1.0)
             inPFrac(m,n,3) = 1.0 - inPFrac(m,n,1) - inPFrac(m,n,2)

          endif
       else

          LUC_EXPT%prim_only(k)=.TRUE.

       endif

       ! don't consider LUC events in desert or tundra
       if (inveg(m,n,1)==14 .OR.  inveg(m,n,1)==8 ) THEN
          LUC_EXPT%prim_only(k) = .TRUE.
          LUC_EXPT%primaryf(k)  = 1.0
          LUC_EXPT%secdf(k)     = 0.0
          LUC_EXPT%grass(k)     = 0.0
          inPFrac(m,n,1)   = 1.0
          inPFrac(m,n,2:3) = 0.0
          inVeg(m,n,2:3)   = 0
       endif

    ENDDO

  END SUBROUTINE LUC_EXPT_SET_TILES


  ! ------------------------------------------------------------------


  SUBROUTINE LUC_EXPT_SET_TILES_BIOS(inVeg, inPfrac, LUC_EXPT )

    USE cable_bios_met_obs_params, ONLY: cable_bios_load_fracC4

    IMPLICIT NONE

    INTEGER,              INTENT(INOUT) :: inVeg(:,:,:)
    REAL,                 INTENT(INOUT) :: inPFrac(:,:,:)
    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT
    
    REAL :: fracC4(mland)
    INTEGER :: k, m, n

    CALL cable_bios_load_fracC4(fracC4)

    DO k=1, mland
       m = landpt(k)%ilon
       n = landpt(k)%ilat

       if (LUC_EXPT%prim_only(k) ) then

          inVeg(m,n,1)     = LUC_EXPT%ivegp(k)
          inVeg(m,n,2:3)   = 0
          inPFrac(m,n,2:3) = 0.0
          inPFrac(m,n,1)   = 1.0
          if ( LUC_EXPT%grass(k) .gt. 0.01 ) then
             if (fracC4(k).gt.0.5) then
                inVeg(m,n,2) = 7 ! C4 grass
             else
                inVeg(m,n,2) = 6 ! C3 grass
             endif
             inPFrac(m,n,1) = min(LUC_EXPT%primaryf(k),1.0)
             inPFrac(m,n,2) = 1.0 - min(LUC_EXPT%primaryf(k),1.0)
          endif

       elseif ((.NOT.LUC_EXPT%prim_only(k)) ) then
          
          inVeg(m,n,1) = LUC_EXPT%ivegp(k)
          inVeg(m,n,2) = LUC_EXPT%ivegp(k)
          if (fracC4(k).gt.0.5) then
             inVeg(m,n,3) = 7 ! C4 grass
          else
             inVeg(m,n,3) = 6 ! C3 grass
          endif
          inPFrac(m,n,1) = min(LUC_EXPT%primaryf(k),1.0)
          inPFrac(m,n,2) = min(LUC_EXPT%secdf(k),1.0)
          inPFrac(m,n,3) = 1.0 - inPFrac(m,n,1) - inPFrac(m,n,2)

       endif
       
    ENDDO

  END SUBROUTINE LUC_EXPT_SET_TILES_BIOS


  ! ------------------------------------------------------------------


  SUBROUTINE READ_ClimateFile(LUC_EXPT)

    use netcdf, only: nf90_open, nf90_nowrite, nf90_inq_varid, nf90_inq_dimid, &
         nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_close, nf90_noerr
#ifdef __MPI__
    use mpi,    only: MPI_Abort
#endif

    IMPLICIT NONE

    TYPE(LUC_EXPT_type), INTENT(INOUT) :: LUC_EXPT ! climate variables

    INTEGER(KIND=4) :: mp4
    INTEGER(KIND=4), parameter :: pmp4 = 0
    INTEGER, parameter :: fmp4 = kind(pmp4)
    INTEGER(KIND=4)    :: STATUS
    INTEGER(KIND=4)    :: FILE_ID, dID, i, land_dim
    CHARACTER :: FNAME*99

    ! 0 dim arrays
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(3) :: A1
    ! 1 dim arrays (integer) (npt )
    CHARACTER(len=20),DIMENSION(2) :: AI1

    ! REAL(r_2), DIMENSION(mland) :: LAT, LON, TMP
    REAL, DIMENSION(mland) :: LAT, LON, TMP
    INTEGER(KIND=4) :: TMPI(mland)
    LOGICAL         :: EXISTFILE
#ifdef __MPI__
    integer :: ierr
#endif

    mp4 = int(mland,fmp4)
    A0(1) = 'nyears'
    A0(2) = 'year'

    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'mtemp_min20'

    AI1(1) = 'iveg'
    AI1(2) = 'biome'

    fname = TRIM(LUC_EXPT%ClimateFile)

    INQUIRE( FILE=TRIM( fname ), EXIST=EXISTFILE )

    IF ( .NOT.EXISTFILE) THEN
       write(*,'(a)') trim(fname)//' does not exist!!'
    ELSE
       write(*,'(a)') 'reading biome from : '//trim(fname)
    ENDIF
    ! Open NetCDF file:
    STATUS = NF90_OPEN(TRIM(fname), NF90_NOWRITE, FILE_ID)
    ! print*, 'OOpen32 ', file_id, TRIM(fname)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! dimensions:
    ! Land (number of points)


    STATUS = NF90_INQ_DIMID(FILE_ID, 'land'   , dID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    IF ( land_dim .NE. mland) THEN
       WRITE(*,*) "Dimension misfit, ", fname
       WRITE(*,*) "land_dim", land_dim
#ifdef __MPI__
       call MPI_Abort(0, 4, ierr) ! Do not know comm nor rank here
#else
       stop 4
#endif
    ENDIF

    ! LAT & LON
    STATUS = NF90_INQ_VARID( FILE_ID, A1(1), dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_GET_VAR( FILE_ID, dID, LAT )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_VARID( FILE_ID, A1(2), dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_GET_VAR( FILE_ID, dID, LON )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


    ! READ 1-dimensional real fields
    DO i = 3, SIZE(A1)
       STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE (TRIM(A1(i)))
       CASE ('mtemp_min20')
          LUC_EXPT%mtemp_min20 = TMP
       END SELECT
    END DO

    ! READ 1-dimensional integer fields
    DO i = 1, SIZE(AI1)
       write(*,*)  TRIM(AI1(i))
       STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(AI1(i)))
       CASE ('iveg')
          LUC_EXPT%ivegp = TMPI
       CASE ('biome')
          LUC_EXPT%biome = TMPI
       END SELECT
    END DO

    ! non-woody potential vegetation not considered to undergo LU change
    WHERE (LUC_EXPT%ivegp .GT. 5)
       LUC_EXPT%prim_only=.TRUE.
    ENDWHERE

    ! Close NetCDF file:
    ! print*, 'OClose32 ', file_id
    STATUS  = NF90_close(FILE_ID)
    file_id = -1
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    write(*,'(a)') "end read_climatefile"

  END SUBROUTINE READ_CLIMATEFILE


  ! ------------------------------------------------------------------

  
  SUBROUTINE READ_LUH2(LUC_EXPT)

    use netcdf, only: nf90_get_var

    IMPLICIT NONE

    TYPE(LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT

    REAL    ::  tmp
    REAL, ALLOCATABLE :: tmparr(:,:)
    INTEGER :: t, i, k
    INTEGER :: xds, yds
    INTEGER :: STATUS

    xds = LUC_EXPT%xdimsize
    yds = LUC_EXPT%ydimsize
    t   = LUC_EXPT%CTSTEP
    IF (.NOT.ALLOCATED(tmparr)) ALLOCATE(tmparr(xds,yds))

    if (t .LE. LUC_EXPT%nrec) then
       DO i=1, LUC_EXPT%nfile
          ! print*, 'ORead30 ', i, LUC_EXPT%nfile, LUC_EXPT%F_ID(i)
          IF ( LUC_EXPT%DirectRead ) THEN
             DO k = 1, mland
                STATUS = NF90_GET_VAR( LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmp, &
                     start=(/land_x(k),land_y(k),t/) )
                CALL HANDLE_ERR(STATUS, "Reading direct from "//LUC_EXPT%TransFile(i) )
                LUC_EXPT%INPUT(i)%VAL(k) = tmp
             END DO
          ELSE
             STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr, &
                  start=(/1,1,t/),count=(/xds,yds,1/) )
             CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
             DO k = 1, mland
                LUC_EXPT%INPUT(i)%VAL(k) = tmparr( land_x(k), land_y(k) )
                if (LUC_EXPT%INPUT(i)%VAL(k) .gt. 1.0) then
                   LUC_EXPT%INPUT(i)%VAL(k) = 0.0
                endif
             END DO
          ENDIF
       ENDDO
    else
       write(*,'(a)') 'warning: past end of LUH2 record'
    endif

    ! Adjust transition areas based on primary wooded fraction
    WHERE (LUC_EXPT%biome .eq. 3 .or. LUC_EXPT%biome .eq. 11)  ! savanna/ xerophytic woods
       LUC_EXPT%INPUT(ptos)%VAL   =  LUC_EXPT%INPUT(ptos)%VAL   * 0.4
       LUC_EXPT%INPUT(ptog)%VAL   =  LUC_EXPT%INPUT(ptog)%VAL   * 0.4
       LUC_EXPT%INPUT(gtos)%VAL   =  LUC_EXPT%INPUT(gtos)%VAL   * 0.4
       LUC_EXPT%INPUT(stog)%VAL   =  LUC_EXPT%INPUT(stog)%VAL   * 0.4
       LUC_EXPT%INPUT(smharv)%VAL =  LUC_EXPT%INPUT(smharv)%VAL * 0.4
       LUC_EXPT%INPUT(syharv)%VAL =  LUC_EXPT%INPUT(syharv)%VAL * 0.4
    ELSEWHERE (LUC_EXPT%biome .eq. 12 .or. LUC_EXPT%biome .eq. 13 &
         .or. LUC_EXPT%biome .eq. 15 .or. LUC_EXPT%biome .eq. 16  ) ! shrub
       LUC_EXPT%INPUT(ptos)%VAL   =  LUC_EXPT%INPUT(ptos)%VAL   * 0.2
       LUC_EXPT%INPUT(ptog)%VAL   =  LUC_EXPT%INPUT(ptog)%VAL   * 0.2
       LUC_EXPT%INPUT(gtos)%VAL   =  LUC_EXPT%INPUT(gtos)%VAL   * 0.2
       LUC_EXPT%INPUT(stog)%VAL   =  LUC_EXPT%INPUT(stog)%VAL   * 0.2
       LUC_EXPT%INPUT(smharv)%VAL =  LUC_EXPT%INPUT(smharv)%VAL * 0.2
       LUC_EXPT%INPUT(syharv)%VAL =  LUC_EXPT%INPUT(syharv)%VAL * 0.2
    ELSEWHERE (LUC_EXPT%biome .eq. 7 .or. LUC_EXPT%biome .eq. 8 &
         .or. LUC_EXPT%biome .eq. 9 .or. LUC_EXPT%biome .eq. 10) ! boreal
       LUC_EXPT%INPUT(ptos)%VAL   =  LUC_EXPT%INPUT(ptos)%VAL   * 0.8
       LUC_EXPT%INPUT(ptog)%VAL   =  LUC_EXPT%INPUT(ptog)%VAL   * 0.8
       LUC_EXPT%INPUT(gtos)%VAL   =  LUC_EXPT%INPUT(gtos)%VAL   * 0.8
       LUC_EXPT%INPUT(stog)%VAL   =  LUC_EXPT%INPUT(stog)%VAL   * 0.8
       LUC_EXPT%INPUT(smharv)%VAL =  LUC_EXPT%INPUT(smharv)%VAL * 0.8
       LUC_EXPT%INPUT(syharv)%VAL =  LUC_EXPT%INPUT(syharv)%VAL * 0.8
    ELSEWHERE (LUC_EXPT%biome .eq. 5 .or. LUC_EXPT%biome .eq. 6 ) ! DBL
       LUC_EXPT%INPUT(ptos)%VAL   =  LUC_EXPT%INPUT(ptos)%VAL   * 0.7
       LUC_EXPT%INPUT(ptog)%VAL   =  LUC_EXPT%INPUT(ptog)%VAL   * 0.7
       LUC_EXPT%INPUT(gtos)%VAL   =  LUC_EXPT%INPUT(gtos)%VAL   * 0.7
       LUC_EXPT%INPUT(stog)%VAL   =  LUC_EXPT%INPUT(stog)%VAL   * 0.7
       LUC_EXPT%INPUT(smharv)%VAL =  LUC_EXPT%INPUT(smharv)%VAL * 0.7
       LUC_EXPT%INPUT(syharv)%VAL =  LUC_EXPT%INPUT(syharv)%VAL * 0.7
    ENDWHERE

  END SUBROUTINE READ_LUH2

  
  ! ------------------------------------------------------------------

  
  SUBROUTINE close_luh2(LUC_EXPT)

    use netcdf, only: nf90_close

    implicit none

    type(luc_expt_type), intent(inout) :: LUC_EXPT

    integer :: i, status

    write(*,*) 'Closing LUH2 files.'
    do i=1, LUC_EXPT%nfile
       if (LUC_EXPT%F_ID(i) > -1) then
          ! print*, 'OClose30 ', i, LUC_EXPT%nfile, LUC_EXPT%F_ID(i)
          status = nf90_close(LUC_EXPT%F_ID(i))
          LUC_EXPT%F_ID(i) = -1
          call handle_err(status, "Closing LUH2 transition file "//trim(LUC_EXPT%TransFile(i)))
       end if
    end do

  END SUBROUTINE close_luh2

  ! ------------------------------------------------------------------

END MODULE CABLE_LUC_EXPT
