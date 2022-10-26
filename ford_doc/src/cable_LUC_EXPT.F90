MODULE CABLE_LUC_EXPT

  USE netcdf
  USE casa_ncdf_module, ONLY: HANDLE_ERR, GET_UNIT
  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, LEAP_DAY

  USE cable_IO_vars_module, ONLY: logn, land_x, land_y, landpt, latitude, longitude
  USE cable_def_types_mod,  ONLY: mland, r_2

  IMPLICIT NONE

  TYPE LUC_INPUT_TYPE
     REAL, DIMENSION(:), ALLOCATABLE :: VAL
     !  INTEGER :: YEAR
  END TYPE LUC_INPUT_TYPE

  TYPE LUC_EXPT_TYPE
     CHARACTER(len=200)::	TransitionFilePath,ClimateFile, Run
     LOGICAL  :: DirectRead, READrst, WRITErst
     LOGICAL, ALLOCATABLE ::  prim_only(:)
     LOGICAL, ALLOCATABLE :: ptos(:), ptog(:), stog(:), gtos(:)
     INTEGER, ALLOCATABLE :: ivegp(:)
     INTEGER, ALLOCATABLE :: biome(:)
     INTEGER :: YearStart, YearEnd, nfile
     INTEGER :: CTSTEP
     REAL, ALLOCATABLE :: primaryf(:), mtemp_min20(:), grass(:), secdf(:)
     CHARACTER(len=200),DIMENSION(9) :: TransFile
     CHARACTER(len=12) ,DIMENSION(9) :: VAR_NAME
     INTEGER,           DIMENSION(9) :: F_ID, V_ID
     TYPE(LUC_INPUT_TYPE), DIMENSION(9):: INPUT
     INTEGER :: YEAR, ydimsize, xdimsize, nrec, FirstYear

  END TYPE LUC_EXPT_TYPE

  TYPE (LUC_EXPT_TYPE), SAVE :: LUC_EXPT

  INTEGER,  PARAMETER :: &
       ptos         =  1, &
       ptog         =  2, &
       stog         =  3, &
       gtos         =  4, &
       grassfrac    =  5, &
       primffrac    =  6, &
       pharv        =  7, &
       smharv       =  8, &
       syharv       =  9


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALBE_LUC_EXPT routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ==============================================================================
  SUBROUTINE LUC_EXPT_INIT( LUC_EXPT)

    IMPLICIT NONE

    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT

    REAL    :: tmparr(720,360), tmp
    INTEGER :: t, i, ii, k, x, y, realk
    INTEGER :: fID, vID, timID,latID, lonID, tdimsize, xdimsize, ydimsize
    INTEGER :: xds, yds, tds
    INTEGER :: STATUS,  iu
    CHARACTER(len=15)    :: Run
    CHARACTER(len=200)   :: TransitionFilePath, ClimateFile
    LOGICAL :: DirectRead
    INTEGER :: YearStart, YearEnd
    REAL, ALLOCATABLE :: tmpvec(:), tmparr3(:,:,:)

    NAMELIST /LUCNML/  TransitionFilePath, ClimateFile, Run, DirectRead, YearStart, YearEnd

    ALLOCATE( LUC_EXPT%prim_only(mland) )
    ALLOCATE( LUC_EXPT%ivegp(mland) )
    ALLOCATE( LUC_EXPT%biome(mland) )
    ALLOCATE( LUC_EXPT%ptos(mland) )
    ALLOCATE( LUC_EXPT%ptog(mland) )
    ALLOCATE( LUC_EXPT%stog(mland) )
    ALLOCATE( LUC_EXPT%gtos(mland) )
    ALLOCATE( LUC_EXPT%primaryf(mland) )
    ALLOCATE( LUC_EXPT%secdf(mland) )
    ALLOCATE( LUC_EXPT%grass(mland) )
    ALLOCATE( LUC_EXPT%mtemp_min20(mland) )



    ! READ LUC_EXPT settings
    CALL GET_UNIT(iu)
    OPEN (iu,FILE="LUC.nml",STATUS='OLD',ACTION='READ')
    READ (iu,NML=LUCNML)
    CLOSE(iu)
    LUC_EXPT%TransitionFilePath = TransitionFilePath
    LUC_EXPT%ClimateFile        = ClimateFile
    LUC_EXPT%DirectRead         = DirectRead
    LUC_EXPT%YearStart          = YearStart
    LUC_EXPT%YearEnd            = YearEnd

    WRITE(*   ,*)"================== LUC_EXPT  ============"
    WRITE(*   ,*)"LUC_EXPT settings chosen:"
    WRITE(*   ,*)" TransitionFilePath: ",TRIM(LUC_EXPT%TransitionFilePath)
    WRITE(*   ,*)" ClimateFile       : ",TRIM(LUC_EXPT%ClimateFile)


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

    LUC_EXPT%VAR_NAME(1) = 'ptos'
    LUC_EXPT%VAR_NAME(2) = 'ptog'
    LUC_EXPT%VAR_NAME(3) = 'stog'
    LUC_EXPT%VAR_NAME(4) = 'gtos'
    LUC_EXPT%VAR_NAME(5) = 'grass'
    LUC_EXPT%VAR_NAME(6) = 'primaryf'
    LUC_EXPT%VAR_NAME(7) = 'pharv'
    LUC_EXPT%VAR_NAME(8) = 'smharv'
    LUC_EXPT%VAR_NAME(9) = 'syharv'

    LUC_EXPT%nfile = 9

    DO x = 1, LUC_EXPT%nfile
       ALLOCATE( LUC_EXPT%INPUT(x)%VAL(mland) )
    END DO

    ! OPEN LUC INPUT FILES
    DO i = 1, LUC_EXPT%nfile

       WRITE(*   ,*) 'LUC input data file: ', LUC_EXPT%TransFile(i)
       WRITE(logn,*) 'LUC input data file: ', LUC_EXPT%TransFile(i)

       STATUS = NF90_OPEN(TRIM(LUC_EXPT%TransFile(i)), NF90_NOWRITE, LUC_EXPT%F_ID(i))
       CALL HANDLE_ERR(STATUS, "Opening LUH2 file "//LUC_EXPT%TransFile(i) )
       STATUS = NF90_INQ_VARID(LUC_EXPT%F_ID(i),TRIM(LUC_EXPT%VAR_NAME(i)), LUC_EXPT%V_ID(i))
       CALL HANDLE_ERR(STATUS, "Inquiring LUC_EXPT var "//TRIM(LUC_EXPT%VAR_NAME(i))// &
            " in "//LUC_EXPT%TransFile(i) )

       ! inquire dimensions
       IF (i.EQ.1) THEN
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
          STATUS = NF90_INQUIRE_DIMENSION(FID,timID,len=tdimsize)
          CALL HANDLE_ERR(STATUS, "Inquiring 'time'"//TRIM(LUC_EXPT%TransFile(i)))
          LUC_EXPT%nrec = tdimsize

!!$          STATUS = NF90_GET_VAR( Luc_expt%f_id(i), timID, tmp, &
!!$               start=(/1,1,1/) )
!!$          CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )



          xds = LUC_EXPT%xdimsize
          yds = LUC_EXPT%ydimsize
       ENDIF
       !write(*,*) 'length LUH2 data: ', tdimsize
    ENDDO



    LUC_EXPT%FirstYEAR = 850
    ! Set internal counter
    LUC_EXPT%CTSTEP = 1 +  LUC_EXPT%YearStart- LUC_EXPT%FirstYEAR

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



    LUC_EXPT%grass = MIN(LUC_EXPT%grass, 1.0)
    LUC_EXPT%primaryf = MIN(LUC_EXPT%primaryf, 1.0- LUC_EXPT%grass)
    LUC_EXPT%secdf = MAX((1.0 -  LUC_EXPT%grass - LUC_EXPT%primaryf), 0.0)

    CALL READ_ClimateFile(LUC_EXPT)
    ! hot desert
    WHERE (LUC_EXPT%biome.EQ.15 )
       LUC_EXPT%ivegp = 14
    ENDWHERE

    WHERE (LUC_EXPT%biome .EQ. 3 .OR. LUC_EXPT%biome .EQ. 11) ! savanna/ xerophytic woods
       LUC_EXPT%grass = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf)*1.0/2.0
       LUC_EXPT%primaryf =  LUC_EXPT%primaryf * 1.0/2.0
       LUC_EXPT%secdf =  LUC_EXPT%secdf * 1.0/2.0
    ELSEWHERE (LUC_EXPT%biome .EQ. 12 .OR. LUC_EXPT%biome .EQ. 13 & ! shrub
         .OR. LUC_EXPT%biome .EQ. 15 .OR. LUC_EXPT%biome .EQ. 16  )
       LUC_EXPT%grass = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf)*4.0/5.0
       LUC_EXPT%primaryf =  LUC_EXPT%primaryf * 1.0/5.0
       LUC_EXPT%secdf =  LUC_EXPT%secdf * 1.0/5.0
    ELSEWHERE (LUC_EXPT%biome .EQ. 7 .OR. LUC_EXPT%biome .EQ. 8 &  ! boreal
         .OR. LUC_EXPT%biome .EQ. 9 .OR. LUC_EXPT%biome .EQ. 10  )
       LUC_EXPT%grass = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf)*1.0/5.0
       LUC_EXPT%primaryf =  LUC_EXPT%primaryf * 4.0/5.0
       LUC_EXPT%secdf =  LUC_EXPT%secdf * 4.0/5.0
    ELSEWHERE (LUC_EXPT%biome .EQ. 5 .OR. LUC_EXPT%biome .EQ. 6 ) ! DBL
       LUC_EXPT%grass = LUC_EXPT%grass + (LUC_EXPT%primaryf+LUC_EXPT%secdf)*0.3
       LUC_EXPT%primaryf =  LUC_EXPT%primaryf *0.7
       LUC_EXPT%secdf =  LUC_EXPT%secdf * 0.7
    END WHERE




    ! READ transitions from primary to see if primary remains primary

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
             IF (SUM(tmpvec).GT.1e-3 ) LUC_EXPT%prim_only(k) = .FALSE.
          END DO

       ELSE

          STATUS = NF90_GET_VAR(LUC_EXPT%F_ID(i), LUC_EXPT%V_ID(i), tmparr3, &
               start=(/1,1,1/),count=(/xds,yds,tdimsize/) )
          CALL HANDLE_ERR(STATUS, "Reading from "//LUC_EXPT%TransFile(i) )
          DO k = 1, mland
             tmpvec = tmparr3(  land_x(k), land_y(k) , :)
             ! IF (sum(tmpvec).gt.1e-3.OR. LUC_EXPT%primaryf(k).lt.0.99) LUC_EXPT%prim_only(k) = .FALSE.
             IF (SUM(tmpvec).GT.1e-3) LUC_EXPT%prim_only(k) = .FALSE.
          END DO

       ENDIF

    END DO

    ! set secondary vegetation area to be zero where land use transitions don't occur
    ! set grass component of primary vegetation cover
    WHERE (LUC_EXPT%prim_only .EQV. .TRUE.)
       LUC_EXPT%secdf = 0.0
       LUC_EXPT%primaryf = 1.0
       LUC_EXPT%grass = 0.0
       WHERE (LUC_EXPT%biome .EQ. 3 .OR. LUC_EXPT%biome .EQ. 11) ! savanna/ xerophytic woods
          LUC_EXPT%grass = LUC_EXPT%primaryf*1.0/2.0
          LUC_EXPT%primaryf =  LUC_EXPT%primaryf * 1.0/2.0
       ELSEWHERE (LUC_EXPT%biome .EQ. 12 .OR. LUC_EXPT%biome .EQ. 13 &
            .OR. LUC_EXPT%biome .EQ. 15 .OR. LUC_EXPT%biome .EQ. 16  ) ! shrub
          LUC_EXPT%grass = LUC_EXPT%primaryf*4.0/5.0
          LUC_EXPT%primaryf =  LUC_EXPT%primaryf * 1.0/5.0
       ELSEWHERE (LUC_EXPT%biome .EQ. 7 .OR. LUC_EXPT%biome .EQ. 8 &
            .OR. LUC_EXPT%biome .EQ. 9 .OR. LUC_EXPT%biome .EQ. 10) ! boreal
          LUC_EXPT%grass = LUC_EXPT%primaryf*1.0/5.0
          LUC_EXPT%primaryf =  LUC_EXPT%primaryf * 4.0/5.0
       ELSEWHERE (LUC_EXPT%biome .EQ. 5 .OR. LUC_EXPT%biome .EQ. 6 ) ! DBL
          LUC_EXPT%grass = LUC_EXPT%primaryf*0.3
          LUC_EXPT%primaryf =  LUC_EXPT%primaryf *0.7
       END WHERE
    END WHERE


!!$    WHERE (LUC_EXPT%ivegp == 14)
!!$       LUC_EXPT%prim_only = .TRUE.
!!$    END WHERE

  END SUBROUTINE LUC_EXPT_INIT

  ! ==============================================================================


  SUBROUTINE LUC_EXPT_SET_TILES(inVeg, inPfrac, LUC_EXPT )

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: inVeg(:,:,:)
    REAL,   INTENT(INOUT) :: inPFrac(:,:,:)
    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT
    INTEGER :: k, m, n




    DO k=1,mland
       m = landpt(k)%ilon
       n = landpt(k)%ilat


       IF (inVeg(m,n,1).LT.11) THEN ! vegetated

          IF (LUC_EXPT%prim_only(k) ) THEN

             inVeg(m,n,1) = LUC_EXPT%ivegp(k)
             inVeg(m,n,2:3) = 0
             inPFrac(m,n,2:3) = 0
             inPFrac(m,n,1) = 1.0
             IF ( LUC_EXPT%grass(k) .GT. 0.01) THEN
                IF (LUC_EXPT%mtemp_min20(k).LE. 15.5) THEN
                   inVeg(m,n,2) = 6 ! C3 grass
                ELSE
                   inVeg(m,n,2) = 7 ! C4 grass
                ENDIF
                inPFrac(m,n,1) =  MIN(LUC_EXPT%primaryf(k),1.0)
                inPFrac(m,n,2) = 1.0 -  MIN(LUC_EXPT%primaryf(k),1.0)
             ENDIF

          ELSEIF ((.NOT.LUC_EXPT%prim_only(k)) ) THEN
             inVeg(m,n,1) = LUC_EXPT%ivegp(k)
             inVeg(m,n,2) = LUC_EXPT%ivegp(k)

             IF (LUC_EXPT%mtemp_min20(k).LE. 15.5) THEN
                inVeg(m,n,3) = 6 ! C3 grass
             ELSE
                inVeg(m,n,3) = 7 ! C4 grass
             ENDIF
             inPFrac(m,n,1) = MIN(LUC_EXPT%primaryf(k),1.0)
             inPFrac(m,n,2) = MIN(LUC_EXPT%secdf(k),1.0)
             inPFrac(m,n,3) = 1.0 - inPFrac(m,n,1) - inPFrac(m,n,2)


          ENDIF
       ELSE
          LUC_EXPT%prim_only(k)=.TRUE.

       ENDIF

       ! don't consider LUC events in desert or tundra
       IF (inveg(m,n,1)==14 .OR.  inveg(m,n,1)==8 ) THEN
          LUC_EXPT%prim_only(k)=.TRUE.
          LUC_EXPT%primaryf(k) = 1.0
          LUC_EXPT%secdf(k) = 0.0
          LUC_EXPT%grass(k) = 0.0
          inPFrac(m,n,1) = 1.0
          inPFrac(m,n,2:3) = 0.0
          inVeg(m,n,2:3) = 0
       ENDIF


    ENDDO



991 FORMAT(1166(e12.4,2x))

  END SUBROUTINE LUC_EXPT_SET_TILES
  ! ==============================================================================

  SUBROUTINE READ_ClimateFile(LUC_EXPT)

    USE netcdf


    IMPLICIT NONE

    TYPE (LUC_EXPT_type), INTENT(INOUT)       :: LUC_EXPT ! climate variables

    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS
    INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID, dID, i, land_dim
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

    ! 0 dim arrays
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(3) :: A1
    ! 1 dim arrays (integer) (npt )
    CHARACTER(len=20),DIMENSION(2) :: AI1

    REAL(r_2), DIMENSION(mland)          :: LAT, LON, TMP
    INTEGER*4 :: TMPI(mland), TMPI0
    LOGICAL            ::  EXISTFILE

    mp4=INT(mland,fmp4)
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
       WRITE(*,*) fname, ' does not exist!!'
    ELSE
       WRITE(*,*) 'reading biome from : ', fname
    ENDIF
    ! Open NetCDF file:
    STATUS = NF90_OPEN(fname, NF90_NOWRITE, FILE_ID)
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
       STOP
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

       SELECT CASE ( TRIM(A1(i)))
       CASE ('mtemp_min20' ) ; LUC_EXPT%mtemp_min20 = TMP
       END SELECT
    END DO

    ! READ 1-dimensional integer fields
    DO i = 1, SIZE(AI1)

       WRITE(*,*)  TRIM(AI1(i))
       STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(AI1(i)))
       CASE ('iveg'      ) ; LUC_EXPT%ivegp     = TMPI
       CASE ('biome'      ) ; LUC_EXPT%biome     = TMPI
       END SELECT
    END DO

    ! non-woody potential vegetation not considered to undergo LU change
    WHERE (LUC_EXPT%ivegp.GT.5)
       LUC_EXPT%prim_only=.TRUE.
    ENDWHERE

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)




  END SUBROUTINE READ_CLIMATEFILE

  ! ==============================================================================
  SUBROUTINE READ_LUH2(LUC_EXPT)

    IMPLICIT NONE

    TYPE (LUC_EXPT_TYPE), INTENT(INOUT) :: LUC_EXPT

    REAL    ::  tmp
    REAL, ALLOCATABLE :: tmparr(:,:)
    INTEGER :: t, i, ii, k, x, y, realk
    INTEGER :: fid, vid, tid
    INTEGER :: xds, yds, tds
    INTEGER :: STATUS,  iu

    yds = LUC_EXPT%ydimsize
    xds = LUC_EXPT%xdimsize
    t = LUC_EXPT%CTSTEP
    IF(.NOT.ALLOCATED(tmparr)) ALLOCATE(tmparr(xds,yds))

    IF (t.LE. LUC_EXPT%nrec) THEN
       DO i=1, LUC_EXPT%nfile
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
                LUC_EXPT%INPUT(i)%VAL(k) =tmparr( land_x(k), land_y(k) )
                IF  (LUC_EXPT%INPUT(i)%VAL(k).GT.1.0) THEN
                   LUC_EXPT%INPUT(i)%VAL(k) = 0.0
                ENDIF
             END DO

          ENDIF
       ENDDO

    ELSE

       WRITE(*,*) 'warning: past end of LUH2 record'

    ENDIF


    ! Adjust transition areas based on primary wooded fraction
    WHERE (LUC_EXPT%biome .EQ. 3 .OR. LUC_EXPT%biome .EQ. 11)  ! savanna/ xerophytic woods
       LUC_EXPT%INPUT(ptos)%VAL =  LUC_EXPT%INPUT(ptos)%VAL * 1.0/2.0
       LUC_EXPT%INPUT(ptog)%VAL =  LUC_EXPT%INPUT(ptog)%VAL * 1.0/2.0
       LUC_EXPT%INPUT(gtos)%VAL =  LUC_EXPT%INPUT(gtos)%VAL * 1.0/2.0
       LUC_EXPT%INPUT(stog)%VAL =  LUC_EXPT%INPUT(stog)%VAL * 1.0/2.0
    ELSEWHERE (LUC_EXPT%biome .EQ. 12 .OR. LUC_EXPT%biome .EQ. 13 &
         .OR. LUC_EXPT%biome .EQ. 15 .OR. LUC_EXPT%biome .EQ. 16  ) ! shrub
       LUC_EXPT%INPUT(ptos)%VAL =  LUC_EXPT%INPUT(ptos)%VAL * 1.0/5.0
       LUC_EXPT%INPUT(ptog)%VAL =  LUC_EXPT%INPUT(ptog)%VAL * 1.0/5.0
       LUC_EXPT%INPUT(gtos)%VAL =  LUC_EXPT%INPUT(gtos)%VAL * 1.0/5.0
       LUC_EXPT%INPUT(stog)%VAL =  LUC_EXPT%INPUT(stog)%VAL * 1.0/5.0
    ELSEWHERE (LUC_EXPT%biome .EQ. 7 .OR. LUC_EXPT%biome .EQ. 8 &
         .OR. LUC_EXPT%biome .EQ. 9 .OR. LUC_EXPT%biome .EQ. 10) ! boreal
       LUC_EXPT%INPUT(ptos)%VAL =  LUC_EXPT%INPUT(ptos)%VAL * 0.8
       LUC_EXPT%INPUT(ptog)%VAL =  LUC_EXPT%INPUT(ptog)%VAL * 0.8
       LUC_EXPT%INPUT(gtos)%VAL =  LUC_EXPT%INPUT(gtos)%VAL * 0.8
       LUC_EXPT%INPUT(stog)%VAL =  LUC_EXPT%INPUT(stog)%VAL * 0.8
    ELSEWHERE (LUC_EXPT%biome .EQ. 5 .OR. LUC_EXPT%biome .EQ. 6 ) ! DBL
       LUC_EXPT%INPUT(ptos)%VAL =  LUC_EXPT%INPUT(ptos)%VAL * 0.7
       LUC_EXPT%INPUT(ptog)%VAL =  LUC_EXPT%INPUT(ptog)%VAL * 0.7
       LUC_EXPT%INPUT(gtos)%VAL =  LUC_EXPT%INPUT(gtos)%VAL * 0.7
       LUC_EXPT%INPUT(stog)%VAL =  LUC_EXPT%INPUT(stog)%VAL * 0.7
    ENDWHERE

  END SUBROUTINE READ_LUH2
  ! ==============================================================================


END MODULE CABLE_LUC_EXPT
