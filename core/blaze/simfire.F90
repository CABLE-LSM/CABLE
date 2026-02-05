MODULE SIMFIRE_MOD

TYPE TYPE_SIMFIRE
   INTEGER, DIMENSION(:), ALLOCATABLE    :: IGBP, BIOME, REGION, NDAY
   REAL,    DIMENSION(:), ALLOCATABLE    :: POPD, MAX_NESTEROV, CNEST, LAT, LON, FLI, FAPAR, POPDENS, AREA
   REAL,    DIMENSION(:,:), ALLOCATABLE  :: SAV_NESTEROV, SAV_FAPAR, BA_MONTHLY_CLIM
   INTEGER   :: SYEAR, EYEAR, NCELLS
   REAL      :: RES, RESF
   CHARACTER :: IGBPFILE*120, HYDEPATH*100, OUTMODE*6, BA_CLIM_FILE*100
   LOGICAL   :: STOCH_AREA
END TYPE TYPE_SIMFIRE

! IGBP2BIOME MAPPING:
! vegetation formation (should be derived from inputs to function)
!  1 Cropland/Urban/Natural Vegetation Mosaic (IGBP 12-14)
!  2 Needleleaf forest (IGBP 1,3): >60% cover, height>2m
!  3 Broadleaf forest (IGBP 2,4): >60% cover, height>2m
!  4 Mixed forest (IGBP 4): >60% cover, height>2m, none >60%
!  5 Shrubland (IGBP 6,7 and latitude<50): >10% woody cover, height<2m
!  6 Savanna or Grassland (IGBP 8-10): herbaceous component present, <60% tree cover
!  7 Tundra (IGBP 6,7,16 and latitude>=50): height<2m
!  8 Barren or Sparsely Vegetated (IGBP 16 and latitude<50): <10% vegetation cover
!
! IGBP:
!  0 Water bodies
!  1 Evergreen Needleleaf Forest % 1: >60% cover, height>2m
!  2 Evergreen Broadleaf Forest % 2: >60% cover, height>2m
!  3 Deciduous Needleleaf Forest % 3: >60% cover, height>2m
!  4 Deciduous Broadleaf Forest % 4: >60% cover, height>2m
!  5 Mixed Forest % 5: >60% cover, height>2m, no forest type>60% cover
!  6 Closed Shrubland % 6: >60% woody cover, height<2m
!  7 Open Shrubland % 7: 10-60% woody cover, height<2m
!  8 Woody Savanna % 8: 30-60% tree cover, height>2m, herbaceous or other understory
!  9 Savanna % 9: 10-30% tree cover, height>2m, herbaceous or other understory
! 10 Grassland % 10: <10% tree and shrub cover
! 11 Permanent Wetland % 11: mixture of water and herbaceous or woody vegetation
! 12 Cropland % 12
! 13 Urban and Built-Up % 13
! 14 Cropland/Natural Vegetation Mosaic % 14 none of forest, shrubland, cropland,
!    grassland >60% cover
! 15 Permanent Snow and Ice % 15
! 16 Barren or Sparsely Vegetated % 16: <10% vegetation cover all year
! 17 Unclassified / No data
!
INTEGER, DIMENSION(16,2), PARAMETER :: IGBP2BIOME = reshape( &
     (/ 2, 3, 2, 3, 4, 5, 5, 6, 6, 6, 0, 1, 1, 1, 0, 8 ,   & ! LAT  < 50
        2, 3, 2, 3, 4, 7, 7, 6, 6, 6, 0, 1, 1, 1, 0, 7 /), & ! LAT >= 50
     (/16,2/) )

INTEGER, PARAMETER :: FAPAR_AVG_INT = 3

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INI_SIMFIRE( NCELLS, SF, modis_igbp )

  USE CABLE_COMMON_MODULE,  ONLY: GET_UNIT, HANDLE_ERR
  USE CABLE_IO_VARS_MODULE, ONLY: LATITUDE, LONGITUDE
  use netcdf

  IMPLICIT NONE

  TYPE (TYPE_SIMFIRE), INTENT(INOUT) :: SF
  INTEGER,             INTENT(IN)    :: NCELLS, modis_igbp(NCELLS)
  CHARACTER(len=400)  :: HydePath, BurnedAreaFile = "", BurnedAreaClimatologyFile, SIMFIRE_REGION
  CHARACTER(len=10)   :: BurnedAreaSource = "SIMFIRE", blazeTStep = "annually"
  INTEGER :: F_ID, V_ID, V_ID_lat, V_ID_lon, ilat,ilon
  INTEGER :: iu
  INTEGER :: i
  REAL, DIMENSION(720):: lon_BA
  REAL, DIMENSION(360):: lat_BA
  integer :: status
  LOGICAL :: STOCHFLAG = .false.

  NAMELIST /SIMFIRENML/ SIMFIRE_REGION, HydePath, BurnedAreaClimatologyFile, STOCHFLAG

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SF%RES    = 0.5 ! CLN should stay 0.5 deg becasue SIMFIRE was trained on it. 
  SF%NCELLS = NCELLS

  ALLOCATE( SF%IGBP        (NCELLS) )
  ALLOCATE( SF%BIOME       (NCELLS) )
  ALLOCATE( SF%REGION      (NCELLS) )
  ALLOCATE( SF%POPD        (NCELLS) )
  ALLOCATE( SF%MAX_NESTEROV(NCELLS) )
  ALLOCATE( SF%CNEST       (NCELLS) )
  ALLOCATE( SF%NDAY        (NCELLS) )
  ALLOCATE( SF%FAPAR       (NCELLS) )
  ALLOCATE( SF%LAT         (NCELLS) )
  ALLOCATE( SF%LON         (NCELLS) )
  ALLOCATE( SF%AREA        (NCELLS) )
  ALLOCATE( SF%SAV_NESTEROV(NCELLS,12) )
  ALLOCATE( SF%SAV_FAPAR(NCELLS,FAPAR_AVG_INT) )
  ALLOCATE( SF%POPDENS     (NCELLS) )
  ALLOCATE( SF%BA_MONTHLY_CLIM(NCELLS,12)) ! fraction of annual burned area in each month

  call zero_simfire(SF)

  SF%IGBP = modis_igbp
  SF%LAT  = LATITUDE
  SF%LON  = LONGITUDE

  !get effective grid cell area - dummy until robust solution determined
  CALL get_grid_areakm2(NCELLS,SF%LAT,SF%LON,SF%AREA)

  !=============================================================================
  ! VEGTYPE from IGBP dataset
  !=============================================================================

  ! Assign SIMFIRE Parameters
  ! IGBP, SF-BIOME and REGION(of optimization)

  !CLN READ modis-IGBP here!!!
  !vh!
  ! inherit modis_igbp from climate variable

  ! READ BLAZE settings
  CALL GET_UNIT(iu)
  OPEN (iu,FILE="blaze.nml",STATUS='OLD',ACTION='READ')
  READ (iu,NML=SIMFIRENML)
  CLOSE(iu)
  
  SF%HYDEPATH = TRIM(HydePath)
  SF%BA_CLIM_FILE = TRIM(BurnedAreaClimatologyFile)
  !WRITE(*,*)"SIMFIRENML :", SIMFIRE_REGION, HydePath, BurnedAreaClimatologyFile
  SF%IGBP = modis_igbp
  SF%STOCH_AREA = STOCHFLAG

  DO i = 1, NCELLS
     IF ( SF%IGBP(i) .LT. 1 .OR. SF%IGBP(i) .GT. 16 ) THEN
        WRITE(*,*) "Pixel i:",i," doesn't have proper IGBP veg:",SF%IGBP(i)
        SF%BIOME(i) = 0
     ELSEIF ( ABS(SF%LAT(i)) .GE. 50. ) THEN
        SF%BIOME(i) = IGBP2BIOME(SF%IGBP(i),2)
     ELSE
        SF%BIOME(i) = IGBP2BIOME(SF%IGBP(i),1)
     ENDIF
     
     !WRITE(*,FMT='(A5,I2,A17,I1,(1X,F7.2))')"IGBP ",SF%IGBP(i)," -> SIMFIREBIOME ", &
     !     SF%BIOME(i),SF%LAT(i), SF%LON(i)

     IF ( TRIM(SIMFIRE_REGION) == "GLOBAL" ) THEN ! GLOBAL
        SF%REGION(i) = 1
     ELSE IF ( TRIM(SIMFIRE_REGION) == "ANZ" ) THEN ! AUSTRALIA
        SF%REGION(i) = 2
     ELSE IF ( TRIM(SIMFIRE_REGION) == "EUROPE" ) THEN ! EUROPE
        SF%REGION(i) = 3
     ELSE
        WRITE(*,*)"Invalid SIMFIRE_REGION in cable.nml:",SIMFIRE_REGION
        STOP
     ENDIF
  END DO

  !WRITE(*,*)"SIMFIRE Optimisation chosen:", TRIM(SIMFIRE_REGION)

  !WRITE(*,*) "reading monthly burned area fraction from: ", TRIM(SF%BA_CLIM_FILE)

  STATUS = NF90_OPEN(TRIM(SF%BA_CLIM_FILE), NF90_NOWRITE, F_ID)
  CALL HANDLE_ERR(STATUS, "Opening BA Clim File "//SF%BA_CLIM_FILE )
  STATUS = NF90_INQ_VARID(F_ID,'monthly_ba', V_ID)
  CALL HANDLE_ERR(STATUS, "Inquiring  var monthly_ba in "//SF%BA_CLIM_FILE )
  STATUS = NF90_INQ_VARID(F_ID,'latitude', V_ID_lat)
  CALL HANDLE_ERR(STATUS, "Inquiring  var latitude in "//SF%BA_CLIM_FILE )
  STATUS = NF90_INQ_VARID(F_ID,'longitude', V_ID_lon)
  CALL HANDLE_ERR(STATUS, "Inquiring  var longitude in "//SF%BA_CLIM_FILE )
  STATUS = NF90_GET_VAR( F_ID, V_ID_lat, lat_BA, &
        start=(/1/)  )
  STATUS = NF90_GET_VAR( F_ID, V_ID_lon, lon_BA, &
           start=(/1/)  )

  DO i = 1, SF%NCELLS

      ilat = MINLOC(ABS(lat_BA - SF%LAT(i)),DIM=1)
      ilon = MINLOC(ABS(lon_BA - SF%LON(i)),DIM=1)

      STATUS = NF90_GET_VAR( F_ID, V_ID, SF%BA_MONTHLY_CLIM(i,:), &
           start=(/1,ilon,ilat/) )
      CALL HANDLE_ERR(STATUS, "Reading direct from "//SF%BA_CLIM_FILE )

  ENDDO

  STATUS = NF90_CLOSE(F_ID)

  WRITE(*,*) "application of stochastic burnt area is", SF%STOCH_AREA

END SUBROUTINE INI_SIMFIRE


SUBROUTINE zero_simfire(SF)

  implicit none

  type(type_simfire), intent(inout) :: sf

  sf%IGBP         = 0
  sf%BIOME        = 0
  sf%REGION       = 0
  sf%POPD         = 0
  sf%MAX_NESTEROV = 0
  sf%CNEST        = 0
  sf%NDAY         = 0
  sf%FAPAR        = 0
  sf%LAT          = 0
  sf%LON          = 0
  ! sf%FLI          = 0
  sf%SAV_NESTEROV = 0
  sf%SAV_FAPAR    = 0
  sf%POPDENS      = 0
  sf%BA_MONTHLY_CLIM = 0

END SUBROUTINE zero_simfire


SUBROUTINE print_simfire(SF)

  implicit none

  type(type_simfire), intent(in) :: sf

  write(*,*) 'IGBP ', sf%IGBP
  write(*,*) 'BIOME ', sf%BIOME
  write(*,*) 'REGION ', sf%REGION
  write(*,*) 'POPD ', sf%POPD
  write(*,*) 'MAX_NESTEROV ', sf%MAX_NESTEROV
  write(*,*) 'CNEST ', sf%CNEST
  write(*,*) 'NDAY ', sf%NDAY
  write(*,*) 'FAPAR ', sf%FAPAR
  write(*,*) 'LAT ', sf%LAT
  write(*,*) 'LON ', sf%LON
  write(*,*) 'SAV_NESTEROV ', sf%SAV_NESTEROV
  write(*,*) 'SAV_FAPAR ', sf%SAV_FAPAR
  write(*,*) 'POPDENS ', sf%POPDENS
  write(*,*) 'BA_MONTHLY_CLIM ', sf%BA_MONTHLY_CLIM

END SUBROUTINE print_simfire


SUBROUTINE GET_POPDENS ( SF, YEAR )

  USE CABLE_COMMON_MODULE, ONLY: GET_UNIT

  IMPLICIT NONE

  TYPE(TYPE_SIMFIRE)  :: SF
  INTEGER, INTENT(IN) :: YEAR

  REAL,    PARAMETER :: HYRES = 1./12. ! [deg]
  INTEGER, PARAMETER :: NLON=4320, NLAT=2160
  INTEGER, PARAMETER :: FINAL_YEAR = 2005

  REAL, DIMENSION(NLON,NLAT) :: RVAL
  INTEGER       :: ISTEP, iu, i ,j, nr, ix0, jy0, ix, jy, dxy, NREAD
  INTEGER       :: RYEAR
  REAL          :: wPOPD, wTOT, tf
  CHARACTER     :: FNAME*100, cYEAR*5, suf*2
  LOGICAL       :: EXISTFILE

  INTEGER, SAVE :: RF
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: X,Y
  REAL,    DIMENSION(:), ALLOCATABLE, SAVE :: SPOPD, EPOPD
  REAL,    SAVE :: LAND_AREA(NLON,NLAT)
  LOGICAL, SAVE :: CALL1 = .TRUE.

  ! CLN Put into cable_in
  !SF%HYDEPATH = "/OSM/CBR/OA_GLOBALCABLE/work/Data_BLAZE/HYDE3.1"

  !=============================================================================
  ! POPDENS Population Density from HYDE 3.1 popd comes at 5' res
  !=============================================================================

  ! Determine current interpolation step
  IF ( YEAR .LT. 0 ) THEN
     ISTEP = 1000
  ELSE IF ( YEAR .LT. 1700 ) THEN
     ISTEP = 100
  ELSE IF ( YEAR .LT. 2000 ) THEN
     ISTEP = 10
  ELSE 
     ISTEP = 5
  END IF
  !CLN WROOONGGG!!!! SF%RES    = HYRES

  !INH need to deal with runs that start after FINAL_YEAR here
  IF ( CALL1 ) THEN
     RF = NINT(SF%RES/HYRES)
     ! Check for Res being an integral multiple of 5' [RES] = fract. deg
     IF ( REAL(RF) .NE. SF%RES/HYRES .OR. SF%RES .LT. HYRES ) THEN
        WRITE(*,*) 'Spatial resolution must be integer multiple of HYDE res. '
        WRITE(*,*) "RES:",SF%RES,"/ HYDE:",HYRES," = ",SF%RES/HYRES
        STOP "get_popdens in simfire_mod.f90"
     END IF

     ! READ Land_area for each gridcell [km^2]
     CALL GET_UNIT(iu)
     FNAME = TRIM(SF%HYDEPATH)//"/mland_cr.asc"
     INQUIRE( FILE=TRIM(FNAME), EXIST=EXISTFILE)
     IF ( .NOT. EXISTFILE ) THEN
        WRITE(*,*)"Hyde 3.1 Pop Dens mland_cr File doesnt exist!",TRIM(FNAME)
        STOP -1
     END IF

     OPEN(iu, FILE=TRIM(FNAME), ACTION="READ", STATUS="OLD")
     ! Skip header
     DO i = 1, 6
        READ(iu,*)
     END DO
     ! Read data backwards ( llcorner is -180E, -90N )
     DO j = NLAT, 1, -1
        READ(iu,*)(LAND_AREA(i,j),i=1,NLON)
     END DO
     CLOSE(iu)

     ! Get years to interpolate between
     IF ( YEAR .LT. 0 ) THEN
        SF%SYEAR = INT(REAL(YEAR)/REAL(ISTEP) - 1. ) * ISTEP
     ELSE
        SF%SYEAR = INT(REAL(YEAR)/REAL(ISTEP)) * ISTEP
     ENDIF
     SF%EYEAR = SF%SYEAR + ISTEP

     !runs that start after FINAL_YEAR
     IF (YEAR .GE. FINAL_YEAR) THEN
        SF%EYEAR = FINAL_YEAR
        SF%SYEAR = FINAL_YEAR
     END IF

     ALLOCATE( SPOPD(SF%NCELLS), EPOPD(SF%NCELLS) )
     ALLOCATE( x(SF%NCELLS),y(SF%NCELLS) )
     DO i = 1, SF%NCELLS
        X(i) = INT( (SF%LON(i) + 180.) / SF%RES ) + 1
        Y(i) = INT( (SF%LAT(i) +  90.) / SF%RES ) + 1
     END DO

     !write(*,*) 'X,Y:', X, Y,SF%RES
     NREAD = 2
     CALL1 = .FALSE.
  ELSE IF ( YEAR .GE. FINAL_YEAR ) THEN
     SF%POPD  = EPOPD
     WRITE(*,*)"No population development after ", FINAL_YEAR
     RETURN
  ELSE IF ( YEAR .EQ. SF%EYEAR ) THEN
     SF%SYEAR = SF%EYEAR
     SPOPD = EPOPD
     SF%EYEAR = SF%SYEAR + ISTEP
     NREAD = 1
  ELSE
     NREAD = 0
  END IF

  !CLN SF%RES = 1./12.
  RF = NINT(SF%RES/HYRES)


  IF ( NREAD .GT. 0 ) THEN
     CALL GET_UNIT(iu)
     DO nr = 1, NREAD
        IF ( nr .EQ. 1 ) THEN
           RYEAR = SF%EYEAR
        ELSE
           RYEAR = SF%SYEAR
        ENDIF

        IF ( ABS(RYEAR) .GE. 10000 ) THEN
           WRITE(cYEAR,FMT="(I5)")ABS(RYEAR)
        ELSE IF ( ABS(RYEAR) .GE. 1000 ) THEN
           WRITE(cYEAR,FMT="(I4,A1)")ABS(RYEAR)," "
        ELSE IF ( ABS(RYEAR) .GE. 100  ) THEN
           WRITE(cYEAR,FMT="(I3,A2)")ABS(RYEAR),"  "
        ELSE
           cYEAR = "0    "
        END IF

        IF ( RYEAR .LT. 0 ) THEN
           suf = "BC"
        ELSE
           suf = "AD"
        ENDIF
        FNAME = TRIM(SF%HYDEPATH)//"/data/popc_"//TRIM(cYEAR)//suf//".asc"
        INQUIRE( FILE=TRIM(FNAME), EXIST=EXISTFILE)
        IF ( .NOT. EXISTFILE ) THEN
           WRITE(*,*)"Hyde 3.1 Pop Dens File doesnt exist!",TRIM(FNAME)
           STOP -1
        ENDIF

        OPEN(iu, FILE=TRIM(FNAME), ACTION="READ", STATUS="OLD")
        ! Skip header
        DO i = 1, 6
           READ(iu,*)
        END DO
        ! Read data backwards ( llcorner is -180E, 90N )
        !write(*,*) NLAT, NLON
        DO j = NLAT, 1, -1
           READ(iu,*)(RVAL(i,j),i=1,NLON)
          ! write(3334,"(4320e16.6)") (RVAL(i,j),i=1,NLON)
        END DO
        CLOSE(iu)


        DO i = 1, SF%NCELLS

           ix0 = RF * (X(i)-1) + 1
           jy0 = RF * (Y(i)-1) + 1
           dxy = RF - 1
           !write(*,*) 'ix0,iy0', RF,ix0, jy0
           ! average over sub-gridcells, weighted by land area of cell
           wPOPD = 0.
           wTOT  = 0.
           DO jy = jy0, jy0+dxy
              DO ix = ix0, ix0+dxy
                 IF ( RVAL(ix,jy) .LT. 0. ) CYCLE
                 !9/2025 popc data is humans per HYDE grid cell 
                 !- need population density humans / km2 
                 wPOPD = wPOPD + RVAL(ix,jy) !* LAND_AREA(ix,jy)
                 wTOT  = wTOT  + LAND_AREA(ix,jy)

                 !write(*,*) 'RVAL: ',   RVAL(ix,jy), ix, jy

              END DO
           END DO

           IF ( wTOT .EQ. 0. ) THEN
              ! There are land-pixels in the CABLE grid that are non-land in HYDE 3.1
              ! therefore assume zero population on these
              WRITE(*,*)"Pixel LAT:",SF%LAT(i)," LON:",SF%LON(i)," does not contain land in HYDE population data"
              !CLN STOP "GET_POPDENS in simfire_mod.f90"
              wPOPD = 0.
              wTOT  = 1.
           ENDIF

           IF ( nr .EQ. 1 ) THEN
              EPOPD(i) = wPOPD / wTOT
           ELSE
              SPOPD(i) = wPOPD / wTOT
           END IF

        END DO
     END DO
  END IF

  ! Finally get time-interpolated population density
  tf = REAL( SF%EYEAR - YEAR ) / REAL(ISTEP)
  
  !INH for runs that start after FINAL_YEAR - SPOPD=EPOPD so tf irrelevant
  SF%POPD = tf * SPOPD + (1.-tf) * EPOPD

  !write(*,*) 'POPD', SPOPD, EPOPD,  SF%POPD , SF%LON, SF%LAT

  !stop

END SUBROUTINE GET_POPDENS

!FUNCTION ANNUAL_BA ( FAPAR, FIRE_IDX, POPDENS, BIOME, REGIO_FLAG )
SUBROUTINE GET_ANNUAL_BA ( FAPAR, FIRE_IDX, POPDENS, BIOME, REGIO_FLAG, ANNUAL_BA, T1, T2, T3, T4 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Annual burned area model by W.Knorr et al. 2014
! Ref.:  http://www.biogeosciences.net/11/1085/2014/bg-11-1085-2014.html
!
! Implementation by L.P. Nieradzik 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE

REAL   , INTENT(IN) :: FAPAR, FIRE_IDX, POPDENS
INTEGER, INTENT(IN) :: BIOME, REGIO_FLAG
!REAL                :: ANNUAL_BA
REAL, INTENT(OUT)   :: ANNUAL_BA, T1, T2, T3, T4

INTEGER :: ai

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PARAMETERS FOR SIMFIRE:
! GLOBAL (1) , AUSTRALIA-NZ-OPTIMISED (2), AND EUROPE (3)
! BIOME COEFFICIENTS a(BIOME)
REAL, DIMENSION(8,3), PARAMETER :: a = reshape( &
     (/0.110,  0.095    ,0.092  ,0.127  ,0.470  ,0.889 ,0.059  ,0.113,     & ! GLOBAL
       0.06974,0.6535   ,0.6341 ,0.6438 ,2.209  ,1.710 ,    0. ,2.572,     & ! ANZ
       0.02589,0.0008087,0.04896,0.06248,0.01966,0.1191,0.01872,0.08873/), & ! EUR
       (/8,3/) )
! Biome:  crop NLfor     BLfor   mixedfor shrub  grass  tundra  barren
!
! fAPAR EXPONENT b
REAL, DIMENSION(3), PARAMETER :: b = (/0.905, & ! GLOBAL
                                       1.297, & ! ANZ
                                       0.9164/) ! EUR

! NESTEROV-INDEX EXPONENT c
REAL, DIMENSION(3), PARAMETER :: c = (/0.860, & ! GLOBAL
                                       1.038, & ! ANZ
                                       0.4876/) ! EUR

! POPULATION COEFFICIENT e
!CLNORIREAL, DIMENSION(3), PARAMETER :: e = (/-0.0168, & ! GLOBAL
!CLNORI                                       -0.2131, & ! ANZ
!CLNORI                                       -0.017 /)  ! EUR
REAL, DIMENSION(3), PARAMETER :: e = (/-0.0168, & ! GLOBAL
                                       -0.05  , & ! ANZ
                                       -0.017 /)  ! EUR

! CORRECTION FACTORS (only for fpar_leafon!)
! REAL, PARAMETER :: fpar_corr1 = 0.428
! REAL, PARAMETER :: fpar_corr2 = 0.148
REAL, PARAMETER :: scalar     = 1e-5

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IF ( REGIO_FLAG .LT. 1 .OR. REGIO_FLAG .GT. 3 )THEN
   WRITE(*,*)"Wrong REGIO_FLAG chosen for SIMFIRE: ",REGIO_FLAG
   STOP "Either global=1, AUS/NZ=2, EUR=3!"
ENDIF

IF ( BIOME .LT. 0 .OR. BIOME .GT. 8 ) THEN
   WRITE(*,*)"Wrong BIOME chosen for SIMFIRE: ",BIOME
   STOP "Should be element of [0,8]! Check simfire_mod.f90!"
ENDIF

IF ( FAPAR .LT. 0. .OR. FIRE_IDX .LT. 0. ) THEN
   WRITE(*,*)"Strange parmeters input to SIMFIRE:"
   WRITE(*,*)"FAPAR: ",FAPAR," FIRE_IDX: ",FIRE_IDX
   STOP "Check simfire.f90!"
ENDIF

ai    = REGIO_FLAG
IF ( BIOME .EQ. 0 ) THEN
   ANNUAL_BA = 0.
   T1 = 0.
   T2 = 0.
   T3 = 0.
   T4 = 0.
ELSE
   ANNUAL_BA = &
        a(BIOME,ai) * FAPAR ** b(ai) * (scalar * FIRE_IDX) ** c(ai) * EXP(e(ai)*POPDENS)
   T1 = a(BIOME,ai)
   T2 = FAPAR ** b(ai) 
   T3 = (scalar * FIRE_IDX) ** c(ai)
   T4 = EXP(e(ai)*POPDENS)
!CLNELSE
!CLN ! W.KNORR: Instead of fpar_corr1 * fpar_leafon + fpar_corr2 * fpar_leafon * fpar_leafon,
!CLN ! simply use FAPAR - the correction takes into account that fpar_leafon has a high bias
!CLN ! compared to observed FAPAR.
!CLN   ABA = a(BIOME,ai) * &
!CLN        (fpar_corr1 * FPAR_LEAFON + fpar_corr2 * FPAR_LEAFON * FPAR_LEAFON) ** b(ai) * &
!CLN        (scalar * FIRE_IDX) ** c(ai) * &
!CLN        EXP(e(ai)*POPDENS)
ENDIF

!END FUNCTION ANNUAL_BA
END SUBROUTINE GET_ANNUAL_BA

SUBROUTINE SIMFIRE ( SF, RAINF, TMAX, TMIN, DOY,MM, YEAR, AB, annAB, climate, FAPARSOURCE, FSTEP, T1,T2,T3,T4 )

  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR
  USE cable_IO_vars_module, ONLY:  landpt, patch

  USE CABLE_DEF_TYPES_MOD, ONLY:  climate_type
  USE netcdf

  IMPLICIT NONE

  TYPE (TYPE_SIMFIRE) :: SF
  REAL,    INTENT(IN) :: RAINF(*), TMAX(*), TMIN(*)
  REAL,    INTENT(OUT):: AB(*)
  REAL,    INTENT(INOUT) :: annAB(*)
  REAL,    INTENT(OUT):: T1(*), T2(*), T3(*), T4(*)
  INTEGER, INTENT(IN) :: YEAR, MM
  CHARACTER(len=10), INTENT(IN) :: FAPARSOURCE
  CHARACTER(len=7), INTENT(IN)  :: FSTEP                !trigger on whether to use daily/annual Nesterov
  TYPE (CLIMATE_TYPE), INTENT(IN)     :: climate
  INTEGER, PARAMETER :: stoch_trig = 10                 !number of days between call to stoch generator

  INTEGER :: i, DOM(12), DOY, p, patch_index, iSTOCH

  DOM = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  IF ( IS_LEAPYEAR(YEAR) ) DOM(2) = 29

  !INH 2025-04 - add case whereby inline evaluation of %fapar_ann_max_last_year used
  ! as the climatolgical value for AvgAnnMaxFAPAR.  Choice controlled by switch %faparsource in BLAZE NML
  ! AvgAnnMaxFAPAR is read in from file %faparfilename and constant through run
  IF (FAPARSOURCE .eq. "inline") THEN
      !evaluate (current estimate) grid-cell averaged max_fapar
      DO i=1,SF%NCELLS
         SF%FAPAR(i) = 0.0
         DO p = 1, landpt(i)%nap  ! loop over number of active patches
            patch_index = landpt(i)%cstart + p - 1 ! patch index in CABLE vector
            SF%FAPAR(i) = SF%FAPAR(i) + real(climate%fapar_ann_max_last_year(patch_index)*patch(patch_index)%frac)
         ENDDO
         !constrain %FAPAR and apply scaling factor to account for mismatch between CABLE-POP and RSensing
         SF%FAPAR(i) = 0.5*MIN(MAX(SF%FAPAR(i),0.0),1.0)
         
      ENDDO

  ELSE  !FAPARSOURCE is "fromfile - original code
     SF%FAPAR = climate%AvgAnnMaxFAPAR(landpt(:)%cstart)
  END IF

  !INH: using BLAZE%FSTEP to switch between daily Nesterov and annual max Nesterov
  ! - defaults to annual max
  IF (TRIM(FSTEP) .eq. "daily") THEN
     !over 1950-2020 current_Nestrov is ~0.4 of annual Max Nesterov
     !100000 is max value allowed in cable_climate
     SF%MAX_NESTEROV = 2.5*climate%Nesterov_Current(landpt(:)%cstart)
     SF%MAX_NESTEROV = MIN(SF%MAX_NESTEROV,100000.0)
  ELSE
     SF%MAX_NESTEROV =  climate%Nesterov_ann_running_max(landpt(:)%cstart)
  END IF
  !SF%BA_CLIM_FILE = "/OSM/CBR/OA_GLOBALCABLE/work/Data_BLAZE/simfire_monthly_ba.nc"

  ! Housekeeping first
  IF ( DOY.EQ. 1 ) THEN
     CALL GET_POPDENS ( SF, YEAR )
  ENDIF

  DO i = 1, SF%NCELLS
     !AB(i) = ANNUAL_BA( SF%FAPAR(i), SF%MAX_NESTEROV(i), SF%POPD(i), SF%BIOME(i), SF%REGION(i) )
     CALL GET_ANNUAL_BA( SF%FAPAR(i), SF%MAX_NESTEROV(i), SF%POPD(i), SF%BIOME(i), SF%REGION(i), &
                         AB(i), T1(i), T2(i), T3(i), T4(i) )

     !catch out-of-bounds issues
     AB(i) = MAX( 0., MIN(AB(i),.99) )

     !write(*,*) 'SF%FAPAR, SF%MAX_NESTEROV, SF%POPD, SF%BIOME, SF%REGION, AB'
     !write(*,"(200e16.6)") ,SF%FAPAR(i), SF%MAX_NESTEROV(i), SF%POPD(i), real(SF%BIOME(i)), real(SF%REGION(i)), AB(i)

!!$     ! convert to daily burned area using GFED climatology
!!$
!!$      STATUS = NF90_OPEN(TRIM(SF%BA_CLIM_FILE), NF90_NOWRITE, F_ID)
!!$      CALL HANDLE_ERR(STATUS, "Opening BA Clim File "//SF%BA_CLIM_FILE )
!!$      STATUS = NF90_INQ_VARID(F_ID,'monthly_ba', V_ID)
!!$      CALL HANDLE_ERR(STATUS, "Inquiring  var monthly_ba in "//SF%BA_CLIM_FILE )
!!$
!!$      STATUS = NF90_INQ_VARID(F_ID,'latitude', V_ID_lat)
!!$      CALL HANDLE_ERR(STATUS, "Inquiring  var latitude in "//SF%BA_CLIM_FILE )
!!$
!!$      STATUS = NF90_INQ_VARID(F_ID,'longitude', V_ID_lon)
!!$      CALL HANDLE_ERR(STATUS, "Inquiring  var longitude in "//SF%BA_CLIM_FILE )
!!$
!!$      STATUS = NF90_GET_VAR( F_ID, V_ID_lat, lat_BA, &
!!$               start=(/1/) )
!!$      STATUS = NF90_GET_VAR( F_ID, V_ID_lon, lon_BA, &
!!$           start=(/1/)  )
!!$
!!$      ilat = MINLOC(ABS(lat_BA - SF%LAT(i)),DIM=1)
!!$      ilon = MINLOC(ABS(lon_BA - SF%LON(i)),DIM=1)
!!$
!!$      STATUS = NF90_GET_VAR( F_ID, V_ID, monthly_ba, &
!!$           start=(/MM,ilon,ilat/) )
!!$      CALL HANDLE_ERR(STATUS, "Reading direct from "//SF%BA_CLIM_FILE )

      !if (year==2003) then
      !   AB(i) = 0.8        ! test vh
      !endif

      !apply randomness to AB annual
      !IF (SF%STOCH_AREA) THEN
      !   CALL STOCH_AREA(AB(i),SF%LAT(i),SF%LON(i),SF%AREA(i),1,YEAR,365)
      !END IF

      ! Monthly Burned Area
      ! Both BLAZE%FSTEP options use monthly_clim information to go from annual->monthly
      AB(i) = AB(i) *  SF%BA_MONTHLY_CLIM(i,MM) 

      !apply randomness to AB monthly 
      IF (SF%STOCH_AREA) THEN
         !iSTOCH increments by 1 each stoch_trig days 
         !- if stoch_trig>1 then allow some fires to last multiple days
         !- to get daily stochaticity set stoch_trig=1, to get annual set to >366
         iSTOCH = INT( REAL(DOY)/REAL(stoch_trig)) + 1
        
         !STOCH_AREA is deterministic, if iSTOCH is the same then the same AB will be produced.
         CALL STOCH_AREA(AB(i),SF%LAT(i),SF%LON(i),SF%AREA(i),iSTOCH,YEAR,DOM(MM))

      END IF 

      !Daily burned area
      AB(i) = AB(i) / DOM(MM)

      !enforce that annAB can't exceed 1.0 through year
      !sumBLAZE accumulation (equiv to annAB) is done in update_sumblaze
      AB(i) = MAX( 0., MIN(AB(i),0.99-annAB(i)) )
      annAB(i) = annAB(i) + AB(i)

   END DO

END SUBROUTINE SIMFIRE

SUBROUTINE STOCH_AREA(AB,lat,lon,gca,iSTOCH,YEAR,NDAY)

   !simple, empirical function to apply some deterministic randomness to the burnt area 
   ! while preserving the time/ensemble average.  this is a surrogate for ignition
   ! needs to be gridcell area dependent since (conceptually) more likely to be 0 or 1 AB for small cells
   !
   !function is:  ABout = ABin ( 1 + tanh( (seed - gamma AB)/ ABx)) / denom
   !
   !interpretation - 2 parts - given a seed on 0-1
   !               a) ABout1 = mx (1 - sign( (seed - gamma ABin) ) ) / 2
   !where (1-gamma ABin)*mx = AB and mx is maxBA permitted per call to randmiser 
   !
   !This function randomises AB so that ABout=0 with a probability of gamma ABin and =mx with a probability of ABin
   !The mean value of ABout (for x in 0,1) is ABin
   !
   !               b) ABout2 = mx (1 + tanh ((seed - gamma ABin)/ABx) ) / denom
   !
   !This blurs the step function nature of ABout1 over a range of seed of scale ABx, centred on seed = gamma AB
   !denom is needed to ensure that mean value of ABout (for x in 0,1) is ABin

      IMPLICIT NONE

   !in/outs
   REAL, INTENT(INOUT) :: AB
   REAL, INTENT(IN)    :: lat, lon, gca
   INTEGER, INTENT(IN) :: iSTOCH, YEAR, NDAY

   !parameters: eps needed to ensure randmomisation occurs globally, beta>=3 scaling within function 
   !max_ba = maxBA possible per day ~200km2  
   REAL, PARAMETER     :: max_ba = 200., eps = 1.0, beta = 3.0

   !working vars
   REAL                :: mx, ABx, seed, denom, gamma, dlta
   INTEGER             :: intseed

   IF ( (AB .le. 0.0) ) THEN
      !satisfy INOUT requirement
      AB = 1.0*AB

   ELSE
      !create a 'random' integer that is deterministic but varies in time, space (and ~randomly)
      !need intseed to differ for each grid cell, year and call (as given by iSTOCH)
      ! - note expected min resolution is 0.05 degrees so YEAR factor ensures cells do not align
      intseed = INT(REAL(iSTOCH*YEAR)*ABS(lat + lon + eps)*ABS(lat - lon + eps)/2.0)

      !use alternate lcg_generator to generate REAL seed on (0,1) - doesn't require the repeated ALLOCATION
      CALL seeded_random(intseed,seed)

      !apply function using seed - mx is max_BA per call to randomiser = nday_fire*max_ba_per_day/grid_cell_area
      mx = MAX(MIN(REAL(MIN(NDAY,5)) * max_ba / gca, 1.0), AB)

      !gamma*AB is positioning of max steepness
      gamma = (1.0-AB/mx)/AB

      ! 16/1/2025 - dlta is a representative grd scale in degrees
      ! allows more variation in AB with larger grid cells, more bimodal (burn/no burn) at small cells
      ! gca in km2, beta~3 required to ensure that AB reaches max/min values for seed in (0,1)
      ! smaller minimum value allows larger ABout at small values of ABin but can break the average
      !
      ! important note: if AB < 1 - 89.5*ABx then single precision numerics fail (cosh(89.5)->infty)
      ! minimum value needs to be >~ 0.011 (resolution and stoch_trig invariant)
      dlta = SQRT(gca/100.0/100.0)
      ABx = max( min(gamma*AB/beta,(1.-gamma*AB)/beta, dlta/beta), 0.0125)

      !normalising factor to ensure mean (over seed=0,1) is AB
      denom = 1.0 + ABx*( log(cosh((gamma*AB-1.0)/ABx))-log(cosh(gamma*AB/ABx)))

      !apply randomisation
      AB = AB*(1.0 + tanh( (seed-gamma*AB)/ABx) ) / denom
   
   END IF

END SUBROUTINE STOCH_AREA

SUBROUTINE seeded_random(x,seed)
   ! Linear congruential generator for single precision random numbers on (0,1)
   ! uses the matlab incarnation of coefficients (random0, Chapman 2015)
   ! https://en.wikipedia.org/wiki/Linear_congruential_generator
   !
   ! note x has to be <= 264432 in order to avoid overflow issues
   ! behaviour is compiler & serial/mpi dependent - so see catch on seed

   real, intent(out) :: seed
   integer, intent(in) :: x
   integer :: work

   !use matlab/fortran coeffs
   work = modulo (8121*x+28411, 134456)
   seed = work / real(134456)
   IF (seed .lt. 0.0) THEN
      seed = seed + 1.0
   END IF
 
END SUBROUTINE seeded_random

SUBROUTINE get_minlatlon_increment(NCELLS,latlon,dlatlon)

   !routine determines the effective resolution of input from SF TYPE 
   !dummy routine while a more robust solution found
   INTEGER, iNTENT(IN) :: NCELLS 
   REAL, INTENT(IN)    :: latlon(NCELLS)
   REAL, INTENT(OUT)   :: dlatlon
   
   INTEGER :: i, j 

   IF ((NCELLS .eq. 1) .or. (MAXVAL(latlon) .eq. MINVAL(latlon)) ) THEN
      dlatlon = 0.05
   ELSE
      !initial value
      dlatlon = 360.0
      DO i = 1,NCELLS 
         DO j = i,NCELLS 
            !loop over cells to find smallest dlat or dlon
            IF ( (abs(latlon(i)-latlon(j)) < dlatlon) .and. (latlon(i) .ne. latlon(j)) ) THEN
               dlatlon = abs(latlon(i)-latlon(j))
            END IF
         END DO
      END DO
   END IF

   !enforce a minimum value (single site runs)
   dlatlon = MAX(dlatlon,0.05)

END SUBROUTINE get_minlatlon_increment

SUBROUTINE get_grid_areakm2(NCELLS,lat,lon,area)

   !NB assumes that lat,lon are given as grid cell centres
   INTEGER, INTENT(IN) :: NCELLS
   REAL, INTENT(IN)    :: lat(NCELLS), lon(NCELLS)
   REAL, INTENT(OUT)   :: area(NCELLS)

   REAL :: dlat, dlon, PI, rEarth=6371.0
   INTEGER :: i

   CALL get_minlatlon_increment(NCELLS,lat,dlat)
   CALL get_minlatlon_increment(NCELLS,lon,dlon)

   WRITE(*,*) "SIMFIRE: resolution of CABLE simulation, dlat: ", dlat, " dlon: ", dlon

   PI = 4.0*ATAN(1.0)
             
   piR2 = PI*rEarth*rEarth/180.0
   deg2rad = PI/180.0
   DO i=1,NCELLS
      area(i) = piR2*dlon*ABS( sin(deg2rad*(lat(i)-dlat/2.0))-sin(deg2rad*(lat(i)+dlat/2.0)))
   END DO

END SUBROUTINE get_grid_areakm2


END MODULE SIMFIRE_MOD
