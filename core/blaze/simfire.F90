MODULE SIMFIRE_MOD

TYPE TYPE_SIMFIRE
   INTEGER, DIMENSION(:), ALLOCATABLE    :: IGBP, BIOME, REGION, NDAY
   REAL,    DIMENSION(:), ALLOCATABLE    :: POPD, MAX_NESTEROV, CNEST, LAT, LON, FLI, FAPAR, POPDENS
   REAL,    DIMENSION(:,:), ALLOCATABLE  :: SAV_NESTEROV, SAV_FAPAR, BA_MONTHLY_CLIM
   INTEGER   :: SYEAR, EYEAR, NCELLS
   REAL      :: RES, RESF
   CHARACTER :: IGBPFILE*120, HYDEPATH*100, OUTMODE*6, BA_CLIM_FILE*100
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
  CHARACTER(len=400)   :: HydePath,  BurnedAreaSource, BurnedAreaFile, &
       BurnedAreaClimatologyFile, SIMFIRE_REGION
  INTEGER :: F_ID, V_ID, V_ID_lat, V_ID_lon, ilat,ilon
  INTEGER :: iu
  INTEGER :: i
  REAL, DIMENSION(720):: lon_BA
  REAL, DIMENSION(360):: lat_BA
  integer :: status
 
  NAMELIST /BLAZENML/ HydePath,  BurnedAreaSource, BurnedAreaFile, BurnedAreaClimatologyFile, &
       SIMFIRE_REGION
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! SF%RES    = 1./12.
 
 ! SF%RES = 0.5
  SF%NCELLS = NCELLS

  ALLOCATE( SF%IGBP        (NCELLS) )
  SF%IGBP = modis_igbp
  ALLOCATE( SF%BIOME       (NCELLS) )
  ALLOCATE( SF%REGION      (NCELLS) )
  ALLOCATE( SF%POPD        (NCELLS) )
  ALLOCATE( SF%MAX_NESTEROV(NCELLS) )
  ALLOCATE( SF%CNEST       (NCELLS) )
  ALLOCATE( SF%NDAY        (NCELLS) )
  ALLOCATE( SF%FAPAR       (NCELLS) )
  ALLOCATE( SF%LAT         (NCELLS) )
  SF%LAT  = LATITUDE
  ALLOCATE( SF%LON         (NCELLS) )
  SF%LON  = LONGITUDE
  ALLOCATE( SF%SAV_NESTEROV(NCELLS,12) )
  ALLOCATE( SF%SAV_FAPAR(NCELLS,FAPAR_AVG_INT) )
  ALLOCATE( SF%POPDENS     (NCELLS) )
  ALLOCATE( SF%BA_MONTHLY_CLIM(NCELLS,12)) ! fraction of annual burned area in each month
  

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
   OPEN (iu,FILE="BLAZE.nml",STATUS='OLD',ACTION='READ')
   READ (iu,NML=BLAZENML)
   CLOSE(iu)

   SF%HYDEPATH = TRIM(HydePath)
   SF%BA_CLIM_FILE = TRIM(BurnedAreaClimatologyFile)

  SF%IGBP = modis_igbp

  DO i = 1, NCELLS
     IF ( SF%IGBP(i) .LT. 1 .OR. SF%IGBP(i) .GT. 16 ) THEN
        WRITE(*,*) "Pixel i:",i," doesn't have proper IGBP veg:",SF%IGBP(i)
        SF%BIOME(i) = 0
     ELSEIF ( ABS(SF%LAT(i)) .GE. 50. ) THEN
        SF%BIOME(i) = IGBP2BIOME(SF%IGBP(i),2)
     ELSE
        SF%BIOME(i) = IGBP2BIOME(SF%IGBP(i),1)
     ENDIF

     WRITE(*,FMT='(A5,I2,A17,I1,(1X,F7.2))')"IGBP ",SF%IGBP(i)," -> SIMFIREBIOME ", &
          SF%BIOME(i),SF%LAT(i), SF%LON(i) 
     
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

  WRITE(*,*)"SIMFIRE Optimisation chosen:", TRIM(SIMFIRE_REGION)

  WRITE(*,*) "reading monthly burned area fraction from: ", TRIM(SF%BA_CLIM_FILE)
  
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
  

  
END SUBROUTINE INI_SIMFIRE

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
 SF%RES    = HYRES
  IF ( CALL1 ) THEN
     RF = NINT(SF%RES/HYRES)
     ! Check for Res being an integral multiple of 5' [RES] = fract. deg
     IF ( REAL(RF) .NE. SF%RES/HYRES .OR. SF%RES .LT. HYRES ) THEN
        WRITE(*,*) 'Spatial resolution must be integral multiple of HYDE res. '
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

     ALLOCATE( SPOPD(SF%NCELLS), EPOPD(SF%NCELLS) )
     ALLOCATE( x(SF%NCELLS),y(SF%NCELLS) )
     DO i = 1, SF%NCELLS
        X(i) = INT( (SF%LON(i) + 180.) / SF%RES ) + 1
        Y(i) = INT( (SF%LAT(i) +  90.) / SF%RES ) + 1
     END DO

     write(*,*) 'X,Y:', X, Y,SF%RES
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

  SF%RES = 1./12.
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
        write(*,*) NLAT, NLON
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
                 wPOPD = wPOPD + RVAL(ix,jy) * LAND_AREA(ix,jy)
                 wTOT  = wTOT  + LAND_AREA(ix,jy)

                 write(*,*) 'RVAL: ',   RVAL(ix,jy), ix, jy
              END DO
           END DO

           IF ( wTOT .EQ. 0. ) THEN
              WRITE(*,*)"Pixel LAT:",SF%LAT(i)," LON:",SF%LON(i)," does not contain land!"
              STOP "GET_POPDENS in simfire_mod.f90"
           ELSEIF ( nr .EQ. 1 ) THEN 
              EPOPD(i) = wPOPD / wTOT
           ELSE
              SPOPD(i) = wPOPD / wTOT
           END IF
           
        END DO
     END DO
  END IF

  ! Finally get time-interpolated population density
  tf = REAL( SF%EYEAR - YEAR ) / REAL(ISTEP)
  
  SF%POPD = tf * SPOPD + (1.-tf) * EPOPD

  !write(*,*) 'POPD', SPOPD, EPOPD,  SF%POPD , SF%LON, SF%LAT

  !stop

END SUBROUTINE GET_POPDENS

FUNCTION ANNUAL_BA ( FAPAR, FIRE_IDX, POPDENS, BIOME, REGIO_FLAG )

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
REAL                :: ANNUAL_BA

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
ELSE
   ANNUAL_BA = &
        a(BIOME,ai) * FAPAR ** b(ai) * (scalar * FIRE_IDX) ** c(ai) * EXP(e(ai)*POPDENS)


    ANNUAL_BA = &
        a(BIOME,ai) * FAPAR ** b(ai) * (scalar * FIRE_IDX) ** c(ai) * EXP(e(ai)*POPDENS) 
!CLNELSE
!CLN ! W.KNORR: Instead of fpar_corr1 * fpar_leafon + fpar_corr2 * fpar_leafon * fpar_leafon, 
!CLN ! simply use FAPAR - the correction takes into account that fpar_leafon has a high bias 
!CLN ! compared to observed FAPAR.
!CLN   ABA = a(BIOME,ai) * &
!CLN        (fpar_corr1 * FPAR_LEAFON + fpar_corr2 * FPAR_LEAFON * FPAR_LEAFON) ** b(ai) * &
!CLN        (scalar * FIRE_IDX) ** c(ai) * &
!CLN        EXP(e(ai)*POPDENS) 
ENDIF

END FUNCTION ANNUAL_BA

SUBROUTINE UPDATE_FIRE_BIOME








  
END SUBROUTINE UPDATE_FIRE_BIOME

SUBROUTINE SIMFIRE ( SF, RAINF, TMAX, TMIN, DOY,MM, YEAR, AB, climate )

  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR
  USE cable_IO_vars_module, ONLY:  landpt
  
  USE CABLE_DEF_TYPES_MOD, ONLY:  climate_type
  USE netcdf

  IMPLICIT NONE

  TYPE (TYPE_SIMFIRE) :: SF
  REAL,    INTENT(IN) :: RAINF(*), TMAX(*), TMIN(*)
  REAL,    INTENT(OUT):: AB(*)
  INTEGER, INTENT(IN) :: YEAR, MM
  TYPE (CLIMATE_TYPE), INTENT(IN)     :: climate
  
  INTEGER :: i, DOM(12), DOY

  DOM = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  IF ( IS_LEAPYEAR(YEAR) ) DOM(2) = 29
  
  SF%FAPAR = climate%AvgAnnMaxFAPAR(landpt(:)%cstart)
  SF%MAX_NESTEROV =  climate%Nesterov_ann_running_max(landpt(:)%cstart)
  !SF%BA_CLIM_FILE = "/OSM/CBR/OA_GLOBALCABLE/work/Data_BLAZE/simfire_monthly_ba.nc"
  
  ! Housekeeping first
  IF ( DOY.EQ. 1 ) THEN
     CALL GET_POPDENS ( SF, YEAR )
  ENDIF

  DO i = 1, SF%NCELLS
     AB(i) = ANNUAL_BA( SF%FAPAR(i), SF%MAX_NESTEROV(i), SF%POPD(i), SF%BIOME(i), SF%REGION(i) )
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

      if (year==2003) then
         AB(i) = 0.8        ! test vh
      endif
      
      ! Daily Burned Area
      AB(i) = AB(i) *  SF%BA_MONTHLY_CLIM(i,MM) / DOM(MM)

   END DO



END SUBROUTINE SIMFIRE

END MODULE SIMFIRE_MOD
