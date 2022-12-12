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
! Purpose: Input and output code for CASA-CNP when run offline
!          ACCESS version may use some of this code but split into different files?
!
! Contact: Yingping.Wang@csiro.au and Bernard.Pak@csiro.au
!
! History: Developed for offline code.  Expect to re-write for MPI and ACCESS
!          versions
!
!
! ==============================================================================
! casa_inout.f90
!
! the following routines are used when "casacnp" is coupled to "cable"
!   casa_readbiome
!   casa_readphen
!   casa_readpoint   (removed, now done in parameter_module)
!   casa_init
!   casa_poolout
!   casa_cnpflux  (zeros casabal quantites on doy 1 and updates casabal at end of biogeochem)
!   biogeochem
!CABLE_LSM:This has to be commented for offline
!#define UM_BUILD YES
MODULE casa_offline_inout_module

USE casavariable, ONLY : casafile

CONTAINS

#ifndef UM_BUILD
  SUBROUTINE WRITE_CASA_RESTART_NC ( casamet, casapool, casaflux, phen, CASAONLY )

    USE casavariable, ONLY : casa_met, casa_pool, casa_flux, icycle, mplant, mlitter, msoil
    USE cable_common_module
  USE casa_ncdf_module, ONLY: HANDLE_ERR
   
    USE cable_def_types_mod, ONLY: met_type, mp
    USE phenvariable
    USE netcdf

    IMPLICIT NONE


    TYPE (casa_met),  INTENT(IN) :: casamet
    TYPE (casa_pool),  INTENT(IN) :: casapool
    TYPE (casa_flux),           INTENT(IN) :: casaflux
    TYPE (phen_variable),       INTENT(IN) :: phen

    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS
    INTEGER*4   :: FILE_ID, land_ID, plnt_ID, litt_ID, soil_ID, i
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

    ! ! 1 dim arrays (npt )
    ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
    !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
    ! ! 2 dim arrays (npt,mplant)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
    ! ! 2 dim arrays (npt,mlitter)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
    ! ! 2 dim arrays (npt,msoil)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)

    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(12) :: A1
    CHARACTER(len=20),DIMENSION(2) :: AI1
    ! 2 dim arrays (npt,mplant)
    CHARACTER(len=20),DIMENSION(3) :: A2
    ! 2 dim arrays (npt,mlitter)
    CHARACTER(len=20),DIMENSION(3) :: A3
    ! 2 dim arrays (npt,msoil)
    CHARACTER(len=20),DIMENSION(3) :: A4
    INTEGER*4 :: VID1(SIZE(A1)), VIDI1(SIZE(AI1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))

    mp4=INT(mp,fmp4)
    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'glai'
    A1(4) = 'clabile'
    A1(5) = 'psoillab'
    A1(6) = 'psoilsorb'
    A1(7) = 'psoilocc'
    A1(8) = 'frac_sapwood'
    A1(9) = 'sapwood_area'
    A1(10) = 'phen'
    A1(11) = 'aphen'
    A1(12) = 'nsoilmin'

    AI1(1) = 'phase'
    AI1(2) = 'doyphase3'


    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'
    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'
    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'

    ! Get File-Name
    WRITE(CYEAR, FMT='(I4)') CurYear + 1

    IF (LEN( TRIM(casafile%cnpepool) ) .GT. 0) THEN
       fname=TRIM(casafile%cnpepool)
    ELSE
       fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
            '_casa_rst.nc'
    ENDIF
    ! Create NetCDF file:
    STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    WRITE(*,*) 'writing casa restart', fname
    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle", icycle  )
    IF ( CASAONLY ) THEN
       dum = 'CASA-ONLY run'
    ELSE
       dum = 'CABLE-CASA coupled run'
    ENDIF
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

    ! Define dimensions:
    ! Land (number of points)
    STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(AI1)
       STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID/),VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A3)
       STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID/),VID3(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A4)
       STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID/),VID4(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT LAT / LON
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), casamet%lat )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), casamet%lon )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT VARS
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), casamet%glai )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), casapool%clabile )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), casaflux%frac_sapwood )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), casaflux%sapwood_area )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), phen%phen )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), phen%aphen )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casapool%Nsoilmin )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), phen%phase )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), phen%doyphase(:,3) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), casapool%psoillab )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), casapool%psoilsorb )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), casapool%psoilocc )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant  )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  END SUBROUTINE WRITE_CASA_RESTART_NC

#ifndef UM_BUILD
  SUBROUTINE READ_CASA_RESTART_NC (  casamet, casapool, casaflux,phen )

    USE CASAVARIABLE
    USE phenvariable
    USE CABLE_COMMON_MODULE
    USE casa_ncdf_module, ONLY: HANDLE_ERR
    USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE, r_2, mp
    USE netcdf

    IMPLICIT NONE

    !INTEGER, INTENT(in)    :: YEAR
    TYPE (casa_met) , INTENT(inout) :: casamet
    TYPE (casa_pool), INTENT(inout) :: casapool
    TYPE (casa_flux), INTENT(inout) :: casaflux
    TYPE (phen_variable),       INTENT(INOUT) :: phen

    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS, i
    INTEGER*4   :: FILE_ID, dID, land_dim, mp_dim, ml_dim, ms_dim, mw_dim
    CHARACTER :: FRST_IN*99, CYEAR*4, CDATE*12, RSTDATE*12, FNAME*99

    ! ! 1 dim arrays (npt )
    ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
    !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
    ! ! 2 dim arrays (npt,mplant)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
    ! ! 2 dim arrays (npt,mlitter)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
    ! ! 2 dim arrays (npt,msoil)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)
    REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
    REAL(r_2)                         :: TMP2(mp,mplant),TMP3(mp,mlitter),TMP4(mp,msoil)

    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(12) :: A1
    CHARACTER(len=20),DIMENSION(2) :: AI1
    ! 2 dim arrays (npt,mplant)
    CHARACTER(len=20),DIMENSION(3) :: A2
    ! 2 dim arrays (npt,mlitter)
    CHARACTER(len=20),DIMENSION(3) :: A3
    ! 2 dim arrays (npt,msoil)
    CHARACTER(len=20),DIMENSION(3) :: A4
    ! 2-d array (npt,mwood)
    CHARACTER(len=20),DIMENSION(3) :: A5

    INTEGER :: VID1(SIZE(A1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))
    LOGICAL            ::  EXISTFILE, EXISTFILE1
    mp4=INT(mp,fmp4)
    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'glai'
    A1(4) = 'clabile'
    A1(5) = 'psoillab'
    A1(6) = 'psoilsorb'
    A1(7) = 'psoilocc'
    A1(8) = 'frac_sapwood'
    A1(9) = 'sapwood_area'
    A1(10) = 'phen'
    A1(11) = 'aphen'
    A1(12) = 'nsoilmin'

    AI1(1) = 'phase'
    AI1(2) = 'doyphase3'

    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'
    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'
    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'

    A5(1) = 'cwoodprod'
    A5(2) = 'nwoodprod'
    A5(3) = 'pwoodprod'



    !fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
    !       '_casa_rst.nc'
    fname =  TRIM(casafile%cnpipool)
    INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
    IF (EXISTFILE) THEN
       STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       PRINT *, 'initial pool from restart file: ', fname
    ELSE
       WRITE(*,*) 'CASA restart file:', TRIM(fname), ' does not exist'
       fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
            '_casa_rst.nc'
       INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE1 )
       IF (EXISTFILE1) THEN
          STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          PRINT *, 'initial pool from restart file: ', fname
       ELSE
          WRITE(*,*) 'CASA restart file:', TRIM(fname), ' does not exist either'
          WRITE(*,*) 'Set cable_user%CASA_fromZero to true to initialise without restart file.'
          WRITE(*,*) 'Otherwise set casafile%cnpipool to netcdf restart file name in cable.nml'
          STOP
       ENDIF
    ENDIF

    ! TIME
    STATUS = NF90_GET_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", RSTDATE )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$
    WRITE(CYEAR, FMT="(I4)") CurYear
    CDATE = '01/01/'//CYEAR
    ! compare current year with restart year (only for non-site type met data)
    IF ( CDATE .NE. RSTDATE .AND. &
         TRIM(cable_user%MetType).NE.'' .AND. TRIM(cable_user%MetType).NE.'site' ) THEN
       WRITE(*,*)"Restart Date in rst file doesn't match start date of Run!"
       WRITE(*,*)"File: "//RSTDATE//' Run: '//CDATE
       ! STOP
    ENDIF

    ! DIMS
    STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'mplant', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mp_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'mlitter', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ml_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'msoil', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ms_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    if(l_landuse) then
       STATUS = NF90_INQ_DIMID( FILE_ID, 'mwood', dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mw_dim )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    endif

    IF ( land_dim .NE. SIZE(casamet%lon) .OR. mp_dim .NE. mplant .OR. &
         ml_dim   .NE. mlitter             .OR. ms_dim .NE. msoil ) THEN
       WRITE(*,*)"Dimension misfit!"
       WRITE(*,*)"Restart file      Run"
       WRITE(*,*)"# points  ",land_dim,"     ",SIZE(casamet%lon)
       WRITE(*,*)"# mplant  ",mp_dim,"     ",mplant
       WRITE(*,*)"# mlitter ",ml_dim,"     ",mlitter
       WRITE(*,*)"# msoil   ",ms_dim,"     ",msoil
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

    ! CHECK FOR VALID LONS

    ! READ 1-dimensional fields
    DO i = 3, SIZE(A1)
       STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A1(i)))
       CASE ('glai'      ) ; casamet%glai       = TMP
       CASE ('clabile'   ) ; casapool%clabile   = TMP
       CASE ('frac_sapwood' ) ; casaflux%frac_sapwood  = TMP
       CASE ( 'sapwood_area' ) ; casaflux%sapwood_area  = TMP
       CASE ( 'phen' ) ; phen%phen  = TMP
       CASE ( 'aphen' ) ; phen%aphen  = TMP
       CASE ( 'nsoilmin' ) ; casapool%Nsoilmin  = TMP
       END SELECT
    END DO
    IF (icycle==3) THEN
       DO i = 3, SIZE(A1)
          STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A1(i)))
          CASE ('psoillab'  ) ; casapool%psoillab  = TMP
          CASE ('psoilsorb' ) ; casapool%psoilsorb = TMP
          CASE ('psoilocc'  ) ; casapool%psoilocc  = TMP
          END SELECT
       END DO
    ENDIF

    DO i = 1, SIZE(AI1)
       STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(AI1(i)))
       CASE ( 'phase' ) ; phen%phase  = TMP
       CASE ( 'doyphase3' ) ; phen%doyphase(:,3)  = TMP
       END SELECT
    END DO

    ! READ 2-dimensional fields (mplant)
    DO i = 1, SIZE(A2)
       STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A2(i)))
       CASE ('cplant' ) ; casapool%cplant = TMP2
       CASE ('nplant' ) ; casapool%nplant = TMP2
       END SELECT
    END DO


    IF (icycle==3) THEN
       DO i = 1, SIZE(A2)
          STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A2(i)))
          CASE ('pplant' ) ; casapool%pplant = TMP2
          END SELECT
       END DO
    ENDIF

    ! READ 2-dimensional fields (mlitter)
    DO i = 1, SIZE(A3)
       STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A3(i)))
       CASE ('clitter' ) ; casapool%clitter = TMP3
       CASE ('nlitter' ) ; casapool%nlitter = TMP3
       END SELECT
    END DO

    IF (icycle==3) THEN

       DO i = 1, SIZE(A3)
          STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A3(i)))
          CASE ('plitter' ) ; casapool%plitter = TMP3
          END SELECT
       END DO


    ENDIF

    ! READ 2-dimensional fields (msoil)
    DO i = 1, SIZE(A4)
       STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       SELECT CASE ( TRIM(A4(i)))
       CASE ('csoil' ) ; casapool%csoil = TMP4
       CASE ('nsoil' ) ; casapool%nsoil = TMP4
       END SELECT
    END DO

    if(l_landuse) then

       DO i = 1, SIZE(A5)
          STATUS = NF90_INQ_VARID( FILE_ID, A5(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          SELECT CASE ( TRIM(A4(i)))
             CASE ('cwoodprod' ) ; casapool%cwoodprod = TMP4
             CASE ('nwoodprod' ) ; casapool%nwoodprod = TMP4
             CASE ('pwoodprod' ) ; casapool%pwoodprod = TMP4
          END SELECT
        END DO

    endif

    IF (icycle==3) THEN
       DO i = 1, SIZE(A4)
          STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A4(i)))
          CASE ('psoil' ) ; casapool%psoil = TMP4
          END SELECT
       END DO
    ENDIF

    STATUS = NF90_CLOSE( FILE_ID )

  END SUBROUTINE READ_CASA_RESTART_NC
#endif
  SUBROUTINE WRITE_CASA_OUTPUT_NC ( veg, casamet, casapool, casabal, casaflux, &
       CASAONLY, ctime, FINAL )

    USE CASAVARIABLE
    USE CABLE_COMMON_MODULE
  USE casa_ncdf_module, ONLY: HANDLE_ERR


    USE cable_def_types_mod, ONLY: veg_parameter_type

    USE netcdf

    IMPLICIT NONE

    TYPE (casa_met) ,   INTENT(in) :: casamet
    TYPE (casa_pool),   INTENT(in) :: casapool
    TYPE (casa_balance),INTENT(in) :: casabal
    TYPE (casa_flux),   INTENT(in) :: casaflux
    TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters

    INTEGER   :: STATUS, ctime
    INTEGER   :: land_ID, plnt_ID, litt_ID, soil_ID, t_ID, i
    LOGICAL   :: CASAONLY, FINAL
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! ! 1 dim arrays (mp )
    ! CHARACTER(len=20),DIMENSION(2), PARAMETER :: A0 = (/ 'latitude', 'longitude' /)
    ! ! 2 dim arrays (mp,t)
    ! CHARACTER(len=20),DIMENSION(44),PARAMETER :: A1 = (/ 'glai', 'clabile',      &
    !      'psoillab','psoilsorb','psoilocc', 'sumcbal','sumnbal','sumpbal','Cgpp',&
    !      'Cnpp','stemnpp','Crp','Crgplant','Nminfix','Plabuptake','Clabloss',    &
    !      'fraclabile','Cnep','Crsoil','Nmindep','Nminloss','Nminleach',          &
    !      'Nupland','Nlittermin','Nsmin','Nsnet','fNMinloss','fNMinleach','Pdep', &
    !      'pwea','Pleach','Ploss','Pupland','Plittermin','Psmin','Psimm','Psnet', &
    !      'fPleach','kPlab','kPsorb','kpocc','kmlabP','Psorbmax','FluxCtoco2'/)
    ! ! 3 dim arrays (mp,mplant,t)
    ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A2 = (/ 'cplant' , 'nplant' ,   &
    !      'pplantc','fracCalloc','fracNalloc','fracPalloc','kplant','Crmplant'/)
    ! ! 3 dim arrays (mp,mlitter,t)
    ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A3 = (/ 'clitter', 'nlitter',   &
    !      'plitter','klitter','fromL2CO2','FluxCtolitter','FluxNtolitter',        &
    !      'FluxPtolitter' /)
    ! ! 3 dim arrays (mp,msoil,t)
    ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A4 = (/ 'csoil','nsoil','psoil',&
    !      'ksoil','fromStoCO2','FluxCtosoil','FluxNtosoil','FluxPxtosoil'/)

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim arrays (mp,t)
    CHARACTER(len=20),DIMENSION(51):: A1
    ! 3 dim arrays (mp,mplant,t)
    CHARACTER(len=20),DIMENSION(8) :: A2
    ! 3 dim arrays (mp,mlitter,t)
    CHARACTER(len=20),DIMENSION(8) :: A3
    ! 3 dim arrays (mp,msoil,t)
    CHARACTER(len=20),DIMENSION(8) :: A4

    ! 4 dim arrays (mp,mlitter,mplant,t)
    CHARACTER(len=20),DIMENSION(1), PARAMETER :: A5 = (/ 'fromPtoL'/)
    ! 4 dim arrays (mp,msoil,mlitter,t)
    CHARACTER(len=20),DIMENSION(1), PARAMETER :: A6 = (/ 'fromLtoS'/)
    ! 4 dim arrays (mp,msoil,msoil,t)
    CHARACTER(len=20),DIMENSION(1), PARAMETER :: A7 = (/ 'fromStoS'/)

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VID2(SIZE(A2)),VID3(SIZE(A3))
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)),VID6(SIZE(A6)),VID7(SIZE(A7))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    LOGICAL   :: EXRST
    CHARACTER(len=50) :: RecordDimName
    ! temporary fixer: ypwang 10-6-2021
    REAL(r_2), DIMENSION(mp,mplant)   :: tempFCrmyear
    REAL(R_2), DIMENSION(mp)          :: tempFCrgyear

    tempFCrmyear(:,1)  = casabal%FCrmleafyear(:)   ! leaf
    tempFCrmyear(:,2)  = casabal%FCrmwoodyear(:)   ! wood
    tempFCrmyear(:,3)  = casabal%FCrmrootyear(:)   ! froot
    tempFCrgyear(:)    = casabal%FCrgrowyear(:)



    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'glai'
    A1(2) = 'clabile'
    A1(3) = 'psoillab'
    A1(4) = 'psoilsorb'
    A1(5) = 'psoilocc'
    A1(6) = 'sumcbal'
    A1(7) = 'sumnbal'
    A1(8) = 'sumpbal'
    A1(9) = 'Cgpp'
    A1(10) = 'Cnpp'
    A1(11) = 'stemnpp'
    A1(12) = 'Crp'
    A1(13) = 'Crgplant'
    A1(14) = 'Nminfix'
    A1(15) = 'Plabuptake'
    A1(16) = 'Clabloss'
    A1(17) = 'fraclabile'
    A1(18) = 'Cnep'
    A1(19) = 'Crsoil'
    A1(20) = 'Nmindep'
    A1(21) = 'Nminloss'
    A1(22) = 'Nminleach'
    A1(23) = 'Nupland'
    A1(24) = 'Nlittermin'
    A1(25) = 'Nsmin'
    A1(26) = 'Nsimm'
    A1(27) = 'Nsnet'
    A1(28) = 'fNMinloss'
    A1(29) = 'Pdep'
    A1(30) = 'pwea'
    A1(31) = 'Pleach'
    A1(32) = 'Ploss'
    A1(33) = 'Pupland'
    A1(34) = 'Plittermin'
    A1(35) = 'Psmin'
    A1(36) = 'Psimm'
    A1(37) = 'Psnet'
    A1(38) = 'fPleach'
    A1(39) = 'kPlab'
    A1(40) = 'kPsorb'
    A1(41) = 'kpocc'
    A1(42) = 'kmlabP'
    A1(43) = 'Psorbmax'
    A1(44) = 'FluxCtoco2'
    A1(45) = 'FCgppyear'
    A1(46) = 'FCrpyear'
    A1(47) = 'FCnppyear'
    A1(48) = 'FCrsyear'
    A1(49) = 'FCNeeyear'
    A1(50) = 'vcmax'
    A1(51) = 'Nsoilmin'

    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'
    A2(4) = 'fracCalloc'
    A2(5) = 'fracNalloc'
    A2(6) = 'fracPalloc'
    A2(7) = 'kplant'
    A2(8) = 'Crmplant'

    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'
    A3(4) = 'klitter'
    A3(5) = 'fromL2CO2'
    A3(6) = 'FluxCtolitter'
    A3(7) = 'FluxNtolitter'
    A3(8) = 'FluxPtolitter'

    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'
    A4(4) = 'ksoil'
    A4(5) = 'fromStoCO2'
    A4(6) = 'FluxCtosoil'
    A4(7) = 'FluxNtosoil'
    A4(8) = 'FluxPxtosoil'


    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       IF (TRIM(cable_user%MetType).NE.'' ) THEN

          WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
          IF (CABLE_USER%YEARSTART.LT.1000.AND.CABLE_USER%YEAREND.LT.1000) THEN
             WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
          ELSEIF (CABLE_USER%YEARSTART.LT.1000) THEN
             WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
          ENDIF
          fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
               TRIM(dum)//'_casa_out.nc'
       ELSE
          ! site data
          fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_casa_out.nc'
       ENDIF
       INQUIRE( FILE=TRIM( fname ), EXIST=EXRST )
       EXRST = .FALSE.
       IF ( EXRST ) THEN
          STATUS = NF90_open(fname, mode=nf90_write, ncid=FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          CALL1 = .FALSE.

          STATUS = nf90_inq_dimid(FILE_ID, 'time', t_id)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          !CRM        status = nf90_inquire_dimension(FILE_ID, t_id,name = RecordDimName, len = CNT)
          !CRM        if (status /= nf90_noerr) call handle_err(status)
          !CRM        CNT = CNT+1

          STATUS = nf90_inq_varid(FILE_ID, 'time', VIDTime)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          DO i = 1, SIZE(A0)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A0(i)),VID0(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A1)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A1(i)), VID1(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A2)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A2(i)) , VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A3)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A3(i)) ,VID3(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A4)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A4(i)) ,VID4(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A5)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A5(i)), VID5(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A6)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A6(i)), VID6(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A7)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A7(i)),VID7(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

       ELSE
          ! Create NetCDF file:
          STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Put the file in define mode:
          STATUS = NF90_redef(FILE_ID)

          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle"   , icycle  )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )
          IF ( CASAONLY ) THEN
             dum = 'CASA-ONLY run'
          ELSE
             dum = 'CABLE-CASA coupled run'
          ENDIF
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

          ! Define dimensions:
          ! Land (number of points)
          STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define variables
          STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),VIDtime )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          DO i = 1, SIZE(A0)
             STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A1)
             STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID,t_ID/),VID1(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A3)
             STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID,t_ID/),VID3(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A4)
             STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID,t_ID/),VID4(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A5)
             STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT, &
                  (/land_ID,litt_ID,plnt_ID,t_ID/),VID5(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A6)
             STATUS = NF90_def_var(FILE_ID,TRIM(A6(i)) ,NF90_FLOAT, &
                  (/land_ID,soil_ID,litt_ID,t_ID/),VID6(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A7)
             STATUS = NF90_def_var(FILE_ID,TRIM(A7(i)) ,NF90_FLOAT, &
                  (/land_ID,soil_ID,soil_ID,t_ID/),VID7(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          ! End define mode:
          STATUS = NF90_enddef(FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


          ! PUT LAT / LON ( mp )
          STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), casamet%lat )
          IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

          STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), casamet%lon )
          IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

          CALL1 = .FALSE.
       ENDIF !( EXRST )
    ENDIF

    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    SELECT CASE(icycle)
    CASE(1)
       ! PUT 2D VARS ( mp, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), REAL(casaflux%Cnpp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), tempFCrgyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casaflux%Crsoil,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), tempFCrmyear,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       ! PUT 3D VARS ( mp, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    CASE(2)


       ! PUT 2D VARS ( mp, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), casapool%psoillab,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), casapool%psoilsorb,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), casapool%psoilocc,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), casabal%sumpbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), casaflux%Cnpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), tempFCrgyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casabal%FNfixyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       ! STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Nminuptake, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       ! IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casabal%FPupyear, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casabal%FCrsyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), casabal%FNdepyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), casabal%FNlossyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), casabal%FNleachyear,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), casabal%FNleachyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), casaflux%Nlittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), casaflux%Nsmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), casaflux%Nsimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), casabal%FNsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), casaflux%fNminloss,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), casabal%FPdustyear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), casabal%FPweayear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), casabal%FPleachyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), casabal%FPlossyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), casaflux%Pupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(34), casaflux%Plittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(35), casaflux%Psmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(36), casaflux%Psimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(37), casabal%FPsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(38), casaflux%fPleach,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(39), casaflux%kplab,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(40), casaflux%kpsorb,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(41), casaflux%kpocc,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(42), casaflux%kmlabP,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(43), casaflux%Psorbmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(44), casaflux%FluxCtoco2,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(51), casapool%Nsoilmin,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), casaflux%fracNalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), casaflux%fracPalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), tempFCrmyear,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(7), casaflux%FluxNtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(8), casaflux%FluxPtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), casaflux%FluxNtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(8), casaflux%FluxPtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    CASE(3)
       ! PUT 2D VARS ( mp, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), casapool%psoillab,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), casapool%psoilsorb,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), casapool%psoilocc,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), casabal%sumpbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), casaflux%Cnpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), tempFCrgyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casabal%FNfixyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       !  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       !  STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Nminuptake, start=(/ 1, CNT /), count=(/ mp, 1 !/) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casabal%FPupyear, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casabal%FCrsyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), casabal%FNdepyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), casabal%FNlossyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), casabal%FNleachyear,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), casabal%FNupyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), casaflux%Nlittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), casaflux%Nsmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), casaflux%Nsimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), casabal%FNsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), casaflux%fNminloss,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), casabal%FPdustyear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), casabal%FPweayear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), casabal%FPleachyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), casabal%FPlossyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), casaflux%Pupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(34), casaflux%Plittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(35), casaflux%Psmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(36), casaflux%Psimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(37), casabal%FPsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(38), casaflux%fPleach,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(39), casaflux%kplab,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(40), casaflux%kpsorb,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(41), casaflux%kpocc,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(42), casaflux%kmlabP,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(43), casaflux%Psorbmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(44), casaflux%FluxCtoCo2,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(51), casapool%Nsoilmin,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), casaflux%fracNalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), casaflux%fracPalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), tempFCrmyear,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(7), casaflux%FluxNtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(8), casaflux%FluxPtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), casaflux%FluxNtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(8), casaflux%FluxPtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    END SELECT

    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " Casa Output written to ",fname
    ENDIF

  END SUBROUTINE WRITE_CASA_OUTPUT_NC
#endif

END MODULE casa_offline_inout_module
