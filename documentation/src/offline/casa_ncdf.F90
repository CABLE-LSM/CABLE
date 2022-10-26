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
! Purpose: handles additional, dynamically decided diagnostic output from model.
!          permanently used for bitwise identical testing. more applications
!          will follow.
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Currently stripped down version of cable_diag here. will be
!          re-implemented in time.
!
! ==============================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ USE this module in any subr. you wish to write vars from.             +++!
!+++ x is typically the number of landpoints(tiles). binary file is        +++!
!+++ then appended every timestep with the new foo(x_i)                    +++!
!+++                                                                       +++!
!+++ CALL syntax:                                                          +++!
!+++                                                                       +++!
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, var1 )     +++!
!+++                                                                       +++!
!+++ output binaries can be interpreted from the command line              +++!
!+++ using a suite of tools. Currently, only zero_diff.ksh is supported.   +++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!CABLE_LSM:This has to be commented for offline
!#define UM_BUILD YES
MODULE casa_ncdf_module
   
   IMPLICIT NONE

#ifndef UM_BUILD
  interface put_var_nc
     module procedure put_var_ncr1, put_var_ncr2, put_var_ncr3
  end interface put_var_nc

  interface get_var_nc
     module procedure get_var_ncr2, get_var_ncr3
  end interface get_var_nc

#endif

CONTAINS

#ifndef UM_BUILD
  subroutine def_dims(nd, ncid, dimID, dim_len, dim_name )
    use netcdf
    implicit none
    integer, intent(in) :: nd, ncid
    character(len=*), dimension(:), intent(in) :: dim_name
    integer, dimension(:), intent(out) :: dimID
    integer, dimension(:), intent(in) :: dim_len
    integer :: j, ncok

    do j=1, nd
       ncok = NF90_DEF_DIM(ncid, trim(dim_name(j)), dim_len(j), dimID(j) )
       if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def dim ', dim_name(j))
    enddo

    return
  end subroutine def_dims




  subroutine def_vars(nv, ncid,  xtype, dimID, var_name,varID )
    use netcdf
    implicit none
    integer, intent(in) :: nv, ncid, xtype
    integer, dimension(:), intent(in) :: dimID
    integer, dimension(:), intent(inout) :: varID
    character(len=*), dimension(:), intent(in) :: var_name
    integer :: j, ncok

    ! lat
    ncok = NF90_DEF_VAR( ncid, trim(var_name(1)), xtype, &
         (/ dimID(1) /), varID(1))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(1))

    ! lon
    ncok = NF90_DEF_VAR(ncid, trim(var_name(2)), xtype, &
         (/ dimID(1) /), varID(2))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(2))

    ! tairk
    ncok = NF90_DEF_VAR(ncid, trim(var_name(3)), xtype, &
         (/ dimID(1), dimID(3) /), varID(3))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(3))

    !tsoil
    ncok = NF90_DEF_VAR(ncid, trim(var_name(4)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(4))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(4))

    ! moist
    ncok = NF90_DEF_VAR(ncid, trim(var_name(5)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(5))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(5))

    !cgpp
    ncok = NF90_DEF_VAR(ncid, trim(var_name(6)), xtype, &
         (/ dimID(1), dimID(3)/), varID(6))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(6))

    !crmplant
    ncok = NF90_DEF_VAR(ncid, trim(var_name(7)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(7))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(7))

    !phenphase
    ncok = NF90_DEF_VAR(ncid, trim(var_name(8)), xtype, &
         (/ dimID(1), dimID(3)/), varID(8))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(8))

    !doyphase1
    ncok = NF90_DEF_VAR(ncid, trim(var_name(9)), xtype, &
         (/ dimID(1), dimID(3)/), varID(9))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(9))

    !doyphase2
    ncok = NF90_DEF_VAR(ncid, trim(var_name(10)), xtype, &
         (/ dimID(1), dimID(3)/), varID(10))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(10))

    !doyphase3
    ncok = NF90_DEF_VAR(ncid, trim(var_name(11)), xtype, &
         (/ dimID(1), dimID(3)/), varID(11))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(11))

    !doyphase4
    ncok = NF90_DEF_VAR(ncid, trim(var_name(12)), xtype, &
         (/ dimID(1), dimID(3)/), varID(12))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(12))


    !mtemp
    ncok = NF90_DEF_VAR(ncid, trim(var_name(13)), xtype, &
         (/ dimID(1),dimID(3)/), varID(13))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(13))

    !Ndep
    ncok = NF90_DEF_VAR(ncid, trim(var_name(14)), xtype, &
         (/ dimID(1),dimID(3)/), varID(14))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(14))

    return
  end subroutine def_vars

  subroutine def_var_atts( ncfile_in, ncid, varID )
    use netcdf
    implicit none
    character(len=*), intent(in) :: ncfile_in
    integer, intent(in):: ncid       ! netcdf file ID
    integer, dimension(:), intent(in) :: varID ! (1) ~ tvair, (2) ~ pmb
    integer :: j, ncok
    character(len=10) dummy

    write(dummy,11) varID(1)
11  format(i2)
    ncok = NF90_PUT_ATT(ncid, nf90_global, "Title", "Forcing for define_air subroutine")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', ncfile_in)
    ncok = NF90_PUT_ATT(ncid, varID(3), "longname", "air temperature within canopy")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', dummy)
    ncok = NF90_PUT_ATT(ncid, varID(3), "units", "K")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', dummy)

    write(dummy,11) varID(2)


    return
  end subroutine def_var_atts


  subroutine put_var_ncr1(ncid, var_name, var )
    use netcdf
    use cable_def_types_mod, only : mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real, dimension(:),intent(in) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID,j

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1/), &
         count=(/mp/) )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

  end subroutine put_var_ncr1


  subroutine put_var_ncr2(ncid, var_name, var, n_call )
    use netcdf
    use cable_def_types_mod, only : r_2, mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real(r_2), dimension(:),intent(in) :: var
    integer, intent(in) :: ncid, n_call
    integer :: ncok, varID

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,n_call /), &
         count=(/mp,1/) )

    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

  end subroutine put_var_ncr2

  !soil vars
  subroutine put_var_ncr3(ncid, var_name, var, n_call, nl)
    use netcdf
    use cable_def_types_mod, only : r_2, mp, ms
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:,:),intent(in) :: var
    integer, intent(in) :: ncid, n_call, nl
    integer :: ncok, varID

    ncok = NF90_INQ_VARID( ncid, var_name, varId )
    IF( ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name )

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,1,n_call /), &
         count=(/mp,nl,1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

    return
  end subroutine put_var_ncr3



  subroutine get_var_ncr2(ncid, var_name, var, n_call )
    use netcdf
    use cable_def_types_mod, only : r_2,mp
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:),intent(out) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID, n_call
    real, dimension(mp) :: temp

    temp = 0.

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)
    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,n_call/), &
         count=(/mp,1/) )

    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'getting var ', var_name)

    var = real( temp, r_2 )
  end subroutine get_var_ncr2

  subroutine get_var_ncr3(ncid, var_name, var, n_call, nl )
    use netcdf
    use cable_def_types_mod, only : r_2, mp, ms
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:,:),intent(out) :: var
    integer, intent(in) :: ncid, n_call, nl
    integer :: ncok, varID
    real, dimension(mp,1:nl) :: temp

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,1,n_call /), &
         count=(/mp, nl, 1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)
    var = real( temp, r_2 )
  end subroutine get_var_ncr3

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




  subroutine stderr_nc(status,message, var)
    use netcdf
    character(len=*), intent(in) :: message, var
    INTEGER, INTENT(IN) :: status
    character(len=7) :: err_mess
    err_mess = 'ERROR:'
    print *, (err_mess//message), var
    PRINT*,NF90_STRERROR(status)
    stop
  end subroutine stderr_nc
#endif
  SUBROUTINE YMDHMS2DOYSOD( YYYY,MM,DD,HOUR,MINUTE,SECOND,DOY,SOD )
USE cable_common_module, ONLY: IS_LEAPYEAR

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
USE cable_common_module, ONLY: IS_LEAPYEAR

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

  FUNCTION IS_CASA_TIME(iotype, yyyy, ktau, kstart, koffset, kend, ktauday, logn)

  USE cable_common_module, ONLY: CABLE_USER 
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



END MODULE casa_ncdf_module



