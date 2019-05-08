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
! Purpose: Writing routines for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: No significant changes since v1.4b except addition of extra variables
!
!
! ==============================================================================
!
! CALLed from:    cable_initialise.F90
!                 cable_input.F90
!
! MODULEs used:   cable_abort_module
!                 cable_def_types_mod
!                 cable_IO_vars_module
!                 netcdf
!
! CALLs:          define_output_variable_r1
!                 define_output_variable_r2
!                 define_output_parameter_r1
!                 define_output_parameter_r2
!                 write_output_variable_r1
!                 write_output_variable_r2
!                 write_output_parameter_r1
!                 write_output_parameter_r1d
!                 write_output_parameter_r2
!                 write_output_parameter_r2d
!
!
!
! Notes: Single precision netcdf writes are forced to single precision here
!        (using REAL(,4)) in case of compilation with -r8

MODULE cable_write_module


  USE cable_abort_module
  USE cable_def_types_mod
  USE cable_IO_vars_module, ONLY: landpt, patch, max_vegpatches, parID_type,           &
       metGrid, land_x, land_y, logn, output,               &
       xdimsize, ydimsize, check, mask
  USE netcdf
  IMPLICIT NONE
  PRIVATE
  PUBLIC define_ovar, write_ovar, otmp1, otmp1l, otmp2lt, otmp2xy, otmp2lp,    &
       otmp2ls, otmp2lpc, otmp2lsc, otmp2lsf, otmp2lr, otmp2lsn, otmp3xyt,   &
       otmp3lpt, otmp3lst, otmp3lsnt, otmp3lrt, otmp3lpct, otmp3lsct,        &
       otmp3xyp, otmp3xys, otmp3xypc, otmp3xysc, otmp3lps, otmp3lppc,        &
       otmp3lpsc, otmp3xysf, otmp3lpr, otmp3lpsn, otmp4xypt, otmp4xyzt,      &
       otmp4xyst, otmp4xysnt, otmp4xyrt, otmp4xypct, otmp4xysct, otmp4lpst,  &
       otmp4lpsnt, otmp4lprt, otmp4lpsct, otmp4lppct, otmp4xyps,             &
       otmp4xyppc, otmp4xypsc, otmp5xypst, otmp5xypsnt, otmp5xyprt,          &
       otmp5xyppct, otmp5xypsct, nullify_write
  INTERFACE define_ovar
     ! Defines an output variable in the output netcdf file. Units, long name,
     ! variable, dimensions etc are created.
     MODULE PROCEDURE define_output_variable_r1
     MODULE PROCEDURE define_output_variable_r2
     MODULE PROCEDURE define_output_parameter_r1
     MODULE PROCEDURE define_output_parameter_r2
  END INTERFACE
  INTERFACE write_ovar
     ! Writes a single time step of an output variable to the output netcdf
     ! file
     MODULE PROCEDURE write_output_variable_r1
     MODULE PROCEDURE write_output_variable_r2
     MODULE PROCEDURE write_output_parameter_r1
     MODULE PROCEDURE write_output_parameter_r1d
     MODULE PROCEDURE write_output_parameter_r2
     MODULE PROCEDURE write_output_parameter_r2d
  END INTERFACE

  INTEGER :: ncmissingi = -9999999
  INTEGER :: ok ! netcdf file read status

  ! Temporary variables of same dimension as variables in netcdf file;
  ! e.g. 'o'utput 'tmp'orary with '2' dimensions: 'l'and and 't'ime -> otmp2lt
  ! Other dimension abbrevs: 'x','y','z','p'atch,'s'oil,'sn'ow,
  ! 'r'adiation,'p'lant 'c'arbon,'s'oil 'c'arbon,'s'urface 'f'raction
  REAL, POINTER, DIMENSION(:) :: otmp1, otmp1l
  REAL, POINTER, DIMENSION(:, :) :: otmp2lt, otmp2xy, otmp2lp, otmp2ls,   &
       otmp2lpc, otmp2lsc, otmp2lsf,         &
       otmp2lr, otmp2lsn
  REAL, POINTER, DIMENSION(:, :, :) :: otmp3xyt, otmp3lpt, otmp3lst,      &
       otmp3lsnt, otmp3lrt, otmp3lpct,    &
       otmp3lsct, otmp3xyp, otmp3xys,     &
       otmp3xypc, otmp3xysc, otmp3lps,    &
       otmp3lppc, otmp3lpsc, otmp3xysf,   &
       otmp3lpr, otmp3lpsn, otmp3xyr
  REAL, POINTER, DIMENSION(:, :, :, :) :: otmp4xypt, otmp4xyzt,           &
       otmp4xyst, otmp4xysnt,          &
       otmp4xyrt, otmp4xypct,          &
       otmp4xysct, otmp4lpst,          &
       otmp4lpsnt, otmp4lprt,          &
       otmp4lpsct, otmp4lppct,         &
       otmp4xyps, otmp4xyppc,          &
       otmp4xypsc, otmp4xypr
  REAL, POINTER, DIMENSION(:, :, :, :, :) :: otmp5xypst, otmp5xypsnt,     &
       otmp5xyprt, otmp5xyppct,     &
       otmp5xypsct
  REAL :: ncmissingr = -1.0e+33

CONTAINS

  ! Nullify all temporary pointers so that one can query associated(pointer)
  SUBROUTINE nullify_write()
    IMPLICIT NONE

    NULLIFY(otmp1)
    NULLIFY(otmp1l)

    NULLIFY(otmp2lt)
    NULLIFY(otmp2xy)
    NULLIFY(otmp2lp)
    NULLIFY(otmp2ls)
    NULLIFY(otmp2lpc)
    NULLIFY(otmp2lsc)
    NULLIFY(otmp2lsf)
    NULLIFY(otmp2lr)
    NULLIFY(otmp2lsn)

    NULLIFY(otmp3xyt)
    NULLIFY(otmp3lpt)
    NULLIFY(otmp3lst)
    NULLIFY(otmp3lsnt)
    NULLIFY(otmp3lrt)
    NULLIFY(otmp3lpct)
    NULLIFY(otmp3lsct)
    NULLIFY(otmp3xyp)
    NULLIFY(otmp3xys)
    NULLIFY(otmp3xypc)
    NULLIFY(otmp3xysc)
    NULLIFY(otmp3lps)
    NULLIFY(otmp3lppc)
    NULLIFY(otmp3lpsc)
    NULLIFY(otmp3xysf)
    NULLIFY(otmp3lpr)
    NULLIFY(otmp3lpsn)
    NULLIFY(otmp3xyr)

    NULLIFY(otmp4xypt)
    NULLIFY(otmp4xyzt)
    NULLIFY(otmp4xyst)
    NULLIFY(otmp4xysnt)
    NULLIFY(otmp4xyrt)
    NULLIFY(otmp4xypct)
    NULLIFY(otmp4xysct)
    NULLIFY(otmp4lpst)
    NULLIFY(otmp4lpsnt)
    NULLIFY(otmp4lprt)
    NULLIFY(otmp4lpsct)
    NULLIFY(otmp4lppct)
    NULLIFY(otmp4xyps)
    NULLIFY(otmp4xyppc)
    NULLIFY(otmp4xypsc)
    NULLIFY(otmp4xypr)

    NULLIFY(otmp5xypst)
    NULLIFY(otmp5xypsnt)
    NULLIFY(otmp5xyprt)
    NULLIFY(otmp5xyppct)
    NULLIFY(otmp5xypsct)

  END SUBROUTINE nullify_write

  SUBROUTINE define_output_variable_r1(ncid, varID, vname,                     &
       vunits, longname, writepatch,           &
       dimswitch, xID, yID, zID, landID,       &
       patchID, tID)
    ! Subroutine for defining a real valued 1D variable
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(OUT) :: varID ! variable's netcdf ID
    ! netcdf dimension IDs
    INTEGER, INTENT(IN) :: xID, yID, zID, landID, patchID, tID
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: vunits ! variable units
    CHARACTER(LEN=*), INTENT(IN) :: longname ! full variable name
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF(output%grid(1:3) == 'mas' .OR.                                          &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.                 &
         output%grid(1:3) == 'ALM') THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND.                                  &
            (.NOT. output%grid(1:3) == 'ALM')) THEN
          WRITE(logn, *) 'Writing '//vname//                                     &
               ' to output file using mask grid with patch-specific info'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/xID, yID, patchID, tID/), &
               varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(INTERFACE define_ovar)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(.NOT.ASSOCIATED(otmp4xypt))                                         &
               ALLOCATE(otmp4xypt(xdimsize, ydimsize, max_vegpatches, 1))
       ELSE ! only grid point values, no patch-specific info
          ! If this is an ALMA 4D surface variable
          ! AND the user has forced the grid type as ALMA:
          IF(dimswitch == 'ALMA' .AND. output%grid(1:3) == 'ALM') THEN
             WRITE(logn, *) 'Writing '//vname//' to output file using mask grid'
             ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/xID, yID, zID, tID/),   &
                  varID)
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF(.NOT.ASSOCIATED(otmp4xyzt))                                       &
                  ALLOCATE(otmp4xyzt(xdimsize, ydimsize, 1, 1))
          ELSE ! normal x-y-t mask grid
             WRITE(logn, *) 'Writing '//vname//' to output file using mask grid'
             ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/xID, yID, tID/), varID)
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF(.NOT.ASSOCIATED(otmp3xyt))ALLOCATE(otmp3xyt(xdimsize, ydimsize, 1))
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(INTERFACE define_ovar)')
       END IF
    ELSE IF(output%grid(1:3) == 'lan'                                          &
         .OR.(output%grid(1:3) == 'def' .AND. metGrid == 'land')) THEN
       ! Should patch-specific info be written for this variable?
       IF(writepatch .OR. output%patch) THEN
          WRITE(logn, *) 'Writing '//vname//                                     &
               ' to output file using land grid with patch-specific info'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/landID, patchID, tID/),   &
               varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok,'Error defining '//vname//' variable in output file. '// &
               '(INTERFACE define_ovar)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF( .NOT. ASSOCIATED(otmp3lpt)) ALLOCATE(otmp3lpt(mland,               &
               max_vegpatches, 1))
       ELSE ! only grid point values, no patch-specific info
          WRITE(logn, *) 'Writing '//vname//' to output file using land grid'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/landID,tID/), varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(INTERFACE define_ovar)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF( .NOT. ASSOCIATED(otmp2lt)) ALLOCATE(otmp2lt(mland, 1))
       END IF
    ELSE
       CALL abort('Unknown grid specification (INTERFACE define_ovar)')
    END IF
    ! Define variable units:
    ok = NF90_PUT_ATT(ncid, varID, 'units', vunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(INTERFACE define_ovar)')
    ! Define long name:
    ok = NF90_PUT_ATT(ncid,varID, 'long_name', longname)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(INTERFACE define_ovar)')
    ! Define missing/fill values:
    ok = NF90_PUT_ATT(ncid, varID, '_FillValue', REAL(ncmissingr, 4))
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(INTERFACE define_ovar)')
    ! Define missing/fill values:
    ok = NF90_PUT_ATT(ncid, varID, 'missing_value', REAL(ncmissingr, 4))
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(INTERFACE define_ovar)')

  END SUBROUTINE define_output_variable_r1
  !=============================================================================
  SUBROUTINE define_output_variable_r2(ncid, varID, vname, vunits, longname,   &
       writepatch, dimswitch, xID, yID, zID,   &
       landID, patchID, othdimID, tID)
    ! Subroutine for defining a real valued 2D variable
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    ! netcdf dimension IDs
    INTEGER, INTENT(IN) :: xID, yID, zID, landID, patchID, tID
    INTEGER, INTENT(IN) :: othdimID ! ID of variable's second dimension
    INTEGER, INTENT(OUT) :: varID ! variable's netcdf ID
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: vunits ! variable units
    CHARACTER(LEN=*), INTENT(IN) :: longname ! full variable name
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF(output%grid(1:3) == 'mas' .OR.                                          &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.                 &
         output%grid(1:3) == 'ALM') THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND.                                  &
            ( .NOT. output%grid(1:3) == 'ALM')) THEN
          WRITE(logn, *) 'Writing '//vname//                                     &
               ' to output file using mask grid with patch-specific info'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/xID, yID, patchID,        &
               othdimID, tID/), varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(INTERFACE define_ovar)')
          IF(dimswitch == 'soil') THEN ! other dim is soil
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp5xypst))                                    &
                  ALLOCATE(otmp5xypst(xdimsize, ydimsize, max_vegpatches, ms, 1))
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp5xypsnt))                                   &
                  ALLOCATE(otmp5xypsnt(xdimsize, ydimsize, max_vegpatches, msn, 1))
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp5xyprt))                                    &
                  ALLOCATE(otmp5xyprt(xdimsize, ydimsize, max_vegpatches, nrb, 1))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp5xyppct))                                   &
                  ALLOCATE(otmp5xyppct(xdimsize, ydimsize, max_vegpatches, ncp, 1))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp5xypsct))                                   &
                  ALLOCATE(otmp5xypsct(xdimsize, ydimsize, max_vegpatches, ncs, 1))
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch// &
                  ' - in SUBROUTINE define_output_variable_r2')
          END IF
       ELSE ! only grid point values, no patch-specific info
          WRITE(logn, *) 'Writing '//vname//' to output file using mask grid'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/xID, yID, othdimID,       &
               tID/), varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(SUBROUTINE define_output_variable_r2)')
          IF(dimswitch == 'soil') THEN ! other dim is soil
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xyst))                                     &
                  ALLOCATE(otmp4xyst(xdimsize, ydimsize, ms, 1))
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xysnt))                                    &
                  ALLOCATE(otmp4xysnt(xdimsize, ydimsize, msn, 1))
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xyrt))                                     &
                  ALLOCATE(otmp4xyrt(xdimsize, ydimsize, nrb, 1))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xypct))                                    &
                  ALLOCATE(otmp4xypct(xdimsize, ydimsize, ncp, 1))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xysct))                                    &
                  ALLOCATE(otmp4xysct(xdimsize, ydimsize, ncs, 1))
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch// &
                  ' - in SUBROUTINE define_output_variable_r2')
          END IF
       END IF
    ELSE IF(output%grid(1:3) == 'lan'                                          &
         .OR. (output%grid(1:3) == 'def' .AND. metGrid == 'land')) THEN
       ! Should patch-specific info be written for this variable?
       IF(writepatch .OR. output%patch) THEN
          WRITE(logn, *) 'Writing '//vname//                                     &
               ' to output file using land grid with patch-specific info'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/landID, patchID,          &
               othdimID, tID/), varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(SUBROUTINE define_output_variable_r2)')
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok,'Error defining '//vname//' variable in output file. '// &
               '(SUBROUTINE define_output_variable_r2)')
          IF(dimswitch == 'soil') THEN ! other dim is soil
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4lpst))                                     &
                  ALLOCATE(otmp4lpst(mland, max_vegpatches, ms, 1))
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xysnt))                                    &
                  ALLOCATE(otmp4xysnt(mland, max_vegpatches, msn, 1))
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xyrt))                                     &
                  ALLOCATE(otmp4xyrt(mland, max_vegpatches, nrb, 1))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xypct))                                    &
                  ALLOCATE(otmp4xypct(mland, max_vegpatches, ncp, 1))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp4xysct))                                    &
                  ALLOCATE(otmp4xysct(mland, max_vegpatches, ncs, 1))
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch// &
                  ' - in SUBROUTINE define_output_variable_r2')
          END IF
       ELSE ! only grid point values, no patch-specific info
          WRITE(logn, *) 'Writing '//vname//' to output file using land grid'
          ok = NF90_DEF_VAR(ncid, vname, NF90_FLOAT, (/landID, othdimID, tID/),  &
               varID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                    &
               (ok, 'Error defining '//vname//' variable in output file. '// &
               '(SUBROUTINE define_output_variable_r2)')
          IF(dimswitch == 'soil') THEN ! other dim is soil
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp3lst)) ALLOCATE(otmp3lst(mland, ms, 1))
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp3lsnt)) ALLOCATE(otmp3lsnt(mland, msn, 1))
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp3lrt)) ALLOCATE(otmp3lrt(mland, nrb, 1))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF(.NOT.ASSOCIATED(otmp3lpct)) ALLOCATE(otmp3lpct(mland, ncp, 1))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             ! If not already allocated, allocate a temporary storage variable
             ! of this dim:
             IF( .NOT. ASSOCIATED(otmp3lsct)) ALLOCATE(otmp3lsct(mland, ncs, 1))
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch// &
                  ' - in SUBROUTINE define_output_variable_r2')
          END IF
       END IF
    ELSE
       CALL abort('Unknown grid specification (SUBROUTINE '//                   &
            'define_output_variable_r2)')
    END IF
    ! Define variable units:
    ok = NF90_PUT_ATT(ncid, varID, 'units', vunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(SUBROUTINE define_output_variable_r2)')
    ! Define long name:
    ok = NF90_PUT_ATT(ncid, varID, 'long_name', longname)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(SUBROUTINE define_output_variable_r2)')
    ! Define missing/fill values:
    ok = NF90_PUT_ATT(ncid, varID, '_FillValue', REAL(ncmissingr, 4))
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok,'Error defining '//vname//' variable attributes in output file. '// &
         '(INTERFACE define_ovar)')
    ! Define missing/fill values:
    ok = NF90_PUT_ATT(ncid, varID, 'missing_value', REAL(ncmissingr, 4))
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//vname//' variable attributes in output file. '// &
         '(INTERFACE define_ovar)')

  END SUBROUTINE define_output_variable_r2
  !=============================================================================
  SUBROUTINE define_output_parameter_r1(ncid, parID, pname, punits, longname,  &
       writepatch, dimswitch, xID, yID, zID,  &
       landID, patchID, restart)
    ! Subroutine for defining a real valued 1D parameter (time invariant)
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: xID, yID, zID, landID, patchID ! netcdf
    ! dimension IDs
    INTEGER, INTENT(OUT) :: parID ! variable's netcdf ID
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    LOGICAL, INTENT(IN), OPTIONAL :: restart ! are we writing to a restart file?                                                          ! dimension IDs
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: punits ! variable units
    CHARACTER(LEN=*), INTENT(IN) :: longname ! full variable name
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimension of parameter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF((output%grid(1:3) == 'mas' .OR.                                         &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.               &
         output%grid(1:3) == 'ALM') .AND. .NOT. PRESENT(restart)) THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND.                                 &
            (.NOT. output%grid(1:3) == 'ALM')) THEN
          WRITE(logn, *) 'Writing '//pname//                                   &
               ' to output file using mask grid with patch-specific info'
          IF(dimswitch(1:1) == 'r') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/xID, yID, patchID/)  &
                  , parID)
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_INT, (/xID, yID, patchID/)    &
                  , parID)
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok, 'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r1)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(.NOT. ASSOCIATED(otmp3xyp))                                       &
               ALLOCATE(otmp3xyp(xdimsize, ydimsize, max_vegpatches))
       ELSE ! only grid point values, no patch-specific info
          WRITE(logn, *) 'Writing '//pname//' to output file using mask grid'
          IF(dimswitch(1:1) == 'r') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/xID, yID/), parID)
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_INT, (/xID, yID/), parID)
          END IF
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(.NOT. ASSOCIATED(otmp2xy)) ALLOCATE(otmp2xy(xdimsize, ydimsize))
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok, 'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r1)')
       END IF
    ELSE IF(output%grid(1:3) == 'lan' .OR. (output%grid(1:3) == 'def' .AND.    &
         metGrid == 'land') .OR. PRESENT(restart)) THEN ! land-only grid
       ! Should patch-specific info be written for this variable?
       ! If this variable has been requested by user with patch-specific info
       ! (writepatch) OR all have been (output%patch) AND we're NOT writing
       ! a restart file (which uses a different technique to store patch info):
       IF((writepatch .OR. output%patch) .AND. .NOT. PRESENT(restart)) THEN
          WRITE(logn, *) 'Writing '//pname//                                   &
               ' to output file using land grid with patch-specific info'
          IF(dimswitch(1:2) == 're') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/landID, patchID/)    &
                  , parID)
          ELSE IF(dimswitch(1:2) == 'r2') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_DOUBLE, (/landID, patchID/)   &
                  , parID)
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_INT, (/landID, patchID/)      &
                  , parID)
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok, 'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r1)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(.NOT. ASSOCIATED(otmp2lp)) ALLOCATE(otmp2lp(mland, max_vegpatches))
       ELSE ! only grid point values without patch-specific info UNLESS a
          ! restart variable
          ! Restart file definitions will be directed to this part of interface.
          ! If not writing a restart file, report variable writing to log file:
          IF(.NOT. PRESENT(restart)) WRITE(logn, *) 'Writing '//pname//        &
               ' to output file using land grid'
          IF(dimswitch(1:2) == 're') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/landID/), parID)
          ELSE IF(dimswitch(1:2) == 'r2') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_DOUBLE, (/landID/), parID)
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_INT, (/landID/), parID)
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok,'Error defining '//pname//' variable in output or '// &
               'restart file. (SUBROUTINE define_output_parameter_r1)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dimension structure:
          IF(.NOT. ASSOCIATED(otmp1l)) ALLOCATE(otmp1l(mland))
       END IF
    ELSE
       CALL abort('Unknown grid specification '//                              &
            '(SUBROUTINE define_output_parameter_r1)')
    END IF
    ! Define variable units:
    ok = NF90_PUT_ATT(ncid, parID, 'units', punits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//pname//' variable attributes in '//             &
         'output file. (SUBROUTINE define_output_parameter_r1)')
    ! Define long name:
    ok = NF90_PUT_ATT(ncid, parID, 'long_name', longname)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//pname//' variable attributes in '//             &
         'output file. (SUBROUTINE define_output_parameter_r1)')
    ! Define missing/fill values:
    IF(dimswitch(1:1) == 'i') THEN
       ok = NF90_PUT_ATT(ncid, parID, '_FillValue', ncmissingi)
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
       ok = NF90_PUT_ATT(ncid, parID, 'missing_value', ncmissingi)
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
    ELSE IF(dimswitch(1:2) == 'r2') THEN
       ok = NF90_PUT_ATT(ncid, parID, '_FillValue', REAL(ncmissingr, 8))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
       ok = NF90_PUT_ATT(ncid, parID, 'missing_value', REAL(ncmissingr, 8))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
    ELSE
       ok = NF90_PUT_ATT(ncid, parID, '_FillValue', REAL(ncmissingr, 4))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
       ok = NF90_PUT_ATT(ncid, parID, 'missing_value', REAL(ncmissingr, 4))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//        &
            'output file. (INTERFACE define_ovar)')
    END IF

  END SUBROUTINE define_output_parameter_r1
  !=============================================================================
  SUBROUTINE define_output_parameter_r2(ncid, parID, pname, punits, longname,  &
       writepatch, othdimID, dimswitch, xID,  &
       yID, zID, landID, patchID, restart)
    ! Subroutine for defining a real valued 2D parameter (time invariant)
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: othdimID ! ID of parameter's second dimension
    INTEGER, INTENT(IN) :: xID, yID, zID, landID, patchID ! netcdf
    ! dimension IDs
    INTEGER, INTENT(OUT) :: parID ! variable's netcdf ID
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    LOGICAL,INTENT(IN),OPTIONAL :: restart ! are we writing to a restart file?
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: punits ! variable units
    CHARACTER(LEN=*), INTENT(IN) :: longname ! full variable name
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF((output%grid(1:3) == 'mas' .OR.                                         &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.               &
         output%grid(1:3) == 'ALM') .AND. .NOT. PRESENT(restart)) THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND. (.NOT. output%grid(1:3)         &
            == 'ALM') .AND.(dimswitch/='surftype')) THEN
          WRITE(logn, *) 'Writing '//pname//                                   &
               ' to output file using mask grid with patch-specific info'
          ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/xID, yID, patchID,      &
               othdimID/),parID)
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok, 'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r2)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(dimswitch == 'soil' .OR. dimswitch == 'r2soil') THEN
             IF(.NOT. ASSOCIATED(otmp4xyps))                                   &
                  ALLOCATE(otmp4xyps(xdimsize, ydimsize, max_vegpatches, ms))
          ELSE IF(dimswitch == 'plantcarbon') THEN
             IF(.NOT. ASSOCIATED(otmp4xyppc))                                  &
                  ALLOCATE(otmp4xyppc(xdimsize, ydimsize, max_vegpatches, ncp))
          ELSE IF(dimswitch == 'soilcarbon') THEN
             IF(.NOT. ASSOCIATED(otmp4xypsc))                                  &
                  ALLOCATE(otmp4xypsc(xdimsize, ydimsize, max_vegpatches, ncs))
          ELSE IF(dimswitch == 'radiation') THEN
             IF(.NOT. ASSOCIATED(otmp4xypr))                                   &
                  ALLOCATE(otmp4xypr(xdimsize, ydimsize, max_vegpatches, nrb))
          END IF
       ELSE ! only grid point values, no patch-specific info
          WRITE(logn, *) 'Writing '//pname//' to output file using mask grid'
          ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/xID, yID, othdimID/)    &
               , parID)
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(dimswitch == 'soil' .OR. dimswitch == 'r2soil') THEN
             IF(.NOT. ASSOCIATED(otmp3xys)) ALLOCATE(otmp3xys(xdimsize,        &
                  ydimsize, ms))
          ELSE IF(dimswitch == 'plantcarbon') THEN
             IF(.NOT. ASSOCIATED(otmp3xypc))                                   &
                  ALLOCATE(otmp3xypc(xdimsize, ydimsize, ncp))
          ELSE IF(dimswitch == 'soilcarbon') THEN
             IF(.NOT. ASSOCIATED(otmp3xysc))                                   &
                  ALLOCATE(otmp3xysc(xdimsize, ydimsize, ncs))
          ELSE IF(dimswitch == 'radiation') THEN
             IF(.NOT. ASSOCIATED(otmp3xyr))                                    &
                  ALLOCATE(otmp3xyr(xdimsize, ydimsize, nrb))
          ELSE IF(dimswitch == 'surftype') THEN
             IF(.NOT. ASSOCIATED(otmp3xysf)) ALLOCATE(otmp3xysf(xdimsize,      &
                  ydimsize, 4))
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok, 'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r2)')
       END IF
    ELSE IF(output%grid(1:3) == 'lan' .OR. (output%grid(1:3) == 'def'          &
         .AND. metGrid=='land') .OR. PRESENT(restart)) THEN
       ! Should patch-specific info be written for this variable?
       ! If this variable has been requested by user with patch-specific info
       ! (writepatch) OR all have been (output%patch) AND we're NOT writing
       ! a restart file (which uses a different technique to store patch info):
       IF((writepatch .OR. output%patch) .AND. (dimswitch /= 'surftype')       &
            .AND. .NOT. PRESENT(restart)) THEN
          WRITE(logn, *) 'Writing '//pname//                                   &
               ' to output file using land grid with patch-specific info'
          ! Define parameter as double precision if required:
          IF(dimswitch(1:2) == 'r2') THEN
             ok = NF90_DEF_VAR(ncid, pname, NF90_DOUBLE, (/landID, patchID,    &
                  othdimID/), parID)
          ELSE
             ok = NF90_DEF_VAR(ncid, pname, NF90_FLOAT, (/landID, patchID,     &
                  othdimID/), parID)
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort                                  &
               (ok, 'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r2)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dim:
          IF(dimswitch == 'soil' .OR. dimswitch == 'r2soil') THEN
             IF(.NOT. ASSOCIATED(otmp3lps)) ALLOCATE(otmp3lps(mland,           &
                  max_vegpatches, ms))
          ELSE IF(dimswitch == 'plantcarbon') THEN
             IF(.NOT. ASSOCIATED(otmp3lppc))                                   &
                  ALLOCATE(otmp3lppc(mland, max_vegpatches, ncp))
          ELSE IF(dimswitch == 'soilcarbon') THEN
             IF(.NOT. ASSOCIATED(otmp3lpsc))                                   &
                  ALLOCATE(otmp3lpsc(mland, max_vegpatches, ncs))
          ELSE IF(dimswitch == 'radiation') THEN
             IF(.NOT. ASSOCIATED(otmp3lpr))                                    &
                  ALLOCATE(otmp3lpr(mland, max_vegpatches, nrb))
          ELSE IF(dimswitch == 'snow') THEN
             IF(.NOT. ASSOCIATED(otmp3lpsn))                                   &
                  ALLOCATE(otmp3lpsn(mland, max_vegpatches, msn))
          END IF
       ELSE ! variable has no explicit patch dimension (incl. restart file)
          ! Restart file definitions will be directed to this part of interface.
          ! If not writing a restart file, report variable writing to log file:
          IF(.NOT.PRESENT(restart)) WRITE(logn,*) 'Writing '//pname// &
               ' to output file using land grid'
          ! Define parameter as double precision if required for restart file:
          IF(dimswitch(1:2)=='r2') THEN
             ok=NF90_DEF_VAR(ncid,pname,NF90_DOUBLE,(/landID,othdimID/),parID)
          ELSE
             ok=NF90_DEF_VAR(ncid,pname,NF90_FLOAT,(/landID,othdimID/),parID)
          END IF
          IF (ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error defining '//pname//' variable in output file. '// &
               '(SUBROUTINE define_output_parameter_r2)')
          ! If not already allocated, allocate a temporary storage variable
          ! of this dimension structure:
          IF(dimswitch=='soil'.OR.dimswitch=='r2soil') THEN
             IF(.NOT.ASSOCIATED(otmp2ls)) ALLOCATE(otmp2ls(mland,ms))
          ELSE IF(dimswitch=='plantcarbon') THEN
             IF(.NOT.ASSOCIATED(otmp2lpc)) ALLOCATE(otmp2lpc(mland,ncp))
          ELSE IF(dimswitch=='soilcarbon') THEN
             IF(.NOT.ASSOCIATED(otmp2lsc)) ALLOCATE(otmp2lsc(mland,ncs))
          ELSE IF(dimswitch=='radiation') THEN
             IF(.NOT.ASSOCIATED(otmp2lr)) ALLOCATE(otmp2lr(mland,nrb))
          ELSE IF(dimswitch=='snow') THEN
             IF(.NOT.ASSOCIATED(otmp2lsn)) ALLOCATE(otmp2lsn(mland,msn))
          ELSE IF(dimswitch=='surftype') THEN
             IF(.NOT.ASSOCIATED(otmp2lsf)) ALLOCATE(otmp2lsf(mland,4))
          END IF
       END IF
    ELSE
       CALL abort('Unknown grid specification '//                              &
            '(SUBROUTINE define_output_parameter_r2)')
    END IF
    ! Define variable units:
    ok = NF90_PUT_ATT(ncid ,parID, 'units', punits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//pname//' variable attributes in '//             &
         'output file. (SUBROUTINE define_output_parameter_r2)')
    ! Define long name:
    ok = NF90_PUT_ATT(ncid, parID, 'long_name', longname)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining '//pname//' variable attributes in '//             &
         'output file. (SUBROUTINE define_output_parameter_r2)')
    ! Define missing/fill values:
    IF(dimswitch(1:1) == 'i') THEN
       ok = NF90_PUT_ATT(ncid, parID, '_FillValue', ncmissingi)
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
       ok = NF90_PUT_ATT(ncid, parID, 'missing_value', ncmissingi)
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
    ELSE IF(dimswitch(1:2) == 'r2') THEN
       ok = NF90_PUT_ATT(ncid, parID, '_FillValue', REAL(ncmissingr, 8))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
       ok = NF90_PUT_ATT(ncid, parID, 'missing_value', REAL(ncmissingr, 8))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
    ELSE
       ok = NF90_PUT_ATT(ncid, parID, '_FillValue', REAL(ncmissingr, 4))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//        &
            'output file. (INTERFACE define_ovar)')
       ok = NF90_PUT_ATT(ncid, parID, 'missing_value', REAL(ncmissingr, 4))
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
            (ok, 'Error defining '//pname//' variable attributes in '//          &
            'output file. (INTERFACE define_ovar)')
    END IF

  END SUBROUTINE define_output_parameter_r2
  !=============================================================================
  SUBROUTINE write_output_variable_r1(ktau, ncid, varID, vname, var_r1,        &
       vrange, writepatch, dimswitch, met)
    ! Subroutine for writing a real valued 1D variable
    INTEGER, INTENT(IN) :: ktau ! current time step #
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: varID ! variable's netcdf ID
    REAL(KIND=4), DIMENSION(:), INTENT(IN) :: var_r1 ! variable values
    REAL, DIMENSION(2), INTENT(IN) :: vrange ! max and min for variable
    ! error checking
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter
    TYPE(met_type), INTENT(IN) :: met  ! met data

    INTEGER :: i,j ! do loop counter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF(output%grid(1:3) == 'mas' .OR.                                          &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.                &
         output%grid(1:3) == 'ALM') THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND. (.NOT. output%grid(1:3)         &
            == 'ALM')) THEN
          DO i = 1, mland ! over all land grid points
             ! First write data for active patches:
             otmp4xypt(land_x(i), land_y(i), 1:landpt(i)%nap, 1)               &
                  = var_r1(landpt(i)%cstart:landpt(i)%cend)
             ! Then write data for inactive patches (if any) as dummy value:
             IF(landpt(i)%nap < max_vegpatches) otmp4xypt(land_x(i),           &
                  land_y(i), (landpt(i)%nap + 1):max_vegpatches, 1) = ncmissingr
             IF(check%ranges) THEN  ! Check ranges:
                DO j = 1, landpt(i)%nap ! only check active patches
                   IF((otmp4xypt(land_x(i), land_y(i), j, 1) < vrange(1)) .OR. &
                        (otmp4xypt(land_x(i), land_y(i), j, 1) > vrange(2)))   &
                        CALL range_abort(vname//' is out of specified ranges!',&
                        ktau, met, otmp4xypt(land_x(i), land_y(i), j, 1),      &
                        vrange, i, land_x(i),land_y(i))
                END DO
             END IF
          END DO
          ! Fill non-land points with dummy value:
          DO j = 1, max_vegpatches
             WHERE(mask /= 1) otmp4xypt(:, :, j, 1) = ncmissingr ! not land
          END DO
          ! write data to file
          ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xypt(:, :, :, 1), 4),       &
               start = (/1, 1, 1, ktau/),                         &
               count = (/xdimsize, ydimsize, max_vegpatches, 1/))
       ELSE ! only grid point values, no patch-specific info
          ! If this is an ALMA 4D surface variable
          ! AND the user has forced the grid type as ALMA:
          IF(dimswitch == 'ALMA' .AND. output%grid(1:3) == 'ALM') THEN
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (area weighted average across all
                ! patches):
                otmp4xyzt(land_x(i), land_y(i), 1, 1) =                           &
                     SUM(var_r1(landpt(i)%cstart: &
                     landpt(i)%cend) * patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                IF(check%ranges) THEN  ! Check ranges:
                   IF((otmp4xyzt(land_x(i), land_y(i), 1, 1) < vrange(1)) .OR.      &
                        (otmp4xyzt(land_x(i), land_y(i), 1, 1) > vrange(2)))          &
                        CALL range_abort(vname//' is out of specified ranges!',       &
                        ktau, met, otmp4xyzt(land_x(i), land_y(i), 1, 1), vrange, i, &
                        land_x(i), land_y(i))
                END IF
             END DO
             ! Fill non-land points with dummy value:
             WHERE(mask /= 1) otmp4xyzt(:, :, 1, 1) = ncmissingr ! not land
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xyzt, 4),                   &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/xdimsize, ydimsize, 1, 1/)) ! write data to file
          ELSE ! normal x-y-t mask grid
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (area weighted average across all
                ! patches):
                otmp3xyt(land_x(i), land_y(i), 1) = SUM(var_r1(landpt(i)%cstart:   &
                     landpt(i)%cend) * patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                IF(check%ranges) THEN  ! Check ranges:
                   IF((otmp3xyt(land_x(i), land_y(i), 1) < vrange(1)) .OR.          &
                        (otmp3xyt(land_x(i), land_y(i), 1) > vrange(2)))              &
                        CALL range_abort(vname//' is out of specified ranges!',       &
                        ktau, met, otmp3xyt(land_x(i), land_y(i), 1), vrange, i, &
                        land_x(i), land_y(i))
                END IF
             END DO
             ! Fill non-land points with dummy value:
             WHERE(mask /= 1) otmp3xyt(:, :, 1) = ncmissingr ! not land
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3xyt, 4),                    &
                  start = (/1,1,ktau/),                                &
                  count = (/xdimsize, ydimsize, 1/)) ! write data to file
          END IF
       END IF
    ELSE IF(output%grid(1:3) == 'lan'                                          &
         .OR. (output%grid(1:3) == 'def' .AND. metGrid == 'land')) THEN
       ! Should patch-specific info be written for this variable?
       IF(writepatch .OR. output%patch) THEN
          DO i = 1, mland ! over all land grid points
             ! First write data for active patches:
             otmp3lpt(i, 1:landpt(i)%nap, 1) =                                   &
                  var_r1(landpt(i)%cstart:landpt(i)%cend)
             ! Then write data for inactive patches as dummy value:
             IF(landpt(i)%nap < max_vegpatches)                                  &
                  otmp3lpt(i, (landpt(i)%nap + 1):max_vegpatches, 1) = ncmissingr
             IF(check%ranges) THEN  ! Check ranges for active patches:
                DO j = 1,landpt(i)%nap
                   IF((otmp3lpt(i, j, 1) < vrange(1)) .OR.                          &
                        (otmp3lpt(i, j, 1) > vrange(2)))                            &
                        CALL range_abort(vname//' is out of specified ranges!',     &
                        ktau, met, otmp3lpt(i, j, 1), vrange, i)
                END DO
             END IF
          END DO
          ! write data to file
          ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3lpt(:, :, 1), 4),             &
               start = (/1, 1, ktau/), count = (/mland, max_vegpatches, 1/))
       ELSE ! only grid point values, no patch-specific info
          DO i = 1, mland ! over all land grid points
             ! Write to temporary variable (area weighted average across all
             ! patches):
             otmp2lt(i, 1) = SUM(var_r1(landpt(i)%cstart:                         &
                  landpt(i)%cend) * patch(landpt(i)%cstart:landpt(i)%cend)%frac)
             IF(check%ranges) THEN  ! Check ranges:
                IF((otmp2lt(i, 1) < vrange(1)) .OR.                                &
                     (otmp2lt(i, 1) > vrange(2)))                                  &
                     CALL range_abort(vname//' is out of specified ranges!',       &
                     ktau, met, otmp2lt(i, 1), vrange, i)
             END IF
          END DO
          ok = NF90_PUT_VAR(ncid, varID, REAL(otmp2lt, 4),                       &
               start = (/1, ktau/), count = (/mland, 1/)) ! write data to file
       END IF
    ELSE
       CALL abort('Unknown grid specification '//                               &
            '(SUBROUTINE write_output_variable_r1)')
    END IF
    ! Check writing was successful:
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//vname//           &
         ' variable to output file (SUBROUTINE write_output_variable_r1)')

  END SUBROUTINE write_output_variable_r1
  !=============================================================================
  SUBROUTINE write_output_variable_r2(ktau, ncid, varID, vname, var_r2,        &
       vrange, writepatch, dimswitch, met)
    ! Subroutine for writing a real valued 2D variable
    INTEGER, INTENT(IN) :: ktau ! current time step #
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: varID ! variable's netcdf ID
    REAL(KIND=4), DIMENSION(:, :), INTENT(IN) :: var_r2 ! variable values
    REAL, DIMENSION(2), INTENT(IN) :: vrange ! max and min for variable
    ! error checking
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    CHARACTER(LEN=*), INTENT(IN) :: vname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter
    TYPE(met_type), INTENT(IN) :: met  ! met data

    INTEGER :: i, j, k ! do loop counter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF(output%grid(1:3) == 'mas' .OR.                                          &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.                 &
         output%grid(1:3) == 'ALM') THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND. (.NOT. output%grid(1:3)          &
            == 'ALM')) THEN
          ! Decide what the second dimension of this variable is:
          IF(dimswitch == 'soil') THEN ! other dim is soil
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp5xypst(land_x(i), land_y(i), 1:landpt(i)%nap, :, 1)           &
                     = var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches (if any) as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp5xypst(land_x(i),          &
                     land_y(i), (landpt(i)%nap+1):max_vegpatches,:,1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, ms
                         IF((otmp5xypst(land_x(i), land_y(i), j, k, 1) < vrange(1))   &
                              .OR. (otmp5xypst(land_x(i), land_y(i), j, k, 1) > vrange(2)))&
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp5xypst(land_x(i), land_y(i), j, k, 1),   &
                              vrange, i, land_x(i), land_y(i))
                      END DO
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, ms
                   WHERE(mask /=1 ) otmp5xypst(:, :, j, k, 1) = ncmissingr ! not land
                END DO
             END DO
             ! Write data to file:
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp5xypst(:, :, :, :,1), 4),    &
                  start = (/1, 1, 1, 1, ktau/),                      &
                  count = (/xdimsize, ydimsize, max_vegpatches, ms, 1/))
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp5xypsnt(land_x(i), land_y(i), 1:landpt(i)%nap, :, 1)          &
                     = var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp5xypsnt(land_x(i),         &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, msn
                         IF((otmp5xypsnt(land_x(i), land_y(i), j, k, 1) < vrange(1))  &
                              .OR. (otmp5xypsnt(land_x(i),land_y(i),j,k,1)>vrange(2)))  &
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp5xypsnt(land_x(i), land_y(i), j, k, 1),  &
                              vrange, i, land_x(i), land_y(i))
                      END DO
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, msn
                   ! not land
                   WHERE(mask /= 1) otmp5xypsnt(:, :, j, k, 1) = ncmissingr
                END DO
             END DO
             ! Write data to file:
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp5xypsnt(:, :, :, :, 1), 4),  &
                  start = (/1, 1, 1, 1, ktau/),                      &
                  count = (/xdimsize, ydimsize, max_vegpatches, msn, 1/))
          ELSE IF(dimswitch=='radiation') THEN ! other dim is radiation bands
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp5xyprt(land_x(i), land_y(i), 1:landpt(i)%nap, :, 1)           &
                     = var_r2(landpt(i)%cstart:landpt(i)%cend,:)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp5xyprt(land_x(i),         &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, nrb
                         IF((otmp5xyprt(land_x(i), land_y(i), j, k, 1) < vrange(1))   &
                              .OR. (otmp5xyprt(land_x(i), land_y(i), j, k, 1) >        &
                              vrange(2)))                                              &
                              CALL range_abort(vname//' is out of specified ranges!',  &
                              ktau, met, otmp5xyprt(land_x(i), land_y(i), j, k, 1),    &
                              vrange, i, land_x(i), land_y(i))
                      END DO
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, nrb
                   ! not land
                   WHERE(mask /= 1) otmp5xyprt(:, :, j, k, 1) = ncmissingr
                END DO
             END DO
             ! Write data to file:
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp5xyprt(:, :, :, :, 1), 4),   &
                  start = (/1, 1, 1, 1, ktau/),                      &
                  count = (/xdimsize, ydimsize, max_vegpatches, nrb, 1/))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp5xyppct(land_x(i), land_y(i), 1:landpt(i)%nap, :, 1)          &
                     = var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches (if any) as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp5xyppct(land_x(i),        &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, ncp
                         IF((otmp5xyppct(land_x(i), land_y(i), j, k, 1) < vrange(1))  &
                              .OR. (otmp5xyppct(land_x(i), land_y(i), j, k, 1) >        &
                              vrange(2)))                                               &
                              CALL range_abort(vname//' is out of specified ranges!',   &
                              ktau, met, otmp5xyppct(land_x(i), land_y(i), j, k, 1),    &
                              vrange, i, land_x(i), land_y(i))
                      END DO
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, ncp
                   ! not land
                   WHERE(mask /= 1) otmp5xyppct(:, :, j, k, 1) = ncmissingr
                END DO
             END DO
             ! Write data to file:
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp5xyppct(:, :, :, :, 1), 4),  &
                  start = (/1, 1, 1, 1, ktau/),                      &
                  count = (/xdimsize, ydimsize, max_vegpatches, ncp, 1/))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp5xypsct(land_x(i), land_y(i), 1:landpt(i)%nap, :, 1)          &
                     = var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp5xypsct(land_x(i),         &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, ncs
                         IF((otmp5xypsct(land_x(i), land_y(i), j, k, 1) < vrange(1))  &
                              .OR. (otmp5xypsct(land_x(i), land_y(i), j, k, 1) >       &
                              vrange(2)))                                              &
                              CALL range_abort(vname//' is out of specified ranges!',  &
                              ktau, met, otmp5xypsct(land_x(i), land_y(i), j, k, 1),   &
                              vrange, i, land_x(i), land_y(i))
                      END DO
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, ncs
                   ! not land
                   WHERE(mask /= 1) otmp5xypsct(:, :, j, k, 1) = ncmissingr
                END DO
             END DO
             ! Write data to file:
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp5xypsct(:, :, :, :, 1), 4),  &
                  start = (/1, 1, 1, 1, ktau/),                      &
                  count = (/xdimsize, ydimsize, max_vegpatches, ncs, 1/))
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in INTERFACE write_ovar')
          END IF
       ELSE ! only grid point values, no patch-specific info
          ! Decide what the second dimension of this variable is:
          IF(dimswitch == 'soil') THEN ! other dim is soil
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by
                ! fraction):
                DO j = 1, ms
                   otmp4xyst(land_x(i), land_y(i), j, 1) = SUM(                    &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) *  &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, ms
                      IF((otmp4xyst(land_x(i), land_y(i), j, 1) < vrange(1)) .OR.    &
                           (otmp4xyst(land_x(i), land_y(i), j, 1) > vrange(2)))      &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp4xyst(land_x(i), land_y(i), j, 1),         &
                           vrange, i, land_x(i), land_y(i))
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, ms
                WHERE(mask /= 1) otmp4xyst(:, :, j, 1) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xyst, 4),                   &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/xdimsize, ydimsize, ms, 1/)) ! write data to file
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by
                ! fraction):
                DO j = 1, msn
                   otmp4xysnt(land_x(i), land_y(i), j, 1) = SUM(                   &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) *  &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, msn
                      IF((otmp4xysnt(land_x(i), land_y(i), j, 1) < vrange(1)) .OR.   &
                           (otmp4xysnt(land_x(i), land_y(i), j, 1) > vrange(2)))     &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp4xysnt(land_x(i), land_y(i), j, 1),        &
                           vrange, i, land_x(i), land_y(i))
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, msn
                WHERE(mask /= 1) otmp4xysnt(:, :, j, 1) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xysnt, 4),                  &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/xdimsize, ydimsize, msn, 1/)) ! write data to file
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by
                ! fraction):
                DO j = 1, nrb
                   otmp4xyrt(land_x(i), land_y(i), j, 1) = SUM(                    &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, nrb
                      IF((otmp4xyrt(land_x(i), land_y(i), j, 1) < vrange(1)) .OR.    &
                           (otmp4xyrt(land_x(i), land_y(i), j, 1) > vrange(2)))      &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp4xyrt(land_x(i), land_y(i), j, 1),         &
                           vrange, i, land_x(i), land_y(i))
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, nrb
                WHERE(mask /= 1) otmp4xyrt(:, :, j, 1) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xyrt, 4),                   &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/xdimsize, ydimsize, nrb, 1/)) ! write data to file
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by fraction):
                DO j = 1, ncp
                   otmp4xypct(land_x(i), land_y(i), j, 1) = SUM(                   &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, ncp
                      IF((otmp4xypct(land_x(i), land_y(i), j, 1) < vrange(1)) .OR.   &
                           (otmp4xypct(land_x(i), land_y(i), j, 1) > vrange(2)))     &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp4xypct(land_x(i), land_y(i), j, 1),        &
                           vrange, i, land_x(i), land_y(i))
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, ncp
                WHERE(mask /= 1) otmp4xypct(:, :, j, 1) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xypct, 4),                  &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/xdimsize, ydimsize, ncp, 1/)) ! write data to file
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by fraction):
                DO j = 1, ncs
                   otmp4xysct(land_x(i), land_y(i), j, 1) = SUM(                   &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, ncs
                      IF((otmp4xysct(land_x(i), land_y(i), j, 1) < vrange(1)) .OR.   &
                           (otmp4xysct(land_x(i), land_y(i), j, 1) > vrange(2)))     &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp4xysct(land_x(i), land_y(i), j, 1),        &
                           vrange, i, land_x(i), land_y(i))
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, ncs
                WHERE(mask /= 1) otmp4xysct(:, :, j, 1) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4xysct, 4),                  &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/xdimsize, ydimsize, ncs, 1/)) ! write data to file
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in INTERFACE write_ovar')
          END IF
       END IF
    ELSE IF(output%grid(1:3) == 'lan'                                          &
         .OR.(output%grid(1:3) == 'def' .AND. metGrid == 'land')) THEN
       ! Should patch-specific info be written for this variable
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND. (.NOT. output%grid(1:3)          &
            == 'ALM')) THEN
          ! Decide what the second dimension of this variable is:
          IF(dimswitch == 'soil') THEN ! other dim is soil
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp4lpst(i, 1:landpt(i)%nap, :, 1)                               &
                     = var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches (if any) as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4lpst(i,                  &
                     (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, ms
                         IF((otmp4lpst(i, j, k, 1) < vrange(1)) .OR.                  &
                              (otmp4lpst(i, j, k, 1) > vrange(2)))                    &
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp4lpst(i, j, k, 1), vrange, i)
                      END DO
                   END DO
                END IF
             END DO
             ! Write data to file:
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4lpst(:, :, :, 1), 4),       &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/mland, max_vegpatches, ms, 1/))
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp4lpsnt(i, 1:landpt(i)%nap, :, 1) =                           &
                     var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4lpsnt(i,                 &
                     (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, msn
                         IF((otmp4lpsnt(i, j, k, 1) < vrange(1)) .OR.                 &
                              (otmp4lpsnt(i, j, k, 1) > vrange(2)))                   &
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp4lpsnt(i, j, k, 1), vrange, i)
                      END DO
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4lpsnt(:, :, :, 1), 4),      &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/mland, max_vegpatches, msn, 1/))
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp4lprt(i, 1:landpt(i)%nap, :, 1) =                             &
                     var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4lprt(i,                  &
                     (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, nrb
                         IF((otmp4lprt(i, j, k, 1) < vrange(1)) .OR.                  &
                              (otmp4lprt(i, j, k, 1) > vrange(2)))                    &
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp4lprt(i, j, k, 1), vrange, i)
                      END DO
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4lprt(:, :, :, 1), 4),       &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/mland, max_vegpatches, nrb, 1/))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp4lppct(i, 1:landpt(i)%nap, :, 1) =                            &
                     var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4lppct(i,                 &
                     (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, ncp
                         IF((otmp4lppct(i, j, k, 1) < vrange(1)) .OR.                 &
                              (otmp4lppct(i, j, k, 1) > vrange(2)))                   &
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp4lppct(i, j, k, 1), vrange, i)
                      END DO
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4lppct(:, :, :, 1), 4),      &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/mland, max_vegpatches, ncp, 1/))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp4lpsct(i, 1:landpt(i)%nap, :, 1) =                            &
                     var_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4lpsct(i,                 &
                     (landpt(i)%nap + 1):max_vegpatches, :, 1) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges for active patches:
                   DO j = 1, landpt(i)%nap
                      DO k = 1, ncs
                         IF((otmp4lpsct(i, j, k, 1) < vrange(1)) .OR.                 &
                              (otmp4lpsct(i, j, k, 1) > vrange(2)))                   &
                              CALL range_abort(vname//' is out of specified ranges!', &
                              ktau, met, otmp4lpsct(i, j, k, 1), vrange, i)
                      END DO
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp4lpsct(:, :, :, 1), 4),      &
                  start = (/1, 1, 1, ktau/),                         &
                  count = (/mland, max_vegpatches, ncs, 1/))
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in INTERFACE write_ovar')
          END IF
       ELSE ! only grid point values, no patch-specific info
          ! Decide what the second dimension of this variable is:
          IF(dimswitch == 'soil') THEN ! other dim is soil
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by
                ! fraction):
                DO j = 1, ms
                   otmp3lst(i, j, 1) = SUM(                                        &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, ms
                      IF((otmp3lst(i, j, 1) < vrange(1)) .OR.                        &
                           (otmp3lst(i, j, 1) > vrange(2)))                          &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp3lst(i, j, 1), vrange, i)
                   END DO
                END IF
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3lst, 4),                    &
                  start = (/1, 1, ktau/),                            &
                  count = (/mland, ms, 1/)) ! write data to file
          ELSE IF(dimswitch == 'snow') THEN ! other dim is snow
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by
                ! fraction):
                DO j = 1, msn
                   otmp3lsnt(i, j, 1) = SUM(                                       &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, msn
                      IF((otmp3lsnt(i, j, 1) < vrange(1)) .OR.                       &
                           (otmp3lsnt(i, j, 1) > vrange(2)))                         &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp3lsnt(i, j, 1), vrange, i)
                   END DO
                END IF
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3lsnt, 4),                   &
                  start = (/1, 1, ktau/),                            &
                  count = (/mland, msn, 1/)) ! write data to file
          ELSE IF(dimswitch == 'radiation') THEN ! other dim is radiation bands
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by fraction):
                DO j = 1, nrb
                   otmp3lrt(i, j, 1) = SUM(                                        &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, nrb
                      IF((otmp3lrt(i, j, 1) < vrange(1)) .OR.                        &
                           (otmp3lrt(i, j, 1) > vrange(2)))                          &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp3lrt(i, j, 1), vrange, i)
                   END DO
                END IF
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3lrt, 4),                    &
                  start = (/1, 1, ktau/),                            &
                  count = (/mland, nrb, 1/)) ! write data to file
          ELSE IF(dimswitch == 'plantcarbon') THEN ! other dim is plant carbon
             ! pools
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by fraction):
                DO j = 1, ncp
                   otmp3lpct(i, j, 1) = SUM(                                       &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, ncp
                      IF((otmp3lpct(i, j, 1) < vrange(1)) .OR.                       &
                           (otmp3lpct(i, j, 1) > vrange(2)))                         &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp3lpct(i, j, 1), vrange, i)
                   END DO
                END IF
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3lpct, 4),                   &
                  start = (/1, 1, ktau/),                            &
                  count = (/mland, ncp, 1/)) ! write data to file
          ELSE IF(dimswitch == 'soilcarbon') THEN ! other dim is soil carbon pools
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (sum over patches & weight by fraction):
                DO j = 1, ncs
                   otmp3lsct(i, j, 1) = SUM(                                       &
                        var_r2(landpt(i)%cstart:landpt(i)%cend, j) * &
                        patch(landpt(i)%cstart:landpt(i)%cend)%frac)
                END DO
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, ncs
                      IF((otmp3lsct(i, j, 1) < vrange(1)) .OR.                       &
                           (otmp3lsct(i, j, 1) > vrange(2)))                         &
                           CALL range_abort(vname//' is out of specified ranges!',   &
                           ktau, met, otmp3lsct(i, j, 1), vrange, i)
                   END DO
                END IF
             END DO
             ok = NF90_PUT_VAR(ncid, varID, REAL(otmp3lsct, 4),                   &
                  start = (/1, 1, ktau/),                            &
                  count = (/mland, ncs, 1/)) ! write data to file
          ELSE
             CALL abort('Variable '//vname//                                      &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in SUBROUTINE write_output_variable_r2')
          END IF
       END IF ! patch info or no patch info
    ELSE
       CALL abort('Unknown grid specification '//                               &
            '(SUBROUTINE write_output_variable_r2)')
    END IF ! grid type

    ! Check writing was successful:
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//vname//           &
         ' variable to output file (SUBROUTINE write_output_variable_r2)')

  END SUBROUTINE write_output_variable_r2
  !=============================================================================
  SUBROUTINE write_output_parameter_r1(ncid, parID, pname, par_r1,             &
       prange, writepatch, dimswitch, restart)
    ! Subroutine for writing a real valued 1D parameter (time invariant)
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: parID ! variable's netcdf ID
    REAL(KIND=4), DIMENSION(:), INTENT(IN) :: par_r1 ! variable values
    REAL, DIMENSION(2), INTENT(IN) :: prange ! max and min for variable
    ! error checking
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    LOGICAL, INTENT(IN), OPTIONAL :: restart ! are we writing to a restart file?
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    INTEGER :: i, j ! do loop counter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF((output%grid(1:3) == 'mas' .OR.                                         &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.                &
         output%grid(1:3) == 'ALM') .AND. .NOT. PRESENT(restart)) THEN
       ! Should patch-specific info be written for this parameter
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND. (.NOT. output%grid(1:3)         &
            == 'ALM')) THEN
          DO i = 1, mland ! over all land grid points
             ! First write data for active patches:
             otmp3xyp(land_x(i), land_y(i), 1:landpt(i)%nap)                   &
                  = par_r1(landpt(i)%cstart:landpt(i)%cend)
             ! Then write data for inactive patches as dummy value:
             IF(dimswitch(1:1) == 'r') THEN
                IF(landpt(i)%nap<max_vegpatches) otmp3xyp(land_x(i), land_y(i),&
                     (landpt(i)%nap + 1):max_vegpatches) = ncmissingr
             ELSE IF(dimswitch(1:1) == 'i') THEN
                IF(landpt(i)%nap < max_vegpatches) otmp3xyp(land_x(i),         &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches) = ncmissingi
             END IF
             IF(check%ranges) THEN  ! Check ranges over active patches:
                DO j = 1, landpt(i)%nap
                   IF((otmp3xyp(land_x(i), land_y(i), j) < prange(1)) .OR.     &
                        (otmp3xyp(land_x(i), land_y(i), j) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                        &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ',i, 'patch #', j
                      WRITE(*, *) 'Value: ', otmp3xyp(land_x(i), land_y(i), j)
                      CALL abort('Aborting.')
                   END IF
                END DO
             END IF
          END DO
          ! Write data to file:
          IF(dimswitch(1:1) == 'r') THEN
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                WHERE(mask /= 1) otmp3xyp(:, :, j) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3xyp(:, :, :), 4),        &
                  start = (/1, 1, 1/),                            &
                  count = (/xdimsize, ydimsize, max_vegpatches/))
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                WHERE(mask /= 1) otmp3xyp(:, :, j) = ncmissingi ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, INT(otmp3xyp(:, :, :)),            &
                  start = (/1, 1, 1/),                            &
                  count = (/xdimsize, ydimsize, max_vegpatches/))
          END IF
       ELSE ! only grid point values, no patch-specific info
          DO i = 1, mland ! over all land grid points
             ! Write to temporary variable. Use dominant patch info only,
             ! as averaging parameters over patches wouldn't make nec sense:
             otmp2xy(land_x(i), land_y(i)) = par_r1(landpt(i)%cstart)
             IF(check%ranges) THEN  ! Check ranges:
                IF((otmp2xy(land_x(i), land_y(i)) < prange(1)) .OR.            &
                     (otmp2xy(land_x(i), land_y(i)) > prange(2))) THEN
                   WRITE(*, *) 'Parameter '//pname//                           &
                        ' is set at a value out of specified ranges!'
                   WRITE(*, *) 'Land point # ',i
                   WRITE(*, *) 'Value: ', otmp2xy(land_x(i),land_y(i))
                   CALL abort('Aborting.')
                END IF
             END IF
          END DO
          ! Write data to file:
          IF(dimswitch(1:1) == 'r') THEN
             ! Fill non-land points with dummy value:
             WHERE(mask /= 1) otmp2xy(:, :) = ncmissingr ! not land
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2xy, 4),                  &
                  start = (/1, 1/),                               &
                  count = (/xdimsize, ydimsize/)) ! write data to file
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ! Fill non-land points with dummy value:
             WHERE(mask /= 1) otmp2xy(:, :) = ncmissingi ! not land
             ok = NF90_PUT_VAR(ncid, parID, INT(otmp2xy), start = (/1, 1/),    &
                  count = (/xdimsize, ydimsize/)) ! write data to file
          END IF
       END IF
    ELSE IF(output%grid(1:3) == 'lan'                                          &
         .OR. (output%grid(1:3) == 'def' .AND. metGrid == 'land')           &
         .OR. PRESENT(restart)) THEN
       ! Is patch-specific info written for this variable?
       ! If this variable has been requested by user with patch-specific info
       ! (writepatch) OR all have been (output%patch) AND we're NOT writing
       ! a restart file (which uses a different technique to store patch info):
       IF((writepatch .OR. output%patch) .AND. .NOT. PRESENT(restart)) THEN
          DO i = 1, mland ! over all land grid points
             ! First write data for active patches:
             otmp2lp(i, 1:landpt(i)%nap) = par_r1(landpt(i)%cstart:            &
                  landpt(i)%cend)
             ! Then write data for inactive patches as dummy value:
             IF(landpt(i)%nap < max_vegpatches)  THEN
                IF(dimswitch(1:1) == 'r') THEN
                   otmp2lp(i, (landpt(i)%nap + 1):max_vegpatches) = ncmissingr
                ELSE IF(dimswitch(1:1)=='i') THEN
                   otmp2lp(i, (landpt(i)%nap + 1):max_vegpatches) = ncmissingi
                END IF
             END IF
             IF(check%ranges) THEN  ! Check ranges over active patches:
                DO j = 1, landpt(i)%nap
                   IF((otmp2lp(i, j) < prange(1)) .OR.                         &
                        (otmp2lp(i, j) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                        &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ',i, 'patch #', j
                      WRITE(*, *) 'Value: ', otmp2lp(i, j)
                      CALL abort('Aborting.')
                   END IF
                END DO
             END IF
          END DO
          ! Write data to file
          IF(dimswitch(1:1) == 'r') THEN
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2lp(:, :), 4),            &
                  start = (/1, 1/), count = (/mland, max_vegpatches/))
          ELSE IF(dimswitch(1:1) == 'i') THEN
             ok = NF90_PUT_VAR(ncid, parID, INT(otmp2lp(:, :)),                &
                  start = (/1, 1/), count = (/mland, max_vegpatches/))
          END IF
       ELSE ! only grid point values without patch-specific info UNLESS restart
          ! file
          ! All 1D single precision restart file variables are written here.
          IF(PRESENT(restart)) THEN ! If writing restart data:
             ! Write output:
             IF(dimswitch(1:1) == 'r') THEN
                ok = NF90_PUT_VAR(ncid, parID, REAL(par_r1, 4),                &
                     start = (/1/), count = (/mp/)) ! write data to file
             ELSE IF(dimswitch(1:1) == 'i') THEN
                ok = NF90_PUT_VAR(ncid, parID, INT(par_r1),                    &
                     start = (/1/), count = (/mp/)) ! write data to file
             END IF
          ELSE
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (use dominant patch info only!):
                otmp1l(i) = par_r1(landpt(i)%cstart)
                IF(check%ranges) THEN  ! Check ranges:
                   IF((otmp1l(i) < prange(1)) .OR.                             &
                        (otmp1l(i) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                        &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ',i
                      WRITE(*, *) 'Value: ', otmp1l(i)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ! Write output:
             IF(dimswitch(1:1) == 'r') THEN
                ok = NF90_PUT_VAR(ncid, parID, REAL(otmp1l, 4),                &
                     start = (/1/), count = (/mland/)) ! write data to file
             ELSE IF(dimswitch(1:1) == 'i') THEN
                ok = NF90_PUT_VAR(ncid, parID, INT(otmp1l),                    &
                     start = (/1/), count = (/mland/)) ! write data to file
             END IF
          END IF ! If writing restart
       END IF ! If writing with a patch dimension in output file
    ELSE
       CALL abort('Unknown grid specification '//                              &
            '(SUBROUTINE write_output_parameter_r1)')
    END IF  ! mask x-y or land-only grid
    ! Check writing was successful:
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//pname//           &
         ' parameter/variable to output file '//                               &
         '(SUBROUTINE write_output_parameter_r1)')

  END SUBROUTINE write_output_parameter_r1
  !=============================================================================
  SUBROUTINE write_output_parameter_r1d(ncid, parID, pname, par_r1d,           &
       prange, writepatch, dimswitch, restart)
    ! Subroutine for writing a double precision 1D parameter (time invariant)
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: parID ! variable's netcdf ID
    REAL(r_2), DIMENSION(:), INTENT(IN) :: par_r1d ! variable values
    REAL, DIMENSION(2), INTENT(IN) :: prange ! max and min for variable
    ! error checking
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    LOGICAL,INTENT(IN),OPTIONAL :: restart ! are we writing to a restart file?
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    INTEGER :: i, j ! do loop counter
    REAL(r_2), POINTER, DIMENSION(:, :) :: tmpout

    IF(PRESENT(restart)) THEN ! If writing to a a restart file
       ! Write parameter data:
       ok = NF90_PUT_VAR(ncid, parID, par_r1d,                                 &
            start = (/1, 1/), count = (/mp/)) ! write data to file
       ! Check writing was successful:
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//pname//        &
            ' parameter to restart file (SUBROUTINE write_output_parameter_r1d)')
    ELSE ! a 1D double precision time invariant parameter for output file
       ALLOCATE(tmpout(mland, max_vegpatches))
       DO i = 1, mland ! over all land grid points
          ! First write data for active patches:
          tmpout(i, 1:landpt(i)%nap) = par_r1d(landpt(i)%cstart:landpt(i)%cend)
          ! Then write data for inactive patches as dummy value:
          IF(landpt(i)%nap < max_vegpatches)                                   &
               tmpout(i, (landpt(i)%nap + 1):max_vegpatches) =                 &
               REAL(ncmissingr, r_2)
          IF(check%ranges) THEN  ! Check ranges over active patches:
             DO j = 1, landpt(i)%nap
                IF((tmpout(i, j) < prange(1)) .OR.                             &
                     (tmpout(i, j) > prange(2))) THEN
                   WRITE(*, *) 'Parameter '//pname//                           &
                        ' is set at a value out of specified ranges!'
                   WRITE(*, *) 'Land point # ',i, 'patch #', j
                   WRITE(*, *) 'Value: ', tmpout(i, j)
                   CALL abort('Aborting.')
                END IF
             END DO
          END IF
       END DO
       ok = NF90_PUT_VAR(ncid, parID, REAL(tmpout(:, :), 8), start = (/1, 1/), &
            count = (/mland, max_vegpatches/)) ! write data to file
       DEALLOCATE(tmpout)
       ! Check writing was successful:
       IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//pname//        &
            ' variable to output file (SUBROUTINE write_output_parameter_r1d)')
    END IF ! If writing to a a restart file

  END SUBROUTINE write_output_parameter_r1d
  !=============================================================================
  SUBROUTINE write_output_parameter_r2(ncid, parID, pname, par_r2, prange,     &
       writepatch, dimswitch, restart)
    ! Subroutine for writing a real valued 2D parameter (time invariant)
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: parID ! variable's netcdf ID
    REAL(KIND=4), DIMENSION(:, :), INTENT(IN) :: par_r2 ! variable values
    REAL, DIMENSION(2), INTENT(IN) :: prange ! max and min for variable
    ! error checking
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    LOGICAL,INTENT(IN),OPTIONAL :: restart ! are we writing to a restart file?
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    INTEGER :: i, j, k ! do loop counter

    ! First, decide which grid to use. If user has forced grid using output%grid
    ! in the namelist file, use this grid. Else use format of met file.
    IF((output%grid(1:3) == 'mas' .OR.                                         &
         (output%grid(1:3) == 'def' .AND. metGrid == 'mask') .OR.                &
         output%grid(1:3) == 'ALM') .AND. .NOT. PRESENT(restart)) THEN
       ! Should patch-specific info be written for this parameter
       ! (no patches in ALMA format)?
       IF((writepatch .OR. output%patch) .AND. (.NOT. output%grid(1:3)          &
            == 'ALM') .AND. (dimswitch /= 'surftype')) THEN
          ! Check the nature of the parameter's second dimension:
          IF(dimswitch == 'soil') THEN ! i.e. spatial and soil
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (all patches for current grid point):
                ! First write data for active patches:
                otmp4xyps(land_x(i), land_y(i), 1:landpt(i)%nap, :)               &
                     = par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4xyps(land_x(i),           &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp4xyps(land_x(i), land_y(i), j, :) < prange(1)) .OR. &
                           ANY(otmp4xyps(land_x(i), land_y(i), j, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp4xyps(land_x(i), land_y(i), j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j =1, max_vegpatches
                DO k = 1, ms
                   WHERE(mask /= 1) otmp4xyps(:, :, j, k) = ncmissingr ! not land
                END DO
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp4xyps(:, :, :, :), 4),       &
                  start = (/1, 1, 1, 1/),                            &
                  count = (/xdimsize, ydimsize, max_vegpatches, ms/))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! i.e. spatial and plant carbon
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (all patches for current grid point):
                ! First write data for active patches:
                otmp4xyppc(land_x(i), land_y(i), 1:landpt(i)%nap, :)              &
                     = par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4xyppc(land_x(i),          &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp4xyppc(land_x(i), land_y(i), j, :) < prange(1))     &
                           .OR. ANY(otmp4xyppc(land_x(i), land_y(i), j, :) >           &
                           prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp4xyppc(land_x(i), land_y(i), j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, ncp
                   WHERE(mask /= 1) otmp4xyppc(:, :, j, k) = ncmissingr ! not land
                END DO
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp4xyppc(:, :, :, :), 4),      &
                  start = (/1, 1, 1, 1/),                            &
                  count = (/xdimsize, ydimsize, max_vegpatches, ncp/))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! i.e. spatial and soil carbon
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (all patches for current grid point):
                ! First write data for active patches:
                otmp4xypsc(land_x(i), land_y(i), 1:landpt(i)%nap, :)              &
                     = par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4xypsc(land_x(i),          &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp4xypsc(land_x(i), land_y(i), j, :) < prange(1))     &
                           .OR. ANY(otmp4xypsc(land_x(i), land_y(i), j, :) >           &
                           prange(2)))  THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp4xypsc(land_x(i), land_y(i), j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, ncs
                   WHERE(mask /= 1) otmp4xypsc(:, :, j, k) = ncmissingr ! not land
                END DO
             END DO
             ! write data to file:
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp4xypsc(:, :, :, :), 4),      &
                  start = (/1, 1, 1, 1/),                            &
                  count = (/xdimsize, ydimsize, max_vegpatches, ncs/))
          ELSE IF(dimswitch == 'radiation') THEN ! i.e. spatial and soil carbon
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (all patches for current grid point):
                ! First write data for active patches:
                otmp4xypr(land_x(i), land_y(i), 1:landpt(i)%nap, :)               &
                     = par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches) otmp4xypr(land_x(i),           &
                     land_y(i), (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp4xypr(land_x(i), land_y(i), j, :) < prange(1))      &
                           .OR. ANY(otmp4xypr(land_x(i), land_y(i), j, :) >            &
                           prange(2)))  THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp4xypr(land_x(i), land_y(i), j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, max_vegpatches
                DO k = 1, nrb
                   WHERE(mask /= 1) otmp4xypr(:, :, j, k) = ncmissingr ! not land
                END DO
             END DO
             ! write data to file:
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp4xypr(:, :, :, :), 4),       &
                  start = (/1, 1, 1, 1/),                            &
                  count = (/xdimsize, ydimsize, max_vegpatches, nrb/))
          ELSE
             CALL abort('Parameter '//pname//                                     &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in INTERFACE write_ovar')
          END IF
       ELSE ! only grid point values, no patch-specific info
          ! Check the nature of the parameter's second dimension:
          IF(dimswitch == 'soil') THEN ! i.e. spatial and soil
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (use dominant patch info only!):
                otmp3xys(land_x(i), land_y(i), :) = par_r2(landpt(i)%cstart, :)
                IF(check%ranges) THEN  ! Check ranges:
                   IF(ANY(otmp3xys(land_x(i), land_y(i), :) < prange(1)) .OR.       &
                        ANY(otmp3xys(land_x(i), land_y(i), :) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                              &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ', i
                      WRITE(*, *) 'Values: ', otmp3xys(land_x(i), land_y(i), :)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, ms
                WHERE(mask /= 1) otmp3xys(:, :, j) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3xys, 4),                    &
                  start = (/1, 1, 1/),                               &
                  count = (/xdimsize, ydimsize, ms/)) ! write data to file
          ELSE IF(dimswitch == 'plantcarbon') THEN ! i.e. spatial and plant carbon
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (use dominant patch info only!):
                otmp3xypc(land_x(i), land_y(i), :) = par_r2(landpt(i)%cstart, :)
                IF(check%ranges) THEN  ! Check ranges:
                   IF(ANY(otmp3xypc(land_x(i), land_y(i), :) < prange(1)) .OR.      &
                        ANY(otmp3xypc(land_x(i), land_y(i), :) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                              &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ', i
                      WRITE(*, *) 'Values: ', otmp3xypc(land_x(i), land_y(i), :)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, ncp
                WHERE(mask /= 1) otmp3xypc(:, :, j) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3xypc, 4),                   &
                  start = (/1, 1, 1/),                               &
                  count = (/xdimsize, ydimsize, ncp/)) ! write data to file
          ELSE IF(dimswitch == 'soilcarbon') THEN ! i.e. spatial and soil carbon
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (use dominant patch info only!):
                otmp3xysc(land_x(i), land_y(i), :) = par_r2(landpt(i)%cstart, :)
                IF(check%ranges) THEN  ! Check ranges:
                   IF(ANY(otmp3xysc(land_x(i), land_y(i), :) < prange(1)) .OR.      &
                        ANY(otmp3xysc(land_x(i), land_y(i), :) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                              &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ', i
                      WRITE(*, *) 'Values: ', otmp3xysc(land_x(i), land_y(i), :)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, ncs
                WHERE(mask /= 1) otmp3xysc(:, :, j) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3xysc, 4),                   &
                  start = (/1, 1, 1/),                               &
                  count = (/xdimsize, ydimsize, ncs/)) ! write data to file
          ELSE IF(dimswitch == 'radiation') THEN ! i.e. spatial and soil carbon
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (use dominant patch info only!):
                otmp3xyr(land_x(i), land_y(i), :) = par_r2(landpt(i)%cstart, :)
                IF(check%ranges) THEN  ! Check ranges:
                   IF(ANY(otmp3xyr(land_x(i), land_y(i), :) < prange(1)) .OR.       &
                        ANY(otmp3xyr(land_x(i), land_y(i), :) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                              &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ', i
                      WRITE(*, *) 'Values: ', otmp3xyr(land_x(i), land_y(i), :)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, nrb
                WHERE(mask /= 1) otmp3xyr(:, :, j) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3xyr, 4),                    &
                  start = (/1, 1, 1/),                               &
                  count = (/xdimsize, ydimsize, nrb/)) ! write data to file
          ELSE IF(dimswitch == 'surftype') THEN ! i.e. surface fraction
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (surf fraction only has mp dimension):
                otmp3xysf(land_x(i), land_y(i), :) = par_r2(i, :)
                IF(check%ranges) THEN  ! Check ranges:
                   IF(ANY(otmp3xysf(land_x(i), land_y(i), :) < prange(1)) .OR.      &
                        ANY(otmp3xysf(land_x(i), land_y(i), :) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                              &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ', i
                      WRITE(*, *) 'Values: ', otmp3xysf(land_x(i), land_y(i), :)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ! Fill non-land points with dummy value:
             DO j = 1, 4
                WHERE(mask /= 1) otmp3xysf(:, :, j) = ncmissingr ! not land
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3xysf, 4),                   &
                  start = (/1, 1, 1/),                               &
                  count = (/xdimsize, ydimsize, 4/)) ! write data to file
          ELSE
             CALL abort('Parameter '//pname//                                     &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in SUBROUTINE write_output_parameter_r2')
          END IF
       END IF
    ELSE IF(output%grid(1:3) == 'lan'                                          &
         .OR.(output%grid(1:3) == 'def' .AND. metGrid == 'land')            &
         .OR. PRESENT(restart)) THEN
       ! Does this variable have a patch dimension (restart does not)?
       IF((writepatch .OR. output%patch) .AND. (dimswitch /= 'surftype')        &
            .AND. .NOT. PRESENT(restart)) THEN
          ! Check the nature of the parameter's second dimension:
          IF(dimswitch == 'soil') THEN ! i.e. spatial and soil
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp3lps(i, 1:landpt(i)%nap, :) =                                 &
                     par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches)                                &
                     otmp3lps(i, (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges over active patches:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp3lps(i, j, :) < prange(1)) .OR.                     &
                           ANY(otmp3lps(i, j, :) > prange(2)))  THEN
                         WRITE(*, *) 'Parameter '//pname// &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp3lps(i, j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3lps(:, :, :), 4),           &
                  start = (/1, 1, 1/), count = (/mland, max_vegpatches, ms/))
          ELSE IF(dimswitch == 'plantcarbon') THEN ! i.e. spatial and plant carbon
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp3lppc(i, 1:landpt(i)%nap, :) =                                &
                     par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches)                                &
                     otmp3lppc(i, (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges over active patches:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp3lppc(i, j, :) < prange(1)) .OR.                    &
                           ANY(otmp3lppc(i, j, :) > prange(2)))  THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp3lppc(i, j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3lppc(:, :, :), 4),          &
                  start = (/1, 1, 1/), count = (/mland, max_vegpatches, ncp/))
          ELSE IF(dimswitch == 'soilcarbon') THEN ! i.e. spatial and soil carbon
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp3lpsc(i, 1:landpt(i)%nap, :) =                                &
                     par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches)                                &
                     otmp3lpsc(i, (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges over active patches:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp3lpsc(i, j, :) < prange(1)) .OR.                    &
                           ANY(otmp3lpsc(i, j, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp3lpsc(i, j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3lpsc(:, :, :), 4),          &
                  start = (/1, 1, 1/), count = (/mland, max_vegpatches, ncs/))
          ELSE IF(dimswitch == 'radiation') THEN ! i.e. spatial and radiation
             ! bands
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp3lpr(i, 1:landpt(i)%nap, :) =                                 &
                     par_r2(landpt(i)%cstart:landpt(i)%cend,:)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches)                                &
                     otmp3lpr(i, (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges over active patches:
                   DO j = 1,landpt(i)%nap
                      IF(ANY(otmp3lpr(i, j, :) < prange(1)) .OR.                     &
                           ANY(otmp3lpr(i, j, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp3lpr(i, j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3lpr(:, :, :), 4),           &
                  start = (/1, 1, 1/), count = (/mland, max_vegpatches, nrb/))
          ELSE IF(dimswitch == 'snow') THEN ! i.e. spatial and radiation bands
             DO i = 1, mland ! over all land grid points
                ! First write data for active patches:
                otmp3lpsn(i, 1:landpt(i)%nap, :) =                                &
                     par_r2(landpt(i)%cstart:landpt(i)%cend, :)
                ! Then write data for inactive patches as dummy value:
                IF(landpt(i)%nap < max_vegpatches)                                &
                     otmp3lpsn(i, (landpt(i)%nap + 1):max_vegpatches, :) = ncmissingr
                IF(check%ranges) THEN  ! Check ranges over active patches:
                   DO j = 1, landpt(i)%nap
                      IF(ANY(otmp3lpsn(i, j, :) < prange(1)) .OR.                    &
                           ANY(otmp3lpsn(i, j, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                            &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i, 'patch #', j
                         WRITE(*, *) 'Values: ', otmp3lpsn(i, j, :)
                         CALL abort('Aborting.')
                      END IF
                   END DO
                END IF
             END DO
             ! write data to file
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp3lpsn(:, :, :), 4),          &
                  start = (/1, 1, 1/), count = (/mland, max_vegpatches, msn/))
          ELSE
             CALL abort('Parameter '//pname//                                     &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in SUBROUTINE write_output_parameter_r2')
          END IF
       ELSE ! Varaible has no patch dimension
          ! Check the nature of the parameter's second dimension:
          IF(dimswitch=='soil') THEN ! i.e. spatial and soil
             IF(PRESENT(restart)) THEN
                ! Write data to restart file
                ok=NF90_PUT_VAR(ncid,parID,REAL(par_r2,4), &
                     start=(/1,1/),count=(/mp,ms/))
             ELSE
                DO i = 1, mland ! over all land grid points
                   ! Write to temporary variable (use dominant patch info only!):
                   otmp2ls(i,:) = par_r2(landpt(i)%cstart,:)
                   IF(check%ranges) THEN  ! Check ranges:
                      IF(ANY(otmp2ls(i,:)<prange(1)).OR. &
                           ANY(otmp2ls(i,:)>prange(2))) THEN
                         WRITE(*,*) 'Parameter '//pname// &
                              ' is set at a value out of specified ranges!'
                         WRITE(*,*) 'Land point # ',i
                         WRITE(*,*) 'Values: ', otmp2ls(i,:)
                         CALL abort('Aborting.')
                      END IF
                   END IF
                END DO
                ok=NF90_PUT_VAR(ncid,parID,REAL(otmp2ls,4), &
                     start=(/1,1/),count=(/mland,ms/)) ! write data to file
             END IF
          ELSE IF(dimswitch == 'plantcarbon') THEN ! i.e. spatial and plant carbon
             IF(PRESENT(restart)) THEN
                ! Write data to restart file
                ok = NF90_PUT_VAR(ncid, parID, REAL(par_r2, 4),                  &
                     start = (/1, 1/), count = (/mp, ncp/))
             ELSE
                DO i = 1, mland ! over all land grid points
                   ! Write to temporary variable (use dominant patch info only!):
                   otmp2lpc(i, :) = par_r2(landpt(i)%cstart, :)
                   IF(check%ranges) THEN  ! Check ranges:
                      IF(ANY(otmp2lpc(i, :) < prange(1)) .OR.                    &
                           ANY(otmp2lpc(i, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                       &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i
                         WRITE(*, *) 'Values: ', otmp2lpc(i, :)
                         CALL abort('Aborting.')
                      END IF
                   END IF
                END DO
                ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2lpc, 4),                &
                     start = (/1, 1/), count = (/mland, ncp/)) ! write data to file
             END IF
          ELSE IF(dimswitch == 'soilcarbon') THEN ! i.e. spatial and soil carbon
             IF(PRESENT(restart)) THEN
                ! Write data to restart file
                ok = NF90_PUT_VAR(ncid, parID, REAL(par_r2, 4),                  &
                     start = (/1, 1/), count = (/mp, ncs/))
             ELSE
                DO i = 1, mland ! over all land grid points
                   ! Write to temporary variable (use dominant patch info only!):
                   otmp2lsc(i, :) = par_r2(landpt(i)%cstart, :)
                   IF(check%ranges) THEN  ! Check ranges:
                      IF(ANY(otmp2lsc(i, :) < prange(1)) .OR.                    &
                           ANY(otmp2lsc(i, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                       &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i
                         WRITE(*, *) 'Values: ', otmp2lsc(i, :)
                         CALL abort('Aborting.')
                      END IF
                   END IF
                END DO
                ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2lsc, 4),                &
                     start = (/1, 1/), count=(/mland, ncs/)) ! write data to file
             END IF
          ELSE IF(dimswitch == 'radiation') THEN ! i.e. spatial and radiation
             ! bands
             IF(PRESENT(restart)) THEN
                ! Write data to restart file
                ok = NF90_PUT_VAR(ncid, parID, REAL(par_r2, 4),                  &
                     start = (/1, 1/), count = (/mp, nrb/))
             ELSE ! writing to output file
                DO i = 1, mland ! over all land grid points
                   ! Write to temporary variable (use dominant patch info only!):
                   otmp2lr(i, :) = par_r2(landpt(i)%cstart, :)
                   IF(check%ranges) THEN  ! Check ranges:
                      IF(ANY(otmp2lr(i, :) < prange(1)) .OR.                     &
                           ANY(otmp2lr(i, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                       &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i
                         WRITE(*, *) 'Values: ', otmp2lr(i, :)
                         CALL abort('Aborting.')
                      END IF
                   END IF
                END DO
                ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2lr, 4),                 &
                     start = (/1, 1/), count = (/mland, nrb/)) ! write data to file
             END IF
          ELSE IF(dimswitch == 'snow') THEN ! i.e. spatial and radiation bands
             IF(PRESENT(restart)) THEN
                ! Write data to restart file
                ok = NF90_PUT_VAR(ncid, parID, REAL(par_r2, 4),                  &
                     start = (/1, 1/), count = (/mp, msn/))
             ELSE ! writing to output file
                DO i = 1, mland ! over all land grid points
                   ! Write to temporary variable (use dominant patch info only!):
                   otmp2lsn(i, :) = par_r2(landpt(i)%cstart, :)
                   IF(check%ranges) THEN  ! Check ranges:
                      IF(ANY(otmp2lsn(i, :) < prange(1)) .OR.                    &
                           ANY(otmp2lsn(i, :) > prange(2))) THEN
                         WRITE(*, *) 'Parameter '//pname//                       &
                              ' is set at a value out of specified ranges!'
                         WRITE(*, *) 'Land point # ', i
                         WRITE(*, *) 'Values: ', otmp2lsn(i, :)
                         CALL abort('Aborting.')
                      END IF
                   END IF
                END DO
                ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2lsn, 4),                &
                     start = (/1, 1/), count = (/mland, msn/)) ! write data to file
             END IF
          ELSE IF(dimswitch == 'surftype') THEN
             DO i = 1, mland ! over all land grid points
                ! Write to temporary variable (use dominant patch info only!):
                otmp2lsf(i, :) = par_r2(i, :)
                IF(check%ranges) THEN  ! Check ranges:
                   IF(ANY(otmp2lsf(i, :) < prange(1)) .OR.                          &
                        ANY(otmp2lsf(i, :) > prange(2))) THEN
                      WRITE(*, *) 'Parameter '//pname//                              &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ', i
                      WRITE(*, *) 'Values: ', otmp2lsf(i, :)
                      CALL abort('Aborting.')
                   END IF
                END IF
             END DO
             ok = NF90_PUT_VAR(ncid, parID, REAL(otmp2lsf, 4),                    &
                  start = (/1, 1/), count = (/mland, 4/)) ! write data to file
          ELSE
             CALL abort('Parameter '//pname//                                     &
                  ' defined with unknown dimension switch - '//dimswitch//  &
                  ' - in SUBROUTINE write_output_parameter_r2')
          END IF
       END IF
    ELSE
       CALL abort('Unknown grid specification '//                               &
            '(SUBROUTINE write_output_parameter_r2)')
    END IF
    ! Check writing was successful:
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//pname//           &
         ' variable to output file (SUBROUTINE write_output_parameter_r2)')

  END SUBROUTINE write_output_parameter_r2
  !==============================================================================
  SUBROUTINE write_output_parameter_r2d(ncid, parID, pname, par_r2d, prange,   &
       writepatch, dimswitch, restart)
    ! Subroutine for writing a double precision 2D parameter (time invariant)
    ! ONLY USED FOR RESTART FILE.
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: parID ! variable's netcdf ID
    REAL(r_2), DIMENSION(:, :), INTENT(IN) :: par_r2d ! variable values
    REAL, DIMENSION(2), INTENT(IN) :: prange ! max and min for variable
    ! error checking
    LOGICAL, INTENT(IN) :: writepatch ! write patch-specific info for this var?
    LOGICAL,INTENT(IN),OPTIONAL :: restart ! are we writing to a restart file?
    CHARACTER(LEN=*), INTENT(IN) :: pname ! name of variable
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    INTEGER :: i,j ! do loop counter
    REAL(r_2),POINTER,DIMENSION(:,:,:) :: tmpout

    ! Check the nature of the parameter's second dimension:
    IF(dimswitch == 'soil') THEN ! i.e. spatial and soil
       IF(PRESENT(restart)) THEN
          ! Write data to restart file
          ok = NF90_PUT_VAR(ncid, parID, par_r2d,                              &
               start = (/1, 1/), count = (/mp, ms/))
       ELSE
          ALLOCATE(tmpout(mland, max_vegpatches, ms))
          DO i = 1, mland ! over all land grid points
             ! First write data for active patches:
             tmpout(i, 1:landpt(i)%nap, :) =                                   &
                  par_r2d(landpt(i)%cstart:landpt(i)%cend, :)
             ! Then write data for inactive patches as dummy value:
             IF(landpt(i)%nap < max_vegpatches)                                &
                  tmpout(i, (landpt(i)%nap + 1):max_vegpatches, :) =           &
                  REAL(ncmissingr, r_2)
             IF(check%ranges) THEN  ! Check ranges over active patches:
                DO j = 1, landpt(i)%nap
                   IF(ANY(tmpout(i, j, :) < prange(1)) .OR.                    &
                        ANY(tmpout(i, j, :) > prange(2)))  THEN
                      WRITE(*, *) 'Parameter '//pname//                        &
                           ' is set at a value out of specified ranges!'
                      WRITE(*, *) 'Land point # ',i, 'patch #', j
                      WRITE(*, *) 'Values: ', tmpout(i, j, :)
                      CALL abort('Aborting.')
                   END IF
                END DO
             END IF
          END DO
          ok = NF90_PUT_VAR(ncid, parID, REAL(tmpout(:, :, :), 8),             &
               start = (/1, 1, 1/),                               &
               count = (/mland, max_vegpatches, ms/)) ! write data to file
          DEALLOCATE(tmpout)
       END IF
    ELSE
       CALL abort('Parameter '//pname//                                        &
            ' defined with unknown dimension switch - '//dimswitch//     &
            ' - in SUBROUTINE write_output_parameter_r2d')
    END IF
    ! Check writing was successful:
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing '//pname//           &
         ' variable to output file (SUBROUTINE write_output_parameter_r2d)')

  END SUBROUTINE write_output_parameter_r2d


END MODULE cable_write_module
