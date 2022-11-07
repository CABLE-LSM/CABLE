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
! Purpose: Read routines for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: New routines (redistr*) to do land use change
!
!
! ==============================================================================
! CALLed from:    cable_driver.f90
! MODULEs used:   cable_abort_module
!                 cable_IO_vars_module
!                 netcdf
!
! CALLs:          readpar_i
!                 readpar_r
!                 readpar_rd
!                 readpar_r2
!                 readpar_r2d
!                 nc_abort
!                 redistr*
!
MODULE cable_read_module


  USE cable_abort_module
  USE cable_def_types_mod, ONLY : ms, ncp, r_2, mland, mp, ncs, nrb, msn
  USE cable_IO_vars_module, ONLY: landpt, exists, land_x, land_y, metGrid
  USE netcdf

  IMPLICIT NONE
  PRIVATE
  PUBLIC readpar, redistr_i, redistr_r, redistr_rd, redistr_r2, redistr_r2d

  INTEGER :: ok ! netcdf error status
  INTERFACE readpar
     ! Loads a parameter from the met file - chooses subroutine
     ! below depending on number/type/dimension of arguments
     MODULE PROCEDURE readpar_i   ! for integer parameter read
     MODULE PROCEDURE readpar_r   ! for real parameter read
     MODULE PROCEDURE readpar_rd  ! for double precision real parameter read
     MODULE PROCEDURE readpar_r2  ! for 2d real parameter read
     MODULE PROCEDURE readpar_r2d ! for double precision 2d real parameter read
  END INTERFACE
  ! INTERFACE redistr
  !   MODULE PROCEDURE redistr_i
  !   MODULE PROCEDURE redistr_r
  !   MODULE PROCEDURE redistr_rd
  !   MODULE PROCEDURE redistr_r2
  !   MODULE PROCEDURE redistr_r2d
  ! END INTERFACE

CONTAINS

  SUBROUTINE readpar_i(ncid, parname, completeSet, var_i, filename,            &
       npatch, dimswitch, from_restart, INpatch)
    ! Subroutine for loading an integer-valued parameter
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: npatch ! # of veg patches in parameter's file
    INTEGER, INTENT(IN),OPTIONAL :: INpatch
    INTEGER, DIMENSION(:), INTENT(INOUT) :: var_i ! returned parameter
    ! values
    LOGICAL, INTENT(IN),OPTIONAL :: from_restart ! reading from restart file?
    LOGICAL, INTENT(INOUT) :: completeSet ! has every parameter been loaded?
    CHARACTER(LEN=*), INTENT(IN) :: parname ! name of parameter
    CHARACTER(LEN=*), INTENT(IN) :: filename ! file containing parameter values
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimension of
    ! parameter

    INTEGER :: parID ! parameter's netcdf ID
    INTEGER :: pardims ! # dimensions of parameter
    INTEGER :: i ! do loop counter
    INTEGER, DIMENSION(1) :: data1i ! temporary for ncdf read in
    INTEGER, DIMENSION(1, 1) :: data2i ! temporary for ncdf read in
    INTEGER, DIMENSION(:, :), POINTER :: tmp2i ! temporary for ncdf read in
    INTEGER, DIMENSION(:, :, :), POINTER :: tmp3i ! temporary for ncdf read
    ! in

    ! Check if parameter exists:
    ok = NF90_INQ_VARID(ncid,parname, parID)
    IF(ok /= NF90_NOERR) THEN ! if it doesn't exist
       completeSet=.FALSE.
       ! If this routine is reading from the restart, abort
       IF(PRESENT(from_restart))  WRITE(*,*) ' Error reading '//parname//' in file ' &
            //TRIM(filename)//' (SUBROUTINE readpar_i)'
    ELSE
       exists%parameters = .TRUE. ! Note that pars were found in file
       ! Check for grid type - restart file uses land type grid
       IF(metGrid == 'land' .OR. PRESENT(from_restart)) THEN
          ! Collect data from land only grid in netcdf file.
          ! First, check whether parameter has patch dimension:
          ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
          IF(pardims == 1) THEN ! no patch dimension, just a single land
             ! dimension
             IF(PRESENT(from_restart)) THEN
                ok = NF90_GET_VAR(ncid, parID, var_i, start=(/1/),             &
                     count=(/INpatch/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok,'Error reading '//parname//' in file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_i)')
             ELSE
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, data1i, start=(/i/),         &
                        count=(/1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok, 'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_i)')
                   ! Write non-patch-specific value to all patches:
                   var_i(landpt(i)%cstart:landpt(i)%cend) = data1i(1)
                END DO
             END IF
          ELSE IF(pardims == 2) THEN ! i.e. parameter has a patch dimension
             ALLOCATE(tmp2i(1, npatch))
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp2i,                          &
                     start=(/i,1/), count=(/1,npatch/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in met data file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_i)')
                ! Set values for this par for the # patches that exist
                var_i(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1)) =      &
                     tmp2i(1, :)
             END DO
             DEALLOCATE(tmp2i)
          ELSE
             CALL abort('Dimension of '//parname//' parameter in '//           &
                  TRIM(filename)//' unknown.')
          END IF
       ELSE IF(metGrid == 'mask') THEN ! Get data from land/sea mask type grid:
          ! First, check whether parameter has patch dimension:
          ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
          IF(pardims == 2) THEN ! i.e. no patch dimension, just x-y grid
             ! dimensions
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, data2i,                         &
                     start=(/land_x(i),land_y(i)/), count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_i)')
                ! Set all patches to have the same value for this par:
                var_i(landpt(i)%cstart:landpt(i)%cend) = data2i(1, 1)
             END DO
          ELSE IF(pardims == 3) THEN ! i.e. parameter has a patch dimension
             ALLOCATE(tmp3i(1, 1, npatch))
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp3i,                          &
                     start=(/land_x(i), land_y(i), 1/),           &
                     count=(/1, 1, npatch/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok,'Error reading '//parname//' in file '                &
                     //TRIM(filename)//' (SUBROUTINE readpar_i)')
                ! Set values for this par for the # patches that exist
                var_i(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1)) =      &
                     tmp3i(1, 1, :)
             END DO
             DEALLOCATE(tmp3i)
          ELSE
             CALL abort('Dimension of '//parname//' parameter in met file '//  &
                  'unknown.')
          END IF
       END IF ! gridtype land or mask
    END IF ! parameter's existence

  END SUBROUTINE readpar_i
  !=============================================================================
  SUBROUTINE readpar_r(ncid, parname, completeSet, var_r, filename,            &
       npatch, dimswitch, from_restart, INpatch)
    ! Subroutine for loading a real-valued parameter
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: npatch ! # of veg patches in parameter's file
    INTEGER, INTENT(IN), OPTIONAL :: INpatch
    REAL(KIND=4), DIMENSION(:), INTENT(INOUT) :: var_r ! returned parameter
    ! values
    LOGICAL, INTENT(IN), OPTIONAL :: from_restart ! reading from restart file?
    LOGICAL, INTENT(INOUT) :: completeSet ! has every parameter been loaded?
    CHARACTER(LEN=*), INTENT(IN) :: parname ! name of parameter
    CHARACTER(LEN=*), INTENT(IN) :: filename ! file containing parameter values
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    INTEGER :: parID ! parameter's netcdf ID
    INTEGER :: pardims ! # dimensions of parameter
    INTEGER :: i ! do loop counter
    REAL(KIND=4), DIMENSION(1) :: data1r ! temporary for ncdf read in
    REAL(KIND=4), DIMENSION(1, 1) :: data2r ! temporary for ncdf read in
    REAL(KIND=4), DIMENSION(:, :), POINTER :: tmp2r ! temporary for ncdf read in
    REAL(KIND=4), DIMENSION(:, :, :), POINTER :: tmp3r ! temporary for ncdf read
    ! in

    ! Check if parameter exists:
    ok = NF90_INQ_VARID(ncid, parname, parID)
    IF(ok /= NF90_NOERR) THEN ! if it doesn't exist
       completeSet = .FALSE.
       ! If this routine is reading from the restart, abort
       IF(PRESENT(from_restart)) CALL nc_abort(ok,'Error reading '//parname//  &
            ' in file '//TRIM(filename)// '(SUBROUTINE readpar_r)')
    ELSE
       exists%parameters = .TRUE. ! Note that pars were found in file
       ! If block to distinguish params with non-spatial dimensions:
       IF(dimswitch == 'def') THEN ! i.e. parameters with one spatial dim
          ! of length mland*maxpatches
          ! Check for grid type - restart file uses land type grid
          IF(metGrid == 'land' .OR. PRESENT(from_restart)) THEN
             ! Collect data from land only grid in netcdf file.
             ! First, check whether parameter has patch dimension:
             ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
             IF(pardims == 1) THEN ! no patch dimension, just a single land
                ! dimension
                IF(PRESENT(from_restart)) THEN
                   ok = NF90_GET_VAR(ncid, parID, var_r, start=(/1/),          &
                        count=(/INpatch/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok,'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_r)')
                ELSE
                   DO i = 1, mland ! over all land points/grid cells
                      ok = NF90_GET_VAR(ncid, parID, data1r, start=(/i/),      &
                           count=(/1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort                       &
                           (ok,'Error reading '//parname//' in file ' &
                           //TRIM(filename)//' (SUBROUTINE readpar_r)')
                      ! All patches set to the same value if no patch info:
                      var_r(landpt(i)%cstart:landpt(i)%cend) = data1r(1)
                   END DO
                END IF
             ELSE IF(pardims == 2) THEN ! i.e. parameter has a patch dimension
                ALLOCATE(tmp2r(1, npatch))
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, tmp2r,                       &
                        start=(/i, 1/), count=(/1, npatch/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok,'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_r)')
                   ! Set values for this par for the # patches that exist
                   var_r(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1)) =   &
                        tmp2r(1, :)
                END DO
                DEALLOCATE(tmp2r)
             ELSE
                CALL abort('Dimension of '//parname//                          &
                     ' parameter in met file unknown.')
             END IF
          ELSE IF(metGrid == 'mask') THEN ! Get data from land/sea mask type
             ! grid:
             ! First, check whether parameter has patch dimension:
             ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
             IF(pardims == 2) THEN ! no patch dimension, just x-y grid
                ! dimensions
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, data2r,                      &
                        start=(/land_x(i),land_y(i)/), count=(/1, 1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok,'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_r)')
                   ! Set all patches to have the same value for this par:
                   var_r(landpt(i)%cstart:landpt(i)%cend) = data2r(1, 1)
                END DO
             ELSE IF(pardims == 3) THEN ! i.e. parameter has a patch dimension
                ALLOCATE(tmp3r(1, 1, npatch))
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, tmp3r,                       &
                        start=(/land_x(i), land_y(i), 1/),        &
                        count=(/1, 1, npatch/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok,'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_r)')
                   ! Set values for this par for the # patches that exist
                   var_r(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1)) =   &
                        tmp3r(1, 1, :)
                END DO
                DEALLOCATE(tmp3r)
             ELSE
                CALL abort('Dimension of '//parname//                          &
                     ' parameter in met file unknown.')
             END IF
          END IF ! gridtype land or mask
       ELSE IF(dimswitch == 'ms') THEN ! ie par has only soil dimension, no
          ! spatial
          ! Load parameter values (e.g. zse):
          ok = NF90_GET_VAR(ncid, parID, var_r, start=(/1/), count=(/ms/))
          IF(ok /= NF90_NOERR) CALL nc_abort                                   &
               (ok,'Error reading '//parname//' in file ' &
               //TRIM(filename)//' (SUBROUTINE readpar_r)')
       ELSE IF(dimswitch == 'ncp') THEN ! ie par has only ncp dimension e.g.
          ! ratecp
          ! Load ratecp parameter values:
          ok = NF90_GET_VAR(ncid, parID, var_r, start=(/1/), count=(/ncp/))
          IF(ok /= NF90_NOERR) CALL nc_abort                                   &
               (ok, 'Error reading '//parname//' in file ' &
               //TRIM(filename)//' (SUBROUTINE readpar_r)')
       ELSE IF(dimswitch == 'ncs') THEN ! ie par has only ncs dimension e.g.
          ! ratecs
          ! Load parameter values:
          ok = NF90_GET_VAR(ncid, parID, var_r, start=(/1/), count=(/ncs/))
          IF(ok /= NF90_NOERR) CALL nc_abort                                   &
               (ok,'Error reading '//parname//' in file ' &
               //TRIM(filename)//' (SUBROUTINE readpar_r)')
       ELSE
          CALL abort('Parameter or initial state '//parname//                  &
               ' called with unknown dimension switch - '//dimswitch//   &
               ' - in INTERFACE readpar')
       END IF ! dimension of parameter i.e. is this zse or ratecp or ratecs
    END IF ! parameter's existence

  END SUBROUTINE readpar_r
  !=============================================================================
  SUBROUTINE readpar_rd(ncid, parname, completeSet, var_rd, filename,          &
       npatch, dimswitch, from_restart, INpatch)
    ! Subroutine for loading a double precision real-valued parameter
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: npatch ! # of veg patches in parameter's file
    REAL(r_2), DIMENSION(:), INTENT(INOUT) :: var_rd ! returned parameter
    ! values
    LOGICAL, INTENT(IN), OPTIONAL :: from_restart ! reading from restart file?
    LOGICAL, INTENT(INOUT) :: completeSet ! has every parameter been loaded?
    CHARACTER(LEN=*), INTENT(IN) :: parname ! name of parameter
    CHARACTER(LEN=*), INTENT(IN) :: filename ! file containing parameter values
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of
    ! parameter
    INTEGER, INTENT(IN),OPTIONAL :: INpatch
    INTEGER :: parID ! parameter's netcdf ID
    INTEGER :: pardims ! # dimensions of parameter
    INTEGER :: i ! do loop counter
    REAL(4), DIMENSION(1) :: data1r ! temporary for ncdf read in
    REAL(4), DIMENSION(1, 1) :: data2r ! temporary for ncdf read in
    REAL(4), DIMENSION(:), POINTER :: tmp1r ! temporary for ncdf read in
    REAL(4), DIMENSION(:, :), POINTER :: tmp2r ! temporary for ncdf read in
    REAL(4), DIMENSION(:, :, :), POINTER :: tmp3r ! temporary for ncdf read in

    ! Check if parameter exists:
    ok = NF90_INQ_VARID(ncid, parname, parID)
    IF(ok /= NF90_NOERR) THEN ! if it doesn't exist
       completeSet = .FALSE.
       ! If this routine is reading from the restart, abort
       IF(PRESENT(from_restart)) CALL nc_abort(ok,'Error reading '//parname//  &
            ' in file '//TRIM(filename)// '(SUBROUTINE readpar_rd)')
    ELSE
       exists%parameters = .TRUE. ! Note that pars were found in file
       ! If block to distinguish params with non-spatial dimensions:
       IF(dimswitch(1:3) == 'def') THEN ! i.e. parameters with one spatial dim
          ! of length mland*maxpatches
          ! Check for grid type - restart file uses land type grid
          IF(metGrid == 'land' .OR. PRESENT(from_restart)) THEN
             ! Collect data from land only grid in netcdf file.
             ! First, check whether parameter has patch dimension:
             ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
             IF(pardims == 1) THEN ! no patch dimension, just a single land
                ! dimension
                IF(PRESENT(from_restart)) THEN
                   !                   IF(dimswitch(1:4) == 'defd') THEN ! ie we're expecting to
                   ! read double prec.
                   !                      ! Read double precision data:
                   !                      ok = NF90_GET_VAR(ncid,parID,var_rd,start=(/1/),        &
                   !                                        count=(/INpatch/))
                   !                   ELSE ! ie we're reading single prec. var in netdf file
                   ! Read single precision data:
                   ALLOCATE(tmp1r(INpatch))
                   ok = NF90_GET_VAR(ncid, parID, tmp1r, start=(/1/),       &
                        count=(/INpatch/))
                   var_rd = REAL(tmp1r, r_2)
                   DEALLOCATE(tmp1r)
                   !                   END IF
                   IF(ok /= NF90_NOERR) CALL nc_abort                       &
                        (ok,'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_rd)')
                ELSE ! reading from met file
                   DO i = 1, mland ! over all land points/grid cells
                      ok = NF90_GET_VAR(ncid, parID, data1r, start=(/i/),      &
                           count=(/1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort                       &
                           (ok,'Error reading '//parname//' in file ' &
                           //TRIM(filename)//' (SUBROUTINE readpar_rd)')
                      ! Give single value to all patches if no patch specific
                      ! info:
                      var_rd(landpt(i)%cstart:landpt(i)%cend) =                &
                           REAL(data1r(1))
                   END DO
                END IF
             ELSE IF(pardims == 2) THEN ! i.e. parameter has a patch dimension
                ! NB restart file will not have a patch dimension, therefore
                ! all reads here are of single precision variables.
                ALLOCATE(tmp2r(1, npatch))
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, tmp2r,                       &
                        start=(/i, 1/), count=(/1, npatch/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok, 'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_rd)')
                   ! Set values for this par for the # patches that exist
                   var_rd(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1)) =  &
                        REAL(tmp2r(1, :), r_2)
                END DO
                DEALLOCATE(tmp2r)
             ELSE
                CALL abort('Dimension of '//parname//                          &
                     ' parameter in met file unknown.')
             END IF
          ELSE IF(metGrid == 'mask') THEN ! Get data from land/sea mask type
             ! grid:
             ! NB restart file will not have mask grid, therefore all reads
             ! here are of single precision variables in the netcdf met file
             ! First, check whether parameter has patch dimension:
             ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
             IF(pardims == 2) THEN ! no patch dimension, just x-y grid
                ! dimensions
                DO i=1, mland ! over all land points/grid cells
                   ok= NF90_GET_VAR(ncid, parID, data2r,                       &
                        start=(/land_x(i), land_y(i)/),            &
                        count=(/1, 1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok, 'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_rd)')
                   ! Set all patches to have the same value for this par:
                   var_rd(landpt(i)%cstart:landpt(i)%cend) =                   &
                        REAL(data2r(1, 1), r_2)
                END DO
             ELSE IF(pardims == 3) THEN ! i.e. parameter has a patch dimension
                ALLOCATE(tmp3r(1, 1, npatch))
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, tmp3r,                       &
                        start=(/land_x(i), land_y(i), 1/),        &
                        count=(/1, 1, npatch/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok, 'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_rd)')
                   ! Set values for this par for the # patches that exist
                   var_rd(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1))    &
                        = REAL(tmp3r(1, 1, :), r_2)
                END DO
                DEALLOCATE(tmp3r)
             ELSE
                CALL abort('Dimension of '//parname//                          &
                     ' parameter in met file unknown.')
             END IF
          ELSE
             CALL abort('Prescribed input grid '//metGrid//' unknown.')
          END IF ! gridtype land or mask

       ELSE IF(dimswitch(1:2) == 'ms') THEN ! ie par has only soil dimension,
          ! no spatial
          ! Load parameter values (e.g. zse):
          DO i = 1, ms
             ok = NF90_GET_VAR(ncid, parID, data1r, start=(/i/), count=(/1/))
             IF(ok /= NF90_NOERR) CALL nc_abort                                &
                  (ok, 'Error reading '//parname//' in file ' &
                  //TRIM(filename)//' (SUBROUTINE readpar_rd)')
             var_rd(i) = REAL(data1r(1), r_2)
          END DO
       ELSE IF(dimswitch(1:3) == 'ncp') THEN ! ie par has only ncp dimension
          ! e.g. ratecp
          ! Load ratecp parameter values:
          DO i = 1, ncp
             ok = NF90_GET_VAR(ncid, parID, data1r, start=(/i/), count=(/1/))
             IF(ok /= NF90_NOERR) CALL nc_abort                                &
                  (ok, 'Error reading '//parname//' in file ' &
                  //TRIM(filename)//' (SUBROUTINE readpar_rd)')
             var_rd(i) = REAL(data1r(1), r_2)
          END DO
       ELSE IF(dimswitch(1:3) == 'ncs') THEN ! ie par has only ncs dimension
          ! e.g. ratecs
          ! Load parameter values:
          DO i = 1, ncs
             ok = NF90_GET_VAR(ncid, parID, data1r, start=(/i/), count=(/1/))
             IF(ok /= NF90_NOERR) CALL nc_abort                                &
                  (ok,'Error reading '//parname//' in file ' &
                  //TRIM(filename)//' (SUBROUTINE readpar_rd)')
             var_rd(i) = REAL(data1r(1), r_2)
          END DO
       ELSE
          CALL abort('Parameter or initial state '//parname//                  &
               ' called with unknown dimension switch - '//dimswitch//   &
               ' - in INTERFACE readpar')
       END IF ! dimension of parameter i.e. is this zse or ratecp or ratecs
    END IF ! parameter's existence

  END SUBROUTINE readpar_rd
  !=============================================================================
  SUBROUTINE readpar_r2(ncid, parname, completeSet, var_r2, filename,          &
       npatch, dimswitch, from_restart, INpatch)
    ! Subroutine for loading a two dimensional real-valued parameter
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: npatch ! number of veg patches in file
    INTEGER, INTENT(IN), OPTIONAL :: INpatch
    REAL(KIND=4), DIMENSION(:,:), INTENT(INOUT) :: var_r2 ! returned parameter
    ! values
    LOGICAL, INTENT(IN), OPTIONAL :: from_restart ! reading from restart file?
    LOGICAL, INTENT(INOUT) :: completeSet ! has every parameter been loaded?
    CHARACTER(LEN=*), INTENT(IN) :: filename ! file containing parameter values
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of
    ! parameter
    CHARACTER(LEN=*), INTENT(IN) :: parname ! name of parameter

    INTEGER :: parID ! parameter's netcdf ID
    INTEGER :: pardims ! # dimensions of parameter
    INTEGER :: dimctr ! size of non-spatial (2nd) dimension of parameter
    INTEGER :: i, j ! do loop counter
    REAL(KIND=4), DIMENSION(:, :), POINTER       :: tmp2r ! temporary for ncdf
    ! read in
    REAL(KIND=4), DIMENSION(:, :, :), POINTER    :: tmp3r ! temporary for ncdf
    ! read in
    REAL(KIND=4), DIMENSION(:, :, :, :), POINTER :: tmp4r ! temporary for ncdf
    ! read in
    REAL :: tmpjh

    ! Check if parameter exists:
    ok = NF90_INQ_VARID(ncid, parname, parID)
    IF(ok /= NF90_NOERR) THEN ! if it doesn't exist
       completeSet = .FALSE.
       ! If this routine is reading from the restart, abort
       IF(PRESENT(from_restart)) CALL nc_abort(ok, 'Error reading '//parname   &
            //' in file '//TRIM(filename)// '(SUBROUTINE readpar_r2)')
    ELSE
       exists%parameters = .TRUE. ! Note that pars were found in file
       ! Decide which 2nd dimension of parameter/init state we're loading:
       IF(dimswitch == 'ms') THEN
          dimctr = ms ! i.e. horizontal spatial and soil
       ELSE IF(dimswitch == 'snow') THEN
          dimctr = msn ! i.e. horizontal spatial and snow
       ELSE IF(dimswitch == 'nrb') THEN
          dimctr = nrb ! i.e. horizontal spatial and radiation bands
       ELSE IF(dimswitch == 'ncp') THEN
          dimctr = ncp ! i.e. horizontal spatial and plant carbon pools
       ELSE IF(dimswitch == 'ncs') THEN
          dimctr = ncs ! i.e. horizontal spatial and soil carbon pools
       ELSE
          CALL abort('Parameter or initial state '//parname//                  &
               ' called with unknown dimension switch - '//dimswitch//   &
               ' - in INTERFACE readpar SUBROUTINE readpar_r2')
       END IF
       ! Check for grid type - restart file uses land type grid
       IF(metGrid == 'land' .OR. PRESENT(from_restart)) THEN
          ! Collect data from land only grid in netcdf file.
          ! First, check whether parameter has patch dimension:
          ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
          IF(pardims == 2) THEN ! no patch dimension, just a land + other
             ! dimension
             IF(PRESENT(from_restart)) THEN
                ok = NF90_GET_VAR(ncid, parID, var_r2, start=(/1, 1/),         &
                     count=(/INpatch, dimctr/))
                !    WRITE(45,*) 'Is tgg read here?'
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in file' &
                     //TRIM(filename)//' (SUBROUTINE readpar_rd)')
             ELSE
                ALLOCATE(tmp2r(1, dimctr))
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, tmp2r,                       &
                        start=(/i, 1/), count=(/1, dimctr/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok, 'Error reading '//parname//' in file ' &
                        //TRIM(filename)//' (SUBROUTINE readpar_r2)')
                   DO j = 1, dimctr
                      ! Set all patches to have the same value:
                      var_r2(landpt(i)%cstart:landpt(i)%cend, j) = tmp2r(1, j)
                   END DO
                END DO
                DEALLOCATE(tmp2r)
             END IF
          ELSE IF(pardims == 3) THEN ! i.e. parameter has a patch dimension
             ALLOCATE(tmp3r(1, npatch, dimctr))
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp3r,                          &
                     start=(/i, 1, 1/), count=(/1, npatch, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2)')
                DO j = 1, dimctr
                   ! Set values for this par for the # patches that exist
                   var_r2(landpt(i)%cstart:(landpt(i)%cstart + npatch          &
                        - 1), j) = tmp3r(1,:,j)
                END DO
             END DO
             DEALLOCATE(tmp3r)
          ELSE
             CALL abort('Dimension of '//parname//' parameter in met file '//  &
                  'unknown.')
          END IF
       ELSEIF(metGrid == 'mask') THEN ! Get data from land/sea mask type grid:
          ! First, check whether parameter has patch dimension:
          ok=NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)

          IF(pardims == 3) THEN ! no patch dimension, just x-y + soil grid
             ! dimension
             ALLOCATE(tmp3r(1, 1, ms))
             DO i = 1, mland ! over all land points/grid cells
                ok =  NF90_GET_VAR(ncid, parID, tmp3r,                         &
                     start=(/land_x(i), land_y(i), 1/),          &
                     count=(/1, 1, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2)')
                ! Set all patches to have the same value for this par:
                DO j = 1, dimctr
                   var_r2(landpt(i)%cstart:landpt(i)%cend,j) = tmp3r(1,1,j)
                END DO
             END DO
             DEALLOCATE(tmp3r)
          ELSE IF(pardims == 4) THEN ! i.e. soil parameter has a patch dimension
             ALLOCATE(tmp4r(1, 1, npatch, dimctr))
             DO  i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp4r,                          &
                     start=(/land_x(i),land_y(i), 1, 1/),         &
                     count=(/1, 1, npatch, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2)')
                DO j = 1, dimctr
                   ! Set values for this par for the # patches that exist
                   var_r2(landpt(i)%cstart:(landpt(i)%cstart + npatch -1 ), j) &
                        = tmp4r(1, 1, :, j)
                END DO
             END DO
             DEALLOCATE(tmp4r)
          ELSE IF(dimswitch == 'nrb') THEN
             PRINT *, 'pardims', pardims
             IF(pardims == 2) THEN ! no patch dimension, just a land + other
                ! dimension
                IF(PRESENT(from_restart)) THEN
                   ok = NF90_GET_VAR(ncid, parID, var_r2, start=(/1, 1/),      &
                        count=(/INpatch, dimctr/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok,'Error reading '//parname//' in file' &
                        //TRIM(filename)//' (SUBROUTINE readpar_rd)')
                ELSE
                   DO i = 1, mland ! over all land points/grid cells
                      ok = NF90_GET_VAR(ncid, parID, tmpjh, start=(/i/))
                      IF(ok /= NF90_NOERR) CALL nc_abort                          &
                           (ok, 'Error reading '//parname//' in file ' &
                           //TRIM(filename)//' (SUBROUTINE readpar_r2)')
                      ! Set all patches to have the same value:
                      PRINT *, 'read albsoil: ', tmpjh
                      var_r2(landpt(i)%cstart:landpt(i)%cend, 1) = tmpjh
                      var_r2(landpt(i)%cstart:landpt(i)%cend, 2) = tmpjh
                   END DO
                END IF
             END IF
          ELSE
             CALL abort('Dimension of '//parname//' parameter in met file '//  &
                  'unknown.')
          END IF
       END IF ! gridtype land or mask
    END IF ! parameter's existence

  END SUBROUTINE readpar_r2
  !=============================================================================
  SUBROUTINE readpar_r2d(ncid, parname, completeSet, var_r2d, filename,        &
       npatch, dimswitch, from_restart, INpatch)
    ! Subroutine for loading a double precision two dimensional real-valued
    ! soil dimensioned parameter
    INTEGER, INTENT(IN) :: ncid ! netcdf file ID
    INTEGER, INTENT(IN) :: npatch ! # of veg patches in parameter's file
    INTEGER, INTENT(IN), OPTIONAL :: INpatch
    REAL(r_2), DIMENSION(:, :), INTENT(INOUT) :: var_r2d ! returned parameter
    ! value
    LOGICAL, INTENT(IN), OPTIONAL :: from_restart ! reading from restart file?
    LOGICAL, INTENT(INOUT) :: completeSet ! has every parameter been loaded?
    CHARACTER(LEN=*), INTENT(IN) :: parname ! name of parameter
    CHARACTER(LEN=*), INTENT(IN) :: filename ! file containing parameter values
    CHARACTER(LEN=*), INTENT(IN) :: dimswitch ! indicates dimesnion of parameter

    INTEGER :: parID ! parameter's netcdf ID
    INTEGER :: pardims ! # dimensions of parameter
    INTEGER :: dimctr ! size of non-spatial (2nd) dimension of parameter
    INTEGER :: i,j ! do loop counter
    REAL(8), DIMENSION(:, :), POINTER       :: tmp2rd ! temporary for ncdf
    ! read in
    REAL(4), DIMENSION(:, :), POINTER       :: tmp2r  ! temporary for ncdf
    ! read in
    REAL(4), DIMENSION(:, :, :), POINTER    :: tmp3r  ! temporary for ncdf
    ! read in
    REAL(4), DIMENSION(:, :, :, :), POINTER :: tmp4r  ! temporary for ncdf
    ! read in

    ! Check if parameter exists:
    ok = NF90_INQ_VARID(ncid,parname,parID)
    IF(ok /= NF90_NOERR) THEN ! if it doesn't exist
       completeSet = .FALSE.
       ! If this routine is reading from the restart, abort
       IF(PRESENT(from_restart))  WRITE(*,*) ' Error reading '//parname//' in file ' &
            //TRIM(filename)//' (SUBROUTINE readpar_r2d)'
    ELSE
       exists%parameters = .TRUE. ! Note that pars were found in file
       ! Decide which 2nd dimension of parameter /init state we're loading:
       IF(dimswitch(1:2) == 'ms') THEN
          dimctr = ms ! i.e. horizontal spatial and soil
       ELSE IF(dimswitch(1:4) == 'snow') THEN
          dimctr = msn ! i.e. horizontal spatial and snow
       ELSE IF(dimswitch(1:3) == 'nrb') THEN
          dimctr = nrb ! i.e. horizontal spatial and radiation bands
       ELSE IF(dimswitch(1:3) == 'ncp') THEN
          dimctr = ncp ! i.e. horizontal spatial and plant carbon pools
       ELSE IF(dimswitch(1:3) == 'ncs') THEN
          dimctr = ncs ! i.e. horizontal spatial and soil carbon pools
       ELSE
          CALL abort('Parameter or initial state '//parname//                  &
               ' called with unknown dimension switch - '//dimswitch//   &
               ' - in INTERFACE readpar')
       END IF
       ! Check for grid type - restart file uses land type grid
       IF(metGrid == 'land' .OR. PRESENT(from_restart)) THEN
          ! Collect data from land only grid in netcdf file.
          ! First, check whether parameter has patch dimension:
          ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
          IF(pardims == 2) THEN ! no patch dimension, just a land+soil
             ! dimensions
             ! If we really are reading a double precision variable
             ! from the netcdf restart file, dimswitch will show this:
             ! equivalent to using "IF(PRESENT(from_restart)) THEN"
             IF(dimswitch == 'msd' .OR. dimswitch == 'snowd' .OR.              &
                  dimswitch == 'nrbd' .OR. dimswitch == 'ncpd'                   &
                  .OR. dimswitch == 'ncsd') THEN
                ALLOCATE(tmp2rd(INpatch, dimctr))
                ok = NF90_GET_VAR(ncid, parID, tmp2rd,                          &
                     start=(/1, 1/), count=(/INpatch, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                              &
                     (ok, 'Error reading '//parname//' in met data file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2d)')
                var_r2d(:, :) = REAL(tmp2rd(:, :), r_2)
                DEALLOCATE(tmp2rd)
                !              ALLOCATE(tmp2rd(1,dimctr))
                !              DO i=1, mland ! over all land points/grid cells
                !                 ok= NF90_GET_VAR(ncid, parID, tmp2rd,                        &
                !                                  start=(/i,1/), count=(/1,dimctr/))
                !                 IF(ok /= NF90_NOERR) CALL nc_abort                           &
                !                    (ok,'Error reading '//parname//' in met data file '       &
                !                      //TRIM(filename)//' (SUBROUTINE readpar_r2d)')
                !                 DO j=1, dimctr
                !                    var_r2d(landpt(i)%cstart:landpt(i)%cend,j) =              &
                !                      REAL(tmp2rd(1,j),r_2)
                !                 END DO
                !              END DO
                !              DEALLOCATE(tmp2rd)
                ! WRITE(45,*) 'After read-in restart values'
                ! WRITE(45,*) '1039_var_r2d = ', var_r2d(1039,:)
                ! WRITE(45,*) '1672_var_r2d = ', var_r2d(1672,:)
             ELSE
                ALLOCATE(tmp2r(1, dimctr))
                DO i = 1, mland ! over all land points/grid cells
                   ok = NF90_GET_VAR(ncid, parID, tmp2r,                       &
                        start=(/i, 1/), count=(/1, dimctr/))
                   IF(ok /= NF90_NOERR) CALL nc_abort                          &
                        (ok, 'Error reading '//parname//' in met data file '    &
                        //TRIM(filename)//' (SUBROUTINE readpar_r2d)')
                   DO j = 1, dimctr
                      var_r2d(landpt(i)%cstart:landpt(i)%cend, j) =            &
                           REAL(tmp2r(1,j))
                   END DO
                END DO
                DEALLOCATE(tmp2r)
             END IF ! reading a d.p. var from netcdf
          ELSE IF(pardims == 3) THEN ! i.e. parameter has a patch dimension
             ! Note that restart file doesn't have a patch dimension,
             ! so that reads below are of single precision vares from met file
             ALLOCATE(tmp3r(1,npatch,dimctr))
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp3r,                          &
                     start=(/i, 1, 1/), count=(/1, npatch, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in met data file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2d)')
                DO j = 1, dimctr
                   ! Set values for this par for the # patches that exist
                   var_r2d(landpt(i)%cstart:(landpt(i)%cstart + npatch - 1),j) &
                        = REAL(tmp3r(1,:,j),r_2)
                END DO
             END DO
             DEALLOCATE(tmp3r)
          ELSE
             CALL abort('Dimension of '//parname//' parameter in met file'//   &
                  'unknown.')
          END IF
       ELSEIF(metGrid == 'mask') THEN ! Get data from land/sea mask type grid:
          ! NB restart file won't have mask grid, therefore below we are
          ! reading single precision variables from the met file
          ! First, check whether parameter has patch dimension:
          ok = NF90_INQUIRE_VARIABLE(ncid, parID, ndims=pardims)
          IF(pardims == 3) THEN ! no patch dimension, just x-y + soil grid
             ! dimension
             ALLOCATE(tmp3r(1, 1, dimctr))
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp3r,                          &
                     start=(/land_x(i), land_y(i), 1/),           &
                     count=(/1, 1, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok, 'Error reading '//parname//' in met data file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2d)')
                ! Set all patches to have the same value for this par:
                DO j = 1, dimctr
                   var_r2d(landpt(i)%cstart:landpt(i)%cend, j) =               &
                        REAL(tmp3r(1, 1, j), r_2)
                END DO
             END DO
             DEALLOCATE(tmp3r)
          ELSE IF(pardims == 4) THEN ! i.e. soil parameter has a patch dimension
             ALLOCATE(tmp4r(1,1,npatch,dimctr))
             DO i = 1, mland ! over all land points/grid cells
                ok = NF90_GET_VAR(ncid, parID, tmp4r,                          &
                     start=(/land_x(i), land_y(i), 1, 1/),        &
                     count=(/1, 1, npatch, dimctr/))
                IF(ok /= NF90_NOERR) CALL nc_abort                             &
                     (ok,'Error reading '//parname//' in met data file ' &
                     //TRIM(filename)//' (SUBROUTINE readpar_r2d)')
                DO j = 1, dimctr
                   ! Set values for this par for the # patches that exist
                   var_r2d(landpt(i)%cstart:(landpt(i)%cstart + npatch         &
                        - 1), j) = REAL(tmp4r(1, 1, :, j), r_2)
                END DO
             END DO
             DEALLOCATE(tmp4r)
          ELSE
             CALL abort('Dimension of '//parname//' parameter in met file'//   &
                  'unknown.')
          END IF
       END IF ! gridtype land or mask
    END IF ! parameter's existence

  END SUBROUTINE readpar_r2d
  !=============================================================================
  SUBROUTINE redistr_i(INpatch, nap, in_i, out_i, parname)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)  :: INpatch
    INTEGER,     INTENT(IN)  :: nap(INpatch)
    INTEGER,     INTENT(IN)  :: in_i(INpatch)
    INTEGER,     INTENT(OUT) :: out_i(mp)
    CHARACTER(LEN=*), INTENT(IN)  :: parname ! name of parameter

    ! local variables
    !    REAL    :: ave_r
    INTEGER :: ii, jj, npt

    npt = 0
    DO ii = 1, mland
       !     ave_r = 0.0
       DO jj = 1, nap(ii)
          npt = npt + 1
          !       ave_r = ave_r + in_i(npt)
       END DO
       !     ave_r = ave_r / FLOAT(nap(ii))
       !     out_i(landpt(ii)%cstart:landpt(ii)%cend) = INT(ave_r)
       ! just take the dominant one for isflag
       out_i(landpt(ii)%cstart:landpt(ii)%cend) = in_i(npt-nap(ii) + 1)
    END DO
    IF (npt /= INpatch) THEN
       PRINT *, parname,' Error: npt /= INpatch, ',npt, INpatch
       STOP
    END IF

  END SUBROUTINE redistr_i

  SUBROUTINE redistr_r(INpatch, nap, in_r, out_r, parname)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)  :: INpatch
    INTEGER,     INTENT(IN)  :: nap(INpatch)
    REAL,        INTENT(IN)  :: in_r(INpatch)
    REAL,        INTENT(OUT) :: out_r(mp)
    CHARACTER(LEN=*), INTENT(IN)  :: parname ! name of parameter

    ! local variables
    REAL    :: ave_r
    INTEGER :: ii, jj, npt

    npt = 0
    DO ii = 1, mland
       ave_r = 0.0
       DO jj = 1, nap(ii)
          npt = npt + 1
          ave_r = ave_r + in_r(npt)
       END DO
       ave_r = ave_r / FLOAT(nap(ii))
       out_r(landpt(ii)%cstart:landpt(ii)%cend) = ave_r
    END DO
    IF (npt /= INpatch) THEN
       PRINT *, parname,' Error: npt /= INpatch, ',npt, INpatch
       STOP
    END IF

  END SUBROUTINE redistr_r

  SUBROUTINE redistr_rd(INpatch,nap,in_rd,out_rd,parname)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)  :: INpatch
    INTEGER,     INTENT(IN)  :: nap(INpatch)
    REAL(r_2),        INTENT(IN)  :: in_rd(INpatch)
    REAL(r_2),        INTENT(OUT) :: out_rd(mp)
    CHARACTER(LEN=*), INTENT(IN)  :: parname ! name of parameter

    ! local variables
    REAL(r_2)    :: ave_rd
    INTEGER :: ii, jj, npt

    npt = 0
    DO ii = 1, mland
       ave_rd = 0.0
       DO jj = 1, nap(ii)
          npt = npt + 1
          ave_rd = ave_rd + in_rd(npt)
       END DO
       ave_rd = ave_rd / FLOAT(nap(ii))
       out_rd(landpt(ii)%cstart:landpt(ii)%cend) = ave_rd
    END DO
    IF (npt /= INpatch) THEN
       PRINT *, parname,' Error: npt /= INpatch, ',npt, INpatch
       STOP
    END IF

  END SUBROUTINE redistr_rd

  SUBROUTINE redistr_r2(INpatch, nap,in_r2, out_r2, parname, dim2)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)  :: INpatch
    INTEGER,     INTENT(IN)  :: dim2
    INTEGER,     INTENT(IN)  :: nap(INpatch)
    REAL,        INTENT(IN)  :: in_r2 (INpatch,dim2)
    REAL,        INTENT(OUT) :: out_r2(mp,dim2)
    CHARACTER(LEN=*), INTENT(IN)  :: parname ! name of parameter

    ! local variables
    REAL    :: ave_r2(dim2)
    INTEGER :: ii, jj, npt

    npt = 0
    DO ii = 1, mland
       ave_r2(:) = 0.0
       DO jj = 1, nap(ii)
          npt = npt + 1
          ave_r2(:) = ave_r2(:) + in_r2(npt,:)
       END DO
       ave_r2(:) = ave_r2(:) / FLOAT(nap(ii))
       DO jj = 1, dim2
          out_r2(landpt(ii)%cstart:landpt(ii)%cend, jj) = ave_r2(jj)
       END DO
    END DO
    IF (npt /= INpatch) THEN
       PRINT *, parname,' Error: npt /= INpatch, ',npt, INpatch
       STOP
    END IF

  END SUBROUTINE redistr_r2

  SUBROUTINE redistr_r2d(INpatch,nap,in_r2d,out_r2d,parname,dim2)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)  :: INpatch
    INTEGER,     INTENT(IN)  :: dim2
    INTEGER,     INTENT(IN)  :: nap(INpatch)
    REAL(r_2),        INTENT(IN)  :: in_r2d (INpatch,dim2)
    REAL(r_2),        INTENT(OUT) :: out_r2d(mp,dim2)
    CHARACTER(LEN=*), INTENT(IN)  :: parname ! name of parameter

    ! local variables
    REAL(r_2)    :: ave_r2d(dim2)
    INTEGER :: ii, jj, npt

    npt = 0
    DO ii = 1, mland
       ave_r2d(:) = 0.0
       DO jj = 1, nap(ii)
          npt = npt + 1
          ave_r2d(:) = ave_r2d(:) + in_r2d(npt,:)
       END DO
       ave_r2d(:) = ave_r2d(:) / FLOAT(nap(ii))
       DO jj = 1, dim2
          out_r2d(landpt(ii)%cstart:landpt(ii)%cend, jj) = ave_r2d(jj)
       END DO
    END DO
    IF (npt /= INpatch) THEN
       PRINT *, parname,' Error: npt /= INpatch, ',npt, INpatch
       STOP
    END IF

  END SUBROUTINE redistr_r2d


END MODULE cable_read_module
