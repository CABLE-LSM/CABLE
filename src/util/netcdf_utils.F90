! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation 
! (CSIRO) ABN 41 687 119 230.

MODULE netcdf_utils
  !! Module for common NetCDF utility procedures.
  USE netcdf
  USE iso_fortran_env, ONLY: int32, real32, real64

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: to_netcdf

  INTERFACE to_netcdf
    !* Overloads for the `to_netcdf` subroutine.
    ! 
    ! Usage:
    ! ```fortran
    ! CALL to_netcdf("data.nc", values)
    ! ```
    ! where `values` is an array and its type is supported by the following overloads.
    ! 
    ! **WARNING:** this subroutine should **not** be used to write standard
    ! CABLE outputs as the resulting NetCDF files do not include the required
    ! metadata. This subroutine is to be used for debugging purposes only. New
    ! outputs for the CABLE model should be implemented in cable_output.F90.
    MODULE PROCEDURE to_netcdf_int32_1d, to_netcdf_int32_2d, to_netcdf_int32_3d
    MODULE PROCEDURE to_netcdf_real32_1d, to_netcdf_real32_2d, to_netcdf_real32_3d
    MODULE PROCEDURE to_netcdf_real64_1d, to_netcdf_real64_2d, to_netcdf_real64_3d
  END INTERFACE to_netcdf

CONTAINS

  SUBROUTINE check(status)
    !! Check NetCDF error status.
    INTEGER, INTENT(IN) :: status !! Error status
    IF (status /= NF90_NOERR) THEN 
      PRINT *, trim(nf90_strerror(status))
      STOP
    END IF
  END SUBROUTINE check

  SUBROUTINE to_netcdf_int32_1d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    INTEGER (kind=int32), DIMENSION(:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(1), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_INT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_int32_1d

  SUBROUTINE to_netcdf_int32_2d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    INTEGER (kind=int32), DIMENSION(:,:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(2), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_INT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_int32_2d

  SUBROUTINE to_netcdf_int32_3d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    INTEGER (kind=int32), DIMENSION(:,:,:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(3), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_INT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_int32_3d

  SUBROUTINE to_netcdf_real32_1d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    REAL (kind=real32), DIMENSION(:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(1), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_FLOAT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_real32_1d

  SUBROUTINE to_netcdf_real32_2d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    REAL (kind=real32), DIMENSION(:,:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(2), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_FLOAT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_real32_2d

  SUBROUTINE to_netcdf_real32_3d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    REAL (kind=real32), DIMENSION(:,:,:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(3), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_FLOAT, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_real32_3d

  SUBROUTINE to_netcdf_real64_1d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    REAL (kind=real64), DIMENSION(:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(1), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_DOUBLE, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_real64_1d

  SUBROUTINE to_netcdf_real64_2d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    REAL (kind=real64), DIMENSION(:,:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(2), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_DOUBLE, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_real64_2d

  SUBROUTINE to_netcdf_real64_3d(filename, values)
    !* Dump values to NetCDF file.
    ! Any existing dataset with the same filename will be overwritten.
    CHARACTER (len=*), INTENT(IN) :: filename
    REAL (kind=real64), DIMENSION(:,:,:), INTENT(IN) :: values
    INTEGER :: ncid, varid, dimids(3), i
    CHARACTER (len=2) :: dimc
    CALL check( nf90_create(filename, NF90_CLOBBER, ncid) )
    DO i = 1, size(shape(values))
      WRITE (dimc, "(I2.2)") i
      CALL check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    END DO
    CALL check( nf90_def_var(ncid, "values", NF90_DOUBLE, dimids, varid) )
    CALL check( nf90_enddef(ncid) )
    CALL check( nf90_put_var(ncid, varid, values) )
    CALL check( nf90_close(ncid) )
  END SUBROUTINE to_netcdf_real64_3d

END MODULE netcdf_utils
