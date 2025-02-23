module netcdf_utils
  use netcdf
  use iso_fortran_env, only: int32, real32, real64

  implicit none

  private

  public :: to_netcdf

  interface to_netcdf
    module procedure to_netcdf_int32_1d
    module procedure to_netcdf_real32_1d, to_netcdf_real32_2d, to_netcdf_real32_3d
    module procedure to_netcdf_real64_1d, to_netcdf_real64_2d
  end interface to_netcdf

contains

  subroutine check(status)
    integer, intent(in) :: status
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop
    end if
  end subroutine check

  subroutine to_netcdf_int32_1d(filename, values)
    character (len=*), intent(in) :: filename
    integer (kind=int32), dimension(:), intent(in) :: values
    integer :: ncid, varid, dimids(1), i
    character (len=2) :: dimc
    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    do i = 1, size(shape(values))
      write (dimc, "(I2.2)") i
      call check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    end do
    call check( nf90_def_var(ncid, "values", NF90_INT, dimids, varid) )
    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, values) )
    call check( nf90_close(ncid) )
  end subroutine to_netcdf_int32_1d

  subroutine to_netcdf_real32_1d(filename, values)
    character (len=*), intent(in) :: filename
    real (kind=real32), dimension(:), intent(in) :: values
    integer :: ncid, varid, dimids(1), i
    character (len=2) :: dimc
    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    do i = 1, size(shape(values))
      write (dimc, "(I2.2)") i
      call check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    end do
    call check( nf90_def_var(ncid, "values", NF90_FLOAT, dimids, varid) )
    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, values) )
    call check( nf90_close(ncid) )
  end subroutine to_netcdf_real32_1d

  subroutine to_netcdf_real32_2d(filename, values)
    character (len=*), intent(in) :: filename
    real (kind=real32), dimension(:,:), intent(in) :: values
    integer :: ncid, varid, dimids(2), i
    character (len=2) :: dimc
    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    do i = 1, size(shape(values))
      write (dimc, "(I2.2)") i
      call check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    end do
    call check( nf90_def_var(ncid, "values", NF90_FLOAT, dimids, varid) )
    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, values) )
    call check( nf90_close(ncid) )
  end subroutine to_netcdf_real32_2d

  subroutine to_netcdf_real32_3d(filename, values)
    character (len=*), intent(in) :: filename
    real (kind=real32), dimension(:,:,:), intent(in) :: values
    integer :: ncid, varid, dimids(3), i
    character (len=2) :: dimc
    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    do i = 1, size(shape(values))
      write (dimc, "(I2.2)") i
      call check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    end do
    call check( nf90_def_var(ncid, "values", NF90_FLOAT, dimids, varid) )
    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, values) )
    call check( nf90_close(ncid) )
  end subroutine to_netcdf_real32_3d

  subroutine to_netcdf_real64_1d(filename, values)
    character (len=*), intent(in) :: filename
    real (kind=real64), dimension(:), intent(in) :: values
    integer :: ncid, varid, dimids(1), i
    character (len=2) :: dimc
    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    do i = 1, size(shape(values))
      write (dimc, "(I2.2)") i
      call check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    end do
    call check( nf90_def_var(ncid, "values", NF90_DOUBLE, dimids, varid) )
    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, values) )
    call check( nf90_close(ncid) )
  end subroutine to_netcdf_real64_1d

  subroutine to_netcdf_real64_2d(filename, values)
    character (len=*), intent(in) :: filename
    real (kind=real64), dimension(:,:), intent(in) :: values
    integer :: ncid, varid, dimids(2), i
    character (len=2) :: dimc
    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    do i = 1, size(shape(values))
      write (dimc, "(I2.2)") i
      call check( nf90_def_dim(ncid, "dim" // dimc, size(values, i), dimids(i)) )
    end do
    call check( nf90_def_var(ncid, "values", NF90_DOUBLE, dimids, varid) )
    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, values) )
    call check( nf90_close(ncid) )
  end subroutine to_netcdf_real64_2d

end module netcdf_utils
