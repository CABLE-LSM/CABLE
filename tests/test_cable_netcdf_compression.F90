! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_netcdf_compression
  !! Tests for the shuffle and deflate options of `def_var` (nf90 implementation).
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_inq_var_deflate, nf90_noerr, nf90_nowrite
  use cable_netcdf_mod, only: cable_netcdf_io_t
  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_NETCDF4C
  use cable_netcdf_nf90_mod, only: cable_netcdf_nf90_io_t
  implicit none
  private

  public :: cable_netcdf_compression_test_list

  character(len=*), parameter :: file_name = "compression_test.nc"

contains

  function cable_netcdf_compression_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_netcdf_compression_settings_are_applied", test_settings_are_applied), &
      test_case("test_netcdf_no_compression_by_default", test_no_compression_by_default), &
      test_case("test_netcdf_compression_ignored_for_scalars", test_ignored_for_scalars) &
    ])
  end function cable_netcdf_compression_test_list

  subroutine write_test_file()
    !! Three variables: compressed with shuffle, uncompressed, and a compressed scalar.
    type(cable_netcdf_nf90_io_t) :: io_handler
    class(cable_netcdf_file_t), allocatable :: file

    call io_handler%init()
    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4C)
    call file%def_dims(["i"], [16])
    call file%def_var("compressed", CABLE_NETCDF_FLOAT, ["i"], shuffle=.true., deflate_level=4)
    call file%def_var("plain", CABLE_NETCDF_FLOAT, ["i"])
    call file%def_var("scalar", CABLE_NETCDF_FLOAT, shuffle=.true., deflate_level=4)
    call file%end_def()
    call file%close()
    call io_handler%finalise()
  end subroutine write_test_file

  subroutine read_deflate_settings(variable_name, shuffle, deflate, level)
    character(len=*), intent(in) :: variable_name
    integer, intent(out) :: shuffle, deflate, level
    integer :: ncid, varid

    call check(nf90_open(file_name, nf90_nowrite, ncid) == nf90_noerr, msg="open test file")
    call check(nf90_inq_varid(ncid, variable_name, varid) == nf90_noerr, msg="find variable " // variable_name)
    call check(nf90_inq_var_deflate(ncid, varid, shuffle, deflate, level) == nf90_noerr, msg="inquire deflate")
    call check(nf90_close(ncid) == nf90_noerr, msg="close test file")
    call delete_test_file()
  end subroutine read_deflate_settings

  subroutine delete_test_file()
    integer :: file_unit

    open(newunit=file_unit, file=file_name, status="old")
    close(file_unit, status="delete")
  end subroutine delete_test_file

  subroutine test_settings_are_applied()
    integer :: shuffle, deflate, level

    call write_test_file()
    call read_deflate_settings("compressed", shuffle, deflate, level)
    call check(shuffle == 1, msg="shuffle should be enabled")
    call check(deflate == 1, msg="deflate should be enabled")
    call check(level == 4, msg="deflate level should be 4")
  end subroutine test_settings_are_applied

  subroutine test_no_compression_by_default()
    integer :: shuffle, deflate, level

    call write_test_file()
    call read_deflate_settings("plain", shuffle, deflate, level)
    call check(shuffle == 0, msg="shuffle should be off")
    call check(deflate == 0, msg="deflate should be off")
  end subroutine test_no_compression_by_default

  subroutine test_ignored_for_scalars()
    !! Writing the file in `write_test_file` would abort if compression were
    !! applied to the scalar, so reaching here shows it was ignored.
    integer :: shuffle, deflate, level

    call write_test_file()
    call read_deflate_settings("scalar", shuffle, deflate, level)
    call check(deflate == 0, msg="scalars are not compressed")
  end subroutine test_ignored_for_scalars

end module test_cable_netcdf_compression
