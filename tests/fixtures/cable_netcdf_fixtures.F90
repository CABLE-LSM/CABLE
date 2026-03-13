module cable_netcdf_fixtures_mod
  !! Fixtures for testing cable_netcdf_mod.

  use fortuno_interface_mod, only: test_case_t
  use fortuno_interface_mod, only: test_item_t
  use fortuno_interface_mod, only: skip
  use fortuno_interface_mod, only: check
  use fortuno_interface_mod, only: check_failed
  use fortuno_interface_mod, only: global_comm
  use fortuno_interface_mod, only: num_ranks
  use fortuno_interface_mod, only: this_rank

  use cable_mpi_mod, only: mpi_grp_t

  use file_utils, only: file_delete

  use cable_netcdf_mod, only: cable_netcdf_io_t
  use cable_netcdf_nf90_mod, only: cable_netcdf_nf90_io_t
  use cable_netcdf_pio_mod, only: cable_netcdf_pio_io_t

  implicit none

  private

  public :: test_case_nf90
  public :: test_case_pio
  public :: io_handler_factory_interface

  character(16), parameter :: nc_file_name = "file.nc"

  type, extends(test_case_t) :: test_case_cable_netcdf_nf90_t
    procedure(cable_netcdf_test_interface), pointer, nopass :: test
  contains
    procedure :: run => run_test_cable_netcdf_nf90
  end type test_case_cable_netcdf_nf90_t

  type, extends(test_case_t) :: test_case_cable_netcdf_pio_t
    procedure(cable_netcdf_test_interface), pointer, nopass :: test
  contains
    procedure :: run => run_test_pio
  end type test_case_cable_netcdf_pio_t

  abstract interface
    function io_handler_factory_interface() result(io_handler)
      import cable_netcdf_io_t
      class(cable_netcdf_io_t), allocatable :: io_handler
    end function
    subroutine cable_netcdf_test_interface(io_handler_factory, file_name)
      import io_handler_factory_interface
      procedure(io_handler_factory_interface) :: io_handler_factory
      character(*), intent(in) :: file_name
    end subroutine
  end interface

contains

  function test_case_nf90(name, test) result(testitem)
    character(*), intent(in) :: name
    procedure(cable_netcdf_test_interface) :: test
    type(test_item_t) :: testitem
    testitem = test_item_t(test_case_cable_netcdf_nf90_t(name=name, test=test))
  end function test_case_nf90

  function nf90_io_handler() result(io_handler)
    class(cable_netcdf_io_t), allocatable :: io_handler
    io_handler = cable_netcdf_nf90_io_t()
  end function

  subroutine run_test_cable_netcdf_nf90(this)
    class(test_case_cable_netcdf_nf90_t), intent(in) :: this

    call check(num_ranks() == 1, msg="NetCDF tests must be run with a single process.")
    if (check_failed()) return

    call this%test(nf90_io_handler, nc_file_name)
    call file_delete(nc_file_name)

  end subroutine run_test_cable_netcdf_nf90

  function test_case_pio(name, test) result(testitem)
    character(*), intent(in) :: name
    procedure(cable_netcdf_test_interface) :: test
    type(test_item_t) :: testitem
    testitem = test_item_t(test_case_cable_netcdf_pio_t(name=name, test=test))
  end function test_case_pio

  function pio_io_handler() result(io_handler)
    class(cable_netcdf_io_t), allocatable :: io_handler
    io_handler = cable_netcdf_pio_io_t(mpi_grp_t(global_comm()))
  end function

  subroutine run_test_pio(this)
    class(test_case_cable_netcdf_pio_t), intent(in) :: this
    type(mpi_grp_t) :: mpi_grp

#ifdef SKIP_PIO_TESTS
    call skip()
    return
#endif

    mpi_grp = mpi_grp_t(global_comm())
    call check(mpi_grp%comm_defined(), msg="MPI communicator must be defined for PIO tests")
    if (check_failed()) return

    call this%test(pio_io_handler, nc_file_name)
    call file_delete(nc_file_name)

  end subroutine run_test_pio

end module cable_netcdf_fixtures_mod
