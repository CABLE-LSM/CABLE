module fixtures_mod
  use fortuno_interface_m, only: test_case_t, test_item, check, check_failed, global_comm, num_ranks, this_rank
  use cable_mpi_mod, only: mpi_grp_t, MPI_COMM_UNDEFINED
  use file_utils, only: file_delete_collective
  use cable_netcdf_mod, only: cable_netcdf_io_t, CABLE_NETCDF_MAX_STR_LEN_FILE
  use cable_netcdf_nf90_mod, only: cable_netcdf_nf90_io_t
  use cable_netcdf_pio_mod, only: cable_netcdf_pio_io_t
  implicit none

  private

  public :: test_case_nf90, test_case_pio, io_handler_factory_interface

  character(len=CABLE_NETCDF_MAX_STR_LEN_FILE), parameter :: nc_file_name = "file.nc"

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
      import cable_netcdf_io_t, io_handler_factory_interface
      procedure(io_handler_factory_interface), pointer, intent(in) :: io_handler_factory
      character(*), intent(in) :: file_name
    end subroutine
  end interface

contains

  function test_case_nf90(name, test) result(testitem)
    character(*), intent(in) :: name
    procedure(cable_netcdf_test_interface) :: test
    type(test_item) :: testitem
    testitem = test_item(test_case_cable_netcdf_nf90_t(name=name, test=test))
  end function test_case_nf90

  function nf90_io_handler() result(io_handler)
    class(cable_netcdf_io_t), allocatable :: io_handler
    io_handler = cable_netcdf_nf90_io_t()
  end function

  subroutine run_test_cable_netcdf_nf90(this)
    class(test_case_cable_netcdf_nf90_t), intent(in) :: this
    class(cable_netcdf_io_t), allocatable :: io_handler

    call check(num_ranks() == 1, msg="NetCDF tests must be run with a single process.")
    if (check_failed()) return

    call this%test(nf90_io_handler, nc_file_name)
    call file_delete_collective(nc_file_name)

  end subroutine run_test_cable_netcdf_nf90

  function test_case_pio(name, test) result(testitem)
    character(*), intent(in) :: name
    procedure(cable_netcdf_test_interface) :: test
    type(test_item) :: testitem
    testitem = test_item(test_case_cable_netcdf_pio_t(name=name, test=test))
  end function test_case_pio

  function pio_io_handler() result(io_handler)
    class(cable_netcdf_io_t), allocatable :: io_handler
    io_handler = cable_netcdf_pio_io_t(mpi_grp_t(global_comm()))
  end function

  subroutine run_test_pio(this)
    class(test_case_cable_netcdf_pio_t), intent(in) :: this
    class(cable_netcdf_io_t), allocatable :: io_handler
    type(mpi_grp_t) :: mpi_grp

    mpi_grp = mpi_grp_t(global_comm())
    call check(&
      mpi_grp%comm%mpi_val /= MPI_COMM_UNDEFINED%mpi_val,&
      msg="MPI communicator must be defined for PIO tests"&
    )
    if (check_failed()) return

    call this%test(pio_io_handler, nc_file_name)
    call file_delete_collective(nc_file_name)

  end subroutine run_test_pio

end module fixtures_mod
