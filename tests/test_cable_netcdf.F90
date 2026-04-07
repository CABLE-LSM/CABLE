! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_netcdf
  !! Tests for cable_netcdf_mod.

  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_suite
  use fortuno_interface_mod, only: num_ranks
  use fortuno_interface_mod, only: this_rank
  use fortuno_interface_mod, only: check
  use fortuno_interface_mod, only: check_failed
  use fortuno_interface_mod, only: all_equal
  use fortuno_interface_mod, only: all_close

  use cable_netcdf_fixtures_mod, only: test_case_nf90
  use cable_netcdf_fixtures_mod, only: test_case_pio
  use cable_netcdf_fixtures_mod, only: io_handler_factory_interface

  use cable_netcdf_mod

  implicit none

  private

  public :: cable_netcdf_test_list

  integer, parameter :: LEN = 16
  integer, parameter :: VAL = 42

contains

  function cable_netcdf_test_list()
    type(test_list_t) :: cable_netcdf_test_list

    cable_netcdf_test_list = test_list_t([&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_int32_1d", test_write_read_darray_int32_1d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_int32_2d", test_write_read_darray_int32_2d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_int32_3d", test_write_read_darray_int32_3d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_real32_1d", test_write_read_darray_real32_1d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_real32_2d", test_write_read_darray_real32_2d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_real32_3d", test_write_read_darray_real32_3d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_real64_1d", test_write_read_darray_real64_1d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_real64_2d", test_write_read_darray_real64_2d),&
      test_case_nf90("test_cable_netcdf_nf90_write_read_darray_real64_3d", test_write_read_darray_real64_3d),&
      test_suite("parallel", test_list_t([&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_int32_1d", test_write_read_darray_int32_1d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_int32_2d", test_write_read_darray_int32_2d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_int32_3d", test_write_read_darray_int32_3d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_real32_1d", test_write_read_darray_real32_1d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_real32_2d", test_write_read_darray_real32_2d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_real32_3d", test_write_read_darray_real32_3d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_real64_1d", test_write_read_darray_real64_1d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_real64_2d", test_write_read_darray_real64_2d),&
        test_case_pio("test_cable_netcdf_pio_write_read_darray_real64_1d", test_write_read_darray_real64_3d)&
      ]))&
    ])

  end function cable_netcdf_test_list

  logical function init_decomp(compmap, start, end, block_per_pe) result(result)
    !* Intialise the `compmap` mapping, `start` and `end` indexes for the current
    ! rank, and the number of elements each rank will process (`block_per_pe`).
    !
    ! Returns `.true.` if any pre-checks fail, `.false.` otherwise.
    integer, allocatable, intent(out) :: compmap(:)
    integer, intent(out) :: start, end, block_per_pe
    integer i

    block_per_pe = LEN / num_ranks()

    result = .false.
    call check(mod(LEN, num_ranks()) == 0, msg="test_cable_netcdf.F90: work not divisible by number of ranks")
    call check(block_per_pe > 0, msg="test_cable_netcdf.F90: not enough work to distribute among pes")
    call check(mod(block_per_pe, 4) == 0, msg="test_cable_netcdf.F90: block_per_pe must be divisible by 4")
    if (check_failed()) then
      result = .true.
      return
    end if

    start = this_rank() * block_per_pe + 1
    end = start + block_per_pe - 1

    allocate(compmap(start:end))
    compmap(start:end) = [(i, i=start, end, 1)]

  end function init_decomp

  subroutine test_write_read_darray_int32_1d(io_handler_factory, file_name)
    !! Test writing and reading 1D integer arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: write_buffer(:), read_buffer(:)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(start:end), source=int(this_rank() + VAL, kind=CABLE_NETCDF_INT32_KIND))
    allocate(read_buffer(start:end), source=int(0, kind=CABLE_NETCDF_INT32_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_INT)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_INT, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_equal(write_buffer, read_buffer))

    call io_handler%finalise()

  end subroutine test_write_read_darray_int32_1d

  subroutine test_write_read_darray_int32_2d(io_handler_factory, file_name)
    !! Test writing and reading 2D integer arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: write_buffer(:, :), read_buffer(:, :)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(block_per_pe / 4, 4), source=int(this_rank() + VAL, kind=CABLE_NETCDF_INT32_KIND))
    allocate(read_buffer(block_per_pe / 4, 4), source=int(0, kind=CABLE_NETCDF_INT32_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_INT)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_INT, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_equal(reshape(write_buffer, [block_per_pe]), reshape(read_buffer, [block_per_pe])))

    call io_handler%finalise()

  end subroutine test_write_read_darray_int32_2d

  subroutine test_write_read_darray_int32_3d(io_handler_factory, file_name)
    !! Test writing and reading 3D integer arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: write_buffer(:, :, :), read_buffer(:, :, :)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(block_per_pe / 4, 2, 2), source=int(this_rank() + VAL, kind=CABLE_NETCDF_INT32_KIND))
    allocate(read_buffer(block_per_pe / 4, 2, 2), source=int(0, kind=CABLE_NETCDF_INT32_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_INT)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_INT, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_equal(reshape(write_buffer, [block_per_pe]), reshape(read_buffer, [block_per_pe])))

    call io_handler%finalise()

  end subroutine test_write_read_darray_int32_3d

  subroutine test_write_read_darray_real32_1d(io_handler_factory, file_name)
    !! Test writing and reading 1D 32-bit real arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: write_buffer(:), read_buffer(:)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(start:end), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL32_KIND))
    allocate(read_buffer(start:end), source=real(0, kind=CABLE_NETCDF_REAL32_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_FLOAT)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_FLOAT, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_close(write_buffer, read_buffer))

    call io_handler%finalise()

  end subroutine test_write_read_darray_real32_1d

  subroutine test_write_read_darray_real32_2d(io_handler_factory, file_name)
    !! Test writing and reading 2D 32-bit real arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: write_buffer(:, :), read_buffer(:, :)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(block_per_pe / 4, 4), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL32_KIND))
    allocate(read_buffer(block_per_pe / 4, 4), source=real(0, kind=CABLE_NETCDF_REAL32_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_FLOAT)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_FLOAT, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_close(reshape(write_buffer, [block_per_pe]), reshape(read_buffer, [block_per_pe])))

    call io_handler%finalise()

  end subroutine test_write_read_darray_real32_2d

  subroutine test_write_read_darray_real32_3d(io_handler_factory, file_name)
    !! Test writing and reading 3D 32-bit real arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: write_buffer(:, :, :), read_buffer(:, :, :)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(block_per_pe / 4, 2, 2), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL32_KIND))
    allocate(read_buffer(block_per_pe / 4, 2, 2), source=real(0, kind=CABLE_NETCDF_REAL32_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_FLOAT)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_FLOAT, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_close(reshape(write_buffer, [block_per_pe]), reshape(read_buffer, [block_per_pe])))

    call io_handler%finalise()

  end subroutine test_write_read_darray_real32_3d

  subroutine test_write_read_darray_real64_1d(io_handler_factory, file_name)
    !! Test writing and reading 1D 64-bit real arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: write_buffer(:), read_buffer(:)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(start:end), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL64_KIND))
    allocate(read_buffer(start:end), source=real(0, kind=CABLE_NETCDF_REAL64_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_DOUBLE)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_DOUBLE, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_close(write_buffer, read_buffer))

    call io_handler%finalise()

  end subroutine test_write_read_darray_real64_1d

  subroutine test_write_read_darray_real64_2d(io_handler_factory, file_name)
    !! Test writing and reading 2D 64-bit real arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: write_buffer(:, :), read_buffer(:, :)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(block_per_pe / 4, 4), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL64_KIND))
    allocate(read_buffer(block_per_pe / 4, 4), source=real(0, kind=CABLE_NETCDF_REAL64_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_DOUBLE)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_DOUBLE, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_close(reshape(write_buffer, [block_per_pe]), reshape(read_buffer, [block_per_pe])))

    call io_handler%finalise()

  end subroutine test_write_read_darray_real64_2d

  subroutine test_write_read_darray_real64_3d(io_handler_factory, file_name)
    !! Test writing and reading 3D 64-bit real arrays using the darray API.
    procedure(io_handler_factory_interface) :: io_handler_factory
      !! Factory procedure to create an IO handler (e.g., NetCDF or PIO)
    character(*), intent(in) :: file_name
      !! Name of the file to create and read from during the test
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: write_buffer(:, :, :), read_buffer(:, :, :)

    if (init_decomp(compmap, start, end, block_per_pe)) return

    allocate(write_buffer(block_per_pe / 4, 2, 2), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL64_KIND))
    allocate(read_buffer(block_per_pe / 4, 2, 2), source=real(0, kind=CABLE_NETCDF_REAL64_KIND))

    io_handler = io_handler_factory()

    call io_handler%init()

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_DOUBLE)

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    call file%def_dims(["i"], [LEN])
    call file%def_var("values", CABLE_NETCDF_DOUBLE, ["i"])
    call file%end_def()
    call file%write_darray("values", write_buffer, decomp)
    call file%sync()
    call file%read_darray("values", read_buffer, decomp)
    call file%close()

    call check(all_close(reshape(write_buffer, [block_per_pe]), reshape(read_buffer, [block_per_pe])))

    call io_handler%finalise()

  end subroutine test_write_read_darray_real64_3d

end module test_cable_netcdf
