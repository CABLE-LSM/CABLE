module test_cable_netcdf
  use fortuno_interface_m, only: check, check_failed, suite, test_list, num_ranks, this_rank, all_equal, all_close
  use fixtures_mod, only: test_case_nf90, test_case_pio, io_handler_factory_interface
  use cable_netcdf_mod
  use cable_netcdf_nf90_mod
  use cable_netcdf_pio_mod
  implicit none

  private

  public :: cable_netcdf_test_list

  integer, parameter :: LEN = 16
  integer, parameter :: VAL = 42

contains

  function cable_netcdf_test_list()
    type(test_list) :: cable_netcdf_test_list

    cable_netcdf_test_list = test_list([&
      test_case_nf90("test_cable_netcdf_nf90_write_darray_int32", test_write_read_darray_int32),&
      test_case_nf90("test_cable_netcdf_nf90_write_darray_real32", test_write_read_darray_real32),&
      test_case_nf90("test_cable_netcdf_nf90_write_darray_real64", test_write_read_darray_real64),&
      suite("parallel", test_list([&
        test_case_pio("test_cable_netcdf_pio_write_darray_int32", test_write_read_darray_int32),&
        test_case_pio("test_cable_netcdf_pio_write_darray_real32", test_write_read_darray_real32),&
        test_case_pio("test_cable_netcdf_pio_write_darray_real64", test_write_read_darray_real64)&
      ]))&
    ])

  end function cable_netcdf_test_list

  logical function check_valid_decomp(compmap, start, end, block_per_pe)
    integer, allocatable, intent(out) :: compmap(:)
    integer, intent(out) :: start, end, block_per_pe
    integer i

    block_per_pe = LEN / num_ranks()

    call check(mod(LEN, num_ranks()) == 0, msg="test_cable_netcdf.F90: work not divisible by number of ranks")
    call check(block_per_pe > 0, msg="test_cable_netcdf.F90: not enough work to distribute among pes")
    call check(mod(block_per_pe, 4) == 0, msg="test_cable_netcdf.F90: block_per_pe must be divisible by 4")

    check_valid_decomp = .not. check_failed()
    if (check_failed()) return

    start = this_rank() * block_per_pe + 1
    end = start + block_per_pe - 1

    allocate(compmap(start:end))
    compmap(start:end) = [(i, i=start, end, 1)]

  end function check_valid_decomp

  subroutine test_write_read_darray_int32(io_handler_factory, file_name)
    procedure(io_handler_factory_interface), pointer, intent(in) :: io_handler_factory
    character(*), intent(in) :: file_name
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end, buffer_shape_2d(2), buffer_shape_3d(3)
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: write_buffer_1d(:), write_buffer_2d(:, :), write_buffer_3d(:, :, :)
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: read_buffer_1d(:), read_buffer_2d(:, :), read_buffer_3d(:, :, :)

    if (.not. check_valid_decomp(compmap, start, end, block_per_pe)) return

    buffer_shape_2d = [block_per_pe / 4, 4]
    buffer_shape_3d = [block_per_pe / 4, 2, 2]

    io_handler = io_handler_factory()

    call io_handler%init()

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    allocate(write_buffer_1d(start:end), source=int(this_rank() + VAL, kind=CABLE_NETCDF_INT32_KIND))
    write_buffer_2d = reshape(write_buffer_1d, buffer_shape_2d)
    write_buffer_3d = reshape(write_buffer_1d, buffer_shape_3d)

    allocate(read_buffer_1d(start:end), source=int(0, kind=CABLE_NETCDF_INT32_KIND))
    read_buffer_2d = reshape(read_buffer_1d, buffer_shape_2d)
    read_buffer_3d = reshape(read_buffer_1d, buffer_shape_3d)

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_INT)

    call file%def_dims(["i"], [LEN])

    call file%def_var("values_1d", ["i"], CABLE_NETCDF_INT)
    call file%def_var("values_2d", ["i"], CABLE_NETCDF_INT)
    call file%def_var("values_3d", ["i"], CABLE_NETCDF_INT)

    call file%end_def()

    call file%write_darray("values_1d", write_buffer_1d, decomp)
    call file%write_darray("values_2d", write_buffer_2d, decomp)
    call file%write_darray("values_3d", write_buffer_3d, decomp)

    call file%sync()

    call file%read_darray("values_1d", read_buffer_1d, decomp)
    call file%read_darray("values_2d", read_buffer_2d, decomp)
    call file%read_darray("values_3d", read_buffer_3d, decomp)

    call check(all_equal(write_buffer_1d, read_buffer_1d))
    call check(all_equal(reshape(write_buffer_2d, [block_per_pe]), reshape(read_buffer_2d, [block_per_pe])))
    call check(all_equal(reshape(write_buffer_3d, [block_per_pe]), reshape(read_buffer_3d, [block_per_pe])))

    call file%close()

    call io_handler%finalise()

  end subroutine test_write_read_darray_int32

  subroutine test_write_read_darray_real32(io_handler_factory, file_name)
    procedure(io_handler_factory_interface), pointer, intent(in) :: io_handler_factory
    character(*), intent(in) :: file_name
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end, buffer_shape_2d(2), buffer_shape_3d(3)
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: write_buffer_1d(:), write_buffer_2d(:, :), write_buffer_3d(:, :, :)
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: read_buffer_1d(:), read_buffer_2d(:, :), read_buffer_3d(:, :, :)

    if (.not. check_valid_decomp(compmap, start, end, block_per_pe)) return

    buffer_shape_2d = [block_per_pe / 4, 4]
    buffer_shape_3d = [block_per_pe / 4, 2, 2]

    io_handler = io_handler_factory()

    call io_handler%init()

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    allocate(write_buffer_1d(start:end), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL32_KIND))
    write_buffer_2d = reshape(write_buffer_1d, buffer_shape_2d)
    write_buffer_3d = reshape(write_buffer_1d, buffer_shape_3d)

    allocate(read_buffer_1d(start:end), source=real(0, kind=CABLE_NETCDF_REAL32_KIND))
    read_buffer_2d = reshape(read_buffer_1d, buffer_shape_2d)
    read_buffer_3d = reshape(read_buffer_1d, buffer_shape_3d)

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_FLOAT)

    call file%def_dims(["i"], [LEN])

    call file%def_var("values_1d", ["i"], CABLE_NETCDF_FLOAT)
    call file%def_var("values_2d", ["i"], CABLE_NETCDF_FLOAT)
    call file%def_var("values_3d", ["i"], CABLE_NETCDF_FLOAT)

    call file%end_def()

    call file%write_darray("values_1d", write_buffer_1d, decomp)
    call file%write_darray("values_2d", write_buffer_2d, decomp)
    call file%write_darray("values_3d", write_buffer_3d, decomp)

    call file%sync()

    call file%read_darray("values_1d", read_buffer_1d, decomp)
    call file%read_darray("values_2d", read_buffer_2d, decomp)
    call file%read_darray("values_3d", read_buffer_3d, decomp)

    call check(all_close(write_buffer_1d, read_buffer_1d))
    call check(all_close(reshape(write_buffer_2d, [block_per_pe]), reshape(read_buffer_2d, [block_per_pe])))
    call check(all_close(reshape(write_buffer_3d, [block_per_pe]), reshape(read_buffer_3d, [block_per_pe])))

    call file%close()

    call io_handler%finalise()

  end subroutine test_write_read_darray_real32

  subroutine test_write_read_darray_real64(io_handler_factory, file_name)
    procedure(io_handler_factory_interface), pointer, intent(in) :: io_handler_factory
    character(*), intent(in) :: file_name
    class(cable_netcdf_io_t), allocatable :: io_handler
    class(cable_netcdf_file_t), allocatable :: file
    class(cable_netcdf_decomp_t), allocatable :: decomp
    integer, allocatable :: compmap(:)
    integer :: block_per_pe, start, end, buffer_shape_2d(2), buffer_shape_3d(3)
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: write_buffer_1d(:), write_buffer_2d(:, :), write_buffer_3d(:, :, :)
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: read_buffer_1d(:), read_buffer_2d(:, :), read_buffer_3d(:, :, :)

    if (.not. check_valid_decomp(compmap, start, end, block_per_pe)) return

    buffer_shape_2d = [block_per_pe / 4, 4]
    buffer_shape_3d = [block_per_pe / 4, 2, 2]

    io_handler = io_handler_factory()

    call io_handler%init()

    file = io_handler%create_file(file_name, iotype=CABLE_NETCDF_IOTYPE_NETCDF4P)

    allocate(write_buffer_1d(start:end), source=real(this_rank() + VAL, kind=CABLE_NETCDF_REAL64_KIND))
    write_buffer_2d = reshape(write_buffer_1d, buffer_shape_2d)
    write_buffer_3d = reshape(write_buffer_1d, buffer_shape_3d)

    allocate(read_buffer_1d(start:end), source=real(0, kind=CABLE_NETCDF_REAL64_KIND))
    read_buffer_2d = reshape(read_buffer_1d, buffer_shape_2d)
    read_buffer_3d = reshape(read_buffer_1d, buffer_shape_3d)

    decomp = io_handler%create_decomp(compmap, dims=[LEN], type=CABLE_NETCDF_DOUBLE)

    call file%def_dims(["i"], [LEN])

    call file%def_var("values_1d", ["i"], CABLE_NETCDF_DOUBLE)
    call file%def_var("values_2d", ["i"], CABLE_NETCDF_DOUBLE)
    call file%def_var("values_3d", ["i"], CABLE_NETCDF_DOUBLE)

    call file%end_def()

    call file%write_darray("values_1d", write_buffer_1d, decomp)
    call file%write_darray("values_2d", write_buffer_2d, decomp)
    call file%write_darray("values_3d", write_buffer_3d, decomp)

    call file%sync()

    call file%read_darray("values_1d", read_buffer_1d, decomp)
    call file%read_darray("values_2d", read_buffer_2d, decomp)
    call file%read_darray("values_3d", read_buffer_3d, decomp)

    call check(all_close(write_buffer_1d, read_buffer_1d))
    call check(all_close(reshape(write_buffer_2d, [block_per_pe]), reshape(read_buffer_2d, [block_per_pe])))
    call check(all_close(reshape(write_buffer_3d, [block_per_pe]), reshape(read_buffer_3d, [block_per_pe])))

    call file%close()

    call io_handler%finalise()

  end subroutine test_write_read_darray_real64

end module test_cable_netcdf
