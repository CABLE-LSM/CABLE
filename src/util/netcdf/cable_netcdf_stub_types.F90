module cable_netcdf_stub_types_mod
  use cable_netcdf_mod
  use cable_mpi_mod, only: mpi_grp_t
  use iso_fortran_env, only: error_unit
  implicit none

  private
  public :: &
    cable_netcdf_stub_decomp_t, &
    cable_netcdf_stub_file_t, &
    cable_netcdf_stub_io_t

  type, extends(cable_netcdf_decomp_t) :: cable_netcdf_stub_decomp_t
  end type

  type, extends(cable_netcdf_io_t) :: cable_netcdf_stub_io_t
  contains
    procedure :: init => cable_netcdf_stub_io_init
    procedure :: finalise => cable_netcdf_stub_io_finalise
    procedure :: create_file => cable_netcdf_stub_io_create_file
    procedure :: open_file => cable_netcdf_stub_io_open_file
    procedure :: create_decomp => cable_netcdf_stub_io_create_decomp
  end type

  type, extends(cable_netcdf_file_t) :: cable_netcdf_stub_file_t
  contains
    procedure :: close => cable_netcdf_stub_file_close
    procedure :: end_def => cable_netcdf_stub_file_end_def
    procedure :: sync => cable_netcdf_stub_file_sync
    procedure :: def_dims => cable_netcdf_stub_file_def_dims
    procedure :: def_var => cable_netcdf_stub_file_def_var
    procedure :: put_att_global_string => cable_netcdf_stub_file_put_att_global_string
    procedure :: put_att_global_int32 => cable_netcdf_stub_file_put_att_global_int32
    procedure :: put_att_global_real32 => cable_netcdf_stub_file_put_att_global_real32
    procedure :: put_att_global_real64 => cable_netcdf_stub_file_put_att_global_real64
    procedure :: put_att_var_string => cable_netcdf_stub_file_put_att_var_string
    procedure :: put_att_var_int32 => cable_netcdf_stub_file_put_att_var_int32
    procedure :: put_att_var_real32 => cable_netcdf_stub_file_put_att_var_real32
    procedure :: put_att_var_real64 => cable_netcdf_stub_file_put_att_var_real64
    procedure :: get_att_global_string => cable_netcdf_stub_file_get_att_global_string
    procedure :: get_att_global_int32 => cable_netcdf_stub_file_get_att_global_int32
    procedure :: get_att_global_real32 => cable_netcdf_stub_file_get_att_global_real32
    procedure :: get_att_global_real64 => cable_netcdf_stub_file_get_att_global_real64
    procedure :: get_att_var_string => cable_netcdf_stub_file_get_att_var_string
    procedure :: get_att_var_int32 => cable_netcdf_stub_file_get_att_var_int32
    procedure :: get_att_var_real32 => cable_netcdf_stub_file_get_att_var_real32
    procedure :: get_att_var_real64 => cable_netcdf_stub_file_get_att_var_real64
    procedure :: inq_dim_len => cable_netcdf_stub_file_inq_dim_len
    procedure :: put_var_int32_0d => cable_netcdf_stub_file_put_var_int32_0d
    procedure :: put_var_int32_1d => cable_netcdf_stub_file_put_var_int32_1d
    procedure :: put_var_int32_2d => cable_netcdf_stub_file_put_var_int32_2d
    procedure :: put_var_int32_3d => cable_netcdf_stub_file_put_var_int32_3d
    procedure :: put_var_real32_0d => cable_netcdf_stub_file_put_var_real32_0d
    procedure :: put_var_real32_1d => cable_netcdf_stub_file_put_var_real32_1d
    procedure :: put_var_real32_2d => cable_netcdf_stub_file_put_var_real32_2d
    procedure :: put_var_real32_3d => cable_netcdf_stub_file_put_var_real32_3d
    procedure :: put_var_real64_0d => cable_netcdf_stub_file_put_var_real64_0d
    procedure :: put_var_real64_1d => cable_netcdf_stub_file_put_var_real64_1d
    procedure :: put_var_real64_2d => cable_netcdf_stub_file_put_var_real64_2d
    procedure :: put_var_real64_3d => cable_netcdf_stub_file_put_var_real64_3d
    procedure :: write_darray_int32_1d => cable_netcdf_stub_file_write_darray_int32_1d
    procedure :: write_darray_int32_2d => cable_netcdf_stub_file_write_darray_int32_2d
    procedure :: write_darray_int32_3d => cable_netcdf_stub_file_write_darray_int32_3d
    procedure :: write_darray_real32_1d => cable_netcdf_stub_file_write_darray_real32_1d
    procedure :: write_darray_real32_2d => cable_netcdf_stub_file_write_darray_real32_2d
    procedure :: write_darray_real32_3d => cable_netcdf_stub_file_write_darray_real32_3d
    procedure :: write_darray_real64_1d => cable_netcdf_stub_file_write_darray_real64_1d
    procedure :: write_darray_real64_2d => cable_netcdf_stub_file_write_darray_real64_2d
    procedure :: write_darray_real64_3d => cable_netcdf_stub_file_write_darray_real64_3d
    procedure :: get_var_int32_0d => cable_netcdf_stub_file_get_var_int32_0d
    procedure :: get_var_int32_1d => cable_netcdf_stub_file_get_var_int32_1d
    procedure :: get_var_int32_2d => cable_netcdf_stub_file_get_var_int32_2d
    procedure :: get_var_int32_3d => cable_netcdf_stub_file_get_var_int32_3d
    procedure :: get_var_real32_0d => cable_netcdf_stub_file_get_var_real32_0d
    procedure :: get_var_real32_1d => cable_netcdf_stub_file_get_var_real32_1d
    procedure :: get_var_real32_2d => cable_netcdf_stub_file_get_var_real32_2d
    procedure :: get_var_real32_3d => cable_netcdf_stub_file_get_var_real32_3d
    procedure :: get_var_real64_0d => cable_netcdf_stub_file_get_var_real64_0d
    procedure :: get_var_real64_1d => cable_netcdf_stub_file_get_var_real64_1d
    procedure :: get_var_real64_2d => cable_netcdf_stub_file_get_var_real64_2d
    procedure :: get_var_real64_3d => cable_netcdf_stub_file_get_var_real64_3d
    procedure :: read_darray_int32_1d => cable_netcdf_stub_file_read_darray_int32_1d
    procedure :: read_darray_int32_2d => cable_netcdf_stub_file_read_darray_int32_2d
    procedure :: read_darray_int32_3d => cable_netcdf_stub_file_read_darray_int32_3d
    procedure :: read_darray_real32_1d => cable_netcdf_stub_file_read_darray_real32_1d
    procedure :: read_darray_real32_2d => cable_netcdf_stub_file_read_darray_real32_2d
    procedure :: read_darray_real32_3d => cable_netcdf_stub_file_read_darray_real32_3d
    procedure :: read_darray_real64_1d => cable_netcdf_stub_file_read_darray_real64_1d
    procedure :: read_darray_real64_2d => cable_netcdf_stub_file_read_darray_real64_2d
    procedure :: read_darray_real64_3d => cable_netcdf_stub_file_read_darray_real64_3d
  end type

contains

  subroutine cable_netcdf_stub_io_init(this)
    class(cable_netcdf_stub_io_t), intent(inout) :: this
  end subroutine

  subroutine cable_netcdf_stub_io_finalise(this)
    class(cable_netcdf_stub_io_t), intent(inout) :: this
  end subroutine

  function cable_netcdf_stub_io_create_file(this, path) result(file)
    class(cable_netcdf_stub_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    class(cable_netcdf_file_t), allocatable :: file
    file = cable_netcdf_stub_file_t()
  end function

  function cable_netcdf_stub_io_open_file(this, path) result(file)
    class(cable_netcdf_stub_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    class(cable_netcdf_file_t), allocatable :: file
    file = cable_netcdf_stub_file_t()
  end function

  function cable_netcdf_stub_io_create_decomp(this, compmap, dims, type) result(decomp)
    class(cable_netcdf_stub_io_t), intent(inout) :: this
    integer, intent(in) :: compmap(:), dims(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp
    decomp = cable_netcdf_stub_decomp_t(compmap, dims, type)
  end function

  subroutine cable_netcdf_stub_file_close(this)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
  end subroutine

  subroutine cable_netcdf_stub_file_end_def(this)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
  end subroutine

  subroutine cable_netcdf_stub_file_sync(this)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
  end subroutine

  subroutine cable_netcdf_stub_file_def_dims(this, dim_names, dim_lens)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: dim_names(:)
    integer, intent(in) :: dim_lens(:)
  end subroutine

  subroutine cable_netcdf_stub_file_def_var(this, var_name, dim_names, type)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, dim_names(:)
    integer, intent(in) :: type
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_global_string(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name, att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_global_int32(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_global_real32(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_global_real64(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_var_string(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name, att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_var_int32(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_var_real32(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
  end subroutine

  subroutine cable_netcdf_stub_file_put_att_var_real64(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_global_string(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    character(len=*), intent(out) :: att_value
    att_value = ""
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_global_int32(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    att_value = 0
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_global_real32(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    att_value = 0.0
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_global_real64(this, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    att_value = 0.0
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_var_string(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    character(len=*), intent(out) :: att_value
    att_value = ""
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_var_int32(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    att_value = 0
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_var_real32(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    att_value = 0.0
  end subroutine

  subroutine cable_netcdf_stub_file_get_att_var_real64(this, var_name, att_name, att_value)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    att_value = 0.0
  end subroutine

  subroutine cable_netcdf_stub_file_inq_dim_len(this, dim_name, dim_len)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: dim_name
    integer, intent(out) :: dim_len
    dim_len = 0
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_int32_0d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_int32_1d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_int32_2d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_int32_3d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real32_0d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real32_1d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real32_2d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real32_3d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real64_0d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real64_1d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real64_2d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_put_var_real64_3d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_int32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_int32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_int32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_real32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_real32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_real32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_real64_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_real64_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_write_darray_real64_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_int32_0d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    values = 0
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_int32_1d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    values = 0
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_int32_2d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    values = 0
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_int32_3d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    values = 0
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real32_0d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real32_1d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real32_2d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real32_3d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real64_0d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real64_1d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real64_2d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_get_var_real64_3d(this, var_name, values, start, count)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_int32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_int32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_int32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_real32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_real32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_real32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_real64_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_real64_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

  subroutine cable_netcdf_stub_file_read_darray_real64_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_stub_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    values = 0.
  end subroutine

end module cable_netcdf_stub_types_mod
