module cable_netcdf_mod
  use iso_fortran_env, only: CABLE_NETCDF_INT32_KIND => int32
  use iso_fortran_env, only: CABLE_NETCDF_REAL32_KIND => real32
  use iso_fortran_env, only: CABLE_NETCDF_REAL64_KIND => real64
  use cable_mpi_mod, only: mpi_grp_t
  implicit none

  private

  public :: &
    cable_netcdf_decomp_t, &
    cable_netcdf_file_t, &
    cable_netcdf_io_t

  public :: &
    cable_netcdf_mod_init, &
    cable_netcdf_mod_end, &
    cable_netcdf_create_file, &
    cable_netcdf_open_file, &
    cable_netcdf_create_decomp

  public :: &
    CABLE_NETCDF_INT32_KIND, &
    CABLE_NETCDF_REAL32_KIND, &
    CABLE_NETCDF_REAL64_KIND, &
    CABLE_NETCDF_INT, &
    CABLE_NETCDF_FLOAT, &
    CABLE_NETCDF_DOUBLE, &
    CABLE_NETCDF_MAX_STR_LEN_FILE, &
    CABLE_NETCDF_MAX_STR_LEN_VAR, &
    CABLE_NETCDF_MAX_STR_LEN_DIM, &
    CABLE_NETCDF_MAX_RANK, &
    CABLE_NETCDF_UNLIMITED

  enum, bind(c)
    enumerator :: &
      CABLE_NETCDF_INT, &
      CABLE_NETCDF_FLOAT, &
      CABLE_NETCDF_DOUBLE
  end enum

  integer, parameter :: CABLE_NETCDF_MAX_STR_LEN_FILE = 200
  integer, parameter :: CABLE_NETCDF_MAX_STR_LEN_VAR = 80
  integer, parameter :: CABLE_NETCDF_MAX_STR_LEN_DIM = 20
  integer, parameter :: CABLE_NETCDF_UNLIMITED = -1
  integer, parameter :: CABLE_NETCDF_MAX_RANK = 3

  type :: cable_netcdf_decomp_t
    integer, allocatable :: compmap(:)
    integer, allocatable :: dims(:)
    integer :: type
  end type

  type, abstract :: cable_netcdf_file_t
  contains
    procedure(cable_netcdf_file_close), deferred :: close
    procedure(cable_netcdf_file_end_def), deferred :: end_def
    procedure(cable_netcdf_file_redef), deferred :: redef
    procedure(cable_netcdf_file_sync), deferred :: sync
    procedure(cable_netcdf_file_def_dims), deferred :: def_dims
    procedure(cable_netcdf_file_def_var), deferred :: def_var
    procedure(cable_netcdf_file_put_att_global_string), deferred :: put_att_global_string
    procedure(cable_netcdf_file_put_att_global_int32), deferred :: put_att_global_int32
    procedure(cable_netcdf_file_put_att_global_real32), deferred :: put_att_global_real32
    procedure(cable_netcdf_file_put_att_global_real64), deferred :: put_att_global_real64
    procedure(cable_netcdf_file_put_att_var_string), deferred :: put_att_var_string
    procedure(cable_netcdf_file_put_att_var_int32), deferred :: put_att_var_int32
    procedure(cable_netcdf_file_put_att_var_real32), deferred :: put_att_var_real32
    procedure(cable_netcdf_file_put_att_var_real64), deferred :: put_att_var_real64
    generic :: put_att => &
      put_att_global_string, put_att_global_int32, put_att_global_real32, put_att_global_real64, &
      put_att_var_string, put_att_var_int32, put_att_var_real32, put_att_var_real64
    procedure(cable_netcdf_file_get_att_global_string), deferred :: get_att_global_string
    procedure(cable_netcdf_file_get_att_global_int32), deferred :: get_att_global_int32
    procedure(cable_netcdf_file_get_att_global_real32), deferred :: get_att_global_real32
    procedure(cable_netcdf_file_get_att_global_real64), deferred :: get_att_global_real64
    procedure(cable_netcdf_file_get_att_var_string), deferred :: get_att_var_string
    procedure(cable_netcdf_file_get_att_var_int32), deferred :: get_att_var_int32
    procedure(cable_netcdf_file_get_att_var_real32), deferred :: get_att_var_real32
    procedure(cable_netcdf_file_get_att_var_real64), deferred :: get_att_var_real64
    generic :: get_att => &
      get_att_global_string, get_att_global_int32, get_att_global_real32, get_att_global_real64, &
      get_att_var_string, get_att_var_int32, get_att_var_real32, get_att_var_real64
    procedure(cable_netcdf_file_inq_dim_len), deferred :: inq_dim_len
    procedure(cable_netcdf_file_put_var_int32_0d), deferred :: put_var_int32_0d
    procedure(cable_netcdf_file_put_var_int32_1d), deferred :: put_var_int32_1d
    procedure(cable_netcdf_file_put_var_int32_2d), deferred :: put_var_int32_2d
    procedure(cable_netcdf_file_put_var_int32_3d), deferred :: put_var_int32_3d
    procedure(cable_netcdf_file_put_var_real32_0d), deferred :: put_var_real32_0d
    procedure(cable_netcdf_file_put_var_real32_1d), deferred :: put_var_real32_1d
    procedure(cable_netcdf_file_put_var_real32_2d), deferred :: put_var_real32_2d
    procedure(cable_netcdf_file_put_var_real32_3d), deferred :: put_var_real32_3d
    procedure(cable_netcdf_file_put_var_real64_0d), deferred :: put_var_real64_0d
    procedure(cable_netcdf_file_put_var_real64_1d), deferred :: put_var_real64_1d
    procedure(cable_netcdf_file_put_var_real64_2d), deferred :: put_var_real64_2d
    procedure(cable_netcdf_file_put_var_real64_3d), deferred :: put_var_real64_3d
    generic :: put_var => &
      put_var_int32_0d, put_var_int32_1d, put_var_int32_2d, put_var_int32_3d, &
      put_var_real32_0d, put_var_real32_1d, put_var_real32_2d, put_var_real32_3d, &
      put_var_real64_0d, put_var_real64_1d, put_var_real64_2d, put_var_real64_3d
    procedure(cable_netcdf_file_write_darray_int32_1d), deferred :: write_darray_int32_1d
    procedure(cable_netcdf_file_write_darray_int32_2d), deferred :: write_darray_int32_2d
    procedure(cable_netcdf_file_write_darray_int32_3d), deferred :: write_darray_int32_3d
    procedure(cable_netcdf_file_write_darray_real32_1d), deferred :: write_darray_real32_1d
    procedure(cable_netcdf_file_write_darray_real32_2d), deferred :: write_darray_real32_2d
    procedure(cable_netcdf_file_write_darray_real32_3d), deferred :: write_darray_real32_3d
    procedure(cable_netcdf_file_write_darray_real64_1d), deferred :: write_darray_real64_1d
    procedure(cable_netcdf_file_write_darray_real64_2d), deferred :: write_darray_real64_2d
    procedure(cable_netcdf_file_write_darray_real64_3d), deferred :: write_darray_real64_3d
    generic :: write_darray => &
      write_darray_int32_1d, write_darray_int32_2d, write_darray_int32_3d, &
      write_darray_real32_1d, write_darray_real32_2d, write_darray_real32_3d, &
      write_darray_real64_1d, write_darray_real64_2d, write_darray_real64_3d
    procedure(cable_netcdf_file_get_var_int32_0d), deferred :: get_var_int32_0d
    procedure(cable_netcdf_file_get_var_int32_1d), deferred :: get_var_int32_1d
    procedure(cable_netcdf_file_get_var_int32_2d), deferred :: get_var_int32_2d
    procedure(cable_netcdf_file_get_var_int32_3d), deferred :: get_var_int32_3d
    procedure(cable_netcdf_file_get_var_real32_0d), deferred :: get_var_real32_0d
    procedure(cable_netcdf_file_get_var_real32_1d), deferred :: get_var_real32_1d
    procedure(cable_netcdf_file_get_var_real32_2d), deferred :: get_var_real32_2d
    procedure(cable_netcdf_file_get_var_real32_3d), deferred :: get_var_real32_3d
    procedure(cable_netcdf_file_get_var_real64_0d), deferred :: get_var_real64_0d
    procedure(cable_netcdf_file_get_var_real64_1d), deferred :: get_var_real64_1d
    procedure(cable_netcdf_file_get_var_real64_2d), deferred :: get_var_real64_2d
    procedure(cable_netcdf_file_get_var_real64_3d), deferred :: get_var_real64_3d
    generic :: get_var => &
      get_var_int32_0d, get_var_int32_1d, get_var_int32_2d, get_var_int32_3d, &
      get_var_real32_0d, get_var_real32_1d, get_var_real32_2d, get_var_real32_3d, &
      get_var_real64_0d, get_var_real64_1d, get_var_real64_2d, get_var_real64_3d
    procedure(cable_netcdf_file_read_darray_int32_1d), deferred :: read_darray_int32_1d
    procedure(cable_netcdf_file_read_darray_int32_2d), deferred :: read_darray_int32_2d
    procedure(cable_netcdf_file_read_darray_int32_3d), deferred :: read_darray_int32_3d
    procedure(cable_netcdf_file_read_darray_real32_1d), deferred :: read_darray_real32_1d
    procedure(cable_netcdf_file_read_darray_real32_2d), deferred :: read_darray_real32_2d
    procedure(cable_netcdf_file_read_darray_real32_3d), deferred :: read_darray_real32_3d
    procedure(cable_netcdf_file_read_darray_real64_1d), deferred :: read_darray_real64_1d
    procedure(cable_netcdf_file_read_darray_real64_2d), deferred :: read_darray_real64_2d
    procedure(cable_netcdf_file_read_darray_real64_3d), deferred :: read_darray_real64_3d
    generic :: read_darray => &
      read_darray_int32_1d, read_darray_int32_2d, read_darray_int32_3d, &
      read_darray_real32_1d, read_darray_real32_2d, read_darray_real32_3d, &
      read_darray_real64_1d, read_darray_real64_2d, read_darray_real64_3d
  end type

  abstract interface
    subroutine cable_netcdf_file_close(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    subroutine cable_netcdf_file_end_def(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    subroutine cable_netcdf_file_redef(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    subroutine cable_netcdf_file_sync(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    subroutine cable_netcdf_file_def_dims(this, dim_names, dim_lens)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: dim_names(:)
      integer, intent(in) :: dim_lens(:)
    end subroutine
    subroutine cable_netcdf_file_def_var(this, var_name, dim_names, type)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in), optional :: dim_names(:)
      integer, intent(in) :: type
    end subroutine
    subroutine cable_netcdf_file_put_att_global_string(this, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name, att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_global_int32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_global_real32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_global_real64(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_var_string(this, var_name, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name, att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_var_int32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_var_real32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
    end subroutine
    subroutine cable_netcdf_file_put_att_var_real64(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_global_string(this, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      character(len=*), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_global_int32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_global_real32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_global_real64(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_var_string(this, var_name, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      character(len=*), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_var_int32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_var_real32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_get_att_var_real64(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name, att_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    end subroutine
    subroutine cable_netcdf_file_inq_dim_len(this, dim_name, dim_len)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: dim_name
      integer, intent(out) :: dim_len
    end subroutine
    subroutine cable_netcdf_file_put_var_int32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_int32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_int32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_int32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real64_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real64_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real64_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_put_var_real64_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_write_darray_int32_1d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_int32_2d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_int32_3d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_real32_1d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_real32_2d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_real32_3d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_real64_1d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_real64_2d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_write_darray_real64_3d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_get_var_int32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_int32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_int32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_int32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real64_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real64_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real64_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_get_var_real64_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
      integer, intent(in), optional :: start(:), count(:)
    end subroutine
    subroutine cable_netcdf_file_read_darray_int32_1d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_int32_2d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_int32_3d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_real32_1d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_real32_2d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_real32_3d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_real64_1d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_real64_2d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
    subroutine cable_netcdf_file_read_darray_real64_3d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
      integer, intent(in), optional :: frame
    end subroutine
  end interface

  type, abstract :: cable_netcdf_io_t
  contains
    procedure(cable_netcdf_io_init), deferred :: init
    procedure(cable_netcdf_io_finalise), deferred :: finalise
    procedure(cable_netcdf_io_create_file), deferred :: create_file
    procedure(cable_netcdf_io_open_file), deferred :: open_file
    procedure(cable_netcdf_io_create_decomp), deferred :: create_decomp
  end type

  abstract interface
    subroutine cable_netcdf_io_init(this)
      import cable_netcdf_io_t
      class(cable_netcdf_io_t), intent(inout) :: this
    end subroutine
    subroutine cable_netcdf_io_finalise(this)
      import cable_netcdf_io_t
      class(cable_netcdf_io_t), intent(inout) :: this
    end subroutine
    function cable_netcdf_io_create_file(this, path) result(file)
      import cable_netcdf_io_t, cable_netcdf_file_t
      class(cable_netcdf_io_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      class(cable_netcdf_file_t), allocatable :: file
    end function
    function cable_netcdf_io_open_file(this, path) result(file)
      import cable_netcdf_io_t, cable_netcdf_file_t
      class(cable_netcdf_io_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      class(cable_netcdf_file_t), allocatable :: file
    end function
    function cable_netcdf_io_create_decomp(this, compmap, dims, type) result(decomp)
      import cable_netcdf_io_t, cable_netcdf_decomp_t
      class(cable_netcdf_io_t), intent(inout) :: this
      integer, intent(in) :: compmap(:), dims(:)
      integer, intent(in) :: type
      class(cable_netcdf_decomp_t), allocatable :: decomp
    end function
  end interface

  interface
    module subroutine cable_netcdf_mod_init(mpi_grp)
      type(mpi_grp_t), intent(in) :: mpi_grp
    end subroutine
    module subroutine cable_netcdf_mod_end()
    end subroutine
    module function cable_netcdf_create_file(path) result(file)
      character(len=*), intent(in) :: path
      class(cable_netcdf_file_t), allocatable :: file
    end function
    module function cable_netcdf_open_file(path) result(file)
      character(len=*), intent(in) :: path
      class(cable_netcdf_file_t), allocatable :: file
    end function
    module function cable_netcdf_create_decomp(compmap, dims, type) result(decomp)
      integer, intent(in) :: compmap(:), dims(:), type
      class(cable_netcdf_decomp_t), allocatable :: decomp
    end function
  end interface

  class(cable_netcdf_io_t), allocatable :: cable_netcdf_io_handler

end module cable_netcdf_mod
