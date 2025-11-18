module cable_netcdf_nf90_mod
  use cable_netcdf_mod

  use cable_mpi_mod, only: mpi_grp_t
  use cable_abort_module, only: cable_abort

  use cable_array_utils_mod, only: array_index

  use netcdf, only: nf90_create
  use netcdf, only: nf90_open
  use netcdf, only: nf90_close
  use netcdf, only: nf90_sync
  use netcdf, only: nf90_strerror
  use netcdf, only: nf90_def_dim
  use netcdf, only: nf90_def_var
  use netcdf, only: nf90_put_att
  use netcdf, only: nf90_get_att
  use netcdf, only: nf90_put_var
  use netcdf, only: nf90_get_var
  use netcdf, only: nf90_inq_dimid
  use netcdf, only: nf90_inq_varid
  use netcdf, only: nf90_inquire
  use netcdf, only: nf90_inquire_dimension
  use netcdf, only: nf90_inquire_variable
  use netcdf, only: nf90_enddef
  use netcdf, only: NF90_NOERR
  use netcdf, only: NF90_NETCDF4
  use netcdf, only: NF90_UNLIMITED
  use netcdf, only: NF90_INT
  use netcdf, only: NF90_FLOAT
  use netcdf, only: NF90_DOUBLE
  use netcdf, only: NF90_FILL_INT
  use netcdf, only: NF90_FILL_FLOAT
  use netcdf, only: NF90_FILL_DOUBLE
  use netcdf, only: NF90_MAX_VAR_DIMS
  use netcdf, only: NF90_GLOBAL

  implicit none

  private

  public :: cable_netcdf_nf90_io_t

  type, extends(cable_netcdf_io_t) :: cable_netcdf_nf90_io_t
  contains
    procedure :: init => cable_netcdf_nf90_io_init
    procedure :: finalise => cable_netcdf_nf90_io_finalise
    procedure :: create_file => cable_netcdf_nf90_io_create_file
    procedure :: open_file => cable_netcdf_nf90_io_open_file
    procedure :: create_decomp => cable_netcdf_nf90_io_create_decomp
  end type

  type, extends(cable_netcdf_file_t) :: cable_netcdf_nf90_file_t
    integer, private :: ncid
  contains
    procedure :: close => cable_netcdf_nf90_file_close
    procedure :: end_def => cable_netcdf_nf90_file_end_def
    procedure :: sync => cable_netcdf_nf90_file_sync
    procedure :: def_dims => cable_netcdf_nf90_file_def_dims
    procedure :: def_var => cable_netcdf_nf90_file_def_var
    procedure :: put_att_global_string => cable_netcdf_nf90_file_put_att_global_string
    procedure :: put_att_global_int32 => cable_netcdf_nf90_file_put_att_global_int32
    procedure :: put_att_global_real32 => cable_netcdf_nf90_file_put_att_global_real32
    procedure :: put_att_global_real64 => cable_netcdf_nf90_file_put_att_global_real64
    procedure :: put_att_var_string => cable_netcdf_nf90_file_put_att_var_string
    procedure :: put_att_var_int32 => cable_netcdf_nf90_file_put_att_var_int32
    procedure :: put_att_var_real32 => cable_netcdf_nf90_file_put_att_var_real32
    procedure :: put_att_var_real64 => cable_netcdf_nf90_file_put_att_var_real64
    procedure :: get_att_global_string => cable_netcdf_nf90_file_get_att_global_string
    procedure :: get_att_global_int32 => cable_netcdf_nf90_file_get_att_global_int32
    procedure :: get_att_global_real32 => cable_netcdf_nf90_file_get_att_global_real32
    procedure :: get_att_global_real64 => cable_netcdf_nf90_file_get_att_global_real64
    procedure :: get_att_var_string => cable_netcdf_nf90_file_get_att_var_string
    procedure :: get_att_var_int32 => cable_netcdf_nf90_file_get_att_var_int32
    procedure :: get_att_var_real32 => cable_netcdf_nf90_file_get_att_var_real32
    procedure :: get_att_var_real64 => cable_netcdf_nf90_file_get_att_var_real64
    procedure :: inq_dim_len => cable_netcdf_nf90_file_inq_dim_len
    procedure :: put_var_int32_0d => cable_netcdf_nf90_file_put_var_int32_0d
    procedure :: put_var_int32_1d => cable_netcdf_nf90_file_put_var_int32_1d
    procedure :: put_var_int32_2d => cable_netcdf_nf90_file_put_var_int32_2d
    procedure :: put_var_int32_3d => cable_netcdf_nf90_file_put_var_int32_3d
    procedure :: put_var_real32_0d => cable_netcdf_nf90_file_put_var_real32_0d
    procedure :: put_var_real32_1d => cable_netcdf_nf90_file_put_var_real32_1d
    procedure :: put_var_real32_2d => cable_netcdf_nf90_file_put_var_real32_2d
    procedure :: put_var_real32_3d => cable_netcdf_nf90_file_put_var_real32_3d
    procedure :: put_var_real64_0d => cable_netcdf_nf90_file_put_var_real64_0d
    procedure :: put_var_real64_1d => cable_netcdf_nf90_file_put_var_real64_1d
    procedure :: put_var_real64_2d => cable_netcdf_nf90_file_put_var_real64_2d
    procedure :: put_var_real64_3d => cable_netcdf_nf90_file_put_var_real64_3d
    procedure :: write_darray_int32_1d => cable_netcdf_nf90_file_write_darray_int32_1d
    procedure :: write_darray_int32_2d => cable_netcdf_nf90_file_write_darray_int32_2d
    procedure :: write_darray_int32_3d => cable_netcdf_nf90_file_write_darray_int32_3d
    procedure :: write_darray_real32_1d => cable_netcdf_nf90_file_write_darray_real32_1d
    procedure :: write_darray_real32_2d => cable_netcdf_nf90_file_write_darray_real32_2d
    procedure :: write_darray_real32_3d => cable_netcdf_nf90_file_write_darray_real32_3d
    procedure :: write_darray_real64_1d => cable_netcdf_nf90_file_write_darray_real64_1d
    procedure :: write_darray_real64_2d => cable_netcdf_nf90_file_write_darray_real64_2d
    procedure :: write_darray_real64_3d => cable_netcdf_nf90_file_write_darray_real64_3d
    procedure :: get_var_int32_0d => cable_netcdf_nf90_file_get_var_int32_0d
    procedure :: get_var_int32_1d => cable_netcdf_nf90_file_get_var_int32_1d
    procedure :: get_var_int32_2d => cable_netcdf_nf90_file_get_var_int32_2d
    procedure :: get_var_int32_3d => cable_netcdf_nf90_file_get_var_int32_3d
    procedure :: get_var_real32_0d => cable_netcdf_nf90_file_get_var_real32_0d
    procedure :: get_var_real32_1d => cable_netcdf_nf90_file_get_var_real32_1d
    procedure :: get_var_real32_2d => cable_netcdf_nf90_file_get_var_real32_2d
    procedure :: get_var_real32_3d => cable_netcdf_nf90_file_get_var_real32_3d
    procedure :: get_var_real64_0d => cable_netcdf_nf90_file_get_var_real64_0d
    procedure :: get_var_real64_1d => cable_netcdf_nf90_file_get_var_real64_1d
    procedure :: get_var_real64_2d => cable_netcdf_nf90_file_get_var_real64_2d
    procedure :: get_var_real64_3d => cable_netcdf_nf90_file_get_var_real64_3d
    procedure :: read_darray_int32_1d => cable_netcdf_nf90_file_read_darray_int32_1d
    procedure :: read_darray_int32_2d => cable_netcdf_nf90_file_read_darray_int32_2d
    procedure :: read_darray_int32_3d => cable_netcdf_nf90_file_read_darray_int32_3d
    procedure :: read_darray_real32_1d => cable_netcdf_nf90_file_read_darray_real32_1d
    procedure :: read_darray_real32_2d => cable_netcdf_nf90_file_read_darray_real32_2d
    procedure :: read_darray_real32_3d => cable_netcdf_nf90_file_read_darray_real32_3d
    procedure :: read_darray_real64_1d => cable_netcdf_nf90_file_read_darray_real64_1d
    procedure :: read_darray_real64_2d => cable_netcdf_nf90_file_read_darray_real64_2d
    procedure :: read_darray_real64_3d => cable_netcdf_nf90_file_read_darray_real64_3d
  end type

contains

  function type_nf90(type)
    integer, intent(in) :: type
    integer :: type_nf90
    select case(type)
    case(CABLE_NETCDF_INT)
      type_nf90 = NF90_INT
    case(CABLE_NETCDF_FLOAT)
      type_nf90 = NF90_FLOAT
    case(CABLE_NETCDF_DOUBLE)
      type_nf90 = NF90_DOUBLE
    case default
      call cable_abort("cable_netcdf_nf90_mod: Error: type not supported")
    end select
  end function type_nf90

  subroutine check_nf90(status)
    integer, intent ( in) :: status
    if(status /= NF90_NOERR) then
      call cable_abort(trim(nf90_strerror(status)), file=__FILE__, line=__LINE__)
    end if
  end subroutine check_nf90

  subroutine cable_netcdf_nf90_io_init(this)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
  end subroutine

  subroutine cable_netcdf_nf90_io_finalise(this)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
  end subroutine

  function cable_netcdf_nf90_io_create_file(this, path) result(file)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    class(cable_netcdf_file_t), allocatable :: file
    integer :: ncid
    call check_nf90(nf90_create(path, NF90_NETCDF4, ncid))
    file = cable_netcdf_nf90_file_t(ncid=ncid)
  end function

  function cable_netcdf_nf90_io_open_file(this, path) result(file)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    class(cable_netcdf_file_t), allocatable :: file
    integer :: ncid
    call check_nf90(nf90_open(path, NF90_NETCDF4, ncid))
    file = cable_netcdf_nf90_file_t(ncid=ncid)
  end function

  function cable_netcdf_nf90_io_create_decomp(this, compmap, dims, type) result(decomp)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
    integer, intent(in) :: compmap(:), dims(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp
    decomp = cable_netcdf_decomp_t(compmap=compmap, dims=dims, type=type)
  end function

  subroutine cable_netcdf_nf90_file_close(this)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    call check_nf90(nf90_close(this%ncid))
  end subroutine

  subroutine cable_netcdf_nf90_file_end_def(this)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    call check_nf90(nf90_enddef(this%ncid))
  end subroutine

  subroutine cable_netcdf_nf90_file_sync(this)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    call check_nf90(nf90_sync(this%ncid))
  end subroutine

  subroutine cable_netcdf_nf90_file_def_dims(this, dim_names, dim_lens)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: dim_names(:)
    integer, intent(in) :: dim_lens(:)
    integer :: i, tmp
    do i = 1, size(dim_names)
      if (dim_lens(i) == CABLE_NETCDF_UNLIMITED) then
        call check_nf90(nf90_def_dim(this%ncid, dim_names(i), NF90_UNLIMITED, tmp))
      else
        call check_nf90(nf90_def_dim(this%ncid, dim_names(i), dim_lens(i), tmp))
      end if
    end do
  end subroutine

  subroutine cable_netcdf_nf90_file_def_var(this, var_name, dim_names, type)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, dim_names(:)
    integer, intent(in) :: type
    integer, allocatable :: dimids(:)
    integer :: i, tmp
    allocate(dimids(size(dim_names)))
    do i = 1, size(dimids)
      call check_nf90(nf90_inq_dimid(this%ncid, dim_names(i), dimids(i)))
    end do
    call check_nf90(nf90_def_var(this%ncid, var_name, type_nf90(type), dimids, tmp))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_global_string(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name, att_value
    call check_nf90(nf90_put_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_global_int32(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
    call check_nf90(nf90_put_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_global_real32(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
    call check_nf90(nf90_put_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_global_real64(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
    call check_nf90(nf90_put_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_var_string(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name, att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_var_int32(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_var_real32(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_att_var_real64(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_global_string(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    character(len=*), intent(out) :: att_value
    call check_nf90(nf90_get_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_global_int32(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    call check_nf90(nf90_get_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_global_real32(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    call check_nf90(nf90_get_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_global_real64(this, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    call check_nf90(nf90_get_att(this%ncid, NF90_GLOBAL, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_var_string(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    character(len=*), intent(out) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_var_int32(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_var_real32(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_att_var_real64(this, var_name, att_name, att_value)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name, att_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_att(this%ncid, varid, att_name, att_value))
  end subroutine

  subroutine cable_netcdf_nf90_file_inq_dim_len(this, dim_name, dim_len)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: dim_name
    integer, intent(out) :: dim_len
    integer :: dimid
    call check_nf90(nf90_inq_dimid(this%ncid, dim_name, dimid))
    call check_nf90(nf90_inquire_dimension(this%ncid, dimid, len=dim_len))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_int32_0d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_int32_1d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_int32_2d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_int32_3d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real32_0d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real32_1d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real32_2d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real32_3d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real64_0d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real64_1d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real64_2d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_put_var_real64_3d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_put_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: values_filled(:)
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    if (present(fill_value)) then
      allocate(values_filled(product(decomp%dims)), source=fill_value)
    else
      allocate(values_filled(product(decomp%dims)), source=NF90_FILL_INT)
    end if
    select rank(values)
    rank(1)
      do i = 1, size(values)
        values_filled(decomp%compmap(i)) = values(i)
      end do
    rank(2)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2))
      end do
    rank(3)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2), mem_index(3))
      end do
    end select
    call check_nf90(nf90_inquire_variable(this%ncid, varid, dimids=dimids, ndims=ndims))
    do i = 1, ndims
      starts(i) = 1
      call check_nf90(nf90_inquire_dimension(this%ncid, dimids(i), len=counts(i)))
    end do
    if (present(frame)) then
      starts(ndims) = frame
      counts(ndims) = 1
    end if
    call check_nf90( &
      nf90_put_var( &
        this%ncid, &
        varid, &
        values_filled, &
        start=starts(:ndims), &
        count=counts(:ndims), &
        map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
      ) &
    )
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_int32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_int32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_int32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_int32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: values_filled(:)
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    if (present(fill_value)) then
      allocate(values_filled(product(decomp%dims)), source=fill_value)
    else
      allocate(values_filled(product(decomp%dims)), source=NF90_FILL_FLOAT)
    end if
    select rank(values)
    rank(1)
      do i = 1, size(values)
        values_filled(decomp%compmap(i)) = values(i)
      end do
    rank(2)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2))
      end do
    rank(3)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2), mem_index(3))
      end do
    end select
    call check_nf90(nf90_inquire_variable(this%ncid, varid, dimids=dimids, ndims=ndims))
    do i = 1, ndims
      starts(i) = 1
      call check_nf90(nf90_inquire_dimension(this%ncid, dimids(i), len=counts(i)))
    end do
    if (present(frame)) then
      starts(ndims) = frame
      counts(ndims) = 1
    end if
    call check_nf90( &
      nf90_put_var( &
        this%ncid, &
        varid, &
        values_filled, &
        start=starts(:ndims), &
        count=counts(:ndims), &
        map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
      ) &
    )
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_real32(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: values_filled(:)
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    if (present(fill_value)) then
      allocate(values_filled(product(decomp%dims)), source=fill_value)
    else
      allocate(values_filled(product(decomp%dims)), source=NF90_FILL_DOUBLE)
    end if
    select rank(values)
    rank(1)
      do i = 1, size(values)
        values_filled(decomp%compmap(i)) = values(i)
      end do
    rank(2)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2))
      end do
    rank(3)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2), mem_index(3))
      end do
    end select
    call check_nf90(nf90_inquire_variable(this%ncid, varid, dimids=dimids, ndims=ndims))
    do i = 1, ndims
      starts(i) = 1
      call check_nf90(nf90_inquire_dimension(this%ncid, dimids(i), len=counts(i)))
    end do
    if (present(frame)) then
      starts(ndims) = frame
      counts(ndims) = 1
    end if
    call check_nf90( &
      nf90_put_var( &
        this%ncid, &
        varid, &
        values_filled, &
        start=starts(:ndims), &
        count=counts(:ndims), &
        map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
      ) &
    )
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real64_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real64_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real64_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    call cable_netcdf_write_darray_real64(this, var_name, values, decomp, fill_value, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_int32_0d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_int32_1d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_int32_2d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_int32_3d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real32_0d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real32_1d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real32_2d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real32_3d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real64_0d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real64_1d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real64_2d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_nf90_file_get_var_real64_3d(this, var_name, values, start, count)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    integer, intent(in), optional :: start(:), count(:)
    integer varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_get_var(this%ncid, varid, values, start=start, count=count))
  end subroutine

  subroutine cable_netcdf_read_darray_int32(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: values_filled(:)
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    allocate(values_filled(product(decomp%dims)))
    call check_nf90(nf90_inquire_variable(this%ncid, varid, dimids=dimids, ndims=ndims))
    do i = 1, ndims
      starts(i) = 1
      call check_nf90(nf90_inquire_dimension(this%ncid, dimids(i), len=counts(i)))
    end do
    if (present(frame)) then
      starts(ndims) = frame
      counts(ndims) = 1
    end if
    call check_nf90( &
      nf90_get_var( &
        this%ncid, &
        varid, &
        values_filled, &
        start=starts(:ndims), &
        count=counts(:ndims), &
        map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
      ) &
    )
    select rank(values)
    rank(1)
      do i = 1, size(values)
        values(i) = values_filled(decomp%compmap(i))
      end do
    rank(2)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values(mem_index(1), mem_index(2)) = values_filled(decomp%compmap(i))
      end do
    rank(3)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values(mem_index(1), mem_index(2), mem_index(3)) = values_filled(decomp%compmap(i))
      end do
    end select
  end subroutine cable_netcdf_read_darray_int32

  subroutine cable_netcdf_nf90_file_read_darray_int32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_int32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_int32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_int32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_int32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_int32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_read_darray_real32(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: values_filled(:)
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    allocate(values_filled(product(decomp%dims)))
    call check_nf90(nf90_inquire_variable(this%ncid, varid, dimids=dimids, ndims=ndims))
    do i = 1, ndims
      starts(i) = 1
      call check_nf90(nf90_inquire_dimension(this%ncid, dimids(i), len=counts(i)))
    end do
    if (present(frame)) then
      starts(ndims) = frame
      counts(ndims) = 1
    end if
    call check_nf90( &
      nf90_get_var( &
        this%ncid, &
        varid, &
        values_filled, &
        start=starts(:ndims), &
        count=counts(:ndims), &
        map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
      ) &
    )
    select rank(values)
    rank(1)
      do i = 1, size(values)
        values(i) = values_filled(decomp%compmap(i))
      end do
    rank(2)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values(mem_index(1), mem_index(2)) = values_filled(decomp%compmap(i))
      end do
    rank(3)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values(mem_index(1), mem_index(2), mem_index(3)) = values_filled(decomp%compmap(i))
      end do
    end select
  end subroutine cable_netcdf_read_darray_real32

  subroutine cable_netcdf_nf90_file_read_darray_real32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_real32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_real32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_real32(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_read_darray_real64(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(..)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: values_filled(:)
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    allocate(values_filled(product(decomp%dims)))
    call check_nf90(nf90_inquire_variable(this%ncid, varid, dimids=dimids, ndims=ndims))
    do i = 1, ndims
      starts(i) = 1
      call check_nf90(nf90_inquire_dimension(this%ncid, dimids(i), len=counts(i)))
    end do
    if (present(frame)) then
      starts(ndims) = frame
      counts(ndims) = 1
    end if
    call check_nf90( &
      nf90_get_var( &
        this%ncid, &
        varid, &
        values_filled, &
        start=starts(:ndims), &
        count=counts(:ndims), &
        map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
      ) &
    )
    select rank(values)
    rank(1)
      do i = 1, size(values)
        values(i) = values_filled(decomp%compmap(i))
      end do
    rank(2)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values(mem_index(1), mem_index(2)) = values_filled(decomp%compmap(i))
      end do
    rank(3)
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values(mem_index(1), mem_index(2), mem_index(3)) = values_filled(decomp%compmap(i))
      end do
    end select
  end subroutine cable_netcdf_read_darray_real64

  subroutine cable_netcdf_nf90_file_read_darray_real64_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_real64(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real64_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_real64(this, var_name, values, decomp, frame)
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real64_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    call cable_netcdf_read_darray_real64(this, var_name, values, decomp, frame)
  end subroutine

end module cable_netcdf_nf90_mod
