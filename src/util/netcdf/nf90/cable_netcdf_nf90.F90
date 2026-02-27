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
  use netcdf, only: nf90_redef
  use netcdf, only: NF90_NOERR
  use netcdf, only: NF90_NETCDF4
  use netcdf, only: NF90_CLASSIC_MODEL
  use netcdf, only: NF90_CLOBBER
  use netcdf, only: NF90_NOCLOBBER
  use netcdf, only: NF90_WRITE
  use netcdf, only: NF90_NOWRITE
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

  type, extends(cable_netcdf_decomp_t) :: cable_netcdf_nf90_decomp_int32_t
    integer(kind=CABLE_NETCDF_INT32_KIND), allocatable :: values_filled(:)
  end type

  type, extends(cable_netcdf_decomp_t) :: cable_netcdf_nf90_decomp_real32_t
    real(kind=CABLE_NETCDF_REAL32_KIND), allocatable :: values_filled(:)
  end type

  type, extends(cable_netcdf_decomp_t) :: cable_netcdf_nf90_decomp_real64_t
    real(kind=CABLE_NETCDF_REAL64_KIND), allocatable :: values_filled(:)
  end type

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
    procedure :: redef => cable_netcdf_nf90_file_redef
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
    procedure :: inq_var_ndims => cable_netcdf_nf90_file_inq_var_ndims
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

  function cmode_nf90(iotype, mode)
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    integer :: cmode_nf90
    select case(iotype)
    case (CABLE_NETCDF_IOTYPE_CLASSIC)
      cmode_nf90 = NF90_CLASSIC_MODEL
    case (CABLE_NETCDF_IOTYPE_NETCDF4C, CABLE_NETCDF_IOTYPE_NETCDF4P)
      cmode_nf90 = NF90_NETCDF4
    case default
      call cable_abort("Error: iotype not supported", __FILE__, __LINE__)
    end select
    if (present(mode)) then
      select case(mode)
      case (CABLE_NETCDF_MODE_NOCLOBBER)
        cmode_nf90 = ior(cmode_nf90, NF90_NOCLOBBER)
      case (CABLE_NETCDF_MODE_CLOBBER)
        cmode_nf90 = ior(cmode_nf90, NF90_CLOBBER)
      case (CABLE_NETCDF_MODE_WRITE)
        cmode_nf90 = ior(cmode_nf90, NF90_WRITE)
      case (CABLE_NETCDF_MODE_NOWRITE)
        cmode_nf90 = ior(cmode_nf90, NF90_NOWRITE)
      case default
        call cable_abort("Error: mode not supported", __FILE__, __LINE__)
      end select
    end if
  end function cmode_nf90

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

  function cable_netcdf_nf90_io_create_file(this, path, iotype, mode) result(file)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    class(cable_netcdf_file_t), allocatable :: file
    integer :: ncid
    call check_nf90(nf90_create(path, cmode_nf90(iotype, mode), ncid))
    file = cable_netcdf_nf90_file_t(ncid=ncid)
  end function

  function cable_netcdf_nf90_io_open_file(this, path, iotype, mode) result(file)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    class(cable_netcdf_file_t), allocatable :: file
    integer :: ncid
    call check_nf90(nf90_open(path, cmode_nf90(iotype, mode), ncid))
    file = cable_netcdf_nf90_file_t(ncid=ncid)
  end function

  function cable_netcdf_nf90_io_create_decomp(this, compmap, dims, type) result(decomp)
    class(cable_netcdf_nf90_io_t), intent(inout) :: this
    integer, intent(in) :: compmap(:), dims(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp
    select case(type)
    case (CABLE_NETCDF_INT)
      block
        type(cable_netcdf_nf90_decomp_int32_t) :: decomp_int32
        decomp_int32 = cable_netcdf_nf90_decomp_int32_t(compmap=compmap, dims=dims, type=type)
        allocate(decomp_int32%values_filled(product(dims)))
        decomp = decomp_int32
      end block
    case (CABLE_NETCDF_FLOAT)
      block
        type(cable_netcdf_nf90_decomp_real32_t) :: decomp_real32
        decomp_real32 = cable_netcdf_nf90_decomp_real32_t(compmap=compmap, dims=dims, type=type)
        allocate(decomp_real32%values_filled(product(dims)))
        decomp = decomp_real32
      end block
    case (CABLE_NETCDF_DOUBLE)
      block
        type(cable_netcdf_nf90_decomp_real64_t) :: decomp_real64
        decomp_real64 = cable_netcdf_nf90_decomp_real64_t(compmap=compmap, dims=dims, type=type)
        allocate(decomp_real64%values_filled(product(dims)))
        decomp = decomp_real64
      end block
    case default
      call cable_abort("cable_netcdf_nf90_mod: Error: type not supported")
    end select
  end function

  subroutine cable_netcdf_nf90_file_close(this)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    call check_nf90(nf90_close(this%ncid))
  end subroutine

  subroutine cable_netcdf_nf90_file_end_def(this)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    call check_nf90(nf90_enddef(this%ncid))
  end subroutine

  subroutine cable_netcdf_nf90_file_redef(this)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    call check_nf90(nf90_redef(this%ncid))
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
    character(len=*), intent(in) :: var_name
    character(len=*), intent(in), optional :: dim_names(:)
    integer, intent(in) :: type
    integer, allocatable :: dimids(:)
    integer :: i, tmp

    if (.not. present(dim_names)) then
      call check_nf90(nf90_def_var(this%ncid, var_name, type_nf90(type), tmp))
      return
    end if

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

  subroutine cable_netcdf_nf90_file_inq_var_ndims(this, var_name, ndims)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer, intent(out) :: ndims
    integer :: varid
    call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
    call check_nf90(nf90_inquire_variable(this%ncid, varid, ndims=ndims))
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

  subroutine cable_netcdf_nf90_file_write_darray_int32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_int32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_INT
      end if
      do i = 1, size(values)
        decomp%values_filled(decomp%compmap(i)) = values(i)
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_int32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_int32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_int32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_INT
      end if
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        decomp%values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2))
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_int32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_int32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_int32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_INT
      end if
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        decomp%values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2), mem_index(3))
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_int32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real32_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_INT
      end if
      do i = 1, size(values)
        decomp%values_filled(decomp%compmap(i)) = values(i)
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real32_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_INT
      end if
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        decomp%values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2))
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real32_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_INT
      end if
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        decomp%values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2), mem_index(3))
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real64_1d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real64_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_DOUBLE
      end if
      do i = 1, size(values)
        decomp%values_filled(decomp%compmap(i)) = values(i)
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real64_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real64_2d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real64_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_DOUBLE
      end if
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        decomp%values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2))
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real64_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_write_darray_real64_3d(this, var_name, values, decomp, fill_value, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real64_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
      if (present(fill_value)) then
        decomp%values_filled = fill_value
      else
        decomp%values_filled = NF90_FILL_DOUBLE
      end if
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        decomp%values_filled(decomp%compmap(i)) = values(mem_index(1), mem_index(2), mem_index(3))
      end do
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real64_t", file=__FILE__, line=__LINE__)
    end select
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

  subroutine cable_netcdf_nf90_file_read_darray_int32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_int32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        values(i) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_int32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_int32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_int32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values(mem_index(1), mem_index(2)) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_int32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_int32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_int32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values(mem_index(1), mem_index(2), mem_index(3)) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_int32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real32_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        values(i) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real32_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values(mem_index(1), mem_index(2)) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real32_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real32_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values(mem_index(1), mem_index(2), mem_index(3)) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real32_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real64_1d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real64_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        values(i) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real64_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real64_2d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real64_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:2))
        values(mem_index(1), mem_index(2)) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real64_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

  subroutine cable_netcdf_nf90_file_read_darray_real64_3d(this, var_name, values, decomp, frame)
    class(cable_netcdf_nf90_file_t), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
    class(cable_netcdf_decomp_t), intent(inout) :: decomp
    integer, intent(in), optional :: frame
    integer :: i, varid, ndims
    integer :: dimids(NF90_MAX_VAR_DIMS), starts(NF90_MAX_VAR_DIMS), counts(NF90_MAX_VAR_DIMS), mem_index(CABLE_NETCDF_MAX_RANK)
    select type (decomp)
    type is (cable_netcdf_nf90_decomp_real64_t)
      call check_nf90(nf90_inq_varid(this%ncid, var_name, varid))
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
          decomp%values_filled, &
          start=starts(:ndims), &
          count=counts(:ndims), &
          map=[1, (product(decomp%dims(:i)), i = 1, size(decomp%dims) - 1)] &
        ) &
      )
      do i = 1, size(values)
        call array_index(i, shape(values), mem_index(:3))
        values(mem_index(1), mem_index(2), mem_index(3)) = decomp%values_filled(decomp%compmap(i))
      end do
    class default
      call cable_abort("Error: decomp must be of type cable_netcdf_nf90_decomp_real64_t", file=__FILE__, line=__LINE__)
    end select
  end subroutine

end module cable_netcdf_nf90_mod
