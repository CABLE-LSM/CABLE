! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_netcdf_mod
  !* Interface for netCDF file handling in CABLE.
  !
  ! This module defines abstract types and interfaces for working with netCDF
  ! files using different underlying libraries (e.g., NetCDF Fortran, ParallelIO).
  ! Concrete implementations should extend the abstract types and implement the
  ! deferred procedures.

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
    CABLE_NETCDF_UNLIMITED, &
    CABLE_NETCDF_IOTYPE_NETCDF, &
    CABLE_NETCDF_IOTYPE_CLASSIC, &
    CABLE_NETCDF_IOTYPE_NETCDF4C, &
    CABLE_NETCDF_IOTYPE_NETCDF4P, &
    CABLE_NETCDF_MODE_CLOBBER, &
    CABLE_NETCDF_MODE_NOCLOBBER, &
    CABLE_NETCDF_MODE_WRITE, &
    CABLE_NETCDF_MODE_NOWRITE

  !> Data type constants for netCDF variables and attributes.
  enum, bind(c)
    enumerator :: CABLE_NETCDF_INT
      !! Data type constant for 32-bit integer variables and attributes in netCDF files.
    enumerator :: CABLE_NETCDF_FLOAT
      !! Data type constant for 32-bit real variables and attributes in netCDF files.
    enumerator :: CABLE_NETCDF_DOUBLE
      !! Data type constant for 64-bit real variables and attributes in netCDF files.
  end enum

  !> I/O type options for opening and creating files.
  !! Follows the I/O type options defined in ParallelIO:
  !! <https://github.com/NCAR/ParallelIO/blob/pio2_6_8/src/flib/pio_types.F90#L102-L106>
  enum, bind(c)
    enumerator :: CABLE_NETCDF_IOTYPE_NETCDF
      !! Serial read/write of NetCDF files.
    enumerator :: CABLE_NETCDF_IOTYPE_CLASSIC
      !! Parallel read/write of pNetCDF files.
    enumerator :: CABLE_NETCDF_IOTYPE_NETCDF4C
      !! NetCDF4 (HDF5 format) file opened for compression (serial write access only).
    enumerator :: CABLE_NETCDF_IOTYPE_NETCDF4P
      !! NetCDF4 (HDF5 format) file opened in parallel.
  end enum

  !> File mode options for opening and creating files
  enum, bind(c)
    enumerator :: CABLE_NETCDF_MODE_CLOBBER
      !! Overwrite existing file.
    enumerator :: CABLE_NETCDF_MODE_NOCLOBBER
      !! Do not overwrite existing file.
    enumerator :: CABLE_NETCDF_MODE_WRITE
      !! Open file for writing.
    enumerator :: CABLE_NETCDF_MODE_NOWRITE
      !! Open file for reading only.
  end enum

  integer, parameter :: CABLE_NETCDF_UNLIMITED = -1
    !! Constant to represent an unlimited dimension length in netCDF files.

  type, abstract :: cable_netcdf_decomp_t
    !* Abstract type for describing netCDF decomposition information.
    ! This type describes the mapping from the local in-memory layout of an array
    ! on the current process to the global layout of a netCDF variable on disk
    ! following the degree of freedom decomposition described in Denis et al.
    ! (2011)
    ! [[10.1177/1094342011428143](https://www.doi.org/10.1177/1094342011428143)].
  end type

  type, abstract :: cable_netcdf_file_t
    !* Abstract type for netCDF file handling.
    ! This type defines the interface for operations on netCDF files, such as
    ! defining dimensions and variables, writing and reading data, and managing
    ! attributes.
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
    procedure(cable_netcdf_file_inq_var_ndims), deferred :: inq_var_ndims
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
    !> Close the netCDF file and release any associated resources.
    subroutine cable_netcdf_file_close(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    !> End the definition mode of the netCDF file, allowing data to be written.
    subroutine cable_netcdf_file_end_def(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    !> Re-enter definition mode for the netCDF file, allowing dimensions and
    !! variables to be defined.
    subroutine cable_netcdf_file_redef(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    !> Synchronize the netCDF file, ensuring that all data is written to disk
    !! and visible to other processes.
    subroutine cable_netcdf_file_sync(this)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
    end subroutine
    !> Define dimensions in the netCDF file.
    subroutine cable_netcdf_file_def_dims(this, dim_names, dim_lens)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: dim_names(:)
        !! Array of dimension names to define.
      integer, intent(in) :: dim_lens(:)
        !* Array of dimension lengths corresponding to dim_names. Use
        ! `CABLE_NETCDF_UNLIMITED` for unlimited dimensions.
    end subroutine
    !> Define a variable in the netCDF file with the specified name, dimensions, and type.
    subroutine cable_netcdf_file_def_var(this, var_name, dim_names, type)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to define.
      character(len=*), intent(in), optional :: dim_names(:)
        !* Array of dimension names for the variable. If not provided, the
        ! variable will be defined as a scalar.
      integer, intent(in) :: type
        !* Data type of the variable, using the `CABLE_NETCDF_*` constants (e.g.,
        ! `CABLE_NETCDF_INT`, `CABLE_NETCDF_REAL32`).
    end subroutine
    !> Define a global attribute with a string value in the netCDF file.
    subroutine cable_netcdf_file_put_att_global_string(this, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to define.
      character(len=*), intent(in) :: att_value
        !! Value of the global attribute to define.
    end subroutine
    !> Define a global attribute with an 32-bit integer value in the netCDF file.
    subroutine cable_netcdf_file_put_att_global_int32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to define.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
        !! Value of the global attribute to define.
    end subroutine
    !> Define a global attribute with a 32-bit real value in the netCDF file.
    subroutine cable_netcdf_file_put_att_global_real32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to define.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
        !! Value of the global attribute to define.
    end subroutine
    !> Define a global attribute with a 64-bit real value in the netCDF file.
    subroutine cable_netcdf_file_put_att_global_real64(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to define.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
        !! Value of the global attribute to define.
    end subroutine
    !> Define a variable attribute with a string value in the netCDF file.
    subroutine cable_netcdf_file_put_att_var_string(this, var_name, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to define the attribute for.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to define.
      character(len=*), intent(in) :: att_value
        !! Value of the variable attribute to define.
    end subroutine
    !> Define a variable attribute with an 32-bit integer value in the netCDF file.
    subroutine cable_netcdf_file_put_att_var_int32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to define the attribute for.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to define.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: att_value
        !! Value of the variable attribute to define.
    end subroutine
    !> Define a variable attribute with a 32-bit real value in the netCDF file.
    subroutine cable_netcdf_file_put_att_var_real32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to define the attribute for.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to define.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: att_value
        !! Value of the variable attribute to define.
    end subroutine
    !> Define a variable attribute with a 64-bit real value in the netCDF file.
    subroutine cable_netcdf_file_put_att_var_real64(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to define the attribute for.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to define.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: att_value
        !! Value of the variable attribute to define.
    end subroutine
    !> Read a global attribute with a string value from the netCDF file.
    subroutine cable_netcdf_file_get_att_global_string(this, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to read.
      character(len=*), intent(out) :: att_value
        !! Value of the global attribute read.
    end subroutine
    !> Read a global attribute with an 32-bit integer value from the netCDF file.
    subroutine cable_netcdf_file_get_att_global_int32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to read.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
        !! Value of the global attribute read.
    end subroutine
    !> Read a global attribute with a 32-bit real value from the netCDF file.
    subroutine cable_netcdf_file_get_att_global_real32(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to read.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
        !! Value of the global attribute read.
    end subroutine
    !> Read a global attribute with a 64-bit real value from the netCDF file.
    subroutine cable_netcdf_file_get_att_global_real64(this, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: att_name
        !! Name of the global attribute to read.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
        !! Value of the global attribute read.
    end subroutine
    !> Read a variable attribute with a string value from the netCDF file.
    subroutine cable_netcdf_file_get_att_var_string(this, var_name, att_name, att_value)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read the attribute from.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to read.
      character(len=*), intent(out) :: att_value
        !! Value of the variable attribute read.
    end subroutine
    !> Read a variable attribute with an 32-bit integer value from the netCDF file.
    subroutine cable_netcdf_file_get_att_var_int32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read the attribute from.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to read.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: att_value
        !! Value of the variable attribute read.
    end subroutine
    !> Read a variable attribute with a 32-bit real value from the netCDF file.
    subroutine cable_netcdf_file_get_att_var_real32(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read the attribute from.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to read.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: att_value
        !! Value of the variable attribute read.
    end subroutine
    !> Read a variable attribute with a 64-bit real value from the netCDF file.
    subroutine cable_netcdf_file_get_att_var_real64(this, var_name, att_name, att_value)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read the attribute from.
      character(len=*), intent(in) :: att_name
        !! Name of the variable attribute to read.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: att_value
        !! Value of the variable attribute read.
    end subroutine
    !> Inquire the length of a dimension in the netCDF file by its name.
    subroutine cable_netcdf_file_inq_dim_len(this, dim_name, dim_len)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: dim_name
        !! Name of the dimension to inquire.
      integer, intent(out) :: dim_len
        !* Length of the dimension returned. For unlimited dimensions, the
        ! number of records will be returned.
    end subroutine
    !> Inquire the number of dimensions of a variable in the netCDF file by its name.
    subroutine cable_netcdf_file_inq_var_ndims(this, var_name, ndims)
      import cable_netcdf_file_t
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to inquire.
      integer, intent(out) :: ndims
        !! Number of dimensions of the variable returned.
    end subroutine
    !> Write a 0-dimensional (scalar) 32-bit integer variable to the netCDF file.
    subroutine cable_netcdf_file_put_var_int32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values
        !! Value to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 1-dimensional array of 32-bit integers to the netCDF file.
    subroutine cable_netcdf_file_put_var_int32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
        !! Array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 2-dimensional array of 32-bit integers to the netCDF file.
    subroutine cable_netcdf_file_put_var_int32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
        !! 2-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 3-dimensional array of 32-bit integers to the netCDF file.
    subroutine cable_netcdf_file_put_var_int32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
        !! 3-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 0-dimensional (scalar) 32-bit real variable to the netCDF file.
    subroutine cable_netcdf_file_put_var_real32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values
        !! Scalar value to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 1-dimensional array of 32-bit real values to the netCDF file.
    subroutine cable_netcdf_file_put_var_real32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
        !! 1-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 2-dimensional array of 32-bit real values to the netCDF file.
    subroutine cable_netcdf_file_put_var_real32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
        !! 2-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 3-dimensional array of 32-bit real values to the netCDF file.
    subroutine cable_netcdf_file_put_var_real32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
        !! 3-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 0-dimensional (scalar) 64-bit real variable to the netCDF file.
    subroutine cable_netcdf_file_put_var_real64_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values
        !! The scalar value to write.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 1-dimensional array of 64-bit real values to the netCDF file.
    subroutine cable_netcdf_file_put_var_real64_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
        !! 1-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 2-dimensional array of 64-bit real values to the netCDF file.
    subroutine cable_netcdf_file_put_var_real64_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
        !! 2-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 3-dimensional array of 64-bit real values to the netCDF file.
    subroutine cable_netcdf_file_put_var_real64_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
        !! 3-dimensional array of values to write to the variable.
      integer, intent(in), optional :: start(:)
        !* Starting indices for writing the variable. If not provided, writing
        ! will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !* Count of elements to write for each dimension. If not provided, the
        ! entire variable will be written.
    end subroutine
    !> Write a 1-dimensional array of 32-bit integers to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_int32_1d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:)
        !! 1-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 2-dimensional array of 32-bit integers to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_int32_2d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :)
        !! 2-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 3-dimensional array of 32-bit integers to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_int32_3d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in) :: values(:, :, :)
        !! 3-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 1-dimensional array of 32-bit real values to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_real32_1d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:)
        !! 1-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 2-dimensional array of 32-bit real values to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_real32_2d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :)
        !! 2-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 3-dimensional array of 32-bit real values to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_real32_3d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in) :: values(:, :, :)
        !! 3-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 2-dimensional array of 32-bit real values to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_real64_1d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:)
        !! 1-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 2-dimensional array of 32-bit real values to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_real64_2d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :)
        !! 2-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Write a 3-dimensional array of 32-bit real values to the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_write_darray_real64_3d(this, var_name, values, decomp, fill_value, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to write to.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in) :: values(:, :, :)
        !! 3-dimensional array of values to write to the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(in), optional :: fill_value
        !! Fill value to use for any missing data. If not provided, the default fill value for the variable's type will be used.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be written to the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 0-dimensional (scalar) 32-bit integer variable from the netCDF file.
    subroutine cable_netcdf_file_get_var_int32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values
        !! The scalar value read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 1-dimensional array of 32-bit integers from the netCDF file.
    subroutine cable_netcdf_file_get_var_int32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
        !! 1-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will
    end subroutine
    !> Read a 2-dimensional array of 32-bit integers from the netCDF file.
    subroutine cable_netcdf_file_get_var_int32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
        !! 2-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 3-dimensional array of 32-bit integers from the netCDF file.
    subroutine cable_netcdf_file_get_var_int32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
        !! 3-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 0-dimensional (scalar) 32-bit real variable from the netCDF file.
    subroutine cable_netcdf_file_get_var_real32_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values
        !! Scalar value to store the value read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 1-dimensional array of 32-bit real values from the netCDF file.
    subroutine cable_netcdf_file_get_var_real32_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
        !! 1-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 2-dimensional array of 32-bit real values from the netCDF file.
    subroutine cable_netcdf_file_get_var_real32_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
        !! 2-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 3-dimensional array of 32-bit real values from the netCDF file.
    subroutine cable_netcdf_file_get_var_real32_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
        !! 3-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 0-dimensional (scalar) 64-bit real variable from the netCDF file.
    subroutine cable_netcdf_file_get_var_real64_0d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values
        !! Scalar value to store the value read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 1-dimensional array of 64-bit real values from the netCDF file.
    subroutine cable_netcdf_file_get_var_real64_1d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
        !! 1-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 2-dimensional array of 64-bit real values from the netCDF file.
    subroutine cable_netcdf_file_get_var_real64_2d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
        !! 2-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 3-dimensional array of 64-bit real values from the netCDF file.
    subroutine cable_netcdf_file_get_var_real64_3d(this, var_name, values, start, count)
      import cable_netcdf_file_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
        !! 3-dimensional array to store the values read from the variable.
      integer, intent(in), optional :: start(:)
        !! Starting indices for reading the variable. If not provided, reading will start at the beginning of the variable.
      integer, intent(in), optional :: count(:)
        !! Count of elements to read for each dimension. If not provided, the entire variable will be read.
    end subroutine
    !> Read a 1-dimensional array of 32-bit integers from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_int32_1d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:)
        !! 1-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 2-dimensional array of 32-bit integers from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_int32_2d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :)
        !! 2-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 3-dimensional array of 32-bit integers from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_int32_3d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_INT32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      integer(kind=CABLE_NETCDF_INT32_KIND), intent(out) :: values(:, :, :)
        !! 3-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 1-dimensional array of 32-bit real values from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_real32_1d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:)
        !! 1-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 2-dimensional array of 32-bit real values from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_real32_2d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :)
        !! 2-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 3-dimensional array of 32-bit real values from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_real32_3d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL32_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL32_KIND), intent(out) :: values(:, :, :)
        !! 3-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 1-dimensional array of 64-bit real values from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_real64_1d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:)
        !! 1-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 2-dimensional array of 64-bit real values from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_real64_2d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :)
        !! 2-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
    !> Read a 3-dimensional array of 64-bit real values from the netCDF file using decomposition information for parallel I/O.
    subroutine cable_netcdf_file_read_darray_real64_3d(this, var_name, values, decomp, frame)
      import cable_netcdf_file_t, cable_netcdf_decomp_t, CABLE_NETCDF_REAL64_KIND
      class(cable_netcdf_file_t), intent(inout) :: this
      character(len=*), intent(in) :: var_name
        !! Name of the variable to read from.
      real(kind=CABLE_NETCDF_REAL64_KIND), intent(out) :: values(:, :, :)
        !! 3-dimensional array to store the values read from the variable.
      class(cable_netcdf_decomp_t), intent(inout) :: decomp
        !! Decomposition information for parallel I/O.
      integer, intent(in), optional :: frame
        !! Optional frame index for record variables. If provided, the data will be read from the specified frame along the unlimited dimension.
    end subroutine
  end interface

  type, abstract :: cable_netcdf_io_t
    !* Abstract type defining the interface for netCDF I/O handlers.
    ! This allows for different implementations (e.g. NetCDF, ParallelIO) to be
    ! used interchangeably within the CABLE code
  contains
    procedure(cable_netcdf_io_init), deferred :: init
    procedure(cable_netcdf_io_finalise), deferred :: finalise
    procedure(cable_netcdf_io_create_file), deferred :: create_file
    procedure(cable_netcdf_io_open_file), deferred :: open_file
    procedure(cable_netcdf_io_create_decomp), deferred :: create_decomp
  end type

  abstract interface
    !> Initialise the netcdf I/O handler.
    !! This procedure is called during module initialization and can be used to set
    !! up any necessary state or resources for the I/O handler.
    subroutine cable_netcdf_io_init(this)
      import cable_netcdf_io_t
      class(cable_netcdf_io_t), intent(inout) :: this
    end subroutine
    !> Finalise the netcdf I/O handler.
    !! This procedure is called during module finalization and can be used to clean
    !! up any resources allocated by the I/O handler.
    subroutine cable_netcdf_io_finalise(this)
      import cable_netcdf_io_t
      class(cable_netcdf_io_t), intent(inout) :: this
    end subroutine
    !> Create a new netCDF file with the specified path and I/O type.
    function cable_netcdf_io_create_file(this, path, iotype, mode) result(file)
      import cable_netcdf_io_t, cable_netcdf_file_t
      class(cable_netcdf_io_t), intent(inout) :: this
      character(len=*), intent(in) :: path
        !! Path to the netCDF file to create.
      integer, intent(in) :: iotype
        !! I/O type to use for the file using the `CABLE_NETCDF_IOTYPE_*` constants.
      integer, intent(in), optional :: mode
        !! Optional mode flags for file creation using the `CABLE_NETCDF_MODE_*` constants.
      class(cable_netcdf_file_t), allocatable :: file
    end function
    !> Open an existing netCDF file with the specified path and I/O type.
    function cable_netcdf_io_open_file(this, path, iotype, mode) result(file)
      import cable_netcdf_io_t, cable_netcdf_file_t
      class(cable_netcdf_io_t), intent(inout) :: this
      character(len=*), intent(in) :: path
        !! Path to the netCDF file to open.
      integer, intent(in) :: iotype
        !! I/O type to use for the file using the `CABLE_NETCDF_IOTYPE_*` constants.
      integer, intent(in), optional :: mode
        !! Optional mode flags for file opening using the `CABLE_NETCDF_MODE_*` constants.
      class(cable_netcdf_file_t), allocatable :: file
    end function
    !> Create a new decomposition for parallel I/O.
    !! This follows the degree of freedom approach for specifying I/O
    !! decompositions as described in Denis et al. (2011)
    !! [[10.1177/1094342011428143](https://www.doi.org/10.1177/1094342011428143)].
    function cable_netcdf_io_create_decomp(this, compmap, dims, type) result(decomp)
      import cable_netcdf_io_t, cable_netcdf_decomp_t
      class(cable_netcdf_io_t), intent(inout) :: this
      integer, intent(in) :: compmap(:)
        !* An array of data offsets describing where each element of the
        ! in-memory array is located in the netCDF variable data on disk.
      integer, intent(in) :: dims(:)
        !! An array of the global dimensions used to describe the shape of the data in the netCDF file
      integer, intent(in) :: type
        !! The data type of the in-memory array using the `CABLE_NETCDF_TYPE_*` constants.
      class(cable_netcdf_decomp_t), allocatable :: decomp
    end function
  end interface

  interface
    !> Module initialization procedure for the cable_netcdf_mod module.
    !! This procedure should be called before any other procedures in this module.
    module subroutine cable_netcdf_mod_init(mpi_grp)
      type(mpi_grp_t), intent(in) :: mpi_grp
        !* The MPI group for the set of processes that will be performing netCDF I/O operations.
        ! All procedures in this module should be called collectively by all
        ! processes in this group.
    end subroutine
    !> Module finalization procedure for the cable_netcdf_mod module.
    module subroutine cable_netcdf_mod_end()
    end subroutine
    !> Create a new netCDF file with the specified path and I/O type.
    module function cable_netcdf_create_file(path, iotype, mode) result(file)
      character(len=*), intent(in) :: path
        !! Path to the netCDF file to create.
      integer, intent(in) :: iotype
        !! I/O type to use for the file using the `CABLE_NETCDF_IOTYPE_*` constants.
      integer, intent(in), optional :: mode
        !! Optional mode flags for file creation using the `CABLE_NETCDF_MODE_*` constants.
      class(cable_netcdf_file_t), allocatable :: file
    end function
    !> Open an existing netCDF file with the specified path and I/O type.
    module function cable_netcdf_open_file(path, iotype, mode) result(file)
      character(len=*), intent(in) :: path
        !! Path to the netCDF file to open.
      integer, intent(in) :: iotype
        !! I/O type to use for the file using the `CABLE_NETCDF_IOTYPE_*` constants.
      integer, intent(in), optional :: mode
        !! Optional mode flags for file opening using the `CABLE_NETCDF_MODE_*` constants.
      class(cable_netcdf_file_t), allocatable :: file
    end function
    !> Create a new decomposition for parallel I/O.
    !! This follows the degree of freedom approach for specifying I/O
    !! decompositions as described in Denis et al. (2011)
    !! [[10.1177/1094342011428143](https://www.doi.org/10.1177/1094342011428143)].
    module function cable_netcdf_create_decomp(compmap, dims, type) result(decomp)
      integer, intent(in) :: compmap(:)
        !* An array of data offsets describing where each element of the
        ! in-memory array is located in the netCDF variable data on disk.
      integer, intent(in) :: dims(:)
        !! An array of the global dimensions used to describe the shape of the data in the netCDF file
      integer, intent(in) :: type
        !! The data type of the in-memory array using the `CABLE_NETCDF_TYPE_*` constants.
      class(cable_netcdf_decomp_t), allocatable :: decomp
    end function
  end interface

  !> The internal I/O handler to use for all netCDF file operations.
  class(cable_netcdf_io_t), allocatable, private :: cable_netcdf_io_handler

end module cable_netcdf_mod
