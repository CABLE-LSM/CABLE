! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_mod
  !* This module provides the interface for interacting with the CABLE output system.
  !
  ! The output system is responsible for writing CABLE output variables to one or
  ! more netCDF output and/or restart files, and includes functionality for
  ! performing parallel I/O in MPI mode, grid cell reductions over sub-grid tiles,
  ! and time aggregations of diagnostic variables.
  !
  ! Using the output system involves the following steps:
  !
  ! 1. [[cable_output_mod_init]] must be called before any
  ! other procedures in this module to initialise the output system.
  ! 2. Diagnostics should be registered with the output system via
  ! [[cable_output_register_output_variables]]. This involves creating an array of
  ! `cable_output_variable_t` instances which describe the available diagnostics
  ! and passing this array to `cable_output_register_output_variables`. For
  ! example, a 1-dimensional diagnostic variable defined on the patch
  ! dimension could be registered as follows:
  ! ```fortran
  ! call cable_output_register_output_variables([ &
  !   cable_output_variable_t( &
  !     field_name="my_diagnostic", &
  !     data_shape=[cable_output_get_dimension("patch")], &
  !     aggregator=new_aggregator(my_diagnostic_working_variable) &
  !   ), &
  !   cable_output_variable_t( &
  !     ...
  !   ) &
  ! ])
  ! ```
  ! Note that registering an output variable does not necessarily mean that the
  ! variable will be written to an output stream - this can depend on whether the
  ! output variable is active, which often depends on the output configuration
  ! file, or if the variable is a restart variable and whether we are writing to a
  ! restart file. There are additional properties which may be specified for each
  ! registered output variable - please see [[cable_output_variable_t]] for more
  ! details. In general, output variables should be registered if their associated
  ! diagnostic working variables are initialised in the model as this can help
  ! provide information on the diagnostics which are available.
  ! 3. Output streams should be initialised via [[cable_output_init_streams]]. This
  ! should be done after registering output variables as the output stream
  ! initialisation involves determining which output variables are active in each
  ! output stream based on the current output configuration. Once an output stream
  ! has been initialised, data can be written to disk.
  ! 4. Typically on the first time step of the simulation,
  ! [[cable_output_write_parameters]] should be called to write out any non-time
  ! varying parameter output variables.
  ! 5. On each time step, [[cable_output_update]] should be called to update the
  ! time aggregation accumulation for any output variables that are active in an
  ! output stream. After `cable_output_update` is called, [[cable_output_write]]
  ! should be called to write out the output variables for any output streams with
  ! a sampling frequency that aligns with the current time step.
  ! 6. If writing a CABLE restart file is required, then
  ! [[cable_output_write_restart]] should be called at the end of the simulation
  ! to write out the restart variables to the CABLE restart file.
  ! 7. Lastly, after all output has been written, [[cable_output_mod_end]] should
  ! be called to close any open files and perform any necessary cleanup of
  ! resources.

  use cable_error_handler_mod, only: cable_abort
  use iso_fortran_env, only: int32, real32, real64
  use aggregator_mod, only: aggregator_t
  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_def_types_mod, only: mp
  use cable_def_types_mod, only: mp_global
  use cable_def_types_mod, only: mland
  use cable_def_types_mod, only: mland_global
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncp
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: met_type
  use cable_io_vars_module, only: xdimsize
  use cable_io_vars_module, only: ydimsize
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: patch_type, land_type

  implicit none
  private

  integer, parameter :: CABLE_OUTPUT_VAR_TYPE_UNDEFINED = -1

  !> List of allowed reduction methods for output variables.
  !! Please refer to [[cable_grid_reductions_mod]] for more details on grid reductions.
  character(32), parameter, public :: allowed_reduction_methods(3) = [ &
    "none                    ", &
    "grid_cell_average       ", &
    "first_patch_in_grid_cell" &
  ]

  !> List of allowed aggregation methods for output variables.
  !! Please refer to [[aggregator_mod]] for more details on aggregation methods.
  character(32), parameter, public :: allowed_aggregation_methods(5) = [ &
      "point", &
      "mean ", &
      "max  ", &
      "min  ", &
      "sum  " &
  ]

  !> List of allowed grid types for an output stream.
  character(32), parameter, public :: allowed_grid_types(3) = [ &
      "mask   ", &
      "land   ", &
      "restart" &
  ]

  integer(kind=int32), parameter, public :: CABLE_OUTPUT_FILL_VALUE_INT32  = -9999999_int32
  real(kind=real32),   parameter, public :: CABLE_OUTPUT_FILL_VALUE_REAL32 = -1.0e+33_real32
  real(kind=real64),   parameter, public :: CABLE_OUTPUT_FILL_VALUE_REAL64 = -1.0e+33_real64

  character(64), parameter :: NATIVE_DIM_NAME_PATCH           = "patch_native"
  character(64), parameter :: NATIVE_DIM_NAME_PATCH_GLOBAL    = "patch_global_native"
  character(64), parameter :: NATIVE_DIM_NAME_PATCH_GRID_CELL = "patch_grid_cell_native"
  character(64), parameter :: NATIVE_DIM_NAME_LAND            = "land_native"
  character(64), parameter :: NATIVE_DIM_NAME_LAND_GLOBAL     = "land_global_native"

  type, public :: cable_output_dim_t
    !* Type for describing both in-memory and netCDF variable dimensions used by
    ! the output module.
    !
    ! Instances of `cable_output_dim_t` are created by
    ! [[cable_output_get_dimension]] and is used to describe the in-memory shape
    ! of the native diagnostic of each output variable in
    ! `cable_output_variable_t`.
    !
    ! Components of this type are private to ensure that dimensions are only
    ! created via `cable_output_get_dimension` as several dimension names are
    ! reserved for special handling by the output module. NetCDF variable
    ! dimensions are handled internally in the output module. For more details on
    ! how netCDF variable dimensions are inferred from `cable_output_dim_t`
    ! instances, please refer to [[native_to_netcdf_dimensions]].
    private
    character(64) :: dim_name !! Dimension name.
    integer :: dim_size !! Dimension size.
  contains
    procedure, public :: name => cable_output_dim_get_name !! Return the dimension name.
    procedure, public :: size => cable_output_dim_get_size !! Return the dimension size.
  end type

  type, public :: cable_output_attribute_t
    !! Type for describing string valued netCDF file attributes.
    character(64) :: name !! Name of the attribute.
    character(256) :: value !! Value of the attribute
  end type

  type, public :: cable_output_variable_t
    !* Type for describing output variables.
    !
    ! This type provides the basis for registering output variables with the
    ! output module via [[cable_output_register_output_variables]], and is used in
    ! the definition and writing of output variables in various output streams.
    character(64) :: field_name
      !* The name of the variable as used in the CABLE code. This name is used
      ! as the netCDF variable name when writing CABLE restart files.
    character(64) :: netcdf_name = ""
      !* The name of the variable as it should appear in netCDF output files. If
      ! not specified, this defaults to `field_name`.
    character(64) :: accumulation_frequency = "all"
      !* The frequency at which the variable is accumulated when computing time
      ! aggregations. Please refer to the [[cable_timing_frequency_matches]]
      ! procedure for more information on the available frequency settings. If not
      ! specified, this defaults to "all", meaning that the variable is
      ! accumulated at every CABLE time step.
    character(64) :: reduction_method = "none"
      !* The grid cell reduction method to apply to the variable. The allowed
      ! reduction methods are specified in [[allowed_reduction_methods]]. Please
      ! refer to [[cable_grid_reductions_mod]] for more details on grid
      ! reductions.
    character(64) :: aggregation_method = "point"
      !* The time aggregation method to apply when sampling a diagnostic. Please refer to
      ! [[allowed_aggregation_methods]] for more details on the available
      ! aggregation methods.
    logical :: active = .true.
      !* A flag indicating whether the variable is active in the default output stream.
    logical :: parameter = .false.
      !* A flag indicating whether the variable is a non-time varying parameter.
      ! Variables with `parameter = .true.` are written once on the first time
      ! step via [[cable_output_write_parameters]].
    logical :: distributed = .true.
      !* A flag indicating whether the variable is distributed across multiple
      ! processes. If `distributed = .true.`, the output module will infer an
      ! appropriate parallel I/O decomposition from `data_shape` to perform a
      ! distributed write to disk. If `distributed = .false.`, it is assumed by
      ! the output module that each processes has a copy of the data, and only the
      ! data on the root process will be written.
    logical :: restart = .false.
      !* A flag indicating whether the variable should be written to the CABLE
      ! restart file at the end of the run. Please see
      ! [[cable_output_write_restart]] for more details on how restart variables
      ! are written.
    logical :: patchout = .false.
      !* A flag indicating whether subgrid patch information should be included
      ! in the output variable output. If `patchout = .true.`, this has the same
      ! effect as setting `reduction_method = "none"`. This is a legacy flag for
      ! backward compatibility with the CABLE output namelist settings.
    integer :: var_type = CABLE_OUTPUT_VAR_TYPE_UNDEFINED
      !* The netCDF variable type using `CABLE_NETCDF_<type>` constants. If not
      ! specified, the output module will use the native type of the data as the
      ! netCDF variable type.
    real :: scale_by = 1.0
      !* A multiplicative factor to apply to the native diagnostic values when
      ! writing output.
    real :: divide_by = 1.0
      !* A divisional factor to apply to the native diagnostic values when
      ! writing output.
    real :: offset_by = 0.0
      !* An additive offset to apply to the native diagnostic values when
      ! writing output.
    real, private :: range_native(2) = [-huge(0.0), huge(0.0)]
      !* The valid range of physical values for the output variable in the units
      ! of the native diagnostic.
    real, allocatable :: range(:)
      !* The valid range of physical values for the output variable. If a unit
      ! conversion is applied to the native diagnostic, the range should be given
      ! in the units of the output variable after applying the unit conversion.
      ! If unspecified, all values are considered valid.
    type(cable_output_dim_t), allocatable :: data_shape(:)
      !* An array of in-memory dimensions describing the shape of the variable
      ! data. The dimensions must be created via [[cable_output_get_dimension]]
      ! to ensure that reserved dimension names are handled correctly by the
      ! output module. If not specified, the data shape is assumed to be a
      ! scalar.
    class(aggregator_t), allocatable :: aggregator
      !* The aggregator object associated with the diagnostic working variable
      ! to be written for this output variable. The aggregator object should not
      ! be initialised when registering output variables as this is done
      ! internally in the output module the output variable is active.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !* NetCDF variable attributes to be written with the variable.
  contains
    procedure, private :: get_netcdf_name => cable_output_variable_get_netcdf_name
      !* Return the netCDF variable name, which defaults to `field_name` if not
      ! specified via `netcdf_name`.
  end type

  type :: cable_output_stream_t
    !* Type for describing a netCDF file output stream.
    real :: previous_write_time = 0.0
      !* The simulation time at which the output stream was last written.
    integer :: frame = 0
      !* The current index along the unlimited time dimension for the output stream.
    character(64) :: sampling_frequency
      !* The frequency at which all output variables in the output stream are
      ! aggregated in time and written to disk. Please refer to the
      ! [[cable_timing_frequency_matches]] procedure for more information on the available
      ! frequency settings.
    character(64) :: grid_type
      !* The grid type of the output stream. This controls the netCDF dimensions
      ! and coordinate variables used to describe non-vertical spatial coordinates
      ! in the netCDF file. Common grid types in CABLE include the compressed land
      ! grid, or the lat-lon mask grid. The allowed grid types are specified in
      ! [[allowed_grid_types]].
    character(256) :: file_name
      !* The name of the netCDF file to which the output stream is written.
    class(cable_netcdf_file_t), allocatable :: output_file
      !* The netCDF file object associated with the output stream.
    type(cable_output_variable_t), allocatable :: coordinate_variables(:)
      !* An array of coordinate variables to be written to the output stream.
    type(cable_output_variable_t), allocatable :: output_variables(:)
      !* An array of output variables to be written to the output stream.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !* Global netCDF file attributes to be written to the output stream.
  end type

  public cable_output_mod_init
  interface cable_output_mod_init
    module subroutine cable_output_impl_init()
      !* Module initialisation procedure for `cable_output_mod`.
      !
      ! This procedure must be called before any other procedures in
      ! `cable_output_mod`.
    end subroutine
  end interface

  public cable_output_mod_end
  interface cable_output_mod_end
    module subroutine cable_output_impl_end()
      !* Module finalization procedure for `cable_output_mod`.
      !
      ! This procedure should be called at the end of the simulation after all
      ! output has been written.
    end subroutine
  end interface

  public cable_output_register_output_variables
  interface cable_output_register_output_variables
    module subroutine cable_output_impl_register_output_variables(output_variables)
      !* Registers output variables with the output module. Note that
      ! registering an output variable does not necessarily mean that the variable
      ! will be written to an output stream - this can depend on whether the
      ! output variable is active, or if it is a restart variable. Output
      ! variables should be registered if their associated diagnostic working
      ! variables are initialised in the model as this can help provide the
      ! information on the diagnostics which are available.
      type(cable_output_variable_t), dimension(:), intent(in) :: output_variables
        !! An array of output variable definitions to be registered.
    end subroutine
  end interface

  public cable_output_init_streams
  interface cable_output_init_streams
    module subroutine cable_output_impl_init_streams(dels)
      !! Initialise output streams based on the current output configuration.
      real, intent(in) :: dels !! The current time step size in seconds.
    end subroutine
  end interface

  public cable_output_update
  interface cable_output_update
    module subroutine cable_output_impl_update(time_index, dels, met)
      !* Updates the time aggregation accumulation for any output variables that
      ! are active in an output stream with an accumulation frequency that matches
      ! the current time step.
      integer, intent(in) :: time_index !! The current time step index in the simulation.
      real, intent(in) :: dels !! The current time step size in seconds.
      type(met_type), intent(in) :: met
        !* Met variables at the current time step to provide informative error
        ! messages for CABLE range checks.
    end subroutine
  end interface

  public cable_output_write
  interface cable_output_write
    module subroutine cable_output_impl_write(time_index, dels, met, patch, landpt)
      !* Writes output variables to disk for any output streams with a sampling
      ! frequency that matches the current time step.
      integer, intent(in) :: time_index !! The current time step index in the simulation.
      real, intent(in) :: dels !! The current time step size in seconds.
      type(met_type), intent(in) :: met
        !* Met variables at the current time step to provide informative error
        ! messages for CABLE range checks.
      type(patch_type), intent(in) :: patch(:)
        !! The patch type instance for performing grid reductions over the patch dimension if required.
      type(land_type), intent(in) :: landpt(:)
        !! The land type instance for performing grid reductions over the patch dimension if required.
    end subroutine
  end interface

  public cable_output_write_parameters
  interface cable_output_write_parameters
    module subroutine cable_output_impl_write_parameters(time_index, patch, landpt)
      !* Writes non-time varying parameter output variables to disk. This is
      ! done on the first time step of the simulation after the output streams
      ! have been initialised.
      integer, intent(in) :: time_index !! The current time step index in the simulation.
      type(patch_type), intent(in) :: patch(:)
        !! The patch type instance for performing grid reductions over the patch dimension if required.
      type(land_type), intent(in) :: landpt(:)
        !! The land type instance for performing grid reductions over the patch dimension if required.
    end subroutine
  end interface

  public cable_output_write_restart
  interface cable_output_write_restart
    module subroutine cable_output_impl_write_restart(current_time)
      !* Writes variables to the CABLE restart file. This is done at the end of
      ! the simulation.
      real, intent(in) :: current_time !! Current simulation time
    end subroutine
  end interface

  public cable_output_get_dimension

contains

  function cable_output_get_dimension(name) result(dim)
    !* Returns an output variable dimension. This function contains the
    ! definitions of all dimensions used to describe the in-memory data shapes
    ! of CABLE variables.
    !
    ! @note "Note on adding new dimensions and shapes for output variables"
    ! Adding new dimensions and shapes for output variables is possible, however
    ! it is currently more involved than adding new output variables and requires
    ! making changes to the output module implementation. The steps to add a new
    ! dimension to the output module are as follows:
    !
    ! 1. Add the new dimension name and size definition to `cable_output_get_dimension`.
    ! 2. If grid cell reductions are required for variables involving the new
    ! dimension, add a new grid reduction buffer allocation in
    ! [[cable_output_reductions]] consistent with the data shape and any
    ! necessary code to associate the buffer with an output variable.
    ! 3. If distributed writes are required for variables involving the new
    ! dimension, add a new decomposition definition in `cable_output_decomp_smod`
    ! consistent with the data shape and any necessary code to associate the
    ! decomposition with an output variable.
    !
    ! In future versions this can be improved by generating the necessary grid
    ! reduction buffers and parallel I/O decompositions based on the active output
    ! variables across all output streams, rather than requiring hard coded
    ! definitions for each dimension and shape in the output module implementation.
    ! @endnote
    character(*), intent(in) :: name
      !* Name of the dimension. Please see the implementation of this
      ! function for the list of allowed dimension names and their meanings.
    type(cable_output_dim_t) :: dim
      !! The output dimension object corresponding to the requested dimension name.

    select case(name)
    case ("patch")
      dim = cable_output_dim_t(NATIVE_DIM_NAME_PATCH, mp)
    case ("patch_global")
      dim = cable_output_dim_t(NATIVE_DIM_NAME_PATCH_GLOBAL, mp_global)
    case ("patch_grid_cell")
      dim = cable_output_dim_t(NATIVE_DIM_NAME_PATCH_GRID_CELL, max_vegpatches)
    case ("land")
      dim = cable_output_dim_t(NATIVE_DIM_NAME_LAND, mland)
    case ("land_global")
      dim = cable_output_dim_t(NATIVE_DIM_NAME_LAND_GLOBAL, mland_global)
    case ("soil")
      dim = cable_output_dim_t("soil", ms)
    case ("snow")
      dim = cable_output_dim_t("snow", msn)
    case ("rad")
      dim = cable_output_dim_t("rad", nrb)
    case ("plant_carbon_pools")
      dim = cable_output_dim_t("plant_carbon_pools", ncp)
    case ("soil_carbon_pools")
      dim = cable_output_dim_t("soil_carbon_pools", ncs)
    case ("x")
      dim = cable_output_dim_t("x", xdimsize)
    case ("y")
      dim = cable_output_dim_t("y", ydimsize)
    case default
      call cable_abort("Invalid dimension requested: " // name, __FILE__, __LINE__)
    end select

  end function cable_output_get_dimension

  elemental function cable_output_dim_get_name(this) result(name)
    !! Return the dimension name.
    class(cable_output_dim_t), intent(in) :: this
    character(64) :: name
    name = this%dim_name
  end function

  elemental function cable_output_dim_get_size(this) result(size)
    !! Return the dimension size.
    class(cable_output_dim_t), intent(in) :: this
    integer :: size
    size = this%dim_size
  end function

  elemental function cable_output_variable_get_netcdf_name(this) result(netcdf_name)
    !* Return the netCDF variable name, which defaults to `field_name` if not
    ! specified via `netcdf_name`.
    class(cable_output_variable_t), intent(in) :: this
    character(64) :: netcdf_name
    if (len_trim(this%netcdf_name) > 0) then
      netcdf_name = this%netcdf_name
    else
      netcdf_name = this%field_name
    end if
  end function

end module
