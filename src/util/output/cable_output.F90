! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_mod
  !* This module provides the interface for interacting with the CABLE output system.
  !
  ! The output system writes CABLE output variables to one or more netCDF
  ! output and/or restart files, and includes functionality for performing
  ! parallel I/O in MPI mode, grid cell reductions over sub-grid tiles, and time
  ! aggregations of diagnostic variables.
  !
  ! What is written is described by an output configuration file (YAML) made of
  ! streams (a file and a write frequency) and the variables directed to them.
  ! See the user guide for the file format.
  !
  ! Using the output system involves the following steps:
  !
  ! 1. [[cable_output_mod_init]] must be called before any other procedures in
  ! this module to initialise the output system.
  !
  ! 2. Output variable definitions are built from the output catalogue and the
  ! bindings to model state (see [[cable_output_catalogue_mod]]) and registered
  ! with [[cable_output_register_output_variables]]. A definition
  ! (`cable_output_variable_definition_t`) says what a variable is; it does not
  ! say whether or how it is written. Registering a variable does not mean it
  ! is written: that is decided by the configuration. Variables that the
  ! current model configuration cannot provide are registered as unavailable.
  !
  ! 3. The output configuration is read and checked (see
  ! [[cable_output_config_mod]]), then [[cable_output_init_streams]] creates
  ! the files, defines their contents and writes their coordinate variables.
  ! Every use of a definition (`cable_output_variable_t`) owns its own
  ! aggregator, so one variable can appear in several streams.
  !
  ! 4. Typically on the first time step of the simulation,
  ! [[cable_output_write_parameters]] should be called to write out any
  ! non-time varying parameter output variables.
  !
  ! 5. On each time step, [[cable_output_update]] should be called to update the
  ! time aggregation accumulation for the output variables, followed by
  ! [[cable_output_write]], which writes the streams whose write frequency
  ! matches the current time step.
  !
  ! 6. If writing a CABLE restart file is required, then
  ! [[cable_output_write_restart]] should be called at the end of the simulation
  ! to write out the definitions that have a restart name, as they are.
  !
  ! 7. Lastly, after all output has been written, [[cable_output_mod_end]] should
  ! be called to close any open output streams and perform any necessary cleanup of
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
  character(32), parameter, public :: allowed_reduction_methods(4) = [character(32) :: &
    "none", &
    "grid_cell_average", &
    "first_tile_on_cell", &
    "dominant_tile" &
  ]

  !> List of allowed aggregation methods for output variables.
  !! Please refer to [[aggregator_mod]] for more details on aggregation methods.
  character(32), parameter, public :: allowed_aggregation_methods(5) = [character(32) :: &
    "instant", &
    "mean", &
    "max", &
    "min", &
    "sum" &
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
    !* Type for describing one use of an output variable in an output stream.
    !
    ! Uses are built from a `cable_output_variable_definition_t` and the
    ! aggregation, reduction and NetCDF name chosen in the output configuration,
    ! and are what the output module defines and writes in each output stream.
    ! The constructor taking the individual settings is used for the coordinate
    ! variables written to every file.
    character(64) :: field_name
      !* The name of the variable as used in the CABLE code. This name is used
      ! as the netCDF variable name when writing CABLE restart files.
    character(64) :: netcdf_name = ""
      !* The name of the variable as it should appear in netCDF output files. If
      ! not specified, this defaults to `field_name`.
    character(64) :: accumulation_frequency = "timestep"
      !* The frequency at which the variable is accumulated when computing time
      ! aggregations. Please refer to the [[cable_timing_frequency_matches]]
      ! procedure for more information on the available frequency settings. If not
      ! specified, this defaults to "timestep", meaning that the variable is
      ! accumulated at every CABLE time step.
    character(64) :: reduction_method = "none"
      !* The grid cell reduction method to apply to the variable. The allowed
      ! reduction methods are specified in `allowed_reduction_methods`. Please
      ! refer to [[cable_grid_reductions_mod]] for more details on grid
      ! reductions.
    character(64) :: aggregation_method = "instant"
      !* The time aggregation method to apply when sampling a diagnostic. Please refer to
      ! `allowed_aggregation_methods` for more details on the available
      ! aggregation methods.
    logical :: parameter = .false.
      !* A flag indicating whether the variable is a non-time varying parameter.
      ! Variables with `parameter = .true.` are written once on the first time
      ! step via [[cable_output_write_parameters]].
    logical :: distributed = .true.
      !* A flag indicating whether the variable is distributed across multiple
      ! processes. If `distributed = .true.`, the output module will infer an
      ! appropriate parallel I/O decomposition from `data_shape` to perform a
      ! distributed write to disk. If `distributed = .false.`, it is assumed by
      ! the output module that each process has a copy of the data, and only the
      ! data on the root process will be written.
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
    real :: range(2) = [-huge(0.0), huge(0.0)]
      !* The valid range of physical values for the output variable. If a unit
      ! conversion is applied to the native diagnostic via the `scale_by`,
      ! `divide_by`, or `offset_by` components, the range should be given in the
      ! units of the output variable after applying the unit conversion. If
      ! unspecified, all values are considered valid.
    type(cable_output_dim_t), allocatable :: data_shape(:)
      !* An array of in-memory dimensions describing the shape of the variable
      ! data. The dimensions must be created via [[cable_output_get_dimension]]
      ! to ensure that reserved dimension names are handled correctly by the
      ! output module. If not specified, the data shape is assumed to be a
      ! scalar.
    class(aggregator_t), allocatable :: aggregator
      !* The aggregator object associated with the diagnostic working variable
      ! to be written for this output variable. The aggregator object should not
      ! be initialised by the caller; this is done internally in the output module
      ! when the stream containing the variable is started.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !* NetCDF variable attributes to be written with the variable.
  contains
    procedure, private :: get_netcdf_name => cable_output_variable_get_netcdf_name
      !* Return the netCDF variable name, which defaults to `field_name` if not
      ! specified via `netcdf_name`.
  end type

  interface cable_output_variable_t
    procedure cable_output_variable_constructor
    procedure cable_output_variable_from_definition
  end interface

  type, public :: cable_output_binding_t
    !* Ties an output catalogue name to the CABLE model state it describes.
    !
    ! The catalogue (`output_catalogue.yaml`) describes what a variable is. A
    ! binding supplies the part that only Fortran can: the model variable
    ! (through an aggregator), the unit conversion, the valid range, and whether
    ! the variable exists in the current model configuration.
    character(64) :: name = ""
      !! Catalogue name of the variable this binding belongs to.
    class(aggregator_t), allocatable :: aggregator
      !! Aggregator wrapping the model working variable. Unallocated if the
      !! variable is not available.
    real :: scale_by = 1.0
      !! Multiply the working variable by this to obtain the output value.
    real :: divide_by = 1.0
      !! Divide the working variable by this to obtain the output value.
    real :: offset_by = 0.0
      !! Add this to the working variable to obtain the output value.
    real, allocatable :: range(:)
      !! Valid range of the variable in output units. Unallocated if unchecked.
    logical :: available = .true.
      !! Whether the variable exists in the current model configuration.
  end type

  interface cable_output_binding_t
    procedure cable_output_binding_constructor
  end interface

  type, public :: cable_output_variable_definition_t
    !* Everything about an output variable that does not depend on how it is
    ! written: what it is, where its data comes from and how it is converted.
    ! Aggregation, reduction, frequency, stream and NetCDF name are chosen per
    ! use by the output configuration.
    character(64) :: field_name = ""
      !! User-facing name, also the default NetCDF variable name.
    character(64) :: restart_name = ""
      !! Name used in restart files. Empty if the variable is not restarted.
    logical :: outputtable = .true.
      !! False for variables that are only written to restart files.
    logical :: available = .true.
      !! Whether the variable exists in the current model configuration.
    type(cable_output_dim_t), allocatable :: data_shape(:)
      !! In-memory shape of the working variable. Empty for scalars.
    integer :: var_type = CABLE_OUTPUT_VAR_TYPE_UNDEFINED
      !! NetCDF type of the output.
    logical :: parameter = .false.
      !! True for non-time-varying variables, written once without a time axis.
    logical :: distributed = .true.
      !! False if every process holds a full copy of the data.
    character(64) :: native_frequency = "timestep"
      !! How often the model itself updates the working variable.
    character(64) :: group_name = ""
      !! Convenience group the variable belongs to. Empty for none.
    character(64) :: module_name = ""
      !! Convenience module the variable belongs to. Empty for none.
    real :: scale_by = 1.0
    real :: divide_by = 1.0
    real :: offset_by = 0.0
    real :: range(2) = [-huge(0.0), huge(0.0)]
      !! Valid range in the units of the working variable.
    class(aggregator_t), allocatable :: aggregator
      !! Aggregator wrapping the working variable. Copied for each use.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !! NetCDF variable attributes, e.g. `units` and `long_name`.
  end type

  integer, parameter, public :: CABLE_OUTPUT_CONFIG_LENGTH = 256
  character(len=*), parameter, public :: CABLE_OUTPUT_DEFAULT_NETCDF_NAME = "{field_name}"

  type, public :: cable_output_stream_config_t
    !! Settings of one output stream.
    integer :: stream_id = 0
    character(len=CABLE_OUTPUT_CONFIG_LENGTH) :: file_name = ""
    character(len=16) :: frequency = ""
    character(len=CABLE_OUTPUT_CONFIG_LENGTH) :: netcdf_name = CABLE_OUTPUT_DEFAULT_NETCDF_NAME
      !! Template for NetCDF variable names in this stream.
    logical :: shuffle = .true.
    integer :: compression_level = 1
      !! Deflate level 0 to 9. Zero means no compression.
    logical :: separate_file_per_variable = .false.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !! Global attributes, with templates filled in.
  end type cable_output_stream_config_t

  type, public :: cable_output_variable_config_t
    !! One requested use of an output variable, after expansion of groups and modules.
    character(len=64) :: field_name = ""
    integer :: stream_id = 0
    character(len=16) :: aggregation = ""
    character(len=32) :: reduction = "none"
    character(len=CABLE_OUTPUT_CONFIG_LENGTH) :: netcdf_name = ""
      !! NetCDF variable name, with templates filled in.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !! Variable attributes from the file, with templates filled in.
  end type cable_output_variable_config_t

  type, public :: cable_output_config_t
    type(cable_output_stream_config_t), allocatable :: streams(:)
    type(cable_output_variable_config_t), allocatable :: variables(:)
  end type cable_output_config_t

  type, public :: cable_output_stream_t
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
      ! `allowed_grid_types`.
    character(256) :: file_name
      !* The name of the netCDF file to which the output stream is written.
    logical :: shuffle = .false.
      !! Whether to apply the shuffle filter to the variables in the file.
    integer :: compression_level = 0
      !! Deflate level 0 to 9 for the variables in the file. Zero means no compression.
    class(cable_netcdf_file_t), allocatable :: output_file
      !* The netCDF file object associated with the output stream.
    type(cable_output_variable_t), allocatable :: coordinate_variables(:)
      !* An array of coordinate variables to be written to the output stream.
    type(cable_output_variable_t), allocatable :: output_variables(:)
      !* An array of output variables to be written to the output stream.
    type(cable_output_attribute_t), allocatable :: metadata(:)
      !* Global netCDF file attributes to be written to the output stream.
  end type

  public :: cable_output_restart_variable
  public :: cable_output_grid_type

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
    module subroutine cable_output_impl_register_definitions(definitions)
      !* Registers output variable definitions with the output module. A
      ! definition says what a variable is; whether and how it is written is
      ! decided by the output configuration passed to [[cable_output_init_streams]].
      type(cable_output_variable_definition_t), dimension(:), intent(in) :: definitions
        !! The output variable definitions to register.
    end subroutine
  end interface

  public cable_output_build_streams
  interface cable_output_build_streams
    module function cable_output_impl_build_streams(config, definitions, grid_type) result(streams)
      !* Works out which files the configuration produces and which variables
      ! go into each, without creating any files.
      type(cable_output_config_t), intent(in) :: config
        !! The output configuration.
      type(cable_output_variable_definition_t), intent(in) :: definitions(:)
        !! The output variable definitions the configuration refers to.
      character(len=*), intent(in) :: grid_type
        !! The grid type of every output file, one of `allowed_grid_types`.
      type(cable_output_stream_t), allocatable :: streams(:)
    end function
  end interface

  public cable_output_init_streams
  interface cable_output_init_streams
    module subroutine cable_output_impl_init_streams_from_config(config, grid_type, dels)
      !* Creates the output files described by the configuration, defines
      ! their contents and writes their coordinate variables. Definitions must
      ! have been registered first.
      type(cable_output_config_t), intent(in) :: config
        !! The output configuration.
      character(len=*), intent(in) :: grid_type
        !! The grid type of every output file, one of `allowed_grid_types`.
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
    ! @note "Adding new dimensions"
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

  function cable_output_variable_constructor(field_name, aggregator, netcdf_name, &
    accumulation_frequency, reduction_method, aggregation_method, &
    parameter, distributed, var_type, scale_by, divide_by, &
    offset_by, range, data_shape, metadata &
  ) result(this)
    !! Constructor for `cable_output_variable_t`.
    character(*), intent(in) :: field_name
    class(aggregator_t), intent(in) :: aggregator
    character(*), intent(in), optional :: netcdf_name
    character(*), intent(in), optional :: accumulation_frequency
    character(*), intent(in), optional :: reduction_method
    character(*), intent(in), optional :: aggregation_method
    logical, intent(in), optional :: parameter
    logical, intent(in), optional :: distributed
    integer, intent(in), optional :: var_type
    real, intent(in), optional :: scale_by
    real, intent(in), optional :: divide_by
    real, intent(in), optional :: offset_by
    real, intent(in), optional :: range(:)
    type(cable_output_dim_t), intent(in), optional :: data_shape(:)
    type(cable_output_attribute_t), intent(in), optional :: metadata(:)
    type(cable_output_variable_t) :: this

    this%field_name = field_name
    allocate(this%aggregator, source=aggregator)
    if (present(netcdf_name)) this%netcdf_name = netcdf_name
    if (present(accumulation_frequency)) this%accumulation_frequency = accumulation_frequency
    if (present(reduction_method)) this%reduction_method = reduction_method
    if (present(aggregation_method)) this%aggregation_method = aggregation_method
    if (present(parameter)) this%parameter = parameter
    if (present(distributed)) this%distributed = distributed
    if (present(var_type)) this%var_type = var_type
    if (present(scale_by)) this%scale_by = scale_by
    if (present(divide_by)) this%divide_by = divide_by
    if (present(offset_by)) this%offset_by = offset_by
    if (present(range)) then
      ! Convert range to native units for comparison with working variable:
      this%range = (range - this%offset_by) * this%divide_by / this%scale_by
    end if
    if (present(data_shape)) this%data_shape = data_shape
    if (present(metadata)) this%metadata = metadata

  end function cable_output_variable_constructor

  function cable_output_binding_constructor(name, aggregator, scale_by, divide_by, offset_by, range, available) result(this)
    !! Constructor for `cable_output_binding_t`.
    character(*), intent(in) :: name
    class(aggregator_t), intent(in), optional :: aggregator
      !! Absent for variables that do not exist in the current model configuration.
    real, intent(in), optional :: scale_by
    real, intent(in), optional :: divide_by
    real, intent(in), optional :: offset_by
    real, intent(in), optional :: range(:)
    logical, intent(in), optional :: available
    type(cable_output_binding_t) :: this

    this%name = name
    ! `source=` makes a copy of the aggregator that was passed in, so the
    ! binding owns its own object. A binding for an unavailable variable is
    ! built without one.
    if (present(aggregator)) allocate(this%aggregator, source=aggregator)
    ! Every other setting is optional and keeps the default declared in the
    ! type (scale 1, divide 1, offset 0, no range check, available) if absent.
    if (present(scale_by)) this%scale_by = scale_by
    if (present(divide_by)) this%divide_by = divide_by
    if (present(offset_by)) this%offset_by = offset_by
    if (present(range)) this%range = range
    if (present(available)) this%available = available
  end function cable_output_binding_constructor

  function cable_output_variable_from_definition(definition, aggregation, reduction, netcdf_name, metadata) result(this)
    !* Constructs the use of a variable definition described by the arguments:
    ! the definition says what the variable is, the arguments say how it is
    ! written. The result owns its own copy of the aggregator, so one definition
    ! can be used in several streams.
    type(cable_output_variable_definition_t), intent(in) :: definition
    character(*), intent(in) :: aggregation !! One of `allowed_aggregation_methods`.
    character(*), intent(in) :: reduction !! One of `allowed_reduction_methods`.
    character(*), intent(in) :: netcdf_name !! Name of the variable in the NetCDF file.
    type(cable_output_attribute_t), intent(in), optional :: metadata(:)
      !! Attributes to add to those of the definition.
    type(cable_output_variable_t) :: this
    type(cable_output_attribute_t), allocatable :: attributes(:)
    integer :: definition_count, user_count

    ! Part 1: settings chosen for this use (the arguments), plus the settings
    ! that are a fixed property of the variable (taken from the definition).
    this%field_name = definition%field_name
    this%netcdf_name = netcdf_name
    ! How often the variable is sampled is a property of the model, not of the
    ! file, so the definition's native frequency becomes the sampling frequency.
    this%accumulation_frequency = definition%native_frequency
    this%reduction_method = reduction
    this%aggregation_method = aggregation
    this%parameter = definition%parameter
    this%distributed = definition%distributed
    this%var_type = definition%var_type
    this%scale_by = definition%scale_by
    this%divide_by = definition%divide_by
    this%offset_by = definition%offset_by
    this%range = definition%range
    if (allocated(definition%data_shape)) this%data_shape = definition%data_shape
    ! Copy (not share) the aggregator. It holds the running totals, and two
    ! uses of one variable in different streams must not add into the same totals.
    ! The copy still points at the same model variable, which is what we want.
    allocate(this%aggregator, source=definition%aggregator)

    ! Part 2: the NetCDF attributes, in a fixed order: the definition's own
    ! (units, long_name), then cell_methods, then any the user added.
    definition_count = 0
    user_count = 0
    if (allocated(definition%metadata)) definition_count = size(definition%metadata)
    if (present(metadata)) user_count = size(metadata)
    allocate(attributes(definition_count))
    if (definition_count > 0) attributes = definition%metadata
    ! cell_methods describes how the data was averaged, so it only applies to
    ! variables that change with time. It is worked out from the aggregation
    ! and reduction chosen, never stored, so it cannot disagree with the data.
    if (.not. definition%parameter) then
      attributes = [attributes, cable_output_attribute_t("cell_methods", cell_methods(aggregation, reduction))]
    end if
    ! Appending builds a longer array from the old one plus the new items.
    if (user_count > 0) attributes = [attributes, metadata]
    this%metadata = attributes
  end function cable_output_variable_from_definition

  function cable_output_restart_variable(definition) result(this)
    !* Constructs the variable that writes a definition to the restart file:
    ! the working variable is written as it is, without aggregation or
    ! reduction, under its restart name.
    type(cable_output_variable_definition_t), intent(in) :: definition
    type(cable_output_variable_t) :: this

    ! The restart file uses the variable's restart name, which may differ from
    ! its user-facing name so that restart files stay readable across versions.
    ! No aggregation or reduction is set: the write routine takes the current
    ! model value as it is when it is told this is a restart write.
    this%field_name = definition%restart_name
    this%distributed = definition%distributed
    this%var_type = definition%var_type
    if (allocated(definition%data_shape)) this%data_shape = definition%data_shape
    if (allocated(definition%metadata)) this%metadata = definition%metadata
    allocate(this%aggregator, source=definition%aggregator)
  end function cable_output_restart_variable

  function cable_output_grid_type(output_grid, met_grid) result(grid_type)
    !* The grid type of the output files, from the `output%grid` setting and
    ! the grid type of the meteorological forcing.
    character(len=*), intent(in) :: output_grid !! `default`, `land`, `mask` or `ALMA`.
    character(len=*), intent(in) :: met_grid !! `land` or `mask`.
    character(32) :: grid_type

    ! "default" means follow the forcing: a land-point file gives compressed
    ! land-point output, a lat/lon file gives lat/lon ("mask") output. "land",
    ! "mask" and "ALMA" force the choice (ALMA is a convention for gridded
    ! output, which is the lat/lon layout).
    if (output_grid == "land" .or. (output_grid == "default" .and. met_grid == "land")) then
      grid_type = "land"
    else if ((output_grid == "default" .and. met_grid == "mask") .or. output_grid == "mask" .or. output_grid == "ALMA") then
      grid_type = "mask"
    else
      ! Unrecognised combination, for example a forcing grid that is neither
      ! land nor mask.
      call cable_abort("Unable to determine output grid type.", __FILE__, __LINE__)
    end if
  end function cable_output_grid_type

  pure function cell_methods(aggregation, reduction) result(text)
    !! The CF `cell_methods` attribute for a variable aggregated and reduced as given.
    character(*), intent(in) :: aggregation
    character(*), intent(in) :: reduction
    character(:), allocatable :: text

    ! CF writes cell_methods as "dimension: method" pairs, spatial method first.
    ! Only an area average is worth stating; picking one tile's value (first or
    ! dominant) is not an averaging method.
    text = ""
    if (reduction == "grid_cell_average") text = "area: mean "
    ! The time part uses CF's own words, so "max" becomes "maximum" and
    ! "instant" becomes "point" (a value at one instant).
    select case (aggregation)
    case ("mean")
      text = text // "time: mean"
    case ("sum")
      text = text // "time: sum"
    case ("max")
      text = text // "time: maximum"
    case ("min")
      text = text // "time: minimum"
    case default
      text = text // "time: point"
    end select
  end function cell_methods

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
