module cable_output_types_mod

  use iso_fortran_env, only: int32, real32, real64

  use aggregator_mod, only: aggregator_t

  use cable_netcdf_mod, only: cable_netcdf_file_t

  use cable_enum_mod, only: cable_enum_t

  implicit none
  private

  type, extends(cable_enum_t), public :: cable_output_dim_t
  end type

  type, public :: cable_output_attribute_t
    character(len=64) :: name
    character(len=256) :: value
  end type

  type, public :: cable_output_variable_t
    character(len=64)  :: name
    character(len=64)  :: accumulation_frequency = "all"
    character(len=64)  :: reduction_method = "none"
    character(len=64)  :: aggregation_method
    logical :: active
    logical :: parameter = .false.
    logical :: distributed = .true.
    logical :: restart = .false.
    logical :: patchout = .false.
    integer :: var_type
    real, dimension(2) :: range
    type(cable_output_dim_t), allocatable :: data_shape(:)
    class(aggregator_t), allocatable :: aggregator
    type(cable_output_attribute_t), allocatable :: metadata(:)
  end type

  type, public :: cable_output_profile_t
    real :: previous_write_time = 0.0
    integer :: frame = 0
    character(len=64) :: sampling_frequency
    character(len=64) :: grid_type
    character(len=256) :: file_name
    class(cable_netcdf_file_t), allocatable :: output_file
    type(cable_output_variable_t), allocatable :: output_variables(:)
    type(cable_output_attribute_t), allocatable :: metadata(:)
  end type

  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_PATCH       = cable_output_dim_t(0)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_SOIL        = cable_output_dim_t(1)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_SNOW        = cable_output_dim_t(2)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_RAD         = cable_output_dim_t(3)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_PLANTCARBON = cable_output_dim_t(4)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_SOILCARBON  = cable_output_dim_t(5)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_LAND        = cable_output_dim_t(6)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_LAND_GLOBAL = cable_output_dim_t(7)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_X           = cable_output_dim_t(8)
  type(cable_output_dim_t), parameter, public :: CABLE_OUTPUT_DIM_Y           = cable_output_dim_t(9)

  integer(kind=int32), parameter, public :: FILL_VALUE_INT32  = -9999999_int32
  real(kind=real32),   parameter, public :: FILL_VALUE_REAL32 = -1.0e+33_real32
  real(kind=real64),   parameter, public :: FILL_VALUE_REAL64 = -1.0e+33_real64

end module
