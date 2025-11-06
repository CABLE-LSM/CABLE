! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

! TODO(Sean): The preprocessor define ENFORCE_SINGLE_PRECISION is enabled
! temporarily to restore bitwise reproducibility with the previous output module
! which enforces double precision variables to be sampled using single precision
! arrays.
#define ENFORCE_SINGLE_PRECISION

module aggregator_mod
  !* This module defines the `aggregator_t` type and its extensions, which are
  ! used to perform various types of time aggregations (e.g., mean, sum, min, max)
  ! on source data arrays.
  !
  ! Aggregators are utilised by first initialising an aggregator via
  ! [[new_aggregator]] which accepts as its argument a source data array to be
  ! sampled. The source data array argument is then associated with a pointer
  ! which is used to sample the source data array to compute the time aggregated
  ! value. The source data array argument must be an allocated pointer or
  ! allocatable array. The source data array argument is typically a CABLE working
  ! variable whos memory remains allocated throughout the duration of the
  ! simulation. Aggregator instances can be represented using a polymorphic
  ! `aggregator_t` type or a concrete type corresponding to the data type and rank
  ! of the source data array (e.g., `aggregator_real32_2d_t` for a 2D 32-bit real
  ! source data array). For example:
  !
  ! ```fortran
  ! real(kind=real32), dimension(:,:), allocatable :: source_data_array
  ! class(aggregator_t), allocatable :: aggregator_polymorphic
  ! type(aggregator_real32_2d_t) :: aggregator_concrete
  !
  ! allocate(source_data_array(42, 42))
  !
  ! aggregator_polymorphic = new_aggregator(source_data_array)
  ! aggregator_concrete = new_aggregator(source_data_array)
  ! ```
  !
  ! Aggregator instances must then be initialised by calling their `init`
  ! method, which will allocate the necessary memory for the aggregated data and
  ! set the accumulation and reset methods according to the specified aggregation
  ! method (e.g., mean, sum, min, max). For example:
  !
  ! ```fortran
  ! call aggregator_polymorphic%init('mean')
  ! call aggregator_concrete%init('sum')
  ! ```
  !
  ! Once initialised, the aggregator can be accumulated any number of times
  ! throughout a simulation by calling its `accumulate` method, which will update
  ! the aggregated data from the source data using the appropriate aggregation
  ! method. Once the aggregator has accumulated enough values over the time
  ! interval of interest, the aggregated data can then be accessed via the
  ! `aggregated_data` component of the specific aggregator type (e.g.,
  ! `aggregator_real32_2d_t%aggregated_data`). To reset the aggregator back to its
  ! initial state for the next time interval, the `reset` method can be called,
  ! which will reset the aggregator according to the specified aggregation method.
  !
  ! Resources allocated by the aggregator on initialisation will automatically be
  ! deallocated when the aggregator instance goes out of scope.

  use iso_fortran_env, only: int32, real32, real64
  use cable_error_handler_mod, only: cable_abort

  implicit none
  private

  public :: aggregator_t
  public :: aggregator_int32_0d_t
  public :: aggregator_int32_1d_t
  public :: aggregator_int32_2d_t
  public :: aggregator_int32_3d_t
  public :: aggregator_real32_0d_t
  public :: aggregator_real32_1d_t
  public :: aggregator_real32_2d_t
  public :: aggregator_real32_3d_t
  public :: aggregator_real64_0d_t
  public :: aggregator_real64_1d_t
  public :: aggregator_real64_2d_t
  public :: aggregator_real64_3d_t
  public :: new_aggregator

  type, abstract :: aggregator_t
    !* The `aggregator_t` type is an abstract base type for performing time
    ! aggregations on source data arrays. Support for different data types and
    ! array ranks is provided through extensions of this type.
    integer :: counter = 0 !! The number of times the source data has been accumulated.
    procedure(accumulate_data), pointer :: accumulate
      !* A procedure pointer to the accumulation method used to compute the time
      ! aggregated quantity. The specific method (e.g., mean, sum, min, max) is
      ! set by the `set_method` procedure.
    procedure(reset_data), pointer :: reset
      !* A procedure pointer to the reset method used to reset the aggregated data. The
      ! specific method (e.g., mean, sum, min, max) is set by the `set_method` procedure.
  contains
    procedure :: init => aggregator_init !! Initialise the aggregator.
    procedure :: type => aggregator_type !! Return a string identifier of the aggregator type.
    procedure :: rank => aggregator_rank !! Return the rank of the aggregator.
    procedure :: shape => aggregator_shape !! Return the shape of the aggregator.
    procedure :: scale => aggregator_scale !! Scale the aggregated data by a specified factor.
    procedure :: div => aggregator_div !! Divide the aggregated data by a specified factor.
    procedure :: offset => aggregator_offset !! Add a specified offset to the aggregated data.
    procedure, private :: set_method => aggregator_set_method !! Set the aggregation method.
  end type aggregator_t

  abstract interface
    !* Interfaces for the procedure pointers in the `aggregator_t` type to be
    ! implemented by the specific aggregation methods (e.g., mean, sum, min, max).
    subroutine accumulate_data(this, scale, div, offset)
      !! Accumulate the aggregated data from the source data.
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
      real, intent(in), optional :: scale
        !* An optional scaling factor to apply to the source data before
        ! accumulation. Defaults to 1.0 if not provided.
      real, intent(in), optional :: div
        !* An optional division factor to apply to the source data before
        ! accumulation. Defaults to 1.0 if not provided.
      real, intent(in), optional :: offset
        !* An optional offset to add to the source data before accumulation.
        ! Defaults to 0.0 if not provided.
    end subroutine accumulate_data
    subroutine reset_data(this)
      !! Reset the aggregated data to its initial state.
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
    end subroutine reset_data
  end interface

  type, extends(aggregator_t) :: aggregator_int32_0d_t
    !! An aggregator for 0-dimensional (scalar) 32-bit integer data.
    integer(kind=int32), allocatable :: aggregated_data
    integer(kind=int32), pointer :: source_data => null()
  end type aggregator_int32_0d_t

  type, extends(aggregator_t) :: aggregator_int32_1d_t
    !! An aggregator for 1-dimensional 32-bit integer data.
    integer(kind=int32), dimension(:), allocatable :: aggregated_data
    integer(kind=int32), dimension(:), pointer :: source_data => null()
  end type aggregator_int32_1d_t

  type, extends(aggregator_t) :: aggregator_int32_2d_t
    !! An aggregator for 2-dimensional 32-bit integer data.
    integer(kind=int32), dimension(:,:), allocatable :: aggregated_data
    integer(kind=int32), dimension(:,:), pointer :: source_data => null()
  end type aggregator_int32_2d_t

  type, extends(aggregator_t) :: aggregator_int32_3d_t
    !! An aggregator for 3-dimensional 32-bit integer data.
    integer(kind=int32), dimension(:,:,:), allocatable :: aggregated_data
    integer(kind=int32), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_int32_3d_t

  type, extends(aggregator_t) :: aggregator_real32_0d_t
    !! An aggregator for 0-dimensional (scalar) 32-bit real data.
    real(kind=real32), allocatable :: aggregated_data
    real(kind=real32), pointer :: source_data => null()
  end type aggregator_real32_0d_t

  type, extends(aggregator_t) :: aggregator_real32_1d_t
    !! An aggregator for 1-dimensional 32-bit real data.
    real(kind=real32), dimension(:), allocatable :: aggregated_data
    real(kind=real32), dimension(:), pointer :: source_data => null()
  end type aggregator_real32_1d_t

  type, extends(aggregator_t) :: aggregator_real32_2d_t
    !! An aggregator for 2-dimensional 32-bit real data.
    real(kind=real32), dimension(:,:), allocatable :: aggregated_data
    real(kind=real32), dimension(:,:), pointer :: source_data => null()
  end type aggregator_real32_2d_t

  type, extends(aggregator_t) :: aggregator_real32_3d_t
    !! An aggregator for 3-dimensional 32-bit real data.
    real(kind=real32), dimension(:,:,:), allocatable :: aggregated_data
    real(kind=real32), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_real32_3d_t

  type, extends(aggregator_t) :: aggregator_real64_0d_t
    !! An aggregator for 0-dimensional (scalar) 64-bit real data.
#ifdef ENFORCE_SINGLE_PRECISION
    real(kind=real32), allocatable :: aggregated_data
#else
    real(kind=real64), allocatable :: aggregated_data
#endif

    real(kind=real64), pointer :: source_data => null()
  end type aggregator_real64_0d_t

  type, extends(aggregator_t) :: aggregator_real64_1d_t
    !! An aggregator for 1-dimensional 64-bit real data.
#ifdef ENFORCE_SINGLE_PRECISION
    real(kind=real32), dimension(:), allocatable :: aggregated_data
#else
    real(kind=real64), dimension(:), allocatable :: aggregated_data
#endif
    real(kind=real64), dimension(:), pointer :: source_data => null()
  end type aggregator_real64_1d_t

  type, extends(aggregator_t) :: aggregator_real64_2d_t
    !! An aggregator for 2-dimensional 64-bit real data.
#ifdef ENFORCE_SINGLE_PRECISION
    real(kind=real32), dimension(:,:), allocatable :: aggregated_data
#else
    real(kind=real64), dimension(:,:), allocatable :: aggregated_data
#endif
    real(kind=real64), dimension(:,:), pointer :: source_data => null()
  end type aggregator_real64_2d_t

  type, extends(aggregator_t) :: aggregator_real64_3d_t
    !! An aggregator for 3-dimensional 64-bit real data.
#ifdef ENFORCE_SINGLE_PRECISION
    real(kind=real32), dimension(:,:,:), allocatable :: aggregated_data
#else
    real(kind=real64), dimension(:,:,:), allocatable :: aggregated_data
#endif
    real(kind=real64), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_real64_3d_t

  interface new_aggregator
    !* Factory interface for creating new aggregator instances. The specific
    ! type of aggregator created is determined by the type of the source data
    ! array provided.
    module procedure new_aggregator_int32_0d_t
    module procedure new_aggregator_int32_1d_t
    module procedure new_aggregator_int32_2d_t
    module procedure new_aggregator_int32_3d_t
    module procedure new_aggregator_real32_0d
    module procedure new_aggregator_real32_1d
    module procedure new_aggregator_real32_2d
    module procedure new_aggregator_real32_3d
    module procedure new_aggregator_real64_0d
    module procedure new_aggregator_real64_1d
    module procedure new_aggregator_real64_2d
    module procedure new_aggregator_real64_3d
  end interface

contains

  subroutine aggregator_init(this, method)
    !* Initialise the aggregator by allocating the aggregated data array and its
    ! aggregation method. The values in the aggregated data array are reset
    ! according to the specified aggregation method.
    class(aggregator_t), intent(inout) :: this
    character(len=*), intent(in) :: method
      !! The aggregation method to use (e.g., "mean", "sum", "point", "min", "max").

    select type (this)
    type is (aggregator_int32_0d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_int32_1d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_int32_2d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_int32_3d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real32_0d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real32_1d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real32_2d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real32_3d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real64_0d_t)
      if (.not. allocated(this%aggregated_data)) allocate( &
#ifdef ENFORCE_SINGLE_PRECISION
        this%aggregated_data, mold=real(this%source_data, kind=real32) &
#else
        this%aggregated_data, mold=real(this%source_data, kind=real64) &
#endif
      )
    type is (aggregator_real64_1d_t)
      if (.not. allocated(this%aggregated_data)) allocate( &
#ifdef ENFORCE_SINGLE_PRECISION
        this%aggregated_data, mold=real(this%source_data, kind=real32) &
#else
        this%aggregated_data, mold=real(this%source_data, kind=real64) &
#endif
      )
    type is (aggregator_real64_2d_t)
      if (.not. allocated(this%aggregated_data)) allocate( &
#ifdef ENFORCE_SINGLE_PRECISION
        this%aggregated_data, mold=real(this%source_data, kind=real32) &
#else
        this%aggregated_data, mold=real(this%source_data, kind=real64) &
#endif
      )
    type is (aggregator_real64_3d_t)
      if (.not. allocated(this%aggregated_data)) allocate( &
#ifdef ENFORCE_SINGLE_PRECISION
        this%aggregated_data, mold=real(this%source_data, kind=real32) &
#else
        this%aggregated_data, mold=real(this%source_data, kind=real64) &
#endif
      )
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    call this%set_method(method)

    call this%reset()

  end subroutine aggregator_init

  subroutine aggregator_set_method(this, method)
    !* Set the aggregation method for the aggregator by assigning the appropriate
    ! accumulation and reset procedures based on the specified method.
    class(aggregator_t), intent(inout) :: this
    character(len=*), intent(in) :: method
      !! The aggregation method to use (e.g., "mean", "sum", "point", "min", "max").

    if (method == "mean") then
      this%accumulate => mean_accumulate
      this%reset => other_reset
    elseif (method == "sum") then
      this%accumulate => sum_accumulate
      this%reset => other_reset
    elseif (method == "point") then
      this%accumulate => point_accumulate
      this%reset => point_reset
    elseif (method == "min") then
      this%accumulate => min_accumulate
      this%reset => min_reset
    elseif (method == "max") then
      this%accumulate => max_accumulate
      this%reset => max_reset
    else
      call cable_abort("Aggregation method "//method//" is invalid.", file=__FILE__, line=__LINE__)
    endif

  end subroutine aggregator_set_method

  character(16) function aggregator_type(this)
    !! Return a string identifier of the aggregator type (e.g., "int32", "real32", "real64").
    class(aggregator_t), intent(in) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      aggregator_type = "int32"
    type is (aggregator_int32_1d_t)
      aggregator_type = "int32"
    type is (aggregator_int32_2d_t)
      aggregator_type = "int32"
    type is (aggregator_int32_3d_t)
      aggregator_type = "int32"
    type is (aggregator_real32_0d_t)
      aggregator_type = "real32"
    type is (aggregator_real32_1d_t)
      aggregator_type = "real32"
    type is (aggregator_real32_2d_t)
      aggregator_type = "real32"
    type is (aggregator_real32_3d_t)
      aggregator_type = "real32"
    type is (aggregator_real64_0d_t)
      aggregator_type = "real64"
    type is (aggregator_real64_1d_t)
      aggregator_type = "real64"
    type is (aggregator_real64_2d_t)
      aggregator_type = "real64"
    type is (aggregator_real64_3d_t)
      aggregator_type = "real64"
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

  end function aggregator_type

  integer function aggregator_rank(this)
    !! Return the rank of the aggregator.
    class(aggregator_t), intent(in) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      aggregator_rank = 0
    type is (aggregator_int32_1d_t)
      aggregator_rank = 1
    type is (aggregator_int32_2d_t)
      aggregator_rank = 2
    type is (aggregator_int32_3d_t)
      aggregator_rank = 3
    type is (aggregator_real32_0d_t)
      aggregator_rank = 0
    type is (aggregator_real32_1d_t)
      aggregator_rank = 1
    type is (aggregator_real32_2d_t)
      aggregator_rank = 2
    type is (aggregator_real32_3d_t)
      aggregator_rank = 3
    type is (aggregator_real64_0d_t)
      aggregator_rank = 0
    type is (aggregator_real64_1d_t)
      aggregator_rank = 1
    type is (aggregator_real64_2d_t)
      aggregator_rank = 2
    type is (aggregator_real64_3d_t)
      aggregator_rank = 3
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

  end function aggregator_rank

  function aggregator_shape(this) result(agg_shape)
    !! Return the shape of the aggregator.
    class(aggregator_t), intent(in) :: this
    integer, allocatable :: agg_shape(:)

    select type (this)
    type is (aggregator_int32_0d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_int32_1d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_int32_2d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_int32_3d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real32_0d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real32_1d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real32_2d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real32_3d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real64_0d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real64_1d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real64_2d_t)
      agg_shape = shape(this%source_data)
    type is (aggregator_real64_3d_t)
      agg_shape = shape(this%source_data)
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

  end function aggregator_shape

  subroutine aggregator_scale(this, scale)
    !! Scale the aggregated data by a specified factor.
    class(aggregator_t), intent(inout) :: this
    real, intent(in) :: scale !! The factor by which to scale the aggregated data.

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_int32_1d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_int32_2d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_int32_3d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data * scale
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data * scale
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

  end subroutine aggregator_scale

  subroutine aggregator_div(this, div)
    !! Divide the aggregated data by a specified factor.
    class(aggregator_t), intent(inout) :: this
    real, intent(in) :: div !! The factor by which to divide the aggregated data.

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_int32_1d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_int32_2d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_int32_3d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data / div
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data / div
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

  end subroutine aggregator_div

  subroutine aggregator_offset(this, offset)
    !! Offset the aggregated data by a specified value.
    class(aggregator_t), intent(inout) :: this
    real, intent(in) :: offset !! The value by which to offset the aggregated data.

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_int32_1d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_int32_2d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_int32_3d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data + offset
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data + offset
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

  end subroutine aggregator_offset

  subroutine get_accumulate_args(scale, div, offset, scale_out, div_out, offset_out)
    !* Helper subroutine to get initialise optional scale, div, and offset arguments for
    ! accumulate procedures.
    real, intent(in), optional :: scale
      !* An optional scaling factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: div
      !* An optional division factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: offset
      !* An optional offset to add to the source data before accumulation.
      ! Defaults to 0.0 if not provided.
    real :: scale_out, div_out, offset_out

    if (present(scale)) then
      scale_out = scale
    else
      scale_out = 1.0
    end if

    if (present(div)) then
      div_out = div
    else
      div_out = 1.0
    end if

    if (present(offset)) then
      offset_out = offset
    else
      offset_out = 0.0
    end if

  end subroutine

  subroutine mean_accumulate(this, scale, div, offset)
    !* Accumulate the aggregated data from the source data using the mean
    ! aggregation method.
    class(aggregator_t), intent(inout) :: this
    real, intent(in), optional :: scale
      !* An optional scaling factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: div
      !* An optional division factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: offset
      !* An optional offset to add to the source data before accumulation.
      ! Defaults to 0.0 if not provided.
    real :: scale_val, div_val, offset_val

    call get_accumulate_args(scale, div, offset, scale_val, div_val, offset_val)

    select type (this)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div + offset_val - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div + offset_val - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div + offset_val - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div + offset_val - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data + ( &
        real( &
#ifdef ENFORCE_SINGLE_PRECISION
          scale_val * this%source_data / div + offset_val, kind=real32 &
#else
          scale_val * this%source_data / div + offset_val, kind=real64 &
#endif
        ) - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data + ( &
        real( &
#ifdef ENFORCE_SINGLE_PRECISION
          scale_val * this%source_data / div + offset_val, kind=real32 &
#else
          scale_val * this%source_data / div + offset_val, kind=real64 &
#endif
        ) - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data + ( &
        real( &
#ifdef ENFORCE_SINGLE_PRECISION
          scale_val * this%source_data / div + offset_val, kind=real32 &
#else
          scale_val * this%source_data / div + offset_val, kind=real64 &
#endif
        ) - this%aggregated_data &
      ) / (this%counter + 1)
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data + ( &
        real( &
#ifdef ENFORCE_SINGLE_PRECISION
          scale_val * this%source_data / div + offset_val, kind=real32 &
#else
          scale_val * this%source_data / div + offset_val, kind=real64 &
#endif
        ) - this%aggregated_data &
      ) / (this%counter + 1)
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = this%counter + 1

  end subroutine mean_accumulate

  subroutine sum_accumulate(this, scale, div, offset)
    !* Accumulate the aggregated data from the source data using the sum
    ! aggregation method.
    class(aggregator_t), intent(inout) :: this
    real, intent(in), optional :: scale
      !* An optional scaling factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: div
      !* An optional division factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: offset
      !* An optional offset to add to the source data before accumulation.
      ! Defaults to 0.0 if not provided.
    real :: scale_val, div_val, offset_val

    call get_accumulate_args(scale, div, offset, scale_val, div_val, offset_val)

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = this%aggregated_data + int( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_int32_1d_t)
      this%aggregated_data = this%aggregated_data + int( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_int32_2d_t)
      this%aggregated_data = this%aggregated_data + int( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_int32_3d_t)
      this%aggregated_data = this%aggregated_data + int( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data + ( &
        scale_val * this%source_data / div_val + offset_val &
      )
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data + real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data + real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data + real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data + real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = this%counter + 1

  end subroutine sum_accumulate

  subroutine point_accumulate(this, scale, div, offset)
    !* Accumulate the aggregated data from the source data using the point
    ! aggregation method.
    class(aggregator_t), intent(inout) :: this
    real, intent(in), optional :: scale
      !* An optional scaling factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: div
      !* An optional division factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: offset
      !* An optional offset to add to the source data before accumulation.
      ! Defaults to 0.0 if not provided.
    real :: scale_val, div_val, offset_val

    call get_accumulate_args(scale, div, offset, scale_val, div_val, offset_val)

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = int(scale_val * this%source_data / div_val + offset_val)
    type is (aggregator_int32_1d_t)
      this%aggregated_data = int(scale_val * this%source_data / div_val + offset_val)
    type is (aggregator_int32_2d_t)
      this%aggregated_data = int(scale_val * this%source_data / div_val + offset_val)
    type is (aggregator_int32_3d_t)
      this%aggregated_data = int(scale_val * this%source_data / div_val + offset_val)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = scale_val * this%source_data / div_val + offset_val
    type is (aggregator_real32_1d_t)
      this%aggregated_data = scale_val * this%source_data / div_val + offset_val
    type is (aggregator_real32_2d_t)
      this%aggregated_data = scale_val * this%source_data / div_val + offset_val
    type is (aggregator_real32_3d_t)
      this%aggregated_data = scale_val * this%source_data / div_val + offset_val
    type is (aggregator_real64_0d_t)
      this%aggregated_data = real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    type is (aggregator_real64_1d_t)
      this%aggregated_data = real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    type is (aggregator_real64_2d_t)
      this%aggregated_data = real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    type is (aggregator_real64_3d_t)
      this%aggregated_data = real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      )
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = this%counter + 1

  end subroutine point_accumulate

  subroutine min_accumulate(this, scale, div, offset)
    !* Accumulate the aggregated data from the source data using the min
    ! aggregation method.
    class(aggregator_t), intent(inout) :: this
    real, intent(in), optional :: scale
      !* An optional scaling factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: div
      !* An optional division factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: offset
      !* An optional offset to add to the source data before accumulation.
      ! Defaults to 0.0 if not provided.
    real :: scale_val, div_val, offset_val

    call get_accumulate_args(scale, div, offset, scale_val, div_val, offset_val)

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = min(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_int32_1d_t)
      this%aggregated_data = min(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_int32_2d_t)
      this%aggregated_data = min(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_int32_3d_t)
      this%aggregated_data = min(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_0d_t)
      this%aggregated_data = min(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_1d_t)
      this%aggregated_data = min(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_2d_t)
      this%aggregated_data = min(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_3d_t)
      this%aggregated_data = min(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real64_0d_t)
      this%aggregated_data = min(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    type is (aggregator_real64_1d_t)
      this%aggregated_data = min(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    type is (aggregator_real64_2d_t)
      this%aggregated_data = min(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    type is (aggregator_real64_3d_t)
      this%aggregated_data = min(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = this%counter + 1

  end subroutine min_accumulate

  subroutine max_accumulate(this, scale, div, offset)
    !* Accumulate the aggregated data from the source data using the max
    ! aggregation method.
    class(aggregator_t), intent(inout) :: this
    real, intent(in), optional :: scale
      !* An optional scaling factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: div
      !* An optional division factor to apply to the source data before
      ! accumulation. Defaults to 1.0 if not provided.
    real, intent(in), optional :: offset
      !* An optional offset to add to the source data before accumulation.
      ! Defaults to 0.0 if not provided.
    real :: scale_val, div_val, offset_val

    call get_accumulate_args(scale, div, offset, scale_val, div_val, offset_val)

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = max(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_int32_1d_t)
      this%aggregated_data = max(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_int32_2d_t)
      this%aggregated_data = max(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_int32_3d_t)
      this%aggregated_data = max(this%aggregated_data, int( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_0d_t)
      this%aggregated_data = max(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_1d_t)
      this%aggregated_data = max(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_2d_t)
      this%aggregated_data = max(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real32_3d_t)
      this%aggregated_data = max(this%aggregated_data, ( &
        scale_val * this%source_data / div_val + offset_val &
      ))
    type is (aggregator_real64_0d_t)
      this%aggregated_data = max(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    type is (aggregator_real64_1d_t)
      this%aggregated_data = max(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    type is (aggregator_real64_2d_t)
      this%aggregated_data = max(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    type is (aggregator_real64_3d_t)
      this%aggregated_data = max(this%aggregated_data, real( &
#ifdef ENFORCE_SINGLE_PRECISION
        scale_val * this%source_data / div_val + offset_val, kind=real32 &
#else
        scale_val * this%source_data / div_val + offset_val, kind=real64 &
#endif
      ))
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = this%counter + 1

  end subroutine max_accumulate

  subroutine point_reset(this)
    !* Reset the aggregated data for the point aggregation method. This is a
    ! no-op since point aggregation always takes the value of the most recent data
    ! point.
    class(aggregator_t), intent(inout) :: this
  end subroutine point_reset

  subroutine min_reset(this)
    !! Reset the aggregated data for the min aggregation method.
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = huge(0_int32)
    type is (aggregator_int32_1d_t)
      this%aggregated_data = huge(0_int32)
    type is (aggregator_int32_2d_t)
      this%aggregated_data = huge(0_int32)
    type is (aggregator_int32_3d_t)
      this%aggregated_data = huge(0_int32)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = huge(0.0_real32)
    type is (aggregator_real32_1d_t)
      this%aggregated_data = huge(0.0_real32)
    type is (aggregator_real32_2d_t)
      this%aggregated_data = huge(0.0_real32)
    type is (aggregator_real32_3d_t)
      this%aggregated_data = huge(0.0_real32)
    type is (aggregator_real64_0d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = huge(0.0_real32)
#else
      this%aggregated_data = huge(0.0_real64)
#endif
    type is (aggregator_real64_1d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = huge(0.0_real32)
#else
      this%aggregated_data = huge(0.0_real64)
#endif
    type is (aggregator_real64_2d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = huge(0.0_real32)
#else
      this%aggregated_data = huge(0.0_real64)
#endif
    type is (aggregator_real64_3d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = huge(0.0_real32)
#else
      this%aggregated_data = huge(0.0_real64)
#endif
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = 0

  end subroutine min_reset

  subroutine max_reset(this)
    !! Reset the aggregated data for the max aggregation method.
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = -huge(0_int32)
    type is (aggregator_int32_1d_t)
      this%aggregated_data = -huge(0_int32)
    type is (aggregator_int32_2d_t)
      this%aggregated_data = -huge(0_int32)
    type is (aggregator_int32_3d_t)
      this%aggregated_data = -huge(0_int32)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = -huge(0.0_real32)
    type is (aggregator_real32_1d_t)
      this%aggregated_data = -huge(0.0_real32)
    type is (aggregator_real32_2d_t)
      this%aggregated_data = -huge(0.0_real32)
    type is (aggregator_real32_3d_t)
      this%aggregated_data = -huge(0.0_real32)
    type is (aggregator_real64_0d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = -huge(0.0_real32)
#else
      this%aggregated_data = -huge(0.0_real64)
#endif
    type is (aggregator_real64_1d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = -huge(0.0_real32)
#else
      this%aggregated_data = -huge(0.0_real64)
#endif
    type is (aggregator_real64_2d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = -huge(0.0_real32)
#else
      this%aggregated_data = -huge(0.0_real64)
#endif
    type is (aggregator_real64_3d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = -huge(0.0_real32)
#else
      this%aggregated_data = -huge(0.0_real64)
#endif
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = 0

  end subroutine max_reset

  subroutine other_reset(this)
    !! Reset the aggregated data for aggregation methods other than point, min, and max.
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = 0_int32
    type is (aggregator_int32_1d_t)
      this%aggregated_data = 0_int32
    type is (aggregator_int32_2d_t)
      this%aggregated_data = 0_int32
    type is (aggregator_int32_3d_t)
      this%aggregated_data = 0_int32
    type is (aggregator_real32_0d_t)
      this%aggregated_data = 0.0_real32
    type is (aggregator_real32_1d_t)
      this%aggregated_data = 0.0_real32
    type is (aggregator_real32_2d_t)
      this%aggregated_data = 0.0_real32
    type is (aggregator_real32_3d_t)
      this%aggregated_data = 0.0_real32
    type is (aggregator_real64_0d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = 0.0_real32
#else
      this%aggregated_data = 0.0_real64
#endif
    type is (aggregator_real64_1d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = 0.0_real32
#else
      this%aggregated_data = 0.0_real64
#endif
    type is (aggregator_real64_2d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = 0.0_real32
#else
      this%aggregated_data = 0.0_real64
#endif
    type is (aggregator_real64_3d_t)
#ifdef ENFORCE_SINGLE_PRECISION
      this%aggregated_data = 0.0_real32
#else
      this%aggregated_data = 0.0_real64
#endif
    class default
      call cable_abort("Unexpected aggregator type.", file=__FILE__, line=__LINE__)
    end select

    this%counter = 0

  end subroutine other_reset

  function new_aggregator_int32_0d_t(source_data) result(agg)
    !! Create a new 0D integer aggregator.
    integer(kind=int32), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_int32_0d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_0d_t

  function new_aggregator_int32_1d_t(source_data) result(agg)
    !! Create a new 1D integer aggregator.
    integer(kind=int32), dimension(:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_int32_1d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_1d_t

  function new_aggregator_int32_2d_t(source_data) result(agg)
    !! Create a new 2D integer aggregator.
    integer(kind=int32), dimension(:,:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_int32_2d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_2d_t

  function new_aggregator_int32_3d_t(source_data) result(agg)
    !! Create a new 3D integer aggregator.
    integer(kind=int32), dimension(:,:,:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_int32_3d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_3d_t

  function new_aggregator_real32_0d(source_data) result(agg)
    !! Create a new 0D 32-bit real aggregator.
    real(kind=real32), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real32_0d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_0d

  function new_aggregator_real32_1d(source_data) result(agg)
    !! Create a new 1D 32-bit real aggregator.
    real(kind=real32), dimension(:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real32_1d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_1d

  function new_aggregator_real32_2d(source_data) result(agg)
    !! Create a new 2D 32-bit real aggregator.
    real(kind=real32), dimension(:,:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real32_2d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_2d

  function new_aggregator_real32_3d(source_data) result(agg)
    !! Create a new 3D 32-bit real aggregator.
    real(kind=real32), dimension(:,:,:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real32_3d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_3d

  function new_aggregator_real64_0d(source_data) result(agg)
    !! Create a new 0D 64-bit real aggregator.
    real(kind=real64), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real64_0d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_0d

  function new_aggregator_real64_1d(source_data) result(agg)
    !! Create a new 1D 64-bit real aggregator.
    real(kind=real64), dimension(:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real64_1d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_1d

  function new_aggregator_real64_2d(source_data) result(agg)
    !! Create a new 2D 64-bit real aggregator.
    real(kind=real64), dimension(:,:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real64_2d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_2d

  function new_aggregator_real64_3d(source_data) result(agg)
    !! Create a new 3D 64-bit real aggregator.
    real(kind=real64), dimension(:,:,:), intent(inout), target :: source_data
      !! The source data array to be sampled by the aggregator.
    type(aggregator_real64_3d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_3d

end module
