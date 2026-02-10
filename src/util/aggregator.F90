module aggregator_mod

  use iso_fortran_env, only: int32, real32, real64
  use cable_abort_module, only: cable_abort

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
    integer :: counter = 0
    procedure(accumulate_data), pointer :: accumulate
    procedure(reset_data), pointer :: reset
  contains
    procedure :: init => aggregator_init
    procedure :: set_method => aggregator_set_method
    procedure :: rank => aggregator_rank
    procedure :: shape => aggregator_shape
    procedure :: scale => aggregator_scale
    procedure :: offset => aggregator_offset
  end type aggregator_t

  abstract interface
    subroutine accumulate_data(this)
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
    end subroutine accumulate_data
    subroutine reset_data(this)
      import aggregator_t
      class(aggregator_t), intent(inout) :: this
    end subroutine reset_data
  end interface

  type, extends(aggregator_t) :: aggregator_int32_0d_t
    integer(kind=int32), allocatable :: aggregated_data
    integer(kind=int32), pointer :: source_data => null()
  end type aggregator_int32_0d_t

  type, extends(aggregator_t) :: aggregator_int32_1d_t
    integer(kind=int32), dimension(:), allocatable :: aggregated_data
    integer(kind=int32), dimension(:), pointer :: source_data => null()
  end type aggregator_int32_1d_t

  type, extends(aggregator_t) :: aggregator_int32_2d_t
    integer(kind=int32), dimension(:,:), allocatable :: aggregated_data
    integer(kind=int32), dimension(:,:), pointer :: source_data => null()
  end type aggregator_int32_2d_t

  type, extends(aggregator_t) :: aggregator_int32_3d_t
    integer(kind=int32), dimension(:,:,:), allocatable :: aggregated_data
    integer(kind=int32), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_int32_3d_t

  type, extends(aggregator_t) :: aggregator_real32_0d_t
    real(kind=real32), allocatable :: aggregated_data
    real(kind=real32), pointer :: source_data => null()
  end type aggregator_real32_0d_t

  type, extends(aggregator_t) :: aggregator_real32_1d_t
    real(kind=real32), dimension(:), allocatable :: aggregated_data
    real(kind=real32), dimension(:), pointer :: source_data => null()
  end type aggregator_real32_1d_t

  type, extends(aggregator_t) :: aggregator_real32_2d_t
    real(kind=real32), dimension(:,:), allocatable :: aggregated_data
    real(kind=real32), dimension(:,:), pointer :: source_data => null()
  end type aggregator_real32_2d_t

  type, extends(aggregator_t) :: aggregator_real32_3d_t
    real(kind=real32), dimension(:,:,:), allocatable :: aggregated_data
    real(kind=real32), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_real32_3d_t

  type, extends(aggregator_t) :: aggregator_real64_0d_t
    real(kind=real64), allocatable :: aggregated_data
    real(kind=real64), pointer :: source_data => null()
  end type aggregator_real64_0d_t

  type, extends(aggregator_t) :: aggregator_real64_1d_t
    real(kind=real64), dimension(:), allocatable :: aggregated_data
    real(kind=real64), dimension(:), pointer :: source_data => null()
  end type aggregator_real64_1d_t

  type, extends(aggregator_t) :: aggregator_real64_2d_t
    real(kind=real64), dimension(:,:), allocatable :: aggregated_data
    real(kind=real64), dimension(:,:), pointer :: source_data => null()
  end type aggregator_real64_2d_t

  type, extends(aggregator_t) :: aggregator_real64_3d_t
    real(kind=real64), dimension(:,:,:), allocatable :: aggregated_data
    real(kind=real64), dimension(:,:,:), pointer :: source_data => null()
  end type aggregator_real64_3d_t

  interface new_aggregator
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
    class(aggregator_t), intent(inout) :: this
    character(len=*), intent(in) :: method

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
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real64_1d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real64_2d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    type is (aggregator_real64_3d_t)
      if (.not. allocated(this%aggregated_data)) allocate(this%aggregated_data, mold=this%source_data)
    end select

    call this%set_method(method)

    call this%reset()

  end subroutine aggregator_init

  subroutine aggregator_set_method(this, method)
    class(aggregator_t), intent(inout) :: this
    character(len=*), intent(in) :: method

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
      call cable_abort("Aggregation method "//method//" is invalid.")
    endif

  end subroutine aggregator_set_method

  integer function aggregator_rank(this)
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
    end select

  end function aggregator_rank

  function aggregator_shape(this) result(agg_shape)
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
    end select

  end function aggregator_shape

  subroutine aggregator_scale(this, scale)
    class(aggregator_t), intent(inout) :: this
    real, intent(in) :: scale

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
    end select

  end subroutine aggregator_scale

  subroutine aggregator_offset(this, offset)
    class(aggregator_t), intent(inout) :: this
    real, intent(in) :: offset

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
    end select

  end subroutine aggregator_offset

  subroutine mean_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data + (this%source_data - this%aggregated_data) / (this%counter + 1)
    end select

    this%counter = this%counter + 1

  end subroutine mean_accumulate

  subroutine sum_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_int32_1d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_int32_2d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_int32_3d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%aggregated_data + this%source_data
    end select

    this%counter = this%counter + 1

  end subroutine sum_accumulate

  subroutine point_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_int32_1d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_int32_2d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_int32_3d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real32_0d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real32_1d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real32_2d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real32_3d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real64_0d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real64_1d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real64_2d_t)
      this%aggregated_data = this%source_data
    type is (aggregator_real64_3d_t)
      this%aggregated_data = this%source_data
    end select

    this%counter = this%counter + 1

  end subroutine point_accumulate

  subroutine min_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_int32_1d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_int32_2d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_int32_3d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real32_1d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real32_2d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real32_3d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real64_0d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real64_1d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real64_2d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    type is (aggregator_real64_3d_t)
      this%aggregated_data = min(this%aggregated_data, this%source_data)
    end select

    this%counter = this%counter + 1

  end subroutine min_accumulate

  subroutine max_accumulate(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_int32_1d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_int32_2d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_int32_3d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real32_0d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real32_1d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real32_2d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real32_3d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real64_0d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real64_1d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real64_2d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    type is (aggregator_real64_3d_t)
      this%aggregated_data = max(this%aggregated_data, this%source_data)
    end select

    this%counter = this%counter + 1

  end subroutine max_accumulate

  subroutine point_reset(this)
    class(aggregator_t), intent(inout) :: this
  end subroutine point_reset

  subroutine min_reset(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = huge(int(0_int32))
    type is (aggregator_int32_1d_t)
      this%aggregated_data = huge(int(0_int32))
    type is (aggregator_int32_2d_t)
      this%aggregated_data = huge(int(0_int32))
    type is (aggregator_int32_3d_t)
      this%aggregated_data = huge(int(0_int32))
    type is (aggregator_real32_0d_t)
      this%aggregated_data = huge(real(0.0_real32))
    type is (aggregator_real32_1d_t)
      this%aggregated_data = huge(real(0.0_real32))
    type is (aggregator_real32_2d_t)
      this%aggregated_data = huge(real(0.0_real32))
    type is (aggregator_real32_3d_t)
      this%aggregated_data = huge(real(0.0_real32))
    type is (aggregator_real64_0d_t)
      this%aggregated_data = huge(real(0.0_real64))
    type is (aggregator_real64_1d_t)
      this%aggregated_data = huge(real(0.0_real64))
    type is (aggregator_real64_2d_t)
      this%aggregated_data = huge(real(0.0_real64))
    type is (aggregator_real64_3d_t)
      this%aggregated_data = huge(real(0.0_real64))
    end select

    this%counter = 0

  end subroutine min_reset

  subroutine max_reset(this)
    class(aggregator_t), intent(inout) :: this

    select type (this)
    type is (aggregator_int32_0d_t)
      this%aggregated_data = -huge(int(0_int32))
    type is (aggregator_int32_1d_t)
      this%aggregated_data = -huge(int(0_int32))
    type is (aggregator_int32_2d_t)
      this%aggregated_data = -huge(int(0_int32))
    type is (aggregator_int32_3d_t)
      this%aggregated_data = -huge(int(0_int32))
    type is (aggregator_real32_0d_t)
      this%aggregated_data = -huge(real(0.0_real32))
    type is (aggregator_real32_1d_t)
      this%aggregated_data = -huge(real(0.0_real32))
    type is (aggregator_real32_2d_t)
      this%aggregated_data = -huge(real(0.0_real32))
    type is (aggregator_real32_3d_t)
      this%aggregated_data = -huge(real(0.0_real32))
    type is (aggregator_real64_0d_t)
      this%aggregated_data = -huge(real(0.0_real64))
    type is (aggregator_real64_1d_t)
      this%aggregated_data = -huge(real(0.0_real64))
    type is (aggregator_real64_2d_t)
      this%aggregated_data = -huge(real(0.0_real64))
    type is (aggregator_real64_3d_t)
      this%aggregated_data = -huge(real(0.0_real64))
    end select

    this%counter = 0

  end subroutine max_reset

  subroutine other_reset(this)
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
      this%aggregated_data = 0.0_real64
    type is (aggregator_real64_1d_t)
      this%aggregated_data = 0.0_real64
    type is (aggregator_real64_2d_t)
      this%aggregated_data = 0.0_real64
    type is (aggregator_real64_3d_t)
      this%aggregated_data = 0.0_real64
    end select

    this%counter = 0

  end subroutine other_reset

  function new_aggregator_int32_0d_t(source_data) result(agg)
    integer(kind=int32), intent(inout), target :: source_data
    type(aggregator_int32_0d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_0d_t

  function new_aggregator_int32_1d_t(source_data) result(agg)
    integer(kind=int32), dimension(:), intent(inout), target :: source_data
    type(aggregator_int32_1d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_1d_t

  function new_aggregator_int32_2d_t(source_data) result(agg)
    integer(kind=int32), dimension(:,:), intent(inout), target :: source_data
    type(aggregator_int32_2d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_2d_t

  function new_aggregator_int32_3d_t(source_data) result(agg)
    integer(kind=int32), dimension(:,:,:), intent(inout), target :: source_data
    type(aggregator_int32_3d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_int32_3d_t

  function new_aggregator_real32_0d(source_data) result(agg)
    real(kind=real32), intent(inout), target :: source_data
    type(aggregator_real32_0d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_0d

  function new_aggregator_real32_1d(source_data) result(agg)
    real(kind=real32), dimension(:), intent(inout), target :: source_data
    type(aggregator_real32_1d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_1d

  function new_aggregator_real32_2d(source_data) result(agg)
    real(kind=real32), dimension(:,:), intent(inout), target :: source_data
    type(aggregator_real32_2d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_2d

  function new_aggregator_real32_3d(source_data) result(agg)
    real(kind=real32), dimension(:,:,:), intent(inout), target :: source_data
    type(aggregator_real32_3d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real32_3d

  function new_aggregator_real64_0d(source_data) result(agg)
    real(kind=real64), intent(inout), target :: source_data
    type(aggregator_real64_0d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_0d

  function new_aggregator_real64_1d(source_data) result(agg)
    real(kind=real64), dimension(:), intent(inout), target :: source_data
    type(aggregator_real64_1d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_1d

  function new_aggregator_real64_2d(source_data) result(agg)
    real(kind=real64), dimension(:,:), intent(inout), target :: source_data
    type(aggregator_real64_2d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_2d

  function new_aggregator_real64_3d(source_data) result(agg)
    real(kind=real64), dimension(:,:,:), intent(inout), target :: source_data
    type(aggregator_real64_3d_t) :: agg

    agg%source_data => source_data

  end function new_aggregator_real64_3d

end module
