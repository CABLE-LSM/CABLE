! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_yaml_mod
  !! Minimal YAML reading interface for CABLE.
  !!
  !! The YAML file is parsed with yaFyaml and immediately copied into a plain
  !! tree of [[cable_yaml_node_t]] values, so no yaFyaml types are visible
  !! outside this module. Every node is a scalar, a mapping or a sequence.
  !! Mapping entries keep their file order, and each child node stores the key
  !! it was found under.
  use, intrinsic :: iso_fortran_env, only: int64
  use yafyaml, only: Parser, TextStream, YAML_Node, NodeIterator
  use yafyaml, only: to_string, to_int, to_float, to_bool
  use yafyaml, only: YAFYAML_SUCCESS
  use cable_error_handler_mod, only: cable_abort
  implicit none
  private

  public :: cable_yaml_load_file
  public :: cable_yaml_parse_text

  integer, parameter, public :: CABLE_YAML_UNDEFINED = 0
  integer, parameter, public :: CABLE_YAML_SCALAR = 1
  integer, parameter, public :: CABLE_YAML_MAPPING = 2
  integer, parameter, public :: CABLE_YAML_SEQUENCE = 3

  type, public :: cable_yaml_node_t
    !! A node of a parsed YAML document.
    integer :: node_kind = CABLE_YAML_UNDEFINED
      !! One of the `CABLE_YAML_*` kind constants.
    character(:), allocatable :: key
      !! The mapping key this node was found under. Unallocated for the root
      !! node and for sequence items.
    character(:), allocatable :: text
      !! The scalar value as text. Only allocated for scalar nodes.
    type(cable_yaml_node_t), allocatable :: children(:)
      !! Mapping values (in file order) or sequence items.
  contains
    procedure :: copy_node
    generic :: assignment(=) => copy_node
    procedure :: is_scalar => node_is_scalar
    procedure :: is_mapping => node_is_mapping
    procedure :: is_sequence => node_is_sequence
    procedure :: size => node_size
    procedure :: has => node_has
    procedure :: get => node_get
    procedure :: at => node_at
    procedure :: as_string => node_as_string
    procedure :: as_integer => node_as_integer
    procedure :: as_real => node_as_real
    procedure :: as_logical => node_as_logical
    procedure :: get_string => node_get_string
    procedure :: get_integer => node_get_integer
    procedure :: get_real => node_get_real
    procedure :: get_logical => node_get_logical
  end type cable_yaml_node_t

contains

  recursive subroutine copy_node(destination, source)
    !! Deep copy of a node and everything below it.
    !!
    !! This replaces intrinsic assignment, which gfortran 13 implements
    !! incorrectly for recursive allocatable components (the copies share
    !! storage and are freed twice). The source is staged into local variables
    !! first so that the destination may be part of the source, as in
    !! `node = node%at(1)`.
    class(cable_yaml_node_t), intent(inout) :: destination
    type(cable_yaml_node_t), intent(in) :: source
    type(cable_yaml_node_t), allocatable :: staged_children(:)
    character(:), allocatable :: staged_key, staged_text
    integer :: staged_kind
    integer :: child_index

    ! Step 1: copy everything out of the source into local variables. If the
    ! destination is part of the source (`node = node%at(1)`), clearing the
    ! destination in step 2 would otherwise destroy the data we still need.
    staged_kind = source%node_kind
    if (allocated(source%key)) staged_key = source%key
    if (allocated(source%text)) staged_text = source%text
    if (allocated(source%children)) then
      allocate(staged_children(size(source%children)))
      do child_index = 1, size(source%children)
        ! This assignment calls copy_node again, so the whole subtree is copied.
        staged_children(child_index) = source%children(child_index)
      end do
    end if

    ! Step 2: replace the destination's contents with the staged copies. The
    ! move_alloc calls transfer ownership of the copies without another copy.
    destination%node_kind = staged_kind
    if (allocated(destination%key)) deallocate(destination%key)
    if (allocated(staged_key)) call move_alloc(staged_key, destination%key)
    if (allocated(destination%text)) deallocate(destination%text)
    if (allocated(staged_text)) call move_alloc(staged_text, destination%text)
    if (allocated(destination%children)) deallocate(destination%children)
    if (allocated(staged_children)) call move_alloc(staged_children, destination%children)
  end subroutine copy_node

  function cable_yaml_load_file(file_name) result(root)
    !! Parse the YAML file at `file_name`. Aborts if the file cannot be read or
    !! is not valid YAML.
    character(len=*), intent(in) :: file_name
    type(cable_yaml_node_t) :: root
    type(Parser) :: yaml_parser
    class(YAML_Node), allocatable :: parsed_root
    logical :: file_exists
    integer :: status

    ! Check the file exists first, so a missing file gives a clear message
    ! rather than a generic parser failure.
    inquire(file=file_name, exist=file_exists)
    if (.not. file_exists) then
      call cable_abort("YAML file does not exist: " // trim(file_name), __FILE__, __LINE__)
    end if
    ! The 'core' schema reads plain scalars as strings, integers, reals or
    ! booleans according to how they are written.
    yaml_parser = Parser('core')
    parsed_root = yaml_parser%load(file_name, rc=status)
    if (status /= YAFYAML_SUCCESS) then
      call cable_abort("Failed to parse YAML file: " // trim(file_name), __FILE__, __LINE__)
    end if
    ! Copy the parsed tree into our own node type so that no yaFyaml types
    ! are visible to the caller (see the module description).
    root = convert_node(parsed_root)
  end function cable_yaml_load_file

  function cable_yaml_parse_text(text) result(root)
    !! Parse YAML held in a string. Aborts if the text is not valid YAML.
    character(len=*), intent(in) :: text
    type(cable_yaml_node_t) :: root
    type(Parser) :: yaml_parser
    class(YAML_Node), allocatable :: parsed_root
    integer :: status

    ! Same as loading a file, but reading from a string held in memory.
    yaml_parser = Parser('core')
    parsed_root = yaml_parser%load(TextStream(text), rc=status)
    if (status /= YAFYAML_SUCCESS) then
      call cable_abort("Failed to parse YAML text", __FILE__, __LINE__)
    end if
    root = convert_node(parsed_root)
  end function cable_yaml_parse_text

  recursive function convert_node(parsed_node) result(converted)
    !! Copy a yaFyaml node and everything below it into a `cable_yaml_node_t`.
    class(YAML_Node), intent(in), target :: parsed_node
    type(cable_yaml_node_t) :: converted
    class(YAML_Node), pointer :: child_node
    class(NodeIterator), allocatable :: iterator
    integer :: child_index, child_count

    ! A node is one of three kinds, and each is copied differently: a mapping
    ! (key: value pairs), a sequence (a list) or a scalar (a single value).
    if (parsed_node%is_mapping()) then
      converted%node_kind = CABLE_YAML_MAPPING
      ! Count the entries first so the children can be allocated once and
      ! filled in place; growing the array by concatenation double-frees
      ! recursive allocatable components in gfortran.
      associate (first => parsed_node%begin(), last => parsed_node%end())
        ! First pass: walk the mapping once just to count its entries.
        iterator = first
        child_count = 0
        do while (iterator /= last)
          child_count = child_count + 1
          call iterator%next()
        end do
        allocate(converted%children(child_count))
        ! Second pass: walk it again, converting each value and recording the
        ! key it was stored under. Entries keep the order of the file.
        iterator = first
        child_index = 0
        do while (iterator /= last)
          child_index = child_index + 1
          child_node => iterator%second()
          converted%children(child_index) = convert_node(child_node)
          ! Keys can be numbers (stream `1:`), so they go through scalar_text
          ! to be stored as text whatever their type in the file.
          converted%children(child_index)%key = scalar_text(iterator%first())
          call iterator%next()
        end do
      end associate
    else if (parsed_node%is_sequence()) then
      ! A sequence knows its own length, so it can be filled in a single pass.
      converted%node_kind = CABLE_YAML_SEQUENCE
      allocate(converted%children(parsed_node%size()))
      do child_index = 1, parsed_node%size()
        child_node => parsed_node%at(child_index)
        converted%children(child_index) = convert_node(child_node)
      end do
    else
      ! Anything else is a single value, kept as its text.
      converted%node_kind = CABLE_YAML_SCALAR
      converted%text = scalar_text(parsed_node)
    end if
  end function convert_node

  function scalar_text(parsed_node) result(text)
    !! The text of a scalar yaFyaml node. Non-string scalars are rendered in the
    !! form they are written in YAML (`true`/`false` for booleans).
    class(YAML_Node), intent(in), target :: parsed_node
    character(:), allocatable :: text
    character(len=64) :: number_buffer

    ! yaFyaml has already decided what type the scalar is. Every type is
    ! turned back into text, so later code converts it once, where it knows
    ! what type it expects (and can report a mistake by name).
    if (parsed_node%is_string()) then
      text = to_string(parsed_node)
    else if (parsed_node%is_int()) then
      write(number_buffer, '(i0)') to_int(parsed_node)
      text = trim(number_buffer)
    else if (parsed_node%is_bool()) then
      ! merge needs both choices to be the same length, hence the padded 'true '
      ! and the trim afterwards.
      text = merge('true ', 'false', to_bool(parsed_node))
      text = trim(text)
    else if (parsed_node%is_float()) then
      write(number_buffer, '(g0)') to_float(parsed_node)
      text = trim(number_buffer)
    else
      call cable_abort("Unsupported YAML scalar type", __FILE__, __LINE__)
    end if
  end function scalar_text

  pure logical function node_is_scalar(this)
    class(cable_yaml_node_t), intent(in) :: this
    node_is_scalar = this%node_kind == CABLE_YAML_SCALAR
  end function node_is_scalar

  pure logical function node_is_mapping(this)
    class(cable_yaml_node_t), intent(in) :: this
    node_is_mapping = this%node_kind == CABLE_YAML_MAPPING
  end function node_is_mapping

  pure logical function node_is_sequence(this)
    class(cable_yaml_node_t), intent(in) :: this
    node_is_sequence = this%node_kind == CABLE_YAML_SEQUENCE
  end function node_is_sequence

  pure integer function node_size(this)
    !! Number of children of a mapping or sequence. Zero for scalars.
    class(cable_yaml_node_t), intent(in) :: this
    node_size = 0
    if (allocated(this%children)) node_size = size(this%children)
  end function node_size

  pure function node_index(this, key) result(child_index)
    !! Position of `key` among the children of a mapping, or 0 if absent.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    integer :: child_index
    integer :: candidate

    child_index = 0
    ! Only a mapping has keys, so anything else simply has no match.
    if (.not. this%is_mapping()) return
    ! A linear search is fine: configuration entries have a handful of keys.
    ! Fortran ignores trailing blanks when comparing strings, so a stored key
    ! matches whether or not it is padded.
    do candidate = 1, size(this%children)
      if (this%children(candidate)%key == key) then
        child_index = candidate
        return
      end if
    end do
  end function node_index

  pure logical function node_has(this, key)
    !! Whether this node is a mapping containing `key`.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    node_has = node_index(this, key) > 0
  end function node_has

  function node_get(this, key) result(child)
    !! The value stored under `key`. Aborts if this node has no such key.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    type(cable_yaml_node_t) :: child
    integer :: child_index

    child_index = node_index(this, key)
    ! Callers that cannot be sure the key is there check with `has` first, so
    ! reaching this point with a missing key is a mistake worth stopping for.
    if (child_index == 0) then
      call cable_abort("YAML key not found: " // trim(key), __FILE__, __LINE__)
    end if
    child = this%children(child_index)
  end function node_get

  function node_at(this, child_index) result(child)
    !! The `child_index`-th child (1-based) of a mapping or sequence.
    class(cable_yaml_node_t), intent(in) :: this
    integer, intent(in) :: child_index
    type(cable_yaml_node_t) :: child

    ! Explicit range check so a bad index gives a message, not a crash
    ! from reading outside the array.
    if (child_index < 1 .or. child_index > this%size()) then
      call cable_abort("YAML child index out of range", __FILE__, __LINE__)
    end if
    child = this%children(child_index)
  end function node_at

  function node_as_string(this) result(text)
    !! The text of a scalar node. Aborts for mappings and sequences.
    class(cable_yaml_node_t), intent(in) :: this
    character(:), allocatable :: text

    if (.not. this%is_scalar()) then
      call cable_abort("YAML node is not a scalar", __FILE__, __LINE__)
    end if
    text = this%text
  end function node_as_string

  function node_as_integer(this) result(value)
    !! The scalar value as an integer. Aborts if it is not one.
    class(cable_yaml_node_t), intent(in) :: this
    integer :: value
    character(:), allocatable :: text
    integer :: read_status

    ! List-directed read into an integer; a non-zero status means the text is
    ! not a whole number. `text` is a variable because Fortran cannot read
    ! from a function result.
    text = this%as_string()
    read(text, *, iostat=read_status) value
    if (read_status /= 0) then
      call cable_abort("YAML value is not an integer: " // text, __FILE__, __LINE__)
    end if
  end function node_as_integer

  function node_as_real(this) result(value)
    !! The scalar value as a real. Aborts if it is not one.
    class(cable_yaml_node_t), intent(in) :: this
    real :: value
    character(:), allocatable :: text
    integer :: read_status

    ! Same approach as for integers. Whole numbers are accepted as reals.
    text = this%as_string()
    read(text, *, iostat=read_status) value
    if (read_status /= 0) then
      call cable_abort("YAML value is not a real number: " // text, __FILE__, __LINE__)
    end if
  end function node_as_real

  function node_as_logical(this) result(value)
    !! The scalar value as a logical (`true` or `false`). Aborts otherwise.
    class(cable_yaml_node_t), intent(in) :: this
    logical :: value
    character(:), allocatable :: text

    ! YAML allows several spellings of true and false. These three each are
    ! accepted; anything else (such as yes or 1) is rejected rather than guessed.
    text = this%as_string()
    select case (text)
    case ("true", "True", "TRUE")
      value = .true.
    case ("false", "False", "FALSE")
      value = .false.
    case default
      call cable_abort("YAML value is not a boolean: " // text, __FILE__, __LINE__)
    end select
  end function node_as_logical

  function node_get_string(this, key) result(text)
    !! Shorthand for `this%get(key)%as_string()`.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    character(:), allocatable :: text
    type(cable_yaml_node_t) :: child

    child = this%get(key)
    text = child%as_string()
  end function node_get_string

  function node_get_integer(this, key) result(value)
    !! Shorthand for `this%get(key)%as_integer()`.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    integer :: value
    type(cable_yaml_node_t) :: child

    child = this%get(key)
    value = child%as_integer()
  end function node_get_integer

  function node_get_real(this, key) result(value)
    !! Shorthand for `this%get(key)%as_real()`.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    real :: value
    type(cable_yaml_node_t) :: child

    child = this%get(key)
    value = child%as_real()
  end function node_get_real

  function node_get_logical(this, key) result(value)
    !! Shorthand for `this%get(key)%as_logical()`.
    class(cable_yaml_node_t), intent(in) :: this
    character(len=*), intent(in) :: key
    logical :: value
    type(cable_yaml_node_t) :: child

    child = this%get(key)
    value = child%as_logical()
  end function node_get_logical

end module cable_yaml_mod
