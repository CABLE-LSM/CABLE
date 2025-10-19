module cable_netcdf_decomp_util_mod
  use cable_netcdf_mod, only: cable_netcdf_decomp_t, cable_netcdf_create_decomp, CABLE_NETCDF_MAX_STR_LEN_DIM
  use cable_array_utils_mod, only: array_index, array_offset
  use cable_abort_module, only: cable_abort
  implicit none

  private

  public &
    io_decomp_land_to_x_y, &
    io_decomp_patch_to_x_y_patch, &
    io_decomp_land_to_land, &
    io_decomp_patch_to_land_patch, &
    io_decomp_patch_to_patch, &
    dim_spec_t

  type dim_spec_t
    character(CABLE_NETCDF_MAX_STR_LEN_DIM) :: name
    integer :: size
  end type

contains

  integer function subscript(shape_spec, name)
    type(dim_spec_t), intent(in) :: shape_spec(:)
    character(*), intent(in) :: name
    integer i
    do i = 1, size(shape_spec)
      if (shape_spec(i)%name == name) then
        subscript = i
        return
      end if
    end do
    call cable_abort("Name '" // name // "' not found in shape_spec.", file=__FILE__, line=__LINE__)
    subscript = -1
  end function

  integer function patch_land_index(cstart, nap, patch_index)
    integer, intent(in) :: cstart(:), nap(:), patch_index
    integer i
    do i = 1, size(cstart)
      if (patch_index >= cstart(i) .and. patch_index <= cstart(i) + nap(i) - 1) then
        patch_land_index = i
        return
      end if
    end do
    call cable_abort("Patch index does not lie on any land point.", file=__FILE__, line=__LINE__)
    patch_land_index = -1
  end function

  function io_decomp_land_to_x_y(land_x, land_y, mem_shape_spec, var_shape_spec, type) result(decomp)
    integer, intent(in) :: land_x(:), land_y(:)
    type(dim_spec_t), intent(in) :: mem_shape_spec(:), var_shape_spec(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp
    
    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: mem_land_subscript
    integer :: i, mem_offset

    mem_land_subscript = subscript(mem_shape_spec, 'land')

    mem_shape = [(mem_shape_spec(i)%size, i = 1, size(mem_shape_spec))]
    grid_shape = [(var_shape_spec(i)%size, i = 1, size(var_shape_spec))]

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape_spec)
        select case (var_shape_spec(i)%name)
        case ('x')
          grid_index(i) = land_x(mem_index(mem_land_subscript))
        case ('y')
          grid_index(i) = land_y(mem_index(mem_land_subscript))
        case default
          grid_index(i) = mem_index(subscript(mem_shape_spec, var_shape_spec(i)%name))
        end select
      end do
      compmap(mem_offset) = array_offset(grid_index, grid_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, grid_shape, type)

  end function

  function io_decomp_patch_to_x_y_patch(land_x, land_y, cstart, nap, mem_shape_spec, var_shape_spec, type) result(decomp)
    integer, intent(in) :: land_x(:), land_y(:)
    integer, intent(in) :: cstart(:), nap(:) !! These are required to (a) get the land index of each patch index, and (b) get the patch offset value relative to cstart
    type(dim_spec_t), intent(in) :: mem_shape_spec(:), var_shape_spec(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: mem_patch_subscript
    integer :: i, mem_offset, land_index

    mem_patch_subscript = subscript(mem_shape_spec, 'patch')

    mem_shape = [(mem_shape_spec(i)%size, i = 1, size(mem_shape_spec))]
    grid_shape = [(var_shape_spec(i)%size, i = 1, size(var_shape_spec))]

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      land_index = patch_land_index(cstart, nap, mem_index(mem_patch_subscript))
      do i = 1, size(var_shape_spec)
        select case (var_shape_spec(i)%name)
        case ('x')
          grid_index(i) = land_x(land_index)
        case ('y')
          grid_index(i) = land_y(land_index)
        case ('patch')
          grid_index(i) = mem_index(mem_patch_subscript) - cstart(land_index) + 1
        case default
          grid_index(i) = mem_index(subscript(mem_shape_spec, var_shape_spec(i)%name))
        end select
      end do
      compmap(mem_offset) = array_offset(grid_index, grid_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, grid_shape, type)

  end function

  function io_decomp_land_to_land(land_decomp_start, mem_shape_spec, var_shape_spec, type) result(decomp)
    integer, intent(in) :: land_decomp_start
    type(dim_spec_t), intent(in) :: mem_shape_spec(:), var_shape_spec(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: mem_land_subscript
    integer :: i, mem_offset

    mem_land_subscript = subscript(mem_shape_spec, 'land')

    mem_shape = [(mem_shape_spec(i)%size, i = 1, size(mem_shape_spec))]
    grid_shape = [(var_shape_spec(i)%size, i = 1, size(var_shape_spec))]

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape_spec)
        select case (var_shape_spec(i)%name)
        case ('land')
          grid_index(i) = land_decomp_start + mem_index(mem_land_subscript) - 1
        case default
          grid_index(i) = mem_index(subscript(mem_shape_spec, var_shape_spec(i)%name))
        end select
      end do
      compmap(mem_offset) = array_offset(grid_index, grid_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, grid_shape, type)

  end function

  function io_decomp_patch_to_land_patch(land_decomp_start, cstart, nap, mem_shape_spec, var_shape_spec, type) result(decomp)
    integer, intent(in) :: land_decomp_start
    integer, intent(in) :: cstart(:), nap(:)
    type(dim_spec_t), intent(in) :: mem_shape_spec(:), var_shape_spec(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: mem_patch_subscript
    integer :: i, mem_offset, land_index

    mem_patch_subscript = subscript(mem_shape_spec, 'patch')

    mem_shape = [(mem_shape_spec(i)%size, i = 1, size(mem_shape_spec))]
    grid_shape = [(var_shape_spec(i)%size, i = 1, size(var_shape_spec))]

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      land_index = patch_land_index(cstart, nap, mem_index(mem_patch_subscript))
      do i = 1, size(var_shape_spec)
        select case (var_shape_spec(i)%name)
        case ('land')
          grid_index(i) = land_decomp_start + land_index - 1
        case ('patch')
          grid_index(i) = mem_index(mem_patch_subscript) - cstart(land_index) + 1
        case default
          grid_index(i) = mem_index(subscript(mem_shape_spec, var_shape_spec(i)%name))
        end select
      end do
      compmap(mem_offset) = array_offset(grid_index, grid_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, grid_shape, type)

  end function

  function io_decomp_patch_to_patch(patch_decomp_start, mem_shape_spec, var_shape_spec, type) result(decomp)
    integer, intent(in) :: patch_decomp_start
    type(dim_spec_t), intent(in) :: mem_shape_spec(:), var_shape_spec(:)
    integer, intent(in) :: type
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: mem_patch_subscript
    integer :: i, mem_offset

    mem_patch_subscript = subscript(mem_shape_spec, 'patch')

    mem_shape = [(mem_shape_spec(i)%size, i = 1, size(mem_shape_spec))]
    grid_shape = [(var_shape_spec(i)%size, i = 1, size(var_shape_spec))]

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape_spec)
        select case (var_shape_spec(i)%name)
        case ('patch')
          grid_index(i) = patch_decomp_start + mem_index(mem_patch_subscript) - 1
        case default
          grid_index(i) = mem_index(subscript(mem_shape_spec, var_shape_spec(i)%name))
        end select
      end do
      compmap(mem_offset) = array_offset(grid_index, grid_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, grid_shape, type)

  end function

end module
