! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_netcdf_decomp_util_mod
  !! Utilities for generating parallel I/O decompositions for grids used by CABLE.

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_create_decomp

  use cable_array_utils_mod, only: array_index
  use cable_array_utils_mod, only: array_offset

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

  integer, parameter :: SUBSCRIPT_NOT_FOUND = -1
  type dim_spec_t
    !* A derived type which associates a name with a dimension.
    ! The `name` component can be used to identify a dimension which is used for
    ! mapping a dimension in the in-memory array to one or more dimensions in the
    ! netCDF variable data. For example, the patch dimension which ranges from `1`
    ! to `mp` can be mapped to a `(land, patch)` coordinate which ranges from `(1,
    ! 1)` to `(mland_global, max_veg_patches)`.
    character(64) :: name !! The name of the dimension.
    integer :: size !! The size of the dimension.
  end type

contains

  integer function subscript(shape_spec, name)
    !* Returns the subscript of the dimension matching `name` in the
    ! `shape_spec` array. If no such dimension exists, the function aborts.
    type(dim_spec_t), intent(in) :: shape_spec(:) !! The shape_spec array to search for the name.
    character(*), intent(in) :: name !! The name of the dimension to find the subscript of.
    integer i
    subscript = SUBSCRIPT_NOT_FOUND
    do i = 1, size(shape_spec)
      if (shape_spec(i)%name == name) then
        subscript = i
        exit
      end if
    end do
    if (subscript == SUBSCRIPT_NOT_FOUND) then
      call cable_abort("Name '" // name // "' not found in shape_spec.", __FILE__, __LINE__)
    end if
  end function

  integer function patch_land_index(cstart, nap, patch_index)
    !* Returns the land index corresponding to the given patch index, using
    ! `cstart` and `nap` to determine the mapping. If the patch index does not lie on
    ! any land point, the function aborts.
    integer, intent(in) :: cstart(:) !! The starting patch index for each land point.
    integer, intent(in) :: nap(:) !! The number of active patches for each land point.
    integer, intent(in) :: patch_index !! The patch index to map to a land index.
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
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! 'land' dimension to a netCDF variable layout with 'x' and 'y' dimensions,
    ! using the provided `land_x` and `land_y` arrays to determine the mapping
    ! from land indexes to x and y indexes.
    integer, intent(in) :: land_x(:) !! An array mapping land indexes to x (longitude) indexes.
    integer, intent(in) :: land_y(:) !! An array mapping land indexes to y (latitude) indexes.
    type(dim_spec_t), intent(in) :: mem_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! in-memory array. `mem_shape_spec` must include a dimension with name
      ! 'land' which is used to map to the `x` and `y` dimensions.
    type(dim_spec_t), intent(in) :: var_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! netCDF variable. `var_shape_spec` must include dimensions with names
      ! `x` and `y` which are used to map from the 'land' dimension described by
      ! `mem_shape_spec`.
    integer, intent(in) :: type
      !* The data type of the variable for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: i, mem_offset, mem_land_subscript

    mem_land_subscript = subscript(mem_shape_spec, 'land')

    mem_shape = mem_shape_spec(:)%size
    grid_shape = var_shape_spec(:)%size

    allocate(mem_index, mold=mem_shape)
    allocate(grid_index, mold=grid_shape)
    allocate(compmap(product(mem_shape)))

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
    !* Returns a parallel I/O decomposition mapping from a memory layout with
    ! a 'patch' dimension to a netCDF variable layout with 'x', 'y', and 'patch'
    ! dimensions, using the provided `land_x` and `land_y` arrays to determine the
    ! mapping from land indexes to x and y indexes, and using `cstart` and `nap`
    ! to determine the mapping from patch indexes to `(land, patch)` coordinates.
    !
    ! Note that the 'patch' dimension in the memory layout represents the index in
    ! the 1-dimensional vector of patches, and in the netCDF variable layout, the
    ! 'patch' dimension represents the index of the patch on a particular land
    ! point. The mapping from patch index to land point is determined by `cstart`
    ! and `nap`, and the mapping from land point to x and y coordinates is
    ! determined by `land_x` and `land_y`.
    integer, intent(in) :: land_x(:) !! An array mapping land indexes to x (longitude) indexes.
    integer, intent(in) :: land_y(:) !! An array mapping land indexes to y (latitude) indexes.
    integer, intent(in) :: cstart(:) !! The starting patch index for each land point.
    integer, intent(in) :: nap(:) !! The number of active patches for each land point.
    type(dim_spec_t), intent(in) :: mem_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! in-memory array. `mem_shape_spec` must include a dimension with name
      ! 'patch' which is used to map to the `x`, `y` and `patch` dimensions.
    type(dim_spec_t), intent(in) :: var_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! netCDF variable. `var_shape_spec` must include dimensions with names
      ! `x`, `y` and `patch` which are used to map from the 'land' dimension
      ! described by `mem_shape_spec`.
    integer, intent(in) :: type
      !* The data type of the variable for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: i, mem_offset, mem_patch_subscript, land_index

    mem_patch_subscript = subscript(mem_shape_spec, 'patch')

    mem_shape = mem_shape_spec(:)%size
    grid_shape = var_shape_spec(:)%size

    allocate(mem_index, mold=mem_shape)
    allocate(grid_index, mold=grid_shape)
    allocate(compmap(product(mem_shape)))

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
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! local 'land' dimension to a netCDF variable layout with a global 'land'
    ! dimension.
    integer, intent(in) :: land_decomp_start
      !! The starting index of the first local 'land' index along global 'land' dimension.
    type(dim_spec_t), intent(in) :: mem_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! in-memory array. `mem_shape_spec` must include a dimension with name
      ! 'land' which is used to map to the global `land` dimension.
    type(dim_spec_t), intent(in) :: var_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! netCDF variable. `var_shape_spec` must include a dimension with name
      ! `land` which is used to map from the 'land' dimension described by
      ! `mem_shape_spec`.
    integer, intent(in) :: type
      !* The data type of the variable for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: i, mem_offset, mem_land_subscript

    mem_land_subscript = subscript(mem_shape_spec, 'land')

    mem_shape = mem_shape_spec(:)%size
    grid_shape = var_shape_spec(:)%size

    allocate(mem_index, mold=mem_shape)
    allocate(grid_index, mold=grid_shape)
    allocate(compmap(product(mem_shape)))

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
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! local 'patch' dimension to a netCDF variable layout with a global 'land'
    ! dimension and a 'patch' dimension.
    !
    ! Note that the 'patch' dimension in the memory layout represents the index in
    ! the 1-dimensional vector of patches, and in the netCDF variable layout, the
    ! 'patch' dimension represents the index of the patch on a particular land
    ! point. The mapping from patch index to land point is determined by `cstart`
    ! and `nap`, and the mapping from land point to global land index is
    ! determined by `land_decomp_start`.
    integer, intent(in) :: land_decomp_start
      !! The starting index of the first local 'land' index along global 'land' dimension.
    integer, intent(in) :: cstart(:) !! The starting patch index for each land point.
    integer, intent(in) :: nap(:) !! The number of active patches for each land point.
    type(dim_spec_t), intent(in) :: mem_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! in-memory array. `mem_shape_spec` must include a dimension with name
      ! 'patch' which is used to map to the global `land` and `patch` dimensions.
    type(dim_spec_t), intent(in) :: var_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! netCDF variable. `var_shape_spec` must include dimensions with names
      ! `land` and `patch` which are used to map from the 'patch' dimension
      ! described by `mem_shape_spec`.
    integer, intent(in) :: type
      !* The data type of the variable for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: i, mem_offset, mem_patch_subscript, land_index

    mem_patch_subscript = subscript(mem_shape_spec, 'patch')

    mem_shape = mem_shape_spec(:)%size
    grid_shape = var_shape_spec(:)%size

    allocate(mem_index, mold=mem_shape)
    allocate(grid_index, mold=grid_shape)
    allocate(compmap(product(mem_shape)))

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
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! local 'patch' dimension to a netCDF variable layout with a global 'patch'
    ! dimension.
    !
    ! Note that the 'patch' dimension in the memory layout represents the local
    ! index in the 1-dimensional vector of patches, and in the netCDF variable
    ! layout, the 'patch' dimension represents the global index in the
    ! 1-dimensional vector of patches.
    integer, intent(in) :: patch_decomp_start
      !! The starting index of the first local 'patch' index along global 'patch' dimension.
    type(dim_spec_t), intent(in) :: mem_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! in-memory array. `mem_shape_spec` must include a dimension with name
      ! 'patch' which is used to map to the global `patch` dimensions.
    type(dim_spec_t), intent(in) :: var_shape_spec(:)
      !* An array of dim_spec_t describing the shape and dimension names of the
      ! netCDF variable. `var_shape_spec` must include a dimension with name
      ! `patch` which is used to map from the local 'patch' dimension described by
      ! `mem_shape_spec`.
    integer, intent(in) :: type
      !* The data type of the variable for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), mem_shape(:)
    integer, allocatable :: grid_index(:), grid_shape(:)
    integer :: i, mem_offset, mem_patch_subscript

    mem_patch_subscript = subscript(mem_shape_spec, 'patch')

    mem_shape = mem_shape_spec(:)%size
    grid_shape = var_shape_spec(:)%size

    allocate(mem_index, mold=mem_shape)
    allocate(grid_index, mold=grid_shape)
    allocate(compmap(product(mem_shape)))

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
