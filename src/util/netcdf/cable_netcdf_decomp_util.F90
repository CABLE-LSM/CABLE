! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_netcdf_decomp_util_mod
  !* Utilities for generating parallel I/O decompositions for grids used by CABLE.
  !
  ! Parallel I/O decompositions describe the mapping from the local in-memory
  ! layout of an array on the current process to the global layout of a netCDF
  ! variable on disk. For example, a common I/O decomposition in cable is the
  ! mapping of a 1-dimensional array of size `mp` to a netCDF variable with
  ! dimensions `(x, y, patch)` where `x` and `y` are the spatial coordinates of
  ! the land point corresponding to the patch index, and `patch` is the relative
  ! patch index on that land point. Please see the documentation of the individual
  ! decomposition generating functions in this module for more details on the
  ! specific mappings they generate.

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_create_decomp

  use cable_array_utils_mod, only: array_index
  use cable_array_utils_mod, only: array_offset

  use cable_error_handler_mod, only: cable_abort

  implicit none

  private

  public io_decomp_land_to_x_y
  public io_decomp_patch_to_x_y_patch
  public io_decomp_land_to_land
  public io_decomp_patch_to_land_patch
  public io_decomp_patch_to_patch

  integer, parameter :: INDEX_NOT_FOUND = -1

contains

  integer function patch_land_index(cstart, nap, patch_index)
    !* Returns the land index corresponding to the given patch index, using
    ! `cstart` and `nap` to determine the mapping. If the patch index does not lie on
    ! any land point, the function aborts.
    integer, intent(in) :: cstart(:) !! The starting patch index for each land point.
    integer, intent(in) :: nap(:) !! The number of active patches for each land point.
    integer, intent(in) :: patch_index !! The patch index to map to a land index.
    integer i
    patch_land_index = INDEX_NOT_FOUND
    do i = 1, size(cstart)
      if (patch_index >= cstart(i) .and. patch_index <= cstart(i) + nap(i) - 1) then
        patch_land_index = i
        exit
      end if
    end do
    if (patch_land_index == INDEX_NOT_FOUND) then
      call cable_abort("Patch index does not lie on any land point.", file=__FILE__, line=__LINE__)
    end if
  end function

  function io_decomp_land_to_x_y(land_x, land_y, mem_shape, var_shape, type, var_x_index, var_y_index, index_map) result(decomp)
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! 'land' dimension to a netCDF variable layout with 'x' and 'y' dimensions.
    !
    ! The ‘land’ dimension is assumed to be the first (i.e. fastest varying) index
    ! in `mem_shape`.
    integer, intent(in) :: land_x(:) !! An array mapping land indexes to x (longitude) indexes.
    integer, intent(in) :: land_y(:) !! An array mapping land indexes to y (latitude) indexes.
    integer, intent(in) :: mem_shape(:) !! The shape of the in-memory array.
    integer, intent(in) :: var_shape(:) !! The shape of the netCDF variable.
    integer, intent(in) :: type
      !* The data type of the in-memory array for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    integer, intent(in), optional :: var_x_index
      !! The index of the 'x' dimension in `var_shape`. Defaults to 1 if not provided.
    integer, intent(in), optional :: var_y_index
      !! The index of the 'y' dimension in `var_shape`. Defaults to 2 if not provided.
    integer, intent(in), optional :: index_map(:)
      !* An optional mapping from the dimension indexes of `var_shape` to the
      ! dimension indexes of `mem_shape`. If not provided, it is assumed that the
      ! dimensions of `var_shape` map to the dimensions of `mem_shape` in order,
      ! with `var_x_index` and `var_y_index` indexes mapped to the first index
      ! of `mem_shape`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), var_index(:), index_map_local(:)
    integer :: i, var_x_index_local, var_y_index_local, mem_offset

    if (size(var_shape) < 2) call cable_abort("var_shape must have at least 2 dimensions.", __FILE__, __LINE__)

    var_x_index_local = 1
    if (present(var_x_index)) var_x_index_local = var_x_index
    var_y_index_local = 2
    if (present(var_y_index)) var_y_index_local = var_y_index
    index_map_local = [1, (i, i = 1, size(var_shape) - 1)]
    if (present(index_map)) index_map_local = index_map

    allocate(mem_index, mold=mem_shape)
    allocate(var_index, mold=var_shape)
    allocate(compmap(product(mem_shape)))

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape)
        if (i == var_x_index_local) then
          var_index(i) = land_x(mem_index(index_map_local(i)))
        else if (i == var_y_index_local) then
          var_index(i) = land_y(mem_index(index_map_local(i)))
        else
          var_index(i) = mem_index(index_map_local(i))
        end if
      end do
      compmap(mem_offset) = array_offset(var_index, var_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, var_shape, type)

  end function

  function io_decomp_patch_to_x_y_patch(land_x, land_y, cstart, nap, mem_shape, var_shape, type, var_x_index, var_y_index, var_patch_index, index_map) result(decomp)
    !* Returns a parallel I/O decomposition mapping from a memory layout with
    ! a 'patch' dimension to a netCDF variable layout with 'x', 'y', and 'patch'
    ! dimensions.
    !
    ! Note that the 'patch' dimension in the memory layout represents the index in
    ! the 1-dimensional vector of patches, and in the netCDF variable layout, the
    ! 'patch' dimension represents the index of the patch on a particular land
    ! point.
    !
    ! The 'patch' dimension is assumed to be the first (i.e. fastest varying)
    ! index in `mem_shape`.
    integer, intent(in) :: land_x(:) !! An array mapping land indexes to x (longitude) indexes.
    integer, intent(in) :: land_y(:) !! An array mapping land indexes to y (latitude) indexes.
    integer, intent(in) :: cstart(:) !! The starting patch index for each land point.
    integer, intent(in) :: nap(:) !! The number of active patches for each land point.
    integer, intent(in) :: mem_shape(:) !! The shape of the in-memory array.
    integer, intent(in) :: var_shape(:) !! The shape of the netCDF variable.
    integer, intent(in) :: type
      !* The data type of the in-memory array for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    integer, intent(in), optional :: var_x_index
      !! The index of the 'x' dimension in `var_shape`. Defaults to 1 if not provided.
    integer, intent(in), optional :: var_y_index
      !! The index of the 'y' dimension in `var_shape`. Defaults to 2 if not provided.
    integer, intent(in), optional :: var_patch_index
      !! The index of the 'patch' dimension in `var_shape`. Defaults to 3 if not provided.
    integer, intent(in), optional :: index_map(:)
      !* An optional mapping from the dimension indexes of `var_shape` to the
      ! dimension indexes of `mem_shape`. If not provided, it is assumed that the
      ! dimensions of `var_shape` map to the dimensions of `mem_shape` in order,
      ! with `var_x_index`, `var_y_index` and `var_patch_index` being mapped to
      ! the first index of `mem_shape`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), var_index(:), index_map_local(:)
    integer :: i, mem_offset, var_x_index_local, var_y_index_local, var_patch_index_local

    if (size(var_shape) < 3) call cable_abort("var_shape must have at least 3 dimensions.", __FILE__, __LINE__)

    var_x_index_local = 1
    if (present(var_x_index)) var_x_index_local = var_x_index
    var_y_index_local = 2
    if (present(var_y_index)) var_y_index_local = var_y_index
    var_patch_index_local = 3
    if (present(var_patch_index)) var_patch_index_local = var_patch_index
    index_map_local = [1, 1, (i, i = 1, size(var_shape) - 2)]
    if (present(index_map)) index_map_local = index_map

    allocate(mem_index, mold=mem_shape)
    allocate(var_index, mold=var_shape)
    allocate(compmap(product(mem_shape)))

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape)
        if (i == var_x_index_local) then
          var_index(i) = land_x(patch_land_index(cstart, nap, mem_index(index_map_local(i))))
        else if (i == var_y_index_local) then
          var_index(i) = land_y(patch_land_index(cstart, nap, mem_index(index_map_local(i))))
        else if (i == var_patch_index_local) then
          var_index(i) = mem_index(index_map_local(i)) - cstart(patch_land_index(cstart, nap, mem_index(index_map_local(i)))) + 1
        else
          var_index(i) = mem_index(index_map_local(i))
        end if
      end do
      compmap(mem_offset) = array_offset(var_index, var_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, var_shape, type)

  end function

  function io_decomp_land_to_land(land_decomp_start, mem_shape, var_shape, type, var_land_index, index_map) result(decomp)
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! local 'land' dimension to a netCDF variable layout with a global 'land'
    ! dimension.
    !
    ! The 'land' dimension is assumed to be the first (i.e. fastest varying)
    ! index in `mem_shape`.
    integer, intent(in) :: land_decomp_start
      !! The starting index of the first local 'land' index along the global 'land' dimension.
    integer, intent(in) :: mem_shape(:) !! The shape of the in-memory array.
    integer, intent(in) :: var_shape(:) !! The shape of the netCDF variable.
    integer, intent(in) :: type
      !* The data type of the in-memory array for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    integer, intent(in), optional :: var_land_index
      !! The index of the 'land' dimension in `var_shape`. Defaults to 1 if not provided.
    integer, intent(in), optional :: index_map(:)
      !* An optional mapping from the dimension indexes of `var_shape` to the
      ! dimension indexes of `mem_shape`. If not provided, it is assumed that the
      ! dimensions of `var_shape` map to the dimensions of `mem_shape` in order.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), var_index(:), index_map_local(:)
    integer :: i, mem_offset, var_land_index_local

    if (size(var_shape) < 1) call cable_abort("var_shape must have at least 1 dimension.", __FILE__, __LINE__)

    var_land_index_local = 1
    if (present(var_land_index)) var_land_index_local = var_land_index
    index_map_local = [(i, i = 1, size(var_shape))]
    if (present(index_map)) index_map_local = index_map

    allocate(mem_index, mold=mem_shape)
    allocate(var_index, mold=var_shape)
    allocate(compmap(product(mem_shape)))

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape)
        if (i == var_land_index_local) then
          var_index(i) = land_decomp_start + mem_index(index_map_local(i)) - 1
        else
          var_index(i) = mem_index(index_map_local(i))
        end if
      end do
      compmap(mem_offset) = array_offset(var_index, var_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, var_shape, type)

  end function

  function io_decomp_patch_to_land_patch(land_decomp_start, cstart, nap, mem_shape, var_shape, type, var_land_index, var_patch_index, index_map) result(decomp)
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! 'patch' dimension to a netCDF variable layout with a 'land' and 'patch'
    ! dimension.
    !
    ! Note that the 'patch' dimension in the memory layout represents the index in
    ! the 1-dimensional vector of patches, and in the netCDF variable layout, the
    ! 'patch' dimension represents the index of the patch on a particular land
    ! point.
    !
    ! The 'patch' dimension is assumed to be the first (i.e. fastest varying)
    ! index in `mem_shape`.
    integer, intent(in) :: land_decomp_start
      !! The starting index of the first local 'land' index along the global 'land' dimension.
    integer, intent(in) :: cstart(:) !! The starting patch index for each land point.
    integer, intent(in) :: nap(:) !! The number of active patches for each land point.
    integer, intent(in) :: mem_shape(:) !! The shape of the in-memory array.
    integer, intent(in) :: var_shape(:) !! The shape of the netCDF variable.
    integer, intent(in) :: type
      !* The data type of the in-memory array for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    integer, intent(in), optional :: var_land_index
      !! The index of the 'land' dimension in `var_shape`. Defaults to 1 if not provided.
    integer, intent(in), optional :: var_patch_index
      !! The index of the 'patch' dimension in `var_shape`. Defaults to 2 if not provided.
    integer, intent(in), optional :: index_map(:)
      !* An optional mapping from the dimension indexes of `var_shape` to the
      ! dimension indexes of `mem_shape`. If not provided, it is assumed that the
      ! dimensions of `var_shape` map to the dimensions of `mem_shape` in order,
      ! with `var_land_index` and `var_patch_index` being mapped to the first
      ! index of `mem_shape`.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), var_index(:), index_map_local(:)
    integer :: i, mem_offset, var_land_index_local, var_patch_index_local

    if (size(var_shape) < 2) call cable_abort("var_shape must have at least 2 dimensions.", __FILE__, __LINE__)

    var_land_index_local = 1
    if (present(var_land_index)) var_land_index_local = var_land_index
    var_patch_index_local = 2
    if (present(var_patch_index)) var_patch_index_local = var_patch_index
    index_map_local = [1, (i, i = 1, size(var_shape) - 1)]
    if (present(index_map)) index_map_local = index_map

    allocate(mem_index, mold=mem_shape)
    allocate(var_index, mold=var_shape)
    allocate(compmap(product(mem_shape)))

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape)
        if (i == var_land_index_local) then
          var_index(i) = land_decomp_start + patch_land_index(cstart, nap, mem_index(index_map_local(i))) - 1
        else if (i == var_patch_index_local) then
          var_index(i) = mem_index(index_map_local(i)) - cstart(patch_land_index(cstart, nap, mem_index(index_map_local(i)))) + 1
        else
          var_index(i) = mem_index(index_map_local(i))
        end if
      end do
      compmap(mem_offset) = array_offset(var_index, var_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, var_shape, type)

  end function

  function io_decomp_patch_to_patch(patch_decomp_start, mem_shape, var_shape, type, var_patch_index, index_map) result(decomp)
    !* Returns a parallel I/O decomposition mapping from a memory layout with a
    ! local 'patch' dimension to a netCDF variable layout with a global 'patch'
    ! dimension.
    !
    ! The 'patch' dimension is assumed to be the first (i.e. fastest varying)
    ! index in `mem_shape`.
    integer, intent(in) :: patch_decomp_start
      !! The starting index of the first local 'patch' index along the global 'patch' dimension.
    integer, intent(in) :: mem_shape(:) !! The shape of the in-memory array.
    integer, intent(in) :: var_shape(:) !! The shape of the netCDF variable.
    integer, intent(in) :: type
      !* The data type of the in-memory array for which the decomposition is being
      ! created using `CABLE_NETCDF_TYPE_*` constants from `cable_netcdf_mod`.
    integer, intent(in), optional :: var_patch_index
      ! The index of the 'patch' dimension in `var_shape`. Defaults to 1 if not provided.
    integer, intent(in), optional :: index_map(:)
      !* An optional mapping from the dimension indexes of `var_shape` to the
      ! dimension indexes of `mem_shape`. If not provided, it is assumed that the
      ! dimensions of `var_shape` map to the dimensions of `mem_shape` in order.
    class(cable_netcdf_decomp_t), allocatable :: decomp

    integer, allocatable :: compmap(:)
    integer, allocatable :: mem_index(:), var_index(:), index_map_local(:)
    integer :: i, mem_offset, var_patch_index_local

    if (size(var_shape) < 1) call cable_abort("var_shape must have at least 1 dimension.", __FILE__, __LINE__)

    var_patch_index_local = 1
    if (present(var_patch_index)) var_patch_index_local = var_patch_index
    index_map_local = [(i, i = 1, size(var_shape))]
    if (present(index_map)) index_map_local = index_map

    allocate(mem_index, mold=mem_shape)
    allocate(var_index, mold=var_shape)
    allocate(compmap(product(mem_shape)))

    do mem_offset = 1, size(compmap)
      call array_index(mem_offset, mem_shape, mem_index)
      do i = 1, size(var_shape)
        if (i == var_patch_index_local) then
          var_index(i) = patch_decomp_start + mem_index(index_map_local(i)) - 1
        else
          var_index(i) = mem_index(index_map_local(i))
        end if
      end do
      compmap(mem_offset) = array_offset(var_index, var_shape)
    end do

    decomp = cable_netcdf_create_decomp(compmap, var_shape, type)

  end function

end module
