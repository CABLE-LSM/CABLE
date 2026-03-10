! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_array_utils_mod
  !! Utility procedures for working with arrays.
  implicit none
  private

  public array_offset
  public array_index

contains

  !> Calculate the memory offset corresponding to a given index in a multi-dimensional array.
  function array_offset(index, shape) result(offset)
    integer, intent(in) :: index(:) !! The multi-dimensional index for which to calculate the offset.
    integer, intent(in) :: shape(:) !! The shape of the multi-dimensional array.
    integer :: offset !! The calculated memory offset corresponding to the given index.
    integer :: i, scale
    offset = 1; scale = 1
    do i = 1, size(index)
      offset = offset + (index(i) - 1) * scale
      scale = scale * shape(i)
    end do
  end function

  !> Calculate the index corresponding to a given memory offset in a multi-dimensional array.
  subroutine array_index(offset_in, shape, index)
    integer, intent(in) :: offset_in !! The memory offset for which to calculate the multi-dimensional index.
    integer, intent(in) :: shape(:) !! The shape of the multi-dimensional array.
    integer, intent(inout) :: index(:) !! The calculated multi-dimensional index corresponding to the given offset.
    integer :: i, offset
    offset = offset_in
    do i = 1, size(shape)
      index(i) = mod(offset - 1, shape(i)) + 1
      offset = (offset - 1) / shape(i) + 1
    end do
  end subroutine

end module cable_array_utils_mod
