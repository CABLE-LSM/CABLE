module cable_array_utils_mod
  implicit none

contains

  function array_offset(index, shape) result(offset)
    integer, intent(in) :: index(:), shape(:)
    integer :: i, offset, shape_factor
    offset = 1
    shape_factor = 1
    do i = 1, size(index)
      offset = offset + (index(i) - 1) * shape_factor
      shape_factor = shape_factor * shape(i)
    end do
  end function

  subroutine array_index(offset_in, shape, index)
    integer, intent(in) :: offset_in, shape(:)
    integer, intent(inout) :: index(:)
    integer :: i, offset
    offset = offset_in
    do i = 1, size(shape)
      index(i) = mod(offset - 1, shape(i)) + 1
      offset = (offset - 1) / shape(i) + 1
    end do
  end subroutine

end module cable_array_utils_mod