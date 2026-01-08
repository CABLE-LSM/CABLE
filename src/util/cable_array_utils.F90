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

  subroutine array_partition(n, k, p, start, count)
    !* Compute start and count for the p'th partition of an array of size n
    ! where p = 0, 1, ... , k - 1.
    !
    ! For k partitions, an array of n elements can be partitioned into r
    ! partitions of length q + 1, and k - r partitions of length q where q and r
    ! are the quotient and remainder of n divided by k (i.e. n = q * k + r).
    ! Note, we assume that the r partitions of length q + 1 precede the k - r
    ! partitions of length q in the array.
    integer, intent(in) :: n, k, p
    integer, intent(out) :: start, count
    integer :: q, r

    q = n / k
    r = mod(n, k)

    if (p < r) then
      count = q + 1
      start = 1 + (q + 1) * p
    else
      count = q
      start = 1 + (q + 1) * r + q * (p - r)
    end if

  end subroutine array_partition

end module cable_array_utils_mod