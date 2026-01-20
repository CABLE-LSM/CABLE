module cable_enum_mod
  implicit none

  type :: cable_enum_t
    integer :: value
  contains
    procedure, private :: cable_enum_eq
    generic :: operator(==) => cable_enum_eq
    procedure, private :: cable_enum_ne
    generic :: operator(/=) => cable_enum_ne
  end type

contains

  elemental logical function cable_enum_eq(this, other)
    class(cable_enum_t), intent(in) :: this, other
    cable_enum_eq = (this%value == other%value)
  end function

  elemental logical function cable_enum_ne(this, other)
    class(cable_enum_t), intent(in) :: this, other
    cable_enum_ne = (this%value /= other%value)
  end function

end module
