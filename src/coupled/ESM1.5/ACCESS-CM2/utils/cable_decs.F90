module cable_decs_mod

  type CABLE_file
    logical :: want = .false. 
    integer :: funit
    character(len=199) :: folder
    character(len=190) :: filename
    logical :: vars = .false.
    integer :: itau=0,ftau=0,mtau=0 
  End type CABLE_file
  
  LOGICAL, allocatable, DIMENSION(:,:), save :: L_tile_pts
  
  real, allocatable, DIMENSION(:,:,:), save ::sw_down_diag
  
  !___UM parameters: water density
  REAL, parameter :: rho_water = 1000.

End module cable_decs_mod
