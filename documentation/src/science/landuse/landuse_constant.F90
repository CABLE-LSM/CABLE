MODULE landuse_constant
  USE cable_def_types_mod,    ONLY: mvtype,mstype,mland,ms,msn,nrb,ncp,ncs,r_2
  USE casadimension,          ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  IMPLICIT NONE
!  integer, parameter                     :: sp =selected_real_kind(8)
!  integer, parameter                     :: dp =selected_real_kind(16)
  !  
  INTEGER,   PARAMETER                   :: mstate   = 12  ! number of land use states
  INTEGER,   PARAMETER                   :: mvmax    = 17
  INTEGER,   PARAMETER                   :: mharvw   = 5           

  real(r_2),    PARAMETER                   :: thresh_frac  = 1.0e-6
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fwoodprod    =(/0.3,0.4,0.4/)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: ratewoodprod =(/1.0,0.1,0.01/)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fracwoodseed =(/0.4,0.3,0.3/)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fracgrassseed=(/0.6,0.0,0.4/)
  REAL(r_2),    PARAMETER                   :: fseedling = 0.001
END MODULE landuse_constant
