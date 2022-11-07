MODULE landuse_constant
  ! this contains landuse-specific constants
  ! written by YP Wang
  USE cable_def_types_mod,    ONLY: mvtype,mstype,mland,ms,msn,nrb,ncp,ncs,r_2
  USE casadimension,          ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  IMPLICIT NONE
  !  integer, parameter                     :: sp =selected_real_kind(8)
  !  integer, parameter                     :: dp =selected_real_kind(16)
  !> number of land use states as in HYDE land use dataset processed by  
  ! [George Hurtt](https://doi.org/10.5194/gmd-13-5425-2020)
  !
  INTEGER,   PARAMETER                   :: mstate   = 12  ! number of land use states
  !> max number of PFT in CABLE, should be replaced by "mvtype" 
  INTEGER,   PARAMETER                   :: mvmax    = 17
  !> number of harvest types
  INTEGER,   PARAMETER                   :: mharvw   = 5           
  !> threshold area fraction for a separate tile within a land cell
  real(r_2),    PARAMETER                   :: thresh_frac  = 1.0e-6
  !> allocation of harvested wood into different wood product pools
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fwoodprod    =(/0.3,0.4,0.4/)
  !> decay rate of different wood product pools
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: ratewoodprod =(/1.0,0.1,0.01/)
  !> biomass fraction in leaf wood and root for woody seedlings in a cleared area
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fracwoodseed =(/0.4,0.3,0.3/)
  !> biomass fraction in leaf wood and root for grassy seedlings in a cleared area
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fracgrassseed=(/0.6,0.0,0.4/)
  !> total biomass of seeldings in g C/m2 in a cleared area
  REAL(r_2),    PARAMETER                   :: fseedling = 0.001
END MODULE landuse_constant
