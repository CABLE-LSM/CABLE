MODULE landuse_constant
  !! The `landuse_constant` module contains landuse-specific constants.
  !!
  !! **WARNING:** Eventually, `mvmax` in the land-use code should be replaced
  !! with `mvtype`.
  !written by YP Wang
  USE cable_def_types_mod,    ONLY: mvtype,mstype,mland,ms,msn,nrb,ncp,ncs,r_2
  USE casadimension,          ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  IMPLICIT NONE
  !  integer, parameter                     :: sp =selected_real_kind(8)
  !  integer, parameter                     :: dp =selected_real_kind(16)
  INTEGER,   PARAMETER                   :: mstate   = 12
  !* Number of land use states as in HYDE land use dataset processed by
  ! [George Hurtt](https://doi.org/10.5194/gmd-13-5425-2020)
  ! \( (-) \)
  INTEGER,   PARAMETER                   :: mvmax    = 17
  !! Maximum number of plant functional types in CABLE  \( (-) \)
  INTEGER,   PARAMETER                   :: mharvw   = 5           
  !! Number of harvest types \( (-) \)
  real(r_2),    PARAMETER                   :: thresh_frac  = 1.0e-6
  !! Threshold area fraction for a separate tile within a land cell \( (-) \)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fwoodprod    =(/0.3,0.4,0.4/)
  !! Allocation of harvested wood into each wood product pools \( (-) \)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: ratewoodprod =(/1.0,0.1,0.01/)
  !! Decay rate of each wood product pools \( ( year^{-1} ) \)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fracwoodseed =(/0.4,0.3,0.3/)
  !! Biomass fraction in leaf, wood and root pools for woody seedlings in a cleared area \( (-) \)
  REAL(r_2),    PARAMETER, DIMENSION(mwood) :: fracgrassseed=(/0.6,0.0,0.4/)
  !! Biomass fraction in leaf wood, and root pools for grassy seedlings in a cleared area \( (-) \)
  REAL(r_2),    PARAMETER                   :: fseedling = 0.001
  !! Total biomass of seedlings in a cleared area \( ( g C \cdot m^{-2} ) \)
END MODULE landuse_constant
