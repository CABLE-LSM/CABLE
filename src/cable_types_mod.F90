MODULE cable_types_mod
#define UM_BUILD yes

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

!calculated length of cable arrays = number of active tiles
INTEGER :: mp          ! # total no of patches/tiles

LOGICAL, ALLOCATABLE :: l_tile_pts(:,:)

!jh!elevate these to namelist status
INTEGER ::        & 
#ifdef UM_BUILD 
  mvtype=17,      & ! total # vegetation types,   from input
  mstype=9,       & ! total # soil types, needs to de defined atCompile TimeForNow
#else       
  mvtype,         & ! total # vegetation types,   from input
  mstype,         & ! total # soil types,         from input
#endif
  mland                           ! # land grid cells

!jh!elevate these to namelist status
INTEGER, PARAMETER ::                                                        &
       i_d  = KIND(9), &
#ifdef UM_BUILD 
       r_2  = KIND(1.0),&!SELECTED_REAL_KIND(12, 50), &
#else       
       r_2  = KIND(1.d0),&!SELECTED_REAL_KIND(12, 50), &
#endif
       n_tiles = 17,  & ! # possible no of different
       ncp = 3,       & ! # vegetation carbon stores
       ncs = 2,       & ! # soil carbon stores
       mf = 2,        & ! # leaves (sunlit, shaded)
!jh!this needs sorting (nrb=4[2*2])
       nrb = 3,       & ! # radiation bands
       msn = 3,       & ! max # snow layers
       swb = 2,       & ! # shortwave bands
       niter = 4,     & ! number of iterations for za/L
                                !      ms = 12          ! # soil layers
       ms = 6         ! # soil layers - standard

  INTEGER, PARAMETER :: n_ktherm = 3

END MODULE cable_types_mod
