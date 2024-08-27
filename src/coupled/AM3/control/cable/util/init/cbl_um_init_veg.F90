MODULE cbl_um_init_veg_mod

IMPLICIT NONE

PUBLIC initialize_veg

CONTAINS

SUBROUTINE initialize_veg( SurfaceType, SoilType, mp,&
                           ms, nrb, npft, nsurft, land_pts, l_tile_pts,        &
                           ICE_SurfaceType, ICE_SoilType, LAI_thresh,          &
                           soil_zse, surft_pts, surft_index, tile_frac,        &
                           veg, soil, pars )

USE params_io_mod_cbl,    ONLY: params_io_data_type
USE cable_def_types_mod,  ONLY: veg_parameter_type
USE cable_def_types_mod,  ONLY: soil_parameter_type, r_2

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: ms
INTEGER, INTENT(IN) :: nrb
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: npft
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: ICE_SurfaceType               ! index ICE surface type
INTEGER, INTENT(IN) :: ICE_SoilType                  ! index soil type
REAL,    INTENT(IN) :: lai_thresh 
REAL,    INTENT(IN) :: soil_zse(ms)                  ! soil depth per layer 
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)

INTEGER, INTENT(IN) :: SurfaceType(mp)              ! surface tile PFT/nveg
INTEGER, INTENT(IN) :: SoilType(mp)                 ! soil type per tile

TYPE(veg_parameter_type),  INTENT(OUT) :: veg        ! vegetation parameters
TYPE(soil_parameter_type), INTENT(OUT) :: soil       ! soil parameters
TYPE(params_io_data_type), INTENT(IN)  :: pars

! local vars
INTEGER :: i 
REAL    :: totdepth
   
veg%meth  = 1               ! review !should be namelist config if even used
veg%ejmax = 2.0 * veg%vcmax
veg%gamma = 3.e-2           ! for Haverd2013 switch ! offline init _parameters
veg%F10   = 0.85            ! offline set init _parameters

! calculate veg%froot from using rootbeta and soil depth
! (Jackson et al. 1996, Oceologica, 108:389-411)
totdepth = 0.0
DO i = 1, ms-1
  totdepth = totdepth + soil_zse(i) * 100.0  ! unit in centimetres
  veg%froot(:, i) = MIN( 1.0_r_2, 1.0-veg%rootbeta(:)**totdepth )
END DO
veg%froot(:, ms) = 1.0 - veg%froot(:, ms-1)
DO i = ms-1, 2, -1
  veg%froot(:, i) = veg%froot(:, i)-veg%froot(:,i-1)
END DO

RETURN
END SUBROUTINE initialize_veg

END MODULE cbl_um_init_veg_mod




