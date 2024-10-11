MODULE cbl_um_update_mod
   
IMPLICIT NONE   
PUBLIC :: update_data
PRIVATE  

CONTAINS

SUBROUTINE update_data( row_length, rows, land_pts, nsurft, npft, ms, msn,     &
                        timestep, timestep_number, mp, nrb, CO2_MMR,           & 
                        HGT_pft_um, lai_pft_um, land_index, surft_pts,         &
                        surft_index, tile_frac, L_tile_pts, cos_zenith_angle,  &
                        latitude, longitude, sw_down_VIS, sw_down_NIR,         &
                        beamFrac_VIS, beamFrac_NIR, beamFrac_TOT, lw_down,     &
                        ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq, &
                        z1_uv, canopy_tile, rad, met, veg, soil, rough,canopy, &
                        ssnow, HGT_pft_cbl, LAI_pft_cbl, reducedLAIdue2snow )
                        
USE cbl_um_update_met_mod,       ONLY: update_met
USE cbl_um_update_radiation_mod, ONLY: update_radiation
USE cbl_um_update_roughness_mod, ONLY: update_roughness
USE cbl_um_update_canopy_mod,    ONLY: update_canopy
USE cbl_um_update_soilsnow_mod,    ONLY: update_soilsnow
USE cbl_LAI_canopy_height_mod,     ONLY: limit_HGT_LAI
USE params_io_mod_cbl,             ONLY: params_io_data_type
USE cable_other_constants_mod,     ONLY: LAI_THRESH
USE cable_def_types_mod,           ONLY: radiation_type, met_type,             &
                                         veg_parameter_type,                   &
                                         soil_parameter_type, roughness_type,  &
                                         canopy_type, soil_snow_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: npft              ! # plant functional types
INTEGER, INTENT(IN) :: ms                ! # soil layers 
INTEGER, INTENT(IN) :: msn               ! # snow layers 
REAL,    INTENT(IN) :: timestep
INTEGER, INTENT(IN) :: timestep_number     
INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: nrb
REAL,    INTENT(IN) :: co2_mmr
REAL,    INTENT(IN) :: HGT_pft_um(land_pts, npft) 
REAL,    INTENT(IN) :: lai_pft_um(land_pts, npft)
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)          ! tangled cell index of land_pt
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts, nsurft)
REAL,    INTENT(IN) :: cos_zenith_angle(row_length,rows)
REAL,    INTENT(IN) :: latitude(row_length,rows)
REAL,    INTENT(IN) :: longitude(row_length,rows)
REAL,    INTENT(IN) :: lw_down(row_length,rows)
REAL,    INTENT(IN) :: ls_rain(row_length,rows) 
REAL,    INTENT(IN) :: ls_snow(row_length,rows) 
REAL,    INTENT(IN) :: tl_1(row_length,rows)
REAL,    INTENT(IN) :: qw_1(row_length,rows)
REAL,    INTENT(IN) :: vshr_land(row_length,rows)
REAL,    INTENT(IN) :: pstar(row_length,rows)
REAL,    INTENT(IN) :: z1_tq(row_length,rows)
REAL,    INTENT(IN) :: z1_uv(row_length,rows)
REAL,    INTENT(IN) :: sw_down_VIS(row_length,rows)
REAL,    INTENT(IN) :: sw_down_NIR(row_length,rows)
REAL,    INTENT(IN) :: beamFrac_VIS(row_length,rows)
REAL,    INTENT(IN) :: beamFrac_NIR(row_length,rows)
REAL,    INTENT(IN) :: beamFrac_TOT(row_length,rows)
REAL,    INTENT(IN) :: canopy_tile(land_pts, nsurft)
REAL,    INTENT(IN) :: reducedLAIdue2snow(mp)
REAL,    INTENT(OUT) :: HGT_pft_cbl(mp)
REAL,    INTENT(OUT) :: LAI_pft_cbl(mp)

TYPE(soil_parameter_type), INTENT(IN)    :: soil       ! soil parameters
TYPE(radiation_type),      INTENT(OUT)   :: rad
TYPE(met_type),            INTENT(OUT)   :: met
TYPE(roughness_type),      INTENT(OUT)   :: rough
TYPE(canopy_type),         INTENT(OUT)   :: canopy
TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow      ! 
TYPE(veg_parameter_type),  INTENT(INOUT) :: veg        ! vegetation parameters

! Even when prescribed LAI, canopy height are seasonal so this cant be done ONLY
! on first call limit IN height, LAI and initialize veg% equivalents 
CALL limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, nsurft, npft,      &
                    surft_pts, surft_index, tile_frac, l_tile_pts, LAI_pft_um, &
                    HGT_pft_um, LAI_thresh )

! def veg% as they are what is USEd used in core science code
veg%vlai  = LAI_pft_cbl
veg%hc    = HGT_pft_cbl

!--- Met. forcing
CALL update_met( mp, row_length, rows, timestep, land_pts, nsurft,             &
                 surft_pts, surft_index,land_index,L_tile_pts, co2_mmr,        &
                 ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar, met )
                     
CALL update_radiation( mp, row_length, rows, timestep, land_pts, nsurft,       &
                       surft_pts, surft_index,land_index, L_tile_pts,          &
                       latitude, longitude, cos_zenith_angle, sw_down_VIS,     &
                       sw_down_NIR, beamFrac_VIS, beamFrac_NIR,                &
                       beamFrac_TOT, lw_down, rad, met )

CALL update_roughness( row_length, rows, mp, land_pts, nsurft, npft,           &
                       lai_thresh, surft_pts, surft_index, land_index,         &
                       L_tile_pts, z1_tq, z1_uv, HGT_pft_um, HGT_pft_cbl,      &
                       rough, veg )

CALL update_canopy( mp, land_pts, nsurft, L_tile_pts, canopy_tile,                &
                    reducedLAIdue2snow, canopy )

CALL update_soilsnow( mp, soil, ssnow, veg%iveg )

RETURN
END SUBROUTINE update_data
                                       
END MODULE cbl_um_update_mod
