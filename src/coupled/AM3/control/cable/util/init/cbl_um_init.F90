MODULE cbl_um_init_mod
   
IMPLICIT NONE   
PUBLIC :: init_data
PRIVATE  

CONTAINS

SUBROUTINE init_data( row_length, rows, land_pts, nsurft, npft, ms, msn,       &
                      soil_zse, mp, nrb, CO2_MMR, tfrz, ICE_SurfaceType,       &
                      ICE_SoilType, land_index, surft_pts, surft_index,        &
                      tile_frac, L_tile_pts, albsoil, bexp, hcon, satcon,      &
                      sathh, smvcst, smvcwt, smvccl, pars, tl_1, snow_tile,    &
                      SoilTemp, SoilMoisture, FrozenSoilFrac,                  &
                      OneLyrSnowDensity, SnowAge, ThreeLayerSnowFlag,          &
                      SnowDensity, SnowDepth, SnowTemp, SnowMass, rad_trad,    &
                      met_tk, veg, soil, canopy, ssnow, bgc, sum_flux,         &
                      SurfaceType, SoilType, npp_pft_acc, resp_w_pft_acc )
! subrs
USE cbl_um_init_veg_mod,           ONLY: initialize_veg
USE cbl_um_init_soil_mod,          ONLY: initialize_soil
USE cbl_um_init_soilsnow_mod,      ONLY: initialize_soilsnow
USE cable_um_init_respiration_mod, ONLY: init_respiration
USE cable_um_init_bgc_mod,         ONLY: init_bgc_vars
USE cable_um_init_sumflux_mod,     ONLY: init_sumflux_zero
USE cable_pack_mod,                ONLY: cable_pack_rr

! data
USE cable_other_constants_mod,     ONLY: LAI_THRESH
USE grid_constants_mod_cbl,        ONLY: nsnl, nsoil_max
USE progs_cnp_vars_mod,            ONLY: nCpool_casa, nNpool_casa, nPPool_casa
USE cable_def_types_mod,           ONLY: veg_parameter_type, canopy_type,      &
                                         soil_parameter_type, soil_snow_type,  &
                                         bgc_pool_type, sum_flux_type
USE params_io_mod_cbl,             ONLY: params_io_data_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: npft              ! # plant functional types
INTEGER, INTENT(IN) :: ms                ! # soil layers 
INTEGER, INTENT(IN) :: msn               ! # snow layers 
REAL,    INTENT(IN) :: soil_zse(ms)      ! soil layer thicknesses 
REAL,    INTENT(IN) :: tfrz 
INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: nrb
INTEGER, INTENT(IN) :: ICE_SoilType
REAL,    INTENT(IN) :: co2_mmr
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)          ! cell index of land_pt
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts, nsurft)
REAL,    INTENT(IN) :: bexp (land_pts, ms)           ! b in Campbell equation
REAL,    INTENT(IN) :: satcon(land_pts, ms)
                                    ! hydraulic conductivity @ saturation [mm/s]
REAL,    INTENT(IN) :: sathh(land_pts, ms)
REAL,    INTENT(IN) :: smvcst(land_pts, ms)
REAL,    INTENT(IN) :: smvcwt(land_pts, ms)
REAL,    INTENT(IN) :: smvccl(land_pts, ms)
REAL,    INTENT(IN) :: hcon(land_pts)      ! Soil thermal conductivity (W/m/K).
REAL,    INTENT(IN) :: albsoil(land_pts)
REAL,    INTENT(IN) :: tl_1(row_length,rows)
REAL,    INTENT(IN) :: snow_tile(land_pts, nsurft)
REAL,    INTENT(IN) :: SoilTemp(land_pts, nsurft, ms)
REAL,    INTENT(IN) :: SoilMoisture(land_pts, nsurft, ms)
REAL,    INTENT(IN) :: FrozenSoilFrac(land_pts, nsurft, ms)
REAL,    INTENT(IN) :: SnowDepth(land_pts, nsurft,nsnl)
REAL,    INTENT(IN) :: SnowTemp(land_pts, nsurft,nsnl)
REAL,    INTENT(IN) :: SnowMass(land_pts, nsurft,nsnl)
REAL,    INTENT(IN) :: SnowDensity(land_pts, nsurft,nsnl)
REAL,    INTENT(IN) :: OneLyrSnowDensity(land_pts, nsurft)
REAL,    INTENT(IN) :: SnowAge(land_pts, nsurft)
REAL,    INTENT(IN) :: npp_pft_acc(land_pts,npft)
REAL,    INTENT(IN) :: resp_w_pft_acc(land_pts,npft)
INTEGER, INTENT(IN) :: ThreeLayerSnowFlag(land_pts, nsurft)
INTEGER, INTENT(IN) :: ICE_SurfaceType !CABLE surface tile PFT/nveg
INTEGER, INTENT(IN)  :: SurfaceType(mp) ! surface tile PFT/nveg
INTEGER, INTENT(IN)  :: SoilType(mp)    ! soil type per tile
REAL,    INTENT(OUT) :: rad_trad(mp)
REAL,    INTENT(OUT) :: met_tk(mp)

TYPE(canopy_type),         INTENT(OUT) :: canopy
TYPE(veg_parameter_type),  INTENT(OUT) :: veg        ! vegetation parameters
TYPE(soil_parameter_type), INTENT(OUT) :: soil       ! soil parameters
TYPE(soil_snow_type),      INTENT(OUT) :: ssnow      ! 
TYPE(params_io_data_type), INTENT(IN)  :: pars
TYPE(bgc_pool_type),       INTENT(OUT) :: bgc
TYPE(sum_flux_type),       INTENT(OUT) :: sum_flux

! only needed to set rad%otrad on the first timestep. 
canopy%ga      = 0.0
canopy%fes_cor = 0.0
canopy%fhs_cor = 0.0
canopy%us      = 0.01
canopy%fwsoil  = 1.0

CALL initialize_veg( SurfaceType, SoilType, mp, ms,  &
                     nrb, npft, nsurft, land_pts, l_tile_pts, ICE_SurfaceType, &
                     ICE_SoilType, LAI_thresh, soil_zse, surft_pts,            & 
                     surft_index, tile_frac,                                   &
                     veg, soil, pars )
                   

CALL initialize_soil( nsurft, land_pts, ms, mp, nsoil_max, ICE_soiltype,       &
                      surft_pts, surft_index, L_tile_pts, soiltype,            &
                      bexp, hcon, satcon, sathh, smvcst, smvcwt,               &
                      smvccl, albsoil, soil_zse, pars, soil )

! def met_tk on first call as is needed in wetfac calc ( which should be 
! replaced with tss) approximating surface temp of lake and to init trad below
CALL cable_pack_rr( met_tk, tl_1, mp, l_tile_pts, row_length, rows, nsurft,    &
                    land_pts, land_index, surft_pts, surft_index )

CALL initialize_soilsnow( mp, msn, ms, TFRZ, land_pts, nsurft, row_length,     &
                          rows, ICE_SoilType, l_tile_pts, surft_pts,           &
                          surft_index, smvcst, SoilTemp, FrozenSoilFrac,       &
                          SoilMoisture, snow_tile, OneLyrSnowDensity, SnowAge, &
                          ThreeLayerSnowFlag, SnowDensity, SnowDepth,          &
                          SnowMass, SnowTemp, soil, ssnow, veg%iveg, met_tk )

CALL init_bgc_vars( pars, bgc, veg ) 
CALL init_sumflux_zero( sum_flux ) 
CALL init_respiration( land_pts, nsurft, npft, L_tile_pts,                   &
                       npp_pft_acc, resp_w_pft_acc, canopy )

rad_trad = met_tk

RETURN
END SUBROUTINE init_data
                                      
END MODULE cbl_um_init_mod
