MODULE cable_explicit_main_mod
  
CONTAINS

SUBROUTINE cable_explicit_main(                                               &
            timestep_width,                                                   &
            cycleno, numcycles,                                               &
            land_pts, ntiles,                                                 &
            npft, sm_levels,                                                  &
            land_index, tile_frac, tile_pts, tile_index,                      &
            snow_tile, lw_down,                                               &
            surf_down_sw,                                                     &
            tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,          &
            Fland, co2_mmr, sthu, canht_ft, lai_ft ,                          &
            sin_theta_latitude,                                                &
            ftl_tile, fqw_tile,                                               &
            tstar_tile, u_s, u_s_std_tile, cd_tile, ch_tile,                  &
            radnet_tile, fraca, rESFS, resft, z0h_tile, z0m_tile,             &
            recip_l_mo_tile, epot_tile, doy,                                  &
            air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                 &
            ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,               &
            soilin, soil_cbl )

  !subrs called 
USE cable_explicit_driv_mod, ONLY: cable_explicit_driver
USE cable_expl_unpack_mod, ONLY: cable_expl_unpack

USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soilin_type
USE cable_params_mod,         ONLY: soil_parameter_type

!data !H!
USE cable_other_constants_mod, ONLY: z0surf_min

  !processor number, timestep number / width, endstep
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY: cable_runtime
!jhan:this looks like I was testing standalone

!JULES5.3
!C!USE parallel_mod,   ONLY : mype => task_id 
USE model_time_mod, ONLY: rtimestep => timestep

USE atm_fields_bounds_mod, ONLY: tdims

USE jules_soil_mod, ONLY: dzsoil

USE model_grid_mod, ONLY: latitude, longitude

USE p_s_parms,      ONLY: bexp   => bexp_soilt,                               &
                           hcon   => hcon_soilt,                              &
                           satcon => satcon_soilt,                            &
                           sathh  => sathh_soilt,                             &
                           smvcst => smvcst_soilt,                            &
                           smvcwt => smvcwt_soilt,                            &
                           smvccl => smvccl_soilt,                            &
                           albsoil =>albsoil_soilt,                           &
                           cosine_zenith_angle =>  cosz_ij 

!Imports for driving and flux variables
USE forcing, ONLY: ls_rain => ls_rain_ij,                                     &
                    ls_snow => ls_snow_ij

!cable progs are set here
USE cable_prognostic_info_mod, ONLY:                                          &
  SoilTemp        =>  SoilTemp_CABLE,                                         &
  SoilMoisture    =>  SoilMoisture_CABLE,                                     &
  FrozenSoilFrac  => FrozenSoilFrac_CABLE,                                    &
  SnowDepth   => SnowDepth_CABLE,                                             &
  SnowMass    => SnowMass_CABLE,                                              &
  SnowTemp    => SnowTemp_CABLE,                                              &
  SnowDensity => SnowDensity_CABLE,                                           &
  SnowAge     => SnowAge_CABLE,                                               &
  ThreeLayerSnowFlag  => ThreeLayerSnowFlag_CABLE,                            &
  OneLyrSnowDensity   => OneLyrSnowDensity_CABLE

!data
USE cable_types_mod, ONLY: L_tile_pts
USE cable_types_mod, ONLY: mp
USE cable_other_constants_mod, ONLY: nrb
IMPLICIT NONE
 
!___ re-decl input args
INTEGER :: endstep, cycleno, numcycles
REAL :: timestep_width
INTEGER :: row_length, rows

INTEGER ::                                                                    &
  land_pts,         & ! # of land points being processed
  ntiles,           & ! # of tiles 
  npft,             & ! # of plant functional types
  sm_levels          ! # of soil layers 

# if defined(UM_JULES)
REAL,  DIMENSION(tdims%i_end,tdims%j_end) :: latitude,  longitude
# endif   
INTEGER, DIMENSION(land_pts) :: land_index  ! index of land points processed

INTEGER,  DIMENSION(ntiles) :: tile_pts ! # of land points on each tile

INTEGER,  DIMENSION(land_pts, ntiles) :: tile_index ! index of tile points 

REAL, DIMENSION(land_pts, ntiles) :: tile_frac
   
REAL,  DIMENSION(land_pts, ntiles) :: snow_tile
   
REAL,  DIMENSION(tdims%i_end,tdims%j_end) ::                                  &
lw_down!,               &

REAL,  DIMENSION(land_pts, ntiles) :: surf_down_sw 
  
REAL,  DIMENSION(tdims%i_end,tdims%j_end) ::                                  &
  tl_1,                                                                       &
  qw_1,                                                                       &
  vshr_land,                                                                  &
  pstar,                                                                      &
  z1_tq,                                                                      &
  z1_uv

REAL, DIMENSION(land_pts, ntiles) :: canopy_tile
  
REAL, DIMENSION(land_pts) :: fland
  
REAL :: co2_mmr

REAL, DIMENSION(land_pts, sm_levels) :: sthu 
   
REAL, DIMENSION(land_pts, npft) :: canht_ft, lai_ft 
  
REAL :: sin_theta_latitude(tdims%i_end,tdims%j_end) 
    
!___return miscelaneous 
REAL, DIMENSION(land_pts,ntiles) ::                                           &
  ftl_tile,    & ! Surface FTL for land tiles     
  fqw_tile,    & ! Surface FQW for land tiles     
  tstar_tile,  & ! radiative temperature of surface 
  u_s,         & ! Surface friction velocity (m/s)
  u_s_std_tile,& ! Surface friction velocity
  cd_tile,     & ! Drag coefficient
  ch_tile,     &  ! Transfer coefficient for heat & moisture
  radnet_tile, & ! Surface net radiation
  fraca,       & ! Fraction of surface moisture
  resfs,       & ! Combined soil, stomatal & aerodynamic resistance
                 ! factor for fraction (1-FRACA) of snow-free land tiles
  resft,       & ! Total resistance factor.
                 ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                 ! 1 for snow.    
  z0h_tile,                                                                   &
  z0m_tile,                                                                   &
  recip_l_mo_tile,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
  epot_tile

INTEGER :: timestep_number  
INTEGER :: doy

TYPE(air_type),       INTENT(INOUT)  :: air_cbl
TYPE(met_type),       INTENT(INOUT)  :: met_cbl
TYPE(radiation_type),       INTENT(INOUT)  :: rad_cbl
TYPE(roughness_type),     INTENT(INOUT)  :: rough_cbl
TYPE(canopy_type),    INTENT(INOUT)  :: canopy_cbl
TYPE(soil_snow_type),     INTENT(INOUT)  :: ssnow_cbl
TYPE(bgc_pool_type),       INTENT(INOUT)  :: bgc_cbl
TYPE(balances_type),       INTENT(INOUT)  :: bal_cbl
TYPE(sum_flux_type),  INTENT(INOUT)  :: sum_flux_cbl
TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(soilin_type),  INTENT(INOUT) ::  soilin  
TYPE(soil_parameter_type),  INTENT(INOUT) ::  soil_cbl

  !___ local vars
INTEGER,  DIMENSION(land_pts, ntiles) :: iThreeLayerSnowFlag
LOGICAL, SAVE :: first_call = .TRUE.
REAL :: radians_degrees

! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_main"

row_length = tdims%i_end
rows = tdims%j_end
  
!--- initialize cable_runtime% switches 
cable_runtime%um =          .false.
cable_runtime%um_explicit = .TRUE.
   
! initialize processor number, timestep width & number, endstep 
! UM overwrites these defaults. Not yet done in StandAlone 
IF ( first_call ) THEN
  knode_gl  = 0; kend_gl=-1
# if defined(UM_JULES)
      !knode_gl  = mype 
  kwidth_gl = INT(timestep_width)
  kend_gl   = endstep   
# endif
!see Ticket #252
    sin_theta_latitude  = sin(latitude)

      !--- Convert lat/long to degrees
  radians_degrees = 180.0 / ( 4.0 * ATAN(1.0) ) ! 180 / PI
  latitude  = latitude * radians_degrees
  longitude  = longitude * radians_degrees
END IF

timestep_number = INT(rtimestep)
ktau_gl   = timestep_number
kwidth_gl = timestep_width
!----------------------------------------------------------------------------
!--- CALL _driver to run specific and necessary components of CABLE with IN -
!--- args PACKED to force CABLE
!----------------------------------------------------------------------------
iThreeLayerSnowFlag= INT(ThreeLayerSnowFlag)

CALL cable_explicit_driver(                                                   &
!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

!Variables to be calculated and returned by CABLE
!------------------------------------------------------------------------------
ftl_tile, & !surface sensible heat flux  [W/m2]? (up=+ve?) -"FTL_tile" in CABLE
fqw_tile, & !surface moisture flux flux  [kg/m^2/s/m2](up=+ve?) units could be changed?
tstar_tile,         & !surface temperature [K] per tile
u_s,                & ! land point surface friction velocity [m/s]
u_s_std_tile,       & ! per tile surface friction velocity [m/s] canopy%us
cd_tile,                                                                      &
ch_tile,                                                                      &
radnet_tile,        & !Net radiation at surface [W/m2]
fraca,& !Fraction of surface moisture flux with only aerodynamic resistance for
                      !snow-free land tiles.
resfs,& !Combined soil, stomatal and aerodynamic resistance factor for fraction
                      !(1-FRACA) of snow-free land tiles.
resft, & !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for snow.
Z0H_tile,           & ! Tile roughness lengths for heat and moisture (m).
Z0M_tile,           & ! OUT Tile roughness lengths for momentum.
RECIP_L_MO_tile,   & ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
EPOT_tile,          & ! Potential evaporation from surface, per tile
!------------------------------------------------------------------------------

!This is an "outlier" and possibly misleading. CABLE does not actually calculate radiation:
!Generally we speak of 4-band radiation. This is actually only 2-bands VIS/NIR in the SW.
!We further split each of these into Direct Beam and Diffuse components. 
!Offline CABLE splits the total SW forcing into VIS/NIR using a Spitter() fuction
!We include this variable here to connect back to JULES toplevel routines because:
!Online the UM radiation scheme DOES compute surf_down_sw using a more sophisticated model 
!than that which we use Offline, however not until AFTER the surface albedos have been 
!calculated which IS what is done AND here and technically does not require knowledge of 
!the downward SW. JULES aggregates this SW and threads this to the LSM as it is called 
!explicitly. 
!------------------------------------------------------------------------------
surf_down_sw,     & ! ShortWave radiation per rad band (row_length,rows,4) 
!------------------------------------------------------------------------------

!Mostly model dimensions and associated
!------------------------------------------------------------------------------
row_length,         & !grid cell x
rows,               & !grid cell y
land_pts,           & !grid cell land points on the x,y grid
ntiles,             & !grid cell number of surface types [] 
sm_levels,          & !grid cell number of soil levels 
npft,               & !grid cell number of PFTs 
tile_pts,           & !Number of land points per PFT [] 
tile_index,         & !Index of land point in (land_pts) array[] 
land_index, & !Index of land points in (x,y) array - see  corresponding *decs.inc
timestep_width,     & !bin width in seconds of timestep
endstep,            & !last timestep of experiment
timestep_number,                                                              &
doy,                                                                          &
mp,                                                                           &
nrb,                                                                          &
!------------------------------------------------------------------------------

!Surface descriptions generally parametrized
!------------------------------------------------------------------------------
Fland,              & !fraction of land per land point (could be coastal) 
tile_frac,         & !fraction of each surface type per land point [frac_surft]
L_tile_pts,       & !Logical mask TRUE where tile frac > 0. used to PACK/UNPACK
LAI_ft,             & !Leaf area index. [LAI_pft/LAI_pft_um in radiation]
canht_ft,           & !Canopy height [canht_pft/HGT_pft_um in radiation]
albsoil,    & !(albsoil)Snow-free, bare soil albedo [albsoil_soilt(:,1) in um ]
z0surf_min,                                                                   &
dzsoil,             & !soil thicknesses in each layer  
bexp,                                                                         &
hcon,                                                                         &
satcon,                                                                       &
sathh,                                                                        &
smvcst,                                                                       &
smvcwt,                                                                       &
smvccl,                                                                       &
!------------------------------------------------------------------------------

!Variables passed from JULES/UM
!------------------------------------------------------------------------------
latitude,           & !latitude
longitude,          & !longitude
cosine_zenith_angle,& ! cosine_zenith_angle [cosz_ij]
sin_theta_latitude,                                                           &
sthu,                                                                         &
lw_down,                                                                      &
ls_rain,                                                                      &
ls_snow,                                                                      &
tl_1,                                                                         &
qw_1,                                                                         &
vshr_land,                                                                    &
pstar,                                                                        &
z1_tq,                                                                        &
z1_uv,                                                                        &
snow_tile,          & !snow depth equivalent (in water?) [snow_surft]
                      !This is the total snow depth per tile. CABLE also has depth per layer
canopy_tile,                                                                  &
co2_mmr,                                                                      &
!------------------------------------------------------------------------------

!CABLE prognostics
!------------------------------------------------------------------------------
iThreeLayerSnowFlag, & ! flag indicating whether enough snow to treat as 3 layers
                      ! [real(ThreeLayerSnowFlag_CABLE)]
OneLyrSnowDensity, & ! density considering snow as 1 layer [OneLyrSnowDensity_CABLE
SnowAge,                                                                      &
SnowDensity,                                                                  &
SnowDepth,                                                                    &
SnowMass,                                                                     &
SnowTemp,                                                                     &
SoilMoisture,                                                                 &
FrozenSoilFrac,                                                               &
SoilTemp,                                                                     &
air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,                             &
ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl, veg_cbl,                           &
soilin, soil_cbl )
  
  !----------------------------------------------------------------------------
  !--- CALL _unpack to unpack variables from CABLE back to UM format to return
  !----------------------------------------------------------------------------
CALL cable_expl_unpack( latitude, longitude, ftl_tile, fqw_tile, tstar_tile,  &
                        u_s, u_s_std_tile,                                    &
                        cd_tile, ch_tile, fland, radnet_tile,                 &
                        fraca, rESFS, resft, z0h_tile, z0m_tile,              &
                        recip_l_mo_tile, epot_tile, l_tile_pts,               &
                        ssnow_cbl%snowd, ssnow_cbl%cls, air_cbl%rlam, air_cbl%rho, &
                        canopy_cbl%fe, canopy_cbl%fh, canopy_cbl%us, canopy_cbl%cdtq, &
                        canopy_cbl%fwet, canopy_cbl%wetfac_cs, canopy_cbl%rnet, &
                        canopy_cbl%zetar, canopy_cbl%epot, met_cbl%ua, rad_cbl%trad, &
                        rad_cbl%transd, rough_cbl%z0m, rough_cbl%zref_tq,     &
                        canopy_cbl%fes, canopy_cbl%fev )

cable_runtime%um_explicit = .FALSE.
first_call = .FALSE.        

RETURN

END SUBROUTINE cable_explicit_main
  
END MODULE cable_explicit_main_mod

