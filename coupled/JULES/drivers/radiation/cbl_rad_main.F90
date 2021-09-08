MODULE cable_rad_main_mod
  
CONTAINS

SUBROUTINE cable_rad_main(                                                    &
!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

!Variables to be calculated and returned by CABLE
!------------------------------------------------------------------------------
land_albedo, & ! GridBoxMean albedo per rad band (row_length,rows,4) [land_albedo_ij]
alb_surft,    & ! albedo per rad band per tile (land_pts, nsurft, 4) [alb_tile]
!------------------------------------------------------------------------------

!Mostly model dimensions and associated
row_length,          & !grid cell x
rows,                & !grid cell y
land_pts,            & !grid cell land points on the x,y grid
nsurft,              & !grid cell number of surface types [nsurft] 
!sm_levels,           & !grid cell number of soil levels 
npft,                & !grid cell number of PFTs 
tile_pts,            & !Number of land points per PFT [surft_pts] 
tile_index,          & !Index of land point in (land_pts) array[surft_index] 
land_index, & !Index of land points in (x,y) array - see  corresponding *decs.inc
!------------------------------------------------------------------------------

!Surface descriptions generally parametrized
!------------------------------------------------------------------------------
!dzsoil,              & !soil thicknesses in each layer  
tile_frac,         & !fraction of each surface type per land point [frac_surft]
LAI_pft_um,          & !Leaf area index. [LAI_pft]
HGT_pft_um,          & !Canopy height [canht_pft]
soil_alb,          & !(albsoil)Snow-free, bare soil albedo [albsoil_soilt(:,1)]
!------------------------------------------------------------------------------

!Variables passed from JULES/UM
!------------------------------------------------------------------------------
!This is the total snow depth per tile. CABLE also has depth per layer
SnowOtile,           & !snow depth equivalent (in water?) Prev. dt 
snow_tile,           & !snow depth equivalent (in water?) [snow_surft]
cosine_zenith_angle, & ! cosine_zenith_angle [cosz_ij]
!------------------------------------------------------------------------------
  !The args below are passed from control() level as they usually do not exist
  !in the JULES rasiation pathway -------------------------------------------------
!Mostly model dimensions and associated!---------------------------------------
sm_levels,           & !grid cell number of soil levels 

!Surface descriptions generally parametrized!----------------------------------
dzsoil,              & !soil thicknesses in each layer  

!CABLE dimensionsi/masks !------------------------------------------------------------
mp,             &!# CABLE vars assume vector of length mp(=Active patch)
msn,          &!# of snow layers. at present=3 
nrb,            &!# of rad. bands(confused b/n VIS/NIR, dir/dif. wrongly=3
jls_standalone, &! Am I in standalone mode 
L_tile_pts,     &!Logical mask. TRUE where tile frac > 0. else = FALSE
veg_mask, sunlit_mask, sunlit_veg_mask, & !Logical masks. TRUE where veg/sunlit/Both
!introduced prognostics. tiled soil on 6 layers. tiled snow on 3 layers etc!---
SoilTemp_CABLE,           &!soil temperature (IN for rad.)
SnowTemp_CABLE,           &!snow temperature (IN for rad.) REDUNDANT
ThreeLayerSnowFlag_CABLE, &!flag signalling 3 layer treatment (binary) IN only
OneLyrSnowDensity_CABLE,                                                      &

!constants!--------------------------------------------------------------------
Cz0surf_min,                                                                  &
Clai_thresh,                                                                  &
Ccoszen_tols,                                                                 &
Cgauss_w,                                                                     &
Cpi,                                                                          &
Cpi180,                                                                       &

!Vegetation parameters!--------------------------------------------------------
SurfaceType,                                                                  &
VegXfang,                                                                     &
VegTaul,                                                                      &
VegRefl,                                                                      &
!Soil parameters!--------------------------------------------------------
SoiliSoilm                                                                    &
) 
 
!subrs
USE cable_rad_driv_mod,         ONLY: cable_rad_driver
USE cable_rad_unpack_mod,       ONLY: cable_rad_unpack
USE cable_pack_mod,             ONLY: cable_pack_lp !packing
USE cable_pack_mod,             ONLY: cable_pack_rr !packing
USE cbl_LAI_canopy_height_mod,  ONLY: limit_HGT_LAI

!eliminate using data we need only 4 of these - can at least veg(soil)in%get through JaC 
USE cable_common_module,    ONLY: cable_runtime

IMPLICIT NONE

!--- IN ARGS FROM surf_couple_radiation() ------------------------------------
INTEGER :: row_length                       !grid cell x
INTEGER :: rows                             !grid cell y
INTEGER :: land_pts                         !grid cell land points on the x,y grid
INTEGER :: nsurft                           !grid cell number of surface types 
INTEGER :: npft                             !grid cell number of PFTs 

!Variables to be calculated and returned by CABLE
REAL :: land_albedo(row_length,rows,4)      
REAL :: alb_surft(Land_pts,nsurft,4)        

REAL :: tile_frac(land_pts,nsurft)          !fraction of each surface type per land point 
INTEGER :: tile_pts(nsurft)                 !Number of land points per PFT 
INTEGER :: tile_index(land_pts,nsurft)      !Index of land point in (land_pts) array
INTEGER :: land_index(land_pts)             !Index of land points in (x,y) array - see below
REAL :: LAI_pft_um(land_pts, npft)          !Leaf area index.
REAL :: HGT_pft_um(land_pts, npft)          !Canopy height
REAL :: fHGT_pft_um(land_pts, npft)          !Canopy height

REAL :: snow_tile(land_pts,nsurft)          ! snow depth equivalent (in water?)
REAL :: SnowOTile(land_pts,nsurft)         ! snow depth equivalent from previous timestep
REAL :: soil_alb(land_pts)                  !(albsoil)Snow-free, bare soil albedo
REAL :: cosine_zenith_angle(row_length,rows)             !cosine_zenith_angle          

!--- IN ARGS FROM control() ------------------------------------
INTEGER :: sm_levels
REAL :: dzsoil(sm_levels)
INTEGER :: mp
INTEGER :: msn
INTEGER :: nrb
LOGICAL :: L_tile_pts(land_pts,nsurft)
REAL :: SoilTemp_CABLE(land_pts, nsurft, sm_levels )
REAL :: SnowTemp_CABLE(land_pts, nsurft, msn)
REAL :: ThreeLayerSnowFlag_CABLE(land_pts, nsurft )
REAL :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
!constants
REAL :: Cz0surf_min                      !the minimum roughness of bare soil
REAL :: Clai_thresh                     !The minimum LAI below which a "cell" is considred NOT vegetated
REAL :: Ccoszen_tols                    !threshold cosine of sun's zenith angle, below which considered SUNLIT
REAL :: Cgauss_w(nrb)               !Gaussian integration weights
REAL :: Cpi                             !PI - describing the ratio of circumference to diameter
REAL :: Cpi180                          !PI in radians
INTEGER:: SurfaceType(mp) 
REAL :: VegXfang(mp)
REAL :: VegTaul(mp, nrb)
REAL :: VegRefl(mp, nrb)
INTEGER :: SoiliSoilm(mp)

INTEGER :: metDoY(mp)          !local dummy Day of the Year [formerly met%doy]

!Convoluted mapping using land_index(array) to get back to the row_length,rows co-ordinate
! J = ( LAND_INDEX(L)-1 ) / ROW_LENGTH + 1
! I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
! FTL_1(I,J) = where ftl_1(row_length,rows) 
!-------------------------------------------------------------------------------
!these local to CABLE and can be flushed every timestep
REAL :: ExtCoeff_beam(mp)
REAL :: ExtCoeff_dif(mp)
REAL :: EffExtCoeff_beam(mp, nrb)
REAL :: EffExtCoeff_dif(mp, nrb)

REAL :: CanopyTransmit_dif(mp, nrb)
REAL :: CanopyTransmit_beam(mp, nrb)
REAL :: CanopyRefl_dif(mp, nrb)
REAL :: CanopyRefl_beam(mp, nrb)

REAL :: EffSurfRefl_dif(mp, nrb)
REAL :: EffSurfRefl_beam(mp, nrb)

REAL :: coszen(mp)
REAL :: RadFbeam(mp, nrb)
REAL :: RadAlbedo(mp, nrb)
REAL :: AlbSnow(mp, nrb)
REAL :: AlbSoil(mp, nrb)

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp, nrb)
REAL :: rhoch(mp, nrb)
REAL :: xk(mp, nrb)
!-------------------------------------------------------------------------------
!highlight need to allocate these early on and thres to here as SAVE is not an
!option
REAL :: SnowDepth(mp)             !Formerly: ssnow%snowd 
REAL :: SnowODepth(mp)            !Formerly: ssnow%osnowd
REAL :: SnowDensity(mp)           !Formerly: ssnow%ssdnn 
!computed from UM HT(LAI)_PFT passed in explicit call - need at rad call
REAL :: LAI_pft_cbl(mp)           !Formerly: ~veg%vlai
REAL :: HGT_pft_cbl(mp)           !Formerly:  ~veg%hc 
!can compute from  z0surf_min, HGT_pft_cbl, SnowDepth, SnowDensity )
REAL :: HeightAboveSnow(mp)       !Formerly: rough%hruff

REAL :: MetTk(mp) 
REAL :: SoilTemp(mp)
REAL :: SnowTemp(mp)
REAL :: SnowAge(mp)
INTEGER:: SnowFlag_3L(mp)

!--- declare vars local to CABLE -------------------------------------------- 
!packed in pack
LOGICAL :: jls_standalone 
LOGICAL :: um_online = .FALSE. 
LOGICAL :: jls_radiation = .TRUE. !um_radiation = jls_radiation
  
!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
REAL :: reducedLAIdue2snow(mp)
!masks
LOGICAL, ALLOCATABLE :: veg_mask(:),  sunlit_mask(:),  sunlit_veg_mask(:) 

!--- declare local vars to subr -------------------------------------------- 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_main"
LOGICAL, SAVE :: first_call = .TRUE.
INTEGER :: i, j, n
LOGICAL :: skip =.TRUE. 
REAL :: SW_down(mp,2)
LOGICAL, SAVE :: albflip = .FALSE.
REAL, SAVE :: ialb_surft
!from rose-app.conf!!canht_ft_io= !!lai_io=
!LAI_pft_um(1, :) = (/4.0,5.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/) 
!HGT_pft_um(1,:) = (/16.38,19.01,0.0,0.0,0.0,0.79,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
!HGT_pft_um(1,14) = 20.0
jls_radiation= .TRUE.

!--- initialize/zero each timestep 
CALL zero_albedo( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,     &
                  EffExtCoeff_dif, CanopyTransmit_dif, CanopyTransmit_beam,   &
                  CanopyRefl_dif,CanopyRefl_beam, EffSurfRefl_dif,            &
                  EffSurfRefl_beam, coszen, RadFbeam, RadAlbedo, AlbSnow,     &
                  AlbSoil, c1, rhoch, xk )

metDoY = 0    !can pass DoY from current_time%
!jhan:ubiquitous confusion between nrb, dir/dif radiation

!JULES does not have seperate albedos per nrb so assume same value here
albsoil(:,3)  = 0.0
CALL cable_pack_lp( soil_alb, soil_alb, albsoil(:,1), SoiliSoilm, skip )
CALL cable_pack_lp( soil_alb, soil_alb, albsoil(:,2), SoiliSoilm, skip )

!Pack zenith this timestep
CALL cable_pack_rr( cosine_zenith_angle, coszen)

!Store Snow Depth from previous timestep. Treat differently on 1st timestep 
SnowODepth = PACK( snowotile, l_tile_pts )
SnowDepth = PACK( snow_tile, l_tile_pts )
SnowOTile = UNPACK( SnowDepth, l_tile_pts, 0.0 )

SnowDensity = PACK( OneLyrSnowDensity_CABLE, l_tile_pts )

!Treat snow depth across 3Layers? Typecasts from Real to integer
SnowFlag_3L = PACK( ThreeLayerSnowFlag_CABLE, l_tile_pts )

!Surface skin/top layer Soil/Snow temperature
SoilTemp =   PACK( SoilTemp_cable(:,:,1), l_tile_pts )
SnowTemp =   PACK( SnowTemp_cable(:,:,1), l_tile_pts )

! limit IN height, LAI  and initialize existing cable % types
CALL limit_HGT_LAI( fHGT_pft_um, LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, nsurft,           &
                    tile_pts, tile_index, tile_frac, L_tile_pts,              &
                    LAI_pft_um, HGT_pft_um, CLAI_thresh )

!------------------------------------------------------------------------------
! Call CABLE_rad_driver to run specific and necessary components of CABLE 
!------------------------------------------------------------------------------
CALL cable_rad_driver( mp, nrb, Clai_thresh, Ccoszen_tols,  CGauss_w, Cpi,    &
                       veg_mask, sunlit_mask, sunlit_veg_mask, Cpi180,        &
                       jls_standalone, jls_radiation, Cz0surf_min, SurfaceType, &
                       LAI_pft_cbl, HGT_pft_cbl, SnowDepth, SnowODepth,       &
                       SnowFlag_3L, SnowDensity, SoilTemp, SnowTemp, SnowAge, &
                       HeightAboveSnow, ExtCoeff_beam, ExtCoeff_dif,          &
                       EffExtCoeff_beam, EffExtCoeff_dif,                     &
                       CanopyRefl_dif, CanopyRefl_beam,                       &
                       CanopyTransmit_dif, CanopyTransmit_beam,               &
                       EffSurfRefl_dif, EffSurfRefl_beam,                     &
                       RadAlbedo, reducedLAIdue2snow,                         &
                       coszen, VegXfang, VegTaul, VegRefl,                    &
                       c1, rhoch, metDoY, SW_down,                            &
                       RadFbeam, xk, AlbSnow, AlbSoil, metTk )

! Unpack variables (CABLE computed albedos) to JULES 
!------------------------------------------------------------------------------
CALL cable_rad_unpack( land_albedo, alb_surft, mp, nrb, row_length, rows,     &
                       land_pts, nsurft, sm_levels, tile_pts, tile_index,     &
                       land_index, tile_frac, L_tile_pts,                     &
                       EffSurfRefl_dif, EffSurfRefl_beam,                     &
                       ThreeLayerSnowFlag_CABLE, SnowFlag_3L )

!flick switches before leaving  
jls_radiation= .FALSE.
first_call = .FALSE.

RETURN

END SUBROUTINE cable_rad_main

SUBROUTINE zero_albedo( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam, &
                        EffExtCoeff_dif, CanopyTransmit_dif,                  &
                        CanopyTransmit_beam, CanopyRefl_dif, CanopyRefl_beam, &
                        EffSurfRefl_dif, EffSurfRefl_beam,                    &
                        coszen, RadFbeam, RadAlbedo, AlbSnow, AlbSoil, c1,    &
                        rhoch, xk )

INTEGER :: mp, nrb
!these local to CABLE and can be flushed every timestep
REAL :: ExtCoeff_beam(mp)
REAL :: ExtCoeff_dif(mp)
REAL :: EffExtCoeff_beam(mp, nrb)
REAL :: EffExtCoeff_dif(mp, nrb)

REAL :: CanopyTransmit_dif(mp, nrb)
REAL :: CanopyTransmit_beam(mp, nrb)
REAL :: CanopyRefl_dif(mp, nrb)
REAL :: CanopyRefl_beam(mp, nrb)

REAL :: EffSurfRefl_dif(mp, nrb)
REAL :: EffSurfRefl_beam(mp, nrb)

REAL :: coszen(mp)
REAL :: RadFbeam(mp, nrb)
REAL :: RadAlbedo(mp, nrb)
REAL :: AlbSnow(mp, nrb)
REAL :: AlbSoil(mp, nrb)

REAL :: c1(mp, nrb)
REAL :: rhoch(mp, nrb)
REAL :: xk(mp, nrb)

ExtCoeff_beam(:) = 0.0
ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0
EffExtCoeff_dif(:,:) = 0.0
CanopyTransmit_dif(:,:) = 0.0
CanopyTransmit_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0
CanopyRefl_beam(:,:) = 0.0
EffSurfRefl_dif(:,:) = 0.0
EffSurfRefl_beam(:,:) = 0.0
coszen(:) = 0.0
RadFbeam(:,:) = 0.0
RadAlbedo(:,:) = 0.0
AlbSnow(:,:) = 0.0
AlbSoil(:,:) = 0.0
c1(:,:) = 0.0
rhoch(:,:) = 0.0
xk(:,:) = 0.0

END SUBROUTINE zero_albedo


END MODULE cable_rad_main_mod

