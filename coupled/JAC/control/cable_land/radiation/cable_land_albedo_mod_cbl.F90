!******************************COPYRIGHT********************************************
! (c) CSIRO 2022.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT********************************************

MODULE cable_land_albedo_mod

!-----------------------------------------------------------------------------
! Description:
!  Computes the albedo for CABLE and returns it to JULES
!
! This MODULE is USEd in:
!     surf_couple_radiation_mod.F90 (JULES)
!
! This MODULE contains 1 public Subroutine:
!     cable_land_albedo
! Other Subroutines:
!     cable_pack_Albsoil,
!     cable_pack_progs
!
! Code owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: cable_land_albedo
PRIVATE

CONTAINS

SUBROUTINE cable_land_albedo (                                                 &
  !OUT: JULES (per rad band) albedos [GridBoxMean & per tile albedo]
  land_albedo , alb_surft,                                                     &
  !IN: JULES dimensions and associated
  row_length, rows, land_pts, nsurft, npft,                                    &
  surft_pts, surft_index, land_index,                                          &
  !IN: JULES Surface descriptions generally parametrized
  tile_frac, LAI_pft_um, HGT_pft_um, soil_alb,                                 &
  !IN: JULES  timestep varying fields
  cosine_zenith_angle, snow_tile,                                              &
  !IN:CABLE dimensions from grid_constants_cbl
  nrb, nrs, mp,                                                                &
  !IN: CABLE specific surface_type indexes
  ICE_Surfacetype, lakes_SurfaceType, ICE_SoilType,                            &
  !IN: CABLE constants
  Cz0surf_min, Clai_thresh, Ccoszen_tols, Cgauss_w, Cpi, Cpi180,               &
  !IN: CABLE Vegetation/Soil parameters. decl in params_io_cbl.F90
  VeginXfang, VeginTaul, VeginRefl,                                            &
  !IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90
  SoilTemp_CABLE, OneLyrSnowDensity_CABLE, SnowAge_CABLE                       &
)
!-------------------------------------------------------------------------------
! Description:
!   Provide (return) albedo(s) to JULES [land_albedo , alb_surft]
!   per rad stream (VIS/NIR, Direct&Diffuse) [GridBoxMean & per tile albedo]
! Three main sections:
!   1. Pack CABLE variables from those passed from surf_couple_radiation()
!   2. Call CABLE's radiation/albedo scheme
!   3. Unpack albedos to send back to JULES
!-------------------------------------------------------------------------------

! USE subroutines
! Route into CABLE & UNPACKING what we have computed for JULES
USE cable_rad_driv_mod,         ONLY: cable_rad_driver
USE cable_rad_unpack_mod,       ONLY: cable_rad_unpack

! PACKING from JULES array dims to CABLE active tile 1-D vector
USE init_active_tile_mask_mod,  ONLY: init_active_tile_mask_cbl
USE cable_pack_mod,             ONLY: cable_pack_rr

! Define CABLE grid, sunlit/veg masks & initialize surface type params
USE init_cable_parms_mod,       ONLY: init_cable_parms_rad
USE alloc_rad_albedo_vars_mod,  ONLY: alloc_local_vars, flush_local_vars
USE cbl_masks_mod,              ONLY: fveg_mask, L_tile_pts

!Compute canopy exposed above (potential) snow
USE cbl_LAI_canopy_height_mod,  ONLY: limit_HGT_LAI
USE hruff_eff_LAI_mod_cbl, ONLY: HgtAboveSnow, LAI_eff

IMPLICIT NONE
! re-decl dims necessary to declare OUT fields
!-- IN: JULES model dimensions
INTEGER, INTENT(IN) :: row_length, rows        !# grid cell x, y
INTEGER, INTENT(IN) :: nsurft                  !# surface types, PFTS
INTEGER, INTENT(IN) :: land_pts                !# land points on x,y grid
!--- IN: CABLE  declared in grid_cell_constants_cbl
INTEGER, INTENT(IN) :: nrs                     !# rad streams
                                               !(:,:,1) direct beam VIS
                                               !(:,:,2) diffuse visible
                                               !(:,:,3) direct beam NIR
                                               !(:,:,4) diffuse NIR
!--- OUT: JULES (per rad band) albedos [GridBoxMean & per tile albedo]
REAL,    INTENT(OUT) :: land_albedo(row_length,rows,nrs)    ! [land_albedo_ij]
REAL,    INTENT(OUT) :: alb_surft(Land_pts,nsurft,nrs)      ! [alb_tile]

!-- IN: JULES model dimensions
INTEGER, INTENT(IN) :: npft                    !# surface types, PFTS

!---IN: JULES model associated
INTEGER, INTENT(IN) :: surft_pts(nsurft)            ! # land points per PFT
INTEGER, INTENT(IN) :: surft_index(land_pts,nsurft) ! Index in land_pts array
INTEGER, INTENT(IN) :: land_index(land_pts)         ! Index in (x,y) array

!-- IN: JULES Surface descriptions generally parametrized
REAL, INTENT(IN) :: tile_frac(land_pts,nsurft)      ! fraction of each surf type
REAL, INTENT(IN) :: LAI_pft_um(land_pts, npft)      ! Leaf area index.
REAL, INTENT(IN) :: HGT_pft_um(land_pts, npft)      ! Canopy height
REAL, INTENT(IN) :: soil_alb(land_pts)              ! Snow-free, soil albedo

!---IN: JULES  timestep varying fields
REAL, INTENT(IN) :: cosine_zenith_angle(row_length,rows)  ! zenith angle of sun
REAL, INTENT(IN) :: snow_tile(land_pts,nsurft)            ! snow depth (units?)

!--- IN: CABLE  declared in grid_cell_constants_cbl
INTEGER, INTENT(IN) :: nrb                     !# rad bands VIS/NIR + Legacy LW
INTEGER, INTENT(OUT) :: mp        ! curr. NOT requ'd OUT, however it likely will

!--- IN: CABLE specific Surface/Soil type indexes
INTEGER, INTENT(IN) :: ICE_SurfaceType
INTEGER, INTENT(IN) :: lakes_SurfaceType
INTEGER, INTENT(IN) :: ICE_SoilType

!---IN: CABLE constants
REAL, INTENT(IN) :: Cz0surf_min     ! the minimum roughness of bare soil
REAL, INTENT(IN) :: Clai_thresh     ! min. LAI signalling a cell is vegetated
REAL, INTENT(IN) :: Cgauss_w(nrb)   ! Gaussian integration weights
REAL, INTENT(IN) :: Cpi             ! PI
REAL, INTENT(IN) :: Cpi180          ! PI in radians
REAL, INTENT(IN) :: Ccoszen_tols    ! sun rise/set threshold for zenith angle
                                    ! signals daylit
!---IN: CABLE Vegetation/Soil parameters. decl in params_io_cbl.F90
REAL, INTENT(IN) :: VeginXfang(nsurft)               ! Leaf Angle
REAL, INTENT(IN) :: VeginTaul(nrb, nsurft )          ! Leaf Transmisivity
REAL, INTENT(IN) :: VeginRefl(nrb, nsurft )          ! Leaf Reflectivity

!---IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90
REAL, INTENT(IN) :: SoilTemp_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: SnowAge_CABLE(land_pts, nsurft )

!--- local vars - Neither IN nor OUT (passed to subrs)

! Albedos req'd by JULES - Effective Surface Relectance as seen by atmosphere
REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)

! vegetated mask required on albedo pathway
LOGICAL, ALLOCATABLE :: veg_mask(:)

!Vegetation/soil parameters
INTEGER, ALLOCATABLE :: SurfaceType(:)     ! veg%iveg
INTEGER, ALLOCATABLE :: SoilType(:)        ! soil%isoilm
REAL, ALLOCATABLE :: VegXfang(:)           ! Leaf Angle [veg%xfang]
REAL, ALLOCATABLE :: VegTaul(:,:)          ! Leaf Transmisivity [veg%taul]
REAL, ALLOCATABLE :: VegRefl(:,:)          ! Leaf Reflectivity [veg%refl]

! Formerly canopy%vlaiw, rough%hruff SAVEd after explicit call to CABLE
REAL, ALLOCATABLE :: reducedLAIdue2snow(:) ! Eff. LAI IF snow [canopy%vlaiw]
REAL, ALLOCATABLE :: HeightAboveSnow(:)    ! Canopy Hgt above snow (rough%hruff)

! arrays to map IN progs to CABLE vector length
REAL, ALLOCATABLE :: SnowDepth(:)     ! Total Snow depth - water eqivalent -
                                      ! ssnow%snowd
REAL, ALLOCATABLE :: SnowDensity(:)   ! Total Snow density (assumes 1 layer
REAL, ALLOCATABLE :: SoilTemp(:)      ! Soil Temperature of top layer (soil%tgg)
REAL, ALLOCATABLE :: SnowAge(:)       ! Snow age (assumes 1 layer describes snow
                                      ! ssnow%isflag

REAL, ALLOCATABLE :: AlbSoil(:,:)     ! BareSoil Albedo [soill%albsoil]
REAL, ALLOCATABLE :: coszen(:)
! corrected version of UM HT(LAI)_PFT passed in explicit call - need at rad call
REAL, ALLOCATABLE  :: LAI_pft_cbl(:)           !Formerly: ~veg%vlai
REAL, ALLOCATABLE  :: HGT_pft_cbl(:)           !Formerly:  ~veg%hc

REAL, ALLOCATABLE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, ALLOCATABLE :: RadFbeam(:,:)
REAL, ALLOCATABLE :: RadAlbedo(:,:)
REAL, ALLOCATABLE :: AlbSnow(:,:)
REAL, ALLOCATABLE :: c1(:,:)
REAL, ALLOCATABLE :: rhoch(:,:)
REAL, ALLOCATABLE :: xk(:,:)
! used in Calc of Beam calculation NOT on rad/albedo path.
! However Needed to fulfill arg list with dummy
REAL, ALLOCATABLE :: SW_down(:, :)  ! NA at surf_couple_rad layer
INTEGER, ALLOCATABLE :: metDoY(:)     ! can pass DoY from current_time

LOGICAL :: um_online = .FALSE.
LOGICAL :: jls_standalone = .TRUE.
LOGICAL :: jls_radiation = .TRUE.     !um_radiation = jls_radiation

INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_main"
! End header

! initialise INTENT(OUT) fields
land_albedo = 0.0; alb_surft = 0.0

! Determine the number of active tiles
mp = SUM(surft_pts)

! Define mapping mask. i.e. l_tile_pts =TRUE (active) , where tile_frac > 0
CALL init_active_tile_mask_cbl(l_tile_pts, land_pts, nsurft, tile_frac )

! alloc/zero each timestep
! metDoY, SW_down, RadFbeaam, RadAlbedo NOT used on rad/albedo path.
! Nevertheless,  need to fulfill later arg list(s) with dumies
CALL alloc_local_vars( EffSurfRefl_beam, EffSurfRefl_dif, mp, nrb,             &
                       reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,           &
                       HeightAboveSnow, coszen, ExtCoeff_beam,                 &
                       ExtCoeff_dif, EffExtCoeff_beam, EffExtCoeff_dif,        &
                       CanopyTransmit_beam, CanopyTransmit_dif,                &
                       CanopyRefl_beam, CanopyRefl_dif, RadFbeam,              &
                       RadAlbedo, AlbSnow, c1, rhoch, xk, metDoY,              &
                       SnowDepth, SnowDensity, SoilTemp, SnowAge,              &
                       AlbSoil, SW_down)
! -----------------------------------------------------------------------------
! 1. PACK CABLE fields
! -----------------------------------------------------------------------------
! map PFT/soil parameters to mp format
CALL init_cable_parms_rad( VegXfang, VegTaul, VegRefl, SurfaceType, SoilType,  &
                           mp, nrb, l_tile_pts, ICE_SurfaceType, ICE_SoilType, &
                           VeginXfang, VeginTaul, VeginRefl,                   &
                           land_pts, nsurft, tile_frac )
! Pack UM spatial (per landpoint) bare soil albedo to mp vector for CABLE
CALL cable_pack_Albsoil( albsoil, soil_alb, mp, nrb, l_tile_pts,               &
                         nsurft, land_pts )

! Pack UM spatial (per cell) zenith angle to mp vector for CABLE
CALL cable_pack_rr( coszen, cosine_zenith_angle, mp, l_tile_pts, row_length,   &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

! Pack UM spatial (per landpoint,ntile) CABLE prognostics to mp vector
CALL cable_pack_progs( SnowDepth, SnowDensity, SoilTemp, SnowAge, mp,          &
                       land_pts, nsurft, l_tile_pts, snow_tile,                &
                       OneLyrSnowDensity_CABLE, SoilTemp_CABLE, SnowAge_CABLE )
! -----------------------------------------------------------------------------

! limit IN height, LAI  and initialize some existing cable % types
CALL limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, nsurft, npft,      &
                    surft_pts, surft_index, tile_frac, l_tile_pts,             &
                    LAI_pft_um, HGT_pft_um, CLAI_thresh )

! set Height of Canopy above snow (rough%hruff)
CALL HgtAboveSnow( HeightAboveSnow, mp, Cz0surf_min, HGT_pft_cbl,              &
                   SnowDepth, SnowDensity )

! set Effective LAI considering potential snow coverage [cabopy%vlaiw]
CALL LAI_eff( mp, LAI_pft_cbl, HGT_pft_cbl, HeightAboveSnow, reducedLAIdue2snow)

! Define logical mask for vegetated mp cell
CALL fveg_mask( veg_mask, mp, Clai_thresh, reducedLAIdue2snow )

!------------------------------------------------------------------------------
! 2. Call CABLE_rad_driver to run specific and necessary components of CABLE
!------------------------------------------------------------------------------
CALL cable_rad_driver( EffSurfRefl_beam, EffSurfRefl_dif,                      &
                       mp, nrb, ICE_SoilType, lakes_SurfaceType, Clai_thresh,  &
                       Ccoszen_tols, CGauss_w, Cpi, Cpi180, Cz0surf_min,       &
                       veg_mask, jls_standalone, jls_radiation, SurfaceType,   &
                       SoilType, LAI_pft_cbl, HGT_pft_cbl, SnowDepth,          &
                       SnowDensity, SoilTemp, SnowAge, AlbSoil ,coszen,        &
                       VegXfang, VegTaul, VegRefl, HeightAboveSnow,            &
                       reducedLAIdue2snow, ExtCoeff_beam, ExtCoeff_dif,        &
                       EffExtCoeff_beam,  EffExtCoeff_dif, CanopyTransmit_beam,&
                       CanopyTransmit_dif, CanopyRefl_beam,CanopyRefl_dif,     &
                       c1, rhoch, xk, AlbSnow, RadFbeam, RadAlbedo, metDoY,    &
                       SW_down )

!------------------------------------------------------------------------------
! 3. Unpack variables (CABLE computed albedos) to JULES
!------------------------------------------------------------------------------
CALL cable_rad_unpack( land_albedo, alb_surft, mp, nrs, row_length, rows,      &
                       land_pts, nsurft, surft_pts, surft_index,               &
                       land_index, tile_frac, l_tile_pts,                      &
                       EffSurfRefl_dif, EffSurfRefl_beam  )

CALL flush_local_vars( EffSurfRefl_beam, EffSurfRefl_dif,SnowDepth,            &
                       SnowDensity, SoilTemp, SnowAge, AlbSoil,                &
                       reducedLAIdue2snow, LAI_pft_cbl,  HGT_pft_cbl,          &
                       HeightAboveSnow, coszen, ExtCoeff_beam, ExtCoeff_dif,   &
                       EffExtCoeff_beam, EffExtCoeff_dif, CanopyTransmit_beam, &
                       CanopyTransmit_dif, CanopyRefl_beam, CanopyRefl_dif,    &
                       RadFbeam, RadAlbedo, AlbSnow, c1, rhoch, xk, metDoY, SW_down )

!flick switches before leaving
jls_radiation= .FALSE.

RETURN

END SUBROUTINE cable_land_albedo

!==============================================================================
SUBROUTINE cable_pack_Albsoil( AlbSoil, soil_alb, mp, nrb, l_tile_pts,         &
                               nsurft, land_pts )

! Description:
!   Pack spatial UM bare soil albedo to CABLE dimensions, HOWEVER limited by
!   JULES typically has one soil type per cell & no distiction b//n rad bands

IMPLICIT NONE

REAL,    INTENT(OUT), ALLOCATABLE :: AlbSoil(:, :)
INTEGER, INTENT(IN) :: mp, nsurft, land_pts, nrb
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)
REAL,    INTENT(IN)    :: soil_alb(land_pts)
!local vars:
REAL    :: fvar(land_pts, nsurft)
INTEGER :: n,l

IF ( .NOT. ALLOCATED(AlbSoil) ) ALLOCATE (AlbSoil(mp, nrb) )

AlbSoil(:,:) = 0.0 ! albsoil(:,3)  stays =  0.0
fvar(:, :) = 0.0

DO n = 1, nsurft
  DO l = 1, land_pts
    fvar(l,n) = soil_alb(l)
  END DO
END DO

AlbSoil(:,1) = PACK(fvar, l_tile_pts)
AlbSoil(:,2) = AlbSoil(:,1)

RETURN
END SUBROUTINE cable_pack_Albsoil

SUBROUTINE cable_pack_progs( SnowDepth, SnowDensity,SoilTemp, SnowAge,         &
                             mp, land_pts, nsurft, l_tile_pts,                 &
                             snow_tile, OneLyrSnowDensity_CABLE,               &
                             SoilTemp_CABLE, SnowAge_CABLE )

! Description:
!   Pack CABLE prognostics n CABLE dimensions from passed JULES vars

IMPLICIT NONE

INTEGER, INTENT(IN)  :: land_pts, nsurft, mp

! map IN progs to CABLE veector length
REAL, INTENT(OUT), ALLOCATABLE :: SnowDepth(:)   ! Tot Snow depth - water eqiv.
REAL, INTENT(OUT), ALLOCATABLE :: SnowDensity(:) ! Snow density-assumes 1 layer
REAL, INTENT(OUT), ALLOCATABLE :: SoilTemp(:)    ! Soil Temp. of top layer
REAL, INTENT(OUT), ALLOCATABLE :: SnowAge(:)     ! Snow age (assumes 1 layer)

LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, nsurft )

!---IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90
REAL, INTENT(IN) :: SoilTemp_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: SnowAge_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: snow_tile(land_pts,nsurft)            ! snow depth (units?)

IF ( .NOT. ALLOCATED(SnowDepth) )         ALLOCATE( SnowDepth(mp) )
IF ( .NOT. ALLOCATED(SnowDensity) )       ALLOCATE( SnowDensity(mp) )
IF ( .NOT. ALLOCATED(SoilTemp) )          ALLOCATE( SoilTemp(mp) )
IF ( .NOT. ALLOCATED(SnowAge) )           ALLOCATE( SnowAge(mp) )

!Store Snow Depth from previous timestep. Treat differently on 1st timestep
SnowDepth   = PACK( snow_tile, l_tile_pts )
SnowDensity = PACK( OneLyrSnowDensity_CABLE, l_tile_pts )

!Surface skin/top layer Soil/Snow temperature
SoilTemp =   PACK( SoilTemp_CABLE(:,:), l_tile_pts )
SnowAge  =   PACK( SnowAge_CABLE(:,:), l_tile_pts )

RETURN

END SUBROUTINE cable_pack_progs

END MODULE cable_land_albedo_mod

