MODULE casa_um_inout_mod

USE landuse_mod!jh!
USE feedback_mod
USE casa_inout_module

IMPLICIT NONE 

CONTAINS 

SUBROUTINE init_casacnp( mp, land_pts, nsurft, row_length, rows, l_tile_pts,   &
                         surft_pts, surft_index, land_index, nsoil_max,        &
                         ICE_soiltype, soiltype, tile_frac, latitude,          &
                         longitude, smvcst, SoilTemp, FrozenSoilFrac,          &
                         veg, soil, canopy, Cpool_tile, Npool_tile,            &
                         Ppool_tile, soil_order, Nidep, Nifix, Pwea, Pdust,    &
                         wood_hvest_C,wood_hvest_N,wood_hvest_P, wresp_c,      &
                         wresp_n, wresp_P, thinning, gLAI, phenphase,          &
                         prev_yr_sfrac, idoy, casapool, casaflux, sum_casapool,&
                         sum_casaflux, casabiome, casamet,  casabal,  phen )

! subrs
USE casa_readbiome_module, ONLY: casa_readbiome

! data
USE grid_constants_mod_cbl, ONLY: nsl
USE cable_def_types_mod,    ONLY: canopy_type 
USE cable_def_types_mod,    ONLY: soil_parameter_type
USE cable_def_types_mod,    ONLY: veg_parameter_type

! TYPE declarations
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_biome_type_mod,   ONLY: casa_biome   => casa_biome_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE phenology_type_mod,    ONLY: phenology_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: row_length, rows
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft 
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)          ! cell index of land_pt
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft )  ! TRUE if active tile
INTEGER, INTENT(IN) :: ICE_SoilType
INTEGER, INTENT(IN) :: nsoil_max
INTEGER, INTENT(IN) :: SoilType(mp)    ! soil type per tile
REAL,    INTENT(IN) :: tile_frac(land_pts, nsurft)
REAL,    INTENT(IN) :: latitude( row_length,rows )
REAL,    INTENT(IN) :: longitude( row_length,rows )
REAL,    INTENT(IN) :: smvcst(land_pts, nsl)
REAL,    INTENT(IN) :: SoilTemp(land_pts, nsurft, nsl )
REAL,    INTENT(IN) :: FrozenSoilFrac(land_pts, nsurft, nsl )

REAL   , INTENT(INOUT) :: Cpool_tile   ( land_pts, nsurft, 10 )
REAL   , INTENT(INOUT) :: Npool_tile   ( land_pts, nsurft, 10 )
REAL   , INTENT(INOUT) :: Ppool_tile   ( land_pts, nsurft, 12 )
REAL   , INTENT(INOUT) :: soil_order   ( land_pts )
REAL   , INTENT(INOUT) :: Nidep        ( land_pts ) 
REAL   , INTENT(INOUT) :: Nifix        ( land_pts )
REAL   , INTENT(INOUT) :: Pwea         ( land_pts )
REAL   , INTENT(INOUT) :: Pdust        ( land_pts )
REAL   , INTENT(INOUT) :: glai         ( land_pts, nsurft )
REAL   , INTENT(INOUT) :: phenphase    ( land_pts, nsurft)
REAL   , INTENT(INOUT) :: prev_yr_sfrac( land_pts, nsurft )
REAL   , INTENT(INOUT) :: wood_hvest_C ( land_pts, nsurft, 3 )
REAL   , INTENT(INOUT) :: wood_hvest_N ( land_pts, nsurft, 3 )
REAL   , INTENT(INOUT) :: wood_hvest_P ( land_pts, nsurft, 3 )
REAL   , INTENT(INOUT) :: wresp_C      ( land_pts, nsurft, 3 )
REAL   , INTENT(INOUT) :: wresp_N      ( land_pts, nsurft, 3 )
REAL   , INTENT(INOUT) :: wresp_P      ( land_pts, nsurft, 3 )
REAL   , INTENT(INOUT) :: thinning     ( land_pts, nsurft )
TYPE (canopy_type),         INTENT(INOUT) :: canopy
TYPE (soil_parameter_type), INTENT(INOUT) :: soil
TYPE (veg_parameter_type),  INTENT(INOUT) :: veg

TYPE (casa_pool),          INTENT(INOUT) :: casapool
TYPE (casa_flux),          INTENT(INOUT) :: casaflux
TYPE (casa_pool),          INTENT(INOUT) :: sum_casapool
TYPE (casa_flux),          INTENT(INOUT) :: sum_casaflux
TYPE (casa_biome),         INTENT(INOUT) :: casabiome
TYPE (casa_met),           INTENT(INOUT) :: casamet
TYPE (casa_balance),       INTENT(INOUT) :: casabal
TYPE (phenology_type),      INTENT(INOUT) :: phen

INTEGER :: idoy

!jh!CALL alloc_casavariable( casabiome, casapool, casaflux, casamet, casabal, mp )

!jh!CALL alloc_phenology_data_type( phen, mp )

CALL casa_readpoint_pk( mp, land_pts, nsurft, row_length, rows, surft_pts,     &
                        surft_index, land_index, l_tile_pts, tile_frac, latitude,          &
                        longitude, veg, soil, casaflux,        &
                        casamet, nidep, nifix, pwea, pdust, soil_order )

! jh:is ntiles used here actually ntiles?
CALL casa_readbiome( mp, nsl, nsurft, veg%iveg, soil%zse, casabiome,     &
                           casapool, casaflux, casamet, phen ) 

!jh!CALL casa_readphen(veg,casamet,phen)

CALL casa_init_pk( mp, land_pts, nsurft, row_length, rows, l_tile_pts, tile_frac, &
                   casabiome, casaflux, casamet, casapool, casabal, veg,       &
                   canopy, phen, cpool_tile, npool_tile, ppool_tile,           &
                   wood_hvest_c, wood_hvest_n, wood_hvest_p, wresp_c, wresp_n, &
                   wresp_p, thinning, GLAI, PHENPHASE, PREV_YR_SFRAC, idoy )
                       
RETURN
END SUBROUTINE init_casacnp

SUBROUTINE casa_readpoint_pk( mp, land_pts, nsurft, row_length, rows, surft_pts,          &
                              surft_index, land_index, l_tile_pts, tile_frac, latitude, longitude, veg,  &
                              soil, casaflux, casamet, nidep, nifix, pwea,     &
                              pdust, soil_order )

USE cable_def_types_mod,    ONLY: soil_parameter_type
USE cable_def_types_mod,    ONLY: veg_parameter_type
    USE casaparm

USE cable_pack_mod, ONLY: cable_pack_rr, pack_landpts2mp
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft 
INTEGER, INTENT(IN) :: row_length, rows
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft )  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts, nsurft)
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)          ! cell index of land_pt
REAL,    INTENT(IN) :: latitude( row_length,rows )
REAL,    INTENT(IN) :: longitude( row_length,rows )

    TYPE (veg_parameter_type)           :: veg
    TYPE (soil_parameter_type)          :: soil
    TYPE (casa_flux), INTENT(INOUT)     :: casaflux
    TYPE (casa_met) , INTENT(INOUT)     :: casamet
    REAL, DIMENSION(land_pts)       :: soil_order
    REAL, DIMENSION(land_pts)       :: nidep,nifix,pwea,pdust
! local variables
    INTEGER :: k,p,i
    !REAL,DIMENSION(mp)             :: annNdep,annNfix,annPwea,annPdust
    REAL(r_2),DIMENSION(mp)             :: annNdep,annNfix,annPwea,annPdust
    !REAL(r_2),DIMENSION(:),ALLOCATABLE :: annNdep,annNfix,annPwea,annPdust
    REAL(r_2)                           :: annNfert,annPfert
    !REAL(r_2),ALLOCATABLE               :: annNfert,annPfert
    LOGICAL                             :: skip =.TRUE.
    INTEGER           :: sorder( land_pts )

! initialise 
    sorder(:)  = 0
    annNdep(:) = 0.0; annNfix(:) = 0.0 
    annPwea(:) = 0.0; annPdust(:)= 0.0
    annPfert = 0.7/365.0
    annNfert = 4.3/365.0

CALL cable_pack_rr( casamet%lat, latitude, mp, l_tile_pts, row_length,    &
                          rows, nsurft, land_pts, land_index, surft_pts,        &
                          surft_index )

CALL cable_pack_rr( casamet%lon, longitude, mp, l_tile_pts, row_length,    &
                          rows, nsurft, land_pts, land_index, surft_pts,        &
                          surft_index )

! Lest Nov2011 - not correct, but areacell not needed/used for UM ?
! Lest Feb2014 - use UM's A_BOXAREAS
    casamet%areacell = PACK( tile_frac, l_tile_pts )

sorder = INT(soil_order)

CALL um2cable_ilp( mp, land_pts, nsurft, surft_pts, surft_index,               &
                   l_tile_pts, sorder, sorder(1:10), casamet%isorder,          &
                   soil%isoilm, skip )

CALL pack_landpts2mp( nsurft, land_pts, mp, surft_pts, surft_index,            &
                      L_tile_pts, Nifix, annNfix )
  
CALL pack_landpts2mp( nsurft, land_pts, mp, surft_pts, surft_index,            &
                      L_tile_pts, Pwea, annPwea )
  
CALL pack_landpts2mp( nsurft, land_pts, mp, surft_pts, surft_index,            &
                      L_tile_pts, Pdust, annPdust )

    !casaflux%Nmindep = annNdep/365.0      ! gN/m2/day
    casaflux%Nminfix = annNfix/365.0      ! gN/m2/day
    casaflux%Pdep    = annPdust/365.0     ! gP/m2/day
    casaflux%Pwea    = annPwea/365.0      ! gP/m2/day

    do p = 1,mp
       if (veg%iveg(p)==cropland .or. veg%iveg(p)==croplnd2) then
       ! P fertilizer =13 Mt P globally in 1994
         casaflux%Pdep(p) = casaflux%Pdep(p)+annPfert
       ! N fertilizer =86 Mt N globally in 1994
         casaflux%Nminfix(p) = casaflux%Nminfix(p)+annNfert
       endif
    enddo

END SUBROUTINE casa_readpoint_pk




!========================================================================
!========================================================================
!========================================================================

SUBROUTINE casa_init_pk( mp, land_pts, nsurft, row_length, rows, l_tile_pts,   &
                   tile_frac, casabiome, casaflux, casamet, casapool, casabal, &
                   veg, canopy, phen, cpool_tile, npool_tile, ppool_tile,      &
                   wood_hvest_c, wood_hvest_n, wood_hvest_p, wresp_c, wresp_n, &
                   wresp_p, thinning, GLAI, PHENPHASE, PREV_YR_SFRAC, idoy )
                   
                       
!  initialize some values in phenology parameters and leaf growth phase

USE cable_def_types_mod,    ONLY: canopy_type 
USE cable_def_types_mod,    ONLY: veg_parameter_type

USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_biome_type_mod,   ONLY: casa_biome   => casa_biome_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE phenology_type_mod,    ONLY: phenology_type
USE progs_cnp_vars_mod,    ONLY: nCpool_casa, nNpool_casa, nPPool_casa

  USE casadimension        ! icycle,mplant,mlitter,msoil
  USE cable_common_module, ONLY : ktau_gl, l_luc
  USE casaparm             !, ONLY : initcasa

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft 
INTEGER, INTENT(IN) :: row_length, rows
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft )  ! TRUE if active tile
REAL,    INTENT(IN) :: tile_frac(land_pts, nsurft)

TYPE (casa_biome),         INTENT(INOUT) :: casabiome
TYPE (casa_flux),          INTENT(INOUT) :: casaflux
TYPE (casa_met),           INTENT(INOUT) :: casamet
TYPE (casa_pool),          INTENT(INOUT) :: casapool
TYPE (casa_balance),       INTENT(INOUT) :: casabal
TYPE (canopy_type),        INTENT(INOUT) :: canopy
TYPE (veg_parameter_type), INTENT(INOUT) :: veg
TYPE (phenology_type),     INTENT(INOUT) :: phen

REAL, INTENT(INOUT) :: Cpool_tile( land_pts, nsurft, nCpool_casa )
REAL, INTENT(INOUT) :: Npool_tile( land_pts, nsurft, nNpool_casa )
REAL, INTENT(INOUT) :: Ppool_tile( land_pts, nsurft, nPPool_casa )
  REAL   , INTENT(INOUT) :: glai          ( land_pts, nsurft ) 
  REAL   , INTENT(INOUT) :: phenphase     ( land_pts, nsurft ) 
  REAL   , INTENT(INOUT) :: prev_yr_sfrac ( land_pts, nsurft ) 
  REAL   , INTENT(INOUT) :: wood_hvest_C  ( land_pts, nsurft, 3 )
  REAL   , INTENT(INOUT) :: wood_hvest_N  ( land_pts, nsurft, 3 )
  REAL   , INTENT(INOUT) :: wood_hvest_P  ( land_pts, nsurft, 3 )
  REAL   , INTENT(INOUT) :: wresp_C       ( land_pts, nsurft, 3 )
  REAL   , INTENT(INOUT) :: wresp_N       ( land_pts, nsurft, 3 )
  REAL   , INTENT(INOUT) :: wresp_P       ( land_pts, nsurft, 3 )
  REAL   , INTENT(INOUT) :: thinning      ( land_pts, nsurft )
  INTEGER :: idoy,mtau

  casamet%tairk(:)   = 0.0
  casamet%tsoil(:,:) = 0.0
  casamet%moist(:,:) = 0.0
  casaflux%cgpp(:)   = 0.0
  casaflux%cnpp(:)   = canopy%fnpp*86400.     !0.0 TZ initilize npp from prognostic
  casaflux%Crsoil(:) = canopy%frs*86400.     !0.0
  casaflux%crgplant(:)   = canopy%frp*86400. !0.0 TZ this might not be corrext crgplant ne frp 
  casaflux%crmplant(:,:) = 0.0
  casaflux%clabloss(:)   = 0.0

  ! Lest 19/2/14 - will work for coupled and amip
  ! mtau is the step number of the day (1,2,..47,0)
  mtau = mod(ktau_gl,int(24.*3600./ktau_gl))

IF (initcasa==1) THEN
  
  IF (l_luc .and. idoy == 1 .and. mtau == 1) THEN

    CALL casa_reinit_pk( land_pts, nsurft, L_tile_pts, tile_frac, casabiome,   &
                         casamet, casapool, casabal, veg, phen, Cpool_tile,    &
                         Npool_tile, Ppool_tile, wood_hvest_C, wood_hvest_N,   &
                         wood_hvest_P, wresp_C, wresp_N, wresp_P, thinning,    &
                         gLAI, phenphase, prev_yr_sfrac )
                         
                         

  ELSE ! (l_luc .and. idoy == 1 .and. mtau == 1)

    CALL pack_cnppool( land_pts, nsurft, L_tile_pts, casamet, casapool, casabal, phen, &
                       Cpool_tile, Npool_tile, Ppool_tile, gLAI, phenphase )
                          
  ENDIF ! (l_luc .and. idoy == 1 .and. mtau == 1)

ENDIF ! initcasa

! reset labile C pool,comment out by Q.Zhang 10/09/2011
  casapool%cplant     = MAX(0.0,casapool%cplant)
  casapool%clitter    = MAX(0.0,casapool%clitter)
  casapool%csoil      = MAX(0.0,casapool%csoil)
  casabal%cplantlast  = casapool%cplant
  casabal%clitterlast = casapool%clitter
  casabal%csoillast   = casapool%csoil
  casabal%clabilelast = casapool%clabile
  casabal%sumcbal     = 0.0
       casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
       casabal%FCnppyear=0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0

  IF (icycle==1) THEN
    casapool%nplant(:,:) = casapool%cplant(:,:) * casapool%rationcplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF

  IF (icycle >1) THEN
    casapool%nplant      = MAX(1.e-6,casapool%nplant)
    casapool%nlitter     = MAX(1.e-6,casapool%nlitter)
    casapool%nsoil       = MAX(1.e-6,casapool%nsoil)
    casapool%nsoilmin    = MAX(1.e-6,casapool%nsoilmin)
    casabal%nplantlast   = casapool%nplant
    casabal%nlitterlast  = casapool%nlitter
    casabal%nsoillast    = casapool%nsoil
    casabal%nsoilminlast = casapool%nsoilmin
    casabal%sumnbal      = 0.0
       casabal%FNdepyear=0.0;casabal%FNfixyear=0.0;casabal%FNsnetyear=0.0
       casabal%FNupyear=0.0;casabal%FNleachyear=0.0;casabal%FNlossyear=0.0
  ENDIF

  IF (icycle >2) THEN
    casapool%pplant       = MAX(1.0e-7,casapool%pplant)
    casapool%plitter      = MAX(1.0e-7,casapool%plitter)
    casapool%psoil        = MAX(1.0e-7,casapool%psoil)
    casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0, changed according to  YP
    casapool%psoilsorb    = MAX(1.0e-7,casapool%psoilsorb) ! was 10.0, -
    casapool%psoilocc     = MAX(1.0e-7,casapool%psoilocc)  ! was 50.0, -
    casabal%pplantlast    = casapool%pplant
    casabal%plitterlast   = casapool%plitter
    casabal%psoillast     = casapool%psoil
    casabal%psoillablast  = casapool%psoillab
    casabal%psoilsorblast = casapool%psoilsorb
    casabal%psoilocclast  = casapool%psoilocc
    casabal%sumpbal       = 0.0
       casabal%FPweayear=0.0;casabal%FPdustyear=0.0; casabal%FPsnetyear=0.0
       casabal%FPupyear=0.0;casabal%FPleachyear=0.0;casabal%FPlossyear=0.0
  ENDIF

END SUBROUTINE casa_init_pk

SUBROUTINE casa_reinit_pk( land_pts, nsurft, L_tile_pts, tile_frac, casabiome, &
                           casamet, casapool, casabal, veg, phen, Cpool_tile,  &
                           Npool_tile, Ppool_tile, woodhvest_C, woodhvest_N,   &
                           woodhvest_p, wresp_C, wresp_N, wresp_P, thinning,   &
                           gLAI, phenphase, prev_yr_sfrac )

   USE cable_def_types_mod ! combines def_dimensions (mp,r_2) and define_types (mland)  
   USE casadimension
   USE casaparm

USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_biome_type_mod,   ONLY: casa_biome   => casa_biome_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE phenology_type_mod,    ONLY: phenology_type
USE progs_cnp_vars_mod,    ONLY: nCpool_casa, nNpool_casa, nPPool_casa
USE cable_common_module, ONLY : ktau_gl, l_thinforest
   
IMPLICIT NONE

INTEGER,                  INTENT(IN) :: land_pts         
INTEGER,                  INTENT(IN) :: nsurft 
LOGICAL,                  INTENT(IN) :: L_tile_pts(land_pts, nsurft )  
REAL,                     INTENT(IN) :: tile_frac(land_pts, nsurft)
TYPE(casa_biome),         INTENT(INOUT) :: casabiome
TYPE(casa_met),           INTENT(INOUT) :: casamet
TYPE(casa_pool),          INTENT(INOUT) :: casapool
TYPE(casa_balance),       INTENT(INOUT) :: casabal    
TYPE(veg_parameter_type), INTENT(INOUT) :: veg
TYPE(phenology_type),     INTENT(INOUT) :: phen

REAL,  INTENT(INOUT) :: Cpool_tile    ( land_pts, nsurft, nCpool_casa )
REAL,  INTENT(INOUT) :: Npool_tile    ( land_pts, nsurft, nNpool_casa )
REAL,  INTENT(INOUT) :: Ppool_tile    ( land_pts, nsurft, nPPool_casa )
REAL,  INTENT(INOUT) :: glai          ( land_pts, nsurft)
REAL,  INTENT(INOUT) :: phenphase     ( land_pts, nsurft)
REAL,  INTENT(INOUT) :: prev_yr_sfrac ( land_pts, nsurft)

! local variables
REAL(r_2) :: clabile_x  ( land_pts, nsurft)
REAL(r_2) :: cplant_x   ( land_pts, nsurft, mplant )
REAL(r_2) :: clitter_x  ( land_pts, nsurft, mlitter )
REAL(r_2) :: csoil_x    ( land_pts, nsurft, msoil )
REAL(r_2) :: nplant_x   ( land_pts, nsurft, mplant )
REAL(r_2) :: nlitter_x  ( land_pts, nsurft, mlitter )
REAL(r_2) :: nsoil_x    ( land_pts, nsurft, msoil )
REAL(r_2) :: nsoilmin_x ( land_pts, nsurft)
REAL(r_2) :: pplant_x   ( land_pts, nsurft, mplant )
REAL(r_2) :: plitter_x  ( land_pts, nsurft, mlitter )
REAL(r_2) :: psoil_x    ( land_pts, nsurft, msoil )
REAL(r_2) :: psoillab_x ( land_pts, nsurft)
REAL(r_2) :: psoilsorb_x( land_pts, nsurft)
REAL(r_2) :: psoilocc_x ( land_pts, nsurft)

REAL(r_2) :: clabile_y  ( land_pts, nsurft )
REAL(r_2) :: cplant_y   ( land_pts, nsurft, mplant )
REAL(r_2) :: clitter_y  ( land_pts, nsurft, mlitter )
REAL(r_2) :: csoil_y    ( land_pts, nsurft, msoil )
REAL(r_2) :: nplant_y   ( land_pts, nsurft, mplant )
REAL(r_2) :: nlitter_y  ( land_pts, nsurft, mlitter )
REAL(r_2) :: nsoil_y    ( land_pts, nsurft, msoil )
REAL(r_2) :: nsoilmin_y ( land_pts, nsurft )
REAL(r_2) :: pplant_y   ( land_pts, nsurft, mplant )
REAL(r_2) :: plitter_y  ( land_pts, nsurft, mlitter )
REAL(r_2) :: psoil_y    ( land_pts, nsurft, msoil )
REAL(r_2) :: psoillab_y ( land_pts, nsurft )
REAL(r_2) :: psoilsorb_y( land_pts, nsurft )
REAL(r_2) :: psoilocc_y ( land_pts, nsurft )

REAL    :: frac_x  ( land_pts, nsurft )
LOGICAL :: ifpre_x ( land_pts, nsurft )
REAL    :: frac_y  ( land_pts, nsurft )
LOGICAL :: ifpre_y ( land_pts, nsurft )

! To be recorded wood log variables.
! wood_flux
REAL(r_2) :: logc ( land_pts, nsurft )
REAL(r_2) :: logn ( land_pts, nsurft )
REAL(r_2) :: logp ( land_pts, nsurft )


REAL(r_2) :: woodhvest_c ( land_pts, nsurft, 3 )
REAL(r_2) :: woodhvest_n ( land_pts, nsurft, 3 )
REAL(r_2) :: woodhvest_p ( land_pts, nsurft, 3 )

REAL(r_2) :: wresp_c ( land_pts, nsurft, 3 )
REAL(r_2) :: wresp_n ( land_pts, nsurft, 3 )
REAL(r_2) :: wresp_p ( land_pts, nsurft, 3 )

REAL(r_2) :: thinning    ( land_pts, nsurft )

REAL, PARAMETER:: POOL_FRAC(3) =(/0.33, 0.33, 0.34/)
REAL, PARAMETER:: POOL_TIME(3) =(/1.00, 0.10, 0.01/)
REAL(r_2) :: cplant_z    ( land_pts, nsurft, mplant )
REAL(r_2) :: nplant_z    ( land_pts, nsurft, mplant )
REAL(r_2) :: pplant_z    ( land_pts, nsurft, mplant )

! check if all of them are required
INTEGER   :: g, p, k, y ! np
REAL(r_2) :: cbal, nbal, pbal

! Initialize temporary variables

Cplant_x = 0.0
Nplant_x = 0.0
Pplant_x = 0.0

Clitter_x = 0.0
Nlitter_x = 0.0
Plitter_x = 0.0

Csoil_x = 0.0
Nsoil_x = 0.0
Psoil_x = 0.0

Clabile_x   = 0.0
Nsoilmin_x  = 0.0
Psoillab_x  = 0.0
Psoilsorb_x = 0.0
Psoilocc_x  = 0.0

Cplant_y = 0.0
Nplant_y = 0.0
Pplant_y = 0.0

Clitter_y = 0.0
Nlitter_y = 0.0
Plitter_y = 0.0

Csoil_y = 0.0
Nsoil_y = 0.0
Psoil_y = 0.0

Clabile_y   = 0.0
Nsoilmin_y  = 0.0
Psoillab_y  = 0.0
Psoilsorb_y = 0.0
Psoilocc_y  = 0.0

frac_x  = 0.0
frac_y  = 0.0
ifpre_x = .FALSE.
ifpre_y = .FALSE.

logC = 0.0
logN = 0.0
logP = 0.0

Cplant_z = 0.0
Nplant_z = 0.0
Pplant_z = 0.0

wresp_C = 0.0
wresp_N = 0.0
wresp_P = 0.0

thinning = 0.0

! assign "old" cnp pool values (from last dump file, initilization)
Clabile_x(:,:)   = Cpool_tile(:,:,1)
Cplant_x(:,:,1)  = Cpool_tile(:,:,2)
Cplant_x(:,:,2)  = Cpool_tile(:,:,3)
Cplant_x(:,:,3)  = Cpool_tile(:,:,4)
Clitter_x(:,:,1) = Cpool_tile(:,:,5)
Clitter_x(:,:,2) = Cpool_tile(:,:,6)
Clitter_x(:,:,3) = Cpool_tile(:,:,7)
Csoil_x(:,:,1)   = Cpool_tile(:,:,8)
Csoil_x(:,:,2)   = Cpool_tile(:,:,9)
Csoil_x(:,:,3)   = Cpool_tile(:,:,10)

IF (icycle>1) THEN
   Nplant_x(:,:,1)  = Npool_tile(:,:,1)
   Nplant_x(:,:,2)  = Npool_tile(:,:,2)
   Nplant_x(:,:,3)  = Npool_tile(:,:,3)
   Nlitter_x(:,:,1) = Npool_tile(:,:,4)
   Nlitter_x(:,:,2) = Npool_tile(:,:,5)
   Nlitter_x(:,:,3) = Npool_tile(:,:,6)
   Nsoil_x(:,:,1)   = Npool_tile(:,:,7)
   Nsoil_x(:,:,2)   = Npool_tile(:,:,8)
   Nsoil_x(:,:,3)   = Npool_tile(:,:,9)
   Nsoilmin_x(:,:)  = Npool_tile(:,:,10)
END IF

IF (icycle>2) THEN
   Pplant_x(:,:,1)  = Ppool_tile(:,:,1)
   Pplant_x(:,:,2)  = Ppool_tile(:,:,2)
   Pplant_x(:,:,3)  = Ppool_tile(:,:,3)
   Plitter_x(:,:,1) = Ppool_tile(:,:,4)
   Plitter_x(:,:,2) = Ppool_tile(:,:,5)
   Plitter_x(:,:,3) = Ppool_tile(:,:,6)
   Psoil_x(:,:,1)   = Ppool_tile(:,:,7)
   Psoil_x(:,:,2)   = Ppool_tile(:,:,8)
   Psoil_x(:,:,3)   = Ppool_tile(:,:,9)
   Psoillab_x(:,:)  = Ppool_tile(:,:,10)
   Psoilsorb_x(:,:) = Ppool_tile(:,:,11)
   Psoilocc_x(:,:)  = Ppool_tile(:,:,12)
END IF

! assign fractions (previous) (need to get fractions from previous year)
frac_x(:,:) = prev_yr_sfrac(:,:)

! assign fractions (current)
frac_y(:,:) = tile_frac(:,:)

! set the ifpre_x and ifpre_y values to true where fractions are .ne. 0
WHERE(frac_x > 0.) ifpre_x = .TRUE.
WHERE(frac_y > 0.) ifpre_y = .TRUE.

  ! start main loop
  DO g = 1, land_pts

     ! Check all glacier tiles
     IF (ifpre_x(g,iceland).and.ifpre_y(g,iceland)) THEN
         IF (abs(frac_x(g,iceland)-1.0)>0.01 .or. abs(frac_y(g,iceland)-1.0)>0.01) THEN
         print *,'Lest Glacier',g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
            STOP "Glacier fraction .ne. 1"
         ELSE
            print *, 'Lest cycle', g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
            cycle
         END IF
     ELSEIF (.not.ifpre_x(g,iceland) .and. .not.ifpre_y(g,iceland)) THEN
         IF (abs(frac_x(g,iceland)-0.0)>0.01 .or. abs(frac_y(g,iceland)-0.0)>0.01) THEN
         print *, 'Lest landpt fracs x', frac_x(g,:)
         print *, 'Lest landpt fracs y', frac_y(g,:)
         print *,'Lest Glacier',g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
            STOP "Glacier fraction .ne. 0"
         END IF
     ELSE
         print *,'Lest Glacier',g,ifpre_x(g,iceland),ifpre_y(g,iceland),&
                                 frac_x(g,iceland), frac_y(g,iceland)
         STOP "Glacier tiles not consistent"
     END IF

     ! For none glacier tiles
     ! Re-calculate plant C, N, P pools
     CALL newplant(cplant_x(g,:,:),frac_x(g,:),ifpre_x(g,:), &
                   cplant_y(g,:,:),frac_y(g,:),ifpre_y(g,:),logc(g,:))
     IF (icycle > 1) CALL newplant(nplant_x(g,:,:),frac_x(g,:),ifpre_x(g,:), &
                                   nplant_y(g,:,:),frac_y(g,:),ifpre_y(g,:),logn(g,:))
     IF (icycle > 2) CALL newplant(pplant_x(g,:,:),frac_x(g,:),ifpre_x(g,:), &
                                   pplant_y(g,:,:),frac_y(g,:),ifpre_y(g,:),logp(g,:))

     DO y = 1,3 ! pools NOT leaf/wood/root
                      woodhvest_c(g,:,y) = woodhvest_c(g,:,y) + pool_frac(y)*logc(g,:)     !slogc
      IF (icycle > 1) woodhvest_n(g,:,y) = woodhvest_n(g,:,y) + pool_frac(y)*logn(g,:)     !slogn
      IF (icycle > 2) woodhvest_p(g,:,y) = woodhvest_p(g,:,y) + pool_frac(y)*logp(g,:)     !slogp
     END DO

     ! Re-calculate litter C, N, P pools
     CALL newlitter(casabiome,frac_x(g,:),ifpre_x(g,:),frac_y(g,:),ifpre_y(g,:), &
                    cplant_x(g,:,:),nplant_x(g,:,:),pplant_x(g,:,:), &
                    cplant_y(g,:,:),nplant_y(g,:,:),pplant_y(g,:,:), &
                    clitter_x(g,:,:),nlitter_x(g,:,:),plitter_x(g,:,:), &
                    clitter_y(g,:,:),nlitter_y(g,:,:),plitter_y(g,:,:))

     ! Re-calculate soil C, N, P pools
     CALL newsoil(msoil,csoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
                  csoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
     CALL newsoil(1,clabile_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                  clabile_y(g,:),frac_y(g,:),ifpre_y(g,:))
     IF (icycle > 1) THEN
        CALL newsoil(msoil,nsoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
                     nsoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,nsoilmin_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     nsoilmin_y(g,:),frac_y(g,:),ifpre_y(g,:))
     ENDIF

     IF (icycle > 2) THEN
        CALL newsoil(msoil,psoil_x(g,:,:),frac_x(g,:),ifpre_x(g,:),&
                     psoil_y(g,:,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,psoillab_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     psoillab_y(g,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,psoilsorb_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     psoilsorb_y(g,:),frac_y(g,:),ifpre_y(g,:))
        CALL newsoil(1,psoilocc_x(g,:),frac_x(g,:),ifpre_x(g,:),&
                     psoilocc_y(g,:),frac_y(g,:),ifpre_y(g,:))
     ENDIF


     ! TEST Lestevens 6june18 - thinning forests after luc ----
     IF (l_thinforest) THEN
                      cplant_z(g,:,:) = cplant_y(g,:,:)
      if (icycle > 1) nplant_z(g,:,:) = nplant_y(g,:,:)
      if (icycle > 2) pplant_z(g,:,:) = pplant_y(g,:,:)
      DO y=1,3 ! pools for whvest
                        woodhvest_c(g,:,y) = woodhvest_c(g,:,y) + &
                                   (1-thinning(g,:)) * pool_frac(y) * cplant_y(g,:,wood)
        if (icycle > 1) woodhvest_n(g,:,y) = woodhvest_n(g,:,y) + &
                                   (1-thinning(g,:)) * pool_frac(y) * nplant_y(g,:,wood)
        if (icycle > 2) woodhvest_p(g,:,y) = woodhvest_p(g,:,y) + &
                                   (1-thinning(g,:)) * pool_frac(y) * pplant_y(g,:,wood)
      END DO
      DO y=1,mplant
                        cplant_z(g,:,y) = thinning(g,:) * cplant_y(g,:,y)
        if (icycle > 1) nplant_z(g,:,y) = thinning(g,:) * nplant_y(g,:,y)
        if (icycle > 2) pplant_z(g,:,y) = thinning(g,:) * pplant_y(g,:,y)
      END DO
      CALL newlitter_thin(casabiome,frac_y(g,:),ifpre_y(g,:),frac_y(g,:),ifpre_y(g,:), &
                     cplant_y(g,:,:),nplant_y(g,:,:),pplant_y(g,:,:), &
                     cplant_z(g,:,:),nplant_z(g,:,:),pplant_z(g,:,:), &
                     clitter_y(g,:,:),nlitter_y(g,:,:),plitter_y(g,:,:), &
                     clitter_y(g,:,:),nlitter_y(g,:,:),plitter_y(g,:,:),thinning(g,:))
                      cplant_y(g,:,:) = cplant_z(g,:,:)
      if (icycle > 1) nplant_y(g,:,:) = nplant_z(g,:,:)
      if (icycle > 2) pplant_y(g,:,:) = pplant_z(g,:,:)
     ENDIF

! Lestevens 24 Nov 2017 - Implement Yingping's flux to state pools
! Lestevens 7 June 2018 - Moved to after thinning
     !DATA pool_frac/0.33,0.33,0.34/
     !DATA pool_time/1.0 ,0.1 ,0.01/
     DO y = 1,3 ! pools NOT leaf/wood/root
                      wresp_c(g,:,y)     = woodhvest_c(g,:,y) * (1.-exp(-1.*pool_time(y))) !
      IF (icycle > 1) wresp_n(g,:,y)     = woodhvest_n(g,:,y) * (1.-exp(-1.*pool_time(y))) !
      IF (icycle > 2) wresp_p(g,:,y)     = woodhvest_p(g,:,y) * (1.-exp(-1.*pool_time(y))) !
                      woodhvest_c(g,:,y) = woodhvest_c(g,:,y) - wresp_c(g,:,y)
      IF (icycle > 1) woodhvest_n(g,:,y) = woodhvest_n(g,:,y) - wresp_n(g,:,y)
      IF (icycle > 2) woodhvest_p(g,:,y) = woodhvest_p(g,:,y) - wresp_p(g,:,y)
     END DO
! Lestevens 24 Nov 2017 - Implement Yingping's flux to state pools


     ! Balance check
     cbal = sum((sum(cplant_x(g,:,:),2) + sum(clitter_x(g,:,:),2)  &
          + sum(csoil_x(g,:,:),2) + clabile_x(g,:)) * frac_x(g,:))  &
          - (sum((sum(cplant_y(g,:,:),2) + sum(clitter_y(g,:,:),2)  &
          + sum(csoil_y(g,:,:),2) + clabile_y(g,:)) * frac_y(g,:)) + sum(logc(g,:)))

     IF (icycle > 1) nbal = sum((sum(nplant_x(g,:,:),2) + sum(nlitter_x(g,:,:),2)  &
                          + sum(nsoil_x(g,:,:),2) + nsoilmin_x(g,:)) * frac_x(g,:))  &
                          - (sum((sum(nplant_y(g,:,:),2) + sum(nlitter_y(g,:,:),2)  &
                          + sum(nsoil_y(g,:,:),2) + nsoilmin_y(g,:)) * frac_y(g,:)) + sum(logn(g,:)))

     IF (icycle > 2) pbal = sum((sum(pplant_x(g,:,:),2) + sum(plitter_x(g,:,:),2)  &
                          + sum(psoil_x(g,:,:),2) + psoillab_x(g,:)  &
                          + psoilsorb_x(g,:) + psoilocc_x(g,:)) * frac_x(g,:))  &
                          - (sum((sum(pplant_y(g,:,:),2) + sum(plitter_y(g,:,:),2)  &
                          + sum(psoil_y(g,:,:),2) + psoillab_y(g,:) + psoilsorb_y(g,:)  &
                          + psoilocc_y(g,:)) * frac_y(g,:)) + sum(logp(g,:)))


     IF(abs(cbal)>1.e-3 .or.abs(nbal)>1.e-4 .or.abs(pbal)>1.e-4) THEN
        print*, 'imbalance on grid:',g,cbal,nbal,pbal
     END IF

  END DO ! end main loop
  !END DO ! end main loop

! write values back into cnp pool files
Cpool_tile(:,:,1)  = Clabile_y(:,:)
Cpool_tile(:,:,2)  = Cplant_y(:,:,1)
Cpool_tile(:,:,3)  = Cplant_y(:,:,2)
Cpool_tile(:,:,4)  = Cplant_y(:,:,3)
Cpool_tile(:,:,5)  = Clitter_y(:,:,1)
Cpool_tile(:,:,6)  = Clitter_y(:,:,2)
Cpool_tile(:,:,7)  = Clitter_y(:,:,3)
Cpool_tile(:,:,8)  = Csoil_y(:,:,1)
Cpool_tile(:,:,9)  = Csoil_y(:,:,2)
Cpool_tile(:,:,10) = Csoil_y(:,:,3)

! if n is switched on
IF (icycle > 1) THEN
   Npool_tile(:,:,1)  = Nplant_y(:,:,1)
   Npool_tile(:,:,2)  = Nplant_y(:,:,2)
   Npool_tile(:,:,3)  = Nplant_y(:,:,3)
   Npool_tile(:,:,4)  = Nlitter_y(:,:,1)
   Npool_tile(:,:,5)  = Nlitter_y(:,:,2)
   Npool_tile(:,:,6)  = Nlitter_y(:,:,3)
   Npool_tile(:,:,7)  = Nsoil_y(:,:,1)
   Npool_tile(:,:,8)  = Nsoil_y(:,:,2)
   Npool_tile(:,:,9)  = Nsoil_y(:,:,3)
   Npool_tile(:,:,10) = Nsoilmin_y(:,:)     
END IF

! if p is switched on
IF (icycle > 2) THEN
   Ppool_tile(:,:,1)  = Pplant_y(:,:,1)
   Ppool_tile(:,:,2)  = Pplant_y(:,:,2)
   Ppool_tile(:,:,3)  = Pplant_y(:,:,3)
   Ppool_tile(:,:,4)  = Plitter_y(:,:,1)
   Ppool_tile(:,:,5)  = Plitter_y(:,:,2)
   Ppool_tile(:,:,6)  = Plitter_y(:,:,3)
   Ppool_tile(:,:,7)  = Psoil_y(:,:,1)
   Ppool_tile(:,:,8)  = Psoil_y(:,:,2)
   Ppool_tile(:,:,9)  = Psoil_y(:,:,3)
   Ppool_tile(:,:,10) = Psoillab_y(:,:)
   Ppool_tile(:,:,11) = Psoilsorb_y(:,:)
   Ppool_tile(:,:,12) = Psoilocc_y(:,:)
END IF

! pack everything
CALL pack_cnppool( land_pts, nsurft, L_tile_pts, casamet, casapool, casabal, phen, &
                   Cpool_tile, Npool_tile, Ppool_tile, gLAI, phenphase )


  ! update LAI, Vcmax
  IF (icycle > 1) call casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)

  DO p=1,mp
     IF (casamet%iveg2(p) == icewater) THEN
        casamet%glai(p)   = 0.0
     ELSE
        casamet%glai(p)   = MIN(casabiome%glaimax(veg%iveg(p)), MAX(casabiome%glaimin(veg%iveg(p)), &
                               casabiome%sla(veg%iveg(p)) * casapool%cplant(p,leaf)))
     ENDIF

  END DO

  

END SUBROUTINE casa_reinit_pk

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE pack_cnppool( land_pts, nsurft, L_tile_pts, casamet, casapool, casabal, phen, Cpool_tile,         &
                         Npool_tile, Ppool_tile,  gLAI, phenphase )

    USE casadimension        ! icycle,mplant,mlitter,msoil

USE phenology_type_mod,    ONLY: phenology_type
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE progs_cnp_vars_mod,    ONLY: nCpool_casa, nNpool_casa, nPPool_casa

IMPLICIT NONE

! passed in
INTEGER,              INTENT(IN) :: land_pts  
INTEGER,              INTENT(IN) :: nsurft 
LOGICAL,              INTENT(IN) :: L_tile_pts(land_pts, nsurft )
TYPE (casa_met),      INTENT(INOUT) :: casamet
TYPE (casa_pool),     INTENT(INOUT) :: casapool
TYPE (casa_balance),  INTENT(INOUT) :: casabal
TYPE (phenology_type), INTENT(INOUT) :: phen

REAL, INTENT(INOUT) :: Cpool_tile( land_pts, nsurft, nCpool_casa )
REAL, INTENT(INOUT) :: Npool_tile( land_pts, nsurft, nNpool_casa )
REAL, INTENT(INOUT) :: Ppool_tile( land_pts, nsurft, nPPool_casa )
REAL, INTENT(INOUT) :: gLAI      ( land_pts, nsurft)
REAL, INTENT(INOUT) :: phenphase ( land_pts, nsurft)

! local vars
INTEGER k
REAL :: Clabile   ( land_pts, nsurft ) 
REAL :: Nsoilmin  ( land_pts, nsurft )
REAL :: Psoillab  ( land_pts, nsurft )
REAL :: Psoilsorb ( land_pts, nsurft ) 
REAL :: Psoilocc  ( land_pts, nsurft ) 

REAL :: Cplant ( land_pts, nsurft, mplant )
REAL :: Nplant ( land_pts, nsurft, mplant )
REAL :: Pplant ( land_pts, nsurft, mplant )

REAL :: Clitter ( land_pts, nsurft, mlitter )
REAL :: Nlitter ( land_pts, nsurft, mlitter )
REAL :: Plitter ( land_pts, nsurft, mlitter )

REAL :: Csoil ( land_pts, nsurft, msoil )
REAL :: Nsoil ( land_pts, nsurft, msoil )
REAL :: Psoil ( land_pts, nsurft, msoil )

INTEGER  :: phenph ( land_pts, nsurft )

! initialise variables
Clabile(:,:)  = 0.0 
Cplant(:,:,:) = 0.0; Clitter(:,:,:) = 0.0; Csoil(:,:,:) = 0.0
Nplant(:,:,:) = 0.0; Nlitter(:,:,:) = 0.0; Nsoil(:,:,:) = 0.0 
Nsoilmin(:,:) = 0.0
Pplant(:,:,:) = 0.0; Plitter(:,:,:) = 0.0; Psoil(:,:,:) = 0.0 
Psoillab(:,:) = 0.0; Psoilsorb(:,:) = 0.0; Psoilocc(:,:)= 0.0
Phenph(:,:)   = 0

! set to appropriate pools
Clabile(:,:)   = Cpool_tile(:,:,1)
Cplant(:,:,1)  = Cpool_tile(:,:,2)
Cplant(:,:,2)  = Cpool_tile(:,:,3)
Cplant(:,:,3)  = Cpool_tile(:,:,4)
Clitter(:,:,1) = Cpool_tile(:,:,5)
Clitter(:,:,2) = Cpool_tile(:,:,6)
Clitter(:,:,3) = Cpool_tile(:,:,7)
Csoil(:,:,1)   = Cpool_tile(:,:,8)
Csoil(:,:,2)   = Cpool_tile(:,:,9)
Csoil(:,:,3)   = Cpool_tile(:,:,10)

IF (icycle >1) THEN
  Nplant(:,:,1)  = Npool_tile(:,:,1)
  Nplant(:,:,2)  = Npool_tile(:,:,2)
  Nplant(:,:,3)  = Npool_tile(:,:,3)
  Nlitter(:,:,1) = Npool_tile(:,:,4)
  Nlitter(:,:,2) = Npool_tile(:,:,5)
  Nlitter(:,:,3) = Npool_tile(:,:,6)
  Nsoil(:,:,1)   = Npool_tile(:,:,7)
  Nsoil(:,:,2)   = Npool_tile(:,:,8)
  Nsoil(:,:,3)   = Npool_tile(:,:,9)
  Nsoilmin(:,:)  = Npool_tile(:,:,10)
ENDIF

IF (icycle >2) THEN
  Pplant(:,:,1)  = Ppool_tile(:,:,1)
  Pplant(:,:,2)  = Ppool_tile(:,:,2)
  Pplant(:,:,3)  = Ppool_tile(:,:,3)
  Plitter(:,:,1) = Ppool_tile(:,:,4)
  Plitter(:,:,2) = Ppool_tile(:,:,5)
  Plitter(:,:,3) = Ppool_tile(:,:,6)
  Psoil(:,:,1)   = Ppool_tile(:,:,7)
  Psoil(:,:,2)   = Ppool_tile(:,:,8)
  Psoil(:,:,3)   = Ppool_tile(:,:,9)
  Psoillab(:,:)  = Ppool_tile(:,:,10)
  Psoilsorb(:,:) = Ppool_tile(:,:,11)
  Psoilocc(:,:)  = Ppool_tile(:,:,12)
ENDIF

! pack variables
casamet%glai = PACK( glai(:,:)  , l_tile_pts )
phenph       = INT( phenphase )
phen%phase   = PACK( phenph(:,:), l_tile_pts )

casapool%Clabile       = pack(Clabile(:,:), l_tile_pts)
DO k=1,3
  casapool%cplant(:,k)  = pack(cplant(:,:,k) ,l_tile_pts)
  casapool%clitter(:,k) = pack(clitter(:,:,k),l_tile_pts)
  casapool%csoil(:,k)   = pack(csoil(:,:,k)  ,l_tile_pts)
ENDDO

IF (icycle>1) THEN
  DO k=1,3
    casapool%nplant(:,k)  = PACK( Nplant(:,:,k),  l_tile_pts )
    casapool%nlitter(:,k) = PACK( Nlitter(:,:,k), l_tile_pts )
    casapool%nsoil(:,k)   = PACK( Nsoil(:,:,k),   l_tile_pts )
  ENDDO
  casapool%nsoilmin     = PACK( Nsoilmin(:,:), l_tile_pts)
ENDIF

IF (icycle>2) THEN
  DO k=1,3
    casapool%pplant(:,k)  = PACK( Pplant(:,:,k) , l_tile_pts )
    casapool%plitter(:,k) = PACK( Plitter(:,:,k), l_tile_pts )
    casapool%psoil(:,k)   = PACK( Psoil(:,:,k)  , l_tile_pts )
  ENDDO
  casapool%Psoillab     = PACK( Psoillab(:,:),  l_tile_pts )
  casapool%Psoilsorb    = PACK( Psoilsorb(:,:), l_tile_pts )
  casapool%Psoilocc     = PACK( Psoilocc(:,:),  l_tile_pts )
ENDIF

END SUBROUTINE pack_cnppool

SUBROUTINE casa_poolout_unpk( land_pts, nsurft, L_tile_pts, casapool, casaflux,&
                              casamet, casabal, phen, Cpool_tile, Npool_tile,  &
                              Ppool_tile, gLAI, phenphase )

  USE cable_def_types_mod
  USE casadimension        ! icycle,mplant,mlitter,msoil
USE phenology_type_mod,    ONLY: phenology_type
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_biome_type_mod,   ONLY: casa_biome   => casa_biome_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE progs_cnp_vars_mod,    ONLY: nCpool_casa, nNpool_casa, nPPool_casa


IMPLICIT NONE

INTEGER,                  INTENT(IN)    :: land_pts       
INTEGER,                  INTENT(IN)    :: nsurft 
LOGICAL,                  INTENT(IN)    :: L_tile_pts(land_pts, nsurft )  
TYPE (casa_pool),         INTENT(INOUT) :: casapool
TYPE (casa_flux),         INTENT(INOUT) :: casaflux
TYPE (casa_met),          INTENT(INOUT) :: casamet
TYPE (casa_balance),      INTENT(INOUT) :: casabal
TYPE (phenology_type),    INTENT(INOUT) :: phen

REAL, INTENT(INOUT) :: Cpool_tile( land_pts, nsurft, nCpool_casa )
REAL, INTENT(INOUT) :: Npool_tile( land_pts, nsurft, nNpool_casa )
REAL, INTENT(INOUT) :: Ppool_tile( land_pts, nsurft, nPPool_casa )
REAL, INTENT(INOUT) :: gLAI      ( land_pts, nsurft)
REAL, INTENT(INOUT) :: phenphase ( land_pts, nsurft)

! local variables
    REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
    REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
    REAL(r_2), DIMENSION(mp)  :: totpsoil
    INTEGER  npt,nso

  ! Soiltype     soilnumber soil P(g P/m2)
  ! Alfisol   1       61.3
  ! Andisol   2       103.9
  ! Aridisol  3       92.8
  ! Entisol   4       136.9
  ! Gellisol  5       98.2
  ! Histosol  6       107.6
  ! Inceptisol        7      84.1
  ! Mollisol  8       110.1
  ! Oxisol    9       35.4
  ! Spodosol  10      41.0
  ! Ultisol   11      51.5
  ! Vertisol  12      190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

!  WRITE(*,91) nyear,cplantsum,clittersum,csoilsum
  casabal%sumcbal=MIN(9999.0,MAX(-9999.0,casabal%sumcbal))
  casabal%sumnbal=MIN(9999.0,MAX(-9999.0,casabal%sumnbal))
  casabal%sumpbal=MIN(9999.0,MAX(-9999.0,casabal%sumpbal))

  DO npt =1, mp
    nso = casamet%isorder(npt)
    totpsoil(npt) = psorder(nso) *xpsoil50(nso)

    IF (icycle<2) THEN
      casapool%nplant(npt,:) = casapool%rationcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0
    ENDIF

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
      casapool%pplant(npt,:) = casapool%ratiopcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
    ENDIF
  ENDDO

    CALL unpack_cnppool( land_pts, nsurft, L_tile_pts, casamet, casapool, casabal, phen,   &
                         Cpool_tile, Npool_tile, Ppool_tile, gLAI, phenphase )
                        


END SUBROUTINE casa_poolout_unpk

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE unpack_cnppool( land_pts, nsurft, L_tile_pts, casamet, casapool, casabal, phen, &
                           Cpool_tile, Npool_tile, Ppool_tile, gLAI, phenphase )

  USE cable_def_types_mod
  USE casadimension        ! icycle,mplant,mlitter,msoil
USE phenology_type_mod,    ONLY: phenology_type
USE casa_pool_type_mod,    ONLY: casa_pool    => casa_pool_type
USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE casa_met_type_mod,     ONLY: casa_met     => casa_met_type
USE casa_balance_type_mod, ONLY: casa_balance => casa_bal_type             
USE progs_cnp_vars_mod,    ONLY: nCpool_casa, nNpool_casa, nPPool_casa

 IMPLICIT NONE

! passed in
INTEGER, INTENT(IN) :: land_pts  
INTEGER, INTENT(IN) :: nsurft 
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft )
TYPE (casa_met)                  :: casamet
TYPE (casa_pool)                 :: casapool
TYPE (casa_balance)              :: casabal
TYPE (phenology_type)            :: phen

REAL, INTENT(INOUT) :: Cpool_tile( land_pts, nsurft, nCpool_casa )
REAL, INTENT(INOUT) :: Npool_tile( land_pts, nsurft, nNpool_casa )
REAL, INTENT(INOUT) :: Ppool_tile( land_pts, nsurft, nPPool_casa )
REAL, INTENT(INOUT) :: gLAI       ( land_pts, nsurft)
REAL, INTENT(INOUT) :: phenphase  ( land_pts, nsurft)

! local vars
INTEGER k
REAL(r_2) :: miss   = 0.0
REAL, DIMENSION( land_pts, nsurft)         :: clab,nsmin,pslab,&
                                                    pssorb,psocc,sumcbal,&
                                                    sumnbal,sumpbal
REAL  :: Cpl( land_pts, nsurft, mplant )
REAL  :: Npl( land_pts, nsurft, mplant )
REAL  :: Ppl( land_pts, nsurft, mplant )

REAL  :: Clit( land_pts, nsurft, mlitter )
REAL  :: Nlit( land_pts, nsurft, mlitter )
REAL  :: Plit( land_pts, nsurft, mlitter )

REAL  :: Cs( land_pts, nsurft, msoil )
REAL  :: Ns( land_pts, nsurft, msoil )
REAL  :: Ps( land_pts, nsurft, msoil )

! initialise variables
    clab(:,:)  = 0.0 
    cpl(:,:,:) = 0.0; clit(:,:,:) = 0.0; cs(:,:,:) = 0.0
    npl(:,:,:) = 0.0; nlit(:,:,:) = 0.0; ns(:,:,:) = 0.0  
    nsmin(:,:) = 0.0
    ppl(:,:,:) = 0.0; plit(:,:,:) = 0.0; ps(:,:,:) = 0.0  
    pslab(:,:) = 0.0; pssorb(:,:) = 0.0; psocc(:,:)= 0.0
    sumcbal(:,:)  = 0.0; sumnbal(:,:)   = 0.0; sumpbal(:,:) = 0.0

! unpack variables

    GLAI      = unpack(casamet%glai,l_tile_pts,miss)
    PHENPHASE = unpack(REAL(phen%phase),l_tile_pts,miss)
!   PHENPHASE = unpack(phen%phase,l_tile_pts,i_miss)

    do k= 1,3
     cpl(:,:,k)  = unpack(casapool%cplant(:,k) ,l_tile_pts,miss)
     clit(:,:,k) = unpack(casapool%clitter(:,k),l_tile_pts,miss)
     cs(:,:,k)   = unpack(casapool%csoil(:,k)  ,l_tile_pts,miss)
     npl(:,:,k)  = unpack(casapool%nplant(:,k) ,l_tile_pts,miss)
     nlit(:,:,k) = unpack(casapool%nlitter(:,k),l_tile_pts,miss)
     ns(:,:,k)   = unpack(casapool%nsoil(:,k)  ,l_tile_pts,miss)
     ppl(:,:,k)  = unpack(casapool%pplant(:,k) ,l_tile_pts,miss)
     plit(:,:,k) = unpack(casapool%plitter(:,k),l_tile_pts,miss)
     ps(:,:,k)   = unpack(casapool%psoil(:,k)  ,l_tile_pts,miss)
    enddo

    clab   = unpack(casapool%clabile     ,l_tile_pts,miss)
    nsmin  = unpack(casapool%nsoilmin    ,l_tile_pts,miss)
    pslab  = unpack(casapool%psoillab    ,l_tile_pts,miss)
    pssorb = unpack(casapool%psoilsorb   ,l_tile_pts,miss)
    psocc  = unpack(casapool%psoilocc    ,l_tile_pts,miss)
    sumcbal   = unpack(casabal%sumcbal   ,l_tile_pts,miss)
    sumnbal   = unpack(casabal%sumnbal   ,l_tile_pts,miss)
    sumpbal   = unpack(casabal%sumpbal   ,l_tile_pts,miss)

! prognostics
    cpool_tile(:,:,1)  =  clab(:,:)    
    cpool_tile(:,:,2)  =  cpl(:,:,1)   
    cpool_tile(:,:,3)  =  cpl(:,:,2)   
    cpool_tile(:,:,4)  =  cpl(:,:,3)   
    cpool_tile(:,:,5)  =  clit(:,:,1)  
    cpool_tile(:,:,6)  =  clit(:,:,2)  
    cpool_tile(:,:,7)  =  clit(:,:,3)  
    cpool_tile(:,:,8)  =  cs(:,:,1)    
    cpool_tile(:,:,9)  =  cs(:,:,2)    
    cpool_tile(:,:,10) =  cs(:,:,3)    
!   sumcbal(:,:)       =  0.0

!   if (icycle > 1) then
     npool_tile(:,:,1)  =  npl(:,:,1)   
     npool_tile(:,:,2)  =  npl(:,:,2)   
     npool_tile(:,:,3)  =  npl(:,:,3)   
     npool_tile(:,:,4)  =  nlit(:,:,1)  
     npool_tile(:,:,5)  =  nlit(:,:,2)  
     npool_tile(:,:,6)  =  nlit(:,:,3)  
     npool_tile(:,:,7)  =  ns(:,:,1)    
     npool_tile(:,:,8)  =  ns(:,:,2)    
     npool_tile(:,:,9)  =  ns(:,:,3)    
     npool_tile(:,:,10) =  nsmin(:,:)   
!    sumnbal(:,:)       =  0.0
!   endif

!   if (icycle > 2) then
     ppool_tile(:,:,1)  =  ppl(:,:,1)   
     ppool_tile(:,:,2)  =  ppl(:,:,2)   
     ppool_tile(:,:,3)  =  ppl(:,:,3)   
     ppool_tile(:,:,4)  =  plit(:,:,1)  
     ppool_tile(:,:,5)  =  plit(:,:,2)  
     ppool_tile(:,:,6)  =  plit(:,:,3)  
     ppool_tile(:,:,7)  =  ps(:,:,1)    
     ppool_tile(:,:,8)  =  ps(:,:,2)    
     ppool_tile(:,:,9)  =  ps(:,:,3)    
     ppool_tile(:,:,10) =  pslab(:,:)   
     ppool_tile(:,:,11) =  pssorb(:,:)  
     ppool_tile(:,:,12) =  psocc(:,:)   
!    sumpbal(:,:)       =  0.0
!   endif

END SUBROUTINE unpack_cnppool

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE unpack_glai( land_pts, nsurft, L_tile_pts, casamet,phen,GLAI,PHENPHASE)

! Lest 15 Jan 2013 - casa_poolout

  USE cable_def_types_mod
  USE casadimension
USE phenology_type_mod, ONLY: phenology_type
USE casa_met_type_mod,  ONLY: casa_met => casa_met_type

 IMPLICIT NONE

! passed in
INTEGER,               INTENT(IN)    :: land_pts          ! # land points being processed
INTEGER,               INTENT(IN)    :: nsurft 
LOGICAL,               INTENT(IN)    :: L_tile_pts(land_pts, nsurft )  ! TRUE if active tile
TYPE (casa_met),       INTENT(INOUT) :: casamet
TYPE (phenology_type), INTENT(INOUT) :: phen
REAL,                  INTENT(INOUT) :: GLAI    ( land_pts, nsurft )
REAL,                  INTENT(INOUT) ::PHENPHASE( land_pts, nsurft )

! local vars
REAL(r_2) :: miss   = 0.0

! unpack variables
gLAI      = UNPACK( casamet%glai, l_tile_pts, miss )
phenphase = UNPACK( REAL(phen%phase), l_tile_pts, miss )

END SUBROUTINE unpack_glai

!========================================================================
!========================================================================
!========================================================================

   !--- Lestevens 23Nov11: Based on Jhan's um2cable_lp but for Integers.
   !--- UM met forcing vars needed by CABLE which have UM dimensions
   !---(land_points)[_lp], which is no good to cable. These have to be
   !--- re-packed in a single vector of active tiles. Hence we use
   !--- conditional "mask" l_tile_pts(land_pts,nsurft) which is .true.
   !--- if the land point is/has an active tile

SUBROUTINE um2cable_ilp( mp, land_pts, nsurft, surft_pts, surft_index, l_tile_pts, umvar,  &
                         defaultin, cablevar, soiltype, skip )

IMPLICIT NONE

INTEGER, INTENT(IN)    :: mp 
INTEGER, INTENT(IN)    :: land_pts          ! # land points being processed
INTEGER, INTENT(IN)    :: nsurft 
INTEGER, INTENT(IN)    :: surft_pts(nsurft) ! # land points per tile
INTEGER, INTENT(IN)    :: surft_index(land_pts, nsurft) ! land_pt index of point
LOGICAL, INTENT(IN)    :: L_tile_pts(land_pts, nsurft )  ! TRUE if active tile
INTEGER, INTENT(IN)    :: umvar( land_pts )
INTEGER, INTENT(IN)    :: defaultin(10)
INTEGER, INTENT(INOUT) :: cablevar(mp)
INTEGER, INTENT(INOUT) :: soiltype(mp)
INTEGER, ALLOCATABLE   :: fvar(:,:)
LOGICAL, OPTIONAL      :: skip
INTEGER                :: n,k,l,i

ALLOCATE( fvar( land_pts, nsurft) )
fvar = 0.0

DO n=1, nsurft
  DO k=1, surft_pts(n)
    l = surft_index(k,n)
    fvar(l,n) = umvar(l)
    IF(.NOT. PRESENT(skip) ) THEN
       IF( n == nsurft ) THEN
          fvar(l,n) =  defaultin(9)
       ENDIF
    ENDIF
  ENDDO
ENDDO

cablevar     =  PACK(fvar,l_tile_pts)

IF(.NOT. PRESENT(skip) ) THEN
  DO i=1,mp
    IF(soiltype(i)==9) cablevar(i) =  defaultin(9)
  ENDDO
ENDIF

DEALLOCATE(fvar)

RETURN
END SUBROUTINE um2cable_ilp

SUBROUTINE redistr_luc( land_pts, nsurft, tile_frac, prev_yr_sfrac,inVar,outVar)

IMPLICIT NONE

INTEGER, INTENT(IN)    :: land_pts          ! # land points being processed
INTEGER, INTENT(IN)    :: nsurft 
REAL,    INTENT(IN)    :: tile_frac(land_pts, nsurft)
REAL,    INTENT(IN)    :: prev_yr_sfrac ( land_pts, nsurft )
REAL,    INTENT(IN)    :: inVar         ( land_pts, nsurft )
REAL,    INTENT(INOUT) :: outVar        ( land_pts, nsurft )
! local variables
REAL    :: dfrac( nsurft )
REAL    :: tmpVar, Rcount
INTEGER :: L,N
  
    DO L = 1, land_pts
      dfrac(:) = tile_frac(L,:) - prev_yr_sfrac(L,:)
      ! Collect all cut out areas from various decreasing tiles for averaging
      tmpVar = 0.0
      Rcount = 0.0
      DO N = 1, nsurft
        IF (dfrac(N) < 0.0) THEN
          tmpVar = tmpVar + ABS(dfrac(N)) * inVar(L,N)
          Rcount = Rcount + ABS(dfrac(N))
        ENDIF
      Enddo
      if (Rcount > 0) tmpVar = tmpVar / Rcount
      
      ! Add the averaged amount to those increasing tiles
      DO N = 1, nsurft
        IF (dfrac(N) > 0.0) THEN   ! those tiles increasing in size
          outVar(L,N) = (dfrac(N) * tmpVar + &
                        prev_yr_sfrac(L,N) * inVar(L,N)) / tile_frac(L,N)
        ELSE   ! those that are decreasing or no change in size
          outVar(L,N) = inVar(L,N)
        ENDIF
      Enddo
    Enddo
  
  END SUBROUTINE redistr_luc

SUBROUTINE casa_ndep_pk( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,        &
                         surft_pts, surft_index, L_tile_pts, Nidep, SoilType,  &
                         casaflux )

USE casa_flux_type_mod,    ONLY: casa_flux    => casa_flux_type
USE cable_pack_mod, ONLY: pack_landpts2mp_ICE

IMPLICIT NONE

! passed in
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: nsoil_max
INTEGER, INTENT(IN) :: ICE_SoilType
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: SoilType(mp)    ! soil type per tile
REAL,    INTENT(IN) :: nidep( land_pts )
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)  ! TRUE if active tile
TYPE (casa_flux),   INTENT(INOUT) :: casaflux

! pack variables gN/m2/day
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, Nidep, SoilType, &
                          Nidep(1:nsoil_max), casaflux % Nmindep )

END SUBROUTINE casa_ndep_pk

END MODULE casa_um_inout_mod

