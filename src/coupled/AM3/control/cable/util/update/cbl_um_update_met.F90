MODULE cbl_um_update_met_mod
   
 IMPLICIT NONE

CONTAINS

SUBROUTINE update_met( mp, row_length, rows, timestep, land_pts, nsurft,       &
                       surft_pts, surft_index, land_index, L_tile_pts, co2_mmr,&
                       ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar, met )
                       !block!CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,
                       !block!L_CO2_INTERACTIVE )   

USE cable_pack_mod,             ONLY: cable_pack_rr
USE cable_phys_constants_mod,   ONLY: UMIN 
USE cable_def_types_mod,        ONLY: met_type

INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
REAL,    INTENT(IN) :: timestep
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # points on each tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! index of tile points
INTEGER, INTENT(IN) :: land_index(land_pts)          ! index of land points
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)

! "forcing"
REAL,    INTENT(IN) :: co2_mmr
REAL,    INTENT(IN) :: ls_rain(row_length,rows)
REAL,    INTENT(IN) :: ls_snow(row_length,rows)
REAL,    INTENT(IN) :: tl_1(row_length,rows)
REAL,    INTENT(IN) :: qw_1(row_length,rows)
REAL,    INTENT(IN) :: vshr_land(row_length,rows)
REAL,    INTENT(IN) :: pstar(row_length,rows)
  
TYPE(met_type), INTENT(OUT) :: met

!local decs 
REAL    :: ls_rain_dt(row_length,rows)
REAL    :: ls_snow_dt(row_length,rows)
REAL    :: precip_dt(row_length,rows)
REAL    :: CO2_3D(row_length,rows)     ! co2 mass mixing ratio
LOGICAL :: L_CO2_INTERACTIVE = .FALSE. ! namelist?

!block!! rml 2/7/13 Extra atmospheric co2 variables
!block!   LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
!block!   INTEGER, INTENT(IN) :: CO2_DIM_LEN, CO2_DIM_ROW
!block!   REAL, INTENT(IN) :: CO2_3D(:,:)  ! co2 mass mixing ratio

met%DoY = 1.0 !jhan: fudged initialization to prevent dangling arg to init_rad

!jhan:implies precip is given as per second
ls_rain_dt = ls_rain * timestep
ls_snow_dt = ls_snow * timestep
precip_dt  = ls_rain_dt + ls_snow_dt
CALL cable_pack_rr( met%precip, precip_dt, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( met%precip_sn, ls_snow_dt, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( met%tk, tl_1, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( met%qv, qw_1, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( met%ua, vshr_land, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

met%ua = MAX( UMIN, met%ua ) !---this is necessary cloberring at present 

CALL cable_pack_rr( met%pmb, pstar, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

met%pmb = 0.01 * met%pmb  ! what is this .01? Units conversion?

!block!IF( .NOT. ALLOCATED( conv_rain_prevstep(mp) ) ) THEN
!block!  ALLOCATE( conv_rain_prevstep(mp )
!block!  ALLOCATE( conv_snow_prevstep(mp) )
!block!  conv_rain_prevstep = 0. 
!block!  conv_snow_prevstep = 0.
!block!ENDIF   
!block!we would need to get conv_*_from prev step for a start
!block!met%precip   =  (met%precip + conv_rain_prevstep) &
!block!               + (met%precip_sn +  conv_snow_prevstep)
!block!               + (met%precip_sn +  conv_rain_prevstep)

!jhan:This might be peculiar to the UM? and if so presents a problem
met%tvair = met%tk
met%tvrad = met%tk

!******************clobbered by a hard-wired number ********************!      
! rml 24/2/11 Set atmospheric CO2 seen by cable to CO2_MMR (value seen 
! by radiation scheme).  Option in future to have cable see interactive 
! (3d) CO2 field Convert CO2 from kg/kg to mol/mol ( m_air, 
! 28.966 taken from include/constant/ccarbon.h file )
! r935 rml 2/7/13 Add in co2_interactive option
IF (L_CO2_INTERACTIVE) THEN
  CALL cable_pack_rr( met%ca, CO2_3D, mp, l_tile_pts, row_length,              &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )
ELSE
  met%ca = CO2_MMR
ENDIF
!jhan:WTF
met%ca = met%ca * 28.966/44.

END SUBROUTINE update_met

END MODULE cbl_um_update_met_mod
