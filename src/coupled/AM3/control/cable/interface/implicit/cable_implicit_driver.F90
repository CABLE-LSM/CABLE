MODULE Cable_implicit_driv_mod
  
CONTAINS

SUBROUTINE cable_implicit_driver( cycleno, numcycles, i_day_number,            &
                                  timestep, timestep_number, row_length,       & 
                                  rows, land_pts, nsurft, npft, sm_levels,     &
                                  dim_cs1, mp, nrb, land_index, surft_pts,     &
                                  surft_index, l_tile_pts, ls_rain, conv_rain, & 
                                  ls_snow, conv_snow, tl_1, qw_1, ftl_1,fqw_1, &
                                  dtl1_1, dqw1_1, ctctq1, rad, met, rough,     &
                                  canopy, veg, soil, ssnow, bal, air, bgc,     &
                                  sum_flux )
                          !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,                  &
                          !GLAI, PHENPHASE)
                          

!subrs called 
USE cable_pack_mod,           ONLY: cable_pack_rr
USE cable_cbm_module,         ONLY: cbm_impl
  
! data: TYPE definitions of passed asarguments
USE cable_def_types_mod, ONLY : met_type, radiation_type, veg_parameter_type,  &
                                soil_parameter_type, roughness_type,           &
                                canopy_type, soil_snow_type, balances_type,    &
                                air_type, bgc_pool_type, sum_flux_type,        &
                                climate_type

! data: Scalars
USE cable_def_types_mod,      ONLY: msn, ncs, ncp
USE cable_phys_constants_mod, ONLY: TFRZ, CAPP 
USE cable_common_module,      ONLY: l_casacnp, l_vcmaxFeedbk, knode_gl,        &
                                    ktau_gl, kend_gl
                                
USE cable_common_module, ONLY : cable_runtime, cable_user !jhan:have to sort this out for JAC
  
!jhan: Leave for reference
!USE casavariable
!USE phenvariable
!USE casa_types_mod
!USE casa_um_inout_mod
!use POP_TYPES, only : pop_type
!CM3-standaloneUSE river_inputs_mod,   ONLY: river_step

!CM3 - updating CABLE-JULES rivers
USE jules_rivers_mod, ONLY: nstep_rivers

IMPLICIT NONE

INTEGER, INTENT(IN) :: cycleno
INTEGER, INTENT(IN) :: numcycles
INTEGER, INTENT(IN) :: i_day_number  ! day of year (1:365) counter for CASA-CNP
REAL   , INTENT(IN) :: timestep ! timestep length  [s]
INTEGER, INTENT(IN) :: timestep_number     
INTEGER, INTENT(IN) :: row_length,rows, land_pts, nsurft, npft, sm_levels
INTEGER, INTENT(IN) :: dim_cs1
INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: nrb
  
INTEGER, INTENT(IN) :: surft_pts(nsurft) ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)  ! tangled cell index of land_pt
LOGICAL, INTENT(IN) :: l_tile_pts(:,:)

REAL, INTENT(IN) :: ls_rain(row_length,rows)   ! Large scale rain
REAL, INTENT(IN) :: ls_snow(row_length,rows)   ! Large scale snow
REAL, INTENT(IN) :: conv_rain(row_length,rows) ! Convective rain
REAL, INTENT(IN) :: conv_snow(row_length,rows) ! Convective snow
REAL, INTENT(IN) :: tl_1(row_length,rows)      !
REAL, INTENT(IN) :: qw_1(row_length,rows)      !
REAL, INTENT(IN) :: dtl1_1(row_length,rows)     ! Level 1 increment to T field 
REAL, INTENT(IN) :: dqw1_1(row_length,rows)     ! Level 1 increment to q field 
REAL, INTENT(IN) :: ctctq1(row_length,rows)    ! information needed for increment to T an q field   
REAL, INTENT(IN) :: ftl_1(row_length,rows)     ! sensible heat flux to layer 1 H.(W/m2)
REAL, INTENT(IN) :: fqw_1(row_length,rows)     ! Moisture flux to layer 1 (kg/m^2/sec).

TYPE(met_type),            INTENT(INOUT) :: met
TYPE(radiation_type),      INTENT(INOUT) :: rad
TYPE(roughness_type),      INTENT(INOUT) :: rough
TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow 
TYPE(balances_type),       INTENT(INOUT) :: bal 
TYPE(canopy_type),         INTENT(INOUT) :: canopy
TYPE(air_type),            INTENT(INOUT) :: air
TYPE(bgc_pool_type),       INTENT(INOUT) :: bgc
TYPE(sum_flux_type),       INTENT(INOUT) :: sum_flux
TYPE(veg_parameter_type),  INTENT(INOUT) :: veg        ! vegetation parameters
TYPE(soil_parameter_type), INTENT(INOUT) :: soil      ! soil parameters
TYPE (climate_type) :: climate     ! climate variables

!jhan: Leave for reference
!REAL, INTENT(OUT) :: CPOOL_TILE(LAND_PTS,nsurft,10)
!REAL, INTENT(OUT) :: NPOOL_TILE(LAND_PTS,nsurft,10)
!REAL, INTENT(OUT) :: PPOOL_TILE(LAND_PTS,nsurft,12)
!REAL, INTENT(OUT) :: GLAI(LAND_PTS,nsurft)
!REAL, INTENT(OUT) :: PHENPHASE(LAND_PTS,nsurft)
    
INTEGER ::     &
  ktauday,      & ! day counter for CASA-CNP
  idoy            ! day of year (1:365) counter for CASA-CNP
INTEGER, SAVE :: &
  kstart = 1

REAL ::    dtl_mp(mp)
REAL ::    dqw_mp(mp)
!Ticket 132 - need ctctq1, incoming values of ftl_1 and fqw_1 on tiles
REAL :: ctctq1_mp(mp)          ! UM boundary layer coefficient
REAL :: ftl1_mp(mp)            ! grid box averaged FTL
REAL :: fqw1_mp(mp)            ! gird box averaged FQW
  
  REAL, DIMENSION(LAND_PTS) ::                               &
    LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
    SUB_SURF_ROFF, & !
    SURF_ROFF,     & !
    TOT_TFALL        !
!this is NA for single-site standalone 
INTEGER, parameter :: river_step = 1

!___ local vars
!jhan: Leave for reference
!This is a quick fix. These can be organised through namelists
LOGICAL :: spinup=.FALSE., spinconv=.FALSE.,                   &
           dump_read=.FALSE., dump_write=.FALSE.
INTEGER :: loy=365, lalloc=0
  
!___ 1st call in RUN (!=ktau_gl -see below) 
REAL, ALLOCATABLE:: fwork(:,:,:)
REAL :: dummy_rr(row_length,rows)
   
LOGICAL, SAVE :: first_cable_call = .TRUE.
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_implicit_driver"

!Prog Bank copies all prognostic and other variables whose
!values need to be retain from UM timestep to UM timestep.
!NOTE that canopy%cansto is a prognostic variable but is handled
!differently through the canopy%oldcansto variable

!-------- Unique subroutine body -----------

dtl_mp = 0. ; dqw_mp = 0.

!--- All these subrs do is pack a CABLE var with a UM var.
!-------------------------------------------------------------------
!--- UM met forcing vars needed by CABLE which have UM dimensions
!---(rowlength,rows)[_rr], which is no good to CABLE. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,nsurft) which is .true.
!--- if the land point is/has an active tile
!--- generic format:
!--- cable_pack_rr( UM var, default value for snow tile, CABLE var, mask )
!--- where mask tells cable_pack_rr whether or not to use default value 
!--- for snow tile 
!-------------------------------------------------------------------

dummy_rr(:,:) = 0.0
dummy_rr(:,:) = ( ls_rain(:,:) + conv_rain(:,:) ) * timestep
CALL cable_pack_rr( met%precip, dummy_rr, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

dummy_rr(:,:) = 0.0
dummy_rr(:,:) = ( ls_snow(:,:) + conv_snow(:,:) ) * timestep
CALL cable_pack_rr( met%precip_sn, dummy_rr, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )
      
CALL cable_pack_rr( met%tk, TL_1, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )

CALL cable_pack_rr( met%qv, QW_1, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )
      
CALL cable_pack_rr( dtl_mp, dtl1_1, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )
      
CALL cable_pack_rr( dqw_mp, dqw1_1, mp, l_tile_pts, row_length,     &
                    rows, nsurft, land_pts, land_index, surft_pts, surft_index )
      
!Ticket #132 implementation 
!dtl_mp, dqw_mp found on tiles above - these are the corrected dtl_mp and dqw_mp 
IF (cable_user%l_revised_coupling) THEN
  
  CALL cable_pack_rr( ctctq1_mp, ctctq1, mp, l_tile_pts, row_length,     &
                      rows, nsurft, land_pts, land_index, surft_pts, surft_index )
  
  CALL cable_pack_rr( ftl1_mp, ftl_1, mp, l_tile_pts, row_length,     &
                      rows, nsurft, land_pts, land_index, surft_pts, surft_index )
  
  CALL cable_pack_rr( fqw1_mp, fqw_1, mp, l_tile_pts, row_length,     &
                      rows, nsurft, land_pts, land_index, surft_pts, surft_index )

  dtl_mp = dtl_mp - ctctq1_mp * ftl1_mp/CAPP  !NB FTL_1 is in W/m2 hence / CAPP
  dqw_mp = dqw_mp - ctctq1_mp * fqw1_mp

ENDIF

met%precip = met%precip + met%precip_sn
met%tk     = met%tk + dtl_mp
met%qv     = met%qv + dqw_mp
met%tvair  = met%tk
met%tvrad  = met%tk

CALL cbm_impl( cycleno, numcycles, mp, nrb, timestep_number, timestep,         &
               air, bgc, canopy, met, bal, rad, rough,                         &
               soil, ssnow, sum_flux, veg, climate )

! Integrate wb_lake over the river timestep.
! Used to scale river flow within ACCESS
! Zeroed each river step in subroutine cable_lakesriver and on restarts.
!  ssnow_wb_lake in kg/m^2
! CM3 - updated to use JULES7.x vars
if (cycleno == numcycles) THEN
  ssnow%totwblake = ssnow%totwblake + ssnow%wb_lake/nstep_rivers
end if
 
  ! Lestevens - temporary ?
  ktauday = int(24.0*3600.0/TIMESTEP)
  idoy=i_day_number
  
  !Jan 2018: Only call carbon cycle prognostics updates on the last call to 
  !cable_implicit per atmospheric time step
  if (cycleno==numcycles) then
    !Call CASA-CNP
    !CM2!if (l_casacnp) & 
    !CM2!  CALL bgcdriver(ktau_gl,kstart,kend_gl,timestep,met,ssnow,canopy,veg,soil, &
    !CM2!                 climate,casabiome,casapool,casaflux,casamet,casabal,phen, &
    !CM2!                 pop, spinConv,spinup, ktauday, idoy,loy, dump_read,   &
    !CM2!                 dump_write, LALLOC)

    !CM2!CALL sumcflux(ktau_gl,kstart,kend_gl,TIMESTEP,bgc,canopy,soil,ssnow,      &
    !CM2!              sum_flux,veg,met,casaflux,l_vcmaxFeedbk)
  endif

!GLAI = 0.0
!PHENPHASE = 0.0
  ! Only call carbon cycle prognostics updates on the last call to 
  ! cable_implicit per atmospheric time step
  ! Call CASA-CNP collect pools
  !block!if (cycleno==numcycles .AND. l_casacnp) & 
  !block!  CALL casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
  !block!                        CPOOL_TILE,NPOOL_TILE,PPOOL_TILE, &
  !block!                        GLAI,PHENPHASE)

  !-------- End Unique subroutine body -----------

!End Testing puroses:

RETURN

END SUBROUTINE cable_implicit_driver

END module cable_implicit_driv_mod
