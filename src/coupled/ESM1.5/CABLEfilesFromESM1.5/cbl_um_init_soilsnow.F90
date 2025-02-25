MODULE cbl_um_init_soilsnow_mod
   
IMPLICIT NONE

PUBLIC initialize_soilsnow

CONTAINS

! ESM1.5 args moderated by necessary AM3 USEs etc to support AM3 algorithm(s)
! ESM1.5 algorithm(s) super indented 
SUBROUTINE initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile,      &
                                snow_tile, snow_rho1l, snage_tile, isnow_flg3l,&
                                snow_rho3l, snow_cond, snow_depth3l,           &
                                snow_mass3l, snow_tmp3l, fland,                &
                                sin_theta_latitude ) 

! subrs
USE cable_init_wetfac_mod,   ONLY: initialize_wetfac
USE cable_um_init_subrs_mod, ONLY: um2cable_lp

USE cable_def_types_mod, ONLY: mp, msn
USE cable_um_tech_mod,   ONLY: um1, soil, ssnow, met, bal, veg
USE cable_common_module, ONLY: cable_runtime, cable_user, frozen_limit

USE cable_phys_constants_mod, ONLY: density_ice, density_liq
USE cable_phys_constants_mod, ONLY: TFRZ
USE cable_phys_constants_mod, ONLY: cgsnow, csice, cswat 
USE grid_constants_mod_cbl,   ONLY: ICE_SoilType
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: smvcst
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles, um1%sm_levels) ::    &
      sthf_tile, &   !
      smcl_tile, &   !
      tsoil_tile     !

   INTEGER, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) :: isnow_flg3l 

   REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles) :: snow_tile

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles) ::                    &
      snow_rho1l, &  !
      snage_tile     !

   REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles,3) :: snow_cond

   REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles,3) ::                  & 
      snow_rho3l,    & !
      snow_depth3l,  & !
      snow_mass3l,   & !
      snow_tmp3l       !
   
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: fland 
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: sin_theta_latitude
   
!local vars
INTEGER :: i,j,k,L,n
REAL    :: zsetot
REAL    :: ice_vol_tmp(um1%land_pts, um1%ntiles, um1%sm_levels)
REAL    :: wbtot_mp(mp, um1%sm_levels)
REAL    :: wbice_mp(mp, um1%sm_levels)
! from ESM1.5
REAL    :: max_snow_depth=50000. 
REAL    :: sfact(mp)
REAL    :: fvar(um1%land_pts )
LOGICAL :: skip =.TRUE. 

ssnow%pudsto       = 0.0 
ssnow%pudsmx       = 0.0
ssnow%wbtot        = 0.0
ssnow%totwblake    = 0.0  ! wb_lake integrated over river timestep
ssnow%tggav        = 0.0
ssnow%qhz          = 0.0
ssnow%qhlev        = 0.0
ssnow%qrecharge    = 0.0
ssnow%rtevap_sat   = 0.0
ssnow%rtevap_unsat = 0.0

! identify module parameters here and recast (NOT ssnow% state variables)
ssnow%t_snwlr      = 0.05
ssnow%rtsoil       = 50.0
ssnow%satfrac      = 0.5
ssnow%wtd          = 1.0

!from esm1.5 snow_tile    = MIN( max_snow_depth, snow_tile )
ssnow%snowd  = PACK( snow_tile, um1%L_TILE_PTS)
ssnow%snage  = PACK( snage_tile, um1%L_TILE_PTS)
ssnow%ssdnn  = PACK( snow_rho1l, um1%L_TILE_PTS)  
ssnow%isflag = PACK( isnow_flg3L, um1%L_TILE_PTS)  

! over snow layers
DO j=1, msn
  ssnow%sdepth(:,j)= PACK( snow_depth3l(:,:,j), um1%L_TILE_PTS)
  ssnow%smass(:,j) = PACK( snow_mass3l(:,:,j), um1%L_TILE_PTS)  
  ssnow%ssdn(:,j)  = PACK( snow_rho3l(:,:,j), um1%L_TILE_PTS)  
  ssnow%tggsn(:,j) = PACK( snow_tmp3l(:,:,j), um1%L_TILE_PTS)  
  ssnow%sconds(:,J)= PACK(SNOW_COND(:,:,J), um1%L_TILE_PTS) ! rm-ed in AM3
ENDDO 
     
! over soil layers
DO j=1, um1%sm_levels
  ssnow%tgg(:,j)   = PACK( tsoil_tile(:,:,j), um1%L_TILE_PTS)
ENDDO 

! AM3 this is broken into init/update
! IF( first_call) THEN 
        
ssnow%osnowd = ssnow%snowd

zsetot = sum(soil%zse)
DO k = 1, um1%sm_levels
  ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
END DO

!jhan:Ian: check this is indeed equiv to previous algorithm
ice_vol_tmp(:,:,:)  = 0.0
DO n=1,um1%ntiles                                                       
  DO k=1,um1%tile_pts(n)                                           
    i = um1%tile_index(k,n)                                      
    DO j = 1, um1%sm_levels
      ice_vol_tmp(i,n,j)  = sthf_tile(i,n,j) * smvcst(i)
    ENDDO ! J
  ENDDO
ENDDO
 
DO j = 1, um1%sm_levels
  !liq volume  from (tot_mass - ice_mass) / (dz*rho_liq)
  wbtot_mp(:,j)    = PACK( smcl_tile(:,:,j), um1%L_TILE_PTS)
  
  ssnow%wbice(:,J) = PACK( ice_vol_tmp(:,:,j), um1%L_TILE_PTS)
  ssnow%wbice(:,J) = MAX ( 0.0 , ssnow%wbice(:,j) )
  wbice_mp(:,j)    = ssnow%wbice(:,j) * ( soil%zse(j) * density_ice )
  
  ssnow%wbliq(:,j) = ( wbtot_mp(:,j) - wbice_mp(:,j) )  / ( soil%zse(j) * density_liq )
  ssnow%wb(:,j)    = ssnow%wbice(:,j) + ssnow%wbliq(:,j) 
END DO

         ! from ESM1.5
         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               L = um1%TILE_INDEX(K,N)
               fvar(L) = real(L)
            ENDDO
         ENDDO
         CALL um2cable_lp( fland, fland, ssnow%fland, soil%isoilm, skip )
         CALL um2cable_lp( fvar, fvar, ssnow%ifland, soil%isoilm, skip )
         
         ssnow%owetfac = MAX( 0., MIN( 1.0,                                    &
                         ( ssnow%wb(:,1) - soil%swilt ) /                      &
                         ( soil%sfc - soil%swilt) ) )

! wetfac initialized here as used to init owetfac on firstimestep in cbm
! add to startdump? includes Temporay fix for accounting for reduction of 
! soil evaporation due to freezing. Also includes specific lakes case  
! Prevents divide by zero at glaciated points where both wb and wbice=0.
CALL initialize_wetfac( mp, ssnow%wetfac, soil%swilt, soil%sfc,            &
                        ssnow%wb(:,1), ssnow%wbice(:,1), ssnow%snowd,      &
                        veg%iveg, met%tk, tfrz ) 

         ! Temporay fix for accounting for reduction of soil evaporation 
         ! due to freezing
         WHERE( ssnow%wbice(:,1) > 0. )                                        &
            ! Prevents divide by zero at glaciated points where both 
            ! wb and wbice=0.
            ssnow%owetfac = ssnow%owetfac * ( 1.0 - ssnow%wbice(:,1) /         &
                            ssnow%wb(:,1) )**2
      
! initialized here on first call 
ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1) 
 
!jhan: do we want to do this before %owetfac is set 
DO J = 1, um1%sm_levels
  !should be removed!!!!!!!! This cannot conserve if there are any
  !dynamics 
  WHERE( soil%isoilm == ICE_SoilType)
    ssnow%wb(:,j)    = 0.95 * soil%ssat
    ssnow%wbice(:,j) = frozen_limit  * ssnow%wb(:,j)
  ENDWHERE
  
  !no not force rho_water==rho_ice==1000.0
  ssnow%wbtot = ssnow%wbtot + soil%zse(j) *                       &
                             ( ssnow%wbliq(:,j) * density_liq +   &
                               ssnow%wbice(:,j) * density_ice )                     
ENDDO

!initialize %gammz (see soil_snow)
ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil                                        &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * density_liq                                     &
            & + ssnow%wbice(:,1) * csice * density_ice, soil%css * soil%rhosoil ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * cgsnow * ssnow%snowd
     
RETURN
END SUBROUTINE initialize_soilsnow

END MODULE cbl_um_init_soilsnow_mod

