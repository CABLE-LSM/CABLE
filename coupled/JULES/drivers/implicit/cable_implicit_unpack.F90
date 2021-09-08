MODULE cable_implicit_unpack_mod
  
CONTAINS

SUBROUTINE Implicit_unpack( cycleno, & ! nucycles
                            row_length,rows, land_pts, ntiles,L_tile_pts,     &
                            npft, sm_levels,                                  &
                            dim_cs1, dim_cs2,                                 &
                            tsoil, tsoil_tile, smcl, smcl_tile,               &
                            smvcst, sthf, sthf_tile, sthu, sthu_tile,         &
                            snow_tile, snow_rho1l ,isnow_flg3l, snow_depth3l, &
                            snow_mass3l, snow_rho3l, snow_tmp3l, snow_cond,   &
                            ftl_1, ftl_tile, fqw_1,  fqw_tile, tstar_tile,    &
                            surf_ht_flux_land, ecan_tile,                     &
                            esoil_tile, ei_tile, radnet_tile, tot_alb,        &
                            snow_age, canopy_tile, gs, gs_tile, t1p5m_tile,   &
                            q1p5m_tile, canopy_gb, fland, melt_tile,          &
                            !NPP, NPP_FT, GPP, GPP_FT, RESP_S,                  &
                            !RESP_S_TOT, RESP_P, RESP_P_FT, G_LEAF,&
                            !TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC,            &
                            surf_htf_tile, dtrad, dtstar_tile,                &
air, met, rad, rough, canopy,                                                 &
ssnow, bgc, bal, sum_flux, veg,                                               &
soil)

  !processor number, timestep number / width, endstep
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY: cable_runtime
USE cable_phys_constants_mod, ONLY : ctfrz => tfrz
  
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
USE cable_params_mod,         ONLY: soil_parameter_type

USE cable_def_types_mod, ONLY: mp
USE cable_um_tech_mod,   ONLY: um1 

IMPLICIT NONE
        
!___ re-decl input args
INTEGER :: cycleno
INTEGER :: row_length,rows, land_pts, ntiles, npft, sm_levels
INTEGER :: dim_cs1, dim_cs2 

LOGICAL,DIMENSION(land_pts, ntiles) ::                                        &
  L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0
  
REAL, DIMENSION(land_pts) ::                                                  &
  gs,         &  ! OUT "Stomatal" conductance to
  smvcst,     &  ! IN Volumetric saturation point
  fland          ! IN Land fraction on land tiles
   
REAL, DIMENSION(row_length,rows) ::                                           &
  !--- Net downward heat flux at surface over land.
  !--- fraction of gridbox (W/m2).
  surf_ht_flux_land,                                                          &
  !--- Moisture flux between layers. (kg/m^2/sec).
  !--- FQW(,1) is total water flux from surface, 'E'.
  fqw_1,                                                                      &
  !--- FTL(,K) =net turbulent sensible heat flux into layer K
  !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
  ftl_1         

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  surf_htf_tile,                                                              &
  !___Surface FTL, FQL for land tiles
  ftl_tile, fqw_tile,                                                         &
  !___(tiled) latent heat flux, melting, stomatatal conductance
  le_tile, melt_tile, gs_tile,                                                &
  radnet_tile, & ! INOUT Surface net radiation on tiles (W/m2)
  tot_alb,     & ! total albedo
  ei_tile,     & ! OUT EI for land tiles.
  ecan_tile,   & ! OUT ECAN for snow-free land tiles
  esoil_tile,  & ! evapotranspiration from soil moisture store (kg/m2/s) 
  resp_p_ft,   & !
  g_leaf,      & !
  gpp_ft,      & !
  npp_ft,      & !
  npp_ft_acc,    & ! sresp for CASA-CNP
  resp_w_ft_acc, & ! presp for CASA-CNP
  snow_tile,     & !
  snow_rho1l,    & ! Mean snow density
  snow_age,      & !
  canopy_tile,   & !
  t1p5m_tile,                                                                 &
  q1p5m_tile,                                                                 &
  tstar_tile,                                                                 &
  resp_s_tile,                                                                &
  transp_tile,                                                                &
  smgw_tile,                                                                  &
  dtstar_tile      !change in tstar_tile over time step

REAL, DIMENSION(land_pts) ::                                                  &
  resp_p,     & ! 
  npp,        & !
  gpp,        & !
  snow_grd,                                                                   &
  canopy_gb,                                                                  &
  t1p5m

REAL :: dtrad(mp)  ! CABLE change in rad%trad over time step

REAL, DIMENSION(land_pts,ntiles,3) ::                                         &
  snow_depth3l,                                                               &
  snow_mass3l,                                                                &
  snow_rho3l,                                                                 &
  snow_tmp3l,                                                                 &
  snow_cond 

!___flag for 3 layer snow pack
INTEGER :: isnow_flg3l(land_pts,ntiles)

!___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
!___ runoff ??
REAL, DIMENSION(land_pts,sm_levels) ::                                        &
  smcl,       & !
  sthf,                                                                       &
  sthu,                                                                       &
  tsoil       

!___(tiled) soil prognostics: as above 
REAL, DIMENSION(land_pts,ntiles,sm_levels) ::                                 &
  smcl_tile,                                                                  &
  sthu_tile,                                                                  &
  tsoil_tile,                                                                 &
  sthf_tile  

!___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
REAL, DIMENSION(land_pts,ntiles) ::                                           &
  nee_tile

REAL ::                                                                       &
  resp_s(land_pts,dim_cs1),     & !
  RESP_S_old(land_pts,dim_cs1), & !
  resp_s_tot(dim_cs2)             !
  
!___ local vars
INTEGER :: i,j,l,k,n,m

REAL, DIMENSION(mp) ::                                                        &
  fe_dlh,    & !
  fes_dlh,   & !
  fev_dlh      !

REAL, DIMENSION(land_pts,ntiles) ::                                           &
  !--- Local buffer surface FTL, FQL @ prev dt
  FTL_TILE_old, FQW_TILE_old,                                                 &
  lpts_ntiles

TYPE(air_type),      INTENT(INOUT)  :: air
TYPE(met_type),      INTENT(INOUT)  :: met
TYPE(radiation_type),INTENT(INOUT)  :: rad
TYPE(roughness_type),INTENT(INOUT)  :: rough
TYPE(canopy_type),   INTENT(INOUT)  :: canopy
TYPE(soil_snow_type),INTENT(INOUT)  :: ssnow
TYPE(bgc_pool_type), INTENT(INOUT)  :: bgc
TYPE(balances_type), INTENT(INOUT)  :: bal
TYPE(sum_flux_type), INTENT(INOUT)  :: sum_flux
TYPE(veg_parameter_type), INTENT(INOUT) :: veg
TYPE(soil_parameter_type),INTENT(INOUT) ::  soil
INTEGER:: i_miss = 0
REAL :: miss = 0.0

! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_implicit_unpack"

!--- set UM vars to zero
smcl_tile = 0.0; sthf_tile = 0.0; sthu_tile = 0.0
tsoil_tile = 0.0; smgw_tile = 0.0

DO j = 1,sm_levels
    
  tsoil_tile(:,:,j)= UNPACK(ssnow%tgg(:,j), l_tile_pts, miss)
  !liquid mass first
  smcl_tile(:,:,j)= UNPACK(REAL(ssnow%wbliq(:,j)), l_tile_pts, miss)
  smcl_tile(:,:,j) = smcl_tile(:,:,j) * soil%zse(j) * um1%rho_water
  !ice volumetric
  sthf_tile(:,:,j)= UNPACK(REAL(ssnow%wbice(:,j)), l_tile_pts, miss)
    
  !calcualte sthu_tilebefore smcl_tile incoudes ice mass
  DO n = 1,ntiles
      
    DO k = 1,um1%tile_pts(n)
         
      i = um1%tile_index(k,n)

      ! Exclude permanent ice 
      IF ( smvcst(i) > 0.0 ) & !liq mass relaative to max
         sthu_tile(i,n,j)= MAX( 0.0, smcl_tile(i,n,j) /                       &
                           (soil%zse(j) * smvcst(i) * um1%rho_water) )
         
      !add ice mass to liq mass
      smcl_tile(i,n,j) = smcl_tile(i,n,j) +                                   &
                         sthf_tile(i,n,j) * soil%zse(j) * um1%rho_ice
       !relative ice vol 
      IF ( smvcst(i) > 0.0 )                                                  &
          sthf_tile(i,n,j)= sthf_tile(i,n,j) / smvcst(i)

    END DO ! TILE_PTS(N)

  END DO ! NTILES

  smcl(:,j) = SUM(um1%tile_frac * smcl_tile(:,:,j),2)
  tsoil(:,j) = SUM(um1%tile_frac * tsoil_tile(:,:,j),2)

  sthf(:,j) = SUM(um1%tile_frac * sthf_tile(:,:,j),2)
  sthu(:,j) = SUM(um1%tile_frac * sthu_tile(:,:,j),2)

END DO ! SM_LEVELS

smgw_tile(:,:) = UNPACK(ssnow%GWwb(:), l_tile_pts, miss)



!--- unpack snow vars 
snow_rho1l  = UNPACK(ssnow%ssdnn, l_tile_pts, miss)
isnow_flg3l = UNPACK(ssnow%isflag, l_tile_pts, i_miss)
melt_tile   = UNPACK(ssnow%smelt, l_tile_pts, miss)
snow_tile= UNPACK(ssnow%snowd, l_tile_pts, miss)
snow_grd=  SUM(um1%tile_frac * snow_tile,2)  ! gridbox snow mass & snow below canopy 
  
!--- unpack layered snow vars 
DO k = 1,3
  snow_tmp3l(:,:,k) = UNPACK(ssnow%tggsn(:,k), l_tile_pts, miss)
  snow_mass3l(:,:,k)= UNPACK(ssnow%smass(:,k), l_tile_pts, miss)
  snow_rho3l(:,:,k) = UNPACK(ssnow%ssdn(:,k), l_tile_pts, miss)
  snow_cond(:,:,k)  = UNPACK(ssnow%sconds(:,k),l_tile_pts,miss)
  snow_depth3l(:,:,k)  = UNPACK(ssnow%sdepth(:,k),l_tile_pts,miss)
END DO

      
canopy%gswx_T = canopy%gswx_T / air%cmolar
gs_tile = UNPACK(canopy%gswx_T,l_tile_pts,miss)
gs =  SUM(um1%tile_frac * gs_tile,2)

!---preserve fluxes from the previous time step for the coastal grids
FTL_TILE_old = ftl_tile
FQW_TILE_old = fqw_tile
  
ftl_tile = UNPACK(canopy%fh,  l_tile_pts, miss)
fes_dlh = canopy%fes / (air%rlam * ssnow%cls)
fev_dlh = canopy%fev / air%rlam
fe_dlh =  fev_dlh + fes_dlh

fqw_tile      = UNPACK(fe_dlh, l_tile_pts, miss)
tstar_tile    = UNPACK(rad%trad, l_tile_pts, miss)
radnet_tile   = UNPACK( canopy%rnet , l_tile_pts, miss)
tot_alb        = UNPACK(rad%albedo_T,l_tile_pts, miss) 
ecan_tile     = UNPACK(fev_dlh,  l_tile_pts, miss)
!ESOIL_TILE    = UNPACK(fes_dlh, L_TILE_PTS, miss)
surf_htf_tile = UNPACK(canopy%ga,l_tile_pts,miss)
  
!EI_TILE = 0.

!Jun 2018 - need to split %fes into evaporation and sublimation
fes_dlh = 0.0
WHERE (ssnow%cls == 1.0)  fes_dlh = canopy%fes / (air%rlam * ssnow%cls)
esoil_tile = UNPACK(fes_dlh, l_tile_pts, miss)
fes_dlh = 0.0
WHERE (ssnow%cls == 1.1335) fes_dlh = canopy%fes / (air%rlam * ssnow%cls)
ei_tile = UNPACK(fes_dlh, l_tile_pts, miss)

!Jun 2018 - dtstar_tile unpacking
dtstar_tile = UNPACK(dtrad, l_tile_pts, miss)

snow_age = UNPACK(ssnow%snage, l_tile_pts, miss) 

transp_tile = UNPACK(canopy%fevc, l_tile_pts, miss) 

!unpack screen level (1.5m) variables - Convert back to K 
t1p5m_tile     = UNPACK(canopy%tscrn + ctfrz, l_tile_pts, miss)
q1p5m_tile     = UNPACK(canopy%qscrn, l_tile_pts, miss)

canopy_tile    = UNPACK(canopy%cansto, l_tile_pts, miss)
canopy_gb      = SUM(um1%tile_frac * canopy_tile,2)

!initialse full land grids and retain coastal grid fluxes
DO n = 1,ntiles

  DO k = 1,um1%tile_pts(n)
  
    l = um1%tile_index(k,n)
    j=(um1%land_index(l) - 1) / row_length + 1
    i = um1%land_index(l) - (j-1) * row_length
      
    IF ( fland(l) == 1.0) THEN 
      ftl_1(i,j) =  0.0
      fqw_1(i,j) =  0.0
    ELSE
      !retain sea/ice contribution and remove land contribution
      ftl_1(i,j) = ftl_1(i,j) - fland(l) * um1%tile_frac(l,n) *               &
                   FTL_TILE_old(l,n)
      fqw_1(i,j) = fqw_1(i,j) - fland(l) * um1%tile_frac(l,n) *               &
                   FQW_TILE_old(l,n)
    END IF
      
    surf_ht_flux_land(i,j) = 0.0

  END DO !tile_pts(n)

END DO !ntiles

DO n = 1,ntiles

  DO k = 1,um1%tile_pts(n)

    l = um1%tile_index(k,n)
    j=(um1%land_index(l) - 1) / row_length + 1
    i = um1%land_index(l) - (j-1) * row_length
    ftl_1(i,j) = ftl_1(i,j) + fland(l) * um1%tile_frac(l,n) * ftl_tile(l,n)
    fqw_1(i,j) = fqw_1(i,j) + fland(l) * um1%tile_frac(l,n) * fqw_tile(l,n)
      
    !Jun 2018 SURF_HT_FLUX_LAND is averaged over land not grid cell
    !SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
    !                         FLAND(L)*um1%TILE_FRAC(L,N) *              &
    !                         SURF_HTF_TILE(L,N)
    surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j) +                         &
                             um1%tile_frac(l,n) * surf_htf_tile(l,n)

  END DO

END DO

! Initialise grid-cell carbon fields that are accumulated over tiles
resp_p = 0.0;  npp = 0.0;  gpp = 0.0;  resp_s = 0.0

resp_s_tile   = UNPACK(canopy%frs, l_tile_pts, miss)
nee_tile      = UNPACK(canopy%fnee, l_tile_pts, miss)
npp_ft        = UNPACK(canopy%fnpp, l_tile_pts, miss)
g_leaf        = UNPACK(canopy%frday,l_tile_pts, miss)
resp_p_ft     = UNPACK(canopy%frp, l_tile_pts, miss)

gpp_ft = UNPACK(canopy%fnpp + canopy%frp + canopy%frday,                      &
                l_tile_pts, miss)

! convert from CABLE units (gC/m2/s) to UM units (kgC/m2/s)
resp_s_tile = resp_s_tile * 1.0e-3
g_leaf      = g_leaf * 1.0e-3
npp_ft      = npp_ft * 1.0e-3
gpp_ft      = gpp_ft * 1.0e-3
resp_p_ft   = resp_p_ft * 1.0e-3

! If CASA-CNP used, plant and soil resp need to be passed into 
! variables that are dumped to restart, because CASA-CNP only run daily
npp_ft_acc = resp_s_tile
resp_w_ft_acc = resp_p_ft

DO n = 1,ntiles 
  DO k = 1,um1%tile_pts(n)
    l = um1%tile_index(k,n)
    npp(l) = npp(l) + fland(l) * um1%tile_frac(l,n) * npp_ft(l,n)
    gpp(l) = gpp(l) + fland(l) * um1%tile_frac(l,n) * gpp_ft(l,n)
    resp_p(l) = resp_p(l) + fland(l) * um1%tile_frac(l,n) * resp_p_ft(l,n)

    !loop for soil resp. although DIM_CS1=1 (not 1 for triffid)
    DO i = 1,dim_cs1
      resp_s(l,i) = resp_s(l,i) +                                             &
                    fland(l) * um1%tile_frac(l,n) * resp_s_tile(l,n)
    END DO
    resp_s_tot(l) = SUM(resp_s(l,:))
    t1p5m(l) = SUM(t1p5m_tile(l,:))
  END DO
END DO

RETURN

END SUBROUTINE Implicit_unpack

END MODULE cable_implicit_unpack_mod
