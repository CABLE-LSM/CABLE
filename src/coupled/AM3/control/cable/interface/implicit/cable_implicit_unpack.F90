MODULE cable_implicit_unpack_mod
  
CONTAINS

SUBROUTINE implicit_unpack( cycleno, row_length, rows, land_pts, nsurft, npft, &
                            sm_levels, dim_cs1, timestep, mp, nsnl, land_index,&
                            surft_pts, surft_index, tile_frac, l_tile_pts,     &
                            smvcst, tsoil, smcl, sthf, sthu, snow_surft,       &
                            ftl_1, ftl_tile, fqw_1, fqw_tile, le_surft,        &
                            tstar_tile, surf_ht_flux_land, ecan_tile,          &
                            esoil_tile, ei_tile, radnet_tile, canopy_tile,     &
                            gs, gs_tile, t1p5m_tile, q1p5m_tile, canopy_gb,    &
                            fland, melt_tile, npp, npp_ft, gpp, gpp_ft,        &
                            resp_s, resp_s_tot, resp_p, resp_p_ft, g_leaf,     &
                            npp_ft_acc, surf_htf_tile, dtstar_tile,            &
                            progs, work, rad, met, rough, canopy, veg, soil,   &
                            ssnow, bal, air, bgc, sum_flux )

USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type ! CABLE requires extra progs
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep

!processor number, timestep number / width, endstep
USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY : cable_runtime
USE cable_common_module!, ONLY : cable_runtime, cable_user, fudge_out,       &
                         !         L_fudge, ktau_gl
USE cable_def_types_mod, ONLY: air_type, bgc_pool_type, met_type,              &              
                               balances_type, radiation_type, roughness_type,  & 
                               soil_parameter_type, soil_snow_type,            &
                               sum_flux_type, veg_parameter_type, canopy_type
  
USE cable_phys_constants_mod,  ONLY: density_liq, density_ice, tfrz

IMPLICIT NONE
        
!___ re-decl input args
INTEGER, INTENT(IN) :: cycleno
INTEGER, INTENT(IN) :: row_length,rows, land_pts, nsurft, npft, sm_levels
INTEGER, INTENT(IN) :: dim_cs1
INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: nsnl 
REAL,    INTENT(IN) :: timestep

INTEGER, INTENT(IN) :: surft_pts(nsurft) ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
INTEGER, INTENT(IN) :: land_index(land_pts)  ! tangled cell index of land_pt
REAL,    INTENT(IN) :: tile_frac(land_pts, nsurft)
LOGICAL, INTENT(IN) :: l_tile_pts(:,:)

TYPE(progs_cbl_vars_type), INTENT(OUT)  :: progs
TYPE(work_vars_type), INTENT(INOUT)       :: work

REAL, INTENT(IN) :: smvcst(land_pts,sm_levels) ! IN Volumetric saturation point
                   
REAL, INTENT(IN)  :: fland(land_pts)         ! IN Land fraction on land tiles
REAL, INTENT(OUT) :: gs(land_pts)            ! OUT "Stomatal" conductance to
   
!--- FQW(,1) is total water flux from surface, 'E'.
!--- Moisture flux between layers. (kg/m^2/sec).
REAL, INTENT(OUT) :: fqw_1(row_length,rows)

!--- FTL(,K) =net turbulent sensible heat flux into layer K
!--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
REAL, INTENT(OUT) :: ftl_1(row_length,rows)
    
!--- Net downward heat flux at surface over land.
!--- fraction of gridbox (W/m2).
REAL :: SURF_HT_FLUX_LAND(ROW_LENGTH,ROWS)
    
REAL, INTENT(OUT) :: melt_tile(land_pts,nsurft)
REAL, INTENT(OUT) :: smcl(land_pts,sm_levels)
REAL, INTENT(OUT) :: sthf(land_pts,sm_levels)
REAL, INTENT(OUT) :: sthu(land_pts,sm_levels)
REAL, INTENT(OUT) :: tsoil(land_pts,sm_levels)
REAL, INTENT(OUT) :: snow_surft(land_pts,nsurft)
REAL, INTENT(OUT) :: le_surft( land_pts, nsurft )       ! latent heat flux

  REAL, DIMENSION(land_pts,nsurft) ::                                 &
    SURF_HTF_TILE,&
    !___Surface FTL, FQL for land tiles
    FTL_TILE, FQW_TILE,                 &  
    !___(tiled) latent heat flux, melting, stomatatal conductance
    GS_TILE,     &  
    RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
    EI_TILE,     & ! OUT EI for land tiles.
    ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
    ESOIL_TILE,  & ! evapotranspiration from soil moisture store (kg/m2/s) 
    RESP_P_FT,   & !
    G_LEAF,      & !
    GPP_FT,      & !
    NPP_FT,      & !
    NPP_FT_ACC,    & ! sresp for CASA-CNP
    CANOPY_TILE,   & !
    T1P5M_TILE,    &
    Q1P5M_TILE,    &
    TSTAR_TILE,    &
    RESP_S_TILE,   & 
    DTSTAR_TILE      !change in tstar_tile over time step

  REAL, dimension(land_pts) ::                                            &
    RESP_P,     & ! 
    NPP,        & !
    GPP,        & !
    SNOW_GRD,   &  
    CANOPY_GB,  &
    T1P5M

  !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
  REAL, dimension(land_pts,nsurft) ::                                 &
    NEE_TILE

  REAL ::                                                                     &
    RESP_S(LAND_PTS,DIM_CS1),     & !
    RESP_S_old(LAND_PTS,DIM_CS1), & !
    RESP_S_TOT(land_pts)             !
  
  !___ local vars
  REAL :: DTRAD(mp) ! CABLE change in rad%trad over time step

  INTEGER :: i,j,l,k,n,m
  !___(tiled) soil prognostics: as above 
  REAL :: smcl_tile(land_pts,nsurft,sm_levels)
  REAL :: sthu_tile(land_pts,nsurft,sm_levels)
  REAL :: sthf_tile(land_pts,nsurft,sm_levels)
  REAL :: tsoil_tile(land_pts,nsurft,sm_levels)
  REAL :: smcl_ln(land_pts,nsurft,sm_levels)
  REAL :: sthu_ln(land_pts,nsurft,sm_levels)
  REAL :: sthf_ln(land_pts,nsurft,sm_levels)
  REAL :: tsoil_ln(land_pts,nsurft,sm_levels)
  REAL :: SNOW_COND(land_pts,nsurft,nsnl)
  
  REAL :: TOT_ALB(land_pts,nsurft) ! total albedo
  REAL :: RESP_W_FT_ACC(land_pts,nsurft) ! presp for CASA-CNP
  REAL :: SURF_CAB_ROFF(land_pts,nsurft) 
  REAL :: canopy_through_UM(land_pts,nsurft) 
  REAL :: TOT_TFALL_TILE(land_pts,nsurft) 

TYPE(met_type)            :: met
TYPE(radiation_type)      :: rad
TYPE(roughness_type)      :: rough
TYPE(soil_snow_type)      :: ssnow 
TYPE(balances_type)       :: bal 
TYPE(canopy_type)         :: canopy
TYPE(air_type)            :: air
TYPE(bgc_pool_type)       :: bgc
TYPE(sum_flux_type)       :: sum_flux
TYPE(veg_parameter_type)  :: veg        ! vegetation parameters
TYPE(soil_parameter_type) :: soil      ! soil parameters

REAL :: fe_dlh(mp) 
REAL :: fes_dlh(mp) 
REAL :: fev_dlh(mp) 
    
!--- Local buffer surface FTL, FQL @ prev dt
REAL :: ftl_tile_old(land_pts,nsurft)
REAL :: fqw_tile_old(land_pts,nsurft)

INTEGER:: i_miss = 0
REAL :: miss = 0.0
INTEGER :: isnow_flg_cable(land_pts,nsurft)
REAL :: dmdA_liq
REAL :: dmdA_ice
REAL :: gpp_ft_mp(mp)  
 
!-------- Unique subroutine body -----------
!jhan: should fland be applied or not?  - maybe we shouldnt even be dealing with
!any i,j fields
smcl        = 0.0 
sthf        = 0.0
sthu        = 0.0
tsoil       = 0.0
smcl_tile   = 0.0 
sthf_tile   = 0.0
sthu_tile   = 0.0
tsoil_tile  = 0.0
snow_cond   = 0.0
!smgw_tile       = 0.0 

DO j = 1,sm_levels
    
  tsoil_tile(:,:,j)= UNPACK( ssnow%tgg(:,j), L_tile_pts, miss)
  !liquid mass first
  smcl_tile(:,:,j) = UNPACK( REAL( ssnow%wbliq(:,j)), L_tile_pts, miss)
  smcl_tile(:,:,j) = smcl_tile(:,:,j) * soil%zse(j) * density_liq  
  !ice volumetric
  sthf_tile(:,:,j) = UNPACK( REAL( ssnow%wbice(:,j)), L_tile_pts, miss)
    
ENDDO ! SM_LEVELS

DO j = 1,sm_levels
  !calcualte sthu_tilebefore smcl_tile incoudes ice mass
  DO N=1,Nsurft
    
    DO K=1,surft_PTS(N)
       
      I = surft_INDEX(K,N)

      ! Exclude permanent ice 
       IF ( SMVCST(I,j) > 0. ) & !liq mass relaative to max
          STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) /                      &
                            (soil%zse(j)*SMVCST(I,j)*density_liq  ) )
       
      !add ice mass to liq mass
      SMCL_TILE(I,N,j) = SMCL_TILE(I,N,j) +                                 &
                         STHF_TILE(I,N,j) * soil%zse(j) * density_ice
       !relative ice vol 
       IF ( SMVCST(I,j) > 0. ) &
           STHF_TILE(I,N,j)= STHF_TILE(I,N,j)/SMVCST(I,j)

    ENDDO ! TILE_PTS(N)

  ENDDO ! NTILES

  SMCL(:,j) = SUM(TILE_FRAC * SMCL_TILE(:,:,j),2)
  TSOIL(:,j) = SUM(TILE_FRAC * TSOIL_TILE(:,:,j),2)

  STHF(:,J) = SUM(TILE_FRAC * STHF_TILE(:,:,J),2)
  STHU(:,J) = SUM(TILE_FRAC * STHU_TILE(:,:,J),2)

ENDDO ! SM_LEVELS


progs%SoilMoisture_CABLE(:,:,:)   = smcl_tile(:,:,:)
progs%FrozenSoilFrac_CABLE(:,:,:) = sthf_tile(:,:,:)
progs%SoilTemp_CABLE(:,:,:)       = tsoil_tile(:,:,:)

isnow_flg_cable                     = UNPACK(ssnow%isflag, L_TILE_PTS, i_miss)
progs%OneLyrSnowDensity_CABLE       = UNPACK(ssnow%ssdnn, l_tile_pts, miss)
progs%ThreeLayerSnowFlag_CABLE(:,:) = REAL( isnow_flg_cable(:,:) )
  
!--- unpack layered snow vars 
DO k = 1,3
  progs%SnowTemp_CABLE(:,:,k)    = UNPACK(ssnow%tggsn(:,k), L_TILE_PTS, miss)
  progs%SnowMass_CABLE(:,:,k)    = UNPACK(ssnow%smass(:,k), L_TILE_PTS, miss)
  progs%SnowDensity_CABLE(:,:,k) = UNPACK(ssnow%ssdn(:,k), L_TILE_PTS, miss)
  progs%SnowDepth_CABLE(:,:,k)   = UNPACK(ssnow%sdepth(:,k),L_TILE_PTS,miss)
  progs%SnowAge_CABLE(:,:)       = UNPACK(ssnow%snage, L_TILE_PTS, miss) 
  !snow_cond(:,:,k)              = UNPACK(ssnow%sconds(:,k),L_TILE_PTS,miss)i !this
  !doesn't go anywhere
END DO
      
!--- unpack snow vars 
melt_tile       = UNPACK(ssnow%smelt, l_tile_pts, miss)
snow_surft      = UNPACK(ssnow%snowd, l_tile_pts, miss)
snow_grd        = SUM(tile_frac * snow_surft,2)  ! gridbox snow mass & snow below canopy 

work%lying_snow(:) = snow_grd(:)  !jh:this sis done twice AND differently 
canopy%gswx_T = canopy%gswx_T/air%cmolar
gs_tile = UNPACK(canopy%gswx_T,L_TILE_PTS,miss)
gs      = SUM(tile_frac * GS_TILE,2)

!---preserve fluxes from the explicit step for the coastal grids
ftl_tile_old = ftl_tile
fqw_tile_old = fqw_tile
  
fes_dlh   = canopy%fes / ( air%rlam * ssnow%cls )
fev_dlh   = canopy%fev / air%rlam
fe_dlh    = fev_dlh + fes_dlh

ftl_tile = UNPACK(canopy%fh,  l_tile_pts, miss)
le_surft = UNPACK(canopy%fe,  l_tile_pts, miss)
fqw_tile = UNPACK(fe_dlh, l_tile_pts, miss)

!retain sea/ice contribution and remove land contribution j
!anymore
DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(K,N)
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    
    IF( fland(l) == 1.0) THEN 
      ftl_1(i,j) =  0.0
      fqw_1(i,j) =  0.0
    ELSE
      ftl_1(i,j) = ftl_1(i,j) - ( fland(l)* tile_frac(l,n) * ftl_tile_old(l,n) )
      fqw_1(i,j) = fqw_1(i,j) - ( fland(l)* tile_frac(l,n) * fqw_tile_old(l,n) )
    ENDIF
   
  ENDDO !surft_pts(n)
ENDDO !nsurft
!update with this ftl
DO n=1,nsurft
  DO k=1,surft_pts(n)
    l = surft_index(K,N)
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    
    ftl_1(i,j) = ftl_1(i,j) + ( fland(l)* tile_frac(l,n) * ftl_tile(l,n) )
    fqw_1(i,j) = fqw_1(i,j) + ( fland(l)* tile_frac(l,n) * fqw_tile(l,n) )
   
  ENDDO !surft_pts(n)
ENDDO !nsurft
  
tstar_tile    = UNPACK(rad%trad, l_tile_pts, miss)
radnet_tile   = UNPACK( canopy%rnet , l_tile_pts, miss)
ecan_tile     = UNPACK(fev_dlh, L_TILE_PTS, miss)

! need to split %fes into evaporation and sublimation
fes_dlh = 0.
WHERE (ssnow%cls==1.)  fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
esoil_tile = UNPACK(fes_dlh, L_tile_pts, miss)
  
fes_dlh = 0.
WHERE (ssnow%cls==1.1335) fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
ei_tile = UNPACK(fes_dlh, L_TILE_PTS, miss)

dtrad = rad%trad - rad%otrad        
dtstar_tile = UNPACK(dtrad, L_tile_pts, miss)

!CM3 - internal diag only  
!CM3 TRANSP_TILE = UNPACK(canopy%fevc, L_TILE_PTS, miss) 
tot_alb       = UNPACK(rad%albedo_T,L_TILE_PTS, miss) 

!unpack screen level (1.5m) variables - Convert back to K 
t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, L_TILE_PTS, miss)
q1p5m_tile     = UNPACK(canopy%qscrn, L_TILE_PTS, miss)

canopy_tile    = UNPACK(canopy%cansto, L_TILE_PTS, miss)
canopy_gb      = SUM(tile_frac * canopy_tile,2) !fland?


!initialse full land grids and retain coastal grid fluxes
!initialse full land grids and retain coastal grid fluxes
   
surf_htf_tile = UNPACK(canopy%ga,L_TILE_PTS,miss)
surf_ht_flux_land(:,:) = 0.
!fland?
DO n=1,nsurft
  DO K=1,surft_pts(N)

    l = surft_index(K,N)
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    
    surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j) +                   &
                             tile_frac(l,n) * surf_htf_tile(l,n)

  ENDDO
ENDDO

! Initialise grid-cell carbon fields that are accumulated over tiles
resp_p = 0.;  npp = 0.;  gpp = 0.;  resp_s = 0.

resp_s_tile   = UNPACK(canopy%frs,  L_tile_pts, miss) !see ISSUE#51
nee_tile      = UNPACK(canopy%fnee, L_tile_pts, miss)
npp_ft        = UNPACK(canopy%fnpp, L_tile_pts, miss)
g_leaf        = UNPACK(canopy%frday,L_tile_pts, miss)
resp_p_ft     = UNPACK(canopy%frp,  L_tile_pts, miss)
gpp_ft_mp     = canopy%fnpp + canopy%frp + canopy%frday
gpp_ft        = UNPACK(gpp_ft_mp, L_tile_pts, miss)

! convert from CABLE units (gC/m2/s) to UM units (kgC/m2/s)
resp_s_tile = resp_s_tile * 1.e-3
g_leaf      = g_leaf      * 1.e-3
npp_ft      = npp_ft      * 1.e-3
gpp_ft      = gpp_ft      * 1.e-3
resp_p_ft   = resp_p_ft   * 1.e-3

! If CASA-CNP used, plant and soil resp need to be passed into 
! variables that are dumped to restart, because CASA-CNP only run daily
npp_ft_acc    = resp_s_tile !see ISSUE#51
resp_w_ft_acc = resp_p_ft !see ISSUE#51
!fland ?
DO n=1,nsurft 
  DO k=1,surft_pts(N)
    l = surft_index(K,N)
    npp(l)=npp(l)+fland(l)*tile_frac(l,n)*npp_ft(l,n)
    gpp(l)=gpp(l)+fland(l)*tile_frac(l,n)*gpp_ft(l,n)
    resp_p(l)=resp_p(l)+fland(l)*tile_frac(l,n)*resp_p_ft(l,n)

    !loop for soil resp. although DIM_CS1=1 (not 1 for triffid)
    DO i=1,dim_cs1
      resp_s(l,i) = resp_s(l,i) + &
                    fland(l)*tile_frac(l,n)*resp_s_tile(l,n)
    ENDDO
    resp_s_tot(l)=SUM(resp_s(l,:)) !in CM2 this doesnt actually do anything !see ISSUE#51
    t1p5m(L)=SUM(t1p5m_tile(L,:))
  ENDDO
ENDDO
!jhan?fland?
work%lying_snow = SUM(TILE_FRAC * snow_surft,2) !gridbox snow mass

surf_cab_roff   = UNPACK(ssnow%rnof1, L_tile_pts, miss)
work%surf_roff  = SUM(tile_frac * surf_cab_roff,2)

surf_cab_roff      = UNPACK(ssnow%rnof2, L_TILE_PTS, miss)
work%sub_surf_roff = SUM(tile_frac * surf_cab_roff,2)

! %through is /dels in UM app. for STASH output  
canopy_through_UM = UNPACK(canopy%through, L_tile_pts, miss)
tot_tfall_tile    = canopy_through_UM  / timestep
work%tot_tfall    = SUM(tile_frac * tot_tfall_tile,2)


RETURN

END SUBROUTINE Implicit_unpack

End module cable_implicit_unpack_mod
