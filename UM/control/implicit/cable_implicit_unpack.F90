!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Updates CABLE variables (as altered by first pass through boundary 
!          layer and convection scheme), calls cbm, passes CABLE variables back 
!          to UM. 'Implicit' is the second call to cbm in each UM timestep.
!
! Called from: UM/JULES cable_implicit_main
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
! ==============================================================================

module cable_implicit_unpack_mod
  
contains

subroutine Implicit_unpack( cycleno, & ! nucycles
                            row_length,rows, land_pts, ntiles, npft, sm_levels,&
                            dim_cs1, dim_cs2,                                  &
                            TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,SMGW_TILE,      &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,                      &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                            SNOW_AGE, CANOPY_TILE, GS, gs_tile, T1P5M_TILE,    &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE,           &
                            NPP, NPP_FT, GPP, GPP_FT, RESP_S,                  &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,&
                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC,            &
                            SURF_HTF_TILE, DTRAD, DTSTAR_TILE )

  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir

  USE cable_diag_module  

  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
  USE cable_common_module!, ONLY : cable_runtime, cable_user, fudge_out,       &
                         !         L_fudge, ktau_gl
  
  USE cable_def_types_mod, ONLY : mp
  USE cable_data_module,   ONLY : PHYS
  USE cable_um_tech_mod,   ONLY : um1 ,canopy, rad, soil, ssnow, air,         &
                                  basic_diag, veg

  USE cable_decs_mod, ONLY : L_tile_pts!, rho_water

  implicit none
        
  !___ re-decl input args
  integer :: cycleno
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
  integer :: dim_cs1, dim_cs2 

  REAL, DIMENSION(land_pts) ::                                            &
    GS,         &  ! OUT "Stomatal" conductance to
    SMVCST,     &  ! IN Volumetric saturation point
    FLAND          ! IN Land fraction on land tiles
   
  real, dimension(ROW_LENGTH,ROWS) ::                                 &
    !--- Net downward heat flux at surface over land.
    !--- fraction of gridbox (W/m2).
    SURF_HT_FLUX_LAND,           &
    !--- Moisture flux between layers. (kg/m^2/sec).
    !--- FQW(,1) is total water flux from surface, 'E'.
    FQW_1,       &  
    !--- FTL(,K) =net turbulent sensible heat flux into layer K
    !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
    FTL_1         

  REAL, DIMENSION(land_pts,ntiles) ::                                 &
    SURF_HTF_TILE,&
    !___Surface FTL, FQL for land tiles
    FTL_TILE, FQW_TILE,                 &  
    !___(tiled) latent heat flux, melting, stomatatal conductance
    LE_TILE, MELT_TILE, GS_TILE,     &  
    RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
    TOT_ALB,     & ! total albedo
    EI_TILE,     & ! OUT EI for land tiles.
    ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
    ESOIL_TILE,  & ! evapotranspiration from soil moisture store (kg/m2/s) 
    RESP_P_FT,   & !
    G_LEAF,      & !
    GPP_FT,      & !
    NPP_FT,      & !
    NPP_FT_ACC,    & ! sresp for CASA-CNP
    RESP_W_FT_ACC, & ! presp for CASA-CNP
    SNOW_TILE,     & !
    SNOW_RHO1L,    & ! Mean snow density
    SNOW_AGE,      & !
    CANOPY_TILE,   & !
    T1P5M_TILE,    &
    Q1P5M_TILE,    &
    TSTAR_TILE,    &
    RESP_S_TILE,   & 
    TRANSP_TILE,   &
    SMGW_TILE,     &
    DTSTAR_TILE      !change in tstar_tile over time step

  REAL, dimension(land_pts) ::                                            &
    RESP_P,     & ! 
    NPP,        & !
    GPP,        & !
    SNOW_GRD,   &  
    CANOPY_GB,  &
    T1P5M,      &
    DTRAD         ! CABLE change in rad%trad over time step

  REAL, DIMENSION(land_pts,ntiles,3) ::                               &
    SNOW_DEPTH3L,  &
    SNOW_MASS3L,   &
    SNOW_RHO3L,    &
    SNOW_TMP3L,    &
    SNOW_COND 

  !___flag for 3 layer snow pack
  INTEGER :: ISNOW_FLG3L(LAND_PTS,NTILES)

  !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
  !___ runoff ??
  REAL, DIMENSION(land_pts,sm_levels) ::                              &
    SMCL,       & !
    STHF,       &
    STHU,       &
    TSOIL       

  !___(tiled) soil prognostics: as above 
  REAL, DIMENSION(land_pts,ntiles,sm_levels) ::                   &
    SMCL_TILE,  & 
    STHU_TILE,  &
    TSOIL_TILE, &
    STHF_TILE  

  !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
  REAL, dimension(land_pts,ntiles) ::                                 &
    NEE_TILE

  REAL ::                                                                     &
    RESP_S(LAND_PTS,DIM_CS1),     & !
    RESP_S_old(LAND_PTS,DIM_CS1), & !
    RESP_S_TOT(DIM_CS2)             !
  
  !___ local vars
  INTEGER :: i,j,l,k,n,m

  REAL, DIMENSION(mp) ::                                                                     &
    fe_dlh,    & !
    fes_dlh,   & !
    fev_dlh      !

  REAL, DIMENSION(land_pts,ntiles) ::                                 &
    !--- Local buffer surface FTL, FQL @ prev dt
    FTL_TILE_old, FQW_TILE_old, &
    lpts_ntiles

  INTEGER:: i_miss = 0
  REAL :: miss = 0.0
  
  REAL, POINTER :: TFRZ

  ! std template args 
  character(len=*), parameter :: subr_name = "cable_implicit_unpack"

# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------

  TFRZ => PHYS%TFRZ
  
  !--- set UM vars to zero
  SMCL_TILE = 0.; STHF_TILE = 0.; STHU_TILE = 0.
  TSOIL_TILE = 0.; SMGW_TILE = 0.

  DO j = 1,SM_LEVELS
    
    TSOIL_TILE(:,:,j)= UNPACK(ssnow%tgg(:,j), L_TILE_PTS, miss)
    !liquid mass first
    SMCL_TILE(:,:,j)= UNPACK(REAL(ssnow%wbliq(:,j)), L_TILE_PTS, miss)
    SMCL_TILE(:,:,j)=SMCL_TILE(:,:,j)*soil%zse(j)*um1%RHO_WATER
    !ice volumetric
    STHF_TILE(:,:,j)= UNPACK(REAL(ssnow%wbice(:,j)), L_TILE_PTS, miss)
    
    !calcualte sthu_tilebefore smcl_tile incoudes ice mass
    DO N=1,NTILES
      
      DO K=1,um1%TILE_PTS(N)
         
        I = um1%TILE_INDEX(K,N)

        ! Exclude permanent ice 
         IF ( SMVCST(I) > 0. ) & !liq mass relaative to max
            STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) /                      &
                              (soil%zse(j)*SMVCST(I)*um1%RHO_WATER) )
         
        !add ice mass to liq mass
        SMCL_TILE(I,N,j) = SMCL_TILE(I,N,j) +                                 &
                           STHF_TILE(I,N,j) * soil%zse(j) * um1%RHO_ICE
         !relative ice vol 
         IF ( SMVCST(I) > 0. ) &
             STHF_TILE(I,N,j)= STHF_TILE(I,N,j)/SMVCST(I)

      ENDDO ! TILE_PTS(N)

    ENDDO ! NTILES

    SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
    TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)

    STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
    STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)

  ENDDO ! SM_LEVELS

  SMGW_TILE(:,:) = UNPACK(ssnow%GWwb(:), L_TILE_PTS, miss)



  !--- unpack snow vars 
  SNOW_RHO1L  = UNPACK(ssnow%ssdnn, L_TILE_PTS, miss)
  ISNOW_FLG3L = UNPACK(ssnow%isflag, L_TILE_PTS, i_miss)
  MELT_TILE   = UNPACK(ssnow%smelt, L_TILE_PTS, miss)
  SNOW_TILE= UNPACK(ssnow%snowd, L_TILE_PTS, miss)
  SNOW_GRD=  SUM(um1%TILE_FRAC * SNOW_TILE,2)  ! gridbox snow mass & snow below canopy 
  
  !--- unpack layered snow vars 
  do k = 1,3
    SNOW_TMP3L(:,:,k) = UNPACK(ssnow%tggsn(:,k), L_TILE_PTS, miss)
    SNOW_MASS3L(:,:,k)= UNPACK(ssnow%smass(:,k), L_TILE_PTS, miss)
    SNOW_RHO3L(:,:,k) = UNPACK(ssnow%ssdn(:,k), L_TILE_PTS, miss)
    SNOW_COND(:,:,k)  = UNPACK(ssnow%sconds(:,k),L_TILE_PTS,miss)
    SNOW_DEPTH3L(:,:,k)  = UNPACK(ssnow%sdepth(:,k),L_TILE_PTS,miss)
  enddo

      
  canopy%gswx_T = canopy%gswx_T/air%cmolar
  GS_TILE = UNPACK(canopy%gswx_T,L_TILE_PTS,miss)
  GS =  SUM(um1%TILE_FRAC * GS_TILE,2)

  !---preserve fluxes from the previous time step for the coastal grids
  FTL_TILE_old = FTL_TILE
  FQW_TILE_old = FQW_TILE
  
  FTL_TILE = UNPACK(canopy%fh,  l_tile_pts, miss)
  fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
  fev_dlh = canopy%fev/air%rlam
  fe_dlh =  fev_dlh + fes_dlh

  FQW_TILE      = UNPACK(fe_dlh, l_tile_pts, miss)
  TSTAR_TILE    = UNPACK(rad%trad, l_tile_pts, miss)
  RADNET_TILE   = unpack( canopy%rnet , l_tile_pts, miss)
  TOT_ALB       =UNPACK(rad%albedo_T,L_TILE_PTS, miss) 
  ECAN_TILE     = UNPACK(fev_dlh,  L_TILE_PTS, miss)
  !ESOIL_TILE    = UNPACK(fes_dlh, L_TILE_PTS, miss)
  SURF_HTF_TILE = UNPACK(canopy%ga,L_TILE_PTS,miss)
  
  !EI_TILE = 0.

  !Jun 2018 - need to split %fes into evaporation and sublimation
  fes_dlh = 0.
  WHERE (ssnow%cls==1.)  fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
  ESOIL_TILE = UNPACK(fes_dlh, L_TILE_PTS, miss)
  fes_dlh = 0.
  WHERE (ssnow%cls==1.1335) fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
  EI_TILE = UNPACK(fes_dlh, L_TILE_PTS, miss)

  !Jun 2018 - dtstar_tile unpacking
  DTSTAR_TILE = UNPACK(DTRAD, L_TILE_PTS, miss)

  SNOW_AGE = UNPACK(ssnow%snage, L_TILE_PTS, miss) 

  TRANSP_TILE = UNPACK(canopy%fevc, L_TILE_PTS, miss) 

  !unpack screen level (1.5m) variables - Convert back to K 
  t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, L_TILE_PTS, miss)
  q1p5m_tile     = UNPACK(canopy%qscrn, L_TILE_PTS, miss)

  CANOPY_TILE    = UNPACK(canopy%cansto, L_TILE_PTS, miss)
  CANOPY_GB      = SUM(um1%TILE_FRAC * CANOPY_TILE,2)

  !initialse full land grids and retain coastal grid fluxes
  DO N=1,NTILES

    DO K=1,um1%TILE_PTS(N)
  
      L = um1%TILE_INDEX(K,N)
      J=(um1%LAND_INDEX(L)-1)/ROW_LENGTH + 1
      I = um1%LAND_INDEX(L) - (J-1)*ROW_LENGTH
      
      IF( FLAND(L) == 1.0) THEN 
        FTL_1(I,J) =  0.0
        FQW_1(I,J) =  0.0
      ELSE
        !retain sea/ice contribution and remove land contribution
        FTL_1(I,J) = FTL_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                     FTL_TILE_old(L,N)
        FQW_1(I,J) = FQW_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                     FQW_TILE_old(L,N)
      ENDIF
      
      SURF_HT_FLUX_LAND(I,J) = 0.

    ENDDO !tile_pts(n)

  ENDDO !ntiles

  DO N=1,NTILES

    DO K=1,um1%TILE_PTS(N)

      L = um1%TILE_INDEX(K,N)
      J=(um1%LAND_INDEX(L)-1)/ROW_LENGTH + 1
      I = um1%LAND_INDEX(L) - (J-1)*ROW_LENGTH
      FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE(L,N)
      FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE(L,N)
      
      !Jun 2018 SURF_HT_FLUX_LAND is averaged over land not grid cell
      !SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
      !                         FLAND(L)*um1%TILE_FRAC(L,N) *              &
      !                         SURF_HTF_TILE(L,N)
      SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
                               um1%TILE_FRAC(L,N) * SURF_HTF_TILE(L,N)

    ENDDO

  ENDDO

  ! Initialise grid-cell carbon fields that are accumulated over tiles
  RESP_P = 0.;  NPP = 0.;  GPP = 0.;  RESP_S = 0.

  RESP_S_TILE   = UNPACK(canopy%frs, L_TILE_PTS, miss)
  NEE_TILE      = UNPACK(canopy%fnee, L_TILE_PTS, miss)
  NPP_FT        = UNPACK(canopy%fnpp, L_TILE_PTS, miss)
  G_LEAF        = UNPACK(canopy%frday,L_TILE_PTS, miss)
  RESP_P_FT     = UNPACK(canopy%frp, L_TILE_PTS, miss)

  GPP_FT = UNPACK(canopy%fnpp+canopy%frp+canopy%frday,  &
                  L_TILE_PTS, miss)

  ! convert from CABLE units (gC/m2/s) to UM units (kgC/m2/s)
  RESP_S_TILE = RESP_S_TILE*1.e-3
  G_LEAF      = G_LEAF*1.e-3
  NPP_FT      = NPP_FT*1.e-3
  GPP_FT      = GPP_FT*1.e-3
  RESP_P_FT   = RESP_P_FT*1.e-3

  ! If CASA-CNP used, plant and soil resp need to be passed into 
  ! variables that are dumped to restart, because CASA-CNP only run daily
  NPP_FT_ACC = RESP_S_TILE
  RESP_W_FT_ACC = RESP_P_FT

  DO N=1,NTILES 
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      NPP(L)=NPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT(L,N)
      GPP(L)=GPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT(L,N)
      RESP_P(L)=RESP_P(L)+FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT(L,N)

      !loop for soil resp. although DIM_CS1=1 (not 1 for triffid)
      DO I=1,DIM_CS1
        RESP_S(L,I) = RESP_S(L,I) + &
                      FLAND(L)*um1%TILE_FRAC(L,N)*RESP_S_TILE(L,N)
      ENDDO
      RESP_S_TOT(L)=sum(RESP_S(L,:))
      t1p5m(L)=sum(t1p5m_tile(L,:))
    ENDDO
  ENDDO

  !-------- End Unique subroutine body -----------

  fprintf_dir=trim(fprintf_dir_root)//trim(unique_subdir)//"/"
  if(L_cable_fprint) then 
    !basics to std output stream
    if (knode_gl == 0 .and. ktau_gl == 1)  call cable_fprintf(subr_name, .true.) 
    !more detailed output
    vname=trim(subr_name//'_')
    call cable_fprintf( cDiag00, vname, knode_gl, ktau_gl, .true. )
  endif

  if(L_cable_Pyfprint) then 
    !vname='canopy_tscrn'; dimx=mp
    !call cable_Pyfprintf( cDiag2, vname, (canopy%tscrn+tfrz), dimx, .true.)
    !vname='tscrn'; dimx=land_pts
    !call cable_Pyfprintf( cDiag2, vname, t1p5m, dimx, .true.)
  endif

return

END SUBROUTINE Implicit_unpack

End module cable_implicit_unpack_mod
