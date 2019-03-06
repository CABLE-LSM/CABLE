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
! Purpose:
!
! Called from: JULES: surf_couple_ pathway
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_implicit_main_mod
  
contains

subroutine cable_implicit_main( cycleno, & ! num_cycles
              row_length, rows, land_pts, ntiles, npft, sm_levels,             &
              dim_cs1, dim_cs2, Fland,                                         &
              !forcing: precip
              ls_rain_cable, conv_rain_cable, ls_snow_cable, conv_snow_cable,  &
              !forcing & increments: temp and humidity 
              tl_1, qw_1, dtl1_1, dqw1_1, ctctq1,                              &
              !prog: canopy water storage [per landpoint, per tile ]
              canopy_gb, canopy,                                               & 
              !prog: UM soi , per land_pt
              T_SOIL, smcl, STHF, STHU, snow_surft,                            & 
              ftl_1, ftl_surft, fqw_1, fqw_surft,                              &
              tstar_surft, surf_ht_flux_land, surf_htf_surft,                  &
              ecan_surft, esoil_surft, ei_surft, radnet_surft,                 &
              gs, gs_surft, t1p5m_surft, q1p5m_surft, melt_surft,              &
              NPP, NPP_FT, GPP, GPP_FT,                                        &
              RESP_S, RESP_S_TOT, & !RESP_S_TILE, !Kathy intro-ed as diag 
              RESP_P, RESP_P_FT, G_LEAF, snow_depth, dtstar_surft )
  !subrs called 
  USE cable_implicit_driv_mod, ONLY : cable_implicit_driver
  USE cable_implicit_unpack_mod, ONLY : implicit_unpack

  ! processor number, timestep number & width
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, cable_runtime
  
# if defined(UM_JULES)
  ! CABLE prognostics declared at top_level
  USE atm_fields_real_mod, ONLY : soil_temp_cable, soil_moist_cable,           &
                                  soil_froz_frac_cable, snow_dpth_cable,       & 
                                  snow_mass_cable, snow_temp_cable,            &
                                  snow_rho_cable, snow_avg_rho_cable,          &   
                                  snow_age_cable, snow_flg_cable,              &
                                  aquifer_moist_cable
  ! CASA prognostics declared at top_level
  USE atm_fields_real_mod, ONLY : C_pool_casa, N_pool_casa, P_pool_casa,       &
                                  SOIL_ORDER_casa, N_DEP_casa, N_FIX_casa,     &
                                  P_DUST_casa, P_weath_casa, LAI_casa,         &
                                  PHENPHASE_casa, NPP_PFT_ACC, RSP_W_PFT_ACC
  !UM: time info 
  USE model_time_mod, ONLY:    target_end_stepim, i_day, i_day_number
#endif

  USE cable_gather_um_data_decs, ONLY : smvcst_cable  
  USE atmos_physics2_alloc_mod, ONLY : resp_s_tile
  
     
  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir
  
  implicit none
  
  !___ re-decl input args
  integer :: cycleno
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
  integer :: dim_cs1, dim_cs2 
  
  REAL,  DIMENSION(land_pts) :: & 
      fland
   
  REAL, DIMENSION(row_length,rows) :: &
    ls_rain_cable,    &!forcing: precip:rain , large scale
    ls_snow_cable,    &!forcing: precip: snow, large scale
    conv_rain_cable,  &!forcing: precip:rain , convective 
    conv_snow_cable    !forcing: precip:snow, convective 

  REAL, DIMENSION(row_length,rows) :: &
   tl_1, qw_1,       &  !forcing: temp and humidity 
   dtl1_1, dqw1_1       !forcing: increments to temp and humidity 

  !Ticket #132 needs ctctq1
  REAL, DIMENSION(row_length,rows) :: &
   ctctq1               !forcing: information for temp and humidity increment
  
  REAL ::                                                                       &
    canopy_gb(land_pts),         & !prog: canopy water store aggregate over tiles 
    canopy(land_pts, ntiles)       !prog:  per tile

  !prog: UM soil quantities, non-tiled  aggregate over CABLE tiles per land_pt
  real, dimension(land_pts,sm_levels) ::                                       &
    T_SOIL, &
    smcl, &
    STHF,&
    STHU
  
  real, dimension( land_pts, ntiles ) ::                                       &
    snow_surft, &     ! snow ammount on tile nee:snow_tile 
    ftl_surft, &
    fqw_surft,  &
    tstar_surft, &
    surf_htf_surft, &
    ecan_surft, esoil_surft,                 &
    ei_surft, radnet_surft, &
    gs_surft, &
    t1p5m_surft, &
    q1p5m_surft, &
    melt_surft, &
    snow_depth, &
    dtstar_surft
    
  real, dimension( land_pts ) ::                           &
    gs, &             ! conductance for use in dust scheme
    dtrad
    
  real, dimension(row_length,rows) ::      &
    ftl_1, &
    fqw_1, &
    surf_ht_flux_land

  real, dimension( land_pts ) ::                                               &
    gpp,   & ! IN Gross primary productivity (kg C/m2/s).
    npp,   & ! IN Net primary productivity   (kg C/m2/s).
    resp_p   ! IN Plant respiration (kg C/m2/s).

  real, dimension( land_pts,ntiles ) ::                           &
    gpp_ft,   & ! IN Gross primary productivity  on PFTs (kg C/m2/s).
    npp_ft,   & ! IN Net primary productivity  on PFTs (kg C/m2/s).
    g_leaf,   & ! IN Leaf turnover rate (/360days).
    resp_p_ft   ! IN Plant respiration on PFTs  (kg C/m2/s).
   
   real :: resp_s(land_pts,dim_cs1) ! IN Soil respiration (kg C/m2/s).
   real :: resp_s_tot(dim_cs2)      ! IN Total soil resp'n (kg C/m2/s).
 
  !___ local vars

  ! UM type vars but no feedback 
  REAL, dimension(land_pts,ntiles,sm_levels) :: STHU_surft
  REAL, dimension(land_pts,ntiles) ::  TOT_ALB, TRANSP_surft
  REAL, dimension(land_pts,ntiles,3) :: SNOW_COND        

  logical, save :: first_call = .true.
  integer,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable

  ! std template args  
  character(len=*), parameter :: subr_name = "cable_implicit_main"

# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
  
  ! FLAGS def. specific call to CABLE from UM
  cable_runtime%um          = .TRUE.
  cable_runtime%um_implicit = .TRUE.
  cable_runtime%um_explicit = .FALSE.
  
  isnow_flg_cable = int(snow_flg_cable)

  call cable_implicit_driver( i_day_number, cycleno, &! num_cycles  
        row_length,rows, land_pts, ntiles, npft, sm_levels,                    &
        dim_cs1, dim_cs2, Fland,                                               &
        LS_RAIN_cable, CONV_RAIN_cable, LS_SNOW_cable, CONV_SNOW_cable,        &
        DTL1_1, DQW1_1, ctctq1, T_SOIL, soil_temp_cable,                       &
        SMCL, soil_moist_cable, aquifer_moist_cable, real(kwidth_gl),          &
        SMVCST_cable, STHF, soil_froz_frac_cable, STHU,                        &
        snow_surft, snow_avg_rho_cable, isnow_flg_cable, snow_dpth_cable,      &
        snow_mass_cable, snow_rho_cable, snow_temp_cable,                      & 
        FTL_1, FTL_surft, FQW_1, FQW_surft,                                    &
        TSTAR_surft, SURF_HT_FLUX_LAND, ECAN_surft, ESOIL_surft,               &
        EI_surft, RADNET_surft, SNOW_AGE_cable, CANOPY, GS, gs_surft,          &
        T1P5M_surft, Q1P5M_surft, CANOPY_GB, MELT_surft,                       &
        NPP, NPP_FT, GPP, GPP_FT, RESP_S, RESP_S_TOT, RESP_S_TILE,             &
        RESP_P, RESP_P_FT, G_LEAF, TL_1, QW_1, SURF_HTF_surft,                 &
        C_pool_casa, N_pool_casa, P_pool_casa, LAI_casa, PHENPHASE_casa,       &
        NPP_PFT_ACC, RSP_W_PFT_ACC, dtrad)
 
   CALL implicit_unpack( cycleno, row_length,rows, land_pts, ntiles, npft,     &
                        sm_levels, dim_cs1, dim_cs2, T_SOIL, soil_temp_cable,  &
                        SMCL, soil_moist_cable, aquifer_moist_cable,           &
                        SMVCST_cable, STHF, &
                        soil_froz_frac_cable, STHU, STHU_surft, snow_surft,    &
                        snow_avg_rho_cable, isnow_flg_cable, snow_dpth_cable,  &
                        snow_mass_cable, snow_rho_cable, snow_temp_cable,      &
                        SNOW_COND, FTL_1, FTL_surft, FQW_1, FQW_surft,         &
                        TSTAR_surft, SURF_HT_FLUX_LAND, ECAN_surft,            &
                        ESOIL_surft, EI_surft, RADNET_surft, TOT_ALB,          &
                        SNOW_AGE_cable, CANOPY, GS,GS_surft,                   &
                        T1P5M_surft, Q1P5M_surft,                              &
                        CANOPY_GB, FLAND, MELT_surft,                          &
                        NPP, NPP_FT, GPP, GPP_FT, RESP_S, RESP_S_TOT,          &
                        RESP_S_tile, RESP_P, RESP_P_FT, G_LEAF, TRANSP_surft,  &
                        NPP_PFT_ACC, RSP_W_PFT_ACC,SURF_HTF_surft,             &
                        dtrad, dtstar_surft)

  snow_flg_cable = real(isnow_flg_cable)
    
  first_call = .false.        

  cable_runtime%um_implicit = .FALSE.

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
    !vname='latitude'; dimx=mp
    !call cable_Pyfprintf( cDiag1, vname, cable%lat, dimx, .true.)
  endif

return

End subroutine cable_implicit_main
  
End module cable_implicit_main_mod


