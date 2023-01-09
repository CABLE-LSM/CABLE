!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
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
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: cable_explicit_driver
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

MODULE cable_expl_unpack_mod

!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!

implicit none

contains

SUBROUTINE cable_expl_unpack( latitude, longitude, FTL_TILE, FQW_TILE,         &
                              TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE, &
                              FLAND, RADNET_TILE, FRACA, rESFS, RESFT,         &
                              Z0H_TILE, Z0M_TILE, RECIP_L_MO_TILE, EPOT_TILE,  &
                              l_tile_pts, ssnow_snowd, ssnow_cls, air_rlam,    &
                              air_rho, canopy_fe, canopy_fh, canopy_us,        &
                              canopy_cdtq, canopy_fwet, canopy_wetfac_cs,      &
                              canopy_rnet, canopy_zetar, canopy_epot, met_ua,  &
                              rad_trad, rad_transd, rough_z0m, rough_zref_tq,  &
                              canopy_fes, canopy_fev )

  !subrs called 
  
  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
                                
  USE cable_common_module!, ONLY : cable_runtime, cable_user, &
                          !         ktau_gl, knode_gl, kend_gl, &
   
  USE cable_def_types_mod, ONLY : mp, NITER 
USE cable_phys_constants_mod,  ONLY: CAPP 
  USE cable_um_tech_mod,   ONLY : um1, basic_diag
  
  implicit none         

  !___ re-decl input args
  REAL,  DIMENSION(um1%row_length,um1%rows) :: latitude, longitude

  !___ return fluxes UM vars recieve unpacked CABLE vars
  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
    TSTAR_TILE,  & ! surface temperature
    FTL_TILE,    & ! Surface FTL for land tiles     
    FQW_TILE,    & ! Surface FQW for land tiles     
    Z0H_TILE,    & ! roughness
    Z0M_TILE,    & ! roughness
    CD_TILE,     & ! Drag coefficient
    CH_TILE,     & ! Transfer coefficient for heat & moisture
    U_S_STD_TILE,& ! Surface friction velocity
    RADNET_TILE,   & ! Surface net radiation
    RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                     ! factor for fraction (1-FRACA) of snow-free land tiles
    RESFT,         & ! Total resistance factor.
                     ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                     ! 1 for snow.    
    FRACA,         & ! Fraction of surface moisture
    RECIP_L_MO_TILE,&! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
    EPOT_TILE

  REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  :: &
    U_S               ! Surface friction velocity (m/s)

  LOGICAL,DIMENSION(um1%land_pts,um1%ntiles) :: l_tile_pts

  !___UM vars used but NOT returned 
  REAL, INTENT(IN), DIMENSION(um1%land_pts) ::   &
    FLAND(um1%land_pts)              ! IN Land fraction on land tiles.

  !___ CABLE variables to be unpacked

  REAL, INTENT(IN), DIMENSION(mp) :: &
    ssnow_snowd,      & ! snow depth (liquid water)
    ssnow_cls,        & ! factor for latent heat
    met_ua,           & ! surface wind speed (m/s)
    air_rlam,         & ! latent heat for water (j/kg)
    air_rho,          & ! dry air density (kg m-3) 
    rad_trad,         & ! frac SW diffuse transmitted thru canopy
    rad_transd,       & !  rad. temp. (soil and veg)
    canopy_fe,        & ! total latent heat (W/m2)
    canopy_fh,        & ! total sensible heat (W/m2) 
    canopy_fes,       & !  
    canopy_fev,       & ! 
    canopy_fwet,      & ! fraction of canopy wet 
    canopy_wetfac_cs, & ! fraction of canopy wet
    canopy_us,        & ! friction velocity 
    canopy_cdtq,      & ! drag coefficient for momentum
    canopy_rnet,      & ! net rad. absorbed by surface (W/m2) 
    canopy_epot,      & ! total potential evaporation        
    rough_z0m,        & ! roughness length 
    rough_zref_tq       ! Reference height for met forcing
  
  REAL, INTENT(IN), DIMENSION(mp,niter) :: canopy_zetar! stability correction
  
  !___ local vars

  !___vars in local calc. of latent heat fluxes
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
     LE_TILE

  !___vars in local calc of Surface friction velocities
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
     U_S_TILE
  REAL, DIMENSION(mp)  :: &
     CDCAB
  
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
     Lpts_nTILE

  !___local miscelaneous
  REAL, DIMENSION(mp)  :: &
    THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
  INTEGER :: i,j,k,N,L
  REAL :: miss = 0.0
  LOGICAL, SAVE :: first_cable_call = .true.
  LOGICAL :: Lunpack = .false.
  LOGICAL :: checks = .false.
  real :: tols = 0.25
  
  ! std template args 
  character(len=*), parameter :: subr_name = "cable_explicit_unpack"

  !-------- Unique subroutine body -----------
  !___return fluxes
  FTL_TILE = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
  FTL_TILE = FTL_TILE / CAPP

  fe_dlh = ( canopy_fes/(air_rlam*ssnow_cls) )  &
         + ( canopy_fev/air_rlam )

  FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

  !___return temp and roughness
  TSTAR_TILE = UNPACK(rad_trad,  um1%l_tile_pts, miss )

  Z0M_TILE = UNPACK(rough_z0m,  um1%l_tile_pts, miss)
  Z0H_TILE = Z0M_TILE
      
  !___return friction velocities/drags/ etc
  U_S_TILE  =  UNPACK(canopy_us, um1%l_tile_pts, miss)
  CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
  CD_TILE =  UNPACK(CDCAB,um1%l_tile_pts, miss)
  CH_TILE =  UNPACK(canopy_cdtq,um1%l_tile_pts, miss)

  U_S_STD_TILE=U_S_TILE

  U_S = 0.
  DO N=1,um1%ntiles
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
      I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
      U_S(I,J) = U_S(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*U_S_TILE(L,N)
    ENDDO
  ENDDO

  !___return miscelaneous 
  fraca_cab = canopy_fwet * (1.-rad_transd)
  WHERE( ssnow_snowd > 1.0 ) fraca_cab = 1.0
  
  rfsfs_cab = MIN( 1., MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /         &
              MAX( 0.01,1. - fraca_cab ) )
  FRACA = UNPACK( fraca_cab, um1%l_tile_pts, miss )
  RESFT = UNPACK( canopy_wetfac_cs,um1%l_tile_pts, miss )
  RESFS = UNPACK( rfsfs_cab , um1%l_tile_pts, miss )

  RADNET_TILE = UNPACK( canopy_rnet , um1%l_tile_pts, miss )

  THETAST = ABS( canopy_fh ) / ( air_rho * capp*canopy_us )
  RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
  RECIP_L_MO_TILE = UNPACK( RECIPLMOTILE, um1%l_tile_pts, miss )
  EPOT_TILE = UNPACK( canopy_epot, um1%l_tile_pts, miss )

  IF(first_cable_call) THEN 
    l_tile_pts = um1%l_tile_pts
    first_cable_call = .FALSE.
  ENDIF
  !-------- End Unique subroutine body -----------

return

End subroutine cable_expl_unpack
    
End module cable_expl_unpack_mod
