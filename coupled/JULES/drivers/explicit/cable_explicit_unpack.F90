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

IMPLICIT NONE

CONTAINS

SUBROUTINE cable_expl_unpack( latitude, longitude, ftl_tile, fqw_tile,        &
                              tstar_tile, u_s, u_s_std_tile, cd_tile, ch_tile,&
                              fland, radnet_tile, fraca, rESFS, resft,        &
                              z0h_tile, z0m_tile, recip_l_mo_tile, epot_tile, &
                              l_tile_pts, ssnow_snowd, ssnow_cls, air_rlam,   &
                              air_rho, canopy_fe, canopy_fh, canopy_us,       &
                              canopy_cdtq, canopy_fwet, canopy_wetfac_cs,     &
                              canopy_rnet, canopy_zetar, canopy_epot, met_ua, &
                              rad_trad, rad_transd, rough_z0m, rough_zref_tq, &
                              canopy_fes, canopy_fev )

 !subrs called 

  !processor number, timestep number / width, endstep
USE cable_common_module, ONLY: knode_gl, ktau_gl, kwidth_gl, kend_gl
USE cable_common_module, ONLY: cable_runtime
USE cable_phys_constants_mod, ONLY : ccapp   => capp
                                
!USE cable_common_module!, ONLY : cable_runtime, cable_user, &
                        !         ktau_gl, knode_gl, kend_gl, &
   
USE cable_def_types_mod, ONLY: mp, niter 
USE cable_um_tech_mod,   ONLY: um1
  
IMPLICIT NONE         

!___ re-decl input args
REAL,  DIMENSION(um1%row_length,um1%rows) :: latitude, longitude

!___ return fluxes UM vars recieve unpacked CABLE vars
REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) ::                      &
  tstar_tile,  & ! surface temperature
  ftl_tile,    & ! Surface FTL for land tiles     
  fqw_tile,    & ! Surface FQW for land tiles     
  z0h_tile,    & ! roughness
  z0m_tile,    & ! roughness
  cd_tile,     & ! Drag coefficient
  ch_tile,     & ! Transfer coefficient for heat & moisture
  u_s_std_tile,& ! Surface friction velocity
  radnet_tile,   & ! Surface net radiation
  resfs,         & ! Combined soil, stomatal & aerodynamic resistance
                   ! factor for fraction (1-FRACA) of snow-free land tiles
  resft,         & ! Total resistance factor.
                   ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                   ! 1 for snow.    
  fraca,         & ! Fraction of surface moisture
  recip_l_mo_tile,&! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
  epot_tile

REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  ::                     &
  u_s               ! Surface friction velocity (m/s)

LOGICAL,DIMENSION(um1%land_pts,um1%ntiles) :: l_tile_pts

!___UM vars used but NOT returned 
REAL, INTENT(IN), DIMENSION(um1%land_pts) ::                                  &
  fland(um1%land_pts)              ! IN Land fraction on land tiles.

!___ CABLE variables to be unpacked

REAL, INTENT(IN), DIMENSION(mp) ::                                            &
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
REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                   &
   le_tile

!___vars in local calc of Surface friction velocities
REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                   &
   u_s_tile
REAL, DIMENSION(mp)  ::                                                       &
   cdcab
  
REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                   &
   Lpts_nTILE

!___local miscelaneous
REAL, DIMENSION(mp)  ::                                                       &
  thetast,fraca_cab,rfsfs_cab, reciplmotile, fe_dlh
INTEGER :: i,j,k,n,l
REAL :: miss = 0.0
LOGICAL, SAVE :: first_cable_call = .TRUE.
LOGICAL :: Lunpack = .FALSE.
LOGICAL :: checks = .FALSE.
REAL :: tols = 0.25
  
! std template args 
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_unpack"

!___return fluxes
ftl_tile = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
ftl_tile = ftl_tile / ccapp

fe_dlh = ( canopy_fes / (air_rlam * ssnow_cls) )                              &
       + ( canopy_fev / air_rlam )

fqw_tile = UNPACK(fe_dlh, um1%l_tile_pts, miss)

!___return temp and roughness
tstar_tile = UNPACK(rad_trad,  um1%l_tile_pts, miss )

z0m_tile = UNPACK(rough_z0m,  um1%l_tile_pts, miss)
z0h_tile = z0m_tile
      
!___return friction velocities/drags/ etc
u_s_tile  =  UNPACK(canopy_us, um1%l_tile_pts, miss)
cdcab = canopy_us**2 / met_ua**2   ! met%ua is always above umin = 0.1m/s
cd_tile =  UNPACK(cdcab,um1%l_tile_pts, miss)
ch_tile =  UNPACK(canopy_cdtq,um1%l_tile_pts, miss)

u_s_std_tile = u_s_tile

u_s = 0.0
DO n = 1,um1%ntiles
  DO k = 1,um1%tile_pts(n)
    l = um1%tile_index(k,n)
    j=(um1%land_index(l) - 1) / um1%row_length + 1
    i = um1%land_index(l) - (j-1) * um1%row_length
    u_s(i,j) = u_s(i,j) + fland(l) * um1%tile_frac(l,n) * u_s_tile(l,n)
  END DO
END DO

!___return miscelaneous 
fraca_cab = canopy_fwet * (1.0 - rad_transd)
WHERE ( ssnow_snowd > 1.0 ) fraca_cab = 1.0
  
rfsfs_cab = MIN( 1.0, MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /             &
            MAX( 0.01,1.0 - fraca_cab ) )
fraca = UNPACK( fraca_cab, um1%l_tile_pts, miss )
resft = UNPACK( canopy_wetfac_cs,um1%l_tile_pts, miss )
resfs = UNPACK( rfsfs_cab , um1%l_tile_pts, miss )

radnet_tile = UNPACK( canopy_rnet , um1%l_tile_pts, miss )

thetast = ABS( canopy_fh ) / ( air_rho * ccapp * canopy_us )
reciplmotile =  canopy_zetar(:,niter) / rough_zref_tq
recip_l_mo_tile = UNPACK( reciplmotile, um1%l_tile_pts, miss )
epot_tile = UNPACK( canopy_epot, um1%l_tile_pts, miss )

IF (first_cable_call) THEN 
  l_tile_pts = um1%l_tile_pts
  first_cable_call = .FALSE.
END IF

RETURN

END SUBROUTINE cable_expl_unpack
    
END MODULE cable_expl_unpack_mod
