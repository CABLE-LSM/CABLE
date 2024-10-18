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

SUBROUTINE cable_expl_unpack( & 

  ! IN: UM/JULES/CABLE model/grid parameters, fields, mappings
  row_length, rows, land_pts, nsurft, npft, mp, land_index, surft_pts,         &
  surft_index, l_tile_pts, fland, tile_frac, latitude, longitude,              &

  !OUT: UM fields to be updated
  ftl_tile, fqw_tile, tstar_tile, dtstar_tile , u_s, u_s_std_tile, cd_tile,    &
  ch_tile, radnet_tile, fraca, resfs, resft, z0h_tile, z0m_tile,               &
  recip_l_mo_tile, epot_tile,                                                  &
 
  !IN: UM fields to be updated FROM these CABLE fields
  canopy_fh, canopy_fes, canopy_fev, canopy_us, canopy_cdtq, canopy_fwet,      &
  canopy_wetfac_cs, canopy_rnet, canopy_zetar, canopy_epot,                    &
  rad_trad, rad_otrad, rad_transd, rough_z0m, rough_zref_tq,                   &

  !IN: CABLE fields used in derivation of fields to be updated
  ssnow_snowd, ssnow_cls, air_rlam, air_rho, met_ua )

USE cable_def_types_mod, ONLY : NITER 
USE cable_phys_constants_mod,  ONLY: CAPP 
  
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: npft              ! # plant functional types
INTEGER, INTENT(IN) :: mp 
INTEGER, INTENT(IN) :: land_index(land_pts)          ! land point indices 
                                                     ! recipe back to (i,j) cell  
INTEGER, INTENT(IN) :: surft_pts(nsurft)             ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! tile points indices
                                                     ! recipe back to land_index
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, nsurft)
REAL, INTENT(IN)    :: fland(land_pts)               ! Land fraction
REAL, INTENT(IN)    :: tile_frac(land_pts, nsurft)
REAL, INTENT(IN)    :: latitude(row_length,rows)
REAL, INTENT(IN)    :: longitude(row_length,rows) 

!___ return fluxes UM vars recieve unpacked CABLE vars
REAL, INTENT(OUT) :: tstar_tile(land_pts,nsurft) ! surface temperature
REAL, INTENT(OUT) :: dtstar_tile(land_pts,nsurft) ! surface temperature
REAL, INTENT(OUT) :: ftl_tile(land_pts,nsurft)   ! Surface FTL for land tiles     
REAL, INTENT(OUT) :: fqw_tile(land_pts,nsurft)   ! Surface FQW for land tiles     
REAL, INTENT(OUT) :: z0h_tile(land_pts,nsurft)   ! roughness
REAL, INTENT(OUT) :: z0m_tile(land_pts,nsurft)   ! roughness
REAL, INTENT(OUT) :: cd_tile(land_pts,nsurft)    ! Drag coefficient
REAL, INTENT(OUT) :: ch_tile(land_pts,nsurft)    ! Transfer coefficient for heat & moisture
REAL, INTENT(OUT) :: u_s_std_tile(land_pts,nsurft)  ! Surface friction velocity
REAL, INTENT(OUT) :: radnet_tile(land_pts,nsurft)   ! Surface net radiation
REAL, INTENT(OUT) :: resfs(land_pts,nsurft) ! Combined soil, stomatal & aero resistance
                   ! factor for fraction (1-FRACA) of snow-free land tiles
REAL, INTENT(OUT) :: RESFT(land_pts,nsurft) ! Total resistance factor.
                   ! FRACA+(1-FRACA)*RESFS where NO snow
                   ! 1 for snow.    
REAL, INTENT(OUT) :: FRACA(land_pts,nsurft)            ! Fraction of surface moisture
REAL, INTENT(OUT) :: RECIP_L_MO_TILE(land_pts,nsurft)  ! Reciprocal Monin-Obukhov len (m^-1)
REAL, INTENT(OUT) :: EPOT_TILE(land_pts,nsurft)
REAL, INTENT(OUT) :: U_S(row_length,rows) ! Surface friction velocity (m/s)

!___ CABLE variables to be unpacked
REAL, INTENT(IN) :: rad_trad(mp)           ! 
REAL, INTENT(IN) :: rad_otrad(mp)          ! 
REAL, INTENT(IN) :: rad_transd(mp)         ! rad. temp. (soil and veg)
REAL, INTENT(IN) :: canopy_fh(mp)          ! total sensible heat (W/m2) 
REAL, INTENT(IN) :: canopy_fes(mp)         !  
REAL, INTENT(IN) :: canopy_fev(mp)         ! 
REAL, INTENT(IN) :: canopy_fwet(mp)        ! fraction of canopy wet 
REAL, INTENT(IN) :: canopy_wetfac_cs(mp)   ! fraction of canopy wet
REAL, INTENT(IN) :: canopy_us(mp)          ! friction velocity 
REAL, INTENT(IN) :: canopy_cdtq(mp)        ! drag coefficient for momentum
REAL, INTENT(IN) :: canopy_rnet(mp)        ! net rad. absorbed by surface (W/m2) 
REAL, INTENT(IN) :: canopy_epot(mp)        ! total potential evaporation        
REAL, INTENT(IN) :: rough_z0m(mp)          ! roughness length 
REAL, INTENT(IN) :: rough_zref_tq(mp)      ! Reference height for met forcing
REAL, INTENT(IN) :: canopy_zetar(mp,niter) ! stability correction

!IN: CABLE fields used in derivation of fields to be updated
REAL, INTENT(IN) :: ssnow_snowd(mp)        ! snow depth (liquid water)
REAL, INTENT(IN) :: ssnow_cls(mp)          ! factor for latent heat
REAL, INTENT(IN) :: met_ua(mp)             ! surface wind speed (m/s)
REAL, INTENT(IN) :: air_rlam(mp)           ! latent heat for water (j/kg)
REAL, INTENT(IN) :: air_rho(mp)            ! dry air density (kg m-3) 

!___ local vars
REAL :: u_s_tile(land_pts,nsurft)
REAL :: dtrad(mp)                          ! change in rad%trad over time step
REAL :: cdcab(mp)
REAL :: thetast(mp)
REAL :: fraca_cab(mp)
REAL :: rfsfs_cab(mp)
REAL :: reciplmotile(mp)
REAL :: fe_dlh(mp)
REAL :: miss = 0.0
REAL :: miss_tiny = 1.0e-9 
INTEGER :: i,j,k,n,L
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_explicit_unpack"

!-------- Unique subroutine body -----------
!___return fluxes
ftl_tile = UNPACK(canopy_fh,  l_tile_pts, miss)
ftl_tile = ftl_tile / CAPP
fe_dlh   = ( canopy_fes/(air_rlam*ssnow_cls) ) + ( canopy_fev/air_rlam )
fqw_tile = UNPACK(fe_dlh, l_tile_pts, miss)

!___return surface temp and roughness
tstar_tile  = UNPACK(rad_trad,  l_tile_pts, miss )
dtrad       = rad_trad - rad_otrad        
dtstar_tile = UNPACK(dtrad,  l_tile_pts, miss )
z0m_tile    = UNPACK(rough_z0m,  l_tile_pts, miss_tiny)
z0h_tile    = z0m_tile
    
!___return friction velocities/drags/ etc
u_s_tile  = UNPACK( canopy_us, l_tile_pts, miss)
cdcab     = canopy_us**2 / met_ua**2   ! met%ua is always above umin = 0.1m/s
cd_tile   = UNPACK( cdcab,l_tile_pts, miss)
ch_tile   = UNPACK( canopy_cdtq,l_tile_pts, miss)

u_s_std_tile=u_s_tile

u_s = 0.0
DO n=1, nsurft
  DO k=1,surft_pts(n)
    l = surft_index(K,N)
    j = (land_index(l)-1) / row_length + 1
    i = land_index(l) - (j-1)*row_length
    u_s(i,j) = u_s(i,j)+fland(l)*tile_frac(l,n)*u_s_tile(l,n)
  ENDDO
ENDDO

!___return miscelaneous 
fraca_cab = canopy_fwet * (1.-rad_transd)
WHERE( ssnow_snowd > 1.0 ) fraca_cab = 1.0

rfsfs_cab       = MIN( 1., MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /         &
                           MAX( 0.01, 1.0 - fraca_cab ) )
fraca           = UNPACK( fraca_cab, l_tile_pts, miss )
resft           = UNPACK( canopy_wetfac_cs,l_tile_pts, miss )
resfs           = UNPACK( rfsfs_cab , l_tile_pts, miss )
radnet_tile     = UNPACK( canopy_rnet , l_tile_pts, miss )
thetast         = ABS( canopy_fh ) / ( air_rho * capp*canopy_us )
reciplmotile    = canopy_zetar(:,niter) / rough_zref_tq
recip_l_mo_tile = UNPACK( reciplmotile, l_tile_pts, miss )
epot_tile       = UNPACK( canopy_epot, l_tile_pts, miss )

RETURN

END SUBROUTINE cable_expl_unpack
    
END MODULE cable_expl_unpack_mod
