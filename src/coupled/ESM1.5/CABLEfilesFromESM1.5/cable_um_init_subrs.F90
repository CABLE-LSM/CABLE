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
! Purpose: Routines to pass UM variables into appropriate CABLE variables and 
!          to map parameters for each surface type to CABLE arrays
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Rewrites of code in v1.8 (ACCESS1.3)
!
!
! ==============================================================================

MODULE cable_um_init_subrs_mod
   IMPLICIT NONE


CONTAINS

!jhan: code under development for future release
!   subroutine initialize_maps(latitude,longitude, tile_index_mp)
!   USE cable_phys_constants_mod, ONLY: 
!      use cable_um_tech_mod, only : um1
!      use define_dimensions, only : mp
!
!      use cable_diag_module, only : cable_diag 
!      use cable_common_module, only : ktau_gl, knode_gl, cable_user 
!         
!      implicit none
!      real, intent(in), dimension(um1%row_length,um1%rows) :: &
!         latitude, longitude
!      integer, intent(in), dimension(um1%land_pts, um1%ntiles) :: &
!         tile_index_mp  ! index of tile
!          
!      logical, save :: first_call = .true.
!      
!      INTEGER :: j
!
!      real :: dlon
!      real, dimension(um1%row_length) :: tlong,acoslong
!      real, dimension(um1%row_length, um1%rows) :: new_longitude
!
!
!           
!            allocate( cable%lat(mp), cable%lon(mp), cable%tile(mp), cable%tile_frac(mp) )
!
!            !-------------------------------------   
!            !---make indexes for tile, lat, lon
!            !-------------------------------------   
!           
!            !--- get latitude index corresponding to cable points
!            call um2cable_rr( (asin(latitude)/const%math%pi180), cable%lat )
!
!            !--- get longitude index corresponding to cable points.
!            !--- this is not so straight forward as UM longitude index 
!            !--- contains ambiguity. thus define "new_longitude" first
!            acoslong =  acos( longitude(:,1) ) /const%math%pi180  
!       
!            tlong(1) = acoslong(1)
!            do j=2, um1%row_length
!               if( acoslong(j) < acoslong(j-1) ) then  
!                  dlon = acoslong(j) - acoslong(j-1)
!                  tlong(j) = tlong(j-1) - dlon   
!               else 
!                  tlong(j) = acoslong(j)
!               endif           
!            enddo
!            
!            do j=1, um1%row_length
!               new_longitude(j,:) = tlong(j)
!            enddo
!            
!            call um2cable_rr( new_longitude, cable%lon )
!         
!         
!            !--- get tile index/fraction  corresponding to cable points
!            cable%tile = pack(tile_index_mp, um1%l_tile_pts)
!            cable%tile_frac = pack(um1%tile_frac, um1%l_tile_pts)
!
!         !--- write all these maps.  cable_user%initialize_mapping can be 
!         !--- set in namelist cable.nml
!         if ( cable_user%initialize_mapping ) then
!            !write indexes for tile, lat, lon
!            call cable_diag( 1, 'latitude', um1%rows, 1, ktau_gl,  & 
!                  knode_gl, 'latitude', ( asin( latitude(1,:) ) /const%math%pi180 ) ) 
!            
!            call cable_diag( 1, 'longitude', um1%row_length, 1, ktau_gl,  & 
!                  knode_gl, 'longitude', ( new_longitude(:,1) ) ) 
!        
!            !write indexes for tile, lat, lon
!            call cable_diag( 1, 'lat_index', mp, 1, ktau_gl,  & 
!                  knode_gl, 'lat', cable%lat )
!            call cable_diag( 1, 'lon_index', mp, 1, ktau_gl,  & 
!                  knode_gl, 'lon', cable%lon )
!            
!            !this should be integer-ed. typecast for now
!            call cable_diag( 1, 'tile_index', mp, 1, ktau_gl,  & 
!                  knode_gl, 'tile', real(cable%tile) )
!            
!            call cable_diag( 1, 'tile_frac', mp, 1, ktau_gl,  & 
!                  knode_gl, 'tile_frac', cable%tile_frac )
!            
!          endif  
!         
!      
!      return
!   end subroutine initialize_maps
        
SUBROUTINE initialize_veg( clobbered_htveg, land_pts, npft, ntiles, ms, mp,      &
                           canht_ft, lai_ft, soil_zse, veg,                  &
                    tile_pts, tile_index, tile_frac, L_tile_pts,                 &
                    CLAI_thresh )

USE cable_common_module, ONLY : knode_gl 
USE cable_def_types_mod,  ONLY: veg_parameter_type
USE cbl_LAI_canopy_height_mod,  ONLY: limit_HGT_LAI
   
INTEGER ::  land_pts
INTEGER ::  ntiles
INTEGER ::  npft
INTEGER ::  ms                ! soil levels
INTEGER ::  mp                ! active pts CABLE
REAL    ::  CLAI_thresh  
INTEGER ::  tile_pts(ntiles)  ! number of land_pts per tile type
REAL    ::  tile_frac(land_pts, ntiles)
INTEGER ::  tile_index(land_pts, ntiles)  ! index of tile 
LOGICAL ::  L_tile_pts(land_pts, ntiles)  ! true IF vegetation (tile) fraction is greater than 0
REAL    ::  clobbered_htveg(land_pts, ntiles)
REAL    :: soil_zse(ms)       !soil layer thickness [dzsoil]
TYPE(veg_parameter_type), INTENT(INOUT) :: veg
REAL, INTENT(IN) :: canht_ft(land_pts, npft)
REAL, INTENT(IN) :: lai_ft(land_pts, npft) 
LOGICAL, SAVE :: first_call= .TRUE. ! defs 1st call to CABLE in this run
INTEGER :: JSurfaceTypeID(land_pts,ntiles)
INTEGER :: i,j

!--- veg params were read from initialize_soil() 
IF(first_call)  THEN
  !local var to pack surface type:
  JSurfaceTypeID = 0
  DO j = 1, land_pts
    DO i = 1,ntiles
      IF ( tile_frac(j,i) > 0.0 ) JSurfaceTypeID(j,i) = i
    END DO
  END DO
  veg%iveg = PACK( JSurfaceTypeID, L_tile_pts)

ENDIF

! limit IN height, LAI  and initialize existing cable % types
CALL limit_HGT_LAI( clobbered_htveg, veg%vlai, veg%hc, mp, land_pts, ntiles,           &
                    tile_pts, tile_index, tile_frac, L_tile_pts,              &
                    LAI_ft, canht_ft, CLAI_thresh )
     
      
      !--- veg params were read from initialize_soil() 
      IF(first_call)  THEN
         CALL init_veg_pars_fr_vegin( soil_zse ) 
         ! Fix in-canopy turbulence scheme globally:
         veg%meth = 1
      ENDIF
      first_call= .FALSE.
     
END SUBROUTINE initialize_veg

!========================================================================
!========================================================================
!========================================================================
   
SUBROUTINE init_veg_pars_fr_vegin( soil_zse ) 

USE cable_um_tech_mod,   ONLY: veg, soil 
USE cable_def_types_mod, ONLY: mp, ms

IMPLICIT NONE

REAL :: soil_zse(ms) 

CALL init_veg_from_vegin(1, mp, veg, soil_zse) 

!rml 20/11/24 needed for fwsoil_switch=Haverd2013
veg%gamma = 3.e-2

RETURN

END SUBROUTINE init_veg_pars_fr_vegin

SUBROUTINE init_veg_from_vegin(ifmp,fmp, veg, soil_zse ) 

USE cable_def_types_mod,  ONLY: veg_parameter_type, ms, mp 
USE cable_pft_params_mod, ONLY: vegin
USE cable_common_module,  ONLY: cable_user

     integer ::  ifmp,  & ! start local mp, # landpoints (jhan:when is this not 1 )      
                 fmp     ! local mp, # landpoints       
     real, dimension(ms) :: soil_zse 
  
     type(veg_parameter_type) :: veg
     
    INTEGER :: is
    REAL :: totdepth
     integer :: k,h
     
!! Prescribe parameters for current gridcell based on veg/soil type (which
!! may have loaded from default value file or met file):
DO h = 1, mp          ! over each patch in current grid

  veg%canst1(h)   = vegin%canst1(veg%iveg(h))
  veg%ejmax(h)    = 2.* vegin%vcmax(veg%iveg(h)) ! instead of read vegin%ejmax!?
  veg%frac4(h)    = vegin%frac4(veg%iveg(h))
  veg%tminvj(h)   = vegin%tminvj(veg%iveg(h))
  veg%tmaxvj(h)   = vegin%tmaxvj(veg%iveg(h))
  veg%vbeta(h)    = vegin%vbeta(veg%iveg(h))
  veg%rp20(h)     = vegin%rp20(veg%iveg(h))
  veg%rpcoef(h)   = vegin%rpcoef(veg%iveg(h))
  veg%shelrb(h)   = vegin%shelrb(veg%iveg(h))
  veg%vegcf(h)    = vegin%vegcf(veg%iveg(h))
  veg%extkn(h)    = vegin%extkn(veg%iveg(h))
  veg%vcmax(h)    = vegin%vcmax(veg%iveg(h))
  veg%xfang(h)    = vegin%xfang(veg%iveg(h))
  veg%dleaf(h)    = vegin%dleaf(veg%iveg(h))
  veg%xalbnir(h)  = vegin%xalbnir(veg%iveg(h))
  veg%rs20(h)     = vegin%rs20(veg%iveg(h))
  veg%wai(h)      = vegin%wai(veg%iveg(h))
  veg%a1gs(h)     = vegin%a1gs(veg%iveg(h))
  veg%d0gs(h)     = vegin%d0gs(veg%iveg(h))
  veg%g0(h)       = vegin%g0(veg%iveg(h))         ! Ticket #56
  veg%g1(h)       = vegin%g1(veg%iveg(h))         ! Ticket #56
  veg%alpha(h)    = vegin%alpha(veg%iveg(h))
  veg%convex(h)   = vegin%convex(veg%iveg(h))
  veg%cfrd(h)     = vegin%cfrd(veg%iveg(h))
  veg%gswmin(h)   = vegin%gswmin(veg%iveg(h))
  veg%conkc0(h)   = vegin%conkc0(veg%iveg(h))
  veg%conko0(h)   = vegin%conko0(veg%iveg(h))
  veg%ekc(h)      = vegin%ekc(veg%iveg(h))
  veg%eko(h)      = vegin%eko(veg%iveg(h))
  veg%rootbeta(h) = vegin%rootbeta(veg%iveg(h))
  veg%zr(h)       = vegin%zr(veg%iveg(h))
  veg%clitt(h)    = vegin%clitt(veg%iveg(h))
  !veg%hc(h)       = vegin%hc(veg%iveg(h))
  veg%taul(h,1)   = vegin%taul1(veg%iveg(h))
  veg%taul(h,2)   = vegin%taul2(veg%iveg(h))
  veg%refl(h,1)   = vegin%refl1(veg%iveg(h))
  veg%refl(h,2)   = vegin%refl2(veg%iveg(h))
            
 
END DO ! over each veg patch in land point
 
IF (cable_user%access13roots) THEN

  DO h = 1, mp          ! over each patch in current grid
    ! prescribed values from vegin%froot so can be set to ESM1.5 values
    ! which were the same for all pfts (0.05, 0.2, 0.2, 0.2, 0.2, 0.15)
    ! or alternate prescribed root fractions could be used with pft dependence
    veg%froot(h,1) = vegin%froot1(veg%iveg(h))
    veg%froot(h,2) = vegin%froot2(veg%iveg(h))
    veg%froot(h,3) = vegin%froot3(veg%iveg(h))
    veg%froot(h,4) = vegin%froot4(veg%iveg(h))
    veg%froot(h,5) = vegin%froot5(veg%iveg(h))
    veg%froot(h,6) = vegin%froot6(veg%iveg(h))
  END DO ! over each veg patch in land point

ELSE

  ! parametrized values used in CMIP6 ACCESS-CM2
  ! calculate vegin%froot from using rootbeta and soil depth
  ! (Jackson et al. 1996, Oceologica, 108:389-411)
  totdepth = 0.0
  DO is = 1, ms-1
    totdepth = totdepth + soil_zse(is) * 100.0  ! unit in centimetres
    veg%froot(:, is) = MIN( 1.0, 1.0-veg%rootbeta(:)**totdepth )
  END DO
  veg%froot(:, ms) = 1.0 - veg%froot(:, ms-1)
  DO is = ms-1, 2, -1
    veg%froot(:, is) = veg%froot(:, is)-veg%froot(:,is-1)
  END DO

ENDIF

RETURN
END SUBROUTINE init_veg_from_vegin



SUBROUTINE init_respiration(NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC)
   ! Lestevens 23apr13 - for reading in prog soil & plant resp
   USE cable_um_tech_mod,   ONLY : um1, canopy
   !USE cable_common_module, ONLY : cable_runtime, cable_user

   REAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: NPP_FT_ACC
   REAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: RESP_W_FT_ACC
   REAL, INTENT(INOUT),DIMENSION(um1%land_pts, um1%ntiles) :: RESP_S_ACC

!   REAL, ALLOCATABLE :: tempvar(:,:), tempvar2(:,:)
   INTEGER :: l,j,n

!   ALLOCATE( tempvar(um1%land_pts,um1%ntiles) )
!   ALLOCATE( tempvar2(um1%land_pts,um1%ntiles) )
!
!      DO N=1,um1%NTILES
!         DO J=1,um1%TILE_PTS(N)
!
!            L = um1%TILE_INDEX(j,N)  ! It must be landpt index
!
!            !IF( um1%TILE_FRAC(L,N) .gt. 0.0 ) THEN
!
!               !IF(N <= 13 ) THEN
!                  tempvar(L,N)  = NPP_FT_ACC(L,N)
!                  tempvar2(L,N) = RESP_W_FT_ACC(L,N)
!               !ELSE IF(N > 13 ) THEN
!               !   tempvar(L,N)  = 0.
!               !   tempvar2(L,N) = 0.
!               !ENDIF
!
!            !ENDIF
!
!         ENDDO
!      ENDDO

      !---set soil & plant respiration (now in dim(land_pts,ntiles))
      canopy%frs = PACK(NPP_FT_ACC   , um1%L_TILE_PTS)
      canopy%frp = PACK(RESP_W_FT_ACC, um1%L_TILE_PTS)
      canopy%fnpp = PACK(RESP_S_ACC, um1%L_TILE_PTS)
!      canopy%frs = PACK(tempvar, um1%L_TILE_PTS)
!      canopy%frp = PACK(tempvar2, um1%L_TILE_PTS)

      !---convert units to g C m-2 s-1
      canopy%frs = canopy%frs * 1000.
      canopy%frp = canopy%frp * 1000.
      canopy%fnpp = canopy%fnpp * 1000.

!   DEALLOCATE( tempvar, tempvar2 )

END SUBROUTINE init_respiration

!========================================================================
!========================================================================
!========================================================================

!========================================================================
!========================================================================
!========================================================================
        
SUBROUTINE initialize_radiation( sw_down, lw_down, cos_zenith_angle,           &
                                 surf_down_sw, sin_theta_latitude, ls_rain,    &
                                 ls_snow, tl_1, qw_1, vshr_land, pstar,        &
! rml 2/7/13 pass 3d co2 through to cable if required
                   CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE )   

   USE cable_def_types_mod, ONLY : mp
   USE cable_other_constants_mod, ONLY: RAD_thresh
   USE cable_um_tech_mod,   ONLY : um1, rad, soil, met,                        & 
                                   conv_rain_prevstep, conv_snow_prevstep
   USE cable_common_module, ONLY : cable_runtime, cable_user

   REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: sw_down
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::                    & 
      lw_down,           &
      sin_theta_latitude
   
   REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: cos_zenith_angle

   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows, 4) :: surf_down_sw 
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::                    & 
      ls_rain,    &
      ls_snow,    &   
      tl_1,       &
      qw_1,       &
      vshr_land,  & 
      pstar
   
   REAL, INTENT(IN) :: co2_mmr
! rml 2/7/13 Extra atmospheric co2 variables
   LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
   INTEGER, INTENT(IN) ::                              &
      CO2_DIM_LEN                                      &
     ,CO2_DIM_ROW
   REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio
             
   !___defs 1st call to CABLE in this run. OK in UM & coupled
   LOGICAL, SAVE :: first_call= .TRUE.
     
      IF( first_call ) THEN
         rad%albedo_T = soil%albsoil(:,1)
         first_call = .FALSE.
         ALLOCATE( conv_rain_prevstep(mp), conv_snow_prevstep(mp) )
         conv_rain_prevstep = 0. 
         conv_snow_prevstep = 0.
         met%doy = 0.0
      ENDIF   
      
      ! re-set UM rad. forcings to suit CABLE. also called in explicit call to 
      ! CABLE from subr cable_um_expl_update() 
      CALL update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
      
      ! set met. and rad. forcings to CABLE. also called in radiation call to 
      ! CABLE from subr cable_rad_() !jhan?
      ! subr.  um2cable_met_rad_alb() USES CABLE types met%, rad%, soil%
      ! and kblum% rad. calculated in  update_kblum_radiation() above 
      CALL um2cable_met_rad( cos_zenith_angle)
         
      ! UM met forcing vars needed by CABLE which have UM dimensions
      !(row_length,rows)[_rr], which is no good to cable. These have to be 
      ! re-packed in a single vector of active tiles. Hence we use 
      ! conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
      ! if the land point is/has an active tile
      ! generic format:
      ! um2cable_rr( UM var, CABLE var)
      
      ! CABLE met type forcings, not set by um2cable_met_rad()
      CALL um2cable_rr( LW_DOWN, met%fld)
      CALL um2cable_rr( (LS_RAIN*um1%TIMESTEP), met%precip)
      CALL um2cable_rr( (LS_SNOW*um1%TIMESTEP), met%precip_sn)
      CALL um2cable_rr( TL_1, met%tk)
      CALL um2cable_rr( QW_1, met%qv)
      CALL um2cable_rr( VSHR_LAND, met%ua)
      CALL um2cable_rr( PSTAR*0.01, met%pmb)
   
      !---re-set some of CABLE's forcing variables
      met%precip   =  met%precip + met%precip_sn 
      !met%precip   =  (met%precip + conv_rain_prevstep) &
      !               + (met%precip_sn +  conv_snow_prevstep)
      !               + (met%precip_sn +  conv_rain_prevstep)
      met%tvair =     met%tk
      met%tvrad =     met%tk
      met%coszen =    max(met%coszen,1e-8) 

      !---this is necessary clobrring at present 
      WHERE(met%ua < 0.001 ) met%ua = 0.001
      
      ! rml 24/2/11 Set atmospheric CO2 seen by cable to CO2_MMR (value seen 
      ! by radiation scheme).  Option in future to have cable see interactive 
      ! (3d) CO2 field Convert CO2 from kg/kg to mol/mol ( m_air, 
      ! 28.966 taken from include/constant/ccarbon.h file )
      ! rml 2/7/13 Add in co2_interactive option
      IF (L_CO2_INTERACTIVE) THEN
        CALL um2cable_rr(CO2_3D, met%ca)
      ELSE
        met%ca = CO2_MMR
        !met%ca = 4.314801E-04 ! hard coded to 1850 level, only for 1% increase run (CO2 increase seen by radiation only)
      ENDIF
      met%ca = met%ca * 28.966/44.

      WHERE (met%coszen < RAD_THRESH ) 
         rad%fbeam(:,1) = REAL(0) 
         rad%fbeam(:,2) = REAL(0) 
         rad%fbeam(:,3) = REAL(0) 
      ENDWHERE

      !--- CABLE radiation type forcings, not set by um2cable_met_rad(
      !--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
      CALL um2cable_rr( um1%longitude,rad%longitude )

END SUBROUTINE initialize_radiation

!========================================================================
!========================================================================
!========================================================================
          
SUBROUTINE initialize_canopy(canopy_tile)
   USE cable_um_tech_mod,   ONLY : um1, canopy 
   USE cable_common_module, ONLY : cable_runtime, cable_user
   
   REAL, INTENT(IN),DIMENSION(um1%land_pts, um1%ntiles) :: canopy_tile
   
   ! defs 1st call to CABLE in this run. OK in UM & coupled
   LOGICAL, SAVE :: first_call= .TRUE.
   
      !--- %ga is computed (on LHS only) in define_canopy and 
      !--- then used in soilsnow() in implicit call, then unpacked
      IF( first_call ) THEN
         canopy%ga = 0.
         canopy%us = 0.01
         canopy%fes_cor = 0.
         canopy%fhs_cor = 0.
         ! rml 20/11/2024 may be needed to use fwsoil_switch=Haverd2013
         canopy%fwsoil = 1.
         first_call = .FALSE.
      ENDIF
         
     !---set canopy storage (already in dim(land_pts,ntiles) ) 
     canopy%cansto = pack(CANOPY_TILE, um1%l_tile_pts)

END SUBROUTINE initialize_canopy

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE initialize_roughness( z1_tq, z1_uv, htveg )  
   USE cable_um_tech_mod,   ONLY : um1, rough, veg
   USE cable_common_module, ONLY : ktau_gl
   USE cable_def_types_mod, ONLY : mp
   USE cable_common_module, ONLY : cable_runtime, cable_user
   
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) ::  z1_tq, z1_uv
   REAL, INTENT(INOUT), DIMENSION(um1%land_pts, um1%ntiles) :: htveg
   INTEGER :: i,j,k,L,n
   REAL, ALLOCATABLE, DIMENSION(:,:) :: jhruff, jhwork

      !--- CABLE roughness type forcings
      CALL um2cable_rr( Z1_TQ, rough%za_tq)
      CALL um2cable_rr( Z1_UV, rough%za_uv)

      ALLOCATE(jhwork (um1%land_pts,um1%ntiles) ) 
      ALLOCATE(jhruff (um1%land_pts,um1%ntiles) ) 

      !Veg height changes seasonally in MOSES hence no updates here due to snow
      jhwork = 0.
      DO N=1,um1%NTILES
        DO K=1,um1%TILE_PTS(N)
          I = um1%TILE_INDEX(K,N)
          jhWORK(I,N) = MAX(.01,HTVEG(I,N))
        ENDDO
      ENDDO

      jHRUFF= 0.01 
      DO l=1,um1%land_pts
        DO n=1,um1%ntiles     
          IF( jHRUFF(L,N) .lt. jhwork(l,n)) jHRUFF(L,:) =  jhwork(l,n)
        ENDDO
      ENDDO
      
      rough%hruff= MAX(0.01,veg%hc)
      rough%hruff_grmx = pack(jHRUFF, um1%l_tile_pts) 

      DEALLOCATE( jhruff, jhwork ) 

END SUBROUTINE initialize_roughness

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
   USE cable_um_tech_mod!, only : um1, um_rad, kblum_rad
  
   REAL, INTENT(INOUT), DIMENSION(um1%row_length, um1%rows) :: sw_down
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: cos_zenith_angle
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows, 4) :: surf_down_sw 

      !jhan: do you really want to be changing sw_down            
      SW_DOWN = ( surf_down_sw(:,:,1)                                          &
                        + surf_down_sw(:,:,2)                                  &
                        + surf_down_sw(:,:,3)                                  &
                        + surf_down_sw(:,:,4) )                                &
                        * cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_DIR = ( surf_down_sw(:,:,1)                            &
                        + surf_down_sw(:,:,3) )                                &
                        * cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_DIF = ( surf_down_sw(:,:,2)                            & 
                              + surf_down_sw(:,:,4) )                          &
                              *cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_VIS = (surf_down_sw(:,:,1)                             & 
                              + surf_down_sw(:,:,2) )                          &
                              * cos_zenith_angle(:,:)

      kblum_rad%SW_DOWN_NIR = ( surf_down_sw(:,:,3)                            &
                              + surf_down_sw(:,:,4) )                          &
                              *cos_zenith_angle(:,:)
      ! fbeam for VIS
      kblum_rad%FBEAM(:,:,1) = surf_down_sw(:,:,1)                             &
                              * cos_zenith_angle(:,:)                          &
                                 / max( 0.1, kblum_rad%SW_DOWN_VIS )
      ! fbeam for NIR
      kblum_rad%FBEAM(:,:,2) = surf_down_sw(:,:,3)                             &
                              * cos_zenith_angle(:,:)                          &
                              / max( 0.1, kblum_rad%SW_DOWN_NIR )
      !---fbeam for all solar 
      kblum_rad%FBEAM(:,:,3) = kblum_rad%SW_DOWN_DIR /                         &
                              MAX( 0.1, SW_DOWN )
       
END SUBROUTINE Update_kblum_radiation

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE  um2cable_met_rad( cos_zenith_angle)
   USE cable_um_tech_mod, ONLY :um1, kblum_rad, rad, met

   !___ from UM, cosine zenith angle and soil albedo
   REAL, INTENT(INOUT) :: cos_zenith_angle(um1%row_length, um1%rows)

      !--- CABLE met type forcings
      CALL um2cable_rr( cos_zenith_angle, met%coszen)
      CALL um2cable_rr( kblum_rad%SW_DOWN_VIS, met%fsd(:,1))
      CALL um2cable_rr( kblum_rad%SW_DOWN_NIR, met%fsd(:,2))
      
      !--- CABLE radiation type forcings
      !--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
      CALL um2cable_rr( kblum_rad%FBEAM(:,:,1), rad%fbeam(:,1))
      CALL um2cable_rr( kblum_rad%FBEAM(:,:,2), rad%fbeam(:,2))
      CALL um2cable_rr( kblum_rad%FBEAM(:,:,3), rad%fbeam(:,3))

END SUBROUTINE  um2cable_met_rad

!========================================================================
!========================================================================
!========================================================================

!--- UM met forcing vars needed by CABLE commonly have UM dimensions
!---(row_length,rows), which is no good to CABLE. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE um2cable_rr(umvar,cablevar)
   USE cable_def_types_mod, ONLY : mp
   USE cable_um_tech_mod,   ONLY :um1
 
   REAL, INTENT(IN), DIMENSION(um1%row_length, um1%rows) :: umvar   
   REAL, INTENT(INOUT), DIMENSION(mp) :: cablevar
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: fvar   
   INTEGER :: n,k,l,j,i

      fvar = 0.0
      DO N=1,um1%NTILES                     
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
            I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
            fvar(L,N) = umvar(I,J)
         ENDDO
      ENDDO
      cablevar =  pack(fvar,um1%l_tile_pts)

END SUBROUTINE um2cable_rr


!========================================================================
!========================================================================
!========================================================================

!--- UM met forcing vars needed by CABLE which have UM dimensions
!---(land_points)[_lp], which is no good to cable. These have to be 
!--- re-packed in a single vector of active tiles. Hence we use 
!--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!--- if the land point is/has an active tile
SUBROUTINE um2cable_lp(umvar, defaultin, cablevar, soiltype, skip )
   USE cable_def_types_mod, ONLY : mp, mstype
   USE cable_um_tech_mod,   ONLY :um1
  
   REAL, INTENT(IN), DIMENSION(um1%land_pts) :: umvar
   REAL, INTENT(IN), DIMENSION(mstype) :: defaultin    
   REAL, INTENT(INOUT), DIMENSION(mp) :: cablevar
   INTEGER, INTENT(INOUT), DIMENSION(mp) :: soiltype
   REAL, DIMENSION(:,:), ALLOCATABLE:: fvar   
   LOGICAL, OPTIONAL :: skip
   INTEGER :: n,k,l,i

         
      ALLOCATE( fvar(um1%land_pts,um1%ntiles) )
      fvar = 0.0

      ! loop over Ntiles
      DO N=1,um1%NTILES
         ! loop over number of points per tile
         DO K=1,um1%TILE_PTS(N)
            ! index of each point per tile in an array of dim=(land_pts,ntiles)
            L = um1%TILE_INDEX(K,N)
            ! at this point fvar=umvar, ELSE=0.0 
            fvar(L,N) = umvar(L)
            ! unless explicitly SKIPPED by including arg in subr call
            IF(.NOT. PRESENT(skip) ) THEN
               ! on perma frost tile, set fvar=defaultin
               IF( N == um1%ntiles ) THEN
                  fvar(L,N) =  defaultin(9)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
     
      cablevar     =  PACK(fvar,um1%l_tile_pts)
  
      ! unless explicitly SKIPPED by including arg in subr call
      IF(.NOT. PRESENT(skip) ) THEN
         DO i=1,mp
            ! soiltype=9 for perma-frost tiles 
            IF(soiltype(i)==9) cablevar(i) =  defaultin(9)         
         ENDDO        
      ENDIF
   
      DEALLOCATE(fvar)

END SUBROUTINE um2cable_lp
 

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE init_bgc_vars() 
  
USE cable_def_types_mod,  ONLY: ncs, ncp, mp  
USE cable_um_tech_mod,    ONLY: bgc, veg   
USE cable_pft_params_mod, ONLY: vegin

IMPLICIT NONE
   
INTEGER :: h

DO h=1,mp
  bgc%cplant(h,1) = vegin%cplant1( veg%iveg(h) )
  bgc%cplant(h,2) = vegin%cplant2( veg%iveg(h) )
  bgc%cplant(h,3) = vegin%cplant3( veg%iveg(h) )

  bgc%csoil(h,1)  = vegin%csoil1( veg%iveg(h) )
  bgc%csoil(h,2)  = vegin%csoil2( veg%iveg(h) )
END DO 

! jhan: previous implementation was very mixed up version of what is already a 
! shortcut intended to assume ratecp, ratecs are the same for all PFTs 
bgc%ratecp(1)   = vegin%ratecp1(1)
bgc%ratecp(2)   = vegin%ratecp2(1)
bgc%ratecp(3)   = vegin%ratecp3(1)

bgc%ratecs(1)   = vegin%ratecs1(1)
bgc%ratecs(2)   = vegin%ratecs2(1)
   
RETURN
END SUBROUTINE init_bgc_vars

!========================================================================
!========================================================================
!========================================================================

subroutine init_sumflux_zero() 
   USE cable_um_tech_mod, ONLY : sum_flux
      sum_flux%sumpn = 0.; sum_flux%sumrp = 0.; sum_flux%sumrpw = 0.
      sum_flux%sumrpr = 0.; sum_flux%sumrs = 0.; sum_flux%sumrd = 0.
      sum_flux%dsumpn = 0.; sum_flux%dsumrp = 0.; sum_flux%dsumrs = 0.
      sum_flux%dsumrd = 0.; sum_flux%sumxrp = 0.;  sum_flux%sumxrs = 0.
END SUBROUTINE init_sumflux_zero 

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE alloc_cable_types()
   USE cable_def_types_mod, ONLY : mp, alloc_cbm_var
   USE cable_um_tech_mod,   ONLY : air, canopy, met, bal, rad, rough,          &
                                   soil, ssnow, sum_flux, veg, bgc
  USE allocate_veg_params_mod,  ONLY : allocate_veg_parameter_type
  USE allocate_soil_params_mod, ONLY : allocate_soil_parameter_type
USE cable_canopy_type_mod, ONLY : alloc_canopy_type
USE cable_soil_snow_type_mod, ONLY : alloc_soil_snow_type

      CALL alloc_cbm_var(air, mp)
      CALL alloc_cbm_var(met, mp)
      CALL alloc_cbm_var(bal, mp)
      CALL alloc_cbm_var(rad, mp)
      CALL alloc_cbm_var(rough, mp)
      CALL alloc_cbm_var(sum_flux, mp)
      CALL alloc_cbm_var(bgc, mp)
    CALL allocate_veg_parameter_type(veg, mp)
    CALL allocate_soil_parameter_type(soil, mp)
    CALL alloc_canopy_type(canopy, mp)
    CALL alloc_soil_snow_type(ssnow, mp)

END SUBROUTINE alloc_cable_types

!========================================================================
!========================================================================
!========================================================================


END MODULE cable_um_init_subrs_mod




