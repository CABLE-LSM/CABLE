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
! Purpose: Routines to read CABLE namelist, check variables, allocate and 
!          deallocate CABLE arrays
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Rewrite of code from v1.8 (ACCESS1.3)
!          soil_snow_type now ssnow (instead of ssoil)
!
!
! ==============================================================================

MODULE cable_um_tech_mod
   
USE cable_def_types_mod, ONLY : air_type, bgc_pool_type, met_type,             &
                 balances_type, radiation_type, roughness_type, sum_flux_type, &
                 soil_snow_type, canopy_type, veg_parameter_type,              &
                 soil_parameter_type, climate_type

   IMPLICIT NONE

   TYPE(air_type), SAVE             :: air
   TYPE(bgc_pool_type), SAVE        :: bgc
   TYPE(met_type), SAVE             :: met
   TYPE(balances_type), SAVE        :: bal
   TYPE(radiation_type), SAVE       :: rad
   TYPE(roughness_type), SAVE       :: rough
   TYPE(soil_parameter_type), SAVE  :: soil       ! soil parameters
   TYPE(soil_snow_type), SAVE       :: ssnow
   TYPE(sum_flux_type), SAVE        :: sum_flux
   TYPE(veg_parameter_type), SAVE   :: veg        ! vegetation parameters
   TYPE(canopy_type), SAVE          :: canopy
   TYPE(climate_type), SAVE         :: climate

   TYPE derived_rad_bands    
      REAL, ALLOCATABLE ::                                                     &
         SW_DOWN_DIR (:,:), & ! Surface downward SW direct radiation (W/m2).
         SW_DOWN_DIF(:,:), & ! Surface downward SW diffuse radiation (W/m2).
         SW_DOWN_VIS(:,:), & ! Surface downward VIS radiation (W/m2).
         SW_DOWN_NIR(:,:), & ! Surface downward NIR radiation (W/m2).
         FBEAM(:,:,:)      ! Surface downward SW radiation (W/m2).
   END TYPE derived_rad_bands
   
   TYPE um_dimensions 
      INTEGER :: row_length, rows, land_pts, ntiles, npft,                     &
                 sm_levels, timestep 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: tile_pts, land_index
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tile_index
      REAL :: rho_water
      REAL,ALLOCATABLE, DIMENSION(:,:) :: tile_frac
      REAL,ALLOCATABLE, DIMENSION(:,:) :: latitude, longitude
      LOGICAL,ALLOCATABLE, DIMENSION(:,:) :: l_tile_pts
   ENDTYPE um_dimensions 

   TYPE derived_veg_pars
      INTEGER, DIMENSION(:,:), POINTER ::                                      &
         ivegt(:,:),    & ! vegetation  types
         isoilm(:,:)      ! soil types
      REAL, DIMENSION(:,:), POINTER ::                                         &
         htveg(:,:),    &
         laift(:,:)       ! hruffmax(:.:)
   END TYPE derived_veg_pars

      TYPE(derived_rad_bands), SAVE :: kblum_rad    
      TYPE(derived_veg_pars),  SAVE :: kblum_veg    
      TYPE(um_dimensions),     SAVE :: um1
      
      REAL,ALLOCATABLE, DIMENSION(:) :: conv_rain_prevstep, conv_snow_prevstep 

CONTAINS

SUBROUTINE alloc_um_interface_types( row_length, rows, land_pts, ntiles,       &
                                     sm_levels )
      USE cable_common_module, ONLY : cable_runtime, cable_user
      
      INTEGER,INTENT(IN) :: row_length, rows, land_pts, ntiles, sm_levels   

         ALLOCATE( um1%land_index(land_pts) )
         ALLOCATE( um1%tile_pts(ntiles) )
         ALLOCATE( um1%tile_frac(land_pts, ntiles) )
         ALLOCATE( um1%tile_index(land_pts, ntiles) )
         ALLOCATE( um1%latitude(row_length, rows) )
         ALLOCATE( um1%longitude(row_length, rows) )
         ALLOCATE( um1%l_tile_pts(land_pts, ntiles) ) 
        !-------------------------------------------------------
         ALLOCATE( kblum_rad%sw_down_dir(row_length,rows) )
         ALLOCATE( kblum_rad%sw_down_dif(row_length,rows) )
         ALLOCATE( kblum_rad%sw_down_vis(row_length,rows) )
         ALLOCATE( kblum_rad%sw_down_nir(row_length,rows) )
         ALLOCATE( kblum_rad%fbeam(row_length,rows,3) )
         ALLOCATE( kblum_veg%htveg(land_pts,ntiles) )
         ALLOCATE( kblum_veg%laift(land_pts,ntiles) )
         ALLOCATE( kblum_veg%ivegt(land_pts,ntiles) )
         ALLOCATE( kblum_veg%isoilm(land_pts,ntiles) ) 
         
END SUBROUTINE alloc_um_interface_types 

!========================================================================
!========================================================================
!========================================================================

END MODULE cable_um_tech_mod




