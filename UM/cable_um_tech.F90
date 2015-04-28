!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
! 
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located 
! in each directory containing CABLE code.
!
! ==============================================================================
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
   
   USE cable_def_types_mod
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

   INTERFACE check_nmlvar 
      MODULE PROCEDURE check_chvar, check_intvar, check_lgvar
   END INTERFACE check_nmlvar 
 
      TYPE(derived_rad_bands), SAVE :: kblum_rad    
      TYPE(derived_veg_pars),  SAVE :: kblum_veg    
      TYPE(um_dimensions),     SAVE :: um1
      
      REAL,ALLOCATABLE, DIMENSION(:) :: conv_rain_prevstep, conv_snow_prevstep 

CONTAINS

!========================================================================
!========================================================================
!========================================================================

SUBROUTINE cable_um_runtime_vars(runtime_vars_file) 
   USE cable_common_module, ONLY : cable_runtime, cable_user, filename,        &
                                   cable_user, knode_gl, redistrb, wiltParam,  &
                                   satuParam, l_casacnp, l_laiFeedbk,          &
                                   l_vcmaxFeedbk
   USE casavariable, ONLY : casafile
   USE casadimension, ONLY : icycle


   CHARACTER(LEN=*), INTENT(IN) :: runtime_vars_file
   INTEGER :: funit=88
   
   !--- namelist for CABLE runtime vars, files, switches 
   NAMELIST/CABLE/filename, l_casacnp, l_laiFeedbk, l_vcmaxFeedbk, icycle,   &
                  casafile, cable_user, redistrb, wiltParam, satuParam

      !--- assume namelist exists. no iostatus check 
      OPEN(unit=funit,FILE= runtime_vars_file)
         READ(funit,NML=CABLE)
         IF( knode_gl==0)  THEN
            PRINT *, '  '; PRINT *, 'CABLE_log:' 
            PRINT *, '  Opened file - '
            PRINT *, '  ', trim(runtime_vars_file)
            PRINT *, '  for reading runtime vars.' 
            PRINT *, 'End CABLE_log:'; PRINT *, '  '
        ENDIF
      CLOSE(funit)

      if (knode_gl==0) then
        print *, '  '; print *, 'CASA_log:'
        print *, '  icycle =',icycle
        print *, '  l_casacnp =',l_casacnp
        print *, '  l_laiFeedbk =',l_laiFeedbk
        print *, '  l_vcmaxFeedbk =',l_vcmaxFeedbk
        print *, 'End CASA_log:'; print *, '  '
      endif
      IF (l_casacnp  .AND. (icycle == 0 .OR. icycle > 3)) &
          STOP 'CASA_log: icycle must be 1 to 3 when using casaCNP'
      IF ((.NOT. l_casacnp)  .AND. (icycle >= 1)) &
          STOP 'CASA_log: icycle must be <=0 when not using casaCNP'
      IF ((l_laiFeedbk .OR. l_vcmaxFeedbk) .AND. (.NOT. l_casacnp)) &
          STOP 'CASA_log: casaCNP required to get prognostic LAI or Vcmax'
      IF (l_vcmaxFeedbk .AND. icycle < 2) &
          STOP 'CASA_log: icycle must be 2 to 3 to get prognostic Vcmax'
   
      !--- check value of variable 
      CALL check_nmlvar('filename%veg', filename%veg)
      CALL check_nmlvar('filename%soil', filename%soil)
      CALL check_nmlvar('cable_user%DIAG_SOIL_RESP', cable_user%DIAG_SOIL_RESP)
      CALL check_nmlvar('cable_user%LEAF_RESPIRATION',                         &
                        cable_user%LEAF_RESPIRATION)
      CALL check_nmlvar('cable_user%FWSOIL_SWITCH', cable_user%FWSOIL_SWITCH)
      CALL check_nmlvar('cable_user%RUN_DIAG_LEVEL', cable_user%RUN_DIAG_LEVEL)
      CALL check_nmlvar('cable_user%l_new_roughness_soil',                     &
                         cable_user%l_new_roughness_soil)
      CALL check_nmlvar('cable_user%l_new_roughness_soil',                     &
                         cable_user%l_new_roughness_soil)
      CALL check_nmlvar('cable_user%l_new_roughness_soil',                     &
                         cable_user%l_new_roughness_soil)

END SUBROUTINE cable_um_runtime_vars

!jhan: also add real, logical, int interfaces
SUBROUTINE check_chvar(this_var, val_var)
   USE cable_common_module, ONLY : knode_gl

   CHARACTER(LEN=*), INTENT(IN) :: this_var, val_var 
   
      IF (knode_gl==0) THEN
         PRINT *, '  '; PRINT *, 'CABLE_log:' 
         PRINT *, '   run time variable - '
         PRINT *, '  ', trim(this_var) 
         PRINT *, '   defined as - '
         PRINT *, '  ', trim(val_var) 
         PRINT *, 'End CABLE_log:'; PRINT *, '  '
      ENDIf

END SUBROUTINE check_chvar

SUBROUTINE check_intvar(this_var, val_var)
   USE cable_common_module, ONLY : knode_gl

   CHARACTER(LEN=*), INTENT(IN) :: this_var
   INTEGER, INTENT(IN) :: val_var 

      IF (knode_gl==0) THEN
         PRINT *, '  '; PRINT *, 'CABLE_log:' 
         PRINT *, '   run time variable - '
         PRINT *, '  ', trim(this_var) 
         PRINT *, '   defined as - '
         PRINT *, '  ', val_var
         PRINT *, 'End CABLE_log:'; PRINT *, '  '
      ENDIF

END SUBROUTINE check_intvar

SUBROUTINE check_lgvar(this_var, val_var)
   USE cable_common_module, ONLY : knode_gl

   CHARACTER(LEN=*), INTENT(IN) :: this_var
   LOGICAL, INTENT(IN) :: val_var

      IF (knode_gl==0) THEN
         PRINT *, '  '; PRINT *, 'CABLE_log:'
         PRINT *, '   run time variable - '
         PRINT *, '  ', trim(this_var)
         PRINT *, '   defined as - '
         PRINT *, '  ', (val_var)
         PRINT *, 'End CABLE_log:'; PRINT *, '  '
      ENDIf

END SUBROUTINE check_lgvar
    
!========================================================================= 
!=========================================================================
!========================================================================= 
 
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

SUBROUTINE dealloc_vegin_soilin()
   USE cable_common_module, ONLY : cable_runtime, cable_user, vegin, soilin
      
      DEALLOCATE(vegin%canst1)
      DEALLOCATE(vegin%dleaf)
      DEALLOCATE(vegin%vcmax)
      DEALLOCATE(vegin%ejmax)
      DEALLOCATE(vegin%hc)
      DEALLOCATE(vegin%xfang)
      DEALLOCATE(vegin%rp20)
      DEALLOCATE(vegin%rpcoef)
      DEALLOCATE(vegin% rs20)
      DEALLOCATE(vegin%shelrb)
      DEALLOCATE(vegin%vegcf)
      DEALLOCATE(vegin%frac4)
      DEALLOCATE(vegin%refl)
      DEALLOCATE(vegin%taul)
      DEALLOCATE(vegin%xalbnir)
      DEALLOCATE(vegin%extkn)
      DEALLOCATE(vegin%froot)
      DEALLOCATE(vegin%tminvj)
      DEALLOCATE(vegin%tmaxvj)
      DEALLOCATE(vegin%vbeta)
      DEALLOCATE(vegin%cplant)
      DEALLOCATE(vegin%csoil)
      DEALLOCATE(vegin%ratecp)
      DEALLOCATE(vegin%ratecs)
     
      DEALLOCATE(soilin%silt)
      DEALLOCATE(soilin%clay)
      DEALLOCATE(soilin%sand)
      DEALLOCATE(soilin%swilt)
      DEALLOCATE(soilin%sfc)
      DEALLOCATE(soilin%ssat)
      DEALLOCATE(soilin%bch)
      DEALLOCATE(soilin%hyds)
      DEALLOCATE(soilin%sucs)
      DEALLOCATE(soilin%rhosoil)
      DEALLOCATE(soilin%css)

END SUBROUTINE dealloc_vegin_soilin


   !========================================================================
   !========================================================================
   !========================================================================


END MODULE cable_um_tech_mod




