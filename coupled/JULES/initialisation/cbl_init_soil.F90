MODULE cable_um_init_soil_mod
   
IMPLICIT NONE

CONTAINS
SUBROUTINE initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,        &
                            smvccl, albsoil, tsoil_tile, sthu, sthu_tile,     &
                            dzsoil, slope_avg, slope_std,                     &
                             dz_gw,aq_perm,drain_dens, soilin, veg_cbl, soil_cbl, ssnow_cbl ) 

USE cable_um_init_subrs_mod, ONLY: um2cable_lp

USE cable_def_types_mod, ONLY: ms, mstype, mp, r_2
USE cable_um_tech_mod,   ONLY: um1
USE cable_common_module, ONLY: cable_runtime, cable_user,                     &
                                knode_gl, gw_params
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soilin_type
USE cable_params_mod,         ONLY: soil_parameter_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
TYPE(veg_parameter_type),   INTENT(INOUT) :: veg_cbl
TYPE(soilin_type),  INTENT(INOUT) ::  soilin  
TYPE(soil_parameter_type),   INTENT(INOUT) :: soil_cbl
TYPE(soil_snow_type),   INTENT(INOUT) :: ssnow_cbl
  
REAL, INTENT(IN), DIMENSION(um1%land_pts) ::                                  &
   bexp,                                                                      &
   hcon,                                                                      &
   satcon,                                                                    &
   sathh,                                                                     &
   smvcst,                                                                    &
   smvcwt,                                                                    &
   smvccl,                                                                    &
   albsoil 
   
REAL, INTENT(IN), DIMENSION(um1%land_pts) ::                                  &
   slope_avg,                                                                 &
   slope_std,                                                                 &
   dz_gw,aq_perm,drain_dens

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%sm_levels) :: sthu
   
REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles, um1%sm_levels) ::       &
   sthu_tile,                                                                 &
   tsoil_tile

REAL, INTENT(IN), DIMENSION(um1%sm_levels) :: dzsoil

!___defs 1st call to CABLE in this run
LOGICAL, SAVE :: first_call= .TRUE.
INTEGER :: i,j,k,l,n
REAL, ALLOCATABLE :: tempvar(:), tempvar2(:),fwork(:,:)
LOGICAL, PARAMETER :: skip =.TRUE. 
REAL, DIMENSION(mstype) :: dummy 
REAL :: tmp_clay, tmp_sand
REAL, ALLOCATABLE :: znode(:), ssat_bounded(:,:),rho_soil_bulk(:,:)

REAL, PARAMETER :: snow_ccnsw = 2.0,                                          &
                   ssat_lo = 0.15,                                            &
                   ssat_hi = 0.65,                                            &
                   rhob_lo = 810.0,                                           &
                   rhob_hi = 2300.0

REAL :: sucs_sign_factor,hyds_unit_factor,sucs_min_magnitude


IF (cable_user%gw_model) THEN
  sucs_sign_factor = 1.0
  hyds_unit_factor = 1.0
  sucs_min_magnitude = 106.0
ELSE
  sucs_sign_factor = 1.0
  hyds_unit_factor = 1.0 / 1000.0
  sucs_min_magnitude = 106.0 / 1000.0
END IF
      
dummy = 0.0 

IF ( first_call ) THEN 

  ssnow_cbl%pudsto = 0.0; ssnow_cbl%pudsmx = 0.0
      
  !--- soil%isoilm defines soiltype. 
  ! currently is either 2 (arbitrarily) or 9.
  ! type 9 -> permanent ice points which are dealt with by CABLE. 
  ! Spatially explicit soil properties are used by 
  ! the UM anyway, and is only really an issue for soil%css & 
  ! soil%rhosoil, which are set to either 2 or 9. 
  ! dealing with this in CASACNP is another issue.
  !--- %isoilm=10 for Lakes
  soil_cbl%isoilm  =  2
            
  ! set soil type for permanent ice based on where permanent ice 
  ! located in vegetation map (in v1.8 set by soil albedo value)
  ! hard-wired numbers to be removed in future release
  WHERE ( veg_cbl%iveg == 17 ) soil_cbl%isoilm = 9
           
  !--- set CABLE-var soil%albsoil from UM var albsoil
  ! (see below ~ um2cable_lp)
  CALL um2cable_lp( albsoil, dummy, soil_cbl%albsoil(:,1),                    &
                    soil_cbl%isoilm, skip )

  !--- defined in soil_thick.h in UM
  soil_cbl%zse = dzsoil
  soil_cbl%zse_vec = SPREAD(dzsoil,1,mp)
         
  ! distance between consecutive layer midpoints
  soil_cbl%zshh(1) = 0.5 * soil_cbl%zse(1) 
  soil_cbl%zshh(ms+1) = 0.5 * soil_cbl%zse(ms)
  soil_cbl%zshh(2:ms) = 0.5 * (soil_cbl%zse(1:ms-1) + soil_cbl%zse(2:ms))

  !node depths
  IF (ALLOCATED(znode)) DEALLOCATE(znode)
  ALLOCATE(znode(ms))

  znode(1) = soil_cbl%zshh(1)
  DO k = 2,ms
    znode(k) = znode(k-1) * 0.5 * (soil_cbl%zse(k-1) + soil_cbl%zse(k))
  END DO

  !-------------------------------------------------------------------
  !--- UM met forcing vars needed by CABLE which have UM dimensions
  !---(land_pts,ntiles)[_lp], which is no good to cable. These have to be 
  !--- re-packed in a single vector of active tiles. Hence we use 
  !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
  !--- if the land point is/has an active tile. generic format:
  !---     um2cable_lp( UM var, 
  !---                 default value for snow tile  where 
  !---                 peranent ice point to be treatedas a snowtile, 
  !---                 CABLE var, 
  !---                 mask )
  !--- where mask tells um2cable_lp whether or not to use default value 
  !--- for snow tile 
  !-------------------------------------------------------------------
         
  ! parameter b in Campbell equation 
  CALL um2cable_lp( bexp, soilin%bch, soil_cbl%bch, soil_cbl%isoilm)
         
  ALLOCATE( tempvar(mstype), tempvar2(mp) )
   
  tempvar = soilin%sand(9) * 0.3  + soilin%clay(9) *0.25 +                    &
            soilin%silt(9) * 0.265
         
  CALL um2cable_lp( hcon, tempvar, tempvar2, soil_cbl%isoilm)
  soil_cbl%cnsd = REAL( tempvar2, r_2 )
  DEALLOCATE( tempvar, tempvar2 )
         
  ! hydraulic conductivity @saturation (satcon[mm/s], soilin%hyds[m/s] )
  CALL um2cable_lp( satcon,soilin%hyds * 1000.0, soil_cbl%hyds, soil_cbl%isoilm)

  CALL um2cable_lp( sathh,  soilin%sucs,  soil_cbl%sucs,  soil_cbl%isoilm)
  CALL um2cable_lp( smvcst, soilin%ssat,  soil_cbl%ssat,  soil_cbl%isoilm)
  CALL um2cable_lp( smvcwt, soilin%swilt, soil_cbl%swilt, soil_cbl%isoilm)
  CALL um2cable_lp( smvccl, soilin%sfc,   soil_cbl%sfc,   soil_cbl%isoilm)
   
    
  !mrd561
  IF (ALLOCATED(fwork)) DEALLOCATE(fwork) 

  ALLOCATE( fwork(um1%land_pts,um1%ntiles) )

  fwork(:,:) = 20.0
  DO n = 1,um1%ntiles
    DO k = 1,um1%tile_pts(n)
      i = um1%tile_index(k,n)
      fwork(i,n) = dz_gw(i)
    END DO
  END DO

  soil_cbl%GWdz(:) = PACK(fwork(:,:),um1%l_tile_pts)

  DO i = 1,mp 
    IF (soil_cbl%GWdz(i)  <   20.0) soil_cbl%GWdz(i) = 20.0
    IF (soil_cbl%GWdz(i)  >   150.0) soil_cbl%GWdz(i) = 150.0
    IF (veg_cbl%iveg(i)  ==  16) soil_cbl%GWdz(i) = 150.0
  END DO

  fwork(:,:) = 0.02
  DO n = 1,um1%ntiles
    DO k = 1,um1%tile_pts(n)
      i = um1%tile_index(k,n)
      fwork(i,n) = slope_avg(i)
    END DO
  END DO
  soil_cbl%slope(:) = PACK(fwork(:,:),um1%l_tile_pts)
  DO i = 1,mp 
    IF (soil_cbl%slope(i)  <   0.0002) soil_cbl%slope(i) = 0.0002
    IF (soil_cbl%slope(i)  >   0.2) soil_cbl%slope(i) = 0.2
    IF (veg_cbl%iveg(i)  ==  16) soil_cbl%slope(i) = 0.0002
  END DO
 
  fwork(:,:) = .005 
  DO n = 1,um1%ntiles
    DO k = 1,um1%tile_pts(n)
      i = um1%tile_index(k,n)
      fwork(i,n) = slope_std(i)
    END DO
  END DO
  soil_cbl%slope_std(:) = PACK(fwork(:,:),um1%l_tile_pts) 
  DO i = 1,mp 
    IF (soil_cbl%slope_std(i)  <   0.00002) soil_cbl%slope_std(i) = 0.00002
    IF (soil_cbl%slope_std(i)  >   0.2) soil_cbl%slope_std(i) = 0.2
    IF (veg_cbl%iveg(i)  ==  16) soil_cbl%slope_std(i) = 0.00002
  END DO
  fwork(:,:) = 0.0008
  soil_cbl%drain_dens(:) = 0.0
  DO n = 1,um1%ntiles
    DO k = 1,um1%tile_pts(n)
      i = um1%tile_index(k,n)
      fwork(i,n) = drain_dens(i)
    END DO
  END DO
  soil_cbl%drain_dens(:) = PACK(fwork(:,:),um1%l_tile_pts)
  DO i = 1,mp 
    IF (soil_cbl%drain_dens(i)  <   1.0e-6) soil_cbl%drain_dens(i) = 1.0e-6
    IF (soil_cbl%drain_dens(i)  >   0.02) soil_cbl%drain_dens(i) = 0.02
    IF (veg_cbl%iveg(i)  ==  16) soil_cbl%drain_dens(i) = 0.02
  END DO

  WRITE(6,*) 'maxval soil%drain_dens',MAXVAL(soil_cbl%drain_dens,dim = 1)
  WRITE(6,*) 'minval soil%drain_dens',MINVAL(soil_cbl%drain_dens,dim = 1)
  IF (ANY(soil_cbl%drain_dens(:)  ==  0.0))                                   &
                    WRITE( * ,*) 'drain_dens has values of zero'

 

  soil_cbl%GWhyds_vec(:) = 0.0
  fwork = 3.0e-6
  DO n = 1,um1%ntiles
    DO k = 1,um1%tile_pts(n)
      i = um1%tile_index(k,n)
      fwork(i,n) = aq_perm(i) / 10.0
    END DO
  END DO
  soil_cbl%GWhyds_vec(:) = PACK(fwork(:,:),um1%l_tile_pts) 
  DO i = 1,mp 
    IF (soil_cbl%GWhyds_vec(i)  <   1.0e-8) soil_cbl%GWhyds_vec(i)  = 1.0e-8
    IF (soil_cbl%GWhyds_vec(i)  >   1.0e-3) soil_cbl%GWhyds_vec(i)  = 1.0e-3
  END DO

  WRITE(6,*) 'maxval soil%GWhyds_vec',MAXVAL(soil_cbl%GWhyds_vec,dim = 1)
  WRITE(6,*) 'minval soil%GWhyds_vec',MINVAL(soil_cbl%GWhyds_vec,dim = 1)
  IF (ANY(soil_cbl%GWhyds_vec(:)  ==  0.0))                                   &
                    WRITE( * ,*) 'hyds_vec has values of zero'

  DEALLOCATE(fwork) 
            
  !--- (re)set values for CABLE
  soil_cbl%ibp2    =  NINT(soil_cbl%bch) + 2
  soil_cbl%i2bp3   =  2 * NINT(soil_cbl%bch) + 3
         
  ! satcon in UM is in mm/s; Cable needs m/s
  soil_cbl%hyds    = hyds_unit_factor * soil_cbl%hyds
  soil_cbl%sucs    = ABS( soil_cbl%sucs) * sucs_sign_factor
  soil_cbl%sucs    =  MAX(sucs_min_magnitude,soil_cbl%sucs)
  soil_cbl%ssat    =  MAX( soil_cbl%ssat, soil_cbl%sfc + 0.01 )

         
  !jhan:coupled runs 
  soil_cbl%hsbh    =  soil_cbl%hyds * ABS(soil_cbl%sucs) * soil_cbl%bch

  WHERE (soil_cbl%ssat > 0.0 )                                                &
     soil_cbl%pwb_min =  (soil_cbl%swilt / soil_cbl%ssat )**soil_cbl%ibp2
           
  !--- these are temporary 
  soil_cbl%rhosoil =   1598.667 ! match to offline gridinfo for Loobos soilin%rhosoil(soil_cbl%isoilm)
  soil_cbl%css     =  803.5883 ! match to offline gridinfo for Loobos soilin%css(soil_cbl%isoilm)

  DO k = 1,ms
    soil_cbl%ssat_vec(:,k)      = REAL(soil_cbl%ssat(:)   ,r_2)    
    soil_cbl%sucs_vec(:,k)      = REAL(soil_cbl%sucs(:)   ,r_2)   
    soil_cbl%hyds_vec(:,k)      = REAL(soil_cbl%hyds(:)   ,r_2)  
    soil_cbl%swilt_vec(:,k)     = REAL(soil_cbl%swilt(:)  ,r_2)  
    soil_cbl%bch_vec(:,k)       = REAL(soil_cbl%bch(:)    ,r_2)
    soil_cbl%sfc_vec(:,k)       = REAL(soil_cbl%sfc(:)    ,r_2)
    soil_cbl%rhosoil_vec(:,k)   = REAL(soil_cbl%rhosoil(:),r_2)   
    soil_cbl%cnsd_vec(:,k)      = REAL(soil_cbl%cnsd      ,r_2)
    soil_cbl%css_vec(:,k)       = REAL(soil_cbl%css       ,r_2)
    soil_cbl%watr(:,k)          = 0.001_r_2
  END DO

  WHERE (soil_cbl%ssat_vec  <=  0.0 .AND. soil_cbl%sfc_vec  >   0.0)
    soil_cbl%ssat_vec = soil_cbl%sfc_vec + 0.05
  END WHERE
  !--- Lestevens 28 Sept 2012 - Fix Init for soil% textures 
  !--- needed for CASA-CNP

  !default values, overwrite if cable_uer%gw_model selected
  soil_cbl%clay = soilin%clay(soil_cbl%isoilm)
  soil_cbl%silt = soilin%silt(soil_cbl%isoilm)
  soil_cbl%sand = soilin%sand(soil_cbl%isoilm)
  DO k = 1,ms
      !should read texture by layer evantually
    soil_cbl%clay_vec(:,k) = soil_cbl%clay(:)
    soil_cbl%sand_vec(:,k) = soil_cbl%sand(:)
    soil_cbl%silt_vec(:,k) = soil_cbl%silt(:)
  END DO
         
  DO k = 1,ms
    DO i = 1,mp

      IF ( (soil_cbl%silt_vec(i,k)  >   0.99) .OR.                            &
           (soil_cbl%silt_vec(i,k)  <   0.01) .OR.                            &
           (soil_cbl%sand_vec(i,k)  >   0.99) .OR.                            &
           (soil_cbl%sand_vec(i,k)  <   0.01) .OR.                            &
           (soil_cbl%clay_vec(i,k)  >   0.99) .OR.                            &
           (soil_cbl%clay_vec(i,k)  <   0.01) ) THEN

           !all bad
        soil_cbl%clay_vec(i,k) = 0.3
        soil_cbl%sand_vec(i,k) = 0.3
        soil_cbl%silt_vec(i,k) = 0.4

      END IF

    END DO

  END DO


  IF (cable_user%gw_model) THEN
           
    DO k = 1,ms

      DO i = 1,mp  !from reversing pedotransfer functions
                 !,ay cause io issues because not passed into um

        IF (soil_cbl%isoilm(i)  /=  9) THEN

          soil_cbl%hyds_vec(i,k) = soil_cbl%hyds_vec(i,k) * &   !change in hyds
                              EXP(-gw_params%hkrz * ( znode(k) - gw_params%zdepth) )
  
        END IF


      END DO
    END DO

    k = 1
    soil_cbl%hyds(:) = soil_cbl%hyds_vec(:,k)

  END IF

  IF (cable_user%soil_thermal_fix) THEN

    IF (ALLOCATED(ssat_bounded)) DEALLOCATE(ssat_bounded)
    IF (ALLOCATED(rho_soil_bulk)) DEALLOCATE(rho_soil_bulk)

    ALLOCATE(ssat_bounded(SIZE(soil_cbl%ssat_vec,dim = 1),                    &
                          SIZE(soil_cbl%ssat_vec,dim = 2) ) )

    ssat_bounded(:,:) = MIN( ssat_hi, MAX(ssat_lo,                            &
                                       soil_cbl%ssat_vec(:,:) ) )

    ALLOCATE(rho_soil_bulk(SIZE(soil_cbl%rhosoil_vec,dim = 1),                &
                           SIZE(soil_cbl%rhosoil_vec,dim = 2) ) )

    rho_soil_bulk(:,:) = MIN(rhob_hi, MAX(rhob_lo ,                           &
                           (2700.0 * (1.0 - ssat_bounded(:,:)) ) ) )


    DO k = 1,ms
      DO i = 1,mp


        IF (soil_cbl%isoilm(i)  /=  9) THEN

          soil_cbl%rhosoil_vec(i,k) = 2700.0
          soil_cbl%rhosoil_vec(i,k) =  1598.667 ! match to offline gridinfo for Loobos soilin%rhosoil(soil_cbl%isoilm)

          soil_cbl%cnsd_vec(i,k) = ( (0.135 * (1.0 - ssat_bounded(i,k))) +    &
                              (64.7 / rho_soil_bulk(i,k)) ) /                 &
                            (1.0-0.947 * (1.0 - ssat_bounded(i,k)))

        END IF

      END DO
    END DO

    k = 1
    DO i = 1,mp
      IF (soil_cbl%isoilm(i)  /=  9) THEN
        soil_cbl%rhosoil(i) = soil_cbl%rhosoil_vec(i,1)
        soil_cbl%cnsd(i)    = soil_cbl%cnsd_vec(i,1)
      END IF
    END DO

    IF (ALLOCATED(ssat_bounded)) DEALLOCATE(ssat_bounded)
    IF (ALLOCATED(rho_soil_bulk)) DEALLOCATE(rho_soil_bulk)

  END IF

  !always set these though not needed unless gw_model - true
  !should read in the values but they need calibration
  !where (soil%GWssat_vec .lt. 0.0) &
  WHERE (soil_cbl%isoilm  ==  9 .OR. veg_cbl%iveg  ==  16)
    !CAN leave these read in, not enough testing for now
    soil_cbl%GWhyds_vec = soil_cbl%hyds_vec(:,ms)
  END WHERE
  soil_cbl%GWssat_vec = soil_cbl%ssat_vec(:,ms)
  soil_cbl%GWsucs_vec = soil_cbl%sucs_vec(:,ms)
  soil_cbl%GWbch_vec  = soil_cbl%bch_vec(:,ms)
  soil_cbl%GWwatr     = 0.0

         
  !for sli   
  soil_cbl%nhorizons = 1 ! use 1 soil horizon globally
            
  first_call= .FALSE.
END IF

END SUBROUTINE initialize_soil
 
END MODULE cable_um_init_soil_mod
