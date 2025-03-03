MODULE cbl_um_init_soil_mod
   
IMPLICIT NONE
PUBLIC initialize_soil

CONTAINS

SUBROUTINE initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,         &
                        smvccl, albsoil,                             & 
                        hcap, soil_clay, soil_silt, soil_sand,       &
                        tsoil_tile, sthu, sthu_tile, dzsoil )
                        
! subrs
USE cable_pack_mod,    ONLY: pack_landpts2mp_ICE, pack_landpts2mp, cable_pack_rr

! data
!USE params_io_mod_cbl,      ONLY: params_io_data_type
USE cable_common_module,     ONLY : cable_user, gw_params
USE cable_def_types_mod,     ONLY: soil_parameter_type, r_2
USE cable_def_types_mod,     ONLY: ms, mstype, mp, r_2
USE cable_um_tech_mod,       ONLY: um1, soil, veg, ssnow 
USE cable_soil_params_mod,   ONLY: soilin
USE grid_constants_mod_cbl,  ONLY: ICE_SoilType, nsoil_max 
USE cable_surface_types_mod, ONLY: ICE_SurfaceType => ice_cable
IMPLICIT NONE

REAL, INTENT(IN), DIMENSION(um1%land_pts) :: &
   bexp, &
   hcon, &
   satcon, & 
   sathh, &
   smvcst, &
   smvcwt, &
   smvccl, &
   albsoil 

REAL, INTENT(IN) :: soil_clay ( um1%row_length, um1%rows )
REAL, INTENT(IN) :: soil_silt ( um1%row_length, um1%rows )
REAL, INTENT(IN) :: soil_sand ( um1%row_length, um1%rows )
REAL, INTENT(IN) :: hcap(um1%land_pts)

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%sm_levels) :: sthu

REAL, INTENT(IN), DIMENSION(um1%land_pts, um1%ntiles, um1%sm_levels) :: &
   sthu_tile,     &
   tsoil_tile

REAL, INTENT(IN), DIMENSION(um1%sm_levels) :: dzsoil

!___defs 1st call to CABLE in this run
LOGICAL, SAVE :: first_call= .TRUE.
INTEGER :: i,j,k,L,n
REAL, ALLOCATABLE :: tempvar(:), tempvar2(:)
LOGICAL, PARAMETER :: skip =.TRUE. 
REAL, DIMENSION(mstype) :: dummy 
!from AM3
REAL    :: hcon_ICE(nsoil_max) 
REAL    :: soilCnsd_real(mp) 
REAL    :: sucs_sign_factor, hyds_unit_factor, sucs_min_magnitude

soil%heat_cap_lower_limit = 0.01

soil%isoilm(:)=2
DO j=1,mp
  IF (veg%iveg(j) == ICE_SurfaceType) THEN
    soil%isoilm(j) = ICE_SoilType
  END IF
END DO

soil%zse = dzsoil  !ESM!

! distance between consecutive layer midpoints
soil%zshh(1)    = 0.5 * soil%zse(1) 
soil%zshh(ms+1) = 0.5 * soil%zse(ms)
soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))

ssnow%pudsto = 0.0; ssnow%pudsmx = 0.0!ESM!

! albsoil ->  soil%albsoil
CALL pack_landpts2mp( um1%ntiles, um1%land_pts, mp, um1%tile_pts, um1%tile_index,            &
                      um1%L_tile_pts, albsoil, soil%albsoil(:,1) )
  
! bexp -> soil%bch: parameter b in Campbell equation 
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, BEXP, soil%isoilm,  &
                          soilin%bch, soil%bch )
!--- Fix Init for soil% textures  needed for CASA-CNP
CALL cable_pack_rr( soil%sand, soil_sand, mp, um1%l_tile_pts, um1%row_length,    &
                    um1%rows, um1%ntiles, um1%land_pts, um1%land_index, um1%tile_pts, &  
                    um1%tile_index )

CALL cable_pack_rr( soil%clay, soil_clay, mp, um1%l_tile_pts, um1%row_length,    &
                    um1%rows, um1%ntiles, um1%land_pts, um1%land_index, um1%tile_pts, &  
                    um1%tile_index )

CALL cable_pack_rr( soil%silt, soil_silt, mp, um1%l_tile_pts, um1%row_length,    &
                    um1%rows, um1%ntiles, um1%land_pts, um1%land_index, um1%tile_pts, &  
                    um1%tile_index )

!jhan: where do these wweighting factors come from
hcon_ICE(:) = ( soilin%sand(ICE_SoilType) * 0.3 )                         &
             + ( soilin%clay(ICE_SoilType) * 0.25 )                       &
             + ( soilin%silt(ICE_SoilType) * 0.265 )

! hcon -> soil%cnsd
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, hcon, soil%isoilm,  &
                          hcon_ICE, soilCnsd_real )

soil%cnsd = REAL( soilCnsd_real, r_2 )

! CABLE soil parameters are in general PACKed from UM spatial fields
! *EXCEPT* for permanent ice points where we overwrite with parametrs read from
! cable_soil_params (pars%)

! soil%hyds = hydraulic conductivity @saturation is PACKed from satcon[mm/s]
! *EXCEPT* at permanent ice points where we overwrite with parametrs read from
! pars%soilin_hyds[m/s] which are already in correct units. Dividing by 1000  
! soil%hyds is correct for tiles 1:8. Would be 1000* too small, hence *1000 first 
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, satcon, soil%isoilm,&
                          soilin%hyds * 1000.0, soil%hyds )
                          !pars%soilin_hyds * 1000.0, soil%hyds )
soil%hyds    = soil%hyds / 1000.0

! sathh -> soil%sucs
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, sathh, soil%isoilm, &
                          soilin%sucs, soil%sucs )

! smvcst -> soil%ssat
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, smvcst, soil%isoilm,&
                          soilin%ssat, soil%ssat )

! smvcwt -> soil%swilt
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, smvcwt, soil%isoilm,&
                          soilin%swilt, soil%swilt )

! smvccl -> soil%sfc  
CALL pack_landpts2mp_ICE( um1%ntiles, um1%land_pts, mp, nsoil_max, ICE_soiltype,       &
                          um1%tile_pts, um1%tile_index, um1%L_tile_pts, smvccl, soil%isoilm,&
                          soilin%sfc, soil%sfc )

!--- (re)set values for CABLE
soil%ibp2    =  NINT(soil%bch)+2
soil%i2bp3   =  2*NINT(soil%bch)+3

soil%ssat    = MAX( soil%ssat, soil%sfc + 0.01 )
soil%sucs    = ABS( soil%sucs )
sucs_min_magnitude = 106.0/1000.0 ! satcon in UM is in mm/s; Cable needs m/s
soil%sucs          = MAX(sucs_min_magnitude,soil%sucs)
soil%hsbh          = soil%hyds * ABS(soil%sucs) * soil%bch
WHERE( soil%ssat > 0. ) 
  soil%pwb_min =  (soil%swilt / soil%ssat )**soil%ibp2
END WHERE

! USED product (rhosoil*css) assumes css is specific heat capacity 
! JULES hcap is already volumetric heat capacity - requires set rhosoil = 1.0
soil%rhosoil = 1.0 
! hcap ->  soil%css
CALL pack_landpts2mp( um1%ntiles, um1%land_pts, mp, um1%tile_pts, um1%tile_index,            &
                      um1%L_tile_pts, hcap, soil%css )

!initializations require for cable_user%soil_thermal_fix and %Haverd2013 
DO k=1,ms
   soil%ssat_vec(:,k)      = real(soil%ssat(:)   ,r_2)    
   soil%sucs_vec(:,k)      = real(soil%sucs(:)   ,r_2)   
   soil%hyds_vec(:,k)      = real(soil%hyds(:)   ,r_2)  
   soil%swilt_vec(:,k)     = real(soil%swilt(:)  ,r_2)  
   soil%bch_vec(:,k)       = real(soil%bch(:)    ,r_2)
   soil%sfc_vec(:,k)       = real(soil%sfc(:)    ,r_2)
   soil%rhosoil_vec(:,k)   = real(soil%rhosoil(:),r_2)   
   soil%cnsd_vec(:,k)      = real(soil%cnsd(:)   ,r_2)
   soil%css_vec(:,k)       = real(soil%css(:)    ,r_2)
   soil%sand_vec(:,k)      = real(soil%sand(:)   ,r_2)
   soil%clay_vec(:,k)      = real(soil%clay(:)   ,r_2)
   soil%silt_vec(:,k)      = real(soil%silt(:)   ,r_2)
   soil%watr(:,k)          = 0.001_r_2
END DO

WHERE (soil%ssat_vec .LE. 0.0 .AND. soil%sfc_vec .GT. 0.0)
  soil%ssat_vec = soil%sfc_vec + 0.05
END WHERE

RETURN
END SUBROUTINE initialize_soil

END MODULE cbl_um_init_soil_mod
