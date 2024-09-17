MODULE cbl_um_init_soil_mod
   
IMPLICIT NONE
PUBLIC initialize_soil

CONTAINS
        
SUBROUTINE initialize_soil( nsurft, land_pts, ms, mp, nsoil_max, ICE_soiltype, &
                            surft_pts, surft_index, L_tile_pts, soiltype,        &
                            bexp, hcon, satcon, sathh, smvcst, smvcwt,         &
                            smvccl, albsoil, dzsoil, pars, soil )

! subrs
USE cable_pack_mod,      ONLY: pack_landpts2mp_ICE, pack_landpts2mp

! data
USE cable_def_types_mod, ONLY: soil_parameter_type, r_2
USE params_io_mod_cbl,   ONLY: params_io_data_type
USE cable_common_module, ONLY : cable_user, gw_params
IMPLICIT NONE
                            
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: ms 
INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: nsoil_max
REAL,    INTENT(IN) :: dzsoil(ms)        ! soil layer thicknesses 
INTEGER, INTENT(IN) :: ICE_soiltype
INTEGER, INTENT(IN) :: surft_pts(nsurft) ! # land points on each tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! index of tile points 
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)
INTEGER, INTENT(IN) :: SoilType(mp)    !CABLE soil type per tile
REAL,    INTENT(IN) :: bexp(land_pts)
REAL,    INTENT(IN) :: hcon(land_pts)
REAL,    INTENT(IN) :: satcon(land_pts)
REAL,    INTENT(IN) :: sathh(land_pts)
REAL,    INTENT(IN) :: smvcst(land_pts)
REAL,    INTENT(IN) :: smvcwt(land_pts)
REAL,    INTENT(IN) :: smvccl(land_pts)
REAL,    INTENT(IN) :: albsoil(land_pts) 
   
TYPE(soil_parameter_type), INTENT(INOUT) :: soil       ! soil parameters
TYPE(params_io_data_type), INTENT(IN)      :: pars

!__ local vars
INTEGER :: i,j,k
REAL    :: hcon_ICE(nsoil_max) 
REAL    :: soilCnsd_real(mp) 
REAL    :: sucs_sign_factor, hyds_unit_factor, sucs_min_magnitude
REAL, PARAMETER   :: ssat_lo = 0.15
REAL, PARAMETER   :: ssat_hi = 0.65
REAL, PARAMETER   :: rhob_lo = 810.0
REAL, PARAMETER   :: rhob_hi = 2300.0
REAL, ALLOCATABLE :: znode(:), ssat_bounded(:,:),rho_soil_bulk(:,:)

! defined in namelist in UM

soil%zse_vec   = spread(dzsoil,1,mp)
soil%watr(:,:) = 0.05
soil%GWwatr(:) = 0.05
  
! distance between consecutive layer midpoints
soil%zshh(1)    = 0.5 * soil%zse(1) 
soil%zshh(ms+1) = 0.5 * soil%zse(ms)
soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))

! albsoil ->  soil%albsoil
CALL pack_landpts2mp( nsurft, land_pts, mp, surft_pts, surft_index,            &
                      L_tile_pts, albsoil, soil%albsoil(:,1) )
  
! bexp -> soil%bch: parameter b in Campbell equation 
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, BEXP, SoilType,  &
                          pars%soilin_bch, soil%bch )

hcon_ICE(:) = ( pars%soilin_sand(ICE_SoilType) * 0.3 )                         &
             + ( pars%soilin_clay(ICE_SoilType) * 0.25 )                       &
             + ( pars%soilin_silt(ICE_SoilType) * 0.265 )

! hcon -> soil%cnsd
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, hcon, soiltype,  &
                          hcon_ICE, soilCnsd_real )

soil%cnsd = REAL( soilCnsd_real, r_2 )

! CABLE soil parameters are in general PACKed from UM spatial fields
! *EXCEPT* for permanent ice points where we overwrite with parametrs read from
! cable_soil_params (pars%)

! soil%hyds = hydraulic conductivity @saturation is PACKed from satcon[mm/s]
! *EXCEPT* at permanent ice points where we overwrite with parametrs read from
! pars%soilin_hyds[m/s] which are already in correct units. Dividing by 1000  
! soil%hyds is correct for tiles 1:8. Would be 1000* too small, hence *1000 first 
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, satcon, soiltype,&
                          pars%soilin_hyds * 1000.0, soil%hyds )
soil%hyds    = soil%hyds / 1000.0

! sathh -> soil%sucs
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, sathh, soiltype, &
                          pars%soilin_sucs, soil%sucs )

! smvcst -> soil%ssat
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, smvcst, soiltype,&
                          pars%soilin_ssat, soil%ssat )

! smvcwt -> soil%swilt
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, smvcwt, soiltype,&
                          pars%soilin_swilt, soil%swilt )

! smvccl -> soil%sfc  
CALL pack_landpts2mp_ICE( nsurft, land_pts, mp, nsoil_max, ICE_soiltype,       &
                          surft_pts, surft_index, L_tile_pts, smvccl, soiltype,&
                          pars%soilin_sfc, soil%sfc )

!--- (re)set values for CABLE
soil%ibp2    =  NINT(soil%bch) + 2
soil%i2bp3   =  2 * NINT(soil%bch) + 3

soil%ssat    = MAX( soil%ssat, soil%sfc + 0.01 )
soil%sucs    = ABS( soil%sucs )
WHERE( soil%ssat > 0. ) 
  soil%pwb_min =  (soil%swilt / soil%ssat )**soil%ibp2
END WHERE
sucs_min_magnitude = 106.0/1000.0
soil%sucs          = MAX(sucs_min_magnitude,soil%sucs)
soil%hsbh          = soil%hyds * ABS(soil%sucs) * soil%bch
       
!from CM2 - BUT i don't think even necessary here. CM2 has a LOT of GW
!stuff here and cable_user%soil_thermal_fix stuff
DO k=1,ms
   soil%ssat_vec(:,k)      = real(soil%ssat(:)   ,r_2)    
   soil%sucs_vec(:,k)      = real(soil%sucs(:)   ,r_2)   
   soil%hyds_vec(:,k)      = real(soil%hyds(:)   ,r_2)  
   soil%swilt_vec(:,k)     = real(soil%swilt(:)  ,r_2)  
   soil%bch_vec(:,k)       = real(soil%bch(:)    ,r_2)
   soil%sfc_vec(:,k)       = real(soil%sfc(:)    ,r_2)
   soil%rhosoil_vec(:,k)   = real(soil%rhosoil(:),r_2)   
   soil%cnsd_vec(:,k)      = real(soil%cnsd      ,r_2)
   soil%css_vec(:,k)       = real(soil%css       ,r_2)
   soil%watr(:,k)          = 0.001_r_2
END DO

WHERE (soil%ssat_vec .LE. 0.0 .AND. soil%sfc_vec .GT. 0.0)
  soil%ssat_vec = soil%sfc_vec + 0.05
END WHERE

! review:: good chance we want some of this is CM3
!jhan!IF (cable_user%soil_thermal_fix) THEN
!jhan!
!jhan!if (allocated(ssat_bounded)) deallocate(ssat_bounded)
!jhan!if (allocated(rho_soil_bulk)) deallocate(rho_soil_bulk)
!jhan!
!jhan!allocate(ssat_bounded(size(soil%ssat_vec,dim=1),&
!jhan!                      size(soil%ssat_vec,dim=2) ) )
!jhan!
!jhan!ssat_bounded(:,:) = min( ssat_hi, max(ssat_lo, &
!jhan!                                   soil%ssat_vec(:,:) ) )
!jhan!
!jhan!allocate(rho_soil_bulk(size(soil%rhosoil_vec,dim=1),&
!jhan!                       size(soil%rhosoil_vec,dim=2) ) )
!jhan!
!jhan!rho_soil_bulk(:,:) = min(rhob_hi, max(rhob_lo , &
!jhan!                       (2700.0*(1.0 - ssat_bounded(:,:)) ) ) )
!jhan!
!jhan!
!jhan!do k=1,ms
!jhan!   do i=1,mp
!jhan!
!jhan!
!jhan!      if (soil%isoilm(i) .ne. 9) then
!jhan!
!jhan!         soil%rhosoil_vec(i,k) = 2700.0
!jhan!
!jhan!         soil%cnsd_vec(i,k) = ( (0.135*(1.0-ssat_bounded(i,k))) +&
!jhan!                             (64.7/rho_soil_bulk(i,k)) ) / &
!jhan!                           (1.0 - 0.947*(1.0-ssat_bounded(i,k)))
!jhan!
!jhan!      end if
!jhan!
!jhan!   end do
!jhan!end do
!jhan!
!jhan!k=1
!jhan!do i=1,mp
!jhan!   if (soil%isoilm(i) .ne. 9) then
!jhan!      soil%rhosoil(i) = soil%rhosoil_vec(i,1)
!jhan!      soil%cnsd(i)    = soil%cnsd_vec(i,1)
!jhan!   end if
!jhan!end do
!jhan!
!jhan!  if (allocated(ssat_bounded)) deallocate(ssat_bounded)
!jhan!  if (allocated(rho_soil_bulk)) deallocate(rho_soil_bulk)
!jhan!
!jhan!END IF
!jhan!
!jhan!!node depths
!jhan!IF (allocated(znode)) deallocate(znode)
!jhan!allocate(znode(ms))
!jhan!
!jhan!znode(1) = soil%zshh(1)
!jhan!do k=2,ms
!jhan!   znode(k) = znode(k-1) * 0.5*(soil%zse(k-1)+soil%zse(k))
!jhan!end do
!jhan!
!jhan!IF (cable_user%gw_model) THEN
!jhan!  
!jhan!   DO k=1,ms
!jhan!
!jhan!      do i=1,mp  !from reversing pedotransfer functions
!jhan!                 !,ay cause io issues because not passed into um
!jhan!
!jhan!         if (soil%isoilm(i) .ne. 9) then
!jhan!
!jhan!               soil%hyds_vec(i,k) = soil%hyds_vec(i,k) * &   !change in hyds
!jhan!                                   exp(-gw_params%hkrz*( znode(k)-gw_params%zdepth) )
!jhan!
!jhan!         end if
!jhan!
!jhan!
!jhan!      end do
!jhan!   end do
!jhan!
!jhan!   k=1
!jhan!   soil%hyds(:) = soil%hyds_vec(:,k)
!jhan!
!jhan!END IF

RETURN
END SUBROUTINE initialize_soil

END MODULE cbl_um_init_soil_mod
