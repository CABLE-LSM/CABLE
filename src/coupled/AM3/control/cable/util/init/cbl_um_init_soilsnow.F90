MODULE cbl_um_init_soilsnow_mod
   
IMPLICIT NONE

CONTAINS

SUBROUTINE initialize_soilsnow( mp, msn, ms, TFRZ, land_pts, nsurft,           &
                                row_length, rows, ICE_SoilType,l_surft_pts,    &
                                surft_pts, surft_index, smvcst, tsoil_tile,    &
                                sthf_tile, smcl_tile, snow_tile, snow_rho1l,   &
                                snow_age, isnow_flg3l, snow_rho3l,             &
                                snow_depth3l, snow_mass3l, snow_tmp3l, soil,   &
                                ssnow, veg_iveg, met_tk )

! subrs
USE cable_init_wetfac_mod,    ONLY: initialize_wetfac

! data
USE cable_def_types_mod,      ONLY : r_2
USE cable_def_types_mod,      ONLY: soil_parameter_type, veg_parameter_type
USE cable_def_types_mod,      ONLY: met_type, soil_snow_type
USE cable_phys_constants_mod, ONLY: density_ice, density_liq

INTEGER, INTENT(IN) :: row_length        ! # columns in spatial grid
INTEGER, INTENT(IN) :: rows              ! # rows in spatial grid
INTEGER, INTENT(IN) :: land_pts          ! # land points being processed
INTEGER, INTENT(IN) :: nsurft            ! # tiles 
INTEGER, INTENT(IN) :: ms                ! # soil layers 
INTEGER, INTENT(IN) :: msn               ! # snow layers 
INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: ICE_SoilType      ! ice soil type index
INTEGER, INTENT(IN) :: surft_pts(nsurft) ! # land points per tile
INTEGER, INTENT(IN) :: surft_index(land_pts, nsurft) ! land_pt index of point
LOGICAL, INTENT(IN) :: l_surft_pts(land_pts, nsurft)
REAL,    INTENT(IN) :: smvcst(land_pts)
REAL,    INTENT(IN) :: sthf_tile(land_pts, nsurft, ms)
REAL,    INTENT(IN) :: smcl_tile(land_pts, nsurft, ms)
REAL,    INTENT(IN) :: tsoil_tile(land_pts, nsurft, ms)
REAL,    INTENT(IN) :: snow_rho1l(land_pts, nsurft)
REAL,    INTENT(IN) :: snow_tile(land_pts, nsurft)
REAL,    INTENT(IN) :: snow_age(land_pts, nsurft)
REAL,    INTENT(IN) :: snow_rho3l(land_pts, nsurft, msn)
REAL,    INTENT(IN) :: snow_depth3l(land_pts, nsurft, msn)
REAL,    INTENT(IN) :: snow_mass3l(land_pts, nsurft, msn)
REAL,    INTENT(IN) :: snow_tmp3l(land_pts, nsurft, msn)
INTEGER, INTENT(IN) :: isnow_flg3l(land_pts, nsurft)
REAL,    INTENT(IN) :: TFRZ 
INTEGER, INTENT(IN) :: veg_iveg(mp)
REAL,    INTENT(IN) :: met_tk(mp)

TYPE(soil_parameter_type), INTENT(IN)  :: soil       ! soil parameters
TYPE(soil_snow_type),      INTENT(OUT) :: ssnow      ! 

!local vars
INTEGER :: i,j,k,L,n
REAL    :: zsetot
REAL    :: ice_vol_tmp(land_pts, nsurft, ms)
REAL    :: wbtot_mp(mp, ms)
REAL    :: wbice_mp(mp, ms)

ssnow%pudsto       = 0.0 
ssnow%pudsmx       = 0.0
ssnow%wbtot        = 0.0
ssnow%totwblake    = 0.0  ! wb_lake integrated over river timestep
ssnow%tggav        = 0.0
ssnow%qhz          = 0.0
ssnow%qhlev        = 0.0
ssnow%qrecharge    = 0.0
ssnow%rtevap_sat   = 0.0
ssnow%rtevap_unsat = 0.0

! identify module parameters here and recast (NOT ssnow% state variables)
ssnow%t_snwlr      = 0.05
ssnow%rtsoil       = 50.0
ssnow%satfrac      = 0.5
ssnow%wtd          = 1.0

ssnow%snowd  = PACK( snow_tile, l_surft_pts )
ssnow%snage  = PACK( snow_age, l_surft_pts )
ssnow%ssdnn  = PACK( snow_rho1l, l_surft_pts )  
ssnow%isflag = PACK( isnow_flg3L, l_surft_pts )  
ssnow%ssdnn  = PACK( snow_rho1l, l_surft_pts )  
ssnow%isflag = PACK( isnow_flg3l, l_surft_pts )  

! over snow layers
DO j=1, msn
  ssnow%sdepth(:,j)= PACK( snow_depth3l(:,:,j), l_surft_pts )
  ssnow%smass(:,j) = PACK( snow_mass3l(:,:,j), l_surft_pts )  
  ssnow%ssdn(:,j)  = PACK( snow_rho3l(:,:,j), l_surft_pts )  
  ssnow%tggsn(:,j) = PACK( snow_tmp3l(:,:,j), l_surft_pts )  
ENDDO 
         
! over soil layers
DO j=1, ms
  ssnow%tgg(:,j)   = PACK( tsoil_tile(:,:,j), l_surft_pts )
ENDDO 

ssnow%osnowd    = ssnow%snowd

zsetot = sum(soil%zse)
DO k = 1, ms
   ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
END DO

ice_vol_tmp(:,:,:)  = 0.0
DO n=1,nsurft                                                       
  DO k=1,surft_pts(n)                                           
    i = surft_index(k,n)                                      
    DO j = 1, ms
      ice_vol_tmp(i,n,j)  = sthf_tile(i,n,j) * smvcst(i)
    ENDDO ! J
  ENDDO
ENDDO
 
DO j = 1, ms
  !liq volume  from (tot_mass - ice_mass) / (dz*rho_liq)
  wbtot_mp(:,j)    = PACK( smcl_tile(:,:,j), l_surft_pts )
  
  ssnow%wbice(:,J) = PACK( ice_vol_tmp(:,:,j), l_surft_pts )
  ssnow%wbice(:,J) = MAX ( 0.0 , ssnow%wbice(:,j) )
  wbice_mp(:,j)    = ssnow%wbice(:,j) * ( soil%zse(j) * density_ice )
  
  ssnow%wbliq(:,j) = ( wbtot_mp(:,j) - wbice_mp(:,j) )  / ( soil%zse(j) * density_liq )
  ssnow%wb(:,j)    = ssnow%wbice(:,j) + ssnow%wbliq(:,j) 
END DO 

! wetfac initialized here as used to init owetfac on firstimestep in cbm
! add to startdump? includes Temporay fix for accounting for reduction of 
! soil evaporation due to freezing. Also includes specific lakes case  
! Prevents divide by zero at glaciated points where both wb and wbice=0.
CALL initialize_wetfac( mp, ssnow%wetfac, soil%swilt, soil%sfc,            &
                        ssnow%wb(:,1), ssnow%wbice(:,1), ssnow%snowd,      &
                        veg_iveg, met_tk, tfrz ) 

! initialized here on first call 
ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1) 
 
!jhan: do we want to do this before %owetfac is set 
DO J = 1, ms
  !should be removed!!!!!!!! This cannot conserve if there are any
  !dynamics 
  WHERE( soil%isoilm == ICE_SoilType)
    ssnow%wb(:,j)    = 0.95 * soil%ssat
    ssnow%wbice(:,j) = 0.85 * ssnow%wb(:,j)
  ENDWHERE
  
  !no not force rho_water==rho_ice==1000.0
  ssnow%wbtot = ssnow%wbtot + soil%zse(j) *                       &
                             ( ssnow%wbliq(:,j) * density_liq +   &
                               ssnow%wbice(:,j) * density_ice )                     
ENDDO
     
RETURN
END SUBROUTINE initialize_soilsnow

END MODULE cbl_um_init_soilsnow_mod

