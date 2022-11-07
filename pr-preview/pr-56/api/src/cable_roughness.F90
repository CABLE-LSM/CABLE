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
! Purpose: Calculate roughness lengths as a function of soil and canopy 
!          parameters
!
! Contact: Eva.Kowalczyk@csiro.au
!
! History: No significant changes since v1.4b except change to cope with 
!          split timestep in ACCESS (zref_uv, zref_tq)
!
!
! ==============================================================================

MODULE cable_roughness_module

!! Calculate roughness lengths as a function of soil and canopyparameters 
USE cable_phys_constants_mod, ONLY : CCSD   => CSD 
USE cable_phys_constants_mod, ONLY : CCRD   => CRD 
USE cable_phys_constants_mod, ONLY : CCCD   => CCD 
USE cable_phys_constants_mod, ONLY : CCCW_C => CCW_C 
USE cable_phys_constants_mod, ONLY : CUSUHM  => USUHM 
USE cable_phys_constants_mod, ONLY : CVONK   => VONK
USE cable_phys_constants_mod, ONLY : CA33    => A33 
USE cable_phys_constants_mod, ONLY : CCTL    =>  CTL  
USE cable_phys_constants_mod, ONLY : CZDLIN  => ZDLIN
USE cable_phys_constants_mod, ONLY : CCSW    => CSW
USE cable_phys_constants_mod, ONLY : CGRAV   => GRAV
USE cable_other_constants_mod, ONLY : CLAI_THRESH => LAI_THRESH 
   
IMPLICIT NONE

real, parameter :: z0soilsn_min = 1.e-7
real, parameter :: z0soilsn_min_PF = 1.e-4
 
PRIVATE
PUBLIC ruff_resist

CONTAINS

SUBROUTINE ruff_resist(veg, rough, ssnow, canopy, LAI_pft, HGT_pft, reducedLAIdue2snow )

! see: Raupach, 1992, BLM 60 375-395
!      MRR notes "Simplified wind model for canopy", 23-oct-92
!      MRR draft paper "Simplified expressions...", dec-92
! modified to include resistance calculations by Ray leuning 19 Jun 1998  

USE cable_common_module, ONLY : cable_user, cable_runtime
USE cable_def_types_mod, ONLY : veg_parameter_type, roughness_type,         &
                                soil_snow_type, canopy_type, mp  
!subrs
USE cbl_hruff_mod, ONLY : HgtAboveSnow
USE cbl_LAI_eff_mod, ONLY : LAI_eff
!data
USE cable_other_constants_mod, ONLY : z0surf_min

implicit none

!result returned from called subr. Avail. in cross_*paths_module - but unsure
!yet
real :: HeightAboveSnow(mp) 
real :: reducedLAIdue2snow(mp) 

TYPE(roughness_type), INTENT(INOUT) :: rough
TYPE (canopy_type),   INTENT(INOUT) :: canopy
TYPE(soil_snow_type), INTENT(IN)    :: ssnow
TYPE (veg_parameter_type),  INTENT(INOUT) :: veg

real :: LAI_pft(mp)
real :: HGT_pft(mp)


REAL, DIMENSION(mp) ::                                                      &
  xx,      & ! =CCCD*LAI; working variable 
  dh         ! d/h where d is zero-plane displacement
integer :: i
   
! Set canopy height above snow level:
call HgtAboveSnow( HeightAboveSnow, mp, z0soilsn_min, veg%hc, ssnow%snowd, &
                   ssnow%ssdnn )
rough%hruff =  HeightAboveSnow

! LAI decreases due to snow: formerly canopy%vlaiw
call LAI_eff( mp, veg%vlai, veg%hc, HeightAboveSnow, &
              reducedLAIdue2snow )
   
canopy%vlaiw  = reducedLAIdue2snow
canopy%rghlai = canopy%vlaiw

IF (cable_user%soil_struc=='default') THEN

  ! Roughness length of bare soil (m): new formulation- E.Kowalczyk 2014
  IF (.NOT.cable_user%l_new_roughness_soil .AND. (.NOT.cable_user%or_evap)) THEN
    rough%z0soil = 0.0009*min(1.0,canopy%vlaiw) + 1.e-4
    rough%z0soilsn = rough%z0soil 
  ELSE
    rough%z0soil = 0.01*min(1.0,canopy%vlaiw) + 0.02*min(canopy%us**2/CGRAV,1.0)
    rough%z0soilsn = max(1.e-7,rough%z0soil)
  ENDIF

  WHERE( ssnow%snowd .GT. 0.01   )  &
    rough%z0soilsn =  MAX(z0soilsn_min, &
    rough%z0soil - rough%z0soil*MIN(ssnow%snowd,10.)/10.)
  WHERE( ssnow%snowd .GT. 0.01 .AND. veg%iveg == 17  )  &
    rough%z0soilsn =  MAX(rough%z0soilsn, z0soilsn_min_PF )

ELSEIF (cable_user%soil_struc=='sli') THEN

  rough%z0soil = 0.01*MIN(1.0,canopy%vlaiw) + 0.02*MIN(canopy%us**2/CGRAV,1.0)
  rough%z0soilsn = MAX(1.e-2,rough%z0soil) ! (1e-2: Mori et al., J Ag Met, 2010)

  WHERE( ssnow%snowd .GT. 0.01   )  &
    rough%z0soilsn =  MAX( 1.e-2,                                            &      
                          rough%z0soil - rough%z0soil*MIN(ssnow%snowd,10.)/10.)

ENDIF

do i=1,mp
  if( canopy%vlaiw(i) .LE. CLAI_THRESH  .OR.                                          &
      rough%hruff(i) .LT. rough%z0soilsn(i) ) then ! BARE SOIL SURFACE

    rough%z0m(i) = rough%z0soilsn(i)
    rough%rt0us(i) = 0.0
    rough%disp(i) = 0.0

    ! Reference height zref is height above the displacement height
    IF( cable_runtime%esm15 ) THEN
      rough%zref_uv(i) = MAX( 3.5, rough%za_uv(i) )
      rough%zref_tq(i) = MAX( 3.5, rough%za_tq(i) )
    ELSE     
      ! Ticket #148: Reference height is height above the displacement height
      ! noting that this must be above the roughness length and rough%hruff-rough%disp
      ! (though second case is unlikely to be attained)
      rough%zref_uv(i) = MAX( 3.5 + rough%z0m(i), rough%za_uv(i) )
      rough%zref_tq(i) = MAX( 3.5 + rough%z0m(i), rough%za_tq(i) )
      rough%zref_uv(i) = MAX( rough%zref_uv(i), rough%hruff(i)-rough%disp(i) )
      rough%zref_tq(i) = MAX( rough%zref_tq(i), rough%hruff(i)-rough%disp(i) )
    END IF
    
    rough%zruffs(i) = 0.0
    rough%rt1usa(i) = 0.0
    rough%rt1usb(i) = 0.0
   
    ! Friction velocity/windspeed at canopy height
    ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
    ! (CUSUHM set in physical_constants module):
    rough%usuh(i) = MIN( SQRT( CCSD + CCRD * ( canopy%vlaiw(i) * 0.5 ) ), CUSUHM )

    xx(i) = SQRT( CCCD * MAX( ( canopy%vlaiw(i) * 0.5 ), 0.0005 ) )

    ! Displacement height/canopy height:
    ! eq.8 Raupach 1994, BLM, vol 71, p211-216
    dh(i) = 1.0 - ( 1.0 - EXP( -xx(i) ) ) / xx(i)

    ! Extinction coefficient for wind profile in canopy:
    ! eq. 3.14, SCAM manual (CSIRO tech report 132)
    rough%coexp(i) = rough%usuh(i) / ( CVONK * CCCW_C * ( 1.0 - dh(i) ) )

  ELSE ! VEGETATED SURFACE

    ! Friction velocity/windspeed at canopy height
    ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
    ! (CUSUHM set in physical_constants module):
    rough%usuh(i) = MIN( SQRT( CCSD + CCRD * ( canopy%rghlai(i) * 0.5 ) ),       &
         CUSUHM )

    xx(i) = SQRT( CCCD * MAX( ( canopy%rghlai(i) * 0.5 ), 0.0005 ) )

    ! eq.8 Raupach 1994, BLM, vol 71, p211-216:
    dh(i) = 1.0 - ( 1.0 - EXP( -xx(i) ) ) / xx(i)

    ! Calculate zero-plane displacement:
    rough%disp(i) = dh(i)* rough%hruff(i)

    ! Calculate roughness length:
    rough%z0m(i) = ( (1.0 - dh(i)) * EXP( LOG( CCCW_C ) - 1. + 1. / CCCW_C       &
                    - CVONK / rough%usuh(i) ) ) * rough%hruff(i)

    ! Reference height zref is height above the displacement height
    IF( cable_runtime%esm15 ) THEN
      rough%zref_uv(i) = MAX( 3.5, rough%za_uv(i) )
      rough%zref_tq(i) = MAX( 3.5, rough%za_tq(i) )
    ELSE
     ! Reference height zref is height above the displacement height
     ! Ticket #148: Reference height is height above the displacement height
     ! noting that this must be above the roughness length and rough%hruff-rough%disp
     rough%zref_uv(i) = MAX( 3.5 + rough%z0m(i), rough%za_uv(i) )
     rough%zref_tq(i) = MAX( 3.5 + rough%z0m(i), rough%za_tq(i) )
     rough%zref_uv(i) = MAX( rough%zref_uv(i), rough%hruff(i)-rough%disp(i) )
     rough%zref_tq(i) = MAX( rough%zref_tq(i), rough%hruff(i)-rough%disp(i) )
    END IF

    ! find coexp: see notes "simplified wind model ..." eq 34a
    ! Extinction coefficient for wind profile in canopy:
    ! eq. 3.14, SCAM manual (CSIRO tech report 132)
    rough%coexp(i) = rough%usuh(i) / ( CVONK * CCCW_C * ( 1.0 - dh(i) ) )

    rough%term2(i)  = EXP( 2 * CCSW * canopy%rghlai(i) *                          &
         ( 1 - rough%disp(i) / rough%hruff(i) ) )
    rough%term3(i)  = CA33**2 * CCTL * 2 * CCSW * canopy%rghlai(i)
    rough%term5(i)  = MAX( ( 2. / 3. ) * rough%hruff(i) / rough%disp(i), 1.0 )
    rough%term6(i) =  EXP( 3. * rough%coexp(i) * ( rough%disp(i) / rough%hruff(i) -1. ) )
      rough%term6a(i) = EXP(rough%coexp(i) * ( 0.1 * rough%hruff(i) / rough%hruff(i) -1. ))

    ! eq. 3.54, SCAM manual (CSIRO tech report 132)
    rough%rt0us(i)  = rough%term5(i) * ( CZDLIN * LOG(                            &
         CZDLIN * rough%disp(i) / rough%z0soilsn(i) ) +                 &
         ( 1 - CZDLIN ) )                                         &
         * ( EXP( 2 * CCSW * canopy%rghlai(i) )  -  rough%term2(i) )    &
         / rough%term3(i)

    ! See CSIRO SCAM, Raupach et al 1997, eq. 3.49:
    rough%zruffs(i) = rough%disp(i) + rough%hruff(i) * CA33**2 * CCTL / CVONK /    &
         rough%term5(i)

    ! See CSIRO SCAM, Raupach et al 1997, eq. 3.51:
    rough%rt1usa(i) = rough%term5(i) * ( rough%term2(i) - 1.0 ) / rough%term3(i)
    rough%rt1usb(i) = rough%term5(i) * ( MIN( rough%zref_tq(i) + rough%disp(i),          &
         rough%zruffs(i) ) - rough%hruff(i) ) /                          &
         ( CA33**2 * CCTL * rough%hruff(i) )

    rough%rt1usb(i) = MAX( rough%rt1usb(i), 0.0 ) ! in case zrufs < rough%hruff

  END IF
END DO    

IF (cable_user%soil_struc.EQ.'sli') THEN
  WHERE( canopy%vlaiw .GE. CLAI_THRESH  .AND.                                          &
         rough%hruff .GE. rough%z0soilsn ) ! VEGETATED SURFACE

    ! vh ! Haverd et al., Biogeosciences 10, 2011-2040, 2013
    rough%rt0us  = LOG( rough%disp / (0.1 * rough%hruff) ) *                    &
                   EXP( 2. * CCSW * canopy%rghlai ) * rough%disp                &
                 / rough%hruff / (Ca33 ** 2 * Cctl) 

  ENDWHERE
ENDIF


END SUBROUTINE ruff_resist

END MODULE cable_roughness_module
