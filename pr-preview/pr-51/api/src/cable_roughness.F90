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

!*# Overview
!
! The procedures contained in this module calculate the roughness parameters
! and the aerodynamic contribution to the resistances controlling the fluxes 
! of momentum, heat and water vapour between the land and atmosphere for each 
! land point.
!
! The formulations take into account vegetation and snow cover. The dependence
! on atmospheric conditions (i.e. the surface heat fluxes) is incorporated
! later within the [[define_canopy]] subroutine.


IMPLICIT NONE

real, parameter :: z0soilsn_min = 1.e-7
  !! Minimum value for the roughness length for bare soil, \(z_{0,min}\) (m)
real, parameter :: z0soilsn_min_PF = 1.e-4
  !! Minimum value for the roughness length for permanent ice on land, 
  !! \(z_{0,minPF}\) (m)
 
PRIVATE
PUBLIC ruff_resist

CONTAINS

SUBROUTINE ruff_resist(veg, rough, ssnow, canopy, LAI_pft, HGT_pft, reducedLAIdue2snow )
!* Calculates the roughness parameters and the aerodynamic contribution to the
! resistances controlling the fluxes between the land and atmosphere for each 
! land point.
!
!## Scientific description
!
! The scientific basis for the formulae is a combination of Localised Near 
! Field theory for aerodynamic transfer within canopies (for the heat and 
! water vapour fluxes) and bulk formulae for the roughness length and 
! displacement height accounting for roughness sublayer effects (for the 
! momentum flux).
!
! The calculations take into account whether:
! 
! 1. the soil model used is SLI or the default (soilsnow)
! 2. the land point is vegetated (LAI > LAI_THRESH) or not
!
! The primary references for these formulations are
!
! * [Raupach MR (1989)](https://doi.org/10.1002/qj.49711548710)
! * [Raupach MR (1992)](https://doi.org/10.1007/BF00155203)
! * [Raupach MR (1994)](https://doi.org/10.1007/BF00709229)
! * [Kowalczyk et al. (2006)](http://www.cmar.csiro.au/e-print/open/kowalczykea_2006a.pdf)
!
!## Outputs
!
! The principal outputs from the MODULE are
!
! * `rough%zref_tq`: reference height above the displacement height for the 
!   air temperature and humidity, where these forcing observations are deemed
!   to have been collected (m) 
! * `rough%zref_uv`: height above the displacement height for the 
!   wind speeds, where the forcing observations are deemed to have been 
!   collected (m) 
! * `rough%hruff`: the canopy height accounting for the presence of snow (m)
! * `canopy%vlaiw`: the leaf area index accounting for the presence of snow 
!   (m\(^2\) m\(^{-2}\))
! * `rough%z0soil`: the aerodynamic roughness length for soil (m)
! * `rough%z0ssoilsn`: the aerodynamic roughness length for snow (m)
! * `rough%z0m`: the aerodynamic roughness length for the surface 
!   (canopy+soil+snow) (m)
! * `rough%disp`: the displacement height of the surface 
!   (=0.0 if not vegetated) (m)
! * `rough%zruffs`: the depth of the roughness sublayer over vegetated
!   surfaces (m)
! * `rough%coexp`: the extinction coefficient for the wind speed profile 
!   within the canopy, \(c_{0}\) (m\(^{-1}\)). 
!    
!    \(c_{0}\) is defined to be the coefficient within an expontial wind 
!    profile, i.e. the wind speed at height \(z\) above the ground within 
!    a canopy of height \(h_c\) is given by 
!    \( U(z) = U_{h} \exp\{ c_{0} (z-h_c) \} \) where \(U_h\) is the wind 
!    speed at canopy top.
!
! * `rough%usuh`: the ratio of the friction velocity, \(u_*\), to the wind 
!   speed at canopy top (-)
! * `rough%rt0us`: *normalised* aerodynamic resistance for the turbulent
!   transfer from the soil/snow surface to the displacement height (-)
! * `rough%rt1us`: *normalised* aerodynamic resistance for the turbulent 
!   transfer from the displacement height to the reference level (-)
!
!    `rough%rt1us` is evaluated in three subparts (`rough%rt1usa`, 
!    `rough%rt1usb`, and `rough%rt1usc`). One of the resistance terms 
!    (`rough%rt1usc`) is evaluated in subroutine [[define_canopy]].
!
! Each of the normalized resistances are given by the theoretical formulae 
! given by the references. The aerodynamic resistances for the current 
! time step are evaluated later by dividing the *normalized resistances* by 
! the current time step's friction velocity `canopy%us`. 



USE cable_common_module, ONLY : cable_user, cable_runtime
USE cable_def_types_mod, ONLY : veg_parameter_type, roughness_type,         &
                                soil_snow_type, canopy_type, mp  
!subrs
USE cbl_hruff_mod, ONLY : HgtAboveSnow
USE cbl_LAI_eff_mod, ONLY : LAI_eff
!data
USE cable_other_constants_mod, ONLY : z0surf_min

implicit none

!!## Structure

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
!* * evaluates the canopy height and leaf area given the presence of snow 
!    (or not) using [[HgtAboveSnow]] and [[LAI_eff]]
rough%hruff =  HeightAboveSnow

! LAI decreases due to snow: formerly canopy%vlaiw
call LAI_eff( mp, veg%vlai, veg%hc, HeightAboveSnow, &
              reducedLAIdue2snow )
   
canopy%vlaiw  = reducedLAIdue2snow
canopy%rghlai = canopy%vlaiw

!* * sets the value of soil and snow roughness lengths
!    (depends on the configuration of CABLE)
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

!| * evaluates the remaining output variables, depending on whether the land 
!    point is vegetated or not.
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

!> * if the SLI soil model is used - updates the evaluated `rough%rt0us`
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
