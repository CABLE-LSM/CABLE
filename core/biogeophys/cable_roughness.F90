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

   USE cable_data_module, ONLY : irough_type, point2constants

   IMPLICIT NONE

   TYPE ( irough_type ) :: C
   PRIVATE
   PUBLIC ruff_resist

CONTAINS

SUBROUTINE ruff_resist(veg, rough, ssnow, canopy)

   ! m.r. raupach, 24-oct-92
   ! see: Raupach, 1992, BLM 60 375-395
   !      MRR notes "Simplified wind model for canopy", 23-oct-92
   !      MRR draft paper "Simplified expressions...", dec-92
   ! modified to include resistance calculations by Ray leuning 19 Jun 1998

   USE cable_common_module, ONLY : cable_user
   USE cable_def_types_mod, ONLY : veg_parameter_type, roughness_type,         &
                                   soil_snow_type, canopy_type, mp

   TYPE(roughness_type), INTENT(INOUT) :: rough
   TYPE (canopy_type),   INTENT(INOUT) :: canopy
   TYPE(soil_snow_type), INTENT(IN)    :: ssnow
   TYPE (veg_parameter_type),  INTENT(INOUT) :: veg

   REAL, DIMENSION(mp) ::                                                      &
      xx,      & ! =C%CCD*LAI; working variable
      dh         ! d/h where d is zero-plane displacement
   REAL :: take_off_disp

   CALL point2constants( C )

   !ASKJK - in NESP but not gm branch: remove ?
   ! JK: temporary bug fix: don't let veg%hc fall to 0
   veg%hc = MAX(veg%hc, 1.0)
   !write(90,*) 'veg%hc:', veg%hc
   !ASKJK - in NESP but not gm branch: remove ?

   ! for site-level runs, subtract displacement height from reference height
   take_off_disp = 0.0
   if (trim(cable_user%MetType) .EQ. 'site') then
     take_off_disp = 1.0
   endif


   ! Set canopy height above snow level:
   rough%hruff = MAX( 1.e-6, veg%hc - 1.2 * ssnow%snowd /                       &
        MAX( ssnow%ssdnn, 100. ) )

   ! LAI decreases due to snow:
   canopy%vlaiw = veg%vlai * rough%hruff / MAX( 0.01, veg%hc )
   canopy%rghlai = canopy%vlaiw

    IF (cable_user%soil_struc=='default') THEN

       ! Roughness length of bare soil (m): new formulation- E.Kowalczyk 2014
       IF (.not.cable_user%l_new_roughness_soil) THEN
          rough%z0soil = 0.0009*min(1.0,canopy%vlaiw) + 1.e-4
          rough%z0soilsn = rough%z0soil
       ELSE
          rough%z0soil = 0.01*min(1.0,canopy%vlaiw) + 0.02*min(canopy%us**2/C%GRAV,1.0)
          rough%z0soilsn = max(1.e-7,rough%z0soil)
       ENDIF

       WHERE( ssnow%snowd .GT. 0.01   )  &
            rough%z0soilsn =  max( 1.e-7, rough%z0soil - rough%z0soil*min(ssnow%snowd,10.)/10.)

    ELSEIF (cable_user%soil_struc=='sli') THEN

       rough%z0soil = 0.01*min(1.0,canopy%vlaiw) + 0.02*min(canopy%us**2/C%GRAV,1.0)
       rough%z0soilsn = max(1.e-2,rough%z0soil) ! (1e-2: Mori et al., J Ag Met, 2010)


       WHERE( ssnow%snowd .GT. 0.01   )  &
            rough%z0soilsn =  max( 1.e-2, rough%z0soil - rough%z0soil*min(ssnow%snowd,10.)/10.)

    ENDIF

   !! vh_js !! use LAI_THRESH here
   WHERE( canopy%vlaiw .LT. C%LAI_THRESH  .OR.                                          &
           rough%hruff .LT. rough%z0soilsn ) ! BARE SOIL SURFACE

      rough%z0m = rough%z0soilsn
      rough%rt0us = 0.0
      rough%disp = 0.0

      ! Reference height zref is height above the displacement height
      rough%zref_uv = MAX( 2.0, rough%za_uv )
      rough%zref_tq = MAX( 2.0, rough%za_tq )

      rough%zruffs = 0.0
      rough%rt1usa = 0.0
      rough%rt1usb = 0.0

      ! Friction velocity/windspeed at canopy height
      ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
      ! (C%USUHM set in physical_constants module):
      rough%usuh = MIN( SQRT( C%CSD + C%CRD * ( canopy%vlaiw * 0.5 ) ), C%USUHM )

      xx = SQRT( C%CCD * MAX( ( canopy%vlaiw * 0.5 ), 0.0005 ) )

      ! Displacement height/canopy height:
      ! eq.8 Raupach 1994, BLM, vol 71, p211-216
      dh = 1.0 - ( 1.0 - EXP( -xx ) ) / xx

      ! Extinction coefficient for wind profile in canopy:
      ! eq. 3.14, SCAM manual (CSIRO tech report 132)
      rough%coexp = rough%usuh / ( C%VONK * C%CCW_C * ( 1.0 - dh ) )

   ELSEWHERE ! VEGETATED SURFACE

      ! Friction velocity/windspeed at canopy height
      ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
      ! (C%USUHM set in physical_constants module):
      rough%usuh = MIN( SQRT( C%CSD + C%CRD * ( canopy%rghlai * 0.5 ) ),       &
                   C%USUHM )

      xx = SQRT( C%CCD * MAX( ( canopy%rghlai * 0.5 ), 0.0005 ) )

      ! eq.8 Raupach 1994, BLM, vol 71, p211-216:
      dh = 1.0 - ( 1.0 - EXP( -xx ) ) / xx

      ! Calculate zero-plane displacement:
      rough%disp = dh * rough%hruff

      ! Reference height zref is height above the displacement height
      rough%zref_uv = MAX(2.0, rough%za_uv - take_off_disp * rough%disp)
      rough%zref_tq = MAX(2.0, rough%za_tq - take_off_disp * rough%disp)
      rough%zref_uv = MAX(rough%zref_uv, rough%hruff - rough%disp)
      rough%zref_tq = MAX(rough%zref_tq, rough%hruff - rough%disp)
      

      ! Calculate roughness length:
      rough%z0m = ( (1.0 - dh) * EXP( LOG( C%CCW_C ) - 1. + 1. / C%CCW_C       &
                  - C%VONK / rough%usuh ) ) * rough%hruff

      ! find coexp: see notes "simplified wind model ..." eq 34a
      ! Extinction coefficient for wind profile in canopy:
      ! eq. 3.14, SCAM manual (CSIRO tech report 132)
      rough%coexp = rough%usuh / ( C%VONK * C%CCW_C * ( 1.0 - dh ) )

      rough%term2  = EXP( 2 * C%CSW * canopy%rghlai *                          &
                     ( 1 - rough%disp / rough%hruff ) )
      rough%term3  = C%A33**2 * C%CTL * 2 * C%CSW * canopy%rghlai
      rough%term5  = MAX( ( 2. / 3. ) * rough%hruff / rough%disp, 1.0 )
      rough%term6 =  EXP( 3. * rough%coexp * ( rough%disp / rough%hruff -1. ) )
      rough%term6a = EXP(rough%coexp * ( 0.1 * rough%hruff / rough%hruff -1. ))

     
         ! eq. 3.54, SCAM manual (CSIRO tech report 132)
         rough%rt0us  = rough%term5 * ( C%ZDLIN * LOG(                            &
              C%ZDLIN * rough%disp / rough%z0soilsn ) +                 &
              ( 1 - C%ZDLIN ) )                                         &
              * ( EXP( 2 * C%CSW * canopy%rghlai )  -  rough%term2 )    &
              / rough%term3

         ! See CSIRO SCAM, Raupach et al 1997, eq. 3.49:
         rough%zruffs = rough%disp + rough%hruff * C%A33**2 * C%CTL / C%VONK /    &
              rough%term5

         ! See CSIRO SCAM, Raupach et al 1997, eq. 3.51:
         rough%rt1usa = rough%term5 * ( rough%term2 - 1.0 ) / rough%term3
         rough%rt1usb = rough%term5 * ( MIN( rough%zref_tq + rough%disp,          &
              rough%zruffs ) - rough%hruff ) /                          &
              ( C%A33**2 * C%CTL * rough%hruff )

         rough%rt1usb = MAX( rough%rt1usb, 0.0 ) ! in case zrufs < rough%hruff

      END WHERE


       IF (cable_user%soil_struc.eq.'sli') THEN
         WHERE( canopy%vlaiw .GE. C%LAI_THRESH  .AND.                                          &
              rough%hruff .GE. rough%z0soilsn ) ! VEGETATED SURFACE

            rough%rt0us  = log(rough%disp/(0.1 * rough%hruff)) * &
                    EXP(2. * C%CSW * canopy%rghlai) * rough%disp &
                    / rough%hruff / (c%a33 ** 2 * c%ctl) ! vh ! Haverd et al., Biogeosciences 10, 2011-2040, 2013

         ENDWHERE
      ENDIF



END SUBROUTINE ruff_resist


END MODULE cable_roughness_module
