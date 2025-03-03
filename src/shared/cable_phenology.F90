
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
! Purpose: climate-dependent phenology
!
! Called from: SUBROUTINE bgcdriver in casa_cable.F90
!
! History: Vanessa Haverd Jan 2015

! ==============================================================================
MODULE cable_phenology_module

  USE cable_def_types_mod, ONLY: met_type, climate_type, canopy_type, veg_parameter_type, &
       mp, r_2
  USE TypeDef,              ONLY: i4b, dp
  USE cable_IO_vars_module, ONLY: patch
  USE casa_ncdf_module, ONLY: HANDLE_ERR
  USE CABLE_COMMON_MODULE, ONLY: CurYear, filename, cable_user

CONTAINS
  ! ==============================================================================

  SUBROUTINE cable_phenology_clim (veg, climate, phen)

    ! sets the following days of year for use in allocation and leaf senescence
    ! algorithm depends on pft and climate
    !phen%doyphase(np,1) ! DOY for greenup
    !phen%doyphase(np,2) ! DOY for steady LAI
    ! phen%doyphase(np,3) ! DOY for leaf senescence
    !phen%doyphase(np,4) ! DOY for minimal LAI season


    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    IMPLICIT NONE

    TYPE (veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
    TYPE (phen_variable),      INTENT(INOUT) :: phen
    TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables
    INTEGER :: np, days
    REAL:: gdd0
    REAL(r_2) :: phen_tmp
    REAL, PARAMETER :: k_chilla = 0, k_chillb = 100, k_chillk = 0.05
    REAL, PARAMETER :: APHEN_MAX = 200.0, mmoisture_min=0.30
    INTEGER, PARAMETER:: COLDEST_DAY_NHEMISPHERE = 355
    INTEGER, PARAMETER:: COLDEST_DAY_SHEMISPHERE = 172
    REAL :: phengdd5ramp

    DO np= 1,mp

       ! evergreen pfts
       IF (veg%iveg(np) == 31 .OR. veg%iveg(np) == 2 .OR. veg%iveg(np) == 5) THEN
          phen%doyphase(np,1) = -50
          phen%doyphase(np,2) = phen%doyphase(np,1) +14
          phen%doyphase(np,3) = 367
          phen%doyphase(np,4) = phen%doyphase(np,3) + 14
          phen%phase(np) = 2
       ENDIF

       ! summergreen woody pfts
       IF (veg%iveg(np) == 3 .OR. veg%iveg(np) == 4) THEN  ! deciduous needleleaf(3) and broadleaf(4)

          ! Calculate GDD0  base value (=gdd to bud burst) for this PFT given
          !  current length of chilling period (Sykes et al 1996, Eqn 1)
          gdd0 = k_chilla + k_chillb*EXP(-k_chillk*REAL(climate%chilldays(np)))
          phengdd5ramp = 200

          IF (climate%gdd5(np).GT.gdd0 .AND. phen%aphen(np).LT. APHEN_MAX) THEN

             phen_tmp = MIN(1.0_r_2, (climate%gdd5(np)-gdd0)/phengdd5ramp)

          ELSE

             phen_tmp = 0.0

          ENDIF

       ENDIF

       ! summergreen grass or crops
       IF (veg%iveg(np).GE.6.AND.veg%iveg(np).LE.10)  THEN     ! grass or crops

          phengdd5ramp = 50
          phen_tmp = MIN(1.0_r_2, climate%gdd5(np)/phengdd5ramp)

       ENDIF

       ! raingreen pfts
       IF (veg%iveg(np).GE.6.AND.veg%iveg(np).LE.10) THEN ! (grass or crops) need to include raingreen savanna trees here too

          IF (climate%dmoist(np).LT. mmoisture_min) phen_tmp = 0.0


       ENDIF

       IF ((veg%iveg(np) == 3 .OR. veg%iveg(np) == 4) .OR. &
            (veg%iveg(np).GE.6.AND.veg%iveg(np).LE.10)) THEN


          IF (phen_tmp.GT.0.0 .AND.( phen%phase(np).EQ.3 .OR. phen%phase(np).EQ.0 )) THEN
             phen%phase(np) = 1 ! greenup
             phen%doyphase(np,1) = climate%doy
          ELSEIF (phen_tmp.GE.1.0_r_2 .AND. phen%phase(np).EQ.1) THEN
             phen%phase(np) = 2 ! steady LAI
             phen%doyphase(np,2) = climate%doy
          ELSEIF (phen_tmp.LT.1.0_r_2 .AND. phen%phase(np).EQ.2) THEN
             phen%phase(np) = 3 ! senescence
             phen%doyphase(np,3) = climate%doy
          ENDIF

          IF (phen%phase(np)==3) THEN
             days =   MIN(climate%doy,365)-phen%doyphase(np,3)
             IF (days < 0) days = days + 365
             IF (days > 14) phen%phase(np) = 0          ! mimimum LAI
          ENDIF

          ! Update annual leaf-on sum
          IF ((patch(np)%latitude>=0.0 .AND. climate%doy==COLDEST_DAY_NHEMISPHERE).OR. &
               (patch(np)%latitude <0.0 .AND. climate%doy==COLDEST_DAY_SHEMISPHERE) ) &
               phen%aphen(np) = 0
          phen%phen(np) = phen_tmp
          phen%aphen(np) = phen%aphen(np) + phen%phen(np)

       ENDIF

    ENDDO  ! end loop over patches


  END SUBROUTINE cable_phenology_clim


  ! ==============================================================================
END MODULE cable_phenology_module
