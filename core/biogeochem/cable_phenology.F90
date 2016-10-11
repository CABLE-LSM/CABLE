
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

 Use cable_def_types_mod, ONLY: met_type, climate_type, canopy_type, veg_parameter_type, &
      mp, r_2
 USE TypeDef,              ONLY: i4b, dp
 USE cable_IO_vars_module, ONLY: patch
 USE CABLE_COMMON_MODULE, ONLY: CurYear, filename, cable_user, HANDLE_ERR

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
   if (veg%iveg(np) == 31 .or. veg%iveg(np) == 2 .or. veg%iveg(np) == 5) then
      phen%doyphase(np,1) = -50
      phen%doyphase(np,2) = phen%doyphase(np,1) +14
      phen%doyphase(np,3) = 367
      phen%doyphase(np,4) = phen%doyphase(np,3) + 14
      phen%phase(np) = 2
   endif

   ! summergreen woody pfts
   if (veg%iveg(np) == 3 .or. veg%iveg(np) == 4) then  ! deciduous needleleaf(3) and broadleaf(4)

      ! Calculate GDD0  base value (=gdd to bud burst) for this PFT given
      !  current length of chilling period (Sykes et al 1996, Eqn 1)
      gdd0 = k_chilla + k_chillb*exp(-k_chillk*real(climate%chilldays(np)))
      phengdd5ramp = 200

      if (climate%gdd5(np).gt.gdd0 .and. phen%aphen(np).lt. APHEN_MAX) then

         phen_tmp = min(1.0_r_2, (climate%gdd5(np)-gdd0)/phengdd5ramp)

      else

         phen_tmp = 0.0

      endif

   endif

   ! summergreen grass or crops
   if (veg%iveg(np).ge.6.and.veg%iveg(np).le.10)  then     ! grass or crops

      phengdd5ramp = 50
      phen_tmp = min(1.0_r_2, climate%gdd5(np)/phengdd5ramp)

   endif

   ! raingreen pfts
   if (veg%iveg(np).ge.6.and.veg%iveg(np).le.10) then ! (grass or crops) need to include raingreen savanna trees here too

      if (climate%dmoist(np).lt. mmoisture_min) phen_tmp = 0.0


   endif

 if ((veg%iveg(np) == 3 .or. veg%iveg(np) == 4) .or. &
      (veg%iveg(np).ge.6.and.veg%iveg(np).le.10)) then


    if (phen_tmp.gt.0.0 .and.( phen%phase(np).eq.3 .or. phen%phase(np).eq.0 )) then
       phen%phase(np) = 1 ! greenup
       phen%doyphase(np,1) = climate%doy
    elseif (phen_tmp.ge.1.0_r_2 .and. phen%phase(np).eq.1) then
       phen%phase(np) = 2 ! steady LAI
       phen%doyphase(np,2) = climate%doy
    elseif (phen_tmp.lt.1.0_r_2 .and. phen%phase(np).eq.2) then
       phen%phase(np) = 3 ! senescence
       phen%doyphase(np,3) = climate%doy
    endif

    if (phen%phase(np)==3) then
       days =   min(climate%doy,365)-phen%doyphase(np,3)
       if (days < 0) days = days + 365
       if (days > 14) phen%phase(np) = 0          ! mimimum LAI
    endif

    ! Update annual leaf-on sum
    IF ((patch(np)%latitude>=0.0 .and. climate%doy==COLDEST_DAY_NHEMISPHERE).OR. &
         (patch(np)%latitude <0.0 .and. climate%doy==COLDEST_DAY_SHEMISPHERE) ) &
         phen%aphen(np) = 0
    phen%phen(np) = phen_tmp
    phen%aphen(np) = phen%aphen(np) + phen%phen(np)

 endif

ENDDO  ! end loop over patches


END SUBROUTINE cable_phenology_clim


! ==============================================================================
END MODULE cable_phenology_module
