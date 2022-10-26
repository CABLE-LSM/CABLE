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
! Purpose: Computes radiation absorbed by canopy and soil surface
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b
!
!
! ==============================================================================

MODULE cbl_radiation_module

  IMPLICIT NONE

  PUBLIC radiation
  PRIVATE

CONTAINS

SUBROUTINE radiation( ssnow, veg, air, met, rad, canopy, sunlit_veg_mask,&
  !constants
  clai_thresh, Csboltz, Cemsoil, Cemleaf, Ccapp &
)

    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
         veg_parameter_type, soil_snow_type,         &
         air_type, mp, mf, r_2

USE cable_other_constants_mod,  ONLY : Crad_thresh => rad_thresh
IMPLICIT NONE
logical :: sunlit_veg_mask(mp)
!constants
real :: CLAI_thresh
real :: CSboltz
real :: Cemsoil
real :: Cemleaf
real :: Ccapp

    TYPE (canopy_type),   INTENT(IN) :: canopy
    TYPE (air_type),      INTENT(IN) :: air
    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow
    TYPE (met_type),      INTENT(INOUT) :: met
    TYPE (radiation_type),INTENT(INOUT) :: rad

    TYPE (veg_parameter_type), INTENT(IN) :: veg

    REAL, DIMENSION(mp) ::                                                      &
         cf1, &      ! (1.0 - rad%transb * cexpkdm) / (extkb + extkdm(:,b))
         cf3, &      ! (1.0 - rad%transb * cexpkbm) / (extkb + extkbm(:,b))
         cf2n, &     ! exp(-extkn * vlai) (nitrogen)
         emair, &    ! air emissivity
         flpwb, &    ! black-body long-wave radiation
         flwv, &     ! vegetation long-wave radiation (isothermal)
         dummy, dummy2


    INTEGER :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave

    INTEGER, SAVE :: call_number =0

    call_number = call_number + 1

    ! Relative leaf nitrogen concentration within canopy:
    cf2n = EXP(-veg%extkn * canopy%vlaiw)

    rad%transd = 1.0

    WHERE (canopy%vlaiw > cLAI_thresh )    ! where vegetation exists....

       ! Diffuse SW transmission fraction ("black" leaves, extinction neglects
       ! leaf SW transmittance and REFLectance);
       ! from Monsi & Saeki 1953, quoted in eq. 18 of Sellers 1985:
       rad%transd = EXP(-rad%extkd * canopy%vlaiw)

    END WHERE

    ! Define fraction of SW beam tranmitted through canopy:
    !C!jhan: check rel. b/n extkb, extkbm,transb,cexpkbm def. cable_albedo, qsabbs
    !! vh_js !!
    dummy2 = MIN(rad%extkb * canopy%vlaiw,30.) ! vh version to avoid floating underflow !
    dummy = EXP(-dummy2)
    ! dummy2 = -rad%extkb * canopy%vlaiw
    ! dummy = EXP(dummy2)
    rad%transb = REAL(dummy)

    ! Define longwave from vegetation:
    flpwb = CSboltz * (met%tvrad) ** 4
    flwv = Cemleaf * flpwb

    rad%flws = CSboltz*Cemsoil* ssnow%tss **4

    ! Define air emissivity:
    emair = met%fld / flpwb

    rad%gradis = 0.0 ! initialise radiative conductance
    rad%qcan = 0.0   ! initialise radiation absorbed by canopy

    WHERE (canopy%vlaiw > CLAI_thresh )

       ! Define radiative conductance (Leuning et al, 1995), eq. D7:
       rad%gradis(:,1) = ( 4.0 * Cemleaf / (Ccapp * air%rho) ) * flpwb        &
            / (met%tvrad) * rad%extkd                              &
            * ( ( 1.0 - rad%transb * rad%transd ) /                &
            ( rad%extkb + rad%extkd )                              &
            + ( rad%transd - rad%transb ) /                        &
            ( rad%extkb - rad%extkd ) )

       rad%gradis(:,2) = ( 8.0 * Cemleaf / ( Ccapp * air%rho ) ) *            &
            flpwb / met%tvrad * rad%extkd *                        &
            ( 1.0 - rad%transd ) / rad%extkd - rad%gradis(:,1)

       ! Longwave radiation absorbed by sunlit canopy fraction:
       rad%qcan(:,1,3) = (rad%flws - flwv ) * rad%extkd *                       &
            ( rad%transd - rad%transb ) / ( rad%extkb - rad%extkd )&
            + ( emair- Cemleaf ) * rad%extkd * flpwb *            &
            ( 1.0 - rad%transd * rad%transb )                      &
            / ( rad%extkb + rad%extkd )

       ! Longwave radiation absorbed by shaded canopy fraction:
       rad%qcan(:,2,3) = ( 1.0 - rad%transd ) *                                 &
            ( rad%flws + met%fld - 2.0 * flwv ) - rad%qcan(:,1,3)

    END WHERE

    ! Convert radiative conductance from m/s to mol/m2/s:
    rad%gradis=SPREAD(air%cmolar, 2, mf)*rad%gradis
    rad%gradis = MAX(1.0e-3_r_2,rad%gradis)

    ! Update extinction coefficients and fractional transmittance for
    ! leaf transmittance and REFLection (ie. NOT black leaves):
    ! Define qcan for short wave (par, nir) for sunlit leaf:
    ! UM recieves met%fsd(:,b) forcing. assumed for offline that USED met%fsd(:,b) = 1/2* INPUT met%fsd
    DO b = 1, 2 ! 1 = visible, 2 = nir radiaition

       WHERE (sunlit_veg_mask) ! i.e. vegetation and sunlight are present

          cf1 = ( 1.0 - rad%transb * rad%cexpkdm(:,b) ) /                       &
               ( rad%extkb + rad%extkdm(:,b) )
          cf3 = (1.0 - rad%transb * rad%cexpkbm(:,b)) /                         &
               ( rad%extkb + rad%extkbm(:,b) )

          ! scale to real sunlit flux
          rad%qcan(:,1,b) = met%fsd(:,b) * (                                    &
               ( 1.0 - rad%fbeam(:,b) ) * ( 1.0 - rad%reffdf(:,b) )&
               * rad%extkdm(:,b) * cf1                             &
               + rad%fbeam(:,b) * ( 1.0-rad%reffbm(:,b) ) *        &
               rad%extkbm(:,b) * cf3                               &
               + rad%fbeam(:,b) * ( 1.0 - veg%taul(:,b)            &
               - veg%refl(:,b) ) * rad%extkb                       &
               * ( ( 1-rad%transb ) / rad%extkb - ( 1 -            &
               rad%transb**2 ) / ( rad%extkb + rad%extkb ) ) )

          ! Define qcan for short wave (par, nir) for shaded leaf:
          rad%qcan(:,2,b) = met%fsd(:,b) * (                                    &
               ( 1.0 - rad%fbeam(:,b) ) * ( 1.0 - rad%reffdf(:,b) )&
               * rad%extkdm(:,b) *                                 &
               ( ( 1.0 - rad%cexpkdm(:,b) ) / rad%extkdm(:,b)      &
               - cf1 ) + rad%fbeam(:,b) * ( 1. - rad%reffbm(:,b) ) &
               * rad%extkbm(:,b) * ( ( 1.0 - rad%cexpkbm(:,b) ) /  &
               rad%extkbm(:,b) - cf3 ) - rad%fbeam(:,b) *          &
               ( 1.0 - veg%taul(:,b) -veg%refl(:,b)) * rad%extkb   &
               * ( ( 1 - rad%transb ) / rad%extkb -                &
               ( 1 - rad%transb**2 ) / ( rad%extkb + rad%extkb ) ) )

       END WHERE

    END DO

    rad%qssabs = 0.

    WHERE (sunlit_veg_mask) ! i.e. vegetation and sunlight are present

       ! Calculate shortwave radiation absorbed by soil:
       ! (av. of transmitted NIR and PAR through canopy)*SWdown
       rad%qssabs = met%fsd(:,1) * (                                            &
            rad%fbeam(:,1) * ( 1. - rad%reffbm(:,1) ) *                 &
            EXP( -MIN(rad%extkbm(:,1) * canopy%vlaiw,20.) ) +           &
            ( 1. - rad%fbeam(:,1) ) * ( 1. - rad%reffdf(:,1) ) *        &
            EXP( -MIN(rad%extkdm(:,1) * canopy%vlaiw,20.) ) )           &
            + met%fsd(:,2) * ( rad%fbeam(:,2) * ( 1. - rad%reffbm(:,2) )&
            * rad%cexpkbm(:,2) + ( 1. - rad%fbeam(:,2) ) *              &
            ( 1. - rad%reffdf(:,2) ) * rad%cexpkdm(:,2) )

       ! Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C:
       rad%scalex(:,1) = ( 1.0 - rad%transb * cf2n ) / ( rad%extkb + veg%extkn )

       ! LAI of big leaf, sunlit, shaded, respectively:
       rad%fvlai(:,1) = ( 1.0 - rad%transb ) / rad%extkb
       rad%fvlai(:,2) = canopy%vlaiw - rad%fvlai(:,1)

    ELSEWHERE ! i.e. either vegetation or sunlight are NOT present

       ! Shortwave absorbed by soil/snow surface:
       rad%qssabs = ( 1.0 - ssnow%albsoilsn(:,1) ) * met%fsd(:,1) +             &
            ( 1.0 - ssnow%albsoilsn(:,2) ) * met%fsd(:,2)

       rad%scalex(:,1) = 0.0
       rad%fvlai(:,1) = 0.0
       rad%fvlai(:,2) = canopy%vlaiw

    END WHERE

    rad%scalex(:,2) = (1.0 - cf2n) / veg%extkn - rad%scalex(:,1)

    ! Total energy absorbed by canopy:
    rad%rniso = SUM(rad%qcan, 3)

  END SUBROUTINE radiation

END MODULE cbl_radiation_module
