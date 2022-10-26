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
! Purpose: Calculates surface albedo, including from snow covered surface
!
! Called from: cbm
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b (but was previously in cable_radiation)
!
!
! ==============================================================================

MODULE cbl_soilColour_albedo_module

  IMPLICIT NONE

  PUBLIC soilcol_albedo
  PRIVATE

CONTAINS

  SUBROUTINE soilcol_albedo(ssnow, soil)

    USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
         r_2, mp, nrb
    ! Arguments
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow      ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil       ! soil parameters

    ! Local Variables
    INTEGER   :: ib
    REAL(r_2), DIMENSION(mp)      :: inc
    REAL(r_2), DIMENSION(mp, nrb) :: albsod,          & ! soil albedo (direct)
         albsoi             ! soil albedo (indirect)

    ! Look-up tables for soil albedo
    ! saturated soil albedos for 20 color classes and 2 wavebands (1=vis, 2=nir)
    REAL(r_2), DIMENSION(20,nrb) ::                                  &
         albsat,                                                                  &
         albdry

    REAL(r_2), PARAMETER, DIMENSION(20*nrb) ::                                  &
         albsat1D = (/ 0.25,0.23,0.21,0.20,0.19,0.18,0.17,0.16,0.15,0.14,0.13,    &
         0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05,0.04 ,             &
         0.50,0.46,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,    &
         0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,              &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                    &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /),                    &
                                ! dry soil albedos for 20 color classes and 2 wavebands (1=vis, 2=nir)
         albdry1D = (/  0.36,0.34,0.32,0.31,0.30,0.29,0.28,0.27,0.26,0.25,0.24,   &
         0.23,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,             &
         0.61,0.57,0.53,0.51,0.49,0.48,0.45,0.43,0.41,0.39,0.37,   &
         0.35,0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.16,             &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                   &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /)

    albsat = RESHAPE( albsat1D, (/20, nrb/) )
    albdry = RESHAPE( albdry1D, (/20, nrb/) )

    DO ib = 1,2 ! Number of wavebands (vis, nir)
       inc = MAX(0.11-0.40*ssnow%wb(:,1), 0._r_2)
       albsod(:,ib) = MIN(albsat(INT(soil%soilcol),ib)+inc, albdry(INT(soil%soilcol),ib))
       albsoi(:,ib) = albsod(:,ib)
    END DO
    ssnow%albsoilsn = REAL(0.5*(albsod + albsoi))

  END SUBROUTINE soilcol_albedo

END MODULE cbl_soilColour_albedo_module
