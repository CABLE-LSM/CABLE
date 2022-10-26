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
! Purpose: Fills CABLE type 'air' with appropriate values calculating
!          temperature dependent physical constants
!
! Called from: cbm, define_canopy
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b
!
!
! ==============================================================================

MODULE cable_air_module

  ! local pointers to global constants defined in

USE cable_phys_constants_mod, ONLY : CTFRZ   => TFRZ
USE cable_phys_constants_mod, ONLY : CRMAIR  => RMAIR
USE cable_phys_constants_mod, ONLY : CRGAS   => RGAS
USE cable_phys_constants_mod, ONLY : CCAPP   => CAPP
USE cable_phys_constants_mod, ONLY : CHL     => HL
USE cable_phys_constants_mod, ONLY : CRMH2O  => RMH2O
USE cable_phys_constants_mod, ONLY : CTETENA => TETENA
USE cable_phys_constants_mod, ONLY : CTETENB => TETENB
USE cable_phys_constants_mod, ONLY : CTETENC => TETENC
USE cable_phys_constants_mod, ONLY : CTETENA_ICE => TETENA_ICE
USE cable_phys_constants_mod, ONLY : CTETENB_ICE => TETENB_ICE
USE cable_phys_constants_mod, ONLY : CTETENC_ICE => TETENC_ICE


  IMPLICIT NONE

  PUBLIC define_air
  PRIVATE


CONTAINS


  SUBROUTINE define_air(met,air)

    USE cable_def_types_mod,          ONLY : air_type, met_type, mp

    TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
    TYPE (met_type), INTENT(IN)    :: met ! meteorological variables

    ! local vatiables
    REAL, DIMENSION(mp)     :: es ! sat vapour pressure (mb)

    ! END header

    ! Calculate saturation vapour pressure
    es = CTETENA * EXP( CTETENB * ( met%tvair - CTFRZ )                     &
         / ( CTETENC + ( met%tvair - CTFRZ ) ) )

    ! Calculate conversion factor from from m/s to mol/m2/s
    air%cmolar = met%pmb * 100.0 / (CRGAS * (met%tvair))

    ! Calculate dry air density:
    air%rho = MIN(1.3,CRMAIR * air%cmolar)

    ! molar volume (m^3/mol)
    air%volm = CRGAS * (met%tvair) / (100.0 * met%pmb)

    ! latent heat for water (j/kg)
    air%rlam= CHL

    ! saturation specific humidity
    air%qsat = (CRMH2O / CRMAIR) * es / met%pmb

    ! d(qsat)/dT ((kg/kg)/K)
    air%epsi = (air%rlam / CCAPP) * (CRMH2O / CRMAIR) * es * CTETENB *     &
         CTETENC / ( CTETENC + (met%tvair - CTFRZ) ) ** 2 / met%pmb

    ! air kinematic viscosity (m^2/s)
    air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met%tvair - CTFRZ) )

    ! psychrometric constant
    air%psyc = met%pmb * 100.0 * CCAPP * CRMAIR / air%rlam / CRMH2O

    ! d(es)/dT (mb/K)
    air%dsatdk = 100.0*(CTETENA*CTETENB*CTETENC)/((met%tvair-CTFRZ) +      &
         CTETENC)**2 * EXP( CTETENB * ( met%tvair-CTFRZ ) /         &
         ( (met%tvair-CTFRZ) + CTETENC) )

  END SUBROUTINE define_air



END MODULE cable_air_module
