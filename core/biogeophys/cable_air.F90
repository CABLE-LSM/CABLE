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
   USE cable_data_module, ONLY : iair_type, point2constants
   
   IMPLICIT NONE

   PUBLIC define_air
   PRIVATE

   TYPE( iair_type ) :: C


CONTAINS


SUBROUTINE define_air(met,air)

   USE cable_def_types_mod,          ONLY : air_type, met_type, mp    
   USE cable_common_module,   ONLY : cable_runtime, cable_user, ktau_gl 

   TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
   TYPE (met_type), INTENT(IN)    :: met ! meteorological variables
   
   ! local vatiables 
   REAL, DIMENSION(mp)     :: es ! sat vapour pressure (mb)   
   INTEGER :: i 
      
   ! END header
   
   ! local ptrs to constants defined in cable_data_module
   call point2constants( C )
   
   ! Calculate saturation vapour pressure
   es = C%TETENA * EXP( C%TETENB * ( met%tvair - C%TFRZ )                     &
        / ( C%TETENC + ( met%tvair - C%TFRZ ) ) )
   
   ! Calculate conversion factor from from m/s to mol/m2/s
   air%cmolar = met%pmb * 100.0 / (C%RGAS * (met%tvair))
   
   ! Calculate dry air density:
   air%rho = MIN(1.3,C%RMAIR * air%cmolar)
   
   ! molar volume (m^3/mol)
   air%volm = C%RGAS * (met%tvair) / (100.0 * met%pmb)
   
   ! latent heat for water (j/kg)
   air%rlam= C%HL
   
   ! saturation specific humidity
   air%qsat = (C%RMH2O / C%RMAIR) * es / met%pmb
   
   ! d(qsat)/dT ((kg/kg)/K)
   air%epsi = (air%rlam / C%CAPP) * (C%RMH2O / C%RMAIR) * es * C%TETENB *     &
              C%TETENC / ( C%TETENC + (met%tvair - C%TFRZ) ) ** 2 / met%pmb
   
   ! air kinematic viscosity (m^2/s)
   air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met%tvair - C%TFRZ) )
   
   ! psychrometric constant
   air%psyc = met%pmb * 100.0 * C%CAPP * C%RMAIR / air%rlam / C%RMH2O
   
   ! d(es)/dT (mb/K)
   air%dsatdk = 100.0*(C%TETENA*C%TETENB*C%TETENC)/((met%tvair-C%TFRZ) +      &
                C%TETENC)**2 * EXP( C%TETENB * ( met%tvair-C%TFRZ ) /         &
                ( (met%tvair-C%TFRZ) + C%TETENC) )
  
END SUBROUTINE define_air



END MODULE cable_air_module
