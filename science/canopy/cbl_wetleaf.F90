MODULE cable_wetleaf_module

  IMPLICIT NONE

  PUBLIC wetleaf
  PRIVATE

CONTAINS
SUBROUTINE wetLeaf( dels, cansat, tlfy,     &
                    gbhu, gbhf, ghwet, &
                    mp, CLAI_thresh, CCAPP, CRmair, & 
                    reducedLAIdue2snow, sum_rad_rniso, sum_rad_gradis,  & 
                    canopy_fevw,canopy_fevw_pot, canopy_fhvw, &
                    canopy_fwet, canopy_cansto, air_rlam, air_dsatdk, &
                    met_tvair, met_tk, met_dva, air_psyc )

USE cable_def_types_mod, ONLY : r_2
   
INTEGER, INTENT(IN) :: mp
REAL, INTENT(IN) :: CLAI_thresh, CCAPP, CRmair
REAL :: reducedLAIdue2snow(mp)
REAL, INTENT(INOUT) :: canopy_fevw(mp)
REAL, INTENT(INOUT) :: canopy_fevw_pot(mp)
REAL, INTENT(INOUT) :: canopy_fhvw(mp)
REAL, INTENT(INOUT) :: canopy_fwet(mp)

REAL, INTENT(IN) :: canopy_cansto(mp)
REAL, INTENT(IN) :: air_rlam(mp)
REAL, INTENT(IN) :: air_dsatdk(mp)
REAL, INTENT(IN) :: met_tvair(mp)
REAL, INTENT(IN) :: met_tk(mp)
REAL, INTENT(IN) :: met_dva(mp)
REAL, INTENT(IN) :: air_psyc(mp)
REAL, INTENT(IN) :: sum_rad_rniso(mp)
REAL, INTENT(IN) :: sum_rad_gradis(mp)

  REAL,INTENT(IN), DIMENSION(:) ::                                            &
       tlfy,          & ! leaf temp (K) - assCUMINg the temperature of
                            ! wet leaf is equal that of dry leaf ="tlfy"
       cansat           ! max canopy intercept. (mm)

  REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
       gbhu,          & ! forcedConvectionBndryLayerCond
       gbhf             ! freeConvectionBndryLayerCond

  REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
       ghwet            ! cond for heat for a wet canopy

  REAL, INTENT(IN)     :: dels ! integration time step (s)

  ! local variables
  REAL, DIMENSION(mp) ::                                                      &
       ccfevw,        & ! limitation term for
       gwwet,         & ! cond for water for a wet canopy
       ghrwet           ! wet canopy cond: heat & thermal rad

  !i sums, terms of convenience/readability
  REAL, DIMENSION(mp) ::                                                      &
     sum_gbh, xx1     ! xx1 not used - rk4417 - phase2

  INTEGER :: j
   
  ! END header

  ghwet = 1.0e-3
  gwwet = 1.0e-3
  ghrwet= 1.0e-3
  canopy_fevw = 0.0
  canopy_fhvw = 0.0
  sum_gbh = SUM((gbhu+gbhf),2)

  DO j=1,mp

      IF(reducedLAIdue2snow(j) > CLAI_THRESH) THEN

        ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
        ! calculate total thermal resistance, rthv in s/m
        ghwet(j) = 2.0   * sum_gbh(j)
        gwwet(j) = 1.075 * sum_gbh(j)
        ghrwet(j) = sum_rad_gradis(j) + ghwet(j)

        ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
        ! to avoid excessive evaporation:
         ccfevw(j) = MIN(canopy_cansto(j) * air_rlam(j) / dels, &
                       2.0 / (1440.0 / (dels/60.0)) * air_rlam(j) )

         canopy_fevw(j) = MIN( canopy_fwet(j) * ( air_dsatdk(j) *              &
                         ( sum_rad_rniso(j)- CCAPP*Crmair*( met_tvair(j)     &
                         - met_tk(j) ) * sum_rad_gradis(j) )                   &
                         + CCAPP * Crmair * met_dva(j) * ghrwet(j) )         &
                         / ( air_dsatdk(j)+air_psyc(j)*ghrwet(j) / gwwet(j) )  &
             , ccfevw(j) )

        !upwards flux density of water (kg/m2/s) - canopy componenet of 
        !potential evapotranspiration 
         canopy_fevw_pot(j) = ( air_dsatdk(j)* (sum_rad_rniso(j) -             &
                              CCAPP * Crmair * ( met_tvair(j) - met_tk(j) )  &
             *sum_rad_gradis(j) )                             &
                              + CCAPP * Crmair * met_dva(j) * ghrwet(j))     &
                              / (air_dsatdk(j)+air_psyc(j)*ghrwet(j)/gwwet(j) )

        ! calculate sens heat from wet canopy:
         canopy_fhvw(j) = canopy_fwet(j) * ( sum_rad_rniso(j) -CCAPP * Crmair&
                          * ( tlfy(j) - met_tk(j) ) * sum_rad_gradis(j) )      &
                           - canopy_fevw(j)

     ENDIF

  ENDDO

END SUBROUTINE wetLeaf

END MODULE cable_wetleaf_module
