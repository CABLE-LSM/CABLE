SUBROUTINE crop_driver( &
      doy, &
      climate, &
      ssnow, &
      soil, &
      veg, &
      canopy, &
      casaflux, &
      casamet, &
      casapool, &
      crop &
      )

  USE crop_def, only: &
      crop_type, &
      nc, &
      baresoil, &
      planted, &
      emergent, &
      growing, &
      Rgcoeff, &
      DMtoC
  USE crop_module
  USE casavariable, only: casa_flux, casa_met, casa_pool
  USE cable_def_types_mod, only: &
      climate_type, &
      soil_snow_type, &
      soil_parameter_type, &
      veg_parameter_type, &
      canopy_type, &
      dp=>r_2

  IMPLICIT NONE

  INTEGER, INTENT (in) :: doy
  TYPE (climate_type), INTENT (in) :: climate
  TYPE (soil_snow_type), INTENT (in) :: ssnow
  TYPE (soil_parameter_type), INTENT (in) :: soil
  TYPE (veg_parameter_type), INTENT (inout) :: veg
  TYPE (canopy_type), INTENT (inout) :: canopy
  TYPE (casa_flux), INTENT (inout) :: casaflux
  TYPE (casa_met), INTENT (inout) :: casamet
  TYPE (casa_pool), INTENT (inout) :: casapool
  TYPE (crop_type), INTENT (inout) :: crop

  ! local
  INTEGER :: ic,sl ! loop counters: crop type, soil layer
  REAL (dp) :: fPHU_day
  REAL (dp) :: SLA_C

  DO ic=1, nc

    WRITE (70,*) 'DOY: ', doy

    IF (crop%state(ic)==baresoil) THEN
      IF (doy>=crop%sowing_doymin(ic) .and. doy<=crop%sowing_doymax(ic)) THEN
         CALL planting(ic, doy, climate, ssnow, soil, crop)
      END IF

    ELSE IF (crop%state(ic)==planted) THEN
      CALL germination(ic, doy, climate, ssnow, soil, crop)
      CALL irrigation(ic, ssnow, soil, canopy)

      ! shouldn't be needed here! Check initialisation
      casapool%Cplant(ic,:) = 0.0_dp
      casamet%glai(ic) = 0.0_dp

    ELSE IF (crop%state(ic)==emergent .or. crop%state(ic)==growing) THEN

      ! update phenological heat units (start at 0 at germination!)
      fPHU_day = heat_units(climate%dtemp(ic), crop%Tbase(ic), crop%Tmax(ic)) &
          /crop%PHU_maturity(ic)
      crop%fPHU(ic) = crop%fPHU(ic) + fPHU_day

      ! calculate specific leaf area in m2 g-1 C
      SLA_C = SLA_development( &
          crop%fPHU(ic), &
          crop%sla_maturity(ic), &
          crop%sla_beta(ic) &
          )/DMtoC

      IF (crop%vernalisation(ic) .and. .not. crop%vacc(ic)) THEN
        CALL vernalisation(ic, climate, crop)
      END IF

      ! calculate C allocation factors
      CALL C_allocation_crops(ic, casaflux, crop)

      ! calculate growth respiration
      casaflux%Crgplant(ic) = &
          sum(Rgcoeff*casaflux%fracCalloc(ic,:))*casaflux%Cgpp(ic)

      WRITE (66,*) 'sum(Rgcoeff * casaflux%fracCalloc(ic,:))', &
          sum(Rgcoeff * casaflux%fracCalloc(ic,:))
      WRITE (66,*) 'casaflux%fracCalloc(ic,:)', casaflux%fracCalloc(ic,:)

      ! calculate maintenance respiration
      ! as in casa at the moment. Evtl calculate Rm first, then
      ! Rg. Discuss if it makes sense to replace Rgcoeff with Ygrowth as
      ! calculated in casa_cnp only need to calculate Rm and Rg here if we
      ! want to do it differently than in casa_rplant in that case, we might
      ! also change it in casa_rplant directly!

      ! calculate NPP
      casaflux%Cnpp(ic) = casaflux%Cgpp(ic) - &
          sum(casaflux%Crmplant(ic,:)) - casaflux%Crgplant(ic)

      IF (crop%state(ic)==emergent) THEN
        CALL emergence( &
            ic, &
            doy, &
            SLA_C, &
            fPHU_day, &
            veg, &
            casaflux, &
            casapool, &
            casamet, &
            crop &
            )
      ELSE IF (crop%state(ic)==growing) THEN

        WRITE (60,*) 'doy:', doy
        WRITE (60,*) '  casaflux%Cgpp(ic):', casaflux%Cgpp(ic)
        WRITE (60,*) '  sum(casaflux%Crmplant(ic,:)): ', &
            sum(casaflux%Crmplant(ic,:))
        WRITE (60,*) '  casaflux%Crgplant(ic): ', casaflux%Crgplant(ic)
        WRITE (60,*) '  casaflux%Cnpp_first(ic):', casaflux%Cnpp(ic)
        WRITE (70,*) 'climate%dtemp: ', climate%dtemp
        WRITE (70,*) 'fPHU_day: ', fPHU_day
        WRITE (70,*) 'crop%fPHU: ', crop%fPHU

        CALL senescence(ic, fPHU_day, casamet, casapool, casaflux, crop)
        CALL growth(ic, SLA_C, veg, casaflux, casapool, casamet, crop)

        WRITE (70,*) 'casaflux%Cgpp: ', casaflux%Cgpp
        WRITE (70,*) 'casaflux%Cnpp: ', casaflux%Cnpp
        WRITE (70,*) 'casamet%glai: ',  casamet%glai

        ! harvest if enough PHU accumulated
        IF (crop%fPHU(ic)>=1.0_dp) THEN
          CALL harvest(ic, doy, casapool, casamet, veg, crop)

          ! reset heat units and vernalisation requirements etc.
          crop%fPHU(ic) = 0.0_dp
          crop%fsenesc(ic) = 0.0_dp
          crop%VU(ic) = 0.0_dp
          crop%fVU(ic) = 0.0_dp
          crop%vacc(ic) = .FALSE.

          ! reset management settings
          canopy%irrig_surface(ic) = 0.0
          canopy%irrig_sprinkler(ic) = 0.0

          ! reset other variables
          crop%state_nrdays(ic,:) = 0
          crop%sl(ic) = 0

        END IF
      END IF ! crop%state == growing

      CALL irrigation(ic, ssnow, soil, canopy)

    END IF
  END DO

END SUBROUTINE crop_driver
