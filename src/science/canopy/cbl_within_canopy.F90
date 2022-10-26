MODULE cable_within_canopy_module
  
  IMPLICIT NONE

  PUBLIC within_canopy
  PRIVATE

CONTAINS

SUBROUTINE within_canopy( mp, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, CLAI_thresh, &
                           CCAPP, CTFRZ, rad,rough, air, met, veg, canopy, ssnow, gbhu, gbhf,    &
                           qstvair, rt0, rhlitt, relitt )
  
  USE cable_common_module, ONLY : cable_user

  USE cable_def_types_mod, ONLY : r_2

  USE cable_def_types_mod, ONLY : canopy_type, air_type, met_type,             &
                                  radiation_type, roughness_type,              &
                                  veg_parameter_type, soil_snow_type

USE cbl_qsat_module, ONLY: qsatfjh, qsatfjh2

  TYPE (radiation_type), INTENT(INOUT) :: rad
  TYPE (roughness_type), INTENT(INOUT) :: rough
  TYPE (air_type),       INTENT(INOUT) :: air
  TYPE (met_type),       INTENT(INOUT) :: met
  TYPE (canopy_type),    INTENT(INOUT) :: canopy
  TYPE(soil_snow_type),  INTENT(INOUT) :: ssnow
  TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

  REAL,  INTENT(INOUT)  :: qstvair(mp)        ! sat spec humidity at leaf temperature
  INTEGER, INTENT(IN) :: mp
  REAL, INTENT(IN) :: CRMH2o, Crmair, CTETENA, CTETENB, CTETENC
  REAL, INTENT(IN) :: CLAI_thresh, CCAPP, CTFRZ

  REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
       gbhu,    &  ! forcedConvectionBndryLayerCond
       gbhf        ! freeConvectionBndryLayerCond

  REAL, INTENT(IN), DIMENSION(mp) ::                                       &
       rt0,     &  ! resistance from ground to canopy air
       rhlitt,  &  ! additional litter resistance for heat (=0 by default)
       relitt      ! additional litter resistance for water

  REAL, DIMENSION(mp) ::                                                      &
       rrsw,             & ! recipr. stomatal resistance for water
       rrbw,             & ! recipr. leaf boundary layer resistance for water
       dmah,             & ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
       dmbh,             & ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
       dmch,             & ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
       dmae,             & ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
       dmbe,             & ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
       dmce                ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132

  REAL  :: lower_limit, upper_limit
  REAL, DIMENSION(mp) :: fix_eqn,fix_eqn2

  INTEGER :: j

  !INH: rhlitt=relitt=0. if litter resistance not active but case included
  !dmah through to dmce are not A_{H} through C_{E} as per Eqn 3.40
  !in SCAM documentation but rt0*((1+esp)/rs + 1/rb)*A_{H} etc.
  !
  !changes from v1.4 for %cls package, litter and Or hydrology

  rrbw = SUM(gbhu+gbhf,2)/air%cmolar  ! MJT

  ! leaf stomatal resistance for water
  rrsw = SUM(canopy%gswx,2)/air%cmolar ! MJT

  IF (cable_user%or_evap) THEN
     fix_eqn(:) = rt0(:)*(REAL(ssnow%satfrac(:))/(rt0(:)+REAL(ssnow%rtevap_sat(:))) + &
          (1-REAL(ssnow%satfrac(:)))/(rt0(:)+REAL(ssnow%rtevap_unsat(:))))
     !lakes/ice rtevap=0 and wetfac is .ne. 1
     fix_eqn(:) = ssnow%wetfac(:) * fix_eqn(:)*ssnow%cls(:)   !INH correction. & M.Dekker +d wetfac

     fix_eqn2(:) = rt0(:) / (rt0(:) + REAL(ssnow%rt_qh_sublayer) )

  ELSE  !with INH corrections for litter and cls

     fix_eqn(:) = ssnow%cls(:)*rt0(:)/(rt0(:)+relitt(:))
     WHERE (ssnow%potev>0.) fix_eqn(:)=fix_eqn(:)*ssnow%wetfac(:)
     fix_eqn2(:) = rt0(:)/(rt0(:)+rhlitt(:))

  END IF

  DO j=1,mp

     IF(veg%meth(j) > 0 .AND. canopy%vlaiw(j) > CLAI_THRESH .AND.              &
          rough%hruff(j) > rough%z0soilsn(j) ) THEN

        !   use the dispersion matrix (DM) to find the air temperature
        !   and specific humidity
        !   (Raupach, Finkele and Zhang 1997, pp 17)
        ! leaf boundary layer resistance for water
        ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
        dmah(j) = (rt0(j)+fix_eqn2(j)*rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j))  &
             + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

        ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
        dmbh(j) = (-air%rlam(j)/CCAPP)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

        ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
        dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
             (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*CCAPP)

        ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
        dmae(j) = (-air%epsi(j)*CCAPP/air%rlam(j))*(rt0(j)*rough%rt1(j)) *   &
             (rrbw(j)*rrsw(j))

        ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
        dmbe(j) = ( rt0(j) + fix_eqn(j) * rough%rt1(j) ) *               &
             ( (1.+air%epsi(j) ) * rrsw(j) + rrbw(j) ) +                 &
             ( rt0(j) * rough%rt1(j) ) * ( rrbw(j) * rrsw(j) )

        ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
        ! INH: includes modifications for %cls
        dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
             (canopy%fev(j) + canopy%fes(j)/ssnow%cls(j)) /                   &
             (air%rho(j)*air%rlam(j))

        ! Within canopy air temperature:
        met%tvair(j) = met%tk(j) + ( dmbe(j) * dmch(j) - dmbh(j) * dmce(j) )  &
             / (dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)

        !---set limits for comparisson
        lower_limit =  MIN( ssnow%tss(j), met%tk(j)) - 5.0
        upper_limit =  MAX( ssnow%tss(j), met%tk(j)) + 5.0

        !--- tvair within these limits
        met%tvair(j) = MAX(met%tvair(j) , lower_limit)
        met%tvair(j) = MIN(met%tvair(j) , upper_limit)

        ! recalculate using canopy within temperature
        met%qvair(j) = met%qv(j) + (dmah(j)*dmce(j)-dmae(j)*dmch(j)) /        &
             ( dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
        met%qvair(j) = MAX(0.0,met%qvair(j))

        !---set limits for comparisson

        lower_limit =  MIN( ssnow%qstss(j), met%qv(j))
        upper_limit =  MAX( ssnow%qstss(j), met%qv(j))

        !--- qvair within these limits
        met%qvair(j) =  MAX(met%qvair(j),lower_limit)
        met%qvair(j) =  MIN(met%qvair(j), upper_limit)

        ! Saturated specific humidity in canopy:
        CALL  qsatfjh2( qstvair(j), CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, met%tvair(j)-Ctfrz,met%pmb(j))

        ! Saturated vapour pressure deficit in canopy:
        met%dva(j) = ( qstvair(j) - met%qvair(j) ) *  Crmair/CRMH2o         &
             * met%pmb(j) * 100.
     ENDIF

  ENDDO

END SUBROUTINE within_canopy

END MODULE cable_within_canopy_module
