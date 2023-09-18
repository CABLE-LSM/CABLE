MODULE snow_processes_soil_thermal_mod

USE cbl_ssnow_data_mod

PUBLIC snow_processes_soil_thermal

CONTAINS

!SUBROUTINE snow_processes_soil_thermal(dels,ssnow,soil,veg,canopy,met,bal,snowmlt) ! replaced by rk4417 - phase2
! snowmlt is not used in cable_gw_hydro_module and as far as I can tell serves no purpose in the code - rk4417

SUBROUTINE snow_processes_soil_thermal(dels,ssnow,soil,veg,canopy,met,bal)
!* calculate snow processes and thermal soil
  
USE snowl_adjust_mod,             ONLY: snowl_adjust                     
USE snowCheck_mod,                ONLY: snowCheck                         
USE snow_melting_mod,             ONLY: snow_melting
USE snow_accum_mod,               ONLY: snow_accum
USE snowdensity_mod,              ONLY: snowDensity
USE GWstempv_mod,                 ONLY: GWstempv

IMPLICIT NONE

    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    TYPE (balances_type), INTENT(INOUT)      :: bal
    REAL, DIMENSION(mp)                      :: snowmlt !track snow melt 
!    REAL, DIMENSION(:),  INTENT(INOUT)       :: snowmlt  ! replaced by rk4417 - phase2
    INTEGER             :: k,i

    snowmlt = 0.0  ! inserted by rk4417 - phase2
    
    CALL snowcheck (dels, ssnow, soil, met )

    CALL snowdensity (dels, ssnow, soil)

    CALL snow_accum (dels, canopy, met, ssnow, soil )

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add snow melt to global snow melt variable:
    ssnow%smelt(:) = snowmlt(:)
    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ssnow, canopy )

!    IF (cable_user%gw_model) CALL GWstempv(dels, canopy, ssnow, soil) ! replaced by rk4417 - phase2
    CALL GWstempv(dels, canopy, ssnow, soil)

    !do the soil and snow melting, freezing prior to water movement
    DO i=1,mp
       ssnow%tss(i) =  (1-ssnow%isflag(i))*ssnow%tgg(i,1) + ssnow%isflag(i)*ssnow%tggsn(i,1)
    END DO

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add new snow melt to global snow melt variable:
    ssnow%smelt(:) = ssnow%smelt(:) + snowmlt(:)

END SUBROUTINE snow_processes_soil_thermal

END MODULE snow_processes_soil_thermal_mod
