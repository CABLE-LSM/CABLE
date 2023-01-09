MODULE cable_climate_type_mod

IMPLICIT NONE

  ! Climate data:
TYPE climate_type

  INTEGER :: nyear_average = 20
  INTEGER :: nday_average  = 31
  !      INTEGER, POINTER ::                                                  &
  INTEGER ::                                                                  &
       nyears, & ! number of years in climate record
       doy ! day of year

  INTEGER, DIMENSION(:), POINTER ::                                           &
       chilldays, &   ! length of chilling period (period with T<5deg)
       iveg, &        ! potential vegetation type based on climatic constraints
       biome

  REAL, DIMENSION(:), POINTER ::                                              &
       dtemp,        & ! daily temperature
       dmoist,        & ! daily moisture availability
       mtemp,       & ! mean temperature over the last 31 days
       qtemp,       & ! mean temperature over the last 91 days
       mmoist,        & ! monthly moisture availability
       mtemp_min,   & ! minimum monthly temperature
       mtemp_max,   & ! maximum monhtly temperature
       qtemp_max,& ! mean temperature of the warmest quarter (so far this year)
       qtemp_max_last_year, & ! mean temperature of the warmest quarter (last calendar year)
       mtemp_min20,   & ! minimum monthly temperature, averaged over 20 y
       mtemp_max20,   & ! maximum monhtly temperature, averaged over 20 y
       atemp_mean,  & ! annual average temperature
       agdd5,                                                                 &
       gdd5,       & ! growing degree day sum relative to 5deg base temperature
       agdd0,        & !
       gdd0,       & ! growing degree day sum relative to 0deg base temperature
       alpha_PT,    & ! ratio of annual evap to annual PT evap
       evap_PT,    & ! annual PT evap [mm]
       aevap , &       ! annual evap [mm]
       alpha_PT20

  REAL, DIMENSION(:,:), POINTER ::                                            &
       mtemp_min_20, & ! mimimum monthly temperatures for the last 20 y
       mtemp_max_20, & ! maximum monthly temperatures for the last 20 y
       dtemp_31 , &    ! daily temperature for the last 31 days
       dmoist_31 , &    ! daily moisture availability for the last 31 days
       alpha_PT_20, &      ! priestley Taylor Coefft for last 20 y
       dtemp_91     ! daily temperature for the last 91 days

END TYPE climate_type

!Instantiation:
TYPE(climate_type) :: climate_cbl

END MODULE cable_climate_type_mod
