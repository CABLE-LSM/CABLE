MODULE snow_aging_mod

  !This routine evaluates the effective age of any snow pack
  !It is called from cbm module and the output (SnowAge)
  !is used in cbl_snow_albedo

PUBLIC snow_aging

CONTAINS

  SUBROUTINE snow_aging(SnowAge,mp,dels,SnowDepth,SnowODepth,SnowTemp, &
       SoilTemp,SnowFlag_3L, surface_type,soil_type)

    !may need to eliminate this USE data statement  
    USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ

    IMPLICIT NONE

    REAL, INTENT(IN OUT)   :: SnowAge(mp)
    INTEGER, INTENT(IN) :: mp
    REAL, INTENT(IN)    :: dels
    REAL, INTENT(IN)    :: SnowDepth(mp)
    REAL, INTENT(IN)    :: SnowODepth(mp)
    REAL, INTENT(IN)    :: SnowTemp(mp)
    REAL, INTENT(IN)    :: SoilTemp(mp)
    INTEGER, INTENT(IN) :: SnowFlag_3L(mp)
    INTEGER, INTENT(IN) :: surface_type(mp) 
    INTEGER, INTENT(IN) :: soil_type(mp) 

    !hard wired index to be eliminated
    INTEGER, PARAMETER :: perm_ice = 9
    !model parameter shared across subroutines -> cable_phys_constants
    REAL, PARAMETER :: snow_depth_thresh = 1.0

    !working variables - converted to scalars
    REAL ::                                                           &
         ar1,     &  ! factor for crystal growth  (-ve)
         ar2,     &  ! factor for freezing of melt water
         ar3,     &  ! factor for accumulation of dirt
         dnsnow,  &  ! depth of new snow albedo
         dtau,    &  !change in effective snow age
         tmp         !local effective surface temperauture

    INTEGER :: i        !looping variable
    
    DO i=1,mp
       IF (SnowDepth(i)>snow_depth_thresh) THEN

          !depth of new snow (in cm H20)
          dnsnow = MIN (1.0, 0.1 * MAX( 0.0, SnowDepth(i) - SnowODepth(i) ) )

          ! Snow age depends on snow crystal growth, freezing of melt water,
          ! accumulation of dirt and amount of new snow.
          tmp = SnowFlag_3L(i) * SnowTemp(i) + ( 1 - SnowFlag_3L(i) ) * SoilTemp(i)
          tmp = MIN( tmp, CTFRZ )
          ar1 = 5000.0 * (1.0 / (CTFRZ-0.01) - 1.0 / tmp) ! crystal growth  (-ve)
          ar2 = 10.0 * ar1                             ! freezing of melt water

          IF (soil_type(i) == perm_ice) THEN
             ! permanent ice case
             ar3 = 0.0000001
             
             !  NB. dsnow =1,assumes pristine snow; ignores soot etc. ALTERNATIVELY,
             !dnsnow = max (dnsnow(i), 0.5) !increase refreshing of snow in Antarctic
             dnsnow = 1.0

          ELSE
             ! accumulation of dirt
             ar3 = 0.1

          END IF

          !update snow age
          dtau = 1.0e-6 * (EXP( ar1 ) + EXP( ar2 ) + ar3 ) * dels
          SnowAge(i) = MAX(0.0,(SnowAge(i)+dtau)*(1.0-dnsnow))

       END IF

    END DO

    RETURN

  END SUBROUTINE snow_aging

END MODULE snow_aging_mod
