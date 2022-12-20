MODULE surfbv_mod

USE cbl_ssnow_data_mod

PUBLIC  surfbv

CONTAINS 

SUBROUTINE surfbv (dels, met, ssnow, soil, veg, canopy )

USE smoisturev_mod,               ONLY: smoisturev
    USE cable_common_module
IMPLICIT NONE

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(canopy_type), INTENT(IN)       :: canopy

    TYPE(met_type),       INTENT(INOUT) :: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow  ! soil+snow variables

    TYPE(veg_parameter_type),  INTENT(IN)     :: veg
    TYPE(soil_parameter_type), INTENT(INOUT)  :: soil  ! soil parameters

    !jhan:cable.nml
    INTEGER            :: nglacier  ! 0 original, 1 off, 2 new Eva

    REAL, DIMENSION(mp) ::                                                      &
         rnof5,      & !
         sfact,      & !
         sgamm,      & !
         smasstot,   & !
         talb,       & ! snow albedo
         tmp           ! temporary value
    REAL(r_2), DIMENSION(mp) :: xxx

    REAL, DIMENSION(mp,0:3) :: smelt1

    REAL :: wb_lake_T, rnof2_T, ratio
    INTEGER :: k,j

    IF( cable_runtime%UM ) THEN
       nglacier = 0
    ELSE
       nglacier = 2
    ENDIF
    IF( cable_runtime%esm15 ) nglacier = 2

    CALL smoisturev( dels, ssnow, soil, veg )

    DO k = 1, ms
       xxx = REAL( soil%ssat,r_2 )
       ssnow%rnof1 = ssnow%rnof1 + REAL( MAX( ssnow%wb(:,k) - xxx, 0.0_r_2 )  &
            * 1000.0 )  * soil%zse(k)
       ssnow%wb(:,k) = MAX( 0.01_r_2, MIN( ssnow%wb(:,k), xxx ) )
    END DO

    ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
    ! account for soil/ice cracking
    ! fracm = MIN(0.2, 1. - MIN(1., ssnow%wb(:,ms) / soil%sfc ) )
    ! ssnow%wb(:,ms) = ssnow%wb(:,ms) &
    !                  + fracm*ssnow%rnof1/(1000.0*soil%zse(ms))
    ! ssnow%rnof1 = (1. - fracm) * ssnow%rnof1

    ! Scaling  runoff to kg/m^2/s to match rest of the model
    !jhan:replace nested wheres

    !---  glacier formation
    rnof5= 0.

    IF (nglacier == 2) THEN

       smelt1=0.
       WHERE( ssnow%snowd > max_glacier_snowd )

          rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

          !---- change local tg to account for energy - clearly not best method
          WHERE( ssnow%isflag == 0 )
             smasstot = 0.0
             ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * CHLF                    &
                  / REAL( ssnow%gammzz(:,1) )
             ssnow%snowd = ssnow%snowd - rnof5
          ELSEWHERE
             smasstot = ssnow%smass(:,1) + ssnow%smass(:,2) + ssnow%smass(:,3)
          END WHERE

       END WHERE

       DO k = 1, 3

          WHERE( ssnow%snowd > max_glacier_snowd  .AND.  ssnow%isflag > 0 )
             sgamm = ssnow%ssdn(:,k) * Ccgsnow * ssnow%sdepth(:,k)
             smelt1(:,k) = MIN( rnof5 * ssnow%smass(:,k) / smasstot,            &
                  0.2 * ssnow%smass(:,k) )
             ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)
             ssnow%snowd = ssnow%snowd - smelt1(:,k)
          END WHERE

       END DO

       WHERE( ssnow%isflag > 0 ) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)

    END IF

    !  Rescale drainage to remove water added to lakes (wb_lake)
    ssnow%sinfil = 0.0
    WHERE( veg%iveg == 16 )
       ssnow%sinfil  = MIN( ssnow%rnof1, ssnow%wb_lake ) ! water that can be extracted from the rnof1
       ssnow%rnof1   = MAX( 0.0, ssnow%rnof1 - ssnow%sinfil )
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
       ssnow%sinfil  = MIN( ssnow%rnof2, ssnow%wb_lake ) ! water that can be extracted from the rnof2
       ssnow%rnof2   = MAX( 0.0, ssnow%rnof2 - ssnow%sinfil )
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
       xxx = MAX(0.0_r_2, (ssnow%wb(:,ms) - REAL(soil%sfc(:),r_2))*soil%zse(ms)*1000.0)
       ssnow%sinfil  = MIN( REAL(xxx), ssnow%wb_lake )
       ssnow%wb(:,ms) = ssnow%wb(:,ms) - REAL(ssnow%sinfil / (soil%zse(ms)*1000.0), r_2)
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
       xxx = MAX(0.0_r_2, (ssnow%wb(:,ms) - 0.5*(soil%sfc + soil%swilt))*soil%zse(ms)*1000.0)
       ssnow%sinfil  = MIN( REAL(xxx), ssnow%wb_lake )
       ssnow%wb(:,ms) = ssnow%wb(:,ms) - ssnow%sinfil / (soil%zse(ms)*1000.0)
       ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
    ENDWHERE

    !wb_lake_T = sum( ssnow%wb_lake )
    !rnof2_T = sum( ssnow%rnof2 )
    !ratio = min( 1., wb_lake_T/max(rnof2_T,1.))
    !ssnow%rnof2 = ssnow%rnof2 - ratio*ssnow%rnof2
    !ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ratio*ssnow%rnof2)

    !  Rescale drainage to remove water added to lakes (wb_lake)
    !wb_lake_T = 0.0
    !rnof2_T = 0.
    !DO j=1,mp
    !   IF( ssnow%wb_lake(j) >  0.0 ) wb_lake_T = wb_lake_T + ssnow%wb_lake(j)
    !   rnof2_T = rnof2_T + ssnow%rnof2(j)
    !END DO
    !ratio = min( 1., wb_lake_T/max(rnof2_T,1.))
    !ssnow%rnof2 = ssnow%rnof2 - ratio*ssnow%rnof2
    !ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ratio*ssnow%rnof2)

    ssnow%rnof1 = ssnow%rnof1 / dels + rnof5/dels
    ssnow%rnof2 = ssnow%rnof2 / dels
    ssnow%runoff = ssnow%rnof1 + ssnow%rnof2

  END SUBROUTINE surfbv

END MODULE surfbv_mod

