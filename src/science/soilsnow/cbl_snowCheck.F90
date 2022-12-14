MODULE snowcheck_mod

USE cbl_ssnow_data_mod

PUBLIC  snowcheck

CONTAINS

SUBROUTINE snowcheck(dels, ssnow, soil, met )

    USE cable_common_module

IMPLICIT NONE
    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE(met_type),       INTENT(INOUT) :: met ! all met forcing

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters

    INTEGER :: k,j

    DO j=1,mp

       IF( ssnow%snowd(j) <= 0.0 ) THEN

          ssnow%isflag(j) = 0
          ssnow%ssdn(j,:) = 120.0
          ssnow%ssdnn(j) = 120.0
          ssnow%tggsn(j,:) = CTFRZ
          ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)

          ssnow%sdepth(j,2) = 0.
          ssnow%sdepth(j,3) = 0.

          ssnow%smass(j,1) = ssnow%snowd(j)
          ssnow%smass(j,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
          ssnow%smass(j,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007

          ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
       ELSEIF( ssnow%snowd(j) < snmin * ssnow%ssdnn(j) ) THEN

          IF( ssnow%isflag(j) == 1 ) THEN
             ssnow%ssdn(j,1) = ssnow%ssdnn(j)
             ssnow%tgg(j,1) = ssnow%tggsn(j,1)
          ENDIF

          ssnow%isflag(j) = 0
          ssnow%ssdnn(j) = MIN( 400.0, MAX( 120.0, ssnow%ssdn(j,1) ) )

          ssnow%tggsn(j,:) = MIN( CTFRZ,ssnow%tgg(j,1) )

          ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)
          ssnow%sdepth(j,2) = 0.0
          ssnow%sdepth(j,3) = 0.0

          ssnow%smass(j,1) = ssnow%snowd(j)
          ssnow%smass(j,2) = 0.0
          ssnow%smass(j,3) = 0.0

          ssnow%ssdn(j,:) = ssnow%ssdnn(j)



       ELSE ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
          ! sufficient snow now for 3 layer snowpack

          IF( ssnow%isflag(j) == 0 ) THEN

             ssnow%tggsn(j,:) = MIN( CTFRZ, ssnow%tgg(j,1) )

             ssnow%ssdn(j,2) = ssnow%ssdn(j,1)
             ssnow%ssdn(j,3) = ssnow%ssdn(j,1)


             ssnow%sdepth(j,1) = ssnow%t_snwlr(j)

             ssnow%smass(j,1)  =  ssnow%t_snwlr(j) * ssnow%ssdn(j,1)

             ssnow%smass(j,2)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.4
             ssnow%smass(j,3)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.6

             ssnow%sdepth(j,2) = ssnow%smass(j,2) / ssnow%ssdn(j,2)
             ssnow%sdepth(j,3) = ssnow%smass(j,3) / ssnow%ssdn(j,3)

             ssnow%ssdnn(j) = ( ssnow%ssdn(j,1) * ssnow%smass(j,1) +            &
                  ssnow%ssdn(j,2) * ssnow%smass(j,2) +             &
                  ssnow%ssdn(j,3) * ssnow%smass(j,3) )             &
                  / ssnow%snowd(j)

          ENDIF

          ssnow%isflag(j) = 1

       ENDIF ! END: IF( ssnow%snowd(j) <= 0.0 ) THEN


    ENDDO ! END: DO j=1,mp

END SUBROUTINE snowcheck 


END MODULE snowcheck_mod
