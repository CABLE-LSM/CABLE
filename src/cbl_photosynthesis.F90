MODULE cbl_photosynthesis_module

IMPLICIT NONE

PUBLIC :: photosynthesis
PRIVATE

CONTAINS
  ! Ticket #56, xleuningz repalced with gs_coeffz
  SUBROUTINE photosynthesis( csxz, cx1z, cx2z, gswminz,                          &
       rdxz, vcmxt3z, vcmxt4z, vx3z,                       &
       vx4z, gs_coeffz, vlaiz, deltlfz, anxz, fwsoilz )
    USE cable_def_types_mod, ONLY : mp, mf, r_2
! maths & other constants
USE cable_other_constants_mod, ONLY : CLAI_THRESH  => LAI_THRESH
USE cable_photo_constants_mod, ONLY : CRGSWC => RGSWC

    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csxz

    REAL, DIMENSION(mp,mf), INTENT(IN) ::                                       &
         cx1z,       & !
         cx2z,       & !
         gswminz,    & !
         rdxz,       & !
         vcmxt3z,    & !
         vcmxt4z,    & !
         vx4z,       & !
         vx3z,       & !
         gs_coeffz,  & ! Ticket #56, xleuningz repalced with gs_coeffz
         vlaiz,      & !
         deltlfz       !

    REAL, DIMENSION(mp,mf), INTENT(INOUT) :: anxz

    ! local variables
    REAL(r_2), DIMENSION(mp,mf) ::                                              &
         coef0z,coef1z,coef2z, ciz,delcxz,                                        &
         anrubiscoz,anrubpz,ansinkz

    REAL, DIMENSION(mp) :: fwsoilz

    REAL, PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
    ! Bonan,LSM version 1.0, p106)

    INTEGER :: i,j

    DO i=1,mp

       IF (SUM(vlaiz(i,:)) .GT. CLAI_THRESH) THEN

          DO j=1,mf

             IF( vlaiz(i,j) .GT. CLAI_THRESH .AND. deltlfz(i,j) .GT. 0.1) THEN

                ! Rubisco limited:
                coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / CRGSWC + gs_coeffz(i,j) * &
                     ( vcmxt3z(i,j) - ( rdxz(i,j)-vcmxt4z(i,j) ) )

                coef1z(i,j) = (1.0-csxz(i,j)*gs_coeffz(i,j)) *                  &
                     (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j))             &
                     + (gswminz(i,j)*fwsoilz(i)/CRGSWC)*(cx1z(i,j)-csxz(i,j)) &
                     - gs_coeffz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0      &
                     + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j) ) )


                coef0z(i,j) = -(1.0-csxz(i,j)*gs_coeffz(i,j)) *                 &
                     (vcmxt3z(i,j)*cx2z(i,j)/2.0                       &
                     + cx1z(i,j)*( rdxz(i,j)-vcmxt4z(i,j ) ) )         &
                     -( gswminz(i,j)*fwsoilz(i)/CRGSWC ) * cx1z(i,j)*csxz(i,j)


                ! kdcorbin,09/10 - new calculations
                IF( ABS(coef2z(i,j)) .GT. 1.0e-9 .AND. &
                     ABS(coef1z(i,j)) .LT. 1.0e-9) THEN

                   ! no solution, give it a huge number as
                   ! quadratic below cannot handle zero denominator
                   ciz(i,j) = 99999.0

                   anrubiscoz(i,j) = 99999.0 ! OR do ciz=0 and calc anrubiscoz

                ENDIF

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                     ABS( coef1z(i,j) ) >= 1e-9 ) THEN

                   ! same reason as above
                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)

                   ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )

                   anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j) / 2.0 ) / &
                        ( ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) -   &
                        rdxz(i,j)

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0 * coef0z(i,j)              &
                        * coef2z(i,j)

                   ciz(i,j) = ( -coef1z(i,j) + SQRT( MAX( 0.0_r_2 ,             &
                        delcxz(i,j) ) ) ) / ( 2.0*coef2z(i,j) )

                   ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )   ! must be positive, why?

                   anrubiscoz(i,j) = vcmxt3z(i,j) * ( ciz(i,j) - cx2z(i,j)      &
                        / 2.0)  / ( ciz(i,j) + cx1z(i,j) ) +       &
                        vcmxt4z(i,j) - rdxz(i,j)

                ENDIF

                ! RuBP limited:
                coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / CRGSWC + gs_coeffz(i,j) &
                     * ( vx3z(i,j) - ( rdxz(i,j) - vx4z(i,j) ) )

                coef1z(i,j) = ( 1.0 - csxz(i,j) * gs_coeffz(i,j) ) *            &
                     ( vx3z(i,j) + vx4z(i,j) - rdxz(i,j) )             &
                     + ( gswminz(i,j)*fwsoilz(i) / CRGSWC ) *          &
                     ( cx2z(i,j) - csxz(i,j) ) - gs_coeffz(i,j)        &
                     * ( vx3z(i,j) * cx2z(i,j) / 2.0 + cx2z(i,j) *     &
                     ( rdxz(i,j) - vx4z(i,j) ) )

                coef0z(i,j) = -(1.0-csxz(i,j)*gs_coeffz(i,j)) *   &
                     (vx3z(i,j)*cx2z(i,j)/2.0                          &
                     + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j)))                &
                     - (gswminz(i,j)*fwsoilz(i)/CRGSWC)*cx2z(i,j)*csxz(i,j)


                !Ticket #117 - initialize at all times
                ciz(i,j) = 99999.0
                anrubpz(i,j)  = 99999.0

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                     ABS( coef1z(i,j) ) >= 1.e-9) THEN

                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)

                   ciz(i,j) = MAX(0.0_r_2,ciz(i,j))

                   anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /          &
                        (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF ( ABS( coef2z(i,j)) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)

                   ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j))))     &
                        /(2.0*coef2z(i,j))

                   ciz(i,j) = MAX(0.0_r_2,ciz(i,j))

                   anrubpz(i,j)  = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /         &
                        (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)

                ENDIF

                ! Sink limited:
                coef2z(i,j) = gs_coeffz(i,j)

                coef1z(i,j) = gswminz(i,j)*fwsoilz(i)/CRGSWC + gs_coeffz(i,j)   &
                     * (rdxz(i,j) - 0.5*vcmxt3z(i,j))                  &
                     + effc4 * vcmxt4z(i,j) - gs_coeffz(i,j)           &
                     * csxz(i,j) * effc4 * vcmxt4z(i,j)

                coef0z(i,j) = -( gswminz(i,j)*fwsoilz(i)/CRGSWC )*csxz(i,j)*effc4 &
                     * vcmxt4z(i,j) + ( rdxz(i,j)                      &
                     - 0.5 * vcmxt3z(i,j)) * gswminz(i,j)*fwsoilz(i)/CRGSWC

                ! no solution, give it a huge number
                IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND.                           &
                     ABS( coef1z(i,j)) < 1.0e-9 ) THEN

                   ciz(i,j) = 99999.0
                   ansinkz(i,j)  = 99999.0

                ENDIF

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                     ABS( coef1z(i,j) ) >= 1.e-9 ) THEN

                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                   ansinkz(i,j)  = ciz(i,j)

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)

                   ciz(i,j) = (-coef1z(i,j)+SQRT (MAX(0.0_r_2,delcxz(i,j)) ) )  &
                        / ( 2.0 * coef2z(i,j) )

                   ansinkz(i,j) = ciz(i,j)

                ENDIF

                ! minimal of three limited rates
                anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))


             ENDIF

          ENDDO

       ENDIF

    ENDDO



  END SUBROUTINE photosynthesis


END MODULE cbl_photosynthesis_module
