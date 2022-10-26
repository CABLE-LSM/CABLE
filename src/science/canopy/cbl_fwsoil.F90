MODULE cbl_fwsoil_module

IMPLICIT NONE

PUBLIC :: fwsoil_calc_std
PUBLIC :: fwsoil_calc_non_linear
PUBLIC :: fwsoil_calc_Lai_Ktaul
PUBLIC :: fwsoil_calc_sli
PRIVATE

CONTAINS

SUBROUTINE fwsoil_calc_std(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    USE cable_common_module, ONLY : cable_user
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    REAL, DIMENSION(mp) :: rwater ! soil water availability

    !note even though swilt_vec is defined in default model it is r_2
    !and even using real(_vec) gives results different from trunk (rounding
    !errors)

    IF (.NOT.cable_user%gw_model) THEN

       rwater = MAX(1.0e-9,                                                    &
            SUM(veg%froot * MAX(1.0e-9,MIN(1.0, REAL(ssnow%wb) -                   &
            SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))

    ELSE
       rwater = MAX(1.0e-9,                                                    &
            SUM(veg%froot * MAX(1.0e-9,MIN(1.0, REAL((ssnow%wbliq -                 &
            soil%swilt_vec)/(soil%sfc_vec-soil%swilt_vec)) )),2) )

    ENDIF

    ! Remove vbeta #56
    IF(cable_user%GS_SWITCH == 'medlyn') THEN
       fwsoil = MAX(1.0e-4,MIN(1.0, rwater))
    ELSE
       fwsoil = MAX(1.0e-9,MIN(1.0, veg%vbeta * rwater))
    ENDIF


  END SUBROUTINE fwsoil_calc_std

  ! ------------------------------------------------------------------------------

  SUBROUTINE fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    REAL, DIMENSION(mp) :: rwater ! soil water availability
    REAL, DIMENSION(mp,3)          :: xi, ti, si
    INTEGER :: j

    rwater = MAX(1.0e-9,                                                    &
         SUM(veg%froot * MAX(0.0,MIN(1.0, REAL(ssnow%wb) -                   &
         SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))

    fwsoil = 1.

    rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)

    xi(:,1) = soil%swilt
    xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
    xi(:,3) = soil%sfc

    ti(:,1) = 0.
    ti(:,2) = 0.9
    ti(:,3) = 1.0

    si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) *                       &
         (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))

    si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) *                       &
         (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))

    si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) *                       &
         (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))

    DO j=1,mp
       IF (rwater(j) < soil%sfc(j) - 0.02)                                      &
            fwsoil(j) = MAX(0.,MIN(1., ti(j,1)*si(j,1) +                          &
            ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))

    ENDDO

  END SUBROUTINE fwsoil_calc_non_linear

  ! ------------------------------------------------------------------------------

  ! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
  SUBROUTINE fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    INTEGER   :: ns
    REAL, PARAMETER ::rootgamma = 0.01   ! (19may2010)
    REAL, DIMENSION(mp)  :: dummy, normFac
    !--- local level dependent rwater
    REAL, DIMENSION(mp,ms)  :: frwater

    fwsoil(:) = 0.0
    normFac(:) = 0.0

    DO ns=1,ms

       dummy(:) = rootgamma/MAX(1.0e-3_r_2,ssnow%wb(:,ns)-soil%swilt(:))

       frwater(:,ns) = MAX(1.0e-4_r_2,((ssnow%wb(:,ns)-soil%swilt(:))/soil%ssat(:)) &
            ** dummy)

       fwsoil(:) = MIN(1.0,MAX(fwsoil(:),frwater(:,ns)))

       normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)

    ENDDO

  END SUBROUTINE fwsoil_calc_Lai_Ktaul

  ! ------------------------------------------------------------------------------
  SUBROUTINE fwsoil_calc_sli(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    REAL, DIMENSION(mp,ms):: tmp2d1, tmp2d2, delta_root, alpha2a_root, alpha2_root
    ! Lai and Katul formulation for root efficiency function  vh 17/07/09
    alpha2a_root = MAX(ssnow%wb-soil%swilt_vec, 0.001_r_2)/(soil%ssat_vec)
    tmp2d1 = ssnow%wb -soil%swilt_vec
    tmp2d2 = SPREAD(veg%gamma,2,ms)/tmp2d1*LOG(alpha2a_root)
    WHERE ((tmp2d1>0.001) .AND. (tmp2d2 > -10.0))
       alpha2_root = EXP(tmp2d2)
    ELSEWHERE
       alpha2_root = 0.0
    ENDWHERE

    WHERE (veg%froot>0.0)
       delta_root = 1.0
    ELSEWHERE
       delta_root = 0.0
    ENDWHERE

    fwsoil  = MAXVAL(alpha2_root*delta_root, 2)
    fwsoil  = MAX(0.0, fwsoil)

  END SUBROUTINE fwsoil_calc_sli

 
END MODULE cbl_fwsoil_module
