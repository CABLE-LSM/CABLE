MODULE cbl_zetar_module

IMPLICIT NONE

PUBLIC update_zetar

PRIVATE

CONTAINS

SUBROUTINE update_zetar( mp, NITER, canopy_zetar, iter, nrb, CVONK, CGRAV, CCAPP,  &
                   CLAI_THRESH, CZETmul, CZETPOS, CZETNEG,          &
                   cable_user_soil_struc, air_rho, met_tk,  met_fsd, &
                   rough_zref_tq, rough_hruff, rough_term6a, rough_z0soilsn,   &
                   canopy_vlaiw, canopy_zetash,  canopy_us, &
                   canopy_fh, canopy_fe, canopy_fhs, canopy_fes ) 
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: NITER
INTEGER, INTENT(IN) :: nrb

REAL, INTENT(OUT) :: canopy_zetar(mp, NITER)
REAL, INTENT(OUT) :: canopy_zetash(mp, NITER)
INTEGER, INTENT(IN) :: iter

! constants
REAL, INTENT(IN) :: CVONK, CGRAV, CCAPP, CLAI_THRESH, CZETmul, CZETPOS, CZETNEG
CHARACTER, INTENT(IN)  :: cable_user_soil_struc

REAL, INTENT(IN) :: air_rho(mp)
REAL, INTENT(IN) :: met_tk(mp)
REAL, INTENT(IN) :: met_fsd(mp,nrb)
REAL, INTENT(IN) :: rough_zref_tq(mp)
REAL, INTENT(IN) :: rough_hruff(mp)
REAL, INTENT(IN) :: rough_term6a(mp)
REAL, INTENT(IN) :: rough_z0soilsn(mp)
REAL, INTENT(IN) :: canopy_vlaiw(mp)
REAL, INTENT(IN) :: canopy_us(mp)
REAL, INTENT(IN) :: canopy_fh(mp)
REAL, INTENT(IN) :: canopy_fe(mp)
REAL, INTENT(IN) :: canopy_fhs(mp)
REAL, INTENT(IN) :: canopy_fes(mp)

!local vars
INTEGER :: iterplus
INTEGER :: j

! monin-obukhov stability parameter zetar=zref/l
! recompute zetar for the next iteration, except on last iteration
IF (iter < NITER) THEN ! dont compute zetar on the last iter

  iterplus = MAX(iter+1,2)
  canopy_zetar(:,iterplus) = -( CVONK * CGRAV * rough_zref_tq *              &
       ( canopy_fh + 0.07 * canopy_fe ) ) /          &
       ( air_rho * CCAPP * met_tk * canopy_us**3 )

  ! stability parameter at shear height: needed for Harman in-canopy stability correction
  IF (cable_user_soil_struc=='sli') THEN
     WHERE (canopy_vlaiw > CLAI_THRESH .AND. rough_hruff > rough_z0soilsn)
        canopy_zetash(:,iterplus) = -(Cvonk*Cgrav*(0.1*rough_hruff)*(canopy_fhs+0.07*REAL(canopy_fes)))/ &
             MAX( (air_rho*Ccapp*met_tk*(canopy_us*rough_term6a)**3), 1.e-12)
     ELSEWHERE
        canopy_zetash(:,iterplus) = canopy_zetash(:,iter)
     ENDWHERE
  ENDIF

  ! case NITER=2: final zetar=CZETmul*zetar(2) (compute only when iter=1)
  IF (NITER == 2) THEN

     canopy_zetar(:,2) = CZETmul * canopy_zetar(:,2)

     ! stability parameter at shear height: needed for Harman in-canopy stability correction
     IF (cable_user_soil_struc=='sli') THEN
        canopy_zetash(:,2) = CZETmul * canopy_zetash(:,2)
     ENDIF


     DO j=1,mp
        IF ( (met_fsd(j,1)+met_fsd(j,2))  ==  0.0 ) &
             canopy_zetar(j,2) = 0.5 * canopy_zetar(j,2)
     ENDDO

  END IF

  ! constrain zeta to CZETPOS and CZETNEG (set in param0)

  ! zetar too +
  canopy_zetar(:,iterplus) = MIN(CZETPOS,canopy_zetar(:,iterplus))
  !jhan: to get past rigorous build - however (:,i) cant be compared
  !if ( canopy_zetash(:,iterplus) .NE. CZETPOS ) &
  IF (cable_user_soil_struc=='sli') &
       canopy_zetash(:,iterplus) = MIN(CZETPOS,canopy_zetash(:,iterplus))

  ! zetar too -
  canopy_zetar(:,iterplus) = MAX(CZETNEG,canopy_zetar(:,iterplus))
  !if ( canopy_zetash(:,iterplus) .NE. CZETNEG ) &
  IF (cable_user_soil_struc=='sli') &
       canopy_zetash(:,iterplus) = MAX(CZETNEG,canopy_zetash(:,iterplus))

END IF ! (iter < NITER)

RETURN
END SUBROUTINE update_zetar

END MODULE cbl_zetar_module
