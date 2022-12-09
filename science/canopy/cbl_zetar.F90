MODULE cbl_zetar_module
  !* This MODULE contains the SUBROUTINE [[update_zetar]] needed to update
  !  the value of the stability parameter `canopy%zetar`=\(\xi\). 

IMPLICIT NONE

PUBLIC update_zetar

PRIVATE

CONTAINS

SUBROUTINE update_zetar( mp, iterplus, NITER, canopy_zetar, iter, nrb, CVONK, CGRAV, CCAPP,  &
                   CLAI_THRESH, CZETmul, CZETPOS, CZETNEG,          &
                   cable_user_soil_struc, air_rho, met_tk,  met_fsd, &
                   rough_zref_tq, rough_hruff, rough_term6a, rough_z0soilsn,   &
                   canopy_vlaiw, canopy_zetash,  canopy_us, &
                   canopy_fh, canopy_fe, canopy_fhs, canopy_fes )
  !*## Purpose
  !
  ! This SUBROUTINE updates the value of the stability parameter \(\xi\)
  ! during the iteration loop of the Monin-Obukhov (MO) similarity theory in [[define_canopy]].  
  ! \(\xi\) quantifies the role that the surface fluxes play in setting the 
  ! efficiency of turbulent transfer from the land to the atmosphere, and hence
  ! the aerodynamic component of the resistance network for those same surface 
  ! fluxes (an implicit problem which requires iteration to solve).
  !
  ! This SUBROUTINE forms part of the codebase to evaluate the surface
  ! energy balance on a sub-diurnal basis (i.e. every CABLE time step).
  ! It resides in the canopy science directory.
  !  
  !## Method
  !
  ! The two outputs of the SUBROUTINE are:
  !
  ! - `canopy_zetar`, the local (in space, time and by iteration counter) value of \(\xi\)
  !  (Equation 9). It is evaluated from the total land (soil+canopy)
  !  surface fluxes of momentum, sensible heat and latent heat.
  ! - `canopy_zetash` is the equivalent variable evaluated from the soil
  !  contribution to those fluxes only. `canopy_zetash` is used in conjunction
  !  with the [[sli_main_mod]] soil model to moderate the fluxes from the soil
  !  underneath a canopy.
  !
  ! `canopy_zetar` and `canopy_zetash` are initialised to `CZETA0`=0 in
  ! [[define_canopy]] and updated `NITER`(>1) times during the calculation 
  ! of the energy balance. The value of the variables at each iteration are 
  ! stored in memory to aid in the diagnosis of convergence.
  !
  ! A special case applies if `NITER`=2.  
  ! `canopy_zetar` and `canopy_zetash` are also bounded by the interval 
  ! `[CZETNEG, CZETPOS]`.  
  !
  ! The outputs `canopy_zetar` and `canopy_zetash` are known as `canopy%zetar` 
  ! and `canopy%zetash` elsewhere in the code. 
  ! `NITER`(=4) is defined in [[cable_types_mod]]; `CZETMUL`, `CZET0`, 
  ! `CZETPOS` and `CZETNEG` in [[cable_phys_constants_mod]].
  !
  !## References
  !
  ! [Kowalczyk et al. (2006)](http://www.cmar.csiro.au/e-print/open/kowalczykea_2006a.pdf)
  ! - section 3.1, equations 1-9. 
  ! 
  !  <br></br>
  !
  ! **WARNING** The INTENT statements for `canopy_zetar` and `canopy_zetash` 
  ! need to be INTENT(INOUT): currently previous values are reset at each call 
  ! of the subroutine. It means the initialisations in [[define_canopy]] are 
  ! useless (the ITER=1 values are lost) and for SLI, on bare soils, 
  !`canopy_zetash` gets crazy values for all iterations. 
  !

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp     !! number of land points (-)
INTEGER, INTENT(IN) :: NITER  !! number of MO-iterations (-)
INTEGER, INTENT(IN) :: nrb    !! number of radiation bands (-)
INTEGER, INTENT(OUT)  :: iterplus

REAL, INTENT(OUT) :: canopy_zetar(mp, NITER)  !!stability parameter \(\xi\) (-)
REAL, INTENT(OUT) :: canopy_zetash(mp, NITER) !!stability parameter for soil (-)
INTEGER, INTENT(IN) :: iter                   !! iteration counter (-)

! constants
REAL, INTENT(IN) :: CVONK, CGRAV, CCAPP, CLAI_THRESH, CZETmul, CZETPOS, CZETNEG
CHARACTER, INTENT(IN)  :: cable_user_soil_struc !! name of soil model used

REAL, INTENT(IN) :: air_rho(mp)        !! air density (kg m\(^{-3}\))
REAL, INTENT(IN) :: met_tk(mp)         !! reference level air temperature (K)
REAL, INTENT(IN) :: met_fsd(mp,nrb)    !! downwelling shortwave (Wm\(^{-2}\))  
REAL, INTENT(IN) :: rough_zref_tq(mp)  !! reference height for T and q (m)
REAL, INTENT(IN) :: rough_hruff(mp)    !! height of canopy (above snow) (m)
REAL, INTENT(IN) :: rough_term6a(mp)   !! term from [[ruff_resist]] (-)
REAL, INTENT(IN) :: rough_z0soilsn(mp) !! roughness length of soil or snow (m)
REAL, INTENT(IN) :: canopy_vlaiw(mp)   !! canopy leaf area (m\(^2\)m\(^{-2}\))
REAL, INTENT(IN) :: canopy_us(mp)      !! friction velocity (ms\(^{-1}\))
REAL, INTENT(IN) :: canopy_fh(mp)      !! sensible heat flux (Wm\(^{-2}\))
REAL, INTENT(IN) :: canopy_fe(mp)      !! latent heat flux (Wm\(^{-2}\))
REAL, INTENT(IN) :: canopy_fhs(mp)     !! soil sensible heat flux (Wm\(^{-2}\)) 
REAL, INTENT(IN) :: canopy_fes(mp)     !! soil latent heat flux (Wm\(^{-2}\))

!local vars
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
