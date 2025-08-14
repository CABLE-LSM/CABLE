SUBROUTINE BLAZE_DRIVER ( NCELLS, BLAZE, SF, casapool,  casaflux, casamet, &
     climate, shootfrac, idoy, curyear, CTLFLAG, POP, veg, call_blaze, call_pop )

  use cable_def_types_mod, only: r_2
  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, DOYSOD2YMDHMS !, Esatf
  USE CASAVARIABLE,        ONLY: casa_pool, casa_flux, casa_met
  USE BLAZE_MOD,           ONLY: RUN_BLAZE, TYPE_TURNOVER, BLAZE_TURNOVER, &
       METB, STR, CWD, LEAF, WOOD, FROOT, TYPE_BLAZE, MLIT, SLIT, CLIT ! , p_surv_OzSavanna
  USE SIMFIRE_MOD,         ONLY: TYPE_SIMFIRE

  USE cable_IO_vars_module, ONLY: landpt, patch
  USE CABLE_DEF_TYPES_MOD, ONLY:  climate_type
  USE POPMODULE,            ONLY: ADJUST_POP_FOR_FIRE
  USE POP_TYPES,            ONLY: POP_TYPE, i4b, dp
  USE cable_def_types_mod, ONLY: veg_parameter_type

  !CLN implement k_tun_XXX here



  IMPLICIT NONE

  TYPE (casa_pool), INTENT(IN)      :: casapool
  TYPE (casa_flux), INTENT(IN)         :: casaflux
  TYPE (casa_met ), INTENT(IN)         :: casamet
  TYPE (climate_type ), INTENT(IN)         :: climate
  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  INTEGER,          INTENT(IN)         :: idoy, CurYear, CTLFLAG,ncells
  REAL, INTENT(IN)    :: shootfrac
  INTEGER, INTENT(IN) :: call_blaze    !=2 if blaze running but not coupled, 3 if coupled
  LOGICAL, INTENT(IN) :: call_pop      !true if cable_user%call_pop = .true.
  CHARACTER(len=3)    :: TOform='old'  !to come in BLAZE TYPE
                                                 

  TYPE(TYPE_TURNOVER)   ,ALLOCATABLE,SAVE :: TO(:,:)
  REAL,   DIMENSION(NCELLS) :: POP_TO

  INTEGER       :: MM, DD, i, np, j, patch_index, p, pidx
  REAL          :: TSTP
  REAL          :: ag_litter_frac, twto, rkill
  REAL          :: CPLANT_g (ncells,3),CPLANT_w (ncells,3)
  REAL          :: CLITTER_g(ncells,3),CLITTER_w(ncells,3)
  LOGICAL       :: EOY

  TYPE(TYPE_BLAZE), INTENT(INOUT)   :: BLAZE
  TYPE(POP_TYPE),             INTENT(INOUT) :: POP
  INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles
  TYPE (TYPE_SIMFIRE) :: SF

  ! INITIALISATION ============================================================

  IF ( BLAZE%BURNMODE .EQ. 0 ) RETURN

  BLAZE%shootfrac = shootfrac
  IF (.not.allocated(TO)) ALLOCATE(TO(BLAZE%NCELLS,7))
  IF ( idoy .EQ. 1) THEN
     DO i = 1, BLAZE%NCELLS
        IF ( ABS(BLAZE%LAT(i)) .GT. 50.) THEN
           BLAZE%K_TUNE_LITTER(i) = BLAZE%K_LITTER_BOREAL
        ELSE IF ( ABS(BLAZE%LAT(i)) .GT. 30.) THEN
           BLAZE%K_TUNE_LITTER(i) = BLAZE%K_LITTER_TEMPERATE
        ELSE IF ( SF%BIOME(i) .EQ. 6 ) THEN ! Savanna
           BLAZE%K_TUNE_LITTER(i) = BLAZE%K_LITTER_SAVANNA
        ELSE ! Tropics
           BLAZE%K_TUNE_LITTER(i) = BLAZE%K_LITTER_TROPICS
        END IF
     END DO
  END IF
  
!!$  IF (CALL1) THEN
!!$     StartYear = CurYear
!!$     CALL1 = .FALSE.
!!$  ENDIF

  BLAZE%CPLANT_g  = 0.
  CLITTER_g = 0.
  BLAZE%CPLANT_w  = 0.
  CLITTER_w = 0.

  !CVH
  ! BLAZE _g and _w pools are  summed over tiles in each gridcell, weighted by patchfrac
  DO i = 1, BLAZE%NCELLS
     DO p = 1, landpt(i)%nap  ! loop over number of active patches
        patch_index = landpt(i)%cstart + p - 1 ! patch index in CABLE vector
        DO j = 1, 3
           IF ( casamet%lnonwood(patch_index) == 1 ) THEN ! Here non-wood
              BLAZE%CPLANT_g (i,j) = BLAZE%CPLANT_g (i,j) + &
                   real(casapool%cplant(patch_index,j)*patch(patch_index)%frac)
              CLITTER_g(i,j) = CLITTER_g(i,j) + &
                   real(casapool%clitter(patch_index,j)*patch(patch_index)%frac)
           ELSEIF ( casamet%lnonwood(patch_index) == 0 ) THEN ! Here woody patches
              BLAZE%CPLANT_w (i,j) = BLAZE%CPLANT_w (i,j) + &
                   real(casapool%cplant(patch_index,j)*patch(patch_index)%frac)
              CLITTER_w(i,j) = CLITTER_w(i,j) + &
                   real(casapool%clitter(patch_index,j)*patch(patch_index)%frac)
          ENDIF
        END DO
     ENDDO
  ENDDO

  BLAZE%CPLANT_w (:,WOOD) =  BLAZE%CPLANT_w (:,WOOD) * BLAZE%shootfrac(:)

  !above ground litter fraction - applies to all PFTs and LUs
  ag_litter_frac = 0.4
  BLAZE%AGLit_w(:, CWD) = CLITTER_w(:, CWD) * BLAZE%shootfrac(:)
  BLAZE%AGLit_w(:,METB) = CLITTER_w(:,METB) * ag_litter_frac
  BLAZE%AGLit_w(:, STR) = CLITTER_w(:, STR) * ag_litter_frac 

  BLAZE%AGLit_g(:, CWD) = 0.
  BLAZE%AGLit_g(:, STR) = CLITTER_g(:, STR) * ag_litter_frac
  BLAZE%AGLit_g(:,METB) = CLITTER_g(:,METB) * ag_litter_frac

  DO i = 1, 3
     BLAZE%BGLit_g(:,i) = CLITTER_g(:,i) - BLAZE%AGLit_g(:,i)
     BLAZE%BGLit_w(:,i) = CLITTER_w(:,i) - BLAZE%AGLit_w(:,i)
  END DO

  CPLANT_g = BLAZE%CPLANT_g
  CPLANT_w = BLAZE%CPLANT_w

  CALL DOYSOD2YMDHMS( CurYear, idoy, 1, MM, DD )

  IF ( idoy .EQ. 366 .OR. ( idoy .EQ. 365 .AND. .NOT. is_leapyear(CurYear)) ) THEN
     EOY = .TRUE.
  ELSE
     EOY = .FALSE.
  END IF

  IF ( IS_LEAPYEAR(CurYear) ) THEN
     TSTP = 1./366.
  ELSE
     TSTP = 1./365.
  END IF

  BLAZE%time =  BLAZE%time + 86400.0

  ! MAIN  ============================================================


  ! BLAZE used to compute ALL or FLI_ONLY
  IF ( BLAZE%BURNMODE .EQ. 1 .OR.  CTLFLAG .EQ. 1 ) THEN

     ! BURNMode 1: BLAZE computes mortality| BURNMode 2: BLAZE computes FLI only 
     !
     ! Ignition 0: GFED derived BA | Ignition 1: SIMFIRE Simulated BA
     ! Ignition is set in BLAZE at the moment.
     ! 
     ! CTLFLAG  -1: Compute Turnovers (using POP-TO) | 1: Compute FLI 

     CALL RUN_BLAZE( BLAZE, SF, CPLANT_g, CPLANT_w, tstp, CurYear, idoy, TO, climate)

     ! compute c-pool turnovers after POP has provided biomass TO 
  ELSE IF ( CTLFLAG .EQ. -1 ) THEN
     !     IF ( .NOT. PRESENT(POP_TO) ) STOP "Provide POP_TO to blaze_casa.f90!"
     DO np = 1, NCELLS
        CALL BLAZE_TURNOVER( BLAZE%AB(np), CPLANT_g(np,:), CPLANT_w(np,:), BLAZE%AGLit_g(np,:), &
             BLAZE%AGLit_w(np,:), BLAZE%BGLit_g(np,:), BLAZE%BGLit_w(np,:),BLAZE%shootfrac(np),TO(np,:), &
             BLAZE%FLUXES(np,:), BLAZE%BURNMODE, BLAZE%IAM,POP_TO(np) )
     END DO
  ELSE
     STOP "Wrong MODE in blaze_driver.f90!"
  ENDIF

  !CVH
  ! set casa fire turnover rates and partitioning of fire losses here!  
  casaflux%kplant_fire   = 0.0_r_2
  casaflux%klitter_fire  = 0.0_r_2
  casaflux%fromPtoL_fire = 0.0_r_2

  ! set burned area and fire line intensity in veg%disturbance_intensity for use in
  ! calculating tree mortality in POP.
  DO i = 1, BLAZE%NCELLS
     DO p = 1, landpt(i)%nap  ! loop over number of active patches
        patch_index = landpt(i)%cstart + p - 1 ! patch index in CABLE vector
        veg%disturbance_intensity(patch_index,1) = BLAZE%AB(i)  ! needed for ADJUST_POP_FOR_FIRE
        veg%disturbance_intensity(patch_index,2) = BLAZE%FLI(i) ! needed for ADJUST_POP_FOR_FIRE
     ENDDO
  ENDDO


  ! get tree mortality rate by applying above height-class mortality to
  ! POP cohorts and then interpolating fire_mortality between POP patches
  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))
  Iw = POP%Iwood

  !do not adjust POP if blaze not coupled to CASA-POP (i.e. call_blaze/=3
  IF (call_blaze==3) THEN
   CALL ADJUST_POP_FOR_FIRE(pop,int(veg%disturbance_interval(Iw,:), i4b), &
        veg%disturbance_intensity(Iw,1), veg%disturbance_intensity(Iw,2)  )
  ENDIF

  ! Apply turn-overs to biomass killed by fire in POP and evaluate fluxes for this call
  BLAZE%FLUXES = 0.0

!   ! pop_grid index
!   pidx=1
!   DO i = 1, BLAZE%NCELLS
!      DO p = 1, landpt(i)%nap  ! loop over number of active patches
!         patch_index = landpt(i)%cstart + p - 1 ! patch index in CABLE vector
 
!         IF ( casamet%lnonwood(patch_index) == 1 ) THEN ! Here non-wood
! !CLN wrong for kplant. Take out cmass here
!            !CLNcasaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i), r_2) &
!            !CLN     * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2)
!            !CLNcasaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i), r_2) &
!            !CLN     * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)

!          !do not set k-terms if blaze not coupled to CASA (i.e. call_blaze/=3)
!          IF (call_blaze==3) THEN
!            casaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i), r_2)
!            casaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i), r_2)
!            casaflux%kplant_fire(patch_index,WOOD)  = 0.0_r_2

!            casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * ag_litter_frac, r_2)
!            casaflux%klitter_fire(patch_index,STR)  = real(BLAZE%AB(i) * ag_litter_frac, r_2)
!            casaflux%klitter_fire(patch_index,CWD)  = 0.0_r_2

!            !for non-woody patches all combustion goes to ATM
!            casaflux%fromPtoL_fire(patch_index,METB,LEAF)  = 0.0_r_2
!            casaflux%fromPtoL_fire(patch_index,METB,FROOT) = 0.0_r_2
!            casaflux%fromPtoL_fire(patch_index,METB,WOOD)  = 0.0_r_2

!            casaflux%fromPtoL_fire(patch_index,STR,LEAF)  = 0.0_r_2
!            casaflux%fromPtoL_fire(patch_index,STR,FROOT) = 0.0_r_2
!            casaflux%fromPtoL_fire(patch_index,STR,WOOD)  = 0.0_r_2

!            casaflux%fromPtoL_fire(patch_index,CWD,LEAF)  = 0.0_r_2
!            casaflux%fromPtoL_fire(patch_index,CWD,FROOT) = 0.0_r_2
!            casaflux%fromPtoL_fire(patch_index,CWD,WOOD)  = 0.0_r_2
!         ENDIF

!            ! BLAZE fluxes - total grid cell plant/litter fluxes due to fire (gC/day)
!            BLAZE%FLUXES(i,11) = BLAZE%FLUXES(i,11) + casaflux%kplant_fire(patch_index,LEAF ) &
!                 * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i,12) = BLAZE%FLUXES(i,12) + casaflux%kplant_fire(patch_index,FROOT) &
!                 * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)

!            BLAZE%FLUXES(i,13) = BLAZE%FLUXES(i,13) + casaflux%klitter_fire(patch_index,METB) &
!                 * real(casapool%clitter(patch_index,METB)*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i,14) = BLAZE%FLUXES(i,14) + casaflux%klitter_fire(patch_index,STR ) &
!                 * real(casapool%clitter(patch_index,STR )*patch(patch_index)%frac, r_2)
           
!         ELSEIF ( casamet%lnonwood(patch_index) == 0 ) THEN ! Here woody patches

!            !CLN increment iwood for pop_grid
           
!            rkill=POP%pop_grid(pidx)%rkill
 
!            pidx = pidx + 1
!            ! Check if there is mortality and COMBUST has only computed non-woody TO
!            ! When POP is involved these fluxes need to sum up to 1, assuming that
!            ! all that is not going to ATM or STR will be going to CWD (DEADWOOD)
!            !IF (CALL_POP) THEN
!            TO(i, WOOD)%TO_CWD = 1. - TO(i, WOOD)%TO_ATM - TO(i, WOOD)%TO_STR
!            TO(i, FROOT)%TO_STR   = MAX( TO(i, WOOD)%TO_ATM + TO(i, WOOD)%TO_STR + TO(i, WOOD)%TO_CWD - &
!                 TO(i, FROOT)%TO_ATM , 0. )
!            ! ENDIF

!            ! Total wood turn-over
!            twto = MAX(TO(i, WOOD)%TO_ATM * 0.7 + TO(i, WOOD )%TO_CWD + TO(i, WOOD )%TO_STR,1.e-7)
           
!            !do not set k-terms if blaze not coupled to CASA (i.e. call_blaze/=3)
!          IF (call_blaze==3) THEN
!            !total turnover for plant and litter pools
!            casaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i) * TO(i, LEAF )%TO_ATM, r_2)
!            casaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i) * (1.-rkill) * TO(i, FROOT)%TO_ATM, r_2)
!            casaflux%kplant_fire(patch_index,WOOD ) = real(rkill/twto  * TO(i, WOOD )%TO_ATM * 0.7, r_2)

!            casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * TO(i, MLIT )%TO_ATM &
!                 * ag_litter_frac, r_2)
!            casaflux%klitter_fire(patch_index,STR)  = real(BLAZE%AB(i) * TO(i, SLIT )%TO_ATM &
!                 * ag_litter_frac, r_2)
!            casaflux%klitter_fire(patch_index,CWD)  = real(BLAZE%AB(i) * TO(i, CLIT )%TO_ATM &
!                 * ag_litter_frac, r_2)

!           ! fraction of total turnover related to conversion of plant to litter
!            casaflux%fromPtoL_fire(patch_index,STR,LEAF) = real(BLAZE%AB(i) * TO(i, LEAF )%TO_STR, r_2)
!            casaflux%fromPtoL_fire(patch_index,STR,FROOT)= real(rkill       * TO(i, FROOT)%TO_STR, r_2)
!            casaflux%fromPtoL_fire(patch_index,STR,WOOD) = real(rkill/twto  * TO(i, WOOD )%TO_STR, r_2)
!            casaflux%fromPtoL_fire(patch_index,CWD,WOOD) = real(rkill/twto  * TO(i, WOOD )%TO_CWD, r_2)

!          ENDIF

!            ! BLAZE fluxes - total grid cell plant/litter fluxes due to fire (gC/day)
!            BLAZE%FLUXES(i, 1) = BLAZE%FLUXES(i, 1) + casaflux%kplant_fire(patch_index,LEAF ) &
!                 * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i, 2) = BLAZE%FLUXES(i, 2) + casaflux%kplant_fire(patch_index,FROOT) &
!                 * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i, 3) = BLAZE%FLUXES(i, 3) + casaflux%kplant_fire(patch_index,WOOD ) &
!                 * real(casapool%cplant(patch_index,WOOD )*patch(patch_index)%frac, r_2) 

!            BLAZE%FLUXES(i, 4) = BLAZE%FLUXES(i, 4) + casaflux%klitter_fire(patch_index,METB) &
!                 * real(casapool%clitter(patch_index,METB)*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i, 5) = BLAZE%FLUXES(i, 5) + casaflux%klitter_fire(patch_index,STR ) &
!                 * real(casapool%clitter(patch_index,STR )*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i, 6) = BLAZE%FLUXES(i, 6) + casaflux%klitter_fire(patch_index,CWD ) &
!                 * real(casapool%clitter(patch_index,CWD )*patch(patch_index)%frac, r_2)
           
!            ! fluxes from plant to litter as a result of fire (gC/day)
!            BLAZE%FLUXES(i, 7) = BLAZE%FLUXES(i, 7) + casaflux%fromPtoL_fire(patch_index,STR,LEAF )  &
!                 * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i, 8) = BLAZE%FLUXES(i, 8) + casaflux%fromPtoL_fire(patch_index,STR,FROOT)  &
!                 * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i, 9) = BLAZE%FLUXES(i, 9) + casaflux%fromPtoL_fire(patch_index,STR,WOOD )  &
!                 * real(casapool%cplant(patch_index,WOOD )*patch(patch_index)%frac, r_2)
!            BLAZE%FLUXES(i,10) = BLAZE%FLUXES(i,10) + casaflux%fromPtoL_fire(patch_index,CWD,WOOD )  &
!                 * real(casapool%cplant(patch_index,WOOD)*patch(patch_index)%frac, r_2)
           
!         ENDIF

!         !IF (MOD(i,50)==0) THEN
!         ! WRITE(*,*) "old calcs", BLAZE%NCELLS
!         ! WRITE(*,*) i, p, pidx, patch_index, BLAZE%AB(i), BLAZE%FLIx(i), BLAZE%w(i)
!         ! WRITE(*,*) "turnovers"
!         ! WRITE(*,*) TO(i,LEAF)%TO_ATM, TO(i,LEAF)%TO_STR
!         ! WRITE(*,*) TO(i,FROOT)%TO_ATM, TO(i,FROOT)%TO_STR
!         ! WRITE(*,*) TO(i, WOOD)%TO_ATM,  TO(i, WOOD)%TO_STR, TO(i, WOOD)%TO_CWD
!         ! WRITE(*,*) "krates normal"
!         ! WRITE(*,*) (casaflux%kplant(patch_index,j), j=1,3)
!         ! WRITE(*,*) (casaflux%klitter(patch_index,j), j=1,3)
!         ! WRITE(*,*) "krates fires"
!         ! WRITE(*,*) (casaflux%kplant_fire(patch_index,j), j=1,3)
!         ! WRITE(*,*) (casaflux%klitter_fire(patch_index,j), j=1,3)
!         ! DO j=1,3
!         !    WRITE(*,*) (casaflux%fromPtoL_fire(patch_index,j,MM), MM=1,3)
!         ! ENDDO
!         ! WRITE(*,*) "fluxes"
!         ! WRITE(*,*) (BLAZE%FLUXES(i,j),j=1,10)
!         ! WRITE(*,*) " "
!       !ENDIF

!       ENDDO ! number of active patches

!   ENDDO ! number of grid cells

!   !07-08-2025 INH - Large scale rewrite of the above which also 
!   !   i) corrects the turnovers passed back to CASA to match expectations there
!   !      including factor of 1/86400.0 assuming daily call to BLAZE.
!   !   ii) allows for BLAZE/SIMFIRE to be used without POP demography active 
!   !       (note POP grid structure is still required)
!   !   iii) provides for an alternate formulation of turonovers

!   ! get tree mortality rate by applying above height-class mortality to
!   ! POP cohorts and then interpolating fire_mortality between POP patches
!   ! do not adjust POP if blaze not coupled to CASA-POP (i.e. call_blaze/=3)
!    if ( (CALL_POP) .and. (call_blaze==3)) then 
!        if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))
!        Iw = POP%Iwood

!        CALL ADJUST_POP_FOR_FIRE(pop,int(veg%disturbance_interval(Iw,:), i4b), &
!             veg%disturbance_intensity(Iw,1), veg%disturbance_intensity(Iw,2)  )
!    ENDIF

!    ! Evaluate turn-overs to biomass (POP dependent) and evaluate fluxes for this call
    BLAZE%FLUXES = 0.0

   !loop over cells
   pidx=1                !counter for pop_grid index
   DO i = 1, BLAZE%NCELLS
      DO p = 1, landpt(i)%nap  ! loop over number of active patches
        patch_index = landpt(i)%cstart + p - 1 ! patch index in CABLE vector
 
        IF ( casamet%lnonwood(patch_index) == 1 ) THEN ! Here non-wood

           !do not set k-terms if blaze not coupled to CASA (i.e. call_blaze/=3)
           IF (call_blaze==3) THEN
              casaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i), r_2)
              casaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i), r_2)
              !casaflux%kplant_fire(patch_index,WOOD)  = 0.0_r_2

              casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * ag_litter_frac, r_2)
              casaflux%klitter_fire(patch_index,STR)  = real(BLAZE%AB(i) * ag_litter_frac, r_2)
              !casaflux%klitter_fire(patch_index,CWD)  = 0.0_r_2

              !for non-woody patches all combustion goes to ATM so no fraction to litter
              !leave casaflux%fromPtoL_fire = 0.0_r_2 for this patch_index)
           ENDIF

           ! BLAZE fluxes - total loss of plant/litter fluxes due to fire activity (all to atm)
           ! units are gC/day
           BLAZE%FLUXES(i,11) = BLAZE%FLUXES(i,11) + casaflux%kplant_fire(patch_index,LEAF ) &
                * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2) 
           BLAZE%FLUXES(i,12) = BLAZE%FLUXES(i,12) + casaflux%kplant_fire(patch_index,FROOT) &
                * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)

           BLAZE%FLUXES(i,13) = BLAZE%FLUXES(i,13) + casaflux%klitter_fire(patch_index,METB) &
                * real(casapool%clitter(patch_index,METB)*patch(patch_index)%frac, r_2)
           BLAZE%FLUXES(i,14) = BLAZE%FLUXES(i,14) + casaflux%klitter_fire(patch_index,STR ) &
                * real(casapool%clitter(patch_index,STR )*patch(patch_index)%frac, r_2)

        ELSEIF( casamet%lnonwood(patch_index) == 0 ) THEN ! Here woody patches
           
           !rkill set to sum of Surawski turnovers if CALL_POP is .false.
           IF (CALL_POP) THEN
              rkill = POP%pop_grid(pidx)%rkill
           ELSE
              rkill = BLAZE%AB(i)*MAX( TO(i, WOOD)%TO_ATM + TO(i, WOOD)%TO_STR + TO(i, WOOD)%TO_CWD,0.0)
           ENDIF
           pidx = pidx + 1
     
           !original formulation for woody turnovers - requires POP%rkill
           !but modified coupling back to CASA
           IF ( (TOform=='old') .and. (CALL_POP) ) THEN

              ! update TO for POP mortality and COMBUST has only computed non-woody TO
              ! when POP is involved these TO need to sum up to 1, assuming that
              ! all that is not going to ATM or STR will be going to CWD (DEADWOOD)
              TO(i, WOOD)%TO_CWD = 1. - TO(i, WOOD)%TO_ATM - TO(i, WOOD)%TO_STR
              TO(i, FROOT)%TO_STR   = MAX( TO(i, WOOD)%TO_ATM + TO(i, WOOD)%TO_STR + TO(i, WOOD)%TO_CWD - &
                    TO(i, FROOT)%TO_ATM, 0._r_2 )
         
              ! Total wood turn-over - should be >0 if AB>0
              twto = MAX(TO(i, WOOD)%TO_ATM * 0.7 + TO(i, WOOD )%TO_CWD + TO(i, WOOD )%TO_STR,1.e-7)

              !do not set k-terms if blaze not coupled to CASA (i.e. call_blaze/=3)
              IF ( (call_blaze==3) .and. (BLAZE%AB(i)>0.0) ) THEN
                 !total turnover for plant and litter pools - units (/day)
                 casaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i) * TO(i, LEAF )%TO_ATM, r_2) + &
                                                        real(BLAZE%AB(i) * TO(i, LEAF )%TO_STR, r_2)
                 casaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i) * (1.-rkill) * TO(i, FROOT)%TO_ATM, r_2) + &
                                                        real(rkill       * TO(i, FROOT)%TO_STR, r_2)
                 casaflux%kplant_fire(patch_index,WOOD ) = real(rkill/twto  * TO(i, WOOD )%TO_ATM * 0.7, r_2) + &
                            real(rkill/twto  * TO(i, WOOD )%TO_STR, r_2) + real(rkill/twto  * TO(i, WOOD )%TO_CWD, r_2)

                 casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * TO(i, MLIT )%TO_ATM &
                    * ag_litter_frac, r_2)
                 casaflux%klitter_fire(patch_index,STR)  = real(BLAZE%AB(i) * TO(i, SLIT )%TO_ATM &
                    * ag_litter_frac, r_2)
                 casaflux%klitter_fire(patch_index,CWD)  = real(BLAZE%AB(i) * TO(i, CLIT )%TO_ATM &
                    * ag_litter_frac, r_2)

                 ! fraction of total turnover related to conversion of plant to litter
                 IF (casaflux%kplant_fire(patch_index,LEAF)>0.0_r_2) THEN
                     casaflux%fromPtoL_fire(patch_index,STR,LEAF) = 1.0_r_2 - & 
                             real(BLAZE%AB(i) * TO(i, LEAF )%TO_ATM, r_2) / casaflux%kplant_fire(patch_index,LEAF)
                 END if
                 IF (casaflux%kplant_fire(patch_index,FROOT)>0.0_r_2) THEN
                     casaflux%fromPtoL_fire(patch_index,STR,FROOT)= 1.0_r_2 - & 
                              real(BLAZE%AB(i) * TO(i, FROOT)%TO_ATM, r_2) / casaflux%kplant_fire(patch_index,FROOT)
                 END IF
                 IF (casaflux%kplant_fire(patch_index,WOOD)>0.0_r_2) THEN
                     casaflux%fromPtoL_fire(patch_index,STR,WOOD) = real(rkill/twto  * TO(i, WOOD )%TO_STR, r_2) / &
                                                            casaflux%kplant_fire(patch_index,WOOD)
                     casaflux%fromPtoL_fire(patch_index,CWD,WOOD) = real(rkill/twto  * TO(i, WOOD )%TO_CWD, r_2) / &
                                                            casaflux%kplant_fire(patch_index,WOOD)
                 END IF
              END IF  !blaze==3 

           !new formulation - rkill set to total wood turnover above (CALL_POP dependent)
           ELSEIF ( (TOform=='new') ) THEN 

              !do not set k-terms if blaze not coupled to CASA (i.e. call_blaze/=3)
              !FLIX gt. 0 precludes case with BA>0 and TO=0 (which should be caught by AVAIL_FUEL as well)
              IF ( (call_blaze==3) .and. (BLAZE%AB(i)>0.0) .and. (BLAZE%FLIx(i) .gt. 0)) THEN
                 ! Total wood turn-over
                 twto = MAX(TO(i, WOOD)%TO_ATM + TO(i, WOOD )%TO_CWD + TO(i, WOOD )%TO_STR,1.e-7)

                 !total turnover for plant and litter pools (units /day)
                 casaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i) * TO(i, LEAF )%TO_ATM, r_2) + &
                                                        real(BLAZE%AB(i) * TO(i, LEAF )%TO_STR, r_2)
                 casaflux%kplant_fire(patch_index,LEAF)  = MAX(casaflux%kplant_fire(patch_index,LEAF),rkill)
                 casaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i) * TO(i, FROOT)%TO_ATM, r_2) + &
                                                        real(BLAZE%AB(i) * TO(i, FROOT)%TO_STR, r_2)
                 casaflux%kplant_fire(patch_index,FROOT) = MAX(casaflux%kplant_fire(patch_index,FROOT),rkill)
                 casaflux%kplant_fire(patch_index,WOOD ) = real(rkill,r_2)

                 casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * TO(i, MLIT )%TO_ATM &
                      * ag_litter_frac, r_2)
                 casaflux%klitter_fire(patch_index,STR)  = real(BLAZE%AB(i) * TO(i, SLIT )%TO_ATM &
                      * ag_litter_frac, r_2)
                 casaflux%klitter_fire(patch_index,CWD)  = real(BLAZE%AB(i) * TO(i, CLIT )%TO_ATM &
                      * ag_litter_frac, r_2)

                 ! fraction of total turnover related to conversion of plant to litter
                 IF (casaflux%kplant_fire(patch_index,LEAF)>0.0_r_2) THEN
                     casaflux%fromPtoL_fire(patch_index,STR,LEAF) = 1.0_r_2 -   &
                           real(BLAZE%AB(i) * TO(i, LEAF )%TO_ATM , r_2)/ casaflux%kplant_fire(patch_index,LEAF)
                 END if
                 IF (casaflux%kplant_fire(patch_index,FROOT)>0.0_r_2) THEN
                     casaflux%fromPtoL_fire(patch_index,STR,FROOT)= 1.0_r_2 -   & 
                           real(BLAZE%AB(i) * TO(i, FROOT)%TO_ATM, r_2) / casaflux%kplant_fire(patch_index,FROOT)
                 END IF
                 IF (twto>0.0) THEN
                   casaflux%fromPtoL_fire(patch_index,STR,WOOD) = real(TO(i, WOOD )%TO_STR / twto, r_2) 
                   casaflux%fromPtoL_fire(patch_index,CWD,WOOD) = real(TO(i, WOOD )%TO_CWD / twto, r_2)
                 END IF
              ENDIF 
           
           ELSE  !TOform is 'new', CALL_POP false - doesn't work
              WRITE(*,*) "BLAZE DRIVER: incorrect combination of CALL_POP and TOform"
              STOP
                 
           ENDIF  !TOform loop

           ! BLAZE fluxes - grid cell plant to atm due to fire - units are gC/day
           BLAZE%FLUXES(i, 1) = BLAZE%FLUXES(i, 1) +                                                        &
                (1.0-casaflux%fromPtoL_fire(patch_index,STR,LEAF ))*casaflux%kplant_fire(patch_index,LEAF ) &
                * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2)
           BLAZE%FLUXES(i, 2) = BLAZE%FLUXES(i, 2) +                                                        &
                (1.0-casaflux%fromPtoL_fire(patch_index,STR,FROOT))*casaflux%kplant_fire(patch_index,FROOT) &
                * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)
           BLAZE%FLUXES(i, 3) = BLAZE%FLUXES(i, 3) +  (1.0-casaflux%fromPtoL_fire(patch_index,STR,WOOD) -   &
                casaflux%fromPtoL_fire(patch_index,CWD,WOOD)) *casaflux%kplant_fire(patch_index,WOOD )      &
                * real(casapool%cplant(patch_index,WOOD )*patch(patch_index)%frac, r_2)

   
            ! BLAZE fluxes - grid cell litter to atm due to fire - units are gC/day
            BLAZE%FLUXES(i, 4) = BLAZE%FLUXES(i, 4) + casaflux%klitter_fire(patch_index,METB) &
              * real(casapool%clitter(patch_index,METB)*patch(patch_index)%frac, r_2)
            BLAZE%FLUXES(i, 5) = BLAZE%FLUXES(i, 5) + casaflux%klitter_fire(patch_index,STR ) &
                * real(casapool%clitter(patch_index,STR )*patch(patch_index)%frac, r_2)
            BLAZE%FLUXES(i, 6) = BLAZE%FLUXES(i, 6) + casaflux%klitter_fire(patch_index,CWD ) &
                * real(casapool%clitter(patch_index,CWD )*patch(patch_index)%frac, r_2)
           
            ! fluxes from plant to litter as a result of fire - units are gC/day
            BLAZE%FLUXES(i, 7) = BLAZE%FLUXES(i, 7) +                                                  &
                casaflux%fromPtoL_fire(patch_index,STR,LEAF )*casaflux%kplant_fire(patch_index,LEAF )  &
                * real(casapool%cplant(patch_index,LEAF )*patch(patch_index)%frac, r_2)
            BLAZE%FLUXES(i, 8) = BLAZE%FLUXES(i, 8) +                                                  &
                casaflux%fromPtoL_fire(patch_index,STR,FROOT)*casaflux%kplant_fire(patch_index,FROOT)  &
                * real(casapool%cplant(patch_index,FROOT)*patch(patch_index)%frac, r_2)
            BLAZE%FLUXES(i, 9) = BLAZE%FLUXES(i, 9) +                                                  &
                casaflux%fromPtoL_fire(patch_index,STR,WOOD )*casaflux%kplant_fire(patch_index,WOOD )  &
                * real(casapool%cplant(patch_index,WOOD )*patch(patch_index)%frac, r_2)
            BLAZE%FLUXES(i,10) = BLAZE%FLUXES(i,10) +                                                  &
                casaflux%fromPtoL_fire(patch_index,CWD,WOOD )*casaflux%kplant_fire(patch_index,WOOD )  &
                * real(casapool%cplant(patch_index,WOOD)*patch(patch_index)%frac, r_2)

        ENDIF  ! casamet%lnonwood = 0/1

        !IF (MOD(i,50)==0) THEN
        ! WRITE(*,*) "new calcs"
        ! WRITE(*,*) i, p, pidx, patch_index, BLAZE%AB(i), BLAZE%FLIx(i), BLAZE%w(i)
        ! WRITE(*,*) "turnovers"
        ! WRITE(*,*) TO(i,LEAF)%TO_ATM, TO(i,LEAF)%TO_STR
        ! WRITE(*,*) TO(i,FROOT)%TO_ATM, TO(i,FROOT)%TO_STR
        ! WRITE(*,*) TO(i, WOOD)%TO_ATM,  TO(i, WOOD)%TO_STR, TO(i, WOOD)%TO_CWD
        ! WRITE(*,*) "krates normal"
        ! WRITE(*,*) (casaflux%kplant(patch_index,j), j=1,3)
        ! WRITE(*,*) (casaflux%klitter(patch_index,j), j=1,3)
        ! WRITE(*,*) "krates fires"
        ! WRITE(*,*) (casaflux%kplant_fire(patch_index,j), j=1,3)
        ! WRITE(*,*) (casaflux%klitter_fire(patch_index,j), j=1,3)
        ! DO j=1,3
        !    WRITE(*,*) (casaflux%fromPtoL_fire(patch_index,j,MM), MM=1,3)
        ! ENDDO
        ! WRITE(*,*) "fluxes"
        ! WRITE(*,*) (BLAZE%FLUXES(i,j),j=1,10)
        ! WRITE(*,*) " "
        ! WRITE(*,*) " "
      !ENDIF

     END DO    ! number of active patches 
  END DO       ! number of grid cells

  !STOP

END SUBROUTINE BLAZE_DRIVER
