SUBROUTINE BLAZE_DRIVER ( NCELLS, BLAZE, SF, casapool,  casaflux, casamet, &
     climate, shootfrac, idoy, curyear, CTLFLAG )

  use cable_def_types_mod, only: r_2
  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, DOYSOD2YMDHMS, kbl_user_switches !, Esatf
  USE CASAVARIABLE,        ONLY: casa_pool, casa_flux, casa_met
  USE BLAZE_MOD,           ONLY: RUN_BLAZE, TYPE_TURNOVER, BLAZE_TURNOVER, NTO, &
       METB, STR, CWD, LEAF, WOOD, FROOT, TYPE_BLAZE, MLIT, SLIT, CLIT, p_surv_OzSavanna
  USE SIMFIRE_MOD,         ONLY: TYPE_SIMFIRE

  USE cable_IO_vars_module, ONLY: landpt, patch
  USE CABLE_DEF_TYPES_MOD, ONLY:  climate_type
  IMPLICIT NONE
 
  TYPE(kbl_user_switches) :: cable_user
  TYPE (casa_pool), INTENT(IN)      :: casapool
  TYPE (casa_flux), INTENT(IN)         :: casaflux
  TYPE (casa_met ), INTENT(IN)         :: casamet
  TYPE (climate_type ), INTENT(IN)         :: climate
  INTEGER,          INTENT(IN)         :: idoy, CurYear, CTLFLAG,ncells
  REAL, INTENT(IN)    :: shootfrac

  INTEGER, PARAMETER        :: NPOOLS = 3

  TYPE(TYPE_TURNOVER)   ,ALLOCATABLE,SAVE :: TO(:,:)
   REAL,   DIMENSION(NCELLS) :: AGL_wo1,AGL_wo2,AGL_wo3
  REAL,   DIMENSION(NCELLS) :: POP_TO,POP_CWD,POP_STR

  INTEGER,DIMENSION(NCELLS) :: modis_igbp
  INTEGER,DIMENSION(NCELLS) :: t1, t2
  REAL,   DIMENSION(NCELLS) :: AB, relhum, U10, FLI, DFLI, FFDI !CRM , popd, mnest
  INTEGER       :: MM, DD, i, np, j, patch_index, p, nbins, nh
  REAL          :: TSTP, C_CHKSUM
  REAL          :: ag_lit, tot_lit, bg_lit, ag_litter_frac
  REAL          :: CPLANT_g (ncells,3),CPLANT_w (ncells,3)
  REAL          :: CLITTER_g(ncells,3),CLITTER_w(ncells,3)
  LOGICAL       :: EOY
  LOGICAL, SAVE :: CALL1 = .TRUE.
  LOGICAL, SAVE :: YEAR1 = .TRUE.

  TYPE(TYPE_BLAZE), INTENT(INOUT)   :: BLAZE
  TYPE (TYPE_SIMFIRE) :: SF

  REAL    :: C_BIOMASS, C_FIRE, hgt
  LOGICAL,PARAMETER :: CLOSURE_TEST = .FALSE.
  REAL,PARAMETER:: BIN_POWER=1.2 ! bins have muscles

  ! INITIALISATION ============================================================
 
  IF ( BLAZE%BURNMODE .EQ. 0 ) RETURN

  BLAZE%shootfrac = shootfrac
  IF (.not.allocated(TO)) ALLOCATE(TO(BLAZE%NCELLS,7))

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
                   real(casapool%clitter(patch_index,i)*patch(patch_index)%frac)
           ELSEIF ( casamet%lnonwood(patch_index) == 0 ) THEN ! Here woody patches
              BLAZE%CPLANT_w (i,j) = BLAZE%CPLANT_w (i,j) + &
                   real(casapool%cplant(patch_index,i)*patch(patch_index)%frac)
              CLITTER_w(i,j) = CLITTER_w(i,j) + &
                   real(casapool%clitter(patch_index,i)*patch(patch_index)%frac)
           ENDIF
        END DO
     ENDDO
  ENDDO
 

  ! set heights at which tree mortality is calculated
  nbins = 30
  DO i=1,nbins
     casaflux%fire_mortality_vs_height(:,i,1) = real(BIN_POWER**i, r_2)
  ENDDO


  ag_litter_frac = 0.4
  BLAZE%AGLit_w(:, CWD) = CLITTER_w(:, CWD) * BLAZE%shootfrac(:)
  BLAZE%AGLit_w(:,METB) = CLITTER_w(:,METB) * ag_litter_frac
  BLAZE%AGLit_w(:,STR) = CLITTER_w(:,STR) * ag_litter_frac
  
  BLAZE%AGLit_g(:,CWD)  = 0.
  BLAZE%AGLit_g(:,STR ) = CLITTER_g(:,STR) * ag_litter_frac
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

  BLAZE%time =  BLAZE%time + 86400

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
             BLAZE%FLUXES(np,:), BLAZE%BURNMODE, POP_TO(np) )
     END DO
  ELSE
     STOP "Wrong MODE in blaze_driver.f90!"
  ENDIF

!CVH
! set casa fire turnover rates and partitioning of fire losses here!  
  casaflux%kplant_fire(:,LEAF) = 0.0_r_2
  casaflux%kplant_fire(:,FROOT) = 0.0_r_2
  casaflux%klitter_fire = 0.0_r_2
  casaflux%fromPtoL_fire = 0.0_r_2
  
  DO i = 1, BLAZE%NCELLS
     DO p = 1, landpt(i)%nap  ! loop over number of active patches
        patch_index = landpt(i)%cstart + p - 1 ! patch index in CABLE vector
        
        IF ( casamet%lnonwood(patch_index) == 1 ) THEN ! Here non-wood
           
           casaflux%kplant_fire(patch_index,LEAF)  = real(BLAZE%AB(i), r_2)
           casaflux%kplant_fire(patch_index,FROOT) =  real(BLAZE%AB(i), r_2)
           casaflux%kplant_fire(patch_index,WOOD)  = 0.0_r_2
           
           casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * ag_litter_frac, r_2)
           casaflux%klitter_fire(patch_index,STR)  =  real(BLAZE%AB(i) * ag_litter_frac, r_2)
           casaflux%klitter_fire(patch_index,CWD)  = 0.0_r_2

           casaflux%fromPtoL_fire(patch_index,METB,LEAF)  = 0.0_r_2
           casaflux%fromPtoL_fire(patch_index,METB,FROOT) = 0.0_r_2
           casaflux%fromPtoL_fire(patch_index,METB,WOOD)  = 0.0_r_2

           casaflux%fromPtoL_fire(patch_index,STR,LEAF)  = 0.0_r_2
           casaflux%fromPtoL_fire(patch_index,STR,FROOT) = 0.0_r_2
           casaflux%fromPtoL_fire(patch_index,STR,WOOD)  = 0.0_r_2

           casaflux%fromPtoL_fire(patch_index,CWD,LEAF)  = 0.0_r_2
           casaflux%fromPtoL_fire(patch_index,CWD,FROOT) = 0.0_r_2
           casaflux%fromPtoL_fire(patch_index,CWD,WOOD)  = 0.0_r_2
           
        ELSEIF ( casamet%lnonwood(patch_index) == 0 ) THEN ! Here woody patches
           casaflux%kplant_fire(patch_index,LEAF) =  real(BLAZE%AB(i) * &
                                    (TO(i, LEAF )%TO_ATM + TO(i, LEAF )%TO_STR), r_2)
           casaflux%kplant_fire(patch_index,FROOT) = real(BLAZE%AB(i) * &
                (TO(i, LEAF )%TO_ATM + TO(i, LEAF )%TO_STR), r_2)
           
           !casaflux%kplant_fire(patch_index,WOOD) = 0.0 ! to be set after call to POP
           
           casaflux%klitter_fire(patch_index,METB) = real(BLAZE%AB(i) * &
                                                     (TO(i, MLIT )%TO_ATM) * ag_litter_frac, r_2)
           casaflux%klitter_fire(patch_index,STR) =  real(BLAZE%AB(i) * &
                                                     (TO(i, SLIT )%TO_ATM) * ag_litter_frac, r_2)
           casaflux%klitter_fire(patch_index,CWD) =  real(BLAZE%AB(i) * &
                (TO(i, CLIT )%TO_ATM) * ag_litter_frac, r_2)

           casaflux%fromPtoL_fire(patch_index,STR,LEAF) = real(TO(i, LEAF )%TO_STR/ &
                MAX((TO(i, LEAF )%TO_STR + TO(i, LEAF )%TO_ATM),1.e-5), r_2)

           casaflux%fromPtoL_fire(patch_index,STR,FROOT) = real(TO(i, FROOT )%TO_STR/ &
                MAX((TO(i, FROOT )%TO_STR + TO(i, FROOT )%TO_ATM),1.e-5), r_2)

           casaflux%fromPtoL_fire(patch_index,CWD,WOOD) = real(TO(i, WOOD )%TO_CWD/ &
                MAX((TO(i, WOOD )%TO_CWD + TO(i, WOOD )%TO_ATM +  TO(i, WOOD )%TO_STR ),1.e-5), r_2)

           casaflux%fromPtoL_fire(patch_index,STR,WOOD) = real(TO(i, WOOD )%TO_STR/ &
                MAX((TO(i, WOOD )%TO_CWD + TO(i, WOOD )%TO_ATM + &
                TO(i, WOOD )%TO_STR ),1e-5), r_2)
           
           ! get tree mortality by height class
           DO nh = 1,30
              hgt = real(casaflux%fire_mortality_vs_height(patch_index,nh,1))
              casaflux%fire_mortality_vs_height(patch_index,nh,2) = &
                   (1._r_2 - real(p_surv_OzSavanna(hgt,BLAZE%FLI(i)),r_2)) * &
                   (1._r_2 - casaflux%fire_mortality_vs_height(patch_index,nh,2))
           ENDDO
               
        ENDIF

     ENDDO ! number of active patches

  ENDDO ! number of grid cells

END SUBROUTINE BLAZE_DRIVER
