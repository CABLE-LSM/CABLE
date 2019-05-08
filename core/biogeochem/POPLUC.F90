!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: module for land-use change which interacts with POP demography
! via secondary forest age-distribution, and updates casa stocks according to land-use transitions
!
! Called from: cable_driver or cable_mpimaster
!
! SUBROUTINES
! ZeroPOPLUC(POPLUC)
! execute_luc_event(from_state,to_state,frac_change_grid,g,POPLUC)
! CALCULATE_WEIGHTS(POPLUC, g)
! INCREMENT_AGE(POPLUC,g)
! POPLUCStep(POPLUC,year)
! POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)
! POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)
! POPLUC_Init(POPLUC,LUC_EXPT, casapool, casaflux, casabiome, veg, POP, np)
! POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
!  POPLUC_set_params(POPLUC,LUC_EXPT)
! alloc_POPLUC(POPLUC, arraysize)
! WRITE_LUC_OUTPUT_NC ( POPLUC, ctime, FINAL )
! WRITE_LUC_RESTART_NC ( POPLUC, ctime )
! READ_LUC_RESTART_NC (POPLUC)
! WRITE_LUC_OUTPUT_GRID_NC ( POPLUC, ctime, FINAL )


! History: Vanessa Haverd July 2016

! ==============================================================================
!*********************************************************************************
MODULE POPLUC_CONSTANTS
  USE TYPEdef, ONLY: dp, i4b

  INTEGER(i4b),PARAMETER ::  LENGTH_SECDF_HISTORY = 4000
  INTEGER(i4b),PARAMETER :: AGE_MAX = 1000
  INTEGER(i4b),PARAMETER :: disturbance_interval = 100
  ! N.B. needs to be the same as veg%disturbance_interval
  LOGICAL, PARAMETER :: IFHARVEST=.FALSE.
  INTEGER(i4b), PARAMETER :: ROTATION=70
  INTEGER(i4b), PARAMETER :: nLU=3 ! number of land-use tiles (pf, sf, grass)
  INTEGER(i4b), PARAMETER :: nTrans=4
  ! number of possible gross transition types (ptog, ptos, stog, gtos)

END MODULE POPLUC_CONSTANTS

!*******************************************************************************
MODULE POPLUC_Types
  USE TYPEdef, ONLY: dp, i4b
  USE POPLUC_Constants, ONLY: LENGTH_SECDF_HISTORY, AGE_MAX

  TYPE POPLUC_TYPE
     INTEGER(i4b),POINTER :: it
     INTEGER(i4b),POINTER :: np
     INTEGER(i4b),POINTER :: firstyear
     INTEGER(i4b),POINTER :: thisyear
     INTEGER(i4b), DIMENSION(:),POINTER :: n_event ! number of secondary forest transitions
     REAL(dp), DIMENSION(:),POINTER :: latitude, longitude
     REAL(dp), DIMENSION(:),POINTER :: primf, secdf, grass,       &  ! land cover types
          ptos, ptog, STOP, stog, gtop, gtos,    & ! transitions
          frac_primf, frac_forest
     REAL(dp), DIMENSION(:,:),POINTER ::  freq_age_primary, freq_age_secondary, &
          biomass_age_primary, biomass_age_secondary
     REAL(dp), DIMENSION(:,:),POINTER :: age_history_secdf, area_history_secdf
     REAL(dp), DIMENSION(:,:),POINTER :: FNEP, Clitt, Csoil, Cbiomass
     REAL(dp), DIMENSION(:,:),POINTER :: FHarvest, FClearance, FTransferNet
     REAL(dp), DIMENSION(:,:),POINTER :: FTransferGross
     REAL(dp), DIMENSION(:),POINTER :: pharv, smharv, syharv
     REAL(dp), DIMENSION(:,:),POINTER :: HarvProd, ClearProd
     REAL(dp), DIMENSION(:,:),POINTER :: fracHarvProd, fracClearProd
     REAL(dp), DIMENSION(:,:),POINTER :: HarvProdLoss, ClearProdLoss
     REAL(dp), DIMENSION(:),POINTER :: fracHarvResid, fracHarvSecResid, fracClearResid

  END TYPE POPLUC_TYPE

END MODULE POPLUC_Types
!*******************************************************************************

MODULE POPLUC_Module

  !-------------------------------------------------------------------------------
  ! * This module contains all subroutines for POPLUC calcs at a single time step.
  !-------------------------------------------------------------------------------
  USE TYPEdef, ONLY: sp, i4b
  USE POPLUC_Types
  USE POPLUC_Constants
  USE casavariable, ONLY: casa_pool, casa_balance, casa_flux, casa_biome
  USE POP_Types, ONLY: POP_TYPE
  USE cable_common_module, ONLY: cable_user
  USE cable_IO_vars_module, ONLY: landpt, patch
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE
  USE POPModule, ONLY: pop_init_single

CONTAINS

  !*******************************************************************************
  SUBROUTINE ZeroPOPLUC(POPLUC)
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER:: g,np

    np = popluc%np

    POPLUC%firstyear = 0
    POPLUC%thisyear = 0
    POPLUC%primf = 0
    POPLUC%secdf = 0
    POPLUC%grass = 0
    POPLUC%ptos = 0
    POPLUC%ptog = 0
    POPLUC%stop = 0
    POPLUC%stog = 0
    POPLUC%gtop = 0
    POPLUC%gtos = 0
    POPLUC%frac_forest = 0
    POPLUC%frac_primf = 0
    POPLUC%area_history_secdf = 0
    POPLUC%age_history_secdf = 0
    POPLUC%n_event = 0
    POPLUC%freq_age_secondary = 0
    POPLUC%freq_age_primary = 0
    POPLUC%biomass_age_primary = 0
    POPLUC%biomass_age_secondary = 0
    POPLUC%FNEP = 0
    POPLUC%Clitt = 0
    POPLUC%Csoil = 0
    POPLUC%Cbiomass = 0
    POPLUC%FHarvest = 0
    POPLUC%FClearance = 0
    POPLUC%FTransferNet = 0
    POPLUC%FTransferGross = 0
    POPLUC%pharv = 0
    POPLUC%smharv = 0
    POPLUC%syharv = 0
    POPLUC%HarvProd = 0
    POPLUC%HarvProdLoss = 0
    POPLUC%ClearProd = 0
    POPLUC%ClearProdLoss = 0


  END SUBROUTINE ZeroPOPLUC
  !*******************************************************************************

  SUBROUTINE execute_luc_event(from_state,to_state,frac_change_grid,g,POPLUC)
    ! Execute a transition between land use types (states)
    !  frac_change_grid = fractional change in unit fraction of grid cell this year

    IMPLICIT NONE

    CHARACTER(5), INTENT(IN) :: from_state, to_state
    REAL(dp), INTENT(INOUT):: frac_change_grid
    TYPE(POPLUC_TYPE),INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: g  ! grid cell index
    REAL :: frac_open_grid, remaining
    INTEGER(i4b) :: n  ! position of new element in POPLUC%SecFor array
    INTEGER(i4b) :: i
    frac_open_grid = 1.0 -POPLUC%frac_forest(g)
    n = POPLUC%n_event(g)

    IF (from_state=='PRIMF') THEN

       IF(frac_change_grid.GT.POPLUC%primf(g)) THEN
          !PRINT*, "Warning: requested reduction in primary forest area &
          !     exceeds primary forest area"
          IF (to_state=='SECDF') POPLUC%ptos(g) = POPLUC%primf(g)
          IF (to_state=='C3ANN') POPLUC%ptog(g) = POPLUC%primf(g)
          frac_change_grid = POPLUC%primf(g)
          POPLUC%primf(g) = 0.0
       ELSE
          POPLUC%primf(g) = POPLUC%primf(g) &
               -frac_change_grid
       ENDIF

       IF (to_state=='SECDF') THEN
          ! Transition from primary -> secondary forest(ptos)
          POPLUC%n_event(g) = POPLUC%n_event(g)+1
          n = POPLUC%n_event(g)
          POPLUC%area_history_secdf(g,n) = frac_change_grid
          POPLUC%age_history_secdf(g,n) = 0
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_change_grid

       ENDIF

    ELSEIF (from_state=='SECDF') THEN

       IF (to_state=='PRIMF') THEN
          PRINT*, "Error: cannot create primary forest from secondary forest"
          STOP

       ELSE
          ! Transition from secondary -> non forest (stog)
          ! Assumption: youngest stands cleared first

          remaining = frac_change_grid
          i = 1
          DO WHILE (remaining > 0.0 .AND. i <= age_max )
             IF (POPLUC%freq_age_secondary(g,i).GE.remaining) THEN
                POPLUC%freq_age_secondary(g,i) =POPLUC%freq_age_secondary(g,i) &
                     - remaining
                remaining = 0.0
             ELSE
                remaining = remaining - POPLUC%freq_age_secondary(g,i)
                POPLUC%freq_age_secondary(g,i) = 0.0
                i = i+1
             ENDIF

          ENDDO
          !if (remaining.gt.frac_change_grid) POPLUC%stog(g) = POPLUC%stog(g)-remaining
          IF (remaining.GT.0.0) POPLUC%stog(g) = POPLUC%stog(g)-remaining
       ENDIF

    ELSEIF (to_state=='SECDF') THEN

       POPLUC%n_event(g) = POPLUC%n_event(g)+1
       n = POPLUC%n_event(g)
       ! Transition from non-forest to secondary forest (gtos)

       IF (frac_change_grid.LE.frac_open_grid) THEN
          POPLUC%area_history_secdf(g,n) = frac_change_grid
          POPLUC%age_history_secdf(g,n) = 0
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_change_grid
       ELSE
          POPLUC%area_history_secdf(g,n) = frac_open_grid
          POPLUC%age_history_secdf(g,n) = 0
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_open_grid
          ! gtos to frac_change_grid here!!
          POPLUC%gtos(g) =  frac_open_grid


       ENDIF

    ELSEIF (to_state=='PRIMF') THEN

       PRINT*, "Error: cannot create primary forest from non-forest"
       STOP

    ENDIF


  ENDSUBROUTINE execute_luc_event


  !*******************************************************************************
  SUBROUTINE CALCULATE_WEIGHTS(POPLUC, g)
    ! Calculates weights (fraction of total forest area on grid cell)
    !for primary and secondary forest stands up to specified maximum stand age


    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: g
    REAL(dp), PARAMETER :: EPS= 1.e-9
    INTEGER(i4b) :: age , i, iage
    REAL(dp) :: fac, disturbance_freq

    ! First get relative weights for primary forest
    disturbance_freq=1.0/REAL(disturbance_interval,dp)

    !fac = POPLUC%frac_primf/POPLUC%frac_forest
    fac = 1.0
    DO iage = 1, age_max
       POPLUC%freq_age_primary(g,iage) =  REALExponential(disturbance_freq,REAL(iage-1,dp))
       POPLUC%freq_age_secondary(g,iage) = 0.0
    END DO
    POPLUC%freq_age_primary(g,:) =  POPLUC%freq_age_primary(g,:) &
         / SUM(POPLUC%freq_age_primary(g,:))*fac


    !  Loop through secondary forest stands to transfer weights
    fac = 1.0
    DO i = 1, POPLUC%n_event(g)

       age = POPLUC%age_history_secdf(g,i)
       POPLUC%freq_age_secondary(g,age+1) = POPLUC%area_history_secdf(g,i)/fac
    ENDDO


  END SUBROUTINE CALCULATE_WEIGHTS
  !*******************************************************************************
  SUBROUTINE INCREMENT_AGE(POPLUC,g)

    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b) :: n_event, i, j
    INTEGER(i4b), INTENT(IN) :: g
    REAL(dp):: area, remaining

    n_event =  POPLUC%n_event(g)

    POPLUC%freq_age_secondary(g,2:age_max)=POPLUC%freq_age_secondary(g,1:age_max-1)
    POPLUC%freq_age_secondary(g,1) = 0.0

    ! adjust secondary age distribution for secondary forest harvest area
    IF (POPLUC%smharv(g)+POPLUC%syharv(g).GT.0) THEN

       remaining = POPLUC%smharv(g)+POPLUC%syharv(g)
       i = age_max
       DO WHILE (remaining.GT.0.0)
          IF (POPLUC%freq_age_secondary(g,i) .GT. remaining) THEN

             POPLUC%freq_age_secondary(g,i) = POPLUC%freq_age_secondary(g,i) - remaining
             POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + remaining
             remaining = 0.0
          ELSE
             POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + &
                  POPLUC%freq_age_secondary(g,i)
             remaining = remaining - POPLUC%freq_age_secondary(g,i)
             POPLUC%freq_age_secondary(g,i) = 0.0
             i = i-1;
          END IF
          IF (i.LT.2) THEN
             POPLUC%smharv(g) = POPLUC%smharv(g)+POPLUC%syharv(g) - remaining
             POPLUC%syharv(g) = 0.0
             remaining = 0.0
          ENDIF
       ENDDO
    ENDIF

    IF (IFHARVEST .AND. POPLUC%freq_age_secondary(g,ROTATION+1).GT.0.0) THEN
       POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + &
            POPLUC%freq_age_secondary(g,ROTATION+1)
       POPLUC%freq_age_secondary(g,ROTATION+1) = 0.0
    ENDIF

    ! adjust secondary age distribution for natural disturbance
    i = age_max
    DO i = age_max, 2 , -1
       POPLUC%freq_age_secondary(g,1) =  POPLUC%freq_age_secondary(g,1) +  &
            POPLUC%freq_age_secondary(g,i)/disturbance_interval
       POPLUC%freq_age_secondary(g,i) = POPLUC%freq_age_secondary(g,i)* &
            (1. - 1./disturbance_interval)
    ENDDO


  END SUBROUTINE INCREMENT_AGE
  !*******************************************************************************
  SUBROUTINE POPLUCStep(POPLUC,year)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: year
    INTEGER(i4b) :: g,j

    POPLUC%it = POPLUC%it + 1


    IF (year.LT.POPLUC%thisyear) THEN

       DO g = 1,POPLUC%np
          ! CALL calculate_weights(POPLUC, g)
          CALL increment_age(POPLUC,g)
       ENDDO

    ELSE

       DO g = 1,POPLUC%np


          IF (POPLUC%ptos(g) .GT. 0.0) &
               CALL execute_luc_event('PRIMF','SECDF',POPLUC%ptos(g),g,POPLUC)

          IF (POPLUC%ptog(g) .GT. 0.0) &
               CALL execute_luc_event('PRIMF','C3ANN',POPLUC%ptog(g),g,POPLUC)

          IF (POPLUC%stop(g) .GT.0.0) &
               CALL execute_luc_event('SECDF','PRIMF',POPLUC%stop(g),g,POPLUC)
          IF (POPLUC%stog(g) .GT.0.0) &
               CALL execute_luc_event('SECDF','C3ANN',POPLUC%stog(g),g,POPLUC)

          IF (POPLUC%gtop(g) .GT.0.0) &
               CALL execute_luc_event('C3ANN','PRIMF',POPLUC%gtop(g),g,POPLUC)
          IF (POPLUC%gtos(g) .GT.0.0) &
               CALL execute_luc_event('C3ANN','SECDF',POPLUC%gtos(g),g,POPLUC)

          POPLUC%frac_forest(g) =  POPLUC%primf(g)+ SUM(POPLUC%freq_age_secondary(g,:))
          CALL increment_age(POPLUC,g)

       ENDDO

    ENDIF

  END SUBROUTINE POPLUCStep
  !*******************************************************************************
  SUBROUTINE POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

    !-------------------------------------------------------------------------------
    ! This subroutine transfers LUC-determined age distributions to POP
    !-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    INTEGER:: g, k, j, l
    REAL(dp), DIMENSION(:), ALLOCATABLE:: freq_age

    DO l = 1, POP%np
       pop%pop_grid(l)%freq_age = 0.0
    ENDDO
    ALLOCATE (freq_age(age_max))

    DO g = 1,POPLUC%np   ! loop over POPLUC gricells (same as CABLE grid-cells)
       IF (.NOT.LUC_EXPT%prim_only(g) .AND.SUM(POPLUC%freq_age_secondary(g,:)).GT.1e-12) THEN
          ! check if tile is secondary forest
          j = landpt(g)%cstart + 1       ! cable index for secondary forest tile
          freq_age = POPLUC%freq_age_secondary(g,:)/SUM(POPLUC%freq_age_secondary(g,:))
          DO l = 1, POP%np   ! for each wooded tile
             IF (POP%Iwood(l).EQ.j) THEN
                ! check if cable index in pop structure (Iwood) is the target cable index
                pop%pop_grid(l)%freq_age = freq_age
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    DEALLOCATE (freq_age)

  END SUBROUTINE POPLUC_weights_transfer
  !*******************************************************************************
  SUBROUTINE POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal,casaflux,ktauday)

    !-------------------------------------------------------------------------------
    ! This subroutine redestributes carbon (nitrogen and phosphorous) amongst
    ! pools according to land-use transtions
    !-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    TYPE(POP_TYPE), INTENT(IN) :: POP
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    TYPE (casa_pool),           INTENT(INOUT) :: casapool
    TYPE (casa_balance),        INTENT(INOUT) :: casabal
    TYPE (casa_flux),        INTENT(INOUT) :: casaflux
    INTEGER, INTENT(IN) :: ktauday
    ! number of cable time-steps in a day (for needed for LUC flux output)
    INTEGER:: g, k, j, l, idp, irp, idlu, irlu, ilu
    INTEGER,  PARAMETER :: &
         ptos         =  1, &
         ptog         =  2, &
         stog         =  3, &
         gtos         =  4
    INTEGER,  PARAMETER :: &
         p         =  1, &
         s         =  2, &
         gr         =  3
    REAL(dp) :: dcplant(nLU,nLU,3), dclitter(nLU,nLU,3), dcsoil(nLU,nLU,3)
    REAL(dp) :: dcplant_r(nLU,3), dclitter_r(nLU,3), dcsoil_r(nLU,3)
    REAL(dp) :: dcplant_d(nLU,3), dclitter_d(nLU,3), dcsoil_d(nLU,3)
    REAL(dp) :: dnplant(nLU,nLU,3), dnlitter(nLU,nLU,3), dnsoil(nLU,nLU,3)
    REAL(dp) :: dnplant_r(nLU,3), dnlitter_r(nLU,3), dnsoil_r(nLU,3)
    REAL(dp) :: dnplant_d(nLU,3), dnlitter_d(nLU,3), dnsoil_d(nLU,3)
    REAL(dp) :: dpplant(nLU,nLU,3), dplitter(nLU,nLU,3), dpsoil(nLU,nLU,3)
    REAL(dp) :: dpplant_r(nLU,3), dplitter_r(nLU,3), dpsoil_r(nLU,3)
    REAL(dp) :: dpplant_d(nLU,3), dplitter_d(nLU,3), dpsoil_d(nLU,3)
    REAL(dp) :: dnsoilmin_r(nLU), dclabile_r(nLU)
    REAL(dp) :: dnsoilmin_d(nLU), dclabile_d(nLU)
    REAL(dp) :: dclabile(nLU, nLU), dnsoilmin(nLU,nLU)
    REAL(dp) :: dA_r(nLU), dA_d(nLU), dA(nLU), deltaA, dwood_transfer
    REAL(dp), ALLOCATABLE :: dcHarvCLear(:), dcHarv(:), dcClear(:), dcExpand(:), FHarvClear(:)
    REAL(dp) :: kHarvProd(3), kClearProd(3)


    ! turnover rates for harvest and clearance products (y-1)
    kHarvProd(1) = 1.0/1.0
    kHarvProd(2) = 1.0/10.0
    kHarvProd(3) = 1.0/100.0
    kClearProd = kHarvProd

    ! zero POPLUC fluxes
    popluc%FtransferNet = 0.0
    popluc%FtransferGross = 0.0
    popluc%FHarvest = 0.0
    popluc%FClearance = 0.0

    ! local variable for storing sum of biomass change due to
    !secondary harvest, clearance and expansion, and secondary forest
    !  harvest and clearance fluxes
    ALLOCATE(FHarvClear(POPLUC%np))
    ALLOCATE(dcHarvClear(POPLUC%np))
    ALLOCATE(dcHarv(POPLUC%np))
    ALLOCATE(dcClear(POPLUC%np))
    ALLOCATE(dcExpand(POPLUC%np))
    FHarvClear = 0.0
    dcHarvClear = 0.0
    dcHarv = 0.0
    dcClear = 0.0
    dcExpand = 0.0
    POPLUC%FHarvest = 0.0
    POPLUC%FClearance = 0.0
    casaflux%FluxCtohwp = 0.0
    casaflux%FluxCtoclear = 0.0
    casaflux%CtransferLUC = 0.0


    ! Transfer POP age-dependent biomass to POLUC variables (diagnostic only) and
    ! catastrophic mortality in secondary forest tiles to changes in biomass associated
    ! with secondary forest harvest and clearance and expansion.
    DO g = 1,POPLUC%np  ! loop over CABLE grid-cells
       j = landpt(g)%cstart ! index of primary veg tile in each grid-cell
       DO l = 1, POP%np   ! loop over all POP tiles (== wooded tiles in CABLE)
          IF (.NOT.LUC_EXPT%prim_only(g)) THEN ! land-use change may occur

             IF (POP%Iwood(l).EQ.j+1 .AND. &
                  (patch(j+1)%frac+ POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g)) .GT. 0.0   ) THEN
                ! if secondary forest and new secondary forest area > 0
                ! set POLUC diagnostic 2o forest age distribution to POP age distribution
                POPLUC%biomass_age_secondary(g,:)=POP%pop_grid(l)%biomass_age
                ! change in biomass area density due to secondary forest expansion [g C m-2]
                dcExpand(g) = -(POPLUC%gtos(g)+POPLUC%ptos(g))*casapool%cplant(j+1,2)/ &
                     (patch(j+1)%frac + POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g))
                IF (POP%pop_grid(l)%cmass_sum_old .GT. 0.001) THEN
                   ! change in  biomass area density due to secondary harvest and clearance
                   ! (== POP 'cataastrophic' mortality minus natural disturbance minus
                   ! density reduction due to expansion
                   dcHarvClear(g) = -MAX(POP%pop_grid(l)%cat_mortality/ &
                        (POP%pop_grid(l)%cmass_sum_old + POP%pop_grid(l)%growth) &
                        - 1.0/ disturbance_interval ,0.0) &
                        *casapool%cplant(j+1,2) - dcExpand(g)
                   IF (( POPLUC%smharv(g) + POPLUC%syharv(g) + POPLUC%stog(g)).GT.1e-10) THEN
                      ! partition dcHarvClear into harvest and clearance components
                      dcHarv(g) = dcHarvClear(g) * (POPLUC%smharv(g)+POPLUC%syharv(g)) &
                           /( POPLUC%smharv(g) +POPLUC%syharv(g) + POPLUC%stog(g))
                      dcClear(g) = dcHarvClear(g) * POPLUC%stog(g) &
                           /( POPLUC%smharv(g) + POPLUC%syharv(g) + POPLUC%stog(g))

                      IF ((casapool%cplant(j+1,2) + dcExpand(g)+dcClear(g)+dcHarv(g)).GT.0.0) THEN
                         ! flux + A0C0 = (A + dA) * (C + dC)
                         ! flux = A0 * dC + dA( C0 + dC)
                         ! harvest+ clearance flux (not yet corrected for carbon remaining &
                         ! in landscape as litter)
                         !  [g C m-2] (grid-cell basis)
                         FHarvClear(g) = -(patch(j+1)%frac *  (dcHarvClear(g) + dcExpand(g)) + &
                              (POPLUC%ptos(g)+POPLUC%gtos(g)-POPLUC%stog(g))* &
                              ( casapool%cplant(j+1,2) + dcHarvClear(g) + dcExpand(g) ) )
                         ! secondary forest harvest flux (goes to harvest wood products) &
                         ! [g C m-2] (grid-cell basis)
                         !  corrected for carbon remaining in landscape as litter
                         POPLUC%FHarvest(g,2) =  (1.-POPLUC%fracHarvSecResid(g))*FHarvClear(g) &
                              * (POPLUC%smharv(g) + POPLUC%syharv(g)) &
                              /( POPLUC%smharv(g) + POPLUC%syharv(g) + POPLUC%stog(g))
                         ! secondary forest clearance flux (goes to clearance pool) &
                         ! [g C m-2] (grid-cell basis)
                         !  corrected for carbon remaining in landscape as litter
                         POPLUC%FClearance(g,2) =  (1.-POPLUC%fracClearResid(g))* FHarvClear(g) * &
                              POPLUC%stog(g) &
                              /( POPLUC%smharv(g) + POPLUC%syharv(g) + POPLUC%stog(g))
                      ELSE
                         dcHarv(g) = 0.0
                         dcClear(g) = 0.0
                         POPLUC%FHarvest(g,2) = 0.0
                         POPLUC%FClearance(g,2) = 0.0
                      ENDIF

                   ENDIF
                ENDIF

             ENDIF
          ENDIF
          IF (POP%Iwood(l).EQ.j ) THEN
             ! set diagnostic POPLUC variables to POP equivalents &
             !(patch biomass and age distributions)
             POPLUC%biomass_age_primary(g,:)=POP%pop_grid(l)%biomass_age
             POPLUC%freq_age_primary(g,:)=POP%pop_grid(l)%freq_age
          ENDIF
       ENDDO
    ENDDO



    ! Calculate Carbon Pool Transfers
    DO g = 1,POPLUC%np ! loop over POPLUC gridcells (== CABLE gridcells)
       j = landpt(g)%cstart   ! start index of CABLE grid-cell tiles
       l = landpt(g)%cend     ! end index of CABLE grid-cell tiles
       ! initialise local variables
       ! _r refers to receiver tile
       ! _d refers to donor tile
       ! dA refers to change in tile area
       ! dC refers to change in carbon [g m-2] (grid-cell basis)
       ! dn refers to change in nitrogen [g m-2] (grid-cell basis)
       ! dp refers to change in phosphorous [g m-2] (grid-cell basis)
       dA_r = 0.0
       dA_d = 0.0
       dA = 0.0
       dcsoil = 0.0
       dcplant = 0.0
       dclitter = 0.0
       dclabile = 0.0
       dcsoil_r = 0.0
       dcplant_r = 0.0
       dclitter_r = 0.0
       dcsoil_r = 0.0
       dcplant_r = 0.0
       dclitter_r = 0.0
       dclabile_r = 0.0

       dcsoil_d = 0.0
       dcplant_d = 0.0
       dclitter_d = 0.0
       dcsoil_d = 0.0
       dcplant_d = 0.0
       dclitter_d = 0.0
       dclabile_d = 0.0

       dnsoil = 0.0
       dnplant = 0.0
       dnlitter = 0.0
       dnsoil_r = 0.0
       dnplant_r = 0.0
       dnlitter_r = 0.0
       dnsoil_r = 0.0
       dnplant_r = 0.0
       dnlitter_r = 0.0
       dnsoilmin_r = 0.0

       dnlitter = 0.0
       dnsoil_d = 0.0
       dnplant_d = 0.0
       dnlitter_d = 0.0
       dnsoil_d = 0.0
       dnplant_d = 0.0
       dnlitter_d = 0.0
       dnsoilmin_d = 0.0

       dpsoil = 0.0
       dpplant = 0.0
       dplitter = 0.0
       dpsoil_r = 0.0
       dpplant_r = 0.0
       dplitter_r = 0.0
       dpsoil_r = 0.0
       dpplant_r = 0.0
       dplitter_r = 0.0
       dwood_transfer = 0.0

       dpsoil_d = 0.0
       dpplant_d = 0.0
       dplitter_d = 0.0
       dpsoil_d = 0.0
       dpplant_d = 0.0
       dplitter_d = 0.0

       IF (.NOT.LUC_EXPT%prim_only(g)) THEN ! only worry about transitions where land-use change is possible
          DO k = 1,nTrans  ! loop over all possible transition types
             IF (k==1) THEN
                deltaA = POPLUC%ptos(g)        ! transition area: primary to seondary transition
                idp = j                        ! donor patch index
                irp = j+1                      ! receiver patch index
                idlu = p                       ! donor land use index (p = primary)
                irlu = s                       ! receiver land use index (s = secondary)
             ELSEIF (k==2) THEN
                deltaA = POPLUC%ptog(g)        ! transition area: primary to grass (open) transition
                idp = j                        ! donor patch index
                irp = j+2                      ! receiver patch index
                idlu = p                       ! donor land use index (p = primary)
                irlu = gr                      ! receiver land use index (gr = grass or open)
             ELSEIF (k==3) THEN
                deltaA = POPLUC%stog(g)        ! transition area: secondary to grass (open) transition
                idp = j+1                      ! donor patch index
                irp = j+2                      ! receiver patch index
                idlu = s                       ! donor land use index
                irlu = gr                       ! receiver land use index
             ELSEIF(k==4) THEN
                deltaA = POPLUC%gtos(g)        ! transition area: grass (open) to secondary transition
                idp = j+2                      ! donor patch index
                irp = j+1                      ! receiver patch index
                idlu = gr                      ! donor land use index
                irlu = s                       ! receiver land use index (s = secondary)
             ENDIF
             dcsoil   = 0
             dclitter = 0
             dcplant  = 0
             dclabile = 0
             dnsoil   = 0
             dnlitter = 0
             dnplant  = 0
             dpsoil   = 0
             dplitter = 0
             dpplant  = 0
             dwood_transfer = 0.0

             ! transfer fluxes : only consider cases where gross transition area is greater than zero
             IF (deltaA.GT.0.0) THEN
                ! all soil pools

                ! change in carbon associated with gross transition
                dcsoil(irlu,idlu,:) = deltaA*casapool%csoil(idp,:)
                ! change in receiver carbon pool (accumulated over all gross transitions)
                dcsoil_r(irlu,:) =  dcsoil_r(irlu,:) + deltaA*casapool%csoil(idp,:)
                ! change in donor carbon pool (accumulated over all gross transitions)
                dcsoil_d(idlu,:) =  dcsoil_d(idlu,:) - deltaA*casapool%csoil(idp,:)

                dnsoil(irlu,idlu,:) = deltaA*casapool%nsoil(idp,:)
                dnsoil_r(irlu,:) =  dnsoil_r(irlu,:) + deltaA*casapool%nsoil(idp,:)
                dnsoil_d(idlu,:) =  dnsoil_d(idlu,:) - deltaA*casapool%nsoil(idp,:)

                dnsoilmin_r(irlu) = dnsoilmin_r(irlu) + deltaA*casapool%nsoilmin(idp)
                dnsoilmin_d(irlu) = dnsoilmin_d(idlu) - deltaA*casapool%nsoilmin(idp)

                dpsoil(irlu,idlu,:) = deltaA*casapool%psoil(idp,:)
                dpsoil_r(irlu,:) =  dpsoil_r(irlu,:) + deltaA*casapool%psoil(idp,:)
                dpsoil_d(idlu,:) =  dpsoil_d(idlu,:) - deltaA*casapool%psoil(idp,:)
                ! need to inlcude other P pools here too: occluded P and labile P)

                ! microbial litter
                dclitter(irlu,idlu,1) = deltaA*casapool%clitter(idp,1)
                dclitter_r(irlu,1) =  dclitter_r(irlu,1) + dclitter(irlu,idlu,1)
                dclitter_d(idlu,1) =  dclitter_d(idlu,1) - dclitter(irlu,idlu,1)

                dnlitter(irlu,idlu,1) = deltaA*casapool%nlitter(idp,1)
                dnlitter_r(irlu,1) =  dnlitter_r(irlu,1) + dnlitter(irlu,idlu,1)
                dnlitter_d(idlu,1) =  dnlitter_d(idlu,1) - dnlitter(irlu,idlu,1)

                dplitter(irlu,idlu,1) = deltaA*casapool%plitter(idp,1)
                dplitter_r(irlu,1) =  dplitter_r(irlu,1) + dplitter(irlu,idlu,1)
                dplitter_d(idlu,1) =  dplitter_d(idlu,1) - dplitter(irlu,idlu,1)

                ! CWD : donor pool inhreits CWD and residues from harvest/clearance
                IF (idlu == s .AND. casapool%cplant(idp,2).GT.1e-5 ) THEN
                   ! seocndary forest clearance: use secondary forest clearance flux from above
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3) + &
                        POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.-POPLUC%fracClearResid(g))
                   dnlitter(irlu,idlu,3) = deltaA*casapool%nlitter(idp,3) + &
                        POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.-POPLUC%fracClearResid(g))  &
                        *casapool%nplant(idp,2)/ casapool%cplant(idp,2)
                   dplitter(irlu,idlu,3) = deltaA*casapool%plitter(idp,3) + &
                        POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.-POPLUC%fracClearResid(g))  &
                        *casapool%pplant(idp,2)/ casapool%cplant(idp,2)
                ELSEIF (idlu == p .AND. irlu == s) THEN
                   ! primary forest harvest: assume harvest is uniform across age-classes
                   dclitter(irlu,idlu,3) = deltaA*(casapool%clitter(idp,3) + &
                        POPLUC%fracHarvResid(g)*casapool%cplant(idp,2) )
                   dnlitter(irlu,idlu,3) = deltaA*(casapool%nlitter(idp,3) + &
                        POPLUC%fracHarvResid(g)*casapool%nplant(idp,2) )
                   dplitter(irlu,idlu,3) = deltaA*(casapool%plitter(idp,3) + &
                        POPLUC%fracHarvResid(g) *casapool%pplant(idp,2) )
                ELSEIF (idlu == p .AND. irlu == gr) THEN
                   ! primary forest clearance: assume clearance is uniform across age-classes
                   dclitter(irlu,idlu,3) = deltaA*(casapool%clitter(idp,3) + &
                        POPLUC%fracClearResid(g)*casapool%cplant(idp,2) )
                   dnlitter(irlu,idlu,3) = deltaA*(casapool%nlitter(idp,3) + &
                        POPLUC%fracClearResid(g)*casapool%nplant(idp,2) )
                   dplitter(irlu,idlu,3) = deltaA*(casapool%plitter(idp,3) + &
                        POPLUC%fracClearResid(g) *casapool%pplant(idp,2) )
                ELSE
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3)
                   dnlitter(irlu,idlu,3) = deltaA*casapool%nlitter(idp,3)
                   dplitter(irlu,idlu,3) = deltaA*casapool%plitter(idp,3)
                ENDIF
                dclitter_r(irlu,3) =  dclitter_r(irlu,3) + dclitter(irlu,idlu,3)
                dclitter_d(idlu,3) =  dclitter_d(idlu,3) - dclitter(irlu,idlu,3)

                dnlitter_r(irlu,3) =  dnlitter_r(irlu,3) + dnlitter(irlu,idlu,3)
                dnlitter_d(idlu,3) =  dnlitter_d(idlu,3) - dnlitter(irlu,idlu,3)

                dplitter_r(irlu,3) =  dplitter_r(irlu,3) + dplitter(irlu,idlu,3)
                dplitter_d(idlu,3) =  dplitter_d(idlu,3) - dplitter(irlu,idlu,3)

                ! fine structural litter: donor pool inherits leaves and fine roots
                dclitter(irlu,idlu,2) = deltaA*(casapool%clitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
                dclitter_r(irlu,2) =  dclitter_r(irlu,2) + dclitter(irlu,idlu,2)
                dclitter_d(idlu,2) =  dclitter_d(idlu,2) - dclitter(irlu,idlu,2)

                dnlitter(irlu,idlu,2) = deltaA*(casapool%nlitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
                dnlitter_r(irlu,2) =  dnlitter_r(irlu,2) + dnlitter(irlu,idlu,2)
                dnlitter_d(idlu,2) =  dnlitter_d(idlu,2) - dnlitter(irlu,idlu,2)

                dplitter(irlu,idlu,2) = deltaA*(casapool%plitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
                dplitter_r(irlu,2) =  dplitter_r(irlu,2) + dplitter(irlu,idlu,2)
                dplitter_d(idlu,2) =  dplitter_d(idlu,2) - dplitter(irlu,idlu,2)

                ! labile carbon
                dclabile(irlu,idlu) = deltaA*casapool%clabile(idp)
                dclabile_r(irlu) =  dclabile_r(irlu) + deltaA*casapool%clabile(idp)
                dclabile_d(irlu) =  dclabile_d(idlu) - deltaA*casapool%clabile(idp)

                ! biomass: no biomass inherited
                dcplant(irlu,idlu,:) = 0.0
                dcplant_r(irlu,:) =  dcplant_r(irlu,:) + dcplant(irlu,idlu,:)
                dcplant_d(idlu,:) =  dcplant_d(idlu,:) - dcplant(irlu,idlu,:)

                dnplant(irlu,idlu,:) = 0.0
                dnplant_r(irlu,:) =  dnplant_r(irlu,:) + dnplant(irlu,idlu,:)
                dnplant_d(idlu,:) =  dnplant_d(idlu,:) - dnplant(irlu,idlu,:)

                dpplant(irlu,idlu,:) = 0.0
                dpplant_r(irlu,:) =  dpplant_r(irlu,:) + dpplant(irlu,idlu,:)
                dpplant_d(idlu,:) =  dpplant_d(idlu,:) - dpplant(irlu,idlu,:)

                ! Gross Transfer Flux (total C-tranfer associated with  k-th tranistion)
                popluc%FtransferGross(g,k) = &
                     SUM(dcsoil(irlu,idlu,:) + dclitter(irlu,idlu,:) +  dcplant(irlu,idlu,:))

                ! Harvest Flux
                ! primary to secondary forest
                ! (note secondary harvest and clearance fluxes are already evaluated &
                !at the top of this subroutine)
                IF (idlu==p .AND. irlu ==s) THEN
                   popluc%FHarvest(g,idlu) =  (1.0 -POPLUC%fracHarvResid(g)) &
                        *casapool%cplant(idp,2)*deltaA
                ENDIF
                ! Clearance Flux
                ! primary forest to grass
                IF ((idlu==p) .AND. irlu ==gr) THEN
                   popluc%FClearance(g,idlu) =  (1.0 -POPLUC%fracClearResid(g)) &
                        * casapool%cplant(idp,2)*deltaA
                ENDIF

                ! transition area
                dA_r(irlu) = dA_r(irlu) + deltaA
                dA_d(idlu) = dA_d(idlu) - deltaA
             ENDIF
          ENDDO  ! ntrans

          DO ilu=1,nLU
             ! update pools
             irp = ilu + j -1
             dwood_transfer = 0.0
             dA(ilu) = dA_r(ilu) + dA_d(ilu)
             ! Net Transfer Flux
             POPLUC%FTransferNet(g,ilu) = SUM(dcsoil_r(ilu,:) + dcsoil_d(ilu,:)) + &
                  SUM(dclitter_r(ilu,:) + dclitter_d(ilu,:))  + &
                  SUM(dcplant_r(ilu,:)  + dcplant_d(ilu,:)) + &
                  dclabile_r(ilu)  + dclabile_d(ilu)
             IF (ilu ==2) THEN
                ! augment CWD pools in secondary veg tiles by harvest residues
                dclitter_r(ilu,3) = dclitter_r(irlu,3) + &
                     POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0 -POPLUC%fracHarvSecResid(g))
                IF (casapool%cplant(irp,2) .GT. 1.e-5) THEN
                   dnlitter_r(ilu,3) = dnlitter_r(irlu,3) +  &
                        POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0 -POPLUC%fracHarvSecResid(g)) &
                        *casapool%nplant(irp,2)/ casapool%cplant(irp,2)
                   dplitter_r(ilu,3) = dplitter_r(irlu,3) +  &
                        POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0 -POPLUC%fracHarvSecResid(g)) &
                        *casapool%pplant(irp,2)/ casapool%cplant(irp,2)
                ENDIF
             ENDIF

             IF ((patch(irp)%frac+dA(ilu)).GT.1.e-6  ) THEN ! avoid fpe's by ensuring finite &
                ! new tile area

                casapool%nsoilmin(irp) = casapool%nsoilmin(irp) +  &
                     (dnsoilmin_r(ilu) - casapool%nsoilmin(irp)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%clabile(irp) = casapool%clabile(irp) +  &
                     (dclabile_r(ilu) - casapool%clabile(irp)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     (dclabile_r(ilu) - casapool%clabile(irp)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%csoil(irp,:)  = casapool%csoil(irp,:) + &
                     (dcsoil_r(ilu,:)-casapool%csoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     SUM((dcsoil_r(ilu,:)-casapool%csoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu)))

                casapool%nsoil(irp,:)  = casapool%nsoil(irp,:) + &
                     (dnsoil_r(ilu,:)-casapool%nsoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%psoil(irp,:)  = casapool%psoil(irp,:) + &
                     (dpsoil_r(ilu,:)-casapool%psoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))


                casapool%clitter(irp,:)  = casapool%clitter(irp,:) + &
                     (dclitter_r(ilu,:)-casapool%clitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     SUM( (dclitter_r(ilu,:)-casapool%clitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu)))

                casapool%nlitter(irp,:)  = casapool%nlitter(irp,:) + &
                     (dnlitter_r(ilu,:)-casapool%nlitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%plitter(irp,:)  = casapool%plitter(irp,:) + &
                     (dplitter_r(ilu,:)-casapool%plitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%cplant(irp,1)  = casapool%cplant(irp,1) + &
                     (dcplant_r(ilu,1)-casapool%cplant(irp,1)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     (dcplant_r(ilu,1)-casapool%cplant(irp,1)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%nplant(irp,1)  = casapool%nplant(irp,1) + &
                     (dnplant_r(ilu,1)-casapool%nplant(irp,1)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,1)  = casapool%pplant(irp,1) + &
                     (dpplant_r(ilu,1)-casapool%pplant(irp,1)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%cplant(irp,3)  = casapool%cplant(irp,3) + &
                     (dcplant_r(ilu,3)-casapool%cplant(irp,3)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     (dcplant_r(ilu,3)-casapool%cplant(irp,3)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%nplant(irp,3)  = casapool%nplant(irp,3) + &
                     (dnplant_r(ilu,3)-casapool%nplant(irp,3)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,3)  = casapool%pplant(irp,3) + &
                     (dpplant_r(ilu,3)-casapool%pplant(irp,3)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%cplant(irp,2)  = casapool%cplant(irp,2) + &
                     (dcplant_r(ilu,2)-casapool%cplant(irp,2)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     (dcplant_r(ilu,2)-casapool%cplant(irp,2)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%nplant(irp,2)  = casapool%nplant(irp,2) + &
                     (dnplant_r(ilu,2)-casapool%nplant(irp,2)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,2)  = casapool%pplant(irp,2) + &
                     (dpplant_r(ilu,2)-casapool%pplant(irp,2)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))


                ! account here for change in secondary forest biomass density due to:
                !  harvest and clearing, as well as increase in (below ground) CWD
                ! where secondary forest harvest occurs
                IF (ilu .EQ.s .AND. (casapool%cplant(irp,2) +(dcHarv(g)+dCClear(g)) ).GT.0.0  &
                     .AND. casapool%cplant(irp,2).GT.1.e-10  ) THEN


                   casapool%nplant(irp,2) = casapool%nplant(irp,2) + (dcHarv(g)+dcClear(g)) &
                        * casapool%nplant(irp,2)/casapool%cplant(irp,2)
                   casapool%pplant(irp,2) = casapool%pplant(irp,2) + (dcHarv(g)+dcClear(g)) &
                        * casapool%pplant(irp,2)/casapool%cplant(irp,2)
                   casapool%cplant(irp,2) = casapool%cplant(irp,2) + (dcHarv(g)+dcClear(g))


                ELSEIF (ilu .EQ.s .AND. (casapool%cplant(irp,2) +(dcHarv(g)+dCClear(g)) ).LE.0.0  ) THEN

                   POPLUC%FHarvest(g,2) = 0.0
                   POPLUC%FClearance(g,2) = 0.0
                ENDIF

             ENDIF

             IF (patch(irp)%frac .GT. 1e-8) THEN
                casaflux%FluxCtohwp(irp) = POPLUC%FHarvest(g,ilu)/ patch(irp)%frac * ktauday
                casaflux%FluxCtoclear(irp) = POPLUC%FClearance(g,ilu)/ patch(irp)%frac * ktauday
             ENDIF

          ENDDO
          ! POPLUC diagnostics
          ! pools in gC per m2 of gridcell
          ! NEP in g C y-1 per m2 of gridcell
          POPLUC%FNEP(g,:) = casabal%Fcneeyear(j:l)*patch(j:l)%frac  ! note NEE = NEP here
          !(+ve into surface)
          DO ilu=1,nLU
             ! update area weights
             irp = ilu + j -1
             ! Update tile area
             IF (ilu == p) THEN
                POPLUC%primf(g)  = MAX(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)
             ELSEIF (ilu == s) THEN
                POPLUC%secdf(g) = MAX(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)
             ELSEIF (ilu == gr) THEN
                POPLUC%grass(g) = MAX(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)
             ENDIF

             POPLUC%csoil(g,ilu) = SUM(casapool%csoil(irp,:))* &
                  MAX(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)
             POPLUC%clitt(g,ilu) = SUM(casapool%clitter(irp,:))* &
                  MAX(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)
             POPLUC%cbiomass(g,ilu) = SUM(casapool%cplant(irp,:))* &
                  MAX(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)

          ENDDO

       ELSE

          POPLUC%csoil(g,:) = 0.0
          POPLUC%clitt(g,:) = 0.0
          POPLUC%cbiomass(g,:) = 0.0
          POPLUC%FNEP(g,:) = 0.0

          POPLUC%csoil(g,1) = SUM(casapool%csoil(j,:))*patch(j)%frac
          POPLUC%clitt(g,1) = SUM(casapool%clitter(j,:))*patch(j)%frac
          POPLUC%cbiomass(g,1) = SUM(casapool%cplant(j,:))*patch(j)%frac
          POPLUC%FNEP(g,1) =  casabal%Fcneeyear(j)*patch(j)%frac
          POPLUC%primf(g) = patch(j)%frac
          POPLUC%secdf(g) = 0.0

          IF (POPLUC%grass(g).GT.0.0 .AND. l.EQ.j+1) THEN
             POPLUC%csoil(g,3) = SUM(casapool%csoil(l,:))*patch(l)%frac
             POPLUC%clitt(g,3) = SUM(casapool%clitter(l,:))*patch(l)%frac
             POPLUC%cbiomass(g,3) = SUM(casapool%cplant(l,:))*patch(l)%frac
             POPLUC%FNEP(g,3) = casabal%Fcneeyear(l)*patch(l)%frac
             POPLUC%grass(g) = patch(l)%frac
          ENDIF

       ENDIF

       POPLUC%HarvProdLoss(g,:) = kHarvProd * POPLUC%HarvProd(g,:)
       POPLUC%ClearProdLoss(g,:) = kClearProd * POPLUC%ClearProd(g,:)

       DO j=1,3
          POPLUC%HarvProd(g,j) = POPLUC%HarvProd(g,j) + &
               POPLUC%fracHarvProd(g,j)*SUM(POPLUC%FHarvest(g,:)) - POPLUC%HarvProdLoss(g,j)

          POPLUC%ClearProd(g,j) = POPLUC%ClearProd(g,j) + &
               POPLUC%fracClearProd(g,j)*SUM(POPLUC%FClearance(g,:)) - POPLUC%ClearProdLoss(g,j)

       ENDDO


    ENDDO

991 FORMAT(1166(e14.7,2x))

    ! update total carbon pools and "last" pool values for use in carbon balance checks.
    casapool%ctot = SUM(casapool%cplant,2)+SUM(casapool%clitter,2)+ &
         SUM(casapool%csoil,2)+casapool%clabile
    casabal%cplantlast  = casapool%cplant
    casabal%clabilelast = casapool%clabile
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil



  END SUBROUTINE POP_LUC_CASA_transfer
  !*******************************************************************************
  SUBROUTINE POPLUC_Init(POPLUC,LUC_EXPT, casapool, casaflux, casabiome, veg, POP, np)
    USE cable_def_types_mod, ONLY : veg_parameter_type
    USE casaparm, ONLY: LEAF, WOOD, FROOT

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_biome),            INTENT(IN) :: casabiome
    TYPE (veg_parameter_type),    INTENT(IN) :: veg  ! vegetation parameters
    TYPE (POP_TYPE), INTENT(INOUT)     :: POP

    INTEGER(i4b), INTENT(IN) :: np
    INTEGER(i4b) :: g, j, k, l

    CALL alloc_POPLUC(POPLUC,np)

    POPLUC%it = 0
    POPLUC%np = np

    CALL ZeroPOPLUC(POPLUC)

    CALL POPLUC_set_params(POPLUC, LUC_EXPT)

    IF (cable_user%POPLUC_RunType .EQ. 'init') THEN
       POPLUC%frac_primf = LUC_EXPT%primaryf
       POPLUC%primf = LUC_EXPT%primaryf
       POPLUC%grass = LUC_EXPT%grass
       WHERE ((POPLUC%primf + POPLUC%grass) > 1.0) &
            POPLUC%grass = 1.0 - POPLUC%primf
       POPLUC%frac_forest = 1- POPLUC%grass
       POPLUC%freq_age_secondary(:,1) =  MAX(POPLUC%frac_forest - POPLUC%primf, 0.0)
       POPLUC%latitude = patch(landpt(:)%cstart)%latitude
       POPLUC%longitude = patch(landpt(:)%cstart)%longitude

       ! zero biomass in secondary forest tiles (both CASA and POP variables)
       WHERE (veg%iLU ==2)
          casapool%cplant(:,leaf) = 0.01
          casapool%nplant(:,leaf)= casabiome%ratioNCplantmin(veg%iveg,leaf)* casapool%cplant(:,leaf)
          casapool%pplant(:,leaf)= casabiome%ratioPCplantmin(veg%iveg,leaf)* casapool%cplant(:,leaf)

          casapool%cplant(:,froot) = 0.01
          casapool%nplant(:,froot)= casabiome%ratioNCplantmin(veg%iveg,froot)* &
               casapool%cplant(:,froot)
          casapool%pplant(:,froot)= casabiome%ratioPCplantmin(veg%iveg,froot)* &
               casapool%cplant(:,froot)

          casapool%cplant(:,wood) = 0.01
          casapool%nplant(:,wood)= casabiome%ratioNCplantmin(veg%iveg,wood)* casapool%cplant(:,wood)
          casapool%pplant(:,wood)= casabiome%ratioPCplantmin(veg%iveg,wood)* casapool%cplant(:,wood)
          casaflux%frac_sapwood = 1.0
       endwhere


       DO k=1,np
          IF (.NOT.LUC_EXPT%prim_only(k)) THEN
             j = landpt(k)%cstart+1
             DO l=1,SIZE(POP%Iwood)
                IF( POP%Iwood(l) == j) THEN
                   CALL POP_init_single(POP,veg%disturbance_interval,l)
                   EXIT
                ENDIF
             ENDDO
          ENDIF
       ENDDO

    ELSEIF (cable_user%POPLUC_RunType .EQ. 'restart') THEN

       CALL  READ_LUC_RESTART_NC (POPLUC)

       POPLUC%frac_primf = POPLUC%primf
       POPLUC%frac_forest = 1- POPLUC%grass
       POPLUC%latitude = patch(landpt(:)%cstart)%latitude
       POPLUC%longitude = patch(landpt(:)%cstart)%longitude

    ENDIF

    ! set landuse index for secondary forest POP landscapes
    ! (not for 'static') run type because secondary forest dynamics are only
    ! simulated with dynamic land-use forcing
    IF (TRIM(cable_user%POPLUC_RunType) .NE. 'static') THEN
       DO k=1,POP%np
          IF (veg%iLU(POP%Iwood(k)).EQ.2) THEN
             POP%pop_grid(k)%LU = 2
          ENDIF
       ENDDO
    ENDIF


  END SUBROUTINE POPLUC_Init
  !*******************************************************************************
  SUBROUTINE POPLUC_set_patchfrac(POPLUC,LUC_EXPT)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    INTEGER(i4b) :: j, k, l

    DO k=1,POPLUC%np
       j = landpt(k)%cstart
       l = landpt(k)%cend
       IF (.NOT.LUC_EXPT%prim_only(k)) THEN
          patch(j)%frac = POPLUC%primf(k)
          patch(l)%frac = POPLUC%grass(k)
          patch(j+1)%frac = 1.0 -  patch(j)%frac - patch(l)%frac
       ENDIF
    ENDDO


  END SUBROUTINE POPLUC_SET_PATCHFRAC
  !*******************************************************************************
  SUBROUTINE POPLUC_set_params(POPLUC,LUC_EXPT)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    INTEGER(i4b) :: g, np

    np = POPLUC%np
    DO g = 1,np
       POPLUC%fracharvProd(g, 1) = 0.9
       POPLUC%fracharvProd(g, 2) = 0.04
       POPLUC%fracharvProd(g, 3) = 0.06

       POPLUC%fracClearProd(g, 1) = 0.597
       POPLUC%fracClearProd(g, 2) = 0.403
       POPLUC%fracClearProd(g, 3) = 0.00
       POPLUC%fracClearResid(g) = 0.33
       POPLUC%fracHarvResid(g) = 0.79
       POPLUC%fracHarvSecResid(g) = 0.81

       IF (LUC_EXPT%biome(g)==1 .OR. LUC_EXPT%biome(g)==2) THEN
          ! Tropical Evergreen and Tropical Deciduous
          POPLUC%fracharvProd(g, 1) = 0.9
          POPLUC%fracharvProd(g, 2) = 0.04
          POPLUC%fracharvProd(g, 3) = 0.06

          POPLUC%fracClearProd(g, 1) = 0.597
          POPLUC%fracClearProd(g, 2) = 0.403
          POPLUC%fracClearProd(g, 3) = 0.00

          IF (LUC_EXPT%biome(g)==1) POPLUC%fracHarvResid(g) = 0.79
          IF (LUC_EXPT%biome(g)==2) POPLUC%fracHarvResid(g) = 0.86
          IF (LUC_EXPT%biome(g)==1) POPLUC%fracHarvSecResid(g) = 0.71
          IF (LUC_EXPT%biome(g)==2) POPLUC%fracHarvSecResid(g) = 0.81

          POPLUC%fracClearResid(g) = 0.33

       ELSEIF (LUC_EXPT%biome(g).GE.4 .OR. LUC_EXPT%biome(g).LE.10) THEN

          ! Other Forest
          POPLUC%fracharvProd(g, 1) = 0.4
          POPLUC%fracharvProd(g, 2) = 0.24
          POPLUC%fracharvProd(g, 3) = 0.36

          POPLUC%fracClearProd(g, 1) = 0.597
          POPLUC%fracClearProd(g, 2) = 0.2985
          POPLUC%fracClearProd(g, 3) = 0.1045

          IF (LUC_EXPT%ivegp(g)==2) POPLUC%fracHarvResid(g) = 0.83
          IF (LUC_EXPT%ivegp(g)==1 .OR. LUC_EXPT%ivegp(g)==3 ) POPLUC%fracHarvResid(g) = 0.87
          IF (LUC_EXPT%ivegp(g)==4) POPLUC%fracHarvResid(g) = 0.78

          IF (LUC_EXPT%ivegp(g)==2) POPLUC%fracHarvSecResid(g) = 0.75
          IF (LUC_EXPT%ivegp(g)==1 .OR. LUC_EXPT%ivegp(g)==3 ) POPLUC%fracHarvSecResid(g) = 0.82
          IF (LUC_EXPT%ivegp(g)==4) POPLUC%fracHarvSecResid(g) = 0.70


          POPLUC%fracClearResid(g) = 0.33

       ELSEIF (LUC_EXPT%biome(g).EQ.3 .OR. LUC_EXPT%biome(g).GE.11) THEN
          ! savanna and shrub
          POPLUC%fracharvProd(g, 1) = 1.0
          POPLUC%fracharvProd(g, 2) = 0.00
          POPLUC%fracharvProd(g, 3) = 0.00

          POPLUC%fracClearProd(g, 1) = 0.8
          POPLUC%fracClearProd(g, 2) = 0.2
          POPLUC%fracClearProd(g, 3) = 0.0


          IF (LUC_EXPT%biome(g).EQ.13 .OR.(LUC_EXPT%biome(g).EQ.14)  ) THEN
             POPLUC%fracHarvResid(g) = 0.78
             POPLUC%fracHarvSecResid(g) = 0.70


          ELSE
             POPLUC%fracHarvResid(g) = 0.86
             POPLUC%fracHarvSecResid(g) = 0.81

          ENDIF

          POPLUC%fracClearResid(g) = 0.5

       ENDIF

    ENDDO

  END SUBROUTINE POPLUC_SET_PARAMS
  !*******************************************************************************
  SUBROUTINE alloc_POPLUC(POPLUC, arraysize)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE),INTENT(INOUT) :: POPLUC
    INTEGER,            INTENT(IN) :: arraysize

    ALLOCATE(POPLUC%it)
    ALLOCATE(POPLUC%np)
    ALLOCATE(POPLUC%firstyear)
    ALLOCATE(POPLUC%thisyear)
    ALLOCATE(POPLUC%latitude(arraysize))
    ALLOCATE(POPLUC%longitude(arraysize))
    ALLOCATE(POPLUC%n_event(arraysize))
    ALLOCATE(POPLUC%primf(arraysize))
    ALLOCATE(POPLUC%secdf(arraysize))
    ALLOCATE(POPLUC%grass(arraysize))
    ALLOCATE(POPLUC%ptos(arraysize))
    ALLOCATE(POPLUC%ptog(arraysize))
    ALLOCATE(POPLUC%stop(arraysize))
    ALLOCATE(POPLUC%stog(arraysize))
    ALLOCATE(POPLUC%gtop(arraysize))
    ALLOCATE(POPLUC%gtos(arraysize))
    ALLOCATE(POPLUC%ptog(arraysize))
    ALLOCATE(POPLUC%pharv(arraysize))
    ALLOCATE(POPLUC%smharv(arraysize))
    ALLOCATE(POPLUC%syharv(arraysize))
    ALLOCATE(POPLUC%frac_primf(arraysize))
    ALLOCATE(POPLUC%frac_forest(arraysize))
    ALLOCATE(POPLUC%freq_age_primary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%freq_age_secondary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%biomass_age_primary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%biomass_age_secondary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%age_history_secdf(arraysize,LENGTH_SECDF_HISTORY))
    ALLOCATE(POPLUC%area_history_secdf(arraysize,LENGTH_SECDF_HISTORY))
    ALLOCATE(POPLUC%FNEP(arraysize,nLU))
    ALLOCATE(POPLUC%Clitt(arraysize,nLU))
    ALLOCATE(POPLUC%Csoil(arraysize,nLU))
    ALLOCATE(POPLUC%Cbiomass(arraysize,nLU))
    ALLOCATE(POPLUC%FHarvest(arraysize,nLU))
    ALLOCATE(POPLUC%FClearance(arraysize,nLU))
    ALLOCATE(POPLUC%FTransferNet(arraysize,nLU))
    ALLOCATE(POPLUC%FTransferGross(arraysize,nTrans))
    ALLOCATE(POPLUC%HarvProd(arraysize,3))
    ALLOCATE(POPLUC%ClearProd(arraysize,3))
    ALLOCATE(POPLUC%HarvProdLoss(arraysize,3))
    ALLOCATE(POPLUC%ClearProdLoss(arraysize,3))
    ALLOCATE(POPLUC%fracHarvProd(arraysize,3))
    ALLOCATE(POPLUC%fracClearProd(arraysize,3))
    ALLOCATE(POPLUC%fracHarvResid(arraysize))
    ALLOCATE(POPLUC%fracHarvSecResid(arraysize))
    ALLOCATE(POPLUC%fracClearResid(arraysize))

  END SUBROUTINE alloc_POPLUC
  !*******************************************************************************
  ! Exponential distribution
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  REAL(dp) FUNCTION REALExponential(lambda, x)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) ::  x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       REALExponential=0.0
    ELSE
       REALExponential=lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION REALExponential


  !*******************************************************************************
  SUBROUTINE WRITE_LUC_OUTPUT_NC ( POPLUC, ctime, FINAL )

    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    LOGICAL, INTENT(IN)    :: FINAL
    INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, hist_ID, t_ID, nLU_ID, nTrans_ID
    INTEGER   :: i, mp, nprod, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp,t)
    CHARACTER(len=20),DIMENSION(13):: A1
    ! 2 dim integer arrays (mp,t)
    CHARACTER(len=20),DIMENSION(1):: AI1
    ! 3 dim real arrays (mp,age_max,t)
    CHARACTER(len=25),DIMENSION(4) :: A2
    ! 3 dim real arrays (mp,LENGTH_SECDF_HISTORY,t)
    CHARACTER(len=20),DIMENSION(1) :: A3
    ! 3 dim integer arrays (mp,LENGTH_SECDF_HISTORY,t)
    CHARACTER(len=20),DIMENSION(1) :: AI3
    ! 3 dim real arrays (mp,nLU,t)
    CHARACTER(len=20),DIMENSION(7) :: A4
    ! 3 dim real arrays (mp,nTrans,t)
    CHARACTER(len=20),DIMENSION(1) :: A5
    ! 3 dim real arrays (mp,nprod,t)
    CHARACTER(len=20),DIMENSION(4) :: A6

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1))
    INTEGER, SAVE :: VID2(SIZE(A2)),VID3(SIZE(A3)),VIDI3(SIZE(AI3))
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)), VID6(SIZE(A6))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    CHARACTER(len=50) :: RecordDimName
    REAL, ALLOCATABLE :: freq_age_secondary(:,:)
    INTEGER :: g
    LOGICAL :: put_age_vars
    mp = POPLUC%np
    nprod = 3
    put_age_vars=.FALSE.
    ALLOCATE(freq_age_secondary(mp,age_max))

    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'
    A1(4) = 'ptos'
    A1(5) = 'ptog'
    A1(6) = 'stog'
    A1(7) = 'gtop'
    A1(8) = 'gtos'
    A1(9) = 'frac_primf'
    A1(10) = 'frac_forest'
    A1(11) = 'pharv'
    A1(12) = 'smharv'
    A1(13) = 'syharv'


    AI1(1) = 'n_event'

    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'
    A2(3) = 'biomass_age_primary'
    A2(4) = 'biomass_age_secondary'

    A3(1) = 'area_history_secdf'

    AI3(1) = 'age_history_secdf'

    A4(1) = 'FHarvest'
    A4(2) = 'FClearance'
    A4(3) = 'FNEP'
    A4(4) = 'CLitt'
    A4(5) = 'CSoil'
    A4(6) = 'CBiomass'
    A4(7) = 'FTransferNet'

    A5(1) = 'FTransferGross'

    A6(1) = 'HarvProd'
    A6(2) = 'ClearProd'
    A6(3) = 'HarvProdLoss'
    A6(4) = 'ClearProdLoss'


    DO g=1,mp
       IF (SUM(POPLUC%freq_age_secondary(g,:)).GT.1e-12 ) THEN
          freq_age_secondary(g,:) = POPLUC%freq_age_secondary(g,:)/SUM(POPLUC%freq_age_secondary(g,:))
       ELSE
          freq_age_secondary(g,:) = 0.0
       ENDIF
    ENDDO


    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       IF (CABLE_USER%YEARSTART.LT.1000.AND.CABLE_USER%YEAREND.LT.1000) THEN
          WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ELSEIF (CABLE_USER%YEARSTART.LT.1000) THEN
          WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ENDIF
       fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
            TRIM(dum)//'_LUC_out.nc'

       ! Create NetCDF file:
       STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ! Put the file in define mode:
       STATUS = NF90_redef(FILE_ID)

       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

       ! Define dimensions:
       ! Land (number of points)
       STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mage' , age_max , age_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mhist',LENGTH_SECDF_HISTORY , hist_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mLU', nLU , nLU_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mTrans',nTrans , nTrans_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'nprod',nprod , nprod_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ! Define variables
       STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),VIDtime )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       DO i = 1, SIZE(A0)
          STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       DO i = 1, SIZE(A1)
          STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID,t_ID/),VID1(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       DO i = 1, SIZE(AI1)
          STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID,t_ID/),VIDI1(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO
       IF(put_age_vars) THEN
          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,age_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO
       ENDIF
!!$
!!$       DO i = 1, SIZE(A3)
!!$          STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,hist_ID,t_ID/),VID3(i))
!!$          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$       END DO
!!$
!!$       DO i = 1, SIZE(A3)
!!$          STATUS = NF90_def_var(FILE_ID,TRIM(AI3(i)) ,NF90_INT,(/land_ID,hist_ID,t_ID/),VIDI3(i))
!!$          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$       END DO

       DO i = 1, SIZE(A4)
          STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,nLU_ID,t_ID/),VID4(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       DO i = 1, SIZE(A5)
          STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT,(/land_ID,ntrans_ID,t_ID/),VID5(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       DO i = 1, SIZE(A6)
          STATUS = NF90_def_var(FILE_ID,TRIM(A6(i)) ,NF90_FLOAT,(/land_ID,nprod_ID,t_ID/),VID6(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       ! End define mode:
       STATUS = NF90_enddef(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), POPLUC%latitude )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), POPLUC%longitude )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       CALL1 = .FALSE.

    ENDIF ! CALL1



    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), POPLUC%primf,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), POPLUC%secdf,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), POPLUC%grass,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), POPLUC%ptos,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), POPLUC%ptog,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), POPLUC%stog,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), POPLUC%gtop,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), POPLUC%gtos,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), POPLUC%frac_primf,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), POPLUC%frac_forest,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 11), POPLUC%pharv,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 12), POPLUC%smharv,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 13), POPLUC%syharv,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)



    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), POPLUC%n_event,   &
         start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    IF (put_age_vars) THEN
       ! PUT 3D VARS ( mp, mage, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), POPLUC%freq_age_primary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2),freq_age_secondary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), POPLUC%biomass_age_primary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), POPLUC%biomass_age_secondary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ENDIF


    ! PUT 3D VARS ( mp, nLU, t )

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), POPLUC%FHarvest,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), POPLUC%FClearance,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), POPLUC%FNEP,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), POPLUC%Clitt,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), POPLUC%CSoil,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), POPLUC%CBiomass,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), POPLUC%FTransferNet,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    ! PUT 3D VARS ( mp, nTrans, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), POPLUC%FTransferGross,   &
         start=(/ 1,1,CNT /), count=(/ mp,nTrans,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, nprod, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), POPLUC%HarvProd,   &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(2), POPLUC%ClearProd,   &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(3), POPLUC%HarvProdLoss,   &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(4), POPLUC%ClearProdLoss,   &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " POPLUC Output written to ",fname
    ENDIF

  END SUBROUTINE WRITE_LUC_OUTPUT_NC
  !************************************************************************************************************************************
  SUBROUTINE WRITE_LUC_RESTART_NC ( POPLUC, ctime )

    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, nLU_ID, nTrans_ID, i, mp, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp)
    CHARACTER(len=20),DIMENSION(3):: A1
    ! 2 dim real arrays (mp,age_max)
    CHARACTER(len=25),DIMENSION(2) :: A2
    ! 2 dim real arrays (mp,nprod)
    CHARACTER(len=25),DIMENSION(2) :: A3


    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1))
    INTEGER, SAVE :: VID2(SIZE(A2)), VID3(SIZE(A3))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    CHARACTER(len=50) :: RecordDimName
    INTEGER :: g, nprod

    mp = POPLUC%np
    nprod = 3 ! number of product pools

    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'

    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'

    A3(1) = 'HarvProd'
    A3(2) = 'ClearProd'



    ! Get File-Name

    WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
    IF (CABLE_USER%YEARSTART.LT.1000.AND.CABLE_USER%YEAREND.LT.1000) THEN
       WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
    ELSEIF (CABLE_USER%YEARSTART.LT.1000) THEN
       WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
    ENDIF
    fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//'LUC_rst.nc'

    ! Create NetCDF file:
    STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

    ! Define dimensions:
    ! Land (number of points)
    STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mage' , age_max , age_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'nprod' , nprod , nprod_ID)


    ! Define variables

    DO i = 1, SIZE(A0)
       STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,age_ID/),VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A3)
       STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,nprod_ID/),VID3(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO


    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT LAT / LON ( mp )
    STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), POPLUC%latitude )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), POPLUC%longitude )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), POPLUC%primf)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), POPLUC%secdf)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), POPLUC%grass)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, mage, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), POPLUC%freq_age_primary)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2),POPLUC%freq_age_secondary)

    ! PUT 3D VARS ( mp, nprod, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), POPLUC%HarvProd)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID3(2),POPLUC%ClearProd)


    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    WRITE(*,*) " POPLUC Restart written to ",fname


  END SUBROUTINE WRITE_LUC_RESTART_NC
  !*******************************************************************************
  SUBROUTINE READ_LUC_RESTART_NC (POPLUC)
    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC


    INTEGER   :: STATUS, land_dim, mage_dim, nprod_dim
    INTEGER   :: land_ID, age_ID, nLU_ID, nTrans_ID, i, mp, FILE_ID, dID, nprod, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.
    REAL , ALLOCATABLE :: TMP(:), TMP2(:,:), TMP3(:,:)

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp)
    CHARACTER(len=20),DIMENSION(3):: A1
    ! 2 dim real arrays (mp,age_max)
    CHARACTER(len=25),DIMENSION(2) :: A2
    ! 2 dim real arrays (mp,nprod)
    CHARACTER(len=25),DIMENSION(2) :: A3

    mp = POPLUC%np
    nprod = 3
    ALLOCATE(tmp(mp))
    ALLOCATE(tmp2(mp,age_max))
    ALLOCATE(tmp3(mp,nprod))
    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'

    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'

    A3(1) = 'HarvProd'
    A3(2) = 'ClearProd'




    !fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
    !     '_LUC_rst.nc'

    fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//'LUC_rst.nc'
    STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


    ! DIMS
    STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'mage', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mage_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


    STATUS = NF90_INQ_DIMID( FILE_ID, 'nprod', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=nprod_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! READ 1-dimensional fields
    DO i = 1, SIZE(A1)
       STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A1(i)))
       CASE ('primf'      ) ; POPLUC%primf       = TMP
       CASE ('secdf'   ) ; POPLUC%secdf   = TMP
       CASE ('grass'  ) ; POPLUC%grass  = TMP

       END SELECT
    END DO

991 FORMAT(1000(e12.4,2x))
    ! READ 2-dimensional fields (mage)
    DO i = 1, SIZE(A2)
       STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A2(i)))
       CASE ('freq_age_primary' ) ; POPLUC%freq_age_primary = TMP2
       CASE ('freq_age_secondary' ) ; POPLUC%freq_age_secondary = TMP2
       END SELECT
    END DO

    ! READ 3-dimensional fields (nprod)
    DO i = 1, SIZE(A3)
       STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A3(i)))
       CASE ('HarvProd' ) ; POPLUC%HarvProd = TMP3
       CASE ('ClearProd' ) ; POPLUC%ClearProd = TMP3
       END SELECT
    END DO

    STATUS = NF90_CLOSE( FILE_ID )

  END SUBROUTINE READ_LUC_RESTART_NC
  !*******************************************************************************
  SUBROUTINE WRITE_LUC_OUTPUT_GRID_NC ( POPLUC, ctime, FINAL )

    USE cable_IO_vars_module, ONLY: mask, xdimsize, ydimsize , lat_all, lon_all
    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    LOGICAL, INTENT(IN)    :: FINAL
    INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, hist_ID, t_ID, nLU_ID, nTrans_ID, xID, yID
    INTEGER   :: xvID, yvID
    INTEGER   :: i, mp, nprod, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp,t)
    CHARACTER(len=20),DIMENSION(13):: A1
    ! 2 dim integer arrays (mp,t)
    CHARACTER(len=20),DIMENSION(1):: AI1
    ! 3 dim real arrays (mp,age_max,t)
    CHARACTER(len=25),DIMENSION(4) :: A2
    ! 3 dim real arrays (mp,LENGTH_SECDF_HISTORY,t)
    CHARACTER(len=20),DIMENSION(1) :: A3
    ! 3 dim integer arrays (mp,LENGTH_SECDF_HISTORY,t)
    CHARACTER(len=20),DIMENSION(1) :: AI3
    ! 3 dim real arrays (mp,nLU,t)
    CHARACTER(len=20),DIMENSION(7) :: A4
    ! 3 dim real arrays (mp,nTrans,t)
    CHARACTER(len=20),DIMENSION(1) :: A5
    ! 3 dim real arrays (mp,nprod,t)
    CHARACTER(len=20),DIMENSION(4) :: A6

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1))
    INTEGER, SAVE :: VID2(SIZE(A2)),VID3(SIZE(A3)),VIDI3(SIZE(AI3))
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)), VID6(SIZE(A6))
    INTEGER, SAVE :: FILE_ID, CNT = 0, latID, lonID
    CHARACTER(len=50) :: RecordDimName
    REAL, ALLOCATABLE :: freq_age_secondary(:,:)
    INTEGER :: g, k
    LOGICAL :: put_age_vars
    REAL(dp) :: ncmissingr = -1.0e+33
    INTEGER :: ncmissingi = -9999999
    ! REAL :: missing_value = -999999.0 ! for netcdf output
    INTEGER :: nx, ny
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: LandMask
    ! Logical landmask, true for land, false for non-land
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: fieldr
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: fieldi
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmparr1
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmparr2
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmparr3
    mp = POPLUC%np
    nx = xdimsize
    ny = ydimsize
    nprod = 3
    put_age_vars=.FALSE.
    ALLOCATE(freq_age_secondary(mp,age_max))
    ALLOCATE( landmask ( xdimsize, ydimsize) )      ! Local use in this routine (integer)
    ALLOCATE( fieldr ( xdimsize, ydimsize) )  ! field for UNPACK command (reals)
    ALLOCATE( fieldi ( xdimsize, ydimsize) )  ! field for UNPACK command (integers)
    ALLOCATE( tmparr1( xdimsize, ydimsize, nLU) )
    ALLOCATE( tmparr2( xdimsize, ydimsize, nTrans) )
    ALLOCATE( tmparr3( xdimsize, ydimsize, nprod) )
    ! Convert the integer 'mask' into the logical 'landmask'
    WHERE (mask .EQ. 1 )
       landmask = .TRUE.
    ELSEWHERE
       landmask = .FALSE.
    END WHERE
    fieldr = ncmissingr
    fieldi = ncmissingi

    A0(1) = 'local_latitude'
    A0(2) = 'local_longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'
    A1(4) = 'ptos'
    A1(5) = 'ptog'
    A1(6) = 'stog'
    A1(7) = 'gtop'
    A1(8) = 'gtos'
    A1(9) = 'frac_primf'
    A1(10) = 'frac_forest'
    A1(11) = 'pharv'
    A1(12) = 'smharv'
    A1(13) = 'syharv'


    AI1(1) = 'n_event'

    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'
    A2(3) = 'biomass_age_primary'
    A2(4) = 'biomass_age_secondary'

    A3(1) = 'area_history_secdf'

    AI3(1) = 'age_history_secdf'

    A4(1) = 'FHarvest'
    A4(2) = 'FClearance'
    A4(3) = 'FNEP'
    A4(4) = 'CLitt'
    A4(5) = 'CSoil'
    A4(6) = 'CBiomass'
    A4(7) = 'FTransferNet'

    A5(1) = 'FTransferGross'

    A6(1) = 'HarvProd'
    A6(2) = 'ClearProd'
    A6(3) = 'HarvProdLoss'
    A6(4) = 'ClearProdLoss'


    DO g=1,mp
       IF (SUM(POPLUC%freq_age_secondary(g,:)).GT.1e-12 ) THEN
          freq_age_secondary(g,:) = POPLUC%freq_age_secondary(g,:)/SUM(POPLUC%freq_age_secondary(g,:))
       ELSE
          freq_age_secondary(g,:) = 0.0
       ENDIF
    ENDDO


    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       IF (CABLE_USER%YEARSTART.LT.1000.AND.CABLE_USER%YEAREND.LT.1000) THEN
          WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ELSEIF (CABLE_USER%YEARSTART.LT.1000) THEN
          WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ENDIF
       fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
            TRIM(dum)//'_LUC_out.nc'

       ! Create NetCDF file:
       STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ! Put the file in define mode:
       STATUS = NF90_redef(FILE_ID)

       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )


       ! Define dimensions:
       STATUS = NF90_DEF_DIM(FILE_ID, 'x', xdimsize, xID)
       IF (STATUS /= NF90_NOERR) CALL handle_err                                       &
            (STATUS, 'Error defining x dimension in output file. '// &
            '(SUBROUTINE open_output_file)')
       STATUS = NF90_DEF_DIM(FILE_ID, 'y', ydimsize, yID)
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error defining y dimension in output file. '// &
            '(SUBROUTINE open_output_file)')


       STATUS = NF90_def_dim(FILE_ID, 'mage' , age_max , age_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mhist',LENGTH_SECDF_HISTORY , hist_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mLU', nLU , nLU_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mTrans',nTrans , nTrans_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'nprod',nprod , nprod_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ! Define variables
       STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),VIDtime )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       STATUS = NF90_DEF_VAR(FILE_ID, 'latitude', NF90_FLOAT, (/xID, yID/), latID)
       STATUS = NF90_DEF_VAR(FILE_ID, 'longitude', NF90_FLOAT, (/xID, yID/), lonID)

       STATUS = NF90_PUT_ATT(FILE_ID, latID, 'units', 'degrees_north')
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error defining latitude variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')

       STATUS = NF90_PUT_ATT(FILE_ID, lonID, 'units', 'degrees_east')
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error defining longitude variable attributes in output file. '// &
            '(SUBROUTINE open_output_file)')
       ! Write "cordinate variables" to enable reading by GrADS:
       STATUS = NF90_DEF_VAR(FILE_ID, 'x', NF90_FLOAT, (/xID/), xvID)
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error defining "x" variable (for GrADS) in output file. '// &
            '(SUBROUTINE open_output_file)')
       STATUS = NF90_PUT_ATT(FILE_ID, xvID, 'units', 'degrees_east')
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error writing x coordinate variable (GrADS) units in output '// &
            'file. (SUBROUTINE open_output_file)')
       STATUS = NF90_PUT_ATT(FILE_ID, xvID, 'comment',                               &
            'x coordinate variable for GrADS compatibility')
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error writing x variables comment in output file. '// &
            '(SUBROUTINE open_output_file)')
       STATUS = NF90_DEF_VAR(FILE_ID, 'y', NF90_FLOAT, (/yID/), yvID)
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error defining "y" variable (for GrADS) in output file. '// &
            '(SUBROUTINE open_output_file)')
       STATUS = NF90_PUT_ATT(FILE_ID, yvID, 'units', 'degrees_north')
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error writing y coordinate variable (GrADS) units in output '//  &
            'file. (SUBROUTINE open_output_file)')
       STATUS = NF90_PUT_ATT(FILE_ID, yvID, 'comment',                               &
            'y coordinate variable for GrADS compatibility')
       IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
            (STATUS, 'Error writing y variables comment in output file. '// &
            '(SUBROUTINE open_output_file)')



       DO i = 1, SIZE(A0)
          STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/xID, yID/),VID0(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID0(i), '_FillValue', REAL(ncmissingr,4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A0(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID0(i), 'missing_value', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A0(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
       END DO

       DO i = 1, SIZE(A1)
          STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/xID, yID,t_ID/),VID1(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID1(i), '_FillValue', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A1(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID1(i), 'missing_value', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A1(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')

       END DO

       DO i = 1, SIZE(AI1)
          STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/xID, yID,t_ID/),VIDI1(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VIDI1(i), '_FillValue',ncmissingi)
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//AI1(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VIDI1(i), 'missing_value', ncmissingi)
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//AI1(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')

       END DO
       IF(put_age_vars) THEN
          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/xID, yID,age_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO
       ENDIF

       DO i = 1, SIZE(A4)
          STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/xID, yID,nLU_ID,t_ID/),VID4(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID4(i), '_FillValue', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A4(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID4(i), 'missing_value', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A4(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')


       END DO

       DO i = 1, SIZE(A5)
          STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT,(/xID, yID,ntrans_ID,t_ID/),VID5(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID5(i), '_FillValue', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A5(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID5(i), 'missing_value', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A5(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')


       END DO

       DO i = 1, SIZE(A6)
          STATUS = NF90_def_var(FILE_ID,TRIM(A6(i)) ,NF90_FLOAT,(/xID, yID ,nprod_ID,t_ID/),VID6(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID6(i), '_FillValue', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A6(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID6(i), 'missing_value', REAL(ncmissingr, 4))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A6(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')


       END DO

       ! End define mode:
       STATUS = NF90_enddef(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)



       ! Write GrADS coordinate variables
       STATUS = NF90_PUT_VAR(FILE_ID, xvID, REAL(lon_all(:, 1), 4))
       IF(STATUS /= NF90_NOERR) CALL handle_err                                         &
            (STATUS, 'Error writing GrADS x coordinate variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       STATUS = NF90_PUT_VAR(FILE_ID, yvID, REAL(lat_all(1, :), 4))
       IF(STATUS /= NF90_NOERR) CALL handle_err                                         &
            (STATUS, 'Error writing GrADS y coordinate variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')


       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), UNPACK(POPLUC%latitude,landmask, fieldr ))
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), UNPACK(POPLUC%longitude, landmask, fieldr) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, latID, lat_all)
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, lonID, lon_all )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       CALL1 = .FALSE.

    ENDIF ! CALL1



    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), UNPACK(POPLUC%primf,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), UNPACK(POPLUC%secdf,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), UNPACK(POPLUC%grass,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), UNPACK(POPLUC%ptos,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), UNPACK(POPLUC%ptog,landmask, fieldr), &
         start=(/ 1,1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), UNPACK(POPLUC%stog, landmask, fieldr), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), UNPACK(POPLUC%gtop,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), UNPACK(POPLUC%gtos,landmask, fieldr), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), UNPACK(POPLUC%frac_primf,landmask, fieldr), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), UNPACK(POPLUC%frac_forest,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 11), UNPACK(POPLUC%pharv, landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 12), UNPACK(POPLUC%smharv,landmask, fieldr), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 13), UNPACK(POPLUC%syharv,landmask, fieldr), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), UNPACK(POPLUC%n_event, landmask, fieldi), &
         start=(/ 1,1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, nLU, t )
    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FHarvest(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FClearance(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FNEP(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%Clitt(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%Csoil(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%CBiomass(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FTransferNet(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), tmparr1,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( nx,ny, nTrans, t )
    DO k=1,nTrans
       tmparr2(:,:,k) = UNPACK(POPLUC%FTransferGross(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), tmparr2,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nTrans,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( nx,ny, nprod, t )
    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%HarvProd(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), tmparr3,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%ClearProd(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(2), tmparr3,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%HarvProdLoss(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(3), tmparr3,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%ClearProdLoss(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(4), tmparr3,   &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " POPLUC Output written to ",fname
    ENDIF

  END SUBROUTINE WRITE_LUC_OUTPUT_GRID_NC
  !************************************************************************************************************************************




END MODULE POPLUC_MODULE
