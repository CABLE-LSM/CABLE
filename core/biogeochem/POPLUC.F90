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
! POPLUC_set_params(POPLUC,LUC_EXPT)
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

  implicit none

  INTEGER(i4b), PARAMETER :: LENGTH_SECDF_HISTORY = 4000
  INTEGER(i4b), PARAMETER :: AGE_MAX = 1000
  ! N.B. needs to be the same as veg%disturbance_interval
  INTEGER(i4b), PARAMETER :: disturbance_interval = 100
  LOGICAL,      PARAMETER :: IFHARVEST=.FALSE.
  INTEGER(i4b), PARAMETER :: ROTATION=70
  INTEGER(i4b), PARAMETER :: nLU=3    ! number of land-use tiles (pf, sf, grass)
  INTEGER(i4b), PARAMETER :: nTrans=4 ! number of possible gross transition types (ptog, ptos, stog, gtos)

END MODULE POPLUC_CONSTANTS

!*******************************************************************************

MODULE POPLUC_Types

  USE TYPEdef,          ONLY: dp, i4b
  USE POPLUC_Constants, ONLY: LENGTH_SECDF_HISTORY, AGE_MAX

  implicit none

  TYPE POPLUC_TYPE
     INTEGER(i4b), POINTER :: it
     INTEGER(i4b), POINTER :: np
     INTEGER(i4b), POINTER :: firstyear
     INTEGER(i4b), POINTER :: thisyear
     INTEGER(i4b), DIMENSION(:),POINTER :: n_event => null() ! number of secondary forest transitions
     REAL(dp), DIMENSION(:),POINTER :: latitude => null(), longitude => null()
     REAL(dp), DIMENSION(:),POINTER :: primf => null(), secdf => null(), grass => null(), &  ! land cover types
          ptos => null(), ptog => null(), stog => null(), gtop => null(), gtos => null(),    & ! transitions
          frac_primf => null(), frac_forest => null()
     REAL(dp), DIMENSION(:),POINTER :: crop => null(), past => null() ! components of managed grass (crop,pasture)
     ! transitions associated with crop (c) and pasture (q)
     REAL(dp), DIMENSION(:),POINTER :: ptoc => null(), ptoq => null(), stoc => null(), stoq => null(), &
          qtos => null(), ctos => null()
     REAL(dp), DIMENSION(:,:),POINTER ::  freq_age_primary => null(), freq_age_secondary => null(), &
          biomass_age_primary => null(), biomass_age_secondary => null()
     REAL(dp), DIMENSION(:,:),POINTER :: age_history_secdf => null(), area_history_secdf => null()
     REAL(dp), DIMENSION(:,:),POINTER :: FNEP => null(), Clitt => null(), Csoil => null(), Cbiomass => null()
     REAL(dp), DIMENSION(:,:),POINTER :: FHarvest => null(), FClearance => null(), FTransferNet => null()
     REAL(dp), DIMENSION(:,:),POINTER :: FTransferGross => null()
     REAL(dp), DIMENSION(:),POINTER :: pharv => null(), smharv => null(), syharv => null()
     ! ag prod pool (grazing + crop harvest) and loss to atm, loss of C from biosphere due to crop/pasture harvest
     REAL(dp), DIMENSION(:),POINTER :: AgProd => null(), AgProdLoss => null(), FAg => null()
     REAL(dp), DIMENSION(:,:),POINTER :: HarvProd => null(), ClearProd => null() ! wood harvest and clearance pools
     REAL(dp), DIMENSION(:,:),POINTER :: fracHarvProd => null(), fracClearProd => null()
     REAL(dp), DIMENSION(:,:),POINTER :: HarvProdLoss => null(), ClearProdLoss => null()
     REAL(dp), DIMENSION(:),POINTER :: fracHarvResid => null(), fracHarvSecResid => null(), fracClearResid => null()
     REAL(dp), DIMENSION(:),POINTER :: kSecHarv => null(), kNatDist => null(), &
          kExpand1 => null(), kExpand2 => null(), kClear => null()
     ! carbon denisty in sec forest harvest area, relative to tile average
     REAL(dp), DIMENSION(:),POINTER :: cRelClear => null()
     ! biomass density loss rates
     ! For 13C
     ! Residual flux to litter from harvesting primary forest
     real(dp), dimension(:), pointer :: FluxPHarvResidtoLitter => null()
     ! Residual flux to litter from harvesting secondary forest
     real(dp), dimension(:), pointer :: FluxSHarvResidtoLitter => null()
     ! Residual flux to litter from clearing primary forest
     real(dp), dimension(:), pointer :: FluxPClearResidtoLitter => null()
     ! Residual flux to litter from clearing secondary forest
     real(dp), dimension(:), pointer :: FluxSClearResidtoLitter => null()
     ! Harvest and Clearance induced change of plant pool of secondary forest
     real(dp), dimension(:), pointer :: dcSHarvClear => null()
  END TYPE POPLUC_TYPE

END MODULE POPLUC_Types

!*******************************************************************************

MODULE POPLUC_Module

  !-------------------------------------------------------------------------------
  ! * This module contains all subroutines for POPLUC calcs at a single time step.
  !-------------------------------------------------------------------------------
  USE TYPEdef,              ONLY: dp, sp, i4b
  USE POPLUC_Types
  USE POPLUC_Constants
  USE casavariable,         ONLY: casa_pool, casa_balance, casa_flux, casa_biome
  USE POP_Types,            ONLY: POP_TYPE
  USE cable_common_module,  ONLY: cable_user
  USE cable_IO_vars_module, ONLY: landpt, patch, wlogn
  USE CABLE_LUC_EXPT,       ONLY: LUC_EXPT_TYPE
  USE POPModule,            ONLY: pop_init_single

  implicit none

  real(dp), dimension(3) :: kHarvProd, kClearProd
  real(dp)               :: kAgProd

CONTAINS

  !*******************************************************************************

  SUBROUTINE ZeroPOPLUC(POPLUC)

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER:: g,np

    np = popluc%np

    POPLUC%firstyear               = 0_i4b
    POPLUC%thisyear                = 0_i4b
    POPLUC%primf                   = 0._dp
    POPLUC%secdf                   = 0._dp
    POPLUC%grass                   = 0._dp
    POPLUC%crop                    = 0._dp
    POPLUC%past                    = 0._dp
    POPLUC%ptos                    = 0._dp
    POPLUC%ptog                    = 0._dp
    POPLUC%stog                    = 0._dp
    POPLUC%gtop                    = 0._dp
    POPLUC%gtos                    = 0._dp
    POPLUC%ptoc                    = 0._dp
    POPLUC%ptoq                    = 0._dp
    POPLUC%stoc                    = 0._dp
    POPLUC%stoq                    = 0._dp
    POPLUC%qtos                    = 0._dp
    POPLUC%ctos                    = 0._dp
    POPLUC%frac_forest             = 0._dp
    POPLUC%frac_primf              = 0._dp
    POPLUC%area_history_secdf      = 0._dp
    POPLUC%age_history_secdf       = 0._dp
    POPLUC%n_event                 = 0._dp
    POPLUC%freq_age_secondary      = 0._dp
    POPLUC%freq_age_primary        = 0._dp
    POPLUC%biomass_age_primary     = 0._dp
    POPLUC%biomass_age_secondary   = 0._dp
    POPLUC%FNEP                    = 0._dp
    POPLUC%Clitt                   = 0._dp
    POPLUC%Csoil                   = 0._dp
    POPLUC%Cbiomass                = 0._dp
    POPLUC%FHarvest                = 0._dp
    POPLUC%FClearance              = 0._dp
    POPLUC%FTransferNet            = 0._dp
    POPLUC%FTransferGross          = 0._dp
    POPLUC%pharv                   = 0._dp
    POPLUC%smharv                  = 0._dp
    POPLUC%syharv                  = 0._dp
    POPLUC%HarvProd                = 0._dp
    POPLUC%HarvProdLoss            = 0._dp
    POPLUC%ClearProd               = 0._dp
    POPLUC%ClearProdLoss           = 0._dp
    POPLUC%kSecHarv                = 0._dp
    POPLUC%kNatDist                = 0._dp
    POPLUC%kExpand1                = 0._dp
    POPLUC%kExpand2                = 0._dp
    POPLUC%kClear                  = 0._dp
    POPLUC%cRelClear               = 0._dp
    POPLUC%AgProd                  = 0._dp
    POPLUC%AgProdLoss              = 0._dp
    POPLUC%FAg                     = 0._dp
    popluc%FluxPHarvResidtoLitter  = 0._dp
    popluc%FluxSHarvResidtoLitter  = 0._dp
    popluc%FluxPClearResidtoLitter = 0._dp
    popluc%FluxSClearResidtoLitter = 0._dp
    popluc%dcSHarvClear            = 0._dp
    
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
    REAL :: tmp, tmp1, tmp2

    ! frac_open_grid = 1.0_dp-POPLUC%frac_forest(g)
    frac_open_grid = POPLUC%grass(g)
    n = POPLUC%n_event(g)
    IF (from_state=='PRIMF') THEN

       IF (frac_change_grid.GT.POPLUC%primf(g)) THEN
          ! PRINT*, "Warning: requested reduction in primary forest area exceeds primary forest area"
          IF (to_state=='SECDF') POPLUC%ptos(g) = POPLUC%primf(g)
          IF (to_state=='C3ANN') POPLUC%ptog(g) = POPLUC%primf(g)
          frac_change_grid = POPLUC%primf(g)
          POPLUC%primf(g) = 0.0_dp
       ELSE
          POPLUC%primf(g) = POPLUC%primf(g) - frac_change_grid
       ENDIF

       IF (to_state=='SECDF') THEN
          ! Transition from primary -> secondary forest(ptos)
          IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
             tmp = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
                  /sum( POPLUC%freq_age_secondary(g,:))
          ELSE
             tmp = 0.0_dp
          ENDIF
          POPLUC%kExpand1(g) = 0.0_dp
          POPLUC%n_event(g) = POPLUC%n_event(g)+1
          n = POPLUC%n_event(g)
          POPLUC%area_history_secdf(g,n) = frac_change_grid
          POPLUC%age_history_secdf(g,n) = 0.0_dp
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_change_grid

          IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
             tmp1 = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
                  /sum( POPLUC%freq_age_secondary(g,:))
          ELSE
             tmp1 = 0.0_dp
          ENDIF
          if (tmp.gt.0.0_dp) then
             tmp2 = (tmp1-tmp)/tmp
             POPLUC%kExpand1(g) = -tmp2
          else
             POPLUC%kExpand1(g) = 0.0_dp
          endif
       ENDIF

    ELSEIF (from_state=='SECDF') THEN

       IF (to_state=='PRIMF') THEN
          write(*,*) "Error: cannot create primary forest from secondary forest"
          STOP
       ELSE
          ! Transition from secondary -> non forest (stog)
          ! Assumption: youngest stands cleared first, if stands> 10 years accommodate transition
          ! Otherwise uniform harvest

          if (sum( POPLUC%freq_age_secondary(g,10:age_max)) .gt. frac_change_grid) then
             IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
                tmp = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
                     /sum( POPLUC%freq_age_secondary(g,:))
             ELSE
                tmp = 0.0_dp
             ENDIF
             POPLUC%kClear(g) = 0.0_dp
             POPLUC%CRelClear(g) = 0.0_dp
             remaining = frac_change_grid
             i = 10
             DO WHILE (remaining > 0.0_dp .and. i <= age_max )
                IF (POPLUC%freq_age_secondary(g,i).GE.remaining) THEN
                   POPLUC%freq_age_secondary(g,i) =POPLUC%freq_age_secondary(g,i) &
                        - remaining
                   POPLUC%CRelClear(g) = POPLUC%biomass_age_secondary(g,i)*remaining
                   remaining = 0.0_dp
                ELSE
                   remaining = remaining - POPLUC%freq_age_secondary(g,i)
                   POPLUC%freq_age_secondary(g,i) = 0.0_dp
                   POPLUC%CRelClear(g) =  POPLUC%CRelClear(g) &
                        + POPLUC%biomass_age_secondary(g,i)*POPLUC%freq_age_secondary(g,i)
                   i = i+1
                ENDIF

             ENDDO
             if (tmp .gt. 0.0_dp .and. frac_change_grid.gt.0.0_dp) then
                POPLUC%CRelClear(g)  =  POPLUC%CRelClear(g) /frac_change_grid/ tmp
             else
                POPLUC%CRelClear(g) = 0.0_dp
             endif

             if (remaining.gt.0.0_dp) POPLUC%stog(g) = POPLUC%stog(g)-remaining

             IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
                tmp1 = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
                     /sum( POPLUC%freq_age_secondary(g,:))
             ELSE
                tmp1 = 0.0_dp
             ENDIF
             if (tmp.gt.0.0_dp) then
                tmp2 = (tmp1-tmp)/tmp
                POPLUC%kClear(g) = -tmp2
             else
                POPLUC%kClear(g) = 0.0_dp
             endif
          else
             ! uniform clearance across age classes

             POPLUC%kClear(g) = 0.0_dp ! no change in biomass density
             frac_change_grid = min(sum( POPLUC%freq_age_secondary(g,:)), frac_change_grid)
             POPLUC%stog(g) = frac_change_grid
             !if (g==3) write(*,*) 'b4 age_sec', frac_change_grid,  sum(POPLUC%freq_age_secondary(g,:))
             POPLUC%freq_age_secondary(g,:) = POPLUC%freq_age_secondary(g,:)*(1.0_dp-frac_change_grid)
             !if (g==3) write(*,*) 'after age_sec',POPLUC%stog(g),  sum(POPLUC%freq_age_secondary(g,:))
             POPLUC%CRelClear(g) = 1.0_dp
          endif

       ENDIF

    ELSEIF (to_state=='SECDF') THEN
       IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
          tmp = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
               /sum( POPLUC%freq_age_secondary(g,:))
       ELSE
          tmp = 0.0_dp
       ENDIF
       POPLUC%kExpand2(g) = 0.0_dp

       POPLUC%n_event(g) = POPLUC%n_event(g)+1
       n = POPLUC%n_event(g)
       ! Transition from non-forest to secondary forest (gtos)
!!$if (g==4) then
!!$write(*,*) 'gtos1', frac_change_grid, frac_open_grid, POPLUC%grass(g) , sum( POPLUC%freq_age_secondary(g,:))
!!$endif
       if (frac_change_grid.LE.frac_open_grid) THEN
          POPLUC%area_history_secdf(g,n) = frac_change_grid
          POPLUC%age_history_secdf(g,n) = 0.0_dp
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_change_grid
       ELSE
          POPLUC%area_history_secdf(g,n) = frac_open_grid
          POPLUC%age_history_secdf(g,n) = 0.0_dp
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_open_grid
          ! gtos to frac_change_grid here!!
          POPLUC%gtos(g) =  frac_open_grid
       ENDIF

       IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
          tmp1 = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
               /sum( POPLUC%freq_age_secondary(g,:))
       ELSE
          tmp1 = 0.0_dp
       ENDIF

         if (tmp.gt.0.0_dp) then
            tmp2 = (tmp1-tmp)/tmp
            POPLUC%kExpand2(g) = -tmp2
         else
            POPLUC%kExpand2(g) = 0.0_dp
         endif

      ELSEIF (to_state=='PRIMF') THEN

       write(*,*) "Error: cannot create primary forest from non-forest"
       STOP

    ENDIF


  ENDSUBROUTINE execute_luc_event


  !*******************************************************************************

  subroutine calculate_weights(POPLUC, g)

    ! Calculates weights (fraction of total forest area on grid cell)
    !for primary and secondary forest stands up to specified maximum stand age

    implicit none

    type(popluc_type), intent(inout) :: POPLUC
    integer(i4b),      intent(in)    :: g
    
    real(dp), parameter :: eps= 1.e-9_dp
    integer(i4b) :: age, i, iage
    real(dp)     :: fac, disturbance_freq

    ! First get relative weights for primary forest
    disturbance_freq = 1.0_dp / real(disturbance_interval,dp)

    !fac = POPLUC%frac_primf/POPLUC%frac_forest
    fac = 1.0_dp
    do iage=1, age_max
       POPLUC%freq_age_primary(g,iage)   = REALExponential(disturbance_freq, real(iage-1,dp))
       POPLUC%freq_age_secondary(g,iage) = 0.0_dp
    end do
    POPLUC%freq_age_primary(g,:) = POPLUC%freq_age_primary(g,:) / sum(POPLUC%freq_age_primary(g,:)) * fac

    !  Loop through secondary forest stands to transfer weights
    fac = 1.0_dp
    do i=1, POPLUC%n_event(g)
       age = POPLUC%age_history_secdf(g,i)
       POPLUC%freq_age_secondary(g,age+1) = POPLUC%area_history_secdf(g,i)/fac
    enddo

  end subroutine calculate_weights

  !*******************************************************************************

  SUBROUTINE INCREMENT_AGE(POPLUC,g)

    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b) :: n_event, i, j
    INTEGER(i4b), INTENT(IN) :: g
    REAL(dp):: area, remaining
    REAL(dp):: tmp, tmp1, tmp2, tmp3
    n_event =  POPLUC%n_event(g)
    POPLUC%kSecHarv(g) = 0.0_dp
    POPLUC%kNatDist(g) = 0.0_dp

    POPLUC%freq_age_secondary(g,2:age_max)=POPLUC%freq_age_secondary(g,1:age_max-1)
    POPLUC%biomass_age_secondary(g,2:age_max)=POPLUC%biomass_age_secondary(g,1:age_max-1)

    POPLUC%freq_age_secondary(g,1) = 0.0_dp
    POPLUC%biomass_age_secondary(g,1) = 0.0_dp
    ! adjust secondary age distribution for secondary forest harvest area

    if (POPLUC%smharv(g)+POPLUC%syharv(g).gt.0) then
      ! if (sum(POPLUC%biomass_age_secondary(g,:)).gt.0.5*1000) then ! only harvest if biomass density > 0.5 kg Cm-2
       if (sum(POPLUC%biomass_age_secondary(g,:)* POPLUC%freq_age_secondary(g,:)).gt.0.5_dp) then
          ! only harvest if biomass density > 0.5 kg Cm-2
          IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
             tmp = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
                  /sum( POPLUC%freq_age_secondary(g,:))
          ELSE
             tmp = 0.0_dp
          ENDIF

          remaining = POPLUC%smharv(g)+POPLUC%syharv(g)
          i = age_max
          do while (remaining.gt.1.e-10_dp)

             if (POPLUC%freq_age_secondary(g,i) .gt. remaining) then
                POPLUC%freq_age_secondary(g,i) = POPLUC%freq_age_secondary(g,i) - remaining
                POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + remaining
                remaining = 0.0_dp
             else


                POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + &
                     POPLUC%freq_age_secondary(g,i)
                remaining = remaining - POPLUC%freq_age_secondary(g,i)
                POPLUC%freq_age_secondary(g,i) = 0.0_dp

                i = i-1;
             end if
             if (i.lt.2) then
                POPLUC%smharv(g) = POPLUC%smharv(g)+POPLUC%syharv(g) - remaining
                POPLUC%syharv(g) = 0.0_dp
                remaining = 0.0_dp
             endif

          enddo


          IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
             tmp1 = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
                  /sum( POPLUC%freq_age_secondary(g,:))
          ELSE
             tmp1 = 0.0_dp
          ENDIF
          if (tmp.gt.0.0_dp) then
             tmp2 = (tmp1-tmp)/tmp
             POPLUC%kSecHarv(g) = -tmp2
          else
             POPLUC%kSecHarv(g) = 0.0_dp
          endif

       endif
    endif

    ! remove IFHARVEST ?
    IF (IFHARVEST .and. POPLUC%freq_age_secondary(g,ROTATION+1).gt.0.0_dp) THEN
       POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + &
            POPLUC%freq_age_secondary(g,ROTATION+1)
       POPLUC%freq_age_secondary(g,ROTATION+1) = 0.0_dp
    ENDIF

    IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
       tmp = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
            /sum( POPLUC%freq_age_secondary(g,:))
    ELSE
       tmp = 0.0_dp
    ENDIF
    ! adjust secondary age distribution for natural disturbance
    i = age_max
    DO i = age_max, 2 , -1
       POPLUC%freq_age_secondary(g,1) =  POPLUC%freq_age_secondary(g,1) +  &
            POPLUC%freq_age_secondary(g,i)/disturbance_interval
       POPLUC%freq_age_secondary(g,i) = POPLUC%freq_age_secondary(g,i)* &
            (1._dp - 1._dp/disturbance_interval)
    ENDDO

    IF (sum( POPLUC%freq_age_secondary(g,:)) .gt. 0.0_dp) THEN
       tmp1 = sum( POPLUC%freq_age_secondary(g,:)*POPLUC%biomass_age_secondary(g,:)) &
            /sum( POPLUC%freq_age_secondary(g,:))
    ELSE
       tmp1 = 0.0_dp
    ENDIF

    if (tmp.gt.0.0_dp) then
       tmp2 = (tmp1-tmp)/tmp
       POPLUC%kNatDist(g) = -tmp2
    else
       POPLUC%kNatDist(g) = 0.0_dp
    endif


  END SUBROUTINE INCREMENT_AGE

  !*******************************************************************************

  SUBROUTINE POPLUCStep(POPLUC,year)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: year
    INTEGER(i4b) :: g,j

    POPLUC%it = POPLUC%it + 1

    ! vh ! remove code relating to (year.lt.POPLUC%thisyear) ?
    IF (year.lt.POPLUC%thisyear) THEN

       DO g = 1,POPLUC%np
          ! CALL calculate_weights(POPLUC, g)
          CALL increment_age(POPLUC,g)
       ENDDO

    ELSE

       DO g = 1,POPLUC%np

!!$if (POPLUC%thisyear.gt.860) then
!!$         POPLUC%smharv = 0.0_dp ! test
!!$         POPLUC%stog = 0.0_dp ! test
!!$         POPLUC%syharv = 0.0_dp ! test
!!$         POPLUC%ptog = 0.0_dp ! test
!!$         POPLUC%ptos = 0.0_dp ! test
!!$         POPLUC%gtos = 0.0_dp ! test
!!$endif
         POPLUC%kClear(g) = 0.0_dp
         POPLUC%kExpand1(g) = 0.0_dp
         POPLUC%kExpand2(g) = 0.0_dp
         POPLUC%kSecHarv(g) = 0.0_dp
         !CALL increment_age(POPLUC,g)

         if (POPLUC%stog(g) .gt.0.0_dp) &
               CALL execute_luc_event('SECDF','C3ANN',POPLUC%stog(g),g,POPLUC)

         if (POPLUC%ptos(g) .gt. 0.0_dp) &
              CALL execute_luc_event('PRIMF','SECDF',POPLUC%ptos(g),g,POPLUC)

         if (POPLUC%ptog(g) .gt. 0.0_dp) &
              CALL execute_luc_event('PRIMF','C3ANN',POPLUC%ptog(g),g,POPLUC)

         ! if (POPLUC%stog(g) .gt.0.0_dp) &
         !      CALL execute_luc_event('SECDF','C3ANN',POPLUC%stog(g),g,POPLUC)

         if (POPLUC%gtop(g) .gt.0.0_dp) &
              CALL execute_luc_event('C3ANN','PRIMF',POPLUC%gtop(g),g,POPLUC)
         if (POPLUC%gtos(g) .gt.0.0_dp) &
              CALL execute_luc_event('C3ANN','SECDF',POPLUC%gtos(g),g,POPLUC)

         POPLUC%frac_forest(g) =  POPLUC%primf(g)+ SUM(POPLUC%freq_age_secondary(g,:))
!!$if (g==4) write(*,*) 'fracfor: ',  POPLUC%frac_forest(g),  SUM(POPLUC%freq_age_secondary(g,:))
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
    
    TYPE(POPLUC_TYPE),   INTENT(IN)    :: POPLUC
    TYPE(POP_TYPE),      INTENT(INOUT) :: POP
    TYPE(LUC_EXPT_TYPE), INTENT(IN)    :: LUC_EXPT
    integer:: g, k, j, l
    REAL(dp), DIMENSION(:), ALLOCATABLE:: freq_age

    DO l = 1, POP%np
       pop%pop_grid(l)%freq_age = 0.0_dp
    ENDDO
    ALLOCATE(freq_age(age_max))

    DO g = 1,POPLUC%np   ! loop over POPLUC gricells (same as CABLE grid-cells)
       IF (.NOT.LUC_EXPT%prim_only(g) .and.sum(POPLUC%freq_age_secondary(g,:)).gt.1.e-12_dp) THEN
          ! check if tile is secondary forest
          j = landpt(g)%cstart + 1       ! cable index for secondary forest tile
          freq_age = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
          DO l = 1, POP%np   ! for each wooded tile
             if (POP%Iwood(l).eq.j) then
                ! check if cable index in pop structure (Iwood) is the target cable index
                pop%pop_grid(l)%freq_age = freq_age
             endif
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
    
    TYPE(POPLUC_TYPE),   INTENT(INOUT) :: POPLUC
    TYPE(POP_TYPE),      INTENT(INOUT) :: POP
    TYPE(LUC_EXPT_TYPE), INTENT(IN)    :: LUC_EXPT
    TYPE(casa_pool),     INTENT(INOUT) :: casapool
    TYPE(casa_balance),  INTENT(INOUT) :: casabal
    TYPE(casa_flux),     INTENT(INOUT) :: casaflux
    INTEGER,             INTENT(IN)    :: ktauday
    
    ! number of cable time-steps in a day (for needed for LUC flux output)
    integer:: g, k, j, l, idp, irp, idlu, irlu, ilu
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
    REAL(dp), ALLOCATABLE :: dcHarvCLear(:), dcHarv(:), dcClear(:), dcExpand(:), &
         dcNat(:),  FHarvClear(:), FDist(:), dcExpand1(:), dcExpand2(:), &
         FNatDist(:), FHarv(:), FClear(:)
    ! 13C REAL(dp) :: kHarvProd(3), kClearProd(3), kAgProd
    REAL(dp) :: NatDist_loss, Expand_Loss, Clear_Loss, SecHarv_Loss , Dist_Loss, &
         Expand1_Loss, Expand2_Loss, scalefac, tmp
#ifdef __C13DEBUG__
    real(dp) :: tmp_dplant(nLU,3), tmp_tplant(nLU,3), tmp_dlit(nLU,3), tmp_slit(nLU,3)
    real(dp) :: tmp_dsoil(nLU,3)
    integer  :: iwtile, iwpool
#endif

    ! turnover rates for harvest and clearance products (y-1)
    kHarvProd(1) = 1.0_dp/1.0_dp
    kHarvProd(2) = 1.0_dp/10.0_dp
    kHarvProd(3) = 1.0_dp/100.0_dp
    kClearProd = kHarvProd
    kAgProd = 1.0_dp/1.0_dp

    ! zero POPLUC fluxes
    popluc%FtransferNet = 0.0_dp
    popluc%FtransferGross = 0.0_dp
    popluc%FHarvest = 0.0_dp
    popluc%FClearance = 0.0_dp

    ! local variable for storing sum of biomass change due to
    !secondary harvest, clearance and expansion, and secondary forest
    !  harvest and clearance fluxes
    Allocate(FDist(POPLUC%np))
    Allocate(FClear(POPLUC%np))
    Allocate(FHarv(POPLUC%np))
    Allocate(FNatDist(POPLUC%np))
    Allocate(FHarvClear(POPLUC%np))
    Allocate(dcHarvClear(POPLUC%np))
    Allocate(dcHarv(POPLUC%np))
    Allocate(dcClear(POPLUC%np))
    Allocate(dcExpand(POPLUC%np))
    Allocate(dcExpand1(POPLUC%np))
    Allocate(dcExpand2(POPLUC%np))
    Allocate(dcNat(POPLUC%np))
    FHarvClear = 0.0_dp
    FDist = 0.0_dp
    FNatDist = 0.0_dp
    FHarv = 0.0_dp
    FClear = 0.0_dp
    dcHarvClear = 0.0_dp
    dcHarv = 0.0_dp
    dcClear = 0.0_dp
    dcExpand = 0.0_dp
    dcExpand1 = 0.0_dp ! change in sec for carbon density by p->s
    dcExpand2 = 0.0_dp ! change in sec for carbon density by g->s
    dcNat = 0.0_dp
    POPLUC%FHarvest = 0.0_dp
    POPLUC%FClearance = 0.0_dp
    casaflux%FluxCtohwp = 0.0_dp
    casaflux%FluxCtoclear = 0.0_dp
    casaflux%CtransferLUC = 0.0_dp

    popluc%FluxPHarvResidtoLitter  = 0.0_dp
    popluc%FluxSHarvResidtoLitter  = 0.0_dp
    popluc%FluxPClearResidtoLitter = 0.0_dp
    popluc%FluxSClearResidtoLitter = 0.0_dp
    popluc%dcSHarvClear            = 0.0_dp
    
    DO g = 1,POPLUC%np  ! loop over CABLE grid-cells
       dcHarv(g) = 0.0_dp
       dcClear(g) = 0.0_dp
       POPLUC%FHarvest(g,2) = 0.0_dp
       POPLUC%FClearance(g,2) = 0.0_dp
       j = landpt(g)%cstart ! index of primary veg tile in each grid-cell
       DO l = 1, POP%np   ! loop over all POP tiles (== wooded tiles in CABLE)
          IF (.NOT.LUC_EXPT%prim_only(g)) THEN ! land-use change may occur

             if (POP%Iwood(l).eq.j+1 .and. &
                  (patch(j+1)%frac+ POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g)) .gt. 0.0_dp) then
                ! fraction biomass density loss to disturbance
                ! includes natural disturbance, expansion, harvest, clearing
                if ((POP%pop_grid(l)%cmass_sum_old+POP%pop_grid(l)%growth).gt.1.e-6_dp) then
                   Dist_loss = POP%pop_grid(l)%cat_mortality &
                        /(POP%pop_grid(l)%cmass_sum_old+POP%pop_grid(l)%growth)
                else
                   Dist_loss = 0.0_dp
                end if
                POPLUC%biomass_age_secondary(g,:)=POP%pop_grid(l)%biomass_age
             endif

             if (POP%Iwood(l).eq.j+1 .and. &
                  (patch(j+1)%frac+ POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g)) .gt. 0.0_dp  ) then
                ! if secondary forest and new secondary forest area > 0
                ! set POPLUC diagnostic 2o forest age distribution to POP age distribution

                ! change in biomass area density due to secondary forest expansion [g C m-2]
                dcExpand(g) = -(POPLUC%gtos(g)+POPLUC%ptos(g))*casapool%cplant(j+1,2)/ &
                     (patch(j+1)%frac + POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g))
                if (abs((-POPLUC%kNatDist(g) - POPLUC%kSecHarv(g) -  POPLUC%kClear(g))) .gt.0.0_dp) then
                   scalefac = (-Dist_loss -   dcExpand(g)/casapool%cplant(j+1,2))/ &
                        (-POPLUC%kNatDist(g) - POPLUC%kSecHarv(g) -  POPLUC%kClear(g))
                else
                   scalefac = 0.0_dp
                endif

                ! relative  changes in biomass density due to sec forest
                ! harvest, clearing, natural disturbance
                if (scalefac.gt.0) then
                   SecHarv_loss =  scalefac * POPLUC%kSecHarv(g)
                   NatDist_Loss =  scalefac *  POPLUC%kNatDist(g)
                   Clear_loss = scalefac * POPLUC%kClear(g)
                elseif ( abs(-POPLUC%kNatDist(g) - POPLUC%kSecHarv(g)).gt.0.0_dp ) then
                   scalefac = (-Dist_loss -   dcExpand(g)/casapool%cplant(j+1,2))/ &
                        (-POPLUC%kNatDist(g) - POPLUC%kSecHarv(g))
                   Clear_Loss = 0.0_dp
                   SecHarv_loss =  scalefac * POPLUC%kSecHarv(g)
                   NatDist_Loss =  scalefac *  POPLUC%kNatDist(g)
                   if (scalefac.lt.0.0_dp) then
                      SecHarv_loss =  0.0_dp
                      NatDist_Loss =  0.0_dp
                   endif
                else
                   Clear_Loss = 0.0_dp
                   SecHarv_loss =  0.0_dp
                   NatDist_Loss =  0.0_dp

                endif

                ! absolute  changes in biomass density due to sec forest
                ! harvest, clearing, natural disturbance
                dcHarvClear(g) = -(Clear_loss+SecHarv_loss) &
                     *casapool%cplant(j+1,2)

                dcHarv(g) = -(SecHarv_loss) &
                     *casapool%cplant(j+1,2)

                dcClear(g) = -(Clear_loss) &
                     *casapool%cplant(j+1,2)

                dcNat(g) = -NatDist_loss &
                     *casapool%cplant(j+1,2)
!!$               if (g==4) write(*,*) 'b4 fdist', &
!!$               scalefac, casapool%cplant(j+1,2) + dcExpand(g)+dcClear(g)+dcHarv(g),  &
!!$               casapool%cplant(j+1,2) , dcExpand(g),dcClear(g),dcHarv(g), dcnat(g)

                if ((casapool%cplant(j+1,2) + dcExpand(g)+dcClear(g)+dcHarv(g)).gt.0.0_dp .and. &
                     dcnat(g).lt.0.0_dp) then
                   ! flux + A0C0 = (A + dA) * (C + dC)
                   ! flux = A0 * dC + dA( C0 + dC)
                   ! harvest+ clearance flux (not yet corrected for carbon remaining &
                   ! in landscape as litter)
                   !  [g C m-2] (grid-cell basis)
                   FDist(g) = -(patch(j+1)%frac *  (dcHarvClear(g) + dcExpand(g)+dcNat(g)) + &
                        (POPLUC%ptos(g)+POPLUC%gtos(g)-POPLUC%stog(g))* &
                        ( casapool%cplant(j+1,2) + dcHarvClear(g) + dcExpand(g)  ) )

                   Fnatdist(g) = -patch(j+1)%frac*dcNat(g)


                   ! partition disturbance flux between naturl dist, harvest and clearing
                   tmp = -patch(j+1)%frac*(dcNat(g) + dcHarv(g)) &
                        + POPLUC%stog(g)*POPLUC%CRelClear(g)*casapool%cplant(j+1,2)

                   Fnatdist(g) = -patch(j+1)%frac*dcNat(g)/tmp * FDist(g)
                   FHarv(g) = -patch(j+1)%frac*dcHarv(g)/tmp * FDist(g)
                   FClear(g) = POPLUC%stog(g)*POPLUC%CRelClear(g)*casapool%cplant(j+1,2)/ &
                        tmp * FDist(g)
                   FHarvClear(g) = FHarv(g)+FClear(g)

!!$                   if (g==4 ) write(*,*) 'fdist', &
!!$                   FDist(g),  Fnatdist(g), tmp, dcharv(g), FHarv(g), FClear(g), &
!!$                   POPLUC%CRelClear(g)
                   ! adjust biomass density changes to be consistent with FClear & FHarv


                   dcHarvClear(g) = -(FHarvCLear(g)+(-POPLUC%stog(g))&
                        *casapool%cplant(j+1,2))/ &
                        (patch(j+1)%frac + POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g))

                   ! POP%pop_grid(l)%cat_mortality = Fnatdist(g)/casapool%cplant(j+1,2)*&
                   !    (POP%pop_grid(l)%cmass_sum_old+POP%pop_grid(l)%growth)

                   ! secondary forest harvest flux (goes to harvest wood products) &
                   ! [g C m-2] (grid-cell basis)
                   !  corrected for carbon remaining in landscape as litter
                   if ((FHarv(g)).gt.0.0_dp) then
                      POPLUC%FHarvest(g,2) =  (1.0_dp-POPLUC%fracHarvSecResid(g))*FHarv(g)
                   endif
                   ! secondary forest clearance flux (goes to clearance pool) &
                   ! [g C m-2] (grid-cell basis)
                   !  corrected for carbon remaining in landscape as litter
                   if ((FClear(g)).gt.0.0_dp) then
                      POPLUC%FClearance(g,2) =  (1.0_dp-POPLUC%fracClearResid(g))* FClear(g)
                   endif
                else
                   dcHarvClear(g) = 0.0_dp
                   POPLUC%stog(g) = 0.0_dp
                endif
             endif
          ENDIF
          if (POP%Iwood(l).eq.j ) then
             ! set diagnostic POPLUC variables to POP equivalents &
             !(patch biomass and age distributions)
             POPLUC%biomass_age_primary(g,:)=POP%pop_grid(l)%biomass_age
             POPLUC%freq_age_primary(g,:)=POP%pop_grid(l)%freq_age
          endif
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
       dA_r       = 0.0_dp
       dA_d       = 0.0_dp
       dA         = 0.0_dp
       dcsoil     = 0.0_dp
       dcplant    = 0.0_dp
       dclitter   = 0.0_dp
       dclabile   = 0.0_dp
       dcsoil_r   = 0.0_dp
       dcplant_r  = 0.0_dp
       dclitter_r = 0.0_dp
       dcsoil_r   = 0.0_dp
       dcplant_r  = 0.0_dp
       dclitter_r = 0.0_dp
       dclabile_r = 0.0_dp

       dcsoil_d   = 0.0_dp
       dcplant_d  = 0.0_dp
       dclitter_d = 0.0_dp
       dcsoil_d   = 0.0_dp
       dcplant_d  = 0.0_dp
       dclitter_d = 0.0_dp
       dclabile_d = 0.0_dp

       dnsoil      = 0.0_dp
       dnplant     = 0.0_dp
       dnlitter    = 0.0_dp
       dnsoil_r    = 0.0_dp
       dnplant_r   = 0.0_dp
       dnlitter_r  = 0.0_dp
       dnsoil_r    = 0.0_dp
       dnplant_r   = 0.0_dp
       dnlitter_r  = 0.0_dp
       dnsoilmin_r = 0.0_dp

       dnlitter    = 0.0_dp
       dnsoil_d    = 0.0_dp
       dnplant_d   = 0.0_dp
       dnlitter_d  = 0.0_dp
       dnsoil_d    = 0.0_dp
       dnplant_d   = 0.0_dp
       dnlitter_d  = 0.0_dp
       dnsoilmin_d = 0.0_dp

       dpsoil         = 0.0_dp
       dpplant        = 0.0_dp
       dplitter       = 0.0_dp
       dpsoil_r       = 0.0_dp
       dpplant_r      = 0.0_dp
       dplitter_r     = 0.0_dp
       dpsoil_r       = 0.0_dp
       dpplant_r      = 0.0_dp
       dplitter_r     = 0.0_dp
       dwood_transfer = 0.0_dp

       dpsoil_d   = 0.0_dp
       dpplant_d  = 0.0_dp
       dplitter_d = 0.0_dp
       dpsoil_d   = 0.0_dp
       dpplant_d  = 0.0_dp
       dplitter_d = 0.0_dp

#ifdef __C13DEBUG__
       tmp_dplant = 0.0_dp
       tmp_tplant = 0.0_dp
       tmp_dlit   = 0.0_dp
       tmp_slit   = 0.0_dp
       tmp_dsoil = 0.0_dp
#endif

       IF (.NOT.LUC_EXPT%prim_only(g)) THEN ! only worry about transitions where land-use change is possible
          DO k = 1,nTrans  ! loop over all possible transition types
             if (k==1) then
                deltaA = POPLUC%ptos(g)        ! transition area: primary to seondary transition
                idp = j                        ! donor patch index
                irp = j+1                      ! receiver patch index
                idlu = p                       ! donor land use index (p = primary)
                irlu = s                       ! receiver land use index (s = secondary)
             elseif (k==2) then
                deltaA = POPLUC%ptog(g)        ! transition area: primary to grass (open) transition
                idp = j                        ! donor patch index
                irp = j+2                      ! receiver patch index
                idlu = p                       ! donor land use index (p = primary)
                irlu = gr                      ! receiver land use index (gr = grass or open)
             elseif (k==3) then
                deltaA = POPLUC%stog(g)        ! transition area: secondary to grass (open) transition
                idp = j+1                      ! donor patch index
                irp = j+2                      ! receiver patch index
                idlu = s                       ! donor land use index
                irlu = gr                       ! receiver land use index
             elseif(k==4) then
                deltaA = POPLUC%gtos(g)        ! transition area: grass (open) to secondary transition
                idp = j+2                      ! donor patch index
                irp = j+1                      ! receiver patch index
                idlu = gr                      ! donor land use index
                irlu = s                       ! receiver land use index (s = secondary)
             endif
             dcsoil   = 0.0_dp
             dclitter = 0.0_dp
             dcplant  = 0.0_dp
             dclabile = 0.0_dp
             dnsoil   = 0.0_dp
             dnlitter = 0.0_dp
             dnplant  = 0.0_dp
             dpsoil   = 0.0_dp
             dplitter = 0.0_dp
             dpplant  = 0.0_dp
             dwood_transfer = 0.0_dp

             ! transfer fluxes : only consider cases where gross transition area is greater than zero
             if (deltaA.gt.0.0_dp) then
                ! all soil pools

                ! change in carbon associated with gross transition
                dcsoil(irlu,idlu,:) = deltaA*casapool%csoil(idp,:)
#ifdef __C13DEBUG__
                tmp_dsoil(irlu,:) = tmp_dsoil(irlu,:) + deltaA*casapool%csoil(idp,:)
#endif
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
#ifdef __C13DEBUG__
                tmp_dlit(irlu,1) = tmp_dlit(irlu,1) + deltaA*casapool%clitter(idp,1)
#endif
                dclitter_r(irlu,1) =  dclitter_r(irlu,1) + dclitter(irlu,idlu,1)
                dclitter_d(idlu,1) =  dclitter_d(idlu,1) - dclitter(irlu,idlu,1)

                dnlitter(irlu,idlu,1) = deltaA*casapool%nlitter(idp,1)
                dnlitter_r(irlu,1) =  dnlitter_r(irlu,1) + dnlitter(irlu,idlu,1)
                dnlitter_d(idlu,1) =  dnlitter_d(idlu,1) - dnlitter(irlu,idlu,1)

                dplitter(irlu,idlu,1) = deltaA*casapool%plitter(idp,1)
                dplitter_r(irlu,1) =  dplitter_r(irlu,1) + dplitter(irlu,idlu,1)
                dplitter_d(idlu,1) =  dplitter_d(idlu,1) - dplitter(irlu,idlu,1)

                ! CWD : donor pool inhreits CWD and residues from harvest/clearance
                tmp = 0.0_dp
                if (idlu == s .and. casapool%cplant(idp,2).gt.1.e-5_dp ) then
                   ! secondary forest clearance: use secondary forest clearance flux from above
                   ! 13C
                   tmp = POPLUC%fracClearResid(g) * POPLUC%Fclearance(g,2) / (1.0_dp-POPLUC%fracClearResid(g))
                   popluc%FluxSClearResidtoLitter(g) = tmp
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3) + tmp
                   dnlitter(irlu,idlu,3) = deltaA*casapool%nlitter(idp,3) + &
                        POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.0_dp-POPLUC%fracClearResid(g))  &
                        *casapool%nplant(idp,2)/ casapool%cplant(idp,2)
                   dplitter(irlu,idlu,3) = deltaA*casapool%plitter(idp,3) + &
                       POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.0_dp-POPLUC%fracClearResid(g))  &
                          *casapool%pplant(idp,2)/ casapool%cplant(idp,2)
                elseif (idlu == p .and. irlu == s) then
                   ! primary forest harvest: assume harvest is uniform across age-classes
                   ! 13C
                   tmp = deltaA * POPLUC%fracHarvResid(g) * casapool%cplant(idp,2)
                   popluc%FluxPHarvResidtoLitter(g) = tmp
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3) + tmp
                   dnlitter(irlu,idlu,3) = deltaA*(casapool%nlitter(idp,3) + &
                        POPLUC%fracHarvResid(g)*casapool%nplant(idp,2) )
                   dplitter(irlu,idlu,3) = deltaA*(casapool%plitter(idp,3) + &
                        POPLUC%fracHarvResid(g) *casapool%pplant(idp,2) )
                elseif (idlu == p .and. irlu == gr) then
                   ! primary forest clearance: assume clearance is uniform across age-classes
                   ! 13C
                   tmp = deltaA * POPLUC%fracClearResid(g) * casapool%cplant(idp,2)
                   popluc%FluxPClearResidtoLitter(g) = tmp
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3) + tmp
                   dnlitter(irlu,idlu,3) = deltaA*(casapool%nlitter(idp,3) + &
                        POPLUC%fracClearResid(g)*casapool%nplant(idp,2) )
                   dplitter(irlu,idlu,3) = deltaA*(casapool%plitter(idp,3) + &
                       POPLUC%fracClearResid(g) *casapool%pplant(idp,2) )
                else
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3)
                   dnlitter(irlu,idlu,3) = deltaA*casapool%nlitter(idp,3)
                   dplitter(irlu,idlu,3) = deltaA*casapool%plitter(idp,3)
                endif
#ifdef __C13DEBUG__
                tmp_dlit(irlu,3) = tmp_dlit(irlu,3) + deltaA*casapool%clitter(idp,3)
                tmp_slit(irlu,3) = tmp_slit(irlu,3) + tmp
#endif
                dclitter_r(irlu,3) =  dclitter_r(irlu,3) + dclitter(irlu,idlu,3)
                dclitter_d(idlu,3) =  dclitter_d(idlu,3) - dclitter(irlu,idlu,3)

                dnlitter_r(irlu,3) =  dnlitter_r(irlu,3) + dnlitter(irlu,idlu,3)
                dnlitter_d(idlu,3) =  dnlitter_d(idlu,3) - dnlitter(irlu,idlu,3)

                dplitter_r(irlu,3) =  dplitter_r(irlu,3) + dplitter(irlu,idlu,3)
                dplitter_d(idlu,3) =  dplitter_d(idlu,3) - dplitter(irlu,idlu,3)

                ! fine structural litter: donor pool inherits leaves and fine roots
                dclitter(irlu,idlu,2) = deltaA*(casapool%clitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
#ifdef __C13DEBUG__
                tmp_dlit(irlu,2) = tmp_dlit(irlu,2) + deltaA*casapool%clitter(idp,2)
                tmp_slit(irlu,2) = tmp_slit(irlu,2) + deltaA*(casapool%cplant(idp,1) + casapool%cplant(idp,3))
#endif
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
                dcplant(irlu,idlu,:) = 0.0_dp
#ifdef __C13DEBUG__
                tmp_dplant(irlu,:) = tmp_dplant(irlu,:) + 0.0_dp
#endif
                dcplant_r(irlu,:) =  dcplant_r(irlu,:) + dcplant(irlu,idlu,:)
                dcplant_d(idlu,:) =  dcplant_d(idlu,:) - dcplant(irlu,idlu,:)

                dnplant(irlu,idlu,:) = 0.0_dp
                dnplant_r(irlu,:) =  dnplant_r(irlu,:) + dnplant(irlu,idlu,:)
                dnplant_d(idlu,:) =  dnplant_d(idlu,:) - dnplant(irlu,idlu,:)

                dpplant(irlu,idlu,:) = 0.0_dp
                dpplant_r(irlu,:) =  dpplant_r(irlu,:) + dpplant(irlu,idlu,:)
                dpplant_d(idlu,:) =  dpplant_d(idlu,:) - dpplant(irlu,idlu,:)

                ! Gross Transfer Flux (total C-tranfer associated with k-th transition)
                popluc%FtransferGross(g,k) = &
                     sum(dcsoil(irlu,idlu,:) + dclitter(irlu,idlu,:) +  dcplant(irlu,idlu,:))

                ! Harvest Flux
                ! primary to secondary forest
                ! (note secondary harvest and clearance fluxes are already evaluated &
                !  at the top of this subroutine)
                if (idlu==p .and. irlu ==s) then
                   popluc%FHarvest(g,idlu) =  (1.0_dp -POPLUC%fracHarvResid(g)) &
                        *casapool%cplant(idp,2)*deltaA
                endif
                ! Clearance Flux
                ! primary forest to grass
                if ((idlu==p) .and. irlu ==gr) then
                   popluc%FClearance(g,idlu) =  (1.0_dp -POPLUC%fracClearResid(g)) & 
                        * casapool%cplant(idp,2)*deltaA
                endif

                ! transition area
                dA_r(irlu) = dA_r(irlu) + deltaA
                dA_d(idlu) = dA_d(idlu) - deltaA
             endif
          ENDDO  ! ntrans


          do ilu=1, nlu
             ! update pools
             irp = ilu + j -1
#ifdef __C13DEBUG__
             iwtile = l-2
             iwpool = 1
#endif
             !if (g==3 .and. POPLUC%thisyear==1837 .and. ilu==2) &
             !     write(*,*) 'c00',casapool%cplant(irp,2),casapool%cplant(irp,2)*patch(j)%frac
             dwood_transfer = 0.0_dp
             dA(ilu) = dA_r(ilu) + dA_d(ilu)
             ! Net Transfer Flux
             POPLUC%FTransferNet(g,ilu) = sum(dcsoil_r(ilu,:) + dcsoil_d(ilu,:)) + &
                  sum(dclitter_r(ilu,:) + dclitter_d(ilu,:))  + &
                  sum(dcplant_r(ilu,:)  + dcplant_d(ilu,:)) + &
                  dclabile_r(ilu)  + dclabile_d(ilu)
             if (ilu ==2) then
                ! augment CWD pools in secondary veg tiles by harvest residues
                ! 13C
                tmp = POPLUC%fracHarvSecResid(g) * POPLUC%FHarvest(g,2) / (1.0_dp -POPLUC%fracHarvSecResid(g))
                popluc%FluxSHarvResidtoLitter(g) = tmp
                dclitter_r(ilu,3) = dclitter_r(irlu,3) + tmp
#ifdef __C13DEBUG__
                tmp_slit(irlu,3) = tmp_slit(irlu,3) + tmp
#endif
                if (casapool%cplant(irp,2) .gt. 1.e-5_dp) then
                   dnlitter_r(ilu,3) = dnlitter_r(irlu,3) +  &
                        POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0_dp -POPLUC%fracHarvSecResid(g)) &
                        *casapool%nplant(irp,2)/ casapool%cplant(irp,2)
                   dplitter_r(ilu,3) = dplitter_r(irlu,3) +  &
                        POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0_dp -POPLUC%fracHarvSecResid(g)) &
                        *casapool%pplant(irp,2)/ casapool%cplant(irp,2)
                endif
             endif

             if ((patch(irp)%frac+dA(ilu)).gt.1.e-5_dp) then ! avoid fpe's by ensuring finite new tile area
                casapool%nsoilmin(irp) = casapool%nsoilmin(irp) +  &
                     (dnsoilmin_r(ilu) - casapool%nsoilmin(irp)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     (dclabile_r(ilu) - casapool%clabile(irp)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%clabile(irp) = casapool%clabile(irp) +  &
                     (dclabile_r(ilu) - casapool%clabile(irp)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp)+ &
                     sum(dcsoil_r(ilu,:) - casapool%csoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

#ifdef __C13DEBUG__
                if (irp==iwtile) then
                   print*, 'PA01 ', patch(irp)%frac, patch(irp)%frac+dA(ilu), ilu
                   print*, 'PA02 ', dA(ilu) - dA_d(ilu)
                   print*, 'PA03 ', dA(ilu), dA_d(ilu), dA_r(ilu)
                endif
                if (irp==iwtile) then
                   print*, 'PS01 ', casapool%csoil(irp,iwpool)
                   print*, 'PS02 ', tmp_dsoil(ilu,iwpool)
                   print*, 'PS03 ', dcsoil_r(ilu,iwpool)
                   print*, 'PS04 ', casapool%csoil(irp,iwpool)*(dA(ilu) - dA_d(ilu))
                endif
#endif
                casapool%csoil(irp,:) = casapool%csoil(irp,:) + &
                     (dcsoil_r(ilu,:) - casapool%csoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))
#ifdef __C13DEBUG__
                if (irp==iwtile) then
                   print*, 'PS05 ', casapool%csoil(irp,iwpool)
                endif
#endif

                casapool%nsoil(irp,:) = casapool%nsoil(irp,:) + &
                     (dnsoil_r(ilu,:) - casapool%nsoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%psoil(irp,:) = casapool%psoil(irp,:) + &
                     (dpsoil_r(ilu,:) - casapool%psoil(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp) + &
                     sum(dclitter_r(ilu,:) - casapool%clitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

#ifdef __C13DEBUG__
                if (irp==iwtile) then
                   print*, 'PL01 ', casapool%clitter(irp,iwpool)
                   print*, 'PL02 ', tmp_dlit(ilu,iwpool)
                   print*, 'PL03 ', tmp_slit(ilu,iwpool)
                   print*, 'PL04 ', dclitter_r(ilu,iwpool)
                   print*, 'PL05 ', casapool%clitter(irp,iwpool)*(dA(ilu) - dA_d(ilu))
                endif
#endif
                casapool%clitter(irp,:) = casapool%clitter(irp,:) + &
                     (dclitter_r(ilu,:) - casapool%clitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))
#ifdef __C13DEBUG__
                if (irp==iwtile) then
                   print*, 'PL06 ', casapool%clitter(irp,iwpool)
                endif
#endif

                casapool%nlitter(irp,:) = casapool%nlitter(irp,:) + &
                     (dnlitter_r(ilu,:) - casapool%nlitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%plitter(irp,:) = casapool%plitter(irp,:) + &
                     (dplitter_r(ilu,:) - casapool%plitter(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casaflux%CtransferLUC(irp) = casaflux%CtransferLUC(irp) + &
                     sum(dcplant_r(ilu,:) - casapool%cplant(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

#ifdef __C13DEBUG__
                if (irp==iwtile) then
                   print*, 'PP01 ', casapool%cplant(irp,iwpool)
                   print*, 'PP02 ', tmp_dplant(ilu,iwpool)
                   print*, 'PP03 ', tmp_tplant(ilu,iwpool)
                   print*, 'PP04 ', dcplant_r(ilu,iwpool)
                   print*, 'PP05 ', casapool%cplant(irp,iwpool)*(dA(ilu) - dA_d(ilu))
                endif
#endif
                casapool%cplant(irp,:) = casapool%cplant(irp,:) + &
                     (dcplant_r(ilu,:) - casapool%cplant(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))
#ifdef __C13DEBUG__
                if (irp==iwtile) then
                   print*, 'PP06 ', casapool%cplant(irp,iwpool)
                endif
#endif

                casapool%nplant(irp,:) = casapool%nplant(irp,:) + &
                     (dnplant_r(ilu,:) - casapool%nplant(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,:) = casapool%pplant(irp,:) + &
                     (dpplant_r(ilu,:) - casapool%pplant(irp,:)*(dA(ilu) - dA_d(ilu))) &
                     /(patch(irp)%frac+dA(ilu))

                ! account here for change in secondary forest biomass density due to:
                !   harvest and clearing, as well as increase in (below ground) CWD
                !   where secondary forest harvest occurs

                if (ilu .eq.s .and. (casapool%cplant(irp,2)+dcHarvClear(g)).gt.0.0_dp  &
                     .and. casapool%cplant(irp,2).gt.1.e-10_dp  ) then
                   popluc%dcSHarvClear(g) = dcHarvClear(g)
                   casapool%nplant(irp,2) = casapool%nplant(irp,2) + (dcHarvClear(g)) &
                        * casapool%nplant(irp,2)/casapool%cplant(irp,2)
                   casapool%pplant(irp,2) = casapool%pplant(irp,2) + (dcHarvClear(g)) &
                        * casapool%pplant(irp,2)/casapool%cplant(irp,2)
                   casapool%cplant(irp,2) = casapool%cplant(irp,2) + dcHarvClear(g)
#ifdef __C13DEBUG__
                   if (irp==iwtile) then
                      print*, 'PP07 ', casapool%cplant(irp,iwpool), dcHarvClear(g)
                   endif
#endif
                elseif (ilu .eq.s .and. (casapool%cplant(irp,2)+dcHarvClear(g)).le.0.0_dp  ) then
                   POPLUC%FHarvest(g,2)   = 0.0_dp
                   POPLUC%FClearance(g,2) = 0.0_dp
                   popluc%FluxSHarvResidtoLitter(g)  = 0.0_dp
                   popluc%FluxSClearResidtoLitter(g) = 0.0_dp
                endif
             endif ! (patch(irp)%frac+dA(ilu)).gt.1.e-5_dp

             if (patch(irp)%frac .gt. 1.e-8_dp) then
                casaflux%FluxCtohwp(irp)   = POPLUC%FHarvest(g,ilu)   / patch(irp)%frac * real(ktauday,dp)
                casaflux%FluxCtoclear(irp) = POPLUC%FClearance(g,ilu) / patch(irp)%frac * real(ktauday,dp)
             endif
             !if (g==3 .and. POPLUC%thisyear==1837 ) write(*,*) 'c01c',  casapool%cplant(irp,2), ilu,dcharvclear(g)
          enddo ! ilu=1, nlu
          
          ! POPLUC diagnostics
          ! pools in gC per m2 of gridcell
          ! NEP in g C y-1 per m2 of gridcell
          POPLUC%FNEP(g,:) = casabal%Fcneeyear(j:l)*patch(j:l)%frac  ! note NEE = NEP here
          !(+ve into surface)
          DO ilu=1,nLU
             ! update area weights
             irp = ilu + j -1
             ! Update tile area
             if (ilu == p) then
                POPLUC%primf(g)  = max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp)
             elseif (ilu == s) then
                POPLUC%secdf(g) = max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp)
                if (POPLUC%secdf(g) .eq. 0.0_dp) then
                   POPLUC%freq_age_secondary(g,:) = 0.0_dp
                endif
             elseif (ilu == gr) then
                POPLUC%grass(g) = max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp)
                POPLUC%past(g) = min(max(POPLUC%past(g) + POPLUC%ptoq(g) + POPLUC%stoq(g) &
                     - POPLUC%qtos(g),0.0_dp), POPLUC%grass(g))
                POPLUC%crop(g) =min( max(POPLUC%crop(g) + POPLUC%ptoc(g) + POPLUC%stoc(g) &
                     - POPLUC%ctos(g),0.0_dp), POPLUC%grass(g) - POPLUC%past(g))
             endif

             POPLUC%csoil(g,ilu) = sum(casapool%csoil(irp,:))* &
                  max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp)
             POPLUC%clitt(g,ilu) = sum(casapool%clitter(irp,:))* &
                  max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp)
             POPLUC%cbiomass(g,ilu) = sum(casapool%cplant(irp,:))* &
                  max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp)
             ! if (g==4  .and. ilu==2) write(*,*) 'c02', casapool%cplant(irp,2), &
             !    casapool%cplant(irp,2)* max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp), &
             !    patch(irp)%frac,  dA_r(ilu) + dA_d(ilu)

             !  if (g==4  .and. ilu==3) write(*,*) 'c03', casapool%csoil(irp,2), &
             !     casapool%csoil(irp,2)* max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0_dp), &
             !     patch(irp)%frac,  dA_r(ilu) + dA_d(ilu)

          ENDDO
       ELSE

          POPLUC%csoil(g,:) = 0.0_dp
          POPLUC%clitt(g,:) = 0.0_dp
          POPLUC%cbiomass(g,:) = 0.0_dp
          POPLUC%FNEP(g,:) = 0.0_dp

          POPLUC%csoil(g,1) = sum(casapool%csoil(j,:))*patch(j)%frac
          POPLUC%clitt(g,1) = sum(casapool%clitter(j,:))*patch(j)%frac
          POPLUC%cbiomass(g,1) = sum(casapool%cplant(j,:))*patch(j)%frac
          POPLUC%FNEP(g,1) =  casabal%Fcneeyear(j)*patch(j)%frac
          POPLUC%primf(g) = patch(j)%frac
          POPLUC%secdf(g) = 0.0_dp
          if (POPLUC%grass(g).gt.0.0_dp .and. l.eq.j+1) then
             POPLUC%csoil(g,3) = sum(casapool%csoil(l,:))*patch(l)%frac
             POPLUC%clitt(g,3) = sum(casapool%clitter(l,:))*patch(l)%frac
             POPLUC%cbiomass(g,3) = sum(casapool%cplant(l,:))*patch(l)%frac
             POPLUC%FNEP(g,3) = casabal%Fcneeyear(l)*patch(l)%frac
             POPLUC%grass(g) = patch(l)%frac
          endif

       ENDIF

       POPLUC%HarvProdLoss(g,:) = kHarvProd * POPLUC%HarvProd(g,:)
       POPLUC%ClearProdLoss(g,:) = kClearProd * POPLUC%ClearProd(g,:)
       POPLUC%AgProdLoss(g) = kAgProd*POPLUC%AgProd(g)

       if (POPLUC%grass(g).gt.0.0_dp .and. l.eq.j+2) then
       POPLUC%FAg(g) =  casaflux%Charvest(l)*patch(l)%frac
       POPLUC%AgProd(g) = POPLUC%AgProd(g) + casaflux%Charvest(l)*patch(l)%frac - POPLUC%AgProdLoss(g)
       casaflux%charvest(l) = 0.0_dp
       casaflux%nharvest(l) = 0.0_dp
       casaflux%fharvest(l) = min(POPLUC%past(g)/POPLUC%grass(g),1.0_dp)*0.5_dp + &
            min(POPLUC%crop(g)/POPLUC%grass(g),1.0_dp)*0.9_dp !  fraction grass AGB to be removed next year
!write(*,*)'harvest', g,   patch(l)%frac,  casaflux%fharvest(l), POPLUC%crop(g),  POPLUC%past(g)
       !casaflux%fharvest(l) = 0; ! test vh!
        casaflux%fcrop(l) = min(POPLUC%crop(g)/POPLUC%grass(g), 1.0_dp)
       endif
       DO j=1, 3
          POPLUC%HarvProd(g,j) = POPLUC%HarvProd(g,j) + &
               POPLUC%fracHarvProd(g,j)*sum(POPLUC%FHarvest(g,:)) - POPLUC%HarvProdLoss(g,j)

          POPLUC%ClearProd(g,j) = POPLUC%ClearProd(g,j) + &
               POPLUC%fracClearProd(g,j)*sum(POPLUC%FClearance(g,:)) - POPLUC%ClearProdLoss(g,j)
       ENDDO

    ENDDO

991 format(1166(e14.7,2x))

! update total carbon pools and "last" pool values for use in carbon balance checks.
    casapool%ctot = sum(casapool%cplant,2)+sum(casapool%clitter,2)+ &
         sum(casapool%csoil,2)+casapool%clabile
    casabal%cplantlast  = casapool%cplant
    casabal%clabilelast = casapool%clabile
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil

  END SUBROUTINE POP_LUC_CASA_transfer

  !*******************************************************************************

  SUBROUTINE POPLUC_Init(POPLUC, LUC_EXPT, casapool, casaflux, casabiome, veg, POP, np)

    USE cable_def_types_mod, ONLY : veg_parameter_type
    USE casaparm, ONLY: LEAF, WOOD, FROOT

    IMPLICIT NONE

    TYPE(POPLUC_TYPE),        INTENT(INOUT) :: POPLUC
    TYPE(LUC_EXPT_TYPE),      INTENT(IN)    :: LUC_EXPT
    TYPE(casa_pool),          INTENT(INOUT) :: casapool
    TYPE(casa_flux),          INTENT(INOUT) :: casaflux
    TYPE(casa_biome),         INTENT(IN)    :: casabiome
    TYPE(veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
    TYPE(POP_TYPE),           INTENT(INOUT) :: POP

    INTEGER(i4b), INTENT(IN) :: np
    INTEGER(i4b) :: g, j, k, l

    CALL alloc_POPLUC(POPLUC,np)

    POPLUC%it = 0
    POPLUC%np = np

    CALL ZeroPOPLUC(POPLUC)

    CALL POPLUC_set_params(POPLUC, LUC_EXPT)

    IF (cable_user%POPLUC_RunType .eq. 'init') THEN
       POPLUC%frac_primf = LUC_EXPT%primaryf
       POPLUC%primf      = LUC_EXPT%primaryf
       POPLUC%grass      = LUC_EXPT%grass
       POPLUC%past       = LUC_EXPT%past
       POPLUC%crop       = LUC_EXPT%crop
       where ((POPLUC%primf + POPLUC%grass) > 1.0_dp) POPLUC%grass = 1.0_dp - POPLUC%primf
       POPLUC%frac_forest = 1.0_dp - POPLUC%grass
       POPLUC%freq_age_secondary(:,1) = max(POPLUC%frac_forest - POPLUC%primf, 0.0_dp)
       POPLUC%latitude  = real(patch(landpt(:)%cstart)%latitude,dp)
       POPLUC%longitude = real(patch(landpt(:)%cstart)%longitude,dp)

       ! zero biomass in secondary forest tiles (both CASA and POP variables)
       where (veg%iLU ==2)
          casapool%cplant(:,leaf) = 0.01_dp
          casapool%nplant(:,leaf)= casabiome%ratioNCplantmin(veg%iveg,leaf)* casapool%cplant(:,leaf)
          casapool%pplant(:,leaf)= casabiome%ratioPCplantmin(veg%iveg,leaf)* casapool%cplant(:,leaf)

          casapool%cplant(:,froot) = 0.01_dp
          casapool%nplant(:,froot)= casabiome%ratioNCplantmin(veg%iveg,froot)* &
               casapool%cplant(:,froot)
          casapool%pplant(:,froot)= casabiome%ratioPCplantmin(veg%iveg,froot)* &
               casapool%cplant(:,froot)

          casapool%cplant(:,wood) = 0.01_dp
          casapool%nplant(:,wood)= casabiome%ratioNCplantmin(veg%iveg,wood)* casapool%cplant(:,wood)
          casapool%pplant(:,wood)= casabiome%ratioPCplantmin(veg%iveg,wood)* casapool%cplant(:,wood)
          casaflux%frac_sapwood = 1.0_dp
       endwhere

       DO k=1, np
          IF (.NOT. LUC_EXPT%prim_only(k)) THEN
             j = landpt(k)%cstart+1
             do l=1, size(POP%Iwood)
                if (POP%Iwood(l) == j) then
                   CALL POP_init_single(POP,veg%disturbance_interval,l)
                   exit
                endif
             enddo
          ENDIF
       ENDDO

    ELSEIF (cable_user%POPLUC_RunType .eq. 'restart') THEN

       CALL READ_LUC_RESTART_NC(POPLUC)

       POPLUC%frac_primf  = POPLUC%primf
       POPLUC%frac_forest = 1.0_dp - POPLUC%grass
       POPLUC%latitude    = real(patch(landpt(:)%cstart)%latitude)
       POPLUC%longitude   = real(patch(landpt(:)%cstart)%longitude)

    ENDIF

    ! set landuse index for secondary forest POP landscapes
    ! (not for 'static') run type because secondary forest dynamics are only
    ! simulated with dynamic land-use forcing
    IF (TRIM(cable_user%POPLUC_RunType) .ne. 'static') THEN
       DO k=1,POP%np
          if (veg%iLU(POP%Iwood(k)).eq.2) then
             POP%pop_grid(k)%LU = 2
          endif
       ENDDO
    ENDIF


  END SUBROUTINE POPLUC_Init

  !*******************************************************************************

  SUBROUTINE POPLUC_set_patchfrac(POPLUC,LUC_EXPT)

    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    INTEGER(i4b) :: j, k, l

    DO k=1, POPLUC%np
       j = landpt(k)%cstart
       l = landpt(k)%cend
       IF (.NOT.LUC_EXPT%prim_only(k)) THEN
          patch(j)%frac = POPLUC%primf(k)
          patch(l)%frac = POPLUC%grass(k)
          patch(j+1)%frac = 1.0_dp -  patch(j)%frac - patch(l)%frac
       ENDIF
    ENDDO

  END SUBROUTINE POPLUC_SET_PATCHFRAC

  !*******************************************************************************

  SUBROUTINE POPLUC_set_params(POPLUC,LUC_EXPT)

    IMPLICIT NONE

    TYPE(POPLUC_TYPE),    INTENT(INOUT) :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(IN)    :: LUC_EXPT

    INTEGER(i4b) :: g, np

    np = POPLUC%np
    DO g=1, np
       POPLUC%fracharvProd(g, 1) = 0.9_dp
       POPLUC%fracharvProd(g, 2) = 0.04_dp
       POPLUC%fracharvProd(g, 3) = 0.06_dp

       POPLUC%fracClearProd(g, 1) = 0.597_dp
       POPLUC%fracClearProd(g, 2) = 0.403_dp
       POPLUC%fracClearProd(g, 3) = 0.0_dp
       POPLUC%fracClearResid(g) = 0.33_dp
       POPLUC%fracHarvResid(g) = 0.79_dp
       POPLUC%fracHarvSecResid(g) = 0.81_dp

       IF (LUC_EXPT%biome(g)==1 .OR. LUC_EXPT%biome(g)==2) THEN
          ! Tropical Evergreen and Tropical Deciduous
          POPLUC%fracharvProd(g, 1) = 0.9_dp
          POPLUC%fracharvProd(g, 2) = 0.04_dp
          POPLUC%fracharvProd(g, 3) = 0.06_dp

          POPLUC%fracClearProd(g, 1) = 0.597_dp
          POPLUC%fracClearProd(g, 2) = 0.403_dp
          POPLUC%fracClearProd(g, 3) = 0.0_dp

          IF (LUC_EXPT%biome(g)==1) POPLUC%fracHarvResid(g) = 0.79_dp
          IF (LUC_EXPT%biome(g)==2) POPLUC%fracHarvResid(g) = 0.86_dp
          IF (LUC_EXPT%biome(g)==1) POPLUC%fracHarvSecResid(g) = 0.71_dp
          IF (LUC_EXPT%biome(g)==2) POPLUC%fracHarvSecResid(g) = 0.81_dp

          POPLUC%fracClearResid(g) = 0.33_dp

       ELSEIF (LUC_EXPT%biome(g).GE.4 .OR. LUC_EXPT%biome(g).LE.10) THEN

          ! Other Forest
          POPLUC%fracharvProd(g, 1) = 0.4_dp
          POPLUC%fracharvProd(g, 2) = 0.24_dp
          POPLUC%fracharvProd(g, 3) = 0.36_dp

          POPLUC%fracClearProd(g, 1) = 0.597_dp
          POPLUC%fracClearProd(g, 2) = 0.2985_dp
          POPLUC%fracClearProd(g, 3) = 0.1045_dp

          IF (LUC_EXPT%ivegp(g)==2) POPLUC%fracHarvResid(g) = 0.83_dp
          IF (LUC_EXPT%ivegp(g)==1 .OR. LUC_EXPT%ivegp(g)==3 ) POPLUC%fracHarvResid(g) = 0.87_dp
          IF (LUC_EXPT%ivegp(g)==4) POPLUC%fracHarvResid(g) = 0.78_dp

          IF (LUC_EXPT%ivegp(g)==2) POPLUC%fracHarvSecResid(g) = 0.75_dp
          IF (LUC_EXPT%ivegp(g)==1 .OR. LUC_EXPT%ivegp(g)==3 ) POPLUC%fracHarvSecResid(g) = 0.82_dp
          IF (LUC_EXPT%ivegp(g)==4) POPLUC%fracHarvSecResid(g) = 0.70_dp

          POPLUC%fracClearResid(g) = 0.33_dp

       ELSEIF (LUC_EXPT%biome(g).EQ.3 .OR. LUC_EXPT%biome(g).GE.11) THEN
          ! savanna and shrub
          POPLUC%fracharvProd(g, 1) = 1.0_dp
          POPLUC%fracharvProd(g, 2) = 0.0_dp
          POPLUC%fracharvProd(g, 3) = 0.0_dp

          POPLUC%fracClearProd(g, 1) = 0.8_dp
          POPLUC%fracClearProd(g, 2) = 0.2_dp
          POPLUC%fracClearProd(g, 3) = 0.0_dp

          IF (LUC_EXPT%biome(g).EQ.13 .OR.(LUC_EXPT%biome(g).EQ.14)  ) THEN
             POPLUC%fracHarvResid(g) = 0.78_dp
             POPLUC%fracHarvSecResid(g) = 0.70_dp
          ELSE
             POPLUC%fracHarvResid(g) = 0.86_dp
             POPLUC%fracHarvSecResid(g) = 0.81_dp
          ENDIF

          POPLUC%fracClearResid(g) = 0.5_dp

       ENDIF
      ! no residue test
      !POPLUC%fracClearResid(g) = 0.0_dp
     ! POPLUC%fracHarvResid(g) = 0.0_dp
     ! POPLUC%fracHarvSecResid(g) = 0.0_dp

    ENDDO

  END SUBROUTINE POPLUC_SET_PARAMS

  !*******************************************************************************

  SUBROUTINE alloc_POPLUC(POPLUC, arraysize)

    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER,           INTENT(IN)    :: arraysize

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
    ALLOCATE(POPLUC%crop(arraysize))
    ALLOCATE(POPLUC%past(arraysize))
    ALLOCATE(POPLUC%ptos(arraysize))
    ALLOCATE(POPLUC%ptog(arraysize))
    ALLOCATE(POPLUC%stog(arraysize))
    ALLOCATE(POPLUC%gtop(arraysize))
    ALLOCATE(POPLUC%gtos(arraysize))
    ALLOCATE(POPLUC%ptog(arraysize))
    ALLOCATE(POPLUC%ptoc(arraysize))
    ALLOCATE(POPLUC%ptoq(arraysize))
    ALLOCATE(POPLUC%stoc(arraysize))
    ALLOCATE(POPLUC%stoq(arraysize))
    ALLOCATE(POPLUC%ctos(arraysize))
    ALLOCATE(POPLUC%qtos(arraysize))
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
    ALLOCATE(POPLUC%AgProd(arraysize))
    ALLOCATE(POPLUC%ClearProd(arraysize,3))
    ALLOCATE(POPLUC%HarvProdLoss(arraysize,3))
    ALLOCATE(POPLUC%ClearProdLoss(arraysize,3))
    ALLOCATE(POPLUC%AgProdLoss(arraysize))
    ALLOCATE(POPLUC%FAg(arraysize))
    ALLOCATE(POPLUC%fracHarvProd(arraysize,3))
    ALLOCATE(POPLUC%fracClearProd(arraysize,3))
    ALLOCATE(POPLUC%fracHarvResid(arraysize))
    ALLOCATE(POPLUC%fracHarvSecResid(arraysize))
    ALLOCATE(POPLUC%fracClearResid(arraysize))
    ALLOCATE(POPLUC%kSecHarv(arraysize))
    ALLOCATE(POPLUC%kNatDist(arraysize))
    ALLOCATE(POPLUC%kExpand1(arraysize))
    ALLOCATE(POPLUC%kExpand2(arraysize))
    ALLOCATE(POPLUC%kClear(arraysize))
    ALLOCATE(POPLUC%cRelClear(arraysize))
    allocate(popluc%FluxPHarvResidtoLitter(arraysize))
    allocate(popluc%FluxSHarvResidtoLitter(arraysize))
    allocate(popluc%FluxPClearResidtoLitter(arraysize)) ! Residual flux to litter from clearing primary forest
    allocate(popluc%FluxSClearResidtoLitter(arraysize))
    allocate(popluc%dcSHarvClear(arraysize))

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

    REAL(dp), INTENT(IN) :: lambda
    REAL(dp), INTENT(IN) :: x

    IF (x.LT.0.0_dp) THEN ! Shouldn't happen but ...
       REALExponential = 0.0_dp
    ELSE
       REALExponential = lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION REALExponential

  !*******************************************************************************

  SUBROUTINE WRITE_LUC_OUTPUT_NC( POPLUC, ctime, FINAL )

    USE CABLE_COMMON_MODULE, ONLY: filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    LOGICAL, INTENT(IN)    :: FINAL
    INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, hist_ID, t_ID, nLU_ID, nTrans_ID
    INTEGER   :: i, mp, nprod, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99, dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp,t)
    CHARACTER(len=20),DIMENSION(24):: A1
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
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)), VID6(size(A6))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    CHARACTER(len=50) :: RecordDimName
    REAL, ALLOCATABLE :: freq_age_secondary(:,:)
    INTEGER :: g
    LOGICAL :: put_age_vars
    mp = POPLUC%np
    nprod = 3
    put_age_vars=.TRUE.
    allocate(freq_age_secondary(mp,age_max))

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
    A1(14) = 'crop'
    A1(15) = 'past'
    A1(16) = 'ptoc'
    A1(17) = 'ptoq'
    A1(18) = 'stoc'
    A1(19) = 'stoq'
    A1(20) = 'qtos'
    A1(21) = 'ctos'
    A1(22) = 'AgProd'
    A1(23) = 'AgProdLoss'
    A1(24) = 'FAg'

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
       if (sum(POPLUC%freq_age_secondary(g,:)).gt.1.e-12_dp ) then
          freq_age_secondary(g,:) = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
       else
          freq_age_secondary(g,:) = 0.0_dp
       endif
    ENDDO

    CNT = CNT + 1

    ! Get File-Name
    if (cable_user%yearstart < 1000) then
       write(dum, fmt="(i3)") cable_user%yearstart
    else
       write(dum, fmt="(i4)") cable_user%yearstart
    endif
    if (cable_user%yearend < 1000) then
       write(dum, fmt="(a,a,i3)") trim(dum), '_', cable_user%yearend
    else
       write(dum, fmt="(a,a,i4)") trim(dum), '_', cable_user%yearend
    endif

    if (len_trim(cable_user%LUC_outfile) .gt. 0) then
       fname = trim(cable_user%LUC_outfile)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_'//trim(dum)//'_LUC_out.nc'
    endif
    
    IF ( CALL1 ) THEN

       ! Create NetCDF file:
       STATUS = NF90_create(trim(fname), ior(nf90_clobber,nf90_64bit_offset), FILE_ID)
       ! print*, 'OCreate65 ', file_id, trim(fname)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

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
       if(put_age_vars) then
          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,age_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO
       endif
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
       STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), real(POPLUC%latitude,sp))
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), real(POPLUC%longitude,sp))
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       CALL1 = .FALSE.
       ! print*, 'OCreated65'

    ENDIF ! CALL1



    ! TIME  ( t )
    ! print*, 'OWrite65 ', file_id
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), real(POPLUC%primf,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), real(POPLUC%secdf,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), real(POPLUC%grass,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), real(POPLUC%ptos,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), real(POPLUC%ptog,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), real(POPLUC%stog,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), real(POPLUC%gtop,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), real(POPLUC%gtos,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), real(POPLUC%frac_primf,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), real(POPLUC%frac_forest,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 11), real(POPLUC%pharv,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 12), real(POPLUC%smharv,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 13), real(POPLUC%syharv,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 14), real(POPLUC%crop,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 15), real(POPLUC%past,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 16), real(POPLUC%ptoc,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 17), real(POPLUC%ptoq,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 18), real(POPLUC%stoc,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 19), real(POPLUC%stoq,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 20), real(POPLUC%qtos,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 21), real(POPLUC%ctos,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 22), real(POPLUC%AgProd,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 23), real(POPLUC%AgProdLoss,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 24), real(POPLUC%FAg,sp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), POPLUC%n_event, start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    if (put_age_vars) then
       ! PUT 3D VARS ( mp, mage, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), real(POPLUC%freq_age_primary,sp), &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2),freq_age_secondary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), real(POPLUC%biomass_age_primary,sp), &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), real(POPLUC%biomass_age_secondary,sp), &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    endif

    ! PUT 3D VARS ( mp, nLU, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), real(POPLUC%FHarvest,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), real(POPLUC%FClearance,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), real(POPLUC%FNEP,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), real(POPLUC%Clitt,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), real(POPLUC%CSoil,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), real(POPLUC%CBiomass,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), real(POPLUC%FTransferNet,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    ! PUT 3D VARS ( mp, nTrans, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), real(POPLUC%FTransferGross,sp),   &
         start=(/ 1,1,CNT /), count=(/ mp,nTrans,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, nprod, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), real(POPLUC%HarvProd,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(2), real(POPLUC%ClearProd,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(3), real(POPLUC%HarvProdLoss,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(4), real(POPLUC%ClearProdLoss,sp), &
         start=(/ 1,1,CNT /), count=(/ mp,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    IF ( FINAL ) THEN
       ! Close NetCDF file:
       ! print*, 'OClose65 ', file_id
       STATUS = NF90_close(FILE_ID)
       file_id = -1
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " POPLUC Output written to ", trim(fname)
    ENDIF

  END SUBROUTINE WRITE_LUC_OUTPUT_NC

  !*********************************************************************************************************************

  SUBROUTINE WRITE_LUC_RESTART_NC(POPLUC, ctime)

    use cable_common_module, only: filename, cable_user, HANDLE_ERR
    use netcdf

    implicit none

    type(popluc_type), intent(in) :: POPLUC
    integer,           intent(in) :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, nLU_ID, nTrans_ID, i, mp, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99, dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20), DIMENSION(2) :: A0
    ! 2 dim real arrays (mp)
    CHARACTER(len=20), DIMENSION(6) :: A1
    ! 2 dim real arrays (mp,age_max)
    CHARACTER(len=25), DIMENSION(2) :: A2
    ! 2 dim real arrays (mp,nprod)
    CHARACTER(len=25), DIMENSION(2) :: A3

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)), VID1(SIZE(A1))
    INTEGER, SAVE :: VID2(SIZE(A2)), VID3(SIZE(A3))
    INTEGER, SAVE :: FILE_ID, CNT=0
    CHARACTER(len=50) :: RecordDimName
    INTEGER :: g, nprod

    mp = POPLUC%np
    nprod = 3 ! number of product pools

    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'
    A1(4) = 'crop'
    A1(5) = 'past'
    A1(6) = 'AgProd'

    A2(1) = 'biomass_age_secondary'
    A2(2) = 'freq_age_secondary'

    A3(1) = 'HarvProd'
    A3(2) = 'ClearProd'

    ! Get File-Name

    write(dum, fmt="(i4,'_',i4)") cable_user%yearstart, cable_user%yearend
    if (cable_user%yearstart.lt.1000 .and. cable_user%yearend.lt.1000) then
       write(dum, fmt="(i3,'_',i3)") cable_user%yearstart, cable_user%yearend
    elseif (cable_user%yearstart.lt.1000) then
       write(dum, fmt="(i3,'_',i4)") cable_user%yearstart, cable_user%yearend
    endif

    if (len_trim(cable_user%luc_restart_out) .gt. 0) then
       fname = trim(cable_user%luc_restart_out)
    else
       fname = trim(filename%path)//'/'//trim(cable_user%RunIden)//'_'//'LUC_rst.nc'
    endif

    ! Create NetCDF file:
    STATUS = NF90_create(trim(fname), ior(nf90_clobber,nf90_64bit_offset), FILE_ID)
    ! print*, 'OCreate66 ', file_id, trim(fname)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

    ! Define dimensions:
    ! Land (number of points)
    status = nf90_def_dim(file_id, 'land',  mp,      land_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(file_id, 'mage',  age_max, age_id)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_dim(file_id, 'nprod', nprod,   nprod_id)

    ! Define variables
    DO i = 1, SIZE(A0)
       STATUS = NF90_def_var(FILE_ID, TRIM(A0(i)), NF90_DOUBLE, (/land_ID/), VID0(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID, TRIM(A1(i)), NF90_DOUBLE, (/land_ID/), VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID, TRIM(A2(i)), NF90_DOUBLE, (/land_ID,age_ID/), VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A3)
       STATUS = NF90_def_var(FILE_ID, TRIM(A3(i)), NF90_DOUBLE, (/land_ID,nprod_ID/), VID3(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT LAT / LON ( mp )
    STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), POPLUC%latitude)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), POPLUC%longitude)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), POPLUC%primf)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), POPLUC%secdf)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), POPLUC%grass)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), POPLUC%crop)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), POPLUC%past)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), POPLUC%AgProd)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, mage, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), POPLUC%biomass_age_secondary)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2),POPLUC%freq_age_secondary)

    ! PUT 3D VARS ( mp, nprod, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), POPLUC%HarvProd)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID3(2),POPLUC%ClearProd)


    ! Close NetCDF file:
    ! print*, 'OClose66 ', file_id
    STATUS = NF90_close(FILE_ID)
    file_id = -1
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    WRITE(*,*) " POPLUC Restart written to ", trim(fname)

  END SUBROUTINE WRITE_LUC_RESTART_NC

  !*******************************************************************************

  SUBROUTINE READ_LUC_RESTART_NC(POPLUC)

    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC

    INTEGER   :: STATUS, land_dim, mage_dim, nprod_dim
    INTEGER   :: land_ID, age_ID, nLU_ID, nTrans_ID, i, mp, FILE_ID, dID, nprod, nprod_ID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.
    REAL(dp), ALLOCATABLE :: TMP(:), TMP2(:,:), TMP3(:,:)

    ! 1 dim arrays (mp )
    CHARACTER(len=20), DIMENSION(2) :: A0
    ! 2 dim real arrays (mp)
    CHARACTER(len=20), DIMENSION(6) :: A1
    ! 2 dim real arrays (mp,age_max)
    CHARACTER(len=25), DIMENSION(2) :: A2
    ! 2 dim real arrays (mp,nprod)
    CHARACTER(len=25), DIMENSION(2) :: A3

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
    A1(4) = 'crop'
    A1(5) = 'past'
    A1(6) = 'AgProd'

    A2(1) = 'biomass_age_secondary'
    A2(2) = 'freq_age_secondary'

    A3(1) = 'HarvProd'
    A3(2) = 'ClearProd'

    IF (LEN_TRIM(cable_user%LUC_restart_in) .gt. 0) THEN
       fname = TRIM(cable_user%LUC_restart_in)
    ELSE
       fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//'LUC_rst.nc'
    ENDIF
    STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
    ! print*, 'OOpen67 ', file_id, trim(fname)
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

       SELECT CASE (TRIM(A1(i)))
       CASE ('primf') ;  POPLUC%primf  = TMP
       CASE ('secdf') ;  POPLUC%secdf  = TMP
       CASE ('grass') ;  POPLUC%grass  = TMP
       CASE ('crop') ;   POPLUC%crop   = TMP
       CASE ('past') ;   POPLUC%past   = TMP
       CASE ('AgProd') ; POPLUC%AgProd = TMP
       END SELECT
    END DO

 991  format(1000(e12.4,2x))
    ! READ 2-dimensional fields (mage)
    DO i = 1, SIZE(A2)
       STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A2(i)))
       CASE ('biomass_age_secondary') ; POPLUC%biomass_age_secondary = TMP2
       CASE ('freq_age_secondary' ) ;   POPLUC%freq_age_secondary    = TMP2
       END SELECT
    END DO

    ! READ 3-dimensional fields (nprod)
    DO i = 1, SIZE(A3)
       STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A3(i)))
       CASE ('HarvProd') ;  POPLUC%HarvProd  = TMP3
       CASE ('ClearProd') ; POPLUC%ClearProd = TMP3
       END SELECT
    END DO

    ! print*, 'OClose67 ', file_id
    STATUS = NF90_CLOSE( FILE_ID )
    file_id = -1

  END SUBROUTINE READ_LUC_RESTART_NC


  !*******************************************************************************

  SUBROUTINE WRITE_LUC_OUTPUT_GRID_NC( POPLUC, ctime, FINAL )

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
    CHARACTER(len=20),DIMENSION(24):: A1
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
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)), VID6(size(A6))
    INTEGER, SAVE :: FILE_ID, CNT = 0, latID, lonID
    CHARACTER(len=50) :: RecordDimName
    REAL, ALLOCATABLE :: freq_age_secondary(:,:)
    INTEGER :: g, k
    LOGICAL :: put_age_vars
    REAL(dp) :: ncmissingr = -1.e+33_dp
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
    allocate(freq_age_secondary(mp,age_max))
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
    A1(14) = 'crop'
    A1(15) = 'past'
    A1(16) = 'ptoc'
    A1(17) = 'ptoq'
    A1(18) = 'stoc'
    A1(19) = 'stoq'
    A1(20) = 'qtos'
    A1(21) = 'ctos'
    A1(22) = 'AgProd'
    A1(23) = 'AgProdLoss'
    A1(24) = 'FAg'


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
       if (sum(POPLUC%freq_age_secondary(g,:)).gt.1.e-12_dp ) then
          freq_age_secondary(g,:) = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
       else
          freq_age_secondary(g,:) = 0.0_dp
       endif
    ENDDO

    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       IF (CABLE_USER%YEARSTART.lt.1000.and.CABLE_USER%YEAREND.lt.1000) THEN
          WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ELSEIF (CABLE_USER%YEARSTART.lt.1000) THEN
          WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ENDIF
       fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
            TRIM(dum)//'_LUC_out.nc'

       ! Create NetCDF file:
       STATUS = NF90_create(trim(fname), ior(nf90_clobber,nf90_64bit_offset), FILE_ID)
       ! print*, 'OCreate68 ', file_id, trim(fname)
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

       if(put_age_vars) then
          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/xID, yID,age_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO
       endif

       DO i = 1, SIZE(A4)
          STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/xID, yID,nLU_ID,t_ID/),VID4(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID4(i), '_FillValue', REAL(ncmissingr, sp))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A4(i)//' variable attributes in output file. '// &
                                                      '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID4(i), 'missing_value', REAL(ncmissingr, sp))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A4(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
       END DO

       DO i = 1, SIZE(A5)
          STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT,(/xID, yID,ntrans_ID,t_ID/),VID5(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID5(i), '_FillValue', REAL(ncmissingr, sp))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A5(i)//' variable attributes in output file. '// &
                                                      '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID5(i), 'missing_value', REAL(ncmissingr, sp))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A5(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
       END DO

       DO i = 1, SIZE(A6)
          STATUS = NF90_def_var(FILE_ID,TRIM(A6(i)) ,NF90_FLOAT,(/xID, yID ,nprod_ID,t_ID/),VID6(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID6(i), '_FillValue', REAL(ncmissingr, sp))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A6(i)//' variable attributes in output file. '// &
                                                      '(INTERFACE define_ovar)')
          ! Define missing/fill values:
          STATUS = NF90_PUT_ATT(FILE_ID, VID6(i), 'missing_value', REAL(ncmissingr, sp))
          IF (STATUS /= NF90_NOERR) CALL handle_err                                        &
               (STATUS, 'Error defining '//A6(i)//' variable attributes in output file. '// &
               '(INTERFACE define_ovar)')
       END DO

       ! End define mode:
       STATUS = NF90_enddef(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ! Write GrADS coordinate variables
       STATUS = NF90_PUT_VAR(FILE_ID, xvID, REAL(lon_all(:, 1), sp))
       IF(STATUS /= NF90_NOERR) CALL handle_err                                         &
            (STATUS, 'Error writing GrADS x coordinate variable to ' &
            //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       STATUS = NF90_PUT_VAR(FILE_ID, yvID, REAL(lat_all(1, :), sp))
       IF(STATUS /= NF90_NOERR) CALL handle_err                                         &
         (STATUS, 'Error writing GrADS y coordinate variable to ' &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')


       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), real(UNPACK(POPLUC%latitude,landmask, fieldr ),sp))
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), real(UNPACK(POPLUC%longitude, landmask, fieldr),sp) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, latID, real(lat_all,sp))
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, lonID, real(lon_all,sp))
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       CALL1 = .FALSE.
       ! print*, 'OCreated68'

    ENDIF ! CALL1

    ! TIME  ( t )
    ! print*, 'OWrite68 ', file_id
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), real(UNPACK(POPLUC%primf,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), real(UNPACK(POPLUC%secdf,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), real(UNPACK(POPLUC%grass,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), real(UNPACK(POPLUC%ptos,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), real(UNPACK(POPLUC%ptog,landmask, fieldr),sp), &
         start=(/ 1,1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), real(UNPACK(POPLUC%stog, landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), real(UNPACK(POPLUC%gtop,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), real(UNPACK(POPLUC%gtos,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), real(UNPACK(POPLUC%frac_primf,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), real(UNPACK(POPLUC%frac_forest,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 11), real(UNPACK(POPLUC%pharv, landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 12), real(UNPACK(POPLUC%smharv,landmask, fieldr),sp), &
         start=(/ 1, 1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 13), real(UNPACK(POPLUC%syharv,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 14), real(UNPACK(POPLUC%crop,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 15), real(UNPACK(POPLUC%past,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 16), real(UNPACK(POPLUC%ptoc,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 17), real(UNPACK(POPLUC%ptoq,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 18), real(UNPACK(POPLUC%stoc,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 19), real(UNPACK(POPLUC%stoq,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 20), real(UNPACK(POPLUC%qtos,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 21), real(UNPACK(POPLUC%ctos,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 22), real(UNPACK(POPLUC%AgProd,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 23), real(UNPACK(POPLUC%AgProdLoss,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 24), real(UNPACK(POPLUC%FAg,landmask, fieldr),sp), &
         start=(/ 1, 1,CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), UNPACK(POPLUC%n_event, landmask, fieldi), &
         start=(/ 1,1, CNT /), count=(/ nx,ny, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, nLU, t )
    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FHarvest(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FClearance(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FNEP(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%Clitt(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%Csoil(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%CBiomass(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nLU
       tmparr1(:,:,k) = UNPACK(POPLUC%FTransferNet(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), real(tmparr1,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( nx,ny, nTrans, t )
    DO k=1,nTrans
       tmparr2(:,:,k) = UNPACK(POPLUC%FTransferGross(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), real(tmparr2,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nTrans,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( nx,ny, nprod, t )
    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%HarvProd(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), real(tmparr3,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%ClearProd(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(2), real(tmparr3,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%HarvProdLoss(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(3), real(tmparr3,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    DO k=1,nprod
       tmparr3(:,:,k) = UNPACK(POPLUC%ClearProdLoss(:,k),landmask, fieldr)
    ENDDO
    STATUS = NF90_PUT_VAR(FILE_ID, VID6(4), real(tmparr3,sp), &
         start=(/ 1,1,1,CNT /), count=(/ nx,ny,nprod,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    IF ( FINAL ) THEN
       ! Close NetCDF file:
       ! print*, 'OClose68 ', file_id
       STATUS = NF90_close(FILE_ID)
       file_id = -1
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " POPLUC Output written to ",fname
    ENDIF

  END SUBROUTINE WRITE_LUC_OUTPUT_GRID_NC

  !********************************************************************************************************************

END MODULE POPLUC_MODULE
