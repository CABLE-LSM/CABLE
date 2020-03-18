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
! Purpose: subroutines for calculating carbon, nitrogen, phosphorus cycle
!          including plant growth
!
! Called from: biogeochem (mostly) or casa_xnp
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed by Yingping Wang (Wang et al., BG, 2011)
!          Current version uses fixed phenology.
!
! Sep 2015: option of climate-driven phenology (V. Haverd)
!           search for cable_user%PHENOLOGY_SWITCH (Ticket #110)
! May 2016: option of acclimation of auttrophic respiration (V. Haverd)
!            search for cable_user%CALL_climate (Ticket#110)
!         : fixes to prevent carbon and nitrogen pools from going negative
!           search for Ticket#108 (V.Haverd)
!         : alternative functional form of vcmax, called when cable_user%vcmax=='Walker2014'
!           (V.Haverd)
!         : alternative allocation switch integer: LALLOC=3. (V.Haverd)
!           leaf:wood allocation set to maintain LA:SA ratio
!           below target value (requires casaflux%sapwood_area
!           inherited from POP demography module. (Ticket#61)
! Feb 2019 : DAMM Enzyme kinetics implemented as a switchable option for soil and litter
! carbon tunover responses to temperature and moisture
! Feb 2020 : additional turnover of plant and litter pools due to fire (inherited from BLAZE)
! ==============================================================================
!
! This module contains the following subroutines:
!   casa_xnp
!   casa_allocation
!   casa_rplant
!   casa_xrateplant,        casa_xratesoil
!   casa_coeffplant,        casa_coeffsoil
!   casa_delplant,          casa_delsoil
!   avgsoil
!   casa_xkN
!   casa_nuptake,           casa_puptake
!   casa_Nrequire,          casa_Prequire
!   casa_cnpcycle
!   casa_poolzero
!   casa_cnpbal
!   casa_ndummy
!   phenology

MODULE casa_cnp_module
  
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_IO_vars_module, ONLY: wlogn
  USE cable_common_module,  only: cable_user ! Custom soil respiration: Ticket #42

  implicit none

  real(r_2), parameter :: zero = 0.0_r_2
  real(r_2), parameter :: one  = 1.0_r_2
  
CONTAINS


SUBROUTINE casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xnplimit
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xNPuptake
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  INTEGER :: np
  REAL(r_2), DIMENSION(mp)        :: xnlimit,xplimit
  REAL(r_2), DIMENSION(mp)        :: xncleaf,xpcleaf
  REAL(r_2), DIMENSION(mp)        :: xnCnpp,xpCnpp
  REAL(r_2), DIMENSION(mp,mplant) :: Nreqmax, Nreqmin, NtransPtoP
  REAL(r_2), DIMENSION(mp)        :: totNreqmax,totNreqmin
  REAL(r_2), DIMENSION(mp)        :: xNuptake,xPuptake
  REAL(r_2), DIMENSION(mp,mplant) :: Preqmax, Preqmin, PtransPtoP
  REAL(r_2), DIMENSION(mp)        :: totPreqmax,totPreqmin

  xnlimit  = 1.0_r_2
  xplimit  = 1.0_r_2
  xnplimit = 1.0_r_2
  casaflux%fracClabile(:) = 0.0_r_2

  SELECT CASE(icycle)
  CASE(2)
    WHERE(casamet%iveg2/=icewater)
      xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10_r_2)
      !xnlimit(:) = xncleaf(:)/(xncleaf(:)+casabiome%KminN(veg%iveg(:)))
      xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01_r_2)
      xplimit(:) = 1.0_r_2
      xnplimit(:) =min(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
    ENDWHERE
  CASE(3)
    WHERE(casamet%iveg2/=icewater)
      xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10_r_2)
      xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01_r_2)
      xpcleaf(:) = casapool%pplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10_r_2)
      !xplimit(:) = xpcleaf(:)/(xpcleaf(:)+casabiome%Kuplabp(veg%iveg(:)))
      xplimit(:) = xpcleaf(:)/(xpcleaf(:)+0.0006_r_2)
      !xnplimit(:) = min(1.0_r_2,casabiome%Kuptake(veg%iveg(:))*min(xnlimit(:),xplimit(:)))
      xnplimit(:) =min(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
    ENDWHERE
  END SELECT

  ! now check if soil nutrient supply can meet the plant uptake,
  ! otherwise reduce NPP
  xNuptake = 1.0_r_2
  xPuptake = 1.0_r_2

  IF(icycle >1) THEN
    Nreqmin(:,:)    = 0.0_r_2
    Nreqmax(:,:)    = 0.0_r_2
    NtransPtoP(:,:) = 0.0_r_2
    totNreqmax = 0.0_r_2
    totNreqmin = 0.0_r_2
    xNuptake   = 1.0_r_2

    xnCnpp = max(0.0_r_2,casaflux%Cnpp)
    call casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
    DO np=1,mp
      IF(casamet%iveg2(np)/=icewater) THEN
        totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
        totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)
        xNuptake(np)   = MAX(0.0_r_2,MIN(1.0_r_2,casapool%Nsoilmin(np) &
                                         /(totNreqmin(np)*deltpool+1.0e-10_r_2)))
      ENDIF
    ENDDO
  ENDIF
  IF(icycle >2) THEN
    Preqmin(:,:)       = 0.0_r_2
    Preqmax(:,:)       = 0.0_r_2
    PtransPtoP(:,:)    = 0.0_r_2
    totPreqmax = 0.0_r_2
    totPreqmin = 0.0_r_2
    xPuptake   = 1.0_r_2
    xpCnpp = max(0.0_r_2,casaflux%Cnpp)
    call casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                       casabiome,casapool,casaflux,casamet)
    DO np=1,mp
      IF(casamet%iveg2(np)/=icewater) THEN
        totPreqmax(np) = Preqmax(np,leaf)+Preqmax(np,wood)+Preqmax(np,froot)
        totPreqmin(np) = Preqmin(np,leaf)+Preqmin(np,wood)+Preqmin(np,froot)
        xPuptake(np)   = MAX(0.0_r_2,MIN(1.0_r_2,casapool%psoillab(np) &
                                         /(totPreqmin(np)*deltpool+1.0e-10_r_2)))
      ENDIF
    ENDDO
  ENDIF

  xnplimit(:)  = 1.0_r_2
  xNPuptake(:)     = min(xnuptake(:), xpuptake(:))
  do np =1, mp
     if(casamet%iveg2(np)/=icewater.and.casaflux%cnpp(np) > 0.0_r_2.and.xNPuptake(np) < 1.0_r_2) then

        if (casaflux%fHarvest(np) .lt. 0.1_r_2) then   ! N limitation on growth not applied to crop/pasture
           casaflux%fracClabile(np) =min(1.0_r_2,max(0.0_r_2,(1.0_r_2- xNPuptake(np)))) * &
                max(0.0_r_2,casaflux%cnpp(np))/(casaflux%cgpp(np) +1.0e-10_r_2)
           casaflux%cnpp(np)    = casaflux%cnpp(np) - casaflux%fracClabile(np) * casaflux%cgpp(np)
        endif
     endif

   enddo

END SUBROUTINE casa_xnp


SUBROUTINE casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
! compute fraction of net photosynthate allocated to leaf, wood and froot
!
! inputs
!   moistavg(mp)           as an argument (volume fraction)
!   tsoilavg(mp)           as an argument (K)
!   btran(mp)              as an argument (dimensionless)
! outputs:
!   fracCalloc(mp,mplant1)
!
! modified Piere's alocation scheme
! input: leaf stage
!        leaf area

  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (phen_variable),       INTENT(INOUT) :: phen
  INTEGER , INTENT(IN) :: LALLOC
  ! local variables
  REAL(r_2), DIMENSION(mp,mplant) :: fracCallocx
  REAL(r_2), DIMENSION(mp)        :: xLalloc,xwsalloc,xTalloc
  REAL(r_2), DIMENSION(mp)        :: xWorNalloc,xNalloc,xWalloc
  REAL(r_2), DIMENSION(mp)        :: totfracCalloc
  REAL(r_2), DIMENSION(mp)        :: newLAI
  
  ! initlization
  casaflux%fracCalloc  = 0.0_r_2
  !casaflux%fracClabile = 0.0_r_2
  fracCallocx = 0.0_r_2
  newLAI = 0.0_r_2
  SELECT CASE (LALLOC)

  CASE(2)   !
    ! calculate the allocation coefficients
    call casa_wolf(veg,casabiome,casaflux,casapool,casamet)

  CASE(1)   ! dynamic allocation
    WHERE(casamet%iveg2/=icewater)
      xLalloc(:) = min(1.0_r_2,max(0.0_r_2,exp(-0.5_r_2*casamet%glai(:))))   ! L limiting
      ! Pseudo-nutrient limitation calculation
      WHERE(casamet%tsoilavg > 0.0_r_2)
        xwsalloc(:) = min( max(casamet%moistavg(:)-soil%swilt(:),0.0_r_2) &
                         /(soil%sfc(:)-soil%swilt(:)), 1.0_r_2)
      ELSE WHERE
        xwsalloc(:) = 0.01_r_2
      END WHERE
      xTalloc(:)    = min(1.0_r_2,max(0.0_r_2,Q10alloc** &
                      ((casamet%tsoilavg(:)-TkzeroC-30.0_r_2)/10.0_r_2) )) !T limiting
      xNalloc(:)    = min(1.0_r_2,max(0.0_r_2,xwsalloc(:)*xTalloc(:)))     !N limiting
      xWalloc(:)    = min(1.0_r_2,max(0.0_r_2,casamet%btran(:)))           !W limiting
      xWorNalloc(:) = min(xWalloc(:),xNalloc(:))
      WHERE(casamet%lnonwood==0)
        casaflux%fracCalloc(:,FROOT) = R0 * 3.0_r_2 * xLalloc(:) &
                                     / (xLalloc(:)+ 2.0_r_2*xWorNalloc(:))
        casaflux%fracCalloc(:,WOOD)  = S0 * 3.0_r_2 * xWorNalloc(:) &
                                     / (2.0_r_2*xLalloc(:)+ xWorNalloc(:))
        casaflux%fracCalloc(:,LEAF)  = 1.0_r_2 - casaflux%fracCalloc(:,FROOT) &
                                     - casaflux%fracCalloc(:,WOOD)
      ELSE WHERE
        casaflux%fracCalloc(:,FROOT) = R0 * 3.0_r_2 * xLalloc(:) &
                                     / (xLalloc(:)+2.0_r_2*xWorNalloc(:))
        casaflux%fracCalloc(:,WOOD)  = 0.0_r_2
        casaflux%fracCalloc(:,LEAF)  = 1.0_r_2 - casaflux%fracCalloc(:,FROOT)
      END WHERE
    END WHERE
  CASE (0)   ! fixed allocation
    casaflux%fracCalloc(:,:) = casabiome%fracnpptop(veg%iveg(:),:)

  CASE (3) ! leaf:wood allocation set to maintain LA:SA ratio
     ! below target value of casabiome%la_to_sa, where phen%phase = 1 or 2
     !(requires casaflux%sapwood_area, which is inherited from the
     ! POP tree demography module. (Ticket #61)
    WHERE(casamet%lnonwood==0)
        casaflux%fracCalloc(:,FROOT) =  casabiome%fracnpptop(veg%iveg(:),FROOT)
        casaflux%fracCalloc(:,WOOD) = 0.01_r_2
        casaflux%fracCalloc(:,LEAF) = 1.0_r_2 - casaflux%fracCalloc(:,FROOT) - &
             casaflux%fracCalloc(:,WOOD)
        newLAI =casamet%glai + (casaflux%fracCalloc(:,LEAF) *casaflux%cnpp- &
             casaflux%kplant_tot(:,leaf) *casapool%cplant(:,LEAF) )*casabiome%sla(veg%iveg(:))
        where (casaflux%sapwood_area.gt.1.e-6_r_2 .and. newLAI.gt.(casabiome%la_to_sa(veg%iveg(:)) &
             *casaflux%sapwood_area) &
             .and. casaflux%cnpp.gt.0.0_r_2)

           casaflux%fracCalloc(:,LEAF) = ((casabiome%la_to_sa(veg%iveg(:))*casaflux%sapwood_area &
                - casamet%glai)/ &
                casabiome%sla(veg%iveg(:)) &
             + casaflux%kplant_tot(:,leaf) *casapool%cplant(:,LEAF)  )/casaflux%cnpp

           casaflux%fracCalloc(:,LEAF) = max(0.0_r_2,  casaflux%fracCalloc(:,LEAF) )
           casaflux%fracCalloc(:,LEAF) = min(1.0_r_2 - casaflux%fracCalloc(:,FROOT) - &
                casaflux%fracCalloc(:,WOOD) ,&
             casaflux%fracCalloc(:,LEAF) )

           casaflux%fracCalloc(:,WOOD) = 1.0_r_2 -  casaflux%fracCalloc(:,FROOT) - &
                casaflux%fracCalloc(:,LEAF)
        end where


     ELSEWHERE

        casaflux%fracCalloc(:,FROOT) =  casabiome%fracnpptop(veg%iveg(:),FROOT)
        casaflux%fracCalloc(:,WOOD) = 0.0_r_2
        casaflux%fracCalloc(:,LEAF) =  casabiome%fracnpptop(veg%iveg(:),LEAF)

     ENDWHERE

  END SELECT

  ! during leaf growth phase 0 or 3, no carbon is allocated to leaf,
  ! during maximal leaf growth phase, all C is allocated to leaf
  ! during steady growth period, C allocation is estimated in such
  ! a way that approach the allometric relationship
  ! the relationships are:(all pools in g C/m2)
  ! for forests
  !   fineroot/totalC C=0.3192-0.0485*(totalC)^0.1755, see mokany et al. (2003)
  !   fineroot = ratiofrootleaf*cleaf
  ! for grassland
  !   root=ratiofinerootleaf*cleaf

! vh edit to avoid overwriting CASE(3) for woody veg
!! vh_js !!
  IF (LALLOC.ne.(3)) THEN

     WHERE(casamet%iveg2/=icewater)
        WHERE(phen%phase==0)
           casaflux%fracCalloc(:,leaf)  = 0.0_r_2
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           casaflux%fracCalloc(:,wood)  = 1.0_r_2 -casaflux%fracCalloc(:,froot)
        END WHERE

        WHERE(phen%phase==1)
           casaflux%fracCalloc(:,leaf)  = 0.95_r_2
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,froot) = 0.5_r_2*(1.0_r_2-casaflux%fracCalloc(:,leaf))
              casaflux%fracCalloc(:,wood)  = 0.5_r_2*(1.0_r_2-casaflux%fracCalloc(:,leaf))
           ELSEWHERE !grassland
              casaflux%fracCalloc(:,froot) = 1.0_r_2-casaflux%fracCalloc(:,leaf)
           ENDWHERE
        END WHERE

        WHERE(phen%phase==3)
           !      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
           casaflux%fracCalloc(:,froot) = 1.0_r_2-casaflux%fracCalloc(:,wood)
           casaflux%fracCalloc(:,leaf)    = 0.0_r_2
        ENDWHERE


        ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
        ! Q.Zhang 17/03/2011
        WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
           casaflux%fracCalloc(:,leaf)  = 0.0_r_2
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           casaflux%fracCalloc(:,wood)  = 1.0_r_2 -casaflux%fracCalloc(:,froot)
        ENDWHERE

        ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
        WHERE(casaflux%Cnpp<0.0_r_2)
           casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)  / sum(casaflux%Crmplant,2)
           casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)  / sum(casaflux%Crmplant,2)
           casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot) / sum(casaflux%Crmplant,2)
        ENDWHERE

        !! vh_js !!
        !! as long as biomass is positive, adjust allocation to be
        !! proportional to stock when NPP -ve (Ticket#108)
        WHERE(casaflux%Cnpp<0.0_r_2 .and. sum(casapool%Cplant,2)>0.0_r_2  )
           casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)  / sum(casapool%Cplant,2)
           casaflux%fracCalloc(:,wood)  = casapool%Cplant(:,wood)  / sum(casapool%Cplant,2)
           casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot) / sum(casapool%Cplant,2)
        ENDWHERE
     ENDWHERE

  ELSE
     WHERE(casamet%iveg2/=icewater)
        WHERE(phen%phase==0)
           casaflux%fracCalloc(:,leaf)  = 0.0_r_2
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           WHERE (casamet%lnonwood==0)
              casaflux%fracCalloc(:,wood)  = 1.0_r_2 -casaflux%fracCalloc(:,froot)
           ELSEWHERE
              casaflux%fracCalloc(:,wood) = 0.0_r_2
           ENDWHERE
        END WHERE

        WHERE(phen%phase==1.and.casamet%lnonwood==1)

           casaflux%fracCalloc(:,leaf)  = 0.80
           casaflux%fracCalloc(:,froot) = 1.0_r_2-casaflux%fracCalloc(:,leaf)
           casaflux%fracCalloc(:,wood) = 0.0_r_2
        ENDWHERE

        WHERE(phen%phase==3)
           !      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
           casaflux%fracCalloc(:,froot) = 1.0_r_2-casaflux%fracCalloc(:,wood)
           casaflux%fracCalloc(:,leaf)  = 0.0_r_2
        ENDWHERE

!! vh !!
!! This fix can lead to over-allocation to roots, in turn bumping up N-uptake
!! , leading to decline in mineral nitrogen availability and spikes in fracCalloc,
!! causing spikes in tree mortality and lack of model convergence in productive
!! regions where LAI is hitting LAImax: for woody pfts, ensure that LAImax is
!! parameter is very high (eg 10)
        ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
        ! Q.Zhang 17/03/2011
        WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
           casaflux%fracCalloc(:,leaf)  = 0.0_r_2
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           WHERE (casamet%lnonwood==0)
              casaflux%fracCalloc(:,wood)  = 1.0_r_2 -casaflux%fracCalloc(:,froot)
           ELSEWHERE
              casaflux%fracCalloc(:,wood) = 0.0_r_2
           ENDWHERE
        ENDWHERE

        WHERE(casamet%glai(:)<casabiome%glaimin(veg%iveg(:)))
           casaflux%fracCalloc(:,leaf)  = 0.8
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,froot) = 0.5*(1.0_r_2-casaflux%fracCalloc(:,leaf))
              casaflux%fracCalloc(:,wood)  = 0.5*(1.0_r_2-casaflux%fracCalloc(:,leaf))
           ELSEWHERE !grassland
              casaflux%fracCalloc(:,froot) = 1.0_r_2-casaflux%fracCalloc(:,leaf)
              casaflux%fracCalloc(:,wood) = 0.0_r_2
           ENDWHERE
        ENDWHERE
        ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
        WHERE(casaflux%Cnpp<0.0_r_2)
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)  / sum(casaflux%Crmplant,2)
              casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)  / sum(casaflux%Crmplant,2)
              casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot) / sum(casaflux%Crmplant,2)
           ELSEWHERE
              casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)  / sum(casaflux%Crmplant,2)
              casaflux%fracCalloc(:,wood)  = 0.0_r_2
              casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot) / sum(casaflux%Crmplant,2)
           ENDWHERE
        ENDWHERE

        !! vh_js !!  Ticket#108
        WHERE(casaflux%Cnpp<0.0_r_2 .and. sum(casapool%Cplant,2)>0.0_r_2  )
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)  / sum(casapool%Cplant,2)
              casaflux%fracCalloc(:,wood)  = casapool%Cplant(:,wood)  / sum(casapool%Cplant,2)
              casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot) / sum(casapool%Cplant,2)
           ELSEWHERE
              casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)  / sum(casapool%Cplant,2)
              casaflux%fracCalloc(:,wood)  = 0.0_r_2
              casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot) / sum(casapool%Cplant,2)
           ENDWHERE
        ENDWHERE


     ENDWHERE
  !   write(*,*) 'alloc2',  casaflux%fracCalloc(1,2), casaflux%Cnpp(1), casapool%Cplant(1,:), &
  !        casamet%lnonwood(1)
!if (ANY(casapool%Cplant(1,:).NE.casapool%Cplant(1,:))) then
!write(*,*) 'cplant', casapool%Cplant(1,:)
!stop
!endif
  ENDIF ! LALLOC=3

!write(*,*) 'alloc 4', casaflux%fracCalloc(4,:),  casapool%Cplant(4,:)
!write(*,*) 'alloc 5', casaflux%fracCalloc(5,:), casapool%Cplant(5,:)
  ! normalization the allocation fraction to ensure they sum up to 1
  totfracCalloc(:) = sum(casaflux%fracCalloc(:,:),2)
  casaflux%fracCalloc(:,leaf)  = casaflux%fracCalloc(:,leaf)  / totfracCalloc(:)
  casaflux%fracCalloc(:,wood)  = casaflux%fracCalloc(:,wood)  / totfracCalloc(:)
  casaflux%fracCalloc(:,froot) = casaflux%fracCalloc(:,froot) / totfracCalloc(:)

END SUBROUTINE casa_allocation


SUBROUTINE casa_wolf(veg,casabiome,casaflux,casapool,casamet)
   ! carbon allocation based on
   ! Wolf, Field and Berry, 2011. Ecological Applications, p1546-1556
   ! Wolf et al. 2011. Global Biogeochemical Cycles, 25, GB3015, doi:10.1019/2010GB003917
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(IN) :: casabiome
  TYPE (casa_met),            INTENT(IN) :: casamet
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux

   real(r_2), parameter :: wolf_alpha1=6.22_r_2
   real(r_2), parameter :: wolf_beta=-1.33_r_2
   real(r_2), parameter :: wolf_c1=-wolf_alpha1/(1._r_2+wolf_beta)
   real(r_2), parameter :: wolf_c2=1.0_r_2/(1.0_r_2+wolf_beta)
   !
   ! local variables
   integer   npt
   real(r_2), dimension(mp) :: totbmdm,ntree,nppdm
   real(r_2), dimension(mp) :: gleaf,gwood,gcroot,gfroot,gtot
   !
   ! input
   !  totleaf, totwood, totcroot, totfroot :    g C m-2
   !  totnpp:                                   g C m-2 d-1
   ! output
   !  fracleaf,fracwood, fraccroot, fracfroot:  fractions
   !

    do npt=1,mp
       IF(casamet%iveg2(npt)==3.and.casaflux%cnpp(npt)>0.0001_r_2) THEN  !forest types
          totbmdm(npt) = sum(casapool%cplant(npt,:)) *10.0_r_2 / fracCbiomass      !10.0 for convert gc/m2 to kg/ha
          totbmdm(npt) = max(30000.0_r_2, totbmdm(npt))
          ! calculate tree stocking density
           ntree(npt) = 10._r_2**(wolf_c1+wolf_c2*log10(totbmdm(npt)))   ! tree ha-1, based on eqn (4) of Wolf et al. 2011, GBC
           ntree(npt) = min(200000.0_r_2,ntree(npt))
           ! changed by ypw 23/april/2012 to avoid negative npp
           nppdm(npt)  = (abs(casaflux%cnpp(npt)) *365.0_r_2*0.001_r_2/fracCbiomass)/(0.0001_r_2*ntree(npt))  ! in kg dm tree-1 yr-1

           gleaf(npt)  = 0.156_r_2*(nppdm(npt)**1.106_r_2)     ! Figure 2a of Wolf, Field and Berry (2011)
           gwood(npt)  = 0.232_r_2*(nppdm(npt)**1.165_r_2)     ! Figure 2b of Wolf, Field and Berry (2011)
           gcroot(npt) = 0.0348_r_2*(nppdm(npt)**1.310_r_2)    ! Figure 2d of Wolf, Field and Berry (2011)
           gfroot(npt) = 0.247_r_2*(nppdm(npt)**0.987_r_2)     ! Figure 2c of Wolf, Field and Berry (2011)
           gtot(npt)   = gleaf(npt) + gwood(npt) + gcroot(npt) + gfroot(npt)

           casaflux%fracCalloc(npt,leaf)  = gleaf(npt)/gtot(npt)
           casaflux%fracCalloc(npt,wood)  = gwood(npt)/gtot(npt)
           casaflux%fracCalloc(npt,froot) = (gcroot(npt)+gfroot(npt))/gtot(npt)
!        write(87,*) 'allocation = ',npt,casamet%iveg2(npt), totbmdm(npt),ntree(npt),nppdm(npt),casaflux%fracCalloc(npt,:)
        ELSE                ! other types
           casaflux%fracCalloc(npt,:) = casabiome%fracnpptop(veg%iveg(npt),:)
        ENDIF
    enddo

END SUBROUTINE casa_wolf


SUBROUTINE casa_rplant(veg, casabiome, casapool, casaflux, casamet, climate)

  ! maintenance respiration of woody tisse and fineroots
  ! see Sitch et al. (2003), GCB, reqn (23)

  implicit none
  
  type(veg_parameter_type), intent(inout) :: veg       ! vegetation parameters
  type(casa_biome),         intent(inout) :: casabiome
  type(casa_pool),          intent(inout) :: casapool
  type(casa_flux),          intent(inout) :: casaflux
  type(casa_met),           intent(inout) :: casamet
  type(climate_type),       intent(in)    :: climate

  integer :: npt, ivt
  real(r_2), dimension(mp)        :: Ygrow        ! growth efficiency Q.Zhang 22/02/2011
  real(r_2), dimension(mp,mplant) :: ratioPNplant ! Q.Zhang 22/02/2011
  real(r_2), dimension(mp)        :: delcrmwood,delcrmfroot    ! reduction in wood and root respiration when NPP <0.0
  real(r_2), dimension(mp)        :: resp_coeff_root, resp_coeff_sapwood
  real(r_2), dimension(mp)        :: nleaf, pleaf, vcmaxmax
  real(r_2) :: c1, c2, c3, c5 ! coefficients for acclimation of maintenance respiration

  resp_coeff_root = 1._r_2
  resp_coeff_sapwood = 1._r_2
  ratioPNplant = 0.0_r_2
  Ygrow        = 0.0_r_2

  WHERE(casapool%Nplant>0.0_r_2)
    ratioPNplant = casapool%Pplant/(casapool%Nplant+ 1.0e-10_r_2)
  ENDWHERE

  Ygrow(:) = 0.65_r_2 + 0.2_r_2*ratioPNplant(:,leaf)/(ratioPNplant(:,leaf)+1.0_r_2/15.0_r_2)

  casaflux%crmplant(:,wood)  = 0.0_r_2
  casaflux%crmplant(:,froot) = 0.0_r_2
  delcrmwood  = 0.0_r_2
  delcrmfroot = 0.0_r_2
  casaflux%crgplant = 0.0_r_2
  casaflux%clabloss = 0.0_r_2

  if (cable_user%CALL_climate) then
     ! coefficients required to implement T-acclimation of autotrophic respiration (Ticket # 110)
     ! adapted from Atkin et al., New Phyt., 2015)
     do npt=1, mp
        ivt=veg%iveg(npt)
        ! max leaf N in g N m-2 leaf
        nleaf(npt) = casabiome%ratioNCplantmax(ivt,leaf)/casabiome%sla(ivt)
        ! max leaf P in g P m-2 leaf
        pleaf(npt) = casabiome%ratioPcplantmax(ivt,leaf)/casabiome%sla(ivt)
        if (ivt.eq.7) then
           ! special for C4 grass: set here to value from  parameter file
           vcmaxmax(npt) = veg%vcmax(npt)
        else
           vcmaxmax(npt) = vcmax_np(real(nleaf(npt)), real(pleaf(npt)))*casabiome%vcmax_scalar(ivt)
        endif

        if (veg%iveg(npt).eq.2 .or. veg%iveg(npt).eq.4) then
           ! broadleaf forestNESP2pt9_BLAZE
           c1 = 1.2818_r_2 * 1.e-6_r_2
           ! c1 = 1.3805e-6_r_2 ! null model
        elseif(veg%iveg(npt).eq.1 .or. veg%iveg(npt).eq.3) then
           ! needleleaf forest
           c1 = 1.2877_r_2 * 1.e-6_r_2
           ! c1 = 1.3247e-6_r_2 ! null model
        elseif (veg%iveg(npt).eq.6 .or. veg%iveg(npt).eq.8 .or. veg%iveg(npt).eq.9) then
           ! C3 grass, tundra, crop
           c1 = 1.6737_r_2 * 1.e-6_r_2
           ! c1 = 1.8904e-6_r_2 ! null model
        else
           ! shrubs and other (C4 grass and crop)
           c1 = 1.5758_r_2 * 1.e-6_r_2
           !c1 = 1.7265e-6_r_2 ! null model
        endif
        c2 = 0.0116_r_2
        c3 = -0.0334_r_2*1.e-6_r_2
        c5 = 1./vcmaxmax(npt)/0.0116_r_2*0.60_r_2
        resp_coeff_root(npt) = casapool%nplant(npt,froot) * c5 * &
             (c1 + c2 *vcmaxmax(npt)*climate%frec(npt)  + c3 * climate%qtemp_max_last_year(npt) )
        
        resp_coeff_sapwood(npt) = casapool%nplant(npt,wood) *casaflux%frac_sapwood(npt)* c5 * &
             (c1 + c2 *vcmaxmax(npt)*climate%frec(npt)  + c3 * climate%qtemp_max_last_year(npt) )
        ! null model (use this with null-model c1 coefft to turn off T-acclimation)
        ! resp_coeff_root(npt) = casapool%nplant(npt,froot) *c1
        ! resp_coeff_sapwood(npt) = casapool%nplant(npt,wood) *casaflux%frac_sapwood(npt)*c1
     enddo ! npt=1,mp
  endif  ! cable_user%CALL_climate - end coefficients for acclimation of autotrophic respiration ticket #110 - 

  if (cable_user%CALL_climate) then
     !  acclimation of autotrophic respiration Ticket #110
     where(casamet%iveg2/=icewater)
        
        where(casamet%tairk >250.0_r_2)
           where(casapool%cplant(:,wood)>1.0e-6_r_2)
              casaflux%crmplant(:,wood)  =   resp_coeff_sapwood * &
                   casabiome%rmplant(veg%iveg(:),wood) &
                   * exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2           &
                   / (casamet%tairk(:)+46.02_r_2-tkzeroc)))

           endwhere
           !vh! prevent floating underflow with this mask
           where (casapool%Clabile(:).gt.1.e-8_r_2)
              casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                   * max(0.0_r_2,casapool%Clabile(:))      &
                   * exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2         &
                   / (casamet%tairk(:)+46.02_r_2-tkzeroc)))
           endwhere
        endwhere
        
        where(casamet%tsoilavg >250.0_r_2 .and. casapool%cplant(:,froot)>1.0e-6_r_2)
           casaflux%crmplant(:,froot) =   resp_coeff_root * &
                casabiome%rmplant(veg%iveg(:),froot) &
                * exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2  &
                / (casamet%tsoilavg(:)+46.02_r_2-tkzeroc)))
        endwhere

        where((casaflux%Cgpp-sum(Casaflux%crmplant,2))>0.0_r_2)
           !casaflux%crgplant(:)  = 0.25* max(0.0_r_2,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
           ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
           casaflux%crgplant(:)  = (1.0_r_2-Ygrow(:))* &
                max(0.0_r_2,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
        elsewhere
           casaflux%crgplant(:) = 0.0_r_2
        endwhere
        
     endwhere ! /= icewater

     casaflux%Cnpp(:) = casaflux%Cgpp(:) - sum(casaflux%Crmplant(:,:),2) - casaflux%Crgplant(:)

  else ! cable_user%CALL_climate

     where(casamet%iveg2/=icewater)
        
        where(casamet%tairk >250.0_r_2)
           where(casapool%cplant(:,wood)>1.0e-6_r_2)
              casaflux%crmplant(:,wood)  =   casaflux%frac_sapwood(:) * &
                   casabiome%rmplant(veg%iveg(:),wood) &
                   * casapool%nplant(:,wood)             &
                   * exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2           &
                   / (casamet%tairk(:)+46.02_r_2-tkzeroc)))
           endwhere
           casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                * max(0.0_r_2,casapool%Clabile(:))      &
                * exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2         &
                / (casamet%tairk(:)+46.02_r_2-tkzeroc)))
        endwhere
        
        where(casamet%tsoilavg >250.0_r_2 .and. casapool%cplant(:,froot)>1.0e-6_r_2)
           casaflux%crmplant(:,froot) =  casabiome%rmplant(veg%iveg(:),froot) &
                * casapool%nplant(:,froot)             &
                * exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2            &
                / (casamet%tsoilavg(:)+46.02_r_2-tkzeroc)))
        endwhere
        !    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + casaflux%clabloss(:)

        where((casaflux%Cgpp-sum(casaflux%crmplant,2))>0.0_r_2)
           !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
           ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
           casaflux%crgplant(:)  = (1.0_r_2-Ygrow(:))* &
                max(0.0_r_2,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
        elsewhere
           casaflux%crgplant(:) = 0.0_r_2
        endwhere
        !casaflux%Cnpp(:) = MAX(0.0_r_2,(casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2) &
        !                 - casaflux%crgplant(:)))
        ! changes made by yp wang 5 april 2013
        casaflux%Cnpp(:) = casaflux%Cgpp(:) - sum(casaflux%Crmplant(:,:),2) - casaflux%Crgplant(:)

     endwhere ! /= icewater

  endif ! cable_user%CALL_climate

END SUBROUTINE casa_rplant


SUBROUTINE casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                           casamet,phen)
! use xleafcold and xleafdry to account for
! cold and drought stress on death rate of leaf
! inputs:
!     ivt(mp)  :       biome type
!     phase(mp):       leaf growth stage
!     tairk(mp)    :   air temperature in K
! outputs
!     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
!     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
!     xkleaf(mp):      set to 0.0 during maximal leaf growth phase

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleafcold
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleafdry
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleaf
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mp)              :: xcoldleaf
  INTEGER :: npt

  xkleafcold(:) = 0.0_r_2
  xkleafdry(:)  = 0.0_r_2
  xkleaf(:)     = 1.0_r_2

  ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
  !    following the formulation of Arora (2005) on the
  !    effect of cold or drought stress on leaf litter fall
  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
    IF(casamet%tairk(npt)>=phen%TKshed(veg%iveg(npt))) THEN
      xcoldleaf(npt) = 1.0_r_2
    ELSE
      IF(casamet%tairk(npt)<=(phen%TKshed(veg%iveg(npt))-5.0)) THEN
        xcoldleaf(npt)=0.0_r_2
      ELSE
        xcoldleaf(npt) = (casamet%tairk(npt)-phen%TKshed(veg%iveg(npt))-5.0)/5.0
      ENDIF
    ENDIF
    xcoldleaf(npt) = min(1.0_r_2,max(0.0_r_2,xcoldleaf(npt)))
    xkleafcold(npt) = casabiome%xkleafcoldmax(veg%iveg(npt)) &
                    * (1.0_r_2-xcoldleaf(npt)) &
                    ** casabiome%xkleafcoldexp(veg%iveg(npt))
    xkleafdry(npt)  = casabiome%xkleafdrymax(veg%iveg(npt)) &
                    * (1.0_r_2-casamet%btran(npt))&
                    ** casabiome%xkleafdryexp(veg%iveg(npt))
    IF (phen%phase(npt)==1) xkleaf(npt)= 0.0_r_2
    ! vh: account for high rate of leaf loss during senescence
    ! vh_js
    if (trim(cable_user%PHENOLOGY_SWITCH)=='climate') then
       ! increases base turnover rate by a factor of 13 (for base turnover time of 1y, this reduces it to 4 weeks)
       IF ((phen%phase(npt)==3.or.phen%phase(npt)==0).and.casamet%lnonwood(npt)==0) &
                xkleaf(npt)= 13.
        IF ((phen%phase(npt)==3.or.phen%phase(npt)==0).and.casamet%lnonwood(npt)==1) &
                xkleaf(npt)= 13.*0.5

    endif
  END IF
  END DO

!  WHERE(casamet%iveg2/=icewater)
!  !    following the formulation of Arora (2005) on the
!  !    effect of cold or drought stress on leaf litter fall
!  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
!    WHERE(casamet%tairk(:)>=phen%TKshed(veg%iveg(:)))
!      xcoldleaf(:) = 1.0_r_2
!    ELSEWHERE
!      WHERE(casamet%tairk(:)<=(phen%TKshed(veg%iveg(:))-5.0))
!        xcoldleaf(:)=0.0_r_2
!      ELSEWHERE
!        xcoldleaf(:) = (casamet%tairk(:)-phen%TKshed(veg%iveg(:))-5.0)/5.0
!      ENDWHERE
!    ENDWHERE
!    xcoldleaf(:) = min(1.0_r_2,max(0.0_r_2,xcoldleaf(:)))
!    xkleafcold(:) = casabiome%xkleafcoldmax(veg%iveg(:)) * (1.0_r_2-xcoldleaf(:)) &
!                 ** casabiome%xkleafcoldexp(veg%iveg(:))
!    xkleafdry(:) = casabiome%xkleafdrymax(veg%iveg(:))*(1.0_r_2-casamet%btran(:))&
!                 ** casabiome%xkleafdryexp(veg%iveg(:))
!    WHERE(phen%phase(:)==1) xkleaf(:)= 0.0_r_2
!  ENDWHERE

END SUBROUTINE casa_xrateplant


SUBROUTINE casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
!  to account for cold and drought stress on death rate of leaf: xleafcold,xleafdry
!  to account for effects of T and W on litter decomposition: xk, xksurf
!  inputs:
!     ivt(mp)  :       biome type
!     tsoilavg(mp):    soil temperature in K
!     moistavg(mp):    volumetric soil moisture
!
!  outputs
!     xk(mp):          modifier of soil litter decomposition rate (dimensionless)
  IMPLICIT NONE
  
  REAL(r_2), DIMENSION(mp),  INTENT(OUT)   :: xklitter, xksoil
  TYPE(veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE(casa_met),            INTENT(INOUT) :: casamet
  TYPE(casa_biome),          INTENT(INOUT) :: casabiome

  ! local variables
  REAL(r_2), parameter :: wfpscoefa=0.55_r_2   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps
  REAL(r_2), parameter :: wfpscoefb=1.70_r_2   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007_r_2 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22_r_2   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481_r_2 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)
  ! ! Kirschbaum function parameters
  ! REAL(r_2), parameter :: xkalpha=-3.764   ! Kirschbaum (1995, SBB)
  ! REAL(r_2), parameter :: xkbeta=0.204
  ! REAL(r_2), parameter :: xktoptc=36.9
  ! Trudinger2016 function parameters (from Trudinger 2016)
  REAL(r_2), parameter ::  wfpswidth1=1.2160310E+00
  REAL(r_2), parameter ::  wfpswidth2=6.4620150E-01
  REAL(r_2), parameter ::  wfpswidth3=4.9625073E-01
  REAL(r_2), parameter :: wfpsscale1=1.6193957E+00
  REAL(r_2), parameter :: wfpswidth0=1.3703876E-02
  REAL(r_2), parameter :: wfpsquad=8.0000000E-01

  ! Trudinger2016 function parameters (corresponds to Haverd 2013)
  ! REAL(r_2), parameter :: wfpswidth1=0.78_r_2
  ! REAL(r_2), parameter :: wfpswidth2=0_r_2
  ! REAL(r_2), parameter :: wfpswidth3=1.5_r_2
  ! REAL(r_2), parameter :: wfpsscale1=10.0_r_2
  ! REAL(r_2), parameter :: wfpswidth0 = 0.0_r_2
  ! REAL(r_2), parameter :: wfpsquad=0._r_2

  ! DAMM temporary variables
  REAL(r_2) :: O2, vol_air_content, Enz, Dliq, Dva

  REAL(r_2), DIMENSION(mp)       :: xkwater,xktemp
  REAL(r_2), DIMENSION(mp)       :: fwps,tsavg
  ! Custom soil respiration - see Ticket #42
  REAL(r_2), DIMENSION(mp)       :: smrf,strf,slopt,wlt,tsoil,sopt
  REAL(r_2) :: f0
  !,tsurfavg  !!, msurfavg
  INTEGER :: npt

 xklitter(:) = 1.0_r_2
 xksoil(:)   = 1.0_r_2
 fwps(:)     =  casamet%moistavg(:)/soil%ssat(:)
 tsavg(:)    =  casamet%tsoilavg(:)

 ! Custom soil respiration - see Ticket #42
 tsoil(:)    =  tsavg(:)-TKzeroC !tsoil in C
 strf(:)     = 1.0_r_2
 smrf(:)     = 1.0_r_2
 slopt(:)    = 1.0_r_2
 sopt(:)     = 1.0_r_2


  ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
    xktemp(npt)  = casabiome%q10soil(veg%iveg(npt))**(0.1_r_2*(tsavg(npt)-TKzeroC-35.0_r_2))


    xkwater(npt) = ((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
               * ((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
    IF (veg%iveg(npt) == cropland .OR. veg%iveg(npt) == croplnd2) &
               xkwater(npt)=1.0_r_2
    xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) * xktemp(npt) * xkwater(npt)

    IF( .NOT. cable_user%SRF) THEN
        ! Use original function, ELSE Ticket #42
       xksoil(npt)   = casabiome%xkoptsoil(veg%iveg(npt))   * xktemp(npt) * xkwater(npt)
    ELSE
    ! Custom soil respiration - see Ticket #42
    ! Implementing alternative parameterizations
      IF(trim(cable_user%SMRF_NAME)=='CASA-CNP') THEN
         smrf(npt)=xkwater(npt)
      ELSE IF (trim(cable_user%SMRF_NAME)=='SOILN') then
         sopt(npt)=0.92_r_2
         slopt(npt)=wlt(npt)+0.1_r_2          !SLOPT is the lower optimum
         IF (fwps(npt)>sopt(npt)) THEN
           smrf(npt)=0.2_r_2+0.8_r_2*(1.0_r_2-fwps(npt))/(1.0_r_2-sopt(npt))
         ELSE IF(slopt(npt)<=fwps(npt) .AND. fwps(npt)<=sopt(npt)) THEN
           smrf(npt) = 1.0_r_2
         ELSE IF (wlt(npt)<=fwps(npt) .AND. fwps(npt) <slopt(npt)) THEN
           smrf(npt)=0.01_r_2+0.99_r_2*(fwps(npt)-wlt(npt))/(slopt(npt)-wlt(npt))
         ELSE IF (fwps(npt)<wlt(npt)) THEN
           smrf(npt) = 0.01_r_2
         END IF
      ELSE IF (trim(cable_user%SMRF_NAME)=='TRIFFID') THEN
         sopt(npt) = 0.5_r_2 * (1._r_2+wlt(npt))
         IF (fwps(npt) > sopt(npt)) THEN
           smrf(npt) =1.0_r_2-0.8_r_2*(fwps(npt)-sopt(npt))
         ELSE IF (wlt(npt)<fwps(npt) .AND. fwps(npt)<=sopt(npt)) THEN
           smrf(npt)=0.01_r_2+0.8_r_2*((fwps(npt)-wlt(npt))/(sopt(npt)-wlt(npt)))
         ELSE IF (fwps(npt)<wlt(npt)) THEN
           smrf(npt) = 0.2_r_2
         END IF
      ELSE IF (trim(cable_user%SMRF_NAME)=='Trudinger2016') THEN

         if (fwps(npt) .le. wfpswidth0) then
            smrf(npt) = wfpsquad*fwps(npt)**2.0_r_2
         else
            if (fwps(npt) .le. (wfpswidth0+wfpswidth1)) then
               f0 = 0.5_r_2*(1.0_r_2-cos(3.1415_r_2*((wfpsscale1-1)*wfpswidth1)/(wfpswidth1*wfpsscale1)))
               smrf(npt) = max((wfpsquad*fwps(npt)**2.0),  &
                    ((0.5_r_2*(1.0_r_2-cos(3.1415_r_2*(fwps(npt)-wfpswidth0+(wfpsscale1-1)*wfpswidth1)/ &
                    (wfpswidth1*wfpsscale1)))-f0)/(1._r_2-f0)))
            else
               if (fwps(npt) .le. wfpswidth0 + wfpswidth1+wfpswidth2) then
                  smrf(npt) = 1.0_r_2
               else
                  smrf(npt) = 0.5_r_2*(1.0_r_2+cos(3.1415_r_2*(fwps(npt)-wfpswidth0-wfpswidth1-wfpswidth2) &
                       /wfpswidth3))
               endif
            endif
         endif

      END IF

      IF(trim(cable_user%STRF_NAME)=='CASA-CNP') THEN
        strf(npt)=xktemp(npt)
      ELSE if (trim(cable_user%STRF_NAME)=='K1995') THEN
      !Kirschbaum from Kirschbaum 1995, eq (4) in SBB, .66 is to collapse smrf
      !to same area
        strf(npt)=exp(-3.764_r_2+0.204_r_2*tsoil(npt)*(1._r_2-0.5_r_2*tsoil(npt)/36.9_r_2))/0.66_r_2
      ELSE IF (trim(cable_user%STRF_NAME)=='PnET-CN') THEN
        strf(npt)=0.68_r_2*exp(0.1_r_2*(tsoil(npt)-7.1_r_2))/12.64_r_2
      ELSE if (trim(cable_user%STRF_NAME)=='LT1994') THEN
         ! Lloyd & Taylor, Func. Ec. 8, 315, 1994, Eq 11.
         ! Normalised to 10 deg C
         strf(npt)= exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2         &
              / (max(tsoil(npt),-20.0_r_2)+46.02_r_2)))/ &
              (exp(308.56_r_2*(1.0_r_2/56.02_r_2-1.0_r_2         &
              / (10.0_r_2 + 46.02_r_2))))
      ELSE if (trim(cable_user%STRF_NAME)=='Q10') THEN
         ! Q10 normalised to 10 deg C
          strf(npt) = casabiome%q10soil(veg%iveg(npt))**(0.1_r_2*(tsavg(npt)-TKzeroC-10.0_r_2))
      END IF

      IF (trim(cable_user%STRF_NAME)=='DAMM' .and. trim(cable_user%SMRF_NAME)=='DAMM') THEN

         ! Implementation follows Sihi et al., AFM 252 (2018) 155-166.
         ! The key difference is that we use the DAMM model to provide rate constant modifiers
         ! for soil and litter carbon turnover, whereas Sihi et al. apply the model directly
         ! to heterotrophic respiration. Also, we account here for T-dependence of O2 diffusion
         ! and diffusivity of Enzyme in the liquid phase (T-dependences of diffusivities taken from
         ! (Haverd and Cuntz 2010, J. Hyd., 388, (2010), 438-455 ).

         vol_air_content = max(soil%ssat(npt) - casamet%moistavg(npt), 0.0_r_2)

         Dva = 1.67_r_2 * 0.209_r_2 * ((casamet%tsoilavg(npt)/TKzeroC) ** 1.88_r_2 )/ &
               (((10.0_r_2  + TKzeroC) / TKzeroC) ** 1.88_r_2 )
         O2 = Dva * (vol_air_content**(4._r_2/3._r_2))

         Dliq = 3.17_r_2 * ( casamet%moistavg(npt)**3) * exp(-577.0_r_2/(casamet%tsoilavg(npt)-145._r_2))/ &
              exp(-577.0_r_2/((10.0_r_2 + 273.15_r_2)-145._r_2))

         Enz=casabiome%DAMM_EnzPool(veg%iveg(npt))* Dliq

         ! Arrhenius temperature response function, normalised to 10 deg C.
         strf(npt) = exp(-casabiome%DAMM_Ea(veg%iveg(npt))/(8.314e-3_r_2 &
              * (casamet%tsoilavg(npt)))) /                   &
         (exp(-casabiome%DAMM_Ea(veg%iveg(npt))/(8.314e-3_r_2 &
              * 283.0_r_2)))

         ! moisture response function
         smrf(npt) =  (Enz/(Enz+casabiome%DAMM_KMcp(veg%iveg(npt)))) &
              * (O2/(casabiome%DAMM_KMO2(veg%iveg(npt))+O2))

         xksoil(npt) =  (10._r_2**casabiome%DAMM_alpha(veg%iveg(npt)))* smrf(npt) * strf(npt)
         xklitter(npt) = xksoil(npt)

      ELSE
         !xksoil(npt) = casabiome%xkoptsoil(veg%iveg(npt))*strf(npt)*smrf(npt)
         !xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) *strf(npt)*smrf(npt)
         xksoil(npt) = strf(npt)*smrf(npt)
         xklitter(npt) = strf(npt)*smrf(npt)
      ENDIF

!write(67,"(i8,18e16.6)") npt, tsoil(npt), fwps(npt), smrf(npt), strf(npt),xkwater(npt),  xktemp(npt), xklitter(npt)
!write(67,"(i8,18e16.6)") npt, tsoil(npt), strf(npt),  xktemp(npt)

   END IF
  END IF
  END DO

!!$ write(63,"(100e16.6)") xklitter
!!$ write(64,"(100e16.6)") xksoil
!!$ write(65,"(100e16.6)") smrf
!!$ write(66,"(100e16.6)") strf
!!$ write(67,"(100e16.6)") fwps
!!$ write(68,"(100e16.6)") tsoil
END SUBROUTINE casa_xratesoil


SUBROUTINE casa_coeffplant(xkleafcold, xkleafdry, xkleaf, veg, casabiome, casapool, &
     casaflux, casamet, phen)
  ! calculate the plant litter fall rate, litter fall and SOM decomposition rate (1/day)
  ! and the transfer coefficients between different pools
  !
  ! inputs:
  !     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
  !     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
  !     xkleaf(mp):      set to 0.0 during maximal leaf growth phase
  !
  ! outputs:
  !     kplant(mp,mplant):        senescence rate of plant pool (1/day)
  !     fromPtoL(mp,mlitter,mplant): fraction of senesced plant biomass to litter pool (fraction)

  implicit none
  
  real(r_2), dimension(mp), intent(in)    :: xkleafcold, xkleafdry, xkleaf
  type(veg_parameter_type), intent(inout) :: veg       ! vegetation parameters
  type(casa_biome),         intent(inout) :: casabiome
  type(casa_pool),          intent(inout) :: casapool
  type(casa_flux),          intent(inout) :: casaflux
  type(casa_met),           intent(inout) :: casamet
  type(phen_variable),      intent(in)    :: phen

  ! local variables
  real(r_2), dimension(mp,mplant) :: ratioligninton
  integer :: npt

  casaflux%fromPtoL(:,:,:) = 0.0_r_2
  casaflux%kplant(:,:)     = 0.0_r_2
  casaflux%kplant_tot(:,:) = 0.0_r_2

  where (casamet%iveg2/=icewater)
     ! using max function to avoid dividing by zero, ypw 14/may/2008
     ratioLignintoN(:,leaf)  = ( casapool%Cplant(:,leaf) &
          / ( max(1.0e-10_r_2, casapool%Nplant(:,leaf)) * casabiome%ftransNPtoL(veg%iveg(:),leaf) ) ) &
          * casabiome%fracLigninplant(veg%iveg(:),leaf)
     ratioLignintoN(:,froot) = ( casapool%Cplant(:,froot) &
          / ( max(1.0e-10_r_2, casapool%Nplant(:,froot)) * casabiome%ftransNPtoL(veg%iveg(:),froot) ) ) &
          * casabiome%fracLigninplant(veg%iveg(:),froot)

     casaflux%fromPtoL(:,metb,leaf)  = max(0.001_r_2, 0.85_r_2 - 0.018_r_2 * ratioLignintoN(:,leaf))
     casaflux%fromPtoL(:,metb,froot) = max(0.001_r_2, 0.85_r_2 - 0.018_r_2 * ratioLignintoN(:,froot))
     casaflux%fromPtoL(:,str,leaf)   = 1.0_r_2 - casaflux%fromPtoL(:,metb,leaf)
     casaflux%fromPtoL(:,str,froot)  = 1.0_r_2 - casaflux%fromPtoL(:,metb,froot)
     casaflux%fromPtoL(:,cwd,wood)   = 1.0_r_2

     casaflux%kplant(:,leaf)  = casabiome%plantrate(veg%iveg(:),leaf) * xkleaf(:) + xkleafcold(:) + xkleafdry(:)
     casaflux%kplant(:,wood)  = casabiome%plantrate(veg%iveg(:),wood)
     casaflux%kplant(:,froot) = casabiome%plantrate(veg%iveg(:),froot)
  endwhere

  ! set leaf turnover to be the same as froot turnover for C3 & C4 grass
  where (veg%iveg==7 .OR. veg%iveg==6) casaflux%kplant(:,froot) = casaflux%kplant(:,leaf)

  where (casamet%iveg2/=icewater)
     ! total turnover rate includes turnover by fire
     casaflux%kplant_tot(:,leaf)  = casaflux%kplant(:,leaf)  + (1.-casaflux%kplant(:,leaf))  * casaflux%kplant_fire(:,leaf)
     casaflux%kplant_tot(:,froot) = casaflux%kplant(:,froot) + (1.-casaflux%kplant(:,froot)) * casaflux%kplant_fire(:,froot)
     casaflux%kplant_tot(:,wood)  = casaflux%kplant(:,wood)  + (1.-casaflux%kplant(:,wood))  * casaflux%kplant_fire(:,wood)
     ! adjust leaf flux to litter for crop/pasture removal
     casaflux%fromPtoL(:,str,leaf)  = casaflux%fromPtoL(:,str,leaf)  * (1.-casaflux%fharvest)
     casaflux%fromPtoL(:,metb,leaf) = casaflux%fromPtoL(:,metb,leaf) * (1.-casaflux%fharvest)
  endwhere

  ! When glai<glaimin, leaf biomass will not decrease anymore. (Q.Zhang 10/03/2011)
  do npt=1, mp
     if (casamet%glai(npt) .le. casabiome%glaimin(veg%iveg(npt))) then
        casaflux%kplant_tot(npt,leaf) = 0.0_r_2
        casaflux%kplant(npt,leaf)     = 0.0_r_2
     endif
  enddo

END SUBROUTINE casa_coeffplant


SUBROUTINE casa_coeffsoil(xklitter, xksoil, veg, soil, casabiome, casaflux, casamet)
  !  calculate the plant litter fall rate, litter fall and sOM decomposition rate (1/day)
  !  and the transfer coefficients between different pools
  !
  ! inputs:
  !     xk(mp):          modifier of soil litter decomposition rate (dimensionless)
  !
  ! outputs:
  !     klitter(mp,mlitter):      decomposition rate of litter pool (1/day)
  !     ksoil(mp,msoil):          decomposition rate of soil pool (1/day)
  !     fromLtoS(mp,mlitter,msoil):  fraction of decomposed litter to soil (fraction)
  !     fromStoS(mp,msoil,msoil):    fraction of decomposed soil C to another soil pool (fraction)
  !     fromLtoCO2(mp,mlitter):      fraction of decomposed litter emitted as CO2 (fraction)
  !     fromStoCO2(mp,msoil):        fraction of decomposed soil C emitted as Co2 (fraction)

  implicit none
  
  real(r_2), dimension(mp),  intent(in)    :: xklitter, xksoil
  type(veg_parameter_type),  intent(inout) :: veg        ! vegetation parameters
  type(soil_parameter_type), intent(inout) :: soil       ! soil parameters
  type(casa_biome),          intent(inout) :: casabiome
  type(casa_flux),           intent(inout) :: casaflux
  type(casa_met),            intent(inout) :: casamet

  ! local variables
  integer :: j, k, kk, nland             ! i: for plant pool, j for litter, k for soil

  casaflux%fromLtoS(:,:,:) = 0.0_r_2
  casaflux%fromStoS(:,:,:) = 0.0_r_2

  do k=1, msoil
     casaflux%fromStoS(:,k,k) = -1.0_r_2     ! flow from soil to soil
  enddo
  casaflux%fromLtoCO2(:,:) = 0.0_r_2         ! flow from L or S to CO2
  casaflux%fromStoCO2(:,:) = 0.0_r_2

  casaflux%klitter(:,:)     = 0.0_r_2        ! initialize klitter (Q.Zhang 03/03/2011)
  casaflux%klitter_tot(:,:) = 0.0_r_2

  where (casamet%iveg2/=icewater)
    casaflux%klitter(:,metb) = xklitter(:) * casabiome%litterrate(veg%iveg(:),metb)
    casaflux%klitter(:,str)  = xklitter(:) * casabiome%litterrate(veg%iveg(:),str) * &
         exp(-3.0*casabiome%fracLigninplant(veg%iveg(:),leaf))
    casaflux%klitter(:,cwd)  = xklitter(:) * casabiome%litterrate(veg%iveg(:),cwd)

    ! add fire turnover of litter
    casaflux%klitter_tot(:,metb) = casaflux%klitter(:,metb) + (1.-casaflux%klitter(:,metb)) * casaflux%klitter_fire(:,metb)
    casaflux%klitter_tot(:,str)  = casaflux%klitter(:,str)  + (1.-casaflux%klitter(:,str))  * casaflux%klitter_fire(:,str)
    casaflux%klitter_tot(:,cwd)  = casaflux%klitter(:,cwd)  + (1.-casaflux%klitter(:,cwd))  * casaflux%klitter_fire(:,cwd)

    ! soil turnover
    casaflux%ksoil(:,mic)  = xksoil(:) * casabiome%soilrate(veg%iveg(:),mic) * &
         (1.0_r_2 - 0.75 *(soil%silt(:)+soil%clay(:)))
    casaflux%ksoil(:,slow) = xksoil(:) * casabiome%soilrate(veg%iveg(:),slow)
    casaflux%ksoil(:,pass) = xksoil(:) * casabiome%soilrate(veg%iveg(:),pass)
    casaflux%kplab(:)      = xksoil(:) * casabiome%xkplab(casamet%isorder(:))
    casaflux%kpsorb(:)     = xksoil(:) * casabiome%xkpsorb(casamet%isorder(:))
    casaflux%kpocc(:)      = xksoil(:) * casabiome%xkpocc(casamet%isorder(:))

    where (veg%iveg==cropland) ! for cultivated land type
       casaflux%ksoil(:,mic)  = casaflux%ksoil(:,mic)  * 1.25
       casaflux%ksoil(:,slow) = casaflux%ksoil(:,slow) * 1.5
       casaflux%ksoil(:,pass) = casaflux%ksoil(:,pass) * 1.5
    endwhere  !

    where (casaflux%fcrop>0.1) ! 50% increase in soil C turnover in cropping areas
       casaflux%ksoil(:,mic)  = casaflux%ksoil(:,mic) * ((1.-casaflux%fcrop) + 1.5*casaflux%fcrop)
       casaflux%ksoil(:,slow) = casaflux%ksoil(:,slow)* ((1.-casaflux%fcrop) + 1.5*casaflux%fcrop)
       casaflux%ksoil(:,pass) = casaflux%ksoil(:,pass)* ((1.-casaflux%fcrop) + 1.5*casaflux%fcrop)
    endwhere
    
    ! flow from litter to soil
    casaflux%fromLtoS(:,mic,metb) = 0.45                                                   ! metb -> mic
    casaflux%fromLtoS(:,mic,str)  = 0.45*(1.0_r_2-casabiome%fracLigninplant(veg%iveg(:),leaf)) ! str  -> mic
    casaflux%fromLtoS(:,slow,str) = 0.7 * casabiome%fracLigninplant(veg%iveg(:),leaf)      ! str  -> slow
    casaflux%fromLtoS(:,mic,cwd)  = 0.40*(1.0_r_2-casabiome%fracLigninplant(veg%iveg(:),wood)) ! cwd  -> fmic
    casaflux%fromLtoS(:,slow,cwd) = 0.7 * casabiome%fracLigninplant(veg%iveg(:),wood)      ! cwd  -> slow

    ! set the following two backflow to set (see Bolker 199x)
    !   casaflux%fromStoS(:,mic,slow)  = 0.45_r_2 * (0.997_r_2 - 0.009_r_2 *soil%clay(:))
    !   casaflux%fromStoS(:,mic,pass)  = 0.45_r_2
    casaflux%fromStoS(:,slow,mic)  = (0.85_r_2 - 0.68_r_2 * (soil%clay(:)+soil%silt(:))) * (0.997_r_2 - 0.032_r_2*soil%clay(:))
    casaflux%fromStoS(:,pass,mic)  = (0.85_r_2 - 0.68_r_2 * (soil%clay(:)+soil%silt(:))) * (0.003_r_2 + 0.032_r_2*soil%clay(:))
    casaflux%fromStoS(:,pass,slow) = 0.45_r_2 * (0.003_r_2 + 0.009_r_2 * soil%clay(:))
  endwhere ! /= icewater

  do nland=1, mp
     if (casamet%iveg2(nland)/=icewater) then
        do j=1, mlitter
           do k=1, msoil
              casaflux%fromLtoCO2(nland,j) = casaflux%fromLtoCO2(nland,j) + casaflux%fromLtoS(nland,k,j)
           enddo
           casaflux%fromLtoCO2(nland,j) = 1.0_r_2 - casaflux%fromLtoCO2(nland,j)
        enddo
        do k=1, msoil
           do kk=1, msoil
              casaflux%fromStoCO2(nland,k) = casaflux%fromStoCO2(nland,k) + casaflux%fromStoS(nland,kk,k)
           enddo
        enddo
        casaflux%fromStoCO2(nland,:) = -casaflux%fromStoCO2(nland,:)
     endif ! /= icewater
  enddo ! nland

END SUBROUTINE casa_coeffsoil


! modified by ypw following Chris Lu 5/nov/2012
SUBROUTINE casa_delplant(veg, casabiome, casapool, casaflux, casamet, &
     cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd,  &
     nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd,  &
     pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd)

  !  calculate the chnage in plant C, N and P pools
  !  uptake of N and P will be computed in casa_uptake
  !  labile C pool will be computed casa_labile

  implicit none

  type(veg_parameter_type), intent(inout) :: veg       ! vegetation parameters
  type(casa_biome),         intent(inout) :: casabiome
  type(casa_pool),          intent(inout) :: casapool
  type(casa_flux),          intent(inout) :: casaflux
  type(casa_met),           intent(inout) :: casamet
  ! added by ypwang following Chris Lu 5/nov/2012
  real(r_2), dimension(mp), intent(out)   :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd, &
       nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd, &
       pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd

  integer   :: npt, nL, nP
  real(r_2) :: Ygrow, ratioPNplant

  ! casa
  casaflux%FluxCtolitter         = 0.0_r_2
  casaflux%FluxNtolitter         = 0.0_r_2
  casaflux%FluxPtolitter         = 0.0_r_2
  ! fire
  casaflux%fluxCtoCO2_plant_fire = 0.0_r_2
  casaflux%fluxNtoAtm_fire       = 0.0_r_2
  ! 13C
  casaflux%FluxFromPtoL          = 0.0_r_2
  casaflux%FluxFromPtoCO2        = 0.0_r_2
  casaflux%FluxFromPtoHarvest    = 0.0_r_2

  ! added by ypwang following Chris Lu 5/nov/2012
  cleaf2met = 0.0_r_2
  cleaf2str = 0.0_r_2
  croot2met = 0.0_r_2
  croot2str = 0.0_r_2
  cwood2cwd = 0.0_r_2

  nleaf2met = 0.0_r_2
  nleaf2str = 0.0_r_2
  nroot2met = 0.0_r_2
  nroot2str = 0.0_r_2
  nwood2cwd = 0.0_r_2

  pleaf2met = 0.0_r_2
  pleaf2str = 0.0_r_2
  proot2met = 0.0_r_2
  proot2str = 0.0_r_2
  pwood2cwd = 0.0_r_2

  !MPI
  do npt=1, mp

     if (casamet%iveg2(npt) /= icewater) then
        casapool%dcplantdt(npt,:)  =  casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)     &
             - casaflux%kplant_tot(npt,:)  * casapool%cplant(npt,:)

        !! vh_js !!
        ! adjust turnover and autotrophic respiration to avoid negative stocks
        ! Ticket#108
        where ( ( (casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) .lt. 0.0_r_2) &
             .or. ( (casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) &
             .lt. 0.5_r_2 * casapool%cplant(npt,2:3) ) )
           casaflux%kplant_tot(npt,2:3) = 0.0_r_2
           casaflux%kplant(npt,2:3)     = 0.0_r_2
           casaflux%crmplant(npt,2:3)   = 0.0_r_2
        endwhere
        if ( any((casapool%dcplantdt(npt,:)*deltpool + casapool%cplant(npt,:)) .lt. 0.0_r_2) ) then
           casaflux%kplant_tot(npt,1) = 0.0_r_2
           casaflux%kplant(npt,1)     = 0.0_r_2
           casaflux%crmplant(npt,1)   = min(casaflux%crmplant(npt,1), 0.5_r_2*casaflux%cgpp(npt))
        endif

        ! revise turnover and NPP and dcplantdt to reflect above adjustments
        casaflux%Cplant_turnover(npt,:) = casaflux%kplant_tot(npt,:)  * casapool%cplant(npt,:)
        if (any((casapool%dcplantdt(npt,:)*deltpool + casapool%cplant(npt,:)).lt. 0.0_r_2) &
          .or. any((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) &
                  .lt. 0.5_r_2 * casapool%cplant(npt,2:3) )) then
           ratioPNplant = 0.0_r_2
           if (casapool%Nplant(npt,leaf)>0.0_r_2) &
                ratioPNplant = casapool%Pplant(npt,leaf)/(casapool%Nplant(npt,leaf)+ 1.0e-10_r_2)

           Ygrow = 0.65_r_2+0.2_r_2*ratioPNplant/(ratioPNplant+1.0_r_2/15.0_r_2)
           IF ((casaflux%Cgpp(npt)-SUM(casaflux%crmplant(npt,:)))>0.0_r_2) THEN
              ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
              casaflux%crgplant(npt)  = (1.0_r_2-Ygrow)* max(0.0_r_2,casaflux%Cgpp(npt)- &
                   SUM(casaflux%crmplant(npt,:)))
           ELSE
              casaflux%crgplant(npt) = 0.0_r_2
           ENDIF
        endif
        ! recalc in any case so that all consistent
        casaflux%Cnpp(npt) = casaflux%Cgpp(npt) - sum(casaflux%Crmplant(npt,:)) - casaflux%Crgplant(npt) - &
             casaflux%fracClabile(npt) * casaflux%Cgpp(npt)
        casapool%dcplantdt(npt,:) = casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:) - &
             casaflux%kplant_tot(npt,:) * casapool%cplant(npt,:)
        !! vh_js !! end of adjustments to avoid negative stocks Ticket#108

        ! change here made by ypw on 26 august 2011
        ! calculate fraction C to labile pool as a fraction of gpp, not npp
        ! casapool%dClabiledt(npt) = casaflux%Cnpp(npt) * casaflux%fracClabile(npt)
        casapool%dClabiledt(npt) = casaflux%Cgpp(npt) * casaflux%fracClabile(npt) - casaflux%clabloss(npt)

        ! added by ypwang 5/nov/2012
        cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
        cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
        croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
        croot2str(npt) = casaflux%fromPtoL(npt,str,froot)  * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
        cwood2cwd(npt) = casaflux%fromPtoL(npt,cwd,wood)   * casaflux%kplant(npt,wood)  * casapool%cplant(npt,wood)

        ! add fire component to above fluxes
        ! print*, 'CN01 ', npt, cleaf2met(npt), cleaf2str(npt), croot2met(npt)
        print*, 'CN02 ', npt, metb, casaflux%fromPtoL_fire(npt,:,:)
        ! print*, 'CN03 ', casaflux%kplant(npt,:)
        ! print*, 'CN04 ', casaflux%kplant_fire(npt,:)
        ! print*, 'CN05 ', casapool%cplant(npt,:)
        ! print*, 'CN06 ', npt, croot2str(npt), cwood2cwd(npt)
        cleaf2met(npt) = cleaf2met(npt) + casaflux%fromPtoL_fire(npt,metb,leaf) &
             * (1.0_r_2 - casaflux%kplant(npt,leaf))  * casaflux%kplant_fire(npt,leaf)  * casapool%cplant(npt,leaf)
        ! print*, 'CN07 ', casaflux%fromPtoL_fire(npt,str,leaf)
        ! print*, 'CN08 ', (1.0_r_2 - casaflux%kplant(npt,leaf))
        ! print*, 'CN09 ', casaflux%fromPtoL_fire(npt,str,leaf) &
        !      * (1.0_r_2 - casaflux%kplant(npt,leaf))
        ! print*, 'CN10 ', casaflux%fromPtoL_fire(npt,str,leaf) &
        !      * (1.0_r_2 - casaflux%kplant(npt,leaf))  * casaflux%kplant_fire(npt,leaf)
        ! print*, 'CN11 ', casaflux%fromPtoL_fire(npt,str,leaf) &
        !      * (1.0_r_2 - casaflux%kplant(npt,leaf))  * casaflux%kplant_fire(npt,leaf)  * casapool%cplant(npt,leaf)
        cleaf2str(npt) = cleaf2str(npt) + casaflux%fromPtoL_fire(npt,str,leaf) &
             * (1.0_r_2 - casaflux%kplant(npt,leaf))  * casaflux%kplant_fire(npt,leaf)  * casapool%cplant(npt,leaf)
        croot2met(npt) = croot2met(npt) + casaflux%fromPtoL_fire(npt,metb,froot) &
             * (1.0_r_2 - casaflux%kplant(npt,froot)) * casaflux%kplant_fire(npt,froot) * casapool%cplant(npt,froot)
        croot2str(npt) = croot2str(npt) + casaflux%fromPtoL_fire(npt,str,froot) &
             * (1.0_r_2 - casaflux%kplant(npt,froot)) * casaflux%kplant_fire(npt,froot) * casapool%cplant(npt,froot)
        cwood2cwd(npt) = cwood2cwd(npt) + casaflux%fromPtoL_fire(npt,cwd,wood) &
             * (1.0_r_2 - casaflux%kplant(npt,wood))  * casaflux%kplant_fire(npt,wood)  * casapool%cplant(npt,wood)

        ! fire flux to atmosphere from burnt plant material
        casaflux%fluxCtoCO2_plant_fire(npt) = 0.0_r_2
        casaflux%FluxFromPtoCO2(npt,:) = 0.0_r_2
        do nP=1,mplant
           casaflux%FluxFromPtoCO2(npt,nP) = &
                (1.0_r_2-casaflux%kplant(npt,nP)) * casaflux%kplant_fire(npt,nP) * &
                casapool%cplant(npt,nP) * &
                ( 1.0_r_2 - SUM(casaflux%fromPtoL_fire(npt,:,nP))) 
           
           casaflux%fluxCtoCO2_plant_fire(npt) = casaflux%fluxCtoCO2_plant_fire(npt) + &
                casaflux%FluxFromPtoCO2(npt,nP)

        enddo

        ! Crop Harvest Flux
        casaflux%Charvest(npt) = casaflux%Charvest(npt) + &
             casaflux%fharvest(npt) * casaflux%kplant(npt,leaf) * casapool%cplant(npt,leaf)

        ! between-pool fluxes for 13CO2
        casaflux%FluxFromPtoL(npt,leaf,metb)  = cleaf2met(npt)
        casaflux%FluxFromPtoL(npt,froot,metb) = croot2met(npt)
        casaflux%FluxFromPtoL(npt,wood,metb)  = 0.0_r_2

        casaflux%FluxFromPtoL(npt,leaf,str)  = cleaf2str(npt)
        casaflux%FluxFromPtoL(npt,froot,str) = croot2str(npt)
        casaflux%FluxFromPtoL(npt,wood,str)  = 0.0_r_2

        casaflux%FluxFromPtoL(npt,leaf,cwd)  = 0.0_r_2
        casaflux%FluxFromPtoL(npt,froot,cwd) = 0.0_r_2
        casaflux%FluxFromPtoL(npt,wood,cwd)  = cwood2cwd(npt)

        casaflux%FluxFromPtoHarvest(npt) = casaflux%FluxFromPtoHarvest(npt) + &
             casaflux%fharvest(npt) * casaflux%kplant(npt,leaf) * casapool%cplant(npt,leaf)
        ! if (casaflux%FluxFromPtoHarvest(npt) > 0.) &
        !      print*, 'We have Harvest 02 ', npt, casaflux%FluxFromPtoHarvest(npt)

        ! Nitrogen
        IF (icycle > 1) THEN

           IF(casaflux%fracNalloc(npt,leaf)==0.0_r_2) THEN
              casapool%dNplantdt(npt,leaf)  = - casaflux%kplant_tot(npt,leaf) * casapool%Nplant(npt,leaf)
           else
              casapool%dNplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf) &
                   * casabiome%ftransNPtoL(veg%iveg(npt),leaf)
              ! add turnover due to fire
              casapool%dNplantdt(npt,leaf)  =  casapool%dNplantdt(npt,leaf) - &
                   (1.0_r_2 - casaflux%kplant(npt,leaf))*casaflux%kplant_fire(npt,leaf) * casapool%Nplant(npt,leaf)

           ENDIF

           IF (casamet%lnonwood(npt)==0) THEN
              casapool%dNplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Nplant(npt,wood) &
                   * casabiome%ftransNPtoL(veg%iveg(npt),wood)

               ! add turnover due to fire
              casapool%dNplantdt(npt,wood)  =  casapool%dNplantdt(npt,wood) - &
                   (1.0_r_2 - casaflux%kplant(npt,wood))*casaflux%kplant_fire(npt,wood) * casapool%Nplant(npt,wood)
           ELSE
              casapool%dNplantdt(npt,wood) = 0.0_r_2
           ENDIF

           casapool%dNplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Nplant(npt,froot) &
                * casabiome%ftransNPtoL(veg%iveg(npt),froot)

           ! add turnover due to fire
           casapool%dNplantdt(npt,froot)  =  casapool%dNplantdt(npt,froot) - &
                (1.0_r_2 - casaflux%kplant(npt,froot))*casaflux%kplant_fire(npt,froot) * casapool%Nplant(npt,froot)

           nleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                * casapool%cplant(npt,leaf)       * ratioNCstrfix
           nroot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                * casapool%cplant(npt,froot)      * ratioNCstrfix

           ! add N flux from leaves and roots to structural litter due to fire
           nleaf2str(npt) = nleaf2str(npt) + casaflux%fromPtoL_fire(npt,str,leaf) * (1. - casaflux%kplant(npt,leaf)) *  &
                casaflux%kplant_fire(npt,leaf)* casapool%cplant(npt,leaf)* ratioNCstrfix
           nroot2str(npt) = nroot2str(npt) + casaflux%fromPtoL_fire(npt,str,froot) * (1. - casaflux%kplant(npt,froot)) *  &
                casaflux%kplant_fire(npt,froot)* casapool%cplant(npt,froot)* ratioNCstrfix

           ! flux of leaf N to met litter: need to deduct harvest losses
           nleaf2met(npt) = - casapool%dNplantdt(npt,leaf)   - nleaf2str(npt) &
                - casaflux%kplant(npt,leaf)  * casapool%nplant(npt,leaf) *  casaflux%fharvest(npt)

           nroot2met(npt) = - casapool%dNplantdt(npt,froot) - nroot2str(npt)

           nwood2cwd(npt) = -casapool%dNplantdt(npt,wood)

           ! Nitrogen lost to harvest
           casaflux%Nharvest(npt) = casaflux%Nharvest(npt) + &
                casaflux%kplant(npt,leaf)  * casapool%nplant(npt,leaf) *  casaflux%fharvest(npt)

           ! Nitrogen lost to fire
           casaflux%fluxNtoAtm_fire(npt) =  (1.0_r_2 - casaflux%kplant(npt,leaf))* &
                casaflux%kplant_fire(npt,leaf) * casapool%Nplant(npt,leaf) &
                * (1.0_r_2 - casaflux%fromPtoL_fire(npt,str,leaf)) + &
                (1.0_r_2 - casaflux%kplant(npt,froot))*casaflux%kplant_fire(npt,froot) * &
                casapool%Nplant(npt,froot) &
                * (1.0_r_2 - casaflux%fromPtoL_fire(npt,str,froot)) + &
                (1.0_r_2 - casaflux%kplant(npt,wood))*casaflux%kplant_fire(npt,wood) * &
                casapool%Nplant(npt,wood) &
                * (1.0_r_2 - casaflux%fromPtoL_fire(npt,str,wood))
        ENDIF

        ! Phosphorus
        IF (icycle >2) THEN

           IF(casaflux%fracPalloc(npt,leaf)==0.0_r_2) THEN
              casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf)
           else
              casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf) &
                   * casabiome%ftransPPtoL(veg%iveg(npt),leaf)
           ENDIF

           casapool%dPplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Pplant(npt,froot) &
                * casabiome%ftransPPtoL(veg%iveg(npt),froot)

           casapool%dPplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Pplant(npt,wood) &
                * casabiome%ftransPPtoL(veg%iveg(npt),wood)

           pleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix
           proot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix


           ! add P loss by fire
           casapool%dPplantdt(npt,leaf) = casapool%dPplantdt(npt,leaf) - &
                (1.0_r_2 - casaflux%kplant(npt,leaf) ) * casaflux%kplant_fire(npt,leaf) * casapool%Pplant(npt,leaf)

           casapool%dPplantdt(npt,froot)  = casapool%dPplantdt(npt,froot) - &
                (1.0_r_2 - casaflux%kplant(npt,froot) ) * casaflux%kplant_fire(npt,froot) * casapool%Pplant(npt,froot)

           casapool%dPplantdt(npt,wood)  = casapool%dPplantdt(npt,wood) - &
                (1.0_r_2 - casaflux%kplant(npt,wood) ) * casaflux%kplant_fire(npt,wood) * casapool%Pplant(npt,wood)

           ! pleaf2str(npt) = pleaf2str(npt) + casaflux%fromPtoL_fire(npt,str,leaf) * casaflux%kplant_fire(npt,leaf)*  &
           !      (1. - casaflux%kplant(npt,leaf)) * ratioNCstrfix/ratioNPstrfix * casapool%cplant(npt,leaf)
           
           ! proot2str(npt) = proot2str(npt) + casaflux%fromPtoL_fire(npt,str,froot) * casaflux%kplant_fire(npt,froot)*  &
           !      (1. - casaflux%kplant(npt,froot)) * ratioNCstrfix/ratioNPstrfix * casapool%cplant(npt,froot)

           ! assume all plant phosphorous release by fire goes to litter
           pleaf2str(npt) = pleaf2str(npt) +  casaflux%kplant_fire(npt,leaf)*  &
                (1.0_r_2 - casaflux%kplant(npt,leaf)) * ratioNCstrfix/ratioNPstrfix * casapool%cplant(npt,leaf)


           proot2str(npt) = proot2str(npt) + casaflux%kplant_fire(npt,froot)*  &
                (1.0_r_2 - casaflux%kplant(npt,froot)) * ratioNCstrfix/ratioNPstrfix * casapool%cplant(npt,froot)

           ! end Ploss by fire

           ! P lost to harvest
           casaflux%Pharvest(npt) = casaflux%kplant(npt,leaf) * casaflux%fharvest(npt) * &
                ratioNCstrfix/ratioNPstrfix * casapool%cplant(npt,froot)


           pleaf2met(npt) = -casapool%dPplantdt(npt,leaf)  - pleaf2str(npt) - casaflux%Pharvest(npt)
           proot2met(npt) = -casapool%dPplantdt(npt,froot) - proot2str(npt)
           pwood2cwd(npt) = -casapool%dPplantdt(npt,wood)
        ENDIF
        
        do nL=1, mlitter
           do nP=1, mplant
              casaflux%FluxCtolitter(npt,nL) = casaflux%FluxCtolitter(npt,nL) &
                   + casaflux%fromPtoL(npt,nL,nP) * casaflux%kplant(npt,nP) * &
                   casapool%cplant(npt,nP) &
                   ! inputs from fire
                   + casaflux%kplant_fire(npt,nP)*(1.0_r_2 - casaflux%kplant(npt,nP)) * &
                   casapool%cplant(npt,nP)* casaflux%fromPtoL_fire(npt,nL,nP)
              
           enddo
        enddo

        
        ! Nitrogen
        IF (icycle > 1) THEN
           ! casaflux%FluxNtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
           !                         * casapool%cplant(npt,leaf)       * ratioNCstrfix              &
           !                         + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
           !                         * casapool%cplant(npt,froot)      * ratioNCstrfix

           !vh! to avoid -ve Nitrogen pools Ticket#108
           ! vh ! need to include fire fluxes here?
           casaflux%FluxNtolitter(npt,str) = min(casaflux%fromPtoL(npt,str,leaf) * &
                casaflux%kplant(npt,leaf)  &
                * casapool%cplant(npt,leaf)       * ratioNCstrfix &
                , -casapool%dNplantdt(npt,leaf))             &
                + min(casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                * casapool%cplant(npt,froot)      * ratioNCstrfix &
                , -casapool%dNplantdt(npt,froot))

           ! casaflux%FluxNtolitter(npt,str) =  min(nleaf2str(npt),  -casapool%dNplantdt(npt,leaf)) &
           !      + min(nroot2str(npt), -casapool%dNplantdt(npt,froot))

           casaflux%FluxNtolitter(npt,metb) = - casapool%dNplantdt(npt,leaf)-casapool%dNplantdt(npt,froot) &
                - casaflux%FluxNtolitter(npt,str)
           casaflux%FluxNtolitter(npt,CWD) = -casapool%dNplantdt(npt,wood)

           ! adding N uptake
           casapool%dNplantdt(npt,:) = casapool%dNplantdt(npt,:) &
                + casaflux%Nminuptake(npt)*casaflux%fracNalloc(npt,:)
           !       casapool%Nsoilmin(npt)    = casapool%Nsoilmin(npt) - casaflux%Nminuptake(npt) *deltpool
        ENDIF !end "icycle >1"


        ! Phosphorus
        IF(icycle>2) THEN
           ! casaflux%FluxPtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
           !      * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix        &
           !      + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
           !      * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix

           casaflux%FluxPtolitter(npt,str) = min(pleaf2str(npt),  -casapool%dPplantdt(npt,leaf)) &
                + min(proot2str(npt), -casapool%dPplantdt(npt,froot))


           casaflux%FluxPtolitter(npt,metb) = -casapool%dPplantdt(npt,leaf)-casapool%dPplantdt(npt,froot) &
                - casaflux%FluxPtolitter(npt,str)
           casaflux%FluxPtolitter(npt,CWD) = -casapool%dPplantdt(npt,wood)
           ! add P uptake
           casapool%dPplantdt(npt,:) = casapool%dPplantdt(npt,:) &
                + casaflux%Plabuptake(npt)*casaflux%fracPalloc(npt,:)
           !       casapool%Psoillab(npt)    = casapool%Psoillab(npt) - casaflux%Plabuptake(npt) * deltpool
        ENDIF ! icycle > 2

     endif ! /= icewater

  enddo ! npt=1,mp
  
  ! npt=2

END SUBROUTINE casa_delplant


SUBROUTINE casa_delsoil(veg, casapool, casaflux, casamet, casabiome)

  ! calculate changes in litter and soil pools

  implicit none
  
  type(veg_parameter_type), intent(inout) :: veg        ! vegetation parameters
  type(casa_pool),          intent(inout) :: casapool
  type(casa_flux),          intent(inout) :: casaflux
  type(casa_met),           intent(inout) :: casamet
  type(casa_biome),         intent(inout) :: casabiome

  ! local variables
  real(r_2), dimension(mp) :: xdplabsorb, fluxptase
  integer   :: j, jj, k, kk, kkk, nL, nS, nSS, nland
  real(r_2) :: iflux

  ! casa
  casaflux%fluxCtoCO2   = 0.0_r_2
  casaflux%fluxCtosoil  = 0.0_r_2
  casaflux%fluxNtosoil  = 0.0_r_2
  casaflux%fluxPtosoil  = 0.0_r_2
  casaflux%Crsoil       = 0.0_r_2 ! initialization added by BP jul2010
  casapool%dClitterdt   = 0.0_r_2
  casapool%dCsoildt     = 0.0_r_2

  casapool%dNlitterdt   = 0.0_r_2
  casapool%dNsoildt     = 0.0_r_2
  casapool%dNsoilmindt  = 0.0_r_2
  casaflux%Nsmin        = 0.0_r_2
  casaflux%Nsimm        = 0.0_r_2
  casaflux%Nsnet        = 0.0_r_2
  casaflux%Nminloss     = 0.0_r_2
  casaflux%Nminleach    = 0.0_r_2
  casaflux%Nlittermin   = 0.0_r_2

  casapool%dPlitterdt   = 0.0_r_2
  casapool%dPsoildt     = 0.0_r_2
  casapool%dPsoillabdt  = 0.0_r_2
  casapool%dPsoilsorbdt = 0.0_r_2
  casapool%dPsoiloccdt  = 0.0_r_2
  casaflux%Psmin        = 0.0_r_2
  casaflux%Psimm        = 0.0_r_2
  casaflux%Psnet        = 0.0_r_2
  casaflux%Pleach       = 0.0_r_2
  casaflux%Ploss        = 0.0_r_2
  casaflux%Plittermin   = 0.0_r_2
  fluxptase             = 0.0_r_2
  
  !13C
  casaflux%fluxCtoCO2_litter_fire = 0.0_r_2
  casaflux%FluxFromLtoS           = 0.0_r_2
  casaflux%FluxFromStoS           = 0.0_r_2
  casaflux%FluxFromStoCO2         = 0.0_r_2
  casaflux%FluxFromLtoCO2         = 0.0_r_2

  do nland=1, mp
     if (casamet%iveg2(nland) /= icewater) then

        if (icycle > 1) then
           !vh! set klitter to zero where Nlitter will go -ve
           ! (occurs occasionally for metabolic litter pool) Ticket#108
           where (casaflux%klitter(nland,:) * max(0.0_r_2,casapool%Nlitter(nland,:)).gt. &
                casapool%Nlitter(nland,:)+casaflux%fluxNtolitter(nland,:)) &
                casaflux%klitter(nland,:) = 0.0_r_2
        endif

        if (icycle > 2) then
           !vh! set klitter to zero where Plitter will go -ve
           ! (occurs occasionally for metabolic litter pool) Ticket#108
           where (casaflux%klitter(nland,:) * max(0.0_r_2,casapool%Plitter(nland,:)).gt. &
                casapool%Plitter(nland,:)+casaflux%fluxPtolitter(nland,:)) &
                casaflux%klitter(nland,:) = 0.0_r_2
        endif

        ! fire flux to atmosphere from burnt litter material
        do nL=1, mlitter
           iflux = (1.-casaflux%klitter(nland,nL)) * casaflux%klitter_fire(nland,nL) * &
                casapool%clitter(nland,nL)
           !if (iflux > 0.) print*, 'We have fire 01 ', nL, iflux,  casaflux%klitter_fire(nland,nL)
           casaflux%fluxCtoCO2_litter_fire(nland) = casaflux%fluxCtoCO2_litter_fire(nland) + iflux
           casaflux%FluxFromLtoCO2(nland,nL)      = casaflux%FluxFromLtoCO2(nland,nL)      + iflux
           
           iflux = casaflux%fromLtoCO2(nland,nL) * casaflux%klitter(nland,nL) * &
                casapool%clitter(nland,nL)
           casaflux%fluxCtoCO2(nland)        = casaflux%fluxCtoCO2(nland)        + iflux
           casaflux%FluxFromLtoCO2(nland,nL) = casaflux%FluxFromLtoCO2(nland,nL) + iflux
        enddo

        do nS=1, msoil
           do nL=1, mlitter
              iflux = casaflux%fromLtoS(nland,nS,nL) * casaflux%klitter(nland,nL) * casapool%clitter(nland,nL)
              casaflux%fluxCtosoil(nland,nS)     = casaflux%fluxCtosoil(nland,nS)     + iflux
              casaflux%FluxFromLtoS(nland,nL,nS) = casaflux%FluxFromLtoS(nland,nL,nS) + iflux
           enddo
           do nSS=1, msoil
              if (nSS /= nS) then
                 iflux = casaflux%fromStoS(nland,nS,nSS) * casaflux%ksoil(nland,nSS) * casapool%csoil(nland,nSS)
                 casaflux%fluxCtosoil(nland,nS)      = casaflux%fluxCtosoil(nland,nS)      + iflux
                 casaflux%FluxFromStoS(nland,nSS,nS) = casaflux%FluxFromStoS(nland,nSS,nS) + iflux
              endif
           enddo
           iflux = casaflux%fromStoCO2(nland,nS) * casaflux%ksoil(nland,nS) * casapool%csoil(nland,nS)
           casaflux%fluxCtoCO2(nland)        = casaflux%fluxCtoCO2(nland)        + iflux
           casaflux%FluxFromStoCO2(nland,nS) = casaflux%FluxFromStoCO2(nland,nS) + iflux
        enddo

        ! Nitrogen
        IF (icycle>1) THEN
           DO j=1,mlitter
              casaflux%Nlittermin(nland) = casaflux%Nlittermin(nland) &
                   + casaflux%klitter(nland,j) * casapool%Nlitter(nland,j)
           ENDDO
           DO k=1,msoil
              casaflux%Nsmin(nland)   = casaflux%Nsmin(nland)   &
                   + casaflux%ksoil(nland,k)   * casapool%Nsoil(nland,k)
           ENDDO    !gross mineralisation

           DO kk=1,msoil
              DO jj=1,mlitter    ! immobilisation from litter to soil
                 casaflux%Nsimm(nland) = casaflux%Nsimm(nland) &
                      - casaflux%fromLtoS(nland,kk,jj) &
                      * casaflux%klitter(nland,jj)     &
                      * casapool%Clitter(nland,jj)     &
                      * casapool%ratioNCsoilnew(nland,kk)
              ENDDO
              DO kkk=1,msoil      ! immobilisation from soil to soil
                 IF(kkk.ne.kk) THEN
                    casaflux%Nsimm(nland) = casaflux%Nsimm(nland) &
                         - casaflux%fromStoS(nland,kk,kkk)  &
                         * casaflux%ksoil(nland,kkk) &
                         * casapool%Csoil(nland,kkk) &
                         * casapool%ratioNCsoilnew(nland,kk)
                 ENDIF
              ENDDO
           ENDDO  ! immobilization
           
           casaflux%Nsnet(nland)=casaflux%Nlittermin(nland) &
                +casaflux%Nsmin(nland)   &
                +casaflux%Nsimm(nland)
           ! net mineralization
           IF(casapool%Nsoilmin(nland)>2.0.AND.casamet%tsoilavg(nland)>273.12) THEN
              casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                   * MAX(0.0_r_2,casaflux%Nsnet(nland))
              casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                   * MAX(0.0_r_2,casapool%Nsoilmin(nland))
           ELSE
              casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                   * MAX(0.0_r_2,casaflux%Nsnet(nland)) &
                   * MAX(0.0_r_2,casapool%Nsoilmin(nland)/2.0)
              casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                   * MAX(0.0_r_2,casapool%Nsoilmin(nland)) &
                   * MAX(0.0_r_2,casapool%Nsoilmin(nland)/2.0)
           ENDIF
           DO k=1,msoil
              DO j=1,mlitter
                 casaflux%FluxNtosoil(nland,k) =  casaflux%FluxNtosoil(nland,k)  &
                      + casaflux%fromLtoS(nland,k,j) &
                      * casaflux%klitter(nland,j)    &
                      * casapool%Clitter(nland,j)    &
                      * casapool%ratioNCsoilnew(nland,k)
              ENDDO  ! end of "j"
              DO kk=1,msoil
                 IF(kk.ne.k) THEN
                    casaflux%FluxNtosoil(nland,k) = casaflux%FluxNtosoil(nland,k)  &
                         + casaflux%fromStoS(nland,k,kk) &
                         * casaflux%ksoil(nland,kk)      &
                         * casapool%Csoil(nland,kk)      &
                         * casapool%ratioNCsoilnew(nland,k)
                 ENDIF
              ENDDO ! end of "kk"
           ENDDO    ! end of "k"
           
        ENDIF ! end of icycle > 1

        ! Phosphorus
        IF (icycle >2) THEN
           DO j=1,mlitter
              casaflux%Plittermin(nland) = casaflux%Plittermin(nland) &
                   + casaflux%klitter(nland,j) * casapool%Plitter(nland,j)
           ENDDO
           DO k=1,msoil
              casaflux%Psmin(nland)   = casaflux%Psmin(nland)   &
                   + casaflux%ksoil(nland,k)   * casapool%Psoil(nland,k)
           ENDDO    !gross mineralisation

           DO kk=1,msoil
              DO jj=1,mlitter    ! immobilisation from litter to soil
                 casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                      - casaflux%fromLtoS(nland,kk,jj) &
                      * casaflux%klitter(nland,jj)     &
                      * casapool%Nlitter(nland,jj)     &
                      /casapool%ratioNPsoil(nland,kk)
                 !                                     * casapool%ratioPCsoil(nland,kk)/casapool%ratioNCsoil(nland,kk)
              ENDDO
              DO kkk=1,msoil      ! immobilisation from soil to soil
                 IF(kkk.ne.kk) THEN
                    casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                         - casaflux%fromStoS(nland,kk,kkk)  &
                         * casaflux%ksoil(nland,kkk) &
                         * casapool%Nsoil(nland,kkk) &
                         /casapool%ratioNPsoil(nland,kk)
                    !                                        * casapool%ratioPCsoil(nland,kk)/casapool%ratioNCsoil(nland,kk)
                 ENDIF
              ENDDO
           ENDDO  ! immobilization
           
           casaflux%Psnet(nland)=casaflux%Plittermin(nland) &
                +casaflux%Psmin(nland)   &
                +casaflux%Psimm(nland)
           ! net mineralization
           !      casaflux%Pleach(nland)  =  (1.0e-4) &
           !                                 * max(0.0_r_2,casapool%Psoillab(nland))           
           casaflux%Pleach(nland)  =  casaflux%fPleach(nland) &
                * max(0.0_r_2,casapool%Psoillab(nland))
           
           DO k=1,msoil
              DO j=1,mlitter
                 casaflux%FluxPtosoil(nland,k) =  casaflux%FluxPtosoil(nland,k)  &
                      + casaflux%fromLtoS(nland,k,j) &
                      * casaflux%klitter(nland,j)    &
                      * casapool%Nlitter(nland,j)    &
                      /casapool%ratioNPsoil(nland,k)
                 !                                 * casapool%ratioPCsoil(nland,k)/casapool%ratioNCsoil(nland,k)
              ENDDO  ! end of "j"
              DO kk=1,msoil
                 IF(kk.ne.k) THEN
                    !               casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
                    !                                    + casaflux%fromStoS(nland,k,kk) &
                    !                                    * casaflux%ksoil(nland,kk)      &
                    !                                    * casapool%Csoil(nland,kk)      &
                    !                                    * casapool%ratioPCsoil(nland,k)
                    casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
                         + casaflux%fromStoS(nland,k,kk) &
                         * casaflux%ksoil(nland,kk)      &
                         * casapool%Nsoil(nland,kk)      &
                         /casapool%ratioNPsoil(nland,k)
                    !                                    * casapool%ratioPCsoil(nland,k)/casapool%ratioNCsoil(nland,k)
                    
                    
                 ENDIF
              ENDDO ! end of "kk"
           ENDDO    ! end of "k"
           ! need to account for flow from sorbed to occluded pool
        ENDIF ! end of icycle > 2
        
     endif  ! end of casamet%iveg2(nland) /= icewater
     
  enddo  ! end of nland

  do nland=1, mp
     if (casamet%iveg2(nland) /= icewater) then
        casapool%dClitterdt(nland,:) = casaflux%fluxCtolitter(nland,:) - &
             casaflux%klitter_tot(nland,:) * casapool%clitter(nland,:) 

        casapool%dCsoildt(nland,:)   = casaflux%fluxCtosoil(nland,:)   - &
             casaflux%ksoil(nland,:)   * casapool%csoil(nland,:)
        casaflux%Crsoil(nland)       = casaflux%fluxCtoCO2(nland)
        casaflux%cnep(nland)         = casaflux%cnpp(nland) - casaflux%Crsoil(nland)

        ! Nitrogen
        IF (icycle > 1) THEN
           casapool%dNlitterdt(nland,:) =  casaflux%fluxNtolitter(nland,:)  &
                - casaflux%klitter(nland,:) &
                * max(0.0_r_2,casapool%Nlitter(nland,:))
           
           casapool%dNsoildt(nland,:) = casaflux%FluxNtosoil(nland,:) &
                - casaflux%ksoil(nland,:) * casapool%Nsoil(nland,:)
           casapool%dNsoilmindt(nland)= casaflux%Nsnet(nland)&
                + casaflux%Nmindep(nland) + casaflux%Nminfix(nland)   &
                - casaflux%Nminloss(nland)   &
                - casaflux%Nminleach(nland)   &
                - casaflux%Nupland(nland)
        ENDIF ! end icycle>1

        ! Phosphorus
        IF (icycle >2) THEN
           fluxptase(nland) =  casabiome%prodptase( veg%iveg(nland) ) * deltcasa    &
                * max( 0.0_r_2, ( casapool%Psoil(nland,2)                   &
                * casaflux%ksoil(nland,2)                &
                + casapool%Psoil(nland,3)                &
                * casaflux%ksoil(nland,3) )              &
                )                                                 &
                * max( 0.0_r_2, ( casabiome%costNPup( veg%iveg(nland) )    &
                - 15.0_r_2 )                                 &
                )                                                 &
                / ( max( 0.0_r_2, ( casabiome%costNPup( veg%iveg(nland) )  &
                - 15.0_r_2 )                               &
                ) + 150.0_r_2                                       &
                )
           xdplabsorb(nland) = 1.0_r_2+ casaflux%Psorbmax(nland)*casaflux%kmlabp(nland) &
                /((casaflux%kmlabp(nland)+casapool%Psoillab(nland))**2)

           !  write(*,*) 'xdplabsorb:',  xdplabsorb(nland), casaflux%Psorbmax(nland), casaflux%kmlabp(nland), &
           !      casaflux%kmlabp(nland), casapool%Psoillab(nland)
           casapool%dPlitterdt(nland,:) = casaflux%fluxPtolitter(nland,:)  &
                - casaflux%klitter(nland,:)                 &
                * max(zero,casapool%Plitter(nland,:))

           casapool%dPsoildt(nland,1) = casaflux%FluxPtosoil(nland,1)                        &
                - casaflux%ksoil(nland,1) * casapool%Psoil(nland,1)

           casapool%dPsoildt(nland,2) = casaflux%FluxPtosoil(nland,2)                        &
                - casaflux%ksoil(nland,2) * casapool%Psoil(nland,2)    &
                - fluxptase(nland) * casaflux%ksoil(nland,2)*casapool%Psoil(nland,2) &
                /(casaflux%ksoil(nland,2)*casapool%Psoil(nland,2)+casaflux%ksoil(nland,3)*casapool%Psoil(nland,3))

           casapool%dPsoildt(nland,3) = casaflux%FluxPtosoil(nland,3)                        &
                - casaflux%ksoil(nland,3) * casapool%Psoil(nland,3)    &
                - fluxptase(nland) * casaflux%ksoil(nland,3)*casapool%Psoil(nland,3) &
                /(casaflux%ksoil(nland,2)*casapool%Psoil(nland,2)+casaflux%ksoil(nland,3)*casapool%Psoil(nland,3))

           casapool%dPsoillabdt(nland)= casaflux%Psnet(nland) + fluxptase(nland)         &
                + casaflux%Pdep(nland) + casaflux%Pwea(nland)      &
                - casaflux%Pleach(nland)-casaflux%pupland(nland)   &
                - casaflux%kpsorb(nland)*casapool%Psoilsorb(nland) &
                + casaflux%kpocc(nland) * casapool%Psoilocc(nland)

           ! here the dPsoillabdt =(dPsoillabdt+dPsoilsorbdt)
           ! dPsoilsorbdt  = xdplabsorb

           !     write(*,*) 'dPsoillabdt: ' ,  casapool%dPsoillabdt,  xdplabsorb(nland)
           casapool%dPsoillabdt(nland)  = casapool%dPsoillabdt(nland)/xdplabsorb(nland)
           !    write(*,*) 'dPsoillabdt: ' ,  casapool%dPsoillabdt
           casapool%dPsoilsorbdt(nland) = 0.0_r_2
           
           casapool%dPsoiloccdt(nland)  = casaflux%kpsorb(nland)* casapool%Psoilsorb(nland) &
                - casaflux%kpocc(nland) * casapool%Psoilocc(nland)
           ! P loss to non-ravailable P pools
           !      casaflux%Ploss(nland)        = casaflux%kpocc(nland) * casapool%Psoilocc(nland)
           
           !      casaflux%Ploss(nland)       = casaflux%fPleach(nland) &
           !                                 * max(zero,casapool%Psoillab(nland))
           casaflux%Ploss(nland)       = 0.0_r_2
        ENDIF ! end icycle>2
        
     endif ! end of casamet%iveg2(nland) /= icewater
     
  enddo ! end nland

END SUBROUTINE casa_delsoil


SUBROUTINE avgsoil(veg, soil, casamet)
  
  ! Get avg soil moisture, avg soil temperature
  ! need to estimate the land cell mean soil temperature and moisture weighted by the area fraction
  ! of each tile within the land cell

  implicit none
  
  type(veg_parameter_type),  intent(inout) :: veg     ! vegetation parameters
  type(soil_parameter_type), intent(inout) :: soil    ! soil parameters
  type(casa_met),            intent(inout) :: casamet

  integer :: ns, nland

  casamet%tsoilavg = 0.0_r_2
  casamet%moistavg = 0.0_r_2
  casamet%btran    = 0.0_r_2

  do ns = 1, ms
     do nland=1,mp
        casamet%tsoilavg(nland)  = casamet%tsoilavg(nland) + veg%froot(nland,ns) * casamet%tsoil(nland,ns)

        if (trim(cable_user%SMRF_NAME)=='Trudinger2016' .or. &
             trim(cable_user%SMRF_NAME)=='DAMM' ) then
           casamet%moistavg(nland)  = casamet%moistavg(nland) + &
                veg%froot(nland,ns) * casamet%moist(nland,ns)
        else
           casamet%moistavg(nland)  = casamet%moistavg(nland) + &
                real(veg%froot(nland,ns),r_2) * min(real(soil%sfc(nland),r_2), casamet%moist(nland,ns))
        endif

        ! Ticket#121
        ! casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
        !         * (min(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland)) &
        !         /(soil%sfc(nland)-soil%swilt(nland))
        casamet%btran(nland) = casamet%btran(nland) + real(veg%froot(nland,ns),r_2)  &
             * ( max( min(real(soil%sfc(nland),r_2), casamet%moist(nland,ns)) - real(soil%swilt(nland),r_2), 0.0_r_2 ) ) &
             / real(soil%sfc(nland)-soil%swilt(nland),r_2)
     enddo ! nland=1,mp
  enddo ! ns=1,ms

END SUBROUTINE avgsoil


SUBROUTINE casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
! computing the reduction in litter and SOM decomposition
! when decomposition rate is N-limiting
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(INOUT) :: xkNlimiting
  TYPE (casa_pool),         INTENT(INOUT) :: casapool
  TYPE (casa_flux),         INTENT(INOUT) :: casaflux
  TYPE (casa_met),          INTENT(INOUT) :: casamet
  TYPE (casa_biome),        INTENT(INOUT) :: casabiome
!
  TYPE (veg_parameter_type),   INTENT(IN) :: veg  ! vegetation parameters

  ! local variables
  INTEGER :: j,k,kk,nland
  REAL(r_2), DIMENSION(mp)         :: xFluxNlittermin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilmin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilimm
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilminnet
! A maximum Clitter set to avoid positive feedback for litter accumulation
! when N mode is activated. (Q.Zhang 23/05/2011)
!  real(r_2), dimension(17)         :: xClitter
!  data xClitter/100.0,100.0,100.0,100.0,50.0,150.0,150.0,100.0,&
!                150.0,150.0,100.0, 20.0,20.0, 20.0, 20.0, 20.0,20.0/

  xkNlimiting  = 1.0_r_2
!  set N mineral N fluxes to zero
  xFluxNlittermin(:)  = 0.0_r_2
  xFluxNsoilmin(:)    = 0.0_r_2
  xFluxNsoilimm(:)    = 0.0_r_2  !negative for microbial upatek and postive for release of mineral N
  xFluxNsoilminnet(:) = 0.0_r_2
!   PRINT *, 'within casa_xkN'

!  calculate gross mineralisation
  DO nland=1,mp
  IF (casamet%iveg2(nland)/=icewater) THEN

    ! calculate C:N ratio of newly formed SOM as function of soil mineral N pool
    IF (casapool%Nsoilmin(nland) < 2.0) THEN
      casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmin(nland,:)  &
                                       + (casapool%ratioNCsoilmax(nland,:) &
                                         -casapool%ratioNCsoilmin(nland,:)) &
                                       * max(0.0_r_2,casapool%Nsoilmin(nland)) / 2.0
    ELSE
      casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmax(nland,:)
    ENDIF

    DO j=1,mlitter
      xFluxNlittermin(nland) = xFluxNlittermin(nland) &
                         + casaflux%klitter(nland,j) * casapool%Nlitter(nland,j)
    ENDDO
    DO k=1,msoil
      xFluxNsoilmin(nland)   = xFluxNsoilmin(nland)   &
                         + casaflux%ksoil(nland,k)   * casapool%Nsoil(nland,k)
    ENDDO

    ! calculate N immobilisation from L to S and S to S
    DO kk=1,msoil
      DO j=1,mlitter    ! immobilisation from litter to soil
        xFluxNsoilimm(nland) = xFluxNsoilimm(nland) &
             - casaflux%fromLtoS(nland,kk,j) * casaflux%klitter(nland,j) &
             * casapool%Clitter(nland,j) * casapool%ratioNCsoilnew(nland,kk)
      ENDDO
      DO k=1,msoil      ! immobilisation from soil to soil
        IF(k.ne.kk) THEN
          xFluxNsoilimm(nland) = xFluxNsoilimm(nland) &
               - casaflux%fromStoS(nland,kk,k) * casaflux%ksoil(nland,k) &
               * casapool%Csoil(nland,k) * casapool%ratioNCsoilnew(nland,kk)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  ENDDO

  ! now check if there is sufficient mineral N
  xFluxNsoilminnet(:) = xFluxNlittermin(:) + xFluxNsoilmin(:) + xFluxNsoilimm(:)
!   PRINT *, 'casamet%iveg2 = ', casamet%iveg2
!   PRINT *, 'deltpool = ',deltpool
!   PRINT *, 'xFluxNsoilminnet = ', xFluxNsoilminnet
! WHERE(casamet%iveg2(:)/=icewater)
!    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0_r_2)
!      xkNlimiting(:) =1.0_r_2
!    ELSEWHERE
!      xkNlimiting(:) =max(0.0_r_2, - (casapool%Nsoilmin(:)-2.0)/(deltpool*xFluxNsoilminnet(:)))
!      xkNlimiting(:) =MIN(1.0_r_2,xkNlimiting(:))
!    ENDWHERE
! ENDWHERE

! Q.Zhang 23/05/2011 test code according to YPW
  WHERE(casamet%iveg2(:)/=icewater)
    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0_r_2)) > 0.0_r_2 &
          .OR. xFluxNsoilminnet(:) .ge. 0.0_r_2)
      xkNlimiting(:) =1.0_r_2
    ELSEWHERE
      xkNlimiting(:) =MAX(0.0_r_2, - (casapool%Nsoilmin(:)-0.5_r_2) &
                                /(deltpool*xFluxNsoilminnet(:)))
      xkNlimiting(:) =MIN(1.0_r_2,xkNlimiting(:))
    ENDWHERE
! Q.Zhang 23/05/2011 test
! If pool size larger than xClitter, turnover rate will not constrained by Nsoilmin.
!    where(casapool%clitter(:,1) > xClitter(veg%iveg(:)))
!     xkNlimiting(:) = 1.0_r_2
!    end where
! end (Q.Zhang 23/05/2011)
    where(sum(casapool%clitter,2) > casabiome%maxfinelitter(veg%iveg(:)) + casabiome%maxcwd(veg%iveg(:)))
     xkNlimiting(:) = 1.0_r_2
    end where
  ENDWHERE


END SUBROUTINE casa_xkN


SUBROUTINE casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
! (1) compute (1)N uptake by plants;
! (2) allocation of uptaken N to plants
!
  IMPLICIT NONE

  TYPE(veg_parameter_type), INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE(casa_biome),         INTENT(INOUT) :: casabiome
  TYPE(casa_pool),          INTENT(INOUT) :: casapool
  TYPE(casa_flux),          INTENT(INOUT) :: casaflux
  TYPE(casa_met),           INTENT(INOUT) :: casamet
  REAL(r_2), DIMENSION(mp), INTENT(IN)    :: xkNlimiting

  ! local variables
  INTEGER :: np
  REAL(r_2), DIMENSION(mp,mplant) :: Nreqmax, Nreqmin, NtransPtoP, xnuptake
  REAL(r_2), DIMENSION(mp)        :: totNreqmax, totNreqmin
  REAL(r_2), DIMENSION(mp)        :: xnCnpp

  Nreqmin(:,:)       = 0.0_r_2
  Nreqmax(:,:)       = 0.0_r_2
  NtransPtoP(:,:)    = 0.0_r_2
  totNreqmax = 0.0_r_2
  totNreqmin = 0.0_r_2

  casaflux%Nminuptake(:)     = 0.0_r_2
  casaflux%fracNalloc(:,:)   = 0.0_r_2
  xnCnpp = max(0.0_r_2,casaflux%Cnpp)
  call casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)

  DO np=1,mp
  IF(casamet%iveg2(np)/=icewater) THEN
    totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
    totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)

    xnuptake(np,leaf) = Nreqmin(np,leaf) + xkNlimiting(np)* (Nreqmax(np,leaf)-Nreqmin(np,leaf))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))
    xnuptake(np,wood) = Nreqmin(np,wood) + xkNlimiting(np)* (Nreqmax(np,wood)-Nreqmin(np,wood))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))
    xnuptake(np,froot) = Nreqmin(np,froot) + xkNlimiting(np)* (Nreqmax(np,froot)-Nreqmin(np,froot))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))

    casaflux%Nminuptake(np) = xnuptake(np,leaf) + xnuptake(np,wood) + xnuptake(np,froot)+1.0e-10_r_2
    casaflux%fracNalloc(np,leaf)  = xnuptake(np,leaf)/casaflux%Nminuptake(np)
    casaflux%fracNalloc(np,wood)  = xnuptake(np,wood)/casaflux%Nminuptake(np)
    casaflux%fracNalloc(np,froot) = xnuptake(np,froot)/casaflux%Nminuptake(np)
  ENDIF
  ENDDO


  casaflux%Nupland = casaflux%Nminuptake

END SUBROUTINE casa_nuptake


SUBROUTINE casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                         casabiome,casapool,casaflux,casamet)
!
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xnCnpp
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Nreqmax, Nreqmin
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: NtransPtoP
  TYPE (veg_parameter_type),             INTENT(INOUT) :: veg
  TYPE (casa_biome),                     INTENT(INOUT) :: casabiome
  TYPE (casa_pool),                      INTENT(INOUT) :: casapool
  TYPE (casa_flux),                      INTENT(INOUT) :: casaflux
  TYPE (casa_met),                       INTENT(INOUT) :: casamet

  ! local variable
  INTEGER :: np
  REAL(r_2), DIMENSION(mp,mplant)     :: ncplantmax

  Nreqmin(:,:)    = 0.0_r_2
  Nreqmax(:,:)    = 0.0_r_2
  NtransPtoP(:,:) = 0.0_r_2

  DO np=1,mp
  IF(casamet%iveg2(np)/=icewater) THEN
    if(casapool%Nsoilmin(np)<2.0) then
       ncplantmax(np,leaf) =casabiome%ratioNCplantmin(veg%iveg(np),leaf)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),leaf)-casabiome%ratioNCplantmin(veg%iveg(np),leaf)) &
                             * min(1.0_r_2,max(0.0_r_2,2.0_r_2**(0.5_r_2*casapool%Nsoilmin(np))-1.0_r_2))
       ncplantmax(np,wood) =casabiome%ratioNCplantmin(veg%iveg(np),wood)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),wood)-casabiome%ratioNCplantmin(veg%iveg(np),wood)) &
                             * min(1.0_r_2,max(0.0_r_2,2.0_r_2**(0.5_r_2*casapool%Nsoilmin(np))-1.0_r_2))
       ncplantmax(np,froot) =casabiome%ratioNCplantmin(veg%iveg(np),froot)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),froot)-casabiome%ratioNCplantmin(veg%iveg(np),froot)) &
                             * min(1.0_r_2,max(0.0_r_2,2.0_r_2**(0.5_r_2*casapool%Nsoilmin(np))-1.0_r_2))
    else
      ncplantmax(np,leaf)  = casabiome%ratioNCplantmax(veg%iveg(np),leaf)
      ncplantmax(np,wood)  = casabiome%ratioNCplantmax(veg%iveg(np),wood)
      ncplantmax(np,froot) = casabiome%ratioNCplantmax(veg%iveg(np),froot)
    endif

    Nreqmax(np,leaf)  = xnCnpp(np)* casaflux%fracCalloc(np,leaf) *ncplantmax(np,leaf)
    Nreqmax(np,wood)  = xnCnpp(np)* casaflux%fracCalloc(np,wood) *ncplantmax(np,wood)
    Nreqmax(np,froot) = xnCnpp(np)* casaflux%fracCalloc(np,froot)*ncplantmax(np,froot)

    Nreqmin(np,leaf) =  xnCnpp(np)* casaflux%fracCalloc(np,leaf) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),leaf)
    Nreqmin(np,wood) =  xnCnpp(np)* casaflux%fracCalloc(np,wood) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),wood)
    Nreqmin(np,froot) =  xnCnpp(np)* casaflux%fracCalloc(np,froot) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),froot)

    NtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Nplant(np,leaf) &
                       * (1.0_r_2-casabiome%ftransNPtoL(veg%iveg(np),leaf))
    NtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Nplant(np,wood) &
                       * (1.0_r_2-casabiome%ftransNPtoL(veg%iveg(np),wood))
    NtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Nplant(np,froot) &
                       * (1.0_r_2-casabiome%ftransNPtoL(veg%iveg(np),froot))

    Nreqmax(np,leaf)  = max(0.0_r_2,Nreqmax(np,leaf) - NtransPtoP(np,leaf))
    Nreqmax(np,wood)  = max(0.0_r_2,Nreqmax(np,wood) - NtransPtoP(np,wood))
    Nreqmax(np,froot) = max(0.0_r_2,Nreqmax(np,froot) - NtransPtoP(np,froot))
    Nreqmin(np,leaf)  = max(0.0_r_2,Nreqmin(np,leaf) - NtransPtoP(np,leaf))
    Nreqmin(np,wood)  = max(0.0_r_2,Nreqmin(np,wood) - NtransPtoP(np,wood))
    Nreqmin(np,froot) = max(0.0_r_2,Nreqmin(np,froot) - NtransPtoP(np,froot))

    if(casapool%nplant(np,leaf)/(casapool%cplant(np,leaf)+1.0e-10_r_2)>casabiome%ratioNCplantmax(veg%iveg(np),leaf)) then
       Nreqmax(np,leaf) = 0.0_r_2
       Nreqmin(np,leaf) =0.0_r_2
    endif
    if(casapool%nplant(np,wood)/(casapool%cplant(np,wood)+1.0e-10_r_2)>casabiome%ratioNCplantmax(veg%iveg(np),wood)) then
       Nreqmax(np,wood) = 0.0_r_2
       Nreqmin(np,wood) =0.0_r_2
    endif
    if(casapool%nplant(np,froot)/(casapool%cplant(np,froot)+1.0e-10_r_2)>casabiome%ratioNCplantmax(veg%iveg(np),froot)) then
       Nreqmax(np,froot) = 0.0_r_2
       Nreqmin(np,froot) =0.0_r_2
    endif

  ENDIF
  ENDDO

END SUBROUTINE casa_Nrequire


SUBROUTINE casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
! (1) compute  P uptake by plants;
! (2) allocation of uptaken P to plants
!
  IMPLICIT NONE
  
  TYPE(veg_parameter_type), INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE(casa_biome),         INTENT(INOUT) :: casabiome
  TYPE(casa_pool),          INTENT(INOUT) :: casapool
  TYPE(casa_flux),          INTENT(INOUT) :: casaflux
  TYPE(casa_met),           INTENT(INOUT) :: casamet
  REAL(r_2), DIMENSION(mp), INTENT(IN)    :: xkNlimiting

  ! local variables
  REAL(r_2), DIMENSION(mp,mplant) :: Preqmax, Preqmin, PtransPtoP, xPuptake
  REAL(r_2), DIMENSION(mp)        :: totPreqmax, totPreqmin
  REAL(r_2), DIMENSION(mp)        :: xpCnpp

  Preqmin(:,:)             = 0.0_r_2
  Preqmax(:,:)             = 0.0_r_2
  PtransPtoP(:,:)          = 0.0_r_2
  casaflux%Plabuptake(:)   = 0.0_r_2
  casaflux%fracPalloc(:,:) = 0.0_r_2
  totPreqmax               = 0.0_r_2
  totPreqmin               = 0.0_r_2

  xpCnpp = max(0.0_r_2,casaflux%cnpp)
  call casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
  WHERE(casamet%iveg2/=icewater)
    totPreqmax(:) = Preqmax(:,leaf)+Preqmax(:,wood)+Preqmax(:,froot)
    totPreqmin(:) = Preqmin(:,leaf)+Preqmin(:,wood)+Preqmin(:,froot)

    xpuptake(:,leaf) = Preqmin(:,leaf) + xkNlimiting(:)* (Preqmax(:,leaf)-Preqmin(:,leaf))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))
    xpuptake(:,wood) = Preqmin(:,wood) + xkNlimiting(:)* (Preqmax(:,wood)-Preqmin(:,wood))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))
    xpuptake(:,froot) = Preqmin(:,froot) + xkNlimiting(:)* (Preqmax(:,froot)-Preqmin(:,froot))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))

    casaflux%Plabuptake(:) = xpuptake(:,leaf) + xpuptake(:,wood) + xpuptake(:,froot)+1.0e-10_r_2
    casaflux%fracPalloc(:,leaf)  = xpuptake(:,leaf)/casaflux%Plabuptake(:)
    casaflux%fracPalloc(:,wood)  = xpuptake(:,wood)/casaflux%Plabuptake(:)
    casaflux%fracPalloc(:,froot) = xpuptake(:,froot)/casaflux%Plabuptake(:)

  ENDWHERE

  casaflux%Pupland = casaflux%Plabuptake

!  ! only used in spinning up the model
!  DO  np=1,mp
!    casaflux%Plabuptake(np) = TotPreqmax(np)
!    casaflux%Pupland(np)    = TotPreqmax(np)
!    casaflux%Pwea(np)       = TotPreqmax(np)
!  ENDDO

END SUBROUTINE casa_puptake


SUBROUTINE casa_Prequire(xpCnpp, Preqmin, Preqmax, PtransPtoP, veg, &
     casabiome, casapool, casaflux, casamet)
  
  IMPLICIT NONE
  
  REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xpCnpp
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Preqmax, Preqmin
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: PtransPtoP
  TYPE(veg_parameter_type),        INTENT(INOUT) :: veg
  TYPE(casa_biome),                INTENT(INOUT) :: casabiome
  TYPE(casa_pool),                 INTENT(INOUT) :: casapool
  TYPE(casa_flux),                 INTENT(INOUT) :: casaflux
  TYPE(casa_met),                  INTENT(INOUT) :: casamet

  ! local variables
  INTEGER :: np

  Preqmin(:,:)    = 0.0_r_2
  Preqmax(:,:)    = 0.0_r_2
  PtransPtoP(:,:) = 0.0_r_2
  
  do np=1,mp
     if (casamet%iveg2(np)/=icewater) then
        Preqmax(np,leaf) = xpCnpp(np)* casaflux%fracCalloc(np,leaf) &
             * (casapool%Nplant(np,leaf)/(casapool%Cplant(np,leaf)+1.0e-10_r_2)) &
             / casabiome%ratioNPplantmin(veg%iveg(np),leaf)
        Preqmax(np,wood) = xpCnpp(np)* casaflux%fracCalloc(np,wood) &
             * (casapool%Nplant(np,wood)/(casapool%Cplant(np,wood)+1.0e-10_r_2)) &
             / casabiome%ratioNPplantmin(veg%iveg(np),wood)
        Preqmax(np,froot) = xpCnpp(np)* casaflux%fracCalloc(np,froot) &
             * (casapool%Nplant(np,froot)/(casapool%Cplant(np,froot)+1.0e-10_r_2)) &
             / casabiome%ratioNPplantmin(veg%iveg(np),froot)

        Preqmin(np,leaf) = xpCnpp(np) * casaflux%fracCalloc(np,leaf) &
             * (casapool%Nplant(np,leaf)/(casapool%Cplant(np,leaf)+1.0e-10_r_2)) &
             / casabiome%ratioNPplantmax(veg%iveg(np),leaf)
        Preqmin(np,wood) = xpCnpp(np) * casaflux%fracCalloc(np,wood) &
             * (casapool%Nplant(np,wood)/(casapool%Cplant(np,wood)+1.0e-10_r_2)) &
             / casabiome%ratioNPplantmax(veg%iveg(np),wood)
        Preqmin(np,froot) = xpCnpp(np) * casaflux%fracCalloc(np,froot) &
             * (casapool%Nplant(np,froot)/(casapool%Cplant(np,froot)+1.0e-10_r_2)) &
             / casabiome%ratioNPplantmax(veg%iveg(np),froot)
        
        PtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Pplant(np,leaf) &
             * (1.0_r_2-casabiome%ftransPPtoL(veg%iveg(np),leaf))
        PtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Pplant(np,wood) &
             * (1.0_r_2-casabiome%ftransPPtoL(veg%iveg(np),wood))
        PtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Pplant(np,froot) &
             * (1.0_r_2-casabiome%ftransPPtoL(veg%iveg(np),froot))
        
        Preqmax(np,leaf)  = max(0.0_r_2,Preqmax(np,leaf) - PtransPtoP(np,leaf))
        Preqmax(np,wood)  = max(0.0_r_2,Preqmax(np,wood) - PtransPtoP(np,wood))
        Preqmax(np,froot) = max(0.0_r_2,Preqmax(np,froot) - PtransPtoP(np,froot))
        
        Preqmin(np,leaf)  = max(0.0_r_2,Preqmin(np,leaf) - PtransPtoP(np,leaf))
        Preqmin(np,wood)  = max(0.0_r_2,Preqmin(np,wood) - PtransPtoP(np,wood))
        Preqmin(np,froot) = max(0.0_r_2,Preqmin(np,froot) - PtransPtoP(np,froot))

        if (casapool%pplant(np,leaf)/(casapool%nplant(np,leaf)+1.0e-10_r_2) > &
             1.0_r_2/casabiome%ratioNPplantmin(veg%iveg(np),leaf)) then
           Preqmax(np,leaf) = 0.0_r_2
           Preqmin(np,leaf) = 0.0_r_2
        endif
        if (casapool%pplant(np,wood)/(casapool%nplant(np,wood)+1.0e-10_r_2) > &
             1.0_r_2/casabiome%ratioNPplantmin(veg%iveg(np),wood)) then
           Preqmax(np,wood) = 0.0_r_2
           Preqmin(np,wood) = 0.0_r_2
        endif
        if(casapool%pplant(np,froot)/(casapool%nplant(np,froot)+1.0e-10_r_2) > &
             1.0_r_2/casabiome%ratioNPplantmin(veg%iveg(np),froot)) then
           Preqmax(np,froot) = 0.0_r_2
           Preqmin(np,froot) = 0.0_r_2
        endif

     endif
     
  enddo

END SUBROUTINE casa_Prequire


SUBROUTINE casa_cnpcycle(veg, casabiome, casapool, casaflux, casamet, LALLOC)

  ! update all pool sizes

  implicit none
  
  type(veg_parameter_type), intent(inout) :: veg       ! vegetation parameters
  type(casa_biome),         intent(inout) :: casabiome
  type(casa_pool),          intent(inout) :: casapool
  type(casa_flux),          intent(inout) :: casaflux
  type(casa_met),           intent(inout) :: casamet
  integer,                  intent(in)    :: lalloc

  ! local variables
  integer :: i, j, k, np

 
  do np=1, mp

     if (casamet%iveg2(np) == icewater) then
        casamet%glai(np) = 0.0_r_2
     else
        casapool%cplant(np,:) = casapool%cplant(np,:) + casapool%dcplantdt(np,:) * deltpool
        casapool%clabile(np)  = casapool%clabile(np)  + casapool%dclabiledt(np)  * deltpool

        if (casapool%cplant(np,leaf) > 0.0_r_2) then
           if (icycle >1) casapool%Nplant(np,:) = casapool%Nplant(np,:) + casapool%dNplantdt(np,:)*deltpool
           if (icycle >2) casapool%Pplant(np,:) = casapool%Pplant(np,:) + casapool%dPplantdt(np,:)*deltpool
        endif

        ! avoid high ratios of N to P in plant material
        casapool%Nplant(np,3) = min( casapool%Nplant(np,3), &
             casabiome%ratioNCplantmax(veg%iveg(np),froot) * casapool%cplant(np,3) )
        casamet%glai(np)      = max( casabiome%glaimin(veg%iveg(np)), &
             casabiome%sla(veg%iveg(np)) * casapool%cplant(np,leaf) )
        ! vh !
        !IF (LALLOC.ne.3) THEN
        casamet%glai(np) = min(casabiome%glaimax(veg%iveg(np)), casamet%glai(np))
        !ENDIF
        casapool%clitter(np,:) = casapool%clitter(np,:) + casapool%dClitterdt(np,:) * deltpool
        casapool%csoil(np,:)   = casapool%csoil(np,:)   + casapool%dCsoildt(np,:)   * deltpool
        
        IF (icycle >1) THEN
           casapool%Nlitter(np,:) = casapool%Nlitter(np,:) + casapool%dNlitterdt(np,:) * deltpool
           casapool%Nsoil(np,:)   = casapool%Nsoil(np,:)   + casapool%dNsoildt(np,:)   * deltpool
           ! vh ! put lower bound of 1.e-3 to prevent Nsoilmin from going negative
           ! Ticket #108
           casapool%Nsoilmin(np)  = max(casapool%Nsoilmin(np) + casapool%dNsoilmindt(np) * deltpool,1.e-3_r_2)
        ENDIF
        
        IF (icycle >2) THEN
           casapool%Plitter(np,:) = casapool%Plitter(np,:) + casapool%dPlitterdt(np,:) * deltpool
           casapool%Psoil(np,:)   = casapool%Psoil(np,:)   + casapool%dPsoildt(np,:)   * deltpool
           ! vh ! put lower bound of 1.e-3 to prevent Psoillab from going negative
           casapool%Psoillab(np)  = max(casapool%Psoillab(np) + casapool%dPsoillabdt(np) * deltpool, 1.e-3_r_2)
           casapool%Psoilsorb(np) = casaflux%Psorbmax(np)*casapool%Psoillab(np) / (casaflux%kmlabp(np)+casapool%Psoillab(np))
           !      casapool%Psoilsorb(np) = casapool%Psoilsorb(np)  &
           !                             + casapool%dPsoilsorbdt(np) * deltpool
           casapool%Psoilocc(np)  = casapool%Psoilocc(np) + casapool%dPsoiloccdt(np) * deltpool
        ENDIF

        do i=1, mplant
           if (casapool%cplant(np,i) < 0.0_r_2)  then
              write(57,*)  'Cpool: np,ivt', np, casamet%lat(np), casamet%lon(np), &
                   casamet%iveg2(np), casapool%cplant(np,:)
              write(*,*)  'Cpool: np,ivt', np, casamet%lat(np), casamet%lon(np), &
                   casamet%iveg2(np), casapool%cplant(np,:)
              call casa_poolzero(np, 1, casapool)
              casapool%cplant(np,i) = max(0.0_r_2, casapool%cplant(np,i))
           endif
        enddo
        
        IF (icycle >1) THEN
           DO i=1, mplant
              IF(casapool%nplant(np,i) < 0.0_r_2) THEN
                 WRITE(57,*) 'Npool:', 'np,ivt,ipool',np,casamet%iveg2(np),casapool%nplant(np,:)
                 call casa_poolzero(np,2,casapool)
                 casapool%nplant(np,i) = max(0.0_r_2, casapool%nplant(np,i))
              ENDIF
           ENDDO
        ENDIF ! end of "icycle >1"

        do j=1, mlitter
           if (casapool%clitter(np,j) < 0.0_r_2) then
              write(57,*) 'Clitter: np,ivt2',np, casamet%iveg2(np), casapool%clitter(np,:)
              write(*,*) 'Clitter: np,ivt2',np, casamet%iveg2(np), casapool%clitter(np,:)
              call casa_poolzero(np, 3, casapool)
              casapool%clitter(np,j) = max(0.0_r_2, casapool%clitter(np,j))
           endif
        enddo

        do k=1, msoil
           if (casapool%csoil(np,k) < 0.0_r_2) then
              write(57,*) 'Csoil: np,ivt2', np, casamet%iveg2(np), casapool%csoil(np,:)
              write(*,*) 'Csoil: np,ivt2', np, casamet%iveg2(np), casapool%csoil(np,:)
              call casa_poolzero(np, 5, casapool)
              casapool%csoil(np,k) = max(0.0_r_2, casapool%csoil(np,k))
           endif
        enddo
        
        !  check if any pool size, and terminate model run if any pool size is negative!!
        IF (icycle >1) THEN
           DO j=1,mlitter
              IF(casapool%nlitter(np,j) < 0.0_r_2)  THEN
                 WRITE(57,*)  'Nlitter: np,ivt2',np,casamet%iveg2(np),casapool%Nlitter(np,:)
                 call casa_poolzero(np,4,casapool)
                 casapool%nlitter(np,j) = max(0.0_r_2, casapool%nlitter(np,j))
              ENDIF
           ENDDO
           DO k=1,msoil
              IF(casapool%nsoil(np,k) < 0.0_r_2) THEN
                 WRITE(57,*)  'Nsoil: np,ivt2',np,casamet%iveg2(np),casapool%nsoil(np,:)
                 call casa_poolzero(np,6,casapool)
                 casapool%nsoil(np,k) = max(0.0_r_2, casapool%nsoil(np,k))
              ENDIF
           ENDDO
        ENDIF  !end of "icycle >1"
        
     endif ! == icewater
     
  enddo ! np=1,mp
  
END SUBROUTINE casa_cnpcycle


SUBROUTINE casa_poolzero(n,ipool,casapool)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ipool
  TYPE (casa_pool), INTENT(INOUT) :: casapool

  WRITE(57,*) ' WARNING: negative pools are reset to ZERO!!'
  SELECT CASE(ipool)
  CASE(1)
     WRITE(57,*) 'plant carbon pool size negative!!'
     WRITE(57,*) 'plant C pools: ', n,casapool%cplant(n,:)
  CASE(2)
     WRITE(57,*) 'plant nitrogen pool size negative!!'
     WRITE(57,*) 'plant C pools: ',n,casapool%cplant(n,:)
     WRITE(57,*) 'plant N pools: ',n,casapool%nplant(n,:)
  CASE(3)
     WRITE(57,*) 'litter carbon pool size negative!!'
     WRITE(57,*) 'litter C pools: ',n,casapool%clitter(n,:)
  CASE(4)
     WRITE(57,*) 'litter nitrogen pool size negative!!'
     WRITE(57,*) 'carbon pool: ',n,casapool%clitter(n,:)
     WRITE(57,*) 'nitrogen pools: ',n,casapool%nlitter(n,:)
  CASE(5)
     WRITE(57,*) 'soil carbon pool size negative!!'
     WRITE(57,*) 'soil C pools: ',n,casapool%csoil(n,:)
  CASE(6)
     WRITE(57,*) 'soil nitrogen pool size negative!!'
     WRITE(57,*) 'soil C pools: ', n,casapool%csoil(n,:)
     WRITE(57,*) 'soil N pools: ', n,casapool%nsoil(n,:)
  END SELECT

END SUBROUTINE casa_poolzero


SUBROUTINE casa_cnpbal(casapool,casaflux,casabal)

  IMPLICIT NONE

  TYPE(casa_pool),    INTENT(INOUT) :: casapool
  TYPE(casa_flux),    INTENT(INOUT) :: casaflux
  TYPE(casa_balance), INTENT(INOUT) :: casabal

  ! local variables
  REAL(r_2), DIMENSION(mp) :: cbalplant, nbalplant, pbalplant
  REAL(r_2), DIMENSION(mp) :: cbalsoil,  nbalsoil,  pbalsoil

  REAL(r_2) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

  cbalplant(:) = 0.0_r_2
  cbalsoil(:)  = 0.0_r_2
  nbalplant(:) = 0.0_r_2
  nbalsoil(:)  = 0.0_r_2
  pbalplant(:) = 0.0_r_2
  pbalsoil(:)  = 0.0_r_2

  casabal%cbalance(:)  = 0.0_r_2
  casabal%nbalance(:)  = 0.0_r_2
  casabal%pbalance(:)  = 0.0_r_2

  !C balance

  ! plant (change in stock = npp minus turnover of plant pools)
   Cbalplant(:)  = sum(casabal%cplantlast,2) -sum(casapool%cplant,2)            &
                 + casabal%Clabilelast(:)-casapool%clabile(:)                   &
                 +(casaflux%Cnpp(:) - SUM((casaflux%kplant_tot*casabal%cplantlast),2))*deltpool &
                 + casapool%dClabiledt(:)* deltpool

   !soil and litter (change in stock = (base plant turnover - heterotrophic resp) +
   ! plant turnover by fire (excluding fraction lost to atm) -
   ! litter loss to atmosphere by fire
   Cbalsoil(:)   = sum(casabal%clitterlast,2) - sum(casapool%clitter,2)         &
                 + sum(casabal%csoillast,2)   - sum(casapool%csoil,2)           &
                 +(SUM((casaflux%kplant*casabal%cplantlast),2)-casaflux%Crsoil(:))*deltpool &
   +(SUM((casaflux%kplant_fire*(1.0_r_2 - casaflux%kplant) * casabal%cplantlast),2) - &
   casaflux%fluxCtoCO2_plant_fire)*deltpool - casaflux%fluxCtoCO2_litter_fire*deltpool
   

   if (abs(Cbalsoil(1)).gt.0.1) then
      print*, 'Cbalsoil(1)', Cbalsoil(1)
      !print*, sum(casabal%clitterlast(1,:)) - sum(casapool%clitter(1,:))
      ! mass bal on litter
      tmp1 = SUM((casaflux%kplant(1,:)*casabal%cplantlast(1,:))) ! plant input
      tmp2 = SUM(casaflux%kplant_fire(1,:) * (1.0_r_2 - casaflux%kplant(1,:)) * casabal%cplantlast(1,:)) - &
           casaflux%fluxCtoCO2_plant_fire(1)
      tmp3 = SUM((casaflux%klitter(1,:)*casabal%clitterlast(1,:)))
      tmp4 = SUM((casaflux%klitter_fire(1,:)*(1.0_r_2 - casaflux%klitter(1,:)) * casabal%clitterlast(1,:)))
      tmp5 = SUM((casaflux%kplant_tot(1,:)*casabal%cplantlast(1,:)))
      tmp6 = casaflux%fluxCtoCO2_plant_fire(1)
      tmp7 = SUM((casaflux%klitter_tot(1,:)*casabal%clitterlast(1,:)))
      tmp8 = SUM((casaflux%kplant_fire(1,:)*(1.0_r_2 - casaflux%kplant(1,:)) * casabal%cplantlast(1,:)))

      print*, 'delta clitt1: ', sum(casabal%clitterlast(1,:)) - sum(casapool%clitter(1,:))
      print*, 'plant input to litter (base turnover): ', tmp1
      print*, 'plant input to litter (fire): ', tmp2
      print*, 'plant input to litter (base turnover + fire): ', tmp1 + tmp2
      print*, 'fluxctolitter: ', sum(casaflux%FluxCtolitter(1,:))
      print*, 'plant loss to atm (fire): ', tmp6
      print*, 'total plant turnover: ', tmp5
      print*, 'litter loss (base turnover): ', tmp3
      print*, 'litter loss (fire): ', tmp4
      print*, 'litter loss (klitter_tot): ',  tmp7
      print*, 'delta clitt2: ', tmp1 + tmp2 - tmp3 - tmp4
      print*, 'fluxCtolitter: ',  SUM(casaflux%fluxCtolitter(1,:))
      print*, 'frac plant fire flux to CO2: ',  casaflux%FluxFromPtoCO2(1,:)
      print*, 'frac plant fire flux to litter: ',  sum(casaflux%fromPtoL_fire(1,:,:),1)
      print*, 'plant loss (fire)', tmp8
  
     ! stop
     ! print*, sum(casabal%csoillast(1,:))   - sum(casapool%csoil(1,:))
     ! print*, SUM((casaflux%kplant(1,:)*casabal%cplantlast(1,:)))-casaflux%Crsoil(1)*deltpool
      !print*, SUM((casaflux%kplant_fire(1,:)*(1.0_r_2 - casaflux%kplant(1,:)) * casabal%cplantlast(1,:)))
     ! print*, casaflux%fluxCtoCO2_plant_fire(1),  casaflux%fluxCtoCO2_litter_fire(1)
      !stop
      endif
   
   
   casabal%cbalance(:) = Cbalplant(:) + Cbalsoil(:)

   ! do npt=1,mp
   !    IF(abs(casabal%cbalance(npt))>1e-10) THEN
   !       write(*,*) 'cbalance',  npt, Cbalplant(npt), Cbalsoil(npt)
   !       write(*,*) 'soil input', SUM((casaflux%kplant_tot(npt,:)*casabal%cplantlast(npt,:))), &
   !          casaflux%kplant_tot(npt,1), casaflux%kplant_tot(npt,3),casaflux%kplant(npt,1), casaflux%kplant(npt,3)
   !       write(*,*)  'soil efflux', casaflux%Crsoil(npt)
   !       write(*,*) 'dclitter', casapool%clitter(npt,:) -  casabal%clitterlast(npt,:)
   !       write(*,*) 'dcsoil', casapool%csoil(npt,:) -  casabal%csoillast(npt,:)
   !       write(*,*) 'cplant', casapool%cplant(npt,:)
   !       write(*,*) 'gpp, npp',casaflux%Cgpp(npt) , &
   !           casaflux%Cnpp(npt)
   !      write(*,*) 'dcplandt',  casapool%dcplantdt(npt,:), sum(casapool%dcplantdt(npt,:))
   !      write(*,*) 'rmplant, rgplant',  casaflux%crmplant(npt,:) , casaflux%crgplant(npt)
   !      write(*,*), 'dclabile',  casapool%dClabiledt(npt)* deltpool
   !       !STOP
   !    ENDIF
   ! ENDDO

   casapool%ctot_0 = sum(casabal%cplantlast,2)+sum(casabal%clitterlast,2) &
        + sum(casabal%csoillast,2)+ casabal%clabilelast
   casapool%ctot = sum(casapool%cplant,2)+sum(casapool%clitter,2) &
        + sum(casapool%csoil,2)+ casapool%clabile
   casabal%cplantlast  = casapool%cplant
   casabal%clabilelast = casapool%clabile
   casabal%clitterlast = casapool%clitter
   casabal%csoillast   = casapool%csoil
   casabal%sumcbal     = casabal%sumcbal + casabal%cbalance


   IF(icycle >1) THEN
      Nbalplant(:) = sum(casabal%nplantlast,2) -sum(casapool%nplant,2)                  &
                    +casaflux%Nminuptake(:) *deltpool
      Nbalsoil(:)  = -sum(casapool%nlitter,2)-sum(casapool%nsoil,2)                     &
                     -casapool%nsoilmin(:)+ casabal%nsoilminlast(:)                     &
                     + sum(casabal%nlitterlast,2)    + sum(casabal%nsoillast,2)         &
                     +(casaflux%Nmindep(:) + casaflux%Nminfix(:)- casaflux%Nminloss(:)                        &
                       -casaflux%Nminleach(:)-casaflux%Nupland(:)) * deltpool

      casabal%nbalance(:) = Nbalplant(:) + Nbalsoil(:)

      casabal%nplantlast  = casapool%nplant
      casabal%nlitterlast = casapool%nlitter
      casabal%nsoillast   = casapool%nsoil
      casabal%nsoilminlast= casapool%nsoilmin
      casabal%sumnbal     = casabal%sumnbal + casabal%nbalance

   ENDIF

   IF(icycle >2) THEN
      Pbalplant(:) = sum(casabal%Pplantlast,2) -sum(casapool%Pplant,2)                      &
                   + casaflux%Plabuptake(:) *deltpool
      Pbalsoil(:)  = -sum(casapool%Plitter,2)        - sum(casapool%Psoil,2)                      &
                     + sum(casabal%Plitterlast,2)    + sum(casabal%Psoillast,2)                   &
                   -casapool%psoillab(:)-casapool%psoilsorb(:)-casapool%psoilocc(:)               &
                   + casabal%psoillablast(:) + casabal%psoilsorblast(:) + casabal%psoilocclast(:) &
                   +(casaflux%Pdep(:) + casaflux%Pwea(:)                                          &
                     -casaflux%Pleach(:)-casaflux%Pupland(:)                                      &
                     -casaflux%Ploss(:)) * deltpool

      casabal%pbalance(:) = pbalplant(:) + pbalsoil(:)

      casabal%pplantlast   = casapool%pplant
      casabal%plitterlast  = casapool%plitter
      casabal%psoillast    = casapool%psoil
      casabal%psoillablast = casapool%psoillab
      casabal%psoilsorblast= casapool%psoilsorb
      casabal%psoilocclast = casapool%psoilocc
      casabal%sumpbal  = casabal%sumpbal + casabal%pbalance
   ENDIF

   !write(6999,"(100(f12.5,2x))"),  casabal%cbalance(:)
   !write(8999,"(100(f12.5,2x))"), - (casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)
   !write(7999,"(100(f12.5,2x))"),  casapool%ctot - casapool%ctot_0
   !write(9999,"(100(f12.5,2x))"), casapool%ctot - casapool%ctot_0 + &
   ! ( casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)

END SUBROUTINE casa_cnpbal


SUBROUTINE casa_ndummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Nplant(:,:) = casapool%Cplant(:,:) * casapool%ratioNCplant(:,:)

END SUBROUTINE casa_ndummy


SUBROUTINE casa_pdummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Pplant(:,:) = casapool%Nplant(:,:) / casapool%ratioNPplant(:,:)

END SUBROUTINE casa_pdummy


SUBROUTINE phenology(iday,veg,phen)
  IMPLICIT NONE
  INTEGER,              INTENT(IN)    :: iday
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (phen_variable), INTENT(INOUT) :: phen

  ! local variables (temprary)
  INTEGER :: np
  INTEGER, DIMENSION(mp)  :: days,days1to2, days2to3, days3to4, days4to1

!  PRINT *, 'Within SUBROUTINE phenology, mp = ', mp
  DO np=1,mp
    days1to2(np) = phen%doyphase(np,2) - phen%doyphase(np,1)
    days2to3(np) = phen%doyphase(np,3) - phen%doyphase(np,2)
    days3to4(np) = phen%doyphase(np,4) - phen%doyphase(np,3)
    days4to1(np) = phen%doyphase(np,1) - phen%doyphase(np,4)
    IF(days1to2(np) < 0) days1to2(np) = days1to2(np) +365
    IF(days2to3(np) < 0) days2to3(np) = days2to3(np) +365
    IF(days3to4(np) < 0) days3to4(np) = days3to4(np) +365
    IF(days4to1(np) < 0) days4to1(np) = days4to1(np) +365
  ENDDO
  ! compute leaf phenology
  DO np=1,mp
    SELECT CASE(phen%phase(np))
      CASE(0)
        days(np) = iday - phen%doyphase(np,4)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days4to1(np)) phen%phase(np) =1
      CASE(1)
        days(np) = iday - phen%doyphase(np,1)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days1to2(np)) phen%phase(np) =2
      CASE(2)
        days(np) = iday - phen%doyphase(np,2)
        IF(days(np) <0) days(np) = days(np) + 365
        IF(days(np) > days2to3(np)) phen%phase(np) =3
      CASE(3)
        days(np) = iday - phen%doyphase(np,3)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days3to4(np)) phen%phase(np) =0
    END SELECT
  ENDDO

  WHERE(veg%iveg==1 .or. veg%iveg ==2 )
       phen%phase = 2
  ENDWHERE

END SUBROUTINE phenology


REAL FUNCTION vcmax_np(nleaf, pleaf)
  
  implicit none

  real, intent(in) :: nleaf ! leaf N in g N m-2 leaf
  real, intent(in) :: pleaf ! leaf P in g P m-2 leaf

  ! Walker, A. P. et al.: The relationship of leaf photosynthetic traits - Vcmax and Jmax -
  !   to leaf nitrogen, leaf phosphorus, and specific leaf area: a meta-analysis and modeling study,
  !   Ecology and Evolution, 4, 3218-3235, 2014.
  vcmax_np = exp(3.946 + 0.921*log(nleaf) + 0.121*log(pleaf) + &
       0.282*log(pleaf)*log(nleaf)) * 1.0e-6 ! units of mol m-2 (leaf)

END FUNCTION vcmax_np


END MODULE casa_cnp_module
