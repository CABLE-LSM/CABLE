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
! ==============================================================================
!
! This module contains the following subroutines:
!   casa_xnp
!   casa_allocation
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
  USE cable_common_module, ONLY: cable_user,l_landuse ! Custom soil respiration: Ticket #42
  USE landuse_constant
  IMPLICIT NONE
  REAL(r_2), PARAMETER :: zero = 0.0_r_2
  REAL(r_2), PARAMETER :: one  = 1.0_r_2
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

    xnlimit  = 1.0
    xplimit  = 1.0
    xnplimit = 1.0
    casaflux%fracClabile(:) = 0.0

    SELECT CASE(icycle)
    CASE(2)
       WHERE(casamet%iveg2/=icewater)
          xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
          !xnlimit(:) = xncleaf(:)/(xncleaf(:)+casabiome%KminN(veg%iveg(:)))
          xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01)
          xplimit(:) = 1.0
          xnplimit(:) =MIN(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
       ENDWHERE
    CASE(3)
       WHERE(casamet%iveg2/=icewater)
          xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
          xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01)
          xpcleaf(:) = casapool%pplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
          !xplimit(:) = xpcleaf(:)/(xpcleaf(:)+casabiome%Kuplabp(veg%iveg(:)))
          xplimit(:) = xpcleaf(:)/(xpcleaf(:)+0.0006)
          !xnplimit(:) = min(1.0,casabiome%Kuptake(veg%iveg(:))*min(xnlimit(:),xplimit(:)))
          xnplimit(:) =MIN(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
       ENDWHERE
    END SELECT

    ! now check if soil nutrient supply can meet the plant uptake,
    ! otherwise reduce NPP
    xNuptake = 1.0
    xPuptake = 1.0

    IF(icycle >1) THEN
       Nreqmin(:,:)    = 0.0
       Nreqmax(:,:)    = 0.0
       NtransPtoP(:,:) = 0.0
       totNreqmax = 0.0
       totNreqmin = 0.0
       xNuptake   = 1.0

       xnCnpp = MAX(0.0,casaflux%Cnpp)
       CALL casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
            casabiome,casapool,casaflux,casamet)
       DO np=1,mp
          IF(casamet%iveg2(np)/=icewater) THEN
             totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
             totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)
             xNuptake(np)   = MAX(0.0,MIN(1.0,casapool%Nsoilmin(np) &
                  /(totNreqmin(np)*deltpool+1.0e-10)))
          ENDIF
       ENDDO
    ENDIF
    IF(icycle >2) THEN
       Preqmin(:,:)       = 0.0
       Preqmax(:,:)       = 0.0
       PtransPtoP(:,:)    = 0.0
       totPreqmax = 0.0
       totPreqmin = 0.0
       xPuptake   = 1.0
       xpCnpp = MAX(0.0,casaflux%Cnpp)
       CALL casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
            casabiome,casapool,casaflux,casamet)
       DO np=1,mp
          IF(casamet%iveg2(np)/=icewater) THEN
             totPreqmax(np) = Preqmax(np,leaf)+Preqmax(np,wood)+Preqmax(np,froot)
             totPreqmin(np) = Preqmin(np,leaf)+Preqmin(np,wood)+Preqmin(np,froot)
             xPuptake(np)   = MAX(0.0,MIN(1.0,casapool%psoillab(np) &
                  /(totPreqmin(np)*deltpool+1.0e-10)))
          ENDIF
       ENDDO
    ENDIF

    xnplimit(:)  = 1.0
    xNPuptake(:)     = MIN(xnuptake(:), xpuptake(:))
    DO np =1, mp
       IF(casamet%iveg2(np)/=icewater.AND.casaflux%cnpp(np) > 0.0.AND.xNPuptake(np) < 1.0) THEN
          casaflux%fracClabile(np) =MIN(1.0,MAX(0.0,(1.0- xNPuptake(np)))) * MAX(0.0,casaflux%cnpp(np))/(casaflux%cgpp(np) +1.0e-10)
          casaflux%cnpp(np)    = casaflux%cnpp(np) - casaflux%fracClabile(np) * casaflux%cgpp(np)
       ENDIF

    ENDDO

    !write(59,91)  xNuptake(1),casapool%Nsoilmin(1), totNreqmin(1)*deltpool
91  FORMAT(20(e12.4,2x))
    !  casaflux%cnpp(:) = xNPuptake(:) * xnplimit(:) * casaflux%cnpp(:)


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
    INTEGER :: npt,ns,is,iv
    REAL(r_2), DIMENSION(mp,mplant) :: fracCallocx
    REAL(r_2), DIMENSION(mp,mplant) :: delc
    REAL(r_2), DIMENSION(mp)        :: ctotal
    REAL(r_2), DIMENSION(mp)        :: xLalloc,xwsalloc,xTalloc
    REAL(r_2), DIMENSION(mp)        :: xWorNalloc,xNalloc,xWalloc
    REAL(r_2), DIMENSION(mp)        :: totfracCalloc
    REAL(r_2), DIMENSION(mp)        :: newLAI
    ! initlization
    casaflux%fracCalloc  = 0.0
    !casaflux%fracClabile = 0.0
    fracCallocx = 0.0
    newLAI = 0.0
    SELECT CASE (LALLOC)

    CASE(2)   !
       ! calculate the allocation coefficients
       CALL casa_wolf(veg,casabiome,casaflux,casapool,casamet)

    CASE(1)   ! dynamic allocation
       WHERE(casamet%iveg2/=icewater)
          xLalloc(:) = MIN(1.0,MAX(0.0,EXP(-0.5*casamet%glai(:))))   ! L limiting
          ! Pseudo-nutrient limitation calculation
          WHERE(casamet%tsoilavg > 0.0)
             xwsalloc(:) = MIN( MAX(casamet%moistavg(:)-soil%swilt(:),0.0) &
                  /(soil%sfc(:)-soil%swilt(:)), 1.0 )
          ELSE WHERE
             xwsalloc(:) = 0.01
          END WHERE
          xTalloc(:)    = MIN(1.0,MAX(0.0,Q10alloc** &
               ((casamet%tsoilavg(:)-TkzeroC-30.0)/10.0) )) !T limiting
          xNalloc(:)    = MIN(1.0,MAX(0.0,xwsalloc(:)*xTalloc(:)))     !N limiting
          xWalloc(:)    = MIN(1.0,MAX(0.0,casamet%btran(:)))           !W limiting
          xWorNalloc(:) = MIN(xWalloc(:),xNalloc(:))
          WHERE(casamet%lnonwood==0)
             casaflux%fracCalloc(:,FROOT) = R0 * 3.0 * xLalloc(:) &
                  / (xLalloc(:)+ 2.0*xWorNalloc(:))
             casaflux%fracCalloc(:,WOOD)  = S0 * 3.0 * xWorNalloc(:) &
                  / (2.0*xLalloc(:)+ xWorNalloc(:))
             casaflux%fracCalloc(:,LEAF)  = 1.0 - casaflux%fracCalloc(:,FROOT) &
                  - casaflux%fracCalloc(:,WOOD)
          ELSE WHERE
             casaflux%fracCalloc(:,FROOT) = R0 * 3.0 * xLalloc(:) &
                  / (xLalloc(:)+2.0*xWorNalloc(:))
             casaflux%fracCalloc(:,WOOD)  = 0.0
             casaflux%fracCalloc(:,LEAF)  = 1.0 - casaflux%fracCalloc(:,FROOT)
          END WHERE
       END WHERE
    CASE (0)   ! fixed allocation
       casaflux%fracCalloc(:,:) = casabiome%fracnpptop(veg%iveg(:),:)

    CASE (3) ! leaf:wood allocation set to maintain LA:SA ratio
       ! below target value of 4000, where phen%phase = 1 or 2
       !(requires casaflux%sapwood_area, which is inherited from the
       ! POP tree demography module. (Ticket #61)
       WHERE(casamet%lnonwood==0)
          casaflux%fracCalloc(:,FROOT) =  casabiome%fracnpptop(veg%iveg(:),FROOT)
          casaflux%fracCalloc(:,WOOD) = 0.01
          casaflux%fracCalloc(:,LEAF) = 1.0 - casaflux%fracCalloc(:,FROOT) - &
               casaflux%fracCalloc(:,WOOD)
          newLAI =casamet%glai + (casaflux%fracCalloc(:,LEAF) *casaflux%cnpp- &
               casaflux%kplant(:,leaf) *casapool%cplant(:,LEAF) )*casabiome%sla(veg%iveg(:))
          WHERE (casaflux%sapwood_area.GT.1.e-6 .AND. newLAI.GT.(4000.*casaflux%sapwood_area) &
               .AND. casaflux%cnpp.GT.0.0)

             casaflux%fracCalloc(:,LEAF) = ((4000.*casaflux%sapwood_area - casamet%glai)/ &
                  casabiome%sla(veg%iveg(:)) &
                  + casaflux%kplant(:,leaf) *casapool%cplant(:,LEAF)  )/casaflux%cnpp

             casaflux%fracCalloc(:,LEAF) = MAX(0.0,  casaflux%fracCalloc(:,LEAF) )
             casaflux%fracCalloc(:,LEAF) = MIN(1.0 - casaflux%fracCalloc(:,FROOT) - &
                  casaflux%fracCalloc(:,WOOD) ,&
                  casaflux%fracCalloc(:,LEAF) )

             casaflux%fracCalloc(:,WOOD) = 1.0 -  casaflux%fracCalloc(:,FROOT) - &
                  casaflux%fracCalloc(:,LEAF)
          END WHERE


       ELSEWHERE

          casaflux%fracCalloc(:,FROOT) =  casabiome%fracnpptop(veg%iveg(:),FROOT)
          casaflux%fracCalloc(:,WOOD) = 0.0
          casaflux%fracCalloc(:,LEAF) =  casabiome%fracnpptop(veg%iveg(:),LEAF)

       ENDWHERE

    END SELECT

991 FORMAT(1166(e14.7,2x))

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
    IF (LALLOC.NE.(3)) THEN

       WHERE(casamet%iveg2/=icewater)
          WHERE(phen%phase==0)
             casaflux%fracCalloc(:,leaf)  = 0.0
             casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                  /(casaflux%fracCalloc(:,froot) &
                  +casaflux%fracCalloc(:,wood))
             casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
          END WHERE

          WHERE(phen%phase==1)
             casaflux%fracCalloc(:,leaf)  = 0.8
             WHERE(casamet%lnonwood==0)  !woodland or forest
                casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
                casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
             ELSEWHERE !grassland
                casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
             ENDWHERE
          END WHERE

          WHERE(phen%phase==3)
             !      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
             casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,wood)
             casaflux%fracCalloc(:,leaf)    = 0.0
          ENDWHERE


          ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
          ! Q.Zhang 17/03/2011
          WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
             casaflux%fracCalloc(:,leaf)  = 0.0
             casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                  /(casaflux%fracCalloc(:,froot) &
                  +casaflux%fracCalloc(:,wood))
             casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
          ENDWHERE

          ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
          WHERE(casaflux%Cnpp<0.0)
             casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/SUM(casaflux%Crmplant,2)
             casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)/SUM(casaflux%Crmplant,2)
             casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/SUM(casaflux%Crmplant,2)
          ENDWHERE

          !! vh_js !!
          !! as long as biomass is positive, adjust allocation to be
          !! proportional to stock when NPP -ve   (Ticket#108)
          WHERE(casaflux%Cnpp<0.0 .AND. SUM(casapool%Cplant,2)>0  )
             casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)/SUM(casapool%Cplant,2)
             casaflux%fracCalloc(:,wood)  = casapool%Cplant(:,wood)/SUM(casapool%Cplant,2)
             casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot)/SUM(casapool%Cplant,2)
          ENDWHERE
       ENDWHERE

    ELSE
       WHERE(casamet%iveg2/=icewater)
          WHERE(phen%phase==0)
             casaflux%fracCalloc(:,leaf)  = 0.0
             casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                  /(casaflux%fracCalloc(:,froot) &
                  +casaflux%fracCalloc(:,wood))
             WHERE (casamet%lnonwood==0)
                casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
             ELSEWHERE
                casaflux%fracCalloc(:,wood) = 0.0
             ENDWHERE
          END WHERE

          WHERE(phen%phase==1.AND.casamet%lnonwood==1)

             casaflux%fracCalloc(:,leaf)  = 0.8
             casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
             casaflux%fracCalloc(:,wood) = 0.0
          ENDWHERE

          WHERE(phen%phase==3)
             !      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
             casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,wood)
             casaflux%fracCalloc(:,leaf)  = 0.0
          ENDWHERE

          !! vh !! don't require this fix for LALLOC = 3 (POP allocation scheme)
          !! Thiss fix can lead to over-allocation to roots, in turn bumping up N-uptake
          !! , leading to decline in mineral nitrogen availability and spikes in fracCalloc,
          !! causing spikes in tree mortality and lack of model convergence in productive
          !! regions where LAI is hitting LAImax.
!!$        ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
!!$        ! Q.Zhang 17/03/2011
!!$        WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
!!$           casaflux%fracCalloc(:,leaf)  = 0.0
!!$           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
!!$                /(casaflux%fracCalloc(:,froot) &
!!$                +casaflux%fracCalloc(:,wood))
!!$           WHERE (casamet%lnonwood==0)
!!$              casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
!!$           ELSEWHERE
!!$              casaflux%fracCalloc(:,wood) = 0.0
!!$           ENDWHERE
!!$        ENDWHERE

          WHERE(casamet%glai(:)<casabiome%glaimin(veg%iveg(:)))
             casaflux%fracCalloc(:,leaf)  = 0.8
             WHERE(casamet%lnonwood==0)  !woodland or forest
                casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
                casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
             ELSEWHERE !grassland
                casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
                casaflux%fracCalloc(:,wood) = 0.0
             ENDWHERE
          ENDWHERE
          ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
          WHERE(casaflux%Cnpp<0.0)
             WHERE(casamet%lnonwood==0)  !woodland or forest
                casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/SUM(casaflux%Crmplant,2)
                casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)/SUM(casaflux%Crmplant,2)
                casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/SUM(casaflux%Crmplant,2)
             ELSEWHERE
                casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/SUM(casaflux%Crmplant,2)
                casaflux%fracCalloc(:,wood)  = 0.0
                casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/SUM(casaflux%Crmplant,2)
             ENDWHERE
          ENDWHERE

          !! vh_js !!  Ticket#108
          WHERE(casaflux%Cnpp<0.0 .AND. SUM(casapool%Cplant,2)>0  )
             WHERE(casamet%lnonwood==0)  !woodland or forest
                casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)/SUM(casapool%Cplant,2)
                casaflux%fracCalloc(:,wood)  = casapool%Cplant(:,wood)/SUM(casapool%Cplant,2)
                casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot)/SUM(casapool%Cplant,2)
             ELSEWHERE
                casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)/SUM(casapool%Cplant,2)
                casaflux%fracCalloc(:,wood)  = 0.0
                casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot)/SUM(casapool%Cplant,2)
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




   ! normalization the allocation fraction to ensure they sum up to 1
   totfracCalloc(:) = SUM(casaflux%fracCalloc(:,:),2)
   casaflux%fracCalloc(:,leaf) = casaflux%fracCalloc(:,leaf)/totfracCalloc(:)
   casaflux%fracCalloc(:,wood) = casaflux%fracCalloc(:,wood)/totfracCalloc(:)
   casaflux%fracCalloc(:,froot) = casaflux%fracCalloc(:,froot)/totfracCalloc(:)

  END SUBROUTINE casa_allocation

  SUBROUTINE casa_wolf(veg,casabiome,casaflux,casapool,casamet)
  ! carbon allocation based on
  ! Wolf,Field and Berry, 2011. Ecological Applications, p1546-1556
  ! Wolf et al. 2011. Global Biogeochemical Cycles, 25, GB3015, doi:10.1019/2010GB003917
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(IN) :: casabiome
  TYPE (casa_met),            INTENT(IN) :: casamet
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux

    REAL, PARAMETER :: wolf_alpha1=6.22
    REAL, PARAMETER :: wolf_beta=-1.33
    REAL, PARAMETER :: wolf_c1=-wolf_alpha1/(1+wolf_beta)
    REAL, PARAMETER :: wolf_c2=1.0/(1.0+wolf_beta)

    REAL  totleaf,totwood,totcroot,totfroot,totnpp
    REAL  fracleaf,fracwood,fraccroot,fracfroot
    !
    ! local variables
    INTEGER   npt
    REAL(r_2), DIMENSION(mp)  ::  totbmdm,ntree,nppdm
    REAL(r_2), DIMENSION(mp)  ::  gleaf,gwood,gcroot,gfroot,gtot
    !
    ! input
    !  totleaf, totwood, totcroot, totfroot :    g C m-2
    !  totnpp:                                   g C m-2 d-1
    ! output
    !  fracleaf,fracwood, fraccroot, fracfroot:  fractions
    !

    DO npt=1,mp
       IF(casamet%iveg2(npt)==3.AND.casaflux%cnpp(npt)>0.0001) THEN  !forest types
          totbmdm(npt) = SUM(casapool%cplant(npt,:)) *10.0 / fracCbiomass      !10.0 for convert gc/m2 to kg/ha
          totbmdm(npt) = MAX(30000.0, totbmdm(npt))
          ! calculate tree stocking density
          ntree(npt) = 10**(wolf_c1+wolf_c2*LOG10(totbmdm(npt)))   ! tree ha-1, based on eqn (4) of Wolf et al. 2011, GBC
          ntree(npt) = MIN(200000.0,ntree(npt))
          ! changed by ypw 23/april/2012 to avoid negative npp
          nppdm(npt)  = (ABS(casaflux%cnpp(npt)) *365.0*0.001/fracCbiomass)/(0.0001*ntree(npt))  ! in kg dm tree-1 yr-1

          gleaf(npt)  = 0.156*(nppdm(npt)**1.106)     ! Figure 2a of Wolf, Field and Berry (2011)
          gwood(npt)  = 0.232*(nppdm(npt)**1.165)     ! Figure 2b of Wolf, Field and Berry (2011)
          gcroot(npt) = 0.0348*(nppdm(npt)**1.310)    ! Figure 2d of Wolf, Field and Berry (2011)
          gfroot(npt) = 0.247*(nppdm(npt)**0.987)     ! Figure 2c of Wolf, Field and Berry (2011)
          gtot(npt)   = gleaf(npt) + gwood(npt) + gcroot(npt) + gfroot(npt)

          casaflux%fracCalloc(npt,leaf)  = gleaf(npt)/gtot(npt)
          casaflux%fracCalloc(npt,wood)  = gwood(npt)/gtot(npt)
          casaflux%fracCalloc(npt,froot) = (gcroot(npt)+gfroot(npt))/gtot(npt)

          !        write(87,*) 'allocation = ',npt,casamet%iveg2(npt), totbmdm(npt),ntree(npt),nppdm(npt),casaflux%fracCalloc(npt,:)

       ELSE                ! other types
          casaflux%fracCalloc(npt,:) = casabiome%fracnpptop(veg%iveg(npt),:)
       ENDIF
    ENDDO

  END SUBROUTINE casa_wolf

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

    xkleafcold(:) = 0.0
    xkleafdry(:)  = 0.0
    xkleaf(:)     = 1.0

    ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
    DO npt=1,mp
       IF(casamet%iveg2(npt)/=icewater) THEN
          !    following the formulation of Arora (2005) on the
          !    effect of cold or drought stress on leaf litter fall
          !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
          IF(casamet%tairk(npt)>=phen%TKshed(veg%iveg(npt))) THEN
             xcoldleaf(npt) = 1.0
          ELSE
             IF(casamet%tairk(npt)<=(phen%TKshed(veg%iveg(npt))-5.0)) THEN
                xcoldleaf(npt)=0.0
             ELSE
                xcoldleaf(npt) = (casamet%tairk(npt)-phen%TKshed(veg%iveg(npt))-5.0)/5.0
             ENDIF
          ENDIF
          xcoldleaf(npt) = MIN(1.0,MAX(0.0,xcoldleaf(npt)))
          xkleafcold(npt) = casabiome%xkleafcoldmax(veg%iveg(npt)) &
               * (1.0-xcoldleaf(npt)) &
               ** casabiome%xkleafcoldexp(veg%iveg(npt))
          xkleafdry(npt)  = casabiome%xkleafdrymax(veg%iveg(npt)) &
               * (1.0-casamet%btran(npt))&
               ** casabiome%xkleafdryexp(veg%iveg(npt))
          IF (phen%phase(npt)==1) xkleaf(npt)= 0.0
          ! vh: account for high rate of leaf loss during senescence
          ! vh_js
          IF (TRIM(cable_user%PHENOLOGY_SWITCH)=='climate') THEN
             IF (phen%phase(npt)==3.OR.phen%phase(npt)==0) xkleaf(npt)= 100.0
          ENDIF
       END IF
    END DO

    !  WHERE(casamet%iveg2/=icewater)
    !  !    following the formulation of Arora (2005) on the
    !  !    effect of cold or drought stress on leaf litter fall
    !  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
    !    WHERE(casamet%tairk(:)>=phen%TKshed(veg%iveg(:)))
    !      xcoldleaf(:) = 1.0
    !    ELSEWHERE
    !      WHERE(casamet%tairk(:)<=(phen%TKshed(veg%iveg(:))-5.0))
    !        xcoldleaf(:)=0.0
    !      ELSEWHERE
    !        xcoldleaf(:) = (casamet%tairk(:)-phen%TKshed(veg%iveg(:))-5.0)/5.0
    !      ENDWHERE
    !    ENDWHERE
    !    xcoldleaf(:) = min(1.0,max(0.0,xcoldleaf(:)))
    !    xkleafcold(:) = casabiome%xkleafcoldmax(veg%iveg(:)) * (1.0-xcoldleaf(:)) &
    !                 ** casabiome%xkleafcoldexp(veg%iveg(:))
    !    xkleafdry(:) = casabiome%xkleafdrymax(veg%iveg(:))*(1.0-casamet%btran(:))&
    !                 ** casabiome%xkleafdryexp(veg%iveg(:))
    !    WHERE(phen%phase(:)==1) xkleaf(:)= 0.0
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
    REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xklitter,xksoil
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome

    ! local variables
    INTEGER nland,np
    REAL(r_2), PARAMETER :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps
    REAL(r_2), PARAMETER :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
    REAL(r_2), PARAMETER :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
    REAL(r_2), PARAMETER :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
    REAL(r_2), PARAMETER :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)
    ! Kirschbaum function parameters
    REAL(r_2), PARAMETER :: xkalpha=-3.764   ! Kirschbaum (1995, SBB)
    REAL(r_2), PARAMETER :: xkbeta=0.204
    REAL(r_2), PARAMETER :: xktoptc=36.9
    REAL(r_2), DIMENSION(mp)       :: xkwater,xktemp
    REAL(r_2), DIMENSION(mp)       :: fwps,tsavg
    ! Custom soil respiration - see Ticket #42
    REAL(r_2), DIMENSION(mp)       :: smrf,strf,slopt,wlt,tsoil,fcap,sopt
    !,tsurfavg  !!, msurfavg
    INTEGER :: npt

    xklitter(:) = 1.0
    xksoil(:)   = 1.0
    fwps(:)     =  casamet%moistavg(:)/soil%ssat(:)
    tsavg(:)    =  casamet%tsoilavg(:)

    ! Custom soil respiration - see Ticket #42
    tsoil(:)    =  tsavg(:)-TKzeroC !tsoil in C
    strf(:)     = 1.0
    smrf(:)     = 1.0
    slopt(:)    = 1.0
    sopt(:)     = 1.0


    ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
    DO npt=1,mp
       IF(casamet%iveg2(npt)/=icewater) THEN
          xktemp(npt)  = casabiome%q10soil(veg%iveg(npt))**(0.1*(tsavg(npt)-TKzeroC-35.0))
          xkwater(npt) = ((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
               * ((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
          IF (veg%iveg(npt) == cropland .OR. veg%iveg(npt) == croplnd2) &
               xkwater(npt)=1.0
          xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) * xktemp(npt) * xkwater(npt)

          IF( .NOT. cable_user%SRF) THEN
             ! Use original function, ELSE Ticket #42
             xksoil(npt)   = casabiome%xkoptsoil(veg%iveg(npt))   * xktemp(npt) * xkwater(npt)
          ELSE
             ! Custom soil respiration - see Ticket #42
             ! Implementing alternative parameterizations
             IF(TRIM(cable_user%SMRF_NAME)=='CASA-CNP') THEN
                smrf(npt)=xkwater(npt)
             ELSE IF (TRIM(cable_user%SMRF_NAME)=='SOILN') THEN
                sopt(npt)=0.92
                slopt(npt)=wlt(npt)+0.1          !SLOPT is the lower optimum
                IF (fwps(npt)>sopt(npt)) THEN
                   smrf(npt)=0.2+0.8*(1.0-fwps(npt))/(1.0-sopt(npt))
                ELSE IF(slopt(npt)<=fwps(npt) .AND. fwps(npt)<=sopt(npt)) THEN
                   smrf(npt) = 1.0
                ELSE IF (wlt(npt)<=fwps(npt) .AND. fwps(npt) <slopt(npt)) THEN
                   smrf(npt)=0.01+0.99*(fwps(npt)-wlt(npt))/(slopt(npt)-wlt(npt))
                ELSE IF (fwps(npt)<wlt(npt)) THEN
                   smrf(npt) = 0.01
                END IF
             ELSE IF (TRIM(cable_user%SMRF_NAME)=='TRIFFID') THEN
                sopt(npt) = 0.5 * (1+wlt(npt))
                IF (fwps(npt) > sopt(npt)) THEN
                   smrf(npt) =1.0-0.8*(fwps(npt)-sopt(npt))
                ELSE IF (wlt(npt)<fwps(npt) .AND. fwps(npt)<=sopt(npt)) THEN
                   smrf(npt)=0.01+0.8*((fwps(npt)-wlt(npt))/(sopt(npt)-wlt(npt)))
                ELSE IF (fwps(npt)<wlt(npt)) THEN
                   smrf(npt) = 0.2
                END IF
             END IF

             IF(TRIM(cable_user%STRF_NAME)=='CASA-CNP') THEN
                strf(npt)=xktemp(npt)
             ELSE IF (TRIM(cable_user%STRF_NAME)=='K1995') THEN
                !Kirschbaum from Kirschbaum 1995, eq (4) in SBB, .66 is to collapse smrf
                !to same area
                strf(npt)=EXP(-3.764+0.204*tsoil(npt)*(1-0.5*tsoil(npt)/36.9))/.66
             ELSE IF (TRIM(cable_user%STRF_NAME)=='PnET-CN') THEN
                strf(npt)=0.68*EXP(0.1*(tsoil(npt)-7.1))/12.64
             END IF
             xksoil(npt) = casabiome%xkoptsoil(veg%iveg(npt))*strf(npt)*smrf(npt)
          END IF
       END IF
    END DO

  END SUBROUTINE casa_xratesoil

  SUBROUTINE casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
       casaflux,casamet,phen)
    ! calculate the plant litter fall rate, litter fall and sOM decomposition rate (1/day)
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

    IMPLICIT NONE
    REAL(r_2), DIMENSION(mp), INTENT(IN)    :: xkleafcold,xkleafdry,xkleaf
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (phen_variable),       INTENT(IN) :: phen

    ! local variables
    REAL(r_2), DIMENSION(mp)  :: xk
    REAL(r_2), DIMENSION(mp,mplant)         :: ratioLignintoN
    INTEGER npt

    casaflux%fromPtoL(:,:,:)      = 0.0
    casaflux%kplant(:,:)          = 0.0   ! (BPjun2010)

    WHERE(casamet%iveg2/=icewater)
       ! using max function to avoid dividing by zero, ypw 14/may/2008
       ratioLignintoN(:,leaf) = (casapool%Cplant(:,leaf) &
            /(MAX(1.0e-10,casapool%Nplant(:,leaf)) *casabiome%ftransNPtoL(veg%iveg(:),leaf))) &
            * casabiome%fracLigninplant(veg%iveg(:),leaf)
       ratioLignintoN(:,froot)= (casapool%Cplant(:,froot)&
            /(MAX(1.0e-10,casapool%Nplant(:,froot))*casabiome%ftransNPtoL(veg%iveg(:),froot))) &
            * casabiome%fracLigninplant(veg%iveg(:),froot)

       casaflux%fromPtoL(:,metb,leaf)    = MAX(0.001, 0.85 - 0.018 *ratioLignintoN(:,leaf))
       casaflux%fromPtoL(:,metb,froot)   = MAX(0.001, 0.85 - 0.018 *ratioLignintoN(:,froot))
       casaflux%fromPtoL(:,str,leaf)    = 1.0 - casaflux%fromPtoL(:,metb,leaf)
       casaflux%fromPtoL(:,str,froot)   = 1.0 - casaflux%fromPtoL(:,metb,froot)
       casaflux%fromPtoL(:,cwd,wood)    = 1.0

       ! calc. of casaflux%kplant drops scaling - see #242
       casaflux%kplant(:,leaf)        = casabiome%plantrate(veg%iveg(:),leaf)

       casaflux%kplant(:,wood)        = casabiome%plantrate(veg%iveg(:),wood)
       casaflux%kplant(:,froot)       = casabiome%plantrate(veg%iveg(:),froot)
    ENDWHERE


    ! When glai<glaimin,leaf biomass will not decrease anymore. (Q.Zhang 10/03/2011)
    DO npt = 1,mp
       IF(casamet%glai(npt).LE.casabiome%glaimin(veg%iveg(npt))) casaflux%kplant(npt,leaf) = 0.0
    ENDDO
    ! end change

  END SUBROUTINE casa_coeffplant

  SUBROUTINE casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
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

    IMPLICIT NONE
    REAL(r_2), DIMENSION(mp), INTENT(IN) :: xklitter,xksoil
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet

    ! local variables
    INTEGER j,k,kk,nland             !i: for plant pool, j for litter, k for soil

    casaflux%fromLtoS(:,:,:)      = 0.0
    casaflux%fromStoS(:,:,:)      = 0.0
    ! flow from soil to soil
    DO k = 1, msoil
       casaflux%fromStoS(:,k,k)   = -1.0
    ENDDO  ! "k"
    casaflux%fromLtoCO2(:,:) = 0.0             ! flow from L or S to CO2
    casaflux%fromStoCO2(:,:) = 0.0

    casaflux%klitter(:,:) = 0.0        !initialize klitter (Q.Zhang 03/03/2011)

    WHERE(casamet%iveg2/=icewater)

       casaflux%klitter(:,metb)   = xklitter(:) * casabiome%litterrate(veg%iveg(:),metb)
       casaflux%klitter(:,str)    = xklitter(:) * casabiome%litterrate(veg%iveg(:),str) &
            * EXP(-3.0*casabiome%fracLigninplant(veg%iveg(:),leaf))
       casaflux%klitter(:,cwd)    = xklitter(:) * casabiome%litterrate(veg%iveg(:),cwd)

       casaflux%ksoil(:,mic)      = xksoil(:) * casabiome%soilrate(veg%iveg(:),mic)   &
            * (1.0 - 0.75 *(soil%silt(:)+soil%clay(:)))
       casaflux%ksoil(:,slow)     = xksoil(:) * casabiome%soilrate(veg%iveg(:),slow)
       casaflux%ksoil(:,pass)     = xksoil(:) * casabiome%soilrate(veg%iveg(:),pass)
       casaflux%kplab(:)          = xksoil(:) * casabiome%xkplab(casamet%isorder(:))
       casaflux%kpsorb(:)         = xksoil(:) * casabiome%xkpsorb(casamet%isorder(:))
       casaflux%kpocc(:)          = xksoil(:) * casabiome%xkpocc(casamet%isorder(:))


       WHERE(veg%iveg==cropland)      ! for cultivated land type
          casaflux%ksoil(:,mic)  = casaflux%ksoil(:,mic) * 1.25
          casaflux%ksoil(:,slow) = casaflux%ksoil(:,slow)* 1.5
          casaflux%ksoil(:,pass) = casaflux%ksoil(:,pass)* 1.5
       ENDWHERE  !

       ! flow from litter to soil
       casaflux%fromLtoS(:,mic,metb)   = 0.45
       ! metb -> mic
       casaflux%fromLtoS(:,mic,str)   = 0.45*(1.0-casabiome%fracLigninplant(veg%iveg(:),leaf))
       ! str -> mic
       casaflux%fromLtoS(:,slow,str)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),leaf)
       ! str -> slow
       casaflux%fromLtoS(:,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(:),wood))
       ! CWD -> fmic
       casaflux%fromLtoS(:,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),wood)
       ! CWD -> slow

       !! set the following two backflow to set (see Bolker 199x)
       !    casaflux%fromStoS(:,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(:))
       !    casaflux%fromStoS(:,mic,pass)  = 0.45

       casaflux%fromStoS(:,slow,mic)  = (0.85 - 0.68 * (soil%clay(:)+soil%silt(:))) &
            * (0.997 - 0.032*soil%clay(:))
       casaflux%fromStoS(:,pass,mic)  = (0.85 - 0.68 * (soil%clay(:)+soil%silt(:))) &
            * (0.003 + 0.032*soil%clay(:))
       casaflux%fromStoS(:,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(:) )

    ENDWHERE

    DO nland=1,mp
       IF(casamet%iveg2(nland)/=icewater) THEN
          DO j=1,mlitter
             DO k=1,msoil
                casaflux%fromLtoCO2(nland,j) = casaflux%fromLtoCO2(nland,j)  &
                     + casaflux%fromLtoS(nland,k,j)
             ENDDO  !"k"
             casaflux%fromLtoCO2(nland,j) = 1.0 - casaflux%fromLtoCO2(nland,j)
          ENDDO !"j"
          DO k=1,msoil
             DO kk=1,msoil
                casaflux%fromStoCO2(nland,k) = casaflux%fromStoCO2(nland,k) &
                     + casaflux%fromStoS(nland,kk,k)
             ENDDO  !"kk"
          ENDDO   !"k"
          casaflux%fromStoCO2(nland,:) = -casaflux%fromStoCO2(nland,:)
       ENDIF
    ENDDO   ! "nland"

  END SUBROUTINE casa_coeffsoil

  ! modified by ypw following Chris Lu 5/nov/2012
  SUBROUTINE casa_delplant(veg,casabiome,casapool,casaflux,casamet,            &
     cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

  !  calculate the chnage in plant C, N and P pools
  !  uptake of N and P will be computed in casa_uptake
  !  labile C pool will be computed casa_labile

    IMPLICIT NONE
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_biome),         INTENT(INOUT) :: casabiome
    TYPE (casa_pool),          INTENT(INOUT) :: casapool
    TYPE (casa_flux),          INTENT(INOUT) :: casaflux
    TYPE (casa_met),           INTENT(INOUT) :: casamet

    ! added by ypwang following Chris Lu 5/nov/2012
    REAL, DIMENSION(mp),INTENT(OUT) :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

    INTEGER  npt,nL,nP,nland
    REAL(r_2)      :: Ygrow, ratioPNplant

    casaflux%FluxCtolitter = 0.0
    casaflux%FluxNtolitter = 0.0
    casaflux%FluxPtolitter = 0.0

    ! added by ypwang following Chris Lu 5/nov/2012

    cleaf2met = 0.0
    cleaf2str = 0.0
    croot2met = 0.0
    croot2str = 0.0
    cwood2cwd = 0.0

    nleaf2met = 0.0
    nleaf2str = 0.0
    nroot2met = 0.0
    nroot2str = 0.0
    nwood2cwd = 0.0

    pleaf2met = 0.0
    pleaf2str = 0.0
    proot2met = 0.0
    proot2str = 0.0
    pwood2cwd = 0.0

    !MPI
    DO npt=1,mp
       IF(casamet%iveg2(npt)/=icewater) THEN
          !    PRINT *, 'npt = ', npt
          !    PRINT *, 'casapool%cplant(npt,:) = ', casapool%cplant(npt,:)
          casapool%dcplantdt(npt,:)  =  casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)     &
               - casaflux%kplant(npt,:)  * casapool%cplant(npt,:)



          !casapool%dcplantdt(npt,2) = casapool%dcplantdt(npt,2)

          !! vh_js !!
          !! adjust turnover and autotrophic respiration to avoid negative stores.
          !! Ticket#108

          WHERE (((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)).LT. 0.0) &
               .OR. ((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) &
               .LT. 0.5 * casapool%cplant(npt,2:3) ))
             casaflux%kplant(npt,2:3) = 0.0
             casaflux%crmplant(npt,2:3)= 0.0
          endwhere
          IF(ANY((casapool%dcplantdt(npt,:)*deltpool + casapool%cplant(npt,:)).LT. 0.0)) THEN
             casaflux%kplant(npt,1) = 0.0
             casaflux%crmplant(npt,1)= MIN(casaflux%crmplant(npt,1),0.5*casaflux%Cgpp(npt))
          ENDIF

          !! revise turnover and NPP and dcplantdt to reflect above adjustments

          casaflux%Cplant_turnover(npt,:) = casaflux%kplant(npt,:)  * casapool%cplant(npt,:)
          IF (ANY((casapool%dcplantdt(npt,:)*deltpool + casapool%cplant(npt,:)).LT. 0.0) &

               .OR. ANY((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) &
               .LT. 0.5 * casapool%cplant(npt,2:3) )) THEN

             ratioPNplant = 0.0
             IF (casapool%Nplant(npt,leaf)>0.0) &
                  ratioPNplant = casapool%Pplant(npt,leaf)/(casapool%Nplant(npt,leaf)+ 1.0e-10)

             Ygrow = 0.65+0.2*ratioPNplant/(ratioPNplant+1.0/15.0)
             IF ((casaflux%Cgpp(npt)-SUM(casaflux%crmplant(npt,:)))>0.0) THEN
                ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
                casaflux%crgplant(npt)  = (1.0-Ygrow)* MAX(0.0,casaflux%Cgpp(npt)- &
                     SUM(casaflux%crmplant(npt,:)))
             ELSE
                casaflux%crgplant(npt) = 0.0
             ENDIF
             casaflux%Cnpp(npt) = casaflux%Cgpp(npt)-SUM(casaflux%crmplant(npt,:)) &
                  - casaflux%crgplant(npt) - casaflux%fracClabile(npt) * casaflux%cgpp(npt)

             casapool%dcplantdt(npt,:)  =  casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)     &
                  - casaflux%kplant(npt,:)  * casapool%cplant(npt,:)

          ENDIF


          !! vh_js !! end of adjustments to avoid negative stores Ticket#108

          ! change here made by ypw on 26august 2011
          ! calculate fraction c to labile pool as a fraction of gpp, not npp
          ! casapool%dClabiledt(npt)   = casaflux%Cnpp(npt)    * casaflux%fracClabile(npt)
          casapool%dClabiledt(npt)   =  casaflux%Cgpp(npt)  * casaflux%fracClabile(npt) &
               - casaflux%clabloss(npt)
          ! added by ypwang 5/nov/2012
          cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          croot2str(npt) = casaflux%fromPtoL(npt,str,froot)  * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          cwood2cwd(npt) = casaflux%fromPtoL(npt,cwd,wood)   * casaflux%kplant(npt,wood)  * casapool%cplant(npt,wood)

          !    PRINT *, 'npt, mp, iveg', npt, mp, veg%iveg(npt)
          IF(icycle > 1) THEN
             !    PRINT *, 'casapool%Nplant(npt,:) = ', casapool%Nplant(npt,:)
             IF(casaflux%fracNalloc(npt,leaf)==0.0) THEN
                casapool%dNplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf)
             ELSE
                casapool%dNplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf) &
                     * casabiome%ftransNPtoL(veg%iveg(npt),leaf)
             ENDIF

             IF (casamet%lnonwood(npt)==0) THEN
                casapool%dNplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Nplant(npt,wood) &
                     * casabiome%ftransNPtoL(veg%iveg(npt),wood)
             ELSE
                casapool%dNplantdt(npt,wood) = 0.0
             ENDIF

             casapool%dNplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Nplant(npt,froot) &
                  * casabiome%ftransNPtoL(veg%iveg(npt),froot)
             ! added by ypwang 5/nov/2012

             nleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                  * casapool%cplant(npt,leaf)       * ratioNCstrfix
             nroot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                  * casapool%cplant(npt,froot)      * ratioNCstrfix

             nleaf2met(npt) = - casapool%dNplantdt(npt,leaf)  - nleaf2str(npt)
             nroot2met(npt) = - casapool%dNplantdt(npt,froot) - nroot2str(npt)

             nwood2cwd(npt) = -casapool%dNplantdt(npt,wood)

          ENDIF


          IF(icycle >2) THEN

             IF(casaflux%fracPalloc(npt,leaf)==0.0) THEN
                casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf)
             ELSE
                casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf) &
                     * casabiome%ftransPPtoL(veg%iveg(npt),leaf)
             ENDIF

             casapool%dPplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Pplant(npt,wood) &
                  * casabiome%ftransPPtoL(veg%iveg(npt),wood)


             casapool%dPplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Pplant(npt,froot) &
                  * casabiome%ftransPPtoL(veg%iveg(npt),froot)
             ! added by ypwang 5/nov/2012

             pleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                  * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix
             proot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                  * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix
             pleaf2met(npt) = -casapool%dPplantdt(npt,leaf)  - pleaf2str(npt)
             proot2met(npt) = -casapool%dPplantdt(npt,froot) - proot2str(npt)
             pwood2cwd(npt) = -casapool%dPplantdt(npt,wood)


          ENDIF

          DO nL=1,mlitter
             DO nP=1,mplant
                casaflux%FluxCtolitter(npt,nL) = casaflux%FluxCtolitter(npt,nL) &
                     + casaflux%fromPtoL(npt,nL,nP) &
                     * casaflux%kplant(npt,nP) &
                     * casapool%cplant(npt,nP)
             ENDDO
          ENDDO

          !    PRINT *, 'before 2nd icycle >1; npt, mp', npt, mp
          IF(icycle > 1) THEN
!!$       casaflux%FluxNtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
!!$                               * casapool%cplant(npt,leaf)       * ratioNCstrfix              &
!!$                               + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
!!$                               * casapool%cplant(npt,froot)      * ratioNCstrfix

             !vh! to avoid -ve Nitrogen pools Ticket#108
             casaflux%FluxNtolitter(npt,str) = MIN(casaflux%fromPtoL(npt,str,leaf) * &
                  casaflux%kplant(npt,leaf)  &
                  * casapool%cplant(npt,leaf)       * ratioNCstrfix &
                  , -casapool%dNplantdt(npt,leaf))             &
                  + MIN(casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                  * casapool%cplant(npt,froot)      * ratioNCstrfix &
                  , -casapool%dNplantdt(npt,froot))

             casaflux%FluxNtolitter(npt,metb) = - casapool%dNplantdt(npt,leaf)-casapool%dNplantdt(npt,froot) &
                  - casaflux%FluxNtolitter(npt,str)
             casaflux%FluxNtolitter(npt,CWD) = -casapool%dNplantdt(npt,wood)

             ! adding N uptake
             casapool%dNplantdt(npt,:) = casapool%dNplantdt(npt,:) &
                  + casaflux%Nminuptake(npt)*casaflux%fracNalloc(npt,:)
             !       casapool%Nsoilmin(npt)    = casapool%Nsoilmin(npt) - casaflux%Nminuptake(npt) *deltpool
          ENDIF !end "icycle >1"


          IF(icycle>2) THEN
             casaflux%FluxPtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                  * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix        &
                  + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                  * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix
             casaflux%FluxPtolitter(npt,metb) = -casapool%dPplantdt(npt,leaf)-casapool%dPplantdt(npt,froot) &
                  - casaflux%FluxPtolitter(npt,str)
             casaflux%FluxPtolitter(npt,CWD) = -casapool%dPplantdt(npt,wood)
             ! add P uptake
             casapool%dPplantdt(npt,:) = casapool%dPplantdt(npt,:) &
                  + casaflux%Plabuptake(npt)*casaflux%fracPalloc(npt,:)
             !       casapool%Psoillab(npt)    = casapool%Psoillab(npt) - casaflux%Plabuptake(npt) * deltpool
          ENDIF  !of "icycle >2"


       ENDIF
    ENDDO
    npt=2
    ! write(911,91) casaflux%kplant(npt,wood), casabiome%ftransNPtoL(veg%iveg(npt),wood), casapool%Nplant(npt,wood), casapool%cplant(npt,wood), &
    !casaflux%Nminuptake(npt)*casaflux%fracNalloc(npt,wood), casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,wood)
91  FORMAT (100(e12.4,2x))

  END SUBROUTINE casa_delplant

  SUBROUTINE casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
    ! calculate changes in litter and soil pools

    IMPLICIT NONE
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome

    ! local variables
    REAL(r_2), DIMENSION(mp)    :: xdplabsorb, fluxptase
    INTEGER i,j,jj,k,kk,kkk,n,iv,npt,nL,nS,nSS,nland
    ! local temporary ypwang
    real latx, lonx
    integer ivegx,isox

    casaflux%fluxCtoCO2    = 0.0
    casaflux%fluxCtosoil   = 0.0
    casaflux%fluxNtosoil   = 0.0
    casaflux%fluxPtosoil   = 0.0
    casaflux%Crsoil        = 0.0  ! initialization added by BP jul2010

    casapool%dClitterdt = 0.0
    casapool%dCsoildt   = 0.0

    casapool%dNlitterdt = 0.0
    casapool%dNsoildt   = 0.0
    casapool%dNsoilmindt= 0.0
    casaflux%Nsmin = 0.0
    casaflux%Nsimm = 0.0
    casaflux%Nsnet = 0.0
    casaflux%Nminloss = 0.0
    casaflux%Nminleach = 0.0
    casaflux%Nlittermin=0.0

    casapool%dPlitterdt   = 0.0
    casapool%dPsoildt     = 0.0
    casapool%dPsoillabdt  = 0.0
    casapool%dPsoilsorbdt = 0.0
    casapool%dPsoiloccdt  = 0.0
    casaflux%Psmin = 0.0
    casaflux%Psimm = 0.0
    casaflux%Psnet = 0.0
    casaflux%Pleach = 0.0
    casaflux%Ploss  = 0.0
    casaflux%Plittermin=0.0
    fluxptase = 0.0

  DO nland=1,mp
       IF(casamet%iveg2(nland)/=icewater) THEN

          IF(icycle > 1) THEN
             !vh! set klitter to zero where Nlitter will go -ve
             !(occurs occasionally for metabolic litter pool) Ticket#108
             WHERE (casaflux%klitter(nland,:) * MAX(0.0,casapool%Nlitter(nland,:)).GT. &
                  casapool%Nlitter(nland,:)+casaflux%fluxNtolitter(nland,:)) casaflux%klitter(nland,:) = 0.0
          ENDIF

          DO nL=1,mlitter
             casaflux%fluxCtoCO2(nland) = casaflux%fluxCtoCO2(nland)  &
                  + casaflux%fromLtoCO2(nland,nL)  &
                  * casaflux%klitter(nland,nL) &
                  * casapool%clitter(nland,nL)
          ENDDO

          DO nS=1,msoil
             DO nL=1,mlitter
                casaflux%fluxCtosoil(nland,nS) = casaflux%fluxCtosoil(nland,nS) &
                     + casaflux%fromLtoS(nland,nS,nL) &
                     * casaflux%klitter(nland,nL) &
                     * casapool%clitter(nland,nL)
             ENDDO
             DO nSS=1,msoil
                IF(nSS/=nS) THEN
                   casaflux%fluxCtosoil(nland,nS) = casaflux%fluxCtosoil(nland,nS) &
                        + casaflux%fromStoS(nland,nS,nSS) &
                        * casaflux%ksoil(nland,nSS) &
                        * casapool%csoil(nland,nSS)
                ENDIF
             ENDDO
             casaflux%fluxCtoCO2(nland) = casaflux%fluxCtoCO2(nland)  &
                  + casaflux%fromStoCO2(nland,nS) &
                  * casaflux%ksoil(nland,nS) &
                  * casapool%csoil(nland,nS)
          ENDDO

          IF(icycle>1) THEN
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
                   IF(kkk.NE.kk) THEN
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
                     * MAX(0.0,casaflux%Nsnet(nland))
                casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                     * MAX(0.0,casapool%Nsoilmin(nland))
             ELSE
                casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                     * MAX(0.0,casaflux%Nsnet(nland)) &
                     * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
                casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                     * MAX(0.0,casapool%Nsoilmin(nland)) &
                     * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
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
                   IF(kk.NE.k) THEN
                      casaflux%FluxNtosoil(nland,k) = casaflux%FluxNtosoil(nland,k)  &
                           + casaflux%fromStoS(nland,k,kk) &
                           * casaflux%ksoil(nland,kk)      &
                           * casapool%Csoil(nland,kk)      &
                           * casapool%ratioNCsoilnew(nland,k)
                   ENDIF
                ENDDO ! end of "kk"
             ENDDO    ! end of "k"

          ENDIF !end of icycle >1

          IF(icycle >2) THEN
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
                        * casapool%Clitter(nland,jj)     &
                        * casapool%ratioPCsoil(nland,kk)
                   !                                     * casapool%ratioPCsoil(nland,kk)/casapool%ratioNCsoil(nland,kk)
                ENDDO
                DO kkk=1,msoil      ! immobilisation from soil to soil
                   IF(kkk.NE.kk) THEN
                      casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                           - casaflux%fromStoS(nland,kk,kkk)  &
                           * casaflux%ksoil(nland,kkk) &
                           * casapool%Csoil(nland,kkk) &
                           * casapool%ratioPCsoil(nland,kk)
                      !                                        * casapool%ratioPCsoil(nland,kk)/casapool%ratioNCsoil(nland,kk)
                   ENDIF
                ENDDO
             ENDDO  ! immobilization

             casaflux%Psnet(nland)=casaflux%Plittermin(nland) &
                  +casaflux%Psmin(nland)   &
                  +casaflux%Psimm(nland)
             ! net mineralization

             !      casaflux%Pleach(nland)  =  (1.0e-4) &
             !                                 * max(0.0,casapool%Psoillab(nland))

             casaflux%Pleach(nland)  =  casaflux%fPleach(nland) &
                  * MAX(0.0,casapool%Psoillab(nland))

             DO k=1,msoil
                DO j=1,mlitter
                   casaflux%FluxPtosoil(nland,k) =  casaflux%FluxPtosoil(nland,k)  &
                        + casaflux%fromLtoS(nland,k,j) &
                        * casaflux%klitter(nland,j)    &
                        * casapool%Clitter(nland,j)    &
                        * casapool%ratioPCsoil(nland,k)
                   !                                 * casapool%ratioPCsoil(nland,k)/casapool%ratioNCsoil(nland,k)
                ENDDO  ! end of "j"
                DO kk=1,msoil
                   IF(kk.NE.k) THEN
                      !               casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
                      !                                    + casaflux%fromStoS(nland,k,kk) &
                      !                                    * casaflux%ksoil(nland,kk)      &
                      !                                    * casapool%Csoil(nland,kk)      &
                      !                                    * casapool%ratioPCsoil(nland,k)
                      casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
                           + casaflux%fromStoS(nland,k,kk) &
                           * casaflux%ksoil(nland,kk)      &
                           * casapool%Csoil(nland,kk)      &
                           * casapool%ratioPCsoil(nland,k)
                   !  changed using PC ratio, ypw 
                   !        * casapool%Csoil(nland,kk)      &
                   !        /casapool%ratioNPsoil(nland,k)
                      !                                    * casapool%ratioPCsoil(nland,k)/casapool%ratioNCsoil(nland,k)
                   ENDIF
                ENDDO ! end of "kk"
             ENDDO    ! end of "k"
             ! need to account for flow from sorbed to occluded pool
          ENDIF
  ENDIF  ! end of /=icewater
    ENDDO  ! end of nland

    DO nland=1,mp
       IF(casamet%iveg2(nland)/=icewater) THEN
          casapool%dClitterdt(nland,:) =  casaflux%fluxCtolitter(nland,:) - casaflux%klitter(nland,:) * casapool%clitter(nland,:)
          casapool%dCsoildt(nland,:)   =  casaflux%fluxCtosoil(nland,:)   - casaflux%ksoil(nland,:)   * casapool%csoil(nland,:)
          casaflux%Crsoil(nland)       =  casaflux%fluxCtoCO2(nland)
          casaflux%cnep(nland)         =  casaflux%cnpp(nland) - casaflux%Crsoil(nland)

          IF(icycle > 1) THEN
             casapool%dNlitterdt(nland,:) =  casaflux%fluxNtolitter(nland,:)  &
                  - casaflux%klitter(nland,:) &
                  * MAX(0.0,casapool%Nlitter(nland,:))

             casapool%dNsoildt(nland,:) = casaflux%FluxNtosoil(nland,:) &
                  - casaflux%ksoil(nland,:) * casapool%Nsoil(nland,:)
             casapool%dNsoilmindt(nland)= casaflux%Nsnet(nland)&
                  + casaflux%Nmindep(nland) + casaflux%Nminfix(nland)   &
                  - casaflux%Nminloss(nland)   &
                  - casaflux%Nminleach(nland)   &
                  - casaflux%Nupland(nland)

          ENDIF
!!$if (nland==1) write(59,91) casaflux%Nsnet(nland) , &
!!$                                 casaflux%Nlittermin(nland),  &
!!$                                 casaflux%Nsmin(nland),   &
!!$                                 casaflux%Nsimm(nland)
!!$                                 , casaflux%Nmindep(nland) ,casaflux%Nminfix(nland)   &
!!$                                 , casaflux%Nminloss(nland)   &
!!$                                 , casaflux%Nminleach(nland)   &
!!$                                 , casaflux%Nupland(nland)

91        FORMAT(20(e12.4,2x))
          IF(icycle >2) THEN

          !  doubling account of "deltcasa" wit "readbiome"   
             fluxptase(nland) =  casabiome%prodptase( veg%iveg(nland) )       &
                  * MAX( 0.0_r_2, ( casapool%Psoil(nland,2)                   &
                  * casaflux%ksoil(nland,2)                &
                  + casapool%Psoil(nland,3)                &
                  * casaflux%ksoil(nland,3) )              &
                  )                                                 &
                  * MAX( 0.0_r_2, ( casabiome%costNPup( veg%iveg(nland) )    &
                  - 15.0 )                                 &
                  )                                                 &
                  / ( MAX( 0.0_r_2, ( casabiome%costNPup( veg%iveg(nland) )  &
                  - 15.0 )                               &
                  ) + 150.0                                       &
                  )

             xdplabsorb(nland) = 1.0+ casaflux%Psorbmax(nland)*casaflux%kmlabp(nland) &
                  /((casaflux%kmlabp(nland)+casapool%Psoillab(nland))**2)
             casapool%dPlitterdt(nland,:) = casaflux%fluxPtolitter(nland,:)  &
                  - casaflux%klitter(nland,:)                 &
                  * MAX(zero,casapool%Plitter(nland,:))

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
             casapool%dPsoillabdt(nland)  = casapool%dPsoillabdt(nland)/xdplabsorb(nland)
             casapool%dPsoilsorbdt(nland) = 0.0

             casapool%dPsoiloccdt(nland)  = casaflux%kpsorb(nland)* casapool%Psoilsorb(nland) &
                  - casaflux%kpocc(nland) * casapool%Psoilocc(nland)
             ! P loss to non-available P pools
             !      casaflux%Ploss(nland)        = casaflux%kpocc(nland) * casapool%Psoilocc(nland)

             !      casaflux%Ploss(nland)       = casaflux%fPleach(nland) &
             !                                 * max(zero,casapool%Psoillab(nland))
             casaflux%Ploss(nland)       = 0.0
          ENDIF
  ENDIF
    ENDDO

  END SUBROUTINE casa_delsoil

  SUBROUTINE avgsoil(veg,soil,casamet)
    ! Get avg soil moisture, avg soil temperature
    ! need to estimate the land cell mean soil temperature and moisture weighted by the area fraction
    ! of each tile within the land cell

    IMPLICIT NONE
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_met),              INTENT(INOUT) :: casamet

    INTEGER                     :: ns,nland

    casamet%tsoilavg   = 0.0
    casamet%moistavg   = 0.0
    casamet%btran      = 0.0

    DO ns = 1, ms
       DO nland=1,mp
          casamet%tsoilavg(nland)  = casamet%tsoilavg(nland)+veg%froot(nland,ns)  &
               * casamet%tsoil(nland,ns)
          casamet%moistavg(nland)  = casamet%moistavg(nland)+ veg%froot(nland,ns) &
               * MIN(soil%sfc(nland),casamet%moist(nland,ns))
          casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
               * (MIN(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland)) &
               /(soil%sfc(nland)-soil%swilt(nland))

          ! Ticket#121

          casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
               * (MAX(MIN(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland),0.0)) &
               /(soil%sfc(nland)-soil%swilt(nland))

       ENDDO
    ENDDO

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
    INTEGER i,j,k,kk,iv,thepoint,nland
    REAL(r_2), DIMENSION(mp)         :: xFluxNlittermin
    REAL(r_2), DIMENSION(mp)         :: xFluxNsoilmin
    REAL(r_2), DIMENSION(mp)         :: xFluxNsoilimm
    REAL(r_2), DIMENSION(mp)         :: xFluxNsoilminnet
    ! A maximum Clitter set to avoid positive feedback for litter accumulation
    ! when N mode is activated. (Q.Zhang 23/05/2011)
    !  real(r_2), dimension(17)         :: xClitter
    !  data xClitter/100.0,100.0,100.0,100.0,50.0,150.0,150.0,100.0,&
    !                150.0,150.0,100.0, 20.0,20.0, 20.0, 20.0, 20.0,20.0/

    xkNlimiting  = 1.0
    !  set N mineral N fluxes to zero
    xFluxNlittermin(:)  = 0.0
    xFluxNsoilmin(:)    = 0.0
    xFluxNsoilimm(:)    = 0.0  !negative for microbial upatek and postive for release of mineral N
    xFluxNsoilminnet(:) = 0.0
    !   PRINT *, 'within casa_xkN'

    !  calculate gross mineralisation
    DO nland=1,mp
       IF (casamet%iveg2(nland)/=icewater) THEN

          ! calculate C:N ratio of newly formed SOM as function of soil mineral N pool
          IF (casapool%Nsoilmin(nland) < 2.0) THEN
             casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmin(nland,:)  &
                  + (casapool%ratioNCsoilmax(nland,:) &
                  -casapool%ratioNCsoilmin(nland,:)) &
                  * MAX(0.0,casapool%Nsoilmin(nland)) / 2.0
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
                IF(k.NE.kk) THEN
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
    !    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0)
    !      xkNlimiting(:) =1.0
    !    ELSEWHERE
    !      xkNlimiting(:) =max(0.0, - (casapool%Nsoilmin(:)-2.0)/(deltpool*xFluxNsoilminnet(:)))
    !      xkNlimiting(:) =MIN(1.0,xkNlimiting(:))
    !    ENDWHERE
    ! ENDWHERE

    ! Q.Zhang 23/05/2011 test code according to YPW
    WHERE(casamet%iveg2(:)/=icewater)
       WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0 &
            .OR. xFluxNsoilminnet(:) .GE. 0.0)
          xkNlimiting(:) =1.0
       ELSEWHERE
          xkNlimiting(:) =MAX(0.0, - (casapool%Nsoilmin(:)-0.5) &
               /(deltpool*xFluxNsoilminnet(:)))
          xkNlimiting(:) =MIN(1.0,xkNlimiting(:))
       ENDWHERE
       ! Q.Zhang 23/05/2011 test
       ! If pool size larger than xClitter, turnover rate will not constrained by Nsoilmin.
       !    where(casapool%clitter(:,1) > xClitter(veg%iveg(:)))
       !     xkNlimiting(:) = 1.0
       !    end where
       ! end (Q.Zhang 23/05/2011)
       WHERE(SUM(casapool%clitter,2) > casabiome%maxfinelitter(veg%iveg(:)) + casabiome%maxcwd(veg%iveg(:)))
          xkNlimiting(:) = 1.0
       END WHERE
    ENDWHERE


  END SUBROUTINE casa_xkN

  SUBROUTINE casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
    ! (1) compute (1)N uptake by plants;
    ! (2) allocation of uptaken N to plants
    !
    IMPLICIT NONE
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    REAL(r_2), DIMENSION(mp),     INTENT(IN)    :: xkNlimiting

    ! local variables
    INTEGER                   nland,np,ip
    REAL(r_2), DIMENSION(mp,mplant)      :: Nreqmax, Nreqmin, NtransPtoP, xnuptake
    REAL(r_2), DIMENSION(mp)             :: totNreqmax,totNreqmin
    REAL(r_2), DIMENSION(mp)             :: xnCnpp

    Nreqmin(:,:)       = 0.0
    Nreqmax(:,:)       = 0.0
    NtransPtoP(:,:)    = 0.0
    totNreqmax = 0.0
    totNreqmin = 0.0

    casaflux%Nminuptake(:)     = 0.0
    casaflux%fracNalloc(:,:)   = 0.0
    xnCnpp = MAX(0.0_r_2,casaflux%Cnpp)
    CALL casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
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

          casaflux%Nminuptake(np) = xnuptake(np,leaf) + xnuptake(np,wood) + xnuptake(np,froot)+1.0e-10
          casaflux%fracNalloc(np,leaf)  = xnuptake(np,leaf)/casaflux%Nminuptake(np)
          casaflux%fracNalloc(np,wood)  = xnuptake(np,wood)/casaflux%Nminuptake(np)
          casaflux%fracNalloc(np,froot) = xnuptake(np,froot)/casaflux%Nminuptake(np)
       ENDIF
    ENDDO

    !  np=1
    !  write(*,911) casapool%nsoilmin(np),casaflux%Nminuptake(np),xnuptake(np,leaf), xnuptake(np,wood), xnuptake(np,froot), &
    !               casaflux%fracNalloc(np,leaf),casaflux%fracNalloc(np,wood),casaflux%fracNalloc(np,froot)

    casaflux%Nupland = casaflux%Nminuptake

911 FORMAT('N uptake:',100(f8.3,2x))
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

    Nreqmin(:,:)    = 0.0
    Nreqmax(:,:)    = 0.0
    NtransPtoP(:,:) = 0.0

    DO np=1,mp
       IF(casamet%iveg2(np)/=icewater) THEN
          IF(casapool%Nsoilmin(np)<2.0) THEN
             ncplantmax(np,leaf) =casabiome%ratioNCplantmin(veg%iveg(np),leaf)  &
                  +(casabiome%ratioNCplantmax(veg%iveg(np),leaf)-casabiome%ratioNCplantmin(veg%iveg(np),leaf)) &
                  * MIN(1.0,MAX(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
             ncplantmax(np,wood) =casabiome%ratioNCplantmin(veg%iveg(np),wood)  &
                  +(casabiome%ratioNCplantmax(veg%iveg(np),wood)-casabiome%ratioNCplantmin(veg%iveg(np),wood)) &
                  * MIN(1.0,MAX(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
             ncplantmax(np,froot) =casabiome%ratioNCplantmin(veg%iveg(np),froot)  &
                  +(casabiome%ratioNCplantmax(veg%iveg(np),froot)-casabiome%ratioNCplantmin(veg%iveg(np),froot)) &
                  * MIN(1.0,MAX(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
          ELSE
             ncplantmax(np,leaf)  = casabiome%ratioNCplantmax(veg%iveg(np),leaf)
             ncplantmax(np,wood)  = casabiome%ratioNCplantmax(veg%iveg(np),wood)
             ncplantmax(np,froot) = casabiome%ratioNCplantmax(veg%iveg(np),froot)
          ENDIF

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
               * (1.0-casabiome%ftransNPtoL(veg%iveg(np),leaf))
          NtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Nplant(np,wood) &
               * (1.0-casabiome%ftransNPtoL(veg%iveg(np),wood))
          NtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Nplant(np,froot) &
               * (1.0-casabiome%ftransNPtoL(veg%iveg(np),froot))

          Nreqmax(np,leaf)  = MAX(0.0,Nreqmax(np,leaf) - NtransPtoP(np,leaf))
          Nreqmax(np,wood)  = MAX(0.0,Nreqmax(np,wood) - NtransPtoP(np,wood))
          Nreqmax(np,froot) = MAX(0.0,Nreqmax(np,froot) - NtransPtoP(np,froot))
          Nreqmin(np,leaf)  = MAX(0.0,Nreqmin(np,leaf) - NtransPtoP(np,leaf))
          Nreqmin(np,wood)  = MAX(0.0,Nreqmin(np,wood) - NtransPtoP(np,wood))
          Nreqmin(np,froot) = MAX(0.0,Nreqmin(np,froot) - NtransPtoP(np,froot))

          IF(casapool%nplant(np,leaf)/(casapool%cplant(np,leaf)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),leaf)) THEN
             Nreqmax(np,leaf) = 0.0
             Nreqmin(np,leaf) =0.0
          ENDIF
          IF(casapool%nplant(np,wood)/(casapool%cplant(np,wood)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),wood)) THEN
             Nreqmax(np,wood) = 0.0
             Nreqmin(np,wood) =0.0
          ENDIF
          IF(casapool%nplant(np,froot)/(casapool%cplant(np,froot)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),froot)) THEN
             Nreqmax(np,froot) = 0.0
             Nreqmin(np,froot) =0.0
          ENDIF

       ENDIF
    ENDDO

  END SUBROUTINE casa_Nrequire

  SUBROUTINE casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
    ! (1) compute  P uptake by plants;
    ! (2) allocation of uptaken P to plants
    !
    IMPLICIT NONE
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    REAL(r_2), DIMENSION(mp),     INTENT(IN)    :: xkNlimiting


    ! local variables
    INTEGER                   nland,np,ip
    REAL(r_2), DIMENSION(mp,mplant) :: Preqmax,Preqmin,PtransPtoP,xPuptake
    REAL(r_2), DIMENSION(mp)        :: totPreqmax,totPreqmin
    REAL(r_2), DIMENSION(mp)        :: xpCnpp

    Preqmin(:,:)             = 0.0
    Preqmax(:,:)             = 0.0
    PtransPtoP(:,:)          = 0.0
    casaflux%Plabuptake(:)   = 0.0
    casaflux%fracPalloc(:,:) = 0.0
    totPreqmax               = 0.0
    totPreqmin               = 0.0

    xpCnpp = MAX(0.0_r_2,casaflux%cnpp)
    CALL casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
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

       casaflux%Plabuptake(:) = xpuptake(:,leaf) + xpuptake(:,wood) + xpuptake(:,froot)+1.0e-10
       casaflux%fracPalloc(:,leaf)  = xpuptake(:,leaf)/casaflux%Plabuptake(:)
       casaflux%fracPalloc(:,wood)  = xpuptake(:,wood)/casaflux%Plabuptake(:)
       casaflux%fracPalloc(:,froot) = xpuptake(:,froot)/casaflux%Plabuptake(:)

    ENDWHERE

    casaflux%Pupland = casaflux%Plabuptake

    !  np=1
    !  write(*,911) casapool%Psoillab(np),casaflux%Plabuptake(np), &
    !               xpuptake(np,leaf), xpuptake(np,wood), xpuptake(np,froot), &
    !               casaflux%fracPalloc(np,leaf),casaflux%fracPalloc(np,wood), &
    !               casaflux%fracPalloc(np,froot)
    !911 format('P uptake:',100(f8.3,2x))

    !  ! only used in spinning up the model
    !  DO  np=1,mp
    !    casaflux%Plabuptake(np) = TotPreqmax(np)
    !    casaflux%Pupland(np)    = TotPreqmax(np)
    !    casaflux%Pwea(np)       = TotPreqmax(np)
    !  ENDDO

  END SUBROUTINE casa_puptake

  SUBROUTINE casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
       casabiome,casapool,casaflux,casamet)
    IMPLICIT NONE
    REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xpCnpp
    REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Preqmax, Preqmin
    REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: PtransPtoP
    TYPE (veg_parameter_type),             INTENT(INOUT) :: veg
    TYPE (casa_biome),                     INTENT(INOUT) :: casabiome
    TYPE (casa_pool),                      INTENT(INOUT) :: casapool
    TYPE (casa_flux),                      INTENT(INOUT) :: casaflux
    TYPE (casa_met),                       INTENT(INOUT) :: casamet

    ! local variables
    INTEGER :: nland,np,ip

    Preqmin(:,:)       = 0.0
    Preqmax(:,:)       = 0.0
    PtransPtoP(:,:)    = 0.0
    DO np=1,mp
       IF(casamet%iveg2(np)/=icewater) THEN
         Preqmax(np,leaf) = xpCnpp(np)* casaflux%fracCalloc(np,leaf) &
                          * casabiome%ratioPCplantmax(veg%iveg(np),leaf)
         Preqmax(np,wood) = xpCnpp(np)* casaflux%fracCalloc(np,wood) &
                          * casabiome%ratioPCplantmax(veg%iveg(np),wood)
         Preqmax(np,froot) = xpCnpp(np)* casaflux%fracCalloc(np,froot) &
                           * casabiome%ratioPCplantmax(veg%iveg(np),froot)

         Preqmin(np,leaf) = xpCnpp(np) * casaflux%fracCalloc(np,leaf) &
                         * casabiome%ratioPCplantmin(veg%iveg(np),leaf)
         Preqmin(np,wood) = xpCnpp(np) * casaflux%fracCalloc(np,wood) &
                         * casabiome%ratioPCplantmin(veg%iveg(np),wood)
         Preqmin(np,froot) = xpCnpp(np) * casaflux%fracCalloc(np,froot) &
                          * casabiome%ratioPCplantmin(veg%iveg(np),froot)

         PtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Pplant(np,leaf) &
                            * (1.0-casabiome%ftransPPtoL(veg%iveg(np),leaf))
         PtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Pplant(np,wood) &
                            * (1.0-casabiome%ftransPPtoL(veg%iveg(np),wood))
         PtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Pplant(np,froot) &
                            * (1.0-casabiome%ftransPPtoL(veg%iveg(np),froot))

         Preqmax(np,leaf)    = max(0.0,Preqmax(np,leaf) - PtransPtoP(np,leaf))
         Preqmax(np,wood)    = max(0.0,Preqmax(np,wood) - PtransPtoP(np,wood))
         Preqmax(np,froot)   = max(0.0,Preqmax(np,froot)- PtransPtoP(np,froot))

         Preqmin(np,leaf)    = max(0.0,Preqmin(np,leaf) - PtransPtoP(np,leaf))
         Preqmin(np,wood)    = max(0.0,Preqmin(np,wood) - PtransPtoP(np,wood))
         Preqmin(np,froot)   = max(0.0,Preqmin(np,froot)- PtransPtoP(np,froot))

          IF(casapool%pplant(np,leaf)/(casapool%cplant(np,leaf)+1.0e-10)> casabiome%ratioPCplantmax(veg%iveg(np),leaf)) THEN
             Preqmax(np,leaf) = 0.0
             Preqmin(np,leaf) =0.0
          ENDIF
          IF(casapool%pplant(np,wood)/(casapool%cplant(np,wood)+1.0e-10)> casabiome%ratioPCplantmax(veg%iveg(np),wood)) THEN
             Preqmax(np,wood) = 0.0
             Preqmin(np,wood) =0.0
          ENDIF
          IF(casapool%pplant(np,froot)/(casapool%cplant(np,froot)+1.0e-10)> casabiome%ratioPCplantmin(veg%iveg(np),froot)) THEN
             Preqmax(np,froot) = 0.0
             Preqmin(np,froot) =0.0
          ENDIF

       ENDIF
    ENDDO

  END SUBROUTINE casa_Prequire


  SUBROUTINE casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet,casabal, LALLOC)
    ! update all pool sizes
    !
    IMPLICIT NONE
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal    
    INTEGER , INTENT(IN) :: LALLOC
    ! local variables
    REAL(r_2), DIMENSION(mp)   :: plabsorb,deltap
    REAL(r_2), DIMENSION(mp,mwood)   :: delcwoodprod,delnwoodprod,delpwoodprod
    INTEGER i,j,k,np,nland

    !  PRINT *, 'mp here is ', mp
    !  print*,'cplant',casapool%cplant(:,leaf)

    ! store the pool size before updating
    !  ypwang 15-6-2021
    casabal%cplantlast  = casapool%cplant
    casabal%clabilelast = casapool%clabile
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil
    if(icycle >1) then
       casabal%nplantlast  = casapool%nplant
       casabal%nlitterlast = casapool%nlitter
       casabal%nsoillast   = casapool%nsoil
       casabal%nsoilminlast= casapool%nsoilmin
       if(icycle >2) then
          casabal%pplantlast   = casapool%pplant
          casabal%plitterlast  = casapool%plitter
          casabal%psoillast    = casapool%psoil
          casabal%psoillablast = casapool%psoillab
          casabal%psoilsorblast= casapool%psoilsorb
          casabal%psoilocclast = casapool%psoilocc
       endif        
    endif
                 
    
    DO np=1,mp
       IF(casamet%iveg2(np) == icewater) THEN
          casamet%glai(np)   = 0.0
       ELSE

          !if (np==2) write(912,91) casapool%cplant(np,:),  casapool%dcplantdt(np,:)  * deltpool
91        FORMAT (100(e12.4,2x))
          casapool%cplant(np,:)  = casapool%cplant(np,:)  &
               + casapool%dcplantdt(np,:)  * deltpool
          casapool%clabile(np)   = casapool%clabile(np)   &
               + casapool%dclabiledt(np)   * deltpool



          IF(casapool%cplant(np,leaf) > 0.0) THEN
             IF(icycle >1) casapool%Nplant(np,:) = casapool%Nplant(np,:) &
                  +casapool%dNplantdt(np,:)*deltpool
             IF(icycle >2) casapool%Pplant(np,:) = casapool%Pplant(np,:) &
                  +casapool%dPplantdt(np,:)*deltpool
          ENDIF
          !    casamet%glai(np)   = MIN(0.0, casabiome%sla(veg%iveg(np))  &
          !                                  * casapool%cplant(np,leaf))
          casamet%glai(np)   = MAX(casabiome%glaimin(veg%iveg(np)), &
               casabiome%sla(veg%iveg(np)) * casapool%cplant(np,leaf))
          ! vh !
          IF (LALLOC.NE.3) THEN
             casamet%glai(np)   = MIN(casabiome%glaimax(veg%iveg(np)), casamet%glai(np))
          ENDIF
          casapool%clitter(np,:) = casapool%clitter(np,:) &
               + casapool%dClitterdt(np,:) * deltpool
          casapool%csoil(np,:)   = casapool%csoil(np,:)   &
               + casapool%dCsoildt(np,:)   * deltpool

          IF(icycle >1) THEN
             casapool%Nlitter(np,:) = casapool%Nlitter(np,:) &
                  + casapool%dNlitterdt(np,:)* deltpool
             casapool%Nsoil(np,:)   = casapool%Nsoil(np,:)   &
                  + casapool%dNsoildt(np,:)  * deltpool
             ! vh ! put lower bound of 1.e-3 to prevent Nsoilmin from going negative
             ! Ticket #108
             casapool%Nsoilmin(np)  = MAX(casapool%Nsoilmin(np)  &
                  + casapool%dNsoilmindt(np) * deltpool,1.e-3)
          ENDIF

          IF(icycle >2) THEN
             casapool%Plitter(np,:) = casapool%Plitter(np,:) &
                  + casapool%dPlitterdt(np,:)* deltpool
             casapool%Psoil(np,:)   = casapool%Psoil(np,:)   &
                  + casapool%dPsoildt(np,:)  * deltpool
             casapool%Psoillab(np)  = casapool%Psoillab(np)  &
                  + casapool%dPsoillabdt(np) * deltpool
             casapool%Psoilsorb(np) = casaflux%Psorbmax(np)*casapool%Psoillab(np) &
                  /(casaflux%kmlabp(np)+casapool%Psoillab(np))
             !      casapool%Psoilsorb(np) = casapool%Psoilsorb(np)  &
             !                             + casapool%dPsoilsorbdt(np) * deltpool
             casapool%Psoilocc(np)   = casapool%Psoilocc(np)  &
                  + casapool%dPsoiloccdt(np)  * deltpool
          ENDIF

          DO i=1,mplant
             IF(casapool%cplant(np,i) < 0.0)  THEN
                WRITE(57,*)  'Cpool: np,ivt',np,casamet%lat(np),casamet%lon(np), &
                     casamet%iveg2(np),casapool%cplant(np,:)
                CALL casa_poolzero(np,1,casapool)
                !stop
                casapool%cplant(np,i) = MAX(0.0, casapool%cplant(np,i))
             ENDIF
          ENDDO
          IF(icycle >1) THEN
             DO i=1,mplant
                IF(casapool%nplant(np,i) < 0.0) THEN
                   WRITE(57,*) 'Npool:', 'np,ivt,ipool',np,casamet%iveg2(np),casapool%nplant(np,:)
                   CALL casa_poolzero(np,2,casapool)
                   casapool%nplant(np,i) = MAX(0.0, casapool%nplant(np,i))
                ENDIF
             ENDDO
          ENDIF ! end of "icycle >1"

          DO j=1,mlitter
             IF(casapool%clitter(np,j) < 0.0)  THEN
                WRITE(57,*)  'Clitter: np,ivt2',np,casamet%iveg2(np),casapool%clitter(np,:)
                CALL casa_poolzero(np,3,casapool)
                casapool%clitter(np,j) = MAX(0.0, casapool%clitter(np,j))
             ENDIF
          ENDDO

          DO k=1,msoil
             IF(casapool%csoil(np,k) < 0.0)    THEN
                WRITE(57,*)  'Csoil: np,ivt2',np,casamet%iveg2(np),casapool%csoil(np,:)
                CALL casa_poolzero(np,5,casapool)
                casapool%csoil(np,k) = MAX(0.0, casapool%csoil(np,k))
             ENDIF
          ENDDO

          !  check if any pool size, and terminate model run if any pool size is negative!!
          IF(icycle >1) THEN
             DO j=1,mlitter
                IF(casapool%nlitter(np,j) < 0.0)  THEN
                   WRITE(57,*)  'Nlitter: np,ivt2',np,casamet%iveg2(np),casapool%Nlitter(np,:)
                   CALL casa_poolzero(np,4,casapool)
                   casapool%nlitter(np,j) = MAX(0.0, casapool%nlitter(np,j))
                ENDIF
             ENDDO
             DO k=1,msoil
                IF(casapool%nsoil(np,k) < 0.0) THEN
                   WRITE(57,*)  'Nsoil: np,ivt2',np,casamet%iveg2(np),casapool%nsoil(np,:)
                   CALL casa_poolzero(np,6,casapool)
                   casapool%nsoil(np,k) = MAX(0.0, casapool%nsoil(np,k))
                ENDIF
             ENDDO
          ENDIF  !end of "icycle >1"
       ENDIF
    ENDDO !end of "np"

    ! calculate the decay of wood product pools
    if(l_landuse) then

       DO np=1,mp
          delcwoodprod(np,:) =0.0        
          delnwoodprod(np,:) =0.0        
          delpwoodprod(np,:) =0.0        
          IF(casamet%iveg2(np) .ne. icewater) THEN
             delcwoodprod(np,:) = - ratewoodprod(:) *casapool%cwoodprod(np,:) * deltcasa
             delnwoodprod(np,:) = - ratewoodprod(:) *casapool%nwoodprod(np,:) * deltcasa
             delpwoodprod(np,:) = - ratewoodprod(:) *casapool%pwoodprod(np,:) * deltcasa
             casapool%cwoodprod(np,:) = casapool%cwoodprod(np,:) + delcwoodprod(np,:)
             casapool%nwoodprod(np,:) = casapool%nwoodprod(np,:) + delnwoodprod(np,:)
             casapool%pwoodprod(np,:) = casapool%pwoodprod(np,:) + delpwoodprod(np,:)
          endif        
       ENDDO   

     endif


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

  SUBROUTINE casa_cnpbal(veg,casamet,casapool,casaflux,casabal)

    IMPLICIT NONE
    TYPE (casa_met),              INTENT(IN)    :: casamet
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_balance),          INTENT(INOUT) :: casabal

    TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters

    ! local variables
    INTEGER :: npt
    REAL(r_2), DIMENSION(mp) :: cbalplant,  nbalplant,  pbalplant
    REAL(r_2), DIMENSION(mp) :: cbalsoil,   nbalsoil,   pbalsoil
    REAL(r_2), DIMENSION(mp) :: cbalplantx, nbalplantx, pbalplantx
    ! temporary variables
    real latx,lonx
    integer ivegx, isox

    cbalplant(:) = 0.0
    cbalsoil(:)  = 0.0
    nbalplant(:) = 0.0
    nbalsoil(:)  = 0.0
    pbalplant(:) = 0.0
    pbalsoil(:)  = 0.0

    casabal%cbalance(:)  = 0.0
    casabal%nbalance(:)  = 0.0
    casabal%pbalance(:)  = 0.0

    !C balance
    Cbalplant(:)  = SUM(casabal%cplantlast,2) -SUM(casapool%cplant,2)            &
         + casabal%Clabilelast(:)-casapool%clabile(:)                   &
         +(casaflux%Cnpp(:) - SUM((casaflux%kplant*casabal%cplantlast),2))*deltpool &
         + casapool%dClabiledt(:)* deltpool
    Cbalsoil(:)   = SUM(casabal%clitterlast,2) - SUM(casapool%clitter,2)         &
         + SUM(casabal%csoillast,2)   - SUM(casapool%csoil,2)           &
         +(SUM((casaflux%kplant*casabal%cplantlast),2)-casaflux%Crsoil(:))*deltpool

    casabal%cbalance(:) = Cbalplant(:) + Cbalsoil(:)


    DO npt=1,mp
       IF(ABS(casabal%cbalance(npt))>1e-10.and.casapool%cplant(npt,1) >1.0e-3) THEN
       
          WRITE(*,*) 'cbalance',   npt, casamet%lat(npt),casamet%lon(npt), &
                                   casamet%iveg2(npt),Cbalplant(npt), Cbalsoil(npt)
          WRITE(*,*) 'gpp, npp',   casaflux%Cgpp(npt), casaflux%Cnpp(npt)
          WRITE(*,*) 'rmplant, rgplant',  casaflux%crmplant(npt,:) , casaflux%crgplant(npt)
          WRITE(*,*) 'DIFF cplant',  SUM(casapool%cplant(npt,:)) -SUM(casabal%cplantlast(npt,:))
          WRITE(*,*) 'dcplandt',   casapool%dcplantdt(npt,:), SUM(casapool%dcplantdt(npt,:))
          WRITE(*,*) 'litterfall', casaflux%kplant(npt,:) * casabal%cplantlast(npt,:)*deltpool, &
                                   SUM((casaflux%kplant(npt,:)*casabal%cplantlast(npt,:)))*deltpool
          write(*,*) 'deplantdt-1', casapool%dcplantdt(npt,1), casapool%cplant(npt,1)-casabal%cplantlast(npt,1), &
                      casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,1) - casaflux%kplant(npt,1)  * casabal%cplantlast(npt,1)
          write(*,*) 'deplantdt-2', casapool%dcplantdt(npt,2), casapool%cplant(npt,2)-casabal%cplantlast(npt,2),  &
                      casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,2) - casaflux%kplant(npt,2)  * casabal%cplantlast(npt,2)
          write(*,*) 'deplantdt-3', casapool%dcplantdt(npt,3),   casapool%cplant(npt,3)-casabal%cplantlast(npt,3),&
                      casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,3) - casaflux%kplant(npt,3)  * casabal%cplantlast(npt,3)

          WRITE(*,*) 'delclabile', casabal%Clabilelast(npt)-casapool%clabile(npt)
          WRITE(*,*), 'dclabile',  casapool%dClabiledt(npt)* deltpool

       ENDIF
    ENDDO




    casapool%ctot_0 = SUM(casabal%cplantlast,2)+SUM(casabal%clitterlast,2) &
         + SUM(casabal%csoillast,2)+ casabal%clabilelast
    casapool%ctot = SUM(casapool%cplant,2)+SUM(casapool%clitter,2) &
         + SUM(casapool%csoil,2)+ casapool%clabile
    casabal%sumcbal     = casabal%sumcbal + casabal%cbalance


    IF(icycle >1) THEN
       Nbalplant(:) = SUM(casabal%nplantlast,2) -SUM(casapool%nplant,2)                  &
            +casaflux%Nminuptake(:) *deltpool
       Nbalsoil(:)  = -SUM(casapool%nlitter,2)-SUM(casapool%nsoil,2)                     &
            -casapool%nsoilmin(:)+ casabal%nsoilminlast(:)                     &
            + SUM(casabal%nlitterlast,2)    + SUM(casabal%nsoillast,2)         &
            +(casaflux%Nmindep(:) + casaflux%Nminfix(:)- casaflux%Nminloss(:)                        &
            -casaflux%Nminleach(:)-casaflux%Nupland(:)) * deltpool

       casabal%nbalance(:) = Nbalplant(:) + Nbalsoil(:)

       casabal%sumnbal     = casabal%sumnbal + casabal%nbalance

    ENDIF

    IF(icycle >2) THEN
       Pbalplant(:) = SUM(casabal%Pplantlast,2) -SUM(casapool%Pplant,2)                      &
            + casaflux%Plabuptake(:) *deltpool
       Pbalsoil(:)  = -SUM(casapool%Plitter,2)        - SUM(casapool%Psoil,2)                      &
            + SUM(casabal%Plitterlast,2)    + SUM(casabal%Psoillast,2)                   &
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





91  FORMAT('balance= ',100(f12.5,2x))

  END SUBROUTINE casa_cnpbal

  SUBROUTINE casa_ndummy(casamet,casabal,casapool)
    IMPLICIT NONE
    TYPE (casa_met),              INTENT(IN)    :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal
    TYPE (casa_pool),             INTENT(INOUT) :: casapool

        casapool%Nplant(:,:) = casapool%ratioNCplant(:,:)  * casapool%cplant(:,:)
        casapool%Nlitter(:,:)= casapool%ratioNClitter(:,:) * casapool%clitter(:,:)
        casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:)   * casapool%Csoil(:,:)
        casapool%nsoilmin(:) = 2.0
        casabal%sumnbal(:)   = 0.0
        WHERE(casamet%iveg2==grass)
              casapool%nplant(:,wood) = 0.0
              casapool%nlitter(:,cwd) = 0.0
        ENDWHERE

   END SUBROUTINE casa_ndummy

  SUBROUTINE casa_pdummy(casamet,casabal,casaflux,casapool)
    IMPLICIT NONE
    TYPE (casa_met),              INTENT(IN)    :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal
    TYPE (casa_flux),             INTENT(IN)    :: casaflux
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    ! ypw: the following data block should be consistent with "casa_poolout"
    ! local variables
    REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
    REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
    REAL(r_2), DIMENSION(mp)  :: totpsoil
    INTEGER  npt,nout,nso

    ! Soiltype     soilnumber soil P(g P/m2)
    ! Alfisol     1       61.3
    ! Andisol     2       103.9
    ! Aridisol    3       92.8
    ! Entisol     4       136.9
    ! Gellisol    5       98.2
    ! Histosol    6       107.6
    ! Inceptisol  7       84.1
    ! Mollisol    8       110.1
    ! Oxisol      9       35.4
    ! Spodosol    10      41.0
    ! Ultisol     11      51.5
    ! Vertisol    12      190.6

    DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
    DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
    DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
    DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
    DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
    DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
    DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/


    totpsoil(:) = psorder(casamet%isorder(:)) * xpsoil50(casamet%isorder(:))
    casabal%sumpbal(:)    = 0.0
    casapool%pplant(:,:)  = casapool%Cplant(:,:)  * casapool%ratioPcplant(:,:)
    casapool%plitter(:,:) = casapool%Clitter(:,:) * casapool%ratioPclitter(:,:)
    casapool%psoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioPcsoil(:,:)

    WHERE(casamet%iveg2==grass)
          casapool%pplant(:,wood) = 0.0
          casapool%plitter(:,cwd)  = 0.0
    ENDWHERE

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

    WHERE(veg%iveg==1 .OR. veg%iveg ==2 )
       phen%phase = 2
    ENDWHERE

  END SUBROUTINE phenology

  REAL FUNCTION vcmax_np(nleaf, pleaf)
    IMPLICIT NONE
    REAL, INTENT(IN) :: nleaf ! leaf N in g N m-2 leaf
    REAL, INTENT(IN) :: pleaf ! leaf P in g P m-2 leaf

    !Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax –
    !to leaf nitrogen, leaf phosphorus, and specific leaf area:
    !a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
    vcmax_np = EXP(3.946 + 0.921*LOG(nleaf) + 0.121*LOG(pleaf) + &
         0.282*LOG(pleaf)*LOG(nleaf)) * 1.0e-6 ! units of mol m-2 (leaf)


  END FUNCTION vcmax_np

END MODULE casa_cnp_module
