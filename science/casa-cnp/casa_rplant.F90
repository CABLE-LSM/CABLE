!#define ESM15 YES
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
!   casa_rplant

MODULE casa_rplant_module
  !jhan:move thesse to subr AND ONLY-ise
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module, ONLY: cable_user,l_landuse ! Custom soil respiration: Ticket #42
#ifndef ESM15 
  USE landuse_constant
#endif
  IMPLICIT NONE

CONTAINS

SUBROUTINE casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)
! maintenance respiration of woody tisse and fineroots
! see Sitch et al. (2003), GCB, reqn (23)

USE casa_cnp_module, ONLY : vcmax_np 
    IMPLICIT NONE
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (casa_biome),          INTENT(INOUT) :: casabiome
    TYPE (casa_pool),           INTENT(INOUT) :: casapool
    TYPE (casa_flux),           INTENT(INOUT) :: casaflux
    TYPE (casa_met),            INTENT(INOUT) :: casamet
    TYPE (climate_type),            INTENT(IN) :: climate
    INTEGER :: npt, ivt

    REAL(r_2), DIMENSION(mp)        :: Ygrow        ! growth efficiency Q.Zhang 22/02/2011
    REAL(r_2), DIMENSION(mp,mplant) :: ratioPNplant ! Q.Zhang 22/02/2011
    REAL(r_2), DIMENSION(mp)        :: delcrmleaf, delcrmwood,delcrmfroot    ! reduction in wood and root respiration when NPP <0.0
    REAL(r_2), DIMENSION(mp)        :: resp_coeff_root, resp_coeff_sapwood, resp_coeff
    REAL,  DIMENSION(mp)        :: nleaf, pleaf, vcmaxmax

    resp_coeff = 1
    resp_coeff_root = 1
    resp_coeff_sapwood = 1
    ratioPNplant = 0.0
    Ygrow        = 0.0

    WHERE(casapool%Nplant>0.0)
       ratioPNplant = casapool%Pplant/(casapool%Nplant+ 1.0e-10)
    ENDWHERE

    Ygrow(:) = 0.65+0.2*ratioPNplant(:,leaf)/(ratioPNplant(:,leaf)+1.0/15.0)
    Ygrow(:) = min(0.85,max(0.65,Ygrow(:)))

    casaflux%crmplant(:,wood) = 0.0
    casaflux%crmplant(:,froot) = 0.0
    delcrmleaf   = 0.0
    delcrmwood   = 0.0
    delcrmfroot  = 0.0
    casaflux%crgplant = 0.0
    casaflux%clabloss = 0.0

    IF (cable_user%CALL_climate) THEN
       ! coefficients required to implement T-acclimation of autotrophic respiration (Ticket # 110)
       ! adapted from Atkin et al., New Phyt., 2015)
       DO npt = 1, mp
          ivt=veg%iveg(npt)
          ! max leaf N in g N m-2 leaf
          nleaf(npt) =  casabiome%ratioNCplantmax(ivt,leaf)/casabiome%sla(ivt)
          ! max leaf P in g P m-2 leaf
          pleaf(npt) = casabiome%ratioPcplantmax(ivt,leaf)/casabiome%sla(ivt)
          IF (ivt .EQ. 7) THEN
             ! special for C4 grass: set here to value from  parameter file
             vcmaxmax(npt) = 1.0e-5
          ELSE
             vcmaxmax(npt) = vcmax_np(nleaf(npt), pleaf(npt))
          ENDIF
          IF (veg%iveg(npt).EQ.2 .OR. veg%iveg(npt).EQ. 4  ) THEN
             ! broadleaf forest

             resp_coeff_root(npt) = (1.2818 * 1.e-6 *casapool%nplant(npt,froot)/ &
                  vcmaxmax(npt)/0.0116   + &
                  casapool%nplant(npt,froot)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                  casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )

             resp_coeff_sapwood(npt) = (1.2818 * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                  casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      )


          ELSEIF (veg%iveg(npt).EQ.1 .OR. veg%iveg(npt).EQ. 3  ) THEN
             ! needleleaf forest

             resp_coeff_root(npt) = (1.2877 * 1.e-6 *casapool%nplant(npt,froot) &
                  /vcmaxmax(npt)/0.0116   + &
                  casapool%nplant(npt,froot)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                  casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )

             resp_coeff_sapwood(npt) = (1.2877 * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                  casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      )



          ELSEIF (veg%iveg(npt).EQ.6 .OR. veg%iveg(npt).EQ.8 .OR. veg%iveg(npt).EQ. 9  ) THEN
             ! C3 grass, tundra, crop

             resp_coeff_root(npt) = (1.6737 * 1.e-6 *casapool%nplant(npt,froot)/ &
                  vcmaxmax(npt)/0.0116   + &
                  casapool%nplant(npt,froot)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                  casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )

             resp_coeff_sapwood(npt) = (1.6737 * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                  casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      )
          ELSE
             ! shrubs and other (C4 grass and crop)
             resp_coeff_root(npt) = (1.5758 * 1.e-6 *casapool%nplant(npt,froot)/ &
                  vcmaxmax(npt)/0.0116   + &
                  casapool%nplant(npt,froot)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                  casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )

             resp_coeff_sapwood(npt) = (1.5758 * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                  casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                  0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                  casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      )
          ENDIF
       ENDDO
       resp_coeff = 0.50
    ENDIF  ! end coefficients for acclimation of autotrophic respiration Ticket #110


    IF (cable_user%CALL_climate) THEN
       !  acclimation of autotrophic respiration Ticket #110
       WHERE(casamet%iveg2/=icewater)
          WHERE(casamet%tairk >250.0)
             WHERE(casapool%cplant(:,wood)>1.0e-6)
                casaflux%crmplant(:,wood)  =  resp_coeff  * resp_coeff_sapwood * &
                     casabiome%rmplant(veg%iveg(:),wood) &
                     * EXP(308.56*(1.0/56.02-1.0           &
                     / (casamet%tairk(:)+46.02-tkzeroc)))

             ENDWHERE
             !vh! prevent floating underflow with this mask
             WHERE (casapool%Clabile(:).GT.1.e-8) &
                  casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                  * MAX(0.0,casapool%Clabile(:))      &
                  * EXP(308.56*(1.0/56.02-1.0         &
                  / (casamet%tairk(:)+46.02-tkzeroc)))


          ENDWHERE

          WHERE(casamet%tsoilavg >250.0.AND.casapool%cplant(:,froot)>1.0e-6)

             casaflux%crmplant(:,froot) =  resp_coeff * resp_coeff_root * &
                  casabiome%rmplant(veg%iveg(:),froot) &
                  * EXP(308.56*(1.0/56.02-1.0            &
                  / (casamet%tsoilavg(:)+46.02-tkzeroc)))

          ENDWHERE

          WHERE((casaflux%Cgpp-SUM(casaflux%crmplant,2))>0.0)
             !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
             ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
             casaflux%crgplant(:)  = (1.0-Ygrow(:))* &
                  MAX(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))

          ELSEWHERE
             casaflux%crgplant(:) = 0.0
          ENDWHERE
       ENDWHERE !(casamet%iveg2/=icewater)

       Casaflux%cnpp(:) = casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2) - casaflux%crgplant(:)

    ELSE !IF (cable_user%CALL_climate) THEN


       WHERE(casamet%iveg2/=icewater)
          WHERE(casamet%tairk >250.0)
             WHERE(casapool%cplant(:,wood)>1.0e-6)

#             ifdef ESM15 
                casaflux%crmplant(:,wood)  = casabiome%rmplant(veg%iveg(:),wood) &
#             else
                casaflux%crmplant(:,wood)  =  resp_coeff * casaflux%frac_sapwood(:) * &
                     casabiome%rmplant(veg%iveg(:),wood) &
#             endif
                     * casapool%nplant(:,wood)             &
                     * EXP(308.56*(1.0/56.02-1.0           &
                     / (casamet%tairk(:)+46.02-tkzeroc)))


             ENDWHERE
             casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                  * MAX(0.0,casapool%Clabile(:))      &
                  * EXP(308.56*(1.0/56.02-1.0         &
                  / (casamet%tairk(:)+46.02-tkzeroc)))
          ENDWHERE
          WHERE(casamet%tsoilavg >250.0.AND.casapool%cplant(:,froot)>1.0e-6)

#         ifdef ESM15 
            casaflux%crmplant(:,froot) = casabiome%rmplant(veg%iveg(:),froot) &
#         else
            casaflux%crmplant(:,froot) =  resp_coeff * casabiome%rmplant(veg%iveg(:),froot) &
#         endif
                  * casapool%nplant(:,froot)             &
                  * EXP(308.56*(1.0/56.02-1.0            &
                  / (casamet%tsoilavg(:)+46.02-tkzeroc)))

          ENDWHERE
          casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + casaflux%clabloss(:)

          WHERE((casaflux%Cgpp-SUM(casaflux%crmplant,2))>0.0)
             !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
             ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
             casaflux%crgplant(:)  = (1.0-Ygrow(:))* &
                  MAX(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))

          ELSEWHERE
             casaflux%crgplant(:) = 0.0
          ENDWHERE


          Casaflux%cnpp(:) = casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2) - casaflux%crgplant(:)

    WHERE(casaflux%Cnpp < 0.0)
! change made here by ypw on 11-7-2016 to include leaf maintenance respiration
      delcrmleaf(:)  = casaflux%Cnpp(:) * casaflux%crmplant(:,leaf) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))
      delcrmwood(:)  = casaflux%Cnpp(:) * casaflux%crmplant(:,wood) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))
      delcrmfroot(:) = casaflux%Cnpp(:) * casaflux%crmplant(:,froot) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))

         casaflux%crmplant(:,leaf)  = casaflux%crmplant(:,leaf)  + delcrmleaf(:)
      casaflux%crmplant(:,wood)  = casaflux%crmplant(:,wood)  + delcrmwood(:)
      casaflux%crmplant(:,froot) = casaflux%crmplant(:,froot) + delcrmfroot(:)
      casaflux%crgplant(:) = 0.0
    ENDWHERE



       ENDWHERE ! (casamet%iveg2/=icewater)
  
  casaflux%Cnpp(:) = casaflux%Cgpp(:) - SUM(casaflux%crmplant(:,:),2) &
                   - casaflux%crgplant(:)

    ENDIF

END SUBROUTINE casa_rplant


SUBROUTINE casa_rplant1(veg,casabiome,casapool,casaflux,casamet)
! maintenance respiration of woody tisse and fineroots
! see Sitch et al. (2003), GCB, reqn (23)

  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  INTEGER :: npt

  real(r_2), dimension(mp)        :: Ygrow        ! growth efficiency Q.Zhang 22/02/2011
  real(r_2), dimension(mp,mplant) :: ratioPNplant ! Q.Zhang 22/02/2011
  real(r_2), dimension(mp)        :: delcrmleaf, delcrmwood,delcrmfroot    ! reduction in wood and root respiration when NPP <0.0

  ratioPNplant(:,:) = 1.0/casabiome%ratioNPplantmin(veg%iveg(:),:)
  WHERE(casapool%Nplant>0.0)
    ratioPNplant = casapool%Pplant/casapool%Nplant
  ENDWHERE

  Ygrow(:) = 0.65+0.2*ratioPNplant(:,leaf)/(ratioPNplant(:,leaf)+1.0/15.0)

  casaflux%crmplant(:,wood) = 0.0
  casaflux%crmplant(:,froot) = 0.0
  delcrmleaf   = 0.0
  delcrmwood   = 0.0
  delcrmfroot  = 0.0
  casaflux%crgplant = 0.0
  casaflux%clabloss = 0.0

  WHERE(casamet%iveg2/=icewater)
    WHERE(casamet%tairk >250.0)
      WHERE(casapool%cplant(:,wood)>1.0e-6)
      casaflux%crmplant(:,wood)  = casabiome%rmplant(veg%iveg(:),wood) &
                                 * casapool%nplant(:,wood)             &
                                 * exp(308.56*(1.0/56.02-1.0           &
                                 / (casamet%tairk(:)+46.02-tkzeroc)))
      ENDWHERE
      casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                            * max(0.0,casapool%Clabile(:))      &
                            * exp(308.56*(1.0/56.02-1.0         &
                            / (casamet%tairk(:)+46.02-tkzeroc)))
    ENDWHERE
    WHERE(casamet%tsoilavg >250.0.and.casapool%cplant(:,froot)>1.0e-6)
      casaflux%crmplant(:,froot) = casabiome%rmplant(veg%iveg(:),froot) &
                                 * casapool%nplant(:,froot)             &
                                 * exp(308.56*(1.0/56.02-1.0            &
                                 / (casamet%tsoilavg(:)+46.02-tkzeroc)))
    ENDWHERE
    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + casaflux%clabloss(:)

    WHERE((casaflux%Cgpp-SUM(casaflux%crmplant,2))>0.0)
    ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
      casaflux%crgplant(:)  = (1.0-Ygrow(:))* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))

    ELSEWHERE
      casaflux%crgplant(:) = 0.0
    ENDWHERE

!!!!!!!!!!!!!!!!!!!!!!! begin from YPW 02/10/17 !!!!!!!!!!!!!!!!!!!!!!!!!!!

    casaflux%Cnpp(:) = casaflux%Cgpp(:) - SUM(casaflux%crmplant(:,:),2) &
                     - casaflux%crgplant(:)

    WHERE(casaflux%Cnpp < 0.0)
! change made here by ypw on 11-7-2016 to include leaf maintenance respiration
      delcrmleaf(:)  = casaflux%Cnpp(:) * casaflux%crmplant(:,leaf) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))
      delcrmwood(:)  = casaflux%Cnpp(:) * casaflux%crmplant(:,wood) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))
      delcrmfroot(:) = casaflux%Cnpp(:) * casaflux%crmplant(:,froot) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))

      casaflux%crmplant(:,leaf)  = casaflux%crmplant(:,leaf)  + delcrmleaf(:)
      casaflux%crmplant(:,wood)  = casaflux%crmplant(:,wood)  + delcrmwood(:)
      casaflux%crmplant(:,froot) = casaflux%crmplant(:,froot) + delcrmfroot(:)
      casaflux%crgplant(:) = 0.0

      ! The logic above can still lead to a negative NPP as
      ! SUM(casaflux%crmplant(:,:),2) can be a tiny number, if this
      ! happens, set NPP to zero
      WHERE(casaflux%Cnpp < 0.0)
         casaflux%cnpp(:) = 0.0
      ENDWHERE

    ENDWHERE

!!!!!!!!!!!!!!!!!!!!!!! end from YPW 02/10/17 !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENDWHERE

  casaflux%Cnpp(:) = casaflux%Cgpp(:) - SUM(casaflux%crmplant(:,:),2) &
                   - casaflux%crgplant(:)

END SUBROUTINE casa_rplant1

END MODULE casa_rplant_module
