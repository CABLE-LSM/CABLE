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
! Purpose: Calculates plant and soil respiration and updates carbon pools.
!          Note: carbon pools reset if ACCESS model run restarted
!          Use CASA-CNP in preference to these routines for carbon fluxes
!
! Called from: cbm
!
! Contact: Rachel.Law@csiro.au
!
! History: Plant respiration code moved from canopy
!          Namelist parameter used to switch between soil respiration options
!
!
! ==============================================================================

MODULE cable_carbon_module

   USE cable_data_module, ONLY : icarbon_type, point2constants

   IMPLICIT NONE

   PUBLIC carbon_pl, soilcarb, plantcarb
   PRIVATE

   TYPE( icarbon_type ) :: C

CONTAINS


SUBROUTINE carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

   USE cable_def_types_mod, ONLY : soil_parameter_type, veg_parameter_type,    &
                                   soil_snow_type, canopy_type, bgc_pool_type, &
                                   mp, mvtype

   USE cable_common_module, ONLY : cable_user

    REAL, INTENT(IN) ::                                                       &
      dels     ! integration time step (s)

    TYPE(soil_snow_type), INTENT(IN)     :: ssnow  ! soil/snow variables
    TYPE(veg_parameter_type), INTENT(IN) :: veg    ! vegetation parameters
    TYPE(canopy_type), INTENT(IN)        :: canopy ! canopy/veg variables
    TYPE(bgc_pool_type), INTENT(INOUT)   :: bgc    ! biogeochemistry variables

    TYPE(soil_parameter_type), INTENT(IN):: soil   ! soil parameters

    REAL, PARAMETER     :: beta = 0.9

    REAL, DIMENSION(mp) ::                                                     &
      cfsf,       & ! fast soil carbon turnover
      cfrts,      & ! roots turnover
      cfwd,       & ! wood turnover
      fcl,        & ! fraction of assimilated carbon that
                     ! goes to the construction of leaves  (eq. 5)
      fr,         & !
      clitt,      & !
      coef_cd,    & ! total stress coeff. for veg (eq. 6)
      coef_cold,  & ! coeff. for cold stress (eq. 7)
      coef_drght, & ! coeff. for drought stress (eq. 8)
      wbav          ! water stress index

    REAL, DIMENSION(:), ALLOCATABLE ::                                         &
      rw,      & !
      tfcl,    & !
      tvclst,  & !
      trnl,    & !
      trnr,    & !
      trnsf,   & !
      trnw


   ALLOCATE( rw(mvtype), tfcl(mvtype), tvclst(mvtype),                        &
              trnl(mvtype), trnr(mvtype), trnsf(mvtype), trnw(mvtype) )

   trnl = 3.17e-8
   trnr = 4.53e-9
   trnsf = 1.057e-10
   trnw = 6.342e-10

   SELECT CASE (mvtype)

      CASE (13)     ! CASA vegetation types

         rw   = (/ 16., 8.7, 12.5, 16., 18., 7.5,                              &
                 6.1, .84, 10.4, 15.1, 9., 5.8, 0.001 /)

         tfcl = (/ 0.248, 0.345, 0.31, 0.42, 0.38, 0.35,                       &
                 0.997, 0.95, 2.4, 0.73, 2.4, 0.55, 0.9500 /)

         tvclst = (/ 283., 278., 278., 235., 268., 278.0,                      &
                   278.0, 278.0, 278.0, 235., 278., 278., 268. /)

      CASE (15)     ! CSIRO types for UM

         rw   = (/ 16., 16., 18., 8.7, 10.4, 6.1, 6.1, 6.1,                    &
                  5.8, 5.8, 0.001, 9.0, 0.001, 0.001, 0.001 /)

         tfcl = (/ 0.42, 0.248, 0.38, 0.345, 2.4, 0.997, 0.997, 0.997,         &
                  0.55, 0.55, 0.9500, 2.4, 0.9500, 0.9500, 0.9500 /)

         tvclst = (/ 235., 283., 268., 278., 278.0, 278.0, 278.0, 278.0,       &
                    278., 278., 278.0, 278., 278., 278., 268. /)

      CASE (16)     ! IGBP vegetation types without water bodies

         rw   = (/ 16., 16., 18., 8.7, 12.5, 15.1, 10.4, 7.5,                  &
                & 6.1, 6.1, 0.001, 5.8, 0.001, 5.8, 0.001, 9.0 /)

         tfcl = (/ 0.42, 0.248, 0.38, 0.345, 0.31, 0.73, 2.4, 0.35,            &
                & 0.997, 0.997, 0.9500, 0.55, 0.9500, 0.55, 0.9500, 2.4 /)

         tvclst = (/ 235., 283., 268., 278., 278., 235., 278.0, 278.0,         &
                  & 278.0, 278.0, 278.0, 278., 278., 278., 268., 278. /)

      CASE (17)     ! IGBP vegetation types with water bodies

         ! rml: may not be the best values for our current 17 types,
         ! but will be superceeded by CASA-CNP anyway
         rw   = (/ 16., 16., 18., 8.7, 12.5, 15.1, 10.4, 7.5,                  &
                 6.1, 6.1, 0.001, 5.8, 0.001, 5.8, 0.001, 9.0, 0.001 /)

         tfcl = (/ 0.42, 0.248, 0.38, 0.345, 0.31, 0.73, 2.4, 0.35, 0.997,     &
                & 0.997, 0.9500, 0.55, 0.9500, 0.55, 0.9500, 2.4, 0.9500 /)

         tvclst = (/ 235., 283., 268., 278., 278., 235., 278.0, 278.0,         &
                   278.0, 278.0, 278.0, 278., 278., 278., 268., 278., 278. /)

      CASE DEFAULT

         PRINT *, 'Error! Dimension not compatible with CASA ',                 &
                  'or CSIRO or IGBP types!'
         PRINT *, 'Dimension =', mvtype
         PRINT *, 'At the rw section.'
         STOP

   END SELECT

   ! Limit size of exponent to avoif overflow when tv is very cold
   coef_cold = EXP( MIN( 1., -( canopy%tv - tvclst( veg%iveg ) ) ) )
   wbav = SUM( veg%froot * real(ssnow%wb), 2)
   wbav = max( 0.01, wbav )  ! EAK Jan2011

   ! drought stress
   coef_drght = EXP( 5.*( MIN( 1., MAX( 1., wbav**( 2 - soil%ibp2 ) - 1.) /    &
                ( soil%swilt**( 2 - soil%ibp2 ) - 1. ) ) - 1. ) )

   coef_cd = ( coef_cold + coef_drght ) * 2.0e-7

   ! CARBON POOLS : fraction of assimilated carbon
   ! that goes to the construction of leaves (eq. 5)
   fcl = EXP( - tfcl( veg%iveg ) * veg%vlai )

   ! LEAF: resp_lfrate is omitted below as fpn represents photosythesis - leaf
   ! transpiration calculated by the CBM
   clitt = ( coef_cd + trnl( veg%iveg ) ) * bgc%cplant(:,1)

   bgc%cplant(:,1) = bgc%cplant(:,1) - dels * ( canopy%fpn * fcl + clitt )

   !    WOOD:
   ! fraction of photosynthate going to roots, (1-fr) to wood, eq. 9
   !!vh!! inserted '0.0001' to avoide floating pt underflow
   fr = MIN( 1., EXP (- rw( veg%iveg ) * beta *0.0001* bgc%cplant(:,3)         &
        / MAX( bgc%cplant(:,2), 0.01 ) ) / beta )

   cfwd = trnw(veg%iveg) * bgc%cplant(:,2)

   bgc%cplant(:,2) = bgc%cplant(:,2) - dels * ( canopy%fpn * ( 1. - fcl )      &
                     * ( 1. - fr ) + canopy%frpw + cfwd )

   ! ROOTS
   cfrts = trnr( veg%iveg ) * bgc%cplant(:,3)

   bgc%cplant(:,3) = bgc%cplant(:,3) - dels * ( canopy%fpn * ( 1. - fcl ) * fr &
                     + cfrts + canopy%frpr )

   ! SOIL
   ! fast carbon
   cfsf = trnsf( veg%iveg ) * bgc%csoil(:,1)
   bgc%csoil(:,1) = bgc%csoil(:,1) + dels * ( 0.98 * clitt + 0.9 * cfrts       &
                    + cfwd - cfsf - 0.98 * canopy%frs )

   ! slow carbon
   bgc%csoil(:,2) = bgc%csoil(:,2) + dels * ( 0.02 * clitt  + 0.1 * cfrts      &
                    + cfsf - 0.02 * canopy%frs)

   !rml 17/1/11 change minimum pool size from 0.001 to 0.
   !(since want 0. for vegtype=ice)
   bgc%cplant(:,:)  = MAX(0.00, bgc%cplant(:,:))

   bgc%csoil(:,:) = MAX(0.00, bgc%csoil(:,:))

   DEALLOCATE( rw, tfcl, tvclst, trnl, trnr, trnsf, trnw )

END SUBROUTINE carbon_pl

! -----------------------------------------------------------------------------

SUBROUTINE soilcarb( soil, ssnow, veg, bgc, met, canopy)

   USE cable_def_types_mod, ONLY : soil_parameter_type, veg_parameter_type,    &
                                   soil_snow_type, canopy_type, bgc_pool_type, &
                                   met_type, mp, mvtype, ms, mstype

   USE cable_common_module

   TYPE (soil_snow_type), INTENT(IN)        :: ssnow
   TYPE (bgc_pool_type), INTENT(IN)         :: bgc
   TYPE (met_type), INTENT(IN)              :: met
   TYPE (canopy_type), INTENT(INOUT)        :: canopy

   TYPE (soil_parameter_type), INTENT(IN)   :: soil
   TYPE (veg_parameter_type), INTENT(IN)    :: veg

   REAL, DIMENSION(mp) ::                                                      &
      den,        & ! sib3
      rswc,       & !
      sss,        & !
      e0rswc,     & !
      ftsoil,     & !
      ftsrs,      & !
      tref,       & !
      tsoil,      & !
      avgtrs,     & ! root weighted mean soil temperature
      avgwrs        !root weighted mean soil moisture

   REAL, DIMENSION(mstype) ::                                                  &
      rswch,         & !
      soilcf           !

   REAL, PARAMETER :: t0 = -46.0

   INTEGER :: k

   IF( cable_user%DIAG_SOIL_RESP == 'off' .OR.                                 &
       cable_user%DIAG_SOIL_RESP == 'OFF' )  THEN

      ! key parameter for this scheme is veg%rs20

      avgwrs = SUM( veg%froot * real(ssnow%wb), 2 )
      avgtrs = MAX( 0.0, SUM( veg%froot * ssnow%tgg, 2 )- C%TFRZ )

      canopy%frs = veg%rs20 * MIN( 1.0, MAX( 0.0, MIN(                         &
                   -0.0178 + 0.2883 * avgwrs + 5.0176 * avgwrs * avgwrs        &
                   -4.5128 * avgwrs * avgwrs * avgwrs,                         &
                   0.3320+22.6726*exp( -5.8184*avgwrs ) ) ) )                  &
                   * MIN( 1.0, MAX( 0.0, MIN( 0.0104 * ( avgtrs**1.3053 ),     &
                   5.5956 - 0.1189 * avgtrs ) ) )

      canopy%frs = canopy%frs * SUM( SPREAD( bgc%ratecs, 1, mp ) * bgc%csoil,  &
                   2 ) / ( 365.0 * 24.0 * 3600.0 )  !convert 1/year to 1/second

      WHERE (ssnow%snowd > 1.)                                                 &
         canopy%frs = canopy%frs / max(0.001,min(100.,ssnow%snowd))

    ELSE

      ! key parameter for this scheme is veg%vegcf

      rswch = 0.16
      soilcf = 1.0

      den = MAX( 0.07,soil%sfc - soil%swilt )

      rswc = MAX(0.0001, veg%froot(:,1)*(REAL(ssnow%wb(:,2)) - soil%swilt))&
         & / den
      tsoil = veg%froot(:,1) * ssnow%tgg(:,2) - C%TFRZ

      tref = MAX( 0., ssnow%tgg(:,ms) - (C%TFRZ-.05 ) )

      DO k = 2,ms

         rswc = rswc + MAX( 0.0001, veg%froot(:,k)                             &
                * ( REAL( ssnow%wb(:,k) ) - soil%swilt ) ) / den

         tsoil = tsoil + veg%froot(:,k) * ssnow%tgg(:,k)

      ENDDO

      rswc = MIN( 1., rswc )
      tsoil = MAX( t0 + 2., tsoil)
      e0rswc = 52.4 + 285. * rswc
      ftsoil = MIN( 0.0015, 1./ ( tref - t0 ) - 1. / ( tsoil - t0 ) )
      sss = MAX( -15., MIN( 1., e0rswc * ftsoil ) )
      ftsrs=EXP( sss )
      canopy%frs = veg%vegcf * ( 144.0 / 44.0e6 )                              &
                   * soilcf(soil%isoilm) * MIN( 1., 1.4 * MAX( .3, .0278 *     &
                   tsoil + .5 ) ) * ftsrs * rswc / ( rswch ( soil%isoilm ) +   &
                   rswc )

   ENDIF

  END SUBROUTINE soilcarb

! -----------------------------------------------------------------------------

! plant respiration subroutine
SUBROUTINE plantcarb(veg, bgc, met, canopy)

   USE cable_def_types_mod, ONLY : veg_parameter_type, met_type,               &
                                   canopy_type, bgc_pool_type, mp , mvtype


   TYPE (veg_parameter_type), INTENT(IN)    :: veg
   TYPE (bgc_pool_type), INTENT(IN)         :: bgc
   TYPE (met_type), INTENT(IN)              :: met
   TYPE (canopy_type), INTENT(INOUT)        :: canopy

   REAL, DIMENSION(mp) ::                                                      &
      poolcoef1,     &! non-leaf carbon turnover rate * non-leaf pool size
      poolcoef1w,    & ! wood carbon turnover rate * wood pool size
      poolcoef1r,    & ! root carbon turnover rate * root pool size
      tmp1,tmp2,tmp3   ! kludge

   REAL, PARAMETER :: sec_per_year  = 365.0*24.0*3600.0

   call point2constants(C)

   poolcoef1=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -                      &
         bgc%ratecp(1)*bgc%cplant(:,1))

   poolcoef1w=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -                     &
        bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(3)*bgc%cplant(:,3))

   poolcoef1r=(sum(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -                     &
        bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(2)*bgc%cplant(:,2))

   tmp1(:) = 3.22 - 0.046 * (met%tk(:)-C%TFRZ)
   tmp2(:) = 0.1 * (met%tk(:)-C%TFRZ-20.0)
   tmp3(:) = tmp1(:) ** tmp2(:)

   canopy%frp  = veg%rp20 * tmp3 * poolcoef1  / sec_per_year
   canopy%frpw = veg%rp20 * tmp3 * poolcoef1w / sec_per_year
   canopy%frpr = veg%rp20 * tmp3 * poolcoef1r / sec_per_year

END SUBROUTINE plantcarb



END MODULE cable_carbon_module
