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
! Purpose: Calculates surface albedo, including from snow covered surface
!
! Called from: cbm
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b (but was previously in cable_radiation)
!
!
! ==============================================================================

MODULE cable_albedo_module

  USE cable_data_module, ONLY : ialbedo_type, point2constants

  IMPLICIT NONE

  PUBLIC surface_albedo
  PRIVATE

  TYPE(ialbedo_type) :: C


CONTAINS


  SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy)

    USE cable_common_module
    USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &
         canopy_type, met_type, radiation_type,      &
         soil_snow_type, mp, r_2, nrb

    TYPE (canopy_type),INTENT(IN)       :: canopy
    TYPE (met_type),INTENT(INOUT)       :: met
    TYPE (radiation_type),INTENT(INOUT) :: rad
    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow

    TYPE (veg_parameter_type),INTENT(INOUT)  :: veg
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL(r_2), DIMENSION(mp)  ::                                                &
         dummy2, & !
         dummy

    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: c1, rhoch

    LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation

    INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave

    ! END header

    CALL point2constants(C)

    IF (.NOT. ALLOCATED(c1)) &
         ALLOCATE( c1(mp,nrb), rhoch(mp,nrb) )


    CALL surface_albedosn(ssnow, veg, met, soil)

    rad%cexpkbm = 0.0
    rad%extkbm  = 0.0
    rad%rhocbm  = 0.0

    ! Initialise effective conopy beam reflectance:
    rad%reffbm = ssnow%albsoilsn
    rad%reffdf = ssnow%albsoilsn
    rad%albedo = ssnow%albsoilsn

    ! Define vegetation mask:
    mask = canopy%vlaiw > C%LAI_THRESH .AND.                                    &
         ( met%fsd(:,1) + met%fsd(:,2) ) > C%RAD_THRESH

    CALL calc_rhoch( veg, c1, rhoch )

    ! Update extinction coefficients and fractional transmittance for
    ! leaf transmittance and reflection (ie. NOT black leaves):
    !---1 = visible, 2 = nir radiaition
    DO b = 1, 2

       rad%extkdm(:,b) = rad%extkd * c1(:,b)

       !--Define canopy diffuse transmittance (fraction):
       rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)

       !---Calculate effective diffuse reflectance (fraction):
       WHERE( canopy%vlaiw > 1e-2 )                                             &
            rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssnow%albsoilsn(:,b)             &
            - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2

       !---where vegetated and sunlit
       WHERE (mask)

          rad%extkbm(:,b) = rad%extkb * c1(:,b)

          ! Canopy reflection (6.21) beam:
          rad%rhocbm(:,b) = 2. * rad%extkb / ( rad%extkb + rad%extkd )          &
               * rhoch(:,b)

          ! Canopy beam transmittance (fraction):
          dummy2 = MIN(rad%extkbm(:,b)*canopy%vlaiw, 20.)
          dummy  = EXP(-dummy2)
          rad%cexpkbm(:,b) = REAL(dummy)

          ! Calculate effective beam reflectance (fraction):
          rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssnow%albsoilsn(:,b)             &
               - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2

       END WHERE

       ! Define albedo:
       WHERE( canopy%vlaiw> C%LAI_THRESH )                                      &
            rad%albedo(:,b) = ( 1. - rad%fbeam(:,b) )*rad%reffdf(:,b) +           &
            rad%fbeam(:,b) * rad%reffbm(:,b)

    END DO


  END SUBROUTINE surface_albedo

  ! ------------------------------------------------------------------------------

  SUBROUTINE surface_albedosn(ssnow, veg, met, soil)

    USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &
         met_type, soil_snow_type, mp
    USE cable_common_module

    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow
    TYPE (met_type),INTENT(INOUT)       :: met

    TYPE (veg_parameter_type),INTENT(INout)  :: veg
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL, DIMENSION(mp) ::                                                      &
         alv,     &  ! Snow albedo for visible
         alir,    &  ! Snow albedo for near infra-red
         ar1,     &  ! crystal growth  (-ve)
         ar2,     &  ! freezing of melt water
         ar3,     &  !
         dnsnow,  &  ! new snow albedo
         dtau,    &  !
         fage,    &  ! age factor
         fzenm,   &  !
         sfact,   &  !
         snr,     &  !
         snrat,   &  !
         talb,    &  ! snow albedo
         tmp         ! temporary value

    REAL, PARAMETER ::                                                          &
         alvo  = 0.95,  &  ! albedo for vis. on a new snow
         aliro = 0.70      ! albedo for near-infr. on a new snow

    soil%albsoilf = soil%albsoil(:,1)

    ! lakes: hard-wired number to be removed in future
    WHERE( veg%iveg == 16 )                                                     &
         soil%albsoilf = -0.022*( MIN( 275., MAX( 260., met%tk ) ) - 260. ) + 0.45

    WHERE(ssnow%snowd > 1. .AND. veg%iveg == 16 ) soil%albsoilf = 0.85

    sfact = 0.68

    WHERE (soil%albsoilf <= 0.14)
       sfact = 0.5
    ELSEWHERE (soil%albsoilf > 0.14 .AND. soil%albsoilf <= 0.20)
       sfact = 0.62
    END WHERE

    ssnow%albsoilsn(:,2) = 2. * soil%albsoilf / (1. + sfact)
    ssnow%albsoilsn(:,1) = sfact * ssnow%albsoilsn(:,2)

    ! calc soil albedo based on colour - Ticket #27
    IF (calcsoilalbedo) THEN
       CALL soilcol_albedo(ssnow, soil)
    END IF

    snrat=0.
    alir =0.
    alv  =0.

    WHERE ( ssnow%snowd > 1. .AND. .NOT. cable_runtime%um_radiation )

       ! new snow (cm H2O)
       dnsnow = MIN ( 1., .1 * MAX( 0., ssnow%snowd - ssnow%osnowd ) )

       ! Snow age depends on snow crystal growth, freezing of melt water,
       ! accumulation of dirt and amount of new snow.
       tmp = ssnow%isflag * ssnow%tggsn(:,1) + ( 1 - ssnow%isflag )            &
            * ssnow%tgg(:,1)
       tmp = MIN( tmp, C%TFRZ )
       ar1 = 5000. * (1. / (C%TFRZ-0.01) - 1. / tmp) ! crystal growth  (-ve)
       ar2 = 10. * ar1 ! freezing of melt water
       snr = ssnow%snowd / MAX (ssnow%ssdnn, 200.)

       WHERE (soil%isoilm == 9)
          ! permanent ice: hard-wired number to be removed in future version
          ar3 = .0000001
          !  NB. dsnow =1,assumes pristine snow; ignores soot etc. ALTERNATIVELY,
          !dnsnow = max (dnsnow, .5) !increase refreshing of snow in Antarctic
          dnsnow = 1.0
          snrat = 1.

       ELSEWHERE

          ! accumulation of dirt
          ar3 = .1
          ! snow covered fraction of the grid
          snrat = MIN (1., snr / (snr + .1) )

       END WHERE

       dtau = 1.e-6 * (EXP( ar1 ) + EXP( ar2 ) + ar3 ) * kwidth_gl

       WHERE (ssnow%snowd <= 1.0)
          ssnow%snage = 0.
       ELSEWHERE
          ssnow%snage = MAX (0.,(ssnow%snage+dtau)*(1.-dnsnow))
       END WHERE

       fage = 1. - 1. / (1. + ssnow%snage ) !age factor

       tmp = MAX( .17365, met%coszen )
       fzenm = MAX( 0.0, MERGE( 0.0,                                           &
            ( 1. + 1./2. ) / ( 1. + 2.*2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)

       ! use dry snow albedo for pernament land ice: hard-wired no to be removed
       WHERE (soil%isoilm == 9)

          tmp = 0.95 * (1.0 - 0.2 * fage)
          alv = .4 * fzenm * (1. - tmp) + tmp
          tmp = 0.75 * (1. - .5 * fage)

       END WHERE

       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo

    ENDWHERE        ! snowd > 0

    ! when it is called from cable_rad_driver (UM)
    ! no need to recalculate snage
    WHERE (ssnow%snowd > 1 .AND. cable_runtime%um_radiation )

       snr = ssnow%snowd / MAX (ssnow%ssdnn, 200.)

       WHERE (soil%isoilm == 9)
          ! permanent ice: hard-wired number to be removed
          snrat = 1.
       ELSEWHERE
          snrat = MIN (1., snr / (snr + .1) )
       END WHERE

       fage = 1. - 1. / (1. + ssnow%snage ) !age factor
       tmp = MAX (.17365, met%coszen )
       fzenm = MAX( 0., MERGE( 0.0,                                             &
            ( 1. + 1./2. ) / ( 1. + 2. * 2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)

       ! use dry snow albedo
       WHERE (soil%isoilm == 9)
          ! permanent ice: hard-wired number to be removed

          tmp = 0.95 * (1.0 - 0.2 * fage)
          alv = .4 * fzenm * (1. - tmp) + tmp
          tmp = 0.75 * (1. - .5 * fage)

       END WHERE

       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo

    ENDWHERE        ! snowd > 0

    IF(cable_user%SOIL_STRUC=='sli') THEN
       WHERE (ssnow%snowd.GT.1.0)
          snrat = 1.0   ! using default parameterisation, albedo is too low,
          ! inhibiting snowpack initiation
       ENDWHERE
    ENDIF
    ssnow%albsoilsn(:,2) = MIN( aliro,                                          &
         ( 1. - snrat ) * ssnow%albsoilsn(:,2) + snrat * alir)

    ssnow%albsoilsn(:,1) = MIN( alvo,                                           &
         ( 1. - snrat ) * ssnow%albsoilsn(:,1) + snrat * alv )

    WHERE (soil%isoilm == 9) ! use dry snow albedo: 1=vis, 2=nir
       ssnow%albsoilsn(:,1) = alvo - 0.05 ! al*o = albedo appropriate for new snow
       ssnow%albsoilsn(:,2) = aliro - 0.05 ! => here al*o LESS arbitrary aging 0.05
    END WHERE
    
!$    WHERE (soil%isoilm == 9)       ! use dry snow albedo   ! This bit appears instead in MYY code --rk4417
!$       ssnow%albsoilsn(:,2) = 0.82
!$       ssnow%albsoilsn(:,1) = 0.82
!$    END WHERE
    
    RETURN

  END SUBROUTINE surface_albedosn

  ! ------------------------------------------------------------------------------

  !jhan:subr was reintroduced here to temporarily resolve issue when
  !creating libcable.a  (repeated in cable_radiation.F90)
  SUBROUTINE calc_rhoch(veg,c1,rhoch)

    USE cable_def_types_mod, ONLY : veg_parameter_type
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg
    REAL, INTENT(INOUT), DIMENSION(:,:) :: c1, rhoch

    c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
    c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
    c1(:,3) = 1.

    ! Canopy reflection black horiz leaves
    ! (eq. 6.19 in Goudriaan and van Laar, 1994):
    rhoch = (1.0 - c1) / (1.0 + c1)

  END SUBROUTINE calc_rhoch

  ! -----------------------------------------------------------------------------
  ! subr to calc soil albedo based on colour - Ticket #27
  SUBROUTINE soilcol_albedo(ssnow, soil)

    USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
         r_2, mp, nrb
    ! Arguments
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow      ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil       ! soil parameters

    ! Local Variables
    INTEGER   :: ib
    REAL(r_2), DIMENSION(mp)      :: inc
    REAL(r_2), DIMENSION(mp, nrb) :: albsod,          & ! soil albedo (direct)
         albsoi             ! soil albedo (indirect)

    ! Look-up tables for soil albedo
    ! saturated soil albedos for 20 color classes and 2 wavebands (1=vis, 2=nir)
    REAL(r_2), DIMENSION(20,nrb) ::                                  &
         albsat,                                                                  &
         albdry

    REAL(r_2), PARAMETER, DIMENSION(20*nrb) ::                                  &
         albsat1D = (/ 0.25,0.23,0.21,0.20,0.19,0.18,0.17,0.16,0.15,0.14,0.13,    &
         0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05,0.04 ,             &
         0.50,0.46,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,    &
         0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,              &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                    &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /),                    &
                                ! dry soil albedos for 20 color classes and 2 wavebands (1=vis, 2=nir)
         albdry1D = (/  0.36,0.34,0.32,0.31,0.30,0.29,0.28,0.27,0.26,0.25,0.24,   &
         0.23,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,             &
         0.61,0.57,0.53,0.51,0.49,0.48,0.45,0.43,0.41,0.39,0.37,   &
         0.35,0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.16,             &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                   &
         0., 0., 0., 0., 0., 0., 0., 0., 0., 0. /)

    albsat = RESHAPE( albsat1D, (/20, nrb/) )
    albdry = RESHAPE( albdry1D, (/20, nrb/) )

    DO ib = 1,2 ! Number of wavebands (vis, nir)
       inc = MAX(0.11-0.40*ssnow%wb(:,1), 0._r_2)
       albsod(:,ib) = MIN(albsat(INT(soil%soilcol),ib)+inc, albdry(INT(soil%soilcol),ib))
       albsoi(:,ib) = albsod(:,ib)
    END DO
    ssnow%albsoilsn = REAL(0.5*(albsod + albsoi))

  END SUBROUTINE soilcol_albedo

END MODULE cable_albedo_module
