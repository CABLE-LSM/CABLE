!******************************************************************************
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) can be found
! at https://github.com/CABLE-LSM/CABLE/blob/main/
!******************************************************************************

MODULE cbl_albedo_mod

!-----------------------------------------------------------------------------
! Description:
!   Computes 4-band (visible/near-infrared, beam/diffuse) reflectances
!   and albedo
!
! IMPORTANT NOTE regarding the masks (Ticket 333)
! Prior to #333, 3 masks were used here - veg_mask, sunlit_mask, veg_mask
! - although passed in sunlit_mask was not used
! For JAC we will not be able to populate the sunlit masks and instead will
! evaluate the EffExtCoeff outside their bounds of applicability here by
! using inclusive masks in their place from the calling routines,
!  ie. JAC will use veg_mask in place of veg_mask
!
! To avoid confusion the mask names here are renamed:
! sunlit_mask now called mask1, veg_mask now called mask2
!
! This MODULE is USEd by:
!      cable_cbm.F90 (ESM1.5),
!      cable_rad_driver.F90 (ESM1.5),
!      cbl_model_driver_offline.F90 (CABLE),
!      rad_driver_cbl.F90 (JULES)
!
! This MODULE contains 1 public Subroutine:
!      albedo
! Other subroutines:
!      CanopyReflectance,
!      CanopyTransmitance,
!      EffectiveSurfaceReflectance,
!      EffectiveReflectance,
!      FbeamRadAlbedo,
!
! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC :: albedo
PRIVATE

CONTAINS

SUBROUTINE Albedo( AlbSnow, AlbSoil, mp, nrb, ICE_SoilType, lakes_cable,       &
                   jls_radiation, veg_mask, Ccoszen_tols, cgauss_w,            &
                   SurfaceType, SoilType, VegRefl, VegTaul,                    &
                   coszen, reducedLAIdue2snow, SnowDepth, SnowDensity,         &
                   SoilTemp, SnowAge, xk, c1, rhoch, RadFbeam, RadAlbedo,      &
                   ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,              &
                   EffExtCoeff_dif, CanopyRefl_beam,CanopyRefl_dif,            &
                   CanopyTransmit_beam, CanopyTransmit_dif,                    &
                   EffSurfRefl_beam, EffSurfRefl_dif)

! Description:
!   Computes the 4-band (visible/near-infrared, beam/diffuse) reflectances
!   and albedo

!subrs called
USE cbl_snow_albedo_module, ONLY: surface_albedosn

IMPLICIT NONE

!model dimensions
INTEGER, INTENT(IN) :: mp              ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb             ! # rad bands: VIS,NIR(,LW)

! Return: Effective Surface Relectance as seen by atmosphere
REAL, INTENT(OUT) :: AlbSnow(mp,nrb)   ! Alb adjusted for snow (ssnow%albsoilsn)
REAL, INTENT(OUT) :: EffSurfRefl_beam(mp,nrb)    ! Direct Beam [SW] (rad%reffbm)
REAL, INTENT(OUT) :: EffSurfRefl_dif(mp,nrb)     ! Diffuse [SW]  (rad%reffdf)
REAL, INTENT(OUT) :: RadAlbedo(mp,nrb)  ! Tot albedo given RadFbeam (rad%albedo)

!--- IN: CABLE specific surface_type indexes
INTEGER, INTENT(IN) :: ICE_SoilType    ! Index of ICE soil type
INTEGER, INTENT(IN) :: lakes_cable     ! Index of lakes surface type

!constants
REAL, INTENT(IN) :: Ccoszen_tols       ! zenith angle threshold for SUNLIT
REAL, INTENT(IN) :: Cgauss_w(nrb)
LOGICAL, INTENT(IN) :: jls_radiation   ! runtime switch = radiation pathway here

!masks
LOGICAL, INTENT(IN) :: veg_mask(mp)    ! TRUE where vegetated

!Vegetation parameters
REAL, INTENT(IN) :: VegTaul(mp,nrb)    ! PARAMETER leaf transmisivity (veg%taul)
REAL, INTENT(IN) :: VegRefl(mp,nrb)    ! PARAMETER leaf reflectivity (veg%refl)
INTEGER, INTENT(IN) :: SurfaceType(mp) ! Integer index of Surface type (veg%iveg)
INTEGER, INTENT(IN) :: SoilType(mp)    ! Integer index of Soil    type (soil%isoilm)

REAL :: reducedLAIdue2snow(mp)         ! Reduced LAI given snow coverage

! Albedos
REAL, INTENT(IN) :: AlbSoil(mp,nrb)    ! Param'ed Bare Soil Albedo(soil%albsoil)

!Forcing
REAL, INTENT(IN)  :: coszen(mp)        ! cosine zenith angle  (met%coszen)

!Prognostics
REAL, INTENT(IN) :: SnowDepth(mp)      ! Total Snow depth - water eqivalent
REAL, INTENT(IN) :: SnowDensity(mp)    ! Total Snow density (assumes 1 layer)
REAL, INTENT(IN) :: SoilTemp(mp)       ! Top layer Soil Temp. (ssnow%tgg)
REAL, INTENT(IN) :: SnowAge(mp)        ! Snow age (assumes 1 layer )

! Variables shared primarily between radiation and albedo and possibly elsewhere
REAL, INTENT(IN) :: RadFbeam(mp,nrb)   ! Beam Fraction of total SW (rad%fbeam)
!common radiation scalings - computed  in init_radiation()
REAL, INTENT(IN) :: xk(mp,nrb)
REAL, INTENT(IN) :: c1(mp,nrb)
REAL, INTENT(IN) :: rhoch(mp,nrb)

!Raw Extinction co-efficients computed in init_radiation()
REAL, INTENT(IN) :: ExtCoeff_beam(mp)  ! Direct Beam component (rad%extkb)
REAL, INTENT(IN) :: ExtCoeff_dif(mp)   ! Diffuse component(rad%extkd)
!Effecctive Extinction co-efficients computed in init_radiation()
REAL, INTENT(IN) :: EffExtCoeff_beam(mp,nrb) ! Direct Beam (rad%extkbm)
REAL, INTENT(IN) :: EffExtCoeff_dif(mp,nrb)  ! Diffuse component (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo()
REAL, INTENT(OUT) :: CanopyRefl_beam(mp,nrb)     ! Beam Canopy refl (rad%rhocbm)
REAL, INTENT(OUT) :: CanopyRefl_dif(mp,nrb)      ! Difuse Canopy refl(rad%rhocdf
REAL, INTENT(OUT) :: CanopyTransmit_beam(mp,nrb) ! Beam Transmit   (rad%cexpkdm)
REAL, INTENT(OUT) :: CanopyTransmit_dif(mp,nrb)  ! Difuse Transmit (rad%cexpkbm)

! local vars
INTEGER :: i

! END header

AlbSnow(:,:) = 0.0
CanopyRefl_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0
!CanopyTransmit_beam(:,:) = 1.0 ! Flushing this =1or0 changes answer - implying?
!CanopyTransmit_dif(:,:) = 1.0  ! MPI (at least inits this = 1.0 at dt=0)

!Modify parametrised soil albedo based on snow coverage
CALL surface_albedosn( AlbSnow, AlbSoil, mp, nrb, ICE_SoilType, lakes_cable,   &
                       SurfaceType, SoilType, SnowDepth, SnowDensity,          &
                       SoilTemp, SnowAge, Coszen )

! Update fractional leaf transmittance and reflection
!---1 = visible, 2 = nir radiaition

! Define canopy Reflectance for diffuse/direct radiation
! Formerly rad%rhocbm, rad%rhocdf
CALL CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif,                       &
                        mp, nrb, CGauss_w, veg_mask,                           &
                        AlbSnow, xk, rhoch,                                    &
                        ExtCoeff_beam, ExtCoeff_dif)

! Define canopy diffuse transmittance
! Formerly rad%cexpkbm, rad%cexpkdm
CALL CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,      &
                        veg_mask, reducedLAIdue2snow,                          &
                        EffExtCoeff_beam, EffExtCoeff_dif)

!---1 = visible, 2 = nir radiaition
! Finally compute Effective 4-band albedo for diffuse/direct radiation-
! In the UM this is the required variable to be passed back on the rad call
! Formerly rad%reffbm, rad%reffdf

! Even when there is no vegetation, albedo is at least snow modified soil albedo
EffSurfRefl_beam = AlbSnow
EffSurfRefl_dif = AlbSnow

CALL EffectiveSurfaceReflectance( EffSurfRefl_beam, EffSurfRefl_dif,           &
                                  mp, nrb, veg_mask, CanopyRefl_beam,          &
                                  CanopyRefl_dif, CanopyTransmit_beam,         &
                                  CanopyTransmit_dif, AlbSnow )

! Compute total albedo to SW given the Effective Surface Reflectance
! (considering Canopy/Soil/Snow contributions)
! we dont need to do this on rad call AND may not haveappropriate RadFbeam
RadAlbedo = AlbSnow
IF (.NOT. jls_radiation)                                                       &
  CALL FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam,                 &
                       EffSurfRefl_beam, EffSurfRefl_dif, AlbSnow )

RETURN
END SUBROUTINE albedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif,                 &
                              mp, nrb, CGauss_w, veg_mask,                     &
                              AlbSnow, xk, rhoch,                              &
                              ExtCoeff_beam, ExtCoeff_dif )
! Description:
!   Computes canopy Reflectance for diffuse/direct radiation

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                    ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb                   ! # rad bands [VIS,NIR(,LW)]
REAL, INTENT(OUT) :: CanopyRefl_dif(mp,nrb)  ! Beam Canopy refl (rad%rhocbm)
REAL, INTENT(OUT) :: CanopyRefl_beam(mp,nrb) ! Difuse Canopy refl(rad%rhocdf
REAL, INTENT(IN) :: Cgauss_w(nrb)
LOGICAL, INTENT(IN) :: veg_mask(mp)          ! TRUE where vegetated
REAL, INTENT(IN) :: AlbSnow(mp,nrb)          ! Ground Albedo with snow
REAL, INTENT(IN) :: xk(mp,nrb)
REAL, INTENT(IN) :: rhoch(mp,nrb)
REAL, INTENT(IN) :: ExtCoeff_beam(mp)        ! Extinction co-efficient
REAL, INTENT(IN) :: ExtCoeff_dif(mp)         ! Extinction co-efficient
!local
INTEGER :: i, b

! Canopy reflection (6.21) beam:
DO i = 1,mp
  DO b = 1, (nrb-1) !because nrb=3 due to legacy
    IF ( veg_mask(i) )                                                         &
      CanopyRefl_beam(i,b) = 2.0 * ExtCoeff_beam(i) /                          &
                            ( ExtCoeff_beam(i) + ExtCoeff_dif(i) )             &
                            * rhoch(i,b)
  END DO
END DO

! Canopy REFLection of diffuse radiation for black leaves:
DO i=1,(nrb-1) !because nrb=3 due to legacy

  CanopyRefl_dif(:,i) = rhoch(:,i) *  2.0 *                                    &
                       ( cgauss_w(1) * xk(:,1) / ( xk(:,1) + ExtCoeff_dif(:) ) &
                       + cgauss_w(2) * xk(:,2) / ( xk(:,2) + ExtCoeff_dif(:) ) &
                       + cgauss_w(3) * xk(:,3) / ( xk(:,3) + ExtCoeff_dif(:) ) )

END DO

RETURN
END SUBROUTINE CanopyReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyTransmitance( CanopyTransmit_beam, CanopyTransmit_dif, mp,    &
                               nrb,mask, reducedLAIdue2snow, EffExtCoeff_beam, &
                               EffExtCoeff_dif )
! Description:
!   Computes canopy diffuse transmittance

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                           ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb                          ! # rad bands [VIS,NIR(,LW)]
REAL, INTENT(OUT) :: CanopyTransmit_dif(mp,nrb)     ! Transmitance (rad%cexpkdm)
REAL, INTENT(OUT) :: CanopyTransmit_beam(mp,nrb)    ! Transmitance (rad%cexpkbm)
LOGICAL, INTENT(IN) :: mask(mp)                     ! TRUE where vegetated
REAL, INTENT(IN) :: reducedLAIdue2snow(mp)          ! Reduced LAI given snow
REAL, INTENT(IN) :: EffExtCoeff_beam(mp,nrb)        ! Extinction co-efficient
REAL, INTENT(IN) :: EffExtCoeff_dif(mp,nrb)         ! Extinction co-efficient
!local vars
INTEGER :: i,b
REAL :: dummy(mp,nrb)

! For beam, compute canopy trasmitance when sunlit (and vegetated)
DO i = 1,mp
  DO b = 1,(nrb-1)
    IF ( mask(i) ) THEN
      dummy(i,b) = MIN( EffExtCoeff_beam(i,b) * reducedLAIdue2snow(i), 20.0 )
      CanopyTransmit_beam(i,b) = EXP( -1.0* dummy(i,b) )
    END IF
  END DO
END DO

! For diffuse rad, always compute canopy trasmitance
DO i = 1,mp
  DO b = 1, (nrb-1) !because nrb=3 due to legacy
    dummy(i,b) = EffExtCoeff_dif(i,b) * reducedLAIdue2snow(i)
    CanopyTransmit_dif(i,b) = EXP( -1.0* dummy(i,b) )
  END DO
END DO

RETURN
END SUBROUTINE CanopyTransmitance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EffectiveSurfaceReflectance(EffSurfRefl_beam, EffSurfRefl_dif,      &
                                       mp, nrb, veg_mask, CanopyRefl_beam,     &
                                       CanopyRefl_dif, CanopyTransmit_beam,    &
                                       CanopyTransmit_dif, AlbSnow )
! Description:
!   Computes the effective 4-band albedo for diffuse/direct radiation.

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                       ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb                      ! # rad bands [VIS,NIR(,LW)
REAL, INTENT(OUT) :: EffSurfRefl_dif(mp,nrb)    ! Diffuse Surf Refl (rad%reffdf)
REAL, INTENT(OUT) :: EffSurfRefl_beam(mp,nrb)   ! Beam Surf Refl (rad%reffbm)
LOGICAL, INTENT(IN) :: veg_mask(mp)             ! TRUE where vegetated
REAL, INTENT(IN) :: CanopyRefl_beam(mp,nrb)     ! Beam Canopy refl (rad%rhocbm)
REAL, INTENT(IN) :: CanopyRefl_dif(mp,nrb)      ! Difuse Canopy refl(rad%rhocdf)
REAL, INTENT(IN) :: CanopyTransmit_dif(mp,nrb)  ! Transmitance (rad%cexpkdm)
REAL, INTENT(IN) :: CanopyTransmit_beam(mp,nrb) ! Transmitance (rad%cexpkbm)
REAL, INTENT(IN) :: AlbSnow(mp,nrb)           ! snow adjustd Alb ssnow%albsoilsn

!Surface reflectance to Difuse Radiation
CALL EffectiveReflectance( EffSurfRefl_dif, mp, nrb, CanopyRefl_dif, AlbSnow,  &
                           CanopyTransmit_dif, veg_mask )

!Surface reflectance to Direct Beam  Radiation
CALL EffectiveReflectance( EffSurfRefl_beam, mp, nrb, CanopyRefl_beam, AlbSnow,&
                           CanopyTransmit_beam, veg_mask )

RETURN
END SUBROUTINE EffectiveSurfaceReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EffectiveReflectance( EffRefl, mp, nrb, CanopyRefl, AlbSnow,        &
                                 CanopyTransmit, mask )
! Description:
!   Computes the Surface reflectance for a given radiation

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                   ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb                  ! # rad bands [VIS,NIR(,LW)
REAL, INTENT(OUT) :: EffRefl(mp,nrb)        ! Effective Surface Reflectance
REAL, INTENT(IN) :: AlbSnow(mp,nrb)         ! snow adjustd Alb ssnow%albsoilsn
REAL, INTENT(IN) :: CanopyRefl(mp,nrb)      ! Beam Canopy refl (rad%rhocbm)
REAL, INTENT(IN) :: CanopyTransmit(mp,nrb)  ! Difuse Canopy refl(rad%rhocdf nrb)
LOGICAL, INTENT(IN) :: mask(mp)            ! TRUE where vegetated
!local vars
INTEGER :: i,b

DO i = 1,mp
  DO b = 1, (nrb-1) !because nrb=3 due to legacy
    IF ( mask(i) ) THEN

       ! Calculate effective beam reflectance (fraction):
      EffRefl(i,b) = CanopyRefl(i,b) + ( AlbSnow(i,b) - CanopyRefl(i,b) )      &
                         * CanopyTransmit(i,b)**2

    END IF
  END DO
END DO

RETURN
END SUBROUTINE EffectiveReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam,             &
                           EffSurfRefl_beam, EffSurfRefl_dif, AlbSnow )
! Description:
!   Computes total albedo to SW given the Effective Surface Reflectance
!   (considering Canopy/Soil/Snow contributions)

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                    ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb                   ! # rad bands [VIS,NIR(,LW)
REAL, INTENT(OUT) :: RadAlbedo(mp,nrb)        ! RadFbeam albedo (rad%albedo)
LOGICAL, INTENT(IN) :: veg_mask(mp)          ! TRUE where vegetated
REAL, INTENT(IN) :: AlbSnow(mp,nrb)          ! snow adjustd Alb (ssnow%albsoilsn
REAL, INTENT(IN) :: RadFbeam(mp,nrb)        ! Beam Fraction of SW (rad%fbeam)
REAL, INTENT(IN) :: EffSurfRefl_dif(mp,nrb)  ! Direct Beam [SW] (rad%reffbm)
REAL, INTENT(IN) :: EffSurfRefl_beam(mp,nrb) ! Diffuse [SW]  (rad%reffdf)
!local vars
INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
INTEGER :: i

! Initialise total albedo:
RadAlbedo = AlbSnow
DO i = 1,mp
  DO b = 1,(nrb-1) !because nrb=3 due to legacy
    ! Define albedo:
    IF ( veg_mask(i) )                                                         &
       RadAlbedo(i,b) = ( 1.0 - radfbeam(i,b) )*EffSurfRefl_dif(i,b)           &
                        + radfbeam(i,b) * EffSurfRefl_beam(i,b)
  END DO
END DO

RETURN
END SUBROUTINE FbeamRadAlbedo

END MODULE cbl_albedo_mod
