!******************************COPYRIGHT********************************************
! (c) CSIRO 2022.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT********************************************
MODULE cable_rad_driv_mod

!-----------------------------------------------------------------------------
! Description:
!    Initialises radiation specific variables and computes the albedo for
!    CABLE.
!
! This MODULE is USEd in:
!     cable_land_albedo_mod_cbl.F90 (JULES)
!
! This MODULE contains 1 public Subroutine:
!     cable_rad_driver
!
! Code owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: cable_rad_driver

CONTAINS

SUBROUTINE cable_rad_driver( EffSurfRefl_beam, EffSurfRefl_dif,                &
                             mp, nrb, ICE_SoilType, lakes_type, Clai_thresh,   &
                             Ccoszen_tols, CGauss_w, Cpi, Cpi180, z0surf_min,  &
                             veg_mask, jls_standalone, jls_radiation,          &
                             SurfaceType,  SoilType,                           &
                             LAI_pft_cbl, HGT_pft_cbl, SnowDepth,              &
                             SnowDensity, SoilTemp, SnowAge, AlbSoil,          &
                             coszen, VegXfang, VegTaul, VegRefl,               &
                             HeightAboveSnow, reducedLAIdue2snow,              &
                             ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,    &
                             EffExtCoeff_dif, CanopyTransmit_beam,             &
                             CanopyTransmit_dif, CanopyRefl_beam,              &
                             CanopyRefl_dif, c1, rhoch, xk, AlbSnow, RadFbeam, &
                             RadAlbedo, metDoY, SW_down )

! Description:
! Nothing further to add to module description.

!subrs:
USE cbl_albedo_mod,             ONLY: albedo
USE cbl_init_radiation_module,  ONLY: init_radiation


IMPLICIT NONE

!model dimensions
INTEGER, INTENT(IN) :: mp       ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb      ! # radiation bands[ 1=VIS,2=NIR,3=LW(legacy)

! Albedos req'd by JULES - Effective Surface Relectance as seen by atmosphere
REAL, INTENT(OUT) :: EffSurfRefl_dif(mp,nrb)  ! Refl to Diffuse component of rad
                                              ! formerly rad%reffdf
REAL, INTENT(OUT) :: EffSurfRefl_beam(mp,nrb) ! Refl to Beam component of rad
                                              ! formerly rad%reffbm

!--- IN: CABLE specific surface_type indexes
INTEGER, INTENT(IN) :: ICE_SoilType
INTEGER, INTENT(IN) :: lakes_type

!constants
REAL, INTENT(IN) :: Ccoszen_tols      ! threshold cosine of sun's zenith angle,
                                      ! below which considered SUNLIT
REAL, INTENT(IN) :: Cgauss_w(nrb)     ! Gaussian integration weights
REAL, INTENT(IN) :: Clai_thresh       ! The minimum LAI below which a "cell" is
                                      ! considred  NOT vegetated
REAL, INTENT(IN) :: Cpi               ! PI
REAL, INTENT(IN) :: Cpi180            ! PI in radians
REAL, INTENT(IN) :: z0surf_min        ! the minimum roughness of bare soil

LOGICAL, INTENT(IN) :: jls_standalone ! local runtime switch for JULES(/UM) run
LOGICAL, INTENT(IN) :: jls_radiation  ! local runtime switch for radiation path

!masks
LOGICAL, INTENT(IN) :: veg_mask(:)         !  vegetated (uses min LAI)

!recieved as spatial maps from the UM. remapped to "mp"
REAL, INTENT(IN) :: LAI_pft_cbl(mp)        ! LAI -  "limited" and remapped
REAL, INTENT(IN) :: HGT_pft_cbl(mp)        ! canopy height -  "limited", remapped
REAL, INTENT(IN) :: coszen(mp)             ! cosine zenith angle  (met%coszen)
REAL, INTENT(IN) :: AlbSoil(mp, nrb)      ! soil%AlbSoil

!computed for CABLE model
REAL, INTENT(IN):: HeightAboveSnow(mp)     ! Height of Canopy above snow
                                           ! (rough%hruff) computed from
                                           ! z0surf_min, HGT_pft_cbl,
                                           ! SnowDepth, SnowDensity
REAL, INTENT(IN) :: reducedLAIdue2snow(mp) ! Reduced LAI given snow coverage

!Prognostics !recieved as spatial maps from the UM. remapped to "mp"
REAL, INTENT(IN) :: SnowDepth(mp)          ! Total Snow depth - water eqivalent -
                                           ! packed from snow_surft (ssnow%snowd)
                                           ! this timestep (ssnow%Osnowd)
REAL, INTENT(IN) :: SnowDensity(mp)        ! Total Snow density (assumes 1 layer
                                           ! describes snow cover) (ssnow%ssdnn)
REAL, INTENT(IN) :: SoilTemp(mp)           ! Soil Temperature of top layer (soil%tgg)
REAL, INTENT(IN) :: SnowAge(mp)            ! Snow age (assumes 1 layer describes snow
                                           ! cover) (ssnow%snage)

!Vegetation parameters !recieved as per PFT params from the UM. remapped to "mp"
INTEGER, INTENT(IN) :: SurfaceType(mp)
INTEGER, INTENT(IN) :: SoilType(mp)
REAL,    INTENT(IN) :: VegXfang(mp)        ! leaf angle PARAMETER (veg%xfang)
REAL,    INTENT(IN) :: VegTaul(mp,nrb)     ! PARAM leaf transmisivity (veg%taul)
REAL,    INTENT(IN) :: VegRefl(mp,nrb)     ! PARAM leaf reflectivity (veg%refl)

!local to Rad/Albedo pathway:
REAL, INTENT(IN OUT) :: ExtCoeff_beam(mp)            ! nee. rad%extkb,
REAL, INTENT(IN OUT) :: ExtCoeff_dif(mp)             ! nee. rad%extkd
REAL, INTENT(IN OUT) :: EffExtCoeff_beam(mp, nrb)    ! nee. rad%extkbm
REAL, INTENT(IN OUT) :: EffExtCoeff_dif(mp, nrb)     ! nee. rad%extkdm,
REAL, INTENT(IN OUT) :: CanopyTransmit_dif(mp, nrb)  ! nee. rad%cexpkdm
REAL, INTENT(IN OUT) :: CanopyTransmit_beam(mp, nrb) ! nee. rad%cexpkbm
REAL, INTENT(IN OUT) :: CanopyRefl_dif(mp, nrb)      ! nee. rad%rhocdf
REAL, INTENT(IN OUT) :: CanopyRefl_beam(mp, nrb)     ! nee. rad%rhocbm
REAL, INTENT(IN OUT) :: RadFbeam(mp, nrb)            ! nee. rad%fbeam
REAL, INTENT(IN OUT) :: RadAlbedo(mp, nrb)           ! nee. rad%albedo
REAL, INTENT(IN OUT) :: AlbSnow(mp, nrb)             ! nee. ssnow%AlbSoilsn
REAL, INTENT(IN OUT) :: c1(mp, nrb)                  ! common rad scalings
REAL, INTENT(IN OUT) :: rhoch(mp, nrb)               ! common rad scalings
REAL, INTENT(IN OUT) :: xk(mp, nrb)                  ! common rad scalings
! used in Calc of Beam calculation NOT on rad/albedo path.
! However Needed to fulfill arg list with dummy
INTEGER, INTENT(IN OUT) :: metDoY(mp)                 ! can pass DoY from current_time
REAL, INTENT(IN OUT)    :: SW_down(mp, nrb)           ! NA at surf_couple_rad layer

CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_driver"
LOGICAL :: cbl_standalone = .FALSE.

!Defines Extinction Coefficients to use in calculation of Canopy
!Reflectance/Transmitance.
CALL init_radiation( ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,            &
                     EffExtCoeff_dif, RadFbeam, c1, rhoch, xk,                 &
                     mp,nrb, Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180, &
                     cbl_standalone, jls_standalone, jls_radiation, subr_name, &
                     veg_mask, VegXfang, VegTaul, VegRefl, coszen, metDoY,     &
                     SW_down, reducedLAIdue2snow )

!Finally call albedo to get what we really need to fulfill contract with JULES
!Defines 4-stream albedos [VIS/NIR bands. direct beam/diffuse components] from
!considering albedo of Ground (snow?) and Canopy Reflectance/Transmitance.
CALL Albedo( AlbSnow, AlbSoil,                                                 &
             mp, nrb, ICE_SoilType, lakes_type, jls_radiation, veg_mask,       &
             Ccoszen_tols, cgauss_w,                                           &
             SurfaceType, SoilType ,VegRefl, VegTaul,                          &
             coszen, reducedLAIdue2snow,                                       &
             SnowDepth, SnowDensity, SoilTemp, SnowAge,                        &
             xk, c1, rhoch,                                                    &
             RadFbeam, RadAlbedo,                                              &
             ExtCoeff_beam, ExtCoeff_dif,                                      &
             EffExtCoeff_beam, EffExtCoeff_dif,                                &
             CanopyRefl_beam,CanopyRefl_dif,                                   &
             CanopyTransmit_beam, CanopyTransmit_dif,                          &
             EffSurfRefl_beam, EffSurfRefl_dif)

RETURN

END SUBROUTINE cable_rad_driver

END MODULE cable_rad_driv_mod

