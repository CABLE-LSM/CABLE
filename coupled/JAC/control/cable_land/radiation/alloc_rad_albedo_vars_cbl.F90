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
MODULE alloc_rad_albedo_vars_mod

!-----------------------------------------------------------------------------
! Description:
!    Allocate and deallocate variables in the rad structure
!
! This MODULE is USEd in:
!     cable_land_albedo_mod_cbl.F90 (JULES)
!
! This MODULE contains 2 public Subroutine:
!     alloc_local_vars,
!     flush_local_vars
!
! Code owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE
PUBLIC :: alloc_local_vars
PUBLIC :: flush_local_vars
PRIVATE

CONTAINS

! Allocate vars in mp format
SUBROUTINE alloc_local_vars( EffSurfRefl_beam, EffSurfRefl_dif, mp, nrb,       &
                             reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,     &
                             HeightAboveSnow, coszen, ExtCoeff_beam,           &
                             ExtCoeff_dif, EffExtCoeff_beam, EffExtCoeff_dif,  &
                             CanopyTransmit_beam, CanopyTransmit_dif,          &
                             CanopyRefl_beam, CanopyRefl_dif, RadFbeam,        &
                             RadAlbedo, AlbSnow, c1, rhoch, xk, metDoY,        &
                             SnowDepth, SnowDensity, SoilTemp, SnowAge,        &
                             AlbSoil, SW_down, veg_mask )
! Description:
! Allocate variables in the rad structure

IMPLICIT NONE

INTEGER :: mp, nrb
REAL, INTENT(OUT), ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: EffSurfRefl_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: SnowDepth(:)
REAL, INTENT(OUT), ALLOCATABLE :: SnowDensity(:)
REAL, INTENT(OUT), ALLOCATABLE :: SoilTemp(:)
REAL, INTENT(OUT), ALLOCATABLE :: SnowAge( :)
REAL, INTENT(OUT), ALLOCATABLE :: AlbSoil(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: LAI_pft_cbl(:)      ! Prescribed LAI
REAL, INTENT(OUT), ALLOCATABLE :: HGT_pft_cbl(:)      ! Prescribed canopy height
REAL, INTENT(OUT), ALLOCATABLE :: reducedLAIdue2snow(:) ! Eff. LAI given snow
REAL, INTENT(OUT), ALLOCATABLE :: HeightAboveSnow(:)  ! Canopy hgt above snow
REAL, INTENT(OUT), ALLOCATABLE :: coszen(:)
REAL, INTENT(OUT), ALLOCATABLE :: ExtCoeff_beam(:)
REAL, INTENT(OUT), ALLOCATABLE :: ExtCoeff_dif(:)
REAL, INTENT(OUT), ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: AlbSnow(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: c1(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: rhoch(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: xk(:,:)
! these are dummies in JULES rad call but req'd to load arg lists
REAL, INTENT(OUT), ALLOCATABLE :: SW_down(:,:)        ! dummy
REAL, INTENT(OUT), ALLOCATABLE :: RadFbeam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: RadAlbedo(:,:)
INTEGER, INTENT(OUT), ALLOCATABLE :: metDoY(:)        ! can pass DoY from current_time
! vegetated mask required on albedo pathway
LOGICAL, INTENT(OUT), ALLOCATABLE :: veg_mask(:)


IF ( .NOT. ALLOCATED(reducedLAIdue2snow))  ALLOCATE(reducedLAIdue2snow(mp) )
IF ( .NOT. ALLOCATED(HeightAboveSnow) )    ALLOCATE(HeightAboveSnow(mp) )
IF ( .NOT. ALLOCATED(LAI_pft_cbl) )        ALLOCATE(LAI_pft_cbl(mp) )
IF ( .NOT. ALLOCATED(HGT_pft_cbl) )        ALLOCATE(HGT_pft_cbl(mp) )
IF ( .NOT. ALLOCATED(AlbSoil) )            ALLOCATE(AlbSoil(mp,nrb) )
IF ( .NOT. ALLOCATED(coszen) )             ALLOCATE(coszen(mp) )
IF ( .NOT. ALLOCATED(EffSurfRefl_dif) )    ALLOCATE(EffSurfRefl_dif(mp, nrb) )
IF ( .NOT. ALLOCATED(EffSurfRefl_beam) )   ALLOCATE(EffSurfRefl_beam(mp, nrb) )
IF ( .NOT. ALLOCATED(SnowDepth) )          ALLOCATE(SnowDepth(mp) )
IF ( .NOT. ALLOCATED(SnowDensity) )        ALLOCATE(SnowDensity(mp) )
IF ( .NOT. ALLOCATED(SoilTemp) )           ALLOCATE(SoilTemp(mp) )
IF ( .NOT. ALLOCATED(SnowAge) )            ALLOCATE(SnowAge(mp) )
IF ( .NOT. ALLOCATED(ExtCoeff_beam) )      ALLOCATE(ExtCoeff_beam(mp) )
IF ( .NOT. ALLOCATED(ExtCoeff_dif) )       ALLOCATE(ExtCoeff_dif(mp) )
IF ( .NOT. ALLOCATED(EffExtCoeff_beam) )   ALLOCATE(EffExtCoeff_beam(mp, nrb) )
IF ( .NOT. ALLOCATED(EffExtCoeff_dif) )    ALLOCATE(EffExtCoeff_dif(mp, nrb) )
IF ( .NOT. ALLOCATED(CanopyTransmit_dif))  ALLOCATE(CanopyTransmit_dif(mp, nrb))
IF ( .NOT. ALLOCATED(CanopyTransmit_beam)) ALLOCATE(CanopyTransmit_beam(mp,nrb))
IF ( .NOT. ALLOCATED(CanopyRefl_dif) )     ALLOCATE(CanopyRefl_dif(mp, nrb) )
IF ( .NOT. ALLOCATED(CanopyRefl_beam) )    ALLOCATE(CanopyRefl_beam(mp, nrb) )
IF ( .NOT. ALLOCATED(AlbSnow) )            ALLOCATE(AlbSnow(mp, nrb) )
IF ( .NOT. ALLOCATED(c1) )                 ALLOCATE(c1(mp, nrb) )
IF ( .NOT. ALLOCATED(rhoch) )              ALLOCATE(rhoch(mp, nrb) )
IF ( .NOT. ALLOCATED(xk) )                 ALLOCATE(xk(mp, nrb) )
IF ( .NOT. ALLOCATED(metDoY) )             ALLOCATE(metDoY(mp) )
IF ( .NOT. ALLOCATED(SW_down) )            ALLOCATE(SW_down(mp,nrb) )
IF ( .NOT. ALLOCATED(RadFbeam) )           ALLOCATE(RadFbeam(mp, nrb) )
IF ( .NOT. ALLOCATED(RadAlbedo) )          ALLOCATE(RadAlbedo(mp, nrb) )
IF (.NOT. ALLOCATED(veg_mask) ) THEN
  ALLOCATE( veg_mask(mp) )
END IF

EffSurfRefl_dif(:,:) = 0.0
EffSurfRefl_beam(:,:) = 0.0
SnowDepth = 0.0
SnowDensity = 0.0
SoilTemp = 0.0
SnowAge = 0.0
LAI_pft_cbl(:) = 0.0
HGT_pft_cbl(:) = 0.0
coszen(:) = 0.0
reducedLAIdue2snow(:) = 0.0
HeightAboveSnow(:) = 0.0
ExtCoeff_beam(:) = 0.0
ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0
EffExtCoeff_dif(:,:) = 0.0
CanopyTransmit_dif(:,:) = 0.0
CanopyTransmit_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0
CanopyRefl_beam(:,:) = 0.0
AlbSnow(:,:) = 0.0
rhoch(:,:) = 0.0
xk(:,:) = 0.0
c1(:,:) = 0.0
RadFbeam(:,:) = 0.0
RadAlbedo(:,:) = 0.0
SW_down(:,:) = 0.0
metDoY(:) = 0  !can pass DoY from current_time%
veg_mask(:) = .FALSE.

RETURN

END SUBROUTINE alloc_local_vars

!flush memory
SUBROUTINE flush_local_vars( EffSurfRefl_beam, EffSurfRefl_dif, SnowDepth,     &
                             SnowDensity, SoilTemp, SnowAge, AlbSoil,          &
                             reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,     &
                             HeightAboveSnow, coszen, ExtCoeff_beam,           &
                             ExtCoeff_dif, EffExtCoeff_beam, EffExtCoeff_dif,  &
                             CanopyTransmit_beam, CanopyTransmit_dif,          &
                             CanopyRefl_beam, CanopyRefl_dif, RadFbeam,        &
                             RadAlbedo, AlbSnow, c1, rhoch, xk, metDoY, SW_down)

! Description:
! Deallocate variables in the rad structure

IMPLICIT NONE

REAL, INTENT(IN OUT), ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: EffSurfRefl_beam(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: SnowDepth(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: SnowDensity(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: SoilTemp(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: SnowAge( :)
REAL, INTENT(IN OUT), ALLOCATABLE :: AlbSoil(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: reducedLAIdue2snow(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: HeightAboveSnow(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: LAI_pft_cbl(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: HGT_pft_cbl(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: coszen(:)
!these local to CABLE and can be flushed every timestep
REAL, INTENT(IN OUT), ALLOCATABLE :: ExtCoeff_beam(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: ExtCoeff_dif(:)
REAL, INTENT(IN OUT), ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: RadFbeam(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: RadAlbedo(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: AlbSnow(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: c1(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: rhoch(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: xk(:,:)
REAL, INTENT(IN OUT), ALLOCATABLE :: SW_down(:,:)        ! dummy
INTEGER, INTENT(IN OUT), ALLOCATABLE :: metDoY(:)  ! pass DoY from current_time

IF ( ALLOCATED(EffSurfRefl_dif)     ) DEALLOCATE ( EffSurfRefl_dif )
IF ( ALLOCATED(EffSurfRefl_beam)    ) DEALLOCATE ( EffSurfRefl_beam )
IF ( ALLOCATED(SnowDepth)           ) DEALLOCATE( SnowDepth )
IF ( ALLOCATED(SnowDensity)         ) DEALLOCATE( SnowDensity )
IF ( ALLOCATED(SoilTemp)            ) DEALLOCATE( SoilTemp )
IF ( ALLOCATED(SnowAge)             ) DEALLOCATE( SnowAge )
IF ( ALLOCATED(AlbSoil)             ) DEALLOCATE( AlbSoil )
IF ( ALLOCATED(reducedLAIdue2snow)  ) DEALLOCATE( reducedLAIdue2snow )
IF ( ALLOCATED(HeightAboveSnow)     ) DEALLOCATE( HeightAboveSnow )
IF ( ALLOCATED(LAI_pft_cbl)         ) DEALLOCATE( LAI_pft_cbl )
IF ( ALLOCATED(HGT_pft_cbl)         ) DEALLOCATE( HGT_pft_cbl )
IF ( ALLOCATED(coszen)              ) DEALLOCATE( coszen )
IF ( ALLOCATED (ExtCoeff_beam)      ) DEALLOCATE (ExtCoeff_beam )
IF ( ALLOCATED (ExtCoeff_dif)       ) DEALLOCATE (ExtCoeff_dif )
IF ( ALLOCATED (EffExtCoeff_beam)   ) DEALLOCATE (EffExtCoeff_beam )
IF ( ALLOCATED (EffExtCoeff_dif)    ) DEALLOCATE (EffExtCoeff_dif )
IF ( ALLOCATED (CanopyTransmit_dif) ) DEALLOCATE (CanopyTransmit_dif )
IF ( ALLOCATED (CanopyTransmit_beam)) DEALLOCATE (CanopyTransmit_beam )
IF ( ALLOCATED (CanopyRefl_dif)     ) DEALLOCATE (CanopyRefl_dif )
IF ( ALLOCATED (CanopyRefl_beam)    ) DEALLOCATE (CanopyRefl_beam )
IF ( ALLOCATED (RadFbeam)           ) DEALLOCATE (RadFbeam )
IF ( ALLOCATED (RadAlbedo)          ) DEALLOCATE (RadAlbedo )
IF ( ALLOCATED (AlbSnow)            ) DEALLOCATE (AlbSnow )
IF ( ALLOCATED (c1)                 ) DEALLOCATE (c1 )
IF ( ALLOCATED (rhoch)              ) DEALLOCATE (rhoch )
IF ( ALLOCATED (xk)                 ) DEALLOCATE (xk )
IF ( ALLOCATED (SW_down)            ) DEALLOCATE (SW_down )
IF ( ALLOCATED (MetDoY)             ) DEALLOCATE (MetDoY )

RETURN
END SUBROUTINE flush_local_vars

END MODULE alloc_rad_albedo_vars_mod
