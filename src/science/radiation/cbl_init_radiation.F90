!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) can be found
! at https://github.com/CABLE-LSM/CABLE/blob/main/
!
! ==============================================================================

MODULE cbl_init_radiation_module

  !*# Overview
  !
  ! This MODULE initialise the radiation parameters.
  !
  ! The MODULE contains one public subroutine [[init_radiation]] and several
  ! private subroutines.

IMPLICIT NONE

PUBLIC :: init_radiation
PUBLIC :: Common_InitRad_Scalings
PRIVATE

CONTAINS

SUBROUTINE init_radiation( ExtCoeff_beam, ExtCoeff_dif,                        &
                           EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,        &
                           c1, rhoch, xk, mp,nrb,                              &
                           Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,   &
                           cbl_standalone, jls_standalone, jls_radiation,      &
                           subr_name, veg_mask, VegXfang, VegTaul, VegRefl,    &
                           coszen, metDoY, SW_down, reducedLAIdue2snow )
! Description:
!   Computes various extinction coefficients for different radiations (visible,
!   diffuse) and various quantities for black leaves.

IMPLICIT NONE

!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"
INTEGER, INTENT(IN) :: nrb                    ! # of rad. bands VIS,NIR(,LW)

!returned variables
REAL, INTENT(OUT) :: ExtCoeff_beam(mp)        ! Extinction co-eff. - Beam SW
REAL, INTENT(OUT) :: ExtCoeff_dif(mp)         ! Extinction co-eff. - Diffuse SW
REAL, INTENT(OUT) :: EffExtCoeff_beam(mp,nrb) ! Effective Ext. co-eff. - Beam SW
REAL, INTENT(OUT) :: EffExtCoeff_dif(mp,nrb)  ! Effective Ext. co-eff. - Dif. SW
REAL, INTENT(OUT) :: RadFbeam(mp,nrb)         ! Beam Fraction of SW [rad%fbeam]
REAL, INTENT(OUT) :: c1(mp,nrb)
REAL, INTENT(OUT) :: rhoch(mp,nrb)            ! REFLection black horiz leaves
REAL, INTENT(OUT) :: xk(mp,nrb)               ! EXT. coef beam SW - black leaves

!constants
REAL, INTENT(IN) :: Clai_thresh               ! min. LAI for vegetated
REAL, INTENT(IN) :: Ccoszen_tols              ! min. cosine zenith for SUNLIT
REAL, INTENT(IN) :: Cgauss_w(nrb)
REAL, INTENT(IN) :: Cpi
REAL, INTENT(IN) :: Cpi180
!what model am i in
LOGICAL, INTENT(IN) :: cbl_standalone         ! TRUE =  cable_standalone
LOGICAL, INTENT(IN) :: jls_standalone         ! TRUE  =  jules_standalone
LOGICAL, INTENT(IN) :: jls_radiation          ! TRUE  =  radiation pathway
CHARACTER(LEN=*), INTENT(IN) :: subr_name     ! where am i called  from
!masks
LOGICAL, INTENT(IN) :: veg_mask(mp)           ! TRUE = vegetated

!vegetation parameters input via namelist
REAL, INTENT(IN) :: VegXfang(mp)              ! Leaf Angle
REAL, INTENT(IN) :: VegTaul(mp,nrb)           ! Leaf Transmisivity
REAL, INTENT(IN) :: VegRefl(mp,nrb)           ! Leaf Reflectivity

REAL, INTENT(IN) :: reducedLAIdue2snow(mp)    ! Effective LAI given snow
REAL, INTENT(IN) :: coszen(mp)                ! cosine zenith angle
REAL, INTENT(IN) :: SW_down(mp,nrb)           ! Downward SW [formerly met%fsd]
INTEGER, INTENT(IN) :: metDoY(mp)             ! Day of the Year [met%doy]

!local_vars - common scaling co-efficients used throughout init_radiation
REAL :: Ccoszen_tols_huge             ! Set constant to avoid Hard-wiring limits
REAL :: Ccoszen_tols_tiny             ! that occur in more than one place
REAL :: xvlai2(mp,nrb)                ! 2D vlai
REAL :: xphi1(mp)                     ! leaf angle parmameter 1
REAL :: xphi2(mp)                     ! leaf angle parmameter 2

!Null Initializations
ExtCoeff_beam(:) = 0.0
ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0
EffExtCoeff_dif(:,:) = 0.0

c1(:,:) = 0.0
rhoch(:,:) = 0.0
xk(:,:) = 0.0

! Compute common scaling co-efficients used throughout init_radiation
CALL Common_InitRad_Scalings( xphi1, xphi2, xk, xvlai2, c1, rhoch,             &
                              mp, nrb, Cpi180,cLAI_thresh, veg_mask,           &
                              reducedLAIdue2snow, VegXfang, VegTaul, VegRefl)

!Limiting Initializations for stability
Ccoszen_tols_huge = Ccoszen_tols * 1e2
Ccoszen_tols_tiny = Ccoszen_tols * 1e-2

! Define Raw extinction co-efficients for direct beam/diffuse radiation
! Largely parametrized per PFT. Does depend on zenith angle and effective LAI
! [Formerly rad%extkb, rad%extkd]
CALL ExtinctionCoeff( ExtCoeff_beam, ExtCoeff_dif, mp, nrb,                    &
                      CGauss_w,Ccoszen_tols_tiny, reducedLAIdue2snow,          &
                      veg_mask, cLAI_thresh, coszen, xphi1, xphi2, xk, xvlai2 )

! Define effective Extinction co-efficient for direct beam/diffuse radiation
! Extincion Co-eff defined by parametrized leaf reflect(transmit)ance - used in
! canopy transmitance calculations (cbl_albeo)
! [Formerly rad%extkbm, rad%extkdm ]
CALL EffectiveExtinctCoeffs( EffExtCoeff_beam, EffExtCoeff_dif,                &
                             mp, nrb, veg_mask, ExtCoeff_beam, ExtCoeff_dif, c1)

! Offline/standalone forcing gives us total downward Shortwave. We have
! previosuly, arbitratily split this into NIR/VIS (50/50). We use
! Spitter function to split these bands into direct beam and diffuse components
IF ( cbl_standalone .OR. jls_standalone .AND. .NOT. jls_radiation )            &
  CALL BeamFraction( RadFbeam, mp, nrb, Cpi, Ccoszen_tols_huge, metDoy,        &
                     coszen, SW_down )
RETURN
END SUBROUTINE init_radiation

!===============================================================================

SUBROUTINE Common_InitRad_Scalings( xphi1, xphi2, xk, xvlai2, c1, rhoch,       &
                            mp, nrb, Cpi180,cLAI_thresh, veg_mask,             &
                            reducedLAIdue2snow,                                &
                            VegXfang, VegTaul, VegRefl)
!* Calculates the extinction coefficients for black leaves. It returns:
!
! * the extinction coefficients for three values of the zenith angle
!  to be used to calculate the real extinction coefficient for the
!  diffuse radiation (Gauss quadrature)
!
! * the \(\phi\) coefficients from Sellers 1985 to calculate the
!  real extinction coefficient for direct beam radaiation in
!  the subroutine [[ExtinctionCoeff_beam]]
!
!### Equations
!
! This subroutine is using the equation B6 from
! [Wang and Leuning 1998](https://www.sciencedirect.com/science/article/abs/pii/S0168192398000616)
! and equation 13 from [Sellers 1985](https://www.tandfonline.com/doi/pdf/10.1080/01431168508948283?needAccess=true):
!
! \[ k_b = \frac{\phi_1}{\cos\theta} + \phi_2 \]
! \[ \phi_1 = 0.5 - 0.633 \chi_L - 0.33 \chi_L^2 \]
! \[ \phi_2 = 0.877 (1.0 - 2.0 \phi_1) \]

!subrs
USE cbl_rhoch_module,   ONLY: calc_rhoch
IMPLICIT NONE
! model dimensions
INTEGER, INTENT(IN) :: mp
  !! Number of tiles
INTEGER, INTENT(IN) :: nrb
  !! Number of radiation bands

REAL, INTENT(OUT) :: xk(mp,nrb)   ! extinct. coef.for beam rad. and black leaves
  !! Extinction coefficients for black leaves at 3 different zenith angles
REAL, INTENT(OUT) :: c1(mp,nrb)
REAL, INTENT(OUT) :: rhoch(mp,nrb)
REAL, INTENT(OUT) :: xphi1(mp)    ! leaf angle parmameter 1
  !! Leaf angle parameter defined by Sellers 1985 \(\phi_1\) 
REAL, INTENT(OUT) :: xphi2(mp)    ! leaf angle parmameter 2
  !! Leaf angle parameter defined by Sellers 1985 \(\phi_2\) 
REAL, INTENT(OUT) :: xvlai2(mp,nrb)  ! 2D vlai
  !! LAI spread over the 3 different zenith angles

REAL, INTENT(IN) :: Cpi180
  !! \(\pi\)

REAL, INTENT(IN) :: cLAI_thresh
  !! Threshold for the LAI under which a tile is considered unvegetated
REAL, INTENT(IN) :: reducedLAIdue2snow(mp)
  !! LAI after the effect of snow
REAL, INTENT(IN) :: VegXfang(mp)
  !! Parameter \(\chi\) in Sellers 1985
REAL, INTENT(IN) :: VegTaul(mp,nrb)
REAL, INTENT(IN) :: VegRefl(mp,nrb)
LOGICAL, INTENT(IN) :: veg_mask(mp)
  !! Mask indicating the presence of vegetation on a tile

!local vars
REAL :: cos3(nrb)      ! cos(15 45 75 degrees)

cos3 = COS(cpi180 * [ 15.0, 45.0, 75.0 ])

xphi1 = 0.0
xphi2 = 0.0
xvlai2 = 0.0
! See Sellers 1985, eq.13 (leaf angle parameters):
WHERE ( veg_mask )
  xphi1 = 0.5 - VegXfang * (0.633 + 0.33 * VegXfang)
  xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
END WHERE

! 2 dimensional LAI
xvlai2 = SPREAD(reducedLAIdue2snow, 2, 3)

! Extinction coefficient for beam radiation and black leaves;
! eq. B6, Wang and Leuning, 1998
WHERE (xvlai2 > cLAI_THRESH) ! vegetated
  xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
ELSE WHERE ! i.e. bare soil
  xk = 0.0
END WHERE

CALL calc_rhoch( c1,rhoch, mp, nrb, VegTaul, VegRefl )

RETURN
END SUBROUTINE Common_InitRad_Scalings

!===============================================================================

SUBROUTINE ExtinctionCoeff( ExtCoeff_beam, ExtCoeff_dif, mp, nrb, CGauss_w,    &
                            Ccoszen_tols_tiny, reducedLAIdue2snow, veg_mask,   &
                            cLAI_thresh, coszen, xphi1, xphi2, xk, xvlai2 )
! Description:
!   Define Raw extinction co-efficients for direct beam/diffuse radiation
!   Largely parametrized per PFT. Does depend on zenith angle and effective LAI

IMPLICIT NONE
!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"
INTEGER, INTENT(IN) :: nrb                    ! # of rad. bands VIS,NIR(,LW)

REAL, INTENT(OUT) :: ExtCoeff_beam(mp)        !extinction co-eff RETURNED
REAL, INTENT(OUT) :: ExtCoeff_dif(mp)         !extinction co-eff RETURNED

LOGICAL, INTENT(IN) :: veg_mask(mp)           !vegetated mask based on a min LAI
REAL, INTENT(IN)  :: Cgauss_w(nrb)
REAL, INTENT(IN)  :: Ccoszen_tols_tiny      ! min threshold of zenith for SUNLIT
REAL, INTENT(IN)  :: cLAI_thresh
REAL, INTENT(IN)  :: coszen(mp)
REAL, INTENT(IN)  :: reducedLAIdue2snow(mp)
REAL, INTENT(IN)  :: xphi1(mp)
REAL, INTENT(IN)  :: xphi2(mp)
REAL, INTENT(IN)  :: xvlai2(mp,nrb)  ! 2D vlai
REAL, INTENT(IN)  :: xk(mp,nrb)      ! ext. coef.for beam rad. and black leaves

! initialization for bare soil
ExtCoeff_beam = 0.5
ExtCoeff_dif = 0.7

! Diffuse component of Extiction  Coeff
WHERE ( veg_mask ) ! vegetated
  ! Extinction coefficient for diffuse radiation for black leaves:
  ExtCoeff_dif = -LOG( SUM(                                                    &
                            SPREAD( cgauss_w, 1, mp )                          &
                            * EXP( -xk * xvlai2 ), 2)                          &
                     )                                                         &
                     / reducedLAIdue2snow
END WHERE


! Direct beam component of Extiction  Coeff

! SW beam extinction coefficient ("black" leaves, extinction neglects
! leaf SW transmittance and REFLectance):
WHERE ( veg_mask .AND. coszen > Ccoszen_tols_tiny )
  ExtCoeff_beam = xphi1 / Coszen + xphi2
END WHERE

! higher value precludes sunlit leaves at night. affects
! nighttime evaporation - Ticket #90
WHERE ( coszen <  Ccoszen_tols_tiny )
  ExtCoeff_beam = 1.0e5
END WHERE

! Seems to be for stability only
WHERE ( ABS(ExtCoeff_beam - ExtCoeff_dif )  < 0.001 )
  ExtCoeff_beam = ExtCoeff_dif + 0.001
END WHERE

RETURN
END SUBROUTINE ExtinctionCoeff

!===============================================================================

SUBROUTINE EffectiveExtinctCoeffs( EffExtCoeff_beam, EffExtCoeff_dif, mp, nrb, &
                                   veg_mask, ExtCoeff_beam, ExtCoeff_dif, c1 )
! Description:
!   Define effective Extinction co-efficient for direct beam/diffuse radiation
!   Extincion Co-eff defined by parametrized leaf reflect(transmit)ance - used in
!   canopy transmitance calculations (cbl_albeo)

IMPLICIT NONE
!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"
INTEGER, INTENT(IN) :: nrb                    ! # of rad. bands VIS,NIR(,LW)

! Effective Extinction co-efficients
REAL, INTENT(OUT) :: EffExtCoeff_beam(mp,nrb)! Direct Beam component to SW
REAL, INTENT(OUT) :: EffExtCoeff_dif(mp,nrb) ! Diffuse component to SW

REAL :: c1(mp,nrb)
LOGICAL :: veg_mask(mp)  !mask -  vegetated
! "raw" Extinction co-efficients
REAL, INTENT(IN) :: ExtCoeff_beam(mp)       ! Direct Beam component to SW
REAL, INTENT(IN) :: ExtCoeff_dif(mp)        ! Diffuse component to SW

EffExtCoeff_beam = 0.0
EffExtCoeff_dif = 0.0

! we re-use basically the same function here for both components
CALL EffectiveExtinctCoeff( EffExtCoeff_beam, mp, ExtCoeff_beam, c1, veg_mask )

CALL EffectiveExtinctCoeff( EffExtCoeff_dif, mp, ExtCoeff_dif, c1 )

RETURN
END SUBROUTINE EffectiveExtinctCoeffs

!===============================================================================

SUBROUTINE EffectiveExtinctCoeff(Eff_ExtCoeff, mp, ExtCoeff, c1, mask )
! Description:
!   modified k diffuse(6.20)(for leaf scattering)

IMPLICIT NONE
!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"

REAL, INTENT(OUT) :: Eff_ExtCoeff(mp,2)

REAL, INTENT(IN) :: ExtCoeff(mp)
REAL, INTENT(IN) :: c1(mp,2)
LOGICAL, INTENT(IN), OPTIONAL :: mask(mp)
!local
INTEGER :: i, b

DO i = 1,mp
  DO b = 1, 2
    !IF mask is present we are doing the beam component then:
    IF ( PRESENT(mask)) THEN
      !then ONLY IF it is sunlit and vegetated -else default
      IF ( mask(i) ) THEN
        Eff_ExtCoeff(i,b) = ExtCoeff(i) * c1(i,b)
      END IF
    ELSE
      Eff_ExtCoeff(i,b) = ExtCoeff(i) * c1(i,b)
    END IF

  END DO
END DO

RETURN
END SUBROUTINE EffectiveExtinctCoeff

!===============================================================================

SUBROUTINE BeamFraction( RadFbeam, mp, nrb, Cpi,Ccoszen_tols_huge, metDoy,     &
                         coszen, SW_down )
! Description:
!   Split the total downward Shortwave into Near-infrared and visible components
!   Based on the Spitter function.

USE cbl_spitter_module, ONLY: Spitter
IMPLICIT NONE
!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"
INTEGER, INTENT(IN) :: nrb                    ! # of rad. bands VIS,NIR(,LW)

REAL, INTENT(IN OUT) :: RadFbeam(mp,nrb)       ! Beam Fraction of SW [rad%fbeam]

REAL, INTENT(IN) :: Cpi                       ! PI
REAL, INTENT(IN) :: Ccoszen_tols_huge         ! min. cos zenith for beam frac >0

INTEGER, INTENT(IN):: metDoY(mp)              ! Day of the Year [met%doy]
REAL, INTENT(IN):: coszen(mp)                 ! Day of the Year [met%coszen]
REAL, INTENT(IN) :: SW_down(mp,nrb)           ! Downward SW [met%fsd]


! Define beam fraction, fbeam:
RadFbeam(:,1) = spitter(mp, cpi, metDoy, coszen, SW_down(:,1))
RadfBeam(:,2) = spitter(mp, cpi, metDoy, coszen, SW_down(:,2))

! coszen is set during met data read in.
WHERE (coszen < Ccoszen_tols_huge )
  RadFbeam(:,1) = 0.0
  RadFbeam(:,2) = 0.0
END WHERE

RETURN
END SUBROUTINE BeamFraction

END MODULE cbl_init_radiation_module
