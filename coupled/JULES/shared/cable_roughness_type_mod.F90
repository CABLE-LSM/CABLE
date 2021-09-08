MODULE cable_roughness_type_mod

IMPLICIT NONE

  ! Roughness variables:
TYPE roughness_type

  REAL, DIMENSION(:), POINTER ::                                              &
       disp,    & ! zero-plane displacement
       hruff,   & ! canopy height above snow level
       hruff_grmx,&! max ht of canopy from tiles on same grid
       rt0us,   & ! eq. 3.54, SCAM manual (CSIRO tech report 132)
       rt1usa,  & ! resistance from disp to hruf
       rt1usb,  & ! resist fr hruf to zruffs (zref if zref<zruffs)
       rt1,     & ! 1/aerodynamic conductance
       za_uv,   & ! level of lowest atmospheric model layer
       za_tq,   & ! level of lowest atmospheric model layer
       z0m,     & ! roughness length
       zref_uv, & ! Reference height for met forcing
       zref_tq, & ! Reference height for met forcing
       zruffs,  & ! SCALAR Roughness sublayer depth (ground=origin)
       z0soilsn,& ! roughness length of bare soil surface
       z0soil     ! roughness length of bare soil surface

  ! "coexp": coefficient in exponential in-canopy wind profile
  ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
  ! canopy and roughness-sublayer U(z) at z=h
  REAL, DIMENSION(:), POINTER ::                                              &
       coexp ! Extinction coef for wind profile in canopy

  ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
  REAL, DIMENSION(:), POINTER ::                                              &
       usuh ! Friction velocity/windspeed at canopy height

  REAL, DIMENSION(:), POINTER ::                                              &
       term2, term3, term5, term6, term6a ! for aerodyn resist. calc.



END TYPE roughness_type

!Instantiation:
TYPE(roughness_type) :: rough_cbl

END MODULE cable_roughness_type_mod
