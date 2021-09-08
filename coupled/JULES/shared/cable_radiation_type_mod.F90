MODULE cable_radiation_type_mod

IMPLICIT NONE

  ! Radiation variables:
TYPE radiation_type

  REAL, DIMENSION(:), POINTER   ::                                            &
       transb,  & ! fraction SW beam tranmitted through canopy
       albedo_T,& ! canopy+soil albedo for VIS+NIR
       longitude,&! longitude
       workp1,  & ! absorbed short-wave radiation for soil
       workp2,  & ! absorbed short-wave radiation for soil
       workp3,  & ! absorbed short-wave radiation for soil
       extkb,   & ! beam radiation extinction coeff
       extkd2,  & ! diffuse 2D radiation extinction coeff
       extkd,   & ! diffuse radiation extinction coeff (-)
       flws,    & ! soil long-wave radiation
       latitude,& ! latitude
       lwabv,   & ! long wave absorbed by vegetation
       qssabs,  & ! absorbed short-wave radiation for soil
       transd,  & ! frac SW diffuse transmitted through canopy
       trad,    & !  radiative temperature (soil and veg)
       otrad      ! radiative temperature on previous timestep (ACCESS)

  REAL, DIMENSION(:,:), POINTER  ::                                           &
       fvlai,   & ! leaf area index of big leaf
       rhocdf,  & ! canopy diffuse reflectance (-)
       rniso,   & ! sum(rad%qcan, 3) total abs by canopy (W/m2)
       scalex,  & ! scaling PARAMETER for big leaf
       albedo,  & ! canopy+soil albedo
       reffdf,  & ! effective conopy diffuse reflectance
       reffbm,  & ! effective conopy beam reflectance
       extkbm,  & ! modified k beam(6.20)(for leaf scattering)
       extkdm,  & ! modified k diffuse(6.20)(for leaf scattering)
       fbeam,   & ! beam fraction
       cexpkbm, & ! canopy beam transmittance
       cexpkdm, & ! canopy diffuse transmittance
       rhocbm,  & ! modified canopy beam reflectance(6.21)
       gradis     ! radiative conductance

  REAL, DIMENSION(:,:,:), POINTER ::                                          &
       qcan ! absorbed radiation for canopy (W/m^2)


END TYPE radiation_type

!Instantiation:
TYPE(radiation_type) :: rad_cbl

END MODULE cable_radiation_type_mod
