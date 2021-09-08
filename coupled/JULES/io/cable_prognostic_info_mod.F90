MODULE cable_prognostic_info_mod

USE cable_types_mod

IMPLICIT NONE

! Tiled soil prognostics to be initialized from IO
REAL, ALLOCATABLE, PUBLIC ::                                                  &
  SoilTemp_CABLE(:,:,:),                                                      &
  SoilMoisture_CABLE(:,:,:),                                                  &
  FrozenSoilFrac_CABLE(:,:,:),                                                &
  SnowDepth_CABLE(:,:,:),                                                     &
  SnowMass_CABLE(:,:,:),                                                      &
  SnowTemp_CABLE(:,:,:),                                                      &
  SnowDensity_CABLE(:,:,:),                                                   &
  ThreeLayerSnowFlag_CABLE(:,:),                                              &
  OneLyrSnowDensity_CABLE(:,:),                                               &
  SnowAge_CABLE(:,:),                                                         &        
  snowOsurft(:,:)

END MODULE cable_prognostic_info_mod

