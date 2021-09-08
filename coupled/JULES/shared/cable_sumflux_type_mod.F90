MODULE cable_sum_flux_type_mod

IMPLICIT NONE
  
TYPE sum_flux_type

  REAL, DIMENSION(:), POINTER ::                                              &
       sumpn,   & ! sum of canopy photosynthesis (g C m-2)
       sumrp,   & ! sum of plant respiration (g C m-2)
       sumrpw,  & ! sum of plant respiration (g C m-2)
       sumrpr,  & ! sum of plant respiration (g C m-2)
       sumrs,   & ! sum of soil respiration (g C m-2)
       sumrd,   & ! sum of daytime respiration (g C m-2)
       dsumpn,  & ! daily sumpn
       dsumrp,  & ! daily sumrp
       dsumrs,  & ! daily sumrs
       dsumrd,  & ! daily sumrd
       sumxrp,  & ! sum plant resp. modifier
       sumxrs     ! sum soil resp. modifier

END TYPE sum_flux_type

!Instantiation:
TYPE(sum_flux_type) :: sum_flux_cbl

END MODULE cable_sum_flux_type_mod

