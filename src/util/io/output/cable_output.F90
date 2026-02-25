module cable_output_prototype_v2_mod ! TODO(Sean): rename to cable_output_mod

  use cable_output_core_mod, only: cable_output_mod_init
  use cable_output_core_mod, only: cable_output_mod_end
  use cable_output_core_mod, only: cable_output_register_output_variables
  use cable_output_core_mod, only: cable_output_profiles_init
  use cable_output_core_mod, only: cable_output_update
  use cable_output_core_mod, only: cable_output_write
  use cable_output_core_mod, only: cable_output_write_parameters
  use cable_output_core_mod, only: cable_output_write_restart

  use cable_output_types_mod, only: cable_output_attribute_t
  use cable_output_types_mod, only: cable_output_variable_t

  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PATCH
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOIL
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SNOW
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_RAD
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PLANTCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOILCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_LAND_GLOBAL

  implicit none
  public

end module
