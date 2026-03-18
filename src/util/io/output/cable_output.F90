! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_mod
  !! This module provides the public interface to the CABLE output system.

  use cable_output_core_mod, only: cable_output_mod_init
  use cable_output_core_mod, only: cable_output_mod_end
  use cable_output_core_mod, only: cable_output_get_dimension
  use cable_output_core_mod, only: cable_output_register_output_variables
  use cable_output_core_mod, only: cable_output_profiles_init
  use cable_output_core_mod, only: cable_output_update
  use cable_output_core_mod, only: cable_output_write
  use cable_output_core_mod, only: cable_output_write_parameters
  use cable_output_core_mod, only: cable_output_write_restart

  use cable_output_types_mod, only: cable_output_attribute_t
  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_dim_t

  implicit none
  public

end module
