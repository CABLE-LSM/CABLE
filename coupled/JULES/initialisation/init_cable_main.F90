! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE init_cable_mod
CONTAINS
SUBROUTINE init_cable()

USE def_cable_grid_mod,       ONLY: def_cable_grid
USE init_cable_params_mod,    ONLY: init_cable_params
USE init_cable_state_mod,     ONLY: init_cable_state
USE init_cable_runtime_opts_mod, ONLY: init_cable_runtime_options
USE cable_types_mod,          ONLY: mp

IMPLICIT NONE

! Declare/allocate/initialize - params,progs,state vars for CABLE
CALL init_cable_runtime_options() ! allocates/initializes veg/soil from vegin/soilin 
CALL def_cable_grid(mp)   ! init. #active tiles,define mask for active tiles
CALL init_cable_params(mp) ! allocates/initializes veg/soil from vegin/soilin 
                        ! params already read from namelist via init_params.inc
CALL init_cable_state(mp)

!CALL set_cable_progs() ! alloc/inits - BUT called from init_ancillaries.inc

RETURN

END SUBROUTINE init_cable
END MODULE init_cable_mod
