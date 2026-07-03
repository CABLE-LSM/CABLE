! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

! TODO(Sean): The preprocessor define ENFORCE_SINGLE_PRECISION is enabled
! temporarily to restore bitwise reproducibility with the previous output module
! which enforces writing both double and single precision data as single
! precision.
#define ENFORCE_SINGLE_PRECISION

submodule (cable_output_mod:cable_output_common_smod) cable_output_decomp_smod
  !* Implementation of procedures for creating and managing I/O decompositions
  ! for the CABLE output system.

  use cable_error_handler_mod, only: cable_abort

  use cable_array_utils_mod, only: array_eq

  use cable_def_types_mod, only: mp
  use cable_def_types_mod, only: mp_global
  use cable_def_types_mod, only: mland
  use cable_def_types_mod, only: mland_global
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: ncp

  use cable_io_vars_module, only: xdimsize
  use cable_io_vars_module, only: ydimsize
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: land_x, land_y
  use cable_io_vars_module, only: landpt
  use cable_io_vars_module, only: land_decomp_start
  use cable_io_vars_module, only: patch_decomp_start

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_x_y
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_x_y_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_land
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_land_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_patch

  implicit none

  type :: cable_output_decomp_t
    !* Data structure for holding the I/O decompositions for each output grid
    ! type and variable type.
    !
    ! Each component represents the in-memory shape of the data being written
    ! (not to be confused with the shape of the netCDF variable on disk).
    class(cable_netcdf_decomp_t), allocatable :: land
      !! I/O decomposition for data with shape `[mland]`
    class(cable_netcdf_decomp_t), allocatable :: land_soil
      !! I/O decomposition for data with shape `[mland, ms]`
    class(cable_netcdf_decomp_t), allocatable :: land_snow
      !! I/O decomposition for data with shape `[mland, msn]`
    class(cable_netcdf_decomp_t), allocatable :: land_rad
      !! I/O decomposition for data with shape `[mland, nrb]`
    class(cable_netcdf_decomp_t), allocatable :: land_plantcarbon
      !! I/O decomposition for data with shape `[mland, ncp]`
    class(cable_netcdf_decomp_t), allocatable :: land_soilcarbon
      !! I/O decomposition for data with shape `[mland, ncs]`
    class(cable_netcdf_decomp_t), allocatable :: patch
      !! I/O decomposition for data with shape `[mp]`
    class(cable_netcdf_decomp_t), allocatable :: patch_soil
      !! I/O decomposition for data with shape `[mp, ms]`
    class(cable_netcdf_decomp_t), allocatable :: patch_snow
      !! I/O decomposition for data with shape `[mp, msn]`
    class(cable_netcdf_decomp_t), allocatable :: patch_rad
      !! I/O decomposition for data with shape `[mp, nrb]`
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon
      !! I/O decomposition for data with shape `[mp, ncp]`
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon
      !! I/O decomposition for data with shape `[mp, ncs]`
  end type

  type(cable_output_decomp_t), target :: decomps_grid_x_y_int32
    !! Decompositions for writing to an x-y grid.
  type(cable_output_decomp_t), target :: decomps_grid_x_y_real32
    !! Decompositions for writing to an x-y grid.
  type(cable_output_decomp_t), target :: decomps_grid_x_y_real64
    !! Decompositions for writing to an x-y grid.
  type(cable_output_decomp_t), target :: decomps_grid_land_int32
    !! Decompositions for writing to a land grid.
  type(cable_output_decomp_t), target :: decomps_grid_land_real32
    !! Decompositions for writing to a land grid.
  type(cable_output_decomp_t), target :: decomps_grid_land_real64
    !! Decompositions for writing to a land grid.
  type(cable_output_decomp_t), target :: decomps_grid_restart_int32
    !! Decompositions for writing to a restart grid.
  type(cable_output_decomp_t), target :: decomps_grid_restart_real32
    !! Decompositions for writing to a restart grid.
  type(cable_output_decomp_t), target :: decomps_grid_restart_real64
    !! Decompositions for writing to a restart grid.

contains

  module subroutine cable_output_decomp_init()
    !! Intialises I/O decompositions used in the output system.

    decomps_grid_x_y_int32%land              = io_decomp_land_to_x_y(land_x, land_y, [mland],      [xdimsize, ydimsize],      CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%land_soil         = io_decomp_land_to_x_y(land_x, land_y, [mland, ms],  [xdimsize, ydimsize, ms],  CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%land_snow         = io_decomp_land_to_x_y(land_x, land_y, [mland, msn], [xdimsize, ydimsize, msn], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%land_rad          = io_decomp_land_to_x_y(land_x, land_y, [mland, nrb], [xdimsize, ydimsize, nrb], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%land_plantcarbon  = io_decomp_land_to_x_y(land_x, land_y, [mland, ncp], [xdimsize, ydimsize, ncp], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%land_soilcarbon   = io_decomp_land_to_x_y(land_x, land_y, [mland, ncs], [xdimsize, ydimsize, ncs], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%patch             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp],      [xdimsize, ydimsize, max_vegpatches],      CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%patch_soil        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [xdimsize, ydimsize, max_vegpatches, ms],  CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%patch_snow        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [xdimsize, ydimsize, max_vegpatches, msn], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%patch_rad         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [xdimsize, ydimsize, max_vegpatches, nrb], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%patch_plantcarbon = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [xdimsize, ydimsize, max_vegpatches, ncp], CABLE_NETCDF_INT)
    decomps_grid_x_y_int32%patch_soilcarbon  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [xdimsize, ydimsize, max_vegpatches, ncs], CABLE_NETCDF_INT)

    decomps_grid_land_int32%land              = io_decomp_land_to_land(land_decomp_start, [mland],      [mland_global],      CABLE_NETCDF_INT)
    decomps_grid_land_int32%land_soil         = io_decomp_land_to_land(land_decomp_start, [mland, ms],  [mland_global, ms],  CABLE_NETCDF_INT)
    decomps_grid_land_int32%land_snow         = io_decomp_land_to_land(land_decomp_start, [mland, msn], [mland_global, msn], CABLE_NETCDF_INT)
    decomps_grid_land_int32%land_rad          = io_decomp_land_to_land(land_decomp_start, [mland, nrb], [mland_global, nrb], CABLE_NETCDF_INT)
    decomps_grid_land_int32%land_plantcarbon  = io_decomp_land_to_land(land_decomp_start, [mland, ncp], [mland_global, ncp], CABLE_NETCDF_INT)
    decomps_grid_land_int32%land_soilcarbon   = io_decomp_land_to_land(land_decomp_start, [mland, ncs], [mland_global, ncs], CABLE_NETCDF_INT)
    decomps_grid_land_int32%patch             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp],      [mland_global, max_vegpatches],      CABLE_NETCDF_INT)
    decomps_grid_land_int32%patch_soil        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [mland_global, max_vegpatches, ms],  CABLE_NETCDF_INT)
    decomps_grid_land_int32%patch_snow        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [mland_global, max_vegpatches, msn], CABLE_NETCDF_INT)
    decomps_grid_land_int32%patch_rad         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [mland_global, max_vegpatches, nrb], CABLE_NETCDF_INT)
    decomps_grid_land_int32%patch_plantcarbon = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [mland_global, max_vegpatches, ncp], CABLE_NETCDF_INT)
    decomps_grid_land_int32%patch_soilcarbon  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [mland_global, max_vegpatches, ncs], CABLE_NETCDF_INT)

    decomps_grid_restart_int32%patch             = io_decomp_patch_to_patch(patch_decomp_start, [mp],      [mp_global],      CABLE_NETCDF_INT)
    decomps_grid_restart_int32%patch_soil        = io_decomp_patch_to_patch(patch_decomp_start, [mp, ms],  [mp_global, ms],  CABLE_NETCDF_INT)
    decomps_grid_restart_int32%patch_snow        = io_decomp_patch_to_patch(patch_decomp_start, [mp, msn], [mp_global, msn], CABLE_NETCDF_INT)
    decomps_grid_restart_int32%patch_rad         = io_decomp_patch_to_patch(patch_decomp_start, [mp, nrb], [mp_global, nrb], CABLE_NETCDF_INT)
    decomps_grid_restart_int32%patch_plantcarbon = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncp], [mp_global, ncp], CABLE_NETCDF_INT)
    decomps_grid_restart_int32%patch_soilcarbon  = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncs], [mp_global, ncs], CABLE_NETCDF_INT)

    decomps_grid_x_y_real32%land              = io_decomp_land_to_x_y(land_x, land_y, [mland],      [xdimsize, ydimsize],      CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%land_soil         = io_decomp_land_to_x_y(land_x, land_y, [mland, ms],  [xdimsize, ydimsize, ms],  CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%land_snow         = io_decomp_land_to_x_y(land_x, land_y, [mland, msn], [xdimsize, ydimsize, msn], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%land_rad          = io_decomp_land_to_x_y(land_x, land_y, [mland, nrb], [xdimsize, ydimsize, nrb], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%land_plantcarbon  = io_decomp_land_to_x_y(land_x, land_y, [mland, ncp], [xdimsize, ydimsize, ncp], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%land_soilcarbon   = io_decomp_land_to_x_y(land_x, land_y, [mland, ncs], [xdimsize, ydimsize, ncs], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%patch             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp],      [xdimsize, ydimsize, max_vegpatches],      CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%patch_soil        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [xdimsize, ydimsize, max_vegpatches, ms],  CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%patch_snow        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [xdimsize, ydimsize, max_vegpatches, msn], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%patch_rad         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [xdimsize, ydimsize, max_vegpatches, nrb], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%patch_plantcarbon = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [xdimsize, ydimsize, max_vegpatches, ncp], CABLE_NETCDF_FLOAT)
    decomps_grid_x_y_real32%patch_soilcarbon  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [xdimsize, ydimsize, max_vegpatches, ncs], CABLE_NETCDF_FLOAT)

    decomps_grid_land_real32%land              = io_decomp_land_to_land(land_decomp_start, [mland],      [mland_global],      CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%land_soil         = io_decomp_land_to_land(land_decomp_start, [mland, ms],  [mland_global, ms],  CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%land_snow         = io_decomp_land_to_land(land_decomp_start, [mland, msn], [mland_global, msn], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%land_rad          = io_decomp_land_to_land(land_decomp_start, [mland, nrb], [mland_global, nrb], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%land_plantcarbon  = io_decomp_land_to_land(land_decomp_start, [mland, ncp], [mland_global, ncp], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%land_soilcarbon   = io_decomp_land_to_land(land_decomp_start, [mland, ncs], [mland_global, ncs], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%patch             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp],      [mland_global, max_vegpatches],      CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%patch_soil        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [mland_global, max_vegpatches, ms],  CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%patch_snow        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [mland_global, max_vegpatches, msn], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%patch_rad         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [mland_global, max_vegpatches, nrb], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%patch_plantcarbon = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [mland_global, max_vegpatches, ncp], CABLE_NETCDF_FLOAT)
    decomps_grid_land_real32%patch_soilcarbon  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [mland_global, max_vegpatches, ncs], CABLE_NETCDF_FLOAT)

    decomps_grid_restart_real32%patch             = io_decomp_patch_to_patch(patch_decomp_start, [mp],      [mp_global],      CABLE_NETCDF_FLOAT)
    decomps_grid_restart_real32%patch_soil        = io_decomp_patch_to_patch(patch_decomp_start, [mp, ms],  [mp_global, ms],  CABLE_NETCDF_FLOAT)
    decomps_grid_restart_real32%patch_snow        = io_decomp_patch_to_patch(patch_decomp_start, [mp, msn], [mp_global, msn], CABLE_NETCDF_FLOAT)
    decomps_grid_restart_real32%patch_rad         = io_decomp_patch_to_patch(patch_decomp_start, [mp, nrb], [mp_global, nrb], CABLE_NETCDF_FLOAT)
    decomps_grid_restart_real32%patch_plantcarbon = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncp], [mp_global, ncp], CABLE_NETCDF_FLOAT)
    decomps_grid_restart_real32%patch_soilcarbon  = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncs], [mp_global, ncs], CABLE_NETCDF_FLOAT)

    decomps_grid_x_y_real64%land              = io_decomp_land_to_x_y(land_x, land_y, [mland],      [xdimsize, ydimsize],      CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%land_soil         = io_decomp_land_to_x_y(land_x, land_y, [mland, ms],  [xdimsize, ydimsize, ms],  CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%land_snow         = io_decomp_land_to_x_y(land_x, land_y, [mland, msn], [xdimsize, ydimsize, msn], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%land_rad          = io_decomp_land_to_x_y(land_x, land_y, [mland, nrb], [xdimsize, ydimsize, nrb], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%land_plantcarbon  = io_decomp_land_to_x_y(land_x, land_y, [mland, ncp], [xdimsize, ydimsize, ncp], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%land_soilcarbon   = io_decomp_land_to_x_y(land_x, land_y, [mland, ncs], [xdimsize, ydimsize, ncs], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%patch             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp],      [xdimsize, ydimsize, max_vegpatches],      CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%patch_soil        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [xdimsize, ydimsize, max_vegpatches, ms],  CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%patch_snow        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [xdimsize, ydimsize, max_vegpatches, msn], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%patch_rad         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [xdimsize, ydimsize, max_vegpatches, nrb], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%patch_plantcarbon = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [xdimsize, ydimsize, max_vegpatches, ncp], CABLE_NETCDF_DOUBLE)
    decomps_grid_x_y_real64%patch_soilcarbon  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [xdimsize, ydimsize, max_vegpatches, ncs], CABLE_NETCDF_DOUBLE)

    decomps_grid_land_real64%land              = io_decomp_land_to_land(land_decomp_start, [mland],      [mland_global],      CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%land_soil         = io_decomp_land_to_land(land_decomp_start, [mland, ms],  [mland_global, ms],  CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%land_snow         = io_decomp_land_to_land(land_decomp_start, [mland, msn], [mland_global, msn], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%land_rad          = io_decomp_land_to_land(land_decomp_start, [mland, nrb], [mland_global, nrb], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%land_plantcarbon  = io_decomp_land_to_land(land_decomp_start, [mland, ncp], [mland_global, ncp], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%land_soilcarbon   = io_decomp_land_to_land(land_decomp_start, [mland, ncs], [mland_global, ncs], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%patch             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp],      [mland_global, max_vegpatches],      CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%patch_soil        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [mland_global, max_vegpatches, ms],  CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%patch_snow        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [mland_global, max_vegpatches, msn], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%patch_rad         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [mland_global, max_vegpatches, nrb], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%patch_plantcarbon = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [mland_global, max_vegpatches, ncp], CABLE_NETCDF_DOUBLE)
    decomps_grid_land_real64%patch_soilcarbon  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [mland_global, max_vegpatches, ncs], CABLE_NETCDF_DOUBLE)

    decomps_grid_restart_real64%patch             = io_decomp_patch_to_patch(patch_decomp_start, [mp],      [mp_global],      CABLE_NETCDF_DOUBLE)
    decomps_grid_restart_real64%patch_soil        = io_decomp_patch_to_patch(patch_decomp_start, [mp, ms],  [mp_global, ms],  CABLE_NETCDF_DOUBLE)
    decomps_grid_restart_real64%patch_snow        = io_decomp_patch_to_patch(patch_decomp_start, [mp, msn], [mp_global, msn], CABLE_NETCDF_DOUBLE)
    decomps_grid_restart_real64%patch_rad         = io_decomp_patch_to_patch(patch_decomp_start, [mp, nrb], [mp_global, nrb], CABLE_NETCDF_DOUBLE)
    decomps_grid_restart_real64%patch_plantcarbon = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncp], [mp_global, ncp], CABLE_NETCDF_DOUBLE)
    decomps_grid_restart_real64%patch_soilcarbon  = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncs], [mp_global, ncs], CABLE_NETCDF_DOUBLE)

  end subroutine cable_output_decomp_init

  module subroutine cable_output_decomp_free()
    !! Deallocates I/O decompositions used in the output system.

    deallocate(decomps_grid_x_y_int32%land)
    deallocate(decomps_grid_x_y_int32%land_soil)
    deallocate(decomps_grid_x_y_int32%land_snow)
    deallocate(decomps_grid_x_y_int32%land_rad)
    deallocate(decomps_grid_x_y_int32%land_plantcarbon)
    deallocate(decomps_grid_x_y_int32%land_soilcarbon)
    deallocate(decomps_grid_x_y_int32%patch)
    deallocate(decomps_grid_x_y_int32%patch_soil)
    deallocate(decomps_grid_x_y_int32%patch_snow)
    deallocate(decomps_grid_x_y_int32%patch_rad)
    deallocate(decomps_grid_x_y_int32%patch_plantcarbon)
    deallocate(decomps_grid_x_y_int32%patch_soilcarbon)

    deallocate(decomps_grid_land_int32%land)
    deallocate(decomps_grid_land_int32%land_soil)
    deallocate(decomps_grid_land_int32%land_snow)
    deallocate(decomps_grid_land_int32%land_rad)
    deallocate(decomps_grid_land_int32%land_plantcarbon)
    deallocate(decomps_grid_land_int32%land_soilcarbon)
    deallocate(decomps_grid_land_int32%patch)
    deallocate(decomps_grid_land_int32%patch_soil)
    deallocate(decomps_grid_land_int32%patch_snow)
    deallocate(decomps_grid_land_int32%patch_rad)
    deallocate(decomps_grid_land_int32%patch_plantcarbon)
    deallocate(decomps_grid_land_int32%patch_soilcarbon)

    deallocate(decomps_grid_restart_int32%patch)
    deallocate(decomps_grid_restart_int32%patch_soil)
    deallocate(decomps_grid_restart_int32%patch_snow)
    deallocate(decomps_grid_restart_int32%patch_rad)
    deallocate(decomps_grid_restart_int32%patch_plantcarbon)
    deallocate(decomps_grid_restart_int32%patch_soilcarbon)

    deallocate(decomps_grid_x_y_real32%land)
    deallocate(decomps_grid_x_y_real32%land_soil)
    deallocate(decomps_grid_x_y_real32%land_snow)
    deallocate(decomps_grid_x_y_real32%land_rad)
    deallocate(decomps_grid_x_y_real32%land_plantcarbon)
    deallocate(decomps_grid_x_y_real32%land_soilcarbon)
    deallocate(decomps_grid_x_y_real32%patch)
    deallocate(decomps_grid_x_y_real32%patch_soil)
    deallocate(decomps_grid_x_y_real32%patch_snow)
    deallocate(decomps_grid_x_y_real32%patch_rad)
    deallocate(decomps_grid_x_y_real32%patch_plantcarbon)
    deallocate(decomps_grid_x_y_real32%patch_soilcarbon)

    deallocate(decomps_grid_land_real32%land)
    deallocate(decomps_grid_land_real32%land_soil)
    deallocate(decomps_grid_land_real32%land_snow)
    deallocate(decomps_grid_land_real32%land_rad)
    deallocate(decomps_grid_land_real32%land_plantcarbon)
    deallocate(decomps_grid_land_real32%land_soilcarbon)
    deallocate(decomps_grid_land_real32%patch)
    deallocate(decomps_grid_land_real32%patch_soil)
    deallocate(decomps_grid_land_real32%patch_snow)
    deallocate(decomps_grid_land_real32%patch_rad)
    deallocate(decomps_grid_land_real32%patch_plantcarbon)
    deallocate(decomps_grid_land_real32%patch_soilcarbon)

    deallocate(decomps_grid_restart_real32%patch)
    deallocate(decomps_grid_restart_real32%patch_soil)
    deallocate(decomps_grid_restart_real32%patch_snow)
    deallocate(decomps_grid_restart_real32%patch_rad)
    deallocate(decomps_grid_restart_real32%patch_plantcarbon)
    deallocate(decomps_grid_restart_real32%patch_soilcarbon)

    deallocate(decomps_grid_x_y_real64%land)
    deallocate(decomps_grid_x_y_real64%land_soil)
    deallocate(decomps_grid_x_y_real64%land_snow)
    deallocate(decomps_grid_x_y_real64%land_rad)
    deallocate(decomps_grid_x_y_real64%land_plantcarbon)
    deallocate(decomps_grid_x_y_real64%land_soilcarbon)
    deallocate(decomps_grid_x_y_real64%patch)
    deallocate(decomps_grid_x_y_real64%patch_soil)
    deallocate(decomps_grid_x_y_real64%patch_snow)
    deallocate(decomps_grid_x_y_real64%patch_rad)
    deallocate(decomps_grid_x_y_real64%patch_plantcarbon)
    deallocate(decomps_grid_x_y_real64%patch_soilcarbon)

    deallocate(decomps_grid_land_real64%land)
    deallocate(decomps_grid_land_real64%land_soil)
    deallocate(decomps_grid_land_real64%land_snow)
    deallocate(decomps_grid_land_real64%land_rad)
    deallocate(decomps_grid_land_real64%land_plantcarbon)
    deallocate(decomps_grid_land_real64%land_soilcarbon)
    deallocate(decomps_grid_land_real64%patch)
    deallocate(decomps_grid_land_real64%patch_soil)
    deallocate(decomps_grid_land_real64%patch_snow)
    deallocate(decomps_grid_land_real64%patch_rad)
    deallocate(decomps_grid_land_real64%patch_plantcarbon)
    deallocate(decomps_grid_land_real64%patch_soilcarbon)

    deallocate(decomps_grid_restart_real64%patch)
    deallocate(decomps_grid_restart_real64%patch_soil)
    deallocate(decomps_grid_restart_real64%patch_snow)
    deallocate(decomps_grid_restart_real64%patch_rad)
    deallocate(decomps_grid_restart_real64%patch_plantcarbon)
    deallocate(decomps_grid_restart_real64%patch_soilcarbon)

  end subroutine cable_output_decomp_free

  module subroutine cable_output_decomp_associate(output_stream, output_var, decomp)
    !* Associates an I/O decomposition pointer with the appropriate I/O
    ! decomposition, taking into account the output variable shape and type, and
    ! the output stream grid type.
    type(cable_output_stream_t), intent(in) :: output_stream
      !! The output stream for which to associate the decomposition.
    type(cable_output_variable_t), intent(in) :: output_var
      !! The output variable for which to associate the decomposition.
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
      !! The decomposition pointer to associate.
    type(cable_output_decomp_t), pointer :: output_decomp

    select case (output_stream%grid_type)
    case ("restart")
      call cable_output_decomp_associate_restart(output_var, decomp)
      return
    case ("mask")
      select case (output_var%aggregator%type())
      case ("int32")
        output_decomp => decomps_grid_x_y_int32
      case ("real32")
        output_decomp => decomps_grid_x_y_real32
      case ("real64")
#ifdef ENFORCE_SINGLE_PRECISION
        output_decomp => decomps_grid_x_y_real32
#else
        output_decomp => decomps_grid_x_y_real64
#endif
      case default
        call cable_abort("Unexpected data type for output variable " // output_var%field_name, __FILE__, __LINE__)
      end select
    case ("land")
      select case (output_var%aggregator%type())
      case ("int32")
        output_decomp => decomps_grid_land_int32
      case ("real32")
        output_decomp => decomps_grid_land_real32
      case ("real64")
#ifdef ENFORCE_SINGLE_PRECISION
        output_decomp => decomps_grid_land_real32
#else
        output_decomp => decomps_grid_land_real64
#endif
      case default
        call cable_abort("Unexpected data type for output variable " // output_var%field_name, __FILE__, __LINE__)
      end select
    case default
      call cable_abort("Unexpected grid type for output profile " // output_stream%file_name, __FILE__, __LINE__)
    end select

    if (array_eq(output_var%data_shape(:)%size(), [mp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch
      else
        decomp => output_decomp%land
      end if
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ms])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil
      else
        decomp => output_decomp%land_soil
      end if
    else if (array_eq(output_var%data_shape(:)%size(), [mp, msn])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow
      else
        decomp => output_decomp%land_snow
      end if
    else if (array_eq(output_var%data_shape(:)%size(), [mp, nrb])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad
      else
        decomp => output_decomp%land_rad
      end if
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon
      else
        decomp => output_decomp%land_plantcarbon
      end if
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncs])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon
      else
        decomp => output_decomp%land_soilcarbon
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine cable_output_decomp_associate

  subroutine cable_output_decomp_associate_restart(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    type(cable_output_decomp_t), pointer :: output_decomp

    select case (output_var%aggregator%type())
    case ("int32")
      output_decomp => decomps_grid_restart_int32
    case ("real32")
      output_decomp => decomps_grid_restart_real32
    case ("real64")
#ifdef ENFORCE_SINGLE_PRECISION
      output_decomp => decomps_grid_restart_real32
#else
      output_decomp => decomps_grid_restart_real64
#endif
    case default
      call cable_abort("Unexpected data type for output variable " // output_var%field_name, __FILE__, __LINE__)
    end select

    if (array_eq(output_var%data_shape(:)%size(), [mp])) then
      decomp => output_decomp%patch
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ms])) then
      decomp => output_decomp%patch_soil
    else if (array_eq(output_var%data_shape(:)%size(), [mp, msn])) then
      decomp => output_decomp%patch_snow
    else if (array_eq(output_var%data_shape(:)%size(), [mp, nrb])) then
      decomp => output_decomp%patch_rad
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncp])) then
      decomp => output_decomp%patch_plantcarbon
    else if (array_eq(output_var%data_shape(:)%size(), [mp, ncs])) then
      decomp => output_decomp%patch_soilcarbon
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine cable_output_decomp_associate_restart

end submodule cable_output_decomp_smod
