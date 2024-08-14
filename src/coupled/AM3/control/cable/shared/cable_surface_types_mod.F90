! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cable_surface_types_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains CABLE surface type information and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE max_dimensions,   ONLY:                                                    &
    elev_tile_max,                                                             &
    ntype_max

USE missing_data_mod, ONLY: imdi

!-----------------------------------------------------------------------------
! Module variables.
!-----------------------------------------------------------------------------
USE jules_surface_types_mod, ONLY:                                             &
  nnpft,                                                                       &
  ncpft,                                                                   &
  npft,                                                                &
  nnvg,                                                                &
  ntype,                                                                       &
  urban,                                                                &
  lake,                                                                &
  soil,                                                                &
  ice 

IMPLICIT NONE

! Index of the various surface types used by CABLE"
INTEGER ::                                                                     &
   evergreen_needleleaf = imdi,                                                &
   evergreen_broadleaf  = imdi,                                                &
   deciduous_needleleaf = imdi,                                                &
   deciduous_broadleaf  = imdi,                                                &
   shrub_cable          = imdi,                                                &
   c3_grassland         = imdi,                                                &
   c4_grassland         = imdi,                                                &
   tundra               = imdi,                                                &
   c3_cropland          = imdi,                                                &
   c4_cropland          = imdi,                                                &
   wetland              = imdi,                                                &
   empty1               = imdi,                                                &
   empty2               = imdi,                                                &
   barren_cable         = imdi,                                                &
   urban_cable          = imdi,                                                & 
   lakes_cable          = imdi,                                                &
   ice_cable            = imdi

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / cable_surface_types/                                               &
  npft, nnvg, & 
  evergreen_needleleaf, evergreen_broadleaf, deciduous_needleleaf,          &
  deciduous_broadleaf, shrub_cable, c3_grassland, c4_grassland,             &
  tundra, c3_cropland, c4_cropland, wetland, empty1, empty2,                &
  barren_cable, urban_cable, lakes_cable, ice_cable

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CABLE_SURFACE_TYPES_MOD'

CONTAINS

SUBROUTINE check_surface_type_value ( surface_type, surface_type_name,         &
   min_value, max_value, RoutineName, errorstatus, nchecks )

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

INTEGER          :: surface_type, min_value, max_value, errorstatus, nchecks
CHARACTER(LEN=*) :: surface_type_name, RoutineName

IF ( surface_type > 0 ) THEN
  nchecks = nchecks + 1
  IF ( surface_type < min_value .OR. surface_type > max_value ) THEN
    errorstatus = 101
    CALL jules_print(RoutineName,                                              &
       TRIM(surface_type_name) // " tile is given but is out of range")
  END IF
END IF

END SUBROUTINE  check_surface_type_value

SUBROUTINE print_nlist_cable_surface_types()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

INTEGER :: i, n ! Loop counter

CHARACTER(LEN=50000) :: lineBuffer


!-----------------------------------------------------------------------------
! This needs to be implemented corresponding to cable_surface_types
CALL jules_print('cable_surface_types',                                        &
                 'Contents of namelist cable_surface_types')

!WRITE(lineBuffer, *) '  npft = ', npft
!CALL jules_print('cable_surface_types', lineBuffer)

!IF ( brd_leaf > 0 ) THEN
!  WRITE(lineBuffer, *) '  brd_leaf = ', brd_leaf
!  CALL jules_print('cable_surface_types', lineBuffer)
!END IF

!WRITE(lineBuffer, *) '   = ', 
!CALL jules_print('cable_surface_types', lineBuffer)

!IF ( > 0 ) THEN
!  WRITE(lineBuffer, *) '   = ', 
!  CALL jules_print('cable_surface_types', lineBuffer)
!END IF

CALL jules_print('cable_surface_types',                                        &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_cable_surface_types



SUBROUTINE set_derived_variables_cable_surface_types()

USE jules_print_mgr, ONLY: jules_print, jules_message

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Derive ntype and nnpft from the namelist values
!-----------------------------------------------------------------------------
ntype = npft + nnvg
nnpft = npft - ncpft

CALL jules_print('set_derived_variables_cable_surface_types',                  &
                 'Derived variables from cable_surface_types')

!WRITE(jules_message, *) '  ntype = ', ntype
!CALL jules_print('set_derived_variables_cable_surface_types', jules_message)
!
!WRITE(jules_message, *) '  nnpft = ', nnpft
!CALL jules_print('set_derived_variables_cable_surface_types', jules_message)

CALL jules_print('set_derived_variables_cable_surface_types',                  &
    '- - - - - - end of derived variables - - - - - -')

RETURN

END SUBROUTINE set_derived_variables_cable_surface_types



SUBROUTINE check_cable_surface_types()

USE max_dimensions, ONLY: npft_max, ncpft_max, nnvg_max

USE ereport_mod, ONLY: ereport
USE jules_print_mgr, ONLY: jules_print, jules_message

!-----------------------------------------------------------------------------
! Description:
!   Checks cable_surface_types namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: i, nchecks ! Loop counter

INTEGER :: errorstatus

CHARACTER(LEN=*), PARAMETER ::                                                 &
   RoutineName = 'CHECK_CABLE_SURFACE_TYPES'
!-----------------------------------------------------------------------------
! Check that the given values are less than the fixed values for IO
!-----------------------------------------------------------------------------
IF ( npft > npft_max ) THEN
  errorstatus = 101
  CALL ereport(RoutineName, errorstatus,                                       &
               "npft > npft_max - increase npft_max and recompile")
END IF

IF ( nnvg > nnvg_max ) THEN
  errorstatus = 101
  CALL ereport(RoutineName, errorstatus,                                       &
               "nnvg > nnvg_max - increase nnvg_max and recompile")
END IF

!-----------------------------------------------------------------------------
! Check values for the specific surface types are sensible
!-----------------------------------------------------------------------------
! PFTs
errorstatus = 0
nchecks = 0

!!CALL check_surface_type_value ( brd_leaf, "brd_leaf", 1, npft,                 &
!!   RoutineName, errorstatus, nchecks )
!!CALL check_surface_type_value ( , "", 1, npft,         &
!!   RoutineName, errorstatus, nchecks )
!
!! PFT surface types must come before non-veg types, so if urban, lake, soil,
!! ice, urban_canyon or urban_roof are given (i.e. > 0) then they must be > npft
!! A soil type is required
!!CALL check_surface_type_value ( urban, "urban", npft+1, ntype,                 &
!!   RoutineName, errorstatus, nchecks )
!!CALL check_surface_type_value ( , "", npft+1, ntype,                   &
!!   RoutineName, errorstatus, nchecks )

!jhan:this will need to be properly implemented
nchecks = 17
! Check that all present surface types have been checked for range compliance
! This check should also ensure that a check is added for each new surface type
IF ( nchecks /= ntype ) THEN
  errorstatus = 101
  CALL jules_print(RoutineName,                                                &
     "At least one surface type in namelist does not have a range check.")
  WRITE(jules_message,'(A,I0,A,I0)')                                           &
     "These should be the same; ntype = ", ntype, ", nchecks = ", nchecks
  CALL jules_print(RoutineName, jules_message)
END IF

! Now that all surface types have been checked issue abort if required
IF ( errorstatus > 0 ) THEN
  CALL ereport(RoutineName, errorstatus,                                       &
     "Error(s) found. Please see job.out for information ")
END IF

END SUBROUTINE check_cable_surface_types

SUBROUTINE read_nml_cable_surface_types (unitnumber)

! Description:
!  Read the CABLE_SURFACE_TYPES namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_CABLE_SURFACE_TYPES'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 19 !+ ntype_max

TYPE my_namelist
  !!SEQUENCE
  INTEGER :: npft
  INTEGER :: nnvg
  INTEGER :: evergreen_needleleaf
  INTEGER :: evergreen_broadleaf
  INTEGER :: deciduous_needleleaf
  INTEGER :: deciduous_broadleaf
  INTEGER :: shrub_cable
  INTEGER :: c3_grassland
  INTEGER :: c4_grassland
  INTEGER :: tundra
  INTEGER :: c3_cropland
  INTEGER :: c4_cropland
  INTEGER :: wetland
  INTEGER :: empty1
  INTEGER :: empty2
  INTEGER :: barren_cable
  INTEGER :: urban_cable
  INTEGER :: lakes_cable
  INTEGER :: ice_cable
  !INTEGER :: tile_map_ids(ntype_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = cable_surface_types, IOSTAT = errorstatus,    &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist cable_surface_types",iomessage)

  my_nml % npft                 = npft
  my_nml % nnvg                 = nnvg
  my_nml % evergreen_needleleaf = evergreen_needleleaf
  my_nml % evergreen_broadleaf  = evergreen_broadleaf
  my_nml % deciduous_needleleaf = deciduous_needleleaf
  my_nml % deciduous_broadleaf  = deciduous_broadleaf
  my_nml % shrub_cable          = shrub_cable
  my_nml % c3_grassland         = c3_grassland
  my_nml % c4_grassland         = c4_grassland
  my_nml % tundra               = tundra
  my_nml % c3_cropland          = c3_cropland
  my_nml % c4_cropland          = c4_cropland
  my_nml % wetland              = wetland
  my_nml % empty1               = empty1
  my_nml % empty2               = empty2
  my_nml % barren_cable         = barren_cable
  my_nml % urban_cable          = urban_cable
  my_nml % lakes_cable          = lakes_cable
  my_nml % ice_cable            = ice_cable

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  npft                  = my_nml % npft
  nnvg                  = my_nml % nnvg
  evergreen_needleleaf  = my_nml % evergreen_needleleaf
  evergreen_broadleaf   = my_nml % evergreen_broadleaf
  deciduous_needleleaf  = my_nml % deciduous_needleleaf
  deciduous_broadleaf   = my_nml % deciduous_broadleaf
  shrub_cable           = my_nml % shrub_cable
  c3_grassland          = my_nml % c3_grassland
  c4_grassland          = my_nml % c4_grassland
  tundra                = my_nml % tundra
  c3_cropland           = my_nml % c3_cropland
  c4_cropland           = my_nml % c4_cropland
  wetland               = my_nml % wetland
  empty1                = my_nml % empty1
  empty2                = my_nml % empty2
  barren_cable          = my_nml %barren_cable
  urban_cable           = my_nml % urban_cable
  lakes_cable           = my_nml % lakes_cable
  ice_cable             = my_nml % ice_cable

END IF

CALL mpl_type_free(mpl_nml_type,icode)

soil         = barren_cable
ice          = ice_cable 
lake         = lakes_cable 
urban        = urban_cable 
 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_cable_surface_types

END MODULE cable_surface_types_mod
