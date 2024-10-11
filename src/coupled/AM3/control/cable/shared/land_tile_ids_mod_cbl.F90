! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting Tile ID numbers, which are used to identify which land
!   surface tiles are present.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE land_tile_ids_mod_cbl

USE max_dimensions,   ONLY:                                                    &
    ntype_max,                                                                 &
    snow_layers_max,                                                           &
    elev_tile_max

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

!-----------------------------------------------------------------------

INTEGER, PRIVATE :: i ! Loop counter

INTEGER :: surface_type_ids(ntype_max) = imdi
                ! Array which maps pseudo levels to tile types

INTEGER :: ml_snow_type_ids(ntype_max * snow_layers_max) = imdi
                ! Array which maps pseudo levels to tile types

INTEGER :: tile_ids_in(ntype_max * snow_layers_max) = imdi
                ! Tile IDs in the input header

CONTAINS

SUBROUTINE set_surface_type_ids_cbl( )
USE ereport_mod,             ONLY: ereport

USE jules_print_mgr, ONLY:                                                     &
   jules_message,                                                              &
   jules_print,                                                                &
   jules_format

USE errormessagelength_mod, ONLY: errormessagelength

USE cable_surface_types_mod, ONLY:                                             & 
   ntype,                                                                      &
   evergreen_needleleaf,                                                       &  
   evergreen_broadleaf,                                                        &  
   deciduous_needleleaf,                                                       &  
   deciduous_broadleaf,                                                        &  
   shrub_cable,                                                                &  
   c3_grassland,                                                               &  
   c4_grassland,                                                               &  
   tundra,                                                                     &  
   c3_cropland,                                                                &  
   c4_cropland,                                                                &  
   wetland,                                                                    &  
   empty1,                                                                     &  
   empty2,                                                                     &  
   barren_cable,                                                               &  
   urban_cable,                                                                &  
   lakes_cable,                                                                &  
   ice_cable  

IMPLICIT NONE

INTEGER      :: i ! Loop counter
INTEGER      :: errorstatus
CHARACTER(LEN=18), PARAMETER      :: routinename='set_tile_id_arrays'

!There is presently no other option for CABLE surface type indexing
surface_type_ids( evergreen_needleleaf ) = 1 
surface_type_ids( evergreen_broadleaf  ) = 2
surface_type_ids( deciduous_needleleaf ) = 3
surface_type_ids( deciduous_broadleaf  ) = 4
surface_type_ids( shrub_cable          ) = 5  
surface_type_ids( c3_grassland         ) = 6
surface_type_ids( c4_grassland         ) = 7
surface_type_ids( tundra               ) = 8
surface_type_ids( c3_cropland          ) = 9
surface_type_ids( c4_cropland          ) = 10
surface_type_ids( wetland              ) = 11
surface_type_ids( empty1               ) = 12
surface_type_ids( empty2               ) = 13
surface_type_ids( barren_cable         ) = 14
surface_type_ids( urban_cable          ) = 15
surface_type_ids( lakes_cable          ) = 16
surface_type_ids( ice_cable            ) = 17

! Print the surface types that are present
WRITE(jules_format,'(a3,i3,a8)') '(a,',ntype,'(1x,i6))'
WRITE(jules_message,jules_format) ' Surface types present =',                  &
   surface_type_ids(1:ntype)
CALL jules_print(routinename, jules_message)

! Check that all surface types have been specified
IF ( ANY( surface_type_ids(1:ntype) == imdi ) ) THEN
  errorstatus = 30
  WRITE(jules_message, '(A,I6)')                                               &
     ' All surface types need to be specified. Please see job.out.'
  CALL ereport ( routinename, errorstatus, jules_message )
END IF

! Check that tile IDs are unique
DO i = 1, ntype
  IF ( COUNT( surface_type_ids(:) == surface_type_ids(i) ) /= 1 ) THEN
    errorstatus = 30
    WRITE(jules_message, '(A,I6)')                                             &
       ' Surface type ID not unique :', surface_type_ids(i)
    CALL ereport ( routinename, errorstatus, jules_message )
  END IF
END DO

END SUBROUTINE set_surface_type_ids_cbl

END MODULE land_tile_ids_mod_cbl
