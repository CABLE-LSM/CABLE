!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
! 
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located 
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Reorders JULES and Unified Model tiled parameters to an order
!          appropriate for CABLE and ACCESS-CM2
!           
!          assumes input from JULES/UM is on 17 tiles with order:
!          (1) dec broadlf, (2) evgn broadlf tropics, (3) evgn broadlf temp,
!          (4) dec needlelf, (5) evgn needlelf, (6) grass c3, (7)crop c3, 
!          (8) pasture c3, (9) grass c4, (10) crop c4, (11) pasture c4,
!          (12) dec shrub, (13) evgn shrub, (14) urban, (15) lakes, 
!          (16) baresoil, (17) permanent ice
!
!          mapping at time of writing is:
!
!  CABLE    egnneedle tile 	1 -> 	ndl_leaf_eg 	tile 5 JULES
!	        egnbroad	    2 ->	brd_leaf_evg_	0.5*(2+3) [special case]
!	        decneedle	    3 -> 	ndl_leaf_dec	4
!	        decbroad 	    4 ->	brd_leaf_dec	1
!	        shrub 		    5 ->	shrub_dec/evg	12 or 13
!	        grass c3	    6 ->	c3_grass	    6
!	        grass c4	    7 ->	c4_grass 	    9
!	        tundra		    8 -> 	c3_grass 6 or evg shrub 13
!	        crop c3		    9 -> 	c3_crop	        7
!	        crop_c4	        10 ->	c4_crop	        10
!	        wetland	        11 -> 	lake 15 or c3 grass  6 or c4 grass 9
!	        empty		    12 ->	c3_grass 	    6
!	        empty		    13 -> 	c3_grass	    6
!	        bare soil	    14 ->	soil		    16
!	        urban		    15 -> 	urban		    14
!	        lakes		    16 ->	lakes		    15
!	        ice		        17 ->	ice		        17
!
!  CABLE shrub map largely coincides with JULES tile 12 map
!  CABLE tundra map largely conincides with JULES tile 13 map
!
!  integer switch shrub permits CABLE shrub to take another value
!
!  integer switch tundra permits CABLE tundra to take another value
!
!  integer switch wetland13/17 permits CABLE wetlands to take another value e.g. 
!  lake (15), c3 grass (6) or c4 grass (9) (ideally it depends on the process)
!  note cannot take a lake value if only veg pfts are available.
!
!  Future development needed to permit external/namelist prescription of the 
!  integer switches for shrub, tundra, wetland13/17  
!
!  Called from: ukca routines surfddr, ddepaer & ddepaer_inc_sedi
!
! ===============================================================================

MODULE cable_jules_links_mod
    
    IMPLICIT NONE
    
    PUBLIC tile_resort_1D, tile_resort_2D
    PRIVATE tile_order
    
CONTAINS
    SUBROUTINE tile_order(ntiles,t_sort)
        !sets new order of tiles
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: ntiles      !number of tiles - should be 4, 13 or 17
        INTEGER, DIMENSION(ntiles), INTENT(OUT) :: t_sort
        
        INTEGER :: shrub                   !switch to specify shrub tile values
        INTEGER :: tundra                  !switch to specify tundra tile values
        INTEGER :: wetland13               !switch to specify wetland tile values
        INTEGER :: wetland17               !switch to specify wetland tile values
        
        !edit shrub, tundra and wetland mapping - other tile mapping is hard-wired
        shrub = 12                         !to be 12 (default) or 13
        tundra = 13                        !to be 6, 12 or 13 (default)
        wetland13 = 6                      !to be 6 (default for 13 tiles) or 9 
        wetland17 = 15                     !to be 6, 9 or 15 (default for 13 tiles) 
 
        !set the tile order -----------------------------------------------------
        IF (ntiles == 17) THEN             !all tiles          
            t_sort = (/5,2,4,1,shrub,6,9,tundra,7,10,wetland17,6,6,16,14,15,17/)
        
        ELSEIF (ntiles == 13) THEN         !veg tiles only
            t_sort = (/5,2,4,1,shrub,6,9,tundra,7,10,wetland13,6,6/)
        
        ELSEIF (ntiles == 4) THEN          !non-veg tiles only
            t_sort = (/3,1,2,4/)
       
        ELSE
            !error message
            t_sort(:) = MAX(1,MIN(ntiles-1,6)) !catch all as c3 grass, bare soil
        ENDIF
               
        RETURN
    END SUBROUTINE

    SUBROUTINE tile_resort_2D(ntiles, nvar, t_params)
        !2D array case
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: ntiles      !number of tiles - should be 13 or 17
        INTEGER, INTENT(IN) :: nvar        !number of variables in array to sort
        
        REAL, DIMENSION(ntiles,nvar), INTENT (INOUT) :: t_params
                                           !variable array to resort
        
        !working variables
        INTEGER :: i                       !looping variable
        INTEGER, DIMENSION(:), ALLOCATABLE :: t_sort
                                           !ordering of tiles
        REAL, DIMENSION(:,:), ALLOCATABLE :: t_params_new
                                           !resorted variable array
        !END header
        ALLOCATE(t_sort(ntiles))
        ALLOCATE(t_params_new(ntiles,nvar))
 
        !set the tile order -----------------------------------------------------
        CALL tile_order(ntiles, t_sort)
        
        !reorder ----------------------------------------------------------------
        
        DO i = 1,ntiles
            t_params_new(i,:) = t_params(t_sort(i),:)
            !special cases
            IF ((ntiles==13) .or. (ntiles==17)) THEN
                t_params_new(2,:) = 0.5*(t_params(2,:)+t_params(3,:))
            END IF
        END DO
        t_params(:,:) = t_params_new(:,:)
        
        DEALLOCATE(t_sort,t_params_new)
        
        RETURN
        
    END SUBROUTINE tile_resort_2D
    
    SUBROUTINE tile_resort_1D(ntiles, t_params)
        !1D array case
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: ntiles      !number of tiles - should be 13 or 17
        
        REAL, DIMENSION(ntiles), INTENT (INOUT) :: t_params
                                           !variable array to resort
        
        !working variables
        INTEGER :: i                       !looping variable
        INTEGER, DIMENSION(:), ALLOCATABLE :: t_sort
                                           !ordering of tiles
        REAL, DIMENSION(:), ALLOCATABLE :: t_params_new
                                           !resorted variable array
        !END header
        ALLOCATE(t_sort(ntiles))
        ALLOCATE(t_params_new(ntiles))
 
        !set the tile order -----------------------------------------------------
        CALL tile_order(ntiles, t_sort)
        
        !reordering -------------------------------------------------------------
        DO i = 1,ntiles
            t_params_new(i) = t_params(t_sort(i))
            !special cases
            IF ((ntiles==13) .or. (ntiles==17)) THEN
                t_params_new(2) = 0.5*(t_params(2)+t_params(3))
            END IF
        END DO
        t_params(:) = t_params_new(:)
        
        DEALLOCATE(t_sort,t_params_new)
        
        RETURN
        
    END SUBROUTINE tile_resort_1D
END MODULE

    