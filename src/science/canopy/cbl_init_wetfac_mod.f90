MODULE cable_init_wetfac_mod
   
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE initialize_wetfac( &
        mp, &
        ssnow_wetfac, &
        soil_swilt, &
        soil_sfc, & 
        ssnow_wb, &
        ssnow_wbice, &
        ssnow_snowd, &
        veg_iveg, &
        met_tk, &
        Ctfrz &
    )

        !! Initialize the surface wetness factor of the soil/snow (ssnow_wetfac) array

        ! Imports
        ! use cable_surface_types_mod, ONLY: lakes => lakes_cable        
        USE grid_constants_mod_cbl, ONLY: lakes_cable
        USE cable_def_types_mod,  ONLY : r_2
        USE cable_other_constants_mod, ONLY : wilt_limitfactor
                                
        IMPLICIT NONE                     

        ! Arguments (See soil params file for units?)
        INTEGER, INTENT(IN) :: mp !! Number of active land points
        REAL, INTENT(INOUT) :: ssnow_wetfac(mp) !! Soil/snow wetness factor
        REAL, INTENT(IN) :: soil_swilt(mp) !! Wilting factor, point at which plants in soil start to wilt (unit?)
        REAL, INTENT(IN) :: soil_sfc(mp) !! Soil sfc array
        REAL(r_2), INTENT(IN) :: ssnow_wb(mp) !! Soil/snow moisture
        REAL(r_2), INTENT(IN) :: ssnow_wbice(mp) !! Frozen component of soil/snow moisture
        REAL, INTENT(IN) :: ssnow_snowd(mp) !! Soil/snow snow depth (mm of water)
        REAL, INTENT(IN) :: met_tk(mp) !! Air temperature (kelvin)
        REAL, INTENT(IN) :: Ctfrz  !! Freezing temperature (CONSTANT)
        integer, INTENT(IN) :: veg_iveg(mp) !! Surface types (one of the 17, integer)

        ! Local variables
        REAL :: wilting_pt(mp) ! Wilting point
        REAL :: wetfac_num(mp) ! Wetness factor numerator
        REAL :: wetfac_den(mp) ! Wetness factor denominator
        REAL :: ice_ratio ! Ice ratio
        REAL :: ice_factor ! Ice factor
        INTEGER :: i ! Index to iterate through

        ! Work out the numerator / denominator
        wilting_pt(:) = soil_swilt(:) / wilt_limitfactor 
        wetfac_num(:) = REAL(ssnow_wb(:)) - wilting_pt(:)
        wetfac_den(:) = REAL(soil_sfc(:) - wilting_pt(:))
        wetfac_den(:) = MAX(0.0830, wetfac_den(:)) !? What is this magical number? It should be a constant and described

        ! Set some "meaningful" defaults
        ssnow_wetfac(:) = wetfac_num(:) / wetfac_den(:)
        ssnow_wetfac(:) = MIN(1.0, ssnow_wetfac(:))
        ssnow_wetfac(:) = MAX(0.0, ssnow_wetfac(:))

        ! Loop through the number of land points
        DO i=1,mp

            ! Ultimately reduces surface wetness considering wetness locked up in ice
            IF (ssnow_wbice(i) > 0.0) THEN
                
                ! ice_ratio = ice moisture / total moisture (** ?)
                ice_ratio  = (ssnow_wbice(i) / ssnow_wb(i))**2 ! Why is this squared?

                !~ 1-Ice_ratio^2 
                ice_factor = 1._r_2 - MIN(0.2_r_2, ice_ratio)
                ice_factor = REAL(MAX(0.5_r_2, ice_factor))
                ssnow_wetfac(i) = ssnow_wetfac(i) * ice_factor

            END IF

            ! If snow depth is greater than 0.1m then 90% wet
            IF (ssnow_snowd(i) > 0.1) THEN
                ssnow_wetfac(i) = 0.9
            END IF

            ! If we are on a lake
            IF (veg_iveg(i) == lakes_cable) THEN

                ! When the air temperature is >= +5 deg above freezing it is 100% wet
                IF ( met_tk(i) >= Ctfrz + 5. ) THEN
                    ssnow_wetfac(i) = 1.0
                END IF

                ! When the air temperature is < +5 deg above freezing it is 70% wet
                IF( met_tk(i) < Ctfrz + 5. ) THEN
                    ssnow_wetfac(i) = 0.7
                END IF

            END IF

        ENDDO

        RETURN
    
    END SUBROUTINE initialize_wetfac  

END MODULE cable_init_wetfac_mod