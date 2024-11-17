MODULE cable_init_wetfac_mod
    !! Module containing subroutine to initialise the surface wetness factor
    !! of the soil/snow (ssnow_wetfac) array
   
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
        !! ## Purpose
        !!
        !! Initialize the surface wetness factor of the soil/snow 
        !! (ssnow_wetfac) array
        !!
        !! ## Method
        !!
        !! **Warning**: The original subroutine from which this was ported
        !! lacks any documented methodolody.
        !!
        !! ## References
        !!
        !! **Warning**: The original subroutine from which this was ported
        !! lacks any literature reference.

        

        ! Imports
        USE grid_constants_mod_cbl, ONLY: lakes_cable
        USE cable_def_types_mod,  ONLY : r_2
        USE cable_other_constants_mod, ONLY : wilt_limitfactor
                                
        IMPLICIT NONE                     

        ! Arguments (See soil params file for units?)
        INTEGER, INTENT(IN) :: mp !! Number of active land points
        REAL, INTENT(INOUT) :: ssnow_wetfac(mp) !! Surface wetness factor at current time step
        REAL, INTENT(IN) :: soil_swilt(mp) !! Wilting factor, point at which plants in soil start to wilt
        REAL, INTENT(IN) :: soil_sfc(mp) !! Volumetric H20 @ field capacity
        REAL(r_2), INTENT(IN) :: ssnow_wb(mp) !! Volumetric soil moisture (solid+liquid)
        REAL(r_2), INTENT(IN) :: ssnow_wbice(mp) !! Soil ice
        REAL, INTENT(IN) :: ssnow_snowd(mp) !! Soil/snow snow depth (mm of water)
        REAL, INTENT(IN) :: met_tk(mp) !! Air temperature (Kelvin)
        REAL, INTENT(IN) :: Ctfrz  !! Freezing temperature (Kelvin)
        INTEGER, INTENT(IN) :: veg_iveg(mp) !! Surface types

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
        wetfac_den(:) = MAX(0.0830, wetfac_den(:)) ! WARNING: CABLE#457

        ! Set some "meaningful" defaults
        ssnow_wetfac(:) = wetfac_num(:) / wetfac_den(:)
        ssnow_wetfac(:) = MIN(1.0, ssnow_wetfac(:))
        ssnow_wetfac(:) = MAX(0.0, ssnow_wetfac(:))

        ! Loop through the number of land points
        DO i=1,mp

            ! Ultimately reduces surface wetness considering wetness locked up in ice
            IF (ssnow_wbice(i) > 0.0) THEN
                
                ! ice_ratio = ice moisture / total moisture (** ?)
                ice_ratio  = (ssnow_wbice(i) / ssnow_wb(i))**2

                !~ 1-Ice_ratio^2 
                ice_factor = 1._r_2 - MIN(0.2_r_2, ice_ratio)
                ice_factor = REAL(MAX(0.5_r_2, ice_factor))
                ssnow_wetfac(i) = ssnow_wetfac(i) * ice_factor

            END IF

            ! If snow depth is greater than 0.1m then soil is at 90% of the total available water
            IF (ssnow_snowd(i) > 0.1) THEN
                ssnow_wetfac(i) = 0.9
            END IF

            ! If we are on a lake
            IF (veg_iveg(i) == lakes_cable) THEN

                ! When the air temperature is >= +5 deg above freezing it at 100% of the total available water
                IF ( met_tk(i) >= Ctfrz + 5. ) THEN
                    ssnow_wetfac(i) = 1.0
                END IF

                ! When the air temperature is < +5 deg above freezing it is at 70% of the total available water
                IF( met_tk(i) < Ctfrz + 5. ) THEN
                    ssnow_wetfac(i) = 0.7
                END IF

            END IF

        ENDDO

        RETURN
    
    END SUBROUTINE initialize_wetfac  

END MODULE cable_init_wetfac_mod