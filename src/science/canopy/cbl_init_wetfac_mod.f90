module cable_init_wetfac_mod
   
    implicit none

    contains

    subroutine initialize_wetfac( &
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
        use grid_constants_mod_cbl, ONLY: lakes_cable
        use cable_def_types_mod,  ONLY : r_2
        use cable_other_constants_mod, ONLY : wilt_limitfactor
                                
        implicit none                     

        ! Arguments (See soil params file for units?)
        integer, intent(in) :: mp !! Number of active land points
        real, intent(inout) :: ssnow_wetfac(mp) !! Soil/snow wetness factor
        real, intent(in) :: soil_swilt(mp) !! Wilting factor, point at which plants in soil start to wilt (unit?)
        real, intent(in) :: soil_sfc(mp) !! Soil sfc array
        real(r_2), intent(in) :: ssnow_wb(mp) !! Soil/snow moisture
        real(r_2), intent(in) :: ssnow_wbice(mp) !! Frozen component of soil/snow moisture
        real, intent(in) :: ssnow_snowd(mp) !! Soil/snow snow depth (mm of water)
        real, intent(in) :: met_tk(mp) !! Air temperature (kelvin)
        real, intent(in) :: Ctfrz  !! Freezing temperature (CONSTANT)
        integer, intent(in) :: veg_iveg(mp) !! Surface types (one of the 17, integer)

        ! Local variables
        real :: wilting_pt(mp) ! Wilting point
        real :: wetfac_num(mp) ! Wetness factor numerator
        real :: wetfac_den(mp) ! Wetness factor denominator
        real :: ice_ratio ! Ice ratio
        real :: ice_factor ! Ice factor
        integer :: i ! Index to iterate through

        ! Work out the numerator / denominator
        wilting_pt(:) = soil_swilt(:) / wilt_limitfactor 
        wetfac_num(:) = real(ssnow_wb(:)) - wilting_pt(:)
        wetfac_den(:) = real(soil_sfc(:) - wilting_pt(:))
        wetfac_den(:) = max(0.0830, wetfac_den(:)) !? What is this magical number? It should be a constant and described

        ! Set some "meaningful" defaults
        ssnow_wetfac(:) = wetfac_num(:) / wetfac_den(:)
        ssnow_wetfac(:) = min(1.0, ssnow_wetfac(:))
        ssnow_wetfac(:) = max(0.0, ssnow_wetfac(:))

        ! Loop through the number of land points
        do i=1,mp

            ! Ultimately reduces surface wetness considering wetness locked up in ice
            if( ssnow_wbice(i) > 0.0 ) then
                
                ! ice_ratio = ice moisture / total moisture (** ?)
                ice_ratio  = (ssnow_wbice(i) / ssnow_wb(i))**2 ! Why is this squared?

                !~ 1-Ice_ratio^2 
                ice_factor = 1._r_2 - min(0.2_r_2, ice_ratio)
                ice_factor = real(max(0.5_r_2, ice_factor))
                ssnow_wetfac(i) = ssnow_wetfac(i) * ice_factor

            end if

            ! This looks like something regarding snow depth?
            if( ssnow_snowd(i) > 0.1) then 
                ssnow_wetfac(i) = 0.9
            end if

            ! If we are on a lake
            if ( veg_iveg(i) == lakes_cable ) then

                ! When the air temperature is >= +5 deg above freezing it is 100% wet
                if ( met_tk(i) >= Ctfrz + 5. ) then
                    ssnow_wetfac(i) = 1.0
                end if

                ! When the air temperature is < +5 deg above freezing it is 70% wet
                if( met_tk(i) < Ctfrz + 5. ) then
                    ssnow_wetfac(i) = 0.7
                end if

            end if

        enddo

        return
    
    end subroutine initialize_wetfac  

end module cable_init_wetfac_mod