SUBROUTINE crop_driver(ktau,ktauday,doy,climate,ssnow,soil,casaflux,crop)

  use crop_def,            only: crop_type, nc
  use crop_module
  use casavariable,        only: casa_flux
  use cable_def_types_mod, only: climate_type, soil_snow_type, soil_parameter_type, &
                                 dp => r_2

  implicit none
  
  integer,                   intent(in)    :: ktau
  integer,                   intent(in)    :: ktauday
  integer,                   intent(in)    :: doy
  type(climate_type),        intent(in)    :: climate
  type(soil_snow_type),      intent(in)    :: ssnow
  type(soil_parameter_type), intent(in)    :: soil
  type(casa_flux),           intent(inout) :: casaflux
  type(crop_type),           intent(inout) :: crop
  
  ! local
  integer :: ic,sl  ! loop counters: crop type, soil layer
  real(dp), dimension(nc) :: fPHU_day

  
  do ic=1, nc   
write(70,*) 'DOY: ', doy    
    select case (crop%state(ic))  ! maybe if-else construct makes more sense
      case (0)  ! fallow
write(70,*) 'crop%sowing_doymin: ', crop%sowing_doymin(1)
        if (doy > crop%sowing_doymin(ic) .AND. doy < crop%sowing_doymax(ic)) then 
          if (mod(ktau,ktauday) == ktauday/2) then  ! midday   
             call sowing(doy,soil,crop)
write(70,*) 'crop%state: ', crop%state
write(70,*) 'crop%Tbase: ', crop%Tbase
          end if
        end if

      case (1)  ! sown
        if (mod(ktau,ktauday) == 0) then  ! end of day 
          call germination(doy,climate,ssnow,soil,crop)
        end if

      case (2)  ! emerging
        call emergence(doy,crop)

      case (3)  ! growing
        if (mod(ktau,ktauday) == 0) then ! end of day
          ! update phenological heat units
          fPHU_day(ic)  = heat_units(crop%Tbase(ic),climate%dtemp(ic)) / crop%PHU_maturity(ic)
          crop%fPHU(ic) = crop%fPHU(ic) + fPHU_day(ic)

write(70,*) 'climate%dtemp: ', climate%dtemp
write(70,*) 'crop%fPHU: ', crop%fPHU

          ! calculate GPP (do this first, calculated every timestep, not every day)


          ! calculate maintenance and growth respiration
          !call casa_rplant(...,...)

          ! calculate carbon allocation
          call C_allocation_crops(casaflux,crop)

          ! update root length and vegetation height
          
          
          ! harvest if enough PHU accumulated
          if (crop%fPHU(ic) >= 1.0_dp) then
            call harvest(doy,crop)
          end if
        end if
          
     end select
   end do
 
END SUBROUTINE crop_driver

