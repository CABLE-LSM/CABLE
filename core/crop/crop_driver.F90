SUBROUTINE crop_driver(ktau,ktauday,doy,climate,ssnow,soil,veg,casaflux,casamet,&
                       casapool,crop)

  use crop_def,            only: crop_type, nc, baresoil, sown, emergent, growing, &
                                 DMtoC
  use crop_module
  use casavariable,        only: casa_flux, casa_met, casa_pool
  use cable_def_types_mod, only: climate_type, soil_snow_type, soil_parameter_type, &
                                 veg_parameter_type, dp => r_2

  implicit none
  
  integer,                   intent(in)    :: ktau
  integer,                   intent(in)    :: ktauday
  integer,                   intent(in)    :: doy
  type(climate_type),        intent(in)    :: climate
  type(soil_snow_type),      intent(in)    :: ssnow
  type(soil_parameter_type), intent(in)    :: soil
  type(veg_parameter_type),  intent(inout) :: veg
  type(casa_flux),           intent(inout) :: casaflux
  type(casa_met),            intent(inout) :: casamet
  type(casa_pool),           intent(inout) :: casapool
  type(crop_type),           intent(inout) :: crop
  
  ! local
  integer :: ic,sl  ! loop counters: crop type, soil layer
  real(dp), dimension(nc) :: fPHU_day
  real(dp), dimension(nc) :: SLA_C
  real(dp), dimension(nc) :: VU  ! vernalisation units
  real(dp), dimension(nc) :: fVU ! vernalisation requirements
  logical,  dimension(nc) :: vt  ! vernalisation been accounted for?
real(dp),dimension(nc) :: NPP_test
  
  
  do ic=1, nc   
write(70,*) 'DOY: ', doy
write(70,*) 'VU(i):', VU(ic)
    if (crop%state(ic) == baresoil) then
      if (doy > crop%sowing_doymin(ic) .AND. doy < crop%sowing_doymax(ic)) then 
         call sowing(doy,soil,crop)
write(70,*) 'crop%state: ', crop%state
write(70,*) 'crop%Tbase: ', crop%Tbase
      end if

    else if (crop%state(ic) == sown) then
      call germination(doy,climate,ssnow,soil,crop)
casapool%Cplant(ic,:) = 0.0_dp ! shouldn't be needed here!! Check initialisation
casamet%glai(ic) = 0.0_dp

    else if (crop%state(ic) == emergent .or. crop%state(ic) == growing) then
       ! update phenological heat units (start at 0 at germination!)
       fPHU_day(ic)  = heat_units(crop%Tbase(ic),climate%dtemp(ic)) / crop%PHU_maturity(ic)
       crop%fPHU(ic) = crop%fPHU(ic) + fPHU_day(ic)

       ! calculate SLA in m2 g-1 C
       SLA_C(ic) = SLA_development(crop%fPHU(ic),crop%sla_maturity(ic),crop%sla_beta(ic)) / DMtoC

       if (crop%vernalisation(ic) .and. .not. crop%vacc(ic)) then
         call vernalisation(climate,crop)
       endif

       if (crop%state(ic) == emergent) then

         call emergence(doy,SLA_C,fPHU_day,veg,casaflux,casapool,casamet,crop)

         
       else if (crop%state(ic) == growing) then



write(70,*) 'climate%dtemp: ', climate%dtemp
write(70,*) 'fPHU_day: ', fPHU_day          
write(70,*) 'crop%fPHU: ', crop%fPHU


         ! calculate maintenance and growth respiration
         !call casa_rplant(...,...) ! no need to call it here again, already done within biogeochem
         ! only need to calculate Rm and Rg here if we want to do it differently than in casa_rplant
         ! in that case, we might also change it in casa_rplant directly!!
         NPP_test(ic) = 0.6_dp * casaflux%Cgpp(ic)
         !casaflux%Crmplant(ic,:) = 0.05_dp * casaflux%Cgpp(ic)  
         !casaflux%Crgplant(ic)   = 0.2_dp * casaflux%Cgpp(ic)
         !casaflux%Cnpp(ic) = casaflux%Cgpp(ic) - sum(casaflux%Crmplant(ic,:)) - casaflux%Crgplant(ic)
         casaflux%Cnpp(ic) = casaflux%Cgpp(ic)
write(70,*) 'casaflux%Cgpp: ', casaflux%Cgpp
         ! calculate carbon allocation
write(60,*) 'doy:', doy
         call growth(SLA_C,veg,casaflux,casapool,casamet,crop)

write(70,*) 'casaflux%Cgpp: ', casaflux%Cgpp
write(70,*) 'casaflux%Cnpp: ', casaflux%Cnpp
write(70,*) 'NPP_test: ',      NPP_test
write(70,*) 'casamet%glai: ',  casamet%glai
         ! update root length and vegetation height
          
          
         ! harvest if enough PHU accumulated
         if (crop%fPHU(ic) >= 1.0_dp) then
            call harvest(doy,casapool,crop)

            ! reset heat units etc. (maybe better in harvest routine)
            crop%fPHU(ic) = 0.0_dp
            crop%VU(ic) = 0.0_dp
            crop%fVU(ic) = 0.0_dp
            crop%vacc(ic) = .FALSE.
         end if
       end if
     end if 
   end do
 
END SUBROUTINE crop_driver

