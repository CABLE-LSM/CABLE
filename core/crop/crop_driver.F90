SUBROUTINE crop_driver(doy,climate,ssnow,soil,veg,canopy,casaflux,casamet, &
                       casapool,crop)

  use crop_def,            only: crop_type, nc, baresoil, planted, emergent, growing, &
                                 Rgcoeff, DMtoC
  use crop_module
  use casavariable,        only: casa_flux, casa_met, casa_pool
  use cable_def_types_mod, only: climate_type, soil_snow_type, soil_parameter_type, &
                                 veg_parameter_type, canopy_type, dp => r_2

  implicit none

  integer,                   intent(in)    :: doy
  type(climate_type),        intent(in)    :: climate
  type(soil_snow_type),      intent(in)    :: ssnow
  type(soil_parameter_type), intent(in)    :: soil
  type(veg_parameter_type),  intent(inout) :: veg
  type(canopy_type),         intent(inout) :: canopy
  type(casa_flux),           intent(inout) :: casaflux
  type(casa_met),            intent(inout) :: casamet
  type(casa_pool),           intent(inout) :: casapool
  type(crop_type),           intent(inout) :: crop
  
  ! local
  integer  :: ic,sl  ! loop counters: crop type, soil layer
  real(dp) :: fPHU_day
  real(dp) :: SLA_C
  
  
  do ic=1, nc
write(70,*) 'DOY: ', doy
    if (crop%state(ic) == baresoil) then
      if (doy >= crop%sowing_doymin(ic) .and. doy <= crop%sowing_doymax(ic)) then 
         call planting(ic,doy,climate,ssnow,soil,crop)

      end if

    else if (crop%state(ic) == planted) then
       call germination(ic,doy,climate,ssnow,soil,crop)
       call irrigation(ic,ssnow,soil,canopy)
casapool%Cplant(ic,:) = 0.0_dp ! shouldn't be needed here!! Check initialisation
casamet%glai(ic) = 0.0_dp

    else if (crop%state(ic) == emergent .or. crop%state(ic) == growing) then

       ! update phenological heat units (start at 0 at germination!)
       fPHU_day = heat_units(climate%dtemp(ic),crop%Tbase(ic),crop%Tmax(ic)) / crop%PHU_maturity(ic)
       crop%fPHU(ic) = crop%fPHU(ic) + fPHU_day

       ! calculate SLA in m2 g-1 C
       SLA_C = SLA_development(crop%fPHU(ic),crop%sla_maturity(ic),crop%sla_beta(ic)) / DMtoC

       if (crop%vernalisation(ic) .and. .not. crop%vacc(ic)) then
         call vernalisation(ic,climate,crop)
       endif

       ! calculate C allocation factors
       call C_allocation_crops(ic,casaflux,crop)

       ! calculate growth respiration
       casaflux%Crgplant(ic) = sum(Rgcoeff * casaflux%fracCalloc(ic,:)) * casaflux%Cgpp(ic)
write(65,*) 'sum(Rgcoeff * casaflux%fracCalloc(ic,:))', sum(Rgcoeff * casaflux%fracCalloc(ic,:))
write(65,*) 'casaflux%fracCalloc(ic,:)', casaflux%fracCalloc(ic,:)      
       ! calculate maintenance respiration
       ! as in casa at the moment. Evtl calculate Rm first, then Rg. Discuss
       ! if it makes sense to replace Rgcoeff with Ygrowth as calculated in casa_cnp
       ! only need to calculate Rm and Rg here if we want to do it differently than in casa_rplant
       ! in that case, we might also change it in casa_rplant directly!!

       ! calculate NPP
       casaflux%Cnpp(ic) = casaflux%Cgpp(ic) - sum(casaflux%Crmplant(ic,:)) - casaflux%Crgplant(ic)
       

       if (crop%state(ic) == emergent) then

         call emergence(ic,doy,SLA_C,fPHU_day,veg,casaflux,casapool,casamet,crop)

       else if (crop%state(ic) == growing) then

write(60,*) 'doy:', doy
write(60,*) '  casaflux%Cgpp(ic):', casaflux%Cgpp(ic)
write(60,*) '  sum(casaflux%Crmplant(ic,:)): ', sum(casaflux%Crmplant(ic,:))
write(60,*) '  casaflux%Crgplant(ic): ', casaflux%Crgplant(ic)
write(60,*) '  casaflux%Cnpp_first(ic):', casaflux%Cnpp(ic)
write(70,*) 'climate%dtemp: ', climate%dtemp
write(70,*) 'fPHU_day: ', fPHU_day          
write(70,*) 'crop%fPHU: ', crop%fPHU


         call senescence(ic,fPHU_day,casamet,casapool,casaflux,crop)

         call growth(ic,SLA_C,veg,casaflux,casapool,casamet,crop)

write(70,*) 'casaflux%Cgpp: ', casaflux%Cgpp
write(70,*) 'casaflux%Cnpp: ', casaflux%Cnpp
write(70,*) 'casamet%glai: ',  casamet%glai
          
         ! harvest if enough PHU accumulated
         if (crop%fPHU(ic) >= 1.0_dp) then
            call harvest(ic,doy,casapool,casamet,veg,crop)

            ! reset heat units and vernalisation requirements etc.
            crop%fPHU(ic)    = 0.0_dp
            crop%fsenesc(ic) = 0.0_dp
            crop%VU(ic)      = 0.0_dp
            crop%fVU(ic)     = 0.0_dp
            crop%vacc(ic)    = .FALSE.

            ! reset management settings
            canopy%irrig_surface(ic)   = 0.0
            canopy%irrig_sprinkler(ic) = 0.0

            ! reset other variables
            crop%state_nrdays(ic,:) = 0
            crop%sl(ic) = 0
            
         end if
       end if ! crop%state == growing

         call irrigation(ic,ssnow,soil,canopy)
      
     end if 
   end do
 
END SUBROUTINE crop_driver

