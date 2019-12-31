MODULE crop_module

  use cable_def_types_mod,  only: climate_type, soil_snow_type, soil_parameter_type, &
                                  dp => r_2
  use crop_def,             only: crop_type, nc, maxdays_ger, Cinit_root, Cinit_stem, &
                                  Cinit_leaf, DMtoC
  use casavariable,         only: casa_flux, casa_met, casa_pool
  use casaparm,             only: froot,wood,leaf,product, &   ! cpool
                                  metb,str,cwd                 ! litter
  
  implicit none

  private
  
  public :: sowing             ! determines sowing date
  public :: germination        ! determines germination date
  public :: emergence          ! plant growth at emergence
  public :: C_allocation_crops ! carbon allocation to roots, stems, leaves, products
  public :: harvest            ! removal of products and redistribution of remaining biomass
                               ! to litter and soil pools
  public :: heat_units         ! calculates heat units per day

  
contains


  subroutine sowing(doy,soil,crop)

    integer,                   intent(in)    :: doy  ! day of year
    type(soil_parameter_type), intent(in)    :: soil
    type(crop_type),           intent(inout) :: crop

    ! local
    integer :: i, sl  ! loop counters: crop type, soil layer
    real, dimension(size(soil%zse)) :: cumzse ! cumulative soil layer depths (bottom) (m)
    
    do i=1,nc
      if (doy == 280) then
         crop%state(i) = 1  ! sown
         crop%sowing_doy(i) = doy ! keep track of sowing doy

         ! if sowing occurred (i.e. crop%state switches to 1),
         ! determine soil layer in which seeds are located
         cumzse(:) = soil%zse(1)
         sl = 1
         do while (cumzse(sl) < crop%sowing_depth(i))
           cumzse(sl+1) = cumzse(sl) + soil%zse(sl+1)
           sl = sl + 1
           if (sl >= size(soil%zse)) then
              STOP "Implausible sowing depth value!"
           end if
         end do
         crop%sl = sl
      end if
    end do

      
  end subroutine sowing


  

  subroutine germination(doy,climate,ssnow,soil,crop)

    integer,                   intent(in)    :: doy ! day of year
    type(climate_type),        intent(in)    :: climate
    type(soil_snow_type),      intent(in)    :: ssnow
    type(soil_parameter_type), intent(in)    :: soil
    type(crop_type),           intent(inout) :: crop

    ! local
    integer  :: i       ! crop type
    integer,  dimension(nc) :: days_since_sowing
    real(dp), dimension(nc) :: fwsger  ! water stress factor for germination (0-1)
    real(dp), dimension(nc) :: fPHUger ! daily heat units for germination
    real(dp), dimension(nc) :: fgerday ! germination requirement for actual day

    ! initialise locals
    days_since_sowing = 0
    fwsger  = 0._dp
    fPHUger = 0._dp
    fgerday = 0._dp
    
    do i=1,nc
      ! calculate days since sowing 
      days_since_sowing(i) = doy - crop%sowing_doy(i) ! implement case year switch!!

      ! calculate heat requirements for germination (based on temp in soil layer crop%sl)
      fPHUger(i) = heat_units(crop%Tbase(i),climate%dtempsoil(i,crop%sl(i))) / &
                   crop%PHU_germination(i)
      
      ! calculate soil water requirements for germination
      fwsger(i) = water_stress_ger(ssnow%wb(i,crop%sl(i)),real(soil%ssat(i),dp), &
                  real(soil%swilt(i),dp))

      ! calculate total germination requirements (0-1)
      fgerday(i) = fPHUger(i) * fwsger(i)
      crop%fgermination(i) = crop%fgermination(i) + fgerday(i)

      ! based on the information calculated above, determine
      ! if germination occurred, if it is still ongoing,
      ! or if it has failed.
      if (days_since_sowing(i) > maxdays_ger) then ! germination assumed to have failed
         crop%state(i) = 0            ! fallow again, no plant growth
         crop%fgermination(i) = 0.0
      else if (crop%fgermination(i) >= 1.0_dp .AND. days_since_sowing(i) <= maxdays_ger) then
         crop%state(i) = 2           ! emerging
         crop%fgermination(i) = 0.0
      end if
write(60,*) 'doy:', doy
write(60,*) '  climate%dtempsoil(:,crop%sl(1)):', climate%dtempsoil(:,crop%sl(1))     
write(60,*) '  fPHUger:', fPHUger
write(60,*) '  fwsger:', fwsger
write(60,*) '  fgerday:', fgerday 
write(60,*) '  crop%fgermination:', crop%fgermination
write(60,*) '  crop%state:', crop%state
    end do
      
  end subroutine germination


  

  subroutine emergence(doy,fPHU_day,casaflux,casapool,casamet,crop)

    integer,  intent(in)           :: doy ! day of year
    real(dp), dimension(nc), intent(in) :: fPHU_day ! fPHU of actual day
    type(casa_flux), intent(inout) :: casaflux
    type(casa_pool), intent(inout) :: casapool
    type(casa_met),  intent(inout) :: casamet
    type(crop_type), intent(inout) :: crop

    ! local
    integer                 :: i ! crop type
    real(dp), dimension(nc) :: SLA_C
    real(dp), dimension(nc) :: LAIday

    do i=1,nc
       ! calculate SLA in m2 g-1 C
       SLA_C(i) = SLA_development(crop%fPHU(i),crop%sla_maturity(i),crop%sla_beta(i)) / DMtoC
write(60,*) 'doy:', doy
write(60,*) '  crop%sla_maturity(i):', crop%sla_maturity(i) 
write(60,*) '  crop%fPHU(i):', crop%fPHU(i)
write(60,*) '  SLA_C(i):', SLA_C(i)      
       ! initial C allocation
       casaflux%fracCalloc(i,froot)   = Cinit_root
       casaflux%fracCalloc(i,wood)    = Cinit_stem
       casaflux%fracCalloc(i,leaf)    = Cinit_leaf
       casaflux%fracCalloc(i,product) = 0.0_dp


       ! at emergence, utilize carbon reserves in the seeds (given by crop%Cseed),
       ! as well as carbon assimilated by first leaves
write(60,*) '  casaflux%Cnpp(i):', casaflux%Cnpp(i)
       casapool%dcplantdt(i,:) = casaflux%fracCalloc(i,:) * crop%Cseed(i) * fPHU_day(i)/crop%fPHU_emergence(i)
write(60,*) '  casapool%dcplantdt(i,:):', casapool%dcplantdt(i,:)
       casapool%dcplantdt(i,:) = casapool%dcplantdt(i,:) + (max(casaflux%Cnpp(i),0.0_dp) &
                                                          * casaflux%fracCalloc(i,:))
write(60,*) '  casapool%dcplantdt(i,:):', casapool%dcplantdt(i,:)
       casapool%Cplant(i,:) = casapool%Cplant(i,:) + casapool%dcplantdt(i,:)
write(60,*) '  casapool%Cplant(i,:):', casapool%Cplant(i,:)      
       
       ! initial LAI
write(60,*) '  casamet%glai:', casamet%glai
       LAIday(i) = casapool%dcplantdt(i,leaf) * SLA_C(i)
       casamet%glai(i) = casamet%glai(i) + LAIday(i)
       

write(60,*) '  casamet%glai:', casamet%glai
    end do
      
  end subroutine emergence





  subroutine C_allocation_crops(casaflux,casapool,casamet,crop) !! TBD: this could also be part of CASA!

    type(casa_flux), intent(inout) :: casaflux
    type(casa_pool), intent(inout) :: casapool
    type(casa_met),  intent(inout) :: casamet
    type(crop_type), intent(inout) :: crop

    ! local
    integer  :: i ! crop type
    real(dp), dimension(nc) :: SLA_C
    real(dp), dimension(nc) :: LAIday
       
    do i=1,nc
       SLA_C(i) = SLA_development(crop%fPHU(i),crop%sla_maturity(i),crop%sla_beta(i)) / DMtoC
       write(60,*) '  SLA_C(i):', SLA_C(i) 

       casaflux%fracCalloc(i,froot)   = 0.33_dp
       casaflux%fracCalloc(i,wood)    = 0.33_dp
       casaflux%fracCalloc(i,leaf)    = 0.33_dp
       casaflux%fracCalloc(i,product) = 0.01_dp

write(60,*) '  casaflux%Cnpp:', casaflux%Cnpp
write(60,*) '  casaflux%Cnpp/casaflux%Cgpp:', casaflux%Cnpp/casaflux%Cgpp
       ! update Cpools
       casapool%dcplantdt(i,:) = casaflux%Cnpp(i) * casaflux%fracCalloc(i,:)
       casapool%Cplant(i,:) = casapool%Cplant(i,:) + casapool%dcplantdt(i,:)

write(60,*) '  casapool%dcplantdt(i,:):', casapool%dcplantdt(i,:)
write(60,*) '  casapool%Cplant:', casapool%Cplant

       ! update LAI
       LAIday(i) = casapool%dcplantdt(i,leaf) * SLA_C(i)
       casamet%glai(i) = casamet%glai(i) + LAIday(i)
write(60,*) '  LAIday(i):', LAIday(i)
write(60,*) '  casamet%glai:', casamet%glai

       
    end do
    
  end subroutine C_allocation_crops



  
  
  subroutine harvest(doy,casapool,crop)

    integer, intent(in)             :: doy  ! day of year
    type(casa_pool),  intent(inout) :: casapool
    type (crop_type), intent(inout) :: crop 

    ! local
    integer                 :: i ! crop type
    real(dp), dimension(nc) :: Cstemleaf_removed 

    Cstemleaf_removed = 0.0_dp
write(60,*) 'doy:', doy    
    do i=1,nc
       ! remove product pool (aka yield)
       crop%yield(i) = casapool%Cplant(i,product)
       crop%harvest_index(i) = casapool%Cplant(i,product) / sum(casapool%Cplant(i,:)) ! diagnostic

write(60,*) '  casapool%Clitter(i,str)_BEFORE:', casapool%Clitter(i,str)
       ! add leaves and stems partly and roots completely to litter pool
       Cstemleaf_removed(i) = sum((1.0_dp - crop%Cplant_remove) * casapool%Cplant(i,leaf:wood)) 
       casapool%Clitter(i,str) = casapool%Clitter(i,str) + Cstemleaf_removed(i)  ! add to dt variable first?     
       casapool%Clitter(i,str) = casapool%Clitter(i,str) + casapool%Cplant(i,froot)

       ! set all plant pools to 0
       casapool%Cplant(i,:) = 0.0_dp
       
       ! change crop state to 0 (bare soil/inactive)
       crop%state(i) = 0  ! fallow

write(60,*) '  crop%yield(i):', crop%yield(i)
write(60,*) '  crop%harvest_index(i):', crop%harvest_index(i)
write(60,*) '  Cstemleaf_removed(i):', Cstemleaf_removed(i)
write(60,*) '  casapool%Clitter(i,str)_AFTER:', casapool%Clitter(i,str)
       
    end do   
    
  end subroutine harvest
    


  

  
  
  ! standard calculation of daily heat units
  function heat_units(basetemp,dtemp) result(HU)

    real, intent(in)  :: basetemp    ! base temperature (degC)
    real, intent(in)  :: dtemp       ! mean daily temperature (degC)
    real(dp)          :: HU          ! daily heat units
    
    HU = max(0.0, dtemp - max(0.0,basetemp)) 

  end function heat_units


  
  ! soil water stress factor for germination
  ! soil moisture (theta) is taken from soil layer at sowing depth
  function water_stress_ger(theta,thetas,thetaw) result(fwsger)

    real(dp), intent(in) :: theta  ! soil moisture at sowing depth
    real(dp), intent(in) :: thetas ! saturation soil moisture
    real(dp), intent(in) :: thetaw ! soil moisture at wilting point
    real(dp)             :: fwsger ! soil water stress factor for germination (0-1)

    fwsger = max(min((theta - thetaw) / (thetas - thetaw),1.0),0.0)
    
  end function water_stress_ger



  function SLA_development(fPHU,alpha,beta) result(SLA)

    real(dp), intent(in) :: fPHU   ! fraction of phenological heat units
    real(dp), intent(in) :: alpha  ! baseline value of SLA (at fPHU=1)
    real(dp), intent(in) :: beta   ! exponent (rate of decrease over developmental stage)
    real(dp)             :: SLA    ! specific leaf area (m2 gDM-1)

    SLA = alpha * fPHU**beta

  end function SLA_development
  
  
END MODULE crop_module

    
