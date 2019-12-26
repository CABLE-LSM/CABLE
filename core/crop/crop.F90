MODULE crop_module

  use cable_def_types_mod,  only: climate_type, soil_snow_type, soil_parameter_type, &
                                  dp => r_2
  use crop_def,             only: crop_type, nc, maxdays_ger
  use casavariable,         only: casa_flux
  
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

      ! calculate heat requirements for germination
      fPHUger(i) = heat_units(crop%Tbase(i),climate%dtemp(i)) / crop%PHU_germination(i)
      !dtemp needs to be soil temp!!!!!
      
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
        crop%state(i) = 0  ! fallow again, no plant growth
        crop%fgermination(i) = 0.0
      else if (crop%fgermination(i) >= 1.0_dp .AND. days_since_sowing(i) <= maxdays_ger) then
        crop%state(i) = 2  ! emerging
     end if
write(60,*) 'doy:', doy     
write(60,*) '  fPHUger:', fPHUger
write(60,*) '  fwsger:', fwsger
write(60,*) '  fgerday:', fgerday 
write(60,*) '  crop%fgermination:', crop%fgermination
    end do
      
  end subroutine germination


  

  subroutine emergence(doy,crop)

    integer, intent(in)             :: doy ! day of year
    type (crop_type), intent(inout) :: crop

    ! local
    integer :: i ! crop type

    do i=1,nc
      ! seed/plant density

      ! initial C allocation

      ! initial LAI

      ! update state
      crop%state(i) = 3  ! growing
    end do
      
  end subroutine emergence





  subroutine C_allocation_crops(casaflux,crop)

    type (casa_flux), intent(inout) :: casaflux
    type (crop_type), intent(inout) :: crop

    ! local
    integer :: i ! crop type

    do i=1,nc
      write(*,*) 'nothing here yet'
    end do
    
  end subroutine C_allocation_crops



  
  
  subroutine harvest(doy,crop)

    integer, intent(in)             :: doy  ! day of year
    type (crop_type), intent(inout) :: crop 

    ! local
    integer :: i ! crop type

    do i=1,nc
      ! re-allocate biomass to soil or remove
      crop%state(i) = 0  ! fallow
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
  
END MODULE crop_module

    
