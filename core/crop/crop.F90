MODULE crop_module

  use cable_def_types_mod,  only: climate_type, soil_snow_type, soil_parameter_type, &
                                  veg_parameter_type, dp => r_2
  use crop_def,             only: crop_type, nc, maxdays_ger, fPHU_flowering, rCsen_Cgr, &
                                  fCleaf_mobile, DMtoC, baresoil, sown, emergent, growing
  use casavariable,         only: casa_flux, casa_met, casa_pool
  use casadimension,        only: mplant
  use casaparm,             only: froot,wood,leaf,product, &   ! cpool
                                  metb,str,cwd                 ! litter
  
  implicit none

  private
  
  public :: sowing             ! determines sowing date
  public :: germination        ! determines germination date
  public :: emergence          ! plant growth at emergence
  public :: growth             ! plant growth in vegetative and reproductive phase
  public :: senescence         ! plant senescence and effects on C pools and LAI + C remobilisation
  public :: harvest            ! removal of products and redistribution of remaining biomass
                               ! to litter and soil pools
  public :: vernalisation      ! calcualtes vernalisation factor and consequences for crop development
  public :: C_allocation_crops ! C allocation to plant pools (crops only)
  public :: heat_units         ! calculates heat units per day
  public :: SLA_development    ! calculates SLA in dependence of crop age/ developmental stage
  public :: irrigation         ! calculates irrigation requirements
  
contains


  subroutine sowing(i,doy,soil,crop)

    integer,                   intent(in)    :: i    ! crop type
    integer,                   intent(in)    :: doy  ! day of year
    type(soil_parameter_type), intent(in)    :: soil
    type(crop_type),           intent(inout) :: crop

    ! local
    integer :: sl  ! loop counters: crop type, soil layer
    real, dimension(size(soil%zse)) :: cumzse ! cumulative soil layer depths (bottom) (m)
    
    if (doy == 280) then
       crop%state(i) = sown
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

   
  end subroutine sowing


  

  subroutine germination(i,doy,climate,ssnow,soil,crop)

    integer,                   intent(in)    :: i   ! crop type
    integer,                   intent(in)    :: doy ! day of year
    type(climate_type),        intent(in)    :: climate
    type(soil_snow_type),      intent(in)    :: ssnow
    type(soil_parameter_type), intent(in)    :: soil
    type(crop_type),           intent(inout) :: crop

    ! local
    integer  :: days_since_sowing
    real(dp) :: fwsger  ! water stress factor for germination (0-1)
    real(dp) :: fPHUger ! daily heat units for germination
    real(dp) :: fgerday ! germination requirement for actual day

    ! initialise locals
    days_since_sowing = 0
    fwsger  = 0._dp
    fPHUger = 0._dp
    fgerday = 0._dp
    

    ! calculate days since sowing 
    days_since_sowing = doy - crop%sowing_doy(i) ! implement case year switch!!

    ! calculate heat requirements for germination (based on temp in soil layer crop%sl)
    fPHUger = heat_units(climate%dtempsoil(i,crop%sl(i)),crop%Tbase(i),crop%Tmax(i)) / &
              crop%PHU_germination(i)
      
    ! calculate soil water requirements for germination
    fwsger = water_stress_ger(ssnow%wb(i,crop%sl(i)),real(soil%ssat(i),dp), &
             real(soil%swilt(i),dp))

    ! calculate total germination requirements (0-1)
    fgerday = fPHUger * fwsger
    crop%fgermination(i) = crop%fgermination(i) + fgerday

    ! based on the information calculated above, determine
    ! if germination occurred, if it is still ongoing,
    ! or if it has failed.
    if (days_since_sowing > maxdays_ger) then ! germination assumed to have failed
       crop%state(i) = baresoil           ! bare soil again, no plant growth
       crop%fgermination(i) = 0.0
    else if (crop%fgermination(i) >= 1.0_dp .AND. days_since_sowing <= maxdays_ger) then
       crop%state(i) = emergent           ! successful germination
       crop%fgermination(i) = 0.0
    end if
write(60,*) 'doy:', doy
write(60,*) '  climate%dtempsoil(:,crop%sl(1)):', climate%dtempsoil(:,crop%sl(1))     
write(60,*) '  fPHUger:', fPHUger
write(60,*) '  fwsger:', fwsger
write(60,*) '  fgerday:', fgerday 
write(60,*) '  crop%fgermination:', crop%fgermination
write(60,*) '  crop%state:', crop%state
      
  end subroutine germination


  

  subroutine emergence(i,doy,SLA_C,fPHU_day,veg,casaflux,casapool,casamet,crop)

    integer,  intent(in)    :: i       ! crop type
    integer,  intent(in)    :: doy      ! day of year
    real(dp), intent(in)    :: SLA_C
    real(dp), intent(inout) :: fPHU_day ! fPHU of actual day
    type(veg_parameter_type),intent(inout) :: veg
    type(casa_flux), intent(inout) :: casaflux
    type(casa_pool), intent(inout) :: casapool
    type(casa_met),  intent(inout) :: casamet
    type(crop_type), intent(inout) :: crop

    ! local
    real(dp) :: LAIday


write(60,*) 'doy:', doy
write(60,*) '  crop%fPHU(i):', crop%fPHU(i)
write(60,*) '  SLA_C:', SLA_C     

    ! at emergence, utilize carbon reserves in the seeds (given by crop%Cseed),
    ! as well as carbon assimilated by first leaves
       
    ! first, make sure C allocation from seeds does not exceed crop%Cseed
    ! also update crop%state if fPHU is high enough
    if (crop%fPHU(i) >= crop%fPHU_emergence(i)) then  
       fPHU_day = fPHU_day - (crop%fPHU(i) - crop%fPHU_emergence(i))
       crop%state(i) = growing
    endif
      
write(60,*) '  casaflux%Cnpp(i):', casaflux%Cnpp(i)
    casapool%dcplantdt(i,:) = casaflux%fracCalloc(i,:) * crop%Cseed(i) * fPHU_day/crop%fPHU_emergence(i)
write(60,*) '  casapool%dcplantdt(i,:):', casapool%dcplantdt(i,:)
    casapool%dcplantdt(i,:) = casapool%dcplantdt(i,:) + (max(casaflux%Cnpp(i),0.0_dp) &
                                                          * casaflux%fracCalloc(i,:))
write(60,*) '  casapool%dcplantdt(i,:):', casapool%dcplantdt(i,:)
    casapool%Cplant(i,:) = casapool%Cplant(i,:) + casapool%dcplantdt(i,:)
write(60,*) '  casapool%Cplant(i,:):', casapool%Cplant(i,:)

    ! update Creserve pool in stems
    crop%Cstem_mobile(i) = crop%Cstem_mobile(i) + crop%fCstem_mobile(i) * casapool%dcplantdt(i,wood) 
       
    ! update initial LAI
write(60,*) '  casamet%glai:', casamet%glai
    LAIday = casapool%dcplantdt(i,leaf) * SLA_C
    casamet%glai(i) = casamet%glai(i) + LAIday
       
write(60,*) '  casamet%glai:', casamet%glai
    ! update vegetation height
    veg%hc(i) = vegheight_Cstem_Cleaf(casapool%Cplant(i,wood), &
                                         casapool%Cplant(i,leaf)) 
write(60,*) '  veg%hc(i):', veg%hc(i)
    ! update rooting depth and root distribution
      
  end subroutine emergence

  

  subroutine growth(i,SLA_C,veg,casaflux,casapool,casamet,crop)

    integer,                   intent(in)    :: i ! crop type
    real(dp),                  intent(in)    :: SLA_C
    type(veg_parameter_type),  intent(inout) :: veg
    type(casa_flux),           intent(inout) :: casaflux
    type(casa_pool),           intent(inout) :: casapool
    type(casa_met),            intent(inout) :: casamet
    type(crop_type),           intent(inout) :: crop

    ! local
    real(dp) :: LAIday

write(60,*) '  crop%fPHU(i):', crop%fPHU(i)
write(60,*) '  casaflux%fracCalloc(i,froot):', casaflux%fracCalloc(i,froot)
write(60,*) '  casaflux%fracCalloc(i,leaf):', casaflux%fracCalloc(i,leaf)
write(60,*) '  casaflux%fracCalloc(i,wood):', casaflux%fracCalloc(i,wood)
write(60,*) '  casaflux%fracCalloc(i,product):', casaflux%fracCalloc(i,product)
write(60,*) '  casaflux%Cgpp:', casaflux%Cgpp
write(60,*) '  casaflux%Cnpp:', casaflux%Cnpp
write(60,*) '  casaflux%Cnpp/casaflux%Cgpp:', casaflux%Cnpp/casaflux%Cgpp

    ! update Cpools
    casapool%dcplantdt(i,:) = casaflux%Cnpp(i) * casaflux%fracCalloc(i,:)
    casapool%Cplant(i,:) = casapool%Cplant(i,:) + casapool%dcplantdt(i,:)

    ! update C reserve pool in stems
    crop%Cstem_mobile(i) = crop%Cstem_mobile(i) + crop%fCstem_mobile(i) * casapool%dcplantdt(i,wood)
       
write(60,*) '  casapool%dcplantdt(i,:):', casapool%dcplantdt(i,:)
write(60,*) '  casapool%Cplant:', casapool%Cplant
write(60,*) '  crop%Cstem_mobile(i):', crop%Cstem_mobile(i)

    if (casapool%dcplantdt(i,leaf) > 0.0_dp) then
      ! update LAI
      LAIday = casapool%dcplantdt(i,leaf) * SLA_C
      casamet%glai(i) = casamet%glai(i) + LAIday
write(60,*) '  LAIday:', LAIday
write(60,*) '  casamet%glai:', casamet%glai
    endif
    if (casapool%dcplantdt(i,leaf) > 0.0_dp .or. casapool%dcplantdt(i,wood) > 0.0_dp) then
      ! update vegetation height
      veg%hc(i) = vegheight_Cstem_Cleaf(casapool%Cplant(i,wood), &
                                        casapool%Cplant(i,leaf)) 
write(60,*) '  veg%hc(i):', veg%hc(i)
    endif

  end subroutine growth


  
  subroutine C_allocation_crops(i,casaflux,crop) !! TBD: this could also be part of CASA!

    integer,         intent(in)    :: i ! crop type
    type(casa_flux), intent(inout) :: casaflux
    type(crop_type), intent(inout) :: crop


    ! determine C allocation fractions
    casaflux%fracCalloc(i,froot)   = Calloc_roots(crop%fPHU(i),crop%fCalloc_root_init(i),crop%Calloc_root_end(i))            
    casaflux%fracCalloc(i,leaf)    = Calloc_leaves(crop%fPHU(i),crop%Calloc_leaf_bpt(i),crop%Calloc_leaf_end(i), &
                                                   crop%fCalloc_leaf_init(i),crop%fCalloc_leaf_bpt(i))
    casaflux%fracCalloc(i,product) = Calloc_products(crop%fPHU(i),fPHU_flowering,crop%Calloc_prod_max(i))
    casaflux%fracCalloc(i,wood)    = 1.0_dp - (casaflux%fracCalloc(i,froot) + casaflux%fracCalloc(i,leaf) + &
                                               casaflux%fracCalloc(i,product))

    !if (crop%dynamic_allocation) then
    !  write(*,*) 'not yet implemented!'
    !endif
    
  end subroutine C_allocation_crops



  subroutine senescence(i,fPHU_day,casamet,casapool,casaflux,crop)
    ! models sequential senescence (constant leaf senescence due to e.g. ageing) and
    ! reproductive senescence (occurs after flowering)

    integer,         intent(in)    :: i ! crop type
    real(dp),        intent(in)    :: fPHU_day
    type(casa_met),  intent(inout) :: casamet
    type(casa_pool), intent(inout) :: casapool
    type(casa_flux), intent(inout) :: casaflux
    type(crop_type), intent(inout) :: crop

    ! local
    integer  :: p            ! plant pool
    real(dp) :: sen_period   ! period (in fPHU) over which senescence occurs
    real(dp) :: rsenesc      ! senescence rate leaves (0-1)
    real(dp) :: C_loss_leaf  ! carbon loss caused by senescence
    real(dp) :: C_loss_mobil ! 'lost' carbon remobilised to product pool
    real(dp) :: C_loss_resp  ! 'lost' respired directly from C pools
    

    ! model reproductive senescence as soon as C allocation to leaves ceases
    ! senescence refers to leaf pool only. Senescence in stems is only modelled
    ! as remobilisation of C to products, senescence in roots is not modelled.
    ! --> at the moment stems and roots do not experience mass loss due to respiration
    ! at senescence, i.e. any respiration C is coming from GPP!
    ! currently assumes that senescence stops at harvest, not before
    sen_period = 1.0_dp - crop%Calloc_leaf_end(i)
    if (crop%fPHU(i) > crop%Calloc_leaf_end(i)) then
       rsenesc = ((crop%fPHU(i) - crop%Calloc_leaf_end(i)) / &
                   sen_period)**crop%drsen(i) - &
                   ((max(crop%fPHU(i) - fPHU_day,crop%Calloc_leaf_end(i)) - crop%Calloc_leaf_end(i)) / &
                   sen_period)**crop%drsen(i)
    else
       rsenesc = 0.0_dp    
    endif
    ! update complete canopy senescence fraction
    ! avoid fsenesc exceeding 1.0
    if (rsenesc + crop%fsenesc(i) > 1.0_dp) then
       rsenesc = rsenesc - ((rsenesc + crop%fsenesc(i)) - 1.0_dp)
    endif
    crop%fsenesc(i) = min(crop%fsenesc(i) + rsenesc,1.0_dp)
write(60,*) '  crop%fsenesc(i):', crop%fsenesc(i)
    ! update dead and green LAI
    crop%LAItot(i)  = max(0.0_dp,casamet%glai(i) + crop%LAIsen(i))
    crop%LAIsen(i)  = max(0.0_dp,crop%fsenesc(i) * crop%LAItot(i))
    casamet%glai(i) = max(0.0_dp,crop%LAItot(i) - crop%LAIsen(i))

    do p=1,mplant
      if (casaflux%fracCalloc(i,p) > 1.0e-9) then
          crop%Cmax(i,p) = casapool%Cplant(i,p)
          if (p==wood) then
            crop%Cmaxstemmob(i) = crop%Cstem_mobile(i)
          endif
       endif
    end do

    ! update Rm and NPP given fsenesc (NPP is updated later in this subroutine)
    ! fsenesc calculated for leaves, but represents overall senescence state of
    ! the crop in this case, to be improved
write(60,*) '  casaflux%Crmplant1(i,:):', casaflux%Crmplant(i,:)
    casaflux%Crmplant(i,:) = casaflux%Crmplant(i,:) * (1.0_dp - crop%fsenesc(i))
    casaflux%Cnpp(i) = casaflux%Cgpp(i) - sum(casaflux%Crmplant(i,:)) - casaflux%Crgplant(i)
write(60,*) '  casaflux%Crmplant2(i,:):', casaflux%Crmplant(i,:)

    ! Calculate C loss in leaves due to senescence 
    ! Part of this C loss is respired, the other part is mobilised to the product 
    C_loss_leaf = rsenesc * (1.0_dp - rCsen_Cgr) * crop%Cmax(i,leaf)
    C_loss_resp = fCleaf_mobile * C_loss_leaf
    if (C_loss_resp > casaflux%Crmplant(i,leaf)) then
       C_loss_resp = casaflux%Crmplant(i,leaf)
    endif
    C_loss_mobil = C_loss_leaf - C_loss_resp
write(60,*) '  C_loss_leaf:', C_loss_leaf
       
    call remobilisation(i,casapool,crop,fPHU_day,C_loss_mobil)

    ! update NPP: usually, C for Rm is taken from GPP, but in senescence it is taken directly from the C
    ! pools in the leaves/stems. Thus updating NPP to conserve C mass balance.
    casaflux%Cnpp(i) = casaflux%Cnpp(i) + (sum(casaflux%Crmplant(i,(/leaf,wood/))) - C_loss_resp)
write(60,*) '  sum(casaflux%Crmplant(i,(/leaf,wood/))):', sum(casaflux%Crmplant(i,(/leaf,wood/)))
write(60,*) '  casaflux%Cnpp2:', casaflux%Cnpp    
    
  end subroutine senescence 



  
  subroutine remobilisation(i,casapool,crop,fPHU_day,Closs_leaf_mobile)

    integer,         intent(in)    :: i       ! crop type
    type(casa_pool), intent(inout) :: casapool
    type(crop_type), intent(inout) :: crop
    real(dp),        intent(in)    :: fPHU_day
    real(dp),        intent(in)    :: Closs_leaf_mobile ! C loss of leaves due to remobilisation

    ! local
    real(dp) :: remob_period   ! period over which remobilisation occurs
    real(dp) :: rCloss_stem    ! rate of Carbon loss of mobile reserves in stems
    real(dp) :: Closs_stem     ! C loss of stem due to remobilisation


    ! --------------------
    !  C remobilisation  !
    ! --------------------   

    ! 1) from stem to products
    !    for stems, it is assumed that all of their mobile reserves are
    !    relocated to the products during the first half of senescence
    remob_period = (1.0_dp - crop%Calloc_prod_max(i)) * 0.5_dp
    if (crop%fPHU(i) > crop%Calloc_prod_max(i) .and. crop%fPHU(i) < &
         (crop%Calloc_prod_max(i) + remob_period)) then          
       rCloss_stem = ((crop%fPHU(i) - crop%Calloc_prod_max(i)) / &
                          remob_period)**0.5_dp - &
                          ((max(crop%fPHU(i) - fPHU_day,crop%Calloc_prod_max(i)) &
                           - crop%Calloc_prod_max(i)) / &
                           remob_period)**0.5_dp
    else
       rCloss_stem = 0.0_dp
    endif
    ! avoid fsenesc exceeding 1.0 
    if (rCloss_stem + crop%fmobilise(i) > 1.0_dp) then
        rCloss_stem = rCloss_stem - ((rCloss_stem + crop%fmobilise(i)) - 1.0_dp)
    endif
    crop%fmobilise(i) = min(crop%fmobilise(i) + rCloss_stem,1.0_dp)

    ! subtract C_loss from mobile reserves (not ideal right now, because
    ! C has to be subtracted from both Cstem_mobile and Cpools(i,wood)...
    Closs_stem = rCloss_stem * crop%Cmaxstemmob(i)
    crop%Cstem_mobile(i)    = crop%Cstem_mobile(i) - Closs_stem
    casapool%Cplant(i,wood) = casapool%Cplant(i,wood) - Closs_stem

    ! add carbon lost from stem to product pool
    casapool%Cplant(i,product) = casapool%Cplant(i,product) + Closs_stem

      
      
    ! 2) from leaves to products
    !    a fraction corresponding to fCleaf_mobile is moved from leaves
    !    to stems during leaf senescence (calculated in senescence subroutine)
    casapool%Cplant(i,leaf) = casapool%Cplant(i,leaf) - Closs_leaf_mobile
    ! add carbon lost from stem to product pool
    casapool%Cplant(i,product) = casapool%Cplant(i,product) + Closs_leaf_mobile

      
    ! --------------------
    !  N remobilisation  !
    ! --------------------

  
    
  end subroutine remobilisation
  

  
  !----------------------------
  ! Management subroutines ----
  !----------------------------  
  subroutine irrigation(i,veg,ssnow,soil,crop)
    ! irrigate whenever soil moisture falls below field capacity

    integer,                   intent(in)    :: i       ! crop type
    type(veg_parameter_type),  intent(in)    :: veg
    type(soil_snow_type),      intent(in)    :: ssnow
    type(soil_parameter_type), intent(in)    :: soil
    type(crop_type),           intent(inout) :: crop

    ! local
    real(dp) :: missing
write(50,*) 'ssnow%wb:', ssnow%wb
write(50,*) 'soil%ssat,soil%sfc,soil%swilt:', soil%ssat,soil%sfc,soil%swilt   

    ! check where roots are growing
    !veg%zr()

! 2) determine missing water for this region
! 3) check results: fwsoil should be close to 1, runoff should not be significantly increased
    missing = (soil%sfc(1) * sum(soil%zse(1:3))*1000._dp) - sum(ssnow%wb(1,1:3))    ! missing water in mm
    
write(50,*) 'missing:', missing      
    
  end subroutine irrigation


  !subroutine fertilisation

  !end subroutine fertilisation

  
  subroutine harvest(i,doy,casapool,casamet,veg,crop)

    integer,                  intent(in)    :: i       ! crop type
    integer,                  intent(in)    :: doy  ! day of year
    type(casa_pool),          intent(inout) :: casapool
    type(casa_met),           intent(inout) :: casamet
    type(veg_parameter_type), intent(inout) :: veg
    type(crop_type),          intent(inout) :: crop

    ! local
    real(dp) :: Cstemleaf_removed 

    Cstemleaf_removed = 0.0_dp
write(60,*) 'doy:', doy    

    ! remove product pool (aka yield)
    crop%yield(i) = casapool%Cplant(i,product)
    crop%harvest_index(i) = casapool%Cplant(i,product) / sum(casapool%Cplant(i,(/leaf,wood,product/))) ! diagnostic

write(60,*) '  casapool%Clitter(i,str)_BEFORE:', casapool%Clitter(i,str)
    ! add leaves and stems partly and roots completely to litter pool
    Cstemleaf_removed = sum((1.0_dp - crop%Cplant_remove) * casapool%Cplant(i,leaf:wood)) 
    casapool%Clitter(i,str) = casapool%Clitter(i,str) + Cstemleaf_removed  ! add to dt variable first?     
    casapool%Clitter(i,str) = casapool%Clitter(i,str) + casapool%Cplant(i,froot)

    ! set all plant pools to 0
    casapool%Cplant(i,:) = 0.0_dp

    ! set LAI to 0
    casamet%glai(i) = 0.0_dp
    crop%LAItot(i)  = 0.0_dp
    crop%LAIsen(i)  = 0.0_dp

    ! set plant height to 0
    veg%hc(i) = 0.0_dp
       
    ! change crop state to 0 (bare soil/inactive)
    crop%state(i) = baresoil

write(60,*) '  crop%yield(i):', crop%yield(i)
write(60,*) '  crop%harvest_index(i):', crop%harvest_index(i)
write(60,*) '  Cstemleaf_removed:', Cstemleaf_removed
write(60,*) '  casapool%Clitter(i,str)_AFTER:', casapool%Clitter(i,str)   
    
  end subroutine harvest
    


  
  subroutine vernalisation(i,climate,crop)

    integer,             intent(in)        :: i       ! crop type
    type (climate_type), intent(in)        :: climate
    type (crop_type),    intent(inout)     :: crop 

    ! local
    real :: VU_day ! daily accumulated vernalisation unit

    VU_day = vernalisation_rate(climate%dtemp(i))
    crop%VU(i) = crop%VU(i) + real(VU_day,dp) 
    crop%fVU(i) = min((crop%VU(i)**5.0_dp) / (22.5_dp**5.0_dp + crop%VU(i)**5.0_dp),1.0_dp)
write(60,*) ' crop%fVU(i):', crop%fVU(i)
      
    if (crop%fPHU(i) >= 0.45_dp) then
       crop%fPHU(i) = crop%fPHU(i) * crop%fVU(i) ! alternative: increase fPHU_maturity.
       ! also affects C allocation to grains!!
       crop%vacc(i) = .TRUE.        
    endif
      
  end subroutine vernalisation


  

  !-----------------------------
  ! FUNCTIONS
  !-----------------------------

  ! calculation of a daily vernalisation rate (required for some crops)
  ! following Streck et al. 2003
  ! note typo in fvn function in Streck et al. 2003. Right formulation
  ! can be found in e.g. Lu et al. 2017.
  function vernalisation_rate(dtemp) result(fvn)

    real, intent(in) :: dtemp   ! mean daily temperature (degC)
    
    real, parameter :: Tmin=-1.3  ! minimum temperature for vernalisation
    real, parameter :: Topt=4.9   ! optimum temperature for vernalisation
    real, parameter :: Tmax=15.7  ! maximum temperature for vernalisation

    real :: alpha  ! used in for fvn calculation
    real :: fvn    ! daily vernalisation rate

    alpha = log(2.0) / log((Tmax - Tmin)/(Topt - Tmin))

    if (dtemp >= Tmin .and. dtemp <= Tmax) then
       fvn = (2.0*(dtemp - Tmin)**alpha * (Topt-Tmin)**alpha - &
             (dtemp - Tmin)**(2.0*alpha)) / &
             ((Topt - Tmin)**(2.0*alpha))
    else
       fvn = 0.0
    endif

  end function vernalisation_rate
  
  
  ! standard calculation of daily heat units
  function heat_units(dtemp,basetemp,maxtemp) result(HU)

    real, intent(in)  :: dtemp       ! mean daily temperature (degC)
    real, intent(in)  :: basetemp    ! base temperature (degC)
    real, intent(in)  :: maxtemp     ! maximum mean daily temp. (degC)
                                     ! no further increases in HU above maxtemp
    real(dp)          :: HU          ! daily heat units
    
    HU = min(max(0.0, dtemp - max(0.0,basetemp)),maxtemp)

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

  

  ! decrease of SLA with crop age
  function SLA_development(fPHU,alpha,beta) result(SLA)

    real(dp), intent(in) :: fPHU   ! fraction of phenological heat units
    real(dp), intent(in) :: alpha  ! baseline value of SLA (at fPHU=1)
    real(dp), intent(in) :: beta   ! exponent (rate of decrease over developmental stage)
    real(dp)             :: SLA    ! specific leaf area (m2 gDM-1)

    SLA = alpha * fPHU**beta

  end function SLA_development



  
  ! Functions describing the C allocation to roots, leaves, and products
  ! as a function of the phenological development
  ! Roots
  function Calloc_roots(fPHU,fCalloc_init,Calloc_end) result(fCalloc)

    real(dp), intent(in) :: fPHU         ! fraction of phenological heat units
    real(dp), intent(in) :: fCalloc_init ! initial fraction of C allocation (at fPHU=0)
    real(dp), intent(in) :: Calloc_end   ! developmental stage (fPHU) at which fCalloc reaches 0
    real(dp)             :: fCalloc

    if (fPHU <= Calloc_end) then
       fCalloc = fCalloc_init - 1.0_dp/Calloc_end*fCalloc_init * fPHU
    else
       fCalloc = 0.0_dp
    endif

  end function Calloc_roots


  ! Leaves
  function Calloc_leaves(fPHU,fPHU_bpt,Calloc_end,fCalloc_init,fCalloc_bpt) result(fCalloc)

    real(dp), intent(in) :: fPHU           ! fraction of phenological heat units
    real(dp), intent(in) :: fPHU_bpt       ! fPHU at breakpoint
    real(dp), intent(in) :: Calloc_end     ! developmental stage (fPHU) at which fCalloc reaches 0
    real(dp), intent(in) :: fCalloc_init   ! initial fraction of C allocation (at fPHU=0)
    real(dp), intent(in) :: fCalloc_bpt    ! fCalloc to leaves at breakpoint
    real(dp)             :: fCalloc
    

    if (fPHU <= fPHU_bpt) then
       fCalloc = fCalloc_init - 1.0_dp/fPHU_bpt*(fCalloc_init - fCalloc_bpt) * fPHU
    else if (fPHU > fPHU_bpt .and. fPHU < Calloc_end) then
       fCalloc = fCalloc_bpt - fCalloc_bpt/(Calloc_end - fPHU_bpt) * (fPHU - fPHU_bpt)
    else
       fCalloc = 0.0_dp
    endif

  end function Calloc_leaves

  
  ! Products
  function Calloc_products(fPHU,fPHU_flowering,Calloc_max) result(fCalloc)

    real(dp), intent(in) :: fPHU              ! fraction of phenological heat units
    real(dp), intent(in) :: fPHU_flowering    ! fPHU at flowering (usually 0.5)
    real(dp), intent(in) :: Calloc_max        ! developmental stage (fPHU) at which fCalloc reaches 1
    real(dp)             :: fCalloc

    if (fPHU <= fPHU_flowering) then
       fCalloc = 0.0_dp
    else if (fPHU > fPHU_flowering .and. fPHU < Calloc_max) then
       fCalloc = 1.0_dp / (Calloc_max - fPHU_flowering) * (fPHU - fPHU_flowering)
    else
       fCalloc = 1.0_dp
    endif

  end function Calloc_products
  

  
  ! Crop height from allometric relationships
  ! Arora and Boer 2005 (Eq. A2): dependence on Cstem and Cleaf
  function vegheight_Cstem_Cleaf(Cstem,Cleaf) result(height)

    real(dp), intent(in) :: Cstem  ! C content in stem (gC m-2)
    real(dp), intent(in) :: Cleaf  ! C content in leaves (gC m-2)
    real(dp)             :: height ! crop height (m)

    height = ((Cstem + Cleaf)/1000.0_dp)**0.385_dp 

  end function vegheight_Cstem_Cleaf
  
  
END MODULE crop_module

    
