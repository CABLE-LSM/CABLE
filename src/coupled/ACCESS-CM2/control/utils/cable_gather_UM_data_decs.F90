module cable_gather_UM_data_decs
   implicit none

   REAL, DIMENSION(:,:), allocatable, save :: &
    ls_rain_cable,    &
    ls_snow_cable,    &
    conv_rain_cable,    &
    conv_snow_cable

   REAL, DIMENSION(:,:,:), allocatable, save :: &
    surf_down_sw_cable
      
   REAL, DIMENSION(:), allocatable, save :: &
    bexp_cable,    &
    satcon_cable,    &
    sathh_cable,    &
    smvcst_cable

!   REAL, DIMENSION(:), allocatable :: &
!    canopy_cable,    &
!    gs_cable

Contains


subroutine cable_set_ls_precip( row_length, rows,                         &
                ls_rain, ls_snow )

  implicit none
  
  integer :: row_length,rows
  
  real, dimension(row_length,rows) :: &
    ls_rain,    &
    ls_snow

    if(.NOT. allocated(ls_rain_cable) ) &
      allocate( ls_rain_cable(row_length,rows),    &
                ls_snow_cable(row_length,rows) )
    
    ls_rain_cable = ls_rain
    ls_snow_cable = ls_snow

End subroutine cable_set_ls_precip

subroutine cable_set_conv_precip( row_length, rows,                         &
                conv_rain, conv_snow )

  implicit none
  
  integer :: row_length,rows
  
  real, dimension(row_length,rows) :: &
    conv_rain,    &
    conv_snow
    
    if(.NOT. allocated(conv_rain_cable) ) &
      allocate(conv_rain_cable(row_length,rows),    &
                conv_snow_cable(row_length,rows) )
    
    conv_rain_cable = conv_rain
    conv_snow_cable = conv_snow

End subroutine cable_set_conv_precip

subroutine cable_set_shortwave( row_length, rows, nbands,                     &
                surf_down_sw, first_call )

  implicit none
  
  integer :: row_length,rows, nbands
  
  real, dimension(row_length,rows,nbands) :: &
    surf_down_sw

  logical :: first_call
  
  if( .NOT. allocated(surf_down_sw_cable) ) & 
    allocate(surf_down_sw_cable(row_length,rows, nbands) )
  
  surf_down_sw_cable = surf_down_sw
  
  return

End subroutine cable_set_shortwave

subroutine cable_set_soil_params( land_pts, satcon, sathh, clapp, smvcst )

  implicit none
  
  integer :: land_pts 
  
  real, dimension(land_pts) :: &
   clapp, &
   satcon, &
   sathh, &
   smvcst 

    if(.NOT. allocated(bexp_cable) )  &
      allocate(                       &
        bexp_cable(land_pts),          &
        satcon_cable(land_pts),        &
        sathh_cable(land_pts),          &
        smvcst_cable(land_pts)          &
      )
    
    bexp_cable      = clapp
    satcon_cable    = satcon
    sathh_cable     = sathh
    smvcst_cable     = smvcst

End subroutine cable_set_soil_params

!subroutine cable_set_canopy( land_pts, canopy_gb )
!
!  implicit none
!  
!  integer :: land_pts 
!  
!  real, dimension(land_pts) :: &
!   canopy_gb
!
!    if(.NOT. allocated(canopy_cable) )   &
!      allocate(                    & 
!        canopy_cable(land_pts),    &
!      )
!    
!    canopy_cable  = canopy_gb
!
!End subroutine cable_set_canopy


End module cable_gather_UM_data_decs



