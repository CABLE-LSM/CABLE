MODULE cable_psm

   USE cable_def_types_mod, only : r_2,ms,mp,air_type,met_type,soil_snow_type,&
                                  canopy_type,soil_parameter_type,veg_parameter_type,&
                                  roughness_type
   USE cable_common_module, only : cable_user

implicit none


   REAL(r_2), parameter :: Dff=2.5e-5, &  !diffusivity water vapor in air
                      lm=1.73e-5, &       !converts units
                      pi = 3.14159265358979324, &  !obvous
                      c2 = 2.0,&                  !params
                      litter_thermal_diff=8.3e-6  !param based on vh thermal diffusivity

   real(r_2), parameter :: rtevap_max = 10000.0

PUBLIC  or_soil_evap_resistance,update_or_soil_resis

contains

  recursive function my_gamma(a) result(g)

   
    real(r_2), intent(in) :: a 
    real(r_2) :: g 

    real(r_2), parameter :: pi = 3.14159265358979324
    integer, parameter :: cg = 7

    ! these precomputed values are taken by the sample code in Wikipedia,
    ! and the sample itself takes them from the GNU Scientific Library
    real(r_2), dimension(0:8), parameter :: p = &
         (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
         771.32342877765313, -176.61502916214059, 12.507343278686905, &
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)

    real(r_2) :: t, w, x 
    integer :: i 

    x = a

    if ( x < 0.5 ) then 
       g = pi / ( sin(pi*x) * my_gamma(1.0-x) )
    else 
       x = x - 1.0
       t = p(0) 
       do i=1, cg+2 
          t = t + p(i-1)/(x+real(i,r_2))
       end do
       w = x + real(cg,r_2) + 0.5
       g = sqrt(2.0*pi) * w**(x+0.5) * exp(-w) * t
    end if
  end function my_gamma

SUBROUTINE or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough)

   TYPE (air_type), INTENT(INOUT)       :: air
   TYPE (met_type), INTENT(INOUT)       :: met
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type), INTENT(INOUT) :: veg
   TYPE (roughness_type), INTENT(INOUT) :: rough



   REAL(r_2), DIMENSION(mp) :: sublayer_dz, eddy_shape,eddy_mod,soil_moisture_mod, &
                          soil_moisture_mod_sat, wb_liq, &
                          pore_size,pore_radius, rel_s,hk_zero,hk_zero_sat,time_scale  !note pore_size in m

   REAL(r_2), DIMENSION(mp) :: litter_dz

   INTEGER, DIMENSION(mp) :: int_eddy_shape

   !

   integer :: i,j,k 

   if (cable_user%litter) then
      litter_dz(:) = veg%clitt*0.003
   else
      litter_dz(:) = 0.0
   endif

   pore_radius(:) = 0.148  / (1000.0*9.81*abs(soil%sucs_vec(:,1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle
   pore_size(:) = pore_radius(:)*sqrt(pi)

      !scale ustar according to the exponential wind profile, assuming we are a mm from the surface
      eddy_shape = 0.3*met%ua/ max(1.0e-4,canopy%us*exp(-rough%coexp*(1.0-canopy%sublayer_dz/max(1e-2,rough%hruff))))
      int_eddy_shape = floor(eddy_shape)
      eddy_mod(:) = 0.0
      do i=1,mp   
         eddy_mod(i) = 2.2*sqrt(112.0*pi) / (2.0**(eddy_shape(i)+1.0) * sqrt(eddy_shape(i)+1.0))

         if (int_eddy_shape(i) .gt. 0) then
            eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
            do k=1,int_eddy_shape(i)
               eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
            end do
         end if
      end do
      canopy%sublayer_dz = max(eddy_mod(:) * air%visc / max(1.0e-4,canopy%us*&
                           exp(-rough%coexp*(1.0-canopy%sublayer_dz/max(1e-2,rough%hruff)))),1e-7) 


   wb_liq(:) = real(max(0.0001,min(pi/4.0, &
                (ssnow%wb(:,1)-ssnow%wbice(:,1) - ssnow%satfrac(:)*soil%ssat_vec(:,1)) / &
                max((1._r_2 - ssnow%satfrac(:)),1e-5) ) ) )

   rel_s = real( max(wb_liq(:)-soil%watr(:,1),0._r_2)/(soil%ssat_vec(:,1)-soil%watr(:,1)) )
   hk_zero = max(0.001*soil%hyds_vec(:,1)*(min(max(rel_s,0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(:,1)+3._r_2) ),1e-8)
   hk_zero_sat = max(0.001*soil%hyds_vec(:,1),1e-8)

   soil_moisture_mod(:)     = 1.0/pi/sqrt(wb_liq)* ( sqrt(pi/(4.0*wb_liq))-1.0)
   soil_moisture_mod_sat(:) = 1.0/pi/sqrt(soil%ssat_vec(:,1))* ( sqrt(pi/(4.0*soil%ssat_vec(:,1)))-1.0)

   where(ssnow%isflag(:) .ne. 0)
      soil_moisture_mod = 0.
      soil_moisture_mod_sat = 0.
   elsewhere
      canopy%sublayer_dz = canopy%sublayer_dz + litter_dz
   endwhere


   where(canopy%sublayer_dz .ge. 1.0e-7) 
      ssnow%rtevap_unsat(:) = min(rtevap_max, &
                               rough%z0soil/canopy%sublayer_dz * (lm/ (4.0*hk_zero) +&
                               (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod) / Dff))
      ssnow%rtevap_sat(:)  = min(rtevap_max, &
                               rough%z0soil/canopy%sublayer_dz * (lm/ (4.0*hk_zero_sat) + &
                              (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod_sat) / Dff))

      ssnow%rt_qh_sublayer = canopy%sublayer_dz / litter_thermal_diff

   elsewhere
      ssnow%rtevap_unsat(:) = min(rtevap_max, &
                           lm/ (4.0*hk_zero) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod) / Dff)
      ssnow%rtevap_sat(:)  = min(rtevap_max, &
                         lm/ (4.0*hk_zero_sat) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod_sat) / Dff)

      ssnow%rt_qh_sublayer = 0.0
   endwhere


   !no additional evap resistane over lakes
   where(veg%iveg .eq. 16) 
      ssnow%rtevap_sat = 0.0
      ssnow%rtevap_unsat = 0.0
      ssnow%rt_qh_sublayer = 0.0
   endwhere


END SUBROUTINE or_soil_evap_resistance


SUBROUTINE update_or_soil_resis(ssnow,canopy,veg,dq,dqu)

   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy
   TYPE (veg_parameter_type), INTENT(INOUT) :: veg
   REAL, DIMENSION(mp), INTENT(IN) :: dq,&
                                      dqu

   INTEGER :: i


   do i=1,mp

      if (veg%iveg(i) .lt. 16 .and. ssnow%snowd(i) .lt. 1e-7) THEN

         if (dq(i) .le. 0.0) THEN
            ssnow%rtevap_sat(i) = min(rtevap_max,canopy%sublayer_dz(i)/Dff)
         end if

         if (dqu(i) .le. 0.0) THEN
            ssnow%rtevap_unsat(i) = min(rtevap_max,canopy%sublayer_dz(i)/Dff)
         end if

      end if
   end do






END SUBROUTINE update_or_soil_resis

END MODULE cable_psm

