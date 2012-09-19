!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Calculates surface exchange fluxes through the solution of surface 
!          energy balance and its interaction with plant physiology. Specific 
!        representation of the transport of scalars within a canopy is included.
!
! Called from: cbm
!
! Contact: Yingping.Wang@csiro.au and Eva.Kowalczyk@csiro.au
!
! History: Revision of canopy temperature calculation (relative to v1.4b) 
!          Reorganisation of code (dryLeaf, wetLeaf, photosynthesis subroutines
!          taken out of define_canopy)
!
!
! ==============================================================================

MODULE cable_canopy_module
   
   USE cable_data_module, ONLY : icanopy_type, point2constants 
   
   IMPLICIT NONE
   
   PUBLIC define_canopy
   PRIVATE
   
   TYPE( icanopy_type ) :: C
  
     
CONTAINS
 

SUBROUTINE define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy)
   USE cable_def_types_mod
   USE cable_radiation_module
   USE cable_air_module
   USE cable_common_module   

   TYPE (balances_type), INTENT(INOUT)  :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type), INTENT(INOUT)       :: air
   TYPE (met_type), INTENT(INOUT)       :: met
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy

   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   
   REAL, INTENT(IN)               :: dels ! integration time setp (s)
   
   INTEGER  ::                                                                 &
      iter,  & ! iteration #
      iterplus !

   REAL, DIMENSION(mp) ::                                                      &
      rt0,           & ! turbulent resistance
      ortsoil,       & ! turb. resist. prev t-step
      rt1usc,        & ! eq. 3.53, SCAM manual, 1997
      tstar,         & ! 
      zscrn,         & !
      qstar,         & !
      rsts,          & !
      qsurf,         & !
      qtgnet,        & !
      tss4,          & ! soil/snow temperature**4
      qstvair,       & ! sat spec humidity at leaf temperature
      xx,            & ! delta-type func 4 sparse canopy limit, p20 SCAM manual
      r_sc,          & !
      zscl,          & !
      pwet,          & !   
      dq,            & ! sat sp
      xx1,           & !
      sum_rad_rniso, & ! 
      sum_rad_gradis   ! 
   
   ! temporary buffers to simplify equations
   REAL, DIMENSION(mp) ::                                                      &
      ftemp,z_eff,psim_arg, psim_1, psim_2, rlower_limit,                      &
      term1, term2, term3, term5 

   REAL, DIMENSION(:), POINTER ::                                              & 
      cansat,        & ! max canopy intercept. (mm)
      dsx,           & ! leaf surface vpd
      fwsoil,        & ! soil water modifier of stom. cond
      tlfx,          & ! leaf temp prev. iter (K)
      tlfy             ! leaf temp (K)

   REAL(r_2), DIMENSION(mp) ::                                                 &
     gbvtop                   ! bnd layer cond. top leaf

   REAL(r_2), DIMENSION(:), POINTER ::                                         &
      ecy,           & ! lat heat fl dry big leaf
      hcy,           & ! veg. sens heat
      rny,           & ! net rad
      ghwet             ! cond for heat for a wet canopy
   
   REAL(r_2), DIMENSION(:,:), POINTER ::                                       &
      gbhu,          & ! forcedConvectionBndryLayerCond
      gbhf,          & ! freeConvectionBndryLayerCond
      csx              ! leaf surface CO2 concentration

   REAL  :: rt_min

   INTEGER :: j
   
   INTEGER, SAVE :: call_number =0
   
   ! END header
   
   call_number = call_number + 1
           
   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

   ! ACCESS version has this statement but elsewhere?
   IF( .NOT. cable_runtime%um)                                                 &
      canopy%cansto =  canopy%oldcansto

   ALLOCATE( cansat(mp), gbhu(mp,mf))
   ALLOCATE( dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp) )
   ALLOCATE( ecy(mp), hcy(mp), rny(mp))
   ALLOCATE( gbhf(mp,mf), csx(mp,mf))
   ALLOCATE( ghwet(mp))

   ! BATS-type canopy saturation proportional to LAI:
   cansat = veg%canst1 * canopy%vlaiw

   !---compute surface wetness factor, update cansto, through
   CALL surf_wetness_fact( cansat, canopy, ssnow,veg,met, soil, dels )

   canopy%fevw_pot = 0.0
   canopy%gswx = 1e-3     ! default stomatal conuctance 
   gbhf = 1e-3     ! default free convection boundary layer conductance
   gbhu = 1e-3     ! default forced convection boundary layer conductance
   ssnow%evapfbl = 0.0

   ! Initialise in-canopy temperatures and humidity:
   csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
   met%tvair = met%tk
   met%qvair = met%qv
   canopy%tv = met%tvair

   CALL define_air (met, air)
   
   CALL qsatfjh(qstvair,met%tvair-C%tfrz,met%pmb)

   met%dva = (qstvair - met%qvair) *  C%rmair/C%rmh2o * met%pmb * 100.0
   dsx = met%dva     ! init. leaf surface vpd
   
   tlfx = met%tk  ! initialise leaf temp iteration memory variable (K)
   tlfy = met%tk  ! initialise current leaf temp (K)
   
   ortsoil = ssnow%rtsoil
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)
   tss4 = ssnow%tss**4
   canopy%fes = 0.
   canopy%fess = 0.
   canopy%fesp = 0.
   ssnow%potev = 0.
   canopy%fevw_pot = 0.

   CALL radiation( ssnow, veg, air, met, rad, canopy )

   canopy%zetar(:,1) = C%ZETA0 ! stability correction terms
   canopy%zetar(:,2) = C%ZETPOS + 1 


   DO iter = 1, NITER

      ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
      ! resistances rt0, rt1 (elements of dispersion matrix):
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
      CALL comp_friction_vel()
      
      ! Turbulent aerodynamic resistance from roughness sublayer depth 
      ! to reference height, x=1 if zref+disp>zruffs, 
      ! 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
      xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)
      
      ! correction  by Ian Harman to the 2nd psis term
      rt1usc = xx * ( LOG( rough%zref_tq/MAX( rough%zruffs-rough%disp,         &
                                              rough%z0soilsn ) )               &
               - psis( canopy%zetar(:,iter) )                                  &
               + psis( canopy%zetar(:,iter) * ( MAX( rough%zruffs-rough%disp,  &
                                                rough%z0soilsn ) )            &
                       / rough%zref_tq ) ) / C%VONK
      
      rt_min = 5.      
      rt0 = max(rt_min,rough%rt0us / canopy%us)
      
      ! Aerodynamic resistance (sum 3 height integrals)/us
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
      rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
      
      DO j=1,mp
     
         IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
            ssnow%rtsoil(j) = rt0(j)
         ELSE
            ssnow%rtsoil(j) = rt0(j) + rough%rt1(j)
         ENDif 
     
      ENDDO 
      
      ssnow%rtsoil = max(rt_min,ssnow%rtsoil)   
      
      DO j=1,mp
      
         IF( ssnow%rtsoil(j) > 2.*ortsoil(j) .OR.                              &
             ssnow%rtsoil(j) < 0.5*ortsoil(j) ) THEN
              
            ssnow%rtsoil(j) = MAX(rt_min,0.5*(ssnow%rtsoil(j) + ortsoil(j)))
        
        ENDIF    
      
      ENDDO 
   
      ! Vegetation boundary-layer conductance (mol/m2/s)
      ! C%prandt = kinematic viscosity/molecular diffusivity
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
      DO j=1,mp

         IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
            gbvtop(j) = air%cmolar(j)*C%APOL * air%visc(j) / C%prandt /        &
                        veg%dleaf(j) * (canopy%us(j) / MAX(rough%usuh(j),1.e-6)&
                        * veg%dleaf(j) / air%visc(j) )**0.5                    &
                        * C%prandt**(1.0/3.0) / veg%shelrb(j)
            gbvtop(j) = MAX (0.05,gbvtop(j) )      ! for testing (BP aug2010)
            
            ! Forced convection boundary layer conductance                     
            ! (see Wang & Leuning 1998, AFM):
            gbhu(j,1) = gbvtop(j)*(1.0-EXP(-canopy%vlaiw(j)                    &
                        *(0.5*rough%coexp(j)+rad%extkb(j) ))) /                &
                        (rad%extkb(j)+0.5*rough%coexp(j))
            
            gbhu(j,2) = (2.0/rough%coexp(j))*gbvtop(j)*  &
                        (1.0-EXP(-0.5*rough%coexp(j)*canopy%vlaiw(j)))         &
                        - gbhu(j,1)
         ENDIF 
      
      ENDDO 
     
      rny = SUM(rad%rniso,2) ! init current estimate net rad
      hcy = 0.0              ! init current estimate lat heat
      ecy = rny - hcy        ! init current estimate lat heat

      sum_rad_rniso = SUM(rad%rniso,2)

      CALL dryLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, soil, ssnow, dsx,                             &
                    fwsoil, tlfx, tlfy, ecy, hcy,                              &
                    rny, gbhu, gbhf, csx, cansat,                              &
                    ghwet,  iter )
     
      CALL wetLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, cansat, tlfy,                                 &
                    gbhu, gbhf, ghwet )

     
      ! Calculate latent heat from vegetation:
      ! Calculate sensible heat from vegetation:
      ! Calculate net rad absorbed by canopy:
      canopy%fev = REAL(canopy%fevc + canopy%fevw)
      ftemp = (1.0 - canopy%fwet) *  REAL(hcy) + canopy%fhvw
      canopy%fhv = real(ftemp) 
      ftemp= (1.0-canopy%fwet)*REAL(rny)+canopy%fevw+canopy%fhvw
      canopy%fnv = real(ftemp)

      ! canopy rad. temperature calc from long-wave rad. balance
      sum_rad_gradis = SUM(rad%gradis,2)

      DO j=1,mp

         IF ( canopy%vlaiw(j) > C%LAI_THRESH .AND.                             &
              rough%hruff(j) > rough%z0soilsn(j) ) THEN

            rad%lwabv(j) = C%CAPP * C%rmair * ( tlfy(j) - met%tk(j) ) *        &
                           sum_rad_gradis(j) 

            canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
                           * C%SBOLTZ*C%EMLEAF)+met%tk(j)**4)**0.25
         
         ELSE! sparse canopy
         
           canopy%tv(j) = met%tk(j)
         
         ENDIF
          
      ENDDO 
     

      ! Calculate net rad to soil:
      canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*canopy%tv**4 - C%EMSOIL*C%SBOLTZ* tss4


      ! Saturation specific humidity at soil/snow surface temperature:
      call qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

      IF(cable_user%ssnow_POTEV== "P-M") THEN
         
         !--- uses %ga from previous timestep    
         ssnow%potev =  Penman_Monteith(canopy%ga) 
      
      ELSE !by default assumes Humidity Deficit Method
      
         dq = ssnow%qstss - met%qv
         ssnow%potev =  Humidity_deficit_method(dq,ssnow%qstss )
          
      ENDIF

      ! Soil latent heat:
      CALL latent_heat_flux()

      ! Calculate soil sensible heat:
      canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tk) /ssnow%rtsoil


      CALL within_canopy( gbhu, gbhf )

      ! Saturation specific humidity at soil/snow surface temperature:
      call qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

      IF(cable_user%ssnow_POTEV== "P-M") THEN
         
         !--- uses %ga from previous timestep    
         ssnow%potev =  Penman_Monteith(canopy%ga) 
      
      ELSE !by default assumes Humidity Deficit Method
      
         dq = ssnow%qstss - met%qvair
         ssnow%potev =  Humidity_deficit_method(dq,ssnow%qstss )
          
      ENDIF

         
      ! Soil latent heat:
      CALL latent_heat_flux()

      ! Soil sensible heat:
      canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
      canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssnow%cls
      
      ! Set total latent heat:
      canopy%fe = canopy%fev + canopy%fes
      
      ! Set total sensible heat:
      canopy%fh = canopy%fhv + canopy%fhs
      
      !---diagnostic purposes
      DO j=1,mp
      
         IF (ssnow%potev(j) .GE. 0.) THEN
            ssnow%potev(j) = max(0.00001,ssnow%potev(j))
         ELSE
            ssnow%potev(j) = min(-0.0002,ssnow%potev(j))
         ENDIF
         
         IF (canopy%fevw_pot(j) .ge. 0.) then
            canopy%fevw_pot(j) = max(0.000001,canopy%fevw_pot(j))
         ELSE
            canopy%fevw_pot(j) = min(-0.002,canopy%fevw_pot(j))
         ENDIF

      ENDDO 


      canopy%rnet = canopy%fnv + canopy%fns  
      canopy%epot = ((1.-rad%transd)*canopy%fevw_pot +                         &
                    rad%transd*ssnow%potev) * dels/air%rlam  

      ! convert to mm/day
      rlower_limit = canopy%epot * air%rlam / dels
      
      canopy%wetfac_cs = max(0., min(1.0,canopy%fe / rlower_limit ))
      
      DO j=1,mp

         IF ( canopy%wetfac_cs(j) .LE. 0. )                                    &
            canopy%wetfac_cs(j) = MAX( 0., MIN( 1.,                            &
                                  MAX( canopy%fev(j)/canopy%fevw_pot(j),       &
                                  canopy%fes(j)/ssnow%potev(j) ) ) )
      
      ENDDO 

      CALL update_zetar()

   END DO           ! do iter = 1, NITER


   canopy%cduv = canopy%us * canopy%us / (max(met%ua,C%UMIN))**2

   !---diagnostic purposes
   canopy%gswx_T = rad%fvlai(:,1)/MAX( C%LAI_THRESH, canopy%vlaiw(:) )         & 
                   * canopy%gswx(:,1) + rad%fvlai(:,2) / MAX(C%LAI_THRESH,     &
                   canopy%vlaiw(:))*canopy%gswx(:,2)

   canopy%gswx_T = max(1.e-05,canopy%gswx_T )
          
   canopy%cdtq = canopy%cduv *( LOG( rough%zref_uv / rough%z0m) -              &
                 psim( canopy%zetar(:,NITER) * rough%zref_uv/rough%zref_tq )   &
                 ) / ( LOG( rough%zref_uv /(0.1*rough%z0m) )                   &
                 - psis( canopy%zetar(:,NITER)) )

   ! Calculate screen temperature: 1) original method from SCAM
   ! screen temp., windspeed and relative humidity at 1.5m
   ! screen temp., windspeed and relative humidity at 2.0m
    tstar = - canopy%fh / ( air%rho*C%CAPP*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
    zscrn = MAX(rough%z0m,2.0-rough%disp)
    ftemp = ( LOG(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) +       &
            psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /C%VONK

    ! Calculate screen temperature:
    canopy%tscrn = met%tk - C%tfrz - tstar * ftemp

   ! Calculate radiative/skin temperature; 
   ! at this stage old soil temperature is used
   ! calculation of screen temepratures for LAI > 0.1 . Method by Ian Harman
    
   term1=0.
   term2=0.
   term5=0.
   term3 = 0. ! Work around for Intel compiler problem with nested whres
   r_sc = 0.
   zscl = MAX(rough%z0soilsn,2.0)

   ! assume screen temp of bareground if all these conditions are not met
   DO j=1,mp
      
      IF ( canopy%vlaiw(j) > C%LAI_THRESH .and. rough%hruff(j) > 0.01) THEN
      
         IF ( rough%disp(j)  > 0.0 ) then
     
            term1(j) = EXP(2*C%CSW*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
            term2(j) = EXP(2*C%CSW*canopy%rghlai(j) *                          &
                       (1-rough%disp(j)/rough%hruff(j)))
            term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)
         
         ENDIF
        
         term3(j) = C%A33**2*C%CTL*2*C%CSW*canopy%rghlai(j)

         IF( zscl(j) < rough%disp(j) ) THEN

            r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
                      ( EXP(2*C%CSW*canopy%rghlai(j)) - term1(j) ) / term3(j)

         ELSEIF( rough%disp(j) <= zscl(j) .AND.                                &
                 zscl(j) < rough%hruff(j) ) THEN
                 
            r_sc(j) = rough%rt0us(j) + term5(j) * ( term2(j) - term1(j) ) /    &
                      term3(j)

         ELSEIF( rough%hruff(j) <= zscl(j) .AND.                               &
                 zscl(j) <  rough%zruffs(j) ) THEN

            r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + term5(j) *            &
                      ( zscl(j) - rough%hruff(j) ) /                           &
                      ( C%A33**2 * C%CTL * rough%hruff(j) )


         ELSEIF( zscl(j) >= rough%zruffs(j) ) THEN
            
            r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) +     &
                      ( LOG( (zscl(j) - rough%disp(j)) /                       &
                      MAX( rough%zruffs(j)-rough%disp(j),                      &
                      rough%z0soilsn(j) ) ) - psis1( (zscl(j)-rough%disp(j))   &
                      / (rough%zref_tq(j)/canopy%zetar(j,iterplus) ) )         &
                      + psis1( (rough%zruffs(j) - rough%disp(j) )              &
                      / (rough%zref_tq(j)/canopy%zetar(j,iterplus ) ) ) )      &
                      / C%VONK

         ENDIF

        canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *          &
                          MIN(1.,r_sc(j) / MAX( 1.,                            &
                          rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)   &
                          + rt1usc(j))) - C%tfrz 
      ENDIF  

   ENDDO  
  
   CALL qsatfjh(rsts,canopy%tscrn,met%pmb)
     
   qtgnet = rsts * ssnow%wetfac - met%qv
   
   DO j=1,mp
      
      IF (qtgnet(j) .GT. 0. ) THEN
         qsurf(j) = rsts(j) * ssnow%wetfac(j)
      ELSE
         qsurf(j) = 0.1*rsts(j)*ssnow%wetfac(j) + 0.9*met%qv(j)
      ENDIF

      canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

      IF( canopy%vlaiw(j) >C%LAI_THRESH .and. rough%hruff(j) > 0.01)           &

            canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) * MIN( 1.,     &
                              r_sc(j) / MAX( 1., rough%rt0us(j) +              &
                              rough%rt1usa(j) + rough%rt1usb(j) + rt1usc(j) ) )

   ENDDO 


   ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
   canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
        dels * 1.0e3 / (C%RHOW*air%rlam)

   ! Add dewfall to canopy water storage:
   canopy%cansto = canopy%cansto + canopy%dewmm
   
   ! Modify canopy water storage for evaporation:
   canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw))*dels &
     *1.0e3/(C%RHOW*air%rlam), 0.0)

   ! Calculate canopy water storage excess:
   canopy%spill=max(0.0, canopy%cansto-cansat)

   ! Move excess canopy water to throughfall:
   canopy%through = canopy%through + canopy%spill
   
   ! Initialise 'throughfall to soil' as 'throughfall from canopy'; 
   ! snow may absorb
   canopy%precis = max(0.,canopy%through)
   
   ! Update canopy storage term:
   canopy%cansto=canopy%cansto - canopy%spill
   
   ! Calculate the total change in canopy water store (mm/dels):
   canopy%delwc = canopy%cansto-canopy%oldcansto
   
   ! calculate dgdtg, derivative of ghflux 3 instances
   ! d(canopy%fns)/d(ssnow%tgg)
   ! d(canopy%fhs)/d(ssnow%tgg)
   ! d(canopy%fes)/d(dq)
   ssnow%dfn_dtg = (-1.)*4.*C%EMSOIL*C%SBOLTZ*tss4/ssnow%tss  
   ssnow%dfh_dtg = air%rho*C%CAPP/ssnow%rtsoil      
   ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam/ssnow%rtsoil  
  
   ssnow%ddq_dtg = (C%rmh2o/C%rmair) /met%pmb * C%TETENA*C%TETENB * C%TETENC   &
                   / ( ( C%TETENC + ssnow%tss-C%tfrz )**2 )*EXP( C%TETENB *       &
                   ( ssnow%tss-C%tfrz ) / ( C%TETENC + ssnow%tss-C%tfrz ) )
   canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%cls*ssnow%dfe_ddq *    &
                  ssnow%ddq_dtg

   bal%drybal = REAL(ecy+hcy) - SUM(rad%rniso,2)                               &
                + C%CAPP*C%rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009

   bal%wetbal = canopy%fevw + canopy%fhvw - SUM(rad%rniso,2) * canopy%fwet      &
                + C%CAPP*C%rmair * (tlfy-met%tk) * SUM(rad%gradis,2) *          &
                canopy%fwet  ! YP nov2009

   DEALLOCATE(cansat,gbhu)
   DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
   DEALLOCATE(ecy, hcy, rny)
   DEALLOCATE(gbhf, csx)
   DEALLOCATE(ghwet)

CONTAINS

! ------------------------------------------------------------------------------

SUBROUTINE comp_friction_vel()
   USE cable_def_types_mod, only : mp
   REAL, DIMENSION(mp)  :: lower_limit, rescale

   psim_1 = psim(canopy%zetar(:,iter)) 
   
      rescale = C%VONK * MAX(met%ua,C%UMIN)
      z_eff = rough%zref_uv / rough%z0m
   
   psim_arg = canopy%zetar(:,iter) / z_eff 
   !---fix for compiler limitation. bitwise reproducable whilst we  
   !---we know it to 11th decimal. psim_arg typically of a few 
   !psim_arg = nint(psim_arg * 1.e11)*1.e-11

   psim_2 = psim( psim_arg )
          
   lower_limit = rescale / ( LOG(z_eff) - psim_1 + psim_2 )

   canopy%us = MAX(1.e-6, lower_limit )

END SUBROUTINE comp_friction_vel

! ------------------------------------------------------------------------------

FUNCTION Penman_Monteith( ground_H_flux ) RESULT(ssnowpotev)
   USE cable_def_types_mod, only : mp
   REAL, INTENT(IN), DIMENSION(mp)  :: ground_H_flux
   REAL, DIMENSION(MP)  ::                                                     &
      ssnowpotev,      & ! returned result of function    
      sss,             & ! var for Penman-Monteith soil evap
      cc1,             & ! var for Penman-Monteith soil evap
      cc2,             & ! var for Penman-Monteith soil evap
      qsatfvar           !
   INTEGER :: j
   
   ! Penman-Monteith formula
   sss=air%dsatdk
   cc1=sss/(sss+air%psyc )
   cc2=air%psyc /(sss+air%psyc )
   
   CALL qsatfjh(qsatfvar,met%tvair-C%tfrz,met%pmb)

   ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
   cc2 * air%rho * air%rlam*(qsatfvar  - met%qvair)/ssnow%rtsoil
 
END FUNCTION Penman_Monteith


! ------------------------------------------------------------------------------
! method alternative to P-M formula above
FUNCTION humidity_deficit_method(dq,qstss ) RESULT(ssnowpotev)

   USE cable_def_types_mod, only : mp
   
   REAL, DIMENSION(mp) ::                                                      &
      ssnowpotev,    & ! 
      dq,            & ! sat spec hum diff.
      qstss             !dummy var for compilation
       
   INTEGER :: j
   
   DO j=1,mp
      !if(ssnow%snowd(j) > 1.0) dq(j) = max( -0.1e-3, dq(j))
      IF( ssnow%snowd(j)>1.0 .OR. ssnow%tgg(j,1).EQ.C%tfrz)                      &
         dq(j) = max( -0.1e-3, dq(j))
   ENDDO 
   
   ssnowpotev =air%rho * air%rlam * dq /ssnow%rtsoil
   
END FUNCTION Humidity_deficit_method

! ------------------------------------------------------------------------------
 
SUBROUTINE Latent_heat_flux() 

   USE cable_common_module
   USE cable_def_types_mod, only : mp

   REAL, DIMENSION(mp) ::                                                      &
      frescale,  flower_limit, fupper_limit

   INTEGER :: j
   
   ! Soil latent heat:
   canopy%fess= ssnow%wetfac * ssnow%potev
   WHERE (ssnow%potev < 0. ) canopy%fess = ssnow%potev
   
   ! Reduce soil evap due to presence of puddle
   pwet = max(0.,min(0.2,ssnow%pudsto/max(1.,ssnow%pudsmx)))
   canopy%fess = canopy%fess * (1.-pwet)

   frescale = soil%zse(1) * 1000. * air%rlam / dels         

   DO j=1,mp
      
      IF(ssnow%snowd(j) < 0.1 .AND. canopy%fess(j) .GT. 0. ) THEN

         flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)/2.0
         fupper_limit(j) = MAX( 0._r_2,                                        &
                           flower_limit(j) * frescale(j)                       &
                           - ssnow%evapfbl(j,1)*air%rlam(j)/dels)

         canopy%fess(j) = MIN(canopy%fess(j), fupper_limit(j))
         
         fupper_limit(j) = REAL(ssnow%wb(j,1)-ssnow%wbice(j,1)) * frescale(j)

         canopy%fess(j) = min(canopy%fess(j), fupper_limit(j))

      END IF

      ssnow%cls(j)=1.
      
      IF (ssnow%snowd(j) >= 0.1 .and. ssnow%potev(j) > 0.) THEN

         ssnow%cls(j) = 1.1335
         canopy%fess(j) = MIN( (ssnow%wetfac(j)*ssnow%potev(j)),               &
                          ssnow%snowd(j)/dels*air%rlam(j)*ssnow%cls(j))
      
      ENDIF

   ENDDO 
   
   ! Evaporation form soil puddle
   canopy%fesp = min(ssnow%pudsto/dels*air%rlam,max(pwet*ssnow%potev,0.))
   canopy%fes = canopy%fess + canopy%fesp

END SUBROUTINE latent_heat_flux

! -----------------------------------------------------------------------------

SUBROUTINE within_canopy( gbhu, gbhf )

   USE cable_def_types_mod, only : mp, r_2

   REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
      gbhu,    &  ! forcedConvectionBndryLayerCond
      gbhf        ! freeConvectionBndryLayerCond
      
   REAL, DIMENSION(mp) ::                                                      &
      rrsw,             & ! recipr. stomatal resistance for water
      rrbw,             & ! recipr. leaf boundary layer resistance for water
      dmah,             & ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
      dmbh,             & ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
      dmch,             & ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
      dmae,             & ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
      dmbe,             & ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
      dmce                ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132

   REAL  :: lower_limit, upper_limit
 
   INTEGER :: j
   
   rrbw = sum(gbhu+gbhf,2)/air%cmolar  ! MJT 
   
   ! leaf stomatal resistance for water
   rrsw = sum(canopy%gswx,2)/air%cmolar ! MJT
   
   DO j=1,mp
   
      IF(veg%meth(j) > 0 .AND. canopy%vlaiw(j) > C%LAI_THRESH .AND.              &
         rough%hruff(j) > rough%z0soilsn(j) ) THEN

         !   use the dispersion matrix (DM) to find the air temperature 
         !   and specific humidity 
         !   (Raupach, Finkele and Zhang 1997, pp 17)
         ! leaf boundary layer resistance for water
         ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmah(j) = (rt0(j)+rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j))  &
                   + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
         
         ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbh(j) = (-air%rlam(j)/C%CAPP)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
         
         ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                   (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*C%CAPP)

         ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmae(j) = (-air%epsi(j)*C%CAPP/air%rlam(j))*(rt0(j)*rough%rt1(j)) *   &
                   (rrbw(j)*rrsw(j))

         ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbe(j) = ( rt0(j) + ssnow%wetfac(j) * rough%rt1(j) ) *               &
                   ( (1.+air%epsi(j) ) * rrsw(j) + rrbw(j) ) +                 &
                   ( rt0(j) * rough%rt1(j) ) * ( rrbw(j) * rrsw(j) )

         ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                   (canopy%fev(j) + canopy%fes(j))/(air%rho(j)*air%rlam(j))
      
         ! Within canopy air temperature:
         met%tvair(j) = met%tk(j) + ( dmbe(j) * dmch(j) - dmbh(j) * dmce(j) )  &
                        / (dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
           
         !---set limits for comparisson
         lower_limit =  MIN( ssnow%tss(j), met%tk(j)) - 5.0
         upper_limit =  MAX( ssnow%tss(j), met%tk(j)) + 5.0
      
         !--- tvair within these limits
         met%tvair(j) = MAX(met%tvair(j) , lower_limit)
         met%tvair(j) = MIN(met%tvair(j) , upper_limit)
   
         ! recalculate using canopy within temperature
         met%qvair(j) = met%qv(j) + (dmah(j)*dmce(j)-dmae(j)*dmch(j)) /        &
                        ( dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
         met%qvair(j) = MAX(0.0,met%qvair(j))
      
         !---set limits for comparisson
         lower_limit =  MIN( ssnow%qstss(j), met%qv(j)) 
         upper_limit =  MAX( ssnow%qstss(j), met%qv(j))
      
            !--- qvair within these limits
         met%qvair(j) =  MAX(met%qvair(j),lower_limit)
         met%qvair(j) =  MIN(met%qvair(j), upper_limit)
      
         ! Saturated specific humidity in canopy:
         CALL qsatfjh2(qstvair(j),met%tvair(j)-C%tfrz,met%pmb(j))
         
         ! Saturated vapour pressure deficit in canopy:
         met%dva(j) = ( qstvair(j) - met%qvair(j) ) *  C%rmair/C%RMH2o         &
                      * met%pmb(j) * 100.
      ENDIF 

   ENDDO 
     
END SUBROUTINE within_canopy

! -----------------------------------------------------------------------------

SUBROUTINE update_zetar()
   
   INTEGER :: j
   
   ! monin-obukhov stability parameter zetar=zref/l
   ! recompute zetar for the next iteration, except on last iteration
   IF (iter < NITER) THEN ! dont compute zetar on the last iter

      iterplus = MAX(iter+1,2)
      canopy%zetar(:,iterplus) = -( C%VONK * C%GRAV * rough%zref_tq *              &
                                 ( canopy%fh + 0.07 * canopy%fe ) ) /          &
                                 ( air%rho * C%CAPP * met%tk * canopy%us**3 )

      ! case NITER=2: final zetar=C%ZETmul*zetar(2) (compute only when iter=1)
      IF (NITER == 2) THEN
    
         canopy%zetar(:,2) = C%ZETmul * canopy%zetar(:,2)
    
         DO j=1,mp
            IF ( (met%fsd(j,1)+met%fsd(j,2))  ==  0.0 ) &
               canopy%zetar(j,2) = 0.5 * canopy%zetar(j,2)
         ENDDO

      END IF

      ! constrain zeta to C%ZETPOS and C%ZETNEG (set in param0)
      
      ! zetar too +
      canopy%zetar(:,iterplus) = MIN(C%ZETPOS,canopy%zetar(:,iterplus))        
      
      ! zetar too -
      canopy%zetar(:,iterplus) = MAX(C%ZETNEG,canopy%zetar(:,iterplus))        
    
   END IF ! (iter < NITER)
      
END SUBROUTINE update_zetar

! -----------------------------------------------------------------------------

FUNCTION qsatf(j,tair,pmb) RESULT(r)
  ! MRR, 1987
  ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
  ! HUMIDITY (KG/KG) FROM TETEN FORMULA
   
   REAL, INTENT(IN) ::                                                         &
      tair,         & ! air temperature (C)
      pmb             ! pressure PMB (mb)
  
   INTEGER, INTENT(IN) :: j 
  
   REAL           :: r    ! result; sat sp humidity
  
   r = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb
END FUNCTION qsatf

! -----------------------------------------------------------------------------

SUBROUTINE qsatfjh(var,tair,pmb) 
   USE cable_def_types_mod, only : mp
   REAL, INTENT(IN), DIMENSION(mp) ::                                          &
      tair,                        & ! air temperature (C)
      pmb                            ! pressure PMB (mb)
  
   REAL, INTENT(OUT), DIMENSION(mp) ::                                         &
      var                            ! result; sat sp humidity
  
  INTEGER :: j
  
  DO j=1,mp

     var(j) = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair(j)/(C%TETENC+tair(j))))    &
              / pmb(j)
  ENDDO   

END SUBROUTINE qsatfjh

! -----------------------------------------------------------------------------

SUBROUTINE qsatfjh2(var,tair,pmb) 
   
   REAL, INTENT(IN) ::                                                         &
      tair,         & ! air temperature (C)
      pmb             ! pressure PMB (mb)
   
   REAL, INTENT(OUT) ::                                                        &
      var             ! result; sat sp humidity
      
      var = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb

END SUBROUTINE qsatfjh2

! -----------------------------------------------------------------------------

FUNCTION psim(zeta) RESULT(r)
   USE cable_def_types_mod, only : mp
   ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
   ! computes integrated stability function psim(z/l) (z/l=zeta)
   ! for momentum, using the businger-dyer form for unstable cases
   ! and the Beljaars and Holtslag (1991) form for stable cases.
   
   
   REAL, INTENT(IN), DIMENSION(mp) ::  zeta       !

   ! function result   
   REAL, DIMENSION(mp) :: r
   
   REAL, DIMENSION(mp) ::                                                      &
      x,       & !
      z,       & !
      stable,  & !
      unstable   !
      
   REAL, PARAMETER ::                                                          &
      gu = 16.0,  & !
      gs = 5.0
  
   REAL, PARAMETER ::                                                          &
      a = 1.0,       & !
      b = 0.667,     & !
      xc = 5.0,       & !
      d = 0.35         !

   z = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
   
   ! Beljaars and Holtslag (1991) for stable
   stable = -a*zeta - b*(zeta - xc/d)*exp( -d*zeta) - b*xc/d
   x      = (1.0 + gu*abs(zeta))**0.25
   unstable = ALOG((1.0+x*x)*(1.0+x)**2/8) - 2.0*atan(x) + C%PI_C*0.5
   r = z*stable + (1.0-z)*unstable
END FUNCTION psim

! -----------------------------------------------------------------------------

ELEMENTAL FUNCTION psis(zeta) RESULT(r)

   ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
   ! computes integrated stability function psis(z/l) (z/l=zeta)
   ! for scalars, using the businger-dyer form for unstable cases
   ! and the webb form for stable cases. see paulson (1970).

   REAL, INTENT(IN)     :: zeta
   
   REAL, PARAMETER      ::                                                     &
      gu = 16.0,        & !
      gs = 5.0,         & !
      a = 1.0,          & !
      b = 0.667,        & !
      c = 5.0,          & !
      d = 0.35
 
   REAL                 ::                                                     &
      r,                & !
      stzeta,           & !
      ustzeta,          & !
      z,                & !
      y,                & !
      stable,           & !
      unstable
   
   z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 

   ! Beljaars and Holtslag (1991) for stable
   stzeta = MAX(0.,zeta)
   stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
             b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.
   y      = (1.0 + gu*abs(zeta))**0.5
   unstable = 2.0 * alog((1+y)*0.5)
   r   = z*stable + (1.0-z)*unstable

END FUNCTION psis

! -----------------------------------------------------------------------------

FUNCTION psis1(zeta) RESULT(r)
   ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
   ! computes integrated stability function psis(z/l) (z/l=zeta)
   ! for scalars, using the businger-dyer form for unstable cases
   ! and the webb form for stable cases. see paulson (1970).
   REAL, INTENT(IN)     :: zeta
   
   REAL, PARAMETER      :: gu = 16.0
   REAL, PARAMETER      :: gs = 5.0
   REAL, PARAMETER      :: a = 1.0
   REAL, PARAMETER      :: b = 0.667
   REAL, PARAMETER      :: c = 5.0
   REAL, PARAMETER      :: d = 0.35
 
   REAL                 :: r
   REAL                 :: stable
   REAL                 :: unstable
   REAL                 :: stzeta
 
   REAL                 :: z
   REAL                 :: y
   !REAL                 :: stable
   !REAL                 :: unstable
 
   z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
   
   ! Beljaars and Holtslag (1991) for stable
   stzeta = max(0.,zeta)
   stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
             b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.
 
   y      = (1.0 + gu*abs(zeta))**0.5
   unstable = 2.0 * alog((1+y)*0.5)
   r   = z*stable + (1.0-z)*unstable

END FUNCTION psis1

! -----------------------------------------------------------------------------

ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
   REAL, INTENT(IN)     :: rpconst
   REAL, INTENT(IN)     :: rpcoef
   REAL, INTENT(IN)     :: tair
   REAL                 :: z
   z = rpconst * exp(rpcoef * tair)
END FUNCTION rplant

! -----------------------------------------------------------------------------

SUBROUTINE wetLeaf( dels, rad, rough, air, met, veg, canopy, cansat, tlfy,     &
                    gbhu, gbhf, ghwet )

   USE cable_def_types_mod

   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (canopy_type),    INTENT(INOUT) :: canopy

   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

   REAL,INTENT(IN), DIMENSION(:) ::                                            &
      tlfy,          & ! leaf temp (K) - assC%UMINg the temperature of 
                       ! wet leaf is equal that of dry leaf ="tlfy"
      cansat           ! max canopy intercept. (mm)

   REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
      gbhu,          & ! forcedConvectionBndryLayerCond
      gbhf             ! freeConvectionBndryLayerCond

   REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
      ghwet            ! cond for heat for a wet canopy

   REAL, INTENT(IN)     :: dels ! integration time step (s)

   ! local variables  
   REAL, DIMENSION(mp) ::                                                      &
     ccfevw,        & ! limitation term for
     gwwet,         & ! cond for water for a wet canopy
     ghrwet           ! wet canopy cond: heat & thermal rad
   
   !i sums, terms of convenience/readability
   REAL, DIMENSION(mp) ::                                                      &
     sum_gbh, sum_rad_rniso, sum_rad_gradis, xx1

   INTEGER :: j
   
   ! END header

   ghwet = 1.0e-3
   gwwet = 1.0e-3
   ghrwet= 1.0e-3
   canopy%fevw = 0.0
   canopy%fhvw = 0.0
   sum_gbh = SUM((gbhu+gbhf),2)
   sum_rad_rniso = SUM(rad%rniso,2)
   sum_rad_gradis = SUM(rad%gradis,2)

   DO j=1,mp

      IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN

         ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
         ! calculate total thermal resistance, rthv in s/m
         ghwet(j) = 2.0   * sum_gbh(j) 
         gwwet(j) = 1.075 * sum_gbh(j) 
         ghrwet(j) = sum_rad_gradis(j) + ghwet(j)
         
         ! Calculate fraction of canopy which is wet:
         canopy%fwet(j) = MAX( 0.0, MIN( 1.0,                                  &
                          0.8 * canopy%cansto(j) / MAX( cansat(j), 0.01 ) ) )
         
         ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
         ! to avoid excessive evaporation:
         ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                       2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )
   
         canopy%fevw(j) = MIN( canopy%fwet(j) * ( air%dsatdk(j) *              &
                         ( sum_rad_rniso(j)- C%CAPP*C%rmair*( met%tvair(j)     &
                         - met%tk(j) ) * sum_rad_gradis(j) )                   &
                         + C%CAPP * C%rmair * met%dva(j) * ghrwet(j) )         &
                         / ( air%dsatdk(j)+air%psyc(j)*ghrwet(j) / gwwet(j) )  &
                         , ccfevw(j) )

         canopy%fevw_pot(j) = ( air%dsatdk(j)* (sum_rad_rniso(j) -             &
                              C%CAPP * C%rmair * ( met%tvair(j) - met%tk(j) )  &
                              *sum_rad_gradis(j) )                             &
                              + C%CAPP * C%rmair * met%dva(j) * ghrwet(j))     &
                              / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )
          
         ! calculate sens heat from wet canopy:
         canopy%fhvw(j) = canopy%fwet(j) * ( sum_rad_rniso(j) -C%CAPP * C%rmair&
                          * ( tlfy(j) - met%tk(j) ) * sum_rad_gradis(j) )      &
                           - canopy%fevw(j)

          xx1(j) = canopy%fhvw(j)

         canopy%fhvw(j) = sum_rad_rniso(j) - canopy%fevc(j) - canopy%fevw(j)   &
                          - ( 1.0 - canopy%fwet(j) ) *  REAL( hcy(j) ) 
          
         canopy%fhvw(j) =  xx1(j)

      ENDIF
       
   ENDDO 
           
END SUBROUTINE wetLeaf

! -----------------------------------------------------------------------------

END SUBROUTINE define_canopy

! -----------------------------------------------------------------------------

SUBROUTINE Surf_wetness_fact( cansat, canopy, ssnow,veg, met, soil, dels )

   USE cable_common_module
   USE cable_def_types_mod

   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   TYPE (soil_snow_type), intent(inout):: ssnow
   TYPE (soil_parameter_type), intent(inout)   :: soil
   TYPE (canopy_type), INTENT(INOUT)   :: canopy
   TYPE (met_type), INTENT(INOUT)   :: met
   
   REAL, INTENT(IN) :: dels ! integration time setp (s)

   REAL,INTENT(IN), DIMENSION(:) :: cansat ! max canopy intercept. (mm)

   !local variables
   REAL, DIMENSION(mp)  :: lower_limit, upper_limit,ftemp
   
   INTEGER :: j

   ! Rainfall variable is limited so canopy interception is limited,
   ! used to stabilise latent fluxes.
   ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
   upper_limit = 4.0 * MIN(dels,1800.0) / (60.0 * 1440.0 ) 
   ftemp =MIN(met%precip-met%precip_sn, upper_limit )
   ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
   lower_limit = cansat - canopy%cansto
   upper_limit = max(lower_limit, 0.0) 
   canopy%wcint = MERGE( MIN( upper_limit, ftemp ), 0.0,                       &
                  ftemp > 0.0  .AND. met%tk > C%tfrz)  !EAK, 09/10

   ! Define canopy throughfall (100% of precip if temp < 0C, see above):
   canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn ,          &
                    MAX( 0.0, met%precip - met%precip_sn - canopy%wcint) ) 

   ! Add canopy interception to canopy storage term:
   canopy%cansto = canopy%cansto + canopy%wcint

   ! Calculate fraction of canopy which is wet:
   canopy%fwet   = MAX( 0.0, MIN( 0.9, 0.8 * canopy%cansto /                   &
                   MAX( cansat, 0.01 ) ) )

   ssnow%wetfac = MAX( 1.e-6, MIN( 1.0,                                        &
                  ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0 )                  &
                  / ( soil%sfc - soil%swilt/2.0 ) ) )
  
   DO j=1,mp
   
      IF( ssnow%wbice(j,1) > 0. )                                              &
         ssnow%wetfac(j) = ssnow%wetfac(j) * MAX( 0.5, 1. - MIN( 0.2,          &
                           ( ssnow%wbice(j,1) / ssnow%wb(j,1) )**2 ) )

      IF( ssnow%snowd(j) > 0.1) ssnow%wetfac(j) = 0.9
      
      IF ( veg%iveg(j) == 16 .and. met%tk(j) >= C%tfrz + 5. )                  &
         ssnow%wetfac(j) = 1.0 ! lakes: hard-wired number to be removed
      
      IF( veg%iveg(j) == 16 .and. met%tk(j) < C%tfrz + 5. )                    &
         ssnow%wetfac(j) = 0.7 ! lakes: hard-wired number to be removed

   ENDDO 
      

   ! owetfac introduced to reduce sharp changes in dry regions,
   ! especially in offline runs in which there may be discrepancies b/n
   ! timing of precip and temperature change (EAK apr2009)
   ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

END SUBROUTINE Surf_wetness_fact

! -----------------------------------------------------------------------------

SUBROUTINE dryLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, soil, ssnow, dsx,                             &
                    fwsoil, tlfx,  tlfy,  ecy, hcy,                            &
                    rny, gbhu, gbhf, csx,                                      &
                    cansat, ghwet, iter )

   USE cable_def_types_mod
   USE cable_common_module

   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow

   TYPE (veg_parameter_type),  INTENT(INOUT)   :: veg
   TYPE (soil_parameter_type), INTENT(inout)   :: soil
   
   REAL, INTENT(INOUT), DIMENSION(:) ::                                        &
      dsx,        & ! leaf surface vpd
      fwsoil,     & ! soil water modifier of stom. cond
      tlfx,       & ! leaf temp prev. iter (K)
      tlfy          ! leaf temp (K)
   
   REAL(R_2),INTENT(INOUT), DIMENSION(:) ::                                    &
      ecy,        & ! lat heat fl dry big leaf
      hcy,        & ! veg. sens heat
      rny         !& !

   REAL(R_2),INTENT(INOUT), DIMENSION(:,:) ::                                  &
      gbhu,       & ! forcedConvectionBndryLayerCond
      gbhf,       & ! freeConvectionBndryLayerCond
      csx           ! leaf surface CO2 concentration

   REAL,INTENT(IN), DIMENSION(:) :: cansat

   REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
      ghwet  ! cond for heat for a wet canopy

   INTEGER,INTENT(IN) :: iter
 
   REAL, INTENT(IN)     :: dels ! integration time step (s)
   
   !local variables
   REAL, PARAMETER  ::                                                         &
      co2cp3 = 0.0,  & ! CO2 compensation pt C3
      jtomol = 4.6e-6  ! Convert from J to Mol for light

   REAL, DIMENSION(mp) ::                                                      &
      conkct,        & ! Michaelis Menton const.
      conkot,        & ! Michaelis Menton const.
      cx1,           & ! "d_{3}" in Wang and Leuning,
      cx2,           & !     1998, appendix E
      tdiff,         & ! leaf air temp diff.
      tlfxx,         & ! leaf temp of current iteration (K)
      abs_deltlf,    & ! ABS(deltlf)
      deltlf,        & ! deltlfy of prev iter.
      deltlfy,       & ! del temp successive iter.
      gras,          & ! Grashof coeff
      evapfb,        & !
      sum_rad_rniso, & !
      sum_rad_gradis,& ! 
      gwwet,         & ! cond for water for a wet canopy
      ghrwet,        & ! wet canopy cond: heat & thermal rad
      sum_gbh,       & !
      ccfevw,        & ! limitation term for
                       ! wet canopy evaporation rate
      temp             !

   REAL(r_2), DIMENSION(mp)  ::                                                &
      ecx,        & ! lat. hflux big leaf
      ecx_t,      & ! lat. hflux big leaf
      hcx,        & ! sens heat fl big leaf prev iteration
      rnx           ! net rad prev timestep

   REAL, DIMENSION(mp,ms)  :: oldevapfbl
      
   REAL, DIMENSION(mp,mf)  ::                                                  &
      gw,         & ! cond for water for a dry canopy
      gh,         & ! cond for heat for a dry canopy
      ghr,        & ! dry canopy cond for heat & thermal rad
      anx,        & ! net photos. prev iteration
      an_y,       & ! net photosynthesis soln
      rdx,        & ! daytime leaf resp rate, prev iteration
      rdy,        & ! daytime leaf resp rate
      ejmax2,     & ! jmax of big leaf
      ejmxt3,     & ! jmax big leaf C3 plants
      vcmxt3,     & ! vcmax big leaf C3
      vcmxt4,     & ! vcmax big leaf C4
      vx3,        & ! carboxylation C3 plants
      vx4,        & ! carboxylation C4 plants
      xleuning,   & ! leuning stomatal coeff
      psycst,     & ! modified pych. constant
      frac42,     & ! 2D frac4
      temp2

   REAL, DIMENSION(:,:), POINTER :: gswmin ! min stomatal conductance
   
   REAL, DIMENSION(mp,2) ::  gsw_term, lower_limit2  ! local temp var 

   INTEGER :: i, j, k, kk  ! iteration count
   
   ! END header


   ALLOCATE( gswmin(mp,mf ))

   ! Soil water limitation on stomatal conductance:
   IF( iter ==1) THEN
   
      IF(cable_user%FWSOIL_SWITCH == 'standard') THEN
         CALL fwsoil_calc_std( fwsoil, soil, ssnow, veg) 
      ELSEIf (cable_user%FWSOIL_SWITCH == 'non-linear extrapolation') THEN
         !EAK, 09/10 - replace linear approx by polynomial fitting
         CALL fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg) 
      ELSEIF(cable_user%FWSOIL_SWITCH == 'Lai and Ktaul 2000') THEN
         CALL fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg) 
      ELSE
         STOP 'fwsoil_switch failed.'
      ENDIF

   ENDIF

   ! weight min stomatal conductance by C3 an C4 plant fractions
   frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants

   gsw_term = C%gsw03 * (1. - frac42) + C%gsw04 * frac42
   lower_limit2 = rad%scalex * (C%gsw03 * (1. - frac42) + C%gsw04 * frac42)
   gswmin = max(1.e-6,lower_limit2)
         

   gw = 1.0e-3 ! default values of conductance
   gh = 1.0e-3
   ghr= 1.0e-3
   rdx = 0.0
   anx = 0.0
   rnx = SUM(rad%rniso,2)
   abs_deltlf = 999.0


   gras = 1.0e-6
   an_y= 0.0
   hcx = 0.0              ! init sens heat iteration memory variable
   hcy = 0.0
   rdy = 0.0
   ecx = SUM(rad%rniso,2) ! init lat heat iteration memory variable
   tlfxx = tlfx
   psycst(:,:) = SPREAD(air%psyc,2,mf)
   canopy%fevc = 0.0
   ssnow%evapfbl = 0.0

   ghwet = 1.0e-3
   gwwet = 1.0e-3
   ghrwet= 1.0e-3
   canopy%fevw = 0.0
   canopy%fhvw = 0.0
   sum_gbh = SUM((gbhu+gbhf),2)
   sum_rad_rniso = SUM(rad%rniso,2)
   sum_rad_gradis = SUM(rad%gradis,2)

   DO kk=1,mp

      IF(canopy%vlaiw(kk) <= C%LAI_THRESH) THEN
         rnx(kk) = 0.0 ! intialise
         ecx(kk) = 0.0 ! intialise
         ecy(kk) = ecx(kk) ! store initial values
         abs_deltlf(kk)=0.0
         rny(kk) = rnx(kk) ! store initial values
         ! calculate total thermal resistance, rthv in s/m
      END IF

   ENDDO
   
   deltlfy = abs_deltlf
   k = 0

   !kdcorbin, 08/10 - doing all points all the time
   DO WHILE (k < C%MAXITER)
      k = k + 1

      DO i=1,mp
         
         IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1) THEN
         
            ghwet(i) = 2.0   * sum_gbh(i)
            gwwet(i) = 1.075 * sum_gbh(i)
            ghrwet(i) = sum_rad_gradis(i) + ghwet(i)
       
            ! Calculate fraction of canopy which is wet:
            canopy%fwet(i) = MAX( 0.0, MIN( 1.0, 0.8 * canopy%cansto(i)/ MAX(  &
                             cansat(i),0.01 ) ) )

            ! Calculate lat heat from wet canopy, may be neg.                  
            ! if dew on wet canopy to avoid excessive evaporation:
            ccfevw(i) = MIN(canopy%cansto(i) * air%rlam(i) / dels,             &
                        2.0 / (1440.0 / (dels/60.0)) * air%rlam(i) )
   
            ! Grashof number (Leuning et al, 1995) eq E4:
            gras(i) = MAX(1.0e-6,                                              &
                      1.595E8* ABS( tlfx(i)-met%tvair(i))* (veg%dleaf(i)**3.0) )

            ! See Appendix E in (Leuning et al, 1995):
            gbhf(i,1) = rad%fvlai(i,1) * air%cmolar(i) * 0.5*C%dheat           &
                       * ( gras(i)**0.25 ) / veg%dleaf(i)
            gbhf(i,2) = rad%fvlai(i,2) * air%cmolar(i) * 0.5 * C%dheat         &
                        * ( gras(i)**0.25 ) / veg%dleaf(i)
            gbhf(i,:) = MAX( 1.e-6, gbhf(i,:) )
      
            ! Conductance for heat:
            gh(i,:) = 2.0 * (gbhu(i,:) + gbhf(i,:))
      
            ! Conductance for heat and longwave radiation:
            ghr(i,:) = rad%gradis(i,:)+gh(i,:)
      
            ! Leuning 2002 (P C & E) equation for temperature response
            ! used for Vcmax for C3 plants:
            temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
            
            vcmxt3(i,1) = rad%scalex(i,1) * temp(i)
            vcmxt3(i,2) = rad%scalex(i,2) * temp(i)
    
            ! Temperature response Vcmax, C4 plants (Collatz et al 1989):
            temp(i) = xvcmxt4(tlfx(i)-C%tfrz) * veg%vcmax(i) * veg%frac4(i)
            vcmxt4(i,1) = rad%scalex(i,1) * temp(i)
            vcmxt4(i,2) = rad%scalex(i,2) * temp(i)
    
            ! Leuning 2002 (P C & E) equation for temperature response
            ! used for Jmax for C3 plants:
            temp(i) = xejmxt3(tlfx(i)) * veg%ejmax(i) * (1.0-veg%frac4(i))
            ejmxt3(i,1) = rad%scalex(i,1) * temp(i)
            ejmxt3(i,2) = rad%scalex(i,2) * temp(i)
            
            ! Difference between leaf temperature and reference temperature:
            tdiff(i) = tlfx(i) - C%TREFK
            
            ! Michaelis menten constant of Rubisco for CO2:
            conkct(i) = C%conkc0 * EXP( (C%ekc / ( C%rgas*C%trefk) ) *         &
                        ( 1.0 - C%trefk/tlfx(i) ) )

            ! Michaelis menten constant of Rubisco for oxygen:
            conkot(i) = C%conko0 * EXP( ( C%eko / (C%rgas*C%trefk) ) *         &
                        ( 1.0 - C%trefk/tlfx(i) ) )
   
            ! Store leaf temperature
            tlfxx(i) = tlfx(i)
   
            ! "d_{3}" in Wang and Leuning, 1998, appendix E:
            cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
            cx2(i) = 2.0 * C%gam0 * ( 1.0 + C%gam1 * tdiff(i) +                    &
                     C%gam2 * tdiff(i) * tdiff(i ))
    
            ! All equations below in appendix E in Wang and Leuning 1998 are
            ! for calculating anx, csx and gswx for Rubisco limited,
            ! RuBP limited, sink limited
            temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
            temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
            vx3(i,1)  = ej3x(temp2(i,1),ejmxt3(i,1))
            vx3(i,2)  = ej3x(temp2(i,2),ejmxt3(i,2))
    
            temp2(i,1) = rad%qcan(i,1,1) * jtomol * veg%frac4(i)
            temp2(i,2) = rad%qcan(i,2,1) * jtomol * veg%frac4(i)
            vx4(i,1)  = ej4x(temp2(i,1),vcmxt4(i,1))
            vx4(i,2)  = ej4x(temp2(i,2),vcmxt4(i,2))
    
            rdx(i,1) = (C%cfrd3*vcmxt3(i,1) + C%cfrd4*vcmxt4(i,1))*fwsoil(i)  
            rdx(i,2) = (C%cfrd3*vcmxt3(i,2) + C%cfrd4*vcmxt4(i,2))*fwsoil(i)
            
            xleuning(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) )              &
                          * ( ( 1.0 - veg%frac4(i) ) * C%A1C3 / ( 1.0 + dsx(i) &
                          / C%d0c3 ) + veg%frac4(i)    * C%A1C4 / (1.0+dsx(i)/ &
                          C%d0c4) )

            xleuning(i,2) = ( fwsoil(i) / ( csx(i,2)-co2cp3 ) )                &
                            * ( (1.0-veg%frac4(i) ) * C%A1C3 / ( 1.0 + dsx(i) /&
                            C%d0c3 ) + veg%frac4(i)    * C%A1C4 / (1.0+ dsx(i)/&
                            C%d0c4) )
    
         
         ENDIF
         
      ENDDO !i=1,mp
   
      CALL photosynthesis( csx(:,:),                                           &
                           SPREAD( cx1(:), 2, mf ),                            &
                           SPREAD( cx2(:), 2, mf ),                            &
                           gswmin(:,:), rdx(:,:), vcmxt3(:,:),                 &
                           vcmxt4(:,:), vx3(:,:), vx4(:,:),                    &
                           xleuning(:,:), rad%fvlai(:,:),                      &
                           SPREAD( abs_deltlf, 2, mf ),                        &
                           anx(:,:) )

      DO i=1,mp
         
         IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1) Then
      
            DO kk=1,mf
               
               IF(rad%fvlai(i,kk)>C%LAI_THRESH) THEN

                  csx(i,kk) = met%ca(i) - C%RGBWC*anx(i,kk) / (                &
                              gbhu(i,kk) + gbhf(i,kk) )
                  csx(i,kk) = MAX( 1.0e-4, csx(i,kk) )

                  canopy%gswx(i,kk) = MAX( 1.e-3, gswmin(i,kk) +               &
                                      MAX( 0.0, C%RGSWC * xleuning(i,kk) *     &
                                      anx(i,kk) ) )

                  !Recalculate conductance for water:
                  gw(i,kk) = 1.0 / ( 1.0 / canopy%gswx(i,kk) +                 &
                             1.0 / ( 1.075 * ( gbhu(i,kk) + gbhf(i,kk) ) ) )

                  gw(i,kk) = MAX( gw(i,kk), 0.00001 )

                  ! Modified psychrometric constant 
                  ! (Monteith and Unsworth, 1990)
                  psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk) / gw(i,kk) )
           
               ENDIF
            
            ENDDO

            ecx(i) = ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair     &
                     * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) )        &
                     + C%capp * C%rmair * met%dva(i) * ghr(i,1) )              &
                     / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i)       &
                     * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) -  &
                     met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair *      &
                     met%dva(i) * ghr(i,2) ) /                                 &
                     ( air%dsatdk(i) + psycst(i,2) ) 


            IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) Then
               evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels      &
                           / air%rlam(i)

               DO kk = 1,ms
                  
                  ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk),      &
                                        MAX( 0.0, REAL( ssnow%wb(i,kk) ) -     &
                                        1.1 * soil%swilt(i) ) *                &
                                        soil%zse(kk) * 1000.0 )

               ENDDO

               canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
    
               ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))

            ENDIF

            ! Update canopy sensible heat flux:
            hcx(i) = (SUM(rad%rniso(i,:))-ecx(i)                               &
               - C%capp*C%rmair*(met%tvair(i)-met%tk(i))                       &
               * SUM(rad%gradis(i,:)))                                         &
               * SUM(gh(i,:))/ SUM(ghr(i,:))

            ! Update leaf temperature:
            tlfx(i)=met%tvair(i)+REAL(hcx(i))/(C%capp*C%rmair*SUM(gh(i,:)))
      
            ! Update net radiation for canopy:
            rnx(i) = SUM( rad%rniso(i,:)) -                                    &
                     C%CAPP * C%rmair *( tlfx(i)-met%tk(i) ) *                 &
                     SUM( rad%gradis(i,:) )

            ! Update leaf surface vapour pressure deficit:
            dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))

            ! Store change in leaf temperature between successive iterations:
            deltlf(i) = tlfxx(i)-tlfx(i)
            abs_deltlf(i) = ABS(deltlf(i))

         ENDIF !lai/abs_deltlf

      ENDDO !i=1,mp

      ! Whhere leaf temp change b/w iterations is significant, and
      ! difference is smaller than the previous iteration, store results:
      DO i=1,mp
      
         IF ( abs_deltlf(i) < ABS( deltlfy(i) ) ) THEN

            deltlfy(i) = deltlf(i)
            tlfy(i) = tlfx(i)
            rny(i) = rnx(i)
            hcy(i) = hcx(i)
            ecy(i) = ecx(i)
            rdy(i,1) = rdx(i,1)
            rdy(i,2) = rdx(i,2)
            an_y(i,1) = anx(i,1)
            an_y(i,2) = anx(i,2)
            
            ! save last values calculated for ssnow%evapfbl
            oldevapfbl(i,1) = ssnow%evapfbl(i,1)
            oldevapfbl(i,2) = ssnow%evapfbl(i,2)
            oldevapfbl(i,3) = ssnow%evapfbl(i,3)
            oldevapfbl(i,4) = ssnow%evapfbl(i,4)
            oldevapfbl(i,5) = ssnow%evapfbl(i,5)
            oldevapfbl(i,6) = ssnow%evapfbl(i,6)

         ENDIF
          
         IF( abs_deltlf(i) > 0.1 )                                             &
            
            ! after 4 iterations, take mean of current & previous estimates
            ! as the next estimate of leaf temperature, to avoid oscillation
            tlfx(i) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *tlfxx(i) + &
                      ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) )   &
                      * tlfx(i)
     
         IF(k==1) THEN
         
            ! take the first iterated estimates as the defaults
            tlfy(i) = tlfx(i)
            rny(i) = rnx(i)
            hcy(i) = hcx(i)
            ecy(i) = ecx(i)
            rdy(i,:) = rdx(i,:)
            an_y(i,:) = anx(i,:)
            ! save last values calculated for ssnow%evapfbl
            oldevapfbl(i,:) = ssnow%evapfbl(i,:)
         
         END IF
      
      END DO !over mp 

   END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < C%MAXITER)


   ! dry canopy flux
   canopy%fevc = (1.0-canopy%fwet) * ecy

   ! Recalculate ssnow%evapfbl as ecy may not be updated with the ecx
   ! calculated in the last iteration.
   ! DO NOT use simple scaling as there are times that ssnow%evapfbl is zero.
   ! ** ssnow%evapfbl(i,:) = ssnow%evapfbl(i,:) * ecy(i) / ecx(i) **
   DO i = 1, mp
      
       IF( ecy(i) > 0.0 .AND. canopy%fwet(i) < 1.0 ) THEN
         
         IF( ABS( ecy(i) - ecx(i) ) > 1.0e-6 ) THEN
            
            IF( ABS( canopy%fevc(i) - ( SUM( oldevapfbl(i,:)) * air%rlam(i)    &
                /dels ) ) > 1.0e-4 ) THEN
               
               PRINT *, 'Error! oldevapfbl not right.', ktau_gl, i
               PRINT *, 'ecx, ecy = ', ecx(i), ecy(i)
               PRINT *, 'or in mm = ', ecx(i) * ( 1.0 - canopy%fwet(i) )       &
                                       / air%rlam(i) * dels,                   &
                                       ecy(i) * ( 1.0 - canopy%fwet(i) ) /     &
                                       air%rlam(i) * dels

               PRINT *,'fevc = ', canopy%fevc(i), SUM( oldevapfbl(i,:) ) *     &
                                  air%rlam(i) / dels
               PRINT *, 'fwet = ', canopy%fwet(i)
               PRINT *, 'oldevapfbl = ', oldevapfbl(i,:)
               PRINT *, 'ssnow%evapfbl before rescaling: ',                    &
                                                           ssnow%evapfbl(i,:)
               STOP
            
            ELSE
            
               ssnow%evapfbl(i,:) = oldevapfbl(i,:)
            
            END IF
         
         END IF

      END IF
   
   END DO

   canopy%frday = 12.0 * SUM(rdy, 2)
   canopy%fpn = -12.0 * SUM(an_y, 2)
   canopy%evapfbl = ssnow%evapfbl
   
   DEALLOCATE( gswmin )

END SUBROUTINE dryLeaf

! -----------------------------------------------------------------------------

SUBROUTINE photosynthesis( csxz, cx1z, cx2z, gswminz,                          &
                           rdxz, vcmxt3z, vcmxt4z, vx3z,                       &
                           vx4z, xleuningz, vlaiz, deltlfz, anxz )
   USE cable_def_types_mod, only : mp, mf, r_2
   
   REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csxz
   
   REAL, DIMENSION(mp,mf), INTENT(IN) ::                                       &
      cx1z,       & !
      cx2z,       & !     
      gswminz,    & !
      rdxz,       & !
      vcmxt3z,    & !
      vcmxt4z,    & !
      vx4z,       & !
      vx3z,       & !
      xleuningz,  & !
      vlaiz,      & !
      deltlfz       ! 

   REAL, DIMENSION(mp,mf), INTENT(INOUT) :: anxz
   
   ! local variables
   REAL(r_2), DIMENSION(mp,mf) ::                                              &
      coef0z,coef1z,coef2z, ciz,delcxz,                                        &
      anrubiscoz,anrubpz,ansinkz

   REAL, PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
                                       ! Bonan,LSM version 1.0, p106)

   INTEGER :: i,j   
  
   
   DO i=1,mp
      
      IF (sum(vlaiz(i,:)) .GT. C%LAI_THRESH) THEN
      
         DO j=1,mf
            
            IF( vlaiz(i,j) .GT. C%LAI_THRESH .AND. deltlfz(i,j) .GT. 0.1) THEN

               ! Rubisco limited:
               coef2z(i,j) = gswminz(i,j) / C%RGSWC + xleuningz(i,j) *           &
                             ( vcmxt3z(i,j) - ( rdxz(i,j)-vcmxt4z(i,j) ) )

               coef1z(i,j) = (1.0-csxz(i,j)*xleuningz(i,j)) *                  &
                             (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j))             &
                             + (gswminz(i,j)/C%RGSWC)*(cx1z(i,j)-csxz(i,j))      &
                             - xleuningz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0      &
                             + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j) ) )
               
                
               coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) *                 &    
                             (vcmxt3z(i,j)*cx2z(i,j)/2.0                       &
                             + cx1z(i,j)*( rdxz(i,j)-vcmxt4z(i,j ) ) )         &
                             -( gswminz(i,j)/C%RGSWC ) * cx1z(i,j)*csxz(i,j)

               ! kdcorbin,09/10 - new calculations
               IF( ABS(coef2z(i,j)) .GT. 1.0e-9 .AND. &
                   ABS(coef1z(i,j)) .LT. 1.0e-9) THEN
                  
                  ! no solution, give it a huge number as 
                  ! quadratic below cannot handle zero denominator
                  ciz(i,j) = 99999.0        
                  
                  anrubiscoz(i,j) = 99999.0 ! OR do ciz=0 and calc anrubiscoz
                   
               ENDIF

               ! solve linearly
               IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                   ABS( coef1z(i,j) ) >= 1e-9 ) THEN
                  
                  ! same reason as above
                  ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)  
          
                  ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )

                  anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j) / 2.0 ) / &
                                    ( ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) -   &
                                    rdxz(i,j)
   
               ENDIF
   
               ! solve quadratic (only take the more positive solution)
               IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN
                  
                  delcxz(i,j) = coef1z(i,j)**2 -4.0 * coef0z(i,j)              &
                                * coef2z(i,j)
          
                  ciz(i,j) = ( -coef1z(i,j) + SQRT( MAX( 0.0_r_2 ,             &
                             delcxz(i,j) ) ) ) / ( 2.0*coef2z(i,j) )
                             
                  ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )   ! must be positive, why?
   
                  anrubiscoz(i,j) = vcmxt3z(i,j) * ( ciz(i,j) - cx2z(i,j)      &
                                    / 2.0)  / ( ciz(i,j) + cx1z(i,j) ) +       &
                                    vcmxt4z(i,j) - rdxz(i,j)
               
               ENDIF
   
               ! RuBP limited:
               coef2z(i,j) = gswminz(i,j) / C%RGSWC + xleuningz(i,j)             &
                             * ( vx3z(i,j) - ( rdxz(i,j) - vx4z(i,j) ) )
   
               coef1z(i,j) = ( 1.0 - csxz(i,j) * xleuningz(i,j) ) *            &
                             ( vx3z(i,j) + vx4z(i,j) - rdxz(i,j) )             &
                             + ( gswminz(i,j) / C%RGSWC ) *                      &
                             ( cx2z(i,j) - csxz(i,j) ) - xleuningz(i,j)        &
                             * ( vx3z(i,j) * cx2z(i,j) / 2.0 + cx2z(i,j) *     &
                             ( rdxz(i,j) - vx4z(i,j) ) )                          
                             
                             coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) *   &
                             (vx3z(i,j)*cx2z(i,j)/2.0                          &
                             + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j)))                &
                             - (gswminz(i,j)/C%RGSWC)*cx2z(i,j)*csxz(i,j)
   
   
               !kdcorbin, 09/10 - new calculations
               ! no solution, give it a huge number
               IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND.                           &
                   ABS( coef1z(i,j) ) < 1.0e-9 ) THEN

                  ciz(i,j) = 99999.0
                  anrubpz(i,j)  = 99999.0

               ENDIF
   
               ! solve linearly
               IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                   ABS( coef1z(i,j) ) >= 1.e-9) THEN
                  
                  ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                  
                  ciz(i,j) = MAX(0.0_r_2,ciz(i,j))
                  
                  anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /          &
                                 (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
               
               ENDIF
   
               ! solve quadratic (only take the more positive solution)
               IF ( ABS( coef2z(i,j)) >= 1.e-9 ) THEN
                  
                  delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
                  
                  ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j))))     &
                             /(2.0*coef2z(i,j))
                   
                  ciz(i,j) = MAX(0.0_r_2,ciz(i,j)) 
                  
                  anrubpz(i,j)  = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /         &
                                  (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
   
               ENDIF
                 
               ! Sink limited:
               coef2z(i,j) = xleuningz(i,j)
               
               coef1z(i,j) = gswminz(i,j)/C%RGSWC + xleuningz(i,j)               &
                             * (rdxz(i,j) - 0.5*vcmxt3z(i,j))                  &
                             + effc4 * vcmxt4z(i,j) - xleuningz(i,j)           &
                             * csxz(i,j) * effc4 * vcmxt4z(i,j)  
                                            
               coef0z(i,j) = -( gswminz(i,j)/C%RGSWC ) * csxz(i,j) * effc4       &
                             * vcmxt4z(i,j) + ( rdxz(i,j)                      &
                             - 0.5 * vcmxt3z(i,j)) * gswminz(i,j)/C%RGSWC
          
               ! no solution, give it a huge number
               IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND.                           &
                   ABS( coef1z(i,j)) < 1.0e-9 ) THEN

                  ciz(i,j) = 99999.0
                  ansinkz(i,j)  = 99999.0

               ENDIF
          
               ! solve linearly
               IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                   ABS( coef1z(i,j) ) >= 1.e-9 ) THEN

                  ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                  ansinkz(i,j)  = ciz(i,j)

               ENDIF
               
               ! solve quadratic (only take the more positive solution)
               IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN
                  
                  delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
                  
                  ciz(i,j) = (-coef1z(i,j)+SQRT (MAX(0.0_r_2,delcxz(i,j)) ) )  &
                             / ( 2.0 * coef2z(i,j) )
   
                  ansinkz(i,j) = ciz(i,j)
   
               ENDIF
          
               ! minimal of three limited rates
               anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))
        
            ENDIF
      
         ENDDO
    
      ENDIF

   ENDDO
     
END SUBROUTINE photosynthesis

! ------------------------------------------------------------------------------

FUNCTION ej3x(parx,x) RESULT(z)
   
   REAL, INTENT(IN)     :: parx
   REAL, INTENT(IN)     :: x
   REAL                 :: z
   
   z = MAX(0.0,                                                                &
       0.25*((C%alpha3*parx+x-sqrt((C%alpha3*parx+x)**2 -                      &
       4.0*C%convx3*C%alpha3*parx*x)) /(2.0*C%convx3)) )

END FUNCTION ej3x

! ------------------------------------------------------------------------------

FUNCTION ej4x(parx,x) RESULT(z)
   
   REAL, INTENT(IN)     :: parx
   REAL, INTENT(IN)     :: x
   REAL                 :: z
 
   z = MAX(0.0,                                                                &
        (C%alpha4*parx+x-sqrt((C%alpha4*parx+x)**2 -                           &
        4.0*C%convx4*C%alpha4*parx*x))/(2.0*C%convx4))
 
END FUNCTION ej4x

! ------------------------------------------------------------------------------

! Explicit array dimensions as temporary work around for NEC inlining problem
FUNCTION xvcmxt4(x) RESULT(z)
   
   REAL, PARAMETER      :: q10c4 = 2.0
   REAL, INTENT(IN) :: x
   REAL :: z
 
   z = q10c4 ** (0.1 * x - 2.5) /                                              &
        ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))
 
END FUNCTION xvcmxt4

! ------------------------------------------------------------------------------

FUNCTION xvcmxt3(x) RESULT(z)
   
   !  leuning 2002 (p c & e) equation for temperature response
   !  used for vcmax for c3 plants
   REAL, INTENT(IN) :: x
   REAL :: xvcnum,xvcden,z
 
   REAL, PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
   REAL, PARAMETER  :: xVccoef = 1.17461 ! derived parameter
                     ! xVccoef=1.0+exp((EntropJx*C%TREFK-EHdJx)/(Rconst*C%TREFK))
 
   xvcnum=xvccoef*exp( ( ehavc / ( C%rgas*C%TREFK ) )* ( 1.-C%TREFK/x ) )
   xvcden=1.0+exp( ( entropvc*x-ehdvc ) / ( C%rgas*x ) )
   z = max( 0.0,xvcnum / xvcden )

END FUNCTION xvcmxt3

! ------------------------------------------------------------------------------

FUNCTION xejmxt3(x) RESULT(z)
   
   !  leuning 2002 (p c & e) equation for temperature response
   !  used for jmax for c3 plants
 
   REAL, INTENT(IN) :: x
   REAL :: xjxnum,xjxden,z   
 
   REAL, PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
   REAL, PARAMETER  :: xjxcoef = 1.16715 ! derived parameter
 
   xjxnum = xjxcoef*exp( ( ehajx / ( C%rgas*C%TREFK ) ) * ( 1.-C%TREFK / x ) )
   xjxden=1.0+exp( ( entropjx*x-ehdjx) / ( C%rgas*x ) )
   z = max(0.0, xjxnum/xjxden)

END FUNCTION xejmxt3

! ------------------------------------------------------------------------------

SUBROUTINE fwsoil_calc_std(fwsoil, soil, ssnow, veg) 
   USE cable_def_types_mod
   TYPE (soil_snow_type), INTENT(INOUT):: ssnow
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
   REAL, DIMENSION(mp) :: rwater ! soil water availability

   rwater = MAX(1.0e-4_r_2,                                                    &
            SUM(veg%froot * MAX(0.024,MIN(1.0_r_2,ssnow%wb -                   &
            SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))
  
   fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))
      
END SUBROUTINE fwsoil_calc_std 

! ------------------------------------------------------------------------------

SUBROUTINE fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg) 
   USE cable_def_types_mod
   TYPE (soil_snow_type), INTENT(INOUT):: ssnow
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
   REAL, DIMENSION(mp) :: rwater ! soil water availability
   REAL, DIMENSION(mp,3)          :: xi, ti, si
   INTEGER :: j

   rwater = MAX(1.0e-4_r_2,                                                    &
            SUM(veg%froot * MAX(0.024,MIN(1.0_r_2,ssnow%wb -                   &
            SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))

   fwsoil = 1.

   rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)
   
   xi(:,1) = soil%swilt
   xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
   xi(:,3) = soil%sfc
   
   ti(:,1) = 0.
   ti(:,2) = 0.9
   ti(:,3) = 1.0
   
   si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) *                       &
             (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))

   si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) *                       &
             (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))

   si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) *                       &
             (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))

   DO j=1,mp
      IF (rwater(j) < soil%sfc(j) - 0.02)                                      &
         fwsoil(j) = max(0.,min(1., ti(j,1)*si(j,1) +                          &
                       ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))

   ENDDO

END SUBROUTINE fwsoil_calc_non_linear 

! ------------------------------------------------------------------------------

! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
SUBROUTINE fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg) 
   USE cable_def_types_mod
   TYPE (soil_snow_type), INTENT(INOUT):: ssnow
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
   REAL, DIMENSION(mp) :: rwater ! soil water availability
   INTEGER   :: ns
   REAL, parameter ::rootgamma = 0.01   ! (19may2010)
   REAL, DIMENSION(mp)  :: dummy, normFac
   !--- local level dependent rwater 
   REAL, DIMENSION(mp,ms)  :: frwater

   fwsoil(:) = 1.0e-4
   normFac(:) = 0.0

   DO ns=1,ms
     
      dummy(:) = rootgamma/max(1.0e-3,ssnow%wb(:,ns)-soil%swilt(:))

      frwater(:,ns) = MAX(1.0e-4,((ssnow%wb(:,ns)-soil%swilt(:))/soil%ssat(:)) &
                      ** dummy)
      
      fwsoil(:) = min(1.0,max(fwsoil(:),frwater(:,ns)))
      
      normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)

   ENDDO

END SUBROUTINE fwsoil_calc_Lai_Ktaul


    
END MODULE cable_canopy_module
