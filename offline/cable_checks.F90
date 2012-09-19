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
! Purpose: defines ranges to verify validity of inputs and outputs
!          checks mass balance and energy balance
!          switched on/off through namelist variables: check%*
!
! Contact: Bernard.Pak@csiro.au
!
! History: Small change to energy balance equation relative to 1.4b
!          Additional variables from 1.4b for range checking
!
!
!==============================================================================

MODULE cable_checks_module
! Ranges_type in the module sets the acceptable ranges for all variables
! coming in or going out of the offline netcdf driver. The mass_balance
! and energy_balance subroutines calculate cumulative and per-timestep 
! balances, as well as allow user to scrutinise balances in
! particular sections of the code - largely for diagnostics/fault finding.
! rh_sh - converts relative to sensible humidity if met file units require it
!
   USE cable_radiation_module, ONLY: sinbet
   USE cable_def_types_mod

   IMPLICIT NONE

   PRIVATE
   PUBLIC ranges_type, ranges, mass_balance, energy_balance, rh_sh

   TYPE units_type
      CHARACTER(LEN=1) :: Rainf ! 's' (mm/s) or 'h' (mm/h)
      CHARACTER(LEN=1) :: PSurf  ! 'h'(hPa or mbar) or 'P'(Pa)
      CHARACTER(LEN=1) :: Tair  ! 'C' or 'K'
      CHARACTER(LEN=1) :: Qair  ! '%' or 'g' (spec hum)
      CHARACTER(LEN=1) :: CO2air ! 'p' (ppmv)
      CHARACTER(LEN=1) :: Wind ! 'm'(m/s)
   END TYPE units_type
   TYPE(units_type) :: units

   TYPE ranges_type 
      REAL, DIMENSION(2) ::                                               &
           nav_lon = (/-360.0,360.0/),         & 
           nav_lat = (/-90.0,90.0/),           &   
           time,                               &     
           timestp,                            &      
           ! possible forcing variables for CABLE
           SWdown = (/0.0,1360.0/),            & ! W/m^2
           LWdown = (/0.0,750.0/),             & ! W/m^2
           Rainf = (/0.0,0.03/),               & ! mm/s
           Snowf = (/0.0,0.0085/),             & ! mm/s
           PSurf = (/500.0,1100.0/),           & ! mbar/hPa
           Tair = (/200.0,333.0/),             & ! K
           Qair = (/0.0,0.04/),                & ! g/g
           CO2air = (/160.0,2000.0/),          & ! ppmv   
           Wind = (/0.0,75.0/),                & ! m/s
           Wind_N = (/-75.0,75.0/),            & ! m/s
           Wind_E = (/-75.0,75.0/),            & ! m/s
           ! possible output variables
           Qh = (/-1000.0,1000.0/),            & ! W/m^2
           Qle = (/-1000.0,1000.0/),           & ! W/m^2
           Qg = (/-1000.0,1000.0/),            & ! W/m^2   
           SWnet = (/0.0,1350.0/),             & ! W/m^2 (YP oct07)
           ! SWnet = (/0.0,1250.0/),            & ! W/m^2
           LWnet = (/-500.0,510.0/),           & ! W/m^2 
           Rnet = (/-500.0,1250.0/),           & ! W/m^2 
           Evap = (/-0.0003,0.00035/),         &      
           Ewater = (/-0.0003,0.0003/),        &
           ESoil = (/-0.0003,0.0003/),         &
           TVeg = (/-0.0003,0.0003/),          &
           ECanop = (/-0.0003,0.0003/),        &
           PotEvap = (/-0.0006,0.0006/),       &
           ACond = (/0.0,1.0/),                &
           SoilWet = (/-0.4,1.2/),             &
           Albedo = (/0.0,1.0/),               &
           VegT = (/213.0,333.0/),             &
           SoilTemp = (/213.0,343.0/),         &
           SoilMoist = (/0.0,2000.0/),         &
           Qs = (/0.0,5.0/),                   &
           Qsb = (/0.0,5.0/),                  &
           DelSoilMoist  = (/-2000.0,2000.0/), & 
           DelSWE  = (/-2000.0,2000.0/),       &
           DelIntercept = (/-100.0,100.0/),    &
           SnowT  = (/213.0,280.0/),           &
           BaresoilT = (/213.0,343.0/),        &
           AvgSurfT = (/213.0,333.0/),         &
           RadT = (/200.0,373.0/),             &
           SWE = (/0.0,2000.0/),               &
           RootMoist = (/0.0,2000.0/),         &
           CanopInt = (/0.0,100.0/),           &
           NEE = (/-70.0,50.0/),               & ! umol/m2/s
           NPP = (/-20.0,75.0/),               & ! umol/m2/s 
           GPP = (/-20.0,100.0/),              & ! umol/m2/s 
           AutoResp = (/-50.0,20.0/),          & ! umol/m2/s
           LeafResp = (/-50.0,20.0/),          & ! umol/m2/s
           HeteroResp = (/-50.0,20.0/),        & ! umol/m2/s
           HSoil = (/-1000.0,1000.0/),         &
           HVeg = (/-1000.0,1000.0/),          &
           SnowDepth = (/0.0,50.0/),           & ! EK nov07
           Wbal = (/-999999.0,999999.0/),      &
           Ebal = (/-999999.0,999999.0/),      &
           ! parameters:
           albsoil = (/0.0,0.9/),              &
           isoil = (/1.0,30.0/),               &
           iveg = (/1.0,30.0/),                &
           bch = (/2.0,15.0/),                 &
           latitude = (/-90.0,90.0/),          &
           c3 = (/0.0,1.0/),                   & ! EK nov07   
           clay = (/0.0,1.0/),                 &
           css = (/700.0,2200.0/),             &
           rhosoil = (/300.0,3000.0/),         &
           hyds = (/5.0E-7,8.5E-4/),           &
           rs20 = (/0.0,10.0/),                &
           sand = (/0.0,1.0/),                 &
           sfc = (/0.1,0.5/),                  & 
           silt = (/0.0,1.0/),                 &
           ssat = (/0.35,0.5/),                &
           sucs = (/-0.8,-0.03/),              &
           swilt = (/0.05,0.4/),               &
           froot = (/0.0,1.0/),                &
           zse = (/0.0,5.0/),                  &
           canst1 = (/0.05,0.15/),             &
           dleaf = (/0.005,0.4/),              &
           ejmax = (/1.0E-5,3.0E-4/),          &
           frac4 = (/0.0,1.0/),                &
           hc = (/0.0,100.0/),                 &
           lai = (/0.0,8.0/),                  &
           rp20 = (/0.0,10.0/),                &
           vbeta =(/-999999.0,999999.0/),      &
           xalbnir = (/0.0,1.5/),              &
           meth = (/0.0,1.0/),                 &
           za =(/0.0,150.0/),                  &
           rpcoef = (/0.05,1.5/),              &
           shelrb = (/1.0,3.0/),               &
           vcmax = (/5.0E-6,1.5E-4/),          &
           xfang = (/-1.0,0.5/),               &
           ratecp = (/0.01,3.0/),              & 
           ratecs = (/0.01,3.0/),              &
           refsbare = (/0.0,0.5/),             &
           taul = (/0.0,0.3/),                 &
           refl = (/0.0,0.5/),                 &
           tauw = (/0.0,0.1/),                 &
           refw = (/0.0,0.5/),                 &
           extkn = (/0.0,10.0/),               & ! YP oct07
           wai = (/0.0,5.0/),                  & ! YP oct07
           vegcf = (/0.0,100.0/),              & ! YP oct07
           tminvj = (/-20.0,15.0/),            &  
           tmaxvj = (/-15.0,30.0/),            &
           rootbeta = (/0.7,1.0/),             & ! YP oct07
           veg_class = (/1.0,20.0/),           &
           soil_class = (/1.0,20.0/)  
   END TYPE ranges_type
   TYPE(ranges_type),SAVE :: ranges

CONTAINS

!==============================================================================
!
! Name: mass_balance
!
! Purpose: Calculate cumulative and per-timestep balance, as well as allow user
!          to scrutinise balance in particular sections of the code - largely 
!          for diagnostics/fault finding.
!
! CALLed from: write_output
!
!
!==============================================================================

SUBROUTINE mass_balance(dels,ktau, ssnow,soil,canopy,met,                            &
                        air,bal)

   ! Input arguments
   REAL,INTENT(IN)                           :: dels        ! time step size
   INTEGER, INTENT(IN)                       :: ktau        ! timestep number  
   TYPE (soil_snow_type),INTENT(IN)          :: ssnow       ! soil data
   TYPE (soil_parameter_type),INTENT(IN)     :: soil        ! soil data
   TYPE (canopy_type),INTENT(IN)             :: canopy      ! canopy variable data
   TYPE(met_type),INTENT(IN)                 :: met         ! met data
   TYPE (air_type),INTENT(IN)                :: air

   ! Local variables
   REAL(r_2), DIMENSION(:,:,:),POINTER, SAVE :: bwb         ! volumetric soil moisture
   REAL(r_2), DIMENSION(mp)                  :: delwb       ! change in soilmoisture
                                                            ! b/w tsteps
   REAL, DIMENSION(mp)                  :: canopy_wbal !canopy water balance
   TYPE (balances_type),INTENT(INOUT)        :: bal 
   INTEGER                              :: j, k        ! do loop counter
    
   IF(ktau==1) THEN
      ALLOCATE( bwb(mp,ms,2) )
      ! initial vlaue of soil moisture
      bwb(:,:,1)=ssnow%wb
   ELSE
      ! Calculate change in soil moisture b/w timesteps:
      IF(MOD(REAL(ktau),2.0)==1.0) THEN         ! if odd timestep
         bwb(:,:,1)=ssnow%wb
         DO k=1,mp           ! current smoist - prev tstep smoist
            delwb(k) = SUM((bwb(k,:,1)                                         &
                  - (bwb(k,:,2)))*soil%zse)*1000.0
         END DO
      ELSE IF(MOD(REAL(ktau),2.0)==0.0) THEN    ! if even timestep
         bwb(:,:,2)=ssnow%wb
         DO k=1,mp           !  current smoist - prev tstep smoist
            delwb(k) = SUM((bwb(k,:,2)                                         &
                 - (bwb(k,:,1)))*soil%zse)*1000.0
         END DO
      END IF
   END IF

   ! IF(ktau==kend) DEALLOCATE(bwb)

   ! net water into soil (precip-(change in canopy water storage) 
   !  - (change in snow depth) - (surface runoff) - (deep drainage)
   !  - (evaporated water from vegetation and soil(excluding fevw, since
   !      it's included in change in canopy storage calculation))
   ! rml 28/2/11 ! BP changed rnof1+rnof2 to ssnow%runoff which also included rnof5
   ! which is used when nglacier=2 in soilsnow routines (BP feb2011)
   bal%wbal = REAL(met%precip - canopy%delwc - ssnow%snowd+ssnow%osnowd        &
        - ssnow%runoff-(canopy%fevw+canopy%fevc                                &
        + canopy%fes/ssnow%cls)*dels/air%rlam - delwb)
!   bal%wbal = REAL(met%precip - canopy%delwc - ssnow%snowd+ssnow%osnowd        & 
!        - ssnow%rnof1-ssnow%rnof2-(canopy%fevw+canopy%fevc                     &
!        + canopy%fes/ssnow%cls)*dels/air%rlam - delwb)
   ! Canopy water balance: precip-change.can.storage-throughfall-evap+dew
   canopy_wbal = REAL(met%precip-canopy%delwc-canopy%through                   &
        - (canopy%fevw+MIN(canopy%fevc,0.0))*dels/air%rlam)

   bal%wbal_tot = 0. 
   IF(ktau>10) THEN
      ! Add current water imbalance to total imbalance
      ! (method 1 for water balance):
      bal%wbal_tot = bal%wbal_tot + bal%wbal
    
      ! Add to accumulation variables:
      bal%precip_tot = bal%precip_tot + met%precip
      bal%rnoff_tot = bal%rnoff_tot + ssnow%rnof1 + ssnow%rnof2
      bal%evap_tot = bal%evap_tot                                              &
           + (canopy%fev+canopy%fes/ssnow%cls) * dels/air%rlam
   END IF

END SUBROUTINE mass_balance

!==============================================================================
!
! Name: energy_balance
!
! Purpose: Calculate cumulative and per-timestep balance, as well as allow user
!          to scrutinise balance in particular sections of the code - largely 
!          for diagnostics/fault finding.
!
! CALLed from: write_output
!
! MODULEs used: cable_data (inherited) 
!
!==============================================================================

SUBROUTINE energy_balance( dels,met,rad,canopy,bal,ssnow,                    &
                             SBOLTZ,EMLEAF, EMSOIL )

   ! Input arguments
   REAL,INTENT(IN)              :: dels   ! time step size
   TYPE (canopy_type),INTENT(IN)     :: canopy ! canopy variable data
   TYPE(met_type),INTENT(IN)         :: met    ! met data
   TYPE(radiation_type),INTENT(IN)   :: rad    ! met data
   TYPE (balances_type),INTENT(INOUT):: bal 
   TYPE (soil_snow_type),INTENT(IN)  :: ssnow  ! soil data
   REAL, INTENT(IN)::                                                         &
      SBOLTZ,  & !Stefan-Bolzman constant
      EMLEAF,  & !leaf emissivity
      EMSOIL     !leaf emissivity

     
 
   ! SW absorbed + LW absorbed - (LH+SH+ghflux) should = 0
   !bal%ebal = SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs         &
   bal%ebal_cncheck = SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs &
        +met%fld-sboltz*emleaf*canopy%tv**4*(1-rad%transd)                     &
        -sboltz*emsoil*ssnow%tss**4*rad%transd -canopy%fev-canopy%fes          &
        !& -rad%flws*rad%transd -canopy%fev-canopy%fes * ssnow%cls &
        * ssnow%cls -canopy%fh -canopy%ga                ! removed bug (EK 1jul08)
   ! Add to cumulative balance:
   !bal%ebal_tot = bal%ebal_tot + bal%ebal
   bal%ebal_tot_cncheck = bal%ebal_tot_cncheck+ bal%ebal_cncheck

      

!   Atmosphere and radiation code in the global model can 'see' only fluxes, surface 
!   temperature and albedo from CABLE. Those variables as calculated at the end of the given 
!   time step need to be used for the evaluation of the SEB
!   
   !bal%ebal_cncheck = (1.0-rad%albedo(:,1))*met%fsd(:,1) + (1.0-rad%albedo(:,2))*met%fsd(:,2)  &
   bal%ebal = (1.0-rad%albedo(:,1))*met%fsd(:,1) + (1.0-rad%albedo(:,2))       &
        *met%fsd(:,2)+met%fld-sboltz*emleaf*canopy%tv**4*(1-rad%transd)        &
        -sboltz*emsoil*ssnow%tss**4*rad%transd -canopy%fev-canopy%fes          &
        * ssnow%cls-canopy%fh                 
   bal%ebal_tot = bal%ebal_tot + bal%ebal
  
      
   !bal%ebal_tot_cncheck = bal%ebal_tot_cncheck+ bal%ebal_cncheck
!    print 11,bal%ebal_cncheck,bal%ebal_tot_cncheck,bal%ebal,bal%ebal_tot
!11  format(1x,'Ebal_cncheck, Ebal',f6.1,f8.0,2x,f6.1,f8.0)
    
END SUBROUTINE energy_balance

!==============================================================================
!
! Name: rh_sh
!
! Purpose: Converts relative humidity to specific humidity
!
! CALLed from: units_in
!              get_met_data
!
! CALLs: svp
!
!==============================================================================

SUBROUTINE rh_sh (relHum,tk,psurf,specHum)

   ! Input arguments
   REAL, INTENT (IN)  ::                                                  &
        psurf,  & ! surface pressure (hPa)
        relHum, & ! relative humidity (%)
        tk        ! air temp (K) 
   REAL, INTENT (OUT) :: specHum ! specific humidity (kg/kg)
   
   ! Local variables
   REAL ::                                                                &
        es,     & ! saturation vapour pressure
        ws        ! specific humidity at saturation

   es = svp (tk) ! saturation vapour pressure
   ws = 0.622 * es / (psurf - es) ! specific humidity at saturation
   specHum = (relHum/100.0) * ws ! specific humidity

END SUBROUTINE rh_sh

!==============================================================================
!
! Name: svp
!
! Purpose: Calculates saturation vapour pressure
!
! CALLed from: rh_sh
!
!==============================================================================

FUNCTION svp(tk) RESULT (F_Result)

   ! Local variables
   REAL ::                  &
        eilog,                   &
        ewlog,                   &
        ewlog2,                  &
        ewlog3,                  &
        ewlog4,                  &
        F_Result,                &
        temp,                    &
        tk,                      &
        toot,                    &
        toto,                    &
        tsot

   temp = tk - 273.15
   IF (temp < -20.0) THEN
      ! ice saturation
      toot = 273.15 / tk
      toto = 1. / toot
      eilog = -9.09718 * (toot-1) - 3.56654 * (LOG (toot) / LOG (10.0))        &
           + 0.876793 * (1-toto) + (LOG (6.1071) / LOG (10.0))
      F_Result = 10.0**eilog
   ELSE
      tsot = 373.15 / tk
      ewlog = -7.90298 * (tsot-1) + 5.02808 * (LOG (tsot) / LOG (10.0))
      ewlog2 = ewlog - 1.3816e-07 * (10**(11.344 * (1 - (1/tsot))) - 1)
      ewlog3 = ewlog2 + 0.0081328 * (10**(-3.49149 * (tsot-1)) - 1)
      ewlog4 = ewlog3 + (LOG (1013.246) / LOG (10.0))
      F_Result = 10.0**ewlog4
   END IF

END FUNCTION svp
!==============================================================================
END MODULE cable_checks_module
