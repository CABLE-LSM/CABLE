!=========================================================
! purpose: all routines for calculate hydraulic processes in soil
!
! contanct: zihan.lu@inrae.fr
!==========================================================

MODULE cable_soil_hydraulics_module
   
   use cable_data_module,    only: icanopy_type, point2constants
   type(icanopy_type) :: C
   REAL, PARAMETER :: PA_2_MPa = 1E-6
   REAL, PARAMETER                   :: gC2DM = 1./0.49
   PUBLIC :: calc_soil_root_resistance, calc_swp, calc_weighted_swp_and_frac_uptake

CONTAINS
   ! ----------------------------------------------------------------------------
   SUBROUTINE calc_soil_root_resistance(ssnow, soil, veg, casapool, casabiome,root_length_density, i, wbpsdo)
      ! Calculate root & soil hydraulic resistance following SPA approach
      ! (Williams et al.)
      !
      ! Root hydraulic resistance declines linearly with increasing root
      ! biomass according to root resistivity (400) [MPA s m2 mmol-1].
      !
      ! Soil hydraulic resistance depends on soil conductivity, root length,
      ! depth of layer and distance between roots.
      !
      ! In units conversion, useful to recall that:
      ! m s-1 = m3 m-1 m-1 s-1
      ! m3 (amount of water) m-1 (per unit length) m-1 (per unit hydraulic head,
      !                                                 measured in meters) s-1
      !
      ! References:
      ! -----------
      ! * Duursma, R. A. 2008. Predicting the decline in daily maximum
      !   transpiration rate of two pine stands during drought based on
      !   constant minimum leaf water potential and plant hydraulic conductance.
      !   Tree Physiology, 28, 265–276.
      ! * Gardner, W.R. 1964. Relation of root distribution to water uptake
      !   and availability. Agron. J. 56:41–45.
      ! * Newman, E.I. 1969. Resistance to water flow in soil and plant. I.
      !   Soil resistance in relation to amounts of root: theoretical
      !   estimates. J. Appl. Ecol. 6:1–12.
      ! * Williams, M. et al. 1996. Modeling the soil–plant–atmosphere continuum
      !   in a Quercus–Acer stand at Harvard Forest: the regulation of stomatal
      !   conductance by light, nitrogen and soil/plant hydraulic properties.
      !   Plant Cell Environ. 19:911–927.
      !
      ! Martin De Kauwe, 3rd June, 2019

      USE cable_def_types_mod
      USE cable_common_module
      USE casavariable
      use mo_constants, only: pi => pi_dp
      IMPLICIT NONE

      TYPE (soil_snow_type), INTENT(INOUT)        :: ssnow
      TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
      TYPE (casa_pool),  INTENT(IN)           :: casapool
      TYPE (casa_biome),  INTENT(IN)           :: casabiome
      REAL, DIMENSION(:), INTENT(INOUT) :: root_length_density
      real(r_2), dimension(:,:), intent(in), optional :: wbpsdo
      INTEGER, INTENT(IN) :: i
      ! All from Williams et al. 2001, Tree phys
      REAL, PARAMETER :: root_radius = 0.0005                 ! m
      REAL, PARAMETER :: root_xsec_area = pi * root_radius**2 ! m2
      REAL, PARAMETER :: root_density = 0.5e6                ! g biomass m-3 root
      REAL, PARAMETER :: root_resistivity = 25.             ! MPa s g mmol-1, Bonan
      REAL, PARAMETER :: root_k = 100.0

      ! unit conv
      REAL, PARAMETER :: head = 0.009807             ! head of pressure  (MPa/m)
      REAL, PARAMETER :: MM_TO_M = 0.001
      
      REAL, PARAMETER :: G_WATER_TO_MOLE = 1.0 / 18.01528
      REAL, PARAMETER :: CUBIC_M_WATER_2_GRAMS = 1E6
      REAL, PARAMETER :: MOL_2_MMOL = 1000.0
      REAL, PARAMETER :: TINY_NUMBER = 1E-35
      REAL, PARAMETER :: HUGE_NUMBER = 1E35
      REAL, PARAMETER :: BIG_NUMBER = 1E9
      REAL, PARAMETER :: SMALL_NUMBER = 1E-9
      REAL, PARAMETER :: SRA = 48.0 * 1e-3 !m2/gC, specific root area
      REAL, PARAMETER :: Bfr_Bl = 0.26!! ratio of fine root dry mass to leaf dry mass

      REAL, DIMENSION(ms) :: depth
      REAL                :: root_mass, rs, Ksoil0, Ksoil, root_biomass, root_depth, root_mass_density, RAI, Lsr
      REAL                :: soil_resist, rsum, conv, shoot_biomass, leaf_biomass
      real(r_2), dimension(mp,ms) :: wbtmp

      INTEGER :: j,rootRM,soilRM

      REAL               :: ht, stem_biomass
      REAL, PARAMETER    :: Kbiometric = 50.0 ! cst in height-diameter relationship
      REAL, PARAMETER    :: WD = 300.0 ! Wood density kgC/m3
      !REAL, PARAMETER    :: root_conduc = 1e-7 ! kg s-1 Mpa-1 m-1(root length)
      INTEGER, PARAMETER :: STEM_INDEX = 2
      rootRM = 1 ! 1: ; 2: no rootR
      soilRM = 1 ! 1: original, 2: ED
      if (present(wbpsdo)) then
         wbtmp = wbpsdo
      else
         wbtmp = ssnow%wb
      endif
      ! stem_biomass = bgc%cplant(i,STEM_INDEX) * gC2DM
      ! ht = (Kbiometric**(3.0/4.0))*(4.*stem_biomass/(WD*PI))**(1.0/4.0)

      !print*, ht
      call point2constants(C)
      ! convert from gC to g biomass, i.e. twice the C content
      !root_biomass = bgc%cplant(i,ROOT_INDEX) * gC2DM
      !root_biomass = 1443.0 * gC2DM ! EBF value
      ! root_biomass = 832.0 * gC2DM ! Eucface value
      ! root_biomass = casapool%cplant(i,3) * gC2DM !  g m-2
      ! print*, 'root_biomass original: ', root_biomass
      ! another method to calculate root biomass
      ! leaf_biomass = veg%vlai(i) / casabiome%sla(veg%iveg(i)) ! gc m-2
      ! shoot_biomass = ((casapool%cplant(i,1)+casapool%cplant(i,2))/casapool%cplant(i,1)) * leaf_biomass
      ! root_biomass = veg%root_shoot(i) * shoot_biomass * gC2DM
      ! print*, 'root_biomass/leaf biomass: ', root_biomass/leaf_biomass
      !root_biomass = 318.9 * gC2DM ! Spruce experiment
      !! new method to calculate fine root biomass
      leaf_biomass = veg%vlai(i) / casabiome%sla(veg%iveg(i)) ! gc m-2
      root_biomass = Bfr_Bl * leaf_biomass * gC2DM

      ! sensitivity experiment values
      !root_biomass = 200. * gC2DM ! Range from Williams 2001, 200-1000
      !root_biomass = 400. * gC2DM ! Range from Williams 2001, 200-1000
      !root_biomass = 600. * gC2DM ! Range from Williams 2001, 200-1000
      !root_biomass = 800. * gC2DM ! Range from Williams 2001, 200-1000
      !root_biomass = 1000. * gC2DM ! Range from Williams 2001, 200-1000


      ! Always provide a minimum root biomass
      root_biomass = MAX(5., root_biomass)

      ! Store each layers resistance, used in LWP calculatons
      rsum = 0.0

      !print*, bgc%cplant
      !stop
      DO j = 1, ms ! Loop over 6 soil layers

         ! Soil Hydraulic conductivity (m s-1), Campbell 1974
         Ksoil0 = soil%hyds(i) * &
            (max(real(wbtmp(i,j)), soil%sres(i)) / soil%ssat(i))**(2.0 * soil%bch(i) + 3.0)

         ! Calculate soil-root hydraulic resistance

         ! prevent floating point error
         !ksoil = max(ksoil, small_number)
         !IF (Ksoil < TINY_NUMBER) THEN
         ! IF (Ksoil < SMALL_NUMBER) THEN
         !    !ssnow%soilR(i,j) = HUGE_NUMBER
         !    ssnow%soilR(i,j) = BIG_NUMBER / C%RHOW
         !    rsum = rsum + ( 1.0 / ssnow%soilR(i,j) )
         ! ELSE

            ! Root biomass density (g biomass m-3 soil)
            ! Divide root mass up by the frac roots in the layer (g m-3)
            ! plant carbon is g m-3
            root_mass_density = root_biomass * veg%froot(i,j) / soil%zse(j)
            ! Root length density (m root m-3 soil)
            root_length_density(j) = root_mass_density / (root_density * root_xsec_area)
            if (rootRM==1) then
            ! kg m-2 s-1 Mpa-1
            ssnow%rootR(i,j) = 1.0 / (veg%root_conduc(i) * root_length_density(j) * soil%zse(j)) 
            elseif (rootRM==2) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  method2 for rootR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ssnow%rootR(i,j) = SMALL_NUMBER
            endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!method 1 for soilR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Conductance of the soil-to-root pathway can be estimated
            ! assuming that the root system consists of one long root that
            ! has access to a surrounding cylinder of soil
            ! (Gardner 1960, Newman 1969)
            if (soilRM == 1) then
            rs = SQRT(1.0 / (root_length_density(j) * pi))
            ! converts from m s-1 to m2 s-1 MPa-1
            Ksoil = Ksoil0 / (C%grav * C%RHOW * PA_2_MPa )
            ! converts from m2 s-1 MPa-1 to kg s-1 Mpa-1 m-1
            Ksoil = Ksoil * C%RHOW
            ! if (j==1) then
            !    print*, 'ksoil kg s-1 Mpa-1 m-1',Ksoil
            ! endif
            ! Soil-to-root resistance (MPa s m2 m-3)
            soil_resist = LOG(rs / root_radius) / &
               (2.0 * pi * root_length_density(j) * soil%zse(j) * Ksoil)
            ! if (j==1) then
            !       print*, 'rs/root_radius',1.0/LOG(rs / root_radius)
            ! endif 
            !! convert from MPa s m2 m-3 to MPa s m2 kg-1
            !soil_resist = soil_resist / C%RHOW
            ! if (j==1) then
            !    print*, 'ksoil kg m-2 Mpa-1 s-1',1.0/soil_resist
            ! endif
            endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! method 2 for soilR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!! reference: Katul, G., Leuning, R., & Oren, R. (2003). Relationship between plant hydraulic and 
            !!biochemical properties derived from a steady-state coupled water and carbon transport model.
            !!! Plant, Cell & Environment, 26(3), 339–350. https://doi.org/10.1046/j.1365-3040.2003.00965.x
            if (soilRM==2) then
               Ksoil = Ksoil0 * C%RHOW ! kg m-2 s-1
               ksoil  = ksoil / (C%grav * C%RHOW * PA_2_MPa )! kg m-1 Mpa-1 s-1
            RAI = root_biomass / gC2DM * SRA * veg%froot(i,j) !! m2/m2
            Lsr = pi * soil%zse(j) / sqrt(RAI)
            soil_resist = Lsr / Ksoil ! kg m-2 Mpa-1 s-1
            ! if (j==1) then
            !    print*, 'ED method: ksoil',1.0/soil_resist
            ! endif
            endif
            ! root_resistance is commented out : don't use root-component of
            ! resistance (is part of plant resistance)
            ! MPa s m2 mmol-1 H2O
            !ssnow%soilR(i,j) = soil_resist !+ root_resist
            ssnow%soilR(i,j) = soil_resist


         ! END IF

         !print*, "DEBUG soilR:", j,  ssnow%soilR(i,j), Ksoil, root_mass, root_biomass , rsum


         ! Need to combine resistances in parallel, but we only want the
         ! soil term as the root component is part of the plant resistance
         rsum = rsum + ( 1.0 / ssnow%soilR(i,j) )

      END DO

      ! rsum calc above is 1/rsum (i.e. conductance(, as we're combining in
      ! parallel, turn back into resistance (MPa s m2 kg-1 H2O)
      ssnow%Rsr(i) = 1.0 / rsum


   END SUBROUTINE calc_soil_root_resistance
   ! ----------------------------------------------------------------------------

   ! ----------------------------------------------------------------------------
   SUBROUTINE calc_swp(ssnow, soil, i, wbpsdo)
      ! Calculate the soil water potential.
      !
      ! Martin De Kauwe, 2019; Manon Sabot, 2022
      !
      ! When the soil water potential is below the wilting point, we apply the
      ! dry soil correction from Webb 2000, as parameterised in Schneider & Goss
      ! 2012.
      !
      ! References:
      !
      ! * Webb, S. W. (2000). A simple extension of two-phase characteristic
      ! curves to include the dry region. Water Resources Research, 36(6),
      ! 1425-1430.
      !
      ! * Schneider, M., & Goss, K. U. (2012). Prediction of water retention
      ! curves for dry soils from an established pedotransfer function:
      ! Evaluation of the Webb model. Water Resources Research, 48(6).

      USE cable_def_types_mod
      USE cable_common_module

      IMPLICIT NONE

      TYPE (soil_snow_type), INTENT(INOUT)        :: ssnow
      TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
      real(r_2), dimension(:,:), intent(in), optional :: wbpsdo

      INTEGER             :: j
      INTEGER, INTENT(IN) :: i
      REAL                :: psi_sat, psi_wilt
      real(r_2), dimension(mp,ms) :: wbtmp
      REAL, PARAMETER :: cmH2O_TO_MPa = 1.0 / (10.0 * 1036)
      if (present(wbpsdo)) then
         wbtmp = wbpsdo
      else
         wbtmp = ssnow%wb
      endif
      call point2constants(C)
      ssnow%psi_soil(i,:) = 0.0 ! MPa

      ! Soil matric potential at saturation (m of head to MPa: 9.8 * KPA_2_MPA)
      psi_sat = soil%sucs(i) * C%grav * C%RHOW *PA_2_MPa

      ! Soil matric potential at wilting point, MPa
      psi_wilt = psi_sat * MAX(1.E-9, MIN(1.0, soil%swilt(i) / soil%ssat(i))) &
         ** (-soil%bch(i))

      DO j = 1, ms ! Loop over 6 soil layers

         ! Below the wilting point, the water potential drops to silly values, due
         ! to the non-linear reln btw swc and swp. This is really more of an
         ! implentation issue. The very thin upper layers would obv get extremely
         ! negative, which in turn would bias the weighted_swp.
         ! The bounding here can be considered as a physical disconnection of the
         ! roots from the soil.

         ssnow%psi_soil(i,j) = psi_sat * MAX(1.E-9, MIN(1.0, real(wbtmp(i,j)) / &
            soil%ssat(i))) ** (-soil%bch(i))

         ! bound psi_soil by the wilting point in upper soil layers
         IF (j < 3) THEN

            ssnow%psi_soil(i,j) = MAX(psi_wilt, ssnow%psi_soil(i,j))

            ! below wilting point, apply dry soil correction from Webb 2000
         ELSE IF (wbtmp(i,j) < soil%swilt(i)) THEN

            ssnow%psi_soil(i,j) = (-10.0 ** ((LOG10(-psi_wilt / cmH2O_TO_MPa) &
               - 6.8) * wbtmp(i,j) / soil%swilt(i) + 6.8)) &
               * cmH2O_TO_MPa

         END IF

      END DO

   END SUBROUTINE calc_swp
   ! ----------------------------------------------------------------------------

   ! ----------------------------------------------------------------------------
   SUBROUTINE calc_weighted_swp_and_frac_uptake(ssnow, soil, canopy, veg, &
      root_length, i)
      !
      ! Determine weighted SWP given the the maximum rate of water supply from
      ! each rooted soil layer and hydraulic resistance of each layer. We are
      ! also calculating a weighting fraction for water extraction. This is
      ! achieved by roughly estimating the maximum rate of water supply from each
      ! rooted soil layer, using SWP and hydraulic resistance of each layer.
      ! Actual water from each layer is determined using the estimated value as a
      ! weighted factor.
      !
      ! Martin De Kauwe, 2019

      USE cable_def_types_mod
      USE cable_common_module

      IMPLICIT NONE

      TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow
      TYPE (soil_parameter_type), INTENT(INOUT) :: soil
      TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables
      TYPE(veg_parameter_type), INTENT(IN)      :: veg

      ! The minimum root water potential (MPa), used in determining fractional
      ! water uptake in soil layers
      REAL, PARAMETER :: min_root_wp = -3.0

      REAL, DIMENSION(ms)            :: swp, est_evap
      REAL, DIMENSION(:), INTENT(IN) :: root_length
      REAL                           :: total_est_evap, swp_diff, depth_sum

      INTEGER, INTENT(IN) :: i
      INTEGER             :: j

      ! SPA method to figure out relative water uptake.
      LOGICAL :: SPA_relative_uptake
      SPA_relative_uptake = .TRUE.

      total_est_evap = 0.0
      est_evap = 0.0


      ! Estimate max transpiration from gradient-gravity / soil resistance
      DO j = 1, ms ! Loop over 6 soil layers

         IF (ssnow%soilR(i,j) .GT. 0.0 .AND. veg%froot(i,j) .GT. 0.0) THEN

            est_evap(j) = MAX(0.0, &
               (ssnow%psi_soil(i,j) - min_root_wp) / ssnow%soilR(i,j))
         ELSE
            est_evap(j) = 0.0 ! when no roots present
         ENDIF

         ! No uptake from frozen soils
         IF ( ssnow%wbice(i,j) .gt. 0.0 ) THEN
            est_evap(j) = 0.0
         ENDIF

    
            ssnow%psi_rootzone(i) = 0.0
            ssnow%fraction_uptake(i,:) = 0.0
            IF (veg%froot(i,j) .GT. 0.0) THEN
               ! Soil water potential weighted by layer Emax (from SPA)
               ssnow%psi_rootzone(i) = ssnow%psi_rootzone(i) + &
                  ssnow%psi_soil(i,j) * est_evap(j)
            ENDIF
         

      END DO
      total_est_evap = SUM(est_evap)


     
         ! calculate the weighted psi_soil, ms8355 2022
         IF (total_est_evap > 0.0) THEN

            ! Soil water potential is weighted by total_est_evap.
            ssnow%psi_rootzone(i) = ssnow%psi_rootzone(i) / total_est_evap

         ELSE

            ssnow%psi_rootzone(i) = 0.0
            depth_sum = 0.0

            DO j = 1, ms ! Loop over 6 soil layers

               IF (veg%froot(i,j) .GT. 0.0) THEN

                  ssnow%psi_rootzone(i) = ssnow%psi_rootzone(i) + &
                     ssnow%psi_soil(i,j) * soil%zse(j)
                  depth_sum = depth_sum + soil%zse(j)

               ENDIF

            END DO

            ssnow%psi_rootzone(i) = ssnow%psi_rootzone(i) / depth_sum

         END IF

         ! SPA method to figure out relative water uptake.
         ! Fraction uptake in each layer by Emax in each layer
         IF (SPA_relative_uptake) THEN

            IF (total_est_evap > 0.0) THEN
               DO j = 1, ms ! Loop over 6 soil layers

                  ! fraction of water taken from layer, I've lower bounded frac
                  ! uptake because when soilR is set to a huge number
                  ! (see calc_soil_root_resistance), then frac_uptake will be so
                  ! small you end up with numerical issues.
                  ssnow%fraction_uptake(i,j) = MAX(1E-09, &
                     est_evap(j) / total_est_evap)

                  IF ((ssnow%fraction_uptake(i,j) > 1.0) .OR. &
                     (ssnow%fraction_uptake(i,j) < 0.0)) THEN
                     PRINT *, 'Problem with the uptake fraction (either >1 or 0<)'
                     STOP
                  END IF

               END DO
            ELSE
               ! No water was evaporated
               ssnow%fraction_uptake(i,:) = 1.0 / FLOAT(ms)
            END IF

            ! Use Taylor-Keppler root water uptake distribution.
         ELSE
            ! Taylor and Keppler: relative water uptake is
            ! proportional to root length density and Psi difference.
            ! See : Taylor, H.M. and B. Keppler. 1975. Water uptake by cotton root
            ! systems: an examination of assumptions in the single root model.
            ! Soil Science. 120:57-67.
            DO j = 1, ms ! Loop over 6 soil layers

               IF (total_est_evap .GT. 0.) THEN
                  swp_diff = MAX(0., (ssnow%psi_soil(i,j) - min_root_wp))
                  ssnow%fraction_uptake(i,j) = root_length(j) * swp_diff
               ELSE
                  ! no water uptake possible
                  ssnow%fraction_uptake(i,j) = 0.0
               END IF
            END DO

            IF (SUM(ssnow%fraction_uptake) .GT. 0) THEN
               ! Make sure that it sums to 1.
               ssnow%fraction_uptake = ssnow%fraction_uptake / &
                  SUM(ssnow%fraction_uptake)
            ELSE
               ssnow%fraction_uptake = 0.0
            ENDIF

         ENDIF
   


   END SUBROUTINE calc_weighted_swp_and_frac_uptake
   ! ----------------------------------------------------------------------------
   SUBROUTINE calc_psix(ssnow, soil, canopy, veg, casapool, ex, psix, kplant, i)
      USE cable_def_types_mod
      USE cable_common_module
      use cable_veg_hydraulics_module, only: get_xylem_vulnerability
      USE casavariable
      use mo_constants, only: pi => pi_sp
      IMPLICIT NONE

      TYPE (soil_snow_type), INTENT(INOUT)        :: ssnow
      TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
      TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables
      TYPE(casa_pool),        INTENT(IN)    :: casapool
      real(r_2), INTENT(IN)    :: ex !kg m-2 s-1 
      real(r_2), INTENT(OUT)    :: psix
      real(r_2), INTENT(OUT)    :: kplant
      INTEGER, INTENT(IN) :: i
      REAL, DIMENSION(ms) ::  layer_depth, zsetmp, froottmp, frcuptmp
      REAL :: SoilMoistPFTtemp

      REAL, DIMENSION(ms) :: a
      INTEGER :: k
      real :: k1, k2, pd, BAI, WD, AGB_pl, DBH, plc, sumpsiksoil, sumksoil
      real:: huber_value, SLA
      !huber_value = 3.65e-4 ! m2 m-2 (sapwood area/leaf area) for Fagus sylvatica
      !SLA = 154 ! cm2 g-1 for Fagus sylvatica
      !real(r_2) :: psi_sr
      !ELSE
      ! zihanlu: calculate psi_soil no matter which soil_sche is used
      !DO i = 1, mp
      !CALL calc_swp(ssnow, soil, i)
      !write(logn,*),'calculate soilMoistPFT'
      call point2constants(C)

      layer_depth(1) = 0.0_r_2
      do k=2, ms
         layer_depth(k) = sum(soil%zse(1:k-1))
      enddo


         ! zsetmp = soil%zse
         ! where (layer_depth > veg%zr(i))
         !    zsetmp = 0.
         ! elsewhere
         !    zsetmp = min(real(veg%zr(i)) - layer_depth, soil%zse)
         ! endwhere


         ! if (cable_user%calSoilMean == 'zW') then

         !    SoilMoistPFTtemp = sum(real(ssnow%wb(i,:)) * zsetmp, 1) / sum(zsetmp)
         !    ssnow%psi_rootzone(i) = soil%sucs(i) * 9.8 * 0.001 * MAX(1.E-9, MIN(1.0, SoilMoistPFTtemp / &
         !       soil%ssat(i))) ** (-soil%bch(i))

         ! elseif (cable_user%calSoilMean == 'frootW') then
         !    froottmp = veg%froot(i,:)
         !    SoilMoistPFTtemp = sum(real(ssnow%wb(i,:)) * froottmp * zsetmp) / sum(froottmp * zsetmp)
         !    ssnow%psi_rootzone(i) = soil%sucs(i) * 9.8 * 0.001 * MAX(1.E-9, MIN(1.0, SoilMoistPFTtemp / &
         !       soil%ssat(i))) ** (-soil%bch(i))
         ! elseif (cable_user%calSoilMean == 'FrcUpW') then
         !    frcuptmp = ssnow%fraction_uptake(i,:)
         !    where (layer_depth > veg%zr(i) )
         !       frcuptmp = 0
         !    endwhere
         !    SoilMoistPFTtemp = sum(real(ssnow%wb(i,:)) * frcuptmp) / sum(frcuptmp)
         !    ssnow%psi_rootzone(i) = soil%sucs(i) * 9.8 * 0.001 * MAX(1.E-9, MIN(1.0, SoilMoistPFTtemp / &
         !       soil%ssat(i))) ** (-soil%bch(i))

         ! endif
         !SoilMoistPFTtemp = real(sum(ssnow%wb * 1000.0_r_2 * real(spread(soil%zse,1,mp),r_2),2),r_2)
         sumksoil = 0.0
         sumpsiksoil = 0.0
         DO k = 1, ms
            sumpsiksoil = sumpsiksoil + ssnow%psi_soil(i,k) / (ssnow%soilR(i,k) + ssnow%rootR(i,k) )
            sumksoil = sumksoil + 1.0 / (ssnow%soilR(i,k) + ssnow%rootR(i,k) )
            ! if (k==1) then
            !    print*, 'calc_psix: psi_soil ',ssnow%psi_soil(i,k)
            !    print*, 'calc_psix: ksoil, rootR ',1.0/ssnow%soilR(i,k), ssnow%rootR(i,k)
            ! endif
         END DO
         psix = (sumpsiksoil - ex) / sumksoil
         ! print*, 'calc_psix: psix ',psix
         !canopy%psix(i) = psi_sr - ex / veg%kmax(i)
         !psix = psi_sr - ex / veg%kmax(i)
         ! Plant hydraulic conductance (mmol m-2 leaf s-1 MPa-1)
         ! k1 = 50.0_r_2
         ! k2 = 1.5_r_2
         k1 = 0.2351 ! coefficient in 
         k2 = 2.3226
         WD = 300.0_r_2 ! kgC m-2
         !Biomass = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7_r_2
         !pd = 4.0_r_2 * real(casapool%cplant(i,2) / 1000.0_r_2,r_2) / (WD * pi * veg%hc * (veg%hc / k1)**(2.0_r_2*k2))
         pd = 0.07_r_2 ! pl m-2
         !DBH = (veg%hc/k1)**k2
         AGB_pl = casapool%cplant(i,2) / 1000.0_r_2 * gC2DM / pd ! kg pl-1
         DBH = (AGB_pl/k1)**(1.0_r_2/k2) ! cm 
         BAI = (DBH/200.0_r_2)**2.0_r_2*pi*pd ! m2 m-2
         ! print*, 'old BAI',BAI
         ! plc = get_xylem_vulnerability(ssnow%psi_rootzone(i), &
         ! veg%b_plant(i), veg%c_plant(i))
         !BAI = casapool%cplant(i,1) * gC2DM * SLA * huber_value !  m2 m-2
         huber_value = veg%huber_value(i)
         BAI = veg%vlai(i) * huber_value !  m2 m-2 
         ! print*, 'new BAI',BAI
         plc = get_xylem_vulnerability(psix, &
         veg%b_plant(i), veg%c_plant(i))
         !canopy%kplant(i) = veg%kmax(i) * BAI / veg%hc(i) * plc
         kplant = veg%kmax(i) * BAI / veg%hc(i) * plc
         !print*,'Entry kplant: ',casapool%cplant(i,2)/1000.0_r_2,AGB_pl,DBH,BAI,plc,canopy%kplant(i)

   
   END SUBROUTINE calc_psix


END MODULE cable_soil_hydraulics_module
