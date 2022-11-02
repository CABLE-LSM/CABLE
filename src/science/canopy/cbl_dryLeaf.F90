MODULE cbl_dryLeaf_module

IMPLICIT NONE

PUBLIC :: dryLeaf
PRIVATE

CONTAINS

  SUBROUTINE dryLeaf( dels, rad, rough, air, met,                                &
       veg, canopy, soil, ssnow, dsx,                             &
       fwsoil, tlfx,  tlfy,  ecy, hcy,                            &
       rny, gbhu, gbhf, csx,                                      &
       cansat, ghwet, iter,climate )

    USE cable_def_types_mod
    USE cable_common_module
USE cbl_photosynthesis_module,  ONLY : photosynthesis
USE cbl_fwsoil_module,        ONLY : fwsoil_calc_std, fwsoil_calc_non_linear,           &
                                     fwsoil_calc_Lai_Ktaul, fwsoil_calc_sli
! maths & other constants
USE cable_other_constants_mod, ONLY : CLAI_THRESH  => LAI_THRESH
! physical constants
USE cable_phys_constants_mod, ONLY : CTFRZ   => TFRZ
USE cable_phys_constants_mod, ONLY : CDHEAT  => DHEAT
USE cable_phys_constants_mod, ONLY : CRGAS   => RGAS
USE cable_phys_constants_mod, ONLY : CCAPP   => CAPP
USE cable_phys_constants_mod, ONLY : CRMAIR  => RMAIR
! photosynthetic constants
USE cable_photo_constants_mod, ONLY : CMAXITER  => MAXITER ! only integer here
USE cable_photo_constants_mod, ONLY : CTREFK => TREFK
USE cable_photo_constants_mod, ONLY : CGAM0  => GAM0
USE cable_photo_constants_mod, ONLY : CGAM1  => GAM1
USE cable_photo_constants_mod, ONLY : CGAM2  => GAM2
USE cable_photo_constants_mod, ONLY : CRGSWC => RGSWC
USE cable_photo_constants_mod, ONLY : CRGBWC => RGBWC

    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (air_type),       INTENT(INOUT) :: air
    TYPE (met_type),       INTENT(INOUT) :: met
    TYPE (canopy_type),    INTENT(INOUT) :: canopy
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow

    TYPE (veg_parameter_type),  INTENT(INOUT)   :: veg
    TYPE (soil_parameter_type), INTENT(inout)   :: soil
    TYPE (climate_type), INTENT(IN)    :: climate

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
         rnx,        & ! net rad prev timestep
         fwsoil_coef   !

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
                                ! Ticket #56, xleuning is no longer used, we replace it with
                                ! gs_coeff,
                                ! which is computed differently based on the new GS_SWITCH. If GS_SWITCH
                                ! is "leuning", it's the same, if "medlyn", then the new Medlyn model
                                ! xleuning,   & ! leuning stomatal coeff
         gs_coeff,   & ! stom coeff, Ticket #56
         psycst,     & ! modified pych. constant
         frac42,     & ! 2D frac4
         temp2

    REAL, DIMENSION(:,:), POINTER :: gswmin ! min stomatal conductance

    REAL, DIMENSION(mp,2) ::  gsw_term, lower_limit2  ! local temp var

    INTEGER :: i, j, k, kk  ! iteration count
    REAL :: vpd, g1 ! Ticket #56
#define VanessasCanopy
#ifdef VanessasCanopy
    REAL, DIMENSION(mp,mf)  ::                                                  &
         xleuning    ! leuning stomatal coeff
#endif

    REAL :: medlyn_lim  !INH 2018: should be a parameter in long-term
    ! END header

    ALLOCATE( gswmin(mp,mf ))

    ! Soil water limitation on stomatal conductance:
    IF( iter ==1) THEN
       IF ((cable_user%soil_struc=='default').AND.(cable_user%FWSOIL_SWITCH.NE.'Haverd2013')) THEN
          IF(cable_user%FWSOIL_SWITCH == 'standard') THEN
             CALL fwsoil_calc_std( fwsoil, soil, ssnow, veg)
          ELSEIF (cable_user%FWSOIL_SWITCH == 'non-linear extrapolation') THEN
             !EAK, 09/10 - replace linear approx by polynomial fitting
             CALL fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg)
          ELSEIF(cable_user%FWSOIL_SWITCH == 'Lai and Ktaul 2000') THEN
             CALL fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg)
          ELSE
             STOP 'fwsoil_switch failed.'
          ENDIF
          canopy%fwsoil = fwsoil
       ELSEIF ((cable_user%soil_struc=='sli').OR.(cable_user%FWSOIL_SWITCH=='Haverd2013')) THEN
          fwsoil = canopy%fwsoil
       ENDIF

    ENDIF

    ! weight min stomatal conductance by C3 an C4 plant fractions
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
    gsw_term = SPREAD(veg%gswmin,2,mf)
    lower_limit2 = rad%scalex * gsw_term
    gswmin = MAX(1.e-6,lower_limit2)


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

       IF(canopy%vlaiw(kk) <= CLAI_THRESH) THEN
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
    DO WHILE (k < CMAXITER)
       k = k + 1
       DO i=1,mp

          IF (canopy%vlaiw(i) > CLAI_THRESH .AND. abs_deltlf(i) > 0.1) THEN

             ghwet(i) = 2.0   * sum_gbh(i)
             gwwet(i) = 1.075 * sum_gbh(i)
             ghrwet(i) = sum_rad_gradis(i) + ghwet(i)

             ! Calculate lat heat from wet canopy, may be neg.
             ! if dew on wet canopy to avoid excessive evaporation:
             ccfevw(i) = MIN(canopy%cansto(i) * air%rlam(i) / dels,             &
                  2.0 / (1440.0 / (dels/60.0)) * air%rlam(i) )

             ! Grashof number (Leuning et al, 1995) eq E4:
             gras(i) = MAX(1.0e-6,                                              &
                  1.595E8* ABS( tlfx(i)-met%tvair(i))* (veg%dleaf(i)**3.0) )

             ! See Appendix E in (Leuning et al, 1995):
             gbhf(i,1) = rad%fvlai(i,1) * air%cmolar(i) * 0.5*Cdheat           &
                  * ( gras(i)**0.25 ) / veg%dleaf(i)
             gbhf(i,2) = rad%fvlai(i,2) * air%cmolar(i) * 0.5 * Cdheat         &
                  * ( gras(i)**0.25 ) / veg%dleaf(i)
             gbhf(i,:) = MAX( 1.e-6_r_2, gbhf(i,:) )

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
             temp(i) = xvcmxt4(tlfx(i)-Ctfrz) * veg%vcmax(i) * veg%frac4(i)
             vcmxt4(i,1) = rad%scalex(i,1) * temp(i)
             vcmxt4(i,2) = rad%scalex(i,2) * temp(i)

             ! Leuning 2002 (P C & E) equation for temperature response
             ! used for Jmax for C3 plants:
             temp(i) = xejmxt3(tlfx(i)) * veg%ejmax(i) * (1.0-veg%frac4(i))
             ejmxt3(i,1) = rad%scalex(i,1) * temp(i)
             ejmxt3(i,2) = rad%scalex(i,2) * temp(i)

             ! Difference between leaf temperature and reference temperature:
             tdiff(i) = tlfx(i) - CTREFK

             ! Michaelis menten constant of Rubisco for CO2:
             conkct(i) = veg%conkc0(i) * EXP( ( veg%ekc(i) / (Crgas*Ctrefk) ) &
                  * ( 1.0 - Ctrefk/tlfx(i) ) )

             ! Michaelis menten constant of Rubisco for oxygen:
             conkot(i) = veg%conko0(i) * EXP( ( veg%eko(i) / (Crgas*Ctrefk) ) &
                  * ( 1.0 - Ctrefk/tlfx(i) ) )

             ! Store leaf temperature
             tlfxx(i) = tlfx(i)

             ! "d_{3}" in Wang and Leuning, 1998, appendix E:
             cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
             cx2(i) = 2.0 * Cgam0 * ( 1.0 + Cgam1 * tdiff(i)                  &
                  + Cgam2 * tdiff(i) * tdiff(i) )

             ! All equations below in appendix E in Wang and Leuning 1998 are
             ! for calculating anx, csx and gswx for Rubisco limited,
             ! RuBP limited, sink limited
             temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
             temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
             vx3(i,1)  = ej3x(temp2(i,1),veg%alpha(i),veg%convex(i),ejmxt3(i,1))
             vx3(i,2)  = ej3x(temp2(i,2),veg%alpha(i),veg%convex(i),ejmxt3(i,2))
             temp2(i,1) = rad%qcan(i,1,1) * jtomol * veg%frac4(i)
             temp2(i,2) = rad%qcan(i,2,1) * jtomol * veg%frac4(i)
             vx4(i,1)  = ej4x(temp2(i,1),veg%alpha(i),veg%convex(i),vcmxt4(i,1))
             vx4(i,2)  = ej4x(temp2(i,2),veg%alpha(i),veg%convex(i),vcmxt4(i,2))

             rdx(i,1) = (veg%cfrd(i)*Vcmxt3(i,1) + veg%cfrd(i)*vcmxt4(i,1))
             rdx(i,2) = (veg%cfrd(i)*vcmxt3(i,2) + veg%cfrd(i)*vcmxt4(i,2))

             !Vanessa - the trunk does not contain xleauning as of Ticket#56 inclusion
             !as well as other inconsistencies here that need further investigation. In the
             !interests of getting this into the trunk ASAP just isolate this code for now
             !default side of this condition is to use trunk version

             !#ifdef VanessasCanopy


             IF (cable_user%CALL_climate) THEN

                ! Atkins et al. 2015, Table S4,
                ! modified by saling factor to reduce leaf respiration to
                ! expected proportion of GPP
                !Broad-leaved trees: Rdark a25 =
                !1.2818 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)
                !C3 herbs/grasses: Rdark,a25 =
                !1.6737 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)
                !Needle-leaved trees: Rdark,a25 =
                !1.2877 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)
                !Shrubs: Rdark,a25 = 1.5758 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)

                IF (veg%iveg(i).EQ.2 .OR. veg%iveg(i).EQ. 4  ) THEN ! broadleaf forest

                   rdx(i,1) = 0.60*(1.2818e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)

                ELSEIF (veg%iveg(i).EQ.1 .OR. veg%iveg(i).EQ. 3  ) THEN ! needleleaf forest
                   rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)

                ELSEIF (veg%iveg(i).EQ.6 .OR. veg%iveg(i).EQ.8 .OR. &
                     veg%iveg(i).EQ. 9  ) THEN ! C3 grass, tundra, crop
                   rdx(i,1) = 0.60*(1.6737e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)

                ELSE  ! shrubs and other (C4 grass and crop)
                   rdx(i,1) = 0.60*(1.5758e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)
                ENDIF


                ! modify for leaf area and instanteous temperature response (Rd25 -> Rd)
                rdx(i,1) = rdx(i,1) * xrdt(tlfx(i)) * rad%scalex(i,1)
                rdx(i,2) = rdx(i,2) * xrdt(tlfx(i)) * rad%scalex(i,2)



                ! reduction of daytime leaf dark-respiration to account for
                !photo-inhibition
                !Mercado, L. M., Huntingford, C., Gash, J. H. C., Cox, P. M.,
                ! and Jogireddy, V.:
                ! Improving the representation of radiation
                !interception and photosynthesis for climate model applications,
                !Tellus B, 59, 553-565, 2007.
                ! Equation 3
                ! (Brooks and Farquhar, 1985, as implemented by Lloyd et al., 1995).
                ! Rc = Rd 0 < Io < 10 μmol quantam−2s−1
                ! Rc = [0.5 − 0.05 ln(Io)] Rd Io > 10μmol quantam−2s−1

                IF (jtomol*1.0e6*rad%qcan(i,1,1).GT.10.0) &
                     rdx(i,1) = rdx(i,1) * &
                     (0.5 - 0.05*LOG(jtomol*1.0e6*rad%qcan(i,1,1)))

                IF (jtomol*1.0e6*rad%qcan(i,1,2).GT.10.0) &
                     rdx(i,2) = rdx(i,2) * &
                     (0.5 - 0.05*LOG(jtomol*1.0e6*rad%qcan(i,1,2)))

!$                xleuning(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) )              &
!$                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))
!$                xleuning(i,2) = ( fwsoil(i) / ( csx(i,2) - co2cp3 ) )              &
!$                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

             ELSE !cable_user%call_climate

!$!Vanessa:note there is no xleuning to go into photosynthesis etc anymore
!$             gs_coeff = xleuning

                !#else
                rdx(i,1) = (veg%cfrd(i)*vcmxt3(i,1) + veg%cfrd(i)*vcmxt4(i,1))
                rdx(i,2) = (veg%cfrd(i)*vcmxt3(i,2) + veg%cfrd(i)*vcmxt4(i,2))

             ENDIF !cable_user%call_climate

             ! Ticket #56 added switch for Belinda Medlyn's model
             IF (cable_user%GS_SWITCH == 'leuning') THEN
                gs_coeff(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) )              &
                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

                gs_coeff(i,2) = ( fwsoil(i) / ( csx(i,2) - co2cp3 ) )              &
                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

                ! Medlyn BE et al (2011) Global Change Biology 17: 2134-2144.
             ELSEIF(cable_user%GS_SWITCH == 'medlyn') THEN

                gswmin = veg%g0(i)

                IF (dsx(i) < 50.0) THEN
                   vpd  = 0.05 ! kPa
                ELSE
                   vpd = dsx(i) * 1E-03 ! Pa -> kPa
                END IF

                g1 = veg%g1(i)

                gs_coeff(i,1) = (1.0 + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                gs_coeff(i,2) = (1.0 + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)

                !INH 2018: enforce gs_coeff to vary proportionally to fwsoil in dry soil conditions
                ! required to avoid transpiration without soil water extraction
                medlyn_lim = 0.05
                IF (fwsoil(i) <= medlyn_lim) THEN
                   gs_coeff(i,1) = (fwsoil(i) / medlyn_lim + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                   gs_coeff(i,2) = (fwsoil(i) / medlyn_lim + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)
                END IF

             ELSE
                STOP 'gs_model_switch failed.'
             ENDIF ! IF (cable_user%GS_SWITCH == 'leuning') THEN
             !#endif

          ENDIF !IF (canopy%vlaiw(i) > CLAI_THRESH .AND. abs_deltlf(i) > 0.1)

       ENDDO !i=1,mp

       CALL photosynthesis( csx(:,:),                                           &
            SPREAD( cx1(:), 2, mf ),                            &
            SPREAD( cx2(:), 2, mf ),                            &
            gswmin(:,:), rdx(:,:), vcmxt3(:,:),                 &
            vcmxt4(:,:), vx3(:,:), vx4(:,:),                    &
                                ! Ticket #56, xleuning replaced with gs_coeff here
            gs_coeff(:,:), rad%fvlai(:,:),&
            SPREAD( abs_deltlf, 2, mf ),                        &
            anx(:,:), fwsoil(:) )

       DO i=1,mp


          IF (canopy%vlaiw(i) > CLAI_THRESH .AND. abs_deltlf(i) > 0.1 ) THEN

             DO kk=1,mf

                IF(rad%fvlai(i,kk)>CLAI_THRESH) THEN

                   csx(i,kk) = met%ca(i) - CRGBWC*anx(i,kk) / (                &
                        gbhu(i,kk) + gbhf(i,kk) )
                   csx(i,kk) = MAX( 1.0e-4_r_2, csx(i,kk) )


                   ! Ticket #56, xleuning replaced with gs_coeff here
                   canopy%gswx(i,kk) = MAX( 1.e-3, gswmin(i,kk)*fwsoil(i) +     &
                        MAX( 0.0, CRGSWC * gs_coeff(i,kk) *     &
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

             ecx(i) = ( air%dsatdk(i) * ( rad%rniso(i,1) - Ccapp * Crmair     &
                  * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) )        &
                  + Ccapp * Crmair * met%dva(i) * ghr(i,1) )              &
                  / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i)       &
                  * ( rad%rniso(i,2) - Ccapp * Crmair * ( met%tvair(i) -  &
                  met%tk(i) ) * rad%gradis(i,2) ) + Ccapp * Crmair *      &
                  met%dva(i) * ghr(i,2) ) /                                 &
                  ( air%dsatdk(i) + psycst(i,2) )

             IF (cable_user%fwsoil_switch=='Haverd2013') THEN
                ! avoid root-water extraction when fwsoil is zero
                IF (fwsoil(i).LT.1e-6) THEN
                   anx(i,:) =  - rdx(i,:)
                   ecx(i) = 0.0
                ENDIF

                canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))

                CALL getrex_1d(ssnow%wbliq(i,:),&
                     ssnow%rex(i,:), &
                     canopy%fwsoil(i), &
                     REAL(veg%froot(i,:),r_2),&
                     soil%ssat_vec(i,:), &
                     soil%swilt_vec(i,:), &
                     MAX(REAL(canopy%fevc(i)/air%rlam(i)/1000_r_2,r_2),0.0_r_2), &
                     REAL(veg%gamma(i),r_2), &
                     REAL(soil%zse,r_2), REAL(dels,r_2), REAL(veg%zr(i),r_2))

                fwsoil(i) = canopy%fwsoil(i)
                ssnow%evapfbl(i,:) = ssnow%rex(i,:)*dels*1000_r_2 ! mm water &
                !(root water extraction) per time step

             ELSE

                IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
                   evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels      &
                        / air%rlam(i)

                   DO kk = 1,ms

                      ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk),      &
                           MAX( 0.0, REAL( ssnow%wb(i,kk) ) -     &
                           1.1 * soil%swilt(i) ) *                &
                           soil%zse(kk) * 1000.0 )

                   ENDDO
                   IF (cable_user%soil_struc=='default') THEN
                      canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
                      ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                   ELSEIF (cable_user%soil_struc=='sli') THEN
                      canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))
                   ENDIF

                ENDIF

             ENDIF
             ! Update canopy sensible heat flux:
             hcx(i) = (SUM(rad%rniso(i,:))-ecx(i)                               &
                  - Ccapp*Crmair*(met%tvair(i)-met%tk(i))                       &
                  * SUM(rad%gradis(i,:)))                                         &
                  * SUM(gh(i,:))/ SUM(ghr(i,:))
             ! Update leaf temperature:
             tlfx(i)=met%tvair(i)+REAL(hcx(i))/(Ccapp*Crmair*SUM(gh(i,:)))

             ! Update net radiation for canopy:
             rnx(i) = SUM( rad%rniso(i,:)) -                                    &
                  CCAPP * Crmair *( tlfx(i)-met%tk(i) ) *                 &
                  SUM( rad%gradis(i,:) )

             ! Update leaf surface vapour pressure deficit:
             dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))

             dsx(i)=  MAX(dsx(i),0.0)

             ! Store change in leaf temperature between successive iterations:
             deltlf(i) = tlfxx(i)-tlfx(i)
             abs_deltlf(i) = ABS(deltlf(i))

          ENDIF !lai/abs_deltlf

       ENDDO !i=1,mp
       ! Where leaf temp change b/w iterations is significant, and
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
             oldevapfbl(i,:) = ssnow%evapfbl(i,:)

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


    END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < CMAXITER)


    ! dry canopy flux
    canopy%fevc = (1.0-canopy%fwet) * ecy

    IF (cable_user%fwsoil_switch.NE.'Haverd2013') THEN

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
                   ! STOP

                ELSE

                   ssnow%evapfbl(i,:) = oldevapfbl(i,:)

                END IF

             END IF

          END IF

       END DO

    ENDIF

    canopy%frday = 12.0 * SUM(rdy, 2)
    ! vh ! inserted min to avoid -ve values of GPP
    canopy%fpn = MIN(-12.0 * SUM(an_y, 2), canopy%frday)
    canopy%evapfbl = ssnow%evapfbl


    DEALLOCATE( gswmin )

  END SUBROUTINE dryLeaf
 
 
   SUBROUTINE getrex_1d(theta, rex, fws, Fs, thetaS, thetaw, Etrans, gamma, dx, dt, zr)

    ! root extraction : Haverd et al. 2013
    USE cable_def_types_mod, ONLY: r_2

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:), INTENT(IN)    :: theta      ! volumetric soil moisture
    REAL(r_2), DIMENSION(:), INTENT(INOUT)   :: rex    ! water extraction per layer
    REAL(r_2),               INTENT(INOUT) :: fws    ! stomatal limitation factor
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: Fs     ! root length density
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaS ! saturation soil moisture
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaw ! soil moisture at wiliting point
    REAL(r_2),               INTENT(IN)    :: Etrans ! total transpiration
    REAL(r_2),               INTENT(IN)    :: gamma  ! skew of Li & Katul alpha2 function
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: dx     ! layer thicknesses (m)
    REAL(r_2),               INTENT(IN)    :: dt
    REAL(r_2),               INTENT(IN)    :: zr

    ! Gets rate of water extraction compatible with CABLE stomatal conductance model
    ! theta(:) - soil moisture(m3 m-3)
    ! rex(:)   - rate of water extraction by roots from layers (cm/h).
    REAL(r_2), DIMENSION(1:SIZE(theta)) ::  lthetar, alpha_root, delta_root, layer_depth
    REAL(r_2)                       :: trex, e3, one, zero
    INTEGER :: k

    e3 = 0.001
    one = 1.0;
    zero = 0.0;

    layer_depth(1) = 0.0
    DO k=2,SIZE(theta)
       layer_depth(k) = SUM(dx(1:k-1))
    ENDDO

    !theta(:)   = S(:)*thetaS(:)
    lthetar(:) = LOG(MAX(theta(:)-thetaw(:),e3)/thetaS(:))

    WHERE ((theta(:)-thetaw(:)) > e3)
       alpha_root(:) = EXP( gamma/MAX(theta(:)-thetaw(:), e3) * lthetar(:) )
    ELSEWHERE
       alpha_root(:) = zero
    endwhere

    WHERE (Fs(:) > zero .AND. layer_depth < zr )  ! where there are roots and we are aobe max rooting depth
       delta_root(:) = one
    ELSEWHERE
       delta_root(:) = zero
    endwhere

    rex(:) = alpha_root(:)*Fs(:)

    trex = SUM(rex(:))
    IF (trex > zero) THEN
       rex(:) = rex(:)/trex
    ELSE
       rex(:) = zero
    ENDIF
    rex(:) = Etrans*rex(:)


    ! reduce extraction efficiency where total extraction depletes soil moisture below wilting point
    WHERE (((rex*dt) > (theta(:)-thetaw(:))*dx(:)) .AND. ((rex*dt) > zero))
       alpha_root = alpha_root*(theta(:)-thetaw(:))*dx(:)/(1.1_r_2*rex*dt)
    endwhere
    rex(:) = alpha_root(:)*Fs(:)

    trex = SUM(rex(:))
    IF (trex > zero) THEN
       rex(:) = rex(:)/trex
    ELSE
       rex(:) = zero
    ENDIF
    rex(:) = Etrans*rex(:)

    ! check that the water available in each layer exceeds the extraction
    !if (any((rex*dt) > (theta(:)-0.01_r_2)*dx(:))) then
    IF (ANY(((rex*dt) > MAX((theta(:)-thetaw(:)),zero)*dx(:)) .AND. (Etrans > zero))) THEN
       fws = zero
       ! distribute extraction according to available water
       !rex(:) = (theta(:)-0.01_r_2)*dx(:)
       rex(:) = MAX((theta(:)-thetaw(:))*dx(:),zero)
       trex = SUM(rex(:))
       IF (trex > zero) THEN
          rex(:) = rex(:)/trex
       ELSE
          rex(:) = zero
       ENDIF
       rex(:) = Etrans*rex(:)
    ELSE
       fws    = MAXVAL(alpha_root(2:)*delta_root(2:))
    ENDIF

  END SUBROUTINE getrex_1d


  FUNCTION ej3x(parx,alpha,convex,x) RESULT(z)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: alpha
    REAL, INTENT(IN)     :: convex
    REAL, INTENT(IN)     :: x
    REAL                 :: z

    z = MAX(0.0,                                                                &
         0.25*((alpha*parx+x-SQRT((alpha*parx+x)**2 -                      &
         4.0*convex*alpha*parx*x)) /(2.0*convex)) )
  END FUNCTION ej3x


  FUNCTION ej4x(parx,alpha,convex,x) RESULT(z)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: alpha
    REAL, INTENT(IN)     :: convex
    REAL, INTENT(IN)     :: x
    REAL                 :: z

    z = MAX(0.0,                                                                &
         (alpha*parx+x-SQRT((alpha*parx+x)**2 -                           &
         4.0*convex*alpha*parx*x))/(2.0*convex))

  END FUNCTION ej4x

  ! Explicit array dimensions as temporary work around for NEC inlining problem
  FUNCTION xvcmxt4(x) RESULT(z)

    REAL, PARAMETER      :: q10c4 = 2.0
    REAL, INTENT(IN) :: x
    REAL :: z

    z = q10c4 ** (0.1 * x - 2.5) /                                              &
         ((1.0 + EXP(0.3 * (13.0 - x))) * (1.0 + EXP(0.3 * (x - 36.0))))

  END FUNCTION xvcmxt4

  ! ------------------------------------------------------------------------------

  FUNCTION xvcmxt3(x) RESULT(z)

USE cable_phys_constants_mod, ONLY : CRGAS   => RGAS
USE cable_photo_constants_mod, ONLY : CTREFK => TREFK
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for vcmax for c3 plants
    REAL, INTENT(IN) :: x
    REAL :: xvcnum,xvcden,z

    REAL, PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
    REAL, PARAMETER  :: xVccoef = 1.17461 ! derived parameter
    ! xVccoef=1.0+exp((EntropJx*CTREFK-EHdJx)/(Rconst*CTREFK))

    xvcnum=xvccoef*EXP( ( ehavc / ( Crgas*CTREFK ) )* ( 1.-CTREFK/x ) )
    xvcden=1.0+EXP( ( entropvc*x-ehdvc ) / ( Crgas*x ) )
    z = MAX( 0.0,xvcnum / xvcden )

  END FUNCTION xvcmxt3

  ! ------------------------------------------------------------------------------
  REAL FUNCTION xrdt(x)

    !  Atkins et al. (Eq 1, New Phytologist (2015) 206: 614–636)
    !variable Q10 temperature of dark respiration
    ! Originally from Tjoelker et al. 2001

    REAL, INTENT(IN) :: x


    xrdt = (3.09 - 0.043*((x-273.15)+25.)/2.0)**((x-273.15 -25.0)/10.0)

  END FUNCTION xrdt

  ! ------------------------------------------------------------------------------

  FUNCTION xejmxt3(x) RESULT(z)
USE cable_phys_constants_mod, ONLY : CRGAS   => RGAS
USE cable_photo_constants_mod, ONLY : CTREFK => TREFK

    !  leuning 2002 (p c & e) equation for temperature response
    !  used for jmax for c3 plants

    REAL, INTENT(IN) :: x
    REAL :: xjxnum,xjxden,z

    REAL, PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
    REAL, PARAMETER  :: xjxcoef = 1.16715 ! derived parameter

    xjxnum = xjxcoef*EXP( ( ehajx / ( Crgas*CTREFK ) ) * ( 1.-CTREFK / x ) )
    xjxden=1.0+EXP( ( entropjx*x-ehdjx) / ( Crgas*x ) )
    z = MAX(0.0, xjxnum/xjxden)

  END FUNCTION xejmxt3



END MODULE cbl_dryLeaf_module
