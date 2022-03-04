MODULE cable_pft_params_mod
  IMPLICIT NONE

  TYPE vegin_type

     REAL, DIMENSION(:),ALLOCATABLE ::                                        &
          canst1,     & !
          dleaf,      & !
          length,     & !
          width,      & !
          vcmax,      & !
          ejmax,      & !
          hc,         & !
          xfang,      & !
          rp20,       & !
          rpcoef,     & !
          rs20,       & !
          wai,        & !
          rootbeta,   & !
          shelrb,     & !
          vegcf,      & !
          frac4,      & !
          xalbnir,    & !
          extkn,      & !
          tminvj,     & !
          tmaxvj,     & !
          vbeta,      &
          a1gs,       &
          d0gs,       &
          alpha,      &
          convex,     &
          cfrd,       &
          gswmin,     &
          conkc0,     &
          conko0,     &
          ekc,        &
          eko,        &
          g0,         & !  Ticket #56
          g1,         & !  Ticket #56
          zr,         &
          clitt,      &
          sf,         & ! mgk576
          psi_f ,     & ! mgk576
          X_hyd,      & ! mgk576
          p50,        & ! mgk576MPa
          Kmax,       & ! mgk576MPa +1
          Kcrit,       & ! mgk576MPa +2
          b_plant,    & ! mgk576MPa +3
          c_plant,    & ! mgk576MPa +4
          s50,        & ! mgk576MPa
          kp_sat,     & ! mgk576MPa
          Cl,         & ! mgk576MPa
          Cs,         & ! mgk576MPa
          gmin          ! mgk576MPa

     REAL, DIMENSION(:,:),ALLOCATABLE ::                                      &
          froot,      & !
          cplant,     & !
          csoil,      & !
          ratecp,     & !
          ratecs,     & !
          refl,     & !
          taul        !

  END TYPE vegin_type

  CHARACTER(LEN=70), DIMENSION(:), POINTER ::                                 &
       veg_desc    ! decriptions of veg type

  TYPE(vegin_type),  SAVE  :: vegin

CONTAINS


  SUBROUTINE cable_pft_params()

    ! Gets parameter values for each vegetation type
    USE cable_def_types_mod, ONLY : mvtype, ms, ncs, ncp, nrb

    INTEGER :: a, jveg ! do loop counter
    LOGICAL, SAVE :: first_call = .TRUE.
    !mvtype=17
    !mvtype=22 ! mgk576
    mvtype=32 ! mgk576, 26 Jun 2021, profit max
    ! Allocate memory for type-specific vegetation parameters:
    IF( first_call ) THEN

       ALLOCATE (                                                               &
            vegin%canst1( mvtype ), vegin%dleaf( mvtype ),                        &
            vegin%length( mvtype ), vegin%width( mvtype ),                        &
            vegin%vcmax( mvtype ),  vegin%ejmax( mvtype ),                        &
            vegin%hc( mvtype ), vegin%xfang( mvtype ),                            &
            vegin%rp20( mvtype ), vegin%rpcoef( mvtype ),                         &
            vegin%rs20( mvtype ), vegin%wai( mvtype ),                            &
            vegin%rootbeta( mvtype ), vegin%shelrb( mvtype ),                     &
            vegin%vegcf( mvtype ), vegin%frac4( mvtype ),                         &
            vegin%xalbnir( mvtype ), vegin%extkn( mvtype ),                       &
            vegin%tminvj( mvtype ), vegin%tmaxvj( mvtype ),                       &
            vegin%vbeta( mvtype ), vegin%froot( ms, mvtype ),                     &
            vegin%cplant( ncp, mvtype ), vegin%csoil( ncs, mvtype ),              &
            vegin%ratecp( ncp, mvtype ), vegin%ratecs( ncs, mvtype ),             &
            vegin%refl( nrb, mvtype ), vegin%taul( nrb, mvtype ),                 &
            veg_desc( mvtype ),                                                   &
            vegin%a1gs(mvtype), vegin%d0gs(mvtype),                               &
            vegin%alpha(mvtype),vegin%convex(mvtype),vegin%cfrd(mvtype),          &
            vegin%gswmin(mvtype),vegin%conkc0(mvtype), vegin%conko0(mvtype),      &
            vegin%ekc(mvtype), vegin%eko(mvtype),                                 &
                                ! Ticket #56
            vegin%g0( mvtype ), vegin%g1( mvtype ),                               &
                                !! vh_veg_params !!
            vegin%zr(mvtype), vegin%clitt(mvtype) ,                               &
            !mgk576
            vegin%sf(mvtype), vegin%psi_f(mvtype) ,                               &
            vegin%X_hyd(mvtype), vegin%p50(mvtype), vegin%Kmax(mvtype),           &
            vegin%Kcrit(mvtype), vegin%b_plant(mvtype),                           &
            vegin%c_plant(mvtype), vegin%s50(mvtype),                             &
            vegin%kp_sat(mvtype), vegin%Cl(mvtype), vegin%Cs(mvtype),           &
            vegin%gmin(mvtype))

       !PFT parameters: description and corresponding variable name in code.
       !PFT parameters are assigned as TYPE vegin% but later used as veg%

       !PFT: evergreen_needleleaf
       !=========================================================

       ! Hydraulics
       vegin%Kmax(1) = 0.8
       vegin%Kcrit(1) = 0.04 ! vegin%Kmax(1) * 0.05
       vegin%b_plant(1) = 4.102123226
       vegin%c_plant(1) = 4.347350048
       vegin%gmin(1) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(1) = 0.000112
       vegin%ejmax(1) = 0.000187

       vegin%canst1(1) =        0.100000
       vegin%length(1) =        0.055000
       vegin%width(1) =        0.001000
       !vegin%vcmax(1) =        0.000040
       !vegin%ejmax(1) =        0.000000
       vegin%hc(1) =       17.000000
       vegin%xfang(1) =        0.010000
       vegin%rp20(1) =        3.000000
       vegin%rpcoef(1) =        0.083200
       vegin%rs20(1) =        1.000000
       vegin%wai(1) =        1.000000
       vegin%rootbeta(1) =        0.943000
       vegin%shelrb(1) =        2.000000
       vegin%vegcf(1) =        9.000000
       vegin%frac4(1) =        0.000000
       vegin%xalbnir(1) =        1.000000
       vegin%extkn(1) =        0.001000
       vegin%tminvj(1) =      -15.000000
       vegin%tmaxvj(1) =      -10.000000
       vegin%vbeta(1) =        2.000000
       vegin%froot(1,1) =        0.050000
       vegin%froot(2,1) =        0.050000
       vegin%froot(3,1) =        0.050000
       vegin%froot(4,1) =        0.050000
       vegin%froot(5,1) =        0.050000
       vegin%froot(6,1) =        0.050000
       vegin%refl(1,1) =        0.090000
       vegin%taul(1,1) =        0.090000
       vegin%refl(2,1) =        0.300000
       vegin%taul(2,1) =        0.300000
       vegin%refl(3,1) =        0.010000
       vegin%taul(3,1) =        0.010000
       vegin%csoil(1,1) =      184.000000
       vegin%ratecs(1,1) =        2.000000
       vegin%csoil(2,1) =      367.000000
       vegin%ratecs(2,1) =        0.500000
       vegin%cplant(1,1) =      200.000000
       vegin%ratecp(1,1) =        1.000000
       vegin%cplant(2,1) =    10217.000000
       vegin%ratecp(2,1) =        0.030000
       vegin%cplant(3,1) =      876.000000
       vegin%ratecp(3,1) =        0.140000
       vegin%a1gs(1) =        9.000000
       vegin%d0gs(1) =     1500.000000
       vegin%alpha(1) =        0.200000
       vegin%convex(1) =        0.700000
       vegin%cfrd(1) =        0.015000
       vegin%gswmin(1) =        0.010000
       vegin%conkc0(1) =        0.000302
       vegin%conko0(1) =        0.256000
       vegin%ekc(1) =    59430.000000
       vegin%eko(1) =    36000.000000
       vegin%g0(1) =        0.000000
       vegin%g1(1) =        2.346064
       vegin%zr(1) =        1.800000
       vegin%clitt(1) =       20.000000


       !PFT: evergreen_broadleaf
       !=========================================================
       vegin%canst1(2) =        0.100000
       vegin%length(2) =        0.100000
       vegin%width(2) =        0.050000
       vegin%vcmax(2) =        0.000055
       vegin%ejmax(2) =        0.000000
       vegin%hc(2) =       35.000000
       vegin%xfang(2) =        0.100000
       vegin%rp20(2) =        0.600000
       vegin%rpcoef(2) =        0.083200
       vegin%rs20(2) =        1.000000
       vegin%wai(2) =        1.000000
       vegin%rootbeta(2) =        0.962000
       vegin%shelrb(2) =        2.000000
       vegin%vegcf(2) =       14.000000
       vegin%frac4(2) =        0.000000
       vegin%xalbnir(2) =        1.000000
       vegin%extkn(2) =        0.001000
       vegin%tminvj(2) =      -15.000000
       vegin%tmaxvj(2) =      -10.000000
       vegin%vbeta(2) =        2.000000
       vegin%froot(1,2) =        0.200000
       vegin%froot(2,2) =        0.200000
       vegin%froot(3,2) =        0.200000
       vegin%froot(4,2) =        0.200000
       vegin%froot(5,2) =        0.200000
       vegin%froot(6,2) =        0.200000
       vegin%refl(1,2) =        0.090000
       vegin%taul(1,2) =        0.090000
       vegin%refl(2,2) =        0.290000
       vegin%taul(2,2) =        0.290000
       vegin%refl(3,2) =        0.010000
       vegin%taul(3,2) =        0.010000
       vegin%csoil(1,2) =      303.000000
       vegin%ratecs(1,2) =        2.000000
       vegin%csoil(2,2) =      606.000000
       vegin%ratecs(2,2) =        0.500000
       vegin%cplant(1,2) =      300.000000
       vegin%ratecp(1,2) =        1.000000
       vegin%cplant(2,2) =    16833.000000
       vegin%ratecp(2,2) =        0.030000
       vegin%cplant(3,2) =     1443.000000
       vegin%ratecp(3,2) =        0.140000
       vegin%a1gs(2) =        9.000000
       vegin%d0gs(2) =     1500.000000
       vegin%alpha(2) =        0.200000
       vegin%convex(2) =        0.700000
       vegin%cfrd(2) =        0.015000
       vegin%gswmin(2) =        0.010000
       vegin%conkc0(2) =        0.000302
       vegin%conko0(2) =        0.256000
       vegin%ekc(2) =    59430.000000
       vegin%eko(2) =    36000.000000
       vegin%g0(2) =        0.000000
       vegin%g1(2) =        4.114762
       vegin%zr(2) =        3.000000
       vegin%clitt(2) =        6.000000

       ! mgk576, hydraulics stuff
       !vegin%g1(2) = 12.0  ! JED 15 for teretoconis seedlings, using Jim's value for Eucface
       !vegin%sf(2) = 8.0
       !vegin%psi_f(2) = -2.0
       !vegin%X_hyd(2) = 50.0
       !vegin%p50(2) = -4.
       !vegin%s50(2) = 30.0
       !vegin%kp_sat(2) = 4.0
       !vegin%Cl(2) = 10000.
       !vegin%Cs(2) = 120000.

       ! mgk576, hydraulics stuff
       vegin%gmin(2) = 1.275279  ! mmol m-2 s-1
       !vegin%vcmax(2) = 0.000085
       !vegin%ejmax(2) = 0.000142
       vegin%g1(2) = 3.154297 ! 4.12
       vegin%sf(2) = 2.000000
       vegin%psi_f(2) = -2.455474
       vegin%X_hyd(2) = 50.0
       vegin%p50(2) = -3.002384
       vegin%s50(2) = 35.26948
       !vegin%kp_sat(2) = 1.686987 ! 4.0
       vegin%Cl(2) = 342.904821
       vegin%Cs(2) = 53266.089926

       vegin%Kmax(2) = 1.5
       vegin%Kcrit(2) = 0.075 ! vegin%Kmax(18) * 0.05

       !PFT: deciduous_needleleaf
       !=========================================================
       vegin%canst1(3) =        0.100000
       vegin%length(3) =        0.040000
       vegin%width(3) =        0.001000
       vegin%vcmax(3) =        0.000040
       vegin%ejmax(3) =        0.000000
       vegin%hc(3) =       15.500000
       vegin%xfang(3) =        0.010000
       vegin%rp20(3) =        3.000000
       vegin%rpcoef(3) =        0.083200
       vegin%rs20(3) =        1.000000
       vegin%wai(3) =        1.000000
       vegin%rootbeta(3) =        0.966000
       vegin%shelrb(3) =        2.000000
       vegin%vegcf(3) =        9.000000
       vegin%frac4(3) =        0.000000
       vegin%xalbnir(3) =        1.000000
       vegin%extkn(3) =        0.001000
       vegin%tminvj(3) =        5.000000
       vegin%tmaxvj(3) =       10.000000
       vegin%vbeta(3) =        2.000000
       vegin%froot(1,3) =        0.200000
       vegin%froot(2,3) =        0.200000
       vegin%froot(3,3) =        0.200000
       vegin%froot(4,3) =        0.200000
       vegin%froot(5,3) =        0.200000
       vegin%froot(6,3) =        0.200000
       vegin%refl(1,3) =        0.075000
       vegin%taul(1,3) =        0.075000
       vegin%refl(2,3) =        0.300000
       vegin%taul(2,3) =        0.300000
       vegin%refl(3,3) =        0.010000
       vegin%taul(3,3) =        0.010000
       vegin%csoil(1,3) =      107.000000
       vegin%ratecs(1,3) =        2.000000
       vegin%csoil(2,3) =      214.000000
       vegin%ratecs(2,3) =        0.500000
       vegin%cplant(1,3) =      200.000000
       vegin%ratecp(1,3) =        1.000000
       vegin%cplant(2,3) =     5967.000000
       vegin%ratecp(2,3) =        0.030000
       vegin%cplant(3,3) =      511.000000
       vegin%ratecp(3,3) =        0.140000
       vegin%a1gs(3) =        9.000000
       vegin%d0gs(3) =     1500.000000
       vegin%alpha(3) =        0.200000
       vegin%convex(3) =        0.700000
       vegin%cfrd(3) =        0.015000
       vegin%gswmin(3) =        0.010000
       vegin%conkc0(3) =        0.000302
       vegin%conko0(3) =        0.256000
       vegin%ekc(3) =    59430.000000
       vegin%eko(3) =    36000.000000
       vegin%g0(3) =        0.000000
       vegin%g1(3) =        2.346064
       vegin%zr(3) =        2.000000
       vegin%clitt(3) =       10.000000

       !PFT: deciduous_broadleaf
       !=========================================================
       vegin%canst1(4) =        0.100000
       vegin%length(4) =        0.150000
       vegin%width(4) =        0.080000
       vegin%vcmax(4) =        0.000060
       vegin%ejmax(4) =        0.000000
       vegin%hc(4) =       20.000000
       vegin%xfang(4) =        0.250000
       vegin%rp20(4) =        2.200000
       vegin%rpcoef(4) =        0.083200
       vegin%rs20(4) =        1.000000
       vegin%wai(4) =        1.000000
       vegin%rootbeta(4) =        0.961000
       vegin%shelrb(4) =        2.000000
       vegin%vegcf(4) =        8.000000
       vegin%frac4(4) =        0.000000
       vegin%xalbnir(4) =        1.000000
       vegin%extkn(4) =        0.001000
       vegin%tminvj(4) =        5.000000
       vegin%tmaxvj(4) =       15.000000
       vegin%vbeta(4) =        2.000000
       vegin%froot(1,4) =        0.200000
       vegin%froot(2,4) =        0.200000
       vegin%froot(3,4) =        0.200000
       vegin%froot(4,4) =        0.200000
       vegin%froot(5,4) =        0.200000
       vegin%froot(6,4) =        0.200000
       vegin%refl(1,4) =        0.090000
       vegin%taul(1,4) =        0.090000
       vegin%refl(2,4) =        0.290000
       vegin%taul(2,4) =        0.290000
       vegin%refl(3,4) =        0.010000
       vegin%taul(3,4) =        0.010000
       vegin%csoil(1,4) =      216.000000
       vegin%ratecs(1,4) =        2.000000
       vegin%csoil(2,4) =      432.000000
       vegin%ratecs(2,4) =        0.500000
       vegin%cplant(1,4) =      300.000000
       vegin%ratecp(1,4) =        1.000000
       vegin%cplant(2,4) =    12000.000000
       vegin%ratecp(2,4) =        0.030000
       vegin%cplant(3,4) =     1029.000000
       vegin%ratecp(3,4) =        0.140000
       vegin%a1gs(4) =        9.000000
       vegin%d0gs(4) =     1500.000000
       vegin%alpha(4) =        0.200000
       vegin%convex(4) =        0.700000
       vegin%cfrd(4) =        0.015000
       vegin%gswmin(4) =        0.010000
       vegin%conkc0(4) =        0.000302
       vegin%conko0(4) =        0.256000
       vegin%ekc(4) =    59430.000000
       vegin%eko(4) =    36000.000000
       vegin%g0(4) =        0.000000
       vegin%g1(4) =        4.447321
       vegin%zr(4) =        2.000000
       vegin%clitt(4) =       13.000000

       !PFT: shrub
       !=========================================================
       vegin%canst1(5) =        0.100000
       vegin%length(5) =        0.100000
       vegin%width(5) =        0.005000
       vegin%vcmax(5) =        0.000040
       vegin%ejmax(5) =        0.000000
       vegin%hc(5) =        0.600000
       vegin%xfang(5) =        0.010000
       vegin%rp20(5) =        1.000000
       vegin%rpcoef(5) =        0.083200
       vegin%rs20(5) =        1.000000
       vegin%wai(5) =        0.000000
       vegin%rootbeta(5) =        0.964000
       vegin%shelrb(5) =        2.000000
       vegin%vegcf(5) =        5.000000
       vegin%frac4(5) =        0.000000
       vegin%xalbnir(5) =        1.000000
       vegin%extkn(5) =        0.001000
       vegin%tminvj(5) =      -15.000000
       vegin%tmaxvj(5) =      -10.000000
       vegin%vbeta(5) =        4.000000
       vegin%froot(1,5) =        0.200000
       vegin%froot(2,5) =        0.200000
       vegin%froot(3,5) =        0.200000
       vegin%froot(4,5) =        0.200000
       vegin%froot(5,5) =        0.200000
       vegin%froot(6,5) =        0.200000
       vegin%refl(1,5) =        0.090000
       vegin%taul(1,5) =        0.090000
       vegin%refl(2,5) =        0.300000
       vegin%taul(2,5) =        0.300000
       vegin%refl(3,5) =        0.010000
       vegin%taul(3,5) =        0.010000
       vegin%csoil(1,5) =      100.000000
       vegin%ratecs(1,5) =        2.000000
       vegin%csoil(2,5) =      250.000000
       vegin%ratecs(2,5) =        0.500000
       vegin%cplant(1,5) =      159.000000
       vegin%ratecp(1,5) =        1.000000
       vegin%cplant(2,5) =     5000.000000
       vegin%ratecp(2,5) =        0.030000
       vegin%cplant(3,5) =      500.000000
       vegin%ratecp(3,5) =        0.140000
       vegin%a1gs(5) =        9.000000
       vegin%d0gs(5) =     1500.000000
       vegin%alpha(5) =        0.200000
       vegin%convex(5) =        0.700000
       vegin%cfrd(5) =        0.015000
       vegin%gswmin(5) =        0.010000
       vegin%conkc0(5) =        0.000302
       vegin%conko0(5) =        0.256000
       vegin%ekc(5) =    59430.000000
       vegin%eko(5) =    36000.000000
       vegin%g0(5) =        0.000000
       vegin%g1(5) =        4.694803
       vegin%zr(5) =        2.500000
       vegin%clitt(5) =        2.000000

       !PFT: C3
       !=========================================================
       vegin%canst1(6) =        0.100000
       vegin%length(6) =        0.300000
       vegin%width(6) =        0.010000
       vegin%vcmax(6) =        0.000060
       vegin%ejmax(6) =        0.000000
       vegin%hc(6) =        0.567000
       vegin%xfang(6) =       -0.300000
       vegin%rp20(6) =        1.500000
       vegin%rpcoef(6) =        0.083200
       vegin%rs20(6) =        1.000000
       vegin%wai(6) =        0.000000
       vegin%rootbeta(6) =        0.943000
       vegin%shelrb(6) =        2.000000
       vegin%vegcf(6) =        7.000000
       vegin%frac4(6) =        0.000000
       vegin%xalbnir(6) =        1.000000
       vegin%extkn(6) =        0.001000
       vegin%tminvj(6) =      -15.000000
       vegin%tmaxvj(6) =      -10.000000
       vegin%vbeta(6) =        4.000000
       vegin%froot(1,6) =        0.150000
       vegin%froot(2,6) =        0.150000
       vegin%froot(3,6) =        0.150000
       vegin%froot(4,6) =        0.150000
       vegin%froot(5,6) =        0.150000
       vegin%froot(6,6) =        0.150000
       vegin%refl(1,6) =        0.110000
       vegin%taul(1,6) =        0.110000
       vegin%refl(2,6) =        0.340000
       vegin%taul(2,6) =        0.340000
       vegin%refl(3,6) =        0.010000
       vegin%taul(3,6) =        0.010000
       vegin%csoil(1,6) =      275.000000
       vegin%ratecs(1,6) =        2.000000
       vegin%csoil(2,6) =      314.000000
       vegin%ratecs(2,6) =        0.500000
       vegin%cplant(1,6) =      250.000000
       vegin%ratecp(1,6) =        1.000000
       vegin%cplant(2,6) =        0.000000
       vegin%ratecp(2,6) =        0.030000
       vegin%cplant(3,6) =      500.000000
       vegin%ratecp(3,6) =        0.140000
       vegin%a1gs(6) =        9.000000
       vegin%d0gs(6) =     1500.000000
       vegin%alpha(6) =        0.200000
       vegin%convex(6) =        0.700000
       vegin%cfrd(6) =        0.015000
       vegin%gswmin(6) =        0.010000
       vegin%conkc0(6) =        0.000302
       vegin%conko0(6) =        0.256000
       vegin%ekc(6) =    59430.000000
       vegin%eko(6) =    36000.000000
       vegin%g0(6) =        0.000000
       vegin%g1(6) =        5.248500
       vegin%zr(6) =        0.500000    !1.5 in Haverd et al. (2016)
       vegin%clitt(6) =        2.000000

       !PFT: C4
       !=========================================================
       vegin%canst1(7) =        0.100000
       vegin%length(7) =        0.300000
       vegin%width(7) =        0.010000
       vegin%vcmax(7) =        0.000010
       vegin%ejmax(7) =        0.000000
       vegin%hc(7) =        0.567000
       vegin%xfang(7) =       -0.300000
       vegin%rp20(7) =        2.800000
       vegin%rpcoef(7) =        0.083200
       vegin%rs20(7) =        1.000000
       vegin%wai(7) =        0.000000
       vegin%rootbeta(7) =        0.943000
       vegin%shelrb(7) =        2.000000
       vegin%vegcf(7) =        7.000000
       vegin%frac4(7) =        1.000000
       vegin%xalbnir(7) =        1.000000
       vegin%extkn(7) =        0.001000
       vegin%tminvj(7) =      -15.000000
       vegin%tmaxvj(7) =      -10.000000
       vegin%vbeta(7) =        4.000000
       vegin%froot(1,7) =        0.000000
       vegin%froot(2,7) =        0.000000
       vegin%froot(3,7) =        0.000000
       vegin%froot(4,7) =        0.000000
       vegin%froot(5,7) =        0.000000
       vegin%froot(6,7) =        0.000000
       vegin%refl(1,7) =        0.110000
       vegin%taul(1,7) =        0.110000
       vegin%refl(2,7) =        0.340000
       vegin%taul(2,7) =        0.340000
       vegin%refl(3,7) =        0.010000
       vegin%taul(3,7) =        0.010000
       vegin%csoil(1,7) =      275.000000
       vegin%ratecs(1,7) =        2.000000
       vegin%csoil(2,7) =      314.000000
       vegin%ratecs(2,7) =        0.500000
       vegin%cplant(1,7) =      250.000000
       vegin%ratecp(1,7) =        1.000000
       vegin%cplant(2,7) =        0.000000
       vegin%ratecp(2,7) =        0.030000
       vegin%cplant(3,7) =      500.000000
       vegin%ratecp(3,7) =        0.140000
       vegin%a1gs(7) =        4.000000
       vegin%d0gs(7) =     1500.000000
       vegin%alpha(7) =        0.050000
       vegin%convex(7) =        0.800000
       vegin%cfrd(7) =        0.025000
       vegin%gswmin(7) =        0.040000
       vegin%conkc0(7) =        0.000302
       vegin%conko0(7) =        0.256000
       vegin%ekc(7) =    59430.000000
       vegin%eko(7) =    36000.000000
       vegin%g0(7) =        0.000000
       vegin%g1(7) =        1.616178
       vegin%zr(7) =        0.500000    !2.4 in Haverd et al. (2016)
       vegin%clitt(7) =        0.300000

       !PFT: Tundra
       !=========================================================
       vegin%canst1(8) =        0.100000
       vegin%length(8) =        0.300000
       vegin%width(8) =        0.010000
       vegin%vcmax(8) =        0.000040
       vegin%ejmax(8) =        0.000000
       vegin%hc(8) =        0.567000
       vegin%xfang(8) =       -0.300000
       vegin%rp20(8) =        2.500000
       vegin%rpcoef(8) =        0.083200
       vegin%rs20(8) =        1.000000
       vegin%wai(8) =        0.000000
       vegin%rootbeta(8) =        0.943000
       vegin%shelrb(8) =        2.000000
       vegin%vegcf(8) =        5.000000
       vegin%frac4(8) =        0.000000
       vegin%xalbnir(8) =        1.000000
       vegin%extkn(8) =        0.001000
       vegin%tminvj(8) =      -15.000000
       vegin%tmaxvj(8) =      -10.000000
       vegin%vbeta(8) =        4.000000
       vegin%froot(1,8) =        0.000000
       vegin%froot(2,8) =        0.000000
       vegin%froot(3,8) =        0.000000
       vegin%froot(4,8) =        0.000000
       vegin%froot(5,8) =        0.000000
       vegin%froot(6,8) =        0.000000
       vegin%refl(1,8) =        0.075000
       vegin%taul(1,8) =        0.075000
       vegin%refl(2,8) =        0.320000
       vegin%taul(2,8) =        0.320000
       vegin%refl(3,8) =        0.010000
       vegin%taul(3,8) =        0.010000
       vegin%csoil(1,8) =      275.000000
       vegin%ratecs(1,8) =        2.000000
       vegin%csoil(2,8) =      314.000000
       vegin%ratecs(2,8) =        0.500000
       vegin%cplant(1,8) =      250.000000
       vegin%ratecp(1,8) =        1.000000
       vegin%cplant(2,8) =        0.000000
       vegin%ratecp(2,8) =        0.030000
       vegin%cplant(3,8) =      500.000000
       vegin%ratecp(3,8) =        0.140000
       vegin%a1gs(8) =        9.000000
       vegin%d0gs(8) =     1500.000000
       vegin%alpha(8) =        0.200000
       vegin%convex(8) =        0.700000
       vegin%cfrd(8) =        0.015000
       vegin%gswmin(8) =        0.010000
       vegin%conkc0(8) =        0.000302
       vegin%conko0(8) =        0.256000
       vegin%ekc(8) =    59430.000000
       vegin%eko(8) =    36000.000000
       vegin%g0(8) =        0.000000
       vegin%g1(8) =        2.222156
       vegin%zr(8) =        0.500000
       vegin%clitt(8) =        0.300000

       !PFT: C3
       !=========================================================
       vegin%canst1(9) =        0.100000
       vegin%length(9) =        0.300000
       vegin%width(9) =        0.010000
       vegin%vcmax(9) =        0.000080
       vegin%ejmax(9) =        0.000000
       vegin%hc(9) =        0.550000
       vegin%xfang(9) =       -0.300000
       vegin%rp20(9) =        1.500000
       vegin%rpcoef(9) =        0.083200
       vegin%rs20(9) =        1.000000
       vegin%wai(9) =        0.000000
       vegin%rootbeta(9) =        0.961000
       vegin%shelrb(9) =        2.000000
       vegin%vegcf(9) =        7.000000
       vegin%frac4(9) =        0.000000
       vegin%xalbnir(9) =        1.000000
       vegin%extkn(9) =        0.001000
       vegin%tminvj(9) =      -15.000000
       vegin%tmaxvj(9) =      -10.000000
       vegin%vbeta(9) =        2.000000
       vegin%froot(1,9) =        0.000000
       vegin%froot(2,9) =        0.000000
       vegin%froot(3,9) =        0.000000
       vegin%froot(4,9) =        0.000000
       vegin%froot(5,9) =        0.000000
       vegin%froot(6,9) =        0.000000
       vegin%refl(1,9) =        0.110000
       vegin%taul(1,9) =        0.110000
       vegin%refl(2,9) =        0.340000
       vegin%taul(2,9) =        0.340000
       vegin%refl(3,9) =        0.010000
       vegin%taul(3,9) =        0.010000
       vegin%csoil(1,9) =      149.000000
       vegin%ratecs(1,9) =        2.000000
       vegin%csoil(2,9) =      300.000000
       vegin%ratecs(2,9) =        0.500000
       vegin%cplant(1,9) =      150.000000
       vegin%ratecp(1,9) =        1.000000
       vegin%cplant(2,9) =        0.000000
       vegin%ratecp(2,9) =        0.030000
       vegin%cplant(3,9) =      607.000000
       vegin%ratecp(3,9) =        0.140000
       vegin%a1gs(9) =        9.000000
       vegin%d0gs(9) =     1500.000000
       vegin%alpha(9) =        0.200000
       vegin%convex(9) =        0.700000
       vegin%cfrd(9) =        0.015000
       vegin%gswmin(9) =        0.010000
       vegin%conkc0(9) =        0.000302
       vegin%conko0(9) =        0.256000
       vegin%ekc(9) =    59430.000000
       vegin%eko(9) =    36000.000000
       vegin%g0(9) =        0.000000
       vegin%g1(9) =        5.789377
       vegin%zr(9) =        0.500000    !1.5 in Haverd et al. (2016)
       vegin%clitt(9) =        0.000000

       !PFT: C4
       !=========================================================
       vegin%canst1(10) =        0.100000
       vegin%length(10) =        0.300000
       vegin%width(10) =        0.010000
       vegin%vcmax(10) =        0.000080
       vegin%ejmax(10) =        0.000000
       vegin%hc(10) =        0.550000
       vegin%xfang(10) =       -0.300000
       vegin%rp20(10) =        1.000000
       vegin%rpcoef(10) =        0.083200
       vegin%rs20(10) =        1.000000
       vegin%wai(10) =        0.000000
       vegin%rootbeta(10) =        0.961000
       vegin%shelrb(10) =        2.000000
       vegin%vegcf(10) =        1.000000
       vegin%frac4(10) =        1.000000
       vegin%xalbnir(10) =        1.000000
       vegin%extkn(10) =        0.001000
       vegin%tminvj(10) =      -15.000000
       vegin%tmaxvj(10) =      -10.000000
       vegin%vbeta(10) =        2.000000
       vegin%froot( 1,10) =        0.000000
       vegin%froot( 2,10) =        0.000000
       vegin%froot( 3,10) =        0.000000
       vegin%froot( 4,10) =        0.000000
       vegin%froot( 5,10) =        0.000000
       vegin%froot( 6,10) =        0.000000
       vegin%refl( 1,10) =        0.110000
       vegin%taul( 1,10) =        0.110000
       vegin%refl( 2,10) =        0.340000
       vegin%taul( 2,10) =        0.340000
       vegin%refl( 3,10) =        0.010000
       vegin%taul( 3,10) =        0.010000
       vegin%csoil( 1,10) =      149.000000
       vegin%ratecs( 1,10) =        2.000000
       vegin%csoil( 2,10) =      300.000000
       vegin%ratecs( 2,10) =        0.500000
       vegin%cplant( 1,10) =      150.000000
       vegin%ratecp( 1,10) =        1.000000
       vegin%cplant( 2,10) =        0.000000
       vegin%ratecp( 2,10) =        0.030000
       vegin%cplant( 3,10) =      607.000000
       vegin%ratecp( 3,10) =        0.140000
       vegin%a1gs(10) =        4.000000
       vegin%d0gs(10) =     1500.000000
       vegin%alpha(10) =        0.050000
       vegin%convex(10) =        0.800000
       vegin%cfrd(10) =        0.025000
       vegin%gswmin(10) =        0.040000
       vegin%conkc0(10) =        0.000302
       vegin%conko0(10) =        0.256000
       vegin%ekc(10) =    59430.000000
       vegin%eko(10) =    36000.000000
       vegin%g0(10) =        0.000000
       vegin%g1(10) =        1.616178
       vegin%zr(10) =        0.500000    !1.5 in Haverd et al. (2016)
       vegin%clitt(10) =        0.000000

       !PFT: wetland
       !=========================================================
       vegin%canst1(11) =        0.100000
       vegin%length(11) =        0.300000
       vegin%width(11) =        0.010000
       vegin%vcmax(11) =        0.000060
       vegin%ejmax(11) =        0.000000
       vegin%hc(11) =        0.567000
       vegin%xfang(11) =       -0.300000
       vegin%rp20(11) =        1.500000
       vegin%rpcoef(11) =        0.083200
       vegin%rs20(11) =        1.000000
       vegin%wai(11) =        0.000000
       vegin%rootbeta(11) =        0.943000
       vegin%shelrb(11) =        2.000000
       vegin%vegcf(11) =        7.000000
       vegin%frac4(11) =        0.000000
       vegin%xalbnir(11) =        1.000000
       vegin%extkn(11) =        0.001000
       vegin%tminvj(11) =      -15.000000
       vegin%tmaxvj(11) =      -10.000000
       vegin%vbeta(11) =        4.000000
       vegin%froot( 1,11) =        0.000000
       vegin%froot( 2,11) =        0.000000
       vegin%froot( 3,11) =        0.000000
       vegin%froot( 4,11) =        0.000000
       vegin%froot( 5,11) =        0.000000
       vegin%froot( 6,11) =        0.000000
       vegin%refl( 1,11) =        0.108000
       vegin%taul( 1,11) =        0.075000
       vegin%refl( 2,11) =        0.343000
       vegin%taul( 2,11) =        0.146000
       vegin%refl( 3,11) =        0.010000
       vegin%taul( 3,11) =        0.010000
       vegin%csoil( 1,11) =      275.000000
       vegin%ratecs( 1,11) =        2.000000
       vegin%csoil( 2,11) =      314.000000
       vegin%ratecs( 2,11) =        0.500000
       vegin%cplant( 1,11) =      250.000000
       vegin%ratecp( 1,11) =        1.000000
       vegin%cplant( 2,11) =        0.000000
       vegin%ratecp( 2,11) =        0.030000
       vegin%cplant( 3,11) =      500.000000
       vegin%ratecp( 3,11) =        0.140000
       vegin%a1gs(11) =        9.000000
       vegin%d0gs(11) =     1500.000000
       vegin%alpha(11) =        0.200000
       vegin%convex(11) =        0.700000
       vegin%cfrd(11) =        0.015000
       vegin%gswmin(11) =        0.010000
       vegin%conkc0(11) =        0.000302
       vegin%conko0(11) =        0.256000
       vegin%ekc(11) =    59430.000000
       vegin%eko(11) =    36000.000000
       vegin%g0(11) =        0.000000
       vegin%g1(11) =        5.248500
       vegin%zr(11) =        1.800000
       vegin%clitt(11) =        2.000000

       !PFT: empty
       !=========================================================
       vegin%canst1(12) =        0.100000
       vegin%length(12) =        0.030000
       vegin%width(12) =        0.003000
       vegin%vcmax(12) =        0.000017
       vegin%ejmax(12) =        0.000000
       vegin%hc(12) =        0.200000
       vegin%xfang(12) =        0.100000
       vegin%rp20(12) =        1.000000
       vegin%rpcoef(12) =        0.083200
       vegin%rs20(12) =        0.000000
       vegin%wai(12) =        0.000000
       vegin%rootbeta(12) =        0.975000
       vegin%shelrb(12) =        2.000000
       vegin%vegcf(12) =        1.000000
       vegin%frac4(12) =        0.000000
       vegin%xalbnir(12) =        1.000000
       vegin%extkn(12) =        0.001000
       vegin%tminvj(12) =      -15.000000
       vegin%tmaxvj(12) =      -10.000000
       vegin%vbeta(12) =        4.000000
       vegin%froot( 1,12) =        0.000000
       vegin%froot( 2,12) =        0.000000
       vegin%froot( 3,12) =        0.000000
       vegin%froot( 4,12) =        0.000000
       vegin%froot( 5,12) =        0.000000
       vegin%froot( 6,12) =        0.000000
       vegin%refl( 1,12) =        0.055000
       vegin%taul( 1,12) =        0.023000
       vegin%refl( 2,12) =        0.190000
       vegin%taul( 2,12) =        0.198000
       vegin%refl( 3,12) =        0.010000
       vegin%taul( 3,12) =        0.010000
       vegin%csoil( 1,12) =        1.000000
       vegin%ratecs( 1,12) =        2.000000
       vegin%csoil( 2,12) =        1.000000
       vegin%ratecs( 2,12) =        0.500000
       vegin%cplant( 1,12) =        1.000000
       vegin%ratecp( 1,12) =        1.000000
       vegin%cplant( 2,12) =        0.000000
       vegin%ratecp( 2,12) =        0.030000
       vegin%cplant( 3,12) =        1.000000
       vegin%ratecp( 3,12) =        0.140000
       vegin%a1gs(12) =        9.000000
       vegin%d0gs(12) =     1500.000000
       vegin%alpha(12) =        0.200000
       vegin%convex(12) =        0.700000
       vegin%cfrd(12) =        0.015000
       vegin%gswmin(12) =        0.010000
       vegin%conkc0(12) =        0.000302
       vegin%conko0(12) =        0.256000
       vegin%ekc(12) =    59430.000000
       vegin%eko(12) =    36000.000000
       vegin%g0(12) =        0.000000
       vegin%g1(12) =        5.248500
       vegin%zr(12) =        3.100000
       vegin%clitt(12) =        2.000000

       !PFT: empty
       !=========================================================
       vegin%canst1(13) =        0.100000
       vegin%length(13) =        0.242000
       vegin%width(13) =        0.015000
       vegin%vcmax(13) =        0.000001
       vegin%ejmax(13) =        0.000000
       vegin%hc(13) =        6.017000
       vegin%xfang(13) =        0.000000
       vegin%rp20(13) =        1.000000
       vegin%rpcoef(13) =        0.083200
       vegin%rs20(13) =        1.000000
       vegin%wai(13) =        0.000000
       vegin%rootbeta(13) =        0.961000
       vegin%shelrb(13) =        2.000000
       vegin%vegcf(13) =        1.000000
       vegin%frac4(13) =        0.000000
       vegin%xalbnir(13) =        1.000000
       vegin%extkn(13) =        0.001000
       vegin%tminvj(13) =      -15.000000
       vegin%tmaxvj(13) =      -10.000000
       vegin%vbeta(13) =        2.000000
       vegin%froot( 1,13) =        0.000000
       vegin%froot( 2,13) =        0.000000
       vegin%froot( 3,13) =        0.000000
       vegin%froot( 4,13) =        0.000000
       vegin%froot( 5,13) =        0.000000
       vegin%froot( 6,13) =        0.000000
       vegin%refl( 1,13) =        0.091000
       vegin%taul( 1,13) =        0.059000
       vegin%refl( 2,13) =        0.310000
       vegin%taul( 2,13) =        0.163000
       vegin%refl( 3,13) =        0.010000
       vegin%taul( 3,13) =        0.010000
       vegin%csoil( 1,13) =        0.100000
       vegin%ratecs( 1,13) =        2.000000
       vegin%csoil( 2,13) =        0.100000
       vegin%ratecs( 2,13) =        0.500000
       vegin%cplant( 1,13) =        0.100000
       vegin%ratecp( 1,13) =        1.000000
       vegin%cplant( 2,13) =        0.000000
       vegin%ratecp( 2,13) =        0.030000
       vegin%cplant( 3,13) =        0.100000
       vegin%ratecp( 3,13) =        0.140000
       vegin%a1gs(13) =        9.000000
       vegin%d0gs(13) =     1500.000000
       vegin%alpha(13) =        0.200000
       vegin%convex(13) =        0.700000
       vegin%cfrd(13) =        0.015000
       vegin%gswmin(13) =        0.010000
       vegin%conkc0(13) =        0.000302
       vegin%conko0(13) =        0.256000
       vegin%ekc(13) =    59430.000000
       vegin%eko(13) =    36000.000000
       vegin%g0(13) =        0.000000
       vegin%g1(13) =        0.000000
       vegin%zr(13) =        3.000000
       vegin%clitt(13) =        0.000000

       !PFT: barren
       !=========================================================
       vegin%canst1(14) =        0.100000
       vegin%length(14) =        0.030000
       vegin%width(14) =        0.001000
       vegin%vcmax(14) =        0.000017
       vegin%ejmax(14) =        0.000000
       vegin%hc(14) =        0.200000
       vegin%xfang(14) =        0.000000
       vegin%rp20(14) =        1.000000
       vegin%rpcoef(14) =        0.083200
       vegin%rs20(14) =        0.000000
       vegin%wai(14) =        0.000000
       vegin%rootbeta(14) =        0.961000
       vegin%shelrb(14) =        2.000000
       vegin%vegcf(14) =        1.000000
       vegin%frac4(14) =        0.000000
       vegin%xalbnir(14) =        1.000000
       vegin%extkn(14) =        0.001000
       vegin%tminvj(14) =      -15.000000
       vegin%tmaxvj(14) =      -10.000000
       vegin%vbeta(14) =        4.000000
       vegin%froot( 1,14) =        0.000000
       vegin%froot( 2,14) =        0.000000
       vegin%froot( 3,14) =        0.000000
       vegin%froot( 4,14) =        0.000000
       vegin%froot( 5,14) =        0.000000
       vegin%froot( 6,14) =        0.000000
       vegin%refl( 1,14) =        0.238000
       vegin%taul( 1,14) =        0.039000
       vegin%refl( 2,14) =        0.457000
       vegin%taul( 2,14) =        0.189000
       vegin%refl( 3,14) =        0.010000
       vegin%taul( 3,14) =        0.010000
       vegin%csoil( 1,14) =        1.000000
       vegin%ratecs( 1,14) =        2.000000
       vegin%csoil( 2,14) =        1.000000
       vegin%ratecs( 2,14) =        0.500000
       vegin%cplant( 1,14) =        0.000000
       vegin%ratecp( 1,14) =        1.000000
       vegin%cplant( 2,14) =        0.000000
       vegin%ratecp( 2,14) =        0.030000
       vegin%cplant( 3,14) =        0.000000
       vegin%ratecp( 3,14) =        0.140000
       vegin%a1gs(14) =        9.000000
       vegin%d0gs(14) =     1500.000000
       vegin%alpha(14) =        0.200000
       vegin%convex(14) =        0.700000
       vegin%cfrd(14) =        0.015000
       vegin%gswmin(14) =        0.010000
       vegin%conkc0(14) =        0.000302
       vegin%conko0(14) =        0.256000
       vegin%ekc(14) =    59430.000000
       vegin%eko(14) =    36000.000000
       vegin%g0(14) =        0.000000
       vegin%g1(14) =        5.248500
       vegin%zr(14) =        1.000000
       vegin%clitt(14) =        0.000000

       !PFT: urban
       !=========================================================
       vegin%canst1(15) =        0.100000
       vegin%length(15) =        0.030000
       vegin%width(15) =        0.001000
       vegin%vcmax(15) =        0.000017
       vegin%ejmax(15) =        0.000000
       vegin%hc(15) =        0.200000
       vegin%xfang(15) =        0.000000
       vegin%rp20(15) =        1.000000
       vegin%rpcoef(15) =        0.083200
       vegin%rs20(15) =        0.000000
       vegin%wai(15) =        0.000000
       vegin%rootbeta(15) =        0.961000
       vegin%shelrb(15) =        2.000000
       vegin%vegcf(15) =        1.000000
       vegin%frac4(15) =        0.000000
       vegin%xalbnir(15) =        1.000000
       vegin%extkn(15) =        0.001000
       vegin%tminvj(15) =      -15.000000
       vegin%tmaxvj(15) =      -10.000000
       vegin%vbeta(15) =        4.000000
       vegin%froot( 1,15) =        0.000000
       vegin%froot( 2,15) =        0.000000
       vegin%froot( 3,15) =        0.000000
       vegin%froot( 4,15) =        0.000000
       vegin%froot( 5,15) =        0.000000
       vegin%froot( 6,15) =        0.000000
       vegin%refl( 1,15) =        0.143000
       vegin%taul( 1,15) =        0.023000
       vegin%refl( 2,15) =        0.275000
       vegin%taul( 2,15) =        0.113000
       vegin%refl( 3,15) =        0.010000
       vegin%taul( 3,15) =        0.010000
       vegin%csoil( 1,15) =        1.000000
       vegin%ratecs( 1,15) =        2.000000
       vegin%csoil( 2,15) =        1.000000
       vegin%ratecs( 2,15) =        0.500000
       vegin%cplant( 1,15) =        1.000000
       vegin%ratecp( 1,15) =        1.000000
       vegin%cplant( 2,15) =        0.000000
       vegin%ratecp( 2,15) =        0.030000
       vegin%cplant( 3,15) =        1.000000
       vegin%ratecp( 3,15) =        0.140000
       vegin%a1gs(15) =        9.000000
       vegin%d0gs(15) =     1500.000000
       vegin%alpha(15) =        0.200000
       vegin%convex(15) =        0.700000
       vegin%cfrd(15) =        0.015000
       vegin%gswmin(15) =        0.010000
       vegin%conkc0(15) =        0.000302
       vegin%conko0(15) =        0.256000
       vegin%ekc(15) =    59430.000000
       vegin%eko(15) =    36000.000000
       vegin%g0(15) =        0.000000
       vegin%g1(15) =        5.248500
       vegin%zr(15) =        1.000000
       vegin%clitt(15) =        0.000000

       !PFT: lakes
       !=========================================================
       vegin%canst1(16) =        0.100000
       vegin%length(16) =        0.030000
       vegin%width(16) =        0.001000
       vegin%vcmax(16) =        0.000017
       vegin%ejmax(16) =        0.000000
       vegin%hc(16) =        0.200000
       vegin%xfang(16) =        0.000000
       vegin%rp20(16) =        1.000000
       vegin%rpcoef(16) =        0.083200
       vegin%rs20(16) =        0.000000
       vegin%wai(16) =        0.000000
       vegin%rootbeta(16) =        0.961000
       vegin%shelrb(16) =        2.000000
       vegin%vegcf(16) =        1.000000
       vegin%frac4(16) =        0.000000
       vegin%xalbnir(16) =        1.000000
       vegin%extkn(16) =        0.001000
       vegin%tminvj(16) =      -15.000000
       vegin%tmaxvj(16) =      -10.000000
       vegin%vbeta(16) =        4.000000
       vegin%froot( 1,16) =        0.000000
       vegin%froot( 2,16) =        0.000000
       vegin%froot( 3,16) =        0.000000
       vegin%froot( 4,16) =        0.000000
       vegin%froot( 5,16) =        0.000000
       vegin%froot( 6,16) =        0.000000
       vegin%refl( 1,16) =        0.143000
       vegin%taul( 1,16) =        0.023000
       vegin%refl( 2,16) =        0.275000
       vegin%taul( 2,16) =        0.113000
       vegin%refl( 3,16) =        0.010000
       vegin%taul( 3,16) =        0.010000
       vegin%csoil( 1,16) =        1.000000
       vegin%ratecs( 1,16) =        2.000000
       vegin%csoil( 2,16) =        1.000000
       vegin%ratecs( 2,16) =        0.500000
       vegin%cplant( 1,16) =        1.000000
       vegin%ratecp( 1,16) =        1.000000
       vegin%cplant( 2,16) =        0.000000
       vegin%ratecp( 2,16) =        0.030000
       vegin%cplant( 3,16) =        1.000000
       vegin%ratecp( 3,16) =        0.140000
       vegin%a1gs(16) =        9.000000
       vegin%d0gs(16) =     1500.000000
       vegin%alpha(16) =        0.200000
       vegin%convex(16) =        0.700000
       vegin%cfrd(16) =        0.015000
       vegin%gswmin(16) =        0.010000
       vegin%conkc0(16) =        0.000302
       vegin%conko0(16) =        0.256000
       vegin%ekc(16) =    59430.000000
       vegin%eko(16) =    36000.000000
       vegin%g0(16) =        0.000000
       vegin%g1(16) =        5.248500
       vegin%zr(16) =        1.000000
       vegin%clitt(16) =        0.000000

       !PFT: ice
       !=========================================================
       vegin%canst1(17) =        0.100000
       vegin%length(17) =        0.030000
       vegin%width(17) =        0.001000
       vegin%vcmax(17) =        0.000017
       vegin%ejmax(17) =        0.000000
       vegin%hc(17) =        0.200000
       vegin%xfang(17) =        0.000000
       vegin%rp20(17) =        1.000000
       vegin%rpcoef(17) =        0.083200
       vegin%rs20(17) =        0.000000
       vegin%wai(17) =        0.000000
       vegin%rootbeta(17) =        0.961000
       vegin%shelrb(17) =        2.000000
       vegin%vegcf(17) =        1.000000
       vegin%frac4(17) =        0.000000
       vegin%xalbnir(17) =        1.000000
       vegin%extkn(17) =        0.001000
       vegin%tminvj(17) =      -15.000000
       vegin%tmaxvj(17) =      -10.000000
       vegin%vbeta(17) =        4.000000
       vegin%froot( 1,17) =        0.000000
       vegin%froot( 2,17) =        0.000000
       vegin%froot( 3,17) =        0.000000
       vegin%froot( 4,17) =        0.000000
       vegin%froot( 5,17) =        0.000000
       vegin%froot( 6,17) =        0.000000
       vegin%refl( 1,17) =        0.159000
       vegin%taul( 1,17) =        0.026000
       vegin%refl( 2,17) =        0.305000
       vegin%taul( 2,17) =        0.113000
       vegin%refl( 3,17) =        0.010000
       vegin%taul( 3,17) =        0.010000
       vegin%csoil( 1,17) =        1.000000
       vegin%ratecs( 1,17) =        2.000000
       vegin%csoil( 2,17) =        1.000000
       vegin%ratecs( 2,17) =        0.500000
       vegin%cplant( 1,17) =        0.000000
       vegin%ratecp( 1,17) =        1.000000
       vegin%cplant( 2,17) =        0.000000
       vegin%ratecp( 2,17) =        0.030000
       vegin%cplant( 3,17) =        0.000000
       vegin%ratecp( 3,17) =        0.140000
       vegin%a1gs(17) =        9.000000
       vegin%d0gs(17) =     1500.000000
       vegin%alpha(17) =        0.200000
       vegin%convex(17) =        0.700000
       vegin%cfrd(17) =        0.015000
       vegin%gswmin(17) =        0.010000
       vegin%conkc0(17) =        0.000302
       vegin%conko0(17) =        0.256000
       vegin%ekc(17) =    59430.000000
       vegin%eko(17) =    36000.000000
       vegin%g0(17) =        0.000000
       vegin%g1(17) =        5.248500
       vegin%zr(17) =        1.000000
       vegin%clitt(17) =        0.000000

       !PFT: New PFT, Eucalyptus blakelyi, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(18) = 0.75 !1.5
       vegin%Kcrit(18) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(18) = 5.030576014
       vegin%c_plant(18) = 3.355198155
       vegin%gmin(18) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(18) = 0.000087
       vegin%ejmax(18) = 0.000145

       vegin%canst1(18) =        0.100000
       vegin%length(18) =        0.100000
       vegin%width(18) =        0.050000
       !vegin%vcmax(18) =        0.000055
       !vegin%ejmax(18) =        0.000000
       vegin%hc(18) =       35.000000
       vegin%xfang(18) =        0.100000
       vegin%rp20(18) =        0.600000
       vegin%rpcoef(18) =        0.083200
       vegin%rs20(18) =        1.000000
       vegin%wai(18) =        1.000000
       vegin%rootbeta(18) =        0.962000
       vegin%shelrb(18) =        2.000000
       vegin%vegcf(18) =       14.000000
       vegin%frac4(18) =        0.000000
       vegin%xalbnir(18) =        1.000000
       vegin%extkn(18) =        0.001000
       vegin%tminvj(18) =      -15.000000
       vegin%tmaxvj(18) =      -10.000000
       vegin%vbeta(18) =        2.000000
       vegin%froot(1,18) =        0.200000
       vegin%froot(2,18) =        0.200000
       vegin%froot(3,18) =        0.200000
       vegin%froot(4,18) =        0.200000
       vegin%froot(5,18) =        0.200000
       vegin%froot(6,18) =        0.200000
       vegin%refl(1,18) =        0.076000
       vegin%taul(1,18) =        0.050000
       vegin%refl(2,18) =        0.350000
       vegin%taul(2,18) =        0.250000
       vegin%refl(3,18) =        0.010000
       vegin%taul(3,18) =        0.010000
       vegin%csoil(1,18) =      303.000000
       vegin%ratecs(1,18) =        2.000000
       vegin%csoil(2,18) =      606.000000
       vegin%ratecs(2,18) =        0.500000
       vegin%cplant(1,18) =      300.000000
       vegin%ratecp(1,18) =        1.000000
       vegin%cplant(2,18) =    16833.000000
       vegin%ratecp(2,18) =        0.030000
       vegin%cplant(3,18) =     1443.000000
       vegin%ratecp(3,18) =        0.140000
       vegin%a1gs(18) =        9.000000
       vegin%d0gs(18) =     1500.000000
       vegin%alpha(18) =        0.200000
       vegin%convex(18) =        0.700000
       vegin%cfrd(18) =        0.015000
       vegin%gswmin(18) =        0.010000
       vegin%conkc0(18) =        0.000302
       vegin%conko0(18) =        0.256000
       vegin%ekc(18) =    59430.000000
       vegin%eko(18) =    36000.000000
       vegin%g0(18) =        1E-09!0.000000
       !vegin%g1(18) =        4.114762

       vegin%zr(18) =        3.000000
       vegin%clitt(18) =        6.000000

       !PFT: New PFT, Eucalyptus camaldulensis, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(19) = 0.75 !1.5
       vegin%Kcrit(19) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(19) = 4.102123226
       vegin%c_plant(19) = 4.347350048
       vegin%gmin(19) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(19) = 0.000112
       vegin%ejmax(19) = 0.000187

       vegin%canst1(19) =        0.100000
       vegin%length(19) =        0.100000
       vegin%width(19) =        0.050000
       !vegin%vcmax(19) =        0.000055
       !vegin%ejmax(19) =        0.000000
       vegin%hc(19) =       35.000000
       vegin%xfang(19) =        0.100000
       vegin%rp20(19) =        0.600000
       vegin%rpcoef(19) =        0.083200
       vegin%rs20(19) =        1.000000
       vegin%wai(19) =        1.000000
       vegin%rootbeta(19) =        0.962000
       vegin%shelrb(19) =        2.000000
       vegin%vegcf(19) =       14.000000
       vegin%frac4(19) =        0.000000
       vegin%xalbnir(19) =        1.000000
       vegin%extkn(19) =        0.001000
       vegin%tminvj(19) =      -15.000000
       vegin%tmaxvj(19) =      -10.000000
       vegin%vbeta(19) =        2.000000
       vegin%froot(1,19) =        0.200000
       vegin%froot(2,19) =        0.200000
       vegin%froot(3,19) =        0.200000
       vegin%froot(4,19) =        0.200000
       vegin%froot(5,19) =        0.200000
       vegin%froot(6,19) =        0.200000
       vegin%refl(1,19) =        0.076000
       vegin%taul(1,19) =        0.050000
       vegin%refl(2,19) =        0.350000
       vegin%taul(2,19) =        0.250000
       vegin%refl(3,19) =        0.010000
       vegin%taul(3,19) =        0.010000
       vegin%csoil(1,19) =      303.000000
       vegin%ratecs(1,19) =        2.000000
       vegin%csoil(2,19) =      606.000000
       vegin%ratecs(2,19) =        0.500000
       vegin%cplant(1,19) =      300.000000
       vegin%ratecp(1,19) =        1.000000
       vegin%cplant(2,19) =    16833.000000
       vegin%ratecp(2,19) =        0.030000
       vegin%cplant(3,19) =     1443.000000
       vegin%ratecp(3,19) =        0.140000
       vegin%a1gs(19) =        9.000000
       vegin%d0gs(19) =     1500.000000
       vegin%alpha(19) =        0.200000
       vegin%convex(19) =        0.700000
       vegin%cfrd(19) =        0.015000
       vegin%gswmin(19) =        0.010000
       vegin%conkc0(19) =        0.000302
       vegin%conko0(19) =        0.256000
       vegin%ekc(19) =    59430.000000
       vegin%eko(19) =    36000.000000
       vegin%g0(19) =        1E-09!0.000000
       !vegin%g1(19) =        4.114762

       vegin%zr(19) =        3.000000
       vegin%clitt(19) =        6.000000

       !PFT: New PFT, Eucalyptus crebra, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(20) = 0.75 !1.5
       vegin%Kcrit(20) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(20) = 5.52013939
       vegin%c_plant(20) = 3.075600898
       vegin%gmin(20) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(20) = 0.000087
       vegin%ejmax(20) = 0.000145

       vegin%canst1(20) =        0.100000
       vegin%length(20) =        0.100000
       vegin%width(20) =        0.050000
       !vegin%vcmax(20) =        0.000055
       !vegin%ejmax(20) =        0.000000
       vegin%hc(20) =       35.000000
       vegin%xfang(20) =        0.100000
       vegin%rp20(20) =        0.600000
       vegin%rpcoef(20) =        0.083200
       vegin%rs20(20) =        1.000000
       vegin%wai(20) =        1.000000
       vegin%rootbeta(20) =        0.962000
       vegin%shelrb(20) =        2.000000
       vegin%vegcf(20) =       14.000000
       vegin%frac4(20) =        0.000000
       vegin%xalbnir(20) =        1.000000
       vegin%extkn(20) =        0.001000
       vegin%tminvj(20) =      -15.000000
       vegin%tmaxvj(20) =      -10.000000
       vegin%vbeta(20) =        2.000000
       vegin%froot(1,20) =        0.200000
       vegin%froot(2,20) =        0.200000
       vegin%froot(3,20) =        0.200000
       vegin%froot(4,20) =        0.200000
       vegin%froot(5,20) =        0.200000
       vegin%froot(6,20) =        0.200000
       vegin%refl(1,20) =        0.076000
       vegin%taul(1,20) =        0.050000
       vegin%refl(2,20) =        0.350000
       vegin%taul(2,20) =        0.250000
       vegin%refl(3,20) =        0.010000
       vegin%taul(3,20) =        0.010000
       vegin%csoil(1,20) =      303.000000
       vegin%ratecs(1,20) =        2.000000
       vegin%csoil(2,20) =      606.000000
       vegin%ratecs(2,20) =        0.500000
       vegin%cplant(1,20) =      300.000000
       vegin%ratecp(1,20) =        1.000000
       vegin%cplant(2,20) =    16833.000000
       vegin%ratecp(2,20) =        0.030000
       vegin%cplant(3,20) =     1443.000000
       vegin%ratecp(3,20) =        0.140000
       vegin%a1gs(20) =        9.000000
       vegin%d0gs(20) =     1500.000000
       vegin%alpha(20) =        0.200000
       vegin%convex(20) =        0.700000
       vegin%cfrd(20) =        0.015000
       vegin%gswmin(20) =        0.010000
       vegin%conkc0(20) =        0.000302
       vegin%conko0(20) =        0.256000
       vegin%ekc(20) =    59430.000000
       vegin%eko(20) =    36000.000000
       vegin%g0(20) =        1E-09!0.000000
       !vegin%g1(20) =        4.114762

       vegin%zr(20) =        3.000000
       vegin%clitt(20) =        6.000000

       !PFT: New PFT, Eucalyptus dunnii, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(21) = 0.75 !1.5
       vegin%Kcrit(21) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(21) = 5.557394056
       vegin%c_plant(21) = 3.059620633
       vegin%gmin(21) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(21) = 0.000087
       vegin%ejmax(21) = 0.000145

       vegin%canst1(21) =        0.100000
       vegin%length(21) =        0.100000
       vegin%width(21) =        0.050000
       !vegin%vcmax(21) =        0.000055
       !vegin%ejmax(21) =        0.000000
       vegin%hc(21) =       35.000000
       vegin%xfang(21) =        0.100000
       vegin%rp20(21) =        0.600000
       vegin%rpcoef(21) =        0.083200
       vegin%rs20(21) =        1.000000
       vegin%wai(21) =        1.000000
       vegin%rootbeta(21) =        0.962000
       vegin%shelrb(21) =        2.000000
       vegin%vegcf(21) =       14.000000
       vegin%frac4(21) =        0.000000
       vegin%xalbnir(21) =        1.000000
       vegin%extkn(21) =        0.001000
       vegin%tminvj(21) =      -15.000000
       vegin%tmaxvj(21) =      -10.000000
       vegin%vbeta(21) =        2.000000
       vegin%froot(1,21) =        0.200000
       vegin%froot(2,21) =        0.200000
       vegin%froot(3,21) =        0.200000
       vegin%froot(4,21) =        0.200000
       vegin%froot(5,21) =        0.200000
       vegin%froot(6,21) =        0.200000
       vegin%refl(1,21) =        0.076000
       vegin%taul(1,21) =        0.050000
       vegin%refl(2,21) =        0.350000
       vegin%taul(2,21) =        0.250000
       vegin%refl(3,21) =        0.010000
       vegin%taul(3,21) =        0.010000
       vegin%csoil(1,21) =      303.000000
       vegin%ratecs(1,21) =        2.000000
       vegin%csoil(2,21) =      606.000000
       vegin%ratecs(2,21) =        0.500000
       vegin%cplant(1,21) =      300.000000
       vegin%ratecp(1,21) =        1.000000
       vegin%cplant(2,21) =    16833.000000
       vegin%ratecp(2,21) =        0.030000
       vegin%cplant(3,21) =     1443.000000
       vegin%ratecp(3,21) =        0.140000
       vegin%a1gs(21) =        9.000000
       vegin%d0gs(21) =     1500.000000
       vegin%alpha(21) =        0.200000
       vegin%convex(21) =        0.700000
       vegin%cfrd(21) =        0.015000
       vegin%gswmin(21) =        0.010000
       vegin%conkc0(21) =        0.000302
       vegin%conko0(21) =        0.256000
       vegin%ekc(21) =    59430.000000
       vegin%eko(21) =    36000.000000
       vegin%g0(21) =        1E-09!0.000000
       !vegin%g1(21) =        4.114762

       vegin%zr(21) =        3.000000
       vegin%clitt(21) =        6.000000

       !PFT: New PFT, Eucalyptus globulus, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(22) = 0.75 !1.5
       vegin%Kcrit(22) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(22) = 2.550186588
       vegin%c_plant(22) = 8.298063366
       vegin%gmin(22) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(22) = 0.000086
       vegin%ejmax(22) = 0.000144

       vegin%canst1(22) =        0.100000
       vegin%length(22) =        0.100000
       vegin%width(22) =        0.050000
       !vegin%vcmax(22) =        0.000055
       !vegin%ejmax(22) =        0.000000
       vegin%hc(22) =       35.000000
       vegin%xfang(22) =        0.100000
       vegin%rp20(22) =        0.600000
       vegin%rpcoef(22) =        0.083200
       vegin%rs20(22) =        1.000000
       vegin%wai(22) =        1.000000
       vegin%rootbeta(22) =        0.962000
       vegin%shelrb(22) =        2.000000
       vegin%vegcf(22) =       14.000000
       vegin%frac4(22) =        0.000000
       vegin%xalbnir(22) =        1.000000
       vegin%extkn(22) =        0.001000
       vegin%tminvj(22) =      -15.000000
       vegin%tmaxvj(22) =      -10.000000
       vegin%vbeta(22) =        2.000000
       vegin%froot(1,22) =        0.200000
       vegin%froot(2,22) =        0.200000
       vegin%froot(3,22) =        0.200000
       vegin%froot(4,22) =        0.200000
       vegin%froot(5,22) =        0.200000
       vegin%froot(6,22) =        0.200000
       vegin%refl(1,22) =        0.076000
       vegin%taul(1,22) =        0.050000
       vegin%refl(2,22) =        0.350000
       vegin%taul(2,22) =        0.250000
       vegin%refl(3,22) =        0.010000
       vegin%taul(3,22) =        0.010000
       vegin%csoil(1,22) =      303.000000
       vegin%ratecs(1,22) =        2.000000
       vegin%csoil(2,22) =      606.000000
       vegin%ratecs(2,22) =        0.500000
       vegin%cplant(1,22) =      300.000000
       vegin%ratecp(1,22) =        1.000000
       vegin%cplant(2,22) =    16833.000000
       vegin%ratecp(2,22) =        0.030000
       vegin%cplant(3,22) =     1443.000000
       vegin%ratecp(3,22) =        0.140000
       vegin%a1gs(22) =        9.000000
       vegin%d0gs(22) =     1500.000000
       vegin%alpha(22) =        0.200000
       vegin%convex(22) =        0.700000
       vegin%cfrd(22) =        0.015000
       vegin%gswmin(22) =        0.010000
       vegin%conkc0(22) =        0.000302
       vegin%conko0(22) =        0.256000
       vegin%ekc(22) =    59430.000000
       vegin%eko(22) =    36000.000000
       vegin%g0(22) =        1E-09!0.000000
       !vegin%g1(22) =        4.114762

       vegin%zr(22) =        3.000000
       vegin%clitt(22) =        6.000000


       !PFT: New PFT, Eucalyptus grandis, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(23) = 0.75 !1.5
       vegin%Kcrit(23) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(23) = 3.579776939
       vegin%c_plant(23) = 5.286522315
       vegin%gmin(23) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(23) = 0.000094
       vegin%ejmax(23) = 0.000157

       vegin%canst1(23) =        0.100000
       vegin%length(23) =        0.100000
       vegin%width(23) =        0.050000
       !vegin%vcmax(23) =        0.000055
       !vegin%ejmax(23) =        0.000000
       vegin%hc(23) =       35.000000
       vegin%xfang(23) =        0.100000
       vegin%rp20(23) =        0.600000
       vegin%rpcoef(23) =        0.083200
       vegin%rs20(23) =        1.000000
       vegin%wai(23) =        1.000000
       vegin%rootbeta(23) =        0.962000
       vegin%shelrb(23) =        2.000000
       vegin%vegcf(23) =       14.000000
       vegin%frac4(23) =        0.000000
       vegin%xalbnir(23) =        1.000000
       vegin%extkn(23) =        0.001000
       vegin%tminvj(23) =      -15.000000
       vegin%tmaxvj(23) =      -10.000000
       vegin%vbeta(23) =        2.000000
       vegin%froot(1,23) =        0.200000
       vegin%froot(2,23) =        0.200000
       vegin%froot(3,23) =        0.200000
       vegin%froot(4,23) =        0.200000
       vegin%froot(5,23) =        0.200000
       vegin%froot(6,23) =        0.200000
       vegin%refl(1,23) =        0.076000
       vegin%taul(1,23) =        0.050000
       vegin%refl(2,23) =        0.350000
       vegin%taul(2,23) =        0.250000
       vegin%refl(3,23) =        0.010000
       vegin%taul(3,23) =        0.010000
       vegin%csoil(1,23) =      303.000000
       vegin%ratecs(1,23) =        2.000000
       vegin%csoil(2,23) =      606.000000
       vegin%ratecs(2,23) =        0.500000
       vegin%cplant(1,23) =      300.000000
       vegin%ratecp(1,23) =        1.000000
       vegin%cplant(2,23) =    16833.000000
       vegin%ratecp(2,23) =        0.030000
       vegin%cplant(3,23) =     1443.000000
       vegin%ratecp(3,23) =        0.140000
       vegin%a1gs(23) =        9.000000
       vegin%d0gs(23) =     1500.000000
       vegin%alpha(23) =        0.200000
       vegin%convex(23) =        0.700000
       vegin%cfrd(23) =        0.015000
       vegin%gswmin(23) =        0.010000
       vegin%conkc0(23) =        0.000302
       vegin%conko0(23) =        0.256000
       vegin%ekc(23) =    59430.000000
       vegin%eko(23) =    36000.000000
       vegin%g0(23) =        1E-09!0.000000
       !vegin%g1(23) =        4.114762

       vegin%zr(23) =        3.000000
       vegin%clitt(23) =        6.000000

       !PFT: New PFT, Eucalyptus largiflorens, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(24) = 0.75 !1.5
       vegin%Kcrit(24) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(24) = 8.28282906
       vegin%c_plant(24) = 3.251978832
       vegin%gmin(24) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(24) = 0.000087
       vegin%ejmax(24) = 0.000145

       vegin%canst1(24) =        0.100000
       vegin%length(24) =        0.100000
       vegin%width(24) =        0.050000
       !vegin%vcmax(24) =        0.000055
       !vegin%ejmax(24) =        0.000000
       vegin%hc(24) =       35.000000
       vegin%xfang(24) =        0.100000
       vegin%rp20(24) =        0.600000
       vegin%rpcoef(24) =        0.083200
       vegin%rs20(24) =        1.000000
       vegin%wai(24) =        1.000000
       vegin%rootbeta(24) =        0.962000
       vegin%shelrb(24) =        2.000000
       vegin%vegcf(24) =       14.000000
       vegin%frac4(24) =        0.000000
       vegin%xalbnir(24) =        1.000000
       vegin%extkn(24) =        0.001000
       vegin%tminvj(24) =      -15.000000
       vegin%tmaxvj(24) =      -10.000000
       vegin%vbeta(24) =        2.000000
       vegin%froot(1,24) =        0.200000
       vegin%froot(2,24) =        0.200000
       vegin%froot(3,24) =        0.200000
       vegin%froot(4,24) =        0.200000
       vegin%froot(5,24) =        0.200000
       vegin%froot(6,24) =        0.200000
       vegin%refl(1,24) =        0.076000
       vegin%taul(1,24) =        0.050000
       vegin%refl(2,24) =        0.350000
       vegin%taul(2,24) =        0.250000
       vegin%refl(3,24) =        0.010000
       vegin%taul(3,24) =        0.010000
       vegin%csoil(1,24) =      303.000000
       vegin%ratecs(1,24) =        2.000000
       vegin%csoil(2,24) =      606.000000
       vegin%ratecs(2,24) =        0.500000
       vegin%cplant(1,24) =      300.000000
       vegin%ratecp(1,24) =        1.000000
       vegin%cplant(2,24) =    16833.000000
       vegin%ratecp(2,24) =        0.030000
       vegin%cplant(3,24) =     1443.000000
       vegin%ratecp(3,24) =        0.140000
       vegin%a1gs(24) =        9.000000
       vegin%d0gs(24) =     1500.000000
       vegin%alpha(24) =        0.200000
       vegin%convex(24) =        0.700000
       vegin%cfrd(24) =        0.015000
       vegin%gswmin(24) =        0.010000
       vegin%conkc0(24) =        0.000302
       vegin%conko0(24) =        0.256000
       vegin%ekc(24) =    59430.000000
       vegin%eko(24) =    36000.000000
       vegin%g0(24) =        1E-09!0.000000
       !vegin%g1(24) =        4.114762

       vegin%zr(24) =        3.000000
       vegin%clitt(24) =        6.000000

       !PFT: New PFT, Eucalyptus macrorhyncha, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(25) = 0.75 !1.5
       vegin%Kcrit(25) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(25) = 4.400036199
       vegin%c_plant(25) = 3.948576733
       vegin%gmin(25) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(25) = 0.000087
       vegin%ejmax(25) = 0.000145

       vegin%canst1(25) =        0.100000
       vegin%length(25) =        0.100000
       vegin%width(25) =        0.050000
       !vegin%vcmax(25) =        0.000055
       !vegin%ejmax(25) =        0.000000
       vegin%hc(25) =       35.000000
       vegin%xfang(25) =        0.100000
       vegin%rp20(25) =        0.600000
       vegin%rpcoef(25) =        0.083200
       vegin%rs20(25) =        1.000000
       vegin%wai(25) =        1.000000
       vegin%rootbeta(25) =        0.962000
       vegin%shelrb(25) =        2.000000
       vegin%vegcf(25) =       14.000000
       vegin%frac4(25) =        0.000000
       vegin%xalbnir(25) =        1.000000
       vegin%extkn(25) =        0.001000
       vegin%tminvj(25) =      -15.000000
       vegin%tmaxvj(25) =      -10.000000
       vegin%vbeta(25) =        2.000000
       vegin%froot(1,25) =        0.200000
       vegin%froot(2,25) =        0.200000
       vegin%froot(3,25) =        0.200000
       vegin%froot(4,25) =        0.200000
       vegin%froot(5,25) =        0.200000
       vegin%froot(6,25) =        0.200000
       vegin%refl(1,25) =        0.076000
       vegin%taul(1,25) =        0.050000
       vegin%refl(2,25) =        0.350000
       vegin%taul(2,25) =        0.250000
       vegin%refl(3,25) =        0.010000
       vegin%taul(3,25) =        0.010000
       vegin%csoil(1,25) =      303.000000
       vegin%ratecs(1,25) =        2.000000
       vegin%csoil(2,25) =      606.000000
       vegin%ratecs(2,25) =        0.500000
       vegin%cplant(1,25) =      300.000000
       vegin%ratecp(1,25) =        1.000000
       vegin%cplant(2,25) =    16833.000000
       vegin%ratecp(2,25) =        0.030000
       vegin%cplant(3,25) =     1443.000000
       vegin%ratecp(3,25) =        0.140000
       vegin%a1gs(25) =        9.000000
       vegin%d0gs(25) =     1500.000000
       vegin%alpha(25) =        0.200000
       vegin%convex(25) =        0.700000
       vegin%cfrd(25) =        0.015000
       vegin%gswmin(25) =        0.010000
       vegin%conkc0(25) =        0.000302
       vegin%conko0(25) =        0.256000
       vegin%ekc(25) =    59430.000000
       vegin%eko(25) =    36000.000000
       vegin%g0(25) =        1E-09!0.000000
       !vegin%g1(25) =        4.114762

       vegin%zr(25) =        3.000000
       vegin%clitt(25) =        6.000000








       !PFT: New PFT, Eucalyptus melliodora, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(26) = 0.75 !1.5
       vegin%Kcrit(26) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(26) = 5.668675883
       vegin%c_plant(26) = 3.015931136
       vegin%gmin(26) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(26) = 0.000087
       vegin%ejmax(26) = 0.000145

       vegin%canst1(26) =        0.100000
       vegin%length(26) =        0.100000
       vegin%width(26) =        0.050000
       !vegin%vcmax(26) =        0.000055
       !vegin%ejmax(26) =        0.000000
       vegin%hc(26) =       35.000000
       vegin%xfang(26) =        0.100000
       vegin%rp20(26) =        0.600000
       vegin%rpcoef(26) =        0.083200
       vegin%rs20(26) =        1.000000
       vegin%wai(26) =        1.000000
       vegin%rootbeta(26) =        0.962000
       vegin%shelrb(26) =        2.000000
       vegin%vegcf(26) =       14.000000
       vegin%frac4(26) =        0.000000
       vegin%xalbnir(26) =        1.000000
       vegin%extkn(26) =        0.001000
       vegin%tminvj(26) =      -15.000000
       vegin%tmaxvj(26) =      -10.000000
       vegin%vbeta(26) =        2.000000
       vegin%froot(1,26) =        0.200000
       vegin%froot(2,26) =        0.200000
       vegin%froot(3,26) =        0.200000
       vegin%froot(4,26) =        0.200000
       vegin%froot(5,26) =        0.200000
       vegin%froot(6,26) =        0.200000
       vegin%refl(1,26) =        0.076000
       vegin%taul(1,26) =        0.050000
       vegin%refl(2,26) =        0.350000
       vegin%taul(2,26) =        0.250000
       vegin%refl(3,26) =        0.010000
       vegin%taul(3,26) =        0.010000
       vegin%csoil(1,26) =      303.000000
       vegin%ratecs(1,26) =        2.000000
       vegin%csoil(2,26) =      606.000000
       vegin%ratecs(2,26) =        0.500000
       vegin%cplant(1,26) =      300.000000
       vegin%ratecp(1,26) =        1.000000
       vegin%cplant(2,26) =    16833.000000
       vegin%ratecp(2,26) =        0.030000
       vegin%cplant(3,26) =     1443.000000
       vegin%ratecp(3,26) =        0.140000
       vegin%a1gs(26) =        9.000000
       vegin%d0gs(26) =     1500.000000
       vegin%alpha(26) =        0.200000
       vegin%convex(26) =        0.700000
       vegin%cfrd(26) =        0.015000
       vegin%gswmin(26) =        0.010000
       vegin%conkc0(26) =        0.000302
       vegin%conko0(26) =        0.256000
       vegin%ekc(26) =    59430.000000
       vegin%eko(26) =    36000.000000
       vegin%g0(26) =        1E-09!0.000000
       !vegin%g1(26) =        4.114762

       vegin%zr(26) =        3.000000
       vegin%clitt(26) =        6.000000

       !PFT: New PFT, Eucalyptus obliqua, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(27) = 0.75 !1.5
       vegin%Kcrit(27) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(27) = 2.725198383
       vegin%c_plant(27) = 7.667745303
       vegin%gmin(27) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(27) = 0.000087
       vegin%ejmax(27) = 0.000145

       vegin%canst1(27) =        0.100000
       vegin%length(27) =        0.100000
       vegin%width(27) =        0.050000
       !vegin%vcmax(27) =        0.000055
       !vegin%ejmax(27) =        0.000000
       vegin%hc(27) =       35.000000
       vegin%xfang(27) =        0.100000
       vegin%rp20(27) =        0.600000
       vegin%rpcoef(27) =        0.083200
       vegin%rs20(27) =        1.000000
       vegin%wai(27) =        1.000000
       vegin%rootbeta(27) =        0.962000
       vegin%shelrb(27) =        2.000000
       vegin%vegcf(27) =       14.000000
       vegin%frac4(27) =        0.000000
       vegin%xalbnir(27) =        1.000000
       vegin%extkn(27) =        0.001000
       vegin%tminvj(27) =      -15.000000
       vegin%tmaxvj(27) =      -10.000000
       vegin%vbeta(27) =        2.000000
       vegin%froot(1,27) =        0.200000
       vegin%froot(2,27) =        0.200000
       vegin%froot(3,27) =        0.200000
       vegin%froot(4,27) =        0.200000
       vegin%froot(5,27) =        0.200000
       vegin%froot(6,27) =        0.200000
       vegin%refl(1,27) =        0.076000
       vegin%taul(1,27) =        0.050000
       vegin%refl(2,27) =        0.350000
       vegin%taul(2,27) =        0.250000
       vegin%refl(3,27) =        0.010000
       vegin%taul(3,27) =        0.010000
       vegin%csoil(1,27) =      303.000000
       vegin%ratecs(1,27) =        2.000000
       vegin%csoil(2,27) =      606.000000
       vegin%ratecs(2,27) =        0.500000
       vegin%cplant(1,27) =      300.000000
       vegin%ratecp(1,27) =        1.000000
       vegin%cplant(2,27) =    16833.000000
       vegin%ratecp(2,27) =        0.030000
       vegin%cplant(3,27) =     1443.000000
       vegin%ratecp(3,27) =        0.140000
       vegin%a1gs(27) =        9.000000
       vegin%d0gs(27) =     1500.000000
       vegin%alpha(27) =        0.200000
       vegin%convex(27) =        0.700000
       vegin%cfrd(27) =        0.015000
       vegin%gswmin(27) =        0.010000
       vegin%conkc0(27) =        0.000302
       vegin%conko0(27) =        0.256000
       vegin%ekc(27) =    59430.000000
       vegin%eko(27) =    36000.000000
       vegin%g0(27) =        1E-09!0.000000
       !vegin%g1(27) =        4.114762

       vegin%zr(27) =        3.000000
       vegin%clitt(27) =        6.000000

       !PFT: New PFT, Eucalyptus populnea, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(28) = 0.75 !1.5
       vegin%Kcrit(28) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(28) = 6.410296129
       vegin%c_plant(28) = 2.8629035
       vegin%gmin(28) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(28) = 0.000087
       vegin%ejmax(28) = 0.000145

       vegin%canst1(28) =        0.100000
       vegin%length(28) =        0.100000
       vegin%width(28) =        0.050000
       !vegin%vcmax(28) =        0.000055
       !vegin%ejmax(28) =        0.000000
       vegin%hc(28) =       35.000000
       vegin%xfang(28) =        0.100000
       vegin%rp20(28) =        0.600000
       vegin%rpcoef(28) =        0.083200
       vegin%rs20(28) =        1.000000
       vegin%wai(28) =        1.000000
       vegin%rootbeta(28) =        0.962000
       vegin%shelrb(28) =        2.000000
       vegin%vegcf(28) =       14.000000
       vegin%frac4(28) =        0.000000
       vegin%xalbnir(28) =        1.000000
       vegin%extkn(28) =        0.001000
       vegin%tminvj(28) =      -15.000000
       vegin%tmaxvj(28) =      -10.000000
       vegin%vbeta(28) =        2.000000
       vegin%froot(1,28) =        0.200000
       vegin%froot(2,28) =        0.200000
       vegin%froot(3,28) =        0.200000
       vegin%froot(4,28) =        0.200000
       vegin%froot(5,28) =        0.200000
       vegin%froot(6,28) =        0.200000
       vegin%refl(1,28) =        0.076000
       vegin%taul(1,28) =        0.050000
       vegin%refl(2,28) =        0.350000
       vegin%taul(2,28) =        0.250000
       vegin%refl(3,28) =        0.010000
       vegin%taul(3,28) =        0.010000
       vegin%csoil(1,28) =      303.000000
       vegin%ratecs(1,28) =        2.000000
       vegin%csoil(2,28) =      606.000000
       vegin%ratecs(2,28) =        0.500000
       vegin%cplant(1,28) =      300.000000
       vegin%ratecp(1,28) =        1.000000
       vegin%cplant(2,28) =    16833.000000
       vegin%ratecp(2,28) =        0.030000
       vegin%cplant(3,28) =     1443.000000
       vegin%ratecp(3,28) =        0.140000
       vegin%a1gs(28) =        9.000000
       vegin%d0gs(28) =     1500.000000
       vegin%alpha(28) =        0.200000
       vegin%convex(28) =        0.700000
       vegin%cfrd(28) =        0.015000
       vegin%gswmin(28) =        0.010000
       vegin%conkc0(28) =        0.000302
       vegin%conko0(28) =        0.256000
       vegin%ekc(28) =    59430.000000
       vegin%eko(28) =    36000.000000
       vegin%g0(28) =        1E-09!0.000000
       !vegin%g1(28) =        4.114762

       vegin%zr(28) =        3.000000
       vegin%clitt(28) =        6.000000

       !PFT: New PFT, Eucalyptus saligna, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(29) = 0.75 !1.5
       vegin%Kcrit(29) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(29) = 3.651439207
       vegin%c_plant(29) = 5.13712216
       vegin%gmin(29) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(29) = 0.000077
       vegin%ejmax(29) = 0.000129

       vegin%canst1(29) =        0.100000
       vegin%length(29) =        0.100000
       vegin%width(29) =        0.050000
       !vegin%vcmax(29) =        0.000055
       !vegin%ejmax(29) =        0.000000
       vegin%hc(29) =       35.000000
       vegin%xfang(29) =        0.100000
       vegin%rp20(29) =        0.600000
       vegin%rpcoef(29) =        0.083200
       vegin%rs20(29) =        1.000000
       vegin%wai(29) =        1.000000
       vegin%rootbeta(29) =        0.962000
       vegin%shelrb(29) =        2.000000
       vegin%vegcf(29) =       14.000000
       vegin%frac4(29) =        0.000000
       vegin%xalbnir(29) =        1.000000
       vegin%extkn(29) =        0.001000
       vegin%tminvj(29) =      -15.000000
       vegin%tmaxvj(29) =      -10.000000
       vegin%vbeta(29) =        2.000000
       vegin%froot(1,29) =        0.200000
       vegin%froot(2,29) =        0.200000
       vegin%froot(3,29) =        0.200000
       vegin%froot(4,29) =        0.200000
       vegin%froot(5,29) =        0.200000
       vegin%froot(6,29) =        0.200000
       vegin%refl(1,29) =        0.076000
       vegin%taul(1,29) =        0.050000
       vegin%refl(2,29) =        0.350000
       vegin%taul(2,29) =        0.250000
       vegin%refl(3,29) =        0.010000
       vegin%taul(3,29) =        0.010000
       vegin%csoil(1,29) =      303.000000
       vegin%ratecs(1,29) =        2.000000
       vegin%csoil(2,29) =      606.000000
       vegin%ratecs(2,29) =        0.500000
       vegin%cplant(1,29) =      300.000000
       vegin%ratecp(1,29) =        1.000000
       vegin%cplant(2,29) =    16833.000000
       vegin%ratecp(2,29) =        0.030000
       vegin%cplant(3,29) =     1443.000000
       vegin%ratecp(3,29) =        0.140000
       vegin%a1gs(29) =        9.000000
       vegin%d0gs(29) =     1500.000000
       vegin%alpha(29) =        0.200000
       vegin%convex(29) =        0.700000
       vegin%cfrd(29) =        0.015000
       vegin%gswmin(29) =        0.010000
       vegin%conkc0(29) =        0.000302
       vegin%conko0(29) =        0.256000
       vegin%ekc(29) =    59430.000000
       vegin%eko(29) =    36000.000000
       vegin%g0(29) =        1E-09!0.000000
       !vegin%g1(29) =        4.114762

       vegin%zr(29) =        3.000000
       vegin%clitt(29) =        6.000000

       !PFT: New PFT, E. sideroxylon, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(30) = 0.75 !1.5
       vegin%Kcrit(30) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(30) = 4.538102944
       vegin%c_plant(30) = 3.791935013
       vegin%gmin(30) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(30) = 0.000100
       vegin%ejmax(30) = 0.000167

       vegin%canst1(30) =        0.100000
       vegin%length(30) =        0.100000
       vegin%width(30) =        0.050000
       !vegin%vcmax(30) =        0.000055
       !vegin%ejmax(30) =        0.000000
       vegin%hc(30) =       35.000000
       vegin%xfang(30) =        0.100000
       vegin%rp20(30) =        0.600000
       vegin%rpcoef(30) =        0.083200
       vegin%rs20(30) =        1.000000
       vegin%wai(30) =        1.000000
       vegin%rootbeta(30) =        0.962000
       vegin%shelrb(30) =        2.000000
       vegin%vegcf(30) =       14.000000
       vegin%frac4(30) =        0.000000
       vegin%xalbnir(30) =        1.000000
       vegin%extkn(30) =        0.001000
       vegin%tminvj(30) =      -15.000000
       vegin%tmaxvj(30) =      -10.000000
       vegin%vbeta(30) =        2.000000
       vegin%froot(1,30) =        0.200000
       vegin%froot(2,30) =        0.200000
       vegin%froot(3,30) =        0.200000
       vegin%froot(4,30) =        0.200000
       vegin%froot(5,30) =        0.200000
       vegin%froot(6,30) =        0.200000
       vegin%refl(1,30) =        0.076000
       vegin%taul(1,30) =        0.050000
       vegin%refl(2,30) =        0.350000
       vegin%taul(2,30) =        0.250000
       vegin%refl(3,30) =        0.010000
       vegin%taul(3,30) =        0.010000
       vegin%csoil(1,30) =      303.000000
       vegin%ratecs(1,30) =        2.000000
       vegin%csoil(2,30) =      606.000000
       vegin%ratecs(2,30) =        0.500000
       vegin%cplant(1,30) =      300.000000
       vegin%ratecp(1,30) =        1.000000
       vegin%cplant(2,30) =    16833.000000
       vegin%ratecp(2,30) =        0.030000
       vegin%cplant(3,30) =     1443.000000
       vegin%ratecp(3,30) =        0.140000
       vegin%a1gs(30) =        9.000000
       vegin%d0gs(30) =     1500.000000
       vegin%alpha(30) =        0.200000
       vegin%convex(30) =        0.700000
       vegin%cfrd(30) =        0.015000
       vegin%gswmin(30) =        0.010000
       vegin%conkc0(30) =        0.000302
       vegin%conko0(30) =        0.256000
       vegin%ekc(30) =    59430.000000
       vegin%eko(30) =    36000.000000
       vegin%g0(30) =        1E-09!0.000000
       !vegin%g1(30) =        4.114762

       vegin%zr(30) =        3.000000
       vegin%clitt(30) =        6.000000


       !PFT: New PFT, Eucalyptus tereticornis, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(31) = 0.75 !1.5
       vegin%Kcrit(31) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(31) = 4.362501154
       vegin%c_plant(31) = 3.994093418
       vegin%gmin(31) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(31) = 0.000087
       vegin%ejmax(31) = 0.000145

       vegin%canst1(31) =        0.100000
       vegin%length(31) =        0.100000
       vegin%width(31) =        0.050000
       !vegin%vcmax(31) =        0.000055
       !vegin%ejmax(31) =        0.000000
       vegin%hc(31) =       35.000000
       vegin%xfang(31) =        0.100000
       vegin%rp20(31) =        0.600000
       vegin%rpcoef(31) =        0.083200
       vegin%rs20(31) =        1.000000
       vegin%wai(31) =        1.000000
       vegin%rootbeta(31) =        0.962000
       vegin%shelrb(31) =        2.000000
       vegin%vegcf(31) =       14.000000
       vegin%frac4(31) =        0.000000
       vegin%xalbnir(31) =        1.000000
       vegin%extkn(31) =        0.001000
       vegin%tminvj(31) =      -15.000000
       vegin%tmaxvj(31) =      -10.000000
       vegin%vbeta(31) =        2.000000
       vegin%froot(1,31) =        0.200000
       vegin%froot(2,31) =        0.200000
       vegin%froot(3,31) =        0.200000
       vegin%froot(4,31) =        0.200000
       vegin%froot(5,31) =        0.200000
       vegin%froot(6,31) =        0.200000
       vegin%refl(1,31) =        0.076000
       vegin%taul(1,31) =        0.050000
       vegin%refl(2,31) =        0.350000
       vegin%taul(2,31) =        0.250000
       vegin%refl(3,31) =        0.010000
       vegin%taul(3,31) =        0.010000
       vegin%csoil(1,31) =      303.000000
       vegin%ratecs(1,31) =        2.000000
       vegin%csoil(2,31) =      606.000000
       vegin%ratecs(2,31) =        0.500000
       vegin%cplant(1,31) =      300.000000
       vegin%ratecp(1,31) =        1.000000
       vegin%cplant(2,31) =    16833.000000
       vegin%ratecp(2,31) =        0.030000
       vegin%cplant(3,31) =     1443.000000
       vegin%ratecp(3,31) =        0.140000
       vegin%a1gs(31) =        9.000000
       vegin%d0gs(31) =     1500.000000
       vegin%alpha(31) =        0.200000
       vegin%convex(31) =        0.700000
       vegin%cfrd(31) =        0.015000
       vegin%gswmin(31) =        0.010000
       vegin%conkc0(31) =        0.000302
       vegin%conko0(31) =        0.256000
       vegin%ekc(31) =    59430.000000
       vegin%eko(31) =    36000.000000
       vegin%g0(31) =        1E-09!0.000000
       !vegin%g1(31) =        4.114762

       vegin%zr(31) =        3.000000
       vegin%clitt(31) =        6.000000

       !PFT: New PFT, Eucalyptus viminalis, mgk576, 25/7/21
       !=========================================================

       vegin%Kmax(32) = 0.75 !1.5
       vegin%Kcrit(32) = 0.0375 !0.075 ! vegin%Kmax * 0.05
       vegin%b_plant(32) = 3.355422803
       vegin%c_plant(32) = 5.801519551
       vegin%gmin(32) = 0.720828 ! single sided, above was double-sided, mmol m-2 s-1
       vegin%vcmax(32) = 0.000078
       vegin%ejmax(32) = 0.000130

       vegin%canst1(32) =        0.100000
       vegin%length(32) =        0.100000
       vegin%width(32) =        0.050000
       !vegin%vcmax(32) =        0.000055
       !vegin%ejmax(32) =        0.000000
       vegin%hc(32) =       35.000000
       vegin%xfang(32) =        0.100000
       vegin%rp20(32) =        0.600000
       vegin%rpcoef(32) =        0.083200
       vegin%rs20(32) =        1.000000
       vegin%wai(32) =        1.000000
       vegin%rootbeta(32) =        0.962000
       vegin%shelrb(32) =        2.000000
       vegin%vegcf(32) =       14.000000
       vegin%frac4(32) =        0.000000
       vegin%xalbnir(32) =        1.000000
       vegin%extkn(32) =        0.001000
       vegin%tminvj(32) =      -15.000000
       vegin%tmaxvj(32) =      -10.000000
       vegin%vbeta(32) =        2.000000
       vegin%froot(1,32) =        0.200000
       vegin%froot(2,32) =        0.200000
       vegin%froot(3,32) =        0.200000
       vegin%froot(4,32) =        0.200000
       vegin%froot(5,32) =        0.200000
       vegin%froot(6,32) =        0.200000
       vegin%refl(1,32) =        0.076000
       vegin%taul(1,32) =        0.050000
       vegin%refl(2,32) =        0.350000
       vegin%taul(2,32) =        0.250000
       vegin%refl(3,32) =        0.010000
       vegin%taul(3,32) =        0.010000
       vegin%csoil(1,32) =      303.000000
       vegin%ratecs(1,32) =        2.000000
       vegin%csoil(2,32) =      606.000000
       vegin%ratecs(2,32) =        0.500000
       vegin%cplant(1,32) =      300.000000
       vegin%ratecp(1,32) =        1.000000
       vegin%cplant(2,32) =    16833.000000
       vegin%ratecp(2,32) =        0.030000
       vegin%cplant(3,32) =     1443.000000
       vegin%ratecp(3,32) =        0.140000
       vegin%a1gs(32) =        9.000000
       vegin%d0gs(32) =     1500.000000
       vegin%alpha(32) =        0.200000
       vegin%convex(32) =        0.700000
       vegin%cfrd(32) =        0.015000
       vegin%gswmin(32) =        0.010000
       vegin%conkc0(32) =        0.000302
       vegin%conko0(32) =        0.256000
       vegin%ekc(32) =    59430.000000
       vegin%eko(32) =    36000.000000
       vegin%g0(32) =        1E-09!0.000000
       !vegin%g1(32) =        4.114762

       vegin%zr(32) =        3.000000
       vegin%clitt(32) =        6.000000


    ENDIF

    first_call = .FALSE.

    ! new calculation dleaf since April 2012 (cable v1.8 did not use width)
    vegin%dleaf = SQRT(vegin%width * vegin%length)

  END SUBROUTINE cable_pft_params

END MODULE cable_pft_params_mod
