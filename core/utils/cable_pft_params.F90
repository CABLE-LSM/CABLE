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
          clitt

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


subroutine cable_pft_params()

   ! Gets parameter values for each vegetation type 
   USE cable_def_types_mod, ONLY : mvtype, ms, ncs, ncp, nrb 

   INTEGER :: a, jveg ! do loop counter
  logical, save :: first_call = .true.
   mvtype=17    

    ! Allocate memory for type-specific vegetation parameters:
  if( first_call ) then
  
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
         vegin%zr(mvtype), vegin%clitt(mvtype) )

 !PFT parameters: description and corresponding variable name in code. 
 !PFT parameters are assigned as TYPE vegin% but later used as veg%
 
 !PFT: evergreen_needleleaf                                                  
 !=========================================================
    vegin%canst1(1) =        0.100000
   vegin%length(1) =        0.055000
    vegin%width(1) =        0.001000
    vegin%vcmax(1) =        0.000040
    vegin%ejmax(1) =        0.000000
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
  endif

  first_call = .false.
      
   ! new calculation dleaf since April 2012 (cable v1.8 did not use width)
   vegin%dleaf = SQRT(vegin%width * vegin%length)
    
End subroutine cable_pft_params

END MODULE cable_pft_params_mod

