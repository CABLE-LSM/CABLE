!> \file mo_c13o2_photosynthesis.f90

!> \brief Photosynthetic isotope discrimination

!> \details Photosynthetic carbon isotope discrimination after Farquhar et al.
!> (Aust J Plant Physiol, 1982), Farquhar (Aust J Plant Physiol, 1983),
!> Evans et al. (Aust J Plant Physiol, 1986), and Lloyd & Farquhar (Oecologia, 1994).
!> It includes extensions similar to Wingate et al. (Plant Cell Environ, 2007)
!> so that discrimination can be calculated over the whole diurnal cycle.

!> \author Matthias Cuntz
!> \date Apr 2019

! License
! -------
! This file is part of the JAMS Fortran package, distributed under the MIT License.
!
! Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

MODULE mo_c13o2_photosynthesis

  use mo_kind, only: dp, i4

  implicit none

  private

  ! Public routines and variables
  
  public :: c13o2_discrimination        ! Photosynthetic 13C discrimination
  public :: c13o2_discrimination_simple ! Simple 13C discrimination = a+(b-a)*ci/ca
  public :: init_starch_pool            ! Initialise starch pool and its isotope ratio
  public :: init_sugar_pools            ! Initialise isotope ratios of leaf carbohydrate pools

  ! Transitory starch concentration in leaf [mol(C)/m2]
  real(dp), dimension(:), allocatable, public :: Vstarch
  ! Isotopic composition if transitory starch
  real(dp), dimension(:), allocatable, public :: Rstarch
  ! Isotopic composition if leaf sucrose
  real(dp), dimension(:), allocatable, public :: Rsucrose
  ! Isotopic composition if pool used for photorespiration
  real(dp), dimension(:), allocatable, public :: Rphoto

  
  ! Private parameters
  
  ! Mesophyll conductance factor for C3
  ! ISOLSM takes 8000 but for C4. John Evans suggests that it should be ca. 3000 for C3.
  ! gm for C4 is about twice that of C3 (Evans & v.Caemmerer 1996).
  ! John Evans is critical about temp dependence of gm from Bernacchi et al. (2002)
  ! and suggests: if any temp depence, take the one of Vcmax.
  ! -> C3: gm(CO2) =   meso_cond_fac*Vcmax [mol(CO2) m-2 s-1]
  !    C4: gm(CO2) = 2*meso_cond_fac*Vcmax
  real(dp), parameter :: meso_cond_fac = 4000.0_dp

  ! fraction of leaf respiration in mesophyll in C4, rest in bundle sheats: Rdm = frdm*Rd
  real(dp), parameter :: frdm = 0.5_dp    ! von Caemmerer (2000), Table 4.1

  ! mesophyll to bundle sheat conductance for CO2 in C4, von Caemmerer (2000) Section 4.3.2.2
  real(dp), parameter :: gbsc = 2.e-3_dp 

  ! Leakage of bundle sheat to mesophyll in C4 plants during the day, von Caemmerer (2000) Section 4.3.3.2
  ! LPJ used 0.4
  real(dp), parameter :: Phi  = 0.2_dp ! von Caemmerer (2000), Fig. 4.11

  ! Fractionation of leaf respiration
  ! Ghashghaie et al. (2003) reviews around -6.
  ! Tcherkez et al. (2004) calculates -5.5
  ! Gessler et al. (2008) determines -2. for Ricinus
  real(dp), parameter :: eps_e = -6.0e-03_dp

  ! Fractionation of photo respiration
  ! Ghashghaie et al. (2003) reviews around +10.
  ! Tcherkez et al. (2004) calculates -9.2
  real(dp), parameter :: eps_f = 10.0e-03_dp

  ! Parameters for photosynthesis discrimination: Lloyd and Farquhar (1994)
  real(dp), parameter :: &
       eps_a    =  4.4e-03_dp, & ! diffusion frac.
       eps_al   =  0.7e-03_dp, & ! diffusion in liquid
       eps_b3   = 30.0e-03_dp, & ! fractionation during RUBISCO carboxylation
       eps_b_ci = 27.0e-03_dp, & ! effective carboxylation fractionation in 'simple' model
       beta     =  0.05_dp       ! fraction of PEP carboxylation in C3

  ! Sucrose content of leaf
  ! 100-500 mol(C6)/gDW guess from figures in Grimmer & Komor (1999) and Komor (2000) for Ricinus
  ! Conversion to mol(CO2)/m2(leaf) based on values from Ricinus experiment with Claudia Keitel
  ! -> mol(C6)/gDW * C/C6 * gDW/leaf * leaf/cm2(leaf) * cm2(leaf)/m2(leaf)
  real(dp), parameter :: Vsucrose = 300.e-6_dp * 6._dp * 2.2_dp * 1.0_dp/500._dp * 1.e4_dp

  ! Photorespiration pool
  ! This is some mixed sugar/enzyme pool.
  ! Because 90% of the sugars in leaves are sucrose, this can be a maximum 10%.
  ! This is definitely too large. It is little, though, compared to the oxygenation flux and hence mixes rapidly.
  ! Mathematically, we only need that the isotope ratio is not equal the assimilation.
  ! Hence, take 10% of sucrose pool.
  real(dp), parameter :: Vphoto = Vsucrose/10._dp

  ! Fractionation of starch synthesis
  ! Tcherkez et al. (2004) models about Rchloroplast/Rinput=1.006, i.e. eps_starch = -6.
  ! The measured equilibrium fractionation is -4.4 (according to Gerd Gleixner)
  real(dp), parameter :: eps_starch = -4.4e-03_dp

  ! Fractionation of sucrose
  ! Tcherkez et al. (2004) models about Rcytoplasm/Rinput=1.0019, i.e. eps_sucrose ~ -2.
  ! Also in this simple model, sucrose and starch are linked (ca.):
  !   A*RA = F_starch*(1-eps_starch)*RA + F_sucrose*(1-eps_sucrose)*RA
  ! so if F_starch = A/3 and F_sucrose=2A/3 then eps_sucrose = -eps_starch/2.
  ! This is the 'wrong' direction, according to Tcherkez et al. (2004).
  ! But the effect on sucrose is rather small so that we rather 'close the mass balance' and take -eps_starch/2.
  real(dp), parameter :: eps_sucrose = -eps_starch * 0.5_dp

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         c13o2_discrimination

  !     PURPOSE
  !>        \brief 13C leaf discrimination

  !>        \details 13C leaf discrimination D after Farquhar et al. (Aust J Plant Physiol, 1982),
  !>        Farquhar (Aust J Plant Physiol, 1983), Evans et al. (Aust J Plant Physiol, 1986), and
  !>        Lloyd & Farquhar (Oecologia, 1994), with the extension of Wingate et al.
  !>        (Plant Cell Environ, 2007) for leaf respiration and an analog extension for photorespiration.
  !>        This is for C3 during daytime: \n
  !>            D = k*ca/(k*(ca-G)-Rd) * (a + (b-a)*cc/ca + G/ca*(df-da-f) + Rd/(k*Ca)*(de-da-e)) \n
  !>        where \n
  !>            k = (An+Rd)/(ci-G) so that the carboxylation rate is Vc = k*ci and the oxygenation rate
  !>                is Vo = 2.*k*G, \n
  !>            G is the CO2 compensation point in absence of dark respiration (Gamma star), \n
  !>            Rd is mitochondrial respiration, \n
  !>            ca is the ambient mixing ratio outside the leaf boundary layer, \n
  !>            cc is the mixing ratio at the site of carboxylation (in chloroplast), \n
  !>            a is the effective diffusion fractionation from the ambient air outside the leaf boundary layer
  !>                into the site of carboxylation (interior of chloroplasts for C3, mesophyll for C4), \n
  !>            b is the effective fractionation of (RuBisCO and PEP-C) carboxylation, \n
  !>            e is the effective fractionation during mitochondrial respiration (-6. permil), \n
  !>            f is the effective fractionation during photorespiration (+10. permil), \n
  !>            da is the isotopic composition of ambient air outside the leaf boundary layer, \n
  !>            de is the isotopic composition of the substrate for mitochondrial respiration, and \n
  !>            df is the isotopic composition of the substrate for photorespiration. \n
  !>        A sugar pool (sucrose) with the size of 300 mol(C6)/gDW is taken as the substrate for
  !>        mitochondrial respiration. The isotopic composition of the sugar pool is renewed by assimilation. \n
  !>        The photorespiration pool is a mixture of sugars and enzymes. 90% of sugars in leaves
  !>        are sucrose so that the substrate pool is taken as 10% of the sugar pool. The pool is also
  !>        turned over isotopically by assimilation. \n
  !>        A starch pool is built up during day by 1/3 of assimilation, with an equilibrium fractionation of -4.4 permil.
  !>        Mitochondrial respiration respires this starch pool at night, with the average
  !>        isotopic composition of the pool at the beginning of the night, without fractionation. \n

  !     CALLING SEQUENCE
  !         call subroutine c13o2_discrimination(dt, isc3, Vcmax, GPP, Rd, Gammastar, &
  !                                              ca, ci, ga, gb, gs, Tl, Rair, &
  !                                              Vstarch, Rstarch, Rsucrose, Rphoto, &
  !                                              Disc, Ass13)

  !>        This is an elemental subroutine and can be called with scalars and arrays. <- ToDo

  !     PARAMETER
  !>        \param[in]    "real(dp) :: dt"         Time step [s]
  !>        \param[in]    "logical  :: isc3"       C3 mask with .true. for C3 and .false. for C4 plants
  !>        \param[in]    "real(dp) :: Vcmax"      Maximal carboxylation rate at leaf temperature Vcmax(Tl) [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: GPP"        Gross assimilation [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: Rd"         Leaf respiration [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: Gammstar"   CO2 compensation point in absence of dark respiration [mol(CO2)/mol(air)]
  !>        \param[in]    "real(dp) :: ca"         Ambient CO2 mixing ratio [mol(CO2)/mol(air)]
  !>        \param[in]    "real(dp) :: ci"         Stomatal CO2 mixing ratio [mol(CO2)/mol(air)]
  !>        \param[in]    "real(dp) :: gsc"        Stomatal conductance for CO2 [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: gac"        Aerodynamic conductance for CO2 [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: gbc"        Boundary layer conductance for CO2 [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: Tl"         Leaf temperature [K]
  !>        \param[in]    "real(dp) :: Rair"       Isotope ratio of ambient CO2 ([13C]/[12C])
  !>        \param[inout] "real(dp) :: Vstarch"    Transitory starch pool [mol(C)/m2]
  !>        \param[inout] "real(dp) :: Rstarch"    Isotope ratio of transitory starch ([13C]/[12C])
  !>        \param[inout] "real(dp) :: Rsucrose"   Isotope ratio of sucrose pool ([13C]/[12C])
  !>        \param[inout] "real(dp) :: Rphoto"     Isotope ratio of pool for photorespiration ([13C]/[12C])
  !>        \param[out]   "real(dp) :: Disc"       Leaf discrimination
  !>        \param[out]   "real(dp) :: Ass13"      13C net assimilation [mol(13CO2)/m^2/s]

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr-Nov 2019
  subroutine c13o2_discrimination( &
       ! -- Input
       dt, isc3, &
       ! Photosynthesis variables
       Vcmax, GPP, Rd, Gammastar, &
       ! CO2 concentrations
       ca, ci, &
       ! Conductances
       gac, gbc, gsc, &
       ! leaf temperature
       Tl, &
       ! Ambient isotope ratio
       Rair, &
       ! -- Inout
       ! Starch pool and isotope ratios of pools for respiration
       Vstarch, Rstarch, Rsucrose, Rphoto, &
       ! -- Output
       ! discrimination
       Disc, &
       ! 13CO2 flux
       Ass13)

    use mo_utils,     only: ne
    use mo_constants, only: onethird_dp, twothird_dp, T0_dp

    implicit none

    real(dp),               intent(in)    :: dt                   ! time step
    logical,  dimension(:), intent(in)    :: isc3                 ! C3 mask
    real(dp), dimension(:), intent(in)    :: Vcmax                ! Vcmax(Tl) [mol(air)/m2s]
    real(dp), dimension(:), intent(in)    :: GPP                  ! A+Rd [mol(co2)/m2s]
    real(dp), dimension(:), intent(in)    :: Rd                   ! Leaf respiration Rd [mol(co2)/m2s]
    real(dp), dimension(:), intent(in)    :: Gammastar            ! CO2 compensation point Gamma* [ppm]
    real(dp), dimension(:), intent(in)    :: ca                   ! Ambient CO2 concentration [ppm]
    real(dp), dimension(:), intent(in)    :: ci                   ! Stomatal CO2 concentration [ppm]
    real(dp), dimension(:), intent(in)    :: gsc                  ! Stomatal conductance for CO2 [mol(CO2)/m2s]
    real(dp), dimension(:), intent(in)    :: gac                  ! Aerodynamic conductance for CO2 [mol(CO2)/m2s]
    real(dp), dimension(:), intent(in)    :: gbc                  ! Boundary layer conductance for CO2 [mol(CO2)/m2s]
    real(dp), dimension(:), intent(in)    :: Tl                   ! Leaf temperature [K]
    real(dp), dimension(:), intent(in)    :: Rair                 ! Isotopic composition of ambient CO2
    real(dp), dimension(:), intent(inout) :: Vstarch              ! Transitory starch pool [mol(C)/m2]
    real(dp), dimension(:), intent(inout) :: Rstarch              ! Isotopic composition of transitory starch
    real(dp), dimension(:), intent(inout) :: Rsucrose             ! Isotopic composition of sucrose pool
    real(dp), dimension(:), intent(inout) :: Rphoto               ! Isotopic composition of pool for photorespiration
    real(dp), dimension(:), intent(out)   :: Disc                 ! Discrimination
    real(dp), dimension(:), intent(out)   :: Ass13                ! 13CO2 flux [mol(13CO2)/m2s]

    ! Local variables
    integer(i4) :: nn  ! number of grid points
    integer(i4) :: jl  ! counter
    real(dp)    :: tmp ! temporary real
    real(dp), dimension(size(isc3)) :: Rdm, Rds, leakage, Vp, Phi1 ! for C4
    real(dp), dimension(size(isc3)) :: Ass                         ! net assimilation
    real(dp), dimension(size(isc3)) :: k                           ! carboxylation efficiency
    real(dp), dimension(size(isc3)) :: Vc                          ! carboxylation rate
    real(dp), dimension(size(isc3)) :: Photo                       ! Photorespiration
    ! resistances and conductances
    real(dp), dimension(size(isc3)) :: rac, rbc, rsc, rmc, rwc, rmlc, rabsmc ! resistances for CO2
    real(dp), dimension(size(isc3)) :: gmc, gabsc, gabsmc                    ! conductances for CO2
    ! CO2 concentrations
    real(dp), dimension(size(isc3)) :: cc, cbs, ctmp ! , cs, cw
    ! fractionations
    real(dp), dimension(size(isc3)) :: eps_es, eps_b4, eps_b, eps_s, eps_a_eff
    real(dp) :: eps_ab
    real(dp) :: eps_night ! fractionation during night
    ! for pools update
    real(dp) :: add_flux  ! new flux to pool
    real(dp) :: add_r     ! isotope ratio of new flux

    nn = size(isc3,1)
    
    ! For C4
    Rdm = frdm*Rd  ! C4: respiration in mesophyll
    Rds = Rd - Rdm ! C4: respiration in bundle sheat
    ! L  - Leakage of bundle sheat cells
    ! Vp - Rate of PEP-Carboxylations
    ! Phi = L/Vp seems to be rather constant in C4 plants = 0.2,
    ! probably because the amount of C3 to C4 photosynthesis is regulated (von Caemmerer 2000, section 4.3.3.2).
    ! However, Vp=0 (and L<>0) at night and Phi->inf.
    ! Phi should be 1 if GPP equals the amount of respiration in the bundle sheats (GPP=Rds).
    ! We calc PhiL=1/Phi because this goes to zero at night. We take a linear relationship between
    ! PhiL and GPP, so that PhiL=1 at Rds and PhiL capped at 5=1/0.2.
    ! PhiL(GPP=0)=0; PhiL(GPP=Rds)=1 -> PhiL = 1/Rds*GPP < 5.
    ! Rds is changing with time but this effect of a changing slope is rather small.
    ! This comes to vp=gpp and leakage=rds if gpp<5*rds otherwise phiL=1/phi.
    tmp = 1.0_dp / (1.0_dp-Phi)
    where (GPP > (Rds/Phi))
       leakage = (GPP-Rds)*(Phi*tmp)
       Vp      = (GPP-Rds)*tmp
    elsewhere
       leakage = Rds
       Vp      = GPP
    end where
    ! New continuous Phi
    where (Vp > 0.0_dp)
       Phi1    = leakage/Vp
    elsewhere
       Phi1    = 0.0_dp ! dummy
    end where

    ! net assimilation
    Ass = GPP - Rd

    !
    !-- Resistances & Conductances
    !

    ! Aerodynamic resistance
    where (gac < huge(1.0_dp)*0.1_dp)
       rac = 1.0_dp / gac ! [mol(CO2) m-2 s-1]^-1
    elsewhere
       rac = tiny(1.0_dp)
    endwhere
    ! Leaf boundary layer resistance
    where (gbc < huge(1.0_dp)*0.1_dp)
       rbc = 1.0_dp / gbc ! [mol(CO2) m-2 s-1]^-1
    elsewhere
       rbc = tiny(1.0_dp)
    endwhere
    ! Stomatal resistance
    where (gsc < huge(1.0_dp)*0.1_dp)
       rsc = 1.0_dp / gsc ! [mol(CO2) m-2 s-1]^-1
    elsewhere
       rsc = tiny(1.0_dp)
    endwhere
    ! Mesophyll conductance [mol(CO2) m-2 s-1]
    ! ISOLSM takes 8000 for C4. John Evans suggests that it should be ca. 3000 for C3.
    ! gm for C4 is about twice that of C3 (Evans & v.Caemmerer 1996).
    ! John Evans is critical about temp dependence of gm from Bernacchi et al. (2002)
    ! and suggests: if any temp depence, take the one of Vcmax.
    where (isc3) ! c3
       gmc =        meso_cond_fac*Vcmax
    elsewhere    ! c4
       gmc = 2.0_dp*meso_cond_fac*Vcmax
    end where
    where (gmc > 0.0_dp)
       rmc = 1.0_dp / gmc ! [mol(co2) m-2 s-1]^-1
    elsewhere
       rmc = 0.0_dp       ! dummy
    end where
    ! Resistance from stoma middle to the chloroplast surface [mol(CO2) m-2 s-1]^-1
    ! 0.25 of mesophyll resistance: Lloyd & Farquhar (1994)
    rwc  = 0.25_dp * rmc
    ! Resistance from the chloroplast surface to inside the chloroplast [mol(CO2) m-2 s-1]^-1
    rmlc = (1.0_dp-0.25_dp) * rmc

    ! Total conductance from canopy to stomatal air space for CO2
    gabsc = 1.0_dp / (rac+rbc+rsc)
    ! Total conductance from canopy to sites of carboxylation for CO2
    rabsmc = rac + rbc + rsc + rwc + rmlc
    gabsmc = 1.0_dp / rabsmc
    
    !
    !-- CO2 concentrations
    !

    ! CO2 at the site of carboxylation, i.e. in the chloroplasts
    where ((gmc > 0.0_dp) .and. (ass > 0.0_dp))
       cc = ci - ass*rmc
       cc = max(cc, 1.1*Gammastar) ! from isolsm
    elsewhere
       cc = ci
    end where
    ! C4: CO2 in the bundle sheat (von Caemmerer 2000: Eq. 4.12 and Table 4.1)
    cbs = cc + leakage/gbsc
    ! ! CO2 at the leaf surface
    ! cs   = ca - ass*(rac+rbc) ! same: cs = ci + ass*rsc
    ! ! CO2 at the chloroplast surface = 1/4 down from Ci->Cc
    ! ! cf. rwc, rmlc above (0.25 from Lloyd & Farquhar 1994)
    ! cw = ci - 0.25_dp*(ci-cc)

    
    !
    !-- Other photosynthesis variables
    !

    ! Total conductance from canopy to sites of carboxylation
    ! Done again to account for possible cropping of Cc.
    ! Do not do rabsmc again so that calc of eps_a_eff is still correct.
    where (ca > cc) gabsmc = gabsc * (ca-ci)/(ca-cc)
    ! if (any(abs(gabsmc-1.0_dp/rabsmc)>(gabsmc*1e-6_dp))) then
    !    print*, 'GABSMC00 ', ass
    !    print*, 'GABSMC01 ', ca
    !    print*, 'GABSMC02 ', ci
    !    print*, 'GABSMC03 ', cc
    !    print*, 'GABSMC04 ', 1.1*Gammastar
    !    print*, 'GABSMC05 ', gabsmc
    !    print*, 'GABSMC06 ', 1.0_dp/rabsmc
    ! endif

    ! Carboxylation efficiency: initial slope of A vs Ci
    ! From Farquhar et al. (1982), eq. B11
    !    k_ci = GPP/(ci-Gammastar)
    ! Same for Cc in C3 and Cbs in C4
    where (isc3)
       ctmp = cc
    elsewhere
       ctmp = cbs
    end where
    where (ctmp > Gammastar)
       k = GPP/(ctmp-Gammastar)
    elsewhere
       k = 0.0_dp
    end where
    ! Carboxylation rate
    Vc = k*ctmp
    ! Photorespiration rate = 0.5*Vo
    photo = k*Gammastar
    
    !
    !-- Fractionations
    !
    
    ! diffusion fractionation through lamina
    eps_ab = 1.0_dp - (1.0_dp-eps_a)**twothird_dp
    ! frac. during CO2 dissolution
    eps_es = (1.18_dp - 0.0041_dp*(Tl-T0_dp))*1.e-3_dp ! Vogel et al. (Z Physik, 1970), Szaran (Chemical Geology, 1998)
    ! discrimination by PEP-c (<0)
    eps_b4 = (26.19_dp - 9483._dp/Tl)*1.e-3_dp         ! Henderson et al. (Aust J Plant Physiol, 1992)         
    ! effective discrimination of carboxylation in C3 plants
    eps_b = eps_b3*(1.0_dp-beta) + eps_b4*beta         ! Brugnoli & Farquhar (2000)
    ! frac. during leakage of bundle sheets in C4
    eps_s = eps_es + eps_al
    ! effective fractionation from canopy canopy to sites of carboxylation (chloroplast interior)
    eps_a_eff = (rac*0.0_dp + rbc*eps_ab + rsc*eps_a + rwc*eps_a + rmlc*eps_s) / (rabsmc)
    
    !
    !-- 13CO2 fluxes
    !
    
    eps_night = 0.0_dp
    do jl=1, nn
       if (k(jl) > 0.0_dp) then
          ! Day
          ! Same for photorespiration as Wingate et al. (2007) for leaf respiration
          if (isc3(jl)) then ! C3
             ! print*, 'DD01 ', jl, k(jl), gabsmc(jl)
             ! print*, 'DD02 ', eps_a_eff(jl), eps_b(jl), eps_f, eps_e
             ! print*, 'DD03 ', ca(jl), ci(jl), cc(jl), Gammastar(jl)
             ! print*, 'DD04 ', GPP(jl), Rd(jl)
             Ass13(jl) = (1.0_dp-eps_a_eff(jl)) * gabsmc(jl) / &
                  ( (1.0_dp-eps_a_eff(jl)) * gabsmc(jl) + (1.0_dp-eps_b(jl)) * k(jl) ) * &
                  ( (1.0_dp-eps_b(jl)) * Rair(jl) * k(jl) * ca(jl) - &
                  (1.0_dp-eps_f) * Rphoto(jl) * k(jl) * Gammastar(jl) - &
                  (1.0_dp-eps_e) * Rsucrose(jl) * Rd(jl) )
             ! print*, 'DD05 ', Ass(jl), Ass13(jl)
             if (ne(Ass(jl),0.0_dp)) then
                Disc(jl) = 1.0_dp - Ass13(jl) / (Ass(jl)*Rair(jl))
                ! print*, 'DD06 ', Ass13(jl)/Ass(jl), Disc(jl)*1000._dp
             else
                Disc(jl) = 0.0_dp
             end if             
             ! print*, 'DD07 ', ass(jl), gmc(jl), Vcmax(jl)
          else               ! C4
             tmp = cc(jl) * ( 1.0_dp - (1.0_dp-eps_a_eff(jl)) / (1.0_dp-eps_b4(jl)) * Ass(jl)/Vp(jl) - &
                  (1.0_dp-eps_s(jl)) / (1.0_dp-eps_b3) * (1.0_dp-eps_a_eff(jl)) / (1.0_dp-eps_b4(jl)) * &
                  Phi1(jl) * Ass(jl) / vc(jl) )
             if (abs(1.0_dp-tmp/ca(jl)) > 0.01_dp) then
                Ass13(jl) = (1.0_dp-eps_a_eff(jl)) * ( Rair(jl)*ca(jl) - & ! no *ass
                     (1.0_dp-eps_e) / (1.0_dp-eps_b4(jl)) * Rdm(jl) / Vp(jl) * Rsucrose(jl) * cc(jl) - &
                     (1.0_dp-eps_s(jl)) / (1.0_dp-eps_b3) / (1.0_dp-eps_b4(jl)) * Phi1(jl) * cc(jl) * &
                     ( (1.0_dp-eps_e) * Rsucrose(jl) * Rd(jl) + &
                     (1.0_dp-eps_f) * Photo(jl) * Rphoto(jl) ) / vc(jl) ) / &
                     (ca(jl)-tmp)
             else
                Ass13(jl) = 0.0_dp
             end if
             Disc(jl)  = 1.0_dp - Ass13(jl)/Rair(jl)
             Ass13(jl) = Ass13(jl) * Ass(jl)                               ! *ass
          endif              ! end C3/C4
       else
          ! Night
          Ass13(jl) = (1.0_dp-eps_night) * Rstarch(jl) ! no *ass
          Disc(jl)  = 1.0_dp - Ass13(jl)/Rair(jl)
          Ass13(jl) = Ass13(jl) * Ass(jl)              ! *ass
       end if ! end day/night

       !
       !-- 13CO2 pools
       !
       
       ! Update substrate pools
       ! 2/3 of assimilation is sugar, 1/3 goes into starch pool
       if (Ass(jl) > 0.0_dp) then
          ! average sucrose pool (continuous)
          add_flux     = dt * twothird_dp * Ass(jl)
          add_r        = (1.0_dp-eps_sucrose) * Ass13(jl)/Ass(jl)
          Rsucrose(jl) = (Vsucrose*Rsucrose(jl) + add_flux*add_r) / (Vsucrose + add_flux)
          ! average photorespiration pool (continuous)
          add_flux    = dt * 2.0_dp * Photo(jl) ! vo=2*photo
          add_r       = (1.0_dp-eps_sucrose) * Ass13(jl)/Ass(jl)
          Rphoto(jl)  = (Vphoto*Rphoto(jl) + add_flux*add_r) / (Vphoto + add_flux)
          ! integrated starch pool (new at every day)
          add_flux    = dt * onethird_dp * Ass(jl)
          add_r       = (1.0_dp-eps_starch) * Ass13(jl)/Ass(jl)
          Rstarch(jl) = (Vstarch(jl)*Rstarch(jl) + add_flux*add_r) / (Vstarch(jl) + add_flux)
          Vstarch(jl) = Vstarch(jl) + add_flux
       else
          ! Rsucrose(jl) = Rsucrose(jl)
          ! Rphoto(jl)   = Rphoto(jl)
          ! Rstarch(jl)  = Rstarch(jl)
          Vstarch(jl) = 0.0_dp ! start new starch integration each day
       end if

    end do

    return

  end subroutine c13o2_discrimination


  ! ------------------------------------------------------------------

  !     NAME
  !         c13o2_discrimination_simple

  !     PURPOSE
  !>        \brief Simple model of 13C leaf discrimination.

  !>        \details 13C leaf discrimination D after Farquhar et al. (Aust J Plant Physiol, 1982)
  !>        and Farquhar (Aust J Plant Physiol, 1983), assuming no day/dark nor photo-respiration
  !>        and using the CO2 mixing ratio in the sub-stomatal cavity.
  !>        This is for C3 during daytime: \n
  !>            D = a + (b-a) * ci/ca \n
  !>        where a is 13C fractionation during diffusion (-4.4 permil),
  !>        b is the effective fractionation during RUBISCO carboxylation (inifinite mesophyll conductance) (27. permil),
  !>        ca is the ambient mixing ratio outside the leaf boundary layer, and
  !>        ci is the mixing ratio inside the sub-stomatal cavity. \n
  !>        The routine returns discrimination D and also the 13CO2 net assimilation rate A': \n
  !>            A' = (1-D) * A * Ra \n
  !>        with the 12CO2 net assimilation rate A and the isotope ratio of the ambient air Ra.
  !>        A starch pool is built up during day by 1/3 of assimilation without any fractionation.
  !>        Mitochondrial respiration respires this starch pool at night, with the average
  !>        isotopic composition of the pool at the beginning of the night, without fractionation. \n

  !     CALLING SEQUENCE
  !         call c13o2_discrimination_simple(dt, isc3, GPP, Rd, ca, ci, Tl, Rair, Vstarch, Rstarch, Disc, Ass13)

  !>        This is an elemental subroutine and can be called with scalars and arrays.

  !     PARAMETER
  !>        \param[in]    "real(dp) :: dt"        Time step [s]
  !>        \param[in]    "logical  :: isc3"      C3 mask with .true. for C3 and .false. for C4 plants
  !>        \param[in]    "real(dp) :: GPP"       Gross assimilation [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: Rd"        Leaf respiration [mol(CO2)/m^2/s]
  !>        \param[in]    "real(dp) :: ca"        Ambient CO2 mixing ratio [mol(CO2)/mol(air)]
  !>        \param[in]    "real(dp) :: ci"        Stomatal CO2 mixing ratio [mol(CO2)/mol(air)]
  !>        \param[in]    "real(dp) :: Tl"        Leaf temperature [K]
  !>        \param[in]    "real(dp) :: Rair"      Isotope ratio of ambient CO2 ([13C]/[12C])
  !>        \param[inout] "real(dp) :: Vstarch"   Transitory starch pool [mol(C)/m2]
  !>        \param[inout] "real(dp) :: Rstarch"   Isotope ratio of transitory starch ([13C]/[12C])
  !>        \param[out]   "real(dp) :: Disc"      Leaf discrimination
  !>        \param[out]   "real(dp) :: Ass13"     13C net assimilation [mol(13CO2)/m^2/s]

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure subroutine c13o2_discrimination_simple( &
       ! -- Input
       dt, isc3, &
       ! GPP and Leaf respiration
       GPP, Rd, &
       ! Ambient and stomatal CO2 mixing ratio
       ca, ci, &
       ! leaf temperature
       Tl, &
       ! Ambient isotope ratio
       Rair, &
       ! -- Inout
       ! Starch pool and its isotope ratio
       Vstarch, Rstarch, &
       ! -- Output
       ! discrimination
       Disc, &
       ! 13CO2 flux
       Ass13)

    use mo_constants, only: onethird_dp, T0_dp

    implicit none

    real(dp), intent(in)    :: dt      ! time step
    logical,  intent(in)    :: isc3    ! C3 mask
    real(dp), intent(in)    :: GPP     ! A+Rd [mol(CO2)/m2s]
    real(dp), intent(in)    :: Rd      ! Leaf respiration Rd [mol(CO2)/m2s]
    real(dp), intent(in)    :: ci      ! Stomatal CO2 concentration [mol(CO2)/mol(air)]
    real(dp), intent(in)    :: ca      ! Ambient CO2 concentration [mol(CO2)/mol(air)]
    real(dp), intent(in)    :: Tl      ! Leaf temperature [K]
    real(dp), intent(in)    :: Rair    ! Isotopic composition of ambient CO2
    real(dp), intent(inout) :: Vstarch ! Transitory starch pool [mol(C)/m2]
    real(dp), intent(inout) :: Rstarch ! Isotopic composition of transitory starch
    real(dp), intent(out)   :: Disc    ! Discrimination
    real(dp), intent(out)   :: Ass13   ! 13CO2 flux [mol(13CO2)/m2s]

    ! Local variables
    real(dp) :: tmp
    real(dp) :: Ass                   ! net assimilation
    ! fractionations
    real(dp) :: eps_es, eps_b4, eps_s ! for C4
    ! for pool update
    real(dp) :: add_flux  ! new flux to pool
    real(dp) :: add_r     ! isotope ratio of new flux

    ! net assimilation
    Ass = GPP - Rd

    ! For C4
    ! frac. during CO2 dissolution
    eps_es = (1.18_dp - 0.0041_dp*(Tl-T0_dp))*1.e-3_dp ! Vogel et al. (Z Physik, 1970), Szaran (Chemical Geology, 1998)
    ! discrimination by PEP-c (<0)
    eps_b4 = (26.19_dp - 9483._dp/Tl)*1.e-3_dp         ! Henderson et al. (Aust J Plant Physiol, 1992)
    ! frac. during leakage of bundle sheets in C4
    eps_s = eps_es + eps_al

    !
    !-- 13CO2 fluxes
    !

    if (ca > ci) then ! day
       if (isc3) then ! C3
          Disc = eps_a + (eps_b_ci-eps_a) * ci/ca
       else           ! C4
          tmp = ci * ( 1.0_dp - (1.0_dp-eps_a) / (1.0_dp-eps_b4) * (1.0_dp-Phi) - &
               (1.0_dp-eps_s) / (1.0_dp-eps_b3) * (1.0_dp-eps_a) / (1.0_dp-eps_b4) * Phi )
          Ass13 = (1.0_dp-eps_a) * ca / (ca-tmp) ! no *A*Ra
          Disc  = 1.0_dp - Ass13                 ! D = 1 - A'/A/Ra
       end if ! end C3/C4
       Ass13 = (1.0_dp-Disc) * Ass * Rair        ! D = 1 - A'/A/Ra
    else              ! night
       Ass13 = Rstarch     ! no *ass
       Disc  = 1.0_dp - Ass13/Rair
       Ass13 = Ass13 * Ass ! *ass
    endif ! end day/night

    !
    !-- Starch pool
    !

    ! Update substrate pools
    ! 2/3 of assimilation is sugar, 1/3 goes into starch pool
    if (ca > ci) then ! day
       ! integrated starch pool (new at every day)
       add_flux = dt * onethird_dp * Ass
       add_r    = Ass13/Ass
       Rstarch  = (Vstarch*Rstarch + add_flux*add_r) / (Vstarch + add_flux)
       Vstarch  = Vstarch + add_flux
    else
       ! Rstarch  = Rstarch
       Vstarch = 0.0_dp ! start new starch integration each day
    end if

    return

  end subroutine c13o2_discrimination_simple


  ! ------------------------------------------------------------------

  !     NAME
  !         init_starch_pool

  !     PURPOSE
  !>        \brief Initialise transitory starch pool and its isotope ratio.

  !>        \details Initialises the transitory starch pool and its isotope ratio
  !>        in the chloroplasts.

  !     CALLING SEQUENCE
  !         call init_starch_pool(isc3, Vstarch, Rstarch, Rinitc3, Rinitc4)

  !>        This is an elemental subroutine and can be called with scalars and arrays.

  !     PARAMETER
  !>        \param[in]    "logical  :: isc3"       C3 mask with .true. for C3 and .false. for C4 plants
  !>        \param[inout] "real(dp) :: Vstarch"    Transitory starch pool
  !>        \param[inout] "real(dp) :: Rstarch"    Isotope ratio of transitory starch
  !>        \param[in]    "real(dp) :: Rinitc3"    Initial isotope ratio if isc3=.true., i.e. of C3 plants
  !>        \param[in]    "real(dp) :: Rinitc4"    Initial isotope ratio if isc3=.false., i.e. of C4 plants

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr-Nov 2019
  elemental pure subroutine init_starch_pool(isc3, Vstarch, Rstarch, Rinitc3, Rinitc4)

    implicit none

    logical,  intent(in)    :: isc3
    real(dp), intent(inout) :: Vstarch
    real(dp), intent(inout) :: Rstarch
    real(dp), intent(in)    :: Rinitc3
    real(dp), intent(in)    :: Rinitc4
    
    ! initialise
    Vstarch  = 0.0_dp
    if (isc3) then
       Rstarch  = Rinitc3
    else
       Rstarch  = Rinitc4
    endif

    return

  end subroutine init_starch_pool


  ! ------------------------------------------------------------------

  !     NAME
  !         init_sugar_pools

  !     PURPOSE
  !>        \brief Initialise isotope ratios of leaf carbohydrate pools.

  !>        \details Initialises the isotope ratios of the leaf carbohydrate pools
  !>        for sucrose and the pool for photorespiration.

  !     CALLING SEQUENCE
  !         call init_sugar_pools(isc3, Rsucrose, Rphoto, Rinitc3, Rinitc4)

  !>        This is an elemental subroutine and can be called with scalars and arrays.

  !     PARAMETER
  !>        \param[in]    "logical  :: isc3"       C3 mask with .true. for C3 and .false. for C4 plants
  !>        \param[inout] "real(dp) :: Rsucrose"   Isotope ratio of leaf sucrose
  !>        \param[inout] "real(dp) :: Rphoto"     Isotope ratio of substrate for photorespiration
  !>        \param[in]    "real(dp) :: Rinitc3"    Initial isotope ratio if isc3=.true., i.e. of C3 plants
  !>        \param[in]    "real(dp) :: Rinitc4"    Initial isotope ratio if isc3=.false., i.e. of C4 plants

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr-Nov 2019
  elemental pure subroutine init_sugar_pools(isc3, Rsucrose, Rphoto, Rinitc3, Rinitc4)

    implicit none

    logical,  intent(in)    :: isc3
    real(dp), intent(inout) :: Rsucrose
    real(dp), intent(inout) :: Rphoto
    real(dp), intent(in)    :: Rinitc3
    real(dp), intent(in)    :: Rinitc4
    
    ! initialise
    if (isc3) then
       Rsucrose = Rinitc3
       Rphoto   = Rinitc3
    else
       Rsucrose = Rinitc4
       Rphoto   = Rinitc4
    endif

    return

  end subroutine init_sugar_pools

END MODULE mo_c13o2_photosynthesis
