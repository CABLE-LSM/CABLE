! UCLEM urban canopy model

! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)

! This file is part of the UCLEM urban canopy model
!
! UCLEM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UCLEM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UCLEM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! This code was originally inspired by the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)
! The snow scheme is based on Douville, Royer and Mahfouf, Climate Dynamics, 12, p21 (1995)
! The in-canyon vegetation is based on Kowalczyk et al, DAR Tech Paper 32 (1994), but simplified by assuming sigmaf=1.

! The important differences from previous models include an alternative formulation of the conduction model by
! Lipson et al (2017), an alternative formulation for in-canyon aerodynamical resistances based on Harman, et al (2004)
! and Kanada et al (2007), combined with a second canyon wall for completeness.  The scheme includes nrefl order
! reflections in the canyon for both longwave and shortwave radiation (in TEB infinite reflections are used for
! shortwave and 1st order reflections in longwave). A big-leaf vegetation tile is included in the canyon using the
! Kowalczyk et al (1994) scheme but with a simplified soil moisture budget and no modelling of the soil temperature
! since sigmaf=1.  Snow is also included in the canyon and on roofs using a single-layer scheme based on
! Douville, et al (1995).  Time dependent traffic heat fluxes are based on Coutts, et al (2007).

! The Urban CLimate Energy Model (UCLEM) results when the exterior canyon model is combined with an interior model
! of energy use from Lipson, et al (2018).  This model includes a representation of energy storage and conduction for the
! interior slab, prognostic internal air temperature, as well as dynamic anthropogenic energy generation as a function
! of internal air temperature and a representation of human behaviours over the diurnal cycle based on Thatcher (2007).


module uclem_science

implicit none

private



public uclem_eval, uclem_most_profile, uclem_canyon_roughness, init_lwcoeff, init_internal, getswcoeff, getnetalbedo
public uclem_calc_alb, uclem_prepinterior, uclem_rebuild_derived
public facetdata, facetparams, hydrodata, vegdata, intldata, fparmdata, pdiagdata
public urbtemp,energytol,resmeth,zohmeth,acmeth,nrefl,                                 &
       scrnmeth,wbrelaxc,wbrelaxr,ncyits,nfgits,tol,                                   &
       zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden,refheight,     &
       zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf,        &
       intairtmeth,intmassmeth,statsmeth,lwintmeth,cvcoeffmeth,infilmeth,acfactor,     &
       ac_heatcap,ac_coolcap,ac_deltat,cyc_prop,cyc_base,cyc_traf,cyc_tran,icmax
public nl, vkar, waterden, pi
! lf is exported so the host can convert pd%u_melt, which is an energy flux, back
! to the mass flux the host reports as melt.  Using UCLEM's own value keeps that
! conversion consistent with the energy balance u_melt was taken from.
public lf


type facetdata
  real(kind=8), dimension(:,:), allocatable :: nodetemp        ! Temperature of node (prognostic)       [K]
  real(kind=8), dimension(:,:), allocatable :: storage ! Facet energy storage (diagnostic)
end type facetdata

type facetparams
  real, dimension(:,:), allocatable :: depth         ! Layer depth                              [m]
  real, dimension(:,:), allocatable :: volcp         ! Layer volumetric heat capacity           [J m^-3 K-1]
  real, dimension(:,:), allocatable :: lambda        ! Layer conductivity                       [W m^-1 K^-1]
  real, dimension(:),   allocatable :: alpha         ! Facet albedo (internal & external)
  real, dimension(:),   allocatable :: emiss         ! Facet emissivity (internal & external)
end type facetparams


type hydrodata
  real, dimension(:), allocatable   :: surfwater
  real, dimension(:), allocatable   :: leafwater
  real, dimension(:), allocatable   :: soilwater
  real, dimension(:), allocatable   :: snow
  real, dimension(:), allocatable   :: den
  real, dimension(:), allocatable   :: snowalpha
end type hydrodata


type vegdata
  real, dimension(:), allocatable :: temp          ! Temperature of veg (prognostic)  [K]
  real, dimension(:), allocatable :: sigma         ! Fraction of veg on roof/canyon
  real, dimension(:), allocatable :: alpha         ! Albedo of veg
  real, dimension(:), allocatable :: emiss         ! Emissivity of veg
  real, dimension(:), allocatable :: lai           ! Leaf area index of veg
  real, dimension(:), allocatable :: zo            ! Roughness of veg
  real, dimension(:), allocatable :: rsmin         ! Minimum stomatal resistance of veg
end type vegdata


type intldata
  real, dimension(:,:,:), allocatable :: psi   ! internal radiation
  real, dimension(:,:,:), allocatable :: viewf ! internal radiation
end type intldata


type fparmdata
  real, dimension(:), allocatable :: hwratio,effhwratio,sigmabld
  real, dimension(:), allocatable :: industryfg,intgains_flr,trafficfg,bldheight,bldwidth
  real, dimension(:), allocatable :: ctime,vangle,hangle,fbeam,weekdayload
  real, dimension(:), allocatable :: bldairtemp
  real, dimension(:), allocatable :: swilt,sfc,ssat,rfvegdepth
  real, dimension(:), allocatable :: infilach,ventilach,heatprop,coolprop
  real, dimension(:), allocatable :: sigmau
  integer, dimension(:), allocatable :: intmassn
end type fparmdata


type pdiagdata
  real, dimension(:), allocatable :: lzom, lzoh, cndzmin, cduv, cdtq, zol
  real, dimension(:), allocatable :: tscrn, qscrn, uscrn, u10, emiss, snowmelt
  real, dimension(:), allocatable :: bldheat, bldcool, traf, intgains_full
  real, dimension(:), allocatable :: irrig,acond_vegw
  real, dimension(:), allocatable :: surfrunoff,soilwetness,soilwater
  real, dimension(:), allocatable :: transveg,soilmoist,delsoilmoist
  real, dimension(:), allocatable :: rootmoistc
  real, dimension(:), allocatable :: frac_sigma
  real, dimension(:), allocatable :: delintercept, snowt, vegt, swe, surfstor
  real, dimension(:), allocatable :: snowfrac, salbedo, calbedo, taircanyon, ulai
  real, dimension(:), allocatable :: delswe
  real, dimension(:), allocatable :: roof_water_runoff, roof_snow_runoff, roof_soil_runoff
  real, dimension(:), allocatable :: road_water_runoff, road_snow_runoff, road_soil_runoff
  real, dimension(:), allocatable :: ac_heat_on, ac_cool_on
  real(kind=8), dimension(:), allocatable :: surferr, atmoserr, surferr_bias, atmoserr_bias
  real(kind=8), dimension(:), allocatable :: storage_flux
  ! real(kind=8), dimension(:,:), allocatable :: storagetot_road, storagetot_walle, storagetot_wallw, storagetot_roof
  real, dimension(:), allocatable :: u_alb, u_melt
  real(kind=8), dimension(:), allocatable :: d_atmosflux
  real, dimension(:), allocatable :: acond_roof, acond_road, acond_walle, acond_wallw
  real, dimension(:), allocatable :: d_tempr, d_canyontemp
  real, dimension(:), allocatable :: rg_roof, rg_road, rg_walle, rg_wallw
  real, dimension(:), allocatable :: rg_vegc, rg_vegr, rg_rfsn, rg_rdsn, d_canyonrgout, d_roofrgout
  real, dimension(:), allocatable :: sg_roof, sg_road, sg_walle, sg_wallw
  real, dimension(:), allocatable :: sg_vegc, sg_vegr, sg_rfsn, sg_rdsn
  real, dimension(:), allocatable :: acond_vegc, acond_vegr, acond_rfsn, acond_rdsn
  real, dimension(:), allocatable :: d_rfsndelta, d_rdsndelta, rfveg_sigma, cnveg_sigma
  real, dimension(:), allocatable :: rdsntemp, rfsntemp, rfveg_temp
  real, dimension(:), allocatable :: d_ac_inside, int_infilflux, ac_coeff, d_totdepth
  real, dimension(:), allocatable :: eg_roof, eg_vegr, eg_road, eg_vegc, eg_rdsn, eg_rfsn

  logical :: first_call
end type pdiagdata

! default namelist here
! model parameters
integer, save      :: resmeth=3            ! Canyon sensible heat transfer (0=Masson, 2=Kusaka, 3=Harman (fixed width))
integer, save      :: zohmeth=1            ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, save      :: acmeth=1             ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, 3=COP of 1.0)
integer, save      :: intairtmeth=0        ! Internal air temperature (0=prescribed, 1=aggregated varying, 2=fractional varying)
integer, save      :: intmassmeth=2        ! Internal thermal mass (0=none, 1=one floor, 2=dynamic floor number)
integer, save      :: cvcoeffmeth=1        ! Internal surface convection heat transfer coefficient (0=DOE, 1=ISO6946, 2=fixed)
integer, save      :: statsmeth=1          ! Use statistically based diurnal QF amendments (0=off, 1=on) from Thatcher 2007
integer, save      :: infilmeth=0          ! Method to calculate infiltration rate (0=constant, 1=EnergyPlus/BLAST, 2=ISO)
integer, save      :: nrefl=3              ! Number of canyon reflections for radiation (default=3)
integer, save      :: scrnmeth=1           ! Screen diagnostic method (0=Slab, 1=Hybrid, 2=Canyon)
integer, save      :: wbrelaxc=0           ! Relax canyon soil moisture for irrigation (0=Off, 1=On)
integer, save      :: wbrelaxr=0           ! Relax roof soil moisture for irrigation (0=Off, 1=On)
integer, save      :: lwintmeth=2          ! Interior longwave flux method (0=off, 1=explicit, 2=with forward step facet temps)
integer, parameter :: nl=4                 ! Number of layers (default 4, must be a multiple of 4)
integer, save      :: iqt=314              ! Diagnostic point (in terms of host grid)
! sectant solver parameters
integer, save      :: ncyits=6             ! Number of iterations for balancing canyon sensible and latent heat fluxes (default=6)
integer, save      :: nfgits=3             ! Number of iterations for balancing veg and snow energy budgets (default=3)
real, save         :: tol=0.001            ! Sectant method tolerance for sensible heat flux (default=0.001)
real, save         :: infilalpha=0.5       ! Weighting for dampening factor when calculating internal infiltration
real(kind=8), save :: energytol=0.005_8    ! Tolerance for acceptable energy closure in each timestep
real, save         :: urbtemp=290.         ! reference temperature to improve precision
! physical parameters
real, parameter    :: waterden=1000.       ! water density (kg m^-3)
real, parameter    :: icelambda=2.22       ! conductance of ice (W m^-1 K^-1)
real, parameter    :: aircp=1004.64        ! Heat capacity of dry air (J kg^-1 K^-1)
real, parameter    :: icecp=2100.          ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter    :: grav=9.80616         ! gravity (m s^-2)
real, parameter    :: vkar=0.4             ! von Karman constant
real, parameter    :: lv=2.501e6           ! Latent heat of vaporisation (J kg^-1)
real, parameter    :: lf=3.337e5           ! Latent heat of fusion (J kg^-1)
real, parameter    :: ls=lv+lf             ! Latent heat of sublimation (J kg^-1)
real, parameter    :: pi=3.14159265        ! pi (must be rounded down for shortwave)
real, parameter    :: rd=287.04            ! Gas constant for dry air
real, parameter    :: rv=461.5             ! Gas constant for water vapour
real, parameter    :: sbconst=5.67e-8      ! Stefan-Boltzmann constant
! snow parameters
real, save         :: zosnow=0.001         ! Roughness length for snow (m)
real, save         :: snowemiss=1.         ! snow emissivity
real, save         :: maxsnowalpha=0.85    ! max snow albedo
real, save         :: minsnowalpha=0.5     ! min snow albedo
real, save         :: maxsnowden=300.      ! max snow density (kg m^-3)
real, save         :: minsnowden=100.      ! min snow density (kg m^-3)
! generic urban parameters
real, save         :: refheight=0.6        ! Displacement height as a fraction of building height (Kanda et al 2007)
real, save         :: zomratio=0.10        ! Ratio of roughness length to building height (default=0.1 or 10%)
real, save         :: zocanyon=0.01        ! Roughness length of in-canyon surfaces (m)
real, save         :: zoroof=0.01          ! Roughness length of roof surfaces (m)
real, save         :: maxrfwater=1.        ! Maximum roof water (kg m^-2)
real, save         :: maxrdwater=1.        ! Maximum road water (kg m^-2)
real, save         :: maxrfsn=1.           ! Maximum roof snow (kg m^-2)
real, save         :: maxrdsn=1.           ! Maximum road snow (kg m^-2)
real, save         :: maxvwatf=0.1         ! Factor multiplied to LAI to predict maximum leaf water (kg m^-2)
real, save         :: acfactor=2.          ! Air conditioning inefficiency factor
real, save         :: ac_heatcap=10.       ! Maximum heating/cooling capacity (W m^-3)
real, save         :: ac_coolcap=10.       ! Maximum heating/cooling capacity (W m^-3)
real, save         :: ac_deltat=1.         ! Comfort range for temperatures (+-K)
! atmosphere stability parameters
integer, save      :: icmax=5              ! number of iterations for stability functions (default=5)
real, save         :: a_1=1.
real, save         :: b_1=2./3.
real, save         :: c_1=5.
real, save         :: d_1=0.35
! diurnal cycle profiles (start & end at 1am)
! traffic diurnal cycle weights approximated from Chapman et al., 2016
real, dimension(25), save :: cyc_traf = (/ 0.17, 0.12, 0.12, 0.17, 0.37, 0.88, &
                                           1.29, 1.48, 1.37, 1.42, 1.5 , 1.52, &
                                           1.50, 1.57, 1.73, 1.84, 1.84, 1.45, &
                                           1.01, 0.77, 0.65, 0.53, 0.41, 0.27, 0.17 /)
! base electricity demand cycle weights from Thatcher (2007), mean for NSW, VIC, QLD, SA
real, dimension(25), save :: cyc_base = (/ 0.92, 0.86, 0.81, 0.78, 0.8 , 0.87, &
                                           0.98, 1.06, 1.08, 1.09, 1.09, 1.09, &
                                           1.08, 1.08, 1.06, 1.06, 1.08, 1.11, &
                                           1.08, 1.06, 1.03, 1.00, 0.98, 0.95, 0.92 /)
! normalised proportion of heating/cooling appliances in use approximated from Thatcher (2007)
real, dimension(25), save :: cyc_prop = (/ 0.43, 0.40, 0.36, 0.33, 0.31, 0.36, &
                                           0.53, 0.66, 0.70, 0.64, 0.60, 0.60, &
                                           0.64, 0.69 ,0.75, 0.81, 0.88, 1.00, &
                                           0.99, 0.91, 0.81, 0.69, 0.57, 0.43, 0.43 /)
! internal temperature translation cycle from Thatcher (2007) (only used for intairtmeth=0)
real, dimension(25), save :: cyc_tran = (/ -1.09, -1.21, -2.12, -2.77, -3.06, -2.34, &
                                           -0.37,  1.03,  1.88,  2.37,  2.44,  2.26, &
                                            1.93,  1.41,  0.74,  0.16,  0.34,  1.48, &
                                            1.03,  0.14, -0.74, -1.17, -1.15, -1.34, -1.09/)




contains



! progcalc = 0   Update prognostic and diagnostic variables
! progcalc = 1   Update diagnostic variables only

subroutine uclem_eval(u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl,ddt,a_sg,a_rg,a_rho,a_temp,a_temproof, &
                      a_mixr,a_mixrroof,a_ps,a_umag,a_umagroof,a_udir,a_rnd,a_snd,a_zmin,a_zroof,fp,  &
                      fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,rdhyd,rfhyd,rfveg,road,roof,    &
                      room,slab,walle,wallw,cnveg,intl,mp,progcalc,diag,qf,u_temp9,u_shf9)

implicit none

integer, intent(in) :: mp, progcalc, diag
real, intent(in) :: ddt
integer, parameter :: nl=4                 ! Number of layers (default 4, must be a multiple of 4)
integer k


real, dimension(mp), intent(in) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(mp), intent(in) :: a_temproof,a_mixrroof,a_umagroof,a_zroof

real, dimension(mp) :: c_a_sg,c_a_rg,c_a_rho,c_a_temp,c_a_mixr,c_a_ps,c_a_umag,c_a_udir,c_a_rnd,c_a_snd,c_a_zmin
real, dimension(mp) :: c_a_temproof,c_a_mixrroof,c_a_umagroof,c_a_zroof


real, dimension(mp), intent(out) :: u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl
real, dimension(mp), intent(out) :: qf
real, dimension(mp,9), intent(out) :: u_temp9
real, dimension(mp,9), intent(out) :: u_shf9

real, dimension(mp) :: ggint_roof,ggint_walle,ggint_wallw,ggint_road,ggint_slab,ggint_intm
real, dimension(mp) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt,garfsn,gardsn
real, dimension(mp) :: wallpsi,roadpsi,fgtop,egtop,qsatr,qsata
real, dimension(mp) :: cu,fgrooftop,egrooftop
real, dimension(mp) :: we,ww,wr,zolog,a,n,zom,zonet,dis
real, dimension(mp) :: roofvegwetfac,roadvegwetfac,d_irrigwater
real, dimension(mp) :: z_on_l, pa, dts, dtt
real, dimension(mp) :: u_alb, u_melt, evspsbltop, evspsblrooftop, sbltop, sblrooftop
real, dimension(mp) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(mp) :: rg_roof,rg_road,rg_walle,rg_wallw,rg_vegc,rg_vegr,rg_rfsn,rg_rdsn
real, dimension(mp) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab,zero_flux
real, dimension(mp) :: fg_roof,fg_road,fg_walle,fg_wallw,fg_vegc,fg_vegr,fg_rfsn,fg_rdsn
real, dimension(mp) :: eg_roof,eg_road,eg_vegc,eg_vegr,eg_rfsn,eg_rdsn
real, dimension(mp) :: acond_roof,acond_road,acond_walle,acond_wallw
real, dimension(mp) :: acond_vegc,acond_vegr,acond_rfsn,acond_rdsn
real, dimension(mp) :: abase_road,abase_walle,abase_wallw,abase_vegc,abase_rdsn
real, dimension(mp) :: d_roofdelta,d_roaddelta,d_vegdeltac,d_vegdeltar,d_rfsndelta,d_rdsndelta
real, dimension(mp) :: d_tempc,d_mixrc,d_tempr,d_mixrr,d_sigd,d_sigr,d_rfdzmin
real, dimension(mp) :: d_ac_canyon,d_canyonrgout,d_roofrgout,d_tranc,d_evapc,d_tranr,d_evapr
real, dimension(mp) :: d_totdepth,d_netemiss,d_netrad,d_topu
real, dimension(mp) :: d_cwa,d_cw0,d_cww,d_cwr,d_cra,d_crr,d_crw
real, dimension(mp) :: d_canyontemp,d_canyonmix,d_traf
real, dimension(mp) :: ggext_roof, ggext_walle, ggext_wallw, ggext_road, ggext_slab, ggext_intm
real, dimension(mp) :: ggext_impl, ggint_impl, tmp
real, dimension(mp) :: ac_coeff, ggint
real, dimension(mp) :: d_ac_inside, d_intgains_bld, int_infilflux
real, dimension(mp) :: cyc_traffic, cyc_basedemand, cyc_proportion, cyc_translation
real, dimension(mp) :: int_infilfg, delintercept, rf_delsnow, rd_delsnow
real, dimension(mp) :: iroomtemp
real, dimension(mp) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(mp) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(mp) :: d_roomstor
real, dimension(mp) :: d_surfwater, d_soilwater, d_leafwater, d_snow, d_den, d_snowalpha
real, dimension(mp) :: fp_coeffbldheight

real, dimension(mp,1:4)  :: iskintemp          ! interim internal skin temperature
real(kind=8), dimension(mp,0:nl) :: inodetemp  ! interim tridiagonal solution



type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(hydrodata), intent(inout) :: rdhyd, rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(vegdata), intent(inout) :: cnveg
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd
!

u_fg = 0
u_eg = 0
u_ts = 0
u_wf = 0
u_rn = 0
u_evspsbl = 0
u_sbl = 0
qf = 0

c_a_temp = a_temp - urbtemp ! 273.15
c_a_temproof = a_temproof - urbtemp ! 273.15
!------------------------------------------------------------------------------
! Calculate fluxes
!------------------------------------------------------------------------------

! calculate water and snow area cover fractions
d_roofdelta = max(rfhyd%surfwater/maxrfwater,0.)**(2./3.)
d_roaddelta = max(rdhyd%surfwater/maxrdwater,0.)**(2./3.)
d_vegdeltac = max(rdhyd%leafwater/max(maxvwatf*cnveg%lai,1.E-8),0.)**(2./3.)
d_vegdeltar = max(rfhyd%leafwater/max(maxvwatf*rfveg%lai,1.E-8),0.)**(2./3.)

d_rfsndelta = rfhyd%snow/(rfhyd%snow+maxrfsn)
d_rdsndelta = rdhyd%snow/(rdhyd%snow+maxrdsn)


! Mar23 - cal d_mixrc and d_mixrr

! canyon level air temp and water vapor (displacement height at refheight*building height)
pa      = a_ps*exp(-grav*a_zmin/(rd*(c_a_temp+urbtemp)))
d_sigd  = a_ps
a       = (d_sigd/pa)**(rd/aircp)

!d_tempc = (a_temp+urbtemp)*a - urbtemp
d_tempc = c_a_temp*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempc,d_sigd)
call getqsat(qsata,c_a_temp,pa)

d_mixrc = a_mixr*qsatr/qsata


! roof level air temperature and water vapor (displacement height at building height)
pa      = a_ps*exp(-grav*a_zroof/(rd*(c_a_temproof+urbtemp)))
d_sigr  = a_ps*exp(-grav*fp%bldheight*(1.-refheight)/(rd*(c_a_temproof+urbtemp)))
a       = (d_sigr/pa)**(rd/aircp)
!d_tempr = (a_temproof+urbtemp)*a - urbtemp
d_tempr = c_a_temproof*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempr,d_sigr)
call getqsat(qsata,c_a_temproof,pa)

d_mixrr = a_mixrroof*qsatr/qsata


! calculate soil data
d_totdepth = sum(fp_road%depth,2)
! variables for canyon geometry with vegetation
fp_coeffbldheight = fp%effhwratio/fp%hwratio





! Mar 24
! calculate shortwave reflections
! Here we modify the effective canyon geometry to account for in-canyon vegetation
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,& ! output
                fp%hwratio,fp%effhwratio,fp%vangle,fp%hangle,fp%fbeam,                            & ! grid info
                cnveg%sigma,fp_road%alpha,cnveg%alpha,fp_wall%alpha, & ! grid info
                rdhyd%snowalpha,d_rdsndelta)



sg_walle = sg_walle*fp_coeffbldheight ! shadow due to in-canyon vegetation
sg_wallw = sg_wallw*fp_coeffbldheight ! shadow due to in-canyon vegetation



call getnetalbedo(u_alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,  &
                  fp%hwratio,fp%sigmabld,rfveg%sigma,fp_roof%alpha,rfveg%alpha,             &
                  cnveg%sigma,fp_road%alpha,fp_wall%alpha,cnveg%alpha,                      &
                  rfhyd%snowalpha,rdhyd%snowalpha,d_rfsndelta,d_rdsndelta)

sg_roof  = (1.-fp_roof%alpha)*sg_roof*a_sg
sg_vegr  = (1.-rfveg%alpha)*sg_vegr*a_sg
sg_walle = (1.-fp_wall%alpha)*sg_walle*a_sg
sg_wallw = (1.-fp_wall%alpha)*sg_wallw*a_sg
sg_road  = (1.-fp_road%alpha)*sg_road*a_sg
sg_vegc  = (1.-cnveg%alpha)*sg_vegc*a_sg
sg_rfsn  = (1.-rfhyd%snowalpha)*sg_rfsn*a_sg
sg_rdsn  = (1.-rdhyd%snowalpha)*sg_rdsn*a_sg



! calculate long wave reflections (pregenerated before canyonflux subroutine)
call getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta, &
                wallpsi,roadpsi,cnveg%sigma,fp_road%emiss,cnveg%emiss,fp_wall%emiss)

pd%emiss = d_rfsndelta*snowemiss+(1.-d_rfsndelta)*((1.-rfveg%sigma)*fp_roof%emiss+rfveg%sigma*rfveg%emiss)
pd%emiss = fp%sigmabld*pd%emiss+(1.-fp%sigmabld)*(2.*fp_wall%emiss*fp%effhwratio*d_cwa+d_netemiss*d_cra) ! diagnostic only

! Estimate the in-canyon and grid-cell urban roughness consistently with the
! coupling-side MOST profile calculation.
call uclem_canyon_roughness(zom, zonet, fp, cnveg, rdhyd, mp)

! Calculate distance from atmosphere to displacement height
d_rfdzmin = max(a_zroof-fp%bldheight,zoroof+0.2,rfveg%zo+0.2) ! distance to roof displacement height
pd%cndzmin = max(a_zmin-refheight*fp%bldheight,1.5,zom+0.2)   ! distance to canyon displacement height
pd%lzom    = log(pd%cndzmin/zom)

! 2
! calculate canyon wind speed and bulk transfer coefficients
! (i.e., acond = 1/(aerodynamic resistance) )
! some terms are updated when calculating canyon air temperature
select case(resmeth)
  case(0) ! Masson (2000)
    cu=exp(-0.25*fp%effhwratio)
    abase_road =cu ! bulk transfer coefficients are updated in canyonflux
    abase_walle=cu
    abase_wallw=cu
    abase_rdsn =cu
    abase_vegc =cu
  case(2) ! Kusaka et al (2001)
    cu=exp(-0.386*fp%effhwratio)
    abase_road =cu ! bulk transfer coefficients are updated in canyonflux
    abase_walle=cu
    abase_wallw=cu
    abase_rdsn =cu
    abase_vegc =cu
  case(3) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    call getincanwindb(we,ww,wr,a_udir,zonet,fp,mp)

    dis=max(0.1*fp_coeffbldheight*fp%bldheight,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_walle=a*we                 ! east wall bulk transfer
    abase_wallw=a*ww                 ! west wall bulk transfer
    dis=max(0.1*fp_coeffbldheight*fp%bldheight,zocanyon+0.2,cnveg%zo+0.2,zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_road=a*wr                  ! road bulk transfer
    zolog=log(dis/cnveg%zo)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_rdsn=a*wr                  ! road snow bulk transfer
  case default
    !write(6,*) "ERROR: Invalid choice for resmeth ",resmeth
    stop
end select


! 3
call getdiurnal(fp%ctime,fp%weekdayload,cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation)

! traffic sensible heat flux
pd%traf = fp%trafficfg*cyc_traffic
d_traf = pd%traf/(1.-fp%sigmabld)
! internal gains sensible heat flux
d_intgains_bld = (fp%intmassn+1.)*fp%intgains_flr*cyc_basedemand ! building internal gains
pd%intgains_full= fp%sigmabld*d_intgains_bld                     ! full domain internal gains


! 4
! calculate canyon fluxes
call solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,abase_road,                          &
                 sg_walle,rg_walle,fg_walle,acond_walle,abase_walle,                             &
                 sg_wallw,rg_wallw,fg_wallw,acond_wallw,abase_wallw,                             &
                 sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,abase_vegc,                          &
                 sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn, &
                 a_umag,a_rho,a_rg,a_rnd,a_snd,we,ww,                                            &
                 d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,                &
                 d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf,d_ac_inside,             &
                 d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,              &
                 d_cwr,d_totdepth,d_intgains_bld,d_roomstor,fgtop,egtop,evspsbltop,sbltop,       &
                 int_infilflux,int_infilfg,cyc_proportion,cyc_translation,ddt,                   &
                 cnveg,fp,fp_intm,fp_road,fp_roof,fp_wall,fp_slab,intm,pd,rdhyd,rfveg,road,      &
                 roof,room,slab,walle,wallw,iroomtemp,cvcoeff_roof,cvcoeff_walle,                &
                 cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,mp,diag)

zero_flux = 0.

! 5
! calculate roof fluxes (fg_roof updated below)
eg_roof = 0. ! For cray compiler
call solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
               sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
               sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
               a_rg,a_umagroof,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,     &
               d_sigr,ddt,fp_roof,rfhyd,rfveg,roof,fp,mp)


! update interior fluxes (first estimate)
if ( intairtmeth==0 .or. lwintmeth==0 ) then
  rgint_slab  = 0.
  rgint_wallw = 0.
  rgint_roof  = 0.
  rgint_walle = 0.
else if ( lwintmeth==1 .or. progcalc>0 ) then
  iskintemp(:,1) = real( slab%nodetemp(:,nl) )
  iskintemp(:,2) = real( wallw%nodetemp(:,nl) )
  iskintemp(:,3) = real( roof%nodetemp(:,nl) )
  iskintemp(:,4) = real( walle%nodetemp(:,nl) )
  ! calculate interior lwfluxes using last timestep temps lw (rg) fluxes
  call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle,    &
                       iskintemp,fp,intl,mp)
end if


! 6 only needed in the impilicit call
!------------------------------------------------------------------------------
! Update prognostic temperatures
!------------------------------------------------------------------------------

if ( progcalc==0 ) then

  ggext_roof = (1.-d_rfsndelta)*(sg_roof+rg_roof-eg_roof+aircp*a_rho*d_tempr*acond_roof) &
                +d_rfsndelta*garfsn
  ggext_walle= sg_walle+rg_walle+aircp*a_rho*d_canyontemp*acond_walle*fp_coeffbldheight
  ggext_wallw= sg_wallw+rg_wallw+aircp*a_rho*d_canyontemp*acond_wallw*fp_coeffbldheight
  ggext_road = (1.-d_rdsndelta)*(sg_road+rg_road-eg_road                        &
               +aircp*a_rho*d_canyontemp*acond_road) + d_rdsndelta*gardsn
  ggext_slab = 0.
  if ( intairtmeth/=0 ) then
    if ( intmassmeth/=0 ) then
      ggext_intm = cvcoeff_intm1*iroomtemp
    else
      ggext_intm = 0.
    end if
  else
    ggext_intm = 0.
  end if


  if ( intairtmeth/=0 .and. lwintmeth==2 ) then
    ! explicit longwave using update skin temperature from intermediate facet temp calc
    rgint_slab  = 0.
    rgint_wallw = 0.
    rgint_roof  = 0.
    rgint_walle = 0.

    !!!! slab !!!!
    ggint_impl = cvcoeff_slab
    ggint = -cvcoeff_slab*iroomtemp
    ggext_impl = 0.  ! zero_flux boundary condition
    call solvetridiag(ggext_slab,ggint,rgint_slab,ggext_impl,ggint_impl,slab%nodetemp,  &
                      inodetemp,fp_slab%depth,fp_slab%volcp,fp_slab%lambda,ddt,mp)
    iskintemp(:,1) = real(inodetemp(:,nl))

    !!!! wallw !!!!
    ggint_impl = cvcoeff_wallw
    ggint = -cvcoeff_wallw*iroomtemp
    ggext_impl = aircp*a_rho*acond_wallw*fp_coeffbldheight ! implicit, later update fg_wallw with final T
    call solvetridiag(ggext_wallw,ggint,rgint_wallw,ggext_impl,ggint_impl,wallw%nodetemp,  &
                      inodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,mp)
    iskintemp(:,2) = real(inodetemp(:,nl))

    !!!! roof !!!!
    ggint_impl = cvcoeff_roof
    ggint = -cvcoeff_roof*iroomtemp
    ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! implicit, later update fg_roof with final T
    call solvetridiag(ggext_roof,ggint,rgint_roof,ggext_impl,ggint_impl,roof%nodetemp,  &
                      inodetemp,fp_roof%depth,fp_roof%volcp,fp_roof%lambda,ddt,mp)
    iskintemp(:,3) = real(inodetemp(:,nl))

    !!!! walle !!!!
    ggint_impl = cvcoeff_walle
    ggint = -cvcoeff_walle*iroomtemp
    ggext_impl = aircp*a_rho*acond_walle*fp_coeffbldheight ! implicit, later update fg_walle with final T
    call solvetridiag(ggext_walle,ggint,rgint_walle,ggext_impl,ggint_impl,walle%nodetemp,  &
                      inodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,mp)
    iskintemp(:,4) = real(inodetemp(:,nl))

    ! calculate interior lwfluxes using iskintemp
    call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle,    &
                         iskintemp,fp,intl,mp)
  end if ! intairtmeth/=0 .and. lwintmeth==2


  ! solve final conduction through facets with explicit interior, implicit exterior boundary
  ggint_impl = cvcoeff_roof
  ggint = -cvcoeff_roof*iroomtemp
  ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! later update fg_roof with final roof skin T
  call solvetridiag(ggext_roof,ggint,rgint_roof,ggext_impl,ggint_impl,roof%nodetemp,     &
                    roof%nodetemp,fp_roof%depth,fp_roof%volcp,fp_roof%lambda,ddt,mp)

  ggint_impl = cvcoeff_walle
  ggint = -cvcoeff_walle*iroomtemp
  ggext_impl = aircp*a_rho*acond_walle*fp_coeffbldheight ! later update fg_walle with final walle skin T
  call solvetridiag(ggext_walle,ggint,rgint_walle,ggext_impl,ggint_impl,walle%nodetemp,  &
                    walle%nodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,mp)

  ggint_impl = cvcoeff_wallw
  ggint = -cvcoeff_wallw*iroomtemp
  ggext_impl = aircp*a_rho*acond_wallw*fp_coeffbldheight ! later update fg_wallw with final wallw skin T
  call solvetridiag(ggext_wallw,ggint,rgint_wallw,ggext_impl,ggint_impl,wallw%nodetemp,  &
                    wallw%nodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,mp)

  ggint_impl = 0.  ! zero_flux boundary condition
  ggint_road = 0.  ! zero_flux boundary condition
  ggext_impl = (1.-d_rdsndelta)*aircp*a_rho*acond_road ! later update fg_road with final road skin T
  call solvetridiag(ggext_road,ggint_road,zero_flux,ggext_impl,ggint_impl,road%nodetemp,   &
                    road%nodetemp,fp_road%depth,fp_road%volcp,fp_road%lambda,ddt,mp)

  ! calculate internal facet conduction
  if ( intairtmeth==0 ) then
    slab%nodetemp = 1._8 ! + urbtemp
  else
    ggint_impl = cvcoeff_slab
    ggint = -cvcoeff_slab*iroomtemp
    ggext_slab = 0.  ! zero_flux boundary condition
    ggext_impl = 0.  ! zero_flux boundary condition
    call solvetridiag(ggext_slab,ggint,rgint_slab,ggext_impl,ggint_impl,slab%nodetemp,    &
                      slab%nodetemp,fp_slab%depth,fp_slab%volcp,fp_slab%lambda,ddt,mp)
    if ( intmassmeth/=0 ) then
      ggint_impl = cvcoeff_intm2
      ggint = -cvcoeff_intm2*iroomtemp
      ggext_impl = cvcoeff_intm1
      call solvetridiag(ggext_intm,ggint,zero_flux,ggext_impl,ggint_impl,intm%nodetemp,    &
                        intm%nodetemp,fp_intm%depth,fp_intm%volcp,fp_intm%lambda,ddt,mp)
    end if
  end if

end if ! if progcalc==0 ...


! 7 everything before energy closure
!------------------------------------------------------------------------------
! Update sensible heat fluxes
!------------------------------------------------------------------------------

fg_roof  = aircp*a_rho*(real(roof%nodetemp(:,0))-d_tempr)*acond_roof
fg_walle = aircp*a_rho*(real(walle%nodetemp(:,0))-d_canyontemp)*acond_walle*fp_coeffbldheight
fg_wallw = aircp*a_rho*(real(wallw%nodetemp(:,0))-d_canyontemp)*acond_wallw*fp_coeffbldheight
fg_road  = aircp*a_rho*(real(road%nodetemp(:,0))-d_canyontemp)*acond_road
if ( intairtmeth/=0 ) then
  if ( intmassmeth/=0 ) then
    ggext_intm = -cvcoeff_intm1*(real(intm%nodetemp(:,0))-iroomtemp)
  else
    ggext_intm = 0.
  end if
else
  ggext_intm = 0.
end if

ggint_roof  = cvcoeff_roof*(real(roof%nodetemp(:,nl))-iroomtemp)
ggint_walle = cvcoeff_walle*(real(walle%nodetemp(:,nl))-iroomtemp)
ggint_wallw = cvcoeff_wallw*(real(wallw%nodetemp(:,nl))-iroomtemp)
if ( intairtmeth/=0 ) then
  ggint_slab = cvcoeff_slab*(real(slab%nodetemp(:,nl))-iroomtemp)
  if ( intmassmeth/=0 ) then
    ggint_intm = cvcoeff_intm2*(real(intm%nodetemp(:,nl))-iroomtemp)
  else
    ggint_intm = 0.
  end if
else
  ggint_slab = 0.
  ggint_intm = 0.
end if

d_ac_inside = -ggint_roof - ggint_slab + (ggext_intm - ggint_intm)*real(fp%intmassn) &
              - (ggint_walle+ggint_wallw)*(fp%bldheight/fp%bldwidth)                 &
              - int_infilflux - d_intgains_bld + d_roomstor


where ( d_ac_inside<0. )
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
  ac_coeff = max( acfactor*(d_canyontemp+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
elsewhere
  ac_coeff = 0.
end where


select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, COP of 1.0)
  case(0) ! unrealistic cooling (buildings act as heat sink)
    d_ac_canyon = 0.
    pd%bldheat = max(0.,d_ac_inside*fp%sigmabld)
    pd%bldcool = max(0.,d_ac_inside*fp%sigmabld)
  case(1) ! d_ac_canyon pumps conducted heat + ac waste heat back into canyon
    d_ac_canyon = max(0.,-d_ac_inside*(1.+ac_coeff)*fp%sigmabld/(1.-fp%sigmabld))  ! canyon domain W/m/m
    pd%bldheat = max(0.,d_ac_inside*fp%sigmabld)            ! entire domain W/m/m
    pd%bldcool = max(0.,-d_ac_inside*ac_coeff*fp%sigmabld)  ! entire domain W/m/m
  case(2) ! reversible heating and cooling (for testing energy conservation)
    d_ac_canyon = -d_ac_inside*fp%sigmabld/(1.-fp%sigmabld)
    pd%bldheat = 0.
    pd%bldcool = 0.
  case DEFAULT
    write(6,*) "ERROR: Unknown acmeth mode ",acmeth
    stop
end select

qf = pd%bldheat + pd%bldcool + pd%traf + fp%industryfg + pd%intgains_full  ! complete AHF = heating + cooling(work) + traffic + industry + internal gains (matches atmos energy balance, line ~1267)


! update canyon flux
fgtop = fp%hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_canyon - int_infilfg

!------------------------------------------------------------------------------
! Road hydrology calculations
!------------------------------------------------------------------------------

! reset inception of leaf water
pd%delintercept = 0.

d_surfwater = rdhyd%surfwater
d_soilwater = rdhyd%soilwater
d_leafwater = rdhyd%leafwater
d_snow = rdhyd%snow
d_den = rdhyd%den
d_snowalpha = rdhyd%snowalpha

! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  d_den = (d_snow*d_den+a_snd*ddt*minsnowden)/(d_snow+ddt*a_snd)
  ! reset snow albedo
  d_snowalpha = maxsnowalpha
end where

! calculate water/snow budgets for road surface
call updatewater(ddt,d_surfwater,d_soilwater,d_leafwater,d_snow,              &
                 d_den,d_snowalpha,rdsnmelt,a_rnd,a_snd,eg_road,              &
                 eg_rdsn,d_tranc,d_evapc,d_totdepth, cnveg%lai,wbrelaxc,      &
                 fp%sfc,fp%swilt,d_irrigwater,delintercept,rd_delsnow,mp)
! record grid irrigation road component [kg m-2 s-1]
pd%irrig = d_irrigwater*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)/ddt
! record leaf intercept component [kg m-2 s-1]
pd%delintercept = pd%delintercept + delintercept*cnveg%sigma

pd%road_water_runoff = max(d_surfwater-maxrdwater,0.)
pd%road_snow_runoff  = max(d_snow-maxrdsn,0.)
pd%road_soil_runoff  = max(d_soilwater-fp%ssat,0.)

!------------------------------------------------------------------------------
! Update road prognostic hydrology variables
!------------------------------------------------------------------------------

if ( progcalc==0 ) then
  rdhyd%surfwater(1:mp) = max(d_surfwater - pd%road_water_runoff,0.)
  rdhyd%soilwater(1:mp) = max(d_soilwater - pd%road_soil_runoff,fp%swilt)
  rdhyd%leafwater          = max(d_leafwater,0.)
  rdhyd%snow(1:mp)      = max(d_snow - pd%road_snow_runoff,0.)
  rdhyd%den                = min(max(d_den,minsnowden),maxsnowden)
  rdhyd%snowalpha          = min(max(d_snowalpha,minsnowalpha),maxsnowalpha)
end if

!------------------------------------------------------------------------------
! Roof hydrology calculations
!------------------------------------------------------------------------------

d_surfwater = rfhyd%surfwater
d_soilwater = rfhyd%soilwater
d_leafwater = rfhyd%leafwater
d_snow = rfhyd%snow
d_den = rfhyd%den
d_snowalpha = rfhyd%snowalpha

! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  d_den = (rfhyd%snow*rfhyd%den+a_snd*ddt*minsnowden)/(rfhyd%snow+ddt*a_snd)
  ! reset snow albedo
  d_snowalpha = maxsnowalpha
end where

! calculate water/snow budgets for roof surface
call updatewater(ddt,d_surfwater,d_soilwater,d_leafwater,d_snow,             &
                 d_den,d_snowalpha,rfsnmelt,a_rnd,a_snd,eg_roof,             &
                 eg_rfsn,d_tranr,d_evapr,fp%rfvegdepth,rfveg%lai,wbrelaxr,   &
                 fp%sfc,fp%swilt,d_irrigwater,delintercept,rf_delsnow,mp)
! add grid irrigation roof component [kg m-2 s-1]
pd%irrig = pd%irrig + d_irrigwater*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld/ddt
! record leaf intercept component [kg m-2 s-1]
pd%delintercept = pd%delintercept + delintercept*rfveg%sigma

pd%roof_water_runoff = max(d_surfwater-maxrfwater,0.)
pd%roof_snow_runoff  = max(d_snow-maxrfsn,0.)
pd%roof_soil_runoff  = max(d_soilwater-fp%ssat,0.)

!------------------------------------------------------------------------------
! Update roof prognostic hydrology variables
!------------------------------------------------------------------------------

if ( progcalc==0 ) then
  rfhyd%surfwater(1:mp) = max(d_surfwater - pd%roof_water_runoff,0.)
  rfhyd%soilwater(1:mp) = max(d_soilwater - pd%roof_soil_runoff,fp%swilt)
  rfhyd%leafwater          = max(d_leafwater,0.)
  rfhyd%snow(1:mp)      = max(d_snow - pd%roof_snow_runoff,0.)
  rfhyd%den                = min(max(d_den,minsnowden),maxsnowden)
  rfhyd%snowalpha          = min(max(d_snowalpha,minsnowalpha),maxsnowalpha)
end if

!------------------------------------------------------------------------------
! Calculate pd diagnostic and u_? output variables
!------------------------------------------------------------------------------

! area weighted vegetation aerodynamic conductance
pd%acond_vegw = (acond_vegc*cnveg%sigma*(1.-fp%sigmabld) + acond_vegr*rfveg%sigma*fp%sigmabld)/  &
               (rfveg%sigma*fp%sigmabld + cnveg%sigma*(1.-fp%sigmabld))

! calculate runoff (leafwater runoff already accounted for in precip reaching canyon floor) [kg m-2]
u_rn = pd%roof_water_runoff*fp%sigmabld*(1.-rfveg%sigma)                     & ! roof surfwater runnoff
      +pd%road_water_runoff*(1.-fp%sigmabld)*(1.-cnveg%sigma)                & ! road surfwater runnoff
      +pd%roof_snow_runoff*fp%sigmabld                                       & ! roof snow removal
      +pd%road_snow_runoff*(1.-fp%sigmabld)                                  & ! road snow removal
      +pd%roof_soil_runoff*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld    & ! roofveg saturation runnoff
      +pd%road_soil_runoff*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)    ! canyveg saturation runnoff

! diagnose surface drainage rate [kg m-2 s-1]
pd%surfrunoff = u_rn/ddt

! diagnose total vegetation transpiration [kg m-2 s-1]
pd%transveg = (d_tranc*cnveg%sigma*(1.-fp%sigmabld) + d_tranr*rfveg%sigma*fp%sigmabld)/lv

! combine snow and snow-free tiles for fluxes
d_roofrgout = a_rg-d_rfsndelta*rg_rfsn-(1.-d_rfsndelta)*((1.-rfveg%sigma)*rg_roof+rfveg%sigma*rg_vegr)
fgrooftop   = d_rfsndelta*fg_rfsn+(1.-d_rfsndelta)*((1.-rfveg%sigma)*fg_roof+rfveg%sigma*fg_vegr)
egrooftop   = d_rfsndelta*eg_rfsn+(1.-d_rfsndelta)*((1.-rfveg%sigma)*eg_roof+rfveg%sigma*eg_vegr)
evspsblrooftop = d_rfsndelta*eg_rfsn/ls+(1.-d_rfsndelta)*((1.-rfveg%sigma)*eg_roof+rfveg%sigma*eg_vegr)/lv
sblrooftop  = d_rfsndelta*eg_rfsn/ls



!fgtop       = d_rdsndelta*fg_rdsn+(1.-d_rdsndelta)*((1.-cnveg%sigma)*fg_road+cnveg%sigma*fg_vegc)   &
!             +fp%hwratio*(fg_walle+fg_wallw)+d_traf+d_ac_canyon
!egtop       = d_rdsndelta*eg_rdsn+(1.-d_rdsndelta)*((1.-cnveg%sigma)*eg_road+cnveg%sigma*eg_vegc)

! calculate wetfac for roof and road vegetation (see sflux.f or cable_canopy.f90)
roofvegwetfac = max(min((rfhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)
roadvegwetfac = max(min((rdhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)

! area weighted soil water (canyon and roof)
pd%soilwater = (rdhyd%soilwater*cnveg%sigma*(1.-fp%sigmabld) + rfhyd%soilwater*rfveg%sigma*fp%sigmabld) / &
               max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)

! record previous soilmoist
tmp = pd%soilmoist
! update soilmoist only in areas which support vegetation (i.e. have soil moisture)
pd%soilmoist = ( rdhyd%soilwater*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)   &  ! canyon integrated soil moisture
               + rfhyd%soilwater*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld ) / &  ! roof integrated soil moisture
               max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)     ! limited to areas with soil
! calculate difference in soilmoist, weighted over entire grid (i.e. include dry area under roads + buildings)
if ( pd%first_call ) then
  pd%delsoilmoist = 0.
else
  pd%delsoilmoist = (pd%soilmoist - tmp)*(rfveg%sigma*fp%sigmabld + cnveg%sigma*(1.-fp%sigmabld))
end if

! rootmoist calculated only for canyon veg subgrid
pd%rootmoistc = rdhyd%soilwater*waterden*d_totdepth

! soilwetness: (soil moisture - wilting) / (saturation - wilting).
tmp = (rdhyd%soilwater*cnveg%sigma*(1.-fp%sigmabld) + rfhyd%soilwater*rfveg%sigma*fp%sigmabld) / &
      max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)
pd%soilwetness = (tmp-fp%swilt)/(fp%ssat-fp%swilt)

! calculate longwave, sensible heat latent heat outputs
! estimate surface temp from outgoing longwave radiation
u_ts = ((fp%sigmabld*d_roofrgout+(1.-fp%sigmabld)*d_canyonrgout)/sbconst)**0.25 - urbtemp

u_fg = fp%sigmabld*fgrooftop+(1.-fp%sigmabld)*fgtop+fp%industryfg

u_eg = fp%sigmabld*egrooftop+(1.-fp%sigmabld)*egtop

u_evspsbl = fp%sigmabld*evspsblrooftop + (1.-fp%sigmabld)*evspsbltop

u_sbl = fp%sigmabld*sblrooftop + (1.-fp%sigmabld)*sbltop

u_temp9(:,1) = real(roof%nodetemp(:,0))
u_temp9(:,2) = real(road%nodetemp(:,0))
u_temp9(:,3) = real(walle%nodetemp(:,0))
u_temp9(:,4) = real(wallw%nodetemp(:,0))
u_temp9(:,5) = real(slab%nodetemp(:,0))
u_temp9(:,6) = real(intm%nodetemp(:,0))
u_temp9(:,7) = real(room%nodetemp(:,1))
u_temp9(:,8) = cnveg%temp
u_temp9(:,9) = rfveg%temp

! 1 fg_roof, 2 fg_road, 3 fg_walle, 4 fg_wallw, 5 fg_cnveg, 6 fg_rfveg, 7 ac/other
! Area-weight components so the sum matches u_fg (including roof/road snow, traffic,
! infiltration, and industry heat).
u_shf9(:,1) = fp%sigmabld * ( d_rfsndelta*fg_rfsn + (1.-d_rfsndelta)*(1.-rfveg%sigma)*fg_roof )
u_shf9(:,2) = (1.-fp%sigmabld) * ( d_rdsndelta*fg_rdsn + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road )
u_shf9(:,3) = (1.-fp%sigmabld) * fp%hwratio * fg_walle
u_shf9(:,4) = (1.-fp%sigmabld) * fp%hwratio * fg_wallw
u_shf9(:,5) = (1.-fp%sigmabld) * (1.-d_rdsndelta) * cnveg%sigma * fg_vegc
u_shf9(:,6) = fp%sigmabld * (1.-d_rfsndelta) * rfveg%sigma * fg_vegr
u_shf9(:,7) = (1.-fp%sigmabld) * (d_ac_canyon + d_traf - int_infilfg) + fp%industryfg
u_shf9(:,8) = pd%bldheat   ! building heating load
u_shf9(:,9) = pd%bldcool   ! building cooling (compressor work)

u_wf = fp%sigmabld*(1.-d_rfsndelta)*((1.-rfveg%sigma)*d_roofdelta       &
      +rfveg%sigma*((1.-d_vegdeltar)*roofvegwetfac+d_vegdeltar))        &
      +(1.-fp%sigmabld)*(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta  &
      +cnveg%sigma*((1.-d_vegdeltac)*roadvegwetfac+d_vegdeltac))



pd%snowmelt = fp%sigmabld*rfsnmelt + (1.-fp%sigmabld)*rdsnmelt
u_melt = lf*(fp%sigmabld*d_rfsndelta*rfsnmelt + (1.-fp%sigmabld)*d_rdsndelta*rdsnmelt)

! calculate average snow temperature and Snow water equivilent
pd%snowt = ( (rfsntemp+urbtemp)*d_rfsndelta*fp%sigmabld + (rdsntemp+urbtemp)*d_rdsndelta*(1.-fp%sigmabld) ) / &
           max( d_rfsndelta*fp%sigmabld + d_rdsndelta*(1.-fp%sigmabld), 1.e-10 )
pd%swe = rfhyd%snow*d_rfsndelta*fp%sigmabld + rdhyd%snow*d_rdsndelta*(1.-fp%sigmabld)
pd%snowfrac = d_rfsndelta*fp%sigmabld + d_rdsndelta*(1.-fp%sigmabld)
pd%salbedo = ( rfhyd%snowalpha*d_rfsndelta*fp%sigmabld + rdhyd%snowalpha*d_rdsndelta*(1.-fp%sigmabld) ) / &
           max( d_rfsndelta*fp%sigmabld + d_rdsndelta*(1.-fp%sigmabld), 1.e-10 )
pd%delswe = rf_delsnow*d_rfsndelta*fp%sigmabld + rd_delsnow*d_rdsndelta*(1.-fp%sigmabld)

! calculate total liquid water storage
! road term: non-building fraction, masked by road snow and canyon veg
! roof term: building fraction, masked by roof snow and roof veg
pd%surfstor = rdhyd%surfwater*(1.-fp%sigmabld)*(1.-d_rdsndelta)*(1.-cnveg%sigma) &
            + rfhyd%surfwater*fp%sigmabld*(1.-d_rfsndelta)*(1.-rfveg%sigma)

! calculate average vegetation temperature
pd%vegt = ( (cnveg%temp+urbtemp)*(1.-fp%sigmabld)*cnveg%sigma + (rfveg%temp+urbtemp)*fp%sigmabld*rfveg%sigma ) / &
          max( (1.-fp%sigmabld)*cnveg%sigma + fp%sigmabld*rfveg%sigma, 1.e-10 )
pd%calbedo = ( cnveg%alpha*(1.-fp%sigmabld)*cnveg%sigma + rfveg%alpha*fp%sigmabld*rfveg%sigma ) / &
          max( (1.-fp%sigmabld)*cnveg%sigma + fp%sigmabld*rfveg%sigma, 1.e-10 )
pd%ulai = ( cnveg%lai*(1.-fp%sigmabld)*cnveg%sigma + rfveg%lai*fp%sigmabld*rfveg%sigma ) / &
          max( (1.-fp%sigmabld)*cnveg%sigma + fp%sigmabld*rfveg%sigma, 1.e-10 )

! Save canyon air temperature
pd%taircanyon = d_canyontemp




! (re)calculate heat roughness length for MOST (diagnostic only)
call getqsat(a,u_ts,d_sigd)
a = a*u_wf + min(a_mixr,a)*(1.-u_wf)
dts = u_ts + (u_ts+urbtemp)*0.61*a
dtt = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
select case(zohmeth)
  case(0) ! Use veg formulation
    pd%lzoh = 2.3+pd%lzom
    call getinvres(pd%cdtq,pd%cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,1)
  case(1) ! Use Kanda parameterisation
    pd%lzoh = 2.3+pd%lzom ! replaced in getlna
    call getinvres(pd%cdtq,pd%cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,2)
  case(2) ! Use Kanda parameterisation
    pd%lzoh = 6.+pd%lzom
    call getinvres(pd%cdtq,pd%cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,4)
end select

! z_on_l is z/L at the UCLEM canyon reference height.  Retain it so the
! UM-facing screen/10-m interpolation uses the same stability as the UCLEM
! roughness and exchange coefficients.
pd%zol = z_on_l

! calculate screen level diagnostics
call scrncalc(a_mixr,a_umag,c_a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,a,rdsntemp,zonet, &
              cnveg,fp,pd,rdhyd,road,mp)


! 8 energy closure
!------------------------------------------------------------------------------
! Check energy conservation
!------------------------------------------------------------------------------
if ( progcalc==0 ) then

  call energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                     sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                     rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                     eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                     a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,                  &
                     ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                     ggext_intm,ggint_slab,ggint_intm,d_intgains_bld,        &
                     int_infilflux,d_ac_inside,fp,ddt,                       &
                     cnveg,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,  &
                     rfveg,road,roof,room,slab,walle,wallw,mp)
end if

  ! Save longwave radiation parameters
  pd%rg_roof = rg_roof
  pd%rg_road = rg_road
  pd%rg_walle = rg_walle
  pd%rg_wallw = rg_wallw
  pd%rg_vegc = rg_vegc
  pd%rg_vegr = rg_vegr
  pd%rg_rfsn = rg_rfsn
  pd%rg_rdsn = rg_rdsn
  pd%d_canyonrgout = d_canyonrgout
  pd%d_roofrgout = d_roofrgout

  ! save shortwave radiation parameters
  pd%sg_roof = sg_roof
  pd%sg_road = sg_road
  pd%sg_walle = sg_walle
  pd%sg_wallw = sg_wallw
  pd%sg_vegc = sg_vegc
  pd%sg_vegr = sg_vegr
  pd%sg_rfsn = sg_rfsn
  pd%sg_rdsn = sg_rdsn


  ! save u_fg budget terms


  pd%acond_roof = acond_roof
  pd%acond_road = acond_road
  pd%acond_walle = acond_walle
  pd%acond_wallw = acond_wallw
  pd%d_tempr = d_tempr
  pd%d_canyontemp = d_canyontemp
  pd%acond_vegc = acond_vegc
  pd%acond_vegr = acond_vegr
  pd%acond_rfsn = acond_rfsn
  pd%acond_rdsn = acond_rdsn

  pd%d_rfsndelta = d_rfsndelta
  pd%d_rdsndelta = d_rdsndelta
  pd%rfveg_sigma = rfveg%sigma
  pd%cnveg_sigma = cnveg%sigma
  pd%rdsntemp = rdsntemp
  pd%rfsntemp = rfsntemp

  pd%d_ac_inside = d_ac_inside
  pd%int_infilflux = int_infilflux
  pd%ac_coeff = ac_coeff

  ! eg
  pd%eg_roof = eg_roof
  pd%eg_vegr = eg_vegr
  pd%eg_road = eg_road
  pd%eg_vegc = eg_vegc
  pd%eg_rdsn = eg_rdsn
  pd%eg_rfsn = eg_rfsn

  ! rn
  pd%d_totdepth = d_totdepth

return
end subroutine uclem_eval


subroutine uclem_canyon_roughness(zom, zonet, fp, cnveg, rdhyd, mp)
! Calculate the urban momentum roughness used by both UCLEM and the
! coupling-side MOST profile.  Keeping this calculation in one place avoids
! applying a nominal 0.1H roughness to the forcing while exporting a different
! snow- and vegetation-blended roughness to the UM.

implicit none

integer, intent(in) :: mp
real, dimension(mp), intent(out) :: zom, zonet
type(fparmdata), intent(in) :: fp
type(vegdata), intent(in) :: cnveg
type(hydrodata), intent(in) :: rdhyd

real, dimension(mp) :: dis, zolog, snowfrac, fp_coeffbldheight, road_snowfrac

fp_coeffbldheight = fp%effhwratio/fp%hwratio
road_snowfrac = rdhyd%snow/(rdhyd%snow+maxrdsn)
dis = max(max(max(0.1*fp_coeffbldheight*fp%bldheight,zocanyon+0.2),             &
      cnveg%zo+0.2), zosnow+0.2)
zolog = 1./sqrt(road_snowfrac/log(dis/zosnow)**2 +                              &
        (1.-road_snowfrac)*(cnveg%sigma/log(dis/cnveg%zo)**2 +                  &
        (1.-cnveg%sigma)/log(dis/zocanyon)**2))
zonet = dis*exp(-zolog)

zom = zomratio*fp%bldheight
where (zom*fp%sigmabld < zonet*(1.-fp%sigmabld))
  zom = zonet
end where
snowfrac = rdhyd%snow/(rdhyd%snow+maxrdsn+0.408*grav*zom)
zom = (1.-snowfrac)*zom + snowfrac*zosnow

end subroutine uclem_canyon_roughness


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tridiagonal solver for temperatures

! This version has an implicit estimate for roof sensible heat flux

! [ ggB ggC         ] [ temp ] = [ ggD ]
! [ ggA ggB ggC     ] [ temp ] = [ ggD ]
! [     ggA ggB ggC ] [ temp ] = [ ggD ]
! [         ggA ggB ] [ temp ] = [ ggD ]

subroutine solvetridiag(ggext,ggint,rgint,ggext_impl,ggint_impl,    &
                        nodetemp,inodetemp,depth,volcp,lambda,ddt,mp)

implicit none

integer, intent(in) :: mp
real, dimension(mp),     intent(in)    :: ggext,ggint,rgint  ! surface energy fluxes
real, dimension(mp),     intent(in)    :: ggext_impl         ! implicit update for stability
real, dimension(mp),     intent(in)    :: ggint_impl         ! implicit update for stability
real(kind=8), dimension(mp,0:nl),intent(in)  :: nodetemp     ! temperature of each node
real(kind=8), dimension(mp,0:nl),intent(out) :: inodetemp    ! updated temperature solution
real, dimension(mp,nl),  intent(in)    :: depth,volcp,lambda ! facet depth, heat capacity, conductivity
real(kind=8), dimension(mp,nl)         :: cap,res            ! layer capacitance & resistance
real(kind=8), dimension(mp,0:nl)       :: ggA,ggB,ggC,ggD    ! tridiagonal matrices
real(kind=8) ggX                                                ! tridiagonal coefficient
real(kind=8), dimension(mp)            :: ans                ! tridiagonal solution
real, intent(in)                          :: ddt                ! timestep
integer k, iqu

! Add numerical stability checks
res = real(depth,8)/max(real(lambda,8), 1.0e-10_8)  ! Prevent division by zero
cap = real(depth,8)*real(volcp,8)


!!!!!!!!! interface conduction !!!!!!!!!!!
ggA(:,1:nl)   = -1._8/res(:,1:nl)
ggB(:,0)      =  1._8/res(:,1) +0.5_8*cap(:,1)/real(ddt,8) + real(ggext_impl,8)
ggB(:,1:nl-1) =  1._8/res(:,1:nl-1) +1._8/res(:,2:nl) +0.5_8*(cap(:,1:nl-1) &
                   +cap(:,2:nl))/real(ddt,8)
ggB(:,nl)     =  1._8/res(:,nl) + 0.5_8*cap(:,nl)/real(ddt,8) + real(ggint_impl,8)
ggC(:,0:nl-1) = -1._8/res(:,1:nl)
ggD(:,0)      = nodetemp(:,0)*0.5_8*cap(:,1)/real(ddt,8) + real(ggext,8)
ggD(:,1:nl-1) = nodetemp(:,1:nl-1)*0.5_8*(cap(:,1:nl-1)+cap(:,2:nl))/real(ddt,8)
ggD(:,nl)     = nodetemp(:,nl)*0.5_8*cap(:,nl)/real(ddt,8) - real(ggint,8) - real(rgint,8)
! tridiagonal solver (Thomas algorithm) to solve node temperatures
do k = 1,nl
  do iqu = 1,mp
    ggX = ggA(iqu,k)/ggB(iqu,k-1)
    ggB(iqu,k) = ggB(iqu,k)-ggX*ggC(iqu,k-1)
    ggD(iqu,k) = ggD(iqu,k)-ggX*ggD(iqu,k-1)
  end do
end do
! Check for zero diagonal elements in back substitution
where (abs(ggB(:,nl)) < 1.0e-15_8)
  ggB(:,nl) = sign(1.0e-15_8, ggB(:,nl))
end where

ans = ggD(:,nl)/ggB(:,nl)
inodetemp(:,nl) = ans
do k = nl-1,0,-1
  where (abs(ggB(:,k)) < 1.0e-15_8)
    ggB(:,k) = sign(1.0e-15_8, ggB(:,k))
  end where
  ans = (ggD(:,k) - ggC(:,k)*ans)/ggB(:,k)
  inodetemp(:,k) = ans
end do

end subroutine solvetridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conservation of energy check

subroutine energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                         sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                         rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                         eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                         a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,                  &
                         ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                         ggext_intm,ggint_slab,ggint_intm,d_intgains_bld,      &
                         int_infilflux,d_ac_inside,fp,ddt,cnveg,fp_intm,fp_road, &
                         fp_roof,fp_slab,fp_wall,intm,pd,rfveg,road,roof,room,   &
                         slab,walle,wallw,mp)

implicit none

integer, intent(in) :: mp
real, intent(in) :: ddt
real, dimension(mp), intent(in) :: sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle
real, dimension(mp), intent(in) :: sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw
real, dimension(mp), intent(in) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab
real, dimension(mp), intent(in) :: eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta
real, dimension(mp), intent(in) :: a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt
real, dimension(mp), intent(in) :: ggint_roof,ggint_road,ggint_walle,ggint_wallw
real, dimension(mp), intent(in) :: ggext_intm,ggint_slab,ggint_intm,d_intgains_bld
real, dimension(mp), intent(in) :: int_infilflux,d_ac_inside
real(kind=8), dimension(mp) :: d_roofflux,d_walleflux,d_wallwflux,d_roadflux,d_slabflux,d_intmflux,d_roomflux
real(kind=8), dimension(mp) :: d_roofstor,d_wallestor,d_wallwstor,d_roadstor,d_slabstor,d_intmstor,d_roomstor
real(kind=8), dimension(mp) :: d_faceterr
real(kind=8), dimension(mp) :: d_atmosflux
real(kind=8), dimension(mp,nl) :: roadstorage_prev, roofstorage_prev, wallestorage_prev, wallwstorage_prev
real(kind=8), dimension(mp,nl) :: slabstorage_prev, intmstorage_prev
real(kind=8), dimension(mp,1) :: roomstorage_prev
logical, dimension(mp) :: is_urban
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(facetdata), intent(inout) :: intm
type(vegdata), intent(in) :: cnveg, rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

! In this host UCLEM is evaluated on the whole CABLE mp vector rather than on a
! packed list of urban points, so every column carries a solution, but only the
! columns with a non-zero urban fraction carry a meaningful one.  The rest run
! with default urban parameters under whatever forcing their gridbox supplies
! and will routinely miss the closure tolerance.  Restricting the checks below
! to urban columns keeps them usable: without the mask they fire on almost every
! timestep and maxval reports a non-urban column, which hides any real urban
! closure failure.  fp%sigmau is already used as the urban point mask in
! cbl_uclem_coupling_mod.
is_urban = fp%sigmau > 0.

! Store previous calculation to determine flux
roofstorage_prev(:,:)  = roof%storage(:,:)
roadstorage_prev(:,:)  = road%storage(:,:)
wallestorage_prev(:,:) = walle%storage(:,:)
wallwstorage_prev(:,:) = wallw%storage(:,:)
slabstorage_prev(:,:)  = slab%storage(:,:)
intmstorage_prev(:,:)  = intm%storage(:,:)
roomstorage_prev(:,:)  = room%storage(:,:)
pd%surferr = 0.


! assume a_rho=1 for interior
room%storage(:,1) = real(fp%bldheight(:),8)*real(aircp,8)*room%nodetemp(:,1)
! Sum heat stored in urban materials from layer 1 to nl
! interface conduction
roof%storage(:,:)  = 0.5_8*real(fp_roof%depth(:,:),8)*real(fp_roof%volcp(:,:),8)                        &
                        *(roof%nodetemp(:,0:nl-1)+roof%nodetemp(:,1:nl))
road%storage(:,:)  = 0.5_8*real(fp_road%depth(:,:),8)*real(fp_road%volcp(:,:),8)                        &
                        *(road%nodetemp(:,0:nl-1)+road%nodetemp(:,1:nl))
walle%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8)                        &
                        *(walle%nodetemp(:,0:nl-1)+walle%nodetemp(:,1:nl))
wallw%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8)                        &
                        *(wallw%nodetemp(:,0:nl-1)+wallw%nodetemp(:,1:nl))
slab%storage(:,:)  = 0.5_8*real(fp_slab%depth(:,:),8)*real(fp_slab%volcp(:,:),8)                        &
                        *(slab%nodetemp(:,0:nl-1)+slab%nodetemp(:,1:nl))
intm%storage(:,:)  = 0.5_8*real(fp_intm%depth(:,:),8)*real(fp_intm%volcp(:,:),8)                        &
                        *(intm%nodetemp(:,0:nl-1)+intm%nodetemp(:,1:nl))

! u_alb and u_melt are copies of intent(in) arguments and have nothing to do with
! the comparison against the previous timestep, so they are stored before the
! first-call return.  Leaving them below it left pd%u_alb at its initialised zero
! for the first prognostic step of every run and every restart, which the host
! now sees because the urban tile takes its net radiation from pd%u_alb.
pd%u_alb  = u_alb
pd%u_melt = u_melt

! first call do not calculate comparison with last timestep (return)
if ( pd%first_call ) then
  pd%first_call = .false.
  return
end if

d_roofstor = sum(roof%storage-roofstorage_prev,dim=2)/real(ddt,8)
d_roofflux = (1._8-real(d_rfsndelta,8))*(real(sg_roof,8)+real(rg_roof,8)-real(fg_roof,8)-real(eg_roof,8))  &
           + real(d_rfsndelta,8)*real(garfsn,8) - real(ggint_roof,8) - real(rgint_roof,8)
d_faceterr  = d_roofstor - d_roofflux
pd%surferr = pd%surferr + d_faceterr
if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
  write(6,*) "aTEB roof facet closure error:", maxval(abs(d_faceterr),mask=is_urban)
d_roadstor = sum(road%storage-roadstorage_prev,dim=2)/real(ddt,8)
d_roadflux = (1._8-real(d_rdsndelta,8))*(real(sg_road,8)+real(rg_road,8)-real(fg_road,8)-real(eg_road,8)) &
           + real(d_rdsndelta,8)*real(gardsn,8) - real(ggint_road,8)
d_faceterr  = d_roadstor - d_roadflux
pd%surferr = pd%surferr + d_faceterr
if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
  write(6,*) "aTEB road facet closure error:", maxval(abs(d_faceterr),mask=is_urban)
d_wallestor= sum(walle%storage-wallestorage_prev,dim=2)/real(ddt,8)
d_walleflux= real(sg_walle,8)+real(rg_walle,8)-real(fg_walle,8) - real(ggint_walle,8) - real(rgint_walle,8)
d_faceterr = d_wallestor - d_walleflux
pd%surferr = pd%surferr + d_faceterr
if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
  write(6,*) "aTEB walle facet closure error:", maxval(abs(d_faceterr),mask=is_urban)
d_wallwstor= sum(wallw%storage-wallwstorage_prev,dim=2)/real(ddt,8)
d_wallwflux= real(sg_wallw,8)+real(rg_wallw,8)-real(fg_wallw,8) - real(ggint_wallw,8) - real(rgint_wallw,8)
d_faceterr = d_wallwstor - d_wallwflux
pd%surferr = pd%surferr + d_faceterr
if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
  write(6,*) "aTEB wallw facet closure error:", maxval(abs(d_faceterr),mask=is_urban)
d_slabstor = sum(slab%storage-slabstorage_prev,dim=2)/real(ddt,8)
d_slabflux = -real(ggint_slab,8) - real(rgint_slab,8)
d_faceterr = d_slabstor - d_slabflux
pd%surferr = pd%surferr + d_faceterr
if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
  write(6,*) "aTEB slab facet closure error:", maxval(abs(d_faceterr),mask=is_urban)
if (intairtmeth/=0) then
  d_intmstor = sum(intm%storage-intmstorage_prev,dim=2)/real(ddt,8)
  d_intmflux = real(ggext_intm,8) - real(ggint_intm,8)
  d_faceterr = d_intmstor - d_intmflux
  pd%surferr = pd%surferr + d_faceterr
  if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
    write(6,*) "aTEB intm facet closure error:", maxval(abs(d_faceterr),mask=is_urban)
  d_roomstor = (room%storage(:,1)-roomstorage_prev(:,1))/real(ddt,8)
  d_roomflux = real(ggint_roof,8)+real(ggint_slab,8)-real(fp%intmassn,8)*d_intmflux                    &
            + (real(fp%bldheight,8)/real(fp%bldwidth,8))*(real(ggint_walle,8) + real(ggint_wallw,8))   &
            + real(int_infilflux,8) + real(d_ac_inside,8) + real(d_intgains_bld,8)
  d_faceterr = d_roomstor - d_roomflux
  pd%surferr = pd%surferr + d_faceterr
  if (any(is_urban .and. abs(d_faceterr)>=energytol)) &
    write(6,*) "aTEB room volume closure error:", maxval(abs(d_faceterr),mask=is_urban)
else
  d_intmstor = 0._8
  d_roomstor = 0._8
end if

pd%storage_flux = d_roofstor*real(fp%sigmabld,8)*(1._8-real(rfveg%sigma,8))           &
                + d_roadstor*(1._8-real(fp%sigmabld,8))*(1._8-real(cnveg%sigma,8))    &
                + d_wallestor*(1._8-real(fp%sigmabld,8))*real(fp%hwratio,8)           &
                + d_wallwstor*(1._8-real(fp%sigmabld,8))*real(fp%hwratio,8)           &
                + d_slabstor*real(fp%sigmabld,8)                                      &
                + d_intmstor*real(fp%sigmabld,8)*real(fp%intmassn,8)                  &
                + d_roomstor*real(fp%sigmabld,8)

! atmosphere energy flux = (SWdown-SWup) + (LWdown-LWup) - Turbulent + Anthropogenic
d_atmosflux = (real(a_sg,8)-real(a_sg,8)*real(u_alb,8)) + (real(a_rg,8)-real(sbconst,8)*(real(u_ts,8)+urbtemp)**4) &
            - (real(u_fg,8)+real(u_eg,8)+real(u_melt,8)) + real(pd%bldheat,8) + real(pd%bldcool,8)                 &
            + real(pd%traf,8) + real(fp%industryfg,8) + real(pd%intgains_full,8)
pd%atmoserr = pd%storage_flux - d_atmosflux

pd%d_atmosflux = d_atmosflux

if ( any(is_urban .and. abs(pd%atmoserr)>=energytol) ) then
  write(6,*) "UCLEM energy not conserved! Atmos. error:", maxval(abs(pd%atmoserr),mask=is_urban)
end if

return
end subroutine energyclosure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rebuild the state that energyclosure carries across a timestep but that the
! host does not hold in its restart dump, then release the first-call guard.
!
! energyclosure copies %storage into its *_prev arrays on entry and then
! recomputes %storage from the current nodetemp, so d_*stor is a difference
! across the step.  pd%delsoilmoist is formed the same way from the previous
! pd%soilmoist.  Neither %storage nor pd%soilmoist needs to be a restart field:
! both are pure functions of quantities a host restart does carry, %storage of
! the facet node temperatures and pd%soilmoist of the roof and road soil water.
! What a restart cannot supply is the value they held at the end of the previous
! run, and instead they hold whatever the host initialised them to.  pd%first_call
! exists to skip the two differences in that situation.
!
! Rebuilding both here from the restored prognostics removes the reason for the
! skip, so pd%first_call is cleared.  Clearing it is not optional: energyclosure
! returns before assigning pd%storage_flux, so without it the rebuilt values are
! never used and the first timestep of every run and every restart reports a zero
! storage flux to the host.
!
! Deliberately placed next to energyclosure.  The expressions below duplicate the
! %storage block in energyclosure and the pd%soilmoist update in uclem_eval, and
! must be changed with them.
!
! The caller is responsible for invoking this only after all geometry and
! material parameters are final, since bldheight, sigmabld, the vegetation
! fractions, the facet depths and the volumetric heat capacities all enter.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine uclem_rebuild_derived(roof,road,walle,wallw,slab,intm,room,          &
                                 fp,fp_roof,fp_road,fp_wall,fp_slab,fp_intm,    &
                                 rdhyd,rfhyd,cnveg,rfveg,pd,mp)

implicit none

integer, intent(in) :: mp
type(facetdata), intent(inout) :: roof, road, walle, wallw, slab, intm, room
type(facetparams), intent(in) :: fp_roof, fp_road, fp_wall, fp_slab, fp_intm
type(fparmdata), intent(in) :: fp
type(hydrodata), intent(in) :: rdhyd, rfhyd
type(vegdata), intent(in) :: cnveg, rfveg
type(pdiagdata), intent(inout) :: pd

real, dimension(mp) :: d_totdepth, veg_area

! With a prescribed internal air temperature the slab is not an active thermal
! facet.  uclem_eval pins slab%nodetemp to 1._8 on every prognostic call, and the
! slab facet closure requires the resulting storage change to be zero because
! both ggint_slab and rgint_slab are held at zero in that configuration.  Apply
! the same value here, otherwise the storage rebuilt from the start dump is
! overwritten on the first prognostic step and the difference is charged to the
! slab and atmosphere energy budgets.
if ( intairtmeth==0 ) then
  slab%nodetemp = 1._8
end if

! assume a_rho=1 for interior, as in energyclosure
room%storage(:,1) = real(fp%bldheight(:),8)*real(aircp,8)*room%nodetemp(:,1)

roof%storage(:,:)  = 0.5_8*real(fp_roof%depth(:,:),8)*real(fp_roof%volcp(:,:),8) &
                        *(roof%nodetemp(:,0:nl-1)+roof%nodetemp(:,1:nl))
road%storage(:,:)  = 0.5_8*real(fp_road%depth(:,:),8)*real(fp_road%volcp(:,:),8) &
                        *(road%nodetemp(:,0:nl-1)+road%nodetemp(:,1:nl))
walle%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8) &
                        *(walle%nodetemp(:,0:nl-1)+walle%nodetemp(:,1:nl))
wallw%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8) &
                        *(wallw%nodetemp(:,0:nl-1)+wallw%nodetemp(:,1:nl))
slab%storage(:,:)  = 0.5_8*real(fp_slab%depth(:,:),8)*real(fp_slab%volcp(:,:),8) &
                        *(slab%nodetemp(:,0:nl-1)+slab%nodetemp(:,1:nl))
intm%storage(:,:)  = 0.5_8*real(fp_intm%depth(:,:),8)*real(fp_intm%volcp(:,:),8) &
                        *(intm%nodetemp(:,0:nl-1)+intm%nodetemp(:,1:nl))

d_totdepth = sum(fp_road%depth,2)
veg_area   = max(cnveg%sigma*(1.-fp%sigmabld) + rfveg%sigma*fp%sigmabld, 1.e-10)
pd%soilmoist = ( rdhyd%soilwater*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld) &
               + rfhyd%soilwater*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld ) &
               / veg_area

pd%first_call = .false.

return
end subroutine uclem_rebuild_derived

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update water prognostic variables for roads and roofs

subroutine updatewater(ddt,surfwater,soilwater,leafwater,snow,den,alpha, &
                       snmelt,a_rnd,a_snd,eg_surf,eg_snow,d_tran,d_evap, &
                       d_totdepth,fp_vegrlai,iwbrelax,                   &
                       fp_sfc,fp_swilt,irrigwater,delintercept,          &
                       delsnow,mp)

implicit none

integer, intent(in) :: mp
integer, intent(in) :: iwbrelax
real, intent(in) :: ddt
real, dimension(mp), intent(inout) :: surfwater,soilwater,leafwater,snow,den,alpha
real, dimension(mp), intent(in) :: snmelt,a_rnd,a_snd,eg_surf,eg_snow
real, dimension(mp), intent(in) :: d_tran,d_evap,d_totdepth,fp_vegrlai
real, dimension(mp) :: modrnd, oldleafwater
real, dimension(mp), intent(in) :: fp_sfc, fp_swilt
real, dimension(mp), intent(out):: irrigwater, delintercept, delsnow

modrnd = max(a_rnd-d_evap/lv-max(maxvwatf*fp_vegrlai-leafwater,0.)/ddt,0.) ! rainfall reaching the soil under vegetation

! soil moisture based on Kowalczyk et al 1994 eq 12
! note that since sigmaf=1, then there is no soil evaporation, only transpiration.
! Evaporation only occurs from water on leafs.
surfwater = surfwater+ddt*(a_rnd-eg_surf/lv+snmelt)                                    ! surface
soilwater = soilwater+ddt*(modrnd+snmelt*den/waterden-d_tran/lv)/(waterden*d_totdepth) ! soil
oldleafwater = leafwater - ddt*d_evap/lv
leafwater = leafwater+ddt*(a_rnd-d_evap/lv)                                            ! leaf
leafwater = min(max(leafwater,0.),maxvwatf*fp_vegrlai)
delintercept = leafwater - oldleafwater

if (iwbrelax==1) then
  ! increase soil moisture for irrigation
  irrigwater = max(0.75*fp_swilt+0.25*fp_sfc-soilwater,0.)/(86400./ddt+1.)    ! 24h e-fold time
  soilwater=soilwater+irrigwater
else
  irrigwater = 0.
end if

! snow fields
snow  = snow + ddt*(a_snd-eg_snow/ls-snmelt)
delsnow = ddt*(a_snd-eg_snow/ls-snmelt)
den   = den + (maxsnowden-den)/(86400./(0.24*ddt)+1.)
alpha = alpha + (minsnowalpha-alpha)/(86400./(0.24*ddt)+1.)

return
end subroutine updatewater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsatr8(qsat,temp,ps)

implicit none

real(kind=8), dimension(:), intent(in) :: temp
real, dimension(:), intent(in) :: ps
real, dimension(:), intent(out) :: qsat
real, dimension(size(temp)) :: tempr4

tempr4 = real(temp)
call getqsat(qsat,tempr4,ps)

return
end subroutine getqsatr8

subroutine getqsat(qsat,temp,ps)

implicit none

real, dimension(:), intent(in) :: temp
real, dimension(:), intent(in) :: ps
real, dimension(:), intent(out) :: qsat
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix

real, dimension(0:220), parameter :: table =                                          &
  (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                               & !-146C
     6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                                            & !-141C
     36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                                         & !-136C
     0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,                 & !-131C
     0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,                 & !-126C
     0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,                  & !-121C
     0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,                      & !-116C
     0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,                        & !-111C
     0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,                           & !-106C
     0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,                            & !-101C
     0.001403, 0.001719, 0.002101, 0.002561, 0.003117, 0.003784,                      & !-95C
     0.004584, 0.005542, 0.006685, 0.008049, 0.009672, 0.01160, 0.01388, 0.01658,     & !-87C
     0.01977, 0.02353, 0.02796, 0.03316, 0.03925, 0.04638, 0.05472, 0.06444, 0.07577, & !-78C
     0.08894, 0.1042, 0.1220, 0.1425, 0.1662, 0.1936, 0.2252, 0.2615, 0.3032,         & !-69C
     0.3511, 0.4060, 0.4688, 0.5406, 0.6225, 0.7159, 0.8223, 0.9432, 1.080,           & !-60C
     1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476,                   & !-51C
     3.935, 4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,                          & !-43C
     10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,                   & !-34C
     27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85,                & !-24C
     77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67,                     & !-16C
     171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78,                 & !-8C
     353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78,                  & !0C
     656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,                  & !8C
     1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,                  & !16C
     1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,                  & !24C
     3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,                  & !32C
     5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,                  & !40C
     7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,                        & !47C
     11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,                   & !54C
     15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,                   & !61C
     21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,                   & !68C
     29845.0, 31169.0 /)


tdiff=min(max( temp+(urbtemp-123.16), 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/max(ps-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for calculating ustar and thetastar

subroutine getinvres(invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,a_umag,mode)

implicit none

integer, intent(in) :: mode
real, dimension(:), intent(in) :: ilzom
real, dimension(:), intent(in) :: zmin,sthetav,thetav
real, dimension(:), intent(in) :: a_umag
real, dimension(:), intent(out) :: invres,cd,z_on_l
real, dimension(:), intent(inout) :: olzoh
real, dimension(size(ilzom)) :: lna,thetavstar,integralh

lna=olzoh-ilzom
call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,zmin,ilzom,lna,mode)
invres=vkar*sqrt(cd)*a_umag/integralh
olzoh=lna+ilzom

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Translate forcing from one effective MOST height to another using
! the same Dyer-Hicks stability functions as UCLEM.

subroutine uclem_most_profile(a_temp_tgt,a_mixr_tgt,a_umag_tgt,cd_tgt,cdtq_tgt,z_on_l_tgt, &
                              a_temp_src,a_mixr_src,a_umag_src,s_temp,s_mixr,z_src,z_tgt,  &
                              z0m,z0h,mode,mp)

implicit none

integer, intent(in) :: mode, mp
real, dimension(mp), intent(in)  :: a_temp_src,a_mixr_src,a_umag_src
real, dimension(mp), intent(in)  :: s_temp,s_mixr,z_src,z_tgt,z0m,z0h
real, dimension(mp), intent(out) :: a_temp_tgt,a_mixr_tgt,a_umag_tgt
real, dimension(mp), intent(out) :: cd_tgt,cdtq_tgt,z_on_l_tgt

real, dimension(mp) :: z_src_eff,z_tgt_eff,z0m_eff,z0h_eff
real, dimension(mp) :: ilzom_src,ilzom_tgt,lna
real, dimension(mp) :: thetav_src,sthetav
real, dimension(mp) :: integralh_src,integralm_tgt,integralh_tgt
real, dimension(mp) :: cd_src,thetavstar,z_on_l_src
real, dimension(mp) :: ustar,tstar,qstar

real, parameter :: minz0m = 1.E-6
real, parameter :: minz0h = 1.E-8
real, parameter :: mineps = 1.E-10

z0m_eff   = max(z0m,minz0m)
z0h_eff   = max(z0h,minz0h)
z_src_eff = max(z_src,z0m_eff+1.E-4,z0h_eff+1.E-4)
z_tgt_eff = max(z_tgt,z0m_eff+1.E-4,z0h_eff+1.E-4)

ilzom_src = log(z_src_eff/z0m_eff)
lna       = log(z0m_eff/z0h_eff)

thetav_src = a_temp_src - urbtemp + a_temp_src*0.61*a_mixr_src
sthetav    = s_temp     - urbtemp + s_temp    *0.61*s_mixr

call dyerhicks(integralh_src,z_on_l_src,cd_src,thetavstar,thetav_src,sthetav, &
               a_umag_src,z_src_eff,ilzom_src,lna,mode)

ustar = sqrt(cd_src)*max(a_umag_src,0.01)
tstar = vkar*(a_temp_src-s_temp)/max(integralh_src,mineps)
qstar = vkar*(a_mixr_src-s_mixr)/max(integralh_src,mineps)

z_on_l_tgt = z_on_l_src*z_tgt_eff/z_src_eff
ilzom_tgt  = log(z_tgt_eff/z0m_eff)

call uclem_most_integrals(integralm_tgt,integralh_tgt,z_on_l_tgt,ilzom_tgt,lna)

a_umag_tgt = max(ustar*integralm_tgt/vkar,0.)
a_temp_tgt = s_temp + tstar*integralh_tgt/vkar
a_mixr_tgt = max(s_mixr + qstar*integralh_tgt/vkar,1.E-8)

cd_tgt   = (vkar/max(integralm_tgt,mineps))**2
cdtq_tgt = vkar*sqrt(cd_tgt)*max(a_umag_tgt,0.01)/max(integralh_tgt,mineps)

return
end subroutine uclem_most_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stability functions using Dyer-Hicks

subroutine dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,umagin,zmin,ilzom,lna,mode)

implicit none

integer, intent(in) :: mode
integer ic
real, dimension(:), intent(in) :: thetav
real, dimension(:), intent(in) :: sthetav,umagin,zmin,ilzom
real, dimension(:), intent(inout) :: lna
real, dimension(:), intent(out) :: cd,thetavstar
real, dimension(:), intent(out) :: integralh,z_on_l
real, dimension(size(thetav)) :: z0_on_l,zt_on_l,olzoh,umag
real, dimension(size(thetav)) :: pm0,ph0,pm1,ph1,integralm
!real, parameter :: aa1 = 3.8
!real, parameter :: bb1 = 0.5
!real, parameter :: cc1 = 0.3

umag = max(umagin, 0.01)
cd=(vkar/ilzom)**2                         ! first guess
call getlna(lna,cd,umag,zmin,ilzom,mode)
olzoh=ilzom+lna
integralh=sqrt(cd)*ilzom*olzoh/vkar        ! first guess
thetavstar=vkar*(thetav-sthetav)/integralh ! first guess

do ic=1,icmax
  z_on_l=vkar*zmin*grav*thetavstar/((thetav+urbtemp)*cd*umag**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-ilzom)
  zt_on_l  = z0_on_l*exp(-lna)
  where (z_on_l<0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph1     = (1.-16.*z_on_l)**(-0.5)
    integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
               +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = olzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5+b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l)+b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
    integralm = ilzom-(pm1-pm0)
    integralh = olzoh-(ph1-ph0)
  endwhere
  integralm = max( integralm, 1.e-10 )
  integralh = max( integralh, 1.e-10 )
  !where (z_on_l<=0.4)
    cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
  !elsewhere
  !  cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*(1.0+cc1* z_on_l**(1.-bb1)) &
  !      -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))/umag)**2
  !endwhere
  thetavstar= vkar*(thetav-sthetav)/integralh
  call getlna(lna,cd,umag,zmin,ilzom,mode)
end do

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Integrated Dyer-Hicks momentum and heat functions for a prescribed z/L.

subroutine uclem_most_integrals(integralm,integralh,z_on_l,ilzom,lna)

implicit none

real, dimension(:), intent(in)  :: z_on_l,ilzom,lna
real, dimension(:), intent(out) :: integralm,integralh
real, dimension(size(z_on_l)) :: z0_on_l,zt_on_l,olzoh
real, dimension(size(z_on_l)) :: pm0,ph0,pm1,ph1

z0_on_l = min(z_on_l*exp(-ilzom),10.)
zt_on_l = min(z0_on_l*exp(-lna),10.)
olzoh   = ilzom + lna

where (z_on_l<0.)
  pm0 = (1.-16.*z0_on_l)**(-0.25)
  ph0 = (1.-16.*zt_on_l)**(-0.5)
  pm1 = (1.-16.*z_on_l)**(-0.25)
  ph1 = (1.-16.*z_on_l)**(-0.5)
  integralm = ilzom - 2.*log((1.+1./pm1)/(1.+1./pm0))                         &
             - log((1.+1./pm1**2)/(1.+1./pm0**2))                              &
             + 2.*(atan(1./pm1)-atan(1./pm0))
  integralh = olzoh - 2.*log((1.+1./ph1)/(1.+1./ph0))
elsewhere
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
  pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
  ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5+b_1*(zt_on_l-(c_1/d_1))                &
        *exp(-d_1*zt_on_l)+b_1*c_1/d_1-1.)
  ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5+b_1*(z_on_l-(c_1/d_1))                  &
        *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
  integralm = ilzom - (pm1-pm0)
  integralh = olzoh - (ph1-ph0)
endwhere

integralm = max(integralm,1.E-10)
integralh = max(integralh,1.E-10)

return
end subroutine uclem_most_integrals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate roughness length for heat
!

pure subroutine getlna(lna,cd,umag,zmin,ilzom,mode)

implicit none

integer, intent(in) :: mode
integer iqu
! intent(inout), not intent(out): mode 4 leaves lna alone by design, so the value supplied
! by the caller has to survive the call.  With intent(out) that value is formally undefined
! on entry and a compiler is free to discard it under copy-in/copy-out.
real, dimension(:), intent(inout) :: lna
real, dimension(:), intent(in) :: cd,umag,zmin,ilzom
real re
real, parameter :: nu = 1.461E-5
!real, parameter :: eta0 = 1.827E-5
!real, parameter :: t0 = 291.15
!real, parameter :: c = 120.
!eta=eta0*((t0+c)/(theta+c))*(theta/t0)**(2./3.)
!nu=eta/rho

select case(mode) ! roughness length for heat
  case(1) ! zot=zom/10.
    lna=2.3
  case(2) ! Kanda et al 2007
    do iqu = 1,size(lna)
      re=max(sqrt(cd(iqu))*umag(iqu)*zmin(iqu)*exp(-ilzom(iqu))/nu,10.)
      !lna(iqu)=2.46*re**0.25-2. !(Brutsaet, 1982)
      lna(iqu)=1.29*re**0.25-2.  !(Kanda et al, 2007)
    end do
  case(3) ! zot=zom (neglect molecular diffusion)
    lna=0.
  case(4) ! user defined
    ! no change
end select

return
end subroutine getlna

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficients (modified to include 2nd wall)

subroutine getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,fp_hwratio,   &
                      fp_effhwratio,fp_vangle,fp_hangle,fp_fbeam,fp_vegsigmac,fp_roadalpha,fp_vegalphac,fp_wallalpha, &
                      ird_alpha,rdsndelta)

implicit none

integer k
real, dimension(:), intent(in) :: rdsndelta
real, dimension(:), intent(in) :: ird_alpha
real, dimension(:), intent(out) :: wallpsi,roadpsi
real, dimension(:), intent(in) :: fp_hwratio, fp_effhwratio
real, dimension(:), intent(in) :: fp_vangle,fp_hangle,fp_fbeam,fp_vegsigmac,fp_roadalpha,fp_vegalphac
real, dimension(:), intent(in) :: fp_wallalpha
real, dimension(:), intent(out) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(size(rdsndelta)) :: thetazero,walles,wallws,roads,ta,tc,xa,ya,roadnetalpha
real, dimension(size(rdsndelta)) :: nwalles,nwallws,nroads

wallpsi=0.5*(fp_effhwratio+1.-sqrt(fp_effhwratio*fp_effhwratio+1.))/fp_effhwratio
roadpsi=sqrt(fp_effhwratio*fp_effhwratio+1.)-fp_effhwratio

! integrate through 180 deg instead of 360 deg.  Hence partitioning to east and west facing walls
where (fp_vangle>=0.5*pi)
  walles=(1.-fp_fbeam)*wallpsi
  wallws=fp_fbeam/fp_hwratio + (1.-fp_fbeam)*wallpsi
  roads=(1.-fp_fbeam)*roadpsi
elsewhere
  ta=tan(fp_vangle)
  thetazero=asin(1./max(fp_hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  xa=min(max(fp_hangle-thetazero,0.),pi)-max(fp_hangle-pi+thetazero,0.)-min(fp_hangle+thetazero,0.)
  ya=cos(max(min(0.,fp_hangle),fp_hangle-pi))-cos(max(min(thetazero,fp_hangle),fp_hangle-pi)) &
    +cos(min(0.,-fp_hangle))-cos(min(thetazero,-fp_hangle)) &
    +cos(max(0.,pi-fp_hangle))-cos(max(thetazero,pi-fp_hangle))
  ! note that these terms now include the azimuth angle
  walles=fp_fbeam*(xa/fp_hwratio+ta*ya)/pi+(1.-fp_fbeam)*wallpsi
  wallws=fp_fbeam*((pi-2.*thetazero-xa)/fp_hwratio+ta*(tc-ya))/pi+(1.-fp_fbeam)*wallpsi
  roads=fp_fbeam*(2.*thetazero-fp_hwratio*ta*tc)/pi+(1.-fp_fbeam)*roadpsi
end where

! Calculate short wave reflections
roadnetalpha=rdsndelta*ird_alpha+(1.-rdsndelta)*((1.-fp_vegsigmac)*fp_roadalpha+fp_vegsigmac*fp_vegalphac)
sg_walle = walles
sg_wallw = wallws
sg_road = roads
do k = 1,nrefl
  nwalles = roadnetalpha*wallpsi*roads + fp_wallalpha*(1.-2.*wallpsi)*wallws
  nwallws = roadnetalpha*wallpsi*roads + fp_wallalpha*(1.-2.*wallpsi)*walles
  nroads = fp_wallalpha*(1.-roadpsi)*0.5*(walles+wallws)
  walles = nwalles
  wallws = nwallws
  roads = nroads
  sg_walle = sg_walle + walles
  sg_wallw = sg_wallw + wallws
  sg_road = sg_road + roads
end do
sg_roof=1.
sg_vegr=1.
sg_rfsn=1.
sg_rdsn=sg_road
sg_vegc=sg_road

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate longwave radiation coefficients (modified to include 2nd wall)

subroutine getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,fp_vegsigmac, &
                      fp_roademiss,fp_vegemissc,fp_wallemiss)

implicit none

integer k
real, dimension(:), intent(inout) :: d_netemiss
real, dimension(:), intent(inout) :: d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta
real, dimension(:), intent(in) :: fp_vegsigmac,fp_roademiss,fp_vegemissc,fp_wallemiss
real, dimension(:), intent(in) :: wallpsi,roadpsi
real, dimension(size(d_netemiss)) :: rcwa,rcra,rcwe,rcww,rcrw,rcrr,rcwr
real, dimension(size(d_netemiss)) :: ncwa,ncra,ncwe,ncww,ncrw,ncrr,ncwr

d_netemiss=d_rdsndelta*snowemiss+(1.-d_rdsndelta)*((1.-fp_vegsigmac)*fp_roademiss+fp_vegsigmac*fp_vegemissc)
d_cwa=wallpsi
d_cra=roadpsi
d_cw0=0.
d_cww=1.-2.*wallpsi
d_crw=0.5*(1.-roadpsi)
d_crr=0.
d_cwr=wallpsi
rcwa=d_cwa
rcra=d_cra
rcwe=d_cw0
rcww=d_cww
rcrw=d_crw
rcrr=d_crr
rcwr=d_cwr
do k = 1,nrefl
  ncwa=(1.-d_netemiss)*wallpsi*rcra+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-fp_wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-d_netemiss)*wallpsi*rcrw+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-d_netemiss)*wallpsi*rcrw+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcwe
  ncrw=(1.-fp_wallemiss)*(1.-roadpsi)*0.5*(rcww+rcwe)
  ncwr=(1.-d_netemiss)*wallpsi*rcrr+(1.-fp_wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-fp_wallemiss)*(1.-roadpsi)*rcwr
  rcwa=ncwa
  rcra=ncra
  rcwe=ncwe
  rcww=ncww
  rcrw=ncrw
  rcrr=ncrr
  rcwr=ncwr
  d_cwa=d_cwa+rcwa
  d_cra=d_cra+rcra
  d_cw0=d_cw0+rcwe
  d_cww=d_cww+rcww
  d_crw=d_crw+rcrw
  d_cwr=d_cwr+rcwr
  d_crr=d_crr+rcrr
end do

end subroutine getlwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes vegetation canopy temperature and canyon temperature)

subroutine solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,abase_road,                          &
                       sg_walle,rg_walle,fg_walle,acond_walle,abase_walle,                             &
                       sg_wallw,rg_wallw,fg_wallw,acond_wallw,abase_wallw,                             &
                       sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,abase_vegc,                          &
                       sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn, &
                       a_umag,a_rho,a_rg,a_rnd,a_snd,we,ww,                                            &
                       d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,                &
                       d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf,d_ac_inside,             &
                       d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,              &
                       d_cwr,d_totdepth,d_intgains_bld,d_roomstor,fgtop,egtop,evspsbltop,sbltop,       &
                       int_infilflux,int_infilfg,cyc_proportion,cyc_translation,ddt,                   &
                       cnveg,fp,fp_intm,fp_road,fp_roof,fp_wall,fp_slab,intm,pd,rdhyd,rfveg,road,      &
                       roof,room,slab,walle,wallw,iroomtemp,cvcoeff_roof,cvcoeff_walle,                &
                       cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,mp,diag)
implicit none

integer, intent(in) :: mp,diag
integer k,l
real, intent(in)    :: ddt
real, dimension(mp), intent(inout) :: rg_road,fg_road,eg_road,abase_road
real, dimension(mp), intent(inout) :: rg_walle,fg_walle,abase_walle
real, dimension(mp), intent(inout) :: rg_wallw,fg_wallw,abase_wallw
real, dimension(mp), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,abase_vegc
real, dimension(mp), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn
real, dimension(mp), intent(in) :: sg_road,sg_walle,sg_wallw,sg_vegc,sg_rdsn
real, dimension(mp), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd,d_sigd,we,ww
real, dimension(mp), intent(out) :: d_canyontemp,d_canyonmix,d_topu,d_netrad
real, dimension(mp), intent(out) :: d_ac_inside
real, dimension(mp), intent(inout) :: d_tempc,d_mixrc
real, dimension(mp), intent(inout) :: d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf
real, dimension(mp), intent(inout) :: d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr
real, dimension(mp), intent(inout) :: d_totdepth,d_intgains_bld
real, dimension(mp), intent(out) :: fgtop,egtop,evspsbltop,sbltop,int_infilflux
real, dimension(mp), intent(out) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(mp), intent(out) :: int_infilfg
real, dimension(mp), intent(out) :: iroomtemp
real, dimension(mp), intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(mp), intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(mp), intent(out) :: d_roomstor
real, dimension(mp) :: newval,sndepth,snlambda,ldratio,roadqsat,vegqsat,rdsnqsat
real, dimension(mp) :: cu,topinvres,dts,dtt,cduv,z_on_l,dumroaddelta,dumvegdelta,res
real, dimension(mp) :: aa,bb,cc,dd,ee,ff
real, dimension(mp) :: cyc_proportion,cyc_translation
real, dimension(mp,2) :: evct,evctx,oldval
real, dimension(mp) :: d_canyontemp_prev, fp_coeffbldheight
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_wall,fp_slab
type(hydrodata), intent(in) :: rdhyd
type(vegdata), intent(inout) :: cnveg, rfveg
type(facetdata), intent(in) :: intm
type(facetdata), intent(in) :: roof, slab
type(facetdata), intent(inout) :: road, room, walle, wallw
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd


! -------Initialise out variables------------------------
d_ac_inside = 0.

! snow conductance
sndepth  = rdhyd%snow*waterden/rdhyd%den
snlambda = icelambda*(rdhyd%den/waterden)**1.88

! first guess for canyon and room air temperature
! also guess for canyon veg, snow temperatures and water vapour mixing ratio
d_canyontemp      = d_tempc
d_canyontemp_prev = d_tempc
d_canyonmix       = d_mixrc
cnveg%temp        = d_tempc
rdsntemp          = real(road%nodetemp(:,1))
iroomtemp         = real(room%nodetemp(:,1))

rdsnmelt        = 0.
dumvegdelta     = 0. ! cray compiler bug


d_netrad=sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                                  &
        +(1.-d_rdsndelta)*(1.-cnveg%sigma)*fp_road%emiss*(real(road%nodetemp(:,0))+urbtemp)**4 &
        +(1.-d_rdsndelta)*cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4)

! Solve for canyon air temperature and water vapour mixing ratio
do l = 1,ncyits
  !  solve for aerodynamical resistance between canyon and atmosphere
  ! assume zoh=zom when coupling to canyon air temperature
  pd%lzoh = pd%lzom
  dts    = d_canyontemp + (d_canyontemp+urbtemp)*0.61*d_canyonmix
  dtt    = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
  call getinvres(topinvres,cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,3)
  call gettopu(d_topu,a_umag,z_on_l,fp%bldheight,cduv,pd%cndzmin,fp%hwratio,mp)

  if ( resmeth==0 ) then


    acond_road  = (11.8+4.2*sqrt((d_topu*abase_road)**2+cduv*a_umag**2))/(aircp*a_rho)  ! From Rowley, et al (1930)
    acond_walle = acond_road
    acond_wallw = acond_road
    acond_rdsn  = acond_road
    acond_vegc  = acond_road

else if ( resmeth==2 ) then

    cu = abase_road*d_topu
    where (cu<=5.)
      acond_road = (6.15+4.18*cu)/(aircp*a_rho)
    elsewhere
      acond_road = (7.51*cu**0.78)/(aircp*a_rho)
    end where
    acond_walle = acond_road
    acond_wallw = acond_road
    acond_rdsn  = acond_road
    acond_vegc  = acond_road

else

    acond_road  = d_topu*abase_road
    acond_walle = d_topu*abase_walle
    acond_wallw = d_topu*abase_wallw
    acond_vegc  = d_topu*abase_vegc
    acond_rdsn  = d_topu*abase_rdsn
  end if

  ! saturated mixing ratio for road
  call getqsatr8(roadqsat,road%nodetemp(:,0),d_sigd)   ! evaluate using pressure at displacement height

 ! correction for dew
  where (roadqsat<d_canyonmix)
    dumroaddelta=1.
  elsewhere
    dumroaddelta=d_roaddelta
  end where


  ! calculate canyon road latent heat flux
  aa=rdhyd%surfwater/ddt+a_rnd+rdsnmelt
  eg_road=lv*min(a_rho*d_roaddelta*(roadqsat-d_canyonmix)*acond_road,aa)

  ! Calculate longwave radiation emitted from the canyon floor
  ! MJT notes - This could be included within the iterative solver for snow and vegetation temperatures.
  ! However, it creates a (weak) coupling between these two variables and therefore could require
  ! a multivariate root finding method (e.g,. Broyden's method). Instead we explicitly solve for d_netrad,
  ! which allows us to decouple the solutions for snow and vegetation temperatures.
  d_netrad=sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                                   &
          +(1.-d_rdsndelta)*((1.-cnveg%sigma)*fp_road%emiss*(real(road%nodetemp(:,0))+urbtemp)**4 &
          +cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4))


  ! solve for road snow and canyon veg temperatures -------------------------------
  ldratio  = 0.5*( sndepth/snlambda + fp_road%depth(:,1)/fp_road%lambda(:,1) )
  cnveg%temp = cnveg%temp + 0.5 ! 1st guess
  rdsntemp   = rdsntemp + 0.5
  oldval(:,1) = cnveg%temp
  oldval(:,2) = rdsntemp
  call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                  sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                  a_rg,a_rho,a_rnd,a_snd,                                                       &
                  d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                     &
                  d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,ldratio,ddt,                         &
                  cnveg,fp,fp_wall,rdhyd,road,walle,wallw,mp)
  cnveg%temp = cnveg%temp - 0.5 ! 2nd guess
  rdsntemp   = rdsntemp - 0.5
  do k = 1,nfgits ! sectant
    evctx = evct
    call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                    sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                    a_rg,a_rho,a_rnd,a_snd,                                                       &
                    d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                     &
                    d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,ldratio,ddt,                         &
                    cnveg,fp,fp_wall,rdhyd,road,walle,wallw,mp)
    evctx = evct - evctx
    where ( abs(evctx(:,1))>tol )
      newval      = max(min(cnveg%temp-evct(:,1)*(cnveg%temp-oldval(:,1))/evctx(:,1),400.-urbtemp),200.-urbtemp)
      oldval(:,1) = cnveg%temp
      cnveg%temp  = newval
    end where
    where ( abs(evctx(:,2))>tol )
      newval      = max(min(rdsntemp-evct(:,2)*(rdsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
      oldval(:,2) = rdsntemp
      rdsntemp    = newval
    end where
  end do

  ! balance canyon latent heat budget
  aa = d_rdsndelta*acond_rdsn
  bb = (1.-d_rdsndelta)*(1.-cnveg%sigma)*dumroaddelta*acond_road
  cc = (1.-d_rdsndelta)*cnveg%sigma*(dumvegdelta*acond_vegc+(1.-dumvegdelta)/(1./max(acond_vegc,1.e-10)+res))
  dd = topinvres
  d_canyonmix = (aa*rdsnqsat+bb*roadqsat+cc*vegqsat+dd*d_mixrc)/(aa+bb+cc+dd)

  !!!!!!!!!!!!!!!!!!!! start interior models !!!!!!!!!!!!!!!!!!!!!!!!!!
  select case(intairtmeth)
    case(0) ! simple internal physics (per Thatcher and Hurley, 2012)
      call interiorflux_simple(mp,l,ddt,                                 &
                               d_intgains_bld,cyc_translation,                   &
                               d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor,  &
                               iroomtemp,int_infilflux,int_infilfg,              &
                               cvcoeff_roof,cvcoeff_walle,                       &
                               cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,         &
                               cvcoeff_intm2,fp,fp_roof,fp_wall,rfveg,roof,      &
                               walle,wallw,slab,intm)

    case(1,2) ! complex internal physics (per Lipson et al., 2018) or with interior frac
      call interiorflux_complex(mp,l,ddt,a_rho,we,ww,                    &
                                d_intgains_bld,cyc_proportion,d_canyontemp_prev, &
                                d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor, &
                                iroomtemp,int_infilflux,int_infilfg,             &
                                cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,        &
                                cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,fp,     &
                                fp_roof,fp_wall,fp_slab,fp_intm,roof,walle,      &
                                wallw,slab,intm,room,pd)
    case default
      write(6,*) "ERROR: Unknown option intairtmeth = ",intairtmeth
      stop
  end select

  d_canyontemp_prev = d_canyontemp

  ! update canyon temperature estimate
  aa = aircp*a_rho*topinvres
  bb = d_rdsndelta*aircp*a_rho*acond_rdsn
  cc = (1.-d_rdsndelta)*(1.-cnveg%sigma)*aircp*a_rho*acond_road
  dd = (1.-d_rdsndelta)*cnveg%sigma*aircp*a_rho*acond_vegc
  ee = fp%effhwratio*aircp*a_rho*acond_walle
  ff = fp%effhwratio*aircp*a_rho*acond_wallw
  d_canyontemp = (aa*d_tempc+bb*rdsntemp+cc*real(road%nodetemp(:,0))+dd*cnveg%temp+ee*real(walle%nodetemp(:,0)) &
                +ff*real(wallw%nodetemp(:,0))+d_traf+d_ac_canyon-int_infilfg)/(aa+bb+cc+dd+ee+ff)

end do


room%nodetemp(:,1) = real(iroomtemp,8)

fp_coeffbldheight = fp%effhwratio/fp%hwratio

! solve for canyon sensible heat flux
fg_walle = aircp*a_rho*(real(walle%nodetemp(:,0))-d_canyontemp)*acond_walle*fp_coeffbldheight

fg_wallw = aircp*a_rho*(real(wallw%nodetemp(:,0))-d_canyontemp)*acond_wallw*fp_coeffbldheight

fg_road  = aircp*a_rho*(real(road%nodetemp(:,0))-d_canyontemp)*acond_road

fg_vegc  = sg_vegc+rg_vegc-eg_vegc

fg_rdsn  = sg_rdsn+rg_rdsn-eg_rdsn-lf*rdsnmelt-gardsn*(1.-cnveg%sigma)

fgtop = fp%hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_canyon - int_infilfg


! solve for canyon latent heat flux
egtop = (1.-d_rdsndelta)*(1.-cnveg%sigma)*eg_road + (1.-d_rdsndelta)*cnveg%sigma*eg_vegc &
      + d_rdsndelta*eg_rdsn
evspsbltop = (1.-d_rdsndelta)*(1.-cnveg%sigma)*eg_road/lv + (1.-d_rdsndelta)*cnveg%sigma*eg_vegc/lv &
      + d_rdsndelta*eg_rdsn/ls
sbltop = d_rdsndelta*eg_rdsn/ls


! calculate longwave radiation
rg_walle = fp_wall%emiss*(a_rg*d_cwa+sbconst*(real(walle%nodetemp(:,0))+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)          &
         + sbconst*(real(wallw%nodetemp(:,0))+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
rg_wallw = fp_wall%emiss*(a_rg*d_cwa+sbconst*(real(wallw%nodetemp(:,0))+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)          &
         + sbconst*(real(walle%nodetemp(:,0))+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
rg_road = fp_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(real(road%nodetemp(:,0))+urbtemp)**4)                    &
        + sbconst*fp_wall%emiss*((real(walle%nodetemp(:,0))+urbtemp)**4+(real(wallw%nodetemp(:,0))+urbtemp)**4)*d_crw)
! outgoing longwave radiation
! note that eff terms are used for outgoing longwave radiation, whereas rg terms are used for heat conduction
d_canyonrgout = a_rg-d_rdsndelta*rg_rdsn-(1.-d_rdsndelta)*((1.-cnveg%sigma)*rg_road+cnveg%sigma*rg_vegc)              &
              - fp%hwratio*(rg_walle+rg_wallw)


return
end subroutine solvecanyon





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon/ internal heat fluxes and air temperatures

subroutine interiorflux_simple(mp,l,ddt,                                              &
                               d_intgains_bld,cyc_translation,                        &
                               d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor,       &
                               iroomtemp,int_infilflux,int_infilfg,                   &
                               cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,              &
                               cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,fp,fp_roof,   &
                               fp_wall,rfveg,roof,walle,wallw,slab,intm)

implicit none

!START------Initialise variables -------------------------------------------------
!--+-+------Canyon/internal air temperature loop (default ncyits=6) -----------+-+
!  | |                                                                         | |
!  | |        (intairtmeth=0) Evaluating with fixed internal air temperature   | |
!  | +---0.1: Estimate internal surface convection coefficients                | |
!  | +---0.2: Calculate internal conducted flux (ggint)                        | |
!  | +---0.3: Calculate cooling/heating fluxes (fixed air temp)                | |
!  | +---0.4: Calculate acflux into canyon ------------------------------------+ |
!---------------------------------------------------------------------------------

! passed variables
integer,                    intent(in) :: mp
real,                       intent(in) :: ddt
integer,                    intent(in) :: l
real, dimension(mp),     intent(in) :: d_intgains_bld
real, dimension(mp),     intent(in) :: cyc_translation
real, dimension(mp),     intent(in) :: d_canyontemp
real, dimension(mp),     intent(out) :: d_ac_canyon
real, dimension(mp),     intent(inout) :: d_ac_inside
real, dimension(mp),     intent(inout) :: iroomtemp
real, dimension(mp),     intent(out) :: int_infilflux
real, dimension(mp),     intent(out) :: int_infilfg
real, dimension(mp),     intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(mp),     intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(mp),     intent(out) :: d_roomstor

type(fparmdata),   intent(in) :: fp
type(facetparams), intent(in) :: fp_roof,fp_wall
type(facetdata),   intent(in) :: roof,walle,wallw,slab,intm
type(vegdata),     intent(in) :: rfveg

! local variables
real, dimension(mp) :: ggint_roof,ggint_walle,ggint_wallw,ggint_slab
real, dimension(mp) :: ac_coeff


iroomtemp = fp%bldairtemp + cyc_translation
d_roomstor = 0.

! ---0.1: Estimate internal surface convection coefficients -------------
call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,  &
                    cvcoeff_intm1,cvcoeff_intm2,iroomtemp,roof,slab,intm,mp)
! ---0.2: Calculate internal conducted flux (ggint) ----------------------
ggint_slab=0.  ! no slab for intairtmeth=0
ggint_roof  = cvcoeff_roof*(real(roof%nodetemp(:,nl))-iroomtemp)
ggint_walle = cvcoeff_walle*(real(walle%nodetemp(:,nl))-iroomtemp)
ggint_wallw = cvcoeff_wallw*(real(wallw%nodetemp(:,nl))-iroomtemp)




! ---0.3: Calculate cooling/heating fluxes (fixed air temp) --------------
! flux into room potentially pumped out into canyon (depends on AC method)
d_ac_inside = -ggint_roof - ggint_slab                         &
              - (ggint_walle+ggint_wallw)*(fp%bldheight/fp%bldwidth)

! ---0.4: Calculate acflux into canyon -----------------------------------
! update heat pumped into canyon
where ( d_ac_inside<0. )
  ! original UCLEM COP calculation (Lipson et al., 2018)
  ! ac_coeff = max(acfactor*(d_canyontemp-iroomtemp)/(iroomtemp+urbtemp), 1./ac_copmax) ! T&H2012 Eq. 10
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
  ac_coeff = max( acfactor*(d_canyontemp+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
elsewhere
  ac_coeff = 0.
end where
select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, COP of 1.0)
  case(0) ! unrealistic cooling (buildings act as heat sink)
    d_ac_canyon  = 0.
  case(1) ! d_ac_canyon pumps conducted heat + ac waste heat back into canyon
    d_ac_canyon = max(0.,-d_ac_inside*(1.+ac_coeff)*fp%sigmabld/(1.-fp%sigmabld))  ! canyon domain W/m/m
  case(2) ! reversible heating and cooling (for testing energy conservation)
    d_ac_canyon  = -d_ac_inside*fp%sigmabld/(1.-fp%sigmabld)
  case DEFAULT
    write(6,*) "ERROR: Unknown acmeth mode ",acmeth
    stop
end select

! immediately eject intgains as infiltration flux
int_infilflux = -d_intgains_bld
int_infilfg = int_infilflux*fp%sigmabld/(1.-fp%sigmabld)

return
end subroutine interiorflux_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon/ internal heat fluxes and air temperatures
subroutine interiorflux_complex(mp,l,ddt,a_rho,d_weu,d_wwu,                            &
                                d_intgains_bld,cyc_proportion,d_canyontemp_prev,       &
                                d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor,       &
                                iroomtemp,int_infilflux,int_infilfg,                   &
                                cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, &
                                cvcoeff_intm1,cvcoeff_intm2,fp,fp_roof,fp_wall,        &
                                fp_slab,fp_intm,roof,walle,wallw,slab,intm,room,pd)

implicit none

!START------Initialise variables -------------------------------------------------
!--+-+------Canyon/internal air temperature loop (default ncyits=6) -----------+-+
!  |          (intairtmeth=1) Evaluating with varying internal air temperature   |
!  +-----1.1: Estimate internal surface convection coefficients                  |
!  +-----1.2: Calculate air exchange from infiltration and open windows          |
!  +-----1.3: Calculate cooling/heating fluxes (varying air temp)                |
!  +-----1.4: Calculate internal air temperature implicitly                      |
!  +--+         Final loop                                                       |
!  |  +-----1.5: Final iteration flux and temperature                            |
!  +-----1.8: Calculate acflux into canyon and update canyon temperature --------+
!---------------------------------------------------------------------------------

! passed variables
integer,                    intent(in) :: mp                   ! tile vector length
integer,                    intent(in) :: l                       ! ncyits loop integer
real,                       intent(in) :: ddt                     ! timestep
real, dimension(mp),     intent(in) :: a_rho                   ! atmospheric air density
real, dimension(mp),     intent(in) :: d_weu,d_wwu             ! wind speed wall east and west
real, dimension(mp),     intent(in) :: d_intgains_bld          ! internal gain flux for building
real, dimension(mp),     intent(in) :: cyc_proportion          ! cyclic proportion of conditioning units
real, dimension(mp),     intent(in) :: d_canyontemp            ! canyon air temperature (diagnostic)
real, dimension(mp),     intent(in) :: d_canyontemp_prev       ! previous canyon air temperature
real, dimension(mp),     intent(out) :: d_ac_canyon            ! conditioning flux into canyon
real, dimension(mp),     intent(inout) :: d_ac_inside          ! conditioning flux inside
real, dimension(mp),     intent(inout) :: iroomtemp            ! update of room air temperature
real, dimension(mp),     intent(out) :: int_infilflux          ! internal infiltration flux
real, dimension(mp),     intent(out) :: int_infilfg            ! external infiltration flux
real, dimension(mp),     intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab
real, dimension(mp),     intent(out) :: cvcoeff_intm1,cvcoeff_intm2
real, dimension(mp),     intent(out) :: d_roomstor

type(fparmdata),   intent(in) :: fp
type(facetparams), intent(in) :: fp_roof,fp_wall,fp_slab,fp_intm
type(facetdata),   intent(in) :: roof,walle,wallw,slab,intm
type(facetdata),   intent(inout) :: room
type(pdiagdata),   intent(in) :: pd

! local variables
integer iqu
real, dimension(mp) :: ggint_roof, ggint_walle, ggint_wallw ! facet convection fluxes
real, dimension(mp) :: ggint_slab, ggext_intm, ggint_intm   ! facet convection fluxes
real, dimension(mp) :: rm,rf,we,ww,sl,im1,im2,infl
real, dimension(mp) :: ac_coeff,d_ac_cool,d_ac_heat,ac_load
real, dimension(mp) :: infl_dynamic,itemp
real, dimension(mp) :: ac_heat_temp,ac_cool_temp
real, dimension(mp) :: int_airden
real, dimension(mp) :: d_openwindows
real, dimension(mp) :: dummy, canyontemp_ave, tempr4

! calculate canyon temperature for fluxes
canyontemp_ave = infilalpha*d_canyontemp + (1.-infilalpha)*d_canyontemp_prev


! ---1.1: Estimate internal surface convection coefficients -------------
tempr4 = real(room%nodetemp(:,1))
call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,       &
                    cvcoeff_intm1,cvcoeff_intm2,tempr4,roof,slab,intm,mp)


ac_heat_temp = fp%bldairtemp - ac_deltat
ac_cool_temp = fp%bldairtemp + ac_deltat


! ---1.2: Calculate air exchange from infiltration and open windows -----
call calc_openwindows(d_openwindows,fp,iroomtemp,canyontemp_ave,roof,walle,wallw,slab,mp)
select case(infilmeth)
  case(0) ! constant
    infl_dynamic = fp%infilach
  case(1) ! EnergyPlus with BLAST coefficients (Coblenz and Achenbach, 1963)
    infl_dynamic = fp%infilach*(0.606 + 0.03636*abs(iroomtemp-canyontemp_ave) + 0.1177*0.5*(d_wwu+d_weu) + 0.)
  !case(2) ! AccuRate (Chen, 2010)
  !  not yet implemented
  case default
    write(6,*) "ERROR: Unknown option infilmeth = ",infilmeth
    stop
end select


! ---1.3: Calculate cooling/heating fluxes (varying air temp)------------
! int_airden = air_pr/(air_rd*(iroomtemp+urbtemp))
int_airden = 1.
infl = aircp*a_rho*fp%bldheight*(infl_dynamic+d_openwindows*fp%ventilach)/3600.
rm = aircp*int_airden*fp%bldheight/ddt
rf = cvcoeff_roof
we = (fp%bldheight/fp%bldwidth)*cvcoeff_walle
ww = (fp%bldheight/fp%bldwidth)*cvcoeff_wallw
sl = cvcoeff_slab
im1 = cvcoeff_intm1*real(fp%intmassn)
im2 = cvcoeff_intm2*real(fp%intmassn)


! ---1.4: Calculate cooling/heating fluxes (varying air temp)------------

! remove acflux component added to implicit formulation of air temp
itemp = (rm*real(room%nodetemp(:,1))              & ! room temperature
       + rf*real(roof%nodetemp(:,nl))             & ! roof conduction
       + we*real(walle%nodetemp(:,nl))            & ! wall conduction east
       + ww*real(wallw%nodetemp(:,nl))            & ! wall conduction west
       + sl*real(slab%nodetemp(:,nl))             & ! slab conduction
       + im1*real(intm%nodetemp(:,0))             & ! mass conduction side 1
       + im2*real(intm%nodetemp(:,nl))            & ! mass conduction side 2
       + infl*canyontemp_ave                      & ! weighted infiltration
       + d_intgains_bld                           & ! internal gains (explicit)
       )/(rm+rf+we+ww+sl+im1+im2+infl)

d_ac_cool = 0.
d_ac_heat = 0.

! heating load
where ( itemp < ac_heat_temp )
  !ac_load = rm*(ac_heat_temp-itemp)
  ac_load = (rm+rf+we+ww+sl+im1+im2+infl)*(ac_heat_temp-itemp)
  d_ac_heat = min(ac_heatcap*fp%bldheight,ac_load)*pd%ac_heat_on
end where
! cooling load
where ( itemp > ac_cool_temp )
  !ac_load = rm*(itemp-ac_cool_temp)
  ac_load = (rm+rf+we+ww+sl+im1+im2+infl)*(itemp-ac_cool_temp)
  d_ac_cool = min(ac_coolcap*fp%bldheight,ac_load)*pd%ac_cool_on
end where

! modulate acfluxes to account for proportional heating/cooling area and time of day
if (intairtmeth==1) then
  d_ac_heat = d_ac_heat*fp%heatprop*cyc_proportion
  d_ac_cool = d_ac_cool*fp%coolprop*cyc_proportion
end if

d_ac_inside = d_ac_heat - d_ac_cool

! ---1.4: Calculate internal air temperature implicitly -----------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fully implicit - old temps with acflux included

iroomtemp = (rm*real(room%nodetemp(:,1))              & ! room temperature
           + rf*real(roof%nodetemp(:,nl))             & ! roof conduction
           + we*real(walle%nodetemp(:,nl))            & ! wall conduction east
           + ww*real(wallw%nodetemp(:,nl))            & ! wall conduction west
           + sl*real(slab%nodetemp(:,nl))             & ! slab conduction
           + im1*real(intm%nodetemp(:,0))             & ! mass conduction side 1
           + im2*real(intm%nodetemp(:,nl))            & ! mass conduction side 2
           + infl*canyontemp_ave                      & ! weighted infiltration
           + d_ac_inside                              & ! ac flux (iterative)
           + d_intgains_bld                           & ! internal gains (explicit)
           )/(rm+rf+we+ww+sl+im1+im2+infl)
!iroomtemp = itemp + d_ac_inside/rm ! time split version


! calculate infiltration flux with dampening factor from previous timestep (reduce oscillation)
int_infilflux = infl*(canyontemp_ave - iroomtemp)
int_infilfg = (fp%sigmabld/(1.-fp%sigmabld))*int_infilflux

d_roomstor = (iroomtemp - real(room%nodetemp(:,1)))*(aircp*int_airden*fp%bldheight)/ddt

! ---1.5 (Depreciated by MJT)

! ---1.6: Calculate acflux into canyon and update canyon temperature-----

where ( d_ac_cool>0. )
  ! original UCLEM COP calculation (Lipson et al., 2018)
  ! ac_coeff = max(acfactor*(d_canyontemp-iroomtemp)/(iroomtemp+urbtemp), 1./ac_copmax) ! T&H2012 Eq. 10
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
  !ac_coeff = max( acfactor*(canyontemp_ave+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
  ac_coeff = max( acfactor*(d_canyontemp+30.-iroomtemp)/(iroomtemp+30.+urbtemp), 0. )
elsewhere
  ac_coeff = 0.
end where
select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible)
  case(0) ! unrealistic cooling (buildings act as heat sink)
    d_ac_canyon = 0.
  case(1) ! d_ac_canyon pumps conducted heat + ac waste heat back into canyon
    d_ac_canyon = max(0.,d_ac_cool*(1. + ac_coeff)*fp%sigmabld/(1.-fp%sigmabld))  ! canyon domain W/m2
  case(2) ! reversible heating and cooling (for testing energy conservation)
    d_ac_canyon = -d_ac_inside*fp%sigmabld/(1.-fp%sigmabld)
  case DEFAULT
    write(6,*) "ERROR: Unknown acmeth mode ",acmeth
    stop
end select

! ! COP printout when canyon air is >35°C
! if ( l==ncyits ) then
!   if (( (d_canyontemp(1) + urbtemp) > (35. + 273.15) ) .and. ( (d_canyontemp(1) + urbtemp) < (36. + 273.15) )) then
!     if ( ac_coeff(1) > 0 ) then
!       write(6,*) 'COP:', 1./ac_coeff
!     end if
!   end if
! end if


return
end subroutine interiorflux_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon veg and snow fluxes

subroutine canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,       &
                      sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat,  &
                      a_rg,a_rho,a_rnd,a_snd,                                                        &
                      d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                      &
                      d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,ldratio,ddt,                          &
                      cnveg,fp,fp_wall,rdhyd,road,walle,wallw,mp)


implicit none

integer, intent(in) :: mp
real, intent(in) :: ddt
real, dimension(mp,2), intent(out) :: evct
real, dimension(mp), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta
real, dimension(mp), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat
real, dimension(mp), intent(in) :: sg_vegc,sg_rdsn
real, dimension(mp), intent(in) :: a_rg,a_rho,a_rnd,a_snd,d_sigd
real, dimension(mp), intent(in) :: ldratio
real, dimension(mp), intent(inout) :: d_canyontemp,d_canyonmix,d_netrad,d_tranc,d_evapc
real, dimension(mp), intent(inout) :: d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac
real, dimension(mp) :: ff,f1,f2,f3,f4
real, dimension(mp) :: snevap
type(vegdata), intent(in) :: cnveg
type(facetparams), intent(in) :: fp_wall
type(hydrodata), intent(in) :: rdhyd
type(facetdata), intent(in) :: road, walle, wallw
type(fparmdata), intent(in) :: fp

! estimate mixing ratio for vegetation and snow
call getqsat(vegqsat,cnveg%temp,d_sigd)
call getqsat(rdsnqsat,rdsntemp,d_sigd)

! correction for dew
where (vegqsat<d_canyonmix)
  dumvegdelta=1.
elsewhere
  dumvegdelta=d_vegdeltac
end where

! vegetation transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (cnveg%zo<0.5)
  ff=1.1*sg_vegc/max(cnveg%lai*150.,1.E-8)
elsewhere
  ff=1.1*sg_vegc/max(cnveg%lai*30.,1.E-8)
end where
f1=(1.+ff)/(ff+cnveg%rsmin*cnveg%lai/5000.)
f2=max(0.5*(fp%sfc-fp%swilt)/max(rdhyd%soilwater-fp%swilt,1.E-9),1.)
f3=max(1.-0.00025*(vegqsat-d_canyonmix)*d_sigd/0.622,0.5) ! increased limit from 0.05 to 0.5 following Mk3.6
f4=max(1.-0.0016*(298.-urbtemp-d_canyontemp)**2,0.05)     ! 0.2 in Mk3.6
res=max(30.,cnveg%rsmin*f1*f2/(f3*f4))

! solve for vegetation and snow sensible heat fluxes
fg_vegc=aircp*a_rho*(cnveg%temp-d_canyontemp)*acond_vegc
fg_rdsn=aircp*a_rho*(rdsntemp-d_canyontemp)*acond_rdsn

! calculate longwave radiation for vegetation and snow
rg_vegc=cnveg%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(cnveg%temp+urbtemp)**4)    &
                  +sbconst*fp_wall%emiss*((real(walle%nodetemp(:,0))+urbtemp)**4          &
                  +(real(wallw%nodetemp(:,0))+urbtemp)**4)*d_crw)
rg_rdsn=snowemiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(rdsntemp+urbtemp)**4)        &
                  +sbconst*fp_wall%emiss*((real(walle%nodetemp(:,0))+urbtemp)**4          &
                  +(real(wallw%nodetemp(:,0))+urbtemp)**4)*d_crw)

! estimate snow melt
rdsnmelt=min(max(0.,rdsntemp+(urbtemp-273.16))*icecp*rdhyd%snow/(ddt*lf),rdhyd%snow/ddt)

! calculate latent heat of transpiration and evaporation of in-canyon vegetation
d_tranc=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_canyonmix)/(1./max(acond_vegc,1.e-10)+res),0.), &
               max((rdhyd%soilwater-fp%swilt)*d_totdepth*waterden/ddt,0.))
d_evapc=lv*min(dumvegdelta*a_rho*(vegqsat-d_canyonmix)*acond_vegc,rdhyd%leafwater/ddt+a_rnd)
eg_vegc=d_evapc+d_tranc

! calculate canyon snow latent heat and ground fluxes
snevap=min(a_rho*max(0.,rdsnqsat-d_canyonmix)*acond_rdsn,rdhyd%snow/ddt+a_snd-rdsnmelt)
eg_rdsn=lv*snevap
rdsnmelt=rdsnmelt+snevap
gardsn=(rdsntemp-real(road%nodetemp(:,0)))/ldratio ! use road temperature to represent canyon bottom surface temperature
                                                   ! (i.e., we have omitted soil under vegetation temperature)

! vegetation energy budget error term
evct(:,1) = sg_vegc+rg_vegc-fg_vegc-eg_vegc

! road snow energy balance error term
evct(:,2) = sg_rdsn+rg_rdsn-fg_rdsn-eg_rdsn-lf*rdsnmelt-gardsn*(1.-cnveg%sigma)

return
end subroutine canyonflux




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for roof fluxes

subroutine solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
                     sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
                     sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
                     a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,         &
                     d_sigr,ddt,fp_roof,rfhyd,rfveg,roof,fp,mp)

implicit none

integer, intent(in) :: mp
integer k
real, intent(in) :: ddt
real, dimension(mp), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn
real, dimension(mp), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(mp), intent(inout) :: rg_roof,eg_roof,acond_roof
real, dimension(mp), intent(in) :: sg_rfsn,sg_vegr,sg_roof
real, dimension(mp), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(mp), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr
real, dimension(mp), intent(inout) :: d_rfsndelta,d_vegdeltar,d_roofdelta
real, dimension(mp) :: lzomroof,lzohroof,qsatr,dts,dtt,cdroof,z_on_l,newval,ldratio
real, dimension(mp) :: aa,dd,ee
real, dimension(mp,2) :: oldval,evctx,evctveg
type(facetparams), intent(in) :: fp_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: roof
type(fparmdata), intent(in) :: fp

lzomroof=log(d_rfdzmin/zoroof)
lzohroof=2.3+lzomroof
call getqsatr8(qsatr,roof%nodetemp(:,0),d_sigr)
dts=real(roof%nodetemp(:,0)) + (real(roof%nodetemp(:,0))+urbtemp)*0.61*d_roofdelta*qsatr
dtt=d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_roof,cdroof,z_on_l,lzohroof,lzomroof,d_rfdzmin,dts,dtt,a_umag,1)

! update green roof and snow temperature
rfveg%temp=d_tempr
rfsntemp  =max(min(real(roof%nodetemp(:,0)),300.-urbtemp),100.-urbtemp)
rg_vegr = fp_roof%emiss*(a_rg-sbconst*(real(roof%nodetemp(:,0))+urbtemp)**4) ! 1st guess
rg_rfsn = fp_roof%emiss*(a_rg-sbconst*(real(roof%nodetemp(:,0))+urbtemp)**4) ! 1st guess
eg_vegr = 0.
eg_rfsn = 0.
rfsnmelt = 0.
garfsn = 0.
d_tranr = 0.
d_evapr = 0.
acond_vegr = acond_roof
acond_rfsn = acond_roof
if ( any( d_rfsndelta>0. .or. rfveg%sigma>0. ) ) then
  evctveg = 0.
  oldval(:,1)=rfveg%temp+0.5
  oldval(:,2)=rfsntemp+0.5
  call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,                &
                d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,mp)
  ! turn off roof snow and roof vegetation if they are not needed
  where ( rfveg%sigma>0. )
    rfveg%temp=rfveg%temp-0.5
  end where
  where ( d_rfsndelta>0. )
    rfsntemp  =rfsntemp-0.5
  end where
  do k=1,nfgits
    evctx=evctveg
    call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                  sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                  d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,                &
                  d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,mp)
    evctx=evctveg-evctx
    where ( abs(evctx(:,1))>tol .and. rfveg%sigma>0. )
      newval=rfveg%temp-evctveg(:,1)*(rfveg%temp-oldval(:,1))/evctx(:,1)
      oldval(:,1)=rfveg%temp
      rfveg%temp=newval
    end where
    where ( abs(evctx(:,2))>tol .and. d_rfsndelta>0. )
      newval=max(min(rfsntemp-evctveg(:,2)*(rfsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
      oldval(:,2)=rfsntemp
      rfsntemp=newval
    end where
  end do
end if
fg_vegr=sg_vegr+rg_vegr-eg_vegr
fg_rfsn=sg_rfsn+rg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)

! estimate roof latent heat flux
where ( qsatr<d_mixrr )
  ! dew
  eg_roof=lv*a_rho*(qsatr-d_mixrr)*acond_roof
elsewhere
  ! evaporation
  aa=rfhyd%surfwater/ddt+a_rnd+rfsnmelt
  eg_roof=lv*min(a_rho*d_roofdelta*(qsatr-d_mixrr)*acond_roof,aa)
end where

! calculate net roof longwave radiation
! (sensible heat flux will be updated in solvetridiag)
rg_roof=fp_roof%emiss*(a_rg-sbconst*(real(roof%nodetemp(:,0))+urbtemp)**4)

return
end subroutine solveroof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof and snow fluxes

subroutine roofflux(evct,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,   &
                    sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,   &
                    d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,               &
                    d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,mp)

implicit none

integer, intent(in) :: mp
real, intent(in) :: ddt
real, dimension(mp,2), intent(inout) :: evct
real, dimension(mp), intent(in) :: rfsntemp,sg_vegr,sg_rfsn
real, dimension(mp), intent(inout) :: rfsnmelt,garfsn
real, dimension(mp), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(mp), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn
real, dimension(mp), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(mp), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr
real, dimension(mp), intent(inout) :: d_vegdeltar,d_rfsndelta
real, dimension(mp) :: lzomvegr,lzohvegr,vwetfac,dts,dtt,z_on_l,ff,f1,f2,f3,f4,cdvegr
real, dimension(mp) :: vegqsat,dumvegdelta,res,sndepth,snlambda,ldratio,lzosnow,rfsnqsat,cdrfsn
real, dimension(mp) :: lzotdum, snevap
type(facetparams), intent(in) :: fp_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(in) :: rfveg
type(facetdata), intent(in) :: roof
type(fparmdata), intent(in) :: fp

call getqsat(vegqsat,rfveg%temp,d_sigr)
where ( vegqsat<d_mixrr )
  dumvegdelta = 1.
elsewhere
  dumvegdelta = d_vegdeltar
end where

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where ( rfveg%zo<0.5 )
  ff = 1.1*sg_vegr/max(rfveg%lai*150.,1.E-8)
elsewhere
  ff = 1.1*sg_vegr/max(rfveg%lai*30.,1.E-8)
end where
f1 = (1.+ff)/(ff+rfveg%rsmin*rfveg%lai/5000.)
f2 = max(0.5*(fp%sfc-fp%swilt)/max(rfhyd%soilwater-fp%swilt,1.E-9),1.)
f3 = max(1.-.00025*(vegqsat-d_mixrr)*d_sigr/0.622,0.5)
f4 = max(1.-0.0016*((298.-urbtemp)-d_tempr)**2,0.05)
res = max(30.,rfveg%rsmin*f1*f2/(f3*f4))

vwetfac = max(min((rfhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
vwetfac = (1.-dumvegdelta)*vwetfac+dumvegdelta
lzomvegr = log(d_rfdzmin/rfveg%zo)
! xe is a dummy variable for lzohvegr
lzohvegr = 2.3+lzomvegr
dts = rfveg%temp + (rfveg%temp+urbtemp)*0.61*vegqsat*vwetfac
dtt = d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_vegr,cdvegr,z_on_l,lzohvegr,lzomvegr,d_rfdzmin,dts,dtt,a_umag,1)
! acond_vegr is multiplied by a_umag

where ( rfveg%sigma>0. )
  ! longwave radiation
  rg_vegr=rfveg%emiss*(a_rg-sbconst*(rfveg%temp+urbtemp)**4)

  ! sensible heat flux
  fg_vegr=aircp*a_rho*(rfveg%temp-d_tempr)*acond_vegr

  ! calculate transpiration and evaporation of roof vegetation
  d_tranr=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_mixrr)/(1./acond_vegr+res),0.), &
                 max((rfhyd%soilwater-fp%swilt)*fp%rfvegdepth*waterden/ddt,0.))
  d_evapr=lv*min(dumvegdelta*a_rho*(vegqsat-d_mixrr)*acond_vegr,rfhyd%leafwater/ddt+a_rnd)
  eg_vegr=d_evapr+d_tranr

  ! balance green roof energy budget
  evct(:,1)=sg_vegr+rg_vegr-fg_vegr-eg_vegr
end where


! snow conductance
sndepth=rfhyd%snow*waterden/rfhyd%den
snlambda=icelambda*(rfhyd%den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+fp_roof%depth(:,1)/fp_roof%lambda(:,1))

! Update roof snow energy budget
lzosnow=log(d_rfdzmin/zosnow)
call getqsat(rfsnqsat,rfsntemp,d_sigr)
lzotdum=2.3+lzosnow
dts=rfsntemp + (rfsntemp+urbtemp)*0.61*rfsnqsat
call getinvres(acond_rfsn,cdrfsn,z_on_l,lzotdum,lzosnow,d_rfdzmin,dts,dtt,a_umag,1)
! acond_rfsn is multiplied by a_umag

where ( d_rfsndelta>0. )
  rfsnmelt=min(max(0.,rfsntemp+(urbtemp-273.16))*icecp*rfhyd%snow/(ddt*lf),rfhyd%snow/ddt)
  rg_rfsn=snowemiss*(a_rg-sbconst*(rfsntemp+urbtemp)**4)
  fg_rfsn=aircp*a_rho*(rfsntemp-d_tempr)*acond_rfsn
  snevap=min(a_rho*max(0.,rfsnqsat-d_mixrr)*acond_rfsn,rfhyd%snow/ddt+a_snd-rfsnmelt)
  eg_rfsn=lv*snevap
  rfsnmelt=rfsnmelt+snevap
  garfsn=(rfsntemp-real(roof%nodetemp(:,0)))/ldratio

  ! balance snow energy budget
  evct(:,2)=sg_rfsn+rg_rfsn-fg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)
end where

return
end subroutine roofflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define weights during the diurnal cycle. Array starts at 1am.

subroutine getdiurnal(fp_ctime,fp_weekdayload,icyc_traffic,icyc_basedemand,icyc_proportion,icyc_translation)

implicit none

real, dimension(:), intent(in) :: fp_ctime,fp_weekdayload
real, dimension(:), intent(out) :: icyc_traffic,icyc_basedemand,icyc_proportion,icyc_translation
real, dimension(size(fp_ctime)) :: real_p
integer, dimension(size(fp_ctime)) :: int_p

! ! base temperature translation cycle optimised for Melbourne
! real, dimension(25), save :: trancycle =(/ -1.2, -1.7, -2.4, -2.9, -3.0, -2.1,  &
!                                            -1.1, -0.6, -0.7, -0.8, -1.2, -1.4,  &
!                                            -1.6, -1.8, -1.6, -1.0, -0.2,  1.1,  &
!                                             2.1,  2.7,  3.0,  2.6,  1.6,  0.2, -1.2 /)

int_p=int(24.*fp_ctime)
real_p=24.*fp_ctime-real(int_p)
where (int_p<1) int_p=int_p+24

if (any(int_p < 1 .or. int_p > 24)) then

icyc_traffic = 0.0
icyc_basedemand = 0.0
icyc_proportion = 0.0
icyc_translation = 0.0

write(6,*) "ERROR in getdiurnal: int_p out of range [1,24]"
write(6,*) "  min(int_p)=", minval(int_p), "  max(int_p)=", maxval(int_p)
write(6,*) "  min(fp_ctime)=", minval(fp_ctime), "  max(fp_ctime)=", maxval(fp_ctime)


return
end if


icyc_traffic     = ((1.-real_p)*cyc_traf(int_p)+real_p*cyc_traf(int_p+1))*fp_weekdayload
icyc_basedemand  = ((1.-real_p)*cyc_base(int_p)+real_p*cyc_base(int_p+1))*fp_weekdayload
icyc_proportion  = ((1.-real_p)*cyc_prop(int_p)+real_p*cyc_prop(int_p+1))*fp_weekdayload
icyc_translation = ((1.-real_p)*cyc_tran(int_p)+real_p*cyc_tran(int_p+1))

! remove statistical energy use diurnal adjustments
if (statsmeth==0) then
  icyc_basedemand=1.
  icyc_proportion=1.
  icyc_translation=0.
end if

return
end subroutine getdiurnal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define heating and cooling for different interior model tiles

subroutine uclem_prepinterior(ifrac,fp,pd,mp)

implicit none

integer, intent(in) :: ifrac, mp
type(fparmdata),   intent(in) :: fp
type(pdiagdata),   intent(inout) :: pd
real, dimension(mp) :: cyc_traffic, cyc_basedemand, cyc_proportion, cyc_translation

! fractional interior spaces
call getdiurnal(fp%ctime,fp%weekdayload,cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation)
select case( intairtmeth + (ifrac-1)*10)
  case(0,1)
    pd%frac_sigma = 1.
    pd%ac_heat_on = 1.
    pd%ac_cool_on = 1.
  case(2)  ! heating and cooling
    pd%frac_sigma = min( fp%heatprop, fp%coolprop )*cyc_proportion
    pd%ac_heat_on = 1.
    pd%ac_cool_on = 1.
  case(12) ! heating or cooling
    pd%frac_sigma = (max( fp%heatprop, fp%coolprop ) - min( fp%heatprop, fp%coolprop ))*cyc_proportion
    where ( fp%heatprop>fp%coolprop )
      pd%ac_heat_on = 1.
      pd%ac_cool_on = 0.
    elsewhere
      pd%ac_heat_on = 0.
      pd%ac_cool_on = 1.
    end where
  case(22) ! no conditioning
    pd%frac_sigma = 1. - max( fp%heatprop, fp%coolprop )*cyc_proportion
    pd%ac_heat_on = 0.
    pd%ac_cool_on = 0.
  case default
    write(6,*) "ERROR: Invalid choice for intairtemth",intairtmeth + (ifrac-1)*10
    stop
end select

return
end subroutine uclem_prepinterior

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version holds the eddy size in the canyon to a constant width
! that allows for an analytic solution to the integral

subroutine getincanwindb(ueast,uwest,ufloor,a_udir,z0,fp,mp)

implicit none

integer, intent(in) :: mp
real, dimension(mp), intent(out) :: ueast,uwest,ufloor
real, dimension(mp), intent(in) :: z0
real, dimension(mp) :: wsuma,wsumb,fsum
real, dimension(mp) :: theta1,wdir,h,w
real, dimension(mp) :: dufa,dura,duva,ntheta
real, dimension(mp) :: dufb,durb,duvb
real, dimension(mp) :: fp_coeffbldheight
real, dimension(mp), intent(in) :: a_udir
type(fparmdata), intent(in) :: fp

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

fp_coeffbldheight = fp%effhwratio/fp%hwratio

h=fp%bldheight*fp_coeffbldheight
w=fp%bldheight/fp%hwratio
theta1=acos(min(w/(3.*h),1.))

call winda(dufa,dura,duva,h,w,z0,mp) ! jet on road
call windb(dufb,durb,duvb,h,w,z0,mp) ! jet on wall
ntheta=2. ! i.e., int_0^pi sin(theta) dtheta = 2.)
where (wdir<theta1.or.wdir>pi-theta1) ! jet on wall
  wsuma=duvb*ntheta
  wsumb=durb*ntheta
  fsum=dufb*ntheta
elsewhere                             ! jet on road
  wsuma=dura*ntheta
  wsumb=duva*ntheta
  fsum=dufa*ntheta
end where

! Correct for rotation of winds at start of subroutine
! 0.5 to adjust for factor of 2 in gettopu
where (a_udir>=0.)
  ueast=0.5*wsuma
  uwest=0.5*wsumb
elsewhere
  ueast=0.5*wsumb
  uwest=0.5*wsuma
end where
ufloor=0.5*fsum      ! floor

return
end subroutine getincanwindb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on road

subroutine winda(uf,ur,uv,h,w,z0,mp)

implicit none

integer, intent(in) :: mp
integer iqu
real, dimension(mp), intent(out) :: uf,ur,uv
real, dimension(mp), intent(in) :: h,w,z0
real a,u0,cuven,zolog

u0=exp(-0.9*sqrt(13./4.))
do iqu = 1,mp
  a=0.15*max(1.,3.*h(iqu)/(2.*w(iqu)))
  zolog=log(max(h(iqu),z0(iqu)+0.2)/z0(iqu))
  cuven=log(max(refheight*h(iqu),z0(iqu)+0.2)/z0(iqu))/log(max(h(iqu),z0(iqu)+0.2)/z0(iqu))
  cuven=max(cuven*max(1.-3.*h(iqu)/w(iqu),0.),(u0/a)*(h(iqu)/w(iqu))*(1.-exp(max(-a*max(w(iqu)/h(iqu)-3.,0.),-40.))))
  uf(iqu)=(u0/a)*(h(iqu)/w(iqu))*(1.-exp(-3.*a))+cuven
  !uf(iqu)=(u0/a)*(h(iqu)/w(iqu))*(2.-exp(-a*3.)-exp(-a*(w(iqu)/h(iqu)-3.)))
  ur(iqu)=(u0/a)*exp(-a*3.)*(1.-exp(-a))
  ! MJT suggestion
  cuven=1.-1./zolog
  uv(iqu)=(u0/a)*exp(max(-a*max(w(iqu)/h(iqu)-3.,0.),-40.))*(1.-exp(-a))
  uv(iqu)=max(cuven,uv(iqu))
  !uv(iqu)=(u0/a)*exp(-a*(w(iqu)/h(iqu)-3.))*(1.-exp(-a))
end do

return
end subroutine winda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on wall

subroutine windb(uf,ur,uv,h,win,z0,mp)

implicit none

integer, intent(in) :: mp
integer iqu
real, dimension(mp), intent(out) :: uf,ur,uv
real, dimension(mp), intent(in) :: h,win,z0
real a,dh,u0,w
real zolog,cuven

do iqu = 1,mp
  w=min(win(iqu),1.5*h(iqu))

  a=0.15*max(1.,3.*h(iqu)/(2.*w))
  dh=max(2.*w/3.-h(iqu),0.)
  u0=exp(-0.9*sqrt(13./4.)*dh/h(iqu))

  zolog=log(max(h(iqu),z0(iqu)+0.2)/z0(iqu))
  ! MJT suggestion (cuven is multiplied by dh to avoid divide by zero)
  cuven=h(iqu)-(h(iqu)-dh)*log(max(h(iqu)-dh,z0(iqu)+0.2)/z0(iqu))/zolog-dh/zolog
  ! MJT cuven is back to the correct units of m/s
  cuven=max(cuven/h(iqu),(u0/a)*(1.-exp(-a*dh/h(iqu))))

  uf(iqu)=(u0/a)*(h(iqu)/w)*exp(-a*(1.-dh/h(iqu)))*(1.-exp(-a*w/h(iqu)))
  ur(iqu)=(u0/a)*exp(-a*(1.-dh/h(iqu)+w/h(iqu)))*(1.-exp(-a))
  uv(iqu)=(u0/a)*(1.-exp(-a*(1.-dh/h(iqu))))+cuven
  !uv(iqu)=(u0/a)*(2.-exp(-a*(1.-dh/h(iqu)))-exp(-a*dh/h(iqu)))
end do

return
end subroutine windb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate wind speed at canyon top
subroutine gettopu(d_topu,a_umag,z_on_l,fp_bldheight,pd_cduv,pd_cndzmin,fp_hwratio,mp)

implicit none

integer, intent(in) :: mp
real, dimension(mp), intent(in) :: z_on_l
real, dimension(mp) :: z0_on_l,bldheight
real, dimension(mp) :: pm0,pm1,integralm
real, dimension(mp) :: ustar,neutral
real, dimension(mp), intent(out) :: d_topu
real, dimension(mp), intent(in) :: a_umag
real, dimension(mp), intent(in) :: fp_bldheight, fp_hwratio
real, dimension(mp), intent(inout) :: pd_cduv, pd_cndzmin

bldheight=fp_bldheight*(1.-refheight)
ustar=sqrt(pd_cduv)*a_umag

z0_on_l=min(bldheight,pd_cndzmin)*z_on_l/pd_cndzmin ! calculate at canyon top
z0_on_l=min(z0_on_l,10.)
neutral = log(pd_cndzmin/min(bldheight,pd_cndzmin))
where (z_on_l<0.)
  pm0     = (1.-16.*z0_on_l)**(-0.25)
  pm1     = (1.-16.*z_on_l)**(-0.25)
  integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                +2.*(atan(1./pm1)-atan(1./pm0))
elsewhere
  !-------Beljaars and Holtslag (1991) heat function
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
  pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
  integralm = neutral-(pm1-pm0)
end where
where (bldheight<pd_cndzmin)
  d_topu=(2./pi)*(a_umag-ustar*integralm/vkar)
elsewhere ! within canyon
  d_topu=(2./pi)*a_umag*exp(0.5*fp_hwratio*(1.-pd_cndzmin/bldheight))
end where
d_topu=max(d_topu,0.1)

return
end subroutine gettopu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine scrncalc(a_mixr,a_umag,c_a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,smixr,rdsntemp,zonet, &
                    cnveg,fp,pd,rdhyd,road,mp)

implicit none

integer, intent(in) :: mp
real, dimension(mp), intent(in) :: smixr,rdsntemp,zonet
real, dimension(mp) :: cd,thetav,sthetav
real, dimension(mp) :: thetavstar,z_on_l,z0_on_l
real, dimension(mp) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(mp) :: ustar,qstar,z10_on_l
real, dimension(mp) :: neutral,neutral10,pm10
real, dimension(mp) :: integralm10,tts,tetp
real, dimension(mp) :: tstar,lna
real, dimension(mp) :: utop,ttop,qtop,wf,tsurf,qsurf,n
real, dimension(mp), intent(in) :: a_mixr,a_umag,c_a_temp
real, dimension(mp), intent(in) :: u_ts
real, dimension(mp), intent(in) :: d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd
type(vegdata), intent(in) :: cnveg
type(hydrodata), intent(in) :: rdhyd
type(facetdata), intent(in) :: road
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

real, parameter :: z0  = 1.5
real, parameter :: z10 = 10.

select case(scrnmeth)
  case(0) ! estimate screen diagnostics (slab at displacement height approach)
    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=u_ts + (u_ts+urbtemp)*0.61*smixr
    lna=pd%lzoh-pd%lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,pd%cndzmin,pd%lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tstar=vkar*(c_a_temp-u_ts)/integralh

    z0_on_l  = z0*z_on_l/pd%cndzmin
    z10_on_l = z10*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(pd%cndzmin/z0)
    neutral10 = log(pd%cndzmin/z10)
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh   = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm   = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))     &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    pd%tscrn = c_a_temp - tstar*integralh/vkar
    pd%qscrn = a_mixr - qstar*integralh/vkar
    pd%uscrn = max(a_umag-ustar*integralm/vkar,0.)
    pd%u10   = max(a_umag-ustar*integralm10/vkar,0.)

  case(1) ! estimate screen diagnostics (two step canopy approach)
    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=u_ts + (u_ts+urbtemp)*0.61*smixr
    lna=pd%lzoh-pd%lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,pd%cndzmin,pd%lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tts=vkar*(thetav-sthetav)/integralh
    tstar=vkar*(c_a_temp-u_ts)/integralh

    z0_on_l  = fp%bldheight*(1.-refheight)*z_on_l/pd%cndzmin ! calculate at canyon top
    z10_on_l = max(z10-fp%bldheight*refheight,1.)*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(pd%cndzmin/(fp%bldheight*(1.-refheight)))
    neutral10 = log(pd%cndzmin/max(z10-fp%bldheight*refheight,1.))
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))   &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    ttop = thetav - tts*integralh/vkar
    tetp = c_a_temp - tstar*integralh/vkar
    qtop = a_mixr - qstar*integralh/vkar
    utop = a_umag - ustar*integralm/vkar

    where (fp%bldheight<=z10) ! above canyon
      pd%u10=max(a_umag-ustar*integralm10/vkar,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf = d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*real(road%nodetemp(:,0))+cnveg%sigma*cnveg%temp)
    n=max(min((rdhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)
    wf = (1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta+cnveg%sigma*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(fp%bldheight/zonet)

    thetav=ttop + (ttop+urbtemp)*0.61*qtop
    sthetav=tsurf + (tsurf+urbtemp)*0.61*qsurf
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,utop,fp%bldheight,n,lna,1)
    ustar=sqrt(cd)*utop
    tstar=vkar*(tetp-tsurf)/integralh
    qstar=vkar*(qtop-qsurf)/integralh

    z0_on_l   = z0*z_on_l/fp%bldheight
    z10_on_l  = max(z10,fp%bldheight)*z_on_l/fp%bldheight
    z0_on_l   = min(z0_on_l,10.)
    z10_on_l  = min(z10_on_l,10.)
    neutral   = log(fp%bldheight/z0)
    neutral10 = log(fp%bldheight/max(z10,fp%bldheight))
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))   &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere

    pd%tscrn = tetp-tstar*integralh/vkar
    pd%qscrn = qtop-qstar*integralh/vkar
    pd%uscrn = max(utop-ustar*integralm/vkar,0.)
    where (fp%bldheight>z10) ! within canyon
      pd%u10 = max(utop-ustar*integralm10/vkar,0.)
    end where

  case(2) ! calculate screen diagnostics from canyon only
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*real(road%nodetemp(:,0))+cnveg%sigma*cnveg%temp)
    n=max(min((rdhyd%soilwater-fp%swilt)/(fp%sfc-fp%swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta+cnveg%sigma*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(fp%bldheight/zonet)

    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=tsurf + (tsurf+urbtemp)*0.61*qsurf
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,pd%cndzmin,n,lna,1)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tstar=vkar*(c_a_temp-tsurf)/integralh

    z0_on_l  = z0*z_on_l/pd%cndzmin
    z10_on_l = z10*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(pd%cndzmin/z0)
    neutral10 = log(pd%cndzmin/z10)
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh   = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm   = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))     &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    pd%tscrn = c_a_temp-tstar*integralh/vkar
    pd%qscrn = a_mixr-qstar*integralh/vkar
    pd%uscrn = max(a_umag-ustar*integralm/vkar,0.)
    pd%u10   = max(a_umag-ustar*integralm10/vkar,0.)

end select
pd%qscrn       = max(pd%qscrn,1.E-4)

return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following subroutines are used for internal varying temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates net longwave radiation flux (flux_rg) at each surface
! longwave flux is temperature dependent, so this subroutine should be run at each timestep
subroutine internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle, &
                           skintemps,fp,intl,mp)

implicit none
integer, intent(in) :: mp
integer :: j
real, dimension(mp,4), intent(in) :: skintemps  ! floor, wall, ceiling, wall temperature array
real, dimension(mp,4) :: epsil     ! floor, wall, ceiling, wall emissivity array
real, dimension(mp,4) :: radnet    ! net flux density on ith surface (+ve leaving)
real, dimension(mp)   :: radtot    ! net leaving flux density (B) on ith surface
real, dimension(mp), intent(out) :: rgint_slab,rgint_wallw,rgint_roof,rgint_walle
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp

radnet = 0.
epsil  = 0.9

!epsil = reshape((/(fp_slab%emiss,fp_wall%emiss,fp_roof%emiss,fp_wall%emiss, &
!                    i=1,mp)/), (/mp,4/))

do j = 2,4
  radnet(:,j) = epsil(:,j)/(1.-epsil(:,j))*((sbconst*(skintemps(:,j)+urbtemp)**4)  &
               - sum(intl%psi(:,j,:)*(sbconst*(skintemps(:,:)+urbtemp)**4),dim=2))
end do

radnet(:,1) = -fp%bldheight*(radnet(:,2)+radnet(:,4))/fp%bldwidth - radnet(:,3)

! energy conservation check
radtot(:) = abs(fp%bldwidth(:)*(radnet(:,1)+radnet(:,3)) + fp%bldheight*(radnet(:,2)+radnet(:,4)))

! Restricted to urban columns for the reason given in energyclosure.  This loop
! is per-point rather than guarded by any(), so on the full mp vector it can emit
! one line per non-urban column per timestep per MPI rank.
do j = 1,mp
  if ( fp%sigmau(j)>0. .and. radtot(j)>energytol ) then
    write(6,*) "error: radiation energy non-closure: ", radtot(j)
  end if
end do

rgint_slab  = radnet(:,1)
rgint_wallw = radnet(:,2)
rgint_roof  = radnet(:,3)
rgint_walle = radnet(:,4)

return
end subroutine internal_lwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets internal surface convective heat transfer coefficients
! Compares temperature of innermost layer temperature with air temperature
! Considers horizontal and vertical orientation
! Based on EnergyPlus: Simple Natural Convection Algorithm [W m^-2 K^-1]

subroutine calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, &
                          cvcoeff_intm1,cvcoeff_intm2,roomtemp,roof,slab,intm,mp)
implicit none

integer, intent(in) :: mp
real, dimension(mp), intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw
real, dimension(mp), intent(out) :: cvcoeff_intm1,cvcoeff_intm2,cvcoeff_slab
real, dimension(mp), intent(in)  :: roomtemp
type(facetdata), intent(in) :: roof
type(facetdata), intent(in) :: slab
type(facetdata), intent(in) :: intm
real, parameter :: tdel = 0.1
real, dimension(mp) :: x

select case(cvcoeffmeth)
  case(0) ! DOE Simple
    cvcoeff_walle = 3.067    ! vertical surface coefficient constant
    cvcoeff_wallw = 3.067    ! vertical surface coefficient constant
    where ( real(roof%nodetemp(:,nl))>=roomtemp+tdel )
      cvcoeff_roof(:)=0.948  ! reduced convection
    elsewhere ( real(roof%nodetemp(:,nl))>roomtemp-tdel )
      x = (real(roof%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_roof(:)=(1.-x)*4.04 + x*0.948 ! blend
    elsewhere
      cvcoeff_roof(:)=4.040  ! enhanced convection
    end where
    cvcoeff_intm1 = 3.067   ! vertical surface coefficient constant
    cvcoeff_intm2 = 3.067   ! vertical surface coefficient constant
    where ( real(slab%nodetemp(:,nl))<=roomtemp-tdel )
      cvcoeff_slab(:)=0.7   ! reduced convection
    elsewhere ( real(slab%nodetemp(:,nl))<roomtemp+tdel )
      x = (real(slab%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_slab(:) = (1.-x)*0.7 + x*4.04 ! blend
    elsewhere
      cvcoeff_slab(:)=4.040  ! enhanced convection
    end where
  case(1) ! dynamic, from international standard ISO6946:2007, Annex A
    cvcoeff_walle = 2.5    ! vertical surface coefficient constant
    cvcoeff_wallw = 2.5    ! vertical surface coefficient constant
    where ( real(roof%nodetemp(:,nl))>=roomtemp+tdel )
      cvcoeff_roof(:)=0.7  ! reduced convection (upper surface)
    elsewhere ( real(roof%nodetemp(:,nl))>roomtemp-tdel )
      x = (real(roof%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_roof(:) = (1.-x)*5.0 + x*0.7 ! blend
    elsewhere
      cvcoeff_roof(:)=5.0  ! enhanced convection (upper surface)
    end where
    where ( real(intm%nodetemp(:,nl))>=roomtemp+tdel )
      cvcoeff_intm2 = 0.7+5.7   ! reduced convection (upper surface) + radiation @20 deg C
    elsewhere ( real(intm%nodetemp(:,nl))>roomtemp-tdel )
      x = (real(intm%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_intm2 = (1.-x)*(5.+5.7) + x*(0.7+5.7) ! blend
    elsewhere
      cvcoeff_intm2 = 5.0+5.7   ! reduced convection (upper surface) + radiation @20 deg C
    end where
    where ( real(intm%nodetemp(:,0))<=roomtemp-tdel )
      cvcoeff_intm1 = 0.7+5.7   ! reduced convection (lower surface) + radiation @20 deg C
    elsewhere ( real(intm%nodetemp(:,0))<roomtemp+tdel )
      x = (real(intm%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_intm1 = (1.-x)*(0.7+5.7) + x*(5.+5.7) ! blend
    elsewhere
      cvcoeff_intm1 = 5.0+5.7   ! reduced convection (lower surface) + radiation @20 deg C
    end where
    where ( real(slab%nodetemp(:,nl))<=roomtemp-tdel )
      cvcoeff_slab(:)=0.7   ! reduced convection (lower surface)
    elsewhere ( real(slab%nodetemp(:,nl))<roomtemp+tdel )
      x = (real(slab%nodetemp(:,nl))-roomtemp+tdel)/(2.*tdel)
      cvcoeff_slab = (1.-x)*0.7 + x*5.
    elsewhere
      cvcoeff_slab(:)=5.0   ! enhanced convection (lower surface)
    end where
  case(2) ! fixed, from international standard IS6946:2007, 5.2
    cvcoeff_roof  = 1./0.10 ! fixed coefficient up
    cvcoeff_walle = 1./0.13 ! fixed coefficient horizontal
    cvcoeff_wallw = 1./0.13 ! fixed coefficient horizontal
    cvcoeff_intm1 = 1./0.10 ! fixed coefficient up
    cvcoeff_intm2 = 1./0.17 ! fixed coefficient down
    cvcoeff_slab  = 1./0.17 ! fixed coefficient down
end select

end subroutine calc_convcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the proportion of open windows for ventilation
subroutine calc_openwindows(d_openwindows,fp,iroomtemp,d_canyontemp, &
                            roof,walle,wallw,slab,mp)
implicit none

integer, intent(in)                 :: mp
real, dimension(mp), intent(in)  :: d_canyontemp,iroomtemp
real, dimension(mp), intent(out) :: d_openwindows
real, dimension(mp)              :: xtemp, mrt
type(facetdata), intent(in) :: roof, walle, wallw, slab
type(fparmdata), intent(in) :: fp

! mean radiant temperature estimation
mrt = 0.5*(fp%bldheight/(fp%bldwidth+fp%bldheight)*(real(walle%nodetemp(:,nl) + wallw%nodetemp(:,nl)))) &
    + 0.5*(fp%bldwidth/ (fp%bldwidth+fp%bldheight)*(real(roof%nodetemp(:,nl) + slab%nodetemp(:,nl))))
! globe temperature approximation (average of mrt and air temperature) [Celsius]
xtemp = 0.5*(iroomtemp + mrt)

d_openwindows = 1./(1. + exp( 0.5*(fp%bldairtemp+ac_deltat+5.-xtemp) )) &
    *1./(1. + exp( max(min(d_canyontemp-iroomtemp,40.),-40.) ))

return
end subroutine calc_openwindows

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getnetalbedo(alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,  &
                        fp_hwratio,fp_sigmabld,fp_vegsigmar,fp_roofalpha,fp_vegalphar,          &
                        fp_vegsigmac,fp_roadalpha,fp_wallalpha,fp_vegalphac,                    &
                        roofalpha,roadalpha,snowdeltar,snowdeltac)

implicit none

integer iqu
real, dimension(:), intent(out) :: alb
real, dimension(:), intent(in) :: sg_roof, sg_vegr, sg_road, sg_walle, sg_wallw, sg_vegc
real, dimension(:), intent(in) :: sg_rfsn, sg_rdsn
real, dimension(:), intent(in) :: fp_hwratio, fp_sigmabld
real, dimension(:), intent(in) :: fp_vegsigmar, fp_roofalpha, fp_vegalphar
real, dimension(:), intent(in) :: fp_vegsigmac, fp_roadalpha, fp_vegalphac, fp_wallalpha
real, dimension(:), intent(in) :: roofalpha, roadalpha, snowdeltar, snowdeltac
real albu, albr

do iqu = 1,size(alb)
  ! canyon
  albu=1.-(fp_hwratio(iqu)*(sg_walle(iqu)+sg_wallw(iqu))*(1.-fp_wallalpha(iqu))         &
      +snowdeltac(iqu)*sg_rdsn(iqu)*(1.-roadalpha(iqu))                                 &
      +(1.-snowdeltac(iqu))*((1.-fp_vegsigmac(iqu))*sg_road(iqu)*(1.-fp_roadalpha(iqu)) &
      +fp_vegsigmac(iqu)*sg_vegc(iqu)*(1.-fp_vegalphac(iqu))))
  ! roof
  albr=(1.-snowdeltar(iqu))*((1.-fp_vegsigmar(iqu))*sg_roof(iqu)*fp_roofalpha(iqu) &
      +fp_vegsigmar(iqu)*sg_vegr(iqu)*fp_vegalphar(iqu))                           &
      +snowdeltar(iqu)*sg_rfsn(iqu)*roofalpha(iqu)
  ! net
  alb(iqu)=fp_sigmabld(iqu)*albr+(1.-fp_sigmabld(iqu))*albu
end do

return
end subroutine getnetalbedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine init_lwcoeff(fp,intl,mp)
! This subroutine calculates longwave reflection coefficients (int_psi) at each surface
! longwave coefficients do not change, so this subroutine should only be run once
! Infinite reflections per Harman et al., (2004) "Radiative Exchange in Urban Street Canyons"
! Per method in "Radiation Heat Transfer, Sparrow & Cess 1978, Ch 3-3"
! array surface order is: (1) floor; (2) wallw; (3) ceiling; (4) walle

! local variables
integer, intent(in) :: mp
real, dimension(mp,4,4) :: chi
real, dimension(4,4)       :: krondelta
real, dimension(mp)     :: h, w
real, dimension(mp,4)   :: epsil   ! floor, wall, ceiling, wall emissivity array
integer :: i, j
integer :: ierr       ! inverse matrix error flag
type(intldata), intent(inout) :: intl
type(fparmdata), intent(in) :: fp


krondelta = 0.
chi = 0.
intl%psi = 0.
h = fp%bldheight
w = fp%sigmabld*(fp%bldheight/fp%hwratio)/(1.-fp%sigmabld)

! set int_vfactors
intl%viewf(:,1,1) = 0.                                    ! floor to self
intl%viewf(:,1,2) = 0.5*(1.+(h/w)-sqrt(1.+(h/w)**2))      ! floor to wallw
intl%viewf(:,1,3) = sqrt(1.+(h/w)**2)-(h/w)               ! floor to ceiling
intl%viewf(:,1,4) = intl%viewf(:,1,2)                     ! floor to walle
intl%viewf(:,2,1) = 0.5*(1.+(w/h)-sqrt(1.+(w/h)**2))      ! wallw to floor
intl%viewf(:,2,2) = 0.                                    ! wallw to self
intl%viewf(:,2,3) = intl%viewf(:,2,1)                     ! wallw to ceiling
intl%viewf(:,2,4) = sqrt(1.+(w/h)**2)-(w/h)               ! wallw to walle
intl%viewf(:,3,1) = intl%viewf(:,1,3)                     ! ceiling to floor
intl%viewf(:,3,2) = intl%viewf(:,1,2)                     ! ceiling to wallw
intl%viewf(:,3,3) = 0.                                    ! ceiling to self
intl%viewf(:,3,4) = intl%viewf(:,1,2)                     ! ceiling walle
intl%viewf(:,4,1) = intl%viewf(:,2,1)                     ! walle to floor
intl%viewf(:,4,2) = intl%viewf(:,2,4)                     ! walle to wallw
intl%viewf(:,4,3) = intl%viewf(:,2,1)                     ! walle to ceiling
intl%viewf(:,4,4) = 0.                                    ! walle to self

!epsil = reshape((/(f_slab%emiss,f_wall%emiss,f_roof%emiss,f_wall%emiss, &
!                    i=1,mp_g)/), (/mp_g,4/))
epsil = 0.9

do i = 1,4
  krondelta(i,i) = 1.
end do
do j = 1,4
  do i = 1,4
    chi(:,i,j) = (krondelta(i,j) - (1.-epsil(:,i))*intl%viewf(:,i,j))/(epsil(:,i))
  end do
end do

! invert matrix
intl%psi = chi
call minverse(intl%psi,ierr)

end subroutine init_lwcoeff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises internal variables;
! building width and number of internal mass floors
subroutine init_internal(fp)

implicit none

type(fparmdata), intent(inout) :: fp
real :: floor_height = 3. ! assumed average floor to floor height (m)

fp%bldwidth = fp%sigmabld*(fp%bldheight/fp%hwratio)/(1.-fp%sigmabld)

! define number of internal mass floors (based on building height)
select case(intmassmeth)
  case(0) ! no internal mass
    fp%intmassn = 0
  case(1) ! one floor of internal mass
    fp%intmassn = 1
  case(2) ! dynamic floors of internal mass
    fp%intmassn = max((nint(fp%bldheight/floor_height)-1),0)
end select

end subroutine init_internal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the inverse of a NxN matrix
! input/output = a, size=s, error flag=ier

subroutine minverse(a_inout,ierr)

implicit none

real, dimension(:,:,:), intent(inout) :: a_inout
real(kind=8), dimension(size(a_inout,1),size(a_inout,2),size(a_inout,3)) :: a
real(kind=8), dimension(size(a_inout,1)) :: det, d, amax
real(kind=8), dimension(size(a_inout,2)) :: x
real(kind=8) :: y
integer, intent(out)  :: ierr
integer s, ns, iq, i, j, nu
integer, dimension(size(a_inout,1),size(a_inout,2)) :: row, col
integer, dimension(size(a_inout,1)) :: prow, pcol
logical, dimension(size(a_inout,1),size(a_inout,2)) :: notpiv

a = real(a_inout,8)

nu = size(a_inout,1)
s = size(a_inout,2)

det = 0.
d = 1.
notpiv = .TRUE.

do ns = 1,s

  amax(:) = 0.
  do j = 1,s
    do i = 1,s
      where ( notpiv(:,j) .and. notpiv(:,i) .and. amax(:)<abs(a(:,i,j)) )
        amax(:) = abs(a(:,i,j))
        prow(:) = i
        pcol(:) = j
      end where
    end do
  end do

  if ( any(amax<0.) ) then
    ierr = 1 ! cannot invert matrix
    return
  end if

  do iq = 1,nu
    notpiv(iq,pcol(iq)) = .FALSE.
    if ( prow(iq)/=pcol(iq) ) then
      d(iq) = -d(iq)
      x(1:s) = a(iq,prow(iq),1:s)
      a(iq,prow(iq),1:s) = a(iq,pcol(iq),1:s)
      a(iq,pcol(iq),1:s) = x(1:s)
    end if
  end do

  row(:,ns) = prow(:)
  col(:,ns) = pcol(:)
  do iq = 1,nu
    amax(iq) = a(iq,pcol(iq),pcol(iq))
  end do
  d(:) = d(:)*amax(:)

  if ( any(abs(d)<=0.) ) then
    ierr = 1 ! cannot invert matrix
    return
  end if

  amax(:) = 1./amax(:)
  do iq = 1,nu
    a(iq,pcol(iq),pcol(iq)) = 1.
    a(iq,pcol(iq),1:s) = a(iq,pcol(iq),1:s)*amax(iq)
  end do

  do i=1,s
    do iq = 1,nu
      if ( i/=pcol(iq) ) then
        y = a(iq,i,pcol(iq))
        a(iq,i,pcol(iq)) = 0.
        do j = 1,s
          a(iq,i,j) = a(iq,i,j) - y*a(iq,pcol(iq),j)
        end do
      end if
    end do
  end do

end do

det(:) = d(:)

do ns = s,1,-1
  prow(:) = row(:,ns)
  pcol(:) = col(:,ns)
  do iq = 1,nu
    if ( prow(iq)/=pcol(iq) ) then
      do i = 1,s
        y = a(iq,i,prow(iq))
        a(iq,i,prow(iq)) = a(iq,i,pcol(iq))
        a(iq,i,pcol(iq)) = y
      end do
    end if
  end do
end do

a_inout = real(a)

ierr = 0

return
end subroutine minverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo
!

subroutine uclem_calc_alb(alb,fp_hwratio,fp_effhwratio,fp_vangle,      &
                          fp_hangle,fp_fbeam,                          &
                          fp_vegsigmac,fp_vegsigmar,fp_sigmabld,       &
                          fp_roadalpha,fp_wallalpha,fp_roofalpha,      &
                          fp_vegalphac,fp_vegalphar,                   &
                          road_snow_alpha,roof_snow_alpha,             &
                          can_snow,roof_snow)

implicit none

real, dimension(:), intent(in) :: fp_hwratio, fp_effhwratio, fp_vangle, fp_hangle, fp_fbeam
real, dimension(:), intent(in) :: fp_roadalpha, fp_wallalpha, fp_roofalpha
real, dimension(:), intent(in) :: fp_vegalphar, fp_vegalphac
real, dimension(:), intent(in) :: road_snow_alpha, roof_snow_alpha
real, dimension(:), intent(in) :: can_snow, roof_snow
real, dimension(:), intent(in) :: fp_sigmabld, fp_vegsigmar, fp_vegsigmac
real, dimension(:), intent(out) :: alb
real, dimension(size(alb)) :: snowdeltar, snowdeltac, wallpsi, roadpsi
real, dimension(size(alb)) :: sg_roof, sg_vegr, sg_road, sg_walle, sg_wallw, sg_vegc, sg_rfsn, sg_rdsn
real, dimension(size(alb)) :: fp_coeffbldheight

fp_coeffbldheight = fp_effhwratio/fp_hwratio

! roof
snowdeltar = roof_snow/(roof_snow+maxrfsn)
! canyon
snowdeltac = can_snow/(can_snow+maxrdsn)

! call the subroutine
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,       &
                sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,         &
                fp_hwratio,fp_effhwratio,fp_vangle,fp_hangle,    &
                fp_fbeam,fp_vegsigmac,fp_roadalpha,fp_vegalphac, &
                fp_wallalpha,road_snow_alpha,snowdeltac)
sg_walle = sg_walle*fp_coeffbldheight
sg_wallw = sg_wallw*fp_coeffbldheight

call getnetalbedo(alb,sg_roof,sg_vegr,sg_road,sg_walle,          &
                  sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,              &
                  fp_hwratio,fp_sigmabld,fp_vegsigmar,           &
                  fp_roofalpha,fp_vegalphar,fp_vegsigmac,        &
                  fp_roadalpha,fp_wallalpha,fp_vegalphac,        &
                  roof_snow_alpha,road_snow_alpha,               &
                  snowdeltar,snowdeltac)

return
end subroutine uclem_calc_alb

end module uclem_science
