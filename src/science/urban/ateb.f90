! UCLEM urban canopy model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! The Urban CLimate Energy Model (UCLEM) results when the exterior canyon model (aTEB) is combined with an interior model
! of energy use from Lipson, et al (2018).  This model includes a representation of energy storage and conduction for the
! interior slab, prognostic internal air temperature, as well as dynamic anthropogenic energy generation as a function
! of internal air temperature and a representation of human behaviours over the diurnal cycle based on Thatcher (2007).


! Usual practice is:
!   call atebinit            ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call atebloadm           ! to load previous state arrays (from tebsavem)
!   call atebtype            ! to define urban type (or use tebfndef to define urban properties at each grid point)
!   ...
!   do t=1,tmax
!     ...
!     call atebnewangle1     ! store current solar zenith and azimuthal angle (use atebccangle for CCAM)
!     call atebalb1(split=1) ! returns urban albedo for direct component of shortwave radiation
!     call atebalb1(split=2) ! returns urban albedo for diffuse component of shortwave radiation
!     ...
!     call atebfbeam         ! store fraction of direct shortwave radiation (or use atebspitter to estimate fraction)
!     call atebalb1          ! returns net urban albedo (i.e., default split=0)
!     ...
!     call atebcalc          ! calculates urban temperatures, fluxes, etc and blends with input
!     call atebcd            ! returns urban drag coefficient (or use atebzo for roughness length)
!     call atebscrnout       ! returns screen level diagnostics
!     ...
!   end do
!   ...
!   call atebsavem           ! to save current state arrays (for use by tebloadm)
!   call atebend             ! to deallocate memory before quitting

! only atebinit and atebcalc are mandatory.  All other subroutine calls are optional.


! DEPENDICES:

! atebalb1(split=1)                   depends on     atebnewangle1 (or atebccangle)
! atebalb1(split=2)                   depends on     atebnewangle1 (or atebccangle)
! atebalb1 (i.e., default split=0)    depends on     atebnewangle1 (or atebccangle) and atebfbeam (or atebspitter)
! atebcalc                            depends on     atebnewangle1 (or atebccangle) and atebfbeam (or atebspitter)  
! atebcd (or atebzo)                  depends on     atebcalc
! atebscrnout                         depends on     atebcalc


! URBAN TYPES:
 
! 1 = Urban               (TAPM 31)
! 2 = Urban (low)         (TAPM 32)
! 3 = Urban (medium)      (TAPM 33)
! 4 = Urban (high)        (TAPM 34)
! 5 = Urban (cbd)         (TAPM 35)
! 6 = Industrial (low)    (TAPM 36)
! 7 = Industrial (medium) (TAPM 37)
! 8 = Industrial (high)   (TAPM 38)

module ateb

!#ifdef CCAM
!use cc_omp, only : imax, ntiles
!#endif

implicit none

private
public atebinit,atebcalc,atebend,atebzo,atebtype,atebalb1,                             &
       atebnewangle1,atebccangle,atebdisable,atebcd,                                   &
       atebdwn,atebscrnout,atebfbeam,atebspitter,atebsigmau,energyrecord,atebdeftype,  &
       atebhydro,atebenergy,atebloadd,atebsaved,atebmisc,atebdeftype_export,           &
       atebavetemp
public atebnmlfile,urbtemp,energytol,resmeth,zohmeth,acmeth,nrefl,vegmode,             &
       soilunder,scrnmeth,wbrelaxc,wbrelaxr,lweff,ncyits,nfgits,tol,alpha,             &
       zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden,refheight,     &
       zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf,        &
       intairtmeth,intmassmeth,statsmeth,lwintmeth,cvcoeffmeth,infilmeth,acfactor,     &
       ac_heatcap,ac_coolcap,ac_deltat,cyc_prop,cyc_base,cyc_traf,cyc_tran


!#ifdef CCAM
!public upack_g,ufull_g,nl,nfrac
!public f_roof,f_wall,f_road,f_slab,f_intm
!public intm_g,rdhyd_g,rfhyd_g,rfveg_g
!public road_g,roof_g,room_g,slab_g,walle_g,wallw_g,cnveg_g,intl_g
!public f_g,p_g
!public facetparams,facetdata,hydrodata,vegdata,intldata
!public fparmdata,pdiagdata
!#endif

! state arrays
integer, save :: ifull, nfrac
integer, save :: ntiles = 1     ! Emulate OMP
integer, save :: imax = 0       ! Emulate OMP
integer, dimension(:), allocatable, save :: ufull_g
logical, save :: ateb_active = .false.
logical, dimension(:,:), allocatable, save :: upack_g
real, dimension(:,:), allocatable, save :: atebdwn ! This array is for CCAM onthefly.f
real, dimension(0:220), save :: table

type facetdata
  real, dimension(:,:), allocatable :: nodetemp        ! Temperature of node (prognostic)       [K]
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
  real, dimension(:), allocatable :: hwratio,coeffbldheight,effhwratio,sigmabld
  real, dimension(:), allocatable :: industryfg,intgains_flr,trafficfg,bldheight,bldwidth
  real, dimension(:), allocatable :: ctime,vangle,hangle,fbeam,weekdayload
  real, dimension(:), allocatable :: bldairtemp
  real, dimension(:), allocatable :: swilt,sfc,ssat,rfvegdepth
  real, dimension(:), allocatable :: infilach,ventilach,heatprop,coolprop
  real, dimension(:), allocatable :: sigmau
  integer, dimension(:), allocatable :: intmassn
end type fparmdata

type pdiagdata
  real, dimension(:), allocatable :: lzom, lzoh, cndzmin, cduv, cdtq
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
  real(kind=8), dimension(:), allocatable :: surferr, atmoserr, surferr_bias, atmoserr_bias
  real(kind=8), dimension(:), allocatable :: storage_flux
  ! real(kind=8), dimension(:,:), allocatable :: storagetot_road, storagetot_walle, storagetot_wallw, storagetot_roof
  logical :: first_call
end type pdiagdata

type(facetdata), dimension(:,:), allocatable, save :: roof_g, road_g, walle_g, wallw_g, slab_g, intm_g, room_g
type(hydrodata), dimension(:,:), allocatable, save :: rfhyd_g, rdhyd_g
type(pdiagdata), dimension(:,:), allocatable, save :: p_g
type(facetparams), dimension(:), allocatable, save :: f_roof, f_road, f_wall, f_slab, f_intm
type(vegdata), dimension(:), allocatable,     save :: cnveg_g, rfveg_g
type(intldata), dimension(:), allocatable,    save :: intl_g
type(fparmdata), dimension(:), allocatable,   save :: f_g


! model parameters
integer, save      :: atebnmlfile=11       ! Read configuration from nml file (0=off, >0 unit number (default=11))
integer, save      :: resmeth=1            ! Canyon sensible heat transfer (0=Masson, 1=Harman (varying width), 2=Kusaka,
                                           ! 3=Harman (fixed width))
integer, save      :: zohmeth=1            ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, save      :: acmeth=1             ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, 3=COP of 1.0)
integer, save      :: intairtmeth=0        ! Internal air temperature (0=prescribed, 1=aggregated varying, 2=fractional varying)
integer, save      :: intmassmeth=2        ! Internal thermal mass (0=none, 1=one floor, 2=dynamic floor number)
integer, save      :: cvcoeffmeth=1        ! Internal surface convection heat transfer coefficient (0=DOE, 1=ISO6946, 2=fixed)
integer, save      :: statsmeth=1          ! Use statistically based diurnal QF amendments (0=off, 1=on) from Thatcher 2007 
integer, save      :: infilmeth=1          ! Method to calculate infiltration rate (0=constant, 1=EnergyPlus/BLAST, 2=ISO)
integer, save      :: nrefl=3              ! Number of canyon reflections for radiation (default=3)
integer, save      :: vegmode=2            ! In-canyon vegetation mode (0=50%/50%, 1=100%/0%, 2=0%/100%, where out/in=X/Y.
                                           ! Negative values are X=abs(vegmode))
integer, save      :: soilunder=1          ! Modify road heat capacity to extend under
                                           ! (0=road only, 1=canveg, 2=bld, 3=canveg & bld)
integer, save      :: scrnmeth=1           ! Screen diagnostic method (0=Slab, 1=Hybrid, 2=Canyon)
integer, save      :: wbrelaxc=0           ! Relax canyon soil moisture for irrigation (0=Off, 1=On)
integer, save      :: wbrelaxr=0           ! Relax roof soil moisture for irrigation (0=Off, 1=On)
integer, save      :: lweff=2              ! Modification of LW flux for effective canyon height (0=insulated, 1=coupled, 2=full)
integer, save      :: lwintmeth=2          ! Interior longwave flux method (0=off, 1=explicit, 2=with forward step facet temps)
integer, parameter :: nl=4                 ! Number of layers (default 4, must be factors of 4)
integer, save      :: iqt=314              ! Diagnostic point (in terms of host grid)
! sectant solver parameters
integer, save      :: ncyits=6             ! Number of iterations for balancing canyon sensible and latent heat fluxes (default=6)
integer, save      :: nfgits=3             ! Number of iterations for balancing veg and snow energy budgets (default=3)
real, save         :: tol=0.001            ! Sectant method tolerance for sensible heat flux (default=0.001)
real, save         :: alpha=1.             ! Weighting for determining the rate of convergence when calculating canyon temperatures
real, save         :: infilalpha=0.5       ! Weighting for dampening factor when calculating internal infiltration
real(kind=8), save :: energytol=0.05_8    ! Tolerance for acceptable energy closure in each timestep
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

interface atebcalc
  module procedure atebcalc_standard, atebcalc_thread
end interface
  
interface atebenergy
  module procedure atebenergy_standard, atebenergy_thread
end interface

interface atebzo
  module procedure atebzo_standard, atebzo_thread
end interface

interface atebcd
  module procedure atebcd_standard, atebcd_thread
end interface

interface atebhydro
  module procedure atebhydro_standard, atebhydro_thread
end interface

interface atebtype
  module procedure atebtype_standard, atebtype_thread
end interface

interface atebmisc
  module procedure atebmisc_standard, atebmisc_thread
end interface










contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepare the arrays used by the aTEB scheme
! This is a compulsory subroutine that must be called during
! model initialisation

subroutine atebinit(ifin,sigu,diag)

implicit none

integer, intent(in) :: ifin,diag
integer, dimension(:), allocatable, save :: utype
integer tile, is, ie, ifrac
real, dimension(ifin), intent(in) :: sigu

if (diag>=1) write(6,*) "Initialising aTEB"

ateb_active = .true.

ifull=ifin

if ( ntiles<1 ) then
  write(6,*) "ERROR: Invalid ntiles ",ntiles
  stop
end if

imax = ifull/ntiles
if ( mod(ifull,ntiles)/=0 ) then
  write(6,*) "ERROR: Invalid ntiles ",ntiles," for ifull ",ifull
  stop
end if

allocate( f_roof(ntiles), f_road(ntiles), f_wall(ntiles), f_slab(ntiles), f_intm(ntiles) )
allocate( cnveg_g(ntiles), rfveg_g(ntiles) )
allocate( f_g(ntiles) )
allocate( intl_g(ntiles) )
allocate( ufull_g(ntiles) )
allocate( upack_g(imax,ntiles) )

allocate( utype(imax) )

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  upack_g(1:imax,tile) = sigu(is:ie)>0.
  ufull_g(tile) = count( upack_g(1:imax,tile) )

  allocate(f_g(tile)%rfvegdepth(ufull_g(tile)))
  allocate(f_g(tile)%ctime(ufull_g(tile)),f_g(tile)%bldairtemp(ufull_g(tile)),f_g(tile)%weekdayload(ufull_g(tile)))
  allocate(f_g(tile)%hangle(ufull_g(tile)),f_g(tile)%vangle(ufull_g(tile)),f_g(tile)%fbeam(ufull_g(tile)))
  allocate(f_g(tile)%hwratio(ufull_g(tile)),f_g(tile)%coeffbldheight(ufull_g(tile)))
  allocate(f_g(tile)%effhwratio(ufull_g(tile)),f_g(tile)%bldheight(ufull_g(tile)))
  allocate(f_g(tile)%sigmabld(ufull_g(tile)),f_g(tile)%industryfg(ufull_g(tile)))
  allocate(f_g(tile)%intgains_flr(ufull_g(tile)),f_g(tile)%trafficfg(ufull_g(tile)))
  allocate(f_g(tile)%swilt(ufull_g(tile)),f_g(tile)%sfc(ufull_g(tile)),f_g(tile)%ssat(ufull_g(tile)))
  allocate(f_g(tile)%heatprop(ufull_g(tile)),f_g(tile)%coolprop(ufull_g(tile)))
  allocate(f_g(tile)%intmassn(ufull_g(tile)),f_g(tile)%bldwidth(ufull_g(tile)))
  allocate(f_g(tile)%infilach(ufull_g(tile)),f_g(tile)%ventilach(ufull_g(tile)))
  allocate(f_g(tile)%sigmau(ufull_g(tile)))
  
  allocate(f_roof(tile)%depth(ufull_g(tile),nl),f_roof(tile)%lambda(ufull_g(tile),nl))
  allocate(f_roof(tile)%volcp(ufull_g(tile),nl))
  allocate(f_wall(tile)%depth(ufull_g(tile),nl),f_wall(tile)%lambda(ufull_g(tile),nl))
  allocate(f_wall(tile)%volcp(ufull_g(tile),nl))
  allocate(f_road(tile)%depth(ufull_g(tile),nl),f_road(tile)%lambda(ufull_g(tile),nl))
  allocate(f_road(tile)%volcp(ufull_g(tile),nl))
  allocate(f_slab(tile)%depth(ufull_g(tile),nl),f_slab(tile)%lambda(ufull_g(tile),nl))
  allocate(f_slab(tile)%volcp(ufull_g(tile),nl))
  allocate(f_intm(tile)%depth(ufull_g(tile),nl),f_intm(tile)%lambda(ufull_g(tile),nl))
  allocate(f_intm(tile)%volcp(ufull_g(tile),nl))
  allocate(f_roof(tile)%emiss(ufull_g(tile)),f_roof(tile)%alpha(ufull_g(tile)))
  allocate(f_wall(tile)%emiss(ufull_g(tile)),f_wall(tile)%alpha(ufull_g(tile)))
  allocate(f_road(tile)%emiss(ufull_g(tile)),f_road(tile)%alpha(ufull_g(tile)))
  allocate(f_slab(tile)%emiss(ufull_g(tile)))
 
  ! veg temperature is recalculated each time-step and hence does not need to be stored for each ifrac
  allocate(cnveg_g(tile)%emiss(ufull_g(tile)),cnveg_g(tile)%sigma(ufull_g(tile)))
  allocate(cnveg_g(tile)%alpha(ufull_g(tile)),cnveg_g(tile)%rsmin(ufull_g(tile)))
  allocate(cnveg_g(tile)%zo(ufull_g(tile)),cnveg_g(tile)%lai(ufull_g(tile)))
  allocate(rfveg_g(tile)%emiss(ufull_g(tile)),rfveg_g(tile)%sigma(ufull_g(tile)))
  allocate(rfveg_g(tile)%alpha(ufull_g(tile)),rfveg_g(tile)%rsmin(ufull_g(tile)))
  allocate(rfveg_g(tile)%zo(ufull_g(tile)),rfveg_g(tile)%lai(ufull_g(tile)))
  allocate(rfveg_g(tile)%temp(ufull_g(tile)),cnveg_g(tile)%temp(ufull_g(tile)))
  
  allocate(intl_g(tile)%viewf(ufull_g(tile),4,4),intl_g(tile)%psi(ufull_g(tile),4,4))
  
  if ( ufull_g(tile)>0 ) then
      
    ! define grid arrays
    f_g(tile)%sigmau = pack(sigu(is:ie),upack_g(1:imax,tile))
    f_g(tile)%rfvegdepth=0.1
    f_g(tile)%hwratio=1.
    f_g(tile)%sigmabld=0.5
    f_g(tile)%industryfg=0.
    f_g(tile)%intgains_flr=0.
    f_g(tile)%trafficfg=0.
    f_g(tile)%bldheight=10.
    f_g(tile)%bldairtemp=1. ! + urbtemp
    f_g(tile)%vangle=0.
    f_g(tile)%hangle=0.
    f_g(tile)%ctime=0.
    f_g(tile)%fbeam=1.
    f_g(tile)%swilt=0.
    f_g(tile)%sfc=0.5
    f_g(tile)%ssat=1.
    f_g(tile)%infilach=0.5
    f_g(tile)%ventilach=2.
    f_g(tile)%weekdayload=1.0
    
    f_roof(tile)%depth=0.1
    f_roof(tile)%volcp=2.E6
    f_roof(tile)%lambda=2.
    f_roof(tile)%alpha=0.2
    f_roof(tile)%emiss=0.97
    f_wall(tile)%depth=0.1
    f_wall(tile)%volcp=2.E6
    f_wall(tile)%lambda=2.
    f_wall(tile)%alpha=0.2
    f_wall(tile)%emiss=0.97
    f_road(tile)%depth=0.1
    f_road(tile)%volcp=2.E6
    f_road(tile)%lambda=2.
    f_road(tile)%alpha=0.2
    f_road(tile)%emiss=0.97
    f_slab(tile)%depth=0.1
    f_slab(tile)%volcp=2.E6
    f_slab(tile)%lambda=2.
    f_slab(tile)%emiss=0.97
    f_intm(tile)%depth=0.1
    f_intm(tile)%lambda=2.
    f_intm(tile)%volcp=2.E6
    
    cnveg_g(tile)%sigma=0.5
    cnveg_g(tile)%alpha=0.2
    cnveg_g(tile)%emiss=0.97
    cnveg_g(tile)%zo=0.1
    cnveg_g(tile)%lai=1.
    cnveg_g(tile)%rsmin=200.
    cnveg_g(tile)%temp=1. ! + urbtemp             ! updated in atebcalc
    rfveg_g(tile)%sigma=0.
    rfveg_g(tile)%alpha=0.2
    rfveg_g(tile)%emiss=0.97
    rfveg_g(tile)%zo=0.1
    rfveg_g(tile)%lai=1.
    rfveg_g(tile)%rsmin=200.
    rfveg_g(tile)%temp=1. ! + urbtemp             ! updated in atebcalc
    
    utype=1 ! default urban
    call atebtype(utype,diag,f_g(tile),cnveg_g(tile),rfveg_g(tile),      &
                  f_roof(tile),f_road(tile),f_wall(tile),f_slab(tile),   &
                  f_intm(tile),intl_g(tile),upack_g(:,tile),             &
                  ufull_g(tile))
    
  end if
    
end do

deallocate( utype )


! moved here so atebtype can read namelist
select case( intairtmeth )
  case(0,1)
    nfrac = 1
  case(2)
    nfrac = 3
  case default
    write(6,*) "ERROR: Unknown intairtmeth option ",intairtmeth
    stop
end select

allocate( roof_g(nfrac,ntiles), road_g(nfrac,ntiles), walle_g(nfrac,ntiles), wallw_g(nfrac,ntiles) )
allocate( slab_g(nfrac,ntiles), intm_g(nfrac,ntiles) )
allocate( room_g(nfrac,ntiles) )
allocate( rfhyd_g(nfrac,ntiles), rdhyd_g(nfrac,ntiles) )
allocate( p_g(nfrac,ntiles) )

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax

  upack_g(1:imax,tile) = sigu(is:ie)>0.
  ufull_g(tile) = count( upack_g(1:imax,tile) )
    
  do ifrac = 1,nfrac
  
    allocate(roof_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl),road_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl))
    allocate(walle_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl),wallw_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl))
    allocate(slab_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl),intm_g(ifrac,tile)%nodetemp(ufull_g(tile),0:nl))
    allocate(room_g(ifrac,tile)%nodetemp(ufull_g(tile),1))
    allocate(road_g(ifrac,tile)%storage(ufull_g(tile),nl),roof_g(ifrac,tile)%storage(ufull_g(tile),nl))
    allocate(walle_g(ifrac,tile)%storage(ufull_g(tile),nl),wallw_g(ifrac,tile)%storage(ufull_g(tile),nl))
    allocate(slab_g(ifrac,tile)%storage(ufull_g(tile),nl),intm_g(ifrac,tile)%storage(ufull_g(tile),nl))
    allocate(room_g(ifrac,tile)%storage(ufull_g(tile),1)) 
  
    allocate(rfhyd_g(ifrac,tile)%surfwater(ufull_g(tile)),rfhyd_g(ifrac,tile)%snow(ufull_g(tile)))
    allocate(rfhyd_g(ifrac,tile)%den(ufull_g(tile)),rfhyd_g(ifrac,tile)%snowalpha(ufull_g(tile)))
    allocate(rfhyd_g(ifrac,tile)%leafwater(ufull_g(tile)),rfhyd_g(ifrac,tile)%soilwater(ufull_g(tile)))
    allocate(rdhyd_g(ifrac,tile)%surfwater(ufull_g(tile)),rdhyd_g(ifrac,tile)%snow(ufull_g(tile)))
    allocate(rdhyd_g(ifrac,tile)%den(ufull_g(tile)),rdhyd_g(ifrac,tile)%snowalpha(ufull_g(tile)))
    allocate(rdhyd_g(ifrac,tile)%leafwater(ufull_g(tile)),rdhyd_g(ifrac,tile)%soilwater(ufull_g(tile)))

    allocate(p_g(ifrac,tile)%lzom(ufull_g(tile)),p_g(ifrac,tile)%lzoh(ufull_g(tile)),p_g(ifrac,tile)%cndzmin(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%cduv(ufull_g(tile)),p_g(ifrac,tile)%cdtq(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%tscrn(ufull_g(tile)),p_g(ifrac,tile)%qscrn(ufull_g(tile)),p_g(ifrac,tile)%uscrn(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%u10(ufull_g(tile)),p_g(ifrac,tile)%emiss(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%bldheat(ufull_g(tile)),p_g(ifrac,tile)%bldcool(ufull_g(tile)),p_g(ifrac,tile)%traf(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surfrunoff(ufull_g(tile)),p_g(ifrac,tile)%irrig(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%transveg(ufull_g(tile)), p_g(ifrac,tile)%acond_vegw(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%soilmoist(ufull_g(tile)),p_g(ifrac,tile)%delsoilmoist(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%rootmoistc(ufull_g(tile)), p_g(ifrac,tile)%ulai(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%intgains_full(ufull_g(tile)),p_g(ifrac,tile)%storage_flux(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surferr(ufull_g(tile)),p_g(ifrac,tile)%atmoserr(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surferr_bias(ufull_g(tile)),p_g(ifrac,tile)%atmoserr_bias(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%soilwetness(ufull_g(tile)),p_g(ifrac,tile)%soilwater(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%snowmelt(ufull_g(tile)),p_g(ifrac,tile)%frac_sigma(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%delintercept(ufull_g(tile)),p_g(ifrac,tile)%snowt(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%vegt(ufull_g(tile)),p_g(ifrac,tile)%swe(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%surfstor(ufull_g(tile)),p_g(ifrac,tile)%snowfrac(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%salbedo(ufull_g(tile)),p_g(ifrac,tile)%calbedo(ufull_g(tile)))
    allocate(p_g(ifrac,tile)%taircanyon(ufull_g(tile)),p_g(ifrac,tile)%delswe(ufull_g(tile)))

    if ( ufull_g(tile)>0 ) then

      ! Initialise state variables
      roof_g(ifrac,tile)%nodetemp=1.  ! + urbtemp
      roof_g(ifrac,tile)%storage =0._8
      road_g(ifrac,tile)%nodetemp=1.  ! + urbtemp
      road_g(ifrac,tile)%storage =0._8
      walle_g(ifrac,tile)%nodetemp=1. ! + urbtemp
      walle_g(ifrac,tile)%storage=0._8
      wallw_g(ifrac,tile)%nodetemp=1. ! + urbtemp
      wallw_g(ifrac,tile)%storage=0._8
      slab_g(ifrac,tile)%nodetemp=1. ! + urbtemp
      slab_g(ifrac,tile)%storage=0._8
      intm_g(ifrac,tile)%nodetemp=1. ! + urbtemp
      intm_g(ifrac,tile)%storage=0._8
      room_g(ifrac,tile)%nodetemp=1.  ! + urbtemp
      room_g(ifrac,tile)%storage=0._8

      rfhyd_g(ifrac,tile)%surfwater=0.
      rfhyd_g(ifrac,tile)%snow=0.
      rfhyd_g(ifrac,tile)%den=minsnowden
      rfhyd_g(ifrac,tile)%snowalpha=maxsnowalpha
      rfhyd_g(ifrac,tile)%leafwater=0.
      rfhyd_g(ifrac,tile)%soilwater=0.
      rdhyd_g(ifrac,tile)%surfwater=0.
      rdhyd_g(ifrac,tile)%snow=0.
      rdhyd_g(ifrac,tile)%den=minsnowden
      rdhyd_g(ifrac,tile)%snowalpha=maxsnowalpha
      rdhyd_g(ifrac,tile)%leafwater=0.
      rdhyd_g(ifrac,tile)%soilwater=0.25
    
      p_g(ifrac,tile)%cndzmin=max(10.,0.1*f_g(tile)%bldheight+2.)                 ! updated in atebcalc
      p_g(ifrac,tile)%lzom=log(p_g(ifrac,tile)%cndzmin/(0.1*f_g(tile)%bldheight)) ! updated in atebcalc
      p_g(ifrac,tile)%lzoh=6.+p_g(ifrac,tile)%lzom ! (Kanda et al 2005)           ! updated in atebcalc
      p_g(ifrac,tile)%cduv=(vkar/p_g(ifrac,tile)%lzom)**2                         ! updated in atebcalc
      p_g(ifrac,tile)%cdtq=vkar**2/(p_g(ifrac,tile)%lzom*p_g(ifrac,tile)%lzoh)    ! updated in atebcalc
      p_g(ifrac,tile)%tscrn=1.      ! + urbtemp                                   ! updated in atebcalc
      p_g(ifrac,tile)%qscrn=0.                                                    ! updated in atebcalc
      p_g(ifrac,tile)%uscrn=0.                                                    ! updated in atebcalc
      p_g(ifrac,tile)%u10=0.                                                      ! updated in atebcalc
      p_g(ifrac,tile)%emiss=0.97                                                  ! updated in atebcalc
      p_g(ifrac,tile)%bldheat=0.
      p_g(ifrac,tile)%bldcool=0.
      p_g(ifrac,tile)%surfrunoff=0.
      p_g(ifrac,tile)%soilwetness=0.
      p_g(ifrac,tile)%soilwater=0.
      p_g(ifrac,tile)%transveg=0.
      p_g(ifrac,tile)%soilmoist=rdhyd_g(ifrac,tile)%soilwater*waterden*4.
      p_g(ifrac,tile)%rootmoistc=rdhyd_g(ifrac,tile)%soilwater*waterden*4.
      p_g(ifrac,tile)%delsoilmoist=0.
      p_g(ifrac,tile)%irrig=0.
      p_g(ifrac,tile)%acond_vegw=0.1
      p_g(ifrac,tile)%traf=0.
      p_g(ifrac,tile)%intgains_full=0.
      if ( ifrac==1 ) then 
        p_g(ifrac,tile)%frac_sigma=1.
      else
        p_g(ifrac,tile)%frac_sigma=0.  
      end if    
      p_g(ifrac,tile)%surferr=0._8
      p_g(ifrac,tile)%atmoserr=0._8
      p_g(ifrac,tile)%surferr_bias=0._8
      p_g(ifrac,tile)%atmoserr_bias=0._8
      p_g(ifrac,tile)%storage_flux=0._8
      p_g(ifrac,tile)%first_call=.true.
      p_g(ifrac,tile)%delintercept=0.
      p_g(ifrac,tile)%snowt=273.
      p_g(ifrac,tile)%vegt=273.
      p_g(ifrac,tile)%swe=0.
      p_g(ifrac,tile)%surfstor=0.
      p_g(ifrac,tile)%snowfrac=0.
      p_g(ifrac,tile)%salbedo=0.
      p_g(ifrac,tile)%calbedo=0.
      p_g(ifrac,tile)%taircanyon=273.
      p_g(ifrac,tile)%delswe=0.
      p_g(ifrac,tile)%ulai=1.

    end if
    
  end do  
  
end do


! for getqsat
table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
table(219:220)=(/ 29845.0, 31169.0 /)

return
end subroutine atebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine atebend(diag)

implicit none

integer, intent(in) :: diag
integer tile, ifrac

if (diag>=1) write(6,*) "Deallocating aTEB arrays"

if ( ateb_active ) then

  do tile = 1,ntiles

    deallocate(f_roof(tile)%depth,f_wall(tile)%depth,f_road(tile)%depth,f_slab(tile)%depth,f_intm(tile)%depth)
    deallocate(f_roof(tile)%volcp,f_wall(tile)%volcp,f_road(tile)%volcp,f_slab(tile)%volcp,f_intm(tile)%volcp)
    deallocate(f_roof(tile)%lambda,f_wall(tile)%lambda,f_road(tile)%lambda,f_slab(tile)%lambda,f_intm(tile)%lambda)
    deallocate(f_roof(tile)%alpha,f_wall(tile)%alpha,f_road(tile)%alpha)
    deallocate(f_roof(tile)%emiss,f_wall(tile)%emiss,f_road(tile)%emiss)
    deallocate(f_slab(tile)%emiss)

    deallocate(cnveg_g(tile)%sigma,cnveg_g(tile)%alpha)
    deallocate(cnveg_g(tile)%emiss,rfveg_g(tile)%sigma,rfveg_g(tile)%alpha,rfveg_g(tile)%emiss)
    deallocate(cnveg_g(tile)%zo,cnveg_g(tile)%lai,cnveg_g(tile)%rsmin,rfveg_g(tile)%zo,rfveg_g(tile)%lai)
    deallocate(rfveg_g(tile)%rsmin,cnveg_g(tile)%temp,rfveg_g(tile)%temp)

    deallocate(intl_g(tile)%viewf,intl_g(tile)%psi)

    deallocate(f_g(tile)%sigmabld,f_g(tile)%hwratio,f_g(tile)%bldheight,f_g(tile)%coeffbldheight)
    deallocate(f_g(tile)%effhwratio)
    deallocate(f_g(tile)%industryfg,f_g(tile)%intgains_flr,f_g(tile)%trafficfg,f_g(tile)%vangle)
    deallocate(f_g(tile)%ctime,f_g(tile)%hangle,f_g(tile)%fbeam,f_g(tile)%weekdayload)
    deallocate(f_g(tile)%swilt,f_g(tile)%sfc,f_g(tile)%ssat)
    deallocate(f_g(tile)%bldairtemp,f_g(tile)%rfvegdepth)
    deallocate(f_g(tile)%intmassn,f_g(tile)%infilach,f_g(tile)%ventilach,f_g(tile)%heatprop,f_g(tile)%coolprop)
    deallocate(f_g(tile)%sigmau)
    
    do ifrac = 1,nfrac  
      
      deallocate(rfhyd_g(ifrac,tile)%surfwater,rfhyd_g(ifrac,tile)%snow,rfhyd_g(ifrac,tile)%den,rfhyd_g(ifrac,tile)%snowalpha)
      deallocate(rdhyd_g(ifrac,tile)%surfwater,rdhyd_g(ifrac,tile)%snow,rdhyd_g(ifrac,tile)%den,rdhyd_g(ifrac,tile)%snowalpha)
      deallocate(rdhyd_g(ifrac,tile)%leafwater,rdhyd_g(ifrac,tile)%soilwater,rfhyd_g(ifrac,tile)%leafwater)
      deallocate(rfhyd_g(ifrac,tile)%soilwater)

      deallocate(roof_g(ifrac,tile)%nodetemp,road_g(ifrac,tile)%nodetemp,walle_g(ifrac,tile)%nodetemp)
      deallocate(wallw_g(ifrac,tile)%nodetemp)
      deallocate(slab_g(ifrac,tile)%nodetemp,intm_g(ifrac,tile)%nodetemp,room_g(ifrac,tile)%nodetemp)
      deallocate(road_g(ifrac,tile)%storage,roof_g(ifrac,tile)%storage,walle_g(ifrac,tile)%storage,wallw_g(ifrac,tile)%storage)
      deallocate(slab_g(ifrac,tile)%storage,intm_g(ifrac,tile)%storage,room_g(ifrac,tile)%storage)
      
      deallocate(p_g(ifrac,tile)%lzom,p_g(ifrac,tile)%lzoh,p_g(ifrac,tile)%cndzmin,p_g(ifrac,tile)%cduv,p_g(ifrac,tile)%cdtq)
      deallocate(p_g(ifrac,tile)%tscrn,p_g(ifrac,tile)%qscrn,p_g(ifrac,tile)%uscrn,p_g(ifrac,tile)%u10,p_g(ifrac,tile)%emiss)
      deallocate(p_g(ifrac,tile)%surferr,p_g(ifrac,tile)%atmoserr,p_g(ifrac,tile)%surferr_bias,p_g(ifrac,tile)%atmoserr_bias)
      deallocate(p_g(ifrac,tile)%bldheat,p_g(ifrac,tile)%bldcool,p_g(ifrac,tile)%traf,p_g(ifrac,tile)%intgains_full)
      deallocate(p_g(ifrac,tile)%surfrunoff,p_g(ifrac,tile)%irrig,p_g(ifrac,tile)%acond_vegw)
      deallocate(p_g(ifrac,tile)%soilwetness,p_g(ifrac,tile)%soilwater)
      deallocate(p_g(ifrac,tile)%transveg,p_g(ifrac,tile)%soilmoist,p_g(ifrac,tile)%delsoilmoist)
      deallocate(p_g(ifrac,tile)%rootmoistc,p_g(ifrac,tile)%ulai)
      deallocate(p_g(ifrac,tile)%storage_flux,p_g(ifrac,tile)%snowmelt)
      deallocate(p_g(ifrac,tile)%delintercept,p_g(ifrac,tile)%snowt)
      deallocate(p_g(ifrac,tile)%vegt,p_g(ifrac,tile)%swe)
      deallocate(p_g(ifrac,tile)%surfstor,p_g(ifrac,tile)%snowfrac)
      deallocate(p_g(ifrac,tile)%salbedo,p_g(ifrac,tile)%calbedo)
      deallocate(p_g(ifrac,tile)%taircanyon,p_g(ifrac,tile)%delswe)
      
    end do  

  end do
  
  deallocate( roof_g, road_g, walle_g, wallw_g, slab_g, intm_g )
  deallocate( room_g )
  deallocate( f_roof, f_road, f_wall, f_slab, f_intm )
  deallocate( rfhyd_g, rdhyd_g )
  deallocate( cnveg_g, rfveg_g )
  deallocate( intl_g )
  deallocate( f_g )
  deallocate( p_g )
  deallocate( ufull_g )
  deallocate( upack_g )

end if
    
return
end subroutine atebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general version of tebload

subroutine atebloadd(urban,mode,ifrac,diag)

implicit none

integer, intent(in) :: ifrac, diag
integer ii, tile, is, ie
real, dimension(ifull), intent(in) :: urban
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if (diag>=1) write(6,*) "Load aTEB state array"
if (.not.ateb_active) return

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        roof_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( roof_g(ifrac,tile)%nodetemp(:,ii)>150. )
          roof_g(ifrac,tile)%nodetemp(:,ii) = roof_g(ifrac,tile)%nodetemp(:,ii) - urbtemp
        end where
      end if
    end do
    return
  end if
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        walle_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( walle_g(ifrac,tile)%nodetemp(:,ii)>150. )
          walle_g(ifrac,tile)%nodetemp(:,ii) = walle_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        wallw_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( wallw_g(ifrac,tile)%nodetemp(:,ii)>150. )
          wallw_g(ifrac,tile)%nodetemp(:,ii) = wallw_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        road_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( road_g(ifrac,tile)%nodetemp(:,ii)>150. )
          road_g(ifrac,tile)%nodetemp(:,ii) = road_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        slab_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( slab_g(ifrac,tile)%nodetemp(:,ii)>150. )
          slab_g(ifrac,tile)%nodetemp(:,ii) = slab_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if  
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        intm_g(ifrac,tile)%nodetemp(:,ii)=pack(urban(is:ie),upack_g(:,tile))
        where ( intm_g(ifrac,tile)%nodetemp(:,ii)>150. )
          intm_g(ifrac,tile)%nodetemp(:,ii) = intm_g(ifrac,tile)%nodetemp(:,ii) - urbtemp  
        end where  
      end if
    end do
    return
  end if   
end do  
  
select case(mode)
  case("roomtemp")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        room_g(ifrac,tile)%nodetemp(:,1)=pack(urban(is:ie),upack_g(:,tile))  
        where ( room_g(ifrac,tile)%nodetemp(:,1)>150. )
          room_g(ifrac,tile)%nodetemp(:,1) = room_g(ifrac,tile)%nodetemp(:,1) - urbtemp  
        end where 
      end if
    end do
    return
  case("canyonsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%soilwater=pack(urban(is:ie),upack_g(:,tile))
      end if
    end do
    return
  case("roofsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%soilwater=pack(urban(is:ie),upack_g(:,tile))
      end if
    end do
    return
  case("roadsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%surfwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%surfwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("canyonleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%leafwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%leafwater=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roadsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%snow=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%snow=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roadsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%den=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%den=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roadsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rdhyd_g(ifrac,tile)%snowalpha=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
  case("roofsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        rfhyd_g(ifrac,tile)%snowalpha=pack(urban(is:ie),upack_g(:,tile))  
      end if
    end do
    return
end select
  
write(6,*) "ERROR: Unknown mode for atebloadd ",trim(mode)
stop
    
return
end subroutine atebloadd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine atebtype_standard(itype,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
integer, dimension(ifull), intent(in) :: itype

if (diag>=1) write(6,*) "Load aTEB building properties"
if (.not.ateb_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebtype_thread(itype(is:ie),diag,f_g(tile),cnveg_g(tile),rfveg_g(tile), &
                         f_roof(tile),f_road(tile),f_wall(tile),f_slab(tile),     &
                         f_intm(tile),intl_g(tile),upack_g(:,tile),               &
                         ufull_g(tile))
  end if
end do

return
end subroutine atebtype_standard

subroutine atebtype_thread(itype,diag,fp,cnveg,rfveg,fp_roof,fp_road,fp_wall,fp_slab, &
                           fp_intm,intl,upack,ufull)

implicit none

integer, intent(in) :: diag, ufull
integer ii,j,ierr,nlp
integer, dimension(imax), intent(in) :: itype
integer, dimension(ufull) :: itmp
integer, parameter :: maxtype = 8
integer conductmeth, useonewall ! depreciated
real x
real, dimension(ufull) :: tsigveg,tsigmabld
! In-canyon vegetation fraction
real, dimension(maxtype) ::    csigvegc=(/ 0.38, 0.45, 0.38, 0.34, 0.05, 0.40, 0.30, 0.20 /)
! Green roof vegetation fraction
real, dimension(maxtype) ::    csigvegr=(/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
! Area fraction occupied by buildings
real, dimension(maxtype) ::   csigmabld=(/ 0.45, 0.40, 0.45, 0.46, 0.65, 0.40, 0.45, 0.50 /)
! Building height (m)
real, dimension(maxtype) ::  cbldheight=(/   6.,   4.,   6.,   8.,  18.,   4.,   8.,  12. /)
! Building height to width ratio
real, dimension(maxtype) ::    chwratio=(/  0.4,  0.2,  0.4,  0.6,   2.,  0.5,   1.,  1.5 /)
! Industrial sensible heat flux (W m^-2)
real, dimension(maxtype) :: cindustryfg=(/   0.,   0.,   0.,   0.,   0.,  10.,  20.,  30. /)
! Internal gains sensible heat flux [floor] (W m^-2)
real, dimension(maxtype) ::   cintgains=(/   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype) ::  ctrafficfg=(/  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5 /)
! Comfort temperature (K)
real, dimension(maxtype) :: cbldtemp=(/ 295.16, 295.16, 295.16, 295.16, 295.16, 295.16, 295.16, 295.16 /)
! Roof albedo
real, dimension(maxtype) ::  croofalpha=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)    ! (Fortuniak 08) Masson = 0.15
! Wall albedo
real, dimension(maxtype) ::  cwallalpha=(/ 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 /)    ! (Fortuniak 08) Masson = 0.25
! Road albedo
real, dimension(maxtype) ::  croadalpha=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)    ! (Fortuniak 08) Masson = 0.08
! Canyon veg albedo
real, dimension(maxtype) ::  cvegalphac=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof veg albedo
real, dimension(maxtype) ::  cvegalphar=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof emissivity
real, dimension(maxtype) ::  croofemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
! Wall emissivity
real, dimension(maxtype) ::  cwallemiss=(/ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 /) 
! Road emissivity
real, dimension(maxtype) ::  croademiss=(/ 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94 /)
! Slab emissivity
real, dimension(maxtype) ::  cslabemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /) 
! Canyon veg emissivity
real, dimension(maxtype) ::  cvegemissc=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Roof veg emissivity
real, dimension(maxtype) ::  cvegemissr=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Green roof soil depth
real, dimension(maxtype) ::   cvegdeptr=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)
! Roughness length of in-canyon vegetation (m)
real, dimension(maxtype) ::    czovegc=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
! In-canyon vegetation LAI
real, dimension(maxtype) ::  cvegrlaic=(/   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0 /)
! Unconstrained canopy stomatal resistance
real, dimension(maxtype) :: cvegrsminc=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
! Roughness length of green roof vegetation (m)
real, dimension(maxtype) ::    czovegr=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
! Green roof vegetation LAI
real, dimension(maxtype) ::  cvegrlair=(/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0 /)
! Unconstrained canopy stomatal resistance
real, dimension(maxtype) :: cvegrsminr=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
! Soil wilting point (m^3 m^-3)
real, dimension(maxtype) ::     cswilt=(/  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18 /)
! Soil field capacity (m^3 m^-3)
real, dimension(maxtype) ::       csfc=(/  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26 /)
! Soil saturation point (m^3 m^-3)
real, dimension(maxtype) ::      cssat=(/  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42 /)
! Infiltration air volume changes per hour (m^3 m^-3)
real, dimension(maxtype) ::  cinfilach=(/  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50 /)
! Ventilation air volume changes per hour (m^3 m^-3)
real, dimension(maxtype) :: cventilach=(/  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /)
! Fraction of spaces with heating devices
real, dimension(maxtype) ::  cheatprop=(/  0.5,  0.5,  0.5,  0.5,  1.0,  0.00,  0.00,  0.00 /)
! Fraction of spaces with cooling devices
real, dimension(maxtype) ::  ccoolprop=(/  0.5,  0.5,  0.5,  0.5,  1.0,  0.00,  0.00,  0.00 /)

real, dimension(maxtype,nl) :: croofdepth
real, dimension(maxtype,nl) :: cwalldepth
real, dimension(maxtype,nl) :: croaddepth
real, dimension(maxtype,nl) :: cslabdepth
real, dimension(maxtype,nl) :: croofcp
real, dimension(maxtype,nl) :: cwallcp
real, dimension(maxtype,nl) :: croadcp
real, dimension(maxtype,nl) :: cslabcp
real, dimension(maxtype,nl) :: crooflambda
real, dimension(maxtype,nl) :: cwalllambda
real, dimension(maxtype,nl) :: croadlambda
real, dimension(maxtype,nl) :: cslablambda

logical, dimension(imax), intent(in) :: upack

type(fparmdata), intent(inout) :: fp
type(vegdata), intent(inout) :: cnveg, rfveg
type(facetparams), intent(inout) :: fp_roof, fp_road, fp_wall, fp_slab, fp_intm
type(intldata), intent(inout) :: intl

namelist /atebnml/  resmeth,useonewall,zohmeth,acmeth,intairtmeth,intmassmeth,nrefl,vegmode,soilunder, &
                    conductmeth,cvcoeffmeth,statsmeth,scrnmeth,wbrelaxc,wbrelaxr,lweff,iqt,            &
                    infilmeth,lwintmeth,infilalpha
namelist /atebsnow/ zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden
namelist /atebgen/  refheight,zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf, &
                    acfactor,ac_heatcap,ac_coolcap,ac_deltat,cyc_prop,cyc_base,cyc_traf,cyc_tran
namelist /atebtile/ czovegc,cvegrlaic,cvegrsminc,czovegr,cvegrlair,cvegrsminr,cswilt,csfc,cssat,       &
                    cvegemissc,cvegemissr,cvegdeptr,cvegalphac,cvegalphar,csigvegc,csigvegr,           &
                    csigmabld,cbldheight,chwratio,cindustryfg,cintgains,ctrafficfg,cbldtemp,           &
                    croofalpha,cwallalpha,croadalpha,croofemiss,cwallemiss,croademiss,croofdepth,      &
                    cwalldepth,croaddepth,croofcp,cwallcp,croadcp,crooflambda,cwalllambda,croadlambda, &
                    cslabdepth,cslabcp,cslablambda,cinfilach,cventilach,cheatprop,ccoolprop

! facet array where: rows=maxtypes (landtypes) and columns=nl (material layers)
nlp=nl/4 ! number of layers in each material segment (over 4 material segments)
! depths (m)
croofdepth= reshape((/ ((0.01/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.09/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.40/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.10/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
cwalldepth= reshape((/ ((0.01/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.04/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.10/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
croaddepth= reshape((/ ((0.01/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.04/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.45/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((3.50/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
cslabdepth=reshape((/  ((0.05/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
! heat capacity (J m^-3 K^-1)
croofcp =   reshape((/ ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((0.28E6, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                       ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
cwallcp =   reshape((/ ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
croadcp =   reshape((/ ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((1.28E6, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                       ((1.28E6, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                       (/maxtype,nl/))
cslabcp=   reshape((/  ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                       (/maxtype,nl/))
! heat conductivity (W m^-1 K^-1)
crooflambda=reshape((/ ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((0.0800, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                       ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
cwalllambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
croadlambda=reshape((/ ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((0.2513, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                       ((0.2513, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                       (/maxtype,nl/))
cslablambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                       (/maxtype,nl/))

itmp=pack(itype,upack)
where ( itmp>maxtype ) ! set default to new urban types as generic
  itmp = 1
end where
if ((minval(itmp)<1).or.(maxval(itmp)>maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

if (atebnmlfile/=0) then
  open(unit=atebnmlfile,file='ateb.nml',action="read",iostat=ierr)
  if (ierr==0) then
    write(6,*) "Reading ateb.nml"
    read(atebnmlfile,nml=atebnml)
    read(atebnmlfile,nml=atebsnow)
    read(atebnmlfile,nml=atebgen)
    read(atebnmlfile,nml=atebtile)
    close(atebnmlfile)
  end if
end if

select case(vegmode)
  case(0)
    tsigveg=0.5*csigvegc(itmp)/(1.-0.5*csigvegc(itmp))
    tsigmabld=csigmabld(itmp)/(1.-0.5*csigvegc(itmp))
    fp%sigmau=fp%sigmau*(1.-0.5*csigvegc(itmp))
  case(1)
    tsigveg=0.
    tsigmabld=csigmabld(itmp)/(1.-csigvegc(itmp))
    fp%sigmau=fp%sigmau*(1.-csigvegc(itmp))
  case(2)
    tsigveg=csigvegc(itmp)
    tsigmabld=csigmabld(itmp)
  case DEFAULT
    if (vegmode<0) then
      x=real(abs(vegmode))/100.
      x=max(min(x,1.),0.)
      tsigveg=x*csigvegc(itmp)/(1.-(1.-x)*csigvegc(itmp))
      tsigmabld=csigmabld(itmp)/(1.-(1.-x)*csigvegc(itmp))
      fp%sigmau=fp%sigmau*(1.-(1.-x)*csigvegc(itmp))
    else
      write(6,*) "ERROR: Unsupported vegmode ",vegmode
      stop
    end if
end select
cnveg%sigma=max(min(tsigveg/(1.-tsigmabld),1.),0.)
rfveg%sigma=max(min(csigvegr(itmp),1.),0.)
fp%sigmabld=max(min(tsigmabld,1.),0.)
!fp%hwratio=chwratio(itmp)*fp%sigmabld/(1.-fp%sigmabld) ! MJT suggested new definition
fp%hwratio=chwratio(itmp)          ! MJL simple definition

fp%industryfg=cindustryfg(itmp)
fp%intgains_flr=cintgains(itmp)
fp%trafficfg=ctrafficfg(itmp)
fp%bldheight=cbldheight(itmp)
fp_roof%alpha=croofalpha(itmp)
fp_wall%alpha=cwallalpha(itmp)
fp_road%alpha=croadalpha(itmp)
cnveg%alpha=cvegalphac(itmp)
rfveg%alpha=cvegalphar(itmp)
fp_roof%emiss=croofemiss(itmp)
fp_wall%emiss=cwallemiss(itmp)
fp_road%emiss=croademiss(itmp)
cnveg%emiss=cvegemissc(itmp)
rfveg%emiss=cvegemissr(itmp)
fp%bldairtemp=cbldtemp(itmp) - urbtemp
fp%rfvegdepth=cvegdeptr(itmp)
do ii=1,nl
  fp_roof%depth(:,ii)=croofdepth(itmp,ii)
  fp_wall%depth(:,ii)=cwalldepth(itmp,ii)
  fp_road%depth(:,ii)=croaddepth(itmp,ii)
  fp_roof%lambda(:,ii)=crooflambda(itmp,ii)
  fp_wall%lambda(:,ii)=cwalllambda(itmp,ii)
  fp_road%lambda(:,ii)=croadlambda(itmp,ii)
  fp_roof%volcp(:,ii)=croofcp(itmp,ii)
  fp_wall%volcp(:,ii)=cwallcp(itmp,ii)
  select case(soilunder)
    case(0) ! storage under road only
      fp_road%volcp(:,ii)=croadcp(itmp,ii)
    case(1) ! storage under road and canveg
      fp_road%volcp(:,ii)=croadcp(itmp,ii)/(1.-cnveg%sigma)
    case(2) ! storage under road and bld
      fp_road%volcp(:,ii)=croadcp(itmp,ii)*(1./(1.-cnveg%sigma)*(1./(1.-fp%sigmabld)-1.) +1.)
    case(3) ! storage under road and canveg and bld (100% of grid point)
      fp_road%volcp(:,ii)=croadcp(itmp,ii)/(1.-cnveg%sigma)/(1.-fp%sigmabld)
    case DEFAULT
      write(6,*) "ERROR: Unknown soilunder mode ",soilunder
      stop
  end select
end do
cnveg%zo=czovegc(itmp)
cnveg%lai=cvegrlaic(itmp)
cnveg%rsmin=cvegrsminc(itmp)/max(cnveg%lai,1.E-8)
rfveg%zo=czovegr(itmp)
rfveg%lai=cvegrlair(itmp)
rfveg%rsmin=cvegrsminr(itmp)/max(rfveg%lai,1.E-8)
fp%swilt=cswilt(itmp)
fp%sfc=csfc(itmp)
fp%ssat=cssat(itmp)

! for varying internal temperature
fp_slab%emiss = cslabemiss(itmp)
fp%infilach = cinfilach(itmp)
fp%ventilach = cventilach(itmp)
fp%heatprop = cheatprop(itmp)
fp%coolprop = ccoolprop(itmp)
do ii=1,nl
  fp_slab%depth(:,ii)=cslabdepth(itmp,ii)
  fp_intm%depth(:,ii)=cslabdepth(itmp,ii)  ! internal mass material same as slab
  fp_slab%lambda(:,ii)=cslablambda(itmp,ii)
  fp_intm%lambda(:,ii)=cslablambda(itmp,ii)
  fp_slab%volcp(:,ii)=cslabcp(itmp,ii)
  fp_intm%volcp(:,ii)=cslabcp(itmp,ii)
end do

! Here we modify the effective canyon geometry to account for in-canyon vegetation tall vegetation
fp%coeffbldheight = max(fp%bldheight-6.*cnveg%zo,0.2)/fp%bldheight
fp%effhwratio   = fp%hwratio*fp%coeffbldheight

call init_internal(fp)
call init_lwcoeff(fp,intl,ufull)

if ( diag>0 ) then
  write(6,*) 'hwratio, eff',fp%hwratio, fp%effhwratio
  write(6,*) 'bldheight, eff',fp%bldheight, fp%coeffbldheight
  write(6,*) 'sigmabld, sigmavegc', fp%sigmabld, cnveg%sigma
  write(6,*) 'roadcp multiple for soilunder:', soilunder,fp_road%volcp(itmp,1)/croadcp(itmp,1)
end if

return
end subroutine atebtype_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine atebdeftype(paramdata,typedata,paramname,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie, i, ifrac, maxtype
integer, dimension(ifull), intent(in) :: typedata
integer, dimension(imax) :: itmp
real, dimension(:), intent(in) :: paramdata
logical found
character(len=*), intent(in) :: paramname
character(len=20) :: vname

if ( diag>=1 ) write(6,*) "Load aTEB parameter: ",trim(paramname)
if ( .not.ateb_active ) return

maxtype = size(paramdata)

select case(paramname)
  case('bldheight')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%bldheight = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('hwratio')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%hwratio = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do
  case('sigvegc')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%sigma = paramdata(itmp(1:ufull_g(tile)))/(1.-f_g(tile)%sigmabld)  
      end if
    end do  
  case('sigmabld')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%sigma = cnveg_g(tile)%sigma*(1.-f_g(tile)%sigmabld)/(1.-paramdata(itmp(1:ufull_g(tile))))
        f_g(tile)%sigmabld = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('industryfg')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%industryfg = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('trafficfg')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%trafficfg = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roofalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_roof(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('wallalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_wall(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roadalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_road(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roofemiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_roof(tile)%emiss = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('wallemiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_wall(tile)%emiss = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('roademiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_road(tile)%emiss = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do   
  case('vegalphac')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%alpha = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('zovegc')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        cnveg_g(tile)%zo = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('infilach')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%infilach = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('intgains')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%intgains_flr = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case('bldairtemp')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%bldairtemp = paramdata(itmp(1:ufull_g(tile))) - urbtemp
        do ifrac = 1,nfrac
          room_g(ifrac,tile)%nodetemp(:,1) = f_g(tile)%bldairtemp
        end do  
      end if
    end do   
  case('heatprop')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%heatprop = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do   
  case('coolprop')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
        if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
          write(6,*) "ERROR: Urban type is out of range"
          stop 
        end if
        f_g(tile)%coolprop = paramdata(itmp(1:ufull_g(tile)))
      end if
    end do  
  case default
    found = .false.
    do i = 1,4
      write(vname,'("roofthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_roof(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roofcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_roof(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roofcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_roof(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_wall(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_wall(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_wall(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_road(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_road(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_road(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_slab(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
            f_intm(tile)%depth(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_slab(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
            f_intm(tile)%volcp(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            itmp(1:ufull_g(tile)) = pack(typedata(is:ie),upack_g(:,tile))
            if ( minval(itmp(1:ufull_g(tile)))<1 .or. maxval(itmp(1:ufull_g(tile)))>maxtype ) then
              write(6,*) "ERROR: Urban type is out of range"
              stop 
            end if
            f_slab(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
            f_intm(tile)%lambda(:,i) = paramdata(itmp(1:ufull_g(tile)))
          end if
        end do 
        found = .true.
        exit
      end if
    end do
    if ( .not.found ) then
      write(6,*) "ERROR: Unknown aTEB parameter name ",trim(paramname)
      stop
    end if  
end select

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    ! Here we modify the effective canyon geometry to account for in-canyon vegetation tall vegetation
    f_g(tile)%coeffbldheight = max(f_g(tile)%bldheight-6.*cnveg_g(tile)%zo,0.2)/f_g(tile)%bldheight
    f_g(tile)%effhwratio     = f_g(tile)%hwratio*f_g(tile)%coeffbldheight   
    
    call init_internal(f_g(tile))
    call init_lwcoeff(f_g(tile),intl_g(tile),ufull_g(tile))
  end if
end do

return
end subroutine atebdeftype

subroutine atebdeftype_export(paramdata,paramname,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie, i
real, dimension(ifull), intent(inout) :: paramdata
logical found
character(len=*), intent(in) :: paramname
character(len=20) :: vname

if ( diag>=1 ) write(6,*) "Export aTEB parameter: ",trim(paramname)
if ( .not.ateb_active ) return

select case(paramname)
  case('bldheight')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%bldheight,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('hwratio')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax 
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%hwratio,upack_g(:,tile),paramdata(is:ie))
      end if
    end do
  case('sigvegc')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(cnveg_g(tile)%sigma*(1.-f_g(tile)%sigmabld),upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('sigmabld')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax 
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%sigmabld,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('industryfg')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%industryfg,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('trafficfg')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%trafficfg,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('roofalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_roof(tile)%alpha,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('wallalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_wall(tile)%alpha,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('roadalpha')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax 
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_road(tile)%alpha,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('roofemiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax 
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_roof(tile)%emiss,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('wallemiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax 
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_wall(tile)%emiss,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('roademiss')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_road(tile)%emiss,upack_g(:,tile),paramdata(is:ie))
      end if
    end do   
  case('vegalphac')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax 
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(cnveg_g(tile)%alpha,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('zovegc')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(cnveg_g(tile)%zo,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('infilach')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%infilach,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case('intgains')
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax  
      if ( ufull_g(tile)>0 ) then
        paramdata(is:ie)=unpack(f_g(tile)%intgains_flr,upack_g(:,tile),paramdata(is:ie))
      end if
    end do  
  case default
    found = .false.
    do i = 1,4
      write(vname,'("roofthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_roof(tile)%depth(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roofcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax 
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_roof(tile)%volcp(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roofcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_roof(tile)%lambda(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax 
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_wall(tile)%depth(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax 
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_wall(tile)%volcp(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("wallcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_wall(tile)%lambda(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_road(tile)%depth(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_road(tile)%volcp(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("roadcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_road(tile)%lambda(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabthick",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_slab(tile)%depth(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabcp",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_slab(tile)%volcp(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
      write(vname,'("slabcond",(I1.1))') i
      if ( trim(paramname)==trim(vname) ) then
        do tile = 1,ntiles
          is = (tile-1)*imax + 1
          ie = tile*imax  
          if ( ufull_g(tile)>0 ) then
            paramdata(is:ie)=unpack(f_slab(tile)%lambda(:,i),upack_g(:,tile),paramdata(is:ie))
          end if
        end do 
        found = .true.
        exit
      end if
    end do
    if ( .not.found ) then
      write(6,*) "ERROR: Unknown aTEB parameter name ",trim(paramname)
      stop
    end if  
end select

return
end subroutine atebdeftype_export

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_lwcoeff(fp,intl,ufull)
! This subroutine calculates longwave reflection coefficients (int_psi) at each surface
! longwave coefficients do not change, so this subroutine should only be run once
! Infinite reflections per Harman et al., (2004) "Radiative Exchange in Urban Street Canyons"
! Per method in "Radiation Heat Transfer, Sparrow & Cess 1978, Ch 3-3"
! array surface order is: (1) floor; (2) wallw; (3) ceiling; (4) walle

! local variables
integer, intent(in) :: ufull
real, dimension(ufull,4,4) :: chi
real, dimension(4,4)       :: krondelta
real, dimension(ufull)     :: h, w
real, dimension(ufull,4)   :: epsil   ! floor, wall, ceiling, wall emissivity array
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
!                    i=1,ufull_g)/), (/ufull_g,4/))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general version of atebsave

subroutine atebsaved(urban,mode,ifrac,diag,rawtemp)

implicit none

integer, intent(in) :: ifrac, diag
integer ii, tile, is, ie
real, dimension(ifull), intent(inout) :: urban
real urbtempadj
logical, intent(in), optional :: rawtemp
logical rawmode
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if (diag>=1) write(6,*) "Load aTEB state array"
if (.not.ateb_active) return

rawmode = .false.
if ( present(rawtemp) ) then
  rawmode = rawtemp
end if

if ( rawmode ) then
  urbtempadj = 0.
else
  urbtempadj = urbtemp
end if

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(roof_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(walle_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(wallw_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(road_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(slab_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if  
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(intm_g(ifrac,tile)%nodetemp(:,ii)+urbtempadj,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  end if   
end do  

select case(mode)
  case("roomtemp")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(room_g(ifrac,tile)%nodetemp(:,1)+urbtempadj,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("canyonsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%soilwater,upack_g(:,tile),urban(is:ie))  
      end if
    end do
    return
  case("roofsoilmoisture")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%soilwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%surfwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofsurfacewater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%surfwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("canyonleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%leafwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofleafwater")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%leafwater,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%snow,upack_g(:,tile),urban(is:ie))
      end if
    end do
    return
  case("roofsnowdepth")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%snow,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%den,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofsnowdensity")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%den,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roadsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rdhyd_g(ifrac,tile)%snowalpha,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
  case("roofsnowalbedo")
    do tile = 1,ntiles
      is = (tile-1)*imax + 1
      ie = tile*imax
      if ( ufull_g(tile)>0 ) then
        urban(is:ie)=unpack(rfhyd_g(ifrac,tile)%snowalpha,upack_g(:,tile),urban(is:ie))    
      end if
    end do
    return
end select

write(6,*) "ERROR: Unknown mode for atebsaved ",trim(mode)
stop

return
end subroutine atebsaved

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine collects and passes energy closure information to atebwrap

subroutine energyrecord(o_atmoserr,o_atmoserr_bias,o_surferr,o_surferr_bias, &
                        o_heating,o_cooling,o_intgains,o_traf,o_bldtemp)

implicit none

integer tile, is, ie, ifrac
real, dimension(ifull), intent(out) :: o_atmoserr,o_atmoserr_bias,o_surferr,o_surferr_bias
real, dimension(ifull), intent(out) :: o_heating,o_cooling,o_intgains,o_traf,o_bldtemp

o_atmoserr = 0.
o_atmoserr_bias = 0.
o_surferr = 0.
o_surferr_bias = 0.
o_heating = 0.
o_cooling = 0.
o_intgains = 0.
o_traf = 0.
o_bldtemp = 0.

if ( .not.ateb_active ) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
      
    do ifrac = 1,nfrac  
      
      p_g(ifrac,tile)%atmoserr_bias = p_g(ifrac,tile)%atmoserr_bias + p_g(ifrac,tile)%atmoserr
      p_g(ifrac,tile)%surferr_bias = p_g(ifrac,tile)%surferr_bias + p_g(ifrac,tile)%surferr

      o_atmoserr(is:ie) = o_atmoserr(is:ie)            &
          + unpack(real(p_g(ifrac,tile)%atmoserr)*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_surferr(is:ie) = o_surferr(is:ie)              &
          + unpack(real(p_g(ifrac,tile)%surferr)*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_atmoserr_bias(is:ie) = o_atmoserr_bias(is:ie)  &
          + unpack(real(p_g(ifrac,tile)%atmoserr_bias)*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_surferr_bias(is:ie) = o_surferr_bias(is:ie)    &
          + unpack(real(p_g(ifrac,tile)%surferr_bias)*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_heating(is:ie) = o_heating(is:ie)              &
          + unpack(p_g(ifrac,tile)%bldheat*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_cooling(is:ie) = o_cooling(is:ie)              &
          + unpack(p_g(ifrac,tile)%bldcool*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_intgains(is:ie) = o_intgains(is:ie)            &
          + unpack(p_g(ifrac,tile)%intgains_full*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_traf(is:ie) = o_traf(is:ie)                    &
          + unpack(p_g(ifrac,tile)%traf*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
      o_bldtemp(is:ie) = o_bldtemp(is:ie)              &
          + unpack((room_g(ifrac,tile)%nodetemp(:,1)+urbtemp)*p_g(ifrac,tile)%frac_sigma,upack_g(:,tile),0.)
    
    end do  
  end if
end do

return
end subroutine energyrecord

subroutine atebenergy_standard(o_data,mode,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: o_data
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Extract energy output"
if (.not.ateb_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebenergy_thread(o_data(is:ie),mode,diag,f_g(tile),p_g(:,tile),        & 
                           rdhyd_g(:,tile),rfhyd_g(:,tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine atebenergy_standard

subroutine atebenergy_thread(o_data,mode,diag,fp,pd,rdhyd,rfhyd,upack,ufull)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: o_data
real, dimension(ufull) :: ctmp, dtmp
character(len=*), intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd
type(hydrodata), dimension(nfrac), intent(in) :: rdhyd, rfhyd

if ( diag>=2 ) write(6,*) "THREAD: Extract energy output"
if ( ufull==0 ) return

select case(mode)
  case("anthropogenic")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + (pd(ifrac)%bldheat+pd(ifrac)%bldcool+pd(ifrac)%traf+fp%industryfg+pd(ifrac)%intgains_full) &
                    *pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("elecgas")
    dtmp = 0.
    do ifrac = 1,nfrac   
      dtmp = dtmp + real(pd(ifrac)%bldheat+pd(ifrac)%bldcool+pd(ifrac)%intgains_full)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("heating")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + real(pd(ifrac)%bldheat)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("cooling")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + real(pd(ifrac)%bldcool)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("storage")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + real(pd(ifrac)%storage_flux)*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case default
    write(6,*) "ERROR: Unknown atebenergy mode ",trim(mode)
    stop
end select    

return
end subroutine atebenergy_thread

subroutine atebmisc_standard(o_data,mode,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: o_data
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Extract energy output"
if (.not.ateb_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebmisc_thread(o_data(is:ie),mode,diag,room_g(:,tile),f_g(tile),p_g(:,tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine atebmisc_standard

subroutine atebmisc_thread(o_data,mode,diag,room,fp,pd,upack,ufull)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: o_data
real, dimension(ufull) :: ctmp, dtmp
character(len=*), intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(facetdata), dimension(nfrac), intent(in) :: room
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if ( diag>=2 ) write(6,*) "THREAD: Extract energy output"
if ( ufull==0 ) return

select case(mode)
  case("emissivity")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%emiss*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("snowt")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%snowt*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(o_data,upack)
    ctmp=(1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data=unpack(ctmp,upack,o_data)    
  case("vegt")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%vegt*pd(ifrac)%frac_sigma  
    end do    
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("tairbuilding")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + (room(ifrac)%nodetemp(:,1)+urbtemp)*pd(ifrac)%frac_sigma  
    end do    
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("salbedo")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%salbedo*pd(ifrac)%frac_sigma  
    end do    
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)      
  case("calbedo")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%calbedo*pd(ifrac)%frac_sigma  
    end do    
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("urbanlai")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%ulai*pd(ifrac)%frac_sigma  
    end do    
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case("taircanyon")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%taircanyon*pd(ifrac)%frac_sigma  
    end do    
    ctmp = pack(o_data, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    o_data = unpack(ctmp+urbtemp, upack, o_data)      
  case default
    write(6,*) "ERROR: Unknown atebmisc mode ",trim(mode)
    stop
end select    

return
end subroutine atebmisc_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine atebzo_standard(zom,zoh,zoq,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: zom, zoh, zoq
logical, intent(in), optional :: raw
logical mode

if ( diag>=1 ) write(6,*) "Calculate urban roughness lengths"
if (.not.ateb_active) return

mode=.false.
if (present(raw)) mode=raw

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebzo_thread(zom(is:ie),zoh(is:ie),zoq(is:ie),diag,p_g(:,tile),f_g(tile),upack_g(:,tile),ufull_g(tile),raw=mode)
  end if
end do

return
end subroutine atebzo_standard
                             
subroutine atebzo_thread(zom,zoh,zoq,diag,pd,fp,upack,ufull,raw)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: zom, zoh, zoq
real, dimension(ufull) :: workb,workc,workd,zmtmp,zhtmp,zqtmp
real, dimension(ufull) :: sumlzom, sumlzoh, sumlzoq, sumcndzmin
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat
logical, intent(in), optional :: raw
logical mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if ( diag>=2 ) write(6,*) "THREAD: Calculate urban roughness length"
if ( ufull==0 ) return

mode=.false.
if (present(raw)) mode=raw

sumlzom = 0.
sumlzoh = 0.
sumlzoq = 0.
sumcndzmin = 0.
do ifrac = 1,nfrac
  sumlzom = sumlzom + pd(ifrac)%frac_sigma/pd(ifrac)%lzom**2
  sumlzoh = sumlzoh + pd(ifrac)%frac_sigma/(pd(ifrac)%lzom*pd(ifrac)%lzoh)
  sumlzoq = sumlzoq + pd(ifrac)%frac_sigma/(pd(ifrac)%lzom*pd(ifrac)%lzoh)
  sumcndzmin = sumcndzmin + pd(ifrac)%frac_sigma*pd(ifrac)%cndzmin
end do

if (mode) then
  sumlzom = sqrt(sumlzom)
  sumlzoh = sumlzoh/sumlzom
  sumlzoq = sumlzoq/sumlzom
  zom=unpack(sumcndzmin*exp(-1./sumlzom),upack,zom)
  zoh=unpack(max(sumcndzmin*exp(-1./sumlzoh),zr),upack,zoh)
  zoq=unpack(max(sumcndzmin*exp(-1./sumlzoq),zr),upack,zoq)
else 
  ! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
  zmtmp=pack(zom,upack)
  zhtmp=pack(zoh,upack)
  zqtmp=pack(zoq,upack)
  workb=sqrt((1.-fp%sigmau)/log(sumcndzmin/zmtmp)**2+fp%sigmau*sumlzom)
  workc=(1.-fp%sigmau)/(log(sumcndzmin/zmtmp)*log(sumcndzmin/zhtmp))+fp%sigmau*sumlzoh
  workc=workc/workb
  workd=(1.-fp%sigmau)/(log(sumcndzmin/zmtmp)*log(sumcndzmin/zqtmp))+fp%sigmau*sumlzoq
  workd=workd/workb
  workb=sumcndzmin*exp(-1./workb)
  workc=max(sumcndzmin*exp(-1./workc),zr)
  workd=max(sumcndzmin*exp(-1./workd),zr)
  zom=unpack(workb,upack,zom)
  zoh=unpack(workc,upack,zoh)
  zoq=unpack(workd,upack,zoq)
  if (minval(workc)<=zr) write(6,*) "WARN: minimum zoh reached"
end if

return
end subroutine atebzo_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends the urban drag coeff
!

subroutine atebcd_standard(cduv,cdtq,diag,raw)
 
implicit none
 
integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: cduv, cdtq
logical, intent(in), optional :: raw
logical outmode

if ( diag>=1 ) write(6,*) "Calculate urban drag coeff"
if (.not.ateb_active) return

outmode=.false.
if (present(raw)) outmode=raw

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebcd_thread(cduv(is:ie),cdtq(is:ie),diag,p_g(:,tile),f_g(tile),upack_g(:,tile),ufull_g(tile),raw=outmode)
  end if
end do

return
end subroutine atebcd_standard

subroutine atebcd_thread(cduv,cdtq,diag,pd,fp,upack,ufull,raw)
 
implicit none
 
integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: cduv, cdtq
real, dimension(ufull) :: ctmp, dtmp
logical, intent(in), optional :: raw
logical outmode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd
 
if (diag>=2) write(6,*) "THREAD: Calculate urban drag coeff"
if ( ufull==0 ) return
 
outmode=.false.
if (present(raw)) outmode=raw
 

dtmp = 0.
do ifrac = 1,nfrac
  dtmp = dtmp + pd(ifrac)%cduv*pd(ifrac)%frac_sigma  
end do
if ( outmode ) then
  ctmp=dtmp
else
  ctmp=pack(cduv,upack)  
  ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
end if
cduv=unpack(ctmp,upack,cduv)
 

dtmp = 0.
do ifrac = 1,nfrac
  dtmp = dtmp + pd(ifrac)%cdtq*pd(ifrac)%frac_sigma  
end do
if ( outmode ) then
  ctmp=dtmp
else
  ctmp=pack(cdtq,upack)  
  ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
end if
cdtq=unpack(ctmp,upack,cdtq)
 
return
end subroutine atebcd_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is for hydrological outputs
!
 
subroutine atebhydro_standard(hydroout,mode,diag)

implicit none
 
integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: hydroout
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Calculate hydrological outputs"
if (.not.ateb_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebhydro_thread(hydroout(is:ie),mode,diag,p_g(:,tile),f_g(tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine atebhydro_standard

subroutine atebhydro_thread(hydroout,mode,diag,pd,fp,upack,ufull)
 
implicit none
 
integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: hydroout
real, dimension(ufull) :: ctmp, dtmp
character(len=*), intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd
 
if ( diag>=2 ) write(6,*) "THREAD: Calculate hydrological outputs"
if ( ufull==0 ) return
 
select case(mode)
  case("snowmelt")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%snowmelt*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("soilmoisture")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%soilmoist*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp =(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("soilwet")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%soilwetness*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp =(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("delsoilmoist")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%delsoilmoist*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("rootmoisture")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%rootmoistc*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp =(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("irrigation")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%irrig*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("transpirationveg")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%transveg*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout = unpack(ctmp, upack, hydroout)
  case("aeroconductionveg")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%acond_vegw*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("delintercept")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%delintercept*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("swe")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%swe*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)    
  case("surfstor")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%surfstor*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)   
  case("snowfrac")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%snowfrac*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout) 
  case("delswe")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%delswe*pd(ifrac)%frac_sigma  
    end do    
    ctmp=pack(hydroout,upack)
    ctmp=(1.-fp%sigmau)*ctmp+fp%sigmau*dtmp
    hydroout=unpack(ctmp,upack,hydroout)
  case("surfrunoff")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%surfrunoff*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout = unpack(ctmp, upack, hydroout)
  case("soilwaterfrac")
    dtmp = 0.
    do ifrac = 1,nfrac
      dtmp = dtmp + pd(ifrac)%soilwater*pd(ifrac)%frac_sigma
    end do  
    ctmp = pack(hydroout, upack)
    ctmp = (1.-fp%sigmau)*ctmp + fp%sigmau*dtmp
    hydroout = unpack(ctmp, upack, hydroout)
  case default
    write(6,*) "ERROR: Unknown atebhydro mode ",trim(mode)
    stop
end select
 
return
end subroutine atebhydro_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is for temperature outputs
!

subroutine atebavetemp(tempout,mode,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(:), intent(inout) :: tempout
character(len=*), intent(in) :: mode

if ( diag>=1 ) write(6,*) "Calculate temperature outputs"
if (.not.ateb_active) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebavetemp_thread(tempout(is:ie),mode,diag,roof_g(:,tile),walle_g(:,tile),wallw_g(:,tile),road_g(:,tile),slab_g(:,tile), &
                            intm_g(:,tile),p_g(:,tile),upack_g(:,tile),ufull_g(tile))
  end if
end do

return
end subroutine atebavetemp

subroutine atebavetemp_thread(tempout,mode,diag,roof,walle,wallw,road,slab, &
                              intm,pd,upack,ufull)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac, ii
real, dimension(:), intent(inout) :: tempout
real, dimension(ufull) :: ctmp
logical, dimension(:), intent(in) :: upack
character(len=*), intent(in) :: mode
character(len=10) :: teststr
type(facetdata), dimension(nfrac), intent(in) :: roof, walle, wallw, road, slab, intm
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if ( diag>=2 ) write(6,*) "THREAD: Calculate temperature outputs"
if ( ufull==0 ) return

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + roof(ifrac)%nodetemp(:,ii)*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if  
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + walle(ifrac)%nodetemp(:,ii)*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if 
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + wallw(ifrac)%nodetemp(:,ii)*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + road(ifrac)%nodetemp(:,ii)*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + slab(ifrac)%nodetemp(:,ii)*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    ctmp = 0.
    do ifrac = 1,nfrac
      ctmp = ctmp + intm(ifrac)%nodetemp(:,ii)*pd(ifrac)%frac_sigma
    end do
    tempout = unpack( ctmp+urbtemp, upack, tempout )
    return
  end if
end do

write(6,*) "ERROR: unknown option for atebavetmp_thread ",trim(mode)
stop
  
return
end subroutine atebavetemp_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Store fraction of direct radiation
!

subroutine atebfbeam(is,ifin,fbeam,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ifinish,ib,ie
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: fbeam

if ( diag>=1 ) write(6,*) "Assign urban direct beam ratio"
if ( .not.ateb_active ) return

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is                 ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is               ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then
        f_g(tile)%fbeam(ib:ie)=pack(fbeam(jstart:jfinish),upack_g(kstart:kfinish,tile))
      end if
    end if
  end if
end do

return
end subroutine atebfbeam


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use Spitter et al (1986) method to estimate fraction of direct

subroutine atebspitter(is,ifin,fjd,sg,cosin,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ib,ie,ifinish
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: sg,cosin
real, dimension(imax) :: tmpr,tmpk,tmprat
real, dimension(imax) :: lsg,lcosin
real, intent(in) :: fjd
real, parameter :: solcon = 1370.

if ( diag>=1 ) write(6,*) "Diagnose urban direct beam ratio"
if ( .not.ateb_active ) return

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
      
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then

        lsg(ib:ie)   =pack(sg(jstart:jfinish),upack_g(kstart:kfinish,tile))
        lcosin(ib:ie)=pack(cosin(jstart:jfinish),upack_g(kstart:kfinish,tile))

        tmpr(ib:ie)=0.847+lcosin(ib:ie)*(1.04*lcosin(ib:ie)-1.61)
        tmpk(ib:ie)=(1.47-tmpr(ib:ie))/1.66
        where (lcosin(ib:ie)>1.0e-10 .and. lsg(ib:ie)>10.)
          tmprat(ib:ie)=lsg(ib:ie)/(solcon*(1.+0.033*cos(2.*pi*(fjd-10.)/365.))*lcosin(ib:ie))
        elsewhere
          tmprat(ib:ie)=0.
        end where
        where (tmprat(ib:ie)>tmpk(ib:ie))
          f_g(tile)%fbeam(ib:ie)=max(1.-tmpr(ib:ie),0.)
        elsewhere (tmprat(ib:ie)>0.35)
          f_g(tile)%fbeam(ib:ie)=min(1.66*tmprat(ib:ie)-0.4728,1.)
        elsewhere (tmprat(ib:ie)>0.22)
          f_g(tile)%fbeam(ib:ie)=6.4*(tmprat(ib:ie)-0.22)**2
        elsewhere
          f_g(tile)%fbeam(ib:ie)=0.
        end where
        
      end if
    end if
  end if
end do

return
end subroutine atebspitter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contribution to albedo.
! (selected grid points only)

! raw   (.false.=blend, .true.=output only)
! split (0=net albedo, 1=direct albedo, 2=diffuse albedo)

subroutine atebalb1(is,ifin,alb,diag,raw,split)

implicit none

integer, intent(in) :: is,ifin,diag
integer ucount,ib,ie,ifinish,albmode
integer tile, js, je, kstart, kfinish, jstart, jfinish
integer ifrac
integer, intent(in), optional :: split
real, dimension(ifin), intent(inout) :: alb
real, dimension(imax) :: ualb,utmp
real, dimension(imax) :: dumfbeam,snowdeltac,snowdeltar
real, dimension(imax) :: tmpalb,wallpsi,roadpsi
real, dimension(imax) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
logical, intent(in), optional :: raw
logical outmode

if ( diag>=1 ) write(6,*) "Calculate urban albedo (broad)"
if ( .not.ateb_active ) return

outmode=.false.
if (present(raw)) outmode=raw

albmode=0 ! net albedo
if (present(split)) albmode=split

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
      
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then

        ucount = ie - ib + 1  

        select case(albmode)
          case default ! net albedo
            dumfbeam(ib:ie)=f_g(tile)%fbeam(ib:ie)
          case(1)      ! direct albedo
            dumfbeam(ib:ie)=1.
          case(2)      ! diffuse albedo
            dumfbeam(ib:ie)=0.
        end select

        ualb(ib:ie) = 0.  
        do ifrac = 1,nfrac  
  
          ! roof
          snowdeltar(ib:ie)=rfhyd_g(ifrac,tile)%snow(ib:ie)/(rfhyd_g(ifrac,tile)%snow(ib:ie)+maxrfsn)
  
          ! canyon
          snowdeltac(ib:ie)=rdhyd_g(ifrac,tile)%snow(ib:ie)/(rdhyd_g(ifrac,tile)%snow(ib:ie)+maxrdsn)
          call getswcoeff(sg_roof(ib:ie),sg_vegr(ib:ie),sg_road(ib:ie),sg_walle(ib:ie),sg_wallw(ib:ie),      &
                          sg_vegc(ib:ie),sg_rfsn(ib:ie),sg_rdsn(ib:ie),wallpsi(ib:ie),roadpsi(ib:ie),        &
                          f_g(tile)%hwratio(ib:ie),f_g(tile)%effhwratio(ib:ie),f_g(tile)%vangle(ib:ie),      &
                          f_g(tile)%hangle(ib:ie),dumfbeam(ib:ie),cnveg_g(tile)%sigma(ib:ie),                &
                          f_road(tile)%alpha(ib:ie),cnveg_g(tile)%alpha(ib:ie),f_wall(tile)%alpha(ib:ie),    &
                          rdhyd_g(ifrac,tile)%snowalpha(ib:ie),snowdeltac(ib:ie))
          sg_walle(ib:ie)=sg_walle(ib:ie)*f_g(tile)%coeffbldheight(ib:ie)
          sg_wallw(ib:ie)=sg_wallw(ib:ie)*f_g(tile)%coeffbldheight(ib:ie)

          call getnetalbedo(tmpalb(ib:ie),sg_roof(ib:ie),sg_vegr(ib:ie),sg_road(ib:ie),sg_walle(ib:ie),       &
                            sg_wallw(ib:ie),sg_vegc(ib:ie),sg_rfsn(ib:ie),sg_rdsn(ib:ie),                     &
                            f_g(tile)%hwratio(ib:ie),f_g(tile)%sigmabld(ib:ie),rfveg_g(tile)%sigma(ib:ie),    &
                            f_roof(tile)%alpha(ib:ie),rfveg_g(tile)%alpha(ib:ie),                             &
                            cnveg_g(tile)%sigma(ib:ie),f_road(tile)%alpha(ib:ie),                             &
                            f_wall(tile)%alpha(ib:ie),cnveg_g(tile)%alpha(ib:ie),                             &
                            rfhyd_g(ifrac,tile)%snowalpha(ib:ie),rdhyd_g(ifrac,tile)%snowalpha(ib:ie),        &
                            snowdeltar(ib:ie),snowdeltac(ib:ie))

          ualb(ib:ie) = ualb(ib:ie) + tmpalb(ib:ie)*p_g(ifrac,tile)%frac_sigma(ib:ie)
        end do  


        if (outmode) then
          alb(jstart:jfinish)=unpack(ualb(ib:ie),upack_g(kstart:kfinish,tile),alb(jstart:jfinish))
        else
          utmp(ib:ie)=pack(alb(jstart:jfinish),upack_g(kstart:kfinish,tile))
          utmp(ib:ie)=(1.-f_g(tile)%sigmau(ib:ie))*utmp(ib:ie)+f_g(tile)%sigmau(ib:ie)*ualb(ib:ie)
          alb(jstart:jfinish)=unpack(utmp(ib:ie),upack_g(kstart:kfinish,tile),alb(jstart:jfinish))
        end if
        
      end if
    end if
  end if
end do  

return
end subroutine atebalb1

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

subroutine atebnewangle1(is,ifin,cosin,azimuthin,ctimein)

implicit none

integer, intent(in) :: is,ifin
integer ifinish,ib,ie
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, dimension(ifin), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifin), intent(in) :: azimuthin ! azimuthal angle
real, dimension(ifin), intent(in) :: ctimein   ! local hour (0<=ctime<=1)

if (.not.ateb_active) return

ifinish = is + ifin - 1

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then
        f_g(tile)%hangle(ib:ie)=0.5*pi-pack(azimuthin(jstart:jfinish),upack_g(kstart:kfinish,tile))
        f_g(tile)%vangle(ib:ie)=acos(pack(cosin(jstart:jfinish),upack_g(kstart:kfinish,tile)))
        f_g(tile)%ctime(ib:ie)=pack(ctimein(jstart:jfinish),upack_g(kstart:kfinish,tile))
      end if
    end if
  end if
end do  

return
end subroutine atebnewangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM and TAPM
!

subroutine atebccangle(is,ifin,cosin,rlon,rlat,fjd,slag,dt,sdlt)

implicit none

integer, intent(in) :: is,ifin
integer ifinish, ib, ie, iqu
integer tile, js, je, kstart, kfinish, jstart, jfinish
real, intent(in) :: fjd,slag,dt,sdlt
real cdlt, x, y
real, dimension(ifin), intent(in) :: cosin,rlon,rlat
real, dimension(imax) :: hloc,lattmp

! cosin = cosine of zenith angle
! rlon = longitude
! rlat = latitude
! fjd = day of year
! slag = sun lag angle
! sdlt = sin declination of sun

if (.not.ateb_active) return

ifinish = is + ifin - 1

cdlt=sqrt(min(max(1.-sdlt**2,0.),1.))

do tile = 1,ntiles
  js = (tile-1)*imax + 1 ! js:je is the tile portion of 1:ifull
  je = tile*imax         ! js:je is the tile portion of 1:ifull
  if ( ufull_g(tile)>0 ) then
      
    kstart = max( is - js + 1, 1)          ! kstart:kfinish is the requested portion of 1:imax
    kfinish = min( ifinish - js + 1, imax) ! kstart:kfinish is the requested portion of 1:imax
    if ( kstart<=kfinish ) then
      jstart = kstart + js - is             ! jstart:jfinish is the tile portion of 1:ifin
      jfinish = kfinish + js - is           ! jstart:jfinish is the tile portion of 1:ifin
      ib = count(upack_g(1:kstart-1,tile))+1
      ie = count(upack_g(kstart:kfinish,tile))+ib-1
      if ( ib<=ie ) then

        lattmp(ib:ie)=pack(rlat(jstart:jfinish),upack_g(kstart:kfinish,tile))

        ! from CCAM zenith.f
        hloc(ib:ie)=2.*pi*fjd+slag+pi+pack(rlon(jstart:jfinish),upack_g(kstart:kfinish,tile))+dt*pi/86400.
        
        ! estimate azimuth angle
        do iqu = ib,ie
          x=sin(-hloc(iqu))*cdlt
          y=-cos(-hloc(iqu))*cdlt*sin(lattmp(iqu))+cos(lattmp(iqu))*sdlt
          !azimuth=atan2(x,y)
          f_g(tile)%hangle(iqu)=0.5*pi-atan2(x,y)
        end do  
        
        f_g(tile)%vangle(ib:ie)=acos(pack(cosin(jstart:jfinish),upack_g(kstart:kfinish,tile)))
        f_g(tile)%ctime(ib:ie)=min(max(mod(0.5*hloc(ib:ie)/pi-0.5,1.),0.),1.)
        ! calculate weekdayload loading
        f_g(tile)%weekdayload(ib:ie)=1.0
        if ( (mod(floor(fjd),7)==3) .or. (mod(floor(fjd),7)==4) ) then
          f_g(tile)%weekdayload(ib:ie)=0.9  ! weekend = 90%
        end if  
      end if
    end if
  end if
end do  

return
end subroutine atebccangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates screen level diagnostics
!

subroutine atebscrnout(tscrn,qscrn,uscrn,u10,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(ifull), intent(inout) :: tscrn,qscrn,uscrn,u10
logical, intent(in), optional :: raw
logical mode

if (diag>=1) write(6,*) "Calculate urban 2m diagnostics"
if (.not.ateb_active) return

mode=.false.
if (present(raw)) mode=raw

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebscrnout_thread(tscrn(is:ie),qscrn(is:ie),uscrn(is:ie),u10(is:ie),f_g(tile),p_g(:,tile), &
                            upack_g(:,tile),ufull_g(tile),diag,mode)
  end if
end do  

return
end subroutine atebscrnout

subroutine atebscrnout_thread(tscrn,qscrn,uscrn,u10,fp,pd,upack,ufull,diag,mode)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, dimension(:), intent(inout) :: tscrn,qscrn,uscrn,u10
real, dimension(ufull) :: tmp
real, dimension(ufull) :: tmp_tscrn, tmp_qscrn, tmp_uscrn, tmp_u10
logical, intent(in) :: mode
logical, dimension(:), intent(in) :: upack
type(fparmdata), intent(in) :: fp
type(pdiagdata), dimension(nfrac), intent(in) :: pd

if (diag>=1) write(6,*) "THREAD: Calculate urban 2m diagnostics"
if (.not.ateb_active) return

tmp_tscrn = 0.
tmp_qscrn = 0.
tmp_uscrn = 0.
tmp_u10 = 0.
do ifrac = 1,nfrac
  tmp_tscrn = tmp_tscrn + (pd(ifrac)%tscrn+urbtemp)*pd(ifrac)%frac_sigma
  tmp_qscrn = tmp_qscrn + pd(ifrac)%qscrn*pd(ifrac)%frac_sigma
  tmp_uscrn = tmp_uscrn + pd(ifrac)%uscrn*pd(ifrac)%frac_sigma
  tmp_u10 = tmp_u10 + pd(ifrac)%u10*pd(ifrac)%frac_sigma
end do    

if (mode) then
  tscrn(:)=unpack(tmp_tscrn,upack,tscrn(:))
  qscrn(:)=unpack(tmp_qscrn,upack,qscrn(:))
  uscrn(:)=unpack(tmp_uscrn,upack,uscrn(:))
  u10(:)  =unpack(tmp_u10,  upack,u10(:)  )
else
  tmp(:)=pack(tscrn(:),upack)
  tmp(:)=fp%sigmau*tmp_tscrn + (1.-fp%sigmau)*tmp(:)
  tscrn(:)=unpack(tmp(:),upack,tscrn(:))
  tmp(:)=pack(qscrn(:),upack)
  tmp(:)=fp%sigmau*tmp_qscrn + (1.-fp%sigmau)*tmp(:)
  qscrn(:)=unpack(tmp(:),upack,qscrn(:))
  tmp(:)=pack(uscrn(:),upack)
  tmp(:)=fp%sigmau*tmp_uscrn + (1.-fp%sigmau)*tmp(:)
  uscrn(:)=unpack(tmp(:),upack,uscrn(:))
  tmp(:)=pack(u10(:),upack)
  tmp(:)=fp%sigmau*tmp_u10 + (1.-fp%sigmau)*tmp(:)
  u10(:)=unpack(tmp(:),upack,u10(:))
end if

return
end subroutine atebscrnout_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract urban fraction
subroutine atebsigmau(sigu,diag)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, dimension(ifull), intent(out) :: sigu

if (diag>=1) write(6,*) "Calculate urban cover fraction"
sigu=0.
if (.not.ateb_active) return
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    sigu(is:ie)=unpack(f_g(tile)%sigmau,upack_g(:,tile),0.)
  end if
end do

return
end subroutine atebsigmau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for calculating urban flux contribution

! ifull = number of horizontal grid points
! dt = model time step (sec)
! zmin = first model level height (m)
! sg = incoming short wave radiation (W/m^2)
! rg = incoming long wave radiation (W/m^2)
! rnd = incoming rainfall/snowfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level (kg/m^3)
! temp = atmospheric temperature at first model level (K)
! mixr = atmospheric mixing ratio at first model level (kg/kg)
! ps = surface pressure (Pa)
! pa = pressure at first model level (Pa)
! uu = U component of wind speed at first model level (m/s)
! vv = V component of wind speed at first model level (m/s)
! umin = minimum wind speed (m/s)
! ofg = Input/Output sensible heat flux (W/m^2)
! oeg = Input/Output latent heat flux (W/m^2)
! ots = Input/Output radiative/skin temperature (K)
! owf = Input/Output wetness fraction/surface water (%)
! orn = runoff
! oevspsbl = evapotranspiration
! osbl = sublimation
! diag = diagnostic message mode (0=off, 1=basic messages, 2=more detailed messages, etc)

subroutine atebcalc_standard(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,sg,rg,rnd,snd,rho, &
                             temp,mixr,ps,uu,vv,umin,diag,raw)

implicit none

integer, intent(in) :: diag
integer tile, is, ie
real, intent(in) :: dt,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf,orn,oevspsbl,osbl
logical, intent(in), optional :: raw
logical mode

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

if ( .not.ateb_active ) return

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  if ( ufull_g(tile)>0 ) then
    call atebcalc_thread(ofg(is:ie),oeg(is:ie),ots(is:ie),owf(is:ie),orn(is:ie),oevspsbl(is:ie),      &
                        osbl(is:ie),dt,zmin(is:ie),sg(is:ie),rg(is:ie),rnd(is:ie),snd(is:ie),         &
                        rho(is:ie),temp(is:ie),mixr(is:ie),ps(is:ie),uu(is:ie),vv(is:ie),umin,        &
                        f_g(tile),f_intm(tile),f_road(tile),f_roof(tile),f_slab(tile),f_wall(tile),   &
                        intm_g(:,tile),p_g(:,tile),rdhyd_g(:,tile),rfhyd_g(:,tile),rfveg_g(tile),     &
                        road_g(:,tile),roof_g(:,tile),room_g(:,tile),slab_g(:,tile),walle_g(:,tile),  &
                        wallw_g(:,tile),cnveg_g(tile),intl_g(tile),upack_g(:,tile),ufull_g(tile),     &
                        diag,raw=mode)
  end if
end do

return
end subroutine atebcalc_standard

subroutine atebcalc_thread(ofg,oeg,ots,owf,orn,oevspsbl,osbl,dt,zmin,sg,rg,rnd,snd,rho,temp,    &
                    mixr,ps,uu,vv,umin,fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,      &
                    rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,intl,               &
                    upack,ufull,diag,raw)

implicit none

integer, intent(in) :: ufull, diag
integer ifrac
real, intent(in) :: dt, umin
real, dimension(:), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(:), intent(inout) :: ofg,oeg,ots,owf,orn,oevspsbl,osbl
real, dimension(ufull) :: tmp
real, dimension(ufull) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull) :: u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl
real, dimension(ufull) :: tmp_fg, tmp_eg, tmp_ts, tmp_wf, tmp_rn, tmp_evspsbl, tmp_sbl
logical, intent(in), optional :: raw
logical mode
logical, dimension(:), intent(in) :: upack
type(hydrodata), dimension(nfrac), intent(inout) :: rdhyd, rfhyd
type(facetdata), dimension(nfrac), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(pdiagdata), dimension(nfrac), intent(inout) :: pd
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(vegdata), intent(inout) :: rfveg
type(vegdata), intent(inout) :: cnveg
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp

if ( ufull==0 ) return ! no urban grid points

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

! Host model meteorological data
a_zmin=pack(zmin,                 upack)
a_sg  =pack(sg,                   upack)
a_rg  =pack(rg,                   upack)
a_rho =pack(rho,                  upack)
a_temp=pack(temp-urbtemp,         upack)
a_mixr=pack(mixr,                 upack)
a_ps  =pack(ps,                   upack)
a_umag=max(pack(sqrt(uu*uu+vv*vv),upack),umin)
a_udir=pack(atan2(vv,uu),         upack)
a_rnd =pack(rnd-snd,              upack)
a_snd =pack(snd,                  upack)

u_fg = 0.
u_eg = 0.
u_ts = 0.
u_wf = 0.
u_rn = 0.
u_evspsbl = 0.
u_sbl = 0.

do ifrac = 1,nfrac
  !WRITE(*,*),'1111'
  ! Update urban prognostic variables
  call atebeval(tmp_fg,tmp_eg,tmp_ts,tmp_wf,tmp_rn,tmp_evspsbl,tmp_sbl,dt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps, &
                a_umag,a_udir,a_rnd,a_snd,a_zmin,fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm(ifrac),      &
                pd(ifrac),rdhyd(ifrac),rfhyd(ifrac),rfveg,road(ifrac),roof(ifrac),room(ifrac),slab(ifrac),    &
                walle(ifrac),wallw(ifrac),cnveg,intl,ufull,ifrac,diag)

  u_fg = u_fg + tmp_fg*pd(ifrac)%frac_sigma
  u_eg = u_eg + tmp_eg*pd(ifrac)%frac_sigma
  u_ts = u_ts + (tmp_ts+urbtemp)**4*pd(ifrac)%frac_sigma ! combine as longwave radiation
  u_wf = u_wf + tmp_wf*pd(ifrac)%frac_sigma
  u_rn = u_rn + tmp_rn*pd(ifrac)%frac_sigma
  u_evspsbl = u_evspsbl + tmp_evspsbl*pd(ifrac)%frac_sigma
  u_sbl = u_sbl + tmp_sbl*pd(ifrac)%frac_sigma
  
end do

u_ts = u_ts**0.25 ! convert back to temperature
! export urban fluxes on host grid
if (mode) then
  ofg = unpack(u_fg,upack,ofg)
  oeg = unpack(u_eg,upack,oeg)
  ots = unpack(u_ts,upack,ots)
  owf = unpack(u_wf,upack,owf)
  orn = unpack(u_rn,upack,orn)
  oevspsbl = unpack(u_evspsbl,upack,oevspsbl)
  osbl = unpack(u_sbl,upack,osbl)
else
  tmp = pack(ofg,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_fg
  ofg = unpack(tmp,upack,ofg)
  tmp = pack(oeg,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_eg
  oeg = unpack(tmp,upack,oeg)
  tmp = pack(ots,upack)
  tmp = ((1.-fp%sigmau)*tmp**4+fp%sigmau*u_ts**4)**0.25
  ots = unpack(tmp,upack,ots)
  tmp = pack(owf,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_wf
  owf = unpack(tmp,upack,owf)
  tmp = pack(orn,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_rn
  orn = unpack(tmp,upack,orn)
  tmp = pack(oevspsbl,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_evspsbl
  oevspsbl = unpack(tmp,upack,oevspsbl)
  tmp = pack(osbl,upack)
  tmp = (1.-fp%sigmau)*tmp + fp%sigmau*u_sbl
  osbl = unpack(tmp,upack,osbl)
end if

return
end subroutine atebcalc_thread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! urban flux calculations

! Basic loop is:
!  Canyon short wave flux (nrefl reflections, depending on snow)
!  Canyon long wave flux (nrefl reflections precomputed)
!  Calculate canyon facet aerodynamic resistances
!  Solve canyon sensible and latent heat budgets
!    Solve for canyon snow and vegetation energy budgets
!      Calculate canyon vegetation temperature
!      Calculate canyon snow temperature
!    End canyon snow and vegetation energy budgets loop
!    Estimate interior fluxes into the room
!    Solve for interior air energy budget and interior air temperature
!    Calculate air conditioning fluxes
!    Calculate canyon air temperature
!    Calculate canyon water vapor mixing ratio
!  End canyon sensible and latent heat budgets loop
!  Canyon longwave, sensible and latent heat fluxes                    
!  Solve roof snow and vegetation energy budget
!    Roof snow temperature
!    Roof vegetation temperature
!  End roof snow and vegetation energy budgets loop
!  Roof longwave, sensible and latent heat fluxes
!  Solve interior longwave fluxes
!  Update urban roof, road and wall and slab temperatures with conduction model
!  Update water on road and roof surfaces
!  Update canyon and roof snow albedo and density
!  Estimate bulk roughness length for heat
!  Estimate bulk long wave flux and surface temperature
!  Estimate bulk sensible and latent heat fluxes

subroutine atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl,ddt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,  &
                    a_udir,a_rnd,a_snd,a_zmin,fp,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,rdhyd,      &
                    rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,intl,ufull,ifrac,diag)

implicit none

integer, intent(in) :: ufull, ifrac
integer, intent(in) :: diag
integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(in) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull), intent(out) :: u_fg,u_eg,u_ts,u_wf,u_rn,u_evspsbl,u_sbl
real, dimension(ufull) :: ggint_roof,ggint_walle,ggint_wallw,ggint_road,ggint_slab,ggint_intm
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt,garfsn,gardsn
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr,qsata
real, dimension(ufull) :: cu,fgrooftop,egrooftop
real, dimension(ufull) :: we,ww,wr,zolog,a,n,zom,zonet,dis
real, dimension(ufull) :: roofvegwetfac,roadvegwetfac,d_irrigwater
real, dimension(ufull) :: z_on_l, pa, dts, dtt
real, dimension(ufull) :: u_alb, u_melt, evspsbltop, evspsblrooftop, sbltop, sblrooftop
real, dimension(ufull) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ufull) :: rg_roof,rg_road,rg_walle,rg_wallw,rg_vegc,rg_vegr,rg_rfsn,rg_rdsn
real, dimension(ufull) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab,zero_flux
real, dimension(ufull) :: fg_roof,fg_road,fg_walle,fg_wallw,fg_vegc,fg_vegr,fg_rfsn,fg_rdsn
real, dimension(ufull) :: eg_roof,eg_road,eg_vegc,eg_vegr,eg_rfsn,eg_rdsn
real, dimension(ufull) :: acond_roof,acond_road,acond_walle,acond_wallw
real, dimension(ufull) :: acond_vegc,acond_vegr,acond_rfsn,acond_rdsn
real, dimension(ufull) :: abase_road,abase_walle,abase_wallw,abase_vegc,abase_rdsn
real, dimension(ufull) :: d_roofdelta,d_roaddelta,d_vegdeltac,d_vegdeltar,d_rfsndelta,d_rdsndelta
real, dimension(ufull) :: d_tempc,d_mixrc,d_tempr,d_mixrr,d_sigd,d_sigr,d_rfdzmin
real, dimension(ufull) :: d_ac_canyon,d_canyonrgout,d_roofrgout,d_tranc,d_evapc,d_tranr,d_evapr
real, dimension(ufull) :: d_totdepth,d_netemiss,d_netrad,d_topu
real, dimension(ufull) :: d_cwa,d_cw0,d_cww,d_cwr,d_cra,d_crr,d_crw
real, dimension(ufull) :: d_canyontemp,d_canyonmix,d_traf
real, dimension(ufull) :: ggext_roof, ggext_walle, ggext_wallw, ggext_road, ggext_slab, ggext_intm
real, dimension(ufull) :: ggext_impl, ggint_impl, tmp
real, dimension(ufull) :: ac_coeff, ggint
real, dimension(ufull) :: d_ac_inside, d_intgains_bld, int_infilflux
real, dimension(ufull) :: cyc_traffic, cyc_basedemand, cyc_proportion, cyc_translation
real, dimension(ufull) :: int_infilfg, delintercept, rf_delsnow, rd_delsnow
real, dimension(ufull) :: iroomtemp
real, dimension(ufull) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull) :: d_roomstor
real, dimension(ufull,0:nl) :: inodetemp      ! interim tridiagonal solution
real, dimension(ufull,1:4)  :: iskintemp      ! interim internal skin temperature
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(hydrodata), intent(inout) :: rdhyd, rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(vegdata), intent(inout) :: cnveg
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

if ( diag>=1 ) write(6,*) "Evaluating aTEB"
! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  rfhyd%den = (rfhyd%snow*rfhyd%den+a_snd*ddt*minsnowden)/(rfhyd%snow+ddt*a_snd)
  rdhyd%den = (rdhyd%snow*rdhyd%den+a_snd*ddt*minsnowden)/(rdhyd%snow+ddt*a_snd)
  ! reset snow albedo
  rfhyd%snowalpha = maxsnowalpha
  rdhyd%snowalpha = maxsnowalpha
end where

! reset inception of leafe water
pd%delintercept = 0.

! calculate water and snow area cover fractions
d_roofdelta = max(rfhyd%surfwater/maxrfwater,0.)**(2./3.)
d_roaddelta = max(rdhyd%surfwater/maxrdwater,0.)**(2./3.)
d_vegdeltac = max(rdhyd%leafwater/max(maxvwatf*cnveg%lai,1.E-8),0.)**(2./3.)
d_vegdeltar = max(rfhyd%leafwater/max(maxvwatf*rfveg%lai,1.E-8),0.)**(2./3.)
d_rfsndelta = rfhyd%snow/(rfhyd%snow+maxrfsn)
d_rdsndelta = rdhyd%snow/(rdhyd%snow+maxrdsn)

! canyon level air temp and water vapor (displacement height at refheight*building height)
pa      = a_ps*exp(-grav*a_zmin/(rd*(a_temp+urbtemp)))
d_sigd  = a_ps
a       = (d_sigd/pa)**(rd/aircp)
d_tempc = a_temp*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempc,d_sigd)
call getqsat(qsata,a_temp,pa)
d_mixrc = a_mixr*qsatr/qsata

! roof level air temperature and water vapor (displacement height at building height)
d_sigr  = a_ps*exp(-grav*fp%bldheight*(1.-refheight)/(rd*(a_temp+urbtemp)))
a       = (d_sigr/pa)**(rd/aircp)
d_tempr = a_temp*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempr,d_sigr)
d_mixrr = a_mixr*qsatr/qsata

! calculate soil data
d_totdepth = sum(fp_road%depth,2)

! calculate shortwave reflections
! Here we modify the effective canyon geometry to account for in-canyon vegetation
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,fp%hwratio,   &
                fp%effhwratio,fp%vangle,fp%hangle,fp%fbeam,cnveg%sigma,fp_road%alpha,cnveg%alpha,fp_wall%alpha, &
                rdhyd%snowalpha,d_rdsndelta)
sg_walle = sg_walle*fp%coeffbldheight ! shadow due to in-canyon vegetation
sg_wallw = sg_wallw*fp%coeffbldheight ! shadow due to in-canyon vegetation
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

! calculate long wave reflections to nrefl order (pregenerated before canyonflux subroutine)
call getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta, &
                wallpsi,roadpsi,cnveg%sigma,fp_road%emiss,cnveg%emiss,fp_wall%emiss)
                
pd%emiss = d_rfsndelta*snowemiss+(1.-d_rfsndelta)*((1.-rfveg%sigma)*fp_roof%emiss+rfveg%sigma*rfveg%emiss)
pd%emiss = fp%sigmabld*pd%emiss+(1.-fp%sigmabld)*(2.*fp_wall%emiss*fp%effhwratio*d_cwa+d_netemiss*d_cra) ! diagnostic only

! estimate bulk in-canyon surface roughness length
dis   = max(max(max(0.1*fp%coeffbldheight*fp%bldheight,zocanyon+0.2),cnveg%zo+0.2),zosnow+0.2)
zolog = 1./sqrt(d_rdsndelta/log(dis/zosnow)**2+(1.-d_rdsndelta)*(cnveg%sigma/log(dis/cnveg%zo)**2  &
       +(1.-cnveg%sigma)/log(dis/zocanyon)**2))
zonet = dis*exp(-zolog)

! estimate overall urban roughness length
zom = zomratio*fp%bldheight
where ( zom*fp%sigmabld<zonet*(1.-fp%sigmabld) ) ! MJT suggestion
  zom = zonet
end where
n   = rdhyd%snow/(rdhyd%snow+maxrdsn+0.408*grav*zom)   ! snow cover for urban roughness calc (Douville, et al 1995)
zom = (1.-n)*zom + n*zosnow                            ! blend urban and snow roughness lengths (i.e., snow fills canyon)

! Calculate distance from atmosphere to displacement height
d_rfdzmin = max(a_zmin-fp%bldheight,zoroof+0.2,rfveg%zo+0.2) ! distance to roof displacement height
pd%cndzmin = max(a_zmin-refheight*fp%bldheight,1.5,zom+0.2)  ! distance to canyon displacement height
pd%lzom    = log(pd%cndzmin/zom)

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
  case(1) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    ! estimate wind speed along canyon surfaces
    call getincanwind(we,ww,wr,a_udir,zonet,fp,ufull)
    dis=max(0.1*fp%coeffbldheight*fp%bldheight,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    ! calculate terms for turbulent fluxes
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_walle=a*we                 ! east wall bulk transfer
    abase_wallw=a*ww                 ! west wall bulk transfer
    dis=max(0.1*fp%coeffbldheight*fp%bldheight,zocanyon+0.2,cnveg%zo+0.2,zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_road=a*wr                  ! road bulk transfer
    zolog=log(dis/cnveg%zo)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_rdsn=a*wr                  ! road snow bulk transfer
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
    call getincanwindb(we,ww,wr,a_udir,zonet,fp,ufull)
    dis=max(0.1*fp%coeffbldheight*fp%bldheight,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_walle=a*we                 ! east wall bulk transfer
    abase_wallw=a*ww                 ! west wall bulk transfer
    dis=max(0.1*fp%coeffbldheight*fp%bldheight,zocanyon+0.2,cnveg%zo+0.2,zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_road=a*wr                  ! road bulk transfer
    zolog=log(dis/cnveg%zo)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_rdsn=a*wr                  ! road snow bulk transfer
end select
  
call getdiurnal(fp%ctime,fp%weekdayload,cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation)
! remove statistical energy use diurnal adjustments
if (statsmeth==0) then
  cyc_basedemand=1.
  cyc_proportion=1.
  cyc_translation=0.
end if
! traffic sensible heat flux
pd%traf = fp%trafficfg*cyc_traffic
d_traf = pd%traf/(1.-fp%sigmabld)
! internal gains sensible heat flux
d_intgains_bld = (fp%intmassn+1.)*fp%intgains_flr*cyc_basedemand ! building internal gains 
pd%intgains_full= fp%sigmabld*d_intgains_bld                     ! full domain internal gains
! fractional interior spaces
select case( intairtmeth + (ifrac-1)*10 )
  case(0,1)
    pd%frac_sigma = 1.
  case(2)  ! heating and cooling
    pd%frac_sigma = min( fp%heatprop, fp%coolprop )*cyc_proportion  
  case(12) ! heating or cooling
    pd%frac_sigma = (max( fp%heatprop, fp%coolprop ) - min( fp%heatprop, fp%coolprop ))*cyc_proportion
  case(22) ! no conditioning
    pd%frac_sigma = 1. - max( fp%heatprop, fp%coolprop )*cyc_proportion
end select      
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
                 cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,ifrac,ufull,diag)
zero_flux = 0.
! calculate roof fluxes (fg_roof updated in solvetridiag)
eg_roof = 0. ! For cray compiler
call solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
               sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
               sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
               a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,         &
               d_sigr,ddt,fp_roof,rfhyd,rfveg,roof,fp,ufull)


ggext_roof = (1.-d_rfsndelta)*(sg_roof+rg_roof-eg_roof+aircp*a_rho*d_tempr*acond_roof) &
              +d_rfsndelta*garfsn
ggext_walle= sg_walle+rg_walle+aircp*a_rho*d_canyontemp*acond_walle*fp%coeffbldheight
ggext_wallw= sg_wallw+rg_wallw+aircp*a_rho*d_canyontemp*acond_wallw*fp%coeffbldheight
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

if ( intairtmeth/=0 ) then
  select case(lwintmeth)
    case(0) ! interior longwave off
      rgint_slab  = 0.
      rgint_wallw = 0.
      rgint_roof  = 0.
      rgint_walle = 0.
    case(1) ! explicit longwave, no intermediate facet temperature calc (can cause oscillation)
      iskintemp(:,1) = slab%nodetemp(:,nl)
      iskintemp(:,2) = wallw%nodetemp(:,nl)
      iskintemp(:,3) = roof%nodetemp(:,nl)
      iskintemp(:,4) = walle%nodetemp(:,nl)
      ! calculate interior lwfluxes using last timestep temps lw (rg) fluxes
      call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle,    &
                       iskintemp,fp,intl,ufull)

    case(2) ! explicit longwave using update skin temperature from intermediate facet temp calc
      rgint_slab  = 0.
      rgint_wallw = 0.
      rgint_roof  = 0.
      rgint_walle = 0.

      !!!! slab !!!!
      ggint_impl = cvcoeff_slab
      ggint = -cvcoeff_slab*iroomtemp
      ggext_impl = 0.  ! zero_flux boundary condition
      call solvetridiag(ggext_slab,ggint,rgint_slab,ggext_impl,ggint_impl,slab%nodetemp,  &
                        inodetemp,fp_slab%depth,fp_slab%volcp,fp_slab%lambda,ddt,ufull)
      iskintemp(:,1) = inodetemp(:,nl)

      !!!! wallw !!!!
      ggint_impl = cvcoeff_wallw
      ggint = -cvcoeff_wallw*iroomtemp
      ggext_impl = aircp*a_rho*acond_wallw*fp%coeffbldheight ! implicit, later update fg_wallw with final T
      call solvetridiag(ggext_wallw,ggint,rgint_wallw,ggext_impl,ggint_impl,wallw%nodetemp,  &
                        inodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)
      iskintemp(:,2) = inodetemp(:,nl)

      !!!! roof !!!!
      ggint_impl = cvcoeff_roof
      ggint = -cvcoeff_roof*iroomtemp
      ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! implicit, later update fg_roof with final T
      call solvetridiag(ggext_roof,ggint,rgint_roof,ggext_impl,ggint_impl,roof%nodetemp,  &
                        inodetemp,fp_roof%depth,fp_roof%volcp,fp_roof%lambda,ddt,ufull)
      iskintemp(:,3) = inodetemp(:,nl)

      !!!! walle !!!!
      ggint_impl = cvcoeff_walle
      ggint = -cvcoeff_walle*iroomtemp
      ggext_impl = aircp*a_rho*acond_walle*fp%coeffbldheight ! implicit, later update fg_walle with final T
      call solvetridiag(ggext_walle,ggint,rgint_walle,ggext_impl,ggint_impl,walle%nodetemp,  &
                        inodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)         
      iskintemp(:,4) = inodetemp(:,nl)

      ! calculate interior lwfluxes using iskintemp
      call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle,    &
                           iskintemp,fp,intl,ufull)
  end select
end if


! solve final conduction through facets with explicit interior, implicit exterior boundary
ggint_impl = cvcoeff_roof
ggint = -cvcoeff_roof*iroomtemp
ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! later update fg_roof with final roof skin T

call solvetridiag(ggext_roof,ggint,rgint_roof,ggext_impl,ggint_impl,roof%nodetemp,     &
                  roof%nodetemp,fp_roof%depth,fp_roof%volcp,fp_roof%lambda,ddt,ufull)

ggint_impl = cvcoeff_walle
ggint = -cvcoeff_walle*iroomtemp
ggext_impl = aircp*a_rho*acond_walle*fp%coeffbldheight ! later update fg_walle with final walle skin T
call solvetridiag(ggext_walle,ggint,rgint_walle,ggext_impl,ggint_impl,walle%nodetemp,  &
                  walle%nodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)

ggint_impl = cvcoeff_wallw
ggint = -cvcoeff_wallw*iroomtemp
ggext_impl = aircp*a_rho*acond_wallw*fp%coeffbldheight ! later update fg_wallw with final wallw skin T
call solvetridiag(ggext_wallw,ggint,rgint_wallw,ggext_impl,ggint_impl,wallw%nodetemp,  &
                  wallw%nodetemp,fp_wall%depth,fp_wall%volcp,fp_wall%lambda,ddt,ufull)

ggint_impl = 0.  ! zero_flux boundary condition
ggint_road = 0.  ! zero_flux boundary condition
ggext_impl = (1.-d_rdsndelta)*aircp*a_rho*acond_road ! later update fg_road with final road skin T
call solvetridiag(ggext_road,ggint_road,zero_flux,ggext_impl,ggint_impl,road%nodetemp,   &
                  road%nodetemp,fp_road%depth,fp_road%volcp,fp_road%lambda,ddt,ufull)

! calculate internal facet conduction
if ( intairtmeth/=0 ) then
  ggint_impl = cvcoeff_slab
  ggint = -cvcoeff_slab*iroomtemp
  ggext_slab = 0.  ! zero_flux boundary condition
  ggext_impl = 0.  ! zero_flux boundary condition
  call solvetridiag(ggext_slab,ggint,rgint_slab,ggext_impl,ggint_impl,slab%nodetemp,    &
                    slab%nodetemp,fp_slab%depth,fp_slab%volcp,fp_slab%lambda,ddt,ufull)
  if ( intmassmeth/=0 ) then
    ggint_impl = cvcoeff_intm2
    ggint = -cvcoeff_intm2*iroomtemp
    ggext_impl = cvcoeff_intm1
    call solvetridiag(ggext_intm,ggint,zero_flux,ggext_impl,ggint_impl,intm%nodetemp,    &
                      intm%nodetemp,fp_intm%depth,fp_intm%volcp,fp_intm%lambda,ddt,ufull)
  end if
end if


! implicit update for fg to improve stability for thin layers
fg_roof  = aircp*a_rho*(roof%nodetemp(:,0)-d_tempr)*acond_roof
fg_walle = aircp*a_rho*(walle%nodetemp(:,0)-d_canyontemp)*acond_walle*fp%coeffbldheight
fg_wallw = aircp*a_rho*(wallw%nodetemp(:,0)-d_canyontemp)*acond_wallw*fp%coeffbldheight
fg_road  = aircp*a_rho*(road%nodetemp(:,0)-d_canyontemp)*acond_road
if ( intairtmeth/=0 ) then
  if ( intmassmeth/=0 ) then  
    ggext_intm = -cvcoeff_intm1*(intm%nodetemp(:,0)-iroomtemp)
  else
    ggext_intm = 0.  
  end if
else
  ggext_intm = 0.  
end if  

! implicit update for interior model to improve stability
ggint_roof  = cvcoeff_roof*(roof%nodetemp(:,nl)-iroomtemp)
ggint_walle = cvcoeff_walle*(walle%nodetemp(:,nl)-iroomtemp)
ggint_wallw = cvcoeff_wallw*(wallw%nodetemp(:,nl)-iroomtemp)
if ( intairtmeth/=0 ) then
  ggint_slab = cvcoeff_slab*(slab%nodetemp(:,nl)-iroomtemp)  
  if ( intmassmeth/=0 ) then
    ggint_intm = cvcoeff_intm2*(intm%nodetemp(:,nl)-iroomtemp) 
  else
    ggint_intm = 0.  
  end if    
else
  ggint_slab = 0.  
  ggint_intm = 0.
end if
d_ac_inside = -ggint_roof - ggint_slab + (ggext_intm - ggint_intm)*real(fp%intmassn) & 
              - (ggint_walle+ggint_wallw)*(fp%bldheight/fp%bldwidth) &
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

! update canyon flux
fgtop = fp%hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_canyon - int_infilfg


! calculate water/snow budgets for road surface
call updatewater(ddt,rdhyd%surfwater,rdhyd%soilwater,rdhyd%leafwater,rdhyd%snow,    &
                     rdhyd%den,rdhyd%snowalpha,rdsnmelt,a_rnd,a_snd,eg_road,        &
                     eg_rdsn,d_tranc,d_evapc,d_totdepth, cnveg%lai,wbrelaxc,        &
                     fp%sfc,fp%swilt,d_irrigwater,delintercept,rd_delsnow,ufull)
! record grid irrigation road component [kg m-2 s-1]
pd%irrig = d_irrigwater*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)/ddt
! record leaf intercept component [kg m-2 s-1]
pd%delintercept = pd%delintercept + delintercept*cnveg%sigma

! calculate water/snow budgets for roof surface
call updatewater(ddt,rfhyd%surfwater,rfhyd%soilwater,rfhyd%leafwater,rfhyd%snow,     &
                     rfhyd%den,rfhyd%snowalpha,rfsnmelt,a_rnd,a_snd,eg_roof,         &
                     eg_rfsn,d_tranr,d_evapr,fp%rfvegdepth,rfveg%lai,wbrelaxr,       &
                     fp%sfc,fp%swilt,d_irrigwater,delintercept,rf_delsnow,ufull)
! add grid irrigation roof component [kg m-2 s-1]
pd%irrig = pd%irrig + d_irrigwater*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld/ddt
! record leaf intercept component [kg m-2 s-1]
pd%delintercept = pd%delintercept + delintercept*rfveg%sigma

! area weighted vegetation aerodynamic conductance
pd%acond_vegw = (acond_vegc*cnveg%sigma*(1.-fp%sigmabld) + acond_vegr*rfveg%sigma*fp%sigmabld)/  & 
               (rfveg%sigma*fp%sigmabld + cnveg%sigma*(1.-fp%sigmabld))

! calculate runoff (leafwater runoff already accounted for in precip reaching canyon floor) [kg m-2]
u_rn = max(rfhyd%surfwater-maxrfwater,0.)*fp%sigmabld*(1.-rfveg%sigma)                   & ! roof surfwater runnoff 
      +max(rdhyd%surfwater-maxrdwater,0.)*(1.-fp%sigmabld)*(1.-cnveg%sigma)              & ! road surfwater runnoff
      +max(rfhyd%snow-maxrfsn,0.)*fp%sigmabld                                            & ! roof snow removal     
      +max(rdhyd%snow-maxrdsn,0.)*(1.-fp%sigmabld)                                       & ! road snow removal     
      +max(rfhyd%soilwater-fp%ssat,0.)*waterden*fp%rfvegdepth*rfveg%sigma*fp%sigmabld    & ! roofveg saturation runnoff
      +max(rdhyd%soilwater-fp%ssat,0.)*waterden*d_totdepth*cnveg%sigma*(1.-fp%sigmabld)    ! canyveg saturation runnoff

! diagnose surface drainage rate [kg m-2 s-1]
pd%surfrunoff = u_rn/ddt

! diagnose total vegetation transpiration [kg m-2 s-1]
pd%transveg = (d_tranc*cnveg%sigma*(1.-fp%sigmabld) + d_tranr*rfveg%sigma*fp%sigmabld)/lv

! remove round-off problems
rdhyd%soilwater(1:ufull) = min(max(rdhyd%soilwater(1:ufull),fp%swilt),fp%ssat)
rfhyd%soilwater(1:ufull) = min(max(rfhyd%soilwater(1:ufull),fp%swilt),fp%ssat)
rfhyd%surfwater(1:ufull) = min(max(rfhyd%surfwater(1:ufull),0.),maxrfwater)
rdhyd%surfwater(1:ufull) = min(max(rdhyd%surfwater(1:ufull),0.),maxrdwater)
rdhyd%leafwater(1:ufull) = min(max(rdhyd%leafwater(1:ufull),0.),maxvwatf*cnveg%lai)
rfhyd%leafwater(1:ufull) = min(max(rfhyd%leafwater(1:ufull),0.),maxvwatf*rfveg%lai)
rfhyd%snow(1:ufull)      = min(max(rfhyd%snow(1:ufull),0.),maxrfsn)
rdhyd%snow(1:ufull)      = min(max(rdhyd%snow(1:ufull),0.),maxrdsn)
rfhyd%den(1:ufull)       = min(max(rfhyd%den(1:ufull),minsnowden),maxsnowden)
rdhyd%den(1:ufull)       = min(max(rdhyd%den(1:ufull),minsnowden),maxsnowden)
rfhyd%snowalpha(1:ufull) = min(max(rfhyd%snowalpha(1:ufull),minsnowalpha),maxsnowalpha)
rdhyd%snowalpha(1:ufull) = min(max(rdhyd%snowalpha(1:ufull),minsnowalpha),maxsnowalpha)

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
u_wf = fp%sigmabld*(1.-d_rfsndelta)*((1.-rfveg%sigma)*d_roofdelta       &
      +rfveg%sigma*((1.-d_vegdeltar)*roofvegwetfac+d_vegdeltar))        &
      +(1.-fp%sigmabld)*(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta  &
      +cnveg%sigma*((1.-d_vegdeltac)*roadvegwetfac+d_vegdeltac))

pd%snowmelt = fp%sigmabld*rfsnmelt + (1.-fp%sigmabld)*rdsnmelt
u_melt = lf*(fp%sigmabld*d_rfsndelta*rfsnmelt + (1.-fp%sigmabld)*d_rdsndelta*rdsnmelt)

! calculate average snow temperature and Snow water equivilent
pd%snowt = ( (rfsntemp+urbtemp)*d_rfsndelta*(1.-fp%sigmabld) + (rdsntemp+urbtemp)*d_rdsndelta*fp%sigmabld ) / &
           max( d_rfsndelta*(1.-fp%sigmabld) + d_rdsndelta*fp%sigmabld, 1.e-10 )
pd%swe = rfhyd%snow*d_rfsndelta*(1.-fp%sigmabld) + rdhyd%snow*d_rdsndelta*fp%sigmabld
pd%snowfrac = d_rfsndelta*(1.-fp%sigmabld) + d_rdsndelta*fp%sigmabld
pd%salbedo = ( rfhyd%snowalpha*d_rfsndelta*(1.-fp%sigmabld) + rdhyd%snowalpha*d_rdsndelta*fp%sigmabld ) / &
           max( d_rfsndelta*(1.-fp%sigmabld) + d_rdsndelta*fp%sigmabld, 1.e-10 ) 
pd%delswe = rf_delsnow*d_rfsndelta*(1.-fp%sigmabld) + rd_delsnow*d_rdsndelta*fp%sigmabld

! calculate total liquid water storage
pd%surfstor = rdhyd%surfwater*(1.-fp%sigmabld)*(1.-d_rfsndelta)*(1.-rfveg%sigma) &
            + rdhyd%surfwater*fp%sigmabld*(1.-d_rdsndelta)*(1.-cnveg%sigma)

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

! calculate screen level diagnostics
call scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,a,rdsntemp,zonet, &
              cnveg,fp,pd,rdhyd,road,ufull)

call energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                   sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                   rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                   eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                   a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,                  &
                   ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                   ggext_intm,ggint_slab,ggint_intm,d_intgains_bld,        &
                   int_infilflux,d_ac_inside,fp,ddt,                       &
                   cnveg,fp_intm,fp_road,fp_roof,fp_slab,fp_wall,intm,pd,  &
                   rfveg,road,roof,room,slab,walle,wallw,ufull)

return
end subroutine atebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tridiagonal solver for temperatures

! This version has an implicit estimate for roof sensible heat flux

! [ ggB ggC         ] [ temp ] = [ ggD ]
! [ ggA ggB ggC     ] [ temp ] = [ ggD ]
! [     ggA ggB ggC ] [ temp ] = [ ggD ]
! [         ggA ggB ] [ temp ] = [ ggD ]

subroutine solvetridiag(ggext,ggint,rgint,ggext_impl,ggint_impl,    &
                        nodetemp,inodetemp,depth,volcp,lambda,ddt,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull),     intent(in)    :: ggext,ggint,rgint  ! surface energy fluxes
real, dimension(ufull),     intent(in)    :: ggext_impl         ! implicit update for stability
real, dimension(ufull),     intent(in)    :: ggint_impl         ! implicit update for stability
real, dimension(ufull,0:nl),intent(in)    :: nodetemp           ! temperature of each node
real, dimension(ufull,0:nl),intent(out)   :: inodetemp          ! updated temperature solution
real, dimension(ufull,nl),  intent(in)    :: depth,volcp,lambda ! facet depth, heat capacity, conductivity
real(kind=8), dimension(ufull,nl)         :: cap,res            ! layer capacitance & resistance
real(kind=8), dimension(ufull,0:nl)       :: ggA,ggB,ggC,ggD    ! tridiagonal matrices
real(kind=8) ggX                                                ! tridiagonal coefficient
real(kind=8), dimension(ufull)            :: ans                ! tridiagonal solution
real, intent(in)                          :: ddt                ! timestep
integer k, iqu

res = real(depth,8)/real(lambda,8)
cap = real(depth,8)*real(volcp,8)
!!!!!!!!! interface conduction !!!!!!!!!!!
ggA(:,1:nl)   = -1._8/res(:,1:nl)
ggB(:,0)      =  1._8/res(:,1) +0.5_8*cap(:,1)/ddt + ggext_impl
ggB(:,1:nl-1) =  1._8/res(:,1:nl-1) +1._8/res(:,2:nl) +0.5_8*(cap(:,1:nl-1) &
                   +cap(:,2:nl))/ddt
ggB(:,nl)     =  1._8/res(:,nl) + 0.5_8*cap(:,nl)/real(ddt,8) + ggint_impl
ggC(:,0:nl-1) = -1._8/res(:,1:nl)
ggD(:,0)      = nodetemp(:,0)*0.5_8*cap(:,1)/real(ddt,8) + ggext
ggD(:,1:nl-1) = nodetemp(:,1:nl-1)*0.5_8*(cap(:,1:nl-1)+cap(:,2:nl))/real(ddt,8)
ggD(:,nl)     = nodetemp(:,nl)*0.5_8*cap(:,nl)/ddt - ggint - rgint
! tridiagonal solver (Thomas algorithm) to solve node temperatures
do k=1,nl
  do iqu = 1,ufull
    ggX = ggA(iqu,k)/ggB(iqu,k-1)
    ggB(iqu,k) = ggB(iqu,k)-ggX*ggC(iqu,k-1)
    ggD(iqu,k) = ggD(iqu,k)-ggX*ggD(iqu,k-1)
  end do  
end do
ans = ggD(:,nl)/ggB(:,nl)
inodetemp(:,nl) = real(ans)
do k=nl-1,0,-1
  ans = (ggD(:,k) - ggC(:,k)*ans)/ggB(:,k)
  inodetemp(:,k) = real(ans)
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
                         slab,walle,wallw,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull), intent(in) :: sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle
real, dimension(ufull), intent(in) :: sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw
real, dimension(ufull), intent(in) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab
real, dimension(ufull), intent(in) :: eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta
real, dimension(ufull), intent(in) :: a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt
real, dimension(ufull), intent(in) :: ggint_roof,ggint_road,ggint_walle,ggint_wallw
real, dimension(ufull), intent(in) :: ggext_intm,ggint_slab,ggint_intm,d_intgains_bld
real, dimension(ufull), intent(in) :: int_infilflux,d_ac_inside
real(kind=8), dimension(ufull) :: d_roofflux,d_walleflux,d_wallwflux,d_roadflux,d_slabflux,d_intmflux,d_roomflux 
real(kind=8), dimension(ufull) :: d_roofstor,d_wallestor,d_wallwstor,d_roadstor,d_slabstor,d_intmstor,d_roomstor
real(kind=8), dimension(ufull) :: d_faceterr
real(kind=8), dimension(ufull) :: d_atmosflux
real(kind=8), dimension(ufull,nl) :: roadstorage_prev, roofstorage_prev, wallestorage_prev, wallwstorage_prev
real(kind=8), dimension(ufull,nl) :: slabstorage_prev, intmstorage_prev
real(kind=8), dimension(ufull,1) :: roomstorage_prev
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_slab, fp_wall
type(facetdata), intent(inout) :: intm
type(vegdata), intent(in) :: cnveg, rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw
type(fparmdata), intent(in) :: fp
type(pdiagdata), intent(inout) :: pd

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
room%storage(:,1) = real(fp%bldheight(:)*aircp*room%nodetemp(:,1),8)
! Sum heat stored in urban materials from layer 1 to nl
! interface conduction
roof%storage(:,:)  = 0.5_8*real(fp_roof%depth(:,:),8)*real(fp_roof%volcp(:,:),8)                        & 
                        *(real(roof%nodetemp(:,0:nl-1),8)+real(roof%nodetemp(:,1:nl),8))
road%storage(:,:)  = 0.5_8*real(fp_road%depth(:,:),8)*real(fp_road%volcp(:,:),8)                        & 
                        *(real(road%nodetemp(:,0:nl-1),8)+real(road%nodetemp(:,1:nl),8))
walle%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8)                        & 
                        *(real(walle%nodetemp(:,0:nl-1),8)+real(walle%nodetemp(:,1:nl),8))
wallw%storage(:,:) = 0.5_8*real(fp_wall%depth(:,:),8)*real(fp_wall%volcp(:,:),8)                        & 
                        *(real(wallw%nodetemp(:,0:nl-1),8)+real(wallw%nodetemp(:,1:nl),8))
slab%storage(:,:)  = 0.5_8*real(fp_slab%depth(:,:),8)*real(fp_slab%volcp(:,:),8)                        & 
                        *(real(slab%nodetemp(:,0:nl-1),8)+real(slab%nodetemp(:,1:nl),8))
intm%storage(:,:)  = 0.5_8*real(fp_intm%depth(:,:),8)*real(fp_intm%volcp(:,:),8)                        & 
                        *(real(intm%nodetemp(:,0:nl-1),8)+real(intm%nodetemp(:,1:nl),8))

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
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB roof facet closure error:", maxval(abs(d_faceterr))
d_roadstor = sum(road%storage-roadstorage_prev,dim=2)/real(ddt,8)
d_roadflux = (1._8-real(d_rdsndelta,8))*(real(sg_road,8)+real(rg_road,8)-real(fg_road,8)-real(eg_road,8)) &
           + real(d_rdsndelta,8)*real(gardsn,8) - real(ggint_road,8)
d_faceterr  = d_roadstor - d_roadflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB road facet closure error:", maxval(abs(d_faceterr))
d_wallestor= sum(walle%storage-wallestorage_prev,dim=2)/real(ddt,8)
d_walleflux= real(sg_walle,8)+real(rg_walle,8)-real(fg_walle,8) - real(ggint_walle,8) - real(rgint_walle,8)
d_faceterr = d_wallestor - d_walleflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB walle facet closure error:", maxval(abs(d_faceterr))
d_wallwstor= sum(wallw%storage-wallwstorage_prev,dim=2)/real(ddt,8)
d_wallwflux= real(sg_wallw,8)+real(rg_wallw,8)-real(fg_wallw,8) - real(ggint_wallw,8) - real(rgint_wallw,8)
d_faceterr = d_wallwstor - d_wallwflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB wallw facet closure error:", maxval(abs(d_faceterr))
d_slabstor = sum(slab%storage-slabstorage_prev,dim=2)/real(ddt,8)
d_slabflux = -real(ggint_slab,8) - real(rgint_slab,8)
d_faceterr = d_slabstor - d_slabflux
pd%surferr = pd%surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB slab facet closure error:", maxval(abs(d_faceterr))
if (intairtmeth/=0) then
  d_intmstor = sum(intm%storage-intmstorage_prev,dim=2)/real(ddt,8)
  d_intmflux = real(ggext_intm,8) - real(ggint_intm,8)
  d_faceterr = d_intmstor - d_intmflux
  pd%surferr = pd%surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB intm facet closure error:", maxval(abs(d_faceterr))
  d_roomstor = (room%storage(:,1)-roomstorage_prev(:,1))/real(ddt,8)
  d_roomflux = real(ggint_roof,8)+real(ggint_slab,8)-real(fp%intmassn,8)*d_intmflux                    & 
            + (real(fp%bldheight,8)/real(fp%bldwidth,8))*(real(ggint_walle,8) + real(ggint_wallw,8))   &
            + real(int_infilflux,8) + real(d_ac_inside,8) + real(d_intgains_bld,8)
  d_faceterr = d_roomstor - d_roomflux
  pd%surferr = pd%surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB room volume closure error:", maxval(abs(d_faceterr))
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

if ( any(abs(pd%atmoserr)>=energytol) ) then
  write(6,*) "aTEB energy not conserved! Atmos. error:", maxval(abs(pd%atmoserr))
end if

return
end subroutine energyclosure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update water prognostic variables for roads and roofs
                            
subroutine updatewater(ddt,surfwater,soilwater,leafwater,snow,den,alpha, &
                       snmelt,a_rnd,a_snd,eg_surf,eg_snow,d_tran,d_evap, &
                       d_totdepth,fp_vegrlai,iwbrelax,                   &
                       fp_sfc,fp_swilt,irrigwater,delintercept,          &
                       delsnow,ufull)

implicit none

integer, intent(in) :: ufull
integer, intent(in) :: iwbrelax
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: surfwater,soilwater,leafwater,snow,den,alpha
real, dimension(ufull), intent(in) :: snmelt,a_rnd,a_snd,eg_surf,eg_snow
real, dimension(ufull), intent(in) :: d_tran,d_evap,d_totdepth,fp_vegrlai
real, dimension(ufull) :: modrnd, oldleafwater
real, dimension(ufull), intent(in) :: fp_sfc, fp_swilt
real, dimension(ufull), intent(out):: irrigwater, delintercept, delsnow

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

subroutine getqsat(qsat,temp,ps)

implicit none

real, dimension(:), intent(in) :: temp
real, dimension(:), intent(in) :: ps
real, dimension(:), intent(out) :: qsat
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix

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
! Estimate roughness length for heat
!

pure subroutine getlna(lna,cd,umag,zmin,ilzom,mode)

implicit none

integer, intent(in) :: mode
integer iqu
real, dimension(:), intent(out) :: lna
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

! Calculate short wave reflections to nrefl order
roadnetalpha=rdsndelta*ird_alpha+(1.-rdsndelta)*((1.-fp_vegsigmac)*fp_roadalpha+fp_vegsigmac*fp_vegalphac)
sg_walle=walles
sg_wallw=wallws
sg_road=roads
do k=1,nrefl
  nwalles=roadnetalpha*wallpsi*roads+fp_wallalpha*(1.-2.*wallpsi)*wallws
  nwallws=roadnetalpha*wallpsi*roads+fp_wallalpha*(1.-2.*wallpsi)*walles
  nroads=fp_wallalpha*(1.-roadpsi)*0.5*(walles+wallws)
  walles=nwalles
  wallws=nwallws
  roads=nroads
  sg_walle=sg_walle+walles
  sg_wallw=sg_wallw+wallws
  sg_road=sg_road+roads
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
                       cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,ifrac,ufull,diag)
implicit none

integer, intent(in) :: ifrac,ufull,diag
integer k,l
real, intent(in)    :: ddt
real, dimension(ufull), intent(inout) :: rg_road,fg_road,eg_road,abase_road
real, dimension(ufull), intent(inout) :: rg_walle,fg_walle,abase_walle
real, dimension(ufull), intent(inout) :: rg_wallw,fg_wallw,abase_wallw
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,abase_vegc
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn
real, dimension(ufull), intent(in) :: sg_road,sg_walle,sg_wallw,sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd,d_sigd,we,ww
real, dimension(ufull), intent(out) :: d_canyontemp,d_canyonmix,d_topu,d_netrad
real, dimension(ufull), intent(out) :: d_ac_inside
real, dimension(ufull), intent(inout) :: d_tempc,d_mixrc
real, dimension(ufull), intent(inout) :: d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_canyon,d_traf
real, dimension(ufull), intent(inout) :: d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr
real, dimension(ufull), intent(inout) :: d_totdepth,d_intgains_bld
real, dimension(ufull), intent(out) :: fgtop,egtop,evspsbltop,sbltop,int_infilflux
real, dimension(ufull), intent(out) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(ufull), intent(out) :: int_infilfg
real, dimension(ufull), intent(out) :: iroomtemp
real, dimension(ufull), intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull), intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull), intent(out) :: d_roomstor
real, dimension(ufull) :: newval,sndepth,snlambda,ldratio,roadqsat,vegqsat,rdsnqsat
real, dimension(ufull) :: cu,topinvres,dts,dtt,cduv,z_on_l,dumroaddelta,dumvegdelta,res
real, dimension(ufull) :: effwalle,effwallw,effroad,effrdsn,effvegc
real, dimension(ufull) :: aa,bb,cc,dd,ee,ff
real, dimension(ufull) :: lwflux_walle_road, lwflux_wallw_road, lwflux_walle_rdsn, lwflux_wallw_rdsn
real, dimension(ufull) :: lwflux_walle_vegc, lwflux_wallw_vegc
real, dimension(ufull) :: cyc_proportion,cyc_translation
real, dimension(ufull,2) :: evct,evctx,oldval
real, dimension(ufull) :: d_canyontemp_prev
type(facetparams), intent(in) :: fp_intm, fp_road, fp_roof, fp_wall,fp_slab
type(facetdata), intent(in) :: intm
type(hydrodata), intent(in) :: rdhyd
type(vegdata), intent(inout) :: cnveg, rfveg
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
d_canyontemp    = d_tempc
d_canyontemp_prev = d_tempc
d_canyonmix     = d_mixrc
cnveg%temp      = d_tempc
rdsntemp        = road%nodetemp(:,1)
iroomtemp       = room%nodetemp(:,1)

rdsnmelt        = 0.
dumvegdelta     = 0. ! cray compiler bug
d_netrad=sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                            &
        +(1.-d_rdsndelta)*(1.-cnveg%sigma)*fp_road%emiss*(road%nodetemp(:,0)+urbtemp)**4 &
        +(1.-d_rdsndelta)*cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4)
! Solve for canyon air temperature and water vapour mixing ratio
do l = 1,ncyits
  !  solve for aerodynamical resistance between canyon and atmosphere  
  ! assume zoh=zom when coupling to canyon air temperature
  pd%lzoh = pd%lzom
  dts    = d_canyontemp + (d_canyontemp+urbtemp)*0.61*d_canyonmix
  dtt    = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
  call getinvres(topinvres,cduv,z_on_l,pd%lzoh,pd%lzom,pd%cndzmin,dts,dtt,a_umag,3)
  call gettopu(d_topu,a_umag,z_on_l,fp%bldheight,cduv,pd%cndzmin,fp%hwratio,ufull)
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
  call getqsat(roadqsat,road%nodetemp(:,0),d_sigd)   ! evaluate using pressure at displacement height
  
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
  d_netrad=sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                             &
          +(1.-d_rdsndelta)*((1.-cnveg%sigma)*fp_road%emiss*(road%nodetemp(:,0)+urbtemp)**4 &
          +cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4))
  
  if ( lweff/=1 ) then
    lwflux_walle_road = 0.
    lwflux_wallw_road = 0.
    lwflux_walle_rdsn = 0.
    lwflux_wallw_rdsn = 0.
    lwflux_walle_vegc = 0.
    lwflux_wallw_vegc = 0.
  else
    lwflux_walle_road = sbconst*(fp_road%emiss*(road%nodetemp(:,0)+urbtemp)**4         &
                       -fp_wall%emiss*(walle%nodetemp(:,0)+urbtemp)**4)*(1.-fp%coeffbldheight)
    lwflux_wallw_road = sbconst*(fp_road%emiss*(road%nodetemp(:,0)+urbtemp)**4         &
                       -fp_wall%emiss*(wallw%nodetemp(:,0)+urbtemp)**4)*(1.-fp%coeffbldheight)
    lwflux_walle_rdsn = sbconst*(snowemiss*(rdsntemp+urbtemp)**4                       &
                       -fp_wall%emiss*(walle%nodetemp(:,0)+urbtemp)**4)*(1.-fp%coeffbldheight)
    lwflux_wallw_rdsn = sbconst*(snowemiss*(rdsntemp+urbtemp)**4                       &
                       -fp_wall%emiss*(wallw%nodetemp(:,0)+urbtemp)**4)*(1.-fp%coeffbldheight)
    lwflux_walle_vegc = sbconst*(cnveg%emiss*(cnveg%temp+urbtemp)**4                   &
                       -fp_wall%emiss*(walle%nodetemp(:,0)+urbtemp)**4)*(1.-fp%coeffbldheight)
    lwflux_wallw_vegc = sbconst*(cnveg%emiss*(cnveg%temp+urbtemp)**4                   &
                       -fp_wall%emiss*(wallw%nodetemp(:,0)+urbtemp)**4)*(1.-fp%coeffbldheight)
  end if
  
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
                  d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,                                     &
                  effvegc,effrdsn,ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,                  &
                  lwflux_walle_vegc,lwflux_wallw_vegc,ddt,                                      &
                  cnveg,fp,fp_wall,rdhyd,road,walle,wallw,ufull)
  cnveg%temp = cnveg%temp - 0.5 ! 2nd guess
  rdsntemp   = rdsntemp - 0.5
  do k = 1,nfgits ! sectant
    evctx = evct
    call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                    sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                    a_rg,a_rho,a_rnd,a_snd,                                                       &
                    d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                     &
                    d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,                                     &
                    effvegc,effrdsn,ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,                  &
                    lwflux_walle_vegc,lwflux_wallw_vegc,ddt,                                      &
                    cnveg,fp,fp_wall,rdhyd,road,walle,wallw,ufull)
    evctx = evct - evctx
    where ( abs(evctx(:,1))>tol )
      newval      = max(min(cnveg%temp-alpha*evct(:,1)*(cnveg%temp-oldval(:,1))/evctx(:,1),400.-urbtemp),200.-urbtemp)
      oldval(:,1) = cnveg%temp
      cnveg%temp  = newval
    end where
    where ( abs(evctx(:,2))>tol )
      newval      = max(min(rdsntemp-alpha*evct(:,2)*(rdsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
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
      call interiorflux_simple(ufull,diag,l,ddt,                                 &
                               d_intgains_bld,cyc_translation,                   &
                               d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor,  &
                               iroomtemp,int_infilflux,int_infilfg,              &
                               cvcoeff_roof,cvcoeff_walle,                       &
                               cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,         &
                               cvcoeff_intm2,fp,fp_roof,fp_wall,rfveg,roof,      &
                               walle,wallw,slab,intm)

    case(1,2) ! complex internal physics (per Lipson et al., 2018) or with interior frac      
      call interiorflux_complex(ifrac,ufull,diag,l,ddt,a_rho,we,ww,              &
                                d_intgains_bld,cyc_proportion,d_canyontemp_prev, &
                                d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor, &
                                iroomtemp,int_infilflux,int_infilfg,             &
                                cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,        &
                                cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2,fp,     &
                                fp_roof,fp_wall,fp_slab,fp_intm,roof,walle,      &
                                wallw,slab,intm,room)
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
  d_canyontemp = (aa*d_tempc+bb*rdsntemp+cc*road%nodetemp(:,0)+dd*cnveg%temp+ee*walle%nodetemp(:,0) & 
                +ff*wallw%nodetemp(:,0)+d_traf+d_ac_canyon-int_infilfg)/(aa+bb+cc+dd+ee+ff)

end do

room%nodetemp(:,1) = iroomtemp

! solve for canyon sensible heat flux
fg_walle = aircp*a_rho*(walle%nodetemp(:,0)-d_canyontemp)*acond_walle*fp%coeffbldheight ! canyon vegetation blocks turbulent flux
fg_wallw = aircp*a_rho*(wallw%nodetemp(:,0)-d_canyontemp)*acond_wallw*fp%coeffbldheight ! canyon vegetation blocks turbulent flux
fg_road  = aircp*a_rho*(road%nodetemp(:,0)-d_canyontemp)*acond_road
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
if ( lweff/=2 ) then
  effwalle=fp_wall%emiss*(a_rg*d_cwa+sbconst*(walle%nodetemp(:,0)+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)                & 
                                  +sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_walle=effwalle*fp%coeffbldheight+lwflux_walle_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/fp%hwratio                  &
                                  +lwflux_walle_vegc*(1.-d_rdsndelta)*cnveg%sigma/fp%hwratio                          &
                                  +lwflux_walle_rdsn*d_rdsndelta/fp%hwratio
  effwallw=fp_wall%emiss*(a_rg*d_cwa+sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)                &
                                  +sbconst*(walle%nodetemp(:,0)+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_wallw=effwallw*fp%coeffbldheight+lwflux_wallw_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/fp%hwratio                  &
                                  +lwflux_wallw_vegc*(1.-d_rdsndelta)*cnveg%sigma/fp%hwratio                          &
                                  +lwflux_wallw_rdsn*d_rdsndelta/fp%hwratio
  effroad=fp_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(road%nodetemp(:,0)+urbtemp)**4)                          &
                    +sbconst*fp_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
  rg_road=effroad-lwflux_walle_road-lwflux_wallw_road
  ! outgoing longwave radiation
  ! note that eff terms are used for outgoing longwave radiation, whereas rg terms are used for heat conduction
  d_canyonrgout=a_rg-d_rdsndelta*effrdsn-(1.-d_rdsndelta)*((1.-cnveg%sigma)*effroad+cnveg%sigma*effvegc)            &
                    -fp%hwratio*fp%coeffbldheight*(effwalle+effwallw)
else
  effwalle=fp_wall%emiss*(a_rg*d_cwa+sbconst*(walle%nodetemp(:,0)+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)                & 
                                  +sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_walle=effwalle
  effwallw=fp_wall%emiss*(a_rg*d_cwa+sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*(fp_wall%emiss*d_cw0-1.)                &
                                  +sbconst*(walle%nodetemp(:,0)+urbtemp)**4*fp_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_wallw=effwallw
  effroad=fp_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(road%nodetemp(:,0)+urbtemp)**4)                          &
                    +sbconst*fp_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
  rg_road=effroad
  ! outgoing longwave radiation
  ! note that eff terms are used for outgoing longwave radiation, whereas rg terms are used for heat conduction
  d_canyonrgout=a_rg-d_rdsndelta*effrdsn-(1.-d_rdsndelta)*((1.-cnveg%sigma)*effroad+cnveg%sigma*effvegc)            &
                    -fp%hwratio*(effwalle+effwallw)
end if
!0. = d_rdsndelta*(lwflux_walle_rdsn+lwflux_wallw_rdsn)                            &
!     +(1.-d_rdsndelta)*((1.-cnveg%sigma)*(lwflux_walle_road+lwflux_wallw_road)    &
!    +cnveg%sigma*(lwflux_walle_vegc+lwflux_wallw_vegc))                           &
!    -fp%hwratio*(lwflux_walle_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/fp%hwratio   &
!    +lwflux_walle_vegc*(1.-d_rdsndelta)*cnveg%sigma/fp%hwratio                    &
!     +lwflux_walle_rdsn*d_rdsndelta/fp%hwratio)                                   &
!    - fp%hwratio*(lwflux_wallw_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/fp%hwratio  &
!    +lwflux_wallw_vegc*(1.-d_rdsndelta)*cnveg%sigma/fp%hwratio                    &
!    +lwflux_wallw_rdsn*d_rdsndelta/fp%hwratio)

return
end subroutine solvecanyon


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon/ internal heat fluxes and air temperatures

subroutine interiorflux_simple(ufull,diag,l,ddt,                                      &
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
integer,                    intent(in) :: diag
integer,                    intent(in) :: ufull
integer,                    intent(in) :: l
real,                       intent(in) :: ddt
real, dimension(ufull),     intent(in) :: d_intgains_bld
real, dimension(ufull),     intent(in) :: cyc_translation
real, dimension(ufull),     intent(in) :: d_canyontemp
real, dimension(ufull),     intent(out) :: d_ac_canyon
real, dimension(ufull),     intent(inout) :: d_ac_inside
real, dimension(ufull),     intent(inout) :: iroomtemp
real, dimension(ufull),     intent(out) :: int_infilflux
real, dimension(ufull),     intent(out) :: int_infilfg
real, dimension(ufull),     intent(out) :: cvcoeff_roof, cvcoeff_walle, cvcoeff_wallw, cvcoeff_slab
real, dimension(ufull),     intent(out) :: cvcoeff_intm1, cvcoeff_intm2
real, dimension(ufull),     intent(out) :: d_roomstor

type(fparmdata),   intent(in) :: fp
type(facetparams), intent(in) :: fp_roof,fp_wall
type(facetdata),   intent(in) :: roof,walle,wallw,slab,intm
type(vegdata),     intent(in) :: rfveg

! local variables
real, dimension(ufull) :: ggint_roof,ggint_walle,ggint_wallw,ggint_slab
real, dimension(ufull) :: ac_coeff

if ( diag>=1 ) write(6,*) "Evaluating with fixed internal air temperature"

iroomtemp = fp%bldairtemp + cyc_translation
d_roomstor = 0.

! ---0.1: Estimate internal surface convection coefficients -------------
call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,  & 
                    cvcoeff_intm1,cvcoeff_intm2,iroomtemp,roof,slab,intm,ufull)
! ---0.2: Calculate internal conducted flux (ggint) ----------------------
ggint_slab=0.  ! no slab for intairtmeth=0
ggint_roof  = cvcoeff_roof*(roof%nodetemp(:,nl)-iroomtemp)
ggint_walle = cvcoeff_walle*(walle%nodetemp(:,nl)-iroomtemp)
ggint_wallw = cvcoeff_wallw*(wallw%nodetemp(:,nl)-iroomtemp)

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

subroutine interiorflux_complex(ifrac,ufull,diag,l,ddt,a_rho,d_weu,d_wwu,              &
                                d_intgains_bld,cyc_proportion,d_canyontemp_prev,       &
                                d_canyontemp,d_ac_canyon,d_ac_inside,d_roomstor,       &
                                iroomtemp,int_infilflux,int_infilfg,                   &
                                cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, &
                                cvcoeff_intm1,cvcoeff_intm2,fp,fp_roof,fp_wall,        &
                                fp_slab,fp_intm,roof,walle,wallw,slab,intm,room)

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
integer,                    intent(in) :: ifrac                   ! fractional conditioning integer
integer,                    intent(in) :: diag                    ! diagnostic integer
integer,                    intent(in) :: ufull                   ! tile vector length
integer,                    intent(in) :: l                       ! ncyits loop integer
real,                       intent(in) :: ddt                     ! timestep
real, dimension(ufull),     intent(in) :: a_rho                   ! atmospheric air density
real, dimension(ufull),     intent(in) :: d_weu,d_wwu             ! wind speed wall east and west
real, dimension(ufull),     intent(in) :: d_intgains_bld          ! internal gain flux for building
real, dimension(ufull),     intent(in) :: cyc_proportion          ! cyclic proportion of conditioning units
real, dimension(ufull),     intent(in) :: d_canyontemp            ! canyon air temperature (diagnostic)
real, dimension(ufull),     intent(in) :: d_canyontemp_prev       ! previous canyon air temperature
real, dimension(ufull),     intent(out) :: d_ac_canyon            ! conditioning flux into canyon
real, dimension(ufull),     intent(inout) :: d_ac_inside          ! conditioning flux inside
real, dimension(ufull),     intent(inout) :: iroomtemp            ! update of room air temperature
real, dimension(ufull),     intent(out) :: int_infilflux          ! internal infiltration flux
real, dimension(ufull),     intent(out) :: int_infilfg            ! external infiltration flux
real, dimension(ufull),     intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab
real, dimension(ufull),     intent(out) :: cvcoeff_intm1,cvcoeff_intm2
real, dimension(ufull),     intent(out) :: d_roomstor

type(fparmdata),   intent(in) :: fp
type(facetparams), intent(in) :: fp_roof,fp_wall,fp_slab,fp_intm
type(facetdata),   intent(in) :: roof,walle,wallw,slab,intm
type(facetdata),   intent(inout) :: room

! local variables
integer iqu
real, dimension(ufull) :: ggint_roof, ggint_walle, ggint_wallw ! facet convection fluxes
real, dimension(ufull) :: ggint_slab, ggext_intm, ggint_intm   ! facet convection fluxes
real, dimension(ufull) :: ac_heat_on,ac_cool_on
real, dimension(ufull) :: rm,rf,we,ww,sl,im1,im2,infl
real, dimension(ufull) :: ac_coeff,d_ac_cool,d_ac_heat,ac_load
real, dimension(ufull) :: infl_dynamic,itemp
real, dimension(ufull) :: ac_heat_temp,ac_cool_temp
real, dimension(ufull) :: int_airden
real, dimension(ufull) :: d_openwindows
real, dimension(ufull) :: dummy, canyontemp_ave


if ( diag>=1 ) write(6,*) "Evaluating with varying internal air temperature"

! calculate canyon temperature for fluxes
canyontemp_ave = infilalpha*d_canyontemp + (1.-infilalpha)*d_canyontemp_prev


! fractional interior conditioning cases (1=on)
select case( ifrac )
  case(1)
    ac_heat_on = 1.
    ac_cool_on = 1.
  case(2)
    where ( fp%heatprop>fp%coolprop )
      ac_heat_on = 1.
      ac_cool_on = 0.
    elsewhere
      ac_heat_on = 0.
      ac_cool_on = 1.
    end where
  case(3)
    ac_heat_on = 0.
    ac_cool_on = 0.
end select


! ---1.1: Estimate internal surface convection coefficients -------------
!call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,       & 
!                    cvcoeff_intm1,cvcoeff_intm2,iroomtemp,roof,slab,intm,ufull)
call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,       & 
                    cvcoeff_intm1,cvcoeff_intm2,room%nodetemp(:,1),roof,slab,intm,ufull)


ac_heat_temp = fp%bldairtemp - ac_deltat
ac_cool_temp = fp%bldairtemp + ac_deltat


! ---1.2: Calculate air exchange from infiltration and open windows -----
call calc_openwindows(d_openwindows,fp,iroomtemp,canyontemp_ave,roof,walle,wallw,slab,ufull)
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
itemp = (rm*room%nodetemp(:,1)                    & ! room temperature
       + rf*roof%nodetemp(:,nl)                   & ! roof conduction
       + we*walle%nodetemp(:,nl)                  & ! wall conduction east
       + ww*wallw%nodetemp(:,nl)                  & ! wall conduction west
       + sl*slab%nodetemp(:,nl)                   & ! slab conduction
       + im1*intm%nodetemp(:,0)                   & ! mass conduction side 1
       + im2*intm%nodetemp(:,nl)                  & ! mass conduction side 2
       + infl*canyontemp_ave                      & ! weighted infiltration
       + d_intgains_bld                           & ! internal gains (explicit)
       )/(rm+rf+we+ww+sl+im1+im2+infl)

d_ac_cool = 0.
d_ac_heat = 0.

! heating load
where ( itemp < ac_heat_temp )
  !ac_load = rm*(ac_heat_temp-itemp)
  ac_load = (rm+rf+we+ww+sl+im1+im2+infl)*(ac_heat_temp-itemp)
  d_ac_heat = min(ac_heatcap*fp%bldheight,ac_load)*ac_heat_on
end where
! cooling load
where ( itemp > ac_cool_temp )
  !ac_load = rm*(itemp-ac_cool_temp)
  ac_load = (rm+rf+we+ww+sl+im1+im2+infl)*(itemp-ac_cool_temp)
  d_ac_cool = min(ac_coolcap*fp%bldheight,ac_load)*ac_cool_on
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

iroomtemp = (rm*room%nodetemp(:,1)                    & ! room temperature
           + rf*roof%nodetemp(:,nl)                   & ! roof conduction
           + we*walle%nodetemp(:,nl)                  & ! wall conduction east
           + ww*wallw%nodetemp(:,nl)                  & ! wall conduction west
           + sl*slab%nodetemp(:,nl)                   & ! slab conduction
           + im1*intm%nodetemp(:,0)                   & ! mass conduction side 1
           + im2*intm%nodetemp(:,nl)                  & ! mass conduction side 2
           + infl*canyontemp_ave                      & ! weighted infiltration
           + d_ac_inside                              & ! ac flux (iterative)
           + d_intgains_bld                           & ! internal gains (explicit)
           )/(rm+rf+we+ww+sl+im1+im2+infl)
!iroomtemp = itemp + d_ac_inside/rm ! time split version


! calculate infiltration flux with dampening factor from previous timestep (reduce oscillation)
int_infilflux = infl*(canyontemp_ave - iroomtemp)
int_infilfg = (fp%sigmabld/(1.-fp%sigmabld))*int_infilflux

d_roomstor = (iroomtemp - room%nodetemp(:,1))*(aircp*int_airden*fp%bldheight)/ddt

! ---1.5: Final iteration flux and temperature
!if ( l==ncyits ) then
!
!  !! diagnose fluxes from implicit results
!  ggint_roof  = cvcoeff_roof*(roof%nodetemp(:,nl)-iroomtemp)
!  ggint_slab  = cvcoeff_slab*(slab%nodetemp(:,nl)-iroomtemp)
!  ggint_walle = cvcoeff_walle*(walle%nodetemp(:,nl)-iroomtemp)
!  ggint_wallw = cvcoeff_wallw*(wallw%nodetemp(:,nl)-iroomtemp)
!  ggext_intm  = -cvcoeff_intm1*(intm%nodetemp(:,0)-iroomtemp) ! negative flux leaving facet
!  ggint_intm  = cvcoeff_intm2*(intm%nodetemp(:,nl)-iroomtemp) ! positive flux leaving facet
!
!  ! explicit air temperature evolution
!  iroomtemp = room%nodetemp(:,1) + (ddt/(aircp*int_airden*fp%bldheight))    & ! previous temp
!            * (ggint_roof  + ggint_slab                                     & ! roof/slab fluxes
!            + (ggint_walle + ggint_wallw)*(fp%bldheight/fp%bldwidth)        & ! wall fluxes
!            + (ggint_intm  - ggext_intm)*real(fp%intmassn)                  & ! mass fluxes (w/ neg for ggext)
!            +  int_infilflux + d_intgains_bld + d_ac_inside)                  ! other fluxes
!  d_roomstor = (iroomtemp - room%nodetemp(:,1))*(aircp*int_airden*fp%bldheight)/ddt
!  
!end if    


! ---1.6: Calculate acflux into canyon and update canyon temperature-----

where ( d_ac_cool>0. )
  ! original UCLEM COP calculation (Lipson et al., 2018)
  ! ac_coeff = max(acfactor*(d_canyontemp-iroomtemp)/(iroomtemp+urbtemp), 1./ac_copmax) ! T&H2012 Eq. 10
  ! modified UCLEM COP calculation (based on "Coefficient of performance for finite speed heat pump" - Blanchard 1980)
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
                      d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac,                                      &
                      effvegc,effrdsn,ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,                   &
                      lwflux_walle_vegc,lwflux_wallw_vegc,ddt,                                       &
                      cnveg,fp,fp_wall,rdhyd,road,walle,wallw,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull,2), intent(out) :: evct
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat
real, dimension(ufull), intent(in) :: sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_rg,a_rho,a_rnd,a_snd,d_sigd
real, dimension(ufull), intent(in) :: ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,lwflux_walle_vegc,lwflux_wallw_vegc
real, dimension(ufull), intent(out) :: effvegc,effrdsn
real, dimension(ufull), intent(inout) :: d_canyontemp,d_canyonmix,d_netrad,d_tranc,d_evapc
real, dimension(ufull), intent(inout) :: d_cra,d_crr,d_crw,d_totdepth,d_vegdeltac
real, dimension(ufull) :: ff,f1,f2,f3,f4
real, dimension(ufull) :: snevap
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
effvegc=cnveg%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(cnveg%temp+urbtemp)**4)    &
                  +sbconst*fp_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4          &
                  +(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
rg_vegc=effvegc-lwflux_walle_vegc-lwflux_wallw_vegc
effrdsn=snowemiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(rdsntemp+urbtemp)**4)        &
                  +sbconst*fp_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4          &
                  +(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
rg_rdsn=effrdsn-lwflux_walle_rdsn-lwflux_wallw_rdsn

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
gardsn=(rdsntemp-road%nodetemp(:,0))/ldratio ! use road temperature to represent canyon bottom surface temperature
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
                     d_sigr,ddt,fp_roof,rfhyd,rfveg,roof,fp,ufull)

implicit none

integer, intent(in) :: ufull
integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn
real, dimension(ufull), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(inout) :: rg_roof,eg_roof,acond_roof
real, dimension(ufull), intent(in) :: sg_rfsn,sg_vegr,sg_roof
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr
real, dimension(ufull), intent(inout) :: d_rfsndelta,d_vegdeltar,d_roofdelta
real, dimension(ufull) :: lzomroof,lzohroof,qsatr,dts,dtt,cdroof,z_on_l,newval,ldratio
real, dimension(ufull) :: aa,dd,ee
real, dimension(ufull,2) :: oldval,evctx,evctveg
type(facetparams), intent(in) :: fp_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: roof
type(fparmdata), intent(in) :: fp

lzomroof=log(d_rfdzmin/zoroof)
lzohroof=2.3+lzomroof
call getqsat(qsatr,roof%nodetemp(:,0),d_sigr)
dts=roof%nodetemp(:,0) + (roof%nodetemp(:,0)+urbtemp)*0.61*d_roofdelta*qsatr
dtt=d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_roof,cdroof,z_on_l,lzohroof,lzomroof,d_rfdzmin,dts,dtt,a_umag,1)

! update green roof and snow temperature
rfveg%temp=d_tempr
rfsntemp  =max(min(roof%nodetemp(:,0),300.-urbtemp),100.-urbtemp)
rg_vegr = fp_roof%emiss*(a_rg-sbconst*(roof%nodetemp(:,0)+urbtemp)**4) ! 1st guess
rg_rfsn = fp_roof%emiss*(a_rg-sbconst*(roof%nodetemp(:,0)+urbtemp)**4) ! 1st guess
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
                d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,ufull)
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
                  d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,ufull)
    evctx=evctveg-evctx
    where ( abs(evctx(:,1))>tol .and. rfveg%sigma>0. )
      newval=rfveg%temp-alpha*evctveg(:,1)*(rfveg%temp-oldval(:,1))/evctx(:,1)
      oldval(:,1)=rfveg%temp
      rfveg%temp=newval
    end where
    where ( abs(evctx(:,2))>tol .and. d_rfsndelta>0. )
      newval=max(min(rfsntemp-alpha*evctveg(:,2)*(rfsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
      oldval(:,2)=rfsntemp
      rfsntemp=newval
    end where
  end do
end if
fg_vegr=sg_vegr+rg_vegr-eg_vegr
fg_rfsn=sg_rfsn+rg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)

! estimate roof latent heat flux (approx roof_skintemp with roof%nodetemp(:,1))
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
rg_roof=fp_roof%emiss*(a_rg-sbconst*(roof%nodetemp(:,0)+urbtemp)**4)

return
end subroutine solveroof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof and snow fluxes

subroutine roofflux(evct,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,   &
                    sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,   &
                    d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr,d_vegdeltar,               &
                    d_rfsndelta,ddt,fp,fp_roof,rfhyd,rfveg,roof,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull,2), intent(inout) :: evct
real, dimension(ufull), intent(in) :: rfsntemp,sg_vegr,sg_rfsn
real, dimension(ufull), intent(inout) :: rfsnmelt,garfsn
real, dimension(ufull), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_sigr
real, dimension(ufull), intent(inout) :: d_vegdeltar,d_rfsndelta
real, dimension(ufull) :: lzomvegr,lzohvegr,vwetfac,dts,dtt,z_on_l,ff,f1,f2,f3,f4,cdvegr
real, dimension(ufull) :: vegqsat,dumvegdelta,res,sndepth,snlambda,ldratio,lzosnow,rfsnqsat,cdrfsn
real, dimension(ufull) :: lzotdum, snevap
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
  garfsn=(rfsntemp-roof%nodetemp(:,0))/ldratio
  
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

icyc_traffic     = ((1.-real_p)*cyc_traf(int_p)+real_p*cyc_traf(int_p+1))*fp_weekdayload
icyc_basedemand  = ((1.-real_p)*cyc_base(int_p)+real_p*cyc_base(int_p+1))*fp_weekdayload
icyc_proportion  = ((1.-real_p)*cyc_prop(int_p)+real_p*cyc_prop(int_p+1))*fp_weekdayload
icyc_translation = ((1.-real_p)*cyc_tran(int_p)+real_p*cyc_tran(int_p+1))

return
end subroutine getdiurnal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version allows the eddy size to change with canyon orientation
! which requires a numerical solution to the integral

subroutine getincanwind(ueast,uwest,ufloor,a_udir,z0,fp,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: a,b,wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull), intent(in) :: a_udir
type(fparmdata), intent(in) :: fp

! rotate wind direction so that all cases are between 0 and pi
! walls are flipped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

h=fp%bldheight*fp%coeffbldheight
w=fp%bldheight/fp%hwratio

theta1=asin(min(w/(3.*h),1.))
wsuma=0.
wsumb=0.
fsum=0.  ! floor

! integrate jet on road, venting side (A)
a=0.
b=max(0.,wdir-pi+theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0,ufull)

! integrate jet on wall, venting side
a=max(0.,wdir-pi+theta1)
b=max(0.,wdir-theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,1,ufull)

! integrate jet on road, venting side (B)
a=max(0.,wdir-theta1)
b=wdir
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0,ufull)

! integrate jet on road, recirculation side (A)
a=wdir
b=min(pi,wdir+theta1)
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0,ufull)

! integrate jet on wall, recirculation side
a=min(pi,wdir+theta1)
b=min(pi,wdir+pi-theta1)
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,1,ufull)

! integrate jet on road, recirculation side (B)
a=min(pi,wdir+pi-theta1)
b=pi
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0,ufull)

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
end subroutine getincanwind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version fixes the eddy size to the canyon width which allows
! for an analytic solution to the integral

subroutine getincanwindb(ueast,uwest,ufloor,a_udir,z0,fp,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull) :: dufa,dura,duva,ntheta
real, dimension(ufull) :: dufb,durb,duvb
real, dimension(ufull), intent(in) :: a_udir
type(fparmdata), intent(in) :: fp

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

h=fp%bldheight*fp%coeffbldheight
w=fp%bldheight/fp%hwratio

theta1=acos(min(w/(3.*h),1.))

call winda(dufa,dura,duva,h,w,z0,ufull) ! jet on road
call windb(dufb,durb,duvb,h,w,z0,ufull) ! jet on wall
ntheta=2. ! i.e., int_0^pi sin(theta) dtheta = 2.)
where (wdir<theta1.or.wdir>pi-theta1) ! jet on wall
  wsuma=duvb*ntheta
  wsumb=durb*ntheta
  fsum=dufb*ntheta
elsewhere                                   ! jet on road
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
! integrate winds

subroutine integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,mode,ufull)

implicit none

integer, intent(in) :: ufull
integer, intent(in) :: mode
integer n
!integer, parameter :: ntot=45
integer, parameter :: ntot=1 ! simplified method
real, dimension(ufull), intent(in) :: a,b,h,w,wdir,z0
real, dimension(ufull), intent(inout) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta,dtheta,st,nw
real, dimension(ufull) :: duf,dur,duv

dtheta=(b-a)/real(ntot)
if (any(dtheta>0.)) then
  select case(mode)
    case(0) ! jet on road
      do n=1,ntot
        theta=dtheta*(real(n)-0.5)+a
        st=abs(sin(theta-wdir))
        nw=max(w/max(st,1.E-9),3.*h)
        call winda(duf,dur,duv,h,nw,z0,ufull)
        wsuma=wsuma+dur*st*dtheta
        wsumb=wsumb+duv*st*dtheta
        fsum=fsum+duf*st*dtheta
      end do
    case(1) ! jet on wall
      do n=1,ntot
        theta=dtheta*(real(n)-0.5)+a
        st=abs(sin(theta-wdir))
        nw=min(w/max(st,1.E-9),3.*h)
        call windb(duf,dur,duv,h,nw,z0,ufull)
        wsuma=wsuma+dur*st*dtheta
        wsumb=wsumb+duv*st*dtheta
        fsum=fsum+duf*st*dtheta
      end do
    case DEFAULT
      write(6,*) "ERROR: Unknown ateb.f90 integratewind mode ",mode
      stop
  end select
end if

return
end subroutine integratewind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on road

subroutine winda(uf,ur,uv,h,w,z0,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,w,z0
real, dimension(ufull) :: a,u0,cuven,zolog

a=0.15*max(1.,3.*h/(2.*w))
u0=exp(-0.9*sqrt(13./4.))

zolog=log(max(h,z0+0.2)/z0)
cuven=log(max(refheight*h,z0+0.2)/z0)/log(max(h,z0+0.2)/z0)
cuven=max(cuven*max(1.-3.*h/w,0.),(u0/a)*(h/w)*(1.-exp(max(-a*max(w/h-3.,0.),-40.))))
uf=(u0/a)*(h/w)*(1.-exp(-3.*a))+cuven
!uf=(u0/a)*(h/w)*(2.-exp(-a*3.)-exp(-a*(w/h-3.)))
ur=(u0/a)*exp(-a*3.)*(1.-exp(-a))
! MJT suggestion
cuven=1.-1./zolog
uv=(u0/a)*exp(max(-a*max(w/h-3.,0.),-40.))*(1.-exp(-a))
uv=max(cuven,uv)
!uv=(u0/a)*exp(-a*(w/h-3.))*(1.-exp(-a))

return
end subroutine winda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on wall

subroutine windb(uf,ur,uv,h,win,z0,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,win,z0
real, dimension(ufull) :: a,dh,u0,w
real, dimension(ufull) :: zolog,cuven

w=min(win,1.5*h)

a=0.15*max(1.,3.*h/(2.*w))
dh=max(2.*w/3.-h,0.)
u0=exp(-0.9*sqrt(13./4.)*dh/h)

zolog=log(max(h,z0+0.2)/z0)
! MJT suggestion (cuven is multiplied by dh to avoid divide by zero)
cuven=h-(h-dh)*log(max(h-dh,z0+0.2)/z0)/zolog-dh/zolog
! MJT cuven is back to the correct units of m/s
cuven=max(cuven/h,(u0/a)*(1.-exp(-a*dh/h)))

uf=(u0/a)*(h/w)*exp(-a*(1.-dh/h))*(1.-exp(-a*w/h))
ur=(u0/a)*exp(-a*(1.-dh/h+w/h))*(1.-exp(-a))
uv=(u0/a)*(1.-exp(-a*(1.-dh/h)))+cuven
!uv=(u0/a)*(2.-exp(-a*(1.-dh/h))-exp(-a*dh/h))

return
end subroutine windb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate wind speed at canyon top
subroutine gettopu(d_topu,a_umag,z_on_l,fp_bldheight,pd_cduv,pd_cndzmin,fp_hwratio,ufull)
      
implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(in) :: z_on_l
real, dimension(ufull) :: z0_on_l,bldheight
real, dimension(ufull) :: pm0,pm1,integralm
real, dimension(ufull) :: ustar,neutral
real, dimension(ufull), intent(out) :: d_topu
real, dimension(ufull), intent(in) :: a_umag
real, dimension(ufull), intent(in) :: fp_bldheight, fp_hwratio
real, dimension(ufull), intent(inout) :: pd_cduv, pd_cndzmin


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

subroutine scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,smixr,rdsntemp,zonet, &
                    cnveg,fp,pd,rdhyd,road,ufull)
      
implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(in) :: smixr,rdsntemp,zonet
real, dimension(ufull) :: cd,thetav,sthetav
real, dimension(ufull) :: thetavstar,z_on_l,z0_on_l
real, dimension(ufull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(ufull) :: ustar,qstar,z10_on_l
real, dimension(ufull) :: neutral,neutral10,pm10
real, dimension(ufull) :: integralm10,tts,tetp
real, dimension(ufull) :: tstar,lna
real, dimension(ufull) :: utop,ttop,qtop,wf,tsurf,qsurf,n
real, dimension(ufull), intent(in) :: a_mixr,a_umag,a_temp
real, dimension(ufull), intent(in) :: u_ts
real, dimension(ufull), intent(in) :: d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd
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
    tstar=vkar*(a_temp-u_ts)/integralh
    
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
    pd%tscrn = a_temp - tstar*integralh/vkar
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
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = fp%bldheight*(1.-refheight)*z_on_l/pd%cndzmin ! calculate at canyon top
    z10_on_l = max(z10-fp%bldheight*refheight,1.)*z_on_l/pd%cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
!write(*,*),pd%cndzmin
!write(*,*),fp%bldheight
!write(*,*),refheight
!write(*,*),z10
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
    tetp = a_temp - tstar*integralh/vkar
    qtop = a_mixr - qstar*integralh/vkar
    utop = a_umag - ustar*integralm/vkar

    where (fp%bldheight<=z10) ! above canyon
      pd%u10=max(a_umag-ustar*integralm10/vkar,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf = d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*road%nodetemp(:,0)+cnveg%sigma*cnveg%temp)
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
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*road%nodetemp(:,0)+cnveg%sigma*cnveg%temp)
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
    tstar=vkar*(a_temp-tsurf)/integralh
    
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
    pd%tscrn = a_temp-tstar*integralh/vkar
    pd%qscrn = a_mixr-qstar*integralh/vkar
    pd%uscrn = max(a_umag-ustar*integralm/vkar,0.)
    pd%u10   = max(a_umag-ustar*integralm10/vkar,0.)
    
end select
pd%qscrn       = max(pd%qscrn,1.E-4)
      
return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disables aTEB so subroutine calls have no effect

subroutine atebdisable(diag)

implicit none

integer, intent(in) :: diag

if ( diag>=1 ) write(6,*) "Disable aTEB"

ateb_active = .false.

return
end subroutine atebdisable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following subroutines are used for internal varying temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates net longwave radiation flux (flux_rg) at each surface
! longwave flux is temperature dependent, so this subroutine should be run at each timestep
subroutine internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle, &
                           skintemps,fp,intl,ufull)

implicit none
integer, intent(in) :: ufull
integer :: j
real, dimension(ufull,4), intent(in) :: skintemps  ! floor, wall, ceiling, wall temperature array
real, dimension(ufull,4) :: epsil     ! floor, wall, ceiling, wall emissivity array
real, dimension(ufull,4) :: radnet    ! net flux density on ith surface (+ve leaving)
real, dimension(ufull)   :: radtot    ! net leaving flux density (B) on ith surface
real, dimension(ufull), intent(out) :: rgint_slab,rgint_wallw,rgint_roof,rgint_walle
type(intldata), intent(in) :: intl
type(fparmdata), intent(in) :: fp

radnet = 0.
epsil  = 0.9

!epsil = reshape((/(fp_slab%emiss,fp_wall%emiss,fp_roof%emiss,fp_wall%emiss, & 
!                    i=1,ufull)/), (/ufull,4/))

do j = 2,4
  radnet(:,j) = epsil(:,j)/(1. -epsil(:,j))*((sbconst*(skintemps(:,j) +urbtemp)**4)  & 
               - sum(intl%psi(:,j,:)*(sbconst*(skintemps(:,:) +urbtemp)**4),dim=2))
end do

radnet(:,1) = -fp%bldheight*(radnet(:,2)+radnet(:,4))/fp%bldwidth - radnet(:,3)

! energy conservation check
radtot(:) = abs(fp%bldwidth(:)*(radnet(:,1)+radnet(:,3)) + fp%bldheight*(radnet(:,2)+radnet(:,4)))

do j = 1,ufull
  if ( radtot(j)>energytol ) write(6,*) "error: radiation energy non-closure: ", radtot(j)
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
                          cvcoeff_intm1,cvcoeff_intm2,roomtemp,roof,slab,intm,ufull)
implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw
real, dimension(ufull), intent(out) :: cvcoeff_intm1,cvcoeff_intm2,cvcoeff_slab
real, dimension(ufull), intent(in)  :: roomtemp
type(facetdata), intent(in) :: roof
type(facetdata), intent(in) :: slab
type(facetdata), intent(in) :: intm
real, parameter :: tdel = 0.1
real, dimension(ufull) :: x

select case(cvcoeffmeth)
  case(0) ! DOE Simple
    cvcoeff_walle = 3.067    ! vertical surface coefficient constant
    cvcoeff_wallw = 3.067    ! vertical surface coefficient constant
    where ( roof%nodetemp(:,nl)>=roomtemp+tdel )
      cvcoeff_roof(:)=0.948  ! reduced convection 
    elsewhere ( roof%nodetemp(:,nl)>roomtemp-tdel )
      x = (roof%nodetemp(:,nl)-roomtemp+tdel)/(2.*tdel)
      cvcoeff_roof(:)=(1.-x)*4.04 + x*0.948 ! blend
    elsewhere
      cvcoeff_roof(:)=4.040  ! enhanced convection
    end where    
    cvcoeff_intm1 = 3.067   ! vertical surface coefficient constant
    cvcoeff_intm2 = 3.067   ! vertical surface coefficient constant
    where ( slab%nodetemp(:,nl)<=roomtemp-tdel )
      cvcoeff_slab(:)=0.7   ! reduced convection
    elsewhere ( slab%nodetemp(:,nl)<roomtemp+tdel )
      x = (slab%nodetemp(:,nl)-roomtemp+tdel)/(2.*tdel)
      cvcoeff_slab(:) = (1.-x)*0.7 + x*4.04 ! blend
    elsewhere
      cvcoeff_slab(:)=4.040  ! enhanced convection
    end where
  case(1) ! dynamic, from international standard ISO6946:2007, Annex A
    cvcoeff_walle = 2.5    ! vertical surface coefficient constant
    cvcoeff_wallw = 2.5    ! vertical surface coefficient constant
    where ( roof%nodetemp(:,nl)>=roomtemp+tdel )
      cvcoeff_roof(:)=0.7  ! reduced convection (upper surface)
    elsewhere ( roof%nodetemp(:,nl)>roomtemp-tdel )
      x = (roof%nodetemp(:,nl)-roomtemp+tdel)/(2.*tdel)
      cvcoeff_roof(:) = (1.-x)*5.0 + x*0.7 ! blend
    elsewhere
      cvcoeff_roof(:)=5.0  ! enhanced convection (upper surface)  
    end where
    where ( intm%nodetemp(:,nl)>=roomtemp+tdel )
      cvcoeff_intm2 = 0.7+5.7   ! reduced convection (upper surface) + radiation @20 deg C
    elsewhere ( intm%nodetemp(:,nl)>roomtemp-tdel )
      x = (intm%nodetemp(:,nl)-roomtemp+tdel)/(2.*tdel)
      cvcoeff_intm2 = (1.-x)*(5.+5.7) + x*(0.7+5.7) ! blend
    elsewhere
      cvcoeff_intm2 = 5.0+5.7   ! reduced convection (upper surface) + radiation @20 deg C
    end where
    where ( intm%nodetemp(:,0)<=roomtemp-tdel )    
      cvcoeff_intm1 = 0.7+5.7   ! reduced convection (lower surface) + radiation @20 deg C
    elsewhere ( intm%nodetemp(:,0)<roomtemp+tdel )
      x = (intm%nodetemp(:,nl)-roomtemp+tdel)/(2.*tdel)
      cvcoeff_intm1 = (1.-x)*(0.7+5.7) + x*(5.+5.7) ! blend
    elsewhere
      cvcoeff_intm1 = 5.0+5.7   ! reduced convection (lower surface) + radiation @20 deg C
    end where
    where ( slab%nodetemp(:,nl)<=roomtemp-tdel )
      cvcoeff_slab(:)=0.7   ! reduced convection (lower surface)
    elsewhere ( slab%nodetemp(:,nl)<roomtemp+tdel )
      x = (slab%nodetemp(:,nl)-roomtemp+tdel)/(2.*tdel)
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
                            roof,walle,wallw,slab,ufull)
implicit none

integer, intent(in)                 :: ufull
real, dimension(ufull), intent(in)  :: d_canyontemp,iroomtemp
real, dimension(ufull), intent(out) :: d_openwindows
real, dimension(ufull)              :: xtemp, mrt
type(facetdata), intent(in) :: roof, walle, wallw, slab
type(fparmdata), intent(in) :: fp

! mean radiant temperature estimation
mrt = 0.5*(fp%bldheight/(fp%bldwidth+fp%bldheight)*(walle%nodetemp(:,nl) + wallw%nodetemp(:,nl))) & 
    + 0.5*(fp%bldwidth/ (fp%bldwidth+fp%bldheight)*(roof%nodetemp(:,nl) + slab%nodetemp(:,nl)))
! globe temperature approximation (average of mrt and air temperature) [Celsius]
xtemp = 0.5*(iroomtemp + mrt)

d_openwindows = 1./(1. + exp( 0.5*(fp%bldairtemp+ac_deltat+5.-xtemp) ))*1./(1. + exp( (d_canyontemp-iroomtemp) ))

return
end subroutine calc_openwindows

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
