! ateb+cable offline wrapper
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
! This code is a wrapper for the ATEB (Thatcher and Hurley) and UCLEM (Lipson 
! et al. 2018) urban land surface models, running offline tiled with CABLE v1.4.
!
! OUTPUTS:
! - ateb_flx.dat:  100% ateb tile surface energy flux components
! - cable_flx.dat: 100% cable tile surface energy flux components
! - grid_flx.dat:  grid-averaged fluxes weighted using urbanfrac
! - ateb_bld.dat:  100% ateb/uclem tile anthropogenic energy variables
! - ateb_hyd.dat:  100% ateb hydrology variables
! - ateb_alma.dat: ateb alma outputs (50+ outputs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module atebcable_module

contains

!program atebcable
SUBROUTINE atebcable()
! ateb
use ateb
use zenith_m

! cable
!use cable_air_module
!use cable_albedo_module
!use cable_canopy_module
!use cable_carbon_module
!use cable_common_module
!use cable_data_module
!use cable_def_types_mod
!use cable_radiation_module
!use cable_roughness_module
!use cable_soil_snow_module

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DECLARATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! diagnostic variables
integer :: diag    = 0    ! diagnostic message mode (0=off, 1=basic, 2=detailed)
integer :: ierr    = 0    ! error code
integer :: mp      = 0
real :: cpu_start  = 0    ! start time
real :: cpu_stop   = 0    ! stop time

integer :: nstep = 0                ! cumulative timestep
integer :: nout               ! number of timesteps in output interval
integer :: ifin = 1                 ! number of horizontal grid points
integer :: i                        ! line looper
integer :: k                        ! layer looper
real, dimension(1) :: sigmu = 1.    ! proportion of urban tile (hardset to 1)
logical :: write_step =.false.      ! =.true. if writing to file this loop
logical :: first_step =.true.       ! for initialisation

! offline namelist variables
real :: dt = 1800.                  ! model time step                   [s]
real :: dt_in = 1800.               ! time step of input data           [s]
real :: dt_out = 1800.              ! time step of output data          [s]
real :: lat = -37.81                ! latitude of gridpoint             [degrees_north]
real :: lon = 144.97                ! longitude of gridpoint            [degrees_east]
real :: urbanfrac = 1.              ! fraction of urban                 [1]
real, dimension(1) :: azmin = 40.   ! first model level height          [m]
character(len=60)  :: forcingfile = 'alpha04.dat'           ! default file name of forcing data
character(len=60)  :: contact =''   ! default contact
integer :: forcingskiprows = 7      ! number of lines to skip in forcing file

real :: utc_offset = 10.            ! default for Melbourne, updated in ateb.nml

! cable namelist variables
real :: soil_albedo = 0.2           ! cable soil albedo                 [1]
integer :: cveg1 = 2                ! cable 1st veg type for vegmap

! model parameters
real, parameter :: vmodmin = 0.01     ! minimum wind speed
real, parameter :: bpyear = 0         ! zenith_m mode (0=present day)
real, parameter :: sbconst = 5.67e-8  ! Stefan-Boltzmann constant
real, parameter :: lv = 2.501e6       ! Latent heat of vaporisation       [J kg^-1]
real, parameter :: rd = 287.04        ! Gas constant for dry air          [J kg^-1 K^-1]

! forcing variables
integer :: iyr     = 0    ! year
integer :: iday    = 0    ! day of year
integer :: itime   = 0    ! hour and minute of day
integer :: itime_out     ! hour and minute at start of averaging period

! atmospheric variables at measurement height
real :: dectime    = 0.   ! day of year and decimal time of day (local)
real :: dectime_out       ! day of year at start of averaging period
real :: swdown     = 0.   ! shortwave flux down           [W m-2]
real :: lwdown     = 0.   ! longwave flux down            [W m-2]
real :: windv      = 0.   ! wind v component              [m s-1]
real :: windu      = 0.   ! wind u component              [m s-1]
real :: psurf      = 0.   ! air pressure                  [Pa]
real :: tair       = 0.   ! air temperature               [K]
real :: qair       = 0.   ! specific humidity             [1]
real :: rainf      = 0.   ! rainfall rate                 [kg s-1]
real :: snowf      = 0.   ! snowfall rate                 [kg s-1]
character(len=80)  :: dum ! to discard

real :: swnet
real :: lwnet
real :: wind

! ateb input variables (atebcalc)
real, dimension(1) :: fg        ! sensible heat flux                        [W m-2]
real, dimension(1) :: eg        ! latent heat flux                          [W m-2]
real, dimension(1) :: tss       ! radiative surface temperature             [K]
real, dimension(1) :: wetfac    ! soil wetness (fraction of soil moisture above wilting point relative to field capacity) [1]
real, dimension(1) :: newrunoff ! surface runnoff                           [kg m-2]
!                     dt        ! model time step                           [s]
!                     azmin     ! first model level height                  [m]
real, dimension(1) :: sg        ! incoming short wave radiation             [W m-2]
real, dimension(1) :: rg        ! incoming long wave radiation              [W m-2]
real, dimension(1) :: rnd       ! incoming rainfall rate                    [kg m-2 s-1]
real, dimension(1) :: evspsbl   ! evapotranspiration
real, dimension(1) :: sbl       ! sublimation
real, dimension(1) :: snd       ! incoming snowfall rate                    [kg m-2 s-1]
real, dimension(1) :: rhoair    ! atmospheric density at first model level  [kg m-3]
real, dimension(1) :: t         ! air temperature at first model level      [K]
real, dimension(1) :: qg        ! atmospheric mixing ratio                  [kg kg-1]
real, dimension(1) :: ps        ! pressure at first model level             [Pa]
real, dimension(1) :: uzon      ! wind U component (zonal)                  [m s-1]
real, dimension(1) :: vmer      ! wind V component (meridional)             [m s-1]
!                     vmodmin   ! minimum wind speed                        [m s-1]

real, dimension(1) :: alb       ! ateb albedo                               [1]
real, dimension(1) :: coszro2   ! cos(zenith angle of sun)                  [-]
real, dimension(1) :: rlongg    ! longitude in radians                      [rad]
real, dimension(1) :: rlatt     ! latitude in radians                       [rad]
real, dimension(1) :: taudar2   ! daylight fraction for spitter             [-]

real, dimension(1) :: tscrn,qscrn,uscrn,u10

! instantaneous output variables for ateb and cable
real :: swup    = 0.     ! upwards short wave radiation    [W m-2]
real :: lwup    = 0.     ! upwards long wave radiation     [W m-2]
real :: fgup    = 0.     ! upwards sensible heat flux      [W m-2]
real :: egup    = 0.     ! upwards sensible heat flux      [W m-2]

! ateb output variables
real ateb_swdown,ateb_lwdown,ateb_swup,ateb_lwup,ateb_fgup,ateb_egup
! cable output variables
real cable_swdown,cable_lwdown,cable_swup,cable_lwup,cable_fgup,cable_egup,cable_storage
! grid weighted output variables
real grid_swdown,grid_lwdown,grid_swup,grid_lwup,grid_fgup,grid_egup,grid_stor,grid_anth

! further ateb variables
real, dimension(3) :: roomtemp  ! room air temperature for three conditions            [K]
real, dimension(1) :: bldheat   ! anthropogenic heat flux from heating buildings       [W m-2]
real, dimension(1) :: bldcool   ! anthropogenic heat flux from cooling buildings       [W m-2]
real, dimension(1) :: intgain   ! anthropogenic heat flux from appliances in buildings [W m-2]
real, dimension(1) :: traffic   ! anthropogenic heat flux from traffic                 [W m-2]
real, dimension(1) :: anthrop   ! net anthropogenic heat flux                          [W m-2]
real, dimension(1) :: storage   ! storage heat flux in canopy                          [W m-2]
real, dimension(1) :: allbldtemp   ! building air temperature                             [K]
! real, dimension(1) :: scrntemp  ! screen level air temperature                         [K]
real ateb_bldheat,ateb_bldcool,ateb_intgain,ateb_traffic,ateb_anthrop
real ateb_storage

real, dimension(4) :: roadtemp_acc, roadtemp_out    ! road layer temperature [K]
real, dimension(4) :: rooftemp_acc, rooftemp_out    ! road layer temperature [K]
real, dimension(4) :: walletemp_acc, walletemp_out  ! road layer temperature [K]
real, dimension(4) :: wallwtemp_acc, wallwtemp_out  ! road layer temperature [K]
real, dimension(4) :: slabtemp_acc, slabtemp_out  ! internal slab layer temperature [K]
real, dimension(4) :: intmtemp_acc, intmtemp_out  ! internal mass layer temperature [K]

real :: surfstor_acc, surfstor_out         ! SurfStor       Surface water storage (excl. soil, snow, intercept) [kg/m2]
real :: swe_acc, swe_out                   ! SWE            Snow water equivalent                               [kg/m2]
real :: qsm_acc, qsm_out                   ! Qsm           Snowmelt (solid to liquid)  [kg/m2/s]
real :: soilmoist_acc, soilmoist_out       ! SoilMoist      Average layer soil moisture (1 layer)               [kg/m2]
real :: rootmoist_acc, rootmoist_out       ! RootMoist      Root zone soil moisture                             [kg/m2]
real :: acondveg_acc, acondveg_out         ! ACond          Aerodynamic conductance vegetation canopy           [m s-1]
real :: snowt_acc, snowt_out               ! SnowT          Snow surface temperature                            [K]
real :: vegt_acc, vegt_out                 ! VegT           Vegetation canopy temperature                       [K]
real :: radt_acc, radt_out                 ! RadT           Surface radiative temperature                       [K]

real :: roadsurft_out                      ! RoadSurfT      Road surface temperature (skin)                     [K]
real :: roofsurft_out                      ! RoadSurfT      Road surface temperature (skin)                     [K]
real :: wallsurft_out                      ! RoadSurfT      Road surface temperature (skin)                     [K]
real :: taircanyon_acc, taircanyon_out     ! TairCanyon     Air temperature in street canyon (bulk)             [K]
real :: tairsurf_acc, tairsurf_out         ! TairSurf       Near surface air temperature (2m)                   [K]
real :: tairbld_acc, tairbld_out           ! TairBuilding   Air temperature in buildings (bulk)                 [K]
real :: albedo_acc, albedo_out             ! Albedo         Surface albedo                                      [1]
real :: salbedo_acc, salbedo_out           ! SAlbedo        Snow albedo                                         [1]
real :: calbedo_acc, calbedo_out           ! CAlbedo        Vegetation canopy albedo                            [1]
real :: lai_acc, lai_out                   ! LAI            Leaf area index                                     [1]
real :: snowfrac_acc, snowfrac_out         ! SnowFrac       Snow covered fraction                               [1]
real :: soilwet_acc, soilwet_out           ! SoilWet        Total soil wetness                                  [1]
real :: qirrig_acc, qirrig_out             ! Qirrig         Anthropogenic water flux from irrigation (increase) [kg/m2/s]
real :: tranveg_acc, tranveg_out           ! TVeg           Vegetation transpiration                            [kg/m2/s]
real :: evap_acc, evap_out                 ! Evap           Total evapotranspiration (upward)                   [kg/m2/s]
real :: runoff_acc,runoff_out              ! Qs             Surface runoff (out of gridcell)                    [kg/m2/s]
real :: soilwaterfrac_acc,soilwaterfrac_out ! SMLiqFrac     Average layer fraction of liquid moisture           [1]

! Change in the simulated vertically integrated soil water volume, averaged over a grid cell, accumulated over the sampling time interval.
real :: delsoilmoist_acc, delsoilmoist_out ! DelSoilMoist   Change in soil moisture                             [kg/m2]  
! DelSWE change in snow water equivalent on ground and conopy voer output timestep  [kg/m2]
real :: delswe_acc, delswe_out             ! DelSWE         Change in snow water equivalent (increase)          [kg/m2]
real :: delintercept_acc,delintercept_out  ! DelIntercept   Change in liquid water storage in the canopy        [kg/m2]

! alma output variables
real, dimension(1) :: cduv_work   = 0. 
real, dimension(1) :: cdtq_work   = 0.

real :: qtau_acc,qtau_out     = 0.  ! Qtau        Momentum flux (downward)                [N m-2]
real :: tmpvar

! cable
!real, dimension(1) :: fbeam
!type (air_type),            save :: air
!type (bgc_pool_type),       save :: bgc
!type (met_type),            save :: met
!type (balances_type),       save :: bal
!type (radiation_type),      save :: rad
!type (roughness_type),      save :: rough
!type (soil_parameter_type), save :: soil
!type (soil_snow_type),      save :: ssnow
!type (sum_flux_type),       save :: sum_flux
!type (veg_parameter_type),  save :: veg
!type (canopy_type),         save :: canopy
!type (physical_constants),  save :: c

! for zenith_m
real :: fjd           ! day of year in universal time
real :: slag          ! apparent sun lag angle (west positive)
real :: dlt           ! declination of sun
real :: dhr           ! timestep in hours
real :: r1            ! radius vector (distance to sun in AU)
real :: alp           ! right ascension of sun

real, dimension(1) :: atmoserr=0          ! ateb energy closure residual in atmosphere
real, dimension(1) :: atmoserr_bias=0     ! ateb cumulative energy closure error in atmosphere
real, dimension(1) :: surferr=0           ! ateb energy closure residual in materials
real, dimension(1) :: surferr_bias=0      ! ateb cumulative energy closure error in materials

! for user defined vegetation                                                - MJT
integer, parameter :: mp_max = 5 ! maximum number of PFTs with mp <= mp_max  - MJT
integer, dimension(mp_max) :: veg_map ! PFT for each CABLE index             - MJT
character(len=1024) :: veg_param_file ! input file for veg PFTs              - MJT

! TO DELETE
real dayload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAMELISTS

namelist /offline/ dt,dt_in,dt_out,lat,lon,azmin,urbanfrac,utc_offset,forcingskiprows,forcingfile,contact
!namelist /cablenml/ soil_albedo,cveg1
namelist /atebnml/ intairtmeth

open(2,file='ateb.nml')
read(2,nml=offline)
!read(2,nml=cablenml)
close(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INITIALISATION

call cpu_time(cpu_start)

call error_check(urbanfrac,lat,lon)

! initialise urban model with sigmu=1. (urbanfrac weighted in this wrapper)
call atebinit(ifin,sigmu,diag)
!WRITE(*,*), 'working here'
mp=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONFIGURATION

rlongg(1)  = lon*3.1415927/180.        ! longitude
rlatt(1)   = lat*3.1415927/180.        ! latitude
dhr        = dt/3600.                  ! timestep in hours (for zenith_m)
nout = int(dt_out/dt)            ! number of timesteps in each output step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! input/output

write(6,'(A15,A50)')     'Forcing file: ', forcingfile
write(6,'(A15,F6.0,A4)') 'timestep (dt):', dt, ' [s]'
write(6,'(A15,F6.0,A4)') 'input dt:', dt_in, ' [s]'
write(6,'(A15,F6.0,A4)') 'output dt:', dt_out, ' [s]'
write(6,'(A15,F6.2,A4)') 'urban fraction:', urbanfrac, ' [1]'

! open files for input and output
open(1,file=forcingfile)
open(11,file="ateb_flx.dat")
open(12,file="ateb_bld.dat")
open(13,file="ateb_hyd.dat")
open(14,file="cable_flx.dat")
open(15,file="ateb_alma.dat")
open(21,file="grid_flx.dat")

if ( first_step ) then 
    ! header: ateb_flx.dat
    write(11,*) "Flux output from ATEB+UCLEM. Contact: ",contact
    write(11,'(A5,A10,A6,8A10)') 'YYYY','DOY.0000','HHMM','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2'
    write(11,'(A5,A10,A6,8A10)') 'Year','Dectime','Time','SWdown','LWdown','SWup','LWup','QHup','QEup','dQS','Qanth'

    ! header: ateb_bld.dat
    write(12,*) "Building output from UCLEM. Contact: ",contact
    write(12,'(A5,A10,A6,8A10)') 'YYYY','DOY.0000','HHMM','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','K','K'
    write(12,'(A5,A10,A6,8A10)') 'Year','Dectime','Time','AtmosErr','SurfErr','BldHeat','BldCool','IntGains','Traffic',   & 
        'ABldTemp','ScrnTemp'

    ! header: ateb_hyd.dat
    write(13,*) "Hyd output from ATEB. Contact: ",contact
    write(13,'(A5,A10,A6,13A10)') 'YYYY','DOY.0000','HHMM','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s',   &
        'kg/m2/s','1','1','m/s','kg/m2','kg/m2','kg/m2'
    write(13,'(A5,A10,A6,13A10)') 'Year','Dectime','Time','Rainf','Snowf','TotEvap','TransVeg','SurfDrain','Snowmelt',    & 
        'Irrig','SoilWet','SoilWFrac','ACond','DelSMoist','DelSWE','SoilMoist'

    ! header: cable_flx.dat
    write(14,*) "Flux output from CABLE. Contact: ",contact
    write(14,'(A5,A10,A6,8A10)') 'YYYY','DOY.0000','HHMM','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2'
    write(14,'(A5,A10,A6,8A10)') 'Year','Dectime','Time','SWdown','LWdown','SWup','LWup','QHup','QEup','dQS','Qanth'

    ! header: grid_flx.dat
    write(21,*) "Weighted flux output from ATEB+UCLEM+CABLE. Contact: ",contact
    write(21,'(A5,A10,A6,8A10)') 'YYYY','DOY.0000','HHMM','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2'
    write(21,'(A5,A10,A6,8A10)') 'Year','Dectime','Time','SWdown','LWdown','SWup','LWup','QHup','QEup','dQS','Qanth'

    ! header: ateb_alma.dat
    write(15,*) "Alma output from ATEB. Contact: ",contact
    write(15,'(A5,A10,A6,8A12,7A12,6A12,6A12,6A12,11A12,8A10)')                &
        'YYYY','DOY.0000','HHMM',                                              & 
        'W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2','W/m2',               &
        'W/m2','W/m2','W/m2','K','m/s','kg/kg','Pa',                           &
        'kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s',           &
        'kg/m2/s','kg/m2','kg/m2','kg/m2','kg/m2','kg/m2',                     &
        'kg/m2','kg/m2','N/m2','m/s','K','K',                                  &
        'K','K','K','K','K','K','K','K','K','K','K',                           &
        '1','1','1','1','1','1','1','1'
        
    write(15,'(A5,A10,A6,8A12,7A12,6A12,6A12,6A12,11A12,8A10)')                &
        'Year','Dectime','Time',                                               & 
        'SWnet','LWnet','Qle','Qh','Qanth','Qstor','SWup','LWup',              & 
        'SWdown','LWdown','Qanth_Qh','Tair','Wind','Qair','PSurf',             &
        'Rainf','Snowf','Evap','Qs','Qsm','Qirrig',                            &
        'TVeg','DelSoilM','DelInter','DelSWE','SWE','SurfStor',                &
        'SoilMoist','RootMoist','Qtau','ACond','SnowT','VegT',                 &
        'RadT','RoofSurfT','WallSurfT','RoadSurfT','TairSurf','TairCan',       &
        'TairBld','SoilTemp1','SoilTemp2','SoilTemp3','SoilTemp4',             &
        'Albedo','SAlbedo','CAlbedo','UAlbedo','LAI','SnowFrac','SoilWet','SMLiqFrac'
end if

! discard header lines
do i=1,forcingskiprows
  read(1,*) dum
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fg        = 0.
eg        = 0.
tss       = 0.
wetfac    = 0.
newrunoff = 0.
sg        = 0.
rg        = 0.
rnd       = 0.
snd       = 0.
rhoair    = 1.
t         = 0.
qg        = 0.
ps        = 0.
uzon      = 0.
vmer      = 0.

ateb_swdown   = 0.
ateb_lwdown   = 0.
ateb_swup     = 0.
ateb_lwup     = 0.
ateb_fgup     = 0.
ateb_egup     = 0.
ateb_storage  = 0.
ateb_anthrop  = 0.

cable_swdown  = 0.
cable_lwdown  = 0.
cable_swup    = 0.
cable_lwup    = 0.
cable_fgup    = 0.
cable_egup    = 0.
cable_storage = 0.

ateb_bldheat  = 0. 
ateb_bldcool  = 0.
ateb_intgain  = 0. 
ateb_traffic  = 0. 
ateb_anthrop  = 0.

roomtemp = 273.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main loop
if ( diag>=1 ) write(6,*) "Beginning main loop"

do ! exit when no more lines to read
    if ( mod( nstep+1,nout ) == 0 ) write_step=.true.

    if ( mod( nstep*int(dt), int(dt_in) ) == 0 ) then ! input step
        if (forcingfile=='alpha04.dat') then
            read(1,*,iostat=ierr) iyr,dectime,iday,itime,swdown,lwdown,windv,windu,psurf,tair,qair,rainf
            snowf = 0.
        else
            read(1,*,iostat=ierr) iyr,dectime,iday,itime,swdown,lwdown,windv,windu,psurf,tair,qair,rainf,snowf

        end if

        if ( ierr /= 0 ) exit ! exits if no more lines to read
        if ( diag>=1 ) write(6,*) "<-- input: nstep, dectime, itime ", nstep, dectime, itime

    else ! not input step (keep forcing constant, interpolat time)

        itime = itime + int(dt/60.)
        dectime = dectime + dt/(24.*60.*60)

        if ( diag>=1 ) write(6,*) "--- between: nstep, dectime, itime ", nstep, dectime, itime
    end if

    if ( mod( nstep,nout ) == 0 ) then ! start of output time averaging interval

        itime_out = itime
        dectime_out = dectime

    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set forcing

    ! set-up solar angles for zenith_m (radiation)
    ! fjd=dectime-rlongg(1)*0.5/3.1415927
    ! fjd = dectime - lon/360.

    ! replace fjd estimate with known utc offset
    fjd=dectime-(utc_offset/24.)

    call solargh( fjd,    & ! (in) decimal day of year
                  bpyear, & ! (in) year before present 
                  r1,     & ! (out) radius vector (distance to sun in a.u.)
                  dlt,    & ! (out) declination of sun
                  alp,    & ! (out) right ascension of sun
                  slag)     ! (out) apparent sun lag angle (west of mean sun is plus)

    call zenith( fjd,     & ! (in) decimal day of year
                 r1,      & ! (in) from solargh
                 dlt,     & ! (in) from solargh
                 slag,    & ! (in) from solargh
                 rlatt,   & ! (in) lat in radians 
                 rlongg,  & ! (in) lon in radians
                 dhr,     & ! (in) timestep (in hours)
                 1,       & ! (in) number of points
                 coszro2, & ! (out) cosine of zenith angle
                 taudar2)   ! (frac) daylight fraction

    call atebccangle( 1,        & ! (in) grid start integer
                      mp,       & ! (in) grid finish integer
                      coszro2,  & ! (in) from zenith
                      rlongg,   & ! (in) lon in radians
                      rlatt,    & ! (in) lat in radians
                      fjd,      & ! (in) decimal day of year
                      slag,     & ! (in) from solargh
                      dt,       & ! (in) timestep
                      sin(dlt))   ! (in) sin(declination) from solargh

    ! Update urban prognostic variables
    t(1)=tair     ! air temperature
    qg(1)=qair    ! mixing ratio of water vapour
    ps(1)=psurf   ! surface pressure
    sg(1)=swdown  ! downwelling shortwave
    rg(1)=lwdown  ! downwelling longwave
    rnd(1)=rainf  ! rainfall
    snd(1)=snowf  ! snowfall
    ! rhoair(1)=1.  ! air density
    uzon(1)=windu ! wind U
    vmer(1)=windv ! wind V

    wind = sqrt(windu**2 + windv**2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! run ateb
    if ( urbanfrac > 0. ) then

        ! use spitter 1986 to estimate diffuse/direct components of sg
        call atebspitter(1,1,fjd,sg,coszro2,diag)

        ! ! set direct SW beam at 100%, diffuse at 0%
        ! fbeam(1) = 1.
        ! call atebfbeam(1,1,fbeam,diag)

        ! if weekend reduce dayload for energy use (not used)
        call get_dayload(iyr,dectime,dayload)

        ! Extract albedo
        call atebalb1(1,1,alb,diag,raw=.true.,split=0)
        !WRITE(*,*), 'working here'
        ! ateb urban model main call
        call atebcalc(fg,eg,tss,wetfac,newrunoff,evspsbl,sbl,dt,azmin,  &
                    sg,rg,rnd,snd,rhoair,t,qg,ps,uzon,vmer,vmodmin,diag)
        call atebenergy(storage,"storage",diag)
        ! roomtemp(1)=fully conditioned, roomtemp(2)=either heating or cooling, roomtemp(3)=no heating/cooling
        call atebsaved(roomtemp(1),"roomtemp",1,0)
        if ( intairtmeth == 2 ) then
            call atebsaved(roomtemp(2),"roomtemp",2,0)
            call atebsaved(roomtemp(3),"roomtemp",3,0)
        else
            roomtemp(2) = -999.
            roomtemp(3) = -999.
        end if


        ! store timpestep flux data                  
        fgup=fg(1)
        egup=eg(1)
        lwup=sbconst*tss(1)**4
        swup=swdown*alb(1)

        ! store energy flux data
        call energyrecord(atmoserr,atmoserr_bias,surferr,surferr_bias,bldheat,bldcool,intgain,traffic,allbldtemp)
        ! screen level diagnostics
        call atebscrnout(tscrn,qscrn,uscrn,u10,diag,raw=.true.)

        ! net anthropogenic heat flux (missing industry)
        anthrop = bldheat + bldcool + intgain + traffic

        ! calculate storage as residual for first step
        if ( nstep == 0 ) storage(1) = swdown + lwdown + anthrop(1) - swup - lwup - fgup - egup

        ! store ateb output fluxes cumulatively
        ateb_swdown = ateb_swdown + swdown
        ateb_lwdown = ateb_lwdown + lwdown
        ! re-weight to make 100% weighted for ateb output
        ateb_swup = ateb_swup + swup
        ateb_lwup = ateb_lwup + lwup
        ateb_fgup = ateb_fgup + fgup
        ateb_egup = ateb_egup + egup
        ateb_storage = ateb_storage + storage(1)

        ateb_bldheat  = ateb_bldheat + bldheat(1)
        ateb_bldcool  = ateb_bldcool + bldcool(1)
        ateb_intgain  = ateb_intgain + intgain(1)
        ateb_traffic  = ateb_traffic + traffic(1)
        ateb_anthrop  = ateb_anthrop + anthrop(1)

        ! (ps=surface pressure in Pa, rdry=287.04 and tss is surface temperature (estimate of air temperature at surface)).
        rhoair = ps/(rd*tair)
        call atebcd(cduv_work,cdtq_work,diag,raw=.true.)
        tmpvar = rhoair(1)*cduv_work(1)*wind**2           ! (windspeed is the wind speed at the atmosphere level (e.g., 40m))
        call set_outvar(tmpvar,qtau_acc,qtau_out,nstep,nout,write_step,diag)

        call set_atebhydro('snowmelt',0.,qsm_acc,qsm_out,nstep,nout,write_step,diag)
        call set_atebhydro('swe',0.,swe_acc,swe_out,nstep,nout,write_step,diag)

        call set_atebhydro('surfstor',0.,surfstor_acc,surfstor_out,nstep,nout,write_step,diag)
        call set_atebhydro('soilmoisture',0.,soilmoist_acc,soilmoist_out,nstep,nout,write_step,diag)
        call set_atebhydro('rootmoisture',0.,rootmoist_acc,rootmoist_out,nstep,nout,write_step,diag)
        call set_atebhydro('aeroconductionveg',0.,acondveg_acc,acondveg_out,nstep,nout,write_step,diag)
        call set_atebhydro('snowfrac',0.,snowfrac_acc,snowfrac_out,nstep,nout,write_step,diag)
        call set_atebhydro('soilwet',0.,soilwet_acc,soilwet_out,nstep,nout,write_step,diag)
        call set_atebhydro('irrigation',0.,qirrig_acc,qirrig_out,nstep,nout,write_step,diag)
        call set_atebhydro('transpirationveg',0.,tranveg_acc,tranveg_out,nstep,nout,write_step,diag)
        call set_atebhydro('surfrunoff',0.,runoff_acc,runoff_out,nstep,nout,write_step,diag)
        call set_atebhydro('soilwaterfrac',0.,soilwaterfrac_acc,soilwaterfrac_out,nstep,nout,write_step,diag)

        call set_atebhydro_accumulate('delsoilmoist',0.,delsoilmoist_acc,delsoilmoist_out,nstep,nout,write_step,diag)
        call set_atebhydro_accumulate('delswe',0.,delswe_acc,delswe_out,nstep,nout,write_step,diag)
        call set_atebhydro_accumulate('delintercept',0.,delintercept_acc,delintercept_out,nstep,nout,write_step,diag)
            
        call set_atebmisc('snowt',273.,snowt_acc,snowt_out,nstep,nout,write_step,diag)
        call set_atebmisc('vegt',273.,vegt_acc,vegt_out,nstep,nout,write_step,diag)
        call set_atebmisc('taircanyon',273.,taircanyon_acc,taircanyon_out,nstep,nout,write_step,diag)
        call set_atebmisc('tairbuilding',273.,tairbld_acc,tairbld_out,nstep,nout,write_step,diag)
        call set_atebmisc('salbedo',0.,salbedo_acc,salbedo_out,nstep,nout,write_step,diag)
        call set_atebmisc('calbedo',0.,calbedo_acc,calbedo_out,nstep,nout,write_step,diag)
        call set_atebmisc('urbanlai',0.,lai_acc,lai_out,nstep,nout,write_step,diag)

        call set_outvar(tss(1),radt_acc,radt_out,nstep,nout,write_step,diag)
        call set_outvar(tscrn(1),tairsurf_acc,tairsurf_out,nstep,nout,write_step,diag)
        call set_outvar(alb(1),albedo_acc,albedo_out,nstep,nout,write_step,diag)
        tmpvar = eg(1)/lv
        call set_outvar(tmpvar,evap_acc,evap_out,nstep,nout,write_step,diag)

        call set_atebavetemp('roadtemp',273.,roadtemp_acc,roadtemp_out,nstep,nout,write_step,diag)
        roadsurft_out = roadtemp_out(1)

        call set_atebavetemp('rooftemp',273.,rooftemp_acc,rooftemp_out,nstep,nout,write_step,diag)
        roofsurft_out = rooftemp_out(1)

        call set_atebavetemp('walletemp',273.,walletemp_acc,walletemp_out,nstep,nout,write_step,diag)
        call set_atebavetemp('wallwtemp',273.,wallwtemp_acc,wallwtemp_out,nstep,nout,write_step,diag)
        wallsurft_out = 0.5 * (walletemp_out(1) + wallwtemp_out(1))

    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! output

    ! timestep for write out
    if ( write_step ) then
        if ( diag>=1 ) write(6,*) "--> output: nstep, dectime, itime, dectime_out, itime_out", & 
            nstep, dectime, itime, dectime_out, itime_out

        ! ateb energy flux averages
        ateb_swup = ateb_swup/real(nout)
        ateb_lwup = ateb_lwup/real(nout)
        ateb_fgup = ateb_fgup/real(nout)
        ateb_egup = ateb_egup/real(nout)
        ateb_swdown = ateb_swdown/real(nout)
        ateb_lwdown = ateb_lwdown/real(nout)
        ateb_storage = ateb_storage/real(nout)

        ! cable energy flux averages
        cable_swdown = cable_swdown/real(nout)
        cable_lwdown = cable_lwdown/real(nout)
        cable_swup = cable_swup/real(nout)
        cable_lwup = cable_lwup/real(nout)
        cable_fgup = cable_fgup/real(nout)
        cable_egup = cable_egup/real(nout)
        ! cable heat storage flux by residual (assuming zero advection)
        cable_storage = cable_swdown + cable_lwdown - cable_swup - cable_lwup - cable_fgup - cable_egup

        ! anthropogenic variables
        ateb_bldheat = ateb_bldheat/real(nout)
        ateb_bldcool = ateb_bldcool/real(nout)
        ateb_intgain = ateb_intgain/real(nout)
        ateb_traffic = ateb_traffic/real(nout)
        ateb_anthrop = ateb_anthrop/real(nout)

        ! ateb fluxes out
        write(11,'(I5,F10.4,I6,8F10.3)') iyr,dectime_out,itime_out,  &
            ateb_swdown,ateb_lwdown,ateb_swup,ateb_lwup,ateb_fgup,ateb_egup,ateb_storage,ateb_anthrop

        ! uclem building out
        write(12,'(I5,F10.4,I6,2F10.5,6F10.2)') iyr,dectime_out,itime_out, &
            atmoserr,surferr,ateb_bldheat,ateb_bldcool,ateb_intgain,ateb_traffic,tairbld_out,tairsurf_out

        ! ateb hydro out
        write(13,'(I5,F10.4,I6,13F10.5)') iyr,dectime_out,itime_out, &
            rainf,snowf,evap_out,tranveg_out,runoff_out,qsm_out,qirrig_out, &
            soilwet_out,soilwaterfrac_out,acondveg_out,delsoilmoist_out,delswe_out,soilmoist_out

        ! cable fluxes out
        write(14,'(I5,F10.4,I6,8F10.3)') iyr,dectime_out,itime_out, &
            cable_swdown,cable_lwdown,cable_swup,cable_lwup,cable_fgup,cable_egup,cable_storage,0.

        swnet     = ateb_swdown - ateb_swup
        lwnet     = ateb_lwdown - ateb_lwup

        ! alma fluxes out for ateb
        write(15,'(I5,F10.4,I6,8F12.4,5F12.4,1F12.6,1F12.2,6E12.4,6E12.4,6F12.4,11F12.4,8F10.4)')   &
            iyr,dectime_out,itime_out,                                                              & 
            swnet,lwnet,ateb_egup,ateb_fgup,ateb_anthrop,ateb_storage,ateb_swup,ateb_lwup,          & 
            ateb_swdown,ateb_lwdown,ateb_anthrop,tair,wind,qair,psurf,                              &
            rainf,snowf,evap_out,runoff_out,qsm_out,qirrig_out,                                     &
            tranveg_out,delsoilmoist_out,delintercept_out,delswe_out,swe_out,surfstor_out,          &
            soilmoist_out,rootmoist_out,qtau_out,acondveg_out,snowt_out,vegt_out,                   &
            radt_out,roofsurft_out,wallsurft_out,roadsurft_out,tairsurf_out,taircanyon_out,         &
            tairbld_out, roadtemp_out(1),roadtemp_out(2),roadtemp_out(3),roadtemp_out(4),           &
            albedo_out,salbedo_out,calbedo_out,albedo_out,lai_out,snowfrac_out,soilwet_out,soilwaterfrac_out

        ! grid weighted fluxes
        grid_swdown = ateb_swdown*urbanfrac + cable_swdown*(1.-urbanfrac)
        grid_lwdown = ateb_lwdown*urbanfrac + cable_lwdown*(1.-urbanfrac)
        grid_swup   = ateb_swup*urbanfrac + cable_swup*(1.-urbanfrac)
        grid_lwup   = ateb_lwup*urbanfrac + cable_lwup*(1.-urbanfrac)
        grid_fgup   = ateb_fgup*urbanfrac + cable_fgup*(1.-urbanfrac)
        grid_egup   = ateb_egup*urbanfrac + cable_egup*(1.-urbanfrac)
        grid_stor   = ateb_storage*urbanfrac + cable_storage*(1.-urbanfrac)
        grid_anth   = ateb_anthrop*urbanfrac + 0.

        write(21,'(I5,F10.4,I6,8F10.3)') iyr,dectime_out,itime_out, &
            grid_swdown,grid_lwdown,grid_swup,grid_lwup,grid_fgup,grid_egup,grid_stor,grid_anth

        ! zero out avg fluxes

        ! ateb
        ateb_swup = 0.  
        ateb_lwup = 0.  
        ateb_fgup = 0.  
        ateb_egup = 0.  
        ateb_swdown = 0.
        ateb_lwdown = 0.
        ateb_storage = 0.
        ateb_anthrop = 0.
        ! cable
        cable_lwdown   = 0.
        cable_swdown   = 0.
        cable_swup = 0.
        cable_lwup = 0.
        cable_fgup = 0.
        cable_egup = 0.
        cable_storage = 0.
        ! uclem
        ateb_bldheat = 0.
        ateb_bldcool = 0.
        ateb_intgain = 0.
        ateb_traffic = 0.

        write_step=.false.

    end if

    nstep = nstep + 1
    first_step=.false.
    
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

close(1)
close(11)
close(12)
close(13)
close(14)
close(15)
close(21)

call atebend(diag)

call cpu_time(cpu_stop)
print *, "Time:", cpu_stop - cpu_start, "seconds"

!end program
end SUBROUTINE atebcable









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_outvar(var,var_acc,var_out,nstep,nout,write_step,diag)

    implicit none

    ! sets accumulation and output values
    integer, intent(in)           :: nstep       ! timestep number
    integer, intent(in)           :: nout        ! timesteps in output period
    integer, intent(in)           :: diag        ! diagnostics (>1 for output)
    real, intent(in)              :: var         ! variable value at current timestep
    real, intent(inout)           :: var_acc     ! accumulated value
    real, intent(out)             :: var_out     ! output value
    logical, intent(in)           :: write_step

    ! zero out accumulation on first loop
    if (nstep==0) var_acc=0.

    ! accumulate variable
    var_acc = var_acc + var

    if (write_step) then
        var_out = var_acc/real(nout)
        var_acc = 0.
    end if

end subroutine set_outvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_atebhydro(mode,var_init,var_acc,var_out, & 
                         nstep,nout,write_step,diag)
    ! calls atebhydro and averages accumulation over output timestep

    use ateb
    implicit none

    character(len=*), intent(in)  :: mode        ! case to pass to ateb
    integer, intent(in)           :: nstep       ! timestep number
    integer, intent(in)           :: nout        ! timesteps in output period
    integer, intent(in)           :: diag        ! diagnostics (>1 for output)
    real, intent(in)              :: var_init    ! initial value to pass to ateb
    real, intent(inout)           :: var_acc     ! accumulated value
    real, intent(out)             :: var_out     ! output value
    logical, intent(in)           :: write_step
    ! local
    real, dimension(1)           :: var_ateb

    ! zero out accumulation on first loop
    if (nstep==0) var_acc=0.

    ! set ateb vector and pass to ateb
    var_ateb = var_init
    call atebhydro(var_ateb,mode,diag)

    ! accumulate variable
    var_acc = var_acc + var_ateb(1)

    if (write_step) then
        var_out = var_acc/real(nout)
        var_acc = 0.
    end if

end subroutine set_atebhydro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_atebhydro_accumulate(mode,var_init,var_acc,var_out, & 
                                    nstep,nout,write_step,diag)
    ! calls atebhydro and accumulates output timestep (not averaged)

    use ateb
    implicit none

    character(len=*), intent(in)  :: mode        ! case to pass to ateb
    integer, intent(in)           :: nstep       ! timestep number
    integer, intent(in)           :: nout        ! timesteps in output period
    integer, intent(in)           :: diag        ! diagnostics (>1 for output)
    real, intent(in)              :: var_init    ! initial value to pass to ateb
    real, intent(inout)           :: var_acc     ! accumulated value
    real, intent(out)             :: var_out     ! output value
    logical, intent(in)           :: write_step
    ! local
    real, dimension(1)           :: var_ateb

    ! zero out accumulation on first loop
    if (nstep==0) var_acc=0.

    ! set ateb vector and pass to ateb
    var_ateb = var_init
    call atebhydro(var_ateb,mode,diag)

    ! accumulate variable and zero out
    var_acc = var_acc + var_ateb(1)
    if (write_step) then
        var_out = var_acc
        var_acc = 0.
    end if

end subroutine set_atebhydro_accumulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_atebmisc(mode,var_init,var_acc,var_out, & 
                         nstep,nout,write_step,diag)
    ! calls atebmisc and averages accumulation over output timestep

    use ateb
    implicit none

    character(len=*), intent(in)  :: mode        ! case to pass to ateb
    integer, intent(in)           :: nstep       ! timestep number
    integer, intent(in)           :: nout        ! timesteps in output period
    integer, intent(in)           :: diag        ! diagnostics (>1 for output)
    real, intent(in)              :: var_init    ! initial value to pass to ateb
    real, intent(inout)           :: var_acc     ! accumulated value
    real, intent(out)             :: var_out     ! output value
    logical, intent(in)           :: write_step
    ! local
    real, dimension(1)           :: var_ateb

    ! zero out accumulation on first loop
    if (nstep==0) var_acc=0.

    ! set ateb vector and pass to ateb
    var_ateb = var_init
    call atebmisc(var_ateb,mode,diag)

    if (mode=='snowt') then
        where ( var_ateb < 100. )
            var_ateb = -999.
        end where
    end if

    ! accumulate variable
    var_acc = var_acc + var_ateb(1)

    if (write_step) then
        ! convert variable to output and zero out
        var_out = var_acc/real(nout)
        var_acc = 0.
    end if

end subroutine set_atebmisc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_atebavetemp(mode,var_init,var_acc,var_out, & 
                           nstep,nout,write_step,diag)
    ! calls atebavetemp and averages accumulation over output timestep

    use ateb
    implicit none

    character(len=*), intent(in)  :: mode        ! case to pass to ateb
    integer, intent(in)           :: nstep       ! timestep number
    integer, intent(in)           :: nout        ! timesteps in output period
    integer, intent(in)           :: diag        ! diagnostics (>1 for output)
    real, intent(in)              :: var_init    ! initial value to pass to ateb
    real, dimension(4), intent(inout)           :: var_acc     ! accumulated value
    real, dimension(4), intent(out)             :: var_out     ! output value
    logical, intent(in)           :: write_step
    ! local
    real, dimension(1)           :: var_ateb
    integer :: ii
    character(len=10) :: callstr

    ! zero out accumulation on first loop
    if (nstep==0) var_acc=0.

    do ii = 1,4
        ! create ateb call string for each layer
        var_ateb = var_init
        write(callstr,'(A9,I1.1)') mode, ii
        callstr = adjustl(callstr)
        call atebavetemp(var_ateb,callstr,diag)
        var_acc(ii) = var_acc(ii) + var_ateb(1)
    end do

    if (write_step) then
        ! convert variable to output and zero out
        var_out = var_acc/real(nout)
        var_acc = 0.
    end if

end subroutine set_atebavetemp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine error_check(urbanfrac,lat,lon)

real, intent(in) :: urbanfrac
real, intent(in) :: lat,lon

if ( (urbanfrac < 0.) .or. (urbanfrac > 1.0) ) then
    write (6,*) 'ERROR: urban fraction (urbanfrac) must be between 0 and 1, exiting',urbanfrac
    stop
end if

if ( (lat < -90.) .or. (lat  > 90.) ) then
    write (6,*) 'ERROR: latitude (lat) must be between -90 and 90 degrees, exiting',lat
    stop
end if

if ( (lon < -180.) .or. (lon > 180.) ) then
    write (6,*) 'ERROR: longitude (lon) must be between -90 and 90 degrees, exiting',lon
    stop
end if

end subroutine error_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! not used

subroutine get_dayload(iyr,dectime,dayload)

integer :: iyr                  ! year
real, intent(in) :: dectime         ! julian day
real, intent(out) :: dayload

! find weekends for day of week loading
dayload=1.
select case(iyr)
    ! 1 (first day of weekend) 
    case(2000,2005,2011)
      if ( (mod(floor(dectime),7)==1) .or. (mod(floor(dectime),7)==2) ) dayload=0.9
    ! 2 
    case(2010,2016)
      if ( (mod(floor(dectime),7)==2) .or. (mod(floor(dectime),7)==3) ) dayload=0.9
    ! 3 
    case(2,2004,2009,2015)
      if ( (mod(floor(dectime),7)==3) .or. (mod(floor(dectime),7)==4) ) dayload=0.9
    ! 4 
    case(0,1,2003,2014,2020)
      if ( (mod(floor(dectime),7)==4) .or. (mod(floor(dectime),7)==5) ) dayload=0.9
    ! 5 
    case(2002,2008,2013,2019)
      if ( (mod(floor(dectime),7)==5) .or. (mod(floor(dectime),7)==6) ) dayload=0.9
    ! 6 
    case(2001,2007,2018)
      if ( (mod(floor(dectime),7)==6) .or. (mod(floor(dectime),7)==7) ) dayload=0.9
    ! 7 
    case(2006,2012,2017)
      if ( (mod(floor(dectime),7)==7) .or. (mod(floor(dectime),7)==1) ) dayload=0.9
    ! none 
    case DEFAULT
      if ( (mod(floor(dectime),7)==1) .or. (mod(floor(dectime),7)==2) ) dayload=0.9
end select

end subroutine get_dayload
end module atebcable_module
