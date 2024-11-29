module zenith_m

!  This module contains vectorisable versions of solargh and zenith
!  for the conformal cubic model. The calculation of hang has been 
!  moved from solargh to zenith so that solargh need be called only 
!  once per time stp, with zenith doing all the geographically varying
!  calculations.

!  Martin Dix 1999-03-12

   public :: zenith, solargh

   real, parameter, private :: pi     = 3.14159265358979, &
                               hpi    = 0.5 * pi,         &
                               tpi    = 2.0 * pi,         &
                               piph   = pi  + hpi,        &
                               radian = pi  / 180.0,      &
                               pid24  = pi / 24.0 

contains

   subroutine solargh(day,bpyear,r,dlt,alp,slag)
!
!  subroutine solar is called by march and computes the radius
!    vector, the declination and right ascension of the sun, the
!    apparent sun lag angle (related to the equation of time), the hour
!    angle of the sun at sunset, the mean cosine of the sun's zenith
!    angle and daylight fraction for n specified latitudes given the
!    day and fraction.  see notes of august 1989.
!
!   This routine comes from Bryant McAvaney via Josef Sytkus. It has an
!   additional argument not mentioned above
!      bpyear =       Years before present
!   This is used to set up the orbital parameters for paleoclimate runs.
!   Allowed values are 0, 6000 and 21000.
!
!   Martin Dix 23/3/94
!--------------
      implicit none
!
      real, intent(in)  :: day       ! day of year
      real, intent(in)  :: bpyear    ! Years before present
      real, intent(out) :: r         ! radius vector (distance to sun in a. u.)
      real, intent(out) :: dlt       ! declination of sun
      real, intent(out) :: alp       ! right ascension of sun
      real, intent(out) :: slag      ! apparent sun lag angle 
                                     ! (west of mean sun is plus)
      real :: ec, peril, oblqty
      real :: sind, ecsq, angin, sni, cyear, crate, sma, en, sunper
      real :: craten, cosphi, beta, peri, velngm, el, elt, rlam
!
!     set fundamental astronomical parameters
      if ( bpyear == 0.0 ) then
         ec  = 0.016724
         peril  = 102.04
         oblqty  = 23.446
      else if ( bpyear == 6000. ) then
         ec  = 0.018682
         peril  =   0.87
         oblqty  = 24.105
      else if ( bpyear == 21000. ) then
         ec  = 0.018994
         peril  = 114.42
         oblqty  = 22.949
      else
         print *,'bad year for astronomy' ,bpyear
         print *,'valid years are 0, 6000 and 21000'
         stop
      endif
      ecsq = ec**2
      angin  = radian * oblqty
      sni    = sin(angin)
!
!   tropical or calendar year in mean solar days is cyear
      cyear = 365.0
!
!   crate is mean daily angular motion of planet in degrees with
!   respect to equinox
      crate = 360.0 / cyear
!
!   semimajor axis of planetary orbit in astronomical units is sma.
!   one a. u. = 149,600,000 km.  for earth sma=1.
      sma = 1.0
!
!  The reference point for the astronomical calendar is the vernal
!  equinox, which is defined to occur at 0.0 hours on March 21. This
!  occurs 79 days after January 1, which is the reference point for
!  the model calendar.
      en = day - 79.0
!
!   longitude of sun at perihelion in degrees is sunper.  this is
!    exactly 180 degrees from that of earth at perihelion.
!
      sunper = peril - 180.0
      craten = crate * en
!
!  computation of longitude of mean sun follows the procedure
!  described in Berger (JAS, vol. 35, 2362-2367), using eq. 2
!  from annex 2 of Berger's paper "Precession, eccentricity, obliquity,
!  insolation and paleoclimates" from the NATO ASI "Long Term
!  Climatic Variations, Data and Modelling" (J. C. Duplessy, ed.).
!  The equation is a variant of that given on p.2365 of the JAS paper.
      cosphi = sqrt(1.0-ecsq)
      beta = (1.0-cosphi)/ec
      peri = radian*sunper
!
!  longitude of mean sun at vernal equinox is velngm
      velngm = -2.0*(beta*(1.0+cosphi)*sin(-peri)) +                &
                beta*beta*(0.5+cosphi)*sin(-2.0*peri) +             &
           beta*beta*beta*(1.0/3.0+cosphi)*sin(-3.0*peri)         
!

!  mean longitude of sun in radians is el
      el = velngm + craten*radian
!
!  constrain el to range of zero to 2 * pi
      el = el - tpi*ifix(el/tpi)
      if (el .lt. 0.0) el= el+tpi
      elt   = el
!
!  ecliptic longitude of sun in radians is rlam
!  Berger (JAS, vol. 35, p. 2365)
      rlam = el + 2.0*ec*sin(el-peri)                            &
                   + 5.0*ecsq*sin(2.0*(el-peri))/4.0             &
                   + 13.0*ec*ecsq*sin(3.0*(el-peri))/12.0
!
!   constrain rlam to range of zero to 2 * pi
      rlam = rlam - tpi*ifix(rlam/tpi)
      if (rlam .lt. 0.0) rlam = rlam + tpi
!
!   compute right ascension of sun (alp) in radians
      alp = atan(cos(angin)*tan(rlam))
!
!   following 6 statements determine proper quadrant for root of
!    equation for alp above and insure that the difference slag
!    below has an absolute magnitude of order pi/10.
!
      if (alp .lt. 0.0) alp = alp + tpi
!
      if ((elt-alp) .gt. piph) elt = elt - tpi
      if ((elt-alp) .lt.-piph) elt = elt + tpi
!
      if (abs(elt-alp) .gt. hpi) then
         if ((elt-alp) .gt. hpi) alp=alp+pi
         if ((elt-alp) .lt.-hpi) alp=alp-pi
      endif
!
!   compute declination of sun (dlt) in radians
      sind = sni*sin(rlam)
      dlt  = asin(sind)
!
!   compute radius vector (distance to sun in a.u.)
      r = (sma*(1.0-ecsq))/(1.0+ec*cos(rlam-peri))
!
!   compute sun lag angle in radians (720*slag/pi = equation
!    of time in minutes)
      slag = elt - alp

   end subroutine solargh

subroutine zenith(fjd,r,dlt,slag,xlat,xlon,dhr,npts,coszrs,frac)

! Version of zenith for the global CC model. Based on standard GCM version
! Revision 1.4  

! This version works with arrays containing any random set of latitiudes,
! longitudes and the corresponding hour angles.

!
!  Zenith computes effective mean cosine of zenith angle and daylight
!  fraction from latitude and parameters given by subroutine solar
!
!  All angles are in radians.
!
   implicit none
   integer, intent(in) :: npts
   real, intent(in) :: fjd, r, dlt, slag, dhr
   real, intent(in), dimension(npts) :: xlat, xlon
   real, intent(out), dimension(npts) :: coszrs, frac

   logical, dimension(npts) :: rise, set, done
   real :: rs, gha, arg
   real, dimension(npts) :: ha, ss, cc, hloc, hs, he, dele, delw, width, mid

!  ha is the hour angle of the sun at that latitude at sunrise/set [0:pi]
!  This has to be adjusted for the local time.

!  Over the period a-b to a+b the average zenith angle is
!  1/2b \int_{a-b}^{a+b} cos(z) dz = cos(a) sin(b) / b
!  so to calculate the average angle need the mid time and the half width.
!  Here b is the half-width, called width in the code here.

!  All times are in range [-pi:pi] with zero at noon.
!  Let the period be [start:end], rise time -ha, set time ha
!  There are 6 cases. The inequalities are tricky because have to allow for
!  the possibility that the averaging period might wrap around midnight.
!  1) Sun is up for the whole period
!     start > rise and end < set
!  2) Sun is down for whole period
!     start > set and end < rise
!  3) Sun rises period
!     start < rise and end > rise
!  4) Sun sets in period
!     start <  set and end > set
!  5) Sun rises and sets in period
!     start < rise and end > set
!  6) Sun sets and rises in period


!  First compute hour angle of sunset at all latitudes (calculation moved
!  from solargh).

!  = acos ( tan(alat) * tan(dlt) )
!  If |tan(alat) * tan(dlt)| > 1 then the sun doesn't set or rise
!  If it's positive then the sun is up all the time, h=pi
!  If negative it's dark and h=0.

   ss = sin(xlat) * sin(dlt)
   cc = cos(xlat) * cos(dlt)
   where ( abs(ss) > abs(cc) ) 
      ha = hpi * ( ss + abs(ss) ) / abs(ss)  ! pi if ss > 0, 0 otherwise
   elsewhere
      ha = pi - acos ( ss / cc ) 
   endwhere

   rs = 1.
   gha = fjd*tpi+slag+pi    ! Add pi because time is zero at local noon
   arg = dhr*pid24
   
   frac   = 0.0
   coszrs = 0.0

!  Local hour angle shifted by half of the averaging period
   hloc = gha+xlon+arg
!  Put hloc into the range -pi, pi
   hloc = mod(hloc,tpi)
   where ( hloc > pi ) 
      hloc = hloc-tpi
   endwhere
   hs = hloc-arg  ! Local hour angle at start of averaging period
   he = hloc+arg  ! Local hour angle at end of averaging period

!  Note that set = T means that the averaging period ends after sunset. 
!  It doesn't mean that the sun actually sets during the period, because the
!  start may also be after sunset. Similarly with rise.
   set  = he > ha
   rise = hs < -ha

!  Note that unless arg > pi both these conditions can not be true
   done = .false.
   where ( hs < -pi .and. .not. (rise.and.set) )
      dele = 0.5*max(he+ha,0.)
      delw = 0.5*max(ha-hs-tpi,0.)
      frac = (dele+delw)/arg
      done = .true.
   endwhere
   where ( he > pi .and. .not. (rise.and.set) ) 
      dele = 0.5*max(he+ha-tpi,0.)
      delw = 0.5*max(ha-hs,0.)
      frac = (dele+delw)/arg
      done = .true.
   endwhere
   where ( frac /= 0.0  )
      coszrs  = rs*(ss+cc*(cos(ha-dele)*sin(dele)+                &
                           cos(ha-delw)*sin(delw))/(dele+delw))
   endwhere

   ! Need with the ifort compiler and IEEE trapping?
   width = 0.
   mid = 0.
   where ( rise .and. set .and. .not. done )
      width = ha
      mid = 0.0
   endwhere
   where ( rise .and. .not. set .and. .not. done ) 
      width = 0.5 * ( he + ha )
      mid = -ha+width
   endwhere
   where ( .not. rise .and. set .and. .not. done ) 
      width = 0.5 * ( ha - hs )
      mid = ha - width
   endwhere
   where ( .not. rise .and. .not. set .and. .not. done ) 
      width = arg
      mid = hloc
   endwhere
   where ( .not. done .and. ha > 0. .and. width > 0 )
      frac = width / arg
      coszrs = rs*(ss+cc*cos(mid)*sin(width)/width)
   endwhere

   coszrs =  min ( 1.0, coszrs )
   frac   =  min ( 1.0, frac )

end subroutine zenith

end module zenith_m
