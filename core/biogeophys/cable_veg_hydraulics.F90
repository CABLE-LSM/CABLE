!=========================================================
! purpose: all routines for calculate hydraulic processes in vegetation
!
! contanct: zihan.lu@inrae.fr
!==========================================================
MODULE cable_veg_hydraulics_module
   USE cable_data_module, ONLY : icanopy_type, point2constants
   implicit none

   PUBLIC :: optimisation  ! does the actual optimization between supply and demand of transpiration
   PUBLIC :: zbrent, solve_psi_leaf, calc_psi_leaf, calc_psi_leaf_fails, &
      get_xylem_vulnerability, get_xylem_vulnerabilityx, calc_plc, integrate_vulnerability, &
      calc_transpiration, assim, QUADP, calc_michaelis_menten_constants, arrh, peaked_arrh
   PRIVATE

   TYPE(icanopy_type) :: C
CONTAINS
   SUBROUTINE optimisation(canopy, rad, vpd, press, tleaf, csx, &
      psi_soil, kcmax, kmax, PLCcrit, b_plant, &
      c_plant, N, vcmxt3, ejmxt3, rdx, vx3, cx1, an_canopy, &
      e_canopy, avg_kcan, gamma_star, p, i)
      ! zihanlu :ejmxt3 is not used, try to delete it

      ! Optimisation wrapper for the ProfitMax model. The Sperry model assumes that
      ! plants maximise the normalised (0-1) difference between relative gain and
      ! relative hydraulic risk at every timestep.
      !
      ! Implementation broadly follows Manon's python code.
      !
      ! References:
      !
      ! * Sperry JS, Venturas MD, Anderegg WRL, Mencuccini M, Mackay DS,  Wang Y,
      !   Love DM. 2017. Predicting stomatal responses to the environment from the
      !   optimization of photosynthetic gain and hydraulic cost. Plant, Cell &
      !   Environment 40: 816-830.
      !
      ! * Sabot, M.E.B., De Kauwe, M.G., Pitman, A.J., Medlyn, B.E., Verhoef, A.,
      !   Ukkola, A.M. and Abramowitz, G. (2020), Plant profit maximization improves
      !   predictions of European forest responses to drought. New Phytol, 226:
      !   1638-1655. doi:10.1111/nph.16376
      !
      ! Martin De Kauwe, 2020; Manon Sabot, 2022

      USE cable_def_types_mod
      USE cable_common_module

      IMPLICIT NONE

      TYPE (canopy_type), INTENT(INOUT) :: canopy
      TYPE (radiation_type), INTENT(INOUT) :: rad

      INTEGER, INTENT(IN) :: i, N
      REAL, INTENT(IN) :: cx1, kmax, PLCcrit, b_plant, c_plant, vpd, press, &
         tleaf, gamma_star
      REAL(r_2), INTENT(IN) :: psi_soil
      REAL, DIMENSION(mp), INTENT(IN) :: kcmax
      REAL, DIMENSION(mp,mf), INTENT(IN) :: vcmxt3, ejmxt3, rdx, vx3
      REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csx

      REAL, INTENT(INOUT) :: e_canopy
      REAL, DIMENSION(mf), INTENT(INOUT) :: an_canopy
      REAL(r_2), DIMENSION(N), INTENT(INOUT) :: p

      REAL, DIMENSION(mp), INTENT(OUT) :: avg_kcan

      ! local variables
      INTEGER :: j, k, idx
      LOGICAL, DIMENSION(N) ::  mask
      REAL :: p_crit, lower, upper, Cs, gsw, Vcmax, Jmax, Rd, Vj, Km
      REAL, DIMENSION(mf) :: fsun, apar, e_leaves, p_leaves, kc_leaves
      REAL, DIMENSION(N) :: Ci, Ac, Aj, A, an_leaf, gsc, Cx, kc, kc_iter, e_leaf, cost, &
         gain, profit

      REAL, PARAMETER :: &
         J_TO_MOL = 4.6E-6, & ! Convert from J to Mol for light
         MOL_TO_UMOL = 1E6, MMOL_2_MOL = 1E-3, MOL_TO_MMOL = 1E3

      ! Michaelis Menten coefficient, umol m-2 s-1
      Km = cx1 * MOL_TO_UMOL

      ! Xylem pressure beyond which there's full dessication, MPa
      p_crit = -b_plant * LOG(100. / (100. - PLCcrit)) ** (1. / c_plant)

      ! Loop over sunlit, shaded parts of the canopy and solve the carbon uptake
      ! and transpiration
      DO j=1, 2

         ! absorbed par for the sunlit or shaded leaf, umol m-2 -s-1
         apar(j) = rad%qcan(i,j,1) * J_TO_MOL * MOL_TO_UMOL

         ! If there is bugger all light, assume there are no fluxes
         IF (apar(j) <= 0.) THEN

            ! load into stores
            an_canopy(j) = -rdx(i,j) * MOL_TO_UMOL ! umol m-2 s-1
            e_leaves(j) = 0. ! mol H2O m-2 s-1

         ELSE

            fsun(j) = rad%fvlai(i,j) / canopy%vlaiw(i)

            ! CO2 concentration at the leaf surface, umol mol-1
            Cs = csx(i,j) * MOL_TO_UMOL

            ! Generate a sequence of Ci's that we will solve the optimisation
            ! model for, range btw >gamma_star and up to <Cs; umol mol-1
            DO k=1, INT(0.95 * N)

               Ci(k)  = gamma_star + FLOAT(k) * (Cs - gamma_star) / FLOAT(N - 1)

            END DO

            ! max rate of rubisco activity, scaled up to sunlit/shaded canopy
            Vcmax = vcmxt3(i,j) * MOL_TO_UMOL

            ! potential rate of electron transport, scaled up to sun/shade
            Jmax = ejmxt3(i,j) * MOL_TO_UMOL

            ! day respiration, scaled up to sunlit/shaded canopy
            Rd = rdx(i,j) * MOL_TO_UMOL

            ! Rate of electron transport (NB this is = J/4) so is a func of apar
            Vj = vx3(i,j) * MOL_TO_UMOL

            ! Calculate the sunlit/shaded A_leaf (i.e. scaled up), umol m-2 s-1
            Ac = assim(Ci, gamma_star, Vcmax, Km) ! umol m-2 s-1
            Aj = assim(Ci, gamma_star, Vj, 2. * gamma_star) ! umol m-2 s-1
            A = -QUADP(1. - 1E-4, Ac + Aj, Ac * Aj) ! umol m-2 s-1
            an_leaf = A - Rd ! Net photosynthesis, umol m-2 s-1

            ! Use an_leaf to infer gsc_sun/sha. NB. An is the scaled up values via
            ! scalex applied to Vcmax/Jmax
            gsc = MAX(0., A / (Cs - Ci)) ! mol CO2 m-2 s-1

            ! Infer E_sun/sha from gsc. NB. as we're iterating, Tleaf will change
            ! and so will leaf surface VPD, maintaining energy balance
            e_leaf = gsc * C%RGSWC / press * vpd ! mol H2O m-2 s-1

            ! Infer leaf water potential, MPa, having rescaled E from big-leaf to
            ! unit leaf first.
            ! N.B: solving with an accuracy of 0.1*kcrit on kc
            p = calc_psi_leaf(psi_soil, &
               e_leaf * MOL_TO_MMOL / rad%scalex(i,j), kmax, &
               b_plant, c_plant, 1E-3 * (100. - PLCcrit) * kmax, N)

            ! Alternative zbrent solver which doesn't work
            !p = calc_psi_leaf_fails(p_crit + 1E-2, p_sat - 1E-2, &
            !                        e_leaf * MOL_TO_MMOL / rad%scalex(i,j), kmax, &
            !                        b_plant, c_plant, 1E-6, N)

            ! Soil-plant hydraulic conductance at canopy xylem pressure,
            ! mmol m-2 s-1 MPa-1
            kc = kmax * get_xylem_vulnerabilityx(p, b_plant, c_plant)

            ! Ensure we don't check for profit in bad search space
            WHERE (p < psi_soil .AND. p > p_crit)

               mask = .TRUE.

            ELSEWHERE  ! not physically possible

               mask = .FALSE.

            END WHERE

            ! normalised gain (-)
            gain = an_leaf / MAXVAL(an_leaf, mask=mask)

            ! normalised cost (-)
            cost = (kcmax(i) - kc) / (kcmax(i) - kmax * (100. - PLCcrit) / 100.)

            ! profit
            profit = gain - cost

            ! Mask instances of very high/low profit. These result from our
            ! optimisation approximation, which doesn't fully ramp up the cost as
            ! would be expected in the full optimiality model. Here, we are
            ! ignoring profit above/below 98/2%, which we can fairly conclude is
            ! noise. For example, profit ~100% is typically due to gsc being
            ! bonkers because of an unrealistic Ci:Cs.
            WHERE (profit < 0.005 .OR. profit > 0.995)

               mask = .FALSE.

            END WHERE

            ! Locate maximum profit
            idx = MAXLOC(profit, 1, mask=mask)

            IF (idx > 1) THEN ! the optimisation conditions are satisfied

               an_canopy(j) = an_leaf(idx) ! umolC m-2 s-1
               e_leaves(j) = e_leaf(idx) ! molH2O m-2 s-1
               canopy%gswx(i,j) = MAX(1E-9, gsc(idx) * C%RGSWC) ! molH2O m-2 s-1
               p_leaves(j) = p(idx)
               kc_leaves(j) = kc(idx)

            ELSE ! the optimisation conditions are not satisfied

               an_canopy(j) = -rdx(i,j) * MOL_TO_UMOL ! umol m-2 s-1
               e_leaves(j) = 0. ! mol H2O m-2 s-1

               ! It's not clear what should happen to the leaf water potential if
               ! there is no light. Setting it to the root zone water potential
               ! on the basis that there is some evidence the leaf and soil water
               ! potential come back into equilibrium overnight.
               p_leaves(j) = psi_soil

               ! Unclear what should happen to the conductance too...
               kc_leaves(j) = kcmax(i)

            END IF

            !IF (vpd < 1.0) THEN
            !   print*, vpd, apar(j), profit(idx), e_leaves(j), e_leaf(idx)
            !END IF

         END IF

      END DO

      ! sunlit/shaded weighted components
      IF (apar(1) > 0. .AND. apar(2) > 0.) THEN

         canopy%psi_can(i) = (p_leaves(1) * fsun(1)) + (p_leaves(2) * fsun(2))
         avg_kcan(i) = SUM(fsun) / (fsun(1) / kc_leaves(1) + fsun(2) / &
            kc_leaves(2))

      ELSE IF (apar(1) > 0.) THEN

         canopy%psi_can(i) = p_leaves(1)
         avg_kcan(i) = kc_leaves(1)

      ELSE IF (apar(2) > 0.) THEN

         canopy%psi_can(i) = p_leaves(2)
         avg_kcan(i) = kc_leaves(2)

      ELSE ! nothing has happened, we're in the roots

         ! It's not clear what should happen to the leaf water potential if
         ! there is no light. Setting it to the root zone water potential on
         ! the
         ! basis that there is some evidence the leaf and soil water potential
         ! come back into equilibrium overnight.
         canopy%psi_can(i) = psi_soil

         ! Unclear what should happen to the conductance too...
         avg_kcan(i) = kcmax(i)

      END IF

      e_canopy = SUM(e_leaves) ! mol H2O m-2 s-1

   END SUBROUTINE optimisation
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION zbrent(func , x1 , x2 , tol, aval, bval, cval, dval, eval) &
      RESULT (xroot)

      ! This is a bisection routine. When ZBRENT is called, we provide a reference
      ! to a particular function and also two values which bound the arguments for
      ! the function of interest. ZBRENT finds a root of the function (i.e. the
      ! point where the function equals zero), that lies between the two bounds.
      ! This implementation of ZBRENT allows for parsing five additional arguments
      ! to the function being optimised, therefore the function being optimised
      ! needs to be written as a six-argument function, even if those arguments
      ! are not used (in which case their value will not matter).
      !
      ! For a full description of ZBRENT see Press et al. (1986).
      !
      ! Code adapted from SPA; Manon Sabot, 2023

      IMPLICIT NONE

      REAL, INTENT(IN) :: x1, x2, tol, aval, bval, cval, dval, eval

      ! Interfaces are the correct way to pass procedures as arguments
      INTERFACE

         REAL FUNCTION func(xval, aval, bval, cval, dval, eval)
            REAL, INTENT(IN) :: xval, aval, bval, cval, dval, eval
         END FUNCTION func

      END INTERFACE

      ! local variables
      INTEGER :: iter
      INTEGER, PARAMETER :: NITER = 100
      REAL :: a, b, c, d, e, p, q, fa, fb, fc, tol1, xm, xroot
      REAL, PARAMETER :: EPS = 3.e-8

      ! Calculations...
      a = x1
      b = x2
      fa = func(a, aval, bval, cval, dval, eval)
      fb = func(b, aval, bval, cval, dval, eval)

      ! Check that we haven't started with the root
      IF ( fa .eq. 0. ) THEN
         xroot = a

      ELSEIF ( fb .eq. 0. ) THEN
         xroot = b

      ELSE
         ! Ensure the supplied x-values yield y-values either side of the root
         IF ( sign(1.,fa) .eq. sign(1.,fb) ) THEN
            fa = func(a, aval, bval, cval, dval, eval)
            fb = func(b, aval, bval, cval, dval, eval)
         END IF

         c = b
         fc = fb

         DO iter = 1 , NITER
            ! Adjust new f(c) if it doesn't bracket the root with f(b)
            IF ( sign(1.,fb) .eq. sign(1.,fc) ) THEN
               c = a
               fc = fa
               d = b - a
               e = d
            END IF

            IF ( abs(fc) .lt. abs(fb) ) THEN
               a = b
               b = c
               c = a
               fa = fb
               fb = fc
               fc = fa
            END IF

            tol1 = 2. * EPS * abs(b) + 0.5 * tol
            xm = 0.5 * ( c - b )

            IF ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0. ) ) THEN
               xroot = b
            END IF

            IF ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) THEN
               IF ( a .eq. c ) THEN
                  p = 2. * xm * fb / fa
                  q = 1. - fb / fa
               ELSE
                  p = fb / fa * ( 2. * xm * fa / fc * ( fa / fc - fb / fc ) - &
                     ( b - a ) * ( fb / fc - 1. ) )
                  q = ( fa / fc - 1. ) * ( fb / fc - 1. ) * ( fb / fa - 1. )
               END IF

               IF ( p .gt. 0. ) q = -q
               p = abs( p )

               IF ( (2.*p) .lt. min( 3.*xm*q-abs(tol1*q) , abs(e*q) ) ) THEN
                  e = d
                  d = p / q
               ELSE
                  d = xm
                  e = d
               END IF

            ELSE
               d = xm
               e = d
            END IF

            a = b
            fa = fb

            IF ( abs(d) .gt. tol1 ) THEN
               b = b + d
            ELSE
               b = b + sign( tol1 , xm )
            END IF

            fb = func(b, aval, bval, cval, dval, eval)

         END DO

         xroot = b

      END IF

   END FUNCTION zbrent
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION solve_psi_leaf(p, psi_soil, e_leaf, kmax, b_plant, c_plant) &
      RESULT(resid)

      ! Calculation the approximate matching leaf water potential for a given
      ! transpiration per unit leaf.
      !
      ! N.B.: the exp() function is the same as in get_xylem_vulnerability,
      ! except that the type of p is different to avoid clashing with the
      ! interface in the zbrent.
      !
      ! Manon Sabot, 2023

      IMPLICIT NONE

      REAL, INTENT(IN) :: p, psi_soil, e_leaf, kmax, b_plant, c_plant

      ! local variables
      REAL :: resid

      ! residual must be ~0 if unit leaf trans matches hydraulics trans
      resid = ABS(e_leaf - kmax * min(1., &
         max(1E-9, exp(-(-p / b_plant)**c_plant))) * (psi_soil - p))

   END FUNCTION solve_psi_leaf
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION calc_psi_leaf_fails(p_crit, p_sat, e_leaf, kmax, b_plant, c_plant, pacc, N) &
      RESULT(psi_leaf)

      ! Calculation the approximate matching psi_leaf from the transpiration array.
      ! Manon Sabot, 2023

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: p_crit, p_sat, kmax, b_plant, c_plant, pacc
      REAL, DIMENSION(N), INTENT(IN) :: e_leaf

      ! local variables
      INTEGER :: iter
      REAL, DIMENSION(N) :: psi_leaf

      ! Infer the matching leaf water potential (MPa)
      DO iter = 1, N
         psi_leaf(iter) = zbrent(solve_psi_leaf, p_crit, p_sat, pacc, p_sat, &
            e_leaf(iter), kmax, b_plant, c_plant)
      END DO

   END FUNCTION calc_psi_leaf_fails
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION calc_psi_leaf(psi_sat, e_leaf, kmax, b_plant, c_plant, acc, N) &
      RESULT(psi_leaf)

      ! Calculation the approximate matching psi_leaf from the transpiration array.
      ! The Newton-Raphson procedure and associated methods fail to converge for
      ! this problem, and a z-brent also fails to secure a valid result for every
      ! index on the transpiration array, thus we use this trick.
      !
      ! Manon Sabot, 2023

      USE cable_def_types_mod

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: kmax, b_plant, c_plant, acc
      REAL(r_2), INTENT(IN) :: psi_sat
      REAL, DIMENSION(N), INTENT(IN) :: e_leaf

      ! local variables
      INTEGER :: iter
      REAL, DIMENSION(N) :: kc, kc_iter
      REAL(r_2), DIMENSION(N) :: psi_leaf

      ! Initialise the soil-plant hydraulic conductance at canopy xylem pressure
      ! mmol m-2 s-1 MPa-1
      kc = kmax * get_xylem_vulnerability(psi_sat, b_plant, c_plant)

      ! Infer the leaf water potentials via convergence of kc; this is not
      ! dissimilar to doing a series expansion
      DO iter = 1, N

         ! Update the leaf water potential estimates
         psi_leaf = psi_sat - e_leaf / kc

         ! Check the convergence on kc
         kc_iter = kmax * get_xylem_vulnerabilityx(psi_leaf, b_plant, c_plant)

         IF (iter > 2 .AND. ALL((ABS(kc - kc_iter)) < acc)) THEN

            EXIT

         END IF

         kc = kc_iter

      END DO

   END FUNCTION calc_psi_leaf
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION get_xylem_vulnerabilityx(psi, b_plant, c_plant) RESULT(weibull)

      ! Calculate the vulnerability to cavitation using a Weibull function
      !
      ! Martin De Kauwe, 27th August, 2020

      USE cable_def_types_mod

      IMPLICIT NONE

      REAL, INTENT(IN) :: b_plant, c_plant
      REAL(r_2), DIMENSION(:), INTENT(IN) :: psi

      ! local variables
      REAL, DIMENSION( SIZE(psi) ) :: weibull

      weibull = min(1., max(1E-9, exp(-(-psi / b_plant)**c_plant)))

   END FUNCTION get_xylem_vulnerabilityx
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION get_xylem_vulnerability(psi, b_plant, c_plant) RESULT(weibull)

      ! Calculate the vulnerability to cavitation using a Weibull function
      !
      ! Martin De Kauwe, 27th August, 2020

      USE cable_def_types_mod

      IMPLICIT NONE

      REAL, INTENT(IN) :: b_plant, c_plant
      REAL(r_2), INTENT(IN) :: psi

      ! local variables
      REAL :: weibull

      weibull = min(1., max(1E-9, exp(-(-psi / b_plant)**c_plant)))

   END FUNCTION get_xylem_vulnerability
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION calc_plc(kplant, kp_sat, PLC_crit) RESULT(plc)
      ! Calculates the percent loss of conductivity, PLC (-)
      !
      ! Manon Sabot, 2022

      IMPLICIT NONE

      REAL, INTENT(IN) :: kplant, kp_sat, PLC_crit

      ! local variables
      REAL :: plc

      ! Percent loss of conductivity (-)
      plc = 100. * (1. - kplant / kp_sat)
      plc = max(1E-9, min(PLC_crit, plc))

   END FUNCTION calc_plc
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION integrate_vulnerability(N, a, b, b_plant, c_plant) RESULT(value2)

      ! Approximate the integration with the mid-point rule
      !
      ! Martin De Kauwe, 27th August, 2020

      USE cable_def_types_mod

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: b_plant, c_plant
      REAL(r_2), INTENT(IN) :: a, b

      ! local variables
      INTEGER :: h
      REAL :: value1, value2

      value1 = 0.
      value2 = 0.

      DO h=1, N+1
         value1 = value1 + get_xylem_vulnerability(a + ((n - 0.5) * ((b - a) /&
            float(N))), b_plant, c_plant)
      END DO

      value2 = ((b - a) / float(N)) * value1

   END FUNCTION integrate_vulnerability
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION calc_transpiration(p, N, kmax, b_plant, c_plant) RESULT(e_leaf)

      ! At steady-state, transpiration is the integral of the plant's vulnerability
      ! curve from zero (no cuticular conductance) to its maximum (e_crit)
      ! (Sperry & Love 2015). By integrating across the soil–plant vulnerability
      ! curve, the relation between transpiration and a given total pressure drop
      ! can be found.
      !
      ! Martin De Kauwe, 27th August, 2020

      USE cable_def_types_mod

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: b_plant, c_plant, kmax
      REAL(r_2), DIMENSION(N), INTENT(IN) :: p

      ! local variables
      INTEGER :: h
      REAL, DIMENSION(N) :: e_leaf
      REAL, PARAMETER :: MMOL_2_MOL = 1E-3

      ! integrate over the full range of water potentials from psi_soil to e_crit
      DO h=1, N
         e_leaf(h) = integrate_vulnerability(N, p(h), p(1), b_plant, c_plant)

         IF (e_leaf(h) > 1E-9) then
            e_leaf(h) = e_leaf(h) * kmax * MMOL_2_MOL ! mol m-2 s-1
         END IF

      END DO

   END FUNCTION calc_transpiration
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION assim(Ci, gamma_star, a1, a2) RESULT(assimilation)

      ! Calculation of photosynthesis with the limitation defined by the variables
      ! passed as a1 and a2, i.e. if we are calculating vcmax or jmax limited
      ! assimilation rates.
      !
      ! Martin De Kauwe, 27th August, 2020

      IMPLICIT NONE

      REAL, INTENT(IN) :: gamma_star, a1, a2
      REAL, DIMENSION(:), INTENT(IN) :: Ci

      ! local variables
      REAL, DIMENSION( SIZE(Ci) ) :: assimilation

      assimilation = a1 * (Ci - gamma_star) / (a2 + Ci)

   END FUNCTION assim
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION QUADP(A,B,C) RESULT(root)

      ! Solves the quadratic equation - finds larger root.

      IMPLICIT NONE

      REAL, DIMENSION(:), INTENT(IN) :: B,C

      ! local variables
      REAL :: A
      REAL, DIMENSION( SIZE(B) ) :: d, root

      d = B*B - 4. * A * C ! discriminant

      where (d < 0.)
         root = 0.
      end where

      root = (- B + SQRT(B*B - 4*A*C)) / (2.*A)
      where (A == 0. .AND. B  == 0.)
         root = 0.
      elsewhere (A == 0. .AND. B > 0.)
         root = -C/B
      end where

   END FUNCTION QUADP
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION calc_michaelis_menten_constants(Tleaf) RESULT(Km)

      REAL, INTENT(IN) :: tleaf

      ! local variables
      REAL :: Kc, Ko, Km, Kc25, Ko25, Ec, Eo, Oi

      ! Michaelis-Menten coeffs for carboxylation by Rubisco at 25 C or 298 K,
      ! umol mol-1
      Kc25 = 404.9

      ! Michaelis-Menten coefficents for oxygenation by Rubisco at 25 C, mmol mol-1
      ! Note value in Bernacchie 2001 is in mmol
      Ko25 = 278.4

      ! Activation energy for carboxylation, J mol-1
      Ec = 79430.

      ! Activation energy for oxygenation, J mol-1
      Eo = 36380.

      ! intercellular concentration of O2, mmol mol-1
      Oi = 210.

      Kc = arrh(Kc25, Ec, Tleaf)
      Ko = arrh(Ko25, Eo, Tleaf)

      Km = Kc * (1. + Oi / Ko)

   END FUNCTION calc_michaelis_menten_constants
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION arrh(k25, Ea, Tk) RESULT(a)

      REAL, INTENT(IN) :: k25, Ea, Tk

      ! local variables
      REAL :: a
      REAL, PARAMETER :: RGAS = 8.314

      a = k25 * exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))

   END FUNCTION arrh
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
   FUNCTION peaked_arrh(k25, Ea, Tk, deltaS, Hd) RESULT(pa)

      REAL, INTENT(IN) :: k25, Ea, Tk, deltaS, Hd

      ! local variables
      REAL :: arg1, arg2, arg3, pa
      REAL, PARAMETER :: RGAS = 8.314

      arg1 = arrh(k25, Ea, Tk)
      arg2 = 1. + exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
      arg3 = 1. + exp((Tk * deltaS - Hd) / (Tk * RGAS))

      pa = arg1 * arg2 / arg3

   END FUNCTION peaked_arrh
! ------------------------------------------------------------------------------


END MODULE cable_veg_hydraulics_module
