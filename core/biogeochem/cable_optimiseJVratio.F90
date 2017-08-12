
!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: acclimation of ratio of Jmax to Vcmax, with corresponding optimisation
! of contributions of electron-tranport and carboxylation rates to rate of net photosynthesis
!
! Called from: SUBROUTINE bgcdriver in casa_cable.F90
!
! History: Vanessa Haverd Aug 2017

! ==============================================================================
MODULE cable_optimise_JV_module

 Use cable_def_types_mod, ONLY: met_type, climate_type, canopy_type, veg_parameter_type, &
      mp, r_2
 
 USE cable_data_module, ONLY : icanopy_type, point2constants
 USE TypeDef,              ONLY: i4b, dp

 TYPE( icanopy_type ) :: C

 ! variables local to module
 REAL, ALLOCATABLE :: APAR(:), Dleaf(:), Tleaf(:), cs(:), scalex(:)
 REAL :: Anet, vcmax00, bjv, g1, fwsoil, Kc0, Ko0, ekc, eko, alpha
 REAL ::     convex, Neff, Rd0
 INTEGER :: nt
 !REAL, PARAMETER :: relcost_J = 1.6 ! Chen et al. Oecologia, 1993, 93: 63-69
 REAL, PARAMETER :: relcost_J = 2.0
CONTAINS
! ==============================================================================


SUBROUTINE optimise_JV (veg, climate, ktauday, bjvref)

  IMPLICIT NONE

  TYPE (veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
  TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables
  INTEGER, INTENT(IN) :: ktauday
  REAL, INTENT(IN) :: bjvref
  INTEGER:: k
  REAL :: Anet_cost, bjv_new
  REAL, PARAMETER :: l_bound = 0.5
  REAL, PARAMETER :: u_bound = 3.0

  nt = ktauday * 5
  ALLOCATE(APAR(nt))
  ALLOCATE(Dleaf(nt))
  ALLOCATE(Tleaf(nt))
  ALLOCATE(cs(nt))
  ALLOCATE(scalex(nt))

  ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)


    
    DO k=1,mp
    !DO k=7,7
       if (veg%frac4(k).lt.0.001) then ! not C4
          vcmax00 = veg%vcmax(k) ! vcmax at standard temperature, and bjvref
          Kc0 = veg%conkc0(k)
          Ko0 = veg%conko0(k)
          ekc = veg%ekc(k)
          eko = veg%eko(k)
          g1 = veg%g1(k)
          Rd0 = veg%cfrd(k) * veg%vcmax(k)
          fwsoil =(climate%dmoist_31(k,31))
          alpha = climate%frec(k)*veg%alpha(k) ! quantum efficiency for
          ! electron transport (adjusted down for spring recovery in conifers)
          convex = veg%convex(k) 
          Neff = vcmax00 + relcost_J*bjvref*vcmax00/4. ! effective nitrogen amount 
          !for distribution between e-limited and c-limited processes
 

          ! optimisation for shade leaves
          APAR = climate%APAR_leaf_shade(k,:)*1e-6;
          Dleaf = max(climate%Dleaf_shade(k,:), 50.0)*1e-3 ! Pa -> kPa
          Tleaf = climate%Tleaf_shade(k,:)
          cs = climate%cs_shade(k,:)*1e-6
          scalex = climate%scalex_shade(k,:)

!write(*,*) 'init guesses',total_photosynthesis_cost(bjvref) , &
!total_photosynthesis_cost(l_bound), total_photosynthesis_cost(u_bound)
!write(*,*)
          if(total_photosynthesis_cost(bjvref).lt.total_photosynthesis_cost(l_bound).and. &
               total_photosynthesis_cost(bjvref).lt.total_photosynthesis_cost(u_bound)) then

             Anet_cost = golden(l_bound,bjvref,u_bound,total_photosynthesis_cost,0.01,bjv_new)
             veg%vcmax_shade(k) = Neff/(1.+relcost_J*bjv_new/4.0)
             veg%ejmax_shade(k) = veg%vcmax_shade(k)*bjv_new
          else
             bjv_new = bjvref
             veg%vcmax_shade(k) = veg%vcmax(k)
             veg%ejmax_shade(k) = veg%ejmax(k)
          endif
          
          !Anet = total_photosynthesis(bjv_new)
!write(*,*) bjv_new, Anet,  total_photosynthesis(bjvref), veg%vcmax_shade(k),  veg%ejmax_shade(k),veg%vcmax(k),   veg%ejmax(k)
   
          ! optimisation for sun leaves
          APAR = climate%APAR_leaf_sun(k,:)*1e-6;
          Dleaf = max(climate%Dleaf_sun(k,:), 50.0)*1e-3 ! Pa -> kPa
          Tleaf = climate%Tleaf_sun(k,:)
          cs = climate%cs_sun(k,:)*1e-6
          scalex = climate%scalex_sun(k,:)

!if (k==7) then
!   write(5333,"(200e16.6)") APAR
!   write(5334,"(200e16.6)") Dleaf
!   write(5335,"(200e16.6)") Tleaf
!   write(5336,"(200e16.6)") cs
!   write(5337,"(200e16.6)") scalex
!endif

          if(total_photosynthesis_cost(bjvref).lt.total_photosynthesis_cost(l_bound).and. &
               total_photosynthesis_cost(bjvref).lt.total_photosynthesis_cost(u_bound)) then

             Anet_cost = golden(l_bound,bjvref,u_bound,total_photosynthesis_cost,0.01,bjv_new)
             veg%vcmax_sun(k) = Neff/(1.+relcost_J*bjv_new/4.0)
             veg%ejmax_sun(k) = veg%vcmax_sun(k)*bjv_new
          else
             bjv_new = bjvref
             veg%vcmax_sun(k) = veg%vcmax(k)
             veg%ejmax_sun(k) = veg%ejmax(k)
          endif

       else
          bjv_new = bjvref
          veg%vcmax_shade(k) = veg%vcmax(k)
          veg%ejmax_shade(k) = veg%ejmax(k)

          veg%vcmax_sun(k) = veg%vcmax(k)
          veg%ejmax_sun(k) = veg%ejmax(k)
       endif
!           if (k==1 ) then
!              write(99,"(200e16.6)") bjv_new, total_photosynthesis(bjv_new), maxval(APAR*1e6), maxval(Tleaf) - 273.15, maxval(Dleaf), maxval(scalex), maxval(cs)*1e6, total_photosynthesis(1.38)
!           endif


       !if (k==7) then
       !   write(6333,"(200e16.6)") APAR
       !   write(6334,"(200e16.6)") Dleaf
       !   write(6335,"(200e16.6)") Tleaf
       !   write(6336,"(200e16.6)") cs
       !   write(6337,"(200e16.6)") scalex
       !endif


!if (k==7) then
!   write(4333,"(200e16.6)") vcmax00, g1, fwsoil, alpha, convex, Rd0, veg%vcmax_shade(k) , &
!        veg%ejmax_shade(k) , veg%vcmax_sun(k) , veg%ejmax_sun(k) 
!endif

    ENDDO

DEALLOCATE(APAR)
DEALLOCATE(Dleaf)
DEALLOCATE(Tleaf)
DEALLOCATE(cs)
DEALLOCATE(scalex)

END SUBROUTINE optimise_JV
 ! ------------------------------------------------------------------------------
subroutine fAn(a,b,c,A2)
  REAL, INTENT(IN) :: a,b,c
  REAL, INTENT(OUT) :: A2
  REAL :: s2
  s2 = b**2 - 4.*a*c
  A2 = (-b - sqrt(s2))/(2.*a)
end subroutine fAn
 ! ------------------------------------------------------------------------------
subroutine fabc(Cs,g0,x,gamma,beta,Gammastar,Rd,a,b,c)
  REAL, INTENT(IN) :: Cs,g0,x,gamma,beta,Gammastar,Rd
  REAL, INTENT(OUT) :: a,b,c
  a = (1.-x)*Cs - x*beta
  b = -g0*Cs**2 + ((1.-x)*(Rd-gamma)-g0*beta)*Cs - x*(gamma*Gammastar+Rd*beta)
  c = -g0*(Rd-gamma)*Cs**2 - g0*(gamma*Gammastar+Rd*beta)*Cs
end subroutine fabc
 ! ------------------------------------------------------------------------------

REAL FUNCTION total_photosynthesis_cost(bjv)
USE cable_canopy_module, ONLY :  xvcmxt3,xejmxt3, ej3x, xrdt 
!TYPE( icanopy_type ) :: C
  REAL, INTENT(IN) :: bjv
INTEGER :: k, j
REAL :: kct, kot, tdiff
REAL :: g0, x, gamma,  beta, gammastar, Rd, a, b, c1, jmaxt
REAL:: Anc, Ane, vcmax0
REAL :: An(nt)

CALL point2constants(C)
!write(*,*) C%TrefK
An = 0.0
j = 1
vcmax0 = Neff/(1.+relcost_J*bjv/4.0);

DO k=1,nt
   if (APAR(k) .gt. 60e-6) then

      g0 = 0.0
      x = 1.0  + (g1 * fwsoil) / SQRT(Dleaf(k))
      gamma =  Vcmax0*scalex(k)*xvcmxt3(Tleaf(k)) 
      tdiff = Tleaf(k) - C%Trefk
      gammastar = C%gam0 * ( 1.0 + C%gam1 * tdiff                  &
                                          + C%gam2 * tdiff * tdiff )
      Rd  =Rd0*scalex(k)*xrdt(Tleaf(k))
      kct = kc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
                                              * ( 1.0 - C%trefk/Tleaf(k) ) )
      kot = ko0 *EXP( ( eko / (C%rgas*C%trefk) ) &
                                              * ( 1.0 - C%trefk/Tleaf(k) ) )

      beta = kct * (1.0+0.21/kot)
      CALL fabc(cs(k), g0, x, gamma, beta, gammastar, Rd, a, b, c1)
      CALL fAn(a,b,c1,Anc) ! rubisco-limited
      jmaxt = bjv*Vcmax0*scalex(k)*xejmxt3(Tleaf(k))
      gamma = ej3x(APAR(k), alpha,convex, jmaxt) 
      beta = 2.0 * gammastar
      CALL fabc(cs(k), g0, x, gamma, beta, gammastar, Rd, a, b, c1)
      CALL fAn(a,b,c1,Ane) ! e-transport limited
      An(j) = min(Anc, Ane)
   else
      An(j) = 0.0
   endif
   j = j+1
ENDDO

if (sum(An) > 0.0) then
   total_photosynthesis_cost = (relcost_J*Vcmax0*bjv/4.0 + Vcmax0)/sum(An)
else
   total_photosynthesis_cost = 0.0
endif

END FUNCTION total_photosynthesis_cost

REAL FUNCTION total_photosynthesis(bjv)
USE cable_canopy_module, ONLY :  xvcmxt3,xejmxt3, ej3x, xrdt
!TYPE( icanopy_type ) :: C
REAL, INTENT(IN) :: bjv
INTEGER :: k, j
REAL :: kct, kot, tdiff
REAL :: g0, x, gamma,  beta, gammastar, Rd, a, b, c1, jmaxt
REAL:: Anc, Ane, vcmax0
REAL :: An(nt)

CALL point2constants(C)
!write(*,*) C%TrefK
An = 0.0
j = 1
vcmax0 = Neff/(1.+relcost_J*bjv/4.0);

DO k=1,nt
   if (APAR(k) .gt. 60e-6) then

      g0 = 0.0
      x = 1.0  + (g1 * fwsoil) / SQRT(Dleaf(k))
      gamma =  Vcmax0*scalex(k)*xvcmxt3(Tleaf(k))
      tdiff = Tleaf(k) - C%Trefk
      gammastar = C%gam0 * ( 1.0 + C%gam1 * tdiff                  &
                                          + C%gam2 * tdiff * tdiff )
      Rd  =Rd0*scalex(k)*xrdt(Tleaf(k))
      kct = kc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
                                              * ( 1.0 - C%trefk/Tleaf(k) ) )
      kot = ko0 *EXP( ( eko / (C%rgas*C%trefk) ) &
                                              * ( 1.0 - C%trefk/Tleaf(k) ) )

      beta = kct * (1.0+0.21/kot)
      CALL fabc(cs(k), g0, x, gamma, beta, gammastar, Rd, a, b, c1)
      CALL fAn(a,b,c1,Anc) ! rubisco-limited
      jmaxt = bjv*Vcmax0*scalex(k)*xejmxt3(Tleaf(k))
      gamma = ej3x(APAR(k), alpha,convex, jmaxt) 
      beta = 2.0 * gammastar
      CALL fabc(cs(k), g0, x, gamma, beta, gammastar, Rd, a, b, c1)
      CALL fAn(a,b,c1,Ane) ! e-transport limited
      An(j) = min(Anc, Ane)
      write(3368,"(200e16.6)") An(j), Ane,Anc,Tleaf(k),APAR(k),x,cs(k),scalex(k)
   else
      An(j) = 0.0
   endif
   j = j+1
ENDDO


   total_photosynthesis = sum(An)


END FUNCTION total_photosynthesis


	FUNCTION golden(ax,bx,cx,func,tol,xmin)
	!USE nrtype
	IMPLICIT NONE
	REAL, INTENT(IN) :: ax,bx,cx,tol
	REAL, INTENT(OUT) :: xmin
	REAL :: golden
	INTERFACE
		FUNCTION func(x)
		!USE nrtype
		IMPLICIT NONE
		REAL, INTENT(IN) :: x
		REAL :: func
		END FUNCTION func
	END INTERFACE
	REAL, PARAMETER :: R=0.61803399,C=1.0-R
	REAL :: f1,f2,x0,x1,x2,x3
	x0=ax
	x3=cx
	if (abs(cx-bx) > abs(bx-ax)) then
		x1=bx
		x2=bx+C*(cx-bx)
	else
		x2=bx
		x1=bx-C*(bx-ax)
	end if
	f1=func(x1)
	f2=func(x2)
	do
		if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
		if (f2 < f1) then
			call shft3(x0,x1,x2,R*x2+C*x3)
			call shft2(f1,f2,func(x2))
		else
			call shft3(x3,x2,x1,R*x1+C*x0)
			call shft2(f2,f1,func(x1))
		end if
	end do
	if (f1 < f2) then
		golden=f1
		xmin=x1
	else
		golden=f2
		xmin=x2
	end if
	CONTAINS
!BL
	SUBROUTINE shft2(a,b,c)
	REAL, INTENT(OUT) :: a
	REAL, INTENT(INOUT) :: b
	REAL, INTENT(IN) :: c
	a=b
	b=c
	END SUBROUTINE shft2
!BL
	SUBROUTINE shft3(a,b,c,d)
	REAL, INTENT(OUT) :: a
	REAL, INTENT(INOUT) :: b,c
	REAL, INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft3
	END FUNCTION golden



! ==============================================================================
END MODULE cable_optimise_JV_module
