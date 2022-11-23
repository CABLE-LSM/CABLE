!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: bgcdriver - interface between casacnp and cable
!          sumcflux  - accumulating carbon fluxes (not required for UM)
!
! Called from: cable_driver for offline version
!              Not currently called/available for ACCESS version
!
! Contact: Yingping.Wang@csiro.au
!
! History: Model development by Yingping Wang, coupling to Mk3L by Bernard Pak
!          ssoil changed to ssnow
!
! ==============================================================================

module feedback_mod

contains

 SUBROUTINE casa_feedback(ktau,veg,casabiome,casapool,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE casa_cnp_module, ONLY: vcmax_np
  USE cable_common_module,  ONLY:  CABLE_USER
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx, pleafx, nleafx ! local variables
  real, dimension(17)                   ::  xnslope
  data xnslope/0.80,1.00,2.00,1.00,1.00,1.00,0.50,1.00,0.34,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/

  ! first initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  npleafx = 14.2

  DO np=1,mp
    ivt=veg%iveg(np)
    IF (casamet%iveg2(np)/=icewater &
        .AND. casamet%glai(np)>casabiome%glaimin(ivt)  &
        .AND. casapool%cplant(np,leaf)>0.0) THEN
      ncleafx(np) = MIN(casabiome%ratioNCplantmax(ivt,leaf), &
                        MAX(casabiome%ratioNCplantmin(ivt,leaf), &
                            casapool%nplant(np,leaf)/casapool%cplant(np,leaf)))
      IF (icycle>2 .AND. casapool%pplant(np,leaf)>0.0) THEN
        npleafx(np) = MIN(30.0,MAX(8.0,real(casapool%nplant(np,leaf) &
                /casapool%pplant(np,leaf))))
      ENDIF
    ENDIF

    IF (TRIM(cable_user%vcmax).eq.'standard') then
       IF (casamet%glai(np) > casabiome%glaimin(ivt)) THEN
          IF (ivt/=2) THEN
             veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                  + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
          ELSE
             IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
                veg%vcmax(np) = ( casabiome%nintercept(ivt)  &
                     + casabiome%nslope(ivt)*(0.4+9.0/npleafx(np)) &
                     * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
             ELSE
                veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                     + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.0e-6
             ENDIF
          ENDIF
          veg%vcmax(np) =veg%vcmax(np)* xnslope(ivt)
       ENDIF

    elseif (TRIM(cable_user%vcmax).eq.'Walker2014') then
       !Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax –
       !to leaf nitrogen, leaf phosphorus, and specific leaf area:
       !a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
       ! veg%vcmax(np) = exp(3.946 + 0.921*log(nleafx(np)) + 0.121*log(pleafx(np)) + &
       !      0.282*log(pleafx(np))*log(nleafx(np))) * 1.0e-6
       nleafx(np) = ncleafx(np)/casabiome%sla(ivt) ! leaf N in g N m-2 leaf
       pleafx(np) = nleafx(np)/npleafx(np) ! leaf P in g P m-2 leaf
       if (ivt .EQ. 7) then
          veg%vcmax(np) = 1.0e-5 ! special for C4 grass: set here to value from  parameter file
       else
          veg%vcmax(np) = vcmax_np(nleafx(np), pleafx(np))
       endif
    else
       stop ('invalid vcmax flag')
    endif

  ENDDO

  veg%ejmax = 2.0 * veg%vcmax

END SUBROUTINE casa_feedback

End module feedback_mod
