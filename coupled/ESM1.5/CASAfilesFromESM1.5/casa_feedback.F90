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
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(IN) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet
  real, dimension(17)                   ::  nintercept,nslope,xnslope
! Will move to look-up table in later version
  data nintercept/6.32,4.19,6.32,5.73,14.71,6.42,2.00,14.71,4.71,14.71,14.71,7.00,14.71,14.71,14.71,14.71,14.71/
!! r29
  data nslope/9.00,28.00,12.00,40.00,10.00,25.00,20.00,4.00,45.00,23.15,23.15,10.00,23.15,23.15,23.15,23.15,23.15/
! Modified further in ACCESS-1.4
  data xnslope/0.80,1.00,2.00,1.00,1.00,1.00,0.50,1.00,0.34,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx  ! local variables

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
        npleafx(np) = MIN(30.0,MAX(8.0,casapool%nplant(np,leaf) &
                                      /casapool%pplant(np,leaf)))
      ENDIF
    ENDIF

    IF (casamet%glai(np) > casabiome%glaimin(ivt)) THEN
      IF (ivt/=2) THEN
        veg%vcmax(np) = ( nintercept(ivt) &
                        + nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
      ELSE
        IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
          veg%vcmax(np) = ( nintercept(ivt)  &
                          + nslope(ivt)*(0.4+9.0/npleafx(np)) &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ELSE
          veg%vcmax(np) = ( nintercept(ivt) &
                          + nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.0e-6
        ENDIF
      ENDIF
      veg%vcmax(np) = veg%vcmax(np) * xnslope(ivt)
    ENDIF

  ENDDO

  veg%ejmax = 2.0 * veg%vcmax
 END SUBROUTINE casa_feedback

End module feedback_mod



