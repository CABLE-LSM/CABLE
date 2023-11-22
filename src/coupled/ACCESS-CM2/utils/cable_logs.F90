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
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: UM code sf_exch
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

module cable_logs_mod
   
  !--- USE CABLE modules 
  USE cable_decs_mod, ONLY : CABLE_file
  
  implicit none
   
  logical, save :: first_CABLE_call = .true.
  
  ! CABLE progress log - file details (from:cable_decs_mod) 
  type(CABLE_file), save :: KBL_progress

contains

SUBROUTINE cable_logs( )

  !USE UM modules 
  USE um_parcore, ONLY : mype                   ! processor number
  USE timestep_mod,   ONLY : timestep_number, & ! number
                             timestep           ! width in S 
 USE model_time_mod, ONLY: &
    target_end_stepim 
      
  !USE CABLE modules 
  USE cable_common_module, ONLY : knode_gl,        & ! processor number
                                  ktau_gl,         & ! number
                                  kwidth_gl,       & ! width in S 
                                  kend_gl

  IMPLICIT NONE
 
  !--- IN ARGS FROM surf_couple_explicit() -----------------------------------------------
  
  !--- End INPUT ARGS FROM sf_exch() -------------------------------------------

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_logs"
   
  !cludge 
  INTEGER :: ioerror=0
  
  !--- namelist for explicit driver checks 
  NAMELIST/cable_logs_nml/KBL_progress
    
  !--- End header
  
  !CABLE uses k****_gl throughout   
  knode_gl  = mype
  kwidth_gl = int(timestep)
  kend_gl   = target_end_stepim(1) !UM uses index atmos_im=1 

  if( knode_gl==0 ) then
    !defaults
    KBL_progress%want = .true.
    KBL_progress%funit = 12832
    KBL_progress%folder   = "./"
    KBL_progress%filename = "CABLE_progress.log"
    KBL_progress%itau = 1
    KBL_progress%mtau = -1
    KBL_progress%ftau = 2 
  endif    

  !jhan:TODO: read namelist possibly overwriting above defaults
  
  first_CABLE_call = .false.

return

END SUBROUTINE cable_logs

subroutine init_nml_prefs(xin, xout, xinout, xname)  
  USE cable_decs_mod, ONLY : CABLE_file
  implicit none
  type(CABLE_file) :: xin
  type(CABLE_file), optional :: xout, xinout
  character(len=*) :: xname
  integer, save :: funit = 31465

  !defaults
  funit = funit+1 
  xin%want = .true.
  xin%funit = funit
  xin%folder   = "./"
  xin%filename =trim( xname //"_IN")
  xin%vars = .true.
  xin%itau = 1
  xin%mtau = -1
  xin%ftau = 2 
  if(present(xout)) then 
    funit = funit+1 
    xout%want = .true.
    xout%funit = funit
    xout%folder   = "./"
    xout%filename =trim( xname //"_OUT")
    xout%vars = .true.
    xout%itau = 1
    xout%mtau = -1
    xout%ftau = 2 
  endif 
  if(present(xinout)) then 
    funit = funit+1 
    xinout%want = .true.
    xinout%funit = funit
    xinout%folder   = "./"
    xinout%filename =trim( xname //"_INOUT")
    xinout%vars = .true.
    xinout%itau = 1
    xinout%mtau = -1
    xinout%ftau = 2 
  endif 
  
  !jhan:TODO - read namelists
    
    
End subroutine init_nml_prefs
 
End module cable_logs_mod
