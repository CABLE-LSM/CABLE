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
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_write_logs_mod

contains

subroutine write_progress_log( tlog, ktau, subr_name, cycleno, vars )
  USE cable_decs_mod, ONLY : CABLE_file
  implicit none
  type(CABLE_file) :: tlog
  integer :: ktau
  character(len=*) :: subr_name 
  integer :: cycleno
  logical, optional :: vars 
  integer :: funit
  logical :: already_open 
  logical, save :: first_call =.true.
  real :: fcputime
  
    funit = tlog%funit
    inquire(unit=funit, opened=already_open)
    call cpu_time(fcputime)
   
    if( tlog%want .AND. .NOT. already_open ) & 
      open( unit=funit,                                                        & 
            file=(trim(tlog%folder)//trim(tlog%filename)),                     &
            access='append', status="unknown", action="write" )
 
 
    if ( subr_name=='cable_logs' .AND. first_call ) then
        write(funit,*) "-------------------------------------------------"
        write(funit,*) "ENDGAME numcycles: 2"!, numcycles  
        write(funit,*) "-------------------------------------------------"
        first_call =.false.
    endif
  

    if ( .NOT. present(vars) ) then
      write(funit,*) "-------------------------------------------------"
      write(funit,*) "timestep:         ", ktau 
      write(funit,*) "subroutine:       ", subr_name 
      write(funit,*) "ENDGAME cycle no: ", cycleno
      write(funit,*) "cputime:          ", fcputime 
      write(funit,*) "-------------------------------------------------"
      write(funit,*) ""
    endif

    if ( tlog%vars ) then
      if ( tlog%itau > 0 ) call generic_xtau( tlog%itau, funit )
      if ( tlog%mtau > 0 ) call generic_xtau( tlog%mtau, funit )
      if ( tlog%ftau > 0 ) call generic_xtau( tlog%ftau, funit )
    endif

End subroutine write_progress_log
    
subroutine generic_xtau( tau, funit )
  integer :: tau, funit 
      if ( tau > 0 ) then
        write(funit,*) "Generating logs for timestep: ", tau
        !jhan:TODO - write vars for itau
      else
        write(funit,*) "Specified *tau is a timestep and must be > 0 "
        write(funit,*) "to write logs "
      endif
End subroutine generic_xtau

End module cable_write_logs_mod
