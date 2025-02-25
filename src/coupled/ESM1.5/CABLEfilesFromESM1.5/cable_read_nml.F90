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
! Purpose: Routines to read CABLE namelist, check variables, allocate and 
!          deallocate CABLE arrays
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Rewrite of code from v1.8 (ACCESS1.3)
!          soil_snow_type now ssnow (instead of ssoil)
!
!
! ==============================================================================

MODULE cable_read_nml_mod
   
CHARACTER(LEN=*), PARAMETER :: runtime_vars_file="cable.nml"

INTERFACE check_nmlvar 
  MODULE PROCEDURE check_chvar, check_intvar, check_lgvar
END INTERFACE check_nmlvar 
 
CONTAINS

SUBROUTINE cable_read_nml() 

USE cable_common_module,    ONLY: cable_runtime, filename,                     & 
                                  redistrb, l_casacnp, l_laiFeedbk,            &
                                  l_vcmaxFeedbk, l_luc, l_thinforest
USE cable_runtime_opts_mod, ONLY: wiltParam, satuParam, snmin, cable_user
USE casavariable,           ONLY: casafile
USE casadimension,          ONLY: icycle

IMPLICIT NONE

INTEGER :: funit=88813

!--- namelist for CABLE runtime vars, switches 
NAMELIST/CABLE/filename, l_thinforest, l_luc, l_laiFeedbk, l_vcmaxFeedbk,      &
               icycle, casafile, cable_user, redistrb, wiltParam, satuParam,   &
               snmin, l_casacnp

!--- assume namelist exists. no iostatus check 
OPEN(unit=funit,FILE= runtime_vars_file)
  READ(funit,NML=CABLE)
  PRINT *, '  '; PRINT *, 'CABLE_log:' 
  PRINT *, '  Opened file - '
  PRINT *, '  ', TRIM(runtime_vars_file)
  PRINT *, '  for reading runtime vars.' 
  PRINT *, 'End CABLE_log:'; PRINT *, '  '
CLOSE(funit)

!jh:revise list, include error traps 
!--- check value of variable 
CALL check_nmlvar('cable_user%DIAG_SOIL_RESP', cable_user%DIAG_SOIL_RESP)
CALL check_nmlvar('cable_user%LEAF_RESPIRATION',                         &
                  cable_user%LEAF_RESPIRATION)
CALL check_nmlvar('cable_user%FWSOIL_SWITCH', cable_user%FWSOIL_SWITCH)
CALL check_nmlvar('cable_user%RUN_DIAG_LEVEL', cable_user%RUN_DIAG_LEVEL)
CALL check_nmlvar('cable_user%l_new_roughness_soil',                     &
                   cable_user%l_new_roughness_soil)
CALL check_nmlvar('cable_user%l_new_roughness_soil',                     &
                   cable_user%l_new_roughness_soil)
CALL check_nmlvar('cable_user%l_new_roughness_soil',                     &
                         cable_user%l_new_roughness_soil)

PRINT *, '  ' 
PRINT *, 'CASA_log:'
PRINT *, '  icycle =        ', icycle
PRINT *, '  l_laiFeedbk =   ', l_laiFeedbk
PRINT *, '  l_vcmaxFeedbk = ', l_vcmaxFeedbk
PRINT *, '  l_luc =         ', l_luc
PRINT *, 'End CASA_log:' 
PRINT *, '  '

IF (l_casacnp  .AND. (icycle == 0 .OR. icycle > 3)) THEN 
  STOP 'CASA_log: icycle must be 1 to 3 when using casaCNP'
END IF  
IF ((.NOT. l_casacnp)  .AND. (icycle >= 1)) THEN
  STOP 'CASA_log: icycle must be <=0 when not using casaCNP'
END IF  
IF ((l_laiFeedbk .OR. l_vcmaxFeedbk) .AND. (.NOT. l_casacnp)) THEN
  STOP 'CASA_log: casaCNP required to get prognostic LAI or Vcmax'
END IF  
IF (l_vcmaxFeedbk .AND. icycle < 2) THEN
  STOP 'CASA_log: icycle must be 2 to 3 to get prognostic Vcmax'
END IF  

RETURN
END SUBROUTINE cable_read_nml

!jhan: also add real, logical, int interfaces
SUBROUTINE check_chvar(this_var, val_var)

CHARACTER(LEN=*), INTENT(IN) :: this_var, val_var 

PRINT *, '  ' 
PRINT *, 'CABLE_log:' 
PRINT *, '   run time variable - '
PRINT *, '  ', TRIM(this_var) 
PRINT *, '   defined as - '
PRINT *, '  ', TRIM(val_var) 
PRINT *, 'End CABLE_log:' 
PRINT *, '  '

RETURN
END SUBROUTINE check_chvar

SUBROUTINE check_intvar(this_var, val_var)

CHARACTER(LEN=*), INTENT(IN) :: this_var
INTEGER, INTENT(IN) :: val_var 

PRINT *, '  ' 
PRINT *, 'CABLE_log:' 
PRINT *, '   run time variable - '
PRINT *, '  ', TRIM(this_var) 
PRINT *, '   defined as - '
PRINT *, '  ', val_var
PRINT *, 'End CABLE_log:'
PRINT *, '  '

RETURN
END SUBROUTINE check_intvar

SUBROUTINE check_lgvar(this_var, val_var)

CHARACTER(LEN=*), INTENT(IN) :: this_var
LOGICAL, INTENT(IN) :: val_var

PRINT *, '  '
PRINT *, 'CABLE_log:'
PRINT *, '   run time variable - '
PRINT *, '  ', TRIM(this_var)
PRINT *, '   defined as - '
PRINT *, '  ', val_var
PRINT *, 'End CABLE_log:' 
PRINT *, '  '

RETURN
END SUBROUTINE check_lgvar

END MODULE cable_read_nml_mod





