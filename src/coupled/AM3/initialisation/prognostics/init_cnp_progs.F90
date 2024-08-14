MODULE init_cnp_progs_mod

IMPLICIT NONE

PUBLIC :: init_cnp_progs
PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_CNP_PROGS_MOD'

CONTAINS

SUBROUTINE init_cnp_progs( land_pts, nsurft, progs_cnp_vars,                   &
                           progs_cnp_vars_data )

USE progs_cnp_vars_mod,           ONLY: progs_cnp_vars_type,                   &
                                        progs_cnp_vars_data_type,              &
                                        progs_cnp_vars_alloc,                  &
                                        progs_cnp_vars_assoc

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nsurft
TYPE(progs_cnp_vars_data_type), INTENT(OUT) :: progs_cnp_vars_data
TYPE(progs_cnp_vars_type),      INTENT(OUT) :: progs_cnp_vars

CALL progs_cnp_vars_alloc(land_pts, nsurft, progs_cnp_vars_data )
CALL progs_cnp_vars_assoc(progs_cnp_vars, progs_cnp_vars_data )

RETURN

END SUBROUTINE init_cnp_progs

END MODULE init_cnp_progs_mod
