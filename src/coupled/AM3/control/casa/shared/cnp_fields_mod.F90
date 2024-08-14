MODULE cnp_fields_mod

! TYPE definitions
USE progs_cnp_vars_mod,       ONLY: progs_cnp_vars_type,                       &
                                    progs_cnp_vars_data_type

USE work_vars_mod_cnp,        ONLY: work_vars_type
                                    
PUBLIC

! instantiated TYPES to hold the data
TYPE(progs_cnp_vars_data_type), TARGET :: progs_cnp_vars_data

! instantiated TYPES we pass around. These happen to be pointers to the data
! types above but this should be transparent
TYPE(progs_cnp_vars_type) :: progs_cnp_vars
TYPE(work_vars_type)      :: work_vars_cnp

END MODULE cnp_fields_mod
