
MODULE cable_fields_mod

! module containing instances of the data types for jules standalone

!type definitions
USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type, progs_cbl_vars_data_type

!TYPES to hold the data
TYPE(progs_cbl_vars_data_type), TARGET :: progs_cbl_vars_data

!TYPES we pass around. These happen to be pointers to the data types above
!but this should be transparent
TYPE(progs_cbl_vars_type) ::progs_cbl_vars

END MODULE cable_fields_mod
