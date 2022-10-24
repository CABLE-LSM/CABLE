MODULE cable_fields_mod

! module containing instances of the data types for CABLE in jules standalone
! analagous to jules_fields_mod. This is the central place where ubiquitous
! fields are distributed and USEd from throughout the model are instantiated.
! Currently only prognostic fields are here but many more are on the way

!type definitions
USE progs_cbl_vars_mod,       ONLY: progs_cbl_vars_type,                       &
                                    progs_cbl_vars_data_type
USE work_vars_mod_cbl,        ONLY: work_vars_data_type,                       &
                                    work_vars_type
!These are dimensionally in per PFT format
USE params_io_mod_cbl,        ONLY: params_io_type,                            &
                                    params_io_data_type

PUBLIC

!TYPES to hold the data
TYPE(progs_cbl_vars_data_type), TARGET :: progs_cbl_vars_data
TYPE(work_vars_data_type),      TARGET :: work_vars_data_cbl
!These are dimensionally in per PFT format
TYPE(params_io_data_type),      TARGET :: pars_io_data_cbl

!TYPES we pass around. These happen to be pointers to the data types above
!but this should be transparent
TYPE(progs_cbl_vars_type) :: progs_cbl_vars
TYPE(work_vars_type)      :: work_vars_cbl
!These are dimensionally in per PFT format
TYPE(params_io_type)      :: pars_io_cbl

END MODULE cable_fields_mod
