MODULE cable_runtime_opts_mod

IMPLICIT NONE

  ! hydraulic_redistribution parameters _soilsnow module
REAL :: wiltParam = 0.0, satuParam = 0.0

! user switches turned on/off by the user thru namelists
! CABLE-2.0 user switches all in single namelist file cable.nml
! clean these up for new namelist(s) format
TYPE kbl_user_switches
  !jhan:make this logical
  CHARACTER(LEN=3) :: diag_soil_resp=''

  CHARACTER(LEN=20) :: fwsoil_switch=''

  ! Ticket #56
  !jhan:options?
  CHARACTER(LEN=20) :: gs_switch='medlyn'

  !INH - new switch for revised coupling on implicit step of ACCESS-CM2 Ticket #132
  LOGICAL :: l_revised_coupling = .FALSE.

  !INH -apply revised sensitvity/correction terms to soilsnow energy balance
  LOGICAL :: l_rev_corr = .FALSE.     !switch to revert to unchanged code

  !ticket#179
  LOGICAL :: soil_thermal_fix = .FALSE.

  !jhan:options?
  CHARACTER(LEN=3) :: ssnow_potev=''
     
     
     
     
     
  !jhan: this is redundant now we all use filename%veg?
  CHARACTER(LEN=200) ::                                                       &
       veg_pars_file  !

  CHARACTER(LEN=20) ::                                                        &
       !H!FWSOIL_SWITCH, &     !
       phenology_switch = 'MODIS'   ! alternative is 'climate'
  !--- LN ------------------------------------------[

  CHARACTER(LEN=10) :: RunIden       = 'STANDARD'  !
  CHARACTER(LEN=6)  :: MetType       = ' ' !
  CHARACTER(LEN=20) :: soil_struc    = "default" ! 'default' or 'sli'
  CHARACTER(LEN=3)  :: POP_out       = 'rst' ! POP output type ('epi' or 'rst')
  CHARACTER(LEN=50) :: POP_rst       = ' ' !
  CHARACTER(LEN=8)  :: casa_out_freq = 'annually' ! 'daily', 'monthly', 'annually'
  CHARACTER(LEN=10)  :: vcmax = 'standard' ! "standard" or "Walker2014"
  CHARACTER(LEN=10)  :: POPLUC_RunType = 'static' ! 'static', 'init', 'restart'

  LOGICAL ::                                                                  &
       call_pop               = .FALSE., & !
       POP_fromZero           = .FALSE.,                                      &
       CALL_Climate           = .FALSE.,                                      &
       Climate_fromZero       = .FALSE.,                                      &
       CASA_fromZero          = .FALSE.,                                      &
       popluc                 = .FALSE.

  INTEGER  ::                                                                 &
       casa_spin_startyear = 1950,                                            &
       casa_spin_endyear   = 1960,                                            &
       yearstart           = 0,                                               &
       yearend             = 0,                                               &
       casa_nrep           = 1
  !--- LN ------------------------------------------]

  CHARACTER(LEN=5) ::                                                         &
       run_diag_level  !

  CHARACTER(LEN=3) ::                                                         &
       !H!DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical)
       leaf_respiration    ! either ON or OFF (jhan:Make Logical)

  ! Custom soil respiration - see Ticket #42
  CHARACTER(LEN=10) ::                                                        &
       smrf_name,   & ! Soil Moist Respiration Function
       strf_name      ! Soil Temp Respiration Function

  LOGICAL ::                                                                  &
       initialize_mapping    = .FALSE., & !
       consistency_check     = .FALSE., & !
       casa_dump_read        = .FALSE., & !
       casa_dump_write       = .FALSE., & !
       cable_runtime_coupled = .TRUE. , & !
       LogWorker             = .TRUE. , & ! Write Output of each worker
                             ! L.Stevens - Test Switches
       l_new_roughness_soil  = .FALSE., & !
       l_new_runoff_speed    = .FALSE., & !
       l_new_reduce_soilevp  = .FALSE., & !

                             ! Switch for customized soil respiration - see Ticket #42
       srf = .FALSE.,                                                         &

                             ! vh_js !
       litter = .FALSE.

  !MD
  LOGICAL :: gw_model = .FALSE.
  LOGICAL :: alt_forcing = .FALSE.

  !using GSWP3 forcing?
  LOGICAL :: gswp3 = .FALSE.
  LOGICAL :: or_evap = .FALSE.
  LOGICAL :: test_new_gw = .FALSE.
  LOGICAL :: sync_nc_file = .FALSE.
  INTEGER :: max_spins = -1
  LOGICAL :: fix_access_roots = .FALSE.  !use pft dependent roots in ACCESS
  !ACCESS roots
  LOGICAL :: access13roots = .FALSE.     !switch to use ACCESS1.3 %froot

  LOGICAL :: l_limit_labile = .FALSE.    ! #237: limit Labile in spinup
  LOGICAL :: NtilesThruMetFile = .FALSE. ! #199: Specify Ntiles thru met file 

END TYPE kbl_user_switches

! instantiate internal switches
TYPE(kbl_user_switches), SAVE :: cable_user



END MODULE cable_runtime_opts_mod

