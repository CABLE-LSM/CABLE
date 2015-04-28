#/bin/csh
      tail -n +20 offline/cable_write.F90 > junk &&         mv junk offline/cable_write.F90 
cat newheader.txt offline/cable_write.F90 > junk &&         mv junk offline/cable_write.F90 
      tail -n +20 offline/cable_read.F90  > junk &&         mv junk offline/cable_read.F90 
 cat newheader.txt offline/cable_read.F90  > junk &&         mv junk offline/cable_read.F90 
      tail -n +20 offline/cable_output.F90  > junk &&       mv junk offline/cable_output.F90 
 cat newheader.txt offline/cable_output.F90  > junk &&       mv junk offline/cable_output.F90 
      tail -n +20 offline/cable_mpidrv.F90  > junk &&       mv junk offline/cable_mpidrv.F90 
 cat newheader.txt offline/cable_mpidrv.F90  > junk &&       mv junk offline/cable_mpidrv.F90 
      tail -n +20 offline/cable_iovars.F90  > junk &&       mv junk offline/cable_iovars.F90 
 cat newheader.txt offline/cable_iovars.F90  > junk &&       mv junk offline/cable_iovars.F90 
      tail -n +20 offline/cable_input.F90  > junk &&        mv junk offline/cable_input.F90 
 cat newheader.txt offline/cable_input.F90  > junk &&        mv junk offline/cable_input.F90 
      tail -n +20 offline/cable_initialise.F90  > junk &&   mv junk offline/cable_initialise.F90 
 cat newheader.txt offline/cable_initialise.F90  > junk &&   mv junk offline/cable_initialise.F90 
      tail -n +20 offline/cable_driver.F90  > junk &&       mv junk offline/cable_driver.F90 
 cat newheader.txt offline/cable_driver.F90  > junk &&       mv junk offline/cable_driver.F90 
      tail -n +20 offline/cable_checks.F90  > junk &&       mv junk offline/cable_checks.F90 
 cat newheader.txt offline/cable_checks.F90  > junk &&       mv junk offline/cable_checks.F90 
      tail -n +20 offline/cable_abort.F90  > junk &&        mv junk offline/cable_abort.F90 
 cat newheader.txt offline/cable_abort.F90  > junk &&        mv junk offline/cable_abort.F90 
      tail -n +20 core/biogeophys/cable_soilsnow.F90  > junk &&  mv junk core/biogeophys/cable_soilsnow.F90 
 cat newheader.txt core/biogeophys/cable_soilsnow.F90  > junk &&  mv junk core/biogeophys/cable_soilsnow.F90 
      tail -n +20 core/biogeophys/cable_roughness.F90  > junk && mv junk core/biogeophys/cable_roughness.F90 
 cat newheader.txt core/biogeophys/cable_roughness.F90  > junk && mv junk core/biogeophys/cable_roughness.F90 
      tail -n +20 core/biogeophys/cable_radiation.F90  > junk && mv junk core/biogeophys/cable_radiation.F90 
 cat newheader.txt core/biogeophys/cable_radiation.F90  > junk && mv junk core/biogeophys/cable_radiation.F90 
      tail -n +20 core/biogeophys/cable_diag.F90  > junk &&      mv junk core/biogeophys/cable_diag.F90 
 cat newheader.txt core/biogeophys/cable_diag.F90  > junk &&      mv junk core/biogeophys/cable_diag.F90 
      tail -n +20 core/biogeophys/cable_cbm.F90  > junk &&       mv junk core/biogeophys/cable_cbm.F90 
 cat newheader.txt core/biogeophys/cable_cbm.F90  > junk &&       mv junk core/biogeophys/cable_cbm.F90 
      tail -n +20 core/biogeophys/cable_carbon.F90  > junk &&    mv junk core/biogeophys/cable_carbon.F90 
 cat newheader.txt core/biogeophys/cable_carbon.F90  > junk &&    mv junk core/biogeophys/cable_carbon.F90 
      tail -n +20 core/biogeophys/cable_albedo.F90  > junk &&    mv junk core/biogeophys/cable_albedo.F90 
 cat newheader.txt core/biogeophys/cable_albedo.F90  > junk &&    mv junk core/biogeophys/cable_albedo.F90 
      tail -n +20 core/biogeophys/cable_air.F90  > junk &&       mv junk core/biogeophys/cable_air.F90 
 cat newheader.txt core/biogeophys/cable_air.F90  > junk &&       mv junk core/biogeophys/cable_air.F90 
      tail -n +20 UM/casa_um_inout.F90  > junk &&        mv junk UM/casa_um_inout.F90 
 cat newheader.txt UM/casa_um_inout.F90  > junk &&        mv junk UM/casa_um_inout.F90 
      tail -n +20 UM/casa_types.F90  > junk &&           mv junk UM/casa_types.F90 
 cat newheader.txt UM/casa_types.F90  > junk &&           mv junk UM/casa_types.F90 
      tail -n +20 UM/cable_um_tech.F90  > junk &&        mv junk UM/cable_um_tech.F90 
 cat newheader.txt UM/cable_um_tech.F90  > junk &&        mv junk UM/cable_um_tech.F90 
      tail -n +20 UM/cable_um_init.F90  > junk &&        mv junk UM/cable_um_init.F90 
 cat newheader.txt UM/cable_um_init.F90  > junk &&        mv junk UM/cable_um_init.F90 
      tail -n +20 UM/cable_rad_driver.F90  > junk &&     mv junk UM/cable_rad_driver.F90 
 cat newheader.txt UM/cable_rad_driver.F90  > junk &&     mv junk UM/cable_rad_driver.F90 
      tail -n +20 UM/cable_implicit_driver.F90  > junk &&   mv junk UM/cable_implicit_driver.F90 
 cat newheader.txt UM/cable_implicit_driver.F90  > junk &&   mv junk UM/cable_implicit_driver.F90 
      tail -n +20 UM/cable_hyd_driver.F90  > junk &&        mv junk UM/cable_hyd_driver.F90 
 cat newheader.txt UM/cable_hyd_driver.F90  > junk &&        mv junk UM/cable_hyd_driver.F90 
      tail -n +20 UM/cable_explicit_driver.F90  > junk &&   mv junk UM/cable_explicit_driver.F90 
 cat newheader.txt UM/cable_explicit_driver.F90  > junk &&   mv junk UM/cable_explicit_driver.F90 
      tail -n +20 core/biogeochem/casa_inout.F90  > junk && mv junk core/biogeochem/casa_inout.F90 
 cat newheader.txt core/biogeochem/casa_inout.F90  > junk && mv junk core/biogeochem/casa_inout.F90 
      tail -n +20 offline/cable_parameters.F90  > junk &&   mv junk offline/cable_parameters.F90 
 cat newheader.txt offline/cable_parameters.F90  > junk &&   mv junk offline/cable_parameters.F90 
      tail -n +20 offline/cable_mpiworker.F90  > junk &&    mv junk offline/cable_mpiworker.F90 
 cat newheader.txt offline/cable_mpiworker.F90  > junk &&    mv junk offline/cable_mpiworker.F90 
      tail -n +20 offline/cable_mpimaster.F90  > junk &&    mv junk offline/cable_mpimaster.F90 
 cat newheader.txt offline/cable_mpimaster.F90  > junk &&    mv junk offline/cable_mpimaster.F90 
      tail -n +20 offline/cable_mpicommon.F90  > junk &&    mv junk offline/cable_mpicommon.F90 
 cat newheader.txt offline/cable_mpicommon.F90  > junk &&    mv junk offline/cable_mpicommon.F90 
      tail -n +20 core/biogeophys/cable_define_types.F90  > junk &&  mv junk core/biogeophys/cable_define_types.F90 
 cat newheader.txt core/biogeophys/cable_define_types.F90  > junk &&  mv junk core/biogeophys/cable_define_types.F90 
      tail -n +20 core/biogeophys/cable_data.F90 > junk &&           mv junk core/biogeophys/cable_data.F90
 cat newheader.txt core/biogeophys/cable_data.F90 > junk &&           mv junk core/biogeophys/cable_data.F90
      tail -n +20 core/biogeophys/cable_common.F90 > junk &&         mv junk core/biogeophys/cable_common.F90
 cat newheader.txt core/biogeophys/cable_common.F90 > junk &&         mv junk core/biogeophys/cable_common.F90
      tail -n +20 core/biogeophys/cable_canopy.F90 > junk &&         mv junk core/biogeophys/cable_canopy.F90
 cat newheader.txt core/biogeophys/cable_canopy.F90 > junk &&         mv junk core/biogeophys/cable_canopy.F90
      tail -n +20 core/biogeochem/casa_variable.F90 > junk &&        mv junk core/biogeochem/casa_variable.F90
 cat newheader.txt core/biogeochem/casa_variable.F90 > junk &&        mv junk core/biogeochem/casa_variable.F90
      tail -n +20 core/biogeochem/casa_cnp.F90 > junk &&             mv junk core/biogeochem/casa_cnp.F90
 cat newheader.txt core/biogeochem/casa_cnp.F90 > junk &&             mv junk core/biogeochem/casa_cnp.F90
      tail -n +20 core/biogeochem/casa_cable.F90 > junk &&           mv junk core/biogeochem/casa_cable.F90
 cat newheader.txt core/biogeochem/casa_cable.F90 > junk &&           mv junk core/biogeochem/casa_cable.F90
      tail -n +20 UM/cable_um_init_subrs.F90 > junk &&               mv junk UM/cable_um_init_subrs.F90
 cat newheader.txt UM/cable_um_init_subrs.F90 > junk &&               mv junk UM/cable_um_init_subrs.F90

