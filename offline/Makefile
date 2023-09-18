#==============================================================================
# This file is part of the
# Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
# This work is licensed under the CABLE Academic User Licence Agreement
# (the "Licence").
# You may not use this file except in compliance with the Licence.
# A copy of the Licence and registration form can be obtained from
# http://www.cawcr.gov.au/projects/access/cable
# You need to register and read the Licence agreement before use.
# Please contact cable_help@nf.nci.org.au for any questions on
# registration and the Licence.
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the Licence is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.
# ==============================================================================
#
# Purpose: builds CABLE executable for offline case
#
# Contact: Bernard.Pak@csiro.au
#
# History: Additions for CASA-CNP

# Additions for SLI (Vanessa Haverd)
# Additions for SpinCASA (Lars Nieradzik)
# Additions for POP, climate and phenology (Vanessa Haverd)
# Additions for Land use change (Vanessa Haverd)
# ==============================================================================
##########################################################################

# these are all the files we are compiling
LSRC =                        \
cable_maths_constants_mod.F90	\
cable_other_constants_mod.F90	\
cable_photo_constants_mod.F90	\
cable_phys_constants_mod.F90	\
grid_constants_cbl.F90				\
cable_climate_type_mod.F90		\
cable_common.F90							\
cable_runtime_opts_mod.F90		\
masks_cbl.F90									\
cbl_snow_albedo.F90						\
cbl_soilColour_albedo.F90			\
cable_canopy.F90							\
cbl_dryLeaf.F90								\
cbl_friction_vel.F90					\
cbl_fwsoil.F90								\
cbl_LAI_canopy_height.F90			\
cbl_latent_heat.F90						\
cbl_photosynthesis.F90				\
cbl_pot_evap_snow.F90					\
cbl_qsat.F90									\
cbl_SurfaceWetness.F90				\
cbl_wetleaf.F90								\
cbl_within_canopy.F90					\
cbl_zetar.F90									\
bgcdriver.F90									\
biogeochem_casa.F90						\
casa_cnp.F90									\
casa_dimension.F90						\
casa_feedback.F90							\
casa_inout.F90								\
casa_param.F90								\
casa_phenology.F90						\
casa_sumcflux.F90							\
casa_variable.F90							\
cable_gw_hydro.F90						\
cable_psm.F90									\
cbl_GW.F90										\
landuse3.F90									\
landuse_constant.F90					\
cable_air.F90									\
cable_carbon.F90							\
cable_climate.F90							\
pop_constants.F90							\
pop_def.F90										\
POPdriver.F90									\
POP.F90												\
pop_io.F90										\
POPLUC.F90										\
pop_types.F90									\
cbl_init_radiation.F90				\
cbl_radiation.F90							\
cbl_rhoch.F90									\
cbl_sinbet.F90								\
cbl_spitter.F90								\
cable_roughness.F90						\
roughnessHGT_effLAI_cbl.F90   \
cable_sli_main.F90						\
cable_sli_numbers.F90					\
cable_sli_roots.F90						\
cable_sli_solve.F90						\
cable_sli_utils.F90						\
cbl_conductivity.F90					\
cbl_GW.F90										\
cbl_hyd_redistrib.F90					\
cbl_Oldconductivity.F90				\
cbl_remove_trans.F90					\
cbl_smoisturev.F90						\
cbl_snowAccum.F90							\
cbl_snow_aging.F90						\
cbl_snowCheck.F90							\
cbl_snowDensity.F90						\
cbl_snowl_adjust.F90					\
cbl_snowMelt.F90							\
cbl_soilfreeze.F90						\
cbl_soilsnow_data.F90					\
cbl_soilsnow_init_special.F90	\
cbl_soilsnow_main.F90					\
cbl_stempv.F90								\
cbl_surfbv.F90								\
cbl_thermal.F90								\
cbl_trimb.F90									\
cable_abort.F90								\
cable_checks.F90							\
cable_cru_TRENDY.F90					\
cable_define_types.F90				\
cable_initialise.F90					\
cable_input.F90								\
cable_iovars.F90							\
cable_LUC_EXPT.F90						\
cable_metutils.F90						\
cable_namelist_input.F90			\
cable_output.F90							\
cable_parameters.F90					\
cable_pft_params.F90					\
cable_phenology.F90						\
cable_plume_mip.F90						\
cable_read.F90								\
cable_site.F90								\
cable_soil_params.F90					\
cable_weathergenerator.F90		\
cable_write.F90								\
casa_cable.F90								\
casa_ncdf.F90									\
CASAONLY_LUC.F90							\
cbl_model_driver_offline.F90	\
landuse_inout.F90							\
spincasacnp.F90							  \
bgcdriver.F90 							  \
biogeochem_casa.F90						\
casa_feedback.F90							\
casa_readbiome.F90						\
casa_rplant.F90							  \
casa_offline_inout.F90				\
casa_sumcflux.F90


# this is the executable we are building
PROG = cable

# from the compiled object files
OBJS	= ${LSRC:.F90=.o}

#compiler switches and flags
CINC = -I$(NCMOD)

#suffixes we use
.SUFFIXES:
.SUFFIXES: .F90 .o
#default rules for these suffixes
.F90.o:
	$(FC) $(CFLAGS) $(CINC) -c $<

# default target by convention is ``all''
all : 

# dependencies, compilation rules for ALL files needed for "all" (LSRC)
#================================================================
#├── params
grid_constants_cbl.o: grid_constants_cbl.F90
cable_maths_constants_mod.o: cable_maths_constants_mod.F90
cable_phys_constants_mod.o: cable_phys_constants_mod.F90
cable_other_constants_mod.o: cable_other_constants_mod.F90 grid_constants_cbl.o
cable_photo_constants_mod.o: cable_photo_constants_mod.F90
cable_runtime_opts_mod.o: cable_runtime_opts_mod.F90
cable_pft_params.o: cable_pft_params.F90 cable_define_types.o
cable_soil_params.o: cable_soil_params.F90 cable_define_types.o
casa_param.o: casa_param.F90

##│   ├── casa-cnp
bgcdriver.o: bgcdriver.F90 biogeochem_casa.o pop_types.o POP.o
biogeochem_casa.o: biogeochem_casa.F90 casa_rplant.o pop_types.o
casa_feedback.o: casa_feedback.F90
casa_readbiome.o: casa_readbiome.F90 casa_variable.o
casa_rplant.o: casa_rplant.F90 casa_cnp.o
casa_sumcflux.o: casa_sumcflux.F90

#var decs/allocs
casa_dimension.o: casa_dimension.F90

#tech files including IO manipulation
masks_cbl.o: masks_cbl.F90
casa_ncdf.o: casa_ncdf.F90 cable_define_types.o cable_common.o
cable_abort.o: cable_abort.F90 cable_iovars.o cable_define_types.o

#science seperated subrs

roughnessHGT_effLAI_cbl.o: roughnessHGT_effLAI_cbl.F90

#science partially seperated subrs
#albedo
cbl_albedo.o: cbl_albedo.F90 cbl_snow_albedo.o cbl_rhoch.o cable_define_types.o cable_common.o
cbl_snow_albedo.o: cbl_snow_albedo.F90 cbl_soilColour_albedo.o cable_phys_constants_mod.o
cbl_soilColour_albedo.o: cbl_soilColour_albedo.F90 
#radiation
cbl_rhoch.o: cbl_rhoch.F90
cbl_spitter.o: cbl_spitter.F90
cbl_sinbet.o: cbl_sinbet.F90 cable_maths_constants_mod.o
cbl_radiation.o: cbl_radiation.F90 cable_define_types.o cable_common.o cbl_init_radiation.o cable_phys_constants_mod.o
cbl_init_radiation.o: cbl_init_radiation.F90 cable_other_constants_mod.o cable_maths_constants_mod.o cable_define_types.o cable_common.o cbl_spitter.o cbl_rhoch.o
#soilsnow
cbl_soilsnow_data.o: cbl_soilsnow_data.F90 cable_phys_constants_mod.o cable_define_types.o cable_common.o 
cbl_trimb.o: cbl_trimb.F90 cbl_soilsnow_data.o 
cbl_soilfreeze.o: cbl_soilfreeze.F90 cbl_soilsnow_data.o cable_phys_constants_mod.o cable_define_types.o cable_common.o 
cbl_snowMelt.o: cbl_snowMelt.F90 cbl_soilsnow_data.o cable_define_types.o
cbl_snowl_adjust.o: cbl_snowl_adjust.F90 cbl_soilsnow_data.o cable_define_types.o
cbl_snowDensity.o: cbl_snowDensity.F90 cbl_soilsnow_data.o cable_define_types.o
cbl_snowCheck.o: cbl_snowCheck.F90 cbl_soilsnow_data.o cable_define_types.o cable_common.o 
cbl_snowAccum.o: cbl_snowAccum.F90 cbl_soilsnow_data.o cable_define_types.o cable_common.o
cbl_snow_aging.o: cbl_snow_aging.F90 cable_phys_constants_mod.o cable_common.o
cbl_remove_trans.o: cbl_remove_trans.F90 cbl_soilsnow_data.o cable_define_types.o cable_common.o 
cbl_hyd_redistrib.o: cbl_hyd_redistrib.F90 cbl_soilsnow_data.o cable_define_types.o cable_common.o 
cbl_conductivity.o: cbl_conductivity.F90 cbl_soilsnow_data.o cable_define_types.o
cbl_Oldconductivity.o: cbl_Oldconductivity.F90 cbl_soilsnow_data.o cable_define_types.o
cbl_GW.o: cbl_GW.F90 cbl_soilsnow_data.o cbl_trimb.o cbl_Oldconductivity.o cbl_conductivity.o cable_define_types.o cable_common.o 
cbl_thermal.o: cbl_thermal.F90 cbl_soilsnow_data.o cbl_snowCheck.o cbl_snowDensity.o cbl_snowAccum.o cbl_snowMelt.o cbl_snowl_adjust.o 
cbl_surfbv.o: cbl_surfbv.F90 cbl_soilsnow_data.o cable_define_types.o cable_common.o cbl_smoisturev.o
cbl_stempv.o: cbl_stempv.F90 cbl_trimb.o cbl_soilsnow_data.o cbl_conductivity.o cbl_Oldconductivity.o
cbl_smoisturev.o: cbl_smoisturev.F90 cbl_trimb.o cbl_soilsnow_data.o cable_define_types.o cable_common.o 
cbl_soilsnow_init_special.o: cbl_soilsnow_init_special.F90 cbl_soilsnow_data.o cable_iovars.o
cbl_soilsnow_main.o: cbl_soilsnow_main.F90 cbl_soilsnow_data.o cbl_soilfreeze.o cbl_snowMelt.o cbl_snowl_adjust.o cbl_snowDensity.o cbl_snowCheck.o cbl_snowAccum.o cbl_remove_trans.o cbl_hyd_redistrib.o cbl_conductivity.o cbl_Oldconductivity.o cbl_GW.o cbl_thermal.o cbl_surfbv.o cbl_stempv.o cbl_smoisturev.o 

#science subrs
casa_phenology.o: casa_phenology.F90

cable_air.o: cable_air.F90 cable_define_types.o cable_common.o

cable_canopy.o: cable_canopy.F90 cbl_friction_vel.o cbl_pot_evap_snow.o cbl_qsat.o cbl_zetar.o cbl_latent_heat.o cbl_wetleaf.o cbl_within_canopy.o cbl_SurfaceWetness.o cbl_dryLeaf.o cbl_fwsoil.o cbl_photosynthesis.o cable_maths_constants_mod.o cable_phys_constants_mod.o cable_other_constants_mod.o cable_photo_constants_mod.o cable_define_types.o cbl_radiation.o cable_roughness.o cable_air.o cable_common.o cable_sli_utils.o cable_sli_main.o cable_gw_hydro.o cable_iovars.o cable_psm.o cable_gw_hydro.o cable_climate_type_mod.o

cbl_friction_vel.o: cbl_friction_vel.F90
cbl_pot_evap_snow.o: cbl_pot_evap_snow.F90 cbl_qsat.o
cbl_qsat.o: cbl_qsat.F90
cbl_zetar.o: cbl_zetar.F90
cbl_latent_heat.o: cbl_latent_heat.F90 cable_define_types.o
cbl_wetleaf.o: cbl_wetleaf.F90 cable_define_types.o
cbl_within_canopy.o: cbl_within_canopy.F90 cable_define_types.o cable_common.o cbl_qsat.o
# added cable_gw_hydro.o to line below - rk4417 - phase2
cbl_SurfaceWetness.o: cbl_SurfaceWetness.F90 cable_define_types.o cable_common.o cable_phys_constants_mod.o cable_gw_hydro.o
cbl_dryLeaf.o: cbl_dryLeaf.F90 cbl_photosynthesis.o cbl_fwsoil.o cable_define_types.o cable_common.o cable_phys_constants_mod.o cable_photo_constants_mod.o 
cbl_fwsoil.o: cbl_fwsoil.F90 cable_define_types.o cable_common.o
cbl_photosynthesis.o: cbl_photosynthesis.F90 cable_define_types.o cable_other_constants_mod.o cable_photo_constants_mod.o 

cable_carbon.o: cable_carbon.F90 cable_define_types.o cable_common.o

cbl_model_driver_offline.o: cbl_model_driver_offline.F90  grid_constants_cbl.o cable_phys_constants_mod.o cable_carbon.o cbl_soilsnow_main.o cable_canopy.o cable_define_types.o cable_roughness.o cbl_radiation.o cbl_albedo.o casa_variable.o cable_common.o cable_gw_hydro.o cable_psm.o masks_cbl.o cbl_snow_aging.o

cable_checks.o: cable_checks.F90 cbl_radiation.o cable_common.o

cable_define_types.o: cable_define_types.F90

cable_iovars.o: cable_iovars.F90 cable_define_types.o

cable_read.o: cable_read.F90 cable_iovars.o cable_abort.o

cable_initialise.o: cable_initialise.F90 cable_read.o cable_iovars.o cable_abort.o  cable_define_types.o

cable_write.o: cable_write.F90 cable_iovars.o cable_abort.o cable_define_types.o casa_variable.o

cable_metutils.o: cable_metutils.F90

cable_input.o: cable_input.F90 pop_types.o POP.o cable_checks.o cable_abort.o cable_read.o cable_initialise.o cbl_sinbet.o \
	       cable_parameters.o cable_common.o cable_LUC_EXPT.o POPLUC.o casa_inout.o cable_metutils.o casa_readbiome.o

cable_output.o: cable_output.F90 cable_abort.o cable_checks.o cable_write.o cable_common.o

cable_parameters.o: cable_parameters.F90 cable_define_types.o cable_iovars.o cable_abort.o casa_variable.o cable_common.o cable_LUC_EXPT.o casa_phenology.o

cable_roughness.o: cable_roughness.F90 cable_phys_constants_mod.o cable_define_types.o cable_common.o roughnessHGT_effLAI_cbl.o cable_other_constants_mod.o


casa_variable.o: casa_variable.F90 POP.o casa_dimension.o casa_param.o

spincasacnp.o: spincasacnp.F90  cable_define_types.o  cable_carbon.o cable_common.o \
	 casa_variable.o casa_phenology.o  casa_cable.o casa_ncdf.o pop_types.o POP.o 

CASAONLY_LUC.o: CASAONLY_LUC.F90 cable_define_types.o  cable_carbon.o cable_common.o casa_variable.o pop_types.o POP.o POPLUC.o cable_LUC_EXPT.o casa_ncdf.o casa_offline_inout.o

casa_cable.o: casa_cable.F90 pop_types.o POP.o cable_define_types.o  cable_carbon.o casa_variable.o \
	      cable_phenology.o casa_cnp.o casa_inout.o casa_ncdf.o biogeochem_casa.o

landuse_constant.o: landuse_constant.F90 cable_define_types.o casa_dimension.o

casa_cnp.o: casa_cnp.F90 cable_define_types.o  casa_variable.o landuse_constant.o

casa_inout.o: casa_inout.F90 POP.o cable_define_types.o  cable_iovars.o casa_variable.o casa_cnp.o casa_ncdf.o  casa_offline_inout.o
casa_offline_inout.o: casa_offline_inout.F90 POP.o cable_define_types.o  cable_iovars.o casa_variable.o casa_cnp.o casa_ncdf.o

cable_common.o: cable_common.F90  cable_iovars.o cable_pft_params.o cable_soil_params.o cable_runtime_opts_mod.o

POP.o: POP.F90 cable_define_types.o pop_def.o pop_constants.o pop_types.o

POPLUC.o: POPLUC.F90  pop_def.o pop_types.o cable_define_types.o cable_iovars.o casa_variable.o cable_common.o cable_runtime_opts_mod.o casa_ncdf.o

cable_LUC_EXPT.o: cable_LUC_EXPT.F90 cable_define_types.o cable_iovars.o casa_variable.o cable_common.o casa_ncdf.o

pop_io.o: pop_io.F90 pop_constants.o pop_types.o casa_variable.o cable_common.o casa_ncdf.o

pop_def.o: pop_def.F90

pop_constants.o: pop_constants.F90 pop_def.o

pop_types.o: pop_types.F90 pop_def.o pop_constants.o

cable_sli_main.o: cable_sli_main.F90  cable_define_types.o cable_sli_numbers.o cable_sli_roots.o cable_sli_solve.o 

cable_sli_numbers.o: cable_sli_numbers.F90  cable_define_types.o

cable_sli_roots.o: cable_sli_roots.F90 cable_sli_numbers.o cable_define_types.o

cable_sli_solve.o: cable_sli_solve.F90 cable_sli_numbers.o cable_sli_roots.o cable_sli_utils.o cable_define_types.o

cable_sli_utils.o: cable_sli_utils.F90 cable_sli_numbers.o  cable_define_types.o

cable_plume_mip.o: cable_plume_mip.F90 cable_iovars.o cable_define_types.o cable_checks.o cable_weathergenerator.o

cable_cru_TRENDY.o: cable_cru_TRENDY.F90 cable_iovars.o cable_define_types.o cable_checks.o \
                    cable_weathergenerator.o casa_variable.o

cable_weathergenerator.o: cable_weathergenerator.F90

cable_climate.o: cable_climate.F90 cable_define_types.o cable_iovars.o  cable_common.o casa_ncdf.o

cable_climate_type_mod.o:  cable_climate_type_mod.F90

cable_phenology.o: cable_phenology.F90 cable_define_types.o cable_iovars.o  cable_common.o  casa_variable.o casa_phenology.o

cable_gw_hydro.o: cable_gw_hydro.F90 cbl_thermal.o cbl_trimb.o cable_iovars.o cable_define_types.o cable_common.o

cable_psm.o: cable_psm.F90 cable_define_types.o cable_common.o cable_air.o 

cable_namelist_input.o: cable_namelist_input.F90

cable_site.o: cable_site.F90 cable_iovars.o cable_common.o


landuse3.o: landuse3.F90 cable_define_types.o casa_dimension.o landuse_constant.o cable_common.o \
            cable_iovars.o casa_variable.o casa_phenology.o

landuse_inout.o: landuse_inout.F90 cable_define_types.o landuse_constant.o landuse3.o cable_common.o cable_iovars.o cable_abort.o

# these are all the files we are compiling
all  : grid_constants_cbl.o    \
  cable_maths_constants_mod.o  \
	cable_phys_constants_mod.o 	 \
	cable_other_constants_mod.o	 \
	cable_photo_constants_mod.o  \
	cable_runtime_opts_mod.o   	 \
	cable_pft_params.o   			 	 \
	cable_soil_params.o        	 \
	masks_cbl.o 								 \
	cable_iovars.o            	 \
	casa_ncdf.o               	 \
	cbl_albedo.o         		 	   \
	cbl_snow_albedo.o					   \
	cbl_soilColour_albedo.o      \
	cbl_rhoch.o         	       \
	cbl_spitter.o       	       \
	cbl_sinbet.o        	       \
	cbl_radiation.o     	       \
	cbl_init_radiation.o	       \
	cable_abort.o   						 \
	cable_parameters.o        	 \
	cable_read.o              	 \
	cable_initialise.o        	 \
	cable_input.o            	   \
	cable_write.o            	   \
	cable_output.o           	   \
	cable_checks.o           	   \
	cable_define_types.o     	   \
	cable_common.o           	   \
	cbl_soilsnow_main.o   	 	   \
	cbl_soilsnow_init_special.o  \
	cbl_soilsnow_init_special.o  \
	cbl_soilfreeze.o             \
	cbl_snowMelt.o               \
	cbl_snowl_adjust.o           \
	cbl_snowDensity.o            \
	cbl_snowCheck.o              \
	cbl_snowAccum.o              \
	cbl_smoisturev.o             \
	cbl_remove_trans.o           \
	cbl_hyd_redistrib.o          \
	cbl_conductivity.o           \
	cbl_Oldconductivity.o        \
	cbl_GW.o                     \
	cbl_thermal.o                \
	cbl_soilsnow_data.o          \
	cbl_trimb.o                  \
	cbl_stempv.o                 \
	cbl_surfbv.o                 \
  cbl_snow_aging.o             \
	cable_air.o               	 \
	cable_roughness.o         	 \
	roughnessHGT_effLAI_cbl.o    \
	cable_carbon.o            	 \
	cable_canopy.o            	 \
	cable_sli_main.o          	 \
	cable_sli_numbers.o       	 \
	cable_sli_roots.o       		 \
	cable_sli_solve.o       		 \
	cable_sli_utils.o       		 \
	POP.o                    	   \
	pop_io.o                 	   \
  POPLUC.o  		 							 \
  pop_constants.o 				     \
  pop_def.o  		 							 \
  pop_types.o  	 							 \
	cable_plume_mip.o        	   \
	cable_weathergenerator.o 	   \
	cable_climate.o          	   \
	cable_climate_type_mod.o     \
	cable_phenology.o        	   \
	cable_cru_TRENDY.o       	   \
	cable_LUC_EXPT.o         	   \
	cable_gw_hydro.o         	   \
	cable_psm.o              	   \
	cable_namelist_input.o   	   \
	cable_site.o             	   \
	casa_variable.o          	   \
	casa_param.o             	   \
	casa_dimension.o         	   \
	casa_phenology.o         	   \
	casa_cnp.o               	   \
	casa_cable.o             	   \
	spincasacnp.o             	 \
	CASAONLY_LUC.o             	 \
	cbl_model_driver_offline.o   \
	casa_inout.o                 \
	casa_offline_inout.o         \
	landuse_constant.o           \
	landuse3.o                   \
	landuse_inout.o              \
  bgcdriver.o   							 \
  biogeochem_casa.o  					 \
  casa_feedback.o  						 \
  casa_readbiome.o  					 \
  casa_rplant.o  							 \
  casa_sumcflux.o  

# make clean option
clean:
	rm -fr .tmp

