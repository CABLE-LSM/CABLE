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
LSRC =                          \
 	cable_math_constants_mod.F90  \
 	cable_phys_constants_mod.F90  \
 	cable_other_constants_mod.F90 \
 	cable_photo_constants_mod.F90 \
	cable_runtime_opts_mod.F90 	  \
	cable_pft_params.F90          \
	cable_soil_params.F90         \
	casa_ncdf.F90               	\
	cbl_masks.F90              \
	cbl_hgtAbove_snow.F90     		\
	cbl_lai_eff.F90    		 				\
	cbl_model_driver_offline.F90 \
  cbl_albedo.F90    					  \
	cbl_snow_albedo.F90					  \
  cbl_soilColour_albedo.F90     \
	cbl_rhoch.F90         	      \
	cbl_spitter.F90       	      \
	cbl_sinbet.F90        	      \
	cbl_radiation.F90     	      \
	cbl_init_radiation.F90	      \
	casa_dimension.F90 casa_phenology.F90 casa_param.F90 											 \
	cable_common.F90 cable_define_types.F90 cable_iovars.F90	   \
	cable_abort.F90 casa_variable.F90 cable_parameters.F90                     \
	cbl_soilsnow_main.F90 cbl_soilsnow_subrs.F90 cbl_soilsnow_init_special.F90 \
	cable_air.F90 cable_roughness.F90 																				 \
	cable_carbon.F90 cable_canopy.F90 cable_checks.F90 												 \
	cable_read.F90 cable_initialise.F90 cable_input.F90 cable_write.F90				 \
	cable_output.F90 casa_cnp.F90 casa_cable.F90 casa_inout.F90								 \
	cable_sli_main.F90 cable_sli_numbers.F90 cable_sli_roots.F90							 \
	cable_sli_solve.F90 cable_sli_utils.F90 																	 \
	cable_diag.F90 POP.F90 pop_io.F90 cable_plume_mip.F90 									 	 \
	cable_weathergenerator.F90 cable_climate.F90 cable_phenology.F90 					 \
	cable_cru_TRENDY.F90 POPLUC.F90 cable_LUC_EXPT.F90												 \
	cable_gw_hydro.F90 cable_psm.F90 cable_namelist_input.F90									 \
	cable_site.F90 cable_climate_type_mod.F90 cable_metutils.F90

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

# dependencies
#
# parameters/constants 
cable_math_constants_mod.o: cable_math_constants_mod.F90
cable_phys_constants_mod.o: cable_phys_constants_mod.F90
cable_other_constants_mod.o: cable_other_constants_mod.F90
cable_photo_constants_mod.o: cable_photo_constants_mod.F90
cable_runtime_opts_mod.o: cable_runtime_opts_mod.F90
cable_pft_params.o: cable_pft_params.F90 cable_define_types.o
cable_soil_params.o: cable_soil_params.F90 cable_define_types.o
casa_param.o: casa_param.F90

#var decs/allocs
casa_dimension.o: casa_dimension.F90

#tech files including IO manipulation
cbl_masks.o: cbl_masks.F90
casa_ncdf.o: casa_ncdf.F90 cable_define_types.o cable_common.o
cable_abort.o: cable_abort.F90 cable_iovars.o cable_define_types.o

#tech file diag - exclude from official relllease but part of CABLE for now
cable_fFile.o: cable_fFile.F90
cable_diag.o: cable_diag.F90 cable_common.o cable_fFile.o

#science seperated subrs
cbl_hgtAbove_snow.o: cbl_hgtAbove_snow.F90
cbl_lai_eff.o: cbl_lai_eff.F90

#science partially seperated subrs
#albedo
cbl_albedo.o: cbl_albedo.F90 cbl_snow_albedo.o cbl_rhoch.o cable_define_types.o cable_common.o
cbl_snow_albedo.o: cbl_snow_albedo.F90 cbl_soilColour_albedo.o cable_phys_constants_mod.o
cbl_soilColour_albedo.o: cbl_soilColour_albedo.F90 
#radiation
cbl_rhoch.o: cbl_rhoch.F90
cbl_spitter.o: cbl_spitter.F90
cbl_sinbet.o: cbl_sinbet.F90 cable_math_constants_mod.o
cbl_radiation.o: cbl_radiation.F90 cable_define_types.o cable_common.o cbl_init_radiation.o cable_phys_constants_mod.o
cbl_init_radiation.o: cbl_init_radiation.F90 cable_other_constants_mod.o cable_math_constants_mod.o cable_define_types.o cable_common.o cbl_spitter.o cbl_rhoch.o
#soilsnow
cbl_soilsnow_main.o: cbl_soilsnow_main.F90 cable_iovars.o cable_define_types.o cable_common.o cbl_soilsnow_subrs.o
cbl_soilsnow_subrs.o: cbl_soilsnow_subrs.F90 cable_iovars.o cable_define_types.o cable_common.o
cbl_soilsnow_init_special.o: cbl_soilsnow_init_special.F90 cable_iovars.o cable_define_types.o cable_common.o

#science subrs
casa_phenology.o: casa_phenology.F90

cable_air.o: cable_air.F90 cable_define_types.o cable_common.o

cable_canopy.o: cable_canopy.F90 cable_math_constants_mod.o cable_phys_constants_mod.o cable_other_constants_mod.o cable_photo_constants_mod.o cable_define_types.o cbl_radiation.o cable_roughness.o cable_air.o cable_common.o cable_sli_utils.o cable_sli_main.o cable_gw_hydro.o cable_iovars.o cable_psm.o cable_gw_hydro.o

cable_carbon.o: cable_carbon.F90 cable_define_types.o cable_common.o

cbl_model_driver_offline.o: cbl_model_driver_offline.F90 cable_phys_constants_mod.o cable_carbon.o cbl_soilsnow_main.o cable_canopy.o cable_define_types.o cable_roughness.o cbl_radiation.o cbl_albedo.o casa_variable.o cable_common.o cable_gw_hydro.o cable_psm.o cbl_masks.o

cable_checks.o: cable_checks.F90 cbl_radiation.o cable_common.o

cable_define_types.o: cable_define_types.F90

cable_iovars.o: cable_iovars.F90 cable_define_types.o

cable_read.o: cable_read.F90 cable_iovars.o cable_abort.o

cable_initialise.o: cable_initialise.F90 cable_read.o cable_iovars.o cable_abort.o  cable_define_types.o

cable_write.o: cable_write.F90 cable_iovars.o cable_abort.o cable_define_types.o casa_variable.o

cable_metutils.o: cable_metutils.F90

cable_input.o: cable_input.F90 POP.o cable_checks.o cable_abort.o cable_read.o cable_initialise.o cbl_sinbet.o \
	       cable_parameters.o cable_common.o cable_LUC_EXPT.o POPLUC.o casa_inout.o cable_metutils.o

cable_output.o: cable_output.F90 cable_abort.o cable_checks.o cable_write.o cable_common.o

cable_parameters.o: cable_parameters.F90 cable_define_types.o cable_iovars.o cable_abort.o casa_variable.o cable_common.o cable_LUC_EXPT.o casa_phenology.o

cable_roughness.o: cable_roughness.F90 cable_phys_constants_mod.o cable_define_types.o cable_common.o cbl_hgtAbove_snow.o cbl_lai_eff.o cable_other_constants_mod.o


casa_variable.o: casa_variable.F90 POP.o casa_dimension.o casa_param.o

spincasacnp.o: spincasacnp.F90  cable_define_types.o  cable_carbon.o cable_common.o \
	 casa_variable.o casa_phenology.o  casa_cable.o casa_ncdf.o

#cable_LUC.o: cable_LUC.F90 cable_carbon.o cable_common.o casa_variable.o POP.o POPLUC.o cable_LUC_EXPT.o

CASAONLY_LUC.o: CASAONLY_LUC.F90 cable_define_types.o  cable_carbon.o cable_common.o casa_variable.o POP.o POPLUC.o cable_LUC_EXPT.o casa_ncdf.o

casa_cable.o: casa_cable.F90 POP.o cable_define_types.o  cable_carbon.o casa_variable.o \
	      cable_diag.o cable_phenology.o casa_cnp.o casa_inout.o casa_ncdf.o

casa_cnp.o: casa_cnp.F90 cable_define_types.o  casa_variable.o

casa_inout.o: casa_inout.F90 POP.o cable_define_types.o  cable_iovars.o casa_variable.o casa_cnp.o casa_ncdf.o

cable_common.o: cable_common.F90  cable_iovars.o cable_pft_params.o cable_soil_params.o cable_runtime_opts_mod.o

POP.o: POP.F90 cable_define_types.o

POPLUC.o: POPLUC.F90 cable_define_types.o cable_iovars.o casa_variable.o cable_common.o cable_runtime_opts_mod.o casa_ncdf.o

cable_LUC_EXPT.o: cable_LUC_EXPT.F90 cable_define_types.o cable_iovars.o casa_variable.o cable_common.o casa_ncdf.o

pop_io.o: pop_io.F90 POP.o casa_variable.o cable_common.o casa_ncdf.o

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

cable_gw_hydro.o: cable_gw_hydro.F90 cable_iovars.o cable_define_types.o cable_common.o cbl_soilsnow_subrs.o 

cable_psm.o: cable_psm.F90 cable_define_types.o cable_common.o cable_air.o 

cable_namelist_input.o: cable_namelist_input.F90

cable_site.o: cable_site.F90 cable_iovars.o cable_common.o


# these are all the files we are compiling
all  : cable_math_constants_mod.o   \
	cable_phys_constants_mod.o 	 \
	cable_other_constants_mod.o	 \
	cable_photo_constants_mod.o  \
	cable_runtime_opts_mod.o   	 \
	cable_pft_params.o   			 	 \
	cable_soil_params.o        	 \
	cbl_masks.o 								 \
	cbl_hgtAbove_snow.o					 \
	cbl_lai_eff.o								 \
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
	cbl_soilsnow_subrs.o   	 	   \
	cbl_soilsnow_init_special.o  \
	cable_air.o               	 \
	cable_roughness.o         	 \
	cable_carbon.o            	 \
	cable_canopy.o            	 \
	cable_sli_main.o          	 \
	cable_sli_numbers.o       	 \
	cable_sli_roots.o       		 \
	cable_sli_solve.o       		 \
	cable_sli_utils.o       		 \
	POP.o                    	   \
	pop_io.o                 	   \
	cable_plume_mip.o        	   \
	cable_weathergenerator.o 	   \
	cable_climate.o          	   \
	cable_climate_type_mod.o           \
	cable_phenology.o        	   \
	cable_cru_TRENDY.o       	   \
	POPLUC.o                 	   \
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
	spincasacnp.o             	   \
	CASAONLY_LUC.o             	   \
	cbl_model_driver_offline.o             	   \
	casa_inout.o

# make clean option
clean:
	rm -f *.o *.mod 

