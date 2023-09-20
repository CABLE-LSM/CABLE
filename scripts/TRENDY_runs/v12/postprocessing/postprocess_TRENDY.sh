#!/bin/bash

#PBS -N TRENDY_postprocessing3
#PBS -P pr09
#PBS -q normal
#PBS -l walltime=03:00:00
#PBS -l mem=64GB
#PBS -l ncpus=1
#PBS -l storage=gdata/pr09+gdata/x45
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

# Skript to postprocess TRENDY simulations to give outputs ready for upload as specified in the TRENDY protocol
# Juergen Knauer, August 2021
#  - July 2022: TRENDY test suite and ILAMB functionalities added
#  - August 2022: casa variables added (using Matthias' Python script)
#  - August 2023: - script adjusted to handle outputs from serial setup
#				  - additional output variables added

# The script loops through all TRENDY experiments and all variables and does all steps of postprocessing, including
# change of names and units, unit conversion, spatial and temporal aggregation, etc.
# The first step is to copy the TRENDY results from the simdir to the updir. If this step should be repeated (e.g. new runs)
# these files must be deleted first (or the entire folder).
# Assumptions to be met for this script to work:
#       - cable output is at monthly resolution
#       - LUC output is at annual resolution
#       - spatial resolution for all outputs is 1 deg
#       - every grid cell has 3 PFTs (1-2 of them may be missing values)
#       - soil has 6 layers
#       - soil C has 3 pools
#       - 0.7 of wood is belowground and added to roots 

# Note: this script currently is set up for gadi and TRENDY runs S0-S3 only

## TODO:
# - do not process patchfrac for every PFT-type variable but skip if it already exists

### Modules
module load cdo/2.0.5
module load nco/5.0.5
export R_LIBS=/g/data/x45/R/libs
export PYTHONPATH=/home/599/jk8585/.local/lib/python3.10/site-packages  # differs from PYTHONPATH used in main script, harmonise in the future


### ------------------------------------------------------------------------------------------------
### Settings 
### ------------------------------------------------------------------------------------------------
# Paths
modelname="CABLE-POP"
startdir=$PWD
simdir="/g/data/pr09/TRENDY_v12"           # Location of TRENDY outputs
auxdir="/g/data/pr09/TRENDY_v12/aux" 
updir="/g/data/pr09/TRENDY_v12/upload"     # Upload directory
gridfile="${auxdir}/trendy.grid.1deg"
varfile="${startdir}/variables_TRENDY_v12.csv"    # file that contains all information regarding output variables

# Run type, experiment and variable
test_suite=0   # run TRENDY test suite on 1000 grid points (1) or the whole globe (0)? 
run_ilamb=0    # run ILAMB and prepare outputs for it?
#exps="S0 S1 S2 S3"  # experiments to postprocess
exps="S3"    # must end with "_1000pts" if test_suite -eq 1
#vars=$(cut -d, -f 2 ${varfile} | tail -n +2)  # variables to postprocess

#vars="tas mrro evapotrans cVeg cSoil npp npppft gpp gpppft ra rh lai tran nbp"  # Part 1 (+ILAMB)
#vars="pr rsds mrso evapotranspft transpft evapo albedopft snow_depthpft shflxpft rnpft cLitter cProduct cVegpft cSoilpft fLuc rhpft nbppft landCoverFrac oceanCoverFrac laipft cLeaf cWood cRoot cCwd tsl msl evspsblveg evspsblsoi tskinpft mslpft theightpft"  # Part 2
vars="fVegLitter fLeafLitter fWoodLitter fRootLitter fVegSoil fNdep fBNF fNup fNnetmin nSoil nOrgSoil nInorgSoil fNloss nVeg nLeaf nWood nRoot nLitter cSoilpools fAllocLeaf fAllocWood fAllocRoot" # Part 3
#vars="cSoil cSoilpft"
new_run=0    # 1/0: do/do not delete all previous post-processed variables. 
			 # Set to 1 when postprocessing new runs and to 0 when postprocessing new variables of an existing run.
			 # Warning: new_run=1 deletes all existing outputs in the upload folder!


# Model properties
startyear=1700
endyear=2022
layer_depths="0.022 0.058 0.154 0.409 1.085 6.872"   # depths of soil layers (m)
ivegs="1 2 3 4 5 6 7 8 14 17"   # existing ivegs in CABLE
PFT_names="0=Evergreen Needleleaf Forest, 1=Evergreen Broadleaf Forest, 2=Deciduous Needleleaf Forest, 3=Deciduous Broadleaf Forest, 4=Shrub, 5=C3 Grass, 6=C4 Grass, 7=Tundra, 8=Barren, 9=Ice"
SoilCPool_names="0=microbial, 1=slow, 2=passive"

# Output
compress_files=1      # 1/0: do/do not compress output file
strip_history=1       # 1/0: do/do not delete history attributes in netcdf files

## global attributes
creation_date=$(date)
source_code="https://trac.nci.org.au/svn/cable/branches/Share/CABLE-POP_TRENDY"
revision="9636"
institution="CSIRO Environment; Western Sydney University"
contact="Juergen Knauer (J.Knauer@westernsydney.edu.au); Peter Briggs (Peter.Briggs@csiro.au)"

## constants
missval=-99999
gtokg=0.001
kgtoPg=1.0e-12
km2tom2=1.0e6
secperyear=$(echo "86400 * 365" | bc)  # Note: no leap years in cru forcing (also 'leaps=FALSE' in cable.nml)

## pool indices
leaf=0
wood=1
froot=2

## Test suite
if [[ $test_suite -eq 1 ]] ; then
   startyear=1901
   ## specify regions to average over (lon1,lon2,lat1,lat2)
   regions="-180.0,180.0,-60.0,90.0  
            -180.0,180.0,-30.0,30.0  
            -180.0,180.0,30.0,90.0"
fi


## ILAMB                                      
if [[ $run_ilamb -eq 1 ]] ; then
   idir="${simdir}/ILAMB"                             # ILAMB directory
   ivarfile="${startdir}/ILAMB/variables_ilamb.csv"   # file with information on ILAMB variables
   iconfig="${startdir}/ILAMB/ilamb_TRENDY.cfg"       # ILAMB config file

   export PATH=${PATH}:${HOME}/.local/bin
   export MPLBACKEND=Agg 
   export ILAMB_ROOT=${idir}

   mkdir -p ${idir}/RESULTS
   for exp in ${exps} ; do
      mkdir -p ${idir}/MODELS/${exp} 
	  if [[ -d ${idir}/RESULTS/${exp} ]] ; then
         rm -r ${idir}/RESULTS/${exp}
      fi
   done

   ivars=$(cut -d, -f 1 ${ivarfile} | tail -n +2)        # all ILAMB variables (ILAMB name)
   ivars_cable=$(cut -d, -f 2 ${ivarfile} | tail -n +2)  # equivalent name of CABLE variable

   if [[ ! -d ${idir}/DATA ]] ; then
      echo "ILAMB evaluation data need to be downloaded first (ilamb-fetch)"
   fi
fi



### ----------------------------------------------------------------------------------------------
### Start Script
### ----------------------------------------------------------------------------------------------

if [[ $test_suite -eq 1 ]] ; then
   mkdir -p ${updir}/results
else 
   mkdir -p ${updir}/zonal_nbp
fi
    
## loop over experiments
for exp in ${exps} ; do

    echo "Processing experiment ${exp}"

	# delete existing postprocessed outputs for new runs
	if [[ ${new_run} -eq 1 ]] ; then
		rm -r ${updir}/${exp}
	fi

    # Experiment directory
    if [[ $test_suite -eq 1 ]] ; then
        mkdir -p ${updir}/results
        mkdir -p ${updir}/${exp}/results
		cd ${updir}/${exp}/results
    else
		mkdir -p ${updir}/${exp}
		cd ${updir}/${exp}
    fi

    ## process output files and copy to experiment directory within the upload folder: mergetime, set time axis, regrid
	if [[ ! -f cable_output.nc ]] ; then  # only if this hasn't been done before

		## CABLE
		# merge cable output files
		if [[ $test_suite -eq 1 ]] ; then
	    	cp ${simdir}/${exp}/output/cru_out_cable_${startyear}_${endyear}.nc cable_tmp1.nc
		else
	    	cdo -O -s mergetime ${simdir}/${exp}/output/cru_out_cable_${startyear}_1900.nc \
	        	                ${simdir}/${exp}/output/cru_out_cable_1901_${endyear}.nc cable_tmp1.nc
		fi
	    
		# set time axis and calendar
		cdo setcalendar,365_day cable_tmp1.nc cable_tmp.nc 
		cdo settaxis,${startyear}-01-15,12:00:00,1mon cable_tmp.nc cable_tmp1.nc  # ncview will likely display months incorrectly.
		cdo setrtomiss,9e33,inf cable_tmp1.nc cable_tmp.nc

		# regrid output (fill missing values <60degS)
        cdo -f nc4 remapnn,${gridfile} cable_tmp.nc cable_output.nc

		

		## CASA
		if [[ $test_suite -eq 1 ]] ; then
	    	cp ${simdir}/${exp}/output/cru_out_casa_${startyear}_${endyear}.nc casa_tmp1.nc
		else
	    	cdo -O -s mergetime ${simdir}/${exp}/output/cru_out_casa_${startyear}_1900.nc \
	        	                ${simdir}/${exp}/output/cru_out_casa_1901_${endyear}.nc casa_tmp1.nc
		fi
		
		# set time axis and calendar
		cdo setcalendar,365_day casa_tmp1.nc casa_tmp.nc 
		cdo settaxis,${startyear}-01-15,12:00:00,1mon casa_tmp.nc casa_tmp1.nc  # ncview will likely display months incorrectly.
		cdo setrtomiss,9e33,inf casa_tmp1.nc casa_tmp.nc

		# regrid output (fill missing values <60degS)
       cdo -f nc4 remapnn,${gridfile} casa_tmp.nc casa_output.nc



		## LUC for S3 upwards
		if [[ "${exp}" == *"S3"* ]] ; then
	    	if [[ $test_suite -eq 1 ]] ; then		
                cp ${simdir}/${exp}/output/cru_out_LUC_${startyear}_${endyear}.nc LUC_tmp1.nc
	    	else
                cdo -O -s mergetime ${simdir}/${exp}/output/cru_out_LUC_${startyear}_1900.nc \
		    		                ${simdir}/${exp}/output/cru_out_LUC_1901_${endyear}.nc LUC_tmp1.nc
	    	fi
	    	cdo settunits,years LUC_tmp1.nc LUC_tmp.nc
	    	cdo setcalendar,365_day LUC_tmp.nc LUC_tmp1.nc
        	cdo settaxis,${startyear}-06-15,12:00:00,1year LUC_tmp1.nc LUC_tmp.nc
			cdo setrtomiss,9e33,inf LUC_tmp.nc LUC_tmp1.nc

			# regrid output (fill missing values <60degS)
	    	cdo -f nc4 remapnn,${gridfile} LUC_tmp1.nc LUC_output.nc
        fi
    fi


    ## loop over variables
    for var in ${vars} ; do

		echo "Processing variable ${var}"

		## name of output file
		outfile=${modelname}_${exp}_${var}.nc
	
		## get variable characteristics from variable file
		simvar=$(awk -F',' -v var=${var} '{ if ($2 == var) print $3}' $varfile)      # variable name in CABLE
		longname=$(awk -F',' -v var=${var} '{ if ($2 == var) print $1}' $varfile)    # full name as per TRENDY protocol
		filetype=$(awk -F',' -v var=${var} '{ if ($2 == var) print $4}' $varfile)    # output file where variable is located: cable, casa, or LUC
		unit=$(awk -F',' -v var=${var} '{ if ($2 == var) print $5}' $varfile)        # final unit of variable
		scalefac=$(awk -F',' -v var=${var} '{ if ($2 == var) print $7}' $varfile)    # conversion constant (format is operator:number, where operator is one of add,sub,mul,div)
		extradim=$(awk -F',' -v var=${var} '{ if ($2 == var) print $9}' $varfile)    # if empty: longitude, latitude, time
		freq=$(awk -F',' -v var=${var} '{ if ($2 == var) print $11}' $varfile)       # time step: monthly or annual
		varcomment=$(awk -F',' -v var=${var} '{ if ($2 == var) print $13}' $varfile) # comments to be added as variable attribute
	
		echo simvar = $simvar
		echo longname = $longname

    	## Skip variables we don't simulate (a similar check happens further below with $simvar)
		if [[ "${filetype}" == "NA" ]] ; then
	    	echo "No model output for variable ${var}. Continuing with next variable."
	    	continue
		fi

	
		## extract variable from model output or calculate from more than one variable
		if [[ "${simvar}" == "calculate" ]] ; then
	    	case $var in
			mrro)
		    	cdo expr,'mrro=Qs+Qsb;patchfrac;' ${filetype}_output.nc tmp.nc
				;;
			mrso|msl)
				cdo -s selname,SoilMoist ${filetype}_output.nc tmp.nc
				n=0
				for layer_depth in $layer_depths ; do
					ncks -O -d soil,${n} tmp.nc tmp1.nc   # m3 m-3

					n=$(echo "${n} + 1" | bc)

					cdo -s mulc,${layer_depth} tmp1.nc tmp2.nc       # m3 m-2
					cdo -s mulc,1000 tmp2.nc SoilMoist_layer${n}.nc  # l m-2 = kg m-2

					if [[ $n -eq 1 ]] ; then
						cp SoilMoist_layer${n}.nc previous.nc
					else
						cdo -s add SoilMoist_layer${n}.nc previous.nc total_layer${n}.nc
						mv total_layer${n}.nc previous.nc
					fi
				done
				mv previous.nc SoilMoist_cumulative.nc
				rm tmp?.nc

				if [[ "${var}" == "mrso" ]] ; then
					ncrename -O -v SoilMoist,${var} SoilMoist_cumulative.nc tmp.nc
				elif [[ "${var}" == "msl" ]] ; then
					for ((layer=1;layer<=6;layer++)) ; do
						ncpdq -O -a soil,time SoilMoist_layer${layer}.nc SoilMoist_layer${layer}_tmp.nc   # make soil record dimension
					done
					ncrcat -O SoilMoist_layer?_tmp.nc tmp1.nc
				    ncpdq -O -a time,soil tmp1.nc tmp.nc  # reset time as record dimension
					# renaming and attributes
					ncrename -O -v SoilMoist,${var} tmp.nc
					ncrename -O -d soil,smlayer tmp.nc
					ncatted -O -a cell_methods,${var},d,c, tmp.nc
					rm SoilMoist_layer?_tmp.nc
				fi
				rm SoilMoist_layer?.nc SoilMoist_cumulative.nc
				;;
			cProduct)
				if [[ "${exp}" == *"S3"* ]] ; then
					cdo selname,ClearProd,HarvProd ${filetype}_output.nc tmp1.nc
					cdo vertsum tmp1.nc tmp2.nc
					cdo expr,'cProduct=ClearProd+HarvProd;' tmp2.nc tmp1.nc
					cdo mulc,${gtokg} tmp1.nc tmp.nc   # gC m-2 -> kgC m-2
				else
					echo "${var} not calculated for this experiment!"
					continue
				fi
				;;
			cRoot)
				cdo expr,'cRoot=(0.3*PlantCarbWood)+PlantCarbFineRoot;patchfrac;' ${filetype}_output.nc tmp.nc
				;;
			oceanCoverFrac)
				cdo selname,patchfrac ${filetype}_output.nc tmp1.nc  # land mask
				cdo seltimestep,1 tmp1.nc tmp2.nc
				cdo vertsum tmp2.nc tmp1.nc
				cdo setmisstoc,2 tmp1.nc tmp2.nc  # missing data (ocean) to 2
				cdo subc,1 tmp2.nc tmp1.nc        # subtract 1 (land = 0, ocean = 1)
				ncrename -O -v patchfrac,${var} tmp1.nc tmp2.nc
				# remove time information (S3 and higher only)
				if [[ "${exp}" == *"S3"* ]] ; then
					ncwa -O -a time tmp2.nc tmp1.nc
					ncks -O -C -x -v time tmp1.nc tmp.nc
				else
					mv tmp2.nc tmp.nc
				fi
				;;
			landCoverFrac)
				cdo -s selname,patchfrac,iveg ${filetype}_output.nc tmp.nc
				;;
			cSoilpools)
				cdo -s selname,SoilCarbFast,SoilCarbSlow,SoilCarbPassive,patchfrac ${filetype}_output.nc SoilCPools_tmp.nc
				;;
			fAllocRoot)
		        cdo expr,'fAllocRoot=(0.3*CallocWood)+CallocFineRoot;patchfrac;' ${filetype}_output.nc tmp.nc
				cdo selname,iveg cable_output.nc iveg_tmp.nc  # TODO: add iveg to casa output
		        cdo -O merge tmp.nc iveg_tmp.nc tmp1.nc
				cdo -O setrtomiss,-inf,-1 tmp1.nc tmp.nc
				;;
			fVegLitter)
				cdo -s selname,FluxCtolitter ${filetype}_output.nc Ctolitter_tmp.nc
				cdo -s vertsum Ctolitter_tmp.nc tmp.nc
				ncrename -O -v FluxCtolitter,${var} tmp.nc
				;;
			fLeafLitter)
				cdo -s selname,FluxCLeaftolitter ${filetype}_output.nc Cleaftolitter_tmp.nc
				cdo -s vertsum Cleaftolitter_tmp.nc tmp.nc
				ncrename -O -v FluxCLeaftolitter,${var} tmp.nc
				;;
			fWoodLitter)
				cdo -s expr,'fWoodLitter=0.7*FluxCWoodtolitter;' ${filetype}_output.nc tmp1.nc
				cdo -s vertsum tmp1.nc tmp.nc
				;;
			fRootLitter) # note FluxCFineRoottolitte seems to exceed maximum length of variable name and is cropped
			    cdo -s expr,'fRootLitter=0.3*FluxCWoodtolitter+FluxCFineRoottolitte;' ${filetype}_output.nc tmp1.nc
				cdo -s vertsum tmp1.nc tmp.nc
				;;
			fLitterSoil)
				echo "not yet implemented. Requires new output from code"
				continue
				## casaflux%ctosoil includes both fluxes from litter to soil as well as between soil pools. The latter are gross fluxes only!
				## Solution is to report the fromLtoS (not the fromStoS) fluxes only.
				#cdo -s selname,FluxCtosoil ${filetype}_output.nc Ctosoil_tmp.nc
				#cdo -s vertsum Ctosoil_tmp.nc tmp.nc
				#ncrename -O -v FluxCtosoil,${var} tmp.nc
				;;
			fVegSoil)
				cdo -s selname,FluxCtolitter ${filetype}_output.nc Ctosoil_tmp.nc
				cdo -s setrtoc2,0,inf,0.0,0.0 Ctosoil_tmp.nc tmp.nc  # 0.0 everywhere on land
				ncrename -O -v FluxCtolitter,${var} tmp.nc
				;;
			nSoil)
				# organic soil N 
				cdo -s selname,nsoil ${filetype}_output.nc nsoilorg_tmp1.nc
				cdo -s vertsum nsoilorg_tmp1.nc nsoilorg_tmp.nc

				# inorganic soil N
				cdo -s selname,Nsoilmin,patchfrac ${filetype}_output.nc nsoilmin_tmp1.nc
				ncwa -O -a patch -w patchfrac -v Nsoilmin nsoilmin_tmp1.nc nsoilmin_tmp.nc
				cdo -s selname,Nsoilmin nsoilmin_tmp.nc nsoilmin_tmp1.nc

				# total soil N = organic + inorganic soil N
				cdo -s add nsoilorg_tmp.nc nsoilmin_tmp1.nc tmp.nc
				ncrename -O -v nsoil,${var} tmp.nc
				;;
			nOrgSoil)
				cdo -s selname,nsoil ${filetype}_output.nc tmp1.nc
				cdo -s vertsum tmp1.nc tmp.nc
				ncrename -O -v nsoil,${var} tmp.nc
				;;
			nVeg)
				cdo -s selname,nplant ${filetype}_output.nc nveg_tmp.nc
				cdo -s vertsum nveg_tmp.nc tmp.nc
				ncrename -O -v nplant,nVeg tmp.nc
			    ;;
			nLeaf)
				cdo -s selname,nplant ${filetype}_output.nc nveg_tmp.nc
				ncks -O -d mplant,${leaf} nveg_tmp.nc tmp.nc
				ncrename -O -v nplant,nLeaf tmp.nc
				;;
			nWood)
				cdo -s selname,nplant ${filetype}_output.nc nveg_tmp.nc
				ncks -O -d mplant,${wood} nveg_tmp.nc tmp.nc
				ncrename -O -v nplant,nWood tmp.nc
				;;
			nRoot)
				cdo -s selname,nplant ${filetype}_output.nc nveg_tmp.nc
				ncks -O -d mplant,${froot} nveg_tmp.nc nfroot_tmp.nc
				ncks -O -d mplant,${wood} nveg_tmp.nc nwood_tmp1.nc
				cdo -s mulc,0.3 nwood_tmp1.nc nwood_tmp.nc 
				cdo -s add nfroot_tmp.nc nwood_tmp.nc tmp.nc
				ncrename -O -v nplant,nRoot tmp.nc
				;;
			nLitter)
				cdo -s selname,nlitter ${filetype}_output.nc nlitter_tmp.nc
				cdo -s vertsum nlitter_tmp.nc tmp.nc
				ncrename -O -v nlitter,nLitter tmp.nc
			    ;;
			fNloss)
				cdo -s expr,'fNloss=Nminloss+Nminleach;patchfrac;' ${filetype}_output.nc tmp.nc
				ncwa -O -a patch -w patchfrac -v ${var} tmp.nc tmp1.nc
				ncatted -O -a cell_methods,${var},d,c, tmp1.nc tmp.nc
				;;
			*)
				echo "Variable ${var} not yet calculated!"
				continue
				;;
			esac    

		elif [[ "${simvar}" == "" ]] ; then
			if [[ "${var}" == "fLUC" ]] ; then
				echo "${var} is calculated after experiment loop."
			else
				echo "${var} not available from model output! Skip Variable..."
			fi
			continue
		else  # just extract simvar from file, together with patchfrac and iveg (PFT-output only)
			if [[ "${filetype}" == "casa" ]] ; then
				cdo selname,${simvar},patchfrac ${filetype}_output.nc tmp.nc
				if [[ "${extradim}" == "PFT" ]] ; then
					cdo selname,iveg cable_output.nc iveg_tmp.nc  # TODO: add iveg to casa output
					cdo -O merge tmp.nc iveg_tmp.nc tmp1.nc
					cdo -O setrtomiss,-inf,-1 tmp1.nc tmp.nc
				fi
			else
				if [[ "${extradim}" == "PFT" && "${var}" != "landCoverFrac" ]] ; then
					cdo selname,${simvar},patchfrac,iveg ${filetype}_output.nc tmp.nc 
				elif [[ "${extradim}" == "stlayer" ]] ; then
					cdo selname,${simvar} ${filetype}_output.nc tmp.nc 
				else
					cdo selname,${simvar},patchfrac ${filetype}_output.nc tmp.nc
				fi
			fi
		fi	


		## rename variable if necessary
		if [[ "${simvar}" != "${var}" && "${simvar}" != "calculate" ]] ; then    # "calculate" variables are named above
			ncrename -O -v ${simvar},${var} tmp.nc
		fi


        ## patch averaging (if necessary)
		if [[ "${extradim}" != "PFT" ]] ; then
			if [[ "${extradim}" == "SoilCPool" ]] ; then
				ncwa -O -a patch -w patchfrac SoilCPools_tmp.nc SoilCPools_tmp1.nc
				ncks -O -5 SoilCPools_tmp1.nc SoilCPools_tmp2.nc # needed for ncap2 command to work
				# add Soil C Pool dimension
				ncap2 -O -s 'defdim("soilCpool",3);cSoilpools[time,soilCpool,latitude,longitude]=-1e+33' SoilCPools_tmp2.nc SoilCPools_tmp1.nc
				ncap2 -O -s 'cSoilpools(:,0,:,:)=SoilCarbFast;cSoilpools(:,1,:,:)=SoilCarbSlow;cSoilpools(:,2,:,:)=SoilCarbPassive' SoilCPools_tmp1.nc SoilCPools_tmp2.nc
				cdo -s setctomiss,-1e+33 SoilCPools_tmp2.nc SoilCPools_tmp1.nc
				cdo -s selname,cSoilpools SoilCPools_tmp1.nc tmp.nc
			else # TODO: maybe safer to say elif extradim == ""
				if [[ "${var}" != "oceanCoverFrac" && "${var}" != "cProduct" && ("${filetype}" != "casa" || "${simvar}" != "calculate") && "${extradim}" != "stlayer" && "${extradim}" != "smlayer" ]] ; then # already averaged
					ncwa -O -a patch -w patchfrac -v ${var} tmp.nc tmp1.nc
					ncatted -O -a cell_methods,${var},d,c, tmp1.nc tmp.nc
				fi
			fi
		else  # PFT-level output
			for iveg in $ivegs ; do
				for ((patch=0;patch<=2;patch++)) ; do
					# select individual patch and its iveg, patchfrac and variable
					ncks -O -d patch,${patch} tmp.nc tmp1.nc
					cdo -s selname,iveg tmp1.nc iveg_tmp.nc
					cdo -s selname,patchfrac tmp1.nc patchfrac${patch}_tmp.nc
					if [[ "${var}" != "landCoverFrac" ]] ; then
						cdo -s selname,${var} tmp1.nc var_patch${patch}_tmp1.nc
					fi
				
					# create iveg mask
					cdo -s setvrange,$iveg,$iveg iveg_tmp.nc iveg_tmp2.nc

					# use current iveg to mask variable and patchfrac
					cdo -s ifthen iveg_tmp2.nc patchfrac${patch}_tmp.nc patchfrac${patch}_tmp1.nc
					if [[ "${var}" != "landCoverFrac" ]] ; then
						cdo -s ifthen iveg_tmp2.nc var_patch${patch}_tmp1.nc var_patch${patch}_tmp.nc
					fi

					# set missval to 0 for sum to work
					cdo -s setmisstoc,0 patchfrac${patch}_tmp1.nc patchfrac${patch}_tmp.nc

					# sum up patch fraction for that vegtype
					if [[ ${patch} -eq 0 ]] ; then
						cp patchfrac${patch}_tmp.nc previous_patchfrac_tmp.nc
					else
						cdo -s add patchfrac${patch}_tmp.nc previous_patchfrac_tmp.nc total_patchfrac_tmp.nc
						mv total_patchfrac_tmp.nc previous_patchfrac_tmp.nc
					fi

				done  # end patch loop
				cdo -O -s setctomiss,0 previous_patchfrac_tmp.nc total_patchfrac${iveg}_tmp.nc
		
				if [[ "${var}" == "landCoverFrac" ]] ; then
					if [[ "${exp}" == *"S3"* ]] ; then
						ncpdq -O -a patch,time total_patchfrac${iveg}_tmp.nc var_iveg${iveg}_tmp.nc  # make PFT record dimension
						ncrename -O -v patchfrac,landCoverFrac var_iveg${iveg}_tmp.nc
					else
						ncks -O --mk_rec_dmn patch total_patchfrac${iveg}_tmp.nc var_iveg${iveg}_tmp.nc  # make PFT record dimension (would be good to use same command for S0-S2 and S3!)
						ncrename -O -v patchfrac,landCoverFrac var_iveg${iveg}_tmp.nc
					fi
				else  # all other variables
					for ((patch=0;patch<=2;patch++)) ; do
						cdo -s div patchfrac${patch}_tmp.nc total_patchfrac${iveg}_tmp.nc patchfrac_weight${patch}_tmp.nc
						cdo -s mul patchfrac_weight${patch}_tmp.nc var_patch${patch}_tmp.nc var_patch${patch}_weighted_tmp.nc

						if [[ ${patch} -eq 0 ]] ; then
							cp var_patch${patch}_weighted_tmp.nc varsum_previous_tmp.nc
						else
							cdo -s add var_patch${patch}_weighted_tmp.nc varsum_previous_tmp.nc varsum_total_tmp.nc
							mv varsum_total_tmp.nc varsum_previous_tmp.nc
						fi
					done  # end patch loop
					if [[ "${exp}" == *"S3"* ]] ; then
						ncrename -O -v patchfrac,${var} varsum_previous_tmp.nc
					fi
					ncpdq -O -a patch,time varsum_previous_tmp.nc var_iveg${iveg}_tmp.nc  # make PFT record dimension
				fi
			done  # end iveg loop

			## concatenate files along patch dimension
			ncrcat -O var_iveg?_tmp.nc var_iveg??_tmp.nc tmp1.nc

			# reset time as record dimension (if existing)
			if [[ "${var}" == "landCoverFrac" && "${exp}" != *"S3"* ]] ; then
				mv tmp1.nc tmp.nc
			else
				ncpdq -O -a time,patch tmp1.nc tmp.nc 
			fi
			
			# rename patch dimension to PFT 
			ncrename -O -d patch,PFT tmp.nc

			# make sure variable has correct name
			if [[ "${simvar}" != "${var}" && "${var}" != "landCoverFrac" && "${exp}" != *"S3"* ]] ; then 
				ncrename -O -v patchfrac,${var} tmp.nc
			fi
			ncatted -O -a cell_methods,${var},d,c, tmp.nc
		fi

	

	
		## Unit conversion
        operator=$(echo $scalefac | cut -d ":" -f 1)
        convconst=$(echo $scalefac | cut -d ":" -f 2)
	
        case $operator in
			add)
				cdo addc,${convconst} tmp.nc tmp1.nc
				;;
			sub)
				cdo subc,${convconst} tmp.nc tmp1.nc
				;;
			mul)
				cdo mulc,${convconst} tmp.nc tmp1.nc
				;;
			div)
				cdo divc,${convconst} tmp.nc tmp1.nc
				;;
			*)
				echo no unit conversion
				if [[ -f tmp.nc ]] ; then
					mv tmp.nc tmp1.nc
				fi
				;;
		esac

        if [[ -f tmp1.nc ]] ; then
	   		mv tmp1.nc tmp.nc
		fi

	
	
        ## subtract harvest fluxes for NBP in S3 and above
		if [[ "${var}" == "nbp" && "${exp}" == *"S3"* ]] ; then
	    	cdo selname,AgProdLoss,ClearProdLoss,HarvProdLoss LUC_output.nc tmp1.nc  # select fluxes
	    	cdo vertsum tmp1.nc tmp2.nc  # sum over product pools

        	convfac=$(echo "scale=28; ${gtokg} / ${secperyear}" | bc)  # g C m-2 yr-1 -> kg C m-2 s-1
	    	cdo -mulc,${convfac} tmp2.nc tmp1.nc
	    	cdo expr,'AllProdLoss=AgProdLoss+ClearProdLoss+HarvProdLoss;' tmp1.nc ProdLossTotal.nc

	    	for ((year=${startyear};year<=${endyear};year++)) ; do
				cdo -s selyear,${year} ProdLossTotal.nc ProdLossTotalYear.nc # same value every month
		
				cdo -s selyear,${year} tmp.nc NBPyear.nc
				cdo -s sub NBPyear.nc ProdLossTotalYear.nc NBP_${year}.nc
	    	done
	    	cdo -O mergetime NBP_????.nc tmp.nc
	    	rm NBPyear.nc NBP_????.nc ProdLossTotalYear.nc
	    	    
		elif [[ "${var}" == "nbppft" && "${exp}" == *"S3"* ]] ; then
	    	echo "NBP at PFT-level not defined in current CABLE version. Skipping ${var}"
	    	continue
		fi


		## temporal averaging from monthly to annual (always do the averaging to keep it simple)
		if [[ "${filetype}" == "cable" || "${filetype}" == "casa" ]] ; then
	    	if [[ ("${var}" == "landCoverFrac" && "${exp}" != *"S3"*) || ("${var}" == "oceanCoverFrac") ]] ; then
				echo no time aggregation for variable ${var} and experiment ${exp}
				mv tmp.nc tmp_annual.nc
	    	else
				mv tmp.nc tmp_monthly.nc  # copy monthly output for later use
				cdo yearmonmean tmp_monthly.nc tmp1.nc
				# delete redundant 'time_bnds' variable
				ncatted -O -a bounds,time,d,c, tmp1.nc tmp2.nc
				ncks -O -C -x -v time_bnds tmp2.nc tmp1.nc
				# correct time information
				cdo settunits,years tmp1.nc tmp2.nc
				cdo settaxis,${startyear}-06-15,12:00:00,1year tmp2.nc tmp_annual.nc
	    	fi
		elif [[ "${filetype}" == "LUC" ]] ; then
	    	mv tmp.nc tmp_annual.nc
		fi


	     
		## set missing value and set/modify attributes
    	if [[ -f tmp_monthly.nc ]] ; then
	   		ncatted -O -a _FillValue,${var},o,f,${missval} tmp_monthly.nc
			ncatted -O -a missing_value,${var},o,f,${missval} tmp_monthly.nc
	   		ncatted -O -a long_name,${var},o,c,"${longname}" tmp_monthly.nc
	   		ncatted -O -a units,${var},o,c,"${unit}" tmp_monthly.nc	
		fi
		ncatted -O -a _FillValue,${var},o,f,${missval} tmp_annual.nc
		ncatted -O -a missing_value,${var},o,f,${missval} tmp_annual.nc
		ncatted -O -a long_name,${var},o,c,"${longname}" tmp_annual.nc
		ncatted -O -a units,${var},o,c,"${unit}" tmp_annual.nc

	
    	## Prepare standard TRENDY output files
        if [[ $test_suite -eq 0 ]] ; then

	   		# make sure output time step is correct
	   		if [[ "${freq}" == "monthly" ]] ; then
	      		cp tmp_monthly.nc tmp.nc
           	elif [[ "${freq}" == "annual" || "${freq}" == "once" ]] ; then
	      		cp tmp_annual.nc tmp.nc
	   		else
	     		echo "Frequency must be 'monthly', 'yearly', or 'once'!"
	   		fi

	   		# set/modify attributes
	   		if [[ "${var}" == "mrso" || "${var}" == "msl" || "${var}" == "mslpft" || "${var}" == "tsl" ]] ; then
	      		ncatted -O -a layer_depths_m,${var},c,c,"${layer_depths}" tmp.nc
				if [[ "${var}" == "tsl" ]] ; then
					ncrename -O -d soil,stlayer tmp.nc
				fi
	   		fi
	   		if [[ "${varcomment}" != "" ]] ; then
	      		ncatted -O -a comments,${var},c,c,"${varcomment}" tmp.nc
	   		fi
	   		if [[ "${extradim}" == "PFT" ]] ; then
	      		ncatted -O -a PFT_key,${var},c,c,"${PFT_names}" tmp.nc
	   		fi
	   		if [[ "${extradim}" == "SoilCPool" ]] ; then
              	ncatted -O -a SoilCPool_key,${var},c,c,"${SoilCPool_names}" tmp.nc
	   		fi

	   		# global attributes
	   		ncatted -O -a source_code,global,c,c,"${source_code}" tmp.nc
	   		ncatted -O -a code_revision,global,c,c,"${revision}" tmp.nc
           	ncatted -O -a institution,global,c,c,"${institution}" tmp.nc
           	ncatted -O -a contact,global,c,c,"${contact}" tmp.nc

	   		# cleanup: delete redundant attributes
	   		ncatted -O -a Output_averaging,global,d,c, tmp.nc  # confusing --> delete
	   		if [[ "${filetype}" == "LUC" ]] ; then
	      		ncatted -O -a StartYear,global,d,c, tmp.nc
	      		ncatted -O -a EndYear,global,d,c, tmp.nc
	   		fi
			if [[ "${filetype}" == "casa" ]] ; then
				ncatted -O -a icycle,global,d,c, tmp.nc
				ncatted -O -a startyear,global,d,c, tmp.nc
				ncatted -O -a endyear,global,d,c, tmp.nc
				ncatted -O -a runiden,global,d,c, tmp.nc
				ncatted -O -a run-type,global,d,c, tmp.nc
			fi
	   		if [[ $strip_history -eq 1 ]] ; then
	      		ncatted -h -a history,global,d,c, tmp.nc
	      		ncatted -h -a history_of_appended_files,global,d,c, tmp.nc
	   		fi
	   
	   		## name file
	   		mv tmp.nc ${outfile}

	   		## compress file and give permissions
	   		if [[ ${compress_files} -eq 1 ]] ; then
	       		if [[ $test_suite -eq 0 && "${var}" == "nbp" ]] ; then
                	echo "NBP: file is compressed later"
	       		else
		   			gzip -f ${outfile}
	           		chmod 775 ${outfile}.gz
	       		fi
	   		else
	       		chmod 775 $outfile
	   		fi
	     
		else   # $test_suite -eq 1
			### Specific for test suite: calculate mean values across latitudinal bands
        	for region in ${regions} ; do
        		lon1=$(echo $region | cut -d "," -f 1)
             	lon2=$(echo $region | cut -d "," -f 2)
	      		lat1=$(echo $region | cut -d "," -f 3)
	      		lat2=$(echo $region | cut -d "," -f 4)

	      		## annual output
              	cdo sellonlatbox,${lon1},${lon2},${lat1},${lat2} tmp_annual.nc tmp_band.nc
              	cdo fldmean tmp_band.nc ${var}_${lon1}_${lon2}_${lat1}_${lat2}_annual.nc

              	## monthly output
	      		cdo sellonlatbox,${lon1},${lon2},${lat1},${lat2} tmp_monthly.nc tmp_band.nc
              	cdo fldmean tmp_band.nc ${var}_${lon1}_${lon2}_${lat1}_${lat2}_monthly.nc
	  
           	done # end region loop
		fi

		# copy monthly files for use in ilamb
		if [[ $run_ilamb -eq 1 ]] ; then
			if [[ "${ivars_cable}" == *"${var}"* || "${var}" == "rh" || "${var}" == "ra" ]] ; then
			   cp tmp_monthly.nc ${var}_tmp_global.nc
			fi
		fi
	
    done # end variable loop



    
    ### Run ILAMB if run_ilamb = True
    if [[ $run_ilamb -eq 1 ]] ; then

		### 1) process data in preparation for ILAMB	
		for ivar in ${ivars} ; do

	    	var=$(awk -F',' -v ivar=${ivar}             '{ if ($1 == ivar) print $2}' $ivarfile)  # variable name in CABLE
            #istartyear=$(awk -F',' -v ivar=${ivar}       '{ if ($1 == ivar) print $3}' $ivarfile)  # startyear of ILAMB data (not used)
            #iendyear=$(awk -F',' -v ivar=${ivar}         '{ if ($1 == ivar) print $4}' $ivarfile)  # endyear of ILAMB data (not used)
            unit_conversion=$(awk -F',' -v ivar=${ivar} '{ if ($1 == ivar) print $5}' $ivarfile)  # unit conversion factor
	         	
	    	if [[ "${var}" == "derive" ]] ; then
	       		case $ivar in
	          		nee)
		     			cdo add ra_tmp_global.nc rh_tmp_global.nc reco_tmp.nc
		     			cdo sub reco_tmp.nc gpp_tmp_global.nc tmp.nc
		     			ncrename -O -v ra,${ivar} tmp.nc
	             	;;
	          		reco)
                     	cdo add ra_tmp_global.nc rh_tmp_global.nc tmp.nc
		     			ncrename -O -v ra,${ivar} tmp.nc
	             	;;
		  			*)
		     			echo "Nothing implemented here!"
		     		;;
 	       		esac
	    	else
               cp ${var}_tmp_global.nc tmp.nc
	    	fi
	      
	    	# select years and change time units
	    	#cdo selyear,${istartyear}/${iendyear} tmp.nc tmp1.nc
	    	cdo settunits,days tmp.nc tmp1.nc
	    	if [[ ${var} == "lai" ]] ; then
				ncatted -O -a units,lai,o,c,"1" tmp1.nc tmp.nc
	    	else
				mv tmp1.nc tmp.nc
	    	fi
	      
            # convert unit if necessary
	    	if [[ "${unit_conversion}" != "" ]] ; then
				echo "Unit needs conversion but not yet implemented...check!"
	    	fi
	      
	    	# rename variable
	    	if [[ "${var}" != "${ivar}" && "${var}" != "derive" ]] ; then
                ncrename -O -v ${var},${ivar} tmp.nc ${idir}/MODELS/${exp}/CABLE-POP_${exp}_${ivar}.nc
	    	else
				mv tmp.nc ${idir}/MODELS/${exp}/CABLE-POP_${exp}_${ivar}.nc
	    	fi

		done # end ivar loop 

        ### 2) Run ILAMB (configuration at the beginning of the script)
		ilamb-run --config ${iconfig} --model_root ${idir}/MODELS/ --build_dir ${idir}/RESULTS/${exp} \
		 		  --models ${exp} --regions global

    fi

    ## cleanup
    rm *tmp*.nc

    
    ### ---------------------------------------------------------------------------------------------------------------
    ### create Ascii files with annual NBP output for different latitude bands
    ###----------------------------------------------------------------------------------------------------------------
    if [[ $test_suite -eq 0 ]] ; then
    	if [[ "${vars}" == *"nbp"* ]] ; then  # if nbp is in the variable list...
	    
          	## 2.1) get area in km^2 (identical across experiments)
          	cdo selname,Area cable_output.nc area_tmp1.nc
          	cdo vertsum area_tmp1.nc area_tmp2.nc
          	cdo mulc,${km2tom2} area_tmp2.nc area.nc   # area in m^2

         	## 2.2) calculate annual NBP (kgC m-2 s-1)
          	cdo yearmonmean ${modelname}_${exp}_nbp.nc global_nbp_annual_tmp.nc
	  
	  		if [[ ${compress_files} -eq 1 ]] ; then  # compress file
	      		gzip -f ${modelname}_${exp}_nbp.nc
	      		chmod 775 ${modelname}_${exp}_nbp.nc.gz
	  		fi
	  
          	## 2.3) calculate absolute NBP per grid cell (kgC s-1 gridcell-1)
          	cdo mul global_nbp_annual_tmp.nc area.nc global_nbp_tmp.nc

          	## 2.4) convert NBP to PgC yr-1 gridcell-1
          	cdo mulc,${secperyear} global_nbp_tmp.nc global_nbp_tmp1.nc
          	cdo mulc,${kgtoPg} global_nbp_tmp1.nc nbp_global.nc
    
          	## 2.5) extract latitudinal bands
          	# Use: sellonlatbox,lon1,lon2,lat1,lat2
          	cdo sellonlatbox,-180.0,180.0,30.0,90.0    nbp_global.nc  nbp_north.nc
          	cdo sellonlatbox,-180.0,180.0,-30.0,30.0   nbp_global.nc  nbp_tropics.nc
          	cdo sellonlatbox,-180.0,180.0,-90.0,-30.0  nbp_global.nc  nbp_south.nc
    
          	## 2.6) calculate annual values and export to Ascii file (in R)
          	# command line arguments: 1) model name, 2) experiment, 3) startyear, 4) endyear, 5) inpath, 6) outpath
          	${startdir}/create_zonal_nbp_tables.r $modelname $exp $startyear $endyear $PWD ${updir}/zonal_nbp
    
          	## 2.7) cleanup
          	rm *tmp*.nc
			rm nbp_north.nc nbp_south.nc nbp_tropics.nc nbp_global.nc 
       fi
    else   # $test_suite -eq 1
	
       # call R script that summarises all the results in one table
       ${startdir}/create_summary_tables.r $exp "$vars" $startyear $endyear $PWD ${updir}/results

    fi   # check $test_suite

    # TODO: check if another cleanup is necessary here 

done # end experiment loop




### -------------------------------------------------------------------------------------------------------------------------------
### Calculations across experiments
### -------------------------------------------------------------------------------------------------------------------------------
# TODO: check if this works for test_suite -eq 1 (maybe define S2_nbp first depending on $test_suite)
cd ${updir}

# fLUC defined as NBP_S2 - NBP_S3 (note that this gives the flux from land to atmosphere due to LUC as asked for in the script)
if [[ ("${exp}" == *"S3"*) && (-f S2/${modelname}_S2_nbp.nc || -f S2/${modelname}_S2_nbp.nc.gz) && ( "${vars}" == *"fLUC"* ) ]] ; then

   if [[ ${compress_files} -eq 1 ]] ; then
    	gunzip S2/${modelname}_S2_nbp.nc.gz
       	gunzip S3/${modelname}_S3_nbp.nc.gz
      	cdo sub S2/${modelname}_S2_nbp.nc S3/${modelname}_S3_nbp.nc S3/tmp.nc
   fi
 
   ncrename -O -v nbp,fLUC S3/tmp.nc
   ncatted -O -a long_name,fLUC,o,c,"CO2 Flux to Atmosphere from Land Use Change" S3/tmp.nc S3/${modelname}_S3_fLUC.nc
   #ncatted -O -a units,fLUC,o,c,"kg m-2 s-1" S3/tmp.nc S3/${modelname}_S3_fLUC.nc
   if [[ ${compress_files} -eq 1 ]] ; then
       gzip -f S3/${modelname}_S3_fLUC.nc
       gzip -f S2/${modelname}_S2_nbp.nc
       gzip -f S3/${modelname}_S3_nbp.nc
       chmod 775 S3/${modelname}_S3_fLUC.nc.gz
   else
       chmod 775 S3/${modelname}_S3_fLUC.nc
   fi
    
    rm S3/tmp.nc
fi


## End of Script!
