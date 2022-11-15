#!/bin/bash

#PBS -N TRENDY_postprocessing
#PBS -P x45
#PBS -q express
#PBS -l walltime=06:00:00
#PBS -l mem=192GB
#PBS -l ncpus=1
#PBS -l storage=gdata/vl59+gdata/x45
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

# Skript to postprocess TRENDY simulations to give outputs ready for upload as specified in the TRENDY protocol
# Juergen Knauer, August 2021
#  - July 2022: TRENDY test suite and ILAMB functionalities added
#  - August 2022: casa variables added (using Matthias' Python script)

# The script loops through all TRENDY experiments and all variables and does all steps of postprocessing, including
# change of names and units, unit conversion, spatial and temporal aggregation, etc.
# The first step is to copy the TRENDY results from the simdir to the updir. If this step should be repeated (e.g. new runs)
# these files must be deleted first (or the entire folder).
# Assumptions to be met for this script to work:
#       - cable output is at monthly resolution
#       - LUC output is at annual resolution
#       - every grid cell has 3 PFTs (1-2 of them may be missing values)
#       - soil has 6 layers
#       - 3 soil C pools

# Note: this script currently is set up for gadi only and for 1deg resolution and for runs S0-S3 only

## TODO:
# - do not process patchfrac for every PFT-type variable but skip if it already exists
# - check if deep_cleanup is needed and when

### Modules
module load cdo/2.0.5
module load nco/5.0.5


### ------------------------------------------------------------------------------------------------
### Settings 
### ------------------------------------------------------------------------------------------------
# Paths
modelname="CABLE-POP"
startdir=$PWD
simdir="/g/data/vl59/TRENDY_v11"           # Location of TRENDY outputs
auxdir="/g/data/vl59/TRENDY_v11/aux" 
updir="/g/data/vl59/TRENDY_v11/upload"     # Upload directory
gridfile="${auxdir}/trendy.grid.1deg"
varfile="${startdir}/variables_TRENDY_v11.csv"    # file that contains all information regarding output variables

# Run type, experiment and variable
test_suite=0   # run TRENDY test suite on 1000 grid points (1) or the whole globe (0)? 
run_ilamb=0  # run ILAMB and prepare outputs for this?
exps="S0 S1 S2 S3"  # experiments to postprocess
#exps="S0"  # must end with "_1000pts" if test_suite -eq 1
vars=$(cut -d, -f 2 ${varfile} | tail -n +2)  # variables to postprocess
#vars="mrro evapotrans cVeg cSoil gpp ra rh lai tran nbp"
vars="nVeg nLitter nSoil fNdep fBNF fNup fNnetmin fNloss fN2O"
# specify variables from casa output files needed (note: not the same as in vars!)
casavars="kplant,cplant,fracCalloc,Cnpp,ksoil,csoil"         # Carbon
casavars="${casavars},nplant,nlitter,nsoil,Nmindep,Nminfix,Nupland,Nsnet,Nminloss,Nminleach"    # Nitrogen
#casavars="${casavars},clabile,Clabloss,Nsoilmin"  # additional diagnostic casa vars


# Model properties
startyear=1700
endyear=2021
layer_depths="0.022 0.058 0.154 0.409 1.085 6.872"   # depths of soil layers (m)
ivegs="1 2 3 4 5 6 7 8 14 17"   # existing ivegs in CABLE
PFT_names="0=Evergreen Needleleaf Forest, 1=Evergreen Broadleaf Forest, 2=Deciduous Needleleaf Forest, 3=Deciduous Broadleaf Forest, 4=Shrub, 5=C3 Grass, 6=C4 Grass, 7=Tundra, 8=Barren, 9=Ice"
SoilCPool_names="0=microbial, 1=slow, 2=passive"

# Output
compress_files=1      # 1/0: do/do not compress output file
strip_history=1       # 1/0: do/do not delete history attributes in netcdf files
deep_cleanup=0        # 1/0: do/do not delete the cable_output.nc and LUC_output files (set to 1 if final runs) 

## global attributes
creation_date=$(date)
source_code="https://trac.nci.org.au/svn/cable/branches/Share/CABLE-POP_TRENDY"
revision="9072"
institution="CSIRO Climate Science Centre; Western Sydney University"
contact="Juergen Knauer (J.Knauer@westernsydney.edu.au); Peter Briggs (Peter.Briggs@csiro.au)"

## constants
missval=-99999
gtokg=0.001
kgtoPg=1.0e-12
km2tom2=1.0e6
secperyear=$(echo "86400 * 365" | bc)  # Note: no leap years in cru forcing (also 'leaps=FALSE' in cable.nml)


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

   mkdir -p ${idir}/MODELS
   mkdir -p ${idir}/RESULTS

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

		# merge cable output files
		if [[ $test_suite -eq 1 ]] ; then
	    	cp ${simdir}/${exp}/outputs/cru_out_cable_${startyear}_${endyear}.nc cable_tmp1.nc
		else
	    	cdo -O -s mergetime ${simdir}/${exp}/outputs/cru_out_cable_${startyear}_1800.nc ${simdir}/${exp}/outputs/cru_out_cable_1801_1900.nc \
	        	                ${simdir}/${exp}/outputs/cru_out_cable_1901_${endyear}.nc cable_tmp1.nc
		fi
	    
		# set time axis and calendar
		cdo setcalendar,365_day cable_tmp1.nc cable_tmp.nc 
		cdo settaxis,${startyear}-01-15,12:00:00,1mon cable_tmp.nc cable_tmp1.nc  # ncview will likely display months incorrectly.

		# regrid output (fill missing values <60degS)
        cdo -f nc4 remapnn,${gridfile} cable_tmp1.nc cable_output.nc

        ## same for LUC output for S3 upwards
		if [[ "${exp}" == *"S3"* ]] ; then
	    	if [[ $test_suite -eq 1 ]] ; then		
                cp ${simdir}/${exp}/outputs/cru_out_LUC_${startyear}_${endyear}.nc LUC_tmp1.nc
	    	else
                cdo -O -s mergetime ${simdir}/${exp}/outputs/cru_out_LUC_${startyear}_1800.nc ${simdir}/${exp}/outputs/cru_out_LUC_1801_1900.nc \
		    		                ${simdir}/${exp}/outputs/cru_out_LUC_1901_${endyear}.nc LUC_tmp1.nc
	    	fi
	    	cdo settunits,years LUC_tmp1.nc LUC_tmp.nc
	    	cdo setcalendar,365_day LUC_tmp.nc LUC_tmp1.nc
        	cdo settaxis,${startyear}-06-15,12:00:00,1year LUC_tmp1.nc LUC_tmp.nc
	    	cdo -f nc4 remapnn,${gridfile} LUC_tmp.nc LUC_output.nc
        fi
    fi

    # now extract 5dim variables that are needed (SoilMoist and SoilTemp only; cdo can't handle them). Only extract values and merge files.
    # set time axis and regridding needs to happen later.
    if [[ ("${vars}" == *"mrso"* || "${vars}" == *"msl"* || "${vars}" == *"tsl"*) && (! -f SoilMoistTemp_output.nc) ]] ; then
		if [[ $test_suite -eq 1 ]] ; then
	    	ncks -v SoilMoist,SoilTemp,patchfrac,iveg ${simdir}/${exp}/outputs/cru_out_cable_${startyear}_${endyear}.nc SoilMoistTemp_output.nc
        else
	    	ncks -v SoilMoist,SoilTemp,patchfrac,iveg ${simdir}/${exp}/outputs/cru_out_cable_${startyear}_1800.nc  5dimvars_${startyear}_1800.nc
	    	ncks -v SoilMoist,SoilTemp,patchfrac,iveg ${simdir}/${exp}/outputs/cru_out_cable_1801_1900.nc          5dimvars_1801_1900.nc
            ncks -v SoilMoist,SoilTemp,patchfrac,iveg ${simdir}/${exp}/outputs/cru_out_cable_1901_${endyear}.nc    5dimvars_1901_${endyear}.nc
	    	#ncrcat 5dimvars_${startyear}_1900.nc 5dimvars_1901_${endyear}.nc SoilMoistTemp_output.nc
	    	ncrcat 5dimvars_${startyear}_1800.nc 5dimvars_1801_1900.nc 5dimvars_1901_${endyear}.nc SoilMoistTemp_output.nc
	    	rm 5dimvars*
		fi
    fi


    ## casa variables
    # preselect variables in the casa output file, then bring to 2D.
    if [[ ! -f casa_output.nc ]] ; then
    	cdo -s selname,${casavars},patchfrac,latitude,longitude ${simdir}/${exp}/outputs/cru_out_casa_${startyear}_1800.nc  casa_tmp_p1.nc
	    cdo -s selname,${casavars},patchfrac,latitude,longitude ${simdir}/${exp}/outputs/cru_out_casa_1801_1900.nc          casa_tmp_p2.nc
    	cdo -s selname,${casavars},patchfrac,latitude,longitude ${simdir}/${exp}/outputs/cru_out_casa_1901_${endyear}.nc    casa_tmp_p3.nc
	
    	cdo -O -s mergetime casa_tmp_p1.nc casa_tmp_p2.nc casa_tmp_p3.nc casa_tmp.nc

		# convert to 2D (lat,lon)
		# Note that this step does also PFT-averaging! OK for now as none of the casa variables requires PFT-level output. 
		${startdir}/sum_patchcasa2d.py -o casa_tmp1.nc casa_tmp.nc

		# set time axis and calendar
		cdo setcalendar,365_day casa_tmp1.nc casa_tmp.nc 
		cdo settaxis,${startyear}-01-15,12:00:00,1mon casa_tmp.nc casa_tmp1.nc

		# regrid output (fill missing values <60degS)
		# for some reason, a grid has to be set first before remapnn can be applied
		cp ${gridfile} casagrid
		sed -i "s/gridsize  = 64800/gridsize  = 54000/" casagrid
		sed -i "s/ysize     = 180/ysize     = 150/" casagrid
		sed -i "s/yfirst    = -89.5/yfirst    = -59.5/" casagrid
	
		cdo -O setgrid,casagrid casa_tmp1.nc casa_tmp.nc
    	cdo -f nc4 remapnn,${gridfile} casa_tmp.nc casa_output.nc
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
	

    	## Skip variables we don't simulate (a similar check happens further below with $simvar)
		if [[ "${filetype}" == "NA" ]] ; then
	    	echo "No model output for variable ${simvar}. Continuing with next variable."
	    	continue
		fi

	
		## extract variable from model output or calculate from more than one variable
		if [[ "${simvar}" == "calculate" ]] ; then
	    	case $var in
			mrro)
		    	cdo expr,'mrro=Qs+Qsb;patchfrac;' ${filetype}_output.nc tmp.nc
				;;
			mrso|msl|mslpft)
				ncks -O -v patchfrac SoilMoistTemp_output.nc pf_tmp.nc
				cdo -s -f nc4 remapnn,${gridfile} pf_tmp.nc pf_tmp1.nc
				if [[ "${exp}" == *"S3"* ]] ; then
					cdo -s setcalendar,365_day pf_tmp1.nc pf_tmp.nc
					cdo -s settaxis,${startyear}-01-15,12:00:00,1mon pf_tmp.nc patchfrac.nc
				else
					mv pf_tmp1.nc patchfrac.nc
				fi
				
				if [[ ! -f SoilMoist_layer6.nc ]] ; then
					ncks -O -v SoilMoist SoilMoistTemp_output.nc SoilMoist.nc
					n=0
					for layer_depth in $layer_depths ; do
						ncks -O -d soil,${n} SoilMoist.nc tmp1.nc   # m3 m-3
						ncwa -O -a soil tmp1.nc tmp2.nc             # needs high memory

						# set time axis and regrid here (since it's 4dim now)
						cdo -s setcalendar,365_day tmp2.nc tmp1.nc
						cdo -s settaxis,${startyear}-01-15,12:00:00,1mon tmp1.nc tmp2.nc
						cdo -s -f nc4 remapnn,${gridfile} tmp2.nc tmp1.nc

						n=$(echo "${n} + 1" | bc)

						cdo -s mulc,${layer_depth} tmp1.nc tmp2.nc       # m3 m-2
						cdo -s mulc,1000 tmp2.nc SoilMoist_layer${n}_tmp.nc  # l m-2 = kg m-2

						if [[ $n -eq 1 ]] ; then
							cp SoilMoist_layer${n}_tmp.nc previous.nc
						else
							cdo -s add SoilMoist_layer${n}_tmp.nc previous.nc total_layer${n}.nc
							mv total_layer${n}.nc previous.nc
						fi
						ncks -A patchfrac.nc SoilMoist_layer${n}_tmp.nc
						cdo -s settaxis,${startyear}-01-15,12:00:00,1mon SoilMoist_layer${n}_tmp.nc SoilMoist_layer${n}.nc
					done
					mv previous.nc SoilMoist_cumulative.nc
					rm SoilMoist.nc SoilMoist_layer?_tmp.nc tmp?.nc
				fi 

				if [[ "${var}" == "mrso" || "${var}" == "mslpft" ]] ; then
					ncks -A patchfrac.nc SoilMoist_cumulative.nc  # append patchfrac
					if [[ "${var}" == "mslpft" ]] ; then
						ncks -O -v iveg SoilMoistTemp_output.nc iveg_tmp.nc
						cdo -s -f nc4 remapnn,${gridfile} iveg_tmp.nc iveg.nc
						ncks -A iveg.nc SoilMoist_cumulative.nc   # append iveg
					fi
					cdo -s settaxis,${startyear}-01-15,12:00:00,1mon SoilMoist_cumulative.nc tmp.nc   # PFT averaging happens below
					ncrename -O -v SoilMoist,${var} tmp.nc
				elif [[ "${var}" == "msl" ]] ; then
					for ((layer=1;layer<=6;layer++)) ; do
						ncrename -O -v SoilMoist,${var} SoilMoist_layer${layer}.nc
					done
				fi
				;;
			tsl)
				ncks -O -v SoilTemp SoilMoistTemp_output.nc SoilTemp.nc
				ncks -O -v patchfrac SoilMoistTemp_output.nc pf_tmp.nc
				cdo -s -f nc4 remapnn,${gridfile} pf_tmp.nc pf_tmp1.nc
				if [[ "${exp}" == *"S3"* ]] ; then
					cdo -s setcalendar,365_day pf_tmp1.nc pf_tmp.nc
					cdo -s settaxis,${startyear}-01-15,12:00:00,1mon pf_tmp.nc patchfrac.nc
				else
					mv pf_tmp1.nc patchfrac.nc
				fi
				ncrename -O -v SoilTemp,${var} SoilTemp.nc
				
				n=0
				for layer_depth in $layer_depths ; do
					ncks -O -d soil,${n} SoilTemp.nc tmp1.nc
					ncwa -O -a soil tmp1.nc tmp2.nc
				
					# set time axis and regrid here (since it's 4dim now)
					cdo -s setcalendar,365_day tmp2.nc tmp1.nc
					#cdo settaxis,${startyear}-01-15,12:00:00,1mon tmp1.nc tmp2.nc
					cdo -s -f nc4 remapnn,${gridfile} tmp1.nc tmp2.nc

					n=$(echo "${n} + 1" | bc)
					ncks -A patchfrac.nc tmp2.nc
					cdo -s settaxis,${startyear}-01-15,12:00:00,1mon tmp2.nc SoilTemp_layer${n}.nc
					rm tmp2.nc
				done
				rm tmp1.nc SoilTemp.nc pf_tmp.nc
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
			nVeg)
				cdo -s selname,nplant ${filetype}_output.nc nveg_tmp.nc
				cdo -s vertsum nveg_tmp.nc tmp.nc
				ncrename -O -v nplant,nVeg tmp.nc
			    ;;
			nLitter)
				cdo -s selname,nlitter ${filetype}_output.nc nlitter_tmp.nc
				cdo -s vertsum nlitter_tmp.nc tmp.nc
				ncrename -O -v nlitter,nLitter tmp.nc
			    ;;
			nSoil)
				cdo -s selname,nsoil ${filetype}_output.nc nsoil_tmp.nc
				cdo -s vertsum nsoil_tmp.nc tmp.nc
				ncrename -O -v nsoil,nSoil tmp.nc
			    ;;		
			*)
				echo "Variable ${simvar} not yet calculated!"
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
				cdo selname,${simvar} ${filetype}_output.nc tmp.nc
			else
				if [[ "${extradim}" == "PFT" && "${var}" != "landCoverFrac" ]] ; then
					cdo selname,${simvar},patchfrac,iveg ${filetype}_output.nc tmp.nc 
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
			# notes ncwa: outfile must be present; patchfrac variable automatically deleted
			if [[ "${extradim}" == "smlayer" ]] ; then
				if [[ "${var}" != "msl" ]] ; then
					echo "WARNING: 'smlayer' as additional dimension only defined for variable 'msl', not for ${var}...Rewrite this code!"
				fi
				n=1
				for layer_depth in $layer_depths ; do
					ncwa -O -a patch -w patchfrac -v ${var} SoilMoist_layer${n}.nc tmp1.nc
					ncap2 -O -s 'defdim("soil",1);msl[time,soil,latitude,longitude]=msl' tmp1.nc tmp2.nc  # add back soil dimension
					ncpdq -O -a soil,time tmp2.nc SoilMoist_layer${n}_pftavg.nc  # make soil record dimension
					n=$(echo "${n} + 1" | bc)
				done
				## concatenate files along soil layer dimension
				ncrcat -O SoilMoist_layer?_pftavg.nc tmp1.nc
				ncpdq -O -a time,soil tmp1.nc tmp.nc  # reset time as record dimension
			elif [[ "${extradim}" == "stlayer" ]] ; then
				if [[ "${var}" != "tsl" ]] ; then
					echo "WARNING: 'stlayer' as additional dimension only defined for variable 'tsl', not for ${var}...Rewrite this code!"
				fi
				n=1
				for layer_depth in $layer_depths ; do
					ncwa -O -a patch -w patchfrac -v ${var} SoilTemp_layer${n}.nc tmp1.nc
					ncap2 -O -s 'defdim("soil",1);tsl[time,soil,latitude,longitude]=tsl' tmp1.nc tmp2.nc  # add back soil dimension
					ncpdq -O -a soil,time tmp2.nc SoilTemp_layer${n}_pftavg.nc  # make soil record dimension
					n=$(echo "${n} + 1" | bc)
				done
				## concatenate files along soil layer dimension
				ncrcat -O SoilTemp_layer?_pftavg.nc tmp1.nc
				ncpdq -O -a time,soil tmp1.nc tmp.nc  # reset time as record dimension
			elif [[ "${extradim}" == "SoilCPool" ]] ; then
				ncwa -O -a patch -w patchfrac SoilCPools_tmp.nc SoilCPools_tmp1.nc
				ncks -O -5 SoilCPools_tmp1.nc SoilCPools_tmp2.nc # needed for ncap2 command to work
				# add Soil C Pool dimension
				ncap2 -O -s 'defdim("soilCpool",3);cSoilpools[time,soilCpool,latitude,longitude]=-1e+33' SoilCPools_tmp2.nc SoilCPools_tmp1.nc
				ncap2 -O -s 'cSoilpools(:,0,:,:)=SoilCarbFast;cSoilpools(:,1,:,:)=SoilCarbSlow;cSoilpools(:,2,:,:)=SoilCarbPassive' SoilCPools_tmp1.nc SoilCPools_tmp2.nc
				cdo -s setctomiss,-1e+33 SoilCPools_tmp2.nc SoilCPools_tmp1.nc
				cdo -s selname,cSoilpools SoilCPools_tmp1.nc tmp.nc
			else # TODO: maybe safer to say elif extradim == ""
				if [[ "${var}" != "oceanCoverFrac" && "${var}" != "cProduct" && "${filetype}" != "casa" ]] ; then # already averaged
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
			if [[ "${simvar}" != "${var}" && "${var}" != "landCoverFrac" ]] ; then 
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
	   		ncatted -O -a _FillValue,${var},d,f, tmp_monthly.nc
	   		ncatted -O -a long_name,${var},o,c,"${longname}" tmp_monthly.nc
	   		ncatted -O -a units,${var},o,c,"${unit}" tmp_monthly.nc	
	   		ncatted -O -a missing_value,${var},o,f,${missval} tmp_monthly.nc
		fi
		ncatted -O -a _FillValue,${var},o,f,${missval} tmp_annual.nc
		ncatted -O -a _FillValue,${var},d,f, tmp_annual.nc  # TODO: check if this is necessary
		ncatted -O -a long_name,${var},o,c,"${longname}" tmp_annual.nc
		ncatted -O -a units,${var},o,c,"${unit}" tmp_annual.nc
		ncatted -O -a missing_value,${var},o,f,${missval} tmp_annual.nc

	
	
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



	   		# rename some variables
           	if [[ "${var}" == "msl" ]] ; then
	      		for ((layer=1;layer<=6;layer++)) ; do
	         		ncrename -O -v ${var},SoilMoist SoilMoist_layer${layer}.nc
	      		done
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
	     
		fi  # end if $test_suite -eq 0



	   
		### Specific for test suite: calculate mean values across latitudinal bands
		if [[ $test_suite -eq 1 || $run_ilamb -eq 1 ]] ; then
        	for region in ${regions} ; do
        		lon1=$(echo $region | cut -d "," -f 1)
             	lon2=$(echo $region | cut -d "," -f 2)
	      		lat1=$(echo $region | cut -d "," -f 3)
	      		lat2=$(echo $region | cut -d "," -f 4)

	      		## annual output
              	cdo sellonlatbox,${lon1},${lon2},${lat1},${lat2} tmp_annual.nc tmp_band.nc
              	cdo fldmean tmp_band.nc ${var}_${lon1}_${lon2}_${lat1}_${lat2}_annual.nc

              	## monthly output
	      		if [[ $run_ilamb -eq 1 && "${region}" == "-180.0,180.0,-60.0,90.0" ]] ; then # retain 2D fields for ILAMB if necessary
	         		toutfile="${var}_tmp_global.nc"
	      		else
	         		toutfile="tmp_band.nc"
	      		fi
	      		cdo sellonlatbox,${lon1},${lon2},${lat1},${lat2} tmp_monthly.nc ${toutfile}
              	cdo fldmean ${toutfile} ${var}_${lon1}_${lon2}_${lat1}_${lat2}_monthly.nc
	  
           	done # end region loop
	   
		fi
	
    done # end variable loop



    
    ### Run ILAMB if run_ilamb = True
    if [[ $run_ilamb -eq 1 ]] ; then

		### 1) process data in preparation for ILAMB
		glob_latlon=${lon1}_${lon2}_${lat1}_${lat2}
        ivars=$(cut -d, -f 1 ${ivarfile} | tail -n +2)  # all ILAMB variables
	
		for ivar in ${ivars} ; do

	    	var=$(awk -F',' -v ivar=${ivar}             '{ if ($1 == ivar) print $2}' $ivarfile)  # variable name in CABLE
            #istartyear=$(awk -F',' -v ivar=${ivar}       '{ if ($1 == ivar) print $3}' $ivarfile)  # startyear of ILAMB data (not used)
            #iendyear=$(awk -F',' -v ivar=${ivar}         '{ if ($1 == ivar) print $4}' $ivarfile)  # endyear of ILAMB data (not used)
            unit_conversion=$(awk -F',' -v ivar=${ivar} '{ if ($1 == ivar) print $5}' $ivarfile)  # unit conversion factor

	    	mkdir -p ${idir}/MODELS/${exp}  # TODO: could do this earlier

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

		if [[ -d ${idir}/RESULTS/${exp} ]] ; then
            rm -r ${idir}/RESULTS/${exp}
        fi

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
          	cdo mulc,${km2tom2} area_tmp2.nc area.nc   # m^2

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
          	#rm cable_output.nc LUC_output.nc
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
