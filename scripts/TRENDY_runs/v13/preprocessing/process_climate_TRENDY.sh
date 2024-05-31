#!/usr/bin/env bash

#PBS -N climate_processing
#PBS -P x45
#PBS -q normal
#PBS -l walltime=01:45:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -l storage=gdata/x45+gdata/pr09+scratch/pr09
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

# Script to process climate forcing data for CABLE-POP for use in TRENDY.
# Juergen Knauer, May/June 2024.

# The script downloads climate data from the Exeter server and processes them.
# The script assumes the following file structure and creates it if not present:
#   basepath/6hourly     <- original data as downloaded from server
#   basepath/daily       <- data aggregated to daily values
#   basepath/daily_1deg  <- daily aggregated to 1 deg (as used by CABLE)

## TODO: revisit aggregation of fd (daymean) and consistency with other tswrf

## set up workspace
module purge
module load nco/5.0.5
module load cdo/2.0.5
module load netcdf/4.9.2

vars="dlwrf fd pre pres spfh tmax tmin tswrf ugrd vgrd"  # climate variables to process
vars="dlwrf"
basepath="/g/data/x45/CRUJRA2022/test"                   # absolute path where forcing files were downloaded to.
gridfile="/g/data/pr09/TRENDY_v12/aux/input.grid.1deg"   # grid file used for TRENDY inputs (excluding everything below 60degS).
cruversion="crujra.v2.2.5d"                              # version name of CRUJRA product to be used
startyear=1901
endyear=2021

mkdir -p ${basepath}/6hourly
mkdir -p ${basepath}/daily
mkdir -p ${basepath}/daily_1deg

for var in ${vars} ; do

    echo "processing variable ${var}"
    
    for year in {${startyear}..${endyear}} ; do
        
        echo "year ${year}"

        # base name of the file
        basename="${cruversion}.${var}.${year}.365d.noc"
        
        # unzip the required file
        gunzip ${basepath}/6hourly/${var}/${basename}.nc.gz

        # set time axis to something sensible
        # TODO: check if this is needed
        # cdo settaxis,$year-01-01,03:00:00,6hours /datasets/work/oa-globalcable/work/CRUJRA2021/6hourly/dlwrf/crujra.v2.2.5d.dlwrf.$year.365d.noc.nc /datasets/work/oa-globalcable/work/CRUJRA2021/6hourly/dlwrf/crujra.v2.2.5d.dlwrf.$year.365d.noc.taxis.nc

        ## aggregate from 6-hourly to daily, change unit attribute if necessary, and aggregate to 1deg
        if [[ "${var}" == "pre" ]] ; then

            # calculate daily sum
            cdo -O daysum ${basepath}/6hourly/${var}/${basename}.nc ${basepath}/daily/${var}/${basename}.daytot.nc
            ncatted -O -a units,pre,m,c,"mm d-1" ${basepath}/daily/${var}/${basename}.daytot.nc
            cdo -O remapcon,${gridfile} ${basepath}/daily/${var}/${basename}.daytot.nc ${basepath}/daily_1deg/${var}/${basename}.daytot.1deg.nc

        elif [[ "${var}" == "tmax" ]] ; then

            # extract daily maximum
            cdo -O -m 9.96921e+36 daymax ${basepath}/6hourly/${var}/${basename}.nc ${basepath}/daily/${var}/${basename}.daymax.nc     
            cdo -O remapcon,${gridfile} ${basepath}/daily/${var}/${basename}.daymax.nc ${basepath}/daily_1deg/${var}/${basename}.daymax.1deg.nc 
        
        elif [[ "${var}" == "tmin" ]] ; then
        
            # extract daily minimum
            cdo -O -m 9.96921e+36 daymin ${basepath}/6hourly/${var}/${basename}.nc ${basepath}/daily/${var}/${basename}.daymin.nc      
            cdo -O remapcon,${gridfile} ${basepath}/daily/${var}/${basename}.daymin.nc ${basepath}/daily_1deg/${var}/${basename}.daymin.1deg.nc 

        else    
        
            # calculate daily mean 
            cdo -O daymean ${basepath}/6hourly/${var}/${basename}.nc ${basepath}/daily/${var}/${basename}.daymean.nc
            cdo -O remapcon,${gridfile} ${basepath}/daily/${var}/${basename}.daymean.nc ${basepath}/daily_1deg/${var}/${basename}.daymean.1deg.nc 

        fi

        ## 3) zip the original 6-hourly file back up again
        gzip ${basepath}/6hourly/${var}/${basename}.nc
    
done
