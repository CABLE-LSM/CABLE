#!/usr/bin/env bash

#PBS -N LUH2_processing
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

# Script to process LUH2 data for CABLE-POP for use in TRENDY.
# Juergen Knauer, May 2024. Gives identical results (+/- numerical precision) 
# compared to Peter Briggs' original scripts.

# The file assumes that LUH2 data are downloaded manually (see TRENDY protocol for location)
# there should be three files: management.nc, states.nc, and transitions.nc.
# management.nc is currently not used

## set up workspace
module purge
module load nco/5.0.5
module load cdo/2.0.5
module load netcdf/4.9.2

#files="management.nc states.nc transitions.nc"
files="states.nc transitions.nc"
filepath="/g/data/x45/LUH2/test"                         # absolute path where files were downloaded to.
outpath="/g/data/x45/LUH2/test/EXTRACT"                  # absolute path for processed files.
gridfile="/g/data/pr09/TRENDY_v12/aux/input.grid.1deg"   # grid file used for TRENDY inputs (excluding everything below 60degS).

startyear=1580
endyear=2023   # note that LUH2 goes one year longer than climate data! (e.g. until 2023 for GCB2023)
#endyear=$(ncap2 -v -O -s 'print(time.size(),"%ld\n");' ${filepath}/states.nc tmptime.nc) # last year of time series = size of time dimension


mkdir -p $outpath
cd $outpath

for file in ${files} ; do
    
    ## 1) select time period of interest, set time axis to years since 1580 and set calendar.
    # latitude and longitude dimensions are renamed to lat and lon (the variables latitude(lat) and longitude(lon) remain). 

    # The transitions file always has one less data record than the others.
    if [[ "${file}" == "transitions.nc" ]] ; then
        endyear_file=$(echo "${endyear} - 1" | bc)
    else 
        endyear_file=${endyear}
    fi

    cp ${filepath}/${file} ${outpath}/${file}
    #ncks -O -d time,${startyear},${endyear_file} ${file} tmp.nc
    cdo -O -f nc4 -z zip5 selyear,${startyear}/${endyear_file} ${file} tmp.nc
    cdo -O -f nc4 -settaxis,${startyear}-01-01,00:00:00,1year -setcalendar,365_day tmp.nc tmp1.nc

    ## 2) Aggregate from 0.25 to 1 degree using conservative remapping
    cdo -O -z zip5 remapcon,${gridfile} tmp1.nc tmp2.nc
    #ncrename -d latitude,lat -d longitude,lon tmp2.nc

    # set optimal netcdf chunking
    # Note: This is now done as part of the ncap2 command below. Chunking can likely be improved!
    #nccopy -d0 -c longitude/64,latitude/64,time/1 tmp2.nc tmp3.nc

    ## 3) Aggregate land use classes into broader categories
    tmpfile="tmp2.nc"

    if [[ "${file}" == "management.nc" ]] ; then

        echo 'no variables extracted from management.nc!'

    elif [[ "${file}" == "states.nc" ]] ; then

        ncap2 -O -s "grass=c3ann+c4ann+c3per+c4per+c3nfx+pastr+range" -v --cnk_map rd1 $tmpfile ${outpath}/grass.nc
        ncap2 -O -s "primaryf=primf+primn"                            -v --cnk_map rd1 $tmpfile ${outpath}/primaryf.nc
        ncap2 -O -s "secondaryf=secdf+secdn"                          -v --cnk_map rd1 $tmpfile ${outpath}/secondaryf.nc

        ncap2 -O -s "crop=c3ann+c4ann+c3per+c4per+c3nfx" -v --cnk_map rd1 $tmpfile ${outpath}/crop.nc
        ncap2 -O -s "past=pastr"                         -v --cnk_map rd1 $tmpfile ${outpath}/past.nc
        ncap2 -O -s "rang=range"                         -v --cnk_map rd1 $tmpfile ${outpath}/rang.nc

    elif [[ "${file}" == "transitions.nc" ]] ; then
        
        ncap2 -O -s "pharv=primf_harv+primn_harv"  -v --cnk_map rd1 $tmpfile ${outpath}/pharv.nc
        ncap2 -O -s "smharv=secmf_harv+secnf_harv" -v --cnk_map rd1 $tmpfile ${outpath}/smharv.nc
        ncap2 -O -s "syharv=secyf_harv"            -v --cnk_map rd1 $tmpfile ${outpath}/syharv.nc

        ncap2 -O -s "ptos=primf_harv+primn_to_secdf+primn_harv" -v --cnk_map rd1 $tmpfile ${outpath}/ptos.nc
        ncap2 -O -s "ptog=primf_to_c3ann+primf_to_c4ann+primf_to_c3per+primf_to_c4per+primf_to_c3nfx+primf_to_pastr+
                          primf_to_range+primn_to_c3ann+primn_to_c4ann+primn_to_c3per+primn_to_c4per+primn_to_c3nfx" -v --cnk_map rd1 $tmpfile ${outpath}/ptog.nc
        ncap2 -O -s "stog=secdf_to_c3ann+secdf_to_c4ann+secdf_to_c3per+secdf_to_c4per+secdf_to_c3nfx+secdf_to_pastr+
                          secdf_to_range+secdn_to_c3ann+secdn_to_c4ann+secdn_to_c3per+secdn_to_c4per+secdn_to_c3nfx" -v --cnk_map rd1 $tmpfile ${outpath}/stog.nc
        ncap2 -O -s "gtos=c3ann_to_secdf+c4ann_to_secdf+c3per_to_secdf+c4per_to_secdf+c3nfx_to_secdf+pastr_to_secdf+
                          range_to_secdf+c3ann_to_secdn+c4ann_to_secdn+c3per_to_secdn+c4per_to_secdn+c3nfx_to_secdn+
                          pastr_to_secdn+range_to_secdn" -v --cnk_map rd1 $tmpfile ${outpath}/gtos.nc
        ncap2 -O -s "ptoc=primf_to_c3ann+primf_to_c4ann+primf_to_c3per+primf_to_c4per+primf_to_c3nfx+primn_to_c3ann+
                          primn_to_c4ann+primn_to_c3per+primn_to_c4per+primn_to_c3nfx" -v --cnk_map rd1 $tmpfile  ${outpath}/ptoc.nc
        ncap2 -O -s "ptoq=primf_to_pastr+primn_to_pastr" -v --cnk_map rd1 $tmpfile  ${outpath}/ptoq.nc
        ncap2 -O -s "stoc=secdf_to_c3ann+secdf_to_c4ann+secdf_to_c3per+secdf_to_c4per+secdf_to_c3nfx+secdn_to_c3ann+
                          secdn_to_c4ann+secdn_to_c3per+secdn_to_c4per+secdn_to_c3nfx" -v $tmpfile --cnk_map rd1 ${outpath}/stoc.nc
        ncap2 -O -s "stoq=secdf_to_pastr + secdn_to_pastr" -v --cnk_map rd1 $tmpfile ${outpath}/stoq.nc
        ncap2 -O -s "ctos=c3ann_to_secdf+c4ann_to_secdf+c3per_to_secdf+c4per_to_secdf+c3nfx_to_secdf+c3ann_to_secdn+
                          c4ann_to_secdn+c3per_to_secdn+c4per_to_secdn+c3nfx_to_secdn" -v --cnk_map rd1 $tmpfile ${outpath}/ctos.nc
        ncap2 -O -s "qtos=pastr_to_secdf+pastr_to_secdn" -v --cnk_map rd1 $tmpfile ${outpath}/qtos.nc

   fi
done


# Give permissions and remove intermediate files.
chmod 775 *.nc
rm tmp.nc tmp*.nc 
rm states.nc transitions.nc