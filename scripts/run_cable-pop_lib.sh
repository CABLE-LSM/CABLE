#!/usr/bin/env bash

#----------------------------------------------------------------
#
# Library with helper bash functions for run_cable-pop.sh script
#
# Current functions:
#     File / Directory
#     ----------------
#         abspath     absolute path name
#         absfile     file name with absolute path name
#         copyid      copy files adding first argument to filenames
#         cpid        copy files adding first argument to filenames
#         imv         mv file/directory if present
#         irm         rm file/directory if present
#         mvid        rename files adding first argument to filenames
#         renameid    rename files adding first argument to filenames
#         saveid      Copy and move output, restart, log, and namelists
#
#     File content
#     ------------
#         applysed    apply file with "a=b" as sed script
#         csed        make sed command from comma-separated list: var1=str1,var2=str2
#         nckslatlon  return argument to extract lat and lon with ncks
#         ncvarlist   get list of variables in netcdf file
#         text2sed    transform a simple file with "a = b" into sed command file
#
#     String
#     ------
#         findex      index of field in comma-separated list
#         getfield    field in comma-separated list corresponding to name in a header line
#         isin        check if first argument string is present in rest of arguments
#         lower       make string all lowercase
#         tolower     make string all lowercase
#         upper       make string all uppercase
#         toupper     make string all uppercase
#
# Written,  Matthias Cuntz, Mar 2021
#
#----------------------------------------------------------------

#----------------------------------------------------------------
#
# Functions, almost alphabetical
#
#----------------------------------------------------------------

#----------------------------------------------------------------
# absolute path
#   p="../path"
#   p=$(abspath ${p})
function abspath()
{
    idir=${PWD}
    cd ${1}
    odir=${PWD}
    cd ${idir}
    echo "${odir}"
}

#----------------------------------------------------------------
# filename with absolute path
#   f="../path/file"
#   f=$(absfile ${f})
function absfile()
{
    f=$(basename ${1})
    d=$(dirname ${1})
    d=$(abspath ${d})
    echo "${d}/${f}"
}

#----------------------------------------------------------------
# transform a simple file with "a = b" into sed command file
#   echo "inputpath = ../input" >> ini
#   text2sed ini
#   sed -E -f ini musical.nml.in > musical.nml
function text2sed()
{
    ff=${1}
    iblank=${2-1}
    if [[ ${iblank} -eq 1 ]] ; then
        subst='\1 = \2|/'
    else
        subst='\1=\2|/'
    fi
    sed -E \
        -e 's/[[:blank:]]*([[:alnum:]_+.-\%]*)[[:blank:]]*=[[:blank:]]*([^|]*)/s|^[[:blank:]]*\1[[:blank:]]*=.*|'"${subst}" \
        ${ff} > ${tmp}/text2sed.${pid}
    mv ${tmp}/text2sed.${pid} ${ff}
}

#----------------------------------------------------------------
# apply a simple file with "a=b" as sed script, writing new file or
# overwriting input file, and remove simple file
# if $#==3, then new file will be written, if $#==2, then input file will be overwritten
#   echo "inputpath = ../input" >> ini
#   applysed ini musical.nml.ori musical.nml
#   echo "outputpath = ../output" >> out
#   applysed out musical.nml
function applysed()
{
    script=${1}
    if [[ ! -f ${script} ]] ; then
        return
    fi
    text2sed ${script}
    nfile=${2}
    if [[ $# -gt 2 ]] ; then
        ofile=${3}
        if [[ -f ${ofile} || -h ${ofile} ]] ; then rm ${ofile} ; fi
        sed -E -f ${script} ${nfile} > ${ofile}
    else
        sed -E -f ${script} ${nfile} > ${nfile}.$$
        mv ${nfile}.$$ ${nfile}
    fi
    rm ${script}
}

#----------------------------------------------------------------
# copy files adding first argument to filenames
#   copyid ${pid} *.nc
function copyid()
{
    rid=${1}
    shift 1
    for i in "$@" ; do
        if [[ -f ${i} ]] ; then
            cp ${i} ${i%.*}_${rid}.${i##*.}
        fi
    done
}
function cpid()
{
    copyid "$@"
}

#----------------------------------------------------------------
# make sed command from comma-separated list: var1=str1,var2=str2
# no , but = allowed in str
#    com=$(csed "basepath=\"${sitepath}\"")
#    com=${com}$(csed "gridinfo_file=\"${GlobalLandMaskFile}\"")
#    sed ${com} filein > fileout
function csed()
{
    com=""
    for i in $(echo ${1} | tr ',' ' ') ; do
        v=$(echo ${i} | cut -d '=' -f 1)
        s=$(echo ${i} | cut -d '=' -f 2-)
        com="${com} -e s|${v}[[:blank:]]*=.*|${v}=${s}|"
    done
    printf "%s" "${com}"
}

#----------------------------------------------------------------
# get index of field in comma-separated list
# return -1 if not present
#   header="time,co2,ch4"
#   dat="1990.5,350.,1800."
#   ii=$(findex "co2" ${header})
#   if [[ ${ii} -gt 0 ]] ; then
#       out=$(echo ${dat} | cut -d ',' -f ${ii})
#   fi
function findex()
{
    toindx=${1}
    shift
    list="$@"
    out=$(echo ${list} | tr ',' '\n' | grep -n "${toindx}" | cut -f 1 -d ':')
    if [[ -z ${out} ]] ; then out="-1" ; fi
    echo ${out}
}

#----------------------------------------------------------------
# return field in comma-separated list corresponding to name in header line
# returns empty if not present
#   header="time,co2,ch4"
#   dat="1990.5,350.,1800."
#   out=$(getfield "co2" ${header} ${dat})
function getfield()
{
    toget=${1}
    shift
    header=${1}
    shift
    list="$@"
    ii=$(findex ${toget} ${header})
    if [[ ${ii} -gt 0 ]] ; then
        out=$(echo ${list} | cut -d ',' -f ${ii})
    else
        out=""
    fi
    echo ${out}
}

#----------------------------------------------------------------
# mv file/directory if present
# do nothing if target missing
#   imv source target
function imv()
{
    if [[ $# -gt 1 ]] ; then
        for ii in "$@" ; do
            odir=${ii}
        done
        for ii in "$@" ; do
            if [[ "${ii}" != "${odir}" ]] ; then
                if [[ -f ${ii} ]] ; then
                    mv "${ii}" "${odir}"
                fi
            fi
        done
    fi
}

#----------------------------------------------------------------
# rm file/directory if present
#   irm file directory
function irm()
{
    for ii in "$@" ; do
        if [[ -f ${ii} || -h ${ii} ]] ; then
            rm ${ii}
        elif [[ -d ${ii} ]] ; then
            rm -r ${ii}
        fi
    done
}

#----------------------------------------------------------------
# check if first argument string is present in rest of arguments
# returns first argument if present otherwise empty string (check with -z)
#   list="co2 ch4 n2o"
#   if [[ -n $(isin co2 ${list}) ]] ; then echo "Got it" ; fi
#   if [[ -z $(isin co2 ${list}) ]] ; then echo "Need it" ; fi
function isin()
{
    tofind=${1}
    shift
    out=""
    for ii in $@ ; do
        if [[ ${ii} == ${tofind} ]] ; then out=${ii} ; fi
    done
    echo ${out}
}

#----------------------------------------------------------------
# get list of variables in netcdf file with nco (cdo omits the dimension variables)
# returns space-separated list of variable names
#   vars=$(ncvarlist netcdf_file)
function ncvarlist()
{
    out=$(ncks --trd -m ${1} | grep -E ': type' | cut -f 1 -d ' ' | sed -e 's/://' | sort)
    echo ${out}
}

#----------------------------------------------------------------
# return a string that can be used as argument to ncks
# to extract a single lat and lon grid cell
# latlon="48.6742166667,7.06461666667"
# for nc in ${GlobalMetPath}/pres/*.nc ; do
#     ff=$(basename ${nc})
#     ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${MetPath}/pres/${ff}
# done
function nckslatlon()
{
    vars=$(ncvarlist ${1})
    if [[ -z $(isin latitude ${vars}) ]] ; then ilat="lat" ; else ilat="latitude" ; fi
    if [[ -z $(isin longitude ${vars}) ]] ; then ilon="lon" ; else ilon="longitude" ; fi
    if [[ -z $(echo ${2} | cut -f 3 -d ",") || -z $(echo ${2} | cut -f 4 -d ",") ]] ; then
        iilat=$(echo ${2} | cut -f 1 -d ",")
        iilon=$(echo ${2} | cut -f 2 -d ",")
        echo "-d ${ilat},${iilat} -d ${ilon},${iilon}"
    else
        iilat1=$(echo ${2} | cut -f 1 -d ",")
        iilat2=$(echo ${2} | cut -f 2 -d ",")
        iilon1=$(echo ${2} | cut -f 3 -d ",")
        iilon2=$(echo ${2} | cut -f 4 -d ",")
        echo "-d ${ilat},${iilat1},${iilat2} -d ${ilon},${iilon1},${iilon2}"
    fi
}

#----------------------------------------------------------------
# rename files adding first argument to filenames
#   renameid ${pid} *.nc
function renameid()
{
    rid=${1}
    shift 1
    for i in "$@" ; do
        if [[ -f ${i} ]] ; then
            mv ${i} ${i%.*}_${rid}.${i##*.}
        fi
    done
}
function mvid()
{
    renameid "$@"
}

#----------------------------------------------------------------
# copy restart files adding first argument to filenames,
# and move output files, log files, and namelists adding first argument to filenames
#   saveid ${rid} ${mettype} ${doc13o2}
function saveid()
{
    rid=${1}
    mettype=${2}
    doc13o2=${3}
    # mv namelists
    for i in ${mettype} luc cable ; do
	renameid ${rid} ${i}.nml
    done
    mv *_${rid}.nml restart/
    # mv logs
    cd logs
    for i in log log_out ; do
	renameid ${rid} ${i}_cable.txt
    done
    cd ..
    # mv outputs
    cd outputs
    for i in cable casa LUC ; do
	renameid ${rid} ${mettype}_out_${i}.nc
    done
    if [[ ${doc13o2} -eq 1 ]] ; then
	renameid ${rid} ${mettype}_out_casa_c13o2.nc
    fi
    cd ..
    # cp restarts
    cd restart
    for i in climate cable LUC casa_biome casa_met casa_pool casa_phen casa_flux casa_bal ; do
	copyid ${rid} ${mettype}_${i}_rst.nc
    done
    copyid ${rid} pop_${mettype}_ini.nc
    if [[ ${doc13o2} -eq 1 ]] ; then
	for i in flux pools luc ; do
	    copyid ${rid} ${mettype}_c13o2_${i}_rst.nc
	done
    fi
    cd ..
}

#----------------------------------------------------------------
# make string all lowercase
# does not preserve spaces
#   str=$(tolower This String)
function lower()
{
    echo "$@" | tr A-Z a-z
}
function tolower()
{
    echo "$@" | tr A-Z a-z
}

#----------------------------------------------------------------
# make string all uppercase
# does not preserve spaces
#   str=$(toupper This String)
function upper()
{
    echo "$@" | tr a-z A-Z
}
function toupper()
{
    echo "$@" | tr a-z A-Z
}
