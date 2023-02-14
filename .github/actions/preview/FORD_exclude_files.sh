#!/bin/bash
# This script excludes all the Fortran files from FORD except those that have
# been modified by the Pull Request
#
# Inputs:
#  SHA from base branch (usually main)
#  SHA from the merged commit from the PR
#  project code. Accepted values: "cable", "pop"
#
# The steps are:
#  1. build a list of all Fortran files under src/ (filenames only)
#  2. build a list of modified Fortran files with git diff and grep
#  3. remove the modified files from the list from 1. using sed
#  4. add "exclude: " to all lines in the list of files from 3. using sed
#  5. Output 1/0 whether FORD needs to run.
#  5. Concatenate the generic config options for FORD, the list of excluded
#     files and the front page description for CABLE using cat.
#  6. Clean up by removing unnecessary files
#
# Output
#  TO_BUILD: 1 if modified files have been found and FORD needs to run, 0 otherwise.

# Get SHA numbers as input
base_SHA=$1
PR_SHA=$2

# Get the project we need to build for and define the correct files and directories accordingly
proj=$3
echo proj: $proj

if [ $proj == "cable" ]
then
    echo "Looking for modified files in CABLE source src/"
    all_conf_file="cable_ford.md"
    conf_file="src_config/cable_config_ford.md"
    desc_file="src_config/cable_desc_ford.md"
    src_dir="src"
fi

if [ $proj == "pop" ]
then
    echo "Looking for modified files in POP source src_pop/"
    all_conf_file="pop_ford.md"
    conf_file="pop_config/pop_config_ford.md"
    desc_file="pop_config/pop_desc_ford.md"
    src_dir="src_pop"
fi

# List of the basenames of all source file pathnames in the trunk
find ${src_dir} -name "*.F90" | rev | cut -d '/' -f1 | rev > exclude_files_ford.txt

# Number of files found. Need input redirection to get number only
total_files=$(wc -l < exclude_files_ford.txt)
echo total number of files found: $total_files

# List of modified files in the last commit
modif=$(git diff --name-only --diff-filter=ACMRT ${base_SHA} ${PR_SHA} | grep .F90$ | xargs)

echo modified files found
echo $modif

# Remove modified files from the list to exclude
for path in $modif
do
    # Check if the modified file is in the correct source directory
    # Necessary because POP and CABLE have the same set of files
    # and FORD only takes the filename without the path.
    if [ $(echo $path | cut -d/ -f1 ) = $src_dir ]
    then
        file=$(basename $path)
        sed -i "/$file/d" exclude_files_ford.txt
    fi
done

unchanged=$(wc -l < exclude_files_ford.txt)
echo number of unchanged files to exclude from FORD: ${unchanged}

# Output 1 if we need to run FORD for src_dir.
# ie. if there are changed files in the correct source directory
# $modif contains all changed files, but we want the changed files 
# under the correct source directory only.
echo "TO_BUILD=$[ (( ${total_files} - ${unchanged} )) > 0 ]" >> $GITHUB_OUTPUT

# Print the list to exclude in FORD
sed "s/^/exclude: /" exclude_files_ford.txt > exclude.txt

# concatenate the various parts for the FORD config file
cat documentation/${conf_file} exclude.txt documentation/${desc_file} > documentation/${all_conf_file}

# Clean up
rm exclude_files_ford.txt exclude.txt
