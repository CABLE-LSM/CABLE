#!/bin/bash
# This script excludes all the Fortran files from FORD except these that have
# been modified by the Pull Request
# The steps are:
#  1. build a list of all Fortran files under src/ (filenames only)
#  2. build a list of modified Fortran files with git diff and grep
#  3. remove the modified files from the list from 1. using sed
#  4. add "exclude: " to all lines in the list of files from 3. using sed
#  5. Concatenate the generic config options for FORD, the list of excluded
#     files and the front page description for CABLE using cat.
#  6. Clean up by removing unnecessary files

# Get SHA numbers as input
base_SHA=$1
PR_SHA=$2

# List of all source files in the trunk
find src -name "*.F90" | rev | cut -d '/' -f1 | rev > exclude_files_ford.txt

# FORD fails if we exclude everything.
# Remove the cable_driver.F90 from the list so that at least this one is built.
sed -i "/cable_driver.F90/d" exclude_files_ford.txt

# Number of files found
echo number of files found:
wc -l exclude_files_ford.txt

# List of modified files in the last commit
echo "git diff --name-only --diff-filter=ACMRT ${base_SHA} ${PR_SHA} | grep .F90$ | xargs"
modif=$(git diff --name-only --diff-filter=ACMRT ${base_SHA} ${PR_SHA} | grep .F90$ | xargs)

echo modified files found
echo $modif

# Remove each modified file from the list to exclude
for path in $modif
do
    file=$(basename $path)
    sed -i "/$file/d" exclude_files_ford.txt
done

echo number of unchanged files to exclude from FORD:
wc -l exclude_files_ford.txt

# Print the list to exclude in FORD
sed "s/^/exclude: /" exclude_files_ford.txt > exclude.txt

# concatenate the various parts for the FORD config file
cat documentation/cable_config_ford.md exclude.txt documentation/cable_desc_ford.md > documentation/cable_ford.md

# Clean up
rm exclude_files_ford.txt exclude.txt
