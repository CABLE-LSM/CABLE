#!/bin/bash

# List of all source files in the trunk
find src -name "*.F90" | rev | cut -d '/' -f1 | rev > exclude_files_ford.txt

# Number of files fund
wc -l exclude_files_ford.txt

# List of modified files in the last commit
modif=$(git diff --name-only --diff-filter=ACMRT f7eeea5954823195437b68c4db0fd9f113aa7774 25d6bb3b45cbb47ff6a03d86125b44858fb0695a)

echo $modif

# Remove each modified file from the list to exclude
for path in $modif
do
    file=$(basename $path)
    /opt/miniconda3/bin/sed -i "/$file/d" exclude_files_ford.txt
done

wc -l exclude_files_ford.txt

# Print the list to exclude in FORD
/opt/miniconda3/bin/sed "s/^/exclude: /" exclude_files_ford.txt > toto.txt

# concatenate the various parts for the FORD config file
cat << EOF_beg > beg_config
---
project: CABLE
author: CABLE community
src_dir: ../src
exclude_dir: ../src/coupled
output_dir: site/api
coloured_edges: true
preprocess: false
display: public
display: private
EOF_beg

cat << EOF_end > end_config
---

The Community Atmosphere Biosphere Land Exchange (CABLE) model is developed by a community of users under the CSIRO - MIT BSD license.

In particular, CABLE is used as the land surface model for the ACCESS coupled climate model supported in Australia for the Australian research community.
EOF_end

cat beg_config toto.txt end_config > documentation/cable_FORD.md

