#!/bin/bash

# Set configuration
WORK_DIR="~/work/library_shuffle"

CABLE_REPO="git@github.com:CABLE-LSM/CABLE.git"
CABLE_BRANCH="538-the-great-library-shuffle"
CABLE_DIR="{WORK_DIR}/cable"
CABLE_SHARED_DIR=$CABLE_DIR/src/shared
CABLE_ESM15_DIR=$CABLE_DIR/src/coupled/ESM1.5
CABLE_ESM16_DIR=$CABLE_DIR/src/coupled/esm16
CABLE_ESM_DIR=$CABLE_DIR/src/coupled/esm

UM7_REPO="git@github.com:ACCESS-NRI/UM7.git"
UM7_BRANCH="64-the-great-library-shuffle"
UM7_DIR="{WORK_DIR}/um7"
UM7_ESM15_DIR=$UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/ESM1.5
UM7_CABLE_DIR=$UM7_DIR/umbase_hg3/src/atmosphere/land_surface/cable

# Delete and remake work dir
rm -rf $WORK_DIR
mkdir -p $WORK_DIR

# Clone everything
cd $WORK_DIR
git clone $CABLE_REPO $CABLE_DIR
git clone $UM7_REPO $UM7_DIR

# Check out the CABLE branch for the shuffle
cd $CABLE_DIR
git fetch origin
git checkout $CABLE_BRANCH

# Check out specific branches
cd $UM7_DIR
git fetch origin
git checkout $UM7_BRANCH

# Create a src/shared dir in CABLE
mkdir -p $CABLE_SHARED_DIR

# Move from offline to shared
mv $CABLE_DIR/src/offline/cable_LUC_EXPT.F90 $CABLE_SHARED_DIR/
mv $CABLE_DIR/src/offline/cable_phenology.F90 $CABLE_SHARED_DIR/
mv $CABLE_DIR/src/offline/casa_ncdf.F90 $CABLE_SHARED_DIR/
mv $CABLE_DIR/src/offline/casa_offline_inout.F90 $CABLE_SHARED_DIR/

# Remove them from offline (done in move above)

# Remove from ESM1.5 (none)

# Remove from esm16
rm -rf $CABLE_ESM16_DIR/cable_LUC_EXPT.F90
rm -rf $CABLE_ESM16_DIR/cable_phenology.F90
rm -rf $CABLE_ESM16_DIR/casa_ncdf.F90
rm -rf $CABLE_ESM16_DIR/casa_offline_inout.F90

# Add shared to CMakeLists.txt (manual)

# Create the land_surface/cable dir
mkdir -p $UM7_CABLE_DIR

# Move source into it (not sure if we need to delete the src location?)
cp $UM7_ESM15_DIR/casa_types.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/pack_mod_cbl.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_rad_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_um_init.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_um_init_subrs.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_um_tech.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/allocate_soil_params_cbl.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/allocate_veg_params_cbl.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_cbm.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_explicit_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_hyd_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cable_implicit_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cbl_um_init_soil.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cbl_um_init_soilsnow.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_DIR/cbl_um_update_soilsnow.F90 $UM7_CABLE_DIR/

# casa_um_inout.F90 (currently in esm16/ but it needs to be deleted from there)
# it is not there?

# Put in $CABLE_DIR/src/coupled/esm/, copied from $UM7_DIR/src/atmosphere/CABLE/src/coupled/ESM1.5/
cp $UM7_ESM15_DIR/cable_pft_params.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/cable_soil_params.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/cable_surface_tyoes.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/cable_landuse.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/cable_define_types.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/cable_iovars.F90 $CABLE_ESM_DIR/

# Perform the git commit/push manually once we verify this is all moved correctly.