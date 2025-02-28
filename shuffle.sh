#!/bin/bash

set -e

# Set configuration
WORK_DIR="${HOME}/work/library_shuffle"

CABLE_REPO="git@github.com:CABLE-LSM/CABLE.git"
CABLE_BRANCH="538-the-great-library-shuffle"
CABLE_DIR="${WORK_DIR}/cable"
CABLE_SHARED_DIR=$CABLE_DIR/src/shared
CABLE_ESM15_DIR=$CABLE_DIR/src/coupled/ESM1.5
CABLE_ESM16_DIR=$CABLE_DIR/src/coupled/esm16
CABLE_ESM_DIR=$CABLE_DIR/src/coupled/esm
CABLE_COUPLED_SHARED_DIR=$CABLE_DIR/src/coupled/shared

UM7_REPO="git@github.com:ACCESS-NRI/UM7.git"
UM7_BRANCH="64-the-great-library-shuffle"
UM7_DIR="${WORK_DIR}/um7"
UM7_ESM15_DIR=$UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/ESM1.5
UM7_COUPLED_SHARED_DIR=$UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/shared
UM7_ESM15_SUBDIR=$UM7_ESM15_DIR/CABLEfilesFromESM1.5
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

echo "check content of UM7_ESM15_DIR"
ls $UM7_ESM15_DIR/

# Create a src/shared dir in CABLE
mkdir -p $CABLE_SHARED_DIR

# Move from offline to shared
mv -v $CABLE_DIR/src/offline/cable_LUC_EXPT.F90 $CABLE_SHARED_DIR/
mv -v $CABLE_DIR/src/offline/cable_phenology.F90 $CABLE_SHARED_DIR/
mv -v $CABLE_DIR/src/offline/casa_ncdf.F90 $CABLE_SHARED_DIR/
mv -v $CABLE_DIR/src/offline/casa_offline_inout.F90 $CABLE_SHARED_DIR/

echo rm $CABLE_DIR/src/offline/cable_pft_params.F90
echo rm $CABLE_DIR/src/offline/cable_soil_params.F90

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
cp $UM7_ESM15_SUBDIR/cable_rad_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_um_init.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_um_init_subrs.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_um_tech.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/allocate_soil_params_cbl.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/allocate_veg_params_cbl.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_cbm.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_explicit_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_hyd_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_implicit_driver.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cbl_um_init_soil.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cbl_um_init_soilsnow.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cbl_um_update_soilsnow.F90 $UM7_CABLE_DIR/
cp $UM7_ESM15_SUBDIR/cable_read_nml.F90 $UM7_CABLE_DIR/

# casa_um_inout.F90 (currently in esm16/ but it needs to be deleted from there)
# it is not there?

# Put in $CABLE_DIR/src/coupled/esm/, copied from $UM7_DIR/src/atmosphere/CABLE/src/coupled/ESM1.5/
mkdir -p $CABLE_ESM_DIR
cp $UM7_ESM15_DIR/cable_surface_types.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/casa_landuse.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_SUBDIR/cable_define_types.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_SUBDIR/cable_iovars.F90 $CABLE_ESM_DIR/
cp $UM7_ESM15_DIR/cable_pft_params_mod.F90 $CABLE_SHARED_DIR/
cp $UM7_ESM15_DIR/cable_soil_params_mod.F90 $CABLE_SHARED_DIR/

# Remove unnecessary directories
rm -r $CABLE_ESM16_DIR
rm -r $CABLE_ESM15_DIR
rm -r $CABLE_DIR/src/coupled/ACCESS-CM2

# Replace CABLE directory in UM7 with the one with the reorganised files. We can't do a clean sweep as the sync between the repos isn't complete.
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/science $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/util $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/params $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/
mv -v $UM7_COUPLED_SHARED_DIR $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/
cp -r $CABLE_DIR/src $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/offline
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/AM3
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/JAC
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/shared
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/params
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/science
rm -r $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/util
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/params $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/science $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/util $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/ 
mv -v $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/shared $UM7_DIR/umbase_hg3/src/atmosphere/CABLE/src/coupled/
# Uncomment the #define lines for the pre-processor flags in the science code manually.

# Perform the git commit/push manually once we verify this is all moved correctly.