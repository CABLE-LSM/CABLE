gridinfo_CSIRO_1x1.nc: from /g/data/x45/CABLE-AUX/offline/gridinfo_CSIRO_1x1.nc
trendy.grid.1deg: same as v10 and v11

/landmasks
/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc
/g/data/x45/TRENDY_test_suite/mask_1000pts/landmask_1000pts.nc
/g/data/x45/TRENDY_test_suite/mask_1000pts/latlonlist_1000points.csv


added May 2024:
input.grid.1deg, which is equal to input.grid.1deg but excludes latitudes below 60degS.
It is used primarily for input preprocessing since CABLE in the TRENDY setup does not 
simulate this region.