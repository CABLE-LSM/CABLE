## script to generate a random number of lonlat points for use in CABLE runs
import pandas as pd
import netCDF4
import sys
import random

# settings
nrpoints = int(sys.argv[1])        # number of points
basepath = '/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_run/parallel_runs/'
# path/name of gridinfo file
gridinfo_file = basepath + 'surface_data/gridinfo_CSIRO_1x1.nc'
outname = basepath + 'landmasks/latlonlist_' + str(nrpoints) + 'points.csv'

# read gridinfo file
gridinfo = netCDF4.Dataset(gridinfo_file, mode='r')
iveg     = gridinfo.variables['iveg']
latmask  = gridinfo.variables['latitude'][:]
lonmask  = gridinfo.variables['longitude'][:]

# convert gridinfo to dataframe (long format)
landmask = pd.DataFrame(iveg[:, :], index=latmask, columns=lonmask)
alllonlat = pd.melt(landmask.reset_index(), id_vars='index', value_name='iveg')
alllonlat.columns = ["latitude", "longitude", "iveg"]

# randomly select points (lats and lons) where land exists
landlonlat = alllonlat.loc[(alllonlat['iveg'] > -1) &
                           (alllonlat['iveg'] < 15)]
points = landlonlat.iloc[random.sample(range(0, landlonlat.shape[0] - 1),
                                       nrpoints), :]
lonlat = points[['latitude', 'longitude']]

# write to csv file
lonlat.to_csv(outname, sep=",", header=False, index=False)
