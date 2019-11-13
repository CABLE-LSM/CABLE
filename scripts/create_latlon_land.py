#/usr/bin/env python
#
# Script to generate a random number of lonlat points for use in CABLE runs
# Same as script generate_latlonlist.py but using variable land (0/1)
# instead of variable iveg (1-17).
#
# Directory and land mask file is hard coded
#
import pandas as pd
import netCDF4
import sys
import random

## settings - Matthias@Explor
nrpoints = int(sys.argv[1])        # number of points
basepath = '/home/oqx29/zzy20/data/crujra/daily_1deg'
gridinfo_file = basepath + '/glob_ipsl_1x1.nc'   # path/name of gridinfo file
outname = './latlonlist_' + str(nrpoints) + 'points.csv'

# read gridinfo file
gridinfo = netCDF4.Dataset(gridinfo_file,mode='r')
land     = gridinfo.variables['land']
latmask  = gridinfo.variables['latitude'][:]
lonmask  = gridinfo.variables['longitude'][:]

# convert gridinfo to dataframe (long format)
landmask=pd.DataFrame(land[:,:],index=latmask,columns=lonmask)
alllatlon=pd.melt(landmask.reset_index(),id_vars='index',value_name='land')
alllatlon.columns = ["latitude","longitude","land"]

# randomly select points (lats and lons) where land exists
landlatlon=alllatlon.loc[alllatlon['land'] > 0]
points=landlatlon.iloc[random.sample(range(0,landlatlon.shape[0]-1),nrpoints),:]
latlon=points[['latitude','longitude']]

# write to csv file
latlon.to_csv(outname,sep=",",header=False,index=False)
