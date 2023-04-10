###############################################################################
#
# Generate global landmask file for CABLE from a list of lat, lon tupels in
# a file. This tool takes care of CABLE specifications, i.e. the flipped
# Latitude axis. If you want a standard-landmask use the cable switch
#
# LN 09/2015
#
# VH 08/2018: slight modification so that script works for grid that doesn't cover whole globe.
################################################################################

import sys
sys.path.insert(1,'../scripts'); from lnutils import latlon2ixjy
import numpy as np
import netCDF4 as nc
import pandas as pd

#===============================================================================
# User modification (enter your settings below this line)
#===============================================================================

# path to output/input if desired
path="/home/599/jk8585/CABLE_run/parallel_runs/TestPoint/mask"

# name of maskfile to be created
maskfname="/home/599/jk8585/CABLE_run/parallel_runs/TestPoint/mask/TestPoint_landmask.nc"
# input file with lat, lon tupels (i.e. comma separated, 1 tupel per row)
latlonfile=path+'latlonlist_3points.csv'
gridinfo_file="/g/data/x45/CABLE-AUX/offline/gridinfo_CSIRO_1x1.nc"


# spatial resolution (in case ny,nx aren't integer multiples
res = 1.0 # spatial resolution (in case ny,nx aren't integer multiples
          # of 360. 180 resectively)


# range of lats needs lons (needs to be compatible with CABLE grid-info)
latmin =  -60.
latmax =   90.
lonmin = -180.
lonmax =  180.

# CABLE switch if set to 'cable' flipping is done,
# if set "standard" standard landmask is generated (case-sensitive!)

cable_switch = 'cable'

# latlon string given on command line
if len(sys.argv) > 1:
    randomPixels = False
    coords = sys.argv[1]
    cc = coords.split(',')
    if len(cc) == 2:
        mlatmin = np.floor(float(cc[0])) + res * 0.5
        mlatmax = mlatmin
        mlonmin = np.floor(float(cc[1])) + res * 0.5
        mlonmax = mlonmin
    elif len(cc) == 4:
        mlatmin = np.floor(float(cc[0])) + res * 0.5
        mlatmax = np.floor(float(cc[1])) + res * 0.5
        mlonmin = np.floor(float(cc[2])) + res * 0.5
        mlonmax = np.floor(float(cc[3])) + res * 0.5
    else:
        print("'coords' must be of length 2 or 4!")
        exit(1)
else:
    randomPixels = True


#===============================================================================
# End user modification (don't touch anything below this line!!!)
#===============================================================================

if cable_switch != 'cable' and cable_switch != 'standard':
    print("Please use either 'cable' or standard' for cable_switch!")
    exit()


nx       = int( (lonmax - lonmin)/res ) # number of lons (fields in W-E direction
ny       = int( (latmax - latmin)/res ) # number of lats (fields in N-S direction)
c_offset = res * 0.5

# geometrical settings (we're doing centerpoints and lats are flipped for CABLE!)
if cable_switch == 'cable':
    latset = np.array([latmax-c_offset, latmin-c_offset, -res])
else:
    latset = np.array([latmin+c_offset, latmax+c_offset,  res])
lonset = np.array([lonmin+c_offset, lonmax+c_offset, res])

#create netCDF Dataset
rootgrp   = nc.Dataset(maskfname, 'w', format='NETCDF4')
latitude  = rootgrp.createDimension('latitude', ny)
longitude = rootgrp.createDimension('longitude', nx)

# define netCDF dimensions
latitude = rootgrp.createVariable('latitude', 'f4', ('latitude',))
latitude.units = 'degrees North'

longitude = rootgrp.createVariable('longitude', 'f4', ('longitude',))
longitude.units = 'degrees East'

# generate dimension values
lats = np.arange(latset[0], latset[1], latset[2])
lons = np.arange(lonset[0], lonset[1], lonset[2])

# populate dimensions
longitude[:] = lons
latitude[:]  = lats

# create variable mask
mask = rootgrp.createVariable('land','i2',('latitude','longitude',))
mask.units = '0:no land, 1:land'

# read pixels from csv file
if randomPixels==True:
    data = np.genfromtxt(latlonfile, delimiter=',')

    # split lats and lons
    if len(data.shape) >= 2:
        ilatlist = data[:,0]
        ilonlist = data[:,1]
    else:
        ilatlist = data[0]
        ilonlist = data[1]
else:
    mlats = np.arange(mlatmin, mlatmax+res, res)   # + res because np.arange does not include endpoint
    mlons = np.arange(mlonmin, mlonmax+res, res)
    latlist = np.array([(x,y) for x in mlats for y in mlons ])[:,0]
    lonlist = np.array([(x,y) for x in mlats for y in mlons ])[:,1]

    latlist = latlist.tolist()
    lonlist = lonlist.tolist()

    # note: the following steps have been taken care of in generate_latlonlist.py if randomPixels==True
    # read gridinfo file
    gridinfo    = nc.Dataset(gridinfo_file, mode='r')
    try:
        iveg = gridinfo.variables['iveg']
        landname='iveg'
    except KeyError:
        iveg = gridinfo.variables['land']
        landname='land'
    latgridinfo = gridinfo.variables['latitude'][:]
    longridinfo = gridinfo.variables['longitude'][:]

    # convert gridinfo to dataframe (long format)
    gridinfo_landmask = pd.DataFrame(iveg[:,:], index=latgridinfo, columns=longridinfo)
    alllonlat = pd.melt(gridinfo_landmask.reset_index(), id_vars='index', value_name=landname)
    alllonlat.columns = ["latitude", "longitude", landname]

    # select points (lats and lons) where land exists according to gridinfo file
    if landname == 'iveg':
        landlonlat = alllonlat.loc[alllonlat['latitude'].isin(latlist) & alllonlat['longitude'].isin(lonlist) & \
                                   alllonlat['iveg'].notna()]
    elif landname == 'land':
        landlonlat = alllonlat.loc[alllonlat['latitude'].isin(latlist) & alllonlat['longitude'].isin(lonlist) & \
                                   alllonlat['land'] > 0]
    else:
        print("invalid 'land' variable in land mask file!")
    ilatlist = landlonlat['latitude'].to_numpy(dtype=float, copy=True)
    ilonlist = landlonlat['longitude'].to_numpy(dtype=float, copy=True)

# get mask
tmask = latlon2ixjy(ilatlist,ilonlist,latmin,latmax,lonmin,lonmax,
                    nx,ny,mtype='mask').astype(int)

# flip mask and copy to netCDF variable
if cable_switch == 'cable':
    mask[:,:] = np.flipud(tmask)
else:
    mask[:,:] = tmask

# close netCDF-file
rootgrp.close()
print("file "+maskfname+" created")
