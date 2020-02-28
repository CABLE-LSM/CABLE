################################################################################
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
import math
from lnutils import latlon2ixjy
import numpy as np
from netCDF4 import Dataset

#===============================================================================
# User modification (enter your settings below this line) 
#===============================================================================

# path to output/input if desired
path='/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_run/parallel_runs/landmasks/'

# name of maskfile to be created
maskfname=path+'test3points_1x1.nc'
# input file with lat, lon tupels (i.e. comma separated, 1 tupel per row)
latlonfile=path+'latlonlist_3points.csv'

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

try:
    sys.argv[1]
except IndexError:
    randomPixels = True
else:
    coords = sys.argv[1]
    randomPixels = False
    mlatmin = math.floor(float(coords.split(',')[0])) + res * 0.5
    mlatmax = math.floor(float(coords.split(',')[1])) + res * 0.5
    mlonmin = math.floor(float(coords.split(',')[2])) + res * 0.5
    mlonmax = math.floor(float(coords.split(',')[3])) + res * 0.5


#===============================================================================
# End user modification (don't touch anything below this line!!!)
#===============================================================================

if cable_switch != 'cable' and cable_switch != 'standard':
    print("Please use either 'cable' or standard' for cable_switch!")
    exit()


nx       = np.int( (lonmax - lonmin)/res ) # number of lons (fields in W-E direction
ny       = np.int( (latmax - latmin)/res ) # number of lats (fields in N-S direction)
c_offset = res * 0.5

# geometrical settings (we're doing centerpoints and lats are flipped for CABLE!) 
if cable_switch == 'cable':
    latset   = np.array([latmax-c_offset, latmin-c_offset, -1.*res])
else:
    latset   = np.array([latmin+c_offset, latmax+c_offset,     res])
lonset   = np.array([lonmin+c_offset, lonmax+c_offset ,    res])

#create netCDF Dataset   
rootgrp   = Dataset(maskfname, 'w', format='NETCDF4')
latitude  = rootgrp.createDimension('latitude',ny)
longitude = rootgrp.createDimension('longitude',nx)

# define netCDF dimensions
latitude = rootgrp.createVariable('latitude','f4',('latitude',))
latitude.units = 'degrees North'

longitude = rootgrp.createVariable('longitude','f4',('longitude',))
longitude.units = 'degrees East'

# generate dimension values
lats = np.arange(latset[0],latset[1],latset[2])
lons = np.arange(lonset[0],lonset[1],lonset[2])

# populate dimensions
longitude[:] = lons
latitude[:] = lats

# create variable mask
mask = rootgrp.createVariable('land','i2',('latitude','longitude',))
mask.units = '0:no land, 1:land'


# read pixels from csv file
if randomPixels==True:
    data = np.genfromtxt(latlonfile,delimiter=',')

    # split lats and lons
    if len(data.shape) >= 2:
        latlist=data[:,0]
        lonlist=data[:,1]
    else:
        latlist=data[0]
        lonlist=data[1]
else:
    mlats=np.arange(mlatmin,mlatmax+res,res)   # + res because np.arange does not include endpoint
    mlons=np.arange(mlonmin,mlonmax+res,res)
    latlist=np.array([(x,y) for x in mlats for y in mlons ])[:,0]
    lonlist=np.array([(x,y) for x in mlats for y in mlons ])[:,1]
    
# get mask 
tmask = latlon2ixjy(latlist,lonlist,latmin,latmax,lonmin,lonmax,
                    nx,ny,mtype='mask').astype(int)

# flip mask and copy to netCDF variable
if cable_switch == 'cable':
    mask[:,:] = np.flipud(tmask)
else:
    mask[:,:] = tmask

# close netCDF-file
rootgrp.close()
print("file "+maskfname+" created")
