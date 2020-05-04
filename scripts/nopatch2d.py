#!/usr/bin/env python
from __future__ import print_function
"""
usage: nopatch2d.py [-h] [-o output_netcdf] [-v] [-z] [input_netcdf]

Transform Cable or Casa output with summed patches to 2D arrays with latitude/longitude.

positional arguments:
  input_netcdf          input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-2d.nc).
  -v, --verbose         Feedback during copy (default: no feedback).
  -z, --zip             Use netCDF4 variable compression (default: same format
                        as input file).


Example
-------
  python nopatch2d.py -z -o cru_out_casa_2009_2011-no_patch-2d.nc cru_out_casa_2009_2011-no_patch.nc


History
-------
Written  Matthias Cuntz, May 2020

ToDo: create 2d fields per time step rather than whole arrays at once.
"""

# -------------------------------------------------------------------------
# Functions
#

# https://github.com/Unidata/python-workshop/blob/fall-2016/notebooks/netcdf-by-coordinates.ipynb
import numpy as np
from scipy.spatial import cKDTree
class llKDTree(object):
    def __init__(self, latvar, lonvar):
        rad_factor = np.pi/180.0 # for trigonometry, need angles in radians
        self.latvals = latvar * rad_factor
        self.lonvals = lonvar * rad_factor
        self.shape   = self.latvals.shape
        clat, clon = np.cos(self.latvals), np.cos(self.lonvals)
        slat, slon = np.sin(self.latvals), np.sin(self.lonvals)
        clat_clon  = clat*clon
        clat_slon  = clat*slon
        triples  = list(zip(np.ravel(clat*clon), np.ravel(clat*slon), np.ravel(slat)))
        self.kdt = cKDTree(triples) # build tree during initialisation

    def query(self, lat0, lon0):
        rad_factor = np.pi/180.0 
        lat0_rad = lat0 * rad_factor
        lon0_rad = lon0 * rad_factor
        clat0, clon0 = np.cos(lat0_rad), np.cos(lon0_rad)
        slat0, slon0 = np.sin(lat0_rad), np.sin(lon0_rad)
        dist_sq_min, minindex_1d = self.kdt.query([clat0*clon0, clat0*slon0, slat0]) # flattened index
        iy_min, ix_min = np.unravel_index(minindex_1d, self.shape) # 2d index
        return iy_min, ix_min


# -------------------------------------------------------------------------
# Command line
#

import argparse

ofile   = None
verbose = False
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Transform Cable or Casa output with summed patches to 2D arrays with latitude/longitude.'''))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-2d.nc).')
parser.add_argument('-v', '--verbose', action='store_true', default=verbose, dest='verbose',
                    help='Feedback during copy (default: no feedback).')
parser.add_argument('-z', '--zip', action='store_true', default=izip, dest='izip',
                    help='Use netCDF4 variable compression (default: same format as input file).')
parser.add_argument('ifile', nargs='?', default=None, metavar='input_netcdf',
                   help='input netcdf file.')
args    = parser.parse_args()
ofile   = args.ofile
verbose = args.verbose
izip    = args.izip
ifile   = args.ifile
del parser, args

if ifile is None:
    raise IOError('Input file must be given.')

import sys
# import numpy as np
import netCDF4 as nc
import cablepop as cp
import time as ptime
if verbose:
    tstart = ptime.time()

# -------------------------------------------------------------------------
# Copy data
#

#
# Open files
#
if ofile is None: # Default output filename
    ofile = cp.set_output_filename(ifile, '-2d')
fi = nc.Dataset(ifile, 'r')
if verbose: print('Input file:  ', ifile)
if verbose: print('Output file: ', ofile)
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    if 'file_format' in dir(fi):
        fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    else:
        fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')

# get latitude/longitude indices
if 'local_lat' in fi.variables:
    ilats = fi.variables['local_lat'][:]
    ilons = fi.variables['local_lon'][:]
else:
    ilats = fi.variables['latitude'][:]
    ilons = fi.variables['longitude'][:]
nland = ilats.size
if 'x' in fi.dimensions: # use existing grid
    olat = fi.variables['y'][:]
    olon = fi.variables['x'][:]
    nlat   = olat.size
    nlon   = olon.size
else:                    # great new global grid -60 to +90 latitudes
    dlat = np.abs(np.diff(np.unique(np.sort(ilats)))).min() # 0.5, 1 degree
    dlon = np.abs(np.diff(np.unique(np.sort(ilons)))).min() # 0.5, 1 degree
    nlat = np.rint(150./dlat).astype(np.int) # 300, 150
    nlon = np.rint(360./dlon).astype(np.int) # 720, 360
    clat = ilats.min() % 1. # 0.0 or 0.25, 0.0 or 0.5
    clon = ilons.min() % 1. # 0.0 or 0.25, 0.0 or 0.5
    olat = -60.  + clat + np.arange(nlat)/float(nlat-1) * (150.-dlat) # new lats
    olon = -180. + clon + np.arange(nlon)/float(nlon-1) * (360.-dlon) # new lons
olon2d, olat2d = np.meshgrid(olon, olat) # new lats, lons in 2D
lltree = llKDTree(olat2d, olon2d) # KD-tree
iidl   = np.arange(nland, dtype=np.int) # indices of land in input grid
oidx   = np.empty(nland, dtype=np.int)  # indices of lon in output grid
oidy   = np.empty(nland, dtype=np.int)  # indices of lat in output grid
for i in range(nland):
    iy, ix = lltree.query(ilats[i], ilons[i])
    oidx[i] = ix
    oidy[i] = iy

# Copy global attributes, adding script
cp.set_global_attributes(fi, fo, add={'history':ptime.asctime()+': '+' '.join(sys.argv)})

# Copy dimensions
cp.create_dimensions(fi, fo, removedim=['land', 'ntile'], adddim={'x':nlon, 'y':nlat})

# create static variables (independent of time)
if 'local_lat' in fi.variables:
    renvar = {} #{'latitude':'nav_lat', 'longitude':'nav_lon'}
else:
    renvar = {}
cp.create_variables(fi, fo, time=False, izip=izip, fill=-1e33, chunksizes=False,
                    renamevar=renvar, replacedim={'land':('y','x'), 'ntile':('y','x')})
# create dynamic variables (time dependent)
cp.create_variables(fi, fo, time=True, izip=izip, fill=-1e33, chunksizes=False,
                    renamevar=renvar, replacedim={'land':('y','x'), 'ntile':('y','x')})

#
# Copy variables from in to out expanding the land dimension to y, x
#

# copy static and dynamic variables
if verbose: print('Copy variables')
# np.seterr(over='ignore', under='ignore', invalid='ignore') # multiply masked elements
nvar = len(fi.variables.keys())
if nvar < 20:
    n10 = 1
else:
    n10 = nvar // 10
n = 0
for ivar in fi.variables.values():
    if verbose and (n > 0) and ((n % n10) == 0): print('  {:d}%'.format(int(n/nvar*100.)))
    n += 1
    ovar  = fo.variables[ivar.name]
    invar = ivar[:] # read whole field, otherwise times increasing sharpely
    if ('land' in ivar.dimensions) or ('ntile' in ivar.dimensions):
        outshape = list(invar.shape)
        outshape[-1] = nlat   # land to y,x
        outshape.append(nlon)
        # fill in memory, then write to disk in one go
        out = np.full(outshape, ovar._FillValue)
        if ivar.name == 'latitude':
            out = olat2d
        elif ivar.name == 'longitude':
            out = olon2d
        else:
            if len(invar.shape) == 1:
                out[oidy, oidx] = invar[iidl]
            else:
                out[..., oidy, oidx] = invar[..., iidl]
    else:
        out = invar
    ovar[:] = out

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

if verbose:
    tstop = ptime.time()
    print('Finished in [s]: {:.2f}'.format(tstop-tstart))
