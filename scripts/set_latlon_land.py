#!/usr/bin/env python
from __future__ import print_function
"""
usage: set_latlon_land.py [-h] [-o output_netcdf] [-z] latlon input_netcdf

Sets the land variable to 1 at given latitude/longitude points, otherwise 0.

positional arguments:
  latlon input_netcdf   lat,lon string or latlon input file & input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-masked.nc).
  -z, --zip             Use netCDF4 variable compression (default: same format as input file).


Example
-------
  python set_latlon_land.py -o ofile.nc 48.234,9.567 ifile.nc
  python set_latlon_land.py -o ofile.nc latlonFile ifile.nc


History
-------
Written  Matthias Cuntz Nov 2019
"""

# -------------------------------------------------------------------------
# Functions
#

# Copied from jams.netcdf4
def getVariableDefinition(ncvar):
    out  = ncvar.filters() if ncvar.filters() else {}
    dims = list(ncvar.dimensions)
    if 'x' in dims: dims[dims.index('x')] = 'lon'
    if 'y' in dims: dims[dims.index('y')] = 'lat'
    chunks = ncvar.chunking() if not isinstance(ncvar.chunking(), str) else None
    out.update({
        "name"       : ncvar.name,
        "dtype"      : ncvar.dtype,
        "dimensions" : dims,
        "chunksizes" : chunks,
        "fill_value" : ncvar._FillValue if "_FillValue" in dir(ncvar) else None,
        })
    return out

# get index in vector vec of value closest to num
def closest(vec, num):
    return np.argmin(np.abs(vec-num))


# -------------------------------------------------------------------------
# Command line
#

import argparse

ofile   = None
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Sets the land variable to 1 at given latitude/longitude points, otherwise 0.'''))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-masked.nc).')
parser.add_argument('-z', '--zip', action='store_true', default=izip, dest='izip',
                    help='Use netCDF4 variable compression (default: same format as input file).')
parser.add_argument('ifiles', nargs='*', default=None, metavar='latlon input_netcdf',
                   help='lat,lon string or latlon input file & input netcdf file.')
args   = parser.parse_args()
ofile  = args.ofile
izip   = args.izip
ifiles = args.ifiles
del parser, args

if len(ifiles) < 2:
    raise IOError('lat,lon string or latlon file, and input file must be given.')
latlon = ifiles[0]
ifile  = ifiles[1]

import numpy as np
import netCDF4 as nc
import os
import time as ptime
tstart = ptime.time()

# -------------------------------------------------------------------------
# Get lat/lon indices
#

if os.path.exists(latlon):
    dat = np.genfromtxt(latlon, delimiter=',')
    if len(dat) > 1:
        lats = dat[:,0]
        lons = dat[:,1]
    else:
        lats = np.array([dat[0]])
        lons = np.array([dat[1]])
else:
    lats, lons = latlon.split(',')
    lats = np.array([float(lats)])
    lons = np.array([float(lons)])
nlatlon = lats.size

# -------------------------------------------------------------------------
# Copy data
#

if ofile is None: # Default output filename
    sifile = ifile.split('.')
    sifile[-2] = sifile[-2]+'-masked'
    ofile = '.'.join(sifile)
fi = nc.Dataset(ifile, 'r')
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    # try:
    fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    # except:
    #     fo = nc.Dataset(ofile, 'w', format='NETCDF4')

# lat/lon indices
if 'latitude' in fi.dimensions.keys():
    latname = 'latitude'
else:
    latname = 'lat'
if 'longitude' in fi.dimensions:
    lonname = 'longitude'
else:
    lonname = 'lon'
ilats = fi.variables[latname][:]
ilons = fi.variables[lonname][:]
iilats = []
iilons = []
for ii in range(nlatlon):
    iilats.append(closest(ilats, lats[ii]))
    iilons.append(closest(ilons, lons[ii]))
    
# Copy dimensions
for d in fi.dimensions.values():
    fo.createDimension(d.name, None if d.isunlimited() else len(d))
    
# Create output variables
# patchfrac = fi.variables['patchfrac']
# create static variables (non-time dependent)
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        invardef = getVariableDefinition(ivar)
        if izip: invardef.update({'zlib':True})
        ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
        for k in ivar.ncattrs():
            if k != '_FillValue':
                iattr = ivar.getncattr(k)
                ovar.setncattr(k, iattr)
# create dynamic variables (time dependent)
for ivar in fi.variables.values():
    if 'time' in ivar.dimensions:
        invardef = getVariableDefinition(ivar)
        if izip: invardef.update({'zlib':True})
        ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
        for k in ivar.ncattrs():
            iattr = ivar.getncattr(k)
            ovar.setncattr(k, iattr)
# Copy variables from in to out
# copy static variables
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        ovar = fo.variables[ivar.name]
        if ivar.name == 'land':
            ovar[:] = 0
            ovar[iilats,iilons] = 1
        else:
            ovar[:] = ivar[:]
# copy dynamic variables
if 'time' in fi.dimensions:
    ntime = fi.dimensions['time'].size
    # verbose: t10 = ntime // 10
    for tt in range(ntime):
        # verbose: if tt%t10 == t10-1: print('{:d}%'.format((tt+1)//t10 * 10))
        for ivar in fi.variables.values():
            if 'time' in ivar.dimensions:
                ovar = fo.variables[ivar.name]
                ovar[tt,...] = ivar[tt,...]

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

tstop = ptime.time()
print('Finished in [s]: {:.2f}'.format(tstop-tstart))
