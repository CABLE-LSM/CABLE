#!/usr/bin/env python
from __future__ import print_function
"""
usage: duplicate_lon.py [-h] [-o output_netcdf] [-z] [input_netcdf]

Copies all variables of a single grid cell file into a secnd grid cell with the same lat/lon.

positional arguments:
  input_netcdf          input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-2.nc).
  -z, --zip             Use netCDF4 variable compression (default: same format
                        as input file).


Example
-------
  python duplicate_lon.py -o ofile.nc ifile.nc


History
-------
Written  Matthias Cuntz Dec 2019
"""

# -------------------------------------------------------------------------
# Functions
#

# Copied from jams.netcdf4
def getVariableDefinition(ncvar):
    out  = ncvar.filters() if ncvar.filters() else {}
    dims = list(ncvar.dimensions)
    # Do not set chunksize because it should be only 2
    # chunks = ncvar.chunking() if not isinstance(ncvar.chunking(), str) else None
    if "missing_value" in dir(ncvar):
        ifill = ncvar.missing_value
    elif "_FillValue" in dir(ncvar):
        ifill = ncvar._FillValue
    else:
        ifill = None
    out.update({
        "name"       : ncvar.name,
        "dtype"      : ncvar.dtype,
        "dimensions" : dims,
        # "chunksizes" : chunks,
        "fill_value" : ifill,
    })
    return out


# -------------------------------------------------------------------------
# Command line
#

import argparse

ofile   = None
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Copies all variables of a single grid cell file into a secnd grid cell with the same lat/lon.''')) # 1 degree to the west.'''))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-2.nc).')
parser.add_argument('-z', '--zip', action='store_true', default=izip, dest='izip',
                    help='Use netCDF4 variable compression (default: same format as input file).')
parser.add_argument('ifile', nargs='?', default=None, metavar='input_netcdf',
                   help='input netcdf file.')
args  = parser.parse_args()
ofile = args.ofile
izip  = args.izip
ifile = args.ifile
del parser, args

import numpy as np
import netCDF4 as nc
import os
import time as ptime
tstart = ptime.time()

# -------------------------------------------------------------------------
# Copy data
#

if ofile is None: # Default output filename
    sifile = ifile.split('.')
    sifile[-2] = sifile[-2]+'-2'
    ofile = '.'.join(sifile)
fi = nc.Dataset(ifile, 'r')
# Check that only one longitude present
for d in fi.dimensions.values():
    l = None if d.isunlimited() else len(d)
    if (d.name=='longitude' or d.name=='lon' or d.name=='x') and (l != 1):
        fi.close()
        raise IOError('More than one longitude present in input file: '+str(l))
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
elif 'lat' in fi.dimensions.keys():
    latname = 'lat'
else:
    latname = 'y'
if 'longitude' in fi.dimensions:
    lonname = 'longitude'
elif 'lon' in fi.dimensions:
    lonname = 'lon'
else:
    lonname = 'x'
    
# Copy dimensions
for d in fi.dimensions.values():
    l = None if d.isunlimited() else len(d)
    if (d.name=='longitude' or d.name=='lon' or d.name=='x'): l = 2
    fo.createDimension(d.name, l)
    
# Create output variables
# patchfrac = fi.variables['patchfrac']
# create static variables (non-time dependent)
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        invardef = getVariableDefinition(ivar)
        if izip: invardef.update({'zlib':True})
        ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
        for k in ivar.ncattrs():
            iattr = ivar.getncattr(k)
            if (k != 'missing_value') and (k != '_FillValue'):
                ovar.setncattr(k, iattr)
# create dynamic variables (time dependent)
for ivar in fi.variables.values():
    if 'time' in ivar.dimensions:
        invardef = getVariableDefinition(ivar)
        if izip: invardef.update({'zlib':True})
        ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
        for k in ivar.ncattrs():
            iattr = ivar.getncattr(k)            
            if (k != 'missing_value') and (k != '_FillValue'):
                ovar.setncattr(k, iattr)
# Copy variables from in to out
# copy static variables
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        ovar = fo.variables[ivar.name]
        if (latname in ivar.dimensions) and (lonname in ivar.dimensions):
            ovar[:,0,...] = ivar[:,0,...]
            ovar[:,1,...] = ivar[:,0,...]
        elif (lonname in ivar.dimensions):
            # if (ivar.name == lonname):
            #     ovar[0,...] = ivar[0,...] - 1. # one degree to the west
            # else:
            #     ovar[0,...] = ivar[0,...]
            ovar[0,...] = ivar[0,...]
            ovar[1,...] = ivar[0,...]
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
                if (latname in ivar.dimensions) and (lonname in ivar.dimensions):
                    ovar[tt,:,0,...] = ivar[tt,:,0,...]
                    ovar[tt,:,1,...] = ivar[tt,:,0,...]
                elif (lonname in ivar.dimensions):
                    ovar[tt,0,...] = ivar[tt,0,...]
                    ovar[tt,1,...] = ivar[tt,0,...]
                else:
                    ovar[tt,...] = ivar[tt,...]

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

tstop = ptime.time()
# print('Finished in [s]: {:.2f}'.format(tstop-tstart))
