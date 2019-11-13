#!/usr/bin/env python
from __future__ import print_function
"""


Example
-------
  python extract_latlon.py -l latlon -o omask.nc imask


History
-------
Written  Matthias Cuntz Nov 2019
"""

# -------------------------------------------------------------------------
# Command line arguments
#

def filebase(f):
    f1 = f
    if f.startswith('..'):
        f1 = f[2:]
    elif f.startswith('.'):
        f1 = f[1:]
    else:
        f1 = f
    if '.' in f1:
        return f[0:f.rfind(".")]
    else:
        return f

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

latlon  = None
ofile   = None
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Extracts latitude/longitude point(s) from input file.'''))
parser.add_argument('-l', '--latlon', action='store',
                    default=latlon, dest='latlon', metavar='latlon_string_or_file',
                        help='comma-separated string of one lat,lon or input file with several lat,lon entries.')
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-extract_latlon.nc).')
parser.add_argument('-z', '--zip', action='store_true', default=izip, dest='izip',
                    help='Use netCDF4 variable compression (default: same format as input file).')
parser.add_argument('ifile', nargs='?', default=None, metavar='input_netcdf',
                   help='input netcdf file.')
args   = parser.parse_args()
latlon = args.latlon
ofile  = args.ofile
izip   = args.izip
ifile  = args.ifile
del parser, args

if ifile is None:
    raise IOError('Input file must be given.')
if latlon is None:
    raise IOError('-l latlon must be either string of the form "lat,lon" of file having one or more lines of "lat,lon" combinations.')

import numpy as np
import netCDF4 as nc
import os
import time as ptime
tstart = ptime.time()

# -------------------------------------------------------------------------
# Determine output dimension
#

if os.path.exists(latlon):
    dat = np.genfromtxt(latlon,delimiter=',')
    if len(dat) >= 2:
        lats = dat[:,0]
        lons = dat[:,1]
    else:
        lats = dat[0]
        lons = dat[1]
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
    sifile[-2] = sifile[-2]+'-extract_latlon'
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
if 'lat' in fi.dimensions:
    latname = 'lat'
else:
    latname = 'latitude'
if 'lin' in fi.dimensions:
    lonname = 'lon'
else:
    lonname = 'longitude'
ilats = fi.variables[latname][:]
ilons = fi.variables[lonname][:]
olats = []
olons = []
for ii in range(nlatlon):
    olats.append(closest(ilats, lats[ii])
    olons.append(closest(ilons, lons[ii])
    
# Copy dimensions
for d in fi.dimensions.values():
    l = None if d.isunlimited() else len(d)
    if d.name=='lat' or d.name=='latitude' or d.name=='lon' or d.name=='longitude': l = nlatlon
    fo.createDimension(d.name, l)
    
# Create output variables
# patchfrac = fi.variables['patchfrac']
ntime = fi.dimensions['time'].size
# create static variables (non-time dependent)
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        invardef = getVariableDefinition(ivar)
        if izip: invardef.update({'zlib':True})
        ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
        for k in ivar.ncattrs():
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
        !!!
        ovar = fo.variables[ivar.name]
        ovar[:] = ivar[:]
# copy dynamic variables
# verbose: t10 = ntime // 10
for tt in range(ntime):
    # verbose: if tt%t10 == t10-1: print('{:d}%'.format((tt+1)//t10 * 10))
    for ivar in fi.variables.values():
        if (ivar.name in variables) and ('time' in ivar.dimensions):
            ovar = fo.variables[ivar.name]
            invar = ivar[tt,...]
            iiv = variables.index(ivar.name)
            ovar[tt,...] = invar * convfac[iiv]

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

tstop = ptime.time()
print('Finished in [s]: {:.2f}'.format(tstop-tstart))
