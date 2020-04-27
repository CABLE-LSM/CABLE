#!/usr/bin/env python
from __future__ import print_function
"""
usage: extract_latlon.py [-h] [-o output_netcdf] [-z] lat,lon input_netcdf

Extracts latitude/longitude point from input file.

positional arguments:
  lat,lon input_netcdf  lat,lon string & input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-extract_latlon.nc).
  -z, --zip             Use netCDF4 variable compression (default: same format as input file).


Example
-------
  python extract_latlon.py -o ofile.nc 48.5,8.1 ifile.nc

same as
  ncks -O -d latitude,48.5 -d longitude,8.1 ifile.nc ofile.nc
or
  cdo -sellonlatbox,7.75,8.45,48.1,48.9 ifile.nc ofile.nc


History
-------
Written  Matthias Cuntz, Nov 2019
Modified Matthias Cuntz, Jan 2020 - copy global attributes and append history
         Matthias Cuntz, Apr 2020 - use create_variable function as in sum_patchfrac.py
"""

# -------------------------------------------------------------------------
# Functions
#

def _get_variable_definition(ncvar):
    out  = ncvar.filters() if ncvar.filters() else {}
    dims = list(ncvar.dimensions)
    if 'x' in dims: dims[dims.index('x')] = 'lon'
    if 'y' in dims: dims[dims.index('y')] = 'lat'
    # Do not set chunksize because it should be only nlatlon
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
        # "fill_value" : ncvar._FillValue if "_FillValue" in dir(ncvar) else None,
        "fill_value" : ifill,
    })
    return out


def _create_output_variables(fi, fo, time=False, izip=False):
    # loop over input variables
    for ivar in fi.variables.values():
        if time:
            itime = 'time' in ivar.dimensions
        else:
            itime = 'time' not in ivar.dimensions
        if itime:
            # create output variable
            invardef = _get_variable_definition(ivar)
            if izip: invardef.update({'zlib':True})
            ovar = fo.createVariable(invardef.pop("name"), invardef.pop("dtype"), **invardef)
            for k in ivar.ncattrs():
                iattr = ivar.getncattr(k)
                if (k != 'missing_value') and (k != '_FillValue'):
                    ovar.setncattr(k, iattr)
    return


# get index in vector vec of value closest to num
import numpy as np
def _closest(vec, num):
    return np.argmin(np.abs(np.array(vec)-num))


# -------------------------------------------------------------------------
# Command line
#

import argparse

ofile   = None
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Extracts latitude/longitude point from input file.'''))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-extract_latlon.nc).')
parser.add_argument('-z', '--zip', action='store_true', default=izip, dest='izip',
                    help='Use netCDF4 variable compression (default: same format as input file).')
parser.add_argument('ifiles', nargs='*', default=None, metavar='lat,lon input_netcdf',
                   help='lat,lon string & input netcdf file.')
args   = parser.parse_args()
ofile  = args.ofile
izip   = args.izip
ifiles = args.ifiles
del parser, args

if len(ifiles) < 2:
    raise IOError('lat,lon string and input file must be given.')
latlon = ifiles[0]
ifile  = ifiles[1]

import numpy as np
import netCDF4 as nc
import os
import sys
import time as ptime
tstart = ptime.time()

# -------------------------------------------------------------------------
# Get lat/lon indices
#

lats, lons = latlon.split(',')
lats    = float(lats)
lons    = float(lons)
nlatlon = 1

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
iilats = _closest(ilats, lats)
iilons = _closest(ilons, lons)

# Copy global attributes
for k in fi.ncattrs():
    iattr = fi.getncattr(k)
    if (k == 'history'):
        iattr = iattr + '\n' + ptime.asctime() + ': ' + ' '.join(sys.argv)
    fo.setncattr(k, iattr)
    
# Copy dimensions
for d in fi.dimensions.values():
    l = None if d.isunlimited() else len(d)
    if d.name=='lat' or d.name=='latitude' or d.name=='lon' or d.name=='longitude': l = nlatlon
    fo.createDimension(d.name, l)

#
# Create output variables
#

# create static variables (independent of time)
_create_output_variables(fi, fo, time=False, izip=izip)

# create dynamic variables (time dependent)
_create_output_variables(fi, fo, time=True, izip=izip)

#
# Copy variables from in to out
#

# copy static variables
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        ovar = fo.variables[ivar.name]
        if (latname in ivar.dimensions) and (lonname in ivar.dimensions):
            ovar[:] = ivar[iilats,iilons,...]
        elif (latname in ivar.dimensions):
            ovar[:] = ivar[iilats,...]
        elif (lonname in ivar.dimensions):
            ovar[:] = ivar[iilons,...]
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
                    ovar[tt,...] = ivar[tt,iilats,iilons,...]
                elif (latname in ivar.dimensions):
                    ovar[t,...] = ivar[tt,iilats,...]
                elif (lonname in ivar.dimensions):
                    ovar[tt,...] = ivar[tt,iilons,...]
                else:
                    ovar[tt,...] = ivar[tt,...]

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

tstop = ptime.time()
print('Finished in [s]: {:.2f}'.format(tstop-tstart))
