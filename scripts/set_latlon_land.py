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
Written  Matthias Cuntz, Nov 2019
Modified Matthias Cuntz, Jan 2020 - copy global attributes and append history
         Matthias Cuntz, Apr 2020 - use create_variable function as in sum_patchfrac.py
                                  - lat,lon given with -l option
"""

# -------------------------------------------------------------------------
# Functions
#

# Copied from jams.netcdf4
def _get_variable_definition(ncvar):
    out  = ncvar.filters() if ncvar.filters() else {}
    dims = list(ncvar.dimensions)
    if 'x' in dims: dims[dims.index('x')] = 'lon'
    if 'y' in dims: dims[dims.index('y')] = 'lat'
    chunks = ncvar.chunking() if not isinstance(ncvar.chunking(), str) else None
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
        "chunksizes" : chunks,
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

latlon = None
ofile  = None
izip   = False
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Sets the land variable to 1 at given latitude/longitude points, otherwise 0.'''))
parser.add_argument('-l', '--latlon', action='store',
                    default=latlon, dest='latlon', metavar='lat,lon',
                        help='latitude,longitude of grid point to mask/extract, or file with several latitude,longitude coordinates.')
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-masked.nc).')
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

if latlon is None:
    raise IOError('latitude,longitude or file with latitudes,longitudes must be given with -l flag.')

if ifile is None:
    raise IOError('Netcdf input file must be given.')

import numpy as np
import netCDF4 as nc
import os
import sys
import time as ptime
# tstart = ptime.time()

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
if 'land' not in fi.variables:
    fi.close()
    print('Land variable must be in input file: '+ifile)
    import sys
    sys.exit()
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
    iilats.append(_closest(ilats, lats[ii]))
    iilons.append(_closest(ilons, lons[ii]))

# Copy global attributes
for k in fi.ncattrs():
    iattr = fi.getncattr(k)
    if (k == 'history'):
        iattr = iattr + '\n' + ptime.asctime() + ': ' + ' '.join(sys.argv)
    fo.setncattr(k, iattr)
    
# Copy dimensions
for d in fi.dimensions.values():
    fo.createDimension(d.name, None if d.isunlimited() else len(d))

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

# tstop = ptime.time()
# print('Finished in [s]: {:.2f}'.format(tstop-tstart))
