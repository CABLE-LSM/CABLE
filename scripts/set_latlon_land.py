#!/usr/bin/env python3
"""
usage: set_latlon_land.py [-h] [-o output_netcdf] [-z] latlon input_netcdf

Sets the land variable to 1 at given latitude/longitude points, otherwise 0.

positional arguments:
  latlon input_netcdf   lat,lon string or latlon input file &
                        input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-masked.nc).
  -z, --zip             Use netCDF4 variable compression
                        (default: same format as input file).


Example
-------
  python set_latlon_land.py -o ofile.nc 48.234,9.567 ifile.nc
  python set_latlon_land.py -o ofile.nc latlonFile ifile.nc


History
-------
Written  Matthias Cuntz, Nov 2019
Modified Matthias Cuntz, Jan 2020
             - copy global attributes and append history
         Matthias Cuntz, Apr 2020
             - use create_variable function as in sum_patchfrac.py
             - lat,lon given with -l option
         Matthias Cuntz, Apr 2020
             - use python module cablepop
         Matthias Cuntz, May 2023
             - flake8 compatible

"""

# -------------------------------------------------------------------------
# Command line
#

import argparse

latlon = None
ofile  = None
izip   = False
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('Sets the land variable to 1 at given latitude/longitude'
                 ' points, otherwise 0.'))
parser.add_argument('-l', '--latlon', action='store',
                    default=latlon, dest='latlon', metavar='lat,lon',
                    help=('latitude,longitude of grid point to mask/extract,'
                          ' or file with several latitude,longitude'
                          ' coordinates.'))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                    help='output netcdf file name (default: input-masked.nc).')
parser.add_argument('-z', '--zip', action='store_true', default=izip,
                    dest='izip',
                    help=('Use netCDF4 variable compression (default:'
                          ' same format as input file).'))
parser.add_argument('ifile', nargs='?', default=None, metavar='input_netcdf',
                    help='input netcdf file.')
args   = parser.parse_args()
latlon = args.latlon
ofile  = args.ofile
izip   = args.izip
ifile  = args.ifile
del parser, args

if latlon is None:
    raise IOError('latitude,longitude or file with latitudes,longitudes'
                  ' must be given with -l flag.')

if ifile is None:
    raise IOError('Netcdf input file must be given.')

import numpy as np
import netCDF4 as nc
import os
import sys
import cablepop as cp
import time as ptime
# tstart = ptime.time()

# -------------------------------------------------------------------------
# Get lat/lon indices
#

if os.path.exists(latlon):
    dat = np.genfromtxt(latlon, delimiter=',')
    if len(dat) > 1:
        lats = dat[:, 0]
        lons = dat[:, 1]
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

if ofile is None:  # Default output filename
    ofile = cp.set_output_filename(ifile, '-masked')
fi = nc.Dataset(ifile, 'r')
if 'land' not in fi.variables:
    fi.close()
    print('Land variable must be in input file: ' + ifile)
    import sys
    sys.exit()
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    if 'file_format' in dir(fi):
        fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    else:
        fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')

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
    iilats.append(cp.closest(ilats, lats[ii]))
    iilons.append(cp.closest(ilons, lons[ii]))

# Copy global attributes, adding script
cp.set_global_attributes(fi, fo,
                         add={'history': ptime.asctime() + ': ' +
                              ' '.join(sys.argv)})

# Copy dimensions
cp.create_dimensions(fi, fo)

# create static variables (independent of time)
cp.create_variables(fi, fo, time=False, izip=izip,
                    renamedim={'x': 'lon', 'y': 'lat'})

# create dynamic variables (time dependent)
cp.create_variables(fi, fo, time=True, izip=izip,
                    renamedim={'x': 'lon', 'y': 'lat'})

#
# Copy variables from in to out
#

# copy static variables
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        ovar = fo.variables[ivar.name]
        if ivar.name == 'land':
            ovar[:] = 0
            ovar[iilats, iilons] = 1
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
                ovar[tt, ...] = ivar[tt, ...]

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

# tstop = ptime.time()
# print('Finished in [s]: {:.2f}'.format(tstop-tstart))
