#!/usr/bin/env python
from __future__ import print_function
"""
usage: extract_restart_point.py [-h] [-l lat,lon] [-z] [netcdf [netcdf ...]]

Create mask for Cable with only latitude,longitude set to 1 on given latitudes,longitudes, .
Extract the latitudes,longitudes from any number of Cable-POP restart files.

Output files will be named input_name-extract.nc.

positional arguments:
  netcdf                input netcdf files (>=1): mask [[restart_file1
                        [restart_file2 ...]]

optional arguments:
  -h, --help            show this help message and exit
  -l lat,lon, --latlon lat,lon
                        latitude,longitude of grid point to mask/extract, or
                        file with several latitude,longitude coordinates.
  -z, --zip             Use netCDF4 variable compression (default: same format
                        as input file).


Examples
--------
  # extract Hawaii
  python extract_restart_point.py glob_ipsl_1x1.nc *_rst_*.nc pop_*_ini_*.nc -l 19.8968,-155.5828

  # extract on Southern hemisphere
  python extract_restart_point.py glob_ipsl_1x1.nc *_rst_*.nc pop_*_ini_*.nc -l='-12.5,-67.5'


History
-------
Written  Matthias Cuntz, Oct 2018
Modified Matthias Cuntz, Apr 2020 - updated to same structure as other netcdf utilities
                                  - allow several lats/lons via input file
         Matthias Cuntz, Apr 2020 - use python module cablepop
"""

# -------------------------------------------------------------------------
# Command line
#

import argparse

latlon = None
izip   = False
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=('''Create mask for Cable with only latitude,longitude set to 1 on given latitudes,longitudes, .\nExtract the latitudes,longitudes from any number of Cable-POP restart files.\n\nOutput files will be named input_name-extract.nc.\n'''))
parser.add_argument('-l', '--latlon', action='store',
                    default=latlon, dest='latlon', metavar='lat,lon',
                        help='latitude,longitude of grid point to mask/extract, or file with several latitude,longitude coordinates.')
parser.add_argument('-z', '--zip', action='store_true', default=izip, dest='izip',
                    help='Use netCDF4 variable compression (default: same format as input file).')
parser.add_argument('ifiles', nargs='*', default=None, metavar='netcdf',
                   help='input netcdf files (>=1): mask [[restart_file1 [restart_file2 ...]]')
args   = parser.parse_args()
latlon = args.latlon
izip   = args.izip
ifiles = args.ifiles
del parser, args

if latlon is None:
    raise IOError('latitude,longitude or file with latitudes,longitudes must be given with -l flag.')

if ifiles is None:
    raise IOError('At least mask input file needed. Input files: mask_file restart_file1 restart_file2 ...')
elif len(ifiles) == 1:
    restartfiles = []
else:
    restartfiles = ifiles[1:]
maskfile = ifiles[0]

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
        lats = dat[:,0]
        lons = dat[:,1]
    else:
        lats = np.array([dat[0]])
        lons = np.array([dat[1]])
else:
    lats, lons = latlon.split(',')
    lats = np.array([float(lats)])
    lons = np.array([float(lons)])
    # for old code
    ilat, ilon = [ float(i) for i in latlon.split(',') ]
nlatlon = lats.size


# -------------------------------------------------------------------------
# Create new mask file
#

ifile = maskfile
ofile = cp.set_output_filename(ifile, '-extract')
print('Create ', ofile)
fi = nc.Dataset(ifile, 'r')
if 'land' not in fi.variables:
    fi.close()
    print('Land variable must be in first input, which must be the mask file: '+ifile)
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
if 'latitude' in fi.variables:
    latname = 'latitude'
else:
    latname = 'lat'
if 'longitude' in fi.variables:
    lonname = 'longitude'
else:
    lonname = 'lon'
inlats = fi.variables[latname][:]
inlons = fi.variables[lonname][:]
iinlats = []
iinlons = []
for ii in range(nlatlon):
    iinlats.append(cp.closest(inlats, lats[ii]))
    iinlons.append(cp.closest(inlons, lons[ii]))
    
# Copy global attributes, adding script
cp.set_global_attributes(fi, fo, add={'history':ptime.asctime()+': '+' '.join(sys.argv)})

# Copy dimensions
cp.create_dimensions(fi, fo)

# create static variables (independent of time)
cp.create_variables(fi, fo, time=False, izip=izip, chunksizes=False)

# create dynamic variables (time dependent)
cp.create_variables(fi, fo, time=True, izip=izip, chunksizes=False)

# copy static variables
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        ovar = fo.variables[ivar.name]
        if ivar.name == 'land':
            ovar[:] = 0
            ovar[iinlats,iinlons] = 1
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
# Extract restart files
#

for ifile in restartfiles:
    ofile = cp.set_output_filename(ifile, '-extract')
    print('Extract ', ofile)
    fi = nc.Dataset(ifile, 'r')

    if 'latitude' in fi.variables:
        latname = 'latitude'
    else:
        latname = 'lat'
    if 'longitude' in fi.variables:
        lonname = 'longitude'
    else:
        lonname = 'lon'
    if 'mp' in fi.dimensions.keys():
        # cable file is different to other restart files casa, LUC, climate, c13o2
        #   lats/lons are given per grid cells but most output is given per tile
        #   -> make lats/lons per tile with the help of number of tiles per grid cell 'nap'
        mlats = fi.variables[latname][:]
        mlons = fi.variables[lonname][:]
        nap   = fi.variables['nap'][:].astype(np.int)
        ntile = fi.dimensions['mp'].size
        ngrid = fi.dimensions['mland'].size
        tlats = np.full(ntile, -999.)
        tlons = np.full(ntile, -999.)
        z = 0
        for i in range(ngrid):
            tlats[z:z+nap[i]] = mlats[i]
            tlons[z:z+nap[i]] = mlons[i]
            z += nap[i]
        # indices of grid cells lats/lons on grid cell basis
        miiland = []
        for ii in range(nlatlon):
            miilat0 = cp.closest(mlats, lats[ii])
            miilon0 = cp.closest(mlons, lons[ii])
            miiland.extend(np.where((np.array(mlats) == mlats[miilat0]) & (np.array(mlons) == mlons[miilon0]))[0])
        if len(miiland) == 0:
            print('-> no land point.')
            fi.close()
            continue
    else:
        tlats = fi.variables[latname][:]
        tlons = fi.variables[lonname][:]
        miiland = []
    # indices of grid cells lats/lons on tile basis
    tiiland = []
    for ii in range(nlatlon):
        tiilat0 = cp.closest(tlats, lats[ii])
        tiilon0 = cp.closest(tlons, lons[ii])
        tiiland.extend(np.where((np.array(tlats) == tlats[tiilat0]) & (np.array(tlons) == tlons[tiilon0]))[0])
    if len(tiiland) == 0:
        print('-> no land point or forested area.')
        fi.close()
        continue

    # Create output file
    if izip:
        fo = nc.Dataset(ofile, 'w', format='NETCDF4')
    else:
        if 'file_format' in dir(fi):
            fo = nc.Dataset(ofile, 'w', format=fi.file_format)
        else:
            fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')
    
    # Copy global attributes, adding script
    cp.set_global_attributes(fi, fo, add={'history':ptime.asctime()+': '+' '.join(sys.argv)})

    # Copy dimensions
    cp.create_dimensions(fi, fo, changedim={'mp':len(tiiland), 'land':len(tiiland), 'nland':len(tiiland), 'mland':len(miiland)})

    # create static variables (independent of time)
    cp.create_variables(fi, fo, time=False, izip=izip, chunksizes=False)

    # create dynamic variables (time dependent)
    cp.create_variables(fi, fo, time=True, izip=izip, chunksizes=False)

    # Copy variables, selecting the land points with latitude,longitude    
    for ivar in fi.variables.values():
        ovar = fo.variables[ivar.name]
        if ('mp' in ivar.dimensions) or ('land' in ivar.dimensions) or ('nland' in ivar.dimensions):
            if len(ivar.shape) == 1:
                ovar[:] = ivar[tiiland]
            else:
                ovar[:] = ivar[...,tiiland]
        elif ('mland' in ivar.dimensions):
            if len(ivar.shape) == 1:
                ovar[:] = ivar[miiland]
            else:
                ovar[:] = ivar[...,miiland]
        else:
            ovar = fo.variables[ivar.name]
            ovar[:] = ivar[:]

    # Close restart files
    fi.close()
    fo.close()

# -------------------------------------------------------------------------
# Finish

# tstop = ptime.time()
# print('Finished in [s]: {:.2f}'.format(tstop-tstart))
