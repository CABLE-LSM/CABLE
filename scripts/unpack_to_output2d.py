#!/usr/bin/env python3
"""
usage: unpack_to_output2d.py [-h] [-o output_netcdf] [-v] [-z] [input_netcdf]

Transform Cable or Casa land format to 2D arrays with latitude/longitude.

positional arguments:
  input_netcdf          input netcdf file.

options:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-2d.nc).
  -v, --verbose         Feedback during copy (default: no feedback).
  -z, --zip             Use netCDF4 variable compression
                        (default: same format as input file).


Example
-------
  python pack_output2d.py -z -o cru_out_casa_2009_2011-packed.nc \
      cru_out_casa_2009_2011.nc
  python unpack_to_output2d.py -z -o cru_out_casa_2009_2011-2d.nc \
      cru_out_casa_2009_2011-packed.nc


History
-------
Written  Matthias Cuntz, May 2020
             - from nopatch2d.py

"""
import argparse
import sys
import numpy as np
import netCDF4 as nc
import cablepop as cp
import pyjams as pj
import time as ptime

# -------------------------------------------------------------------------
# Command line
#

ofile   = None
verbose = False
izip    = False
parser  = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('Transform Cable or Casa land format to 2D arrays with'
                 ' latitude/longitude.') )
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                    help='output netcdf file name (default: input-2d.nc).')
parser.add_argument('-v', '--verbose', action='store_true', default=verbose,
                    dest='verbose',
                    help='Feedback during copy (default: no feedback).')
parser.add_argument('-z', '--zip', action='store_true', default=izip,
                    dest='izip',
                    help=('Use netCDF4 variable compression (default:'
                          ' same format as input file).'))
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

if verbose:
    tstart = ptime.time()

# -------------------------------------------------------------------------
# Copy data
#

# Input file
fi = nc.Dataset(ifile, 'r')
if verbose:
    print('Input file:  ', ifile)

# Output file
if ofile is None:  # Default output filename
    ofile = pj.ncio.set_output_filename(ifile, '-2d')
if verbose:
    print('Output file: ', ofile)
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
if 'x' in fi.dimensions:  # use existing grid
    olat = fi.variables['y'][:]
    olon = fi.variables['x'][:]
    nlat   = olat.size
    nlon   = olon.size
else:                     # create new global grid -60 to +90 latitudes
    # 0.5, 1 degree
    dlon = np.abs(np.diff(np.unique(np.sort(ilons)))).min()
    if np.all(ilats == ilats[0]):
        # if only few cells on one latitude band
        dlat = dlon
    else:
        # 0.5, 1 degree
        dlat = np.abs(np.diff(np.unique(np.sort(ilats)))).min()
    nlat = np.rint(150. / dlat).astype(int)  # 300, 150
    nlon = np.rint(360. / dlon).astype(int)  # 720, 360
    clat = ilats.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
    clon = ilons.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
    # new lats
    olat = -60.  + clat + np.arange(nlat) / float(nlat - 1) * (150. - dlat)
    olat = olat[::-1]
    # new lons
    olon = -180. + clon + np.arange(nlon) / float(nlon - 1) * (360. - dlon)
olon2d, olat2d = np.meshgrid(olon, olat)  # new lats, lons in 2D
lltree = cp.llKDTree(olat2d, olon2d)  # KD-tree
iidl   = np.arange(nland, dtype=int)  # indices of land in input grid
oidx   = np.empty(nland, dtype=int)   # indices of lon in output grid
oidy   = np.empty(nland, dtype=int)   # indices of lat in output grid
for i in range(nland):
    iy, ix = lltree.query(ilats[i], ilons[i])
    oidx[i] = ix
    oidy[i] = iy

# Copy global attributes, adding this script
pj.ncio.copy_global_attributes(fi, fo,
                               add={'history': ptime.asctime() + ': ' +
                                    ' '.join(sys.argv)})

# Copy dimensions
pj.ncio.copy_dimensions(fi, fo,
                        removedim=['land', 'ntile'],
                        adddim={'x': nlon, 'y': nlat})

# Create static variables (independent of time)
if 'local_lat' in fi.variables:
    renvar = {}  # {'latitude': 'nav_lat', 'longitude': 'nav_lon'}
else:
    renvar = {}
pj.ncio.create_variables(fi, fo, time=False, izip=izip, fill=True,
                         chunksizes=False, renamevar=renvar,
                         replacedim={'land': ('y', 'x'), 'ntile': ('y', 'x')})
# create dynamic variables (time dependent)
pj.ncio.create_variables(fi, fo, time=True, izip=izip, fill=True,
                         chunksizes=False, renamevar=renvar,
                         replacedim={'land': ('y', 'x'), 'ntile': ('y', 'x')})

if 'x' not in fi.variables:
    print('Create x')
    nvar = {'name': 'x',
            'dtype': ilons.dtype,
            'dimensions': ('x'),
            'units': 'degrees_east'}
    ovar = pj.ncio.create_new_variable(nvar, fo)
    ovar[:] = olon

if 'y' not in fi.variables:
    print('Create y')
    nvar = {'name': 'y',
            'dtype': ilats.dtype,
            'dimensions': ('y'),
            'units': 'degrees_north'}
    ovar = pj.ncio.create_new_variable(nvar, fo)
    ovar[:] = olat

#
# Copy variables from in to out expanding the land dimension to y, x
#

# copy static and dynamic variables
if verbose:
    print('Copy variables')
nvar = len(fi.variables.keys())
if nvar < 20:
    n10 = 1
else:
    n10 = nvar // 10
n = 0
for ivar in fi.variables.values():
    if verbose:
        print(f'    {ivar.name}')
    n += 1
    ovar = fo.variables[ivar.name]
    invar = ivar[:]  # read whole field, otherwise times increasing sharpely
    if ('land' in ivar.dimensions) or ('ntile' in ivar.dimensions):
        outshape = list(invar.shape)
        outshape[-1] = nlat   # land to y, x
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
    print('Finished in [s]: {:.2f}'.format(tstop - tstart))
