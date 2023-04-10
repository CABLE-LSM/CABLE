#!/usr/bin/env python
"""
usage: sum_patchcasa2d.py [-h] [-o output_netcdf] [-v] [-z] [input_netcdf]

Copy Casa output summing the patches on the same grid point and make lat/lon output.

positional arguments:
  input_netcdf          input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-
                        no_patch2d.nc).
  -v, --verbose         Feedback during copy (default: no feedback).
  -z, --zip             Use netCDF4 variable compression (default: same format
                        as input file).


Example
-------
  python sum_patchcasa2d.py -o cru_out_casa_2009_2011-no_patch2d.nc \
      cru_out_casa_2009_2011.nc


History
-------
Written  Matthias Cuntz, May 2020
             - from sum_patchcasa.py and nopatch2d.py
Modified Matthias Cuntz, Jul 2022
             - flake8 compatible
             - Take into account patchfrac, which is Casa output since rev 9053
         Matthias Cuntz, Apr 2023
             - add x and y variables

"""
import sys
import argparse
import numpy as np
import netCDF4 as nc
import cablepop as cp
import time as ptime

# -------------------------------------------------------------------------
# Command line
#

ofile   = None
verbose = False
izip    = False
parser  = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=('Copy Casa output summing the patches on the same grid point'
                 ' and make lat/lon output.'))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                    help=('output netcdf file name (default:'
                          ' input-no_patch2d.nc).'))
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

#
# Open files
#
fi = nc.Dataset(ifile, 'r')
if verbose:
    print('Input file:  ', ifile)
# check for patchfrac
if 'patchfrac' not in fi.variables:
    fi.close()
    raise ValueError('variable patchfrac not in input file')

if ofile is None:  # Default output filename
    ofile = cp.set_output_filename(ifile, '-no_patch2d')
if verbose:
    print('Output file: ', ofile)
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    if 'file_format' in dir(fi):
        fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    else:
        fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')

# patchfrac with different dimensions
# time is always the first and land the last dimension
patchfrac = fi.variables['patchfrac'][:]
patchfrac1d = patchfrac[0, :]
patchfrac2d = patchfrac[:, :]
patchfrac3d = patchfrac[:, np.newaxis, :]
patchfrac4d = patchfrac[:, np.newaxis, np.newaxis, :]

# Determine number of land points
lats = fi.variables['latitude'][:]  # dim land
lons = fi.variables['longitude'][:]
sidx = np.zeros(lats.size, dtype=int)  # index of new grid cells
lidx = np.ones(lats.size, dtype=int)   # no. of patches per grid cell
nidx = 0
for i in range(1, lons.size):
    if (lons[i] == lons[i - 1]) and (lats[i] == lats[i - 1]):
        # next patch of same grid cell
        lidx[nidx] += 1
    else:
        # new grid cell
        nidx += 1
        sidx[nidx] = i
        lidx[nidx] = 1
sidx = sidx[:nidx]
lidx = lidx[:nidx]

# Get latitude/longitude indices
# input grid
ilats = fi.variables['latitude'][sidx]  # dim no. of grid cells
ilons = fi.variables['longitude'][sidx]
nland = ilats.size
# output grid
dlat = np.abs(np.diff(np.unique(np.sort(ilats)))).min()  # 0.5, 1 degree
dlon = np.abs(np.diff(np.unique(np.sort(ilons)))).min()  # 0.5, 1 degree
nlat = np.rint(150. / dlat).astype(int)  # 300, 150
nlon = np.rint(360. / dlon).astype(int)  # 720, 360
clat = ilats.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
clon = ilons.min() % 1.  # 0.0 or 0.25, 0.0 or 0.5
olat = (-60. + clat + np.arange(nlat) / float(nlat - 1) *
        (150. - dlat))  # new lats
olon = (-180. + clon + np.arange(nlon) / float(nlon - 1) *
        (360. - dlon))  # new lons
olon2d, olat2d = np.meshgrid(olon, olat)  # new lats, lons in 2D
lltree = cp.llKDTree(olat2d, olon2d)  # KD-tree
iidl = np.arange(nland, dtype=int)  # indices of land in input grid
oidx = np.empty(nland, dtype=int)   # indices of lon in output grid
oidy = np.empty(nland, dtype=int)   # indices of lat in output grid
for i in range(nland):
    iy, ix = lltree.query(ilats[i], ilons[i])
    oidx[i] = ix
    oidy[i] = iy

# Copy global attributes, adding script
cp.set_global_attributes(
    fi, fo, add={'history': (ptime.asctime() + ': ' + ' '.join(sys.argv))})

# Copy dimensions
cp.create_dimensions(fi, fo, removedim=['land', 'ntile'],
                     adddim={'x': nlon, 'y': nlat})

# create static variables (independent of time)
cp.create_variables(fi, fo, time=False, izip=izip, fill=True,
                    chunksizes=False,
                    replacedim={'land': ('y', 'x'), 'ntile': ('y', 'x')})

# create x and y variables
if 'x' not in fi.variables:
    xattr = {'standard_name': 'longitude',
             'long_name': 'longitude',
             'units': 'degrees_east',
             'axis': 'X'}
    x = fo.createVariable('x', 'f4', ('x',))
    for a in xattr:
        x.setncattr(a, xattr[a])

if 'y' not in fi.variables:
    yattr = {'standard_name': 'latitude',
             'long_name': 'latitude',
             'units': 'degrees_north',
             'axis': 'Y'}
    y = fo.createVariable('y', 'f4', ('y',))
    for a in yattr:
        y.setncattr(a, yattr[a])

# create dynamic variables (time dependent)
cp.create_variables(fi, fo, time=True, izip=izip, fill=True,
                    chunksizes=False,
                    replacedim={'land': ('y', 'x'), 'ntile': ('y', 'x')})

#
# Copy variables from in to out, summing the patches and make lat/lon
#

# set x and y
if 'x' not in fi.variables:
    fo.variables['x'][:] = olon
if 'y' not in fi.variables:
    fo.variables['y'][:] = olat

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
    if verbose and (n > 0) and ((n % n10) == 0):
        print('  {:d}%'.format(int(n / nvar * 100.)))
    n += 1
    ovar  = fo.variables[ivar.name]
    # read whole field, otherwise execution time increases sharpely
    invar = ivar[:]
    if ('land' in ivar.dimensions) or ('ntile' in ivar.dimensions):
        # sum patch
        if len(invar.shape) == 1:
            ipatchfrac = patchfrac1d
        elif len(invar.shape) == 2:
            ipatchfrac = patchfrac2d
        elif len(invar.shape) == 3:
            ipatchfrac = patchfrac3d
        elif len(invar.shape) == 4:
            ipatchfrac = patchfrac4d
        outshape     = list(invar.shape)
        outshape[-1] = nlat   # patch to y,x
        outshape.append(nlon)
        # fill in memory, then write to disk in one go
        out = np.full(outshape, ovar._FillValue)
        if ivar.name == 'latitude':
            out = olat2d
        elif ivar.name == 'longitude':
            out = olon2d
        elif ivar.name == 'area_gridcell':
            if len(invar.shape) == 1:
                for i in range(nidx):
                    out[oidy[i], oidx[i]] = np.sum(
                        invar[sidx[i]:sidx[i] + lidx[i]])
            else:
                for i in range(nidx):
                    out[..., oidy[i], oidx[i]] = np.sum(
                        invar[..., sidx[i]:sidx[i] + lidx[i]], axis=-1)
        else:
            if len(invar.shape) == 1:
                for i in range(nidx):
                    out[oidy[i], oidx[i]] = np.sum(
                        invar[sidx[i]:sidx[i] + lidx[i]] *
                        ipatchfrac[sidx[i]:sidx[i] + lidx[i]])
            else:
                for i in range(nidx):
                    out[..., oidy[i], oidx[i]] = np.sum(
                        invar[..., sidx[i]:sidx[i] + lidx[i]] *
                        ipatchfrac[..., sidx[i]:sidx[i] + lidx[i]],
                        axis=-1)
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
