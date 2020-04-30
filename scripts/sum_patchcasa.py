#!/usr/bin/env python
from __future__ import print_function
"""
usage: sum_patchcasa.py [-h] [-o output_netcdf] [-v] [-z] [input_netcdf]

Copy Casa output summing the patches on the same grid point.

positional arguments:
  input_netcdf          input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o output_netcdf, --outfile output_netcdf
                        output netcdf file name (default: input-no_patch.nc).
  -v, --verbose         Feedback during copy (default: no feedback).
  -z, --zip             Use netCDF4 variable compression (default: same format
                        as input file).


Example
-------
  python sum_patchcasa.py -o cru_out_casa_2009_2011-no_patch.nc cru_out_casa_2009_2011.nc


History
-------
Written  Matthias Cuntz, Apr 2020 - from sum_patchfrac.py and casa2d.py
"""

# -------------------------------------------------------------------------
# Command line
#

import argparse

ofile   = None
verbose = False
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Copy Casa output summing the patches on the same grid point.'''))
parser.add_argument('-o', '--outfile', action='store',
                    default=ofile, dest='ofile', metavar='output_netcdf',
                        help='output netcdf file name (default: input-no_patch.nc).')
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
import numpy as np
import netCDF4 as nc
import cablepop as cp
if verbose:
    import time as ptime
    tstart = ptime.time()

# -------------------------------------------------------------------------
# Copy data
#

#
# Open files
#
if ofile is None: # Default output filename
    ofile = cp.set_output_filename(ifile, '-no_patch')
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

# Determine number of land points
lats = fi.variables['latitude'][:]
lons = fi.variables['longitude'][:]
sidx = np.zeros(lats.size, dtype=np.int)
lidx = np.ones(lats.size,  dtype=np.int)
nidx = 0
for i in range(1,lons.size): 
   if (lons[i] == lons[i-1]) and (lats[i] == lats[i-1]): 
       lidx[nidx] += 1
   else:
       nidx += 1
       sidx[nidx] = i
       lidx[nidx] = 1
sidx = sidx[:nidx]
lidx = lidx[:nidx]
    
# Copy global attributes, adding script
cp.set_global_attributes(fi, fo, add={'history':ptime.asctime()+': '+' '.join(sys.argv)})

# Copy dimensions
cp.create_dimensions(fi, fo, change={'land':nidx, 'ntile':nidx})

# create static variables (independent of time)
cp.create_variables(fi, fo, time=False, izip=izip, fill=-1e33, chunksizes=False)

# create dynamic variables (time dependent)
cp.create_variables(fi, fo, time=True, izip=izip, fill=-1e33, chunksizes=False)

#
# Copy variables from in to out, summing the patches
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
        outshape[-1] = nidx
        # fill in memory, then write to disk in one go
        out = np.full(outshape, ovar._FillValue)
        if len(outshape) == 1:
            for i in range(nidx):
                out[i] = invar[sidx[i]:sidx[i]+lidx[i]].sum()
        else:
            for i in range(nidx):
                out[...,i] = invar[...,sidx[i]:sidx[i]+lidx[i]].sum(axis=-1)
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
