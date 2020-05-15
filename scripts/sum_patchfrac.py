#!/usr/bin/env python
from __future__ import print_function
"""
usage: sum_patchfrac.py [-h] [-o output_netcdf] [-v] [-z] [input_netcdf]

Copy Cable output summing variables weighted with patchfrac.

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
  python sum_patchfrac.py -o cru_out_cable_2009_2011-no_patch.nc cru_out_cable_2009_2011.nc


History
-------
Written  Matthias Cuntz, Oct 2018
Modified Matthias Cuntz, Apr 2020 - set _FillValue when creating variables
                                  - script more general for all Cable-POP output files
         Matthias Cuntz, Apr 2020 - use python module cablepop
"""

# -------------------------------------------------------------------------
# Command line
#

import argparse

ofile   = None
verbose = False
izip    = False
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''Copy Cable output summing variables weighted with patchfrac.'''))
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
import time as ptime
if verbose:
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
if 'patchfrac' not in fi.variables:
    fi.close()
    print('No variable patchfrac in input file '+ifile)
    import sys
    sys.exit()
# Check patch dimensions 1
patchfrac = fi.variables['patchfrac']
if 'patch' in patchfrac.dimensions:
    pname = 'patch'
elif 'mp' in patchfrac.dimensions:
    pname = 'mp'
else:
    fi.close()
    print("Did not understand patchfrac: neither 'patch' nor 'mp' in dimensions - "+ifile)
    import sys
    sys.exit()

if verbose: print('Input file:  ', ifile)
if verbose: print('Output file: ', ofile)
if izip:
    fo = nc.Dataset(ofile, 'w', format='NETCDF4')
else:
    if 'file_format' in dir(fi):
        fo = nc.Dataset(ofile, 'w', format=fi.file_format)
    else:
        fo = nc.Dataset(ofile, 'w', format='NETCDF3_64BIT_OFFSET')

# Check patch dimensions 2
pidx   = patchfrac.dimensions.index(pname)   # positive index from front
npatch = patchfrac.shape[pidx]
if 'time' in patchfrac.dimensions:
    pidxt = True
else:
    pidxt = False
ntime  = fi.dimensions['time'].size

# Copy global attributes, adding script
cp.set_global_attributes(fi, fo, add={'history':ptime.asctime()+': '+' '.join(sys.argv)})

# Copy dimensions
cp.create_dimensions(fi, fo, removedim=[pname])

# create static variables (independent of time)
cp.create_variables(fi, fo, time=False, izip=izip, removedim=pname)

# create dynamic variables (time dependent)
cp.create_variables(fi, fo, time=True, izip=izip, removedim=pname)

#
# Copy variables from in to out, summing the patches
#

# copy static variables
if verbose: print('Copy static variables')
if pidxt:
    pfrac  = patchfrac[0,...]
else:
    pfrac  = patchfrac
np.seterr(over='ignore', under='ignore', invalid='ignore') # multiply masked elements
for ivar in fi.variables.values():
    if 'time' not in ivar.dimensions:
        ovar  = fo.variables[ivar.name]
        invar = ivar[:]
        if pname in ivar.dimensions:
            idx = ivar.dimensions.index(pname)
            if ivar.name == 'iveg':
                # take first patch, could also be patch with maximum frac
                out = np.delete(invar, range(1,npatch), axis=idx).squeeze()
            else:
                out = np.sum(invar*pfrac, axis=idx) # numpy broadcasting
        else:
            out = invar
        ovar[:] = out

# copy dynamic variables
if verbose:
    print('Copy dynamic variables')
if ntime < 20:
    t10 = 1
else:
    t10 = ntime // 10
for tt in range(ntime):
    if verbose and (tt > 0) and ((tt % t10) == 0): print('  {:d}%'.format(int(tt/ntime*100.)))
    if pidxt:
        pfrac  = patchfrac[tt,...]
    else:
        pfrac  = patchfrac
    for ivar in fi.variables.values():
        if 'time' in ivar.dimensions:
            ovar = fo.variables[ivar.name]
            invar = ivar[tt,...]
            if pname in ivar.dimensions:
                idx = ivar.dimensions.index(pname) - 1 # -1 because time removed at beginning
                if ivar.name == 'iveg':
                    # take first patch, could also be patch with maximum frac
                    out = np.delete(invar, range(1,npatch), axis=idx).squeeze()
                else:
                    out = np.sum(invar*pfrac, axis=idx) # numpy broadcasting
            else:
                out = invar
            ovar[tt,...] = out

fi.close()
fo.close()

# -------------------------------------------------------------------------
# Finish

if verbose:
    tstop = ptime.time()
    print('Finished in [s]: {:.2f}'.format(tstop-tstart))
