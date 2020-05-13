#!/usr/bin/env python
"""
Accessing netCDF Data by Coordinates

This is from a blog post by Unidata how to quickly access netCDF data based on
geospatial coordinates instead of array indices:
    https://www.unidata.ucar.edu/blogs/developer/en/entry/accessing_netcdf_data_by_coordinates
The code is take from the accompanying notebook at:
    https://github.com/Unidata/python-workshop/blob/fall-2016/notebooks/netcdf-by-coordinates.ipynb

Imagine Longitude and Latitude are 2-D netCDF variables indexed by Y and X
dimensions.

There's no direct way to compute grid indexes from coordinates via a coordinate
system and projection parameters. Instead, we have to rely on the latitude and
longitude auxiliary coordinate variables, as required by the CF conventions for
data not on a simple lat,lon grid.

For a given lat/lon, e.g. 50N/140W, we need to find Y and X indexes iy and ix
such that (Longitude[iy, ix], Latitude[iy, ix]) is "close" to (50.0, -140.0).

For a flat Earth, this would be minimizing the distance:
(Latitude[iy, ix] - lat0)**2 + (Longitude[iy, ix] - lon0)**2

This code takes a spherical Earth and the "tunnel distance". It also flexible
if longitudes are from 0 to 360 or -180 to 180.

The code uses nearest-neighbour queries: (the Cython version of) KDTree:
scipy.spatial.cKDTree

Examples
--------
We have n land points with lats and lons.

We want to bring them on a regular 2-D grid of 1.0 degree.
nlat = 180
nlon = 360
dlat = 0.5
dlon = 0.5
olat = -90.  + dlat + np.arange(nlat)/float(nlat-1) * (180.-dlat) # new lats
olon = -180. + dlon + np.arange(nlon)/float(nlon-1) * (360.-dlon) # new lons

One builds the tree with the target lat/lon grid (2-D grids):
olon2d, olat2d = np.meshgrid(olon, olat) # new lats, lons in 2D
lltree = llKDTree(olat2d, olon2d) # KD-tree

The input indices are simply 0, ..., n-1:
iidl = np.arange(n, dtype=np.int) # indices of land in input grid

Now we can walk through the input indices and query the tree for the
corresponding output indices:
oidx = np.empty(n, dtype=np.int)  # indices of lon in output grid
oidy = np.empty(n, dtype=np.int)  # indices of lat in output grid
for i in iidl:
    iy, ix = lltree.query(lats[i], lons[i])
    oidx[i] = ix
    oidy[i] = iy

A 1-D field invar with land points can now be transformed to 2-D with:
out = np.full((nlat,nlon), np.nan)
out[oidy, oidx] = invar[iidl]

A n-D field can be transformed like:
outshape = list(invar.shape)
outshape[-1] = nlat   # land to y,x
outshape.append(nlon)
out = np.full(outshape, np.nan)
out[..., oidy, oidx] = invar[..., iidl]

This module was by Matthias Cuntz while at Institut National de Recherche pour
l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy, France.

Copyright (c) 2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written May 2020 by Matthias Cuntz (mc (at) macu (dot) de)

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   llKDTree
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from scipy.spatial import cKDTree

class llKDTree(object):
    """
    Lat/Lon-KD-Tree

    Accessing data by geospatial coordinates instead of array indices.

    This is from a blog post by Unidata how to quickly access netCDF data based on
    geospatial coordinates instead of array indices:
        https://www.unidata.ucar.edu/blogs/developer/en/entry/accessing_netcdf_data_by_coordinates
    The code is take from the accompanying notebook at:
        https://github.com/Unidata/python-workshop/blob/fall-2016/notebooks/netcdf-by-coordinates.ipynb

    For a given lat/lon, e.g. 50N/140W, we need to find Y and X indexes iy and ix
    such that (Longitude[iy, ix], Latitude[iy, ix]) is "close" to (50.0, -140.0).

    This code takes a spherical Earth and the "tunnel distance". It also flexible
    if longitudes are from 0 to 360 or -180 to 180.

    The code uses nearest-neighbour queries: (the Cython version of) KDTree:
    scipy.spatial.cKDTree

    Examples
    --------
    We have n land points with lats and lons.

    We want to bring them on a regular 2-D grid of 1.0 degree.
    nlat = 180
    nlon = 360
    dlat = 0.5
    dlon = 0.5
    olat = -90.  + dlat + np.arange(nlat)/float(nlat-1) * (180.-dlat) # new lats
    olon = -180. + dlon + np.arange(nlon)/float(nlon-1) * (360.-dlon) # new lons

    One builds the tree with the target lat/lon grid (2-D grids):
    olon2d, olat2d = np.meshgrid(olon, olat) # new lats, lons in 2D
    lltree = llKDTree(olat2d, olon2d) # KD-tree

    The input indices are simply 0, ..., n-1:
    iidl = np.arange(n, dtype=np.int) # indices of land in input grid

    Now we can walk through the input indices and query the tree for the
    corresponding output indices:
    oidx = np.empty(n, dtype=np.int)  # indices of lon in output grid
    oidy = np.empty(n, dtype=np.int)  # indices of lat in output grid
    for i in iidl:
        iy, ix = lltree.query(lats[i], lons[i])
        oidx[i] = ix
        oidy[i] = iy

    A 1-D field invar with land points can now be transformed to 2-D with:
    out = np.full((nlat,nlon), np.nan)
    out[oidy, oidx] = invar[iidl]

    A n-D field can be transformed like:
    outshape = list(invar.shape)
    outshape[-1] = nlat   # land to y,x
    outshape.append(nlon)
    out = np.full(outshape, np.nan)
    out[..., oidy, oidx] = invar[..., iidl]


    History
    -------
    Written,  Unidata,        Aug 2013
    Modified, Matthias Cuntz, May 2020 - module
    """
    def __init__(self, latvar, lonvar):
        rad_factor = np.pi/180.0 # for trigonometry, need angles in radians
        self.latvals = latvar * rad_factor
        self.lonvals = lonvar * rad_factor
        self.shape   = self.latvals.shape
        clat, clon = np.cos(self.latvals), np.cos(self.lonvals)
        slat, slon = np.sin(self.latvals), np.sin(self.lonvals)
        clat_clon  = clat*clon
        clat_slon  = clat*slon
        triples  = list(zip(np.ravel(clat*clon), np.ravel(clat*slon), np.ravel(slat)))
        self.kdt = cKDTree(triples) # build tree during initialisation

    def query(self, lat0, lon0):
        rad_factor = np.pi/180.0 
        lat0_rad = lat0 * rad_factor
        lon0_rad = lon0 * rad_factor
        clat0, clon0 = np.cos(lat0_rad), np.cos(lon0_rad)
        slat0, slon0 = np.sin(lat0_rad), np.sin(lon0_rad)
        dist_sq_min, minindex_1d = self.kdt.query([clat0*clon0, clat0*slon0, slat0]) # flattened index
        iy_min, ix_min = np.unravel_index(minindex_1d, self.shape) # 2d index
        return iy_min, ix_min


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
