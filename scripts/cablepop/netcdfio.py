#!/usr/bin/env python
"""
NetCDF4 helper functions to copy one netcdf file to another while
doing some transformations.

This module was written by Matthias Cuntz while at Institut National de Recherche
pour l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy, France.

It borrows some concepts from the netcdf4 thin layer of David Schaefer.

Copyright (c) 2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Apr 2020 by Matthias Cuntz (mc (at) macu (dot) de)

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   create_dimensions
   create_variables
   set_global_attributes
   set_output_filename
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import netCDF4 as nc


__all__ = ['create_dimensions', 'create_variables',
           'set_global_attributes', 'set_output_filename']


def _tolist(arg):
    """
    Assure that `arg` is a list, e.g. if string or None are given.

    Parameters
    ----------
    arg :
        argument to make list

    Returns
    -------
    list
        list(arg)

    Examples
    --------
    >>> _tolist('string')
    ['string']
    >>> _tolist([1,2,3])
    [1, 2, 3]
    >>> _tolist(None)
    [None]
    """
    if isinstance(arg, str):
        return [arg]
    try:
        return list(arg)
    except TypeError:
        return [arg]


def _get_variable_definition(ncvar):
    """
    Collect information on input variable.

    Parameters
    ----------
    ncvar : netcdf4 variable
        variable of input file

    Returns
    -------
    dict
        containing information on input variable withkey/value pairs.

        The following keys are returned: 'name', 'dtype', 'dimensions', 'fill_vallue', 'chunksizes'

    Examples
    --------
    _get_variable_definition(fi.variables['GPP'])
    """
    out = ncvar.filters() if ncvar.filters() else {}
    # chunksizes
    chunks = list(ncvar.chunking()) if not isinstance(ncvar.chunking(), str) else None
    # missing value
    if "missing_value" in dir(ncvar):
        ifill = ncvar.missing_value
    elif "_FillValue" in dir(ncvar):
        ifill = ncvar._FillValue
    else:
        ifill = None
    # output variable
    out.update({
        "name"       : ncvar.name,
        "dtype"      : ncvar.dtype,
        "dimensions" : list(ncvar.dimensions),
        "fill_value" : ifill,
        "chunksizes" : chunks,
        })
    return out


def create_dimensions(fi, fo, removedim=[], renamedim={}, changedim={}, adddim={}):
    """
    Create dimensions in output from dimensions in input file.

    Parameters
    ----------
    fi : file_handle
        file_handle of opened netcdf input file
    fo : file_handle
        file_handle of opened netcdf output file
    removedim : list of str, optional
        Do not create dimensions given in `removedim` in output file.
    renamedim : dict, optional
        Rename dimensions in output file compared to input file.
        Dimension names in input file are given as dictionary keys,
        corresponding dimension names of output file are give as dictionary values.
    add : dict, optional
        Add dimension to output file.
        New dimension names are given as dictionary keys and new dimension sizes
        are give as dictionary values.
    replacedim : dict, optional
        Replace dimensions in output file compared to input file.
        Dimension names in input file are given as dictionary keys,
        corresponding dimension names of output file are give as dictionary values.

    Returns
    -------
    nothing
        the output file will have the altered or unaltered dimensions
        of the input file.

    Examples
    --------
    create_dimensions(fi, fo, removedim=['patch'], renamedim={'x':'lon', 'y':'lat'}, changedim={'mland':1})
    """
    removedim = _tolist(removedim)
    for d in fi.dimensions.values():
        # remove dimension if in removedim
        if d.name not in removedim:
            # change dimension size if in changedim
            if d.name in changedim.keys():
                nd = changedim[d.name]
            elif d.isunlimited():
                nd = None
            else:
                nd = len(d)
            # rename dimension if in renamedim
            if d.name in renamedim.keys():
                oname = renamedim[d.name]
            else:
                oname = d.name
            # create dimension
            fo.createDimension(oname, nd)
    # add new dimensions
    for d in adddim.keys():
        if d not in fo.dimensions:
            fo.createDimension(d, adddim[d])



def create_variables(fi, fo, time=None, izip=False, fill=None, chunksizes=True,
                     removevar=[], renamevar={}, removedim=[], renamedim={}, replacedim={}):
    """
    Create variables in output from variables in input file.

    Parameters
    ----------
    fi : file_handle
        file handle of opened netcdf input file
    fo : file_handle
        file handle of opened netcdf output file
    time : None or bool, optional
        None:  create all variables (default).

        True:  copy only variables having dimension 'time'.

        False: copy only variables that do not have dimension 'time'.
    izip : bool, optional
        True: the data will be compressed in the netCDF file using gzip compression (default: False).
    fill : float, bool or None, optional
        Determine the behaviour if variable have no _FillValue or missing_value.

        If None or False: no _FillValue will be set.

        If True: _FillValue will be set to default value of the Python packege netCDF4 for this type.

        If number: _FillValue will be set to number.
    chunksizes : bool, optional
        True: include possible chunksizes in output file (default).

        False: do not include chunksize information from input file in output file.
        Set to False, for example, if dimension size gets changed because the chunksize
        on a dimension can not be greater than the dimension size.
    removevar : list of str, optional
        do not create variables given in `removevar` in output file.
    renamevar : dict, optional
        rename variables in output file compared to input file.
        Variable names in input file are given as dictionary keys,
        corresponding variable names of output file are give as dictionary values.
    removedim : list of str, optional
        Remove dimensions from variable definitions in output file.
    renamedim : dict, optional
        Rename dimensions for variables in output file.
        Dimension names in input file are given as dictionary keys,
        corresponding dimension names of output file are give as dictionary values.
    replacedim : dict, optional
        Replace dimensions for variables in output file.
        Dimension names in input file are given as dictionary keys,
        corresponding dimension names of output file are give as dictionary values.

        The output names can be tuples or lists to extend dimensions of a variable.

    Returns
    -------
    nothing
        the output file will have the altered or unaltered variables
        of the input file defined.

    Examples
    --------
    create_variable(fi, fo, fill=True, izip=True, removedim=['patch'], renamevar={'lon':'longitude'}, replacedim={'land':('y','x')})
    """
    removevar = _tolist(removevar)
    removedim = _tolist(removedim)
    for ivar in fi.variables.values():
        # remove variable if in removevar
        if ivar.name not in removevar:
            if time is None:
                itime = True
            else:
                if time:
                    itime = 'time' in ivar.dimensions
                else:
                    itime = 'time' not in ivar.dimensions
            if itime:
                invardef = _get_variable_definition(ivar)
                if izip: invardef.update({'zlib':True})
                # rename variable if in renamevar
                if ivar.name in renamevar.keys():
                    invardef['name'] = renamevar[ivar.name]
                # remove dimension if in removedim
                dims   = invardef['dimensions']
                chunks = invardef['chunksizes']
                for dd in removedim:
                    if dd in dims:
                        ip = dims.index(dd)
                        dims.remove(dd)
                        if chunks:
                            _ = chunks.pop(ip)
                # rename dimension if in renamedim
                for dd in renamedim.keys():
                    if dd in dims:
                        dims[dims.index(dd)] = renamedim[dd]
                # replace dimensions if in replacedim
                for dd in replacedim.keys():
                    if dd in dims:
                        rdim = _tolist(replacedim[dd])
                        ip = dims.index(dd)
                        dims = dims[:ip] + rdim + dims[ip+1:]
                        rchunk = []
                        for cc in rdim:
                            rchunk.append(len(fo.dimensions[cc])) # use fo not fi because new dims perhaps not yet in fi
                        chunks = chunks[:ip] + rchunk + chunks[ip+1:]
                invardef['dimensions'] = dims
                invardef['chunksizes'] = chunks
                # set missing value if None
                if invardef['fill_value'] is None:
                    if fill:
                        if isinstance(fill, bool):
                            if invardef['dtype'] == np.dtype(np.int8):
                                invardef['fill_value'] = nc.default_fillvals['i1']
                            elif invardef['dtype'] == np.dtype(np.int16):
                                invardef['fill_value'] = nc.default_fillvals['i2']
                            elif invardef['dtype'] == np.dtype(np.int32):
                                invardef['fill_value'] = nc.default_fillvals['i4']
                            elif invardef['dtype'] == np.dtype(np.int64):
                                invardef['fill_value'] = nc.default_fillvals['i8']
                            elif invardef['dtype'] == np.dtype(np.uint8):
                                invardef['fill_value'] = nc.default_fillvals['u1']
                            elif invardef['dtype'] == np.dtype(np.uint16):
                                invardef['fill_value'] = nc.default_fillvals['u2']
                            elif invardef['dtype'] == np.dtype(np.uint32):
                                invardef['fill_value'] = nc.default_fillvals['u4']
                            elif invardef['dtype'] == np.dtype(np.uint64):
                                invardef['fill_value'] = nc.default_fillvals['u8']
                            elif invardef['dtype'] == np.dtype(np.float32):
                                invardef['fill_value'] = nc.default_fillvals['f4']
                            elif invardef['dtype'] == np.dtype(np.float64):
                                invardef['fill_value'] = nc.default_fillvals['f8']
                            else:
                                import warnings
                                warnings.warn("dtype of variable "+invardef['name']+" unknown: "+invardef['dtype'])
                        else:
                            invardef['fill_value'] = fill
                # exclude chunksizes
                if not chunksizes:
                    _ = invardef.pop('chunksizes')
                oname = invardef.pop("name")
                otype = invardef.pop("dtype")
                ovar = fo.createVariable(oname, otype, **invardef)
                for k in ivar.ncattrs():
                    iattr = ivar.getncattr(k)
                    if (k != 'missing_value') and (k != '_FillValue'):
                        ovar.setncattr(k, iattr)


def set_global_attributes(fi, fo, add={}):
    """
    Create global output file attributes from input global file attributes.

    Parameters
    ----------
    fi : file_handle
        file_handle of opened netcdf input file
    fo : file_handle
        file_handle of opened netcdf output file
    add : dict, optional
        dict values will be given to attributes given in dict keys.

        Attributes will be created if they do not exist yet.

    Returns
    -------
    nothing
        output will have global file attributes

    Examples
    --------
    set_global_attributes(fi, fo, add={'history':time.asctime()+': '+' '.join(sys.argv)})
    """
    for k in fi.ncattrs():
        iattr = fi.getncattr(k)
        # add to existing global attribute
        if k in add.keys():
            iattr += '\n'+add[k]
        fo.setncattr(k, iattr)
    # add if global attribute does not exist yet
    for k in add.keys():
        if k not in fi.ncattrs():
            fo.setncattr(k, add[k])


def set_output_filename(ifile, ext):
    """
    Create output file name from input file name by adding `ext` before
    the file suffix.

    Parameters
    ----------
    ifile : str
        input file name
    ext : str
        string to add before file suffix

    Returns
    -------
    str
        output filename with ext before file suffix

    Examples
    --------
    >>> set_output_filename('in.nc', '-no_patch')
    in-no_patch.nc
    >>> set_output_filename('in.nc', '.nop'')
    in.nop.nc
    """
    sifile = ifile.split('.')
    sifile[-2] = sifile[-2]+ext
    ofile = '.'.join(sifile)
    return ofile


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
