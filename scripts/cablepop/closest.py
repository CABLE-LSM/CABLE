#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
closest : Index of an array (arr) at which the entry is closest to a given number (num),
          which is `argmin(abs(arr-num))`.

This module was by Matthias Cuntz while at Department of Computational
Hydrosystems, Helmholtz Centre for Environmental Research - UFZ,
Leipzig, Germany, and continued while at Institut National de Recherche
pour l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy, France.

Copyright (c) 2012-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jan 2012 by Matthias Cuntz (mc (at) macu (dot) de)
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Make numpy doctsring format, Apr 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   closest
"""
import numpy as np

def closest(arr, num, value=False):
    """
    Index of an array (arr) at which the entry is closest to a given number (num),
    which is `argmin(abs(arr-num))`.

    Parameters
    ----------
    arr : array_like
        Array to search closest entry
    num : number
        Number to which the closest entry is search in arr
    value : bool, optional
        If true, give closest array element instead of index (default: False)

    Returns
    -------
    index : int
        Index of element closest to given number in flattend array.
        Use np.unravel_index to get index tuple.

    Examples
    --------
    >>> arr = np.arange(100)/99.*5.
    >>> print(closest(arr, 3.125))
    62
    >>> out = closest(arr, 3.125, value=True)
    >>> print('{:.3f}'.format(out))
    3.131

    >>> arr = np.arange(100).reshape((10,10))/99.*5.
    >>> out = closest(arr, 3.125, value=True)
    >>> print('{:.3f}'.format(out))
    3.131
    >>> print(closest(arr, 3.125))
    62
    >>> print(np.unravel_index(closest(arr, 3.125), arr.shape))
    (6, 2)
    >>> out = arr[np.unravel_index(closest(arr, 3.125), arr.shape)]
    >>> print('{:.3f}'.format(out))
    3.131

    History
    -------
    Written,  MC, Jan 2012
    Modified, MC, Feb 2013 - ported to Python 3
              MC, Apr 2020 - numpy docstring format
    """
    out = np.ma.argmin(np.ma.abs(np.ma.array(arr)-num))
    if value:
      return arr.flat[out]
    else:
      return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
