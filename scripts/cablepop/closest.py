#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np

def closest(vec, num, value=False):
    """
        Index in array at which the entry is closest to a given number.


        Definition
        ----------
        def closest(vec, num, value=False):


        Input
        -----
        vec        ND-array
        num        number to compare


        Optional Input
        --------------
        value      If true, give closest array element instead of index


        Output
        ------
        Flat index of element closest to given number num.
        Use np.unravel_index to get index tuple.


        Restrictions
        ------------
        None.


        Examples
        --------
        >>> vec = np.arange(100)/99.*5.
        >>> print(closest(vec, 3.125))
        62
        >>> from autostring import astr
        >>> print(astr(closest(vec, 3.125, value=True),3,pp=True))
        3.131

        >>> vec = np.arange(100).reshape((10,10))/99.*5.
        >>> print(astr(closest(vec, 3.125, value=True),3,pp=True))
        3.131
        >>> print(closest(vec, 3.125))
        62
        >>> print(np.unravel_index(closest(vec, 3.125), vec.shape))
        (6, 2)
        >>> print(astr(vec[np.unravel_index(closest(vec, 3.125), vec.shape)],3,pp=True))
        3.131


        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT
        License. The JAMS Python package originates from the former UFZ Python library,
        Department of Computational Hydrosystems, Helmholtz Centre for Environmental
        Research - UFZ, Leipzig, Germany.

        Copyright (c) 2012-2015 Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written,  MC, Jan 2012
        Modified, MC, Feb 2013 - ported to Python 3
                  MC, Jul 2015 - 
    """
    out = np.ma.argmin(np.ma.abs(np.ma.array(vec)-num))
    if value:
      return vec.flat[out]
    else:
      return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
