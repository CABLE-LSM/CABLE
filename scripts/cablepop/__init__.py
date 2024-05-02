#!/usr/bin/env python
"""
Purpose
=======

cablepop provides functions to use with the CSIRO Atmosphere Biosphere
Land Exchange model (CABLE), coupled to the Populations-Order-Physiology
model (POP) and the Soil-Litter-Iso model (SLI).

The package uses several functions of the JAMS Python package
https://github.com/mcuntz/jams_python
The JAMS package and cablepop are synchronised irregularily.

:copyright: Copyright 2020 Matthias Cuntz, see AUTHORS.md for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
   closest
   llkdtree
"""
from .closest  import closest
from .llkdtree import llKDTree

from .version import __version__, __author__
