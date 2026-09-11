########################
#   HJCFIT computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

""" Likelihood sub-package. """
from .likelihood import *
from ._methods import *
from . import random
from . import optimization

# The fitting layer. It imports nothing beyond numpy, scipy and the likelihood
# itself -- it is duck-typed on a mechanism rather than importing one -- so
# exposing it here costs the bare install nothing. What it needs from
# hjcfit[fitting] is a mechanism to be handed, not a module to import.
from . import fitting
