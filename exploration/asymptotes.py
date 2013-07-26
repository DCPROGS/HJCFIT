########################
#   DCProgs computes missed-events likelihood as described in
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

from numpy import array
from scipy.optimize import brentq
from dcprogs.likelihood import DeterminantEq

matrix = array([[ -3050,        50,  3000,      0,    0 ], 
                [ 2./3., -1502./3.,     0,    500,    0 ],  
                [    15,         0, -2065,     50, 2000 ],  
                [     0,     15000,  4000, -19000,    0 ],  
                [     0,         0,    10,      0,  -10 ] ])
determinant_equation = DeterminantEq(matrix, 2, 1e-4, True)

print brentq(determinant_equation, -20000, -1000, full_output=True)
print brentq(determinant_equation, -1000, 1000, full_output=True)

