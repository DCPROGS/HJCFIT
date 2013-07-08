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

