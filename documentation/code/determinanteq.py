from numpy import abs, all, array
from dcprogs import internal_dtype
from dcprogs.likelihood import DeterminantEq, QMatrix

# Define parameters.
matrix = [ [-3050,        50,  3000,      0,    0], 
           [2./3., -1502./3.,     0,    500,    0], 
           [   15,         0, -2065,     50, 2000], 
           [    0,     15000,  4000, -19000,    0], 
           [    0,         0,    10,      0,  -10] ]
qmatrix = QMatrix(matrix, 2)
tau = 1e-4

det0 = DeterminantEq(qmatrix, 1e-4);
det1 = DeterminantEq(matrix, 2, 1e-4);

x = [0, -1, -1e2]
assert all(abs(det0(x) - det1(x)) < 1e-8)

roots = array([-3045.285776037674, -162.92946543451328])
assert all(abs(det0(roots)) < 1e-6 * abs(roots))

transpose = det0.transpose()
roots = array( [-17090.192769236815, -2058.0812921673496, -0.24356535498785126],
               dtype=internal_dtype )
assert all(abs(transpose(roots)) < 1e-6 * abs(roots))

print("  * H({0}):\n{1}\n".format(0, det0.H(0)))
print("  * H({0}):\n{1}\n".format(-1e2, det0.H(-1e2)))
print("  * d[sI-H(s)]/ds for s={0}:\n{1}\n".format(0, det0.s_derivative(0)))
print("  * d[sI-H(s)]/ds for s={0}:\n{1}\n".format(-1e2, det0.s_derivative(-1e2)))
