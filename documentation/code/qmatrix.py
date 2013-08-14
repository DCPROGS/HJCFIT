from numpy import sum, abs, all
from dcprogs.likelihood import QMatrix

matrix = [ [-3050,        50,  3000,      0,    0], 
           [2./3., -1502./3.,     0,    500,    0], 
           [   15,         0, -2065,     50, 2000], 
           [    0,     15000,  4000, -19000,    0], 
           [    0,         0,    10,      0,  -10] ]

qmatrix = QMatrix(matrix, 2);
 
print(qmatrix)

assert qmatrix.nopen == 2
assert qmatrix.nshut == 3
assert all(abs(qmatrix.matrix - matrix) < 1e-12)

  
# Compute sum over rows, row by row.
for i, row in enumerate(qmatrix.matrix): print("sum(row[{0}]) = {1}".format(i, sum(row)))

# Compute sum over rows, but let numpy do it.
print("Sum over rows: {0}".format(sum(qmatrix.matrix, axis=1)))
  
# Check that sum over rows returns zero.
assert all(abs(sum(qmatrix.matrix, axis=1)) < 1e-12)
  
# Gets AA block.
print("AA block:\n---------\n{0}".format(qmatrix.aa))
                
# aa is simply wrapper around the matrix block 
qmatrix.aa[0, 0] = 42
assert abs(qmatrix[0, 0] - 42) < 1e-12
