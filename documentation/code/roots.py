from numpy import all
from dcprogs.likelihood import eig
from dcprogs.likelihood import find_upper_bound_for_roots, find_lower_bound_for_roots,        \
                               find_root_intervals, brentq, find_roots, QMatrix, DeterminantEq

qmatrix = QMatrix([ [-3050,        50,  3000,      0,    0], 
                    [2./3., -1502./3.,     0,    500,    0], 
                    [   15,         0, -2065,     50, 2000], 
                    [    0,     15000,  4000, -19000,    0], 
                    [    0,         0,    10,      0,  -10] ], 2)
det = DeterminantEq(qmatrix, 1e-4);


upper_bound = find_upper_bound_for_roots(det);
lower_bound = find_lower_bound_for_roots(det);

get_eigenvalues = lambda s: eig(det.H(s))[0].T
assert all(get_eigenvalues(lower_bound) > lower_bound) 
assert all(get_eigenvalues(upper_bound) < upper_bound) 


print("Root Determination\n"                                         \
      "==================\n\n"                                       \
      "  * Interval containing roots: {lower}, {upper}\n"            \
      "  * Eigenvalues of H at lower bound: {eiglow}\n"              \
      "  * Eigenvalues of H at upper bound: {eigup}\n"               \
      .format( lower=lower_bound, upper=upper_bound,
               eiglow = get_eigenvalues(lower_bound),
               eigup  = get_eigenvalues(upper_bound) ))

intervals = find_root_intervals(det, lower_bound, upper_bound)
for (start, end), multiplicity in intervals:
  root, iterations, function_calls = brentq(det, start, end)
  print("  * Root interval: [{0}, {1}]\n"
        "    Corresponding root: {2}\n".format(start, end, root))


roots = find_roots(det);
print("  * All roots: {0}\n".format([root for root, multiplicity in roots]))
