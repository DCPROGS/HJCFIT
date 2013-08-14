from numpy import all, abs, arange
from dcprogs.likelihood import QMatrix, DeterminantEq, MissedEventsG

# Define parameters.
qmatrix = QMatrix([ [-3050,        50,  3000,      0,    0], 
                    [2./3., -1502./3.,     0,    500,    0], 
                    [   15,         0, -2065,     50, 2000], 
                    [    0,     15000,  4000, -19000,    0], 
                    [    0,         0,    10,      0,  -10] ], 2)
tau = 1e-4
 
# Create eG from prior knowledge of roots
determinant_eq = DeterminantEq(qmatrix, tau)
af_roots = [( -3045.285776037674, 1), (-162.92946543451328, 1)]
fa_roots = [(-17090.192769236815, 1), (-2058.0812921673496, 1), (-0.24356535498785126, 1)]
eG_from_roots = MissedEventsG(determinant_eq, af_roots, determinant_eq.transpose(), fa_roots)

# Create eG automaticallye
eG_automatic = MissedEventsG(qmatrix, tau)

# Checks the three initialization are equivalent at tau
assert all(abs(eG_from_roots.af(tau) - eG_automatic.af(tau)) < 1e-8)
assert all(abs(eG_from_roots.fa(tau) - eG_automatic.fa(tau)) < 1e-8)

# Checks the three initialization are equivalent at different times
# The functions can be applied to arrays. 
times = arange(tau, 10*tau, 0.1*tau)
assert eG_from_roots.af(times).shape == (len(times), 2, 3)
assert eG_from_roots.fa(times).shape == (len(times), 3, 2)
assert all(abs(eG_from_roots.af(times) - eG_automatic.af(times)) < 1e-8)
assert all(abs(eG_from_roots.fa(times) - eG_automatic.fa(times)) < 1e-8)
