from numpy import all, abs, NaN
from dcprogs.likelihood import Log10Likelihood

bursts = [  [0.1, 0.2, 0.1],                  # 1st burst 
            [0.2],                            # 2nd burst
            [0.15, 0.16, 0.18, 0.05, 0.1] ]   # 3rd burst
""" List of bursts.  

    Each burst is a list of observed open and shut intervals. 
    There should always be an odd number of intervals, since bursts end in a shut states.
"""

likelihood = Log10Likelihood(bursts, nopen=2, tau=0.01, tcritical=None)

print(likelihood)

matrix = [[ -3050,        50,  3000,      0,    0 ], 
          [ 2./3., -1502./3.,     0,    500,    0 ],  
          [    15,         0, -2065,     50, 2000 ],  
          [     0,     15000,  4000, -19000,    0 ],  
          [     0,         0,    10,      0,  -10 ] ]

result = likelihood(matrix)

print("Computation: {0}".format(result))

# Get the number of bursts
assert len(likelihood) == 3
# Check the second burst is what we expect
assert all(abs(likelihood[2] - [0.15, 0.16, 0.18, 0.05, 0.1]) < 1e-12)
# Modify, and check the second burst
likelihood[2] = [0.15, 0.16, 0.18]
assert all(abs(likelihood[2] - [0.15, 0.16, 0.18]) < 1e-12)
# Add an extra burst, and check that it is there
likelihood.append([0.25, 0.013, 0.013])
assert len(likelihood) == 4
assert all(abs(likelihood[3] - [0.25, 0.013, 0.013]) < 1e-12)

# some attributes act as switches and get converted back to None when given special values
assert likelihood.tcritical is None
likelihood.tcritical = NaN
assert likelihood.tcritical is None
likelihood.tcritical = -1
assert likelihood.tcritical is None
likelihood.tcritical = 0.5
assert likelihood.tcritical is not None and abs(likelihood.tcritical - 0.5) < 1e-12

assert likelihood.lower_bound is None
likelihood.lower_bound = NaN
assert likelihood.lower_bound is None
likelihood.lower_bound = -1e6
assert likelihood.lower_bound is not None and abs(likelihood.lower_bound + 1e6) < 1e-12
likelihood.lower_bound = None
assert likelihood.lower_bound is None
