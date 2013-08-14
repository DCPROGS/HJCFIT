from dcprogs.likelihood import QMatrix, ApproxSurvivor

# Define parameters.
qmatrix = QMatrix([ [-3050,        50,  3000,      0,    0], 
                    [2./3., -1502./3.,     0,    500,    0], 
                    [   15,         0, -2065,     50, 2000], 
                    [    0,     15000,  4000, -19000,    0], 
                    [    0,         0,    10,      0,  -10] ], 2)
tau = 1e-4

survivor = ApproxSurvivor(qmatrix, 1e-4);

print(survivor)

print("AF values\n"  \
      "---------\n")
times = [1e-4, 1e-5, 2.0e-5, 2.5e-5]
af_values = survivor.af(times)
for t, v in zip(times, af_values): print("  * at time t={0}:\n{1}\n".format(t, v))

print("")

print("  * Exponents: {0}".format([root for matrix, root in survivor.af_components]))
