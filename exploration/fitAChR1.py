#! /usr/bin/python
"""
Maximum likelihood fit demo.
"""

import sys
import time
import math
import numpy as np
from scipy.optimize import minimize

from dcpyps.samples import samples
from dcpyps import dataset
from dcpyps import mechanism

from dcprogs.likelihood import Log10Likelihood

def dcprogslik(x):
    mec.theta_unsqueeze(np.exp(x))
    mec.set_eff('c', conc)
    return -likelihood(mec.Q) * logfac

def printiter(theta):
    global iternum
    iternum += 1
    if iternum % 100 == 0:
        lik = dcprogslik(theta)
        print("iteration # {0:d}; log-lik = {1:.6f}".format(iternum, -lik))
        print(np.exp(theta))

print('\nTesting CH82 fit to single channel data:')

tres = 0.0001
tcrit = 0.004
conc = 100e-9
# LOAD DATA.
filename = "../../DCPYPS/dcpyps/samples/scn/CH82.scn"
rec1 = dataset.SCRecord([filename], conc, tres, tcrit)
rec1.record_type = 'recorded'
rec1.printout()

# LOAD DEMO MECHANISM (C&H82 numerical example).
mec = samples.CH82()
# PREPARE RATE CONSTANTS.
# Fixed rates.
mec.Rates[7].fixed = True
# Constrained rates.
mec.Rates[5].is_constrained = True
mec.Rates[5].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[5].constrain_args = [4, 2]
mec.Rates[6].is_constrained = True
mec.Rates[6].constrain_func = mechanism.constrain_rate_multiple
mec.Rates[6].constrain_args = [8, 2]
mec.update_constrains()

# Initial guesses. 
rates = mec.unit_rates()
rates = [100, 3000, 10000, 100, 1000, 1000, 1e+7, 5e+7, 6e+7, 10]
#    rates = [6.5, 14800, 3640, 362, 1220, 2440, 1e+7, 5e+8, 2.5e+8, 55]
mec.set_rateconstants(rates)
mec.update_mr()
mec.printout(sys.stdout)
theta = mec.theta()
print ('\ntheta=', theta)

#######   DCPROGS likelihood
bursts = rec1.bursts.intervals()
logfac = math.log(10)
likelihood = Log10Likelihood(bursts, mec.kA, tres, tcrit)

iternum = 0
start = time.clock()
start_wall = time.time()
res = minimize(dcprogslik, np.log(theta), method='Nelder-Mead', callback=printiter,)
t3 = time.clock() - start
t3_wall = time.time() - start_wall
print ("\n\n\nScyPy.minimize (Nelder-Mead) Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
    %time.localtime()[0:6])
print ('CPU time in ScyPy.minimize (Nelder-Mead)=', t3)
print ('Wall clock time in ScyPy.minimize (Nelder-Mead)=', t3_wall)
print ('xout', res.x)
mec.theta_unsqueeze(np.exp(res.x))
print ("\n Final rate constants:")
mec.printout(sys.stdout)
lik = dcprogslik(res.x)
print ("\nFinal likelihood = {0:.6f}".format(-lik))
