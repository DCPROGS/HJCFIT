import time
import math
import os
import sys
import pickle
import numpy as np
from scipy.optimize import minimize

from dcpyps.samples import samples
from dcpyps import dataset
from dcpyps import mechanism
from dcpyps.sccalc import scsim
from dcprogs.likelihood import Log10Likelihood
from mpi4py import MPI
from dcprogs.mpihelpers import MPILikelihoodSolver


def simulate_bursts(conc, mec, tr, inst, nmax):
    bursts = []
    for c in conc:
        mec.set_eff('c', c)  # Get Q-matrix for given concentration
        tints, ampls, flags, ntrans = scsim.simulate_intervals(mec,
                                                               tr,
                                                               inst,
                                                               nintmax=nmax)
        rec = dataset.SCRecord(conc=c, tres=tr,
                               itint=tints, iampl=ampls, iprops=flags)
        bursts.append(rec.bursts.intervals())
    return bursts


def constrain(mec):

    for i in range(len(mec.Rates)):
        mec.Rates[i].fixed = False
    # Constrained rates.
    mec.Rates[21].is_constrained = True
    mec.Rates[21].constrain_func = mechanism.constrain_rate_multiple
    mec.Rates[21].constrain_args = [20, 1.5]
    mec.Rates[18].is_constrained = True
    mec.Rates[18].constrain_func = mechanism.constrain_rate_multiple
    mec.Rates[18].constrain_args = [19, 2]
    mec.Rates[14].is_constrained = True
    mec.Rates[14].constrain_func = mechanism.constrain_rate_multiple
    mec.Rates[14].constrain_args = [12, 3]
    mec.Rates[13].is_constrained = True
    mec.Rates[13].constrain_func = mechanism.constrain_rate_multiple
    mec.Rates[13].constrain_args = [12, 2]
    mec.Rates[15].is_constrained = True
    mec.Rates[15].constrain_func = mechanism.constrain_rate_multiple
    mec.Rates[15].constrain_args = [17, 3]
    mec.Rates[16].is_constrained = True
    mec.Rates[16].constrain_func = mechanism.constrain_rate_multiple
    mec.Rates[16].constrain_args = [17, 2]
    mec.set_mr(True, 9, 0)
    mec.set_mr(True, 11, 1)
    mec.update_constrains()
    return mec

mysolver = MPILikelihoodSolver()

# LOAD Burzomato 2004 mechanism (GlyR, glycine, WT)
mec_true = samples.GlyR_flip()
ig = [4200, 28000, 130000, 3400, 2100, 6700, 180, 6800, 22000,
      29266, 18000, 948, 302, 604, 906, 1.77e6, 1.18e6, 0.59e6, 300e6, 150e6,
      2500, 3750]
mec_true.set_rateconstants(ig)
mec_true = constrain(mec_true)
if mysolver.rank == 0:
    theta_true = np.log(mec_true.theta())
mec = samples.GlyR_flip()
tr = 0.000030
tres = [0.000030, 0.000030, 0.000030, 0.000030]
tcrit = [0.004, -1, -0.06, -0.02]
conc = [10e-6, 30e-6, 100e-6, 1000e-6]
nintmax = 15000  # Number of intervals to be simulated
inst = mec.k - 1  # initial state
ig1 = [500, 2000, 15000, 5000, 2700, 800, 1200, 12000, 120000,
       7500, 1500, 300, 1000, 2000, 3000, 13.5e6, 9e6, 4.5e6, 900e6, 450e6,
       4000, 6000]
ig2 = [8000, 20000, 160000, 5000, 2700, 4000, 100, 4000, 30000,
       6250, 10000, 3000, 500, 1000, 1500, 3e6, 2e6, 1e6, 200e6, 100e6,
       4000, 6000]
kwargs = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12, 'itermax': 100,
          'lower_bound': -1e6, 'upper_bound': 0}
simplex_options = {'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': 5000,
                   'maxfev': 10000, 'disp': True}

if mysolver.rank == 0:
    start = time.clock()
    burstdata = simulate_bursts(conc, mec_true, tr, inst, nintmax)
    data = {'bursts': burstdata}
    end = time.clock()
    print('CPU time in simulation=', end - start)
else:
    data = None

data = mysolver.comm.bcast(data, root=0)
bursts = data['bursts']
# Set initial guesses
mec.set_rateconstants(ig2)
mec = constrain(mec)

mysolver = MPILikelihoodSolver()
mysolver.mec = mec
mysolver.bursts = bursts
mysolver.conc = conc
mysolver.tres = tres
mysolver.tcrit = tcrit

mysolver.set_likelihood_func(kwargs)

if mysolver.rank == 0:
    likelihood_true = mysolver.complete_likelihood(theta_true)
    mysolver.mec.set_rateconstants(ig2)
    print("Generated likelihood {}".format(-likelihood_true))

mysolver.run_optimizer()

if mysolver.rank == 0:
    directory = 'results'
    if not os.path.exists(directory):
        os.makedirs(directory)
    fname = directory + '/' + time.strftime("%Y%m%d-%H%M%S") + '.result'
    f = open(fname, 'wb')
    pickle.dump([mysolver.result,
                 mysolver.cpu_time,
                 mysolver.wallclock_time,
                 likelihood_true],
                f,
                pickle.HIGHEST_PROTOCOL)
