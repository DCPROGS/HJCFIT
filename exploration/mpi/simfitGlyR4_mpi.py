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

if not MPI.Is_initialized():
    MPI.Init()

# getting basic info
comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()


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


# LOAD Burzomato 2004 mechanism (GlyR, glycine, WT)
mec_true = samples.GlyR_flip()
ig = [4200, 28000, 130000, 3400, 2100, 6700, 180, 6800, 22000,
      29266, 18000, 948, 302, 604, 906, 1.77e6, 1.18e6, 0.59e6, 300e6, 150e6,
      2500, 3750]
mec_true.set_rateconstants(ig)
mrc_true = constrain(mec_true)
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

start = time.clock()
bursts = simulate_bursts(conc, mec_true, tr, inst, nintmax)
end = time.clock()
print('CPU time in simulation=', end - start)
# Set initial guesses
mec.set_rateconstants(ig2)
mec = constrain(mec)

theta = np.log(mec.theta())
likelihood = []
for i in range(len(bursts)):
    likelihood.append(Log10Likelihood(bursts[i], mec.kA,
                                      tres[i], tcrit[i], **kwargs))


def dcprogslik(x, args=None):
    mec.theta_unsqueeze(np.exp(x))
    lik = 0
    for i in range(len(conc)):
        mec.set_eff('c', conc[i])
        lik += -likelihood[i](mec.Q) * math.log(10)
    return lik


def mpidcprogslik(x, args=None):
    comm.Bcast([mpi_status, MPI.INT], root=0)
    comm.Bcast([x, MPI.DOUBLE], root=0)
    mec.theta_unsqueeze(np.exp(x))
    lik = np.array(0.0, 'd')
    like = np.array(0.0, 'd')
    mec.set_eff('c', conc[rank])
    lik += -likelihood[rank](mec.Q) * math.log(10)
    comm.Reduce([lik, MPI.DOUBLE], [like, MPI.DOUBLE], op=MPI.SUM, root=0)
    return like


def mpislavedcprogslik():
    comm.Bcast([mpi_status, MPI.INT], root=0)
    if not mpi_status:
        return mpi_status
    x = np.empty(14, dtype='d')
    comm.Bcast([x, MPI.DOUBLE], root=0)
    mec.theta_unsqueeze(np.exp(x))
    mec.set_eff('c', conc[rank])
    lik = np.array(0.0, 'd')
    lik += -likelihood[rank](mec.Q) * math.log(10)
    comm.Reduce([lik, MPI.DOUBLE], None, op=MPI.SUM, root=0)
    return mpi_status


mpi_status = np.array(1, 'int')
if rank == 0:
    start = time.clock()
    wallclock_start = time.time()
    success = False
    result = None
    result = minimize(mpidcprogslik, theta, method='Nelder-Mead',
                      options=simplex_options)
    # Signal slaves to stop
    mpi_status = np.array(0, 'int')
    comm.Bcast([mpi_status, MPI.INT], root=0)
else:
    while mpi_status:
        mpi_status = mpislavedcprogslik()

if rank == 0:
    end = time.clock()
    wallclock_end = time.time()
    print("\nDCPROGS Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
          %time.localtime()[0:6])
    print('CPU time in simplex=', end - start)
    print('Wallclock time in simplex=', wallclock_end - wallclock_start)

    directory = 'results'
    if not os.path.exists(directory):
        os.makedirs(directory)
    fname = directory + '/' + time.strftime("%Y%m%d-%H%M%S") + '.result'
    f = open(fname, 'wb')
    pickle.dump([result, end - start], f, pickle.HIGHEST_PROTOCOL)
