import time
import math
import sys
import numpy as np
from scipy.optimize import minimize

from dcpyps import dcio
from dcpyps import dataset
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood
from mpi4py import MPI
from mpihelpers import MPIHelper

# LOAD DATA: Burzomato 2004 example set.
scnfiles = [["../../../DCPYPS/dcpyps/samples/glydemo/A-10.scn"],
            ["../../../DCPYPS/dcpyps/samples/glydemo/B-30.scn"],
            ["../../../DCPYPS/dcpyps/samples/glydemo/C-100.scn"],
            ["../../../DCPYPS/dcpyps/samples/glydemo/D-1000.scn"]]
tres = [0.000030, 0.000030, 0.000030, 0.000030]
tcrit = [0.004, -1, -0.06, -0.02]
conc = [10e-6, 30e-6, 100e-6, 1000e-6]

mympi = MPIHelper()
mympi.load_data(scnfiles, tres, tcrit, conc)

recs = mympi.recs
bursts = mympi.bursts

# LOAD FLIP MECHANISM USED Burzomato et al 2004
mecfn = "../../../DCPYPS/dcpyps/samples/mec/demomec.mec"
# PREPARE RATE CONSTANTS.
rates = [4500.0, 700.0, 2500.0, 1800.0, 900.0, 18000.0, 200.0,
         0.1100E+06, 4900.0, 0.4000E+09, 1850.0, 10000.0, 5000.0,
         0.7500E+09, 8500.0, 1050.0, 3500.0, 0.5000E+07, 2300.0,
         0.9500E+07, 1950, 0.130000E+08]
mympi.load_mec(mecfn, rates)

for i in range(len(mympi.mec.Rates)):
    mympi.mec.Rates[i].fixed = False

# Constrained rates.
mympi.mec.Rates[21].is_constrained = True
mympi.mec.Rates[21].constrain_func = mechanism.constrain_rate_multiple
mympi.mec.Rates[21].constrain_args = [17, 3]
mympi.mec.Rates[19].is_constrained = True
mympi.mec.Rates[19].constrain_func = mechanism.constrain_rate_multiple
mympi.mec.Rates[19].constrain_args = [17, 2]
mympi.mec.Rates[16].is_constrained = True
mympi.mec.Rates[16].constrain_func = mechanism.constrain_rate_multiple
mympi.mec.Rates[16].constrain_args = [20, 3]
mympi.mec.Rates[18].is_constrained = True
mympi.mec.Rates[18].constrain_func = mechanism.constrain_rate_multiple
mympi.mec.Rates[18].constrain_args = [20, 2]
mympi.mec.Rates[8].is_constrained = True
mympi.mec.Rates[8].constrain_func = mechanism.constrain_rate_multiple
mympi.mec.Rates[8].constrain_args = [12, 1.5]
mympi.mec.Rates[13].is_constrained = True
mympi.mec.Rates[13].constrain_func = mechanism.constrain_rate_multiple
mympi.mec.Rates[13].constrain_args = [9, 2]
mympi.mec.update_constrains()

mympi.mec.set_mr(True, 7, 0)
mympi.mec.set_mr(True, 15, 1)

mympi.mec_printout()

theta = np.log(mympi.mec.theta())
likelihood_kwargs = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12, 'itermax': 100,
          'lower_bound': -1e6, 'upper_bound': 0}

mympi.set_likelihood_func(likelihood_kwargs)


def dcprogslik(x, args=None):
    mympi.mec.theta_unsqueeze(np.exp(x))
    lik = 0
    for i in range(len(conc)):
        mympi.mec.set_eff('c', conc[i])
        lik += -mympi.likelihood[i](mympi.mec.Q) * math.log(10)
    return lik


# def mpidcprogslik(x, args=None):
#     mympi.comm.Bcast([mpi_status, MPI.INT], root=0)
#     mympi.comm.Bcast([x, MPI.DOUBLE], root=0)
#     mympi.mec.theta_unsqueeze(np.exp(x))
#     lik = np.array(0.0, 'd')
#     like = np.array(0.0, 'd')
#     mympi.mec.set_eff('c', conc[mympi.rank])
#     lik += -mympi.likelihood[mympi.rank](mympi.mec.Q) * math.log(10)
#     mympi.comm.Reduce([lik, MPI.DOUBLE], [like, MPI.DOUBLE], op=MPI.SUM, root=0)
#     return like


def mpislavedcprogslik():
    mympi.comm.Bcast([mympi.mpi_status, MPI.INT], root=0)
    if not mympi.mpi_status:
        return mympi.mpi_status
    x = np.empty(14, dtype='d')
    mympi.comm.Bcast([x, MPI.DOUBLE], root=0)
    mympi.mec.theta_unsqueeze(np.exp(x))
    mympi.mec.set_eff('c', conc[mympi.rank])
    lik = np.array(0.0, 'd')
    lik += -mympi.likelihood[mympi.rank](mympi.mec.Q) * math.log(10)
    mympi.comm.Reduce([lik, MPI.DOUBLE], None, op=MPI.SUM, root=0)
    return mympi.mpi_status


def printiter(theta):
    global iternum
    iternum += 1
    if iternum % 100 == 0:
        lik = dcprogslik(theta)
        print("iteration # {0:d}; log-lik = {1:.6f}".format(iternum, -lik))
        print(np.exp(theta))


iternum = 0
lik = dcprogslik(theta)
mympi.mpi_status = np.array(1, 'int')
if mympi.rank == 0:
    print("\nStarting likelihood (DCprogs)= {0:.6f} on {1}".format(-lik, mympi.rank))
    start = time.clock()
    wallclock_start = time.time()
    success = False
    result = None
    options = {'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': 5000,
               'maxfev': 10000, 'disp': True}
    result = minimize(mympi.mpi_master_likelihood, theta, method='Nelder-Mead',
                      callback=printiter, options=options)
    # Signal slaves to stop
    mympi.mpi_status = np.array(0, 'int')
    mympi.comm.Bcast([mympi.mpi_status, MPI.INT], root=0)
else:
    while mympi.mpi_status:
        mympi.mpi_status = mpislavedcprogslik()

if mympi.rank == 0:
    end = time.clock()
    wallclock_end = time.time()
    print("\nDCPROGS Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
          %time.localtime()[0:6])
    print('CPU time in simplex=', end - start)
    print('Wallclock time in simplex=', wallclock_end - wallclock_start)
    print('\n\nresult=')
    print(result)

    print('\n Final log-likelihood = {0:.6f}'.format(-result.fun))
    print('\n Number of iterations = {0:d}'.format(result.nit))
    print('\n Number of evaluations = {0:d}'.format(result.nfev))
    mympi.mec.theta_unsqueeze(np.exp(result.x))
    print("\n Final rate constants:")
    mympi.mec.printout(sys.stdout)
    print('\n\n')
