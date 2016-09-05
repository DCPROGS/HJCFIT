from mpi4py import MPI
import sys
from dcpyps import dataset
from dcpyps import dcio
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood
import math
import numpy as np
from scipy.optimize import minimize
import time


class MPILikelihoodSolver:
    def __init__(self, print_freq=100):
        if not MPI.Is_initialized():
            MPI.Init()
        self.comm = MPI.COMM_WORLD
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.size = MPI.COMM_WORLD.Get_size()
        self.mpi_status = np.array(1, 'int')
        self.iternum = 0
        self.print_freq = print_freq

    def load_data(self, scnfiles, tres, tcrit, conc, verbose=True):
        self.recs = []
        self.bursts = []
        self.tcrit = tcrit
        self.tres = tres
        self.conc = conc
        for i in range(len(scnfiles)):
            rec = dataset.SCRecord(scnfiles[i], conc[i], tres[i], tcrit[i])
            rec.record_type = 'recorded'
            self.recs.append(rec)
            self.bursts.append(rec.bursts.intervals())
            if self.rank == 0 and verbose:
                rec.printout()
            if self.size != len(conc):
                outputstring = ("Number of MPI processes much match number of"
                                "concentrations. Got {} MPI processes "
                                "and {} concentrations.".format(self.size,
                                                                len(conc)))
                raise RuntimeError(outputstring)

    def load_mec(self, mecfn, rates):
        version, meclist, max_mecnum = dcio.mec_get_list(mecfn)
        self.mec = dcio.mec_load(mecfn, meclist[2][0])
        self.mec.set_rateconstants(rates)

    def mec_printout(self, target=sys.stdout):
        if self.rank == 0:
            self.mec.printout(target)

    def set_likelihood_func(self, likelihood_kwargs=None):
        self.likelihood = []
        assert len(self.tres) == len(self.tcrit) == len(self.conc)
        if likelihood_kwargs is None:
            likelihood_kwargs = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12,
                                 'itermax': 100, 'lower_bound': -1e6,
                                 'upper_bound': 0}
        for i in range(len(self.tres)):
            self.likelihood.append(Log10Likelihood(self.bursts[i],
                                                   self.mec.kA,
                                                   self.tres[i],
                                                   self.tcrit[i],
                                                   **likelihood_kwargs))

    def mpi_master_likelihood(self, x, args=None):
        self.comm.Bcast([self.mpi_status, MPI.INT], root=0)
        self.comm.Bcast([x, MPI.DOUBLE], root=0)
        self.mec.theta_unsqueeze(np.exp(x))
        lik = np.array(0.0, 'd')
        like = np.array(0.0, 'd')
        self.mec.set_eff('c', self.conc[self.rank])
        lik += -self.likelihood[self.rank](self.mec.Q) * math.log(10)
        self.comm.Reduce([lik, MPI.DOUBLE],
                         [like, MPI.DOUBLE],
                         op=MPI.SUM,
                         root=0)
        return like

    def mpi_slave_likelihood(self):
        self.comm.Bcast([self.mpi_status, MPI.INT], root=0)
        if not self.mpi_status:
            return
        x = np.empty(14, dtype='d')
        self.comm.Bcast([x, MPI.DOUBLE], root=0)
        self.mec.theta_unsqueeze(np.exp(x))
        self.mec.set_eff('c', self.conc[self.rank])
        lik = np.array(0.0, 'd')
        lik += -self.likelihood[self.rank](self.mec.Q) * math.log(10)
        self.comm.Reduce([lik, MPI.DOUBLE], None, op=MPI.SUM, root=0)
        return

    def complete_likelihood(self, x, args=None):
        self.mec.theta_unsqueeze(np.exp(x))
        lik = 0
        for i in range(len(self.conc)):
            self.mec.set_eff('c', self.conc[i])
            lik += -self.likelihood[i](self.mec.Q) * math.log(10)
        return lik

    def print_likelihood_status(self, theta):
        self.iternum += 1
        if self.iternum % self.print_freq == 0:
            lik = self.complete_likelihood(theta)
            print("iteration # {0:d}; log-lik = {1:.6f}".format(self.iternum,
                                                                -lik))
            print(np.exp(theta))

    def run_optimizer(self, options=None, method='Nelder-Mead'):
        if self.rank == 0:
            theta = np.log(self.mec.theta())
            lik = self.complete_likelihood(theta)
            print("\nStarting likelihood (DCprogs)= "
                  "{0:.6f}".format(-lik))
            start = time.clock()
            wallclock_start = time.time()
            self.result = None
            if options is None:
                options = {'xtol': 1e-4, 'ftol': 1e-4, 'maxiter': 5000,
                           'maxfev': 10000, 'disp': True}
            self.result = minimize(self.mpi_master_likelihood,
                                   theta,
                                   method=method,
                                   callback=self.print_likelihood_status,
                                   options=options)
            # Signal slaves to stop
            self.mpi_status = np.array(0, 'int')
            self.comm.Bcast([self.mpi_status, MPI.INT], root=0)
        else:
            while self.mpi_status:
                self.mpi_slave_likelihood()
        if self.rank == 0:
            end = time.clock()
            wallclock_end = time.time()
            print("\nDCPROGS Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
                  %time.localtime()[0:6])
            self.cpu_time = end - start
            self.wallclock_time = wallclock_end - wallclock_start
            print('CPU time in simplex=', self.cpu_time)
            print('Wallclock time in simplex=', self.wallclock_time)
            print('\n\nresult=')
            print(self.result)

            print('\n Final log-likelihood = {0:.6f}'.format(-self.result.fun))
            print('\n Number of iterations = {0:d}'.format(self.result.nit))
            print('\n Number of evaluations = {0:d}'.format(self.result.nfev))
            self.mec.theta_unsqueeze(np.exp(self.result.x))
            print("\n Final rate constants:")
            self.mec.printout(sys.stdout)
            print('\n\n')
