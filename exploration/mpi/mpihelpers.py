from mpi4py import MPI
import sys
from dcpyps import dataset
from dcpyps import dcio
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood
import math
import numpy as np


class MPIHelper:
    def __init__(self):
        if not MPI.Is_initialized():
            MPI.Init()
        self.comm = MPI.COMM_WORLD
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.size = MPI.COMM_WORLD.Get_size()
        self.mpi_status = np.array(0, 'int')

    def load_data(self, scnfiles, tres, tcrit, conc, verbose=True):
        self.recs = []
        self.bursts = []
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

    def set_likelihood_func(self, kwargs):
        self.likelihood = []

        for i in range(len(self.recs)):
            self.likelihood.append(Log10Likelihood(self.bursts[i],
                                                   self.mec.kA,
                                                   self.recs[i].tres,
                                                   self.recs[i].tcrit,
                                                   **kwargs))

    def mpi_master_likelihood(self, x, args=None):
        self.comm.Bcast([self.mpi_status, MPI.INT], root=0)
        self.comm.Bcast([x, MPI.DOUBLE], root=0)
        self.mec.theta_unsqueeze(np.exp(x))
        lik = np.array(0.0, 'd')
        like = np.array(0.0, 'd')
        self.mec.set_eff('c', self.conc[self.rank])
        lik += -self.likelihood[self.rank](self.mec.Q) * math.log(10)
        self.comm.Reduce([lik, MPI.DOUBLE], [like, MPI.DOUBLE], op=MPI.SUM, root=0)
        return like
