from mpi4py import MPI
from dcpyps import dataset


class MPIHelper:
    def __init__(self):
        if not MPI.Is_initialized():
            MPI.Init()
        self.comm = MPI.COMM_WORLD
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.size = MPI.COMM_WORLD.Get_size()

    def load_data(self, scnfiles, tres, tcrit, conc, verbose=True):
        self.recs = []
        self.bursts = []
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
    # 
    # def load_mec(self, mecfile, rates)
    # # @property
    # # def
