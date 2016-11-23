###########################################
Running likelihood calculations in parallel
###########################################




OpenMP
======

HJCFIT will by default be compiled with openmp support. The parallelisation
is by default over either the number of bursts or over the individual open
close transitions within the burst. Typically experiments either have many
short bursts or a few long bursts so it really only makes sense to parallelise
over one of these axes. The code takes care of detecting which axis to
parallelise over automatically. The number of threads can be controlled by
the usual environmental variable ``OMP_NUM_THREADS``. Running on a PC this
will probably be set to the number of cores in the computer which is probably
the optimal solution in most cases.

CMake takes care of identifying the correct compiler flags and enables OpenMP
automatically on all supported platforms. Currently (2016) Clang on OSX does not
support OpenMP (but the code can be compiled on OSX using gcc from homebrew
or similar)

OpenMP can be disabled explicitly by setting the CMake variable ``openmp`` to
off.



MPI
===

The MPI parallelisation over experiments (antagonist concentrations) is
implemented in the Python layer. This means that in an example such as
``exploration/fitGly$4.py`` the MPI code would be implemented directly in the
example complicating the individual examples. In order to simplify this a
wrapper python class has been implemented in ``mpihelpers.MPILikelihoodSolver``.
An example of using this for the same purpose can be seen in
``exploration/mpi/fitGlyR4_mpi.py``. To launch this example you should run
something like ``mpiexec np -4 python fitGlyR4_mpi.py``. This runes 4 MPI
processes (matching the 4 experiments in the fitGlyR4 example) Each MPI
processes may in addition use OpenMP as detailed above to parallelize the
computations of the likelihood for the individual simulations. I.e on a 24 core
Archer node you would most likely want to use 4 MPI processes who in turn  runs
6 OpenMP threads each. The syntax for running a MPI job will depend on the
specific cluster that you are running on so it's recommended to check the
clusters documentation to see how to launch a job. 
