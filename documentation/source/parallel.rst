###########################################
Running likelihood calculations in parallel
###########################################




OpenMP
======

HJCFIT will by default be complied with openmp support. The parallelisation
is by default over either the number of bursts or over the individual open
close transitions within the burst. Typically experiments either have many
short bursts or a few long bursts so it really only makes sense to parallelise
over one of these axes. The code takes care of detecting which axis to
parallelise over automatically. The number of threads can be controlled by
the usual environmental variable ```OMP_NUM_THREADS``. Running on a PC this
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

Write something about MPI
