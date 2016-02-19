#!/bin/bash

# This scripts loads all the necessary modules and environment configuration
# to build HJCFIT

# Use the right version of CMake
# Note only CMake 2.6 seems to be available in the compute nodes
module unload cmake
module load cmake/3.2.3

# Use GNU compiler suite environment
module unload PrgEnv-cray
module load PrgEnv-gnu/5.2.56

# Swig is needed by the python bindings
module load swig

# From the login nodes, we can use anaconda
# Elsewhere, python-compute is recommended

# Load Anaconda for compute nodes with Python2
# Note there is an anaconda-compute module with Python3 available in the
# login nodes, but that only the Python2 one is available in the MOM nodes
# The module load messes with LO_LIBRARY_PATH we may have to unset this but
# leave it for the moment.
# export OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
# export OLD_LIBRARY_PATH=$LIBRARY_PATH
module load anaconda-compute/2.2.0-python2
# export LD_LIBRARY_PATH=$OLD_LD_LIBRARY_PATH
# export LIBRARY_PATH=$OLD_LIBRARY_PATH
# The anaconda module sets PYTHONHOME Which is wrong and breaks conda envs so remove it
unset PYTHONHOME

# Aparently that anaconda-compute module makes git fail, and that is needed by CookOff, etc, so:
module load git

# Environmental variables needed for the "dcprogs" virtual environment to work
export CONDA_ENVS_PATH=$WORK/.conda/envs

# Needed to dynamically link libraries
export CRAYPE_LINK_TYPE=dynamic
