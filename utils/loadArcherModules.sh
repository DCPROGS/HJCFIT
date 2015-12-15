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

# Load Python3
module load anaconda-compute/2.2.0-python3

# Python tests are run with behave
# Behave won't work unless we add its install location to $PATH
pip install --user behave
export PATH=~/.local/bin:$PATH

# Needed to dynamically link libraries
export CRAYPE_LINK_TYPE=dynamic
