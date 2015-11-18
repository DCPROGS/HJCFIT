########################
#   DCProgs computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4 FATAL_ERROR)
PROJECT(DCProgs)

option(tests          "Enable testing."                         on)
option(pythonBindings "Enable python bindings."                 on)
option(compileDocs    "Compile c++11 documentation examples."   on)

# Provides backwards compatibility and prevents MACOS_RPATH warning
cmake_policy(VERSION 3.0)

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)
# Set default ExternalProject root directory
set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)

set(CMAKE_SCRIPTS ${PROJECT_SOURCE_DIR}/cmake_modules)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

enable_language(CXX)
include(${CMAKE_SCRIPTS}/AddCpp11Flags.cmake)
include(${CMAKE_SCRIPTS}/CheckCXX11Features.cmake)
include(${CMAKE_SCRIPTS}/CheckIsNaNStatus.cmake)

# Downloads and installs GreatCMakeCookOff
# It contains a number of cmake recipes
include(LookUp-GreatCMakeCookOff)
include(PackageLookup)
# Take advantage of GreatCMakeCookOff to install Eigen
lookup_package(Eigen)

#find_package(eigen REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${PROJECT_BINARY_DIR})

if(tests) 
  enable_testing()
  include(AddGTest)
endif(tests)
if(pythonBindings) 
  include(${CMAKE_SCRIPTS}/AllPythonBindings.cmake)
endif(pythonBindings)

# if(NOT HAS_CXX11_NOEXCEPT)
#   set(noexcept TRUE)
# endif()
# if(NOT DCPROGS_LONG_DOUBLE) 
#   set(DCPROGS_LONG_DOUBLE False CACHE BOOL
#       "If True, will use long doubles rather than simple doubles.")
# endif(NOT DCPROGS_LONG_DOUBLE)

# configure_file (
#  "${PROJECT_SOURCE_DIR}/DCProgsConfig.h.in"
#  "${PROJECT_BINARY_DIR}/DCProgsConfig.h"
# )
# install(FILES ${PROJECT_BINARY_DIR}/DCProgsConfig.h DESTINATION include/dcprogs)

# include(${CMAKE_MODULE_PATH}/documentation.cmake)

# add_subdirectory(likelihood)
# add_subdirectory(data)
# add_subdirectory(python)
# add_subdirectory(documentation)

