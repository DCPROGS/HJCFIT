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

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)

# Downloads and installs GreatCMakeCookOff
# It contains a number of cmake recipes
include(LookUp-GreatCMakeCookOff)
include(CMakeDependentOption)
# Hunter helps us download stuff automatically
include(AddHunter)
set(HUNTER_ROOT ${CMAKE_BINARY_DIR}/hunter_root)
initialize_hunter(HUNTER_ROOT ${HUNTER_ROOT} LOCAL)

PROJECT(HJCFIT)

option(tests          "Enable testing."                         on)
option(pythonBindings "Enable python bindings."                 on)
option(compileDocs    "Compile c++11 documentation examples."   on)
option(openmp         "Enable OpenMP parallelization."          on)
option(executenotebooks "Execute example Jupyter notebooks."    on)
CMAKE_DEPENDENT_OPTION(DCPROGS_USE_MPFR "Enable fallback to MPFR Multi precision" ON "UNIX" OFF)

set(JUPYTER_KERNEL "" CACHE STRING "Jupyter kernel to run notebooks with defaults to python2 or python3 depending on python version")
# Provides backwards compatibility and prevents MACOS_RPATH warning
if((CMAKE_MAJOR_VERSION GREATER 3) OR (CMAKE_MAJOR_VERSION EQUAL 3))
      set(CMAKE_MACOSX_RPATH 1)
endif()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# Check if our link directory is already an implicit link directory
# Only add it to RPath if not.
if(NOT WIN32)
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif()
endif()
# Set a default build type to Release if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)
# Set default ExternalProject root directory
set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)
set(CMAKE_SCRIPTS ${PROJECT_SOURCE_DIR}/cmake_modules)

include(PackageLookup)

if(pythonBindings)
  include(${CMAKE_SCRIPTS}/AllPythonBindings.cmake)
  # Unless explicitly initialized from command line or installation occurs in a
  # virtual environment, set install prefix to CookOff-set path to python
  # interpreter.
  if (EXISTS $ENV{VIRTUAL_ENV})
    set(CMAKE_INSTALL_PREFIX "$ENV{VIRTUAL_ENV}" CACHE string "" FORCE)
  elseif (EXISTS $ENV{CONDA_ENV_PATH})
    set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_ENV_PATH}" CACHE string "" FORCE)
  elseif (WIN32 AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    FILE(TO_CMAKE_PATH ${PYTHON_INTERP_PREFIX} PYTHON_INTERP_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${PYTHON_INTERP_PREFIX}" CACHE string "" FORCE)
  endif()
endif(pythonBindings)

# Install Eigen the CookOff way
# Need at least eigen 3.2 to support mpfr
lookup_package(Eigen3 REQUIRED 3.2.0)

# Adds flags to compiler to launch it into c++11 land
include(AddCPP11Flags)

# The following will print out all available features.
include(CheckCXX11Features)

# Check features (before REQUIRED they are optional)
cxx11_feature_check(constexpr unique_ptr random_device constructor_delegate
                    trivial_type_traits noexcept
                    REQUIRED auto lambda static_assert rvalue_references decltype
                    enable_if type_traits shared_ptr initialization)
include(CheckIsNaN)
if(openmp)
  FIND_PACKAGE(OpenMP)
endif()

if(OPENMP_FOUND)
  if(NOT (${CMAKE_C_FLAGS} MATCHES ${OpenMP_C_FLAGS}))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  endif()
  if(NOT (${CMAKE_CXX_FLAGS} MATCHES ${OpenMP_CXX_FLAGS}))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  endif()
  if(NOT ("${CMAKE_EXE_LINKER_FLAGS}" MATCHES "${OpenMP_EXE_LINKER_FLAGS}"))
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" CACHE STRING "" FORCE)
  endif()
elseif(openmp)
  MESSAGE(STATUS "Code will be compiled without OpenMP support.")
endif(OPENMP_FOUND)
if ( CMAKE_COMPILER_IS_GNUCC )
  set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
endif ( CMAKE_COMPILER_IS_GNUCC )
if(tests)
  enable_testing()
  include(${CMAKE_SCRIPTS}/CTestCustomPreWorkAround.cmake)
  # Use the CookOff way of including GTest
  include(AddGTest)
endif(tests)

if(NOT HAS_CXX11_NOEXCEPT)
  set(noexcept TRUE)
endif()

if(NOT DCPROGS_LONG_DOUBLE)
  set(DCPROGS_LONG_DOUBLE False CACHE BOOL
      "If True, will use long doubles rather than simple doubles.")
endif(NOT DCPROGS_LONG_DOUBLE)

if(NOT DCPROGS_USE_MPFR)
  set(DCPROGS_USE_MPFR False CACHE BOOL
  "If True, will use MPFR arbitrary precision floats as fall back if regular double calculations fail.")
endif(NOT DCPROGS_USE_MPFR)

if(DCPROGS_USE_MPFR)
  find_or_add_hunter_package(GMP)
  find_or_add_hunter_package(MPFR)
  include_directories(${PROJECT_SOURCE_DIR})
endif(DCPROGS_USE_MPFR)
configure_file (
 "${PROJECT_SOURCE_DIR}/DCProgsConfig.h.in"
 "${PROJECT_BINARY_DIR}/DCProgsConfig.h"
)

# Save all DCPROGS headers in the same place in Windows
install(FILES ${PROJECT_BINARY_DIR}/DCProgsConfig.h DESTINATION include/dcprogs)

include(${CMAKE_SCRIPTS}/documentation.cmake)

add_subdirectory(likelihood)
add_subdirectory(data)
add_subdirectory(python)
add_subdirectory(documentation)
