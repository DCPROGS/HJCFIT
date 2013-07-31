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
PROJECT(DCprogs)

option(tests          "Enable testing."               on)
option(pythonBindings "Enable python bindings."       on)

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)
# Set default ExternalProject root directory
set(EXTERNAL_ROOT ${CMAKE_BINARY_DIR}/external)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

enable_language(CXX)
include(${CMAKE_MODULE_PATH}/AddCpp11Flags.cmake)
include(${CMAKE_MODULE_PATH}/CheckCXX11Features.cmake)
include(${CMAKE_MODULE_PATH}/CheckIsNaNStatus.cmake)


find_package(eigen REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR})

if(tests) 
  find_package(GTest)
  enable_testing()
endif(tests)
if(pythonBindings) 
  include(${CMAKE_MODULE_PATH}/AllPythonBindings.cmake)
endif(pythonBindings)

if(NOT HAS_CXX11_NOEXCEPT)
  set(noexcept TRUE)
endif()
if(NOT DCPROGS_LONG_DOUBLE) 
  set(DCPROGS_LONG_DOUBLE False CACHE BOOL
      "If True, will use long doubles rather than simple doubles.")
endif(NOT DCPROGS_LONG_DOUBLE)

configure_file (
 "${PROJECT_SOURCE_DIR}/DCProgsConfig.h.in"
 "${PROJECT_BINARY_DIR}/DCProgsConfig.h"
)

add_subdirectory(likelihood)
add_subdirectory(data)
add_subdirectory(python)

include(${CMAKE_MODULE_PATH}/documentation.cmake)
