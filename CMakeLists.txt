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
