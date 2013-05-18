CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4 FATAL_ERROR)
PROJECT(DCprogs)

option(tests       "Enable testing."       on)

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)
# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/external)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

enable_language(CXX)
include(${CMAKE_MODULE_PATH}/AddCpp11Flags.cmake)
include(${CMAKE_MODULE_PATH}/CheckCXX11Features.cmake)

if(NOT HAS_CXX11_NOEXCEPT)
  set(noexcept TRUE)
endif()
configure_file (
 "${PROJECT_SOURCE_DIR}/DCProgsConfig.h.in"
 "${PROJECT_BINARY_DIR}/DCProgsConfig.h"
)

find_package(eigen REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR})

if(tests) 
  find_package(GTest)
  enable_testing()
endif(tests)

add_subdirectory(likelihood)
