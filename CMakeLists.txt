CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4 FATAL_ERROR)
PROJECT(DCprogs)

# Enable ExternalProject CMake module
INCLUDE(ExternalProject)
# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/external)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

enable_language(CXX)

option(tests       "Enable testing."       on)


find_package(eigen)

if(tests) 
  find_package(GTest)
  enable_testing()
  add_subdirectory(tests)
endif(tests)
