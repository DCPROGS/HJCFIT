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


find_package(eigen REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR})

if(tests) 
  find_package(GTest)
  enable_testing()
endif(tests)
if(pythonBindings) 
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  find_package(PythonLibs REQUIRED)
  find_package(numpy REQUIRED) 
  if(NOT PYINSTALL_DIR)
    execute_process( 
      COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
      OUTPUT_VARIABLE PYINSTALL_DIR
    )
    if(PYINSTALL_DIR )
      string (STRIP ${PYINSTALL_DIR} PYINSTALL_DIR)
      set(PYINSTALL_DIR ${PYINSTALL_DIR} CACHE STRING "Version of the Python interpreter.")
    else()
      set( PYINSTALL_DIR lib/python${PYTHON_VERSION}/site-packages/ 
           CACHE PATH "Python modules will be installed here." )
    endif(PYINSTALL_DIR)
    mark_as_advanced(PYINSTALL_DIR)
    MESSAGE(STATUS "[Python] installation directory: ${PYINSTALL_DIR}")
  endif(NOT PYINSTALL_DIR)

  set(DCPROGS_PYTHON_BINDINGS True)
endif(pythonBindings)

if(NOT HAS_CXX11_NOEXCEPT)
  set(noexcept TRUE)
endif()
configure_file (
 "${PROJECT_SOURCE_DIR}/DCProgsConfig.h.in"
 "${PROJECT_BINARY_DIR}/DCProgsConfig.h"
)



add_subdirectory(likelihood)
