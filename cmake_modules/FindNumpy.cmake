# +-----------------------------------------------------------------------------+
# | $Id:: FindNumPy.cmake 3140 2009-09-29 13:49:36Z baehren                   $ |
# +-----------------------------------------------------------------------------+
# |   Copyright (C) 2007                                                        |
# |   Lars B"ahren (bahren@astron.nl)                                           |
# |                                                                             |
# |   This program is free software; you can redistribute it and/or modify      |
# |   it under the terms of the GNU General Public License as published by      |
# |   the Free Software Foundation; either version 2 of the License, or         |
# |   (at your option) any later version.                                       |
# |                                                                             |
# |   This program is distributed in the hope that it will be useful,           |
# |   but WITHOUT ANY WARRANTY; without even the implied warranty of            |
# |   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             |
# |   GNU General Public License for more details.                              |
# |                                                                             |
# |   You should have received a copy of the GNU General Public License         |
# |   along with this program; if not, write to the                             |
# |   Free Software Foundation, Inc.,                                           |
# |   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 |
# +-----------------------------------------------------------------------------+

# - Check for the presence of NUMPY
#
# The following variables are set when NUMPY is found:
#  NUMPY_FOUND              = Set to true, if all components of NUMPY have been
#                             found.
#  NUMPY_INCLUDES           = Include path for the header files of NUMPY
#  NUMPY_MULTIARRAY_LIBRARY = Path to the multiarray shared library
#  NUMPY_SCALARMATH_LIBRARY = Path to the scalarmath shared library
#  NUMPY_LIBRARIES          = Link these to use NUMPY
#  NUMPY_LFLAGS             = Linker flags (optional)
#  NUMPY_API_VERSION        = API version of the installed and available NumPy
#                             package

## -----------------------------------------------------------------------------
## As the shared libraries of a Python module typically do not contain the 
## usual prefix, we need to remove it while searching for the NumPy libraries.
## In order however not to affect other CMake modules we need to swap back in the
## original prefixes once the end of this module is reached.
if(NOT PYTHONLIBS_FOUND)
  find_package(PythonLibs REQUIRED)
endif(NOT PYTHONLIBS_FOUND)
if(NOT PYTHONINTERP_FOUND)
  find_package(PythonInterp REQUIRED)
endif(NOT PYTHONINTERP_FOUND)
if(NUMPY_FOUND AND NUMPY_LIBRARIES AND NUMPY_INCLUDES)
  set(NUMPY_FIND_QUIETLY TRUE)
endif(NUMPY_FOUND AND NUMPY_LIBRARIES AND NUMPY_INCLUDES)


set (TMP_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES})

set (CMAKE_FIND_LIBRARY_PREFIXES "" CACHE STRING
  "Library prefixes"
  FORCE
  )

## -----------------------------------------------------------------------------
## Let's assume the Python executable is smarter about finding NumPy than we
## are, and try asking it before searching ourselves.
## This is necessary to e.g. pick up the MacPorts NumPy installation, which
## ends up in /opt/local/Library/Frameworks/Python.framework ...

execute_process (
  COMMAND ${PYTHON_EXECUTABLE} -c "import numpy, os; print os.path.dirname(numpy.__file__)"
  OUTPUT_VARIABLE numpy_path
  )
if (numpy_path)
  string (STRIP ${numpy_path} numpy_search_path)
else (numpy_path)
  set (numpy_search_path ${lib_locations})
endif (numpy_path)

## -----------------------------------------------------------------------------
## Check for the header files

find_path (NUMPY_ARRAYOBJECT_H numpy/arrayobject.h
  PATHS
  ${numpy_search_path}
  PATH_SUFFIXES
  python
  core/include
  python/numpy/core/include
  NO_DEFAULT_PATH
)

find_path (NUMPY_NDARRAYOBJECT_H numpy/ndarrayobject.h
  PATHS
  ${numpy_search_path}
  PATH_SUFFIXES
  python
  core/include
  python/numpy/core/include
  NO_DEFAULT_PATH
)

find_path (NUMPY_INCLUDES numpy/__multiarray_api.h numpy/multiarray_api.txt
  PATHS
  ${numpy_search_path}
  PATH_SUFFIXES
  python
  core/include
  python/numpy/core/include
  python${PYTHON_VERSION}
  python${PYTHON_VERSION}/site-packages/numpy
  python${PYTHON_VERSION}/site-packages/numpy/core/include
  NO_DEFAULT_PATH
  )

## -----------------------------------------------------------------------------
## Check for the library

find_library (NUMPY_MULTIARRAY_LIBRARY multiarray
  PATHS
  ${numpy_search_path}
  PATH_SUFFIXES
  python
  core
  python/numpy/core
  python${PYTHON_VERSION}/site-packages/numpy/core
  NO_DEFAULT_PATH
  )
if (NUMPY_MULTIARRAY_LIBRARY)
  get_filename_component (NUMPY_MULTIARRAY_LIBRARY ${NUMPY_MULTIARRAY_LIBRARY} ABSOLUTE)
  list (APPEND NUMPY_LIBRARIES ${NUMPY_MULTIARRAY_LIBRARY})
endif (NUMPY_MULTIARRAY_LIBRARY)

find_library (NUMPY_SCALARMATH_LIBRARY scalarmath
  PATHS
  ${numpy_search_path}
  PATH_SUFFIXES
  python
  core
  python/numpy/core
  python${PYTHON_VERSION}/site-packages/numpy/core
  NO_DEFAULT_PATH
  )
if (NUMPY_SCALARMATH_LIBRARY)
  list (APPEND NUMPY_LIBRARIES ${NUMPY_SCALARMATH_LIBRARY})
endif (NUMPY_SCALARMATH_LIBRARY)

## -----------------------------------------------------------------------------
## Try to determine the API version

## Try to locate numpy/version.py first; if this script is not available, loading
## the module in order to query its version does not make any sense.

find_file (NUMPY_VERSION_PY version.py
  PATHS
  ${numpy_search_path}
  /usr/lib64
  /usr/local/lib64
  PATH_SUFFIXES
  python/numpy
  python${PYTHON_VERSION}/site-packages/numpy
  NO_DEFAULT_PATH
  )

if (NUMPY_VERSION_PY AND PYTHON_EXECUTABLE)
  ## some basic feedback
  if (NOT NUMPY_FIND_QUIETLY)
    message (STATUS "[NumPy] Found version.py - running Python to import module numpy.")
  endif (NOT NUMPY_FIND_QUIETLY)
  ## Run Python to import module numpy and print the version information
  execute_process (
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.__version__"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    RESULT_VARIABLE numpy_version_test_result
    OUTPUT_VARIABLE numpy_version_test_output
    ERROR_VARIABLE numpy_version_test_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else (NUMPY_VERSION_PY AND PYTHON_EXECUTABLE)
  message (STATUS "[NumPy] Unable to process version.py script!")
endif (NUMPY_VERSION_PY AND PYTHON_EXECUTABLE)

if (numpy_version_test_output)

  set (NUMPY_API_VERSION ${numpy_version_test_output})
  if (NOT NUMPY_FIND_QUIETLY)
    MESSAGE(STATUS "[NumPy] version ${NUMPY_API_VERSION}")
  endif (NOT NUMPY_FIND_QUIETLY)
  
else (numpy_version_test_output)
  
  if (NUMPY_NDARRAYOBJECT_H)
    file (STRINGS ${NUMPY_NDARRAYOBJECT_H} NPY_VERSION
      REGEX "NPY_VERSION"
      )
    if (NPY_VERSION)
      string (REGEX REPLACE "#define NPY_VERSION 0x" "" NUMPY_VERSION ${NPY_VERSION})
    else (NPY_VERSION)
      message (STATUS "[NumPy] Unable to extract version from header file.")
    endif (NPY_VERSION)
  endif (NUMPY_NDARRAYOBJECT_H)
  
  find_file (NUMPY_TEST_PROGRAM TestNumPyVersion.cc
    PATHS ${CMAKE_MODULE_PATH} 
    )
  
  if (NUMPY_TEST_PROGRAM AND PYTHON_INCLUDE_DIRS)
    ## try to compile and run
    try_run (
      NUMPY_TEST_RUN_RESULT
      NUMPY_TEST_COMPILE_RESULT
      ${CMAKE_BINARY_DIR}
      ${NUMPY_TEST_PROGRAM}
      COMPILE_DEFINITIONS -I${PYTHON_INCLUDE_DIRS} -I${NUMPY_INCLUDES}
      COMPILE_OUTPUT_VARIABLE NUMPY_TEST_COMPILE
      RUN_OUTPUT_VARIABLE NUMPY_TEST_RUN
      )
    ## display results
    if (NOT NUMPY_FIND_QUIETLY)
      message (STATUS "[NumPy] NUMPY_TEST_RUN_RESULT     = ${NUMPY_TEST_RUN_RESULT}")
      message (STATUS "[NumPy] NUMPY_TEST_COMPILE_RESULT = ${NUMPY_TEST_COMPILE_RESULT}")
      message (STATUS "[NumPy] NUMPY_TEST_RUN            = ${NUMPY_TEST_RUN}")
    endif (NOT NUMPY_FIND_QUIETLY)
  else (NUMPY_TEST_PROGRAM AND PYTHON_INCLUDE_DIRS)
    message (STATUS "[NumPy] Unable to locate test program!")
  endif (NUMPY_TEST_PROGRAM AND PYTHON_INCLUDE_DIRS)
  
endif (numpy_version_test_output)

macro(numpy_feature_test OUTVARNAME testfilename testname)
  if (NUMPY_INCLUDES AND NUMPY_LIBRARIES AND NOT ${OUTVARNAME}) # only if numpy found.
    find_file(NUMPY_TESTFILE ${testfilename} PATHS ${CMAKE_MODULE_PATH})
    if (NOT NUMPY_TESTFILE)
      message (FATAL "[NumPy] ${testname} -- Unable to locate test program ${testfilename}!")
    endif(NOT NUMPY_TESTFILE)
    
    if (NUMPY_TESTFILE AND PYTHON_INCLUDE_DIRS)
      ## try to compile and run
      try_compile (
        ${OUTVARNAME}
        ${CMAKE_BINARY_DIR}
        ${NUMPY_TESTFILE}
        COMPILE_DEFINITIONS -I${PYTHON_INCLUDE_DIRS}  -I${NUMPY_INCLUDES} 
        CMAKE_FLAGS -DLINK_LIBRARIES:STRING=${PYTHON_LIBRARIES}
        OUTPUT_VARIABLE NUMPY_TESTCOMPILE
        )
      ## display results
      if (NOT NUMPY_FIND_QUIETLY)
        message (STATUS "[NumPy] ${testname} = ${OUTVARNAME} ${${OUTVARNAME}}")
      endif (NOT NUMPY_FIND_QUIETLY)
    endif (NUMPY_TESTFILE AND PYTHON_INCLUDE_DIRS)
    set(NUMPY_TESTFILE)
    set(NUMPY_TESTCOMPILE)
    mark_as_advanced(${OUTVARNAME})
  endif (NUMPY_INCLUDES AND NUMPY_LIBRARIES AND NOT ${OUTVARNAME})
endmacro()

numpy_feature_test(DCPROGS_NPY_LONG_DOUBLE
                   test_numpy_long_double.cc 
                   "Long double exists")
numpy_feature_test(DCPROGS_NPY_BOOL
                   test_numpy_ubyte.cc 
                   "Bool is a separate type")
numpy_feature_test(DCPROGS_NPY_ARRAY
                   test_numpy_is_noarray.cc 
                   "NPY_ARRAY_* macros exist")
numpy_feature_test(DCPROGS_NPY_ENABLEFLAGS
                   test_numpy_has_enableflags.cc 
                   "PyArray_ENABLEFLAGS exists")

## -----------------------------------------------------------------------------
## Actions taken when all components have been found

if (NUMPY_INCLUDES AND NUMPY_LIBRARIES)
  set (NUMPY_FOUND TRUE CACHE PATH "Says wether Numpy Python was found.")
  set (NUMPY_INCLUDES ${NUMPY_INCLUDES} CACHE PATH "Path to numpy C-API headers.")
  set (NUMPY_LIBRARIES ${NUMPY_LIBRARIES} CACHE PATH "Path to numpy C-API libraries.")
else (NUMPY_INCLUDES AND NUMPY_LIBRARIES)
  set (NUMPY_FOUND FALSE CACHE PATH "Says wether Numpy Python was found.")
  if (NOT NUMPY_FIND_QUIETLY)
    if (NOT NUMPY_INCLUDES)
      message (STATUS "[NumPy] Unable to find NUMPY header files!")
    endif (NOT NUMPY_INCLUDES)
    if (NOT NUMPY_LIBRARIES)
      message (STATUS "[NumPy] Unable to find NUMPY library files!")
    endif (NOT NUMPY_LIBRARIES)
  endif (NOT NUMPY_FIND_QUIETLY)
endif (NUMPY_INCLUDES AND NUMPY_LIBRARIES)

if (NUMPY_FOUND)
  if (NOT NUMPY_FIND_QUIETLY)
    message (STATUS "[NumPy] Found components for NUMPY")
    message (STATUS "[NumPy] NUMPY_INCLUDES  = ${NUMPY_INCLUDES}")
    message (STATUS "[NumPy] NUMPY_LIBRARIES = ${NUMPY_LIBRARIES}")
  endif (NOT NUMPY_FIND_QUIETLY)
else (NUMPY_FOUND)
  if (NUMPY_FIND_REQUIRED)
    message (FATAL_ERROR "[NumPy] Could not find NUMPY!")
  endif (NUMPY_FIND_REQUIRED)
endif (NUMPY_FOUND)

## -----------------------------------------------------------------------------
## Mark advanced variables

unset(NUMPY_ARRAYOBJECT_H CACHE)
unset(NUMPY_BOOLTEST_FILE CACHE)
unset(NUMPY_LONGTEST_FILE CACHE)
unset(NUMPY_NDARRAYOBJECT_H CACHE)
unset(NUMPY_VERSION_PY CACHE)
mark_as_advanced (
  NUMPY_INCLUDES
  NUMPY_LIBRARIES
  NUMPY_MULTIARRAY_LIBRARY
  NUMPY_SCALARMATH_LIBRARY
  NUMPY_API_VERSION)

## -----------------------------------------------------------------------------
## Reinstate the original value of CMAKE_FIND_LIBRARY_PREFIXES

set (CMAKE_FIND_LIBRARY_PREFIXES ${TMP_FIND_LIBRARY_PREFIXES} CACHE STRING
  "Library prefixes"
  FORCE
  )

# Get the major and minor version numbers
string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST "${NUMPY_API_VERSION}")
list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)
