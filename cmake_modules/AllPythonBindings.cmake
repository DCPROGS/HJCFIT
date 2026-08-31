# Python bindings are a bit messy, so done here rather than main file.
#
# Previously this used GreatCMakeCookOff's CoherentPython and Numpy find modules
# and read the install directory out of distutils. CMake's own FindPython3 does
# all of it, is maintained, and does not need a network fetch at configure time.
# distutils was removed from the standard library in Python 3.12, which made the
# old path fail silently: execute_process does not check its result, so
# PYTHON_PKG_DIR simply came back empty and the install paths collapsed without
# any error (F4).

include(PythonPackage)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# The rest of the build, and the SWIG glue in particular, was written against
# the older PYTHON_* spellings. Map them rather than churn every call site.
set(PYTHON_EXECUTABLE   ${Python3_EXECUTABLE})
set(PYTHON_LIBRARIES    ${Python3_LIBRARIES})
set(PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
set(PYTHON_INCLUDE_PATH ${Python3_INCLUDE_DIRS})
set(PYTHON_VERSION      "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
set(NUMPY_INCLUDE_DIRS  ${Python3_NumPy_INCLUDE_DIRS})
message(STATUS "[Python] ${Python3_VERSION} at ${Python3_EXECUTABLE}")
message(STATUS "[NumPy] ${Python3_NumPy_VERSION}")

# sys.prefix, which the Windows branch of the install-prefix logic in the top-
# level CMakeLists uses. CookOff's CoherentPython supplied this as
# PYTHON_INTERP_PREFIX; nothing replaced it when CookOff was removed, and the
# reference is inside a WIN32 guard so no Linux or macOS build could reach it.
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.prefix)"
  OUTPUT_VARIABLE PYTHON_INTERP_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _prefix_result)
if(NOT _prefix_result EQUAL 0)
  message(FATAL_ERROR "Could not determine the python interpreter prefix.")
endif()

# HJCFITConfig.h.in consumes these. NUMPY_VERSION_MINOR alone is not enough to
# decide anything now that NumPy is on 2.x - see the gate in likelihood.i (F3).
string(REPLACE "." ";" _npy_version_list "${Python3_NumPy_VERSION}")
list(GET _npy_version_list 0 NUMPY_VERSION_MAJOR)
list(GET _npy_version_list 1 NUMPY_VERSION_MINOR)
unset(_npy_version_list)

# NPY_ARRAY_* and PyArray_ENABLEFLAGS arrived in NumPy 1.7 (2013). The minimum
# NumPy that works with this project is far beyond that, so these are constants.
set(NUMPY_NPY_ARRAY TRUE)
set(NUMPY_NPY_ENABLEFLAGS TRUE)

# These two are not "does the macro exist" questions, despite the names. They
# ask whether npy_longdouble and npy_bool are types *distinct from* npy_double
# and npy_ubyte. numpy_eigen.h builds a compile-time type -> NPY_* table, so if
# a pair collapses to the same C++ type the specialisations collide:
#
#   numpy_eigen.h:59: error: redefinition of 'type<double>'
#   numpy_eigen.h:375: error: duplicate case value
#
# Both answers are platform-dependent and neither can be assumed. numpy
# typedefs npy_longdouble to double where long double carries no extra
# precision, which is the case on Apple Silicon; and npy_bool is unsigned char
# on every platform seen so far, which is why numpy_eigen.h:69 has an #else
# branch. Probe rather than guess.
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_INCLUDES ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})

check_cxx_source_compiles("
#include <Python.h>
#include <numpy/arrayobject.h>
#include <type_traits>
int main() {
  static_assert(!std::is_same<npy_double, npy_longdouble>::value, \"same type\");
  return 0;
}" NUMPY_NPY_LONG_DOUBLE)

check_cxx_source_compiles("
#include <Python.h>
#include <numpy/arrayobject.h>
#include <type_traits>
int main() {
  static_assert(!std::is_same<npy_ubyte, npy_bool>::value, \"same type\");
  return 0;
}" NUMPY_NPY_BOOL)

unset(CMAKE_REQUIRED_INCLUDES)
message(STATUS "[NumPy] long double distinct from double = ${NUMPY_NPY_LONG_DOUBLE}")
message(STATUS "[NumPy] bool distinct from unsigned char = ${NUMPY_NPY_BOOL}")

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(HJCFIT_PYTHON3 TRUE)

# Where a plain "pip install" would put a pure-python package.
if(NOT DEFINED PYTHON_PKG_DIR)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
    OUTPUT_VARIABLE PYTHON_PKG_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pkg_dir_result)
  if(NOT _pkg_dir_result EQUAL 0)
    message(FATAL_ERROR "Could not determine the python package directory.")
  endif()
  set(PYTHON_PKG_DIR ${PYTHON_PKG_DIR} CACHE PATH "Main python package repository.")
  mark_as_advanced(PYTHON_PKG_DIR)
endif()

# There is an issue on Windows where pyconfig.h defines a macro hypot that screws up swig+c++11
# Test for issue and add -include cmath otherwise
if(MSYS)
  file(WRITE  ${CMAKE_BINARY_DIR}/test_cmath_python.cc
       "#include <Python.h>\n"
       "#include <cmath>\n"
       "int main() { return 0; }" )
  try_compile(
    NEED_CMATH_INCLUDE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/test_cmath_python.cc
    COMPILE_DEFINITIONS -I${PYTHON_INCLUDE_DIRS}
                        -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
    CMAKE_FLAGS -DLINK_LIBRARIES:STRING=${PYTHON_LIBRARIES}
                -DCMAKE_CXX_FLAGS_DEBUG:STRING="${CMAKE_CXX_FLAGS_RELEASE}"
                -DCMAKE_C_FLAGS_DEBUG:STRING="${CMAKE_C_FLAGS_RELEASE}"
                -DCMAKE_EXE_LINKER_FLAGS_DEBUG:STRING="${CMAKE_EXE_LINKER_FLAGS_RELEASE}"
    OUTPUT_VARIABLE OUTVAR
  )
  file(REMOVE ${CMAKE_BINARY}/test_cmath_python.cc)
  if(NOT NEED_CMATH_INCLUDE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include cmath")
  else(NOT NEED_CMATH_INCLUDE)
    message(FATAL_ERROR "[CXX PYTHON] checks this does not get tripped outside MSYS ${NEED_CMATH_INCLUDE}")
  endif(NOT NEED_CMATH_INCLUDE)
  unset(NEED_CMATH_INCLUDE)
endif(MSYS)

set(HJCFIT_PYTHON_BINDINGS True)

if(tests)
  if(NOT DEFINED TEST_INSTALL_DIRECTORY)
    set(TEST_INSTALL_DIRECTORY test-results/install)
    set(TEST_INSTALL_ABSPATH ${CMAKE_BINARY_DIR}/${TEST_INSTALL_DIRECTORY})
  endif(NOT DEFINED TEST_INSTALL_DIRECTORY)


  # A macro to run tests via python or behave.
  function(_python_test name filename thiscommand)
    if(WIN32)
      set(WORKINGDIR ${TEST_INSTALL_ABSPATH}/lib/site-packages)
      set(ADD_TO_PATH ${TEST_INSTALL_ABSPATH}/Library/bin)
    else()
      set(WORKINGDIR ${TEST_INSTALL_ABSPATH}/lib/python${PYTHON_VERSION}/site-packages)
      set(ADD_TO_PATH ${TEST_INSTALL_ABSPATH}/lib)
    endif(WIN32)

    add_test(NAME python_${name}
             WORKING_DIRECTORY ${WORKINGDIR}
             COMMAND ${thiscommand} ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${ARGN})
    if(MSVC OR MSYS)
      set_tests_properties(python_${name} PROPERTIES CONFIGURATIONS Release)
      set(PATH_STRING "${ADD_TO_PATH};$ENV{PATH}")
      STRING(REPLACE "\\;" ";" PATH_STRING "${PATH_STRING}")
      STRING(REPLACE ";" "\\;" PATH_STRING "${PATH_STRING}")
      file(TO_NATIVE_PATH "${WORKINGDIR}" PYTHON_PATH)
      STRING(REPLACE "\\;" ";" PYTHONPATH "${PYTHON_PATH};$ENV{PYTHONPATH}")
      STRING(REPLACE ";" "\\;" PYTHON_PATH "${PYTHON_PATH}")
      set_tests_properties(python_${name} PROPERTIES ENVIRONMENT
                           "PYTHONPATH=${PYTHON_PATH};PATH=${PATH_STRING}")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set_tests_properties(python_${name} PROPERTIES ENVIRONMENT
                           "PYTHONPATH=${WORKINGDIR}:$ENV{PYTHONPATH};DYLD_LIBRARY_PATH=${ADD_TO_PATH}:$ENV{DYLD_LIBRARY_PATH} ")
    elseif(UNIX)
      set_tests_properties(python_${name} PROPERTIES ENVIRONMENT
                           "LD_LIBRARY_PATH=${ADD_TO_PATH}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${WORKINGDIR}:$ENV{PYTHONPATH}")
    endif(MSVC OR MSYS)
  endfunction(_python_test)

  # Look for behave
  if(NOT BEHAVE_EXECUTABLE)
    find_program(BEHAVE_EXECUTABLE behave DOC "Path to the behave executable")
    if(NOT BEHAVE_EXECUTABLE)
      message(FATAL_ERROR "[behave] Not found. Cannot run python tests.")
    endif(NOT BEHAVE_EXECUTABLE)
    message(STATUS "[behave] ${BEHAVE_EXECUTABLE}")
  endif(NOT BEHAVE_EXECUTABLE)
  # A macro to run tests via behave.
  function(feature_test name filename)
    _python_test(${name} ${filename} ${BEHAVE_EXECUTABLE} --junit --junit-directory
                 ${CMAKE_BINARY_DIR}/test-results/ -q ${ARGN})

  endfunction(feature_test)

  function(python_test name filename)
    _python_test(${name} ${filename} ${PYTHON_EXECUTABLE} ${ARGN})
  endfunction(python_test)
endif(tests)


# Follow conda convention for Python installation in Windows
if(WIN32)
  set(PYINSTALL_DIRECTORY lib/site-packages)
else()
  set(PYINSTALL_DIRECTORY lib/python${PYTHON_VERSION}/site-packages)
endif(WIN32)

