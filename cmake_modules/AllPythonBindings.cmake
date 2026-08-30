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

# HJCFITConfig.h.in consumes these. NUMPY_VERSION_MINOR alone is not enough to
# decide anything now that NumPy is on 2.x - see the gate in likelihood.i (F3).
string(REPLACE "." ";" _npy_version_list "${Python3_NumPy_VERSION}")
list(GET _npy_version_list 0 NUMPY_VERSION_MAJOR)
list(GET _npy_version_list 1 NUMPY_VERSION_MINOR)
unset(_npy_version_list)

# NPY_ARRAY_* and PyArray_ENABLEFLAGS arrived in NumPy 1.7 (2013) and the bool
# and long double type numbers are older still. The minimum NumPy that works
# with this project is far beyond all of them, so probing is pointless; the
# defines stay because the SWIG layer reads them.
set(NUMPY_NPY_ARRAY TRUE)
set(NUMPY_NPY_ENABLEFLAGS TRUE)
set(NUMPY_NPY_LONG_DOUBLE TRUE)

# NUMPY_NPY_BOOL is deliberately left undefined. Despite the name it does not
# ask whether NPY_BOOL exists, but whether npy_bool is a type distinct from
# unsigned char. numpy typedefs it to unsigned char, so defining both
# numpy::type<npy_bool> and numpy::type<npy_ubyte> is a redefinition and gives
# duplicate case labels in numpy_eigen.h. CookOff probed this and reported
# "Bool is a separate type = FALSE"; numpy_eigen.h:69 has an #else branch for
# exactly this case.

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

