# Python bindings are a bit messy, so done here rather than main file.

find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c
    "from sys import version_info; print(version_info.major==3)"
  OUTPUT_VARIABLE PYTHON_VERSION_MAJOR
)

if(NOT DEFINED PYTHON_PKG_DIR)
  execute_process( 
    COMMAND ${PYTHON_EXECUTABLE} -c 
              "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
              OUTPUT_VARIABLE PYTHON_PKG_DIR
  )
  if(PYTHON_PKG_DIR )
    string (STRIP ${PYTHON_PKG_DIR} PYTHON_PKG_DIR)
    set(PYTHON_PKG_DIR ${PYTHON_PKG_DIR} CACHE PATH
        "Python modules will be installed here." )
  else(PYTHON_PKG_DIR)
    set( PYTHON_PKG_DIR lib/python${PYTHON_VERSION}/site-packages
         CACHE PATH "Python modules will be installed here." )
  endif(PYTHON_PKG_DIR)
  mark_as_advanced(PYTHON_PKG_DIR)
  MESSAGE(STATUS "[Python] installation directory: ${PYTHON_PKG_DIR}")
endif(NOT DEFINED PYTHON_PKG_DIR)
if(NOT DEFINED CMAKE_PYINSTALL_PREFIX)
  if(WIN32)
    set(CMAKE_PYINSTALL_PREFIX python-pkg/dcprogs)
  else(WIN32)
    set(CMAKE_PYINSTALL_PREFIX lib/python${PYTHON_VERSION}/site-packages/dcprogs)
  endif(WIN32)
  set(CMAKE_PYINSTALL_PREFIX ${CMAKE_PYINSTALL_PREFIX} CACHE
      PATH "Installation path of the python package")
endif(NOT DEFINED CMAKE_PYINSTALL_PREFIX)

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

find_package(numpy REQUIRED) 

set(DCPROGS_PYTHON_BINDINGS True)

if(tests AND pythonBindings)
  if(NOT DEFINED TEST_INSTALL_DIRECTORY)
    set(TEST_INSTALL_DIRECTORY ${CMAKE_BINARY_DIR}/tests/install 
        CACHE PATH "Path of a fake install for testing purposes")
    mark_as_advanced(TEST_INSTALL_DIRECTORY)
  endif(NOT DEFINED TEST_INSTALL_DIRECTORY)

  file(TO_NATIVE_PATH ${TEST_INSTALL_DIRECTORY} TESTNATDIR)
  file(WRITE ${CMAKE_BINARY_DIR}/CTestCustom.cmake
       "set(CTEST_CUSTOM_PRE_TEST \"\\\"${CMAKE_COMMAND}\\\""
             "  -DCMAKE_INSTALL_PREFIX=${TESTNATDIR}"
             "  -P \\\"${CMAKE_BINARY_DIR}/cmake_install.cmake\\\"\")\n" )
  unset(TESTNATDIR)

  # A macro to run tests via behave.
  function(feature_test name filename)
    add_test(NAME python_${name} 
             WORKING_DIRECTORY ${TEST_INSTALL_DIRECTORY}/${CMAKE_PYINSTALL_PREFIX}/..
             COMMAND behave ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -q --summary ${ARGN})
  endfunction(feature_test)
endif(tests AND pythonBindings)
