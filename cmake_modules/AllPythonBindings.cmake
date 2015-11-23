# Python bindings are a bit messy, so done here rather than main file.

# Install Python with GreatCMakeCookOff, as well as numpy and behave
include(PythonPackage)
include(PythonPackageLookup)
find_package(CoherentPython)
find_python_package(numpy)
find_package(Numpy REQUIRED)
#find_python_package(behave)

if(NOT PYTHON_VERSION AND PYTHONINTERP_FOUND)
  execute_process( 
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(\"%i.%i\" % sys.version_info[:2])"
    OUTPUT_VARIABLE PYTHON_VERSION
  )
  if( PYTHON_VERSION )
    string (STRIP ${PYTHON_VERSION} PYTHON_VERSION)
    set(PYTHON_VERSION ${PYTHON_VERSION} CACHE STRING "Version of the Python interpreter.")
    mark_as_advanced(PYTHON_VERSION)
    MESSAGE(STATUS "[Python] Version: ${PYTHON_VERSION}")
  else( PYTHON_VERSION )
    MESSAGE(STATUS "Could not determine python version.")
  endif( PYTHON_VERSION )
endif(NOT PYTHON_VERSION AND PYTHONINTERP_FOUND)

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
  if(PYTHON_PKG_DIR)
    string (STRIP ${PYTHON_PKG_DIR} PYTHON_PKG_DIR)
    set(PYTHON_PKG_DIR ${PYTHON_PKG_DIR} CACHE PATH "Main python package repository.")
    mark_as_advanced(PYTHON_PKG_DIR)
  endif(PYTHON_PKG_DIR)
endif(NOT DEFINED PYTHON_PKG_DIR)

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

set(DCPROGS_PYTHON_BINDINGS True)

if(tests)
  if(NOT DEFINED TEST_INSTALL_DIRECTORY)
    set(TEST_INSTALL_DIRECTORY test-results/install)
    set(TEST_INSTALL_ABSPATH ${CMAKE_BINARY_DIR}/${TEST_INSTALL_DIRECTORY})
  endif(NOT DEFINED TEST_INSTALL_DIRECTORY)

  include(${CMAKE_SCRIPTS}/CTestCustomPreWorkAround.cmake)

  # A macro to run tests via python or behave.
  function(_python_test name filename thiscommand)
    if(WIN32)
      set(WORKINGDIR ${TEST_INSTALL_ABSPATH}/dcprogs/python-pkg/)
      set(ADD_TO_PATH ${TEST_INSTALL_ABSPATH}/dcprogs/DLLs)
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
                           "PATH=${PATH_STRING};PYTHONPATH=${PYTHON_PATH}")
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



if(WIN32)
  set(PYINSTALL_DIRECTORY dcprogs/python-pkg)
else()
  set(PYINSTALL_DIRECTORY lib/python${PYTHON_VERSION}/site-packages)
endif(WIN32)

if(NOT PYTHON_VERSION VERSION_LESS "3.0.0")
  set(DCPROGS_PYTHON3 TRUE)
endif(NOT PYTHON_VERSION VERSION_LESS "3.0.0")

