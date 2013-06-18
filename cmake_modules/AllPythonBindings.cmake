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

if(NOT PYINSTALL_DIR)
  execute_process( 
    COMMAND ${PYTHON_EXECUTABLE} -c 
              "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE PYINSTALL_DIR
  )
  if(PYINSTALL_DIR )
    string (STRIP ${PYINSTALL_DIR} PYINSTALL_DIR)
    set(PYINSTALL_DIR ${PYINSTALL_DIR} CACHE PATH "Version of the Python interpreter.")
  else()
    set( PYINSTALL_DIR lib/python${PYTHON_VERSION}/site-packages/ 
         CACHE PATH "Python modules will be installed here." )
  endif(PYINSTALL_DIR)
  mark_as_advanced(PYINSTALL_DIR)
  MESSAGE(STATUS "[Python] installation directory: ${PYINSTALL_DIR}")
endif(NOT PYINSTALL_DIR)

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
