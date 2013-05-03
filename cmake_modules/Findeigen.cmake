if(eigen_INCLUDE_DIR)
  set( eigen_FIND_QUIETLY True)
endif(eigen_INCLUDE_DIR)

find_path(eigen_INCLUDE_DIR
  NAMES
  Eigen/Core
  Eigen/LU
  Eigen/Geometry
  Eigen/Cholesky
  PATHS
  $ENV{eigen_INCLUDE_DIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  eigen3 eigen2 eigen
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(eigen DEFAULT_MSG eigen_INCLUDE_DIR)

# gets wether eigen2 or eigen3
if(EIGEN_FOUND) 
  file(WRITE ${PROJECT_BINARY_DIR}/pylada_dummy.cc
        "#include <Eigen/Core>\n"  
        "#if not EIGEN_VERSION_AT_LEAST(3,0,0)\n"
        "#  error\n" 
        "#endif\n" )
  try_compile(is_eigen3  ${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/pylada_dummy.cc
              COMPILE_DEFINITIONS -I${eigen_INCLUDE_DIR} 
              CMAKE_FLAGS -DCMAKE_CXX_LINK_EXECUTABLE="echo" 
              OUTPUT_VARIABLES _is_eigen3)
  file(REMOVE ${PROBJECT_BINARY_DIR}/pylada_dummy.cc)
endif(EIGEN_FOUND)
