include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
  "#include <cmath>\nint main() { bool a = std::isnan(0e0); return 0; }\n" 
  CXX_HAS_STD_ISNAN)
if(NOT CXX_HAS_STD_ISNAN)
  CHECK_CXX_SOURCE_COMPILES(
    "#include <math.h>\nint main() { bool a = ::isnan(0e0); return 0; }\n" 
    CXX_HAS_ISNAN)
endif(NOT CXX_HAS_STD_ISNAN)
if(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)
  CHECK_CXX_SOURCE_COMPILES(
    "#include <math.h>\nint main() { bool a = _isnan(0e0); return 0; }\n" 
    CXX_HAS___ISNAN)
endif(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)
# if(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)
#   CHECK_CXX_SOURCE_COMPILES(
#     "# include <float.h>\nint main() { bool a = _isnan(0e0); return 0; }\n" 
#     CXX_HAS__ISNAN)
# endif(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN)
if(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN AND NOT CXX_HAS___ISNAN)
  message(FATAL_ERROR "[isnan] could not find standard function on this OS.")
endif(NOT CXX_HAS_STD_ISNAN AND NOT CXX_HAS_ISNAN AND NOT CXX_HAS___ISNAN)
#
