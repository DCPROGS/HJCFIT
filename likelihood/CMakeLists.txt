########################
#   DCProgs computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

set(LIKELIHOOD_SOURCE qmatrix.cc idealG.cc occupancies.cc time_filter.cc exact_survivor.cc
                      root_finder.cc determinant_equation.cc asymptotes.cc approx_survivor.cc
                      missed_eventsG.cc laplace_survivor.cc brentq.cc likelihood.cc)
set(LIKELIHOOD_HEADERS qmatrix.h errors.h idealG.h occupancies.h time_filter.h recursion_formula.h
                       exact_survivor.h root_finder.h determinant_equation.h asymptotes.h
                       approx_survivor.h missed_eventsG.h laplace_survivor.h likelihood.h brentq.h)

include_directories(${PROJECT_BINARY_DIR})
if(EIGEN3_INCLUDE_DIR)
  include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

add_library(likelihood SHARED ${LIKELIHOOD_SOURCE} ${LIKELIHOOD_HEADERS})
add_dependencies(likelihood lookup_dependencies)

if(DCPROGS_USE_MPFR)
  target_link_libraries(likelihood ${GMP_LIBRARIES})
  target_link_libraries(likelihood ${MPFR_LIBRARIES})
endif(DCPROGS_USE_MPFR)

# Windows + Visual studio specifics.
if(MSVC)
  set_target_properties(likelihood PROPERTIES DEFINE_SYMBOL DCPROGS_LIKELIHOOD_DLLEXPORT)
endif(MSVC)

# Put all libraries and output in same place.
# Otherwise windobe cannot find the dll it generate...
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)

# Follow conda convention for Python installation in Windows
if(WIN32)
  set(LIB_DIRECTORY Library/bin)
else()
  set(LIB_DIRECTORY lib)
endif(WIN32)
install(TARGETS likelihood DESTINATION ${LIB_DIRECTORY})
install(FILES ${LIKELIHOOD_HEADERS} DESTINATION include/dcprogs/likelihood)


# Test subdirectory.
add_subdirectory(tests)
add_subdirectory(python)
