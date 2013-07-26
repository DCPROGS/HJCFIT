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

add_library(likelihood SHARED qmatrix.h qmatrix.cc errors.h idealG.h idealG.cc 
            occupancies.h occupancies.cc time_filter.h time_filter.cc recursion_formula.h
            exact_survivor.h exact_survivor.cc root_finder.h root_finder.cc determinant_equation.h
            determinant_equation.cc asymptotes.h asymptotes.cc approx_survivor.h approx_survivor.cc
            missed_eventsG.h missed_eventsG.cc laplace_survivor.h laplace_survivor.cc likelihood.h
            brentq.cc brentq.h)

# Windows + Visual studio specifics.
if(MSVC)
  set_target_properties(likelihood PROPERTIES DEFINE_SYMBOL DCPROGS_LIKELIHOOD_DLLEXPORT)
endif(MSVC)

# Put all libraries and output in same place.
# Otherwise windobe cannot find the dll it generate...
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)

add_dependencies(likelihood eigen)

if(WIN32)
  set(LIB_DIRECTORY dcprogs/DLLs)
else()
  set(LIB_DIRECTORY lib)
endif(WIN32)
install(TARGETS likelihood DESTINATION ${LIB_DIRECTORY}) 

# Test subdirectory.
add_subdirectory(tests)
add_subdirectory(python)
