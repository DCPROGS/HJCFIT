add_library(likelihood SHARED qmatrix.h qmatrix.cc errors.h idealG.h idealG.cc 
            occupancies.h occupancies.cc time_filter.h time_filter.cc recursion_formula.h
            exact_survivor.h exact_survivor.cc root_finder.h root_finder.cc determinant_equation.h
            determinant_equation.cc asymptotes.h asymptotes.cc approx_survivor.h approx_survivor.cc
            missed_eventsG.h missed_eventsG.cc laplace_survivor.h laplace_survivor.cc )

# Windows + Visual studio specifics.
if(MSVC)
  set_target_properties(likelihood PROPERTIES DEFINE_SYMBOL DCPROGS_LIKELIHOOD_DLLEXPORT)
endif(MSVC)

# Put all libraries and output in same place.
# Otherwise windobe cannot find the dll it generate...
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)

add_dependencies(likelihood eigen)

install(TARGETS likelihood DESTINATION lib) 

# Test subdirectory.
add_subdirectory(tests)
add_subdirectory(python)
