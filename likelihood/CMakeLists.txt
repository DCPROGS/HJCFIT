add_library(likelihood SHARED state_matrix.h errors.h idealG.h idealG.cc 
            equilibrium.h equilibrium.cc time_filter.h time_filter.cc recursion_formula.h
            state_matrix.cc exact_survivor.h exact_survivor.cc root_finder.h root_finder.cc
            determinant_equation.h determinant_equation.cc 
            asymptotes.h asymptotes.cc approx_survivor.h approx_survivor.cc)

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
