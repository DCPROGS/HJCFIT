add_library(likelihood SHARED state_matrix.h errors.h idealG.h idealG.cc equilibrium.h equilibrium.cc)

# Windows + Visual studio specifics.
if(MSVC)
  set_target_properties(likelihood PROPERTIES DEFINE_SYMBOL DCPROGS_LIKELIHOOD_DLLEXPORT)
endif(MSVC)

# Put all libraries and output in same place.
# Otherwise windobe cannot find the dll it generate...
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/likelihood)

add_dependencies(likelihood eigen)

# Test subdirectory.
add_subdirectory(tests)
add_subdirectory(python)
