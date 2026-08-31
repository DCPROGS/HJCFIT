########################
#   HJCFIT computes missed-events likelihood as described in
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

set(LIKELIHOOD_SOURCE qmatrix.cc idealG.cc vectors.cc time_filter.cc exact_survivor.cc
                      root_finder.cc determinant_equation.cc asymptotes.cc approx_survivor.cc
                      missed_eventsG.cc laplace_survivor.cc brentq.cc likelihood.cc)
set(LIKELIHOOD_HEADERS qmatrix.h errors.h idealG.h vectors.h time_filter.h recursion_formula.h
                       exact_survivor.h root_finder.h determinant_equation.h asymptotes.h
                       approx_survivor.h missed_eventsG.h laplace_survivor.h likelihood.h brentq.h)

include_directories(${PROJECT_BINARY_DIR})
if(EIGEN3_INCLUDE_DIR)
  include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

add_library(likelihood SHARED ${LIKELIHOOD_SOURCE} ${LIKELIHOOD_HEADERS})

if(HJCFIT_USE_MPFR)
  # Link the libraries found in the top-level CMakeLists by full path. The bare
  # name "mpfr" only resolved because Hunter or the system happened to put it on
  # the default link path, which is not true of Homebrew's prefix.
  target_link_libraries(likelihood ${MPFR_LIBRARY} ${GMP_LIBRARY})
endif(HJCFIT_USE_MPFR)

# Windows + Visual studio specifics.
if(MSVC)
  set_target_properties(likelihood PROPERTIES DEFINE_SYMBOL HJCFIT_LIKELIHOOD_DLLEXPORT)
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
if(DEFINED SKBUILD_PLATLIB_DIR)
  # In a wheel the library ships inside the package on every platform, since a
  # wheel cannot install to lib/. The RPATH set at the top level makes the
  # extension look beside itself for it.
  install(TARGETS likelihood
          RUNTIME DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood
          LIBRARY DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood
          ARCHIVE DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood)
elseif(WIN32 AND pythonBindings)
  # A .pyd resolves its dependent DLLs from its own directory and from whatever
  # os.add_dll_directory() was given -- since Python 3.8 it does not search
  # PATH. The conda layout put likelihood.dll in Library/bin and the extension
  # in site-packages, so importing HJCFIT.likelihood failed with
  #
  #   ImportError: DLL load failed while importing _likelihood
  #
  # even though both files installed correctly. The test harness compensated by
  # adding Library/bin to PATH, which stopped working at Python 3.8 and is why
  # nothing caught it.
  #
  # Putting the DLL beside the extension is also what a wheel requires, since a
  # wheel cannot install anything outside its own package directory.
  install(TARGETS likelihood
          RUNTIME DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood
          ARCHIVE DESTINATION lib)
else()
  install(TARGETS likelihood DESTINATION ${LIB_DIRECTORY})
endif()
install(FILES ${LIKELIHOOD_HEADERS} DESTINATION include/HJCFIT/likelihood)


# Test subdirectory.
add_subdirectory(tests)
add_subdirectory(python)
