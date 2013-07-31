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

if(pythonBindings)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_SWIG_FLAGS "")

  set_source_files_properties(likelihood.i PROPERTIES CPLUSPLUS ON)
  if(PYTHON_VERSION VERSION_LESS "3.0.0")
    set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall" )
  else()
    set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall;-py3" )
  endif(PYTHON_VERSION VERSION_LESS "3.0.0")
  set(SWIG_MODULE_likelihood_EXTRA_DEPS determinant_equation.swg qmatrix.swg root_finder.swg
                                        asymptotes.swg exact_survivor.swg approx_survivor.swg
                                        missed_eventsG.swg idealg.swg chained.swg time_filter.swg
                                        math.swg brentq.swg log10likelihood.swg )
  swig_add_module(likelihood python likelihood.i)
  set_target_properties(_likelihood PROPERTIES LINKER_LANGUAGE CXX)
  swig_link_libraries(likelihood ${PYTHON_LIBRARIES} likelihood)

  
  install(TARGETS _likelihood DESTINATION ${PYINSTALL_DIRECTORY}/dcprogs/likelihood)
  install(FILES  random.py _methods.py optimization.py
                 ${CMAKE_CURRENT_BINARY_DIR}/likelihood.py __init__.py
          DESTINATION ${PYINSTALL_DIRECTORY}/dcprogs/likelihood)
endif(pythonBindings)

if(pythonBindings AND tests)

  feature_test(qmatrix features/qmatrix.feature)
  feature_test(determinant_equation features/determinant_equation.feature)
  feature_test(find_roots features/find_roots.feature)
  feature_test(asymptotes features/asymptotes.feature)
  feature_test(idealg features/idealg.feature)
  feature_test(missed_eventsG features/missed_eventsG.feature)
  feature_test(approx_survivor features/approx_survivor.feature)
  feature_test(reduce_likelihood features/reduce_likelihood.feature)
  feature_test(time_filter features/time_filter.feature)
  
endif(pythonBindings AND tests)
