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

if(pythonBindings)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
  set(CMAKE_SWIG_FLAGS "")

  set_source_files_properties(likelihood.i PROPERTIES CPLUSPLUS ON)
  if(PYTHON_VERSION VERSION_LESS "3.0.0")
    set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall" )
  else()
    set_source_files_properties(likelihood.i PROPERTIES SWIG_FLAGS "-includeall" )
  endif(PYTHON_VERSION VERSION_LESS "3.0.0")
  set(SWIG_MODULE_likelihood_python_EXTRA_DEPS
                                        determinant_equation.swg qmatrix.swg root_finder.swg
                                        asymptotes.swg exact_survivor.swg approx_survivor.swg
                                        missed_eventsG.swg idealg.swg chained.swg time_filter.swg
                                        math.swg brentq.swg log10likelihood.swg )
  # The C++ shared library target is already called "likelihood". swig_add_module used to
  # prefix an underscore and so produced "_likelihood"; swig_add_library uses the name
  # verbatim and would collide. Hence a distinct target name. OUTPUT_NAME carries the SWIG
  # module name and must NOT include the leading underscore: UseSWIG passes it to swig as
  # -module, and swig's python backend adds the underscore itself when emitting the
  # "import _likelihood" in the generated likelihood.py.
  swig_add_library(likelihood_python LANGUAGE python SOURCES likelihood.i)
  set_target_properties(likelihood_python PROPERTIES
                        OUTPUT_NAME likelihood
                        LINKER_LANGUAGE CXX)
  target_link_libraries(likelihood_python ${PYTHON_LIBRARIES} likelihood)

  
  install(TARGETS likelihood_python DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood)
  install(FILES  random.py _methods.py optimization.py
                 ${CMAKE_CURRENT_BINARY_DIR}/likelihood.py __init__.py
          DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood)
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
