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
  # PREFIX is set explicitly rather than left to UseSWIG. The legacy
  # LIBRARY_OUTPUT_PATH in likelihood/CMakeLists.txt puts every target in one
  # directory, and under the Visual Studio generator the module came out as
  # likelihood rather than _likelihood, colliding with the C++ library it links
  # against:
  #
  #   LNK1149: output filename matches input filename
  #
  # OUTPUT_NAME stays the swig module name, which UseSWIG passes to swig as
  # -module; the leading underscore belongs to the file, not the module.
  set_target_properties(likelihood_python PROPERTIES
                        OUTPUT_NAME likelihood
                        PREFIX "_"
                        IMPORT_PREFIX "_"
                        LINKER_LANGUAGE CXX)
  # Python3::Module, not ${PYTHON_LIBRARIES}. On Unix an extension module must
  # NOT link libpython: its symbols are resolved from the interpreter that
  # loads it. Linking it made delocate treat the framework as an external
  # dependency and copy it into the wheel, so importing loaded a second Python
  # runtime into a process that already had one, and segfaulted:
  #
  #   Copying library .../Python.framework/Versions/3.12/Python to HJCFIT/.dylibs/Python
  #   import HJCFIT.likelihood -> failed with code -11
  #
  # Python3::Module links the import library on Windows, where that is required,
  # and nothing on Unix, which is what each platform actually wants.
  target_link_libraries(likelihood_python Python3::Module likelihood)

  
  install(TARGETS likelihood_python DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood)
  # Named rather than globbed, because CMake will not re-run on a new file
  # appearing and a glob would make the wheel depend on when cmake last ran.
  set(PYTHON_MODULES __init__.py random.py _methods.py optimization.py
                     simplex_hjc.py fitting.py fitspec.py runner.py cli.py)

  install(FILES  ${PYTHON_MODULES}
                 ${CMAKE_CURRENT_BINARY_DIR}/likelihood.py
          DESTINATION ${PYINSTALL_DIRECTORY}/HJCFIT/likelihood)

  # And checked, because a hand-maintained list falls behind silently. fitting.py
  # was added to this directory without being added above, so it was absent from
  # the wheel; the build succeeded, every C++ and behave test passed, and the
  # only thing that noticed was an integration test that imported it. A module
  # missing from the package is not a build failure unless something asks for it.
  file(GLOB PYTHON_MODULES_PRESENT RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
  foreach(module ${PYTHON_MODULES_PRESENT})
    list(FIND PYTHON_MODULES ${module} module_listed)
    if(module_listed EQUAL -1)
      message(FATAL_ERROR
        "${module} is in likelihood/python but not in PYTHON_MODULES, so it "
        "would not be installed. Add it to that list in "
        "likelihood/python/CMakeLists.txt.")
    endif()
  endforeach()
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
