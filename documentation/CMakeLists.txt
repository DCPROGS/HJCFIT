if (DOXYGEN_FOUND)
  add_custom_target ( doxydocs 
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
    SOURCES ${CMAKE_BINARY_DIR}/Doxyfile )
endif(DOXYGEN_FOUND)

if(SPHINX_FOUND)

  if(WIN32)
    set(WORKINGDIR ${TEST_INSTALL_ABSPATH}/dcprogs/python-pkg/)
    set(ADD_TO_PATH ${TEST_INSTALL_ABSPATH}/dcprogs/DLLs)
  else()
    set(WORKINGDIR ${TEST_INSTALL_ABSPATH}/lib/python${PYTHON_VERSION}/site-packages)
    set(ADD_TO_PATH ${TEST_INSTALL_ABSPATH}/lib)
  endif(WIN32)

  set(SOURCE_DOC_FILES source/bibliography.bib source/api/cpp/exceptions.rst
                       source/api/cpp/idealG.rst source/api/cpp/log10likelihood.rst
                       source/api/cpp/missed_eventsG.rst source/api/cpp/occupancies.rst
                       source/api/cpp/qmatrix.rst source/api/cpp/typehierarchy.rst
                       source/api/cpp.rst source/api/python/idealG.rst
                       source/api/python/log10likelihood.rst source/api/python/missed_eventsG.rst
                       source/api/python/qmatrix.rst source/api/python.rst
                       source/api.rst source/index.rst source/manual/exact_survivor.rst
                       source/manual/idealG.rst source/manual/likelihood_of_Q.rst
                       source/manual/missedeventsG.rst source/manual/occupancies.rst
                       source/manual/qmatrix.rst source/manual.rst)

  add_custom_target(documentation 
    ${SPHINX_EXECUTABLE}
    -q -b ${SPHINX_BUILDER}
    -c ${CMAKE_BINARY_DIR}
    "${CMAKE_CURRENT_SOURCE_DIR}/source"
    "${CMAKE_BINARY_DIR}/documentation/sphinx"
    COMMENT "Building HTML documentation with Sphinx"
    WORKING_DIRECTORY ${WORKINGDIR} 
    DEPENDS )
#   DEPENDS doxydocs )

endif(SPHINX_FOUND)

if(compileDocs)

  file(GLOB CXX_SOURCE_FILES code/*.cc)
  file(GLOB PY_SOURCE_FILES code/*.py)

  foreach(pathfile ${CXX_SOURCE_FILES})
    get_filename_component(pathname ${pathfile} NAME_WE)
    add_executable(doc_cxx_${pathname} ${pathfile})
    target_link_libraries(doc_cxx_${pathname} likelihood)
    include_directories(${PROJECT_SOURCE_DIR})
  endforeach()

  if(tests)
    foreach(pathfile ${CXX_SOURCE_FILES})
      get_filename_component(pathname ${pathfile} NAME_WE)
      add_test(cxx_doctest_${pathname} doc_cxx_${pathname})
    endforeach()

    foreach(pathfile ${PY_SOURCE_FILES})
      get_filename_component(pathname ${pathfile} NAME_WE)
      file(RELATIVE_PATH pathfile ${CMAKE_CURRENT_SOURCE_DIR} ${pathfile})
      python_test(doctest_py_${pathname} ${pathfile})
    endforeach()
  endif(tests)
endif(compileDocs)
