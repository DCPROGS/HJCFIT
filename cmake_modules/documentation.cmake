#-- Add an Option to toggle the generation of the API documentation
include(PythonPackage)
find_package(Doxygen)
if(pythonBindings)
  find_package(Sphinx)
  find_package(Pandoc)
  find_package(NbConvert)
  find_python_package(ipykernel)
endif(pythonBindings)
set(GENERATE_XML NO)
set(SPHINX_SOURCE_DIR ${PROJECT_SOURCE_DIR}/documentation/source)
if(SPHINX_FOUND)
  set(GENERATE_XML YES)
  if(DOXYGEN_FOUND)
    set(SPHINX_EXTENSIONS "'sphinxcontrib.bibtex', 'breathe'")
  else()
    set(SPHINX_EXTENSIONS "'sphinxcontrib.bibtex'")
  endif(DOXYGEN_FOUND)
  if(DCPROGS_USE_MPFR)
    set(SPHINX_DCPROGS_USE_MPFR "True")
  else()
    set(SPHINX_DCPROGS_USE_MPFR "False")
  endif(DCPROGS_USE_MPFR)
endif(SPHINX_FOUND)

if (DOXYGEN_FOUND)
  configure_file( ${PROJECT_SOURCE_DIR}/documentation/doxygen.in
                  ${CMAKE_BINARY_DIR}/Doxyfile
                  @ONLY IMMEDIATE )
  if(NOT SPHINX_FOUND)
    message(STATUS "[Doxygen] found. C++ Documentation target is \"doxydocs\".")
  endif(NOT SPHINX_FOUND)
else()
  message(STATUS "[Doxygen] not found. Cannot build documentation.")
endif(DOXYGEN_FOUND)

if(SPHINX_FOUND)

  if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME default)
  endif()
  if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME_DIR)
  endif()
  if(NOT DEFINED SPHINX_BUILDER)
    set(SPHINX_BUILDER html)
  endif(NOT DEFINED SPHINX_BUILDER)

  file(TO_NATIVE_PATH "${SPHINX_SOURCE_DIR}" SPHINX_SOURCE_DIR) 


  configure_file(
    "${PROJECT_SOURCE_DIR}/documentation/source/conf.py"
    "${CMAKE_BINARY_DIR}/conf.py"
    @ONLY)

  message(STATUS "[Sphinx] found. Documentation target is \"documentation\".")
  if(PANDOC_FOUND AND NBCONVERT_FOUND)
    message(STATUS "Pandoc and nbconvert found. Will integrate Jupyer notebooks in docs")
    if(NOT IPYKERNEL_FOUND)
      message(STATUS "IPyKernel not found cannot execute notebooks. Will include static notebooks")
    endif()
  else()
    message(STATUS "Pandoc and nbconvert not found. Will not integrate Jupyer notebooks in docs")
  endif(PANDOC_FOUND AND NBCONVERT_FOUND)
else()
  message(STATUS "[Sphinx] not found. Documentation cannot be built.")
endif(SPHINX_FOUND)
