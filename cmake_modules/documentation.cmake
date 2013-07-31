#-- Add an Option to toggle the generation of the API documentation
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
  configure_file( ${PROJECT_SOURCE_DIR}/documentation/doxygen.in
                  ${PROJECT_BINARY_DIR}/Doxyfile
                  @ONLY IMMEDIATE )
  add_custom_target ( doxydocs 
                      COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
                      SOURCES ${PROJECT_BINARY_DIR}/Doxyfile )
  message(STATUS "[Doxygen] found. \">make doxydocs\" will create c++ documentation.")
else()
  message(STATUS "[Doxygen] not found. Cannot build documentation.")
endif(DOXYGEN_FOUND)
