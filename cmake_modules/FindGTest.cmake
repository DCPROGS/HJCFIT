# Add gtest
ExternalProject_Add(
    googletest
    SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
    TIMEOUT 10
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CMAKE_ARGS 
      -DBUILD_SHARED_LIBS=OFF
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

find_package(Threads)

macro(cxx_test name source)
  add_executable(test_${name} ${source})

  ExternalProject_Get_Property(googletest source_dir)
  include_directories(${source_dir}/include)
  # Better, but only works on CMake 2.8.6?
  # get_target_property(THISTEST_INCLUDE test_${name} INCLUDE_DIRECTORIES)
  # set_target_properties(test_${name} PROPERTIES INCLUDE_DIRECTORIES
  #                       "${source_dir}/include;${THISTEST_INCLUDE}") 

  ExternalProject_Get_Property(googletest binary_dir)
  target_link_libraries(test_${name} ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
  target_link_libraries(test_${name} ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)
  target_link_libraries(test_${name} ${CMAKE_THREAD_LIBS_INIT})

  add_dependencies(test_${name} googletest)
  if(NOT "${ARGN}" STREQUAL "")
    target_link_libraries(test_${name} ${ARGN})
  endif(NOT "${ARGN}" STREQUAL "")

  add_test(test_${name} test_${name} --gtest_output=xml:${CMAKE_BINARY_DIR}/test-results/test_${name}.xml)
endmacro()
