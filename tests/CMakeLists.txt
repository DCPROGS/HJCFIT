# Enable ExternalProject CMake module
INCLUDE(ExternalProject)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/external)

# Add gtest
ExternalProject_Add(
    googletest
    SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
    TIMEOUT 10
    # Force separate output paths for debug and release builds to allow easy
    # identification of correct lib in subsequent TARGET_LINK_LIBRARIES commands
    CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

function(cxx_test name source)
  add_executable(test_${name} ${source})

  ExternalProject_Get_Property(googletest source_dir)
  set_target_properties(test_${name} PROPERTIES INCLUDE_DIRECTORIES "${source_dir}/include")

  ExternalProject_Get_Property(googletest binary_dir)
  target_link_libraries(test_${name} ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
  target_link_libraries(test_${name} ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a)

  add_dependencies(test_${name} googletest)

  add_test(test_${name} test_${name} --gtest_output=xml:test_${name}.xml)
endfunction()

cxx_test(hello hello.cc)
