if(NOT WIN32)
  # Newer cmakes can't cope with CTEST_CUSTOM_PRE_TEST having arguments, at
  # least on Mac and Linux. Hence the BS below.
  file(WRITE ${CMAKE_BINARY_DIR}/fakeinstall.cc
       "// Workaround for CTEST_CUSTOM_POST_TEST not allowing arguments \n"
       "#include <stdlib.h> \n"
       "int main() \n"
       "{ \n"
       "  return system(\""
            "${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${TEST_INSTALL_DIRECTORY} "
                             "-P ${CMAKE_BINARY_DIR}/cmake_install.cmake \");"
       "}\n")

  add_executable(fake_test_install ${CMAKE_BINARY_DIR}/fakeinstall.cc)

  file(WRITE ${CMAKE_BINARY_DIR}/CTestCustom.cmake
    "set(CTEST_CUSTOM_PRE_TEST ${CMAKE_BINARY_DIR}/fake_test_install)\n")
else(NOT WIN32)
  file(WRITE ${CMAKE_BINARY_DIR}/CTestCustom.cmake
       "set(CTEST_CUSTOM_PRE_TEST \"\\\"${CMAKE_COMMAND}\\\""
  	   "  -DCMAKE_INSTALL_PREFIX=${TEST_INSTALL_DIRECTORY}"
  	   "  -P \\\"${CMAKE_BINARY_DIR}/cmake_install.cmake\\\"\")\n")
endif(NOT WIN32)
