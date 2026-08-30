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
#
# Local replacement for GreatCMakeCookOff's AddGTest.
#
# Originally written to shadow GreatCMakeCookOff's module of the same name.
# CookOff has since been removed entirely, so this is simply the project's
# googletest integration.
#
# Why it exists
# -------------
# CookOff's LookUpGTest.cmake calls ExternalProject_Add for googletest with a
# repository but no GIT_TAG. ExternalProject defaults GIT_TAG to "master".
# googletest renamed its default branch to "main" and no longer has a "master"
# branch at all, so the download step fails:
#
#     [4/63] Performing download step (git clone) for 'Lookup-GTest'
#     FAILED: external/src/Lookup-GTest-stamp/Lookup-GTest-download
#     fatal: Remote branch master not found in upstream origin
#
# That broke the build without a single commit landing in HJCFIT, which is the
# hazard of resolving dependencies against a moving default branch. The fix is
# to pin.
#
# Why v1.18.0
# -----------
# This was pinned to release-1.12.1, the last googletest supporting C++11, for
# as long as HJCFIT compiled as C++11. The project now targets C++17, so the
# constraint is gone and the pin moves to a current release.
#
# Pinned to a commit SHA rather than the tag, because tags can be moved and a
# SHA cannot.

include(FetchContent)

if(CMAKE_VERSION VERSION_LESS 3.14)
  message(FATAL_ERROR
    "FetchContent_MakeAvailable requires CMake 3.14 or newer (found "
    "${CMAKE_VERSION}). Either upgrade CMake or set -Dtests=OFF.")
endif()

set(HJCFIT_GOOGLETEST_REPOSITORY "https://github.com/google/googletest.git"
    CACHE STRING "googletest repository to fetch")
set(HJCFIT_GOOGLETEST_GIT_TAG "063de7e9578f82b369302001269680b4b1553359"
    CACHE STRING "googletest commit to build against (v1.18.0)")
mark_as_advanced(HJCFIT_GOOGLETEST_REPOSITORY HJCFIT_GOOGLETEST_GIT_TAG)

# googletest build knobs. Must be set before FetchContent_MakeAvailable, since
# that is when googletest's own CMakeLists reads them.
#   BUILD_GMOCK        the tests use gtest only; skipping gmock halves the build
#   INSTALL_GTEST      do not let a test dependency add itself to `make install`
#   gtest_force_shared_crt  MSVC: match the CRT the rest of the project uses
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# GIT_SHALLOW is deliberately left off: a shallow fetch of a bare commit SHA is
# not supported by every server configuration, and reliability matters more here
# than the few seconds saved.
FetchContent_Declare(
  googletest
  GIT_REPOSITORY ${HJCFIT_GOOGLETEST_REPOSITORY}
  GIT_TAG        ${HJCFIT_GOOGLETEST_GIT_TAG}
)

# For an offline or air-gapped build, point CMake at an existing checkout:
#   cmake -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/path/to/googletest ..
FetchContent_MakeAvailable(googletest)

find_package(Threads)

# add_gtest(<name> <sources> [libraries...])
#
# <sources> may be a single file or a semicolon-separated list, e.g.
#   add_gtest(qmatrix qmatrix.cc likelihood)
#   add_gtest(brentq "brentq.cc;random_matrix.h;random_matrix.cc" likelihood)
#
# Every test in likelihood/tests/ defines its own main() calling
# InitGoogleTest and RUN_ALL_TESTS, so these link gtest and not gtest_main.
macro(add_gtest name source)
  add_executable(test_${name} ${source})
  target_link_libraries(test_${name} PRIVATE gtest)

  if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(test_${name} PRIVATE ${CMAKE_THREAD_LIBS_INIT})
  endif()

  if(NOT "${ARGN}" STREQUAL "")
    target_link_libraries(test_${name} PRIVATE ${ARGN})
  endif()

  # Referring to the target rather than a path under EXECUTABLE_OUTPUT_PATH lets
  # CMake resolve the location itself, including per-config subdirectories on
  # multi-config generators.
  add_test(NAME cxx_${name}
           COMMAND test_${name}
           --gtest_output=xml:${CMAKE_BINARY_DIR}/test-results/test_${name}.xml)
  set_tests_properties(cxx_${name} PROPERTIES LABELS "gtest")
endmacro()
