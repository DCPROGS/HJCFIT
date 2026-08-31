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

cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)

include(CMakeDependentOption)

PROJECT(HJCFIT)

option(tests          "Enable testing."                         on)
option(pythonBindings "Enable python bindings."                 on)
option(compileDocs    "Compile C++ documentation examples."      on)
option(openmp         "Enable OpenMP parallelization."          on)
option(executenotebooks "Execute example Jupyter notebooks."    on)
# MPFR gives an arbitrary-precision fallback for root bracketing when double
# precision fails (likelihood/root_finder.cc). Default ON for Linux only.
#
# Not Apple: HJCFITConfig.h.in declares std::complex<mpfr::mpreal>, and the
# standard only specifies std::complex for float, double and long double.
# Instantiating it with anything else is undefined behaviour. libstdc++ has
# always tolerated it; libc++ enforces the rule, so on macOS the build fails
# inside <complex> and <cmath> with "floating point classification requires
# argument of floating point type (passed in 'mpfr::mpreal')".
#
# Turning it on there needs a complex type that does not rely on that
# guarantee, which is real work rather than a flag. Until then macOS builds
# without the MPFR fallback, which is a supported configuration -- this has
# always been an option, and Windows has never had it.
CMAKE_DEPENDENT_OPTION(HJCFIT_USE_MPFR "Enable fallback to MPFR multiple precision"
                       ON "UNIX;NOT APPLE" OFF)

set(JUPYTER_KERNEL "" CACHE STRING "Jupyter kernel to run notebooks with defaults to python2 or python3 depending on python version")
# Provides backwards compatibility and prevents MACOS_RPATH warning
if((CMAKE_MAJOR_VERSION GREATER 3) OR (CMAKE_MAJOR_VERSION EQUAL 3))
      set(CMAKE_MACOSX_RPATH 1)
endif()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# Check if our link directory is already an implicit link directory
# Only add it to RPath if not.
if(DEFINED SKBUILD_PLATLIB_DIR)
  # Wheel build: the extension and the library it needs sit in the same
  # directory inside the package, so look next to the loading binary rather
  # than at an absolute prefix that will not exist on the installing machine.
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
  else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
  endif()
elseif(NOT WIN32)
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif()
endif()
# Set a default build type to Release if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_SCRIPTS ${PROJECT_SOURCE_DIR}/cmake_modules)

if(pythonBindings)
  include(${CMAKE_SCRIPTS}/AllPythonBindings.cmake)
  # Unless explicitly initialized from command line or installation occurs in a
  # virtual environment, set install prefix to path to the python
  # interpreter.
  if (EXISTS $ENV{VIRTUAL_ENV})
    set(CMAKE_INSTALL_PREFIX "$ENV{VIRTUAL_ENV}" CACHE string "" FORCE)
  elseif (EXISTS $ENV{CONDA_ENV_PATH})
    set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_ENV_PATH}" CACHE string "" FORCE)
  elseif (WIN32 AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    FILE(TO_CMAKE_PATH "${PYTHON_INTERP_PREFIX}" PYTHON_INTERP_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${PYTHON_INTERP_PREFIX}" CACHE string "" FORCE)
  endif()
endif(pythonBindings)

# Eigen. Config mode: every current distribution and Homebrew ships
# Eigen3Config.cmake, which defines Eigen3::Eigen. Point CMAKE_PREFIX_PATH or
# Eigen3_DIR at a private build if you have one.
#
# Eigen3Config declares SameMajorVersion compatibility, so asking for 3.3
# *rejects* Eigen 4 and 5 rather than accepting them as newer. Homebrew now
# ships Eigen 5. Prefer a version this project has actually been built against,
# then accept a newer one with a warning rather than refusing to configure.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  find_package(Eigen3 REQUIRED NO_MODULE)
  message(WARNING
    "Found Eigen ${Eigen3_VERSION}. HJCFIT is tested against Eigen 3.3 and 3.4; "
    "other major versions are unverified.")
endif()
message(STATUS "[Eigen3] ${Eigen3_VERSION}")

# Eigen 3.x's config sets EIGEN3_INCLUDE_DIR, which likelihood/CMakeLists.txt
# uses. Newer configs may only define the target, so derive it if absent.
if(NOT EIGEN3_INCLUDE_DIR AND TARGET Eigen3::Eigen)
  get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
endif()

if(MSVC)
  # Two things CookOff's AddCPP11Flags did besides set C++ flags, and which
  # went missing when that include was replaced by CMAKE_CXX_STANDARD.
  #
  # MSWINDOBE gates the dllimport/dllexport machinery in HJCFITConfig.h.in.
  # Without it the 36 MSWINDOBE annotations in the headers expand to nothing,
  # the DLL exports no symbols, MSVC emits no import library, and linking the
  # SWIG module fails with
  #
  #   LNK1181: cannot open input file 'likelihood.lib'
  #
  # /wd4251 silences the warning MSVC issues for every STL member of an exported
  # class. It is noise here because the library and its consumers are built
  # together; it would matter if we shipped headers for a differently-configured
  # compiler.
  set(MSWINDOBE TRUE)
  add_compile_options(/EHsc /wd4251)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# These were probed one by one when C++11 support was still uneven across
# compilers. The project now requires C++11 outright, and every feature in the
# list is mandatory in that standard, so probing for them can only ever return
# true. They stay as defines because HJCFITConfig.h.in and the documentation
# examples consume them.
set(HAS_CXX11_TYPE_TRAITS TRUE)
set(HAS_CXX11_TRIVIAL_TYPE_TRAITS TRUE)
set(HAS_CXX11_CONSTEXPR TRUE)
set(HAS_CXX11_RANDOM_DEVICE TRUE)
set(HAS_CXX11_CONSTRUCTOR_DELEGATE TRUE)
set(HAS_CXX11_NOEXCEPT TRUE)
set(HAS_CXX11_INITIALIZATION TRUE)

# std::isnan is <cmath> in C++11; the alternative spellings existed for
# compilers that predate it.
set(CXX_HAS_STD_ISNAN TRUE)
if(openmp)
  FIND_PACKAGE(OpenMP)
endif()

if(OPENMP_FOUND)
  if(NOT (${CMAKE_C_FLAGS} MATCHES ${OpenMP_C_FLAGS}))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  endif()
  if(NOT (${CMAKE_CXX_FLAGS} MATCHES ${OpenMP_CXX_FLAGS}))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  endif()
  if(NOT ("${CMAKE_EXE_LINKER_FLAGS}" MATCHES "${OpenMP_EXE_LINKER_FLAGS}"))
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" CACHE STRING "" FORCE)
  endif()
elseif(openmp)
  MESSAGE(STATUS "Code will be compiled without OpenMP support.")
endif(OPENMP_FOUND)
if ( CMAKE_COMPILER_IS_GNUCC )
  set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
endif ( CMAKE_COMPILER_IS_GNUCC )
if(tests)
  enable_testing()
  include(${CMAKE_SCRIPTS}/CTestCustomPreWorkAround.cmake)
  # Resolves to cmake_modules/AddGTest.cmake, which shadows the CookOff module
  # of the same name: CookOff fetched googletest from a moving default branch
  # that no longer exists. Ours pins a commit. See that file for details.
  include(AddGTest)
endif(tests)

if(NOT HAS_CXX11_NOEXCEPT)
  set(noexcept TRUE)
endif()

if(NOT HJCFIT_LONG_DOUBLE)
  set(HJCFIT_LONG_DOUBLE False CACHE BOOL
      "If True, will use long doubles rather than simple doubles.")
endif(NOT HJCFIT_LONG_DOUBLE)

if(NOT HJCFIT_USE_MPFR)
  set(HJCFIT_USE_MPFR False CACHE BOOL
  "If True, will use MPFR arbitrary precision floats as fall back if regular double calculations fail.")
endif(NOT HJCFIT_USE_MPFR)

if(HJCFIT_USE_MPFR)
  find_path(GMP_INCLUDE_DIR gmp.h)
  find_library(GMP_LIBRARY NAMES gmp)
  find_path(MPFR_INCLUDE_DIR mpfr.h)
  find_library(MPFR_LIBRARY NAMES mpfr)
  if(NOT GMP_LIBRARY OR NOT MPFR_LIBRARY)
    message(FATAL_ERROR
      "HJCFIT_USE_MPFR is ON but GMP or MPFR was not found. Install them "
      "(libgmp-dev and libmpfr-dev on Debian and Ubuntu, gmp and mpfr on "
      "Homebrew) or configure with -DHJCFIT_USE_MPFR=OFF.")
  endif()
  mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY MPFR_INCLUDE_DIR MPFR_LIBRARY)
  include_directories(${PROJECT_SOURCE_DIR}/mpfr ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
endif(HJCFIT_USE_MPFR)
configure_file (
 "${PROJECT_SOURCE_DIR}/HJCFITConfig.h.in"
 "${PROJECT_BINARY_DIR}/HJCFITConfig.h"
)

# Save all HJCFIT headers in the same place in Windows
install(FILES ${PROJECT_BINARY_DIR}/HJCFITConfig.h DESTINATION include/HJCFIT)

include(${CMAKE_SCRIPTS}/documentation.cmake)

add_subdirectory(likelihood)
add_subdirectory(data)
add_subdirectory(python)
add_subdirectory(documentation)
