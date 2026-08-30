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
# Local replacement for GreatCMakeCookOff's PythonPackage module.
#
# Only find_python_package was ever used, and only to answer "is this module
# importable". CookOff's version could also download and install the package;
# that behaviour is deliberately not reproduced. A build system should report a
# missing dependency, not reach out and install one.

if(_HJCFIT_PYTHON_PACKAGE_INCLUDED)
  return()
endif()
set(_HJCFIT_PYTHON_PACKAGE_INCLUDED TRUE)

# find_python_package(<module> [REQUIRED] [QUIET])
#
# Sets <MODULE>_FOUND, and <MODULE>_VERSION when the module exposes
# __version__. Uses PYTHON_EXECUTABLE, which AllPythonBindings.cmake maps from
# Python3_EXECUTABLE.
function(find_python_package module)
  cmake_parse_arguments(fpp "REQUIRED;QUIET" "" "" ${ARGN})

  string(TOUPPER "${module}" MODULE)

  if(NOT PYTHON_EXECUTABLE)
    if(fpp_REQUIRED)
      message(FATAL_ERROR "[${module}] no python interpreter to check against.")
    endif()
    set(${MODULE}_FOUND FALSE PARENT_SCOPE)
    return()
  endif()

  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
            "import ${module}; print(getattr(${module}, '__version__', ''))"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _version
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(_result EQUAL 0)
    set(${MODULE}_FOUND TRUE PARENT_SCOPE)
    set(${MODULE}_VERSION "${_version}" PARENT_SCOPE)
    if(NOT fpp_QUIET)
      if(_version)
        message(STATUS "[${module}] ${_version}")
      else()
        message(STATUS "[${module}] found")
      endif()
    endif()
  else()
    set(${MODULE}_FOUND FALSE PARENT_SCOPE)
    if(fpp_REQUIRED)
      message(FATAL_ERROR "[${module}] not importable by ${PYTHON_EXECUTABLE}.")
    elseif(NOT fpp_QUIET)
      message(STATUS "[${module}] not found.")
    endif()
  endif()
endfunction()
