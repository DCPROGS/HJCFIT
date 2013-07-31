/***********************
    DCProgs computes missed-events likelihood as described in
    Hawkes, Jalali and Colquhoun (1990, 1992)

    Copyright (C) 2013  University College London

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
************************/

// Starts the likelihood sub-package
%module(package="dcprogs") likelihood
// C++ definitions that are needed to compile the python bindings.
%{
#define SWIG_FILE_WITH_INIT
#  include <DCProgsConfig.h>
#  include <iostream>
#  include <sstream>
#  include <memory>

#  if NUMPY_VERSION_MINOR >= 7 
#    define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#  endif
#  include <numpy/arrayobject.h>
#  ifndef DCPROGS_NPY_ARRAY
#    define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS 
#    define NPY_ARRAY_WRITEABLE    NPY_WRITEABLE
#  endif
#  ifndef DCPROGS_NPY_ENABLEFLAGS
#    define PyArray_ENABLEFLAGS(ARRAY, FLAGS)  (ARRAY)->flags |= FLAGS
#    define PyArray_CLEARFLAGS(ARRAY, FLAGS)   (ARRAY)->flags &= (!FLAGS)
#    define PyArray_SetBaseObject(ARRAY, BASE) (ARRAY)->base   = BASE
#  endif
#  include "object.h"
#  include "numpy_eigen.h"
#  include "helpers.h"

#  include "../qmatrix.h"
#  include "../idealG.h"
#  include "../occupancies.h"
#  include "../determinant_equation.h"
#  include "../root_finder.h"
#  include "../asymptotes.h"
#  include "../exact_survivor.h"
#  include "../approx_survivor.h"
#  include "../missed_eventsG.h"
#  include "../likelihood.h"
#  include "../time_filter.h"


#ifdef DCPROGS_CATCH
# error DCPROGS_CATCH already defined.
#endif 
#define DCPROGS_CATCH(ONERROR)                                                 \
    catch (DCProgs::errors::ComplexEigenvalues &_e) {                          \
      PyErr_SetString(PyExc_ArithmeticError, _e.what());                       \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::NaN &_e) {                                       \
      PyErr_SetString(PyExc_ArithmeticError, _e.what());                       \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::Mass &_e) {                                      \
      PyErr_SetString(PyExc_ArithmeticError, _e.what());                       \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::Domain &_e) {                                    \
      PyErr_SetString(PyExc_ArithmeticError, _e.what());                       \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::Math &_e) {                                      \
      PyErr_SetString( PyExc_ArithmeticError,                                  \
                       ( std::string("Math error in dcprogs: ")                \
                         + _e.what()).c_str() );                               \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::Index &_e) {                                     \
      PyErr_SetString(PyExc_IndexError, _e.what());                            \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::PythonTypeError &_e) {                           \
      PyErr_SetString(PyExc_TypeError, _e.what());                             \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::PythonValueError &_e) {                          \
      PyErr_SetString(PyExc_ValueError, _e.what());                            \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::PythonErrorAlreadyThrown &) {                    \
      ONERROR;                                                                 \
    } catch (DCProgs::errors::Python &_e) {                                    \
      PyErr_SetString(PyExc_RuntimeError, _e.what());                          \
      ONERROR;                                                                 \
    } catch(DCProgs::errors::NotImplemented &_e) {                             \
      PyErr_SetString(PyExc_NotImplementedError, _e.what());                   \
      ONERROR;                                                                 \
    } catch(DCProgs::errors::Root &_e) {                                       \
      PyErr_SetString( PyExc_RuntimeError,                                     \
                       (std::string("Encountered unknonw error in DCProgs\n")  \
                        + _e.what()).c_str() );                                \
      ONERROR;                                                                 \
    } catch(std::exception &_e) {                                              \
      PyErr_SetString(PyExc_RuntimeError, _e.what());                          \
      ONERROR;                                                                 \
    } catch(...) {                                                             \
      PyErr_SetString(PyExc_RuntimeError, "Caught unknown exception.");        \
      ONERROR;                                                                 \
    }
%}


// Tells swig that we will deal with exceptions.
%include "exception.i"
%init %{ import_array();  %}

%exception {
  try { $function }
  DCPROGS_CATCH(return NULL;);
}

// Tells swig about our type hierarchy. 
// These types should make it easier to go from one system to another, but they do make it slightly
// more difficult for swig to understand our code.
%apply int { DCProgs::t_int }; 
%apply double { DCProgs::t_real }; 
// unsigned should be checked for negative values. Otherwise swig comes out with a truly unhelpfull
// message. Also checks for numpy scalar at the same type.
%typemap(in) DCProgs::t_uint {
  try {
    using namespace DCProgs;
    long integer;
    if(PyLong_Check($input)) integer = PyLong_AsLong($input);
#   ifndef DCPROGS_PYTHON3
      else if(PyInt_Check($input)) integer = PyInt_AsLong($input);
#   endif
    else if(PyArray_CheckScalar($input))  {
      PyArray_CastScalarToCtype($input, static_cast<void*>(&integer), 
                                PyArray_DescrFromType(numpy::type<long>::value));
      if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
    } else throw errors::PythonTypeError("Unexpected type in burst.");
    if(integer < 0) throw errors::PythonValueError("Input integer cannot be negative.");
    $1 = static_cast<t_uint>(integer);
  } DCPROGS_CATCH(SWIG_fail);
}
%typemap(out) DCProgs::t_uint = unsigned int;
%typemap(typecheck) DCProgs::t_uint = int;

// Adds some standard converters for eigen, 
// + bindings for svd, eig, ... in case we are compiling with t_real > double.
%include "math.swg"
// These macros help us translate from C++ exceptions to python exceptions
//! General namespace for all things DCProgs.
namespace DCProgs {

  %include "qmatrix.swg"
  %include "idealg.swg"
  %include "determinant_equation.swg"
  %include "root_finder.swg"
  %include "asymptotes.swg"
  %include "exact_survivor.swg"
  %include "approx_survivor.swg"
  %include "missed_eventsG.swg"
  %include "log10likelihood.swg"

}
%include "time_filter.swg"
%include "chained.swg"

#undef DCPROGS_CATCH
