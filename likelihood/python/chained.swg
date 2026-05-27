/***********************
    HJCFIT computes missed-events likelihood as described in
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

%{
  namespace {
    //! Computes likelihood of a time series in a numpy array.
    //! \param[in] _intervals: Time intervals, starting and ending with an "open" interval.
    //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
    //!                member function, where the argument is the length of an open or shut interval.
    //! \param[in] _initial: initial vectors.
    //! \param[in] _final: final vectors.
    template<class T_G> 
      HJCFIT::t_real py_chained_likelihood(T_G const & _g, PyObject *_intervals, 
                                            PyObject *_initial, PyObject *_final) {
        using namespace HJCFIT;
        if( (not PyArray_Check(_intervals))
            or (not (PyArray_FLAGS((PyArrayObject*)_intervals) & NPY_ARRAY_C_CONTIGUOUS)) ) {
          Object<PyObject> convert = steal_ref( (PyObject*)
            PyArray_FromObject(_intervals, numpy::type<t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
          return py_chained_likelihood(_g, ~convert, _initial, _final);
        }
        t_initvec const initial = numpy::map_to_initvec(_initial).array();
        t_rvector const final   = numpy::map_to_rvector(_final).array();
        int const type = PyArray_TYPE((PyArrayObject*)_intervals);
        npy_intp const size = PyArray_SIZE((PyArrayObject*)_intervals);
        void * const first = PyArray_DATA((PyArrayObject*)_intervals);
  
#       ifdef HJCFIT_MACRO
#         error HJCFIT_MACRO is already defined.
#       endif
#       define HJCFIT_MACRO(TYPENAME)                                                        \
           if(type == numpy::type<TYPENAME>::value)                                           \
             return chained_likelihood( _g, static_cast<TYPENAME*>(first),                    \
                                            static_cast<TYPENAME*>(first) + size,             \
                                        initial, final ); 
          HJCFIT_MACRO(npy_double);
          HJCFIT_MACRO(npy_float);
          HJCFIT_MACRO(npy_longlong);
          HJCFIT_MACRO(npy_ulonglong);
          HJCFIT_MACRO(npy_long);
          HJCFIT_MACRO(npy_ulong);
          HJCFIT_MACRO(npy_int);
          HJCFIT_MACRO(npy_uint);
          HJCFIT_MACRO(npy_short);
          HJCFIT_MACRO(npy_ushort);
          HJCFIT_MACRO(npy_byte);
          HJCFIT_MACRO(npy_ubyte);
#         ifdef NUMPY_NPY_LONG_DOUBLE
            HJCFIT_MACRO(npy_longdouble);
#         endif
#         ifdef NUMPY_NPY_LONG_DOUBLE
            HJCFIT_MACRO(npy_bool);
#         endif
#       undef HJCFIT_MACRO
        throw errors::PythonTypeError("Unexpect numpy array type");
        return t_real(0);
      }
    HJCFIT::t_real chained_likelihood(HJCFIT::IdealG *_idealg, PyObject* _intervals, 
                                       PyObject* _initial, PyObject* _final) {
      return py_chained_likelihood(*_idealg, _intervals, _initial, _final);
    } 
    HJCFIT::t_real chained_likelihood(HJCFIT::MissedEventsG *_idealg, PyObject* _intervals, 
                                       PyObject* _initial, PyObject* _final) {
      return py_chained_likelihood(*_idealg, _intervals, _initial, _final);
    } 
  }
%}

%feature("autodoc", "Computes likelihood of a series of time intervals.\n\n"
         ":param firstArgument: \n"
         "   Can be either an :class:`IdealG` or and :class:`MissedEventsG` instance.\n"
         ":param intervals:\n"
         "   A list of lists of time intervals. Each inner list should be a burst.\n"
         ":param initial:\n"
         "   The initial occupancy vector.\n" 
         ":param final:\n"
         "   The final occupancy vector.\n")
        chained_likelihood(HJCFIT::IdealG *_idealg, PyObject *_intervals, 
                           PyObject* _initial, PyObject* _final);
HJCFIT::t_real chained_likelihood(HJCFIT::IdealG *idealg, PyObject* intervals, 
                                   PyObject* initial, PyObject* final);
HJCFIT::t_real chained_likelihood(HJCFIT::MissedEventsG *g, PyObject* intervals, 
                                   PyObject* initial, PyObject* final);
