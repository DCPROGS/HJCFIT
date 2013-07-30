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

%{
  namespace {
    //! Computes likelihood of a time series in a numpy array.
    //! \param[in] _intervals: Time intervals, starting and ending with an "open" interval.
    //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
    //!                member function, where the argument is the length of an open or shut interval.
    //! \param[in] _initial: initial occupancies.
    //! \param[in] _final: final occupancies.
    template<class T_G> 
      DCProgs::t_real py_chained_likelihood(T_G const & _g, PyObject *_intervals, 
                                            PyObject *_initial, PyObject *_final) {
        using namespace DCProgs;
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
  
#       ifdef DCPROGS_MACRO
#         error DCPROGS_MACRO is already defined.
#       endif
#       define DCPROGS_MACRO(TYPENAME)                                                        \
           if(type == numpy::type<TYPENAME>::value)                                           \
             return chained_likelihood( _g, static_cast<TYPENAME*>(first),                    \
                                            static_cast<TYPENAME*>(first) + size,             \
                                        initial, final ); 
          DCPROGS_MACRO(npy_double);
          DCPROGS_MACRO(npy_float);
          DCPROGS_MACRO(npy_longlong);
          DCPROGS_MACRO(npy_ulonglong);
          DCPROGS_MACRO(npy_long);
          DCPROGS_MACRO(npy_ulong);
          DCPROGS_MACRO(npy_int);
          DCPROGS_MACRO(npy_uint);
          DCPROGS_MACRO(npy_short);
          DCPROGS_MACRO(npy_ushort);
          DCPROGS_MACRO(npy_byte);
          DCPROGS_MACRO(npy_ubyte);
#         ifdef DCPROGS_NPY_HAS_LONG_DOUBLE
            DCPROGS_MACRO(npy_longdouble);
#         endif
#         ifdef DCPROGS_NPY_HAS_BOOL
            DCPROGS_MACRO(npy_bool);
#         endif
#       undef DCPROGS_MACRO
        throw errors::PythonTypeError("Unexpect numpy array type");
        return t_real(0);
      }
    DCProgs::t_real chained_likelihood(DCProgs::IdealG *_idealg, PyObject* _intervals, 
                                       PyObject* _initial, PyObject* _final) {
      return py_chained_likelihood(*_idealg, _intervals, _initial, _final);
    } 
    DCProgs::t_real chained_likelihood(DCProgs::MissedEventsG *_idealg, PyObject* _intervals, 
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
        chained_likelihood(IdealG *_idealg, PyObject *_intervals, 
                           PyObject* _initial, PyObject* _final);
DCProgs::t_real chained_likelihood(DCProgs::IdealG *idealg, PyObject* intervals, 
                                   PyObject* initial, PyObject* final);
DCProgs::t_real chained_likelihood(DCProgs::MissedEventsG *g, PyObject* intervals, 
                                   PyObject* initial, PyObject* final);
