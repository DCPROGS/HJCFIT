%{
  //! Computes likelihood of a time series in a numpy array.
  //! \param[in] _intervals: Time intervals, starting and ending with an "open" interval.
  //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
  //!                member function, where the argument is the length of an open or shut interval.
  //! \param[in] _initial: initial occupancies.
  //! \param[in] _final: final occupancies.
  template<class T_G> 
    DCProgs::t_real py_chained_likelihood(T_G const & _g, PyObject *_intervals, 
                                          PyObject *_initial, PyObject *_final) {
      if( (not PyArray_Check(_intervals))
          or (not (PyArray_FLAGS((PyArrayObject*)_intervals) & NPY_ARRAY_C_CONTIGUOUS)) ) {
        Object<PyObject> convert = steal_ref( (PyObject*)
          PyArray_FromObject(_intervals, DCProgs::numpy::type<DCProgs::t_real>::value, 0, 0)
        );
        if(PyErr_Occurred()) throw DCProgs::errors::PythonErrorAlreadyThrown();
        return py_chained_likelihood(_g, ~convert, _initial, _final);
      }
      DCProgs::t_rvector const initial = DCProgs::numpy::map_to_rmatrix(_initial).array();
      DCProgs::t_rvector const final   = DCProgs::numpy::map_to_rmatrix(_final).array();
      int const type = PyArray_TYPE((PyArrayObject*)_intervals);
      npy_intp const size = PyArray_SIZE((PyArrayObject*)_intervals);
      void * const first = PyArray_DATA((PyArrayObject*)_intervals);

#     ifdef DCPROGS_MACRO
#       error DCPROGS_MACRO is already defined.
#     endif
#     define DCPROGS_MACRO(TYPENAME)                                                             \
         if(type == DCProgs::numpy::type<TYPENAME>::value)                                       \
           return DCProgs::chained_likelihood( static_cast<TYPENAME*>(first),                    \
                                               static_cast<TYPENAME*>(first) + size,             \
                                               _g, initial, final ); 
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
#       ifdef DCPROGS_NPY_HAS_LONG_DOUBLE
          DCPROGS_MACRO(npy_longdouble);
#       endif
#       ifdef DCPROGS_NPY_HAS_BOOL
          DCPROGS_MACRO(npy_bool);
#       endif
#     undef DCPROGS_MACRO
      throw DCProgs::errors::PythonTypeError("Unexpect numpy array type");
      return DCProgs::t_real(0);
    }
  DCProgs::t_real chained_likelihood(PyObject* _intervals, DCProgs::IdealG *_idealg, 
                                     PyObject* _initial, PyObject* _final) {
    return py_chained_likelihood(*_idealg, _intervals, _initial, _final);
  } 
  DCProgs::t_real chained_likelihood(PyObject* _intervals, DCProgs::MissedEventsG *_idealg, 
                                     PyObject* _initial, PyObject* _final) {
    return py_chained_likelihood(*_idealg, _intervals, _initial, _final);
  } 
%}

%feature("autodoc", "Computes likelihood of a series of time intervals.")
        chained_likelihood(IdealG *_idealg, PyObject *_intervals, 
                           t_rvector const &_initial, t_rvector const &_final);
DCProgs::t_real chained_likelihood(DCProgs::IdealG *idealg, PyObject* intervals, 
                                   PyObject* initial, PyObject* final);
DCProgs::t_real chained_likelihood(DCProgs::MissedEventsG *g, PyObject* intervals, 
                                   PyObject* initial, PyObject* final);
