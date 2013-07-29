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



%{

  namespace {
    DCProgs::t_Burst convert_pyburst(PyObject *_burst) {
      using namespace DCProgs;
      if(not PySequence_Check(_burst))
        throw errors::PythonTypeError("Expected a sequence on input.");
      t_Burst result;
      Object<> py_iterator = steal_ref(PyObject_GetIter(_burst));
      if(not py_iterator) throw errors::PythonErrorAlreadyThrown();
      while(Object<> item = steal_ref(PyIter_Next(~py_iterator))) {
          if(PyFloat_Check(~item)) 
            result.emplace_back(static_cast<t_real>(PyFloat_AS_DOUBLE(~item)));
          else if(PyInt_Check(~item))
            result.emplace_back(static_cast<t_real>(PyInt_AS_LONG(~item)));
          else if(PyLong_Check(~item))
            result.emplace_back(static_cast<t_real>(PyLong_AsDouble(~item)));
          else throw errors::PythonTypeError("Unexpected type in burst.");
      }
      if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
      return result;
    }
    DCProgs::t_Bursts convert_pybursts(PyObject *_bursts) {
      using namespace DCProgs;
      if(not PySequence_Check(_bursts))
        throw errors::PythonTypeError("Expected a sequence on input.");
      t_Bursts result;
      Object<> py_iterator = steal_ref(PyObject_GetIter(_bursts));
      if(not py_iterator) throw errors::PythonErrorAlreadyThrown();
      while(Object<> item = steal_ref(PyIter_Next(~py_iterator))) {
        if(not PySequence_Check(~item))
          throw errors::PythonTypeError("Expected a sequence of sequences on input.");

        t_Burst burst = convert_pyburst(~item);
        if(burst.size() == 0) continue;
        if(burst.size() % 2 != 1)  
          throw errors::PythonValueError("Each burst must have an odd number of intervals.");
        result.emplace_back(std::move(burst));
      }
  
      if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
      return result;
    }
  }
%}

namespace DCProgs {
  class Log10Likelihood {
    public:
      //! Number of open states.
      DCProgs::t_uint nopen;
      //! Max length of missed events
      DCProgs::t_real tau;
      //! \brief tcrit. 
      //! \detail If negative or null, will use equilibrium occupancies rather than CHS occupancies.
      DCProgs::t_real tcritical;
      //! Number of intervals for which to compute exact result.
      DCProgs::t_uint nmax;
      //! Tolerance for root finding.
      DCProgs::t_real xtol;
      //! Tolerance for root finding.
      DCProgs::t_real rtol;
      //! Maximum number of iterations for root finding.
      DCProgs::t_uint itermax;

      // The $self argument below is quite important. It tells SWIG that the memory for the matrix is
      // owned by c++, and that the c++ object should stay alive as long as the matrix lives in python.
      %typemap(out) DCProgs::t_Burst & { 
        try {
          Eigen::Map<DCProgs::t_rvector> map(result->data(), result->size());
          $result = DCProgs::numpy::wrap_to_numpy(map, $self); 
        } DCPROGS_CATCH(SWIG_fail);
      }
      //! Constructor
      %extend {
        Log10Likelihood   ( PyObject *_bursts, DCProgs::t_uint _nopen,
                            DCProgs::t_real _tau, DCProgs::t_real _tcritical=-1e0,
                            DCProgs::t_uint _nmax=2, DCProgs::t_real _xtol=1e-10,
                            DCProgs::t_real _rtol=1e-10, DCProgs::t_uint _itermax=100 )  {
          using namespace DCProgs;
          return new Log10Likelihood( convert_pybursts(_bursts), _nopen, _tau, _tcritical,
                                      _nmax, _xtol, _rtol, _itermax );
        }

        DCProgs::t_Burst & __getitem__(int i) {
          using namespace DCProgs;
          if(i < 0) i += $self->bursts.size();
          if(i < 0 or i >= $self->bursts.size())
            throw errors::Index("Index to burst is out of range.");
          return $self->bursts[i];
        }
        void __setitem__(int i, PyObject* _burst) {
          using namespace DCProgs;
          if(i < 0) i += $self->bursts.size();
          if(i < 0 or i >= $self->bursts.size())
            throw errors::Index("Index to burst is out of range.");
          t_Burst burst = convert_pyburst(_burst);
          if(burst.size() == 0) return;
          if(burst.size() % 2 != 1)  
            throw errors::PythonValueError("Each burst must have an odd number of intervals.");
          $self->bursts.insert($self->bursts.begin() + i, burst);
        }
        PyObject* __delitem__(int i) {
          using namespace DCProgs;
          if(i < 0) i += $self->bursts.size();
          if(i < 0 or i >= $self->bursts.size())
            throw errors::Index("Index to burst is out of range.");
          t_Burst out = *$self->bursts.erase($self->bursts.begin() + i);
          Eigen::Map<t_rvector> map(out.data(), out.size());
          return numpy::wrap_to_numpy(map);
        }
        int __len__() { return $self->bursts.size(); }
        void append(PyObject* _in) {
          using namespace DCProgs;
          
          t_Burst burst = convert_pyburst(_in);
          if(burst.size() == 0) return;
          if(burst.size() % 2 != 1)  
            throw errors::PythonValueError("Each burst must have an odd number of intervals.");
          $self->bursts.emplace_back(std::move(burst));
        }
        PyObject* __str__() {
          std::ostringstream sstr; sstr << *($self);
          return PyString_FromString(sstr.str().c_str());
        }
        PyObject* __repr__() {
          std::ostringstream sstr; sstr << *($self);
          return PyString_FromString(sstr.str().c_str());
        }
      }
     
      //! Computes likelihood for each burst in separate value.
      DCProgs::t_rvector vector(DCProgs::t_rmatrix const &_Q) const;
      //! Computes likelihood for each burst in separate value.
      DCProgs::t_rvector vector(QMatrix const &_Q) const;
      //! Log-likelihood 
      DCProgs::t_real operator()(DCProgs::t_rmatrix const &_Q) const;
      //! Log-likelihood 
      DCProgs::t_real operator()(QMatrix const &_Q) const;
  };
}
