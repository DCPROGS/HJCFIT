%{
  namespace {
    //! Special type for bursts
    typedef std::vector<DCProgs::t_rvector> t_Bursts;
    //! \brief Creates a pre-translated (to cpp) lists of bursts.
    //! \details This could make it a bit faster to call the likelihood when optimizing.
    //!          It should be used in conjunction with compute_bursts.
    t_Bursts* create_bursts(PyObject *_bursts) {
      using namespace DCProgs;
      if(not PySequence_Check(_bursts))
        throw errors::PythonTypeError("Expected a sequence on input.");
      t_Bursts* result(new t_Bursts);
      if(not result) errors::Runtime("Could not create burst.");
      Object<> py_iterator = steal_ref(PyObject_GetIter(_bursts));
      if(not py_iterator) throw errors::PythonErrorAlreadyThrown();
      while(Object<> item = steal_ref(PyIter_Next(~py_iterator))) {
        if(not PySequence_Check(~item))
          throw errors::PythonTypeError("Expected a sequence of sequences on input.");
        result->emplace_back(numpy::map_to_rvector(~item));
      }
  
      if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
      return result;
    }
  
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
          std::cout << "AM HERE" << std::endl;
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
    
    //! Computes likelihood of a time series in a pre-translated vector of bursts
    //! \param[in] _bursts: Vector of bursts.
    //! \param[in] _g: The likelihood functor. It should have an `af(t_real)` and an `fa(t_real)`
    //!                member function, where the argument is the length of an open or shut interval.
    //! \param[in] _initial: initial occupancies.
    //! \param[in] _final: final occupancies.
    template<class T_G> 
      PyObject* py_compute_bursts(T_G const & _g, t_Bursts *_bursts, 
                                  PyObject *_initial, PyObject *_final) {
        using namespace DCProgs;
        t_initvec const initial = numpy::map_to_initvec(_initial).array();
        t_rvector const final   = numpy::map_to_rvector(_final).array();
        
        npy_intp dims[1] = { static_cast<npy_intp>(_bursts->size()) };
        Object<> result = steal_ref(PyArray_SimpleNew(1, dims, numpy::type<t_real>::value));
        if(not result) throw errors::PythonErrorAlreadyThrown();
        t_real *i_result
            = static_cast<t_real*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(~result)));
        for(t_Bursts::const_reference burst: *_bursts) {
          *i_result = chained_likelihood(_g, &burst(0), &burst(0) + burst.size(), initial, final); 
          ++i_result;
        }
        return result.release();
      }
    PyObject* compute_bursts(DCProgs::IdealG *_idealg, t_Bursts* _bursts, 
                             PyObject* _initial, PyObject* _final) {
      return py_compute_bursts(*_idealg, _bursts, _initial, _final);
    } 
    PyObject* compute_bursts(DCProgs::MissedEventsG *_idealg, t_Bursts* _bursts, 
                             PyObject* _initial, PyObject* _final) {
      return py_compute_bursts(*_idealg, _bursts, _initial, _final);
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

%feature("autodoc", "Computes likelihood of an opaque of bursts.\n\n"
         "A slightly faster form of computing the likelihood over a list of bursts. "
         "The list is pre-translated to a c++ object over which we can loop in c++, "
         "without going back to python. We may then be able to apply parallelization " 
         "methods. Without parallelization, it does not seem to be better than 10\% faster.\n\n"
         ":param firstArgument: \n"
         "   Can be either an :class:`IdealG` or and :class:`MissedEventsG` instance.\n"
         ":param bursts:\n"
         "   This must be the opaque object returned by :meth:`create_bursts`. "
            "It contains the bursts in c++ form.\n"
         ":param initial:\n"
         "   The initial occupancy vector.\n" 
         ":param final:\n"
         "   The final occupancy vector.\n")
        compute_bursts(DCProgs::IdealG *_idealg, t_Bursts* _bursts, 
                       PyObject* _initial, PyObject* _final);
PyObject* compute_bursts(DCProgs::IdealG *_idealg, t_Bursts* _bursts, 
                         PyObject* _initial, PyObject* _final);
PyObject* compute_bursts(DCProgs::MissedEventsG *_idealg, t_Bursts* _bursts, 
                         PyObject* _initial, PyObject* _final);


%newobject create_bursts(PyObject*);
t_Bursts* create_bursts(PyObject *_bursts);
