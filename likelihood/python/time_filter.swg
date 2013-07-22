%feature("autodoc", "Filters out time series for intervals smaller that :math:`\\tau`\n\n"
         "This procedure erases small intervals from the time series.\n\n"
         ":param times:\n"
         "   A timer series to filter. It can also arrays of time series. In that case, "
             "The filters are applied independently to all but the innermost sequence.\n\n"
             ":returns: A numpy array, or a list of numpy arrays.")
         time_filter;
%feature("autodoc", "Filters out  intervals smaller that :math:`\\tau`\n\n"
         "This procedure erases small intervals from the time series that can be created "
         "from the input. \n\n"
         ":param times:\n"
         "   A timer series to filter. It can also arrays of time series. In that case, "
             "The filters are applied independently to all but the innermost sequence.\n\n"
             ":returns: A numpy array, or a list of numpy arrays.") 
        interval_filter;


%{
  namespace {

    //! Calls either time_filter or interval filter.
    template<class T_FUNC>
      PyObject* filters_impl_(PyObject* _input, DCProgs::t_real _tau, T_FUNC const &_func) {
        using namespace DCProgs;
        if(not PySequence_Check(_input))
          throw errors::PythonTypeError("Expected a sequence on input.");

        Py_ssize_t const N = PySequence_Size(_input);
        if(N == 0) {
          npy_intp dims[1] = {0};
          return PyArray_SimpleNew(1, dims, numpy::type<t_real>::value);
        }

        Object<> first_item = steal_ref(PySequence_GetItem(_input, 0));
        if(PySequence_Check(~first_item)) {
          Object<> result = steal_ref(PyList_New(N));
          for(Py_ssize_t i(0); i < N; ++i) {
            Object<> item_in = steal_ref(PySequence_GetItem(_input, i));
            PyList_SET_ITEM(~result, i, filters_impl_(~item_in, _tau, _func)); 
          }
          return result.release();
        }

        return numpy::wrap_to_numpy(_func(numpy::map_to_rmatrix(_input), _tau)); 
      }

    //! Bindings for time_filter.
    PyObject* time_filter(PyObject* _input, DCProgs::t_real _tau) {
      return filters_impl_(_input, _tau, DCProgs::time_filter);
    }
    //! Bindings for interval_filter.
    PyObject* interval_filter(PyObject* _input, DCProgs::t_real _tau) {
      return filters_impl_(_input, _tau, DCProgs::interval_filter);
    }

  }
%}

PyObject* time_filter(PyObject* _input, DCProgs::t_real _tau);
PyObject* interval_filter(PyObject* _input, DCProgs::t_real _tau);
