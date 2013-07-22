%feature("autodoc", "Filters out time series for intervals smaller that :math:`\\tau`\n\n"
         ":param times:\n"
         "   A timer series to filter. It can also arrays of time series. In that case, "
             "The filters are applied independently to all but the innermost sequence.\n\n"
             ":returns: A numpy array, or a list of numpy arrays.")
         time_filter;


%{
  PyObject* time_filter(PyObject* _input, DCProgs::t_real _tau) {
    using namespace DCProgs;
    typedef ::Eigen::Map<t_rvector, 0, ::Eigen::InnerStride<>> t_Map;
    if(not PySequence_Check(_input))
      throw errors::PythonTypeError("Expected a sequence on input.");
    Py_ssize_t const N = PySequence_Size(_input);
    if(N == 0) return numpy::wrap_to_numpy(t_rvector::Zero(0));
    Object<> first_item = steal_ref(PySequence_GetItem(_input, 0));
    if(PySequence_Check(~first_item)) {
      Object<> result = steal_ref(PyList_New(N));
      for(Py_ssize_t i(0); i < N; ++i) {
        Object<> item_in = steal_ref(PySequence_GetItem(_input, i));
        Object<> item_out = steal_ref(time_filter(~item_in, _tau)); 
        PyList_SET_ITEM(~result, i, item_out.release());
      }
      return ~result;
    }
    return numpy::wrap_to_numpy(time_filter(numpy::map_to_rmatrix(_input), _tau)); 
  }
%}
PyObject* time_filter(PyObject* _input, DCProgs::t_real _tau);
