%{
  // Helper function to convert output of root finding objects to cpp
  std::vector<DCProgs::Root> convert_to_roots(PyObject *_in) {
    Object<> iterator = steal_ref(PyObject_GetIter(_in)); 
    if(not iterator) throw DCProgs::errors::PythonErrorAlreadyThrown();

    std::vector<DCProgs::Root> result;
    while(Object<> item = steal_ref(PyIter_Next(~iterator))) {

      double root;
      int multiplicity;

      if(not PyArg_ParseTuple(~item, "di", &root, &multiplicity))
        throw DCProgs::errors::PythonErrorAlreadyThrown();
      result.emplace_back(DCProgs::t_real(root), DCProgs::t_int(multiplicity));
    }
    if(PyErr_Occurred() != NULL) throw DCProgs::errors::PythonErrorAlreadyThrown();
    
    return result;
  }
%}
%feature("autodoc", "Computes approximate missed-event survivor function.") ApproxSurvivor;
%feature("autodoc", "Open to close transitions") ApproxSurvivor::af;
%feature("autodoc", "Closed to open transitions") ApproxSurvivor::fa;
class ApproxSurvivor {
  public:
  
  // Other way to create this object.
  %extend {
    ApproxSurvivor(DCProgs::DeterminantEq const &_af, PyObject * _roots_af, 
                   DCProgs::DeterminantEq const &_fa, PyObject * _roots_fa ) {
      std::vector<DCProgs::Root> roots_af = convert_to_roots(_roots_af);
      std::vector<DCProgs::Root> roots_fa = convert_to_roots(_roots_fa);
      return new DCProgs::ApproxSurvivor(_af, roots_af, _fa, roots_fa);
    }
  }

  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);

  %extend {
    PyObject* const af_components;
    PyObject* const fa_components;
  }
  %{ namespace { 
    template<class T> 
      PyObject* affa_components_get(DCProgs::ApproxSurvivor *_in,
                                    Py_ssize_t const _N, T const &_getstuff) {
        Object<> result = steal_ref( PyList_New(_N) );
        if(not result) throw DCProgs::errors::PythonErrorAlreadyThrown();
        for(Py_ssize_t i(0); i < _N; ++i) {
          DCProgs::Asymptotes::t_MatrixAndRoot const &item = _getstuff(i);
          Object<> py_item = steal_ref( PyTuple_New(2) );
          if(not py_item) throw DCProgs::errors::PythonErrorAlreadyThrown();
          Object<> py_root = steal_ref(PyFloat_FromDouble(std::get<1>(item)));
          if(not py_root) throw DCProgs::errors::PythonErrorAlreadyThrown();
          DCProgs::t_rmatrix const & matrix = std::get<0>(item);
          Object<> py_matrix = steal_ref(DCProgs::numpy::wrap_to_numpy(matrix));
          if(not py_matrix) throw DCProgs::errors::PythonErrorAlreadyThrown();
          PyArray_CLEARFLAGS((PyArrayObject*) ~py_matrix, NPY_ARRAY_WRITEABLE);
          PyTuple_SET_ITEM(~py_item, 0, py_matrix.release());
          PyTuple_SET_ITEM(~py_item, 1, py_root.release());
          PyList_SET_ITEM(~result, i, py_item.release());
        }
        return result.release();
      }
    PyObject* DCProgs_ApproxSurvivor_af_components_get(DCProgs::ApproxSurvivor *_in) {
        return affa_components_get( _in, (Py_ssize_t)_in->nb_af_components(),
                                    [&_in](DCProgs::t_int i)
                                      { return _in->get_af_components(i); } );
    }
    PyObject* DCProgs_ApproxSurvivor_fa_components_get(DCProgs::ApproxSurvivor *_in) {
        return affa_components_get( _in, (Py_ssize_t)_in->nb_fa_components(),
                                    [&_in](DCProgs::t_int i)
                                       { return _in->get_fa_components(i); } );
    }
 
  } %}

  // Add some printout functions.
  %extend {
    PyObject* __str__() {
      std::ostringstream sstr; sstr << *($self);
      return PyString_FromString(sstr.str().c_str());
    }
    PyObject* __repr__() {
      std::ostringstream sstr; sstr << *($self);
      return PyString_FromString(sstr.str().c_str());
    }
  }
};
%clear ApproxSurvivor;

%pythoncode %{
  def create_approx_survivor(rates, tau):

    """ Creates a ApproxSurvivor function from knowledge of rate matrix. """
    determinant_af = DeterminantEq(rates, tau) 
    determinant_fa = DeterminantEq(rates.transpose(), tau) 
    roots_af = find_roots(determinant_af)
    roots_fa = find_roots(determinant_fa)
    return ApproxSurvivor(determinant_af, roots_af, determinant_fa, roots_fa)
%}
