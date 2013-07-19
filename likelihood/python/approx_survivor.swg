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
%feature("autodoc", "Computes approximate missed-event survivor function.\n\n"
         "This object can be created  explicitely from the determinantal equations and the "
         "roots of the AF and FA blocks:\n\n"
         ">>> approx = ApproxSurvivor(determinant_af, roots_af, determinant_fa, roots_fa)\n\n"
         "However, it can also be created using the helper function "
         ":meth:`create_approx_survivor`.\n\n"
         ":param DeterminantEq determinant_af:\n"
         "  The determinantal equation for the af block.\n"
         ":param roots_af:\n"
         "  The roots of the af determinantal equation. The should come in the format "
            "`[(root, multiplicity), (root, multiplicity), ...]`.\n" 
         ":param DeterminantEq determinant_fa:\n"
         "  The determinantal equation for the fa block. It should likely be "
           "the transpose of determinant_fa.\n"
         ":param roots_fa:\n"
         "  The roots of the fa determinantal equation. The should come in the format "
            "`[(root, multiplicity), (root, multiplicity), ...]`.\n" 
         ":param float tau:\n"
         "  The maximum length of the missed events.\n") ApproxSurvivor;
%feature("autodoc", "Open to close transitions") ApproxSurvivor::af;
%feature("autodoc", "Closed to open transitions") ApproxSurvivor::fa;
class ApproxSurvivor {
  public:
  
  %extend {
    // Creates object from determinants and roots.
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
                                    [&_in](Py_ssize_t i) {
                 return _in->get_af_components(static_cast<DCProgs::t_int>(i)); } );
    }
    PyObject* DCProgs_ApproxSurvivor_fa_components_get(DCProgs::ApproxSurvivor *_in) {
        return affa_components_get( _in, (Py_ssize_t)_in->nb_fa_components(),
                                    [&_in](Py_ssize_t i) {
                 return _in->get_fa_components(static_cast<DCProgs::t_int>(i)); } );
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

