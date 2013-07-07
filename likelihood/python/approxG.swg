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
%feature("autodoc", "Computes approximate missed-event G.") ApproxG;
%feature("autodoc", "Open to close transitions") ApproxG::af;
%feature("autodoc", "Closed to open transitions") ApproxG::fa;
class ApproxG {
  public:
  
  // Other way to create this object.
  %extend {
    ApproxG(DCProgs::DeterminantEq const &_af, PyObject * _roots_af, 
            DCProgs::DeterminantEq const &_fa, PyObject * _roots_fa ) {
      std::vector<DCProgs::Root> roots_af = convert_to_roots(_roots_af);
      std::vector<DCProgs::Root> roots_fa = convert_to_roots(_roots_fa);
      return new DCProgs::ApproxG(_af, roots_af, _fa, roots_fa);
    }
  }

  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);
  DCProgs::t_rmatrix R_af(DCProgs::t_real t);
  DCProgs::t_rmatrix R_fa(DCProgs::t_real t);
};
%clear ApproxG;

%pythoncode %{
  def create_approxg(rates, tau):

    """ Creates a ApproxG matrix from knowledge of rate matrix. """
    determinant_af = DeterminantEq(rates, tau) 
    determinant_fa = DeterminantEq(rates.transpose(), tau) 
    roots_af = find_roots(determinant_af)
    roots_fa = find_roots(determinant_fa)
    return ApproxG(determinant_af, roots_af, determinant_fa, roots_fa)
%}
