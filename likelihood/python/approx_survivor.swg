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
