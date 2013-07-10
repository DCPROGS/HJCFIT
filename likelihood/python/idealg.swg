%feature("autodoc", "Ideal Likelihood.\n\n"
         "This object can be instantiated one of several way:\n\n"
         "   - With a matrix and an integer\n\n"
         "     >>> idealg = IdealG(array([...]), 2)\n\n"
         "   - With a StateMatrix\n\n"
         "     >>> matrix = StateMatrix(array([...]), 2)\n"
         "     >>> idealg = IdealG(matrix)\n\n") IdealG;
%feature("autodoc", "AF transitions with respect to time: :math:`e^{tQ_{FF}}Q_{FA}`.") IdealG::af;
%feature("autodoc", "FA transitions with respect to time: :math:`e^{tQ_{AA}}Q_{AF}$.") IdealG::fa;
%feature("autodoc", "AF transitions with respect to scale: :math:`(sI - Q_{AA})^{-1}Q_{AF}`.")
        IdealG::laplace_af;
%feature("autodoc", "FA transitions with respect to scale: :math:`(sI - Q_{FF})^{-1}Q_{FA}$`.")
        IdealG::laplace_fa;
// Now we define the bindings for the StateMatrix
// \brief State matrix that can  be partitioned into open/shut states.
// \details In practice, this is a two tuple with some helper functions to get corners.
class IdealG {
  public:

  // Some constructors.
  %extend { 
    IdealG(DCProgs::StateMatrix const &states) {
      return new DCProgs::IdealG(states.matrix, states.nopen); 
    }
    IdealG(PyObject* _py_matrix, int _nopen) {
      if(_nopen < 0)
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative.");
      if(not PyArray_Check(_py_matrix))
        throw DCProgs::errors::PythonTypeError("Expected a numpy array on input.");
      DCProgs::t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix((PyArrayObject*)_py_matrix);
      if(_nopen > std::max(matrix.rows(), matrix.cols()) )
        throw DCProgs::errors::PythonValueError(
                "Number of open states cannot be larger than the number states.");
      return new DCProgs::IdealG(std::move(matrix), _nopen); 
    }
  }

  %extend {
    DCProgs::t_initvec _occupancies_initial() { return DCProgs::occupancies(*$self); }
    DCProgs::t_initvec _occupancies_final() { return DCProgs::occupancies(*$self, false); }
  }
  %pythoncode %{
    __swig_getmethods__["occupancies_initial"] = _occupancies_initial
    __swig_getmethods__["occupancies_final"] = _occupancies_final
    if _newclass:
      occupancies_initialf = property(_occupancies_initial, doc="Occupancies for initial states.")
      occupancies_fa = property(_occupancies_final, doc="Occupancies for final states.")
  %}

  DCProgs::t_rmatrix af(DCProgs::t_real t) const;
  DCProgs::t_rmatrix fa(DCProgs::t_real t) const;
  DCProgs::t_rmatrix laplace_fa(DCProgs::t_real s) const;
  DCProgs::t_rmatrix laplace_af(DCProgs::t_real s) const;
};
