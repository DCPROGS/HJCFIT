%feature("autodoc", "Ideal Likelihood.\n\n"
         "This object can be instantiated one of several way:\n\n"
         "   - With a matrix and an integer\n\n"
         "     >>> idealg = IdealG(array([...]), 2)\n\n"
         "   - With a QMatrix\n\n"
         "     >>> matrix = QMatrix(array([...]), 2)\n"
         "     >>> idealg = IdealG(matrix)\n\n") IdealG;
%feature("autodoc", "AF transitions with respect to time: :math:`e^{tQ_{FF}}Q_{FA}`.") IdealG::af;
%feature("autodoc", "FA transitions with respect to time: :math:`e^{tQ_{AA}}Q_{AF}$.") IdealG::fa;
%feature("autodoc", "AF transitions with respect to scale: :math:`(sI - Q_{AA})^{-1}Q_{AF}`.")
        IdealG::laplace_af;
%feature("autodoc", "FA transitions with respect to scale: :math:`(sI - Q_{FF})^{-1}Q_{FA}$`.")
        IdealG::laplace_fa;
// Now we define the bindings for the QMatrix
// \brief State matrix that can  be partitioned into open/shut states.
// \details In practice, this is a two tuple with some helper functions to get corners.
class IdealG {
  public:

  // Some constructors.
  %extend { 
    IdealG(DCProgs::QMatrix const &states) {
      return new DCProgs::IdealG(states.matrix, states.nopen); 
    }
    IdealG(PyObject* _py_matrix, int _nopen) {
      if(_nopen < 0)
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative.");
      DCProgs::t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix(_py_matrix);
      if(_nopen > std::max(matrix.rows(), matrix.cols()) )
        throw DCProgs::errors::PythonValueError(
                "Number of open states cannot be larger than the number states.");
      return new DCProgs::IdealG(std::move(matrix), _nopen); 
    }
  }

  %extend {
    DCProgs::t_initvec _initial_occupancies() { return DCProgs::occupancies(*$self); }
    DCProgs::t_initvec _final_occupancies() { return DCProgs::occupancies(*$self, false); }
  }
  %pythoncode %{
    __swig_getmethods__["initial_occupancies"] = _initial_occupancies
    __swig_getmethods__["final_occupancies"]   = _final_occupancies
    if _newclass:
      initial_occupancies = property(_initial_occupancies, doc="Occupancies for initial states.")
      final_occupancies = property(_final_occupancies, doc="Occupancies for final states.")
  %}

  DCProgs::t_rmatrix af(DCProgs::t_real t) const;
  DCProgs::t_rmatrix fa(DCProgs::t_real t) const;
  DCProgs::t_rmatrix laplace_fa(DCProgs::t_real s) const;
  DCProgs::t_rmatrix laplace_af(DCProgs::t_real s) const;

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
