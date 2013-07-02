// More translating exception from c++ to python
%exception DeterminantEq::DeterminantEq(PyObject *_in, int _nopen, double _s, bool _doopen=true) {
  try { $action }
  DCPROGS_CATCH(SWIG_fail);
}
%exception DeterminantEq::DeterminantEq(StateMatrix *_in, double _s, bool _doopen=true) {
  try { $action }
  DCPROGS_CATCH(SWIG_fail);
}

%feature("autodoc", "Compute determinant W needed to approximate missed event G\n\n"
                    "There are several ways to instanciate this object:\n\n"
                    ">>> DeterminantEq(matrix, nopen, tau, doopen=True)\n\n"
                    "or,\n\n"
                    ">>> DeterminantEq(statematrix, tau, doopen=True)\n\n"
                    "Parameters:\n"
                    "   matrix: numpy array\n"
                    "      Transition matrix, where the upper left corner"
                         " contain open-open transitions\n"  
                    "   nopen: int\n"
                    "      Number of open states in the transition matrix.\n"
                    "   statematrix: StateMatrix\n"
                    "      State matrix created from a transition matrix "
                          "and the number of open states\n" 
                    "   tau: number\n"
                    "      Max length of missed events.\n"
                    "   doopen: Boolean\n"
                    "      Whether to compute AF (True, default) or AF determinant.\n")
         DeterminantEq;
//! Adds determinant functor
class DeterminantEq {

  public: 

  // Rename get_, set_tau so they are hidden (private) in python.
  %rename(_get_tau) get_tau() const; 
  %rename(_set_tau) set_tau(DCProgs::t_real); 

  // Extends the constructors so they can take numpy arrays as input.
  %extend {
    DeterminantEq(PyObject *_in, int _nopen, double _tau, bool _doopen=true) {
      if(_nopen <= 0)
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative or zero.");
      if(not PyArray_Check(_in))
        throw DCProgs::errors::PythonTypeError("Expected a numpy array on input.");
      DCProgs::t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix((PyArrayObject*)_in);
      if(matrix.rows() != matrix.cols())
        throw DCProgs::errors::PythonValueError("Expected a square matrix on input.");
      if(matrix.rows() == 0)
        throw DCProgs::errors::PythonValueError("Expected a non-empty square matrix on input.");
      if(_nopen >= matrix.rows())
        throw DCProgs::errors::PythonValueError("Number of closed states cannot be zero.");
      return new DCProgs::DeterminantEq( DCProgs::StateMatrix(std::move(matrix), _nopen),
                                         _tau, _doopen );
    }
    DeterminantEq(StateMatrix *_in, double _tau, bool _doopen=true) {
      if(_in->nopen <= 0) 
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative or zero.");
      if(_in->matrix.rows() != _in->matrix.cols()) 
        throw DCProgs::errors::PythonValueError("Expected a square state matrix on input.");
      if(_in->matrix.rows() == 0)
        throw DCProgs::errors::PythonValueError("Expected a non-empty square matrix on input.");
      if(_in->nopen >= _in->matrix.rows()) 
        throw DCProgs::errors::PythonValueError("Number of closed states cannot be zero.");
      return new DCProgs::DeterminantEq( *_in, _tau, _doopen );
    }
  }

  // Converts get_, set_tau to property.
  DCProgs::t_real get_tau() const;
  void set_tau(DCProgs::t_real);
  %pythoncode %{
    __swig_getmethods__["tau"] = _get_tau
    __swig_setmethods__["tau"] = _set_tau
    if _newclass: tau = property(_get_tau, _set_tau, doc="Max length of mixed events.")
  %}

  // Defines determinant function for a scalar
  DCProgs::t_real operator()(DCProgs::t_real _s);
  // Defines determinant function for a scalar + changing tau
  DCProgs::t_real operator()(DCProgs::t_real _s, DCProgs::t_real _tau);
  // The following two functions enhance the behavior of the function in python by allowing users to
  // input numpy arrays and get the determinants for each.
  %extend { 
    // Defines determinant function for a numpy array
    PyObject* operator()(PyObject * _s) {
      try { return apply_real(_s, *$self); }
      DCPROGS_CATCH(return NULL;);
    }
    // Defines determinant function for a numpy array + changing tau
    PyObject* operator()(PyObject * _s, DCProgs::t_real _tau) {
      try {
        DCProgs::t_real const oldtau = $self->get_tau();
        $self->set_tau(_tau);
        PyObject * result;
        try { result = apply_real(_s, *$self); }
        catch(...) {
          $self->set_tau(oldtau); 
          throw;
        }
        $self->set_tau(oldtau); 
        return result;
      } DCPROGS_CATCH(return NULL;);
    }
  }
  %pythoncode %{
    __call__.__doc__ = "Computes determinant W\n\n"                                    \
                       "Parameters:\n"                                                 \
                       "  s: scalar, tuple, list, array\n"                             \
                       "    The laplace scale.\n"                                      \
                       "  tau: *optional* number\n"                                    \
                       "    If present, it is the max length of missed events.\n\n"    \
                       "Returns: If a scalar, returns a scalar. "                      \
                                "Otherwise returns a numpy array."
  %}
                     
  // The following functions allows us to get H for scalar values.
  %extend { 
    // Defines determinant function for a numpy array
    PyObject* H(DCProgs::t_real _s) {
      try { return DCProgs::numpy::wrap_to_numpy($self->H(_s)); }
      DCPROGS_CATCH(return NULL;);
    }
    // Defines determinant function for a numpy array + changing tau
    PyObject* H(DCProgs::t_real _s, DCProgs::t_real _tau) {
      try { return DCProgs::numpy::wrap_to_numpy($self->H(_s)); }
      DCPROGS_CATCH(return NULL;);
    }
  }
  %pythoncode %{
    H.__doc__ = "Computes the matrix H\n\n"                                            \
                "H is defined as  :math:`"                                             \
                "Q_{AA} + Q_{AF}\ \int_0^\tau e^{-st}e^{Q_{FF}t}\partial\,t\ Q_{FA}`." \
                "\n\nParameters:\n"                                                    \
                "  s: scalar\n"                                                        \
                "    The laplace scale.\n"                                             \
                "  tau: *optional* number\n"                                           \
                "    If present, it is the max length of missed events.\n\n"           \
                "Returns: a numpy array."
  %}

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
