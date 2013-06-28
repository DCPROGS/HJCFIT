// More translating exception from c++ to python
%exception StateMatrix::StateMatrix(PyObject *_in, int _nopen) {
  try { $action }
  DCPROGS_CATCH(SWIG_fail);
}

// Now we define the bindings for the StateMatrix
//! \brief State matrix that can  be partitioned into open/shut states.
//! \details In practice, this is a two tuple with some helper functions to get corners.
class StateMatrix {
  public:


  //! Number of open states.
  DCProgs::t_int nopen; 
  %typemap(in) DCProgs::t_rmatrix matrix { 
    try { $1 = DCProgs::numpy::map_to_rmatrix((PyArrayObject*)$input); }
    DCPROGS_CATCH(SWIG_fail);
  }
  %typemap(out) DCProgs::t_rmatrix { 
    try { $result = DCProgs::numpy::wrap_to_numpy(arg1->matrix, $self); }
    DCPROGS_CATCH(SWIG_fail);
  }
  //! The matrix itself.
  DCProgs::t_rmatrix matrix; 

  %clear DCProgs::t_rmatrix;

  StateMatrix();

  // Extends the constructors so they can take numpy arrays as input.
  %extend {
    StateMatrix(PyObject *_in, int _nopen) {
      if(_nopen < 0)
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative.");
      if(not PyArray_Check(_in))
        throw DCProgs::errors::PythonTypeError("Expected a numpy array on input.");
      DCProgs::t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix((PyArrayObject*)_in);
      if(_nopen > std::max(matrix.rows(), matrix.cols()) )
        throw DCProgs::errors::PythonValueError(
                "Number of open states cannot be larger than the number states.");
      return new DCProgs::StateMatrix(std::move(matrix), _nopen);
    }
  }

  // Add some printout functions.
  %extend {
    PyObject* __str__() {
      if($self->matrix.rows() == 0 or $self->matrix.cols() == 0) 
        return PyString_FromString("Un-initialized transition matrix.");
      std::ostringstream sstr;
      sstr << "Transition matrix with " << $self->nopen
           << " open states:\n" << DCProgs::numpy_io($self->matrix) << "\n";
      return PyString_FromString(sstr.str().c_str());
    }
    PyObject* __repr__() {
      if($self->matrix.rows() == 0 or $self->matrix.cols() == 0)
        return PyString_FromString("StateMatrix()");
      std::ostringstream sstr;
      sstr << "DeterminantEq( array(" << DCProgs::numpy_io($self->matrix) << "), "
           << $self->nopen << ")";
      return PyString_FromString(sstr.str().c_str());
    }
  }
};
// Translate the convenience functions aa/af/fa/ff to properties in python. 
// This means we have a natural access to the matrix blocks as "matrix.af", without parenthesis. 
%pythoncode %{
  StateMatrix.aa = property(lambda self: self.matrix[:self.nopen, :self.nopen],
                            doc=""" Open to open transitions. """)
  StateMatrix.af = property(lambda self: self.matrix[:self.nopen, self.nopen:],
                            doc=""" Open to close transitions. """)
  StateMatrix.fa = property(lambda self: self.matrix[self.nopen:, :self.nopen],
                            doc=""" Open to close transitions. """)
  StateMatrix.ff = property(lambda self: self.matrix[self.nopen:, self.nopen:],
                            doc=""" Open to close transitions. """)
%}

