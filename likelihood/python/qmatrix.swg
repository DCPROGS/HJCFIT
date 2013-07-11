%feature("autodoc", "Transition rate matrix.\n\n"
         "Open-open transitions should be located in the top-left corner.") QMatrix;
%feature("autodoc", "Swaps AA block with FF, AF block with FA, nopen with nclose.\n\n")
         QMatrix::transpose;
%feature("autodoc", "Transition matrix itself.") QMatrix::matrix;
%feature("autodoc", "Number of open states.") QMatrix::nopen;
// Now we define the bindings for the QMatrix
// \brief State matrix that can  be partitioned into open/shut states.
// \details In practice, this is a two tuple with some helper functions to get corners.
class QMatrix {
  public:


  // Number of open states.
  DCProgs::t_int nopen; 
  // Tells SWIG how to convert eigen matrices to numpy matrices.
  %typemap(in) DCProgs::t_rmatrix matrix { 
    try { $1 = DCProgs::numpy::map_to_rmatrix($input); }
    DCPROGS_CATCH(SWIG_fail);
  }
  // The $self argument below is quite important. It tells SWIG that the memory for the matrix is
  // owned by c++, and that the c++ object should stay alive as long as the matrix lives in python.
  %typemap(out) DCProgs::t_rmatrix matrix { 
    try { $result = DCProgs::numpy::wrap_to_numpy(arg1->matrix, $self); }
    DCPROGS_CATCH(SWIG_fail);
  }
  //! The matrix itself.
  DCProgs::t_rmatrix matrix; 

  // This is the end of automatic conversion from eigen to numpy.
  // It is removed because elsewhere we don't know who will own the memory to the eigen matrix.
  %clear DCProgs::t_rmatrix matrix;

  QMatrix();

  // Extends the constructors so they can take numpy arrays as input.
  %extend {
    QMatrix(PyObject *_py_matrix, int _nopen) {
      if(_nopen < 0)
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative.");
      DCProgs::t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix(_py_matrix);
      if(_nopen > std::max(matrix.rows(), matrix.cols()) )
        throw DCProgs::errors::PythonValueError(
                "Number of open states cannot be larger than the number states.");
      return new DCProgs::QMatrix(std::move(matrix), _nopen);
    }
  }

  // Add some printout functions.
  %extend {
    PyObject* __str__() {
      if($self->matrix.rows() == 0 or $self->matrix.cols() == 0) 
        return PyString_FromString("Un-initialized transition matrix.");
      std::ostringstream sstr; sstr << *($self);
      return PyString_FromString(sstr.str().c_str());
    }
    PyObject* __repr__() {
      if($self->matrix.rows() == 0 or $self->matrix.cols() == 0)
        return PyString_FromString("QMatrix()");
      std::ostringstream sstr;
      sstr << "QMatrix( array(" << DCProgs::numpy_io($self->matrix) << "), "
           << $self->nopen << ")";
      return PyString_FromString(sstr.str().c_str());
    }
  }

  QMatrix transpose() const;
};
// Translate the convenience functions aa/af/fa/ff to properties in python. 
// This means we have a natural access to the matrix blocks as "matrix.af", without parenthesis. 
%pythoncode %{
  QMatrix.aa = property(lambda self: self.matrix[:self.nopen, :self.nopen],
                        doc=""" Open to open transitions. """)
  QMatrix.af = property(lambda self: self.matrix[:self.nopen, self.nopen:],
                        doc=""" Open to close transitions. """)
  QMatrix.fa = property(lambda self: self.matrix[self.nopen:, :self.nopen],
                        doc=""" Open to close transitions. """)
  QMatrix.ff = property(lambda self: self.matrix[self.nopen:, self.nopen:],
                        doc=""" Open to close transitions. """)
  QMatrix.nshut = property(lambda self: self.matrix.shape[0] - self.nopen,
                           lambda self, value: \
                               setattr(self, "nopen", self.matrix.shape[0] - value),
                           doc=""" Number of shut states. """);
%}

