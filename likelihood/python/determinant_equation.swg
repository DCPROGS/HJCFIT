/***********************
    DCProgs computes missed-events likelihood as described in
    Hawkes, Jalali and Colquhoun (1990, 1992)

    Copyright (C) 2013  University College London

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
************************/

%feature("autodoc", "Computes the matrix H\n\n" 
         "H is defined as  :math:`"
         "Q_{AA} + Q_{AF}\ \\int_0^\\tau e^{-st}e^{Q_{FF}t}\\partial\,t\ Q_{FA}`."
         ":param number s:\n"
         "    The laplace scale.\n"
         ":param number tau:\n"
         "    *Optional*. If present, it is the max length of missed events.\n\n"
         ":returns: a numpy array.") DeterminantEq::H;
%feature("autodoc", "Computes the derivative of H versus s\n\n" 
         "H is defined as  :math:`"
         "Q_{AA} + Q_{AF}\ \\int_0^\\tau e^{-st}e^{Q_{FF}t}\\partial\,t\ Q_{FA}`."
         ":param number s:\n"
         "    The laplace scale.\n"
         ":param number tau:\n"
         "    *Optional*. If present, it is the max length of missed events.\n\n"
         ":returns: a numpy array.") DeterminantEq::s_derivative(DeterminantEq const *, DCProgs::t_real s);
%feature("autodoc", "Compute determinant W needed to approximate missed event G\n\n"
                    "There are several ways to instanciate this object:\n\n"
                    ">>> DeterminantEq(matrix, nopen, tau)\n\n"
                    "or,\n\n"
                    ">>> DeterminantEq(qmatrix, tau)\n\n"
                    "Parameters:\n"
                    "   matrix: numpy array\n"
                    "      Transition matrix, where the upper left corner"
                         " contain open-open transitions\n"  
                    "   nopen: int\n"
                    "      Number of open states in the transition matrix.\n"
                    "   qmatrix: QMatrix\n"
                    "      State matrix created from a transition matrix "
                          "and the number of open states\n" 
                    "   tau: number\n"
                    "      Max length of missed events.\n")
         DeterminantEq;
//! Adds determinant functor
class DeterminantEq {

  public: 

  // Extends the constructors so they can take numpy arrays as input.
  %extend {
    DeterminantEq(PyObject *_in, int _nopen, double _tau) {
      if(_nopen <= 0)
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative or zero.");
      DCProgs::t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix(_in);
      if(matrix.rows() != matrix.cols())
        throw DCProgs::errors::PythonValueError("Expected a square matrix on input.");
      if(matrix.rows() == 0)
        throw DCProgs::errors::PythonValueError("Expected a non-empty square matrix on input.");
      if(_nopen >= matrix.rows())
        throw DCProgs::errors::PythonValueError("Number of closed states cannot be zero.");
      return new DCProgs::DeterminantEq( DCProgs::QMatrix(std::move(matrix), _nopen), _tau );
    }
    DeterminantEq(QMatrix *_in, double _tau) {
      if(_in->nopen <= 0) 
        throw DCProgs::errors::PythonValueError("Number of open states cannot be negative or zero.");
      if(_in->matrix.rows() != _in->matrix.cols()) 
        throw DCProgs::errors::PythonValueError("Expected a square state matrix on input.");
      if(_in->matrix.rows() == 0)
        throw DCProgs::errors::PythonValueError("Expected a non-empty square matrix on input.");
      if(_in->nopen >= _in->matrix.rows()) 
        throw DCProgs::errors::PythonValueError("Number of closed states cannot be zero.");
      return new DCProgs::DeterminantEq( *_in, _tau );
    }
  }

  DCProgs::DeterminantEq transpose();
  // Rename get_, set_tau so they are hidden (private) in python.
  %rename(_get_tau) get_tau() const; 
  %rename(_set_tau) set_tau(DCProgs::t_real); 
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
  // The following functions enhance the behavior of the function in python by allowing users to
  // input numpy arrays and get the determinants for each value. 
  // Unfortunately, It doesn't seem to work with overload for tau.
  %extend { 
    // Defines determinant function for a numpy array
    PyObject* operator()(PyObject * _s) { return apply_numpy(_s, *$self); }
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
                     
  DCProgs::t_rmatrix H(DCProgs::t_real s);
  DCProgs::t_rmatrix H(DCProgs::t_real s, DCProgs::t_real tau);
  DCProgs::t_rmatrix s_derivative(DCProgs::t_real s);
  DCProgs::t_rmatrix s_derivative(DCProgs::t_real s, DCProgs::t_real tau);
  // Add the ability to pass in numpy arrays.
  // Unfortunately, no tau overload.
  %extend {
    PyObject* H(PyObject *_s) {
      return apply_numpy(_s, [$self](DCProgs::t_real s) { return $self->H(s); }); 
    }
    PyObject* s_derivative(PyObject *_s) {
      return apply_numpy(_s, [$self](DCProgs::t_real s) { return $self->s_derivative(s); }); 
    }
  }

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

// adds some conversion function for PyArg_Parse... functions
%{
  namespace {
    Py_ssize_t convert_determinant(PyObject *_object, void* _result) {
      try {
        int const tryconv = SWIG_ConvertPtr( _object, &_result, 
                                             SWIGTYPE_p_DCProgs__DeterminantEq, 0 );
        if (!SWIG_IsOK(tryconv)) {
          PyErr_SetString(PyExc_TypeError, "Input is not a DeterminantEq.");
          return 0;
        }
      } DCPROGS_CATCH(return 0;);
      return 1;
    }
  }
%}
