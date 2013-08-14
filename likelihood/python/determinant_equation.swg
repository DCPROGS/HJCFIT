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
         "\\mathcal{Q}_{AA} + \\mathcal{Q}_{AF}\ "
         "\\int_0^\\\\tau e^{-st}e^{\\mathcal{Q}_{FF}t}\\partial\,t\ \\mathcal{Q}_{FA}`.\n\n"
         ":param number s:\n"
         "    The laplace scale. It can be a scalar or a numpy array of any shape.\n"
         "    In the latter case, the return is a array of the same shape where\n"
         "    each element is a matrix corresponding to the element in the input\n"
         "    array.\n"
         ":param number tau:\n"
         "    *Optional*. If present, it is the max length of missed events.\n\n"
         ":returns: a numpy array.\n")  DeterminantEq::H(PyObject*);
%feature("autodoc", "Computes the derivative of H versus s\n\n" 
         "H is defined as  :math:`"
         "\\mathcal{Q}_{AA} + \\mathcal{Q}_{AF}\ \\int_0^\\\\tau "
         "e^{-st}e^{\\mathcal{Q}_{FF}t}\\partial\,t\ \\mathcal{Q}_{FA}`.\n\n"
         ":param number s:\n"
         "    The laplace scale. It can be a scalar or a numpy array of any shape.\n"
         "    In the latter case, the return is a array of the same shape where\n"
         "    each element is a matrix corresponding to the element in the input\n"
         "    array.\n"
         ":param number tau:\n"
         "    *Optional*. If present, it is the max length of missed events.\n\n"
         ":returns: a numpy array.\n") DeterminantEq::s_derivative(PyObject*);
%feature("autodoc", "Compute determinant W needed to approximate missed event G\n\n"
         "This object can be instantiated from a square matrix, the number of open states, and the\n"
         "resolution, or maximum length of missed events :math:`\\\\tau`:\n\n"
         ">>> DeterminantEq(matrix, nopen, tau)\n\n"
         "or, it can be instantiated from a :py:class:`~dcprogs.likelihood.QMatrix` instance \n"
         "and :math:`\\\\tau`\n\n"
         ">>> DeterminantEq(qmatrix, tau)\n\n"
         ":param matrix:\n"
         "      Transition rate matrix, where the upper left corner contain open-open transitions\n"  
         ":param integer nopen:\n"
         "      Number of open states in the transition matrix.\n"
         ":param qmatrix:\n"
         "      :py:class:`QMatrix` instance\n" 
         ":param number tau:\n"
         "      Max length of missed events.\n")
         DeterminantEq;
%feature("autodoc", "Computes the determinant of :math:`W(s)`.\n\n"
         ":param number s:\n"
         "    The laplace scale. It can be a scalar or a numpy array of any shape.\n"
         "    In the latter case, the return is a array of the same shape where\n"
         "    each element is a matrix corresponding to the element in the input\n"
         "    array.\n"
         ":param number tau:\n"
         "    *Optional*. If present, it is the max length of missed events.\n\n"
         ":returns: a numpy array.\n") DeterminantEq::operator()(PyObject*);


// Adds determinant functor
class DeterminantEq {

  public: 

  // Extends the constructors so they can take numpy arrays as input.
  %extend {
    DeterminantEq(PyObject *_in, DCProgs::t_uint _nopen, DCProgs::t_real _tau) {
      using namespace DCProgs;
      if(_nopen == 0)
        throw errors::PythonValueError("Number of open states cannot be negative or zero.");
      t_rmatrix const matrix = DCProgs::numpy::map_to_rmatrix(_in);
      if(matrix.rows() != matrix.cols())
        throw errors::PythonValueError("Expected a square matrix on input.");
      if(matrix.rows() == 0)
        throw errors::PythonValueError("Expected a non-empty square matrix on input.");
      if(_nopen >= static_cast<t_uint>(matrix.rows()))
        throw errors::PythonValueError("Number of closed states cannot be zero.");
      return new DeterminantEq( DCProgs::QMatrix(std::move(matrix), _nopen), _tau );
    }
    DeterminantEq(QMatrix *_in, DCProgs::t_real _tau) {
      using namespace DCProgs;
      if(_in->nopen == 0) 
        throw errors::PythonValueError("Number of open states cannot be negative or zero.");
      if(_in->matrix.rows() != _in->matrix.cols()) 
        throw errors::PythonValueError("Expected a square state matrix on input.");
      if(_in->matrix.rows() == 0)
        throw errors::PythonValueError("Expected a non-empty square matrix on input.");
      if(_in->nopen >= static_cast<t_uint>(_in->matrix.rows())) 
        throw errors::PythonValueError("Number of closed states cannot be zero.");
      return new DeterminantEq( *_in, _tau );
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
                     
  // Add the ability to pass in numpy arrays.
  // Unfortunately, no tau overload.
  DCProgs::t_rmatrix H(DCProgs::t_real s);
  DCProgs::t_rmatrix H(DCProgs::t_real s, DCProgs::t_real tau);
  DCProgs::t_rmatrix s_derivative(DCProgs::t_real s);
  DCProgs::t_rmatrix s_derivative(DCProgs::t_real s, DCProgs::t_real tau);
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
