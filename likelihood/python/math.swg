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
// Some functions to convert back and forth between python and t_real.
%{ namespace {
  // Converts input to real
  DCProgs::t_real convert_to_real(PyObject* _in) {
    using namespace DCProgs;
    if(PyLong_Check(_in)) return static_cast<t_real>(PyLong_AsLong(_in));
#   ifndef DCPROGS_PYTHON3
      else if(PyInt_Check(_in)) return static_cast<t_real>(PyInt_AsLong(_in));
#   endif
    else if(PyFloat_Check(_in)) return static_cast<t_real>(PyFloat_AsDouble(_in));
    else if(PyArray_CheckScalar(_in))  {
      t_real real_result;
      PyArray_CastScalarToCtype(_in, static_cast<void*>(&real_result), 
                                PyArray_DescrFromType(numpy::type<t_real>::value));
      if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
      return real_result;
    } else throw errors::PythonTypeError("Unexpected type is not a real number.");
  }
  bool convertible_to_real(PyObject* _in) {
    using namespace DCProgs;
    if(PyFloat_Check(_in)) return true;
    if(PyLong_Check(_in)) return true;
#   ifndef DCPROGS_PYTHON3
      if(PyInt_Check(_in)) return true;
#   endif
    if(PyArray_CheckScalar(_in)) {
     Object<PyArray_Descr> descr = steal_ref<PyArray_Descr>(PyArray_DescrFromScalar(_in));
     if(not descr) throw errors::PythonErrorAlreadyThrown();
     return PyArray_CanCastSafely(numpy::type<t_real>::value, (~descr)->type_num) != 0;
    }
    return false;
  }
  int convert_to_real(PyObject* _in, void *_out) {
    using namespace DCProgs;
    try {
      t_real result = convert_to_real(_in);
      *static_cast<t_real*>(_out) = std::move(result);
      return 1;
    } DCPROGS_CATCH(return 0);
  }
  PyObject* convert_real_to_python(DCProgs::t_real const & _in, PyObject *base=NULL) {
    using namespace DCProgs;
#   ifndef DCPROGS_LONG_DOUBLE
      return PyFloat_FromDouble(_in);
#   else
      t_real copyme(_in);
      return PyArray_Scalar( static_cast<void*>(&copyme), 
                             PyArray_DescrFromType( numpy::type<t_real>::value),
                             base );
#   endif
  }
}%}

// Tells swig about our type hierarchy. 
// These types should make it easier to go from one system to another, but they do make it slightly
// more difficult for swig to understand our code.
%apply int { DCProgs::t_int }; 
// unsigned should be checked for negative values. Otherwise swig comes out with a truly unhelpfull
// message. Also checks for numpy scalar at the same type.
%typemap(in) DCProgs::t_uint {
  try {
    using namespace DCProgs;
    long integer;
    if(PyLong_Check($input)) integer = PyLong_AsLong($input);
#   ifndef DCPROGS_PYTHON3
      else if(PyInt_Check($input)) integer = PyInt_AsLong($input);
#   endif
    else if(PyArray_CheckScalar($input))  {
      PyArray_CastScalarToCtype($input, static_cast<void*>(&integer), 
                                PyArray_DescrFromType(numpy::type<long>::value));
      if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
    } else throw errors::PythonTypeError("Unexpected type in burst.");
    if(integer < 0) throw errors::PythonValueError("Input integer cannot be negative.");
    $1 = static_cast<t_uint>(integer);
  } DCPROGS_CATCH(SWIG_fail);
}
%typemap(out) DCProgs::t_uint = unsigned int;
%typemap(typecheck) DCProgs::t_uint = int;
%typemap(out) DCProgs::t_real { 
  try{ $result = convert_real_to_python($1); }
  DCPROGS_CATCH(SWIG_fail);
}
%typemap(in) DCProgs::t_real {
  try { $1 = convert_to_real($input); }
  DCPROGS_CATCH(SWIG_fail);
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_DOUBLE) DCProgs::t_real {
  $1 = convertible_to_real($input) ? 1: 0;
}
%typemap(out) DCProgs::t_rvector { 
  try { $result = DCProgs::numpy::wrap_to_numpy($1); }
  DCPROGS_CATCH(SWIG_fail);
}
%typemap(out) DCProgs::t_initvec { 
  try { $result = DCProgs::numpy::wrap_to_numpy($1); }
  DCPROGS_CATCH(SWIG_fail);
}
%typemap(out) DCProgs::t_rmatrix { 
  try { $result = DCProgs::numpy::wrap_to_numpy($1, NULL, true); }
  DCPROGS_CATCH(SWIG_fail);
};


%{
  namespace {
    PyObject* _dcprogs_dtype() {
      return PyArray_TypeObjectFromType(DCProgs::numpy::type<DCProgs::t_real>::value);
    }
  }
%}

// Bindings for eig, svd, and linalg stuff in general.
// Useful if t_real > double since numpy does not provide these functions only for 64-bit reals.
%{
  namespace {
    PyObject* eig(PyObject* _matrix) {
  
      using namespace DCProgs;
      DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
      Eigen::EigenSolver<t_rmatrix> eigsolver(matrix);
      if(eigsolver.info() != Eigen::Success)  {
        std::ostringstream sstr("Could not solve eigenvalue problem.");
        sstr << numpy_io(matrix) << "\n";
        throw errors::Mass(sstr.str());
      }
      t_cvector const eigenvalues = eigsolver.eigenvalues();
      t_cmatrix const eigenvectors = eigsolver.eigenvectors();
      if( (eigenvalues.imag().array().abs() > 1e-12).all() )  {
        Object<> py_eigs = steal_ref(numpy::wrap_to_numpy(eigenvalues));
        if(not py_eigs) throw errors::PythonErrorAlreadyThrown();
        Object<> py_vecs = steal_ref(numpy::wrap_to_numpy(eigenvectors));
        if(not py_vecs) throw errors::PythonErrorAlreadyThrown();
        return PyTuple_Pack(2, py_eigs.release(), py_vecs.release());
      }
      Object<> py_eigs = steal_ref(numpy::wrap_to_numpy(eigenvalues.real().eval()));
      if(not py_eigs) throw errors::PythonErrorAlreadyThrown();
      Object<> py_vecs = steal_ref(numpy::wrap_to_numpy(eigenvectors.real().eval()));
      if(not py_vecs) throw errors::PythonErrorAlreadyThrown();
      return PyTuple_Pack(2, py_eigs.release(), py_vecs.release());
    }
    PyObject* inv(PyObject* _matrix) {
      using namespace DCProgs;
      DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
      return numpy::wrap_to_numpy(matrix.inverse().eval()); 
    }
    double det(PyObject* _matrix) {
      using namespace DCProgs;
      DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
      return static_cast<t_real>(matrix.determinant());
    }
    PyObject* svd(PyObject* _matrix) {
  
      using namespace DCProgs;
      DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
      Eigen::JacobiSVD<t_rmatrix> svd(matrix, Eigen::ComputeThinU|Eigen::ComputeThinV);
      Object<> py_eigs = steal_ref(numpy::wrap_to_numpy(svd.singularValues().eval()));
      if(not py_eigs) throw errors::PythonErrorAlreadyThrown();
      Object<> pyU = steal_ref(numpy::wrap_to_numpy(svd.matrixU().eval()));
      if(not pyU) throw errors::PythonErrorAlreadyThrown();
      Object<> pyV = steal_ref(numpy::wrap_to_numpy(svd.matrixV().transpose().eval()));
      if(not pyV) throw errors::PythonErrorAlreadyThrown();
      return PyTuple_Pack(3, pyU.release(), py_eigs.release(), pyV.release());
    }
    PyObject* expm(PyObject* _matrix) {
      using namespace DCProgs;
      DCProgs::t_rmatrix const matrix = numpy::map_to_rmatrix(_matrix);
      return numpy::wrap_to_numpy(matrix.exp().eval());
    }
  }
%}
PyObject* _dcprogs_dtype();
%feature("autodoc", "Computes eigenvalues\n\n"
         "Only have this so we can compute eigenvalues of 128-bit reals, if dcprogs is so "
         "compiled.") eig;
PyObject* eig(PyObject* _matrix);
%feature("autodoc", "Computes inverse of a matrix\n\n"
         "Only have this so we can compute inverse of 128-bit reals, if dcprogs is so "
         "compiled.") inv;
PyObject* inv(PyObject* _matrix);
DCProgs::t_real det(PyObject* _matrix);
PyObject* svd(PyObject* _matrix);
PyObject* expm(PyObject* _matrix);

