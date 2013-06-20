%module likelihood
%{
#define SWIG_FILE_WITH_INIT
#  include <DCProgsConfig.h>
#  include <iostream>
#  include <sstream>

#  if NUMPY_VERSION_MINOR >= 7 
#    define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#  endif
#  include <numpy/arrayobject.h>
#  ifndef DCPROGS_NPY_ARRAY
#    define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS 
#    define NPY_ARRAY_WRITEABLE    NPY_WRITEABLE
#  endif
#  ifndef DCPROGS_NPY_ENABLEFLAGS
#    define PyArray_ENABLEFLAGS(ARRAY, FLAGS)  (ARRAY)->flags |= FLAGS
#    define PyArray_CLEARFLAGS(ARRAY, FLAGS)   (ARRAY)->flags &= (!FLAGS)
#    define PyArray_SetBaseObject(ARRAY, BASE) (ARRAY)->base   = BASE
#  endif
#  include "numpy_eigen.h"

#  include "../state_matrix.h"
%}

%include "exception.i"
%init %{ import_array();  %}

%apply int { t_int }; 
%apply int { DCProgs::t_int }; 

#ifdef DCPROGS_CATCH
# error DCPROGS_CATCH already defined.
#endif 
#define DCPROGS_CATCH                                                      \
    catch (DCProgs::errors::PythonTypeError &_e) {                         \
      PyErr_SetString(PyExc_TypeError, _e.what());                         \
      SWIG_fail;                                                           \
    }                                                                      \
    catch (DCProgs::errors::PythonValueError &_e) {                        \
      PyErr_SetString(PyExc_ValueError, _e.what());                        \
      SWIG_fail;                                                           \
    } catch(...) {                                                         \
      PyErr_SetString(PyExc_RuntimeError, "Caught unknown exception.");    \
      SWIG_fail;                                                           \
    }
//! General namespace for all things DCProgs.
namespace DCProgs {

  %exception StateMatrix::StateMatrix(PyObject *_in, int _nopen) {
    try { $action }
    DCPROGS_CATCH;
  }

  //! \brief State matrix that can  be partitioned into open/shut states.
  //! \details In practice, this is a two tuple with some helper functions to get corners.
  class StateMatrix {
    public:
 
 
    //! Number of open states.
    DCProgs::t_int nopen; 
    %typemap(in) DCProgs::t_rmatrix matrix { 
      try { $1 = DCProgs::numpy::map_to_rmatrix((PyArrayObject*)$input); }
      DCPROGS_CATCH;
    }
    %typemap(out) DCProgs::t_rmatrix { 
      try { $result = DCProgs::numpy::wrap_to_numpy(arg1->matrix, $self); }
      DCPROGS_CATCH;
    }
    //! The matrix itself.
    DCProgs::t_rmatrix matrix; 

    %clear DCProgs::t_rmatrix;

    StateMatrix();

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

      PyObject* __str__() {
        if($self->matrix.rows() == 0 or $self->matrix.cols() == 0) 
          return PyString_FromString("Un-initialized transition matrix.");
        std::ostringstream sstr;
        sstr << "Transition matrix with " << $self->nopen
             << " open states:\n" << $self->matrix << "\n";
        return PyString_FromString(sstr.str().c_str());
      }
    }
  };
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
}
#undef DCPROGS_CATCH
