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
#  include "../asymptotes.h"

#  include "helpers.h"
%}

%include "exception.i"
%init %{ import_array();  %}

%apply int { t_int }; 
%apply int { DCProgs::t_int }; 
%apply double { t_real }; 
%apply double { DCProgs::t_real }; 

#ifdef DCPROGS_CATCH
# error DCPROGS_CATCH already defined.
#endif 
#define DCPROGS_CATCH(ONERROR)                                                \
    catch (DCProgs::errors::PythonTypeError &_e) {                            \
      PyErr_SetString(PyExc_TypeError, _e.what());                            \
      ONERROR;                                                                \
    } catch (DCProgs::errors::PythonValueError &_e) {                         \
      PyErr_SetString(PyExc_ValueError, _e.what());                           \
      ONERROR;                                                                \
    } catch (DCProgs::errors::PythonErrorAlreadyThrown &_e) {                 \
      ONERROR;                                                                \
    } catch (DCProgs::errors::Python &_e) {                                   \
      PyErr_SetString(PyExc_RuntimeError, "Caught unspecified exception.");   \
      ONERROR;                                                                \
    } catch(...) {                                                            \
      PyErr_SetString(PyExc_RuntimeError, "Caught unknown exception.");       \
      ONERROR;                                                                \
    }
//! General namespace for all things DCProgs.
namespace DCProgs {

  %exception StateMatrix::StateMatrix(PyObject *_in, int _nopen) {
    try { $action }
    DCPROGS_CATCH(SWIG_fail);
  }

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


  %exception DeterminantEq::DeterminantEq(PyObject *_in, int _nopen, double _s, bool _doopen=true) {
    try { $action }
    DCPROGS_CATCH(SWIG_fail);
  }
  %exception DeterminantEq::DeterminantEq(StateMatrix *_in, double _s, bool _doopen=true) {
    try { $action }
    DCPROGS_CATCH(SWIG_fail);
  }

  //! Adds equilibrium functor. 
  class DeterminantEq {

    public: 

    %rename(_get_tau) get_tau() const; 
    %rename(_set_tau) set_tau(DCProgs::t_real); 

    %extend {
      DeterminantEq(PyObject *_in, int _nopen, double _s, bool _doopen=true) {
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
                                           _s, _doopen );
      }
      DeterminantEq(StateMatrix *_in, double _s, bool _doopen=true) {
        if(_in->nopen <= 0) 
          throw DCProgs::errors::PythonValueError("Number of open states cannot be negative or zero.");
        if(_in->matrix.rows() != _in->matrix.cols()) 
          throw DCProgs::errors::PythonValueError("Expected a square state matrix on input.");
        if(_in->matrix.rows() == 0)
          throw DCProgs::errors::PythonValueError("Expected a non-empty square matrix on input.");
        if(_in->nopen >= _in->matrix.rows()) 
          throw DCProgs::errors::PythonValueError("Number of closed states cannot be zero.");
        return new DCProgs::DeterminantEq( *_in, _s, _doopen );
      }
    }

    DCProgs::t_real get_tau() const;
    void set_tau(DCProgs::t_real);
    %pythoncode %{
      __swig_getmethods__["tau"] = _get_tau
      __swig_setmethods__["tau"] = _set_tau
      if _newclass: tau = property(_get_tau, _set_tau)
    %}

    DCProgs::t_real operator()(DCProgs::t_real _s);
    DCProgs::t_real operator()(DCProgs::t_real _s, DCProgs::t_real _tau);
    %extend { 
      PyObject* operator()(PyObject * _s) {
        try { return apply_real(_s, *$self); }
        DCPROGS_CATCH(return NULL;);
      }
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

//   %extend {
//     PyObject* H(DCProgs::t_real _s) const {
//       return DCProgs::numpy::wrap_to_numpy($self->DCProgs::DeterminantEq::H(_s));
//     }
//   }
  };
}
#undef DCPROGS_CATCH
