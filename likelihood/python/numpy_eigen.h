/***********************
    HJCFIT computes missed-events likelihood as described in
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

#ifndef HJCFIT_NUMPY_EIGEN
#define HJCFIT_NUMPY_EIGEN
#include <type_traits>
#include "../errors.h"
#include "object.h"

namespace HJCFIT {
  namespace numpy {
    //! An mpl integer defining the type.
    template<class T> class type;
    
#   ifdef HJCFIT_MACRO
#     error HJCFIT_MACRO already defined
#   endif
#   define HJCFIT_MACRO(TYPE_NAME, TYPE_NUMBER)                     \
    template<> struct type<TYPE_NAME> {                              \
      /*! Original Type */                                           \
      typedef TYPE_NAME np_type;                                     \
      /*! Associated  number */                                      \
      HJCFIT_INIT_CONSTEXPR(int value, TYPE_NUMBER);                \
    };                                                               \
    HJCFIT_DECL_CONSTEXPR(int type<TYPE_NAME>::value, TYPE_NUMBER)

    HJCFIT_MACRO(npy_cdouble,     NPY_CDOUBLE);
    HJCFIT_MACRO(npy_double,      NPY_DOUBLE);
    HJCFIT_MACRO(npy_float,       NPY_FLOAT);
    HJCFIT_MACRO(npy_longlong,    NPY_LONGLONG);
    HJCFIT_MACRO(npy_ulonglong,   NPY_ULONGLONG);
    HJCFIT_MACRO(npy_long,        NPY_LONG);
    HJCFIT_MACRO(npy_ulong,       NPY_ULONG);
    HJCFIT_MACRO(npy_int,         NPY_INT);
    HJCFIT_MACRO(npy_uint,        NPY_UINT);
    HJCFIT_MACRO(npy_short,       NPY_SHORT);
    HJCFIT_MACRO(npy_ushort,      NPY_USHORT);
    HJCFIT_MACRO(npy_byte,        NPY_BYTE);
    HJCFIT_MACRO(npy_ubyte,       NPY_UBYTE);

#   ifdef NUMPY_NPY_LONG_DOUBLE
      HJCFIT_MACRO(npy_longdouble, NPY_LONGDOUBLE);
      HJCFIT_MACRO(npy_clongdouble, NPY_CLONGDOUBLE);
      template<> struct type< std::complex<npy_longdouble> > {
        /*! Original Type */
        typedef npy_clongdouble np_type;
        /*! Associated  number */
        HJCFIT_INIT_CONSTEXPR(int value, NPY_CLONGDOUBLE); 
      };
      HJCFIT_DECL_CONSTEXPR(int type< std::complex<npy_longdouble> >::value, NPY_CLONGDOUBLE);
#   endif
#   ifdef NUMPY_NPY_BOOL
      HJCFIT_MACRO(npy_bool, NPY_BOOL);
#   else
      template<> struct type<bool> {
        /*! Original Type */
        typedef npy_bool np_type;
        /*! Associated  number */
        HJCFIT_INIT_CONSTEXPR(int value, NPY_BOOL); 
      };
      HJCFIT_DECL_CONSTEXPR(int type<bool>::value, NPY_BOOL);
#   endif 
      template<> struct type< std::complex<npy_double> > {
        /*! Original Type */
        typedef npy_cdouble np_type;
        /*! Associated  number */
        HJCFIT_INIT_CONSTEXPR(int value, NPY_CDOUBLE); 
      };
      HJCFIT_DECL_CONSTEXPR(int type< std::complex<npy_double> >::value, NPY_CDOUBLE);
      

#   undef HJCFIT_MACRO
    
    //! Convert/wrap a matrix to numpy.
    template<class T_DERIVED>
      PyObject* wrap_to_numpy_(Eigen::DenseBase<T_DERIVED> const &_in, PyObject *_parent = NULL,
                               bool _keep2d = false)
      {
        typedef type<typename Eigen::DenseBase<T_DERIVED>::Scalar> t_ScalarType;
        npy_intp dims[2] = { _in.rows(), _in.cols() };
        npy_int const ndims = (_in.cols() > 1 or _keep2d) ? 2: 1;
        if(_in.rows() == 0 or _in.cols() == 0)
          return PyArray_ZEROS(ndims, dims, t_ScalarType::value, _in.IsRowMajor ? 0: 1 );
        PyArrayObject *result = _parent == NULL ?
          (PyArrayObject*) PyArray_ZEROS(ndims, dims, t_ScalarType::value, _in.IsRowMajor ? 0: 1):
          (PyArrayObject*) PyArray_SimpleNewFromData( ndims, dims, t_ScalarType::value,
                                                      (void*)(&_in(0,0)) );
        if(result == NULL) return NULL;
        // If has a parent, do not copy data, just incref it as base.
        if(_parent != NULL) 
        {
          // For some reason, eigen is column major, whereas c++ is generally row major.
          if(PyArray_FLAGS(result) & NPY_ARRAY_C_CONTIGUOUS and not _in.IsRowMajor) 
            PyArray_CLEARFLAGS(result, NPY_ARRAY_C_CONTIGUOUS);
          else if((not (PyArray_FLAGS(result) & NPY_ARRAY_C_CONTIGUOUS)) and _in.IsRowMajor) 
            PyArray_ENABLEFLAGS(result, NPY_ARRAY_C_CONTIGUOUS);
          if(_in.cols() == 1)
            PyArray_STRIDES(result)[0] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
          else if(_in.IsRowMajor) 
          {
            PyArray_STRIDES(result)[0] = _in.outerStride() * sizeof(typename t_ScalarType::np_type);
            PyArray_STRIDES(result)[1] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
          }
          else 
          {
            PyArray_STRIDES(result)[0] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
            PyArray_STRIDES(result)[1] = _in.outerStride() * sizeof(typename t_ScalarType::np_type);
          }
          Py_INCREF(_parent);
          PyArray_SetBaseObject(result, _parent);
        }
        // otherwise, copy data.
        else
        {
          for(int i(0); i < _in.rows(); ++i)
            for(int j(0); j < _in.cols(); ++j)
              *((typename t_ScalarType::np_type*) PyArray_GETPTR2(result, i, j)) = _in(i, j);
        }
        return (PyObject*)result;
      }


    //! Convert/wrap a matrix to numpy.
    template<class T_DERIVED>
      PyObject* wrap_to_numpy(Eigen::DenseBase<T_DERIVED> const &_in, bool _keep2d=false) {
        PyObject* const result = wrap_to_numpy_(_in, NULL, _keep2d);
        PyArray_CLEARFLAGS((PyArrayObject*)result, NPY_ARRAY_WRITEABLE);
        return result;
      }
    //! Convert/wrap a matrix to numpy.
    template<class T_DERIVED>
      PyObject* wrap_to_numpy(Eigen::DenseBase<T_DERIVED> &_in, PyObject *_parent = NULL,
                              bool _keep2d=false) {
        PyObject* const result = wrap_to_numpy_(_in, _parent, _keep2d);
        PyArray_ENABLEFLAGS((PyArrayObject*)result, NPY_ARRAY_WRITEABLE);
        return result;
      }

    PyObject* wrap_to_numpy(t_cvector const &_in) {

      typedef type<t_cvector::Scalar> t_ScalarType;
      npy_intp dims[1] = { _in.size() };
      Object<PyArrayObject> result = steal_ref(reinterpret_cast<PyArrayObject*>(
          PyArray_SimpleNew(1, dims, t_ScalarType::value)
      )); 
      if(not result) throw errors::PythonErrorAlreadyThrown();
      
      t_real const *i_data = reinterpret_cast<const t_real(&)[2]>(_in(0));
      std::copy(i_data, i_data + 2 * _in.size(), static_cast<t_real*>(PyArray_DATA(~result)));
      if(PyArray_FLAGS(~result) & NPY_ARRAY_C_CONTIGUOUS and not _in.IsRowMajor) 
        PyArray_CLEARFLAGS(~result, NPY_ARRAY_C_CONTIGUOUS);
      return reinterpret_cast<PyObject*>(result.release());
    }
    PyObject* wrap_to_numpy(t_cmatrix const &_in) {

      typedef type<t_rmatrix::Scalar> t_ScalarType;
      npy_intp dims[2] = { _in.rows(), _in.cols() };
      Object<PyArrayObject> result = steal_ref(reinterpret_cast<PyArrayObject*>(
          PyArray_ZEROS(2, dims, t_ScalarType::value, _in.IsRowMajor ? 0: 1)
      )); 
      if(not result) throw errors::PythonErrorAlreadyThrown();

      t_real const *i_data = reinterpret_cast<const t_real(&)[2]>(_in(0));
      std::copy(i_data, i_data + 2 * _in.size(), static_cast<t_real*>(PyArray_DATA(~result)));
      if(PyArray_FLAGS(~result) & NPY_ARRAY_C_CONTIGUOUS and not _in.IsRowMajor) 
        PyArray_CLEARFLAGS(~result, NPY_ARRAY_C_CONTIGUOUS);

      return reinterpret_cast<PyObject*>(result.release());
    }

    namespace { namespace details {
      template<class T>
        Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, 0,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >
          wrap_to_eigen(PyArrayObject *_in) {
         
            typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> t_Stride;
            typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> t_Matrix;
            typedef Eigen::Map<t_Matrix, 0, t_Stride> t_Map;
         
            int const ndim = PyArray_NDIM(_in);
            npy_intp const * const strides = PyArray_STRIDES(_in);
            npy_intp const * const dims = PyArray_DIMS(_in);
            t_rvector::Index realdims[2] = { static_cast<t_rvector::Index>(dims[0]),
                                             ndim == 2 ? 
                                                static_cast<t_rvector::Index>(dims[1]): 1 };
            t_rvector::Index realstrides[2] = {
              static_cast<t_rvector::Index>(strides[0]),
              static_cast<t_rvector::Index>(ndim == 2 ? strides[1]: strides[0] * dims[0]) 
            };
            
            t_Map result( (T*)PyArray_DATA(_in), realdims[0], realdims[1], 
                          t_Stride( realstrides[1] * sizeof(char) / sizeof(T),
                                    realstrides[0] * sizeof(char) / sizeof(T) ) );
            return result;
         }
    }}


    //! \brief Converts numpy to an eigen matrix.
    //! \details It is best to check PyErr_Occurred after a call to this function.
    //! \param[in] _in a numpy array. 
    //! \return An eigen object which is a copy of the numpy input.
    HJCFIT::t_rmatrix map_to_rmatrix(PyObject *_in) {
       if(not PyArray_Check(_in)) {
          Object<> convert = steal_ref( 
            PyArray_FromObject(_in, HJCFIT::numpy::type<HJCFIT::t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw HJCFIT::errors::PythonErrorAlreadyThrown();
          return map_to_rmatrix(~convert);
       }
       int const type = PyArray_TYPE((PyArrayObject*)_in);
#      ifdef HJCFIT_MACRO
#        error HJCFIT_MACRO is already defined.
#      endif
#      define HJCFIT_MACRO(TYPE, NUM_TYPE)                                                        \
         if(type == NUM_TYPE)                                                                      \
           return details::wrap_to_eigen<TYPE>((PyArrayObject*)_in).cast<t_rmatrix::Scalar>(); 
        
       HJCFIT_MACRO( npy_float,      NPY_FLOAT)      
       else HJCFIT_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else HJCFIT_MACRO( npy_double,     NPY_DOUBLE     )
       else HJCFIT_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else HJCFIT_MACRO( npy_int,        NPY_INT        )
       else HJCFIT_MACRO( npy_uint,       NPY_UINT       )
       else HJCFIT_MACRO( npy_long,       NPY_LONG       )
       else HJCFIT_MACRO( npy_longlong,   NPY_LONGLONG   )
       else HJCFIT_MACRO( npy_ulonglong,  NPY_ULONGLONG  )
       else HJCFIT_MACRO( npy_ubyte,      NPY_BYTE       )
       else HJCFIT_MACRO( npy_short,      NPY_SHORT      )
       else HJCFIT_MACRO( npy_ushort,     NPY_USHORT     )
#      undef HJCFIT_MACRO
       throw HJCFIT::errors::PythonTypeError("Unexpect numpy array type");
       return t_rmatrix();
    }

    //! \brief Converts numpy to an initvec
    //! \param[in] _in a numpy array. 
    //! \return An eigen object which is a copy of the numpy input.
    HJCFIT::t_initvec map_to_initvec(PyObject *_in) {
       if(not PyArray_Check(_in)) {
          Object<> convert = steal_ref( 
            PyArray_FromObject(_in, HJCFIT::numpy::type<HJCFIT::t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw HJCFIT::errors::PythonErrorAlreadyThrown();
          return map_to_rmatrix(~convert);
       }
       // Check dimensionality
       npy_intp const N = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(_in));
       if(N == 0) throw HJCFIT::errors::PythonValueError("Input array is empty or a scalar.");
       npy_intp * dims = PyArray_DIMS(reinterpret_cast<PyArrayObject*>(_in));
       for(npy_intp i(0); i < N - 2; i++, ++dims) 
          if(*dims > 1)
            throw HJCFIT::errors::PythonValueError("Input array is not a (row) vector.");
          else if(*dims == 0) 
            throw HJCFIT::errors::PythonValueError("Input array is empty.");
       if(*dims == 0) throw HJCFIT::errors::PythonValueError("Input array is empty.");


       int const type = PyArray_TYPE((PyArrayObject*)_in);
#      ifdef HJCFIT_MACRO
#        error HJCFIT_MACRO is already defined.
#      endif
#      define HJCFIT_MACRO(TYPE, NUM_TYPE)                                                        \
         if(type == NUM_TYPE)                                                                      \
           return details::wrap_to_eigen<TYPE>((PyArrayObject*)_in).cast<t_rmatrix::Scalar>(); 
        
       HJCFIT_MACRO( npy_float,      NPY_FLOAT)      
       else HJCFIT_MACRO( npy_double,     NPY_DOUBLE     )
       else HJCFIT_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else HJCFIT_MACRO( npy_int,        NPY_INT        )
       else HJCFIT_MACRO( npy_uint,       NPY_UINT       )
       else HJCFIT_MACRO( npy_long,       NPY_LONG       )
       else HJCFIT_MACRO( npy_longlong,   NPY_LONGLONG   )
       else HJCFIT_MACRO( npy_ulonglong,  NPY_ULONGLONG  )
       else HJCFIT_MACRO( npy_ubyte,      NPY_BYTE       )
       else HJCFIT_MACRO( npy_short,      NPY_SHORT      )
       else HJCFIT_MACRO( npy_ushort,     NPY_USHORT     )
#      undef HJCFIT_MACRO
       throw HJCFIT::errors::PythonTypeError("Unexpect numpy array type");
       return t_initvec();
    }
    
    //! \brief Converts numpy to an rvector
    //! \param[in] _in a numpy array. 
    //! \return An eigen object which is a copy of the numpy input.
    HJCFIT::t_rvector map_to_rvector(PyObject *_in) {
       if(not PyArray_Check(_in)) {
          Object<> convert = steal_ref( 
            PyArray_FromObject(_in, HJCFIT::numpy::type<HJCFIT::t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw HJCFIT::errors::PythonErrorAlreadyThrown();
          return map_to_rmatrix(~convert);
       }
       
       // Check dimensionality
       npy_intp const N = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(_in));
       if(N == 0) throw HJCFIT::errors::PythonValueError("Input array is empty or a scalar.");
       npy_intp * dims = PyArray_DIMS(reinterpret_cast<PyArrayObject*>(_in));
       if(N > 2) {
         for(npy_intp i(0); i < N - 2; i++, ++dims) 
            if(*dims > 1)
              throw HJCFIT::errors::PythonValueError("Input array is not a (row) vector.");
            else if(*dims == 0) 
              throw HJCFIT::errors::PythonValueError("Input array is empty.");
       } 
       if(*dims == 0) throw HJCFIT::errors::PythonValueError("Input array is empty.");
       if(N >= 2) {
         npy_intp const dima = *dims;
         npy_intp const dimb = *(dims+1);
         if((dima == 1 and dimb != 1) or (dima != 1 and dimb == 1))
           throw HJCFIT::errors::PythonValueError("Input array is not a (column) vector.");
       }


       int const type = PyArray_TYPE((PyArrayObject*)_in);
#      ifdef HJCFIT_MACRO
#        error HJCFIT_MACRO is already defined.
#      endif
#      define HJCFIT_MACRO(TYPE, NUM_TYPE)                                                        \
         if(type == NUM_TYPE)                                                                      \
           return details::wrap_to_eigen<TYPE>((PyArrayObject*)_in).cast<t_rvector::Scalar>(); 
        
       HJCFIT_MACRO( npy_float,      NPY_FLOAT)      
       else HJCFIT_MACRO( npy_double,     NPY_DOUBLE     )
       else HJCFIT_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else HJCFIT_MACRO( npy_int,        NPY_INT        )
       else HJCFIT_MACRO( npy_uint,       NPY_UINT       )
       else HJCFIT_MACRO( npy_long,       NPY_LONG       )
       else HJCFIT_MACRO( npy_longlong,   NPY_LONGLONG   )
       else HJCFIT_MACRO( npy_ulonglong,  NPY_ULONGLONG  )
       else HJCFIT_MACRO( npy_ubyte,      NPY_BYTE       )
       else HJCFIT_MACRO( npy_short,      NPY_SHORT      )
       else HJCFIT_MACRO( npy_ushort,     NPY_USHORT     )
#      undef HJCFIT_MACRO
       throw HJCFIT::errors::PythonTypeError("Unexpect numpy array type");
       return t_rvector();
    }

    //! Cast data to given type from array type.
    template<class T> T cast(void *_data, int _type) {
      switch(_type) {
#       define HJCFIT_MACRO(TYPENAME)                                                             \
          case type<TYPENAME>::value: return static_cast<T>(*(static_cast<TYPENAME*>(_data)));
          HJCFIT_MACRO(npy_double);
          HJCFIT_MACRO(npy_float);
          HJCFIT_MACRO(npy_longlong);
          HJCFIT_MACRO(npy_ulonglong);
          HJCFIT_MACRO(npy_long);
          HJCFIT_MACRO(npy_ulong);
          HJCFIT_MACRO(npy_int);
          HJCFIT_MACRO(npy_uint);
          HJCFIT_MACRO(npy_short);
          HJCFIT_MACRO(npy_ushort);
          HJCFIT_MACRO(npy_byte);
          HJCFIT_MACRO(npy_ubyte);
#         ifdef NUMPY_NPY_LONG_DOUBLE
            HJCFIT_MACRO(npy_longdouble);
#         endif
#         ifdef NUMPY_NPY_BOOL
            HJCFIT_MACRO(npy_bool);
#         endif
#      undef HJCFIT_MACRO
      }
      throw HJCFIT::errors::PythonTypeError("Unexpect numpy array type");
      return T(0);
    }

  }
}
#endif

