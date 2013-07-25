#ifndef DCPROGS_NUMPY_EIGEN
#define DCPROGS_NUMPY_EIGEN
#include <type_traits>
#include "../errors.h"
#include "object.h"

namespace DCProgs {
  namespace numpy {
    //! An mpl integer defining the type.
    template<class T> class type;
    
#   ifdef DCPROGS_MACRO
#     error DCPROGS_MACRO already defined
#   endif
#   define DCPROGS_MACRO(TYPE_NAME, TYPE_NUMBER)                 \
    template<> struct type<TYPE_NAME> {                          \
      /*! Original Type */                                       \
      typedef TYPE_NAME np_type;                                 \
      /*! Associated  number */                                  \
      DCPROGS_INIT_CONSTEXPR(int value, TYPE_NUMBER);            \
    };                                                           \
    DCPROGS_DECL_CONSTEXPR(int type<TYPE_NAME>::value, TYPE_NUMBER);

    DCPROGS_MACRO(npy_longdouble,    NPY_LONGDOUBLE);
    DCPROGS_MACRO(npy_double,    NPY_DOUBLE);
    DCPROGS_MACRO(npy_float,     NPY_FLOAT);
    DCPROGS_MACRO(npy_longlong,  NPY_LONGLONG);
    DCPROGS_MACRO(npy_ulonglong, NPY_ULONGLONG);
    DCPROGS_MACRO(npy_long,      NPY_LONG);
    DCPROGS_MACRO(npy_ulong,     NPY_ULONG);
    DCPROGS_MACRO(npy_int,       NPY_INT);
    DCPROGS_MACRO(npy_uint,      NPY_UINT);
    DCPROGS_MACRO(npy_short,     NPY_SHORT);
    DCPROGS_MACRO(npy_ushort,    NPY_USHORT);
    DCPROGS_MACRO(npy_byte,      NPY_BYTE);
    DCPROGS_MACRO(npy_ubyte,     NPY_UBYTE);

#   ifdef DCPROGS_NPY_HAS_LONG_DOUBLE
      DCPROGS_MACRO(npy_longdouble, NPY_LONGDOUBLE);
#   endif
#   ifdef DCPROGS_NPY_HAS_BOOL
      DCPROGS_MACRO(npy_bool, NPY_BOOL);
#   else
      template<> struct type<bool> {
        /*! Original Type */
        typedef npy_bool np_type;
        /*! Associated  number */
        DCPROGS_INIT_CONSTEXPR(int value, NPY_BOOL); 
      };
      DCPROGS_DECL_CONSTEXPR(int type<bool>::value, NPY_BOOL);
#   endif 

#   undef DCPROGS_MACRO
    
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
            t_int realdims[2] = { static_cast<t_int>(dims[0]),
                                  ndim == 2 ? static_cast<t_int>(dims[1]): 1 };
            t_int realstrides[2] = {
              static_cast<t_int>(strides[0]),
              static_cast<t_int>(ndim == 2 ? strides[1]: strides[0] * dims[0]) 
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
    DCProgs::t_rmatrix map_to_rmatrix(PyObject *_in) {
       if(not PyArray_Check(_in)) {
          Object<> convert = steal_ref( 
            PyArray_FromObject(_in, DCProgs::numpy::type<DCProgs::t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw DCProgs::errors::PythonErrorAlreadyThrown();
          return map_to_rmatrix(~convert);
       }
       int const type = PyArray_TYPE((PyArrayObject*)_in);
#      ifdef DCPROGS_MACRO
#        error DCPROGS_MACRO is already defined.
#      endif
#      define DCPROGS_MACRO(TYPE, NUM_TYPE)                                                        \
         if(type == NUM_TYPE)                                                                      \
           return details::wrap_to_eigen<TYPE>((PyArrayObject*)_in).cast<t_rmatrix::Scalar>(); 
        
       DCPROGS_MACRO( npy_float,      NPY_FLOAT)      
       else DCPROGS_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else DCPROGS_MACRO( npy_double,     NPY_DOUBLE     )
       else DCPROGS_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else DCPROGS_MACRO( npy_int,        NPY_INT        )
       else DCPROGS_MACRO( npy_uint,       NPY_UINT       )
       else DCPROGS_MACRO( npy_long,       NPY_LONG       )
       else DCPROGS_MACRO( npy_longlong,   NPY_LONGLONG   )
       else DCPROGS_MACRO( npy_ulonglong,  NPY_ULONGLONG  )
       else DCPROGS_MACRO( npy_ubyte,      NPY_BYTE       )
       else DCPROGS_MACRO( npy_short,      NPY_SHORT      )
       else DCPROGS_MACRO( npy_ushort,     NPY_USHORT     )
#      undef DCPROGS_MACRO
       throw DCProgs::errors::PythonTypeError("Unexpect numpy array type");
       return t_rmatrix();
    }

    //! \brief Converts numpy to an initvec
    //! \param[in] _in a numpy array. 
    //! \return An eigen object which is a copy of the numpy input.
    DCProgs::t_initvec map_to_initvec(PyObject *_in) {
       if(not PyArray_Check(_in)) {
          Object<> convert = steal_ref( 
            PyArray_FromObject(_in, DCProgs::numpy::type<DCProgs::t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw DCProgs::errors::PythonErrorAlreadyThrown();
          return map_to_rmatrix(~convert);
       }
       // Check dimensionality
       npy_intp const N = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(_in));
       if(N == 0) throw DCProgs::errors::PythonValueError("Input array is empty or a scalar.");
       npy_intp * dims = PyArray_DIMS(reinterpret_cast<PyArrayObject*>(_in));
       for(npy_intp i(0); i < N - 2; i++, ++dims) 
          if(*dims > 1)
            throw DCProgs::errors::PythonValueError("Input array is not a (row) vector.");
          else if(*dims == 0) 
            throw DCProgs::errors::PythonValueError("Input array is empty.");
       if(*dims == 0) throw DCProgs::errors::PythonValueError("Input array is empty.");


       int const type = PyArray_TYPE((PyArrayObject*)_in);
#      ifdef DCPROGS_MACRO
#        error DCPROGS_MACRO is already defined.
#      endif
#      define DCPROGS_MACRO(TYPE, NUM_TYPE)                                                        \
         if(type == NUM_TYPE)                                                                      \
           return details::wrap_to_eigen<TYPE>((PyArrayObject*)_in).cast<t_rmatrix::Scalar>(); 
        
       DCPROGS_MACRO( npy_float,      NPY_FLOAT)      
       else DCPROGS_MACRO( npy_double,     NPY_DOUBLE     )
       else DCPROGS_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else DCPROGS_MACRO( npy_int,        NPY_INT        )
       else DCPROGS_MACRO( npy_uint,       NPY_UINT       )
       else DCPROGS_MACRO( npy_long,       NPY_LONG       )
       else DCPROGS_MACRO( npy_longlong,   NPY_LONGLONG   )
       else DCPROGS_MACRO( npy_ulonglong,  NPY_ULONGLONG  )
       else DCPROGS_MACRO( npy_ubyte,      NPY_BYTE       )
       else DCPROGS_MACRO( npy_short,      NPY_SHORT      )
       else DCPROGS_MACRO( npy_ushort,     NPY_USHORT     )
#      undef DCPROGS_MACRO
       throw DCProgs::errors::PythonTypeError("Unexpect numpy array type");
       return t_initvec();
    }
    
    //! \brief Converts numpy to an rvector
    //! \param[in] _in a numpy array. 
    //! \return An eigen object which is a copy of the numpy input.
    DCProgs::t_rvector map_to_rvector(PyObject *_in) {
       if(not PyArray_Check(_in)) {
          Object<> convert = steal_ref( 
            PyArray_FromObject(_in, DCProgs::numpy::type<DCProgs::t_real>::value, 0, 0)
          );
          if(PyErr_Occurred()) throw DCProgs::errors::PythonErrorAlreadyThrown();
          return map_to_rmatrix(~convert);
       }
       
       // Check dimensionality
       npy_intp const N = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(_in));
       if(N == 0) throw DCProgs::errors::PythonValueError("Input array is empty or a scalar.");
       npy_intp * dims = PyArray_DIMS(reinterpret_cast<PyArrayObject*>(_in));
       if(N > 2) {
         for(npy_intp i(0); i < N - 2; i++, ++dims) 
            if(*dims > 1)
              throw DCProgs::errors::PythonValueError("Input array is not a (row) vector.");
            else if(*dims == 0) 
              throw DCProgs::errors::PythonValueError("Input array is empty.");
       } 
       if(*dims == 0) throw DCProgs::errors::PythonValueError("Input array is empty.");
       if(N >= 2) {
         npy_intp const dima = *dims;
         npy_intp const dimb = *(dims+1);
         if((dima == 1 and dimb != 1) or (dima != 1 and dimb == 1))
           throw DCProgs::errors::PythonValueError("Input array is not a (column) vector.");
       }


       int const type = PyArray_TYPE((PyArrayObject*)_in);
#      ifdef DCPROGS_MACRO
#        error DCPROGS_MACRO is already defined.
#      endif
#      define DCPROGS_MACRO(TYPE, NUM_TYPE)                                                        \
         if(type == NUM_TYPE)                                                                      \
           return details::wrap_to_eigen<TYPE>((PyArrayObject*)_in).cast<t_rvector::Scalar>(); 
        
       DCPROGS_MACRO( npy_float,      NPY_FLOAT)      
       else DCPROGS_MACRO( npy_double,     NPY_DOUBLE     )
       else DCPROGS_MACRO( npy_longdouble, NPY_LONGDOUBLE )
       else DCPROGS_MACRO( npy_int,        NPY_INT        )
       else DCPROGS_MACRO( npy_uint,       NPY_UINT       )
       else DCPROGS_MACRO( npy_long,       NPY_LONG       )
       else DCPROGS_MACRO( npy_longlong,   NPY_LONGLONG   )
       else DCPROGS_MACRO( npy_ulonglong,  NPY_ULONGLONG  )
       else DCPROGS_MACRO( npy_ubyte,      NPY_BYTE       )
       else DCPROGS_MACRO( npy_short,      NPY_SHORT      )
       else DCPROGS_MACRO( npy_ushort,     NPY_USHORT     )
#      undef DCPROGS_MACRO
       throw DCProgs::errors::PythonTypeError("Unexpect numpy array type");
       return t_rvector();
    }

    //! Cast data to given type from array type.
    template<class T> T cast(void *_data, int _type) {
      switch(_type) {
#       define DCPROGS_MACRO(TYPENAME)                                                             \
          case type<TYPENAME>::value: return static_cast<T>(*(static_cast<TYPENAME*>(_data)));
          DCPROGS_MACRO(npy_longdouble);
          DCPROGS_MACRO(npy_double);
          DCPROGS_MACRO(npy_float);
          DCPROGS_MACRO(npy_longlong);
          DCPROGS_MACRO(npy_ulonglong);
          DCPROGS_MACRO(npy_long);
          DCPROGS_MACRO(npy_ulong);
          DCPROGS_MACRO(npy_int);
          DCPROGS_MACRO(npy_uint);
          DCPROGS_MACRO(npy_short);
          DCPROGS_MACRO(npy_ushort);
          DCPROGS_MACRO(npy_byte);
          DCPROGS_MACRO(npy_ubyte);
#         ifdef DCPROGS_NPY_HAS_LONG_DOUBLE
            DCPROGS_MACRO(npy_longdouble);
#         endif
#         ifdef DCPROGS_NPY_HAS_BOOL
            DCPROGS_MACRO(npy_bool);
#         endif
#      undef DCPROGS_MACRO
      }
      throw DCProgs::errors::PythonTypeError("Unexpect numpy array type");
      return T(0);
    }

  }
}
#endif

