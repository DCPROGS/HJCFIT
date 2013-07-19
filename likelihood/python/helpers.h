namespace {

    // Applies functor to an input array
    // Output array will have same shape as input. It will always be of type real.
    template<class T_FUNC> PyObject* apply_real(PyObject* _in, T_FUNC const & _func) {

      if(not PyArray_Check(_in)) {
        Object<PyObject> convert = steal_ref( (PyObject*)
          PyArray_FromObject(_in, DCProgs::numpy::type<DCProgs::t_real>::value, 0, 0)
        );
        if(PyErr_Occurred()) throw DCProgs::errors::PythonErrorAlreadyThrown();
        return apply_real(~convert, _func);
      }
      PyArrayObject * const array = (PyArrayObject*)_in;
      Object<PyArrayObject> result = steal_ref( (PyArrayObject*)
          PyArray_SimpleNew( PyArray_NDIM(array), PyArray_DIMS(array),
                             DCProgs::numpy::type<DCProgs::t_real>::value )
      );
      if(not result) return NULL;

      auto in_iter = steal_ref(PyArray_IterNew(_in));
      if(not in_iter) return NULL; 
      
      int const typenum = PyArray_TYPE(array);
      
      DCProgs::t_real * i_out = (DCProgs::t_real*) PyArray_DATA((~result));
      while( PyArray_ITER_NOTDONE(~in_iter) ) {
        DCProgs::t_real const s =
          DCProgs::numpy::cast<DCProgs::t_real>(PyArray_ITER_DATA(~in_iter), typenum);
        *i_out = _func(s);
        ++i_out;
        PyArray_ITER_NEXT(~in_iter);
      }
      return reinterpret_cast<PyObject*>(result.release());
    }        

    // Applies functor to an input array. The functor should return a matrix which always has the
    // same size.
    // Output array will have same shape as input + 2 extra dimensions. It will always be of type real.
    template<class T_FUNC> PyObject* apply_rmatrix(PyObject* _in, T_FUNC const & _func) {

      using namespace DCProgs;
      if(not PyArray_Check(_in)) {
        Object<PyObject> convert = steal_ref( static_cast<PyObject*>(
          PyArray_FromObject(_in, numpy::type<t_real>::value, 0, 0)
        ));
        if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
        return apply_rmatrix(~convert, _func);
      }
      PyArrayObject * const array = (PyArrayObject*)_in;
      // First, gets iterator and computes first result:
      // We need the size of the resulting matrix.
      auto in_iter = steal_ref(PyArray_IterNew(_in));
      if(not in_iter) return NULL; 
      // Empty array, returns None.
      if(not PyArray_ITER_NOTDONE(~in_iter)) Py_RETURN_NONE;


      // Compute first result
      int const typenum = PyArray_TYPE(array);
      t_rmatrix const first_result = _func(
        numpy::cast<t_real>(PyArray_ITER_DATA(~in_iter), typenum)
      );
      PyArray_ITER_NEXT(~in_iter);
      t_rmatrix::Index const nrows = first_result.rows();
      t_rmatrix::Index const ncols = first_result.cols();

      // Now construct output array.
      npy_intp const in_N = PyArray_NDIM(array);
      std::vector<npy_intp> dims(in_N + 2);
      std::copy(PyArray_DIMS(array), PyArray_DIMS(array) + in_N, dims.begin());
      dims[in_N] = static_cast<npy_intp>(nrows);
      dims.back() = static_cast<npy_intp>(ncols);

      Object<PyArrayObject> result = steal_ref( reinterpret_cast<PyArrayObject*> (
            PyArray_SimpleNew(in_N + 2, &(dims[0]), numpy::type<t_real>::value)
      ));
      if(not result) return NULL;

      // Stuff we need to create eigen maps.
      typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> t_Stride;
      typedef Eigen::Map<t_rmatrix, 0, t_Stride> t_Map;
      static_assert(sizeof(t_real) % sizeof(char) == 0, "That really shouldn't happen");
      size_t const nbytes = sizeof(t_real) / sizeof(char);
      npy_intp const row_stride = PyArray_STRIDES(~result)[in_N] / nbytes;
      npy_intp const col_stride = PyArray_STRIDES(~result)[in_N+1] / nbytes;
      t_Stride const strides(col_stride, row_stride); 
      
      // finally, loop over array
      // element_stride: to increment i_out from result to result.
      assert(PyArray_STRIDES(array)[in_N-1] % nbytes == 0);
      npy_intp const element_stride = PyArray_STRIDES(~result)[in_N-1] / nbytes;
      t_real * i_out = static_cast<t_real*>(PyArray_DATA((~result)));
      t_Map(i_out, nrows, ncols, t_Stride(col_stride, row_stride)) = first_result;
      while( PyArray_ITER_NOTDONE(~in_iter) ) {
        i_out += element_stride;
        t_real const s =
          numpy::cast<t_real>(PyArray_ITER_DATA(~in_iter), typenum);
        t_Map(i_out, nrows, ncols, strides) = _func(s);
        PyArray_ITER_NEXT(~in_iter);
      }
      return reinterpret_cast<PyObject*>(result.release());
    }        
}
