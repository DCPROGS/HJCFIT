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

namespace {

    //! \brief Helps out applying functor to a function.
    //! \details Applying a functor to a numpy array is mostly the same, whether that functor
    //!          returns a scalar or an array/matrix. `NumpyApplySpecialized` contain the snippet of
    //!          code that do change depending on what the functor returns.
    //! \tparam T Type of the ouput numpy array.
    //! \tparam reduce How many dimensions are gobled up by functor.
    template<class T, size_t reduce=0>
      class NumpyApplySpecialized : public DCProgs::numpy::type<T> {
        public:
          using DCProgs::numpy::type<T>::value;
          using typename DCProgs::numpy::type<T>::np_type;
       
          //! Creates functor.
          NumpyApplySpecialized(PyArrayObject *_in) : input_(acquire_ref(_in)) {
            in_N_ = PyArray_NDIM(_in);
            out_stride_ = 0;
          };
       
          //! Creates output array.
          PyArrayObject* output(T const &) {
            PyArrayObject * result = reinterpret_cast<PyArrayObject*>(
              PyArray_SimpleNew(in_N_ - reduce, PyArray_DIMS(~input_), value)
            );
            size_t const nbytes = sizeof(T) / sizeof(char);
            out_stride_ = PyArray_STRIDES(result)[in_N_-reduce-1] / nbytes;
            return result;
          }
          template<class T_DATA_ITERATOR>
            void set_element(T_DATA_ITERATOR _data_itr, T const &_item) const { *_data_itr = _item; }
          template<class T_DATA>
            void increment_data_pointer(T_DATA & _data_itr) const { _data_itr += out_stride_; }
       
        private:
          //! \brief Array that is the input to the functor
          //! \details Used to figure out size of output array.
          Object<PyArrayObject> input_;
          //! Size of input array.
          npy_intp in_N_;
          //! Stride to move to next element.
          npy_intp out_stride_;
     };

    template<class T, int I, int N, int IS, int NS, size_t reduce>
      struct NumpyApplySpecialized<Eigen::Matrix<T, I, N, 0, IS, NS>, reduce> : public DCProgs::numpy::type<T> {
        public:
          using DCProgs::numpy::type<T>::value;
          using typename DCProgs::numpy::type<T>::np_type;
          //! Number of dims of functor output.
          npy_intp const static add_dims;
  
          //! Creates functor.
          NumpyApplySpecialized(PyArrayObject *_in) : input_(acquire_ref(_in)) {
            in_N_ = PyArray_NDIM(_in);
          };
  
          //! Creates output array.
          PyArrayObject* output(Eigen::Matrix<T, I, N, 0, IS, NS> const &_first) {
            std::vector<npy_intp> dims(in_N_ + add_dims - reduce);
            std::copy(PyArray_DIMS(~input_),
                      PyArray_DIMS(~input_) + in_N_ - reduce, dims.begin());
            if(I != 1 and N != 1) {
              out_sizes_[0] = _first.rows();
              out_sizes_[1] = _first.cols();
              dims[in_N_] = static_cast<npy_intp>(out_sizes_[0]);
              dims.back() = static_cast<npy_intp>(out_sizes_[1]);
            } else if (I == 1 or N == 1) {
              out_sizes_[0] = I == 1 ? 1: _first.size();
              out_sizes_[1] = N == 1 ? 1: _first.size();
              dims.back() = _first.size(); 
            }
            
            PyArrayObject * const result = reinterpret_cast<PyArrayObject*>(
                PyArray_SimpleNew(dims.size(), &(dims[0]), value)
            );
            size_t const nbytes = sizeof(T) / sizeof(char);
            npy_intp const * const ptr_strides = PyArray_STRIDES(result) + in_N_ - reduce;
            data_strides_[0] = I == 1 ? 1: (ptr_strides[0] / nbytes);
            data_strides_[1] = N == 1 ? 1: (ptr_strides[I == 1 ? 0: 1] / nbytes);
            out_stride_ = ptr_strides[-1] / nbytes;
            return result;
          }
          template<class T_DATA_ITERATOR, class T_DERIVED>
            void set_element( T_DATA_ITERATOR  _data_itr,
                              Eigen::DenseBase<T_DERIVED> const &_item ) const {
              // Strides are inverted because Eigen does fortran style and numpy does C style by
              // default. 
              t_Map( _data_itr, out_sizes_[0], out_sizes_[1],
                     t_Stride(data_strides_[1], data_strides_[0]) ) = _item;
            }
          template<class T_DATA>
            void increment_data_pointer(T_DATA & _data_itr) const { 
              _data_itr += out_stride_; 
            }
       
        private:
          typedef Eigen::Stride< ::Eigen::Dynamic, ::Eigen::Dynamic> t_Stride;
          typedef Eigen::Map<DCProgs::t_rmatrix, 0, t_Stride> t_Map;
          //! \brief Array that is the input to the functor
          //! \details Used to figure out size of output array.
          Object<PyArrayObject> input_;
          //! Size of input array.
          npy_intp in_N_;
          //! Stride to move to next element.
          npy_intp out_stride_;
          //! Stride to move to next element.
          typename Eigen::Matrix<T, I, N, 0, IS, NS>::Index out_sizes_[2];
          //! Creates stride for copying to output.
          typename Eigen::Matrix<T, I, N, 0, IS, NS>::Index data_strides_[2];
      };
    template<class T, int I, int N, int IS, int NS, size_t reduce>
      npy_intp const NumpyApplySpecialized<Eigen::Matrix<T, I, N, 0, IS, NS>, reduce>
                     :: add_dims = static_cast<DCProgs::t_int>(I != 1)
                                   + static_cast<DCProgs::t_int>(N != 1);

    // \brief Applies functor to an input array.
    // \details The functor should return a scalar or a matrix which always has the same size.
    //          Output array will have same shape as input (+ extra dimensions if functor returns a
    //          vector or matrix). 
    template<class T_FUNC> PyObject* apply_numpy(PyObject* _in, T_FUNC const & _func) {

      using namespace DCProgs;
      typedef typename NumpyApplySpecialized<decltype(_func(0))>::np_type t_numpy_type;
      npy_intp const numpy_value = NumpyApplySpecialized<decltype(_func(0))>::value;

      if(not PyArray_Check(_in)) {
        Object<PyObject> convert = steal_ref( static_cast<PyObject*>(
          PyArray_FromObject(_in, numpy_value, 0, 0)
        ));
        if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
        return apply_numpy(~convert, _func);
      }

      PyArrayObject * const array = reinterpret_cast<PyArrayObject*>(_in);
      NumpyApplySpecialized<decltype(_func(0))> versatile_code(array); 
      // First, gets iterator and computes first result:
      // We need the size of the resulting matrix.
      auto in_iter = steal_ref(PyArray_IterNew(_in));
      if(not in_iter) return NULL; 
      // If empty array, returns empty array.
      if(not PyArray_ITER_NOTDONE(~in_iter))
        return PyArray_SimpleNew(PyArray_NDIM(array), PyArray_DIMS(array), numpy_value);


      // Compute first result
      int const in_numpy_value = PyArray_TYPE(array);
      auto const first_result = _func(
        numpy::cast<t_numpy_type>(PyArray_ITER_DATA(~in_iter), in_numpy_value)
      );
      PyArray_ITER_NEXT(~in_iter);

      // Now construct output array.
      Object<PyArrayObject> result = steal_ref( versatile_code.output(first_result) );
      if(not result) return NULL;

      // Add first element to array.
      t_real * i_out = static_cast<t_real*>(PyArray_DATA((~result)));
      versatile_code.set_element(i_out, first_result);

      // finally, loop over array
      // element_stride: to increment i_out from result to result.
      while( PyArray_ITER_NOTDONE(~in_iter) ) {
        versatile_code.increment_data_pointer(i_out);
        t_real const s = numpy::cast<t_numpy_type>(PyArray_ITER_DATA(~in_iter), in_numpy_value);
        versatile_code.set_element(i_out, _func(s));
        PyArray_ITER_NEXT(~in_iter);
      }
      return reinterpret_cast<PyObject*>(result.release());
    }        

}
