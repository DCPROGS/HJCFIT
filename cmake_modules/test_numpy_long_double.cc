#include <Python.h>
#include <numpy/npy_common.h>
#include <boost/mpl/int.hpp>

template<class T> class type;
      
template<> struct type<npy_longdouble> : public boost::mpl::int_<0> 
{
  typedef npy_longdouble np_type;
};
template<> struct type<npy_double> : public boost::mpl::int_<1> 
{
  typedef npy_double np_type;
};
int main() {return 0;}
