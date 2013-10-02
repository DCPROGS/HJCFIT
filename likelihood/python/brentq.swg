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

%{
  // Wraps a python object and makes it a c++ functor that brentq understands
  struct CallBackForBrentQ : Object<PyObject> {
    // Constructor acquires input object 
    CallBackForBrentQ(PyObject * const _c) : Object<PyObject>(acquire_ref(_c)) {}
    // Operator calls function object
    DCProgs::t_real operator()(DCProgs::t_real _x) {
      Object<> x = steal_ref(convert_real_to_python(_x));
      if(not x) throw DCProgs::errors::PythonErrorAlreadyThrown();
      Object<> py_result = steal_ref(PyObject_CallFunctionObjArgs(object_, ~x, NULL));
      if(not py_result) throw DCProgs::errors::PythonErrorAlreadyThrown();
      return convert_to_real(~py_result);
    }
  };
  // Wraps brentq so we can transform input and output
  PyObject* brentq(PyObject *function, DCProgs::t_real a, DCProgs::t_real b,
                   DCProgs::t_real xtol, DCProgs::t_real rtol, DCProgs::t_uint itermax) {

    using namespace DCProgs;
    CallBackForBrentQ callback(function);

    // call brentq
    auto c_result = DCProgs::brentq(callback, a, b, xtol, rtol, itermax);

    // convert result back to python
    Object<> py_x = steal_ref( convert_real_to_python(std::get<0>(c_result)) );
    if(not py_x) throw errors::PythonErrorAlreadyThrown();
    Object<> py_iter = steal_ref( convert_real_to_python(std::get<1>(c_result)) );
    if(not py_iter) throw errors::PythonErrorAlreadyThrown();
    Object<> py_fc = steal_ref( convert_real_to_python(std::get<2>(c_result)) );
    if(not py_fc) throw errors::PythonErrorAlreadyThrown();

    Object<> py_result = steal_ref(PyTuple_New(3));
    if(not py_result) throw errors::PythonErrorAlreadyThrown();
    PyTuple_SET_ITEM(~py_result, 0, py_x.release());
    PyTuple_SET_ITEM(~py_result, 1, py_iter.release());
    PyTuple_SET_ITEM(~py_result, 2, py_fc.release());
    return py_result.release();
  }
%}

%rename("_brentq") brentq;
PyObject* brentq(PyObject *function, DCProgs::t_real a, DCProgs::t_real b, 
                 DCProgs::t_real xtol, DCProgs::t_real rtol, DCProgs::t_uint itermax);
%pythoncode %{
  def brentq(function, a, b, xtol=1e-12, rtol=4.440892098e-16, itermax=100):
    """ Computes zeros via brentq.
         
        This is the exact same algorithm as scipy.optimize.brentq. Only, the errors and c types
        have been changed to accomodate DCProgs and protect the innocents.

        :param callable function:
           Function for which to solve :math:`f(x) = 0`.
        :param float xstart:
           Beginning of the interval
        :param float xend:
           End of the interval
        :param float xtol: 
           Tolerance for interval size. 
        :param float rtol: 
           Tolerance for interval size. The convergence criteria is an affine function of the
           root:
           :math:`x_{\\mathrm{tol}}+r_{\\mathrm{tol}} x_{\\mathrm{current}} = \\frac{|x_a-x_b|}{2}`.
        :param int itermax: 
           maximum number of iterations.

        :returns: the tuple (x, iterations, times the function was called)
    """
    if xtol <= 0: raise ValueError("xtol must be strictly positive");
    if rtol <= 0: raise ValueError("rtol must be strictly positive");
    if itermax <= 0: raise ValueError("itermax must be strictly positive");
    if b <= a: raise ValueError("Interval is not valid: {0} >= {1}.".format(a, b))
    return _likelihood._brentq(function, a, b, xtol, rtol, itermax)
%}

