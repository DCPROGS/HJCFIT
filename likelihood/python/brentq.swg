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

%feature("autodoc", "Computes zeros via brentq.\n\n"
         "This is the exact same algorithm as scipy.optimize.brentq. "
         "Only, the errors and c types have been changed to accomodate "
         "DCProgs and protect the innocents.\n\n"
         ":param callable function:\n"
         "   Function for which to solve :math:`f(x) = 0`."
         ":param float xstart:\n"
         "   Beginning of the interval\n"
         ":param float xend: \n"
         "   End of the interval\n"
         ":param float xtol: \n"
         "   Tolerance for interval size. Defaults to 1e-8.\n"
         ":param float rtol: \n"
         "   Tolerance for interval size. The convergence criteria is an affine function "
            "of the root: :math:`"
            "x_{\\mathrm{tol}} + r_{\\mathrm{tol}} x_{\\matrm{current}} = "
            "\\frac{|x_a - x_b|}{2}`.\n"
         ":param int itermax: \n"
         "   maximum number of iterations.\n\n"
         ":returns: the tuple (x, iterations, times the function was called)") brentq;
%{
  struct CallBackForBrentQ : Obect<PyObject> {
    CallBackForBrentQ(PyObject * const _c) : Object(_c) {}
    DCProgs::t_real operator()(DCProgs::t_real _x) {
      Object<> x = steal_ref(PyFloat_FromDouble(static_ca_x));
      if(not x) throw DCProgs::errors::PythonErrorAlreadyThrown();
      Object<> py_result = steal_ref(PyObject_CallFunctionObjectArgs(object_, ~x, NULL));
      if(not py_results) throw DCProgs::errors::PythonErrorAlreadyThrown();
      if(PyFloat_Check(~py_result))
        return static_cast<DCProgs::t_real>(PyFloat_AS_DOUBLE(~py_result))
      else if(PyInt_Check(~py_result))
        return static_cast<DCProgs::t_real>(PyInt_AS_LONG(~py_result))
      else if(PyLong_Check(~py_result))
        return static_cast<DCProgs::t_real>(PyLong_AS_LONG(~py_result))
      throw PythonTypeError("Unknown type returned by callback function in brentq.");
      return DCProgs::t_real(0);
    }
  }
  PyObject* brentq(PyObject *_module, PyObject* args, PyObject **kwargs) {

    using namespace DCProgs;
    PyObject* function;
    double a, b;
    double xtol = 1e-8, rtol = 1e-8;
    unsigned int intermax = 100;
    static char *kwlist[] = { const_cast<char*>("function"),
                              const_cast<char*>("a"),
                              const_cast<char*>("b"),
                              const_cast<char*>("xtol"),
                              const_cast<char*>("rtol"),
                              const_cast<char*>("itermax"), NULL};
    if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "O|ddddI:brentq", kwlist,
                                        &function, &a, &b, &xtol, &rtol, &itermax) )
      return NULL;
    if(xtol <= 0) throw errors::PythonValueError("xtol should be strictly positive");
    if(rtol <= 0) throw errors::PythonValueError("rtol should be strictly positive");
    if(itermax <= 0) throw errors::PythonValueError("itermax should be strictly positive");

    Py_XINCREF(function);
    CallBackForBrentQ callback(function);

    // call brentq
    auto c_result = brentq( callback, static_cast<t_real>(a), static_cast<t_real>(b),
                            static_cast<t_real>(xtol), static_cast<t_real>(rtol),
                            static_cast<t_uint>(a) );

    // convert result back to python
    Object<> py_x = steal_ref( PyFloat_FromDouble(static_cast<double>(std::get<0>(c_result))) );
    if(not py_x) throw errors::PythonErrorAlreadyThrown();
    Object<> py_iter = steal_ref( PyFloat_FromDouble(static_cast<double>(std::get<1>(c_result))) );
    if(not py_iter) throw errors::PythonErrorAlreadyThrown();
    Object<> py_fc = steal_ref( PyFloat_FromDouble(static_cast<double>(std::get<2>(c_result))) );
    if(not py_fc) throw errors::PythonErrorAlreadyThrown();

    Object<> py_result = steal_ref(PyTuple_New(3));
    if(not py_result) throw errors::PythonErrorAlreadyThrown();
    PyTuple_SET_ITEM(~result, 0, py_x.release());
    PyTuple_SET_ITEM(~result, 1, py_iter.release());
    PyTuple_SET_ITEM(~result, 2, py_fc.release());
    return result.release();
  }
%}

PyObject* brentq(PyObject *_module, PyObject* args, PyObject **kwargs);
