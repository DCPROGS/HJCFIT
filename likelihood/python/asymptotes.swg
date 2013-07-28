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

%feature("autodoc", "Computes asymptotic values of missed-event G.") Asymptotes;
class Asymptotes {
  public:

  %extend {
    Asymptotes(DCProgs::DeterminantEq const &_det, PyObject *_array) {
       
       using namespace DCProgs;
       Object<> iterator = steal_ref(PyObject_GetIter(_array));
       if(not iterator) throw errors::PythonErrorAlreadyThrown();
       std::vector<Root> roots;
       while(Object<> item = steal_ref(PyIter_Next(~iterator))) {
         if(not PySequence_Check(~item))
           throw errors::PythonTypeError(
                    "Roots should be a list [(root, multiplicity)].");

         Object<> py_root = steal_ref(PySequence_GetItem(~item, 0));
         if(not py_root) throw errors::PythonErrorAlreadyThrown();
         t_real c_root;
         if(PyFloat_Check(~py_root)) c_root = (t_real)PyFloat_AsDouble(~py_root);
         else if(PyLong_Check(~py_root)) c_root = (t_real)PyLong_AsLong(~py_root);
         else if(PyInt_Check(~py_root)) c_root = (t_real)PyInt_AsLong(~py_root);
         else throw errors::PythonTypeError("Roots should be a real number.");

         Object<> py_mult = steal_ref(PySequence_GetItem(~item, 1));
         if(not py_mult) throw errors::PythonErrorAlreadyThrown();
         t_uint c_mult;
         if(PyLong_Check(~py_mult)) c_mult = static_cast<t_uint>(PyLong_AsLong(~py_mult));
         else if(PyInt_Check(~py_mult)) c_mult = static_cast<t_uint>(PyInt_AsLong(~py_mult));
         else throw errors::PythonTypeError("Multiplicity should an integer.");

         roots.emplace_back(c_root, c_mult);
       }
       if(PyErr_Occurred()) throw errors::PythonErrorAlreadyThrown();
       return new Asymptotes(_det, roots);
    };
  }

  DCProgs::t_rmatrix operator()(DCProgs::t_real _s) const;
  %extend {
    PyObject* operator()(PyObject *t) { return apply_numpy(t, *$self); }
  }
};

%clear Asymptotes;
