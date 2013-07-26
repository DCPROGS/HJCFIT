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
       
       Object<> iterator = steal_ref(PyObject_GetIter(_array));
       if(not iterator) throw DCProgs::errors::PythonErrorAlreadyThrown();
       std::vector<DCProgs::Root> roots;
       while(Object<> item = steal_ref(PyIter_Next(~iterator))) {
         if(not PySequence_Check(~item))
           throw DCProgs::errors::PythonTypeError(
                    "Roots should be a list [(root, multiplicity)].");

         Object<> py_root = steal_ref(PySequence_GetItem(~item, 0));
         if(not py_root) throw DCProgs::errors::PythonErrorAlreadyThrown();
         DCProgs::t_real c_root;
         if(PyFloat_Check(~py_root)) c_root = (DCProgs::t_real)PyFloat_AsDouble(~py_root);
         else if(PyLong_Check(~py_root)) c_root = (DCProgs::t_real)PyLong_AsLong(~py_root);
         else if(PyInt_Check(~py_root)) c_root = (DCProgs::t_real)PyInt_AsLong(~py_root);
         else throw DCProgs::errors::PythonTypeError("Roots should be a real number.");

         Object<> py_mult = steal_ref(PySequence_GetItem(~item, 1));
         if(not py_mult) throw DCProgs::errors::PythonErrorAlreadyThrown();
         DCProgs::t_int c_mult;
         if(PyLong_Check(~py_mult)) c_mult = (DCProgs::t_int)PyLong_AsLong(~py_mult);
         else if(PyInt_Check(~py_mult)) c_mult = (DCProgs::t_int)PyInt_AsLong(~py_mult);
         else throw DCProgs::errors::PythonTypeError("Multiplicity should an integer.");

         roots.emplace_back(c_root, c_mult);
       }
       if(PyErr_Occurred()) throw DCProgs::errors::PythonErrorAlreadyThrown();
       return new DCProgs::Asymptotes(_det, roots);
    };
  }

  DCProgs::t_rmatrix operator()(DCProgs::t_real _s) const;
  %extend {
    PyObject* operator()(PyObject *t) { return apply_numpy(t, *$self); }
  }
};

%clear Asymptotes;
