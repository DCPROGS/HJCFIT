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

%feature("autodoc", "Computes approximate missed-event survivor function.\n\n"
         "There are two ways to create this object: (i) the roots of two determinantal "
         "equations can be given explicitely, or (ii) they can be determined automatically.\n\n"
         ">>> approx = ApproxSurvivor(determinant_af, roots_af, determinant_fa, roots_fa)\n\n"
         ">>> approx = ApproxSurvivor(qmatrix, tau [, xtol, rtol, itermax])\n\n"
         "The parameters in bracket are optional.\n"
         ":param DeterminantEq determinant_af:\n"
         "  The determinantal equation for the af block.\n"
         ":param roots_af:\n"
         "  The roots of the af determinantal equation. The should come in the format "
            "`[(root, multiplicity), (root, multiplicity), ...]`.\n" 
         ":param DeterminantEq determinant_fa:\n"
         "  The determinantal equation for the fa block. It should likely be "
           "the transpose of determinant_fa.\n"
         ":param roots_fa:\n"
         "  The roots of the fa determinantal equation. The should come in the format "
            "`[(root, multiplicity), (root, multiplicity), ...]`.\n" 
         ":param float tau:\n"
         "  The maximum length of the missed events.\n"
         ":param qmatrix: \n"
         "   A :class:`QMatrix` instance. \n"
         ":param int nmax: \n"
         "   The exact missed event likelihood will be computed for times "
            ":math:`t \\in [0, n_{\\mathrm{max}} \\tau`. It defaults to 2.\n"
         ":param float xtol: \n"
         "   Tolerance criteria when computing roots using "
            "py:meth:`dcprogs.likelihood.brentq`. Defaults to 1e-10.\n"
         ":param float rtol: \n"
         "   Tolerance criteria when computing roots using "
            "py:meth:`dcprogs.likelihood.brentq`. Defaults to 1e-10\n"
         ":param float itermax: \n"
         "   Maximum number of iterations when computing roots using "
            "py:meth:`dcprogs.likelihood.brentq`. Defaults to 100.\n") ApproxSurvivor;
%feature("autodoc", "Open to close transitions") ApproxSurvivor::af;
%feature("autodoc", "Closed to open transitions") ApproxSurvivor::fa;
class ApproxSurvivor {
  public:
  
  %extend {
    // Creates object from determinants and roots.
    ApproxSurvivor(DCProgs::DeterminantEq const &_af, PyObject * _roots_af, 
                   DCProgs::DeterminantEq const &_fa, PyObject * _roots_fa ) {
      using namespace DCProgs;
      std::vector<Root> roots_af = convert_to_roots(_roots_af);
      std::vector<Root> roots_fa = convert_to_roots(_roots_fa);
      return new ApproxSurvivor(_af, roots_af, _fa, roots_fa);
    }
    ApproxSurvivor( QMatrix const * const _matrix, DCProgs::t_real _tau, 
                    DCProgs::t_int _nmax=2, DCProgs::t_real _xtol = 1e-10,
                    DCProgs::t_real _rtol = 1e-10, DCProgs::t_int _itermax = 100 ) {
      using namespace DCProgs;
      DeterminantEq const determinant_af(*_matrix, _tau);
      DeterminantEq const determinant_fa = determinant_af.transpose(); 
      std::vector<Root> const roots_af = find_roots(determinant_af, _xtol, _rtol, _itermax);
      std::vector<Root> const roots_fa = find_roots(determinant_fa, _xtol, _rtol, _itermax);
      return new DCProgs::ApproxSurvivor(determinant_af, roots_af, determinant_fa, roots_fa);
    }
  }

  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);
  // Add the ability to pass in numpy arrays.
  %extend {
    PyObject* af(PyObject *_s) {
      return apply_numpy(_s, [$self](DCProgs::t_real s) { return $self->af(s); }); 
    }
    PyObject* fa(PyObject *_s) {
      return apply_numpy(_s, [$self](DCProgs::t_real s) { return $self->fa(s); }); 
    }
  }

  %extend {
    PyObject* const af_components;
    PyObject* const fa_components;
  }
  %{ namespace { 
    template<class T> 
      PyObject* affa_components_get(Py_ssize_t const _N, T const &_getstuff) {
        Object<> result = steal_ref( PyList_New(_N) );
        if(not result) throw DCProgs::errors::PythonErrorAlreadyThrown();
        for(Py_ssize_t i(0); i < _N; ++i) {
          DCProgs::Asymptotes::t_MatrixAndRoot const &item = _getstuff(i);
          Object<> py_item = steal_ref( PyTuple_New(2) );
          if(not py_item) throw DCProgs::errors::PythonErrorAlreadyThrown();
          Object<> py_root = steal_ref(PyFloat_FromDouble(std::get<1>(item)));
          if(not py_root) throw DCProgs::errors::PythonErrorAlreadyThrown();
          DCProgs::t_rmatrix const & matrix = std::get<0>(item);
          Object<> py_matrix = steal_ref(DCProgs::numpy::wrap_to_numpy(matrix));
          if(not py_matrix) throw DCProgs::errors::PythonErrorAlreadyThrown();
          PyArray_CLEARFLAGS((PyArrayObject*) ~py_matrix, NPY_ARRAY_WRITEABLE);
          PyTuple_SET_ITEM(~py_item, 0, py_matrix.release());
          PyTuple_SET_ITEM(~py_item, 1, py_root.release());
          PyList_SET_ITEM(~result, i, py_item.release());
        }
        return result.release();
      }
    PyObject* DCProgs_ApproxSurvivor_af_components_get(DCProgs::ApproxSurvivor *_in) {
        return affa_components_get( (Py_ssize_t)_in->nb_af_components(),
                                    [&_in](Py_ssize_t i) {
                 return _in->get_af_components(static_cast<DCProgs::t_int>(i)); } );
    }
    PyObject* DCProgs_ApproxSurvivor_fa_components_get(DCProgs::ApproxSurvivor *_in) {
        return affa_components_get( (Py_ssize_t)_in->nb_fa_components(),
                                    [&_in](Py_ssize_t i) {
                 return _in->get_fa_components(static_cast<DCProgs::t_int>(i)); } );
    }
 
  } %}

  // Add some printout functions.
  %extend {
    PyObject* __str__() {
      std::ostringstream sstr; sstr << *($self);
      return PyString_FromString(sstr.str().c_str());
    }
    PyObject* __repr__() {
      std::ostringstream sstr; sstr << *($self);
      return PyString_FromString(sstr.str().c_str());
    }
  }
};
%clear ApproxSurvivor;

