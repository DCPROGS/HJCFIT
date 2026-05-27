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

%feature("autodoc", "Ideal Likelihood.\n\n"
         "This object can be instantiated one of several way:\n\n"
         "- With a matrix and an integer\n\n"
         "  >>> idealg = IdealG(array([...]), 2)\n\n"
         "- With a QMatrix\n\n"
         "  >>> matrix = QMatrix(array([...]), 2)\n"
         "  >>> idealg = IdealG(matrix)\n\n") IdealG;
%feature("autodoc", "AF transitions with respect to time\n\n"
         "Implements the ideal likelihood:\n\n"
         ".. math::\n\n  e^{t\\mathcal{Q}_{FF}}\\mathcal{Q}_{FA}.\n\n"
         ":param t:\n"
         "   A scalar or something to a numpy array. In the latter case, the return is an array\n"
         "   of matrices.\n") IdealG::af(PyObject* _s);
%feature("autodoc", "FA transitions with respect to time\n\n"
         "Implements the ideal likelihood:\n\n"
         ".. math::\n\n  e^{t\\mathcal{Q}_{AA}}\\mathcal{Q}_{AF}.\n\n"
         ":param t:\n"
         "   A scalar or something to a numpy array. In the latter case, the return is an array\n"
         "   of matrices.\n") IdealG::fa(PyObject* _s);
%feature("autodoc", "AF transitions with respect to scale\n\n" 
         "Implements the laplace transform of the likelihood:\n\n"
         ".. math::\n\n  (sI - \\mathcal{Q}_{AA})^{-1}\\mathcal{Q}_{AF}.\n\n"
         ":param s:\n"
         "   A scalar or something to a numpy array. In the latter case, the return is an array\n"
         "   of matrices.\n") IdealG::laplace_af(PyObject*);
%feature("autodoc", "FA transitions with respect to scale\n\n"
         "Implements the laplace transform of the likelihood:\n\n"
         ".. math::\n\n  (sI - \\mathcal{Q}_{FF})^{-1}\\mathcal{Q}_{FA}.\n\n"
         ":param s:\n"
         "   A scalar or something to a numpy array. In the latter case, the return is an array\n"
         "   of matrices.\n") IdealG::laplace_fa(PyObject*);
%feature("autodoc", "Equilibrium vectors for initial states.\n\n"
         "Computes the left eigenvector of :math:`\\mathcal{G}_{AF}\\mathcal{G}_{FA}`, where\n"
         ":math:`\\mathcal{G}_{AF}` is the laplacian for :math:`s=0` of the likelihood.")
         IdealG::initial_vectors;
%feature("autodoc", "Equilibrium vectors for final states.\n\n"
         "Computes the right eigenvector of :math:`\\mathcal{G}_{FA}\\mathcal{G}_{AF}`, where\n"
         ":math:`\\mathcal{G}_{FA}` is the laplacian for :math:`s=0` of the likelihood.")
         IdealG::final_vectors;
%feature("autodoc", "Number of open-states.") IdealG::nopen;
%feature("autodoc", "Number of shut-states.") IdealG::nshut;

// Now we define the bindings for the QMatrix
// \brief State matrix that can  be partitioned into open/shut states.
// \details In practice, this is a two tuple with some helper functions to get corners.
class IdealG {
  public:

  // Some constructors.
  %extend { 
    IdealG(PyObject* _py_matrix, HJCFIT::t_uint _nopen) {
      using namespace HJCFIT;
      if(_nopen == 0) 
        throw errors::PythonValueError("Number of open states cannot be negative or zero.");
      t_rmatrix const matrix = HJCFIT::numpy::map_to_rmatrix(_py_matrix);
      if(_nopen > static_cast<t_uint>(std::max(matrix.rows(), matrix.cols())) )
        throw errors::PythonValueError(
                "Number of open states cannot be larger than the number states.");
      return new IdealG(std::move(matrix), _nopen); 
    }
  }
  IdealG(HJCFIT::QMatrix const &states);

  // Adds some quantities as properties
  %extend {
    HJCFIT::t_initvec const initial_vectors;
    HJCFIT::t_rvector const final_vectors;
    HJCFIT::t_uint const nopen;
    HJCFIT::t_uint const nshut;
  }
  %{
    HJCFIT::t_initvec HJCFIT_IdealG_initial_vectors_get(HJCFIT::IdealG *_in) {
      return HJCFIT::vectors(*_in, true);
    }
    HJCFIT::t_rvector HJCFIT_IdealG_final_vectors_get(HJCFIT::IdealG *_in) {
      return HJCFIT::vectors(*_in, false).transpose();
    }
    HJCFIT::t_uint HJCFIT_IdealG_nopen_get(HJCFIT::IdealG *_in) { return _in->get_nopen(); }
    HJCFIT::t_uint HJCFIT_IdealG_nshut_get(HJCFIT::IdealG *_in) { return _in->get_nshut(); }
  %}

  HJCFIT::t_rmatrix af(HJCFIT::t_real t) const;
  HJCFIT::t_rmatrix fa(HJCFIT::t_real t) const;
  %extend {
    PyObject* af(PyObject *t) {
      auto lambda = [$self](HJCFIT::t_real _t) -> HJCFIT::t_rmatrix { return $self->af(_t); };
      return apply_numpy(t, lambda); 
    }
    PyObject* fa(PyObject *t) {
      auto lambda = [$self](HJCFIT::t_real _t) -> HJCFIT::t_rmatrix { return $self->fa(_t); };
      return apply_numpy(t, lambda); 
    }
  }
  HJCFIT::t_rmatrix laplace_fa(HJCFIT::t_real s) const;
  HJCFIT::t_rmatrix laplace_af(HJCFIT::t_real s) const;
  %extend {
    PyObject* laplace_af(PyObject *t) {
      auto lambda = [$self](HJCFIT::t_real _t) -> HJCFIT::t_rmatrix {
        return $self->laplace_af(_t); 
      };
      return apply_numpy(t, lambda); 
    }
    PyObject* laplace_fa(PyObject *t) {
      auto lambda = [$self](HJCFIT::t_real _t) -> HJCFIT::t_rmatrix {
        return $self->laplace_fa(_t); 
      };
      return apply_numpy(t, lambda); 
    }
  }

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
