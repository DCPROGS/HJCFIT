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

%feature("autodoc", "Computes exact missed-event survivor function.\n\n"
         ":param qmatrix:\n"
         "   The Q-matrix for which to compute the exact likelihood. \n"
         ":param float double:\n"
         "   The maximum length of missed events.\n" ) ExactSurvivor;
%feature("autodoc", "Open to close transitions") ExactSurvivor::af;
%feature("autodoc", "Closed to open transitions") ExactSurvivor::fa;
%feature("autodoc", "Probability of no shut times detected between 0 and t.") ExactSurvivor::R_af;
%feature("autodoc", "Probability of no open times detected between 0 and t.") ExactSurvivor::R_fa;
%feature("autodoc", ":math:`C_{iml}` recursion matrices for af.") ExactSurvivor::recursion_af;
%feature("autodoc", ":math:`C_{iml}` recursion matrices for fa.") ExactSurvivor::recursion_fa;
%feature("autodoc", ":math:`D_i` matrices for af.") ExactSurvivor::D_af;
%feature("autodoc", ":math:`D_i` matrices for fa.") ExactSurvivor::D_fa;
%feature("autodoc", "Eigenvalues in of the transition matrix.\n\n"
         "This is the real part of the eigenvalues only.\n") ExactSurvivor::eigenvalues_af;
%feature("autodoc", "Eigenvalues in of the transpose transition matrix.\n\n"
         "This is the real part of the eigenvalues only.\n"
         "For numerical reasons, the order of the eigenvalues for AF and FA recursion\n"
         "may be different") ExactSurvivor::eigenvalues_fa;
%feature("autodoc", "Resolution or maximum length of the missed events") ExactSurvivor::tau;

class ExactSurvivor {
  public:

  // Constructor.
  ExactSurvivor(DCProgs::QMatrix const &_matrix, DCProgs::t_real _tau);

  %typemap(check) DCProgs::t_real t {
    if($1 > 11e0 * arg1->get_tau()) {
      std::ostringstream sstr;
      sstr << "Input value for time (" << $1
           << ") is more than 11 times larger that tau (" << arg1->get_tau() << "). "
           << "Could take a very long-time... Aborting\n"
              "If you want to disable this safeguard, take a look at"
              " exact_survivor.swg, modify and recompile.\n";
      PyErr_SetString(PyExc_ValueError, sstr.str().c_str());
      return NULL;
    }
  }
  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);
  %clear DCProgs::t_real t;
  %extend {
    PyObject* af(PyObject *t) {
      return apply_numpy(t, [$self](DCProgs::t_real _t) { return $self->af(_t); }); 
    }
    PyObject* fa(PyObject *t) {
      return apply_numpy(t, [$self](DCProgs::t_real _t) { return $self->fa(_t); }); 
    }
  }

  DCProgs::t_rmatrix recursion_af(DCProgs::t_uint i, DCProgs::t_uint m, DCProgs::t_uint l) const;
  DCProgs::t_rmatrix recursion_fa(DCProgs::t_uint i, DCProgs::t_uint m, DCProgs::t_uint l) const;
  DCProgs::t_rmatrix D_af(DCProgs::t_uint _i) const;
  DCProgs::t_rmatrix D_fa(DCProgs::t_uint _i) const;

  // Add extensions to compute properties
  %extend {
    DCProgs::t_real const tau;
    DCProgs::t_rvector const eigenvalues_af;
    DCProgs::t_rvector const eigenvalues_fa;
  }
  %{
    DCProgs::t_real DCProgs_ExactSurvivor_tau_get(DCProgs::ExactSurvivor *_in) {
      return _in->get_tau();
    }
    DCProgs::t_rvector const DCProgs_ExactSurvivor_eigenvalues_af_get(DCProgs::ExactSurvivor *_in) {
      return _in->eigenvalues_af();
    }
    DCProgs::t_rvector const DCProgs_ExactSurvivor_eigenvalues_fa_get(DCProgs::ExactSurvivor *_in) {
      return _in->eigenvalues_fa();
    }
  %}

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

%clear ExactSurvivor;
