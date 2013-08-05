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

%feature("autodoc", "Computes missed-events likelihood.\n\n"
         "Exact calculations take place for times smaller than `nmax * tau`.\n"
         "Asymptotic calculations take over for larger times.\n\n"
         "There are two possible ways to instantiate this object:\n\n"
         ">>> MissedEvents(determinant_af, roots_af, determinant_fa, roots_fa[, nmax])\n"
         ">>> MissedEvents(qmatrix, tau, [nmax, xtol, rtol, itermax])\n\n"
         "The parameters bewteen brackets are optional. The last two versions will try\n"
         "and calculate the roots of the determinantal equations automatically.\n\n"
         ":param determinant_af:\n"
         "   A :class:`DeterminantEq` instance, specifically for the af block.\n"
         ":param roots_af:\n"
         "  The roots of the af determinantal equation. The should come in the format "
            "`[(root, multiplicity), (root, multiplicity), ...]`.\n" 
         ":param determinant_fa:\n"
         "   A :class:`DeterminantEq` instance, specifically for the fa block. "
            "It should the transpose of `determinant_af`. It is required so "
            "it need not be recomputed, since it most likely already exists.\n"
         ":param roots_fa:\n"
         "  The roots of the fa determinantal equation. The should come in the format "
            "`[(root, multiplicity), (root, multiplicity), ...]`.\n" 
         ":param qmatrix: \n"
         "   A :class:`QMatrix` instance. \n"
         ":param int nmax: \n"
         "   The exact missed event likelihood will be computed for times "
            ":math:`t \\in [0, n_{\\mathrm{max}} \\tau`. It defaults to 2.\n"
         ":param float xtol: \n"
         "   Tolerance criteria when computing roots using "
            ":py:func:`dcprogs.likelihood.brentq`. Defaults to 1e-10.\n"
         ":param float rtol: \n"
         "   Tolerance criteria when computing roots using "
            ":py:func:`dcprogs.likelihood.brentq`. Defaults to 1e-10\n"
         ":param float itermax: \n"
         "   Maximum number of iterations when computing roots using "
            ":py:func:`dcprogs.likelihood.brentq`. Defaults to 100.\n"
        ) MissedEventsG;


%feature("autodoc", "Open to close transitions\n\n") MissedEventsG::af(PyObject*);
%feature("autodoc", "Closed to open transitions") MissedEventsG::fa(PyObject*);
%feature("autodoc", "Exact missed-events G in Laplace space.\n\n"
         "The exact expression is "
         ":math:`^{e}G_{AF}(s) = {}^AR(s) e^{-s\\\\tau}Q_{AF}e^{Q_{FF}\\\\tau}`, "
         "with :math:`{}^AR(s) = [sI - Q_{AA} - Q_{AF} \\\\int_0^\\\\tau e^{-st}"
         "e^{Q_{FF}t}\\\\partial t Q_{FA}]^{-1}`.\n\n"  
         ":param s:\n"
         "  The laplace scale. A real scalar or something convertible to a numpy array.\n\n"
         ":returns:\n"
         "  A matrix if the input is scalar, an array of matrices otherwise, \n"
         "  with the shape of the input.") MissedEventsG::laplace_af(PyObject*);
%feature("autodoc", "Exact missed-events G in Laplace space.\n\n"
         "The exact expression is "
         ":math:`^{e}G_{FA}(s) = {}^FR(s) e^{-s\\\\tau}Q_{FA}e^{Q_{AA}\\\\tau}`, "
         "with :math:`{}^FR(s) = [sI - Q_{FF} - Q_{FA} \\\\int_0^\\\\tau e^{-st}"
         "e^{Q_{AA}t}\\\\partial t Q_{AF}]^{-1}`.\n\n"
         ":param s:\n"
         "  The laplace scale. A real scalar or something convertible to a numpy array.\n\n"
         ":returns:\n"
         "  A matrix if the input is scalar, an array of matrices otherwise, \n"
         "  with the shape of the input.") MissedEventsG::laplace_fa(PyObject*);
%feature("autodoc", "Computes initial CHS occupancies.") MissedEventsG::initial_CHS_occupancies;
%feature("autodoc", "Computes final CHS occupancies.") MissedEventsG::final_CHS_occupancies;

class MissedEventsG {
  public:

  %extend {
    MissedEventsG(DCProgs::DeterminantEq *_af, PyObject * _roots_af, 
                  DCProgs::DeterminantEq *_fa, PyObject * _roots_fa,
                  DCProgs::t_uint _nmax=2 ) {
      using namespace DCProgs;
      if(_nmax < 0) throw errors::Domain("nmax cannot be negative.");
      if(_nmax > 10) throw errors::Domain("nmax cannot be larger than 10 (that's just unwise).");
      std::vector<Root> roots_af = convert_to_roots(_roots_af);
      std::vector<Root> roots_fa = convert_to_roots(_roots_fa);
      if(roots_af.size() == 0) {
        std::ostringstream sstr;
        sstr << "No AF roots given  on input.\n\n" << *_af << "\n";
        throw errors::Domain(sstr.str());
      }
      if(roots_fa.size() == 0) {
        std::ostringstream sstr; 
        sstr << "No FA roots given  on input.\n\n" << *_fa << "\n";
        throw errors::Domain(sstr.str());
      }
      return new DCProgs::MissedEventsG(*_af, roots_af, *_fa, roots_fa, _nmax);
    }
  }
  MissedEventsG( QMatrix const & _qmatrix, DCProgs::t_real _tau, 
                 DCProgs::t_uint _nmax=2, DCProgs::t_real _xtol = 1e-12,
                 DCProgs::t_real _rtol = 1e-12, DCProgs::t_uint _itermax = 100 );

  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);
  %extend {
    PyObject* af(PyObject *t) {
      return apply_numpy(t, [$self](DCProgs::t_real t) { return $self->af(t); }); 
    }
    PyObject* fa(PyObject *t) {
      return apply_numpy(t, [$self](DCProgs::t_real t) { return $self->fa(t); }); 
    }
  }
  DCProgs::t_rmatrix laplace_af(DCProgs::t_real s);
  DCProgs::t_rmatrix laplace_fa(DCProgs::t_real s);
  %extend {
    PyObject* laplace_af(PyObject *t) {
      return apply_numpy(t, [$self](DCProgs::t_real s) { return $self->laplace_af(s); }); 
    }
    PyObject* laplace_fa(PyObject *t) {
      return apply_numpy(t, [$self](DCProgs::t_real s) { return $self->laplace_fa(s); }); 
    }
  }
  // Make CHS_occupancies member methods of MissedEventsG. 
  // Python is a glue language, so it makes sense to keep both MissedEventsG functionality and
  // CHS vector functionality together, even when it doesn't in c++.
  // Basically, c++ stuff is kept conceptual -- and separate concepts are kept separate.
  // On the other hand, python is just about being practical.
  %extend {
    DCProgs::t_initvec initial_CHS_occupancies(DCProgs::t_real t_critical) {
      return DCProgs::CHS_occupancies(*$self, t_critical, true); 
    }
    DCProgs::t_rvector final_CHS_occupancies(DCProgs::t_real t_critical) { 
      return DCProgs::CHS_occupancies(*$self, t_critical, false); 
    }
  }

  // Renames some values so we can make them properties
  %rename(_get_tau) get_tau() const; 
  %rename(_get_af_factor) get_af_factor() const; 
  %rename(_get_fa_factor) get_fa_factor() const; 
  DCProgs::t_real get_tau() const;
  DCProgs::t_rmatrix get_af_factor() const;
  DCProgs::t_rmatrix get_fa_factor() const;
  // Add extensions to compute occupancies.
  %extend {
    DCProgs::t_initvec _initial_occupancies() { return DCProgs::occupancies(*$self); }
    DCProgs::t_rvector _final_occupancies() { return DCProgs::occupancies(*$self, false); }
  }
  // Now go to python and add the properties.
  %pythoncode %{
    __swig_getmethods__["tau"] = _get_tau
    __swig_getmethods__["af_factor"] = _get_af_factor
    __swig_getmethods__["fa_factor"] = _get_fa_factor
    __swig_getmethods__["initial_occupancies"] = _initial_occupancies
    __swig_getmethods__["final_occupancies"]   = _final_occupancies
    if _newclass:
      tau = property(_get_tau, doc="Max length of mixed events.")
      af_factor = property(_get_af_factor, doc=":math:`Q_{AF}e^{-Q_{FF}\\tau}`")
      fa_factor = property(_get_fa_factor, doc=":math:`Q_{FA}e^{-Q_{AA}\\tau}`")
      initial_occupancies = property(_initial_occupancies, doc="Occupancies for initial states.")
      final_occupancies = property(_final_occupancies, doc="Occupancies for final states.")
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

  // Some more properties
  %extend { 
    DCProgs::t_real const tmax; 
    DCProgs::t_uint nmax; 
    DCProgs::t_uint const nopen;
    DCProgs::t_uint const nshut;
  }
  %{ namespace {
    DCProgs::t_real DCProgs_MissedEventsG_tmax_get(DCProgs::MissedEventsG *_in) {
      return _in->get_tmax();
    }
    DCProgs::t_uint DCProgs_MissedEventsG_nmax_get(DCProgs::MissedEventsG *_in) {
      return _in->get_nmax();
    }
    void DCProgs_MissedEventsG_nmax_set(DCProgs::MissedEventsG *_in, DCProgs::t_uint _n) {
      _in->set_nmax(_n);
    }
    DCProgs::t_uint DCProgs_MissedEventsG_nopen_get(DCProgs::MissedEventsG *_in) {
      return _in->get_qmatrix().nopen;
    }
    DCProgs::t_uint DCProgs_MissedEventsG_nshut_get(DCProgs::MissedEventsG *_in) {
      return _in->get_qmatrix().nshut();
    }
  }%}
};

%clear MissedEventsG;

