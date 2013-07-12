%feature("autodoc", "Computes missed-events likelihood.\n\n"
         "Exact calculations take place for times smaller than `nmax * tau`. "
         "Asymptotic calculations take over for larger times." ) MissedEventsG;
%feature("autodoc", "Open to close transitions\n\n") MissedEventsG::af;
%feature("autodoc", "Closed to open transitions") MissedEventsG::fa;
%feature("autodoc", "Exact missed-events G in Laplace space.\n\\"
         "The exact expression is "
         ":math:`^{e}G_{AF}(s) = {}^AR(s) e^{-s\\tau}Q_{AF}e^{Q_{FF}\\tau}`, "
         "with :math:{}^AR(s) = [sI - Q_{AA} - Q_{AF} \\int_0^\\taue^{-st}"
         "e^{Q_{FF}t}\\partial t Q_{FA}]^{-1}`.") MissedEventsG::laplace_af;
%feature("autodoc", "Exact missed-events G in Laplace space.\n\\"
         "The exact expression is "
         ":math:`^{e}G_{FA}(s) = {}^FR(s) e^{-s\\tau}Q_{FA}e^{Q_{AA}\\tau}`, "
         "with :math:{}^FR(s) = [sI - Q_{FF} - Q_{FA} \\int_0^\\taue^{-st}"
         "e^{Q_{AA}t}\\partial t Q_{AF}]^{-1}`.") MissedEventsG::laplace_af;
%feature("autodoc", "Computes initial CHS occupancies.") MissedEventsG::_initial_CHS_occupancies;
%feature("autodoc", "Computes final CHS occupancies.") MissedEventsG::_final_CHS_occupancies;
class MissedEventsG {
  public:

  %extend {
    MissedEventsG(DCProgs::DeterminantEq *_af, PyObject * _roots_af, 
                  DCProgs::DeterminantEq *_fa, PyObject * _roots_fa,
                  DCProgs::t_int _nmax=2 ) {
      if(_nmax < 0) throw DCProgs::errors::Domain("nmax cannot be negative.");
      if(_nmax > 10)
        throw DCProgs::errors::Domain("nmax cannot be larger than 10 (that's just unwise).");
      std::vector<DCProgs::Root> roots_af = convert_to_roots(_roots_af);
      std::vector<DCProgs::Root> roots_fa = convert_to_roots(_roots_fa);
      return new DCProgs::MissedEventsG(*_af, roots_af, *_fa, roots_fa, _nmax);
    }
  }

  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);
  DCProgs::t_rmatrix laplace_af(DCProgs::t_real t);
  DCProgs::t_rmatrix laplace_fa(DCProgs::t_real t);
  %extend {
    DCProgs::t_initvec _initial_CHS_occupancies(DCProgs::t_real t_critical) {
      return DCProgs::CHS_occupancies(*$self, $self->get_tau(), t_critical); 
    }
    DCProgs::t_initvec _final_CHS_occupancies(DCProgs::t_real t_critical) { 
      return DCProgs::CHS_occupancies(*$self, $self->get_tau(), t_critical); 
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
    DCProgs::t_initvec _final_occupancies() { return DCProgs::occupancies(*$self, false); }
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
  %extend { DCProgs::t_real const tmax; }
  %{
    DCProgs::t_real DCProgs_MissedEventsG_tmax_get(DCProgs::MissedEventsG *_in) {
      return _in->get_tmax();
    }
  %}
};

%clear MissedEventsG;

%pythoncode %{
  def create_missed_eventsG(qmatrix, tau, nmax=2):

    """ Creates a MissedEvents function from knowledge of rate matrix. """
    determinant_af = DeterminantEq(qmatrix, tau) 
    determinant_fa = determinant_af.transpose()
    roots_af = find_roots(determinant_af)
    roots_fa = find_roots(determinant_fa)
    return MissedEventsG(determinant_af, roots_af, determinant_fa, roots_fa, nmax)
%}
