%feature("autodoc", "Computes exact missed-event G function.") MissedEventsG;
%feature("autodoc", "Open to close transitions") MissedEventsG::af;
%feature("autodoc", "Closed to open transitions") MissedEventsG::fa;
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

  // Renames some values so we can make it a property.
  %rename(_get_tau) get_tau() const; 
  %rename(_get_af_factor) get_af_factor() const; 
  %rename(_get_fa_factor) get_fa_factor() const; 
  DCProgs::t_real get_tau() const;
  DCProgs::t_rmatrix get_af_factor() const;
  DCProgs::t_rmatrix get_fa_factor() const;
  %pythoncode %{
    __swig_getmethods__["tau"] = _get_tau
    __swig_getmethods__["af_factor"] = _get_af_factor
    __swig_getmethods__["fa_factor"] = _get_fa_factor
    if _newclass:
      tau = property(_get_tau, doc="Max length of mixed events.")
      af_factor = property(_get_af_factor, doc=":math:`Q_{AF}e^{-Q_{FF}\\tau}`")
      fa_factor = property(_get_fa_factor, doc=":math:`Q_{FA}e^{-Q_{AA}\\tau}`")
  %}
};

%clear MissedEventsG;

%pythoncode %{
  def create_missed_eventsG(rates, tau, nmax=2):

    """ Creates a MissedEvents function from knowledge of rate matrix. """
    determinant_af = DeterminantEq(rates, tau) 
    determinant_fa = DeterminantEq(rates.transpose(), tau) 
    roots_af = find_roots(determinant_af)
    roots_fa = find_roots(determinant_fa)
    return MissedEventsG(determinant_af, roots_af, determinant_fa, roots_fa, nmax)
%}
