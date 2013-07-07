%feature("autodoc", "Computes exact missed-event survivor function.") ExactSurvivor;
%feature("autodoc", "Open to close transitions") ExactSurvivor::af;
%feature("autodoc", "Closed to open transitions") ExactSurvivor::fa;
%feature("autodoc", "Probability of no shut times detected between 0 and t.") ExactSurvivor::R_af;
%feature("autodoc", "Probability of no open times detected between 0 and t.") ExactSurvivor::R_fa;
%feature("autodoc", ":math:`C_{iml}` recursion matrices for af.") ExactSurvivor::recursion_af;
%feature("autodoc", ":math:`C_{iml}` recursion matrices for fa.") ExactSurvivor::recursion_fa;
%feature("autodoc", ":math:`D_i` matrices for af.") ExactSurvivor::D_af;
%feature("autodoc", ":math:`D_i` matrices for fa.") ExactSurvivor::D_fa;
%feature("autodoc", "Eigenvalues in af basis.") ExactSurvivor::eigenvalues_af;
%feature("autodoc", "Eigenvalues in fa basis.") ExactSurvivor::eigenvalues_fa;
class ExactSurvivor {
  public:

  // Constructor.
  ExactSurvivor(DCProgs::StateMatrix const &_matrix, DCProgs::t_real _tau);

  %typemap(check) DCProgs::t_real t {
    if($1 > 11e0 * arg1->get_tau()) {
      std::ostringstream sstr;
      sstr << "Input value for time (" << $1
           << ") is more than 11 times larger that tau (" << arg1->get_tau() 
           << ". Could take a very long-time... Aborting\n"
              "If you want to disable this safeguard, take a look at"
              " exact_survivor.swg, modify and recompile.\n";
      PyErr_SetString(PyExc_ValueError, sstr.str().c_str());
      return NULL;
    }
  }
  DCProgs::t_rmatrix af(DCProgs::t_real t);
  DCProgs::t_rmatrix fa(DCProgs::t_real t);

  %typemap(check) DCProgs::t_real t {
    if($1 > 10e0 * arg1->get_tau()) {
      std::ostringstream sstr;
      sstr << "Input value for time (" << $1
           << ") is more than 10 times larger that tau (" << arg1->get_tau() 
           << ". Could take a very long-time... Aborting\n"
              "If you want to disable this safeguard, take a look at"
              " exact_survivor.swg, modify and recompile.\n";
      PyErr_SetString(PyExc_ValueError, sstr.str().c_str());
      return NULL;
    }
  } 
  %clear DCProgs::t_real t;

  DCProgs::t_rmatrix recursion_af(DCProgs::t_int i, DCProgs::t_int m, DCProgs::t_int l) const;
  DCProgs::t_rmatrix recursion_fa(DCProgs::t_int i, DCProgs::t_int m, DCProgs::t_int l) const;
  DCProgs::t_rmatrix D_af(DCProgs::t_int _i) const;
  DCProgs::t_rmatrix D_fa(DCProgs::t_int _i) const;

  // Renames some values so we can make it a property.
  %rename(_get_tau) get_tau() const; 
  %rename(_eigenvalues_af) eigenvalues_af() const;
  %rename(_eigenvalues_fa) eigenvalues_fa() const;
  DCProgs::t_real get_tau() const;
  DCProgs::t_rvector eigenvalues_af() const;
  DCProgs::t_rvector eigenvalues_fa() const;
  %pythoncode %{
    __swig_getmethods__["tau"] = _get_tau
    __swig_getmethods__["eigenvalues_af"] = _eigenvalues_af
    __swig_getmethods__["eigenvalues_fa"] = _eigenvalues_fa
    if _newclass:
      tau = property(_get_tau, doc="Max length of mixed events.")
      eigenvalues_af = property(_eigenvalues_af,
                                doc="Eigenvalues of the transition matrix.\n\n"  \
                                    "This is the real part of eigenvalues only.\n")
      eigenvalues_fa = property(_eigenvalues_fa,
                                doc="Eigenvalues of the transition matrix.\n\n"  \
                                    "This is the real part of eigenvalues only.\n")
  %}
};

%clear ExactSurvivor;
