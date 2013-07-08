%feature("autodoc", "Ideal Likelihood.\n\n"
         "This object can be instantiated one of several way:\n\n"
         "   - With a matrix and an integer\n\n"
         "     >>> idealg = IdealG(array([...]), 2)\n\n"
         "   - With a StateMatrix\n\n"
         "     >>> matrix = StateMatrix(array([...]), 2)\n"
         "     >>> idealg = IdealG(matrix)\n\n") IdealG;
%feature("autodoc", "AF transitions with respect to time: :math:`e^{tQ_{FF}}Q_{FA}`.") IdealG::af;
%feature("autodoc", "FA transitions with respect to time: :math:`e^{tQ_{AA}}Q_{AF}$.") IdealG::fa;
%feature("autodoc", "AF transitions with respect to scale: :math:`(sI - Q_{AA})^{-1}Q_{AF}`.")
        IdealG::laplace_af;
%feature("autodoc", "FA transitions with respect to scale: :math:`(sI - Q_{FF})^{-1}Q_{FA}$`.")
        IdealG::laplace_fa;
// Now we define the bindings for the StateMatrix
// \brief State matrix that can  be partitioned into open/shut states.
// \details In practice, this is a two tuple with some helper functions to get corners.
class IdealG {
  public:

  // Some constructors.
  %extend { 
    IdealG(DCProgs::StateMatrix const &states) {
      return new DCProgs::IdealG(states.matrix, states.nopen); 
    }
    IdealG(DCProgs::t_rmatrix matrix, DCProgs::t_int nopen) {
      return new DCProgs::IdealG(matrix, nopen); 
    }
  }

  %extend {
    DCProgs::t_initvec _equilibrium_a() { return DCProgs::equilibrium(*$self); }
    DCProgs::t_initvec _equilibrium_f() { return DCProgs::equilibrium(*$self, false); }
  }
  %pythoncode %{
    __swig_getmethods__["equilibrium_a"] = _equilibrium_a
    __swig_getmethods__["equilibrium_f"] = _equilibrium_f
    if _newclass:
      equilibrium_af = property(_equilibrium_a, doc="Equilibrium for ideal AF process.")
      equilibrium_fa = property(_equilibrium_f, doc="Equilibrium for ideal AF process.")
  %}

  DCProgs::t_rmatrix af(DCProgs::t_real t) const;
  DCProgs::t_rmatrix fa(DCProgs::t_real t) const;
  DCProgs::t_rmatrix laplace_fa(DCProgs::t_real s) const;
  DCProgs::t_rmatrix laplace_af(DCProgs::t_real s) const;
};
