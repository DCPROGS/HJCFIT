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

%feature("autodoc", "Computes missed-events likelihood.\n\n"
         "Exact calculations take place for times smaller than :math:`n_{\\mathrm{max}}\\\\tau`. "
         "Asymptotic calculations take over for larger times.\n\n") MissedEventsG;


%feature("autodoc", "Likelihood of an observed open time of length ``t``\n\n"
         ":param t:\n"
         "   A scalar or something to a numpy array. In the latter case, the return is an array\n"
         "   of matrices.\n") MissedEventsG::af(PyObject*);
%feature("autodoc", "Likelihood of a shut time of length ``t``\n\n"
         ":param t:\n"
         "   A scalar or something to a numpy array. In the latter case, the return is an array\n"
         "   of matrices.\n") MissedEventsG::fa(PyObject*);
%feature("autodoc", "Exact missed-events G in Laplace space.\n\n"
         "The exact expression is "
         ":math:`^{e}\\mathcal{G}_{AF}(s) = {}^AR(s) e^{-s\\\\tau}\\mathcal{Q}_{AF}"
                 "e^{\\mathcal{Q}_{FF}\\\\tau}`, "
         "with :math:`{}^AR(s) = [sI - \\mathcal{Q}_{AA} - \\mathcal{Q}_{AF}"
                      " \\\\int_0^\\\\tau e^{-st}"
         "e^{\\mathcal{Q}_{FF}t}\\\\partial t \\mathcal{Q}_{FA}]^{-1}`.\n\n"  
         ":param s:\n"
         "  The laplace scale. A real scalar or something convertible to a numpy array.\n\n"
         ":returns:\n"
         "  A matrix if the input is scalar, an array of matrices otherwise, \n"
         "  with the shape of the input.") MissedEventsG::laplace_af(PyObject*);
%feature("autodoc", "Exact missed-events G in Laplace space.\n\n"
         "The exact expression is "
         ":math:`^{e}\\mathcal{G}_{FA}(s) = {}^FR(s) e^{-s\\\\tau}\\mathcal{Q}_{FA}"
                   "e^{\\mathcal{Q}_{AA}\\\\tau}`, "
         "with :math:`{}^FR(s) = [sI - \\mathcal{Q}_{FF} - \\mathcal{Q}_{FA} "
                     "\\\\int_0^\\\\tau e^{-st}"
         "e^{\\mathcal{Q}_{AA}t}\\\\partial t \\mathcal{Q}_{AF}]^{-1}`.\n\n"
         ":param s:\n"
         "  The laplace scale. A real scalar or something convertible to a numpy array.\n\n"
         ":returns:\n"
         "  A matrix if the input is scalar, an array of matrices otherwise, \n"
         "  with the shape of the input.") MissedEventsG::laplace_fa(PyObject*);
%feature("autodoc", "Computes initial CHS vectors.") MissedEventsG::initial_CHS_vectors;
%feature("autodoc", "Computes final CHS vectors.") MissedEventsG::final_CHS_vectors;
%feature("autodoc", "Resolution or maximum length of the missed events.") MissedEventsG::tau;
%feature("autodoc", "Factor accounting for minimum shut time\n\n"
         "It is the likelihood :math:`\\mathcal{Q}_{AF}e^{-\\mathcal{Q}_{FF}\\tau}` of a shut\n"
         "time of length :math:`\\tau`.") MissedEventsG::af_factor;
%feature("autodoc", "Factor accounting for minimum open time\n\n"
         "It is the likelihood :math:`\\mathcal{Q}_{FA}e^{-\\mathcal{Q}_{AA}\\tau}` of an open\n"
         "time of length :math:`\\tau`.") MissedEventsG::fa_factor;
%feature("autodoc", "Equilibrium vectors for initial states.\n\n"
         "Computes the left eigenvector of :math:`{}^e\\mathcal{G}_{AF}{}^e\\mathcal{G}_{FA}`, \n"
         "where :math:`{}^e\\mathcal{G}_{AF}` is the laplacian for :math:`s=0` of the likelihood.")
         MissedEventsG::initial_vectors;
%feature("autodoc", "Equilibrium vectors for final states.\n\n"
         "Computes the right eigenvector of :math:`{}^e\\mathcal{G}_{FA}{}^e\\mathcal{G}_{AF}`, \n"
         "where :math:`{}^e\\mathcal{G}_{FA}` is the laplacian for :math:`s=0` of the likelihood.")
         MissedEventsG::final_vectors;
%feature("autodoc", "Cut-off time of exact calculations "
         ":math:`t_{\\mathrm{max}} = (n_{\\mathrm{max}}-1)\\\\tau`.\n\n"
         "For practical reasons, the minimum observation time has already"
         "been removed here.") MissedEventsG::tmax;
%feature("autodoc", "Cut-off time of exact calculations in units of :math:`\\\\tau`.")
         MissedEventsG::nmax;
%feature("autodoc", "Number of open-states.") MissedEventsG::nopen;
%feature("autodoc", "Number of shut-states.") MissedEventsG::nshut;

%feature("shadow") MissedEventsG::MissedEventsG %{
  def __init__(self, *args, **kwargs):
    """ Initializes the missed-events likelihood. 

        Exact calculations take place for times smaller than `nmax * tau`.
        Asymptotic calculations take over for larger times.

        There are three possible ways to instantiate this object:

        >>> MissedEvents(determinant_af, roots_af, determinant_fa, roots_fa[, nmax=2])
        >>> MissedEvents(matrix, nopen, tau, **kwargs)
        >>> MissedEvents(qmatrix, tau, **kwargs)

        The parameters between brackets are optional. The last two versions will try 
        and calculate the roots of the determinant equations automatically. A number of parameters
        can be given to control this process.

        :param determinant_af:
            A :class:`DeterminantEq` instance, specifically for the af block.
        :param roots_af:
           The roots of the af determinant equation. The should come in the format
           `[(root, multiplicity), (root, multiplicity), ...]`.
        :param determinant_fa:
           A :class:`DeterminantEq` instance, specifically for the fa block.
           It should the transpose of `determinant_af`. It is required so 
           it need not be recomputed, since it most likely already exists.
        :param roots_fa:
           The roots of the fa determinant equation. The should come in the format 
           `[(root, multiplicity), (root, multiplicity), ...]`.
        :param matrix:
           An object convertible to a square matrix.
        :param integer nopen:
           Number of open states.
        :param qmatrix:
           A :class:`QMatrix` instance.
        :param int nmax: 
           The exact missed event likelihood will be computed for times 
           :math:`t \in [0, n_{\mathrm{max}} \tau`. Defaults to 3.
        :param float xtol: 
           Tolerance criteria when computing roots using :func:`~HJCFIT.likelihood.brentq`.
           Defaults to 1e-12.
        :param float rtol:
           Tolerance criteria when computing roots using :func:`~HJCFIT.likelihood.brentq`.
           Defaults to 1e-12.
        :param float itermax:
           Maximum number of iterations when computing roots using
           :func:`~HJCFIT.likelihood.brentq`.
           Defaults to 100.
        :param float lower_bound:
           Lower bound for *all* roots. Defaults to None, in which the case the lower bound is
           computed from :py:func:`find_lower_bound_for_roots`.
        :param float upper_bound:
           Upper bound for *all* roots. Defaults to None, in which the case the upper bound is
           computed from :py:func:`find_upper_bound_for_roots`.
    """
    from numpy import NaN, array
    from .. import internal_dtype
    from . import QMatrix
    # Try first initalization mode
    det_af   = args[0] if len(args) > 0 else kwargs.get("det_af", None)
    if isinstance(det_af, DeterminantEq):
      roots_af = args[1] if len(args) > 1 else kwargs.get("roots_af", None)
      if roots_af is None: raise TypeError('Missing roots_af argument')
      det_fa   = args[2] if len(args) > 2 else kwargs.get("det_fa", None)
      if det_fa is None: raise TypeError('Missing det_fa argument')
      roots_fa = args[3] if len(args) > 3 else kwargs.get("roots_fa", None)
      if roots_fa is None: raise TypeError('Missing roots_fa argument')
      nmax = int(args[4] if len(args) > 4 else kwargs.get('nmax', 3))
      this = $action(det_af, roots_af, det_fa, roots_fa, nmax)
      try: self.this.append(this)
      except: self.this = this
      return

    # Try second and third instanciation
    matrix = args[0] if len(args) > 0 else kwargs.get('matrix', kwargs.get('qmatrix', None))
    if matrix is None: 
      raise ValueError("Expected the q-matrix or a matrix + number of open states on input.")
    index = 1
    if not isinstance(matrix, QMatrix):
      matrix = array(matrix, internal_dtype)
      nopen = int(args[index]) if len(args) > index else kwargs.get('nopen', None)
      if nopen is None: 
        raise ValueError("A matrix was given on input, but without a nopen parameter.")
      matrix = QMatrix(matrix, nopen)
      index += 1

    tau = args[index] if len(args) > index else kwargs.get('tau', None)
    if tau is None: raise ValueError("No tau given on input.")
    nmax        = int(args[index+1] if len(args) > index+1 else kwargs.get('nmax', 3))
    xtol        = args[index+2] if len(args) > index+2 else kwargs.get('xtol', 1e-12)
    rtol        = args[index+3] if len(args) > index+3 else kwargs.get('rtol', 1e-12)
    itermax     = int(args[index+4] if len(args) > index+4 else kwargs.get('itermax', 100))
    lower_bound = args[index+5] if len(args) > index+5 else kwargs.get('lower_bound', None)
    upper_bound = args[index+6] if len(args) > index+6 else kwargs.get('upper_bound', None)
    if lower_bound is None: lower_bound = NaN
    if upper_bound is None: upper_bound = NaN
    
    this = $action(matrix, tau, nmax, xtol, rtol, itermax, lower_bound, upper_bound)
    try: self.this.append(this)
    except: self.this = this
%}

class MissedEventsG {
  public:

  %extend {
    MissedEventsG(HJCFIT::DeterminantEq *_af, PyObject * _roots_af, 
                  HJCFIT::DeterminantEq *_fa, PyObject * _roots_fa,
                  HJCFIT::t_uint _nmax=3 ) {
      using namespace HJCFIT;
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
      return new HJCFIT::MissedEventsG(*_af, roots_af, *_fa, roots_fa, _nmax);
    }
  }
  MissedEventsG( QMatrix const & _qmatrix, HJCFIT::t_real _tau, 
                 HJCFIT::t_uint _nmax, HJCFIT::t_real _xtol,
                 HJCFIT::t_real _rtol, HJCFIT::t_uint _itermax,
                 HJCFIT::t_real _lowerbound, HJCFIT::t_real _upperbound );

  HJCFIT::t_rmatrix af(HJCFIT::t_real t);
  HJCFIT::t_rmatrix fa(HJCFIT::t_real t);
  %extend {
    PyObject* af(PyObject *t) {
      return apply_numpy(t, [$self](HJCFIT::t_real t) { return $self->af(t); }); 
    }
    PyObject* fa(PyObject *t) {
      return apply_numpy(t, [$self](HJCFIT::t_real t) { return $self->fa(t); }); 
    }
  }
  HJCFIT::t_rmatrix laplace_af(HJCFIT::t_real s);
  HJCFIT::t_rmatrix laplace_fa(HJCFIT::t_real s);
  %extend {
    PyObject* laplace_af(PyObject *t) {
      return apply_numpy(t, [$self](HJCFIT::t_real s) { return $self->laplace_af(s); }); 
    }
    PyObject* laplace_fa(PyObject *t) {
      return apply_numpy(t, [$self](HJCFIT::t_real s) { return $self->laplace_fa(s); }); 
    }
  }
  // Make CHS_vectors member methods of MissedEventsG. 
  // Python is a glue language, so it makes sense to keep both MissedEventsG functionality and
  // CHS vector functionality together, even when it doesn't in c++.
  // Basically, c++ stuff is kept conceptual -- and separate concepts are kept separate.
  // On the other hand, python is just about being practical.
  %extend {
    HJCFIT::t_initvec initial_CHS_vectors(HJCFIT::t_real t_critical) {
      return HJCFIT::CHS_vectors(*$self, t_critical, true); 
    }
    HJCFIT::t_rvector final_CHS_vectors(HJCFIT::t_real t_critical) { 
      return HJCFIT::CHS_vectors(*$self, t_critical, false); 
    }
  }

  // Add extensions to compute properties
  %extend {
    HJCFIT::t_real const tau;
    PyObject * const af_factor;
    PyObject * const fa_factor;
    HJCFIT::t_initvec const initial_vectors;
    HJCFIT::t_rvector const final_vectors;
  }
  %{
    HJCFIT::t_real HJCFIT_MissedEventsG_tau_get(HJCFIT::MissedEventsG *_in) {
      return _in->get_tau();
    }
    PyObject* HJCFIT_MissedEventsG_af_factor_get(HJCFIT::MissedEventsG *_in) {
      PyObject* result = HJCFIT::numpy::wrap_to_numpy(_in->get_af_factor());
      PyArray_CLEARFLAGS((PyArrayObject*)result, NPY_ARRAY_WRITEABLE);
      return result;
    }
    PyObject* HJCFIT_MissedEventsG_fa_factor_get(HJCFIT::MissedEventsG *_in) {
      PyObject* result = HJCFIT::numpy::wrap_to_numpy(_in->get_fa_factor());
      PyArray_CLEARFLAGS((PyArrayObject*)result, NPY_ARRAY_WRITEABLE);
      return result;
    }
    HJCFIT::t_initvec HJCFIT_MissedEventsG_initial_vectors_get(HJCFIT::MissedEventsG *_in) {
      return HJCFIT::vectors(*_in, true);
    }
    HJCFIT::t_rvector HJCFIT_MissedEventsG_final_vectors_get(HJCFIT::MissedEventsG *_in) {
      return HJCFIT::vectors(*_in, false).transpose();
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

  // Some more properties
  %extend { 
    HJCFIT::t_real const tmax; 
    HJCFIT::t_uint nmax; 
    HJCFIT::t_uint const nopen;
    HJCFIT::t_uint const nshut;
  }
  %{ namespace {
    HJCFIT::t_real HJCFIT_MissedEventsG_tmax_get(HJCFIT::MissedEventsG *_in) {
      return _in->get_tmax();
    }
    HJCFIT::t_uint HJCFIT_MissedEventsG_nmax_get(HJCFIT::MissedEventsG *_in) {
      return _in->get_nmax();
    }
    void HJCFIT_MissedEventsG_nmax_set(HJCFIT::MissedEventsG *_in, HJCFIT::t_uint _n) {
      _in->set_nmax(_n);
    }
    HJCFIT::t_uint HJCFIT_MissedEventsG_nopen_get(HJCFIT::MissedEventsG *_in) {
      return _in->get_qmatrix().nopen;
    }
    HJCFIT::t_uint HJCFIT_MissedEventsG_nshut_get(HJCFIT::MissedEventsG *_in) {
      return _in->get_qmatrix().nshut();
    }
  }%}
};

%clear MissedEventsG;

