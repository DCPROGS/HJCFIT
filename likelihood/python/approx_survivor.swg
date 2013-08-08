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

%feature("autodoc", "Computes approximate missed-event survivor function.\n\n") ApproxSurvivor;
%feature("autodoc", "Open to close transitions") ApproxSurvivor::af;
%feature("autodoc", "Closed to open transitions") ApproxSurvivor::fa;
%feature("shadow") ApproxSurvivor::ApproxSurvivor %{
  def __init__(self, *args, **kwargs):
    """ Initializes the approximate survivor function.

        There are three possible ways to instantiate this object:

        >>> ApproxSurvivor(determinant_af, roots_af, determinant_fa, roots_fa)
        >>> ApproxSurvivor(matrix, nopen, tau, **kwargs)
        >>> ApproxSurvivor(qmatrix, tau, **kwargs)

        The parameters bewteen brackets are optional. The last two versions will try and calculate
        the roots of the determinantal equations automatically. A number of parameters can be given
        to control this process. The keyword arguments are listed below.

        :param determinant_af:
            A :class:`DeterminantEq` instance, specifically for the af block.
        :param roots_af:
           The roots of the af determinantal equation. The should come in the format
           `[(root, multiplicity), (root, multiplicity), ...]`.
        :param determinant_fa:
           A :class:`DeterminantEq` instance, specifically for the fa block.
           It should the transpose of `determinant_af`. It is required so 
           it need not be recomputed, since it most likely already exists.
        :param roots_fa:
           The roots of the fa determinantal equation. The should come in the format 
           `[(root, multiplicity), (root, multiplicity), ...]`.
        :param matrix:
           An object convertible to a square matrix.
        :param integer nopen:
           Number of open states.
        :param qmatrix:
           A :class:`QMatrix` instance.
        :param float xtol: 
           Tolerance criteria when computing roots using py:meth:`dcprogs.likelihood.brentq`.
           Defaults to 1e-12.
        :param float rtol:
           Tolerance criteria when computing roots using py:meth:`dcprogs.likelihood.brentq`.
           Defaults to 1e-12.
        :param float itermax: \n"
           Maximum number of iterations when computing roots using
           py:meth:`dcprogs.likelihood.brentq`.
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
    det_af   = args[0] if len(args) > 0 else kargs.get("det_af", None)
    if isinstance(det_af, DeterminantEq):
      roots_af = args[1] if len(args) > 1 else kargs.get("roots_af", None)
      if roots_af is None: raise TypeError('Missing roots_af argument')
      det_fa   = args[2] if len(args) > 2 else kargs.get("det_fa", None)
      if det_fa is None: raise TypeError('Missing det_fa argument')
      roots_fa = args[3] if len(args) > 3 else kargs.get("roots_fa", None)
      if roots_fa is None: raise TypeError('Missing roots_fa argument')
      this = $action(det_af, roots_af, det_fa, roots_fa)
      try: self.this.append(this)
      except: self.this = this
      return

    # Try second and third instanciation
    matrix = args[0] if len(args) > 0 else kwargs.get('matrix', kwargs.get('qmatrix', None))
    if matrix is None: 
      raise TypeError("Expected the q-matrix or a matrix + number of open states on input.")
    index = 1
    if not isinstance(matrix, QMatrix):
      matrix = array(matrix, internal_dtype)
      nopen = int(args[index]) if len(args) > index else kwargs.get('nopen', None)
      if nopen is None: 
        raise TypeError("A matrix was given on input, but without a nopen parameter.")
      matrix = QMatrix(matrix, nopen)
      index += 1

    tau = args[index] if len(args) > index else kwargs.get('tau', None)
    if tau is None: raise TypeError("No tau given on input.")
    xtol        = args[index+1] if len(args) > index+1 else kwargs.get('xtol', 1e-12)
    rtol        = args[index+3] if len(args) > index+2 else kwargs.get('rtol', 1e-12)
    itermax     = int(args[index+3] if len(args) > index+3 else kwargs.get('itermax', 100))
    lower_bound = args[index+4] if len(args) > index+4 else kwargs.get('lower_bound', None)
    upper_bound = args[index+5] if len(args) > index+5 else kwargs.get('upper_bound', None)
    if lower_bound is None: lower_bound = NaN
    if upper_bound is None: upper_bound = NaN
    
    this = $action(matrix, tau, xtol, rtol, itermax, lower_bound, upper_bound)
    try: self.this.append(this)
    except: self.this = this
%}


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
                    DCProgs::t_real _xtol, DCProgs::t_real _rtol, DCProgs::t_uint _itermax,
                    DCProgs::t_real _lowerbound, DCProgs::t_real _upperbound ) {
      using namespace DCProgs;
      DeterminantEq const determinant_af(*_matrix, _tau);
      DeterminantEq const determinant_fa = determinant_af.transpose(); 
      std::vector<Root> const roots_af = find_roots( determinant_af, _xtol, _rtol, _itermax,
                                                     _lowerbound, _upperbound );
      std::vector<Root> const roots_fa = find_roots( determinant_fa, _xtol, _rtol, _itermax, 
                                                     _lowerbound, _upperbound );
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

