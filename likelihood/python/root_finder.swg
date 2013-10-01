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

%rename("_find_lower_bound_for_roots") find_lower_bound_for_roots;
DCProgs::t_real find_lower_bound_for_roots( DCProgs::DeterminantEq const &_det,
                                            DCProgs::t_real _start,
                                            DCProgs::t_real _alpha,
                                            DCProgs::t_uint _itermax );
%pythoncode %{
def find_lower_bound_for_roots(determinant, start=0e0, alpha=2e0, itermax=100):
  """ Figures out lower bound for root.

       The lower bound is obtained by iteratively checking the lowest eigenvalue
       :math:`\\epsilon_i^s` of :math:`H(s_i)`, where :math:`s_i` is the guess at iteration
       :math:`i`. If the lowest eigenvalue is smaller than :math:`s_i`, then
       :math:`s_{i+1} = \\epsilon_i^s + \\alpha(\\epsilon_i^s - s_i)` is created.

       :param DeterminantEq det:
            Function for which to guess bound for lower root.
       :param Number start:
            Starting point. 
       :param Number alpha: 
            Factor from which to determine next best guess. 
       :param Integer itermax:
            Maximum number of iterations before bailing out.
  """
  return _likelihood._find_lower_bound_for_roots(determinant, start, alpha, itermax)
%}

%rename("_find_upper_bound_for_roots") find_upper_bound_for_roots;
DCProgs::t_real find_upper_bound_for_roots( DCProgs::DeterminantEq const &_det,
                                            DCProgs::t_real _start,
                                            DCProgs::t_real _alpha,
                                            DCProgs::t_uint _itermax );
%pythoncode %{
def find_upper_bound_for_roots(determinant, start=0e0, alpha=2e0, itermax=100):
  """ Figures out upper bound for root.

      The upper bound is obtained by iteratively checking the largest eigenvalue
      :math:`\\epsilon_i^s` of :math:`H(s_i)`, where :math:`s_i` is the guess at iteration
      :math:`i`. If the largest eigenvalue is larger than :math:`s_i`, then
      :math:`s_{i+1} = \\epsilon_i^s + \\alpha(\\epsilon_i^s - s_i)` is created.


      :param DeterminantEq det: 
           Function for which to guess bound for lower root.
      :param Number start:
           Starting point. 
      :param Number alpha:
           Factor from which to determine next best guess.
      :param Integer itermax:
           Maximum number of iterations before bailing out.
  """
  return _likelihood._find_upper_bound_for_roots(determinant, start, alpha, itermax)
%}


// Creates a typemap to convert std::vector<RootInterval> to a python list.
// Each item is a "((a, b), multiplicity)".
%typemap(out) std::vector<DCProgs::RootInterval> { 
  Object<PyObject> py_result = steal_ref(PyList_New($1.size()));
  if(not py_result) return NULL;
  
  Py_ssize_t i(0);
  for(DCProgs::RootInterval const &c_interval: $1) {
  
    Object<> left_limit = steal_ref(convert_real_to_python(c_interval.start));
    if(not left_limit) return NULL;
    Object<> right_limit = steal_ref(convert_real_to_python(c_interval.end));
    if(not right_limit) return NULL;
  
    Object<> py_interval = steal_ref(PyTuple_Pack(2, ~left_limit, ~right_limit));
    if(not py_interval) return NULL;
  
    Object<> multiplicity = steal_ref(PyLong_FromLong(c_interval.multiplicity));
    if(not multiplicity) return NULL;
    
    Object<> item = steal_ref(PyTuple_Pack(2, ~py_interval, ~multiplicity));
    if(not item) return NULL;
  
    PyList_SET_ITEM(~py_result, (i++), item.release());
  }
  $result = py_result.release();
}

%rename("_find_root_intervals") find_root_intervals;
std::vector<DCProgs::RootInterval>
   find_root_intervals( DCProgs::DeterminantEq const &det, DCProgs::t_real mins, 
                        DCProgs::t_real maxs, DCProgs::t_real tolerance );
%pythoncode %{
  def find_root_intervals(determinant, mins=None, maxs=None, tolerance=1e-8):
    """ Returns intervals for searching roots.

        Finds a set of intervals between mins and maxs, each containing a single root.

        :param determinant:
             :class:`DeterminantEq` instance for which to guess bound for lower root.
        :param mins: Number
             Lower bound of the interval within which to search for roots. If None, the lower bound
             is determined using :py:func:`find_lower_bound_for_roots`. 
        :param Number maxs:
             Upper bound of the interval within which to search for roots. If None, the upper bound
             is determined using :py:func:`find_upper_bound_for_roots`. 
        :param Number tolerance:
             Size of the smallest possible intervals. Intervals smaller than this are likely to
             have multiple roots.
        :raises: ArithmeticError when NaN is encountered or when eigenvalue solver fails.
        :returns:
             A list `[(a, b), multiplicity]`, where `(a, b)` denotes an interval, and
             `multiplicity` the multiplicity of the root. All roots with `multiplicity`:math:`> 1`
             will have a size of `tolerance` or smaller.
    """
    from numpy import NaN
    if mins is None: mins = NaN
    if maxs is None: maxs = NaN
    return _likelihood._find_root_intervals(determinant, mins, maxs, tolerance)
%}

%rename("_find_root_intervals_brute_force") find_root_intervals_brute_force;
std::vector<DCProgs::RootInterval> 
  find_root_intervals_brute_force( DCProgs::DeterminantEq const &det, DCProgs::t_real resolution,
                                   DCProgs::t_real mins, DCProgs::t_real maxs, 
                                   DCProgs::t_real tolerance );
%pythoncode %{
  def find_root_intervals_brute_force( determinant, resolution=1e-1, mins=None, maxs=None,
                                       tolerance=1e-8 ):
    """ Find intervals for roots using brute force.

        Computes all values between mins and maxs, for a given resolution. If the determinant
        changes sign between two values, or if it comes to within tolerance of zero, then the
        eigenvalues of H are used to determine possible multiplicity.

        :param determinant:
             Instance of :py:class:`DeterminantEq` for which to find root intervals
        :param Number resolution:
             Resolution at which to sample interval.
        :param Number mins: 
             Lower bound of interval. If None, then :py:func:`find_lower_bound_for_roots` is called.
        :param Number maxs:
             Upper bound of interval. If None, then :py:func:`find_upper_bound_for_roots` is called.
        :param Number tolerance:
             Tolerance below which the value of the determinant is considered "_close to zero_".

        :raises: ArithmeticError when NaN is encountered or when eigenvalue solver fails.\n\n"

        :returns:
           A list `[(a, b), multiplicity]`, where `(a, b)` denotes an interval, and `multiplicity`
           the multiplicity of the root. All roots with `multiplicity > 1` will have a size of
           `tolerance` or smaller.
    """
    from numpy import NaN
    if mins is None: mins = NaN
    if maxs is None: maxs = NaN
    return _likelihood._find_root_intervals_brute_force(determinant, resolution, mins, maxs, tolerance)
%}
%{
  // Helper function to convert output of root finding objects to cpp
  std::vector<DCProgs::Root> convert_to_roots(PyObject *_in) {
    Object<> iterator = steal_ref(PyObject_GetIter(_in)); 
    if(not iterator) throw DCProgs::errors::PythonErrorAlreadyThrown();

    std::vector<DCProgs::Root> result;
    while(Object<> item = steal_ref(PyIter_Next(~iterator))) {

      double root;
      int multiplicity;

      if(not PyArg_ParseTuple(~item, "di", &root, &multiplicity))
        throw DCProgs::errors::PythonErrorAlreadyThrown();
      result.emplace_back(DCProgs::t_real(root), DCProgs::t_int(multiplicity));
    }
    if(PyErr_Occurred() != NULL) throw DCProgs::errors::PythonErrorAlreadyThrown();
    
    return result;
  }
%}
