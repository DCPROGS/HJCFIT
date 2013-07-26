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

%feature("autodoc", "Figures out lower bound for root.\n\n"
         "Lower bound is obtained by iteratively checking the "
         "lowest eigenvalue :math:`\\epsilon_i^s` of :math:`H(s_i)`, "
         "where :math:`s_i` is the guess at iteration :math:`i`. If the "
         "lower eigenvalue is lower than :math:`s_i`, than "
         ":math:`s_{i+1} = \\epsilon_i^s + \\alpha(\\epsilon_i^s - s_i)` "
         "is created.\n\n"
         ":param DeterminantEq det: \n"
         "     Function for which to guess bound for lower root.\n"
         ":param Number start:\n"
         "     Starting point. Defaults to zero.\n"
         ":param Number alpha:\n"
         "     Factor from which to determine next best guess.\n"
         ":param Integer itermax:\n"
         "     Maximum number of iterations before bailing out.\n")
         find_lower_bound_for_roots(DeterminantEq const &det);

DCProgs::t_real find_lower_bound_for_roots(DeterminantEq const &det,
                                           DCProgs::t_real start = 1e8,
                                           DCProgs::t_real alpha = 2e0, 
                                           DCProgs::t_int itermax = 100);

%feature("autodoc", "Figures out upper bound for root.\n\n"
         "Lower bound is obtained by iteratively checking the "
         "largest eigenvalue :math:`\\epsilon_i^s` of :math:`H(s_i)`, "
         "where :math:`s_i` is the guess at iteration :math:`i`. If the "
         "lower eigenvalue is lower than :math:`s_i`, than "
         ":math:`s_{i+1} = \\epsilon_i^s + \\alpha(\\epsilon_i^s - s_i)` "
         "is created.\n\n"
         ":param DeterminantEq det: \n"
         "     Function for which to guess bound for lower root.\n"
         ":param Number start:\n"
         "     Starting point. Defaults to zero.\n"
         ":param Number alpha:\n"
         "     Factor from which to determine next best guess.\n"
         ":param Integer itermax:\n"
         "     Maximum number of iterations before bailing out.\n")
         find_upper_bound_for_roots(DeterminantEq const &det);

DCProgs::t_real find_upper_bound_for_roots(DeterminantEq const &det,
                                           DCProgs::t_real start = 1e8,
                                           DCProgs::t_real alpha = 2e0, 
                                           DCProgs::t_int itermax = 100);


// Creates a typemap to convert std::vector<RootInterval> to a python list.
// Each item is a "((a, b), multiplicity)".
%typemap(out) std::vector<DCProgs::RootInterval> { 
  Object<PyObject> py_result = steal_ref(PyList_New($1.size()));
  if(not py_result) return NULL;
  
  Py_ssize_t i(0);
  for(DCProgs::RootInterval const &c_interval: $1) {
  
    Object<> left_limit = steal_ref(PyFloat_FromDouble(c_interval.start));
    if(not left_limit) return NULL;
    Object<> right_limit = steal_ref(PyFloat_FromDouble(c_interval.end));
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

%feature("autodoc", "Returns intervals for searching roots.\n\n"
         "Lower bound is obtained by iteratively checking the "
         "lowest eigenvalue :math:`\\epsilon_i^s` of :math:`H(s_i)`, "
         "where :math:`s_i` is the guess at iteration :math:`i`. If the "
         "lower eigenvalue is lower than :math:`s_i`, than "
         ":math:`s_{i+1} = \\epsilon_i^s + \\alpha(\\epsilon_i^s - s_i)` "
         "is created.\n\n"
         ":param DeterminantEq det:\n"
         "     Function for which to guess bound for lower root.\n"
         ":param mins: Number\n"
         "     Starting point. Defaults to a large number, larger than `maxs`."
             " In that case, the lower bound is determined using"
             " :py:meth:`find_lower_bound_for_roots`. "
             "If it is smaller than `maxs`, than all roots should be larger than this number."
         ":param Number maxs:\n"
         "     Roots should all be smaller than this number.\n"
         ":param Number tolerance:\n"
         "     Size of the smallest possible intervals. Intervals smaller than this are likely to"
             " have multiple roots.\n\n"
         ":raises: ArithmeticError when NaN is encountered or when eigenvalue solver fails.\n\n"
         ":returns: A list `[(a, b), multiplicity]`, where `(a, b)` denotes an interval, "
                   "and `multiplicity` the multiplicity of the root. All roots with "
                   "`multiplicity > 1` will have a size of `tolerance` or smaller.")
         find_root_intervals;
// Gives the arguments a default.
// Despite what it advertises, SWIG doesn't seem to deal well with multiple default arguments in
// c++.
%typemap("default") DCProgs::t_real mins { $1 = 1e8; };
%typemap("default") DCProgs::t_real maxs { $1 = 1e1; };
%typemap("default") DCProgs::t_real tolerance { $1 = 1e-8; };
std::vector<DCProgs::RootInterval>
   find_root_intervals(DeterminantEq const &det, DCProgs::t_real mins, 
                       DCProgs::t_real maxs, DCProgs::t_real tolerance);


%feature("autodoc", "Find intervals for roots using brute force.\n\n"
         "Computes all values between mins and maxs, for a given resolution. "
         "If determinant changes sign between two values, or if it comes to within tolerance of "
         "zero, then computes eigenvalues of H to determine possible multiplicity.\n\n" 
         ":param DeterminantEq det:\n"
         "     Function for which to guess bound for lower root.\n"
         ":param Number resolution:\n"
         "     Resolution at which to sample interval.\n"
         ":param Number mins: \n"
         "     Starting point. Defaults to a large number, larger than `maxs`."
             " In that case, the lower bound is determined using"
             " :py:meth:`find_lower_bound_for_roots`. "
             "If it is smaller than `maxs`, than all roots should be larger than this number."
         ":param Number maxs:\n"
         "     Roots should all be smaller than this number.\n"
         ":param Number tolerance:\n"
         "     Tolerance below which the value of the determinant is considered "
               "_close to zero_.\n\n"
         ":raises: ArithmeticError when NaN is encountered or when eigenvalue solver fails.\n\n"
         ":returns: A list `[(a, b), multiplicity]`, where `(a, b)` denotes an interval, "
                   "and `multiplicity` the multiplicity of the root. All roots with "
                   "`multiplicity > 1` will have a size of `tolerance` or smaller.")
         find_root_intervals_brute_force; 
%typemap("default") DCProgs::t_real resolution { $1 = 1e1; };
std::vector<DCProgs::RootInterval> 
  find_root_intervals_brute_force(DCProgs::DeterminantEq const &det, DCProgs::t_real resolution,
                                  DCProgs::t_real mins, DCProgs::t_real maxs, 
                                  DCProgs::t_real tolerance);
%clear DCProgs::t_real resolution;
%clear DCProgs::t_real mins;
%clear DCProgs::t_real maxs;
%clear DCProgs::t_real tolerance;


