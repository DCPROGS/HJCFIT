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
%typemap("default") DCProgs::t_real maxs { $1 = 0e0; };
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


// Rather than creating a wrapper, we figure roots directly in python. 
// This is easier and more sensible. 
%pythoncode %{
  def find_roots(determinant, intervals=None, tolerance=1e-8):
     """ Computes roots for each interval. 
     
         :param determinant: 
           A function or functor of a single variable.
         :param intervals:
           A list of items `[(a0, b0), ..., (a1, b1)]`, where `(a, b)` is the interval over which to
           look for roots. 

           If this object is None (default), then uses :py:meth:`find_root_intervals` to figure out
           the  intervals.
         :param tolerance:
           Tolerance criteria. Used to determine multiplicity.
         :returns: A list of items `(root, multiplicity)`.
     """
     from scipy.optimize import brentq
     from numpy import abs, sum
     from numpy.linalg import eig
     if intervals is None:
       intervals = [u[0] for u in find_root_intervals(determinant)]
     result = []
     for interval in intervals:
       # left, right: limit of interval.
       left, right = determinant(interval)
       if left * right < 0: root = brentq(determinant, *interval)
       elif left < 0:
         root, value, ierr, numfunc = fminboud(lambda x: -determinant(x), left, right)
         if abs(value) > tolerance: continue
       else:
         root, value, ierr, numfunc = fminboud(determinant, left, right)
         if abs(value) > tolerance: continue

       eigenvalues = eig(determinant.H(root))[0]
       multiplicity = sum(abs(eigenvalues - root) < tolerance)
       if left * right < 0:
         if multiplicity == 0 or multiplicity % 2 != 1: multiplicity = 1
       else:
         if multiplicity == 0 or multiplicity % 2 != 0: multiplicity = 2
       result.append((root, multiplicity));
     return result;

  def plot_roots(determinant, intervals=None, figure=None, main=None, lines=None, size=1000,
                 tolerance=1e-8): 
     """ Computes and plots roots. 

         :param determinant: 
           A function or functor of a single variable.
         :param intervals:
           A list of items `[(a0, b0), ..., (a1, b1)]`, where `(a, b)` is the interval over which to
           look for roots. 

           If this object is None (default), then uses :py:meth:`find_root_intervals` to figure out
           the  intervals.
         :param main:
           A dictionary of values with which to plot the determinant.
         :param lines:
           A dictionary of values with which to plot the roots.
           
         :returns: A figure
     """
     from matplotlib import pyplot as plt
     from numpy import arange, min, max, array

     if intervals is None:
       intervals = [u[0] for u in find_root_intervals(determinant)]
     intervals = array(intervals)
     if main is None: main = {}
     if lines is None: lines = {}

     mini = min(intervals)
     maxi = max(intervals)
     
     roots = find_roots(determinant, intervals, tolerance)
     x = arange(mini, maxi+(maxi-mini)/float(size)*0.5, (maxi-mini)/float(size))
     y = determinant(x)

     figure = plt.figure()
     ax = figure.add_subplot(111)
     ax.plot(x, y, **main)
     ax.set_xlim((mini, maxi))

     ymin, ymax = min(y), max(y)
     ymin = ymin - (ymax - ymin) * 0.05
     ymax = ymax + (ymax - ymin) * 0.05
     ax.set_ylim((ymin, ymax))

     for root in roots:
       ax.plot([root[0], root[0]], [ymin, 0], **lines)
     ax.plot([x[0], x[-1]], [0, 0])

     return figure
%}
