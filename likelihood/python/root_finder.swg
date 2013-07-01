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
         ":raises: ArithmeticError when complex eigenvalues or NaN are encountered.\n\n"
         ":returns: A list `[(a, b), multiplicity]`, where `(a, b)` denotes an interval, "
                   "and `multiplicity` the multiplicity of the root. All roots with "
                   "`multiplicity > 1` will have a size of `tolerance` or smaller.")
         py_find_root_intervals(DeterminantEq const &det);
%{
   namespace DCProgs { 
     // Wrapper transforms input to a more sane python format.
     // This way, we avoid having special wrappers for RootInterval and Root. It doesn't make sense
     // in a dynamically typed environment.
     PyObject* py_find_root_intervals(DeterminantEq const &det, t_real mins=1e8, t_real maxs=0e0,
                                      t_real tolerance=1e-8) {
       try {
         std::vector<RootInterval> intervals = find_root_intervals(det, mins, maxs, tolerance);
         // Now converts it to a list of (2+1)-tuples.
         Object<PyObject> result = steal_ref(PyList_New(intervals.size()));
         if(not result) return NULL;
  
         for(Py_ssize_t i(0); i < intervals.size(); ++i) {
  
           Object<PyObject> left_limit = steal_ref(PyFloat_FromDouble(intervals[i].start));
           if(not left_limit) return NULL;
           Object<PyObject> right_limit = steal_ref(PyFloat_FromDouble(intervals[i].end));
           if(not right_limit) return NULL;
  
           Object<PyObject> interval = steal_ref(PyTuple_Pack(2, ~left_limit, ~right_limit));
           if(not interval) return NULL;
  
           Object<PyObject> multiplicity = steal_ref(PyLong_FromLong(intervals[i].multiplicity));
           if(not multiplicity) return NULL;
           
           Object<PyObject> item = steal_ref(PyTuple_Pack(2, ~interval, ~multiplicity));
           if(not item) return NULL;
  
           PyList_SET_ITEM(~result, i, item.release());
         }
         return result.release();
       } DCPROGS_CATCH(return NULL);
     }
   }
%}


%rename("find_root_intervals") py_find_root_intervals;
PyObject* py_find_root_intervals(DeterminantEq const &det, DCProgs::t_real mins, 
                                 DCProgs::t_real maxs,
                                 DCProgs::t_real tolerance=1e8);
// SWIG fails when trying to use default parameters for maxs and mins... Not sure why.
PyObject* py_find_root_intervals(DeterminantEq const &det, DCProgs::t_real mins);
PyObject* py_find_root_intervals(DeterminantEq const &det);


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
     from numpy.linalg import eig
     if intervals is None:
       intervals = [u[0] for u in find_root_intervals(determinant)]
     print intervals
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
       multiplicity = len(abs(eigenvalues - root) < tolerance)
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
