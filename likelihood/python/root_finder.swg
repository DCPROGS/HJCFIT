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
  def find_roots(determinant, intervals=None):
     """ Computes roots for each interval. 
     
         :param determinant: 
           A function or functor of a single variable.
         :param intervals:
           A list of items `[(a, b), multiplicity]`, where `(a, b)` is the interval over which to
           look for roots. If multiplicity is larger than 1, then the root is set to `(a+b) * 0.5`.
           If it is one, then the root is found using scipy.optimize.brentq.
           If it is zero, then inteval is ignored. 

           If this object is None (default), then uses :py:meth:`find_root_intervals` to figure out
           the  intervals.
         :returns: A list of items `(root, multiplicity)`.
     """
     from scipy.optimize import brentq
     if intervals is None: intervals = find_root_intervals(determinant) 
     result = []
     for interval in intervals:
       if interval[1] == 1: 
         root = brentq(determinant, *interval[0])
         result.append((root, interval[1]));
       elif interval[1] > 1: 
         result.append(((interval[0][0]+interval[0][1]) * 0.5, interval[1]));
     return result;
  def plot_roots(determinant, intervals=None, figure=None, main=None, lines=None): 
     """ Computes and plots roots. 

         :param determinant: 
           A function or functor of a single variable.
         :param intervals:
           A list of items `[(a, b), multiplicity]`, where `(a, b)` is the interval over which to
           look for roots. If multiplicity is larger than 1, then the root is set to `(a+b) * 0.5`.
           If it is one, then the root is found using scipy.optimize.brentq.
           If it is zero, then inteval is ignored. 

           If this object is None (default), then uses :py:meth:`find_root_intervals` to figure out
           the  intervals.
         :param main:
           A dictionary of values with which to plot the determinant.
         :param lines:
           A dictionary of values with which to plot the roots.
           
         :returns: A figure
     """
     from matplotlib.pyplot import Figure
     from numpy import arange, min, max

     if intervals is None: intervals = find_root_intervals(determinant) 
     if main is None: main = {}
     if lines is None: lines = {}

     mini = min(intervals[0][0])
     maxi = min(intervals[0][0])
     for interval in intervals[1:]:
       mini = min([mini] + list(interval[0]))
       maxi = max([maxi] + list(interval[0]))
     
     roots = find_roots(determinant, intervals)
     x = arange(mini, maxi)
     y = determinant(x)

     figure = Figure()
     ax = figure.add_subplot(111)
     ax.plot(x, y, **main)
     ax.set_xlim((mini, maxi))
     ymin, ymax = min(y), max(y)
     ymin = ymin - (ymax - ymin) * 0.05
     ymax = ymax + (ymax - ymin) * 0.05
     ax.set_ylim((ymin, ymax))

     for root in roots:
       ax.plot([root[0], root[0]], [ymin, 0], **lines)

     return figure

%}