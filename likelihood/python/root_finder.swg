
%feature("autodoc", "Figures out lower bound for root.\n\n"
         "Lower bound is obtained by iteratively checking the "
         "lowest eigenvalue :math:`\\epsilon_i^s` of :math:`H(s_i)`, "
         "where :math:`s_i` is the guess at iteration :math:`i`. If the "
         "lower eigenvalue is lower than :math:`s_i`, than "
         ":math:`s_{i+1} = \\epsilon_i^s + \\alpha(\\epsilon_i^s - s_i)` "
         "is created.\n\n"
         ".. Parameters:\n"
         "   det: DeterminantEq\n"
         "     Function for which to guess bound for lower root.\n"
         "   start: Number\n"
         "     Starting point. Defaults to zero.\n"
         "   alpha: Number\n"
         "     Factor from which to determine next best guess.\n"
         "   itermax: Integer\n"
         "     Maximum number of iterations before bailing out.\n") find_lower_bound_for_roots;

DCProgs::t_real find_lower_bound_for_roots(DeterminantEq const &det,
                                           DCProgs::t_real start,
                                           DCProgs::t_real alpha = 2e0, 
                                           DCProgs::t_int itermax = 100);

%{
   namespace DCProgs { 
     PyObject* py_find_root_intervals(DeterminantEq const &_det, t_real _mins, t_real _maxs,
                                      t_real _tolerance) {
       try {
         std::vector<RootInterval> intervals = find_root_intervals(_det, _mins, _maxs, _tolerance);
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
PyObject* py_find_root_intervals(DeterminantEq const &_det, DCProgs::t_real _mins, DCProgs::t_real _maxs,
                                 DCProgs::t_real _tolerance);

