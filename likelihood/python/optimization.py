""" Subpackage for likelihood optimization. """

class Likelihood(object):
  """ A function that can be optimized via scipy or other. """
  def __init__(self, nopen, graph_matrix, intervals, tau, tcrit=None):
    """ Creates the likelihood functor. 
    
        This functor does not incorporate constraints per se. There are quite a few ways that such
        constraints could be defined. It is not clear which should be used a-priori. 


        :param nopen: 
          Number of open states. 
        :param graph_matrix:
          Defines the reaction graph. This is a list of list (matrix) where each element is either
          "V" (variable), 0 (no direct reaction between these two states), or a number (Fixed
          reaction rate). The open-states should be in the top-left corner.
          
          .. code-bloc:: python
          
             [ ["V", 0, 0.1], 
               [0, "V", "V"],
               ["V", "V", "V"] ]

           In the exampl above, there are no direct reactions transforming state (1) into state
           (2), and the reaction rate from state (1) to state (3) is fixed. All others will be
           optimized, subject to the intrinsic constraints (sum over each row is zero, diagonal
           elements are negative) and extrinsic contraints (defined below).

        :param intervals: 
          List of shut and open time intervals. If a list of list is provided, then each defines a
          single patch clamp acquisition. The first interval of each list should correspond to an
          open-time event.

        :param Float tau:
          Max length of missed events.

        :param tcrit: 
          If None, uses equilbrium occupancies. Otherwise, computes and uses CHS vectors.
    """
    from copy import deepcopy
    from numpy import array, zeros, count_nonzero
    from .likelihood import QMatrix



    nstates = len(graph_matrix)
    if nstates < 2: raise ValueError('graph_matrix should contain two or more states.')

    is_fixed = lambda u: (not isinstance(u, str)) and abs(u) > 1e-12
    self.fixed_mask = array([ [is_fixed(graph_matrix(i, j)) for j in xrange(nstates)] 
                              for i in xrange(nstates) ])
    """ Components are True if rates are fixed and non-zero.
    
        If this is None, then there are no fixed components.
    """
    if count_nonzero(self.fixed_mask) == 0: self.fixed_mask = None
    

    if len(intervals) == 0:
      raise ValueError('intervals should a non-empty list of list of intervals.')
    if not hasattr(intervals[1], '__len__'): intervals = [intervals]
    self.intervals = [array(u, dtype='float64') for u in intervals]
    """ List of list of open/shut intervals. """

    self.tau = tau
    """ Max length of missed events. """

    self.loglikelihood = True
    """ Whether to compute the log-likelihood. """

    if not hasattr(tcrit, '__len__'): tcrit = [tcrit] * len(self.intervals)
    if len(tcrit) != len(self.intervals):
      raise ValueError('Number of tcrit and intervals is different.')
    self.tcrit = deepcopy(tcrit)
    """ Critical time t.

        list of critical times t. These times can be None, in which case the equilibrium occupancies
        will be used.
    """





    # sanity checks.
    if any(len(u) != nstates for u in graph_matrix):
      raise ValueError('graph_matrix should be square')
    if self.fixed_mask is not None:
      for i in xrange(nstates):
        if self.fixed_mask[i, i]:
          raise ValueError('Rate from {0} to {0} cannot be fixed.'.format(i))
    for i in xrange(nstates):
      for j in xrange(i+1, nstates):
        if self.disconnect_mask[i, j]: 
          raise ValueError('No reversability for {0} to {1} transformations.'.format(i, j))

 

    # Private members:
    # A qmatrix object with fixed values already set is created once and for all.
    qmatrix = QMatrix(zeros(nstates, nstates), nopen)
    if self.fixed_mask is not None:
      get_fixed = lambda u: (u if is_fixed(u) else 0)
      fixed = array([ [get_fixed(graph_matrix(i, j)) for j in xrange(nstates)] 
                       for i in xrange(nstates) ])
      fixed = fixed[self.fixed_mask].flatten()
      qmatrix.matrix[self.fixed_mask] = fixed

    self._current_qmatrix = qmatrix
    """ Current QMatrix object. """

  @property
  def nopen(self):
    """ Number of open states. """ 
    return self._current_qmatrix.nopen
  @property
  def nshut(self):
    """ Number of open states. """ 
    return self._current_qmatrix.nshut
  @property
  def nstates(self):
    """ Number of states. """ 
    return self._current_qmatrix.matrix.shape[0]


  def vector(self, x):
    """ Computes likelihood for each interval. """
    from numpy import array, zeros, count_nonzero
    from .likelihood import create_missed_eventsG, chained_likelihood

    x = array(x)

    # Figures out whether fixed components have been removed from x or not
    remove_fixed = False
    if self.fixed_mask is not None and x.size != self.qmatrix.size:
      if x.size != self.qmatrix.size - count_nonzero(self.fixed_mask):
        raise ValueError('x input has incorrect size.')
      remove_fixed = True

    # set x values to matrix
    qmatrix = self._current_qmatrix
    # two possibilities here, depending on whether x includes fixed components or not.
    if self.fixed_mask is None: qmatrix.matrix.flat = array(x).flat
    elif remove_fixed: qmatrix.matrix[not self.fixed_mask].flat = array(x).flat
    else:
      qmatrix.matrix[not self.fixed_mask] = x.reshape(qmatrix.matrix.shape)[not self.fixed_mask]

    # create missed events G function
    missed_eventsG = create_missed_eventsG(qmatrix, self.tau)

    # computes likelihood for each state
    initial_occ, final_occ = None, None
    results = zeros(len(self.intervals), dtype='float64')
    for i, (intervals, tcrit) in enumerate(zip(self.intervals, self.tcrit)):

      # figure out initial and final states
      if tcrit is None: 
        if initial_occ is not None: 
          initial, final = initial_occ, final_occ
        else:
          initial = missed_eventsG.initial_occupancies
          final = missed_eventsG.final_occupancies
      else:
        initial = missed_eventsG.initial_CHS_occupancies(tcrit)
        final = missed_eventsG.final_CHS_occupancies(tcrit)

      results[i] = chained_likelihood(missed_eventsG, intervals, initial, final)
    
    # return result
    return results

  def __call__(self, x): 
    """ Computes likelihood for q given input vector x. """
    from numpy import sum, prod, log

    results = self.vector(x)
    return sum(log(results)) if self.loglikelihood else prod(results)

  def intrinsic_linear_equalities(self, with_fixed_points=False):
    """ A matrix defining intrinsic linear equality constraints. 
          
        The intrinsic linear constraints enforce that rows sum to zero.
        The fixed components are not removed from the list of variables (otherwise rows would not
        sum to zero)

        :param with_fixed_points:
          Whether to include rows in the return matrix for each null component. This makes it easier
          to define fixed components as linear constraints.
    """
    from numpy import zeros, count_nonzero, nonzero
    if self.fixed_mask is None: with_fixed_points = False
    nrows = self.nstates * self.nstates
    ncols = self.nstates 
    if with_fixed_points: ncols += count_nonzero(self.fixed_mask)
    result = zeros((nrows, ncols), dtype='float64')

    for i in xrange(self.nstates):
      result[i, i*self.nstates: (i+1)*self.nstates] = 1

    if with_fixed_points:  
      for j, value in enumerate(nonzero(self.fixed_mask)): result[i+j, value] = 1

    return result


def kernel_projection(likelihood, x0, A=None):
  """ Returns a functor that operates in the kernel of linear constraints of the likelihood. """
  from numpy import array, dot, concatenate, zeros, argsort
  from numpy.linalg import svd, norm
  # get linear constraints
  if isinstance(A, Likelihood): A = likelihood.intrinsic_linear_constraints(True)
  x0 = array(x0).copy()
  
  # Add missing vectors to get nullspace. 
  N = x0.size - A.shape[0] 
  A = concatenate((A, zeros((N, x0.size))))

  # Now get kernel
  left, sings, right = svd(A)
  # it is defined as the smallest N singular values.
  indices = argsort( abs(sings) )
  kernel = right[indices[:N], :] 
  # we now normalize the kernel
  for r in kernel: r[:] = r[:] / norm(r)

  # and create to affine transforms to go back and forth between one set of coordinate and the other
  to_kernel_coords = lambda x: dot(kernel, x - x0)
  from_kernel_coords = lambda x: dot(x, kernel) + x0

  def kernel_likelihood(x):
    whole_x = from_kernel_coords(x) 
    return likelihood(whole_x)

  kernel_likelihood.to_kernel_coords = to_kernel_coords
  kernel_likelihood.to_kernel_coords.__doc__  \
      = """ Mapping from whole coordinate space to kernel. """
  kernel_likelihood.from_kernel_coords = from_kernel_coords
  kernel_likelihood.from_kernel_coords.__doc__ \
      = """ Mapping from kernel back to whole coordinate space. """

  return kernel_likelihood
