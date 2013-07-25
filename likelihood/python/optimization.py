""" Subpackage for likelihood optimization. """
__docformat__ = "restructuredtext en"
__all__ = ['Likelihood', 'kernel_projection', 'minimizer_constraints']

class Likelihood(object):
  """ A function that can be optimized via scipy or other. """
  def __init__(self, nopen, graph_matrix, bursts, tau, tcrit=None):
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

        :param bursts: 
          List of shut and open time bursts. If a list of list is provided, then each defines a
          single patch clamp acquisition. The first interval of each list should correspond to an
          open-time event.

        :param Float tau:
          Max length of missed events.

        :param Float tcrit: 
          If None, uses equilbrium occupancies. Otherwise, computes and uses CHS vectors.
    """
    from numpy import array, zeros, count_nonzero
    from .likelihood import QMatrix, create_bursts
    from .. import internal_dtype
  


    nstates = len(graph_matrix)
    if nstates < 2: raise ValueError('graph_matrix should contain two or more states.')

    is_fixed = lambda u: not isinstance(u, str) 
    self.fixed_mask = array([ [is_fixed(graph_matrix[i][j]) for j in xrange(nstates)] 
                              for i in xrange(nstates) ])
    """ Components are True if rates are fixed and non-zero.
    
        If this is None, then there are no fixed components.
    """
    if count_nonzero(self.fixed_mask) == 0: self.fixed_mask = None
    

    if len(bursts) == 0:
      raise ValueError('bursts should a non-empty list of list of intervals.')
    if not hasattr(bursts[0], '__len__'): bursts = [bursts]
    self.bursts = create_bursts(bursts)
    """ List of bursts. """

    self.tau = tau
    """ Max length of missed events. """

    self.loglikelihood = True
    """ Whether to compute the log-likelihood. """

    self.tcrit = tcrit
    """ Critical time t.

        If None, uses equilbrium occupancies. Otherwise, computes and uses CHS vectors.
    """





    # sanity checks.
    if any(len(u) != nstates for u in graph_matrix):
      raise ValueError('graph_matrix should be square')
    if self.fixed_mask is not None:
      for i in xrange(nstates):
        if self.fixed_mask[i, i]:
          raise ValueError('Rate from {0} to {0} cannot be fixed.'.format(i))
    
    is_disconnect = lambda i, j: (not isinstance(graph_matrix[i][j], str)) \
                                 and abs(graph_matrix[i][j]) < 1e-12
    for i in xrange(nstates):
      for j in xrange(i+1, nstates):
        if is_disconnect(i, j) != is_disconnect(j, i):
          raise ValueError( 'No reversability for transformations from states {0} to {1}.'         \
                            .format(i, j) )

 

    # Private members:
    # A qmatrix object with fixed values already set is created once and for all.
    qmatrix = QMatrix(zeros((nstates, nstates), dtype=internal_dtype), nopen)
    if self.fixed_mask is not None:
      get_fixed = lambda u: (u if is_fixed(u) else 0)
      fixed = array([ [get_fixed(graph_matrix[i][j]) for j in xrange(nstates)] 
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
    from numpy import array, count_nonzero, bitwise_not
    from .likelihood import create_missed_eventsG, compute_bursts

    x = array(x)

    # Figures out whether fixed components have been removed from x or not
    remove_fixed = False
    if self.fixed_mask is not None and x.size != self._current_qmatrix.matrix.size:
      if x.size != self._current_qmatrix.matrix.size - count_nonzero(self.fixed_mask):
        raise ValueError('x input has incorrect size.')
      remove_fixed = True

    # set x values to matrix
    qmatrix = self._current_qmatrix
    # two possibilities here, depending on whether x includes fixed components or not.
    if self.fixed_mask is None: qmatrix.matrix.flat = array(x).flat
    elif remove_fixed: qmatrix.matrix[bitwise_not(self.fixed_mask)].flat = array(x).flat
    else:
      qmatrix.matrix[bitwise_not(self.fixed_mask)] =                                               \
        x.reshape(qmatrix.matrix.shape)[bitwise_not(self.fixed_mask)]

    # create missed events G function
    missed_eventsG = create_missed_eventsG(qmatrix, self.tau)

    # figure out initial and final states
    if self.tcrit is None: 
      initial = missed_eventsG.initial_occupancies
      final = missed_eventsG.final_occupancies
    else:
      initial = missed_eventsG.initial_CHS_occupancies(self.tcrit)
      final = missed_eventsG.final_CHS_occupancies(self.tcrit)

    return compute_bursts(missed_eventsG, self.bursts, initial, final)

  def __call__(self, x): 
    """ Computes likelihood for q given input vector x. """
    from numpy import sum, prod, log, NaN, any

    results = self.vector(x)
    if self.loglikelihood:
      if any(results <= 0e0): return NaN
      return sum(log(results))
    return prod(results)

  def intrinsic_linear_equalities(self, with_fixed=False):
    """ A matrix defining intrinsic linear equality constraints. 
          
        The intrinsic linear constraints enforce that rows sum to zero.
        The fixed components are not removed from the list of variables (otherwise rows would not
        sum to zero)

        :param with_fixed:
            Whether to include rows in the return matrix for each fixed component. This makes it
            easier to define fixed components as linear constraints. 

        :returns:
           `(A, b)` with `A` a matrix and `b` a vector, such that intrinsic equality constraints are
           satisfied by a parameter vector `x` when :math:`A\cdot x - b = 0`
    """
    from numpy import zeros, count_nonzero, nonzero
    from .. import internal_dtype

    if self.fixed_mask is None: with_fixed = False
    ncols = self.nstates * self.nstates
    nrows = self.nstates 
    if with_fixed: nrows += count_nonzero(self.fixed_mask)
    A = zeros((nrows, ncols), dtype=internal_dtype)
    b = zeros(nrows, dtype=internal_dtype)

    for i in xrange(self.nstates):
      A[i, i*self.nstates: (i+1)*self.nstates] = 1

    if with_fixed:  
      for j, value in enumerate(nonzero(self.fixed_mask.flat)[0]): A[i+j+1, value] = 1
      b[i+1:] = self._current_qmatrix[self.fixed_mask].flat

    return A, b

  def intrinsic_linear_inequalities(self, with_fixed=False):
    """ A matrix defining intrinsic linear inequality constraints. 
          
        The intrinsic linear constraints enforce that non-diagonal elements should be positive. In
        conjunction with the intrinsic linear constraints, this implies that the diagonal elements
        of the Q-matrix must be negative. However, that constraint is not explicitely included here.

        :param with_fixed:
            Whether to include rows for fixed off-diagonal components.
    """
    from numpy import identity, arange, bitwise_not
    from .. import internal_dtype

    N = self.nstates*self.nstates
    keep = arange(N) % self.nstates != arange(N) // self.nstates
    if not with_fixed: keep &=  bitwise_not(self.fixed_mask.flat)
    return identity(N, dtype=internal_dtype)[keep]


  def random_starting_point(self, with_fixed=False):
    """ Returns a random starting point for the likelihood.

        :param bool with_fixed:
           If True, the return is the full 2d-matrix. Otherwise, it is a vector corresponding to the
           flattened full matrix with the fixed components removed.

        :returns: a random starting points with all linear intrinsic constraints satisfied.
    """
    from numpy import sum, bitwise_not
    from . import QMatrix
    from .random import rate_matrix as random_rate_matrix, qmatrix as random_qmatrix

    def get_qmatrix():
      result = random_rate_matrix(N=self.nstates, zeroprob=0)
      result[self.fixed_mask] = self._current_qmatrix[self.fixed_mask]
      for i in xrange(self.nstates): result[i, i] -= sum(result[i])
      return QMatrix(result, self.nopen)

    # random_qmatrix tries to ensure that the starting point is not pathological.
    # *tries* being the operative word.
    matrix = random_qmatrix(get_qmatrix=get_qmatrix).matrix
    return matrix if with_fixed else matrix[bitwise_not(self.fixed_mask)]


def kernel_projection(likelihood, x0=None, A=None):
  """ Creates functor that operates in the kernel of linear constraints of the likelihood.
  
  
      One possibility to hold linear equality constraints true is to optimize only over the 
      space perpendicular to the constraints, e.g. the nullspace or kernel of the linear equality
      matrix. This routine creates a functor that reduces the input arguments to likelihood to only
      those in the kernel.

      :param likelihood: 
         A likelihood functor, likely an instance of :class:`Likelihood`
      :param x0:
         A starting point that satisfies the constraints. If None, then it is initialized by calling
         the `random_starting_point` of `likelihood`.
      :param A: 
         A matrix of linear constraints. If None, then it is initialized by calling the
         `intrinsic_linear_equalities` of `likelihood` with an argument of `True` (e.g. includes
         fixed components a linear equality constraints).
      
      :result:

         A functor object that takes as  argument a numpy array with the dimensions of the
         nullspace and returns the likelihood.  Additionally, the functor contains two attributes:

         - `to_kernel_coords` takes a numpy array as input and maps the configuration space to the
            reduced reduced dimension space
         - `from_kernel_coods` does the opposite, taking a reduced-dimension vector and returns the
            whole vector, with the missing components initialized to those of `x0`.
  """
  from numpy import array, dot, concatenate, zeros, argsort, count_nonzero, arange
  from numpy.linalg import svd, norm
  # get linear constraints
  if A is None: A, b = likelihood.intrinsic_linear_equalities(True)
  if x0 is None: x0 = likelihood.random_starting_point(True)
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
  for r in kernel:
    r[:] = r[:] / norm(r)
    # check sign: we want as many (all?) diagonal elements negative as possible 
    N = likelihood.nstates
    diagonals = arange(N*N) % N == arange(N*N) // N
    if 2 * count_nonzero(r[diagonals] >= 0e0) >= N: r[:] *= -1

  # and create to affine transforms to go back and forth between one set of coordinate and the other
  non_kernel_vec = (x0.flat - dot(kernel.T, dot(kernel, x0.flat))).reshape(x0.shape)
  to_kernel_coords = lambda x: dot(kernel, x.flat)
  from_kernel_coords = lambda x: dot(x, kernel).reshape(x0.shape) + non_kernel_vec

  def kernel_likelihood(x): return likelihood(from_kernel_coords(x))

  kernel_likelihood.to_kernel_coords = to_kernel_coords
  kernel_likelihood.to_kernel_coords.__doc__  \
      = """ Mapping from whole coordinate space to kernel. """
  kernel_likelihood.from_kernel_coords = from_kernel_coords
  kernel_likelihood.from_kernel_coords.__doc__ \
      = """ Mapping from kernel back to whole coordinate space. """
  kernel_likelihood.kernel = kernel
  kernel_likelihood.start = x0
  kernel_likelihood.__doc__ = """ Likelihood in kernel coordinates. """

  return kernel_likelihood


def minimizer_constraints(likelihood, A=None, b=None, infinity=1e5):
  """ Constructs constraints for scipy.optimize.minimize. 
  
      Converts equality constraints and bounds constraints to something scipy.optimize.minimize
      understands.

      :param likelihood:
         An instance of :class:`likelihood`
      :param eq:
         A matrix giving additional equality constraints. It should be of dimension n by m, with n
         the number of constraints and m the size of the problem (ie number of elements in
         Q-matrix). Ignored if None.
      :param b:
         A vector defining additional equality constraints. It should sport the same number of
         elements as there are rows in `A`. It requires `A`. However, if `A` is given and `b` is
         not, `b` will default to vector of zeros of appropriate size.
      :param float infinity: 
         Upper bound for components of the Q-matrix. Can also be numpy.infty (or equivalently,
         None), although that will most likely have violent consequences.
      :returns: `(equalities, bounds)`, where equalities is a dictionary identifying the linear
                equality constraints, and 
  """
  from numpy import zeros, infty
  from itertools import chain
  if infinity is None: infinity = infty
  equalities = []
  if A is None and b is not None:
    raise ValueError('b requires A on input to minimizer_contraints.')
  if A is not None and b is None: b = zeros(A.shape[0])


  Aintrinsic, bintrinsic = likelihood.intrinsic_linear_equalities(True)
  if A is not None and A.shape[1] != Aintrinsic.shape[1]:
    raise ValueError('A has incorrect number of columns')
    
  def eq_function(x, vector, value, *args, **kwargs):
    from numpy import dot
    return dot(vector, x) - value
  def eq_jacobian(x, vector, *args, **kwargs): return vector

  iterator = zip(Aintrinsic, bintrinsic)
  if A is not None: iterator = chain(zip(A, b), iterator)
  for row, value in iterator:
    equalities.append({ 'type':'eq', 
                        'fun':eq_function, 
                        'jac':eq_jacobian, 
                        'args':[row.copy(), value] }) 

  N = likelihood.nstates
  bounds = []
  for i, fixed in enumerate(likelihood.fixed_mask.flat):
    if fixed: bounds.append((-infinity, infinity))
    elif i % N == i // N: bounds.append((-infinity, 0))
    else: bounds.append((-infinity, infinity))
  return equalities, bounds
  


