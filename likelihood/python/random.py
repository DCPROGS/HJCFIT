""" Methods for creating random objects. 

    The difficulty is to create objects that have some structure, and that are not too pathological.
"""
__docformat__ = "restructuredtext en"
__all__ = ['rate_matrix', 'random_idealg', 'time_intervals', 'time_series']
def rate_matrix(N=(5, 10), zeroprob=0.7, large=0.5, factor=1e4, nonsingular=True, realeigs=True,
                tolerance=1e-8):
  """ Creates a random matrix with some structure to it.

      The structure given to the random matrix is as follows:

        - The matrix is not unlikely to have zero-valued components in it. The zero-valued
          components always come in pairs :math:`matrix[i, j] == matrix[j, i] == 0`.
        - Components are chosen randomly using one of two scales: between 0 and 1, and between 0 and
          `factor`.
        - Each row sums to zero.
        - Off-diagonal elements are positive and real-valued.
        - Off-diagonal elements are negative and real-valued (and non-zero if `nonsingular=True`).

      Furthermore, we can specify that the matrix should not be singular, and that it should have
      real-valued eigenvalues. These garantees are insured by generating random matrices until the
      conditions are satisfied: it might take a while for smaller matrices...

      :param N:
        Size of the resulting square matrix.
        If a 2-tuple `(a, b)` is given, then a random number is chosen in :math:`a \leq N \leq b`. 
      :param zeroprob:
        Probability that values in the matrix are zero. Note that if a value is zero, its transpose
        is zero as well.
      :param large:
        Probability that a non-zero value is large. 
      :param factor:
        Factor by which large values are larger. 
      :param nonsingular:
        If True, matrix will be non-singular. 
      :param realeigs: 
        If True, matrix will have real-valued eigenvalues (imaginary components of eigenvalues will
        be zero).
      :param tolerance: 
        Criteria for non-singularity and real-valued eigenvalues.
  """
  from collections import Sequence
  from numpy import abs, any
  from numpy.random import random, randint
  from numpy.linalg import det, eig
  if isinstance(N, Sequence): N = randint(*N)

  def generate_matrix():
    matrix = random((N, N))
    for i in xrange(N):
      for j in xrange(i+1, N):
        if random(1) < zeroprob: 
          matrix[i, j] = 0
          matrix[j, i] = 0
        elif random(1) > large:
          matrix[(i, j) if random(1) > 0.5  else (j, i)] *= factor 
    for i in xrange(matrix.shape[0]): matrix[i, i] -= sum(matrix[i, :])
    return matrix

  def zero_eig(result): 
    from numpy.linalg import svd
    from numpy import sum, abs
    try: U, sing, V = svd(result)
    except: return False
    else: return sum(abs(sing) < tolerance) == 1

  result = generate_matrix()
  # Loop until conditions are satisfied.
  while (nonsingular and abs(det(result)) < tolerance)                                             \
        or (realeigs and any(abs(eig(result)[0].imag) > tolerance)                                 \
        or not zero_eig(result) ):
    result = generate_matrix()

  return result

def qmatrix(*args, **kwargs):
  """ Creates a random state matrix with some structure to it. """
  from numpy.random import randint
  from .likelihood import QMatrix
  
  def zero_eig(result):
    """ Qff and Qaa cannot be singular. """
    from numpy.linalg import svd
    from numpy import all, any, abs
    try: singular = svd(result.aa)[1]
    except: return False
    if any(abs(singular) < 1e-8): return False
    try: singular = svd(result.ff)[1]
    except: return False
    return all(abs(singular) > 1e-8)


  if 'get_qmatrix' in kwargs: get_qmatrix = kwargs['get_qmatrix']
  else:
    def get_qmatrix():
      matrix = rate_matrix(*args, **kwargs)
      nopen = randint(2, matrix.shape[0]-2)
      return QMatrix(matrix, nopen)

  result = get_qmatrix()
  while not zero_eig(result): result = get_qmatrix()

  return result
# Adds description of parameters and function from the rate_matrix docstring.
qmatrix.__doc__ = "\n".join(qmatrix.__doc__.splitlines()
                                 + rate_matrix.__doc__.splitlines()[1:])


def random_idealg(*args, **kwargs):
  """ Creates a random state matrix with some structure to it. """
  from .likelihood import IdealG
  return IdealG(qmatrix(*args, **kwargs))


def time_intervals(N=100, n=100, nsmall=3, maxsize=10, tau=1e-4):
  """ Create a set of random time intervals.
 
      :param integer N: 
        Number of time intervals once intervals smaller than `tau` have been filtered out.
      :param integer n: 
        Number of time intervals strictly smaller than :math:`\\tau`
      :param integer nsmall:
        When splitting intervals, and in so far as it is possible, will create a cluster of
        :math:`n \in [0, \mathrm{nsmall}]` small intervals.
      :param integer maxsize:
        Maximum size of the large intervals, in units of :math:`\\tau`.
 
      :result: 
        A 2-tuple, where the first item is the list of intervals after filter, and the second item
        the intervals including those smaller than :math:`\\tau`.
  """ 
  from numpy.random import uniform, randint, shuffle
  from numpy import where, sum, concatenate

  perfect = uniform(3*tau + 1e-4, maxsize*tau, N)  
  result = perfect.copy()
  
  while len(result) < N + n:
    # Figure out where small intervals can go
    indices = where(result > 2 * tau + 1e-8)[0]
    if len(indices) == 0: raise Exception('No intervals large enough to split...')
    # Pick one large interval to split
    index = indices[randint(0, len(indices)-1)]
    # Compute number of small intervals (< tau ) that will be added
    nb_small_intervals = min(nsmall, int((result[index] - tau - 1e-8)/tau))
    nb_small_intervals = min(nb_small_intervals, n - len(result) + N)
    if nb_small_intervals in [0, 1]: nb_small_intervals = 1
    else: nb_small_intervals = randint(1, nb_small_intervals) 
    if nb_small_intervals % 2 == 0: nb_small_intervals += 1
    # Create split intervals
    t0 = uniform(tau + 1e-8, result[index] - tau)
    current_size = result[index] - t0
    small_intervals = []
    while len(small_intervals) < nb_small_intervals:
      small_intervals.append(uniform(0, min(tau, current_size)))
      current_size -= small_intervals[-1]
    t1 = result[index] - sum(small_intervals) - t0
    shuffle(small_intervals)
    # Put it all back together 
    result = concatenate(( result[:index],
                           [t0], 
                           small_intervals, 
                           [t1],
                           result[index+1:] ))

  return perfect, result

def time_series(*args, **kwargs): 
  """ Create a random time series. """
  from ._likelihood_methods import intervals_to_series
  start = kwargs.pop('start', 0)
  perfect, result = time_intervals(*args, **kwargs)
  return intervals_to_series(perfect, start), intervals_to_series(result, start)

time_series.__doc__ = "\n".join( time_series.__doc__.splitlines() 
                                 + time_intervals.__doc__
                                                 .replace('intervals', 'times')
                                                 .splitlines()[1:] )
