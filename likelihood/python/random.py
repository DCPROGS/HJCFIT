def qmatrix(N=(5, 10), zeroprob=0.7, large=0.5, factor=1e4, nonsingular=True, realeigs=True,
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

  def equilibrium(result): 
    from numpy.linalg import svd
    from numpy import sum, abs
    try: U, sing, V = svd(result)
    except: return False
    else: return sum(abs(sing) < tolerance) == 1

  result = generate_matrix()
  # Loop until conditions are satisfied.
  while (nonsingular and abs(det(result)) < tolerance)                                             \
        or (realeigs and any(abs(eig(result)[0].imag) > tolerance)                                 \
        or not equilibrium(result) ):
    result = generate_matrix()

  return result

def state_matrix(*args, **kwargs):
  """ Creates a random state matrix with some structure to it. """
  from numpy.random import randint
  from .likelihood import StateMatrix
  matrix = qmatrix(*args, **kwargs)
  nopen = randint(2, matrix.shape[0]-2)
  return StateMatrix(matrix, nopen)

def random_idealg(*args, **kwargs):
  """ Creates a random state matrix with some structure to it. """
  from .likelihood import IdealG
  return IdealG(state_matrix(*args, **kwargs))

# Adds description of parameters and function from the qmatrix docstring.
state_matrix.__doc__ = "\n".join(state_matrix.__doc__.splitlines()
                                 + qmatrix.__doc__.splitlines()[1:])


