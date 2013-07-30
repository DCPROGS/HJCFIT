########################
#   DCProgs computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

""" Subpackage for likelihood optimization. """
__docformat__ = "restructuredtext en"
__all__ = ['reduce_likelihood']

def reduce_likelihood(likelihood, graph_matrix):
  """ Maps likelihood to a set of variable components.
  
      The goal is a callable that takes on input a numpy array with only variable components.
      It hides from the input the components that are fixed or can be obtained as a result of an
      expression.

      :param likelihood:
         This should be a callable object that takes a :class:`QMatrix` or numpy matrix on input.

      :param graph_matrix:
        Defines the reaction graph. This is a list of list (matrix) where each element is either
        "V" (variable), 0 (no direct reaction between these two states), a number (Fixed
        reaction rate), or a string with a valid python expression (that `eval` understands). In the
        latter case, `q` will be replaced with the value of the qmatrix, `i` is set to the current
        row index, and `j` to the current column index. The open-states should be in the top-left
        corner.
        
        .. code-bloc:: python
        
           [ ["V",   0,           0.1], 
             [  0, "V", "2e0*q[2, 1]"],
             ["V", "V",           "V"] ]

         In the exampl above, there are no direct reactions transforming state (1) into state
         (2), and the reaction rate from state (1) to state (3) is fixed. All others will be
         optimized, subject to the intrinsic constraints (sum over each row is zero, diagonal
         elements are negative) and extrinsic contraints (defined below).
        
         Note that diagonal elements are always obtained from the condition that rows of the matrix
         should be equal to one (e.g. as though they were set to "k[i, i] - sum(k[i])"). This
         constraint is the last one imposed on the qmatrix.

         The expression mechanism above is not sufficient to handle cycles, e.g `q[1, 2]` is an
         expression that depends on `q[2, 1]`, which itself is an expression that depends on 
         `q[1, 2]`. Such cases will not produce an error. Howerver, their result is undefined.
         Furthermore, it should not explicitely depend on the diagonal components. Those components
         have not yet been constrained to have rows sum to zero.

     :returns:
     
        A callable from which the fixed components have been abstracted.
        It takes on input a numpy vector with as many components as there are variable components in
        graph_matrix.

        For convenience, the callable has a `to_reduced_coords` method which takes a numpy matrix
        and returns a vector with only the variable components.
        It also sports a  `to_full_coords` coords that maps back to the whole space.
  """
  import numpy
  from numpy import array, zeros, sum
  from .likelihood import QMatrix
  
  nstates = len(graph_matrix)
  """ Number of states in mechanism. """
  if any(len(u) != nstates for u in graph_matrix):
    raise ValueError('graph_matrix should be square')
  
  fixed = zeros((nstates, nstates), dtype='bool')
  """ Array indexing fixed components of the mechanism. """
  variable = zeros((nstates, nstates), dtype='bool')
  """ Array indexing variable components of the mechanism. """
  expressions = []
  """ List of expression used to set matrix components. """
  for i, row in enumerate(graph_matrix):
    for j, value in enumerate(row):
      if i == j: continue
      if not isinstance(value, str): fixed[i, j] = True
      elif value.lower().rstrip().lstrip() == 'v': variable[i, j] = True
      else: expressions.append((i, j, value))
  if expressions is None: expressions = None

  # sanity check
  for i, row in enumerate(graph_matrix):
    for j, value in enumerate(row):
      if fixed[i, j] == False: continue
      if abs(value) > 1e-8: continue
      if fixed[j, i] == False or abs(graph_matrix[j][i]) > 1e-8:
        raise ValueError( 'No reversability for transformations from states {0} to {1}.'         \
                          .format(i, j) )

  # This matrix will hold fixed components.
  # During optimization, it will receive the variable components
  qmatrix = QMatrix(zeros((nstates, nstates)), likelihood.nopen)
  for i, row in enumerate(graph_matrix):
    for j, value in enumerate(row):
      if i == j: qmatrix[i, i] = 0
      elif fixed[i, j]: qmatrix[i, j] = value


  # Finally, creates actual functions.
  def to_full_coords(vector): 
    qmatrix.matrix[variable] = vector
    if expressions is not None:
      local_dict = numpy.__dict__.copy()
      global_dict = globals()
      local_dict['q'] = qmatrix
      for i, j, expression in expressions: 
        local_dict['i'] = i
        local_dict['j'] = j
        qmatrix.matrix[i, j] = eval(expression, global_dict, local_dict)
    for i, row in enumerate(qmatrix.matrix):
      row[i] = 0
      row[i] = -sum(row)
    return qmatrix.matrix

  def reduced_likelihood(vector):
    to_full_coords(vector)
    return likelihood(qmatrix)
  
  reduced_likelihood.to_reduced_coords = lambda x: array(x)[variable]
  reduced_likelihood.to_full_coords = lambda x: to_full_coords(x).copy()
  reduced_likelihood.likelihood = likelihood

  return reduced_likelihood
