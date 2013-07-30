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

from behave import when
from test_setup import register_type
register_type()

@given('the graph matrix below')
def steps(context):
  from numpy import zeros
  context.graph_matrix = eval(context.text)
  N = len(context.graph_matrix)
  context.fixed = zeros((N, N), dtype='bool')
  context.fixed_components = []
  context.variable = zeros((N, N), dtype='bool')
  context.expressions = []
  for i, row in enumerate(context.graph_matrix):
    for j, value in enumerate(row):
      if i == j: continue
      if not isinstance(value, str):
        context.fixed[i, j] = True
        context.fixed_components.append(value)
      elif value.lower() == 'v': context.variable[i, j] = True
      else: context.expressions.append((i, j, value))

@when('the reduced likelihood is created')
def step(context):
  from dcprogs.likelihood.optimization import reduce_likelihood

  class MockLikelihood(object):
    def __init__(self, nopen=2): 
      self.nopen = nopen
    def __call__(self, qmatrix):
      self.qmatrix = qmatrix 

  context.reduced_likelihood = reduce_likelihood(MockLikelihood(), context.graph_matrix)
  
@when('the reduced likelihood called with a random vector')
def step(context):
  from numpy.random import uniform
  from numpy import count_nonzero
  context.vector = uniform(size=(count_nonzero(context.variable)))
  context.reduced_likelihood(context.vector)

@then('the sum of rows of the resulting matrix is zero')
def step(context):
  from numpy import sum, abs
  qmatrix = context.reduced_likelihood.likelihood.qmatrix
  assert qmatrix.nopen == 2
  for i, row in enumerate(qmatrix.matrix):
    assert abs(sum(row)) < 1e-8 

@then('the fixed components of the resulting matrix are those above')
def step(context):
  from numpy import abs
  qmatrix = context.reduced_likelihood.likelihood.qmatrix
  assert all(abs(qmatrix.matrix[context.fixed] - context.fixed_components) < 1e-8)

@then('the variable components of the resulting matrix are from the random matrix')
def step(context):
  from numpy import abs, all
  qmatrix = context.reduced_likelihood.likelihood.qmatrix
  assert all(abs(qmatrix.matrix[context.variable] - context.vector) < 1e-8)

@then('the two expression components come out as expected')
def step(context):
  from numpy import abs, cos
  qmatrix = context.reduced_likelihood.likelihood.qmatrix
  assert abs(qmatrix[1, 3] - 2e0 * qmatrix.matrix[0, 2]) < 1e-8
  assert abs(qmatrix[4, 2] - qmatrix.matrix[0, 2] - 3.0 * cos(qmatrix[2, 4])) < 1e-8
