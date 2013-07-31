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

from behave import given, when, then
from test_setup import register_type
register_type()

@given('a {matrix:Matrix}, {nopen:Integer}, and {tau:Float}')
def step(context, matrix, nopen, tau):
  context.matrix = matrix
  context.nopen = nopen
  context.tau = tau

@given('a determinantal equation instantiated with the {model} Q matrix')
def setp(context, model):
  from test_setup import DetModel
  context.determinant = DetModel(model)

@given('the output matrices below')
def step(context):
  from re import compile
  pattern = compile('\s*q(\d+)')
  if not hasattr(context, 'output'): context.output = {}
  exec(pattern.sub('\nq\\1', context.text), globals(), context.output)

@when('a determinantal equation is instantiated')
def step(context):
  from sys import exc_info
  from dcprogs.likelihood import DeterminantEq

  print(context.matrix)
  try: context.determinant = DeterminantEq(context.matrix, context.nopen, context.tau)
  except: context.initialization_exception = exc_info() 

@then("equation's tau is {tau:Float}")
def step(context, tau):
  assert abs(context.determinant.tau - tau) < 1e-8

@when("a determinantal equation is instantiated from a state matrix and {tau:Float}")
def step(context, tau):
  from sys import exc_info
  from dcprogs.likelihood import DeterminantEq

  try: context.determinant = DeterminantEq(context.qmatrix, tau)
  except: context.initialization_exception = exc_info() 

@when("The determinantal equation is computed for {s:Float}")
def step(context, s):
  from dcprogs.likelihood import DeterminantEq
  context.result = DeterminantEq(context.matrix, context.nopen, context.tau)(s)

@given("the transition matrix below with {nopen:Integer} open states")
def step(context, nopen):
  from numpy import array
  from dcprogs.likelihood import QMatrix
  from dcprogs import internal_dtype

  matrix = context.text.lstrip().rstrip().splitlines()
  matrix = array([[eval(v) for v in u.split(',')] for u in matrix], dtype=internal_dtype)
  context.qmatrix = QMatrix(matrix, nopen)


@when("the {event}-state determinant is computed for s=({s}) and tau={tau:Float}")
def step(context, event, s, tau):
  import numpy
  from dcprogs.likelihood import DeterminantEq
  s = eval(s, globals().copy(), numpy.__dict__.copy())
  if event == "open": context.result = DeterminantEq(context.qmatrix, tau)(s)
  else:               context.result = DeterminantEq(context.qmatrix.transpose(), tau)(s)

@when('H is computed from the scalar or sequence {input:Eval}')
def step(context, input):
  context.H = context.determinant.H(input)

@then("The result is close to zero ({convergence:Float})")
def step(context, convergence): 
  from numpy import abs, any, array
  if any(abs(array(context.result)) > convergence):
    raise Exception("Result is not zero ({0}).".format(context.result))

@then('H returns the matrix or sequence of matrices {output}')
def step(context, output):
  from numpy import all, abs, array
  output = array(eval(output, globals(), context.output))

  assert all(context.H.shape == output.shape)
  print(abs(context.H - output) <= 1e-8 * abs(output) + 1e-8)
  print(context.H)
  assert all(abs(context.H - output) <= 1e-8 * abs(output) + 1e-8) 
