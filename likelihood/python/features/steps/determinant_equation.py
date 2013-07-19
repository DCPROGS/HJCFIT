from behave import given, when, then
from test_setup import register_type
register_type()

@given('a {matrix:Matrix}, {nopen:Integer}, and {tau:Float}')
def step(context, matrix, nopen, tau):
  context.matrix = matrix
  context.nopen = nopen
  context.tau = tau

@when('a determinantal equation is instantiated')
def step(context):
  from sys import exc_info
  from dcprogs.likelihood import DeterminantEq

  print context.matrix
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

  matrix = context.text.lstrip().rstrip().splitlines()
  matrix = array([[eval(v) for v in u.split(',')] for u in matrix], dtype='float64')
  context.qmatrix = QMatrix(matrix, nopen)


@when("the {event}-state determinant is computed for s=({s}) and tau={tau:Float}")
def step(context, event, s, tau):
  import numpy
  from dcprogs.likelihood import DeterminantEq
  s = eval(s, globals().copy(), numpy.__dict__.copy())
  if event == "open": context.result = DeterminantEq(context.qmatrix, tau)(s)
  else:               context.result = DeterminantEq(context.qmatrix.transpose(), tau)(s)

@then("The result is close to zero ({convergence:Float})")
def step(context, convergence): 
  from numpy import abs, any, array
  if any(abs(array(context.result)) > convergence):
    raise Exception("Result is not zero ({0}).".format(context.result))
