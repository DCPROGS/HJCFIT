from behave import given, when, then
from test_setup import register_type
register_type()

@given("StateMatrix is accessible")
def step(context):
  from dcprogs.likelihood import StateMatrix
  context.StateMatrix = StateMatrix

@when("we instantiate StateMatrix without arguments")
def step(context):
  context.statematrix = context.StateMatrix()

@when('we instantiate StateMatrix with {matrix:Matrix} and {nopen:Integer}')
def step(context, matrix, nopen):
  from sys import exc_info
  try: context.statematrix = context.StateMatrix(matrix, nopen)
  except: context.initialization_exception = exc_info()

@given('a StateMatrix instantiated with {matrix:Matrix} and {nopen:Integer}')
def step(context, matrix, nopen):
  from dcprogs.likelihood import StateMatrix
  context.statematrix = StateMatrix(matrix, nopen)

@then('nopen is {nopen:Integer}')
def step(context, nopen):
  assert context.statematrix.nopen == nopen

@then("matrix is a numpy array")
def step(context):
  from numpy import ndarray
  assert isinstance(context.statematrix.matrix, ndarray)

@then("matrix is {matrix:Matrix}")
def step(context, matrix):
  from numpy import all, abs
  assert len(context.statematrix.matrix.shape) == len(matrix.shape)
  assert context.statematrix.matrix.shape == matrix.shape
  assert all(abs(context.statematrix.matrix - matrix) < 1e-8)

@then('{block} is given by rows="{rows}" and cols="{cols}" of {matrix:Matrix}, ' \
      'and nopen={nopen:Integer}')
def step(context, block, rows, cols, matrix, nopen):
  from numpy import all, abs
  block = getattr(context.statematrix, block)
  rows = rows.replace('end', str(matrix.shape[0]))
  cols = cols.replace('end', str(matrix.shape[1]))
  check = eval("input[{rows}, {cols}]".format(rows=rows, cols=cols),
               {'input': matrix, 'nopen': nopen, 'begin': 0})
  assert all(abs(block - check) < 1e-8)
  

@When("item {item:Eval} is set to {value:Float}")
def step(context, item, value):
  context.statematrix.matrix[item] = value
@Then("item {item:Eval} is {value:Float}")
def step(context, item, value):
  from numpy import abs
  assert abs(context.statematrix.matrix[item] - value) < 1e-8
@When("nopen is set to {value:Integer}")
def step(context, value):
  context.statematrix.nopen = value

@given('a list of {n:Integer} random state matrices')
def step(context, n):
  from dcprogs.random import state_matrix as random_state_matrix
  context.matrices = [random_state_matrix() for u in range(n)]
