from behave import given, when, then
from test_setup import register_type
register_type()

@given("QMatrix is accessible")
def step(context):
  from dcprogs.likelihood import QMatrix
  context.QMatrix = QMatrix

@when("we instantiate QMatrix without arguments")
def step(context):
  context.qmatrix = context.QMatrix()

@when('we instantiate QMatrix with {matrix:Matrix} and {nopen:Integer}')
def step(context, matrix, nopen):
  from sys import exc_info
  try: context.qmatrix = context.QMatrix(matrix, nopen)
  except: context.initialization_exception = exc_info()

@given('a QMatrix instantiated with {matrix:Matrix} and {nopen:Integer}')
def step(context, matrix, nopen):
  from dcprogs.likelihood import QMatrix
  context.qmatrix = QMatrix(matrix, nopen)

@then('nopen is {nopen:Integer}')
def step(context, nopen):
  assert context.qmatrix.nopen == nopen

@then("matrix is a numpy array")
def step(context):
  from numpy import ndarray
  assert isinstance(context.qmatrix.matrix, ndarray)

@then("matrix is {matrix:Matrix}")
def step(context, matrix):
  from numpy import all, abs
  assert len(context.qmatrix.matrix.shape) == len(matrix.shape)
  assert context.qmatrix.matrix.shape == matrix.shape
  assert all(abs(context.qmatrix.matrix - matrix) < 1e-8)

@then('{block} is given by rows="{rows}" and cols="{cols}" of {matrix:Matrix}, ' \
      'and nopen={nopen:Integer}')
def step(context, block, rows, cols, matrix, nopen):
  from numpy import all, abs
  block = getattr(context.qmatrix, block)
  rows = rows.replace('end', str(matrix.shape[0]))
  cols = cols.replace('end', str(matrix.shape[1]))
  check = eval("input[{rows}, {cols}]".format(rows=rows, cols=cols),
               {'input': matrix, 'nopen': nopen, 'begin': 0})
  assert all(abs(block - check) < 1e-8)
  

@When("item {item:Eval} is set to {value:Float}")
def step(context, item, value):
  context.qmatrix.matrix[item] = value
@Then("item {item:Eval} is {value:Float}")
def step(context, item, value):
  from numpy import abs
  assert abs(context.qmatrix.matrix[item] - value) < 1e-8
@When("nopen is set to {value:Integer}")
def step(context, value):
  context.qmatrix.nopen = value

@given('a list of {n:Integer} random q-matrices')
def step(context, n):
  from dcprogs.random import qmatrix as random_qmatrix
  context.matrices = [random_qmatrix() for u in range(n)]
