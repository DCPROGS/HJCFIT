from behave import given, when, then, matchers

# First add some parsers.
def Matrix(string): 
  """ Creates matrices from specific strings """
  from numpy import array, identity
  if string == "Qmatrix":
    return array([[ -3050,        50,  3000,      0,    0 ], 
                  [ 2./3., -1502./3.,     0,    500,    0 ],  
                  [    15,         0, -2065,     50, 2000 ],  
                  [     0,     15000,  4000, -19000,    0 ],  
                  [     0,         0,    10,      0,  -10 ] ])
  if string == "empty": return array([])
  if string == "spam": return identity(5).tolist()
  if string == "numpy_spam": return array([['a', 'b', 'c']*3])

  
matchers.register_type(Integer=lambda x: int(x))
matchers.register_type(Float=lambda x: float(x))
matchers.register_type(Eval=lambda x: eval(x))
matchers.register_type(Matrix=Matrix)

@given("StateMatrix is accessible")
def step(context):
  from dcprogs.likelihood import StateMatrix
  context.StateMatrix = StateMatrix

@when("we instantiate StateMatrix without arguments")
def step(context):
  context.statematrix = context.StateMatrix()

@when('we instantiate StateMatrix with {matrix:Matrix} and {nopen:Integer}')
def step(context, matrix, nopen):
  try: context.statematrix = context.StateMatrix(matrix, nopen)
  except (TypeError, ValueError) as e: context.initialization_exception = e

@given('a StateMatrix instantiated with {matrix:Matrix} and {nopen:Integer}')
def step(context, matrix, nopen):
  from dcprogs.likelihood import StateMatrix
  context.statematrix = StateMatrix(matrix, nopen)

@then('instantiation did not throw')
def step(context):
  assert not hasattr(context, 'initialization_exception')

@then('instantiation threw {type}')
def step(context, type):
  assert hasattr(context, 'initialization_exception')
  exception = getattr(context, 'initialization_exception')
  type = eval(type)
  assert isinstance(exception, type)


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
