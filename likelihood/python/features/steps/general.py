from behave import then

@given('a list of random of {n:Integer} {name} between {start:Float} and {end:Float}')
def step(context, n, name, start, end):
  from numpy import random
  setattr(context, name, random.sample(n) * (end - start) + start)
@given('a parameter {name}={value:Eval}')
def step(context, name, value):
  setattr(context, name, value)

@then('instantiation did not throw')
def step(context):
  if hasattr(context, 'initialization_exception'):
    raise context.initialization_exception[1]

@then('instantiation threw {type}')
def step(context, type):
  assert hasattr(context, 'initialization_exception')
  exception = getattr(context, 'initialization_exception')
  type = eval(type)
  if not isinstance(exception[1], type): raise exception[1]

@then("no exception was thrown")
def step(context):
  if hasattr(context, 'exception'): raise context.exception[1]

