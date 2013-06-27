from behave import then

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

