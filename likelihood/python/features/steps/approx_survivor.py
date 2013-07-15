from behave import when
from test_setup import register_type
register_type()

@when('ApproxSurvivor objects are instantiated with the q-matrices and tau={tau:Float}')
def steps(context, tau):
  from dcprogs.likelihood import create_approx_survivor
  if not hasattr(context, "approx_survivors"): context.approx_survivors = []
  for i, qmatrix in enumerate(context.qmatrices):
    if qmatrix is None: context.approx_survivors.append(None); continue
    try: 
      context.approx_survivors.append(create_approx_survivor(qmatrix, tau))
    except ArithmeticError: 
      context.approx_survivors.append(None)
      context.qmatrices[i] = None
      continue
    except:
      print qmatrix
      raise

@given('a list of {n:Integer} random approximate survivor functions with tau={tau:Float}')
def step(context, n, tau):
  from dcprogs.random import qmatrix as random_qmatrix
  from dcprogs.likelihood import create_approx_survivor
  qmatrices, survivors, i = [], [], 10*n
  while len(survivors) != n:
    i -= 1
    if i == 0: raise AssertionError('Could not instanciate enough survivor functions.')
    qmatrix = random_qmatrix()
    try: survivor = create_approx_survivor(qmatrix, tau)
    except: continue
    else: survivors.append(survivor)
  if not hasattr(context, 'qmatrices'): context.qmatrices = []
  if not hasattr(context, 'survivors'): context.survivors = []
  context.qmatrices.extend(qmatrices)
  context.survivors.extend(survivors)

@when('the {name} block of the survivor function is computed for t={t:Float}')
def step(context, name, t):
  if not hasattr(context, 'survivor_values'): context.survivor_values = []
  for survivor in context.survivors:
    function = getattr(survivor, name) # gets function for particular block. 
    context.survivor_values.append( function(t) )

@then('the value is the identity')
def step(context):
  from numpy import all, abs, identity
  for survivor, qmatrix, value in zip(context.survivors, context.qmatrices,
                                      context.survivor_values):
    try: 
      assert all(abs(value - identity(value.shape[0])) < context.tolerance)
    except:
      print survivor
      print qmatrix
      raise
