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
