from behave import when
from test_setup import register_type
register_type()
@when('ExactSurvirvor objects are instantiated with the q-matrices and tau={tau:Float}')
def steps(context, tau):
  from dcprogs.likelihood import ExactSurvivor
  if not hasattr(context, "exact_survivors"): context.exact_survivors = []
  for i, qmatrix in enumerate(context.qmatrices):
    if qmatrix is None: context.exact_survivors.append(None); continue
    try: 
      context.exact_survivors.append(ExactSurvivor(qmatrix, tau))
    except ArithmeticError: 
      context.exact_survivors.append(None)
      context.qmatrices[i] = None
      continue
    except:
      print qmatrix
      raise
