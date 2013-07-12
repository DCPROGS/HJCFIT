from behave import when, then
from test_setup import register_type
register_type()
@when('MissedEventsG objects are instantiated with the q-matrices and tau={tau:Float}' \
      'and nmax={nmax:Integer}')
def step(context, tau, nmax):
  from dcprogs.likelihood import create_missed_eventsG
  if not hasattr(context, "missed_events_Gs"): context.missed_events_Gs = []
  for i, qmatrix in enumerate(context.qmatrices):
    if qmatrix is None: context.missed_events_Gs.append(None); continue
    try: 
      context.missed_events_Gs.append(create_missed_eventsG(qmatrix, tau, nmax))
    except ArithmeticError: 
      context.missed_events_Gs.append(None)
      context.qmatrices[i] = None
      continue
    except:
      print qmatrix
      raise

@then('MissedEventsG.{name}(t) is zero if t is between {start:Float} and {end:Float}')
def step(context, name, start, end): 
  from numpy import abs, any
  times = context.times
  times = times[times < end] 
  times = times[times >= start]
  for missed_events_G in context.missed_events_Gs:
    if missed_events_G is None: continue
    for t in times:
      value = getattr(missed_events_G, name)(t)
      if any(abs(value) > context.tolerance): 
        print(missed_events_G)
        print("{name}({t}) = {value}".format(name=name, t=t, value=value))
        raise AssertionError()

@then('MissedEventsG.{name}(t) can be found from ExactSurvivor if t is between '  \
      '{start:Float} and {end:Float}')
def step(context, name, start, end): 
  from numpy import abs, any, dot
  from scipy.linalg import expm
  times = context.times
  times = times[times < end] 
  times = times[times > start]
  for exact, missed_events_G, qmatrix in zip(context.exact_survivors,
                                             context.missed_events_Gs,
                                             context.qmatrices):
    if missed_events_G is None: continue
    tau = missed_events_G.tau
    if name == "af": factor = dot(qmatrix.af, expm(tau * qmatrix.ff))
    else: factor = dot(qmatrix.fa, expm(tau * qmatrix.aa))
    for t in times:
      try: 
        value = getattr(missed_events_G, name)(t)
        check = dot(getattr(exact, name)(t-tau), factor)
        assert any(abs(value - check) < context.tolerance)
      except:
        print(missed_events_G)
        print("  * tmax: {0}".format(missed_events_G.tmax))
        print("{name}({t}*tau): \n{value}".format(name=name, t=t/tau, value=value))
        print("exact.af({t}*tau) * factor: \n{value}".format(name=name, t=(t-tau)/tau, value=check))
        print("factor:\n{factor}".format(factor=factor))
        raise AssertionError()

@then('MissedEventsG.{name}(t) can be found from ApproxSurvivor if t is larger than '  \
      '{start:Float}')
def step(context, name, start): 
  from numpy import abs, any, dot
  from scipy.linalg import expm
  times = context.times
  times = times[times > start]
  for approx, missed_events_G, qmatrix in zip(context.approx_survivors,
                                             context.missed_events_Gs,
                                             context.qmatrices):
    if missed_events_G is None: continue
    tau = missed_events_G.tau
    if name == "af": factor = dot(qmatrix.af, expm(tau * qmatrix.ff))
    else: factor = dot(qmatrix.fa, expm(tau * qmatrix.aa))
    for t in times:
      try: 
        value = getattr(missed_events_G, name)(t)
        check = dot(getattr(approx, name)(t-tau), factor)
        assert any(abs(value - check) < context.tolerance)
      except:
        print(missed_events_G)
        print("  * tmax: {0}".format(missed_events_G.tmax))
        print("{name}({t}*tau): \n{value}".format(name=name, t=t/tau, value=value))
        print("approx.af({t}*tau) * factor: \n{value}".format(name=name, t=(t-tau)/tau, value=check))
        print("factor:\n{factor}".format(factor=factor))
        raise AssertionError()
