from behave import given, when, then
from test_setup import register_type
register_type()

@given('a list of {n:Integer} random missed-events likelihoods with tau={tau:Float} '             \
       'and nmax={nmax:Integer}')
def step(context, n, tau, nmax):
  from dcprogs.random import qmatrix as random_qmatrix
  from dcprogs.likelihood import create_missed_eventsG
  qmatrices, Gs, i = [], [], 10*n
  while len(Gs) != n:
    i -= 1
    if i == 0: raise AssertionError('Could not instanciate enough likelihoods.')
    qmatrix = random_qmatrix()
    try: G = create_missed_eventsG(qmatrix, tau, nmax)
    except: continue
    else:
      Gs.append(G)
      qmatrices.append(qmatrix)
  if not hasattr(context, 'qmatrices'): context.qmatrices = []
  if not hasattr(context, 'likelihoods'): context.likelihoods = []
  context.qmatrices.extend(qmatrices)
  context.likelihoods.extend(Gs)

@when('MissedEventsG objects are instantiated with the q-matrices and tau={tau:Float}'             \
      'and nmax={nmax:Integer}')
def step(context, tau, nmax):
  from dcprogs.likelihood import create_missed_eventsG
  if not hasattr(context, "likelihoods"): context.likelihoods = []
  for i, qmatrix in enumerate(context.qmatrices):
    if qmatrix is None: context.likelihoods.append(None); continue
    try: 
      context.likelihoods.append(create_missed_eventsG(qmatrix, tau, nmax))
    except ArithmeticError: 
      context.likelihoods.append(None)
      context.qmatrices[i] = None
      continue
    except:
      print qmatrix
      raise

@when('the {name} CHS occupancies are computed')
def step(context, name): 
  if not hasattr(context, 'occupancies'): context.occupancies = []
  equname = '{0}_CHS_occupancies'.format(name)
  for G in context.likelihoods:
    function = getattr(G, equname)
    context.occupancies.append([function(t) for t in context.times])
    


@then('{name} is zero if t is between {start:Float} and {end:Float}')
def step(context, name, start, end): 
  from numpy import abs, any
  times = context.times
  times = times[times < end] 
  times = times[times >= start]
  for missed_events_G in context.likelihoods:
    if missed_events_G is None: continue
    for t in times:
      try:
        value = getattr(missed_events_G, name)(t)
        assert any(abs(value) < context.tolerance)
      except:
        print(missed_events_G)
        print("{name}({t}) = {value}".format(name=name, t=t, value=value))
        raise 

@then('{name} can be found from ExactSurvivor if t is between {start:Float} and {end:Float}')
def step(context, name, start, end): 
  from numpy import abs, any, dot
  from scipy.linalg import expm
  times = context.times
  times = times[times < end] 
  times = times[times > start]
  for exact, missed_events_G, qmatrix in zip(context.exact_survivors,
                                             context.likelihoods,
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

@then('{name} can be found from ApproxSurvivor if t is larger than {start:Float}')
def step(context, name, start): 
  from numpy import abs, any, dot
  from scipy.linalg import expm
  times = context.times
  times = times[times > start]
  for approx, missed_events_G, qmatrix in zip(context.approx_survivors,
                                             context.likelihoods,
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

def compute_Hfa(qmatrix, tau, tcrit):
  from dcprogs.likelihood import create_approx_survivor
  from scipy.linalg import expm
  from numpy import exp, dot

  approx = create_approx_survivor(qmatrix, tau)
  result = None
  # This sums ^{F}R_i tau_i e^(-(tcrit - tau) / tau_i )
  for matrix, root in approx.fa_components: 
    if result is None: result = -matrix / root * exp( root * (tcrit - tau) )
    else: result += -matrix / root * exp( root * (tcrit - tau) )
  # multiply by Q_{FA} e^{Q_AA} 
  print "!!!! \n", result
  return dot(result, dot(qmatrix.fa, expm(tau * qmatrix.aa)))


@then('the initial CHS occupancies are the solutions to the CHS equations')
def step(context):
  from numpy import dot, abs, all, any, sum

  didloop = False
  for qmatrix, G, occupancies in zip(context.qmatrices, context.likelihoods, context.occupancies):
    phif = G.final_occupancies

    for t, occ in zip(context.times, occupancies):
      try:  
        Hfa = compute_Hfa(qmatrix, G.tau, t)
        phif_Hfa = dot(phif, Hfa)
        check = phif_Hfa / sum(phif_Hfa[:G.nopen], axis=1)
        assert all(abs(occ - check) < context.tolerance)
        assert any(abs(occ - 2e0*check) > context.tolerance)
      except:
        print(G)
        print("  * occupancies: {0}".format(occ))
        print("  * check: {0}".format(check))
        print("  * Hfa shape: {0}".format(Hfa))
        raise

      didloop = True
  assert didloop

@then('the final CHS occupancies are the solutions to the CHS equations')
def step(context):
  from numpy import abs, all, any, sum

  for qmatrix, G, occupancies in zip(context.qmatrices, context.likelihoods, context.occupancies):
    for t, occ in zip(context.times, occupancies):
      try:  
        Hfa = compute_Hfa(qmatrix, G.tau, t)
        check = sum(Hfa[:, :G.nopen], axis=1) 
        assert all(abs(occ - check) < context.tolerance)
        assert any(abs(occ - 2e0*check) > context.tolerance)
      except:
        print(G)
        print("  * occupancies: {0}".format(occ))
        print("  * check: {0}".format(check))
        print("  * Hfa shape: {0}".format(Hfa.shape))
        raise
