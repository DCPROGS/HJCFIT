########################
#   DCProgs computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

from behave import given, when, then
from test_setup import register_type
register_type()

@given('a list of {n:Integer} random missed-events likelihoods with tau={tau:Float} '             \
       'and nmax={nmax:Integer}')
def step(context, n, tau, nmax):
  from dcprogs.likelihood.random import qmatrix as random_qmatrix
  from dcprogs.likelihood import MissedEventsG
  qmatrices, Gs, i = [], [], 10*n
  while len(Gs) != n:
    i -= 1
    if i == 0: raise AssertionError('Could not instanciate enough likelihoods.')
    qmatrix = random_qmatrix()
    G = MissedEventsG(qmatrix, tau, nmax)
    try: G = MissedEventsG(qmatrix, tau, nmax)
    except: continue
    else:
      Gs.append(G)
      qmatrices.append(qmatrix)
  if not hasattr(context, 'qmatrices'): context.qmatrices = []
  if not hasattr(context, 'likelihoods'): context.likelihoods = []
  context.qmatrices.extend(qmatrices)
  context.likelihoods.extend(Gs)

@given('the {model} missed-events likelihood') 
def step(context, model):
  from test_setup import eG
  context.eG = eG(model)


@when('MissedEventsG objects are instantiated with the q-matrices and tau={tau:Float}'             \
      'and nmax={nmax:Integer}')
def step(context, tau, nmax):
  from dcprogs.likelihood import MissedEventsG
  if not hasattr(context, "likelihoods"): context.likelihoods = []
  for i, qmatrix in enumerate(context.qmatrices):
    if qmatrix is None: context.likelihoods.append(None); continue
    try: 
      context.likelihoods.append(MissedEventsG(qmatrix, tau, nmax))
    except ArithmeticError: 
      context.likelihoods.append(None)
      context.qmatrices[i] = None
      continue
    except:
      print(qmatrix)
      raise

@when('the {name} CHS occupancies are computed')
def step(context, name): 
  if not hasattr(context, 'occupancies'): context.occupancies = []
  equname = '{0}_CHS_occupancies'.format(name)
  for G in context.likelihoods:
    function = getattr(G, equname)
    context.occupancies.append([function(t) for t in context.times])
    
@when('we compute the {which} CHS occupancies with tcrit={tcrit:Float}')
def step(context, which, tcrit):
  context.chs = getattr(context.eG, '{0}_CHS_occupancies'.format(which))(tcrit)





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
  from dcprogs.likelihood import expm
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
  from dcprogs.likelihood import expm
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
        raise 

def compute_Hfa(qmatrix, tau, tcrit):
  from dcprogs.likelihood import ApproxSurvivor
  from dcprogs.likelihood import expm
  from numpy import exp, dot

  approx = ApproxSurvivor(qmatrix, tau)
  result = None
  # This sums ^{F}R_i tau_i e^(-(tcrit - tau) / tau_i )
  for matrix, root in approx.fa_components: 
    if result is None: result = -matrix / root * exp( root * (tcrit - tau) )
    else: result += -matrix / root * exp( root * (tcrit - tau) )
  # multiply by Q_{FA} e^{Q_AA} 
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
        check = phif_Hfa / sum(phif_Hfa[:G.nopen])
        assert all(abs(occ - check) < context.tolerance)
        assert any(abs(occ - 2e0*check) > context.tolerance)
        assert abs(sum(occ) - 1e0) < context.tolerance
      except:
        print(G)
        print("  * occupancies: {0}".format(occ))
        print("  * check: {0}".format(check))
        print("  * Hfa shape: {0}".format(Hfa.shape))
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
        if any(abs(check) > 5e1 * context.tolerance): 
          assert any(abs(occ - 2e0*check) > context.tolerance)
      except:
        print(G)
        print("  * occupancies: {0}".format(occ))
        print("  * check: {0}".format(check))
        print("  * Hfa shape: {0}".format(Hfa.shape))
        print("  * time: {0}".format(t))
        raise

@then('the {which} CHS occupancies compare to {prior:Eval}')
def step(context, which, prior):
  from numpy import all, abs

  assert all(abs(context.chs - prior) < context.tolerance)
