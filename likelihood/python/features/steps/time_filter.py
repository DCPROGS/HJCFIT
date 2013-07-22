from behave import given, when, then
from test_setup import register_type
register_type()

@given('{n:Integer} random time series with tau={tau:Float}')
def steps(context, n, tau):
  from dcprogs.random import time_series
  from numpy.random import randint

  if not hasattr(context, 'series'): context.series = []
  Ns = randint(10, 150, n)
  ns = [randint(5, u+10) for u in Ns]

  context.series.extend([(time_series(N=N, n=n, tau=tau), tau) for N, n in zip(Ns, ns)])

@when('each time series is filtered')
def step(context):
  from dcprogs.likelihood import time_filter
  
  if not hasattr(context, 'filtered'): context.filtered = []
  for (filtered, series), tau in context.series: 
    context.filtered.append(time_filter(series, tau))

@then('expected filtered time series is obtained')
def step(context):
  from numpy import all, abs  

  for ((check, series), tau), filtered in zip(context.series, context.filtered):
    try: assert all(abs(check - filtered) < 1e-8)
    except: 
      print repr(check)
      print repr(series)
      raise
