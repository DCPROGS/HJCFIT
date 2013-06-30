from behave import given, when, then
from test_setup import register_type
register_type()

@given('the {doopen}-states determinantal equation \"{statmat:StatMat}\" with tau={tau:Float}')
def step(context, doopen, statmat, tau):
  from dcprogs.likelihood import DeterminantEq
  context.equation = DeterminantEq(statmat, tau, doopen == "open")

@when("the root intervals are computed")
def step(context):
  from sys import exc_info
  from dcprogs.likelihood import find_root_intervals
  try: context.root_intervals = find_root_intervals(context.equation)
  except: context.exception = exc_info() 


@then("there are {n:Integer} intervals")
def step(context, n):
  assert len(context.root_intervals) == 2

@then("interval {n:Integer} contains {x:Float}")
def step(context, n, x):
  assert x >= context.root_intervals[n][0][0]
  assert x <= context.root_intervals[n][0][1]
