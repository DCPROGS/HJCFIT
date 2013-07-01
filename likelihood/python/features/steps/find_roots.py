from behave import given, when, then
from test_setup import register_type
register_type()

@given('the {doopen}-states determinantal equation \"{statmat:StatMat}\" with tau={tau:Float}')
def step(context, doopen, statmat, tau):
  from dcprogs.likelihood import DeterminantEq
  if statmat == "random":
    # Avoid pb when eigenvalue problem cannot be solved.
    while True:
      try: context.equation = DeterminantEq(statmat, tau, doopen == "open")
      except: continue
  else: context.equation = DeterminantEq(statmat, tau, doopen == "open")
  print context.equation

@when("the root intervals are computed")
def step(context):
  from dcprogs.likelihood import find_root_intervals
  if hasattr(context, "exception"): raise context.exception
  context.intervals = find_root_intervals(context.equation)

@when("a brute force search for roots is perfomed with resolution={resolution:Float}")
def step(context, resolution):
  from numpy import min
  from dcprogs.likelihood import find_root_intervals_brute_force, find_root_intervals
  if hasattr(context, "exception"): raise context.exception
  roots = find_root_intervals(context.equation)
  mini = min([r[0] for r in roots])
  context.intervals_brute_force = find_root_intervals_brute_force(context.equation, resolution, 2*mini)
  context.resolution = resolution

@then("there are {n:Integer} intervals")
def step(context, n):
  assert len(context.intervals) == n

@then("interval {n:Integer} contains {x:Float}")
def step(context, n, x):
  assert x >= context.intervals[n][0][0]
  assert x <= context.intervals[n][0][1]

@then("the intervals larger than {resolution:Float} do overlap")
def step(context, resolution):

  eigsearch = [r[0] for r in context.intervals if r[0][1]- r[0][0] > resolution]
  brusearch = [r[0] for r in context.intervals_brute_force]

  eigsearch = sorted(eigsearch)
  brusearch = sorted(eigsearch)

  if len(eigsearch) != len(brusearch):
    raise AssertionError( "Different number of intervals found be methods.\n"
                          + str(context.equation) )
  
  for a, b in zip(eigsearch, brusearch):
    if not (b[0] < a[1] and a[0] < b[1]):
      raise AssertionError("Intervals do not overlap\n" + str(context.equation))

@then("roots are found within each overlap for which sign changes")
def step(context):
  from dcprogs.likelihood import find_roots

  eigsearch = [r[0] for r in context.intervals if r[0][1]- r[0][0] > context.resolution]
  brusearch = [r[0] for r in context.intervals_brute_force]

  eigsearch = sorted(eigsearch)
  brusearch = sorted(eigsearch)

  if len(eigsearch) != len(brusearch):
    raise AssertionError( "Different number of intervals found be methods.\n"
                          + str(context.equation) )
  
  for a, b in zip(eigsearch, brusearch):
    interval = [max(a[0], b[0]), min(a[1], b[1])]
    roots = find_roots(context.equation, [interval])
    if len(roots) == 0                                                         \
       and context.equation(interval[0]) * context.equation(interval[1]) < 0e0 :
        raise AssertionError("Could not find root in interval\n"
                             + str(context.equation));
    root = roots[0][0]
    if root < interval[0] or root > interval[1]:                  
      raise AssertionError("Root is outside of expected bound.\n"
                           + str(context.equation))

@then("roots are found within each interval for which sign changes")
def step(context):
  from dcprogs.likelihood import find_roots

  intervals = [r[0] for r in context.intervals]
  intervals = sorted(intervals)

  for interval in intervals:
    roots = find_roots(context.equation, [interval])
    if len(roots) == 0                                                         \
       and context.equation(interval[0]) * context.equation(interval[1]) < 0e0 :
      raise AssertionError("Could not find root in interval");
    root = roots[0][0]
    if root < interval[0] or root > interval[1]:                  
      raise AssertionError("Root is outside of expected bound.")
