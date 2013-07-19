from behave import given, when, then
from test_setup import register_type
register_type()

@given('the {doopen}-states determinantal equation \"{qmatrix:QMatrix}\" with tau={tau:Float}')
def step(context, doopen, qmatrix, tau):
  from dcprogs.likelihood import DeterminantEq
  context.equation = DeterminantEq(qmatrix, tau)
  if doopen != "open": context.equation = context.equation.transpose()
  print context.equation

@given('a list of {n:Integer} random determinant equations')
def step(context, n):
  from dcprogs.likelihood import DeterminantEq
  from dcprogs.random import qmatrix as random_qmatrix
  context.matrices = []
  context.equations = []
  while len(context.equations) < n:
    try:
      matrix = random_qmatrix()
      equation =  DeterminantEq(matrix, 1e-4)
    except: continue
    else:
      context.matrices.append(matrix)
      context.equations.append(equation)

@given('allowing for {n:Integer}% failure in tests below')
def step(context, n):
  context.allowance = float(n) / 100.0

@when("the root intervals are computed")
def step(context):
  from dcprogs.likelihood import find_root_intervals
  if hasattr(context, "exception"): raise context.exception
  context.intervals = find_root_intervals(context.equation)

@when("the roots are computed for each")
def step(context):
  from dcprogs.likelihood import find_roots
  context.roots = []
  for equation in context.equations:
    try: roots = find_roots(equation)
    except: roots = None
    context.roots.append(roots)

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

@then('roots are roots indeed, to within tolerance={tolerance:Float} * variation of det W')
def step(context, tolerance): 
  from numpy import max
  isOK = 0
  for (roots, matrix, equation) in zip(context.roots, context.matrices, context.equations): 
    if roots is None: continue # eigenvalue problem. No to be tested here.
    # tries to figure out reasonnable tolerance. 
    # the function changes very rapidly, so that makes it tricky
    variations = [(equation(root[0]+tolerance), equation(root[0]-tolerance)) for root in roots]
    variations = [max([abs(u[0]), abs(u[1]), abs(u[0] - u[1])]) for u in variations]
    false_roots = [(root[0], abs(equation(root[0]))) for variation, root in zip(variations, roots)
                   if abs(equation(root[0])) > tolerance * variation]
    if len(false_roots) > 0:
      isOK += 1 
      print("(s, error)={0} are not roots of {1}.\n"
            .format(false_roots, matrix)) 
  if isOK * getattr(context, 'allowance', 1) >= 1:
    raise AssertionError("Found {0}/{1} systems with incorrect roots."\
                         .format(isOK, len(context.equations)))

@then('the multiplicity of the real roots add up to the number of open states')
def step(context): 
  isOK = 0
  for (roots, matrix) in zip(context.roots, context.matrices): 
    if roots is None: continue # eigenvalue problem. No to be tested here.
    nroots = sum([r[1] for r in roots])
    if nroots != matrix.aa.shape[0]:
      isOK += 1 
      print("Incorrect number of roots ({0}/{1}) for {2}.\n"
            "Roots are {3}."
            .format(len(roots), matrix.aa.shape[0], matrix, roots)) 
  if isOK * getattr(context, 'allowance', 1) >= 1:
    raise AssertionError("Found {0}/{1} systems with incorrect number of roots."\
                         .format(isOK, len(context.equations)))
