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

@given('the {doopen}-states determinant equation \"{qmatrix:QMatrix}\" with tau={tau:Float}')
def step(context, doopen, qmatrix, tau):
  from dcprogs.likelihood import DeterminantEq
  context.equation = DeterminantEq(qmatrix, tau)
  if doopen != "open": context.equation = context.equation.transpose()
  print(context.equation)

@given('a list of {n:Integer} random determinant equations')
def step(context, n):
  from dcprogs.likelihood import DeterminantEq
  from dcprogs.likelihood.random import qmatrix as random_qmatrix
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

@given('a determinantal equation we know has large roots')
def step(context):
  from dcprogs.likelihood import QMatrix, DeterminantEq
  qmatrix = QMatrix([[ -1.89907444e+02,   0.00000000e+00,   2.17917781e+01,
                        1.68115666e+02,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00],
                     [  0.00000000e+00,  -5.59197168e+02,   0.00000000e+00,
                        0.00000000e+00,   2.28813670e+01,   5.36315801e+02,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00],
                     [  1.70657017e+04,   0.00000000e+00,  -1.70657017e+04,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00],
                     [  6.93510518e+04,   0.00000000e+00,   0.00000000e+00,
                       -7.70252433e+04,   0.00000000e+00,   0.00000000e+00,
                        7.67419154e+03,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00],
                     [  0.00000000e+00,   1.62124166e+04,   0.00000000e+00,
                        0.00000000e+00,  -1.62124166e+04,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00],
                     [  0.00000000e+00,   7.15538758e+04,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,  -1.06690747e+06,
                        9.95292553e+05,   0.00000000e+00,   6.10438859e+01,
                        0.00000000e+00,   0.00000000e+00],
                     [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        1.45242879e+04,   0.00000000e+00,   3.69900307e+05,
                       -3.87202655e+05,   2.77806009e+03,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00],
                     [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        2.01908451e+03,  -2.28998927e+03,   2.70904758e+02,
                        0.00000000e+00,   0.00000000e+00],
                     [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   2.87990000e+04,
                        0.00000000e+00,   1.09766807e+02,  -2.90893700e+04,
                        1.80603172e+02,   0.00000000e+00],
                     [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   2.19533614e+02,
                       -3.09835200e+02,   9.03015859e+01],
                     [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                        3.29300422e+02,  -3.29300422e+02]], 3)
  context.determinant = DeterminantEq(qmatrix, 4e-5).transpose()

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

@when("searching for intervals for each root")
def step(context):
  from dcprogs.likelihood import find_root_intervals
  try: find_root_intervals(context.determinant, **context.parameters)
  except: context.noerror = False
  else: context.noerror = True

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
  if isOK / float(len(context.roots)) > getattr(context, 'allowance', 1):
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

@then("the call does not stack-overflow")
def step(context):
  assert context.noerror
