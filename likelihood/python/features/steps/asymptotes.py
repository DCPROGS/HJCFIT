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

from behave import then
from test_setup import register_type
register_type()

@then("for each root, applying H - root * I to asymptote at t={time:Float} "  \
      "yields zero (<{tolerance:Float})")
def step(context, time, tolerance):
  from dcprogs.likelihood import Asymptotes
  from numpy import identity, all, abs, dot
  for roots, equation in zip(context.roots, context.equations): 
    if roots is None: continue
    for root, multiplicity in roots: 
      H = equation.H(root) 
      HminusRoot = H - root * identity(H.shape[0])
      asymptotes = Asymptotes(equation, [(root, multiplicity)])
      apply = abs(dot(HminusRoot, asymptotes(time)) < tolerance)
      if not all(apply): 
        print("H - e * I\n{0}\n".format(HminusRoot))
        print("asymptote({0})\n{1}\n".format(time, asymptotes(time)))
        print("right apply asymptotes: {0}\n".format(apply))
        raise Exception("Incorrect asymptotic matrix for root {0} of\n{1}\n"
                        .format(root, equation))
      apply = abs(dot(asymptotes(time), HminusRoot) < tolerance)
      if not all(apply): 
        print("H - e * I\n{0}\n".format(HminusRoot))
        print("asymptote({0})\n{1}\n".format(time, asymptotes(time)))
        print("right apply asymptotes: {0}\n".format(apply))
        raise Exception("Incorrect asymptotic matrix for root {0} of\n{1}\n"
                        .format(root, equation))
