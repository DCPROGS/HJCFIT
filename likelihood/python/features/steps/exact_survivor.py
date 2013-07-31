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
      print(qmatrix)
      raise
