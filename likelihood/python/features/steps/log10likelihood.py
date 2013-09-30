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

@given('a list of bursts {bursts:Eval}')
def step(context, bursts):
  context.bursts = bursts
@given('a dictionary of {name}')
def step(context, name):
  setattr(context, name, eval(context.text))

@when('The log10-likelihood of the classic matrix is computed')
def setp(context):
  from test_setup import Matrix
  from dcprogs.likelihood import Log10Likelihood
  matrix = Matrix("classic")
  context.result = Log10Likelihood(bursts=context.bursts, **context.parameters)(matrix)

@then('the result is defined and finite')
def step(context):
  from numpy import isnan, isinf
  assert not isnan(context.result)
  assert not isinf(context.result)
