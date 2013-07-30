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

@given('a list of {n:Integer} random {name} between {start:Float} and {end:Float}')
def step(context, n, name, start, end):
  from numpy import random, concatenate
  result = random.sample(n) * (end - start) + start
  if hasattr(context, name): result = concatenate((getattr(context, name), result))
  setattr(context, name, result)
  
@given('a parameter {name}={value:Eval}')
def step(context, name, value):
  setattr(context, name, value)

@then('instantiation did not throw')
def step(context):
  if hasattr(context, 'initialization_exception'):
    raise context.initialization_exception[1]

@then('instantiation threw {type}')
def step(context, type):
  assert hasattr(context, 'initialization_exception')
  exception = getattr(context, 'initialization_exception')
  type = eval(type)
  if not isinstance(exception[1], type): raise exception[1]

@then("no exception was thrown")
def step(context):
  if hasattr(context, 'exception'): raise context.exception[1]

