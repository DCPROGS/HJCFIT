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

""" Some pure python methods used to access/complement the c++ bindings. """
__docformat__ = "restructuredtext en"
__all__ = ['network', 'create_approx_survivor', 'create_missed_eventsG', 'find_roots', 'plot_roots',
           'missed_events_pdf', 'ideal_pdf', 'intervals_to_series', 'series_to_intervals',
           'plot_time_series', 'plot_time_intervals' ]
 
def network(qmatrix): 
  """ Creates networkx graph object from a :class:`QMatrix` object.
  
      Vertices have an "open" attribute to indicate whether they are open or shut. Edges have a
      "k+" and "k-" attribute containing the transition rates for the node with smaller index to
      the node with larger index, and vice-versa. 

      :param qmatrix: 
        A :class:`QMatrix` instance for which to construct a graph
      :return: A `networkx.Graph` object which contains all the information held by the qmatrix
               object, but in a graph format.
  """
  from networkx import Graph

  graph = Graph()
  for i in xrange(qmatrix.nopen): graph.add_node(i, open=True)
  for j in xrange(qmatrix.nshut): graph.add_node(i+j, open=False)

  for i in xrange(qmatrix.matrix.shape[0]):
    for j in xrange(i, qmatrix.matrix.shape[1]):
      if abs(qmatrix.matrix[i,j]) > 1e-8:
        graph.add_edge(i, j)
        graph[i][j]['k+'] = qmatrix.matrix[i, j]
        graph[i][j]['k-'] = qmatrix.matrix[j, i]
  return graph

def _input_to_cpp_wrappers(qmatrix, tau): 
  """ Creates the input we need to create MissedEvents and ApproxSurvivor. """
  from .likelihood import DeterminantEq
  determinant_af = DeterminantEq(qmatrix, tau) 
  determinant_fa = DeterminantEq(qmatrix.transpose(), tau) 
  roots_af = find_roots(determinant_af)
  roots_fa = find_roots(determinant_fa)
  return determinant_af, roots_af, determinant_fa, roots_fa

def create_approx_survivor(qmatrix, tau):
  """ Creates a ApproxSurvivor function from knowledge of rate matrix. """
  from .likelihood import ApproxSurvivor
  return ApproxSurvivor(*_input_to_cpp_wrappers(qmatrix, tau))

def create_missed_eventsG(qmatrix, tau, nmax=2):
  """ Creates a MissedEvents function from knowledge of rate matrix. """
  from .likelihood import MissedEventsG
  args = list(_input_to_cpp_wrappers(qmatrix, tau)) + [nmax]
  return MissedEventsG(*args)

def find_roots(determinant, intervals=None, tolerance=1e-8):
   """ Computes roots for each interval. 
   
       :param determinant: 
         A function or functor of a single variable.
       :param intervals:
         A list of items `[(a0, b0), ..., (a1, b1)]`, where `(a, b)` is the interval over which to
         look for roots. 

         If this object is None (default), then uses :py:meth:`find_root_intervals` to figure out
         the  intervals.
       :param tolerance:
         Tolerance criteria. Used to determine multiplicity.
       :returns: A list of items `(root, multiplicity)`.
   """
   from scipy.optimize import brentq, fminbound
   from numpy import abs, count_nonzero
   from .likelihood import find_root_intervals, eig

   if intervals is None:
     intervals = [u[0] for u in find_root_intervals(determinant)]

   result = []
   for interval in intervals:
     # left, right: limit of interval.
     left, right = determinant(interval)
     if left * right < 0: root = brentq(determinant, *interval)
     elif left < 0:
       root, value, ierr, numfunc = fminbound(lambda x: -determinant(x), left, right)
       if abs(value) > tolerance: continue
     else:
       root, value, ierr, numfunc = fminbound(determinant, left, right, full_output=True)
       if abs(value) > tolerance: continue

     H = determinant.H(root)
     if len(H) > 1:
       # Use Eigen's eigenvalue pb so that we can do 128 bit reals. 
       eigenvalues = eig(determinant.H(root))[0]
       multiplicity = count_nonzero(abs(eigenvalues - root) < tolerance)
     else: multiplicity = 1
     if left * right < 0:
       if multiplicity == 0 or multiplicity % 2 != 1: multiplicity = 1
     else:
       if multiplicity == 0 or multiplicity % 2 != 0: multiplicity = 2
     result.append((root, multiplicity));
   return result;

def plot_roots(determinant, intervals=None, figure=None, main=None, lines=None, size=1000,
               tolerance=1e-8, ax=None): 
   """ Computes and plots roots. 

       :param determinant: 
         A function or functor of a single variable.
       :param intervals:
         A list of items `[(a0, b0), ..., (a1, b1)]`, where `(a, b)` is the interval over which to
         look for roots. 

         If this object is None (default), then uses :py:meth:`find_root_intervals` to figure out
         the  intervals.
       :param main:
         A dictionary of values with which to plot the determinant.
       :param lines:
         A dictionary of values with which to plot the roots.
         
       :returns: A figure
   """
   from matplotlib import pyplot as plt
   from numpy import arange, min, max, array
   from .likelihood import find_root_intervals

   if intervals is None:
     intervals = [u[0] for u in find_root_intervals(determinant)]
   intervals = array(intervals)
   if main is None: main = {}
   if lines is None: lines = {}

   roots = find_roots(determinant, intervals, tolerance)
   mini = min([u[0] for u in roots])
   maxi = max([u[0] for u in roots])
   diff = maxi - mini
   maxi += diff  * 0.05
   mini -= diff  * 0.05
   
   x = arange(mini, maxi+(maxi-mini)/float(size)*0.5, (maxi-mini)/float(size))
   y = determinant(x)

   if ax is None: 
     figure = plt.figure()
     ax = figure.add_subplot(111)
   ax.plot(x, y, **main)
   ax.set_xlim((mini, maxi))

   ymin, ymax = min(y), max(y)
   ymin = ymin - (ymax - ymin) * 0.05
   ymax = ymax + (ymax - ymin) * 0.05
   ax.set_ylim((ymin, ymax))

   for root in roots:
     ax.plot([root[0], root[0]], [ymin, 0], **lines)
   ax.plot([x[0], x[-1]], [0, 0])

   return figure

def _create_pdf(phi, g, shut):
  """ Creates pdf from knowledge of phi, g and whether open or shut."""
  from numpy import dot, sum, zeros_like
  if shut: 
    def missed_events_pdf(t):
      result = zeros_like(t)
      for i, u in enumerate(t.flat):
        result.itemset(i, sum(dot(phi, g.fa(float(u)))))
      return result
  else:
    def missed_events_pdf(t):
      result = zeros_like(t)
      for i, u in enumerate(t.flat):
        result.itemset(i, sum(dot(phi, g.af(float(u)))))
      return result
  return missed_events_pdf

def missed_events_pdf(qmatrix, tau, nmax=2, shut=False, tcrit=None):
  """ A function to compute missed-events pdf """
  g = create_missed_eventsG(qmatrix, tau, nmax)

  if tcrit is not None:
    phi = g.final_CHS_occupancies(tcrit) if shut else g.initial_CHS_occupancies(tcrit) 
  else: 
    phi = g.final_occupancies if shut else g.initial_occupancies
  return _create_pdf(phi, g, shut)
    
def ideal_pdf(qmatrix, shut=False):
  """ A function to compute ideal pdf """
  from .likelihood import IdealG

  g = IdealG(qmatrix)

  phi = g.final_occupancies if shut else g.initial_occupancies
  return _create_pdf(phi, g, shut)
      
def exponential_pdfs(qmatrix, tau, shut=False, tcrit=None):
  """ Returns a list of function that make up the asymptotic pdf. 
  
      This is mostly to re-plot Hawkes et al (1992).
  """
  from operator import itemgetter
  from functools import partial
  from numpy import dot, sum, exp
  from dcprogs.likelihood import create_missed_eventsG, create_approx_survivor
    

  g = create_missed_eventsG(qmatrix, tau) 
  if tcrit is not None:
    phi = g.final_CHS_occupancies(tcrit) if shut else g.initial_CHS_occupancies(tcrit) 
  else: 
    phi = g.final_occupancies if shut else g.initial_occupancies

  if shut: 
    components = create_approx_survivor(qmatrix, tau).fa_components
  else:
    components = create_approx_survivor(qmatrix, tau).af_components
  components = sorted(components, key=itemgetter(1))

  def function(coef, root, t): return coef * exp(root * (t-tau))
        
  results = []
  if shut:
    for matrix, root in components:
      coef = sum(dot(phi, dot(matrix, g.fa_factor)))
      results.append(partial(function, coef, root))
  else:
    for matrix, root in components:
      coef = sum(dot(phi, dot(matrix, g.af_factor)))
      results.append(partial(function, coef, root))
  return results

def intervals_to_series(intervals, start=0):
  """ Converts time intervals to time series. """
  from numpy import zeros
  from .. import internal_dtype
  result = zeros(len(intervals)+1, dtype=internal_dtype)
  result[0] = start
  for i, z in enumerate(intervals): result[i+1] = result[i] + z
  return result

def series_to_intervals(series, start=0):
  """ Converts time intervals to time series. """
  return series[1:] - series[:-1]

def plot_time_series(series, ax=None, **kwargs):
  """ Plots time series """
  from numpy import array, min, max, arange
  from matplotlib.pylab import figure
  x = [series[0]] + [series[i/2] for i in range(2, 2*len(series)-1)]
  y = (arange(0, len(x)) % 4 >= 2).astype('int')
  if ax is None:
    fig = figure()
    ax = fig.add_subplot(111)
  ax.plot(x, 0.1 + array(y), **kwargs)
  ax.set_xlabel("time")
  ax.set_ylim((0, 1.2))
  ax.set_xlim((min(x)-0.1, max(x)+0.1))
  ax.yaxis.set_visible(False)

def plot_time_intervals(series, start=0, ax = None):
  """ Plots time intervals """
  return plot_time_series(series_to_intervals(series, start), ax=ax)

