#######################
#   HJCFIT computes missed-events likelihood as described in
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
__all__ = ['network', 'find_roots', 'plot_roots',
           'missed_events_pdf', 'ideal_pdf', 'intervals_to_series', 'series_to_intervals',
           'plot_time_series', 'plot_time_intervals',
           'log_bin_edges', 'ideal_pdf_scale_factor', 'dwell_time_histogram' ]
 
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
  for i in range(qmatrix.nopen): graph.add_node(i, open=True)
  for j in range(qmatrix.nshut): graph.add_node(i+j, open=False)

  for i in range(qmatrix.matrix.shape[0]):
    for j in range(i, qmatrix.matrix.shape[1]):
      if abs(qmatrix.matrix[i,j]) > 1e-8:
        graph.add_edge(i, j)
        graph[i][j]['k+'] = qmatrix.matrix[i, j]
        graph[i][j]['k-'] = qmatrix.matrix[j, i]
  return graph

def find_roots(determinant, intervals=None, tolerance=1e-12, **kwargs):
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

       :param kwargs:
         Passed on to :py:meth:`brentq` and :py:meth:`find_root_intervals`.

       :returns: A list of items `(root, multiplicity)`.
   """
   from scipy.optimize import fminbound
   from numpy import abs, count_nonzero
   from inspect import getfullargspec
   from .likelihood import find_root_intervals, eig, brentq

   if intervals is None:
     # Create dictionary of keyword arguments
     spec = getfullargspec(find_root_intervals)
     names, defaults = spec.args, spec.defaults
     intervals_kwargs = {'tolerance': tolerance} 
     for name in names[-len(defaults):]: 
       if name in kwargs: intervals_kwargs[name] = kwargs[name]
     intervals = [u[0] for u in find_root_intervals(determinant, **intervals_kwargs)]

   spec = getfullargspec(brentq)
   names, defaults = spec.args, spec.defaults
   brentq_kwargs = {} 
   for name in names[-len(defaults):]: 
     if name in kwargs: brentq_kwargs[name] = kwargs[name]

   result = []
   for interval in intervals:
     # left, right: limit of interval.
     left, right = determinant(interval)
     if left * right < 0: root = brentq(determinant, *interval, **brentq_kwargs)[0]
     elif left < 0:
       root, value, ierr, numfunc = fminbound(lambda x: -determinant(x), left, right)
       if abs(value) > tolerance: continue
     else:
       root, value, ierr, numfunc = fminbound(determinant, left, right, full_output=True)
       if abs(value) > tolerance: continue

     H = determinant.H(root)
     if len(H) > 1:
       # Use Eigen's eigenvalue pb so that we can do long doubles.
       eigenvalues = eig(determinant.H(root))[0]
       multiplicity = int(count_nonzero(abs(eigenvalues - root) < tolerance))
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
  """ Creates pdf from knowledge of phi, g and whether open or shut.

      `result.flat[i] = v` rather than `result.itemset(i, v)`: the method was
      removed in numpy 2.0. `i` indexes `t.flat`, so the flat assignment is the
      exact equivalent and stays correct for a multi-dimensional `t`.
  """
  from numpy import dot, sum, zeros_like
  if shut: 
    def missed_events_pdf(t):
      result = zeros_like(t)
      for i, u in enumerate(t.flat):
        result.flat[i] = sum(dot(phi, g.fa(float(u))))
      return result
  else:
    def missed_events_pdf(t):
      result = zeros_like(t)
      for i, u in enumerate(t.flat):
        result.flat[i] = sum(dot(phi, g.af(float(u))))
      return result
  return missed_events_pdf

def missed_events_pdf(qmatrix, tau, nmax=2, shut=False, tcrit=None):
  """ A function to compute missed-events pdf """
  from .likelihood import MissedEventsG
  g = MissedEventsG(qmatrix, tau, nmax)

  if tcrit is not None:
    phi = g.final_CHS_vectors(tcrit) if shut else g.initial_CHS_vectors(tcrit) 
  else: 
    phi = g.final_vectors if shut else g.initial_vectors
  return _create_pdf(phi, g, shut)
    
def ideal_pdf(qmatrix, shut=False):
  """ A function to compute ideal pdf """
  from .likelihood import IdealG

  g = IdealG(qmatrix)

  phi = g.final_vectors if shut else g.initial_vectors
  return _create_pdf(phi, g, shut)
      
def exponential_pdfs(qmatrix, tau, shut=False, tcrit=None):
  """ Returns a list of function that make up the asymptotic pdf. 
  
      This is mostly to re-plot Hawkes et al (1992).
  """
  from operator import itemgetter
  from functools import partial
  from numpy import dot, sum, exp
  from .likelihood import MissedEventsG, ApproxSurvivor
    

  g = MissedEventsG(qmatrix, tau) 
  if tcrit is not None:
    phi = g.final_CHS_vectors(tcrit) if shut else g.initial_CHS_vectors(tcrit) 
  else: 
    phi = g.final_vectors if shut else g.initial_vectors

  if shut: 
    components = ApproxSurvivor(qmatrix, tau).fa_components
  else:
    components = ApproxSurvivor(qmatrix, tau).af_components
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
  x = [series[0]] + [series[i//2] for i in range(2, 2*len(series)-1)]
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

def log_bin_edges(intervals, tres, nbdec=None):
  """ Geometric bin edges for a dwell-time histogram.

      Dwell times span decades, so bins are uniform in log time rather than in
      time: each is a fixed ratio wider than the last, and the first starts at
      the resolution. The number of bins per decade follows the DCprogs
      convention, which widens the bins for smaller samples so that the counts
      stay usable.

      :param intervals: Observed dwell times, in seconds.
      :param float tres: Resolution. The histogram starts here.
      :param int nbdec:
        Bins per decade. If None, chosen from the sample size: 5 below 300
        intervals, 8 to 1000, 10 to 3000, 12 above.

      :returns: (edges, nbdec)
  """
  from numpy import asarray, log, log10, ceil, arange

  intervals = asarray(intervals, dtype=float)
  if nbdec is None:
    n = len(intervals)
    nbdec = 5 if n <= 300 else 8 if n <= 1000 else 10 if n <= 3000 else 12

  ratio = 10.0 ** (1.0 / nbdec)
  tmax = 10.0 ** ceil(log10(intervals.max()))
  nbin = int(log(tmax / tres) / log(ratio)) + 1
  return tres * ratio ** arange(nbin + 1), nbdec


def ideal_pdf_scale_factor(tres, aa, initial_vectors):
  r""" Renormalises an ideal dwell-time pdf onto the resolved intervals.

      An ideal pdf describes every sojourn; a record contains only those
      longer than the resolution. Comparing one with the other means dividing
      by the fraction that survives,

      .. math:: P(t \ge t_{res}) = \phi_A e^{Q_{AA} t_{res}} u_A

      so that the pdf integrates to one over the intervals actually observed.
      An apparent (missed-events) pdf needs no such factor: it is already a
      density over the resolved intervals.

      :param float tres: Resolution, in seconds.
      :param aa: The :math:`Q_{AA}` block, open-open or shut-shut to match.
      :param initial_vectors: The equilibrium occupancies :math:`\phi_A`.

      :returns: :math:`1 / P(t \ge t_{res})`
  """
  from numpy import asarray, dot, ones
  from scipy.linalg import expm

  aa = asarray(aa, dtype=float)
  phi = asarray(initial_vectors, dtype=float).reshape(-1)
  survival = float(dot(phi, dot(expm(aa * tres), ones(aa.shape[0]))))
  return 1.0 / survival


def dwell_time_histogram(intervals, tres, ax=None, pdf=None, ideal=None,
                         ideal_scale=1.0, tcrit=None, nbdec=None,
                         xlabel='Dwell time (s)'):
  r""" Dwell-time histogram on log time and square-root frequency.

      The display convention of Sigworth and Sine: bins uniform in log time so
      that components spanning decades are all visible, and a square-root
      ordinate so that the scatter is roughly constant down the tail, which
      makes a fitted curve easy to judge by eye.

      **Scaling.** With bins of fixed ratio :math:`dx = 10^{1/nbdec}`, the
      count expected in the bin starting at *t* from a density *f* conditional
      on :math:`t \ge t_{res}` is

      .. math:: N \int_t^{t\,dx} f(u)\,du \;pprox\; N f(t)\, t \ln dx
                = N f(t)\, t rac{\ln 10}{nbdec}

      because the bin is narrow in log space. That is what is drawn, under the
      same square root as the counts. An ideal pdf is not conditional on the
      resolution, so it carries *ideal_scale* as well -- see
      :py:func:`ideal_pdf_scale_factor`.

      :param intervals: Observed dwell times, in seconds.
      :param float tres: Resolution, in seconds.
      :param ax: Axes to draw on. A new figure is created if None.
      :param pdf:
        Apparent (missed-events) density, a callable of t. Drawn solid.
      :param ideal: Ideal density, a callable of t. Drawn dashed.
      :param float ideal_scale: Factor for *ideal*, from
        :py:func:`ideal_pdf_scale_factor`.
      :param tcrit: Critical time, or several, marked with a vertical line.
      :param int nbdec: Bins per decade; see :py:func:`log_bin_edges`.

      :returns: The axes drawn on.
  """
  from numpy import asarray, histogram, sqrt, log, log10, logspace, atleast_1d

  intervals = asarray(intervals, dtype=float)
  edges, nbdec = log_bin_edges(intervals, tres, nbdec)
  counts, edges = histogram(intervals, bins=edges)

  if ax is None:
    from matplotlib import pyplot as plt
    ax = plt.subplots(1, 1)[1]

  # Draw the bins as an outline rather than bars: with a fitted curve on top,
  # filled bars hide it.
  x = [v for pair in zip(edges[:-1], edges[1:]) for v in pair]
  y = [v for pair in zip(counts, counts) for v in pair]
  ax.semilogx(x, sqrt(y), '-k')

  scale = len(intervals) * log(10.0) / nbdec
  if pdf is not None or ideal is not None:
    t = logspace(log10(tres), log10(edges[-1]), 512)
    if pdf is not None:
      ax.plot(t, sqrt(scale * t * asarray(pdf(t), dtype=float)), '-b')
    if ideal is not None:
      ax.plot(t, sqrt(scale * ideal_scale * t * asarray(ideal(t), dtype=float)),
              '--r')

  if tcrit is not None:
    for value in atleast_1d(tcrit):
      if value > 0: ax.axvline(x=value, color='g')

  ax.set_xlabel(xlabel)
  ax.set_ylabel('sqrt(frequency)')
  return ax
