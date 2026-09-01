########################
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

""" HJCFIT python library. """
__docformat__ = "restructuredtext en"
__all__ = ['likelihood', 'random', 'read_idealized_bursts', 'internal_dtype']
import numpy
from . import likelihood
from .likelihood._likelihood import _HJCFIT_dtype

internal_dtype = _HJCFIT_dtype()
""" Type of the numpy array used internally.

    Using this type should make some conversion from python to c++ faster or
    even unnecessary.

    >>> array([0, 1], dtype=internal_dtype)
"""


def read_idealized_bursts(filename, tau, tcrit):
  """ Reads bursts data from an *.scn file.

      Idealised single-channel records are read and dead-time corrected with
      :mod:`dcio`, then segmented into bursts at shut intervals of at least
      *tcrit* with :func:`scalcs.scsim.extract_burst_intervals`. Neither
      package is a dependency of HJCFIT -- they are imported here and only
      here, so the rest of the library works without them.

      This used to be a wrapper around dc-pyps, which is deprecated and has
      not been installable for years. dcio supplies what dcpyps.dataset and
      dcpyps.dcio did; scalcs supplies what dcpyps.mechanism did.

      :param string filename:
        Path to an *.scn file. It can also be the name of one of the sample
        records shipped with HJCFIT ("CH82", "CO", "CCO"), with or without
        the .scn extension.
      :param float tau:
        Resolution/Maximum length of the missed events, in seconds.
      :param float tcrit:
        Critical time, in seconds. Only its magnitude is used here. The sign
        is a flag to :py:class:`~HJCFIT.likelihood.Log10Likelihood`, where a
        negative value selects equilibrium vectors (Colquhoun & Hawkes 1982)
        instead of CHS vectors (Colquhoun, Hawkes & Srodzinski 1996); pass
        the same value to both.

      :returns:
        A list of arrays of intervals **in seconds**. Each array is one burst,
        alternating open and shut and beginning and ending with an opening, so
        every array has odd length. This is what
        :py:class:`~HJCFIT.likelihood.Log10Likelihood` expects, and it takes
        *tau* and *tcrit* in the same units.
  """
  from os.path import exists, dirname, join, abspath, splitext
  from numpy import array

  try:
    from dcio.analysis import from_scn
    from dcio.formats import scn
    from scalcs.scsim import extract_burst_intervals
  except ImportError as e:
    raise ImportError(
        "read_idealized_bursts needs dcio and scalcs, which are not "
        "dependencies of HJCFIT and are not on PyPI. Install them from "
        "https://github.com/remislp/dcio and "
        "https://github.com/remislp/SCALCS. (Original error: {0})".format(e))

  if not exists(filename):
    # Not a path: check whether it names one of the sample records, with or
    # without the extension.
    stem = splitext(filename)[0]
    sample = join(dirname(abspath(__file__)), 'data', '{0}.scn'.format(stem))
    if not exists(sample):
      raise IOError('Could not find file or sample record {0}.'.format(filename))
    filename = sample

  record = from_scn(scn.read(filename), tres=tau)
  # The sign of tcrit is a flag to Log10Likelihood -- negative selects
  # equilibrium vectors over CHS vectors -- while its magnitude is the
  # critical time. Segmentation uses the magnitude.
  # Flags matter: time-course fitting leaves an interval it could not measure
  # with no defined length, flagged unusable. It ends the burst before it, and
  # its nominal duration must never be compared with tcrit -- in one of the
  # Burzomato records that duration is 46 us where tcrit is a second.
  bursts = extract_burst_intervals(record.resolved_intervals,
                                   record.resolved_amplitudes, abs(tcrit),
                                   flags=record.resolved_flags)

  return [array(u, dtype=internal_dtype) for u in bursts]
