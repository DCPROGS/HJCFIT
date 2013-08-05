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

""" DCProgs python library. """
__docformat__ = "restructuredtext en"
__all__ = ['likelihood', 'random', 'read_idealized_bursts', 'internal_dtype']
import numpy
from . import likelihood
from .likelihood._likelihood import _dcprogs_dtype

internal_dtype = _dcprogs_dtype()
""" Type of the numpy array used internally. 

    Using this type should make some conversion from python to c++ faster or even unnecessary.

    >>> array([0, 1], dtype=internal_dtype)
"""

def read_idealized_bursts(filename, tau, tcrit): 
  """ Reads bursts data from *.scn file. 
 
      This functions is a wrapper around dc-pyps.

      :param string filename:
        Path to an *.scn file. It can also the name of a standard model.
      :param float tau:
        Resolution/Maximum length of the missed events, in seconds.
      :param float tcrit:
        Critical time, in seconds

      :returns: 
        A list of lists of intervals in milliseconds. Inner lists represent bursts.
  """
  from glob import iglob
  from os.path import exists, dirname, join, abspath, basename, splitext
  from numpy import array, all, abs
  from dcpyps.dataset import SCRecord

  if not exists(filename):
    # Check that we are not trying to read a sample data
    # First, finds all data files in directory.
    module_data_dir = join(dirname(abspath(__file__)), 'data')
    data_files = [basename(u) for u in iglob(join(module_data_dir, '*.scn'))]
    data_files = [splitext(u)[0] for u in data_files]
    # Now check if name of model compares to something in data directory.
    if filename not in data_files:
      raise IOError('Could not find file or model {0}.'.format(filename))
    # Creates full filename  and move on to reading it
    filename = join(module_data_dir, '{0}.scn'.format(filename))


  time_series = SCRecord()
  time_series.load_from_file(filename)
  time_series.impose_resolution(tau)
  time_series.get_open_shut_periods()
  time_series.get_bursts(tcrit)

  return [array(u, dtype=internal_dtype) for u in time_series.bursts]
 


