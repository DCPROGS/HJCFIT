""" DCProgs python library. """
__docformat__ = "restructuredtext en"
__all__ = ['likelihood', 'random', 'read_idealized_bursts']
import numpy
import likelihood

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
  from numpy import array
  from dcpyps.dataset import TimeSeries

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


  time_series = TimeSeries()
  time_series.load_from_file(filename)
  time_series.impose_resolution(tau)
  time_series.get_open_shut_periods()
  time_series.get_bursts(tcrit)

  return [array(u, dtype='float64') * 1e-3 for u in time_series.bursts.itervalues()]
