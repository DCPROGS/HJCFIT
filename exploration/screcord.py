"""A single-channel record for the example notebooks, on dcio and scalcs.

This replaces `dcpyps.dataset.SCRecord`, which the examples used to read an
*.scn file, impose a resolution and cut the result into bursts. dcpyps is
deprecated and has not been installable for years; everything it did here is
available from `dcio` (reading and dead-time resolution) and `scalcs`
(burst segmentation), which is the same pairing `HJCFIT.read_idealized_bursts`
uses.

Neither package is a dependency of HJCFIT. They are imported here, in an
example helper, and nowhere in the library.
"""

import numpy as np

from dcio.analysis import from_scn
from dcio.formats import scn
from scalcs.scsim import extract_burst_intervals


class _Burst:
    """One burst: the interval sequence, plus the few queries the examples make.

    A burst starts and ends on an opening, so it holds an odd number of
    intervals and (n + 1) // 2 of them are openings.
    """

    def __init__(self, intervals):
        self.intervals = np.asarray(intervals)

    def get_openings_number(self):
        return (len(self.intervals) + 1) // 2

    def get_length(self):
        return float(self.intervals.sum())

    def __len__(self):
        return len(self.intervals)

    def __iter__(self):
        return iter(self.intervals)

    def __array__(self, dtype=None, copy=None):
        return self.intervals if dtype is None else self.intervals.astype(dtype)

    def __repr__(self):
        return ('Burst({0} intervals, {1} openings, {2:.4g} s)'
                .format(len(self.intervals), self.get_openings_number(),
                        self.get_length()))


class _Bursts:
    """The interval sequences of each burst.

    Named for `SCRecord.bursts.intervals()`, which is how the notebooks reach
    them.
    """

    def __init__(self, intervals):
        self._intervals = intervals

    def intervals(self):
        return self._intervals

    def count(self):
        """Number of bursts."""
        return len(self._intervals)

    def all(self):
        """Every burst, as :class:`_Burst` objects.

        `intervals()` returns the raw arrays, which is what the likelihood
        wants; `all()` returns objects that can also be asked how many
        openings they contain.
        """
        return [_Burst(b) for b in self._intervals]

    def __len__(self):
        return len(self._intervals)

    def __iter__(self):
        return iter(self._intervals)

    def __getitem__(self, i):
        return self._intervals[i]

    def __repr__(self):
        if not self._intervals:
            return 'Bursts(0)'
        lengths = [len(b) for b in self._intervals]
        return ('Bursts({0}, intervals per burst {1}-{2})'
                .format(len(self._intervals), min(lengths), max(lengths)))


class SCRecord:
    """Idealised single-channel record, resolved and cut into bursts.

    :param filenames:
        One *.scn path, or several to be pooled at the same concentration.
    :param conc:
        Agonist concentration in M. Carried for the caller's convenience; it
        does not affect the record.
    :param tres:
        Resolution, in seconds.
    :param tcrit:
        Critical shut time separating bursts, in seconds. Only its magnitude
        is used here: the sign is a flag to
        :py:class:`~HJCFIT.likelihood.Log10Likelihood`, where a negative value
        selects equilibrium vectors over CHS vectors. Pass the same value to
        both.
    """

    def __init__(self, filenames, conc=None, tres=0.0, tcrit=None,
                 record_type='recorded'):
        if isinstance(filenames, str):
            filenames = [filenames]

        self.filenames = list(filenames)
        self.conc = conc
        self.tres = tres
        self.tcrit = tcrit
        self.record_type = record_type

        self.records, bursts = [], []
        for filename in self.filenames:
            # Bursts are extracted per file and pooled, rather than pooling
            # the raw records: concatenating two patches would invent an
            # interval boundary at the join.
            record = from_scn(scn.read(filename), tres=tres)
            self.records.append(record)
            if tcrit is not None:
                # flags: an unusable interval has no defined length, ends
                # the burst before it, and is never compared with tcrit
                bursts.extend(extract_burst_intervals(record.resolved_intervals,
                                                      record.resolved_amplitudes,
                                                      abs(tcrit),
                                                      flags=record.resolved_flags))
        self.bursts = _Bursts(bursts)

    @property
    def opint(self):
        """Open period durations, in seconds, pooled over the files."""
        return np.concatenate([r.open_periods for r in self.records])

    @property
    def shint(self):
        """Shut period durations, in seconds, pooled over the files."""
        return np.concatenate([r.shut_periods for r in self.records])

    def printout(self, output=None):
        text = str(self)
        if output is None:
            print(text)
        else:
            output.write(text + '\n')
        return text

    def __str__(self):
        lines = ['Data loaded from: ' + ', '.join(self.filenames)]
        if self.conc is not None:
            lines.append('Concentration of agonist = {0:.3f} microMolar'
                         .format(self.conc * 1e6))
        lines.append('Resolution for HJC calculations = {0:.1f} microseconds'
                     .format(self.tres * 1e6))
        if self.tcrit is not None:
            lines.append('Critical gap length to define end of group (tcrit) '
                         '= {0:.3f} milliseconds'.format(abs(self.tcrit) * 1e3))
        lines.append('Number of resolved intervals = {0:d}'
                     .format(sum(r.n_resolved for r in self.records)))
        op, sh = self.opint, self.shint
        lines.append('Number of resolved periods = {0:d}'.format(len(op) + len(sh)))
        lines.append('Number of open periods = {0:d}'.format(len(op)))
        if len(op):
            lines.append('Mean and SD of open periods = {0:.9f} +/- {1:.9f} ms'
                         .format(op.mean() * 1e3, op.std() * 1e3))
            lines.append('Range of open periods from {0:.9f} ms to {1:.9f} ms'
                         .format(op.min() * 1e3, op.max() * 1e3))
        if len(sh):
            lines.append('Mean and SD of shut periods = {0:.9f} +/- {1:.9f} ms'
                         .format(sh.mean() * 1e3, sh.std() * 1e3))
            lines.append('Range of shut periods from {0:.9f} ms to {1:.9f} ms'
                         .format(sh.min() * 1e3, sh.max() * 1e3))
        lines.append('Number of bursts = {0:d}'.format(len(self.bursts)))
        return '\n'.join(lines)

    __repr__ = __str__
