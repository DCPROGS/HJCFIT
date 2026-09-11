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

""" Maximising the likelihood over a mechanism.

    :py:class:`~HJCFIT.likelihood.Log10Likelihood` takes a Q matrix and returns
    a number. Between that and a person sits everything this module is: a
    record in the form the likelihood wants, a search over the free rate
    constants of a mechanism, and a result worth printing.

    This code was written for the reproduction of Colquhoun, Hatton & Hawkes
    (2003) and is lifted from it with the comments that took the longest to
    earn. Where one of those records a measurement -- four fits in 250, 248
    log\\ :sub:`10` units -- that is a measurement, not an estimate.

    **It does not import scalcs.** The fitter is duck-typed on a mechanism:
    anything offering ``theta()``, ``theta_unsqueeze()``, ``Rates``, ``kA`` and
    ``set_eff()`` will do, which in practice means
    :py:class:`scalcs.mechanism.Mechanism`. That keeps the likelihood itself
    free of the dependency -- ``pip install hjcfit`` still wants nothing but
    numpy and scipy, and ``hjcfit[fitting]`` is what brings a mechanism
    library.

    Two conventions the likelihood enforces, both of which bite:

    * every group must have an **odd** number of intervals -- it must start and
      end with an opening. :py:func:`trim_to_openings` does that to a whole
      record fitted as one group; burst segmentation in ``dcio`` already gives
      odd groups.
    * a Q matrix it cannot handle raises ``ArithmeticError`` or, more
      dangerously because it is silent, returns ``nan``. Both are caught in
      :py:meth:`HJCFitter.cost`.
"""
__docformat__ = "restructuredtext en"
__all__ = ['Record', 'FitResult', 'HJCFitter', 'FAILURE_COST',
           'trim_to_openings']

import time
from dataclasses import dataclass, field

import numpy as np



#: Returned by the cost function when the likelihood cannot be computed. Large,
#: finite and constant -- see :py:meth:`HJCFitter.cost` for why it is a penalty
#: rather than the perturbation HJCFIT itself used.
FAILURE_COST = 1.0e10


def trim_to_openings(intervals, amplitudes):
    """ Trim an alternating record so it starts and ends with an opening.

        The likelihood requires an odd number of intervals per group. A record
        divided into bursts already satisfies that; a whole record fitted as a
        single group generally does not.

        :param intervals: Durations, alternating open and shut.
        :param amplitudes: Matching amplitudes; zero means shut.
        :returns: A view of *intervals* beginning and ending on an opening.
    """
    amplitudes = np.asarray(amplitudes)
    start = 0 if amplitudes[0] != 0.0 else 1
    stop = len(intervals) if amplitudes[-1] != 0.0 else len(intervals) - 1
    return np.asarray(intervals)[start:stop]


@dataclass(frozen=True)
class Record:
    """ One idealised experiment, in the form the likelihood wants.

        :param conc: Agonist concentration [M].
        :param groups:
          Alternating open/shut intervals, each group starting and ending with
          an opening. One group per burst or cluster; a single group when the
          whole record is treated as coming from one channel.
        :param tres: Dead time already imposed [s].
        :param tcrit:
          What the likelihood is given as ``tcritical``. A number asks for CHS
          vectors (Colquhoun, Hawkes & Srodzinski 1996); None asks for
          equilibrium vectors (Colquhoun & Hawkes 1982). Note that this is the
          *magnitude*: elsewhere in the stack a negative ``tcrit`` is a flag
          selecting equilibrium vectors, and segmentation takes ``abs``.
        :param n_raw:
          Intervals before the dead time was imposed, if known.
        :param n_apparent:
          Intervals after it. Defaults to the number actually held.
    """

    conc: float
    groups: tuple
    tres: float
    tcrit: float = None
    n_raw: int = None
    n_apparent: int = None

    @property
    def n_intervals(self):
        """ Intervals across every group. """
        return sum(len(g) for g in self.groups)

    @property
    def n_openings(self):
        """ Openings across every group; each group starts and ends on one. """
        return sum((len(g) + 1) // 2 for g in self.groups)

    def as_lists(self):
        """ The form :py:class:`Log10Likelihood` wants: a list of lists. """
        return [list(g) for g in self.groups]

    def check(self):
        """ Raise if the groups are not what the likelihood requires.

            Called by :py:class:`HJCFitter` on construction, because an even
            group does not fail loudly -- it silently asks the likelihood for a
            product of matrices that does not alternate.

            :raises ValueError: on an empty record or an even-length group.
        """
        if not self.groups:
            raise ValueError("record has no groups")
        bad = [i for i, g in enumerate(self.groups) if len(g) % 2 == 0]
        if bad:
            raise ValueError(
                "groups must start and end with an opening, so must have an "
                "odd number of intervals; groups {0} do not (see "
                "trim_to_openings)".format(bad[:5]))
        return self

    def __str__(self):
        conc = ("{0:g} nM".format(self.conc * 1e9) if self.conc < 1e-6
                else "{0:g} uM".format(self.conc * 1e6))
        vectors = ("CHS vectors" if self.tcrit is not None
                   else "equilibrium vectors")
        raw = "" if self.n_raw is None else "{0} -> ".format(self.n_raw)
        return ("{0}: {1}{2} intervals at {3:g} us -> {4} group{5}, {6} "
                "openings, {7}".format(conc, raw, self.n_intervals,
                                       self.tres * 1e6, len(self.groups),
                                       "" if len(self.groups) == 1 else "s",
                                       self.n_openings, vectors))


@dataclass
class FitResult:
    """ What one fit produced.

        :param rates: Every rate constant by name, constrained ones included.
        :param free_names: Names of the free parameters, in order.
        :param free_values: Their fitted values, as **rates** rather than logs.
        :param log10_likelihood: :math:`\\log_{10} L_{max}`.
        :param nevals: Likelihood evaluations used.
        :param niter: Iterations the optimiser reported.
        :param nfailures:
          Evaluations that could not be computed and cost
          :py:data:`FAILURE_COST`. Reported with every fit rather than
          swallowed, so that a mechanism the likelihood struggles with is
          visible.
        :param seconds: Wall-clock time of the search.
        :param success: The optimiser's own verdict.
        :param message: The optimiser's own words.
        :param path: Best point per iteration, when asked for.
        :param evaluations: Every point evaluated, when asked for.
    """

    rates: dict
    free_names: tuple
    free_values: np.ndarray
    log10_likelihood: float
    nevals: int
    niter: int
    nfailures: int
    seconds: float
    success: bool
    message: str
    path: list = field(default_factory=list, repr=False)
    evaluations: list = field(default_factory=list, repr=False)

    def __str__(self):
        head = ("log10(L) = {0:.4f}   {1} evaluations in {2:.1f} s"
                .format(self.log10_likelihood, self.nevals, self.seconds))
        if self.nfailures:
            head += "   {0} could not be computed".format(self.nfailures)
        if not self.success:
            head += "   [{0}]".format(self.message)
        lines = [head]
        for name, value in zip(self.free_names, self.free_values):
            lines.append("  {0:9s} {1:12.4g}".format(name, value))
        return "\n".join(lines)


class HJCFitter:
    """ Maximise the HJC likelihood of one or more records over a mechanism.

        :param mec:
          A mechanism already carrying its constraints and its initial guess.
          Duck-typed: it needs ``theta()``, ``theta_unsqueeze()``, ``Rates``,
          ``kA`` and ``set_eff()``. ``scalcs.mechanism.Mechanism`` provides
          them.
        :param records:
          A sequence of :py:class:`Record`, fitted **simultaneously**. The
          likelihood takes one Q matrix per call, so each record is evaluated
          at its own concentration and the logarithms added.
        :param bool log_params:
          Search the logarithms of the rate constants rather than the rates.
          This is HJCFIT's own default and three to four times faster
          (Colquhoun, Hatton & Hawkes 2003, p. 702). It also cannot produce a
          negative rate. The fits of that paper's Figures 2-5 and 12-13 were
          made over the rates themselves, which is what the resetting below is
          for.
        :param bool store_path:
          Keep the best vertex at each iteration. Cheap.
        :param bool store_evaluations:
          Keep every point evaluated. Not cheap across many fits.
    """

    def __init__(self, mec, records, log_params=True,
                 store_path=False, store_evaluations=False):
        self.mec = mec
        self.records = [r.check() for r in records]
        if not self.records:
            raise ValueError("no records to fit")
        self.log_params = log_params
        self.store_path = store_path
        self.store_evaluations = store_evaluations

        # A rate constant that leaves its allowed range is reset before the
        # likelihood sees it, rather than the search being stopped from going
        # there. That is what HJCFIT did (p. 702): an upper limit "to prevent
        # physically unrealistic values", and a floor because "if a value of a
        # rate constant should go negative during the fitting process, it can
        # be reset to a value near zero". Searching the rates themselves
        # without it produced four fits in 250 with negative rate constants.
        #
        # HJCFIT.likelihood.optimization.reset_out_of_range is the same idea
        # for a caller driving the simplex directly. Here it is done through
        # the mechanism, which knows each rate's own limits.
        free = [r for r in mec.Rates if r.is_free]
        self._lower = np.array([r.limits[0][0] if r.limits else 1e-12
                                for r in free])
        self._upper = np.array([r.limits[0][1] if r.limits else np.inf
                                for r in free])

        # Imported here rather than at module scope. Taking it at module
        # scope and exposing this module from the package __init__ made
        # HJCFIT.likelihood import itself while partially initialised; every
        # CI job failed on it, and a comment claiming there was no
        # circularity was wrong. This module is imported on demand --
        # `from HJCFIT.likelihood.fitting import HJCFitter` -- and reaches the
        # likelihood only once a fitter is actually built.
        from .likelihood import Log10Likelihood

        self.likelihoods = [
            Log10Likelihood(r.as_lists(), nopen=mec.kA, tau=r.tres,
                            tcritical=r.tcrit)
            for r in self.records
        ]
        self.nevals = 0
        self.nfailures = 0
        self.path = []
        self.evaluations = []

    # -- the likelihood ----------------------------------------------------

    def _apply(self, x):
        """ Put the parameters onto the mechanism, resetting out-of-range. """
        rates = np.exp(x) if self.log_params else np.asarray(x, dtype=float)
        self.mec.theta_unsqueeze(np.clip(rates, self._lower, self._upper))

    def log10_likelihood(self, x=None):
        """ Summed :math:`\\log_{10}` likelihood of every record.

            The likelihood takes one Q matrix per call, so each record is
            evaluated at its own concentration and the logarithms added.

            Raises rather than returning a sentinel; :py:meth:`cost` catches.

            :param x: Parameters, in the space being searched. None uses the
                mechanism as it stands.
            :raises ArithmeticError: if any record's likelihood is not finite.
        """
        if x is not None:
            self._apply(x)
        total = 0.0
        for lik, record in zip(self.likelihoods, self.records):
            self.mec.set_eff('c', record.conc)
            value = lik(self.mec.Q)
            if not np.isfinite(value):
                raise ArithmeticError("likelihood is {0}".format(value))
            total += value
        return total

    def ln_likelihood(self, x=None):
        """ The same thing in natural logarithms.

            Anything treating the log likelihood as a statistical quantity --
            a Hessian, and so the covariance matrix, the standard deviations
            and the likelihood intervals -- needs natural logarithms. Getting
            it wrong is not obvious in the result: it inflates every standard
            deviation by exactly :math:`\\sqrt{\\ln 10} = 1.517`, which looks
            like a badly behaved fit rather than a units error.
        """
        return self.log10_likelihood(x) * np.log(10.0)

    def cost(self, x):
        """ Negative summed :math:`\\log_{10}` likelihood, for a minimiser.

            A Q matrix the likelihood cannot handle raises ``ArithmeticError``,
            or else returns ``nan`` silently. Both become
            :py:data:`FAILURE_COST`.

            This is **not** what HJCFIT itself did. It kept the best parameters
            so far and, on a failure, replaced the current ones with those plus
            a bounded random perturbation (p. 702) -- which its own simplex
            allows, because it owns its vertices. SciPy's Nelder-Mead does not
            expose them, so a large constant penalty is used instead. The
            failure count is reported with every fit, so it is visible if this
            ever matters; the paper had two numerical failures in nearly 50 000
            fits.
        """
        self.nevals += 1
        try:
            value = -self.log10_likelihood(x)
        except (ArithmeticError, ValueError, RuntimeError, FloatingPointError):
            self.nfailures += 1
            value = FAILURE_COST
        if self.store_evaluations:
            self.evaluations.append((np.array(x, dtype=float), value))
        return value

    # -- the search --------------------------------------------------------

    def fit(self, x0=None, search="simplex", restarts=0,
            xatol=1e-4, fatol=1e-4, maxfev=20000, maxiter=20000,
            options=None):
        """ Run the search and return a :py:class:`FitResult`.

            :param x0:
              Starting point, in the space being searched. Defaults to the
              mechanism's current free parameters.
            :param search:
              ``"simplex"`` is
              :py:func:`~HJCFIT.likelihood.optimization.simplex_hjc`, HJCFIT's
              own -- the search that produced every published result, and the
              default here for that reason. ``"scipy"`` is SciPy's
              Nelder-Mead with the restart loop below; reach for it when the
              starting point is poor, because a regular simplex needs every
              one of its vertices to be evaluable and a random Q matrix
              usually cannot offer that.
            :param options:
              Passed to whichever search is chosen. The simplex's own defaults
              are HJCFIT's, so this is for deliberate departures.
            :param restarts:
              ``"scipy"`` only: extra runs from the previous solution,
              stopping early when one gains less than *fatol*. The simplex has
              its own restart rule, capped at ``nresmax``.
            :param xatol: ``"scipy"`` only.
            :param fatol: ``"scipy"`` only.
            :param maxfev: Evaluation budget, both searches.
            :param maxiter: ``"scipy"`` only.
        """
        if search not in ("simplex", "scipy"):
            raise ValueError("search must be 'simplex' or 'scipy', not "
                             + repr(search))
        if x0 is None:
            theta = np.asarray(self.mec.theta(), dtype=float)
            x0 = np.log(theta) if self.log_params else theta
        x0 = np.asarray(x0, dtype=float)

        self.nevals = self.nfailures = 0
        self.path = []
        self.evaluations = []

        callback = None
        if self.store_path:
            def callback(xk, *_):
                self.path.append(np.array(xk, dtype=float))

        started = time.perf_counter()
        if search == "simplex":
            from .optimization import simplex_hjc
            result = simplex_hjc(self.cost, x0, logfit=self.log_params,
                                 maxfev=maxfev, callback=callback,
                                 **(options or {}))
            niter = int(result.nit)
        else:
            from scipy.optimize import minimize
            opts = dict(maxiter=maxiter, maxfev=maxfev,
                        xatol=xatol, fatol=fatol)
            opts.update(options or {})
            result = minimize(self.cost, x0, method="Nelder-Mead",
                              options=opts, callback=callback)
            for _ in range(restarts):
                previous = result.fun
                result = minimize(self.cost, result.x, method="Nelder-Mead",
                                  options=opts, callback=callback)
                if previous - result.fun < fatol:
                    break
            niter = int(result.nit)
        seconds = time.perf_counter() - started

        # Read the rates back off the mechanism rather than off the optimiser,
        # so that a reset applied in _apply is reflected in what is reported
        # and not only in what was fitted.
        self._apply(result.x)
        values = np.array([r.unit_rate() for r in self.mec.Rates if r.is_free])
        return FitResult(
            rates={r.name: float(np.ravel(r.rateconstants)[0])
                   for r in self.mec.Rates},
            free_names=tuple(self.mec.get_free_parameter_names()),
            free_values=values,
            log10_likelihood=-float(result.fun),
            nevals=self.nevals,
            niter=niter,
            nfailures=self.nfailures,
            seconds=seconds,
            success=bool(result.success),
            message=str(result.message),
            path=list(self.path),
            evaluations=list(self.evaluations),
        )
