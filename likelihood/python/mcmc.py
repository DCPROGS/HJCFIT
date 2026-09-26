########################
#   HJCFIT computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2026  University College London
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

""" Bayesian inference over a mechanism.

    :py:class:`~HJCFIT.likelihood.fitting.HJCFitter` looks for the rate
    constants that maximise the HJC likelihood. This module treats them as
    random variables instead, and describes their posterior distribution: the
    likelihood of the records multiplied by a prior. It follows Epstein,
    Calderhead, Girolami & Sivilotti (2016) *Biophys J* 111:333-348, who
    sampled that posterior by Markov chain Monte Carlo with the same exact
    missed-events likelihood.

    Like the fitter, **it does not import scalcs.** A mechanism is duck-typed:
    ``theta()``, ``theta_unsqueeze()``, ``Rates``, ``kA``, ``set_eff()`` and
    ``get_free_parameter_names()``, which
    :py:class:`scalcs.mechanism.Mechanism` provides.

    Two things are deliberately different from fitting:

    * **Nothing is clipped.** The fitter resets a rate that leaves its limits
      before the likelihood sees it, which suits a search. Here, a point
      outside the prior has zero posterior density and is rejected. Clipping
      it to the boundary would instead pile probability onto the boundary.
    * **Failure is a value, not a penalty.** Where the likelihood cannot be
      computed, the log posterior is ``-inf``. A sampler rejects the move, and
      :py:attr:`LogPosterior.nfailures` counts how often that happened.

    Everything is in natural logarithms, because that is what probability
    densities combine in.
"""
__docformat__ = "restructuredtext en"
__all__ = ['UniformPrior', 'LogUniformPrior', 'LogPosterior']

import numpy as np


def _free_rate_limits(mec):
    """ (names, lower, upper) of the free rates, from their limits. """
    names, lower, upper = [], [], []
    missing = []
    for rate in mec.Rates:
        if not rate.is_free:
            continue
        names.append(rate.name)
        if not rate.limits:
            missing.append(rate.name)
            continue
        lo, hi = rate.limits[0]
        lower.append(float(lo))
        upper.append(float(hi))
    if missing:
        raise ValueError("free rate(s) {0} have no limits, so a prior cannot "
                         "be taken from the mechanism; give the bounds "
                         "explicitly".format(missing))
    return tuple(names), np.array(lower), np.array(upper)


class UniformPrior:
    """ Independent uniform distributions on the free rate constants.

        This is the prior of Epstein et al. (2016). It is flat between the
        limits a maximum-likelihood fit would also impose, so the mode of the
        posterior is the maximum-likelihood estimate. The paper uses
        U(1e-2, 1e6) s\\ :sup:`-1` for every rate and U(1e-2, 1e10)
        M\\ :sup:`-1` s\\ :sup:`-1` for association rates, which are the limits
        scalcs gives its rates by default.

        :param lower: Lower bound of each free rate.
        :param upper: Upper bound of each free rate.
        :param names: Names of the free rates, in the same order. Optional.
    """

    def __init__(self, lower, upper, names=None):
        self.lower = np.atleast_1d(np.asarray(lower, dtype=float))
        self.upper = np.atleast_1d(np.asarray(upper, dtype=float))
        if self.lower.shape != self.upper.shape or self.lower.ndim != 1:
            raise ValueError("lower and upper must be matching 1-d arrays")
        if not np.all(np.isfinite(self.lower) & np.isfinite(self.upper)):
            raise ValueError("bounds must be finite for a proper prior")
        if np.any(self.lower >= self.upper):
            raise ValueError("every lower bound must be below its upper bound")
        self.names = None if names is None else tuple(names)
        if self.names is not None and len(self.names) != self.k:
            raise ValueError("{0} names for {1} bounds"
                             .format(len(self.names), self.k))
        self._check_bounds()

    def _check_bounds(self):
        pass

    @classmethod
    def from_mechanism(cls, mec):
        """ The prior bounded by each free rate's own limits.

            :raises ValueError: if a free rate has no limits.
        """
        names, lower, upper = _free_rate_limits(mec)
        return cls(lower, upper, names=names)

    @property
    def k(self):
        """ Number of parameters. """
        return len(self.lower)

    def _rates(self, rates):
        rates = np.asarray(rates, dtype=float)
        if rates.shape != (self.k,):
            raise ValueError("expected {0} rates, got shape {1}"
                             .format(self.k, rates.shape))
        return rates

    def contains(self, rates):
        """ Whether *rates* lie inside the bounds, boundaries included. """
        rates = self._rates(rates)
        return bool(np.all((rates >= self.lower) & (rates <= self.upper)))

    def logpdf(self, rates):
        """ Log density at *rates*: ``-inf`` outside the bounds. """
        if not self.contains(rates):
            return -np.inf
        return float(-np.sum(np.log(self.upper - self.lower)))

    def sample(self, rng, size=None):
        """ Draws from the prior.

            :param rng: A :py:class:`numpy.random.Generator`.
            :param size: Number of draws; None for a single vector.
        """
        shape = (self.k,) if size is None else (size, self.k)
        return rng.uniform(self.lower, self.upper, size=shape)


class LogUniformPrior(UniformPrior):
    """ Independent log-uniform distributions on the free rate constants.

        Uniform in the logarithm of each rate between the bounds: every decade
        is equally likely a priori. This is the scale-free alternative to
        :py:class:`UniformPrior`, which puts almost all of its mass in the top
        decade. The bounds must be positive.
    """

    def _check_bounds(self):
        if np.any(self.lower <= 0):
            raise ValueError("a log-uniform prior needs positive lower bounds")

    def logpdf(self, rates):
        if not self.contains(rates):
            return -np.inf
        rates = self._rates(rates)
        return float(-np.sum(np.log(rates))
                     - np.sum(np.log(np.log(self.upper / self.lower))))

    def sample(self, rng, size=None):
        shape = (self.k,) if size is None else (size, self.k)
        return np.exp(rng.uniform(np.log(self.lower), np.log(self.upper),
                                  size=shape))


class LogPosterior:
    """ The log posterior density of a mechanism's free rate constants.

        Calling it with a vector of free rates returns
        :math:`\\ln p(\\theta) + \\ln L(\\theta)`. The prior is evaluated
        first, and the likelihood is not computed at all where the prior is
        zero.

        :param mec:
          A mechanism carrying its constraints. It is **modified**: every
          evaluation puts its rates on the mechanism, as fitting does.
        :param records:
          A sequence of :py:class:`~HJCFIT.likelihood.fitting.Record`, fitted
          simultaneously, each at its own concentration.
        :param prior:
          Anything with ``logpdf(rates)`` and ``k``. Defaults to
          :py:meth:`UniformPrior.from_mechanism`.
        :param dict solver:
          Root-finding options for the likelihood; see
          :py:data:`~HJCFIT.likelihood.fitting.SOLVER_OPTIONS`.
    """

    def __init__(self, mec, records, prior=None, solver=None):
        from .fitting import HJCFitter

        # The fitter builds and holds one likelihood per record, and sums
        # them at the mechanism's current rates. Its clipping lives in
        # HJCFitter._apply, which is only reached through a parameter vector;
        # evaluating it with none, below, is the unclipped path.
        self.fitter = HJCFitter(mec, records, solver=solver)
        self.mec = mec
        self.names = tuple(mec.get_free_parameter_names())
        self.prior = UniformPrior.from_mechanism(mec) if prior is None else prior
        if self.prior.k != len(self.names):
            raise ValueError("the prior has {0} parameters but the mechanism "
                             "has {1} free rates".format(self.prior.k,
                                                         len(self.names)))
        prior_names = getattr(self.prior, "names", None)
        if prior_names is not None and tuple(prior_names) != self.names:
            raise ValueError("the prior's parameters {0} are not the "
                             "mechanism's free rates {1}, in that order"
                             .format(list(prior_names), list(self.names)))
        #: Likelihood evaluations requested, failures included.
        self.nevals = 0
        #: Likelihood evaluations that could not be computed.
        self.nfailures = 0

    @property
    def k(self):
        """ Number of free parameters. """
        return len(self.names)

    @property
    def records(self):
        return self.fitter.records

    def _rates(self, rates):
        rates = np.asarray(rates, dtype=float)
        if rates.shape != (self.k,):
            raise ValueError("expected {0} rates, got shape {1}"
                             .format(self.k, rates.shape))
        return rates

    def log_prior(self, rates):
        """ :math:`\\ln p(\\theta)`. """
        return float(self.prior.logpdf(self._rates(rates)))

    def log_likelihood(self, rates):
        """ :math:`\\ln L(\\theta)` summed over the records; ``-inf`` if it
            cannot be computed. Rates are used as given, never clipped. """
        rates = self._rates(rates)
        self.nevals += 1
        try:
            if not np.all(np.isfinite(rates)):
                raise ArithmeticError("rates are not finite")
            self.mec.theta_unsqueeze(rates)
            return float(self.fitter.ln_likelihood())
        except (ArithmeticError, ValueError, RuntimeError, FloatingPointError):
            self.nfailures += 1
            return -np.inf

    def __call__(self, rates):
        """ :math:`\\ln p(\\theta) + \\ln L(\\theta)`. """
        log_prior = self.log_prior(rates)
        if log_prior == -np.inf:
            return -np.inf
        return log_prior + self.log_likelihood(rates)
