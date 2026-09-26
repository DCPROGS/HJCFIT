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

    The samplers are the two stages of that paper:

    * :py:func:`mwg_sample`, Metropolis-within-Gibbs: a random walk on the
      logarithm of one rate at a time, with each step size tuned during burn-in.
      Robust from a poor start; it finds the mode, but mixes slowly along
      correlated directions.
    * :py:func:`adaptive_sample`, adaptive Metropolis (Haario, Saksman &
      Tamminen 2001; Roberts & Rosenthal 2009): a block random walk whose
      covariance is learned from the chain's own history. Started at the mode
      the first stage found, it follows correlations such as the
      :math:`\\alpha_2`-:math:`\\beta_2` ridge.

    Both return a :py:class:`Chain`. :py:func:`effective_sample_size`
    measures what a chain is worth.

    The samplers take any callable returning a log density, not only a
    :py:class:`LogPosterior`, which is how they are tested on targets whose
    answers are known.
"""
__docformat__ = "restructuredtext en"
__all__ = ['UniformPrior', 'LogUniformPrior', 'LogPosterior', 'Chain',
           'mwg_sample', 'adaptive_sample', 'autocorrelation',
           'significant_lags', 'effective_sample_size']

import json
import time
from dataclasses import dataclass, field

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


# --------------------------------------------------------------------------
# Chains
# --------------------------------------------------------------------------

@dataclass
class Chain:
    """ What a sampler produced.

        :param samples: ``(n, k)`` parameter values, one row per iteration, as
          rates even when the sampler walked their logarithms.
        :param log_posterior: ``(n,)`` log density at each row of *samples*.
        :param proposals: ``(n, k)`` the points proposed, accepted or not, as
          rates. Plotting them against *samples* shows how well a proposal
          fits the posterior (Epstein et al. 2016, Fig. 3).
        :param accepted: Whether each proposal was accepted: ``(n,)`` for a
          block sampler, ``(n, k)`` for one parameter at a time.
        :param scales: The step-size multipliers in use at each iteration,
          ``(n,)`` or ``(n, k)`` like *accepted*.
        :param burnin: Iterations during which the sampler tuned itself. They
          are not draws from the posterior; :py:meth:`kept` drops them.
        :param names: The parameters' names, if the target had them.
        :param sampler: Which sampler produced the chain.
        :param settings: The sampler's settings, for provenance.
        :param seconds: Wall-clock time of the run.
        :param nevals: Calls to the log density.
        :param nfailures: Of those, how many the target counted as failures.
        :param proposal_covariance: The proposal covariance at the end of an
          adaptive run, in the space it walked.
    """

    samples: np.ndarray
    log_posterior: np.ndarray
    proposals: np.ndarray
    accepted: np.ndarray
    scales: np.ndarray
    burnin: int
    names: tuple = None
    sampler: str = ""
    settings: dict = field(default_factory=dict)
    seconds: float = 0.0
    nevals: int = 0
    nfailures: int = 0
    proposal_covariance: np.ndarray = None

    @property
    def n(self):
        """ Iterations, burn-in included. """
        return len(self.samples)

    def kept(self):
        """ The samples after burn-in: the draws from the posterior. """
        return self.samples[self.burnin:]

    def acceptance_rate(self, after_burnin=True):
        """ Fraction of proposals accepted; per parameter for one-at-a-time
            sampling. """
        accepted = self.accepted[self.burnin:] if after_burnin else self.accepted
        return accepted.mean(axis=0)

    def mode(self):
        """ The sample with the highest log posterior, and that value. This is
            where Epstein et al. (2016) start the adaptive sampler. """
        best = int(np.argmax(self.log_posterior))
        return self.samples[best].copy(), float(self.log_posterior[best])

    def save(self, path):
        """ Write the chain to a NumPy ``.npz`` file. """
        arrays = dict(samples=self.samples, log_posterior=self.log_posterior,
                      proposals=self.proposals, accepted=self.accepted,
                      scales=self.scales)
        if self.proposal_covariance is not None:
            arrays['proposal_covariance'] = self.proposal_covariance
        meta = dict(burnin=self.burnin, names=self.names,
                    sampler=self.sampler, settings=self.settings,
                    seconds=self.seconds, nevals=self.nevals,
                    nfailures=self.nfailures)
        np.savez_compressed(path, meta=np.array(json.dumps(meta)), **arrays)

    @classmethod
    def load(cls, path):
        """ Read a chain written by :py:meth:`save`. """
        with np.load(path, allow_pickle=False) as data:
            meta = json.loads(str(data['meta']))
            names = meta.pop('names')
            return cls(samples=data['samples'],
                       log_posterior=data['log_posterior'],
                       proposals=data['proposals'],
                       accepted=data['accepted'], scales=data['scales'],
                       proposal_covariance=(data['proposal_covariance']
                                            if 'proposal_covariance' in data
                                            else None),
                       names=None if names is None else tuple(names),
                       **meta)


# --------------------------------------------------------------------------
# Samplers
# --------------------------------------------------------------------------

class _Target:
    """ The log density a sampler walks, in the space it walks it.

        Walking the logarithm of the rates changes the density by the
        Jacobian of the transformation: a density :math:`\\pi(\\theta)` over
        rates is :math:`\\pi(e^y)\\prod_i e^{y_i}` over :math:`y = \\ln\\theta`,
        so :math:`\\sum_i y_i` is added. Leaving it out samples a different
        distribution -- one that favours small rates.
    """

    def __init__(self, log_density, log_space):
        self.log_density = log_density
        self.log_space = log_space
        self.nevals = 0

    def to_state(self, rates):
        rates = np.asarray(rates, dtype=float)
        if self.log_space:
            if np.any(rates <= 0):
                raise ValueError("a sampler over log rates needs positive "
                                 "starting rates")
            return np.log(rates)
        return rates.copy()

    def to_rates(self, state):
        return np.exp(state) if self.log_space else state

    def __call__(self, state):
        self.nevals += 1
        value = float(self.log_density(self.to_rates(state)))
        if np.isnan(value):
            return -np.inf
        if self.log_space and value > -np.inf:
            value += float(np.sum(state))
        return value


def _start(log_density, x0, log_space):
    target = _Target(log_density, log_space)
    state = target.to_state(x0)
    if state.ndim != 1:
        raise ValueError("the starting point must be a 1-d vector")
    log_p = target(state)
    if not np.isfinite(log_p):
        raise ValueError("the log density at the starting point is {0}; a "
                         "chain must start where the posterior is positive"
                         .format(log_p))
    return target, state, log_p


def _accept(rng, log_p_new, log_p):
    """ The Metropolis test for a symmetric proposal. """
    if log_p_new == -np.inf:
        return False
    return bool(np.log(rng.random()) < log_p_new - log_p)


def _tune(scale, window, low, high, step):
    """ Shrink a scale when too few proposals are accepted, grow it when too
        many are. Used during burn-in only, so that the chain that is kept is
        a Markov chain with a fixed kernel. """
    rate = window.mean()
    if rate < low:
        return scale * (1.0 - step)
    if rate > high:
        return scale * (1.0 + step)
    return scale


def _check_run(n, burnin):
    if not 0 <= burnin <= n:
        raise ValueError("burnin must lie between 0 and n")


def _finish(target, log_density, sampler, settings, started, **arrays):
    return Chain(names=getattr(log_density, 'names', None), sampler=sampler,
                 settings=settings, seconds=time.perf_counter() - started,
                 nevals=target.nevals,
                 nfailures=int(getattr(log_density, 'nfailures', 0)),
                 **arrays)


def mwg_sample(log_density, x0, n, burnin, rng=None, *, initial_scale=1.0,
               tune_every=50, acceptance=(0.1, 0.5), tune_step=0.1,
               log_space=True, callback=None, callback_every=1000):
    """ Metropolis-within-Gibbs: a random walk on one parameter at a time.

        Each iteration visits every parameter in turn and proposes
        :math:`y_i' = y_i + s_i z`, :math:`z \\sim N(0, 1)`, on the logarithm
        of the rate (a multiplicative step on the rate itself), accepting it by
        the Metropolis rule. This is the pilot sampler of Epstein et al.
        (2016): multiplicative steps suit rates spanning several decades, and
        it locates the posterior mode reliably from a poor starting point.

        During burn-in, every *tune_every* iterations, each :math:`s_i` is
        multiplied by ``1 - tune_step`` if that parameter's acceptance over the
        window is below ``acceptance[0]``, or by ``1 + tune_step`` if above
        ``acceptance[1]``.

        :param log_density: Callable returning the log posterior at a vector
          of rates, such as a :py:class:`LogPosterior`.
        :param x0: Starting rates; the density there must be positive.
        :param int n: Iterations, each a sweep over all parameters.
        :param int burnin: Iterations during which the step sizes are tuned.
        :param rng: A :py:class:`numpy.random.Generator`, or a seed.
        :param initial_scale: Starting :math:`s_i`, a scalar or one per
          parameter. On log rates, 1 is a step of about a factor of e.
        :param bool log_space: Walk the logarithms of the rates (default) or
          the rates themselves.
        :param callback: Called as ``callback(iteration, chain)`` every
          *callback_every* iterations with the chain so far.
        :returns: A :py:class:`Chain` with per-parameter ``accepted`` and
          ``scales``.
    """
    _check_run(n, burnin)
    rng = np.random.default_rng(rng)
    settings = dict(n=n, burnin=burnin, initial_scale=initial_scale,
                    tune_every=tune_every, acceptance=list(acceptance),
                    tune_step=tune_step, log_space=log_space)
    started = time.perf_counter()
    target, state, log_p = _start(log_density, x0, log_space)
    k = len(state)
    scale = np.broadcast_to(np.asarray(initial_scale, dtype=float), (k,)).copy()
    settings['initial_scale'] = scale.tolist()
    low, high = acceptance

    samples = np.empty((n, k))
    proposals = np.empty((n, k))
    accepted = np.zeros((n, k), dtype=bool)
    scales = np.empty((n, k))
    log_posterior = np.empty(n)

    def partial(upto):
        return _finish(target, log_density, 'mwg', settings, started,
                       samples=target.to_rates(samples[:upto]),
                       log_posterior=log_posterior[:upto],
                       proposals=target.to_rates(proposals[:upto]),
                       accepted=accepted[:upto], scales=scales[:upto],
                       burnin=min(burnin, upto))

    for it in range(n):
        scales[it] = scale
        for i in range(k):
            proposal = state.copy()
            proposal[i] += scale[i] * rng.standard_normal()
            proposals[it, i] = proposal[i]
            log_p_new = target(proposal)
            if _accept(rng, log_p_new, log_p):
                state, log_p = proposal, log_p_new
                accepted[it, i] = True
        samples[it] = state
        # The log density of the rates, without the Jacobian of the walk.
        log_posterior[it] = log_p - (np.sum(state) if log_space else 0.0)

        done = it + 1
        if done <= burnin and done % tune_every == 0:
            window = accepted[done - tune_every:done]
            for i in range(k):
                scale[i] = _tune(scale[i], window[:, i], low, high, tune_step)
        if callback is not None and done % callback_every == 0:
            callback(done, partial(done))

    return partial(n)


class _RunningCovariance:
    """ Mean and covariance of a growing sample, one row at a time (Welford).
        Divides by n - 1, as ``numpy.cov`` does. """

    def __init__(self, k):
        self.count = 0
        self.mean = np.zeros(k)
        self._m2 = np.zeros((k, k))

    def add(self, x):
        self.count += 1
        delta = x - self.mean
        self.mean += delta / self.count
        self._m2 += np.outer(delta, x - self.mean)

    @property
    def covariance(self):
        return self._m2 / (self.count - 1)


def _cholesky(matrix):
    """ Lower Cholesky factor, adding the smallest jitter that makes a
        near-singular sample covariance usable. """
    jitter = 0.0
    size = float(np.mean(np.abs(np.diag(matrix)))) or 1.0
    for _ in range(12):
        try:
            return np.linalg.cholesky(matrix + jitter * np.eye(len(matrix)))
        except np.linalg.LinAlgError:
            jitter = size * 1e-12 if jitter == 0.0 else jitter * 10.0
    raise np.linalg.LinAlgError("proposal covariance is not positive definite")


def adaptive_sample(log_density, x0, n, burnin, rng=None, *, adapt_start=100,
                    beta=0.05, initial_step=0.1, initial_covariance=None,
                    mixture='sum', tune_every=50, acceptance=(0.1, 0.5),
                    tune_step=0.1, log_space=False, callback=None,
                    callback_every=1000):
    """ Adaptive Metropolis: a block random walk that learns its covariance.

        For the first *adapt_start* iterations the proposal is
        :math:`x + (\\epsilon/\\sqrt{k})\\,L_0 z`, where :math:`\\epsilon` is
        *initial_step* and :math:`L_0 L_0^T` is *initial_covariance*, or the
        identity if none is given. From then on :math:`\\hat\\Sigma`, the
        covariance of every sample so far, takes over. It is scaled by
        :math:`2.38^2/k`, which is optimal for a Gaussian target (Roberts &
        Rosenthal 2001). It is also mixed with the small step, so that the
        chain cannot lock onto a degenerate :math:`\\hat\\Sigma`:

        * ``mixture='sum'``, as in Epstein et al. (2016):
          :math:`x' = x + (1-\\beta)\\,L z_1 + \\beta\\,(\\epsilon/\\sqrt{k})\\,L_0 z_2`,
          where :math:`L L^T = (2.38^2/k)\\,s\\,\\hat\\Sigma`.
        * ``mixture='choice'``, as in Roberts & Rosenthal (2009): the small
          step with probability :math:`\\beta`, otherwise the learned one.

        Both are symmetric. During burn-in, a global multiplier :math:`s` on
        :math:`\\hat\\Sigma` is tuned from the acceptance rate every
        *tune_every* iterations, as in :py:func:`mwg_sample`.

        By default the walk is on the rates themselves, as in that paper; start
        it at the mode :py:func:`mwg_sample` found. With ``log_space=True`` it
        walks their logarithms instead.

        :param log_density: Callable returning the log posterior at a vector
          of rates, such as a :py:class:`LogPosterior`.
        :param x0: Starting rates; the density there must be positive.
        :param int n: Iterations.
        :param int burnin: Iterations during which *s* is tuned. The
          covariance keeps adapting afterwards. The weight of each new sample
          in it falls as :math:`1/n`, which is what keeps the chain ergodic
          (Haario et al. 2001).
        :param rng: A :py:class:`numpy.random.Generator`, or a seed.
        :param int adapt_start: Iterations before the learned covariance is
          used.
        :param float beta: Weight of the small step.
        :param initial_covariance: :math:`L_0 L_0^T`. The inverse Hessian at a
          maximum-likelihood fit is a good choice. It shapes the small step
          throughout.
        :param callback: Called as ``callback(iteration, chain)`` every
          *callback_every* iterations with the chain so far.
        :returns: A :py:class:`Chain`; ``proposal_covariance`` is the learned
          covariance at the end, scaled as it was being used.
    """
    if mixture not in ('sum', 'choice'):
        raise ValueError("mixture must be 'sum' or 'choice', not "
                         + repr(mixture))
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must lie between 0 and 1")
    _check_run(n, burnin)
    rng = np.random.default_rng(rng)
    settings = dict(n=n, burnin=burnin, adapt_start=adapt_start, beta=beta,
                    initial_step=initial_step, mixture=mixture,
                    tune_every=tune_every, acceptance=list(acceptance),
                    tune_step=tune_step, log_space=log_space,
                    initial_covariance=initial_covariance is not None)
    started = time.perf_counter()
    target, state, log_p = _start(log_density, x0, log_space)
    k = len(state)
    if initial_covariance is None:
        small = initial_step / np.sqrt(k) * np.eye(k)
    else:
        c0 = np.asarray(initial_covariance, dtype=float)
        if c0.shape != (k, k):
            raise ValueError("initial_covariance must be {0} x {0}".format(k))
        small = initial_step / np.sqrt(k) * _cholesky(c0)
    optimal = 2.38 ** 2 / k
    s = 1.0
    low, high = acceptance

    running = _RunningCovariance(k)
    samples = np.empty((n, k))
    proposals = np.empty((n, k))
    accepted = np.zeros(n, dtype=bool)
    scales = np.empty(n)
    log_posterior = np.empty(n)
    learned_factor = None

    def partial(upto):
        return _finish(target, log_density, 'adaptive', settings, started,
                       samples=target.to_rates(samples[:upto]),
                       log_posterior=log_posterior[:upto],
                       proposals=target.to_rates(proposals[:upto]),
                       accepted=accepted[:upto], scales=scales[:upto],
                       burnin=min(burnin, upto),
                       proposal_covariance=(
                           None if learned_factor is None
                           else learned_factor @ learned_factor.T))

    for it in range(n):
        scales[it] = s
        if it < adapt_start or running.count < 2:
            proposal = state + small @ rng.standard_normal(k)
        else:
            learned_factor = _cholesky(optimal * s * running.covariance)
            learned = learned_factor @ rng.standard_normal(k)
            isotropic = small @ rng.standard_normal(k)
            if mixture == 'sum':
                proposal = state + (1.0 - beta) * learned + beta * isotropic
            elif rng.random() < beta:
                proposal = state + isotropic
            else:
                proposal = state + learned
        proposals[it] = proposal
        log_p_new = target(proposal)
        if _accept(rng, log_p_new, log_p):
            state, log_p = proposal, log_p_new
            accepted[it] = True
        samples[it] = state
        log_posterior[it] = log_p - (np.sum(state) if log_space else 0.0)
        running.add(state)

        done = it + 1
        if done <= burnin and done % tune_every == 0:
            s = _tune(s, accepted[done - tune_every:done], low, high,
                      tune_step)
        if callback is not None and done % callback_every == 0:
            callback(done, partial(done))

    return partial(n)


# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------

def autocorrelation(x, max_lag):
    """ Sample autocorrelation of a series at lags ``0 .. max_lag``.

        The biased estimator: autocovariances divided by the length of the
        series, then by the variance. It is the usual estimator for MCMC
        output, and the one MATLAB's ``autocorr`` computes.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("autocorrelation takes a 1-d series")
    n = len(x)
    max_lag = min(int(max_lag), n - 1)
    centred = x - x.mean()
    variance = centred @ centred / n
    if variance == 0.0:
        raise ValueError("the series is constant")
    return np.array([centred[:n - lag] @ centred[lag:] / n / variance
                     for lag in range(max_lag + 1)])


def significant_lags(x, max_lag=200, nsd=2.0):
    """ The first lag whose autocorrelation falls below
        :math:`\\text{nsd}/\\sqrt{N}`, the approximate bound for a series with
        no autocorrelation; *max_lag* if none does.

        Epstein et al. (2016, Table 2) report this as the number of
        significant lags, and truncate the effective sample size there.
    """
    x = np.asarray(x, dtype=float)
    rho = autocorrelation(x, max_lag)
    below = np.nonzero(rho < nsd / np.sqrt(len(x)))[0]
    return int(below[0]) if len(below) else len(rho) - 1


def effective_sample_size(x, max_lag=None):
    """ Effective sample size by Geyer's (1992) initial monotone sequence.

        Autocorrelations are summed in adjacent pairs,
        :math:`\\Gamma_m = \\rho_{2m} + \\rho_{2m+1}`. The pairs are made
        non-increasing and kept while positive, and then
        :math:`\\text{ESS} = N / (-1 + 2\\sum_m \\Gamma_m)`. The denominator
        is floored at 1, so the estimate never exceeds *N*. This is the
        estimator of Epstein et al. (2016).

        :param x: A series, or an ``(n, k)`` array of them.
        :param max_lag: Autocorrelations beyond this lag are not used. None
          takes each series' own :py:func:`significant_lags`. That paper
          truncates every parameter at the significant lags of
          :math:`\\alpha_2`; pass that number to do the same.
        :returns: A float, or one per column.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim == 2:
        return np.array([effective_sample_size(col, max_lag) for col in x.T])
    lag = significant_lags(x) if max_lag is None else int(max_lag)
    rho = autocorrelation(x, lag)
    npairs = len(rho) // 2
    gamma = rho[0:2 * npairs:2] + rho[1:2 * npairs:2]
    gamma = np.minimum.accumulate(gamma)
    tau = max(1.0, -rho[0] + 2.0 * gamma[gamma > 0].sum())
    return len(x) / tau
