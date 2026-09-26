"""Bayesian layer: priors and the log posterior.

As in test_fitting.py, nothing here imports scalcs. The posterior is
duck-typed on a mechanism, so the fake one from test_fitting.py is enough;
the real-mechanism checks live in test_dcio_integration.py.
"""

import subprocess
import sys

import numpy as np
import pytest

from HJCFIT.likelihood.fitting import SOLVER_OPTIONS, HJCFitter, Record
from HJCFIT.likelihood.mcmc import LogPosterior, LogUniformPrior, UniformPrior

from test_fitting import FakeMechanism


@pytest.fixture()
def record():
    """A short alternating record from the fake mechanism's rough timescales,
    as in test_fitting.py."""
    rng = np.random.default_rng(4)
    n = 41
    intervals = np.empty(n)
    intervals[0::2] = rng.exponential(1.0 / 1000.0, len(intervals[0::2]))
    intervals[1::2] = rng.exponential(1.0 / 100.0, len(intervals[1::2]))
    return Record(conc=0.0, groups=(tuple(intervals),), tres=0.0)


LIMITS = [[1.0, 1e5]]


def test_mcmc_does_not_import_scalcs():
    out = subprocess.run(
        [sys.executable, "-c",
         "import sys, HJCFIT.likelihood.mcmc; print('scalcs' in sys.modules)"],
        capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"


# --------------------------------------------------------------------------
# Priors
# --------------------------------------------------------------------------

class TestUniformPrior:

    def test_density_inside_is_the_product_of_reciprocal_widths(self):
        prior = UniformPrior([0.0, 10.0], [2.0, 30.0])
        assert prior.logpdf([1.0, 20.0]) == pytest.approx(-np.log(2 * 20))

    def test_boundaries_are_inside(self):
        prior = UniformPrior([0.0, 10.0], [2.0, 30.0])
        assert np.isfinite(prior.logpdf([0.0, 30.0]))

    def test_outside_is_minus_infinity(self):
        prior = UniformPrior([0.0, 10.0], [2.0, 30.0])
        assert prior.logpdf([1.0, 31.0]) == -np.inf
        assert prior.logpdf([-1e-9, 20.0]) == -np.inf

    def test_from_mechanism_reads_the_free_rates_limits(self):
        prior = UniformPrior.from_mechanism(FakeMechanism(limits=LIMITS))
        assert prior.names == ("alpha", "beta")
        np.testing.assert_array_equal(prior.lower, [1.0, 1.0])
        np.testing.assert_array_equal(prior.upper, [1e5, 1e5])

    def test_from_mechanism_needs_limits(self):
        with pytest.raises(ValueError, match="no limits"):
            UniformPrior.from_mechanism(FakeMechanism())

    def test_rejects_inverted_or_infinite_bounds(self):
        with pytest.raises(ValueError, match="below"):
            UniformPrior([1.0], [1.0])
        with pytest.raises(ValueError, match="finite"):
            UniformPrior([0.0], [np.inf])

    def test_rejects_the_wrong_number_of_rates(self):
        with pytest.raises(ValueError, match="expected 2"):
            UniformPrior([0.0, 0.0], [1.0, 1.0]).logpdf([0.5])

    def test_samples_stay_inside(self):
        prior = UniformPrior([0.0, 10.0], [2.0, 30.0])
        draws = prior.sample(np.random.default_rng(1), size=1000)
        assert draws.shape == (1000, 2)
        assert all(prior.contains(d) for d in draws)


class TestLogUniformPrior:

    def test_integrates_to_one(self):
        prior = LogUniformPrior([1e-2], [1e6])
        x = np.logspace(-2, 6, 200001)
        density = np.exp([prior.logpdf([v]) for v in x])
        assert np.trapezoid(density, x) == pytest.approx(1.0, rel=1e-4)

    def test_every_decade_equally_likely(self):
        prior = LogUniformPrior([1.0], [1e4])
        draws = prior.sample(np.random.default_rng(2), size=40000)[:, 0]
        counts = np.histogram(np.log10(draws), bins=4, range=(0, 4))[0]
        assert np.all(np.abs(counts / 10000 - 1) < 0.05)

    def test_needs_positive_bounds(self):
        with pytest.raises(ValueError, match="positive"):
            LogUniformPrior([0.0], [1.0])


# --------------------------------------------------------------------------
# Solver options on the fitter
# --------------------------------------------------------------------------

class TestSolverOptions:

    def test_reach_every_likelihood(self, record):
        fitter = HJCFitter(FakeMechanism(), [record, record],
                           solver=dict(nmax=2, xtol=1e-12, rtol=1e-11,
                                       itermax=50))
        for lik in fitter.likelihoods:
            assert (lik.nmax, lik.xtol, lik.rtol, lik.itermax) == (
                2, 1e-12, 1e-11, 50)

    def test_defaults_are_left_alone(self, record):
        lik = HJCFitter(FakeMechanism(), [record]).likelihoods[0]
        assert (lik.nmax, lik.xtol, lik.rtol) == (3, 1e-10, 1e-10)

    def test_unknown_option_is_refused(self, record):
        with pytest.raises(ValueError, match="unknown solver option"):
            HJCFitter(FakeMechanism(), [record], solver=dict(tol=1e-9))

    def test_the_names_are_what_log10likelihood_takes(self):
        assert SOLVER_OPTIONS == ('nmax', 'xtol', 'rtol', 'itermax',
                                  'lower_bound', 'upper_bound')


# --------------------------------------------------------------------------
# LogPosterior
# --------------------------------------------------------------------------

class TestLogPosterior:

    def test_is_prior_plus_likelihood(self, record):
        mec = FakeMechanism(limits=LIMITS)
        post = LogPosterior(mec, [record])
        rates = np.array([900.0, 120.0])
        fitter = HJCFitter(FakeMechanism(limits=LIMITS), [record])
        expected = fitter.ln_likelihood(np.log(rates)) - 2 * np.log(1e5 - 1.0)
        assert post(rates) == pytest.approx(expected, rel=1e-12)

    def test_names_and_size_come_from_the_mechanism(self, record):
        post = LogPosterior(FakeMechanism(limits=LIMITS), [record])
        assert post.names == ("alpha", "beta") and post.k == 2

    def test_outside_the_prior_skips_the_likelihood(self, record):
        post = LogPosterior(FakeMechanism(limits=LIMITS), [record])
        assert post([2e5, 100.0]) == -np.inf
        assert post.nevals == 0

    def test_the_likelihood_is_never_clipped(self, record):
        """The fitter would reset alpha = 3000 to its limit, 2000. The
        posterior's likelihood must see 3000 itself."""
        limits = [[1.0, 2000.0]]
        post = LogPosterior(FakeMechanism(limits=limits), [record])
        unclipped = HJCFitter(FakeMechanism(), [record])
        clipped = HJCFitter(FakeMechanism(limits=limits), [record])
        rates = np.array([3000.0, 100.0])
        value = post.log_likelihood(rates)
        assert value == pytest.approx(
            unclipped.ln_likelihood(np.log(rates)), rel=1e-12)
        assert value != pytest.approx(
            clipped.ln_likelihood(np.log(rates)), rel=1e-6)

    def test_a_failure_is_minus_infinity_and_counted(self, record):
        post = LogPosterior(FakeMechanism(limits=LIMITS), [record],
                            prior=UniformPrior([-np.float64(1e9)] * 2,
                                               [np.float64(1e9)] * 2))
        assert post.log_likelihood([np.nan, 100.0]) == -np.inf
        assert (post.nevals, post.nfailures) == (1, 1)

    def test_a_zero_rate_does_not_raise(self, record):
        """A rate of zero gives a singular Q; the posterior must return
        -inf or a number, never an exception."""
        post = LogPosterior(FakeMechanism(limits=[[0.0, 1e5]]), [record])
        value = post([0.0, 100.0])
        assert value == -np.inf or np.isfinite(value)

    def test_solver_options_reach_the_likelihood(self, record):
        post = LogPosterior(FakeMechanism(limits=LIMITS), [record],
                            solver=dict(nmax=2))
        assert post.fitter.likelihoods[0].nmax == 2

    def test_explicit_prior_must_match_the_mechanism(self, record):
        with pytest.raises(ValueError, match="free rates"):
            LogPosterior(FakeMechanism(limits=LIMITS), [record],
                         prior=UniformPrior([0.0], [1.0]))
        with pytest.raises(ValueError, match="in that order"):
            LogPosterior(FakeMechanism(limits=LIMITS), [record],
                         prior=UniformPrior([0, 0], [1, 1],
                                            names=("beta", "alpha")))

    def test_wrong_number_of_rates_is_an_error_not_minus_infinity(
            self, record):
        post = LogPosterior(FakeMechanism(limits=LIMITS), [record])
        with pytest.raises(ValueError, match="expected 2"):
            post([100.0])
