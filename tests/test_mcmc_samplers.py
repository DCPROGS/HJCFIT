"""Samplers, chains and diagnostics, on targets whose answers are known.

Nothing here computes an HJC likelihood. The samplers take any log density,
so they are checked against distributions with known moments:

* a correlated Gaussian, the shape of the alpha2-beta2 ridge that
  motivates the adaptive sampler;
* a log-normal, which the log-space walk samples correctly only if it
  applies the Jacobian of the transformation;
* a box, whose edges return -inf.

Every run is seeded, so these are deterministic. The tolerances are several
Monte Carlo standard errors wide, so that a different seed would pass too.
"""

import numpy as np
import pytest

from HJCFIT.likelihood.mcmc import (
    Chain, adaptive_sample, autocorrelation, effective_sample_size,
    mwg_sample, significant_lags)


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------

class Gaussian:
    """Correlated 2-d Gaussian, well away from zero so both walks apply."""

    names = ("a", "b")

    def __init__(self, mean=(10.0, 20.0), sd=(1.0, 2.0), rho=0.9):
        self.mean = np.array(mean)
        sd = np.array(sd)
        self.cov = np.array([[sd[0] ** 2, rho * sd[0] * sd[1]],
                             [rho * sd[0] * sd[1], sd[1] ** 2]])
        self.precision = np.linalg.inv(self.cov)
        self.calls = 0

    def __call__(self, x):
        self.calls += 1
        d = np.asarray(x) - self.mean
        return -0.5 * d @ self.precision @ d


class LogNormal:
    """Independent log-normals: ln(theta) ~ N(mu, sigma**2).

    Written as a density over theta, 1/theta included, which is what a
    posterior over rates is. Walking ln(theta) without the Jacobian would
    sample N(mu - sigma**2, sigma**2) instead.
    """

    def __init__(self, mu=(0.0, 2.0), sigma=(0.5, 0.3)):
        self.mu = np.array(mu)
        self.sigma = np.array(sigma)

    def __call__(self, theta):
        theta = np.asarray(theta)
        if np.any(theta <= 0):
            return -np.inf
        y = np.log(theta)
        return float(np.sum(-0.5 * ((y - self.mu) / self.sigma) ** 2 - y))


def box(x):
    x = np.asarray(x)
    return 0.0 if np.all((x >= 0.0) & (x <= 1.0)) else -np.inf


# --------------------------------------------------------------------------
# Metropolis-within-Gibbs
# --------------------------------------------------------------------------

class TestMWG:

    def test_log_space_walk_applies_the_jacobian(self):
        target = LogNormal()
        chain = mwg_sample(target, np.exp(target.mu), n=20000, burnin=2000,
                           rng=1)
        y = np.log(chain.kept())
        np.testing.assert_allclose(y.mean(axis=0), target.mu, atol=0.03)
        np.testing.assert_allclose(y.std(axis=0), target.sigma, rtol=0.06)

    def test_rate_space_walk_recovers_a_gaussian(self):
        target = Gaussian()
        chain = mwg_sample(target, target.mean, n=20000, burnin=2000, rng=2,
                           log_space=False)
        kept = chain.kept()
        np.testing.assert_allclose(kept.mean(axis=0), target.mean, atol=0.15)
        np.testing.assert_allclose(np.cov(kept.T), target.cov, rtol=0.15)

    def test_tuning_shrinks_an_oversized_step(self):
        target = Gaussian()
        chain = mwg_sample(target, target.mean, n=3000, burnin=2000, rng=3,
                           log_space=False, initial_scale=50.0)
        assert np.all(chain.scales[-1] < 10.0)
        assert np.all(chain.scales[-1] == chain.scales[2000])
        rate = chain.acceptance_rate()
        assert np.all((rate > 0.1) & (rate < 0.6)), rate

    def test_scales_are_frozen_after_burnin(self):
        chain = mwg_sample(Gaussian(), (10.0, 20.0), n=600, burnin=300,
                           rng=4, log_space=False, initial_scale=20.0)
        assert np.all(chain.scales[300:] == chain.scales[300])

    def test_per_parameter_records(self):
        chain = mwg_sample(Gaussian(), (10.0, 20.0), n=100, burnin=50, rng=5)
        assert chain.accepted.shape == (100, 2)
        assert chain.scales.shape == (100, 2)
        assert chain.proposals.shape == (100, 2)
        assert chain.acceptance_rate().shape == (2,)

    def test_one_evaluation_per_parameter_per_iteration(self):
        target = Gaussian()
        chain = mwg_sample(target, target.mean, n=100, burnin=0, rng=6)
        assert chain.nevals == target.calls == 100 * 2 + 1

    def test_log_posterior_is_of_the_rates_not_the_walk(self):
        """The Jacobian is the walk's business; the chain reports the
        target's own density at each sample."""
        target = LogNormal()
        chain = mwg_sample(target, np.exp(target.mu), n=50, burnin=0, rng=7)
        for x, lp in zip(chain.samples, chain.log_posterior):
            assert lp == pytest.approx(target(x))

    def test_needs_positive_rates_to_walk_logs(self):
        with pytest.raises(ValueError, match="positive"):
            mwg_sample(Gaussian(), (-1.0, 20.0), n=10, burnin=0)


# --------------------------------------------------------------------------
# Adaptive Metropolis
# --------------------------------------------------------------------------

class TestAdaptive:

    @pytest.mark.parametrize("mixture", ["sum", "choice"])
    def test_recovers_a_correlated_gaussian(self, mixture):
        target = Gaussian()
        chain = adaptive_sample(target, target.mean, n=30000, burnin=5000,
                                rng=11, mixture=mixture)
        kept = chain.kept()
        np.testing.assert_allclose(kept.mean(axis=0), target.mean, atol=0.15)
        np.testing.assert_allclose(np.cov(kept.T), target.cov, rtol=0.12)

    def test_learns_the_correlation(self):
        """The point of the second stage: its proposals line up with the
        posterior ridge (Epstein et al. 2016, Fig. 3 F)."""
        target = Gaussian(rho=0.95)
        chain = adaptive_sample(target, target.mean, n=10000, burnin=5000,
                                rng=12)
        c = chain.proposal_covariance
        assert c[0, 1] / np.sqrt(c[0, 0] * c[1, 1]) == pytest.approx(
            0.95, abs=0.03)

    def test_mixes_better_than_mwg_on_a_ridge(self):
        """Table 2 of that paper, in miniature: on a correlated target the
        adaptive sampler gives more effective samples per iteration."""
        target = Gaussian(rho=0.98)
        mwg = mwg_sample(target, target.mean, n=10000, burnin=2000, rng=13,
                         log_space=False)
        ada = adaptive_sample(target, target.mean, n=10000, burnin=2000,
                              rng=13)
        ess_mwg = effective_sample_size(mwg.kept()[:, 0]) / len(mwg.kept())
        ess_ada = effective_sample_size(ada.kept()[:, 0]) / len(ada.kept())
        assert ess_ada > 2 * ess_mwg, (ess_ada, ess_mwg)

    def test_log_space_walk_applies_the_jacobian(self):
        target = LogNormal()
        chain = adaptive_sample(target, np.exp(target.mu), n=30000,
                                burnin=5000, rng=14, log_space=True)
        y = np.log(chain.kept())
        np.testing.assert_allclose(y.mean(axis=0), target.mu, atol=0.03)
        np.testing.assert_allclose(y.std(axis=0), target.sigma, rtol=0.06)

    def test_initial_covariance_shapes_the_first_steps(self):
        """Before adapt_start the step is (initial_step / sqrt(k)) L0 z. With
        L0 L0^T = diag(100, 1e-4), the steps in the first coordinate are
        10**3 times the size of those in the second."""
        chain = adaptive_sample(lambda x: 0.0, (0.0, 0.0), n=200, burnin=0,
                                rng=15, adapt_start=200,
                                initial_covariance=np.diag([100.0, 1e-4]))
        steps = np.diff(chain.samples, axis=0)
        ratio = steps[:, 0].std() / steps[:, 1].std()
        assert ratio == pytest.approx(1e3, rel=0.2)

    def test_one_evaluation_per_iteration(self):
        target = Gaussian()
        chain = adaptive_sample(target, target.mean, n=300, burnin=100,
                                rng=16)
        assert chain.nevals == target.calls == 301
        assert chain.accepted.shape == (300,)
        assert chain.scales.shape == (300,)

    def test_rejects_an_unknown_mixture(self):
        with pytest.raises(ValueError, match="mixture"):
            adaptive_sample(Gaussian(), (10.0, 20.0), n=10, burnin=0,
                            mixture="both")


# --------------------------------------------------------------------------
# Both samplers
# --------------------------------------------------------------------------

SAMPLERS = [
    pytest.param(lambda f, x0, **kw: mwg_sample(f, x0, log_space=False, **kw),
                 id="mwg"),
    pytest.param(adaptive_sample, id="adaptive"),
]


@pytest.mark.parametrize("sample", SAMPLERS)
class TestBothSamplers:

    def test_a_seed_reproduces_the_chain(self, sample):
        a = sample(Gaussian(), (10.0, 20.0), n=500, burnin=200, rng=21)
        b = sample(Gaussian(), (10.0, 20.0), n=500, burnin=200, rng=21)
        np.testing.assert_array_equal(a.samples, b.samples)
        np.testing.assert_array_equal(a.proposals, b.proposals)

    def test_never_leaves_the_support(self, sample):
        chain = sample(box, (0.5, 0.5), n=3000, burnin=500, rng=22)
        assert np.all((chain.samples >= 0) & (chain.samples <= 1))
        assert np.any((chain.proposals < 0) | (chain.proposals > 1)), (
            "the test needs proposals outside the box to mean anything")
        np.testing.assert_allclose(chain.kept().mean(axis=0), 0.5, atol=0.05)

    def test_refuses_a_start_with_zero_density(self, sample):
        with pytest.raises(ValueError, match="starting point"):
            sample(box, (2.0, 0.5), n=10, burnin=0)

    def test_nan_is_a_rejection(self, sample):
        def target(x):
            return np.nan if x[0] > 0.6 else box(x)
        chain = sample(target, (0.5, 0.5), n=500, burnin=0, rng=23)
        assert np.all(chain.samples[:, 0] <= 0.6)

    def test_names_come_from_the_target(self, sample):
        chain = sample(Gaussian(), (10.0, 20.0), n=10, burnin=0, rng=24)
        assert chain.names == ("a", "b")

    def test_burnin_cannot_exceed_the_run(self, sample):
        with pytest.raises(ValueError, match="burnin"):
            sample(Gaussian(), (10.0, 20.0), n=10, burnin=11)

    def test_callback_sees_the_chain_so_far(self, sample):
        seen = []
        sample(Gaussian(), (10.0, 20.0), n=250, burnin=100, rng=25,
               callback=lambda i, c: seen.append((i, c.n, c.burnin)),
               callback_every=100)
        assert seen == [(100, 100, 100), (200, 200, 100)]


# --------------------------------------------------------------------------
# Chain
# --------------------------------------------------------------------------

class TestChain:

    @pytest.fixture()
    def chain(self):
        return adaptive_sample(Gaussian(), (10.0, 20.0), n=400, burnin=150,
                               rng=31)

    def test_kept_drops_the_burnin(self, chain):
        assert chain.kept().shape == (250, 2)
        np.testing.assert_array_equal(chain.kept(), chain.samples[150:])

    def test_mode_is_the_best_sample(self, chain):
        x, lp = chain.mode()
        assert lp == chain.log_posterior.max()
        assert Gaussian()(x) == pytest.approx(lp)

    def test_save_and_load_round_trip(self, chain, tmp_path):
        path = tmp_path / "chain.npz"
        chain.save(path)
        back = Chain.load(path)
        for name in ("samples", "log_posterior", "proposals", "accepted",
                     "scales", "proposal_covariance"):
            np.testing.assert_array_equal(getattr(back, name),
                                          getattr(chain, name))
        for name in ("burnin", "names", "sampler", "settings", "nevals",
                     "nfailures"):
            assert getattr(back, name) == getattr(chain, name), name


# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------

def ar1(phi, n, seed):
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.standard_normal() / np.sqrt(1 - phi ** 2)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + rng.standard_normal()
    return x


class TestDiagnostics:

    def test_autocorrelation_starts_at_one(self):
        rho = autocorrelation(ar1(0.5, 1000, 1), 5)
        assert rho[0] == 1.0 and len(rho) == 6

    def test_autocorrelation_of_ar1_decays_geometrically(self):
        rho = autocorrelation(ar1(0.8, 200000, 2), 3)
        np.testing.assert_allclose(rho, 0.8 ** np.arange(4), atol=0.01)

    def test_constant_series_is_refused(self):
        with pytest.raises(ValueError, match="constant"):
            autocorrelation(np.ones(10), 2)

    def test_independent_draws_are_worth_their_number(self):
        x = np.random.default_rng(3).standard_normal(20000)
        assert effective_sample_size(x) == pytest.approx(20000, rel=0.1)

    def test_ar1_matches_theory(self):
        """For AR(1), ESS / N = (1 - phi) / (1 + phi)."""
        phi = 0.9
        x = ar1(phi, 100000, 4)
        ess = effective_sample_size(x, max_lag=200)
        assert ess / len(x) == pytest.approx((1 - phi) / (1 + phi), rel=0.15)

    def test_never_more_than_n(self):
        """An anticorrelated series would otherwise claim more samples than
        it has."""
        x = ar1(-0.5, 20000, 5)
        assert effective_sample_size(x) <= len(x)

    def test_significant_lags_grow_with_correlation(self):
        assert significant_lags(ar1(0.0, 20000, 6)) <= 3
        assert significant_lags(ar1(0.95, 20000, 7)) > 20

    def test_one_value_per_column(self):
        x = np.column_stack([ar1(0.0, 5000, 8), ar1(0.9, 5000, 9)])
        ess = effective_sample_size(x)
        assert ess.shape == (2,) and ess[0] > 3 * ess[1]
