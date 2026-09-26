"""The Hessian and the Gaussian approximation at a mode.

Checked on functions whose second derivatives are known exactly. The
scaled case matters most: rate constants span eight decades, so one
parameter may be ~1 and another ~1e8 in the same Hessian.
"""

import numpy as np
import pytest

from HJCFIT.likelihood.mcmc import (
    GaussianApproximation, LogPosterior, gaussian_approximation, hessian)


def gaussian_log_density(mean, cov):
    precision = np.linalg.inv(cov)

    def f(x):
        d = np.asarray(x) - mean
        return -0.5 * d @ precision @ d
    return f


def scaled_cov(sd, rho):
    sd = np.asarray(sd, dtype=float)
    corr = np.full((len(sd), len(sd)), rho)
    np.fill_diagonal(corr, 1.0)
    return corr * np.outer(sd, sd)


class TestHessian:

    def test_exact_on_a_quadratic(self):
        cov = scaled_cov([1.0, 2.0, 0.5], 0.6)
        f = gaussian_log_density(np.array([3.0, 5.0, 7.0]), cov)
        h, _ = hessian(f, [3.0, 5.0, 7.0])
        np.testing.assert_allclose(h, -np.linalg.inv(cov), rtol=1e-7)

    def test_fourth_order_on_a_curved_function(self):
        """f = x0^4 + x0 x1^2 + exp(x1), whose Hessian is known exactly.
        Richardson extrapolation should do much better than the plain
        central difference at the same step."""
        def f(x):
            return x[0] ** 4 + x[0] * x[1] ** 2 + np.exp(x[1])
        x = np.array([1.3, 0.7])
        exact = np.array([[12 * x[0] ** 2, 2 * x[1]],
                          [2 * x[1], 2 * x[0] + np.exp(x[1])]])
        h, err = hessian(f, x, rel_step=1e-2)
        np.testing.assert_allclose(h, exact, rtol=1e-7)
        assert np.all(err < 1e-2 * np.abs(exact))

    def test_rates_eight_decades_apart(self):
        """One rate ~1, one ~1e8, correlated. Relative steps handle both."""
        mean = np.array([2.0, 4e7, 5e4])
        cov = scaled_cov([0.2, 4e6, 3e3], 0.9)
        h, _ = hessian(gaussian_log_density(mean, cov), mean)
        np.testing.assert_allclose(h, -np.linalg.inv(cov), rtol=1e-6)

    def test_absolute_steps_where_a_coordinate_is_zero(self):
        f = gaussian_log_density(np.zeros(2), np.eye(2))
        with pytest.raises(ValueError, match="absolute steps"):
            hessian(f, [0.0, 0.0])
        h, _ = hessian(f, [0.0, 0.0], steps=1e-3)
        np.testing.assert_allclose(h, -np.eye(2), rtol=1e-8)

    def test_a_step_across_the_boundary_is_an_error(self):
        """A difference across a prior edge gives -inf, never a number."""
        def f(x):
            return -np.sum(x ** 2) if np.all(x < 1.0) else -np.inf
        with pytest.raises(ValueError, match="parameter"):
            hessian(f, [0.9995, 0.5])

    def test_evaluations(self):
        calls = [0]

        def f(x):
            calls[0] += 1
            return -np.sum(np.asarray(x) ** 2)
        hessian(f, np.ones(4))
        k = 4
        assert calls[0] == 1 + 2 * (2 * k + 4 * k * (k - 1) // 2)


class TestGaussianApproximation:

    def test_recovers_the_covariance_of_a_gaussian(self):
        mean = np.array([2.0, 4e7, 5e4])
        cov = scaled_cov([0.2, 4e6, 3e3], 0.9)
        g = gaussian_approximation(gaussian_log_density(mean, cov), mean)
        np.testing.assert_allclose(g.covariance, cov, rtol=1e-6)
        np.testing.assert_allclose(g.sd, [0.2, 4e6, 3e3], rtol=1e-6)
        np.testing.assert_allclose(g.correlation[0, 1], 0.9, rtol=1e-6)
        np.testing.assert_array_equal(g.mode, mean)
        assert g.log_density == pytest.approx(0.0)

    def test_covariance_is_symmetric_positive_definite(self):
        cov = scaled_cov([1.0, 3.0, 0.1, 50.0], 0.95)
        g = gaussian_approximation(
            gaussian_log_density(np.array([1.0, 2.0, 3.0, 4.0]), cov),
            [1.0, 2.0, 3.0, 4.0])
        np.testing.assert_array_equal(g.covariance, g.covariance.T)
        assert np.all(np.linalg.eigvalsh(g.covariance) > 0)

    def test_refuses_a_point_that_is_not_a_maximum(self):
        def saddle(x):
            return -x[0] ** 2 + x[1] ** 2
        with pytest.raises(ValueError, match="not a maximum"):
            gaussian_approximation(saddle, [1.0, 1.0])

    def test_marginal_pdf_integrates_to_one(self):
        g = gaussian_approximation(
            gaussian_log_density(np.array([10.0]), np.array([[4.0]])), [10.0])
        x = np.linspace(-10, 30, 20001)
        assert np.trapezoid(g.marginal_pdf(0, x), x) == pytest.approx(1.0)

    def test_names_come_from_the_density(self):
        f = gaussian_log_density(np.array([1.0, 2.0]), np.eye(2))
        f_named = lambda x: f(x)  # noqa: E731
        f_named.names = ("alpha", "beta")
        g = gaussian_approximation(f_named, [1.0, 2.0])
        assert g.names == ("alpha", "beta")
        assert g.marginal_pdf("beta", 2.0) == pytest.approx(
            1 / np.sqrt(2 * np.pi))

    def test_samples_have_its_moments(self):
        cov = scaled_cov([1.0, 2.0], 0.8)
        g = gaussian_approximation(
            gaussian_log_density(np.array([5.0, 6.0]), cov), [5.0, 6.0])
        draws = g.sample(np.random.default_rng(1), size=40000)
        np.testing.assert_allclose(np.cov(draws.T), cov, rtol=0.05)

    def test_is_a_dataclass_worth_keeping(self):
        g = gaussian_approximation(
            gaussian_log_density(np.array([5.0]), np.eye(1)), [5.0])
        assert isinstance(g, GaussianApproximation)
        assert g.hessian_error.shape == (1, 1)


def test_on_a_log_posterior():
    """With the fake two-state mechanism: the approximation at a fitted
    maximum is a proper covariance."""
    from HJCFIT.likelihood.fitting import HJCFitter, Record
    from test_fitting import FakeMechanism

    rng = np.random.default_rng(4)
    intervals = np.empty(401)
    intervals[0::2] = rng.exponential(1.0 / 1000.0, 201)
    intervals[1::2] = rng.exponential(1.0 / 100.0, 200)
    record = Record(conc=0.0, groups=(tuple(intervals),), tres=0.0)
    fit = HJCFitter(FakeMechanism(), [record]).fit(search="scipy")
    post = LogPosterior(FakeMechanism(limits=[[1.0, 1e5]]), [record])
    g = gaussian_approximation(post, fit.free_values)
    assert g.names == ("alpha", "beta")
    assert np.all(g.sd > 0)
    # For an exponential rate from n intervals, the standard error is about
    # rate / sqrt(n); here n = 200 openings and 200 shuttings.
    np.testing.assert_allclose(g.sd / fit.free_values, 1 / np.sqrt(200),
                               rtol=0.25)


def test_relative_error_reports_how_well_the_hessian_is_determined():
    """-cosh(x - 2) has its maximum at 2 with H = -1. A huge step leaves a
    visible truncation error; the default step leaves almost none."""
    def f(x):
        return -np.sum(np.cosh(np.asarray(x) - 2.0))
    coarse = gaussian_approximation(f, [2.0, 2.0], rel_step=0.3)
    fine = gaussian_approximation(f, [2.0, 2.0])
    assert fine.relative_error < 1e-6
    assert coarse.relative_error > 100 * fine.relative_error
    np.testing.assert_allclose(fine.sd, 1.0, rtol=1e-6)
