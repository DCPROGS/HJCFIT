"""HJCFIT's optional dcio-backed entry points.

read_idealized_bursts and log_bin_edges import dcio (and scalcs) at call time,
so that the likelihood library itself needs nothing but numpy. Nothing
exercised either of them: the delegation added in #180 and the burst reading
reworked in #170 both shipped with no coverage, and CI could not have caught a
break in either because dcio was never installed there.

Run with dcio and scalcs installed; the CI job that does so is dcio-integration.
"""

import sys

import numpy as np
import pytest

import HJCFIT
from HJCFIT.likelihood._methods import (
    ideal_pdf_scale_factor,
    log_bin_edges,
)

pytest.importorskip("dcio", reason="dcio-backed entry points need dcio")


class TestSoftDependency:
    """The whole point of the lazy imports: dcio must not be pulled in just by
    importing the library."""

    def test_importing_hjcfit_does_not_import_dcio(self):
        out = __import__("subprocess").run(
            [sys.executable, "-c",
             "import sys, HJCFIT, HJCFIT.likelihood; "
             "print('dcio' in sys.modules or 'scalcs' in sys.modules)"],
            capture_output=True, text=True, check=True)
        assert out.stdout.strip() == "False", (
            "importing HJCFIT pulled in dcio or scalcs; they must stay lazy")


class TestReadIdealizedBursts:

    def test_ch82_sample(self):
        """Golden values for the shipped CH82 record.

        Any change to dcio's resolution imposition or burst convention moves
        these; that is what they are for."""
        bursts = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
        assert len(bursts) == 572
        assert sum(len(b) for b in bursts) == 1100
        assert sum(b.sum() for b in bursts) == pytest.approx(4.819128, abs=1e-5)

    def test_every_burst_alternates_open_shut(self):
        """Each burst starts and ends on an opening, so its interval count is
        odd. The missed-events likelihood is a product of matrices alternating
        A->F and F->A; an even-length burst would end on a shut."""
        bursts = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
        assert all(len(b) % 2 == 1 for b in bursts)

    def test_tcrit_sign_does_not_change_segmentation(self):
        """A negative tcrit is a flag to Log10Likelihood selecting equilibrium
        vectors over CHS vectors; the magnitude is the critical time, and the
        same number is passed to both."""
        pos = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
        neg = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=-4e-3)
        assert len(pos) == len(neg)
        for a, b in zip(pos, neg):
            np.testing.assert_allclose(a, b)

    def test_named_samples_all_readable(self):
        for name in ("CH82", "CO", "CCO"):
            bursts = HJCFIT.read_idealized_bursts(name, tau=1e-4, tcrit=4e-3)
            assert len(bursts) > 0


class TestLogBinEdges:
    """Delegates to dcio.analysis.histogram (#180)."""

    def test_edges_cover_the_longest_interval(self):
        x = np.random.default_rng(5).exponential(0.005, 2000) + 2.5e-5
        edges, nbdec = log_bin_edges(x, 2.5e-5)
        assert x.max() <= edges[-1]
        assert edges[0] == pytest.approx(2.5e-5)

    def test_nbdec_from_sample_size(self):
        for n, expected in ((200, 5), (800, 8), (2000, 10), (8000, 12)):
            x = np.full(n, 0.01)
            _, nbdec = log_bin_edges(x, 1e-5)
            assert nbdec == expected

    def test_nbdec_bins_span_a_decade(self):
        edges, nbdec = log_bin_edges(np.full(2000, 0.01), 1e-5)
        assert edges[nbdec] / edges[0] == pytest.approx(10.0)


class TestIdealPdfScaleFactor:
    """Stays in HJCFIT: it renormalises an ideal pdf onto the resolved
    intervals from the Q-matrix survival, which is not what EKDIST's
    similarly-named helper computes."""

    def test_unit_resolution_gives_unit_factor(self):
        aa = np.array([[-1000.0]])
        phi = np.array([1.0])
        assert ideal_pdf_scale_factor(0.0, aa, phi) == pytest.approx(1.0)

    def test_factor_exceeds_one_for_positive_resolution(self):
        aa = np.array([[-1000.0]])
        phi = np.array([1.0])
        assert ideal_pdf_scale_factor(1e-4, aa, phi) > 1.0


class TestDwellTimeHistogram:
    """The consumer of log_bin_edges. Needs matplotlib, which HJCFIT does not
    depend on either."""

    def test_returns_axes(self):
        pytest.importorskip("matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        from HJCFIT.likelihood._methods import dwell_time_histogram

        x = np.random.default_rng(6).exponential(0.005, 2000) + 2.5e-5
        ax = dwell_time_histogram(x, 2.5e-5)
        assert ax is not None

    def test_bars_account_for_every_interval(self):
        pytest.importorskip("matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        from HJCFIT.likelihood._methods import dwell_time_histogram

        x = np.random.default_rng(7).exponential(0.005, 1500) + 2.5e-5
        ax = dwell_time_histogram(x, 2.5e-5)
        # the outline is drawn as sqrt(counts); squaring recovers them, and
        # each bar contributes its count twice
        ydata = ax.lines[0].get_ydata()
        assert round((ydata ** 2).sum() / 2) == len(x)


class TestSimplexOnARealLikelihood:
    """simplex_hjc against HJCFIT's own likelihood, not a test function.

    Its 26 unit tests pin the Fortran semantics on quadratics and Rosenbrock.
    None of them touches Log10Likelihood, so nothing would notice if the two
    stopped working together -- which is the gap that let the documented entry
    point sit broken for years.
    """

    @pytest.fixture()
    def bursts(self):
        return HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)

    @staticmethod
    def _q(rates):
        """CH82 with five of its rates free, the rest at their usual values."""
        b1, f1, b2, f2, alpha = rates
        import numpy as np
        Q = np.array([
            [-(b1 + 3000.0), b1, 3000.0, 0.0, 0.0],
            [f1, -(f1 + 500.0), 0.0, 500.0, 0.0],
            [15.0, 0.0, -(15.0 + b2 + 2000.0), b2, 2000.0],
            [0.0, 15000.0, f2, -(15000.0 + f2), 0.0],
            [0.0, 0.0, alpha, 0.0, -alpha],
        ])
        return Q

    def test_it_improves_a_real_likelihood(self, bursts):
        """The point is not the maximum reached -- it is that the search moves
        and returns something finite through the C++ likelihood."""
        from HJCFIT.likelihood import Log10Likelihood
        from HJCFIT.likelihood.optimization import (
            reset_out_of_range, simplex_hjc)

        lik = Log10Likelihood(bursts, nopen=2, tau=1e-4, tcritical=4e-3)
        guess = np.array([50.0, 2.0 / 3.0, 50.0, 15000.0, 10.0])

        def cost(x):
            try:
                value = -lik(self._q(np.exp(x)))
            except Exception:
                return 1e10
            return value if np.isfinite(value) else 1e10

        start = cost(np.log(guess))
        assert np.isfinite(start)

        res = simplex_hjc(reset_out_of_range(cost, lower=1e-12, upper=1e6),
                          np.log(guess), maxfev=600)
        assert res.fun <= start, "the search should not end worse than it began"
        assert res.nfev > len(guess), "it should have evaluated a simplex"
        from HJCFIT.likelihood.simplex_hjc import ICONV
        assert res.iconv in ICONV, res.iconv
        assert ICONV[res.iconv] in res.message, res.message

    def test_it_reports_scipy_field_names(self, bursts):
        """So a script can swap minimize for simplex_hjc and keep reading."""
        from HJCFIT.likelihood.optimization import simplex_hjc

        res = simplex_hjc(lambda x: float((np.asarray(x) ** 2).sum()),
                          np.array([1.0, 2.0]))
        for field in ("x", "fun", "nfev", "nit", "success", "message"):
            assert hasattr(res, field), field


class TestFittingLayerWithARealMechanism:
    """HJCFitter against scalcs' CH82 and a shipped record.

    test_fitting.py exercises the whole search against a mechanism faked in
    that file, which is what proves the fitter needs no mechanism library.
    This is the other half: that the duck-typing actually matches the library
    it was written for.
    """

    @pytest.fixture()
    def fitter(self):
        from scalcs.samples import samples

        from HJCFIT.likelihood.fitting import HJCFitter, Record

        bursts = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
        record = Record(conc=100e-9, tres=1e-4, tcrit=4e-3,
                        groups=tuple(tuple(b) for b in bursts))
        mec = samples.CH82()
        mec.set_eff("c", 100e-9)
        return HJCFitter(mec, [record])

    def test_scalcs_satisfies_the_duck_type(self, fitter):
        """Every attribute the fitter asks of a mechanism."""
        mec = fitter.mec
        for name in ("theta", "theta_unsqueeze", "Rates", "kA", "set_eff",
                     "get_free_parameter_names", "Q"):
            assert hasattr(mec, name), name
        for name in ("is_free", "limits", "name", "rateconstants", "unit_rate"):
            assert hasattr(mec.Rates[0], name), name

    def test_likelihood_at_the_shipped_guess(self, fitter):
        """A golden value. CH82 as scalcs ships it, against CH82.scn as HJCFIT
        ships it, at tres = 100 us and tcrit = 4 ms."""
        assert fitter.log10_likelihood() == pytest.approx(2286.9746, abs=1e-3)

    def test_the_record_is_what_the_likelihood_requires(self, fitter):
        record = fitter.records[0]
        assert len(record.groups) == 572
        assert record.n_intervals == 1100
        assert all(len(g) % 2 == 1 for g in record.groups)

    @pytest.mark.parametrize("search", ["simplex", "scipy"])
    def test_a_fit_improves_the_likelihood(self, fitter, search):
        start = fitter.log10_likelihood()
        result = fitter.fit(search=search, maxfev=1500)
        assert result.log10_likelihood > start
        assert result.nfailures == 0, "CH82 should not defeat the likelihood"
        assert set(result.rates) == set(r.name for r in fitter.mec.Rates)

    def test_both_searches_find_the_same_maximum(self):
        """Two independent optimisers over the same real likelihood. They
        agree on the maximum while disagreeing on where it is, which is the
        alpha-beta ridge the reproduction measured at r = 0.92 -- so the
        likelihood is asserted and the parameters are not."""
        from scalcs.samples import samples

        from HJCFIT.likelihood.fitting import HJCFitter, Record

        bursts = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
        groups = tuple(tuple(b) for b in bursts)

        reached = []
        for search in ("simplex", "scipy"):
            mec = samples.CH82()
            mec.set_eff("c", 100e-9)
            record = Record(conc=100e-9, groups=groups, tres=1e-4, tcrit=4e-3)
            result = HJCFitter(mec, [record]).fit(search=search, maxfev=1500)
            reached.append(result.log10_likelihood)

        assert reached[0] == pytest.approx(reached[1], abs=1e-2), reached
