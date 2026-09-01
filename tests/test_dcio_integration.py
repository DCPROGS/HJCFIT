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
