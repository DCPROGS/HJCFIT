"""The fitting layer: Record, FitResult, HJCFitter.

Lifted from the reproduction of Colquhoun, Hatton & Hawkes (2003), where it
was exercised over ten scenarios and 250 Monte-Carlo sets.

Nothing here imports scalcs, and that is deliberate rather than incidental:
HJCFitter is duck-typed on a mechanism, so a fake one built in this file
exercises the whole search. If these tests ever need scalcs, the fitter has
grown a dependency the likelihood was meant not to have -- see
test_dcio_integration.py for the real-mechanism fit, which runs in the job
that installs it.
"""

import numpy as np
import pytest

from HJCFIT.likelihood.fitting import (
    FAILURE_COST, FitResult, HJCFitter, Record, trim_to_openings)


# --------------------------------------------------------------------------
# Record
# --------------------------------------------------------------------------

class TestRecord:

    def test_counts(self):
        r = Record(conc=1e-7, groups=((0.1, 0.01, 0.2), (0.3,)), tres=1e-4)
        assert r.n_intervals == 4
        assert r.n_openings == 3          # 2 in the first group, 1 in the second

    def test_as_lists_is_what_the_likelihood_wants(self):
        r = Record(conc=1e-7, groups=((0.1, 0.01, 0.2),), tres=1e-4)
        out = r.as_lists()
        assert isinstance(out, list) and isinstance(out[0], list)
        assert out == [[0.1, 0.01, 0.2]]

    def test_check_rejects_an_even_group(self):
        """An even group does not fail loudly: it silently asks the likelihood
        for a product of matrices that does not alternate."""
        with pytest.raises(ValueError, match="odd number of intervals"):
            Record(conc=1e-7, groups=((0.1, 0.01),), tres=1e-4).check()

    def test_check_names_the_offending_groups(self):
        bad = Record(conc=1e-7, groups=((0.1,), (0.1, 0.01), (0.2,)),
                     tres=1e-4)
        with pytest.raises(ValueError, match=r"\[1\]"):
            bad.check()

    def test_check_rejects_an_empty_record(self):
        with pytest.raises(ValueError, match="no groups"):
            Record(conc=1e-7, groups=(), tres=1e-4).check()

    def test_check_returns_self_so_it_can_be_chained(self):
        r = Record(conc=1e-7, groups=((0.1,),), tres=1e-4)
        assert r.check() is r

    def test_str_says_which_vectors(self):
        chs = Record(conc=1e-7, groups=((0.1,),), tres=1e-4, tcrit=4e-3)
        eq = Record(conc=1e-7, groups=((0.1,),), tres=1e-4, tcrit=None)
        assert "CHS vectors" in str(chs)
        assert "equilibrium vectors" in str(eq)

    def test_str_scales_the_concentration(self):
        assert "nM" in str(Record(conc=1e-7, groups=((0.1,),), tres=1e-4))
        assert "uM" in str(Record(conc=1e-5, groups=((0.1,),), tres=1e-4))


class TestTrimToOpenings:

    def test_drops_a_leading_and_trailing_shut(self):
        got = trim_to_openings([0.05, 0.1, 0.01, 0.2, 0.3],
                               [0.0, 5.0, 0.0, 5.0, 0.0])
        np.testing.assert_allclose(got, [0.1, 0.01, 0.2])

    def test_leaves_an_already_odd_record_alone(self):
        got = trim_to_openings([0.1, 0.01, 0.2], [5.0, 0.0, 5.0])
        np.testing.assert_allclose(got, [0.1, 0.01, 0.2])

    def test_result_is_always_odd(self):
        rng = np.random.default_rng(0)
        for _ in range(50):
            n = int(rng.integers(3, 30))
            amps = np.where(np.arange(n) % 2 == int(rng.integers(0, 2)), 5.0, 0.0)
            assert len(trim_to_openings(rng.random(n), amps)) % 2 == 1


# --------------------------------------------------------------------------
# A mechanism, faked
# --------------------------------------------------------------------------
#
# Everything HJCFitter asks of a mechanism, and nothing else. The point is
# that this is enough: the fitter never imports a mechanism library.

class FakeRate:
    def __init__(self, name, value, free=True, limits=None):
        self.name = name
        self.rateconstants = np.array([value])
        self.is_free = free
        self.limits = limits

    def unit_rate(self):
        return float(self.rateconstants[0])


class FakeMechanism:
    """Two states, one open: Q = [[-alpha, alpha], [beta, -beta]].

    Concentration-independent, so ``set_eff`` does nothing -- which is itself
    worth exercising, since the fitter calls it once per record per evaluation.
    """

    kA = 1

    def __init__(self, alpha=1000.0, beta=100.0, limits=None):
        self.Rates = [FakeRate("alpha", alpha, limits=limits),
                       FakeRate("beta", beta, limits=limits)]
        self.set_eff_calls = 0

    def theta(self):
        return np.array([r.unit_rate() for r in self.Rates if r.is_free])

    def theta_unsqueeze(self, values):
        for rate, value in zip([r for r in self.Rates if r.is_free], values):
            rate.rateconstants = np.array([float(value)])

    def get_free_parameter_names(self):
        return [r.name for r in self.Rates if r.is_free]

    def set_eff(self, name, value):
        self.set_eff_calls += 1

    @property
    def Q(self):
        alpha, beta = (r.unit_rate() for r in self.Rates)
        return np.array([[-alpha, alpha], [beta, -beta]])


@pytest.fixture()
def record():
    """A short alternating record, odd-length, from the fake mechanism's
    rough timescales. Its statistical fidelity does not matter: both the
    likelihood and the search see the same intervals."""
    rng = np.random.default_rng(4)
    n = 41
    intervals = np.empty(n)
    intervals[0::2] = rng.exponential(1.0 / 1000.0, len(intervals[0::2]))
    intervals[1::2] = rng.exponential(1.0 / 100.0, len(intervals[1::2]))
    return Record(conc=0.0, groups=(tuple(intervals),), tres=0.0)


# --------------------------------------------------------------------------
# HJCFitter
# --------------------------------------------------------------------------

class TestHJCFitterWithoutAMechanismLibrary:

    def test_it_needs_no_mechanism_library(self, record):
        """The whole search runs against a mechanism defined in this file."""
        fitter = HJCFitter(FakeMechanism(), [record])
        assert np.isfinite(fitter.log10_likelihood())

    def test_scalcs_is_not_imported_by_the_fitting_module(self):
        import subprocess
        import sys

        out = subprocess.run(
            [sys.executable, "-c",
             "import sys, HJCFIT.likelihood.fitting; "
             "print('scalcs' in sys.modules)"],
            capture_output=True, text=True, check=True)
        assert out.stdout.strip() == "False", (
            "importing the fitting layer pulled in scalcs; it is duck-typed "
            "on a mechanism and must stay that way")

    def test_rejects_an_empty_record_list(self):
        with pytest.raises(ValueError, match="no records"):
            HJCFitter(FakeMechanism(), [])

    def test_checks_its_records_on_construction(self):
        even = Record(conc=0.0, groups=((0.1, 0.01),), tres=0.0)
        with pytest.raises(ValueError, match="odd number"):
            HJCFitter(FakeMechanism(), [even])

    def test_ln_and_log10_differ_by_ln10(self, record):
        fitter = HJCFitter(FakeMechanism(), [record])
        assert fitter.ln_likelihood() == pytest.approx(
            fitter.log10_likelihood() * np.log(10.0))

    def test_set_eff_is_called_once_per_record_per_evaluation(self, record):
        mec = FakeMechanism()
        fitter = HJCFitter(mec, [record, record])
        before = mec.set_eff_calls
        fitter.log10_likelihood()
        assert mec.set_eff_calls - before == 2

    def test_records_are_summed_not_averaged(self, record):
        """Two copies of one record should give twice one record's value."""
        one = HJCFitter(FakeMechanism(), [record]).log10_likelihood()
        two = HJCFitter(FakeMechanism(), [record, record]).log10_likelihood()
        assert two == pytest.approx(2.0 * one)


class TestCost:

    def test_cost_is_the_negative_log10_likelihood(self, record):
        fitter = HJCFitter(FakeMechanism(), [record])
        x = np.log(fitter.mec.theta())
        assert fitter.cost(x) == pytest.approx(-fitter.log10_likelihood(x))

    def test_an_unusable_q_matrix_costs_a_finite_penalty(self, record):
        """A Q matrix the likelihood cannot handle raises, or returns nan
        silently. Either way the search must get a number it can walk on."""
        fitter = HJCFitter(FakeMechanism(), [record])
        value = fitter.cost(np.log([1e-300, 1e-300]))
        assert np.isfinite(value)
        if fitter.nfailures:
            assert value == FAILURE_COST

    def test_failures_are_counted_not_swallowed(self, record):
        fitter = HJCFitter(FakeMechanism(), [record])
        fitter.cost(np.log([1e-300, 1e-300]))
        # whatever happened, nfailures and nevals agree with each other
        assert fitter.nevals == 1
        assert fitter.nfailures in (0, 1)

    def test_evaluations_are_stored_only_when_asked(self, record):
        quiet = HJCFitter(FakeMechanism(), [record])
        quiet.cost(np.log(quiet.mec.theta()))
        assert quiet.evaluations == []

        loud = HJCFitter(FakeMechanism(), [record], store_evaluations=True)
        loud.cost(np.log(loud.mec.theta()))
        assert len(loud.evaluations) == 1


class TestLimits:

    def test_a_rate_is_reset_rather_than_searched_out_of_range(self, record):
        """HJCFIT reset a rate that left its range before evaluating (p. 702).
        Searching the rates themselves without it produced four fits in 250
        with negative rate constants."""
        mec = FakeMechanism(limits=[[10.0, 2000.0]])
        fitter = HJCFitter(mec, [record], log_params=False)
        fitter.log10_likelihood(np.array([1e9, 1e9]))
        assert all(r.unit_rate() <= 2000.0 for r in mec.Rates)

        fitter.log10_likelihood(np.array([-5.0, -5.0]))
        assert all(r.unit_rate() >= 10.0 for r in mec.Rates)

    def test_no_limits_means_a_floor_but_no_ceiling(self, record):
        mec = FakeMechanism(limits=None)
        fitter = HJCFitter(mec, [record], log_params=False)
        assert (fitter._upper == np.inf).all()
        assert (fitter._lower > 0).all(), "a floor is always applied"


class TestFit:

    def test_simplex_is_the_default_search(self, record):
        """Asserted by the ending it reports, not by reading the signature.
        simplex_hjc ends on one of four documented candidates and says which;
        SciPy's Nelder-Mead reports something else entirely."""
        from HJCFIT.likelihood.simplex_hjc import ICONV

        default = HJCFitter(FakeMechanism(), [record]).fit(maxfev=300)
        assert any(w in default.message for w in ICONV.values()), default.message

        scipy_run = HJCFitter(FakeMechanism(), [record]).fit(search="scipy",
                                                            maxfev=300)
        assert not any(w in scipy_run.message for w in ICONV.values())

    def test_the_default_matches_asking_for_simplex_explicitly(self, record):
        """The simplex is deterministic given its start, so these are equal to
        the last digit -- which a default that quietly differed would break."""
        a = HJCFitter(FakeMechanism(), [record]).fit(maxfev=300)
        b = HJCFitter(FakeMechanism(), [record]).fit(search="simplex",
                                                     maxfev=300)
        assert a.log10_likelihood == b.log10_likelihood
        assert a.nevals == b.nevals
        np.testing.assert_array_equal(a.free_values, b.free_values)

    @pytest.mark.parametrize("search", ["simplex", "scipy"])
    def test_a_fit_does_not_end_worse_than_it_began(self, record, search):
        fitter = HJCFitter(FakeMechanism(), [record])
        start = fitter.log10_likelihood()
        result = fitter.fit(search=search, maxfev=400)
        assert result.log10_likelihood >= start - 1e-9

    def test_both_searches_reach_the_same_maximum(self, record):
        """Two independent optimisers over the same likelihood. Agreement is
        the point: it is evidence about the likelihood, not the search."""
        a = HJCFitter(FakeMechanism(), [record]).fit(search="simplex",
                                                     maxfev=2000)
        b = HJCFitter(FakeMechanism(), [record]).fit(search="scipy",
                                                     maxfev=2000)
        assert a.log10_likelihood == pytest.approx(b.log10_likelihood, abs=1e-3)

    def test_an_unknown_search_is_rejected(self, record):
        with pytest.raises(ValueError, match="must be 'simplex' or 'scipy'"):
            HJCFitter(FakeMechanism(), [record]).fit(search="powell")

    def test_the_result_reports_rates_not_logarithms(self, record):
        fitter = HJCFitter(FakeMechanism(), [record], log_params=True)
        result = fitter.fit(maxfev=200)
        assert (result.free_values > 0).all(), "these should be rates"
        assert result.free_names == ("alpha", "beta")

    def test_the_result_reads_the_rates_off_the_mechanism(self, record):
        """Not off the optimiser, so that a reset is reflected in what is
        reported and not only in what was fitted."""
        mec = FakeMechanism(limits=[[10.0, 500.0]])
        result = HJCFitter(mec, [record], log_params=False).fit(maxfev=300)
        assert (result.free_values <= 500.0).all()

    def test_path_is_captured_only_when_asked(self, record):
        quiet = HJCFitter(FakeMechanism(), [record]).fit(maxfev=200)
        assert quiet.path == []
        loud = HJCFitter(FakeMechanism(), [record],
                         store_path=True).fit(maxfev=200)
        assert len(loud.path) > 0

    def test_x0_defaults_to_the_mechanism_as_it_stands(self, record):
        mec = FakeMechanism(alpha=1234.0, beta=56.0)
        fitter = HJCFitter(mec, [record], store_evaluations=True)
        fitter.fit(maxfev=5)
        first = np.exp(fitter.evaluations[0][0])
        np.testing.assert_allclose(first, [1234.0, 56.0], rtol=1e-9)

    def test_rates_includes_every_rate_not_only_the_free_ones(self, record):
        mec = FakeMechanism()
        mec.Rates[1].is_free = False
        result = HJCFitter(mec, [record]).fit(maxfev=100)
        assert set(result.rates) == {"alpha", "beta"}
        assert result.free_names == ("alpha",)

    def test_str_is_worth_printing(self, record):
        result = HJCFitter(FakeMechanism(), [record]).fit(maxfev=200)
        text = str(result)
        assert "log10(L)" in text
        assert "alpha" in text and "beta" in text
        assert "evaluations" in text
