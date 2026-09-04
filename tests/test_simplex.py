"""The restored Fortran simplex, HJCFIT.likelihood.optimization.simplex.

A simplex that converges to the wrong place still looks like it worked, so
these check it against problems whose answer is known independently, and check
each of the structural features that make it *not* SciPy's Nelder-Mead: the
regular starting simplex built from a fixed step, the vertex-counting
acceptance rule, convergence on the simplex extent alone, and the four-way
ending with restarts.

The reproduced local-search bug is pinned by a test of its own, so that
correcting it cannot happen silently: the published HJCFIT numbers came from
the code with that bug in it, and `local_search="fortran"` exists to reproduce
them.

Pure Python -- nothing here needs the compiled extension.
"""

import numpy as np
import numpy.testing as npt
import pytest

from HJCFIT.likelihood.optimization import (
    SIMPLEX_DEFAULTS,
    SimplexResult,
    simplex,
)


def quadratic(centre):
    c = np.asarray(centre, float)

    def f(x):
        return float(np.sum((np.asarray(x, float) - c) ** 2))
    return f


def rosenbrock(x):
    x = np.asarray(x, float)
    return float(np.sum(100.0 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2))


# ------------------------------------------------------------- it minimises

def test_finds_the_minimum_of_a_quadratic():
    r = simplex(quadratic([1.0, -2.0, 0.5]), [0.0, 0.0, 0.0], errfac=1e-6)
    assert r.success
    npt.assert_allclose(r.x, [1.0, -2.0, 0.5], atol=1e-3)
    assert r.fun < 1e-6


def test_finds_the_rosenbrock_minimum():
    r = simplex(rosenbrock, [-1.2, 1.0], stpfac=np.e, errfac=1e-7,
                maxfev=100000)
    assert r.success
    npt.assert_allclose(r.x, [1.0, 1.0], atol=1e-2)


def test_linear_mode_scales_step_and_tolerance_with_the_guess():
    r = simplex(quadratic([10.0, 20.0]), [8.0, 25.0], logfit=False,
                errfac=1e-8)
    assert r.success
    npt.assert_allclose(r.x, [10.0, 20.0], rtol=1e-3)


def test_log_mode_recovers_rate_constants():
    """The way HJCFIT uses it: search in log space, exponentiate after."""
    true = np.array([50.0, 6000.0, 2.0e8])

    def f(logx):
        return float(np.sum(np.log(np.exp(logx) / true) ** 2))

    r = simplex(f, np.log([10.0, 1000.0, 5.0e7]), errfac=1e-7)
    npt.assert_allclose(np.exp(r.x), true, rtol=1e-3)


def test_args_are_passed_through():
    def f(x, target):
        return float(np.sum((np.asarray(x, float) - target) ** 2))

    r = simplex(f, [0.0, 0.0], args=(np.array([3.0, -1.0]),), errfac=1e-6)
    npt.assert_allclose(r.x, [3.0, -1.0], atol=1e-3)


# --------------------------------------------------- it is not SciPy's NM

def test_starting_simplex_is_regular_and_fixed_in_log_space():
    """Every parameter moves by the same ln(stpfac), whatever its size.

    SciPy perturbs each coordinate by 5 per cent of its own value, which for
    rate constants spanning eight orders of magnitude is a different search.
    """
    seen = []
    f = quadratic([0.0, 0.0])
    simplex(lambda x: (seen.append(np.array(x, float)), f(x))[1],
            [1.0, 1000.0], stpfac=5.0, maxfev=3)
    start = seen[0]
    step = np.log(5.0)
    expected = step * ((np.sqrt(3) - 1) / (2 * np.sqrt(2)) + 1 / np.sqrt(2))
    for offset in (s - start for s in seen[1:3]):
        # identical displacement for both coordinates despite the 1000-fold
        # difference between them
        npt.assert_allclose(max(abs(offset)), expected, rtol=1e-9)


def test_convergence_is_on_simplex_extent_alone():
    """No function tolerance: a tighter errfac really does work harder."""
    f = quadratic([0.0, 0.0])
    tight = simplex(f, [1.0, 1.0], errfac=1e-9)
    loose = simplex(f, [1.0, 1.0], errfac=1e-2)
    assert tight.nfev > loose.nfev
    assert abs(tight.fun) < abs(loose.fun)


def test_the_starting_simplex_costs_k_plus_one_evaluations():
    f = quadratic(np.zeros(4))
    r = simplex(f, np.ones(4), maxfev=4)
    assert r.nfev >= 5


# --------------------------------------------------------- the four endings

def test_returns_one_of_the_documented_endings():
    r = simplex(quadratic([1.0, 1.0]), [0.0, 0.0])
    assert r.iconv in (2, 3, 4, 5)
    assert r.message.startswith("return with ")


def test_budget_exhaustion_is_reported_not_hidden():
    r = simplex(rosenbrock, [-5.0, 5.0], errfac=1e-14, maxfev=60)
    assert r.iconv == 6
    assert not r.success
    assert "no convergence" in r.message
    assert r.nfev > 60


def test_restarts_are_counted_and_capped():
    r = simplex(rosenbrock, [-1.2, 1.0], stpfac=np.e, errfac=1e-9,
                nresmax=2, maxfev=100000)
    assert r.nrestarts <= 2


def test_no_restarts_when_nresmax_is_zero():
    r = simplex(rosenbrock, [-1.2, 1.0], stpfac=np.e, errfac=1e-6, nresmax=0)
    assert r.nrestarts == 0


# ------------------------------------------------------ the reproduced quirk

def test_fortran_local_search_evaluates_the_unchanged_vertex():
    """Pin the bug the published results were produced with.

    SIMPHJC.FOR line 663 evaluates FUNC(kt, temp, ...) where it has just set
    pnew1(j) = temp(j) + crtstp(j). If this test starts failing, someone has
    quietly corrected the port and its numbers are no longer those the Fortran
    would have given.
    """
    seen = []

    def f(x):
        seen.append(np.array(x, float))
        return float(np.sum(np.asarray(x, float) ** 2))

    r = simplex(f, [1.0, 1.0], errfac=1e-3, local_search="fortran",
                nresmax=0)
    duplicates = sum(1 for s in seen if np.allclose(s, r.x, atol=1e-12))
    assert duplicates >= 2, "the converged vertex should be re-measured"


def test_corrected_local_search_actually_perturbs():
    seen = []

    def f(x):
        seen.append(np.array(x, float))
        return float(np.sum((np.asarray(x, float) - 0.3) ** 2))

    r = simplex(f, [1.0, 1.0], errfac=1e-3, local_search="corrected",
                nresmax=0)
    crt = 1e-3
    moved = [s for s in seen
             if np.isclose(abs(s[0] - r.x[0]), crt, atol=1e-9)
             or np.isclose(abs(s[1] - r.x[1]), crt, atol=1e-9)]
    assert moved, "corrected search never stepped a coordinate by crtstp"


def test_local_search_off_is_cheaper_and_lands_in_the_same_place():
    f = quadratic([0.4, -0.7])
    with_ls = simplex(f, [0.0, 0.0], errfac=1e-4, local_search="corrected",
                      nresmax=0)
    without = simplex(f, [0.0, 0.0], errfac=1e-4, local_search="off",
                      nresmax=0)
    assert without.nfev < with_ls.nfev
    npt.assert_allclose(without.fun, with_ls.fun, atol=1e-6)


def test_local_search_mode_is_validated():
    with pytest.raises(ValueError):
        simplex(quadratic([0.0]), [1.0], local_search="nonsense")


def test_empty_x0_is_rejected():
    with pytest.raises(ValueError):
        simplex(quadratic([]), [])


# ---------------------------------------------------------- the result type

def test_result_uses_scipy_field_names():
    """So a script can swap minimize for simplex and keep reading res.x."""
    r = simplex(quadratic([1.0, 2.0]), [0.0, 0.0], errfac=1e-6)
    for name in ("x", "fun", "nfev", "nit", "success", "message"):
        assert hasattr(r, name), name
    assert isinstance(r, SimplexResult)
    npt.assert_allclose(r.fun, quadratic([1.0, 2.0])(r.x), rtol=1e-9)
    assert r.nfev > 0 and r.nit > 0


def test_callback_sees_the_best_vertex():
    seen = []
    simplex(quadratic([2.0, 2.0]), [0.0, 0.0], callback=seen.append,
            errfac=1e-4)
    assert len(seen) > 1
    npt.assert_allclose(seen[-1], [2.0, 2.0], atol=1e-2)


def test_defaults_are_the_fortran_ones():
    assert SIMPLEX_DEFAULTS["stpfac_log"] == 5.0
    assert SIMPLEX_DEFAULTS["stpfac_lin"] == 0.2
    assert SIMPLEX_DEFAULTS["confac"] == 0.5
    assert SIMPLEX_DEFAULTS["resfac"] == 10.0
    assert SIMPLEX_DEFAULTS["nresmax"] == 3
    assert SIMPLEX_DEFAULTS["errfac"] == 1e-3
    assert SIMPLEX_DEFAULTS["maxfev"] == 20000
