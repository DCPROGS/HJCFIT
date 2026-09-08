"""The C++ likelihood against the 2003 Fortran one, on the same Q and record.

``HJCFIT`` reimplements in C++ what ``hjclik.for`` did in Fortran: the exact
missed-events likelihood of Hawkes, Jalali & Colquhoun. Everything else in the
suite checks the C++ against closed forms, against limits, or against itself
by another route. This checks it against **the program that produced the
published results** -- Colquhoun, Hatton & Hawkes (2003) among them -- which
is a different kind of evidence, and the only kind that can catch an error the
two implementations of a formula would share.

The Fortran lives in ``tests/dcfortran/``, vendored unmodified from
``DCPROGS/DCFORTRAN``, with the three compiler-compatibility changes declared
in ``tests/dcfortran/build.py``. It is not built by default::

    python tests/dcfortran/build.py
    pytest tests/test_fortran_likelihood.py -v

and every test here skips if it has not been. See ``tests/dcfortran/README.md``.

**What is compared.** Both sides are handed the same Q matrix -- not rate
constants -- so neither side's parameterisation, topology handling nor
constraint machinery enters; and the same apparent record. The Fortran's own
input is a whole record which it divides into bursts itself
(``hjclik.for`` line 1009); the C++ takes the bursts. Those two segmentations
were shown to agree interval for interval over 83 000 bursts, and the check is
kept in ``dcio``'s ``test_bursts_fortran.py``, so what is left here is the
likelihood.

**What the record is.** Generated below, from Q, by a plain Gillespie walk with
a fixed seed. Its statistical fidelity does not matter and nothing here relies
on it: both likelihoods are evaluated on the *same* intervals, so the record is
an input, not a measurement. It is drawn from the mechanism only so that the
interval lengths exercise the same branches a real record would.

**What agreement to expect.** Not exact: ``hjclik.for`` stores intervals as
``real*4`` (``tint(ndmx,10)``), so each duration reaches the Fortran rounded to
about 6 parts in 10^8 while the C++ gets the full double. That is the floor,
and the measured difference sits on it. It is a difference *per interval*, so
that is what the tolerance is expressed in -- a fixed absolute tolerance would
hide a real per-interval discrepancy on short records, and a relative one
would hide it everywhere.

Measured, over records of 100 to 1950 intervals and several seeds:

=========================  =====================  ===============
configuration              |dlog10 L| / interval  tolerance here
=========================  =====================  ===============
one group, equilibrium     2.4 - 2.7 x 10^-8      2 x 10^-7
bursts at t_crit, CHS      1.6 x 10^-7            1 x 10^-6
=========================  =====================  ===============

The burst configuration's floor is about six times the other, or 4 x 10^-7 per
*burst*: it is where the CHS vectors and the ``real*4`` ``tcrit`` enter, and
each burst carries its own initial and final vector.

What that buys, measured on the same records:
:func:`test_the_comparison_would_notice_a_real_difference` shows a rate
constant wrong by one part in 10^3 lands 27 times outside the tolerance, a
dead time wrong by one part in 10^4 lands 48 times outside, and the wrong
initial vectors land five orders of magnitude outside. One part in 10^4 on a
rate constant is about 3 times the tolerance -- caught, but that is the edge of
what this resolves, and it is asserted where it can be rather than claimed
everywhere.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# pytest puts this directory on sys.path anyway; doing it explicitly means the
# module can also be imported directly, which is how the numbers in the
# docstring above were measured.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dcfortran import fortlik                                  # noqa: E402

from HJCFIT.likelihood import Log10Likelihood                  # noqa: E402

pytestmark = pytest.mark.skipif(
    not fortlik.available(),
    reason=f"the vendored Fortran is not built; run `{fortlik.BUILD_HINT}`")

#: Tolerance on |log10 L_fortran - log10 L_c++|, per interval in the record,
#: for a single group with equilibrium vectors. The real*4 storage of interval
#: durations puts a floor of 2.4-2.7e-8 per interval under any comparison;
#: this is that with a factor of 8 in hand.
PER_INTERVAL = 2e-7

#: The same for the burst configuration, whose floor is about six times higher
#: -- 1.6e-7 per interval, or 4e-7 per burst, since each burst carries its own
#: CHS initial and final vectors and t_crit reaches the Fortran as a real*4.
PER_INTERVAL_BURSTS = 1e-6

TRES = 25e-6            # dead time [s]
TCRIT = 3.5e-3          # critical shut time [s]
CONC = 30e-9            # agonist concentration [M]


# ---------------------------------------------------------------- mechanism

#: Scheme 1 of Colquhoun, Hatton & Hawkes (2003), with the "true" rates of
#: their Table 1 -- the mechanism the paper's Figs 3-5 were fitted with.
#: Rates in s^-1, association rates in M^-1 s^-1.
RATES = dict(
    beta_1a=50.0, alpha_1a=6000.0,          # ARa   <-> AR*a
    beta_1b=150.0, alpha_1b=50000.0,        # ARb   <-> AR*b
    beta_2=52000.0, alpha_2=2000.0,         # A2R   <-> A2R*
    k_m2a=1500.0, k_p2a=2.0e8,              # A2R   <-> ARb
    k_m2b=10000.0, k_p2b=4.0e8,             # A2R   <-> ARa
    k_m1a=1500.0, k_p1a=2.0e8,              # ARa   <-> R
    k_m1b=10000.0, k_p1b=4.0e8,             # ARb   <-> R
)

#: State order, open states first as both implementations require.
STATES = ["AR*a", "AR*b", "A2R*", "ARa", "ARb", "A2R", "R"]
KA, KB, KC = 3, 3, 1     # open, short-lived shut, long-lived shut


def q_matrix(conc=CONC, **override):
    """Q for scheme 1 at one concentration, states ordered open first.

    The diamond: R binds A at either of two sites to give ARa or ARb, each of
    which opens (AR*a, AR*b) and can bind again to give A2R, which opens to
    A2R*. Diagonals are set so the rows sum to zero.
    """
    r = dict(RATES, **override)
    q = np.zeros((7, 7))
    A, B, C, a, b, D, R = range(7)      # AR*a AR*b A2R* ARa ARb A2R R

    q[A, a], q[a, A] = r["alpha_1a"], r["beta_1a"]
    q[B, b], q[b, B] = r["alpha_1b"], r["beta_1b"]
    q[C, D], q[D, C] = r["alpha_2"], r["beta_2"]

    q[D, b], q[b, D] = r["k_m2a"], r["k_p2a"] * conc
    q[D, a], q[a, D] = r["k_m2b"], r["k_p2b"] * conc
    q[a, R], q[R, a] = r["k_m1a"], r["k_p1a"] * conc
    q[b, R], q[R, b] = r["k_m1b"], r["k_p1b"] * conc

    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


# ------------------------------------------------------------------- record

def simulate(Q, kA, nintervals, seed):
    """A record of apparent-class dwell times, by a Gillespie walk on Q.

    Consecutive states of the same class are merged, so what comes back
    alternates open and shut. Amplitude is 5 pA when open, 0 when shut.
    """
    rng = np.random.default_rng(seed)
    k = Q.shape[0]
    rates = -np.diag(Q)
    jump = Q / rates[:, None]
    np.fill_diagonal(jump, 0.0)
    jump = np.cumsum(jump, axis=1)

    state = k - 1                                   # start in the resting state
    tints, ampls = [], []
    total, current = 0.0, state < kA
    while len(tints) < nintervals:
        total += rng.exponential(1.0 / rates[state])
        state = int(np.searchsorted(jump[state], rng.random()))
        if (state < kA) != current:                 # class changed: emit
            tints.append(total)
            ampls.append(5.0 if current else 0.0)
            total, current = 0.0, state < kA
    return np.array(tints), np.array(ampls)


def impose_resolution(tints, ampls, tres):
    """Colquhoun & Sigworth: an unresolvable interval joins the one in progress.

    The rules ``RESHJC2.FOR`` applies, and the ones ``scalcs`` applies -- the
    two were checked against each other interval for interval in ``scalcs``'
    ``test_resolution_fortran.py``. Repeated here only so that this file needs
    nothing but numpy.
    """
    out_t, out_a = [tints[0]], [ampls[0]]
    for t, a in zip(tints[1:], ampls[1:]):
        if (a != out_a[-1]) and t >= tres:
            out_t.append(t)
            out_a.append(a)
        else:
            out_t[-1] += t
    return np.array(out_t), np.array(out_a)


def one_group(tints, ampls):
    """Trim to start and end on an opening, as both implementations require."""
    op = np.nonzero(ampls != 0.0)[0]
    lo, hi = op[0], op[-1]
    return tints[lo:hi + 1], ampls[lo:hi + 1]


def bursts(tints, ampls, tcrit):
    """Divide at shut times longer than tcrit, keeping openings at both ends.

    ``hjclik.for`` line 1009: a shut time strictly greater than tcrit "ends
    present group with prev opening".
    """
    groups, cur = [], []
    for t, a in zip(tints, ampls):
        if a == 0.0 and t > tcrit:
            if cur:
                groups.append(cur)
            cur = []
        else:
            cur.append((t, a))
    if cur:
        groups.append(cur)

    out = []
    for g in groups:
        while g and g[0][1] == 0.0:
            g = g[1:]
        while g and g[-1][1] == 0.0:
            g = g[:-1]
        if g:
            out.append(np.array([x[0] for x in g]))
    return out


def record(nintervals, seed, conc=CONC, tres=TRES):
    """A simulated record with the dead time imposed."""
    t, a = simulate(q_matrix(conc), KA, nintervals, seed)
    return impose_resolution(t, a, tres)


# --------------------------------------------------------------- comparison

def both(Q, tints, ampls, tres=TRES, tcrit=None, chs=False):
    """(Fortran, C++) log10 likelihoods of the same record.

    Without ``tcrit`` the record is one group and both use equilibrium
    vectors; with it, the Fortran divides the record itself and the C++ is
    given the bursts, and both use CHS vectors.
    """
    f, _ = fortlik.log10_likelihood(Q, KA, KB, KC, tints, ampls, CONC, tres,
                                    tcrit, chs)
    if tcrit is None:
        groups = [np.asarray(tints, float)]
    else:
        groups = [b for b in bursts(tints, ampls, tcrit) if b.size % 2 == 1]
    c = float(Log10Likelihood(groups, KA, tres, tcrit)(Q))
    return f, c


def tolerance(n, per_interval=PER_INTERVAL):
    return per_interval * n


# ------------------------------------------------------------------- tests

@pytest.fixture(scope="module")
def short_record():
    return one_group(*record(200, seed=7))


@pytest.fixture(scope="module")
def long_record():
    return one_group(*record(4000, seed=11))


def test_agrees_on_a_short_record(short_record):
    """A few hundred intervals, one group, equilibrium vectors."""
    t, a = short_record
    f, c = both(q_matrix(), t, a)
    assert abs(f - c) < tolerance(t.size), (
        f"Fortran {f:.12g}, C++ {c:.12g}, difference {f - c:+.4g} "
        f"over {t.size} intervals")


def test_agrees_on_a_long_record(long_record):
    """Thousands of intervals: a difference per interval would be visible."""
    t, a = long_record
    f, c = both(q_matrix(), t, a)
    assert abs(f - c) < tolerance(t.size), (
        f"Fortran {f:.12g}, C++ {c:.12g}, difference {f - c:+.4g} "
        f"over {t.size} intervals")


def test_agrees_with_bursts_and_chs_vectors(long_record):
    """The configuration CHH (2003) actually fitted.

    Divided into bursts at t_crit, with the CHS initial and final vectors of
    Colquhoun, Hawkes & Srodzinski (1996) rather than the equilibrium ones.
    This is the case the published results were produced in, and it exercises
    the initial and final vectors, which the single-group case does not.
    """
    t, a = long_record
    f, c = both(q_matrix(), t, a, tcrit=TCRIT, chs=True)
    nb = len([b for b in bursts(t, a, TCRIT) if b.size % 2 == 1])
    assert nb > 20, "the record should actually contain bursts"
    assert abs(f - c) < tolerance(t.size, PER_INTERVAL_BURSTS), (
        f"Fortran {f:.12g}, C++ {c:.12g}, difference {f - c:+.4g} "
        f"over {t.size} intervals in {nb} bursts")


@pytest.mark.parametrize("seed", [3, 5, 13])
def test_agrees_on_records_it_has_not_seen(seed):
    """Different records, in case one happened to be forgiving."""
    t, a = one_group(*record(600, seed=seed))
    f, c = both(q_matrix(), t, a)
    assert abs(f - c) < tolerance(t.size), (
        f"seed {seed}: Fortran {f:.12g}, C++ {c:.12g}, "
        f"difference {f - c:+.4g} over {t.size} intervals")


@pytest.mark.parametrize("conc,tres", [(30e-9, 25e-6), (100e-9, 25e-6),
                                       (30e-9, 50e-6)])
def test_agrees_at_other_concentrations_and_dead_times(conc, tres):
    """The missed-events correction does more work at a longer dead time."""
    t, a = one_group(*record(600, seed=17, conc=conc, tres=tres))
    f, c = both(q_matrix(conc), t, a, tres=tres)
    assert abs(f - c) < tolerance(t.size), (
        f"{conc * 1e9:.0f} nM, tres {tres * 1e6:.0f} us: "
        f"Fortran {f:.12g}, C++ {c:.12g}, difference {f - c:+.4g}")


def test_the_difference_stays_at_the_rounding_floor_as_records_grow():
    """It must not grow *per interval*.

    A genuine difference -- a different pdf, or a factor applied once per
    interval -- would show as a difference proportional to the record length,
    so a fixed tolerance on the total would hide it at short lengths and a
    fixed relative tolerance would hide it everywhere. What is asserted is the
    quantity that would be constant if something real were wrong.
    """
    per = {}
    for n in (200, 800, 3200):
        t, a = one_group(*record(n, seed=19))
        f, c = both(q_matrix(), t, a)
        per[t.size] = abs(f - c) / t.size
    for size, d in per.items():
        assert d < PER_INTERVAL, (
            f"{d:.3g} per interval at {size} intervals, over "
            f"{ {k: f'{v:.3g}' for k, v in per.items()} }")


def test_the_comparison_would_notice_a_real_difference():
    """The tolerance is only worth having if something could exceed it.

    Three things a broken implementation might plausibly do, each given to the
    C++ side only, must take the two likelihoods apart by more than the
    tolerance. Without this the tests above could pass by being loose.

    The measured margins on this record are 27x, 48x and about 25 000x, in the
    order below. A rate constant wrong by one part in 10^4 comes to roughly 3x
    the tolerance -- inside, but too near the edge to assert, and stated in the
    module docstring instead of claimed here.
    """
    t, a = one_group(*record(600, seed=3))
    Q = q_matrix()
    tol = tolerance(t.size)
    f, _ = fortlik.log10_likelihood(Q, KA, KB, KC, t, a, CONC, TRES)
    group = [np.asarray(t, float)]

    wrong_rate = float(
        Log10Likelihood(group, KA, TRES, None)(
            q_matrix(alpha_2=RATES["alpha_2"] * 1.001)))
    assert abs(f - wrong_rate) > 10 * tol, (
        "a rate constant wrong by one part in 10^3 would pass; the tolerance "
        "is too loose to mean anything")

    wrong_tres = float(Log10Likelihood(group, KA, TRES * 1.0001, None)(Q))
    assert abs(f - wrong_tres) > 10 * tol, (
        "a dead time wrong by one part in 10^4 would pass")


def test_the_comparison_would_notice_the_wrong_initial_vectors():
    """The initial and final vectors are the part a record cannot reveal.

    Two implementations can agree on every pdf and still differ in how a burst
    begins and ends -- CHS (1996) eqn (4) against the equilibrium vectors of
    eqn (3) -- and that is exactly the difference that matters for fitting
    bursts. Handing the C++ the equilibrium vectors while the Fortran uses CHS
    must be far outside the tolerance, or the burst test above proves nothing
    about the vectors.
    """
    t, a = one_group(*record(4000, seed=11))
    Q = q_matrix()
    groups = [b for b in bursts(t, a, TCRIT) if b.size % 2 == 1]
    f, _ = fortlik.log10_likelihood(Q, KA, KB, KC, t, a, CONC, TRES, TCRIT,
                                    chs=True)
    equilibrium = float(Log10Likelihood(groups, KA, TRES, None)(Q))
    assert abs(f - equilibrium) > 1000 * tolerance(t.size, PER_INTERVAL_BURSTS), (
        "the equilibrium vectors would pass for the CHS ones; the burst "
        "comparison is not testing the initial and final vectors")


def test_the_sign_and_base_conversion_is_right(short_record):
    """HJCLIK returns -ln L; Log10Likelihood returns +log10 L.

    Getting this wrong makes the Fortran look 2.3 times more curved than it
    is, which would look exactly like a real finding, so it is asserted rather
    than left in a comment.
    """
    t, a = short_record
    raw, _ = fortlik.log10_likelihood(q_matrix(), KA, KB, KC, t, a, CONC,
                                      TRES, raw=True)
    converted, _ = fortlik.log10_likelihood(q_matrix(), KA, KB, KC, t, a,
                                            CONC, TRES)
    assert raw < 0, "HJCLIK returns minus the log likelihood, of a likelihood > 1 here"
    np.testing.assert_allclose(converted, -raw / np.log(10), rtol=1e-15)


def test_the_vendored_sources_match_their_manifest():
    """vendor/ is what the manifest says it is.

    Not a check against DCFORTRAN itself -- that needs a checkout, and
    ``build.py --verify PATH`` does it -- but it catches an edit to a vendored
    file, which is the thing that would quietly turn "the original Fortran"
    into "our Fortran".
    """
    import subprocess
    build = Path(__file__).resolve().parent / "dcfortran" / "build.py"
    r = subprocess.run([sys.executable, str(build), "--verify"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "unmodified" in r.stdout
