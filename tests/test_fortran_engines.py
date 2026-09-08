"""The substituted numerical kernels against the ones DCFORTRAN actually used.

Eight of the sources ``hjclik.for`` needs are not DCPROGS' own work, and carry
a third party's copyright, so they are not vendored (see
``tests/dcfortran/README.md``):

* five are NAG Mark 2 (1972) -- ``F02AGF`` and the ``F01AKF`` / ``F01APF`` /
  ``F02AQF`` / ``A02ACF`` chain it calls, which is every eigendecomposition of
  a Q matrix, via ``QMAT5.FOR``;
* ``LUDCMPD.FOR`` and ``LUBKSBD.FOR`` say in their own comments that they are
  ``real*8`` versions of the Numerical Recipes routines, and ``MATINV2.FOR``
  calls them for every matrix inversion on the likelihood path;
* ``determ2.for`` holds DCPROGS' own ``DETERM2`` *and* an appended copy of the
  Numerical Recipes ``LUDCMP``, which cannot be separated without editing a
  vendored file.

``tests/dcfortran/replacements/`` supplies working equivalents: netlib EISPACK
for the eigen chain, ordinary partial-pivoting Gaussian elimination for the
LU, and ``DETERM2`` written out again. That is the ``free`` engine, and it is
what builds by default and what runs in CI.

**Substituting the numerical core of the program you are validating against is
not a free move.** If the replacement computed something slightly different,
every agreement measured in ``test_fortran_likelihood.py`` would be an
agreement with a modified program, and would not say what it claims to say.
So this file builds both and runs them against each other.

It needs a DCFORTRAN checkout, since the ``original`` engine cannot be built
without one, and skips otherwise::

    python tests/dcfortran/build.py --engine original --dcfortran ../DCFORTRAN
    pytest tests/test_fortran_engines.py -v

Measured, over the same records ``test_fortran_likelihood.py`` uses:

===============================  =========================
comparison                       |dlog10 L| per interval
===============================  =========================
free engine vs original          1.8 x 10^-14
either engine vs the C++         2.4 - 2.7 x 10^-8
===============================  =========================

The substitution is therefore about a million times smaller than the
difference it sits inside, which is itself only the ``real*4`` storage of the
intervals. NAG's chain and EISPACK are both Fortran translations of the same
ALGOL from Wilkinson & Reinsch (1971) -- NAG's headers name the procedures,
DIRHES, DIRTRANS, HQR2, CDIV -- so this is close to the best case, and the
measurement says so rather than the reasoning.

One difference is real and is not hidden: the Numerical Recipes LU chooses its
pivot by implicit scaling, dividing each candidate by the largest element of
its row, where the replacement compares the candidates themselves. On a badly
scaled matrix the two can pivot differently. The measurement above is what
that is worth here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from dcfortran import fortlik                                  # noqa: E402

import test_fortran_likelihood as base                         # noqa: E402

pytestmark = pytest.mark.skipif(
    not (fortlik.available("free") and fortlik.available("original")),
    reason="needs both engines built; the original one needs a DCFORTRAN "
           "checkout: build.py --engine original --dcfortran PATH")

#: Tolerance on |log10 L_free - log10 L_original|, per interval. The measured
#: difference is 1.8e-14; this is that with room for a different compiler, and
#: is still three orders of magnitude below the real*4 floor that separates
#: either engine from the C++.
PER_INTERVAL = 1e-11


def engines(Q, tints, ampls, tres=base.TRES, tcrit=None, chs=False):
    """(free, original) log10 likelihoods of the same record."""
    out = []
    for engine in ("free", "original"):
        v, _ = fortlik.log10_likelihood(
            Q, base.KA, base.KB, base.KC, tints, ampls, base.CONC, tres,
            tcrit, chs, engine=engine)
        out.append(v)
    return out


@pytest.mark.parametrize("n,seed", [(200, 7), (600, 3), (600, 13), (3200, 19)])
def test_the_replacements_compute_the_same_likelihood(n, seed):
    """Records of a few hundred to a few thousand intervals, one group."""
    t, a = base.one_group(*base.record(n, seed=seed))
    free, orig = engines(base.q_matrix(), t, a)
    assert abs(free - orig) < PER_INTERVAL * t.size, (
        f"free {free:.14g}, original {orig:.14g}, "
        f"difference {free - orig:+.4g} over {t.size} intervals")


def test_the_same_with_bursts_and_chs_vectors():
    """The configuration CHH (2003) fitted, which uses the initial vectors."""
    t, a = base.one_group(*base.record(4000, seed=11))
    free, orig = engines(base.q_matrix(), t, a, tcrit=base.TCRIT, chs=True)
    assert abs(free - orig) < PER_INTERVAL * t.size, (
        f"free {free:.14g}, original {orig:.14g}")


@pytest.mark.parametrize("conc,tres", [(100e-9, 25e-6), (30e-9, 50e-6)])
def test_the_same_at_other_concentrations_and_dead_times(conc, tres):
    t, a = base.one_group(*base.record(600, seed=17, conc=conc, tres=tres))
    free, orig = engines(base.q_matrix(conc), t, a, tres=tres)
    assert abs(free - orig) < PER_INTERVAL * t.size, (
        f"{conc * 1e9:.0f} nM, tres {tres * 1e6:.0f} us: "
        f"free {free:.14g}, original {orig:.14g}")


def test_the_substitution_is_far_smaller_than_the_difference_it_sits_inside():
    """The claim this file exists to support, stated as a number.

    If the two engines were as far apart as either is from the C++, the
    comparison in ``test_fortran_likelihood.py`` would be measuring the
    replacement rather than the original, and the tolerances there would have
    been tuned to a program DCFORTRAN never ran.
    """
    t, a = base.one_group(*base.record(3200, seed=19))
    free, orig = engines(base.q_matrix(), t, a)
    between = abs(free - orig) / t.size
    to_cpp = abs(orig - base.both(base.q_matrix(), t, a)[1]) / t.size
    assert between < to_cpp / 100, (
        f"the engines differ by {between:.3g} per interval and the C++ by "
        f"{to_cpp:.3g}; the substitution is no longer negligible against the "
        f"comparison it is inside")


def test_the_comparison_would_notice_a_broken_replacement():
    """Both engines are being run, and a real difference would show.

    A tolerance this tight is worth nothing if the two programs were somehow
    the same binary, or if the numbers were not really being recomputed. A
    rate constant moved by one part in 10^6 -- far below anything a kernel
    substitution should cause -- must exceed it.
    """
    t, a = base.one_group(*base.record(600, seed=3))
    free, _ = engines(base.q_matrix(), t, a)
    moved, _ = fortlik.log10_likelihood(
        base.q_matrix(alpha_2=base.RATES["alpha_2"] * 1.000001),
        base.KA, base.KB, base.KC, t, a, base.CONC, base.TRES, engine="free")
    assert abs(free - moved) > 10 * PER_INTERVAL * t.size, (
        "a rate constant wrong by one part in 10^6 would pass; this "
        "comparison cannot see anything")


def test_the_two_engines_are_not_the_same_program():
    """Cheap, and it would have caught a build that silently reused objects."""
    free = fortlik.exe_for("free")
    orig = fortlik.exe_for("original")
    assert free != orig and free.exists() and orig.exists()
    assert free.read_bytes() != orig.read_bytes(), (
        "the two engines produced identical binaries, so one of them was not "
        "rebuilt and the comparison is vacuous")
