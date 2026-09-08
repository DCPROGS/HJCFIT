# What is in here, and why

Eight of the sources `hjclik.for` needs are not DCPROGS' own work and carry a
third party's copyright, so they are not vendored in `../vendor/`. This
directory holds working equivalents. Building these is the `free` engine, which
is the default and the one CI runs; `build.py --engine original --dcfortran
PATH` builds the real ones out of a DCFORTRAN checkout instead.

`../../test_fortran_engines.py` runs both and measures the difference. It is
**1.8 × 10⁻¹⁴ per interval** in log₁₀ likelihood — double-precision rounding,
about a millionth of the `real*4` floor that separates either engine from the
C++.

## What was replaced

| not vendored | why | replaced by |
|---|---|---|
| `F02AGF.FOR` | NAG Mark 2, 1972 | `f02agf.f` over EISPACK |
| `F01AKF.FOR` | NAG Mark 2 — its header says `DIRHES` | `eispack/elmhes.f` |
| `F01APF.FOR` | NAG Mark 2 — `DIRTRANS` | `eispack/eltran.f` |
| `F02AQF.FOR` | NAG Mark 2 — `HQR2` | `eispack/hqr2.f` |
| `A02ACF.FOR` | NAG Mark 2 — `CDIV` | `eispack/cdiv.f` |
| `LUDCMPD.FOR` | Numerical Recipes `ludcmp`, `real*8` | `lu.f` |
| `LUBKSBD.FOR` | Numerical Recipes `lubksb`, `real*8` | `lu.f` |
| `determ2.for` | DCPROGS' `DETERM2` **plus** an appended NR `LUDCMP` | `determ2.f` |

The NAG files say so themselves — `C    MARK 2 RELEASE. NAG COPYRIGHT 1972` —
and the Numerical Recipes ones say so in DCPROGS' own comments (*"Double
precision version of Numerical Recipes routine"*). `determ2.for` is the awkward
one: `DETERM2` is DCPROGS', but it shares a file with a second copy of the NR
`LUDCMP`, and separating them would mean editing a vendored file, which is
exactly what `vendor/` exists not to do. So the file goes and `DETERM2` is
written out again here, unchanged in behaviour.

Neither group is incidental. `F02AGF` is reached through `QMAT5.FOR` from
`HJCEXACT`, `HJCASYMP`, `GFUNCA`, `GFUNCF` and `HJCMEAN` — every
eigendecomposition of a Q matrix. `LUDCMPD` is reached through `MATINV2.FOR`
from essentially every matrix inversion on the likelihood path. Stubbing them
was never an option.

## `eispack/` — verbatim, public domain

Downloaded from <https://www.netlib.org/eispack/>, unmodified:

| file | bytes | sha256 |
|---|---|---|
| `cdiv.f` | 390 | `3c2ec3c8f8bfa239175bf673d3085f869b25242cb4b7e9ed8f53e860abb973bd` |
| `elmhes.f` | 2727 | `699c7740e93097b28af396fea60a4563f30ae9a9c812e466ef79214f13f8a972` |
| `eltran.f` | 2263 | `e999eabada3c5b979c65dddfee42c98c83a8e4803b0e861063b9294e946e4576` |
| `hqr2.f` | 14212 | `4c56daecd3e92202a424961de53f799c31259032526b4c2086540ef58bc9dd97` |

EISPACK came out of Argonne National Laboratory and is distributed by netlib
without restriction.

This is the closest substitution available rather than merely a convenient one.
NAG's Mark 2 chain and EISPACK are both Fortran translations of the *same* ALGOL
procedures from Wilkinson & Reinsch, *Handbook for Automatic Computation* vol II
(1971) — NAG's own comment lines name them, `DIRHES`, `DIRTRANS`, `HQR2`,
`CDIV`, which are `elmhes`, `eltran`, `hqr2` and `cdiv` here. Same algorithm,
same era, different translator.

## `f02agf.f`, `lu.f`, `determ2.f` — ours

Written for this, and GPLv3 with the rest of HJCFIT.

- **`f02agf.f`** keeps NAG's interface, because `QMAT5.FOR` calls it, and
  computes it with the EISPACK chain — `elmhes`, `eltran`, `hqr2`, with
  `low = 1` and `igh = n`, i.e. no balancing, which is what NAG's `F02AGF`
  does too. It unpacks `hqr2`'s packed eigenvectors into the separate real and
  imaginary parts NAG's interface promises.

  Scaling does not matter here and the file says why: `QMAT5` rescales every
  column so its top element is 1, then forms the spectral matrices as
  `A(m) = EM(:,m) EN(m,:)` with `EN = inv(EM)`. Scaling a column of `EM` scales
  the matching row of its inverse by the reciprocal, so it cancels exactly.

- **`lu.f`** is Gaussian elimination with partial pivoting — the right-looking
  form, what LAPACK calls `dgetf2` — with the interface `MATINV2.FOR` and
  `DETERM2` use.

  **One deliberate difference.** The Numerical Recipes routine pivots by
  *implicit scaling*: it divides each candidate by the largest element of its
  row before comparing. This one compares the candidates themselves. Both are
  backward stable and both give an exact LU of a permutation of the matrix, but
  on a badly scaled matrix they can choose different pivots and round
  differently in the last bits. That is the one place the two engines can
  disagree, and the 1.8 × 10⁻¹⁴ above is what it costs on these matrices.

- **`determ2.f`** is `DETERM2` again: LU, then the product of the diagonal of
  `U`, taking out a factor of 10¹⁰ into `ndscale` whenever the next
  multiplication would leave the range of a `real*8`. `DETWA.FOR` and
  `DETWF.FOR` use it for the determinant of *W(s)*, whose roots are the
  asymptotic time constants, and those determinants run over many orders of
  magnitude, so the scaling is not decoration.

## The thing to be careful about

Substituting the numerical core of the program you are validating against is
not a free move. If any of this computed something slightly different, the
agreement measured in `../../test_fortran_likelihood.py` would be agreement
with a modified program and would not mean what it says.

That is why `../../test_fortran_engines.py` exists, why it runs both engines
rather than reasoning about them, and why it asserts that the gap between them
is at least a hundred times smaller than the gap it sits inside.
