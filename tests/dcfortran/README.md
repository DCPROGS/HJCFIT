# The 2003 Fortran HJCFIT likelihood, vendored

`hjclik.for` and its callees are the likelihood that produced every published
HJCFIT result — Colquhoun, Hatton & Hawkes (2003), *J Physiol* **547**:699–728
among them. HJCFIT's C++ `Log10Likelihood` reimplements the same theory
(Hawkes, Jalali & Colquhoun's exact treatment of missed events). This directory
builds the Fortran so the two can be evaluated on the same Q matrix and the
same record, which is what `tests/test_fortran_likelihood.py` does.

It is worth having because it is a different *kind* of evidence. The rest of
the suite checks the C++ against closed forms, against limiting cases, and
against itself computed another way. None of that can catch an error that both
implementations of a formula would make, or a misreading of the papers that
looks self-consistent. An independent program written twenty years earlier by
the authors of the theory can.

## Building it

Not built by default, and the tests skip if it is missing:

```bash
python tests/dcfortran/build.py
pytest tests/test_fortran_likelihood.py -v
```

It needs a Fortran compiler, found through `$FC` or on `PATH`:

| platform | |
|---|---|
| Debian/Ubuntu | `apt install gfortran` |
| macOS | `brew install gcc` |
| Windows | `conda install -c conda-forge m2w64-gcc-fortran` |

On Windows the compiler's own `bin` directory has to reach it in Windows form —
from Git Bash a POSIX-style `PATH` breaks gfortran's DLL lookup and it
segfaults on hello-world. `build.py` arranges that itself; it is only a problem
if you invoke gfortran by hand.

Everything the build writes goes under `build/`, which is ignored by git.
`vendor/` is only ever read.

## What is vendored, and that it is unmodified

`vendor/` holds 60 files copied byte for byte out of
[DCPROGS/DCFORTRAN](https://github.com/DCPROGS/DCFORTRAN). `MANIFEST.txt`
records, for each one, where in `Fort90/` it came from and its SHA-256.

```bash
python tests/dcfortran/build.py --verify                     # against the manifest
python tests/dcfortran/build.py --verify /path/to/DCFORTRAN  # against the original
python tests/dcfortran/build.py --vendor  /path/to/DCFORTRAN  # re-copy
```

The first of those runs as a test
(`test_the_vendored_sources_match_their_manifest`), so an edit to a vendored
file fails the suite rather than quietly turning "the original Fortran" into
"our Fortran". The second is the one that actually proves the claim, and needs
a checkout.

The code is deliberately **not** cleaned up. It is fixed-form Fortran 77 with
COMMON blocks, `real*4` storage, computed GOTOs and DOS-era conventions, and
every one of those is part of what is being compared. Tidying it would produce
something that agrees with the C++ because it had been made to; leaving it
alone is the whole point.

## The three changes, and why each is safe

`build.py` copies `vendor/` into `build/src/` and patches the copies. Every
change is in `PATCHES` and `STRIP_BYTES` there, with its reason, and the build
prints them with counts on every run. They are:

**1. `call TIMER(` → `call DCTIMER(`** (12 occurrences)
`TIMER` is both a COMMON block and a subroutine name in this code. Lahey
Fortran allowed that; gfortran does not. Renaming the *calls* leaves
`common/timer/` untouched, and the routine only reads a clock for a debug
printout — `driver/stubs.f` supplies a `DCTIMER` that returns zero. It cannot
reach a likelihood.

**2. `int4(` → `int(`** (5 occurrences)
`INT4` is an external truncation function from the Spindrift utility library,
which is not part of DCFORTRAN. The intrinsic `int()` does the same thing and,
being generic, accepts both the `real*8` and the `real*4` call sites — an
external `INT4` could only have one argument type and would misread the other.
It is used for a bisection step count and for a sign.

**3. Stray control bytes** (1 occurrence)
A `0x1B` after the last statement of `DETWA.FOR`, left by a DOS editor. Bytes
outside any statement; `0x1A` and `0x0C` are filtered too, and none occurs.

Nothing else is altered. A change that affected arithmetic would have to be
declared in the same place, and there is none.

## What is *not* vendored: `driver/`

Three files here are ours, and are not copies of anything:

- **`likdrv.f`** — a non-interactive main program. It reads a Q matrix, a
  record and the settings from a file, calls `HJCLIK` once, and prints the
  number. It exists because the original entry point is a full interactive DOS
  application with prompts, graphics, a simplex and topology machinery, none of
  which is being tested. The COMMON blocks it fills are the ones `HJCLIK`
  reads; the comments in it record the two things that are easy to get wrong
  (`km` is Q's *declared* dimension, 100, not the state count; and `tcrit = 0`
  means "split at every shut time", not "do not split").
- **`qset_stub.f`** — replaces `QSET_HJC`, which builds Q from a parameter
  vector through HJCFIT's constraint, micro-reversibility and concentration
  machinery. This one copies a Q matrix supplied directly, and sets the
  diagonals the way `SETDIAG` does. So what the comparison tests is *the
  likelihood given a generator*, with neither side's parameterisation in it.
- **`stubs.f`** — the DOS/Lahey library routines that are linked but never
  reached on this path: a clock, a bell, a keyboard poll, an underflow mode
  switch, a random number generator used only by the reset-and-perturb path,
  and three string routines from Spindrift used only for printed titles.
  `NBLANK` is implemented properly rather than stubbed, since a wrong value
  could index outside a string.

`QSET_HJC.FOR` and `RANPERT.FOR` are vendored anyway, so that what was replaced
can be read, but they are not compiled.

## Two things found while building it

Recorded because they cost time and would cost it again.

**`hjcasym1.for` calls `ROOT_FB` with nine arguments where `ROOT_FB.FOR`
declares ten.** The missing `nerr` writes over the caller's frame; it silently
zeroed `kAm`, and the likelihood then indexed a zero-length array. The build
uses `HJCASYMP.FOR` instead — both carry the same "Modified 03/13/03" and the
same signature, and only the `ROOT_FB` call differs. Version-matching inside
DCFORTRAN cannot be assumed, and `build.py --check` (bounds checking) is what
found it.

**`HJCLIK` returns minus the *natural* log likelihood**, rescaled by
`nscal * 230.2585093` = ln(10¹⁰⁰) for underflow (`hjclik.for` lines 1306, 1318,
1386). The C++ returns `+log10 L`. Comparing them without both corrections
makes the Fortran look 2.3 times more curved than it is —
`test_the_sign_and_base_conversion_is_right` asserts the conversion rather than
leaving it in a comment.

## Licence

HJCFIT is GPLv3. **DCFORTRAN carries no licence file at all** — it is published
by DCPROGS, the same organisation, and is the direct ancestor of this project,
but that is not the same as a grant. Worth settling: a `LICENSE` in DCFORTRAN,
or a note from its authors, would put this directory on the same footing as the
rest of the repository. Until then, treat these 60 files as vendored on the
strength of common ownership rather than on an explicit licence, and do not
redistribute them separately.

Nothing here ships in a wheel: `pyproject.toml` sets `wheel.packages = []` and
builds with `tests = OFF`, so `tests/` is in the source repository and the
sdist only.
