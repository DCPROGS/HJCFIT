"""Build the 2003 Fortran HJCFIT likelihood as a callable program.

The Fortran in ``vendor/`` is the likelihood that produced every published
HJCFIT result, Colquhoun, Hatton & Hawkes (2003) among them. It is vendored
**unmodified**: ``vendor/`` is byte-for-byte what is in ``DCPROGS/DCFORTRAN``,
and :func:`verify` will prove that against a checkout. The changes gfortran
needs are applied to *copies* under ``build/src/``, and every one of them is
in :data:`PATCHES` with its reason, so the difference from the original is a
short list rather than a diff nobody reads.

    python tests/dcfortran/build.py                    # patch, compile, link
    python tests/dcfortran/build.py --clean
    python tests/dcfortran/build.py --verify PATH      # vendor/ vs DCFORTRAN
    python tests/dcfortran/build.py --vendor PATH      # re-copy from DCFORTRAN

Needs a Fortran compiler: ``apt install gfortran``, ``brew install gcc``, or
``conda install -c conda-forge m2w64-gcc-fortran`` on Windows. It is found
through ``$FC`` or on ``PATH``. On Windows the compiler's own ``bin``
directory has to be on ``PATH`` in *Windows* form for its DLL lookup to work
-- from Git Bash a POSIX-style ``PATH`` makes it segfault on hello-world --
which is what :func:`_env` arranges.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VENDOR = HERE / "vendor"
BUILD = HERE / "build"
#: The patched copies live *inside* ``build/``, which is ignored by git, so
#: that ``vendor/`` cannot be edited in place by accident and "the vendored
#: code is unmodified" stays true without anyone having to remember it.
SRC = BUILD / "src"
MANIFEST = VENDOR / "MANIFEST.txt"

FLAGS = ["-ffixed-form", "-ffixed-line-length-none", "-fno-automatic",
         "-std=legacy", "-w", "-O2"]

#: Added by ``--debug``: symbols and a backtrace, for locating a fault.
DEBUG_FLAGS = ["-g", "-fbacktrace", "-O0"]

#: Added by ``--check``. Old Fortran sometimes indexes past a declared bound
#: deliberately, so a report here is a lead rather than a verdict -- but it is
#: how the ``hjcasym1.for`` / ``ROOT_FB`` argument mismatch below was found.
CHECK_FLAGS = ["-fcheck=bounds"]

#: Fortran sources, as paths relative to ``DCFORTRAN/Fort90``. Order is the
#: order they were added: the likelihood and its immediate callees first, then
#: three rounds of resolving undefined symbols from the linker.
SOURCES = [
    "HJCFIT/hjclik.for",
    "HJCFIT/HJCEXACT.FOR",
    # HJCASYMP.FOR, not hjcasym1.for: the latter calls ROOT_FB with nine
    # arguments where ROOT_FB.FOR declares ten, so the missing nerr writes
    # over the caller's frame -- it silently zeroed kAm, and the likelihood
    # then indexed a zero-length array. Both carry the same "Modified
    # 03/13/03" and the same HJCASYMP signature; only the ROOT_FB call
    # differs. Version-matching inside DCFORTRAN cannot be assumed.
    "HJCFIT/HJCASYMP.FOR",
    "HJCFIT/QSET_HJC.FOR",
    "HJCFIT/eqoc_hjc.for",
    "HJCFIT/HJCMEAN.FOR",
    "HJCFIT/RANPERT.FOR",
    "CALC/rac3d.for",
    "CALC/SUBMAT.FOR",
    "CALC/MATMUL.FOR",
    "CALC/VECMUL.FOR",
    "CALC/PHIO1.FOR",
    "CALC/MATSCL2.FOR",
    # resolved from the first link's undefined symbols
    "HJCFIT/ACONVD.FOR",
    "HJCFIT/CEGAF.FOR",
    "HJCFIT/DENOMA.FOR",
    "HJCFIT/DENOMF.FOR",
    "HJCFIT/DETWA.FOR",
    "HJCFIT/DETWF.FOR",
    "HJCFIT/FINDGRP.FOR",
    "HJCFIT/GFUNCA.FOR",
    "HJCFIT/GFUNCF.FOR",
    "HJCFIT/HJCPHI.FOR",
    "HJCFIT/ROOT_FB.FOR",
    "BISECD.FOR",
    "DEXP1.FOR",
    "CALC/ATYPD.FOR",
    "CALC/ATYPD3.FOR",
    "CALC/DATYP.FOR",
    "CALC/MATADD.FOR",
    "CALC/MATIMIN.FOR",
    "CALC/MATINV.FOR",
    "CALC/MATINV2.FOR",
    "CALC/MATMUL3.FOR",
    "CALC/mattran1.for",
    "CALC/PDFOPEN.FOR",
    "CALC/PDFOUTD.FOR",
    "CALC/PDFOUTS.FOR",
    "CALC/old/PDFSHUT.FOR",
    "CALC/QMAT5.FOR",
    "CALC/submat3d.for",
    "CVFIT/GETGROUP.FOR",
    # resolved from the second link
    "ARRAYD.FOR",
    "SORT2D.FOR",
    "SORTD.FOR",
    "CALC/determ2.for",
    "CALC/F02AGF.FOR",
    "CALC/GMAT1.FOR",
    "CALC/LUBKSBD.FOR",
    "CALC/LUDCMPD.FOR",
    "CALC/MATMUL2.FOR",
    "CALC/MATSCL3.FOR",
    "CALC/MINVD.FOR",
    # the eigenvalue chain F02AGF calls: NAG-compatible routines, local
    # implementations rather than the commercial library
    "CALC/F01AKF.FOR",
    "CALC/F01APF.FOR",
    "CALC/F02AQF.FOR",
    "CALC/A02ACF.FOR",
    # needed by HJCASYMP.FOR
    "HJCFIT/BISECHJC.FOR",
    "HJCFIT/CHECKRW.FOR",
    "HJCFIT/EQOC_RED.FOR",
]

#: Every textual change made to the vendored code, as ``(old, new, why)``.
#: All three are compiler-compatibility changes; none of them touches
#: arithmetic. Anything that changed behaviour would have to say so here.
PATCHES = [
    ("call TIMER(", "call DCTIMER(",
     "TIMER is both a COMMON block and a subroutine. Lahey allowed it, "
     "gfortran does not. Renaming the calls leaves common/timer/ alone, and "
     "the routine only reads a clock for debug printout, so this cannot "
     "change a likelihood."),
    ("call timer(", "call DCTIMER(", "the same, lower case"),
    ("int4(", "int(",
     "INT4 is an external truncation function from the Spindrift utility "
     "library, which is not in DCFORTRAN. The intrinsic int() does the same "
     "thing and, being generic, accepts both the real*8 and the real*4 call "
     "sites -- an external INT4 could only have one argument type and would "
     "misread the other. Used for a bisection step count and a sign."),
    ("INT4(", "int(", "the same, upper case"),
]

#: Bytes that are not valid in fixed-form Fortran but appear in some files --
#: a stray ESC (0x1B) after the last statement of DETWA.FOR, DOS end-of-file
#: marks, a form feed. Dropping them changes no code: they lie outside any
#: statement. Counted and reported like a patch, because a silent byte-level
#: edit is exactly the kind of thing this file exists to make visible.
STRIP_BYTES = bytes([0x1A, 0x1B, 0x0C])

#: Our own sources, in ``driver/``: a non-interactive main program, a
#: QSET_HJC that takes a Q matrix directly instead of building one from
#: HJCFIT's topology machinery, and stubs for the DOS library routines. These
#: are not copies of anything and are not vendored code.
OWN = ["likdrv.f", "qset_stub.f", "stubs.f"]

#: Vendored sources that ``OWN`` replaces. They are vendored anyway, so that
#: what was replaced can be read, but they are not compiled.
REPLACED = {"QSET_HJC.FOR", "RANPERT.FOR"}

EXE = "likdrv.exe" if platform.system() == "Windows" else "likdrv"


def compiler():
    """The Fortran compiler: ``$FC``, or ``gfortran`` on ``PATH``."""
    exe = os.environ.get("FC") or shutil.which("gfortran")
    if not exe:
        raise SystemExit(
            "no Fortran compiler found. Set $FC, or put gfortran on PATH:\n"
            "  apt install gfortran / brew install gcc\n"
            "  conda install -c conda-forge m2w64-gcc-fortran   (Windows)")
    found = shutil.which(exe)
    if not found:
        raise SystemExit(f"$FC is set to {exe!r}, which is not executable")
    return Path(found).resolve()


def _env(fc):
    """Environment with the compiler's own bin directory first on PATH.

    On Windows gfortran finds its DLLs relative to PATH, and a POSIX-style
    PATH inherited from Git Bash makes it fail in ways that look like a
    miscompilation. Harmless elsewhere.
    """
    env = dict(os.environ)
    env["PATH"] = f"{fc.parent}{os.pathsep}{env.get('PATH', '')}"
    return env


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def vendor(dcfortran):
    """Copy the sources out of a DCFORTRAN checkout and write the manifest."""
    root = Path(dcfortran)
    root = root / "Fort90" if (root / "Fort90").is_dir() else root
    VENDOR.mkdir(parents=True, exist_ok=True)
    rows = []
    for rel in SOURCES:
        src = root / rel
        if not src.exists():
            raise SystemExit(f"missing source: {src}")
        dst = VENDOR / Path(rel).name
        shutil.copyfile(src, dst)
        rows.append((Path(rel).name, rel, sha256(dst)))
    MANIFEST.write_text(
        "# Vendored from DCPROGS/DCFORTRAN, unmodified.\n"
        "# name  path-within-Fort90  sha256\n"
        + "".join(f"{n}  {r}  {h}\n" for n, r, h in sorted(rows)))
    print(f"vendored {len(rows)} files into {VENDOR}")
    return 0


def verify(dcfortran=None):
    """Check ``vendor/`` against the manifest, and against a checkout if given.

    The manifest check is the one that runs anywhere; the checkout check is
    the one that means something, and is how the claim that these files are
    unmodified can be tested by a reader rather than believed.
    """
    if not MANIFEST.exists():
        raise SystemExit(f"no manifest at {MANIFEST}")
    want = {}
    for line in MANIFEST.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        name, rel, digest = line.split()
        want[name] = (rel, digest)

    bad = []
    for name, (rel, digest) in sorted(want.items()):
        got = VENDOR / name
        if not got.exists():
            bad.append(f"{name}: missing from vendor/")
        elif sha256(got) != digest:
            bad.append(f"{name}: differs from the manifest")
    extra = {p.name for p in VENDOR.glob("*")} - set(want) - {MANIFEST.name}
    for name in sorted(extra):
        bad.append(f"{name}: in vendor/ but not in the manifest")

    if dcfortran:
        root = Path(dcfortran)
        root = root / "Fort90" if (root / "Fort90").is_dir() else root
        for name, (rel, digest) in sorted(want.items()):
            up = root / rel
            if not up.exists():
                bad.append(f"{name}: not found upstream at {rel}")
            elif sha256(up) != digest:
                bad.append(f"{name}: differs from upstream {rel}")
        print(f"checked {len(want)} files against {root}")
    else:
        print(f"checked {len(want)} files against the manifest "
              f"(pass --verify PATH to check against DCFORTRAN itself)")

    for line in bad:
        print(f"  {line}")
    print("unmodified" if not bad else f"{len(bad)} discrepancies")
    return 1 if bad else 0


def prepare():
    """Copy ``vendor/`` into ``build/src/`` and patch the copies."""
    SRC.mkdir(parents=True, exist_ok=True)
    copied = []
    for rel in SOURCES:
        source = VENDOR / Path(rel).name
        if not source.exists():
            raise SystemExit(
                f"missing vendored source: {source}\n"
                f"run `python {Path(__file__).name} --vendor PATH_TO_DCFORTRAN`")
        target = SRC / Path(rel).name
        shutil.copyfile(source, target)
        copied.append(target)

    counts = {old: 0 for old, _, _ in PATCHES}
    stripped = 0
    for target in copied:
        text = target.read_text(encoding="latin-1")
        for old, new, _ in PATCHES:
            n = text.count(old)
            if n:
                counts[old] += n
                text = text.replace(old, new)
        clean = "".join(c for c in text if ord(c) not in STRIP_BYTES)
        stripped += len(text) - len(clean)
        target.write_text(clean, encoding="latin-1")

    for old, new, why in PATCHES:
        print(f"  {counts[old]:3d} x {old!r} -> {new!r}")
        print(f"      {why}")
    print(f"  {stripped:3d} stray control bytes removed")
    return copied


def compile_all(files, fc, flags):
    """Compile each source; return ``(objects, failures)``."""
    objects, failures = [], []
    for f in files:
        obj = BUILD / (f.stem + ".o")
        r = subprocess.run([str(fc), "-c", *flags, str(f), "-o", str(obj)],
                           capture_output=True, text=True, env=_env(fc))
        if r.returncode == 0:
            objects.append(obj)
        else:
            failures.append((f.name, r.stderr))
    return objects, failures


def link(objects, fc):
    exe = BUILD / EXE
    extra = []
    if platform.system() == "Windows":
        # The program is run from Python without the compiler's environment
        # activated, and a dynamically linked one dies with 0xC0000135
        # (STATUS_DLL_NOT_FOUND) looking for libgfortran. macOS has no static
        # libc at all, and on Linux the runtime is there anyway.
        extra = ["-static", "-static-libgfortran", "-static-libgcc"]
    r = subprocess.run([str(fc), *[str(o) for o in objects], *extra,
                        "-o", str(exe)],
                       capture_output=True, text=True, env=_env(fc))
    return exe, r


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--clean", action="store_true",
                   help="remove build/; vendor/ is never touched")
    p.add_argument("--verify", nargs="?", const="", metavar="DCFORTRAN",
                   help="check vendor/ against the manifest, and against a "
                        "DCFORTRAN checkout if one is given")
    p.add_argument("--vendor", metavar="DCFORTRAN",
                   help="re-copy the sources from a DCFORTRAN checkout")
    p.add_argument("--check", action="store_true",
                   help="array bounds checking (implies --debug)")
    p.add_argument("--debug", action="store_true",
                   help="symbols and a backtrace, for locating a fault")
    p.add_argument("--show-errors", type=int, default=12,
                   help="lines of compiler output to show per failure")
    a = p.parse_args()

    if a.clean:
        if BUILD.exists():
            shutil.rmtree(BUILD)
        print("cleaned")
        return 0
    if a.vendor:
        return vendor(a.vendor)
    if a.verify is not None:
        return verify(a.verify or None)

    fc = compiler()
    r = subprocess.run([str(fc), "--version"], capture_output=True, text=True,
                       env=_env(fc))
    print(r.stdout.splitlines()[0] if r.stdout else str(fc))
    print()

    BUILD.mkdir(parents=True, exist_ok=True)
    print("patching copies of the vendored sources")
    files = [f for f in prepare() if f.name not in REPLACED]
    files += [HERE / "driver" / n for n in OWN]

    flags = list(FLAGS)
    if a.debug or a.check:
        flags += DEBUG_FLAGS
    if a.check:
        flags += CHECK_FLAGS
    print(f"\ncompiling {len(files)} files"
          + (" (debug)" if a.debug or a.check else ""))
    objects, failures = compile_all(files, fc, flags)
    print(f"  {len(objects)} of {len(files)} compiled")
    for name, err in failures:
        print(f"\n--- {name}")
        lines = [ln for ln in err.splitlines() if ln.strip()]
        for ln in lines[:a.show_errors]:
            print(f"    {ln}")
        if len(lines) > a.show_errors:
            print(f"    ... {len(lines) - a.show_errors} more lines")
    if failures:
        return 1

    exe, r = link(objects, fc)
    if r.returncode == 0:
        print(f"\nlinked {exe}")
        return 0
    print("\nlink failed:")
    missing = sorted({m.split("`")[1].split("'")[0]
                      for m in r.stderr.splitlines()
                      if "undefined reference" in m})
    if missing:
        print(f"  {len(missing)} undefined symbols:")
        for m in missing:
            print(f"    {m}")
    else:
        for ln in r.stderr.splitlines()[:a.show_errors]:
            print(f"  {ln}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
