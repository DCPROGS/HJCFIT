"""Call the vendored 2003 Fortran HJC likelihood from Python.

Writes ``likdrv``'s input file, runs it, and parses the number back. The point
is to put the original Fortran ``HJCLIK`` and HJCFIT's C++
``Log10Likelihood`` in front of the same Q matrix and the same record.

The Q matrix goes across as a matrix, not as rate constants, so nothing about
either side's parameterisation or constraint handling enters the comparison.
What is compared is the likelihood *given a generator*.

Units at the boundary, since they differ between the two sides and a mistake
here would look exactly like a finding: this module takes **seconds** and
**M** throughout, matching HJCFIT's Python API, and ``likdrv`` converts to the
milliseconds HJCFIT's Fortran holds internally.
"""

from __future__ import annotations

import platform
import subprocess
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_EXE = "likdrv.exe" if platform.system() == "Windows" else "likdrv"

#: How to get it, quoted whenever it is missing so that a skip says what to do.
BUILD_HINT = f"python {Path(__file__).parent.name}/build.py"


def exe_for(engine="free"):
    """The program built by one engine.

    ``free`` is the default build: DCPROGS' own code, with public-domain
    EISPACK and a plain LU in place of the eight sources that carry NAG's or
    Numerical Recipes' copyright. ``original`` is the program exactly as it
    was, and needs a DCFORTRAN checkout to build. See ``build.py``.
    """
    return HERE / "build" / engine / _EXE


#: The default engine's program. Kept as a module attribute because most
#: callers want only this one.
EXE = exe_for("free")


def available(engine="free"):
    """Has this engine been built?"""
    return exe_for(engine).exists()


def write_input(path, Q, kA, kB, kC, intervals, amplitudes, conc, tres,
                tcrit=None, chs=False):
    """Write one ``likdrv`` input file.

    Parameters
    ----------
    Q : (k, k) array
        Generator with the concentration already in it, in s^-1, states
        ordered open first. Diagonals are ignored -- ``qset_stub`` recomputes
        them, as ``SETDIAG`` does in the original.
    kA, kB, kC : int
        Open, short-shut and long-shut state counts. ``kD`` is taken as 0.
    intervals, amplitudes : array_like
        The apparent record, in seconds; amplitude 0 means shut.
    conc : float
        Agonist concentration [M]. Only carried for the record's own header.
    tres : float
        Dead time [s].
    tcrit : float or None
        Critical shut time [s]. None fits the record as one group.
    chs : bool
        Use CHS vectors (Colquhoun, Hawkes & Srodzinski 1996, eqn 4) rather
        than equilibrium vectors (eqn 3).
    """
    Q = np.asarray(Q, float)
    k = Q.shape[0]
    if Q.shape != (k, k):
        raise ValueError("Q must be square")
    if kA + kB + kC != k:
        raise ValueError(f"kA+kB+kC = {kA + kB + kC} but Q is {k}x{k}")
    t = np.asarray(intervals, float)
    a = np.asarray(amplitudes, float)
    if t.size != a.size:
        raise ValueError("intervals and amplitudes differ in length")

    lines = [f"{k} {kA} {kB} {kC} 0",
             repr(float(conc)),
             repr(float(tres)),
             repr(float(tcrit) if tcrit else 0.0),
             "1" if chs else "0"]
    lines += [repr(float(v)) for v in Q.ravel()]
    lines.append(str(t.size))
    lines += [f"{ti!r} {ai!r}" for ti, ai in zip(t.tolist(), a.tolist())]
    Path(path).write_text("\n".join(lines) + "\n")


#: HJCLIK works in natural logarithms and returns the value a minimiser
#: wants. ``hjclik.for`` line 1306 accumulates ``dlog(OLIK(i))`` and line 1318
#: adds ``nscal * 230.2585093``, which is ln(10^100), the underflow
#: rescaling; line 1386 then negates the total. HJCFIT's C++
#: ``Log10Likelihood`` returns ``+log10 L``. So the two differ by a sign and
#: by ln(10), and comparing them without both corrections makes the Fortran
#: look 2.3 times more curved than it is.
LN10 = 2.302585092994046


def log10_likelihood(Q, kA, kB, kC, intervals, amplitudes, conc, tres,
                     tcrit=None, chs=False, keep=None, timeout=600,
                     raw=False, engine="free"):
    """Run the Fortran likelihood once and return it as **+log10 L**.

    Parameters
    ----------
    raw : bool
        Return ``HJCLIK``'s own value instead -- minus the natural log
        likelihood -- for checking the conversion.
    engine : str
        ``"free"`` or ``"original"``; see :func:`exe_for`.

    Returns
    -------
    value : float
        ``+log10 L``, directly comparable with HJCFIT's C++
        ``Log10Likelihood``, unless ``raw``.
    stdout : str
        Everything the program printed, so that a diagnostic or a warning is
        not silently discarded.
    """
    exe = exe_for(engine)
    if not exe.exists():
        raise FileNotFoundError(
            f"{exe} -- run `{BUILD_HINT} --engine {engine}` first")

    if keep is not None:
        path = Path(keep)
        path.parent.mkdir(parents=True, exist_ok=True)
    else:
        fd = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
        fd.close()
        path = Path(fd.name)

    write_input(path, Q, kA, kB, kC, intervals, amplitudes, conc, tres,
                tcrit, chs)
    r = subprocess.run([str(exe), str(path)], capture_output=True, text=True,
                       cwd=str(exe.parent), timeout=timeout)
    if keep is None:
        path.unlink(missing_ok=True)

    value = None
    for line in r.stdout.splitlines():
        if line.startswith("HJCLIK"):
            value = float(line.split()[1])
    if value is None:
        raise RuntimeError(
            f"likdrv produced no HJCLIK line (exit {r.returncode}).\n"
            f"stdout:\n{r.stdout[-2000:]}\nstderr:\n{r.stderr[-2000:]}")
    if raw:
        return value, r.stdout
    return -value / LN10, r.stdout
