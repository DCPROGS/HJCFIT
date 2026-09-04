########################
#   HJCFIT computes missed-events likelihood as described in
#   Hawkes, Jalali and Colquhoun (1990, 1992)
#
#   Copyright (C) 2013  University College London
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#########################

""" The simplex the Fortran HJCFIT used, restored.

    Every published HJCFIT result, including Colquhoun, Hatton & Hawkes (2003),
    was obtained with the simplex in ``SIMPHJC.FOR`` (``DCPROGS/DCFORTRAN``,
    ``Fort90/HJCFIT/``). This is a port of it, with that program's own defaults.

    It is **not** Nelder--Mead as :func:`scipy.optimize.minimize` implements it,
    and the differences are not cosmetic:

    * the starting simplex is the Spendley regular one, built from a *fixed*
      step, so in log space every parameter starts the same factor away from
      the guess however large the rate is. SciPy perturbs each coordinate by
      5 per cent of its own value;
    * the acceptance rule counts how many vertices the reflection beats and
      branches on ``L >= 2``, ``L == 1``, ``L == 0``; at ``L == 1`` it accepts
      the reflection *and then contracts as well*. SciPy compares the
      reflection with the second-worst vertex;
    * shrinkage uses the contraction factor rather than always halving;
    * convergence is on the extent of the simplex alone, per parameter, with no
      function-value tolerance;
    * ending is a choice between four candidates -- the best vertex, the
      average of the vertices, the best point seen anywhere during the run, and
      the result of a local search -- and if one of the last two wins the search
      restarts there with the step reset to ``resfac * crtstp``.

    On the fits of Colquhoun, Hatton & Hawkes (2003) it reaches the same maximum
    as SciPy's Nelder--Mead on all of 250 records, in a little over half the
    likelihood evaluations.

    A bug in the local search, reproduced on purpose
    ------------------------------------------------

    ``SIMPHJC.FOR`` lines 657--691 step each parameter up and then down by
    ``crtstp`` and keep any improvement. The "up" branch sets
    ``pnew1(j) = temp(j) + crtstp(j)`` and then evaluates ``FUNC(kt, temp, ...)``
    -- the *unaltered* best vertex, not the perturbed point -- so it re-measures
    the point it is stepping away from and can never beat it. ``pnew1`` is also
    never reinitialised, so the "down" branch evaluates a vector whose untested
    coordinates still hold the last extension or contraction the iteration
    happened to try. The local search is close to inert in practice.

    ``local_search="fortran"`` reproduces this and is the default, because the
    published numbers came from the code as it stands. ``"corrected"`` does the
    coordinate search as evidently intended and ``"off"`` skips it; having all
    three is what lets the difference be measured rather than assumed.
"""
__docformat__ = "restructuredtext en"
__all__ = ['simplex', 'SimplexResult', 'SIMPLEX_DEFAULTS']

from collections import namedtuple

import numpy as np

#: HJCFIT's own defaults, read from ``Hjcfit1-09122003.for`` lines 2749--2830.
#: ``stpfac`` is the value the program prompts for; in log mode it passes the
#: logarithm of it to the simplex (line 2794, ``stpsav = dlog(stpfac)``).
SIMPLEX_DEFAULTS = dict(
    stpfac_log=5.0,      # line 2772, "if(logsav) stpfac=5.d0"
    stpfac_lin=0.2,      # line 2771
    confac=0.5,          # line 2770
    resfac=10.0,         # line 2773
    nresmax=3,           # line 2808, "irestrt=3"
    errfac=1.0e-3,       # line 2750
    maxfev=20000,        # line 2855
)

_REFFAC = 1.0            # reflection coefficient, SIMPHJC.FOR line 251
_EXTFAC = 2.0            # extension factor, line 252

#: What ``iconv`` means, as SIMPHJC's own header documents it.
ICONV = {2: "best vertex", 3: "average of the vertices", 4: "best point seen",
         5: "local search", 6: "no convergence"}

_FIELDS = ['x', 'fun', 'nfev', 'nit', 'nrestarts', 'iconv', 'success',
           'message']


class SimplexResult(namedtuple('SimplexResult', _FIELDS)):
    """ The result of a :func:`simplex` run.

        The field names are SciPy's where SciPy has one, so that a script can
        swap :func:`scipy.optimize.minimize` for :func:`simplex` and keep
        reading ``res.x``, ``res.fun``, ``res.nfev``, ``res.nit``,
        ``res.success`` and ``res.message``.

        :param x: the point returned, in the space searched
        :param fun: the function value there
        :param nfev: likelihood evaluations used
        :param nit: iterations
        :param nrestarts: restarts taken, at most ``nresmax``
        :param iconv: which of the four endings was taken; see :data:`ICONV`.
            6 means the evaluation budget ran out
        :param success: False only for ``iconv == 6``
        :param message: ``iconv`` in words
    """
    __slots__ = ()

    def __repr__(self):
        return ("SimplexResult(fun={0.fun!r}, nfev={0.nfev}, nit={0.nit}, "
                "nrestarts={0.nrestarts}, success={0.success}, "
                "message={0.message!r})".format(self))


def simplex(fun, x0, args=(), logfit=True, stpfac=None, confac=None,
            resfac=None, nresmax=None, errfac=None, maxfev=None,
            local_search="fortran", callback=None):
    """ Minimise ``fun`` the way the Fortran HJCFIT did.

        :param fun:
           Callable, ``fun(x, *args) -> float``. In HJCFIT this is minus the
           log-likelihood.
        :param x0:
           Starting point. **In log mode this is ``log(rates)``**: the caller
           takes logarithms before the call and exponentiates after, which is
           what the Fortran does around its own call to SIMPHJC.
        :param args:
           Extra arguments passed to ``fun``.
        :param logfit:
           Use the log-space step and convergence rules. This does **not**
           transform ``x0``; it says what ``x0`` already is, so that the step
           can be the same for every parameter and the convergence test can be
           absolute in ``log(rate)``. HJCFIT's own prompt, "Use log(rate
           constant) in Simplex", defaults to yes in every version of the source
           in ``DCPROGS/DCFORTRAN``, so True is the faithful default.
        :param stpfac:
           The step *factor* as the program prompts for it: 5.0 in log mode,
           0.2 in linear. In log mode its logarithm becomes the step, so every
           parameter starts a factor of ``stpfac`` from the guess.
        :param confac:
           Contraction factor, also used for shrinkage.
        :param resfac:
           On a restart the step is reset to ``resfac * crtstp``.
        :param nresmax:
           Most restarts allowed. 0 disables them.
        :param errfac:
           Converged when the simplex spans less than this in every parameter;
           in log mode it is absolute in ``log(rate)``, so 1e-3 is a tenth of a
           per cent of the rate.
        :param maxfev:
           Evaluation budget. Exhausting it returns ``success=False``.
        :param local_search:
           ``"fortran"``, ``"corrected"`` or ``"off"``; see the module
           docstring, which explains why the default reproduces a bug.
        :param callback:
           Called as ``callback(x)`` with the best vertex after each iteration,
           as :func:`scipy.optimize.minimize` calls it.

        :returns: a :class:`SimplexResult`

        .. code-block:: python

            from HJCFIT.likelihood.optimization import simplex

            result = simplex(lambda x: -likelihood(np.exp(x)), np.log(theta))
            rates = np.exp(result.x)
    """
    if local_search not in ("fortran", "corrected", "off"):
        raise ValueError("local_search must be 'fortran', 'corrected' or 'off'")

    d = SIMPLEX_DEFAULTS
    if stpfac is None:
        stpfac = d["stpfac_log"] if logfit else d["stpfac_lin"]
    confac = d["confac"] if confac is None else confac
    resfac = d["resfac"] if resfac is None else resfac
    nresmax = d["nresmax"] if nresmax is None else nresmax
    errfac = d["errfac"] if errfac is None else errfac
    maxfev = d["maxfev"] if maxfev is None else maxfev

    theta = np.array(x0, dtype=float).ravel()
    k = theta.size
    if k < 1:
        raise ValueError("nothing to minimise: x0 is empty")
    n = k + 1

    # step and tolerance, SIMPHJC.FOR lines 213-234. In log mode both are the
    # same for every parameter; in linear mode both scale with the guess.
    if logfit:
        step = np.full(k, float(np.log(stpfac)))
        crtstp = np.full(k, float(errfac))
    else:
        step = stpfac * theta
        crtstp = errfac * theta

    count = [0]

    def evaluate(x):
        count[0] += 1
        return float(fun(x, *args))

    nit = 0
    nrestarts = 0
    # pnew1 persists across iterations exactly as the Fortran array does; the
    # local search depends on what the last extension or contraction left in it
    pnew1 = theta.copy()

    while True:                                   # label 2001, restart target
        simp = np.empty((n, k))
        fval = np.empty(n)
        simp[0] = theta
        fval[0] = evaluate(theta)
        absmin, thmin = fval[0], theta.copy()

        # the Spendley regular simplex, lines 283-291
        fac = (np.sqrt(n) - 1.0) / (k * np.sqrt(2.0))
        for i in range(1, n):
            simp[i] = simp[0] + step * fac
            simp[i][i - 1] = simp[0][i - 1] + step[i - 1] * (fac + 1.0 / np.sqrt(2.0))
            fval[i] = evaluate(simp[i])
            if fval[i] < absmin:
                absmin, thmin = fval[i], simp[i].copy()

        spent = False
        while True:                               # label 2000, iteration
            nit += 1
            ilo, ihi = int(np.argmin(fval)), int(np.argmax(fval))
            flo = fval[ilo]
            if callback is not None:
                callback(simp[ilo].copy())
            if count[0] > maxfev:                 # line 343
                spent = True
                break

            # centroid of every vertex but the worst, divided by k (line 365)
            centre = (simp.sum(axis=0) - simp[ihi]) / k
            pnew = centre - _REFFAC * (simp[ihi] - centre)
            fnew = evaluate(pnew)
            if fnew < absmin:
                absmin, thmin = fnew, pnew.copy()

            if fnew < flo:
                # better than the best vertex, so try extending (line 389)
                pnew1 = centre + _EXTFAC * (pnew - centre)
                fnew1 = evaluate(pnew1)
                if fnew1 < absmin:
                    absmin, thmin = fnew1, pnew1.copy()
                if fnew1 < fnew:
                    simp[ihi], fval[ihi] = pnew1, fnew1
                else:
                    simp[ihi], fval[ihi] = pnew, fnew
            else:
                # L = how many vertices the reflection beats (line 428)
                L = int((fval > fnew).sum())
                if L >= 2:
                    simp[ihi], fval[ihi] = pnew, fnew
                else:
                    if L == 1:
                        # accept the reflection, then contract as well
                        simp[ihi], fval[ihi] = pnew, fnew
                    # contract on the original side of the centroid (line 452)
                    pnew1 = centre + confac * (simp[ihi] - centre)
                    fnew1 = evaluate(pnew1)
                    if fnew1 < absmin:
                        absmin, thmin = fnew1, pnew1.copy()
                    if fnew1 <= fval[ihi]:
                        simp[ihi], fval[ihi] = pnew1, fnew1
                    else:
                        # shrink towards the best vertex, by confac rather than
                        # a hard half (line 478)
                        for i in range(n):
                            simp[i] = simp[ilo] + confac * (simp[i] - simp[ilo])
                            fval[i] = evaluate(simp[i])
                            if fval[i] < absmin:
                                absmin, thmin = fval[i], simp[i].copy()

            # convergence: the extent of the simplex in every coordinate
            # (lines 537-582)
            if np.all(simp.max(axis=0) - simp.min(axis=0) <= np.abs(crtstp)):
                break

        if spent:
            i = int(np.argmin(fval))
            return SimplexResult(
                x=simp[i].copy(), fun=float(fval[i]), nfev=count[0], nit=nit,
                nrestarts=nrestarts, iconv=6, success=False,
                message="no convergence after {0} evaluations".format(maxfev))

        # ---- converged: choose between four candidates (line 604 onwards) ---
        ihi = int(np.argmin(fval))                # Hill's correction, line 614
        temp, fnewlo = simp[ihi].copy(), fval[ihi]

        cand_f = [fnewlo, None, absmin, fnewlo]
        cand_x = [temp, None, thmin.copy(), temp]

        average = simp.mean(axis=0)               # line 633
        cand_f[1], cand_x[1] = evaluate(average), average

        if local_search != "off":                 # line 657
            for j in range(k):
                if local_search == "fortran":
                    # the bug: pnew1[j] is set but *temp* is evaluated
                    pnew1[j] = temp[j] + crtstp[j]
                    f_up = evaluate(temp)
                else:
                    trial = temp.copy()
                    trial[j] = temp[j] + crtstp[j]
                    f_up, pnew1 = evaluate(trial), trial
                if f_up < fnewlo:
                    cand_f[3], cand_x[3] = f_up, pnew1.copy()
                    break

                if local_search == "fortran":
                    pnew1[j] = temp[j] - crtstp[j]
                    f_dn = evaluate(pnew1)
                else:
                    trial = temp.copy()
                    trial[j] = temp[j] - crtstp[j]
                    f_dn, pnew1 = evaluate(trial), trial
                if f_dn < fnewlo:
                    cand_f[3], cand_x[3] = f_dn, pnew1.copy()
                    break

        il = int(np.argmin(cand_f))               # line 703
        if nrestarts < nresmax and il in (2, 3):
            nrestarts += 1
            theta = np.asarray(cand_x[il], dtype=float).copy()
            step = resfac * crtstp                # line 741
            continue                              # goto 2001

        return SimplexResult(
            x=np.asarray(cand_x[il], dtype=float).copy(), fun=float(cand_f[il]),
            nfev=count[0], nit=nit, nrestarts=nrestarts, iconv=il + 2,
            success=True,
            message="return with " + ICONV[il + 2])
