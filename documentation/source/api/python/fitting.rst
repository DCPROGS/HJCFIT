.. _python_fitting_api:

Fitting a mechanism
-------------------

.. currentmodule:: HJCFIT.likelihood.fitting

:py:class:`~HJCFIT.likelihood.Log10Likelihood` takes a Q matrix and returns a
number. Between that and a person sits this module: a record in the form the
likelihood wants, a search over the free rate constants of a mechanism, and a
result worth printing.

It needs a mechanism, and HJCFIT does not define one — mechanisms, rate
constants and constraints live in SCALCS. Install it with the extra::

    pip install hjcfit[fitting]

This module is imported on demand rather than exposed from
HJCFIT.likelihood itself::

    from HJCFIT.likelihood.fitting import HJCFitter, Record

That is not a stylistic choice. Exposing it from the package __init__ made
HJCFIT.likelihood import itself while partially initialised, and every CI
job failed on it.

**Nothing here imports scalcs.** The fitter is duck-typed on a mechanism:
anything offering ``theta()``, ``theta_unsqueeze()``, ``Rates``, ``kA``,
``set_eff()`` and ``Q`` will do, which in practice means
:class:`scalcs.mechanism.Mechanism`. That is what keeps ``pip install hjcfit``
wanting nothing but numpy and scipy.

A whole fit
"""""""""""

.. code-block:: python

    import HJCFIT
    from HJCFIT.likelihood.fitting import HJCFitter, Record
    from scalcs.samples import samples

    bursts = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
    record = Record(conc=100e-9, tres=1e-4, tcrit=4e-3,
                    groups=tuple(tuple(b) for b in bursts))

    mechanism = samples.CH82()
    mechanism.set_eff("c", 100e-9)

    result = HJCFitter(mechanism, [record]).fit()
    print(result)

.. code-block:: text

    log10(L) = 2289.1326   1313 evaluations in 1.1 s
      beta1            11.31
      beta2        1.328e+04
      alpha1            3668
      alpha2           424.4
      k(-1)             1605
      2k(-2)            4410
      2k(+1)       1.471e+04
      k(+2)        5.288e+08

Several concentrations at once
""""""""""""""""""""""""""""""

The likelihood takes **one Q matrix per call**, so a fit to several
concentrations builds one instance per record and adds the logarithms. Pass
them all and :py:class:`HJCFitter` does that::

    fitter = HJCFitter(mechanism, [rec_10uM, rec_30uM, rec_100uM])

The mechanism is set to each record's own concentration before that record is
evaluated, so the rate constants are shared and the concentrations are not.

The record
""""""""""

.. autoclass:: Record
   :members: n_intervals, n_openings, as_lists, check

.. autofunction:: trim_to_openings

Every group must have an **odd** number of intervals — it must start and end
with an opening — because the likelihood is a product of matrices alternating
:math:`A \rightarrow F` and :math:`F \rightarrow A`. Burst segmentation in
``dcio`` already gives odd groups; a whole record fitted as one group generally
does not, which is what :py:func:`trim_to_openings` is for.

An even group does not fail loudly. :py:meth:`Record.check` is called on
construction of the fitter for that reason.

The fitter
""""""""""

.. autoclass:: HJCFitter
   :members: log10_likelihood, ln_likelihood, cost, fit

.. autodata:: FAILURE_COST

The result
""""""""""

.. autoclass:: FitResult

Three things about a fit that are easy to get wrong
""""""""""""""""""""""""""""""""""""""""""""""""""""

**Log space is the default, and not only for speed.** Searching the logarithms
of the rate constants is three to four times faster (Colquhoun, Hatton & Hawkes
2003, p. 702) and cannot produce a negative rate. Searching the rates
themselves can: on that paper's AChR mechanism it gave **four fits in 250 with
negative rate constants** until the out-of-range reset was added.

**A rate that ends against its limit is not a fit.** It is a statement that the
likelihood wanted to go somewhere the model forbids. :py:class:`HJCFitter`
resets rather than constrains, as HJCFIT did, so nothing tells you this
happened — compare the result against the limits you set.

**Standard errors need natural logarithms.** Anything treating the log
likelihood as a statistical quantity — a Hessian, and so the covariance matrix,
the standard deviations and the likelihood intervals — needs
:py:meth:`HJCFitter.ln_likelihood`, not the base-10 one. Using
log\ :sub:`10` inflates every standard deviation by exactly
:math:`\sqrt{\ln 10} = 1.517`, which looks like a badly behaved fit rather than
a units error.

Which search
""""""""""""

``fit(search="simplex")``, the default, is
:py:func:`~HJCFIT.likelihood.optimization.simplex_hjc` — HJCFIT's own, the
search that produced every published result. ``fit(search="scipy")`` is SciPy's
Nelder–Mead with a restart loop.

Reach for SciPy when the starting point is poor. A regular simplex needs every
one of its vertices to be evaluable, and a random Q matrix usually cannot offer
that: on CH82 with eight free parameters, only 6 of 200 random reduced
coordinate vectors give a finite likelihood. See
:ref:`python_optimization_api` for that in more detail.

On a real fit from a sensible guess the two agree. Over CH82 at 100 nM both
reach log\ :sub:`10` L = 2289.13 from 2286.97, by different paths and to
slightly different parameters — the :math:`\alpha`–:math:`\beta` ridge, which
that paper measures at :math:`r = 0.92`. Agreement on the maximum while
disagreeing on where it sits is evidence about the likelihood rather than about
either search.
