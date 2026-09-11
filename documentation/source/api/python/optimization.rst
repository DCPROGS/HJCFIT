.. _python_optimization_api:

Optimization
------------

.. currentmodule:: HJCFIT.likelihood.optimization

HJCFIT computes a likelihood; maximising it is the caller's business. Two
things here help with that: the simplex the Fortran HJCFIT used, and the
mapping from a reaction graph to the free parameters of a fit.

The simplex HJCFIT used
"""""""""""""""""""""""

Every published HJCFIT result, Colquhoun, Hatton & Hawkes (2003) among them,
was obtained with the simplex in ``SIMPHJC.FOR``. It is **not**
:func:`scipy.optimize.minimize` with ``method='Nelder-Mead'``, and the
differences are not cosmetic — see :func:`simplex_hjc` for what they are.

.. autofunction:: simplex_hjc

.. autoclass:: SimplexHJCResult

.. autodata:: SIMPLEX_HJC_DEFAULTS

Keeping rate constants in range
"""""""""""""""""""""""""""""""

``SIMPHJC.FOR`` has no notion of bounds, and neither does the port: the limits
belonged to the *program* around the subroutine. If you search the rate
constants themselves you have to supply that yourself, or the search will
happily propose a negative rate constant.

.. autofunction:: reset_out_of_range

A worked shape, for one record and a mechanism whose free rates you can get at
as a vector:

.. code-block:: python

    from numpy import log, exp
    from HJCFIT.likelihood import Log10Likelihood
    from HJCFIT.likelihood.optimization import simplex_hjc, reset_out_of_range

    likelihood = Log10Likelihood(bursts, nopen=nopen, tau=tau, tcritical=tcrit)

    def cost(x):
        """Negative log10 likelihood at log rate constants ``x``."""
        mechanism.set_free_rates(exp(x))
        return -likelihood(mechanism.Q)

    result = simplex_hjc(reset_out_of_range(cost, lower=1e-12, upper=1e6),
                         log(guess))
    print(result.fun, result.nfev, result.message)
    rates = exp(result.x)

Three things about that are worth stating rather than leaving to be
rediscovered.

The search is over ``log(rate)`` because that is HJCFIT's own default, it is
three to four times faster (Colquhoun, Hatton & Hawkes 2003, p. 702), and a
logarithm cannot go negative. The fits of that paper's Figures 2–5 and 12–13
were made over the rates themselves, which is why the resetting exists at all.

``result.fun`` is the value of whatever you minimised. If that is the negative
log\ :sub:`10` likelihood, as above, then the log\ :sub:`10` likelihood is
``-result.fun``. Mixing the two conventions between the objective and the
running best is an easy mistake and does not announce itself.

Anything that treats the log likelihood as a statistical quantity — a Hessian,
and so the covariance matrix, the standard deviations and the likelihood
intervals — needs **natural** logarithms, not log\ :sub:`10`. Getting that
wrong inflates every standard deviation by exactly :math:`\sqrt{\ln 10} =
1.517`, which looks like a badly behaved fit rather than a units error.

What it is not for
""""""""""""""""""

``simplex_hjc`` is a local refiner started from a guess a person believes in,
searching the logarithms of the rate constants. That is how HJCFIT used it, and
it is where it works.

It is not a global search from a random Q matrix. Tried that way on CH82 with
eight free parameters it does not converge at all, and the reason is not the
simplex: of 200 random reduced coordinate vectors, **6 give a finite
likelihood**. The starting simplex has nine vertices, so almost every one of
them lands somewhere the likelihood cannot be evaluated and there is nothing to
walk on. A derivative-free simplex needs a valid starting simplex, not merely a
valid starting point.

For a random multi-start over a mechanism you have no guess for, a method that
handles constraints directly — ``COBYLA`` or ``SLSQP`` through
:func:`scipy.optimize.minimize` — will get you into a plausible region.
``simplex_hjc`` is what to finish with, and what to quote.

Reducing a likelihood to its free parameters
""""""""""""""""""""""""""""""""""""""""""""

.. autofunction:: reduce_likelihood
