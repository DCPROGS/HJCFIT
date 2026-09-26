.. _python_mcmc_api:

Bayesian inference
------------------

.. currentmodule:: HJCFIT.likelihood.mcmc

:ref:`python_fitting_api` finds the rate constants that maximise the HJC
likelihood. This module treats them as random variables instead, and describes
their posterior distribution: the likelihood of the records multiplied by a
prior. The approach follows Epstein, Calderhead, Girolami & Sivilotti (2016)
*Biophys J* 111:333-348, who sampled that posterior by Markov chain Monte
Carlo using the same exact missed-events likelihood, and it uses their two
samplers.

Like the fitter, it needs a mechanism from SCALCS (``pip install
hjcfit[fitting]``), and it is duck-typed on one, so it does not import scalcs
itself.

The posterior
"""""""""""""

.. code-block:: python

    import HJCFIT
    from HJCFIT.likelihood.fitting import Record
    from HJCFIT.likelihood.mcmc import LogPosterior
    from scalcs.samples import samples

    bursts = HJCFIT.read_idealized_bursts("CH82", tau=1e-4, tcrit=4e-3)
    record = Record(conc=100e-9, tres=1e-4, tcrit=4e-3,
                    groups=tuple(tuple(b) for b in bursts))

    mechanism = samples.CH82()
    posterior = LogPosterior(mechanism, [record])

    posterior(mechanism.theta())      # ln prior + ln likelihood

.. autoclass:: LogPosterior
   :members: log_prior, log_likelihood, __call__, k, nevals, nfailures

Two things differ from fitting on purpose:

* **Nothing is clipped.** :py:class:`~HJCFIT.likelihood.fitting.HJCFitter`
  resets a rate that leaves its limits, which suits a search. For a
  posterior, a point outside the prior has zero density. Moving it to the
  boundary instead would pile probability there.
* **A likelihood that cannot be computed gives** ``-inf``, **not a penalty.**
  A sampler rejects the move. :py:attr:`LogPosterior.nfailures` counts these,
  so a mechanism that the likelihood struggles with is visible.

Priors
""""""

.. autoclass:: UniformPrior
   :members: from_mechanism, logpdf, contains, sample, k

.. autoclass:: LogUniformPrior

The default prior is :py:meth:`UniformPrior.from_mechanism`: flat between the
limits each free rate already carries. Those are the limits a
maximum-likelihood fit resets against, so the mode of the posterior is the
maximum-likelihood estimate. They are also the prior of Epstein et al.
(2016): scalcs' default limits are U(1e-2, 1e6) s\ :sup:`-1`, and
U(1e-2, 1e10) M\ :sup:`-1` s\ :sup:`-1` for association rates.

Solver settings matter
""""""""""""""""""""""

Pass ``solver=`` to reproduce a published value. On that paper's three AChR
records, its root-finding settings (``nmax=2``, tolerances 1e-12) give a
natural log-likelihood 0.38 higher than the defaults at the same rates. See
:py:data:`~HJCFIT.likelihood.fitting.SOLVER_OPTIONS`.

Sampling
""""""""

Two stages, as in that paper. A Metropolis-within-Gibbs pilot walks the
logarithm of one rate at a time and finds the mode. An adaptive Metropolis
chain, started at that mode, learns the posterior covariance and samples
along correlations that the pilot crosses only slowly:

.. code-block:: python

    from HJCFIT.likelihood.mcmc import (
        adaptive_sample, effective_sample_size, mwg_sample)

    pilot = mwg_sample(posterior, mechanism.theta(), n=10000, burnin=5000,
                       rng=1)
    start, _ = pilot.mode()
    chain = adaptive_sample(posterior, start, n=100000, burnin=50000, rng=2)

    draws = chain.kept()                 # (50000, k) rates, burn-in dropped
    effective_sample_size(draws)         # one value per rate

Step sizes are tuned only during burn-in, so the samples that are kept come
from a chain with a fixed kernel. The adaptive covariance keeps learning
afterwards, but with a weight that falls as :math:`1/n`, which is what keeps
the chain ergodic (Haario et al. 2001).

A good first check of any chain is to plot ``chain.log_posterior`` against
the iteration number, and ``chain.proposals`` against ``chain.samples`` for a
correlated pair. Proposals falling well outside the cloud of accepted samples
mean the proposal does not fit the posterior (Epstein et al. 2016, Fig. 3).

.. autofunction:: mwg_sample

.. autofunction:: adaptive_sample

.. autoclass:: Chain
   :members: n, kept, acceptance_rate, mode, save, load

How much a chain is worth
"""""""""""""""""""""""""

Samples from a Markov chain are correlated, so :math:`N` of them carry less
information than :math:`N` independent draws. The effective sample size is how
many independent draws they are worth.

.. autofunction:: effective_sample_size

.. autofunction:: significant_lags

.. autofunction:: autocorrelation
