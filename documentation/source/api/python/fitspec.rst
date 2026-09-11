.. _python_fitspec_api:

A fit from a file
-----------------

:ref:`python_fitting_api` is the fitting layer as Python objects. This is the
same fit written down: a **fit specification**, which says which records, which
mechanism and how to search, and runs nothing.

That separation is what lets the same description drive a notebook, a command,
a batch runner and -- if one is ever built -- a desktop interface, because a
user interface is a way of editing a fit specification. It is also what makes a
fit reproducible by somebody else: a file they can read, diff and keep beside
the result.

::

    hjcfit template -o my-fit.toml     a specification to edit
    hjcfit check my-fit.toml           say what it would do; run nothing
    hjcfit fit my-fit.toml -o out.json run it, and keep the result

``examples/fit_template.ipynb`` is the same thing as a notebook, over
``examples/CH82.toml``.

Everything here needs the ``[fitting]`` extra::

    pip install 'hjcfit[fitting]'

The file
""""""""

.. code-block:: toml

    title = "CH82 sample record at 100 nM"

    [[data]]
    record = "CH82"        # a sample record, or a path to an .scn file
    conc = 1e-07           # M
    tres = 0.0001          # dead time to impose, s
    tcrit = 0.004          # critical shut time dividing the record, s
    vectors = "chs"        # or "equilibrium"

    [mechanism]
    sample = "CH82"        # a factory in scalcs.samples.samples
    nfree = 8              # free parameters expected; asserted, not applied

    [mechanism.rates]      # the initial guess, by rate name
    beta1 = 15.0
    beta2 = 15000.0

    [search]
    method = "simplex"
    log_params = true
    maxfev = 20000

One ``[[data]]`` section per record. Several means several records contributing
to a single likelihood, so the rate constants are shared and the
concentrations are not.

``hjcfit template`` writes this file with every option in it, commented.

Why TOML and not YAML
"""""""""""""""""""""

A specification is mostly rate constants, and rate constants get written
``1e8``. PyYAML implements YAML 1.1, whose resolver requires an exponent to
carry a sign, so ``yaml.safe_load`` reads ``1e8`` as the **string** ``'1e8'``
-- and so does ``1.0e8``, and so does ``1e+8``. Only ``1.0e+8`` becomes a
float. A rate constant arriving silently as a string is the worst failure this
file could have.

TOML has one number syntax, accepts every form above, and is in the standard
library from Python 3.11; ``tomli`` covers 3.10 and is in the ``[fitting]``
extra. The ``.yaml`` file in ``scalcs/samples`` is a third thing again -- a
pickled Python object graph needing ``unsafe_load`` -- so YAML in this stack
already means something other than a document a person edits.

Check it before you fit it
""""""""""""""""""""""""""

``hjcfit check`` does everything a fit does except the search: it finds the
records, builds the mechanism with every constraint applied, and evaluates the
likelihood once at the initial guess.

.. code-block:: text

    $ hjcfit check examples/CH82.toml
    CH82 sample record at 100 nM
      CH82 at 100 nM: 100 us dead time, groups at tcrit 4 ms, chs vectors
      mechanism: CH82
        guess: beta1 = 15, beta2 = 15000
      search: simplex, over logarithms, up to 20000 evaluations

    100 nM: 4312 -> 1100 intervals at 100 us -> 572 groups, 836 openings, CHS vectors

    CH82: 5 states, 2 open, 10 rate constants
    8 free: beta1, beta2, alpha1, alpha2, k(-1), 2k(-2), 2k(+1), k(+2)
    2 fixed or constrained: k*(+2), 2k*(-2)

    log10L at the guess: 2286.9746
    nothing was fitted; run: hjcfit fit examples/CH82.toml

A misspelled rate name, a critical time below the dead time, or a
free-parameter count nobody expected all surface here in a second, rather than
twenty minutes into a search or in a set of estimates that look plausible.

Two fields rather than one signed number
""""""""""""""""""""""""""""""""""""""""

Elsewhere in the stack a negative critical time is a flag selecting
equilibrium vectors while its magnitude is still the time that divides the
record -- see :py:func:`HJCFIT.read_idealized_bursts`. Here they are
``tcrit`` and ``vectors``, because the critical time does two jobs and only
one of them is usually being changed:

* **``tcrit`` with ``vectors = "chs"``** -- groups cut at a critical time
  chosen to separate the activations of one channel, started and ended with
  the CHS vectors of Colquhoun, Hawkes & Srodzinski (1996). The usual case.
* **``tcrit`` with ``vectors = "equilibrium"``** -- still groups, but
  equilibrium vectors: clusters at a concentration high enough to desensitise.
* **no ``tcrit``** -- the whole record as a single group, which assumes one
  channel throughout. ``vectors`` must then be ``"equilibrium"``, because CHS
  vectors are defined between groups.

Rates are named, never numbered
"""""""""""""""""""""""""""""""

Every reference to a rate constant in a specification is by name. SCALCS
addresses rates by index into ``mec.Rates`` -- ``set_mr(True, 5, 0)`` -- and an
index is what users get wrong silently, because rate 5 of a mechanism is
whatever the sample happened to list fifth.
:py:func:`~HJCFIT.likelihood.runner.build_mechanism` is the one place a name
becomes an index, and it lists the names that exist when one does not.

It also refuses two things SCALCS allows:

* **``mr`` naming a rate outside the chosen cycle.** SCALCS warns on stderr
  and carries on, which is worse than nothing: the rate leaves the
  free-parameter list while the cycle's constraint stays where it was, so the
  rate is neither fitted nor computed -- it is frozen at its initial guess and
  the fit has one parameter fewer than the file says.
* **A second rate of a cycle left flagged.** A cycle determines exactly one
  rate, and ``set_mr`` does not clear the previous one; its own ``update_mr``
  carries the comment ``TODO: check for consistency between cycle.mrconstr and
  rate.mr``. On CH82, whose sample already has ``2k*(-2)`` under microscopic
  reversibility, asking for ``mr = "beta1"`` gave seven free parameters rather
  than eight with ``2k*(-2)`` stuck at whatever the sample carried.
  ``build_mechanism`` releases the stale one.

Rates that end on a limit
"""""""""""""""""""""""""

A rate driven outside its limits is reset before the likelihood sees it, which
is what HJCFIT did, and nothing announces it. Such a rate is not an estimate:
it is a statement that the likelihood wanted to go somewhere the model forbids.
:py:func:`~HJCFIT.likelihood.runner.against_limits` finds them and they are
printed apart from the fitted values.

.. code-block:: text

    log10(L) = 2288.7861   1154 evaluations in 1.3 s
      beta1                5
      ...

    Rates that ended on a limit -- these are not estimates:
      beta1                5  at its upper limit of 5

Note that SCALCS gives every rate default limits -- 10\ :sup:`-15` to
10\ :sup:`9` for a concentration-dependent rate and 10\ :sup:`-15` to
10\ :sup:`6` for the rest -- so this can happen in a specification that sets
no limits at all.

The specification
"""""""""""""""""

.. currentmodule:: HJCFIT.likelihood.fitspec

.. autoclass:: FitSpec
   :members: from_toml, from_dict, to_toml, write_toml, as_dict, validate

.. autoclass:: DataSpec
   :members: from_dict, validate

.. autoclass:: MechanismSpec
   :members: from_dict, validate

.. autoclass:: SearchSpec
   :members: from_dict, validate

.. autoexception:: SpecError

.. autofunction:: load_toml_bytes

.. autodata:: TEMPLATE
   :annotation:

Carrying it out
"""""""""""""""

.. currentmodule:: HJCFIT.likelihood.runner

.. autofunction:: run

.. autoclass:: Outcome

.. autofunction:: load_records

.. autofunction:: build_mechanism

.. autofunction:: against_limits

.. autofunction:: provenance

.. autofunction:: result_as_dict

.. autofunction:: write_result

.. autoexception:: RunnerError

The command
"""""""""""

.. currentmodule:: HJCFIT.likelihood.cli

.. autofunction:: main

Also reachable as ``python -m HJCFIT.likelihood.cli`` when the console script
is not on the path, which is the usual state of affairs inside a conda
environment on Windows.

Exit status is 0 on success, 1 on anything the command can explain -- one line
beginning with ``hjcfit:`` -- and 2 from ``fit`` when the search did not
converge. The estimates are still printed in that case; it is the status that
says so, because that is what a script reads.
