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

""" Carrying out what a :py:class:`~HJCFIT.likelihood.fitspec.FitSpec` says.

    :py:mod:`~HJCFIT.likelihood.fitspec` describes a fit and executes nothing.
    This module is the half that executes: it reads the records through
    ``dcio``, builds the mechanism through ``scalcs``, runs the search through
    :py:class:`~HJCFIT.likelihood.fitting.HJCFitter`, and writes a result with
    enough beside it to be believed a year later.

    ::

        from HJCFIT.likelihood.fitspec import FitSpec
        from HJCFIT.likelihood.runner import run

        outcome = run(FitSpec.from_toml("ch82.toml"))
        print(outcome.result)

    Needs the ``[fitting]`` extra: ``pip install 'hjcfit[fitting]'``.

    Three things this module does that a hand-written script usually does not,
    each because it has gone wrong:

    * **It names rates rather than numbering them.** SCALCS addresses rates by
      index into ``mec.Rates``. :py:func:`build_mechanism` is the one place a
      name becomes an index, and it lists the names that exist when one does
      not.
    * **It reports rates that end against a limit.** A rate outside its limits
      is reset before the likelihood sees it and nothing says so, which means
      a fit can end with a rate sitting exactly on a bound and look like an
      estimate. Note that SCALCS gives every rate default limits --
      1e-15 to 1e9 for a concentration-dependent rate and 1e-15 to 1e6 for the
      rest -- so this can happen in a specification that sets no limits at all.
    * **It asserts the free-parameter count** when the specification states
      one, before the search starts rather than after it.
"""
__docformat__ = "restructuredtext en"
__all__ = ['Outcome', 'build_mechanism', 'load_records', 'run',
           'against_limits', 'provenance', 'result_as_dict', 'write_result',
           'RunnerError', 'SampleOutcome', 'sample', 'samples_as_dict',
           'write_samples']

import json
import os
import platform
import socket
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np

from .fitting import HJCFitter, Record, trim_to_openings


class RunnerError(RuntimeError):
    """ A specification that cannot be carried out.

        Distinct from :py:exc:`~HJCFIT.likelihood.fitspec.SpecError`, which is
        a specification that cannot be *read*. A name that no rate of the
        chosen mechanism carries is well-formed and still impossible, and it
        is this.
    """


def _need(module, what):
    """ Import an extra, and say what to install when it is not there. """
    try:
        __import__(module)
    except ImportError as error:
        raise RunnerError(
            "{0} needs {1}, which is not a dependency of HJCFIT -- the "
            "likelihood itself wants only numpy and scipy. Install the extra: "
            "pip install 'hjcfit[fitting]'. (Original error: {2})"
            .format(what, module.split('.')[0], error))
    return sys.modules[module]


# ------------------------------------------------------------------ records

def _resolve_record(name):
    """ A path, or the name of a sample record shipped with HJCFIT.

        The same lookup :py:func:`HJCFIT.read_idealized_bursts` does, so a
        specification can say ``record = "CH82"`` and run anywhere the wheel
        is installed.
    """
    if os.path.exists(name):
        return name
    from os.path import abspath, dirname, join, splitext
    # HJCFIT/likelihood/runner.py -> HJCFIT/data
    sample = join(dirname(dirname(abspath(__file__))), 'data',
                  '{0}.scn'.format(splitext(name)[0]))
    if os.path.exists(sample):
        return sample
    raise RunnerError(
        "no file or sample record {0!r}. The samples that ship with HJCFIT "
        "are CH82, CO and CCO; anything else has to be a path to an .scn "
        "file.".format(name))


def load_records(spec):
    """ Read every record a specification names.

        :param spec: A :py:class:`~HJCFIT.likelihood.fitspec.FitSpec`.
        :returns: A list of :py:class:`~HJCFIT.likelihood.fitting.Record`, in
            the order the specification gives them.
        :raises RunnerError: on a record that cannot be found or read.

        Two paths, because a record divided into groups and a record fitted
        whole are genuinely different things:

        * with a critical time, the record's **periods** are segmented by
          :py:func:`dcio.analysis.bursts_from_record`. Periods rather than
          resolved intervals because imposing a dead time emits a fresh open
          interval at every change of fitted amplitude, so a record idealised
          with sub-conductance levels does not alternate open and shut, and
          the likelihood is a product of matrices that must.
          :py:func:`HJCFIT.read_idealized_bursts` is the other caller of that
          pairing and carries the argument in full; a test here requires the
          two to produce identical groups rather than similar ones.
        * without one, the whole record is a single group, trimmed to start
          and end on an opening. That assumes one channel for the length of
          the record, which is why it is the exception.

        The critical time is used twice and means two different things, which
        is why the specification separates them: its magnitude divides the
        record, and the choice of vectors decides what the likelihood is told.
        Groups cut with ``vectors = "equilibrium"`` -- clusters at a
        desensitising concentration -- pass None as *tcritical* while still
        being groups.
    """
    scn = _need('dcio.formats.scn', 'reading a record')
    analysis = _need('dcio.analysis', 'reading a record')

    records = []
    for i, data in enumerate(spec.data):
        path = _resolve_record(data.record)
        try:
            raw = scn.read(path)
            screc = analysis.from_scn(raw, tres=data.tres)
        except Exception as error:
            raise RunnerError("data[{0}]: cannot read {1}: {2}: {3}".format(
                i, path, type(error).__name__, error))

        if data.tcrit is None:
            periods = screc.periods
            groups = (tuple(trim_to_openings(periods.intervals,
                                             periods.amplitudes)),)
        else:
            groups = tuple(
                tuple(np.asarray(g, dtype=float))
                for g in analysis.bursts_from_record(screc, data.tcrit,
                                                     intervals_only=True))
        if not groups or not groups[0]:
            raise RunnerError(
                "data[{0}]: {1} yields no intervals at a dead time of "
                "{2:g} us{3}".format(i, path, data.tres * 1e6,
                                     "" if data.tcrit is None else
                                     " and a critical time of {0:g} ms"
                                     .format(data.tcrit * 1e3)))

        records.append(Record(
            conc=data.conc,
            groups=groups,
            tres=data.tres,
            # The magnitude divided the record above; here it selects the
            # vectors, and None is what asks for equilibrium ones.
            tcrit=data.tcrit if data.vectors == 'chs' else None,
            n_raw=int(len(raw.intervals)),
        ).check())
    return records


# --------------------------------------------------------------- mechanism

def _indices(mec, name, where):
    """ Every index in ``mec.Rates`` whose rate carries *name*. """
    found = [i for i, rate in enumerate(mec.Rates) if rate.name == name]
    if not found:
        raise RunnerError(
            "{0}: this mechanism has no rate called {1!r}. It has: {2}"
            .format(where, name,
                    ", ".join(r.name for r in mec.Rates)))
    return found


def _one_index(mec, name, where):
    """ The single index of *name*, for an interface that takes one. """
    found = _indices(mec, name, where)
    if len(found) > 1:
        raise RunnerError(
            "{0}: {1} names {2} rates of this mechanism, and this setting "
            "applies to one".format(where, name, len(found)))
    return found[0]


def _build_bare(mech):
    """ The mechanism named by the specification, before any constraint. """
    if mech.sample is not None:
        samples = _need('scalcs.samples.samples', 'building a mechanism')
        factory = getattr(samples, mech.sample, None)
        if factory is None or not callable(factory):
            available = sorted(
                name for name, value in vars(samples).items()
                if callable(value) and not name.startswith('_')
                and name not in ('apply_rates',))
            raise RunnerError(
                "mechanism.sample: scalcs.samples.samples has no {0!r}. It "
                "has: {1}".format(mech.sample, ", ".join(available)))
        built = factory()
        # Some factories -- the independent-binding loaders -- return the
        # mechanism and the constraints they applied. The mechanism is first.
        return built[0] if isinstance(built, tuple) else built

    scalcsio = _need('scalcs.scalcsio', 'reading a .mec file')
    if not os.path.exists(mech.mec_file):
        raise RunnerError("mechanism.mec_file: no such file {0!r}"
                          .format(mech.mec_file))
    _version, meclist, _max = scalcsio.mec_get_list(mech.mec_file)
    if not meclist:
        raise RunnerError("mechanism.mec_file: {0} holds no mechanisms"
                          .format(mech.mec_file))
    if mech.mec_number is None:
        if len(meclist) > 1:
            listing = "\n".join(
                "    {0}  {1}  {2}".format(entry[1], entry[2], entry[3])
                for entry in meclist)
            raise RunnerError(
                "mechanism.mec_number: {0} holds {1} mechanisms, so which one "
                "has to be said:\n{2}".format(mech.mec_file, len(meclist),
                                              listing))
        chosen = meclist[0]
    else:
        matching = [e for e in meclist if e[1] == mech.mec_number]
        if not matching:
            raise RunnerError(
                "mechanism.mec_number: {0} has no mechanism {1}; it has {2}"
                .format(mech.mec_file, mech.mec_number,
                        ", ".join(str(e[1]) for e in meclist)))
        chosen = matching[-1]        # the last rate set saved for it
    return scalcsio.mec_load(mech.mec_file, chosen[0])


def build_mechanism(spec):
    """ Build the mechanism a specification names, with its constraints.

        :param spec: A :py:class:`~HJCFIT.likelihood.fitspec.FitSpec`.
        :returns: A ``scalcs.mechanism.Mechanism``, carrying the initial
            guess, ready to be fitted.
        :raises RunnerError: on a name no rate has, or a free-parameter count
            that does not match ``mechanism.nfree``.

        The order matters and is the reproduction's: the guess, then the fixed
        rates, then the limits, then microscopic reversibility, then the
        constraints are updated, and only then the EC50 -- which is computed
        *from* the other rates and so has to come last.
    """
    mech = spec.mechanism
    mec = _build_bare(mech)

    for name, value in mech.rates.items():
        for i in _indices(mec, name, 'mechanism.rates'):
            mec.Rates[i].rateconstants = value

    for name, value in mech.fixed.items():
        for i in _indices(mec, name, 'mechanism.fixed'):
            mec.Rates[i].rateconstants = value
            mec.Rates[i].fixed = True

    for name, (lower, upper) in mech.limits.items():
        for i in _indices(mec, name, 'mechanism.limits'):
            mec.Rates[i].limits = [[lower, upper]]

    if mech.mr is not None:
        if not getattr(mec, 'Cycles', None):
            raise RunnerError(
                "mechanism.mr: this mechanism has no cycles, so no rate of it "
                "is determined by microscopic reversibility")
        if not 0 <= mech.mr_cycle < len(mec.Cycles):
            raise RunnerError(
                "mechanism.mr_cycle: {0} is not a cycle of this mechanism, "
                "which has {1}".format(mech.mr_cycle, len(mec.Cycles)))
        index = _one_index(mec, mech.mr, 'mechanism.mr')
        rate = mec.Rates[index]
        cycle = mec.Cycles[mech.mr_cycle]
        # Checked here rather than left to SCALCS, which warns on stderr and
        # carries on. What it does then is worse than nothing: it marks the
        # rate constrained, so it leaves the free-parameter list, while the
        # cycle's own constraint stays on whichever pair it was already on.
        # The named rate is not computed from anything -- it is silently
        # frozen at its initial guess, and the fit has one parameter fewer than
        # the specification says. Measured on CH82 with mr = "k(-1)", which is
        # not in its cycle: eight free parameters became seven, k(-1) stayed at
        # the guess, and the microscopic-reversibility constraint remained on
        # 2k*(-2).
        states = list(getattr(cycle, 'states', []))
        outside = [s for s in (rate.State1.name, rate.State2.name)
                   if s not in states]
        if outside:
            raise RunnerError(
                "mechanism.mr: {0} connects {1} and {2}, and {3} not in cycle "
                "{4} ({5}). A rate outside the cycle cannot be determined by "
                "microscopic reversibility round it.".format(
                    mech.mr, rate.State1.name, rate.State2.name,
                    "{0} is".format(outside[0]) if len(outside) == 1
                    else "neither is",
                    mech.mr_cycle, ", ".join(states)))

        # A cycle determines exactly one rate, so any other rate of this cycle
        # already flagged is stale and has to be released. SCALCS does not do
        # this: set_mr moves the cycle's own constraint and leaves every
        # previous rate's mr flag set, and its update_mr carries the comment
        # "TODO: check for consistency between cycle.mrconstr and rate.mr".
        # The consequence is the same silent freeze as above, in the supported
        # path rather than a misuse of it. On CH82, whose sample already has
        # 2k*(-2) under microscopic reversibility, asking for mr = "beta1"
        # left seven free parameters rather than eight, with 2k*(-2) neither
        # fitted nor computed from the cycle -- stuck at 0.411, the value the
        # sample happened to carry.
        for i, other in enumerate(mec.Rates):
            if i == index or not getattr(other, 'mr', False):
                continue
            if (other.State1.name in states and other.State2.name in states):
                other.mr = False

        mec.set_mr(True, index, mech.mr_cycle)

    mec.update_constrains()

    if mech.ec50 is not None:
        mec.set_EC50_constraint(
            _one_index(mec, mech.ec50['rate'], 'mechanism.ec50.rate'),
            mech.ec50['value'])

    free = list(mec.get_free_parameter_names())
    if not free:
        raise RunnerError(
            "every rate of this mechanism is fixed or constrained, so there "
            "is nothing to fit")
    if mech.nfree is not None and len(free) != mech.nfree:
        raise RunnerError(
            "mechanism.nfree says {0} free parameters; the mechanism has {1}: "
            "{2}".format(mech.nfree, len(free), ", ".join(free)))
    return mec


# ------------------------------------------------------------------ the fit

def against_limits(mec, tolerance=1e-6):
    """ Free rates that ended on one of their limits.

        :param mec: A mechanism, after a fit.
        :param tolerance: Relative closeness that counts as "on".
        :returns: A list of ``(name, value, 'lower'|'upper', limit)``.

        A rate driven outside its limits is reset to the limit before the
        likelihood sees it, which is what HJCFIT did and what
        :py:class:`~HJCFIT.likelihood.fitting.HJCFitter` does. Nothing
        announces it. A rate sitting on a bound is not an estimate of that
        rate -- it is a statement that the likelihood wanted to go somewhere
        the model forbids -- so it has to be reported separately from the
        fitted values, and this is what does that.
    """
    on = []
    for rate in mec.Rates:
        if not rate.is_free or not rate.limits:
            continue
        value = float(np.ravel(rate.rateconstants)[0])
        lower, upper = rate.limits[0][0], rate.limits[0][1]
        for which, limit in (('lower', lower), ('upper', upper)):
            if limit and abs(value - limit) <= tolerance * abs(limit):
                on.append((rate.name, value, which, float(limit)))
    return on


@dataclass
class Outcome:
    """ Everything one run produced.

        :param spec: The specification it was run from.
        :param result: The :py:class:`~HJCFIT.likelihood.fitting.FitResult`.
        :param records: The records as they were fitted.
        :param mec: The mechanism, carrying the **fitted** rate constants --
            which is what predicted distributions have to be drawn from.
        :param limited:
          What :py:func:`against_limits` found. Empty is the good case.
        :param provenance: What :py:func:`provenance` recorded.
    """

    spec: object
    result: object
    records: list
    mec: object = field(repr=False, default=None)
    limited: list = field(default_factory=list)
    provenance: dict = field(default_factory=dict, repr=False)

    def __str__(self):
        lines = [str(r) for r in self.records]
        lines.append(str(self.result))
        if self.limited:
            lines.append("")
            lines.append("Rates that ended on a limit -- these are not "
                         "estimates:")
            for name, value, which, limit in self.limited:
                lines.append("  {0:9s} {1:12.4g}  at its {2} limit of {3:g}"
                             .format(name, value, which, limit))
        return "\n".join(lines)


def run(spec, verbose=False, store_path=False, x0=None):
    """ Read the records, build the mechanism, fit, and report.

        :param spec: A :py:class:`~HJCFIT.likelihood.fitspec.FitSpec`.
        :param verbose: Print the records and the starting likelihood first.
        :param store_path: Keep the best point of every iteration.
        :param x0: Starting point, in the space being searched. Defaults to
            the mechanism's own guess, which is what the specification says.
        :returns: An :py:class:`Outcome`.
    """
    records = load_records(spec)
    mec = build_mechanism(spec)
    search = spec.search

    fitter = HJCFitter(mec, records, log_params=search.log_params,
                       store_path=store_path,
                       solver=spec.likelihood.solver() or None)
    if verbose:
        # Not the records: Outcome.__str__ prints those, so that a run with
        # verbose off still says what was fitted.
        print("{0} free parameters: {1}".format(
            len(mec.get_free_parameter_names()),
            ", ".join(mec.get_free_parameter_names())))
        print("log10L at the guess: {0:.4f}".format(fitter.log10_likelihood()))

    result = fitter.fit(x0=x0, search=search.method,
                        restarts=search.restarts, maxfev=search.maxfev,
                        maxiter=search.maxiter, xatol=search.xatol,
                        fatol=search.fatol, options=search.options or None)
    return Outcome(spec=spec, result=result, records=records, mec=fitter.mec,
                   limited=against_limits(fitter.mec),
                   provenance=provenance())


# ---------------------------------------------------------------- the record

def provenance():
    """ Where and with what this was computed.

        Package versions, interpreter, host and time. Cheap, and the
        difference between a number in a paper that can be traced and one that
        cannot: the reproduction of Colquhoun, Hatton & Hawkes (2003) attaches
        one of these to every cached result, and it is the reason a
        disagreement there can be chased to a version rather than argued about.

        Deliberately not a git description of the caller's working tree.
        Shelling out to git costs the best part of a second, needs a
        repository to be there, and says nothing at all about a fit run from
        an installed wheel -- which is how the people this is for will run it.
    """
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version as _version

    versions = {}
    for package in ('hjcfit', 'dcio', 'scalcs', 'numpy', 'scipy',
                    'matplotlib'):
        try:
            versions[package] = _version(package)
        except PackageNotFoundError:
            versions[package] = None
    return {
        'when': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'versions': versions,
        'python': sys.version.split()[0],
        'platform': platform.platform(),
        'host': socket.gethostname(),
    }


def result_as_dict(outcome):
    """ An outcome as plain JSON-ready data.

        The specification is written back into it, so the file is a complete
        account: what was asked for, what came out, and what it was computed
        with. Nothing here holds the records -- they are thousands of floats
        and they are in the .scn file the specification names.
    """
    result = outcome.result
    return {
        'hjcfit_result': 1,
        'spec': outcome.spec.as_dict(),
        'fit': {
            'log10_likelihood': result.log10_likelihood,
            'ln_likelihood': result.log10_likelihood * float(np.log(10.0)),
            'free_names': list(result.free_names),
            'free_values': [float(v) for v in result.free_values],
            'rates': {k: float(v) for k, v in result.rates.items()},
            'nevals': int(result.nevals),
            'niter': int(result.niter),
            'nfailures': int(result.nfailures),
            'seconds': float(result.seconds),
            'success': bool(result.success),
            'message': str(result.message),
        },
        'records': [
            {'conc': r.conc, 'tres': r.tres, 'tcrit': r.tcrit,
             'groups': len(r.groups), 'intervals': r.n_intervals,
             'openings': r.n_openings, 'raw_intervals': r.n_raw}
            for r in outcome.records
        ],
        'at_limits': [
            {'rate': name, 'value': value, 'which': which, 'limit': limit}
            for name, value, which, limit in outcome.limited
        ],
        'provenance': outcome.provenance,
    }


def write_result(path, outcome):
    """ Write :py:func:`result_as_dict` to *path* as JSON.

        JSON rather than TOML: this is written and read by programs, and TOML
        has no null, which every absent version and every equilibrium-vector
        critical time needs.
    """
    with open(path, 'w', encoding='utf-8', newline='\n') as handle:
        json.dump(result_as_dict(outcome), handle, indent=2, sort_keys=False)
        handle.write("\n")
    return path


# --------------------------------------------------------------- sampling

def _posterior(spec):
    """ Records, mechanism and log posterior, built from the specification.

        Every chain process builds its own. The likelihood objects wrap C++
        state and cannot be pickled, and a specification can.
    """
    from .mcmc import LogPosterior, LogUniformPrior, UniformPrior

    records = load_records(spec)
    mec = build_mechanism(spec)
    prior_class = (UniformPrior if spec.mcmc.prior == 'uniform'
                   else LogUniformPrior)
    post = LogPosterior(mec, records, prior=prior_class.from_mechanism(mec),
                        solver=spec.likelihood.solver() or None)
    return records, mec, post


def _run_chain(spec, index, start, covariance):
    """ One chain of ``spec.mcmc``: the unit of work of one process. """
    from .mcmc import adaptive_sample, mwg_sample

    _, _, post = _posterior(spec)
    m = spec.mcmc
    start = np.asarray(start, dtype=float)
    kwargs = {}
    if m.log_space is not None:
        kwargs['log_space'] = m.log_space
    log_space = (m.log_space if m.log_space is not None
                 else m.sampler == 'mwg')
    if covariance is not None:
        covariance = np.asarray(covariance, dtype=float)
        if log_space:
            # The same Gaussian, to first order, on the logarithms.
            covariance = covariance / np.outer(start, start)
    if m.sampler == 'adaptive':
        if covariance is not None:
            # Proposals shaped like the posterior at the optimal scale from
            # the first iteration, instead of the paper's small isotropic step
            # while the covariance is learned.
            kwargs.update(initial_covariance=covariance, initial_step=2.38)
        return adaptive_sample(post, start, n=m.n, burnin=m.burnin,
                               rng=m.seed + index, mixture=m.mixture,
                               **kwargs)
    if covariance is not None:
        kwargs['initial_scale'] = 2.38 * np.sqrt(np.diag(covariance))
    return mwg_sample(post, start, n=m.n, burnin=m.burnin,
                      rng=m.seed + index, **kwargs)


def _starts(post, start, covariance, chains, seed):
    """ The first chain at *start*; the others at draws from the Gaussian
        approximation that the prior allows, so that chains which agree have
        had a chance not to. """
    starts = [np.asarray(start, dtype=float)]
    rng = np.random.default_rng(seed + 10007)
    for _ in range(1, chains):
        point = starts[0]
        if covariance is not None:
            for _ in range(100):
                draw = rng.multivariate_normal(starts[0], covariance)
                if np.isfinite(post(draw)):
                    point = draw
                    break
        starts.append(point)
    return starts


def _chain_summary(names, chains):
    """ Pooled posterior summary per free rate, over the kept samples. """
    from .mcmc import effective_sample_size, potential_scale_reduction

    rows = []
    for i, name in enumerate(names):
        series = [c.kept()[:, i] for c in chains]
        pooled = np.concatenate(series)
        ess = 0.0
        for s in series:
            try:
                ess += effective_sample_size(s)
            except ValueError:          # a constant series: nothing moved
                pass
        try:
            rhat = potential_scale_reduction(series)
        except ValueError:
            rhat = float('nan')
        q = np.percentile(pooled, [2.5, 50.0, 97.5])
        rows.append({'name': name, 'mean': float(pooled.mean()),
                     'sd': float(pooled.std(ddof=1)),
                     'q2.5': float(q[0]), 'median': float(q[1]),
                     'q97.5': float(q[2]), 'ess': float(ess),
                     'rhat': float(rhat)})
    return rows


@dataclass
class SampleOutcome:
    """ Everything one sampling run produced.

        :param spec: The specification.
        :param records: The records as they were sampled against.
        :param start: Where the first chain started, as rates.
        :param start_from: ``'fit'`` or ``'pilot'``.
        :param fit: The :py:class:`Outcome` of the fit, for ``start = "fit"``.
        :param pilot: The pilot :py:class:`~HJCFIT.likelihood.mcmc.Chain`, for
          ``start = "guess"``.
        :param approximation: The
          :py:class:`~HJCFIT.likelihood.mcmc.GaussianApproximation` at the
          start, when the Hessian could be used.
        :param chains: One :py:class:`~HJCFIT.likelihood.mcmc.Chain` each.
        :param summary: Per free rate: mean, sd, quantiles, ESS summed over
          chains and split-R-hat across them.
        :param notes: Anything the run had to decide that the reader should
          know, such as falling back from the Hessian.
        :param seconds: Wall-clock time of the whole run.
    """

    spec: object
    records: list
    names: tuple
    start: np.ndarray
    start_from: str
    chains: list
    summary: list
    fit: object = None
    pilot: object = None
    approximation: object = None
    notes: list = field(default_factory=list)
    seconds: float = 0.0
    provenance: dict = field(default_factory=dict, repr=False)

    @property
    def max_rhat(self):
        values = [row['rhat'] for row in self.summary
                  if np.isfinite(row['rhat'])]
        return max(values) if values else float('nan')

    def __str__(self):
        lines = [str(r) for r in self.records]
        lines.append("")
        lines.append("{0} chain{1} of {2} from the {3}; {4:.1f} min".format(
            len(self.chains), "" if len(self.chains) == 1 else "s",
            self.spec.mcmc.n, "maximum-likelihood fit"
            if self.start_from == 'fit' else "pilot's best sample",
            self.seconds / 60.0))
        for note in self.notes:
            lines.append("  note: " + note)
        lines.append("")
        lines.append("  {0:9s} {1:>12s} {2:>11s} {3:>12s} {4:>12s} {5:>7s} "
                     "{6:>6s}".format("rate", "mean", "sd", "2.5%", "97.5%",
                                      "ESS", "R-hat"))
        for row in self.summary:
            lines.append("  {0:9s} {1:12.4g} {2:11.3g} {3:12.4g} {4:12.4g} "
                         "{5:7.0f} {6:6.3f}".format(
                             row['name'], row['mean'], row['sd'], row['q2.5'],
                             row['q97.5'], row['ess'], row['rhat']))
        failures = sum(c.nfailures for c in self.chains)
        acceptance = [float(np.mean(c.acceptance_rate()))
                      for c in self.chains]
        lines.append("")
        lines.append("acceptance after burn-in: {0}; likelihood failures: {1}"
                     .format(", ".join("{0:.2f}".format(a)
                                       for a in acceptance), failures))
        if self.max_rhat > 1.01:
            lines.append("R-hat above 1.01: the chains have not agreed yet. "
                         "Run them longer before using these numbers.")
        return "\n".join(lines)


def sample(spec, processes=None, verbose=False):
    """ Sample the posterior a specification describes (``[mcmc]``).

        :param spec: A :py:class:`~HJCFIT.likelihood.fitspec.FitSpec`.
        :param processes: Worker processes for the chains. None uses one per
          chain, up to the number of CPUs; 1 runs them one after another in
          this process.
        :param verbose: Print the stages as they happen.
        :returns: A :py:class:`SampleOutcome`.

        The C++ likelihood already spreads each evaluation over the cores
        with OpenMP, but past two or three threads it gains little. Measured
        on the AChR records of Epstein et al. (2016): 57 ms on one thread,
        36 ms on two, 29 ms on four. So parallel chains are run as separate
        processes, each limited to its share of the logical CPUs.

        What that buys depends on the machine. On a 4-core, 8-thread laptop,
        four CH82 chains ran 1.3-1.6 times faster in four processes than one
        after another in one. Fewer threads per process (one each) was
        slower there, because hyperthreading helps. More physical cores, or a
        likelihood that costs more per evaluation than CH82's 1-2 ms, gain
        more.

        Worker processes are *spawned*. As with any spawned process on
        Windows, a script that calls this with ``processes > 1`` must do so
        under ``if __name__ == '__main__':``. The ``hjcfit`` command and
        notebooks need nothing.
    """
    import time
    from .mcmc import gaussian_approximation, mwg_sample

    started = time.perf_counter()
    m = spec.mcmc
    records, mec, post = _posterior(spec)
    names = tuple(mec.get_free_parameter_names())
    notes = []
    fit = pilot = approximation = None

    if m.start == 'fit':
        if verbose:
            print("fitting, to start the chains at the maximum ...")
        fit = run(spec)
        start = np.asarray(fit.result.free_values, dtype=float)
        start_from = 'fit'
        if fit.limited:
            notes.append("the fit ended with {0} on a limit, which is also "
                         "the prior's edge".format(", ".join(
                             name for name, *_ in fit.limited)))
    else:
        if verbose:
            print("running the pilot from the guess ...")
        pilot = mwg_sample(post, np.asarray(mec.theta(), dtype=float),
                           n=m.pilot_n, burnin=m.pilot_n // 2,
                           rng=m.seed - 1)
        start, _ = pilot.mode()
        start_from = 'pilot'
    if not np.isfinite(post(start)):
        raise RunnerError("the posterior is zero at the starting point")

    covariance = None
    if m.covariance == 'hessian':
        try:
            approximation = gaussian_approximation(post, start)
            covariance = approximation.covariance
            if approximation.relative_error > 1e-2:
                notes.append("the Hessian at the start is poorly determined "
                             "(relative error {0:.1g}); a direction the "
                             "records barely constrain".format(
                                 approximation.relative_error))
        except ValueError as error:
            notes.append("no Hessian at the start ({0}); the chains begin "
                         "with the small isotropic step instead".format(
                             str(error).split(';')[0].split('.')[0]))

    starts = _starts(post, start, covariance, m.chains, m.seed)
    workers = 1 if processes == 1 else min(m.chains, processes
                                           or os.cpu_count() or 1)
    if verbose:
        print("{0} chain{1} of {2} iterations in {3} process{4} ...".format(
            m.chains, "" if m.chains == 1 else "s", m.n, workers,
            "" if workers == 1 else "es"))
    if workers == 1:
        chains = [_run_chain(spec, i, s, covariance)
                  for i, s in enumerate(starts)]
    else:
        chains = _parallel_chains(spec, starts, covariance, workers)

    return SampleOutcome(
        spec=spec, records=records, names=names, start=start,
        start_from=start_from, chains=chains,
        summary=_chain_summary(names, chains), fit=fit, pilot=pilot,
        approximation=approximation, notes=notes,
        seconds=time.perf_counter() - started, provenance=provenance())


def _parallel_chains(spec, starts, covariance, workers):
    """ Run the chains in *workers* spawned processes, each limited to its
        share of the OpenMP threads. The limit is set in the environment the
        processes inherit, before any of them loads the likelihood library,
        and restored afterwards. """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing

    threads = max(1, (os.cpu_count() or 1) // workers)
    previous = os.environ.get('OMP_NUM_THREADS')
    os.environ['OMP_NUM_THREADS'] = str(threads)
    try:
        context = multiprocessing.get_context('spawn')
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=context) as pool:
            futures = [pool.submit(_run_chain, spec, i, s, covariance)
                       for i, s in enumerate(starts)]
            return [f.result() for f in futures]
    finally:
        if previous is None:
            os.environ.pop('OMP_NUM_THREADS', None)
        else:
            os.environ['OMP_NUM_THREADS'] = previous


def samples_as_dict(outcome):
    """ A sampling outcome as JSON-ready data: the specification, the start,
        the summary and the provenance. The chains themselves go in ``.npz``
        files beside it (:py:func:`write_samples`). """
    approx = outcome.approximation
    return {
        'hjcfit_samples': 1,
        'spec': outcome.spec.as_dict(),
        'free_names': list(outcome.names),
        'start': {'from': outcome.start_from,
                  'rates': [float(v) for v in outcome.start]},
        'fit': None if outcome.fit is None else result_as_dict(
            outcome.fit)['fit'],
        'hessian_sd': None if approx is None else [float(v)
                                                   for v in approx.sd],
        'summary': outcome.summary,
        'chains': [{'n': c.n, 'burnin': c.burnin, 'seconds': c.seconds,
                    'acceptance': np.atleast_1d(c.acceptance_rate()).tolist(),
                    'nfailures': c.nfailures} for c in outcome.chains],
        'notes': list(outcome.notes),
        'seconds': outcome.seconds,
        'provenance': outcome.provenance,
    }


def write_samples(prefix, outcome):
    """ Write ``<prefix>.json`` (:py:func:`samples_as_dict`) and one
        ``<prefix>_chain<i>.npz`` per chain
        (:py:meth:`~HJCFIT.likelihood.mcmc.Chain.save`).

        :returns: The paths written.
    """
    prefix = str(prefix)
    if prefix.endswith('.json'):
        prefix = prefix[:-5]
    paths = [prefix + '.json']
    with open(paths[0], 'w', encoding='utf-8', newline='\n') as handle:
        json.dump(samples_as_dict(outcome), handle, indent=2)
        handle.write("\n")
    for i, chain in enumerate(outcome.chains):
        path = '{0}_chain{1}.npz'.format(prefix, i)
        chain.save(path)
        paths.append(path)
    return paths
