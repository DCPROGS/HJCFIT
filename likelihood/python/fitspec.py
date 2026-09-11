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

""" What a fit is, written down.

    A fit needs three things said about it: which records, which mechanism and
    how to search. This module is those three things as data, and **nothing
    else**. It reads no records, builds no mechanism and runs no search, and it
    imports nothing but the standard library.

    That separation is the point. The same description drives
    :py:mod:`~HJCFIT.likelihood.runner`, the ``hjcfit`` command, the notebook
    template in ``examples/`` and -- if one is ever built -- a desktop
    interface, because a user interface is a way of editing a fit
    specification. It is also what makes a fit reproducible by somebody else:
    a file they can read, diff and keep beside the result.

    A specification is a file::

        title = "CH82 sample record at 100 nM"

        [[data]]
        record = "CH82"
        conc   = 100e-9
        tres   = 100e-6
        tcrit  = 4e-3

        [mechanism]
        sample = "CH82"

        [search]
        method = "simplex"

    and :py:meth:`FitSpec.from_toml` reads it back.

    **Why TOML rather than YAML.** A specification is mostly rate constants,
    and rate constants get written ``1e8``. PyYAML implements YAML 1.1, whose
    resolver requires an exponent to carry a sign: ``yaml.safe_load`` reads
    ``1e8`` as the *string* ``'1e8'``, and so does ``1.0e8``, and so does
    ``1e+8``. Only ``1.0e+8`` becomes a float. A rate constant silently
    arriving as a string is the worst failure this file could have. TOML has
    one number syntax, accepts every form above, and is in the standard
    library from Python 3.11 -- ``tomli`` covers 3.10 and is in the
    ``[fitting]`` extra. The ``.yaml`` file in ``scalcs/samples`` is a third
    thing again, a pickled Python object graph needing ``unsafe_load``, so
    YAML in this stack already means something other than a document a person
    edits.

    **Rates are named, never numbered.** Every reference to a rate constant
    here is by name. SCALCS addresses rates by index into ``mec.Rates`` --
    ``set_mr(True, 5, 0)`` -- and an index is what users get wrong silently,
    because rate 5 of a mechanism is whatever the sample happened to list
    fifth. :py:mod:`~HJCFIT.likelihood.runner` resolves names to indices in
    one place, and says which names exist when one does not.
"""
__docformat__ = "restructuredtext en"
__all__ = ['DataSpec', 'MechanismSpec', 'SearchSpec', 'FitSpec',
           'SpecError', 'TEMPLATE', 'load_toml_bytes']

from dataclasses import dataclass, field


class SpecError(ValueError):
    """ A specification that cannot mean anything.

        Raised by the ``from_dict`` constructors and by ``validate``. Every
        message names the field it is about.
    """


#: Accepted values of :py:attr:`DataSpec.vectors`.
VECTORS = ('chs', 'equilibrium')

#: Accepted values of :py:attr:`SearchSpec.method`.
METHODS = ('simplex', 'scipy')


def _number(value, where):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SpecError("{0}: expected a number, got {1!r}".format(where, value))
    return float(value)


def _named_numbers(mapping, where):
    """ A ``{rate name: value}`` table, checked and copied. """
    if mapping is None:
        return {}
    if not isinstance(mapping, dict):
        raise SpecError("{0}: expected a table of 'rate name' = value, got "
                        "{1!r}".format(where, mapping))
    return {str(k): _number(v, "{0}.{1}".format(where, k))
            for k, v in mapping.items()}


def _unexpected(given, known, where):
    """ Refuse a key nobody reads.

        A misspelled key in a specification is worse than a missing one: the
        fit runs, silently, without the setting the file says it has. Every
        table here is closed for that reason.
    """
    extra = [k for k in given if k not in known]
    if extra:
        raise SpecError("{0}: unknown key{1} {2}; known keys are {3}".format(
            where, "" if len(extra) == 1 else "s",
            ", ".join(repr(k) for k in sorted(extra)),
            ", ".join(sorted(known))))


@dataclass(frozen=True)
class DataSpec:
    """ One record, and how it is presented to the likelihood.

        :param record:
          An ``*.scn`` file, or the name of a sample record shipped with
          HJCFIT (``"CH82"``, ``"CO"``, ``"CCO"``).
        :param conc: Agonist concentration [M].
        :param tres: Dead time to impose [s].
        :param tcrit:
          Critical shut time dividing the record into groups [s]. Omitted, the
          whole record is fitted as a single group, which assumes one channel
          throughout.
        :param vectors:
          ``"chs"`` for the CHS vectors of Colquhoun, Hawkes & Srodzinski
          (1996), ``"equilibrium"`` for the equilibrium vectors of Colquhoun &
          Hawkes (1982). Groups cut at a critical time chosen to separate
          activations of one channel want CHS; clusters at a concentration
          high enough to desensitise want equilibrium, and so does a whole
          record. Defaults to CHS when *tcrit* is given and to equilibrium
          when it is not.

        *tcrit* and *vectors* are two fields rather than one signed number.
        Elsewhere in the stack a negative critical time is a flag selecting
        equilibrium vectors while its magnitude is still the time that divides
        the record -- :py:func:`HJCFIT.read_idealized_bursts` documents that,
        and :py:class:`~HJCFIT.likelihood.fitting.Record` takes the magnitude
        or None. Two fields say the same thing without a sign to remember, and
        they make the combination that means nothing -- one group, but vectors
        defined *between* groups -- something :py:meth:`validate` can refuse.
    """

    record: str
    conc: float
    tres: float
    tcrit: float = None
    vectors: str = 'chs'

    KEYS = ('record', 'conc', 'tres', 'tcrit', 'vectors')

    @classmethod
    def from_dict(cls, d, where='data'):
        """ Build one from a parsed table, checking as it goes. """
        _unexpected(d, cls.KEYS, where)
        for required in ('record', 'conc', 'tres'):
            if required not in d:
                raise SpecError("{0}: {1} is required".format(where, required))
        tcrit = d.get('tcrit')
        return cls(
            record=str(d['record']),
            conc=_number(d['conc'], where + '.conc'),
            tres=_number(d['tres'], where + '.tres'),
            tcrit=None if tcrit is None else _number(tcrit, where + '.tcrit'),
            vectors=str(d.get('vectors',
                              'chs' if tcrit is not None else 'equilibrium')),
        ).validate(where)

    def validate(self, where='data'):
        """ Raise :py:exc:`SpecError` unless this record can be fitted. """
        if self.conc < 0.0:
            raise SpecError("{0}.conc: {1} is negative".format(where, self.conc))
        if self.tres <= 0.0:
            raise SpecError("{0}.tres: must be positive, not {1}"
                            .format(where, self.tres))
        if self.vectors not in VECTORS:
            raise SpecError("{0}.vectors: must be one of {1}, not {2!r}".format(
                where, " or ".join(VECTORS), self.vectors))
        if self.tcrit is None:
            if self.vectors != 'equilibrium':
                raise SpecError(
                    "{0}: with no tcrit the whole record is one group, and "
                    "CHS vectors are defined between groups; either give "
                    "tcrit or set vectors to 'equilibrium'".format(where))
        else:
            if self.tcrit <= 0.0:
                raise SpecError(
                    "{0}.tcrit: must be positive. The sign is not a flag "
                    "here -- vectors is.".format(where))
            if self.tcrit <= self.tres:
                raise SpecError(
                    "{0}: tcrit ({1:g} s) must exceed tres ({2:g} s); a "
                    "critical shut time below the dead time divides nothing"
                    .format(where, self.tcrit, self.tres))
        return self

    def as_dict(self):
        d = {'record': self.record, 'conc': self.conc, 'tres': self.tres}
        if self.tcrit is not None:
            d['tcrit'] = self.tcrit
        d['vectors'] = self.vectors
        return d

    def __str__(self):
        conc = ("{0:g} nM".format(self.conc * 1e9) if self.conc < 1e-6
                else "{0:g} uM".format(self.conc * 1e6))
        groups = ("whole record as one group" if self.tcrit is None else
                  "groups at tcrit {0:g} ms".format(self.tcrit * 1e3))
        return "{0} at {1}: {2:g} us dead time, {3}, {4} vectors".format(
            self.record, conc, self.tres * 1e6, groups, self.vectors)


@dataclass(frozen=True)
class MechanismSpec:
    """ Which mechanism, with which constraints and which initial guess.

        :param sample:
          Name of a factory in :py:mod:`scalcs.samples.samples`, for example
          ``"CH82"`` or ``"load_AChR_diamond_independent_binding"``. Exactly
          one of this and *mec_file*.
        :param mec_file:
          A DCprogs ``*.mec`` file, read with
          :py:func:`scalcs.scalcsio.mec_load`.
        :param mec_number:
          Which mechanism in that file, by the sequence number
          ``scalcsio.mec_get_list`` reports. Needed only when the file holds
          more than one.
        :param rates:
          Initial guess, by rate name. Names left out keep the value the
          sample or the file carries, so a guess can be one line.
        :param fixed:
          Rate constants held at a stated value throughout the fit.
        :param limits:
          ``'name' = [lower, upper]``. A rate driven outside its limits is
          **reset** to the limit before the likelihood sees it, which is what
          HJCFIT did; see
          :py:class:`~HJCFIT.likelihood.fitting.HJCFitter`. Nothing announces
          that it happened, so a rate that ends against a limit has to be
          noticed by comparing the result with what is set here --
          :py:mod:`~HJCFIT.likelihood.runner` does that comparison and says so.
        :param mr:
          Name of the rate constant computed from microscopic reversibility
          rather than fitted.
        :param mr_cycle:
          Which cycle, when the mechanism has more than one. Cycles have no
          names in SCALCS, so this one is an index.
        :param ec50:
          ``{rate = 'name', value = EC50 in M}``. The option of Colquhoun,
          Hatton & Hawkes (2003) p. 702: an independently measured EC50 is
          supplied and one rate constant is computed from it and all the
          others at each iteration, which removes a free parameter. An EC50
          given wrongly biases the estimates rather than failing -- that
          paper's Figures 9 and 10 are what a factor of two either way does --
          so what was used belongs in the specification, and in the result.
        :param nfree:
          Free parameters expected once every constraint above is applied.
          Not an input: the runner asserts it. A change in how SCALCS handles
          constraints then shows up as a refusal to start rather than as a
          puzzling set of estimates a day later.
    """

    sample: str = None
    mec_file: str = None
    mec_number: int = None
    rates: dict = field(default_factory=dict)
    fixed: dict = field(default_factory=dict)
    limits: dict = field(default_factory=dict)
    mr: str = None
    mr_cycle: int = 0
    ec50: dict = None
    nfree: int = None

    KEYS = ('sample', 'mec_file', 'mec_number', 'rates', 'fixed', 'limits',
            'mr', 'mr_cycle', 'ec50', 'nfree')

    @classmethod
    def from_dict(cls, d, where='mechanism'):
        """ Build one from a parsed table, checking as it goes. """
        _unexpected(d, cls.KEYS, where)

        limits = {}
        for name, pair in (d.get('limits') or {}).items():
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                raise SpecError("{0}.limits.{1}: expected [lower, upper], got "
                                "{2!r}".format(where, name, pair))
            lower = _number(pair[0], "{0}.limits.{1}[0]".format(where, name))
            upper = _number(pair[1], "{0}.limits.{1}[1]".format(where, name))
            if not lower < upper:
                raise SpecError("{0}.limits.{1}: lower ({2:g}) must be below "
                                "upper ({3:g})".format(where, name, lower,
                                                       upper))
            limits[str(name)] = (lower, upper)

        ec50 = d.get('ec50')
        if ec50 is not None:
            if not isinstance(ec50, dict):
                raise SpecError("{0}.ec50: expected a table with rate and "
                                "value, got {1!r}".format(where, ec50))
            _unexpected(ec50, ('rate', 'value'), where + '.ec50')
            for required in ('rate', 'value'):
                if required not in ec50:
                    raise SpecError("{0}.ec50: {1} is required"
                                    .format(where, required))
            ec50 = {'rate': str(ec50['rate']),
                    'value': _number(ec50['value'], where + '.ec50.value')}

        return cls(
            sample=None if d.get('sample') is None else str(d['sample']),
            mec_file=None if d.get('mec_file') is None else str(d['mec_file']),
            mec_number=(None if d.get('mec_number') is None
                        else int(d['mec_number'])),
            rates=_named_numbers(d.get('rates'), where + '.rates'),
            fixed=_named_numbers(d.get('fixed'), where + '.fixed'),
            limits=limits,
            mr=None if d.get('mr') is None else str(d['mr']),
            mr_cycle=int(d.get('mr_cycle', 0)),
            ec50=ec50,
            nfree=None if d.get('nfree') is None else int(d['nfree']),
        ).validate(where)

    def validate(self, where='mechanism'):
        """ Raise :py:exc:`SpecError` unless this mechanism can be built. """
        if (self.sample is None) == (self.mec_file is None):
            raise SpecError("{0}: give exactly one of sample and mec_file"
                            .format(where))
        if self.mec_file is None and self.mec_number is not None:
            raise SpecError("{0}.mec_number: means nothing without mec_file"
                            .format(where))
        both = sorted(set(self.rates) & set(self.fixed))
        if both:
            raise SpecError(
                "{0}: {1} {2} in both rates and fixed. A fixed rate takes its "
                "value from fixed; naming it in rates as well says the same "
                "thing twice and invites the two to disagree.".format(
                    where, ", ".join(both),
                    "is" if len(both) == 1 else "are"))
        if self.ec50 is not None:
            if self.ec50['value'] <= 0.0:
                raise SpecError("{0}.ec50.value: must be positive, not {1}"
                                .format(where, self.ec50['value']))
            if self.ec50['rate'] in self.fixed:
                raise SpecError(
                    "{0}: {1} is both fixed and computed from the EC50; it "
                    "cannot be both".format(where, self.ec50['rate']))
        if self.nfree is not None and self.nfree < 1:
            raise SpecError("{0}.nfree: must be at least 1, not {1}"
                            .format(where, self.nfree))
        return self

    def as_dict(self):
        d = {}
        for name in ('sample', 'mec_file', 'mec_number', 'mr', 'nfree'):
            value = getattr(self, name)
            if value is not None:
                d[name] = value
        if self.mr is not None and self.mr_cycle:
            d['mr_cycle'] = self.mr_cycle
        for name in ('rates', 'fixed'):
            if getattr(self, name):
                d[name] = dict(getattr(self, name))
        if self.limits:
            d['limits'] = {k: list(v) for k, v in self.limits.items()}
        if self.ec50 is not None:
            d['ec50'] = dict(self.ec50)
        return d


@dataclass(frozen=True)
class SearchSpec:
    """ How to search, and how hard.

        :param method:
          ``"simplex"`` is
          :py:func:`~HJCFIT.likelihood.optimization.simplex_hjc`, HJCFIT's own
          and the search that produced every published result. ``"scipy"`` is
          SciPy's Nelder-Mead with the restart loop below; reach for it when
          the starting point is poor, because a regular simplex needs every
          one of its vertices to be evaluable.
        :param log_params:
          Search the logarithms of the rate constants. HJCFIT's default, three
          to four times faster, and it cannot produce a negative rate. It
          belongs in the specification rather than in the driver because it
          changes which maximum a search falls into: Colquhoun, Hatton &
          Hawkes (2003) say of their Figure 2 fits that the rate constants
          themselves were the free parameters, and the same set of fits run
          over logarithms reaches the second solution at a quite different
          rate.
        :param restarts:
          ``"scipy"`` only: extra runs from the previous solution, stopping
          once one gains less than *fatol*. A single SciPy pass leaves fits
          stranded between maxima. Measured over the twelve stranded fits of
          one reproduction scenario, restarting gained 106.3
          log\\ :sub:`10` units in total, resolved five of them onto a real
          maximum, and never lost any.
        :param maxfev: Evaluation budget; both searches.
        :param maxiter: ``"scipy"`` only.
        :param xatol: ``"scipy"`` only.
        :param fatol: ``"scipy"`` only.
        :param options:
          Passed straight to whichever search is chosen, for deliberate
          departures from its own defaults.
    """

    method: str = 'simplex'
    log_params: bool = True
    restarts: int = 0
    maxfev: int = 20000
    maxiter: int = 20000
    xatol: float = 1e-4
    fatol: float = 1e-4
    options: dict = field(default_factory=dict)

    KEYS = ('method', 'log_params', 'restarts', 'maxfev', 'maxiter',
            'xatol', 'fatol', 'options')

    @classmethod
    def from_dict(cls, d, where='search'):
        """ Build one from a parsed table, checking as it goes. """
        _unexpected(d, cls.KEYS, where)
        options = d.get('options') or {}
        if not isinstance(options, dict):
            raise SpecError("{0}.options: expected a table, got {1!r}"
                            .format(where, options))
        return cls(
            method=str(d.get('method', 'simplex')),
            log_params=bool(d.get('log_params', True)),
            restarts=int(d.get('restarts', 0)),
            maxfev=int(d.get('maxfev', 20000)),
            maxiter=int(d.get('maxiter', 20000)),
            xatol=_number(d.get('xatol', 1e-4), where + '.xatol'),
            fatol=_number(d.get('fatol', 1e-4), where + '.fatol'),
            options=dict(options),
        ).validate(where)

    def validate(self, where='search'):
        """ Raise :py:exc:`SpecError` unless this search can be run. """
        if self.method not in METHODS:
            raise SpecError("{0}.method: must be one of {1}, not {2!r}".format(
                where, " or ".join(METHODS), self.method))
        if self.restarts and self.method != 'scipy':
            raise SpecError(
                "{0}.restarts: only the scipy search takes restarts. "
                "simplex_hjc has its own restart rule, capped by nresmax, "
                "which goes in options.".format(where))
        if self.restarts < 0:
            raise SpecError("{0}.restarts: cannot be negative, got {1}"
                            .format(where, self.restarts))
        for name in ('maxfev', 'maxiter'):
            if getattr(self, name) < 1:
                raise SpecError("{0}.{1}: must be at least 1, not {2}".format(
                    where, name, getattr(self, name)))
        return self

    def as_dict(self):
        d = {'method': self.method, 'log_params': self.log_params,
             'maxfev': self.maxfev}
        if self.method == 'scipy':
            d.update(restarts=self.restarts, maxiter=self.maxiter,
                     xatol=self.xatol, fatol=self.fatol)
        if self.options:
            d['options'] = dict(self.options)
        return d


@dataclass(frozen=True)
class FitSpec:
    """ A whole fit, as data.

        :param data:
          One :py:class:`DataSpec` per record. More than one means several
          records contributing to a single likelihood: the rate constants are
          shared and the concentrations are not.
        :param mechanism: A :py:class:`MechanismSpec`.
        :param search: A :py:class:`SearchSpec`.
        :param title: One line for the person reading the result.

        Read one with :py:meth:`from_toml`, write one with
        :py:meth:`to_toml`, and run one with
        :py:func:`HJCFIT.likelihood.runner.run`.
    """

    data: tuple
    mechanism: MechanismSpec = field(default_factory=MechanismSpec)
    search: SearchSpec = field(default_factory=SearchSpec)
    title: str = ''

    KEYS = ('data', 'mechanism', 'search', 'title')

    # -- reading -----------------------------------------------------------

    @classmethod
    def from_dict(cls, d, where='spec'):
        """ Build one from a parsed document.

            :raises SpecError: naming the field, on anything unusable.
        """
        if not isinstance(d, dict):
            raise SpecError("{0}: expected a table, got {1!r}".format(where, d))
        _unexpected(d, cls.KEYS, where)
        records = d.get('data')
        if not records:
            raise SpecError(
                "{0}: at least one [[data]] section is needed -- a fit "
                "without a record is not a fit".format(where))
        if isinstance(records, dict):       # a single [data] rather than [[data]]
            records = [records]
        if not isinstance(records, (list, tuple)):
            raise SpecError("{0}.data: expected one or more [[data]] "
                            "sections, got {1!r}".format(where, records))
        return cls(
            data=tuple(DataSpec.from_dict(r, "data[{0}]".format(i))
                       for i, r in enumerate(records)),
            mechanism=MechanismSpec.from_dict(d.get('mechanism') or {}),
            search=SearchSpec.from_dict(d.get('search') or {}),
            title=str(d.get('title', '')),
        ).validate(where)

    @classmethod
    def from_toml(cls, path):
        """ Read a specification from a TOML file.

            :param path: The file. Read as bytes, which is what TOML wants:
                the format is UTF-8 by definition, so letting the platform
                choose an encoding could only get it wrong.
            :raises SpecError: on a file TOML cannot parse, or a field this
                module cannot use. Both carry the file name.
        """
        with open(path, 'rb') as handle:
            raw = handle.read()
        try:
            return cls.from_dict(load_toml_bytes(raw))
        except SpecError as error:
            raise SpecError("{0}: {1}".format(path, error))

    def validate(self, where='spec'):
        """ Raise :py:exc:`SpecError` unless the whole thing hangs together. """
        for i, record in enumerate(self.data):
            record.validate("data[{0}]".format(i))
        self.mechanism.validate()
        self.search.validate()
        # Two records at the same concentration are not an error -- two patches
        # at one concentration is an ordinary experiment, and the reproduction
        # of Colquhoun, Hatton & Hawkes (2003) fits three concentrations of
        # which none repeats. What would be an error is a mechanism with no
        # concentration dependence fitted to several concentrations, and that
        # cannot be seen from here: it needs the mechanism built.
        return self

    # -- writing -----------------------------------------------------------

    def as_dict(self):
        """ A plain dictionary, ready for :py:func:`to_toml` or for JSON.

            Defaults are left out, so a specification written back is the
            shortest file that means the same thing.
        """
        d = {}
        if self.title:
            d['title'] = self.title
        d['data'] = [r.as_dict() for r in self.data]
        d['mechanism'] = self.mechanism.as_dict()
        d['search'] = self.search.as_dict()
        return d

    def to_toml(self):
        """ The specification as TOML text.

            Written by hand rather than with a library: ``tomllib`` reads and
            does not write, and the alternative is a dependency to emit forty
            lines of a format with four value types in it.
        """
        lines = []
        if self.title:
            lines.append('title = {0}'.format(_toml_value(self.title)))
            lines.append('')
        for record in self.data:
            lines.append('[[data]]')
            lines.extend(_toml_table(record.as_dict()))
            lines.append('')
        for name in ('mechanism', 'search'):
            table = getattr(self, name).as_dict()
            lines.append('[{0}]'.format(name))
            lines.extend(_toml_table(table, name))
            lines.append('')
        return "\n".join(lines).rstrip() + "\n"

    def write_toml(self, path):
        """ Write :py:meth:`to_toml` to *path*, as UTF-8. """
        with open(path, 'w', encoding='utf-8', newline='\n') as handle:
            handle.write(self.to_toml())

    # -- looking at it -----------------------------------------------------

    def __str__(self):
        lines = [self.title or 'a fit']
        for record in self.data:
            lines.append('  ' + str(record))
        mech = self.mechanism
        source = (mech.sample if mech.sample is not None else
                  '{0}{1}'.format(mech.mec_file,
                                  '' if mech.mec_number is None
                                  else ' #{0}'.format(mech.mec_number)))
        lines.append('  mechanism: {0}'.format(source))
        for label, table in (('guess', mech.rates), ('fixed', mech.fixed)):
            if table:
                lines.append('    {0}: {1}'.format(label, ", ".join(
                    '{0} = {1:g}'.format(k, v)
                    for k, v in sorted(table.items()))))
        if mech.mr is not None:
            lines.append('    {0} from microscopic reversibility'
                         .format(mech.mr))
        if mech.ec50 is not None:
            lines.append('    {0} from an EC50 of {1:g} uM'.format(
                mech.ec50['rate'], mech.ec50['value'] * 1e6))
        if mech.limits:
            lines.append('    limits on {0}'.format(
                ", ".join(sorted(mech.limits))))
        lines.append('  search: {0}, over {1}, up to {2} evaluations'.format(
            self.search.method,
            'logarithms' if self.search.log_params else 'rate constants',
            self.search.maxfev))
        return "\n".join(lines)


# ---------------------------------------------------------------- TOML I/O

def load_toml_bytes(raw):
    """ Parse TOML from bytes, whichever parser this Python has.

        ``tomllib`` is in the standard library from 3.11; ``tomli`` is the
        same parser by the same author and covers 3.10, and is in the
        ``[fitting]`` extra for that reason.

        :raises SpecError: on a file TOML cannot parse, or when neither
            parser is importable.
    """
    try:
        import tomllib as toml
    except ImportError:                      # Python 3.10
        try:
            import tomli as toml
        except ImportError:
            raise SpecError(
                "reading a specification needs a TOML parser. Python 3.11 and "
                "later have one in the standard library; on 3.10, install "
                "tomli -- pip install 'hjcfit[fitting]' does that.")
    try:
        return toml.loads(raw.decode('utf-8'))
    except UnicodeDecodeError as error:
        raise SpecError("not UTF-8, which TOML must be: {0}".format(error))
    except toml.TOMLDecodeError as error:
        raise SpecError("not valid TOML: {0}".format(error))


def _toml_value(value):
    """ One TOML scalar, or an array of them. """
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, int):
        return repr(value)
    if isinstance(value, float):
        # repr round-trips, and TOML accepts every form it produces except the
        # bare infinities and nan, which cannot be a rate constant or a time.
        text = repr(value)
        if text in ('inf', '-inf', 'nan'):
            raise SpecError("{0} cannot go in a specification".format(text))
        return text
    if isinstance(value, (list, tuple)):
        return '[' + ", ".join(_toml_value(v) for v in value) + ']'
    text = str(value)
    return '"' + text.replace('\\', '\\\\').replace('"', '\\"') + '"'


def _toml_key(name):
    """ A bare key where TOML allows one, a quoted key otherwise.

        Rate constants are called things like ``k(-1)`` and ``2k(+2)``, which
        are not bare keys, so most of these come back quoted.
    """
    if name and all(c.isalnum() or c in '-_' for c in name):
        return name
    return _toml_value(name)


def _toml_table(table, prefix=''):
    """ One table as lines, with its sub-tables after its own keys.

        TOML reads every key after a ``[header]`` as belonging to it, so a
        scalar written below a sub-table would silently join the sub-table.
        Scalars therefore come first, which is not a matter of taste.
    """
    lines = [
        '{0} = {1}'.format(_toml_key(k), _toml_value(v))
        for k, v in table.items() if not isinstance(v, dict)
    ]
    for key, value in table.items():
        if isinstance(value, dict):
            if not value:
                continue
            header = '{0}.{1}'.format(prefix, key) if prefix else key
            lines.append('')
            lines.append('[{0}]'.format(header))
            lines.extend(_toml_table(value, header))
    return lines


#: A specification to start from, as written by ``hjcfit template``. It fits
#: the CH82 sample record that ships with HJCFIT, so it runs as it stands.
TEMPLATE = '''\
# A fit specification. Everything a fit needs, and nothing that runs.
#
#   hjcfit check  this-file.toml     say what it would do, and do nothing
#   hjcfit fit    this-file.toml     run it
#
# As it stands this fits the CH82 sample record that ships with HJCFIT, so it
# works before you change anything. Then point `record` at your own .scn file.

title = "CH82 sample record at 100 nM"

# One [[data]] section per record. Several means several records contributing
# to one likelihood: the rate constants are shared, the concentrations are not.
[[data]]
record = "CH82"        # an .scn file, or a sample: CH82, CO, CCO
conc = 1e-07           # M
tres = 0.0001          # dead time to impose, s
tcrit = 0.004          # critical shut time dividing the record into groups, s
vectors = "chs"        # "chs" between groups of one channel's activations;
                       # "equilibrium" for clusters, or for a whole record
                       # (omit tcrit for that, and vectors must be equilibrium)

[mechanism]
sample = "CH82"        # a factory in scalcs.samples.samples
# mec_file = "mechs.mec"   # or a DCprogs .mec file, instead of sample
# mec_number = 1           # which mechanism in it
# mr = "k(-1)"             # a rate computed from microscopic reversibility
# nfree = 8                # free parameters expected; asserted, not applied

# The initial guess, by rate name. Names left out keep the sample's own value,
# so this can be one line. `hjcfit check` prints every name the mechanism has.
[mechanism.rates]
beta1 = 15.0
beta2 = 15000.0

# Held at this value throughout.
# [mechanism.fixed]
# "k(-1)" = 2000.0

# A rate driven outside its limits is reset to the limit before the likelihood
# sees it, and nothing announces it. `hjcfit fit` reports any rate that ends
# against one, because such a rate is not an estimate.
# [mechanism.limits]
# beta1 = [0.1, 100000.0]

# One rate computed from an independently measured EC50 at every iteration,
# instead of being fitted (Colquhoun, Hatton & Hawkes 2003, p. 702).
# [mechanism.ec50]
# rate = "k(+2)"
# value = 3.3e-06

[search]
method = "simplex"     # HJCFIT's own; "scipy" is Nelder-Mead with restarts
log_params = true      # search the logarithms: faster, and cannot go negative
maxfev = 20000
'''
