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

""" The ``hjcfit`` command: a fit from a file, with no Python written.

    ::

        hjcfit template -o ch82.toml     a specification to edit
        hjcfit check ch82.toml           say what it would do; run nothing
        hjcfit fit ch82.toml -o out.json run it, and keep the result

    Also reachable as ``python -m HJCFIT.likelihood.cli`` when the console
    script is not on the path, which is the usual state of affairs inside a
    conda environment on Windows.

    ``check`` is the one worth the habit. It reads the specification, finds the
    records, builds the mechanism with every constraint applied, and evaluates
    the likelihood once at the initial guess -- everything a fit does except
    the search. A misspelled rate name, a critical time below the dead time, a
    mechanism with a free-parameter count nobody expected: all of them surface
    in under a second, rather than twenty minutes into a search or, worse, in
    a set of estimates that look plausible.

    Every error this prints is one line beginning with ``hjcfit:`` and the exit
    status is 1, because this is a command that belongs in scripts.
"""
__docformat__ = "restructuredtext en"
__all__ = ['main']

import argparse
import sys

from .fitspec import TEMPLATE, FitSpec, SpecError


def _add_spec_argument(parser):
    parser.add_argument('spec', metavar='SPEC.toml',
                        help='the fit specification')


def _parser():
    parser = argparse.ArgumentParser(
        prog='hjcfit',
        description='Maximum likelihood fitting of ion-channel mechanisms, '
                    'driven by a specification file.',
        epilog='Start with: hjcfit template -o my-fit.toml')
    parser.add_argument('--version', action='store_true',
                        help='print the versions of everything involved')
    sub = parser.add_subparsers(dest='command', metavar='COMMAND')

    template = sub.add_parser(
        'template', help='write a specification to start from')
    template.add_argument('-o', '--output', metavar='FILE', default=None,
                          help='write here instead of to standard output')
    template.add_argument('-f', '--force', action='store_true',
                          help='overwrite FILE if it exists')

    check = sub.add_parser(
        'check', help='say what a specification would do, and do nothing')
    _add_spec_argument(check)

    fit = sub.add_parser('fit', help='run a specification')
    _add_spec_argument(fit)
    fit.add_argument('-o', '--output', metavar='FILE', default=None,
                     help='write the result, with its provenance, as JSON')
    fit.add_argument('-q', '--quiet', action='store_true',
                     help='print the result and nothing before it')

    show = sub.add_parser(
        'show', help='read a specification and write it back, normalised')
    _add_spec_argument(show)
    return parser


def _fail(message):
    """ One line on stderr, after whatever has already been printed.

        The flush is not decoration. Python buffers stdout when it is not a
        terminal, so without it a message written to the unbuffered stderr
        arrives before the output it is about, and a command that printed
        three useful lines and then failed reads as though it failed first.
    """
    sys.stdout.flush()
    sys.stderr.write('hjcfit: {0}\n'.format(message))


def _load(path):
    """ Read a specification, or fail with one line. """
    try:
        return FitSpec.from_toml(path)
    except FileNotFoundError:
        raise SpecError('no such file: {0}'.format(path))
    except OSError as error:
        raise SpecError('cannot read {0}: {1}'.format(path, error))


def _versions():
    from .runner import provenance
    prov = provenance()
    lines = ['python {0} on {1}'.format(prov['python'], prov['platform'])]
    for package, version in prov['versions'].items():
        lines.append('{0:12s} {1}'.format(
            package, version if version else 'not installed'))
    return "\n".join(lines)


def _template(args):
    if args.output is None:
        sys.stdout.write(TEMPLATE)
        return 0
    import os
    if os.path.exists(args.output) and not args.force:
        raise SpecError('{0} exists; pass --force to overwrite it'
                        .format(args.output))
    with open(args.output, 'w', encoding='utf-8', newline='\n') as handle:
        handle.write(TEMPLATE)
    print('wrote {0}'.format(args.output))
    print('edit it, then: hjcfit check {0}'.format(args.output))
    return 0


def _check(args):
    """ Everything a fit does except the search. """
    from .runner import build_mechanism, load_records
    from .fitting import HJCFitter

    spec = _load(args.spec)
    print(spec)
    print('')

    records = load_records(spec)
    for record in records:
        print(record)

    mec = build_mechanism(spec)
    free = list(mec.get_free_parameter_names())
    print('')
    print('{0}: {1} states, {2} open, {3} rate constants'.format(
        getattr(mec, 'mtitle', '') or 'mechanism', mec.k, mec.kA,
        len(mec.Rates)))
    print('{0} free: {1}'.format(len(free), ", ".join(free)))
    held = [r.name for r in mec.Rates if not r.is_free]
    if held:
        print('{0} fixed or constrained: {1}'.format(
            len(held), ", ".join(held)))

    # Evaluating the likelihood once is the whole point of check: it is the
    # only thing here that proves the records, the mechanism and the dead time
    # are mutually possible. A guess at which the likelihood cannot be
    # computed is a guess no search can start from.
    fitter = HJCFitter(mec, records, log_params=spec.search.log_params)
    print('')
    print('log10L at the guess: {0:.4f}'.format(fitter.log10_likelihood()))
    print('nothing was fitted; run: hjcfit fit {0}'.format(args.spec))
    return 0


def _fit(args):
    from .runner import run, write_result

    spec = _load(args.spec)
    if not args.quiet:
        print(spec)
        print('')
    outcome = run(spec, verbose=not args.quiet)
    if not args.quiet:
        print('')
    print(outcome)
    if args.output:
        write_result(args.output, outcome)
        print('')
        print('wrote {0}'.format(args.output))
    # A fit whose optimiser did not converge has still produced a number, and
    # the number is printed; the exit status is what a script reads, so it says
    # so. Rates sitting on a limit are not a failure of the search -- they are
    # a statement about the model -- and are reported without changing it.
    return 0 if outcome.result.success else 2


def _show(args):
    sys.stdout.write(_load(args.spec).to_toml())
    return 0


def main(argv=None):
    """ The ``hjcfit`` command.

        :param argv: Arguments, without the program name. None takes
            ``sys.argv``.
        :returns: 0 on success, 1 on anything this module can explain, and 2
            from ``fit`` when the search did not converge.
    """
    parser = _parser()
    args = parser.parse_args(argv)

    if args.version:
        print(_versions())
        return 0
    if args.command is None:
        parser.print_help()
        return 0

    handlers = {'template': _template, 'check': _check, 'fit': _fit,
                'show': _show}
    try:
        return handlers[args.command](args)
    except SpecError as error:
        _fail(error)
        return 1
    except Exception as error:
        # RunnerError and anything dcio or scalcs raises. One line, named, and
        # no traceback: a traceback out of a command is noise to the person who
        # mistyped a rate name. --debug is deliberately not offered, because
        # `python -c "from HJCFIT.likelihood.runner import run"` is the same
        # thing with a real traceback.
        from .runner import RunnerError
        if not isinstance(error, (RunnerError, OSError, ValueError,
                                  ArithmeticError, KeyError, TypeError)):
            raise
        _fail('{0}: {1}'.format(type(error).__name__, error))
        return 1


if __name__ == '__main__':
    sys.exit(main())
