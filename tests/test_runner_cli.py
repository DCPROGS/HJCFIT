"""Carrying out a specification: the runner, the hjcfit command, the notebook.

Needs the [fitting] extra -- dcio for the records, scalcs for the mechanism --
so it runs in the dcio-integration job rather than the build matrix, and skips
rather than fails where those are absent.

The fit here is the CH82 sample record that ships with HJCFIT, which is the
same fit documented in documentation/source/api/python/fitting.rst: log10 L
goes from 2286.97 at the guess to 2289.13.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

dcio = pytest.importorskip('dcio', reason='needs the [fitting] extra')
scalcs = pytest.importorskip('scalcs', reason='needs the [fitting] extra')

import HJCFIT
from HJCFIT.likelihood import cli
from HJCFIT.likelihood.fitspec import FitSpec, load_toml_bytes
from HJCFIT.likelihood.runner import (RunnerError, build_mechanism,
                                      load_records, provenance,
                                      result_as_dict, run, write_result)

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent / 'examples'

#: The specification every test below starts from, as a dictionary so that
#: individual tests can vary one field of it.
CH82 = {
    'title': 'CH82 at 100 nM',
    'data': [{'record': 'CH82', 'conc': 1e-7, 'tres': 1e-4, 'tcrit': 4e-3}],
    'mechanism': {'sample': 'CH82', 'rates': {'beta1': 15.0, 'beta2': 15000.0}},
}


def spec_for(**changes):
    return FitSpec.from_dict(dict(CH82, **changes))


def write(tmp_path, spec, name='spec.toml'):
    path = tmp_path / name
    spec.write_toml(str(path))
    return str(path)


# --------------------------------------------------------------------------
# records
# --------------------------------------------------------------------------

class TestLoadRecords:

    def test_a_sample_name_resolves(self):
        record = load_records(spec_for())[0]
        assert record.n_openings == 836
        assert len(record.groups) == 572
        assert record.n_raw == 4312

    def test_groups_are_identical_to_read_idealized_bursts(self):
        """Two callers of the same pairing of dcio functions, so the check is
        identity rather than similarity.

        read_idealized_bursts carries the argument for segmenting the record's
        *periods* rather than its resolved intervals -- imposing a dead time
        emits a fresh open interval at every change of fitted amplitude, so a
        record with sub-conductance levels does not alternate open and shut.
        The runner has to do the same thing, and this is what says it does.
        """
        mine = load_records(spec_for())[0].groups
        theirs = HJCFIT.read_idealized_bursts('CH82', tau=1e-4, tcrit=4e-3)
        assert len(mine) == len(theirs)
        for a, b in zip(mine, theirs):
            assert np.array_equal(np.asarray(a, dtype=float),
                                  np.asarray(b, dtype=float))

    def test_no_tcrit_gives_one_group_starting_and_ending_open(self):
        record = load_records(spec_for(data=[
            {'record': 'CH82', 'conc': 1e-7, 'tres': 1e-4}]))[0]
        assert len(record.groups) == 1
        assert len(record.groups[0]) % 2 == 1       # the likelihood requires it
        assert record.tcrit is None                 # equilibrium vectors
        assert record.n_openings == 836             # every opening, as in bursts

    def test_equilibrium_vectors_on_groups_still_gives_groups(self):
        """The cluster case: divided at a critical time, but started and ended
        with equilibrium vectors rather than CHS ones. tcrit does two jobs and
        only one of them is being turned off."""
        record = load_records(spec_for(data=[
            {'record': 'CH82', 'conc': 1e-7, 'tres': 1e-4, 'tcrit': 4e-3,
             'vectors': 'equilibrium'}]))[0]
        assert len(record.groups) == 572
        assert record.tcrit is None

    def test_a_record_that_is_not_there_says_which_samples_are(self):
        with pytest.raises(RunnerError, match='CH82, CO and CCO'):
            load_records(spec_for(data=[
                {'record': 'not-a-file.scn', 'conc': 1e-7, 'tres': 1e-4}]))

    def test_several_records_come_back_in_order(self):
        records = load_records(spec_for(data=[
            {'record': 'CH82', 'conc': 1e-8, 'tres': 1e-4, 'tcrit': 4e-3},
            {'record': 'CH82', 'conc': 1e-7, 'tres': 1e-4, 'tcrit': 4e-3}]))
        assert [r.conc for r in records] == [1e-8, 1e-7]


# --------------------------------------------------------------------------
# the mechanism
# --------------------------------------------------------------------------

class TestBuildMechanism:

    def test_the_guess_is_applied(self):
        mec = build_mechanism(spec_for())
        rates = {r.name: float(np.ravel(r.rateconstants)[0]) for r in mec.Rates}
        assert rates['beta1'] == 15.0 and rates['beta2'] == 15000.0

    def test_eight_free_parameters_for_ch82(self):
        assert build_mechanism(spec_for()).get_free_parameter_names() == [
            'beta1', 'beta2', 'alpha1', 'alpha2', 'k(-1)', '2k(-2)',
            '2k(+1)', 'k(+2)']

    def test_a_fixed_rate_leaves_the_free_list(self):
        mec = build_mechanism(spec_for(mechanism={
            'sample': 'CH82', 'fixed': {'k(-1)': 2000.0}}))
        assert 'k(-1)' not in mec.get_free_parameter_names()
        assert len(mec.get_free_parameter_names()) == 7

    def test_limits_reach_the_rate(self):
        mec = build_mechanism(spec_for(mechanism={
            'sample': 'CH82', 'limits': {'beta1': [0.5, 50.0]}}))
        beta1 = next(r for r in mec.Rates if r.name == 'beta1')
        assert list(beta1.limits[0]) == [0.5, 50.0]

    def test_an_unknown_rate_name_lists_the_names_there_are(self):
        """The interface problem this module exists to fix: SCALCS addresses
        rates by index into mec.Rates, and an index is what gets got wrong
        silently."""
        with pytest.raises(RunnerError) as error:
            build_mechanism(spec_for(mechanism={
                'sample': 'CH82', 'rates': {'beta9': 1.0}}))
        assert "no rate called 'beta9'" in str(error.value)
        assert 'beta1' in str(error.value) and '2k(+1)' in str(error.value)

    def test_an_unknown_sample_lists_the_samples_there_are(self):
        with pytest.raises(RunnerError) as error:
            build_mechanism(spec_for(mechanism={'sample': 'NotAMechanism'}))
        assert 'CH82' in str(error.value) and 'GlyR_flip' in str(error.value)

    def test_nfree_is_asserted_before_the_search(self):
        with pytest.raises(RunnerError, match='nfree says 7'):
            build_mechanism(spec_for(mechanism={'sample': 'CH82', 'nfree': 7}))
        build_mechanism(spec_for(mechanism={'sample': 'CH82', 'nfree': 8}))

    def test_fixing_everything_is_refused(self):
        every = {r.name: 1.0 for r in build_mechanism(spec_for()).Rates}
        with pytest.raises(RunnerError, match='nothing to fit'):
            build_mechanism(spec_for(mechanism={'sample': 'CH82',
                                                'fixed': every}))

    def test_mr_constrains_a_rate_that_is_in_the_cycle(self):
        mec = build_mechanism(spec_for(mechanism={
            'sample': 'CH82', 'mr': 'beta1'}))
        assert 'beta1' not in mec.get_free_parameter_names()
        assert len(mec.get_free_parameter_names()) == 8   # 2k*(-2) freed

    def test_mr_on_a_rate_outside_the_cycle_is_refused(self):
        """SCALCS warns on stderr and carries on, and what it does then is
        worse than nothing: the rate leaves the free-parameter list while the
        cycle's constraint stays where it was, so the rate is silently frozen
        at its initial guess and the fit has one parameter fewer than the
        specification says. On CH82, mr = "k(-1)" turned eight free parameters
        into seven with k(-1) stuck at the guess.
        """
        with pytest.raises(RunnerError) as error:
            build_mechanism(spec_for(mechanism={'sample': 'CH82',
                                                'mr': 'k(-1)'}))
        assert 'not in cycle 0' in str(error.value)
        assert 'A2R*' in str(error.value)       # the states it could have used

    def test_mr_cycle_out_of_range_is_refused(self):
        with pytest.raises(RunnerError, match='not a cycle of this mechanism'):
            build_mechanism(spec_for(mechanism={
                'sample': 'CH82', 'mr': 'beta1', 'mr_cycle': 3}))


# --------------------------------------------------------------------------
# the fit
# --------------------------------------------------------------------------

class TestRun:

    @pytest.fixture(scope='class')
    def outcome(self):
        return run(spec_for())

    def test_it_reaches_the_documented_maximum(self, outcome):
        """The fit in fitting.rst: 2286.97 at the guess, 2289.13 at the end."""
        assert outcome.result.log10_likelihood > 2289.0
        assert outcome.result.success
        assert outcome.result.nfailures == 0

    def test_the_free_values_are_rates_not_logarithms(self, outcome):
        """log_params searches the logarithms; a result reporting them would
        be a units error that looks like an implausible mechanism."""
        values = dict(zip(outcome.result.free_names, outcome.result.free_values))
        assert 1e4 < values['k(+2)'] < 1e10
        assert 1.0 < values['beta1'] < 100.0

    def test_the_mechanism_carries_the_fitted_rates(self, outcome):
        """Predicted distributions are drawn from this, so it must not still
        hold the guess."""
        rates = {r.name: float(np.ravel(r.rateconstants)[0])
                 for r in outcome.mec.Rates}
        assert rates['beta1'] != 15.0
        assert rates['beta1'] == pytest.approx(
            dict(zip(outcome.result.free_names,
                     outcome.result.free_values))['beta1'])

    def test_nothing_ends_on_a_limit_here(self, outcome):
        assert outcome.limited == []

    def test_a_rate_driven_onto_a_limit_is_reported(self):
        """A rate outside its limits is reset before the likelihood sees it,
        and nothing announces it, so the fitted value reads like an estimate
        when it is a statement that the model forbids where the likelihood
        wanted to go."""
        outcome = run(spec_for(
            mechanism={'sample': 'CH82', 'rates': {'beta1': 4.0},
                       'limits': {'beta1': [1.0, 5.0]}},
            search={'maxfev': 3000}))
        assert [name for name, _, _, _ in outcome.limited] == ['beta1']
        name, value, which, limit = outcome.limited[0]
        assert which == 'upper' and value == pytest.approx(5.0)
        assert 'not estimates' in str(outcome)

    def test_the_scipy_search_reaches_the_same_maximum(self):
        """Different paths, same maximum. Agreement on the height while
        disagreeing on where it sits is evidence about the likelihood rather
        than about either search."""
        outcome = run(spec_for(search={'method': 'scipy', 'restarts': 4}))
        assert outcome.result.log10_likelihood > 2289.0

    def test_str_says_what_was_fitted_as_well_as_what_came_out(self, outcome):
        text = str(outcome)
        assert '572 groups' in text and '836 openings' in text
        assert 'log10(L)' in text and 'beta1' in text

    def test_ln_and_log10_differ_by_ln_10(self, outcome):
        """Anything treating the likelihood as a statistical quantity needs
        natural logarithms; using log10 inflates every standard deviation by
        sqrt(ln 10) = 1.517, which looks like a bad fit rather than a units
        error."""
        d = result_as_dict(outcome)
        assert d['fit']['ln_likelihood'] == pytest.approx(
            d['fit']['log10_likelihood'] * np.log(10.0))


# --------------------------------------------------------------------------
# the written result
# --------------------------------------------------------------------------

class TestWrittenResult:

    def test_it_carries_the_specification_that_produced_it(self, tmp_path):
        outcome = run(spec_for(search={'maxfev': 200}))
        path = tmp_path / 'result.json'
        write_result(str(path), outcome)
        d = json.loads(path.read_text(encoding='utf-8'))
        assert FitSpec.from_dict(d['spec']) == outcome.spec

    def test_it_carries_the_versions_it_was_computed_with(self, tmp_path):
        outcome = run(spec_for(search={'maxfev': 200}))
        d = result_as_dict(outcome)
        assert d['provenance']['versions']['hjcfit']
        assert d['provenance']['versions']['scalcs']
        assert d['provenance']['python'] == sys.version.split()[0]

    def test_it_does_not_carry_the_records(self, tmp_path):
        """Thousands of floats that are already in the .scn file the
        specification names."""
        outcome = run(spec_for(search={'maxfev': 200}))
        path = tmp_path / 'result.json'
        write_result(str(path), outcome)
        assert path.stat().st_size < 8000
        d = json.loads(path.read_text(encoding='utf-8'))
        assert d['records'][0]['groups'] == 572
        assert d['records'][0]['intervals'] == 1100

    def test_provenance_has_no_hole_where_a_missing_package_would_be(self):
        prov = provenance()
        assert set(prov['versions']) >= {'hjcfit', 'dcio', 'scalcs', 'numpy'}


# --------------------------------------------------------------------------
# the command
# --------------------------------------------------------------------------

class TestCommand:

    def test_template_then_check_then_fit(self, tmp_path, capsys, monkeypatch):
        """The path a new user takes, end to end."""
        monkeypatch.chdir(tmp_path)
        path = str(tmp_path / 'mine.toml')

        assert cli.main(['template', '-o', path]) == 0
        assert 'wrote' in capsys.readouterr().out

        assert cli.main(['check', path]) == 0
        out = capsys.readouterr().out
        assert '572 groups' in out
        assert '8 free:' in out
        assert 'log10L at the guess: 2286.9746' in out
        assert 'nothing was fitted' in out

        assert cli.main(['fit', path, '-o', 'out.json']) == 0
        out = capsys.readouterr().out
        assert 'log10(L) = 2289.13' in out
        d = json.loads((tmp_path / 'out.json').read_text(encoding='utf-8'))
        assert d['fit']['log10_likelihood'] > 2289.0

    def test_template_refuses_to_overwrite_without_force(self, tmp_path,
                                                         capsys):
        path = str(tmp_path / 'mine.toml')
        assert cli.main(['template', '-o', path]) == 0
        assert cli.main(['template', '-o', path]) == 1
        assert 'pass --force' in capsys.readouterr().err
        assert cli.main(['template', '-o', path, '--force']) == 0

    def test_show_normalises_a_specification(self, tmp_path, capsys):
        path = write(tmp_path, spec_for())
        assert cli.main(['show', path]) == 0
        out = capsys.readouterr().out
        assert FitSpec.from_dict(load_toml_bytes(out.encode())) == spec_for()

    def test_a_search_that_does_not_converge_exits_2(self, tmp_path, capsys):
        """The number is still printed -- it is the exit status that says so,
        because that is what a script reads."""
        path = write(tmp_path, spec_for(search={'maxfev': 50}))
        assert cli.main(['fit', path, '-q']) == 2
        assert 'log10(L)' in capsys.readouterr().out

    def test_a_bad_specification_is_one_line_on_stderr(self, tmp_path, capsys):
        path = tmp_path / 'bad.toml'
        path.write_text('[[data]]\nrecord = "CH82"\nconc = 1e-7\n',
                        encoding='utf-8')
        assert cli.main(['check', str(path)]) == 1
        err = capsys.readouterr().err
        assert err.startswith('hjcfit: ') and err.count('\n') == 1
        assert 'tres is required' in err

    def test_a_missing_file_is_not_a_traceback(self, tmp_path, capsys):
        assert cli.main(['check', str(tmp_path / 'nope.toml')]) == 1
        assert capsys.readouterr().err == (
            'hjcfit: no such file: {0}\n'.format(tmp_path / 'nope.toml'))

    def test_an_impossible_rate_name_is_one_line_too(self, tmp_path, capsys):
        path = write(tmp_path, spec_for(
            mechanism={'sample': 'CH82', 'rates': {'beta9': 1.0}}))
        assert cli.main(['check', path]) == 1
        err = capsys.readouterr().err
        assert err.count('\n') == 1 and 'beta9' in err

    def test_version_lists_everything_involved(self, capsys):
        assert cli.main(['--version']) == 0
        out = capsys.readouterr().out
        for package in ('hjcfit', 'dcio', 'scalcs', 'numpy', 'scipy'):
            assert package in out

    def test_no_arguments_prints_help_rather_than_failing(self, capsys):
        assert cli.main([]) == 0
        assert 'hjcfit' in capsys.readouterr().out

    def test_the_console_script_entry_point_resolves(self):
        """pyproject declares hjcfit = HJCFIT.likelihood.cli:main. A declared
        entry point nothing imports is untested and will one day be wrong."""
        from importlib.metadata import entry_points
        scripts = {e.name: e for e in entry_points(group='console_scripts')}
        assert 'hjcfit' in scripts
        assert scripts['hjcfit'].load() is cli.main


# --------------------------------------------------------------------------
# the shipped example, and the notebook over it
# --------------------------------------------------------------------------

@pytest.mark.skipif(not EXAMPLES.is_dir(),
                    reason='examples/ is not beside the tests')
class TestExamples:

    def test_the_example_specification_is_valid_and_its_nfree_holds(self):
        """It states nfree = 8, which is an assertion about how SCALCS handles
        constraints. If that changes, this is where it shows."""
        spec = FitSpec.from_toml(str(EXAMPLES / 'CH82.toml'))
        assert spec.mechanism.nfree == 8
        assert len(build_mechanism(spec).get_free_parameter_names()) == 8

    def test_the_notebook_template_runs(self, tmp_path, monkeypatch):
        """A notebook is the de facto interface to this package, and an
        example that has stopped working is worse than none. This is the one
        thing a desktop interface could not have: it is testable.
        """
        nbformat = pytest.importorskip('nbformat')
        nbclient = pytest.importorskip('nbclient')
        pytest.importorskip('matplotlib')

        monkeypatch.setenv('MPLBACKEND', 'Agg')
        # Copied, because the notebook writes result.json beside itself and a
        # test must not leave anything in the source tree.
        import shutil
        for name in ('fit_template.ipynb', 'CH82.toml'):
            shutil.copy(EXAMPLES / name, tmp_path / name)

        nb = nbformat.read(str(tmp_path / 'fit_template.ipynb'), as_version=4)
        nbclient.NotebookClient(
            nb, timeout=900, kernel_name='python3',
            resources={'metadata': {'path': str(tmp_path)}}).execute()

        text = "\n".join(
            output.get('text', '')
            for cell in nb.cells for output in cell.get('outputs', [])
            if output.output_type == 'stream')
        assert 'log10(L) = 2289.13' in text
        # The two checks the notebook makes on its own curves. A predicted
        # curve over an observed histogram is the one plot where a wrong curve
        # looks like a finding about the fit rather than a bug.
        assert 'pdf integrates to 1.00000' in text
        assert os.path.exists(tmp_path / 'result.json')
