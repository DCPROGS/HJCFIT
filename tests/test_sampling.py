"""Sampling from a specification: the [likelihood] and [mcmc] tables, split-
R-hat, the runner's sample() and the ``hjcfit sample`` command.

The specification tests need nothing beyond HJCFIT. The runner and command
tests need the [fitting] extra and skip without it, as in
test_runner_cli.py. They sample the CH82 sample record, which is quick, and
whose one concentration leaves one direction barely constrained. That makes
it a good test of how the run reports a poor Hessian and chains that have
not yet agreed.
"""

import json

import numpy as np
import pytest

from HJCFIT.likelihood.fitspec import (
    TEMPLATE, FitSpec, LikelihoodSpec, MCMCSpec, SpecError, load_toml_bytes)
from HJCFIT.likelihood.mcmc import Chain, potential_scale_reduction

CH82 = {
    'title': 'CH82 at 100 nM',
    'data': [{'record': 'CH82', 'conc': 1e-7, 'tres': 1e-4, 'tcrit': 4e-3}],
    'mechanism': {'sample': 'CH82', 'rates': {'beta1': 15.0,
                                              'beta2': 15000.0}},
}


def spec_for(**changes):
    return FitSpec.from_dict(dict(CH82, **changes))


# --------------------------------------------------------------------------
# [likelihood]
# --------------------------------------------------------------------------

class TestLikelihoodSpec:

    def test_its_names_are_the_fitters(self):
        from HJCFIT.likelihood.fitspec import SOLVER_KEYS
        from HJCFIT.likelihood.fitting import SOLVER_OPTIONS
        assert SOLVER_KEYS == SOLVER_OPTIONS

    def test_only_what_is_given_reaches_the_solver(self):
        spec = LikelihoodSpec.from_dict({'nmax': 2, 'xtol': 1e-12})
        assert spec.solver() == {'nmax': 2, 'xtol': 1e-12}
        assert LikelihoodSpec().solver() == {}

    @pytest.mark.parametrize("table, message", [
        ({'nmax': 0}, "nmax"),
        ({'nmax': 2.5}, "integer"),
        ({'nmax': True}, "integer"),
        ({'xtol': -1e-12}, "positive"),
        ({'lower_bound': 0.0, 'upper_bound': -1.0}, "below"),
        ({'tol': 1e-9}, "tol"),
    ])
    def test_nonsense_is_refused(self, table, message):
        with pytest.raises(SpecError, match=message):
            LikelihoodSpec.from_dict(table)


# --------------------------------------------------------------------------
# [mcmc]
# --------------------------------------------------------------------------

class TestMCMCSpec:

    def test_defaults(self):
        m = MCMCSpec.from_dict({})
        assert (m.start, m.sampler, m.n, m.burnin, m.chains) == (
            'fit', 'adaptive', 20000, 5000, 1)
        assert m.as_dict() == {}

    def test_only_departures_are_written_back(self):
        m = MCMCSpec.from_dict({'n': 5000, 'burnin': 1000, 'chains': 4})
        assert m.as_dict() == {'n': 5000, 'burnin': 1000, 'chains': 4}

    @pytest.mark.parametrize("table, message", [
        ({'start': 'mode'}, "start"),
        ({'sampler': 'hmc'}, "sampler"),
        ({'prior': 'flat'}, "prior"),
        ({'covariance': 'fisher'}, "covariance"),
        ({'n': 1000, 'burnin': 1000}, "keeps nothing"),
        ({'chains': 0}, "chains"),
        ({'n': 1e4}, "integer"),
        ({'log_space': 1}, "true or false"),
        ({'iterations': 100}, "iterations"),
    ])
    def test_nonsense_is_refused(self, table, message):
        with pytest.raises(SpecError, match=message):
            MCMCSpec.from_dict(table)

    def test_a_whole_specification_round_trips(self):
        spec = spec_for(likelihood={'nmax': 2, 'xtol': 1e-12},
                        mcmc={'n': 5000, 'burnin': 1000, 'chains': 3,
                              'start': 'guess', 'log_space': True})
        back = FitSpec.from_dict(load_toml_bytes(spec.to_toml().encode()))
        assert back == spec

    def test_without_the_tables_nothing_is_written(self):
        text = spec_for().to_toml()
        assert '[likelihood]' not in text and '[mcmc]' not in text

    def test_str_says_how_it_would_sample(self):
        text = str(spec_for(mcmc={'chains': 4},
                            likelihood={'nmax': 2}))
        assert 'sampling: 4 x 20000 adaptive' in text
        assert 'likelihood: nmax = 2' in text

    def test_the_templates_example_sections_are_valid(self):
        """The commented [likelihood] and [mcmc] blocks, uncommented."""
        lines, inside = [], False
        for line in TEMPLATE.splitlines():
            if line in ('# [likelihood]', '# [mcmc]'):
                inside = True
            elif not line.strip():
                inside = False
            if inside and line.startswith('# ') and (
                    line.startswith('# [') or '=' in line):
                line = line[2:]
            lines.append(line)
        text = "\n".join(lines)
        assert '\n[mcmc]' in text and '\n[likelihood]' in text
        spec = FitSpec.from_dict(load_toml_bytes(text.encode()))
        assert spec.mcmc.chains == 4
        assert spec.likelihood.solver() == {'nmax': 2, 'xtol': 1e-12,
                                            'rtol': 1e-12}


# --------------------------------------------------------------------------
# split-R-hat
# --------------------------------------------------------------------------

class TestPotentialScaleReduction:

    def test_chains_of_one_distribution_give_one(self):
        rng = np.random.default_rng(0)
        chains = [rng.standard_normal(4000) for _ in range(4)]
        assert potential_scale_reduction(chains) == pytest.approx(1.0,
                                                                  abs=0.01)

    def test_chains_that_disagree_give_more(self):
        rng = np.random.default_rng(1)
        chains = [rng.standard_normal(2000) + shift for shift in (0, 0, 3)]
        assert potential_scale_reduction(chains) > 1.5

    def test_a_single_chain_that_drifts_is_caught(self):
        drift = np.linspace(0, 10, 4000) + np.random.default_rng(2).normal(
            size=4000)
        assert potential_scale_reduction([drift]) > 1.5

    def test_too_short_is_an_error(self):
        with pytest.raises(ValueError, match="four"):
            potential_scale_reduction([np.ones(3)])


# --------------------------------------------------------------------------
# The runner and the command: need the [fitting] extra
# --------------------------------------------------------------------------

@pytest.fixture(scope='module')
def fitting_extra():
    pytest.importorskip('dcio', reason='needs the [fitting] extra')
    pytest.importorskip('scalcs', reason='needs the [fitting] extra')


@pytest.fixture(scope='module')
def outcome(fitting_extra):
    from HJCFIT.likelihood.runner import sample
    spec = spec_for(mcmc={'n': 400, 'burnin': 100, 'chains': 2})
    return sample(spec, processes=1)


class TestSample:

    def test_chains_and_summary(self, outcome):
        assert len(outcome.chains) == 2
        assert all(c.n == 400 and c.burnin == 100 for c in outcome.chains)
        names = [row['name'] for row in outcome.summary]
        assert names == list(outcome.names) and len(names) == 8
        for row in outcome.summary:
            assert row['q2.5'] <= row['median'] <= row['q97.5']
            assert row['sd'] > 0 and row['ess'] > 0

    def test_the_first_chain_starts_at_the_fit(self, outcome):
        assert outcome.start_from == 'fit'
        np.testing.assert_allclose(outcome.chains[0].samples[0],
                                   outcome.fit.result.free_values,
                                   rtol=0.2)
        np.testing.assert_array_equal(outcome.start,
                                      outcome.fit.result.free_values)

    def test_a_poor_hessian_is_reported_not_hidden(self, outcome):
        """CH82 at one concentration has a direction the record barely
        constrains; the run says so."""
        assert outcome.approximation is not None or outcome.notes
        assert any('Hessian' in note for note in outcome.notes)

    def test_str_is_a_table_with_rhat(self, outcome):
        text = str(outcome)
        assert 'R-hat' in text and 'acceptance after burn-in' in text

    def test_written_and_read_back(self, outcome, tmp_path):
        from HJCFIT.likelihood.runner import write_samples
        paths = write_samples(tmp_path / 'post', outcome)
        assert [p.rsplit('/', 1)[-1].rsplit('\\', 1)[-1] for p in paths] == [
            'post.json', 'post_chain0.npz', 'post_chain1.npz']
        with open(paths[0]) as handle:
            data = json.load(handle)
        assert data['free_names'] == list(outcome.names)
        assert data['spec']['mcmc'] == {'n': 400, 'burnin': 100, 'chains': 2}
        assert data['provenance']['versions']['hjcfit']
        back = Chain.load(paths[1])
        np.testing.assert_array_equal(back.samples, outcome.chains[0].samples)

    def test_the_paper_scheme_starts_from_a_pilot(self, fitting_extra):
        from HJCFIT.likelihood.runner import sample
        spec = spec_for(mcmc={'start': 'guess', 'pilot_n': 60, 'n': 200,
                              'burnin': 50, 'covariance': 'identity'})
        out = sample(spec, processes=1)
        assert out.start_from == 'pilot' and out.fit is None
        assert out.pilot.n == 60
        np.testing.assert_array_equal(out.start, out.pilot.mode()[0])

    def test_the_mwg_sampler_can_run_the_chains(self, fitting_extra):
        from HJCFIT.likelihood.runner import sample
        spec = spec_for(mcmc={'sampler': 'mwg', 'n': 60, 'burnin': 20})
        out = sample(spec, processes=1)
        assert out.chains[0].accepted.shape == (60, 8)

    def test_chains_in_parallel_processes(self, fitting_extra):
        """Spawned processes, each rebuilding the posterior from the
        specification; seeded, so the result does not depend on where a
        chain ran."""
        from HJCFIT.likelihood.runner import sample
        spec = spec_for(mcmc={'n': 150, 'burnin': 50, 'chains': 2})
        here = sample(spec, processes=1)
        there = sample(spec, processes=2)
        for a, b in zip(here.chains, there.chains):
            np.testing.assert_allclose(a.samples, b.samples, rtol=1e-9)

    def test_the_likelihood_table_reaches_the_posterior(self, fitting_extra):
        from HJCFIT.likelihood.runner import _posterior
        _, _, post = _posterior(spec_for(likelihood={'nmax': 2,
                                                     'xtol': 1e-12}))
        lik = post.fitter.likelihoods[0]
        assert (lik.nmax, lik.xtol) == (2, 1e-12)

    def test_and_the_fit(self, fitting_extra, monkeypatch):
        from HJCFIT.likelihood import runner
        seen = {}
        original = runner.HJCFitter

        class Spy(original):
            def __init__(self, *args, **kwargs):
                seen.update(kwargs)
                super().__init__(*args, **kwargs)

        monkeypatch.setattr(runner, 'HJCFitter', Spy)
        runner.run(spec_for(likelihood={'nmax': 2},
                            search={'maxfev': 50}))
        assert seen['solver'] == {'nmax': 2}


class TestSampleCommand:

    def test_sample_writes_its_files(self, fitting_extra, tmp_path, capsys):
        from HJCFIT.likelihood import cli
        path = tmp_path / 'ch82.toml'
        spec_for(mcmc={'n': 200, 'burnin': 50}).write_toml(str(path))
        status = cli.main(['sample', str(path), '-o',
                           str(tmp_path / 'post'), '--processes', '1', '-q'])
        out = capsys.readouterr().out
        assert status in (0, 2)          # 2: chains this short may disagree
        assert 'R-hat' in out
        assert (tmp_path / 'post.json').exists()
        assert (tmp_path / 'post_chain0.npz').exists()

    def test_chains_on_the_command_line_override_the_file(
            self, fitting_extra, tmp_path, capsys):
        from HJCFIT.likelihood import cli
        path = tmp_path / 'ch82.toml'
        spec_for(mcmc={'n': 120, 'burnin': 20}).write_toml(str(path))
        cli.main(['sample', str(path), '--chains', '2', '--processes', '1',
                  '-o', str(tmp_path / 'p'), '-q'])
        assert (tmp_path / 'p_chain1.npz').exists()

    def test_zero_chains_is_one_line_on_stderr(self, fitting_extra,
                                              tmp_path, capsys):
        from HJCFIT.likelihood import cli
        path = tmp_path / 'ch82.toml'
        spec_for().write_toml(str(path))
        assert cli.main(['sample', str(path), '--chains', '0']) == 1
        assert capsys.readouterr().err.startswith('hjcfit: --chains')


def test_a_plain_fit_says_nothing_about_sampling():
    """Without an [mcmc] section, `hjcfit fit` and `check` should not talk
    about sampling."""
    assert 'sampling' not in str(spec_for())
