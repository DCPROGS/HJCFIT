"""The fit specification: parsing, validation and round-tripping.

Pure Python. Nothing here imports dcio, scalcs or even the likelihood, because
fitspec does not either -- a specification describes a fit and executes
nothing, and this file is the check that it stays that way.

test_runner_cli.py is the other half: carrying one out.
"""

import subprocess
import sys

import pytest

from HJCFIT.likelihood.fitspec import (
    TEMPLATE, DataSpec, FitSpec, MechanismSpec, SearchSpec, SpecError,
    load_toml_bytes)


MINIMAL = {
    'data': [{'record': 'CH82', 'conc': 1e-7, 'tres': 1e-4, 'tcrit': 4e-3}],
    'mechanism': {'sample': 'CH82'},
}


# --------------------------------------------------------------------------
# the module's one structural promise
# --------------------------------------------------------------------------

def test_fitspec_imports_nothing_but_the_standard_library():
    """A specification must be readable without the [fitting] extra.

    The whole point of separating the description from the run is that
    anything can read one -- a batch runner, a user interface, a script that
    only wants to list what a directory of specifications contains. A stray
    import of scalcs here would make reading a file cost a mechanism library.
    """
    code = (
        "import sys;"
        "import HJCFIT.likelihood.fitspec;"
        "bad = [m for m in ('scalcs', 'dcio', 'matplotlib', 'scipy')"
        "       if m in sys.modules];"
        "print(bad)"
    )
    out = subprocess.run([sys.executable, '-c', code], check=True,
                         capture_output=True, text=True)
    assert out.stdout.strip() == '[]', out.stdout


# --------------------------------------------------------------------------
# DataSpec
# --------------------------------------------------------------------------

class TestDataSpec:

    def test_vectors_default_follows_tcrit(self):
        """CHS vectors are defined between groups, so no groups means no CHS."""
        assert DataSpec.from_dict(
            {'record': 'r', 'conc': 1e-7, 'tres': 1e-4,
             'tcrit': 4e-3}).vectors == 'chs'
        assert DataSpec.from_dict(
            {'record': 'r', 'conc': 1e-7, 'tres': 1e-4}).vectors == 'equilibrium'

    def test_chs_without_groups_is_refused(self):
        with pytest.raises(SpecError, match='CHS vectors are defined between'):
            DataSpec.from_dict({'record': 'r', 'conc': 1e-7, 'tres': 1e-4,
                                'vectors': 'chs'})

    def test_tcrit_below_tres_is_refused(self):
        with pytest.raises(SpecError, match='divides nothing'):
            DataSpec.from_dict({'record': 'r', 'conc': 1e-7, 'tres': 1e-4,
                                'tcrit': 5e-5})

    def test_a_negative_tcrit_is_not_a_flag_here(self):
        """Elsewhere in the stack the sign selects the vectors. Not here."""
        with pytest.raises(SpecError, match='sign is not a flag'):
            DataSpec.from_dict({'record': 'r', 'conc': 1e-7, 'tres': 1e-4,
                                'tcrit': -4e-3, 'vectors': 'equilibrium'})

    def test_zero_concentration_is_allowed(self):
        """A mechanism with no agonist is an experiment, not a mistake."""
        assert DataSpec.from_dict(
            {'record': 'r', 'conc': 0.0, 'tres': 1e-4}).conc == 0.0

    def test_a_missing_required_field_names_itself(self):
        with pytest.raises(SpecError, match='tres is required'):
            DataSpec.from_dict({'record': 'r', 'conc': 1e-7})

    @pytest.mark.parametrize('bad', ['vector', 'tcritical', 'resolution'])
    def test_an_unknown_key_is_refused_with_the_known_ones(self, bad):
        """A misspelled key is worse than a missing one: the fit would run,
        silently, without the setting the file says it has."""
        d = {'record': 'r', 'conc': 1e-7, 'tres': 1e-4, bad: 1.0}
        with pytest.raises(SpecError) as error:
            DataSpec.from_dict(d)
        assert bad in str(error.value) and 'known keys are' in str(error.value)

    def test_a_string_where_a_number_belongs_is_refused(self):
        """The failure PyYAML would hand us silently. See the module docstring:
        yaml.safe_load reads 1e8 as the string '1e8'."""
        with pytest.raises(SpecError, match='expected a number'):
            DataSpec.from_dict({'record': 'r', 'conc': '1e-7', 'tres': 1e-4})

    def test_a_boolean_is_not_a_number(self):
        with pytest.raises(SpecError, match='expected a number'):
            DataSpec.from_dict({'record': 'r', 'conc': True, 'tres': 1e-4})

    def test_str_is_readable(self):
        text = str(DataSpec.from_dict(MINIMAL['data'][0]))
        assert '100 nM' in text and '100 us' in text and '4 ms' in text


# --------------------------------------------------------------------------
# MechanismSpec
# --------------------------------------------------------------------------

class TestMechanismSpec:

    def test_exactly_one_source(self):
        with pytest.raises(SpecError, match='exactly one of sample and'):
            MechanismSpec.from_dict({})
        with pytest.raises(SpecError, match='exactly one of sample and'):
            MechanismSpec.from_dict({'sample': 'CH82', 'mec_file': 'a.mec'})

    def test_mec_number_without_a_file_means_nothing(self):
        with pytest.raises(SpecError, match='means nothing without mec_file'):
            MechanismSpec.from_dict({'sample': 'CH82', 'mec_number': 2})

    def test_a_rate_cannot_be_both_a_guess_and_fixed(self):
        with pytest.raises(SpecError, match='both rates and fixed'):
            MechanismSpec.from_dict({'sample': 'CH82',
                                     'rates': {'beta1': 15.0},
                                     'fixed': {'beta1': 15.0}})

    def test_a_rate_cannot_be_both_fixed_and_from_the_ec50(self):
        with pytest.raises(SpecError, match='cannot be both'):
            MechanismSpec.from_dict({
                'sample': 'CH82', 'fixed': {'k(+2)': 1e8},
                'ec50': {'rate': 'k(+2)', 'value': 3.3e-6}})

    def test_limits_must_be_a_pair_in_order(self):
        with pytest.raises(SpecError, match=r'expected \[lower, upper\]'):
            MechanismSpec.from_dict({'sample': 'CH82',
                                     'limits': {'beta1': [1.0]}})
        with pytest.raises(SpecError, match='must be below'):
            MechanismSpec.from_dict({'sample': 'CH82',
                                     'limits': {'beta1': [10.0, 1.0]}})

    def test_ec50_needs_both_fields(self):
        with pytest.raises(SpecError, match='value is required'):
            MechanismSpec.from_dict({'sample': 'CH82',
                                     'ec50': {'rate': 'k(+2)'}})

    def test_rate_names_survive_a_round_trip(self):
        """Rate constants are called things like 2k(+1), which TOML will not
        take as a bare key. The writer has to quote them."""
        names = {'k(-1)': 1e3, '2k(+1)': 1e7, "k*(+2)": 5e8, 'beta1': 15.0}
        spec = FitSpec.from_dict(dict(MINIMAL, mechanism={
            'sample': 'CH82', 'rates': names}))
        back = FitSpec.from_dict(load_toml_bytes(spec.to_toml().encode()))
        assert back.mechanism.rates == names


# --------------------------------------------------------------------------
# SearchSpec
# --------------------------------------------------------------------------

class TestSearchSpec:

    def test_defaults_are_hjcfits_own(self):
        s = SearchSpec.from_dict({})
        assert s.method == 'simplex' and s.log_params is True

    def test_an_unknown_method_is_refused(self):
        with pytest.raises(SpecError, match='must be one of simplex or scipy'):
            SearchSpec.from_dict({'method': 'nelder-mead'})

    def test_restarts_belong_to_the_scipy_search_alone(self):
        """simplex_hjc has its own restart rule, capped by nresmax. Accepting
        restarts here would silently do nothing."""
        with pytest.raises(SpecError, match='only the scipy search'):
            SearchSpec.from_dict({'method': 'simplex', 'restarts': 8})
        assert SearchSpec.from_dict(
            {'method': 'scipy', 'restarts': 8}).restarts == 8

    def test_a_zero_budget_is_refused(self):
        with pytest.raises(SpecError, match='at least 1'):
            SearchSpec.from_dict({'maxfev': 0})


# --------------------------------------------------------------------------
# FitSpec
# --------------------------------------------------------------------------

class TestFitSpec:

    def test_a_fit_needs_a_record(self):
        with pytest.raises(SpecError, match='at least one'):
            FitSpec.from_dict({'mechanism': {'sample': 'CH82'}})

    def test_a_single_data_table_is_taken_as_one_record(self):
        """[data] rather than [[data]] is the easiest TOML mistake to make and
        means exactly one thing, so it is read rather than refused."""
        spec = FitSpec.from_dict(dict(MINIMAL, data=MINIMAL['data'][0]))
        assert len(spec.data) == 1

    def test_several_records_keep_their_order(self):
        spec = FitSpec.from_dict({
            'data': [{'record': 'a', 'conc': 1e-8, 'tres': 1e-5},
                     {'record': 'b', 'conc': 1e-7, 'tres': 1e-5},
                     {'record': 'c', 'conc': 1e-6, 'tres': 1e-5}],
            'mechanism': {'sample': 'CH82'}})
        assert [d.record for d in spec.data] == ['a', 'b', 'c']

    def test_the_index_of_a_bad_record_is_in_the_message(self):
        with pytest.raises(SpecError, match=r'data\[1\]'):
            FitSpec.from_dict({
                'data': [{'record': 'a', 'conc': 1e-8, 'tres': 1e-5},
                         {'record': 'b', 'conc': 1e-8, 'tres': -1.0}],
                'mechanism': {'sample': 'CH82'}})

    def test_round_trip_through_toml(self):
        spec = FitSpec.from_dict(dict(
            MINIMAL, title='a fit',
            mechanism={'sample': 'CH82', 'rates': {'beta1': 15.0},
                       'fixed': {'k(-1)': 2000.0},
                       'limits': {'beta1': [0.1, 1e5]},
                       'mr': 'k(-1)', 'nfree': 7},
            search={'method': 'scipy', 'restarts': 4, 'log_params': False}))
        assert FitSpec.from_dict(
            load_toml_bytes(spec.to_toml().encode())) == spec

    def test_to_toml_writes_scalars_before_sub_tables(self):
        """TOML reads every key after a [header] as belonging to it, so a
        scalar written below a sub-table would silently join it. This is a
        correctness property of the writer, not a matter of taste."""
        text = FitSpec.from_dict(dict(MINIMAL, mechanism={
            'sample': 'CH82', 'nfree': 8,
            'rates': {'beta1': 15.0}})).to_toml()
        lines = [l.strip() for l in text.splitlines()]
        assert lines.index('sample = "CH82"') < lines.index('[mechanism.rates]')
        assert lines.index('nfree = 8') < lines.index('[mechanism.rates]')

    def test_from_toml_names_the_file_in_its_error(self, tmp_path):
        path = tmp_path / 'broken.toml'
        path.write_text('[[data]]\nrecord = "r"\nconc = 1e-7\n',
                        encoding='utf-8')
        with pytest.raises(SpecError) as error:
            FitSpec.from_toml(str(path))
        assert 'broken.toml' in str(error.value)
        assert 'tres is required' in str(error.value)

    def test_unparseable_toml_says_so(self, tmp_path):
        path = tmp_path / 'nonsense.toml'
        path.write_text('this is not toml = = =\n', encoding='utf-8')
        with pytest.raises(SpecError, match='not valid TOML'):
            FitSpec.from_toml(str(path))

    def test_str_mentions_every_constraint(self):
        text = str(FitSpec.from_dict(dict(MINIMAL, mechanism={
            'sample': 'CH82', 'fixed': {'k(-1)': 2000.0}, 'mr': 'beta1',
            'limits': {'beta2': [1.0, 1e5]},
            'ec50': {'rate': 'k(+2)', 'value': 3.3e-6}})))
        assert 'k(-1) = 2000' in text
        assert 'microscopic reversibility' in text
        assert 'EC50 of 3.3 uM' in text
        assert 'limits on beta2' in text


# --------------------------------------------------------------------------
# the template, and the numbers TOML gets right
# --------------------------------------------------------------------------

class TestTemplate:

    def test_the_template_is_a_valid_specification(self):
        """It is what `hjcfit template` hands a new user, so it has to parse
        and to describe a fit that runs."""
        spec = FitSpec.from_dict(load_toml_bytes(TEMPLATE.encode()))
        assert spec.data[0].record == 'CH82'
        assert spec.mechanism.sample == 'CH82'
        assert spec.search.method == 'simplex'

    def test_the_template_round_trips(self):
        spec = FitSpec.from_dict(load_toml_bytes(TEMPLATE.encode()))
        assert FitSpec.from_dict(
            load_toml_bytes(spec.to_toml().encode())) == spec

    @pytest.mark.parametrize('text', ['1e8', '1.0e8', '1e+8', '1.0e+8',
                                      '100e-9', '1E8'])
    def test_every_way_of_writing_a_rate_constant_is_a_number(self, text):
        """The reason this format is TOML. PyYAML implements YAML 1.1, whose
        resolver wants a signed exponent: yaml.safe_load reads three of these
        six as strings, and a rate constant arriving as a string is the worst
        failure a specification could have.
        """
        value = load_toml_bytes('x = {0}\n'.format(text).encode())['x']
        assert isinstance(value, float)
