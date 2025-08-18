import logging
from types import SimpleNamespace
from electrofit.infra.step_logging import _extract, log_relevant_config


def test__extract_dict_and_object_paths():
    obj = SimpleNamespace(a=SimpleNamespace(b=2), d={'x': {'y': 5}})
    assert _extract(obj, 'a.b') == 2
    assert _extract(obj, 'd.x.y') == 5
    assert _extract(obj, 'missing.path') is None
    assert _extract(None, 'a') is None


def test_log_relevant_config_order_and_values(caplog):
    caplog.set_level(logging.INFO)
    obj = SimpleNamespace(m=1, n=2)
    summary = log_relevant_config('stepX', obj, ['m','n'])
    assert summary == {'m':1,'n':2}
    # Ensure order preserved (m before n) in log text
    log_line = next((r.message for r in caplog.records if '[stepX][cfg]' in r.message), '')
    first_idx = log_line.find('m=')
    second_idx = log_line.find('n=')
    assert first_idx != -1 and second_idx != -1 and first_idx < second_idx
