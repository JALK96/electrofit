import importlib.util as u

def test_imports():
    assert u.find_spec("electrofit")
    assert u.find_spec("electrofit.cli.app")
    assert u.find_spec("electrofit.workflows")