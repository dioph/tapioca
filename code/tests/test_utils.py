import pytest
from ..utils import get_lc_kepler

def test_get_lc_kepler():
    assert len(get_lc_kepler(9654627)) == 15
