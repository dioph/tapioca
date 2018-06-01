import pytest
from ..utils import get_lc_kepler

def test_get_lc_k2():
    assert len(get_lc_kepler(206078331)) == 1
