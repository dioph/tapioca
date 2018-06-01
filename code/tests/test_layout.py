import pytest
from ..tapioca import make_url

def test_make_url():
    assert make_url() == "https://exoplanetarchive.ipac.caltech.edu/" \
                         "cgi-bin/nstedAPI/nph-nstedAPI?table=q1_q17_dr25_stellar"
    with pytest.raises(AssertionError) as exc:
        make_url(select='eggs')
    assert "property" in str(exc)