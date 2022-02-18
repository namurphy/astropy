# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains tests for the name resolve convenience module.
"""


import time
import urllib.request

import pytest
import numpy as np

from astropy.coordinates.name_resolve import (get_icrs_coordinates,
                                              NameResolveError,
                                              sesame_database, _parse_response,
                                              sesame_url)
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.config import paths
from astropy import units as u

from pytest_remotedata.disable_internet import no_internet

_cached_ngc3642 = {
    'simbad': '# NGC 3642    #Q22523669\x1f#=S=Simbad (via url):    1\x1f%@ 503952\x1f%I.0 NGC 3642\x1f%C.0 LIN\x1f%C.N0 15.15.01.00\x1f%J 170.5750583 +59.0742417 = 11:22:18.01 +59:04:27.2\x1f%V z 1593 0.005327 [0.000060] D 2002LEDA.........0P\x1f%D 1.673 1.657 75 (32767) (I) C 2006AJ....131.1163S\x1f%T 5 =32800000 D 2011A&A...532A..74B\x1f%#B 140\x1f\x1f\x1f#====Done (2013-Feb-12,16:37:11z)====',
    'vizier': '# NGC 3642    #Q22523677\x1f#=V=VizieR (local):    1\x1f%J 170.56 +59.08 = 11:22.2     +59:05\x1f%I.0 {NGC} 3642\x1f\x1f\x1f\x1f#====Done (2013-Feb-12,16:37:42z)====',
    'all': '# ngc3642    #Q22523722\x1f#=S=Simbad (via url):    1\x1f%@ 503952\x1f%I.0 NGC 3642\x1f%C.0 LIN\x1f%C.N0 15.15.01.00\x1f%J 170.5750583 +59.0742417 = 11:22:18.01 +59:04:27.2\x1f%V z 1593 0.005327 [0.000060] D 2002LEDA.........0P\x1f%D 1.673 1.657 75 (32767) (I) C 2006AJ....131.1163S\x1f%T 5 =32800000 D 2011A&A...532A..74B\x1f%#B 140\x1f\x1f\x1f#=V=VizieR (local):    1\x1f%J 170.56 +59.08 = 11:22.2     +59:05\x1f%I.0 {NGC} 3642\x1f\x1f\x1f#!N=NED : *** Could not access the server ***\x1f\x1f#====Done (2013-Feb-12,16:39:48z)====',
}

_cached_castor = {
    'all': '# castor    #Q22524249\x1f#=S=Simbad (via url):    1\x1f%@ 983633\x1f%I.0 NAME CASTOR\x1f%C.0 **\x1f%C.N0 12.13.00.00\x1f%J 113.649471640 +31.888282216 = 07:34:35.87 +31:53:17.8\x1f%J.E [34.72 25.95 0] A 2007A&A...474..653V\x1f%P -191.45 -145.19 [3.95 2.95 0] A 2007A&A...474..653V\x1f%X 64.12 [3.75] A 2007A&A...474..653V\x1f%S A1V+A2Vm =0.0000D200.0030.0110000000100000 C 2001AJ....122.3466M\x1f%#B 179\x1f\x1f#!V=VizieR (local): No table found for: castor\x1f\x1f#!N=NED: ****object name not recognized by NED name interpreter\x1f#!N=NED: ***Not recognized by NED: castor\x1f\x1f\x1f\x1f#====Done (2013-Feb-12,16:52:02z)====',
    'simbad': '# castor    #Q22524495\x1f#=S=Simbad (via url):    1\x1f%@ 983633\x1f%I.0 NAME CASTOR\x1f%C.0 **\x1f%C.N0 12.13.00.00\x1f%J 113.649471640 +31.888282216 = 07:34:35.87 +31:53:17.8\x1f%J.E [34.72 25.95 0] A 2007A&A...474..653V\x1f%P -191.45 -145.19 [3.95 2.95 0] A 2007A&A...474..653V\x1f%X 64.12 [3.75] A 2007A&A...474..653V\x1f%S A1V+A2Vm =0.0000D200.0030.0110000000100000 C 2001AJ....122.3466M\x1f%#B 179\x1f\x1f\x1f#====Done (2013-Feb-12,17:00:39z)====',
}


@pytest.mark.remote_data
def test_names():

    # First check that sesame is up
    if urllib.request.urlopen("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame").getcode() != 200:
        pytest.skip("SESAME appears to be down, skipping test_name_resolve.py:test_names()...")

    with pytest.raises(NameResolveError):
        get_icrs_coordinates("m87h34hhh")

    try:
        icrs = get_icrs_coordinates("NGC 3642")
    except NameResolveError:
        ra, dec = _parse_response(_cached_ngc3642["all"])
        icrs = SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree)

    icrs_true = SkyCoord(ra="11h 22m 18.014s", dec="59d 04m 27.27s")

    # use precision of only 1 decimal here and below because the result can
    # change due to Sesame server-side changes.
    np.testing.assert_almost_equal(icrs.ra.degree, icrs_true.ra.degree, 1)
    np.testing.assert_almost_equal(icrs.dec.degree, icrs_true.dec.degree, 1)

    try:
        icrs = get_icrs_coordinates("castor")
    except NameResolveError:
        ra, dec = _parse_response(_cached_castor["all"])
        icrs = SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree)

    icrs_true = SkyCoord(ra="07h 34m 35.87s", dec="+31d 53m 17.8s")
    np.testing.assert_almost_equal(icrs.ra.degree, icrs_true.ra.degree, 1)
    np.testing.assert_almost_equal(icrs.dec.degree, icrs_true.dec.degree, 1)


@pytest.mark.remote_data
def test_name_resolve_cache(tmpdir):
    from astropy.utils.data import get_cached_urls

    target_name = "castor"

    temp_cache_dir = str(tmpdir.mkdir('cache'))
    with paths.set_temp_cache(temp_cache_dir, delete=True):
        assert len(get_cached_urls()) == 0

        icrs1 = get_icrs_coordinates(target_name, cache=True)

        urls = get_cached_urls()
        assert len(urls) == 1
        assert 'http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/' in urls[0]

        # Try reloading coordinates, now should just reload cached data:
        with no_internet():
            icrs2 = get_icrs_coordinates(target_name, cache=True)

        assert len(get_cached_urls()) == 1

        assert u.allclose(icrs1.ra, icrs2.ra)
        assert u.allclose(icrs1.dec, icrs2.dec)


def test_names_parse():
    # a few test cases for parsing embedded coordinates from object name
    test_names = ['CRTS SSS100805 J194428-420209',
                  'MASTER OT J061451.7-272535.5',
                  '2MASS J06495091-0737408',
                  '1RXS J042555.8-194534',
                  'SDSS J132411.57+032050.5',
                  'DENIS-P J203137.5-000511',
                  '2QZ J142438.9-022739',
                  'CXOU J141312.3-652013']
    for name in test_names:
        sc = get_icrs_coordinates(name, parse=True)


@pytest.mark.remote_data
@pytest.mark.parametrize(("name", "db_dict"), [('NGC 3642', _cached_ngc3642),
                                               ('castor', _cached_castor)])
def test_database_specify(name, db_dict):
    # First check that at least some sesame mirror is up
    for url in sesame_url.get():
        if urllib.request.urlopen(url).getcode() == 200:
            break
    else:
        pytest.skip("All SESAME mirrors appear to be down, skipping "
                    "test_name_resolve.py:test_database_specify()...")

    for db in db_dict.keys():
        with sesame_database.set(db):
            icrs = SkyCoord.from_name(name)

        time.sleep(1)
