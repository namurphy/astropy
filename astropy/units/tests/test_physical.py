# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for the physical_type support in the units package
"""

import pytest

from astropy import units as u
from astropy.units import physical
from astropy.constants import hbar
from astropy.tests.helper import raises


def test_simple():
    assert u.m.physical_type == 'length'


def test_power():
    assert (u.cm ** 3).physical_type == 'volume'


def test_speed():
    assert (u.km / u.h).physical_type == 'speed'


def test_unknown():
    assert (u.m * u.s).physical_type == 'unknown'


def test_dimensionless():
    assert (u.m / u.m).physical_type == 'dimensionless'


def test_angular_momentum():
    assert hbar.unit.physical_type == 'angular momentum'


def test_flam():
    flam = u.erg / (u.cm**2 * u.s * u.AA)
    assert flam.physical_type == 'spectral flux density wav'


def test_photlam():
    photlam = u.photon / (u.cm ** 2 * u.s * u.AA)
    assert photlam.physical_type == 'photon flux density wav'


def test_photnu():
    photnu = u.photon / (u.cm ** 2 * u.s * u.Hz)
    assert photnu.physical_type == 'photon flux density'


@raises(ValueError)
def test_redundant_physical_type():
    physical.def_physical_type(u.m, 'utter craziness')


def test_data_quantity():
    assert u.byte.physical_type == 'data quantity'
    assert u.bit.physical_type == 'data quantity'


@pytest.fixture
def length():
    return physical._PhysicalType(u.m)


@pytest.fixture
def time():
    return physical._PhysicalType(u.s)


@pytest.fixture
def velocity():
    return physical._PhysicalType(u.m / u.s)


@pytest.fixture
def area():
    return physical._PhysicalType(u.m ** 2)


def test_physical_type_equality():
    length1 = physical._PhysicalType(u.m)
    length2 = physical._PhysicalType(u.m)
    assert (length1 == length2) is True
    assert (length1 != length2) is False


def test_physical_type_inequality(length, time):
    assert (length != time) is True
    assert (length == time) is False



def test_physical_type_as_set(length):
    assert length.as_set == {"length"}


def test_physical_type_contains(length):
    assert "length" in length


def test_physical_type_multiplication(length, time, velocity):
    assert velocity * time == length


def test_physical_type_division(length, time, velocity):
    assert length / time == velocity


def test_physical_type_power(length, area):
    assert length ** 2 == area
