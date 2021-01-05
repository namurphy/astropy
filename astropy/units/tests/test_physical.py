# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Unit tests for the physical_type support in the units package
"""

import pytest

from astropy import units as u
from astropy.units import physical
from astropy.constants import hbar


unit_physical_type_pairs = [
    (u.m, "length"),
    (u.cm ** 3, "volume"),
    (u.km / u.h, "speed"),
    (u.barn * u.Mpc, "volume"),
    (u.m * u.s ** 8, "unknown"),
    (u.m / u.m, "dimensionless"),
    (hbar.unit, "angular momentum"),
    (hbar.unit, "angular momentum"),
    (u.erg / (u.cm ** 2 * u.s * u.AA), "spectral flux density wav"),  # flam
    (u.photon / (u.cm ** 2 * u.s * u.AA), "photon flux density wav"),  # photlam
    (u.photon / (u.cm ** 2 * u.s * u.Hz), "photon flux density"),  # photnu
    (u.byte, "data quantity"),
    (u.bit, "data quantity"),
    (u.imperial.mi / u.week, "speed"),
    (u.erg / u.s, "power"),
    (u.C / u.s, "electrical current"),
    (u.C / u.s / u.cm ** 2, "electrical current density"),
    (u.T * u.m ** 2, "magnetic flux"),
    (u.N * u.m, "energy"),
    (u.rad / u.ms, "angular speed"),
    (u.Unit(1), "dimensionless"),
    (u.m ** 2, "area"),
    (u.s, "time"),
    (u.rad, "angle"),
    (u.sr, "solid angle"),
    (u.m / u.s ** 2, "acceleration"),
    (u.Hz, "frequency"),
    (u.g, "mass"),
    (u.mol, "amount of substance"),
    (u.K, "temperature"),
    (u.deg_C, "temperature"),
    (u.imperial.deg_F, "temperature"),
    (u.N, "force"),
    (u.J, "energy"),
    (u.Pa, "pressure"),
    (u.W, "power"),
    (u.kg / u.m ** 3, "mass density"),
    (u.m ** 3 / u.kg, "specific volume"),
    (u.mol / u.m ** 3, "molar concentration"),
    (u.kg * u.m / u.s, "momentum/impulse"),  # needs to be split up
    (u.kg * u.m ** 2 / u.s, "angular momentum"),
    (u.rad / u.s, "angular speed"),
    (u.rad / u.s ** 2, "angular acceleration"),
    (u.g / (u.m * u.s), "dynamic viscosity"),
    (u.m ** 2 / u.s, "kinematic viscosity"),
    (u.m ** -1, "wavenumber"),
    (u.A, "electrical current"),
    (u.C, "electrical charge"),
    (u.V, "electrical potential"),
    (u.Ohm, "electrical resistance"),
    (u.S, "electrical conductance"),
    (u.F, "electrical capacitance"),
    (u.C * u.m, "electrical dipole moment"),
    (u.A / u.m ** 2, "electrical current density"),
    (u.V / u.m, "electrical field strength"),
    (u.C / u.m ** 2, "electrical flux density"),
    (u.C / u.m ** 3, "electrical charge density"),
    (u.F / u.m, "permittivity"),
    (u.Wb, "magnetic flux"),
    (u.T, "magnetic flux density"),
    (u.A / u.m, "magnetic field strength"),
    (u.H / u.m, "electromagnetic field strength"),
    (u.H, "inductance"),
    (u.cd, "luminous intensity"),
    (u.lm, "luminous flux"),
    (u.lx, "luminous emittance/illuminance"),  # needs to be split up
    (u.W / u.sr, "radiant intensity"),
    (u.cd / u.m ** 2, "luminance"),
    (u.astrophys.Jy, "spectral flux density"),
    (u.astrophys.R, "photon flux"),
    (u.misc.bit, "data quantity"),
    (u.misc.bit / u.s, "bandwidth"),
    (u.cgs.Franklin, "electrical charge (ESU)"),
    (u.cgs.statampere, "electrical current (ESU)"),
    (u.cgs.Biot, "electrical current (EMU)"),
    (u.cgs.abcoulomb, "electrical charge (EMU)"),
]


@pytest.mark.parametrize("unit, physical_type", unit_physical_type_pairs)
def test_physical_types(unit, physical_type):
    """
    Test that the `physical_type` attribute of `u.Unit` objects provides
    the expected physical type for various units.

    Many of these tests are used to test backwards compatibility.
    """
    if unit.physical_type != physical_type:
        pytest.fail(
            f"{repr(unit)}.physical_type was expected to return "
            f"{repr(physical_type)}, but instead returned "
            f"{unit.physical_type}."
        )


@pytest.fixture
def length():
    return u.m.physical_type


@pytest.fixture
def time():
    return u.s.physical_type


@pytest.fixture
def speed():
    return (u.m / u.s).physical_type


@pytest.fixture
def area():
    return (u.m ** 2).physical_type


@pytest.fixture
def pressure():
    return u.Pa.physical_type


@pytest.fixture
def dimensionless():
    return u.dimensionless_unscaled.physical_type


def test_physical_type_as_set(length, pressure):
    """Test the `as_set` attribute of `physical.PhysicalType`."""
    assert length.as_set == {"length"}
    assert pressure.as_set == {"energy density", "pressure"}


def test_physical_type_id(time):
    """
    Test that the physical type ID of a `physical.PhysicalType`
    instance matches that of a unit of that physical type.
    """
    assert time._get_physical_type_id() == u.hr._get_physical_type_id()


equivalent_unit_pairs = [
    (u.m, u.m),
    (u.m, u.cm),
    (u.m ** 18, u.pc ** 18),
    (u.N, u.kg * u.m * u.s ** -2),
    (u.barn * u.Mpc, u.cm ** 3),
]


@pytest.mark.parametrize("unit1, unit2", equivalent_unit_pairs)
def test_physical_type_id(unit1, unit2):
    """Test that the physical type IDs of equivalent units match."""
    id1 = physical.get_physical_type(unit1)._get_physical_type_id()
    id2 = physical.get_physical_type(unit2)._get_physical_type_id()
    assert id1 == id2


@pytest.mark.parametrize("unit1, unit2", equivalent_unit_pairs)
def test_physical_type_instance_equality(unit1, unit2):
    """
    Test that `physical.PhysicalType` instances for units of the same
    dimensionality are equal.
    """
    physical_type1 = physical.PhysicalType(unit1, "ptype1")
    physical_type2 = physical.PhysicalType(unit2, "ptype2")
    assert (physical_type1 == physical_type2) is True
    assert (physical_type1 != physical_type2) is False


def test_physical_type_subset_equality(pressure):
    """
    Test that equality of a `physical.PhysicalType` instance works for
    when the other object is a subset of the
    """
    # This test could be updated later for a case where there are three
    # or more physical type names associated with a particular set of
    # units.
    assert pressure == {"energy density"}
    assert pressure == {"pressure"}


nonequivalent_unit_pairs = [
    (u.m, u.s),
    (u.m ** 18, u.m ** 19),
    (u.N, u.J),
    (u.barn, u.imperial.deg_F),
]


@pytest.mark.parametrize("unit1, unit2", nonequivalent_unit_pairs)
def test_physical_type_instance_inequality(unit1, unit2):
    """
    Test that `physical.PhysicalType` instances for units with different
    dimensionality are considered unequal
    """
    physical_type1 = physical.PhysicalType(unit1, "ptype1")
    physical_type2 = physical.PhysicalType(unit2, "ptype2")
    assert (physical_type1 != physical_type2) is True
    assert (physical_type1 == physical_type2) is False


def test_physical_type_inequality_other_type():
    """
    Test that a physical type does not turn out to be equal to an object
    of another type.
    """
    assert u.m.physical_type != u.m


def test_physical_type_contains(length):
    """Test PhysicalType.__contains__."""
    assert "length" in length


def test_physical_type_multiplication(length, time, speed):
    """Test PhysicalType.__mul__."""
    assert speed * time == length
    assert speed * u.s == length


def test_physical_type_division(length, time, speed):
    """Test PhysicalType.__div__."""
    assert length / time == speed
    assert length / u.s == speed
    assert u.m / time == speed


def test_physical_type_power(length, area):
    """Test taking `physical.PhysicalType` instances to different powers."""
    assert length ** 2 == area
    assert area ** (1 / 2) == length


def test_physical_type_power_typeerror(length):
    """
    Test that exponentiation with an invalid power raises the
    appropriate exception.
    """
    with pytest.raises(TypeError):
        length ** length


def test_dimensionless(dimensionless):
    """Test a dimensionless instance of `physical.PhysicalType`."""
    assert dimensionless == "dimensionless"
    assert dimensionless == {"dimensionless"}


def test_str_for_unique_physical_type(length, speed):
    """
    Test that a `physical.PhysicalType` instance gets converted to a
    string with `str` correctly.
    """
    assert str(length) == "length"
    assert str(speed) == "speed"


def test_repr_for_unique_physical_type(length, speed):
    """
    Test that units with only one element in the set of physical types
    gets represented correctly as a string with `repr`.  This behavior
    should mimic a string to match the previous API.
    """
    assert repr(length) == "'length'"
    assert repr(speed) == "'speed'"


def test_str_and_repr_for_multiple_physical_types(pressure):
    """
    Test that units with multiple physical types get represented
    correctly as strings with both `str` and `repr`.
    """
    expected = "{'energy density', 'pressure'}"
    assert str(pressure) == expected
    assert repr(pressure) == expected


def test_physical_type_attribute_of_unit():
    """
    Test that a `u.Unit` instance has an attribute named `physical_type`
    which is an instance of `physical.PhysicalType`.
    """
    assert isinstance(u.m.physical_type, physical.PhysicalType)


def test_unknown_unit_physical_type_is_PhysicalType():
    """
    Test that a set of units with an unknown physical type has an
    attribute named `physical_type` that is a `physical.PhysicalType`
    instance.
    """
    unknown_unit = u.s ** 19
    assert isinstance(unknown_unit.physical_type, physical.PhysicalType)


@pytest.mark.parametrize("unit_with_physical_type", [u.m, u.C, u.m / u.s ** 2])
def test_redundant_physical_type(unit_with_physical_type):
    """
    Test that trying to redefine a physical type with
    `physical.def_physical_type` raises the appropriate exception.
    """
    with pytest.raises(ValueError):
        physical.def_physical_type(unit_with_physical_type, "redundant")


@pytest.mark.parametrize("not_unit", [u.m.physical_type, 1, "1"])
def test_get_physical_type_error(not_unit):
    """
    Test that attempting to get a physical type with
    `physical.get_physical_type` for something that is not a unit raises
    the appropriate exception.
    """
    with pytest.raises(TypeError):
        physical.get_physical_type(not_unit)


@pytest.mark.parametrize(
    "temperature_unit1, temperature_unit2",
    [
        (u.K, u.deg_C),
        (u.K, u.imperial.deg_F),
        (u.deg_C, u.imperial.deg_F),
    ],
)
def test_temperature_kelvin_celsius(temperature_unit1, temperature_unit2):
    """
    Because K, °C, & °F have different physical type IDs, test that
    different measurements of temperature are treated as equivalent.
    """
    assert temperature_unit1.physical_type == "temperature"
    assert temperature_unit2.physical_type == "temperature"
    assert temperature_unit1.physical_type == temperature_unit2.physical_type
