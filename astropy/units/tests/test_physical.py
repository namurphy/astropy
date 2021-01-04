# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Unit tests for the physical_type support in the units package
"""

import pytest

from astropy import units as u
from astropy.units import physical
from astropy.constants import hbar
from astropy.tests.helper import raises


unit_physical_type_pairs = [
    (u.m, "length"),
    (u.cm ** 3, "volume"),
    (u.km / u.h, "speed"),
    (u.barn * u.Mpc, "volume"),
    (u.m * u.s ** 8, "unknown"),
    (u.m / u.m, "dimensionless"),
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
    (u.kg * u.m / u.s, "momentum/impulse"),
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
    (u.lx, "luminous emittance/illuminance"),
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
    Test that the `physical_type` attribute of `Unit` objects provides
    the expected physical type for various units.
    """
    if unit.physical_type != physical_type:
        pytest.fail(
            f"{repr(unit)}.physical_type was expected to return "
            f"{repr(physical_type)}, but instead returned "
            f"{unit.physical_type}."
        )


@raises(ValueError)
def test_redundant_physical_type():
    physical.def_physical_type(u.m, "redundant")
