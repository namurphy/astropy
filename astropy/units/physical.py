# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Defines the physical types that correspond to different units.

This module is not intended for use by user code directly.  Instead,
the physical unit name of a `Unit` can be obtained using its
`physical_type` property.
"""

import numbers
from typing import Union, Set

from . import core
from . import si
from . import astrophys
from . import cgs
from . import imperial
from . import misc


__all__ = ["def_physical_type", "get_physical_type", "PhysicalType"]


_physical_unit_mapping = {}
_unit_physical_mapping = {}


class PhysicalType:
    """
    Represents and provides information on the physical type(s) that are
    associated with a set of units.
    """

    def __init__(self, unit, physical_types: Union[str, Set[str]]):
        self._unit = unit
        self._physical_type_id = unit._get_physical_type_id()
        if isinstance(physical_types, str):
            self._as_set = {physical_types}
        elif isinstance(physical_types, set):
            self._as_set = physical_types

    @property
    def as_set(self) -> Set[str]:
        """Return a `set` of all physical types that represent a `Unit`."""
        return self._as_set

    def _get_physical_type_id(self):
        return self._physical_type_id

    def __eq__(self, other):
        """Return `True` if `other` represents a physical type"""
        if isinstance(other, str):
            return other in self.as_set
        elif isinstance(other, PhysicalType):
            return self._get_physical_type_id() == other._get_physical_type_id()
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, item):
        return item in self.as_set

    def __repr__(self):
        if len(self.as_set) == 1:
            return [p for p in self.as_set][0]
        else:
            return str(self.as_set)

    def __str__(self):
        return self.__repr__()

    @staticmethod
    def _identify_unit_from_unit_or_physical_type(obj):
        if isinstance(obj, core.UnitBase):
            return obj
        elif isinstance(obj, _PhysicalType):
            return obj._unit
        else:
            raise TypeError("Expecting a unit or a physical type")

    def __mul__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        new_unit = self._unit * other_unit
        return new_unit.physical_type

    def __truediv__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        new_unit = self._unit / other_unit
        return new_unit.physical_type

    def __pow__(self, power):
        if not isinstance(power, numbers.Real):
            raise TypeError(f"{power} is not a real number")
        return (self._unit ** power).physical_type


def def_physical_type(unit, name):
    """
    Adds a new physical unit mapping.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance
        The unit to map from.

    name : str
        The physical type of the unit.
    """
    r = unit._get_physical_type_id()
    if r in _physical_unit_mapping:
        raise ValueError(
            f"{r!r} ({name!r}) already defined as {_physical_unit_mapping[r]!r}"
        )
    _physical_unit_mapping[r] = name
    _unit_physical_mapping[name] = r


def get_physical_type(unit):
    """
    Given a unit, returns the name of the physical type it represents.
    If it represents an unknown physical type, `"unknown"` is returned.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance
        The unit to get the physical type for

    Returns
    -------
    physical_type : str
        The name of the physical type, or `"unknown"` if not known.
    """
    r = unit._get_physical_type_id()
    return _physical_unit_mapping.get(r, "unknown")


for unit, name in [
    (core.Unit(1), "dimensionless"),
    (si.m, "length"),
    (si.m ** 2, "area"),
    (si.m ** 3, "volume"),
    (si.s, "time"),
    (si.rad, "angle"),
    (si.sr, "solid angle"),
    (si.m / si.s, "speed"),
    (si.m / si.s ** 2, "acceleration"),
    (si.Hz, "frequency"),
    (si.g, "mass"),
    (si.mol, "amount of substance"),
    (si.K, "temperature"),
    (si.deg_C, "temperature"),
    (imperial.deg_F, "temperature"),
    (si.N, "force"),
    (si.J, "energy"),
    (si.Pa, "pressure"),
    (si.W, "power"),
    (si.kg / si.m ** 3, "mass density"),
    (si.m ** 3 / si.kg, "specific volume"),
    (si.mol / si.m ** 3, "molar concentration"),
    (si.kg * si.m / si.s, "momentum/impulse"),
    (si.kg * si.m ** 2 / si.s, "angular momentum"),
    (si.rad / si.s, "angular speed"),
    (si.rad / si.s ** 2, "angular acceleration"),
    (si.g / (si.m * si.s), "dynamic viscosity"),
    (si.m ** 2 / si.s, "kinematic viscosity"),
    (si.m ** -1, "wavenumber"),
    (si.A, "electrical current"),
    (si.C, "electrical charge"),
    (si.V, "electrical potential"),
    (si.Ohm, "electrical resistance"),
    (si.S, "electrical conductance"),
    (si.F, "electrical capacitance"),
    (si.C * si.m, "electrical dipole moment"),
    (si.A / si.m ** 2, "electrical current density"),
    (si.V / si.m, "electrical field strength"),
    (si.C / si.m ** 2, "electrical flux density"),
    (si.C / si.m ** 3, "electrical charge density"),
    (si.F / si.m, "permittivity"),
    (si.Wb, "magnetic flux"),
    (si.T, "magnetic flux density"),
    (si.A / si.m, "magnetic field strength"),
    (si.H / si.m, "electromagnetic field strength"),
    (si.H, "inductance"),
    (si.cd, "luminous intensity"),
    (si.lm, "luminous flux"),
    (si.lx, "luminous emittance/illuminance"),
    (si.W / si.sr, "radiant intensity"),
    (si.cd / si.m ** 2, "luminance"),
    (astrophys.Jy, "spectral flux density"),
    (cgs.erg / si.angstrom / si.cm ** 2 / si.s, "spectral flux density wav"),
    (astrophys.photon / si.Hz / si.cm ** 2 / si.s, "photon flux density"),
    (astrophys.photon / si.AA / si.cm ** 2 / si.s, "photon flux density wav"),
    (astrophys.R, "photon flux"),
    (misc.bit, "data quantity"),
    (misc.bit / si.s, "bandwidth"),
    (cgs.Franklin, "electrical charge (ESU)"),
    (cgs.statampere, "electrical current (ESU)"),
    (cgs.Biot, "electrical current (EMU)"),
    (cgs.abcoulomb, "electrical charge (EMU)"),
]:
    def_physical_type(unit, name)
