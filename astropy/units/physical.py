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

    Parameters
    ----------
    unit : u.Unit
        The unit to be represented by the physical type.

    physical_types : `str` or `set` of `str`
        A `str` representing the name of the physical type of the unit,
        or a `set` containing strings that represent one or more names
        of physical types.

    Examples
    --------
    The preferred method to access a physical type is by accessing the
    ``physical_type`` attribute of a unit.

    >>> import astropy.units as u
    >>> u.meter.physical_type
    'length'

    Occasionally, a set of units will correspond to more than one
    physical type.

    >>> u.Pa.physical_type
    {'energy density', 'pressure'}
    >>> u.Pa.physical_type == 'pressure'
    True
    >>> 'energy density' in u.Pa.physical_type
    True

    Physical types may be used for dimensional analysis.

    >>> length = u.pc.physical_type
    >>> area = (u.cm ** 2).physical_type
    >>> length * area
    'volume'

    Multiplcation, division, and exponentiation work with physical types.

    >>> length * area
    'volume'
    >>> area / length
    'length'
    >>> length ** 3
    'volume'

    Unknown physical types will be labelled as ``"unknown"``.

    >>> (u.s ** 13).physical_type
    'unknown'

    Dimensional analysis may be performed for unknown physical types too.

    >>> length_to_19th_power = (u.m ** 19).physical_type
    >>> length_to_20th_power = (u.m ** 20).physical_type

    >>> length_to_19th_power, length_to_20th_power
    ('unknown', 'unknown')

    >>> length_to_20th_power / length_to_19th_power
    'length'

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

    @property
    def as_string(self) -> str:
        """
        Return a string representation of the physical type(s).
        """
        if len(self.as_set) == 1:
            return list(self.as_set)[0]
        else:
            return "{" + str(sorted(self.as_set))[1:-1] + "}"

    def _get_physical_type_id(self):
        # To match unit API
        return self._physical_type_id

    def __eq__(self, other):
        """Return `True` if `other` represents a physical type"""
        if isinstance(other, str):
            return other in self.as_set or other.replace("_", " ") in self.as_set
        elif isinstance(other, set):
            return other.issubset(self.as_set)
        elif isinstance(other, PhysicalType):
            if self._get_physical_type_id() == other._get_physical_type_id():
                return True
            elif self == "temperature" and other == "temperature":
                # Because K, °C, & °F have different physical type IDs
                return True
            else:
                return False
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, item):
        return item in self.as_set

    def __repr__(self):
        if len(self.as_set) == 1:
            return repr(self.as_string)
        else:
            return self.as_string

    def __str__(self):
        return self.as_string

    @staticmethod
    def _identify_unit_from_unit_or_physical_type(obj):
        """
        If a unit is passed in, return that unit.  If a physical type
        is passed in, return a unit that corresponds to that physical
        type.
        """
        if isinstance(obj, core.UnitBase):
            return obj
        elif isinstance(obj, PhysicalType):
            return obj._unit
        else:
            raise TypeError("Expecting a unit or a physical type")

    def __mul__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        new_unit = self._unit * other_unit
        return new_unit.physical_type

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        new_unit = self._unit / other_unit
        return new_unit.physical_type

    def __rtruediv__(self, other):
        other_unit = self._identify_unit_from_unit_or_physical_type(other)
        new_unit = other_unit / self._unit
        return new_unit.physical_type

    def __pow__(self, power):
        if not isinstance(power, numbers.Real):
            raise TypeError(f"{power} is not a real number")
        return (self._unit ** power).physical_type


def def_physical_type(unit: core.UnitBase, name: Union[str, Set[str]]):
    """
    Add a mapping between a unit and the corresponding physical type(s).

    Parameters
    ----------
    unit : u.Unit
        The unit to be represented by the physical type.

    name : `str` or `set` of `str`
        A `str` representing the name of the physical type of the unit,
        or a `set` containing strings that represent one or more names
        of physical types.
    """

    if name in ("unknown", {"unknown"}):
        raise ValueError("Unable to uniquely define an unknown physical type")

    physical_type_id = unit._get_physical_type_id()
    if physical_type_id in _physical_unit_mapping:
        raise ValueError(
            f"{physical_type_id!r} ({name!r}) already defined as "
            f"{_physical_unit_mapping[physical_type_id]!r}"
        )

    if isinstance(name, str):
        name = {name}

    _physical_unit_mapping[physical_type_id] = PhysicalType(unit, name)

    for physical_type in name:
        _unit_physical_mapping[physical_type] = physical_type_id


def get_physical_type(unit: core.UnitBase) -> PhysicalType:
    """
    Return a representation of the physical type of a unit.

    Parameters
    ----------
    unit : `~astropy.units.UnitBase` instance
        The unit to look up

    Returns
    -------
    physical : PhysicalType
        A representation of the physical type(s) of the unit.
    """
    if not isinstance(unit, core.UnitBase):
        raise TypeError("The input to get_physical_type must be a unit.")

    physical_type_id = unit._get_physical_type_id()

    known_unit = physical_type_id in _physical_unit_mapping

    if known_unit:
        return _physical_unit_mapping[physical_type_id]
    else:
        return PhysicalType(unit, "unknown")


units_and_physical_types = [
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
    (si.Pa, {"pressure", "energy density"}),
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
]

for unit, physical_type in units_and_physical_types:
    def_physical_type(unit, physical_type)
