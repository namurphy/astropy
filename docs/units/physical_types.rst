.. |PhysicalType| replace:: :class:`~astropy.units.PhysicalType`
.. |Quantity| replace:: :class:`~astropy.units.Quantity`

Physical Types
**************

A physical type is..

Accessing Physical Types
========================

.. EXAMPLE START: Getting Physical Types

There are two ways to access |PhysicalType| instances. The physical
type of a unit can be accessed via the `~astropy.units.UnitBase.physical_type`
attribute.

  >>> import astropy.units as u
  >>> u.coulomb.physical_type
  PhysicalType('electrical charge')
  >>> (u.meter ** 3).physical_type

Physical types may be inferred from a wider variety of objects by using
`~astropy.units.get_physical_type`.

  >>> u.get_physical_type()

Ways to use ``get_physcial_type``.sjkldfklfdskl;as

* unit
* quantity
* object that can become a quantity (dimensionless)
* real number (dimensionless)
* string 



.. EXAMPLE END

Physical Type Operations
========================


Dimensional Analysis
====================

.. EXAMPLE START: Dimensional Analysis With Physical Types

`~astropy.units.PhysicalType` instances are useful for dimensional
analysis.

  >>> from astropy import units as u
  >>> length = u.get_physical_type("length")
  >>> time = u.get_physical_type("time")
  >>> length * time


.. EXAMPLE END


Available Physical Types
========================

Defining New Physical Types
===========================

On occasion


