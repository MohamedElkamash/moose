# ThermalFunctionSolidProperties

!syntax description /Modules/SolidProperties/ThermalFunctionSolidProperties

## Description

The `ThermalFunctionSolidProperties` class provides
thermal properties for an arbitrary solid a function of temperature.
Parsed function inputs are provided for density, thermal conductivity, and
specific heat by parameterizing the time variable `t` as temperature.
This user object can also be used to specify constant properties by
excluding any dependence on temperature.

## Range of Validity

The range of validity of this material depends on the correlations provided
by the user. Note that arbitrary units can be specified with this material.

## Example Input File Syntax

!listing modules/solid_properties/test/tests/functional/test.i block=Modules

!syntax parameters /Modules/SolidProperties/ThermalFunctionSolidProperties

!syntax inputs /Modules/SolidProperties/ThermalFunctionSolidProperties

!syntax children /Modules/SolidProperties/ThermalFunctionSolidProperties

!bibtex bibliography
