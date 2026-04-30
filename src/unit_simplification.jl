using Revise 
using FlexUnits, .UnitRegistry
import FlexUnits: _pretty_unit_pwr, usymbol, uscale, dimtype



#================================================================================================================
# Test Code
================================================================================================================#
display_simplified_units(true)


q = 1u"kg/L"*9.81u"m/s^2"*20u"cm"
q = 5u"A^2"*0.01u"ohm"
q = 5u"W"*sqrt(2u"V")