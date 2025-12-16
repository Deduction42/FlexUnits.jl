using Unitful
import FlexUnits
import FlexUnits.UnitRegistry
import FlexUnits.uconvert

x = UnitRegistry.qparse.(["5.0 km/hr", "2.0 N", "10 °C"])
velocity = uconvert(u"km/hr", x[1])
force = uconvert(u"N", x[2])
temperature = uconvert(u"°F", x[3])
x_out = [FlexUnits.Quantity(velocity), FlexUnits.Quantity(force), FlexUnits.Quantity(temperature)]