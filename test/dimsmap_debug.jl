using FlexUnits, .UnitRegistry

Xraw = randn(50,5)*rand(5,5)
Uraw = UnitMap(u_out = u"", u_in = inv.([u"m^3", u"kg", u"kJ", u"kPa", u"rad/s"]))
X = Xraw * Uraw
