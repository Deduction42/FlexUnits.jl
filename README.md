# FlexUnits.jl
The flexible units package. FlexUnits was heavily inspired by DynamicQuantities.jl and Unitful.jl but aims to be more flexible than both. Much like DynamicQuantities.jl, it allows for a single concrete unit type to represent many different units. This means multiple quantities with different units can be stored inside a concretely-typed Array or Dict (which is not possible with Unitful.jl). However, unlike DynamicQuantities.jl, this package fully supports affine units (like °C and °F).

The main differences between FlexUnits.jl and DynamicQuantities.jl
1. Fully supports affine units (like °C and °F) and can potentially support logarithmic units (like dB)
2. The string macro `u_str` and parsing function `uparse` are not automatically exported (allowing users to export their own registries)
3. Extremely easy to build your own unit registry (allowing differnet behaviour for `u_str` and `uparse`)
4. While units are pretty-printed by default, you can enable parseable unit outputs by setting `pretty_print_units(false)` which is useful for outputting unit information in JSON
5. More closely resembles the Unitful.jl API in many ways

The main differences between FlexUnits.jl and Unitufl.jl
1. Units are not specialized (u"m/s" retgurns the same concrete type as u"°C") allowing much more flexible containers with only a slight run-time penalty
2. Units are automatically converted to dimensional form (SI) for any mathematical manipulations and tracked at the dimension level
3. Operations on affine units do not produce errors (due to automatic conversion to dimensional form). **This may may yield unituitive (but more consistent) results**.
4. Unit registries are much simpler. A registry is simply a dict of units, all of the same type, living inside a module.

