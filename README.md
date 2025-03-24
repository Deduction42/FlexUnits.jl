# FlexUnits.jl
The flexible units package. FlexUnits was heavily inspired by DynamicQuantities.jl and Unitful.jl but aims to be more flexible than both. Much like DynamicQuantities.jl, it allows for a single concrete unit type to represent many different units. This means multiple quantities with different units can be stored inside a concretely-typed Array or Dict (which is not possible with Unitful.jl). However, unlike DynamicQuantities.jl, this package fully supports affine units (like °C and °F).

#### Major differences between FlexUnits.jl and DynamicQuantities.jl
1. Fully supports affine units (like °C and °F) and can potentially support logarithmic units (like dB)
2. The string macro `u_str` and parsing function `uparse` are not automatically exported (allowing users to export their own registries)
3. Easier to build your own unit registry (allowing differnet behaviour for `u_str` and `uparse`)
4. While units are pretty-printed by default, you can enable parseable unit outputs by setting `pretty_print_units(false)` which is useful for outputting unit information in JSON
5. More closely resembles the Unitful.jl API in many ways

#### Major differences between FlexUnits.jl and Unitufl.jl
1. Units are not specialized (u"m/s" retgurns the same concrete type as u"°C") which is much faster when units cannot be inferred
2. Units are automatically converted to dimensional form (SI) for any mathematical manipulations and tracked at the dimension level
3. Operations on affine units do not produce errors (due to automatic conversion to dimensional form). **This may may yield unituitive (but more consistent) results**.
4. Unit registries are much simpler. A registry is simply a dict of units, all of the same type, living inside a module.

## Benchmarks
FlexUnits.jl and DynamicQuantities.jl both greatly outperform Unitful.jl when the compiler cannot infer the units.
```
using FlexUnits
using .UnitRegistry
import DynamicQuantities
import Unitful
using BenchmarkTools

v1flex = ubase.([1.0u"m/s", 1.0u"J/kg", 1.0u"A/V"])
v1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"J/kg", 1.0*Unitful.u"A/V"]
v1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"J/kg", 1.0*DynamicQuantities.u"A/V"]

@btime sum(x->x^0.0, v1uni)
  7.375 μs (86 allocations: 3.92 KiB)
@btime sum(x->x^0.0, v1flex)
  135.057 ns (1 allocation: 48 bytes)
@btime sum(x->x^0.0, v1dyn)
  107.281 ns (1 allocation: 48 bytes)
```
Notice the 'μ' instead of the 'n' on the Unitful result, FlexUnits offers a 50x speedup in this case (DynamicQuantities does a bit better, with a 68x speedup). In the case where all types can be inferred, performance is more or less the same, with FlexUnits being slightly worse than the others.
```
t1flex = ubase.((1.0u"m/s", 1.0u"J/kg", 1.0u"A/V"))
t1uni  = (1.0*Unitful.u"m/s", 1.0*Unitful.u"J/kg", 1.0*Unitful.u"A/V")
t1dyn  = (1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"J/kg", 1.0*DynamicQuantities.u"A/V")

@btime sum(x->x^0, t1uni)
  88.518 ns (2 allocations: 64 bytes)
@btime sum(x->x^0, t1flex)
  131.771 ns (3 allocations: 304 bytes)
@btime sum(x->x^0, t1dyn)
  85.328 ns (2 allocations: 256 bytes)
```
