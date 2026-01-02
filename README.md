[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/Deduction42/FlexUnits.jl/badge.svg?branch=main)](https://coveralls.io/github/Deduction42/FlexUnits.jl?branch=main)

# FlexUnits.jl
FlexUnits.jl is a rewrite of Unitful.jl that maintains similar performance to Unitful.jl when units can be inferred, but leverages techniques in DynamicQuantities.jl to eliminate many of Unitful's performance pitfalls when units cannot be statically inferred. Thus, this package blends concepts from Unitful and DynamicQuantities to achieve the best of both worlds. This is achieved through two major design decisions:

1. Promotion rules that convert Unitful-like `StaticDims` into DynamicQuantities-like `Dimensions` when units are mismatched. Operations on `Dimensions` quantities are slower than `StaticDims` quantities because dimension operations need to be performed at runtime; however, dynamic dimension ops are type-stable, compilable, and still much faster than dynamic dispatch that Unitful falls back to.

2. Quantities with units are converted to quantities with dimensions (i.e. SI units) before any calculation is performed. This greatly simplifies many calculation operations. When applied to `StaticDims`, this greatly reduces over-specialization, because there is only one unit for every dimension. For example, velocity is always in "m/s", temperature is always in "K", and pressure is always in "Pa"="kg/(m s²)". This decision tends to result in FlexUnits exhibiting better performance than Unitful in cases where variables are re-assigned during iteration (a common pattern in performance-sensitive code), but this comes at the cost of slower unit conversions. However, unit conversions are often only neccessary in high-level code; eager conversion to dimensional quantities removes the need for conversions in low-level code. 

#### Major differences between FlexUnits.jl and Unitufl.jl
1. Units are not specialized (u"m/s" returns the same concrete type as u"°C") which is much faster when units cannot be inferred
2. The string macro `u_str` and parsing function `uparse` are not automatically exported (allowing users to export their own registries)
3. Units are not dynamically tracked, quantities are displayed as though `upreferred` was called on them
4. The function `upreferred` is replaced by `ubase` which converts to raw dimensions (like SI) which are not configurable
5. Operations on affine units do not produce errors (due to automatic conversion to dimensional form). **This may may yield unituitive (but more consistent) results**.
6. Unit registries are much simpler; a registry is simply a dict of units, all of the same type, living inside a module. Custom registries inside user-defined modules are not neccessary, but are still supported.
7. `Quantity` in FlexUnits.jl does not subtype to `Number` in order to support more value types (such as a Distribution or Array)


## General Use
Much like other unit packages, you can use string macros to build units and quantities. Unlike other packages, you must manually "use" the default registry `UnitRegistry`, this is done so as to not be overly opinionated as to what registry to use (users can create and use their own registries instead).
```
julia> using FlexUnits, .UnitRegistry

julia> u = u"J/(mol*K)"
J/(mol*K)

julia> R = 8.314*u
8.314 J/(mol*K)

julia> v_satp = R*(25u"°C")/(101.3u"kPa") #Temperature is auto-converted to Kelvin
0.024470079960513324 m³/mol
```
You can register units using other units or quantities as follows:
```
julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:
  :Ω      => Ω
  :μs     => μs
  :μV     => μV
```
However, due to the nature of macros, these dictionaries are permanent. You can re-register units with the same values (so that you can re-run scripts) but changing them is not allowed.
```
julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:
  :Ω      => Ω
  :μs     => μs
  :μV     => μV

julia> register_unit("bbl" => 22.5*u"m^3")
ERROR: PermanentDictError: Key bbl already exists. Cannot assign a different value.
```

## Interfacing with Unitful.jl
FlexUnits.jl and Unitful.jl focus on different use cases and can be considered complementary. Unitful.jl shines when explicitly referring to units inside the code (especially at low-level operations) while FlexUnits.jl is much better at high-level operations when units tend to be unknown (for example, when parsing strings). Because of this, FlexUnits (as of version v0.2.10) provides an extension to Unitful that allows for converting between types. A useful example is converting `FlexUnits.Quanity` to a Unitful equivalent through `uconvert`. Note that because both packages have significant overlap in their function/object names, you will have to use `import` on at least one of the packages.
```
julia> using Unitful
julia> import FlexUnits
julia> import FlexUnits.UnitRegistry
julia> import FlexUnits.uconvert

julia> x = UnitRegistry.qparse.(["5.0 km/hr", "2.0 N", "10 °C"])

julia> velocity = uconvert(u"km/hr", x[1])
18.0 km hr^-1

julia> force = uconvert(u"N", x[2])
2.0 N

julia> temperature = uconvert(u"°F", x[3])
49.99999999999994 °F
```
This pattern would be useful when performing low-level calculations on `force`, `velocity` and `temperature` inside a function that takes a mixed-unit vector. Similarly, if one wishes to collect results of dissimilar units, one can simply output them as a type-stable FlexUnit vector
```
julia> x_out = [FlexUnits.Quantity(velocity), FlexUnits.Quantity(force), FlexUnits.Quantity(temperature)]
3-element Vector{FlexUnits.Quantity{Float64, FlexUnits.Dimensions{FlexUnits.FixedRational{25200, Int32}}}}:
 1.3888888888888888 m/s
 2.0 (m kg)/s²
 283.15 K
```
Using both packages together should feel natural due to their similar API and can provide the best of both worlds. However, these similarities also means that care must be taken to manually import any functions that needed from both packages.

## Benchmarks
FlexUnits.jl and DynamicQuantities.jl both greatly outperform Unitful.jl when the compiler cannot infer the units.
```
using FlexUnits
using .UnitRegistry
import DynamicQuantities
import Unitful
using BenchmarkTools

v1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"J/kg", 1.0*Unitful.u"A/V"]
v1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"J/kg", 1.0*DynamicQuantities.u"A/V"]
v1flex = ubase.([1.0u"m/s", 1.0u"J/kg", 1.0u"A/V"])

@btime sum(x->x^0.0, v1uni)
  7.850 μs (86 allocations: 3.92 KiB)
@btime sum(x->x^0.0, v1dyn)
  105.173 ns (1 allocation: 48 bytes)
@btime sum(x->x^0.0, v1flex)
  106.882 ns (1 allocation: 48 bytes)

```
Notice the 'μ' instead of the 'n' on the Unitful result, FlexUnits and DynamicQuantities both offer a ~75x speedup in this case (where unit type cannot be inferred). In the case where all types *can* be inferred, performance is more or less the same in terms of execution time (but Unitful allocates fewer bytes).
```
t1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s"]
t1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s"]
t1flex = ubase.([1.0u"m/s", 1.0u"m/s", 1.0u"m/s"])

@btime sum(x->x*x, t1uni)
  86.902 ns (1 allocation: 16 bytes)
@btime sum(x->x*x, t1dyn)
  86.472 ns (1 allocation: 48 bytes)
@btime sum(x->x*x, t1flex)
  86.260 ns (1 allocation: 48 bytes)
```