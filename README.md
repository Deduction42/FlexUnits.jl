[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/Deduction42/FlexUnits.jl/badge.svg?branch=main)](https://coveralls.io/github/Deduction42/FlexUnits.jl?branch=main)

# FlexUnits.jl
FlexUnits.jl is a rewrite of Unitful.jl that maintains similar performance to Unitful.jl when units can be statically inferred, but leverages techniques in DynamicQuantities.jl to eliminate many of Unitful's performance pitfalls when units are uninferrable. FlexUnits seamlessly blends concepts from Unitful and DynamicQuantities to achieve the best of both worlds through two major design decisions:

1. Promotion rules that convert Unitful-like `StaticDims` into DynamicQuantities-like `Dimensions` when units are mismatched or uninferrable. Operations on `Dimensions` are slower than `StaticDims` because dimension operations need to be performed at runtime; however, dynamic dimension ops are type-stable, compilable, and still much faster than dynamic dispatch that Unitful falls back to.

2. Quantities with units are converted to quantities with dimensions (i.e. SI units) before any calculation is performed. This greatly simplifies many calculation operations. When applied to `StaticDims`, this greatly reduces over-specialization, because there is only one unit for every dimension. For example, velocity is always in "m/s", temperature is always in "K", and pressure is always in "Pa"="kg/(m s²)". This decision tends to result in FlexUnits exhibiting better performance than Unitful in cases where variables are re-assigned during iteration (a common pattern in performance-sensitive code), but this comes at the cost of slower unit conversions (which is more common in less performance-sensitive code).

In addition to these design changes, there are a number of other notable differences.

#### Notable differences between FlexUnits.jl and Unitful.jl
1. The string macro `u_str` and parsing function `uparse` are not automatically exported (allowing users to export their own registries)
2. Units are not dynamically tracked through calculations, ony dimensions; calculation results are displayed as though `upreferred` was called on them.
3. The function `upreferred` is replaced by `ubase` which converts to raw dimensions (like SI)
4. Operations on affine units do not produce errors (due to automatic conversion to dimensionas). **This the correct action for the vast majority of cases, but care must be taken to make sure that affine differences such as ***temperature differences*** are in absolute units.** For example, try running follwing commands: 
    - ```(5u"°C" - 2u"°C") == 3u"°C"``` 
    - ```(5u"°C" - 2u"°C") == 3u"K"```
5. FlexUnits registries are much simpler; a registry a module that contains a dict of uniformly typed dynamic units, and exports string macros. Custom registries inside user-defined modules are not neccessary, but are still supported (in case of naming convention conflicts).
6. `Quantity` in FlexUnits.jl does not subtype to `Number` in order to support more value types (such as a Distribution or Array)

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
v1flex = [1.0u"m/s", 1.0u"J/kg", 1.0u"A/V"]

@btime sum(x->x^0.0, $v1uni)
  7.575 μs (86 allocations: 3.92 KiB)
@btime sum(x->x^0.0, $v1dyn)
  41.667 ns (0 allocations: 0 bytes)
@btime sum(x->x^0.0, $v1flex)
  27.209 ns (0 allocations: 0 bytes)

```
Notice the 'μ' instead of the 'n' on the Unitful result. In such uninferrable cases, FlexUnits and DynamicQuantities both more than a 175x speedup. In the case where all types *can* be inferred, Unitful and FlexUnits perform better than DynamicQuantities.
```
t1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s"]
t1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s"]
t1flex = [1.0u"m/s", 1.0u"m/s", 1.0u"m/s"]

@btime sum(x->x*x, $t1uni)
  2.900 ns (0 allocations: 0 bytes)
@btime sum(x->x*x, $t1dyn)
  7.800 ns (0 allocations: 0 bytes)
@btime sum(x->x*x, $t1flex)
  2.900 ns (0 allocations: 0 bytes)
```
In this case, the performance boost from static inference is only ~2.5x but in more demanding cases, the boosts can be greater (rouighly 5x). While DynamicQuantities can definitely smooth out performance issues in worst-case scenarios, FlexUnits can match it while also being able to leverage performance tricks for the statically inferrable cases. In most benchmarks, FlexUnits performance will tie with the better option of DynamicQuantities and Unitful with two notable exceptions:

1. FlexUnits performance is between Unitful and DynamicQuantities in the area of unit conversion (as FlexUnits doesn't support static unit conversion, only static dimension tracking)
2. FlexUnits outperforms both Unitful and DynamicQuantities in cases where units are staically inferrable but internal variables are repeatedly re-assigned (for example, iterative solvers that re-assign variables, as FlexUnits doesn't over-specialize on units)


## Interfacing with Unitful.jl
Previous versions of FlexUnits did not support static units, so an interface was provided to work with Unitful through uconvert to provide that performance boost where units could be statically inferred. However, now that FlexUnits supports static units with equivalent or better performance (including intelligent promotion to dynamic units), it is recommended to simply use FlexUnits (especially since similar method names between the two packages can lead to confusion). Nevertheless, this interface still exists to support legacy applications.
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