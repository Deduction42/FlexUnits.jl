[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/Deduction42/FlexUnits.jl/badge.svg?branch=main)](https://coveralls.io/github/Deduction42/FlexUnits.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://deduction42.github.io/FlexUnits.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://deduction42.github.io/FlexUnits.jl/dev)

# FlexUnits.jl
FlexUnits.jl is a rewrite of Unitful.jl that maintains similar performance to Unitful.jl when units can be statically inferred, but leverages techniques in DynamicQuantities.jl to eliminate many of Unitful's performance pitfalls when units are uninferrable. In addition, FlexUnits provides shortcut unit-inference methods for mixed-unit linear algebra operations; this allows us to use existing high-performance linear algebra operations on raw numbers with separate, low-overhead unit inference. This allows FlexUnits to be used with many more Julia packages including Statistics and DifferentialEquations.jl (refer to examples in the documentation). 

Through four major design decisions, FlexUnits seamlessly blends concepts from Unitful and DynamicQuantities to achieve the best of both worlds, and surpasses both packages in terms of linear algebra capability:

1. Promotion rules that convert Unitful-like `StaticDims` into DynamicQuantities-like `Dimensions` when units are mismatched or uninferrable. Operations on `Dimensions` are slower than `StaticDims` because dimension operations need to be performed at runtime; however, dynamic dimension ops are type-stable, compilable, and still much faster than dynamic dispatch that Unitful falls back to.

2. Quantities with units are converted to quantities with dimensions (i.e. SI units) before any calculation is performed. This greatly simplifies many calculation operations. When applied to `StaticDims`, this greatly reduces over-specialization, because there is only one unit for every dimension. For example, velocity is always in "m/s", temperature is always in "K", and pressure is always in "Pa"="kg/(m s²)". This decision tends to result in FlexUnits exhibiting better performance than Unitful in cases where variables are re-assigned during iteration (a common pattern in performance-sensitive code), but this comes at the cost of slightly slower unit conversions (which is more common in less performance-sensitive code).

3. The use of a sentinel value to denote an `unknown` dimension. This helps funcitonality in many cases where `zero(T)` is called (such as initializing matrices, or indexing sparse or diagonal matrices) where we want to produce a zero value with unknown units dynamic units. This value sets all dimensions to the exponent's `typemax`, and is displayed as `?/?`. This alone has been able to make certain functions like `mean` and `cov` work out of the box, where other unit packages failed.

4. Intorudcing a special array type called a `LinmapQuant`, a special matrix type of `Quantity` that is intended for linear algebra operations. Not all matrices of `Quantity` can be multiplied; if not carefully structured, multiplication will fail when summing results, as these must have the same dimension. `Quantity` matrices that can be multiplied are quantified linear mappings (hence `LinmapQuant`) that map input dimensions to output dimensions. Enforcing this structure results in a dimension-matrix that can be summarized by two vectors and a scalar, resulting in unit inference techniques that are much simpler than the linear algebra itself. For example, matrix multiplication is an O(n³) operation, but its unit inference is only O(n); matrix inversion is also O(n³) but its unit inference is only O(1). Many two-argument linear algebra functions guarantee a `LinmapQuant` result, so these matrices tend to propagate, resulting in better performance if even one matrix is a `LinmapQuant` (broadcasting is an exception). 

In addition to these design changes, there are a number of other notable differences.

#### Notable differences between FlexUnits.jl and Unitful.jl
1. The string macro `u_str` and parsing function `uparse` are not automatically exported (allowing users to export their own registries)
2. Units are not dynamically tracked through calculations, ony dimensions; calculation results are displayed as though `upreferred` was called on them.
3. The function `upreferred` is replaced by `ubase` which converts to raw dimensions (like SI)
4. Operations on affine units do not produce errors (due to automatic conversion to dimensionas). **This the correct action for the vast majority of cases, but care must be taken to make sure that affine differences such as ***temperature differences*** are in absolute units.** For example, try running follwing commands: 
    - ```(5u"°C" - 2u"°C") == 3u"°C"``` 
    - ```(5u"°C" - 2u"°C") == 3u"K"```
5. FlexUnits registries are somewhat simpler; much like Unitful, a registry is a module that contains units. However, because dynamic units are type-stable, they can all be stored efficiently inside a single dictionary. A FlexUnits registry exports string macros which, at parse time, looks up the units inside its own internal dictionary and substitutes them into the string expression. 
6. `Quantity` in FlexUnits.jl does not subtype to `Number` in order to support more value types (such as a Distribution or Array)

## General Use

### Basic Usage
Much like other unit packages, you can use string macros to build units and quantities. Unlike other packages, you must manually "use" the default registry `UnitRegistry`, this is done so as to not be overly opinionated as to what registry to use (allowing users to easily create and use their own registries instead).
```julia
julia> using FlexUnits, .UnitRegistry

julia> u = u"J/(mol*K)"
J/(mol*K)

julia> R = 8.314*u
8.314 J/(mol*K)

julia> v_satp = R*(25u"°C")/(101.3u"kPa") #Temperature is auto-converted to Kelvin
0.024470079960513324 m³/mol
```
The string macro `@u_str` produces static units, while `@ud_str` produces dynamic units; generally users want static units from string macros as promotion rules will usually promote to dynamic when dynamic units are more peformant. All mathematical operations auto-convert to SI units, including multiplication of units. Use the `quantity` funciton to bypass this behaviour
```julia
julia> 212u"°F"
373.15000000000003 K

julia> quantity(212, u"°F")
212 °F

julia> 212ud"°F"
373.15000000000003 K

julia> quantity(212, ud"°F")
212 °F
```
The uconvert function will always result in the desired units. Note that much like DynamicQuantities, you can use the `|>` operator for unit conversions.
```julia
julia> uconvert(u"°F", 373.15*u"K")
212.0 °F

julia> 9u"μm/(m*K)" |> u"μm/(m*Ra)"
5.0 μm/(m*Ra)
```

### Linear algebra
Linear algebra is accellerated through `LinmapQuant` objects that define a linear mapping from input units to output units. To attach these units, simply multiply a matrix times a `UnitMap` constructer that specifies an example of the input and output units expected by a multiplication.
```julia
u = [u"kg/s", u"kW", u"rad/s", u"N/m"]
julia> M = randn(4,4) * UnitMap(u_in = u, u_out=u)
4×4 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
      0.503317        8.70546e-6 s²/m²        -0.793573 kg     0.202484 s
 -822.909 m²/s²            -0.0408469   1006.26 (m² kg)/s²  -920.258 m²/s
  0.964538 1/kg  0.00108277 s²/(m² kg)          -0.637226   0.507414 s/kg
  -0.168001 1/s      -0.000630687 s/m²     1.90672e-5 kg/s     -0.183613
```

One common pattern is to produce a matrix where each column has similar units. This can be achieved with a unit map having output units of `u""` and input units being the *inverse* of the desired units
```julia
u = [u"kg/s", u"kW", u"rad/s", u"N/m"]
julia> X = randn(300, 4) * UnitMap(u_in = inv.(u), u_out = u"")
300×4 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  -0.527674 kg/s   224.616 (m² kg)/s³   0.693502 1/s   -1.79476 kg/s²
   -1.23491 kg/s    452.61 (m² kg)/s³   -0.77325 1/s   0.683541 kg/s²
   -1.05964 kg/s   1391.06 (m² kg)/s³   -0.64679 1/s   0.828654 kg/s²
               ⋮
```
This constructor is more efficient that producing a matrix of quantities, because it is *lazy*. Only a single scalar and two vectors of dimensions are stored. Quantities are only materialized by indexing.

These matrices can be multiplied and divided throiugh typical linear algebra syntax
```julia
julia> M\X'
4×300 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
     -7.39372 kg/s         5.00905 kg/s        6.09836 kg/s       -5.85272 kg/s  …         15.8498 kg/s        -10.082 kg/s       -4.99927 kg/s      
 4074.2 (m² kg)/s³  -2533.63 (m² kg)/s³  -2455.3 (m² kg)/s³  1216.91 (m² kg)/s³     -4952.89 (m² kg)/s³  2032.26 (m² kg)/s³  1296.76 (m² kg)/s³      
       -3.3304 1/s          4.80668 1/s         4.75294 1/s        -2.29171 1/s              12.253 1/s        -5.75445 1/s        -1.69429 1/s      
     2.54502 kg/s²       0.397339 kg/s²      -1.65875 kg/s²       3.10623 kg/s²          -0.58694 kg/s²       2.04295 kg/s²       4.18943 kg/s² 
```

### Registering new units
The default unit registry exports a function `register_unit` (and by following the template, user-defined registries can do the same). With this function, you can register units using other units or quantities as follows:
```julia
using FlexUnits, .UnitRegistry

julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:
```
However, due to the nature of macros, these dictionaries are permanent. You can re-register units with the same values (so that you can re-run scripts) but changing them is not allowed.
```julia
julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:

julia> register_unit("bbl" => 22.5*u"m^3")
ERROR: PermanentDictError: Key bbl already exists. Cannot assign a different value.
```


## Benchmarks
### Static vs dynamic units
FlexUnits.jl and DynamicQuantities.jl both greatly outperform Unitful.jl when the compiler cannot infer the units.
```julia
using FlexUnits
using .UnitRegistry
import DynamicQuantities
import Unitful
using BenchmarkTools

v1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"J/kg", 1.0*Unitful.u"A/V"]
v1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"J/kg", 1.0*DynamicQuantities.u"A/V"]
v1flex = [1.0u"m/s", 1.0u"J/kg", 1.0u"A/V"]

@btime sum(x->x^0.0, $v1uni)
  8.100 μs (86 allocations: 3.92 KiB)
@btime sum(x->x^0.0, $v1dyn)
  41.717 ns (0 allocations: 0 bytes)
@btime sum(x->x^0.0, $v1flex)
  5.300 ns (0 allocations: 0 bytes)

```
Notice the 'μ' instead of the 'n' on the Unitful result. In such uninferrable cases, FlexUnits and DynamicQuantities both offer more than a 175x speedup. However, in the case where all types *can* be inferred, Unitful and FlexUnits perform better than DynamicQuantities.
```julia
t1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s"]
t1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s"]
t1flex = [1.0u"m/s", 1.0u"m/s", 1.0u"m/s"]

@btime sum(x->x^2, $t1uni)
  3.000 ns (0 allocations: 0 bytes)
@btime sum(x->x^2, $t1dyn)
  7.407 ns (0 allocations: 0 bytes)
@btime sum(x->x^2, $t1flex)
  3.000 ns (0 allocations: 0 bytes)
```
In this case, the performance boost from static inference is only ~2.5x but in more demanding cases, the boosts can be somewhat greater (rouighly 5x). While DynamicQuantities works much better than Unitful in worst-case scenarios, FlexUnits can match performance of both packages in their respective strengths. In most benchmarks, FlexUnits performance will tie with the better option of DynamicQuantities and Unitful with two notable exceptions:

1. FlexUnits performance is between Unitful and DynamicQuantities in the area of unit conversion (as FlexUnits doesn't support static unit conversion, only static dimension tracking)
2. FlexUnits outperforms both Unitful and DynamicQuantities in cases where units are staically inferrable but internal variables are repeatedly re-assigned (for example, iterative solvers that re-assign variables, as FlexUnits doesn't over-specialize on units)

More benchmarks can be accessed through the "benchmarks.jl" file in the "test" folder of this repo.

### Linear algebra
The first example consists of multiplying a 200x4 matrix by a 4x4 matrix
```julia
#Use unitless matrices as a benchmark
Nr = 200
X = randn(Nr, 4)
M = rand(4,4)

#Construct unitful matrices
uu = [Unitful.u"kg/s", Unitful.u"kW", Unitful.u"rad/s", Unitful.u"N/m"]
ut = reshape(uu, 1, :)
Xu = X.*ut
Mu = inv.(uu) .* M .* inv.(ut)

#Construct DynamicQuantity matrices
udq = [DynamicQuantities.u"kg/s", DynamicQuantities.u"kW", DynamicQuantities.u"rad/s", DynamicQuantities.u"N/m"]
udqt = reshape(udq, 1, :)
Xdq = X.*udqt
Mdq = inv.(udq) .* M .* inv.(udqt)

#Construct LinmapQuant matrices
ufq = [UnitRegistry.u"kg/s", UnitRegistry.u"kW", UnitRegistry.u"rad/s", UnitRegistry.u"N/m"]
Xfq = X*UnitMap(u_out = UnitRegistry.u"", u_in = inv.(ufq)))
Mfq = M*UnitMap(u_out = inv.(ufq), u_in=ufq))


julia> @btime X*M #No units
  700.000 ns (3 allocations: 6.35 KiB)

julia> @btime Xu*Mu #Unitful, more than 500x slower
  395.700 μs (5603 allocations: 93.83 KiB)

julia> @btime Xdq*Mdq #DynamicQuantities, about 8x slower
  5.700 μs (3 allocations: 31.34 KiB)

julia> @btime Xfq*Mfq #LinmapQuant, almsot no overhead
  710.000 ns (4 allocations: 6.41 KiB)
```
The main reason why FlexUnits.jl has nearly no overhead is that only the inner product of the units between matrices is considered. Only the first 4-element row of X and the first column of M need to be compared. Unit inference does not touch the other 199 rows of X or the other 3 colums of M.


## Interfacing with Unitful.jl
Previous versions of FlexUnits did not support static units, so an interface was provided to work with Unitful through uconvert to provide that performance boost where units could be statically inferred. However, now that FlexUnits supports static units with equivalent or better performance (including intelligent promotion to dynamic units), it is recommended to simply use FlexUnits (especially since similar method names between the two packages can lead to confusion). Nevertheless, this interface still exists to support legacy applications.
```julia
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
