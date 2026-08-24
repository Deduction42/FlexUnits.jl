[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Build Status](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Deduction42/FlexUnits.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/Deduction42/FlexUnits.jl/badge.svg?branch=main)](https://coveralls.io/github/Deduction42/FlexUnits.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://deduction42.github.io/FlexUnits.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://deduction42.github.io/FlexUnits.jl/dev)

# FlexUnits.jl
FlexUnits.jl is the Julia units package for demanding users who want the low-level static performance of Unitful.jl, the type-stable high-level dynamic performance of DynamicQuantities.jl. Additionally, FlexUnits also enables high-performance linear algebra operations on mixed-unit matrices and extends mixed-unit compatibility for Julia packages like Statistics.jl and DifferentialEquations.jl (something that other unit packages struggle with). 

To get started, simply run
```
import Pkg; Pkg.add("FlexUnits")
using FlexUnits, .UnitRegistry
```
Note that the unit registry isn't exported by default; it must be manually imported (this allows users to create their own custom registries and import them instead). 

## Design Philosophy
This packages is built around the following major design decisions:

1. Mixed functionality between Unitful-like `StaticDims` and DynamicQuantities-like `AbstractDimensions`. Promotion rules convert StaticDims to their dynamic counterparts when units are mismatched or uninferrable allowing for a type-stable dynamic fallback. Parameterization and conversion rules can validate dynamic `AbstractDimensions` to catch dimension errors early and convert them to the static counterparts.

2. Quantities are converted to base units before any calculations are performed. This achieves the run-time performance of DynamicQuantities.jl and greatly improves compile-time performance over Unitful.jl.

3. Unit registries consist of a type-stable dictionary living inside a module that exports standard functionality, including the ability to *register new units at runtime*. A default registry is supplied but FlexUnits but registries are designed to be interchangeable; advanced functionality, such as defining new dimensions or transformations requires building a new registry.

4. Introducing a special array type called a `LinmapQuant`, a special matrix type of `Quantity` that is intended for linear algebra operations on mixed-dimensions. If a mixed-dimensional `Quantity` matrix is valid for matrix multiplication, it is a *linear mapping* (hence `LinmapQuant`) which maps input dimensions to output dimensions. This structure is used to enforce validity of matrix operations and reduce unit inference overhead to linear time or better.

5. The use of a sentinel value to denote an `unknown` dimension displayed as `?/?`. This enhances compatibility for Julia codebases where `zero(T)` is called without foreknowledge of units (such as initializing mixed-unit matrices, or indexing sparse/diagonal matrices). 

In addition to these design changes, there are a number of other notable differences from Unitful.

#### Notable differences between FlexUnits.jl and Unitful.jl
1. FlexUnits registries are designed to be interchangeable and consist of an dictionary of dynamic units living inside a module that exports string macros and parsing functions. String macros can produce static units while string parsing functions produce type-stable dynamic units.
2. The string macro `u_str` and parsing function `uparse` are not automatically exported, but must be exported by a chosen unit registry (allowing users to export their own registries)
3. Only dimensions are tracked through calculations and results are displayed as though `upreferred` was called on them. More intuitive representations can be obtained using `simplify(q)` or setting `display_simplified_units(true)`.
4. The function `upreferred` is replaced by `ubase` which converts quantities to base units.
5. Operations on affine units do not produce errors (due to automatic conversion to dimensions). **This the correct action for the vast majority of cases, but care must be taken to make sure that affine differences such as ***temperature differences*** are in absolute units.** For example, try running following commands:
    - ```(5u"°C" - 2u"°C") == 3u"°C"```
    - ```(5u"°C" - 2u"°C") == 3u"K"```
6. Much like Unitful, `Quantity` subtypes to number, but an additional type `FlexQuant` can support any value type (such as a Distribution or Array). The function `quantity(q, u)` selects the appropriate output type based on the arguments.
7. FlexUnits uses the concept of a `LogQuant` to allow taking the `log` of a `Quantity`, and handling logarithmic units like decibels

## Basic Usage
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
The string macro `@u_str` produces static units, while `@ud_str` produces dynamic units; generally users want static units from string macros as promotion rules will usually promote to dynamic when dynamic units are more performant. All mathematical operations auto-convert to SI units, including multiplication of units. Use the `quantity` function to bypass this behaviour.
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
The `uconvert` function will always result in the desired units. Note that much like DynamicQuantities, you can use the `|>` operator for unit conversions.
```julia
julia> uconvert(u"°F", 373.15*u"K")
212.0 °F

julia> 9u"μm/(m*K)" |> u"μm/(m*Ra)"
5.0 μm/(m*Ra)
```

## Registering New Units
The default unit registry exports a function `register_unit` (and by following the template, user-defined registries can do the same). The unit registration process is much simpler than Unitful.jl as you only need to call this function; *no additional registries or modules are required* and you can even register units at run-time.
```julia
using FlexUnits, .UnitRegistry

julia> register_unit("bbl" => 0.158987*u"m^3") # US barrel of oil
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:
```
However, due to the nature of macros, these dictionaries are permanent. You can re-register units with the same values (so that you can re-run scripts) but changing them is not allowed.
```julia
julia> register_unit("bbl" => 0.158987*u"m^3") # US barrel of oil
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:

julia> register_unit("bbl" => 0.164*u"m^3") # UK barrel of beer
ERROR: PermanentDictError: Key bbl already exists. Cannot assign a different value.
```

## Dimensional Enforcement and Dispatch
FlexUnits simplifies dispatch and dimensional enforcement with the `@D_str` macro exported from the unit registry; this macro produces a static dimension type given input units.
```julia
julia> D"kPa"
StaticDims{kg/(m s²)}
```
This allows you to parameterize quantities, for example:
```julia
julia> const PressureQuant{T} = Quantity{T, D"Pa"}
Quantity{T, StaticDims{kg/(m s²)}} where T
```
This practice comes with two major benefits:
1.  It enforces unit consistency early on in the process, catching errors near the source 
```julia
julia> PressureQuant{Float64}(5*uparse("kg"))
ERROR: ConversionError: Cannot convert unit 'kg' to target unit 'kg/(m s²)' due to a dimension mismatch of '(m s²)' 
```
2. Static knowledge of the dimensions improves run-time performance
```julia
julia> typeof(5*uparse("kPa")) #Slower dynamic quantity
Quantity{Float64, Dimensions{FixRat32}}

julia> typeof(PressureQuant{Float64}(5*uparse("kPa"))) #Faster static quantity
Quantity{Float64, StaticDims{kg/(m s²)}}
```
You can also easily create dispatch patterns based off different dimensions.
```
mass_flow(vf::Quantity{<:Any, D"m^3/s"}, ρ::Quantity{<:Any, D"kg/m^3"}) = vf*ρ
mass_flow(vf::Quantity{<:Any, D"m^3/s"}, vs::Quantity{<:Any, D"m^3/kg"}) = vf/vs
```

## Unit Simplification
While operations always convert to SI units, FlexUnits contains a `simplify` function which applies a greedy algorithm to reduce the number of symbols. 

```julia
julia> r = 5u"V" / 10u"mA" #Electrical resistance is opaque
500.0 (m² kg)/(s³ A²)

julia> r = 5u"V" / 10u"mA" |> simplify #This is much better
500.0 Ω
```

This simplification can have much better success than Unitful for long operations that yield common dimensional quantities.
```julia
julia> p = 1u"kg/L"*9.18u"m/s^2"*1u"ft" |> simplify  #Hydraulic pressure
2798.0640000000003 Pa
```
This is easier to understand than the result you would get from Unitful
```julia
using Unitful
julia> p = 1u"kg/L"*9.18u"m/s^2"*1u"ft" #Hydraulic pressure
9.18 ft kg m L^-1 s^-2
```

### Switching simplification units
FlexUnits simplification works by attempting to fit a predefined set of units (prioritizing high complexity) to a dimensional value. You can easily change the units in this set to display the units you like.

```julia
set_preferred_unit(u"psi") #Pressures are now in PSI

julia> p = 1u"kg/L"*9.18u"m/s^2"*1u"ft" |> simplify  #Hydraulic pressure
0.4058247132605051 psi
```

### Simplified view by default
It may be cumbersome to constantly use `|> simplify` after every interactive operation. FlexUnits has a configuration function that allows you to view results as though `simplify` was applied to them.
```julia
display_simplified_units(true) #Show simplified results

julia> p = 1u"kg/L"*9.18u"m/s^2"*1u"ft" #This is convenient
2798.0640000000003 Pa
```
This begs the question: *Why is this not enabled by default?*. The main reason is ***because it can lie***; it displays results *as though they were simplified*, but does not actually simplify the results. This can mess up `ustrip` ***Ye be warned***
```julia
display_simplified_units(true)
set_preferred_unit(u"psi")

julia> p = 1u"kg/L"*9.18u"m/s^2"*1u"ft" #This is convenient
0.4058247132605051 psi

julia> ustrip(p) #I have been deceived
2798.0640000000003
```
Another reason not to enable by default is that ***if you wish to create custom dimensions, simplification will not work out of the box***. You will need to do some extra work to get simplification to run (see the section below or advanced examples section in the documentation for a guide).

## Custom Dimensions and Registries
While the default registry is usually sufficient for most cases, FlexUnits lets you interchangeably define new registries if the default doesn't suit your needs. One such case is if you need to define new dimensions not present in `Dimensions{FixRat32}`. For example, let us consider the case where we want to extend the dimensions to include currency in Euros. We first need to define a new `AbstractDimensions` type:

```julia
using FlexUnits
@kwdef struct MoneyDimensions{P} <: AbstractDimensions{P}
    m   ::P = zero(FixRat32)
    kg  ::P = zero(FixRat32)
    s   ::P = zero(FixRat32)
    A   ::P = zero(FixRat32)
    K   ::P = zero(FixRat32)
    cd  ::P = zero(FixRat32)
    mol ::P = zero(FixRat32)
    €   ::P = zero(FixRat32)
end
MoneyDimensions(args::Real...) = MoneyDimensions{FixRat32}(args...)
```
We then need to build out a unit registry that uses `MoneyDimensions`. Additionally, because simplification only supports `Dimensions{FixRat32}` out of the box, we need to define simplification routines for `MoneyDimensions` in our new registry.
```julia
module CurrencyUnits
    using FlexUnits.RegistryTools
    import ..MoneyDimensions

    # Define your unit registry dict and auto-populate it
    const UNITS = PermanentDict{Symbol,Units{MoneyDimensions{FixRat32},AffineTransform{Float64}}}()
    registry_defaults!(UNITS) # registry_defaults! is compatible with MoneyDimensions because it is a superset of Dimensions
    register_unit!(UNITS, "EUR" => UNITS[:€]) # additional registrations can happen here

    # Optional support for simplification but recommended if dimensions are not Dimensions{FixRat32}
    const PREFERRED_UNITS = [UNITS[u] for u in [:F, :H, :T, :Ω, :V, :W, :J, :Pa, :N, :C, :L]] # typical simplification basis you can change at will
    @generate_unit_simplifier(PREFERRED_UNITS) # takes care of simplification boilerplate code

    @generate_registry_exports(UNITS) # takes care of typical string macro and unit parsing exports
end
```
We can now use this new module to support string macro parsing and simplification.
```julia
using .CurrencyUnits

julia> electricity_price = 0.195u"€/(kW*hr)"
5.416666666666666e-8 (s² €)/(m² kg)

julia> electricity_price = 0.195u"€/(kW*hr)" |> simplify
5.416666666666666e-8 €/J
```
Admittedly, this process is more involved than registering a new dimension in Unitful.jl, but it is on par with the effort required for registering new units in the same package. FlexUnits makes a rare workflow (registering dimensions) more difficult in exchange for more performance and simplicity in common workflows (type-stable unit parsing and registering units).

## Mixed-unit linear algebra
Linear algebra is accelerated through `LinmapQuant` objects that define a linear mapping from input units to output units. To attach these units, simply multiply a matrix times a `UnitMap` constructor that specifies an example of the input and output units expected by a multiplication.
```julia
u = [u"kg/s", u"kW", u"rad/s", u"N/m"]

julia> M = randn(4,4) * UnitMap(u_in = u, u_out=u)
4×4 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
     -0.723157        0.000463328 s²/m²        0.0841925 kg     0.641013 s
 -29.6711 m²/s²             -0.0910066   -308.33 (m² kg)/s²  -354.274 m²/s
  0.190549 1/kg  0.000153441 s²/(m² kg)            0.91125   -1.04493 s/kg
  -0.940626 1/s         0.00182779 s/m²       0.472231 kg/s     -0.226755 
```

One common pattern is to produce a matrix where each column has similar units. This can be achieved the same way as above or by horizontally concatenating vectors multiplied by units (do NOT broadcast the multiplication, that will create a regular matrix)
```julia
N = 300

julia> X = [randn(N)*u"kg/s" randn(N)*u"kW" randn(N)*u"rad/s" rand(N)*u"N/m"]
300×4 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, StaticArraysCore.SVector{4, Dimensions{FixRat32}}, Vector{StaticDims{}}}}:
  -1.06923 kg/s  -1551.15 (m² kg)/s³    0.37945 1/s   0.39719 kg/s²
 -0.277549 kg/s  -705.525 (m² kg)/s³   0.337886 1/s  0.655497 kg/s²
   1.04352 kg/s  -604.187 (m² kg)/s³   -1.35279 1/s  0.981266 kg/s²
              ⋮    
```

These matrices can be multiplied and divided through typical linear algebra syntax
```julia
julia> M\X'
4×300 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{StaticDims{}}, Vector{Dimensions{FixRat32}}}}:
       5.36796 kg/s        2.10869 kg/s       0.482899 kg/s      0.0537144 kg/s  …        3.73209 kg/s        5.55815 kg/s        3.39604 kg/s       1.10409 kg/s
 2925.74 (m² kg)/s³  1373.05 (m² kg)/s³  984.784 (m² kg)/s³  522.156 (m² kg)/s³     2727.62 (m² kg)/s³  3370.72 (m² kg)/s³  2145.29 (m² kg)/s³  1529.4 (m² kg)/s³
        1.22354 1/s        0.688339 1/s      -0.0661644 1/s         -1.5527 1/s           -1.10673 1/s       -0.709414 1/s        0.142878 1/s       -1.29612 1/s
      2.11238 kg/s²      0.863076 kg/s²       1.46959 kg/s²      -1.03138 kg/s²          1.94117 kg/s²       1.02301 kg/s²       1.68966 kg/s²     0.699505 kg/s²
```

## Logarithmic Quantities and Units
Handling logarithmic units such as decibels comes with some level of controversy. The question is whether a quantity such as `1 dB` is a logarithm of a quantity or merely the *logarithmic representation* of a linear quantity.
1. If `1 dB` is a logarithmic representation: `1 dB + 1 dB = 4.0103 dB`  
2. If `1 dB` is a logarithmic quantity: `1 dB + 1 dB = 2 dB`

FlexUnits adopts the philosophy of the second camp, where a decibel represents the *logarithm of a quantity* and supplies algebraic tools to manipulate logarithms of quantities (referred to as a `LogQuant`). 

### Producing logarithmic quantities
One way to produce a `LogQuant` is by taking a log of a quantity.
```julia
julia> q = log(2u"W")
log(2.0 (m² kg)/s³)
```
Another way is to multiply a number by a logarithmic unit. For example, `dB` is a `LogScale` object that can be imported to construct a logarithmic unit
```julia
import FlexUnits.dB
julia> q = 30dB(u"W")
log(1000.0000000000016 (m² kg)/s³)
```
As you can see here, `30 dB(W)` is equivalent to `1000 W` but it displayed as its logarithm. This helps reinforce how operations are performed based on logarithmic identities. While their logarithms are displayed (to emphasize this algebra), the actual numerical value stored is the logarithmic form
```julia
julia> ustrip(log(2u"W"))
0.6931471805599453
```

### Operations on logarithmic quantities
The algebraic rules `LogQuant` are centered around logarithmic identities.
```julia
julia> log(4u"m") + log(4u"s") # log(x) + log(y) = log(x*y)
log(15.999999999999998 (m s))

julia> log(4u"m") - log(4u"s") # log(x) - log(y) = log(x/y)
log(1.0 m/s)

julia> 2log(4u"m") # nlog(x) = log(x^n)
log(15.999999999999998 m²)
```
We also make use of the ⊕ and ⊖ operators that, in this context, commonly refers to adding/subtracting linearized values and transforming back to log space. It's not exported by default as this symbol could be used by other packages to mean something else.
```julia
import FlexUnits: ⊕, ⊖
julia> log(8u"m") ⊕ log(4u"m") #Observe that the linear addition happened
log(12.0 m)

julia> log(8u"m") ⊖ log(4u"m") #Observe that linear subtraction happened
log(3.9999999999999982 m)
```
Logarithmic quantities can be converted back to regular quantities using `quantity`, `linquant`, or `exp`
```julia
julia> (linquant(log(4u"m")), quantity(log(4u"m")), exp(log(4u"m")))
(4.0 m, 4.0 m, 4.0 m)
```

### Registering logarithmic units
The default unit registry can only register affine units, this is not an issue if you are always constructing logarithmic units with the `dB` or `Np` wrappers:
```julia
using FlexUnits, .UnitRegistry
import FlexUnits: dB, Np

julia> dB_kPa = dB(u"kPa")
dB(kPa)
```
However, if you need to call `uparse` on strings that can represent logarithmic units, you will need to register them in the `LogUnitRegistry` instead. This registry can hold both affine and logarithmic units, but `uparse` can introduce performance issues because the output is a `Union`. ***WARNING, because multiplying `uparse` outputs can produce a Quantity or a LogQuant, based on the string value, it's recommended that you use explicit constructors like `quantity` or `ubase` to always produce linear quantities, or `logquant` or `logubase` always produce logarithmic units.***
```julia
using FlexUnits, .LogUnitRegistry
import FlexUnits: dB, Np

register_unit("dB_V" => dB(u"V"))

julia> 10uparse("dB_V")
log(10.000000000000002 (m² kg)/(s³ A))

julia> 10uparse("V")
10.0 (m² kg)/(s³ A)

julia> ubase(10, uparse("dB_V"))
10.000000000000002 (m² kg)/(s³ A)

julia> logubase(10, uparse("V"))
log(22026.465794806718 (m² kg)/(s³ A))

julia> 10u"dB_V"
log(10.000000000000002 (m² kg)/(s³ A))
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
In this case, the performance boost from static inference is only ~2.5× but in more demanding cases, the boosts can be somewhat greater (roughly 5×). While DynamicQuantities works much better than Unitful in worst-case scenarios, FlexUnits can match performance of both packages in their respective strengths. In most benchmarks, FlexUnits performance will tie with the better option of DynamicQuantities and Unitful with one notable exception: ***unit conversion***.

-  ***Unitful is fastest at static unit conversions***. Because it compiles both dimensions and conversion factors, Unitful outperforms FlexUnits (~10×) which compiles only dimensions and DynamicQuantities (~500×) which relies on the inefficient `SymbolicDimensions{T}` object for conversion
-  ***FlexUnits is fastest at dynamic unit conversions***. Because `FlexUnits.Units{Dimensions{T}}` is type-stable and efficient, FlexUnits outperforms Unitful (~35×) which is not dynamically type-stable, and DynamicQuantities (~40×) due to its use of `SymbolicDimensions{T}`

Dynamic unit conversion is much more useful for repeatable applications as you often don't know beforehand what units your data will be in, or what units your users will want the results in (although dimensions are often known). If input datasets are large, the performance differences can be substantial.

More benchmarks can be accessed through the "benchmarks.jl" file in the "test" folder of this repo.

### Mixed-unit linear algebra
The first example consists of multiplying a 200x4 matrix by a 4x4 matrix with mixed units
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

julia> @btime Xfq*Mfq #LinmapQuant, almost no overhead
  710.000 ns (4 allocations: 6.41 KiB)
```
The main reason why FlexUnits.jl has nearly no overhead is that only the inner product of the units between matrices is considered. Only the first 4-element rows of X and the first column of M need to be compared. Unit inference does not touch the other 199 rows of X or the other 3 columns of M.


## Interfacing with Unitful.jl
Previous versions of FlexUnits did not support static units, so an interface was provided to work with Unitful through `uconvert` to provide that performance boost where units could be statically inferred. However, now that FlexUnits supports static units with equivalent or better performance (including intelligent promotion to dynamic units), it is recommended to simply use FlexUnits (especially since similar method names between the two packages can lead to confusion). Nevertheless, this interface still exists to support legacy applications.
```julia
using Unitful
import FlexUnits
import FlexUnits.UnitRegistry
import FlexUnits.uconvert

julia> x = UnitRegistry.qparse.(["5.0 km/hr", "2.0 N", "10 °C"])

julia> velocity = uconvert(u"km/hr", x[1])
5.0 km hr^-1

julia> force = uconvert(u"N", x[2])
2.0 N

julia> temperature = uconvert(u"°F", x[3])
49.99999999999994 °F
```
