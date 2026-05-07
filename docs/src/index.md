```@meta
DocTestSetup = quote
    using FlexUnits, .UnitRegistry
end
```
# FlexUnits.jl
FlexUnits.jl is a unit package designed to resemble Unitful.jl, with similar performance when units can be statically inferred, but leverages techniques in DynamicQuantities.jl to eliminate many of Unitful's performance pitfalls when units are uninferrable. In addition, this package introduces special vector and matrix types that can infer output units for linear algebra operations at O(n) speeds, even if the numerical matrix operations are higher order (such as matrix inversions). When combining these techniques, this package allows the user to employ the following useful pattern:
1. Define high-level linear-algebra operations on quantities with dynamic units that are fast and type-stable even with mixed-unit types
2. Convert quantities to statically-inferred high-performance units in lower-level performance-sensitive code

## Quick start examples
The FlexUnits API is designed to resemble Unitful in a number of ways. One major difference is that string macros and parsing functions are not exported by default, instead they are exported by a unit registry. This allows users to define their own registries and use those as a basis for unit parsing. FlexUnits provides a unit registry (UnitRegistry) that can be imported to provide default unit parsing functionality.
```julia
import Pkg; Pkg.add("FlexUnits")
using FlexUnits, .UnitRegistry
```

Units are accessed primarily through string macros. The `@u_str` macro produces static units by default (for maximum performance when injected into low-level code), but if dynamic units are desired (particularly for interactive workflows) one can use the `@ud_str` macro.
```julia
julia> 1u"°C"   #@u_str produces static units (multiplication converts to base units)
274.15 K

julia> 1ud"°C"  #@ud_str produces dynamic units (multiplication converts to base units)
274.15 K

julia> quantity(1, u"°C") #Use 'quantity' to avoid eager conversion to base units
1 °C
```

Unlike macros, units that need to be parsed from strings cannot be statically inferred, so the `uparse` function will always produce dynamic units in order to make it type-stable.
```julia
julia> uparse("km/hr")  #Uparse always produces dynamic units
km/hr

julia> typeof(uparse("km/hr")) == typeof(uparse("kPa")) #Uparse is type-stable
true
```

The function `uconvert` pays attention to the unit type converted to. If the unit is static, the resulting quantity will be static, if it's dynamic, the result will be dynamic.
```julia
julia> typeof(uconvert(u"°F", 1u"°C")) #Converting to static units produces a static quantity
Quantity{Float64, Units{StaticDims{K}, AffineTransform{Float64}}}

julia> typeof(uconvert(ud"°F", 1u"°C")) #Converting to dynamic units produces a dynamic quantity
Quantity{Float64, Units{Dimensions{FixRat32}, AffineTransform{Float64}}}

julia> 1u"°C" |> u"°F"  #The "pipe" operator is syntactic sugar for unit conversion
33.7999999999999 °F
```

The function `uconvert` can also be used on two unit types. This produces an `AbstractUnitTransform` that when called on a number, applies the unit conversion conversion formula produced by the two units.
```julia
julia> (u"°C" |> u"°F")(0)  #uconvert between two units produces a callable conversion formula
31.999999999999886
```

Because operations on dimensional quantities are more efficient than ones with units, FlexUnits also provides a `dconvert` function that will convert a quantity to the dimension of the units provided. This is particularly useful when converting dynamic quantities to static-dimension quantities (the most performant type).
```julia
julia> q = dconvert(u"mi/hr", 1ud"km/hr")  #dconvert converts to dimensions of the target unit
0.2777777777777778 m/s

julia> typeof(q) #Using dconvert produces high-performance quantities (especially with static units)
Quantity{Float64, StaticDims{m/s}}
```

Much like Unitful, equivalence implies conversion (especially since there is a lot of eager conversion to dimensional quantities in this package).
```julia
julia> 1u"kg" == 1000ud"g" #Equivalence implies conversion (even if using different static/dynamic modes)
true
```

FlexUnits parsing tends to be more flexible than other unit packages. For example, the empty string will parse to a unitless quantity. String macros also produce units whose symbol matches the input exactly. This information is erased however, if one tries to manipulate units directly. Because this packages converts to dimensional quantities for all mathematical operations, only dimensions need to be tracked, which is much more straightforward than tracking units.
```julia
julia> 1u""  #Empty string produces unitless value
1.0

julia> R = quantity(8.314, u"kJ/(mol K)")  #String macros produce the unit's symbols exactly
8.314 kJ/(mol K)

julia> R = quantity(8.314, u"kJ/(K*mol)")  #This displays a different result
8.314 kJ/(K*mol)

julia> u"kJ"/(u"K"*u"mol")  #Math operations on units delete symbol information (don't do this)
Units{Dimensions{FixRat32}, AffineTransform{Float64}}((m² kg)/(s² K mol), AffineTransform{Float64}(1000.0, 0.0), :_)
```

## Unit Simplification
Unitful.jl tracks units when performing operations, applying only very simple unit cancellations. While not ideal, this usually works well enough, but for long expressions, the results can be painful to look at (which is why Unitful has `upreferred`).
```julia
using Unitful

julia> r = 5u"V"/2u"A" #Ohms is better but this will do
2.5 V A^-1

julia> p = (1.0u"kg/m^3" * 9.81u"m/s^2" * 25u"cm") #This is ghastly
245.25 kg cm m^-2 s^-2
```

FlexUnits foregoes attempts to track units entirely by converting everything to dimensional quantities. This is because these operations can incur significant performance penalties when dimensions are unknown at runtime. While advantageous over long expressions, it can be difficult to interpret results, especially if they're in electrical units, which can have opaque dimensional expressions.
```julia
using FlexUnits, .UnitRegistry

julia> r = 5u"V"/2u"A" #Now this is ghastly
2.5 (m² kg)/(s³ A²)

julia> p = (1.0u"kg/m^3" * 9.81u"m/s^2" * 25u"cm") #This is a bit better
2.4525 kg/(m s²)
```

Because of this, FlexUnits includes a `simplify` function that applies a greedy algorithm over a set of common units to express the dimensions with fewer symbols.
```julia
julia> r = 5u"V"/2u"A" |> simplify #Exactly what I wanted
2.5 Ω

julia> p = (1.0u"kg/m^3" * 9.81u"m/s^2" * 25u"cm") |> simplify #And so is this
2.4525 Pa
```

It also converts logarithmic quantities to decibels
```julia
import FlexUnits: dB
julia> log_power = 20dB(u"W")
log(100.00000000000004 (m² kg)/s³)

julia> log_power = 20dB(u"W") |> simplify
20.0 dB(W)
```

Typing out "simplify" all the time can be annoying, so FlexUnits provides a `display_simplified_units` function to set the setting.
```julia
display_simplified_units(true)

julia> r = 5u"V"/2u"A"
2.5 Ω

julia> p = (1.0u"kg/m^3" * 9.81u"m/s^2" * 25u"cm")
2.4525 Pa
```

You can also change the units that the algorithm tries to simplify things into with `set_preferred_unit`. For example, to resolve pressure in PSI, simply use
```julia
set_preferred_unit(u"psi") 
julia> p = (1.0u"kg/m^3" * 9.81u"m/s^2" * 25u"cm")
0.000355704912136173 psi
```

This might beg the question, "***Why isn't simplified unit desplay turned on by default?***". There are a few reasons:

#### Reason 1: It lies. 
It only displays results _as if_ they were simplified. It does mpt simplify them. This can mess up `ustrip` if you are using units that don't convert directly to dimensional units.
```julia
julia> p = (1.0u"kg/m^3" * 9.81u"m/s^2" * 25u"cm")
0.000355704912136173 psi

julia> ustrip(p) #Results are not im PSI
2.4525
```

#### Reason 2: It is somewhat expensive.
While snappy for exploration, if you're displaying a lot of intermediate results, this may impact performance.

#### Reason 3: It may not work out-of-the box for all cases.
The preferred units are all of a certain type:
```
Units{Dimensions{FixRat32}, AffineTransform{Float64}}
```
which is the same as the default unit registry. If you create a new registry that has different dimensions (say SI dimensiosn but with angles too), you will need to point `simplify` to another set of units (sorted in order of descending complexity).

## Promotion Rules
Promotion rules will convert static dimensional quantities to dynamic ones if the dimensions are different
```julia
julia> [1u"m/s", 2u"m/s", 3u"km/hr"]  #Same static dimensions produces an array with static units
3-element Vector{Quantity{Float64, StaticDims{m/s}}}:
 1.0 m/s
 2.0 m/s
 0.8333333333333334 m/s

julia> [1u"m/s", 2u"m/s", 3u"lb/s"]   #Different static dimensions promote to a dynamic dimension
3-element Vector{Quantity{Float64, Dimensions{FixRat32}}}:
 1.0 m/s
 2.0 m/s
 1.360776 kg/s

julia> [1ud"m/s", 2ud"m/s", 3ud"lb/s"]  #Mathematical operations tend to convert to SI units for performance
3-element Vector{Quantity{Float64, Dimensions{FixRat32}}}:
 1.0 m/s
 2.0 m/s
 1.360776 kg/s
```

## Registering units
You can register units to UnitRegistry (or your own registry) using the `register_unit` function. You can simply provide it a single pair argument of type `Pair{String, AbstractUnitLike}` or `Pair{String, Quantity}` as follows:

```julia
julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:
  :Ω      => Ω
  :μs     => μs
  :μV     => μV
```

However, due to the nature of macros, these dictionaries are permanent. You can re-register units with the same values (so that you can re-run scripts) but changing them is not allowed.
```julia
julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:
  :Ω      => Ω
  :μs     => μs
  :μV     => μV

julia> register_unit("bbl" => 22.5*u"m^3")
ERROR: PermanentDictError: Key bbl already exists. Cannot assign a different value.
```

It is possible for users to create their own registries (for example, if they wanted different dimension types or different classes of units). The default UnitRegistry (see the UnitRegistry.jl file in the source code) was constructed with less than 50 lines of code, and users can use that as a template. Most custom registries only need to modify two of these lines of code:
```julia
using .RegistryTools
const UNITS = PermanentDict{Symbol, Units{Dimensions{FixRat32}, AffineTransform}}()

#Fill the UNITS registry with default values
registry_defaults!(UNITS)
```
The `UNITS` constant is the dictionary where all the units live. One can customize it so that the registry has a different dimension type or even a different transform type. The `registry_defaults!` function fills a dictionary with the default body of units (users can build their own if they like).

## Linear Algebra
This package defines a `LinmapQuant`, a matrix that is intended to be a linear mapping that takes a vector with units `u_in` and produces a vector with units `u_out` thus, the dimensions of all elements can be inferred by these two vectors of units. While in general, the units of matrix elements can be arbitrary, in order to support operations like matrix multiplication, the units must adhere to this structure, and the simplicity of this structure allows for shortcuts for inference.

If we want to construct a matrix where all of the columns have the same units, we let `u_in` be the inverse of the units we desire, and `u_out` be dimensionless.
```julia
julia> Z = randn(200, 5)*randn(5,5)
julia> Zu = LinmapQuant(Z, UnitMap(u_in=inv.([u"lb", u"ft", u"W", u"L", u"mol"]), u_out=u""))
200×5 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  -0.501901 kg    0.137162 m  -0.689416 (m² kg)/s³   -0.00173674 m³    4.09169 mol
    1.84433 kg    0.335457 m   0.987017 (m² kg)/s³   -0.00204628 m³    3.95951 mol
  -0.470226 kg   0.0676007 m   -2.44338 (m² kg)/s³     0.0013448 m³   -2.37904 mol
  -0.373749 kg    0.157584 m   -2.97929 (m² kg)/s³    0.00555965 m³   -2.68827 mol
```
We could have also built the quantity matrix first and then used `LinmapQuant(m::AbstractMatrix{<:Quantity})`. Now let us suppose we wanted to do linear regression to predict the last two columns given the first three.
```julia
julia> Xu = [Zu[:,1:3] fill(1.0u"", size(Zu,1))]
200×4 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  -0.501901 kg    0.137162 m  -0.689416 (m² kg)/s³  1.0
    1.84433 kg    0.335457 m   0.987017 (m² kg)/s³  1.0
  -0.470226 kg   0.0676007 m   -2.44338 (m² kg)/s³  1.0
  -0.373749 kg    0.157584 m   -2.97929 (m² kg)/s³  1.0

julia> Yu = Zu[:,4:5]
200×2 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  -0.00173674 m³    4.09169 mol
  -0.00204628 m³    3.95951 mol
    0.0013448 m³   -2.37904 mol
   0.00555965 m³   -2.68827 mol

julia> Bu = (Xu'*Xu)\(Xu'*Yu)
4×2 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
     -0.000206539 m³/kg          -0.109595 mol/kg
          -0.0032453 m²             2.41413 mol/m
 -0.000997411 (m s³)/kg  0.75641 (s³ mol)/(m² kg)
        -0.000152683 m³              0.159627 mol
```
We can verify that the output of the prediction matches `Yu`
```julia
julia> Xu*Bu
200×2 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  0.000193479 m³  0.0242788 mol
  -0.00260673 m³    1.51392 mol
    0.0021621 m³   -1.47384 mol
   0.00238468 m³   -1.67255 mol
    -0.001745 m³   0.810835 mol
```
