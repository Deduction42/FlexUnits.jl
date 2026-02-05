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
```
add FlexUnits
using FlexUnits, .UnitRegistry
```

Units are accessed primarily through string macros. The `@u_str` macro produces static units by default (for maximum performance when injected into low-level code), but if dynamic units are desired (particularly for interactive workflows) one can use the `@ud_str` macro. 
```
julia> 1u"°C"   #@u_str produces static units (multiplication converts to base units)
274.15 K

julia> 1ud"°C"  #@ud_str produces dynamic units (multiplication converts to base units)
274.15 K

julia> quantity(1, u"°C") #Use "quantity" to avoid eager conversion to base units
1 °C
```

Unlike macros, units that need to be parsed from strings cannot be statically inferred, so the `uparse` function will always produce dynamic units in order to make it type-stable.
```
julia> uparse("km/hr")  #Uparse always produces dynamic units
km/hr

julia> typeof(uparse("km/hr")) == typeof(uparse("kPa")) #Uparse is type-stable
true
```

The function `uconvert` pays attention to the unit type converted to. If the unit is static, the resulting quantity will be static, if it's dynamic, the result will be dynamic.
```
julia> typeof(uconvert(u"°F", 1u"°C")) #Converting to static units produces a static quantity
Quantity{Float64, StaticUnits{K, AffineTransform}}

julia> typeof(uconvert(ud"°F", 1u"°C")) #Converting to dynamic units produces a dynamic quantity
Quantity{Float64, Units{Dimensions{FixRat32}, AffineTransform}}

julia> 1u"°C" |> u"°F"  #The "pipe" operator is syntctic sugar for unit conversion
33.799999999999955 °F
```

The function `uconvert` can also be used on two unit types. This produces an `AbstractUnitTransform` that when called on a number, applies the unit conversion conversion formula produced by the two units.
```
julia> (u"°C" |> u"°F")(0)  #uconvert between two units produces a callable conversion formula
31.999999999999943
```

Because operatins on dimensional quantities are more efficient than ones with units, FlexUnits also provides a `dconvert` function that will convert a quantity to the dimension of the units provided. This is particularly useful when converting dynanmic quantities to static-dimension quantities (the most peformant type).
```
julia> q = dconvert(u"mi/hr", 1ud"km/hr")  #dconvert converts to dimensions of the target unit
0.2777777777777778 m/s

julia> typeof(q) #Using dconvert produces high-performance quantities (expecially with static units)
Quantity{Float64, StaticDims{m/s}}
```

Much like Unitful, equivalence implies conversion (especially since there is a lot of eager conversion to dimensional quantities in this package).
```
julia> 1u"kg" == 1000ud"g" #Equivalence implies conversion (even if using different static/dynamic modes)
true
```

FlexUnits parsing tends to be more flexible than other unit packages. For example, the empty string will parse to a unitless quantity. String macros also produce units whose symbol matches the input exactly. This information is erased however, if one tries to manipulate units directly. Because this packages converts to dimensional quantities for all mathematical operations, only dimensions need to be tracked, which is much more straightforward than tracking units.
```
julia> 1u""  #Empty string produces unitless value
1.0

julia> R = 8.314ud"kJ/(mol K)"  #String macros produce the unit's symbols exactly
8.314 kJ/(mol K)

julia> R = 8.314ud"kJ/(K*mol)"  #This displays a different result
8.314 kJ/(K*mol)

julia> u"kJ"/(u"K"*u"mol")  #Math operations on units delete symbol information (don't do this)
Units{Dimensions{FixRat32}, AffineTransform}((m² kg)/(s² K mol), AffineTransform(1000.0, 0.0), :_)
```

Promotion rules will convert static dimensional quantities to dynamic ones if the dimensions are different
```
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

julia> [1ud"m/s", 2ud"m/s", 3ud"lb/s"]  #Dynamic units are untouched
3-element Vector{Quantity{Int64, Units{Dimensions{FixRat32}, AffineTransform}}}:
 1 m/s
 2 m/s
 3 lb/s

julia> [1ud"m/s", 2ud"m/s", 3ud"lb/s"].*1  #Mathematical operations tend to convert to SI units for performance
3-element Vector{Quantity{Float64, Dimensions{FixRat32}}}:
 1.0 m/s
 2.0 m/s
 1.360776 kg/s
```

## Registering units
You can register units to UnitRegistry (or your own registry) using the `register_unit` function. You can simply provide it a single pair argument of type `Pair{String, AbstractUnitLike}` or `Pair{String, Quantity}` as follows: 

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

It is possible for users to create their own registries (for example, if they wanted different dimension types or different classes of units). The default UnitRegistry (see the UnitRegistry.jl file in the source code) was constructed with less than 50 lines of code, and users can use that as a template. Most custom registries only need to modify two of these lines of code:
```
const UNITS = PermanentDict{Symbol, Units{Dimensions{FixRat32}, AffineTransform}}()

#Fill the UNITS registry with default values
registry_defaults!(UNITS)
```
The `UNITS` constant is the dictionary where all the units live. One can customize it so that the registry has a different dimension type or even a different transform type. The `registry_defaults!` function fills a dictionary with the default body of units (users can build their own if they like).


## Linear Algebra
This package devines a `LinmapQuant`, a matrix that is intended to be a linear mapping that takes a vector with units `u_in` and produces a vector with units `u_out` thus, the dimensiosn of all elements can be inferred by these two vectors of units. While in general, the units of matrix elements can be arbitrary, in order to support operations like matrix multiplication, the units must adhere to this structure, and the simplicity of this structure allows for shortcuts for inference.

If we want to construct a matrix where all of the columns have the same units, we let `u_in` be the inverse of the units we desire, and `u_out` be dimensionless.
```
julia> Z = randn(200, 5)*randn(5,5)
julia> Zu = LinmapQuant(Z, UnitMap(u_in=inv.([u"lb", u"ft", u"W", u"L", u"mol"]), u_out=u""))
200×5 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  -0.501901 kg    0.137162 m  -0.689416 (m² kg)/s³   -0.00173674 m³    4.09169 mol
    1.84433 kg    0.335457 m   0.987017 (m² kg)/s³   -0.00204628 m³    3.95951 mol
  -0.470226 kg   0.0676007 m   -2.44338 (m² kg)/s³     0.0013448 m³   -2.37904 mol
  -0.373749 kg    0.157584 m   -2.97929 (m² kg)/s³    0.00555965 m³   -2.68827 mol
```
We could have also built the quantity matrix first and then used `LinmapQuant(m::AbstractMatrix{<:Quantity})`. Now let us suppose we wanted to do linear regression to predict the last two columsn given the first three.
```
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
```
julia> Xu*Bu
200×2 LinmapQuant{Float64, Dimensions{FixRat32}, Matrix{Float64}, DimsMap{Dimensions{FixRat32}, Vector{Dimensions{FixRat32}}, Vector{Dimensions{FixRat32}}}}:
  0.000193479 m³  0.0242788 mol
  -0.00260673 m³    1.51392 mol
    0.0021621 m³   -1.47384 mol
   0.00238468 m³   -1.67255 mol
    -0.001745 m³   0.810835 mol
```