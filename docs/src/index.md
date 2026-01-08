```@meta
DocTestSetup = quote
    using FlexUnits, .UnitRegistry
end
```
# FlexUnits.jl
FlexUnits.jl is a unit package designed to resemble Unitful.jl, with similar performance when units can be statically inferred, but leverages techniques in DynamicQuantities.jl to eliminate many of Unitful's performance pitfalls when units are uninferrable. This package allows the user to employ the following useful pattern:
1. Define high-level operations on quantities with dynamic units that are type-stable even with mixed unit types
2. Convert quantities to statically-inferred high-performance units in lower-level performance-sensitive code
3. Return values in dynamic units with confidence that function output is type-stable

## Quick start examples
The FlexUnits API is designed to resemble Unitful in a number of ways. One major difference is that string macros and parsing functions are not exported by default, instead they are exported by a unit registry. This way, users can define their own registries and use those as a basis for unit parsing. FlexUnits provides a unit registry (UnitRegistry) that can be imported to provide default unit parsing functionality.
- This package requires Julia 1.10. Older versions will not be supported.
- `] add FlexUnits`
- `using FlexUnits, .UnitRegistry`

```@jldoctest
julia> 1u"°C"   #@u_str produces static units which automatically convert quantities to SI on multiplication
274.15 K

julia> 1ud"°C"   #@ud_str produces dynamic units which don't eagerly convert
1 °C

julia> uparse("km/hr")  #Uparse always produces dynamic units
km/hr

julia> typeof(uparse("km/hr")) == typeof(uparse("kPa")) #Uparse is type-stable
true

julia> typeof(uconvert(u"°F", 1u"°C")) #Converting to static units produces a static quantity
Quantity{Float64, StaticUnits{K, AffineTransform}}

julia> typeof(uconvert(ud"°F", 1u"°C")) #Converting to dynamic units produces a dynamic quantity
Quantity{Float64, Units{Dimensions{FixRat32}, AffineTransform}}

julia> 1u"°C" |> u"°F"  #The "pipe" operator is syntctic sugar for unit conversion
33.799999999999955 °F

julia> (u"°C" |> u"°F")(0)  #uconvert between two units produces a callable conversion formula
31.999999999999943

julia> q = dconvert(u"mi/hr", 1ud"km/hr")  #dconvert converts to dimensions of the target unit
0.2777777777777778 m/s

julia> typeof(q) #Using dconvert produces high-performance quantities (expecially with static units)
Quantity{Float64, StaticDims{m/s}}

julia> 1u"kg" == 1000ud"g" #Equivalence implies conversion (even if using different static/dynamic modes)
true

julia> 1u""  #Empty string produces unitless value
1.0

julia> R = 8.314ud"kJ/(mol K)"  #String macros produce the unit's symbols exactly
8.314 kJ/(mol K)
julia> R = 8.314ud"kJ/(K*mol)"  #This displays a different result
8.314 kJ/(K*mol)
julia> u"kJ"/(u"K"*u"mol")  #Math operations on units delete symbol information (don't do this)
Units{Dimensions{FixRat32}, AffineTransform}((m² kg)/(s² K mol), AffineTransform(1000.0, 0.0), :_)

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