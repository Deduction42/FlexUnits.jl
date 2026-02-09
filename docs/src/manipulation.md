# Unit Manipulation and Conversion

## Retrieving units from a registry
Unit registries are modules that export parsing functions and string macros. FlexUnits contains a default unit registry module named `UnitRegistry` that contains scalar/affine units. To use string macros, simply call the string macros and parsing functions from the registry.
```julia
using FlexUnits
distance1 = 1*UnitRegistry.u"km"
distance2 = 1*UnitRegistry.uparse("m")
distance3 = 1*UnitRegistry.qparse("5 cm")

julia> distance1 + distance2 + distance3
1001.05 m
```

You can export these functions and macros directly from the module (this is not done by default so that users can use their own custom modules if they like).
```julia
using FlexUnits, .UnitRegistry
distance1 = 1*u"km"
distance2 = 1*uparse("m")
distance3 = 1*qparse("5 cm")

julia> distance1 + distance2 + distance3
1001.05 m
```

### String macro and parsing behaviour
There are two types of string macros:
1. One that produces staticly-typed units `@u_str` such as `1u"m/s"`
2. One that produces dynamic units `@ud_str` such as `1ud"m/s"`

Most of the time users will want to use the static unit macro inside code, because the Julia compiler can reason about the units with constants, resulting in better runtime because unit checking is all done at compile-time. One can check the types produced by these macros.

```julia
julia> typeof(1u"m/s")
Quantity{Float64, StaticDims{m/s}}

julia> typeof(1ud"m/s")
Quantity{Float64, Dimensions{FixRat32}}
```

There are also two types of parsing functions:
1. One that parses units `uparse` such as `uparse("N")`
2. One that parses quantities `qparse` such as `qparse("5 lbf")`

Parsing functions always produce dynamic units, as they are primarily used to parse free text as units (such as csv spreadsheets); in such cases, producing dynamic units resolves type-stability issues from not knowing the units beforehand. 

It should also be notced that using string macros and parsing functions will cause units to be displayed exactly how they will be parsed (use `quantity` to prevent eagerly converting to dimensional units).
```julia
julia> quantity(5, u"J/(K*mol)")
5 J/(K*mol)

julia> quantity(5, u"J/(mol K)")
5 J/(mol K)

julia> quantity(5, u"J*(mol K)^(-1)")
5 J*(mol K)^(-1)
```

Dimensional units, however, are displayed in a standard format because their symbols are well-defined.
```julia
julia> 5u"J/(K*mol)"
5.0 (m² kg)/(s² K mol)

julia> 5*u"J*(mol K)^(-1)"
5.0 (m² kg)/(s² K mol)
```

### Registering units
You can register units to a registry using other units or quantities using the `register_unit` function:
```
julia> register_unit("bbl" => 0.158987*u"m^3")
```
However, because macros only look up units at compile time, *changing* these values won't update macro outputs in functions that have already been compiled. Becuase of this, we use the `PermanentDict` object to produce errors when overwriting units with different values. You can re-register units with the same values (so that you can re-run scripts) but overwriting them as different values is not allowed.
```
julia> register_unit("bbl" => 0.158987*u"m^3")
FlexUnits.RegistryTools.PermanentDict{Symbol, AffineUnits{Dimensions{FixedRational{Int32, 25200}}}} with 150 entries:

julia> register_unit("bbl" => 22.5*u"m^3")
ERROR: PermanentDictError: Key bbl already exists. Cannot assign a different value.
```

## Converting units
FlexUnits mimics Unitful where possible, so the `uconvert` function is used to convert quantities to desired units. However, internally, most math is done using dimensional quantities (i.e. SI units) so that no internal conversion factors are neccessary. Because dimensional quantities are higher-performacne, the `dconvert` function is also used, mainly to verify the dimensions of the input units and remove any scaling/offsets before calculation.

### Use `uconvert` for displaying results
The function `uconvert` is most often used for converting results of a calculation (usually SI) to the desired units.
```@docs
uconvert(::AbstractUnitLike, ::QuantUnion)
```
Note that using `uconvert` between two unit object will produce a unit conversion formula that can be called directly. Additionally, pipe operator `|>` can also be used as a shorthand for `uconvert`
```@docs
uconvert(::AbstractUnitLike, ::AbstractUnitLike)
|>(u1::AbstractUnitLike, u2::Union{AbstractUnitLike, QuantUnion})
```

### Use `dconvert` before performance-sensitive code
The funciton `dconvert` converts a quantity to the *dimensions* of the desired unit. This is a convenience function since string macros produce units, not dimensions. This may produce unintuitive results if you forget this behaviour, for example,
```julia
julia> dconvert(u"km/hr", 25u"km/hr")
6.944444444444445 m/s
```
This is useful for validating input unit dimensions and converting to a high-performance object all in one step before performing performance-sensitive code. In such cases, make sure you use the static version (i.e. use `5u"m/s"` not `5ud"m/s"`).
```@docs
dconvert(::AbstractUnitLike, ::FlexUnits.QuantUnion)
```

## Extracting units and values from quantities
The following functions can be used to extract/inspect different attributes of a quantity
- `ustrip(q::Quantity)` extracts the raw numerical value from a quantity
- `ustrip(u::AbstractUnitLike, q::Quantity)` returns the value of a quantity in the desired units
- `dstrip(q::Quantity)` converts a quantitty to dimensional units and returns the numerical value
- `unit(q::Quantity)` extracts the unit from a quantity (`Quantity{Float64, <:AbstractDimLike}` will return a dimension)
- `dimension(q::Quantity)` will return the dimensions of a quantity (if dimensions are static, it returns the static value)

