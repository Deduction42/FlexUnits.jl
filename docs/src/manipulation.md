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

## Extracting units and values from quantities
The following functions can be used to extract/inspect different attributes of a quantity
- `ustrip(q::Quantity)` extracts the raw numerical value from a quantity
- `ustrip(u::AbstractUnitLike, q::Quantity)` returns the value of a quantity in the desired units
- `dstrip(q::Quantity)` converts a quantity to dimensional units and returns the numerical value
- `unit(q::Quantity)` extracts the unit from a quantity (`Quantity{Float64, <:AbstractDimLike}` will return a dimension)
- `dimension(q::Quantity)` will return the dimensions of a quantity (if dimensions are static, it returns the static value)

## Converting units
FlexUnits mimics Unitful where possible, so the `uconvert` function is used to convert quantities to desired units. However, internally, most math is done using dimensional quantities (i.e. SI units) so that no internal conversion factors are necessary. Because dimensional quantities are higher-performance, the `dconvert` function is also used, mainly to verify the dimensions of the input units and remove any scaling/offsets before calculation.

### Use `uconvert` for displaying results
The function `uconvert` is most often used for converting results of a calculation (usually SI) to the desired units.
```@docs
uconvert(::AbstractUnitLike, ::QuantUnion)
```
Note that using `uconvert` between two unit object will produce a unit conversion formula that can be called directly. 
```@docs
uconvert(::AbstractUnitLike, ::AbstractUnitLike)
```
Additionally, pipe operator `|>` can also be used as a shorthand for `uconvert`
```julia
5u"m/s" |> u"km/hr"
18.0 km/hr
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

### String macro and parsing behaviour
There are two types of string macros for units:
1. One that produces statically-typed units `@u_str` such as `1u"m/s"`
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

It should also be noticed that using string macros and parsing functions will cause units to be displayed exactly how they will be parsed (use `quantity` to prevent eagerly converting to dimensional units).
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

### String macros for static unit types
This package also exports type macros allowing users to easily constrain and dispatch on static dimension types. The two macros are
1. `@D_str` which produces the static dimension type of the string expression
2. `@U_str` which produces the static unit type of the string expression 
Generally the use of `@D_str` is preferred as this will produce more performant objects. Lets consider an example where these macros could be used. Consider a process where flow and pressure are measured. One could create an object that stipulates these measurement types.
```julia 
struct MyProcess
    flow :: Quantity{Float64, U"kg/hr"}
    pressure :: Quantity{Float64, U"kPa"}
end

julia> MyProcess(quantity(5, u"lb/hr"), quantity(6, u"kPa"))
MyProcess(5.0 lb/hr, 6.0 kPa)
```

Using incompatible dimensions will desirably produce an error
```julia
julia> MyProcess(quantity(5, u"ft^3/hr"), quantity(6, u"kPa"))
ERROR: ConversionError: Cannot convert unit 'm³/s' to target unit 'kg/s'. Consider multiplying 'kg/s' by 'm³/kg' or similar.
```

Similarly, one can also use this pattern for dimensions as well, that will automatically convert fields into high-performance static dimension quantities
```julia
struct MySiProcess
    flow :: Quantity{Float64, D"kg/hr"}
    pressure :: Quantity{Float64, D"kPa"}
end

julia> MySiProcess(quantity(5, u"lb/hr"), quantity(6, u"kPa"))
MySiProcess(0.0006299888888888889 kg/s, 6000.0 kg/(m s²))
```

Another useful trick with the `@D_str` macro is using the pipe operator as a shorthand equivalent of `dconvert`, which validates dimensions and produces a high-performance quantity with static SI base units.
```julia
#The first function validates and converts, while the second function has no dynamic overhead
pressure(density::Quantity, head::Quantity) = pressure(density |> D"kg/m^3", head |> D"m")
pressure(density::Quantity{<:Any, D"kg/m^3"}, head::Quantity{<:Any, D"m"}) = (9.81u"m/s^2")*density*head 
```
While the example above won't improve performance due to the simplicity of the calculation, more complicated functions benefit from not having to dynamically validate units at all.


## Logarithmic Quantities and Units
Handling logarithmic units such as decibels comes with some level of controversy. The question is whether a quantity such as `1 dB` is a logarithm of a quantity or merely the *logarithmic representation* of a linear quantity.
1. If `1 dB` is a logarithmic representation: `1 dB + 1 dB = 4.0103 dB`  
2. If `1 dB` is a logarithmic quantity: `1 dB + 1 dB = 2 dB`

FlexUnits adopts the philosophy of the ***second camp***, where a decibel represents the *logarithm of a quantity* and supplies algebraic tools to manipulate logarithms of quantities (referred to as a `LogQuant`). 

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
Finally, converting to a logarithmic unit also produces a logarithmic quantity
```julia
julia> 100u"m/s" |> dB(u"m/s")
19.999999999999996 dB(m/s)
```

### Operations on logarithmic quantities
The algebraic rules `LogQuant` are centered around logarithmic identities. Addition will add the numeric values and multiply the units.
```julia
julia> log(4u"m") + log(4u"s") # log(x) + log(y) = log(x*y)
log(15.999999999999998 (m s))
```
Subtraction will subtract the numeric values and divide the units.
```julia
julia> log(4u"m") - log(4u"s") # log(x) - log(y) = log(x/y)
log(1.0 m/s)
```
Multiplication will multiply the numeric values and raise the units to a power.
```julia
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

### Logarithmic units
FlexUnits also contains support for logarithmic unit scales such as decibels `dB` and Nepers `Np` which are callable `LogScale` objects. There aren't exported by default as they could be fairly common symbols. To use them, simply call them as a function on units.
```julia
julia> dB(u"V")
dB(V)
```
This operation produces a logarithmic unit type `Units{<:AbstractDimLike, <:ExpAffTransform}`. Multiply a number by such a unit will produce a logarithmic quantity.
```julia
julia> 5dB(u"kW")
log(3162.2776601683804 (m² kg)/s³)
```
Calling `dB` on a logarithmic quantity (piping will also work) will convert it to decibels with SI base units as reference. 
```julia
julia> 5dB(u"kW") |> dB
34.99999999999999 dB((m² kg)/s³)
```
Converting to a quantity to a logarithmic unit will also result in a logarithmic quantity
```julia
julia> 5u"hp" |> dB(u"kW")
5.715340722972715 dB(kW)
```

## Registering Units

### Basic units
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

### Logarithmic units
The default unit registry can only register affine units, which is sufficient for logarithmic units *unless you need to parse strings to produce logarithmic units*. In such cases, you will need to register logarithmic units with the `LogUnitRegistry` instead. This registry can hold both affine and logarithmic units, but `uparse` can introduce performance issues because the output is a `Union`. ***WARNING, because multiplying `uparse` outputs can produce a Quantity or a LogQuant, based on the string value, it's recommended that you use explicit constructors like `quantity` or `ubase` to always produce linear quantities, or `logquant` or `logubase` always produce logarithmic units.***
```julia
using FlexUnits, .LogUnitRegistry
import FlexUnits.dB

register_unit("dB_V" => dB(u"V"))

julia> 10uparse("dB_V")
log(10.000000000000002 (m² kg)/(s³ A))

julia> 10uparse("V")
10.0 (m² kg)/(s³ A)

julia> ubase(10, uparse("dB_V"))
10.000000000000002 (m² kg)/(s³ A)

julia> logubase(10, uparse("V"))
log(10.000000000000002 (m² kg)/(s³ A))

julia> 10u"dB_V"
log(10.000000000000002 (m² kg)/(s³ A))
```




