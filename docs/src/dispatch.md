# Dispatch Patterns

The most defining differences between FlexUnits and Unitful lie in their respective type systems and dispatch strategy:
1. FlexUnits quantities do not need any information about dimensions and units to be concrete
2. FlexUnits can specialize at the level of dimensions and still be concrete
3. FlexUnits cannot specialize beyond the level of dimensions (does not specialize to units and this is generally unnecessary)

## Parameterization Example: Sensor Model
Let's consider a simple example of a sensor model with a value, and a tolerance. The job of this model is to assign units as context to numerical values coming from a raw dataset. The first thing to note about FlexUnits is that Quantity can have units or dimensions
1. `Quantity{T,U} where {T<:Number, U<:AbstractUnitLike}` represents a quantity with any units
2. `Quantity{T,D} where {T<:Number, D<:AbstractDimLike}` represents a quantity with base units (preferred)

Dimensions are a privileged unit type where mathematical operations between dimensions can occur *without any need for scaling and offsetting*. Recall that a Unit contains three pieces of information:
1. Dimensions (can be static or dynamic)
2. A transformation formula to convert to base units (usually AffineTransform{Float64}, but logarithmic units have ExpAffTransform{Float64})
3. A symbol (often directly obtained from a string macro, or a dimensional display if converting from dimensions)

If the quantity's unit is actually raw dimensions, the transformation is a no-op, and the symbol is inferred directly from the dimensional value itself. This makes `Quantity{T,<:AbstractDimLike}` simpler and more efficient than `Quantity{T,<:AbstractUnitLike}` which is the reason why numerical operations will produce quantities with dimensions rather than units.

```julia
abstract type AbstractSensor{T, D<:AbstractDimLike} end

#Correct, but complicated and slightly slower over repeated calls
mutable struct UnitSensor{T<:Real, D<:AbstractDimLike} <: AbstractSensor{T,D}
    value :: Quantity{T, Units{D, AffineTransform{Float64}}}
    tolerance :: Quantity{T, Units{D, AffineTransform{Float64}}}
    units :: Units{D, AffineTransform{Float64}}
end
UnitSensor{T,D}(u::Units) where {T,D} = UnitSensor{T,D}(NaN*u, NaN*u, u)

#Correct, simpler and more performant (dimensional quantities don't require any conversions)
mutable struct DimSensor{T<:Real, D<:AbstractDimLike} <: AbstractSensor{T,D}
    value :: Quantity{T, D}
    tolerance :: Quantity{T, D}
    units :: Units{D, AffineTransform{Float64}}
end
DimSensor{T,D}(u::Units) where {T,D} = DimSensor{T,D}(NaN*u, NaN*u, u)
```

## Dispatch Example: Reading and Writing Values
Conversions between dimensional and unit quantities should be seamless, so it is possible to write generic functions regardless of these unit types. While you can construct quantities through multiplication, using `quantity` is more explicit and lazier as multiplication will convert to base units. In the example below, if a raw numeric value is written, it's interpreted as having the sensor's units. Otherwise, if it's a quantity, that quantity is written directly.

```julia
writeval!(s::AbstractSensor, v::Quantity) = setproperty!(s, :value, v)
writetol!(s::AbstractSensor, v::Quantity) = setproperty!(s, :tolerance, v)
writeval!(s::AbstractSensor, v::Real) = writeval!(s, quantity(v, s.units))
writetol!(s::AbstractSensor, v::Real) = writetol!(s, quantity(v, s.units))
readval(s::AbstractSensor) = s.value
readrawval(s::AbstractSensor) = s.value |> s.units
```

## Dispatch Example: Generic Sensors
If we don't know what kind of sensor we are dealing with (only that the values are Float64), we can create a generic parameterization of the sensor types that can hold any unit.
```julia
const GenericUnitSensor = UnitSensor{Float64, Dimensions{FixRat32}}
const GenericDimSensor = DimSensor{Float64, Dimensions{FixRat32}}
```
This is advantageous if you need to have multiple sensor types in a collection (such as a dictionary or a vector). Any kind of sensor type can be *concretely* represented by this generic type.
```julia
julia> isconcretetype(GenericUnitSensor)
true

julia> isconcretetype(GenericDimSensor)
true
```

The only problem is that there is no unit validation, so it is possible to write values with conflicting units in this case:
```julia
#These units will be conflicting
s = GenericUnitSensor(u"kg")
writeval!(s, 10u"km")
writetol!(s, quantity(5, u"°C"))

julia> s
GenericUnitSensor(10000.0 m, 5.0 °C, kg)

#These dimensions will be conflicting
s = GenericDimSensor(u"kg")
writeval!(s, 10u"km")
writetol!(s, quantity(5, u"°C"))

julia> s
GenericDimSensor(10000.0 m, 278.15 K, kg)
```

This means that a lot of your functions will need dimensional validation
```julia
writeval!(s::UnitSensor{<:Any, <:Dimensions}, v::Quantity) = setproperty!(s, :value, v |> s.units)
writetol!(s::UnitSensor{<:Any, <:Dimensions}, v::Quantity) = setproperty!(s, :tolerance, v |> s.units)
writeval!(s::DimSensor{<:Any, <:Dimensions}, v::Quantity) = setproperty!(s, :value, v |> dimension(s.units))
writetol!(s::DimSensor{<:Any, <:Dimensions}, v::Quantity) = setproperty!(s, :tolerance, v |> dimension(s.units))

#Writing a value in the wrong dimensions produces an error
s = GenericUnitSensor(u"g")

julia> writeval!(s, 10u"km")
ERROR: ConversionError: Cannot convert unit 'm' to target unit 'g' due to a dimension mismatch of 'm/kg'

s = GenericDimSensor(u"g")

julia> writeval!(s, 10u"km")
ERROR: ConversionError: Cannot convert unit 'm' to target unit 'kg' due to a dimension mismatch of 'm/kg'
```

## Dispatch Example: Pressure Sensors
While there is appeal of having generic sensors that can concretely represent any dimension, there are a couple of notable downsides
1. This `Dimensions{FixRat32}` dimension type, while concrete, is still slower than static alternatives
2. Objects with this type will need to do manual dimensional validation to avoid nonsense values

FlexUnits allows you to have static dimensional types that force consistency and improve speed. This type of static specialization is facilitated by the `@D_str` macro.

```julia
#We focus on DimSesnor here because it is simpler and more performant
const PressureSensor = DimSensor{Float64, D"psi"} #Unit inside D is irrelevant, static dimensions are returned
PressureSensor (alias for DimSensor{Float64, StaticDims{kg/(m s²)}})

julia> ps = PressureSensor(u"K") #You can't even construct an inconsistent object
ERROR: ConversionError: Cannot convert unit 'K' to target unit 'kg/(m s²)' due to a dimension mismatch of '(m s² K)/kg'

julia> ps = PressureSensor(5u"kg", 3u"Pa", u"Pa") #Another failure on creating an inconsistent object
ERROR: ConversionError: Cannot convert unit 'kg' to target unit 'kg/(m s²)' due to a dimension mismatch of '(m s²)'

julia> ps = PressureSensor(u"kPa");

julia> writeval!(ps, 10u"psi") #This works
68947.6 kg/(m s²)

julia> writeval!(ps, 5) #This also works (remember units are in kPa)
5000.0 kg/(m s²)

julia> writeval!(ps, 5u"m") #This throws an error without need of validation code
ERROR: ConversionError: Cannot convert unit 'm' to target unit 'kg/(m s²)' due to a dimension mismatch of '(m² s²)/kg'

julia> readval(ps) #This produces high-performance quantities without any effort
5000.0 kg/(m s²)

julia> readrawval(ps) #This produces values in sensor's units (useful for displaying raw values for auditing)
5.0 kPa
```
