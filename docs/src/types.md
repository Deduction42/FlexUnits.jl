
# Types

## Quantities
FlexUnits.jl has two types of quantity values. FlexQuant, which can contain any object as a value, and Quantity which subtypes to Number (following the convention of most of Julia's unit packages). For convenience, they are often referred to together as `QuantUnion`. Mathematical operations on a `QuantUnion` with a Unit will convert to a dimensional quantity (SI by default) having an `AbstractDimension` instead of a `AbstractUnit`.
```@docs
Quantity
FlexQuant
QuantUnion
```

## Dimensions
Dimensions are the core object used to define units, and `Quantity` values tend to have some form of `AbstractDimension` in its unit field after calculations due to computational simplicity. 
```@docs
AbstractDimensions
Dimensions
StaticDims
```

## Units
Units contain an `AbstractDimensions` and an `AbstractUnitTransform` which represents a conversion formula to convert the Quantity to its pure dimensional form. Typically, the unit transform is called before a calcualtion to convert the quantity to dimensional form, and `uconvert` applies an invers of that transform to convert the result back to the desired unit. Because dimensions have no transforms, it is most efficient to perform calculations directly on dimensional quantities.
```@docs
AbstractUnits
Units
StaticUnits
```

## Unit Transforms
One unique feature to the FlexUnits design is the `AbstractUnitTransform` object. This is a callable object that contains a conversion formula to convert the unit into its dimensional form. The default transforms are `NoTransform` (a property of all dimensions) and `AffineTransform` which can deal with all common linear units. This design could potentially support other transforms like `LogTransform` for logarithmic units likd `dB` and `pH` in the future.
```@docs
AbstractUnitTransform
NoTransform
AffineTransform
```

## Linear Algebra
Linear algebra functionality is achieved by observing that any matrix of quantities that supports multiplication is a special kind of quantity matrix that functions as a *linear mapping*. Such matrices have a special structure for units. These special unit structures are *unit mappings* that contain an input vector of units, and an output vector of units. Accellerated linear algebra operations are achived by keeping numerical matrices and unit mappings separate in a `LinmapQuant` and performing the linear algebra and unit inference separately (not unlike how a `Quantity` operates).

### Unit Maps
```@docs
UnitMap
DimsMap
```

### Objects with Unit Maps
```@docs
LinmapQuant
FunctionQuant
FactorQuant
```
