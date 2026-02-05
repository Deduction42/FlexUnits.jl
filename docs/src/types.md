
# Types

## Quantities
FlexUnits.jl only has one quantity type Quantity containing a value and a unit-like object (units or dimensions). Unlike Unitful, it doesn't subtype into Number, so its value can take any type. Mathematical operations on a Quantity with a Unit will convert to a dimensional quantity (SI by default) having an `AbstractDimension` instead of a `AbstractUnit`.
```@docs
Quantity
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