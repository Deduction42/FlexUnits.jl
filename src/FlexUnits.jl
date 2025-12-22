module FlexUnits

# Write your package code here.
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")

export AbstractUnitLike, AbstractDimensions, AbstractUnits, AbstractAffineUnits, AbstractUnitTransform
export ConversionError, DimensionError, NotScalarError, NotDimensionError, FixRat32
export Dimensions, AffineUnits, Quantity, AbstractQuantity, AffineTransform, MirrorDims, MirrorUnion
export RegistryTools, UnitRegistry
export static_fieldnames, uscale, uoffset, dimension, pretty_print_units
export assert_scalar, assert_dimension, assert_dimensionless
export with_ubase, ustrip, unit
export ubase, uconvert, ustrip_base, ustrip_dimensionless


end
