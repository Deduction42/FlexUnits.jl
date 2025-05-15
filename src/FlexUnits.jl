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
export ConversionError, DimensionError, NotScalarError, NotDimensionError
export Dimensions, NoDims, AffineUnits, Quantity, UnitfulCallable, AbstractQuantity, AffineTransform
export RegistryTools, UnitRegistry
export static_fieldnames, uscale, uoffset, dimension, pretty_print_units
export assert_scalar, assert_dimension, assert_dimensionless
export dimensionalize, quantity, ustrip, unit, unitful_call
export ubase, uconvert, ustrip_base, ustrip_dimensionless, transform


end
