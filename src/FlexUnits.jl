module FlexUnits

# Write your package code here.
include("fixed_rational.jl")
include("types.jl")
include("linalg_types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("linear_algebra.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")

export AbstractUnitLike, AbstractDimensions, AbstractUnits, AbstractUnitTransform
export ConversionError, DimensionError, NotScalarError, NotDimensionError, FixRat32
export Dimensions, Units, Quantity, FlexQuant, AffineTransform, NoTransform 
export StaticDims, StaticUnits
export RegistryTools, UnitRegistry
export static_fieldnames, uscale, uoffset, dimension, pretty_print_units
export assert_scalar, assert_dimension, assert_dimensionless
export with_ubase, ustrip, dstrip, ustrip_base, unit, quantity
export ubase, uconvert, dconvert, ustrip_dimensionless, udynamic
export LinmapQuant, FactorQuant, FunctionQuant, UnitMap, DimsMap, RepDimsMap, SymDimsMap
export uinput, uoutput


end
