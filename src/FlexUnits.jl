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
include("LogUnitRegistry.jl")

export AbstractUnitLike, AbstractDimLike, AbstractDimensions, AbstractUnits, AbstractUnitTransform
export ConversionError, DimensionError, NotScalarError, NotDimensionError, LogLinearError, FixRat32
export Dimensions, Units, StaticDims, Quantity, FlexQuant, LogQuant, QuantUnion, AffineTransform, NoTransform 
export RegistryTools, UnitRegistry
export static_fieldnames, uscale, uoffset, dimension, pretty_print_units
export assert_scalar, assert_dimension, assert_dimensionless
export with_ubase, ustrip, dstrip, ustrip_base, unit, quantity
export ubase, uconvert, dconvert, ustrip_dimensionless, udynamic, ustatic, logubase
export LinmapQuant, VectorQuant, FactorQuant, FunctionQuant, UnitMap, DimsMap
export uinput, uoutput, ufactor


end
