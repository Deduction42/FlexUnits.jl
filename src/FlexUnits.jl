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

# Set default units for Simplification
for ustr in ["F", "H", "T", "Ω", "V", "W", "J", "Pa", "N", "C", "L", "(m/s)"]
    set_preferred_unit(UnitRegistry.uparse(ustr))
end

export AbstractUnitLike, AbstractDimLike, AbstractDimensions, AbstractUnits, AbstractUnitTransform
export ConversionError, DimensionError, NotScalarError, NotDimensionError, FixRat32
export Dimensions, Units, StaticDims, LogScale, Quantity, FlexQuant, LogQuant, QuantUnion, AffineTransform, NoTransform, ExpAffTransform
export RegistryTools, UnitRegistry, LogUnitRegistry
export static_fieldnames, uscale, uoffset, dimension, pretty_print_units
export simplify, set_preferred_unit, set_preferred_logscale, display_simplified_units
export with_ubase, ustrip, dstrip, ustrip_base, unit, quantity, logquant, linquant
export ubase, uconvert, dconvert, scalar, udynamic, ustatic, logubase
export LinmapQuant, VectorQuant, FactorQuant, FunctionQuant, UnitMap, DimsMap
export uinput, uoutput, ufactor


end
