const DEFAULT_PWR_TYPE = FixedRational{DEFAULT_NUMERATOR_TYPE,DEFAULT_DENOM}

abstract type AbstractUnitLike end
abstract type AbstractDimensions{P} <: AbstractUnitLike end
abstract type AbstractUnits{D<:AbstractDimensions} end
abstract type AbstractAffineUnits{D<:AbstractDimensions} <: AbstractUnits{D} end 
abstract type AbstractScalarUnits{D<:AbstractDimensions} <: AbstractAffineUnits{D} end

@kwdef struct SIDims{P} <: AbstractDimensions{P}
    length::P = 0
    mass::P = 0
    time::P = 0
    current::P = 0
    temperature::P = 0
    luminosity::P = 0
    amount::P = 0
end

SIDims(;kwargs...) = SIDims{DEFAULT_PWR_TYPE}(;kwargs...)
SIDims(args...) = SIDims{DEFAULT_PWR_TYPE}(args...)
DEFAULT_DIMENSONS = SIDims{DEFAULT_PWR_TYPE}

uscale(u::AbstractDimensions) = 1
uoffset(u::AbstractDimensions) = 0
dimension(u::AbstractDimensions) = u
usymbol(u::AbstractDimensions) = :nothing



Base.@pure static_fieldnames(t::Type) = Base.fieldnames(t)

"""
    dimension_names(::Type{<:AbstractDimensions})

Overloadable feieldnames for a dimension (returns a tuple of the dimensions field names)
This should be static so that it can be hardcoded during compilation. 
Can use this to overload the default "fieldnames" behaviour
"""
@inline function dimension_names(::Type{D}) where {D<:AbstractDimensions}
    return static_fieldnames(D)
end


@kwdef struct ScalarUnits{D} <: AbstractScalarUnits{D}
    scale::Float64
    dims::D
    symbol::Symbol=:nothing
end

ScalarUnits(scale, dims::D, symbol=:nothing) where {D} = ScalarUnits{D}(scale, dims, symbol)

uscale(u::ScalarUnits) = u.scale
uoffset(u::AbstractScalarUnits) = 0
dimension(u::ScalarUnits) = u.dims
usymbol(u::ScalarUnits) = u.symbol


@kwdef struct AffineUnits{D} <: AbstractAffineUnits{D}
    scale::Float64
    offset::Float64
    dims::D
    symbol::Symbol=:nothing 
end

AffineUnits(scale, offset, dims::D, symbol=:nothing) where {D<:AbstractDimensions} = AffineUnits{D}(scale, offset, dims, symbol)

uscale(u::AffineUnits) = u.scale
uoffset(u::AffineUnits) = u.offset 
dimension(u::ScalarUnits) = u.dims 
usymbol(u::ScalarUnits) = u.symbol


abstract type AbstractQuantity{T,U<:AbstractUnits} <: Number end

struct Quantity{T<:Number,U<:AbstractUnitLike}
    value :: T
    units :: U 
end

ustrip(q::Quantity) = q.value
unit(q::Quantity) = q.units
dimension(q::Quantity) = dimension(unit(q))


"""
    DimensionError{D} <: Exception

Error thrown when an operation is dimensionally invalid given the arguments
"""
struct DimensionError{Q1,Q2} <: Exception
    q1::Q1
    q2::Q2
    DimensionError(q1::Q1, q2::Q2) where {Q1,Q2} = new{Q1,Q2}(q1, q2)
    DimensionError(q1) = DimensionError(q1, nothing)
end

Base.showerror(io::IO, e::DimensionError) = print(io, "DimensionError: ", e.q1, " and ", e.q2, " have incompatible dimensions")
Base.showerror(io::IO, e::DimensionError{<:Any,Nothing}) = print(io, "DimensionError: ", e.q1, " is not dimensionless")

"""
    AffineOffsetError{D} <: Exception

Error thrown when an operation is invalid on affine units (with an offset)
"""
struct AffineOffsetError{D} <: Exception
    dim::D
    AffineOffsetError(dim) = new{typeof(dim)}(dim)
end

Base.showerror(io::IO, e::AffineOffsetError) = print(io, "AffineOffsetError: ", e.dim, " has a non-zero offset, operation not allowed on affine units. Consider using `baseunits(x)` to explicitly convert")
