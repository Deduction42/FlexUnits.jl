const DEFAULT_PWR_TYPE = FixedRational{DEFAULT_NUMERATOR_TYPE,DEFAULT_DENOM}

abstract type AbstractUnits end
abstract type AbstractDimensions{P} <: AbstractUnits end

@kwdef struct Dimensions{P} <: AbstractDimensions{P}
    length::P = 0
    mass::P = 0
    time::P = 0
    current::P = 0
    temperature::P = 0
    luminosity::P = 0
    amount::P = 0
end

Dimensions(;kwargs...) = Dimensions{DEFAULT_PWR_TYPE}(;kwargs...)
Dimensions(args...) = Dimensions{DEFAULT_PWR_TYPE}(args...)
DEFAULT_DIMENSONS = Dimensions{DEFAULT_PWR_TYPE}

uscale(u::AbstractDimensions) = 1
uoffset(u::AbstractDimensions) = 0
udims(u::AbstractDimensions) = u
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



abstract type AbstractScalarUnits{D<:AbstractDimensions} <: AbstractUnits end

@kwdef struct ScalarUnits{D} <: AbstractScalarUnits{D}
    scale::Float64
    dims::D
    symbol::Symbol=:nothing
end

ScalarUnits(scale, dims::D, symbol=:nothing) where {D} = ScalarUnits{D}(scale, dims, symbol)

uscale(u::ScalarUnits) = u.scale
uoffset(u::AbstractScalarUnits) = 0
udims(u::ScalarUnits) = u.dims
usymbol(u::ScalarUnits) = u.symbol


abstract type AbstractAffineUnits{D<:AbstractDimensions} <: AbstractUnits end 

@kwdef struct AffineUnits{D} <: AbstractAffineUnits{D}
    scale::Float64
    offset::Float64
    dims::D
    symbol::Symbol=:nothing 
end

AffineUnits(scale, offset, dims::D, symbol=:nothing) where {D<:AbstractDimensions} = AffineUnits{D}(scale, offset, dims, symbol)

uscale(u::AffineUnits) = u.scale
uoffset(u::AffineUnits) = u.offset 
udims(u::ScalarUnits) = u.dims 
usymbol(u::ScalarUnits) = u.symbol


abstract type AbstractQuantity{T,U<:AbstractUnits} end

struct Quantity{T,U}
    value :: T
    units :: U 
end

ustrip(q::Quantity) = q.value
vstrip(q::Quantity) = q.units


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
