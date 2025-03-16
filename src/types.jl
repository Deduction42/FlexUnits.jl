const DEFAULT_PWR_TYPE = FixedRational{DEFAULT_NUMERATOR_TYPE,DEFAULT_DENOM}
const DEFAULT_SYMBOL = :_

abstract type AbstractUnitLike end
abstract type AbstractDimensions{P} <: AbstractUnitLike end
abstract type AbstractUnits{D<:AbstractDimensions} <: AbstractUnitLike end
abstract type AbstractAffineUnits{D<:AbstractDimensions} <: AbstractUnits{D} end 
abstract type AbstractScalarUnits{D<:AbstractDimensions} <: AbstractAffineUnits{D} end

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
dimension(u::AbstractDimensions) = u
usymbol(u::AbstractDimensions) = DEFAULT_SYMBOL



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
    scale::Float64 = 1
    dims::D
    symbol::Symbol=DEFAULT_USYMBOL
end

ScalarUnits(scale, dims::D, symbol=DEFAULT_SYMBOL) where {D} = ScalarUnits{D}(scale, dims, symbol)

uscale(u::ScalarUnits) = u.scale
uoffset(u::AbstractScalarUnits) = 0
dimension(u::ScalarUnits) = u.dims
usymbol(u::ScalarUnits) = u.symbol


@kwdef struct AffineUnits{D} <: AbstractAffineUnits{D}
    scale::Float64 = 1
    offset::Float64 = 0
    dims::D
    symbol::Symbol=DEFAULT_SYMBOL 
end

AffineUnits(scale, offset, dims::D, symbol=DEFAULT_SYMBOL) where {D<:AbstractDimensions} = AffineUnits{D}(scale, offset, dims, symbol)

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
    constructorof(::Type{T}) where T = Base.typename(T).wrapper

Return the constructor of a type T{PS...} by default it only returns T (i.e. removes type parameters)
This function can be overloaded if custom behaviour is needed
"""
constructorof(::Type{T}) where T = Base.typename(T).wrapper
constructorof(::Type{<:Dimensions}) = Dimensions
constructorof(::Type{<:ScalarUnits}) = ScalarUnits
constructorof(::Type{<:AffineUnits}) = AffineUnits
constructorof(::Type{<:Quantity}) = Quantity


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
    ScalarUnitError{D} <: Exception

Error thrown for non-scalar units when the operation is only valid for scalar units
"""
struct ScalarUnitError{D} <: Exception
    dim::D
    ScalarUnitError(dim) = new{typeof(dim)}(dim)
end

Base.showerror(io::IO, e::ScalarUnitError) = print(io, "ScalarUnitError: ", e.dim, " cannot be treated as scalar, operation only valid for scalar units")
