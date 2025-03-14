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


"""
    dimension_names(::Type{<:AbstractDimensions})

Overloadable feieldnames for a dimension (returns a tuple of the dimensions field names)
This should be static so that it can be hardcoded during compilation. 
Can use this to overload the default "fieldnames" behaviour
"""
@inline function dimension_names(::Type{D}) where {D<:AbstractDimensions}
    return static_fieldnames(D)
end


abstract type AbstractScalarUnits{D<:AbstractDimensions,T<:Number} <: AbstractUnits end

@kwdef struct ScalarUnits{D,T} <: AbstractScalarUnits{D,T}
    scale::T
    dims::D
    symbol::Symbol=:nothing
end

ScalarUnits(scale::T, dims::D, symbol=:nothing) where {T, D<:AbstractDimensions} = ScalarUnits{D,promote_type(T,Float64)}(scale, dims, symbol)

abstract type AbstractAffineUnits{D<:AbstractDimensions, T<:Number} <: AbstractUnits end 

@kwdef struct AffineUnits{D,T} <: AbstractScalarUnits{D,T}
    scale::T 
    offset::T 
    dims::D 
    symbol::Symbol=:nothing 
end

AffineUnits(scale::T1, offset::T2, dims::D, symbol=:nothing) where {T1, T2, D<:AbstractDimensions} = ScalarUnits{D,promote_type(T1,T2,Float64)}(scale, offset, dims, symbol)

"""
    DimensionError{D} <: Exception

Error thrown when an operation is dimensionally invalid given the arguments
"""
struct DimensionError{Q1,Q2} <: Exception
    q1::Q1
    q2::Q2
    DimensionError(q1, q2) = new{typeof(q1),typeof(q2)}(q1, q2)
    DimensionError(q1) = DimensionError(q1, nothing)
end

