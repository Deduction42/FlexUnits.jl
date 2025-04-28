const DEFAULT_RATIONAL = FixedRational{DEFAULT_DENOM, DEFAULT_NUMERATOR_TYPE}
const DEFAULT_USYMBOL = :_

abstract type AbstractUnitLike end
abstract type AbstractDimensions{P} <: AbstractUnitLike end
abstract type AbstractUnits{D<:AbstractDimensions} <: AbstractUnitLike end
abstract type AbstractAffineUnits{D<:AbstractDimensions} <: AbstractUnits{D} end 

const AbstractAffineLike{D} = Union{D, AbstractAffineUnits{D}} where D <: AbstractDimensions
Base.@pure static_fieldnames(t::Type) = Base.fieldnames(t)

#Dimension constructor from other types of dimensions
function (::Type{D})(x::AbstractDimensions) where {P, D<:AbstractDimensions{P}}
    return D(map(Base.Fix1(getproperty, x), static_fieldnames(D))...)
end

@kwdef struct Dimensions{P} <: AbstractDimensions{P}
    length::P = FixedRational(0)
    mass::P = FixedRational(0)
    time::P = FixedRational(0)
    current::P = FixedRational(0)
    temperature::P = FixedRational(0)
    luminosity::P = FixedRational(0)
    amount::P = FixedRational(0)
end

Dimensions(args...) = Dimensions{DEFAULT_RATIONAL}(args...)
Dimensions(d::AbstractDimensions) = Dimensions{DEFAULT_RATIONAL}(d)
DEFAULT_DIMENSONS = Dimensions{DEFAULT_RATIONAL}

uscale(u::AbstractDimensions) = 1 # All AbstractDimensions have unity scale
uoffset(u::AbstractDimensions) = 0 # All AbstractDimensions have no offset
dimension(u::AbstractDimensions) = u
usymbol(u::AbstractDimensions) = DEFAULT_USYMBOL
Base.getindex(d::AbstractDimensions, k::Symbol) = getproperty(d, k)
dimtype(::Type{<:AbstractUnits{D}}) where D = D
dimtype(::Type{D}) where D<:AbstractDimensions = D

function unit_symbols(::Type{<:Dimensions})
    return Dimensions{Symbol}(
        length=:m, mass=:kg, time=:s, current=:A, temperature=:K, luminosity=:cd, amount=:mol
    )
end

"""
    NoDims{P}

A unitless dimension, compatible with any other dimension

Calling `getproperty` will always return `zero(P)`
promote(Type{<:NoDims}, D<:AbstractDimension) will return D 
convert(Type{D}, NoDims) where D<:AbstractDimensions will return D()
"""
struct NoDims{P} <: AbstractDimensions{P} end
NoDims() = NoDims{DEFAULT_RATIONAL}()
Base.getproperty(::NoDims{P}, ::Symbol) where {P} = zero(P)
unit_symbols(::Type{<:NoDims}) = NoDims{Symbol}()


"""
    dimension_names(::Type{<:AbstractDimensions})

Overloadable feieldnames for a dimension (returns a tuple of the dimensions field names)
This should be static so that it can be hardcoded during compilation. 
Can use this to overload the default "fieldnames" behaviour
"""
@inline function dimension_names(::Type{D}) where {D<:AbstractDimensions}
    return static_fieldnames(D)
end

"""
    AffineUnits{D<:AbstractDimensions}(scale::Float64, offset::Float64, dims::D, symbol::Symbol)

Affine-dimensional unit (treated as a scalar when offset=0). Quantities with this unit are eagerly
converted to dimmensional quantities for any operation, WHICH MAY BE UNITUITIVE because operations
do not happen directly on values if there is an offset. If you want operations on the quantity 
values directly, simply use "ustrip" and convert back.

julia> 1*(5u"°C") #Operations convert to Kelvin
278.15 K

julia> 1*(5u"°C") |> u"°C" #Converts operation results back to Celsius
5.0 °C

julia> (5u"°C" + 2u"°C") |> u"°C" #Operation adds values in Kelvin, results converted back to Celsius
280.15 °C

julia> (ustrip(5u"°C") + ustrip(2u"°C"))*u"°C" #Strips, adds raw quantity values, converts raw number to Celsius
7 °C
"""
@kwdef struct AffineUnits{D<:AbstractDimensions} <: AbstractAffineUnits{D}
    scale::Float64 = 1
    offset::Float64 = 0
    dims::D
    symbol::Symbol=DEFAULT_USYMBOL 
end

AffineUnits(scale, offset, dims::D, symbol=DEFAULT_USYMBOL) where {D<:AbstractDimensions} = AffineUnits{D}(scale, offset, dims, symbol)
AffineUnits(scale, offset, dims::U, symbol=DEFAULT_USYMBOL) where {D,U<:AbstractUnits{D}} = AffineUnits{D}(scale, offset, dims, symbol)
#AffineUnits(u::AffineUnits) = u

uscale(u::AffineUnits) = u.scale
uoffset(u::AffineUnits) = u.offset 
dimension(u::AffineUnits) = u.dims 
usymbol(u::AffineUnits) = u.symbol
remove_offset(u::U) where U<:AbstractAffineUnits = constructorof(U)(scale=uscale(u), offset=0, dims=dimension(u))

function Base.show(io::IO, u::AffineUnits; pretty=PRETTY_DIM_OUTPUT[])
    if usymbol(u) != DEFAULT_USYMBOL
        return print(io, usymbol(u))
    else
        print(io, "AffineUnits(scale=", uscale(u), ", offset=", uoffset(u), ", dims=")
        if pretty
            show(io, dimension(u), pretty=true)
        else
            show(io, dimension(u), pretty=false)
        end
        return print(io, ")")
    end
end


#=================================================================================================
Quantity types
The basic type is Quantity, which belongs to <:Any (hence it has no real hierarchy)
Other types are "narrower" in order to slot into different parts of the number hierarchy
=================================================================================================#

struct Quantity{T<:Any,U<:AbstractUnitLike}
    value :: T
    units :: U
end
narrowest_quantity(::Type{<:Any}) = Quantity

struct NumberQuantity{T<:Number,U<:AbstractUnitLike} <: Number
    value :: T
    units :: U
end
narrowest_quantity(::Type{<:Number}) = NumberQuantity

struct RealQuantity{T<:Real,U<:AbstractUnitLike} <: Real
    value :: T
    units :: U
end
narrowest_quantity(::Type{<:Real}) = RealQuantity

narrowest_quantity(x::Any) = narrowest_quantity(typeof(x))

#=================================================================================================
# Generic unions of quantities and fallbacks
=================================================================================================#
const UnionQuantity{T,U} = Union{Quantity{T,U}, NumberQuantity{T,U}, RealQuantity{T,U}}

#Generic fallback constructors
(::Type{Q})(v0, u) where {T, Q<:UnionQuantity{T}} = constructorof(Q)(convert(T, v0), u)
(::Type{Q})(q::UnionQuantity) where Q<:UnionQuantity = Q(ustrip(q), unit(q))

ustrip(q::UnionQuantity) = q.value
unit(q::UnionQuantity) = q.units
dimension(q::UnionQuantity) = dimension(unit(q))


"""
    quantity(x, u::AbstractUnitLike)

Constructs a quantity based on the narrowest quantity type that accepts x as an agument.
If "NoDims" is passed, only x is returned
"""
quantity(x, u::AbstractUnitLike) = narrowest_quantity(x)(x, u)
quantity(x, u::NoDims) = x

"""
    narrowest(q::UnionQuantity)

Returns the narrowest quantity type of `q`
"""
narrowest(q::UnionQuantity) = quantity(ustrip(q), unit(q))

"""
    constructorof(::Type{T}) where T = Base.typename(T).wrapper

Return the constructor of a type T{PS...} by default it only returns T (i.e. removes type parameters)
This function can be overloaded if custom behaviour is needed
"""
constructorof(::Type{T}) where T = Base.typename(T).wrapper
constructorof(::Type{<:Dimensions})  = Dimensions
constructorof(::Type{<:AffineUnits}) = AffineUnits
constructorof(::Type{<:Quantity}) = Quantity
constructorof(::Type{<:RealQuantity}) = RealQuantity
constructorof(::Type{<:NumberQuantity}) = NumberQuantity

#=============================================================================================
Errors and assertion functions
=============================================================================================#

"""
    DimensionError{D} <: Exception

Error thrown when an operation is dimensionally invalid given the arguments
"""
struct DimensionError{T} <: Exception
    items :: T
end
Base.showerror(io::IO, e::DimensionError{<:Tuple}) = print(io, "DimensionError: ", e.items, " have incompatible dimensions")
Base.showerror(io::IO, e::DimensionError{<:UnionQuantity}) = print(io, "DimensionError: ", e.items, " is not dimensionless")
Base.showerror(io::IO, e::DimensionError{<:AbstractUnitLike}) = print(io, "DimensionError: ", e.items, " is not dimensionless")

"""
    ConversionError{U, U0} <: Exception 

Error thrown when trying to convert u0 to u, offers hint on how to make u compatible
"""
struct ConversionError{U<:AbstractUnitLike, U0<:AbstractUnitLike} <: Exception
    u::U
    u0::U0
end

function Base.showerror(io::IO, e::ConversionError{<:AbstractUnitLike, <:AbstractUnitLike})
    io_tmp = IOBuffer()
    pretty_str(x) = (_print_pretty_unit(io_tmp, x); String(take!(io_tmp)))

    uΔ = dimension(e.u0)/dimension(e.u)
    return print(io, "ConversionError: Cannot convert unit '", pretty_str(e.u0), "' to target unit '", pretty_str(e.u), "'. Consider multiplying '", pretty_str(e.u), "' by '", pretty_str(uΔ), "' or similar.")
end

"""
    NotScalarError{D} <: Exception

Error thrown for non-scalar units when the operation is only valid for scalar units
"""
struct NotScalarError{D} <: Exception
    dim::D
    NotScalarError(dim) = new{typeof(dim)}(dim)
end
Base.showerror(io::IO, e::NotScalarError) = print(io, "NotScalarError: ", e.dim, " cannot be treated as scalar, operation only valid for scalar units")


"""
    NotDimensionError{D} <: Exception

Error thrown for non-dimensional units (scaled or affine) when the operation is only valid for dimensional units
"""
struct NotDimensionError{D} <: Exception
    dim::D
    NotDimensionError(dim) = new{typeof(dim)}(dim)
end
Base.showerror(io::IO, e::NotDimensionError) = print(io, "NotDimensionError: ", e.dim, " cannot be treated as dimension, operation only valid for dimension units")


assert_scalar(u::AbstractDimensions)  = u
assert_scalar(u::AbstractAffineUnits) = iszero(uoffset(u)) ? u : throw(NotScalarError(u))
scalar_dimension(u::AbstractUnitLike) = dimension(assert_scalar(u))

assert_dimension(u::AbstractDimensions) =  u
assert_dimension(u::AbstractAffineUnits) = isone(uscale(u)) & iszero(uoffset(u)) ? u : throw(NotDimensionError(u))

assert_dimensionless(u::AbstractUnitLike) = isdimensionless(u) ? u : throw(DimensionError(u))
assert_dimensionless(q::UnionQuantity) = isdimensionless(unit(q)) ? q : throw(DimensionError(q))
dimensionless(u::AbstractUnitLike) = (assert_dimensionless(u); NoDims())
dimensionless(q::UnionQuantity) = ustrip(assert_dimensionless(ubase(q)))
dimensionless(n::Number) = n

function Base.iszero(u::U) where U<:AbstractDimensions
    zero_dimension(obj::AbstractDimensions, fn::Symbol) = iszero(getproperty(obj, fn))
    return all(Base.Fix1(zero_dimension, u), dimension_names(U)) 
end
isdimensionless(u::AbstractUnitLike) = iszero(dimension(u))