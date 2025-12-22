const DEFAULT_RATIONAL = FixRat32
const DEFAULT_USYMBOL = :_

abstract type AbstractUnitLike end
abstract type AbstractDimLike <: AbstractUnitLike end
abstract type AbstractDimensions{P} <: AbstractDimLike end
abstract type AbstractUnits{D<:AbstractDimensions} <: AbstractUnitLike end
abstract type AbstractAffineUnits{D<:AbstractDimensions} <: AbstractUnits{D} end 

const AbstractAffineLike{D} = Union{D, AbstractAffineUnits{D}} where D <: AbstractDimensions
#Base.@assume_effects :consistent static_fieldnames(t::Type) = Base.fieldnames(t)
static_fieldnames(t::Type) = Base.fieldnames(t)
Base.eltype(::Type{<:AbstractDimensions{P}}) where P = P

#=======================================================================================
AbstractDimensions API
=======================================================================================#
#Dimension constructor from other types of dimensions
function (::Type{D})(x::AbstractDimensions) where {P, D<:AbstractDimensions{P}}
    return D(map(Base.Fix1(getproperty, x), static_fieldnames(D))...)
end

#Assign a single value to all dimensions
function (::Type{D})(x::Union{Real,Missing}) where {P, D<:AbstractDimensions{P}}
    f(anything) = x
    return D(map(f, static_fieldnames(D))...)
end

uscale(u::AbstractDimLike) = 1 # All AbstractDimensions have unity scale
uoffset(u::AbstractDimLike) = 0 # All AbstractDimensions have no offset
dimension(u::AbstractDimLike) = u
usymbol(u::AbstractDimLike) = DEFAULT_USYMBOL
dimtype(::Type{<:AbstractUnits{D}}) where D = D
dimtype(::Type{D}) where D<:AbstractDimLike = D
dimpowtype(::Type{D}) where {P, D<:AbstractDimensions{P}} = P
Base.getindex(d::AbstractDimensions, k::Symbol) = getproperty(d, k)
dimpowtype(::Type{U}) where {U<:AbstractUnitLike} = dimpowtype(dimtype(U))

#=======================================================================================
Basic SI dimensions
=======================================================================================#
"""
    Dimensions{P}

Basic SI dimensions:
    length = m, 
    mass = kg, 
    time = s, 
    current = A, 
    temperature = K, 
    luminosity = cd, 
    amount = mol
"""
@kwdef struct Dimensions{P} <: AbstractDimensions{P}
    length::P = FixedRational(0)
    mass::P = FixedRational(0)
    time::P = FixedRational(0)
    current::P = FixedRational(0)
    temperature::P = FixedRational(0)
    luminosity::P = FixedRational(0)
    amount::P = FixedRational(0)
end
const DEFAULT_DIMENSONS = Dimensions{DEFAULT_RATIONAL}
Dimensions(args...) = Dimensions{DEFAULT_RATIONAL}(args...)
#Dimensions(d::AbstractDimensions) = Dimensions{DEFAULT_RATIONAL}(d)

function unit_symbols(::Type{<:Dimensions})
    return Dimensions{Symbol}(
        length=:m, mass=:kg, time=:s, current=:A, temperature=:K, luminosity=:cd, amount=:mol
    )
end


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
AffineUnits(scale, offset, dims::AbstractUnits{D}, symbol=DEFAULT_USYMBOL) where {D<:AbstractDimensions} = AffineUnits(scale, offset, convert(D, dims), symbol)

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
        show(io, dimension(u); pretty)
        return print(io, ")")
    end
end


#=================================================================================================
Quantity types
The basic type is Quantity, which belongs to <:Any (hence it has no real hierarchy)
=================================================================================================#
abstract type AbstractQuantity{T,U} end

"""
    Quantity{T<:Any,U<:AbstractUnitLike}

Generic quantity type with fields `value` and `unit`
"""
struct Quantity{T<:Any,U<:AbstractUnitLike} <: AbstractQuantity{T,U}
    value :: T
    unit  :: U
end
Quantity{T}(x, u::AbstractUnitLike) where T = Quantity{T, typeof(u)}(x, u)

ustrip(q::Quantity) = q.value
unit(q::Quantity) = q.unit
dimension(q::Quantity) = dimension(unit(q))
unittype(::Type{<:AbstractQuantity{T,U}}) where {T,U} = U
dimtype(::Type{<:AbstractQuantity{T,U}}) where {T,U} = dimtype(U)

AffineUnits(scale, offset::Quantity, dims::AbstractDimensions, symbol=DEFAULT_USYMBOL) = AffineUnits(scale, ustrip(dims, offset), dims, symbol)
AffineUnits(scale, offset::Quantity, dims::AbstractUnits, symbol=DEFAULT_USYMBOL) = AffineUnits(scale, ustrip(dims, offset), dims, symbol)

"""
    quantity(x, u::AbstractUnitLike)

Overloading constructor for various quantities
"""
quantity(x, u::AbstractUnitLike) = Quantity(x, u)
quantity(x::Tuple, u::Tuple) = map(quantity, x, u)

"""
    constructorof(::Type{T}) where T = Base.typename(T).wrapper

Return the constructor of a type T{PS...} by default it only returns T (i.e. removes type parameters)
This function can be overloaded if custom behaviour is needed
"""
constructorof(::Type{T}) where T = Base.typename(T).wrapper
constructorof(::Type{<:Dimensions}) = Dimensions
constructorof(::Type{<:AffineUnits}) = AffineUnits 
constructorof(::Type{<:Quantity}) = Quantity


"""
    MirrorDims

A dimension that represents a placeholder value that mirrors any dimension that is combined
with it (useful for initialization when units are unknown). For example 

julia> 1u"m/s" + 0*MirrorDims()
1 m/s

julia> max(1u"m/s", -Inf*MirrorDims())
1 m/s
"""
struct MirrorDims{D<:AbstractDimensions} <: AbstractDimLike end
MirrorDims() = MirrorDims{FixRat32, Dimensions{FixRat32}}()
MirrorDims(::Type{D}) where {D<:AbstractDimensions} = MirrorDims{D}()


const MirrorUnion{D} = Union{D, MirrorDims{D}}
promote_rule(::Type{D}, ::Type{<:MirrorDims}) where {D<:AbstractDimensions} = MirrorUnion{D}
function nomirror(x::Quantity)
    u = unit(x)
    return (u isa MirrorDims) ? throw(ArgumentError("Mirror dimensions found, cannot convert to non-mirror version")) : Quantity(ustrip(x), u)
end

#Quantities with mirror dimensions should include a union
Quantity(x::T, u::MirrorDims{D}) where {T,D<:AbstractDimensions} = Quantity{T, MirrorUnion{D}}(x, u)
Quantity{<:Any, <:MirrorDims}(x, u) = error("MirrorDims should not be a type parameter in a Quantity constructor. Use Quantity{T, MirrorUnion{D}}")

function Base.show(io::IO, d::MirrorDims{D}; pretty=PRETTY_DIM_OUTPUT[]) where {D<:AbstractDimensions}
    return print(io, "MirrorDims{$(D)}()")
end

function Base.show(io::IO, ::Type{MirrorDims{D}}; pretty=PRETTY_DIM_OUTPUT[]) where {D<:AbstractDimensions}
    return print(io, "MirrorDims{$(D)}")
end

function Base.show(io::IO, ::Type{MirrorUnion{D}}; pretty=PRETTY_DIM_OUTPUT[]) where {D<:AbstractDimensions}
    return print(io, "MirrorUnion{$(D)}")
end


"""
    UnitfulCallable{T<:Any, UI<:Any, UO<:Any}

An object representing a callable object that requires inputs to be units UI and produces outputs of UO
This is useful for wrapping functions/models that aren't generic enough to support `Quantity` types.

Example1: Applying quantities to a strictly Real function

    angle_coords(θ::Real, r::Real) = r.*(cos(θ), sin(θ))
    unitful_angle_coords = UnitfulCallable(angle_coords, (u"", u"m") => dimension(u"m"))

    c = unitful_angle_coords(30u"deg", 6u"cm")

Example2: Applying quantities to a strictly Real model

    struct RotatingArm
        len :: Float64
    end

    angle_coords(arm::RotatingArm, θ::Real) = arm.len.*(cos(θ), sin(θ))
    angle_coords(arm::Quantity{<:RotatingArm}, θ::Quantity{<:Real}) = UnitfulCallable(Base.Fix1(angle_coords, ustrip(arm)), u""=>unit(arm))(θ)

    c = angle_coords(Quantity(RotatingArm(1.0), u"m"), 30u"deg")
"""
struct UnitfulCallable{T, UI, UO}
    caller :: T 
    units  :: Pair{UI, UO}
end
UnitfulCallable(p::Pair) = UnitfulCallable(nothing, p)

(bb::UnitfulCallable{F})(x...) where F = _apply_unit_pair(bb.caller, bb.units, x...)
unitful_call(f, bb::UnitfulCallable, x...)  = _apply_unit_pair(f, bb.units, x...)

function _apply_unit_pair(f, u::Pair)
    return quantity(f(), u[2])
end

function _apply_unit_pair(f, u::Pair, x)
    (ui, uo) = u
    raw_args = ustrip(ui, x)
    return quantity(f(raw_args), uo)
end

function _apply_unit_pair(f, u::Pair, x1, xs...)
    (ui, uo) = u
    raw_args = map(ustrip, ui, (x1, xs...))
    return quantity(f(raw_args...), uo)
end



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
DimensionError(arg1, arg2, args...) = DimensionError((arg1, arg2, args...))
Base.showerror(io::IO, e::DimensionError{<:Tuple}) = print(io, "DimensionError: ", e.items, " have incompatible dimensions")
Base.showerror(io::IO, e::DimensionError{<:AbstractQuantity}) = print(io, "DimensionError: ", e.items, " is not dimensionless")
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
assert_dimensionless(q::AbstractQuantity) = isdimensionless(unit(q)) ? q : throw(DimensionError(q))
dimensionless(u::AbstractUnitLike) = dimension(assert_dimensionless(u))
dimensionless(q::AbstractQuantity) = ustrip(assert_dimensionless(ubase(q)))
dimensionless(n::Number) = n

isdimensionless(u::AbstractUnitLike) = iszero(dimension(u))
Base.iszero(u::D) where D<:AbstractDimensions = (u == D(0))
