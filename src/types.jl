const DEFAULT_RATIONAL = FixRat32
const DEFAULT_USYMBOL = :_

"""
    AbstractUnitLike

A broad class representing anything that can be interpreted as a unit or dimension
"""
abstract type AbstractUnitLike end

"""
    AbstractUnitTransform

An abstract object representing a unit conversion formula. 
Any object that subtypes this is made callable.

```julia
# Callable form 
utrans = uconvert(u"°C", u"°F")
utrans(0.0)
31.999999999999986

# Shorthand callable form (syntactic sugar)
(u"°C" |> u"°F")(0.0)
31.999999999999986
```
"""
abstract type AbstractUnitTransform end

"""
    AbstractUnitLike

A broad class representing anything that can be interpreted as a unit
"""
abstract type AbstractDimLike <: AbstractUnitLike end

"""
    AbstractDimensions

A class that represents a specific dimensional schema, by default, its
only child is Dimensions (SI units), but users can build their own versions 
(even based on imperial measurements, or where angles are a dimension)
"""
abstract type AbstractDimensions{P} <: AbstractDimLike end

"""
    AbstractUnits

All units in FlexDims contain two entities: 
(1) AbstractDimLike (anything that can be interpreted as a dimension)
(2) AbstractTransform (contains the formula to convert unit value to dimensions)
"""
abstract type AbstractUnits{D, T<:AbstractUnitTransform} <: AbstractUnitLike end


const UnitOrDims{D} = Union{D, AbstractUnits{D}} where D<:AbstractDimensions
const ScalarOrVec{T} = Union{T, AbstractVector{T}} where T
abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}} <: AbstractMatrix{U} end

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

dimension(u::AbstractDimLike) = u
dimval(u::AbstractDimensions) = u
usymbol(u::AbstractDimLike) = DEFAULT_USYMBOL
uscale(u::AbstractUnitLike) = uscale(todims(u))
uoffset(u::AbstractUnitLike) = uoffset(todims(u))
dimtype(::Type{<:AbstractUnits{D}}) where D = D
dimtype(::Type{D}) where D<:AbstractDimLike = D
dimpowtype(::Type{D}) where {P, D<:AbstractDimensions{P}} = P
Base.getindex(d::AbstractDimensions, k::Symbol) = getproperty(d, k)
dimpowtype(::Type{U}) where {U<:AbstractUnitLike} = dimpowtype(dimtype(U))
is_dimension(u::AbstractUnitLike) = is_identity(todims(u))
is_scalar(u::AbstractUnitLike) = is_scalar(todims(u))
udynamic(u::AbstractDimensions) = u

#=======================================================================================
Basic SI dimensions and transforms
=======================================================================================#
"""
NoTransform object, the default transform returned by todims(x::AbstractDimensionLike). Calling it results in 
an identity.
```julia
t = NoTransform()
t("anything")
"anything"
```
"""
struct NoTransform <: AbstractUnitTransform end 
Base.broadcastable(utrans::AbstractUnitTransform) = Ref(utrans)
(t2::AbstractUnitTransform)(t1::AbstractUnitTransform) = t2 ∘ t1
(t::NoTransform)(x) = x
(t::NoTransform)(t0::AbstractUnitTransform) = t0

Base.:∘(t1::NoTransform, t2::AbstractUnitTransform) = t2 
Base.:∘(t1::AbstractUnitTransform, t2::NoTransform) = t1 
Base.:∘(t1::NoTransform, t2::NoTransform) = t1
Base.inv(t::NoTransform) = t
uscale(t::NoTransform)  = 1
uoffset(t::NoTransform) = 0
todims(u::AbstractDimLike) = NoTransform()
is_identity(t::NoTransform) = true
is_scalar(t::NoTransform) = true

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
Dimensions(args...) = Dimensions{FixRat32}(args...)

function unit_symbols(::Type{<:Dimensions})
    return Dimensions{Symbol}(
        length=:m, mass=:kg, time=:s, current=:A, temperature=:K, luminosity=:cd, amount=:mol
    )
end

"""
    dimension_names(::Type{<:AbstractDimensions})

Overloadable fieldnames for a dimension (returns a tuple of the dimension's field names)
This should be static so that it can be hardcoded during compilation. 
Can use this to overload the default "fieldnames" behaviour
"""
@inline function dimension_names(::Type{D}) where {D<:AbstractDimensions}
    return static_fieldnames(D)
end


"""
    StaticDims{D}

Static dimensions where the `D`` is the dimension value. This improves performance when dimensions are
statically inferrable.
"""
struct StaticDims{D} <: AbstractDimLike
    function StaticDims{D}() where D
        return (D isa AbstractDimensions) ? new{D}() : error("Type parameter must be a dimension")
    end
end 
StaticDims(D::AbstractDimensions) = StaticDims{D}()
StaticDims{D}(d::AbstractDimensions) where D = (D == d) ? StaticDims{D} : throw(ArgumentError("Dimesion $(d) must be equal to $(D)"))
dimtype(::Type{<:StaticDims{D}}) where D = typeof(D)
dimtype(d::StaticDims) = dimtype(typeof(d))
dimval(::Type{<:StaticDims{D}}) where D = D
dimval(d::StaticDims) = dimval(typeof(d))
udynamic(u::StaticDims{D}) where D = D


#=======================================================================================
Affine Units and Transforms
=======================================================================================#
"""
    AffineTransform

A type representing an affine transfomration formula that can be
used to convert values from one affine unit to another. This object is callable.

# Fields
- scale :: Float64
- offset :: Float64

# Constructors
```
AffineTransform(scale::Real, offset::Real)
AffineTransform(; scale, offset)
```
"""
@kwdef struct AffineTransform <: AbstractUnitTransform
    scale  :: Float64 = 1
    offset :: Float64 = 0
end
(t::AffineTransform)(x) = muladd(x, t.scale, t.offset)
(t::AffineTransform)(x::AbstractArray) = t.(x)
(t::AffineTransform)(x::Tuple) = map(t, x)
(t::AffineTransform)(t0::AbstractUnitTransform) = t ∘ t0

function Base.:∘(t2::AffineTransform, t1::AffineTransform)
    return AffineTransform(
        scale  = t1.scale*t2.scale,
        offset = t2.offset + t2.scale*t1.offset 
    )
end
Base.inv(t::AffineTransform) = AffineTransform(scale=inv(t.scale), offset=-t.offset/t.scale)

uscale(t::AffineTransform) = t.scale 
uoffset(t::AffineTransform) = t.offset
is_identity(t::AffineTransform) = isone(t.scale) & iszero(t.offset)
is_scalar(t::AffineTransform) = iszero(t.offset)
remove_offset(t::AffineTransform) = AffineTransform(scale=t.scale, offset=0)

"""
    Units{D<:AbstractDimensions, T<:AbstractUnitTransform}(dims::D, todims::T, symbol::Symbol)

A dynamic unit object that contains dimensions (dims) and its conversion formula to said dimensions (todims). The conversion 
formula determines what kind of unit is referred to. An AffineTransform implies affine units, a NoTransform implies dimensions.
Dynamic units can generated through the `@ud_str` macro.

```julia
julia> 1*(5ud"°C") #Operations on units eagerly convert to dimensions
278.15 K

julia> 1*(5du"°C") |> ud"°C" #Converts operation results back to Celsius
5.0 °C

julia> (5ud"°C" + 2ud"°C") |> ud"°C" #Operation adds values in Kelvin, results converted back to Celsius
280.15 °C

julia> (ustrip(5ud"°C") + ustrip(2ud"°C"))*u"°C" #Strips, adds raw quantity values, converts raw number to Celsius
7 °C
```
"""
@kwdef struct Units{D<:AbstractDimensions, T<:AbstractUnitTransform} <: AbstractUnits{D, T}
    dims   :: D
    todims :: T
    symbol :: Symbol = DEFAULT_USYMBOL
end
Units{D}(units, todims::T, symbol=DEFAULT_USYMBOL) where {D,T<:AbstractUnitTransform} = Units{D,T}(units, todims, symbol)
Units(dims::D, todims::AbstractUnitTransform=NoTransform(), symbol=DEFAULT_USYMBOL) where D<:AbstractDimensions = Units(dims, todims, symbol)
Units(units::D, todims::AbstractUnitTransform, symbol=DEFAULT_USYMBOL) where D<:AbstractUnits = Units(dimension(assert_dimension(units)), todims, symbol)

todims(u::Units) = u.todims
dimension(u::Units) = u.dims 
usymbol(u::Units) = u.symbol
remove_offset(u::U) where U<:AbstractUnits = constructorof(U)(dimension(u), remove_offset(u.todims))
is_scalar(u::AbstractUnits) = is_scalar(todims(u))
is_dimension(u::AbstractUnits) = is_identity(todims(u))
#affine_units(;dims, scale=1, offset=0, symbol=DEFAULT_USYMBOL) = Units(dims=dims, todims=AffineTransform(scale=scale, offset=offset), symbol=symbol)
udynamic(u::Units) = u
dimtype(::Type{Units{D,C}}) where {D,C} = D
dimtype(d::Units) = dimtype(typeof(d))


"""
    StaticUnits{D, T<:AbstractUnitTransform}(todims::T, symbol::Symbol)

A static version of units, where the value of dimensions `D` is a a parameter.
Static units can generated through the `@u_str` macro. This improves performance when
dimensions are statically inferrable.
"""
@kwdef struct StaticUnits{D, T<:AbstractUnitTransform} <: AbstractUnits{D,T}
    todims :: T
    symbol :: Symbol
    function StaticUnits{D,C}(conv::AbstractUnitTransform, symb=DEFAULT_USYMBOL::Symbol) where {D, C<:AbstractUnitTransform}
        return (D isa AbstractDimensions) ? new{D,C}(conv, symb) : error("Type parameter must be a dimension")
    end
    function StaticUnits{D}(conv::C, symb=DEFAULT_USYMBOL::Symbol) where {D, C<:AbstractUnitTransform}
        return (D isa AbstractDimensions) ? new{D,C}(conv, symb) : error("Type parameter must be a dimension")
    end
end
StaticUnits(u::Units) = StaticUnits{dimension(u)}(todims(u), usymbol(u))
StaticUnits(d::AbstractDimensions, todims::AbstractUnitTransform, symb=DEFAULT_USYMBOL) = StaticUnits{d}(todims, symb)
StaticUnits(d::StaticDims{D}, todims::AbstractUnitTransform, symb=DEFAULT_USYMBOL) where D = StaticUnits{D}(todims, symb)

constructorof(::Type{<:StaticUnits}) = StaticUnits
Units(u::StaticUnits) = Units{dimtype(u)}(dimval(u), todims(u), usymbol(u))
udynamic(u::StaticUnits) = Units(u)
todims(u::StaticUnits) = u.todims
dimtype(::Type{StaticUnits{D,C}}) where {D,C} = typeof(D)
dimtype(d::StaticUnits) = dimtype(typeof(d))
dimval(::Type{StaticUnits{D,C}}) where {D,C} = D
dimval(d::StaticUnits) = dimval(typeof(d))
dimension(::Type{StaticUnits{D,T}}) where {D,T} = StaticDims{D}()
dimension(d::StaticUnits) = dimension(typeof(d))
usymbol(u::StaticUnits) = u.symbol

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
Quantity(x::T, u::StaticUnits{D}) where {T,D} = Quantity(u.todims(x), StaticDims{D}())
Quantity{T}(x, u::StaticUnits{D}) where {T,D} = Quantity{T}(convert(T, u.todims(x)), StaticDims{D}())
Quantity{T}(q::AbstractQuantity) where T = Quantity{T}(ustrip(q), unit(q))
Quantity{T,U}(q::AbstractQuantity) where {T,U} = Quantity{T,U}(ustrip(q), unit(q))

ustrip(q::Quantity) = q.value
unit(q::Quantity) = q.unit
dimension(q::Quantity) = dimension(unit(q))
unittype(::Type{<:AbstractQuantity{T,U}}) where {T,U} = U
dimtype(::Type{<:AbstractQuantity{T,U}}) where {T,U} = dimtype(U)
dimtype(q::Quantity) = dimtype(unit(q))
udynamic(q::Quantity) = Quantity(ustrip(q), udynamic(unit(q)))

AffineTransform(scale::Real, offset::Quantity) = AffineTransform(scale=scale, offset=dstrip(offset))


"""
    constructorof(::Type{T}) where T = Base.typename(T).wrapper

Return the constructor of a type `T{PS...}` by default it only returns `T` (i.e. removes type parameters)
This function can be overloaded if custom behaviour is needed
"""
constructorof(::Type{T}) where T = Base.typename(T).wrapper
constructorof(::Type{<:Dimensions}) = Dimensions
constructorof(::Type{<:Units}) = Units
constructorof(::Type{<:Quantity}) = Quantity


"""
    MirrorDims

A dimension that represents a placeholder value that mirrors any dimension that is combined
with it (useful for initialization when units are unknown). For example 

```julia
julia> 1u"m/s" + 0*MirrorDims()
1 m/s

julia> max(1u"m/s", -Inf*MirrorDims())
1 m/s
```
"""
struct MirrorDims{D<:AbstractDimensions} <: AbstractDimLike end
MirrorDims() = MirrorDims{Dimensions{FixRat32}}()
MirrorDims(::Type{D}) where {D<:AbstractDimensions} = MirrorDims{D}()


const MirrorUnion{D} = Union{D, MirrorDims{D}} where D<:AbstractDimensions
Base.promote_rule(::Type{D}, ::Type{<:MirrorDims}) where {D<:AbstractDimensions} = MirrorUnion{D}
function nomirror(x::Quantity)
    u = unit(x)
    return (u isa MirrorDims) ? throw(ArgumentError("Mirror dimensions found, cannot convert to non-mirror version")) : Quantity(ustrip(x), u)
end

#Quantities with mirror dimensions should include a union
Quantity(x::T, u::MirrorDims{D}) where {T, D} = Quantity{T, MirrorUnion{D}}(x, u)
Quantity{T, MirrorDims{D}}(x, u) where {T, D} = error("MirrorDims should not be a type parameter in a Quantity constructor. Use Quantity{T, MirrorUnion{D}}")



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
    pretty_str(x) = (ushow(io_tmp, x, pretty=true); String(take!(io_tmp)))

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


assert_scalar(u::AbstractDimLike) = u
assert_scalar(u::AbstractUnits) = is_scalar(u) ? u : throw(NotScalarError(u))
scalar_dimension(u::AbstractUnitLike) = dimension(assert_scalar(u))

assert_dimension(u::AbstractDimLike) =  u
assert_dimension(u::AbstractUnits) = is_dimension(u) ? u : throw(NotDimensionError(u))

assert_dimensionless(u::AbstractUnitLike) = isdimensionless(u) ? u : throw(DimensionError(u))
assert_dimensionless(q::AbstractQuantity) = isdimensionless(unit(q)) ? q : throw(DimensionError(q))
dimensionless(u::AbstractUnitLike) = dimension(assert_dimensionless(u))
dimensionless(q::AbstractQuantity) = ustrip(assert_dimensionless(ubase(q)))
dimensionless(n::Number) = n

isdimensionless(u::AbstractUnitLike) = iszero(dimension(u))
Base.iszero(u::D) where D<:AbstractDimensions = (u == D(0))
Base.iszero(u::StaticDims{d}) where d = iszero(d)

# This is deprecated in favor of QuantMapping, make sure these cases work
#=
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
    return Quantity(f(), u[2])
end

function _apply_unit_pair(f, u::Pair, x)
    (ui, uo) = u
    raw_args = ustrip(ui, x)
    return Quantity(f(raw_args), uo)
end

function _apply_unit_pair(f, u::Pair, x1, xs...)
    (ui, uo) = u
    raw_args = map(ustrip, ui, (x1, xs...))
    return Quantity(f(raw_args...), uo)
end
=#
